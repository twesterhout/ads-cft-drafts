module EinsteinEquations.Tensor
  ( Tensor (..),
    Device (..),
    flatten,
    reshape,
    tensorLength,
    newTensor,
    tensorSimilar,
    reversePrimArray,
    rowMajorStrides,
    generate1,
    tensorMap,
    tensorMult,
  )
where

import Control.Monad.Primitive (PrimBase, PrimMonad, touch, unsafeIOToPrim, unsafePrimToIO)
import Control.Monad.ST
-- import Foreign.Storable

import DLPack
import Data.Kind
import Data.Maybe (fromJust)
import Data.Primitive.PrimArray
import qualified Data.Primitive.Ptr as P
import Data.Primitive.Types (Prim, alignment, sizeOf)
import qualified Data.Text as T
import Foreign.ForeignPtr
import Foreign.Ptr (Ptr, plusPtr)
import qualified GHC.Exts as GHC
import GHC.ForeignPtr (mallocPlainForeignPtrAlignedBytes)
import GHC.TypeLits
import Halide
import Numeric.TBLIS
import System.IO.Unsafe

natToInt :: forall n. KnownNat n => Int
natToInt = fromIntegral $ GHC.TypeLits.natVal (Proxy @n)
{-# INLINE natToInt #-}

symbolToText :: forall s. KnownSymbol s => Text
symbolToText = toText $ GHC.TypeLits.symbolVal (Proxy @s)

data Device = CPU
  deriving stock (Read, Show, Eq)

data Tensor (device :: Device) (rank :: Nat) (a :: Type) = Tensor
  { tensorData :: {-# UNPACK #-} !(ForeignPtr a),
    tensorShape :: {-# UNPACK #-} !(PrimArray Int),
    tensorStrides :: {-# UNPACK #-} !(PrimArray Int)
  }
  deriving stock (Show)

indexSlice :: (HasCallStack, KnownNat r, 2 <= r) => Int -> Int -> Tensor d r a -> Tensor d (r - 1) a
indexSlice dim i t
  | dim < tensorRank t = undefined
  | otherwise =
    error $
      "invalid dim=" <> show dim <> "; tensor is " <> show (tensorRank t) <> "-dimensional"

withForeignPtr' :: forall a b m. PrimBase m => ForeignPtr a -> (Ptr a -> m b) -> m b
withForeignPtr' fp action = unsafeIOToPrim $ withForeignPtr fp action'
  where
    action' p = unsafePrimToIO (action p)

productPrimArray :: (Prim a, Num a) => PrimArray a -> a
productPrimArray = foldlPrimArray' (*) 1
{-# INLINE productPrimArray #-}

tensorLength :: Tensor device rank a -> Int
tensorLength t = productPrimArray (tensorShape t)
{-# INLINE tensorLength #-}

tensorExtent :: (HasCallStack, KnownNat r) => Int -> Tensor d r a -> Int
tensorExtent dim t
  | dim < tensorRank t = indexPrimArray (tensorShape t) dim
  | otherwise =
    error $
      "index out of bounds: " <> show dim <> "; tensor is " <> show (tensorRank t)
        <> "-dimensional"

tensorRank :: forall device rank a. KnownNat rank => Tensor device rank a -> Int
tensorRank _ = natToInt @rank
{-# INLINE tensorRank #-}

tensorToList1D :: forall a. Prim a => Tensor 'CPU 1 a -> [a]
tensorToList1D t = runST $
  withForeignPtr' (tensorData t) $ \p ->
    go p [] (stride * (extent - 1))
  where
    !stride = indexPrimArray (tensorStrides t) 0
    !extent = indexPrimArray (tensorShape t) 0
    go :: PrimMonad m => Ptr a -> [a] -> Int -> m [a]
    go !p acc !i
      | i >= 0 = do
        !x <- P.readOffPtr p i
        go p (x : acc) (i - stride)
      | otherwise = touch t >> pure acc

tensorFromList1D :: forall a. (HasCallStack, Prim a) => Int -> [a] -> Tensor 'CPU 1 a
tensorFromList1D n xs = runST $ do
  t <- newTensor [n]
  withForeignPtr' (tensorData t) $ \p -> do
    let go !i []
          | i == n = pure ()
          | otherwise = error $ "list is shorter than expected"
        go !i (y : ys)
          | i < n = P.writeOffPtr p i y >> go (i + 1) ys
          | otherwise = error "list is longer than expected"
    go 0 xs
  pure t

rowMajorStride :: forall device rank a. KnownNat rank => Tensor device rank a -> Maybe Int
rowMajorStride t = go 0
  where
    !rank = natToInt @rank
    !shape = tensorShape t
    !stride = tensorStrides t
    go !i
      | i < rank - 1 =
        let !n₂ = indexPrimArray shape (i + 1)
            !s₁ = indexPrimArray stride i
            !s₂ = indexPrimArray stride (i + 1)
         in if s₁ `div` s₂ == n₂ then go (i + 1) else Nothing
      | rank /= 0 = Just $ indexPrimArray stride (rank - 1)
      | otherwise = Just 1

flatten :: KnownNat rank => Tensor device rank a -> Maybe (Tensor device 1 a)
flatten t = case rowMajorStride t of
  Just stride ->
    Just $
      t
        { tensorShape = primArrayFromListN 1 [tensorLength t],
          tensorStrides = primArrayFromListN 1 [stride]
        }
  Nothing -> Nothing

reshape ::
  forall r2 a r1 d.
  (KnownNat r1, KnownNat r2) =>
  [Int] ->
  Tensor d r1 a ->
  Maybe (Tensor d r2 a)
reshape shape t
  | product shape == tensorLength t = case rowMajorStride t of
    Just s₀ ->
      let !shape' = primArrayFromListN (natToInt @r2) shape
          !strides' = rowMajorStrides s₀ shape'
       in Just $ t {tensorShape = shape', tensorStrides = strides'}
    Nothing -> Nothing
  | otherwise = Nothing

splitAt' :: HasCallStack => Int -> [a] -> ([a], [a])
splitAt' = go
  where
    go 0 [] = ([], [])
    go 1 (x : xs) = ([x], xs)
    go m (x : xs) = (x : xs', xs'')
      where
        (xs', xs'') = go (m - 1) xs
    go _ [] = error "wrong list length"

chunksOf :: Int -> [a] -> [[a]]
chunksOf n = go
  where
    go [] = []
    go xs@(_ : _) = let (ys, zs) = splitAt' n xs in ys : go zs

listShape2D :: [[a]] -> [Int]
listShape2D [] = [0, 0]
listShape2D xs@(x : _) = [length xs, length x]
{-# INLINE listShape2D #-}

listShape3D :: [[[a]]] -> [Int]
listShape3D [] = [0, 0, 0]
listShape3D xs@(x : _) = length xs : listShape2D x
{-# INLINE listShape3D #-}

listShape4D :: [[[[a]]]] -> [Int]
listShape4D [] = [0, 0, 0, 0]
listShape4D xs@(x : _) = length xs : listShape3D x
{-# INLINE listShape4D #-}

--
-- listShape5D :: [[[[a]]]] -> [Int]
-- listShape5D [] = [0, 0, 0, 0, 0]
-- listShape5D xs@(x : _) = length xs : listShape4D x

instance Prim a => GHC.IsList (Tensor 'CPU 1 a) where
  type Item (Tensor 'CPU 1 a) = a
  toList = tensorToList1D
  fromList xs = tensorFromList1D (length xs) xs
  fromListN n xs = tensorFromList1D n xs

instance Prim a => GHC.IsList (Tensor 'CPU 2 a) where
  type Item (Tensor 'CPU 2 a) = [a]
  toList t = case flatten t of
    Just t' -> chunksOf d₁ (tensorToList1D t')
      where
        !d₁ = indexPrimArray (tensorShape t) 1
    Nothing -> error "toList does not work with strided tensors"
  fromList xs = fromJust . reshape @2 [d₀, d₁] $ tensorFromList1D (d₀ * d₁) (mconcat xs)
    where
      [!d₀, !d₁] = listShape2D xs

instance Prim a => GHC.IsList (Tensor 'CPU 3 a) where
  type Item (Tensor 'CPU 3 a) = [[a]]
  toList t = case flatten t of
    Just t' -> chunksOf d₁ . chunksOf d₂ $ tensorToList1D t'
      where
        !d₁ = indexPrimArray (tensorShape t) 1
        !d₂ = indexPrimArray (tensorShape t) 2
    Nothing -> error "toList does not work with strided tensors"
  fromList xs =
    fromJust . reshape @3 [d₀, d₁, d₂] $
      tensorFromList1D (d₀ * d₁ * d₂) (mconcat . mconcat $ xs)
    where
      [!d₀, !d₁, !d₂] = listShape3D xs

instance Prim a => GHC.IsList (Tensor 'CPU 4 a) where
  type Item (Tensor 'CPU 4 a) = [[[a]]]
  toList t = case flatten t of
    Just t' -> chunksOf d₁ . chunksOf d₂ . chunksOf d₃ $ tensorToList1D t'
      where
        !d₁ = indexPrimArray (tensorShape t) 1
        !d₂ = indexPrimArray (tensorShape t) 2
        !d₃ = indexPrimArray (tensorShape t) 3
    Nothing -> error "toList does not work with strided tensors"
  fromList xs =
    fromJust . reshape @4 [d₀, d₁, d₂, d₃] $
      tensorFromList1D (d₀ * d₁ * d₂ * d₃) (mconcat . mconcat . mconcat $ xs)
    where
      [!d₀, !d₁, !d₂, !d₃] = listShape4D xs

-- indexSlow :: (Show a, Storable a) => Tensor a -> [Int] -> a
-- indexSlow t i
--   | any (< 0) i = error $ "invalid index: " <> show i
--   | otherwise =
--     let !linearIndex = sum $ zipWith (*) (primArrayToList (tensorStrides t)) i
--      in unsafePerformIO $ do
--           print (i, linearIndex)
--           withForeignPtr (tensorData t) $ \p -> peekElemOff p linearIndex
-- {-# NOINLINE indexSlow #-}

-- tensorToList :: (Show a, Storable a) => Tensor a -> [a]
-- tensorToList t =
--   fmap (indexSlow t)
--     . sequence
--     . fmap (\n -> [0 .. n - 1])
--     . primArrayToList
--     . tensorShape
--     $ t

reversePrimArray :: Prim a => PrimArray a -> PrimArray a
reversePrimArray xs = runST $ do
  let !n = sizeofPrimArray xs
  ys <- newPrimArray n
  let go !i
        | i < n = writePrimArray ys i (indexPrimArray xs (n - 1 - i)) >> go (i + 1)
        | otherwise = pure ()
  go 0
  unsafeFreezePrimArray ys
{-# INLINE reversePrimArray #-}

rowMajorStrides :: Int -> PrimArray Int -> PrimArray Int
rowMajorStrides s₀ shape
  | n > 1 = runST $ do
    strides <- newPrimArray n
    let go !s !i
          | i >= 0 =
            let !s' = s * indexPrimArray shape (i + 1)
             in writePrimArray strides i s' >> go s' (i - 1)
          | otherwise = pure ()
    writePrimArray strides (n - 1) s₀
    go s₀ (n - 2)
    unsafeFreezePrimArray strides
  | otherwise = replicatePrimArray n s₀
  where
    n = sizeofPrimArray shape

-- rowMajorStrides :: Integral a => [a] -> [a]
-- rowMajorStrides = drop 1 . scanr (*) 1

newCpuBuffer :: forall a m. (Prim a, PrimBase m) => Int -> m (ForeignPtr a)
newCpuBuffer extent =
  unsafeIOToPrim $
    mallocPlainForeignPtrAlignedBytes (extent * sizeOf element) (max 64 (alignment element))
  where
    element = undefined :: a

newTensor :: forall r a m. (KnownNat r, Prim a, PrimBase m) => [Int] -> m (Tensor 'CPU r a)
newTensor shape = do
  let !rank = natToInt @r
      !shape' = fromListN rank shape
  p <- newCpuBuffer (productPrimArray shape')
  pure $ Tensor p shape' (rowMajorStrides 1 shape')

tensorSimilar :: (Prim a, PrimBase m) => Tensor 'CPU rank a -> m (Tensor 'CPU rank a)
tensorSimilar t = do
  p <- newCpuBuffer (tensorLength t)
  pure $ Tensor p (tensorShape t) (rowMajorStrides 1 (tensorShape t))

instance (IsDLDataType a, PrimBase m) => IsDLTensor m (Tensor 'CPU r a) where
  withDLTensor t action = withForeignPtr' (tensorData t) $ \p ->
    viaContiguousBuffer
      p
      (primArrayToList $ tensorShape t)
      (primArrayToList $ tensorStrides t)
      action

instance IsHalideType a => IsHalideBuffer (Tensor 'CPU rank a) where
  withHalideBuffer t action = withForeignPtr (tensorData t) $ \p ->
    bufferFromPtrShapeStrides
      p
      (primArrayToList . reversePrimArray $ tensorShape t)
      (primArrayToList . reversePrimArray $ tensorStrides t)
      action

generate1 :: Prim a => Int -> (Int -> a) -> Tensor 'CPU 1 a
generate1 n f = runST $ do
  t <- newTensor [n]
  withForeignPtr' (tensorData t) $ \p ->
    let go !i
          | i < n = P.writeOffPtr p i (f i) >> go (i + 1)
          | otherwise = pure ()
     in go 0
  pure t

loopM :: Monad m => i -> (i -> Bool) -> (i -> i) -> (i -> m ()) -> m ()
loopM i₀ cond inc action = go i₀
  where
    go !i
      | cond i = do () <- action i; go (inc i)
      | otherwise = pure ()
{-# INLINE loopM #-}

makeLoopsC :: forall m. Monad m => PrimArray Int -> PrimArray Int -> (Int -> m ()) -> m ()
makeLoopsC shape strides action
  | rank == 0 = pure ()
  | otherwise = loop 0 0
  where
    !rank = sizeofPrimArray shape
    loop !dim !i₀ =
      loopM 0 (< count) (+ 1) $ \k ->
        let !i = i₀ + k * stride in action' i
      where
        !count = indexPrimArray shape dim
        !stride = indexPrimArray strides dim
        action'
          | dim < rank - 1 = loop (dim + 1)
          | otherwise = action

tensorMap ::
  (KnownNat r, IsTblisType a, Prim a, Prim b) =>
  (a -> b) ->
  Tensor 'CPU r a ->
  Tensor 'CPU r b
tensorMap f a
  | isContiguous a = runST $ do
    b <- newTensor (primArrayToList (tensorShape a))
    withForeignPtr' (tensorData a) $ \src ->
      withForeignPtr' (tensorData b) $ \dest ->
        let !n = tensorLength a
            go !i
              | i < n = P.writeOffPtr dest i (f (P.indexOffPtr src i)) >> go (i + 1)
              | otherwise = pure ()
         in go 0
    pure b
  | otherwise = tensorMap f (toContiguous a)

isContiguous :: KnownNat r => Tensor d r a -> Bool
isContiguous t = case rowMajorStride t of
  Just s -> s == 1
  Nothing -> False

toContiguous ::
  forall r a.
  (HasCallStack, KnownNat r, IsTblisType a, Prim a) =>
  Tensor 'CPU r a ->
  Tensor 'CPU r a
toContiguous a = runST $ do
  b <- tensorSimilar a
  r <- withDLTensor a $ \aBuffer -> withDLTensor b $ \bBuffer ->
    runExceptT $ tblisAdd (1 :: a) aBuffer index (0 :: a) bBuffer index
  case r of
    Right () -> pure b
    Left (TblisError e) -> error e
  where
    index = take (tensorRank a) ['a' ..]

tensorMult ::
  forall s r₃ r₁ r₂ a.
  (HasCallStack, KnownSymbol s, KnownNat r₁, KnownNat r₂, KnownNat r₃, Prim a, IsTblisType a) =>
  Tensor 'CPU r₁ a ->
  Tensor 'CPU r₂ a ->
  Tensor 'CPU r₃ a
tensorMult a b = runST $ do
  trace (show indexA) $ pure ()
  trace (show indexB) $ pure ()
  trace (show indexC) $ pure ()
  c <- newTensor shapeC
  r <-
    withDLTensor a $ \aBuffer ->
      withDLTensor b $ \bBuffer ->
        withDLTensor c $ \cBuffer ->
          runExceptT $
            tblisMult
              (1 :: a)
              aBuffer
              (toString indexA)
              (1 :: a)
              bBuffer
              (toString indexB)
              (0 :: a)
              cBuffer
              (toString indexC)
  case r of
    Right () -> pure c
    Left (TblisError e) -> error e
  where
    expr = symbolToText @s
    (indexA, expr') = T.break (== ',') expr
    (indexB, expr'') = T.breakOn "->" (T.drop 1 expr')
    indexC = T.drop 2 expr''
    getExtent !ch =
      case T.findIndex (== ch) indexA of
        Just i -> tensorExtent i a
        Nothing -> case T.findIndex (== ch) indexB of
          Just i -> tensorExtent i b
          Nothing -> error $ "invalid index expression: " <> expr
    shapeC = T.foldr (\ch ns -> getExtent ch : ns) [] indexC
