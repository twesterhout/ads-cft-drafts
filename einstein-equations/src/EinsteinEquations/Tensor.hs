module EinsteinEquations.Tensor
  ( generate1,
    generate2,
    generateM1,
    generateM2,
    loopM,
    foldLoopM,
    linalgSolve,
    strides,
    norm,
  )
where

-- Tensor (..),
-- Device (..),
-- read,
-- flatten,
-- reshape,
-- tensorLength,
-- newTensor,
-- tensorSimilar,
-- blockFoldM,
-- reversePrimArray,
-- rowMajorStrides,
-- tensorIndices,
-- tensorElements,
-- generate,
-- generateM,
-- tensorMap,
-- tensorMult,
-- SupportedRank,
-- Index (..),

import Control.Foldl (FoldM (..))
import Control.Monad.Primitive (PrimBase, PrimMonad, touch, unsafeIOToPrim, unsafePrimToIO)
import Control.Monad.ST
import DLPack
import Data.Complex
import Data.Kind
import Data.Maybe (fromJust)
import Data.Primitive.PrimArray
import qualified Data.Primitive.Ptr as P
import Data.Primitive.Types (Prim, alignment, sizeOf)
import Data.Reflection
import qualified Data.Text as T
import Data.Type.Equality ((:~:) (..))
import qualified Foreign.ForeignPtr as ForeignPtr
import Foreign.Ptr (Ptr, castPtr, plusPtr)
import qualified GHC.Exts as GHC
import GHC.ForeignPtr (ForeignPtr (..), mallocPlainForeignPtrAlignedBytes, unsafeForeignPtrToPtr)
import qualified GHC.ForeignPtr as ForeignPtr
import qualified GHC.Show as GHC
import GHC.TypeLits
import Halide
import qualified ListT as L
import Numeric.TBLIS
import Streamly.Prelude (IsStream, MonadAsync, Serial, SerialT)
import qualified Streamly.Prelude as S
import System.IO.Unsafe
import Torch
import qualified Torch.Functional.Internal
import Torch.Internal.Cast (cast1)
import qualified Torch.Internal.Managed.Type.Tensor as ATen

natToInt :: forall n. KnownNat n => Int
natToInt = fromIntegral $ GHC.TypeLits.natVal (Proxy @n)
{-# INLINE natToInt #-}

symbolToText :: forall s. KnownSymbol s => Text
symbolToText = toText $ GHC.TypeLits.symbolVal (Proxy @s)

strides :: Tensor -> [Int]
strides t = unsafePerformIO $ (cast1 ATen.tensor_strides) t

linalgSolve ::
  -- | A
  Tensor ->
  -- | B
  Tensor ->
  -- | x = A⁻¹B
  Tensor
linalgSolve = Torch.Functional.Internal.linalg_solve

norm :: Tensor -> Float -> Tensor
norm = Torch.Functional.Internal.normAll

instance IsHalideBuffer Tensor where
  withHalideBuffer :: forall b. Tensor -> (Ptr HalideBuffer -> IO b) -> IO b
  withHalideBuffer t action = withTensor t $ \p ->
    let action' :: IsHalideType a => Ptr a -> IO b
        action' p' = bufferFromPtrShapeStrides p' (reverse $ shape t) (reverse $ strides t) action
     in case dtype t of
          Bool -> action' (castPtr p :: Ptr Bool)
          UInt8 -> action' (castPtr p :: Ptr Word8)
          Int8 -> action' (castPtr p :: Ptr Int8)
          Int16 -> action' (castPtr p :: Ptr Int16)
          Int32 -> action' (castPtr p :: Ptr Int32)
          Int64 -> action' (castPtr p :: Ptr Int64)
          Float -> action' (castPtr p :: Ptr Float)
          Double -> action' (castPtr p :: Ptr Double)
          -- ComplexFloat -> action' (castPtr p :: Ptr (Complex Float))
          -- ComplexDouble -> action' (castPtr p :: Ptr (Complex Double))
          ty -> error $ show ty <> " is not supported"

loopM :: Monad m => i -> (i -> Bool) -> (i -> i) -> (i -> m ()) -> m ()
loopM i₀ cond inc action = go i₀
  where
    go !i
      | cond i = do () <- action i; go (inc i)
      | otherwise = pure ()
{-# INLINE loopM #-}

foldLoopM :: Monad m => i -> (i -> Bool) -> (i -> i) -> a -> (a -> i -> m a) -> m a
foldLoopM i₀ cond inc z₀ action = go z₀ i₀
  where
    go !z !i
      | cond i = do !z' <- action z i; go z' (inc i)
      | otherwise = pure z
{-# INLINE foldLoopM #-}

generateM1 :: forall a. (Prim a, Reifies a DType) => Int -> (Int -> IO a) -> IO Tensor
generateM1 d₀ f = do
  let t = zeros [d₀] $ withDType (reflect (Proxy @a)) defaultOpts
  withTensor t $ \p ->
    loopM 0 (< d₀) (+ 1) $ \ !i₀ ->
      f i₀ >>= P.writeOffPtr (castPtr p :: Ptr a) i₀
  pure t

generate1 :: (Prim a, Reifies a DType) => Int -> (Int -> a) -> Tensor
generate1 d₀ f = unsafePerformIO $ generateM1 d₀ (pure . f)

generateM2 :: forall a. (Prim a, Reifies a DType) => (Int, Int) -> (Int -> Int -> IO a) -> IO Tensor
generateM2 (d₀, d₁) f = do
  let t = zeros [d₀, d₁] $ withDType (reflect (Proxy @a)) defaultOpts
  withTensor t $ \p ->
    loopM 0 (< d₀) (+ 1) $ \ !i₀ ->
      loopM 0 (< d₁) (+ 1) $ \ !i₁ ->
        f i₀ i₁ >>= P.writeOffPtr (castPtr p :: Ptr a) (i₀ * d₁ + i₁)
  pure t

generate2 :: (Prim a, Reifies a DType) => (Int, Int) -> (Int -> Int -> a) -> Tensor
generate2 (d₀, d₁) f = unsafePerformIO $ generateM2 (d₀, d₁) (\i₀ i₁ -> pure (f i₀ i₁))

-- data Device = CPU
--   deriving stock (Read, Show, Eq)

-- data Tensor (device :: Device) (rank :: Nat) (a :: Type) = Tensor
--   { tensorData :: {-# UNPACK #-} !(ForeignPtr a),
--     tensorShape :: {-# UNPACK #-} !(PrimArray Int),
--     tensorStrides :: {-# UNPACK #-} !(PrimArray Int)
--   }
--   deriving stock (Show)

-- infixl 6 :.

-- data Index (r :: Nat) where
--   I0 :: {-# UNPACK #-} !Int -> Index 1
--   (:.) :: (KnownNat r, KnownNat (r + 1)) => {-# UNPACK #-} !(Index r) -> {-# UNPACK #-} !Int -> Index (r + 1)
--
-- select ::
--   (HasCallStack, KnownNat r, KnownNat (r - 1)) =>
--   -- | dim
--   Int ->
--   -- | index
--   Int ->
--   -- | input tensor
--   Tensor d r a ->
--   Tensor d (r - 1) a
-- select dim i t
--   | dim < tensorRank t = undefined
--   | otherwise =
--     error $
--       "invalid dim=" <> show dim <> "; tensor is " <> show (tensorRank t) <> "-dimensional"
--
-- touchForeignPtr :: forall a m. PrimMonad m => ForeignPtr a -> m ()
-- touchForeignPtr (ForeignPtr _ r) = touch r
-- {-# INLINE touchForeignPtr #-}
--
-- unsafeWithForeignPtr :: forall a b m. PrimMonad m => ForeignPtr a -> (Ptr a -> m b) -> m b
-- unsafeWithForeignPtr fp f = do
--   r <- f (unsafeForeignPtrToPtr fp)
--   touchForeignPtr fp
--   pure r
-- {-# INLINE unsafeWithForeignPtr #-}
--
-- unsafeWithTensor :: forall r a b m. PrimMonad m => Tensor 'CPU r a -> (Ptr a -> m b) -> m b
-- unsafeWithTensor t f = unsafeWithForeignPtr (tensorData t) f
-- {-# INLINE unsafeWithTensor #-}
--
-- -- withForeignPtr' :: forall a b m. PrimBase m => ForeignPtr a -> (Ptr a -> m b) -> m b
-- -- withForeignPtr' fp action = unsafeIOToPrim $ ForeignPtr.withForeignPtr fp action'
-- --   where
-- --     action' p = unsafePrimToIO (action p)
--
-- productPrimArray :: (Prim a, Num a) => PrimArray a -> a
-- productPrimArray = foldlPrimArray' (*) 1
-- {-# INLINE productPrimArray #-}
--
-- tensorLength :: Tensor device rank a -> Int
-- tensorLength t = productPrimArray (tensorShape t)
-- {-# INLINE tensorLength #-}
--
-- unsafeStride :: (HasCallStack, KnownNat r) => Int -> Tensor d r a -> Int
-- unsafeStride dim t = indexPrimArray (tensorStrides t) dim
-- {-# INLINE unsafeStride #-}
--
-- unsafeExtent :: (HasCallStack, KnownNat r) => Int -> Tensor d r a -> Int
-- unsafeExtent dim t = indexPrimArray (tensorShape t) dim
-- {-# INLINE unsafeExtent #-}
--
-- tensorExtent :: (HasCallStack, KnownNat r) => Int -> Tensor d r a -> Int
-- tensorExtent dim t
--   | dim < tensorRank t = indexPrimArray (tensorShape t) dim
--   | otherwise =
--     error $
--       "index out of bounds: " <> show dim <> "; tensor is " <> show (tensorRank t)
--         <> "-dimensional"
--
-- tensorShape1 :: Tensor d 1 a -> Index 1
-- tensorShape1 t = I0 (unsafeExtent 0 t)
-- {-# INLINE tensorShape1 #-}
--
-- tensorShape2 :: Tensor d 2 a -> Index 2
-- tensorShape2 t = I0 (unsafeExtent 0 t) :. (unsafeExtent 1 t)
-- {-# INLINE tensorShape2 #-}
--
-- tensorShape3 :: Tensor d 3 a -> Index 3
-- tensorShape3 t = I0 (unsafeExtent 0 t) :. (unsafeExtent 1 t) :. (unsafeExtent 2 t)
-- {-# INLINE tensorShape3 #-}
--
-- tensorShape4 :: Tensor d 4 a -> Index 4
-- tensorShape4 t = I0 (unsafeExtent 0 t) :. (unsafeExtent 1 t) :. (unsafeExtent 2 t) :. (unsafeExtent 3 t)
-- {-# INLINE tensorShape4 #-}
--
-- tensorShape5 :: Tensor d 5 a -> Index 5
-- tensorShape5 t = I0 (unsafeExtent 0 t) :. (unsafeExtent 1 t) :. (unsafeExtent 2 t) :. (unsafeExtent 3 t) :. (unsafeExtent 4 t)
-- {-# INLINE tensorShape5 #-}
--
-- tensorRank :: forall device rank a. KnownNat rank => Tensor device rank a -> Int
-- tensorRank _ = natToInt @rank
-- {-# INLINE tensorRank #-}
--
-- toLinearIndex :: forall r. (KnownNat r) => PrimArray Int -> Index r -> Int
-- toLinearIndex strides index = go (natToInt @r - 1) index
--   where
--     go :: (HasCallStack, KnownNat r') => Int -> Index r' -> Int
--     go !n (I0 i) = i * indexPrimArray strides n
--     go !n (rest :. i) = i * indexPrimArray strides n + go (n - 1) rest
--
-- read :: (KnownNat r, Prim a, PrimMonad m) => Tensor 'CPU r a -> Index r -> m a
-- read t i = unsafeWithTensor t $ \p -> P.readOffPtr p (toLinearIndex (tensorStrides t) i)
--
-- write :: (KnownNat r, Prim a, PrimMonad m) => Tensor 'CPU r a -> Index r -> a -> m ()
-- write t i x = unsafeWithTensor t $ \p -> P.writeOffPtr p (toLinearIndex (tensorStrides t) i) x
--
-- -- linearIndex1 :: Tensor d 1 a -> Int -> Int
-- -- linearIndex1 t i = i * unsafeStride 0 t
-- -- {-# INLINE linearIndex1 #-}
--
-- -- linearIndex2 :: Tensor d 2 a -> (Int, Int) -> Int
-- -- linearIndex2 t (i, j) = i * unsafeStride 0 t + j * unsafeStride 1 t
-- -- {-# INLINE linearIndex2 #-}
--
-- -- read1 :: (Prim a, PrimMonad m) => Tensor 'CPU 1 a -> Int -> m a
-- -- read1 t i = unsafeWithTensor t $ \p -> P.readOffPtr p (linearIndex1 t i)
-- -- {-# INLINE read1 #-}
--
-- -- read2 :: (Prim a, PrimMonad m) => Tensor 'CPU 2 a -> (Int, Int) -> m a
-- -- read2 t i = unsafeWithTensor t $ \p -> P.readOffPtr p (linearIndex2 t i)
-- -- {-# INLINE read2 #-}
--
-- -- write1 :: (Prim a, PrimMonad m) => Tensor 'CPU 1 a -> Int -> a -> m ()
-- -- write1 t i x = unsafeWithTensor t $ \p -> P.writeOffPtr p (linearIndex1 t i) x
-- -- {-# INLINE write1 #-}
--
-- -- write2 :: (Prim a, PrimMonad m) => Tensor 'CPU 2 a -> (Int, Int) -> a -> m ()
-- -- write2 t i x = unsafeWithTensor t $ \p -> P.writeOffPtr p (linearIndex2 t i) x
-- -- {-# INLINE write2 #-}
--
-- tensorToFlatList :: forall r a. (SupportedRank r, Prim a) => Tensor 'CPU r a -> [a]
-- tensorToFlatList t = reverse (tensorToReversedFlatList t)
--
-- tensorToReversedFlatList :: forall r a. (SupportedRank r, Prim a) => Tensor 'CPU r a -> [a]
-- tensorToReversedFlatList t =
--   runST $
--     tensorIndexFoldM (FoldM step initial pure) (tensorShape t)
--   where
--     step :: PrimMonad m => [a] -> Index r -> m [a]
--     step acc i = (: acc) <$> read t i
--     initial = pure []
--
-- tensorFromFlatList :: forall a. (Prim a) => Int -> [a] -> Tensor 'CPU 1 a
-- tensorFromFlatList n xs = runST $ do
--   t <- newTensor [n]
--   let step !i !x = write t (I0 i) x >> pure (i + 1)
--   _ <- foldlM step 0 xs
--   pure t
--
-- tensorToList1D :: forall a. Prim a => Tensor 'CPU 1 a -> [a]
-- tensorToList1D t = tensorToFlatList t {-runST $
--                                       unsafeWithForeignPtr (tensorData t) $ \p ->
--                                         go p [] (stride * (extent - 1))
--                                       where
--                                         !stride = indexPrimArray (tensorStrides t) 0
--                                         !extent = indexPrimArray (tensorShape t) 0
--                                         go :: PrimMonad m => Ptr a -> [a] -> Int -> m [a]
--                                         go !p acc !i
--                                           | i >= 0 = do
--                                             !x <- P.readOffPtr p i
--                                             go p (x : acc) (i - stride)
--                                           | otherwise = touch t >> pure acc-}
--
-- tensorFromList1D :: forall a. (HasCallStack, Prim a) => Int -> [a] -> Tensor 'CPU 1 a
-- tensorFromList1D n xs = tensorFromFlatList n xs {-runST $ do
--                                                 t <- newTensor [n]
--                                                 unsafeWithForeignPtr (tensorData t) $ \p -> do
--                                                   let go !i []
--                                                         | i == n = pure ()
--                                                         | otherwise = error $ "list is shorter than expected"
--                                                       go !i (y : ys)
--                                                         | i < n = P.writeOffPtr p i y >> go (i + 1) ys
--                                                         | otherwise = error "list is longer than expected"
--                                                   go 0 xs
--                                                 pure t-}
--
-- rowMajorStride :: forall device rank a. KnownNat rank => Tensor device rank a -> Maybe Int
-- rowMajorStride t
--   | rank == 0 = Just 1
--   | otherwise = go (indexPrimArray stride 0) 1
--   where
--     !rank = natToInt @rank
--     !shape = tensorShape t
--     !stride = tensorStrides t
--     go !s !i
--       | i < rank =
--         let !n = indexPrimArray shape i
--             !s' = indexPrimArray stride i
--          in if s `div` s' == n then go s' (i + 1) else Nothing
--       | otherwise = Just s
--       | rank /= 0 = Just $! indexPrimArray stride (rank - 1)
--       | otherwise = Just 1
--
-- flatten :: KnownNat rank => Tensor device rank a -> Maybe (Tensor device 1 a)
-- flatten t = case rowMajorStride t of
--   Just stride ->
--     Just $
--       t
--         { tensorShape = primArrayFromListN 1 [tensorLength t],
--           tensorStrides = primArrayFromListN 1 [stride]
--         }
--   Nothing -> Nothing
--
-- reshape ::
--   forall r2 a r1 d.
--   (KnownNat r1, KnownNat r2) =>
--   [Int] ->
--   Tensor d r1 a ->
--   Maybe (Tensor d r2 a)
-- reshape shape t
--   | product shape == tensorLength t = case rowMajorStride t of
--     Just s₀ ->
--       let !shape' = primArrayFromListN (natToInt @r2) shape
--           !strides' = rowMajorStrides s₀ shape'
--        in Just $ t {tensorShape = shape', tensorStrides = strides'}
--     Nothing -> Nothing
--   | otherwise = Nothing
--
-- splitAt' :: HasCallStack => Int -> [a] -> ([a], [a])
-- splitAt' = go
--   where
--     go 0 [] = ([], [])
--     go 1 (x : xs) = ([x], xs)
--     go m (x : xs) = (x : xs', xs'')
--       where
--         (xs', xs'') = go (m - 1) xs
--     go _ [] = error "wrong list length"
--
-- chunksOf :: Int -> [a] -> [[a]]
-- chunksOf n = go
--   where
--     go [] = []
--     go xs@(_ : _) = let (ys, zs) = splitAt' n xs in ys : go zs
--
-- listShape2D :: [[a]] -> [Int]
-- listShape2D [] = [0, 0]
-- listShape2D xs@(x : _) = [length xs, length x]
-- {-# INLINE listShape2D #-}
--
-- listShape3D :: [[[a]]] -> [Int]
-- listShape3D [] = [0, 0, 0]
-- listShape3D xs@(x : _) = length xs : listShape2D x
-- {-# INLINE listShape3D #-}
--
-- listShape4D :: [[[[a]]]] -> [Int]
-- listShape4D [] = [0, 0, 0, 0]
-- listShape4D xs@(x : _) = length xs : listShape3D x
-- {-# INLINE listShape4D #-}
--
-- --
-- -- listShape5D :: [[[[a]]]] -> [Int]
-- -- listShape5D [] = [0, 0, 0, 0, 0]
-- -- listShape5D xs@(x : _) = length xs : listShape4D x
--
-- instance Prim a => GHC.IsList (Tensor 'CPU 1 a) where
--   type Item (Tensor 'CPU 1 a) = a
--   toList = tensorToList1D
--   fromList xs = tensorFromList1D (length xs) xs
--   fromListN n xs = tensorFromList1D n xs
--
-- instance Prim a => GHC.IsList (Tensor 'CPU 2 a) where
--   type Item (Tensor 'CPU 2 a) = [a]
--   toList t = case flatten t of
--     Just t' -> chunksOf d₁ (tensorToList1D t')
--       where
--         !d₁ = indexPrimArray (tensorShape t) 1
--     Nothing -> error "toList does not work with strided tensors"
--   fromList xs = fromJust . reshape @2 [d₀, d₁] $ tensorFromList1D (d₀ * d₁) (mconcat xs)
--     where
--       [!d₀, !d₁] = listShape2D xs
--
-- instance Prim a => GHC.IsList (Tensor 'CPU 3 a) where
--   type Item (Tensor 'CPU 3 a) = [[a]]
--   toList t = case flatten t of
--     Just t' -> chunksOf d₁ . chunksOf d₂ $ tensorToList1D t'
--       where
--         !d₁ = indexPrimArray (tensorShape t) 1
--         !d₂ = indexPrimArray (tensorShape t) 2
--     Nothing -> error "toList does not work with strided tensors"
--   fromList xs =
--     fromJust . reshape @3 [d₀, d₁, d₂] $
--       tensorFromList1D (d₀ * d₁ * d₂) (mconcat . mconcat $ xs)
--     where
--       [!d₀, !d₁, !d₂] = listShape3D xs
--
-- instance Prim a => GHC.IsList (Tensor 'CPU 4 a) where
--   type Item (Tensor 'CPU 4 a) = [[[a]]]
--   toList t = case flatten t of
--     Just t' -> chunksOf d₁ . chunksOf d₂ . chunksOf d₃ $ tensorToList1D t'
--       where
--         !d₁ = indexPrimArray (tensorShape t) 1
--         !d₂ = indexPrimArray (tensorShape t) 2
--         !d₃ = indexPrimArray (tensorShape t) 3
--     Nothing -> error "toList does not work with strided tensors"
--   fromList xs =
--     fromJust . reshape @4 [d₀, d₁, d₂, d₃] $
--       tensorFromList1D (d₀ * d₁ * d₂ * d₃) (mconcat . mconcat . mconcat $ xs)
--     where
--       [!d₀, !d₁, !d₂, !d₃] = listShape4D xs
--
-- -- indexSlow :: (Show a, Storable a) => Tensor a -> [Int] -> a
-- -- indexSlow t i
-- --   | any (< 0) i = error $ "invalid index: " <> show i
-- --   | otherwise =
-- --     let !linearIndex = sum $ zipWith (*) (primArrayToList (tensorStrides t)) i
-- --      in unsafePerformIO $ do
-- --           print (i, linearIndex)
-- --           withForeignPtr (tensorData t) $ \p -> peekElemOff p linearIndex
-- -- {-# NOINLINE indexSlow #-}
--
-- -- tensorToList :: (Show a, Storable a) => Tensor a -> [a]
-- -- tensorToList t =
-- --   fmap (indexSlow t)
-- --     . sequence
-- --     . fmap (\n -> [0 .. n - 1])
-- --     . primArrayToList
-- --     . tensorShape
-- --     $ t
--
-- reversePrimArray :: Prim a => PrimArray a -> PrimArray a
-- reversePrimArray xs = runST $ do
--   let !n = sizeofPrimArray xs
--   ys <- newPrimArray n
--   let go !i
--         | i < n = writePrimArray ys i (indexPrimArray xs (n - 1 - i)) >> go (i + 1)
--         | otherwise = pure ()
--   go 0
--   unsafeFreezePrimArray ys
-- {-# INLINE reversePrimArray #-}
--
-- rowMajorStrides :: Int -> PrimArray Int -> PrimArray Int
-- rowMajorStrides s₀ shape
--   | n > 1 = runST $ do
--     strides <- newPrimArray n
--     let go !s !i
--           | i >= 0 =
--             let !s' = s * indexPrimArray shape (i + 1)
--              in writePrimArray strides i s' >> go s' (i - 1)
--           | otherwise = pure ()
--     writePrimArray strides (n - 1) s₀
--     go s₀ (n - 2)
--     unsafeFreezePrimArray strides
--   | otherwise = replicatePrimArray n s₀
--   where
--     n = sizeofPrimArray shape
--
-- -- rowMajorStrides :: Integral a => [a] -> [a]
-- -- rowMajorStrides = drop 1 . scanr (*) 1
--
-- newCpuBuffer :: forall a m. (Prim a, PrimMonad m) => Int -> m (ForeignPtr a)
-- newCpuBuffer extent =
--   unsafeIOToPrim $
--     mallocPlainForeignPtrAlignedBytes (extent * sizeOf element) (max 64 (alignment element))
--   where
--     element = undefined :: a
--
-- newTensor :: forall r a m. (KnownNat r, Prim a, PrimMonad m) => [Int] -> m (Tensor 'CPU r a)
-- newTensor shape = do
--   let !rank = natToInt @r
--       !shape' = fromListN rank shape
--   p <- newCpuBuffer (productPrimArray shape')
--   pure $ Tensor p shape' (rowMajorStrides 1 shape')
--
-- tensorSimilar :: (Prim a, PrimBase m) => Tensor 'CPU rank a -> m (Tensor 'CPU rank a)
-- tensorSimilar t = do
--   p <- newCpuBuffer (tensorLength t)
--   pure $ Tensor p (tensorShape t) (rowMajorStrides 1 (tensorShape t))
--
-- instance (IsDLDataType a, PrimBase m) => IsDLTensor m (Tensor 'CPU r a) where
--   withDLTensor t action = unsafeWithForeignPtr (tensorData t) $ \p ->
--     viaContiguousBuffer
--       p
--       (primArrayToList $ tensorShape t)
--       (primArrayToList $ tensorStrides t)
--       action
--
-- instance IsHalideType a => IsHalideBuffer (Tensor 'CPU rank a) where
--   withHalideBuffer t action = ForeignPtr.withForeignPtr (tensorData t) $ \p ->
--     bufferFromPtrShapeStrides
--       p
--       (primArrayToList . reversePrimArray $ tensorShape t)
--       (primArrayToList . reversePrimArray $ tensorStrides t)
--       action
--
-- -- generateM1 :: (Prim a, PrimMonad m) => Int -> (Int -> m a) -> m (Tensor 'CPU 1 a)
-- -- generateM1 d₀ f = do
-- --   t <- newTensor [d₀]
-- --   loopM 0 (< d₀) (+ 1) $ \ !i₀ ->
-- --     f i₀ >>= write t (I0 i₀)
-- --   pure t
--
-- -- generate1 :: (Prim a) => Int -> (Int -> a) -> Tensor 'CPU 1 a
-- -- generate1 d₀ f = runST $ generateM1 d₀ (pure . f)
--
-- -- generateM2 :: (Prim a, PrimMonad m) => (Int, Int) -> (Int -> Int -> m a) -> m (Tensor 'CPU 2 a)
-- -- generateM2 (d₀, d₁) f = do
-- --   t <- newTensor [d₀, d₁]
-- --   loopM 0 (< d₀) (+ 1) $ \ !i₀ ->
-- --     loopM 0 (< d₁) (+ 1) $ \ !i₁ ->
-- --       f i₀ i₁ >>= write t (I0 i₀ :. i₁)
-- --   pure t
--
-- -- generate2 :: Prim a => (Int, Int) -> (Int -> Int -> a) -> Tensor 'CPU 2 a
-- -- generate2 (d₀, d₁) f = runST $ generateM2 (d₀, d₁) (\i₀ i₁ -> pure (f i₀ i₁))
--
-- loopM :: Monad m => i -> (i -> Bool) -> (i -> i) -> (i -> m ()) -> m ()
-- loopM i₀ cond inc action = go i₀
--   where
--     go !i
--       | cond i = do () <- action i; go (inc i)
--       | otherwise = pure ()
-- {-# INLINE loopM #-}
--
-- foldLoopM :: Monad m => i -> (i -> Bool) -> (i -> i) -> a -> (a -> i -> m a) -> m a
-- foldLoopM i₀ cond inc z₀ action = go z₀ i₀
--   where
--     go !z !i
--       | cond i = do !z' <- action z i; go z' (inc i)
--       | otherwise = pure z
-- {-# INLINE foldLoopM #-}
--
-- _iterate0 :: forall r a d. (KnownNat r, 1 <= r) => Tensor d r a -> [Index 1]
-- _iterate0 t = I0 <$> enumFromTo 0 (tensorExtent 0 t - 1)
-- {-# INLINE _iterate0 #-}
--
-- _iterateN ::
--   forall r o a d.
--   (KnownNat r, KnownNat o, KnownNat (o + 1), o + 1 <= r) =>
--   Tensor d r a ->
--   [Index o] ->
--   [Index (o + 1)]
-- _iterateN t s = do
--   !i <- s
--   !j <- enumFromTo 0 (tensorExtent (natToInt @(o + 1)) t - 1)
--   pure $ i :. j
-- {-# INLINE _iterateN #-}
--
-- tensorIndexFoldM1 :: Monad m => FoldM m (Index 1) b -> PrimArray Int -> m b
-- tensorIndexFoldM1 (FoldM step initial extract) shape = do
--   !x0 <- initial
--   let !n0 = indexPrimArray shape 0
--   !r <- foldLoopM 0 (< n0) (+ 1) x0 $ \ !x1 !i0 ->
--     step x1 (I0 i0)
--   extract r
--
-- tensorIndexFoldM2 :: Monad m => FoldM m (Index 2) b -> PrimArray Int -> m b
-- tensorIndexFoldM2 (FoldM step initial extract) shape = do
--   !x0 <- initial
--   let !n0 = indexPrimArray shape 0
--       !n1 = indexPrimArray shape 1
--   !r <-
--     foldLoopM 0 (< n0) (+ 1) x0 $ \ !x1 !i0 ->
--       foldLoopM 0 (< n1) (+ 1) x1 $ \ !x2 !i1 ->
--         step x2 (I0 i0 :. i1)
--   extract r
--
-- tensorIndexFoldM3 :: Monad m => FoldM m (Index 3) b -> PrimArray Int -> m b
-- tensorIndexFoldM3 (FoldM step initial extract) shape = do
--   !x0 <- initial
--   let !n0 = indexPrimArray shape 0
--       !n1 = indexPrimArray shape 1
--       !n2 = indexPrimArray shape 2
--   !r <-
--     foldLoopM 0 (< n0) (+ 1) x0 $ \ !x1 !i0 ->
--       foldLoopM 0 (< n1) (+ 1) x1 $ \ !x2 !i1 ->
--         foldLoopM 0 (< n2) (+ 1) x2 $ \ !x3 !i2 ->
--           step x3 (I0 i0 :. i1 :. i2)
--   extract r
--
-- tensorIndexFoldM4 :: Monad m => FoldM m (Index 4) b -> PrimArray Int -> m b
-- tensorIndexFoldM4 (FoldM step initial extract) shape = do
--   !x0 <- initial
--   let !n0 = indexPrimArray shape 0
--       !n1 = indexPrimArray shape 1
--       !n2 = indexPrimArray shape 2
--       !n3 = indexPrimArray shape 3
--   !r <-
--     foldLoopM 0 (< n0) (+ 1) x0 $ \ !x1 !i0 ->
--       foldLoopM 0 (< n1) (+ 1) x1 $ \ !x2 !i1 ->
--         foldLoopM 0 (< n2) (+ 1) x2 $ \ !x3 !i2 ->
--           foldLoopM 0 (< n2) (+ 1) x3 $ \ !x4 !i3 ->
--             step x4 (I0 i0 :. i1 :. i2 :. i3)
--   extract r
--
-- blockFoldM ::
--   forall m r b.
--   (HasCallStack, Monad m, SupportedRank r) =>
--   FoldM m (Index r) b ->
--   PrimArray Int ->
--   m b
-- blockFoldM f shape
--   | natToInt @r == sizeofPrimArray shape = tensorIndexFoldM f shape
--   | otherwise = error "rank mismatch"
--
-- indexToList :: KnownNat r => Index r -> [Int]
-- indexToList i = go i []
--   where
--     go :: KnownNat k => Index k -> [Int] -> [Int]
--     go (ks :. k) acc = go ks (k : acc)
--     go (I0 k) acc = (k : acc)
--
-- generateM ::
--   (SupportedRank r, Prim a, PrimMonad m) =>
--   Index r ->
--   (Index r -> m a) ->
--   m (Tensor 'CPU r a)
-- generateM shape f = do
--   t <- newTensor (indexToList shape)
--   let initial = pure ()
--       extract _ = pure ()
--       step () !i = do !x <- f i; write t i x
--   blockFoldM (FoldM step initial extract) (tensorShape t)
--   pure t
--
-- generate ::
--   (SupportedRank r, Prim a) =>
--   Index r ->
--   (Index r -> a) ->
--   Tensor 'CPU r a
-- generate shape f = runST $ generateM shape (pure . f)
--
-- class (KnownNat r, KnownNat (r - 1)) => SupportedRank r where
--   tensorIndexFoldM :: Monad m => FoldM m (Index r) b -> PrimArray Int -> m b
--   tensorIndices :: Tensor d r a -> [Index r]
--
-- instance SupportedRank 1 where
--   tensorIndexFoldM = tensorIndexFoldM1
--   tensorIndices = _iterate0
--
-- instance SupportedRank 2 where
--   tensorIndexFoldM = tensorIndexFoldM2
--   tensorIndices t = _iterateN t $ _iterate0 t
--
-- instance SupportedRank 3 where
--   tensorIndexFoldM = tensorIndexFoldM3
--   tensorIndices t = _iterateN t $ _iterateN t $ _iterate0 t
--
-- instance SupportedRank 4 where
--   tensorIndexFoldM = tensorIndexFoldM4
--   tensorIndices t = _iterateN t $ _iterateN t $ _iterateN t $ _iterate0 t
--
-- tensorElements :: (SupportedRank r, Prim a, PrimMonad m) => Tensor 'CPU r a -> L.ListT m a
-- tensorElements t = L.traverse (read t) $ L.fromFoldable (tensorIndices t)
--
-- -- equalWith ::
-- --   (SupportedRank r, Prim a) =>
-- --   (a -> a -> Bool) ->
-- --   Tensor 'CPU r a ->
-- --   Tensor 'CPU r a ->
-- --   Bool
-- -- equalWith eq a b
-- --   | tensorShape a == tensorShape b =
-- --     unsafePerformIO
-- --       $! S.all id
-- --       $ S.zipWith eq (tensorElements a) (tensorElements b)
-- --   | otherwise = False
--
-- tensorIndices2 :: (IsStream t, Monad m, Monad (t m)) => Tensor d 2 a -> t m (Index 2)
-- tensorIndices2 t = do
--   i₀ <- S.enumerateFromTo 0 (tensorExtent 0 t)
--   i₁ <- S.enumerateFromTo 0 (tensorExtent 1 t)
--   pure $ I0 i₀ :. i₁
--
-- makeLoopsC :: forall m. Monad m => PrimArray Int -> PrimArray Int -> (Int -> m ()) -> m ()
-- makeLoopsC shape strides action
--   | rank == 0 = pure ()
--   | otherwise = loop 0 0
--   where
--     !rank = sizeofPrimArray shape
--     loop !dim !i₀ =
--       loopM 0 (< count) (+ 1) $ \k ->
--         let !i = i₀ + k * stride in action' i
--       where
--         !count = indexPrimArray shape dim
--         !stride = indexPrimArray strides dim
--         action'
--           | dim < rank - 1 = loop (dim + 1)
--           | otherwise = action
--
-- tensorMap ::
--   (KnownNat r, IsTblisType a, Prim a, Prim b) =>
--   (a -> b) ->
--   Tensor 'CPU r a ->
--   Tensor 'CPU r b
-- tensorMap f a
--   | isContiguous a = runST $ do
--     b <- newTensor (primArrayToList (tensorShape a))
--     unsafeWithForeignPtr (tensorData a) $ \src ->
--       unsafeWithForeignPtr (tensorData b) $ \dest ->
--         let !n = tensorLength a
--             go !i
--               | i < n = P.writeOffPtr dest i (f (P.indexOffPtr src i)) >> go (i + 1)
--               | otherwise = pure ()
--          in go 0
--     pure b
--   | otherwise = tensorMap f (toContiguous a)
--
-- isContiguous :: KnownNat r => Tensor d r a -> Bool
-- isContiguous t = case rowMajorStride t of
--   Just s -> s == 1
--   Nothing -> False
--
-- toContiguous ::
--   forall r a.
--   (HasCallStack, KnownNat r, IsTblisType a, Prim a) =>
--   Tensor 'CPU r a ->
--   Tensor 'CPU r a
-- toContiguous a = runST $ do
--   b <- tensorSimilar a
--   r <- withDLTensor a $ \aBuffer -> withDLTensor b $ \bBuffer ->
--     runExceptT $ tblisAdd (1 :: a) aBuffer index (0 :: a) bBuffer index
--   case r of
--     Right () -> pure b
--     Left (TblisError e) -> error e
--   where
--     index = take (tensorRank a) ['a' ..]
--
-- tensorMult ::
--   forall s r₃ r₁ r₂ a.
--   (HasCallStack, KnownSymbol s, KnownNat r₁, KnownNat r₂, KnownNat r₃, Prim a, IsTblisType a) =>
--   Tensor 'CPU r₁ a ->
--   Tensor 'CPU r₂ a ->
--   Tensor 'CPU r₃ a
-- tensorMult a b = runST $ do
--   -- trace (show indexA) $ pure ()
--   -- trace (show indexB) $ pure ()
--   -- trace (show indexC) $ pure ()
--   c <- newTensor shapeC
--   r <-
--     withDLTensor a $ \aBuffer ->
--       withDLTensor b $ \bBuffer ->
--         withDLTensor c $ \cBuffer ->
--           runExceptT $
--             tblisMult
--               (1 :: a)
--               aBuffer
--               (toString indexA)
--               (1 :: a)
--               bBuffer
--               (toString indexB)
--               (0 :: a)
--               cBuffer
--               (toString indexC)
--   case r of
--     Right () -> pure c
--     Left (TblisError e) -> error e
--   where
--     expr = symbolToText @s
--     (indexA, expr') = T.break (== ',') expr
--     (indexB, expr'') = T.breakOn "->" (T.drop 1 expr')
--     indexC = T.drop 2 expr''
--     getExtent !ch =
--       case T.findIndex (== ch) indexA of
--         Just i -> tensorExtent i a
--         Nothing -> case T.findIndex (== ch) indexB of
--           Just i -> tensorExtent i b
--           Nothing -> error $ "invalid index expression: " <> expr
--     shapeC = T.foldr (\ch ns -> getExtent ch : ns) [] indexC
