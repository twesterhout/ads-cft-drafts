module EinsteinEquations.Tensor
  ( Tensor (..),
    Device (..),
    flatten,
    newTensor,
  )
where

import Control.Monad.Primitive (PrimBase, PrimMonad, touch, unsafeIOToPrim, unsafePrimToIO)
import Control.Monad.ST
import Data.Kind
import Data.Primitive.PrimArray
import qualified Data.Primitive.Ptr as P
import Data.Primitive.Types (Prim, alignment, sizeOf)
import Foreign.ForeignPtr
import Foreign.Ptr (Ptr, plusPtr)
-- import Foreign.Storable
import qualified GHC.Exts as GHC
import GHC.ForeignPtr (mallocPlainForeignPtrAlignedBytes)
import GHC.TypeLits
import Halide
import System.IO.Unsafe

natToInt :: forall n. KnownNat n => Int
natToInt = fromIntegral $ GHC.TypeLits.natVal (Proxy @n)
{-# INLINE natToInt #-}

data Device = CPU
  deriving stock (Read, Show, Eq)

data Tensor (device :: Device) (rank :: Nat) (a :: Type) = Tensor
  { tensorData :: {-# UNPACK #-} !(ForeignPtr a),
    tensorShape :: {-# UNPACK #-} !(PrimArray Int),
    tensorStrides :: {-# UNPACK #-} !(PrimArray Int)
  }
  deriving stock (Show)

withForeignPtr' :: forall a b m. PrimBase m => ForeignPtr a -> (Ptr a -> m b) -> m b
withForeignPtr' fp action = unsafeIOToPrim $ withForeignPtr fp action'
  where
    action' p = unsafePrimToIO (action p)

tensorLength :: Tensor device rank a -> Int
tensorLength t = foldlPrimArray' (*) 1 (tensorShape t)
{-# INLINE tensorLength #-}

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
tensorFromList1D n xs = unsafePerformIO $ do
  t <- newTensor [n]
  withForeignPtr (tensorData t) $ \p -> do
    let go !i []
          | i == n = pure ()
          | otherwise = error $ "list is shorter than expected"
        go !i (y : ys)
          | i < n = P.writeOffPtr p i y >> go (i + 1) ys
          | otherwise = error "list is longer than expected"
    go 0 xs
  pure t

rowMajorStride :: forall device rank a. KnownNat rank => Tensor device rank a -> Maybe Int
rowMajorStride t = let !f = go 0 in f
  where
    rank = natToInt @rank
    shape = tensorShape t
    stride = tensorStrides t
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

instance Prim a => GHC.IsList (Tensor 'CPU 1 a) where
  type Item (Tensor 'CPU 1 a) = a
  toList = tensorToList1D
  fromList xs = tensorFromList1D (length xs) xs
  fromListN n xs = tensorFromList1D n xs

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

instance Prim a => GHC.IsList (Tensor 'CPU 2 a) where
  type Item (Tensor 'CPU 2 a) = [a]
  toList t = case flatten t of
    Just t' -> chunksOf d₁ (tensorToList1D t')
      where
        !d₁ = indexPrimArray (tensorShape t) 1
    Nothing -> error "toList does not work with strided tensors"

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

rowMajorStrides :: Integral a => [a] -> [a]
rowMajorStrides = drop 1 . scanr (*) 1

newCpuBuffer :: forall a. Prim a => Int -> IO (ForeignPtr a)
newCpuBuffer extent =
  mallocPlainForeignPtrAlignedBytes (extent * sizeOf element) (max 64 (alignment element))
  where
    element = undefined :: a

newTensor :: forall rank a. (KnownNat rank, Prim a) => [Int] -> IO (Tensor 'CPU rank a)
newTensor shape = do
  let rank = natToInt @rank
  p <- newCpuBuffer (product shape)
  return $
    Tensor
      p
      (fromListN rank shape)
      (fromListN rank . rowMajorStrides $ shape)

instance IsHalideType a => IsHalideBuffer (Tensor 'CPU rank a) where
  withHalideBuffer t action = withForeignPtr (tensorData t) $ \p ->
    bufferFromPtrShapeStrides
      p
      (reverse . primArrayToList $ tensorShape t)
      (reverse . primArrayToList $ tensorStrides t)
      action
