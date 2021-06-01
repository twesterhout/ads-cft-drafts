module EinsteinEquations.Tensor
  ( Tensor (..),
    Device (..),
    newTensor,
    tensorToList,
  )
where

import Data.Primitive.PrimArray
import Foreign.ForeignPtr
import Foreign.Ptr (plusPtr)
import Foreign.Storable
import GHC.ForeignPtr (mallocPlainForeignPtrAlignedBytes)
import Halide
import System.IO.Unsafe

data Device = CPU
  deriving stock (Read, Show, Eq)

data Tensor a = Tensor
  { tensorData :: {-# UNPACK #-} !(ForeignPtr a),
    tensorShape :: {-# UNPACK #-} !(PrimArray Int),
    tensorStrides :: {-# UNPACK #-} !(PrimArray Int),
    tensorDevice :: !Device
  }
  deriving stock (Show)

indexSlow :: (Show a, Storable a) => Tensor a -> [Int] -> a
indexSlow t i
  | any (< 0) i = error $ "invalid index: " <> show i
  | otherwise =
    let !linearIndex = sum $ zipWith (*) (primArrayToList (tensorStrides t)) i
     in unsafePerformIO $ do
          print (i, linearIndex)
          withForeignPtr (tensorData t) $ \p -> peekElemOff p linearIndex
{-# NOINLINE indexSlow #-}

tensorToList :: (Show a, Storable a) => Tensor a -> [a]
tensorToList t =
  fmap (indexSlow t)
    . sequence
    . fmap (\n -> [0 .. n - 1])
    . primArrayToList
    . tensorShape
    $ t

rowMajorStrides :: Integral a => [a] -> [a]
rowMajorStrides = drop 1 . scanr (*) 1

newTensor :: forall a. Storable a => Device -> [Int] -> IO (Tensor a)
newTensor device@CPU shape
  | any (< 0) shape = error $ "invalid shape: " <> show shape
  | otherwise = do
    let x = x :: a
    p <- mallocPlainForeignPtrAlignedBytes (product shape * sizeOf x) (max 64 (alignment x))
    return $
      Tensor
        p
        (fromList shape)
        (fromList . rowMajorStrides $ shape)
        device

instance IsHalideType a => IsHalideBuffer (Tensor a) where
  withHalideBuffer t action = withForeignPtr (tensorData t) $ \p ->
    bufferFromPtrShapeStrides
      p
      (primArrayToList $ tensorShape t)
      (primArrayToList $ tensorStrides t)
      action
