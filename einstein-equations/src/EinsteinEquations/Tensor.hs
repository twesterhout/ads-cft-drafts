module EinsteinEquations.Tensor
  ( Tensor (..),
    newTensor,
  )
where

import Data.Primitive.PrimArray
import Foreign.Storable
import GHC.ForeignPtr

data Device = CPU
  deriving stock (Read, Show, Eq)

data Tensor a = Tensor
  { tensorData :: {-# UNPACK #-} !(ForeignPtr a),
    tensorShape :: {-# UNPACK #-} !(PrimArray Int64),
    tensorStrides :: {-# UNPACK #-} !(PrimArray Int64),
    tensorDevice :: !Device
  }
  deriving stock (Show)

rowMajorStrides :: [Int64] -> [Int64]
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
        (fromList $ fromIntegral <$> shape)
        (fromList . rowMajorStrides $ fromIntegral <$> shape)
        device
