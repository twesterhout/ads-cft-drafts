module EinsteinEquations.Halide
  (
  )
where

import Foreign.Ptr (Ptr)

-- import Data.Bits

data HalideDimension = HalideDimension
  { halideDimensionMin :: !Int32,
    halideDimensionExtent :: !Int32,
    halideDimensionStride :: !Int32,
    halideDimensionFlags :: !Word32
  }
  deriving stock (Read, Show, Eq)

data HalideDeviceInterface

newtype HalideType = HalideType Word32
  deriving stock (Read, Show, Eq)

data HalideBuffer = HalideBuffer
  { halideBufferDevice :: !Word64,
    halideBufferDeviceInterface :: !(Ptr HalideDeviceInterface),
    halideBufferHost :: !(Ptr Word8),
    halideBufferFlags :: !Word64,
    halideBufferType :: !HalideType,
    halideBufferDimensions :: !Int32,
    halideBufferDim :: !(Ptr HalideDimension),
    halideBufferPadding :: !(Ptr ())
  }
  deriving stock (Show, Eq)
