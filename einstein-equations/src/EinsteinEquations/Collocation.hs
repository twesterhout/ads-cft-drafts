module EinsteinEquations.Collocation
  ( differentiationMatrixPeriodic,
    differentiationMatrixBounded,
  )
where

import Data.Bits (toIntegralSized)
import EinsteinEquations.Tensor
import Foreign.C.Types (CDouble (..), CInt (..))
import Foreign.Ptr (Ptr)
import Foreign.Storable
import Halide
import qualified System.IO.Unsafe

foreign import ccall unsafe "differentiation_matrix_periodic"
  differentiation_matrix_periodic :: Double -> CInt -> Ptr HalideBuffer -> IO ()

foreign import ccall unsafe "differentiation_matrix_bounded"
  differentiation_matrix_bounded :: Double -> Double -> CInt -> Ptr HalideBuffer -> IO ()

differentiationMatrixPeriodic :: HasCallStack => Double -> Int -> Tensor 'CPU 2 Double
differentiationMatrixPeriodic period n = System.IO.Unsafe.unsafePerformIO $ do
  t <- newTensor [n, n]
  let !c_n = case toIntegralSized n of
        Just n' -> n'
        Nothing -> error $ "overflow when casting " <> show n <> " to int"
  withHalideBuffer t $ \buffer -> differentiation_matrix_periodic period c_n buffer
  pure t

differentiationMatrixBounded :: HasCallStack => Double -> Double -> Int -> Tensor 'CPU 2 Double
differentiationMatrixBounded lower upper n = System.IO.Unsafe.unsafePerformIO $ do
  t <- newTensor [n, n]
  let !c_n = case toIntegralSized n of
        Just n' -> n'
        Nothing -> error $ "overflow when casting " <> show n <> " to int"
  withHalideBuffer t $ \buffer -> differentiation_matrix_bounded lower upper c_n buffer
  pure t
