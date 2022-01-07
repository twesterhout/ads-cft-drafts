module EinsteinEquations.Collocation
  ( differentiationMatrixPeriodic,
    differentiationMatrixBounded,
    gridPointsForPeriodic,
    gridPointsForBounded,
    differentiateX,
  )
where

import Data.Bits (toIntegralSized)
import EinsteinEquations.Tensor
import Foreign.C.Types (CDouble (..), CInt (..))
import Foreign.Ptr (Ptr)
import Foreign.Storable
import Halide
import Numeric.TBLIS
import qualified System.IO.Unsafe

foreign import ccall unsafe "differentiation_matrix_periodic"
  differentiation_matrix_periodic :: Double -> CInt -> Ptr HalideBuffer -> IO ()

foreign import ccall unsafe "differentiation_matrix_bounded"
  differentiation_matrix_bounded :: Double -> Double -> CInt -> Ptr HalideBuffer -> IO ()

foreign import ccall unsafe "differentiate_x"
  differentiate_x :: Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> IO ()

differentiationMatrixPeriodic :: HasCallStack => Double -> Int -> Tensor 'CPU 2 Double
differentiationMatrixPeriodic period n = System.IO.Unsafe.unsafePerformIO $ do
  t <- newTensor [n, n]
  let !c_n = case toIntegralSized n of
        Just n' -> n'
        Nothing -> error $ "overflow when casting " <> show n <> " to int"
  withHalideBuffer t $ \buffer -> differentiation_matrix_periodic period c_n buffer
  pure t

gridPointsForPeriodic :: HasCallStack => Double -> Int -> Tensor 'CPU 1 Double
gridPointsForPeriodic period n
  | even n = generate1 n (\i -> (period / fromIntegral n) * fromIntegral i)
  | otherwise = error "currently only even n is supported"

differentiationMatrixBounded :: HasCallStack => Double -> Double -> Int -> Tensor 'CPU 2 Double
differentiationMatrixBounded lower upper n = System.IO.Unsafe.unsafePerformIO $ do
  t <- newTensor [n, n]
  let !c_n = case toIntegralSized n of
        Just n' -> n'
        Nothing -> error $ "overflow when casting " <> show n <> " to int"
  withHalideBuffer t $ \buffer -> differentiation_matrix_bounded lower upper c_n buffer
  pure t

gridPointsForBounded :: Double -> Double -> Int -> Tensor 'CPU 1 Double
gridPointsForBounded = undefined

differentiateX ::
  HasCallStack =>
  Tensor 'CPU 2 Double ->
  Tensor 'CPU 4 Double ->
  Tensor 'CPU 4 Double
differentiateX matrix function = tensorMult @"aj,jbcd->abcd" matrix function

-- System.IO.Unsafe.unsafePerformIO $ do
--   derivative <- tensorSimilar function
--   withHalideBuffer matrix $ \matrixBuffer ->
--     withHalideBuffer function $ \functionBuffer ->
--       withHalideBuffer derivative $ \derivativeBuffer ->
--         differentiate_x matrixBuffer functionBuffer derivativeBuffer
--   pure derivative
