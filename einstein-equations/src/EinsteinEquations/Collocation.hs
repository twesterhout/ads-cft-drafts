module EinsteinEquations.Collocation
  ( differentiationMatrixPeriodic,
    differentiationMatrixBounded,
    gridPointsForPeriodic,
    gridPointsForBounded,
    differentiateX,
    differentiateY,
    differentiateZ,
  )
where

import Control.Foldl (FoldM (..))
import Control.Monad.ST
import Data.Bits (toIntegralSized)
import EinsteinEquations.Tensor
import Foreign.C.Types (CDouble (..), CInt (..))
import Foreign.Ptr (Ptr)
import Foreign.Storable
import Halide
import Numeric.TBLIS
import qualified System.IO.Unsafe
import Torch

-- foreign import ccall unsafe "differentiation_matrix_periodic"
--   differentiation_matrix_periodic :: Double -> CInt -> Ptr HalideBuffer -> IO ()

foreign import ccall unsafe "differentiation_matrix_bounded"
  differentiation_matrix_bounded :: Double -> Double -> CInt -> Ptr HalideBuffer -> IO ()

-- foreign import ccall unsafe "differentiate_x"
--   differentiate_x :: Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> IO ()

differentiationMatrixPeriodic :: HasCallStack => Int -> Tensor
differentiationMatrixPeriodic n
  | even n = generate2 (n, n) f
  | otherwise = error "currently only even n is supported"
  where
    !δx = 2 * pi / fromIntegral n
    f :: Int -> Int -> Double
    f !k !j
      | k /= j = 0.5 * (-1) ^^ (k - j) / Prelude.tan (fromIntegral (k - j) * δx / 2)
      | otherwise = 0

gridPointsForPeriodic :: HasCallStack => Double -> Int -> Tensor
gridPointsForPeriodic period n
  | even n = generate1 n (\ !i -> (period / fromIntegral n) * fromIntegral (i + 1))
  | otherwise = error "currently only even n is supported"

differentiationMatrixBounded :: HasCallStack => Double -> Double -> Int -> Tensor
differentiationMatrixBounded lower upper n =
  System.IO.Unsafe.unsafePerformIO $ do
    withHalideBuffer t $ \buffer -> differentiation_matrix_bounded lower upper c_n buffer
    pure t
  where
    t = ones [n, n] (withDType Double defaultOpts)
    !c_n = case toIntegralSized n of
      Just n' -> n'
      Nothing -> error $ "overflow when casting " <> show n <> " to int"

-- differentiationMatrixBounded :: HasCallStack => Double -> Double -> Int -> Tensor 'CPU 2 Double
-- differentiationMatrixBounded l u n = undefined
--   where
--     x :: Tensor 'CPU 1 Double
--     x = gridPointsForBounded l u n
--     productExcept :: forall a. Num a => Tensor 'CPU 1 a -> Index 1 -> a
--     productExcept t j = runST $ do
--       xj <- read t j
--       let step !z !k
--             | k == j = do xk <- read t k; pure $ z * (xk - xj)
--             | otherwise = pure z
--           initial = pure 1
--           extract = pure
--       blockFoldM (FoldM step initial extract) (tensorShape t)
--     a :: Tensor 'CPU 1 Double
--     a = generate (I0 n) (productExcept x)

gridPointsForBounded :: Double -> Double -> Int -> Tensor
gridPointsForBounded a b n
  | n >= 2 = generate1 n f
  | otherwise = error $ "invalid n: " <> show n
  where
    f !j = (a + b) / 2 + (a - b) / 2 * Prelude.cos (fromIntegral j * pi / fromIntegral (n - 1))

differentiateX ::
  HasCallStack =>
  Tensor ->
  Tensor ->
  Tensor
differentiateX matrix function = einsum "aj,jbcd->abcd" [matrix, function]

differentiateY ::
  HasCallStack =>
  Tensor ->
  Tensor ->
  Tensor
differentiateY matrix function = einsum "bj,ajcd->abcd" [matrix, function]

differentiateZ ::
  HasCallStack =>
  Tensor ->
  Tensor ->
  Tensor
differentiateZ matrix function = einsum "cj,acjd->abcd" [matrix, function]

-- System.IO.Unsafe.unsafePerformIO $ do
--   derivative <- tensorSimilar function
--   withHalideBuffer matrix $ \matrixBuffer ->
--     withHalideBuffer function $ \functionBuffer ->
--       withHalideBuffer derivative $ \derivativeBuffer ->
--         differentiate_x matrixBuffer functionBuffer derivativeBuffer
--   pure derivative
