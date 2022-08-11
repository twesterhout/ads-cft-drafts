module MetricSpec (spec) where

import ArrayFire (AFType, Array, scalar)
import qualified ArrayFire as AF
import Metric
import Test.Hspec

combineGrids :: AFType a => Array a -> Array a -> Array a -> (Array a, Array a, Array a)
combineGrids gridX gridY gridZ = (x, y, z)
  where
    (numX, _, _, _) = AF.getDims gridX
    (numY, _, _, _) = AF.getDims gridY
    (numZ, _, _, _) = AF.getDims gridZ
    !x = AF.tile gridX [1, numY, numZ]
    !y = AF.tile (AF.moddims gridY [1, numY]) [numX, 1, numZ]
    !z = AF.tile (AF.moddims gridZ [1, 1, numZ]) [numX, numY, 1]

allclose :: (AFType a, Num a) => Array a -> Array a -> a -> a -> Bool
allclose a b rtol atol =
  (== 1) . fst . AF.allTrueAll $
    AF.abs (a - b) `AF.le` (scalar atol + scalar rtol * AF.abs b)

spec :: Spec
spec = do
  describe "Metric" $ do
    it "constructs differentiation matrices for bounded domains" $ do
      let n = 6
          xs = gridPointsForBounded (0, 2) n
          dX = differentiationMatrixBounded (0, 2) n
          f = xs ^^ (3 :: Int) - 4 * xs
          dF = differentiateX dX f
          dFExpected = 3 * xs ^^ (2 :: Int) - 4
      allclose dF dFExpected 1.0e-10 (1.0e-12 :: Double) `shouldBe` True
    it "constructs differentiation matrices for periodic domains" $ do
      let n = 18
          k = 1.32
          t = scalar (2 * pi / k)
          xs = gridPointsForPeriodic k n
          dX = differentiationMatrixPeriodic k n
          f = AF.cos (t * 8 * xs - (pi / 3)) + AF.sin (t * xs)
          dF = differentiateX dX f
          dFExpected =
            -AF.sin (t * 8 * xs - (pi / 3)) * t * 8
              + AF.cos (t * xs) * t
      allclose dF dFExpected 1.0e-10 (1.0e-12 :: Double) `shouldBe` True
    it "computes ∂x, ∂y, ∂z" $ do
      let gridX = gridPointsForPeriodic (2 * pi) 6
          dX = differentiationMatrixPeriodic (2 * pi) 6
          gridY = gridPointsForBounded (1, 3) 40
          dY = differentiationMatrixBounded (1, 3) 40
          gridZ = gridPointsForBounded (0, 1) 40
          dZ = differentiationMatrixBounded (0, 1) 40
          (x, y, z) = combineGrids gridX gridY gridZ
          f = AF.sin (x - y * y + 2 * z) + z
          dF_xExpected = AF.cos (x - y * y + 2 * z)
          dF_yExpected = AF.cos (x - y * y + 2 * z) * (-2 * y)
          dF_zExpected = AF.cos (x - y * y + 2 * z) * 2 + 1
      -- print $ gridY
      -- print $ (differentiateY dY f)
      -- print $ dF_yExpected
      allclose (differentiateX dX f) dF_xExpected 1.0e-10 (1.0e-12 :: Double) `shouldBe` True
      allclose (differentiateY dY f) dF_yExpected 1.0e-10 (1.0e-12 :: Double) `shouldBe` True
      allclose (differentiateZ dZ f) dF_zExpected 1.0e-10 (1.0e-12 :: Double) `shouldBe` True
