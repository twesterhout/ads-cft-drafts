module Main (main) where

import Data.Vector.Storable (Storable, Vector)
import qualified Data.Vector.Storable as V
import EinsteinEquations.NewtonRaphson
import Test.Hspec

roundTo :: RealFrac a => a -> Int -> a
roundTo x n = (fromInteger . round $ x * (10 ^ n)) / (10.0 ^^ n)

shouldApproxEqual ::
  (Storable a, Show a, RealFrac a) =>
  Vector a ->
  Vector a ->
  IO ()
shouldApproxEqual x y = V.zipWithM_ (\a b -> roundTo a n `shouldBe` roundTo b n) x y
  where
    n = 7

main :: IO ()
main = hspec $ do
  describe "nrm2" $ do
    it "computes L₂ norm of a vector" $ do
      let x :: (Storable a, Floating a) => Vector a
          x =
            fromList
              [ -0.4495389460471282,
                0.1507946864613462,
                -0.2582198045631925,
                0.19783960251972765,
                -0.08502606567957727,
                -0.027089038952040734,
                -0.45259419418301516,
                -0.34808112568663063,
                -0.2492391171933619,
                -0.03652065590737674
              ]
      nrm2 (x :: Vector Double) `shouldBe` 0.8532651379629458
      nrm2 (x :: Vector Float) `shouldBe` 0.8532651
  describe "solveDense" $ do
    it "solves linear systems of equations" $ do
      let a₁ =
            fromList
              [ [-0.12569571, -0.2457441, -0.0917917],
                [0.4542093, 0.23042252, -0.2675304],
                [0.26157499, -0.38497729, -0.2527988]
              ]
          b₁ = fromList [(0.38405567 :: Double), 0.05975013, 0.31492095]
          x₁ = solveDense a₁ b₁
      x₁ `shouldBe` (fromList [-1.172922576677295, -0.10262275813807724, -2.3030992629111187])
      let a₂ =
            fromList
              [ [0.13702166080474854, -0.19305512309074402],
                [0.036100149154663086, 0.4381260871887207]
              ]
          b₂ = fromList [-0.3771229088306427 :: Float, 0.23093461990356445]
          x₂ = solveDense a₂ b₂
      x₂ `shouldBe` (fromList [-1.8006048202514648, 0.675460159778595])
  describe "partialDerivative" $ do
    it "computes partial derivatives numerically" $ do
      let f :: Monad m => Vector Double -> m (Vector Double)
          f v =
            let x = v V.! 0
                y = v V.! 1
             in return . fromList $ [sin (x + y), cos (x - y)]
          x₀ = 0.1
          y₀ = 0.25
      df₀ <- partialDerivative f Nothing (fromList [x₀, y₀]) 0
      df₀ `shouldApproxEqual` (fromList [cos (x₀ + y₀), - sin (x₀ - y₀)])
      df₁ <- partialDerivative f Nothing (fromList [x₀, y₀]) 1
      df₁ `shouldApproxEqual` (fromList [cos (x₀ + y₀), sin (x₀ - y₀)])
