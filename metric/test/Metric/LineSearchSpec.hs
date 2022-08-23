{-# LANGUAGE ScopedTypeVariables #-}

module Metric.LineSearchSpec (spec) where

import ArrayFire (AFType, Array, scalar)
import qualified ArrayFire as AF
import Metric.LineSearch
import Test.Hspec

isConverged :: LineSearchResult -> Bool
isConverged (LineSearchConverged _ _ _) = True
isConverged _ = False

isClose :: (Num a, Ord a) => a -> a -> a -> a -> Bool
isClose a b rtol atol = abs (a - b) <= (atol + rtol * abs b)

spec :: Spec
spec = do
  describe "lineSearch" $ do
    it "optimizes y = x² starting from [-1, -1]" $ do
      let f x = AF.getScalar (AF.dot x x AF.None AF.None) :: Double
          x₀ = AF.vector 2 [-1.0, -1.0] :: Array Double
          y₀ = f x₀
          g = 2 * x₀
          p = -g
          initSlope = AF.getScalar (AF.dot p g AF.None AF.None)
          f' λ = pure $ f (x₀ + AF.scalar λ * p)
          opts = defaultOptions
      r <- lineSearch opts f' y₀ initSlope
      r `shouldSatisfy` isConverged
      case r of
        LineSearchConverged λ y i -> do
          λ `shouldSatisfy` (\x -> isClose x 0.5 1.0e-8 1.0e-10)
          y `shouldSatisfy` (\x -> isClose x 0.0 1.0e-8 1.0e-10)
          i `shouldBe` 2
        _ -> pure ()
    it "optimizes the Himmelblau function" $ do
      let himmelblau x = (x₁ * x₁ + x₂ - 11) ^ (2 :: Int) + (x₁ + x₂ * x₂ - 7) ^ (2 :: Int)
            where
              [x₁, x₂] = AF.toList x
          himmelblauGradient x = AF.vector 2 [g₁, g₂]
            where
              [x₁, x₂] = AF.toList x
              g₁ = 4.0 * x₁ ^ (3 :: Int) + 4.0 * x₁ * x₂ - 44.0 * x₁ + 2.0 * x₁ + 2.0 * x₂ ^ (2 :: Int) - 14.0
              g₂ = 2.0 * x₁ ^ (2 :: Int) + 2.0 * x₂ - 22.0 + 4.0 * x₁ * x₂ + 4.0 * x₂ ^ (3 :: Int) - 28.0 * x₂
          x₀ = AF.vector 2 [2.0, 2.0] :: Array Double
          y₀ = himmelblau x₀
          p = AF.vector 2 [42.0, 18.0]
          g = himmelblauGradient x₀
          f' λ = himmelblau (x₀ + AF.scalar λ * p)
          initSlope = AF.getScalar (AF.dot p g AF.None AF.None)
      r <- lineSearch defaultOptions (pure . f') y₀ initSlope
      r `shouldSatisfy` isConverged
      case r of
        LineSearchConverged λ y i -> do
          λ `shouldSatisfy` (\x -> isClose x 0.020545340808876406 1.0e-6 1.0e-8)
          y `shouldSatisfy` (\x -> isClose x (f' 0.020545340808876406) 1.0e-6 1.0e-8)
        _ -> pure ()
