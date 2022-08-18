{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE ScopedTypeVariables #-}

-- |
-- Copyright: (c) 2022 Tom Westerhout
-- SPDX-License-Identifier: BSD-3-Clause
-- Maintainer: Tom Westerhout <14264576+twesterhout@users.noreply.github.com>
--
-- See README for more info
module Metric.NewtonRaphson
  ( LineSearchOptions,
    RootOptions (..),
    RootResult (..),
    newtonRaphson,
  )
where

import ArrayFire (AFType, Array)
import qualified ArrayFire as AF
import Control.Monad.IO.Unlift
import Metric.IDR
import Prelude hiding (state)

data LineSearchOptions

data RootOptions = RootOptions
  { rootOptionsCriterion :: !(Double -> Double -> Bool),
    rootOptionsMaxIter :: !Int,
    rootOptionsLineSearch :: !(Maybe LineSearchOptions)
  }

data RootState a = RootState
  { rootStateCurrent :: !(Array a),
    rootStateValue :: !(Array a),
    rootStateResidual :: !Double
  }
  deriving stock (Show)

data RootResult a = RootResult
  { rootResultCurrent :: !(Array a),
    rootResultHistory :: [Double]
  }
  deriving stock (Show)

explicitJacobian ::
  forall a m.
  (AFType a, RealFloat a, MonadUnliftIO m) =>
  (Array a -> m (Array a)) ->
  Array a ->
  Array a ->
  m (Array a)
explicitJacobian f x₀ y₀ = do
  let jacobian v = do
        let h :: Array a
            h = AF.scalar $ 1.0e-8
        -- realToFrac (sqrt (1 + norm x₀)) * 1.0e-8 / realToFrac (norm a)
        -- print (AF.getScalar (AF.cast h :: Array Double) :: Double, norm a)
        y <- f (x₀ + h * v)
        -- print . AF.maxAll . AF.abs $ y - y₀
        pure $ AF.transpose ((y - y₀) / h) False
      n = AF.getElements x₀
      zeros i
        | i == 0 = AF.constant [] 0
        | otherwise = AF.constant [i] 0
      a i = AF.join 0 (AF.join 0 (zeros i) (AF.scalar 1)) (zeros (n - 1 - i))
  rows <- sequence [jacobian (a i) | i <- [0 .. n - 1]]
  -- print rows
  case rows of
    (r : rs) -> pure $ AF.transpose (foldl' (AF.join 0) r rs) False

invertApproximateJacobian ::
  forall a m.
  (AFType a, RealFloat a, MonadUnliftIO m) =>
  (Array a -> m (Array a)) ->
  Array a ->
  Array a ->
  m (Array a)
invertApproximateJacobian f x₀ y₀ = do
  j <- explicitJacobian f x₀ y₀
  let r = AF.solve j y₀ AF.None
  liftIO $ print r
  pure r

-- let jacobian a = do
--       let h :: Array a
--           h = AF.scalar $ realToFrac (sqrt (1 + norm x₀)) * 1.0e-8 / realToFrac (norm a)
--       -- print (AF.getScalar (AF.cast h :: Array Double) :: Double, norm a)
--       y <- f (x₀ + h * a)
--       -- print . AF.maxAll . AF.abs $ y - y₀
--       pure $ (y - y₀) / h
--     n = AF.getElements x₀
--     params = defaultIDRParams {idrParamsS = 128, idrParamsMaxIters = 400}
-- guess <- liftIO $ (0.01 *) <$> AF.randu [n]
-- (IDRResult δx r i isConverged) <- idrs params jacobian y₀ guess
-- unless isConverged $
--   error $ "IDR(s) failed to converge after " <> show i <> " iterations; residual norm is " <> show r
-- -- liftIO . putStrLn $ "IDR(s) converged to " <> show r <> " in " <> show i <> " iterations"
-- -- liftIO . putStrLn $ "J⁻¹(x₀) = " <> show (AF.toList (AF.cast δx :: Array Double))
-- pure δx

norm :: AFType a => Array a -> Double
norm v = AF.norm v AF.NormVector2 0 0

newtonRaphsonStep ::
  forall m a.
  (AFType a, RealFloat a, MonadUnliftIO m) =>
  (Array a -> m (Array a)) ->
  RootState a ->
  m (RootState a)
newtonRaphsonStep f (RootState x₀ y₀ _) = do
  δx <- invertApproximateJacobian f x₀ y₀
  let x = x₀ - δx
  y <- f x
  pure $ RootState x y (norm y)

newtonRaphson ::
  forall a m.
  (AFType a, RealFloat a, MonadUnliftIO m) =>
  RootOptions ->
  (Array a -> m (Array a)) ->
  Array a ->
  m (RootResult a)
newtonRaphson options f x₀ = do
  y₀ <- f x₀
  let iₘₐₓ = rootOptionsMaxIter options
      shouldStop = rootOptionsCriterion options r₀ . rootStateResidual
      go acc !i !s
        | i >= iₘₐₓ || shouldStop s = pure (rootStateCurrent s, acc)
        | otherwise = do
          s' <- newtonRaphsonStep f s
          go (rootStateResidual s' : acc) (i + 1) s'
      r₀ = norm y₀
  (solution, history) <- go [r₀] 0 $ RootState x₀ y₀ r₀
  pure $ RootResult solution (reverse history)
