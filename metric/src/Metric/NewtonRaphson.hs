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

invertApproximateJacobian ::
  (AFType a, RealFloat a, MonadUnliftIO m) =>
  (Array a -> m (Array a)) ->
  Array a ->
  Array a ->
  m (Array a)
invertApproximateJacobian f x₀ y₀ = do
  let jacobian a = do
        let h = AF.scalar $ realToFrac (sqrt (1 + norm x₀)) * sqrt epsilon / realToFrac (norm a)
        y <- f (x₀ + h * a)
        pure $ (y - y₀) / h
  (IDRResult δx r i isConverged) <- idrs defaultIDRParams jacobian y₀ x₀
  unless isConverged $
    error $ "IDR(s) failed to converge after " <> show i <> " iterations; residual norm is " <> show r
  -- liftIO . putStrLn $ "J⁻¹(x₀) = " <> show (AF.toList (AF.cast δx :: Array Double))
  pure δx

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
