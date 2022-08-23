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
module Metric.Jacobian
  ( computeExplicitJacobian,
    approximateJacobianVectorProduct,
  )
where

import ArrayFire (AFType, Array)
import qualified ArrayFire as AF
import Control.Monad.IO.Unlift
import Metric.IDR (epsilon, norm)
import Prelude hiding (state)

computeExplicitJacobian ::
  forall a m.
  (AFType a, RealFloat a, MonadUnliftIO m) =>
  -- | Function @f@
  (Array a -> m (Array a)) ->
  -- | Point around which to linearize @x₀@
  Array a ->
  -- | Function value at point @x₀@, i.e. @y₀ = f(x₀)@
  Array a ->
  m (Array a)
computeExplicitJacobian f x₀ y₀ = do
  let h = realToFrac (sqrt (1 + norm x₀)) * sqrt epsilon
      v = AF.join 0 (AF.constant [1] h) (AF.constant [n - 1] 0)
      n = AF.getElements x₀
      column !i = do
        y <- f (x₀ + AF.shift v i 0 0 0)
        pure $ (y - y₀) / AF.scalar h
  if n == 0
    then pure $ AF.constant [] 0
    else do
      columns <- sequence [column i | i <- [0 .. n - 1]]
      case columns of
        (c : cs) -> pure $ foldl' (AF.join 1) c cs

approximateJacobianVectorProduct ::
  forall a m.
  (AFType a, RealFloat a, MonadUnliftIO m) =>
  -- | Function @f@
  (Array a -> m (Array a)) ->
  -- | Point around which to linearize @x₀@
  Array a ->
  -- | Function value at point @x₀@, i.e. @y₀ = f(x₀)@
  Array a ->
  -- | Jacobian-vector product
  (Array a -> m (Array a))
approximateJacobianVectorProduct f x₀ y₀ = \a -> do
  let h = AF.scalar $ numerator / realToFrac (norm a)
  y <- f (x₀ + h * a)
  pure $ (y - y₀) / h
  where
    numerator :: a
    numerator = realToFrac (sqrt (1 + norm x₀)) * sqrt epsilon
