{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE StandaloneDeriving #-}

-- |
-- Copyright: (c) 2022 Tom Westerhout
-- SPDX-License-Identifier: BSD-3-Clause
-- Maintainer: Tom Westerhout <14264576+twesterhout@users.noreply.github.com>
--
-- See README for more info
module Metric.LineSearch where

import ArrayFire (AFType, Array)
import qualified ArrayFire as AF
import Control.Monad.Except (MonadError (..))
import Control.Monad.IO.Unlift
import Metric.Jacobian (approximateJacobianVectorProduct)
import Prelude hiding (State, state)

data LineSearchOrder = LineSearchQuadratic | LineSearchCubic
  deriving (Show, Eq)

data LineSearchOptions = LineSearchOptions
  { lineSearchAlpha :: !Double,
    lineSearchDamping :: !Double,
    lineSearchMaxIters :: !Int,
    lineSearchMinStep :: !Double,
    lineSearchMaxStep :: !Double,
    lineSearchOrder :: !LineSearchOrder
  }

newtype Function m a = Function {unFunction :: Array a -> m (Array a)}

-- | Jacobian-vector product
newtype Jacobian m a = Jacobian {unJacobian :: Array a -> m (Array a)}

-- getInitialSlope :: forall m a. (Monad m, AFType a, Real a) => Jacobian m a -> Array a -> m a
-- getInitialSlope jac p = do
--   g <- unJacobian jac p
--   let initSlope :: a
--       initSlope = AF.getScalar (AF.dot g p AF.None AF.None)
--   unless (initSlope < 0) $
--     error $ "positive initial slope: " <> show (realToFrac initSlope :: Double)
--   pure initSlope

-- clampStepSize :: (AFType a, Fractional a) => LineSearchOptions -> Array a -> (Array a, a)
-- clampStepSize opts p
--   | pNorm > maxStep = (AF.scalar (realToFrac (maxStep / pNorm)) * p, realToFrac maxStep)
--   where
--     pNorm = AF.norm p AF.NormVector2 0 0
--     maxStep = lineSearchMaxStep opts

getInitialState ::
  Monad m =>
  (Double -> m Double) ->
  Double ->
  LineSearchT m State
getInitialState func = go 0
  where
    go !i !λ = do
      maxIters <- reader lineSearchMaxIters
      if i >= maxIters
        then stopEarly $ LineSearchExceededMaxIter λ
        else do
          y <- lift $ func λ
          if isNaN y || isInfinite y
            then go (i + 1) (λ / 2)
            else pure $ State (i + 1) ((λ, y) :| [])

getNextStep :: Double -> Double -> NonEmpty (Double, Double) -> Double
getNextStep y₀ initSlope ((λ, y) :| []) =
  let λ' = -0.5 * initSlope / (y - y₀ - initSlope)
   in max λ' (0.1 * λ)
getNextStep y₀ initSlope ((λ, y) :| (λprev, yprev) : _) =
  let v₁ = y - y₀ - λ * initSlope
      v₂ = yprev - y₀ - λprev * initSlope
      m₁₁ = 1 / (λ * λ)
      m₁₂ = -1 / (λprev * λprev)
      m₂₁ = -λprev / (λ * λ)
      m₂₂ = λ / (λprev * λprev)
      a = (m₁₁ * v₁ + m₁₂ * v₂) / (λ - λprev)
      b = (m₂₁ * v₁ + m₂₂ * v₂) / (λ - λprev)
      d = max (b * b - 3 * a * initSlope) 0
      λtemp = if abs a < 1.0e-8 then -initSlope / (2 * b) else (-b + sqrt d) / (3 * a)
   in max (min λtemp (0.5 * λ)) (0.1 * λ)

data State = State
  { stateIter :: !Int,
    stateHistory :: NonEmpty (Double, Double)
  }

data LineSearchResult
  = LineSearchConverged !Double !Double !Int
  | LineSearchExceededMaxIter !Double
  | LineSearchReachedMinStep !Double
  | LineSearchNaNOrInfinity !Double
  deriving stock (Show, Eq)

newtype LineSearchT m a = LineSearchT {unLineSearchT :: ReaderT LineSearchOptions (ExceptT LineSearchResult m) a}
  deriving newtype (Functor, Applicative, Monad, MonadIO)

deriving newtype instance Monad m => MonadReader LineSearchOptions (LineSearchT m)

instance MonadTrans LineSearchT where
  lift = LineSearchT . lift . lift

stopEarly :: Monad m => LineSearchResult -> LineSearchT m a
stopEarly r = LineSearchT . lift . ExceptT . pure $ Left r

runLineSearchT :: Monad m => LineSearchT m () -> LineSearchOptions -> m LineSearchResult
runLineSearchT action opts = do
  r <- runExceptT $ runReaderT (unLineSearchT action) opts
  case r of
    Left e -> pure e
    Right () -> error "this should never happen"

searchLoop :: Monad m => (Double -> m Double) -> Double -> Double -> LineSearchT m ()
searchLoop func y₀ initSlope = do
  s <- getInitialState func 1
  go s
  where
    go s@(State i ((λ, y) :| history)) = do
      (α, maxIters, λₘᵢₙ) <- reader (\r -> (lineSearchAlpha r, lineSearchMaxIters r, lineSearchMinStep r))
      when (y < y₀ + α * λ * initSlope) . stopEarly $ LineSearchConverged λ y i
      when (i >= maxIters) . stopEarly $ LineSearchExceededMaxIter λ
      when (λ <= λₘᵢₙ) . stopEarly $ LineSearchReachedMinStep y
      let λ' = getNextStep y₀ initSlope (stateHistory s)
      y' <- lift $ func λ'
      go $ State (i + 1) ((λ', y') :| (λ, y) : history)

clampInitSlope :: Double -> Double
clampInitSlope initSlope
  | initSlope > 0 = -initSlope
  | initSlope == 0 = -1
  | otherwise = initSlope

defaultOptions :: LineSearchOptions
defaultOptions =
  LineSearchOptions
    { lineSearchAlpha = 1.0e-4,
      lineSearchDamping = 1.0,
      lineSearchMaxIters = 40,
      lineSearchMinStep = 1.0e-12,
      lineSearchMaxStep = 1.0,
      lineSearchOrder = LineSearchCubic
    }

lineSearch :: Monad m => LineSearchOptions -> (Double -> m Double) -> Double -> Double -> m LineSearchResult
lineSearch opts func y₀ initSlope = runLineSearchT (searchLoop func y₀ (clampInitSlope initSlope)) opts

lineSearchForNorm ::
  (MonadUnliftIO m, AFType a, RealFloat a) =>
  LineSearchOptions ->
  -- | Objective function @f@
  (Array a -> m (Array a)) ->
  -- | Initial point @x₀@
  Array a ->
  -- | Initial function value @y₀@
  Maybe (Array a) ->
  -- | Approximate search direction
  Array a ->
  m LineSearchResult
lineSearchForNorm opts func x₀ maybeValue p = do
  y₀ <- case maybeValue of
    Just y₀ -> pure y₀
    Nothing -> func x₀
  g <- approximateJacobianVectorProduct func x₀ y₀ p
  let !objective₀ =
        let !yNorm = AF.norm y₀ AF.NormVector2 0 0
         in yNorm * yNorm
      !initSlope = AF.getScalar (AF.dot y₀ g AF.None AF.None)
      scale = if initSlope > 0 then (-1) else 1 :: Double
      objective !λ = do
        y <- func (x₀ + AF.scalar (realToFrac (scale * λ)) * p)
        let !yNorm = AF.norm y AF.NormVector2 0 0
        liftIO . putStrLn $ "objective " <> show λ <> " = " <> show (yNorm * yNorm)
        pure $ yNorm * yNorm
  liftIO . putStrLn $ "objective₀ = " <> show objective₀
  liftIO . putStrLn $ "initSlope = " <> show initSlope
  r <- lineSearch opts objective objective₀ initSlope
  case r of
    LineSearchConverged λ y i -> pure $ LineSearchConverged (scale * λ) y i
    _ -> pure r
