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
module Metric.LineSearch
  ( LineSearchOptions,
    lineSearch,
  )
where

import ArrayFire (AFType, Array)
import qualified ArrayFire as AF
import Control.Monad.IO.Unlift
import Metric.IDR
import Prelude hiding (state)

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

newtype Jacobian m a = Jacobian {unJacobian :: Array a -> m (Array a)}

getInitialSlope :: Monad m => Function m a -> Jacobian m a -> Array a -> Array a -> m a
getInitialSlope func jac p y₀ = do
  w <- unJacobian jac p
  undefined

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

lineSearch :: a
lineSearch = undefined
