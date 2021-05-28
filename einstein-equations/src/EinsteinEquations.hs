{- |
Copyright: (c) 2021 Tom Westerhout
SPDX-License-Identifier: BSD-3-Clause
Maintainer: Tom Westerhout <14264576+twesterhout@users.noreply.github.com>

See README for more info
-}

module EinsteinEquations
       ( someFunc
       ) where

import Data.Primitive


someFunc :: IO ()
someFunc = putStrLn ("someFunc" :: String)

data Tensor (rank :: Nat) (a :: Type)

gdd ::
  a -> -- ^ x
  a -> -- ^ y
  a -> -- ^ z
  PrimArray a -> -- ^ other parameters (e.g. Q)
  Tensor 2 a
gdd = undefined
