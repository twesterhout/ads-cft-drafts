-- |
-- Copyright: (c) 2022 Tom Westerhout
-- SPDX-License-Identifier: BSD-3-Clause
-- Maintainer: Tom Westerhout <14264576+twesterhout@users.noreply.github.com>
--
-- See README for more info
module Metric
  ( someFunc,
  )
where

import qualified ArrayFire as A

someFunc :: IO ()
someFunc = do
  putStrLn ("someFunc" :: String)

-- print $ A.matrix @Double (2, 2) [[1, 2], [3, 4]]
