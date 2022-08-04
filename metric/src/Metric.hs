{-# LANGUAGE OverloadedStrings #-}

-- |
-- Copyright: (c) 2022 Tom Westerhout
-- SPDX-License-Identifier: BSD-3-Clause
-- Maintainer: Tom Westerhout <14264576+twesterhout@users.noreply.github.com>
--
-- See README for more info
module Metric
  ( someFunc,
    gridPointsForPeriodic,
    gridPointsForBounded,
    differentiationMatrixPeriodic,
    differentiationMatrixBounded,
    differentiateX,
    differentiateY,
    differentiateZ,
  )
where

import ArrayFire (AFType, Array, scalar)
import qualified ArrayFire as AF

someFunc :: IO ()
someFunc = do
  putStrLn ("someFunc" :: String)

gridPointsForPeriodic :: Double -> Int -> Array Double
gridPointsForPeriodic period n
  | even n = AF.scalar (period / fromIntegral n) * (AF.iota @Double [n] [] + AF.scalar 1)

gridPointsForBounded :: Double -> Double -> Int -> Array Double
gridPointsForBounded a b n
  | n >= 2 = AF.scalar ((b + a) / 2) + AF.scalar ((b - a) / 2) * AF.cos (scale * js)
  | otherwise = error $ "invalid n: " <> show n
  where
    scale = AF.scalar (pi / fromIntegral (n - 1))
    js = AF.iota [n] []

differentiationMatrixPeriodic :: Int -> Array Double
differentiationMatrixPeriodic n
  | even n =
    AF.select
      isDiag
      (AF.scalar @Double 0)
      ( AF.scalar 0.5 * AF.cast (AF.pow (AF.scalar (-1)) δi)
          / AF.tan (scale * AF.cast δi)
      )
  | otherwise = error "currently only even n is supported"
  where
    scale = AF.scalar @Double (pi / fromIntegral n)
    rowIndices = AF.iota @Int [n] [1, n]
    colIndices = AF.transpose rowIndices False
    δi = rowIndices - colIndices
    isDiag = AF.eq rowIndices colIndices

differentiationMatrixBounded :: Double -> Double -> Int -> Array Double
differentiationMatrixBounded l u n = AF.select isDiag diag offDiag
  where
    isDiag = AF.identity [n, n]
    xᵢ = AF.tile (gridPointsForBounded l u n) [1, n]
    xⱼ = AF.transpose xᵢ False
    δx = xᵢ - xⱼ
    diag = flip AF.sum 1 $ AF.select isDiag (scalar 0) (scalar 1 / δx)
    a = flip AF.product 1 $ AF.select isDiag (scalar 1) δx
    aᵢ = AF.tile a [1, n]
    aⱼ = AF.transpose aᵢ False
    offDiag = aᵢ / (aⱼ * δx)

differentiateX ::
  (HasCallStack, AFType a) =>
  Array a ->
  Array a ->
  Array a
differentiateX d f
  | n₂ == 1 && n₃ == 1 && n₁ == m₀ = df
  | otherwise = error $ "incompatible dimensions: " <> show dDims <> " and " <> show fDims
  where
    dDims@(n₀, n₁, n₂, n₃) = AF.getDims d
    fDims@(m₀, m₁, m₂, m₃) = AF.getDims f
    f' = AF.moddims f [m₀, m₁ * m₂ * m₃]
    df' = AF.matmul d f' AF.None AF.None
    df = AF.moddims df' [n₀, m₁, m₂, m₃]

differentiateY ::
  (HasCallStack, AFType a) =>
  Array a ->
  Array a ->
  Array a
differentiateY d f
  | n₂ == 1 && n₃ == 1 && n₁ == m₁ = df
  | otherwise = error $ "incompatible dimensions: " <> show dDims <> " and " <> show fDims
  where
    dDims@(n₀, n₁, n₂, n₃) = AF.getDims d
    fDims@(m₀, m₁, m₂, m₃) = AF.getDims f
    f' = AF.moddims (AF.reorder f [1, 0, 2, 3]) [m₁, m₀ * m₂ * m₃]
    df' = AF.matmul d f' AF.None AF.None
    df = AF.reorder (AF.moddims df' [n₀, m₀, m₂, m₃]) [1, 0, 2, 3]

differentiateZ ::
  (HasCallStack, AFType a) =>
  Array a ->
  Array a ->
  Array a
differentiateZ d f
  | n₂ == 1 && n₃ == 1 && n₁ == m₂ = df
  | otherwise = error $ "incompatible dimensions: " <> show dDims <> " and " <> show fDims
  where
    dDims@(n₀, n₁, n₂, n₃) = AF.getDims d
    fDims@(m₀, m₁, m₂, m₃) = AF.getDims f
    f' = AF.moddims (AF.reorder f [2, 1, 0, 3]) [m₂, m₁ * m₀ * m₃]
    df' = AF.matmul d f' AF.None AF.None
    df = AF.reorder (AF.moddims df' [n₀, m₁, m₀, m₃]) [2, 1, 0, 3]
