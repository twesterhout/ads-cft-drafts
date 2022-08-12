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
module Metric.IDR where

import ArrayFire (AFType, Array, Seq (..), scalar)
import qualified ArrayFire as AF
import Data.Foldable (foldl1)
import System.IO.Unsafe (unsafePerformIO)

data IDRParams a = IDRParams
  { idrParamsMaxIters :: !Int,
    idrParamsTol :: !Double
  }

testMatrix1 :: (Array Double, Array Double, Array Double)
testMatrix1 = (m, c, x)
  where
    m =
      flip AF.transpose False $
        AF.matrix (4, 4) $
          [ [0.95261881, 0.03465117, 0.87306725, 0.61793314],
            [0.2121299, 0.12656886, 0.94639018, 0.78809681],
            [0.82831116, 0.66228212, 0.48960417, 0.4083788],
            [0.08738027, 0.63494534, 0.72405378, 0.50726822]
          ]
    c = AF.vector 4 [0.62064746, 0.62386813, 0.73991291, 0.35641101]
    x = AF.vector 4 [0.6186891104674163, 0.9996017127630863, 0.036766805575082606, 0.28460881338224864]

testState1 :: (Array Double -> Array Double, IDRState Double)
testState1 = (apply, IDRState 1 (norm r) ω p x r m f u g)
  where
    (matrix, b, x) = testMatrix1
    apply v = AF.matmul matrix v AF.None AF.None
    s = 2
    ω = 1
    p =
      AF.matrix
        (4, s)
        [ [0.020593093173895904, 0.009953819642658623, 0.9320396946592179, 0.02338647451286391],
          [0.9743651736396345, 0.6473466046531297, 0.6410814296533377, 0.002716544045039182]
        ]
    r = b - apply x
    m = AF.identity [s, s]
    f = computeF p r
    u = AF.constant [4, s] 0
    g = AF.constant [4, s] 0

data IDRState a = IDRState
  { idrIter :: !Int,
    idrNormR :: !Double,
    idrOmega :: !a,
    idrP :: !(Array a),
    idrX :: !(Array a),
    idrR :: !(Array a),
    idrM :: !(Array a),
    idrF :: !(Array a),
    idrU :: !(Array a),
    idrG :: !(Array a)
  }
  deriving stock (Show)

keepGoing :: IDRParams a -> IDRState a -> Bool
keepGoing params state =
  idrNormR state > idrParamsTol params
    && idrIter state < idrParamsMaxIters params

idrsStep ::
  forall a.
  (AFType a, Fractional a, Show a) =>
  (Array a -> Array a) ->
  Int ->
  IDRState a ->
  IDRState a
idrsStep apply !k (IDRState !i !rNorm !ω !p !x !r !m !f !u !g) =
  IDRState (i + 1) rNorm' ω p x' r' m' f' u' g'
  where
    (_, !s, _, _) = AF.getDims p
    !k' = fromIntegral k
    !c = computeC k m f
    (!v, !q) = computeVandQ k c g u
    (!uₖ, !gₖ) = computeUₖandGₖ k ω apply (r - v) q
    !αs = computeAlphas k p m gₖ
    (!uₖ', !gₖ') = orthogonalizeUₖandGₖ k αs uₖ gₖ u g
    !u' = AF.assign u [Seq 0 (-1) 1, Seq k' k' 1] uₖ'
    !g' = AF.assign g [Seq 0 (-1) 1, Seq k' k' 1] gₖ'
    !m' = AF.assign m [Seq k' (-1) 1, Seq k' k' 1] (computeNextColumnOfM k p gₖ')
    !β = AF.index f [Seq k' k' 1] / AF.index m' [Seq k' k' 1, Seq k' k' 1]
    !r' = r - β * gₖ'
    !x' = x + β * uₖ'
    !rNorm' = norm r'
    !f' =
      if k < s - 1
        then AF.assign f [Seq (k' + 1) (-1) 1] (AF.index (f - β * col m' k) [Seq (k' + 1) (-1) 1])
        else f

idrsIter ::
  forall a.
  (AFType a, Real a, Fractional a, Show a) =>
  (Array a -> Array a) ->
  IDRParams a ->
  IDRState a ->
  IDRState a
idrsIter apply params state₀ = go 0 $ state₀ {idrF = computeF (idrP state₀) (idrR state₀)}
  where
    (_, s, _, _) = AF.getDims (idrP state₀)
    finalize !state =
      let r = idrR state
          v = r
          q = apply v
          ω' = computeOmega q r
          r' = (idrR state) - scalar ω' * q
       in state
            { idrIter = idrIter state + 1,
              idrNormR = norm r',
              idrOmega = ω',
              idrX = (idrX state) + scalar ω' * v,
              idrR = r'
            }
    go !k !state
      | k < s =
        let !state' = idrsStep apply k state
         in if keepGoing params state' then go (k + 1) state' else finalize state'
      | otherwise = finalize state

idrs ::
  (AFType a, Real a, Fractional a, Show a) =>
  IDRParams a ->
  -- Matrix
  (Array a -> Array a) ->
  -- Initial vector
  IDRState a ->
  -- Solution
  IDRState a
idrs params apply state₀ = go state₀
  where
    go !state
      | keepGoing params state = go (idrsIter apply params state)
      | otherwise = state

computeC :: AFType a => Int -> Array a -> Array a -> Array a
computeC k m f =
  AF.solve
    (AF.index m [Seq k' (-1) 1, Seq k' (-1) 1])
    (AF.index f [Seq k' (-1) 1])
    AF.Lower
  where
    k' = fromIntegral k

-- # Now we have sufficient vectors in G_j to compute residual in G_j+1
-- # Note: r is already perpendicular to P so v = r
-- copyto!(V, R)

-- # Preconditioning
-- ldiv!(Pl, V)

-- mul!(Q, A, V)
-- om = omega(Q, R)
-- R .-= om .* Q
-- X .+= om .* V

-- normR = norm(R)
-- if smoothing
--     T_s .= R_s .- R

--     gamma = dot(R_s, T_s)/dot(T_s, T_s)

--     R_s .-= gamma .* T_s
--     X_s .-= gamma .* (X_s .- X)

--     normR = norm(R_s)
-- end
-- iter += 1
-- # nextiter!(log, mvps=1)
-- # push!(log, :resnorm, normR)

norm :: AFType a => Array a -> Double
norm v = AF.norm v AF.NormVector2 0 0

inner :: AFType a => Array a -> Array a -> a
inner a b = AF.getScalar $ AF.dot a b AF.None AF.None

computeOmega :: (AFType a, Real a, Fractional a) => Array a -> Array a -> a
computeOmega t s = realToFrac $ if ρ < θ then ω * θ / ρ else ω
  where
    θ = 1 / sqrt 2
    sNorm = norm s
    tNorm = norm t
    ts = realToFrac $ inner t s
    ρ = abs $ ts / (sNorm * tNorm)
    ω = ts / (tNorm * tNorm)

col :: (HasCallStack, AFType a) => Array a -> Int -> Array a
col arr i
  | AF.getNumDims arr <= 2 = AF.index arr [Seq 0 (-1) 1, Seq (fromIntegral i) (fromIntegral i) 1]
  | otherwise =
    error $
      "cannot get the "
        <> show i
        <> "'th column of an array of shape "
        <> show (AF.getDims arr)
        <> "; array is not a matrix"

computeF :: AFType a => Array a -> Array a -> Array a
computeF p r = AF.vector s [inner (p `col` i) r | i <- [0 .. s - 1]]
  where
    (_, !s, _, _) = AF.getDims p

computeVandQ :: (HasCallStack, AFType a, Num a, Show a) => Int -> Array a -> Array a -> Array a -> (Array a, Array a)
computeVandQ k c g u = case terms of
  ((c₀, g₀, u₀) : others) -> foldl' combine (c₀ * g₀, c₀ * u₀) others
  _ -> error $ "s = 0 which should never happen"
  where
    (!size, _, _, _) = AF.getDims c
    !terms =
      [ (AF.index c [Seq (fromIntegral i) (fromIntegral i) 1], col g (k + i), col u (k + i))
        | i <- [0 .. size - 1]
      ]
    combine (!v, !q) (cᵢ, gᵢ, uᵢ) = (v + cᵢ * gᵢ, q + cᵢ * uᵢ)

computeUₖandGₖ :: (AFType a, Num a) => Int -> a -> (Array a -> Array a) -> Array a -> Array a -> (Array a, Array a)
computeUₖandGₖ k ω apply v q = (u, g)
  where
    u = q + scalar ω * v
    g = apply u

computeAlphas :: (HasCallStack, AFType a, Fractional a) => Int -> Array a -> Array a -> Array a -> Array a
computeAlphas k p m gₖ
  | k == 0 = AF.constant [] (0 / 0)
  | otherwise = AF.moddims (overlaps / diag) [1, k, 1, 1]
  where
    overlaps = AF.vector k [inner (col p i) gₖ | i <- [0 .. k - 1]]
    diag = AF.index (AF.diagExtract m 0) [Seq 0 (fromIntegral (k - 1)) 1]

orthogonalizeUₖandGₖ ::
  (HasCallStack, AFType a, Num a) =>
  Int ->
  Array a ->
  Array a ->
  Array a ->
  Array a ->
  Array a ->
  (Array a, Array a)
orthogonalizeUₖandGₖ k α uₖ gₖ u g
  | k == 0 = (uₖ, gₖ)
  | otherwise = (uₖ', gₖ')
  where
    relevant arr = AF.index arr [Seq 0 (-1) 1, Seq 0 (fromIntegral k - 1) 1]
    uₖ' = uₖ - AF.sum (α * relevant u) 1
    gₖ' = gₖ - AF.sum (α * relevant g) 1

computeNextColumnOfM :: (HasCallStack, AFType a, Num a) => Int -> Array a -> Array a -> Array a
computeNextColumnOfM k p gₖ = AF.vector (s - k) [inner (col p i) gₖ | i <- [k .. s - 1]]
  where
    (_, !s, _, _) = AF.getDims p

generateP :: AFType a => Int -> Int -> IO (Array a)
generateP s n = AF.randu [n, s]

--   proc computeOmega(t, s) {
--       const angle = sqrt(2.0) / 2;
--       const ns = norm(s);
--       const nt = norm(t);
--       const ts = inner(t, s);
--       const rho = abs(ts / (nt * ns));
--       var om = ts / (nt * nt);
--       if rho < angle then
--           om = om * angle:om.type / rho;
--       return om;
--   }
