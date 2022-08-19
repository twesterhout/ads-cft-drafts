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
module Metric.IDR
  ( IDRParams (..),
    IDRResult (..),
    defaultIDRParams,
    epsilon,
    idrs,
    idrs',
  )
where

import ArrayFire (AFType, Array, Seq (..), scalar)
import qualified ArrayFire as AF
import Control.Monad.IO.Unlift
import Prelude hiding (state)

data IDRParams a = IDRParams
  { -- | Maximal number of iterations
    idrParamsMaxIters :: !Int,
    -- | Desired relative tolerance
    idrParamsTol :: !Double,
    -- | @s@ parameter in IDR(s), i.e. size of the basis
    idrParamsS :: !Int
  }
  deriving stock (Show, Eq, Generic)
  deriving anyclass (NFData)

epsilon :: forall a. RealFloat a => a
epsilon = encodeFloat 1 (1 - floatDigits (undefined :: a))

defaultIDRParams :: forall a. RealFloat a => IDRParams a
defaultIDRParams =
  IDRParams
    { idrParamsMaxIters = -1,
      idrParamsTol = realToFrac (sqrt (epsilon :: a)),
      idrParamsS = 8
    }

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

data IDRResult a = IDRResult
  { idrResultX :: !(Array a),
    idrResultResidual :: !Double,
    idrResultIter :: !Int,
    idrResultConverged :: !Bool
  }
  deriving stock (Show, Eq, Generic)

idrsStep ::
  forall a m.
  (AFType a, Fractional a, Monad m) =>
  (Array a -> m (Array a)) ->
  Int ->
  IDRState a ->
  m (IDRState a)
idrsStep apply k (IDRState i _ ω p x r m f u g) = do
  let (_, s, _, _) = AF.getDims p
      k' = fromIntegral k
  (uₖ', gₖ') <-
    orthogonalizeUandG k p m u g
      <$> computeUandG ω apply r (computeVandQ k (computeC k m f) g u)
  let u' = assignColumn u k uₖ'
      g' = assignColumn g k gₖ'
      m' = AF.assign m [Seq k' (-1) 1, Seq k' k' 1] (computeNextColumnOfM k p gₖ')
      β = AF.index f [Seq k' k' 1] / AF.index m' [Seq k' k' 1, Seq k' k' 1]
      r' = r - β * gₖ'
      x' = x + β * uₖ'
      rNorm' = norm r'
      f' =
        if k < s - 1
          then AF.assign f [Seq (k' + 1) (-1) 1] (AF.index (f - β * col m' k) [Seq (k' + 1) (-1) 1])
          else f
  pure $! IDRState (i + 1) rNorm' ω p x' r' m' f' u' g'

idrsIter ::
  forall a m.
  (AFType a, Real a, Fractional a, Monad m) =>
  IDRParams a ->
  (Array a -> m (Array a)) ->
  IDRState a ->
  m (IDRState a)
idrsIter params apply state₀ =
  go 0 $
    state₀ {idrF = computeF (idrP state₀) (idrR state₀)}
  where
    (_, s, _, _) = AF.getDims (idrP state₀)
    go :: Int -> IDRState a -> m (IDRState a)
    go !k !state
      | k < s = do
        !state' <- idrsStep apply k state
        if keepGoing params state'
          then go (k + 1) state'
          else finalize state'
      | otherwise = finalize state
    finalize :: IDRState a -> m (IDRState a)
    finalize !state = do
      let r = idrR state
          v = r
      q <- apply v
      let ω' = computeOmega q r
          r' = (idrR state) - scalar ω' * q
      pure $
        state
          { idrIter = idrIter state + 1,
            idrNormR = norm r',
            idrOmega = ω',
            idrX = (idrX state) + scalar ω' * v,
            idrR = r'
          }

mkState ::
  (HasCallStack, AFType a, Num a, MonadUnliftIO m) =>
  IDRParams a ->
  Array a ->
  (Array a -> m (Array a)) ->
  Array a ->
  Array a ->
  m (IDRParams a, IDRState a)
mkState params p apply b x
  | AF.getDims b == AF.getDims x && AF.getDims p == (n, s, 1, 1) = do
    r <- (b -) <$> apply x
    let ω = 1
        rNorm = norm r
        m = AF.identity [s, s]
        f = computeF p r
        u = AF.constant [n, s] 0
        g = AF.constant [n, s] 0
        maxIters =
          let k = idrParamsMaxIters params
           in if k == -1 then n else k
        tol = rNorm * idrParamsTol params
    liftIO . putStrLn $ "Initial rNorm " <> show rNorm
    pure (IDRParams maxIters tol s, IDRState 0 rNorm ω p x r m f u g)
  | otherwise =
    error $
      "incompatible dimensions: b has shape "
        <> show (AF.getDims b)
        <> ", x has shape "
        <> show (AF.getDims x)
        <> ", p has shape "
        <> show (AF.getDims p)
        <> ", and s = "
        <> show s
  where
    s = idrParamsS params
    n = case AF.getDims x of
      (size, 1, 1, 1) -> size
      dims -> error $ "x has wrong shape: " <> show dims <> ", but expected a vector"

mkResult :: IDRParams a -> IDRState a -> IDRResult a
mkResult params state = IDRResult (idrX state) (idrNormR state) (idrIter state) converged
  where
    converged = idrNormR state < idrParamsTol params

idrs ::
  (MonadUnliftIO m, AFType a, Real a, Fractional a) =>
  IDRParams a ->
  -- | A
  (Array a -> m (Array a)) ->
  -- | b
  Array a ->
  -- | x₀
  Array a ->
  -- Solution
  m (IDRResult a)
idrs params₀ apply b x₀ = do
  let (n, _, _, _) = AF.getDims b
      s = idrParamsS params₀
  p <- liftIO $ AF.randu [n, s] -- <*> pure 0.5
  idrs' params₀ p apply b x₀

idrs' ::
  (MonadUnliftIO m, AFType a, Real a, Fractional a) =>
  IDRParams a ->
  -- | p
  Array a ->
  -- | A
  (Array a -> m (Array a)) ->
  -- | b
  Array a ->
  -- | x₀
  Array a ->
  -- Solution
  m (IDRResult a)
idrs' params₀ p apply b x₀ = do
  (params, state₀) <- mkState params₀ p apply b x₀
  let go state
        | keepGoing params state = idrsIter params apply state >>= go
        | otherwise = pure state
  mkResult params <$> go state₀

computeC :: AFType a => Int -> Array a -> Array a -> Array a
computeC k m f =
  AF.solve
    (AF.index m [Seq k' (-1) 1, Seq k' (-1) 1])
    (AF.index f [Seq k' (-1) 1])
    AF.Lower
  where
    k' = fromIntegral k

keepGoing :: IDRParams a -> IDRState a -> Bool
keepGoing params state =
  idrNormR state > idrParamsTol params
    && idrIter state < idrParamsMaxIters params

assignColumn :: AFType a => Array a -> Int -> Array a -> Array a
assignColumn arr i rhs =
  AF.assign arr [Seq 0 (-1) 1, Seq (fromIntegral i) (fromIntegral i) 1] rhs

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

computeVandQ :: (AFType a, Num a) => Int -> Array a -> Array a -> Array a -> (Array a, Array a)
computeVandQ k c g u = case terms of
  ((c₀, g₀, u₀) : others) -> foldl' combine (c₀ * g₀, c₀ * u₀) others
  -- Actually, this branch is never taken because s > 0
  [] -> let (n, _, _, _) = AF.getDims g in (AF.constant [n] 0, AF.constant [n] 0)
  where
    (size, _, _, _) = AF.getDims c
    terms =
      [ (AF.index c [Seq (fromIntegral i) (fromIntegral i) 1], col g (k + i), col u (k + i))
        | i <- [0 .. size - 1]
      ]
    combine (!v, !q) (cᵢ, gᵢ, uᵢ) = (v + cᵢ * gᵢ, q + cᵢ * uᵢ)

computeUandG ::
  (Monad m, AFType a, Num a) =>
  a ->
  (Array a -> m (Array a)) ->
  Array a ->
  (Array a, Array a) ->
  m (Array a, Array a)
computeUandG ω apply r (v, q) = do
  let u = q + scalar ω * (r - v)
  g <- apply u
  pure (u, g)

-- computeAlphas :: (HasCallStack, AFType a, Fractional a) => Int -> Array a -> Array a -> Array a -> Array a
-- computeAlphas k p m gₖ
--   | k == 0 = AF.constant [] (0 / 0)
--   | otherwise = AF.moddims (overlaps / diag) [1, k, 1, 1]
--   where
--     overlaps = AF.vector k [inner (col p i) gₖ | i <- [0 .. k - 1]]
--     diag = AF.index (AF.diagExtract m 0) [Seq 0 (fromIntegral (k - 1)) 1]

orthogonalizeUandG ::
  (AFType a, Fractional a) =>
  Int ->
  Array a ->
  Array a ->
  Array a ->
  Array a ->
  (Array a, Array a) ->
  (Array a, Array a)
orthogonalizeUandG k p m u g = go 0
  where
    computeAlpha !j !gₖ = overlap / diag
      where
        j' = fromIntegral j
        overlap = inner (col p j) gₖ
        diag = AF.getScalar $ AF.index m [Seq j' j' 1, Seq j' j' 1]
    go !j (!uₖ, !gₖ)
      | j < k =
        let α = scalar (computeAlpha j gₖ)
            uₖ' = uₖ - α * col u j
            gₖ' = gₖ - α * col g j
         in go (j + 1) (uₖ', gₖ')
      | otherwise = (uₖ, gₖ)

computeNextColumnOfM :: AFType a => Int -> Array a -> Array a -> Array a
computeNextColumnOfM k p gₖ = AF.vector (s - k) [inner (col p i) gₖ | i <- [k .. s - 1]]
  where
    (_, !s, _, _) = AF.getDims p
