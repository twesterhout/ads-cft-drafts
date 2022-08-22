{-# LANGUAGE DeriveAnyClass #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}

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
    finiteDifferencesMatrix,
    differentiateX,
    differentiateY,
    differentiateZ,
    importFields,
    evalEquationsFromInput,
    importExpectedOutputs,
    reissnerNordstromFields,
    evalEquations,
    visualization,
  )
where

import ArrayFire (AFType, Array, Seq (..), scalar)
import qualified ArrayFire as AF
import Control.Monad.Primitive (touch)
import Control.Monad.Trans.Resource (MonadResource, register, release)
import Data.HDF5 (ArrayView' (..))
import qualified Data.HDF5 as H5
import Foreign.C.Types (CBool (..), CInt (..))
import Foreign.ForeignPtr (newForeignPtr_, withForeignPtr)
import Foreign.Ptr (Ptr)
import qualified GHC.ForeignPtr as GHC
import Graphics.Vega.VegaLite hiding (Axis, Fields, select)
import Language.Halide
import Metric.NewtonRaphson
import System.IO.Unsafe (unsafePerformIO)
import Prelude hiding (Seq)

someFunc :: IO ()
someFunc = do
  putStrLn ("someFunc" :: String)

ultimateTest = do
  let properShape v = AF.moddims v [size * size, 8]
      flattenBack v = AF.moddims v [size * size * 8]
      params = askarsParameters {pV0 = 0}
      size = 10
      x = mkPeriodicAxis (2 * pi / pk0 params) size
      y = mkPeriodicAxis (0 / 0) 1
      z = mkBoundedAxis (0, 1) size
      grid = SpaceGrid x y z
  -- (grid, qs, dqs) <- importFields "../test_data.h5"
  noise <- AF.randu [1] -- AF.randu @Double [size * size, 8]
  let q = reissnerNordstromFields params grid
      q' = Fields $ unFields q + 0.2 * (noise - 0.5)
  print $ evalEquations params grid (unFields q')
  newtonRaphson
    (RootOptions (\_ r -> r < 1.0e-5) 10 Nothing)
    (pure . flattenBack . evalEquations params grid . properShape)
    (flattenBack $ unFields q')

{- ORMOLU_DISABLE -}
foreign import ccall unsafe "kernels_metric.h ads_cft_halide_compute_metric"
  c_compute_metric
    :: Double -> Double -> -- length, μ
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- x, y, z
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- Q, ∂Q, ∂∂Q
    Ptr HalideBuffer -> Ptr HalideBuffer -> -- g₁₂, g¹²
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- ∂₁g₂₃, ∂₁g²³, ∂₁₂g₃₄
    IO CInt
{- ORMOLU_ENABLE -}

computeMetric ::
  Double ->
  Double ->
  FlatSpaceGrid Double ->
  Fields Double ->
  FieldsDerivatives Double ->
  (Metric Double, MetricDerivatives Double)
computeMetric l μ grid fields fieldsDerivatives =
  unsafePerformIO $ do
    let n = gridNumPoints grid
        metric = Metric (AF.constant [n, 4, 4] 0) (AF.constant [n, 4, 4] 0)
        metricDerivatives =
          MetricDerivatives
            (AF.constant [n, 4, 4, 4] 0)
            (AF.constant [n, 4, 4, 4] 0)
            (AF.constant [n, 4 * 4 * 4 * 4] 0)
    code <-
      withFlatSpaceGrid grid $ \c_x c_y c_z ->
        withFields fields $ \c_q ->
          withFieldsDerivatives fieldsDerivatives $ \c_dq c_ddq ->
            withMetric metric $ \c_g_dd c_g_UU ->
              withMetricDerivatives metricDerivatives $ \c_dg_ddd c_dg_dUU c_ddg_dddd ->
                c_compute_metric
                  l
                  μ
                  c_x
                  c_y
                  c_z
                  c_q
                  c_dq
                  c_ddq
                  c_g_dd
                  c_g_UU
                  c_dg_ddd
                  c_dg_dUU
                  c_ddg_dddd
    unless (code == 0) $
      error $ "Halide failed with error code: " <> show code
    pure (metric, metricDerivatives)

{- ORMOLU_DISABLE -}
foreign import ccall unsafe "kernels_metric.h ads_cft_halide_compute_christoffel"
  c_compute_christoffel
    :: Ptr HalideBuffer -> -- g¹²
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- ∂₁g₂₃, ∂₁g²³, ∂₁₂g₃₄
    Ptr HalideBuffer -> Ptr HalideBuffer -> -- Γ¹₂₃, ∂₁Γ²₃₄
    IO CInt
{- ORMOLU_ENABLE -}

computeChristoffel ::
  Metric Double ->
  MetricDerivatives Double ->
  Christoffel Double
computeChristoffel metric metricDerivatives =
  unsafePerformIO $ do
    let n = metricNumPoints metric
        christoffel = Christoffel (AF.constant [n, 4, 4, 4] 0) (AF.constant [n, 4 * 4 * 4 * 4] 0)
    code <-
      withMetric metric $ \c_g_dd c_g_UU ->
        withMetricDerivatives metricDerivatives $ \c_dg_ddd c_dg_dUU c_ddg_dddd ->
          withChristoffel christoffel $ \c_Γ_Udd c_dΓ_dUdd ->
            c_compute_christoffel
              c_g_UU
              c_dg_ddd
              c_dg_dUU
              c_ddg_dddd
              c_Γ_Udd
              c_dΓ_dUdd
    unless (code == 0) $
      error $ "Halide failed with error code: " <> show code
    pure christoffel

{- ORMOLU_DISABLE -}
foreign import ccall unsafe "kernels_metric.h ads_cft_halide_compute_maxwell"
  c_compute_maxwell
    :: Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- x, y, z
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- Q, ∂Q, ∂∂Q
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- g¹², ∂₁g₂₃, ∂₁g²³
    Ptr HalideBuffer -> -- Γ¹₂₃
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- F₁₂, F¹₂, F¹²
    Ptr HalideBuffer -> -- ∇F₁
    IO CInt
{- ORMOLU_ENABLE -}

computeMaxwell ::
  FlatSpaceGrid Double ->
  Fields Double ->
  FieldsDerivatives Double ->
  Metric Double ->
  MetricDerivatives Double ->
  Christoffel Double ->
  Maxwell Double
computeMaxwell grid fields fieldsDerivatives metric metricDerivatives christoffel =
  unsafePerformIO $ do
    let n = gridNumPoints grid
        maxwell =
          Maxwell
            (AF.constant [n, 4, 4] 0)
            (AF.constant [n, 4, 4] 0)
            (AF.constant [n, 4, 4] 0)
            (AF.constant [n, 4] 0)
    code <-
      withFlatSpaceGrid grid $ \c_x c_y c_z ->
        withFields fields $ \c_q ->
          withFieldsDerivatives fieldsDerivatives $ \c_dq c_ddq ->
            withMetric metric $ \c_g_dd c_g_UU ->
              withMetricDerivatives metricDerivatives $ \c_dg_ddd c_dg_dUU c_ddg_dddd ->
                withChristoffel christoffel $ \c_Γ_Udd c_dΓ_dUdd ->
                  withMaxwell maxwell $ \c_f_dd c_f_Ud c_f_UU c_divf_d ->
                    c_compute_maxwell
                      c_x
                      c_y
                      c_z
                      c_q
                      c_dq
                      c_ddq
                      c_g_UU
                      c_dg_ddd
                      c_dg_dUU
                      c_Γ_Udd
                      c_f_dd
                      c_f_Ud
                      c_f_UU
                      c_divf_d
    unless (code == 0) $
      error $ "Halide failed with error code: " <> show code
    pure maxwell

{- ORMOLU_DISABLE -}
foreign import ccall unsafe "kernels_metric.h ads_cft_halide_compute_deturck"
  c_compute_deturck
    :: Double -> Double -> -- L, μ
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- x, y, z
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- Q, ∂Q, ∂∂Q
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- g₁₂, g¹², ∂₁g₁₂, ∂₁g¹²
    Ptr HalideBuffer -> Ptr HalideBuffer -> -- Γ¹₂₃, ∂₁Γ²₃₄
    Ptr HalideBuffer -> -- ∇ξ₁₂
    IO CInt
{- ORMOLU_ENABLE -}

computeDeTurck ::
  Double ->
  Double ->
  FlatSpaceGrid Double ->
  Fields Double ->
  FieldsDerivatives Double ->
  Metric Double ->
  MetricDerivatives Double ->
  Christoffel Double ->
  DeTurck Double
computeDeTurck l μ grid fields fieldsDerivatives metric metricDerivatives christoffel =
  unsafePerformIO $ do
    let deturck =
          DeTurck
            (AF.constant [fieldsNumPoints fields, 4, 4] 0)
    code <-
      withFlatSpaceGrid grid $ \c_x c_y c_z ->
        withFields fields $ \c_q ->
          withFieldsDerivatives fieldsDerivatives $ \c_dq c_ddq ->
            withMetric metric $ \c_g_dd c_g_UU ->
              withMetricDerivatives metricDerivatives $ \c_dg_ddd c_dg_dUU c_ddg_dddd ->
                withChristoffel christoffel $ \c_Γ_Udd c_dΓ_dUdd ->
                  withDeTurck deturck $ \c_divξ_dd ->
                    c_compute_deturck
                      l
                      μ
                      c_x
                      c_y
                      c_z
                      c_q
                      c_dq
                      c_ddq
                      c_g_dd
                      c_g_UU
                      c_dg_ddd
                      c_dg_dUU
                      c_Γ_Udd
                      c_dΓ_dUdd
                      c_divξ_dd
    unless (code == 0) $
      error $ "Halide failed with error code: " <> show code
    pure deturck

{- ORMOLU_DISABLE -}
foreign import ccall unsafe "kernels_metric.h ads_cft_halide_evaluate_equations"
  c_evaluate_equations
    :: Double -> -- length
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- x, y, z
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- Q, ∂Q, ∂∂Q
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- g₁₂, g¹², ∂₁g¹²
    Ptr HalideBuffer -> Ptr HalideBuffer -> -- Γ¹₂₃, ∂₁Γ²₃₄
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- F₁₂, F¹₂, F¹²
    Ptr HalideBuffer -> -- ∇F₁
    Ptr HalideBuffer -> -- ∇ξ₁₂
    Ptr HalideBuffer -> -- equations
    IO CInt
{- ORMOLU_ENABLE -}

evaluateEquations ::
  Double ->
  FlatSpaceGrid Double ->
  Fields Double ->
  FieldsDerivatives Double ->
  Metric Double ->
  MetricDerivatives Double ->
  Christoffel Double ->
  Maxwell Double ->
  DeTurck Double ->
  Equations Double
evaluateEquations l grid fields fieldsDerivatives metric metricDerivatives christoffel maxwell deturck =
  unsafePerformIO $ do
    let equations =
          Equations
            (AF.constant [fieldsNumPoints fields, fieldsNumParams fields] 0)
    code <-
      withFlatSpaceGrid grid $ \c_x c_y c_z ->
        withFields fields $ \c_q ->
          withFieldsDerivatives fieldsDerivatives $ \c_dq c_ddq ->
            withMetric metric $ \c_g_dd c_g_UU ->
              withMetricDerivatives metricDerivatives $ \c_dg_ddd c_dg_dUU c_ddg_dddd ->
                withChristoffel christoffel $ \c_Γ_Udd c_dΓ_dUdd ->
                  withMaxwell maxwell $ \c_f_dd c_f_Ud c_f_UU c_divf_d ->
                    withDeTurck deturck $ \c_divξ_dd ->
                      withEquations equations $ \c_eq ->
                        c_evaluate_equations
                          l
                          c_x
                          c_y
                          c_z
                          c_q
                          c_dq
                          c_ddq
                          c_g_dd
                          c_g_UU
                          c_dg_dUU
                          c_Γ_Udd
                          c_dΓ_dUdd
                          c_f_dd
                          c_f_Ud
                          c_f_UU
                          c_divf_d
                          c_divξ_dd
                          c_eq
    unless (code == 0) $
      error $ "Halide failed with error code: " <> show code
    pure equations

{- ORMOLU_DISABLE -}
foreign import ccall unsafe "kernels_metric.h ads_cft_halide_horizon_boundary_conditions"
  c_horizon_boundary_conditions
    :: Double -> Double -> -- L, μ
    Ptr HalideBuffer -> Ptr HalideBuffer -> -- x, y
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- Q, ∂Q, ∂∂Q
    Ptr HalideBuffer -> -- horizon
    IO CInt
{- ORMOLU_ENABLE -}

horizonBoundaryConditions ::
  Double ->
  Double ->
  FlatSpaceGrid Double ->
  Fields Double ->
  FieldsDerivatives Double ->
  Array Double
horizonBoundaryConditions l μ grid fields fieldsDerivatives =
  unsafePerformIO $ do
    let horizon = AF.constant [fieldsNumPoints fields, fieldsNumParams fields] 0
    code <-
      withFlatSpaceGrid grid $ \c_x c_y c_z ->
        withFields fields $ \c_q ->
          withFieldsDerivatives fieldsDerivatives $ \c_dq c_ddq ->
            withHalideBuffer horizon $ \c_horizon ->
              c_horizon_boundary_conditions
                l
                μ
                c_x
                c_y
                c_q
                c_dq
                c_ddq
                c_horizon
    unless (code == 0) $
      error $ "Halide failed with error code: " <> show code
    pure horizon

{- ORMOLU_DISABLE -}
foreign import ccall unsafe "kernels_metric.h ads_cft_halide_conformal_boundary_conditions"
  c_conformal_boundary_conditions
    :: Double -> Double -> Double -> Double -> Double -> -- L, μ, V0, k0, θ
    Ptr HalideBuffer -> Ptr HalideBuffer -> -- x, y
    Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> -- Q, ∂Q, ∂∂Q
    Ptr HalideBuffer -> -- conformal
    IO CInt
{- ORMOLU_ENABLE -}

conformalBoundaryConditions ::
  Double ->
  Double ->
  Double ->
  Double ->
  Double ->
  FlatSpaceGrid Double ->
  Fields Double ->
  FieldsDerivatives Double ->
  Array Double
conformalBoundaryConditions l μ v0 k0 θ grid fields fieldsDerivatives =
  unsafePerformIO $ do
    let conformal = AF.constant [fieldsNumPoints fields, fieldsNumParams fields] 0
    code <-
      withFlatSpaceGrid grid $ \c_x c_y c_z ->
        withFields fields $ \c_q ->
          withFieldsDerivatives fieldsDerivatives $ \c_dq c_ddq ->
            withHalideBuffer conformal $ \c_conformal ->
              c_conformal_boundary_conditions
                l
                μ
                v0
                k0
                θ
                c_x
                c_y
                c_q
                c_dq
                c_ddq
                c_conformal
    unless (code == 0) $
      error $ "Halide failed with error code: " <> show code
    pure conformal

data Axis a = Axis {axisGrid :: !(Array a), axisDiff :: !(Array a)}
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

mkPeriodicAxis :: (AFType a, Floating a) => a -> Int -> Axis a
mkPeriodicAxis period numPoints =
  Axis (gridPointsForPeriodic period numPoints) (differentiationMatrixPeriodic period numPoints)

mkBoundedAxis :: (AFType a, Floating a) => (a, a) -> Int -> Axis a
mkBoundedAxis range numPoints =
  Axis (gridPointsForBounded range numPoints) (differentiationMatrixBounded range numPoints)

axisNumPoints :: AFType a => Axis a -> Int
axisNumPoints (Axis grid _) = AF.getElements grid

data SpaceGrid a = SpaceGrid
  { gridX :: !(Axis a),
    gridY :: !(Axis a),
    gridZ :: !(Axis a)
  }
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

gridShape :: AFType a => SpaceGrid a -> (Int, Int, Int)
gridShape (SpaceGrid x y z) = (axisNumPoints x, axisNumPoints y, axisNumPoints z)

data FlatSpaceGrid a = FlatSpaceGrid
  { flatGridX :: !(Array a),
    flatGridY :: !(Array a),
    flatGridZ :: !(Array a)
  }
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

gridNumPoints :: AFType a => FlatSpaceGrid a -> Int
gridNumPoints (FlatSpaceGrid x _ _) = AF.getElements x

newtype Fields a = Fields {unFields :: Array a}
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

newtype Equations a = Equations {unEquations :: Array a}
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

newtype DeTurck a = DeTurck {unDeTurck :: Array a}
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

data FieldsDerivatives a = FieldsDerivatives
  { fieldsDerivative :: !(Array a),
    fieldsSecondDerivative :: !(Array a)
  }
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

fieldsNumPoints :: AFType a => Fields a -> Int
fieldsNumPoints fields = n
  where
    (n, _, _, _) = AF.getDims (unFields fields)

fieldsNumParams :: AFType a => Fields a -> Int
fieldsNumParams fields = n
  where
    (_, n, _, _) = AF.getDims (unFields fields)

data Metric a = Metric
  { metricDD :: !(Array a),
    metricUU :: !(Array a)
  }
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

metricNumPoints :: AFType a => Metric a -> Int
metricNumPoints metric = n
  where
    (n, _, _, _) = AF.getDims (metricDD metric)

data MetricDerivatives a = MetricDerivatives
  { metricDerivativeDDD :: !(Array a),
    metricDerivativeDUU :: !(Array a),
    metricSecondDerivativeDDDD :: !(Array a)
  }
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

data Maxwell a = Maxwell
  { maxwellDD :: !(Array a),
    maxwellUD :: !(Array a),
    maxwellUU :: !(Array a),
    maxwellDivD :: !(Array a)
  }
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

data Christoffel a = Christoffel
  { christoffelUDD :: !(Array a),
    christoffelDerivativeDUDD :: !(Array a)
  }
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

data SystemParameters a = SystemParameters
  { pL :: !a,
    pμ :: !a,
    pV0 :: !a,
    pk0 :: !a,
    pθ :: !a
  }
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

askarsParameters :: SystemParameters Double
askarsParameters =
  SystemParameters
    { pL = 1,
      pμ = 2.3,
      pV0 = 5,
      pk0 = 2.1,
      pθ = 0.785 -- pi / 4
    }

flattenGrid :: AFType a => SpaceGrid a -> FlatSpaceGrid a
flattenGrid (SpaceGrid (Axis gridX _) (Axis gridY _) (Axis gridZ _)) = FlatSpaceGrid x y z
  where
    (numX, _, _, _) = AF.getDims gridX
    (numY, _, _, _) = AF.getDims gridY
    (numZ, _, _, _) = AF.getDims gridZ
    !x = AF.flat $ AF.tile gridX [1, numY, numZ]
    !y = AF.flat $ AF.tile (AF.moddims gridY [1, numY]) [numX, 1, numZ]
    !z = AF.flat $ AF.tile (AF.moddims gridZ [1, 1, numZ]) [numX, numY, 1]

instance NFData Seq where
  rnf (Seq b e s) = b `deepseq` e `deepseq` (rnf s)

data Selector = Selector
  { selectorIndex :: !Seq,
    selectorShape :: !(Int, Int, Int)
  }
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

runSelector :: AFType a => Selector -> Array a -> Array a
runSelector (Selector !i (n₀, n₁, n₂)) arr
  | n₀ * n₁ * n₂ == d₀ = reshapeBack . index . reshape $ arr
  | otherwise =
    error $
      "runSelector failed: array has shape " <> show shape <> ", but the grid has shape " <> show (n₀, n₁, n₂)
  where
    shape@(d₀, d₁, d₂, d₃) = AF.getDims arr
    reshape = flip AF.moddims [n₀, n₁, n₂, d₁ * d₂ * d₃]
    index = flip AF.index [Seq 0 (-1) 1, Seq 0 (-1) 1, i, Seq 0 (-1) 1]
    reshapeBack t = AF.moddims t [n₀ * n₁ * k, d₁, d₂, d₃]
      where
        (_, _, k, _) = AF.getDims t

bulkSelector :: (AFType a, Num a, Eq a) => SpaceGrid a -> Selector
bulkSelector (SpaceGrid (Axis gridX _) (Axis gridY _) (Axis gridZ _))
  -- Check that we have a valid z grid
  | AF.index gridZ [Seq 0 0 1] == AF.scalar 1
      && AF.index gridZ [Seq (-1) (-1) 1] == AF.scalar 0 =
    Selector (Seq 1 (-2) 1) (AF.getElements gridX, AF.getElements gridY, AF.getElements gridZ)
  | otherwise = error "bulkSelector: invalid z grid"

conformalSelector :: (AFType a, Num a, Eq a) => SpaceGrid a -> Selector
conformalSelector (SpaceGrid (Axis gridX _) (Axis gridY _) (Axis gridZ _))
  -- Check that we have a valid z grid
  | AF.index gridZ [Seq 0 0 1] == AF.scalar 1
      && AF.index gridZ [Seq (-1) (-1) 1] == AF.scalar 0 =
    Selector (Seq (-1) (-1) 1) (AF.getElements gridX, AF.getElements gridY, AF.getElements gridZ)
  | otherwise = error "conformalSelector: invalid z grid"

horizonSelector :: (AFType a, Num a, Eq a) => SpaceGrid a -> Selector
horizonSelector (SpaceGrid (Axis gridX _) (Axis gridY _) (Axis gridZ _))
  -- Check that we have a valid z grid
  | AF.index gridZ [Seq 0 0 1] == AF.scalar 1
      && AF.index gridZ [Seq (-1) (-1) 1] == AF.scalar 0 =
    Selector (Seq 0 0 1) (AF.getElements gridX, AF.getElements gridY, AF.getElements gridZ)
  | otherwise = error "horizonSelector: invalid z grid"

-- conformalBoundarySelector :: (AFType a, Num a) => SpaceGrid a -> Selector a
-- conformalBoundarySelector (SpaceGrid _ _ z) = Selector isBoundary boundaryIndices
--   where
--     isBoundary = AF.eq z (scalar 0)
--     boundaryIndices = AF.where' $ AF.cast isBoundary
--
-- horizonBoundarySelector :: (AFType a, Num a) => SpaceGrid a -> Selector a
-- horizonBoundarySelector (SpaceGrid _ _ z) = Selector isBoundary boundaryIndices
--   where
--     isBoundary = AF.eq z (scalar 1)
--     boundaryIndices = AF.where' $ AF.cast isBoundary

class IsSelectable a where
  select :: Selector -> a -> a

instance AFType a => IsSelectable (FlatSpaceGrid a) where
  select s (FlatSpaceGrid x y z) = FlatSpaceGrid (runSelector s x) (runSelector s y) (runSelector s z)

instance AFType a => IsSelectable (Fields a) where
  select s (Fields qs) = Fields (runSelector s qs)

instance AFType a => IsSelectable (FieldsDerivatives a) where
  select s (FieldsDerivatives dqs ddqs) = FieldsDerivatives (runSelector s dqs) (runSelector s ddqs)

withMaxwell ::
  (AFType a, IsHalideType a) =>
  Maxwell a ->
  (Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> IO b) ->
  IO b
withMaxwell (Maxwell f_dd f_Ud f_UU divf_d) action =
  withHalideBuffer f_dd $ \c_f_dd ->
    withHalideBuffer f_Ud $ \c_f_Ud ->
      withHalideBuffer f_UU $ \c_f_UU ->
        withHalideBuffer divf_d $ \c_divf_d ->
          action c_f_dd c_f_Ud c_f_UU c_divf_d

withChristoffel ::
  (AFType a, IsHalideType a) =>
  Christoffel a ->
  (Ptr HalideBuffer -> Ptr HalideBuffer -> IO b) ->
  IO b
withChristoffel (Christoffel γ_Udd dγ_dUdd) action =
  withHalideBuffer γ_Udd $ \c_Γ_Udd ->
    withHalideBufferOfShape dγ_dUdd [n, 4, 4, 4, 4] $ \c_dΓ_dUdd ->
      action c_Γ_Udd c_dΓ_dUdd
  where
    (n, _, _, _) = AF.getDims γ_Udd

rowMajorStrides :: Integral a => [a] -> [a]
rowMajorStrides = drop 1 . scanr (*) 1

colMajorStrides :: Integral a => [a] -> [a]
colMajorStrides = go 1
  where
    go _ [] = []
    go !s (!d : ds) = s : go (s * d) ds

withHalideBufferOfShape ::
  (HasCallStack, IsHalideType a, AFType a) =>
  Array a ->
  [Int] ->
  (Ptr HalideBuffer -> IO b) ->
  IO b
withHalideBufferOfShape arr shape action
  | d₀ * d₁ * d₂ * d₃ == product shape && backend == AF.CPU =
    AF.withDevicePtr arr $ \devicePtr ->
      bufferFromPtrShapeStrides devicePtr shape (colMajorStrides shape) action
  | backend /= AF.CPU = error $ "unsupported backend: " <> show backend
  | otherwise = error $ "cannot reshape Array of size " <> show currentShape <> " into " <> show shape
  where
    backend = AF.getBackend arr
    currentShape@(d₀, d₁, d₂, d₃) = AF.getDims arr

instance (AFType a, IsHalideType a) => IsHalideBuffer (Array a) where
  withHalideBuffer arr action
    | backend == AF.CPU =
      AF.withDevicePtr arr $ \devicePtr ->
        bufferFromPtrShapeStrides devicePtr shape (colMajorStrides shape) action
    | otherwise = error $ "unsupported backend: " <> show backend
    where
      backend = AF.getBackend arr
      currentShape@(d₀, d₁, d₂, d₃) = AF.getDims arr
      shape = take (AF.getNumDims arr) [d₀, d₁, d₂, d₃]

withFlatSpaceGrid ::
  (AFType a, IsHalideType a) =>
  FlatSpaceGrid a ->
  (Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> IO b) ->
  IO b
withFlatSpaceGrid (FlatSpaceGrid x y z) action =
  withHalideBuffer x $ \c_x ->
    withHalideBuffer y $ \c_y ->
      withHalideBuffer z $ \c_z ->
        action c_x c_y c_z

-- withSpaceGrid ::
--   (AFType a, IsHalideType a) =>
--   SpaceGrid a ->
--   (Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> IO b) ->
--   IO b
-- withSpaceGrid (SpaceGrid x y z) action =
--   withHalideBuffer x $ \c_x ->
--     withHalideBuffer y $ \c_y ->
--       withHalideBuffer z $ \c_z ->
--         action c_x c_y c_z

withFields ::
  (AFType a, IsHalideType a) =>
  Fields a ->
  (Ptr HalideBuffer -> IO b) ->
  IO b
withFields (Fields q) action = withHalideBuffer q action

withFieldsDerivatives ::
  (AFType a, IsHalideType a) =>
  FieldsDerivatives a ->
  (Ptr HalideBuffer -> Ptr HalideBuffer -> IO b) ->
  IO b
withFieldsDerivatives (FieldsDerivatives dq ddq) action =
  withHalideBuffer dq $ \c_dq ->
    withHalideBuffer ddq $ \c_ddq ->
      action c_dq c_ddq

withEquations ::
  (AFType a, IsHalideType a) =>
  Equations a ->
  (Ptr HalideBuffer -> IO b) ->
  IO b
withEquations (Equations eq) action = withHalideBuffer eq action

withDeTurck ::
  (AFType a, IsHalideType a) =>
  DeTurck a ->
  (Ptr HalideBuffer -> IO b) ->
  IO b
withDeTurck (DeTurck divξ_dd) action = withHalideBuffer divξ_dd action

withMetric ::
  (AFType a, IsHalideType a) =>
  Metric a ->
  (Ptr HalideBuffer -> Ptr HalideBuffer -> IO b) ->
  IO b
withMetric (Metric g_dd g_UU) action =
  withHalideBuffer g_dd $ \c_g_dd ->
    withHalideBuffer g_UU $ \c_g_UU ->
      action c_g_dd c_g_UU

withMetricDerivatives ::
  (AFType a, IsHalideType a) =>
  MetricDerivatives a ->
  (Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> IO b) ->
  IO b
withMetricDerivatives (MetricDerivatives dg_ddd dg_dUU ddg_dddd) action =
  withHalideBuffer dg_ddd $ \c_dg_ddd ->
    withHalideBuffer dg_dUU $ \c_dg_dUU ->
      withHalideBufferOfShape ddg_dddd [n, 4, 4, 4, 4] $ \c_ddg_dddd ->
        action c_dg_ddd c_dg_dUU c_ddg_dddd
  where
    (n, _, _, _) = AF.getDims dg_ddd

-- foo :: (AFType a, IsHalideType a, Num a) => [Int] -> Array a
-- foo shape = unsafePerformIO $ do
--   let arr = AF.mkArray shape (repeat 0)
--   code <- withArrayOfShape arr shape c_foo
--   unless (code == 0) $
--     error $ "Halide failed with error code: " <> show code
--   pure arr

withArrayViewOfShape ::
  (HasCallStack, MonadResource m, AFType a) =>
  Array a ->
  [Int] ->
  (ArrayView' a -> m b) ->
  m b
withArrayViewOfShape arr shape action
  | AF.getElements arr == size && backend == AF.CPU = do
    -- This is pretty ugly, and hdf5-hs should be updated to support
    -- MonadUnliftIO instead of just MonadResource
    fp <- liftIO $ AF.withDevicePtr arr newForeignPtr_
    key <- register (touch arr)
    r <- action $ ArrayView' fp (reverse shape) (rowMajorStrides (reverse shape))
    release key
    pure r
  | backend /= AF.CPU = error $ "unsupported backend: " <> show backend
  | otherwise =
    error $
      "cannot reshape Array of size " <> show (AF.getDims arr) <> " into " <> show shape
  where
    backend = AF.getBackend arr
    size = product shape

fromArrayViewOfShape ::
  (HasCallStack, MonadResource m, AFType a) =>
  ArrayView' a ->
  [Int] ->
  m (Array a)
fromArrayViewOfShape (ArrayView' fptr shape strides) shape'
  | strides == colMajorStrides shape =
    liftIO $
      withForeignPtr fptr $ \hostPtr ->
        AF.unsafeFromHostPtrShape hostPtr shape'
  | strides == rowMajorStrides shape =
    fmap rowToColMajor $
      liftIO $
        withForeignPtr fptr $ \hostPtr ->
          AF.unsafeFromHostPtrShape hostPtr (reverse shape')
  | otherwise = error "cannot create an Array from a strided ArrayView"
  where
    rowToColMajor arr =
      let ndims = AF.getNumDims arr
       in AF.reorder arr (reverse [0 .. ndims - 1] <> [ndims .. 3])

fromArrayView ::
  (HasCallStack, MonadResource m, AFType a) =>
  ArrayView' a ->
  m (Array a)
fromArrayView arr@(ArrayView' _ shape _) =
  fromArrayViewOfShape arr shape

withArrayView ::
  (HasCallStack, MonadResource m, AFType a) =>
  Array a ->
  (ArrayView' a -> m b) ->
  m b
withArrayView arr = withArrayViewOfShape arr (take (AF.getNumDims arr) [d₀, d₁, d₂, d₃])
  where
    (d₀, d₁, d₂, d₃) = AF.getDims arr

type instance H5.ElementOf (Array a) = a

instance (AFType a, H5.KnownDatatype a) => H5.KnownDataset' (Array a) where
  withArrayView' = withArrayView
  fromArrayView' = fromArrayView

instance NFData (Array a) where
  rnf (AF.Array (GHC.ForeignPtr !addr !contents)) = ()

importFields :: Text -> IO (SpaceGrid Double, Fields Double, FieldsDerivatives Double)
importFields filename = H5.withFile filename H5.ReadOnly $ \file -> do
  -- xGrid <- fromArrayView =<< H5.readDataset =<< H5.open file "x"
  -- yGrid <- fromArrayView =<< H5.readDataset =<< H5.open file "y"
  -- zGrid <- fromArrayView =<< H5.readDataset =<< H5.open file "z"
  let params = askarsParameters
      x = mkPeriodicAxis (2 * pi / pk0 params) 6
      y = mkPeriodicAxis (0 / 0) 1
      z = mkBoundedAxis (0, 1) 6
      grid = SpaceGrid x y z
  qs <- fromArrayView =<< H5.readDataset =<< H5.open file "Qs"
  dqs <- fromArrayView =<< H5.readDataset =<< H5.open file "DQs"
  ddqs <- fromArrayView =<< H5.readDataset =<< H5.open file "DDQs"
  pure (grid, Fields qs, FieldsDerivatives dqs ddqs)

evalEquationsFromAskarsInputs :: Text -> IO (Array Double)
evalEquationsFromAskarsInputs filename = H5.withFile filename H5.ReadOnly $ \file -> do
  let params = askarsParameters
      x = mkPeriodicAxis (2 * pi / pk0 params) 6
      y = mkPeriodicAxis undefined 1
      z = mkBoundedAxis (0, 1) 6
      grid = SpaceGrid x y z
  qs <- fmap Fields $ fromArrayView =<< H5.readDataset =<< H5.open file "Qs"
  let dqs = computeDerivatives grid qs
      horizon = horizonEquations params grid qs dqs
      bulk = bulkEquations params grid qs dqs
      conformal = conformalEquations params grid qs dqs
      !eq = stackAlongZ grid (horizon :| [bulk, conformal])
  pure eq

importAskarsFields :: Text -> IO (SpaceGrid Double, Fields Double, FieldsDerivatives Double, Array Double)
importAskarsFields filename = H5.withFile filename H5.ReadOnly $ \file -> do
  let params = askarsParameters
      x = mkPeriodicAxis (2 * pi / pk0 params) 6
      y = mkPeriodicAxis undefined 1
      z = mkBoundedAxis (0, 1) 6
      grid = SpaceGrid x y z
  qs <- fromArrayView =<< H5.readDataset =<< H5.open file "Qs"
  dqs <- fromArrayView =<< H5.readDataset =<< H5.open file "DQs"
  ddqs <- fromArrayView =<< H5.readDataset =<< H5.open file "DDQs"
  eq <- fromArrayView =<< H5.readDataset =<< H5.open file "equations"
  pure (grid, Fields qs, FieldsDerivatives dqs ddqs, eq)

computeDerivatives :: (AFType a, Num a) => SpaceGrid a -> Fields a -> FieldsDerivatives a
computeDerivatives grid@(SpaceGrid (Axis gridX dX) (Axis _ dY) (Axis _ dZ)) (Fields qs) =
  FieldsDerivatives dqs ddqs
  where
    (d₀, d₁, d₂) = gridShape grid
    n = d₀ * d₁ * d₂
    diff f = foldl' (AF.join 1) d₀f [d₁f, d₂f, d₃f]
      where
        (_, _, _, m) = AF.getDims f
        d₀f = AF.constant [n, 1, m, 1] 0
        d₁f = AF.moddims (differentiateX dX f) [n, 1, m, 1]
        d₂f = AF.moddims (differentiateY dY f) [n, 1, m, 1]
        d₃f = AF.moddims (differentiateZ dZ f) [n, 1, m, 1]
    !dqs = diff qs'
      where
        (_, m, _, _) = AF.getDims qs
        qs' = AF.moddims qs [d₀, d₁, d₂, m]
    !ddqs = AF.moddims (diff dqs') [n, 4, 4, m]
      where
        (_, _, m, _) = AF.getDims dqs
        dqs' = AF.moddims dqs [d₀, d₁, d₂, 4 * m]

bulkEquations ::
  SystemParameters Double ->
  SpaceGrid Double ->
  Fields Double ->
  FieldsDerivatives Double ->
  Array Double
bulkEquations (SystemParameters l μ _ _ _) fullGrid fullFields fullDFields = unEquations eq
  where
    bulk = bulkSelector fullGrid
    grid = select bulk (flattenGrid fullGrid)
    fields = select bulk fullFields
    dfields = select bulk fullDFields
    (metric, dmetric) = computeMetric l μ grid fields dfields
    christoffel = computeChristoffel metric dmetric
    maxwell = computeMaxwell grid fields dfields metric dmetric christoffel
    deturck = computeDeTurck l μ grid fields dfields metric dmetric christoffel
    !eq = evaluateEquations l grid fields dfields metric dmetric christoffel maxwell deturck

horizonEquations ::
  SystemParameters Double ->
  SpaceGrid Double ->
  Fields Double ->
  FieldsDerivatives Double ->
  Array Double
horizonEquations (SystemParameters l μ _ _ _) fullGrid fullFields fullDFields = eq
  where
    horizon = horizonSelector fullGrid
    grid = select horizon (flattenGrid fullGrid)
    fields = select horizon fullFields
    dfields = select horizon fullDFields
    !eq = horizonBoundaryConditions l μ grid fields dfields

conformalEquations ::
  SystemParameters Double ->
  SpaceGrid Double ->
  Fields Double ->
  FieldsDerivatives Double ->
  Array Double
conformalEquations (SystemParameters l μ v0 k0 θ) fullGrid fullFields fullDFields = eq
  where
    conformal = conformalSelector fullGrid
    grid = select conformal (flattenGrid fullGrid)
    fields = select conformal fullFields
    dfields = select conformal fullDFields
    !eq = conformalBoundaryConditions l μ v0 k0 θ grid fields dfields

stackAlongZ :: AFType a => SpaceGrid a -> NonEmpty (Array a) -> Array a
stackAlongZ grid (arr :| arrs) = foldl' stackTwo arr arrs
  where
    stackTwo !a !b
      | (a₁, a₂, a₃) == (b₁, b₂, b₃)
          && a₀ `mod` (n₀ * n₁) == 0
          && b₀ `mod` (n₀ * n₁) == 0 =
        AF.moddims (AF.join 2 a' b') [a₀ + b₀, a₁, a₂, a₃]
      | otherwise =
        error $
          "could not stack arrays of shape "
            <> show (AF.getDims a)
            <> " and "
            <> show (AF.getDims b)
            <> " along z axis; grid size is "
            <> show (n₀, n₁, n₂)
      where
        shapeA@(a₀, a₁, a₂, a₃) = AF.getDims a
        shapeB@(b₀, b₁, b₂, b₃) = AF.getDims b
        (n₀, n₁, n₂) = gridShape grid
        a' = AF.moddims a [n₀, n₁, a₀ `div` (n₀ * n₁), a₁ * a₂ * a₃]
        b' = AF.moddims b [n₀, n₁, b₀ `div` (n₀ * n₁), b₁ * b₂ * b₃]

evalEquations :: SystemParameters Double -> SpaceGrid Double -> Array Double -> Array Double
evalEquations params grid q = stackAlongZ grid (fromList [horizon, bulk, conformal])
  where
    fields = Fields q
    dfields = computeDerivatives grid fields
    !bulk = bulkEquations params grid fields dfields
    !horizon = horizonEquations params grid fields dfields
    !conformal = conformalEquations params grid fields dfields

reissnerNordstromFields :: (AFType a, Num a) => SystemParameters a -> SpaceGrid a -> Fields a
reissnerNordstromFields params grid = Fields q
  where
    μ = pμ params
    q = AF.tile (AF.moddims (AF.vector 8 [1, 1, 1, 1, 0, μ, 0, 0]) [1, 8]) [d₀ * d₁ * d₂]
    (d₀, d₁, d₂) = gridShape grid

evalEquationsFromInput = do
  (grid, fields, dfields) <- importFields "../test_data.h5"
  let params = askarsParameters
      !bulk = bulkEquations params grid fields dfields
      !horizon = horizonEquations params grid fields dfields
      !conformal = conformalEquations params grid fields dfields
  print conformal
  pure $ stackAlongZ grid (fromList [horizon, bulk, conformal])

importExpectedOutputs :: IO (Array Double)
importExpectedOutputs = H5.withFile "../test_data.h5" H5.ReadOnly $ \file -> do
  (eq :: Array Double) <- fromArrayView =<< H5.readDataset =<< H5.open file "equations"
  pure eq

gridPointsForPeriodic :: (AFType a, Floating a) => a -> Int -> Array a
gridPointsForPeriodic period n
  | n == 0 = AF.constant [] 0
  | n == 1 = AF.scalar (0 / 0)
  | even n = AF.scalar (period / fromIntegral n) * (AF.iota [n] [] + AF.scalar 1)
  | otherwise = error $ "invalid n: " <> show n

gridPointsForBounded :: (AFType a, Floating a) => (a, a) -> Int -> Array a
gridPointsForBounded (a, b) n
  | n >= 2 = AF.scalar ((b + a) / 2) + AF.scalar ((b - a) / 2) * AF.cos (scale * js)
  | otherwise = error $ "invalid n: " <> show n
  where
    scale = AF.scalar (pi / fromIntegral (n - 1))
    js = AF.iota [n] []

differentiationMatrixPeriodic :: forall a. (AFType a, Floating a) => a -> Int -> Array a
differentiationMatrixPeriodic period n
  | n == 0 = AF.constant [] 0
  | n == 1 = scalar 0
  | even n =
    AF.select
      isDiag
      (scalar @a 0)
      ( scalar (0.5 * (2 * pi) / period)
          * AF.cast (AF.pow (AF.scalar (-1)) δi)
          / AF.tan (scale * AF.cast δi)
      )
  | otherwise = error "currently only even n is supported"
  where
    scale = scalar $ pi / fromIntegral n
    rowIndices = AF.iota @Int [n] [1, n]
    colIndices = AF.transpose rowIndices False
    δi = rowIndices - colIndices
    isDiag = AF.eq rowIndices colIndices

differentiationMatrixBounded :: (AFType a, Floating a) => (a, a) -> Int -> Array a
differentiationMatrixBounded (l, u) n = AF.select isDiag diag offDiag
  where
    isDiag = AF.identity [n, n]
    xᵢ = AF.tile (gridPointsForBounded (l, u) n) [1, n]
    xⱼ = AF.transpose xᵢ False
    δx = xᵢ - xⱼ
    diag = flip AF.sum 1 $ AF.select isDiag (scalar 0) (scalar 1 / δx)
    a = flip AF.product 1 $ AF.select isDiag (scalar 1) δx
    aᵢ = AF.tile a [1, n]
    aⱼ = AF.transpose aᵢ False
    offDiag = aᵢ / (aⱼ * δx)

visualization :: IO ()
visualization = do
  let cars = dataFromUrl "https://vega.github.io/vega-datasets/data/cars.json" []
      enc =
        encoding
          . position X [PName "Horsepower", PmType Quantitative]
          . position Y [PName "Miles_per_Gallon", PmType Quantitative, PTitle "Miles per Gallon"]
      bkg = background "rgba(0, 0, 0, 0.05)"
  toHtmlFile "myplot.html" $
    toVegaLite [bkg, width 1024, height 920, cars, mark Circle [MTooltip TTEncoding], enc []]

{- ORMOLU_DISABLE -}
foreign import ccall unsafe "kernels_finite_differences.h ads_cft_halide_finite_differences_matrix"
  c_ads_cft_halide_finite_differences_matrix
    :: Ptr HalideBuffer -> -- x grid
    Ptr HalideBuffer -> -- Dₓ
    IO CInt
{- ORMOLU_ENABLE -}

finiteDifferencesMatrix :: Array Double -> Array Double
finiteDifferencesMatrix x =
  unsafePerformIO $ do
    let n = AF.getElements x
        matrix = AF.constant [n, n] 0
    code <-
      withHalideBuffer x $ \c_x ->
        withHalideBuffer matrix $ \c_matrix ->
          c_ads_cft_halide_finite_differences_matrix c_x c_matrix
    unless (code == 0) $
      error $ "Halide failed with error code: " <> show code
    pure matrix

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

differentiateY :: (HasCallStack, AFType a) => Array a -> Array a -> Array a
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
