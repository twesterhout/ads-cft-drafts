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
import Control.Monad.Primitive (touch)
import Control.Monad.Trans.Resource (MonadResource, register, release)
import Data.HDF5 (ArrayView' (..))
import qualified Data.HDF5 as H5
import Foreign.C.Types (CBool (..), CInt (..))
import Foreign.ForeignPtr (newForeignPtr_, withForeignPtr)
import Foreign.Ptr (Ptr)
import qualified GHC.ForeignPtr as GHC
import Language.Halide
import System.IO.Unsafe (unsafePerformIO)

someFunc :: IO ()
someFunc = do
  putStrLn ("someFunc" :: String)

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
  SpaceGrid Double ->
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
      withSpaceGrid grid $ \c_x c_y c_z ->
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
  SpaceGrid Double ->
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
      withSpaceGrid grid $ \c_x c_y c_z ->
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
  SpaceGrid Double ->
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
      withSpaceGrid grid $ \c_x c_y c_z ->
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
  SpaceGrid Double ->
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
      withSpaceGrid grid $ \c_x c_y c_z ->
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
  SpaceGrid Double ->
  Fields Double ->
  FieldsDerivatives Double ->
  Array Double
horizonBoundaryConditions l μ grid fields fieldsDerivatives =
  unsafePerformIO $ do
    let z = gridZ grid
        isBoundary = AF.eq z (scalar 1)
        boundaryIndices :: Array Int
        boundaryIndices = AF.where' $ AF.cast isBoundary
        grid' =
          SpaceGrid
            (AF.lookup (gridX grid) boundaryIndices 0)
            (AF.lookup (gridY grid) boundaryIndices 0)
            (AF.lookup (gridZ grid) boundaryIndices 0)
        fields' = Fields $ AF.lookup (unFields fields) boundaryIndices 0
        fieldsDerivatives' =
          FieldsDerivatives
            (AF.lookup (fieldsDerivative fieldsDerivatives) boundaryIndices 0)
            (AF.lookup (fieldsSecondDerivative fieldsDerivatives) boundaryIndices 0)
        horizon = AF.constant [fieldsNumPoints fields', fieldsNumParams fields'] 0
    code <-
      withSpaceGrid grid' $ \c_x c_y c_z ->
        withFields fields' $ \c_q ->
          withFieldsDerivatives fieldsDerivatives' $ \c_dq c_ddq ->
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

--   Input<double> _length{"length"};
--   Input<double> _chemical_potential{"chemical_potential"};
--   Input<Buffer<double>> _x{"x", 1};
--   Input<Buffer<double>> _y{"y", 1};
--   Input<Buffer<double>> _Q{"Q", 2};
--   Input<Buffer<double>> _DQ{"DQ", 3};
--   Input<Buffer<double>> _DDQ{"DDQ", 4};
--
--   Output<Buffer<double>> _out_horizon{"horizon", 2};

data SpaceGrid a = SpaceGrid
  { gridX :: !(Array a),
    gridY :: !(Array a),
    gridZ :: !(Array a)
  }
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

gridNumPoints :: AFType a => SpaceGrid a -> Int
gridNumPoints = AF.getElements . gridX

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

mkSpaceGrid :: AFType a => Array a -> Array a -> Array a -> SpaceGrid a
mkSpaceGrid gridX gridY gridZ = SpaceGrid x y z
  where
    (numX, _, _, _) = AF.getDims gridX
    (numY, _, _, _) = AF.getDims gridY
    (numZ, _, _, _) = AF.getDims gridZ
    !x = AF.flat $ AF.tile gridX [1, numY, numZ]
    !y = AF.flat $ AF.tile (AF.moddims gridY [1, numY]) [numX, 1, numZ]
    !z = AF.flat $ AF.tile (AF.moddims gridZ [1, 1, numZ]) [numX, numY, 1]

data Selector a = Selector
  { selectorMask :: !(Array CBool),
    selectorIndices :: !(Array Int)
  }
  deriving stock (Show, Generic)
  deriving anyclass (NFData)

bulkSelector :: (AFType a, Num a) => SpaceGrid a -> Selector a
bulkSelector (SpaceGrid _ _ z) = Selector isBulk bulkIndices
  where
    isBulk = AF.not $ AF.or (AF.eq z (scalar 0)) (AF.eq z (scalar 1))
    bulkIndices = AF.where' $ AF.cast isBulk

conformalBoundarySelector :: (AFType a, Num a) => SpaceGrid a -> Selector a
conformalBoundarySelector (SpaceGrid _ _ z) = Selector isBoundary boundaryIndices
  where
    isBoundary = AF.eq z (scalar 0)
    boundaryIndices = AF.where' $ AF.cast isBoundary

horizonBoundarySelector :: (AFType a, Num a) => SpaceGrid a -> Selector a
horizonBoundarySelector (SpaceGrid _ _ z) = Selector isBoundary boundaryIndices
  where
    isBoundary = AF.eq z (scalar 1)
    boundaryIndices = AF.where' $ AF.cast isBoundary

class IsSelectable f where
  select :: AFType a => Selector a -> f a -> f a

instance IsSelectable SpaceGrid where
  select (Selector _ indices) (SpaceGrid x y z) =
    SpaceGrid (AF.lookup x indices 0) (AF.lookup y indices 0) (AF.lookup z indices 0)

instance IsSelectable Fields where
  select (Selector _ indices) (Fields qs) = Fields (AF.lookup qs indices 0)

instance IsSelectable FieldsDerivatives where
  select (Selector _ indices) (FieldsDerivatives dqs ddqs) =
    FieldsDerivatives (AF.lookup dqs indices 0) (AF.lookup ddqs indices 0)

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

withSpaceGrid ::
  (AFType a, IsHalideType a) =>
  SpaceGrid a ->
  (Ptr HalideBuffer -> Ptr HalideBuffer -> Ptr HalideBuffer -> IO b) ->
  IO b
withSpaceGrid (SpaceGrid x y z) action =
  withHalideBuffer x $ \c_x ->
    withHalideBuffer y $ \c_y ->
      withHalideBuffer z $ \c_z ->
        action c_x c_y c_z

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
    r <- action $ ArrayView' fp shape (colMajorStrides shape)
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

instance NFData (Array a) where
  rnf (AF.Array (GHC.ForeignPtr !addr !contents)) = ()

importFields :: Text -> IO (SpaceGrid Double, Fields Double, FieldsDerivatives Double)
importFields filename = H5.withFile filename H5.ReadOnly $ \file -> do
  xGrid <- fromArrayView =<< H5.readDataset =<< H5.open file "x"
  yGrid <- fromArrayView =<< H5.readDataset =<< H5.open file "y"
  zGrid <- fromArrayView =<< H5.readDataset =<< H5.open file "z"
  let grid = mkSpaceGrid xGrid yGrid zGrid
  qs <- fromArrayView =<< H5.readDataset =<< H5.open file "Qs"
  dqs <- fromArrayView =<< H5.readDataset =<< H5.open file "DQs"
  ddqs <- fromArrayView =<< H5.readDataset =<< H5.open file "DDQs"
  let bulk = bulkSelector grid
  pure (select bulk grid, select bulk (Fields qs), select bulk (FieldsDerivatives dqs ddqs))

evalEquationsFromInput :: IO (Equations Double)
evalEquationsFromInput = do
  (grid, fields, dfields) <- importFields "../test_data.h5"
  let (!metric, !dmetric) = computeMetric 1.0 2.3 grid fields dfields
      !christoffel = computeChristoffel metric dmetric
      !maxwell = computeMaxwell grid fields dfields metric dmetric christoffel
      !deturck = computeDeTurck 1.0 2.3 grid fields dfields metric dmetric christoffel
      !eq = evaluateEquations 1.0 grid fields dfields metric dmetric christoffel maxwell deturck
  pure eq

importExpectedOutputs filename = H5.withFile filename H5.ReadOnly $ \file -> do
  (g_dd :: Array Double) <- fromArrayView =<< H5.readDataset =<< H5.open file "g_dd"
  (g_UU :: Array Double) <- fromArrayView =<< H5.readDataset =<< H5.open file "g_UU"
  (dg_ddd :: Array Double) <- fromArrayView =<< H5.readDataset =<< H5.open file "Dg_ddd"
  (dg_dUU :: Array Double) <- fromArrayView =<< H5.readDataset =<< H5.open file "Dg_dUU"
  (divf_d :: Array Double) <- fromArrayView =<< H5.readDataset =<< H5.open file "divF_d"
  (eq :: Array Double) <- fromArrayView =<< H5.readDataset =<< H5.open file "equations"
  pure (g_dd, g_UU, dg_ddd, dg_dUU, divf_d, eq)

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
