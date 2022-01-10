module EinsteinEquations.NewtonRaphson
  ( JacobianResult (..),
    -- FortranMatrix (..),
    numericalJacobian,
    -- nrm2,
    -- axpy,
    -- scal,
    -- solveDense,
    partialDerivative,
    RootOptions (..),
    RootResult (..),
    newtonRaphson,
  )
where

import Control.Monad.Primitive
import Control.Monad.ST
import Data.Complex
-- import Data.Vector.Storable (MVector (..), Vector)
-- import qualified Data.Vector.Storable as V
-- import qualified Data.Vector.Storable.Mutable as MV

import EinsteinEquations.Tensor
import Foreign.ForeignPtr
import Foreign.Marshal.Alloc
import Foreign.Marshal.Utils
import Foreign.Ptr
import Foreign.Storable
import GHC.Exts (IsList (..))
import GHC.Float
import System.IO.Unsafe (unsafePerformIO)
import Torch
import qualified Torch.Functional.Internal
import Prelude hiding (init)

-- data FortranMatrix a = FortranMatrix {fortranMatrixData :: !(ForeignPtr a), fortranMatrixShape :: !(Int, Int)}

-- instance Storable a => IsList (FortranMatrix a) where
--   type Item (FortranMatrix a) = [a]
--   fromList elements
--     | V.length v == n * m = FortranMatrix (fst (V.unsafeToForeignPtr0 v)) (n, m)
--     | otherwise = error "not a matrix"
--     where
--       n = length elements
--       m = case elements of
--         (x : _) -> length x
--         _ -> 0
--       v = V.fromList . concat . transpose $ elements
--   toList (FortranMatrix p (n, m)) =
--     transpose . chunksOf n . V.toList $ V.unsafeFromForeignPtr0 p (n * m)
--     where
--       chunksOf :: Int -> [a] -> [[a]]
--       chunksOf _ [] = []
--       chunksOf k l
--         | k > 0 = (take k l) : (chunksOf k (drop k l))
--         | otherwise = error $ "invalid k: " <> show k

data CSRMatrix = CSRMatrix
  deriving stock (Show)

data JacobianResult
  = DenseJacobian {-# UNPACK #-} !Tensor
  | SparseJacobian {-# UNPACK #-} !CSRMatrix
  deriving stock (Show)

-- | A GADT which allows to dispatch between BLAS types at runtime.
-- data BlasDatatypeTag :: (Type -> Type) where
--   FloatTag :: BlasDatatypeTag Float
--   DoubleTag :: BlasDatatypeTag Double
--   ComplexFloatTag :: BlasDatatypeTag (Complex Float)
--   ComplexDoubleTag :: BlasDatatypeTag (Complex Double)

-- deriving stock instance Show (BlasDatatypeTag a)

-- type family BlasRealPart a where
--   BlasRealPart Float = Float
--   BlasRealPart Double = Double
--   BlasRealPart (Complex Float) = Float
--   BlasRealPart (Complex Double) = Double

-- | BLAS datatype.
-- class
--   ( Storable a,
--     Storable (BlasRealPart a),
--     Floating a,
--     RealFloat (BlasRealPart a),
--     Show a,
--     Show (BlasRealPart a),
--     Eq (BlasRealPart a),
--     Ord (BlasRealPart a)
--   ) =>
--   BlasDatatype a
--   where
--   blasTag :: proxy a -> BlasDatatypeTag a

-- instance BlasDatatype Float where blasTag _ = FloatTag

-- instance BlasDatatype Double where blasTag _ = DoubleTag

-- instance BlasDatatype (Complex Float) where blasTag _ = ComplexFloatTag

-- instance BlasDatatype (Complex Double) where blasTag _ = ComplexDoubleTag

-- type BlasInt = Int32

-- extern void dgesv( int* n, int* nrhs, double* a, int* lda, int* ipiv,
--                    double* b, int* ldb, int* info );

-- type Cgesv a =
--   Ptr BlasInt ->
--   Ptr BlasInt ->
--   Ptr a ->
--   Ptr BlasInt ->
--   Ptr BlasInt ->
--   Ptr a ->
--   Ptr BlasInt ->
--   Ptr BlasInt ->
--   IO ()

-- foreign import ccall "sgesv_" sgesv :: Cgesv Float

-- foreign import ccall "dgesv_" dgesv :: Cgesv Double

-- solveDense :: forall a. BlasDatatype a => FortranMatrix a -> Vector a -> Vector a
-- solveDense (FortranMatrix p (n, m)) b
--   | n /= m = error $ "matrix A has wrong shape: " <> show (n, m) <> "; expected a square matrix"
--   | n < 0 = error $ "matrix A has invalid shape: " <> show (n, m)
--   | n /= V.length b = error $ "dimensions of A and b do not match: " <> show n <> " != " <> show (V.length b)
--   | otherwise = unsafePerformIO $ do
--     x <- V.thaw b
--     ipiv <- MV.new (max 1 n)
--     info <- with (fromIntegral n) $ \nPtr ->
--       with 1 $ \nrhsPtr ->
--         withForeignPtr p $ \aPtr ->
--           with (fromIntegral $ max 1 n) $ \ldaPtr ->
--             MV.unsafeWith ipiv $ \ipivPtr ->
--               MV.unsafeWith x $ \xPtr ->
--                 with (fromIntegral $ max 1 n) $ \ldxPtr ->
--                   alloca $ \infoPtr ->
--                     gesv nPtr nrhsPtr aPtr ldaPtr ipivPtr xPtr ldxPtr infoPtr >> peek infoPtr
--     case compare info 0 of
--       LT -> error $ "parameter №" <> show (-info) <> " (counting from 1) is invalid"
--       GT -> error $ "matrix A is singular: U" <> show (info - 1, info - 1) <> " is exactly zero"
--       EQ -> V.unsafeFreeze x
--   where
--     gesv :: Cgesv a
--     gesv = case blasTag (Proxy @a) of
--       DoubleTag -> dgesv
--       FloatTag -> sgesv
--       tag -> error $ "GESV is not (yet) implemented for " <> show tag
-- {-# NOINLINE solveDense #-}

invertJacobian :: JacobianResult -> Tensor -> Tensor
invertJacobian (DenseJacobian j) δq = linalgSolve j δq
invertJacobian _ _ = error "sparse Jacobians are not yet supported"

-- type Cnrm2 a =
--   Ptr BlasInt ->
--   Ptr a ->
--   Ptr BlasInt ->
--   IO (BlasRealPart a)

-- foreign import ccall unsafe "snrm2_" snrm2 :: Cnrm2 Float

-- foreign import ccall unsafe "dnrm2_" dnrm2 :: Cnrm2 Double

-- nrm2 :: forall a. BlasDatatype a => Vector a -> BlasRealPart a
-- nrm2 x
--   | V.length x == 0 = 0
--   | otherwise = unsafePerformIO $ do
--     with (fromIntegral (V.length x)) $ \nPtr ->
--       V.unsafeWith x $ \xPtr ->
--         with 1 $ \incxPtr ->
--           go nPtr xPtr incxPtr
--   where
--     go :: Cnrm2 a
--     go = case blasTag (Proxy @a) of
--       FloatTag -> snrm2
--       DoubleTag -> dnrm2
--       tag -> error $ "NRM2 is not (yet) implemented for " <> show tag
-- {-# NOINLINE nrm2 #-}

-- type Caxpy a =
--   Ptr BlasInt ->
--   Ptr a ->
--   Ptr a ->
--   Ptr BlasInt ->
--   Ptr a ->
--   Ptr BlasInt ->
--   IO ()

-- foreign import ccall unsafe "saxpy_" saxpy :: Caxpy Float

-- foreign import ccall unsafe "daxpy_" daxpy :: Caxpy Double

-- axpy :: forall a m. (BlasDatatype a, PrimMonad m) => a -> Vector a -> MVector (PrimState m) a -> m ()
-- axpy α x y@(MVector _ p)
--   | V.length x /= MV.length y =
--     error $
--       "vectors x and y have different lengths: " <> show (V.length x) <> " != " <> show (MV.length y)
--   | otherwise = unsafeIOToPrim $
--     with (fromIntegral (V.length x)) $ \nPtr ->
--       with α $ \αPtr ->
--         V.unsafeWith x $ \xPtr ->
--           with 1 $ \incxPtr ->
--             withForeignPtr p $ \yPtr ->
--               with 1 $ \incyPtr ->
--                 go nPtr αPtr xPtr incxPtr yPtr incyPtr
--   where
--     go :: Caxpy a
--     go = case blasTag (Proxy @a) of
--       FloatTag -> saxpy
--       DoubleTag -> daxpy
--       tag -> error $ "AXPY is not (yet) implemented for " <> show tag
-- {-# NOINLINE axpy #-}

-- type Cscal a =
--   Ptr BlasInt ->
--   Ptr a ->
--   Ptr a ->
--   Ptr BlasInt ->
--   IO ()

-- foreign import ccall unsafe "sscal_" sscal :: Cscal Float

-- foreign import ccall unsafe "dscal_" dscal :: Cscal Double

-- scal :: forall a m. (BlasDatatype a, PrimMonad m) => a -> MVector (PrimState m) a -> m ()
-- scal α (MVector n p) = unsafeIOToPrim $
--   with (fromIntegral n) $ \nPtr ->
--     with α $ \αPtr ->
--       withForeignPtr p $ \xPtr ->
--         with 1 $ \incxPtr ->
--           go nPtr αPtr xPtr incxPtr
--   where
--     go :: Cscal a
--     go = case blasTag (Proxy @a) of
--       FloatTag -> sscal
--       DoubleTag -> dscal
--       tag -> error $ "SCAL is not (yet) implemented for " <> show tag
-- {-# NOINLINE scal #-}

partialDerivative ::
  forall m.
  (Monad m) =>
  (Tensor -> m Tensor) ->
  Maybe Tensor ->
  Tensor ->
  Int ->
  m Tensor
partialDerivative f _y₀ x₀ i = do
  y₀ <- case _y₀ of
    Just y -> return y
    Nothing -> f x₀
  y₁ <- f x₁
  -- Prelude.trace ("y₀ = " <> show y₀) $ pure ()
  -- Prelude.trace ("y₁ = " <> show y₁) $ pure ()
  pure $! (1 / ε) * (y₁ - y₀)
  where
    δx = maskedFill (zerosLike x₀) i ε
    x₁ = x₀ + δx
    ε = scale * precision
    scale = if norm > 0 then 2 ^^ (exponent norm) else 1
    norm = asValue (Torch.Functional.Internal.normCastAll x₀ 2 Double) :: Double
    precision = case dtype x₀ of
      Float -> 0.00048828125
      Double -> 2.9802322387695312e-8
      ty -> error $ "partialDerivative is not implemented for " <> show ty

numericalJacobian ::
  forall m.
  (Monad m) =>
  (Tensor -> m Tensor) ->
  Maybe Tensor ->
  Tensor ->
  m JacobianResult
numericalJacobian f _y₀ x₀ = do
  y₀ <- case _y₀ of
    Just y -> return y
    Nothing -> f x₀
  derivatives <- forM [0 .. numel x₀ - 1] $ partialDerivative f (Just y₀) x₀
  return . DenseJacobian $ stack (Dim 1) derivatives

data LineSearchOptions

data RootOptions = RootOptions
  { rootOptionsCriterion :: !(Double -> Double -> Bool),
    rootOptionsMaxIter :: !Int,
    rootOptionsLineSearch :: !(Maybe LineSearchOptions)
  }

data RootState = RootState
  { rootStateCurrent :: Tensor,
    rootStateValue :: Tensor,
    rootStateResidual :: Double
  }
  deriving stock (Show)

data RootResult = RootResult
  { rootResultCurrent :: Tensor,
    rootResultHistory :: [Double]
  }
  deriving stock (Show)

newtonRaphsonStep ::
  forall m.
  (Monad m) =>
  (Tensor -> m Tensor) ->
  (Tensor -> m JacobianResult) ->
  RootState ->
  m RootState
newtonRaphsonStep f df (RootState x₀ y₀ _) = do
  jacobian <- df x₀
  let !δx = invertJacobian jacobian y₀
      !x = x₀ - δx
  y <- f x
  -- Prelude.trace ("newtonRaphsonStep: x₀ = " <> show x₀ <> ", y₀ = " <> show y₀ <> ", jacobian = " <> show jacobian <> ", δx = " <> show δx <> ", x = " <> show x <> ", y = " <> show y) $ pure ()
  let r = asValue $ Torch.Functional.Internal.normCastAll y 2 Double
  return $ RootState x y r

newtonRaphson ::
  forall m.
  (Monad m) =>
  RootOptions ->
  (Tensor -> m Tensor) ->
  (Tensor -> m JacobianResult) ->
  Tensor ->
  m RootResult
newtonRaphson options f df x₀ = do
  y₀ <- f x₀
  let iₘₐₓ = rootOptionsMaxIter options
      shouldStop = rootOptionsCriterion options r₀ . rootStateResidual
      go acc !i !s
        | i >= iₘₐₓ || shouldStop s = return (rootStateCurrent s, acc)
        | otherwise = do
          s' <- newtonRaphsonStep f df s
          let !h = rootStateResidual s'
          go (acc ++ [h]) (i + 1) s'
      r₀ = asValue (Torch.Functional.Internal.normCastAll y₀ 2 Double) :: Double
      init = RootState x₀ y₀ r₀
  (solution, history) <- go [r₀] 0 init
  return $ RootResult solution (fromList history)
