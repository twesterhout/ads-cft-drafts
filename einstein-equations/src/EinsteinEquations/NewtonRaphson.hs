module EinsteinEquations.NewtonRaphson
  ( JacobianResult (..),
    FortranMatrix (..),
    numericalJacobian,
    nrm2,
    axpy,
    scal,
    solveDense,
    partialDerivative,
  )
where

import Control.Monad.Primitive
import Data.Complex
import Data.Vector.Storable (MVector (..), Vector)
import qualified Data.Vector.Storable as V
import qualified Data.Vector.Storable.Mutable as MV
import Foreign.ForeignPtr
import Foreign.Marshal.Alloc
import Foreign.Marshal.Utils
import Foreign.Ptr
import Foreign.Storable
import GHC.Exts (IsList (..))
import GHC.Float
import System.IO.Unsafe (unsafePerformIO)

data FortranMatrix a = FortranMatrix {fortranMatrixData :: !(ForeignPtr a), fortranMatrixShape :: !(Int, Int)}

instance Storable a => IsList (FortranMatrix a) where
  type Item (FortranMatrix a) = [a]
  fromList elements
    | V.length v == n * m = FortranMatrix (fst (V.unsafeToForeignPtr0 v)) (n, m)
    | otherwise = error "not a matrix"
    where
      n = length elements
      m = case elements of
        (x : _) -> length x
        _ -> 0
      v = V.fromList . concat . transpose $ elements
  toList (FortranMatrix p (n, m)) =
    transpose . chunksOf n . V.toList $ V.unsafeFromForeignPtr0 p (n * m)
    where
      chunksOf :: Int -> [a] -> [[a]]
      chunksOf _ [] = []
      chunksOf k l
        | k > 0 = (take k l) : (chunksOf k (drop k l))
        | otherwise = error $ "invalid k: " <> show k

data CSRMatrix a = CSRMatrix

data JacobianResult a
  = DenseJacobian {-# UNPACK #-} !(FortranMatrix a)
  | SparseJacobian {-# UNPACK #-} !(CSRMatrix a)

-- | A GADT which allows to dispatch between BLAS types at runtime.
data BlasDatatypeTag :: (Type -> Type) where
  FloatTag :: BlasDatatypeTag Float
  DoubleTag :: BlasDatatypeTag Double
  ComplexFloatTag :: BlasDatatypeTag (Complex Float)
  ComplexDoubleTag :: BlasDatatypeTag (Complex Double)

deriving stock instance Show (BlasDatatypeTag a)

type family BlasRealPart a where
  BlasRealPart Float = Float
  BlasRealPart Double = Double
  BlasRealPart (Complex Float) = Float
  BlasRealPart (Complex Double) = Double

-- | BLAS datatype.
class
  ( Storable a,
    Storable (BlasRealPart a),
    Floating a,
    RealFloat (BlasRealPart a),
    Show a,
    Show (BlasRealPart a),
    Eq (BlasRealPart a),
    Ord (BlasRealPart a)
  ) =>
  BlasDatatype a
  where
  blasTag :: proxy a -> BlasDatatypeTag a

instance BlasDatatype Float where blasTag _ = FloatTag

instance BlasDatatype Double where blasTag _ = DoubleTag

instance BlasDatatype (Complex Float) where blasTag _ = ComplexFloatTag

instance BlasDatatype (Complex Double) where blasTag _ = ComplexDoubleTag

type BlasInt = Int32

-- extern void dgesv( int* n, int* nrhs, double* a, int* lda, int* ipiv,
--                    double* b, int* ldb, int* info );

type Cgesv a =
  Ptr BlasInt ->
  Ptr BlasInt ->
  Ptr a ->
  Ptr BlasInt ->
  Ptr BlasInt ->
  Ptr a ->
  Ptr BlasInt ->
  Ptr BlasInt ->
  IO ()

foreign import ccall "sgesv_" sgesv :: Cgesv Float

foreign import ccall "dgesv_" dgesv :: Cgesv Double

solveDense :: forall a. BlasDatatype a => FortranMatrix a -> Vector a -> Vector a
solveDense (FortranMatrix p (n, m)) b
  | n /= m = error $ "matrix A has wrong shape: " <> show (n, m) <> "; expected a square matrix"
  | n < 0 = error $ "matrix A has invalid shape: " <> show (n, m)
  | n /= V.length b = error $ "dimensions of A and b do not match: " <> show n <> " != " <> show (V.length b)
  | otherwise = unsafePerformIO $ do
    x <- V.thaw b
    ipiv <- MV.new (max 1 n)
    info <- with (fromIntegral n) $ \nPtr ->
      with 1 $ \nrhsPtr ->
        withForeignPtr p $ \aPtr ->
          with (fromIntegral $ max 1 n) $ \ldaPtr ->
            MV.unsafeWith ipiv $ \ipivPtr ->
              MV.unsafeWith x $ \xPtr ->
                with (fromIntegral $ max 1 n) $ \ldxPtr ->
                  alloca $ \infoPtr ->
                    gesv nPtr nrhsPtr aPtr ldaPtr ipivPtr xPtr ldxPtr infoPtr >> peek infoPtr
    case compare info 0 of
      LT -> error $ "parameter №" <> show (- info) <> " (counting from 1) is invalid"
      GT -> error $ "matrix A is singular: U" <> show (info - 1, info - 1) <> " is exactly zero"
      EQ -> V.unsafeFreeze x
  where
    gesv :: Cgesv a
    gesv = case blasTag (Proxy @a) of
      DoubleTag -> dgesv
      FloatTag -> sgesv
      tag -> error $ "GESV is not (yet) implemented for " <> show tag
{-# NOINLINE solveDense #-}

invertJacobian :: forall a. BlasDatatype a => JacobianResult a -> Vector a -> Vector a
invertJacobian = undefined

type Cnrm2 a =
  Ptr BlasInt ->
  Ptr a ->
  Ptr BlasInt ->
  IO (BlasRealPart a)

foreign import ccall unsafe "snrm2_" snrm2 :: Cnrm2 Float

foreign import ccall unsafe "dnrm2_" dnrm2 :: Cnrm2 Double

nrm2 :: forall a. BlasDatatype a => Vector a -> BlasRealPart a
nrm2 x
  | V.length x == 0 = 0
  | otherwise = unsafePerformIO $ do
    with (fromIntegral (V.length x)) $ \nPtr ->
      V.unsafeWith x $ \xPtr ->
        with 1 $ \incxPtr ->
          go nPtr xPtr incxPtr
  where
    go :: Cnrm2 a
    go = case blasTag (Proxy @a) of
      FloatTag -> snrm2
      DoubleTag -> dnrm2
      tag -> error $ "NRM2 is not (yet) implemented for " <> show tag
{-# NOINLINE nrm2 #-}

type Caxpy a =
  Ptr BlasInt ->
  Ptr a ->
  Ptr a ->
  Ptr BlasInt ->
  Ptr a ->
  Ptr BlasInt ->
  IO ()

foreign import ccall unsafe "saxpy_" saxpy :: Caxpy Float

foreign import ccall unsafe "daxpy_" daxpy :: Caxpy Double

axpy :: forall a m. (BlasDatatype a, PrimMonad m) => a -> Vector a -> MVector (PrimState m) a -> m ()
axpy α x y@(MVector _ p)
  | V.length x /= MV.length y =
    error $
      "vectors x and y have different lengths: " <> show (V.length x) <> " != " <> show (MV.length y)
  | otherwise = unsafeIOToPrim $
    with (fromIntegral (V.length x)) $ \nPtr ->
      with α $ \αPtr ->
        V.unsafeWith x $ \xPtr ->
          with 1 $ \incxPtr ->
            withForeignPtr p $ \yPtr ->
              with 1 $ \incyPtr ->
                go nPtr αPtr xPtr incxPtr yPtr incyPtr
  where
    go :: Caxpy a
    go = case blasTag (Proxy @a) of
      FloatTag -> saxpy
      DoubleTag -> daxpy
      tag -> error $ "AXPY is not (yet) implemented for " <> show tag
{-# NOINLINE axpy #-}

type Cscal a =
  Ptr BlasInt ->
  Ptr a ->
  Ptr a ->
  Ptr BlasInt ->
  IO ()

foreign import ccall unsafe "sscal_" sscal :: Cscal Float

foreign import ccall unsafe "dscal_" dscal :: Cscal Double

scal :: forall a m. (BlasDatatype a, PrimMonad m) => a -> MVector (PrimState m) a -> m ()
scal α (MVector n p) = unsafeIOToPrim $
  with (fromIntegral n) $ \nPtr ->
    with α $ \αPtr ->
      withForeignPtr p $ \xPtr ->
        with 1 $ \incxPtr ->
          go nPtr αPtr xPtr incxPtr
  where
    go :: Cscal a
    go = case blasTag (Proxy @a) of
      FloatTag -> sscal
      DoubleTag -> dscal
      tag -> error $ "SCAL is not (yet) implemented for " <> show tag
{-# NOINLINE scal #-}

partialDerivative ::
  forall a m.
  (BlasDatatype a, PrimMonad m) =>
  (Vector a -> m (Vector a)) ->
  Maybe (Vector a) ->
  Vector a ->
  Int ->
  m (Vector a)
partialDerivative f _f₀ x₀ i = do
  f₀ <- case _f₀ of
    Just y -> return y
    Nothing -> f x₀
  x₁ <- do
    t <- V.thaw x₀
    MV.modify t (+ ε) i
    V.unsafeFreeze t
  f₁ <- V.unsafeThaw =<< f x₁
  axpy (-1) f₀ f₁
  scal (1 / ε) f₁
  V.unsafeFreeze f₁
  where
    xNorm = nrm2 x₀
    scale =
      if xNorm > 0
        then 2 ^^ (exponent xNorm)
        else 1
    precision = case blasTag (Proxy @a) of
      FloatTag -> 0.00048828125
      DoubleTag -> 2.9802322387695312e-8
      tag -> error $ "partialDerivative is not implemented for " <> show tag
    ε = scale * precision

numericalJacobian ::
  forall a m.
  (BlasDatatype a, PrimMonad m) =>
  (Vector a -> m (Vector a)) ->
  Maybe (Vector a) ->
  Vector a ->
  m (JacobianResult a)
numericalJacobian f _f₀ x₀ = do
  f₀ <- case _f₀ of
    Just y -> return y
    Nothing -> f x₀
  derivatives <- forM [0 .. V.length x₀ - 1] $ partialDerivative f (Just f₀) x₀
  return . DenseJacobian $
    FortranMatrix
      (fst . V.unsafeToForeignPtr0 . V.concat $ derivatives)
      (V.length f₀, V.length x₀)

data LineSearchOptions a

data RootOptions a = RootOptions
  { rootOptionsCriterion :: !(a -> a -> Bool),
    rootOptionsMaxIter :: !Int,
    rootOptionsLineSearch :: !(Maybe (LineSearchOptions a))
  }

data RootState a = RootState
  { rootStateCurrent :: Vector a,
    rootStateValue :: Vector a,
    rootStateResidual :: BlasRealPart a,
  }

newtonRaphsonStep ::
  forall a m.
  (BlasDatatype a, PrimMonad m) =>
  (Vector a -> m (Vector a)) ->
  (Vector a -> m (JacobianResult a)) ->
  RootState a
  -> m (RootState a)
newtonRaphsonStep f df (RootState x₀ y₀ r₀) = do
  jacobian <- df x₀
  let δx = invertJacobian jacobian y₀
  x <- do
    t <- MV.thaw x₀
    axpy (-1) δx t
    V.unsafeFreeze t
  y <- f x
  return $ RootState x y (nrm2 y)

newtonRaphson ::
  forall a m.
  (BlasDatatype a, PrimMonad m) =>
  RootOptions a ->
  (Vector a -> m (Vector a)) ->
  (Vector a -> m (JacobianResult a)) ->
  Vector a ->
  m (Vector a)
newtonRaphson options f df _x₀ = do
  let iₘₐₓ = rootOptionsMaxIter options
      go x₀ r₀ i
        | i >= iₘₐₓ = return x₀
        | otherwise = do
          undefined
  f₀ <- f _x₀
  let r₀ = nrm2 f₀
  go _x₀
