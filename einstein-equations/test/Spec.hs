module Main (main) where

import Control.Foldl (FoldM (..))
import Data.Primitive.Types (Prim)
import Data.Vector.Storable (Storable, Vector)
import qualified Data.Vector.Storable as V
import EinsteinEquations.Collocation
import EinsteinEquations.NewtonRaphson
-- import EinsteinEquations.Tensor
import GHC.Exts (IsList (..))
import qualified Streamly.Prelude as S
import Test.Hspec
import Torch (Tensor, allclose, asTensor, asValue, (!))
import qualified Torch
import Prelude hiding (toList)

roundTo :: RealFrac a => a -> Int -> a
roundTo x n = (fromInteger . round $ x * (10 ^ n)) / (10.0 ^^ n)

shouldApproxEqual ::
  (Storable a, Show a, RealFrac a) =>
  Vector a ->
  Vector a ->
  IO ()
shouldApproxEqual x y = V.zipWithM_ (\a b -> roundTo a n `shouldBe` roundTo b n) x y
  where
    n = 7

-- tensorShouldEqual ::
--   (SupportedRank r, Prim a, Show a, RealFrac a) =>
--   Tensor 'CPU r a ->
--   Tensor 'CPU r a ->
--   IO ()
-- tensorShouldEqual a b = do
--   (tensorShape a) `shouldBe` (tensorShape b)
--   blockFoldM (FoldM step initial extract) (tensorShape a)
--   where
--     n = 10
--     eq x y = roundTo x n `shouldBe` roundTo y n
--     initial = pure ()
--     extract _ = pure ()
--     step () i = do
--       a' <- read a i
--       b' <- read b i
--       eq a' b'

main :: IO ()
main = hspec $ do
  -- describe "nrm2" $ do
  --   it "computes L₂ norm of a vector" $ do
  --     let x :: (Storable a, Floating a) => Vector a
  --         x =
  --           fromList
  --             [ -0.4495389460471282,
  --               0.1507946864613462,
  --               -0.2582198045631925,
  --               0.19783960251972765,
  --               -0.08502606567957727,
  --               -0.027089038952040734,
  --               -0.45259419418301516,
  --               -0.34808112568663063,
  --               -0.2492391171933619,
  --               -0.03652065590737674
  --             ]
  --     nrm2 (x :: Vector Double) `shouldBe` 0.8532651379629458
  --     nrm2 (x :: Vector Float) `shouldBe` 0.8532651
  -- describe "solveDense" $ do
  --   it "solves linear systems of equations" $ do
  --     let a₁ =
  --           fromList
  --             [ [-0.12569571, -0.2457441, -0.0917917],
  --               [0.4542093, 0.23042252, -0.2675304],
  --               [0.26157499, -0.38497729, -0.2527988]
  --             ]
  --         b₁ = fromList [(0.38405567 :: Double), 0.05975013, 0.31492095]
  --         x₁ = solveDense a₁ b₁
  --     x₁ `shouldBe` (fromList [-1.1729225766772948, -0.10262275813807731, -2.303099262911118])
  --     let a₂ =
  --           fromList
  --             [ [0.13702166080474854, -0.19305512309074402],
  --               [0.036100149154663086, 0.4381260871887207]
  --             ]
  --         b₂ = fromList [-0.3771229088306427 :: Float, 0.23093461990356445]
  --         x₂ = solveDense a₂ b₂
  --     x₂ `shouldBe` (fromList [-1.8006048202514648, 0.675460159778595])
  describe "partialDerivative" $ do
    it "computes partial derivatives numerically" $ do
      let f :: Monad m => Tensor -> m Tensor
          f v =
            let x = asValue (v ! (0 :: Int)) :: Float
                y = asValue (v ! (1 :: Int)) :: Float
             in pure . asTensor $ [sin (x + y), cos (x - y)]
          x₀ = 0.1 :: Float
          y₀ = 0.25 :: Float
      df₀ <- partialDerivative f Nothing (asTensor [x₀, y₀]) 0
      print $ df₀
      print $ asTensor [cos (x₀ + y₀), -sin (x₀ - y₀)]
      Torch.allclose
        df₀
        (asTensor [cos (x₀ + y₀), -sin (x₀ - y₀)])
        1.0e-3
        1.0e-5
        False
        `shouldBe` True
      df₁ <- partialDerivative f Nothing (asTensor [x₀, y₀]) 1
      print $ df₁
      print $ asTensor [cos (x₀ + y₀), sin (x₀ - y₀)]
      Torch.allclose
        df₁
        (asTensor [cos (x₀ + y₀), sin (x₀ - y₀)])
        2.0e-3
        1.0e-5
        False
        `shouldBe` True
  describe "newtonRaphson" $ do
    it "solves tanh(x) == 0" $ do
      let f :: Monad m => Tensor -> m Tensor
          f = pure . Torch.tanh
          df :: Monad m => Tensor -> m JacobianResult
          df v = pure . DenseJacobian $ asTensor [[1 / cosh x ^^ (2 :: Int)]]
            where
              x = asValue (v ! (0 :: Int)) :: Float
      (RootResult x _) <-
        newtonRaphson
          (RootOptions (\_ r -> r < 1.0e-8) 10 Nothing)
          f
          df
          (asTensor [0.5 :: Float])
      allclose x (asTensor [0.0 :: Float]) 1.0e-5 1.0e-7 False `shouldBe` True
  it "solves x (x - 2) (x + 2) == 0 using numerical differentiation" $ do
    let f :: Tensor -> IO Tensor
        f x = pure $ x ^^ 3 - 4 * x
        df :: Tensor -> IO JacobianResult
        df x = do
          r <- numericalJacobian f Nothing x
          -- print r
          pure r
    -- let x = v V.! 0
    --  in return . DenseJacobian $ fromList [[3 * x ^^ 2 - 4]]
    (RootResult x _) <-
      newtonRaphson
        (RootOptions (\_ r -> r < 1.0e-8) 10 Nothing)
        f
        df
        (asTensor [15 :: Float])
    print x
    allclose x (asTensor [2.0 :: Float]) 1.0e-5 1.0e-7 False `shouldBe` True
  -- describe "tensorToList" $ do
  --   it "converts Tensor to List" $ do
  --     let v = V.generate 12 id
  --         shape = fromList [2, 3]
  --         strides = fromList [6, 2]
  --         (t :: Tensor 'CPU 2 Int) = Tensor (fst $ V.unsafeToForeignPtr0 v) shape strides
  --     toList t `shouldBe` [[0, 2, 4], [6, 8, 10]]
  --   it "converts List to Tensor" $ do
  --     let v = V.generate 8 id
  --         (t :: Tensor 'CPU 2 Int) = fromList [[0, 1, 2, 3], [4, 5, 6, 7]]
  --     V.unsafeFromForeignPtr0 (tensorData t) (tensorLength t) `shouldBe` v
  describe "differentiationMatrixPeriodic" $ do
    it "constructs differentiation matrix for periodic functions" $ do
      print $ differentiationMatrixBounded 0 1 4
      True `shouldBe` True
  -- describe "reversePrimArray" $ do
  --   it ".." $ do
  --     reversePrimArray (fromList [1, 2, 3 :: Int]) `shouldBe` fromList [3, 2, 1]
  --     reversePrimArray (fromList [1 :: Int]) `shouldBe` fromList [1]
  --     reversePrimArray (fromList ([] :: [Int])) `shouldBe` fromList []
  -- describe "rowMajorStrides" $ do
  --   it ".." $ do
  --     rowMajorStrides 1 (fromList []) `shouldBe` fromList []
  --     rowMajorStrides 2 (fromList [3]) `shouldBe` fromList [2]
  --     rowMajorStrides 1 (fromList [3, 2]) `shouldBe` fromList [2, 1]
  --     rowMajorStrides 1 (fromList [1, 3, 2]) `shouldBe` fromList [6, 2, 1]
  --     rowMajorStrides 1 (fromList [2, 1, 3, 2]) `shouldBe` fromList [6, 6, 2, 1]
  -- describe "reshape" $ do
  --   it ".." $ do
  --     let (t :: Tensor 'CPU 2 Int) = fromList [[1, 2, 3], [4, 5, 6]]
  --         t' = reshape @2 [3, 2] t
  --     (toList <$> t') `shouldBe` (Just [[1, 2], [3, 4], [5, 6]])
  describe "differentiateX" $
    it "computes derivative w.r.t. the first coordinate" $ do
      let n = 6
          xs = gridPointsForPeriodic (pi) n
          dX = differentiationMatrixPeriodic n
          f = Torch.reshape [n, 1, 1, 1] $ Torch.sin (2 * xs)
          -- reshape [n, 1, 1, 1] $ tensorMap (\x -> sin (2 * x)) xs
          dF = differentiateX dX f
          dFExpected = Torch.reshape [n, 1, 1, 1] $ Torch.cos (2 * xs)
      print $ dX
      print $ dFExpected
      print $ dF
      Torch.allclose dF dFExpected 1.0e-10 1.0e-12 False `shouldBe` True

-- t = generate1 n (\i -> sin (
