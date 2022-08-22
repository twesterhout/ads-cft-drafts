{-# LANGUAGE ScopedTypeVariables #-}

module Metric.IDRSpec (spec) where

import ArrayFire (AFType, Array, scalar)
import qualified ArrayFire as AF
import qualified Data.HDF5 as H5
import Metric
import Metric.IDR
import Test.Hspec

allclose :: (AFType a, Num a) => Array a -> Array a -> a -> a -> Bool
allclose a b rtol atol =
  (== 1) . fst . AF.allTrueAll $
    AF.abs (a - b) `AF.le` (scalar atol + scalar rtol * AF.abs b)

testMatrix1 :: (Array Double, Array Double, Array Double, Array Double)
testMatrix1 = (m, c, x, p)
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
    p =
      AF.matrix (4, 8) $
        [ [0.18915290913774319, 0.6439047071509303, 0.33343270739238606, 0.5781210730058675],
          [0.0139400430554667, 0.1744577914801967, 0.8018687143461121, 0.7068193846816763],
          [0.5566447844028717, 0.5152617812512884, 0.2216499363643485, 0.9509623456487204],
          [0.2432638218373302, 0.3424892271425941, 0.5678265568504112, 0.579763352328501],
          [0.1644860312447063, 0.2663956206720568, 0.5961616465156822, 0.22386058610554715],
          [0.1964898778275529, 0.5932308749874771, 0.6269099275997124, 0.7974440963227224],
          [0.4769595479104116, 0.318453698862897, 0.5924848883479175, 0.5229227191543668],
          [0.40020774939254444, 0.5232379234182464, 0.4293031583551137, 0.5432088622777586]
        ]

spec :: Spec
spec = do
  describe "idrs'" $ do
    it "solves a small linear system" $ do
      let (m, b, x, p) = testMatrix1
          apply v = pure $ AF.matmul m v AF.None AF.None
          solution = AF.vector 4 [0.39736576702218857, 0.2181228116472028, -1.2736601900991644, 2.179106537095057]
      r1 <- idrs' defaultIDRParams p apply b x :: IO (IDRResult Double)
      r1 `shouldSatisfy` idrResultConverged
      idrResultX r1 `shouldSatisfy` (\v -> allclose v solution 1.0e-7 1.0e-8)
    it "solves a bigger linear system" $ do
      let n = 100
          s = 16
      -- m <- (-) <$> AF.randu [n, n] <*> pure 0.5
      -- b <- (-) <$> AF.randu [n] <*> pure 0.5
      -- x <- (-) <$> AF.randu [n] <*> pure 0.5
      -- p <- AF.randu [n, s]
      -- H5.withFile "idr_inputs.h5" H5.WriteTruncate $ \file -> do
      --   H5.createDataset file "m" m
      --   H5.createDataset file "b" b
      --   H5.createDataset file "x" x
      --   H5.createDataset file "p" p
      (m, b, x, p :: Array Double) <-
        H5.withFile "idr_inputs.h5" H5.ReadOnly $ \file -> do
          m <- H5.open file "m" >>= H5.readDataset
          b <- H5.open file "b" >>= H5.readDataset
          x <- H5.open file "x" >>= H5.readDataset
          p <- H5.open file "p" >>= H5.readDataset
          pure (AF.transpose m False, b, x, AF.transpose p False)
      -- print m
      let apply v = pure $ AF.matmul m v AF.None AF.None
          solution = AF.solve m b AF.None
      r <- idrs (defaultIDRParams {idrParamsS = s, idrParamsMaxIters = 2 * n}) apply b x :: IO (IDRResult Double)
      print (idrResultIter r)
      print (idrResultResidual r)
      print . AF.maxAll . AF.abs $ idrResultX r - solution
      r `shouldSatisfy` idrResultConverged
      idrResultX r `shouldSatisfy` (\v -> allclose v solution 1.0e-6 1.0e-8)
