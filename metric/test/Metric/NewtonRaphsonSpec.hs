module Metric.NewtonRaphsonSpec (spec) where

import ArrayFire (AFType, Array, scalar)
import qualified ArrayFire as AF
import Metric.NewtonRaphson
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
  describe "newtonRaphson" $ do
    it "solves a system of 2 non-linear equations" $ do
      let f v =
            let [a, b] = AF.toList v
             in pure $ AF.vector 2 [-a * a + a + 0.5 - b, a * a - 5 * a * b - b]
          solution = AF.vector 2 [1.233317793, 0.2122450145]
          v₀ = AF.vector 2 [1.2, 1.2]
          ε = 1.0e-6
          opts = RootOptions (\_ r -> r < ε) 10 Nothing
      (RootResult r _) <- newtonRaphson opts f v₀
      r `shouldSatisfy` (\v -> allclose v solution ε 1.0e-8)
    it "solves a system of 3 non-linear equations" $ do
      let f v =
            let [x, y, z] = AF.toList v
             in pure $
                  AF.vector
                    3
                    [ 9 * x * x + 36 * y * y + 4 * z * z - 36,
                      x * x - 2 * y * y - 20 * z,
                      x * x - y * y + z * z
                    ]
          solution = AF.vector 3 [0.893628, 0.894527, -0.0400893]
          v₀ = AF.vector 3 [1.0, 1.0, 0.0]
          ε = 1.0e-6
          opts = RootOptions (\_ r -> r < ε) 10 Nothing
      (RootResult r _) <- newtonRaphson opts f v₀
      r `shouldSatisfy` (\v -> allclose v solution ε 1.0e-8)
