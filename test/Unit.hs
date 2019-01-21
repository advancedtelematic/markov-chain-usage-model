{-# LANGUAGE DataKinds  #-}
{-# LANGUAGE LambdaCase #-}

module Unit where

import           Data.Proxy
                   (Proxy(Proxy))
import           Numeric.LinearAlgebra.Static
                   (L, R, Sq, ℝ, matrix, norm_2, unwrap, vector)
import           Prelude                      hiding
                   (pi, (<>))
import           Test.DocTest
                   (doctest)
import           Test.Tasty.HUnit
                   (Assertion, (@?), (@?=))

import           MarkovChain

------------------------------------------------------------------------

-- cf. A Simpler and More Direct Derivation of System Reliability Using Markov Chain Usage Models (2017)
-- by Lan Lin, Yufeng Xue and Fengguang Song

-- Transient usage model (transitions from sink omitted).
q :: L 2 3
q = matrix
  [ 0, 0.5, 0.5
  , 0, 0.5, 0.5
  ]

successes :: L 2 3
successes = matrix
  [ 1, 2, 3
  , 4, 5, 6
  ]

failures :: L 2 3
failures = matrix
  [ 1, 0, 1
  , 0, 1, 0
  ]

-- Observed transient success rate.
r :: L 2 3
r = successRate Nothing (successes, failures)

unit_expectSuccessRate :: Assertion
unit_expectSuccessRate = unwrap r @?= unwrap expected
  where
    expected :: L 2 3
    expected = matrix
      [ 2/4, 3/4, 4/6
      , 5/6, 6/8, 7/8
      ]

-- Transient reliability matrix.
tr :: L 2 1
tr = transientReliability q Nothing (successes, failures)

unit_expectTransientReliability :: Assertion
unit_expectTransientReliability =
  (norm_2 (tr - expected)) <= 1.0e-6 @? "differs from expected"
  where
    a = (4/6)/2
    b = (3/4)*(7/8)/4
    c = (7/8)/2
    d = 1 - (6/8)/2

    expected :: L 2 1
    expected = matrix
      [ a + b/d
      , c/d
      ]

-- Single use reliability mean.
sur :: ℝ
sur = singleUseReliability Proxy q Nothing (successes, failures)

------------------------------------------------------------------------

-- cf. Computations for Markov Chain Usage Models (2000)
-- by S. J. Prowell

-- Usage model.
p :: Sq 5
p = matrix
  [ 0, 1,    0,   0,    0
  , 0, 0,    0.5, 0.5,  0
  , 0, 0,    0.5, 0.25, 0.25
  , 0, 0.25, 0,   0,    0.75
  , 1, 0,    0,   0,    0
  ]

q' :: Sq 4
q' = reduceCol $ reduceRow p

-- Stimulus matrix.
s :: L 4 5
s = matrix
--         a     b    c     e     f
{- E -}  [ 1,    0,   0,    0,    0
{- A -}  , 0,    0.5, 0.5,  0,    0
{- B -}  , 0,    0.5, 0.25, 0.25, 0
{- C -}  , 0.25, 0,   0,    0.5,  0.25
         ]

-- State occurences.
unit_occurenceMean :: Assertion
unit_occurenceMean =
  norm_2 (occurrenceMean q' - expected) <= 1.0e-3 @? "differs from expected"
  where
    expected :: R 4
    expected = vector
      [ 1.0, 1.231, 1.231, 0.9231 ]

unit_occurenceVar :: Assertion
unit_occurenceVar =
  norm_2 (occurrenceVar q' - expected) <= 1.0e-3 @? "differs from expected"
  where
    expected :: R 4
    expected = vector
      [ 0, 0.284, 2.556, 0.497 ]

-- Long-run occupancies.
unit_longRunOccupancies :: Assertion
unit_longRunOccupancies =
  norm_2 (pi p - expected) <= 1.0e-4 @? "differs from expected"
  where
    expected :: R 5
    expected = vector
      [ 0.1857, 0.2286, 0.2286, 0.1714, 0.1857 ]

-- Stimulus execution.
et :: L 4 5
et = matrix
  [ 4, 0, 0, 0, 0
  , 0, 3, 2, 0, 0
  , 0, 4, 1, 2, 0
  , 1, 0, 0, 1, 1
  ]

-- Vector of state visitations.
es :: R 5
es = vector [4, 5, 7, 3, 4]

unit_kullbackLeibler :: Assertion
unit_kullbackLeibler =
  kullbackLeibler p s es et - expected <= 1.0e-6 @? "differs from expected"
  where
    expected :: ℝ
    expected = 0.03441

unit_docTest :: IO ()
unit_docTest = doctest ["src/MarkovChain.hs"]
