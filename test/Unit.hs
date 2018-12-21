{-# LANGUAGE DataKinds #-}

module Unit where

import           Numeric.LinearAlgebra.Static
                   (L, R, Sq, ℝ, matrix, norm_2, unwrap, vector)
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
  (norm_2 (tr - expected)) <= 1.0e-3 @? "differs from expected"
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
sur = singleUseReliability q Nothing (successes, failures)

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
q' = reduced p

unit_reducedQ :: Assertion
unit_reducedQ = unwrap q' @?= unwrap expected
  where
    expected :: Sq 4
    expected = matrix
      [ 0, 1,    0,   0
      , 0, 0,    0.5, 0.5
      , 0, 0,    0.5, 0.25
      , 0, 0.25, 0,   0
      ]

-- Stimulus matrix.
s :: L 4 5
s = matrix
--         a     b    c     e     f
{- E -}  [ 1,    0,   0,    0,    0
{- A -}  , 0,    0.5, 0.5,  0,    0
{- B -}  , 0,    0.5, 0.25, 0.25, 0
{- C -}  , 0.25, 0,   0,    0.5,  0.25
         ]

unit_occurenceMean :: Assertion
unit_occurenceMean =
  (norm_2 ((occurrenceMean q') - expected)) <= 1.0e-3 @? "differs from expected"
  where
    expected :: R 4
    expected = vector
      [ 1.0, 1.231, 1.231, 0.9231 ]

unit_occurenceVar :: Assertion
unit_occurenceVar =
  (norm_2 ((occurrenceVar q') - expected)) <= 1.0e-3 @? "differs from expected"
  where
    expected :: R 4
    expected = vector
      [ 0, 0.284, 2.556, 0.497 ]
