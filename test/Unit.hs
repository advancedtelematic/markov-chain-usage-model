{-# LANGUAGE DataKinds #-}

module Unit where

import           GHC.TypeLits
import           Numeric.LinearAlgebra.Static
import           Prelude                      hiding ((<>))
import           Test.Tasty.HUnit

import           MarkovChain

------------------------------------------------------------------------

-- cf. A Simpler and More Direct Derivation of System Reliability Using Markov Chain Usage Models (2017)
-- by Lan Lin, Yufeng Xue and Fengguang Song

-- Usage model.
p :: Sq 3
p = matrix
  [ 0, 0.5, 0.5
  , 0, 0.5, 0.5
  , 0, 0,   0
  ]

-- Transient usage model.
q :: L 2 3
q = reduce p

qdot :: Sq 2
qdot = reduce p

-- Expected transient reliabilities.
r1 :: L 2 3
r1 = matrix
  [ 0.1, 0.2, 0.3
  , 0.4, 0.5, 0.6
  ]

r1dot :: Sq 2
r1dot = reduce r1

-- Observed transient reliabilities.
fancyR1 :: L 2 3
fancyR1 = q * r1

fancyR1dot :: Sq 2
fancyR1dot = qdot * r1dot

-- Last column vector of fancyR.
w :: L 2 1
w = restrict fancyR1

unit_expectW :: Assertion
unit_expectW = unwrap w @?= unwrap expected
  where
    expected :: L 2 1
    expected = matrix
      [ 0.3/2
      , 0.6/2
      ]

-- Success rate from transient state to sink.
rstar :: L 2 1
rstar = inv (eye - fancyR1dot) <> w

unit_expectRstar :: Assertion
unit_expectRstar = unwrap rstar @?= unwrap expected
  where
    a = 0.3/2
    b = 0.2*0.6/4
    c = 0.6/2
    d = 1 - 0.5/2

    expected :: L 2 1
    expected = matrix
      [ a + b/d
      , c/d
      ]

sur :: Double
sur = fst $ headTail (uncol rstar)

------------------------------------------------------------------------

-- cf. Computations for Markov Chain Usage Models (2000)
-- by S. J. Prowell

-- p :: Sq 5
-- p = matrix
--   [ 0, 1,    0,   0,    0
--   , 0, 0,    0.5, 0.5,  0
--   , 0, 0,    0.5, 0.25, 0.25
--   , 0, 0.25, 0,   0,    0.75
--   , 1, 0,    0,   0,    0
--   ]

-- Stimulus matrix.
-- s :: L 4 5
-- s = matrix
-- --         a     b    c     e     f
-- {- E -}  [ 1,    0,   0,    0,    0
-- {- A -}  , 0,    0.5, 0.5,  0,    0
-- {- B -}  , 0,    0.5, 0.25, 0.25, 0
-- {- C -}  , 0.25, 0,   0,    0.5,  0.25
--          ]

-- unit_reduced :: Assertion
-- unit_reduced = unwrap q @?= unwrap expected
--   where
--     expected :: Sq 4
--     expected = matrix
--       [ 0, 1,    0,   0
--       , 0, 0,    0.5, 0.5
--       , 0, 0,    0.5, 0.25
--       , 0, 0.25, 0,   0
--       ]
