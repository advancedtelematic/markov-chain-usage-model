{-# LANGUAGE DataKinds #-}

module Unit where

import           GHC.TypeLits
import           Numeric.LinearAlgebra.Static
import           Test.Tasty.HUnit

import           MarkovChain

------------------------------------------------------------------------

p :: Sq 5
p = matrix
  [ 0, 1,    0,   0,    0
  , 0, 0,    0.5, 0.5,  0
  , 0, 0,    0.5, 0.25, 0.25
  , 0, 0.25, 0,   0,    0.75
  , 1, 0,    0,   0,    0
  ]

q :: Sq 4
q = reduced p

unit_reduced :: Assertion
unit_reduced = unwrap q @?= unwrap expected
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

