{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE ExplicitNamespaces  #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeOperators       #-}

module MarkovChain
  ( reduceRow
  , reduceCol
  , unsafeTransform
  , fundamental
  , variance
  , expectedLength
  , onFundamental
  , occurrenceMean
  , occurrenceVar
  , pi
  , successRate
  , transientReliability
  , singleUseReliability
  , singleUseReliabilityIO
  , loadStaticMatrix
  , load2StaticMatrices
  , saveStaticMatrix
  , kullbackLeibler
  )
  where

import           Control.Applicative
                   (liftA2)
import           Control.Monad
                   (join, liftM2)
import           Data.Maybe
                   (fromJust, fromMaybe)
import           GHC.TypeLits
                   (type (+), type (-), type (<=), KnownNat)
import qualified Numeric.LinearAlgebra.Data   as D
import           Numeric.LinearAlgebra.Static
import           Prelude                      hiding
                   (pi, (<>))

------------------------------------------------------------------------

{- $setup
   >>> :set -XDataKinds
   >>> :{
   let
     p :: Sq 5
     p = matrix
       [ 0, 1,    0,   0,    0
       , 0, 0,    0.5, 0.5,  0
       , 0, 0,    0.5, 0.25, 0.25
       , 0, 0.25, 0,   0,    0.75
       , 1, 0,    0,   0,    0
       ]
:}
-}

-- |
-- >>> reduceRow p :: L 4 5
-- (matrix
--  [ 0.0,  1.0, 0.0,  0.0,  0.0
--  , 0.0,  0.0, 0.5,  0.5,  0.0
--  , 0.0,  0.0, 0.5, 0.25, 0.25
--  , 0.0, 0.25, 0.0,  0.0, 0.75 ] :: L 4 5)
reduceRow :: forall m n. (KnownNat m, KnownNat (m - 1))
          => (KnownNat n, m - 1 <= m)

          => L m n        -- ^ Matrix to reduce.
          -> L (m - 1) n
reduceRow = fst . splitRows

-- |
-- >>> reduceCol p :: L 5 4
-- (matrix
--  [ 0.0,  1.0, 0.0,  0.0
--  , 0.0,  0.0, 0.5,  0.5
--  , 0.0,  0.0, 0.5, 0.25
--  , 0.0, 0.25, 0.0,  0.0
--  , 1.0,  0.0, 0.0,  0.0 ] :: L 5 4)
reduceCol :: forall m n. (KnownNat m, KnownNat n)
          => (KnownNat (n - 1), KnownNat (n - (n - 1)), n - 1 <= n)

          => L m n        -- ^ Matrix to reduce.
          -> L m (n - 1)
reduceCol = fst . splitCols

unsafeTransform :: (Sized t2 a d2, Sized t1 b d1)
                => (d2 t2 -> d1 t1) -> a -> b
unsafeTransform f = fromJust . create . f . unwrap

------------------------------------------------------------------------

-- | The fundamental matrix for absorbing chains. Its (i, j)-th entry is
-- the expected number of occurrences of state j prior to absorption at
-- the sink, given that one starts in state i. So the first row
-- indicates the expected occurence of each state starting from the
-- start state.
--
-- >>> fundamental (reduceCol (reduceRow p) :: Sq 4)
-- (matrix
--  [ 1.0,  1.2307692307692306, 1.2307692307692308, 0.9230769230769231
--  , 0.0,  1.2307692307692306, 1.2307692307692308, 0.9230769230769231
--  , 0.0, 0.15384615384615385, 2.1538461538461537, 0.6153846153846154
--  , 0.0,  0.3076923076923077, 0.3076923076923077, 1.2307692307692308 ] :: L 4 4)
fundamental :: KnownNat n
            => Sq n -- ^ Reduced matrix.
            -> Sq n
fundamental q = inv (eye - q)

-- | Expected variance of the occurrence for each state. The (i, j)-th
-- entry should be read in the same away as that of the fundamental
-- matrix above.
--
-- >>> variance (fundamental (reduceCol (reduceRow p) :: Sq 4))
-- (matrix
--  [ 0.0, 0.28402366863905315, 2.5562130177514795, 0.4970414201183433
--  , 0.0, 0.28402366863905315, 2.5562130177514795, 0.4970414201183433
--  , 0.0, 0.20118343195266267, 2.4852071005917153, 0.5207100591715977
--  , 0.0,  0.3550295857988165, 0.9230769230769231, 0.2840236686390534 ] :: L 4 4)
variance :: KnownNat n
         => Sq n        -- ^ Reduced matrix.
         -> Sq n
variance n = n <> (2 * diag (takeDiag n) - eye) - (n * n)

-- | Expected test case length.
--
-- >>> expectedLength (fundamental (reduceCol (reduceRow p) :: Sq 4))
-- 4.384615384615385
expectedLength :: KnownNat n => Sq n -- ^ Fundamental matrix.
               -> Double
expectedLength n = sumV (toRows n !! 0)
  where
    sumV :: KnownNat n => R n -> Double
    sumV v = v <.> 1

------------------------------------------------------------------------

onFundamental :: (KnownNat n, 1 <= n)
              => (Sq n -> Sq n)
              -> Sq n
              -> R n
onFundamental f = unrow . fst . splitRows . f . fundamental

occurrenceMean, occurrenceVar :: (KnownNat n, 1 <= n)
                              => Sq n -> R n
occurrenceMean = (id `onFundamental`)
occurrenceVar  = (variance `onFundamental`)

-- | The Perron eigenvector (long-run occupancies/probabilities of states).
--
-- >>> unwrap (pi p)
-- [0.1857142857142857,0.22857142857142854,0.22857142857142856,0.17142857142857143,0.1857142857142857]
--
-- __note__: The use of 'unwrap' is needed because hmatrix's pretty printing of
-- vectors seems broken.
pi :: forall n. (KnownNat n, KnownNat (n - 1), KnownNat (n - (n - 1)))
   => (1 <= n - 1, n - 1 <= n, ((n - 1) + 1) ~ n)
   => Sq n -- ^ Transition matrix.
   -> R n
pi p = n1 / linspace (l, l)
  where
    n1 :: R n
    n1 = occurrenceMean (reduceCol $ reduceRow p) & 1

    l :: ℝ
    l  = takeDiag eye <.> n1

------------------------------------------------------------------------

-- |
--
-- >>> successRate Nothing (matrix [10,10,10,10] :: Sq 2, matrix [0,1,2,3])
-- (matrix
--  [ 0.9166666666666666, 0.8461538461538461
--  , 0.7857142857142857, 0.7333333333333333 ] :: L 2 2)
--
-- >>> :{
--   successRate (Just (matrix [10,10,10,10] :: Sq 2, matrix [1,1,1,1]))
--               (matrix [10,10,10,10], matrix [0,1,2,3])
-- :}
-- (matrix
--  [ 0.9523809523809523, 0.9090909090909091
--  , 0.8695652173913043, 0.8333333333333334 ] :: L 2 2)
successRate
  :: forall m n
   . (KnownNat m, KnownNat n)
  => Maybe (L m n, L m n)
  -> (L m n, L m n)
  -> L m n
successRate mprior (obsSuccs, obsFails) = succs / (succs + fails)
  where
    priorSuccs, priorFails :: L m n
    (priorSuccs, priorFails) =
      fromMaybe (build (const . const 1), build (const . const 1)) mprior

    succs = priorSuccs + obsSuccs
    fails = priorFails + obsFails

transientReliability :: forall n. (KnownNat (n - 1), KnownNat n)
  => (((n - (n - 1)) ~ 1), (n - 1) <= n)

  => L (n - 1) n                      -- ^ Reduced transition matrix.
  -> Maybe (L (n - 1) n, L (n - 1) n) -- ^ Prior successes and failures.
  -> (L (n - 1) n, L (n - 1) n)       -- ^ Observed successes and failures.
  -> L (n - 1) 1
transientReliability q mprior obs = rstar
  where
    r :: L (n - 1) n
    r = successRate mprior obs

    fancyR :: L (n - 1) n
    fancyR = q * r

    fancyRdot :: L (n - 1) (n - 1)
    w         :: L (n - 1) 1
    (fancyRdot, w) = splitCols fancyR

    rstar :: L (n - 1) 1
    rstar = fundamental fancyRdot <> w

singleUseReliability :: (KnownNat (n - 1), KnownNat n)
  => (1 <= (n - 1), (n - (n - 1)) ~ 1, (n - 1) <= n)

  => L (n - 1) n                       -- ^ Reduced transition matrix.
  -> Maybe (L (n - 1) n, L (n - 1) n)  -- ^ Prior successes and failures.
  -> (L (n - 1) n, L (n - 1) n)        -- ^ Observed successes and failures
  -> Double
singleUseReliability q mprior obs =
  fst $ headTail $ uncol $ transientReliability q mprior obs

singleUseReliabilityIO :: (KnownNat n, KnownNat (n - 1))
  => (1 <= n - 1, n - 1 <= n, (n - (n - 1)) ~ 1)

  => L (n - 1) n                -- ^ Transition matrix without transitions to
                                -- the sink.

  -> FilePath                   -- ^ Filepath where prior successes are/will be
                                -- saved.

  -> FilePath                   -- ^ Filepath where prior failures are/will be
                                -- saved.

  -> (L (n - 1) n, L (n - 1) n) -- ^ Observed successes and failures.
  -> IO Double
singleUseReliabilityIO q fps fpf (s, f) = do
  mprior <- load2StaticMatrices fps fpf
  case mprior of
    Nothing -> do
      -- We need the elementwise @max 1@ below because otherwise we get division
      -- by zero in 'successRate'. For details see the 2017 paper by by L. Lin,
      -- Y. Xue and F. Song in the README.
      saveStaticMatrix fps "%.1f" (dmmap (max 1) s)
      saveStaticMatrix fpf "%.1f" (dmmap (max 1) f)
    Just (ps, pf) -> do
      saveStaticMatrix fps "%.1f" (ps + s)
      saveStaticMatrix fpf "%.1f" (pf + f)
  return (singleUseReliability q mprior (s, f))

------------------------------------------------------------------------

loadStaticMatrix :: (KnownNat m, KnownNat n) => FilePath -> IO (Maybe (L m n))
loadStaticMatrix = fmap (join . fmap create) . D.loadMatrix'

load2StaticMatrices :: (KnownNat m, KnownNat n, KnownNat o, KnownNat p)
                    => FilePath -> FilePath -> IO (Maybe (L m n, L o p))
load2StaticMatrices fp1 fp2 = liftMA2 (loadStaticMatrix fp1) (loadStaticMatrix fp2)
  where
    liftMA2 :: (Monad m, Applicative f) => m (f a) -> m (f b) -> m (f (a, b))
    liftMA2 = liftM2 (liftA2 (,))

saveStaticMatrix :: (KnownNat m, KnownNat n)
                 => FilePath
                 -> String   -- ^ "printf" format (e.g. "%.2f", "%g", etc.)
                 -> L m n -> IO ()
saveStaticMatrix fp format = D.saveMatrix fp format . unwrap

------------------------------------------------------------------------

kullbackLeibler :: forall m n. (KnownNat m, KnownNat n)
  => (KnownNat (n - 1), KnownNat (n - (n - 1)))
  => ((n - 1) <= n, m <= n, 1 <= (n - 1), ((n - 1) + 1) ~ n)

  => Sq n -> L m n -> R n -> L m n -> ℝ
kullbackLeibler p s es et = k' <.> linspace (1, 1)
  where
    pi' :: R m
    pi' = fst $ split (pi p)

    es' :: R m
    es' = fst (split es)

    m :: Int
    m = size es

    es'' :: L m n
    es'' = unsafeTransform (\mat -> D.repmat mat 1 m) (col es')

    t :: L m n
    t = es'' / et

    k :: R m
    k  = vector $ map sum' $ toRows $ s * dmmap log' (s * t) / log 2

    k' :: R m
    k' = pi' * k

    log' :: ℝ -> ℝ
    log' x
      | x == 0    = 0
      | isNaN x   = 0
      | otherwise = log x

    sum' :: R n -> ℝ
    sum' = (<.> linspace (1, 1))
