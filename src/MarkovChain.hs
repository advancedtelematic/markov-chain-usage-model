{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeOperators       #-}

module MarkovChain where

import           Control.Monad (liftM2, join)
import           Control.Applicative (liftA2)
import           Data.Maybe
                   (fromJust, fromMaybe)
import           GHC.TypeLits
import qualified Numeric.LinearAlgebra.Data   as D
import           Numeric.LinearAlgebra.Static
import           Prelude                      hiding
                   (pi, (<>))

------------------------------------------------------------------------

-- $setup
-- >>> :set -XDataKinds
-- >>> :{
-- let
--   p :: Sq 5
--   p = matrix
--     [ 0, 1,    0,   0,    0
--     , 0, 0,    0.5, 0.5,  0
--     , 0, 0,    0.5, 0.25, 0.25
--     , 0, 0.25, 0,   0,    0.75
--     , 1, 0,    0,   0,    0
--     ]
-- :}

-- |
-- >>> reduceRow p :: L 4 5
-- (matrix
--  [ 0.0,  1.0, 0.0,  0.0,  0.0
--  , 0.0,  0.0, 0.5,  0.5,  0.0
--  , 0.0,  0.0, 0.5, 0.25, 0.25
--  , 0.0, 0.25, 0.0,  0.0, 0.75 ] :: L 4 5)
reduceRow
  :: forall m n
  .  (KnownNat m, KnownNat (m - 1), KnownNat (m - (m - 1)))
  => KnownNat n
  => m - 1 <= m
  => L m n -> L (m - 1) n
reduceRow = fst . splitRows

-- |
-- >>> reduceCol p :: L 5 4
-- (matrix
--  [ 0.0,  1.0, 0.0,  0.0
--  , 0.0,  0.0, 0.5,  0.5
--  , 0.0,  0.0, 0.5, 0.25
--  , 0.0, 0.25, 0.0,  0.0
--  , 1.0,  0.0, 0.0,  0.0 ] :: L 5 4)
reduceCol
  :: forall m n p
  .  KnownNat m
  => (KnownNat n, KnownNat (n - 1), KnownNat (n - (n - 1)))
  => n - 1 <= n
  => L m n -> L m (n - 1)
reduceCol = fst . splitCols

unsafeTransform
  :: (Sized t2 a d2, Sized t1 b d1) => (d2 t2 -> d1 t1) -> a -> b
unsafeTransform f = fromJust . create . f . unwrap

------------------------------------------------------------------------

fundamental :: KnownNat n => Sq n -- ^ Reduced matrix.
            -> Sq n
fundamental q = inv (eye - q)

variance :: KnownNat n => Sq n -> Sq n
variance n = n <> (2 * diag (takeDiag n) - eye) - (n * n)

-- |
-- >>> expectedLength (fundamental (reduceCol (reduceRow p) :: Sq 4))
-- 4.384615384615385
expectedLength :: KnownNat n => Sq n -- ^ Fundamental matrix.
               -> Double
expectedLength n = sumV (toRows n !! 0)
  where
    sumV :: KnownNat n => R n -> Double
    sumV v = v <.> 1

------------------------------------------------------------------------

onFundamental
  :: KnownNat n
  => 1 <= n
  => (Sq n -> Sq n) -> Sq n
  -> R n
onFundamental f = unrow . fst . splitRows . f . fundamental

occurrenceMean, occurrenceVar
  :: KnownNat n
  => 1 <= n
  => Sq n -> R n
occurrenceMean = (id `onFundamental`)
occurrenceVar  = (variance `onFundamental`)

pi
  :: forall n
   . (KnownNat n)
  => KnownNat (n - 1)
  => KnownNat (n - (n - 1))
  => 1 <= n - 1
  => n - 1 <= n
  => ((n - 1) + 1) ~ n
  => Sq n -> R n
pi p = n1 / linspace (l, l)
  where
    n1 :: R n
    n1 = occurrenceMean (reduceCol $ reduceRow p) & 1

    l :: ℝ
    l  = takeDiag eye <.> n1

------------------------------------------------------------------------

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

transientReliability
  :: forall n
   . (KnownNat (n - 1), KnownNat n)
  => (KnownNat ((n - 1) - (n - 1)))
  => (n - (n - 1)) ~ 1
  => (n - 1) <= n
  => L (n - 1) n
  -> Maybe (L (n - 1) n, L (n - 1) n)
  -> (L (n - 1) n, L (n - 1) n)
  -> L (n - 1) 1
transientReliability q mprior obs = rstar
  where
    r :: L (n - 1) n
    r = successRate mprior obs

    fancyR :: L (n - 1) n
    fancyR = q * r

    fancyRdot :: L (n - 1) (n - 1)
    w :: L (n - 1) 1
    (fancyRdot, w) = splitCols fancyR

    rstar :: L (n - 1) 1
    rstar = fundamental fancyRdot <> w

singleUseReliability
  :: (KnownNat (n - 1), KnownNat n)
  => (KnownNat ((n - 1) - (n - 1)))
  => 1 <= (n - 1)
  => (n - (n - 1)) ~ 1
  => (n - 1) <= n

  => L (n - 1) n
  -> Maybe (L (n - 1) n, L (n - 1) n)
  -> (L (n - 1) n, L (n - 1) n)
  -> Double
singleUseReliability q mprior obs =
  fst $ headTail $ uncol $ transientReliability q mprior obs

singleUseReliabilityIO :: (KnownNat n, KnownNat (n - 1), KnownNat ((n - 1 - (n - 1))))
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
  mprior <- load2StaticMatrix fps fpf
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

load2StaticMatrix :: (KnownNat m, KnownNat n, KnownNat o, KnownNat p)
                  => FilePath -> FilePath -> IO (Maybe (L m n, L o p))
load2StaticMatrix fp1 fp2 = liftMA2 (loadStaticMatrix fp1) (loadStaticMatrix fp2)
  where
    liftMA2 :: (Monad m, Applicative f) => m (f a) -> m (f b) -> m (f (a, b))
    liftMA2 = liftM2 (liftA2 (,))

saveStaticMatrix :: (KnownNat m, KnownNat n)
                 => FilePath
                 -> String   -- ^ "printf" format (e.g. "%.2f", "%g", etc.)
                 -> L m n -> IO ()
saveStaticMatrix fp format = D.saveMatrix fp format . unwrap

------------------------------------------------------------------------

kullbackLeibler
  :: forall m n
   . (KnownNat m, KnownNat n)
  => (KnownNat (n - 1), KnownNat (n - (n - 1)))
  => (n - 1) <= n
  => (m <= n)
  => (1 <= (n - 1))
  => ((n - 1) + 1) ~ n
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
