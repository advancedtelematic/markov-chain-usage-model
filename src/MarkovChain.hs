{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeOperators       #-}

module MarkovChain where

import           Data.Maybe
                   (fromMaybe)
import           GHC.TypeLits
import           Numeric.LinearAlgebra.Static
import           Prelude                      hiding
                   ((<>))

------------------------------------------------------------------------

reduced :: (KnownNat m, KnownNat (m - p), KnownNat (m - (m - p)))
        => (KnownNat n, KnownNat (n - q), KnownNat (n - (n - q)))
        => m - p <= m
        => n - q <= n
        => L m n -> L (m - p) (n - q)
reduced = fst . splitCols . fst . splitRows

fundamental :: KnownNat n => Sq n -> Sq n
fundamental q = inv (eye - q)

variance :: KnownNat n => Sq n -> Sq n
variance n = n <> (2 * diag (takeDiag n) - eye) - (n * n)

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
    n1 = occurrenceMean (reduced p) & 1

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
  -> ℝ
singleUseReliability q mprior obs =
  fst $ headTail $ uncol $ transientReliability q mprior obs
