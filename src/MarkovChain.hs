{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeOperators       #-}

module MarkovChain where

import           GHC.TypeLits
import           Numeric.LinearAlgebra.Static
import           Prelude                      hiding
                   ((<>))

------------------------------------------------------------------------

reduce :: (KnownNat m, KnownNat (m - p), KnownNat (m - (m - p)))
       => (KnownNat n, KnownNat (n - q), KnownNat (n - (n - q)))
       => m - p <= m
       => n - q <= n
       => L m n -> L (m - p) (n - q)
reduce = fst . splitCols . fst . splitRows

restrict :: (KnownNat m, KnownNat (m - p), KnownNat p)
         => (KnownNat n, KnownNat (n - q), KnownNat q)
         => p <= m
         => q <= n
         => L m n -> L (m - p) (n - q)
restrict = snd . splitCols . snd . splitRows

fundamental :: KnownNat n => Sq n -> Sq n
fundamental q = inv (eye - q)

variance :: KnownNat n => Sq n -> Sq n
variance n = n <> (2 <> diag (takeDiag n) - eye) - (n * n)

pi :: KnownNat n => Sq n -> Sq n -> R n
pi p n = undefined
