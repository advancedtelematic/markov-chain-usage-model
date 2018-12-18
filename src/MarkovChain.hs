{-# LANGUAGE DataKinds        #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE GADTs            #-}
{-# LANGUAGE TypeOperators    #-}

module MarkovChain where

import           Prelude hiding ((<>))
import           GHC.TypeLits
import           Numeric.LinearAlgebra.Static

------------------------------------------------------------------------

reduced :: (KnownNat n, KnownNat (n - 1), KnownNat (n - (n - 1)))
        => n - 1 <= n
        => Sq n -> Sq (n - 1)
reduced = fst . splitCols . fst . splitRows

restriction :: Sq 5 -> Sq 5
restriction = undefined

fundamental :: KnownNat n => Sq n -> Sq n
fundamental q = inv (eye - q)

variance :: KnownNat n => Sq n -> Sq n
variance n = n <> (2 <> diag (takeDiag n) - eye) - (n * n)

pi :: KnownNat n => Sq n -> Sq n -> R n
pi p n = undefined
