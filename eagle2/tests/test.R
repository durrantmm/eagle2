
require(eagle2)

N=20
K=3
Ti=2

ys=array( rpois( N*K*Ti, 20 ), c(N,Ti,K) )
ns=ys + array( rpois( N*K*Ti, 20 ), c(N,Ti,K) )

bb( ys, ns )

require(rstan)
require(abind)
require(foreach)

concShape=1.0001
concRate=1e-4
allow_homozyg=F
