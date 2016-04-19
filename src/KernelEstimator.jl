VERSION >= v"0.4" && __precompile__()

module KernelEstimator 
using Distributions
using Optim
using Yeppp
using StatsBase
using Cubature
import StatsBase: RealVector, RealMatrix
using StatsFuns
import StatsFuns: invsqrt2π, log2π, sqrt2, invsqrt2
# Compat.@irrational invsqrt2π 0.398942280401432677939946 big(1.)/sqrt(big(2.)*π)



export  kde, kerneldensity, gaussiankernel, ekernel, gammakernel, betakernel,
        localconstant, locallinear, npr, lp0, lp1,
        bwkd, bwnormal, bwlscv, bwlcv, bwlocalconstant, bwlocallinear, bwreg,
#         Sind,
        bootstrapCB
#         BootstrapGoodness


include("kernel.jl")
include("regression.jl")
include("density.jl")
include("bandwidth.jl")
# include("semiparametric.jl")
include("confidenceband.jl")
# include("goodnessoffit.jl")


end
