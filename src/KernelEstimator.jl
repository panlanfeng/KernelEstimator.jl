# require("Distributions")
# require("Cubature")
# require("Optim")
# require("NumericExtensions")
# require("StatsBase")


module KernelEstimator 
using Compat
using Distributions
# using Cubature
using Optim
#using Yeppp
using StatsBase
using Cubature
import StatsBase: RealVector, RealMatrix

 Compat.@irrational invsqrt2pi 0.398942280401432677939946 big(1.)/sqrt(big(2.)*π)



export  kde, kerneldensity, gaussiankernel, ekernel,gammakernel, betakernel,
        lp0, lp1, npr,
        bwkd, bwnormal, bwlscv, bwlcv, bwlp0, bwlp1,bwreg,
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
