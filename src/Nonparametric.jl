require("Distributions")
require("Cubature")
require("Optim")
require("NumericExtensions")

module Nonparametric

using Distributions
using Cubature
using Optim
using NumericExtensions








export  GaussianKernel, KernelDensity, EKernel, Gaussian, Epanechnikov, KernelType,
        LP0, LP1, 
        BandwidthNormalReference, BandwidthLSCV, BandwidthLSCVReg, 
        Sind,
        BootstrapCB, 
        BootstrapGoodness


include("kernel.jl")
include("regression.jl")
include("density.jl")
include("bandwidth.jl")
include("semiparametric.jl")
include("confidenceband.jl")
include("goodnessoffit.jl")


end
