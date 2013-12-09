# Nonparametric
The Julia package for nonparametric density estimate and regression.

-----------------------------------------

This package provides the follwoing functions:
 - `KernelDensity(xeval, xdata, kernel::Function, h)` do kernel density estimate
 - `LP0(xeval, xdata, ydata::Vector, kernel::Function,h)` do local constant regression (or Nadaraya-Watson)
 - `LP1(xeval, xdata, ydata::Vector, kernel,h)` do local linear regression
 - `Sind(xdata::Matrix, ydata::Vector, kernel::Function)` fit single index model (which is semiparametric model actually)
 - `BandwidthNormalReference(xdata::Vector, kernel::Function)` select bandwidth for density estimate by rule of thumb
 - `BandwidthLSCV(xdata, kernel::Function)` select bandwidth for density estiamte by leave-one-out
 - `BandwidthLSCVReg(xdata, ydata::Vector, reg::Function, kernel::Function)` select bandwidth for regression using leave-one-out
 - `BootstrapCB(B::Int64, xeval, xdata, ydata::Vector, reg::Function, kernel::Function, h)` compute 95% wild bootstrap confidence band at xeval for LP0 
 - `BootstrapGoodness(B::Int64, xdata, ydata::Vector, reg::Function, kernel::Function, h, testmodel::Function)`: do the goodness of fit test for `testmodel` based on wild bootstrap, where `testmodel` should take parameters like `testmodel(xdata, ydata::Vector)` and return predicted values at `xdata`. For example, 

    function SemiPredict(xdata::Matrix{Float64}, ydata::Vector{Float64})
      xb=xdata*Sind(xdata,ydata,GaussianKernel)
      LP0(xb,xb,y,GaussianKernel,1.0)
    end
  

In the above functions, 
 - `xeval` is the point where the density or fitted value is calculated; 
 - `xdata` is covariate(s), namely X, which can be one dimensional or multi-dimensional; 
 - `ydata` is the response vector y; 
 - `reg` is the regression function, can be LP0 or LP1. LP1 can only be used for univariate regression currently; 
 - `kernel` is the kernel function used. Only `GaussianKernel` is implemented; 
 - `h` is the bandwidth, can be a scalar or a vector, depending on whether `xdata` is univariate or multivariate.




