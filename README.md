# Nonparametric
The Julia package for nonparametric density estimate and regression.

[![Build Status](https://travis-ci.org/panlanfeng/Nonparametric.jl.png)](https://travis-ci.org/panlanfeng/Nonparametric.jl)

## Functions
This package provides the following functions:	
 - `KernelDensity(xeval, xdata, kernel::Function, h)` do kernel density estimate  

 - `LP0(xeval, xdata, ydata::Vector, kernel::Function,h)` do local constant regression (or Nadaraya-Watson)  

 - `LP1(xeval, xdata, ydata::Vector, kernel,h)` do local linear regression  

 - `Sind(xdata::Matrix, ydata::Vector, kernel::Function)` fit single index model (which is semiparametric model actually)  

 - `BandwidthNormalReference(xdata::Vector)` select bandwidth for density estimate by rule of thumb  

 - `BandwidthLSCV(xdata, kernel::Function)` select bandwidth for density estiamte by leave-one-out  

 - `BandwidthLSCVReg(xdata, ydata::Vector, reg::Function, kernel::Function)` select bandwidth for regression using leave-one-out  

 - `BootstrapCB(B::Int64,xeval,xdata,ydata,testmodel::Function,reg,kernel,h)` compute 95% wild bootstrap confidence band at `xeval` for LP0  

 - `BootstrapGoodness(B::Int64, xdata, ydata, testmodel:Function, reg, kernel, h)`: do the goodness of fit test for `testmodel` based on wild bootstrap, where `testmodel` should be a function taking parameters like `testmodel(xdata, ydata::Vector)` and returning predicted values at `xdata`. For example, 
  
        function SemiPredict(xdata::Matrix{Float64}, ydata::Vector{Float64})
          xb=xdata*Sind(xdata,ydata,GaussianKernel)
          LP0(xb,xb,ydata,GaussianKernel,1.0)
        end
	  

In the above functions, 
 - `xeval` is the point(s) where the density or fitted value is calculated  

 - `xdata` is covariate(s), namely X, which can be one dimensional or multi-dimensional; should be of same dimension with `xeval` and `h`.   

 - `ydata` is the response vector y; should have same length as the rows of `xdata`  

 - `reg` is the regression function, can be `LP0` or `LP1`. `LP1` can only be used for univariate regression currently; default to be `LP0`  

 - `kernel` is the kernel function used. Only `GaussianKernel` is implemented; default to be `GaussianKernel`  

 - `h` is the bandwidth, can be a scalar or a vector, depending on whether `xdata` is univariate or multivariate; should have same length as the columns of `xeval` and `xdata`; default to be the value chosen by LSCV  

`xeval`, `xdata`, `ydata` and `h` should have element type Float64.  

##Demos

 - Kernel density estimate
        
        x=rand(100,4)
        y=x * ones(4) + x .^ 2 * ones(4)
        xeval=rand(50,4)
		den=KernelDensity(xeval, x)
        
 - Local constant regression 
       
        yfit0=LP0(xeval,x,y)

 - Confidence Band

         cb=BootstrapCB(100, xeval, x, y)

 - Goodness of fit test for some model (can be time consuming)

        function SemiPredict(xdata::Matrix{Float64}, ydata::Vector{Float64})
          xb=xdata*Sind(xdata,ydata,GaussianKernel)
          LP0(xb,xb,ydata,GaussianKernel,1.0)
        end

        BootstrapGoodness(10, x, y, SemiPredict) 

 - Univariate kernel density and regression

        using Distributions
        x=rand(Normal(), 100)
        xeval=linspace(minimum(x), maximum(x), 50)
        den=KernelDensity(xeval,x) 
              
        using Gadfly
        plot(layer(x=x,y=y, Geom.point), layer(x=xeval,y=den, Geom.line))

 - Univariate local constant regression and local linear regression
         
        y=2 .* x + rand(Normal(), 100)
        yfit0=LP0(xeval, x, y)
        yfit1=LP1(xeval, x, y)
        plot(layer(x=x, y=y, Geom.point), layer(x=xeval, y=yfit1,Geom.line))



##Reference
 - Lecture notes from Dr. Song Xi Chen  

 - W. Hardle and J. S. Marron (1991). Bootstrap Simultaneous Error Bars for Nonparametric Regression. _The Annals of Statistics_. Vol. 19, No. 2 (Jun., 1991), pp. 778-796  

 - W.Hardle and E. Mammen (1993). Comparing Nonparametric Versus Parametric Regression Fits. _The Annals of Statistics_. Vol. 21, No. 4 (Dec., 1993), pp. 1926-1947  





