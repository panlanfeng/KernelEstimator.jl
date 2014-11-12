# Nonparametric
The Julia package for nonparametric density estimate and regression. Currently includes univariate kernel density estimate, local constant regression (Nadaraya-watson estimator) and local linear regression. The Bootstrap confidence band [4] is provided for both of the two regression functions. The previous support for multivariate density estimation and regression are temporarily removed because of the low efficiency.

This package provides Gamma and Beta kernel to deal with bounded density estimation and regression problem. These two kernels are free of boundary bias for one side bounded and two sides bounded problems respectively, see [2, 3]. In particular, least square cross validation (LSCV) bandwidth selection functions are implemented. The bandwidth selection for these kernels is hard. Numeric integration is used here so it may be slow for large datasets.

Bandwidth selection is critical in kernel estimation. But no Julia packages provide reliable bandwidth selection function (up to my knowledge on Nov 12, 2014). LSCV is always recommended for kernel density estimation. Likelihood cross validation is provided but should be avoided because of known drawbacks. For regression problem, the bandwidth of local constant regression is selected using LSCV while that for local linear regression is chosen by AIC [6].


[![Build Status](https://travis-ci.org/panlanfeng/Nonparametric.jl.png)](https://travis-ci.org/panlanfeng/Nonparametric.jl)

## Functions
This package provides the following functions:

 - `kde(xdata::RealVector, xeval::Real or RealVector, kernel::Function, h::Real)`,  kernel density estimation

 - `kerneldensity(xdata::RealVector; xeval::RealVector=xdata, lb::Real=-Inf, ub::Real=Inf, kernel::Function=gaussiankernel, h::Real=-Inf)`, the easy function for kernel density estimation

 - `lp0(xdata::RealVector, ydata::RealVector; xeval::RealVector=xdata, kernel::Function=gaussiankernel, h::Real=bwlp0(xdata,ydata,kernel))`, local constant regression (or Nadaraya-Watson)

 - `lp1(xdata::RealVector, ydata::RealVector; xeval::RealVector=xdata, kernel::Function=gaussiankernel, h::Real=bwlp0(xdata,ydata,kernel))`,  local linear regression

 - `npr(xdata::RealVector, ydata::RealVector; xeval::RealVector=xdata,reg::Function=LP1, lb::Real=-Inf, ub::Real=Inf, kernel::Function=gaussiankernel, h::Real=-Inf) `, the easy function for regression

 - `bwnormal(xdata::Vector)`, bandwidth selection for density estimate by rule of thumb

 - `bwlscv(xdata::RealVector, kernel::Function)`, bandwidthselection for density estiamte by least square cross validation

  - `bwlcv(xdata::RealVector, kernel::Function)`, bandwidth selection for density estiamte by likelihood cross validation

 - `bwlp0(xdata, ydata::Vector, kernel)` select bandwidth for local constant regression using LSCV

 - `bwlp1(xdata, ydata::Vector, kernel)` select bandwidth for local linear regression using corrected AIC. See reference [6].

 - `bootstrapCB(xdata::RealVector, ydata::RealVector; xeval::RealVector=xdata, B::Int64=500, reg::Function=lp1, lb::Real=-Inf, ub::Real=Inf, kernel::Function=gaussiankernel, h::Real=-Inf)` compute 95% wild bootstrap confidence band at `xeval` for `lp0` and `lp1`



In the above functions,

 - `xeval` is the point(s) where the density or fitted value is calculated

 - `xdata` is the input X

 - `ydata` is the response vector y; should have same length with `xdata`

 - `reg` is the regression function `LP0` or `LP1`

 - `kernel` defaults to be `gaussiankernel`; should be a function

 - `h` is the bandwidth, should be a real scalar

 - `lb`, 'ub` are the boundary for x. Must be provide if you want to use Beta or Gamma kernel


##Demos

 - Kernel density estimate

        using Distributions
        x=rand(Normal(), 500)
        xeval=linspace(minimum(x), maximum(x), 100)
        den=kerneldensity(x, xeval=xeval)


 - Local regression

        y=2 .* x.^2 + rand(Normal(), 500)
        yfit0=lp0(x, y, xeval=xeval)
        yfit1=lp1(x, y, xeval=xeval)
        yfit0=npr(x, y, xeval=xeval, reg=lp0)
        yfit1=npr(x, y, xeval=xeval, reg=lp1)

 - Confidence Band

        cb=bootstrapCB(x, y, xeval=xeval)
        using Gadfly
        plot(layer(x=x, y=y, Geom.point), layer(x=xeval, y=yfit1, Geom.line, Theme(default_color=color("black"))), layer(x=xeval, y=cb[1,:], Geom.line, Theme(default_color=color("red"))), layer(x=xeval, y=cb[2,:], Geom.line, Theme(default_color=color("red"))))


##To be done

 - Implement multivariate version of kde and local regression

 - Add Gamma and Beta kernel which are useful to remove boundary biase


##Reference

 - [1] Lecture notes from Dr. Song Xi Chen

 - [2] Chen, Song Xi. "Beta kernel estimators for density functions." Computational Statistics & Data Analysis 31, no. 2 (1999): 131-145.

 - [3] Chen, Song Xi. "Probability density function estimation using gamma kernels." Annals of the Institute of Statistical Mathematics 52, no. 3 (2000): 471-480.

 - [4] W. Hardle and J. S. Marron (1991). Bootstrap Simultaneous Error Bars for Nonparametric Regression. _The Annals of Statistics_. Vol. 19, No. 2 (Jun., 1991), pp. 778-796

 - [5] W.Hardle and E. Mammen (1993). Comparing Nonparametric Versus Parametric Regression Fits. _The Annals of Statistics_. Vol. 21, No. 4 (Dec., 1993), pp. 1926-1947

 -  [6] Clifford M. Hurvich, Jeffrey S. Simonoff and Chih-Ling Tsai (1998). Smoothing Parameter Selection in Nonparametric Regression Using an Improved Akaike Information Criterion. _Journal of the Royal Statistical Society. Series B (Statistical Methodology)_, Vol. 60, No. 2 (1998), pp. 271-293





