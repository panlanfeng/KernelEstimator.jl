# KernelEstimator 

[![KernelEstimator](http://pkg.julialang.org/badges/KernelEstimator_0.3.svg)](http://pkg.julialang.org/?pkg=KernelEstimator&ver=0.3)
[![KernelEstimator](http://pkg.julialang.org/badges/KernelEstimator_0.4.svg)](http://pkg.julialang.org/?pkg=KernelEstimator&ver=0.4)

[![Build Status](https://travis-ci.org/panlanfeng/KernelEstimator.jl.svg?branch=master)](https://travis-ci.org/panlanfeng/KernelEstimator.jl)
[![Coverage Status](https://coveralls.io/repos/panlanfeng/KernelEstimator.jl/badge.svg?branch=master)](https://coveralls.io/r/panlanfeng/KernelEstimator.jl?branch=master)


The Julia package for nonparametric kernel density estimate and regression. This package currently includes univariate kernel density estimate, local constant regression (Nadaraya-Watson regression) and local linear regression. It can also compute the Bootstrap confidence band [4]. 

This package provides Gamma and Beta kernel to deal with bounded density estimation and regression. These two kernels are free of boundary bias for one side and two sides bounded data respectively, see [2, 3]. In particular, least square cross validation (LSCV) bandwidth selection functions are implemented. 
Bandwidth selection is critical in kernel estimation. LSCV is recommended for kernel density estimation. Likelihood cross validation is provided but should be avoided because of known drawbacks. For regression problem, the bandwidth of local constant regression is selected using LSCV while that for local linear regression is chosen by AIC [6].

To install and use this package in Julia, 

	Pkg.add("KernelEstimator")
	using KernelEstimator

This package calculate densities via direct approach, i.e. adding kernel functions together. To add new kernel, you need to define a new function takes the same arguments as `gaussiankernel` and output the kernel weights at given point. If no bandwidth selection function is provide for that kernel, lscv with numeric integration will be used by default. 

## Functions
This package provides two major functions, `kde` for kernel density estimation and `npr` for nonparametric regression. For kernel density, you can simply use 

	xdata = randn(1000)
	kde(xdata)
	
or specify some options

	xeval = linspace(-3, 3, 100)
	bw = bwlscv(xdata, gaussiankernel)
	kde(xdata, xeval=xeval, lb=-Inf, ub=Inf, kernel=gaussiankernel,h = bw)

`xeval` specifies the position you want to evaluate the density at. Default to be the same as `xdata`. `lb` and `ub` means lower bound and upper bound of the data. If you specify either of them to be some finite value, user choice of kernel function will be suppressed and `gammakernel` will be used with a warning. If you specify both, `betakernel` is used with a warning if user's choice is different. 

For kernel regression, you can use

	x = rand(Beta(4,2), 500) * 10
	y=2 .* x.^2 + x .* rand(Normal(0, 5), 500)
	npr(x, y)
	
or change the default by

	npr(x, y, xeval=x, reg=lp1, kernel=betakernel,lb=0.0, ub=10.0)
	
`reg` specifies the order of local polynomial regression. You can choose `lp0`, local constant regression or `lp1`, local linear regression. `lp1` has better theoretical properties in prediction `y` and is used by default but is more computing intensive. 

There is also a function computing the bootstrap confidence interval for regression. 

	bootstrapCB(x, y; xeval=x, B=500, reg=lp1, lb=-Inf, ub=Inf, kernel=gaussiankernel)
	
`B` specifies the number of bootstrap sampling. 

 The following functions are also provided:

 - `lp0(xdata::RealVector, ydata::RealVector; xeval::RealVector=xdata, kernel::Function=gaussiankernel, h::Real=bwlp0(xdata,ydata,kernel))`, local constant regression (or Nadaraya-Watson)

 - `lp1(xdata::RealVector, ydata::RealVector; xeval::RealVector=xdata, kernel::Function=gaussiankernel, h::Real=bwlp0(xdata,ydata,kernel))`,  local linear regression

and bandwidth selection functions:

 - `bwnormal(xdata::Vector)`, bandwidth selection for density estimate by referencing to normal distribution

 - `bwlscv(xdata::RealVector, kernel::Function)`, bandwidth selection for density estimate by least square cross validation

 - `bwlcv(xdata::RealVector, kernel::Function)`, bandwidth selection for density estimate by likelihood cross validation

 - `bwlp0(xdata, ydata::Vector, kernel)`, bandwidth selection for local constant regression using LSCV

 - `bwlp1(xdata, ydata::Vector, kernel)`, bandwidth selection for local linear regression using corrected AIC. See reference [6]



The meaning of arguments:

 - `xeval` is the point(s) where the density or fitted value is calculated

 - `xdata` is the input X

 - `ydata` is the response vector y; should have same length with `xdata`

 - `reg` is the regression function, `lp0` or `lp1`

 - `kernel` defaults to be `gaussiankernel`; should be a function

 - `h` is the bandwidth, should be a real scalar; If negative, the default bandwidth selection method will be used to find the bandwidth and replace it

 - `lb`, `ub` are the boundary for x. Must provide if use Beta or Gamma kernel


## Demos

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
        plot(layer(x=x, y=y, Geom.point), layer(x=xeval, y=yfit1, Geom.line, Theme(default_color=color("black"))),
          layer(x = xeval, y = cb[1,:], Geom.line, Theme(default_color=color("red"))),
          layer(x=xeval, y=cb[2,:], Geom.line, Theme(default_color=color("red"))))




## Reference

 - [1] Lecture notes from Dr. Song Xi Chen

 - [2] Chen, Song Xi. "Beta kernel estimators for density functions." _Computational Statistics & Data Analysis_ 31, no. 2 (1999): 131-145.

 - [3] Chen, Song Xi. "Probability density function estimation using gamma kernels." _Annals of the Institute of Statistical Mathematics_ 52, no. 3 (2000): 471-480.

 - [4] W. Hardle and J. S. Marron (1991). Bootstrap Simultaneous Error Bars for Nonparametric Regression. _The Annals of Statistics_. Vol. 19, No. 2 (Jun., 1991), pp. 778-796

 - [5] W.Hardle and E. Mammen (1993). Comparing Nonparametric Versus Parametric Regression Fits. _The Annals of Statistics_. Vol. 21, No. 4 (Dec., 1993), pp. 1926-1947

 -  [6] Clifford M. Hurvich, Jeffrey S. Simonoff and Chih-Ling Tsai (1998). Smoothing Parameter Selection in Nonparametric Regression Using an Improved Akaike Information Criterion. _Journal of the Royal Statistical Society. Series B (Statistical Methodology)_, Vol. 60, No. 2 (1998), pp. 271-293
