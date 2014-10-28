# Nonparametric
The Julia package for nonparametric density estimate and regression. Currently includes univariate kernel density estimate, local constant regression (Nadaraya-watson estimator), local linear regression. Bootstrap confidence band is also provided for the regression estimation. The previous support for multivariate density estimation and regression are temporarily removed because of low efficiency.

[![Build Status](https://travis-ci.org/panlanfeng/Nonparametric.jl.png)](https://travis-ci.org/panlanfeng/Nonparametric.jl)

## Functions
This package provides the following functions:

 - `kde(xeval, xdata, kernel::Functor{3}=Gkernel(), h=bwcv(xdata,kernel))` do kernel density estimate

 - `LP0(xeval, xdata, ydata::Vector; kernel::Functor{3}=Gkernel(),h=bwlp0(xdata,ydata,kernel))` do local constant regression (or Nadaraya-Watson)

 - `LP1(xeval, xdata, ydata::Vector; kernel::Functor{3}=Gkernel(),h=bwlp1(xdata,ydata,kernel))` do local linear regression

 - `bwnormal(xdata::Vector)` select bandwidth for density estimate by rule of thumb

 - `bwcv(xdata, kernel)` select bandwidth for density estiamte by leave-one-out. Currently only support for Gaussian Kernel

 - `bwlp0(xdata, ydata::Vector, kernel)` select bandwidth for local constant regression using leave-one-out

 - `bwlp1(xdata, ydata::Vector, kernel)` select bandwidth for local linear regression using corrected AIC. See reference.

 - `BootstrapCB(B::Int64,xeval,xdata,ydata,reg,kernel,h)` compute 95% wild bootstrap confidence band at `xeval` for `LP0` and `LP1`



In the above functions,

 - `xeval` is the point(s) where the density or fitted value is calculated

 - `xdata` is the input X

 - `ydata` is the response vector y; should have same length with `xdata`

 - `reg` is the regression function `LP0` or `LP1`

 - `kernel` defaults to be `Gkernel()`; should be of type `Functor{3}` provided in package `NumericExtensions`

 - `h` is the bandwidth, should be a real scalar


##Demos

 - Kernel density estimate

        using Distributions
        x=rand(Normal(), 100)
        xeval=linspace(minimum(x), maximum(x), 50)
        den=kde(xeval,x)

 - Local constant regression

        y=2 .* x.^2 + rand(Normal(), 100)
        yfit0=LP0(xeval, x, y)
        yfit1=LP1(xeval, x, y)

 - Confidence Band

        cb=BootstrapCB(100, xeval, x, y, LP1, bwlp1(x, y))

        using Gadfly
        plot(layer(x=x, y=y, Geom.point), layer(x=xeval, y=yfit1, Geom.line, Theme(default_color=color("black"))), layer(x=xeval, y=cb[1,:], Geom.line, Theme(default_color=color("red"))), layer(x=xeval, y=cb[2,:], Geom.line, Theme(default_color=color("red"))))




##Reference

 - Lecture notes from Dr. Song Xi Chen

 - W. Hardle and J. S. Marron (1991). Bootstrap Simultaneous Error Bars for Nonparametric Regression. _The Annals of Statistics_. Vol. 19, No. 2 (Jun., 1991), pp. 778-796

 - W.Hardle and E. Mammen (1993). Comparing Nonparametric Versus Parametric Regression Fits. _The Annals of Statistics_. Vol. 21, No. 4 (Dec., 1993), pp. 1926-1947

 -  Clifford M. Hurvich, Jeffrey S. Simonoff and Chih-Ling Tsai (1998). Smoothing Parameter Selection in Nonparametric Regression Using an Improved Akaike Information Criterion. _Journal of the Royal Statistical Society. Series B (Statistical Methodology)_, Vol. 60, No. 2 (1998), pp. 271-293





