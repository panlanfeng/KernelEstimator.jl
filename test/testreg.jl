#Multivariate kde and regression
x=rand(100,4)
y=x * ones(4) + x .^ 2 * ones(4)
xeval=rand(50,4)
den=KernelDensity(xeval, x)

yfit0=LP0(xeval,x,y)

#confidence band
cb=BootstrapCB(100, xeval, x, y)


#Goodness of fit
function SemiPredict(xdata::Matrix{Float64}, ydata::Vector{Float64})
   xb=xdata*Sind(xdata,ydata,GaussianKernel)
   LP0(xb,xb,ydata,GaussianKernel,1.0)
end
BootstrapGoodness(10, x, y, SemiPredict) 



##Univariate kde and regression
using Distributions
x=rand(Normal(), 100)
xeval=linspace(minimum(x), maximum(x), 50)
den=KernelDensity(xeval,x) 



y=2 .* x + rand(Normal(), 100)
yfit0=LP0(xeval, x, y)
yfit1=LP1(xeval, x, y)
