# #Multivariate kde and regression
# x=rand(100,4)
# y=x * ones(4) + x .^ 2 * ones(4)
# xeval=rand(50,4)
# den=kde(xeval, x)

# yfit0=LP0(xeval,x,y)

#confidence band
# cb=BootstrapCB(100, xeval, x, y)


#Goodness of fit
# function SemiPredict(xdata::Matrix{Float64}, ydata::Vector{Float64})
#    xb=xdata*Sind(xdata,ydata,GaussianKernel)
#    LP0(xb,xb,ydata,GaussianKernel,1.0)
# end
# BootstrapGoodness(10, x, y, SemiPredict)



##Univariate kde and regression
using Distributions
x=rand(Normal(), 500)
xeval=linspace(minimum(x), maximum(x), 100)
den=kerneldensity(x, xeval=xeval)


y=2 .* x.^2 + rand(Normal(), 500)
yfit0=lp0(x, y, xeval=xeval)
yfit1=lp1(x, y, xeval=xeval)
yfit0=npr(x, y, xeval=xeval, reg=lp0)
yfit1=npr(x, y, xeval=xeval, reg=lp1)
#confidence band
cb=bootstrapCB(x, y, xeval=xeval)


###Bounded kernel density

x = rand(Beta(4,2), 500) * 10
xeval = linspace(0, 10, 100)
den=kerneldensity(x, xeval=xeval, kernel=betakernel, lb=0.0,ub=10.0)
#plot(layer(x=xeval, y=pdf(Beta(4,2), xeval/10), Geom.line), layer(x=xeval, y=den, Geom.line, Theme(default_color=color("black"))))

##Bounded local regression
y=2 .* x.^2 + x.*rand(Normal(0, 5), 500)
yfit0=npr(x, y, xeval=xeval, reg=lp0, kernel=betakernel,lb=0.0, ub=10.0)
yfit1=npr(x, y, xeval=xeval, reg=lp1, kernel=betakernel, lb=0.0,ub=10.0)
cb=bootstrapCB(x, y, xeval=xeval,reg=lp1, kernel=betakernel, lb=0.0,ub=10.0)
#plot(layer(x=x, y=y, Geom.point), layer(x=xeval, y=yfit1, Geom.line, Theme(default_color=color("black"))))

# using Gadfly
# plot(layer(x=x, y=y, Geom.point), layer(x=xeval, y=yfit1, Geom.line, Theme(default_color=color("black"))), layer(x=xeval, y=cb[1,:], Geom.line, Theme(default_color=color("red"))), layer(x=xeval, y=cb[2,:], Geom.line, Theme(default_color=color("red"))))
