
##Univariate kde and regression
using Distributions
x=rand(Normal(10), 500)
xeval=linspace(minimum(x), maximum(x), 100)
h = bwlscv(x, gaussiankernel)
@test h>0
denvalues=kde(x, xeval=xeval)
@test all(denvalues .>= 0)

y=2 .* x.^2 .+ rand(Normal(0,5),500)
b0, b1 = linreg(x, y)
regfit = b0 .+ b1 .* x
yfit0=lp0(x, y, xeval=x)
@test sumabs2(regfit .- y) > sumabs2(yfit0.-y)

yfit1=lp1(x, y, xeval=x)
@test sumabs2(regfit .- y) > sumabs2(yfit1.-y)

yfit0=npr(x, y, xeval=x, reg=lp0)
@test sumabs2(regfit .- y) > sumabs2(yfit0.-y)

yfit1=npr(x, y, xeval=x, reg=lp1)
@test sumabs2(regfit .- y) > sumabs2(yfit1.-y)

#confidence band
yfit1=npr(x, y, xeval=xeval, reg=lp1)
cb=bootstrapCB(x, y, xeval=xeval)
@test mean(vec(cb[1,:]) .<= yfit1 .<= vec(cb[2,:])) > .9




###Bounded kernel density

x = rand(Beta(4,2), 500) * 10
xeval = linspace(0, 10, 100)
h = bwlscv(x./10, betakernel)
@test h>0
denvalues=kde(x, xeval=xeval, kernel=betakernel,h=h, lb=0.0,ub=10.0)
@test all(denvalues .> 0)


##Bounded local regression
y=2 .* x.^2 + x.*rand(Normal(0, 5), 500)
b0, b1 = linreg(x, y)
regfit = b0 .+ b1 .* x
yfit0=npr(x, y, xeval=x, reg=lp0, kernel=betakernel,lb=0.0, ub=10.0)
@test sumabs2(regfit .- y) > sumabs2(yfit0.-y)

yfit1=npr(x, y, xeval=x, reg=lp1, kernel=betakernel, lb=0.0,ub=10.0)
@test sumabs2(regfit .- y) > sumabs2(yfit1.-y)

yfit1=npr(x, y, xeval=xeval, reg=lp1, kernel=betakernel, lb=0.0,ub=10.0)
cb=bootstrapCB(x, y, xeval=xeval,reg=lp1, kernel=betakernel, lb=0.0,ub=10.0)
@test mean(vec(cb[1,:]) .<= yfit1 .<= vec(cb[2,:])) > .9
