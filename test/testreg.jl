
##Univariate kerneldensity and regression
using Distributions
srand(2017);
x=rand(Normal(10), 500)
xeval=linspace(minimum(x), maximum(x), 100)
h = bwlscv(x, gaussiankernel)
@test h>0
denvalues=kerneldensity(x, xeval=xeval)
@test all(denvalues .>= 0)

y=2 .* x.^2 .+ rand(Normal(0,5),500)
b0, b1 = linreg(x, y)
regfit = b0 .+ b1 .* x
yfit0=localconstant(x, y, xeval=x)
@test sum(abs2, regfit .- y) > sum(abs2, yfit0.-y)

yfit1=locallinear(x, y, xeval=x)
@test sum(abs2, regfit .- y) > sum(abs2, yfit1.-y)

yfit0=npr(x, y, xeval=x, reg=localconstant)
@test sum(abs2, regfit .- y) > sum(abs2, yfit0.-y)

yfit1=npr(x, y, xeval=x, reg=locallinear)
@test sum(abs2, regfit .- y) > sum(abs2, yfit1.-y)

#confidence band
yfit1=npr(x, y, xeval=xeval, reg=locallinear)
cb=bootstrapCB(x, y, xeval=xeval)
@test mean(vec(cb[1,:]) .<= yfit1 .<= vec(cb[2,:])) > .8

#multivariate density estimation
x = rand(Normal(10), 500, 3)
denvalues = kerneldensity(x, h = [1.0, 1.0, 1.0])
@test all(denvalues .>= 0)

#multivariate regression
x = rand(Normal(10), 500, 3)
y = x * ones(3) .+ x.^2 *ones(3)
yfit1 = localconstant(x, y, xeval=x)
regfit = x*inv(x'*x)*x'*y
@test sum(abs2, regfit .- y) > sum(abs2, yfit1.-y)


###Bounded gamma kernel density and regression
x = rand(Gamma(4,2), 500)
xeval = linspace(0.01,20, 100)
h = bwlscv(x, gammakernel)
@test h>0
denvalues = kerneldensity(x, xeval=xeval, kernel=gammakernel, lb=0.0)
@test all(denvalues .> 0)

y=2 .* x.^2 + x.*rand(Normal(0, 5), 500)
b0, b1 = linreg(x, y)
regfit = b0 .+ b1 .* x
yfit0=npr(x, y, xeval=x, reg=localconstant, kernel=gammakernel,lb=0.0)
@test sum(abs2, regfit .- y) > sum(abs2, yfit0.-y)

yfit1=npr(x, y, xeval=x, reg=locallinear, kernel=gammakernel, lb=0.0)
@test sum(abs2, regfit .- y) > sum(abs2, yfit1.-y)

#bounded beta kernel density and regression
x = rand(Beta(4,2), 500) * 10
xeval = linspace(0, 10, 100)
h = bwlscv(x./10, betakernel)
@test h>0
denvalues=kerneldensity(x, xeval=xeval, kernel=betakernel,h=h, lb=0.0,ub=10.0)
@test all(denvalues .> 0)


##Bounded local regression
y=2 .* x.^2 + x.*rand(Normal(0, 5), 500)
b0, b1 = linreg(x, y)
regfit = b0 .+ b1 .* x
yfit0=npr(x, y, xeval=x, reg=localconstant, kernel=betakernel,lb=0.0, ub=10.0)
@test sum(abs2, regfit .- y) > sum(abs2, yfit0.-y)

yfit1=npr(x, y, xeval=x, reg=locallinear, kernel=betakernel, lb=0.0,ub=10.0)
@test sum(abs2, regfit .- y) > sum(abs2, yfit1.-y)

yfit1=npr(x, y, xeval=xeval, reg=locallinear, kernel=betakernel, lb=0.0,ub=10.0)
cb=bootstrapCB(x, y, xeval=xeval,reg=locallinear, kernel=betakernel, lb=0.0,ub=10.0)
@test mean(vec(cb[1,:]) .<= yfit1 .<= vec(cb[2,:])) > .8
