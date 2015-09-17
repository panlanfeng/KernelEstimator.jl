
function kde(xdata::RealVector, xeval::RealVector, kernel::Function, h::Real)
    h > 0.0 || error("h < 0!")
    n = length(xdata)
    w = zeros(n)
    den = zeros(length(xeval))
    for i in 1:length(xeval)
        kernel(xeval[i], xdata, h, w, n)
        den[i] = mean(w)
    end
    return den
end
kde(xdata::RealVector, xeval::Real, kernel::Function, h::Real) = kde(xdaa, [xeval;], kernel, h)

function kde(xdata::RealVector; xeval::RealVector=xdata, lb::Real=-Inf, ub::Real=Inf, kernel::Function=gaussiankernel, h::Real=-Inf)

    xdata, xeval = boundit(xdata, xeval, kernel, lb, ub)
    if h <= 0
        h = bwlscv(xdata, kernel)
    end
    return kde(xdata, xeval, kernel, h)
end

kerneldensity(xdata::RealVector; xeval::RealVector=xdata, lb::Real=-Inf, ub::Real=Inf, kernel::Function=gaussiankernel, h::Real=-Inf) = kde(xdata, xeval=xeval, lb=lb, ub=ub, kernel=kernel, h=h)


# #univariate kernel density
# #
# function kde(xeval::Real, xdata::RealVector; kernel::Functor{3}=Gkernel(), h::Real=bwkd(xdata, kernel))
#     h > 0.0 || error("Bandwidth should be positive")
#     mean(kernel, xdata, xeval, h)
# end
# function kde(xeval::RealVector, xdata::RealVector; kernel::Function=gkernel, h::Real=bwkd(xdata, kernel))
#     h > 0.0 || error("Bandwidth should be positive")
#     Float64[mean(map(kernel, xdata, xeval[i], h)) for i=1:length(xeval)]
# end
# kde(xdata::RealVector;kernel::Function=gkernel, h::Real=bwcv(xdata, kernel)) = kde(xdata, xdata, h=h, kernel=kernel)


# function kerneldensity(xdata::RealVector; xeval::RealVector=xdata, lb::Real=-Inf, ub::Real=Inf, kernel::Function=gkernel, h::Real=-Inf)
#     if (lb == -Inf) & (ub == Inf)
#         if h==-Inf
#             return kde(xeval, xdata, kernel=kernel, h=bwkd(xdata, kernel))
#         else
#             return kde(xeval, xdata, kernel=kernel, h=h)
#         end
#     elseif (lb > -Inf) & (ub < Inf)
#         if h == -Inf
#             all(lb .<= xeval .<= ub) & all(lb .<= xdata .<= ub) || error("Your data are not in [lb,ub]")
#             kde((xeval .- lb)/(ub - lb), (xdata .- lb)./(ub - lb), kernel=betakernel, h=bwkd(xdata, betakernel))
#         else
#             all(lb .<= xeval .<= ub) & all(lb .<= xdata .<= ub) || error("Your data are not in [lb,ub]")
#             kde((xeval .- lb)/(ub - lb), (xdata .- lb)./(ub - lb), kernel=betakernel, h=h)
#         end
#     elseif (lb > -Inf) & (ub == Inf)
#         if h == -Inf
#             all(xeval .>= lb) & all(xdata .>= lb) || error("lb should be less than your data")
#             kde(xeval .- lb, xdata .- lb, kernel=gammakernel, h=bwkd(xdata, gammakernel))
#         else
#             all(xeval .>= lb) & all(xdata .>= lb) || error("lb should be less than your data")
#             kde(xeval .- lb, xdata .- lb, kernel=gammakernel, h=h)
#         end
#     else (lb == -Inf) & (ub < Inf)
#         if h == -Inf
#             all(xeval .<= ub) & all(xdata .<= ub) || error("ub should be greater than your data")
#             kde(ub .- xeval, ub .- xdata, kernel=gammakernel, h=bwkd(xdata, Gammakernel ))
#         else
#             all(xeval .<= ub) & all(xdata .<= ub) || error("ub should be greater than your data")
#             kde(ub .- xeval, ub .- xdata, kernel=gammakernel, h=h)
#         end
#     end
# end

# # #one side bound kde
# # function lkde(xeval::Real, xdata::RealVector, lb::Real; h::Real=bwcv(xdata, kernel))
# #     (xeval >= lb) & all(xdata .>= lb) || error("lb should be less than your data")
# #     kde(xeval - lb, xdata .- lb, kernel=GammaKernel(), h=h)
# # end
# # function lkde(xeval::RealVector, xdata::RealVector, lb::Real; h::Real=bwcv(xdata, kernel))
# #     all(xeval .>= lb) & all(xdata .>= lb) || error("lb should be less than your data")
# #     kde(xeval .- lb, xdata .- lb, kernel=GammaKernel(), h=h)
# # end
# # lkde(xdata::RealVector, lb::Real; h::Real=bwcv(xdata, kernel)) = lkde(xdata,xdata,lb,h=h)

# # function rkde(xeval::Real, xdata::RealVector, ub::Real; h::Real=bwcv(xdata, kernel))
# #     (xeval <= ub) & all(xdata .<= ub) || error("ub should be greater than your data")
# #     kde(ub - xeval, ub .- xdata, kernel=GammaKernel(), h=h)
# # end
# # function rkde(xeval::RealVector, xdata::RealVector, ub::Real; h::Real=bwcv(xdata, kernel))
# #     all(xeval .<= ub) & all(xdata .<= ub) || error("ub should be greater than your data")
# #     kde(ub .- xeval, ub .- xdata, kernel=GammaKernel(), h=h)
# # end
# # rkde(xdata::RealVector, ub::Real; h::Real=bwcv(xdata, kernel)) = rkde(xdata,xdata,ub,h=h)

# # #two side bound kde
# # function bkde(xeval::Real, xdata::RealVector, lb::Real, ub::Real; h::Real=bwcv(xdata, kernel))
# #     (lb <= xeval <= ub) & (lb .<= xdata .<= ub) || error("Your data are not in [lb,ub]")
# #     kde((xeval - lb)/(ub - lb), (xdata .- lb)./(ub - lb), kernel=BetaKernel(), h=h)
# # end
# # function bkde(xeval::RealVector, xdata::RealVector, lb::Real, ub::Real; h::Real=bwcv(xdata, kernel))
# #     all(lb .<= xeval .<= ub) & (lb .<= xdata .<= ub) || error("Your data are not in [lb,ub]")
# #     kde((xeval .- lb)./(ub - lb), (xdata .- lb)./(ub - lb), kernel=BetaKernel(), h=h)
# # end
# # bkde(xdata::RealVector, lb::Real, ub::Real; h::Real=bwcv(xdata, kernel))=bkde(xdata,xdata,lb,ub,h=h)




function kde(xdata::RealMatrix, xeval::RealMatrix, 
    kernel::Array{Function, 1}, h::RealVector)
    
    if any(h .<= 0)
        error("h < 0!")
    end
    m, p=size(xeval)
    n, p1 = size(xdata)
    if length(h) != p || p1 != p
        error("xeval and h should have same dimension as xdata")
    end
    den=zeros(m)
    for i=1:m
        wtmp = zeros(n)
        w = ones(n)
        for j=1:p
            kernel[j](xeval[i, j], xdata[:, j], h[j], wtmp, n)
            w .*= wtmp
        end
        den[i] = mean(w)
    end
    den
end
kde(xdata::RealMatrix, xeval::RealMatrix, kernel::Function, h::RealVector) = kde(xdata, xeval, [kernel for i in 1:size(xdata)[2]], h)
