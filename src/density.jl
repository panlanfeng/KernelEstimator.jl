
function kde{T<:Real}(x::Vector{T}, xdata::Vector{T}, kernel::Function, h::T)
    h > 0.0 || error("h < 0!")
    n = length(xdata)
    w = zeros(n)
    den = zeros(length(x))
    for i in 1:length(x)
        kernel(x[i], xdata, h, w, n)
        den[i] = mean(w)
    end
    return den
end
function kde{T<:Real}(x::T, xdata::Vector{T}, kernel::Function, h::T)
    h > 0.0 || error("h < 0!")
    n = length(xdata)
    w = zeros(n)
    kernel(x, xdata, h, w, n)
    return mean(w)
end

function kerneldensity{T<:Real}(xdata::Vector{T}; xeval::Vector{T}=xdata, lb::T=-Inf, ub::T=Inf, kernel::Function=gaussiankernel, h::T=-Inf)

    if (lb > -Inf) & (ub < Inf)
        all(lb .<= xeval .<= ub) & all(lb .<= xdata .<= ub) || error("Your data are not in [lb,ub]")
        xeval = (xeval .- lb)./(ub - lb)
        xdata = (xdata .- lb)./(ub - lb)
        if kernel != betakernel
            kernel = betakernel
            warn("kernel is set to be beta kernel")
        end
    elseif (lb > -Inf) & (ub == Inf)
        all(xeval .>= lb) & all(xdata .>= lb) || error("lb should be less than your data")
        xeval = xeval .- lb
        xdata = xdata .- lb
        if kernel != gammakernel
            kernel = gammakernel
            warn("kernel is set to be beta kernel")
        end
    elseif (lb == -Inf) & (ub < Inf)
        all(xeval .<= ub) & all(xdata .<= ub) || error("ub should be greater than your data")
        xeval = ub .- xeval
        xdata = ub .- xdata
        if kernel != gammakernel
            kernel = gammakernel
            warn("kernel is set to be beta kernel")
        end
    end
    if h <= 0
        h = bwlscv(xdata, kernel)
    end
    return kde(xeval, xdata, kernel, h)
end

function kerneldensity(xdata::RealVector; xeval::RealVector=xdata, lb::Real=-Inf, ub::Real=Inf, kernel::Function=gaussiankernel, h::Real=-Inf)
    xdata, xeval, lb, ub, h = promote(xdata, xeval, lb, ub, h)
    kerneldensity(xdata, xeval = xeval, lb=lb,ub=ub,kernel=kernel,h=h)
end

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


# # #MultiVariate kernel density
# # function KernelDensity(xeval::Vector{Float64}, xdata::Matrix{Float64},
# #   kernel::KernelType=Gaussian, h::Vector{Float64}=BandwidthLSCV(xdata,kernel))


# #   if any(h .<= 0)
# #     error("Bandwidth should be positive")
# #   end
# #   (n, p)=size(xdata)
# #   if length(h) == 1 && p==1
# #     return KernelDensity(xeval, xdata, kernel, h)
# #   end

# #   if length(xeval) != p || length(h) != p
# #     error("xeval should have same dimension as xdata")
# #   end


# #   s0=0.0
# #   xi=zeros(p)
# #   for i = 1:n
# #     for j=1:p
# #       xi[j]=xdata[i,j]
# #     end
# #     s0 = s0 + kernel.Density(xeval, xi, h)
# #   end
# #   s0 / length(xdata)
# # end

# # function KernelDensity(xeval::Matrix{Float64}, xdata::Matrix{Float64},
# #   kernel::KernelType=Gaussian, h::Vector{Float64}=BandwidthLSCV(xdata,kernel))

# #   if any(h .<= 0)
# #     error("xeval should have same dimension as xdata")
# #   end
# #   (m, p)=size(xeval)
# #   if length(h) != p || size(xdata)[2] != p
# #     error("xeval should have same dimension as xdata")
# #   end
# #   xi_eval=zeros(p)
# #   den=zeros(m)
# #   for i=1:m
# #     for j=1:p
# #       xi_eval[j]=xeval[i, j]
# #     end

# #     den[i]=KernelDensity(xi_eval, xdata, kernel, h)
# #   end
# #   den
# # end
