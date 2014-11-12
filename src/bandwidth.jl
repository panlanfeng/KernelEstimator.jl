####bandwidth selector for kernel density

#rule of thumb
function bwnormal(xdata::RealVector)
    0.9 * min((quantile(xdata, .75) - quantile(xdata, .25)) / 1.34, std(xdata)) * length(xdata) ^ (-0.2)
end
function midrange(x::RealVector)
    lq, uq = quantile(x, [.25, .75])
    uq - lq
end

# J(h)=∑ᵢⱼK'((xᵢ - xⱼ)/h) /(n²h) + 2K(0)/(nh) =
# where K'(u) = invsqrt2pi(exp(-0.25u*u)/sqrt(2) - 2exp(-0.5u*u))
#J(h) = invsqrt2pi/(n²h) ∑ᵢⱼ (exp(-0.25u*u)/sqrt(2) - 2exp(-0.5u*u)) + 2 * invsqrt2pi /nh
#J(h) = 2*invsqrt2pi/(n²h) ∑{i<j} (exp(-0.25u*u)/sqrt(2) - 2exp(-0.5u*u)) + invsqrt2pi/sqrt(2)nh
#For normal kernel
function Jh{T<:FloatingPoint}(xdata::RealVector{T}, h::T, n::Int)
    tmp = 0.0
    @inbounds for i in 1:(n-1)
        for j in (i+1):n
            u = (xdata[i] - xdata[j])/h
            u = exp(-0.25*u*u)
            tmp += u/sqrt(2) - 2*u*u
        end
    end
    2*tmp / (n*n*h) + 1/(sqrt(2)*n*h)
end
Jh(xdata::RealVector, h::Real, n::Int)=Jh(float(xdata), float(h), n)


#leave one out
function leaveoneout(xdata::RealVector, kernel::Function, h::Real, w::Vector, n::Int)
    ind = 1
    ind_end = 1+n
    ll = 0.0
    @inbounds while ind < ind_end
        kernel(xdata[ind], xdata, h, w, n)
        w[ind] = 0.0
        ll += mean(w)
        ind += 1
    end
    ll * 2 / (n-1)
end

# function Jh_betakernel{T<:FloatingPoint}(xdata::Vector{T}, h::T, w::Vector, n::Int, xlb::T, xub::T)
#     pquadrature(x->begin betakernel(x, xdata,h,w,n); mean(w)^2; end, xlb, xub, maxevals=100)[1] - leaveoneout(xdata, betakernel, h, w, n)
# end
# function Jh_gammakernel{T<:FloatingPoint}(xdata::Vector{T}, h::T, w::Vector, n::Int, xlb::T, xub::T)
#     pquadrature(x->begin gammakernel(x, xdata,h,w,n); mean(w)^2; end, xlb, xub, maxevals=100)[1] - leaveoneout(xdata, gammakernel, h, w, n)
# end
function Jh{T<:FloatingPoint}(xdata::RealVector{T}, kernel::Function, h::T, w::Vector, n::Int, xlb::T, xub::T)
    pquadrature(x->begin kernel(x, xdata,h,w,n); mean(w)^2; end, xlb, xub, maxevals=100)[1] - leaveoneout(xdata, kernel, h, w, n)
end

#Least Squares cross validation for Gaussian Kernel
#May fail to work if there are multiple equial x_i
# Silverman suggest search interval be (0.25, 1.5)n^(-0.2)σ
function bwlscv(xdata::RealVector, kernel::Function)
    n=length(xdata)
    h0=bwnormal(xdata)
    if kernel == gaussiankernel
        return Optim.optimize(h -> Jh(xdata, h, n), 0.01*h0, 10*h0, iterations=100, abs_tol=h0/n).minimum
    end

    xlb, xub = extrema(xdata)
    h0=midrange(xdata)
    hlb = h0/n
    hub = h0
    w  = zeros(n)
    if kernel == betakernel
        xlb = 0.0
        xub = 1.0
        hub = 0.25
    elseif kernel == gammakernel
        xlb = 0.0
    end
    return Optim.optimize(h -> Jh(xdata, kernel, h, w, n, xlb,xub), hlb, hub, iterations=100,abs_tol=h0/n^2).minimum
end

# likelihood cross validation for beta and gamma kernel
#there seems no other easy way; least square cross validation can be formidable because their convolution have no close form
#may also work for other kernels, but likelihood cv has some known disadvantages.
#Not recommended, but easy to implement for arbitrary kernel
function lcv(xdata::RealVector, kernel::Function, h::Real, w::Vector, n::Int)
#     -mean(kde(xdata,xdata,kernel,h)) + mean(map(kernel, xdata, xdata, h))
    ind = 1
    ind_end = 1+n
    ll = 1.0
    @inbounds while ind < ind_end
        kernel(xdata[ind], xdata, h, w, n)
        w[ind] = 0.0
        ll += log(mean(w))
        ind += 1
    end
    -ll # *n /(n-1)
end
function bwlcv(xdata::RealVector, kernel::Function)
    n = length(xdata)
    w = zeros(n)
    h0=midrange(xdata)
    hlb = h0/n^2
    hub = h0
    if kernel==betakernel
        hub = 0.25
    end
    return Optim.optimize(h->lcv(xdata,kernel,h,w,n), hlb, hub, iterations=100,abs_tol=h0/n^2).minimum
end


# #to be implemented
# function bwkd(xdata::RealVector, kernel::Function)
#     n=length(xdata)
#     h0=bwnormal(xdata)
#     if kernel==gaussiankernel
#         return Optim.optimize(h -> Jh(xdata, h), h0/n, n*h0).minimum
#     elseif kernel==betakernel
#         return Optim.optimize(h->lcv(xdata, h, kernel), h0^2/n^2, 0.25).minimum
#     elseif kernel==gammakernel
#         return Optim.optimize(h->lcv(xdata, h, kernel), h0^2/n^2, n*h0).minimum
#     else
#         warn("No bandwidth selector for this kernel, likelihood cross validation is used")
#         return Optim.optimize(h->lcv(xdata, h, kernel), h0/n^2, n*h0).minimum
#     end
# end


#multivariate
# function BandwidthLSCV(xdata::Matrix{Float64}, kernel::KernelType=Gaussian)
#     (n, p)=size(xdata)
#     h0=BandwidthNormalReference(reshape(xdata[:,1], n))
#     function res(h::Vector{Float64})
#       tmp1=0.0
#       tmp2=0.0
#       for i in 1:n
#         for j in 1:n
#           if j == i
#             continue
#           end
#           #tmp += kernel(xdata[i], xdata[j], sqrt(2)*h) - 2*kernel(xdata[i], xdata[j], h)
#           #xdiff = ((xdata[i,:] .- xdata[j,:]) ./ h)^2
#           #tmp += (2^(-1/2)*exp(-xdiff/4) - 2*exp(-xdiff/2))
#           xdiff=[(xdata[i,k] - xdata[j,k])::Float64 for k in 1:p]
#           tmp1 += kernel.Convolution(xdiff, h)
#           tmp2 += kernel.Density(xdiff,zeros(p),h)
#         end
#       end
#       tmp1 / (n ^ 2) - tmp2 / (n * (n - 1)) * 2 + kernel.Convolution(zeros(p), h)/n

#     end
#    return optimize(res, [h0 for i in 1:p], iterations=100).minimum .+ .1/n
# end

#leave-one-out LSCV. 1/n \sum((yi - yihat)/(1-wi))
function lscvlp0(xdata::RealVector, ydata::RealVector, kernel::Function, h::Real, w::Vector, n::Int)
    tmp = 0.0
    ind = 1
    ind_end = 1 + n
    @inbounds while ind < ind_end
#         map!(kernel, w, xdata, xdata[i], h)
        kernel(xdata[ind], xdata, h, w, n)
        divide!(w, sum(w))
        tmp += abs2((wsum(w, ydata) - ydata[ind])/(1-w[ind]))
        ind += 1
    end
    tmp/n
end

function bwlp0(xdata::RealVector, ydata::RealVector, kernel::Function=gaussiankernel)
    n=length(xdata)
    length(ydata)==n || error("length(ydata) != length(xdata)")
    w = ones(n)
    if kernel == gaussiankernel
        h0= bwnormal(xdata)
        hlb = 0.1*h0
        hub = 10*h0
    elseif kernel == betakernel
        h0 = midrange(xdata)
        hlb = h0/n
        hub = 0.25
    elseif kernel == gammakernel
        h0 = midrange(xdata)
        hlb = h0/n
        hub = h0
    end
    Optim.optimize(h->lscvlp0(xdata, ydata, kernel, h, w, n), hlb, hub).minimum
end

#see reference:Smoothing Parameter Selection in Nonparametric Regression Using an Improved Akaike Information Criterion
# Clifford M. Hurvich, Jeffrey S. Simonoff and Chih-Ling Tsai
# Journal of the Royal Statistical Society. Series B (Statistical Methodology), Vol. 60, No. 2 (1998), pp. 271-293
#http://www.jstor.org/stable/2985940
function AIClp1(xdata::RealVector, ydata::RealVector, kernel::Function, h::Real, w::Vector, n::Int)
    tmp = 0.0
    traceH = 0.0
    ind = 1
    ind_end = 1+n
    @inbounds while ind < ind_end
#           map!(kernel, w, xdata, xdata[i], h)
          kernel(xdata[ind], xdata, h, w, n)
          s0 = sum(w)
          s1 = s0*xdata[ind] - wsum(w, xdata)
          s2 = wsumsqdiff(w, xdata, xdata[ind])
          sy0 = wsum(w, ydata)
          sy1 = NumericExtensions.wsum(w, YXdiff(), xdata, xdata[ind], ydata)
          tmp+=abs2((s2 * sy0 - s1 * sy1) /(s2 * s0 - s1 * s1) - ydata[ind])
          traceH += s0*w[ind]/(s2 * s0 - s1 * s1)
          ind += 1
    end
    tmp/n  + 2*(traceH+1)/(n-traceH-2)
end

function bwlp1(xdata::RealVector, ydata::RealVector, kernel::Function=gaussiankernel)
    n=length(xdata)
    length(ydata)==n || error("length(ydata) != length(xdata)")
    w = ones(n)
    if kernel == gaussiankernel
        h0= bwnormal(xdata)
        hlb = 0.1*h0
        hub = 10*h0
    elseif kernel == betakernel
        h0 = midrange(xdata)
        hlb = h0/n
        hub = 0.25
    elseif kernel == gammakernel
        h0 = midrange(xdata)
        hlb = h0/n
        hub = h0
    end
    Optim.optimize(h->AIClp1(xdata, ydata, kernel, w, n), hlb, hub).minimum
end

function bwreg(xdata::RealVector, ydata::RealVector, reg::Function, kernel::Function=gaussiankernel)

    if reg == LP1
        return bwlp1(xdata, ydata, kernel)
    else
        return bwlp0(xdata, ydata, kernel)
    end
end

# #leave-one-out LSCV for multivariate NW
# function BandwidthLSCVReg(xdata::Matrix{Float64}, ydata::Vector{Float64}, reg::Function=LP0, kernel::Function=GaussianKernel)

#   (n,p)=size(xdata)
#   h0=BandwidthNormalReference(reshape(xdata[:,1], n))

#   function res(h::Vector{Float64})
#       ls=0.0
#       for i in 1:n
#           # xdata[i,:] is actually a vector, so the returning value is still a vector
#           ls += abs2(ydata[i] - reg(xdata[i,:],xdata[[1:(i-1),(i+1):end],:], ydata[[1:(i-1),(i+1):end]],kernel,h)[1])

#       end
#       ls/n
#   end
#   optimize(res, [h0 for i in 1:p], iterations=100).minimum .+ .1 / n

# end

