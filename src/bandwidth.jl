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
# where K'(u) = invsqrt2π(exp(-0.25u*u)/sqrt(2) - 2exp(-0.5u*u))
#J(h) = invsqrt2π/(n²h) ∑ᵢⱼ (exp(-0.25u*u)/sqrt(2) - 2exp(-0.5u*u)) + 2 * invsqrt2π/nh
#J(h) = 2*invsqrt2π/(n²h) ∑{i<j} (exp(-0.25u*u)/sqrt(2) - 2exp(-0.5u*u)) + invsqrt2π/sqrt(2)nh
#For normal kernel
# function Jh(xdata::RealVector, h::Real, n::Int)
#     tmp = 0.0
#     @inbounds for i in 1:(n-1)
#         for j in (i+1):n
#             u = (xdata[i] - xdata[j])/h
#             # u = exp(-0.25*u*u)
#             u = -0.25*u*u
#             u = exp(u)
#             tmp += u*invsqrt2 - 2*u*u
#         end
#     end
#     2*tmp / (n*n*h) + 1/(sqrt2*n*h)
# end
#For gaussiankernel, equivalent to the above one. The above version requires
#less computing but is inefficient in practice due to the exp function.
function Jh(xdata::RealVector, h::Real, w::Vector, n::Int)
    ll = 0.0
    @inbounds for ind in 1:n
        gaussiankernel(xdata[ind], xdata, sqrt2*h, w, n)
        ll += mean(w)
        # gaussiankernel(xdata[ind], xdata, h, w, n)
        # w[ind] = 0.0
        # ll -= 2*mean(w)
    end
    ll / n - leaveoneout(xdata, gaussiankernel, h, w, n)
end


#for general kernel
function Jh(xdata::RealVector, kernel::Function, h::Real, w::Vector, n::Int, xlb::Real, xub::Real)
    pquadrature(x->begin kernel(x, xdata,h,w,n); mean(w)^2; end, xlb, xub, maxevals=200)[1] - leaveoneout(xdata, kernel, h, w, n)
end
function leaveoneout(xdata::RealVector, kernel::Function, h::Real, w::Vector, n::Int)

    ll = 0.0
    @inbounds for ind in 1:n
        kernel(xdata[ind], xdata, h, w, n)
        w[ind] = 0.0
        ll += mean(w)
    end
    ll * 2 / (n-1)
end
#For betakernel
function Jh(xdata::RealVector, logxdata::RealVector,log1_xdata::RealVector, kernel::Function, h::Real, w::Vector, n::Int, xlb::Real, xub::Real)
    pquadrature(x->begin kernel(x, logxdata, log1_xdata, h,w,n); mean(w)^2; end, xlb, xub, maxevals=200)[1] - leaveoneout(xdata, logxdata, log1_xdata, kernel, h, w, n)
end
function leaveoneout(xdata::RealVector, logxdata::RealVector, log1_xdata::RealVector, kernel::Function, h::Real, w::Vector, n::Int)

    ll = 0.0
    @inbounds for ind in 1:n
        kernel(xdata[ind], logxdata, log1_xdata, h, w, n)
        w[ind] = 0.0
        ll += mean(w)
    end
    ll * 2 / (n-1)
end
#For gammakernel
function Jh(xdata::RealVector, logxdata::RealVector, kernel::Function, h::Real, w::Vector, n::Int, xlb::Real, xub::Real)
    pquadrature(x->begin kernel(x, xdata, logxdata, h,w,n); mean(w)^2; end, xlb, xub, maxevals=200)[1] - leaveoneout(xdata, logxdata, kernel, h, w, n)
end
function leaveoneout(xdata::RealVector, logxdata::RealVector, kernel::Function, h::Real, w::Vector, n::Int)

    ll = 0.0
    @inbounds for ind in 1:n
        kernel(xdata[ind], xdata, logxdata, h, w, n)
        w[ind] = 0.0
        ll += mean(w)
    end
    ll * 2 / (n-1)
end
#Least Squares cross validation for Gaussian Kernel
#May fail to work if there are multiple equial x_i
# Silverman suggest search interval be (0.25, 1.5)n^(-0.2)σ
function bwlscv(xdata::RealVector, kernel::Function)
    n=length(xdata)
    w  = zeros(n)
    h0=bwnormal(xdata)
    #when n is large, leaveoneout is more expensive than numeric integration
    if kernel == gaussiankernel && n<200
        return Optim.minimizer(Optim.optimize(h -> Jh(xdata, h, w, n), 0.01*h0, 10*h0, iterations=200, abs_tol=h0/n))
    end

    xlb, xub = extrema(xdata)
    h0=midrange(xdata)
    hlb = h0/n
    hub = h0
    if kernel == betakernel
        xlb = 0.0
        xub = 1.0
        hub = 0.25
        logxdata = Yeppp.log(xdata)
        log1_xdata = Yeppp.log(1.0 .- xdata)
        return Optim.minimizer(Optim.optimize(h -> Jh(xdata, logxdata, log1_xdata, kernel, h, w, n, xlb,xub), hlb, hub, iterations=200,abs_tol=h0/n^2))
    elseif kernel == gammakernel
        xlb = 0.0
        logxdata = Yeppp.log(xdata)
        return Optim.minimizer(Optim.optimize(h -> Jh(xdata, logxdata, kernel, h, w, n, xlb,xub), hlb, hub, iterations=200,abs_tol=h0/n^2))
    end
    return Optim.minimizer(Optim.optimize(h -> Jh(xdata, kernel, h, w, n, xlb,xub), hlb, hub, iterations=200,abs_tol=h0/n^2))
end

# likelihood cross validation for beta and gamma kernel
#there seems no other easy way; least square cross validation can be formidable because their convolution have no close form
#may also work for other kernels, but likelihood cv has some known disadvantages.
#Not recommended, but easy to implement for arbitrary kernel
function lcv(xdata::RealVector, kernel::Function, h::Real, w::Vector, n::Int)
#     -mean(kerneldensity(xdata,xdata,kernel,h)) + mean(map(kernel, xdata, xdata, h))
    ind = 1
    ind_end = 1+n
    ll = 0.0
    @inbounds while ind < ind_end
        kernel(xdata[ind], xdata, h, w, n)
        w[ind] = 0.0
        ll += log(mean(w))
        ind += 1
    end
    -ll
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
    return Optim.minimizer(Optim.optimize(h->lcv(xdata,kernel,h,w,n), hlb, hub, iterations=200,abs_tol=h0/n^2))
end

function lcv(xdata::RealMatrix, kernel::Array{Function, 1}, h::RealVector, w::Vector, n::Int)
#     -mean(kerneldensity(xdata,xdata,kernel,h)) + mean(map(kernel, xdata, xdata, h))
    if any(h .<= 0.0)
        return Inf
    end
    ind = 1
    ind_end = 1+n
    ll = 0.0
    p = size(xdata)[2]
    @inbounds while ind < ind_end
        fill!(w, 1.0)
        wtmp = ones(n)
        for j in 1:p
            kernel[j](xdata[ind, j], xdata[:, j], h[j], wtmp, n)
            for k in 1:n
                w[k] = w[k] * wtmp[k]
            end
        end
        divide!(w, sum(w))

        w[ind] = 0.0
        ll += log(mean(w))
        ind += 1
    end
    -ll
end
function bwlcv(xdata::RealMatrix, kernel::Array{Function, 1})
    n, p = size(xdata)
    w = ones(n)
    h0 = zeros(p)
    hlb = zeros(p)
    hub = zeros(p)
    for j in 1:p
        if kernel[j] == gaussiankernel
            h0[j] = bwnormal(xdata[:, j])
            hlb[j] = 0.1*h0[j]
            hub[j] = 10*h0[j]
        elseif kernel[j] == betakernel
            h0[j] = midrange(xdata[:, j])
            hlb[j] = h0[j]/n
            hub[j] = 0.25
        elseif kernel[j] == gammakernel
            h0[j] = midrange(xdata[:, j])
            hlb[j] = h0[j]/n
            hub[j] = h0[j]
        end
    end
    Optim.minimizer(Optim.optimize(h->lcv(xdata, kernel, h, w, n), h0))
end


#leave-one-out LSCV. 1/n \sum((yi - yihat)/(1-wi))
function lscvlocalconstant(xdata::RealVector, ydata::RealVector, kernel::Function, h::Real, w::Vector, n::Int)
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

function bwlocalconstant(xdata::RealVector, ydata::RealVector, kernel::Function=gaussiankernel)
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
    Optim.minimizer(Optim.optimize(h->lscvlocalconstant(xdata, ydata, kernel, h, w, n), hlb, hub))
end

#see reference:Smoothing Parameter Selection in Nonparametric Regression Using an Improved Akaike Information Criterion
# Clifford M. Hurvich, Jeffrey S. Simonoff and Chih-Ling Tsai
# Journal of the Royal Statistical Society. Series B (Statistical Methodology), Vol. 60, No. 2 (1998), pp. 271-293
#http://www.jstor.org/stable/2985940
function AIClocallinear(xdata::RealVector, ydata::RealVector, kernel::Function, h::Real, w::Vector, n::Int)
    tmp = 0.0
    traceH = 0.0
    ind = 1
    ind_end = 1+n
    @inbounds while ind < ind_end
#           map!(kernel, w, xdata, xdata[i], h)
          kernel(xdata[ind], xdata, h, w, n)
          s0 = sum(w)
          s1 = s0*xdata[ind] - wsum(w, xdata)
          s2 = wsumsqdiff(w, xdata, xdata[ind], n)
          sy0 = wsum(w, ydata)
          sy1 = wsumyxdiff(w, xdata, xdata[ind], ydata, n)
          tmp+=abs2((s2 * sy0 - s1 * sy1) /(s2 * s0 - s1 * s1) - ydata[ind])
          traceH += s0*w[ind]/(s2 * s0 - s1 * s1)
          ind += 1
    end
    tmp/n  + 2*(traceH+1)/(n-traceH-2)
end

function bwlocallinear(xdata::RealVector, ydata::RealVector, kernel::Function=gaussiankernel)
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
    Optim.minimizer(Optim.optimize(h->AIClocallinear(xdata, ydata, kernel,h, w, n), hlb, hub))
end

function bwreg(xdata::RealVector, ydata::RealVector, reg::Function, kernel::Function=gaussiankernel)

    if (reg == locallinear) || (reg == lp1)
        return bwlocallinear(xdata, ydata, kernel)
    elseif (reg == localconstant) || (reg == lp0)
        return bwlocalconstant(xdata, ydata, kernel)
    else
        error("I don't know what is $reg")
    end
end



# #leave-one-out LSCV for multivariate NW
function lscvlocalconstant(xdata::RealMatrix, ydata::RealVector, kernel::Array{Function, 1}, h::RealVector, w::Vector, n::Int)
    if any(h .<= 0.0)
        return Inf
    end
    tmp = 0.0
    ind = 1
    ind_end = 1+n
    p=size(xdata)[2]
    @inbounds while ind < ind_end
        fill!(w, 1.0)
        wtmp = ones(n)
        for j in 1:p
            kernel[j](xdata[ind, j], xdata[:, j], h[j], wtmp, n)
            for k in 1:n
                w[k] = w[k] * wtmp[k]
            end
        end
        divide!(w, sum(w))
        tmp += abs2((wsum(w, ydata) - ydata[ind])/(1-w[ind]))
        ind += 1
    end
    tmp/n
end

function bwlocalconstant(xdata::RealMatrix, ydata::RealVector, kernel::Array{Function, 1} = Function[gaussiankernel for i in 1:size(xdata)[2]])
    n, p = size(xdata)
    w = ones(n)
    h0 = zeros(p)
    hlb = zeros(p)
    hub = zeros(p)
    for j in 1:p
        if kernel[j] == gaussiankernel
            h0[j] = bwnormal(xdata[:, j])
            hlb[j] = 0.1*h0[j]
            hub[j] = 10*h0[j]
        elseif kernel[j] == betakernel
            h0[j] = midrange(xdata[:, j])
            hlb[j] = h0[j]/n
            hub[j] = 0.25
        elseif kernel[j] == gammakernel
            h0[j] = midrange(xdata[:, j])
            hlb[j] = h0[j]/n
            hub[j] = h0[j]
        end
    end
    h_output=Optim.minimizer(Optim.optimize(h->lscvlocalconstant(xdata, ydata, kernel, h, w, n), h0))
    if any(h_output .<= 0.0)
        for j in 1:p
            if h_output[j] .<= 0.0
                h_output[j] = 2.* h0[j]
            end
        end
        h_output = Optim.minimizer(Optim.optimize(h->lscvlocalconstant(xdata, ydata, kernel, h, w, n), h_output))
    end
    h_output
end
