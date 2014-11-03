####bandwidth selector for kernel density

#rule of thumb
function bwnormal(xdata::RealVector)
    0.9 * min((quantile(xdata, .75) - quantile(xdata, .25)) / 1.34, std(xdata)) * length(xdata) ^ (-0.2)
end


# J(h)=∑ᵢⱼK'((xᵢ - xⱼ)/h) /(n²h) + 2K(0)/(nh) =
# where K'(u) = invsqrt2pi(exp(-0.25u*u)/sqrt(2) - 2exp(-0.5u*u))
#J(h) = invsqrt2pi/(n²h) ∑ᵢⱼ (exp(-0.25u*u)/sqrt(2) - 2exp(-0.5u*u)) + 2 * invsqrt2pi /nh
#J(h) = 2*invsqrt2pi/(n²h) ∑{i<j} (exp(-0.25u*u)/sqrt(2) - 2exp(-0.5u*u)) + invsqrt2pi/sqrt(2)nh
#For normal kernel
function Jh{T<:FloatingPoint}(xdata::Vector{T}, h::T)
    n=length(xdata)
    tmp = 0.0
    for i in 1:(n-1)
        for j in (i+1):n
            u = (xdata[i] - xdata[j])/h
            u = exp(-0.25*u*u)
            tmp += u/sqrt(2) - 2*u*u
        end
    end
    2*tmp / (n*n*h) + 1/(sqrt(2)*n*h)
end
Jh(xdata::RealVector, h::Real)=Jh(float(xdata), float(h))


#Leave-one-out cross validation for Gaussian Kernel
#May fail to work if there are multiple equial x_i
# Silverman suggest search interval be (0.25, 1.5)n^(-0.2)σ
function bwcv(xdata::RealVector)
    n=length(xdata)
    h0=bwnormal(xdata)
    return Optim.optimize(h -> Jh(xdata, h), 0.25*h0, 10*h0).minimum  # add a lower order item to avoid 0 bandwidth
end

# likelihood cross validation for beta and gamma kernel
#there seems no other easy way; least square cross validation can be formidable because their convolution have no close form
#may also work for other kernels, but likelihood cv has some known disadvantages.
function lcv(xdata::RealVector, h::Real, kernel::Functor{3})
    n = length(xdata)
    -mean(kde(xdata,xdata,kernel=kernel,h=h)) + mean(kernel, xdata, xdata, h)
end

#to be implemented
function bwkd(xdata::RealVector, kernel::Functor{3})
    n=length(xdata)
    h0=bwnormal(xdata)
    if kernel==Gkernel()
        return Optim.optimize(h -> Jh(xdata, h), h0/n, n*h0).minimum
    elseif (kernel==Betakernel()) | (kernel==Gammakernel())
        return Optim.optimize(h->lcv(xdata, h, kernel), h0^2/n^2, n*h0^2).minimum
    else
        warning("No bandwidth selector for this kernel, likelihood cross validation is used")
        return Optim.optimize(h->lcv(xdata, h, kernel), h0/n^2, n*h0).minimum
    end
end


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
function cvlp0(xdata::RealVector, ydata::RealVector, h::Real, kernel::Functor{3})
    n = length(ydata)
    tmp = 0.0
    w = ones(n)
    for i in 1:n
        map!(kernel, w, xdata, xdata[i], h)
        divide!(w, sum(w))
        tmp += abs2((wsum(w, ydata) - ydata[i])/(1-w[i]))
    end
    tmp/n
end

function bwlp0(xdata::RealVector, ydata::RealVector, kernel::Functor{3}=Gkernel())

  n=length(xdata)
  length(ydata)==n || error("length(ydata) != length(xdata)")
  h0= bwnormal(xdata)

  Optim.optimize(h->cvlp0(xdata, ydata, h, kernel), 0.25*h0, 10*h0).minimum
end

#see reference:Smoothing Parameter Selection in Nonparametric Regression Using an Improved Akaike Information Criterion
# Clifford M. Hurvich, Jeffrey S. Simonoff and Chih-Ling Tsai
# Journal of the Royal Statistical Society. Series B (Statistical Methodology), Vol. 60, No. 2 (1998), pp. 271-293
#http://www.jstor.org/stable/2985940
function AIClp1(xdata::RealVector, ydata::RealVector, h::Real, kernel::Functor{3})
    n = length(ydata)
    tmp = 0.0
    traceH = 0.0
    w=ones(n)
    for i in 1:n
          map!(kernel, w, xdata, xdata[i], h)
          s0 = sum(w)
          s1 = s0*xdata[i] - wsum(w, xdata)
          s2 = wsumsqdiff(w, xdata, xdata[i])
          sy0 = wsum(w, ydata)
          sy1 = NumericExtensions.wsum(w, YXdiff(), xdata, xdata[i], ydata)
          tmp+=abs2((s2 * sy0 - s1 * sy1) /(s2 * s0 - s1 * s1) - ydata[i])
          traceH += s0*w[i]/(s2 * s0 - s1 * s1)
    end
    tmp/n  + 2*(traceH+1)/(n-traceH-2)
end

function bwlp1(xdata::RealVector, ydata::RealVector, kernel::Functor{3}=Gkernel())
    n=length(xdata)
    length(ydata)==n || error("length(ydata) != length(xdata)")
    h0= bwnormal(xdata)

    Optim.optimize(h->AIClp1(xdata, ydata, h, kernel), 0.25*h0, 10*h0).minimum
end

function bwreg(xdata::RealVector, ydata::RealVector, reg::Function, kernel::Functor{3}=Gkernel())
    n=length(xdata)
    length(ydata)==n || error("length(ydata) != length(xdata)")
    h0= bwnormal(xdata)
    if reg == LP1
        return Optim.optimize(h->AIClp1(xdata, ydata, h, kernel), 0.25*h0, 10*h0).minimum
    else
        return Optim.optimize(h->cvlp0(xdata, ydata, h, kernel), 0.25*h0, 10*h0).minimum
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

