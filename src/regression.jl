#univariate nadaraya-watson estimate

function localconstant(xdata::AbstractVector{<:Real}, ydata::AbstractVector{<:Real}; xeval::AbstractVector{<:Real}=xdata, kernel::Function=gaussiankernel, h::Real=bwlocalconstant(xdata,ydata,kernel))

    n=length(xdata)
    length(ydata) == n || error("length(ydata) != length(xdata)")
    w=ones(n)
    pre = zeros(length(xeval))
    for i in 1:length(xeval)
        kernel(xeval[i], xdata, h, w, n)
        pre[i] = wsum(w, ydata)/sum(w)
    end
    pre
end
localconstant(xdata::AbstractVector{<:Real}, ydata::AbstractVector{<:Real}, xeval::Real; kernel::Function = gaussiankernel, h::Real = bwlocalconstant(xdata,ydata,kernel)) = localconstant(xdata, ydata, xeval = [xeval;], kernel=kernel, h = h)

function wsumsqdiff(w::AbstractVector{<:Real}, xdata::AbstractVector{<:Real}, xeval::Real, n::Int)
    res = 0.0
    for i in 1:n
        @inbounds res += w[i]*(xdata[i]-xeval).^2
    end
    res
end
function wsumyxdiff(w::AbstractVector{<:Real}, xdata::AbstractVector{<:Real}, xeval::Real, ydata::AbstractVector{<:Real}, n::Int)
    res = 0.0
    for i in 1:n
        @inbounds res += w[i]*ydata[i]*(xeval-xdata[i])
    end
    res
end

##univariate local linear
function locallinear(xdata::AbstractVector{<:Real}, ydata::AbstractVector{<:Real}; xeval::AbstractVector{<:Real}=xdata, kernel::Function=gaussiankernel, h::Real=bwlocallinear(xdata, ydata, kernel))
    n=length(xdata)
    length(ydata) == n || error("length of ydata not the same with xdata")
    w = ones(n)
    pre = zeros(length(xeval))
    for i in 1:length(xeval)
        kernel(xeval[i], xdata, h, w, n)
        s0 = sum(w)
        s1 = s0*xeval[i] - wsum(w, xdata)
        s2 = wsumsqdiff(w, xdata, xeval[i], n)
        sy0 = wsum(w, ydata)
        sy1 = wsumyxdiff(w, xdata, xeval[i], ydata, n)
        pre[i] = (s2 * sy0 - s1 * sy1) /(s2 * s0 - s1 * s1)
    end
    pre
end
locallinear(xdata::AbstractVector{<:Real}, ydata::AbstractVector{<:Real}, xeval::Real; kernel::Function = gaussiankernel, h::Real = bwlocallinear(xdata,ydata,kernel)) = locallinear(xdata, ydata, xeval=[xeval;], kernel=kernel, h=h)



function boundit(xdata::AbstractVector{<:Real}, xeval::AbstractVector{<:Real}, kernel::Function, lb::Real, ub::Real)
    if (lb == -Inf) && (ub == Inf)
        return (xdata, xeval, kernel)
    elseif (lb > -Inf) && (ub < Inf)
        (all(lb .<= xeval .<= ub) && all(lb .<= xdata .<= ub)) || error("Your data are not in [lb,ub]")
        xeval = (xeval .- lb)./(ub - lb)
        xdata = (xdata .- lb)./(ub - lb)
        if kernel != betakernel
            warn("Kernel is set to betakernel")
            kernel = betakernel
        end
    elseif (lb > -Inf) && (ub == Inf)
        (all(xeval .>= lb) && all(xdata .>= lb)) || error("lb should be less than your data")
        xeval = xeval .- lb
        xdata = xdata .- lb
        if kernel != gammakernel
            warn("Kernel is set to gammakernel")
            kernel = gammakernel
        end
    elseif (lb == -Inf) && (ub < Inf)
        (all(xeval .<= ub) && all(xdata .<= ub)) || error("ub should be greater than your data")
        xeval = ub .- xeval
        xdata = ub .- xdata
        if kernel != gammakernel
             warn("Kernel is set to gamma kernel")
             kernel = gammakernel
        end
    end
    (xdata, xeval, kernel)
end

function npr(xdata::AbstractVector{<:Real}, ydata::AbstractVector{<:Real}; xeval::AbstractVector{<:Real}=xdata,
        reg::Function=locallinear, lb::Real=-Inf, ub::Real=Inf, kernel::Function=gaussiankernel, h::Real=-Inf)

    xdata, xeval, kernel = boundit(xdata, xeval, kernel, lb, ub)
    if h <= 0
        h = bwreg(xdata, ydata, reg, kernel)
    end
    reg(xdata, ydata, xeval=xeval, kernel=kernel, h=h)
end


# #multi-variate nadaraya-watson regression or local linear
function localconstant(xdata::AbstractMatrix{<:Real}, ydata::AbstractVector{<:Real}; kernel::Array{Function, 1}=Function[gaussiankernel for i in 1:size(xdata)[2]], xeval::AbstractMatrix{<:Real}=xdata,  h::AbstractVector{<:Real}=bwlocalconstant(xdata, ydata, kernel))

    m, p = size(xeval)
    n, p1 = size(xdata)
    if p1 != p || length(h) !=p
        error("xeval, xdata and h should have same dimension!")
    end
    if any(h.<=0.0)
        error("Bandwidth should be positive")
    end
    pre=zeros(m)

    for i=1:m
        w = ones(n)
        wtmp = ones(n)
        for j in 1:p
            kernel[j](xeval[i, j], xdata[:, j], h[j], wtmp, n)
            for k in 1:n
                w[k] = w[k] * wtmp[k]
            end
        end
        divide!(w, sum(w))
        pre[i] = wsum(w, ydata)
    end
    pre
end

lp0=localconstant
lp1=locallinear
