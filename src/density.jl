
function kerneldensity(xdata::RealVector; xeval::RealVector=xdata, lb::Real=-Inf, ub::Real=Inf, kernel::Function=gaussiankernel, h::Real=-Inf)

    xdata, xeval, kernel = boundit(xdata, xeval, kernel, lb, ub)
    if h <= 0
        h = bwlscv(xdata, kernel)
    end
    n = length(xdata)
    w = zeros(n)
    den = zeros(length(xeval))
    for i in 1:length(xeval)
        kernel(xeval[i], xdata, h, w, n)
        den[i] = mean(w)
    end
    return den
end

function kerneldensity(xdata::RealMatrix; xeval::RealMatrix=xdata, 
    kernel::Array{Function, 1}=[gaussiankernel for i in 1:size(xdata)[2]], h::RealVector=bwlcv(xdata, kernel))
    
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
            for k in 1:n
                w[k] = w[k] * wtmp[k]
            end
        end
        den[i] = mean(w)
    end
    den
end

kde=kerneldensity
