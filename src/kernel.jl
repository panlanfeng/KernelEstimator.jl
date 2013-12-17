type KernelType
    Density::Function
    Convolution::Function
end
#univariate normal kernel
function GaussianKernel(xeval::Float64, xi::Float64, h::Float64)
  if h <= 0.0
    return Inf
  end
  exp(-((xeval - xi)/h)^2 / 2 )  / sqrt(2 * pi) / h
end

#MultiVariate Normal Kernel
function GaussianKernel(xeval::Vector{Float64}, xi::Vector{Float64}, h::Vector{Float64})
  if any(h .<= 0.0)
    return Inf
  end
  p=length(xi)

  exp(-wsumsqdiff(1./(h.^2), xeval, xi)) / (2*pi)^(p/2) / prod(h)

end
GaussianKernel2(x::Float64) = GaussianKernel(x, 0.0, sqrt(2))
GaussianKernel2(x::Vector{Float64}) = GaussianKernel(x, zeros(length(x)), sqrt(2).*ones(length(x)))

Gaussian=KernelType(GaussianKernel, GaussianKernel2)

function EKernel(xeval::Float64, xi::Float64,h::Float64)
    if h <= 0.0
      return Inf
    end
    u=(xeval - xi) / h
    if abs(u) > 1.0
        return 0.0
    end
    return (1.0-u*u)*3.0/4.0 
end


function EKernel(xeval::Vector{Float64}, xi::Vector{Float64}, h::Vector{Float64})

    if any(h .<= 0.0)
        return Inf
    end
    p=length(xi)
    tmp=1.0
    for i in 1:p
        tmp *= EKernel(xeval[i],xi[i],h[i])
    end
    tmp
end

function EKernel2(x::Float64)
    ax=abs(x)
    if ax > 2.0
        return 0.0
    end

    (2 - ax)^3*(ax^2+6*ax+4) * 3 / 160

end  
function EKernel2(x::Vector{Float64})
    ax=abs(x)
    if any(ax .> 2.0)
        return 0.0
    end
    p=length(ax)
    tmp=1.0
    for i in 1:p
        tmp *= EKernel2(ax[i])
    end
    tmp
end  
Epanechnikov=KernelType(EKernel, EKernel2)


