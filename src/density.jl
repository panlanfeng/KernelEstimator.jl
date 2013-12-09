#univariate kernel density
function KernelDensity{T<:Float64}(xeval::Float64, xdata::Vector{T}, kernel::Function, h::Float64)
  if h <= 0.0
    return Inf
  end
  s0=0.0
  for i = 1:length(xdata)
    s0 += kernel(xeval, xdata[i], h)
  end
  fhat=s0 / length(xdata)
  fhat::Float64
end

KernelDensity{R<:Float64, T<:Float64}(xeval::Vector{R}, xdata::Vector{T}, 
kernel::Function, h::Float64)=[KernelDensity(xeval[i], xdata, kernel, h)::Float64 for i=1:length(xeval)]

KernelDensity{T<:Float64}(xeval::Float64, xdata::Vector{T}, kernel::Function, 
BWSelector::Function)=KernelDensity(xeval, xdata,kernel, BWSelector(xdata, kernel))
KernelDensity{R<:Float64,T<:Float64}(xeval::Vector{R}, xdata::Vector{T}, kernel::Function, 
BWSelector::Function)=KernelDensity(xeval, xdata,kernel, BWSelector(xdata, kernel))

#MultiVariate kernel density
function KernelDensity{R<: Float64,T<:Float64}(xeval::Vector{R}, xdata::Matrix{T}, 
  kernel::Function, h::Vector{Float64})
  
  if any(h .<= 0)
    return Inf
  end
  (n, p)=size(xdata)
  if length(h) == 1 && p==1
    return KernelDensity(xeval, xdata, kernel, h)
  end
  
  if length(xeval) != p || length(h) != p
    error("xeval should have same dimension as xdata")
  end

  
  s0=0.0
  xi=zeros(p)  
  for i = 1:n
    for j=1:p
      xi[j]=xdata[i,j]
    end
    s0 = s0 + kernel(xeval, xi, h)
  end
  s0 / length(xdata)
end

function KernelDensity{R<: Float64,T<:Float64}(xeval::Matrix{R}, xdata::Matrix{T}, 
  kernel::Function, h::Vector{Float64})
  
  if any(h .<= 0)
    return Inf
  end
  (m, p)=size(xeval)
  if length(h) != p || size(xdata)[2] != p
    error("xeval should have same dimension as xdata")
  end
  xi_eval=zeros(p)
  den=zeros(m)
  for i=1:m
    for j=1:p
      xi_eval[j]=xeval[i, j]
    end
    
    den[i]=KernelDensity(xi_eval, xdata, GaussianKernel, h)
  end
  den  
end

KernelDensity{R<: Float64,T<:Float64}(xeval::Matrix{R}, xdata::Matrix{T}, 
    kernel::Function, BandwidthSelector::Function)=KernelDensity(xeval, xdata, kernel, 
    BandwidthSelector(xdata,kernel))
KernelDensity{R<: Float64,T<:Float64}(xeval::Vector{R}, xdata::Matrix{T}, 
    kernel::Function, BandwidthSelector::Function)=KernelDensity(xeval, xdata, kernel, 
    BandwidthSelector(xdata,kernel))

