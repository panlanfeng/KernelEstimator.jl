#univariate kernel density
function KernelDensity(xeval::Float64, xdata::Vector{Float64}, kernel::KernelType=Gaussian, h::Float64=BandwidthLSCV(xdata, kernel))
  if h <= 0.0
    error("Bandwidth should be positive")
  end
  s0=0.0
  for i = 1:length(xdata)
    s0 += kernel.Density(xeval, xdata[i], h)
  end
  return s0 / length(xdata)
end

KernelDensity(xeval::Vector{Float64}, xdata::Vector{Float64}, kernel::KernelType=Gaussian, 
    h::Float64=BandwidthLSCV(xdata,kernel))=[KernelDensity(xeval[i], xdata, kernel, h)::Float64 for i=1:length(xeval)]


#MultiVariate kernel density
function KernelDensity(xeval::Vector{Float64}, xdata::Matrix{Float64}, 
  kernel::KernelType=Gaussian, h::Vector{Float64}=BandwidthLSCV(xdata,kernel))
  
  if any(h .<= 0)
    error("Bandwidth should be positive")
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
    s0 = s0 + kernel.Density(xeval, xi, h)
  end
  s0 / length(xdata)
end

function KernelDensity(xeval::Matrix{Float64}, xdata::Matrix{Float64}, 
  kernel::KernelType=Gaussian, h::Vector{Float64}=BandwidthLSCV(xdata,kernel))
  
  if any(h .<= 0)
    error("xeval should have same dimension as xdata")
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
    
    den[i]=KernelDensity(xi_eval, xdata, kernel, h)
  end
  den  
end


