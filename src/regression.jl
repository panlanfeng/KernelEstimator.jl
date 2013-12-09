#univariate nadaraya-watson estimate
function LP0{T <: Float64, S<:Float64}(xeval::Float64, xdata::Array{T, 1}, 
  ydata::Array{S, 1}, kernel::Function, h::Float64)
  n=length(xdata)
  s0=0.0
  sy0=0.0
  for i in 1:n
    tmp=kernel(xeval, xdata[i], h)::Float64
    s0 += tmp
    sy0 += tmp * ydata[i]
  end
  fhat=sy0 / s0
  fhat::Float64
end

LP0{R<:Float64, T <: Float64, S<:Float64}(xeval::Vector{R}, xdata::Vector{T}, 
ydata::Vector{S}, kernel::Function, h::Float64)=[LP0(xeval[i], xdata,ydata,kernel,
  h)::Float64 for i=1:length(xeval)]

LP0{T <: Float64, S<:Float64}(xeval::Float64, xdata::Vector{T}, ydata::Vector{S}, 
kernel::Function, BWSelector::Function)=LP0(xeval, xdata, ydata, kernel, 
BWSelector(xdata,ydata,LP0,kernel)::Float64)
LP0{R<:Float64,T <: Float64, S<:Float64}(xeval::Vector{R}, xdata::Vector{T}, ydata::Vector{S}, 
kernel::Function, BWSelector::Function)=LP0(xeval, xdata, ydata, kernel, 
BWSelector(xdata,ydata,LP0, kernel)::Float64)


##univariate local linear
function LP1{T <: Float64, S<:Float64}(xeval::Float64, xdata::Vector{T}, ydata::Vector{S}, 
kernel::Function, h::Float64)
  
  n=length(xdata)
  s0=0.0
  s1=0.0
  s2=0.0
  sy0=0.0
  sy1=0.0
  for i in 1:n
    tmp=kernel(xeval, xdata[i], h)::Float64
    s0 += tmp
    s1 += tmp * (xeval - xdata[i])
    s2 += tmp * (xeval - xdata[i])^2
    sy0 += tmp * ydata[i]
    sy1 += tmp * (xeval - xdata[i]) * ydata[i]
  end
    
  fhat=(s2 * sy0 - s1 * sy1) /(s2 * s0 - s1 * s1)
  fhat::Float64
end
LP1{R<:Float64,T <: Float64, S<:Float64}(xeval::Vector{R}, xdata::Array{T, 1}, ydata::Array{S, 1}, 
kernel::Function, h::Float64)=[LP1(xeval[i], xdata,ydata,kernel,h)::Float64 for i=1:length(xeval)]

LP1{T <: Float64, S<:Float64}(xeval::Float64, xdata::Array{T, 1}, ydata::Array{S, 1}, 
kernel::Function, BWSelector::Function)=LP1(xeval,xdata,ydata,kernel,BWSelector(xdata,ydata,LP1,kernel))
LP1{R<:Float64,T <: Float64, S<:Float64}(xeval::Vector{R}, xdata::Array{T, 1}, ydata::Array{S, 1}, 
kernel::Function, BWSelector::Function)=LP1(xeval,xdata,ydata,kernel,BWSelector(xdata,ydata,LP1,kernel))


#multi-variate nadaraya-watson
function LP0(xeval::Vector{Float64}, xdata::Matrix{Float64}, 
  ydata::Vector{Float64}, kernel::Function, h::Vector{Float64})
    
  (n,p)=size(xdata)
  if length(xeval) != p || length(h) !=p
    error("xeval, xdata and h should have same dimension!")
  end

  tmp=zeros(n)
  for i in 1:n
    tmp[i]=prod([GaussianKernel(xeval[j], xdata[i,j],h[j])::Float64 for j in 1:p])
  end

  s0 = sum(tmp)
  sy0 = sum(tmp .* ydata)
  sy0 / s0
end

#
function LP0(xeval::Matrix{Float64}, xdata::Matrix{Float64}, 
  ydata::Vector{Float64}, kernel::Function, h::Vector{Float64})

  (m,p)=size(xeval)
  den=zeros(m)
  xi=zeros(p)
  for i=1:m
    for k in 1:p
      xi[k] = xeval[i,k]
    end
    den[i] = LP0(xi, xdata, ydata, kernel, h)::Float64
  end
  den
end


LP0{R<:Float64, T<:Float64, S<:Float64}(xeval::Vector{R}, xdata::Matrix{T}, 
    ydata::Vector{S}, kernel::Function, BandwidthSelector::Function)=LP0(xeval, xdata, ydata, kernel, 
BandwidthSelector(xdata, ydata, LP0, kernel))
LP0{R<:Float64, T<:Float64, S<:Float64}(xeval::Matrix{R}, xdata::Matrix{T}, 
    ydata::Vector{S}, kernel::Function, BandwidthSelector::Function)=LP0(xeval, xdata, ydata, kernel, 
    BandwidthSelector(xdata, ydata, LP0,kernel))

