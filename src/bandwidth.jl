#bandwidth selector
function BandwidthNormalReference{T<:Float64}(xdata::Vector{T}, kernel::Function)
  h=1.06 * min((quantile(xdata, .75) - quantile(xdata, .25)) / 1.34, std(xdata)) * length(xdata) ^ (-0.2)
  h::Float64
end


#Leave-one-out cross validation. Currently only work for GaussianKernel
#\sum_{i,j}K'((x_i - x_j)/h) /(n^2*h) + 2K(0)/(nh)
#where K' = K^{(2)} - 2K
#When Kernel is normal it is 
function BandwidthLSCV{T<:Float64}(xdata::Vector{T}, kernel::Function)
  
  n=length(xdata)
  h0=BandwidthNormalReference(xdata, kernel)
  
  function res(hv::Vector{Float64})  
    h=hv[1]
    tmp=0.0
    for i in 1:n
      for j in 1:n
        if j == i
          continue
        end
        #tmp += kernel(xdata[i], xdata[j], sqrt(2)*h) - 2*kernel(xdata[i], xdata[j], h)
        xdiff = ((xdata[i] - xdata[j]) / h)^2
        tmp += (2^(-1/2)*exp(-xdiff/4) - 2*exp(-xdiff/2))::Float64
      end
    end
    (tmp/(n^2) + 2^(-.5)/n)/h
  end
  optimize(res, [h0], iterations=100).minimum[1]  
end

#Leave-one-out cross validation. Currently only work for GaussianKernel
#\sum_{i,j}K'((x_i - x_j)/h) /(n^2*h) + 2K(0)/(nh)
#where K' = K^{(2)} - 2K
#When Kernel is normal it is 
function BandwidthLSCV{T<:Float64}(xdata::Matrix{T}, kernel::Function) 

  (n, p)=size(xdata)  
  h0=BandwidthNormalReference(reshape(xdata[:,1], n), GaussianKernel)  
  function res(h::Vector{Float64})
    tmp=0.0
    xi=zeros(p)
    xj=zeros(p)
    for i in 1:n      
      for j in 1:n
        if j == i
          continue
        end        
        xdiff=sum([(((xdata[i,k]-xdata[j,k]) / h[k]) ^ 2)::Float64 for k in 1:p]) 
        tmp += exp(-xdiff / 4) - 2^(1+p/2) * exp(-xdiff/2)
      end
    end
    (tmp / n^2  + 1/n) / prod(h)
  end  
  optimize(res, [h0 for i in 1:p], iterations=100).minimum
end

function leave_one_out{T <: Float64, S<:Float64}(xdata::Vector{T}, ydata::Vector{S},reg::Function, 
  kernel::Function, hv::Vector{Float64})
  h=hv[1]
  n=length(ydata) 
  ls=0.0
  for i in 1:n
    ls += (ydata[i] - reg(xdata[i], xdata[[1:(i-1), (i+1):end]], 
  ydata[[1:(i-1), (i+1):end]], kernel, h)[1])^2
  end
  ls = ls / n
  ls::Float64
end

#leavee-one-out LSCV. Using 
function BandwidthLSCVReg{T <: Float64, S<:Float64}(xdata::Vector{T}, ydata::Vector{S}, 
  reg::Function, kernel::Function)

  n=length(xdata)
  h0= BandwidthNormalReference(xdata, kernel)
  function res(x::Vector{Float64})
    leave_one_out(xdata,ydata,reg,kernel,x)
  end
  optimize(res, [h0], iterations=100).minimum[1]
end


function leave_one_out{T <: Float64, S<:Float64, R<:Float64}(xdata::Matrix{T}, ydata::Vector{S},reg::Function, 
  kernel::Function, h::Vector{R})
  (n,p)=size(xdata) 
  
  yhat=[reg(xdata[i, :], xdata[[1:(i-1), (i+1):end], :], 
    ydata[[1:(i-1), (i+1):end]], kernel, h)[1]::Float64 for i in 1:n]
  NumericExtensions.vnormdiff(ydata,yhat)  
end

#leave-one-out LSCV for multivariate NW
function BandwidthLSCVReg{T <: Float64, S<:Float64}(xdata::Matrix{T}, ydata::Vector{S}, reg::Function, 
  kernel::Function)

  (n,p)=size(xdata)  
  h0=BandwidthNormalReference(reshape(xdata[:,1], n), GaussianKernel)
  
  function res(x::Vector{Float64})
    leave_one_out(xdata,ydata,reg,kernel,x)  
  end
  optimize(res, [h0 for i in 1:p], iterations=100).minimum
end

