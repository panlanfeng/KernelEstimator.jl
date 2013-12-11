#bandwidth selector
function BandwidthNormalReference(xdata::Vector{Float64})
  1.06 * min((quantile(xdata, .75) - quantile(xdata, .25)) / 1.34, std(xdata)) * length(xdata) ^ (-0.2)
end


#Leave-one-out cross validation. Currently only work for GaussianKernel
#\sum_{i,j}K'((x_i - x_j)/h) /(n^2*h) + 2K(0)/(nh)
#where K' = K^{(2)} - 2K
#When Kernel is normal it is 
function BandwidthLSCV(xdata::Vector{Float64}, kernel::Function=GaussianKernel)
  
  n=length(xdata)
  h0=BandwidthNormalReference(xdata)
  
  if kernel==GaussianKernel
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
          tmp += (2^(-1/2)*exp(-xdiff/4) - 2*exp(-xdiff/2))
        end
      end
      (tmp/(n^2) + 2^(-.5)/n)/h
    end
    return optimize(res, [h0], iterations=100).minimum[1] + .1/n  # add a lower order item to avoid 0 bandwidth
  else
      error("Currently only GaussianKernel is supported for bandwidth selecton.")
  end
end

#Leave-one-out cross validation. Currently only work for GaussianKernel
#\sum_{i,j}K'((x_i - x_j)/h) /(n^2*h) + 2K(0)/(nh)
#where K' = K^{(2)} - 2K
#When Kernel is normal it is 
function BandwidthLSCV(xdata::Matrix{Float64}, kernel::Function=GaussianKernel) 

  (n, p)=size(xdata)  
  h0=BandwidthNormalReference(reshape(xdata[:,1], n))  
  if kernel==GaussianKernel
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
    return optimize(res, [h0 for i in 1:p], iterations=100).minimum .+ .1/n
  else 
      error("Currently only GaussianKernel is supported for bandwidth selection.")
  end
end

function leave_one_out(xdata::Vector{Float64}, ydata::Vector{Float64},reg::Function, 
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
function BandwidthLSCVReg(xdata::Vector{Float64}, ydata::Vector{Float64}, reg::Function=LP0, kernel::Function=GaussianKernel)

  n=length(xdata)
  h0= BandwidthNormalReference(xdata)
  function res(x::Vector{Float64})
    leave_one_out(xdata,ydata,reg,kernel,x)
  end
  optimize(res, [h0], iterations=100).minimum[1] + .1/n
end


function leave_one_out(xdata::Matrix{Float64}, ydata::Vector{Float64},reg::Function, 
  kernel::Function, h::Vector{Float64})
  (n,p)=size(xdata) 
  
  yhat=[reg(xdata[i, :], xdata[[1:(i-1), (i+1):end], :], 
    ydata[[1:(i-1), (i+1):end]], kernel, h)[1]::Float64 for i in 1:n]
  NumericExtensions.vnormdiff(ydata,yhat)  
end

#leave-one-out LSCV for multivariate NW
function BandwidthLSCVReg(xdata::Matrix{Float64}, ydata::Vector{Float64}, reg::Function=LP0, kernel::Function=GaussianKernel)

  (n,p)=size(xdata)  
  h0=BandwidthNormalReference(reshape(xdata[:,1], n))
  
  function res(x::Vector{Float64})
    leave_one_out(xdata,ydata,reg,kernel,x)  
  end
  optimize(res, [h0 for i in 1:p], iterations=100).minimum .+ .1 / n

end

