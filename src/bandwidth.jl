#bandwidth selector
function BandwidthNormalReference(xdata::Vector{Float64})
  1.06 * min((quantile(xdata, .75) - quantile(xdata, .25)) / 1.34, std(xdata)) * length(xdata) ^ (-0.2)
end


#Leave-one-out cross validation. Currently only work for GaussianKernel
#\sum_{i,j}K'((x_i - x_j)/h) /(n^2*h) + 2K(0)/(nh)
#where K' = K^{(2)} - 2K
#When Kernel is normal it is 
function BandwidthLSCV(xdata::Vector{Float64}, kernel::KernelType=Gaussian)
  
    n=length(xdata)
    h0=BandwidthNormalReference(xdata)
    function res(hv::Vector{Float64})  
      h=hv[1]
      tmp1=0.0
      tmp2=0.0
      for i in 1:n
        for j in 1:n
          if j == i
            continue
          end
          xdiff=(xdata[i] - xdata[j]) / h
          tmp1 += kernel.Convolution(xdiff)
          tmp2 += kernel.Density(xdiff,0.0,1.0)
        end
      end

      (tmp1 / (n ^ 2) - tmp2 / (n * (n - 1)) * 2 + kernel.Convolution(0.0)/n)/h
    end
    return optimize(res, [h0], iterations=100).minimum[1] + .1/n  # add a lower order item to avoid 0 bandwidth
end

#Leave-one-out cross validation. Currently only work for GaussianKernel
#\sum_{i,j}K'((x_i - x_j)/h) /(n^2*h) + 2K(0)/(nh)
#where K' = K^{(2)} - 2K
#When Kernel is normal it is 
function BandwidthLSCV(xdata::Matrix{Float64}, kernel::KernelType=Gaussian) 
    (n, p)=size(xdata)  
    h0=BandwidthNormalReference(reshape(xdata[:,1], n))  
    function res(h::Vector{Float64})  
      tmp1=0.0
      tmp2=0.0
      for i in 1:n
        for j in 1:n
          if j == i
            continue
          end
          #tmp += kernel(xdata[i], xdata[j], sqrt(2)*h) - 2*kernel(xdata[i], xdata[j], h)
          #xdiff = ((xdata[i,:] .- xdata[j,:]) ./ h)^2
          #tmp += (2^(-1/2)*exp(-xdiff/4) - 2*exp(-xdiff/2))
          xdiff=[((xdata[i,k] - xdata[j,k]) / h[k])::Float64 for k in 1:p]
          tmp1 += kernel.Convolution(xdiff)
          tmp2 += kernel.Density(xdiff,zeros(p),ones(p))
        end
      end
      (tmp1 / (n ^ 2) - tmp2 / (n * (n - 1)) * 2 + kernel.Convolution(zeros(p))/n)/prod(h)

    end
   return optimize(res, [h0 for i in 1:p], iterations=100).minimum .+ .1/n
end

#leave-one-out LSCV. Using 
function BandwidthLSCVReg(xdata::Vector{Float64}, ydata::Vector{Float64}, reg::Function=LP0, kernel::Function=GaussianKernel)

  n=length(xdata)
  h0= BandwidthNormalReference(xdata)
  function res(h::Vector{Float64})
    ls=0.0
    for i in 1:n
      ls += (ydata[i]-reg(xdata[i], xdata[[1:(i-1), (i+1):end]],ydata[[1:(i-1), (i+1):end]], kernel, h[1]))^2
    end
    ls / n
  end
  optimize(res, [h0], iterations=100).minimum[1] + .1/n
end


#leave-one-out LSCV for multivariate NW
function BandwidthLSCVReg(xdata::Matrix{Float64}, ydata::Vector{Float64}, reg::Function=LP0, kernel::Function=GaussianKernel)

  (n,p)=size(xdata)  
  h0=BandwidthNormalReference(reshape(xdata[:,1], n))
  
  function res(x::Vector{Float64})
      ls=0.0
      for i in 1:n
          ls += (ydata[i]-reg(xdata[i,:],xdata[[1:(i-1),(i+1):end],:], ydata[[1:(i-1),(i+1):end]],kernel,x))^2
      end
      ls/n
  end
  optimize(res, [h0 for i in 1:p], iterations=100).minimum .+ .1 / n

end

