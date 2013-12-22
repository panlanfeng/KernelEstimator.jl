
#Single index model
function Sind(xdata::Matrix{Float64}, ydata::Vector{Float64}, kernel::Function=GaussianKernel)
  (n,p)=size(xdata)
  beta0=zeros(p)
  beta0=linreg(xdata,ydata)[2:(p+1)]
#   beta0[p+1]=BandwidthLSCVReg(xdata * beta0[1:p],y,LP0,kernel)
  
  function ss(x::Vector{Float64})
    xb=xdata*x
    ls=0.0
    for i in 1:n
      ls += (ydata[i] - LP0(xb[i], xb[[1:(i-1), (i+1):end]], ydata[[1:(i-1), (i+1):end]], kernel, [1.0]))^2
    end
    ls / n
  end
  optimize(ss, beta0, iterations=100).minimum  
end

