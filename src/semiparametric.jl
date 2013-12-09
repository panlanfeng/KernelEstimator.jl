
#Single index model
function Sind(xdata::Matrix{Float64}, ydata::Vector{Float64}, kernel::Function)
  (n,p)=size(xdata)
  beta0=zeros(p)
  beta0=linreg(xdata,ydata)[2:(p+1)]
#   beta0[p+1]=BandwidthLSCVReg(xdata * beta0[1:p],y,LP0,kernel)
  
  function ss(x::Vector{Float64})
    xb=xdata*x
    leave_one_out(xb,ydata,LP0,kernel,[1.0])  
  end
  optimize(ss, beta0, iterations=100).minimum  
end

