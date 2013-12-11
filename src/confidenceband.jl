#Bootstrap confidence band for univariate nonparametric regression
function BootstrapCB(B::Int64, xeval::Vector{Float64}, 
  xdata::Vector{Float64}, ydata::Vector{Float64}, reg::Function=LP0, kernel::Function=GaussianKernel, h::Float64=BandwidthLSCVReg(xdata,ydata,reg,kernel))
  
  n=length(xdata)
  y_matrix=zeros(B, length(xeval))
  cb=zeros(2, length(xeval))
  
  yhat=reg(xdata, xdata, ydata, kernel, h)
  mhat=reg(xdata, xdata, ydata, kernel, h)
  e = ydata .- yhat
  coef=[-1,0,2]

  for b in 1:B
    boot_coef=coef[rand(Categorical([1/3, 1/2, 1/6]), n)]
    boot_e=boot_coef .* e
    boot_y=mhat .+ boot_e
    y_matrix[b, :]=reg(xeval, xdata, boot_y,kernel, h)
  end
  for i in 1:length(xeval)
    cb[:, i]=quantile(y_matrix[:, i], [.975, .025])
  end
  cb
end

#Bootstrap confidence band for multivariate nonparametric regression
function BootstrapCB(B::Int64, xeval::Matrix{Float64}, 
  xdata::Matrix{Float64}, ydata::Vector{Float64}, reg::Function=LP0, kernel::Function=GaussianKernel, h::Vector{Float64}=BandwidthLSCVReg(xdata,ydata,reg,kernel))
  
  (n, p)=size(xdata)
  (m,p1 )= size(xeval)
  if p != p1
    error("xeval should have same dimension as xdata")
  end
  
  y_matrix=zeros(B, m)
  cb=zeros(2, m)
  
  yhat=reg(xdata, xdata, ydata, kernel, h)
  mhat=reg(xdata, xdata, ydata, kernel, h)
  e = ydata .- yhat
  coef=[-1,0,2]

  for b in 1:B
#     print(b, "\t")
    boot_coef=coef[rand(Categorical([1/3, 1/2, 1/6]), n)]
    boot_e=boot_coef .* e
    boot_y=mhat .+ boot_e
    y_matrix[b, :]=reg(xeval, xdata, boot_y,kernel, h)
  end
  for i in 1:m
    cb[:, i]=quantile(y_matrix[:, i], [.975, .025])
  end
  cb
end

