#Bootstrap confidence band for univariate nonparametric regression
function BootstrapCB{R<:Float64, T<:Float64, S<:Float64}(B::Int64, xeval::Vector{R}, 
  xdata::Vector{S}, ydata::Vector{T}, reg::Function, kernel::Function, h::Real)
  
  n=length(xdata)
  y_matrix=zeros(B, length(xeval))
  cb=zeros(2, length(xeval))
  
  yhat=reg(xdata, xdata, ydata, kernel, h)
  mhat=reg(xdata, xdata, ydata, kernel, 1.2*h)
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
function BootstrapCB{R<:Float64, T<:Float64, S<:Float64}(B::Int64, xeval::Matrix{R}, 
  xdata::Matrix{S}, ydata::Vector{T}, reg::Function, kernel::Function, h::Vector{Float64})
  
  (n, p)=size(xdata)
  (m,p1 )= size(xeval)
  if p != p1
    error("xeval should have same dimension as xdata")
  end
  
  y_matrix=zeros(B, m)
  cb=zeros(2, m)
  
  yhat=reg(xdata, xdata, ydata, kernel, h)
  mhat=reg(xdata, xdata, ydata, kernel, 1.5*h)
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

