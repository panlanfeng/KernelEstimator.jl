

#Univariate goodness of test
function DiffSquare(xeval::Float64,xdata::Vector{Float64}, ydata::Vector{Float64}, 
  yParametric::Vector{Float64}, reg::Function, kernel::Function, h::Float64)
  
  reg(xeval, xdata, yParametric .- ydata, kernel, h) ^ 2
end

function Tn(xdata::Vector{Float64}, ydata::Vector{Float64}, 
  yParametric::Vector{Float64}, reg::Function, kernel::Function, h::Float64)
  
  function res(x::Vector{Float64})
    DiffSquare(x[1], xdata,ydata, yParametric, reg, kernel, h)
  end  

  pcubature(res, minimum(xdata), maximum(xdata), abstol=1e-2)[1]
end
  

#Univariate goodness of test
function BootstrapGoodness(B::Int64, xdata::Vector{Float64}, ydata::Vector{Float64}, testmodel::Function,
  reg::Function=LP0, kernel::Function=GaussianKernel, h::Float64=BandwidthLSCVReg(xdata,ydata,reg,kernel))
  
  n=length(xdata)
  TB=zeros(B)
  mhat = testmodel(xdata,ydata)
  T_value=Tn(xdata,ydata,mhat,reg,kernel,h)
  
  yhat=reg(xdata, xdata, ydata, kernel, h)

  e = ydata .- yhat
  coef=[-1,0,2]

  for b in 1:B
    boot_coef=coef[rand(Categorical([1/3, 1/2, 1/6]), n)]
    boot_e=boot_coef .* e
    boot_y=mhat .+ boot_e
    TB[b]=Tn(xdata, boot_y, testmodel(xdata, boot_y), reg, kernel, h)
  end
  
  #sum(T_value .< TB) / B
   (T_value, TB,  sum(TB .> T_value)/B)
end


#MultiDimensional Goodness of fit test 
function DiffSquare(xeval::Vector{Float64},xdata::Matrix{Float64}, ydata::Vector{Float64}, 
  yParametric::Vector{Float64}, reg::Function, kernel::Function, h::Vector{Float64})

  reg(xeval, xdata, yParametric .- ydata, kernel, h)^2
#   (reg(xeval, xdata, yParametric, kernel, h) - reg(xeval, xdata, y, kernel,h)) ^ 2
end

function Tn(xdata::Matrix{Float64}, ydata::Vector{Float64}, 
  yParametric::Vector{Float64}, reg::Function, kernel::Function, h::Vector{Float64})
  function res(x::Vector{Float64})
    DiffSquare(x, xdata,ydata, yParametric, reg, kernel, h)
  end  

  pcubature(res, minimum(xdata, 1), maximum(xdata, 1), abstol=1e-2)[1]
end
  


function BootstrapGoodness(B::Int64, xdata::Matrix{Float64}, ydata::Vector{Float64},
  testmodel::Function, reg::Function=LP0, kernel::Function=GaussianKernel, h::Vector{Float64}=BandwidthLSCVReg(xdata,ydata,reg,kernel))
  
  (n,p)=size(xdata)
  TB=zeros(B)
  mhat=testmodel(xdata, ydata)
  T_value=Tn(xdata,ydata,mhat,reg,kernel,h)
  
  # yhat to produce \epsilon; mhat to produce a less biased estimate of m(x)
  yhat=reg(xdata, xdata, ydata, kernel, h)

  e = ydata .- yhat
  coef=[-1.0,0,2]

  for b in 1:B
    print(b, "\t")
    boot_coef=coef[rand(Categorical([1/3, 1/2, 1/6]), n)]
    boot_e=boot_coef .* e
    boot_y=mhat .+ boot_e
    TB[b]=Tn(xdata, boot_y, testmodel(xdata, boot_y), reg, kernel, h)
  end
  
  (T_value, TB,  sum(TB .> T_value)/B)
end


