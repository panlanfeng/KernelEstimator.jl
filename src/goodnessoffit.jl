

#Univariate goodness of test
function DiffSquare(xeval::Float64,xdata::Vector{Float64}, y::Vector{Float64}, 
  yParametric::Vector{Float64}, reg::Function, kernel::Function, h::Real)
  
  reg(xeval, xdata, yParametric .- y, kernel, h) ^ 2
end

function Tn(xdata::Vector{Float64}, y::Vector{Float64}, 
  yParametric::Vector{Float64}, reg::Function, kernel::Function, h::Real)
  
  function res(x::Float64)
    DiffSquare(x, xdata,y, yParametric, reg, kernel, h)
  end  

  pcubature(res, minimum(xdata), maximum(xdata), abstol=1e-2)[1]
end
  

#Univariate goodness of test
function BootstrapGoodness(B::Int64, xdata::Vector{Float64}, ydata::Vector{Float64},
  reg::Function, kernel::Function, h::Real,
  ParaModel::Function)
  
  n=length(xdata)
  TB=zeros(B)
  mhat = ParaModel(xdata,ydata)
  T_value=Tn(xdata,y,mhat,reg,kernel,h)
  
  yhat=reg(xdata, xdata, ydata, kernel, h)

  e = ydata .- yhat
  coef=[-1,0,2]

  for b in 1:B
    boot_coef=coef[rand(Categorical([1/3, 1/2, 1/6]), n)]
    boot_e=boot_coef .* e
    boot_y=mhat .+ boot_e
    TB[b]=Tn(xdata, boot_y, ParaModel(xdata, boot_y), reg, kernel, h)
  end
  
  #sum(T_value .< TB) / B
  (T_value, TB)
end


#MultiDimensional Goodness of fit test 
function DiffSquare(xeval::Vector{Float64},xdata::Matrix{Float64}, y::Vector{Float64}, 
  yParametric::Vector{Float64}, reg::Function, kernel::Function, h::Vector{Float64})

  reg(xeval, xdata, yParametric .- y, kernel, h)^2
#   (reg(xeval, xdata, yParametric, kernel, h) - reg(xeval, xdata, y, kernel,h)) ^ 2
end

function Tn(xdata::Matrix{Float64}, y::Vector{Float64}, 
  yParametric::Vector{Float64}, reg::Function, kernel::Function, h::Vector{Float64})
  function res(x::Vector{Float64})
    DiffSquare(x, xdata,y, yParametric, reg, kernel, h)
  end  

  pcubature(res, minimum(xdata, 1), maximum(xdata, 1), abstol=1e-2)[1]
end
  


function BootstrapGoodness(B::Int64, xdata::Matrix{Float64}, ydata::Vector{Float64},
  reg::Function, kernel::Function, h::Vector{Float64}, testmodel::Function)
  
  (n,p)=size(xdata)
  TB=zeros(B)
  mhat=testmodel(xdata, ydata)
  T_value=Tn(xdata,y,mhat,reg,kernel,h)
  
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
  
  #sum(T_value .< TB) / B
  (T_value, TB)
end


