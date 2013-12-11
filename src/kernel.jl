#univariate normal kernel
function GaussianKernel(xeval::Float64, xi::Float64, h::Float64)
  if h <= 0.0
    return Inf
  end
  exp(-((xeval - xi)/h)^2 / 2 )  / (sqrt(2pi) * h)
end

#MultiVariate Normal Kernel
function GaussianKernel(xeval::Vector{Float64}, xi::Vector{Float64}, h::Vector{Float64})
  if any(h .<= 0)
    return Inf
  end
  p=length(xi)
  tmp=1.0
  for i in 1:p
      tmp *= GaussianKernel(xeval[i], xi[i], h[i])
  end
  tmp
end

