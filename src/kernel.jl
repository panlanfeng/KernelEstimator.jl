#univariate normal kernel
function GaussianKernel(xeval::Float64, xi::Float64, h::Float64)
  if h <= 0.0
    return Inf
  end
  fhat=exp(-((xeval - xi)/h)^2 / 2 )  / (sqrt(2pi) * h)
  fhat::Float64
end

#MultiVariate Normal Kernel
function GaussianKernel(xeval::Vector{Float64}, xi::Vector{Float64}, h::Vector{Float64})
  if any(h .<= 0)
    return Inf
  end
  p=length(xi)
  prod([GaussianKernel(xeval[i], xi[i],h[i])::Float64 for i in 1:p])

end

