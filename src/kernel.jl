#univariate normal kernel
function GaussianKernel(xeval::Float64, xi::Float64, h::Float64)
  if h <= 0.0
    return Inf
  end
  exp(-((xeval - xi)/h)^2 / 2 )  / (sqrt(2*pi) * h)
end

#MultiVariate Normal Kernel
function GaussianKernel(xeval::Vector{Float64}, xi::Vector{Float64}, h::Vector{Float64})
  if any(h .<= 0)
    return Inf
  end
  p=length(xi)
  exp(-sumsq((xeval.-xi)./h)/2) / (2*pi)^(p/2) / prod(h)
end

