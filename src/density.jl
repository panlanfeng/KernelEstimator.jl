


#univariate kernel density
function kde(xeval::Real, xdata::RealVector; kernel::Functor{3}=Gkernel(), h::Real=bwcv(xdata, kernel))
    h > 0.0 || error("Bandwidth should be positive")
    mean(kernel, xdata, xeval, h)
end
# kde{T<:Real}(xeval::Real, xdata::Vector{T}, h::Real;kernel::Functor{3}=Gkernel())=kde(float(xeval), float(xdata), float(h),kernel)

kde(xeval::RealVector, xdata::RealVector; kernel::Functor{3}=Gkernel(), h::Real=bwcv(xdata, kernel))=Float64[kde(xeval[i], xdata, h=h, kernel=kernel) for i=1:length(xeval)]
kde(xdata::RealVector;kernel::Functor{3}=Gkernel(), h::Real=bwcv(xdata, kernel)) = kde(xdata, xdata, h=h, kernel=kernel)

# #MultiVariate kernel density
# function KernelDensity(xeval::Vector{Float64}, xdata::Matrix{Float64},
#   kernel::KernelType=Gaussian, h::Vector{Float64}=BandwidthLSCV(xdata,kernel))


#   if any(h .<= 0)
#     error("Bandwidth should be positive")
#   end
#   (n, p)=size(xdata)
#   if length(h) == 1 && p==1
#     return KernelDensity(xeval, xdata, kernel, h)
#   end

#   if length(xeval) != p || length(h) != p
#     error("xeval should have same dimension as xdata")
#   end


#   s0=0.0
#   xi=zeros(p)
#   for i = 1:n
#     for j=1:p
#       xi[j]=xdata[i,j]
#     end
#     s0 = s0 + kernel.Density(xeval, xi, h)
#   end
#   s0 / length(xdata)
# end

# function KernelDensity(xeval::Matrix{Float64}, xdata::Matrix{Float64},
#   kernel::KernelType=Gaussian, h::Vector{Float64}=BandwidthLSCV(xdata,kernel))

#   if any(h .<= 0)
#     error("xeval should have same dimension as xdata")
#   end
#   (m, p)=size(xeval)
#   if length(h) != p || size(xdata)[2] != p
#     error("xeval should have same dimension as xdata")
#   end
#   xi_eval=zeros(p)
#   den=zeros(m)
#   for i=1:m
#     for j=1:p
#       xi_eval[j]=xeval[i, j]
#     end

#     den[i]=KernelDensity(xi_eval, xdata, kernel, h)
#   end
#   den
# end
