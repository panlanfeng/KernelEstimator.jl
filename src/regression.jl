#univariate nadaraya-watson estimate
function LP0(xeval::Real, xdata::RealVector, ydata::RealVector; kernel::Functor{3}=Gkernel(), h::Real=bwlp0(xdata,ydata,kernel))
    n=length(xdata)
    length(ydata) == n || error("length of ydata not the same with xdata")
    w=ones(n)
    map!(kernel,w, xdata, xeval, h)
    wsum(w, ydata)/sum(w)
end
# LP0(xeval::Real, xdata::RealVector, ydata::RealVector; kernel::Functor{3}=Gkernel(), h::Real=bwlp0(xdata,ydata,kernel)) = LP0(float(xeval), float(xdata), float(ydata), h=float(h), kernel=kernel)
LP0(xeval::RealVector, xdata::RealVector, ydata::RealVector; kernel::Functor{3}=Gkernel(), h::Real=bwlp0(xdata,ydata,kernel))=Float64[LP0(xeval[i], xdata,ydata,h=h, kernel=kernel) for i=1:length(xeval)]
LP0(xdata::RealVector, ydata::RealVector; kernel::Functor{3}=Gkernel(), h::Real=bwlp0(xdata,ydata,kernel))=LP0(xdata, xdata, ydata, kernel=kernel,h=h)

yxdiff{T<:FloatingPoint}(xi::T, xeval::T, y::T)=y*(xeval - xi)
yxdiff{T<:Real}(xi::T, xeval::T, y::T)=yxdiff(float(xi), float(xeval), float(y))
yxdiff(xi::Real, xeval::Real, y::Real)=yxdiff(promote(xi, xeval, y)...)
type YXdiff <: Functor{3} end
NumericExtensions.evaluate(::YXdiff, xi, xeval, y) = yxdiff(xi, xeval, y)


##univariate local linear
function LP1(xeval::Real, xdata::RealVector, ydata::RealVector; kernel::Functor{3}=Gkernel(), h::Real=bwlp1(xdata, ydata, kernel))
    n=length(xdata)
    length(ydata) == n || error("length of ydata not the same with xdata")
    w = ones(n)
    map!(kernel, w, xdata, xeval, h)
    s0 = sum(w)
    s1 = s0*xeval - wsum(w, xdata)
    s2 = wsumsqdiff(w, xdata, xeval)
    sy0 = wsum(w, ydata)
    sy1 = NumericExtensions.wsum(w, YXdiff(), xdata, xeval, ydata)

    (s2 * sy0 - s1 * sy1) /(s2 * s0 - s1 * s1)

end
LP1(xeval::RealVector, xdata::RealVector, ydata::RealVector;kernel::Functor{3}=Gkernel(), h::Real=bwlp1(xdata, ydata, kernel))=Float64[LP1(xeval[i], xdata,ydata,h=h, kernel=kernel) for i=1:length(xeval)]
LP1(xdata::RealVector, ydata::RealVector;kernel::Functor{3}=Gkernel(), h::Real=bwlp1(xdata, ydata, kernel))=LP1(xdata, xdata,ydata,kernel=kernel, h=h)



# #multi-variate nadaraya-watson
# function LP0(xeval::Vector{Float64}, xdata::Matrix{Float64}, ydata::Vector{Float64}, kernel::Function=GaussianKernel, h::Vector{Float64}=BandwidthLSCVReg(xdata,ydata,LP0,kernel))

#   (n,p)=size(xdata)
#   if length(xeval) != p || length(h) !=p
#     error("xeval, xdata and h should have same dimension!")
#   end

#   tmp=zeros(n)
#   for i in 1:n
#     tmp[i]=prod([GaussianKernel(xeval[j], xdata[i,j],h[j])::Float64 for j in 1:p])
#   end

#   s0 = sum(tmp)
#   sy0 = sum(tmp .* ydata)
#   sy0 / s0
# end

# #
# function LP0(xeval::Matrix{Float64}, xdata::Matrix{Float64},
#   ydata::Vector{Float64}, kernel::Function=GaussianKernel, h::Vector{Float64}=BandwidthLSCVReg(xdata,ydata,LP0,kernel))

#   (m,p)=size(xeval)
#   den=zeros(m)
#   xi=zeros(p)
#   for i=1:m
#     for k in 1:p
#       xi[k] = xeval[i,k]
#     end
#     den[i] = LP0(xi, xdata, ydata, kernel, h)::Float64
#   end
#   den
# end


