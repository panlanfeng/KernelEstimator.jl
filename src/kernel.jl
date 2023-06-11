rhoxb(x::Real, b::Real) = 2*b*b + 2.5 - sqrt(4*b^4 + 6*b*b+2.25 - x*x - x/b)

function multiply!(des::AbstractVector{<:Real}, x::AbstractVector{<:Real}, y::Real, n::Int=length(x))
    for i in 1:n
        @inbounds des[i] = x[i]*y
    end
end
multiply!(x::AbstractVector{<:Real}, y::Real) = multiply!(x, x, y)
function divide!(des::AbstractVector{<:Real}, x::AbstractVector{<:Real}, y::Real, n::Int=length(x))
    for i in 1:n
        @inbounds des[i] = x[i]/y
    end
end
divide!(x::AbstractVector{<:Real}, y::Real) = divide!(x, x, y)
function minus!(des::AbstractVector{<:Real}, y::Float64, x::AbstractVector{<:Real}, n::Int64=length(x))
   for i in 1:n
       @inbounds des[i] = y - x[i]
   end
   nothing
end
function add!(x::Vector{Float64}, y::Float64, n::Int64=length(x))
   for i in 1:n
       @inbounds x[i] = x[i] + y
   end
   nothing
end
function abs2!(des::AbstractVector{<:Real}, x::AbstractVector{<:Real}, n::Int64=length(x))
   for i in 1:n
       @inbounds des[i] = abs2(x[i])
   end
   nothing
end
function betakernel(x::Real, xdata::AbstractVector{<:Real}, h::Real, w::Vector, n::Int)
    a = x / h - 1
    b = (1 - x) / h - 1
    if (x < 0) || (x > 1)
        fill!(w, 0.0)
        return nothing
    elseif x < 2h
        a =  rhoxb(x, h) - 1
    elseif x>1-2*h
        b = rhoxb(1-x, h) - 1
    end

    minus!(w, 1.0, xdata, n)
    w .= log.(w)
    wtmp = log.(xdata)
    multiply!(w, b)
    multiply!(wtmp, a)
    w .= w .+ wtmp

    # for ind in 1:n
    #     @inbounds w[ind] = a * log(xdata[ind]) + b * log(1 - xdata[ind])
    # end

    add!(w, -lbeta(a+1, b+1))
    w .= exp.(w)
    nothing
end
function betakernel(x::Real, logxdata::AbstractVector{<:Real}, log1_xdata::AbstractVector{<:Real}, h::Real, w::Vector, n::Int)
    a = x / h - 1
    b = (1 - x) / h - 1
    if (x < 0) || (x > 1)
        fill!(w, 0.0)
        return nothing
    elseif x < 2h
        a =  rhoxb(x, h) - 1
    elseif x>1-2*h
        b = rhoxb(1-x, h) - 1
    end

    for ind in 1:n
        @inbounds w[ind] = a * logxdata[ind] + b * log1_xdata[ind]
    end
    add!(w, -lbeta(a+1, b+1))
    w .= exp.(w)
    nothing
end
#f̂(x) = 1/n ∑ᵢ K(xᵢ;x /b+1, b )
#xdata should be positive, or domain error will be raised.
function gammakernel(x::Real, xdata::AbstractVector{<:Real}, h::Real, w::Vector, n::Int)
    rhob = x/h
    if x <= 0
        fill!(w, 0.0)
        return nothing
    elseif x < 2*h
        rhob = 0.25 * rhob * rhob + 1.0
    end

    w .= log.(xdata)
    multiply!(w, rhob-1.0)
    tmp = -rhob*log(h)-lgamma(rhob)
    add!(w, tmp)
    h1 = 1/h
    for ind in 1:n
        @inbounds w[ind] -= xdata[ind] * h1
    end
    w .= exp.(w)
    nothing
end
function gammakernel(x::Real, xdata::AbstractVector{<:Real}, logxdata::AbstractVector{<:Real}, h::Real, w::Vector, n::Int)
    rhob = x/h
    if x <= 0
        fill!(w, 0.0)
        return nothing
    elseif x < 2*h
        rhob = 0.25 * rhob * rhob + 1.0
    end

    # Yeppp.log!(w, xdata)
    multiply!(w, logxdata, rhob-1.0)
    tmp = -rhob*log(h)-lgamma(rhob)
    add!(w, tmp)
    h1 = 1/h
    for ind in 1:n
        @inbounds w[ind] -= xdata[ind] * h1
    end
    w .= exp.(w)
    nothing
end

function gaussiankernel(x::Real, xdata::AbstractVector{<:Real}, h::Real, w::Vector, n::Int)

    # for ind in 1:n
    #     @inbounds w[ind] = normlogpdf(xdata[ind], h, x)
    # end
    h1= 1.0/h
    # minus!(w, x, xdata, n)
    # multiply!(w, w, h1, n)
    # abs2!(w, w, n)
    # multiply!(w, w, -0.5, n)
    tmp = log(h) + log2π/2
    for ind in 1:n
        @inbounds w[ind]=-0.5*abs2((x - xdata[ind])*h1) - tmp
    end
    # add!(w, tmp, n)
    w .= exp.(w)

    nothing
end
function ekernel(x::Real, xdata::AbstractVector{<:Real}, h::Real, w::Vector, n::Int)
    ind = 1
    ind_end = 1+n
    @inbounds while ind < ind_end
        u = (x - xdata[ind]) / h
        w[ind] = ifelse(abs(u)>=1.0, 0.0, 1-u*u)
        ind += 1
    end
    multiply!(w, 0.75 / h)
    nothing
end

# #univariate normal kernel

# gaussiankernel{T<:FloatingPoint}(xi::T, x::T, h::T)=exp(-0.5*abs2((x - xi)/h)) / h * invsqrt2π
# gaussiankernel{T<:Real}(xi::T, x::T, h::T) = gkernel(float(xi), float(x),float(h))
# gaussiankernel(xi::Real, x::Real, h::Real)=gkernel(promote(xi, x, h)...)
# # type Gkernel <: Functor{3} end
# # NumericExtensions.evaluate(::Gkernel, xi, x, h) = gkernel(xi, x, h)

# function ekernel{T<:FloatingPoint}(xi::T, xeval::T, h::T)
#     u=(xeval - xi) / h
#     abs(u)>=1.0 ? 0.0 : 0.75*(1-u*u)
# end
# ekernel{T<:Real}(xi::T,xeval::T, h::T) = ekernel(float(xi),float(xeval), float(h))
# ekernel(xi::Real, xeval::Real, h::Real)=ekernel(promote(xi, xeval, h)...)
# # type Ekernel <: Functor{3} end
# # NumericExtensions.evaluate(::Ekernel, xi, xeval, h) = ekernel(xi, xeval, h)

# # # data need to be transformed on [0, ∞]
# # # K_{ρₕ(x), h}(xᵢ) = (t/h)^(ρₕ(x)-1) exp(-t/h)/gamma(ρ)/h
# function gammakernel{T<:FloatingPoint}(xi::T, x::T, b::T)
#     rhob = x/b
#     if x < 0
#         return(0)
#     elseif x < 2*b
#         rhob = 0.25 * rhob * rhob + 1
#     else
#         nothing
#     end
#     xi_b = xi/b
#     xi_b^(rhob-1) * exp(-xi_b) / b / gamma(rhob)
# end
# gammakernel{T<:Real}(xi::T,xeval::T, h::T) = gammakernel(float(xi),float(xeval), float(h))
# gammakernel(xi::Real, xeval::Real, h::Real)=gammakernel(promote(xi, xeval, h)...)
# # type Gammakernel <: Functor{3} end
# # NumericExtensions.evaluate(::Gammakernel, xi, xeval, h) = gammakernel(xi, xeval, h)



# # #data need be transformed on [0, 1]
# function betakernel{T<:FloatingPoint}(xi::T, x::T, h::T)
#     a = x / h
#     b = (1 - x) / h
#     if (x < 0) | (x > 1)
#         return(0)
#     elseif x < 2h
#         a =  rhoxb(x, h)
#     elseif x>1-2*h
#         b = rhoxb(1-x, h)
#     else
#         nothing
#     end
#     xi ^ (a-1) * (1 - xi) ^ (b-1) #/ beta(a, b)
# end
# betakernel{T<:Real}(xi::T,xeval::T, h::T) = betakernel(float(xi),float(xeval), float(h))
# betakernel(xi::Real, xeval::Real, h::Real)=betakernel(promote(xi, xeval, h)...)
# # type Betakernel <: Functor{3} end
# NumericExtensions.evaluate(::Betakernel, xi, xeval, h) = betakernel(xi, xeval, h)

# # #MultiVariate Normal Kernel
# # function GaussianKernel(xeval::Vector{Float64}, xi::Vector{Float64}, h::Vector{Float64})
# #   if any(h .<= 0.0)
# #     return Inf
# #   end
# #   p=length(xi)

# #   exp(-sumsq((xeval .- xi) ./ h)) / (2*pi)^(p/2) / prod(h)

# # end
# # # GaussianKernel2(xdiff::Float64, h::Float64) = GaussianKernel(xdiff, 0.0, sqrt(2)*h)
# # # GaussianKernel2(xdiff::Vector{Float64}, h::Vector{Float64}) = GaussianKernel(xdiff, zeros(length(xdiff)), sqrt(2).* h)

# # # Gaussian=KernelType(GaussianKernel, GaussianKernel2)

# # function EKernel(xeval::Float64, xi::Float64,h::Float64)
# #     if h <= 0.0
# #       return Inf
# #     end
# #     u=(xeval - xi) / h
# #     if abs(u) > 1.0
# #         return 0.0
# #     end
# #     return (1.0-u*u)*3.0/4.0
# # end


# # function EKernel(xeval::Vector{Float64}, xi::Vector{Float64}, h::Vector{Float64})

# #     if any(h .<= 0.0)
# #         return Inf
# #     end
# #     p=length(xi)
# #     tmp=1.0
# #     for i in 1:p
# #         tmp *= EKernel(xeval[i],xi[i],h[i])
# #     end
# #     tmp
# # end

# # function EKernel2(xdiff::Float64, h::Float64)
# #     ax=abs(xdiff / h)
# #     if ax > 2.0
# #         return 0.0
# #     end

# #     (2 - ax)^3*(ax^2+6*ax+4) * 3 / 160 / h

# # end
# # function EKernel2(xdiff::Vector{Float64}, h::Vector{Float64})
# #     ax=abs(xdiff ./ h)
# #     if any(ax .> 2.0)
# #         return 0.0
# #     end
# #     p=length(ax)
# #     tmp=1.0
# #     for i in 1:p
# #         tmp *= EKernel2(ax[i], h[i])
# #     end
# #     tmp
# # end
# # Epanechnikov=KernelType(EKernel, EKernel2)
