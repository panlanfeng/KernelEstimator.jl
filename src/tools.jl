#expand.grid(x, y)
function ExpandGrid{T}(x::Vector{T}, y::Vector{T})
  n1=length(x)
  n2=length(y)
  grid=Array(eltype(x), n1*n2,2)
  for i =1:n1
    for j=1:n2
      grid[n2 * (i - 1)+j,1] = x[i]
      grid[n2 * (i - 1)+j,2] = y[j]
    end
  end
  grid
end
