immutable modwt <: WaveletTransform
  
  original::Array{Float64}
  L::Int
  series::Int
  level::Int
  filter::waveletFilter
  boundary::ASCIIString
  W::Array{Float64, 3}
  V::Array{Float64, 3}
  aligned::Bool

  function modwt(X::Array{Float64}, filter::ASCIIString, nLevels::Int, boundary::ASCIIString)
    filter = eval(Expr(:call, symbol(string(filter,"Filter")), 1, true))

    if isa(X, Vector)
      (N,) = size(X)
      nSeries = 1
    else
      (N, nSeries) = size(X)
    end

    (W, V) = modwtDo(X, filter, nLevels, boundary, nSeries, N)
    
    new(X, N, nSeries, nLevels, filter, boundary, W, V, false)

  end

end

# Old version which is about 2 times slower
# function modwtForward(V::Vector{Float64}, filter::waveletFilter, j::Int)
#   (N,) = size(V)
#   Wj = fill(NaN, N)
#   Vj = fill(NaN, N)
  
#   for t = 0:(N-1) 
#     k = t
#     Wjt = filter.h[1]*V[k+1]
#     Vjt = filter.g[1]*V[k+1]
#     for n = 1:(filter.L-1) 
#       k = k - 2^(j-1)
#       if (k < 0)
#         k = k + ceil(-k/N)*N
#       end
#       Wjt = Wjt + filter.h[n+1]*V[k+1]
#       Vjt = Vjt + filter.g[n+1]*V[k+1]
#     end
#     Wj[t+1] = Wjt
#     Vj[t+1] = Vjt    
#   end
  
#   return (Wj, Vj)
# end

function modwtForward(V::Vector{Float64}, filter::waveletFilter, j::Int)
  (N,) = size(V)

  indeces = fill(0, N, filter.L)
  for i=1:N
    indeces[i, 1] = i
    for k=2:filter.L
      indeces[i, k] = indeces[i, k-1] - 2^(j-1)
      if indeces[i, k] <= 0
        indeces[i, k] = N + indeces[i, k]
      end
    end
  end

  return (filter.h' * V[indeces]', filter.g' * V[indeces]')
end

function modwtNew(X::Array{Float64}, filter::waveletFilter, nLevels::Int, boundary::ASCIIString, nSeries::Int, N::Int)
  # reflect X for reflection method
  if (boundary == "reflection")
    X = extendSeries(X, boundary, "double")
    N = 2*N
  end

  # initialize variables for pyramid algorithm
  nBoundary = fill(0.0, nLevels)
  WCoefs = fill(0.0, N, nLevels, nSeries)
  VCoefs = fill(0.0, N, nLevels, nSeries)

  for j=1:nLevels
    Lj = (2^j-1)*(filter.L-1)+1
    nBoundary[j] = min(Lj,N)

    indeces = fill(0, N, filter.L)
    for l=1:N
      indeces[l, 1] = l
      for k=2:filter.L
        indeces[l, k] = indeces[l, k-1] - 2^(j-1)
        if indeces[l, k] <= 0
          indeces[l, k] = N + indeces[l, k]
        end
      end
    end
    for i=1:nSeries
      if j==1 
        WCoefs[:,j,i], VCoefs[:,j,i] = filter.h' * X[:,i][indeces]', filter.g' * X[:,i][indeces]'
      else
        WCoefs[:,j,i], VCoefs[:,j,i] = filter.h' * VCoefs[:,j-1,i][indeces]', filter.g' * VCoefs[:,j-1,i][indeces]'
      end
    end
  end
  return (WCoefs, VCoefs)
end

function modwtBackward(W::Vector{Float64}, V::Vector{Float64}, filter::waveletFilter, j::Int)
  (N,) = size(V)
  (L,) = size(filter.h)
  Vj = fill(NaN, N)

  for t = 0:(N-1)
    k = t
    Vjt = filter.h[1]*W[k+1] + filter.g[1]*V[k+1]
    for n = 1:(L-1)
      k = k + 2^(j-1)
      if k >= N
        k = k - floor(k/N)*N
      end
      Vjt = Vjt + filter.h[n+1]*W[k+1] + filter.g[n+1]*V[k+1]
    end
    Vj[t+1] = Vjt
  end

  return Vj
end

function modwtDo(X::Array{Float64}, filter::waveletFilter, nLevels::Int, boundary::ASCIIString, nSeries::Int, N::Int)
  # reflect X for reflection method
  if (boundary == "reflection")
    X = extendSeries(X, boundary, "double")
    N = 2*N
  end

  # initialize variables for pyramid algorithm
  nBoundary = fill(0.0, nLevels)
  WCoefs = fill(0.0, N, nLevels, nSeries)
  VCoefs = fill(0.0, N, nLevels, nSeries)
  
  # implement the pyramid algorithm
  for i=1:nSeries
    Vj = X[:,i]
    for j=1:nLevels
      (WCoefs[:,j,i], VCoefs[:,j,i]) = modwtForward(Vj,filter,j)
      Vj = convert(Vector, VCoefs[:,j,i])
      Lj = (2^j-1)*(filter.L-1)+1
      nBoundary[j] = min(Lj,N)
    end
  end
  
  return (WCoefs, VCoefs)
end

function imodwtDo(wt::modwt)
  output = fill(0.0, wt.L, wt.series)

  for i=1:wt.series
    Vj = convert(Vector, wt.V[:,end,i])
    for j=wt.level:-1:1
      Vj = modwtBackward(convert(Vector, wt.W[:,j,i]), Vj, wt.filter, j)
      Vj = convert(Vector, Vj)
    end
    output[:,i] = Vj
  end

  return output

end