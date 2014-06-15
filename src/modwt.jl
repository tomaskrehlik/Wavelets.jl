type modwt
  
  L::Int
  series::Int
  level::Int
  filter::waveletFilter
  boundary::ASCIIString
  W::Array{Float64, 3}
  V::Array{Float64, 3}

  function modwt(X::Array{Float64}, filter::ASCIIString, nLevels::Int, boundary::ASCIIString)
    filter = eval(Expr(:call, symbol(string(filter,"Filter")), 1, true))
    if isa(X, Vector)
      (N,) = size(X)
      nSeries = 1
    else
      (N, nSeries) = size(X)
    end
    (W, V) = modwtDo(X, filter, nLevels, boundary, nSeries, N)
    new(L, series, nLevels, filter, boundary, W, V)
  end

end

function modwtForward(V::Vector{Float64}, filter::waveletFilter, j::Int)
  (N,) = size(V)
  Wj = fill(NaN, N)
  Vj = fill(NaN, N)
  
  for t = 0:(N-1) 
    k = t
    Wjt = filter.h[1]*V[k+1]
    Vjt = filter.g[1]*V[k+1]
    for n = 1:(filter.L-1) 
      k = k - 2^(j-1)
      if (k < 0)
        k = k + ceil(-k/N)*N
      end
      Wjt = Wjt + filter.h[n+1]*V[k+1]
      Vjt = Vjt + filter.g[n+1]*V[k+1]
    end
    Wj[t+1] = Wjt
    Vj[t+1] = Vjt    
  end
  
  return (Wj, Vj)
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
    X = extendSeries(X, boundary)
    N = 2*N
  end

  # initialize variables for pyramid algorithm
  nBoundary = fill(0, nLevels)
  WCoefs = fill(0, N, nLevels, nSeries)
  VCoefs = fill(0, N, nLevels, nSeries)
  

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