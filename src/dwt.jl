immutable dwt <: WaveletTransform
  
  original::Array{Float64}
  L::Int
  series::Int
  level::Int
  filter::waveletFilter
  boundary::ASCIIString
  W::Array{Float64, 3}
  V::Array{Float64, 3}

  function dwt(X::Array{Float64}, filter::ASCIIString, nLevels::Int, boundary::ASCIIString)
    filter = eval(Expr(:call, symbol(string(filter,"Filter")), 1, false))

    if isa(X, Vector)
      (N,) = size(X)
      nSeries = 1
    else
      (N, nSeries) = size(X)
    end

    (W, V) = dwtDo(X, filter, nLevels, boundary, nSeries, N)
    
    new(X, N, nSeries, nLevels, filter, boundary, W, V)

  end

end


function dwtBackwards(W::Vector{Float64}, V::Vector{Float64}, filter::waveletFilter)
  (M, ) = size(V)
  (L, ) = size(filter.h)
  Vj = fill(NaN, 2*M)
  l = -2
  m = -1
  for t = 0:(M-1)
    l = l+2
    m = m+2
    u = t
    i = 1
    k = 0

    Vj[l+1] = filter.h[i+1]*W[u+1] + filter.g[i+1]*V[u+1]
    Vj[m+1] = filter.h[k+1]*W[u+1] + filter.g[k+1]*V[u+1]

    if (L > 2)
      for n = 1:((L/2)-1)
        u = u+1
        if (u >= M)
        	u = 0
        end
        i = i+2
        k = k+2
        Vj[l+1] = Vj[l+1] + filter.h[i+1]*W[u+1] + filter.g[i+1]*V[u+1]
        Vj[m+1] = Vj[m+1] + filter.h[k+1]*W[u+1] + filter.g[k+1]*V[u+1]
      end
    end
  end

  return Vj
end

# Defined because of the MRA

function dwtBackwards(W::Vector{Float64}, V::Vector{Float64}, filter::waveletFilter, j::Int)
	return dwtBackwards(W, V, filter)
end

function dwtForward(V::Vector{Float64}, filter::waveletFilter)
  (M,) = size(V)
  Wj = fill(NaN, convert(Int64,M/2))
  Vj = fill(NaN, convert(Int64,M/2))
  for (t = 0:(M/2 - 1))
    u = 2*t + 1
    Wjt = filter.h[1]*V[u+1]
    Vjt = filter.g[1]*V[u+1]
    for (n = 1:(filter.L-1))
      u = u - 1
      if (u < 0)
       		u = M - 1
   	  end
      Wjt = Wjt + filter.h[n+1]*V[u+1]
      Vjt = Vjt + filter.g[n+1]*V[u+1]
    end
    Wj[t+1] = Wjt
    Vj[t+1] = Vjt
  end
  return (Wj, Vj)
end

function idwtDo(wt::dwt)
  output = fill(0.0, wt.L, wt.series)

  for i=1:wt.series
    Vj = convert(Vector, wt.V[:,end,i])
    for j=wt.level:-1:1
      Vj = modwtBackward(Vj, wt.filter)
      Vj = convert(Vector, Vj)
    end
    output[:,i] = Vj
  end

  return output
  
end

function dwtDo(X::Array{Float64}, filter::ASCIIString, nLevels::Int, boundary::ASCIIString, nSeries::Int, N::Int)
  # reflect X for reflection method
  if (boundary == "reflection")
    X = extendSeries(X, boundary)
    N = 2*N
  end

  # initialize variables for pyramid algorithm
  nBoundary = fill(0.0, nLevels)
  WCoefs = fill(NaN, convert(Int64,N/2), nLevels, nSeries)
  VCoefs = fill(NaN, convert(Int64,N/2), nLevels, nSeries)

  # implement the pyramid algorithm
  for i=1:nSeries
    Vj = X[:,i]
    for j=1:nLevels
      (WCoefs[1:(N/(2^j)),j,i], VCoefs[1:(N/(2^j)),j,i]) = dwtForward(Vj,filter)
      Vj = convert(Vector, VCoefs[1:(N/(2^j)),j,i])
      Lj = ceil((filter.L-2)*(1-(1/(2^j))))
      Nj = N/(2^j)
      nBoundary[j] = min(Lj,Nj)
    end
  end

  return (WCoefs, VCoefs)
end