immutable dwt <: WaveletTransform
  
  original::Array{Float64}
  L::Int
  series::Int
  level::Int
  filter::waveletFilter
  boundary::ASCIIString
  W::Array{Float64, 3}
  V::Array{Float64, 3}
  aligned::Bool

  function dwt(X::Array{Float64}, filter::ASCIIString, nLevels::Int, boundary::ASCIIString)
    if isa(X, Vector)
      (N,) = size(X)
      nSeries = 1
    else
      (N, nSeries) = size(X)
    end

    #@assert log(N)%log(2)<eps() "The length of the series should be 2^J."

    filter = eval(Expr(:call, symbol(string(filter,"Filter")), 1, false))

    (W, V) = dwtDo(X, filter, nLevels, boundary, nSeries, N)
    
    new(X, N, nSeries, nLevels, filter, boundary, W, V, false)

  end

end


function dwtBackward(W::Vector{Float64}, V::Vector{Float64}, filter::waveletFilter)
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

function dwtBackward(W::Vector{Float64}, V::Vector{Float64}, filter::waveletFilter, j::Int)
	return dwtBackward(W, V, filter)
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

function dwtForward!(orig::Array{Float64}, W::Array{Float64}, V::Array{Float64}, filter::waveletFilter, levels::Int, series::Int)
  @inbounds begin
    for i=1:series
      for l=1:levels
        if l==1
          M = length(orig)
          for (t = 1:(M/2))
            u = 2*(t-1) + 1
            W[t, l, i] = filter.h[1]*orig[u+1,i]
            V[t, l, i] = filter.g[1]*orig[u+1,i]
            for (n = 2:filter.L)
              u <= 0 ? (u = M - 1) : (u -= 1)
              W[t, l, i] += filter.h[n]*orig[u+1,i]
              V[t, l, i] += filter.g[n]*orig[u+1,i]
            end
          end  
        else
          M = convert(Int, length(orig)/(2^l))
          for (t = 1:(M/2))
            u = 2*(t-1) + 1
            W[t, l, i] = filter.h[1]*V[u+1,l-1,i]
            V[t, l, i] = filter.g[1]*V[u+1,l-1,i]
            for (n = 2:filter.L)
              u <= 0 ? (u = M - 1) : (u -= 1)
              W[t, l, i] += filter.h[n]*V[u+1,l-1,i]
              V[t, l, i] += filter.g[n]*V[u+1,l-1,i]
            end
          end  
        end
      end
    end
  end
  return nothing
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

function dwtDo(X::Array{Float64}, filter::waveletFilter, nLevels::Int, boundary::ASCIIString, nSeries::Int, N::Int)
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
  dwtForward!(X, WCoefs, VCoefs, filter, nLevels, nSeries)
  for i=1:nSeries
    for j=1:nLevels
      nBoundary[j] = min(ceil((filter.L-2)*(1-(1/(2^j)))), N/(2^j))
    end
  end

  return (WCoefs, VCoefs)
end