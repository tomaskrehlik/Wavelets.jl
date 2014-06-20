function extendSeries(X::Array{Float64}, method::String, length::String, n::Int)
  
  (N, ) = size(X)
    
  # determine final length 'n' of series after extension
  if (length == "arbitrary")
    if (n <= N)
     error("Invalid argument: 'n' must be greater than length of series when length='arbitrary'.")
    end
  elseif (length == "powerof2")
    k = N/(2^n)
    if (round(k) == k)
      error("Invalid argument: length of series should not be multiple of 2^j when length='powerof2'.")
    else
      n = ceil(k)*2^n
    end
  elseif (length == "double") 
    error("The integer parameter is arbitrary, you probably made some error.")
  end

  if isa(X, Vector)
    # extend the series to length 'n'
    (method == "periodic") ? X = repeat(X, outer=[convert(Int,ceil(n/N))])[1:n] : ()
    (method == "reflection") ? X = repeat([X,X[end:-1:1]], outer=[convert(Int,ceil(n/N))])[1:n] : ()
    (method == "zeros") ? X = [X, zeros(N-n)] : ()
    (method == "mean") ? X = [X, mean(X).*ones(N-n)] : ()
    # if(method == "reflection.inverse") X = apply(X, 2, function(x,n,N){rep(c(x,2*x[N]-rev(x)),length=n)}, n=n, N=N)
  else
    # extend the series to length 'n'
    (method == "periodic") ? X = mapslices(x -> repeat(x, outer=[convert(Int,ceil(n/N))])[1:n], X, 1) : ()
    (method == "reflection") ? X = mapslices(x -> repeat(x, outer=[convert(Int,ceil(n/N))])[1:n], [X,X[end:-1:1]], 2) : ()
    (method == "zeros") ? X = mapslices(x -> [x, zeros(N-n)], X, 2) : ()
    (method == "mean") ? X = mapslices(x -> [x, mean(x).*ones(N-n)], X, 2) : ()
    # if(method == "reflection.inverse") X = apply(X, 2, function(x,n,N){rep(c(x,2*x[N]-rev(x)),length=n)}, n=n, N=N)
  end

  return X
end

function extendSeries(X::Array{Float64}, method::String, length::String)
  (N, ) = size(X)

  if length == "double"
    n = 2*N
  else
    error("The selected method does not correspond, you probably forgot to include another parameter.")
  end

  if isa(X, Vector)
    # extend the series to length 'n'
    (method == "periodic") ? X = repeat(X, outer=[convert(Int,ceil(n/N))])[1:n] : ()
    (method == "reflection") ? X = repeat([X,X[end:-1:1]], outer=[convert(Int,ceil(n/N))])[1:n] : ()
    (method == "zeros") ? X = [X, zeros(N-n)] : ()
    (method == "mean") ? X = [X, mean(X).*ones(N-n)] : ()
    # if(method == "reflection.inverse") X = apply(X, 2, function(x,n,N){rep(c(x,2*x[N]-rev(x)),length=n)}, n=n, N=N)
  else
    # extend the series to length 'n'
    (method == "periodic") ? X = mapslices(x -> repeat(x, outer=[convert(Int,ceil(n/N))])[1:n], X, 1) : ()
    (method == "reflection") ? X = mapslices(x -> repeat(x, outer=[convert(Int,ceil(n/N))])[1:n], [X,X[end:-1:1]], 2) : ()
    (method == "zeros") ? X = mapslices(x -> [x, zeros(N-n)], X, 2) : ()
    (method == "mean") ? X = mapslices(x -> [x, mean(x).*ones(N-n)], X, 2) : ()
    # if(method == "reflection.inverse") X = apply(X, 2, function(x,n,N){rep(c(x,2*x[N]-rev(x)),length=n)}, n=n, N=N)
  end

  return X
end