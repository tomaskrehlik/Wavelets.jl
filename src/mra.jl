function mra(X::Array{Float64}, filter::ASCIIString, nLevels::Int, boundary::ASCIIString, method::ASCIIString)
  # convert X to a matrix
  if isa(X, Vector)
    (N,) = size(X)
    nSeries = 1
  else
    (N, nSeries) = size(X)
  end
  # initializations
  
  D = fill(0.0, N, nLevels, nSeries)
  S = fill(0.0, N, nLevels, nSeries)

  # compute wavelet transform and other necessary values
  (wtW, wtV) = eval(Expr(:call, symbol(method), X, filter, 1, boundary))
  filter = eval(Expr(:call, symbol(string(filter,"Filter")), 1, method=="modwt" ? true : false))

  (N, nLevels, nSeries) = size(wtW)

  # compute the details and smooths
  for j = 1:nLevels
    DWj = wtW
    SVj = wtV
    SWj = fill(0.0, N, nLevels, nSeries)
    DVj = fill(0.0, N, nLevels, nSeries)

    for k=j:-1:1
      for i=1:nSeries
        DVj[:,j,i] = eval(Expr(:call, symbol(string(method,"Backward")), DWj[:,j,i], DVj[:,j,i], filter, k))
        SVj[:,j,i] = eval(Expr(:call, symbol(string(method,"Backward")), SWj[:,j,i], SVj[:,j,i], filter, k))
      end
    
      DWj = fill(0.0, N, nLevels, nSeries)
      SWj = DWj
    end

    D[:,j,:] = DVj[:,j,:]
    S[:,j,:] = SVj[:,j,:]
  end
  
  return (D, S)
end
