function mra(wt::WaveletTransform)
  D = fill(0.0, wt.L, wt.level, wt.series)
  S = fill(0.0, wt.L, wt.level, wt.series)

  # compute the details and smooths
  for j = 1:wt.level
    DWj = wt.W
    SVj = wt.V
    SWj = fill(0.0, wt.L, wt.level, wt.series)
    DVj = fill(0.0, wt.L, wt.level, wt.series)

    for k=j:-1:1
      for i=1:wt.series
        DVj[:,j,i] = eval(Expr(:call, symbol(string(typeof(wt),"Backward")), DWj[:,j,i], DVj[:,j,i], wt.filter, k))
        SVj[:,j,i] = eval(Expr(:call, symbol(string(typeof(wt),"Backward")), SWj[:,j,i], SVj[:,j,i], wt.filter, k))
      end
    
      DWj = fill(0.0, wt.L, wt.level, wt.series)
      SWj = DWj
    end

    D[:,j,:] = DVj[:,j,:]
    S[:,j,:] = SVj[:,j,:]
  end
  
  return (D, S)
end
