function alignCoef(coef, wavelet, coe, modwt, inverse, filter)
	(N, J, ser) = size(coef)
    for j=1:J
      shift = wtFilterShift(filter, j, wavelet, coe, modwt)

      shift >= N ? shift = shift - floor(shift/N)*N : shift = shift
      
      if shift == 0
        coef[:,j,:] = coef[:,j,:]
      end
      
      if inverse
      	coef[:,j,:] = coef[[(end-shift+1):end,1:(end-shift)],j,:]
      else
      	coef[:,j,:] = coef[[(shift+1):end,1:shift],j,:]
      end
    end
    return coef
end

function align(wt::WaveletTransform, coe::Boolean, inverse::Boolean)
  # Do the alignment only if neccessary
  @assume !inverse == wt.aligned
  
  # Check what kind of transform do I have.
  wt.filter.transform == "modwt" ? (modwt = true) : (modwt = false)
  
  # Make a copy of the original object
  wtShifted = deepcopy(wt)

  # Adjust the coefficients
  wtShifted.W = alignCoef(wt.W, true, coe, modwt, inverse, wt.filter)
  wtShifted.V = alignCoef(wt.V, false, coe, modwt, inverse, wt.filter)
  
  inverse ? wtShifted.aligned = false : wtShifted.aligned = true

  return wtShifted
end