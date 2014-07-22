module Wavelets

	export 	modwt,
			mra,
			dwt,
			inverse,
			dwtDo,
			modwtDo,
			modwtNew,
			haarFilter

	abstract WaveletTransform

	include("waveletFilters.jl")
	include("extendSeries.jl")
	include("dwt.jl")
	include("modwt.jl")
	include("mra.jl")

	# Defining global inverse function for any wavelet transform

	function inverse(wt::WaveletTransform)
		return eval(Expr(:call, symbol(string("i", typeof(wt),"Do")), wt))
	end

end # module
