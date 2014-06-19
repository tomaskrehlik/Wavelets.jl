module Wavelets

	export 	modwt,
			mra,
			dwt

	abstract WaveletTransform

	include("waveletFilters.jl")
	include("extendSeries.jl")
	include("dwt.jl")
	include("modwt.jl")
	include("mra.jl")

end # module
