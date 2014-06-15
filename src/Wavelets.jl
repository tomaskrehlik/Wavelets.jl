module Wavelets

	export 	modwt,
			mra,
			dwt

	include("waveletFilters.jl")
	include("extendSeries.jl")
	include("dwt.jl")
	include("modwt.jl")
	include("mra.jl")

end # module
