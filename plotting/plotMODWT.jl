function plot(data::modwt)
	if data.series!=1
		error("Specify which series decomposition you want to plot.")
	else
		plot(data, 1)
	end
end

function plot(data::modwt, ser::Int)
	using Gadfly
	using DataFrames

	import Gadfly.plot
	
	df = convert(DataFrame, [reshape(data.W[:,:,ser], data.L, data.level) data.V[:,data.level,ser] data.original[:,ser] [1:data.L]])
	names!(df, [[symbol("W$i") for i in 1:data.level], symbol("V$levels"), symbol("Originalseries"), symbol("Time")])
	plotdata = stack(df, [1:(data.level+2)], [data.level+3])
	plot(plotdata, ygroup="variable", y="value", x="Time", Geom.subplot_grid(Geom.line))
end