# Read in the data and load packages

using Wavelets
using DataFrames
using Gadfly

# set_default_plot_size(25cm,25cm)

path = Pkg.dir("Wavelets")
data = readcsv(joinpath(path,"tests","data","testdata.csv"))

# Let us first plot the time series and only concentrate on the first one
(N,k) = size(data)
plot(x=[1:N], y=data[:,1], Geom.line)

# Now do the wavelet decomposition using MODWT on four levels using LA8 wavelet
levels = 4
decomposition = modwt(data[:,1], "la8", 4, "periodic")

# Look what do we have inside the dcompostiion object, should be quite strightforward to guess what means what
names(decomposition)

# Now plot the various wavelet coefficients with the last scaling coefficient and the original time series

# Still dont know how to reorder the things in the picture, if anyone has a clue, please help!!!
df = convert(DataFrame, [reshape(decomposition.W[:,:,1], decomposition.L, decomposition.level) decomposition.V[:,decomposition.level,1] decomposition.original[:,1] [1:decomposition.L]])
names!(df, [[symbol("W$i") for i in 1:decomposition.level], symbol("V$levels"), symbol("Originalseries"), symbol("Time")])
plotdata = stack(df, [1:(decomposition.level+2)], [decomposition.level+3])
plot(plotdata, ygroup="variable", y="value", x="Time", Geom.subplot_grid(Geom.line))

# You can do the similar thing with MRA just use the MRA method
(D, S) = mra(decomposition)
df = convert(DataFrame, [reshape(D[:,:,1], decomposition.L, decomposition.level) S[:,decomposition.level,1] decomposition.original[:,1] [1:decomposition.L]])
names!(df, [[symbol("D$i") for i in 1:decomposition.level], symbol("S$levels"), symbol("Originalseries"), symbol("Time")])
plotdata = stack(df, [1:(decomposition.level+2)], [decomposition.level+3])
plot(plotdata, ygroup="variable", y="value", x="Time", Geom.subplot_grid(Geom.line))