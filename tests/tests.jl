include("/Users/tomaskrehlik/Documents/Julia/wavelets/waveletFilters.jl")
include("/Users/tomaskrehlik/Documents/Julia/wavelets/dwt.jl")
include("/Users/tomaskrehlik/Documents/Julia/wavelets/modwt.jl")
include("/Users/tomaskrehlik/Documents/Julia/wavelets/mra.jl")

data = readcsv("data/testdata.csv")

# MODWT tests

boundaries = ["periodic"]
filters = ["haar", map(x->string("la",x),[8:2:20]), map(x->string("d",x),[4:2:20]), map(x->string("c",x),[6:6:30])]
levels = [1:4]

for i in boundaries, j in filters, l in levels
	run(`Rscript modwt.R $j $l $i`)
	W = readcsv("W.csv")
	V = readcsv("V.csv")
	decomp = modwt(data[:,1], j, l, i)
	(Wjl, Vjl) = (decomp.W, decomp.V)
	@assert all(abs(Wjl - W) .< 2.220446049250313e-12)
	@assert all(abs(Vjl - V) .< 2.220446049250313e-12)
end

run(`rm W.csv`)
run(`rm V.csv`)

# DWT tests

boundaries = ["periodic"]
filters = ["haar", map(x->string("la",x),[8:2:20]), map(x->string("d",x),[4:2:20]), map(x->string("c",x),[6:6:30])]
levels = [1:4]

for i in boundaries, j in filters, l in levels
	run(`Rscript dwt.R $j $l $i`)
	W = readcsv("W.csv")
	V = readcsv("V.csv")
	decomp = dwt(data[1:512,1], j, l, i)
	(Wjl, Vjl) = (decomp.W, decomp.V)
	@assert all((abs(Wjl - W) .< 2.220446049250313e-12) $ (isnan(abs(Wjl - W))))
	@assert all((abs(Vjl - V) .< 2.220446049250313e-12) $ (isnan(abs(Vjl - V))))
end


#  MRA tests for modwt

boundaries = ["periodic"]
filters = ["haar", map(x->string("la",x),[8:2:20]), map(x->string("d",x),[4:2:20]), map(x->string("c",x),[6:6:30])]
levels = [1:4]

for i in boundaries, j in filters, l in levels
	run(`Rscript mra.R $j $l $i modwt`)
	D = readcsv("D.csv")
	S = readcsv("S.csv")
	decomp = = modwt(data[:,1], j, l, i)
	(Djl, Sjl) = mra(decomp)
	@assert all(abs(Djl - D) .< 2.220446049250313e-12)
	@assert all(abs(Sjl - S) .< 2.220446049250313e-12)
end

run(`rm D.csv`)
run(`rm S.csv`)