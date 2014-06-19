# Wavelets
This package should eventually end up replicating the good parts of various implementations of wavelet methods in R and MATLAB that I have used.

## Done so far ##

Most of the things work fine, now should be tested.

- DWT works
- MODWT works
- MRA works
- Test suite against R wavelets package.


## To be done ##

Apart from adding some features, see below, some checks for when the user is not cautious engouh should be added. (For example, it is possible to feed arbitrary length series to DWT, which is obviously wrong.)

Features currently planned:

- Plotting techniques
- Better testing techniques
- Packet wavelet transform
- Coherences, phase differences estimates
 

Other features will be added based on issues. Plotting techniques are so far delayed due to packages not being mature enough. I will, however, add some basic plotting templates into examples using Gadfly.

## Inspired by: ##
So far, the code is insipired by [wavelets package from R](http://cran.r-project.org/web/packages/wavelets/index.html)

## License: ##
As the work is mainly derived from GPL packages, this package and code is licensed under the same terms, i.e. GPL. However, I requested a permission to promote the package to MIT license. If the permission will be granted I will do that instantly.
