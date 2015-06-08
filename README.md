# rnetcarto

[![Build Status](https://travis-ci.org/geeklhem/rnetcarto.svg?branch=master)](https://travis-ci.org/geeklhem/rnetcarto) [![codecov.io](http://codecov.io/github/geeklhem/rnetcarto/coverage.svg?branch=master)](http://codecov.io/github/geeklhem/rnetcarto?branch=master) 

Fast network modularity and roles computation by simulated annealing  ([rgraph C library](https://github.com/geeklhem/rgraph) wrapper for R).

### Example
``` R
> library(rnetcarto)

> # A general network.
> library(igraphdata)
> data(karate)
> web = get.adjacency(karate,sparse=FALSE)
> netcarto(web)
[[1]]
       name module     z-score participation
1     Mr Hi      2  2.13570766  5.390625e-01
2   Actor 2      2  1.40155815  1.975309e-01
3   Actor 3      2  0.30033389  6.200000e-01
4   Actor 4      2  0.66740864  1.110223e-16
5   Actor 5      1 -0.81649658  4.444444e-01
(...) 

[[2]]
[1] 0.4197896


> # A bipartite network
> library(bipartite)
> data(vazarr)
> bipartmod(vazarr)
[[1]]
                       name module    z-score participation
1     Aristotelia chilensis      1  0.0000000     0.4444444
2        Alstroemeria aurea      5  0.0000000     0.6805556
3       Schinus patagonicus      1  0.0000000     0.6111111
4         Berberis darwinii      2  1.4142136     0.5000000
5           Vicia nigricans      5  0.0000000     0.5000000
6           Rosa eglanteria      2 -0.7071068     0.4489796
7          Cynanchum diemii      3  0.0000000     0.0000000
8        Ribes magellanicum      2 -0.7071068     0.0000000
9         Mutisia decurrens      4  0.0000000     0.0000000
10 Calceolaria crenatiflora      6  0.0000000     0.0000000

[[2]]
[1] 0.2357462
```
