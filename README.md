# Better spatial geometries for use in clustering algorithms
Contains implementations of alternative geometries for spaces which can be 
used in clustering algorithms. Using the spaces based on optimal sphere 
packing should result in slightly better affinity measurements with only 
minor extra computational requirements. Comparison is here:

![Space type comparisons](https://github.com/bramcohen/benchmark_spaces/blob/main/compare_everything.png?raw=true)

This metric is a good measure of how much unnecessary noise is added to 
distances and how easily points can slide past each other when being 
annealed. 

The lattice implementation benchmarks look a little bumpy because of 
strangeness in sphere packing in different numbers of dimensions. It is 
optimal in dimensions 2, 3, 4, 6, and 8, but there are improvements which 
could be implemented with a bunch of work in 5, 7, and 9 dimensions. Above 
that the best option is an open mathematical problem but those are large 
numbers of dimensions to use for clustering anyway.

For visualization 2d lattice should be used. It 'crystallizes' nicely as 
shown here:

![Crystallizing in 2d](https://github.com/bramcohen/benchmark_spaces/blob/main/movie_twisted_torus_6.gif?raw=true)

For applications which need antipodes (opposite points) Torus is best.
