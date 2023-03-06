# benchmark_spaces
Implements different two dimensional spaces for clustering algorithms

Contains code for making animations of simulated annealing of points 
repelling each other and also for measuring the average deviation in 
distances between random points (lower is better).

The upshort of this is:

If you're doing two dimensional visualtizations you should use the 
geometry called 'twisted torus'. 

Otherwise if having exact opposites isn't part of expected/wanted 
semantics of your data you should use Projective Sphere.

In the unusual case of polar opposites being expected you should 
use Sphere.

More detailed thoughts are in Commentary.txt
