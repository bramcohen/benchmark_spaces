# benchmark_spaces
Benchmarks of how well different two dimensional spaces work for clustering algorithms

This may be useful for guiding anyone optimizing using clustering for data science or machine learning problems. It also could be used to make cool animations.

Benchmarking is done by putting different numbers of points in the space,
letting them repel each other for long enough to achieve equilibrium, and
measuring how close they get to everything being equidistant. The benchmark
number measures the standard deviation of the logs of distances. Lower is
better, zero is perfect.

To get metrics run benchmark.py from the command line.

There's some highly opinionated interpretation and commentary in Commentary.txt

Current output is:
```
Square
 3 0.03
 4 0.18
 5 0.27
 6 0.32
 7 0.34
 8 0.36
 9 0.36
10 0.39
11 0.40
12 0.41
13 0.42
14 0.43
Circle
 3 0.00
 4 0.20
 5 0.27
 6 0.29
 7 0.31
 8 0.33
 9 0.36
10 0.38
11 0.39
12 0.40
13 0.41
14 0.41
Torus
 3 0.00
 4 0.00
 5 0.00
 6 0.16
 7 0.21
 8 0.25
 9 0.21
10 0.26
11 0.29
12 0.30
13 0.29
14 0.30
Twisted Torus
 3 0.00
 4 0.00
 5 0.12
 6 0.00
 7 0.00
 8 0.18
 9 0.24
10 0.25
11 0.29
12 0.30
13 0.28
14 0.29
Klein Bottle
 3 0.00
 4 0.02
 5 0.00
 6 0.08
 7 0.14
 8 0.21
 9 0.21
10 0.26
11 0.28
12 0.29
13 0.29
14 0.30
Sphere Surface
 3 0.00
 4 0.00
 5 0.26
 6 0.29
 7 0.28
 8 0.30
 9 0.31
10 0.33
11 0.36
12 0.36
13 0.37
14 0.37
Projective Plane
 3 0.00
 4 0.00
 5 0.00
 6 0.00
 7 0.19
 8 0.25
 9 0.27
10 0.29
11 0.29
12 0.29
13 0.30
14 0.31
```
