from random import seed, uniform, shuffle
from statistics import stdev, mean
from math import log, acos

def magnitude(*vals):
    return sum(x*x for x in vals) ** 0.5

class point_in_square:
    """
    A point bounded by a square. When it gets out of the square it gets pushed
    back to the edge.
    """
    def __init__(self):
        self.x, self.y = uniform(0, 1), uniform(0, 1)

    def distance(self, other):
        return magnitude(self.x-other.x, self.y-other.y)

    def move_from(self, other, d):
        dx, dy = self.x - other.x, self.y - other.y
        m = d/magnitude(dx, dy)
        self.x += dx*m
        self.y += dy*m
        if self.x < 0:
            self.x = 0
        if self.y < 0:
            self.y = 0
        if self.x > 1:
            self.x = 1
        if self.y > 1:
            self.y = 1

class point_in_circle:
    """
    A point bounded by a circle. When it gets out of the circle it gets pushed
    back to the edge.
    """
    def __init__(self):
        while True:
            self.x, self.y = uniform(-1, 1), uniform(-1, 1)
            if magnitude(self.x, self.y) < 1:
                break

    def distance(self, other):
        return magnitude(self.x-other.x, self.y-other.y)

    def move_from(self, other, d):
        dx, dy = self.x - other.x, self.y - other.y
        m = d/magnitude(dx, dy)
        self.x += dx*m
        self.y += dy*m
        m2 =  magnitude(self.x, self.y)
        if m2 > 1:
            self.x /= m2
            self.y /= m2

class point_on_torus:
    """
    A point in a Torus represented by a planar square. When it gets outside of
    the square it wraps around to the other side.
    """
    def __init__(self):
        self.x, self.y = uniform(0, 1), uniform(0, 1)

    def distance(self, other):
        return self._distance(other)[0]

    def _distance(self, other):
        cs = []
        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                cs.append((other.x + i, other.y + j))
        ds = []
        for x, y in cs:
            dx = self.x - x
            dy = self.y - y
            ds.append((magnitude(dx, dy), dx, dy))
        return min(ds)

    def move_from(self, other, d):
        mag, dx, dy = self._distance(other)
        m = d/mag
        self.x += dx*m
        self.y += dy*m
        self.x %= 1
        self.y %= 1

class point_on_klein_bottle:
    """
    A point on a Klein Bottle represented by a planar square. When it gets
    outside of the square it wraps around to the other side. In one of the
    dimensions wrapping also flips it upside down in the other dimension.
    """
    def __init__(self):
        self.x, self.y = uniform(0, 1), uniform(0, 1)

    def distance(self, other):
        return self._distance(other)[0]

    def _distance(self, other):
        cs = []
        for j in (-1, 0, 1):
            cs.append((other.x, other.y + j))
        for i in (-1, 1):
            for j in (-1, 0, 1):
                cs.append((other.x + i, 1 - other.y + j))
        ds = []
        for x, y in cs:
            dx = self.x - x
            dy = self.y - y
            ds.append((magnitude(dx, dy), dx, dy))
        return min(ds)

    def move_from(self, other, d):
        mag, dx, dy = self._distance(other)
        m = d/mag
        self.x += dx*m
        self.y += dy*m
        self.y %= 1
        if self.x < 0 or self.x > 1:
            self.y = 1 - self.y
            self.x %= 1

isoconst = (3**0.5)/2

class point_on_twisted_torus:
    """
    A point in a Torus with a half turn in one of the two dimensions.
    The aspect ratio (3 ** 0.5)/2 is used to allow seven points to be exactly
    equidistant.
    """
    def __init__(self):
        self.x, self.y = uniform(0, isoconst), uniform(0, 1)

    def distance(self, other):
        return self._distance(other)[0]

    def _distance(self, other):
        cs = []
        x, y = other.x, other.y
        for j in (-1, 0, 1):
            cs.append((x, y + j))
        cs.append((x + isoconst, y + 0.5))
        cs.append((x + isoconst, y - 0.5))
        cs.append((x - isoconst, y + 0.5))
        cs.append((x - isoconst, y - 0.5))
        ds = []
        for x, y in cs:
            dx = self.x - x
            dy = self.y - y
            ds.append((magnitude(dx, dy), dx, dy))
        return min(ds)

    def move_from(self, other, d):
        mag, dx, dy = self._distance(other)
        m = d/mag
        self.x += dx*m
        self.y += dy*m
        if self.x > isoconst or self.x < 0:
            self.y += 0.5
        self.y %= 1
        self.x %= isoconst

class point_on_sphere:
    """
    A point on the surface of a sphere
    """
    def __init__(self):
        while True:
            self.x, self.y, self.z = uniform(-1, 1), uniform(-1, 1), uniform(-1, 1)
            if magnitude(self.x, self.y, self.z) < 1:
                break

    def distance(self, other):
        return acos(self.x*other.x+self.y*other.y+self.z*other.z)

    def move_from(self, other, d):
        dx, dy, dz = self.x-other.x, self.y-other.y, self.z-other.z
        m = d/magnitude(dx, dy, dz)
        self.x += dx*m
        self.y += dy*m
        self.z += dz*m
        m = magnitude(self.x, self.y, self.z)
        self.x /= m
        self.y /= m
        self.z /= m

class point_on_projective_plane:
    """
    A point on the surface of a Projective Plane, implemented by gluing
    together opposite points on the surface of a Sphere.
    """
    def __init__(self):
        while True:
            self.x, self.y, self.z = uniform(-1, 1), uniform(-1, 1), uniform(-1, 1)
            if magnitude(self.x, self.y, self.z) < 1:
                break

    def distance(self, other):
        return acos(abs(self.x*other.x+self.y*other.y+self.z*other.z))

    def move_from(self, other, d):
        x, y, z = other.x, other.y, other.z
        if self.x*x+self.y*y+self.z*z < 0:
            x, y, z = -x, -y, -z
        dx, dy, dz = self.x-x, self.y-y, self.z-z
        m = d/magnitude(dx, dy, dz)
        self.x += dx*m
        self.y += dy*m
        self.z += dz*m
        m = magnitude(self.x, self.y, self.z)
        self.x /= m
        self.y /= m
        self.z /= m

def anneal(points, iterations = 500):
    """
    Repeatedly move each point away from the one closest to it following an
    exponentially decreasing annealing schedule
    """
    d = 0.1
    for i in range(iterations, 0, -1):
        ps = [x for x in points]
        shuffle(ps)
        for p in ps:
            m = 10
            mo = None
            for o in points:
                if o is not p:
                    mn = p.distance(o)
                    if mn < m:
                        m = mn
                        mo = o
            p.move_from(mo, d)
        d *= 0.99

def benchmark_factory(factory):
    """
    For the given factory print out how close to equidistant it gets different
    numbers of points when they're all trying to repel each other. Run each one
    repeatedly and average to compensate for different local maxima.
    """
    for k in range(3, 15):
        results = []
        for q in range(100):
            points = [factory() for i in range(k)]
            anneal(points)
            vals = []
            for i in range(len(points)):
                for j in range(i):
                    vals.append(log(points[i].distance(points[j])))
            results.append(stdev(vals))
        print('%2i %1.2f' % (k, mean(results)))

from sys import argv

def main():
    if len(argv) > 1:
        seed(int(argv[1]))
    print('Square')
    benchmark_factory(point_in_square)
    print('Circle')
    benchmark_factory(point_in_circle)
    print('Torus')
    benchmark_factory(point_on_torus)
    print('Twisted Torus')
    benchmark_factory(point_on_twisted_torus)
    print('Klein Bottle')
    benchmark_factory(point_on_klein_bottle)
    print('Sphere Surface')
    benchmark_factory(point_on_sphere)
    print('Projective Plane')
    benchmark_factory(point_on_projective_plane)

main()
