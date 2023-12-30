from random import seed, uniform, gauss, shuffle, randrange
from statistics import stdev, mean
from math import log, acos

from PIL import Image, ImageDraw
import imageio
from colorsys import hsv_to_rgb
import matplotlib.pyplot as plt

def make_colors(num):
    r = [hsv_to_rgb(i/num, 1, 1) for i in range(num)]
    return [tuple(int(256*y) for y in x) for x in r]

def draw_image_torus(points, circlesize, boxsize = 100, numboxes = 3):
    img = Image.new('RGB', (numboxes*boxsize, numboxes*boxsize), tuple([192] * 3))
    window = ImageDraw.Draw(img)
    cs = circlesize/2
    colors = make_colors(len(points))
    for i in range(-1, numboxes+1):
        for j in range(-1, numboxes+1):
            for color, point in zip(colors, points):
                cx = point.x * boxsize + i * boxsize
                cy = point.y * boxsize + j * boxsize
                window.ellipse([(cx-cs, cy-cs), (cx+cs, cy+cs)], fill=color)
    return img

def draw_image_twisted_torus(points, circlesize = 20, boxsize = 100, numboxes = 3):
    img = Image.new('RGB', (numboxes*boxsize, numboxes*boxsize), tuple([192] * 3))
    window = ImageDraw.Draw(img)
    cs = circlesize/2
    colors = make_colors(len(points))
    for i in range(-1, numboxes * 2):
        for j in range(-1, numboxes * 2):
            for color, point in zip(colors, points):
                py, px = point.get_coords()
                cy = px * boxsize + (3**0.5)/2 * i * boxsize
                cx = py * boxsize + j * boxsize + boxsize * (i % 2)/2
                window.ellipse([(cx-cs, cy-cs), (cx+cs, cy+cs)], fill=color)
    return img

def magnitude(vals):
    return sum(x*x for x in vals) ** 0.5

class point_in_cube:
    """
    A point bounded by a hypercube. When it gets out of bounds it gets pushed back inside.

    People intuitively assume that when packing into cubes the different 
    dimensions will correspond to real phenomena which provide insights into the 
    data but in reality hypercube is a spiky ball, with many things getting stuck in 
    random corners.
    """
    def __init__(self, dimensions):
        self.coords = [uniform(0, 1) for d in range(dimensions)]

    def distance(self, other):
        return magnitude([x-y for x, y in zip(self.coords, other.coords)])

    def move_from(self, other, dist):
        deltas = [x-y for x, y in zip(self.coords, other.coords)]
        m = dist/magnitude(deltas)
        for d, c in enumerate(self.coords):
            c += deltas[d] * m
            if c < 0:
                c = -c/2
            if c > 1:
                c = 1-(c-1)/2
            self.coords[d] = c

class multipoint:
    """
    Used to implement a higher dimensional point using lower dimensional points. Never optimal 
    but can come close and performs well and is easy to implement.

    Treating the internal points as a black box causes unnecessary square rooting and re-squaring 
    which should be elided in any port to a performance critical system.
    """
    def __init__(self, points):
        self.points = points

    def clone(self):
        return multipoint([x.clone() for x in self.points])
    
    def _canonicalize(self):
        for p in self.points:
            p._canonicalize()

    def distance(self, other):
        return magnitude([a.distance(b) for a, b in zip(self.points, other.points)])

    def move_from(self, other, distance):
        ivals = [a._closest(b) for a, b in zip(self.points, other.points)]
        mag = distance / magnitude([x for x, y in ivals])
        for p, (d, c) in zip(self.points, ivals):
            p.coords = [tocoord + (tocoord-fromcoord)*mag for (tocoord, fromcoord) in zip(p.coords, c)]
            p._canonicalize()

class point_on_torus:
    """
    A point in a Torus represented by a hypercube. Wraps around from one side to the other.
    """
    def __init__(self, dimensions = 2, other = None):
        if other:
            self.coords = [x for x in other.coords]
        else:
            self.coords = [uniform(0, 1) for i in range(dimensions)]

    def clone(self):
        return point_on_torus(0, self)

    def distance(self, other):
        return self._distance(other)[0]

    def _canonicalize(self):
        self.coords = [x % 1 for x in self.coords]

    def get_coords(self):
        return self.coords

    def distance(self, other):
        return self._closest(other)[0]

    def move_from(self, other, distance):
        d, c = self._closest(other)
        mag = distance/d
        self.coords = [tocoord + (tocoord-fromcoord)*mag for (tocoord, fromcoord) in zip(self.coords, c)]
        self._canonicalize()

    def _closest(self, other):
        r = []
        for me, oth in zip(self.coords, other.coords):
            if oth > me:
                if oth - me > 0.5:
                    r.append(oth - 1)
                else:
                    r.append(oth)
            else:
                if me-oth > 0.5:
                    r.append(oth + 1)
                else:
                    r.append(oth)
        return (sum((a-b)**2 for (a, b) in zip(self.coords, r))**(1/2), r)

last_dimensions_to_8 = [-1, -1, -1, -1, -1, 4**(1/2), 3**(1/2), 2**(1/2), 1**(1/2)]

class point_to_8d:
    """
    An efficient implementation of space which wraps around according to an optimal sphere 
    packing in 6 and 8 dimensions. 

    These are constructed by starting with a torus, removing one color from checkerboard coloring 
    it, and adding in antipodal points by offsetting everything by half. In 6 dimensions the 
    last dimension is stretched to make the major and minor diagonals the same length.
    
    Doesn't work in odd numbers because going up to the antipode twice hits an eliminated point.

    If you want to put these in proximity buckets you need to do one final step of 
    canonicalization so every point has one representation instead of 2.
    """
    def __init__(self, d, other_coords = None, other_parity = None):
        if other_coords:
            self.coords = [x for x in other_coords]
            self.parity = other_parity
        else:
            self.coords = [uniform(0, 1) for x in range(d-1)]
            self.coords.append(uniform(0, last_dimensions_to_8[d]))
            self.parity = randrange(2)

    def clone(self):
        return point_to_8d(0, self.coords, self.parity)

    def _canonicalize(self):
        for i in range(len(self.coords) - 1):
            while self.coords[i] > 1:
                self.coords[i] -= 1
                self.parity ^= 1
            while self.coords[i] < 0:
                self.coords[i] += 1
                self.parity ^= 1
        lastlen = last_dimensions_to_8[len(self.coords)]
        while self.coords[-1] > lastlen:
            self.coords[-1] -= lastlen
            self.parity ^= 1
        while self.coords[-1] < 0:
            self.coords[-1] += lastlen
            self.parity ^= 1

    def distance(self, other):
        return self._closest(other)[0]

    def move_from(self, other, distance):
        d, c = self._closest(other)
        mag = distance/d
        self.coords = [tocoord + (tocoord-fromcoord)*mag for (tocoord, fromcoord) in zip(self.coords, c)]
        self._canonicalize()

    def _alternate(self):
        other = self.clone()
        for i in range(len(other.coords) - 1):
            other.coords[i] += 0.5
        other.coords[-1] += last_dimensions_to_8[len(other.coords)]/2
        other._canonicalize()
        return other

    def _closest(self, other):
        d1, coords1 = self._closest2(other)
        d2, coords2 = self._closest2(other._alternate())
        if d1 < d2:
            return (d1, coords1)
        else:
            return (d2, coords2)

    def _closest2(self, other):
        parity = self.parity ^ other.parity
        d = len(self.coords)
        lowest_dim = -1
        lowest_amount = 2
        r = []
        for i in range(d - 1):
            oc = other.coords[i]
            mc = self.coords[i]
            if oc > mc:
                diff = 2 * (oc - mc) - 1
                if diff > 0:
                    r.append(oc - 1)
                    parity ^= 1
                    if diff < lowest_amount:
                        lowest_dim = i
                        lowest_amount = diff
                else:
                    r.append(oc)
                    if -diff < lowest_amount:
                        lowest_dim = i
                        lowest_amount = -diff
            else:
                diff = 2 * (mc - oc) - 1
                if diff > 0:
                    r.append(oc + 1)
                    parity ^= 1
                    if diff < lowest_amount:
                        lowest_dim = i
                        lowest_amount = diff
                else:
                    r.append(oc)
                    if -diff < lowest_amount:
                        lowest_dim = i
                        lowest_amount = -diff
        oc = other.coords[-1]
        mc = self.coords[-1]
        lastval = last_dimensions_to_8[d]
        if oc > mc:
            diff = 2 * (oc - mc) - lastval
            if diff > 0:
                r.append(oc - lastval)
                parity ^= 1
                if lastval * diff < lowest_amount:
                    lowest_dim = d-1
            else:
                r.append(oc)
                if -lastval * diff < lowest_amount:
                    lowest_dim = d-1
        else:
            diff = 2 * (mc - oc) - lastval
            if diff > 0:
                r.append(oc + lastval)
                parity ^= 1
                if lastval * diff < lowest_amount:
                    lowest_dim = d-1
            else:
                r.append(oc)
                if -lastval * diff < lowest_amount:
                    lowest_dim = d-1
        if parity:
            if lowest_dim == d-1:
                if r[-1] > self.coords[-1]:
                    r[-1] -= lastval
                else:
                    r[-1] += lastval
            else:
                if r[lowest_dim] > self.coords[lowest_dim]:
                    r[lowest_dim] -= 1
                else:
                    r[lowest_dim] += 1
        return (sum((a-b)**2 for (a, b) in zip(self.coords, r))**(1/2), r)

last_dimensions = [-1, -1, 2*(3/4)**(1/2), 2*(2/4)**(1/2), 2*(1/4)**(1/2)]

class point_to_4d:
    """
    An efficient implementation of a space which wraps around corresponding to optimal 
    sphere packing in 2, 3, or 4 dimensions.
    """
    def __init__(self, d, original = None):
        if original:
            self.coords = [x for x in original.coords]
        else:
            self.coords = [uniform(0, 1) for x in range(d-1)]
            self.coords.append(uniform(0, last_dimensions[d]))

    def clone(self):
        return point_to_4d(0, self)

    def _canonicalize(self):
        r = [x % 1 for x in self.coords[:-1]]
        r.append(self.coords[-1] % last_dimensions[len(self.coords)])
        self.coords = r

    def get_coords(self):
        lastval = last_dimensions[len(self.coords)]
        if self.coords[-1] > lastval / 2:
            r = [(x + 1/2) % 1 for x in self.coords[:-1]]
            r.append(self.coords[-1] - lastval / 2)
            return r
        return self.coords

    def distance(self, other):
        return self._closest(other)[0]

    def move_from(self, other, distance):
        d, c = self._closest(other)
        mag = distance/d
        self.coords = [tocoord + (tocoord-fromcoord)*mag for (tocoord, fromcoord) in zip(self.coords, c)]
        self._canonicalize()

    def _closest(self, other):
        d1, coords1 = self._closest2(other.coords)
        alternate = [(x + 0.5) % 1 for x in other.coords[:-1]]
        lval = last_dimensions[len(self.coords)]
        alternate.append((other.coords[-1] + lval / 2) % lval)
        d2, coords2 = self._closest2(alternate)
        if d1 < d2:
            return (d1, coords1)
        else:
            return (d2, coords2)

    def _closest2(self, other):
        r = []
        for me, oth in zip(self.coords[:-1], other[:-1]):
            if oth > me:
                if oth - me > 0.5:
                    r.append(oth - 1)
                else:
                    r.append(oth)
            else:
                if me-oth > 0.5:
                    r.append(oth + 1)
                else:
                    r.append(oth)
        lval = last_dimensions[len(self.coords)]
        me, oth = self.coords[-1], other[-1]
        if oth > me:
            if oth - me > lval/2:
                r.append(oth - lval)
            else:
                r.append(oth)
        else:
            if me - oth > lval/2:
                r.append(oth + lval)
            else:
                r.append(oth)
        return (sum((a-b)**2 for (a, b) in zip(self.coords, r))**(1/2), r)

class point_on_sphere:
    """
    A point on the surface of a sphere
    """
    def __init__(self, dimensions):
        coords = [gauss(0, 1) for d in range(dimensions+1)]
        m = 1/magnitude(coords)
        self.coords = [x*m for x in coords]

    def distance(self, other):
        return acos(max(sum([x*y for x, y in zip(self.coords, other.coords)]), -1))

    def move_from(self, other, dist):
        deltas = [x-y for x, y in zip(self.coords, other.coords)]
        m = dist/magnitude(deltas)
        for d in range(len(self.coords)):
            self.coords[d] += deltas[d] * m
        m = magnitude(self.coords)
        self.coords = [x/m for x in self.coords]

class point_on_projective_sphere:
    """
    A point on the surface of a Projective Sphere, implemented by gluing
    together opposite points on the surface of a Sphere.
    """
    def __init__(self, dimensions):
        coords = [gauss(0, 1) for d in range(dimensions+1)]
        m = 1/magnitude(coords)
        self.coords = [x*m for x in coords]

    def distance(self, other):
        return acos(abs(sum([x*y for x, y in zip(self.coords, other.coords)])))

    def move_from(self, other, dist):
        ocoords = other.coords
        if sum([x*y for x, y in zip(self.coords, ocoords)]) < 0:
            ocoords = [-x for x in ocoords]
        deltas = [x - y for x, y in zip(self.coords, ocoords)]
        m = dist/magnitude(deltas)
        for d in range(len(self.coords)):
            self.coords[d] += deltas[d] * m
        m = magnitude(self.coords)
        self.coords = [x/m for x in self.coords]

def show_anneal(points, draw):
    images = []
    def ni(points, generation):
        if generation % 10 != 0:
            return
        images.append(draw(points))
    anneal(points, callback = ni)
    imageio.mimsave('movie.gif', images)

def anneal(points, iterations = 500, callback = None):
    """
    Repeatedly move each point away from the one closest to it following an
    exponentially decreasing annealing schedule
    """
    d = 0.1
    for i in range(iterations, 0, -1):
        if callback:
            callback(points, i)
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

# An example use of the show_anneal function to make an image illustrating annealing
#show_anneal([point_to_4d(2) for i in range(6)], draw_image_twisted_torus)

def benchmark_factory(factory):
    """
    For the given factory print out how close to equidistant it gets different
    numbers of points when they're all trying to repel each other. Runs each one
    repeatedly and averages to compensate for different local maxima.

    This can show artifacts of a space not captured by the other benchmark.
    """
    for k in range(3, 40):
        results = []
        for q in range(100):
            points = [factory() for i in range(k)]
            anneal(points)
            vals = []
            for i in range(len(points)):
                for j in range(i):
                    vals.append(log(points[i].distance(points[j])))
            results.append(stdev(vals))
        print('%2i %1.2f %1.2f' % (k, mean(results), stdev(results)))

lattice_factories = [lambda: point_to_4d(2), 
        lambda: point_to_4d(3), 
        lambda: point_to_4d(4), 
        lambda: multipoint([point_to_4d(4), point_on_torus(1)]),
        lambda: point_to_8d(6),
        lambda: multipoint([point_to_8d(6), point_on_torus(1)]),
        lambda: point_to_8d(8),
        lambda: multipoint([point_to_8d(8), point_on_torus(1)]),
        lambda: multipoint([point_to_8d(8), point_to_4d(2)]),
        lambda: multipoint([point_to_8d(8), point_to_4d(3)]),
        lambda: multipoint([point_to_8d(8), point_to_4d(4)]),
        lambda: multipoint([point_to_8d(8), point_to_4d(4), point_on_torus(1)]),
        lambda: multipoint([point_to_8d(8), point_to_8d(6)]),
        lambda: multipoint([point_to_8d(8), point_to_8d(6), point_on_torus(1)]),
        lambda: multipoint([point_to_8d(8), point_to_8d(8)])]

#benchmark_factory(lattice_factories[8 - 2])

def find_average_distance(factory):
    distances = []
    for i in range(1000000):
        distances.append(log(factory().distance(factory())))
    return stdev(distances)

def show_distances():
    dimensions = range(2, 17)
    cubevals = [find_average_distance(lambda: point_in_cube(dim)) for dim in dimensions]
    torusvals = [find_average_distance(lambda: point_on_torus(dim)) for dim in dimensions]
    spherevals = [find_average_distance(lambda: point_on_sphere(dim)) for dim in dimensions]
    pspherevals = [find_average_distance(lambda: point_on_projective_sphere(dim)) for dim in dimensions]
    packevals = [find_average_distance(lp) for lp in lattice_factories]
    fig, ax = plt.subplots()
    ax.set_xticks(dimensions)
    ax.set_ylim(ymin=0, ymax=0.7)
    ax.set_title('Standard deviation of log of distance between random points\n(lower is better)')
    ax.plot(dimensions, spherevals, 'o-', color='blue', label='Sphere')
    ax.plot(dimensions, pspherevals, 'o-', color='red', label='Projective Sphere')
    ax.plot(dimensions, cubevals, 'o-', color='pink', label='Cube')
    ax.plot(dimensions, torusvals, 'o-', color='cyan', label='Torus')
    ax.plot(dimensions, packevals, 'o-', color='green', label='Lattice')
    ax.set_xlabel('Dimensions')
    ax.set_ylabel('deviation')
    ax.legend()
    plt.show()

#show_distances()

def close(a, b):
    assert a-b < .00001
    assert b-a < .00001

# For debugging, doesn't work on sphere-based geometries
def test_points(factory):
    for i in range(10000):
        p1 = factory()
        p2 = factory()
        close(p1.distance(p2), p2.distance(p1))
        p3 = p2.clone()
        p2.move_from(p1, 0.1)
        close(0.1, p2.distance(p3))

#test_points(make_multipoint_3_2)
