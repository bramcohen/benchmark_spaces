from random import seed, uniform, gauss, shuffle
from statistics import stdev, mean
from math import log, acos

from PIL import Image, ImageDraw
import imageio
from colorsys import hsv_to_rgb
import matplotlib.pyplot as plt

def make_colors(num):
    r = [hsv_to_rgb(i/num, 1, 1) for i in range(num)]
    return [tuple(int(256*y) for y in x) for x in r]

def draw_image_disc(points, circlesize = 20, size = 100):
    img = Image.new('RGB', (numboxes*boxsize, numboxes*boxsize), tuple([192] * 3))
    window = ImageDraw.Draw(img)
    cs = circlesize/2
    colors = make_colors(len(points))
    for color, point in zip(colors, points):
        cx = point.coords[0] * size * 0.4 + size/2
        cy = point.coords[1] * size * 0.4 + size/2
        window.ellipse([(cx-cs, cy-cs), (cx+cs, cy+cs)], fill=color)
    return img

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
                cy = point.x * boxsize + isoconst * i * boxsize
                cx = point.y * boxsize + j * boxsize + boxsize * (i % 2)/2
                window.ellipse([(cx-cs, cy-cs), (cx+cs, cy+cs)], fill=color)
    return img

def magnitude(vals):
    return sum(x*x for x in vals) ** 0.5

class point_in_square:
    """
    A point bounded by a square. When it gets out of the square it gets pushed
    back to the edge.
    """
    def __init__(self, dimensions = 2):
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

class point_on_disc:
    """
    A point bounded by a disc. When it gets out of the disc it gets pushed
    back to the edge.
    """
    def __init__(self, dimensions = 2):
        coords = [gauss(0, 1) for d in range(dimensions)]
        m = (1 - uniform(0, 1) ** dimensions)/magnitude(coords)
        self.coords = [x*m for x in coords]

    def distance(self, other):
        return magnitude([x-y for x, y in zip(self.coords, other.coords)])

    def move_from(self, other, dist):
        deltas = [x-y for x, y in zip(self.coords, other.coords)]
        m = dist/magnitude(deltas)
        for d in range(len(self.coords)):
            self.coords[d] += deltas[d] * m
        m = magnitude(self.coords)
        if m > 1:
            mult = (1-(m-1)/2)/m
            self.coords = [x*mult for x in self.coords]

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
            ds.append((magnitude([dx, dy]), dx, dy))
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
            ds.append((magnitude([dx, dy]), dx, dy))
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
        x, y = other.x, other.y
        cs = [(x, y)]
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
            ds.append((magnitude([dx, dy]), dx, dy))
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
    def __init__(self, dimensions = 3):
        coords = [gauss(0, 1) for d in range(dimensions)]
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
    def __init__(self, dimensions = 3):
        coords = [gauss(0, 1) for d in range(dimensions)]
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

# An example use of the show_anneal function to make an image illustrating annealing
#show_anneal([point_on_twisted_torus() for i in range(6)], draw_image_twisted_torus)

def anneal(points, iterations = 500, callback = None):
    """
    Repeatedly move each point away from the one closest to it following an
    exponentially decreasing annealing schedule
    """
    d = 0.1
    for i in range(iterations, 0, -1):
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

def find_average_distance(factory):
    refpoint = factory()
    distances = []
    for i in range(1000000):
        distances.append(log(refpoint.distance(factory())))
    return stdev(distances)

def show_distances():
    dimensions = range(3, 21)
    spherevals = [find_average_distance(lambda: point_on_sphere(dim)) for dim in dimensions]
    pspherevals = [find_average_distance(lambda: point_on_projective_sphere(dim)) for dim in dimensions]
    fig, ax = plt.subplots()
    ax.set_xticks(dimensions)
    ax.set_ylim(ymin=0)
    ax.set_title('Standard deviation of log of distance between random points')
    ax.plot(dimensions, spherevals, 'o-', color='blue', label='Sphere')
    ax.plot(dimensions, pspherevals, 'o-', color='red', label='Projective Sphere')
    ax.set_xlabel('Dimensions')
    ax.set_ylabel('deviation')
    ax.legend()
    plt.show()

show_distances()
