import csv
import numpy
import sys
import math
import datetime
import geohash
import bitstring

# longitude, latitude
world = ((-180.0, -90.0), (180.0, 90.0))

def load_data(filename = "viirs_full_data.csv"):
    rows = -1 # Ignore header row...
    with open(filename) as f:
        for row in f:
            rows += 1

    points = numpy.zeros((rows, 2))

    with open(filename) as f:
        for idx, row in enumerate(csv.DictReader(f)):
            points[idx][0] = float(row['longitude'])
            points[idx][1] = float(row['latitude'])

    return points


def spacepartition(min, max):
    size = (max[0] - min[0], max[1] - min[1])
    return (
        ((min[0],               min[1]),               (min[0] + size[0] / 2, min[1] + size[1] / 2)), # NW
        ((min[0] + size[0] / 2, min[1]),               (max[0],               min[1] + size[1] / 2)), # NE
        ((min[0] + size[0] / 2, min[1] + size[1] / 2), (max[0],               max[1])),     # SE
        ((min[0],               min[1] + size[1] / 2), (min[0] + size[0] / 2, max[1]))      # SW
        )

def create_quadtree(points):
    start_time = datetime.datetime.now()

    def split(start, end, min, max):
        replaceidx = start
        for idx in xrange(start, end):
            if (    points[idx][0] >= min[0]
                and points[idx][1] >= min[1]
                and points[idx][0] <= max[0]
                and points[idx][1] <= max[1]):
                if idx != replaceidx:
                    lon, lat = points[idx]
                    points[idx] = points[replaceidx]
                    points[replaceidx] = lon, lat
                replaceidx += 1
        return replaceidx

    def check_inside(start, end, min, max):
        for idx in xrange(start, end):
            if (   points[idx][0] < min[0]
                or points[idx][1] < min[1]
                or points[idx][0] > max[0]
                or points[idx][1] > max[1]):
                return False
        return True

    def check_different(start, end):
        for idx in xrange(start, end):
            if points[start][0] != points[idx][0] or points[start][1] != points[idx][1]:
                return True
        return False

    # Format: (start, end, child1, child2, child3, child4)
    # where the start/end are indexes into points and childx:es are an indexes into quadtree
    quadtree = numpy.ones((points.shape[0] * points.itemsize * 8 * 2, 6)) # *8 since 8 bits to a byte, *2 to account for 4 way branching...
    quadtree *= -1 # Just make sure we have "null pointers" everywhere to start with...
    quadtreepos = [0]

    def partition(start, end, min, max):
        res = quadtreepos[0]
        quadtreepos[0] += 1
        quadtree[res] = (start, end, -1, -1, -1, -1)
        if end - start > 1 and check_different(start, end):
            partitions = spacepartition(min, max)

            pstart = start
            for pidx in xrange(0, 4):
                if pidx < 3:
                    pend = split(pstart, end, *partitions[pidx])
                else:
                    pend = end
                if pend - pstart:
                    quadtree[res][2 + pidx] = partition(pstart, pend, *partitions[pidx])
                pstart = pend
        return res
    partition(0, points.shape[0], *world)

    end_time = datetime.datetime.now()
    print "Exec time:", end_time - start_time

    return quadtree

def print_quadtree(data, tree, maxlevel=5, pos=0, prefix='', min=world[0], max=world[1]):
    if tree[pos][2] != -1 or tree[pos][3] != -1 or tree[pos][4] != -1 or tree[pos][5] != -1:
        print prefix + "[%s,%s-%s,%s]: %s" % (min[0], min[1], max[0], max[1], int(tree[pos][1] - tree[pos][0]))
        partitions = spacepartition(min, max)
        for pidx in xrange(0, 4):
            if tree[pos][2 + pidx] != -1:
                print_quadtree(data, tree, maxlevel, tree[pos][2 + pidx], prefix + '  ', *partitions[pidx])
    else:
        for pidx in xrange(int(tree[pos][0]), int(tree[pos][1])):
            point = data[pidx]
            print prefix + "%s,%s" % (point[0], point[1])


def load_hashes(filename = "viirs_full_data.csv", precision = 12):
    with open(filename) as f:
        res = [
            geohash.encode(float(row['latitude']), float(row['longitude']), precision)
            for idx, row in enumerate(csv.DictReader(f))]
        res.sort()
        return res


def geohash_to_binary(hashcode):
    res = bitstring.BitArray()
    pos = 0
    for i in hashcode:
        res.append(bitstring.BitArray(uint=geohash._base32_map[i], length=5))
    return res

def binary_to_geohash(bincode):
    res = []
    for i in bincode.cut(5):
        res.append(geohash._base32[i.uint])
    return ''.join(res)

def binary_inc(str):
    str.reverse()
    for pos in xrange(0, len(str) + 1):
        if pos == len(str):
            str += '0b1'
            break
        if not str[pos]:
            str[pos] = 1
            break
        str[pos] = 0
    str.reverse()
    return str

def binary_search(key, lst):
    """Binary search a sorted list, returning the first item equal to or greater than key."""
    size = len(lst)
    step = size / 2 or 1
    pos = size / 2 or 1

    while True:
        step = (step / 2) or 1

        if step <= 0:
            raise Exception("Step <= 0")
        elif pos > size:
            raise Exception("pos > size")
        elif pos < 0:
            raise Exception("pos < 0")
        assert step > 0 and pos <= size and pos >= 0

        if pos == size:
            if lst[pos - 1] < key:
                return pos
            else:
                pos -= step
        elif lst[pos] < key:
            pos += step
        elif lst[pos] >= key:
            if pos < 1 or lst[pos - 1] < key:
                return pos
            else:
                pos -= step
        else:
            return pos

def find_range(prefix, lst):
    """Binary search a sorted list of geohashes
    Returns the range of indexes [lower, upper] that share the prefix
    """

    lower = binary_search(prefix, lst)
    binprefix = geohash_to_binary(prefix)
    if binprefix.all(1):
        upper = len(lst)
    else:
        upper = binary_search(
            binary_to_geohash(
                binary_inc(binprefix)),
            lst)

    return (lower, upper)

def minmaxdistsqr(hash1, hash2):
    """Calculates the square of the mindist and maxdist between two geohashes.
    Note: Treats lat/lon as cartesian."""
    hash1 = geohash.bbox(hash1)
    hash2 = geohash.bbox(hash2)
    distx = (abs(hash1['e'] - hash2['w']), abs(hash1['w'] - hash2['e']), abs(hash1['e'] - hash2['e']), abs(hash1['w'] - hash2['w']))
    disty = (abs(hash1['n'] - hash2['s']), abs(hash1['s'] - hash2['n']), abs(hash1['n'] - hash2['n']), abs(hash1['s'] - hash2['s']))
    print distx
    print disty

    mindist = min(*distx)**2 + min(*disty)**2
    maxdist = max(*distx)**2 + max(*disty)**2

    return mindist, maxdist

# geohash.bbox('0m1g52jnur27')
# {'s': -61.34645814076066, 'e': -166.1391907185316, 'w': -166.13919105380774, 'n': -61.3464579731226}


def build_locality(q, S):
    pass

