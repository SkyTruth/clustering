import csv
import numpy
import sys
import math
import datetime

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

    def spacepartition(min, max):
        size = (max[0] - min[0], max[1] - min[1])
        return (
            ((min[0],               min[1]),               (min[0] + size[0] / 2, min[1] + size[1] / 2)), # NW
            ((min[0] + size[0] / 2, min[1]),               (max[0],               min[1] + size[1] / 2)), # NE
            ((min[0] + size[0] / 2, min[1] + size[1] / 2), (max[0],               max[1])),     # SE
            ((min[0],               min[1] + size[1] / 2), (min[0] + size[0] / 2, max[1]))      # SW
            )

    # Format: (start, end, child1, child2, child3, child4)
    # where the start/end are indexes into points and childx:es are an indexes into quadtree
    quadtree = numpy.ones((points.shape[0] * points.itemsize * 8 * 2, 6)) # *8 since 8 bits to a byte, *2 to account for 4 way branching...
    quadtree *= -1 # Just make sure we have "null pointers" everywhere to start with...
    quadtreepos = [0]

    def partition(start, end, min, max, prefix = ""):
        # print prefix + "PART", start, end, min, max
        res = quadtreepos[0]
        quadtreepos[0] += 1
        quadtree[res] = (start, end, -1, -1, -1, -1)
        if end - start > 1:
            partitions = spacepartition(min, max)

            pstart = start
            for pidx in xrange(0, 4):
                if pidx < 3:
                    pend = split(pstart, end, *partitions[pidx])
                else:
                    pend = end
                different = check_different(pstart, pend)
                # if not check_inside(pstart, pend, *partitions[pidx]):
                #     print prefix + "  NOT INSIDE", pstart, pend, partitions[pidx]
                #     raise Exception
                if pend - pstart and different:
                    quadtree[res][2 + pidx] = partition(pstart, pend, prefix=prefix+" ", *partitions[pidx])
                # if pend == pstart:
                #     print prefix + "  NULL ", pstart, pend, partitions[pidx]
                # elif not different:
                #     print prefix + "  SAME ", pstart, pend, partitions[pidx]
                pstart = pend
        return res
    partition(0, points.shape[0], *world)

    end_time = datetime.datetime.now()
    print "Exec time:", end_time - start_time

    return quadtree
    
