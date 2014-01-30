import csv
import numpy
import sys
import math
import datetime
import geohash
import bitstring

# longitude, latitude
world = ((-180.0, -90.0), (180.0, 90.0))

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

def count(prefix, lst):
    lower, upper = find_range(prefix, lst)
    return max(upper - lower, 0)

def maxdistsqr(hash1, hash2):
    """Calculates the square of the maxdist between two geohashes.
    Note: Treats lat/lon as cartesian."""
    hash1 = geohash.bbox(hash1)
    hash2 = geohash.bbox(hash2)

    distx = max(hash1['e'], hash2['e']) - min(hash1['w'], hash2['w'])
    disty = max(hash1['n'], hash2['n']) - min(hash1['s'], hash2['s'])
    return distx**2 + disty**2

def mindistsqr(hash1, hash2):
    """Calculates the square of the mindist between two geohashes.
    Note: Treats lat/lon as cartesian."""
    hash1 = geohash.bbox(hash1)
    hash2 = geohash.bbox(hash2)
    if hash1['w'] < hash2['w']:
        distx = max(hash2['w'] - hash1['e'], 0) # Return 0 if there's an overlap of the boxes.
    else:
        distx = max(hash1['w'] - hash2['e'], 0) # Return 0 if there's an overlap of the boxes.

    if hash1['s'] < hash2['s']:
        disty = max(hash2['s'] - hash1['n'], 0) # Return 0 if there's an overlap of the boxes.
    else:
        disty = max(hash1['s'] - hash2['n'], 0) # Return 0 if there's an overlap of the boxes.
    
    return distx**2 + disty**2

# geohash.bbox('0m1g52jnur27')
# {'s': -61.34645814076066, 'e': -166.1391907185316, 'w': -166.13919105380774, 'n': -61.3464579731226}



def build_locality(query_block, Q, k):
    locality = set()
    pruned_list = set()
    prune_dist = 0
    total = k
    # Calculate mindists and maxdists for all the blocks
    Q = [(hash, mindistsqr(query_block, hash), maxdistsqr(query_block, hash)) for hash in Q]
    # Sort by decreasing maxdist order, so we can use pop() to get the
    # one with lowest maxdist
    Q.sort(lambda a, b: cmp(b[2], a[2]))
    while total >= 0:
        b, mindist, maxdist = Q.pop()
        prune_dist = maxdist
        total -= count(b)
        locality.add(b)
    # Sort by decreasing mindist order, so we can use pop() to get the
    # one with lowest mindist
    Q.sort(lambda a, b: cmp(b[1], a[1]))
    while Q:
        b, mindist, maxdist  = Q.pop()
        if mindist <= prune_dist:
            locality.add(b)
        else:
            pruned_list.add(b)
    return locality, pruned_list, math.sqrt(prune_dist)
