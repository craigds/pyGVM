"""
Simplified port of GVM for python.

GVM homepage: http://www.tomgibara.com/clustering/fast-spatial/java-library
"""
import bisect
import collections
from pygvm.libs import heapq
import sys
from array import array

VERSION = (0, 2, 5)

try:
    MAX_FLOAT = sys.float_info.max
except AttributeError:
    # python <2.6 didn't have float_info
    MAX_FLOAT = 9999999999999999999999999.0


class Cluster(object):
    def __init__(self, clusters):
        self.removed = False
        self.mass = 1.0

        # mass-weighted coordinate sum
        self.m1 = [0] * clusters.dimension
        # mass-weighted coordinate-square sum
        self.m2 = [0] * clusters.dimension

        self.variance = 0.0
        self.pairs = []
        self.set_members([])

    @property
    def center(self):
        if self.mass:
            return [coord / self.mass for coord in self.m1]
        else:
            return self.m1[:]

    def __len__(self):
        return len(self.members)

    def clear(self):
        """
        Completely clears this cluster. All points and their associated mass is removed
        """
        self.m0 = 0
        self.m1 = [0] * len(self.m1)
        self.m2 = [0] * len(self.m2)
        self.variance = 0.0
        self.set_members([])

    def set(self, mass, coords, key):
        """
        Sets this cluster equal to a single point.
        """
        self.m1 = [mass * coord for coord in coords]
        self.m2 = [mass * coord * coord for coord in coords]
        self.mass = mass
        self.variance = 0.0
        self.set_members([])
        if key is not None:
            self.add_members([key])

    def add(self, mass, coords, key):
        """
        Adds a point to the cluster.
        """
        self.add_bulk([(mass, coords, key)])

    def add_bulk(self, items):
        """
        Adds a number of points to the cluster.
        This is significantly faster than calling add() a bunch of times.
        Expects an iterable of tuples, each of which would form the arguments to add().
        i.e.: [(mass, coords, key), ...]
        """
        keys = []
        for mass, coords, key in items:
            if not self.mass:
                self.set(mass, coords, key)
            else:
                if mass != 0:
                    self.mass += mass
                    for i, coord in enumerate(coords):
                        self.m1[i] += coord * mass
                        self.m2[i] += coord * coord * mass
                    if key is not None:
                        keys.append(key)
        if keys:
            self.add_members(keys)
        self.update()

    def set_members(self, keys, already_sorted=False):
        self.members = set(keys)

    def add_members(self, keys, already_sorted=False):
        self.members.update(keys)

    def remove_member(self, key):
        self.members.remove(key)

    def remove_members(self, keys):
        self.members.difference_update(keys)

    def __contains__(self, key):
        return key in self.members

    def add_cluster(self, cluster):
        """
        Adds the specified cluster to this cluster.
        """
        self.add(cluster.mass, cluster.center, None)
        self.add_members(cluster.members, already_sorted=True)

    def remove(self, mass, coords, key):
        """
        Removes a point from the cluster.
        This is not strictly needed for GVM, as members are never removed,
        but it may be useful for creating derivative algorithms.
        """
        self.mass -= mass
        for i, coord in enumerate(coords):
            self.m1[i] -= coord * mass
            self.m2[i] -= coord * mass * mass
        self.remove_member(key)
        self.update()

    def test(self, mass, coords):
        """
        Computes the change in this clusters variance if it were to have a new point added to it.
        """
        new_mass = self.mass + mass
        if new_mass == 0:
            variance = 0
        else:
            total = 0.0
            for (coord, m1, m2) in zip(coords, self.m1, self.m2):
                m1 += (mass * coord)
                m2 += (mass * coord * coord)
                total += (m2 * new_mass) - (m1 * m1)
            variance = Clusters.correct(total / new_mass)
        return variance - self.variance

    def test_cluster(self, cluster):
        """
        Computes the variance of a cluster that aggregated this cluster with the supplied cluster
        """
        new_mass = self.mass + cluster.mass
        if new_mass == 0:
            return 0.0
        else:
            total = 0
            for (tm1, tm2, cm1, cm2) in zip(self.m1, self.m2, cluster.m1, cluster.m2):
                m1 = tm1 + cm1
                m2 = tm2 + cm2
                total += (m2 * new_mass) - (m1 * m1)
            return Clusters.correct(total / new_mass)

    def update(self):
        """
        Recompute this cluster's variance.
        """
        if self.mass == 0:
            self.variance = 0
        else:
            total = 0.0
            for m1, m2 in zip(self.m1, self.m2):
                total += ((m2 * self.mass) - (m1 * m1))
            self.variance = Clusters.correct(total / self.mass)

    def __repr__(self):
        return u'<Cluster: %d nodes>' % len(self)


class SortedArrayCluster(Cluster):
    """
    This type of cluster uses a python array instead of a set to store members.
    This is much more memory efficient than standard Clusters.

    Insertion complexity:
      * If no keys are provided, insertion is O(1).
      * If keys are inserted in sorted order, insertion is O(log n)
      * If keys are inserted in unknown order, insertion is O(n)

    Merge complexity:
      * Merge runtime is O(n log n)
    """
    def __init__(self, clusters, array_type='i'):
        self._array_type = array_type
        super(SortedArrayCluster, self).__init__(clusters)

    def set_members(self, keys, already_sorted=False):
        if not already_sorted:
            keys = sorted(keys)
        self.members = array(self._array_type, keys)

    def add_members(self, keys, already_sorted=False):
        if not keys:
            return
        if not already_sorted:
            keys = sorted(keys)
        if (not self.members) or keys[0] > self.members[-1]:
            self.members.extend(keys)
        else:
            # need to merge sorted arrays
            old_members = self.members
            self.set_members(heapq.merge(old_members, keys), already_sorted=True)

    def remove_member(self, key):
        pos = bisect.bisect_left(self.members, key)
        if pos < len(self.members) and self.members[pos] == key:
            self.members.pop(pos)

    def remove_members(self, keys):
        if not keys:
            return
        elif len(keys) == 1:
            self.remove_member(keys[0])
        else:
            keys = set(keys)
            self.set_members([member for member in self.members if member not in keys], already_sorted=True)

    def __contains__(self, key):
        pos = bisect.bisect_left(self.members, key)
        if pos >= len(self.members):
            return False
        return self.members[pos] == key


class ClusterPair(object):
    def __init__(self, c1, c2):
        self.c1 = c1
        self.c2 = c2
        self.value = None
        self.update()

    def __lt__(self, other):
        # not a typo! makes the heap algorithm have large stuff first
        return self.value > other.value

    def update(self):
        self.value = self.c1.test_cluster(self.c2) - self.c1.variance - self.c2.variance


class ClusterPairs(object):
    def __init__(self):
        self.pairs = []

    def add(self, pair):
        heapq.heappush(self.pairs, pair)

    def peek(self):
        if self.pairs:
            return self.pairs[0]
        else:
            return None

    def reprioritize(self, pair):
        pair.update()
        heapq.heapify(self.pairs)


class Clusters(object):
    @classmethod
    def correct(cls, variance):
        if variance >= 0.0:
            return variance
        else:
            return 0.0

    def __init__(self, dimension, capacity, cluster_factory=Cluster):
        self.dimension = dimension
        self.capacity = capacity
        self.clusters = []
        self.pairs = ClusterPairs()
        self.cluster_factory = cluster_factory
        for i in range(capacity):
            self.clusters.append(self.cluster_factory(self))
            self._add_pairs()

    def clear(self):
        self.clusters = []
        self.pairs = {}

    def add(self, mass, coords, key):
        self.add_bulk([(mass, coords, key)])

    def add_bulk(self, items, step=1000):
        """
        Takes a iterable of tuples:
            [(mass, (x, y), key), ...]

        This method can be much faster than calling add() repeatedly.
        """

        # doing individual add/append for each key added is slow,
        # so we delay it until a bunch of things have been added instead.
        cluster_keys = collections.defaultdict(list)

        def _add_members():
            for c, new_members in cluster_keys.items():
                c.add_members(new_members)
                del cluster_keys[c]
        try:
            for i, (mass, coords, key) in enumerate(items):
                if mass == 0:
                    continue

                c = None
                for c in self.clusters:
                    if not (cluster_keys[c] or c.members):
                        # *always* populate empty clusters first
                        break
                else:
                    c = None
                if c:
                    c.add(mass, coords, None)
                    if key is not None:
                        cluster_keys[c].append(key)
                    self._update_pairs(c)
                    continue

                #identify cheapest merge
                merge_pair = self.pairs.peek()
                merge_t = merge_pair and merge_pair.value or MAX_FLOAT

                # find cheapest addition
                addition_c = None
                addition_t = MAX_FLOAT
                for cluster in self.clusters:
                    t = cluster.test(mass, coords)
                    if t < addition_t:
                        addition_c = cluster
                        addition_t = t

                if addition_t <= merge_t:
                    # chose addition
                    addition_c.add(mass, coords, None)
                    if key is not None:
                        cluster_keys[addition_c].append(key)
                    self._update_pairs(addition_c)
                else:
                    # choose merge
                    c1 = merge_pair.c1
                    c2 = merge_pair.c2
                    if c1.mass < c2.mass:
                        (c1, c2) = (c2, c1)
                    c1.add_cluster(c2)
                    cluster_keys[c1].extend(cluster_keys.pop(c2, []))
                    if key is not None:
                        cluster_keys[c2] = [key]
                    c2.set(mass, coords, None)
                if i % step == 0:
                    _add_members()
        finally:
            _add_members()

    def results(self):
        return self.clusters[:]

    def reduce(self, max_variance, min_clusters):
        if len(self.clusters) <= min_clusters:
            return

        count = len(self.clusters)
        total_variance = sum([c.variance for c in self.clusters])
        total_mass = sum([c.mass for c in self.clusters])

        while count > min_clusters:
            if count == 1:
                for c in self.clusters:
                    if not c.removed:
                        c.removed = True
            else:
                merge_pair = self.pairs.peek()
                c1 = merge_pair.c1
                c2 = merge_pair.c2
                if c1.mass < c2.mass:
                    c1, c2 = c2, c1
                if max_variance >= 0:
                    diff = c1.test_cluster(c2) - c1.variance - c2.variance
                    total_variance += diff
                    if total_variance / total_mass > max_variance:
                        # stop here, we are going to exceed maximum
                        break
                c1.add_cluster(c2)
                self._update_pairs(c1)
                self._remove_pairs(c2)
                c2.removed = True
            count -= 1

        # remove dead clusters
        for i, c in reversed(list(enumerate(self.clusters))):
            if c.removed:
                self.clusters.pop(i)

        # remove dead pairs
        for c in self.clusters:
            for i, pair in reversed(list(enumerate(c.pairs))):
                if pair.c1.removed or pair.c2.removed:
                    c.pairs.pop(i)

    def _add_pairs(self):
        cj = self.clusters[-1]
        for ci in self.clusters:
            if ci == cj:
                continue
            pair = ClusterPair(ci, cj)
            ci.pairs.append(pair)
            cj.pairs.append(pair)
            self.pairs.add(pair)

    def _update_pairs(self, cluster):
        pairs = cluster.pairs
        limit = len(self.clusters) - 1
        for i, pair in enumerate(pairs):
            if i >= limit:
                break
            if pair.c1.removed or pair.c2.removed:
                continue
            self.pairs.reprioritize(pair)

    def _remove_pairs(self, cluster):
        for pair in cluster.pairs:
            if pair.c1.removed or pair.c2.removed:
                continue
            self.pairs.remove(pair)
