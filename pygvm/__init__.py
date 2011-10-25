"""
Simplified port of GVM for python.

GVM homepage: http://www.tomgibara.com/clustering/fast-spatial/java-library
"""
import heapq
import sys

VERSION = (0, 0, 'alpha')

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
        self.members = []
        self.pairs = []

    def center(self):
        if self.mass:
            return [coord / self.mass for coord in self.m1]
        else:
            return self.m1[:]

    def clear(self):
        """
        Completely clears this cluster. All points and their associated mass is removed
        """
        self.m0 = 0
        self.m1 = [0] * len(self.m1)
        self.m2 = [0] * len(self.m2)
        self.variance = 0.0
        self.members = []

    def set(self, mass, coords, members):
        """
        Sets this cluster equal to a single point.
        """
        self.m1 = [mass * coord for coord in coords]
        self.m2 = [mass * coord * coord for coord in coords]
        self.mass = mass
        self.variance = 0.0
        self.members = members

    def add(self, mass, coords, members):
        """
        Adds a point to the cluster.
        """
        if not self.mass:
            self.set(mass, coords, members)
        else:
            if mass != 0:
                self.mass += mass
                for i, coord in enumerate(coords):
                    self.m1[i] += coord * mass
                    self.m2[i] += coord * coord * mass
                self.members.extend(members)
                self.update()

    def add_cluster(self, cluster):
        """
        Adds the specified cluster to this cluster.
        """
        self.add(cluster.mass, cluster.center(), cluster.members)

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
        return u'<Cluster: %d nodes>' % len(self.members)


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

    def __init__(self, dimension, capacity):
        self.dimension = dimension
        self.capacity = capacity
        self.clusters = []
        self.pairs = ClusterPairs()

    def clear(self):
        self.clusters = []
        self.pairs = {}

    def create_cluster(self):
        return Cluster(self)

    def add(self, mass, coords, members):
        if mass == 0:
            return
        if len(self.clusters) < self.capacity:
            cluster = self.create_cluster()
            self.clusters.append(cluster)
            cluster.set(mass, coords, members)
            self._add_pairs()
        else:
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
                addition_c.add(mass, coords, members)
                self._update_pairs(addition_c)
            else:
                # choose merge
                c1 = merge_pair.c1
                c2 = merge_pair.c2
                if c1.mass < c2.mass:
                    (c1, c2) = (c2, c1)
                c1.add_cluster(c2)
                c2.members = []

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
