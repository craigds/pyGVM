#!/usr/bin/env python
"""
Simplified port of GVM for python.

GVM homepage: http://www.tomgibara.com/clustering/fast-spatial/java-library
"""
import heapq
import sys

VERSION = (0, 0, 'alpha')

MAX_FLOAT = sys.float_info.max


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

    def clear(self):
        """
        Completely clears this cluster. All points and their associated mass is removed
        """
        self.m0 = 0
        self.m1 = [0] * len(self.m1)
        self.m2 = [0] * len(self.m2)
        self.variance = 0.0
        self.members = []

    def set(self, mass, coords):
        """
        Sets this cluster equal to a single point.
        """
        if mass == 0:
            if len(self.members):
                self.m1 = [0] * len(self.m1)
                self.m2 = [0] * len(self.m2)
        else:
            self.m1 = [mass * coord for coord in coords]
            self.m2 = [mass * coord * coord for coord in coords]
        self.mass = mass
        self.variance = 0.0

    def add(self, mass, coords):
        """
        Adds a point to the cluster.
        """
        if not self.members:
            self.set(mass, coords)
        else:
            if mass != 0:
                self.mass += mass
                for i, coord in enumerate(coords):
                    self.m1[i] += coord * mass
                    self.m2[i] += coord * coord * mass
                self.update()

    def add_cluster(self, cluster):
        """
        Adds the specified cluster to this cluster.
        """
        self.mass += cluster.mass
        self.m1 = [sum(*row) for row in zip(self.m1, cluster.m1)]
        self.m2 = [sum(*row) for row in zip(self.m2, cluster.m2)]
        self.update()

    def test(self, mass, coords):
        """
        Computes this clusters variance if it were to have a new point added to it.
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
            for tm1, tm2, cm1, cm2 in zip(self.m1, self.m2, cluster.m1, cluster.m2):
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
        self.index = None

    def __lt__(self, other):
        return self.value < other.value

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

    def add(self, mass, coords, members):
        if mass == 0:
            return
        if len(self.clusters) < self.capacity:
            cluster = Cluster(self)
            self.clusters.append(cluster)
            cluster.set(mass, coords)
            self._add_pairs()
            cluster.members = members
        else:
            #identify cheapest merge
            merge_pair = self.pairs.peek()
            merge_t = merge_pair.value if merge_pair else MAX_FLOAT

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
                addition_c.add(mass, coords)
                self._update_pairs(addition_c)
                addition_c.members.extend(members)
            else:
                # choose merge
                c1 = merge_pair.c1
                c2 = merge_pair.c2
                if c1.mass < c2.mass:
                    (c1, c2) = (c2, c1)
                c1.members.extend(c2.members)
                c1.add_cluster(c2)
                c2.members = []

    def results(self):
        return self.clusters[:]

    def _add_pairs(self):
        cj = self.clusters[-1]
        for ci in self.clusters:
            if ci == cj:
                continue
            pair = ClusterPair(ci, cj)
            ci.pairs.append(pair)
            cj.pairs.append(pair)

    def _update_pairs(self, cluster):
        pairs = cluster.pairs
        limit = len(self.clusters) - 1
        for i, pair in enumerate(pairs):
            if i >= limit:
                break
            self.pairs.reprioritize(pair)
