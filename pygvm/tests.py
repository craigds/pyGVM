#/usr/bin/env python
import unittest

from pygvm import Clusters

class TestClusters(unittest.TestCase):
    def test_empty_clusters(self):
        cs = Clusters(dimension=2, capacity=2)

        self.assertEqual(len(cs.clusters), 0)

    def test_first_point(self):
        cs = Clusters(dimension=2, capacity=2)

        cs.add(1, (5, 5), 'bert')
        self.assertEqual(len(cs.clusters), 1)
        self.assertEqual(len(cs.clusters[0]), 1)

        self.assertEqual(cs.clusters[0].mass, 1)
        self.assertEqual(cs.clusters[0].center, (5, 5))

    def test_same_point_multiple_times(self):
        # NOTE: this is stupid. If you seed both clusters with the same
        # point, your clustering will basically just be random.
        cs = Clusters(dimension=2, capacity=2)

        cs.add(1, (5, 5), 'tom')
        self.assertEqual(len(cs.clusters), 1)
        self.assertEqual(len(cs.clusters[0]), 1)

        cs.add(1, (5, 5), 'dick')
        self.assertEqual(len(cs.clusters), 2)
        self.assertEqual(len(cs.clusters[0]), 1)
        self.assertEqual(len(cs.clusters[1]), 1)

        cs.add(1, (5, 5), 'harry')
        self.assertEqual(len(cs.clusters), 2)
        self.assertEqual(len(cs.clusters[0]), 2)
        self.assertEqual(len(cs.clusters[1]), 1)

        self.assertEqual(cs.clusters[0].mass, 2)
        self.assertEqual(cs.clusters[1].mass, 1)

if __name__ == '__main__':
    unittest.main()
