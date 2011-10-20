
Simplified port of GVM for python.

GVM homepage: http://www.tomgibara.com/clustering/fast-spatial/java-library

Limitations:

 * No implementation of Clusters.reduce() just yet.
 * Not using a proper heap implementation for ClusterPairs, which means it's slower than it could be for large numbers of clusters.
