from __future__ import absolute_import
import sys
from heapq import *

if sys.version < '2.6':
    def merge(*iterables):
        '''
        Lifted from Python 2.6 source for compatibility with Python 2.5.

        Merge multiple sorted inputs into a single sorted output.

        Similar to sorted(itertools.chain(*iterables)) but returns a generator,
        does not pull the data into memory all at once, and assumes that each of
        the input streams is already sorted (smallest to largest).

        >>> list(merge([1,3,5,7], [0,2,4,8], [5,10,15,20], [], [25]))
        [0, 1, 2, 3, 4, 5, 5, 7, 8, 10, 15, 20, 25]

        '''
        _heappop, _heapreplace, _StopIteration = heappop, heapreplace, StopIteration

        h = []
        h_append = h.append
        for itnum, it in enumerate(map(iter, iterables)):
            try:
                next = it.next
                h_append([next(), itnum, next])
            except _StopIteration:
                pass
        heapify(h)

        while 1:
            try:
                while 1:
                    v, itnum, next = s = h[0]   # raises IndexError when h is empty
                    yield v
                    s[0] = next()               # raises StopIteration when exhausted
                    _heapreplace(h, s)          # restore heap condition
            except _StopIteration:
                _heappop(h)                     # remove empty iterator
            except IndexError:
                return
