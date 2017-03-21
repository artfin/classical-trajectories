import primes
from time import time
from pprint import pprint

__SIZE__ = [1e4]

for size in __SIZE__:
    start = time()
    p = primes.primes(size)
    print 'Time needed: {0}s'.format(time() - start)



