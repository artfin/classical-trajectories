from adding_types import var2
from random import random
import timeit

v = [random() for _ in range(10000)]

def test():
    """Stupid test function."""
    res = var2(v)

number = 10000

if __name__ == '__main__': 
    time = timeit.timeit("test()", setup = 'from __main__ import test', number = number)
    print 'number of calls: {0}; evaluated mean time: {1} microseconds'.format(number, 1000000 * time / number)
