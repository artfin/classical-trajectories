from without_sum import var3
from random import random
import timeit

v = [random() for _ in range(10000)]

def test():
    res = var3(v)

number = 10000

if __name__ == '__main__':
    time = timeit.timeit("test()", setup = 'from __main__ import test', number = number)
    print 'number of calls: {0}; evaluated mean time: {1} microseconds'.format(number, 1e6 * time / number)
