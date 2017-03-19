# file parker.pyx

import numpy as np
import vegas

cdef extern double potential(double R, double theta)

@vegas.batchintegrand
def f(
