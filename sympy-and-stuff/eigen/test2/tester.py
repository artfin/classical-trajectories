from ham import __rhs__ as rhs
import numpy as np
from time import time

# preparing np.array() for right-hand sides to be put in
out = np.array([0,0,0,0,0,0], dtype = 'double')

attempts = int(1e5)

start = time()
for i in range(attempts):
    rhs(out, 1.0, 0.55, 10.0, 2.0, 10.0, 0.55, 0.1)
print 'Time needed: {0} microseconds'.format((time() - start) / attempts * 10e6)


