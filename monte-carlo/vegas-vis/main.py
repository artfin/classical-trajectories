from __future__ import print_function

import vegas
import numpy as np

import matplotlib.pyplot as plt

def f(x):
    return np.exp(-x[0]**2 - x[1]**2) 

integ = vegas.Integrator([[0, 1], [0, 1]])

#m = vegas.AdaptiveMap([[0, 1], [0, 1]], ninc = 5)

#ny = 1000
#y = np.random.uniform(0., 1., (ny, 2)) # 1000 random y's

#x = np.empty(y.shape, float)
#jac = np.empty(y.shape[0], float)
#f2 = np.empty(y.shape[0], float)

#print('initial grid:')
#print(m.settings())

#for itn in range(5):
    #m.map(y, x, jac)
    ##m.show_grid()

    #for j in range(ny):
        #f2[j] = (jac[j] * f(x[j])) ** 2

    #m.add_training_data(y, f2)
    #m.adapt(alpha = 1.5)

    #print('iteration: {0}'.format(itn))
    #print(m.settings())


ax = plt.subplot(1, 3, 1)
integ(f, nitn = 1, neval = 1e3)
for x, wgt in integ.random():
    #print('x[0]: {0}; x[1]: {1}; wgt: {2}'.format(x[0], x[1], wgt))
    ax.scatter(x[0], x[1], color = 'r', marker = '*')
ax.grid()

ax = plt.subplot(1, 3, 2)
integ(f, nitn = 40, neval = 1e3)
for x, wgt in integ.random():
    #print('x[0]: {0}; x[1]: {1}; wgt: {2}'.format(x[0], x[1], wgt))
    ax.scatter(x[0], x[1], color = 'b', marker = '*')
ax.grid()

ax = plt.subplot(1, 3, 3)
integ(f, nitn = 100, neval = 1e3)
for x, wgt in integ.random():
    ax.scatter(x[0], x[1], color = 'g', marker = '*')
ax.grid()

plt.show()
