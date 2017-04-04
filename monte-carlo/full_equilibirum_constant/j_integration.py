import numpy as np
import vegas 
from functools import partial 
from time import time

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

k = 1.38064852 * 10**(-23) # J/K
htoj = 4.35974417 * 10**(-18) # hartree to joules
mu = 36440.

angular_part = lambda J, R: J**2 / (2 * mu * R**2)

def integrand(x, temperature):
    # x = [J, R]
    hamiltonian_value = angular_part(*x)
    return np.exp(- hamiltonian_value * htoj / (k * temperature))

limits = [[-1000, 1000], [-100, 100]]

def cycle(temperature):
    _integrand = partial(integrand, temperature = temperature)

    integ = vegas.Integrator(limits)

    start = time()
    result = integ(_integrand, nitn = 20, neval = 10**4)
    print 'Time neeeded: {0}'.format(time() - start)
    print 'Integral (atomic) = %s Q = %.2f' % (result, result.Q)

cycle(200)

#fig = plt.figure()
#ax = Axes3D(fig)

#X = np.linspace(-1000, 1000, 1000) # J 
#Y = np.linspace(-1000, 1000, 1000) # R
#X, Y = np.meshgrid(X, Y)

#Z = X**2 / (2 * mu * Y**2)

#ax.plot_surface(X, Y, Z, rstride = 1, cstride=1, cmap='hot')

#plt.show()






