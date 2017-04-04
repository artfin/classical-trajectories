import numpy as np
import vegas
from functools import partial
from time import time

import matplotlib.pyplot as plt

import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
from potential_wrapper import potential

k = 1.38064852 * 10**(-23) # J/K
htoj = 4.35974417 * 10**(-23) # hartree to joules
mu2 = 36440.

kinetic_part = lambda pR: pR**2 / (2 * mu2)

print kinetic_part(3.)

def integrand(x, temperature):
    # x = [pR]
    hamiltonian_value = kinetic_part(*x)
    return np.exp(- hamiltonian_value * htoj / (k * temperature))

limits = [[-5000, 5000]]

def cycle(temperature):
    _integrand = partial(integrand, temperature = temperature)

    integ = vegas.Integrator(limits)

    start = time()
    result = integ(_integrand, nitn = 50, neval = 10**3)
    print 'Time needed: {0}'.format(time() - start)
    print 'result = %s Q = %.2f' % (result, result.Q)
    return result.mean

cycle(200)

predicted = np.sqrt(2 * np.pi * mu2 * k * 200)
print 'predicted: {0}'.format(predicted)

#plt.plot(x, y3)
#x = np.linspace(-1e4, 1e4, 1e3)
#y0 = [np.exp(-kinetic_part(_x) * htoj / k / 100) for _x in x]
#y1 = [np.exp(-kinetic_part(_x) * htoj / k / 200) for _x in x]
#y2 = [np.exp(-kinetic_part(_x) * htoj / k / 300) for _x in x]
#y3 = [np.exp(-kinetic_part(_x) * htoj / k / 400) for _x in x]
#plt.plot(x, y0, color = 'magenta')
#plt.plot(x, y1, color = 'blue')
#plt.plot(x, y2, color = 'red')
#plt.plot(x, y3, color = 'cyan')
#plt.show()
