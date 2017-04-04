import numpy as np
import vegas
from functools import partial
from time import time

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
from potential_wrapper import potential

k = 1.38064852 * 10**(-23) # J/K
htoj = 4.35974417 * 10**(-18) # hartree to joules
mu2 = 38193.86 

h = 6.626 * 10**(-34)
h_bar = h / (2 * np.pi)

atomic_mass_unit = 9.109 * 10**(-31) # kg
atomic_length_unit = 5.2917 * 10**(-11) # m
atomic_time_unit = 2.418 * 10**(-17) # s

pr_unit = h_bar / atomic_length_unit
print 'pr_unit (1): {0}'.format(pr_unit)

pr_unit = atomic_mass_unit * atomic_length_unit / atomic_time_unit 
print 'pr_unit (2): {0}'.format(pr_unit)

kinetic_part = lambda pR: pR**2 / (2 * mu2)

def integrand(x, temperature):
    # x = [pR]
    hamiltonian_value = kinetic_part(*x)
    return np.exp(- hamiltonian_value * htoj / (k * temperature))

limits = [[-20, 20]]

def cycle(temperature):
    _integrand = partial(integrand, temperature = temperature)

    integ = vegas.Integrator(limits)

    start = time()
    result = integ(_integrand, nitn = 20, neval = 10**4)
    print 'Time needed: {0}'.format(time() - start)
    print 'integral (atomic) = %s Q = %.2f' % (result, result.Q)
    print 'integral (CI) = %s' % (result * pr_unit)
    return result.mean

cycle(100)
cycle(200)
cycle(300)
cycle(400)

predicted = np.sqrt(2 * np.pi * mu2 * k * atomic_mass_unit* 100)
print 'predicted 100K: {0}'.format(predicted)
predicted = np.sqrt(2 * np.pi * mu2 * k * atomic_mass_unit* 200)
print 'predicted 200K: {0}'.format(predicted)
predicted = np.sqrt(2 * np.pi * mu2 * k * atomic_mass_unit* 300)
print 'predicted 300K: {0}'.format(predicted)
predicted = np.sqrt(2 * np.pi * mu2 * k * atomic_mass_unit* 400)
print 'predicted 400K: {0}'.format(predicted)

#sigma = lambda temperature: np.sqrt(mu2 * k * temperature / htoj)
#print 'FWHM 100K: {0}; 3xsigma: {1}'.format(2.35482 * sigma(100), 3 * sigma(100))
#print 'FWHM 200K: {0}; 3xsigma: {1}'.format(2.35482 * sigma(200), 3 * sigma(200))
#print 'FWHM 300K: {0}; 3xsigma: {1}'.format(2.35482 * sigma(300), 3 * sigma(300))
#print 'FWHM 400K: {0}; 3xsigma: {1}'.format(2.35482 * sigma(400), 3 * sigma(400))

#x = np.linspace(-1e1, 1e1, 1e3)
#y0 = [np.exp(-kinetic_part(_x) * htoj / k / 100) for _x in x]
#y1 = [np.exp(-kinetic_part(_x) * htoj / k / 200) for _x in x]
#y2 = [np.exp(-kinetic_part(_x) * htoj / k / 300) for _x in x]
#y3 = [np.exp(-kinetic_part(_x) * htoj / k / 400) for _x in x]
#plt.plot(x, y0, color = 'magenta')
#plt.plot(x, y1, color = 'blue')
#plt.plot(x, y2, color = 'red')
#plt.plot(x, y3, color = 'cyan')

#magenta_patch = mpatches.Patch(color = 'magenta', label = '100K')
#blue_patch = mpatches.Patch(color = 'blue', label = '200K')
#red_patch = mpatches.Patch(color = 'red', label = '300K')
#cyan_patch = mpatches.Patch(color = 'cyan', label = '400K')

#plt.legend(handles = [magenta_patch, blue_patch, red_patch, cyan_patch])

#plt.savefig('impulse_plot.png')
#plt.show()
