from __future__ import print_function
import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
from potential_wrapper import potential

import numpy as np
import scipy.special as sp
import vegas

temperature = 100
k = 1.38064852 * 10**(-23) # J / K
htoj = 4.35974417 * 10**(-18) # hartree to J

print('potential / kT: {0}'.format(potential(6.0, np.pi/2) * htoj / (k * temperature)))

u_kt = potential(6.0, np.pi/2) * htoj / (k * temperature)

def integrand(x):
    # x = [x1, x2, x3]
    summ = np.sum(x**2)
    if summ + u_kt < 0:
        return np.exp(-summ)
    else:
        return 0

integ = vegas.Integrator([[-np.sqrt(-u_kt),np.sqrt(-u_kt)], [-np.sqrt(-u_kt),np.sqrt(-u_kt)], [-np.sqrt(-u_kt),np.sqrt(-u_kt)]]) 
integral = integ(integrand, nitn = 20, neval = 10**4)
print('integral: {0}'.format(integral.mean))

estimate = np.pi**(1.5) * sp.gammainc(1.5, -u_kt)
print('estimate: {0}'.format(estimate))


