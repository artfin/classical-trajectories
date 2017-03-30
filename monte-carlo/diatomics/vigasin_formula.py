import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
from potential_wrapper import potential
import scipy.special as sp

import numpy as np
import vegas
from functools import partial

from time import time

k = 1.38064852 * 10**(-23) # J / K
htoj = 4.35974417 * 10**(-18) # hartree to J
avogadro = 6.022 * 10**(23) # 1 / mol
R = 8.314 # J / mol / K
pressure_coeff = 9.869 * 10**(-6) # pascals-1 to atm-1
length_unit = 5.291772 * 10**(-11) # alu to m

def integrand(x, Temperature):
    # x = [R]
    potential_value = potential(x[0], np.pi/2) * htoj
    if potential_value < 0:
        return sp.gammainc(1.5, -potential_value / (k * Temperature)) * np.exp(-potential_value / (k * Temperature)) * x[0]**2
    else:
        return 0.0

def determine_sigma():
    x_space = np.linspace(5, 10, 10000)

    for x1, x2 in zip(x_space, x_space[1:]):
        if potential(x1, np.pi/2) > 0 and potential(x2, np.pi/2) < 0:
            return x1

sigma = determine_sigma()
print 'sigma: {0}'.format(sigma)

def cycle(Temperature):
    _integrand = partial(integrand, Temperature = Temperature)
    
    start = time()
    integ = vegas.Integrator([[sigma, 100]])
    result = integ(_integrand, nitn = 50, neval = 10**5)
    print 'Time needed: {0}'.format(time() - start)
    
    print 'result = %s Q = %.2f' % (result, result.Q)
    constant = 2 * np.pi * avogadro / R / Temperature * result.mean * pressure_coeff * length_unit**2 
    print 'Constant %.5f' % constant
    return constant

def save_constants(temperatures, constants):
    with open('vigasin_diatomics.dat', mode = 'w') as out:
        for temperature, constant in zip(temperatures, constants):
            out.write(str(temperature) + ' ' + str(constant) + '\n')

temperatures = [200 + 5 * i for i in range(0, 50)]  
constants = [cycle(temperature) for temperature in temperatures]

save_constants(temperatures, constants)