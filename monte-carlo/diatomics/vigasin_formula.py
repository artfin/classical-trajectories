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
patoatm = 101325 # pa^-1 to atm^-1
length_unit = 5.291772 * 10**(-11) # alu to m

def integrand(x, Temperature):
    # x = [R]
    potential_value = potential(x[0], np.pi/2) * htoj
    if potential_value < 0:
        return sp.gammainc(1.5, -potential_value / (k * Temperature)) * np.exp(-potential_value / (k * Temperature)) * x[0]**2
    else:
        return 0.0

#def determine_sigma():
    #x_space = np.linspace(5, 10, 10000)

    #for x1, x2 in zip(x_space, x_space[1:]):
        #if potential(x1, np.pi/2) > 0 and potential(x2, np.pi/2) < 0:
            #return x1

#sigma = determine_sigma()
#print 'sigma: {0}'.format(sigma)

def cycle(Temperature):
    print 'Temperature: {0}'.format(Temperature)

    _integrand = partial(integrand, Temperature = Temperature)
    
    start = time()
    integ = vegas.Integrator([[0, 50]])
    result = integ(_integrand, nitn = 10, neval = 10**4)
    print 'Time needed: {0}'.format(time() - start)
    
    print 'result = %s Q = %.2f' % (result, result.Q)
    constant = 4 * np.pi * avogadro / R / Temperature * result.mean * patoatm * length_unit**3 
    print 'Constant %.5f' % constant
    return constant

def save_constants(temperatures, constants):
    with open('vigasin_diatomics.dat', mode = 'a') as out:
        for temperature, constant in zip(temperatures, constants):
            out.write(str(temperature) + ' ' + str(constant) + '\n')

temperatures = range(100, 150)
constants = [cycle(temperature) for temperature in temperatures]

#for temperature, constant in zip(temperatures, constants):
    #print 'temperature: {0}; constant: {1}'.format(temperature, constant)

#save_constants(temperatures, constants)
