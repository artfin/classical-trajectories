import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potentials/')

import hutson
import numpy as np
import vegas
import scipy.special as sp
from functools import partial

from time import time

# loading parameters in potdat4 subroutine
hutson.potdat4()

k = 1.38064852 * 10**(-23) # J/k
htoj = 4.35974417 * 10**(-18) # hartree to Joules
avogadro = 6.022 * 10**(23)
length_unit = 5.291772 * 10**(-11)
R = 8.314
pressure_coeff = 9.869 * 10**(-6) # between pascals and atmospheres
energy_coeff = 1. / 4184. # J to kcal

def integrand(x, Temperature):
    # x = [R, theta]
    potential_value = hutson.extpot(x[0], np.cos(x[1])) * htoj
    return (1 - np.exp(- potential_value / (k * Temperature))) * np.sin(x[1]) * x[0] ** 2

def initialization(T):
    _integrand = partial(integrand, Temperature = T)

    integ = vegas.Integrator([[0., 100.], [0., np.pi]])
    result = integ(_integrand, nitn = 100, neval = 1000)
    print 'First integration. result = %s Q = %.2f' % (result, result.Q)

def cycle(T):
    _integrand = partial(integrand, Temperature = T)

    integ = vegas.Integrator([[0., 100.], [0., np.pi]])
    result = integ(_integrand, nitn = 50, neval = 10**4)
    print 'result = %s Q = %.2f' % (result, result.Q)

    SVC = np.pi * avogadro * result.mean * length_unit**2 * energy_coeff
        
    print 'Temperature: %d; SVC %.5f' % (T, SVC)
    
    return SVC

def save_data(temperatures, svcs):
    with open('data/SVC_hutson.dat', mode = 'w') as out:
        for temperature, svc in zip(temperatures, svcs):
            out.write(str(temperature) + ' ' + str(svc) + '\n')

initialization(T = 100)

temperatures = [100 + i * 10 for i in range(71)]
svcs = [cycle(temperature) for temperature in temperatures]

save_data(temperatures, svcs)


