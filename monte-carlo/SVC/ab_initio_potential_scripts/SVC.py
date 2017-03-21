import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')

from potential_wrapper import potential

import numpy as np
import vegas
from time import time
from functools import partial

k = 1.38064852 * 10**(-23) # J/k
htoj = 4.35974417 * 10**(-18) # hartree to Joules
avogadro = 6.022 * 10**(23)
length_unit = 5.291772 * 10**(-11)

def integrand(x, Temperature):
    # x = [R, theta]
    potential_value = potential(x[0], x[1]) * htoj
    value =  (1 - np.exp(- potential_value / (k * Temperature))) * x[0]**2 * np.sin(x[1])
    return value

def initialization(T):
    _integrand = partial(integrand, Temperature = T)

    integ = vegas.Integrator([[3., 50.], [0., np.pi]])
    result = integ(_integrand, nitn = 50, neval = 10**4)
    print 'First integration. result = %s Q = %.2f' % (result, result.Q)

def cycle(T):
    _integrand = partial(integrand, Temperature = T)

    integ = vegas.Integrator([[3., 100.], [0., np.pi]])
    result = integ(_integrand, nitn = 50, neval = 10**4)
    print 'result = %s Q = %.2f' % (result, result.Q)

    SVC = np.pi * avogadro * result.mean * length_unit**3 * 10**6 # 10^6 to convert from m3/mol to cm3/mol
    
    print 'Temperature: %d; SVC: %.5f' % (T, SVC)
    
    return SVC

def save_data(temperatures, svcs):
    file_path = '/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/SVC/data/ab-initio/SVC.dat'
    with open(filename, mode = 'w') as out:
        for temperature, svc in zip(temperatures, svcs):
            out.write(str(temperature) + ' ' + str(svc) + '\n')

initialization(T = 100)

temperatures = [100 + i * 10 for i in range(71)]
svcs = [cycle(temperature) for temperature in temperatures]

save_data(temperatures, svcs)

