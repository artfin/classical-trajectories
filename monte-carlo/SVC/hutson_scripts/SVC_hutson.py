import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers')

from hutson import potdat4, extpot
import numpy as np
import vegas
from functools import partial

from time import time

# loading parameters in potdat4 subroutine
potdat4()

k = 1.38064852 * 10**(-23) # J/k
htoj = 4.35974417 * 10**(-18) # hartree to Joules
avogadro = 6.022140857 # * 10**(23)
length_unit = 5.291772 # * 10**(-11)

def integrand(x, Temperature):
    # x = [R, theta]
    potential_value = extpot(x[0], np.cos(x[1])) * htoj
    return (1 - np.exp(- potential_value / (k * Temperature))) * np.sin(x[1]) * x[0]**2

def initialization(T):
    _integrand = partial(integrand, Temperature = T)

    integ = vegas.Integrator([[0., 100.], [0., np.pi]])
    result = integ(_integrand, nitn = 50, neval = 10**4)
    print 'First integration. result = %s Q = %.2f' % (result, result.Q)

def cycle(T):
    _integrand = partial(integrand, Temperature = T)
    
    start = time()
    integ = vegas.Integrator([[0., 100.], [0., np.pi]])
    result = integ(_integrand, nitn = 50, neval = 10**4)
    print 'result = %s Q = %.2f' % (result, result.Q)

    SVC = np.pi * avogadro * result.mean * length_unit**3 * 10**(-4) # 10**6 to convert from m3/mol to cm3/mol
    # 10**6 * 10**23 * (10**(-11))**3 gives 10**(-4)
        
    print 'Temperature: %d; SVC %.5f' % (T, SVC)
    print 'Time needed: {0}'.format(time() - start)
    
    return SVC

def save_data(temperatures, svcs):
    filename = '/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/SVC/data/hutson/SVC_from0.dat'
    with open(filename, mode = 'a') as out:
        for temperature, svc in zip(temperatures, svcs):
            out.write(str(temperature) + ' ' + str(svc) + '\n')

initialization(T = 100)

temperatures = [100 + i for i in range(400)]
svcs = [cycle(temperature) for temperature in temperatures]

save_data(temperatures, svcs)


