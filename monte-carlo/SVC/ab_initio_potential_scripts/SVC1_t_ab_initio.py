import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/derivative_wrapper')
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')

from derivatives_wrapper import derivative_R, derivative_Theta
from potential_wrapper import potential

import numpy as np
import vegas
from time import time
from functools import partial

atomic_mass_unit = 1.660539040 * 10**(-27) # kg
complex_mass = (40. * 44.) / 84. * atomic_mass_unit # kg

planck_constant = 1.054571800 * 10**(-34) # joules * s
k = 1.38064852 * 10**(-23) # J / k
htoj = 4.35974417 * 10**(-18) # hartree to Joules
avogadro = 6.022 * 10**(23)
length_unit = 5.291772 * 10**(-11) # atomic units to m

def integrand(x, Temperature):
    # x = [R, theta]
    u = potential(x[0], x[1]) * htoj
    du_dr = derivative_R(x[0], x[1])
    
    value = np.exp(-u / (k * Temperature)) * (du_dr)**2 * x[0]**2 * np.sin(x[1])

    return value

def initialization(T):
    _integrand = partial(integrand, Temperature = T)

    integ = vegas.Integrator([[3., 20.], [0., np.pi]])
    start = time()
    result = integ(_integrand, nitn = 20, neval = 40000)
    print 'Time needed: {0}'.format(time() - start)
    print 'First integration, result = %s Q = %.2f' % (result, result.Q)

def cycle(T):
    _integrand = partial(integrand, Temperature = T)

    integ = vegas.Integrator([[0., 100.], [0., np.pi]])
    start = time()
    result = integ(_integrand, nitn = 50, neval = 10**4)
	
    print 'result = %s Q = %.2f' % (result, result.Q)

    SVC_correction = np.pi * avogadro * result.mean * length_unit / 12 / (k * T)**3 * htoj**2 * planck_constant**2 / complex_mass * 10**6 # 10**6 to convert m3/mol to cm3/mol
    print 'Temperature: %d; SVC correction: %.8f' % (T, SVC_correction)
    print 'Time needed: {0}'.format(time() - start)
    
    return SVC_correction

def save_data(temperatures, svc_corrections):
    filename = '/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/SVC/data/ab-initio/SVC1t.dat'
    with open(filename, mode = 'w') as out:
        for temperature, svc_correction in zip(temperatures, svc_corrections):
            out.write(str(temperature) + ' ' + str(svc_correction) + '\n')

# initialization(200)
#temperatures = [100 + i for i in range(101)]
temperatures = [213., 223., 242., 248.2, 262., 273.2, 276., 288.2, 290., 295., 296., 296.15, 300., 303.15, 303.2, 310.0, 313.2, 320.0, 322.85, 323.1, 330.0, 333.15, 363.15, 365., 400., 425., 450., 475.]

svc_corrections = [cycle(temperature) for temperature in temperatures]

save_data(temperatures, svc_corrections)


