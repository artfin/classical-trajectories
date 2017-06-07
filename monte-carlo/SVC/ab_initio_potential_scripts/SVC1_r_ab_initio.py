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
complex_mass = (40. * 12.) / 52. * atomic_mass_unit # kg
c_o_distance = 1.16 * 10**(-10) # m

planck_constant = 1.054571800 * 10**(-34) # joules * s
k = 1.38064852 * 10**(-23) # J / k
htoj = 4.35974417 * 10**(-18) # hartree to Joules
avogadro = 6.022 * 10**(23)
length_unit = 5.291772 * 10**(-11) # atomic units to m
inertia_tensor = 2 * 16 * atomic_mass_unit * c_o_distance**2

def integrand(x, Temperature):
    # x = [R, theta]
    u = potential(x[0], x[1]) * htoj
    du_dtheta = derivative_Theta(x[0], x[1])   
    
    value = np.exp(-u / (k * Temperature)) * (du_dtheta)**2 * x[0]**2 * np.sin(x[1])
    return value

def cycle(T):
    _integrand = partial(integrand, Temperature = T)

    integ = vegas.Integrator([[0., 50.], [0., np.pi]])
    start = time()
    result = integ(_integrand, nitn = 10, neval = 10**4)
	
    print 'result = %s Q = %.2f' % (result, result.Q)

    SVC_correction = np.pi * avogadro * result.mean * length_unit**3 / 24 / (k * T)**3 * htoj**2 * planck_constant**2 / inertia_tensor * 10**6 # 10**6 to convert m3/mol to cm3/mol
    print 'Temperature: %d; SVC correction: %.8f' % (T, SVC_correction)
    print 'Time needed: {0}'.format(time() - start)
    
    return SVC_correction

def save_data(temperatures, svc_corrections):
    filename = '/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/SVC/data/ab-initio/SVC1r_long.dat'

    with open(filename, mode = 'a') as out:
        for temperature, svc_correction in zip(temperatures, svc_corrections):
            out.write(str(temperature) + ' ' + str(svc_correction) + '\n')

# initialization(200)
temperatures = [350 + i for i in range(1, 150)]
#temperatures = [213., 223., 242., 248.2, 262., 273.2, 276., 288.2, 290., 295., 296., 296.15, 300., 303.15, 303.2, 310.0, 313.2, 320.0, 322.85, 323.1, 330.0, 333.15, 363.15, 365., 400., 425., 450., 475.]
svc_corrections = [cycle(temperature) for temperature in temperatures]

save_data(temperatures, svc_corrections)


