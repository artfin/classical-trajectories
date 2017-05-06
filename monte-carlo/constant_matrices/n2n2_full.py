from __future__ import print_function

import sys
from n2n2_hamiltonian.ham import ke

sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/N2N2/')
from n2n2 import potinit, potn2n2

# initializing procedure
potinit()

# converts cm^-1 to hartree
def potential(r, theta1, theta2, phi):
    return potn2n2(r, theta1, theta2, phi) * 4.55633 * 10**(-6)

import numpy as np
import vegas
from time import time
from functools import partial

# CONSTANTS
# ---------------------------------------
k = 1.38064852 * 10**(-23) # J / K
htoj = 4.35974417 * 10**(-18) # hartree to J
h = 6.626070040*10**(-34)
da = 1.660539040 * 10**(-27) # da to kg
avogadro = 6.022 * 10**(23)
R = 8.314
patoatm = 101325

atomic_mass_unit = 9.1093826 * 10**(-31) # kg
atomic_length_unit = 5.291772 * 10**(-11) # m
# ---------------------------------------

def integrand(x, temperature):
    # x = [q1, q2, q3, q4, p1, p2, p3, p4, Jx, Jy, Jz]
    theta1 = x[0] # q1
    theta2 = x[2] # q3
    phi = x[1] # q2 
    r = x[3] # q4
    
    kinetic_energy = ke(*x)
    hamiltonian_value = kinetic_energy + potential(r, theta1, theta2, phi)
    print('ham value: {0}; exp: {1}'.format(hamiltonian_value, np.exp(-hamiltonian_value * htoj / (k * temperature))))

    if hamiltonian_value < 0:
        print('*'*30)
        print('kinetic_energy: {0}: x: {1}'.format(ke(*x), x))
        print('*'*30)
        return np.exp(- hamiltonian_value * htoj / (k * temperature))
    else:
        return 0

def constant(temperature):
    _integrand = partial(integrand, temperature = temperature)

    print('*'*30)
    print('Temperature: {0}'.format(temperature))

    start = time()
    integ = vegas.Integrator([[0.0, np.pi], [0.0, 2 * np.pi], [0.0, np.pi], [4.4, 45.0], [-50.0, 50.0], [-50.0, 50.0], [-50.0, 50.0], [-50.0, 50.0], [-50.0, 50.0], [-50.0, 50.0], [-50.0, 50.0]])
    integral = integ(_integrand, nitn = 20, neval = 10**5)
    print(integral.summary())
    print('integral (atomic units): {0}'.format(integral.mean))
    print('Time needed: {0}'.format(time() - start))

temperatures = range(200, 201, 1)
constants = [constant(_) for _ in temperatures]















