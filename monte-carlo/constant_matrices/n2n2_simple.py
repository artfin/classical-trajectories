from __future__ import print_function

import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/N2N2')

from n2n2 import potinit, potn2n2

# initializing procedure
potinit()

x = potn2n2(5.0, 0.5, 0.5, 0.5)
print("potential value: {0}".format(x))

# converts cm^-1 to hartree
def potential(r, theta1, theta2, phi):
    return potn2n2(r, theta1, theta2, phi) * 4.55633 * 10**(-6)

import numpy as np
import scipy.special as sp
import vegas
from time import time
from functools import partial

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# CONSTANTS
# ---------------------------------------
k = 1.38064852 * 10**(-23) # J/k
htoj = 4.35974417 * 10**(-18) # hartree to Joules
avogadro = 6.022140 * 10**(23) 
atomic_length_unit = 5.291772 * 10**(-11)
patoatm = 101325.0 # Pa to Atm
R = 8.314 # J / K
# ---------------------------------------

def integrand(x, temperature):
    # x = [r, theta1, theta2, phi]
    r = x[0]
    theta1 = x[1]
    theta2 = x[2]

    potential_value = potential(*x) * htoj

    if potential_value < 0:
        u_kt = potential_value / (k * temperature)
        return r**2 * np.sin(theta1) * np.sin(theta2) * sp.gammainc(3.5, -u_kt) * np.exp(-u_kt)
    else:
        return 0

def simple_constant(temperature):
    _integrand = partial(integrand, temperature = temperature)

    print('*'*30)
    print('Temperature: {0}'.format(temperature))

    start = time()
    integ = vegas.Integrator([[4.4, 45.0], [0.0, np.pi], [0.0, np.pi], [0.0, 2 * np.pi]])
    integral = integ(_integrand, nitn = 10, neval = 3 * 10**4)
    print('integral (atomic units): {0}'.format(integral))
    print('Time needed: {0}'.format(time() - start))

    eq_const = avogadro / (4 * R * temperature) * integral.mean * atomic_length_unit**3 * patoatm
    print('Equilibrium constant: {0}'.format(eq_const))
    print('*'*30)

    return eq_const

def save_constants(filename, temperatures, constants):
    with open(filename, mode = 'w') as out:
        for t, c in zip(temperatures, constants):
            out.write(str(t) + ' ' + str(c) + '\n')

def read_constants(filename):
    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    temperatures = []
    constants = []
    for line in lines:
        if len(line) > 1:
            data = line.split()

            temperatures.append(float(data[0]))
            constants.append(float(data[1]))

    return temperatures, constants



#temperatures = np.linspace(100, 500, 50)
#constants = [simple_constant(_) for _ in temperatures]

#save_constants('simple_n2n2.dat', temperatures, constants)

#temperatures, constants = read_constants('simple_n2n2.dat')

#plt.plot(temperatures, constants, '--', color = 'k')

#patch = mpatches.Patch(color = 'k', label = 'Simple')
#plt.legend(handles = [patch])
#plt.grid()
#plt.show()
