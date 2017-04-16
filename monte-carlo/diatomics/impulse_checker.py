import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
from potential_wrapper import potential

import numpy as np
import vegas
from functools import partial

from time import time

# ---------------------
# CONSTANTS
amu = 9.10938291 * 10**(-31) # amu to kg
da = 1.66053904 * 10**(-27) # da to kg
k = 1.38064852 * 10**(-23) # J / K
htoj = 4.35974417 * 10**(-18) # hartree to J
avogadro = 6.022 * 10**(23) # mol^-1
gas_constant = 8.314
alu = 6.29 * 10**(-11) # m to atomic units

mu = (40. * 44.) / 84. * da / amu
print 'mu: {0}'.format(mu)
# ---------------------

def hamiltonian(J, R, pR):
    return J**2 / (2 * mu * R**2) + pR**2 / (2 * mu) + potential(R, np.pi/2)

def kinetic_part(J, R, pR):
    return pR**2 / (2 * mu)

def impulse_integrand(x, temperature):
    # x = [J, R, pR]
    kinetic_part_value = kinetic_part(*x)
    hamiltonian_value = hamiltonian(*x)

    if hamiltonian_value < 0:
        return np.exp(-kinetic_part_value * htoj / (k * temperature))
    else:
        return 0

limits = [[0, 100], # J
          [5, 100], # R
          [-50, 50], # pR
         ]

def integration(temperature):
    _impulse_integrand = partial(impulse_integrand, temperature = temperature)

    integ = vegas.Integrator(limits)
    impulse_integral = integ(_impulse_integrand, nitn = 20, neval = 3 * 10**4)
    print 'impulse_integral (atomic units) = %s Q = %.2f' % (impulse_integral, impulse_integral.Q)
    #print impulse_integral.summary()
    return impulse_integral.mean / (2 * np.pi)**3

integral = integration(temperature = 200)
print '1/h^3 integral: {0}'.format(integral)

def stat_summ(temperature):
    return (mu * k * temperature / htoj / (2 * np.pi))**(1.5)

sm = stat_summ(temperature = 200)
print 'stat summ: {0}'.format(sm)

print 'relation: {0}'.format(integral / sm)




