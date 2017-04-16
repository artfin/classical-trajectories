from __future__ import print_function

import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
from potential_wrapper import potential

import numpy as np
import scipy.special as sp
import vegas
from functools import partial

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
print('mu: {0}'.format(mu))
# ---------------------

temperature = 100
R = 6.0
pR = -3.0
u_kt = potential(6.0, np.pi/2) * htoj / (k * temperature)

def hamiltonian(Jx, Jy):
    return Jx**2 / (2 * mu * R**2) + Jy**2 / (2 * mu * R**2) + potential(R, np.pi/2)

def integrand(x):
    # x = [x1, x2, x3]
    summ = np.sum(x**2)
    if summ + u_kt < 0:
        return np.exp(-summ)
    else:
        return 0

integ = vegas.Integrator([[-np.sqrt(-u_kt),np.sqrt(-u_kt)], [-np.sqrt(-u_kt),np.sqrt(-u_kt)], [-np.sqrt(-u_kt),np.sqrt(-u_kt)]]) 
integral = integ(integrand, nitn = 20, neval = 10**4)
print('integral: {0}'.format(integral.mean))

def full_integrand(x):
    # x = [pR, Jx, Jy]
    hamiltonian_value = hamiltonian(x[1], x[2])
    if hamiltonian_value < 0:
        return np.exp(-x[0]**2 / (2 * mu) - x[1]**2 / (2 * mu * R**2) - x[2]**2 / (2 * mu * R**2))
    else: 
        return 0

print('-2*u_kt_mu: {0}'.format(-2*u_kt*mu))
print('-2*u_kt_mu*R**2: {0}'.format(-2*u_kt*mu*R**2))

integ = vegas.Integrator([[-np.sqrt(- 2*u_kt*mu), np.sqrt(-2*u_kt*mu)],
                          [-np.sqrt(-2*u_kt*mu*R**2), np.sqrt(-2*u_kt*mu*R**2)],
                          [-np.sqrt(-2*u_kt*mu*R**2), np.sqrt(-2*u_kt*mu*R**2)]])
integral = integ(integrand, nitn = 20, neval = 10**6)
print(integral.summary())
print('integral: {0}'.format(integral.mean))

