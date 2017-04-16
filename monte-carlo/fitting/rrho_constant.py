from __future__ import division
from __future__ import print_function

import numpy as np

import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
from potential_wrapper import potential

import matplotlib.pyplot as plt

alu = 5.291 * 10**(-11) # atomic length unit to m

c = 2.998 * 10**10 # cm/s
h = 6.626070040*10**(-34)
k = 1.38064852 * 10**(-23) # J/K
da = 1.66053904 * 10**(-27) # da to kg

mu = 44 * 40 / 84 * da # mu in kg
r = 6.4966 * alu # m -- equilibrium r-value

patoatm = 101325 # atm to pa == pa^-1 to atm^-1

rotational_constant_hz = h / (8 * np.pi**2 * mu * r**2) # (= B * c) in Hz 
print('rotational_constant Hz: {0}'.format(rotational_constant_hz))
rotational_constant_cm = h / (8 * np.pi**2 * mu * c * r**2) 
print('rotational_constant cm^-1: {0}'.format(rotational_constant_cm))

_lambda = lambda temperature: h / (2 * np.pi * mu * k * temperature)**(0.5)
print('lambda 200K: {0}'.format(_lambda(200)))

De = 195.64 * 1.98630 * 10**(-23) # j
omega_hz = 8.08088 * 10**(11) # Hz
omega_cm = omega_hz / c
print('omega Hz: {0}'.format(omega_hz))
print('omega cm^-1: {0}'.format(omega_cm))

braces = lambda temperature: np.exp(De / (k * temperature))
print('braces 200K: {0}'.format(braces(200)))

print('kt/hcb: {0}'.format(k * 200 / (h * rotational_constant_hz)))

equilibrium_constant = lambda temperature: _lambda(temperature)**3 * (k * temperature) / ( (h * omega_hz) * (h * rotational_constant_hz) ) * braces(temperature) * patoatm

for temperature in [t*100 for t in range(1, 6)]:
    print('equilibrium_constant({0}K): {1}'.format(temperature, equilibrium_constant(temperature)))

# SUPER SIMPLE FORMULA
mu = 44 * 40 / 84 # in amu
equilibrium_constant = lambda temperature: 18.86 / (omega_cm * rotational_constant_cm * temperature**(0.5) * mu**(1.5)) * (np.exp(De / (k * temperature)) - 1 - De / (k * temperature))

for t in [100, 200, 300, 400]:
    print('equilibrium_constant {0}: {1}'.format(t, equilibrium_constant(t)))

multiplier = h / ( (2 * np.pi)**(1.5) * k**0.5 * c**2 * da**(1.5)) * 101325
print(multiplier)

#with open('rrho_constant.dat', mode = 'w') as out:
    #for temperature in range(100, 500):
        #out.write(str(temperature) + ' ' + str(equilibrium_constant(temperature)) + '\n')

