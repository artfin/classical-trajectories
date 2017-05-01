from __future__ import print_function
import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
from potential_wrapper import potential

import numpy as np
import scipy.special as sp
import vegas

from time import time
from functools import partial

k = 1.38064852 * 10**(-23) # J / K
htoj = 4.35974417 * 10**(-18) # hartree to J
h = 6.626070040*10**(-34)
da = 1.66 * 10**(-27) # da to kg
avogadro = 6.022 * 10**(23)
R = 8.314
patoatm = 101325

atomic_mass_unit = 9.1093826 * 10**(-31) # kg
atomic_length_unit = 5.291772 * 10**(-11) # m

m1 = 16 * da # oxygen mass in kg

mu1 = 14579.0
mu2 = 36440.0
l = 4.398
complex_mass = 84 * da / atomic_mass_unit
print('complex_mass: {0}'.format(complex_mass))
co2_mass = 44 * da / atomic_mass_unit
print('co2_mass: {0}'.format(co2_mass))
ar_mass = 40 * da / atomic_mass_unit
print('ar_mass: {0}'.format(ar_mass))

def hamiltonian(x):
    # x = [R, theta, pR, pT, Jx, Jy, Jz]
    R = x[0]
    theta = x[1]
    pR = x[2]
    pT = x[3]
    Jx = x[4]
    Jy = x[5]
    Jz = x[6]
    return pR**2 / (2 * mu2) + (1 / (2 * mu2 * R**2) + 1 / (2 * mu1 * l**2)) * pT**2 - pT * Jy/ (mu2 * R**2) + Jy**2 / (2 * mu2 * R**2) + Jx**2 / (2 * mu2 * R**2) + Jz**2 / (2 * np.sin(theta)**2) * (np.cos(theta)**2 / (mu2 * R**2) + 1 / (mu1 * l**2)) + Jx * Jz/ (mu2 * R**2 * np.tan(theta)) + potential(R, theta)

def simple_integrand(x, temperature):
    # x = [R, theta]
    R = x[0]
    theta = x[1]
    potential_value = potential(R, theta)

    if potential_value < 0:
        u_kt = potential_value * htoj / (k * temperature)
        return R**2 * np.sin(theta) * sp.gammainc(2.5, -u_kt) * np.exp(-u_kt)
    else:
        return 0

def adv_integrand(x, temperature):
    # x = [R, theta, pR, pT, Jx, Jy, Jz]
    hamiltonian_value = hamiltonian(x)
    if hamiltonian_value < 0:
        return np.exp(-hamiltonian_value * htoj / (k * temperature))
    else:
        return 0

# simple form of hamtiltonian
def simple_constant(temp):
    _simple_integrand = partial(simple_integrand, temperature = temp)
    
    print('*'*30)
    print('Temperature: {0}'.format(temp))
    
    start = time()
    integ = vegas.Integrator([[3, 20], [0, np.pi]])
    integral = integ(_simple_integrand, nitn = 20, neval = 1 * 10**4)
    print('integral (atomic units): {0}'.format(integral))
    interim = integral.mean * atomic_length_unit**3
    print('Time needed: {0}'.format(time() - start))

    eq_const = 2 * np.pi * avogadro / (R * temperature) * interim * patoatm 
    print('eq_const: {0}'.format(eq_const))
    print('*'*30 + '\n')

    return eq_const

# full phase integral
def full_phase_constant(temperature):
    _adv_integrand = partial(adv_integrand, temperature = temperature)

    print('*'*30)
    print('Temperature: {0}'.format(temperature))

    start = time()
    integ = vegas.Integrator([[3, 20], [0, np.pi], [-50, 50], [-50, 50], [-100, 100], [-100, 100], [-100, 100]])
    integral = integ(_adv_integrand, nitn = 20, neval = 4 * 10**5)
    print(integral.summary())
    print('integral (atomic units): {0}'.format(integral.mean))
    print('Time needed: {0}'.format(time() - start))

    #interim = 180529672.14
    interim = integral.mean

    qtr_complex = (2 * np.pi * 84 * da * k * temperature / h**2)**(1.5)
    q_ar = (2 * np.pi * 40 * da * k * temperature / h**2)**(1.5)
    q_co2 = (2 * np.pi * 44 * da * k * temperature / h**2)**(1.5) * 8 * np.pi**2 * k * temperature * 16 * da * (116.3 * 10**(-12))**2 / h**2

    print('qtr complex: {0}'.format(qtr_complex))
    print('q_ar: {0}'.format(q_ar))
    print('q_co2: {0}'.format(q_co2))

    eq_const = avogadro / (R * temperature) * qtr_complex / q_ar / q_co2 * interim * patoatm * 8 * np.pi**2 / (2 * np.pi)**5 / 1.8857
    print('eq_const: {0}'.format(eq_const))
    print('*'*30 + '\n')

    return eq_const

def save_constants(filename, temperatures, constants):
    with open('output/' + filename, mode = 'w') as out:
        for temperature, constant in zip(temperatures, constants):
            out.write(str(temperature) + ' ' + str(constant) + '\n')

temperatures = range(200, 201, 1)
# calculating simple constants
simple_constants = [simple_constant(temperature) for temperature in temperatures]
#save_constants('simple_constants.dat', temperatures, simple_constants)

 #calculating super full phase constant
phase_constants = [full_phase_constant(temperature) for temperature in temperatures]
#save_constants('phase_constants.dat', temperatures, phase_constants)



