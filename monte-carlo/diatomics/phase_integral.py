import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
from potential_wrapper import potential

import numpy as np
import vegas
from functools import partial

from time import time

amu = 9.10938291 * 10**(-31) # amu to kg
da = 1.66053904 * 10**(-27) # da to kg

co2_mass = 44 * da / amu
ar_mass = 40 * da / amu
complex_mass = 84 * da / amu

reduced_mass_amu = co2_mass * ar_mass / (co2_mass + ar_mass) 
print 'reduced_mass_amu: {0}'.format(reduced_mass_amu)

kinetic_part = lambda J, R, pR: pR**2 / (2 * reduced_mass_amu)
angular_part = lambda J, R, pR: J**2 / (2 * reduced_mass_amu * R**2) + potential(R, np.pi/2)
hamiltonian = lambda J, R, pR: kinetic_part(J, R, pR) + angular_part(J, R, pR)

k = 1.38064852 * 10**(-23) # J / K
htoj = 4.35974417 * 10**(-18) # hartree to J
avogadro = 6.022 * 10**(23) # mol^-1
gas_constant = 8.314

alu = 6.29 * 10**(-11) # m to atomic units

def impulse_integrand(x, temperature):
    # x = [J, R, pR]
    kinetic_part_value = kinetic_part(*x)
    hamiltonian_value = hamiltonian(*x)

    if hamiltonian_value < 0:
        return np.exp(-kinetic_part_value * htoj / (k * temperature))
    else:
        return 0

def angular_integrand(x, temperature):
    # x = [J, R, pR]
    angular_part_value = angular_part(*x)
    hamiltonian_value = hamiltonian(*x)
    
    if hamiltonian_value < 0:
        return x[0]**2 * np.exp(-angular_part_value * htoj / (k * temperature))
    else:
        return 0

limits = [[0, 100], # J
          [5, 100], # R
          [-20, 20], # pR
         ]

def cycle(temperature):
    _impulse_integrand = partial(impulse_integrand, temperature = temperature)
    _angular_integrand = partial(angular_integrand, temperature = temperature)

    integ = vegas.Integrator(limits)
    impulse_integral = integ(_impulse_integrand, nitn = 20, neval = 5 * 10**4)
    print 'impulse integral (atomic units) = %s Q = %.2f' % (impulse_integral, impulse_integral.Q)

    angular_integral = integ(_angular_integrand, nitn = 20, neval = 5 * 10**4)
    print 'angular integral (atomic units) = %s Q = %.2f' % (angular_integral, angular_integral.Q)
    
    return impulse_integral.mean, angular_integral.mean

def eval_constant(temperature, impulse_part, angular_part):

    print 'Temperature: {0}'.format(temperature)
    Q_Ar = (ar_mass * k * temperature / htoj / (2 * np.pi))**(1.5)
    print 'Q Ar: {0}'.format(Q_Ar)

    Q_CO2 = (co2_mass * k * temperature / htoj / (2 * np.pi))**(1.5) 
    print 'Q CO2: {0}'.format(Q_CO2)

    Q_complex = (complex_mass * k * temperature / htoj / (2 * np.pi))**(1.5)
    print 'Q translational complex: {0}'.format(Q_complex)
 
    print 'Q_reduced (impulse_part): {0}'.format(impulse_part)

    stat_part = Q_complex * impulse_part / (Q_Ar * Q_CO2)
    print 'statistical sums relation: {0}'.format(stat_part)

    constant = stat_part / (gas_constant * temperature) * avogadro * angular_part * 2 / np.pi * alu**3 
    print 'Constant: {0}'.format(constant)

    print '*'*30 + '\n'
    return constant

def save_constants(temperatures, constants):
    with open('phase_diatomics.dat', mode = 'a') as out:
        for temperature, constant in zip(temperatures, constants):
            out.write(str(temperature) + ' ' + str(constant) + '\n')

temperatures = range(400, 510, 10) 
constants = []
for temperature in temperatures:
    impulse_part, ang_part = cycle(temperature)
    constants.append(eval_constant(temperature, impulse_part, ang_part))

save_constants(temperatures, constants)
