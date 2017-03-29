import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
from potential_wrapper import potential

from sympy import symbols
from sympy import sin, cos, lambdify

import numpy as np
import vegas
from functools import partial

from time import time

amu = 9.10938291 * 10**(-31) # amu to kg

da = 1.66053904 * 10**(-27) # da to kg
co2_mass = 44 * da
ar_mass = 40 * da
complex_mass = 84 * da

reduced_mass = co2_mass * ar_mass / (co2_mass + ar_mass) # in kg
reduced_mass_amu = reduced_mass / amu # reduced mass in amu

J, beta = symbols('J beta')
R, pR = symbols('R pR')
mu = symbols('mu')

kinetic_energy = pR**2 / (2 * mu) + J**2 * sin(beta)**2 / (2 * mu * R**2)
kinetic_energy = kinetic_energy.subs({'mu': reduced_mass_amu}) 
_kinetic_energy = lambdify((J, beta, R, pR), kinetic_energy)

hamiltonian = lambda J, beta, R, pR: _kinetic_energy(J, beta, R, pR) + potential(R, np.pi/2)

k = 1.38064852 * 10**(-23) # J / K
htoj = 4.35974417 * 10**(-18) # hartree to J
avogadro = 6.022 * 10**(23) # mol^-1
pressure_coeff = 9.869 * 10**(-6) # pascals^-1 to atm^-1
gas_constant = 8.314

def integrand(x, temperature):
    # x = [J, beta, R, pR]
    hamiltonian_value = hamiltonian(*x)

    if hamiltonian_value < 0:
        return x[0]**2 * np.sin(x[1]) * np.exp(-hamiltonian_value * htoj / (k * temperature))
    else:
        return 0

limits = [[0, 100], # J
          [0, np.pi], # beta
          [0, 100], # R
          [-5000, 5000], # pR
         ]

h = 6.626070040 * 10**(-34) 

def cycle(temperature):
    _integrand = partial(integrand, temperature = temperature)

    integ = vegas.Integrator(limits)

    start = time()
    result = integ(_integrand, nitn = 50, neval = 3 * 10**5)
    print 'Time needed: {0}'.format(time() - start)
    print 'result = %s Q = %.2f' % (result, result.Q)
    return result.mean

def eval_constant(temperature, integral):

    print 'Temperature: {0}'.format(temperature)
    Q_Ar = (2 * np.pi * ar_mass * k * temperature / h**2)**(1.5) / avogadro
    print 'Q Ar: {0}'.format(Q_Ar)

    Q_CO2 = (2 * np.pi * co2_mass * k * temperature / h**2)**(1.5) / avogadro
    print 'Q CO2: {0}'.format(Q_CO2)

    Q_complex = (2 * np.pi * complex_mass * k * temperature / h**2)**(1.5) / avogadro
    print 'Q translational complex: {0}'.format(Q_complex)

    pre_constant = Q_complex / Q_Ar / Q_CO2 / (gas_constant * temperature)
    print 'Pre constant: {0}'.format(pre_constant)

    constant = pre_constant * integral / (4 * np.pi**3) * pressure_coeff
    print 'Constant: {0}'.format(constant)

    print '*'*30 + '\n'
    return constant

def save_constants(temperatures, constants):
    with open('phase_diatomics.dat', mode = 'w') as out:
        for temperature, constant in zip(temperatures, constants):
            out.write(str(temperature) + ' ' + str(constant) + '\n')

temperatures = [200 + 5 * i for i in range(50)]
constants = []

for temperature in temperatures:
    integral = cycle(temperature)
    constants.append(eval_constant(temperature, integral))

save_constants(temperatures, constants)



