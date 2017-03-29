from sympy import symbols
from sympy import sin, cos, tan, lambdify
import numpy as np
import vegas
from functools import partial
from time import time

import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
from potential_wrapper import potential


# alpha == varphi, beta == theta
# Jx = J * cos(alpha) * sin(beta)
# Jy = J * sin(alpha) * sin(beta)
# Jz = J * cos(beta)

J, alpha, beta = symbols('J alpha beta')
R, theta = symbols('R theta')
mu1, mu2 = symbols('mu1 mu2')
pR, pT = symbols('p_R p_theta')
l = symbols('l')

pseudo_kinetic_energy = pR**2 / (2 * mu2)

kinetic_energy = pR**2 / (2 * mu2) + (1 / (2 * mu2 * R**2) + 1 / (2 * mu1 * l**2)) * pT**2 - pT * J * sin(alpha) * sin(beta) / (mu2 * R**2) + (J * sin(alpha) * sin(beta))**2 / (2 * mu2 * R**2) + (J * cos(alpha) * sin(beta))**2 / (2 * mu2 * R**2) + (J * cos(beta))**2 / (2 * sin(theta)**2) * (cos(theta)**2 / (mu2 * R**2) + 1 / (mu1 * l**2)) + J * cos(alpha) * sin(beta) * cos(beta) / (mu2 * R**2 * tan(theta))
kinetic_energy = kinetic_energy.subs({'mu1': 14579., 'mu2': 36440., 'l': 4.398})
_kinetic_energy = lambdify((J, alpha, beta, R, pR, theta, pT), kinetic_energy)

pseudo_kinetic_energy = pseudo_kinetic_energy.subs({'mu2': 36440.})
_pseudo_kinetic_energy = lambdify((J, alpha, beta, R, pR, theta, pT), pseudo_kinetic_energy)
pseudo_hamiltonian = lambda J, alpha, beta, R, pR, theta, pT: \
        _pseudo_kinetic_energy(J, alpha, beta, R, pR, theta, pT)

hamiltonian = lambda J, alpha, beta, R, pR, theta, pT: \
    _kinetic_energy(J, alpha, beta, R, pR, theta, pT) + potential(R, theta)

k = 1.38064852 * 10**(-23) # J/K
htoj = 4.35974417 * 10**(-18) # hartree to Joules
R = 8.314 # J / mol / K
avogadro = 6.022 * 10**(23) # mol^-1

def integrand(x, Temperature):
    # x = [J, alpha, beta, R, pR, theta, pT]
    hamiltonian_value = hamiltonian(*x)
    pseudo_hamiltonian_value = pseudo_hamiltonian(*x)

    if hamiltonian_value < 0:
        return x[0]**2 * np.sin(x[2]) * np.exp(- hamiltonian_value * htoj / (k * Temperature))
    else:
        return 0.0

limits = [[0., 50], # J
	  [0., 2 * np.pi], # alpha (J varphi)
	  [0, np.pi], # beta (J theta)
	  [0, 50.], # R
          [-100., 100.], # pR
	  [0, np.pi], # theta 
	  [-50., 50.], # pT
]

h = 6.626070040*10**(-34)
h_bar = h / (2 * np.pi)

atomic_mass_unit = 9.1093826 * 10**(-31) # kg
atomic_length_unit = 5.291772 * 10**(-11) # m

amu = 1.66053904 * 10**(-27) # amu to kg
complex_mass = 84 * amu
co2_mass = 44 * amu
ar_mass = 40 * amu
m1 = 16 * amu
r0 = 116.3 * 10**(-12) # c=o distance in m
pressure_coeff = 9.869 * 10**(-6) # between pascals and atmospheres

def cycle(T):
    _integrand = partial(integrand, Temperature = T)

    integ = vegas.Integrator(limits)
    
    start = time()
    result = integ(_integrand, nitn = 50, neval = 2 * 10**5)
    print 'Time needed: {0}'.format(time() - start)
    print 'result = %s Q = %.2f' % (result, result.Q)
    return result.mean

def eval_constant(Temperature, integral):

    print 'Temperature: {0}'.format(Temperature)
    Q_Ar = (2 * np.pi * ar_mass * k * Temperature / h**2)**(1.5) / avogadro
    print 'Q Ar: {0}'.format(Q_Ar)
    
    # symmetry number = 2
    Q_CO2 = 16 * np.pi**2 * k * Temperature / h**2 * m1 * r0**2 * (2 * np.pi * co2_mass * k * Temperature / h**2)**(1.5) / avogadro
    print 'Q_CO2: {0}'.format(Q_CO2)

    Q_complex = (2 * np.pi * complex_mass * k * Temperature / h**2)**(1.5) / avogadro
    print 'Q complex: {0}'.format(Q_complex)

    pre_constant = Q_complex / Q_Ar / Q_CO2 / (R * Temperature)
    print 'pre_constant: {0}'.format(pre_constant)

    # 8 * np.pi**2 comes from variable change (Jx, Jy, Jz) -> (J, theta, varphi)
    # (h_bar)**5 comes from expressions in integral
    # h_bar**5 / h**5 gives 1/(2 * pi)**5
    constant = pre_constant * integral / (4 * np.pi**3) * pressure_coeff
    print 'Constant: {0}'.format(constant)

    print '*'*30 + '\n'
    return constant

def save_constants(temperatures, constants):
    with open('full_constants.dat', mode = 'w') as out:
        for temperature, constant in zip(temperatures, constants):
            out.write(str(temperature) + ' ' + str(constant) + '\n')

temperatures = [200 + 5 * i for i in range(20)]
constants = []

for temperature in temperatures:
    integral = cycle(T = temperature)
    constants.append(eval_constant(temperature, integral))

save_constants(temperatures, constants)

