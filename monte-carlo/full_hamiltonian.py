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

kinetic_energy = pR**2 / (2 * mu2) + (1 / (2 * mu2 * R**2) + 1 / (2 * mu1 * l**2)) * pT**2 - pT * J * sin(alpha) * sin(beta) / (mu2 * R**2) + (J * sin(alpha) * sin(beta))**2 / (2 * mu2 * R**2) + (J * cos(alpha) * sin(beta))**2 / (2 * mu2 * R**2) + (J * cos(beta))**2 / (2 * sin(theta)**2) * (cos(theta)**2 / (mu2 * R**2) + 1 / (mu1 * l**2)) + J * cos(alpha) * sin(beta) * cos(beta) / (mu2 * R**2 * tan(theta))

kinetic_energy = kinetic_energy.subs({'mu1': 14579., 'mu2': 36440., 'l': 4.398})

_kinetic_energy = lambdify((J, alpha, beta, R, theta, pR, pT), kinetic_energy)
hamiltonian = lambda J, alpha, beta, R, theta, pR, pT: \
    _kinetic_energy(J, alpha, beta, R, theta, pR, pT) + potential(R, theta)

k = 1.38064852 * 10**(-23) # J/K
htoj = 4.35974417 * 10**(-18) # hartree to Joules
R = 8.314

def integrand(x, Temperature):
    # x = [J, alpha, beta, R, theta, pR, pT]
    hamiltonian_value = hamiltonian(*x)

    if hamiltonian_value < 0:
        return x[0]**2 * np.sin(x[2]) * np.exp(- hamiltonian_value * htoj / (k * Temperature))
    else:
	return 0.0

limits = [[0., 100], # J
	  [0., 2 * np.pi], # alpha (J varphi)
	  [0, np.pi], # beta (J theta)
	  [0., 50.], # R 
	  [0, np.pi / 2], # theta 
	  [-100., 100.], # pR,
	  [-100., 100.], # pT
]

h = 6.626070040*10**(-34)

atomic_mass_unit = 9.1093826 * 10**(-31) # kg
atomic_time_unit = 2.418884326505 * 10**(-17) # s
atomic_length_unit = 5.291772 * 10**(-11) # m
atomic_momentum_unit = atomic_mass_unit * atomic_length_unit**2 / atomic_time_unit # kg * m**2 / s equals h/(2 * pi)
pR_unit = atomic_mass_unit * atomic_length_unit / atomic_time_unit  # kg * m / s
pT_unit = atomic_mass_unit * atomic_length_unit**2 / atomic_time_unit # kg * m**2 / s
print 'atomic mometum unit: {0}'.format(atomic_momentum_unit)
print 'pR_unit: {0}'.format(pR_unit)
print 'pT_unit: {0}'.format(pT_unit)

amu = 1.66053904 * 10**(-27) # amu to kg
complex_mass = 84 * amu
co2_mass = 44 * amu
ar_mass = 40 * amu
m1 = 16 * amu
r0 = 116.3 * 10**(-12) # c=o distance in m
pressure_coeff = 9.869 * 10**(-6) # between pascals and atmospheres
R = 8.314

def cycle(T):
    _integrand = partial(integrand, Temperature = T)

    integ = vegas.Integrator(limits)
    
    # turns out that neval = 3*10**5 is too small
    start = time()
    result = integ(_integrand, nitn = 50, neval = 3 * 10**5)
    print 'Time needed: {0}'.format(time() - start)
    print 'result = %s Q = %.2f' % (result, result.Q)
    return result.mean

def eval_constant(Temperature, integral):

    print 'Temperature: {0}'.format(Temperature)
    Q_Ar = (2 * np.pi * ar_mass * k * Temperature / h**2)**(1.5)
    print 'Q Ar: {0}'.format(Q_Ar)
    
    # symmetry number = 2
    Q_CO2 = 8 * np.pi**2 * k * Temperature / h**2 * m1 * r0**2 * (2 * np.pi * co2_mass * k * Temperature / h**2)**(1.5)
    print 'Q_CO2: {0}'.format(Q_CO2)

    Q_complex = (2 * np.pi * complex_mass * k * Temperature / h**2)**(1.5)
    print 'Q complex: {0}'.format(Q_complex)

    pre_constant = Q_complex / Q_Ar / Q_CO2 / h**4 * atomic_length_unit * atomic_momentum_unit * pR_unit * pT_unit / (R * Temperature)
    print 'pre_constant: {0}'.format(pre_constant)

    # 8 * np.pi**2 comes from variable change (Jx, Jy, Jz) -> (J, theta, varphi)
    # 4 comes from theta: [0, 2 * np.pi] -> [0, np.pi/2]
    constant = pre_constant * integral * pressure_coeff * 8 * np.pi**2 * 4
    print 'Constant: {0}'.format(constant)

    # h**2 comes from integrand (J**2 in jacobian)
    
    print '*'*30 + '\n'
    return constant

def save_constants(temperatures, constants):
    with open('full_constants.dat', mode = 'w') as out:
        for temperature, constant in zip(temperatures, constants):
            out.write(str(temperature) + ' ' + str(constant) + '\n')

temperatures = [280 + 5 * i for i in range(30)]
constants = []

for temperature in temperatures:
    integral = cycle(T = temperature)
    constants.append(eval_constant(temperature, integral))

save_constants(temperatures, constants)

