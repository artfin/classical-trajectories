from __future__ import print_function
import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
from potential_wrapper import potential

import numpy as np
import scipy.special as sp
import vegas

temperature = 100
k = 1.38064852 * 10**(-23) # J / K
htoj = 4.35974417 * 10**(-18) # hartree to J
h = 6.626070040*10**(-34)
da = 1.66 * 10**(-27) # da to kg
avogadro = 6.022 * 10**(23)
R = 8.314
patoatm = 101325

m1 = 16 * da # oxygen mass in kg

mu1 = 14579.0
mu2 = 36440.0
l = 4.398

atomic_mass_unit = 9.1093826 * 10**(-31) # kg
atomic_length_unit = 5.291772 * 10**(-11) # m

def simple_integrand(x):
    # x = [R, theta]
    R = x[0]
    theta = x[1]
    potential_value = potential(R, theta)

    if potential_value < 0:
        u_kt = potential_value * htoj / (k * temperature)
        return R**2 * sp.gammainc(2.5, -u_kt) * np.exp(-u_kt)
    else:
        return 0

def full_integrand(x):
    # x = [R, theta]
    R = x[0]
    theta = x[1]
    potential_value = potential(R, theta)

    if potential_value < 0:
        u_kt = potential_value * htoj / (k * temperature)
        radical_1 = np.sqrt(k * temperature / htoj / (1 / (2 * mu2 * R**2) + 1 / (mu1 * l**2)))
        radical_2 = np.sqrt(2 * np.sin(theta)**2 * k * temperature / htoj / (np.cos(theta)**2 / (mu2 * R**2) + 1 / (mu1 * l**2)))
        radical_3 = sp.gammainc(2.5, -u_kt)
        return R**2 * radical_1 * radical_2 * radical_3 * np.exp(-u_kt)
    else:
        return 0

def adv_integrand(x):
    # x = [R, theta]
    R = x[0]
    theta = x[1]
    potential_value = potential(R, theta)

    if potential_value < 0:
        u_kt = potential_value * htoj / (k * temperature)
        return R**2 * np.sin(theta) * np.exp(-u_kt) 
    else:
        return 0

co2_rot_summ = 8 * np.pi**2 * k * temperature / h**2 * m1 * (l * atomic_length_unit / 2)**2  

# simple form of hamtiltonian
integ = vegas.Integrator([[3, 20], [0, np.pi]])
integral = integ(simple_integrand, nitn = 20, neval = 10**4)
print('integral (atomic units): {0}'.format(integral))
# multiplying integral by rotational stat summ this value
# this value we compare with the next interim
interim = integral.mean * co2_rot_summ * atomic_length_unit**3
print('integral multiplied by rot stat summ: {0}'.format(interim))
# for constant
interim = integral.mean * atomic_length_unit**3

eq_const = 2 * np.pi * avogadro / (R * temperature) * interim * patoatm 
print('eq_const: {0}'.format(eq_const))

print('*'*30)

# more general form of hamiltonian
integ = vegas.Integrator([[3, 20], [0, np.pi]])
integral = integ(full_integrand, nitn = 20, neval = 10**4)
print('integral (atomic units): {0}'.format(integral.mean))
interim = integral.mean * atomic_length_unit**3
print('integral (CI) {0}'.format(interim))
# for constant

# translational summs cancel out, but rotational summ stays
eq_const = 2 * np.pi * avogadro / (R * temperature) / co2_rot_summ * interim * patoatm
print('eq_const: {0}'.format(eq_const))


# fullu advanced form of hamiltonian
integ = vegas.Integrator([[3, 20], [0, np.pi]])
integral = integ(adv_integrand, nitn = 20, neval = 10**4)
interim = integral.mean * atomic_length_unit ** 3


