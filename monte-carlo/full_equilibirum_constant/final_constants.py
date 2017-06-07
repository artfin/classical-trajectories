from __future__ import print_function
import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
from potential_wrapper import potential

sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers')
from hutson import potdat4, extpot

# initializing Hutson potential
potdat4()

def hutpotential(R, theta):
    return extpot(R, np.cos(theta))

import numpy as np
import scipy.special as sp
import vegas
from functools import partial
from time import time

k = 1.38064852 * 10**(-23) # J / K
htoj = 4.35974417 * 10**(-18) # hartree to J
h = 6.626070040*10**(-34)
da = 1.66 * 10**(-27) # da to kg
avogadro = 6.022 * 10**(23)
R = 8.314
patoatm = 101325

atomic_mass_unit = 9.1093826 * 10**(-31) # kg
atomic_length_unit = 5.291772 * 10**(-11) # m

A_n1 = [.224247*10**2, .635744*10**2, .991128*10**2, .318652*10**3, .332826*10**3, .435837*10**3]
A_n2 = [-.716288*10**0, -.811806*10**0, -.117577*10, -.188135*10, -.214596*10, -.244616*10] 
A_n3 = [-.869136*10**(-1), -.753131*10**(-1), -.477138*10**(-1), 0.*10**0, 0.*10**0, 0.*10**0]
B_n1 = [.599056*10**0, .207391*10**0, .523497*10**(-1), .374994*10**(-1), .137376*10**(-1), .108283*10**0]
B_n2 = [-.650479*10**0, -.620877*10**0, -.735340*10**0, -.105547*10, -.103942*10, -.204765*10]
B_n3 = [-.320299*10**(-1), -.310717*10**(-1), -.296750*10**(-1), -.182219*10**(-1), -.494781*10**(-1), 0.*10**0]
C_6 = [114.5, 26.6, 0.]
C_8 =[2380.0, 2080.0, 410.0] 
Rn = [6.32925, 6.90026, 6.96450]

def exponent_1(R, number):
	return A_n1[number] * np.exp(A_n2[number] * R + A_n3[number] * R**2)

def exponent_2(R, number):
	return B_n1[number] * np.exp(B_n2[number] * R + B_n3[number] * R**2)

def vandervaals(R, number):
	return C_6[number] / R**6 + C_8[number] / R**8

def v(R, number):
	if number < 3:
		if R > Rn[number]: 
			return exponent_1(R, number = number) - vandervaals(R, number = number)
	
	return exponent_1(R, number = number) - exponent_2(R, number = number)

LEG = [sp.legendre(2 * s) for s in range(6)] 

def parpotential(R, theta):
    return sum([v(R, number = n) * LEG[n](np.cos(theta)) for n in range(6)])

def par_integrand(x, temperature):
    # x = [R, theta]
    R = x[0]
    theta = x[1]
    potential_value = parpotential(R, theta)
    
    if potential_value < 0:
        u_kt = potential_value * htoj / (k * temperature)
        res = R**2 * np.sin(theta) * sp.gammainc(2.5, -u_kt) * np.exp(-u_kt)
        return res
    else:
        return 0

def ab_integrand(x, temperature):
    # x = [R, theta]
    R = x[0]
    theta = x[1]
    potential_value = potential(R, theta)

    if potential_value < 0:
        u_kt = potential_value * htoj / (k * temperature)
        res = R**2 * np.sin(theta) * sp.gammainc(2.5, -u_kt) * np.exp(-u_kt)
        #print('potential_value: {0}; res: {1}'.format(potential_value, res))
        return res
    else:
        return 0

def hut_integrand(x, temperature):
    # x = [R, theta]
    R = x[0]
    theta = x[1]
    potential_value = hutpotential(R, theta)

    if potential_value < 0:
        u_kt = potential_value * htoj / (k * temperature)
        return R**2 * np.sin(theta) * sp.gammainc(2.5, -u_kt) * np.exp(-u_kt)
    else:
        return 0

def constant(temperature):
    _ = partial(par_integrand, temperature = temperature)

    print('*'*30)
    print('Temperature: {0}'.format(temperature))

    start = time()
    integ = vegas.Integrator([[5.0, 30], [0, np.pi]])
    integral = integ(_, nitn = 20, neval = 10**3)
    print('integral (atomic units): {0}'.format(integral))
    print('Time needed: {0}'.format(time() - start))

    eq_const = 2 * np.pi * avogadro / (R * temperature) * integral.mean * atomic_length_unit**3 * patoatm
    print('Constant: {0}'.format(eq_const))
    print('*'*30 + '\n')

    return eq_const

def save_constants(filename, temperatures, constants):
    with open('output/' + filename, mode = 'w') as out:
        for t, c in zip(temperatures, constants):
            out.write(str(t) + ' ' + str(c) + '\n')

temperatures = range(100, 500, 1)
ab_constants = [constant(t) for t in temperatures]
save_constants('final.parker.dat', temperatures, ab_constants)
