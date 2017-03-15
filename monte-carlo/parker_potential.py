import numpy as np
import scipy.special as sp
import vegas
from time import time
from functools import partial

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

L = [sp.legendre(2 * k) for k in range(6)] 

def potential(R, theta):
	return sum([v(R, number = n) * L[n](np.cos(theta)) for n in range(6)])

k = 1.38064852 * 10**(-23) # J/k
htoj = 4.35974417 * 10**(-18) # hartree to Joules
avogadro = 6.022 * 10**(23)
length_unit = 5.291772 * 10**(-11)
pressure_coeff = 9.869 * 10**(-6) # between pascals and atmospheres
R = 8.314

def integrand(x, Temperature):
	# x = [R, theta]
	potential_value = potential(x[0], x[1]) * htoj
	if potential_value < 0:
		return sp.gammainc(2.5, - potential_value / (k * Temperature)) * np.exp(- potential_value / \
			(k * Temperature)) * x[0]**2
	else:
		return 0.0

def initialization(T):
	_integrand = partial(integrand, Temperature = T)

	integ = vegas.Integrator([[3., 20.], [0., np.pi]])
	result = integ(_integrand, nitn = 100, neval = 1000)
	print 'First integration. result = %s Q = %.2f' % (result, result.Q)

def cycle(T):
	_integrand = partial(integrand, Temperature = T)

	integ = vegas.Integrator([[3., 20.], [0., np.pi]])
	result = integ(_integrand, nitn = 10, neval = 10**4)
	print 'result = %s Q = %.2f' % (result, result.Q)
	constant = 4. * np.pi * avogadro / (R * T) * result.mean * length_unit**2 * pressure_coeff
	print 'Constant %.5f' % constant
	return constant

def save_constants(temperatures, constants):
	with open('data/parker.dat', mode = 'w') as out:
		for temperature, constant in zip(temperatures, constants):
			out.write(str(temperature) + ' ' + str(constant) + '\n')

initialization(T = 50)
temperatures = [50 + 5 * i for i in range(0, 61)] # K
constants = [cycle(temperature) for temperature in temperatures]

save_constants(temperatures, constants)