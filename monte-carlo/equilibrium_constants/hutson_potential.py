import hutson
import numpy as np
import vegas
import scipy.special as sp
from functools import partial

from time import time

# loading parameters in potdat4 subroutine
hutson.potdat4()

k = 1.38064852 * 10**(-23) # J/k
htoj = 4.35974417 * 10**(-18) # hartree to Joules
avogadro = 6.022 * 10**(23)
length_unit = 5.291772 * 10**(-11)
R = 8.314
patoatm = 101325 

def integrand(x, Temperature):
	# x = [R, theta]
	potential_value = hutson.extpot(x[0], np.cos(x[1])) * htoj
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
	result = integ(_integrand, nitn = 50, neval = 10**4)
	print 'result = %s Q = %.2f' % (result, result.Q)
	constant = 2. * np.pi * avogadro / (R * T) * result.mean * patoatm * length_unit**3
	print 'Temperature: %d Constant %.5f' % (T, constant)
	return constant

def save_constants(temperatures, constants):
	with open('data/hutson.dat', mode = 'w') as out:
		for temperature, constant in zip(temperatures, constants):
			out.write(str(temperature) + ' ' + str(constant) + '\n')

#initialization(T = 50)
temperatures = [100 + 5 * i for i in range(0, 61)] # K
constants = [cycle(temperature) for temperature in temperatures]

#save_constants(temperatures, constants)


