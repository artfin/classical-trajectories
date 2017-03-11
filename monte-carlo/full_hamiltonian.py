from sympy import symbols
from sympy import sin, cos, lambdify
import numpy as np
import scipy.special as sp
import vegas
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

Jx, Jy, Jz = symbols('Jx Jy Jz')
R, theta = symbols('R theta')
mu1, mu2 = symbols('mu1 mu2')
pR, pT = symbols('p_R p_theta')
l = symbols('l')

psi, varphi, alpha = symbols('psi varphi alpha')
p_psi, p_varphi, p_alpha = symbols('p_psi p_varphi p_alpha')

kinetic_energy = Jx * (Jx / (mu2 * R**2) + Jz * cos(theta) / (sin(theta) * mu2 * R**2)) / 2 + \
				 Jy**2 / (2 * mu2 * R**2) + \
				 Jz * (Jx * cos(theta) / (sin(theta) * mu2 * R**2) + \
				 	   Jz * (mu1 * l**2 * cos(theta)**2 + mu2 * R**2) / (mu1 * l**2 * sin(theta)**2 * mu2 * R**2))/2 + \
				 pT**2 * (mu1 * l**2 + mu2 * R**2) / (2 * mu2 * R**2 * mu1 * l**2) + \
				 pR**2 / (2 * mu2) - pT * Jy / (mu2 * R**2)

_Jx = (sin(psi) * p_varphi + cos(psi) * p_alpha * sin(alpha) - cos(alpha) * sin(psi) * p_psi) / sin(alpha)
_Jy = - (-cos(psi) * p_varphi + sin(psi) * p_alpha * sin(alpha) + cos(alpha) * cos(psi) * p_psi) / sin(alpha)
_Jz = p_psi

_jx_func = lambdify((psi, varphi, alpha, p_psi, p_varphi, p_alpha), _Jx)
_jy_func = lambdify((psi, varphi, alpha, p_psi, p_varphi, p_alpha), _Jy)
_jz_func = lambdify((psi, varphi, alpha, p_psi, p_varphi, p_alpha), _Jz)

# subs euler angles
kinetic_energy = kinetic_energy.subs({'Jx': _Jx, 'Jy': _Jy, 'Jz': _Jz})

# subs constants
kinetic_energy = kinetic_energy.subs({'mu1': 14579., 'mu2': 36440., 'l': 4.398})

_hamiltonian = lambdify((psi, varphi, alpha, p_psi, p_varphi, p_alpha, R, theta, pR, pT), kinetic_energy)
hamiltonian = lambda psi, varphi, alpha, p_psi, p_varphi, p_alpha, R, theta, pR, pT: \
		_hamiltonian(psi, varphi, alpha, p_psi, p_varphi, p_alpha, R, theta, pR, pT) + potential(R, theta)
# print hamiltonian(1., 1., 1., 0.1, 0.1, 0.1, 3., np.pi/2, -10., 0.01)

k = 1.38064852 * 10**(-23) # J/K
htoj = 4.35974417 * 10**(-18) # hartree to Joules
avogadro = 1. #6.022 * 10**(23)
R = 8.314

def integrand(x, Temperature):
	# x = [psi, varphi, alpha, p_psi, p_varphi, p_alpha, R, theta, pR, pT]
	hamiltonian_value = hamiltonian(*x)
	# print 'hamiltonian value: ' + str(hamiltonian_value)

	if hamiltonian_value < 0:
		try:
			res = np.exp(- hamiltonian_value * htoj / (k * Temperature))
			#print 'type res: {0}; res: {1}'.format(type(res), res)
			return res
		except OverflowError:
			return 0.0
	else:
		return 0.0

limits = [[0., 2 * np.pi], # psi
		  [0., 2 * np.pi], # varphi
		  [0., np.pi], # alpha / theta
		  [-10., 10.], # p_psi
		  [-10., 10.], # p_varphi
		  [-10., 10.], # p_alpha / p_thea
		  [3., 20.], # R 
		  [0., np.pi], # theta 
		  [-100., 100.], # pR,
		  [-100., 100.], # pT
		  ]

length_unit = 5.291772 * 10**(-11)
amu = 1.66053904 * 10**(-27) # amu to kg
complex_mass = (40 + 12 + 32) * amu
co2_mass = (12 + 32) * amu
ar_mass = 40 * amu
m1 = 16 * amu
r0 = 116.3 * 10**(-12) # c=0 distance in m
h = 6.626070040*10**(-34)
pressure_coeff = 9.869 * 10**(-6) # between pascals and atmospheres

def initialization(T):
	_integrand = partial(integrand, Temperature = T)

	integ = vegas.Integrator(limits)
	result = integ(_integrand, nitn = 10, neval = 10**6)
	print result.summary()
	print 'First integration. result = %s Q = %.2f' % (result, result.Q)


def pre_constant(temperature):
	return (84./40.)**1.5 / (4. * np.pi**2 * k * temperature / h**2 * m1 * r0**2 * \
		   (2.*np.pi*co2_mass * k *temperature / h**2)**(1.5))

def cycle(T):
	_integrand = partial(integrand, Temperature = T)

	integ = vegas.Integrator(limits)
	result = integ(_integrand, nitn = 10, neval = 10**4)
	print 'result = %s Q = %.2f' % (result, result.Q)
	print 'preconstant: {0}'.format(pre_constant(T))
	constant = pre_constant(T) / (R * T) * result.mean * pressure_coeff
	print 'Constant %.5f' % constant
	return constant

initialization(T = 250)

cycle(T = 250)



