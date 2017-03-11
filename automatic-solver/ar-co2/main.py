import sys
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/lib') # path on ubuntu
sys.path.append('/Users/mac/repos/sympy_project/sympy/automatic-solver') # path on mac

from __particle__ import __particle__
from automatic_solver import AutomaticSolver
import autograd.numpy as np
import scipy.special as sp
from time import time

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

m1 = 29392. # mass of 16O in atomic units
m2 = 73480. # mass of 40Ar in atomic units
m3 = 22044. # mass of 12C in atomic units
M = 2 * m1 + m2 + m3
r0 = 4.398

particle_1 = __particle__(m = m1, __x__ = lambda R, theta: -m2/M * R - r0 / 2 * np.cos(theta), 
								 __y__ = lambda R, theta: r0/2 * np.sin(theta), 
								 __z__ = lambda R, theta: 0 * R,
								 vars = ['R', 'theta'])
particle_2 = __particle__(m = m2, __x__ = lambda R, theta: (2*m1 + m3) / M * R,
								 __y__ = lambda R, theta: 0 * R,
								 __z__ = lambda R, theta: 0 * R,
								 vars = ['R', 'theta'])
particle_3 = __particle__(m = m1, __x__ = lambda R, theta: -m2/M * R + r0 / 2 * np.cos(theta),
								  __y__ = lambda R, theta: - r0 / 2  * np.sin(theta),
								  __z__ = lambda R, theta: 0 * R,
								  vars = ['R', 'theta'])
particle_4 = __particle__(m = m3, __x__ = lambda R, theta: -m2 / M * R,
								  __y__ = lambda R, theta: 0 * R,
								  __z__ = lambda R, theta: 0 * R,
								  vars = ['R', 'theta'])

particles = [particle_1, particle_2, particle_3, particle_4]
__degrees__ = 2

J = 0.15
varphi0 = 0.5
psi0 = 0.850505507928

def potential(R, theta):
	return sum([v(R, number = n) * L[n](np.cos(theta)) for n in range(6)])

AS = AutomaticSolver(particles = particles, __degrees__ = __degrees__, potential = potential)
# setting angular momentum in AS!
AS.J = J

# q = [R, theta]
q = np.array([3., np.pi / 2])
p = np.array([-10., 0.01])

print AS.hamiltonian(q, p, psi0, varphi0) - potential(*q)
print potential(*q)