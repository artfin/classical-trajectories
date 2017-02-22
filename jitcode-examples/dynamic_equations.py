from sympy import *
from sympy.physics.vector import dynamicsymbols

import numpy as np
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt

from pprint import pprint 

class DynamicEquations(object):
	def __init__(self, hamiltonian, angular_momentum, angles, freedom_degrees, conjugate_momentum):
		theta, varphi = symbols('theta varphi')

		self.rhs = []
		self.rhs.append((diff(hamiltonian, angular_momentum[0]) * cos(angles[0]) + 
				    diff(hamiltonian, angular_momentum[1]) * sin(angles[0])) / tan(angles[1]) -
				 	diff(hamiltonian, angular_momentum[2]))
		self.rhs.append( diff(hamiltonian, angular_momentum[0]) * sin(angles[0]) - 
					diff(hamiltonian, angular_momentum[1]) * cos(angles[0]) )
		self.rhs.append( diff(hamiltonian, conjugate_momentum[0]))
		self.rhs.append(-diff(hamiltonian, freedom_degrees[0]))

		J = 10
		Vm = 624 
		Vp = 224

		substitutions = {'m': 1,
						 'r0': 2,
						 'Jx': J * cos(varphi) * sin(theta),
						 'Jy': J * sin(varphi) * sin(theta),
						 'Jz': J * cos(theta),
						 'Vm': Vm,
						 'Vp': Vp}

		self.rhs_w_angles = []	
		for equation in self.rhs:
			self.rhs_w_angles.append(equation.subs(substitutions))


		init = [.555, 0, 1.766, 14.8]
		t = np.linspace(0, 500, 10)
		sol = odeint(self.g, init, t)

		plt.plot(t, sol[:,2], 'b', label = 'q(t)')
		plt.plot(t, sol[:,3], 'g', label = 'p(t)')
		plt.legend(loc = 'best')
		plt.xlabel('t')
		plt.grid()
		plt.show()

	def g(self, y, t):
		print 'here you are: {0}'.format(t)
		theta, varphi, q, p = y
		res = []
		vals = {'theta': theta, 'varphi': varphi, 'q': q, 'p': p}
		return [eq.evalf(subs = vals) for eq in self.rhs_w_angles]

if __name__ == '__main__':
	Jx, Jy, Jz = symbols('Jx Jy Jz')
	q, p = symbols('q p')
	phi, theta = symbols('varphi theta')
	m, r0 = symbols('m r0')
	Vm, Vp = symbols('Vm Vp')

	angles = [phi, theta]
	angular_momentum = [Jx, Jy, Jz]
	freedom_degrees = [q]
	conjugate_momentum = [p]
	hamiltonian = Jx**2 / (2*m*r0**2 * (1 - cos(q))) + Jy**2 / (2*m*r0**2) + Jz**2 / (2*m*r0**2 * (1 + cos(q))) + p**2 / (m * r0**2) + \
				Vm / (2*m*r0**2*(1 - cos(q))) + Vp / (2*m*r0**2*(1 + cos(q)))
	
	dynamicEquations = DynamicEquations(hamiltonian = hamiltonian, angular_momentum = angular_momentum, angles = angles,
		freedom_degrees = freedom_degrees, conjugate_momentum = conjugate_momentum)

	#dynamicEquations.g(vars = [.1, .2, 3, 1])

