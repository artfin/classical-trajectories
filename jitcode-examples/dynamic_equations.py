from sympy import *
from sympy.physics.vector import dynamicsymbols

from pprint import pprint 

class DynamicEquations(object):
	def __init__(self, hamiltonian, angular_momentum, angles, freedom_degrees, conjugate_momentum):
		self.rhs = []
		self.rhs.append((diff(hamiltonian, angular_momentum[0]) * cos(angles[0]) + 
				    diff(hamiltonian, angular_momentum[1]) * sin(angles[0])) / tan(angles[1]) -
				 	diff(hamiltonian, angular_momentum[2]))
		self.rhs.append( diff(hamiltonian, angular_momentum[0]) * sin(angles[0]) - 
					diff(hamiltonian, angular_momentum[1]) * cos(angles[0]) )
		self.rhs.append( diff(hamiltonian, conjugate_momentum[0]))
		self.rhs.append(-diff(hamiltonian, freedom_degrees[0]))

if __name__ == '__main__':
	Jx, Jy, Jz = symbols('Jx Jy Jz')
	q, p = dynamicsymbols('q p')
	phi, theta = dynamicsymbols('varphi theta')
	m, r0 = symbols('m r0')

	angles = [phi, theta]
	angular_momentum = [Jx, Jy, Jz]
	freedom_degrees = [q]
	conjugate_momentum = [p]
	hamiltonian = Jx**2 / (2*m*r0**2 * (1 - cos(q))) + Jy**2 / (2*m*r0**2) + Jz**2 / (2*m*r0**2 * (1 + cos(q))) + p**2 / (m * r0**2)
	
	dynamicEquations = DynamicEquations(hamiltonian = hamiltonian, angular_momentum = angular_momentum, angles = angles,
		freedom_degrees = freedom_degrees, conjugate_momentum = conjugate_momentum)

	pprint(dynamicEquations.rhs)
