import sys
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/lib')
sys.path.append('/Users/mac/repos/sympy_project/sympy/lib')

from __particle__ import __particle__
from hamiltonian import Hamiltonian

import sympy as sp
from sympy import lambdify, symbols, sympify
from sympy.physics.vector import dynamicsymbols
from sympy import oo

import autograd.numpy as np
from autograd.core import primitive
from autograd import elementwise_grad, grad

m, r0, t = symbols('m r0 t')

q = dynamicsymbols('q')

freedom_degrees = [q]

particle_1 = __particle__(m = 1.2, x = -r0 * sp.cos(q/2), y = 0., z = - r0 * sp.sin(q/2))
particle_2 = __particle__(m = 1.2, x = -r0 * sp.cos(q/2), y = 0., z = r0 * sp.sin(q/2))
particle_3 = __particle__(m = oo, x = 0., y = 0., z = 0.)
particles = [particle_1, particle_2, particle_3]

J = 10.
p = 1.0
r0_val = 1.0
theta = 0.1
varphi = 0.5

g = lambdify([r0, q], particle_1.__vec__, modules = "numpy")
g_grad = lambdify([r0, q], particle_1.diff(q), modules = "numpy")

@primitive
def __vector__(q):
	return g(2., q)

def __vectorgrad__(x, ans):
	return g_grad(2., x) * ans

__vector__.defvjp(__vectorgrad__)

print __vector__(1.)


# def function(q_val):
# 	print 'type of q_val: {0}'.format(type(q_val))
# 	print 'Cleaning numerical vectors of particles...'
# 	for particle in particles:
# 		particle.clear_numerical_vector()

# 	return g(1, q_val)
# 	# print 'Adding numerical vectors for all particles...'
# 	# subs = {r0: r0_val, q: q_val}
# 	# for particle in particles:
# 	# 	particle.evaluate(subs)

# 	# subs = {r0: r0_val, q: q_val, 'theta': theta, 'varphi': varphi, 'J': J, 'p': p}
# 	# print 'Initializing hamiltonian object...'
# 	# hamiltonian = Hamiltonian(particles = particles, freedom_degrees = freedom_degrees, subs = subs)

# 	# inertia_tensor = hamiltonian.calculate_inertia_tensor()
# 	# a_matrix = hamiltonian.calculate_a_matrix()
# 	# A_matrix = hamiltonian.calculate_A_matrix()

# 	# G11, G12, G21, G22 = hamiltonian.calculate_G_matrices(inertia_tensor = inertia_tensor, a_matrix = a_matrix, A_matrix = A_matrix)

# 	# H = hamiltonian.calculate_hamiltonian(G_matrices = [G11, G12, G21, G22])
# 	# return H

# print function(1.)
# grad_function = elementwise_grad(function)
# print grad_function(1.)












