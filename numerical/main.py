import sys

from __particle__ import __particle__
from hamiltonian import Hamiltonian

from sympy import *
from sympy.physics.vector import dynamicsymbols
from sympy import oo

import numpy as np

m, r0, t = symbols('m r0 t')

q = dynamicsymbols('q')

freedom_degrees = [q]

particle_1 = __particle__(m = 1.2, x = -r0 * cos(q/2), y = 0., z = - r0 * sin(q/2))
particle_2 = __particle__(m = 1.2, x = -r0 * cos(q/2), y = 0., z = r0 * sin(q/2))
particle_3 = __particle__(m = oo, x = 0., y = 0., z = 0.)
particles = [particle_1, particle_2, particle_3]

J = 10.
r0_val = 1.0

def hamiltonian(q_val, p_val, theta, varphi):
	# print 'Cleaning numerical vectors of particles...'
	for particle in particles:
		particle.clear_numerical_vector()

	# print 'Adding numerical vectors for all particles...'
	subs = {r0: r0_val, q: q_val}
	for particle in particles:
		particle.evaluate(subs)

	subs = {r0: r0_val, q: q_val, 'theta': theta, 'varphi': varphi, 'J': J, 'p': p_val}
	# print 'Initializing hamiltonian object...'
	hamiltonian = Hamiltonian(particles = particles, freedom_degrees = freedom_degrees, subs = subs)

	inertia_tensor = hamiltonian.calculate_inertia_tensor()
	a_matrix = hamiltonian.calculate_a_matrix()
	A_matrix = hamiltonian.calculate_A_matrix()

	G11, G12, G21, G22 = hamiltonian.calculate_G_matrices(inertia_tensor = inertia_tensor, a_matrix = a_matrix, A_matrix = A_matrix)

	h = hamiltonian.calculate_hamiltonian(G_matrices = [G11, G12, G21, G22])
	return h

# first central difference approximation of f'(x)
def fcdiff(func, val, step = 0.0001, **kwargs):
	return (func(val + step) - func(val - step)) / (2 * step)

# second central difference approximation of f'(x)
def scdiff(func, val, step = 0.0001):
	return (func(val + 2 * step) - 4 * func(val + step) + 3 * func(val)) / (2 * step)

# rhs = []
# rhs.append(scdiff(hamiltonian, ))








