import sys
sys.path.append('/Users/mac/repos/sympy_project/sympy/') # path on mac
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/') # path on ubuntu

from lib.main import Particle, COM, Lagrange, Hamilton, LatexOutput

from sympy import *
from sympy.physics.vector import dynamicsymbols
from sympy import oo

mu1, mu2, mu3, l1, l2, t = symbols('mu1 mu2 mu3 l1 l2 t')

q1, q2, q3, q4 = dynamicsymbols('q1 q2 q3 q4')
p1, p2, p3, p4 = dynamicsymbols('p1 p2 p3 p4')

# q1 == theta1, q2 == phi, q3 == theta2, q4 == R

particle1 = Particle(m = mu1, x = l1 * sin(q1), y = 0, z = l1 * cos(q1))
particle2 = Particle(m = mu2, x = l2 * cos(q2) * sin(q3), y = l2 * sin(q2) * sin(q3), z = l2 * cos(q3))
particle3 = Particle(m = mu3, x = 0, y = 0, z = q4)

particles = [particle1, particle2, particle3]

freedom_degrees = [q1, q2, q3, q4]
freedom_degrees_derivatives = [diff(degree, t) for degree in freedom_degrees]
conjugate_momentum = [p1, p2, p3, p4]

lagrange = Lagrange(particles = particles, freedom_degrees = freedom_degrees, freedom_degrees_derivatives = freedom_degrees_derivatives)
hamilton = Hamilton(lagrange = lagrange, freedom_degrees = freedom_degrees, conjugate_momentum = conjugate_momentum)
latex_output = LatexOutput(lagrangian = lagrange.lagrangian, hamiltonian = hamilton.hamiltonian, name = 'N2N2')
