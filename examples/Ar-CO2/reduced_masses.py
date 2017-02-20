import sys
sys.path.append('/Users/mac/repos/sympy_project/sympy/') # path on mac
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/') # path on ubuntu

from lib.main import Particle, COM, Lagrange, Hamilton, LatexOutput

from sympy import *
from sympy.physics.vector import dynamicsymbols
from sympy import oo

mu1, mu2, l, t = symbols('mu1 mu2 l t')

q1, q2 = dynamicsymbols('q1 q2')
p1, p2 = dynamicsymbols('p1 p2')

# q1 == theta, q2 == R

particle1 = Particle(m = mu1, x = l * sin(q1), y = 0, z = l * cos(q1))
particle2 = Particle(m = mu2, x = 0, y = 0, z = q2)

particles = [particle1, particle2]

freedom_degrees = [q1, q2]
freedom_degrees_derivatives = [diff(degree, t) for degree in freedom_degrees]
conjugate_momentum = [p1, p2]

lagrange = Lagrange(particles = particles, freedom_degrees = freedom_degrees, freedom_degrees_derivatives = freedom_degrees_derivatives)
hamilton = Hamilton(lagrange = lagrange, freedom_degrees = freedom_degrees, conjugate_momentum = conjugate_momentum)
# latex_output = LatexOutput(lagrangian = lagrange.lagrangian, hamiltonian = hamilton.hamiltonian, name = 'Ar-CO2')
