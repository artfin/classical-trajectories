import sys
sys.path.append('/Users/mac/repos/sympy_project/sympy/')  # path on mac
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/') # path on ubuntu

from lib.main import Particle, COM, Lagrange, Hamilton, LatexOutput

from sympy import *
from sympy.physics.vector import dynamicsymbols
from sympy import oo

m, M, t = symbols('m M t')

q, r1, r2 = dynamicsymbols('q r1 r2')
p, p1, p2 = dynamicsymbols('p p1 p2')

particle1 = Particle(m = m, x = -r1 * cos(q * Rational(1, 2)), y = 0, z = - r1 * sin(q * Rational(1, 2)))
particle2 = Particle(m = m, x = -r2 * cos(q * Rational(1, 2)), y = 0, z = r2 * sin(q * Rational(1, 2)))
particle3 = Particle(m = M, x = 0, y = 0, z = 0)
particles = [particle1, particle2, particle3]

# list of freedom degrees
freedom_degrees = [q, r1, r2]
freedom_degrees_derivatives = [diff(degree, t) for degree in freedom_degrees]
conjugate_momentum = [p, p1, p2]

com = COM(particles, freedom_degrees = freedom_degrees)
print com
particles = com.particles

lagrange = Lagrange(particles = particles, freedom_degrees = freedom_degrees, freedom_degrees_derivatives = freedom_degrees_derivatives)
hamilton = Hamilton(lagrange = lagrange, freedom_degrees = freedom_degrees, conjugate_momentum = conjugate_momentum)
# latex_output = LatexOutput(lagrangian = lagrange.lagrangian, hamiltonian = hamilton.hamiltonian, name = 'FullWater')