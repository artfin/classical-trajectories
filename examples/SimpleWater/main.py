import sys
sys.path.append('/Users/mac/repos/sympy_project/sympy/')
from lib.main import Particle, COM, Lagrange, Hamilton, LatexOutput

from sympy import *
from sympy.physics.vector import dynamicsymbols
from sympy import oo

m, r0, t = symbols('m r0 t')

q = dynamicsymbols('q')
p = dynamicsymbols('p')

particle1 = Particle(m = m, x = -r0 * cos(q/2), y = 0, z = - r0 * sin(q/2))
particle2 = Particle(m = m, x = -r0 * cos(q/2), y = 0, z = r0 * sin(q/2))
particle3 = Particle(m = oo, x = 0, y = 0, z = 0)

particles = [particle1, particle2, particle3]
com = COM(particles)
particles = com.particles

freedom_degrees = [q]
freedom_degrees_derivatives = [diff(degree, t) for degree in freedom_degrees]
conjugate_momentum = [p]


lagrange = Lagrange(particles = particles, freedom_degrees = freedom_degrees, freedom_degrees_derivatives = freedom_degrees_derivatives)
hamilton = Hamilton(lagrange = lagrange, conjugate_momentum = conjugate_momentum)
latex_output = LatexOutput(lagrangian = lagrange.lagrangian, hamiltonian = hamilton.hamiltonian, name = 'SimpleWater')