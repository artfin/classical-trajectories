import sys
sys.path.append('/Users/mac/repos/sympy_project/sympy/') # path on mac
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/') # path on ubuntu

from lib.main import Particle, COM, Lagrange, Hamilton, LatexOutput

from sympy import *
from sympy.physics.vector import dynamicsymbols
from sympy import oo

m1, m2, r0, t = symbols('m1 m2 r0 t')

particle1 = Particle(m = m1, x = -r0, y = 0, z = 0)
particle2 = Particle(m = m1, x =  r0, y = 0, z = 0)
particle3 = Particle(m = m2, x = 0,   y = 0, z = 0)

particles = [particle1, particle2, particle3]
com = COM(particles)
print com
particles = com.particles

freedom_degrees = []
freedom_degrees_derivatives = [diff(degree, t) for degree in freedom_degrees]
conjugate_momentum = []

lagrange = Lagrange(particles = particles, freedom_degrees = freedom_degrees, freedom_degrees_derivatives = freedom_degrees_derivatives)
hamilton = Hamilton(lagrange = lagrange, freedom_degrees = freedom_degrees, conjugate_momentum = conjugate_momentum)
# latex_output = LatexOutput(lagrangian = lagrange.lagrangian, hamiltonian = hamilton.hamiltonian, name = 'Ar-CO2')
