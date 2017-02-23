import sys
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/lib')

from __particle__ import __particle__
from sympy import *
from sympy.physics.vector import dynamicsymbols
from sympy import oo

import autograd.numpy as np
from autograd import elementwise_grad, grad

m, r0, t = symbols('m r0 t')

q = dynamicsymbols('q')

particle_1 = __particle__(m = m, x = -r0 * cos(q/2), y = 0., z = - r0 * sin(q/2))
particle_2 = __particle__(m = m, x = -r0 * cos(q/2), y = 0., z = r0 * sin(q/2))
particle_3 = __particle__(m = oo, x = 0., y = 0., z = 0.)
particles = [particle_1, particle_2, particle_3]

# inserting r0 = 1 for all particles
subs = {r0: 1.}
for particle in particles:
	particle.evaluate(subs)

# creating a function out of particle_1 vector
f0 = lambdify(q, particle_1.__vec__()[0], modules = np)
f1 = lambdify(q, particle_1.__vec__()[1], modules = np)
f2 = lambdify(q, particle_1.__vec__()[2], modules = np)

# automatically differentiate it
grad_f0 = elementwise_grad(f0)
grad_f1 = elementwise_grad(f1)
grad_f2 = elementwise_grad(f2)
print [grad_f0(10.), grad_f1(10.), grad_f2(10.)]

# check the answer with symbolic differentiation
subs = {q: 10.0}
print particle_1.diff(var = q, subs = subs)



