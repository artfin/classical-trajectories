import sys
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/lib')

from __particle__ import __particle__
from sympy import *
from sympy.physics.vector import dynamicsymbols
from sympy import oo

import autograd.numpy as np
from autograd import elementwise_grad

m, r0, t = symbols('m r0 t')

q = dynamicsymbols('q')
p = dynamicsymbols('p')

particle_1 = __particle__(m = m, x = -r0 * cos(q/2), y = S(0), z = - r0 * sin(q/2))
particle_2 = __particle__(m = m, x = -r0 * cos(1/2), y = S(0), z = r0 * sin(q/2))
particle_3 = __particle__(m = oo, x = 0, y = 0, z = 0)
particles = [particle_1, particle_2, particle_3]

subs = {r0: 1., q: 2.0}

# for particle in particles:
# 	print particle.diff(var = q, subs = subs)

def R(_q):
	subs = {r0: 1, q: _q}
	return [element.evalf(subs = subs) for element in particle_1.__vec__()]

print R(_q = 1)

d_fun = elementwise_grad(R)

values = range(10)
for value in values:
	print 'value: {0}; d_fun: {1}'.format(value, d_fun(value))









