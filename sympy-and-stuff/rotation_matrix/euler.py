from __future__ import print_function

from sympy import symbols, solve
from sympy import sin, cos
from sympy.solvers import solve

phi, theta, psi = symbols('phi theta psi')

jx, jy, jz = 1, 2, 4
Jx, Jy, Jz = 1, 4, 2

eq1 = ((cos(psi) * cos(phi) - cos(theta) * sin(phi) * sin(psi)) * jx + (cos(psi) * sin(phi) + cos(theta) * cos(phi) * sin(psi)) * jy + sin(psi) * sin(theta) * jz - Jx).expand(trig = True)
eq2 = ((-sin(psi) * cos(phi) - cos(theta) * sin(phi) * cos(psi)) * jx + (-sin(psi) * sin(phi) + cos(theta) * cos(phi) * cos(psi)) * jy + cos(psi) * sin(theta) * jz - Jy).expand(trig = True)
eq3 = (sin(theta) * sin(phi) * jx - sin(theta) * cos(phi) * jy + cos(theta) * jz - Jz).expand(trig = True) 

print('Solving system...')
sol = solve([eq1, eq2, eq3], [phi, theta, psi])
print(sol)



