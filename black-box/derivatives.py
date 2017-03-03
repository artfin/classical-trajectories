import sympy as sp
import numpy as np

from itertools import product
from pprint import pprint

J, theta, varphi = sp.symbols('J theta varphi')

M = sp.Matrix([
	[J * sp.cos(varphi) * sp.cos(theta), J * sp.sin(varphi) * sp.cos(theta), - J * sp.sin(theta)],
	[-J * sp.sin(varphi) * sp.sin(theta), J * sp.cos(varphi) * sp.sin(theta), 0],
	[sp.cos(varphi) * sp.sin(theta), sp.sin(varphi) * sp.sin(theta), sp.cos(theta)]
	])

pprint(M)

M_inv = M.inv()

M_inv = sp.Matrix(map(sp.trigsimp, M_inv)).reshape(3,3)
pprint(M_inv)

test = M * M_inv
test = sp.Matrix(map(sp.trigsimp, test)).reshape(3, 3)
pprint(test)



