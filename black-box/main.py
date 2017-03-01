import sys
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/lib') # path on ubuntu
sys.path.append('/Users/mac/repos/sympy_project/sympy/lib') # path on mac

from pprint import pprint
from itertools import product

from __particle__ import __particle__

import autograd.numpy as np
from autograd import grad, grad_named

# from scipy.integrate import odeint
import scipy.integrate as spi

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def extract_column(data, column):
	res = []
	for row in data:
		res.append(row[column])
	return res

# np.dot() instead of np.array().dot(np.array()) !

h = 1.
n = 0.
print 'Deformational Energy level: {0}'.format(n)

q0 = 1.82387
omega0 = 0.00663829
I0 = 6021.0735
 
r0 = 1.810357468
m = I0 / r0**2

J = 25.

####################################################################

particle_1 = __particle__(m = m, __x__ = lambda q: -r0 * np.cos(q/2), 
								 __y__ = lambda q: 0*q, 
								 __z__ = lambda q: -r0 * np.sin(q/2))
particle_2 = __particle__(m = m, __x__ = lambda q: -r0 * np.cos(q/2),
								 __y__ = lambda q: 0*q,
								 __z__ = lambda q:  r0 * np.sin(q/2))
particles = [particle_1, particle_2]

def hamiltonian(q = None, p = None, jx = None, jy = None, jz = None, effective_potential = False):
	"""
	q -- np.array object with all degrees of freedom
	p -- np.array object with conjugate momentum
	"""

	J_vector = np.array([jx, jy, jz])
	p_vector = np.array([p])

	Ixx = sum([particle.m * (particle.__y__(_deg)**2 + particle.__z__(_deg)**2) for _deg in q for particle in particles])
	Iyy = sum([particle.m * (particle.__x__(_deg)**2 + particle.__z__(_deg)**2) for _deg in q for particle in particles])
	Izz = sum([particle.m * (particle.__x__(_deg)**2 + particle.__y__(_deg)**2) for _deg in q for particle in particles])
	Ixy = - sum([particle.m * particle.__x__(_deg) * particle.__y__(_deg) for _deg in q for particle in particles])
	Ixz = - sum([particle.m * particle.__x__(_deg) * particle.__z__(_deg) for _deg in q for particle in particles])
	Iyz = - sum([particle.m * particle.__y__(_deg) * particle.__z__(_deg) for _deg in q for particle in particles])
	
	inertia_tensor = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

	pprint(inertia_tensor)

	# a = np.zeros((q.shape[0], q.shape[0]))
	# for j, k in product(range(q.shape[0]), range(q.shape[0])):
	# 	a[j, k] = sum([particle.m * (particle.__dx__(q[j]) * particle.__dx__(q[k]) + \
	# 							 	 particle.__dy__(q[j]) * particle.__dy__(q[k]) + \
	# 								 particle.__dz__(q[j]) * particle.__dz__(q[k])) for particle in particles])

	a = np.array([sum([particle.m * (particle.__dx__(q)**2 + particle.__dy__(q)**2 + particle.__dz__(q)**2) for particle in particles])])

	pprint(a)

	G11 = np.linalg.inv(inertia_tensor)
	G22 = 1/a
	G12 = 0
	G21 = 0

	angular_component = 0.5 * np.dot(np.dot(J_vector, G11), J_vector)
	potential = Vm / (2 * I0 * (1 - np.cos(q))) + Vp / (2 * I0 * (1 + np.cos(q)))
	
	if not effective_potential:
		kinetic_component = 0.5 * p_vector * G22 * p_vector 
		return angular_component + kinetic_component + potential
	else:
		return angular_component + potential


######################################################################

Vp = 1./4 * I0**2 * omega0**2 * (1 + np.cos(q0))**2 
Vm = 1./4 * I0**2 * omega0**2 * (1 - np.cos(q0))**2
print 'Vm: {0}; Vp: {1}'.format(Vm, Vp)

varphi0 = 0.01
theta0 = 0.15
Jx0 = J * np.cos(varphi0) * np.sin(theta0)
Jy0 = J * np.sin(varphi0) * np.sin(theta0)
Jz0 = J * np.cos(theta0)
print 'Jx0: {0}; Jy0: {1}; Jz0: {2}'.format(Jx0, Jy0, Jz0)

E = h**2 * (n + (np.sqrt(Jz0**2 + Vp) + np.sqrt(Jx0**2 + Vm))/(2 * h)) * \
		   (n + 1 + (np.sqrt(Jz0**2 + Vp) + np.sqrt(Jx0**2 + Vm))/(2 * h)) / I0
print 'Energy: {0}'.format(E)

qe = 1.503583924

effective_potential = hamiltonian(q = np.array([qe]), jx = Jx0, jy = Jy0, jz = Jz0, effective_potential = True)
print 'effective_potential: {0}'.format(effective_potential)

pini = np.sqrt(I0 * (E - effective_potential))
print 'pini: {0}'.format(pini)

print 'hamiltonian: {0}'.format(hamiltonian(q = np.array([qe]), p = pini, jx = Jx0, jy = Jy0, jz = Jz0))

dham_dq = grad_named(hamiltonian, argname = 'q')
dham_dp = grad_named(hamiltonian, argname = 'p')
dham_djx = grad_named(hamiltonian, argname = 'jx')
dham_djy = grad_named(hamiltonian, argname = 'jy')
dham_djz = grad_named(hamiltonian, argname = 'jz')

def fcdiff(func, vals, argnum, step = 0.00001):
	_vals_p, _vals_m = vals[:], vals[:]
	_vals_p[argnum] = _vals_p[argnum] + step
	_vals_m[argnum] = _vals_m[argnum] - step
	return (func(*_vals_p) - func(*_vals_m)) / (2 * step)

print 'automatic dH/dq: {0}'.format(dham_dq(np.array([qe]), pini, Jx0, Jy0, Jz0))
print 'numeric dH/dq: {0}'.format(fcdiff(func = hamiltonian, vals = [np.array([qe]), pini, Jx0, Jy0, Jz0], argnum = 0))

print 'automatic dH/dp: {0}'.format(dham_dp(qe, pini, Jx0, Jy0, Jz0))
print 'numeric dH/dp: {0}'.format(fcdiff(func = hamiltonian, vals = [qe, pini, Jx0, Jy0, Jz0], argnum = 1))

print 'automatic dH/djx: {0}'.format(dham_djx(qe, pini, Jx0, Jy0, Jz0))
print 'numeric dH/djx: {0}'.format(fcdiff(func = hamiltonian, vals = [qe, pini, Jx0, Jy0, Jz0], argnum = 2))

print 'automatic dH/djy: {0}'.format(dham_djy(qe, pini, Jx0, Jy0, Jz0))
print 'numeric dH/djy: {0}'.format(fcdiff(func = hamiltonian, vals = [qe, pini, Jx0, Jy0, Jz0], argnum = 3))

print 'automatic dH/djz: {0}'.format(dham_djz(qe, pini, Jx0, Jy0, Jz0))
print 'numeric dH/djz: {0}'.format(fcdiff(func = hamiltonian, vals = [qe, pini, Jx0, Jy0, Jz0], argnum = 4))

rhs = []

rhs = [lambda q, p, jx, jy, jz:   dham_dp(q, p, jx, jy, jz),
	   lambda q, p, jx, jy, jz: - dham_dq(q, p, jx, jy, jz),
	   lambda q, p, jx, jy, jz:   dham_djz(q, p, jx, jy, jz) * jy - dham_djy(q, p, jx, jy, jz) * jz,
	   lambda q, p, jx, jy, jz:   dham_djx(q, p, jx, jy, jz) * jz - dham_djz(q, p, jx, jy, jz) * jx,
	   lambda q, p, jx, jy, jz:   dham_djy(q, p, jx, jy, jz) * jx - dham_djx(q, p, jx, jy, jz) * jy
	   ]

# rhs should return derivatives in the same order: q_dot, p_dot, theta_dot, varphi_dot
def derivatives(t, y):
	#print t, y
	print t, y[3] 
	return [eq(y[0], y[1], y[2], y[3], y[4]) for eq in rhs]

init = [qe, pini, Jx0, Jy0, Jz0]

t_start = 0.
t_end = 100.
t_step = 1.

ode = spi.ode(derivatives)

# ode.set_integrator('lsoda', nsteps = 500, method = 'bdf', atol = 1e-6)
# ode.set_initial_value(init, t_start)

# sol = []
# t = []
# while ode.successful() and ode.t < t_end:
# 	ode.integrate(ode.t + t_step)
# 	t.append(ode.t)
# 	sol.append(ode.y)


# q_list = extract_column(data = sol, column = 0)
# p_list = extract_column(data = sol, column = 1)
# jx_list = extract_column(data = sol, column = 2)
# jy_list = extract_column(data = sol, column = 3)
# jz_list = extract_column(data = sol, column = 4)

# saving it just in case
# np.savetxt("simple-water.dat", jy_list)

# plt.plot(t, q_list, 'b', label = 'q(t)')
# plt.plot(t, p_list, 'g', label = 'p(t)')

# plt.plot(t, jx_list, 'g', label = 'jx(t)')
# plt.plot(t, jy_list, 'b', label = 'jy(t)')
# plt.plot(t, jz_list, 'r', label = 'jz(t)')
# plt.legend(loc = 'best')
# plt.xlabel('t')
# plt.grid()
# plt.show()

# phi = np.linspace(0, 2 * np.pi, 100)
# theta = np.linspace(0, np.pi, 100)
# xm = r * np.outer(np.cos(phi), np.sin(theta))
# ym = r * np.outer(np.sin(phi), np.sin(theta))
# zm = r * np.outer(np.ones(np.size(phi)), np.cos(theta))

# xx = J * np.cos(varphi_list) * np.sin(theta_list)
# yy = J * np.sin(varphi_list) * np.sin(theta_list)
# zz = J * np.cos(theta_list)

# fig = plt.figure()
# ax = plt.axes(projection='3d')

# ax.plot(jx_list, jy_list, jz_list, '-b')
# plt.tight_layout()
# plt.show()













