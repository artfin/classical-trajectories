import sys
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/lib') # path on ubuntu
sys.path.append('/Users/mac/repos/sympy_project/sympy/lib') # path on mac
sys.path.append('/Users/mac/repos/sympy_project/sympy/black-box') # path on mac

import __builtin__

from __particle__ import __particle__

import autograd.numpy as np
from itertools import chain

from autograd import grad, grad_named

import scipy.integrate as spi

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from time import time

I0 = 6021.0735
r0 = 1.810357468
m = I0 / r0**2

particle_1 = __particle__(m = m, __x__ = lambda q, r1, r2: -r1 * np.cos(q/2), 
								 __y__ = lambda q, r1, r2: 0 * q, 
								 __z__ = lambda q, r1, r2: -r1 * np.sin(q/2),
								 vars = ['q', 'r1', 'r2'])
particle_2 = __particle__(m = m, __x__ = lambda q, r1, r2: -r2 * np.cos(q/2),
								 __y__ = lambda q, r1, r2: 0 * q,
								 __z__ = lambda q, r1, r2:  r2 * np.sin(q/2),
								 vars = ['q', 'r1', 'r2'])
particles = [particle_1, particle_2]

__degrees__ = 3

q0 = 1.82387
omega0 = 0.00663829
Vp = 1./4 * I0**2 * omega0**2 * (1 + np.cos(q0))**2 
Vm = 1./4 * I0**2 * omega0**2 * (1 - np.cos(q0))**2

J = 25.

# adding variables to builtin scope, so that they could be accessed in hamiltonian function
__builtin__.particles = particles
__builtin__.Vm = Vm
__builtin__.Vp = Vp
__builtin__.I0 = I0
__builtin__.J = J
from black_box import hamiltonian, extract_column

h = 1.
n = 0.
n1 = 0.
n2 = 0.
print 'Deformational Energy level: {0}'.format(n)
print 'Vibrational energy levels: {0}'.format(n1, n2)

k = 0.524
omega1 = np.sqrt(k / m)

qe = 1.503583924
r1ini = r0 + np.sqrt(h * omega1 * (n1 + 0.5) / k)
r2ini = r0 + np.sqrt(h * omega1 * (n2 + 0.5) / k)
p1ini = 0.
p2ini = 0.
q = np.array([qe, r1ini, r2ini]).reshape((3,1))

varphi0 = 0.01
theta0 = 0.15
Jx0 = np.array([J * np.cos(varphi0) * np.sin(theta0)]).reshape((1,1))
Jy0 = np.array([J * np.sin(varphi0) * np.sin(theta0)]).reshape((1,1))
Jz0 = np.array([J * np.cos(theta0)]).reshape((1,1))
print 'Jx0: {0}; Jy0: {1}; Jz0: {2}'.format(Jx0, Jy0, Jz0)

E = h**2 * (n + (np.sqrt(Jz0**2 + Vp) + np.sqrt(Jx0**2 + Vm))/(2 * h)) * \
		   (n + 1 + (np.sqrt(Jz0**2 + Vp) + np.sqrt(Jx0**2 + Vm))/(2 * h)) / I0 + \
 	+ h * omega1 * (n1 + 0.5) + h * omega1 * (n2 + 0.5)
print 'Energy: {0}'.format(E)

effective_potential = hamiltonian(q = q, theta = theta0, varphi = varphi0, effective_potential = True)
print 'effective_potential: {0}'.format(effective_potential)

pini = np.sqrt(I0 * (E - effective_potential)).flatten()[0]
print 'pini shape: {0}'.format(pini.shape)
print 'pini: {0}'.format(pini)

p = np.array([pini, p1ini, p2ini]).reshape((3,1))
print 'p: {0}'.format(p)

print 'hamiltonian: {0}'.format(hamiltonian(q = q, p = p, theta = theta0, varphi = varphi0))

dham_dq = grad_named(hamiltonian, argname = 'q')
dham_dp = grad_named(hamiltonian, argname = 'p')
dham_dtheta = grad_named(hamiltonian, argname = 'theta')
dham_dvarphi = grad_named(hamiltonian, argname = 'varphi')

print dham_dq(q, p, theta0, varphi0)

def dham_djx(q, p, theta, varphi):
	return (1. / J) * np.cos(theta) * np.cos(varphi) * dham_dtheta(q, p, theta, varphi) - \
		   (1. / J) * np.sin(varphi) / np.sin(theta) * dham_dvarphi(q, p, theta, varphi)

def dham_djy(q, p, theta, varphi):
	return (1. / J) * np.sin(varphi) * np.cos(theta) * dham_dtheta(q, p, theta, varphi) + \
		   (1. / J) * np.cos(varphi) / np.sin(theta) * dham_dvarphi(q, p, theta, varphi)

def dham_djz(q, p, theta, varphi):
	return - (1. / J) * np.sin(theta) * dham_dtheta(q, p, theta, varphi) 

def flatten_nested_array(arr):
	arr_flatten = []
	for val in arr:
		if 'array' in str(type(val)):
			arr_flatten.extend(list(chain(*val.tolist())))
		else:
			arr_flatten.append(val)
	return arr_flatten

def rhs(t, y):
	print t, y

	eqs = [
	   lambda q, p, theta, varphi:  dham_dp(np.array([q]).reshape((__degrees__, 1)), 
										   np.array([p]).reshape((__degrees__, 1)), 
										   theta, 
										   varphi),

	   lambda q, p, theta, varphi: - dham_dq(np.array([q]).reshape((__degrees__, 1)), 
	   									  np.array([p]).reshape((__degrees__, 1)), 
	   									  theta, 
	   									  varphi),

	   lambda q, p, theta, varphi: dham_djx(np.array([q]).reshape((__degrees__, 1)),
	   										np.array([p]).reshape((__degrees__, 1)),
	   										theta, 
	   										varphi) * np.sin(varphi) - \
	   							   dham_djy(np.array([q]).reshape((__degrees__, 1)),
	   							   			np.array([p]).reshape((__degrees__, 1)),
	   							   			theta, 
	   							   			varphi) * np.cos(varphi),

	   lambda q, p, theta, varphi: (dham_djx(np.array([q]).reshape((__degrees__, 1)), 
	   									  	 np.array([p]).reshape((__degrees__, 1)), 
	   									  	 theta, 
	   									  	 varphi) * np.cos(varphi) + \
	   								dham_djy(np.array([q]).reshape((__degrees__, 1)),
	   										 np.array([p]).reshape((__degrees__, 1)),
	   										 theta, 
	   										 varphi) * np.sin(varphi)) * (1. / np.tan(theta)) - \
	   								dham_djz(np.array([q]).reshape((__degrees__, 1)), 
	   										 np.array([p]).reshape((__degrees__, 1)),
	   										 theta, varphi),
	 ]

	q = y[0:__degrees__]
	p = y[__degrees__ : 2 * __degrees__]
	theta = y[2 * __degrees__]
	varphi = y[2 * __degrees__ + 1]

	#print 'q: {0}; p: {1}; theta: {2}; varphi: {3}'.format(q, p, theta, varphi)

	derivatives = [eq(q, p, theta, varphi) for eq in eqs]

	return flatten_nested_array(derivatives)


start = time()

init = [q, p, theta0, varphi0]
print 'init: {0}'.format(init)
init_flatten = flatten_nested_array(init)

t_start = 0.
t_end = 100.
t_step = 1.

ode = spi.ode(rhs)

ode.set_integrator('lsoda', nsteps = 500, method = 'bdf', atol = 1e-6)
ode.set_initial_value(init_flatten, t_start)

sol = []
t = []
while ode.successful() and ode.t < t_end:
	ode.integrate(ode.t + t_step)
	t.append(ode.t)
	sol.append(ode.y)

print 'needed: {0}s'.format(time() - start)

# q_list = extract_column(data = sol, column = 0)
# p_list = extract_column(data = sol, column = 1)
# jx_list = extract_column(data = sol, column = 2)
# jy_list = extract_column(data = sol, column = 3)
# jz_list = extract_column(data = sol, column = 4)

# q_list = [val.reshape((__degrees__,)) for val in q_list]
# p_list = [val.reshape((__degrees__,)) for val in p_list]
# jx_list = [val.reshape((1,)) for val in jx_list]
# jy_list = [val.reshape((1,)) for val in jy_list]
# jz_list = [val.reshape((1,)) for val in jz_list]

# np.savetxt("simple-water.dat", (q_list, p_list, jx_list, jy_list, jz_list))

# plt.plot(t, q_list, 'b', label = 'q(t)')
# plt.plot(t, p_list, 'g', label = 'p(t)')
# plt.legend(loc = 'best')
# plt.xlabel('t')
# plt.grid()
# plt.show()


# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot(jx_list, jy_list, jz_list, '-r')
# plt.tight_layout()
# plt.savefig('n=0,J=25,theta=0.15.png', facecolor = 'w', edgecolor = 'b', orientation = 'portrait',
			# transparent = False, bbox_inches = None, pad_inches = 0.1, frameon = None)
