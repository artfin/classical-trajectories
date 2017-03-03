import sys
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/lib') # path on ubuntu
sys.path.append('/Users/mac/repos/sympy_project/sympy/lib') # path on mac
sys.path.append('/Users/mac/repos/sympy_project/sympy/black-box') # path on mac

import __builtin__

from __particle__ import __particle__

import autograd.numpy as np
from autograd import grad, grad_named

import scipy.integrate as spi

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from time import time

h = 1.
n = 0.
print 'Deformational Energy level: {0}'.format(n)

q0 = 1.82387
omega0 = 0.00663829
I0 = 6021.0735

Vp = 1./4 * I0**2 * omega0**2 * (1 + np.cos(q0))**2 
Vm = 1./4 * I0**2 * omega0**2 * (1 - np.cos(q0))**2
 
r0 = 1.810357468
m = I0 / r0**2

J = 25.

particle_1 = __particle__(m = m, __x__ = lambda q: -r0 * np.cos(q/2), 
								 __y__ = lambda q: 0*q, 
								 __z__ = lambda q: -r0 * np.sin(q/2),
								 vars = 'q')
particle_2 = __particle__(m = m, __x__ = lambda q: -r0 * np.cos(q/2),
								 __y__ = lambda q: 0*q,
								 __z__ = lambda q:  r0 * np.sin(q/2),
								 vars = 'q')
particles = [particle_1, particle_2]

# adding variables to builtin scope, so that they could be accessed in hamiltonian function
__builtin__.particles = particles
__builtin__.Vm = Vm
__builtin__.Vp = Vp
__builtin__.I0 = I0
__builtin__.J = J
from black_box import hamiltonian, extract_column

__degrees__ = 1

varphi0 = 0.01
theta0 = 0.15

Jx0 = np.array([J * np.cos(varphi0) * np.sin(theta0)]).reshape((1,1))
Jy0 = np.array([J * np.sin(varphi0) * np.sin(theta0)]).reshape((1,1))
Jz0 = np.array([J * np.cos(theta0)]).reshape((1,1))
print 'Jx0: {0}; Jy0: {1}; Jz0: {2}'.format(Jx0, Jy0, Jz0)

E = h**2 * (n + (np.sqrt(Jz0**2 + Vp) + np.sqrt(Jx0**2 + Vm))/(2 * h)) * \
		   (n + 1 + (np.sqrt(Jz0**2 + Vp) + np.sqrt(Jx0**2 + Vm))/(2 * h)) / I0
print 'Energy: {0}'.format(E)

qe = 1.503583924

q = np.array([qe]).reshape((1,1))

effective_potential = hamiltonian(q = q, theta = theta0, varphi = varphi0, effective_potential = True)
print 'effective_potential: {0}'.format(effective_potential)

pini = np.sqrt(I0 * (E - effective_potential))
print 'pini: {0}'.format(pini)

p = np.array([pini]).reshape((1,1))

print 'hamiltonian: {0}'.format(hamiltonian(q = q, p = p, theta = theta0, varphi = varphi0))

dham_dq = grad_named(hamiltonian, argname = 'q')
dham_dp = grad_named(hamiltonian, argname = 'p')
dham_dtheta = grad_named(hamiltonian, argname = 'theta')
dham_dvarphi = grad_named(hamiltonian, argname = 'varphi')

def dham_djx(q, p, theta, varphi):
	return (1. / J) * np.cos(theta) * np.cos(varphi) * dham_dtheta(q, p, theta, varphi) - \
		   (1. / J) * np.sin(varphi) / np.sin(theta) * dham_dvarphi(q, p, theta, varphi)

def dham_djy(q, p, theta, varphi):
	return (1. / J) * np.sin(varphi) * np.cos(theta) * dham_dtheta(q, p, theta, varphi) + \
		   (1. / J) * np.cos(varphi) / np.sin(theta) * dham_dvarphi(q, p, theta, varphi)

def dham_djz(q, p, theta, varphi):
	return - (1. / J) * np.sin(theta) * dham_dtheta(q, p, theta, varphi) 

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
	return [eq(*y) for eq in eqs]

start = time()

init = [q, p, theta0, varphi0]

t_start = 0.
t_end = 100.
t_step = 1.

ode = spi.ode(rhs)

ode.set_integrator('lsoda', nsteps = 500, method = 'bdf', atol = 1e-6)
ode.set_initial_value(init, t_start)

sol = []
t = []
while ode.successful() and ode.t < t_end:
	ode.integrate(ode.t + t_step)
	t.append(ode.t)
	sol.append(ode.y)

print 'needed: {0}s'.format(time() - start)

q_list = extract_column(data = sol, column = 0)
p_list = extract_column(data = sol, column = 1)
theta_list = extract_column(data = sol, column = 2)
varphi_list = extract_column(data = sol, column = 3)

q_list = [val.reshape((__degrees__,)) for val in q_list]
p_list = [val.reshape((__degrees__,)) for val in p_list]
theta_list = [val.reshape((1,)) for val in theta_list]
varphi_list = [val.reshape((1,)) for val in varphi_list]

jx_list = J * np.cos(varphi_list) * np.sin(theta_list)
jy_list = J * np.sin(varphi_list) * np.sin(theta_list)
jz_list = J * np.cos(theta_list)

np.savetxt("simple-water.dat", (jy_list))

# plt.plot(t, J_list, 'b', label = 'q(t)')
# plt.plot(t, p_list, 'g', label = 'p(t)')
plt.plot(t, jy_list, 'g', label = 'jx(t)')
plt.legend(loc = 'best')
plt.xlabel('t')
plt.grid()
plt.show()


# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot(jx_list, jy_list, jz_list, '-r')
# plt.tight_layout()
# plt.savefig('n=0,J=25,theta=0.15.png', facecolor = 'w', edgecolor = 'b', orientation = 'portrait',
			# transparent = False, bbox_inches = None, pad_inches = 0.1, frameon = None)
