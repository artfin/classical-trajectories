import sys
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/lib') # path on ubuntu
sys.path.append('/Users/mac/repos/sympy_project/sympy/lib') # path on mac

from pprint import pprint

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
m = m = I0 / r0**2

J = 25.

####################################################################

particle_1 = __particle__(m = m, __x__ = lambda q: -r0 * np.cos(q/2), 
								 __y__ = lambda q: 0*q, 
								 __z__ = lambda q: -r0 * np.sin(q/2))
particle_2 = __particle__(m = m, __x__ = lambda q: -r0 * np.cos(q/2),
								 __y__ = lambda q: 0*q,
								 __z__ = lambda q:  r0 * np.sin(q/2))
particles = [particle_1, particle_2]

def hamiltonian(q = None, p = None, varphi = None, theta = None, effective_potential = False):
	J_vector = np.array([
						 J * np.sin(theta) * np.cos(varphi), 
			    		 J * np.sin(theta) * np.sin(varphi), 
			   			 J * np.cos(theta)
			   			 ])
	p_vector = np.array([p])

	Ixx = sum([particle.m * (particle.__y__(q)**2 + particle.__z__(q)**2) for particle in particles])
	Iyy = sum([particle.m * (particle.__x__(q)**2 + particle.__z__(q)**2) for particle in particles])
	Izz = sum([particle.m * (particle.__x__(q)**2 + particle.__y__(q)**2) for particle in particles])
	Ixy = - sum([particle.m * particle.__x__(q) * particle.__y__(q) for particle in particles])
	Ixz = - sum([particle.m * particle.__x__(q) * particle.__z__(q) for particle in particles])
	Iyz = - sum([particle.m * particle.__y__(q) * particle.__z__(q) for particle in particles])
	
	inertia_tensor = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

	# print 'inertia_tensor:'
	# pprint(inertia_tensor)
	
	a = np.array([sum([particle.m * (particle.__dx__(q)**2 + particle.__dy__(q)**2 + particle.__dz__(q)**2) for particle in particles])])

	G11 = np.linalg.inv(inertia_tensor)
	G22 = 1/a
	G12 = 0
	G21 = 0

	# print 'G11: {0}'.format(G11)
	# print 'G22: {0}'.format(G22)

	angular_component = 0.5 * np.dot(np.dot(J_vector, G11), J_vector)
	potential = Vm / (2 * I0 * (1 - np.cos(q))) + Vp / (2 * I0 * (1 + np.cos(q)))

	# print 'angular component: {0}'.format(angular_component)
	# print 'potential: {0}'.format(potential)
	
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

effective_potential = hamiltonian(q = qe, theta = theta0, varphi = varphi0, effective_potential = True)
print 'effective_potential: {0}'.format(effective_potential)

pini = np.sqrt(I0 * (E - effective_potential))
print 'pini: {0}'.format(pini)

print 'hamiltonian: {0}'.format(hamiltonian(q = qe, p = pini, theta = theta0, varphi = varphi0))

dham_dq = grad_named(hamiltonian, argname = 'q')
dham_dtheta = grad_named(hamiltonian, argname = 'theta')
dham_dvarphi = grad_named(hamiltonian, argname = 'varphi')
dham_dp = grad_named(hamiltonian, argname = 'p')

# print dham_dq(2.0, 0.5, 0.1, 0.2)
# print dham_dp(2.0, 0.5, 0.1, 0.2)
# print dham_dvarphi(2.0, 0.5, 0.1, 0.2)
# print dham_dtheta(2.0, 0.5, 0.1, 0.2)

rhs = []

def dham_djx(q, p, theta, varphi):
	return 1/(J * np.cos(theta) * np.cos(varphi)) * dham_dtheta(q, p, theta, varphi) - \
		   1/(J * np.sin(theta) * np.sin(varphi)) * dham_dvarphi(q, p, theta, varphi)

def dham_djy(q, p, theta, varphi):
	return 1/(J * np.cos(theta) * np.sin(varphi)) * dham_dtheta(q, p, theta, varphi) + \
		   1/(J * np.sin(theta) * np.cos(varphi)) * dham_dvarphi(q, p, theta, varphi)

def dham_djz(q, p, theta, varphi):
	return - 1/(J * np.sin(theta)) * dham_dtheta(q, p, theta, varphi)

rhs = [lambda q, p, theta, varphi:   dham_dp(q, p, theta, varphi),
	   lambda q, p, theta, varphi: - dham_dq(q, p, theta, varphi),
	   lambda q, p, theta, varphi:   dham_djx(q, p, theta, varphi) * np.sin(varphi) - dham_djy(q, p, theta, varphi) * np.cos(varphi),
	   lambda q, p, theta, varphi:  (dham_djx(q, p, theta, varphi) * np.cos(varphi) + \
									 dham_djy(q, p, theta, varphi) * np.sin(varphi)) * (1./np.tan(theta)) - \
									 dham_djz(q, p, theta, varphi)
	   ]

# y = [q, p, theta, varphi]
# rhs should return derivatives in the same order: q_dot, p_dot, theta_dot, varphi_dot
def derivatives(t, y):
	print t, y
	#print dham_dtheta(y[0], y[1], y[2], y[3])
	#print t, dham_djz(y[0], y[1], y[2], y[3])
	# print dham_dq(y[0], y[1], y[2], y[3])
	return [eq(y[0], y[1], y[2], y[3]) for eq in rhs]

init = [qe, pini, theta0, varphi0]
# print dham_djz(init[0], init[1], init[2], init[3])
# t = np.linspace(0, 500, 500)
# sol = odeint(derivatives, init, t, atol = 10**(-6))

t_start = 0.
t_end = 700.
t_step = 1.

ode = spi.ode(derivatives)

ode.set_integrator('lsoda', nsteps = 500, method = 'bdf', atol = 1e-6)
ode.set_initial_value(init, t_start)

sol = []
t = []
while ode.successful() and ode.t < t_end:
	ode.integrate(ode.t + t_step)
	t.append(ode.t)
	sol.append(ode.y)


# saving it just in case
np.savetxt("simple-water.dat", sol)

q_list = extract_column(data = sol, column = 0)
p_list = extract_column(data = sol, column = 1)
theta_list = extract_column(data = sol, column = 2)
varphi_list = extract_column(data = sol, column = 3)

# plt.plot(t, q_list, 'b', label = 'q(t)')
# plt.plot(t, p_list, 'g', label = 'p(t)')

# plt.plot(t, theta_list, 'b', label = 'theta(t)')
plt.plot(t, varphi_list, 'g', label = 'varphi(t)')
plt.legend(loc = 'best')
plt.xlabel('t')
plt.grid()
plt.show()

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

# ax.plot(xx, yy, zz, '-b')
# plt.tight_layout()
# plt.show()













