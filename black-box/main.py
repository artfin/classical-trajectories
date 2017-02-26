import sys
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/lib') # path on ubuntu
from __particle__ import __particle__

import autograd.numpy as np
from autograd import grad, grad_named

from scipy.integrate import odeint

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
	J_vector = np.array([J * np.sin(theta) * np.cos(varphi), 
			    		 J * np.sin(theta) * np.sin(varphi), 
			   			 J * np.cos(theta)])
	p_vector = np.array([p])

	Ixx = sum([particle.m * particle.__y__(q)**2 + particle.__z__(q)**2 for particle in particles])
	Iyy = sum([particle.m * particle.__x__(q)**2 + particle.__z__(q)**2 for particle in particles])
	Izz = sum([particle.m * particle.__x__(q)**2 + particle.__y__(q)**2 for particle in particles])
	Ixy = - sum([particle.m * particle.__x__(q) * particle.__y__(q)])
	Ixz = - sum([particle.m * particle.__x__(q) * particle.__z__(q)])
	Iyz = - sum([particle.m * particle.__y__(q) * particle.__z__(q)])
	
	inertia_tensor = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
	
	a = np.array([sum([particle.m * (particle.__dx__(q)**2 + particle.__dy__(q)**2 + particle.__dz__(q)**2) for particle in particles])])

	G11 = np.linalg.inv(inertia_tensor)
	G22 = 1/a
	G12 = 0
	G21 = 0

	# print 'G11: {0}'.format(G11)
	# print 'G22: {0}'.format(G22)

	angular_component = 0.5 * np.dot(np.dot(J_vector, G11), J_vector)
	potential = Vm / (1 - np.cos(q)) / I0 + Vp / (1 + np.cos(q)) / I0
	
	if not effective_potential:
		kinetic_component = 0.5 * p_vector * G22 * p_vector 
		return angular_component + kinetic_component + potential
	else:
		return angular_component + potential


######################################################################

Vm = 1./4 * I0**2 * omega0**2 * (1 + np.cos(q0))**2 
Vp = 1./4 * I0**2 * omega0**2 * (1 - np.cos(q0))**2
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

# print hamiltonian(q = 2.0, p = 0.5, varphi = 0.1, theta = 0.2)
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
	return 1/J * ((1/(np.cos(theta) * np.cos(varphi))) * dham_dtheta(q, p, theta, varphi) - \
				   1/(np.sin(theta) * np.sin(varphi))  * dham_dvarphi(q, p, theta, varphi)) 

def dham_djy(q, p, theta, varphi):
	return 1/J * ((1/(np.cos(theta) * np.sin(varphi)) * dham_dtheta(q, p, theta, varphi) - \
				   1/(np.sin(theta) * np.sin(varphi)) * dham_dvarphi(q, p, theta, varphi)))

def dham_djz(q, p, theta, varphi):
	return - 1/(J * np.sin(theta)) * dham_dtheta(q, p, theta, varphi)

rhs = [lambda q, p, theta, varphi:  dham_dp(q, p, theta, varphi),
	   lambda q, p, theta, varphi: -dham_dq(q, p, theta, varphi),
	   lambda q, p, theta, varphi:  dham_djx(q, p, theta, varphi) * np.cos(varphi) + \
									dham_djy(q, p, theta, varphi) * np.sin(varphi) * np.tan(theta) - \
									dham_djz(q, p, theta, varphi),
	   lambda q, p, theta, varphi:  dham_djx(q, p, theta, varphi) * np.sin(varphi) - \
		   							dham_djy(q, p, theta, varphi) * np.cos(varphi)]

# y = [q, p, theta, varphi]
# rhs should return derivatives in the same order: q_dot, p_dot, theta_dot, varphi_dot
def derivatives(y, t):
	print y, t
	return [eq(y[0], y[1], y[2], y[3]) for eq in rhs]

init = [1., 0.01, 0.1, 0.5]
t = np.linspace(0, 1)
# sol = odeint(derivatives, init, t)












