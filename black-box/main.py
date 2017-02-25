import autograd.numpy as np
from autograd import grad, grad_named

from scipy.integrate import odeint

r0 = 1.
m1 = 1.2
m2 = 1.2

Vm = 624. 
Vp = 224.

J = 30.

# particle1
def __x1__(q):
	return - r0 * np.cos(q/2)

def __y1__(q):
	return 0 * q

def __z1__(q):
	return - r0 * np.sin(q/2)

# particle2
def __x2__(q):
	return - r0 * np.cos(q/2)

def __y2__(q):
	return 0 * q

def __z2__(q):
	return r0 / 2 * np.cos(q/2)

__dx1__ = grad(__x1__)
__dy1__ = grad(__y1__)
__dz1__ = grad(__z1__)

__dx2__ = grad(__x2__)
__dy2__ = grad(__y2__)
__dz2__ = grad(__z2__)

def hamiltonian(q, p, varphi, theta):
	J_vector = np.array([J * np.sin(theta) * np.cos(varphi), 
			    		 J * np.sin(theta) * np.sin(varphi), 
			   			 J * np.cos(theta)])
	p_vector = np.array([p])

	Ixx = m1 * (__y1__(q) ** 2 + __z1__(q) ** 2) + m2 * (__y2__(q) ** 2 + __z2__(q) ** 2)
	Iyy = m1 * (__x1__(q) ** 2 + __z1__(q) ** 2) + m2 * (__x2__(q) ** 2 + __z2__(q) ** 2)
	Izz = m1 * (__x1__(q) ** 2 + __y1__(q) ** 2) + m2 * (__x2__(q) ** 2 + __y2__(q) ** 2)
	Ixy = - m1 * __x1__(q) * __y1__(q) - m2 * __x2__(q) * __y2__(q)
	Ixz = - m1 * __x1__(q) * __z1__(q) - m2 * __x2__(q) * __z2__(q)
	Iyz = - m1 * __y1__(q) * __z1__(q) - m2 * __y2__(q) * __z2__(q)

	inertia_tensor = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

	a = np.array([m1 * (__dx1__(q) * __dx1__(q) + __dy1__(q) * __dy1__(q) + __dz1__(q) * __dz1__(q)) + \
				  m2 * (__dx2__(q) * __dx2__(q) + __dy2__(q) * __dy2__(q) + __dz2__(q) * __dz2__(q))])

	G11 = np.linalg.inv(inertia_tensor)
	G22 = 1/a
	G12 = 0
	G21 = 0

	angular_component = 0.5 * np.dot(np.dot(J_vector, G11), J_vector)
	potential = Vm / (1 - np.cos(q)) + Vp / (1 + np.cos(q))
	kinetic_component = 0.5 * p_vector * G22 * p_vector 
	
	return angular_component + kinetic_component + potential

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

rhs = [lambda q, p, theta, varphi:  dham_djx(q, p, theta, varphi) * np.cos(varphi) + \
									dham_djy(q, p, theta, varphi) * np.sin(varphi) * np.tan(theta) - \
									dham_djz(q, p, theta, varphi),
	   lambda q, p, theta, varphi:  dham_djx(q, p, theta, varphi) * np.sin(varphi) - \
		   							dham_djy(q, p, theta, varphi) * np.cos(varphi),
	   lambda q, p, theta, varphi:  dham_dp(q, p, theta, varphi),
	   lambda q, p, theta, varphi: -dham_dq(q, p, theta, varphi)]

# y = [q, p, theta, varphi]
def derivatives(y, t):
	print y, t
	return [eq(y[0], y[1], y[2], y[3]) for eq in rhs]

init = [1., 0.01, 0.1, 0.5]
t = np.linspace(0, 1)
sol = odeint(derivatives, init, t)












