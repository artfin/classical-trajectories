import sys
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/automatic-solver') # ubuntu
sys.path.append('/Users/mac/repos/sympy_project/sympy/automatic-solver') # path on mac

from __particle__ import __particle__
from automatic_solver import AutomaticSolver
import autograd.numpy as np
# import matplotlib.pyplot as plt
from time import time

h = 1.
q0 = 1.82387
qe = 1.503583924
omega0 = 0.00663829
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
harmonic_constant = 0.524

J = 25.

h = 1.
n = 0.
n1 = 4.
n2 = 1.
print 'Deformational Energy level: {0}'.format(n)
print 'Vibrational energy levels: {0}'.format(n1, n2)

omega1 = np.sqrt(harmonic_constant / m)

qe = 1.503583924
r1ini = r0 + np.sqrt(h * omega1 * (n1 + 0.5) / harmonic_constant)
r2ini = r0 + np.sqrt(h * omega1 * (n2 + 0.5) / harmonic_constant)
p1ini = 0.
p2ini = 0.
q = np.array([qe, r1ini, r2ini]).reshape((3,1))

varphi0 = 0.0
theta0 = 0.15
Jx0 = np.array([J * np.cos(varphi0) * np.sin(theta0)]).reshape((1,1))
Jy0 = np.array([J * np.sin(varphi0) * np.sin(theta0)]).reshape((1,1))
Jz0 = np.array([J * np.cos(theta0)]).reshape((1,1))

E = h**2 * (n + (np.sqrt(Jz0**2 + Vp) + np.sqrt(Jx0**2 + Vm))/(2 * h)) * \
		   (n + 1 + (np.sqrt(Jz0**2 + Vp) + np.sqrt(Jx0**2 + Vm))/(2 * h)) / I0 + \
 	+ h * omega1 * (n1 + 0.5) + h * omega1 * (n2 + 0.5)
print 'Energy: {0}'.format(E)

q = np.array([qe, r1ini, r2ini]).reshape((3,1))

######################################################################################
def potential(q):
	return Vm / 2 / I0 / (1 - np.cos(q[0])) + Vp / 2 / I0 / (1 + np.cos(q[0])) + \
		   harmonic_constant * (q[1] - r0) ** 2 + harmonic_constant * (q[2] - r0)**2
######################################################################################

AS = AutomaticSolver(particles = particles, __degrees__ = __degrees__, potential = potential)
# setting angular momentum in AS!
AS.J = J

effective_potential = AS.hamiltonian(q = q, theta = theta0, varphi = varphi0, effective_potential = True)
print 'effective_potential: {0}'.format(effective_potential)

I_modified = 2 * m * r1ini**2 * r2ini**2 / (r1ini**2 + r2ini**2)
pini = np.sqrt(I_modified * (E - effective_potential)).flatten()[0]
print 'pini shape: {0}'.format(pini.shape)
print 'pini: {0}'.format(pini)

p = np.array([pini, p1ini, p2ini]).reshape((3,1))
print 'p: {0}'.format(p)

print 'hamiltonian: {0}'.format(AS.hamiltonian(q = q, p = p, theta = theta0, varphi = varphi0))

####################################################################################

init = [q, p, theta0, varphi0]
t_start = 0.
t_end = 10000.
t_step = 1.
t = np.linspace(t_start, t_end, t_end / t_step)

start = time()
solution = AS.integrate(initial_conditions = init, t_start = t_start, t_end = t_end, t_step = t_step)
print 'Time needed: {0}s'.format(time() - start)

# q_list = AS.extract_column(data = solution, column = 0)
# p_list = AS.extract_column(data = solution, column = 1)
theta_list = AS.extract_column(data = solution, column = 2 * __degrees__)
varphi_list = AS.extract_column(data = solution, column = 2 * __degrees__ + 1)

jx_list = J * np.cos(varphi_list) * np.sin(theta_list)
jy_list = J * np.sin(varphi_list) * np.sin(theta_list)
jz_list = J * np.cos(theta_list)

AS.save_file(filename = 'angular_trajectory_J=25.dat', jx = jx_list, jy = jy_list, jz = jz_list)

# plt.plot(t, theta_list, 'b', label = 'q(t)')
# plt.plot(t, p_list, 'g', label = 'p(t)')
# plt.legend(loc = 'best')
# plt.xlabel('t')
# plt.grid()
# plt.show()
		
