import sys
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/lib') # path on ubuntu
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

__degrees__ = 1

n = 0
varphi0 = 0.01
theta0 = 0.15

Jx0 = np.array([J * np.cos(varphi0) * np.sin(theta0)]).reshape((1,1))
Jy0 = np.array([J * np.sin(varphi0) * np.sin(theta0)]).reshape((1,1))
Jz0 = np.array([J * np.cos(theta0)]).reshape((1,1))

E = h**2 * (n + (np.sqrt(Jz0**2 + Vp) + np.sqrt(Jx0**2 + Vm))/(2 * h)) * \
		   (n + 1 + (np.sqrt(Jz0**2 + Vp) + np.sqrt(Jx0**2 + Vm))/(2 * h)) / I0
print 'Energy: {0}'.format(E)

def potential(q):
	return Vm / 2 / I0 / (1 - np.cos(q)) + Vp / 2 / I0 / (1 + np.cos(q))

AS = AutomaticSolver(particles = particles, __degrees__ = __degrees__, potential = potential)
# setting angular momentum in AS!
AS.J = J

q = np.array([qe]).reshape((1,1))

effective_potential = AS.hamiltonian(q = q, theta = theta0, varphi = varphi0, effective_potential = True)
print 'Effective potential: {0}'.format(effective_potential)

pini = np.sqrt(I0 * (E - effective_potential))
print 'pini: {0}'.format(pini)

p = np.array([pini]).reshape((1,1))

init = [q, p, theta0, varphi0]
t_start = 0.
t_end = 500.
t_step = 1.
t = np.linspace(t_start, t_end, t_end / t_step)

start = time()
solution = AS.integrate(initial_conditions = init, t_start = t_start, t_end = t_end, t_step = t_step)
print 'Time needed: {0}s'.format(time() - start)

q_list = AS.extract_column(data = solution, column = 0)
p_list = AS.extract_column(data = solution, column = 1)
theta_list = AS.extract_column(data = solution, column = 2)
varphi_list = AS.extract_column(data = solution, column = 3)

jx_list = J * np.cos(varphi_list) * np.sin(theta_list)
jy_list = J * np.sin(varphi_list) * np.sin(theta_list)
jz_list = J * np.cos(theta_list)

AS.save_file(filename = 'angular_trajectory.dat', jx = jx_list, jy = jy_list, jz = jz_list)

# plt.plot(t, q_list, 'b', label = 'q(t)')
# plt.plot(t, p_list, 'g', label = 'p(t)')
# plt.legend(loc = 'best')
# plt.xlabel('t')
# plt.grid()
# plt.show()
		
