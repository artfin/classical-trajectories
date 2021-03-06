from sympy import *
from sympy.physics.vector import dynamicsymbols

from jitcode import jitcode, provide_basic_symbols

import numpy as np
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D

from pprint import pprint 
from time import time as _time

class DynamicEquations(object):
	def __init__(self, hamiltonian, angular_momentum, angles, freedom_degrees, conjugate_momentum):
		theta, varphi = symbols('theta varphi')

		self.rhs = []
		self.rhs.append((diff(hamiltonian, angular_momentum[0]) * cos(angles[0]) + 
				    	 diff(hamiltonian, angular_momentum[1]) * sin(angles[0])) / tan(angles[1]) -
				 		 diff(hamiltonian, angular_momentum[2]))
		self.rhs.append( diff(hamiltonian, angular_momentum[0]) * sin(angles[0]) - 
						 diff(hamiltonian, angular_momentum[1]) * cos(angles[0]) )
		self.rhs.append( diff(hamiltonian, conjugate_momentum[0]))
		self.rhs.append(-diff(hamiltonian, freedom_degrees[0]))

		J = 30
		Vm = 624 
		Vp = 224

		t, y = provide_basic_symbols()

		substitutions = {'m': 1,
						 'r0': 2,
						 'Jx': J * cos(varphi) * sin(theta),
						 'Jy': J * sin(varphi) * sin(theta),
						 'Jz': J * cos(theta),
						 'Vm': Vm,
						 'Vp': Vp,
						 'varphi': y(0),
						 'theta': y(1),
						 'q': y(2),
						 'p': y(3)}

		self.rhs_w_angles = []	
		for equation in self.rhs:
			self.rhs_w_angles.append(equation.subs(substitutions))

		# varphi, theta, q, p
		init = [.1, 0.1, 3, 1]

		t = np.linspace(1, 10000, 10000)

		ODE = jitcode(self.rhs_w_angles)
		ODE.set_integrator("dopri5", nsteps = 10000)
		ODE.set_initial_value(init, 0.0)

		start = _time()
		data = []
		for time in t:
			data.append(ODE.integrate(time))
		print 'Time needed: {0}'.format(_time() - start)

		# q_list = self.extract_column(data, column = 2)
		# p_list = self.extract_column(data, column = 3)

		# plt.plot(t, q_list, 'b', label = 'q(t)')
		# plt.plot(t, p_list, 'g', label = 'p(t)')
		# plt.legend(loc = 'best')
		# plt.xlabel('t')
		# plt.grid()
		# plt.show()

		np.savetxt("timeseries.dat", data)

		# sol = odeint(self.g, init, t)

	@staticmethod
	def extract_column(data, column):
		res = []
		for row in data:
			res.append(row[column])
		return res

	def g(self, y):
		#print 'here you are: {0}'.format(t)
		vals = {'y(0)': y[0], 'y(1)': y[1], 'y(2)': y[2], 'y(3)': y[3]}
		return [eq.evalf(subs = vals) for eq in self.rhs_w_angles]


def load_data():
	with open("timeseries.dat", 'r') as inputfile:
		data = inputfile.readlines()

	print data[0]
	data = [map(float, row.split(' ')) for row in data]
	print data[0]
	return data


if __name__ == '__main__':
	Jx, Jy, Jz = symbols('Jx Jy Jz')
	q, p = symbols('q p')
	phi, theta = symbols('varphi theta')
	m, r0 = symbols('m r0')
	Vm, Vp = symbols('Vm Vp')

	angles = [phi, theta]
	angular_momentum = [Jx, Jy, Jz]
	freedom_degrees = [q]
	conjugate_momentum = [p]
	hamiltonian = Jx**2 / (2*m*r0**2 * (1 - cos(q))) + Jy**2 / (2*m*r0**2) + Jz**2 / (2*m*r0**2 * (1 + cos(q))) + p**2 / (m * r0**2) + \
				Vm / (2*m*r0**2*(1 - cos(q))) + Vp / (2*m*r0**2*(1 + cos(q)))
	
	dynamicEquations = DynamicEquations(hamiltonian = hamiltonian, angular_momentum = angular_momentum, angles = angles,
		freedom_degrees = freedom_degrees, conjugate_momentum = conjugate_momentum)

	data = load_data()

	varphi_list = DynamicEquations.extract_column(data = data, column = 0)
	theta_list = DynamicEquations.extract_column(data = data, column = 1)

	# create a sphere
	r = 30
	phi = np.linspace(0, 2 * np.pi, 100)
	theta = np.linspace(0, np.pi, 100)
	xm = r * np.outer(np.cos(phi), np.sin(theta))
	ym = r * np.outer(np.sin(phi), np.sin(theta))
	zm = r * np.outer(np.ones(np.size(phi)), np.cos(theta))

	xx = r * np.cos(varphi_list) * np.sin(theta_list)
	yy = r * np.sin(varphi_list) * np.sin(theta_list)
	zz = r * np.cos(theta_list)

	fig = plt.figure()
	ax = plt.axes(projection='3d')

	ax.plot(xx, yy, zz, '-b')

	# fig = plt.figure()
	# ax = fig.add_subplot(111, projection = '3d')
	# #ax.plot_surface(xm, ym, zm, color = 'y')	

	# ax.plot(xx, yy, zz, color = 'k', s = 20)
	# ax.set_xlim([-1, 1])
	# ax.set_ylim([-1, 1])
	# ax.set_zlim([-1, 1])
	# ax.set_aspect("equal")

	plt.tight_layout()
	plt.show()











