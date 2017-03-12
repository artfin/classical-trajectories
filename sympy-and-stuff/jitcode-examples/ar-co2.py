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

class DynamicEquations(object):
	def __init__(self, hamiltonian, angular_momentum, angles, freedom_degrees, conjugate_momentum):
		theta, varphi = symbols('theta varphi')

		J = 30.

		self.rhs = []
		self.rhs.append((diff(hamiltonian, angular_momentum[0]) * cos(angles[0]) + 
				    	 diff(hamiltonian, angular_momentum[1]) * sin(angles[0])) / tan(angles[1]) -
				 		 diff(hamiltonian, angular_momentum[2]))
		self.rhs.append( diff(hamiltonian, angular_momentum[0]) * sin(angles[0]) - 
						 diff(hamiltonian, angular_momentum[1]) * cos(angles[0]) )
		self.rhs.append( diff(hamiltonian, conjugate_momentum[0]))
		self.rhs.append(-diff(hamiltonian, freedom_degrees[0]))

		t, y = provide_basic_symbols()

		substitutions = {
						 'Jx': J * cos(varphi) * sin(theta),
						 'Jy': J * sin(varphi) * sin(theta),
						 'Jz': J * cos(theta),
						 'varphi': y(0),
						 'theta': y(1),
						 'q': y(2),
						 'p': y(3)
						 }

		self.rhs_w_angles = []	
		for equation in self.rhs:
			self.rhs_w_angles.append(equation.subs(substitutions))

		# varphi, theta, q, p
		init = [.1, 0.1, 3, 1]

		t = np.linspace(1, 1000, 1000)

		ODE = jitcode(self.rhs_w_angles)
		ODE.set_integrator("dopri5", nsteps = 10000)
		ODE.set_initial_value(init, 0.0)

		data = []
		for time in t:
			data.append(ODE.integrate(time))

		np.savetxt("timeseries.dat", data)

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


Jx, Jy, Jz = symbols('Jx Jy Jz')
R, theta = symbols('R theta')
pR, pT = symbols('pR pT')

angular_momentum = [Jx, Jy, Jz]
freedom_degrees = [R, theta]
conjugate_momentum = [pR, pT]

hamiltonian = Jx**2 / (2*m*r0**2 * (1 - cos(q))) + Jy**2 / (2*m*r0**2) + Jz**2 / (2*m*r0**2 * (1 + cos(q))) + p**2 / (m * r0**2) + \
			Vm / (2*m*r0**2*(1 - cos(q))) + Vp / (2*m*r0**2*(1 + cos(q)))

dynamicEquations = DynamicEquations(hamiltonian = hamiltonian, angular_momentum = angular_momentum, angles = angles,
	freedom_degrees = freedom_degrees, conjugate_momentum = conjugate_momentum)

data = load_data()

varphi_list = DynamicEquations.extract_column(data = data, column = 0)
theta_list = DynamicEquations.extract_column(data = data, column = 1)

