from __future__ import division
from __future__ import print_function

import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
from potential_wrapper import potential

from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

x_all = np.linspace(5, 7, 1000)
y_r_all = np.array([potential(x, np.pi/2) for x in x_all])

u_min = min(y_all) # potential well in hartree

x_min = x_all[np.where(y_all == min(y_all))] # in atomic units
x_min_au = x_min[0]
x_min = x_min_au * 5.2917721092 * 10**(-11) # alu to m
print("x_min: {0}; x_min: {1} A".format(x_min_au, x_min * 10**(10)))

x_r = np.linspace(x_min_au + 0.1, x_min_au - 0.1, 100)
y_r = np.array([potential(x_min_au, np.pi/2) for x in x_r])

x_theta = np.linspace(np.pi/2 - 0.01, np.pi/2 + 0.01, 100)
y_theta = np.array([potential(x_min_au, theta) for theta in x_theta])

theta_min = x_plot[np.where(y_plot == min(y_plot))]
print("theta_min: {0}".format(theta_min * 360 / 2 / np.pi))

fit_r = np.polyfit(x_r, y_r, 2)
der2_r = fit_r[0]
print(fit_r)

fit_theta = np.polyfit(x_theta, y_theta, 2)
der2_theta = fit_theta[0]

print('derivative r: {0}'.format(der2_r * 15.569))
print('derivative theta: {0}'.format(der2_theta * 4.3598))

#axes = plt.gca()
#axes.set_ylim([-8.925e-4, -8.9e-4])
#plt.plot(x_red, y_red, '.', x_red, z(x_red), '-')
#plt.show()


#der2_r = der2 * 2625500 # hartree to J
#mu = 40 * 44 / 84 * 1.660 * 10**(-27) # amu to kg
#freq = np.sqrt(der2 / (2 * mu))
#freq = freq * 3.33565 * 10**(-11) # Hz to cm^-1

#print("freq: {0}".format(freq))

u_min = u_min * 219474.63 # hartree to cm^-1
print("potential well: {0}".format(u_min))




