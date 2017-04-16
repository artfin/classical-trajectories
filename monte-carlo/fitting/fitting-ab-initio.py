from __future__ import division
from __future__ import print_function

import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/derivative_wrapper')
from derivatives_wrapper import derivative_R as derR
from derivatives_wrapper import derivative_Theta as derT
from potential_wrapper import potential

from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(6, 7, 10000)
dh_dr = np.array([derR(_x, np.pi/2) for _x in x])
dh_dt = np.array([derT(_x, np.pi/2) for _x in x])

nums = np.where(dh_dr == min(abs(dh_dr))), np.where(dh_dr == -min(abs(dh_dr)))
print('minimum: {0}, {1}'.format(x[nums[0]], x[nums[1]]))
# derivative goes through zero at around 6.4966 bohr == potential minimum

htocm = 219474.63 # hartree to cm^-1
r_well = x[nums[0]][0] # r_min in bohr 
print('r-well: {0}'.format(r_well * 5.29 * 10**(-11)))
well = potential(r_well, np.pi/2) * htocm 
print('potential well: {0}'.format(well))

x = np.linspace(6.48, 6.50, 1000)
y = np.array([derR(_x, np.pi/2) for _x in x])

# linear fit
k_r, b = np.polyfit(x, y, 1)
print('k_r: {0}'.format(k_r))

# k is the d2u/dr2 in atomic units (hartree / alu**2)
# so translating it into CI units

x = np.linspace(np.pi/2 - 0.001, np.pi/2 + 0.001, 1000)
y = np.array([derT(6.4966, _x) for _x in x])

# linear fit
k_t, b = np.polyfit(x, y, 1)
print ('k_t: {0}'.format(k_t)) 

htoj = 4.359 * 10**(-20) # hartree to joules
alu = 5.291 * 10**(-11) # atomic length unit to m
da = 1.660 * 10**(-27) # dalton to kg
hztocm = 3.335 * 10**(-11) # Hz to cm**(-1)
mu = (40 * 44) / (40 + 44) * da # mass of ar-co2 complex in kg
freqau = 6.57968 * 10**(15) # freq from atomic to s^-1

k_r = k_r * htoj / alu**2
freq_r = np.sqrt(k_r / mu) # frequency in Hz
print('freq_r: {0} Hz'.format(freq_r))
freq_r = freq_r * hztocm
print('freq_r: {0} cm^-1'.format(freq_r))

#freq_t = np.sqrt(k_t / ((40 * 44) / (40 + 44) * 1836.42)) * freqau # in atomic units
#print('freq_t: {0} Hz'.format(freq_t))
#freq_t = freq_t * hztocm
#print('freq_t: {0} cm^-1'.format(freq_t))

d0 = well + freq_r
print('D0: {0}'.format(d0))

#axes = plt.gca()
#axes.set_ylim([-1, 1])
#plt.plot(x, dh_dr)
#plt.show()

#x_all = np.linspace(5, 7, 1000)
#y_r_all = np.array([potential(x, np.pi/2) for x in x_all])

#u_min = min(y_all) # potential well in hartree

#x_min = x_all[np.where(y_all == min(y_all))] # in atomic units
#x_min_au = x_min[0]
#x_min = x_min_au * 5.2917721092 * 10**(-11) # alu to m
#print("x_min: {0}; x_min: {1} A".format(x_min_au, x_min * 10**(10)))

#x_r = np.linspace(x_min_au + 0.1, x_min_au - 0.1, 100)
#y_r = np.array([potential(x_min_au, np.pi/2) for x in x_r])

#x_theta = np.linspace(np.pi/2 - 0.01, np.pi/2 + 0.01, 100)
#y_theta = np.array([potential(x_min_au, theta) for theta in x_theta])

#theta_min = x_plot[np.where(y_plot == min(y_plot))]
#print("theta_min: {0}".format(theta_min * 360 / 2 / np.pi))

#fit_r = np.polyfit(x_r, y_r, 2)
#der2_r = fit_r[0]
#print(fit_r)

#fit_theta = np.polyfit(x_theta, y_theta, 2)
#der2_theta = fit_theta[0]

#print('derivative r: {0}'.format(der2_r * 15.569))
#print('derivative theta: {0}'.format(der2_theta * 4.3598))

#axes = plt.gca()
#axes.set_ylim([-8.925e-4, -8.9e-4])
#plt.plot(x_red, y_red, '.', x_red, z(x_red), '-')
#plt.show()


#der2_r = der2 * 2625500 # hartree to J
#mu = 40 * 44 / 84 * 1.660 * 10**(-27) # amu to kg
#freq = np.sqrt(der2 / (2 * mu))
#freq = freq * 3.33565 * 10**(-11) # Hz to cm^-1

#print("freq: {0}".format(freq))

#u_min = u_min * 219474.63 # hartree to cm^-1
#print("potential well: {0}".format(u_min))




