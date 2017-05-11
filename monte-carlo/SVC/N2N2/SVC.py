from __future__ import print_function

import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/N2N2/')

from n2n2 import potinit, potn2n2

# initializing procedure
potinit()

# converts cm^-1 to hartree
def potential(r, theta1, theta2, phi):
    return potn2n2(r, theta1, theta2, phi) * 4.55633 * 10**(-6)

import numpy as np
import vegas
from time import time
from functools import partial

import matplotlib.pyplot as plt

# CONSTANTS
# ---------------------------------------
k = 1.38064852 * 10**(-23) # J/k
htoj = 4.35974417 * 10**(-18) # hartree to Joules
avogadro = 6.022140 * 10**(23) 
length_unit = 5.291772 * 10**(-11)
# ---------------------------------------

def integrand(x, temperature):
    # x = [r, theta1, theta2, phi]
    r = x[0]
    theta1 = x[1]
    theta2 = x[2]
    potential_value = potential(*x) * htoj
    
    if r > 4.4:
        return (1 - np.exp(-potential_value / (k * temperature))) * r**2 * np.sin(theta1) * np.sin(theta2)
    else:
        return r**2 * np.sin(theta1) * np.sin(theta2) 

def cycle(temperature):
    _integrand = partial(integrand, temperature = temperature)

    start = time()
    integ = vegas.Integrator([[0.0, 45.0], [0, np.pi], [0, np.pi], [0.0, 2 * np.pi]])
    result = integ(_integrand, nitn = 10, neval = 5 * 10**4)
    print('result: {0}; Q: {1}'.format(result, result.Q))

    SVC = avogadro / 4 * result.mean * length_unit**3 * 10**6
    print('temperature: {0}; SVC: {1}'.format(temperature, SVC))
    print('Time needed: {0}'.format(time() - start))
    print('*'*30)

    return SVC

def save_data(temperatures, svcs):
    file_path = '/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/SVC/N2N2/SVC.dat'
    with open(file_path, mode = 'w') as out:
        for temperature, svc in zip(temperatures, svcs):
            out.write(str(temperature) + ' ' + str(svc) + '\n')

def read_data(filename):
    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    temperatures = []
    svcs = []

    for line in lines:
        if len(line) > 1:
            data = line.split()

            temperatures.append(float(data[0]))
            svcs.append(float(data[1]))

    return temperatures, svcs

temperatures = [75.0 + i * 5.0 for i in range(50)]

#svcs = [cycle(temperature) for temperature in temperatures]
#save_data(temperatures, svcs)

article_temperatures = [75.0, 80.0, 90.0, 100.0, 110.0, 125.0, 150.0, 200.0, 250.0, 300.0]
article_svcs = [-280.2, -247.0, -196.8, -161.0, -134.1, -104.4, -71.77, -35.79, -16.58, -4.75]

temperatures, svcs = read_data('SVC.dat')

plt.plot(temperatures, svcs, '--', color = 'k')
plt.scatter(article_temperatures, article_svcs, marker = '*', color = 'r')
plt.grid()
plt.show()






