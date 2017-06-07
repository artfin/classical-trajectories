from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

def read_constants(filename):
    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    temperatures = []
    constants = []

    for line in lines:
        temperatures.append(int(line.split()[0]))
        constants.append(np.log(float(line.split()[1])))

    return temperatures, constants

temperatures, ab_constants = read_constants('final.ab_initio.dat')
temperatures, hut_constants = read_constants('final.hutson.dat')
temperatures, par_constants = read_constants('final.parker.dat')

fig = plt.figure()
plt.rc('text', usetex = True)

lw = 1.75
l1, = plt.plot(temperatures, par_constants, color = '0.4', linestyle = 'solid', linewidth = lw)
l2, = plt.plot(temperatures, ab_constants, color = '0.6', linestyle = 'dashed', linewidth = lw)
l3, = plt.plot(temperatures, hut_constants, color = '0.5', linestyle = 'dotted', linewidth = lw)

plt.xlabel(r'\textbf{T}, (K)')
plt.ylabel(r'ln \textbf{K}$_p$')

fig.legend((l1, l2, l3), ('Parker', 'Ab-initio', 'Hutson'), 'lower center', ncol = 3, fancybox = True, shadow = True, prop = {'size': 'large'})

plt.grid(linestyle = ':', alpha = 0.7)
plt.show()
