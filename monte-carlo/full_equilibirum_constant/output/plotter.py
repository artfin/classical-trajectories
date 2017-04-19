from __future__ import print_function

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def read_constants(filename):
    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    temperatures = []
    simple_constants = []

    for line in lines:
        temperatures.append(int(line.split()[0]))
        simple_constants.append(float(line.split()[1]))

    return temperatures, simple_constants

temperatures, simple_constants = read_constants('simple_constants.dat')
temperatures, general_constants = read_constants('general_constants.dat')

plt.plot(temperatures, simple_constants)
plt.plot(temperatures, general_constants)
plt.grid()

blue_patch = mpatches.Patch(color = 'blue', label = 'Simple formula')
orange_patch = mpatches.Patch(color = 'orange', label = 'General formula')

plt.legend(handles = [blue_patch, orange_patch])
#plt.show()
plt.savefig('plot.png')
