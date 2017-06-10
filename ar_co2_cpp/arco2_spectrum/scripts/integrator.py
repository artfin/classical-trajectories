from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import simps

def read_file(filename):
       
    with open(filename, mode = 'r') as inputfile:
        for index, line in enumerate(inputfile):

            if index % 1000 == 0:
                print(index)

            if len(line.split()) > 0:
                if index == 0:
                    data = np.array([float(s) for s in line.split()])
                else:
                    data = np.vstack((data, np.array([float(s) for s in line.split()])))

    return data

#traj_file = read_file('../output/trajectory.dat')
dip_file = read_file('../output/dipole.dat')

time = dip_file[:,0]
phi_dot = dip_file[:,1]
phi = [] 

for i in range(1, len(time)):
    rt = time[:i]
    rp = phi_dot[:i]

    # integrating phi_dot using simpson method
    _phi = simps(rp, rt)
    phi.append(_phi)

    if i % 1000 == 0:
        print(i)

plt.plot(time[1:], phi)
plt.show()

