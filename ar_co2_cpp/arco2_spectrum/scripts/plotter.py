import matplotlib.pyplot as plt
import numpy as np

def read_dipole(filename):
    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    t0 = np.array([])
    dipx = np.array([])
    dipy = np.array([])
    dipz = np.array([])

    for line in lines:
        
        data = line.split()
        if len(data) > 1:
            
            t0 = np.append(t0, float(data[0]))
            dipx = np.append(dipx, float(data[1]))
            dipy = np.append(dipy, float(data[2]))
            dipz = np.append(dipz, float(data[3]))

    return t0, dipx, dipy, dipz

t0, dipx, dipy, dipz = read_dipole('../_dip.txt')

lw = 1.5
plt.plot(t0, dipx, color = 'k', linewidth = lw, linestyle = 'solid')
plt.plot(t0, dipy, color = 'r', linewidth = lw, linestyle = 'solid')
plt.plot(t0, dipz, color = 'y', linewidth = lw, linestyle = 'solid')
plt.show()

