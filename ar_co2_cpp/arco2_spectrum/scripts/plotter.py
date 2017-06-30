import matplotlib.pyplot as plt
import numpy as np

def read_traj(filename):
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

def read_dip(filename):
    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    dip = np.array([])
    dip1 = np.array([])

    for line in lines:
        if len(line) > 1:
            data = line.split()

            dip = np.append(dip, float(data[0]))
            dip1 = np.append(dip1, float(data[1]))
    
    return dip, dip1

dip, dip1 = read_dip('../dipfft.txt')

lw = 1.5
plt.plot(dip, color = 'k', linewidth = lw)
plt.plot(dip1, color = 'r', linewidth = lw)
plt.show()

#t0, dipx, dipy, dipz = read_dipole('../_dip.txt')

#lw = 1.5
#plt.plot(t0, dipx, color = 'k', linewidth = lw, linestyle = 'solid')
#plt.plot(t0, dipy, color = 'r', linewidth = lw, linestyle = 'solid')
#plt.plot(t0, dipz, color = 'y', linewidth = lw, linestyle = 'solid')
#plt.show()

