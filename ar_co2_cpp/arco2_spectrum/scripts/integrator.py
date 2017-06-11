from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import simps
from scipy.fftpack import fft

def read_file(filename):
       
    with open(filename, mode = 'r') as inputfile:
        for index, line in enumerate(inputfile):

            if index % 1000 == 0:
                print(index)

            #if (index + 1) % 15000 == 0:
                #break

            if len(line.split()) > 0:
                if index == 0:
                    data = np.array([float(s) for s in line.split()])
                else:
                    data = np.vstack((data, np.array([float(s) for s in line.split()])))
    
    print('-'*30)
    return data

traj_file = read_file('../output/trajectory.dat')
dip_file = read_file('../output/dipole.dat')

alpha = traj_file[:,5]
beta = traj_file[:,6]
J = traj_file[0,7]

jx = J * np.cos(alpha) * np.sin(beta)
jy = J * np.sin(alpha) * np.sin(beta)
jz = J * np.cos(beta)

time = dip_file[:,0]
phi_dot = dip_file[:,1]

mux = dip_file[:,2]
muy = dip_file[:,3]
muz = dip_file[:,4]

phi = 0
for i in range(1, len(time)):
    rt = time[(i - 1):i]
    rp = phi_dot[(i - 1):i]

    # integrating phi_dot using simpson method
    # phi -- is the value of phi at the moment of time == time[i]
    phi += simps(rp, rt)
    
    _jx = jx[i]
    _jy = jy[i]
    _jz = jz[i]

    S = np.matrix([
        [_jy * np.cos(phi) - 1 / J * _jx * _jz * np.sin(phi), _jx * np.cos(phi) + 1 / J * _jy * _jz * np.sin(phi), (_jx**2 + _jy**2) / J * np.sin(phi)],
        [_jy * np.sin(phi) + 1 /J * _jx * _jz * np.cos(phi), - _jx * np.sin(phi) + 1 / J * _jy * _jz * np.cos(phi), - (_jx**2 + _jy**2) / J * np.cos(phi)],
        [_jx / J * np.sqrt(_jx**2 + _jy**2), _jy / J * np.sqrt(_jx**2 + _jy**2), _jz / J * np.sqrt(_jx**2 + _jy**2)]
    ])

    temp = S * np.array([mux[i], muy[i], muz[i]]).reshape((3, 1))
    temp = temp.reshape((1, 3))

    if i == 1:
        lmu = np.array(temp)
    else:
        lmu = np.vstack((lmu, temp))

    if i  % 1000 == 0:
        print(i)

def sq_mod(arr):
    return arr.real**2 + arr.imag**2

# number of sample points
N = len(time)
# sample spacing
T = time[1] - time[0]
tf = np.linspace(0.0, 1.0 / (2.0 * T), N/2)

lmux = lmu[:,0].flatten() 
lmuy = lmu[:,1].flatten()
lmuz = lmu[:,2].flatten()

lmuxf = fft(lmux).flatten()
lmuyf = fft(lmuy).flatten()
lmuzf = fft(lmuz).flatten()


# N//2 -- the whole part of N/2
lmuf = (2.0 / N)**2 * (np.abs(lmuxf[:N//2])**2 + np.abs(lmuyf[:N//2])**2 + np.abs(lmuzf[:N//2])**2)

fig, ax = plt.subplots()
ax.plot(tf, lmuf)
plt.show()


