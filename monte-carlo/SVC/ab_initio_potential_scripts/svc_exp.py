import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers')

from potential_wrapper import potential
from hutson import potdat4, extpot

import numpy as np
import vegas
from time import time
from functools import partial

import matplotlib.pyplot as plt

# loading parameters
potdat4()

k = 1.38064852 * 10**(-23) # J/k
htoj = 4.35974417 * 10**(-18) # hartree to Joules
avogadro = 6.022 # * 10**(23) 
length_unit = 5.291772 # * 10**(-11)

def hpot(R, theta):
    return extpot(R, np.cos(theta))

def integrand(x, Temperature):
    # x = [R, theta]
    potential_value = potential(x[0], x[1]) * htoj
    value =  (1 - np.exp(- potential_value / (k * Temperature))) * x[0]**2 * np.sin(x[1])
    return value

def hint(x, Temperature):
    # x = [R, theta]
    potential_value = hpot(x[0], x[1]) * htoj
    value = (1 - np.exp(- potential_value / (k * Temperature))) * x[0]**2 * np.sin(x[1])
    return value

def cycle(lb, hb, T):
    _integrand = partial(hint, Temperature = T)
    
    start = time()
    integ = vegas.Integrator([[lb, hb], [0., np.pi]])
    result = integ(_integrand, nitn = 10, neval = 10**3)
    print 'result = %s Q = %.2f' % (result, result.Q)

    SVC = np.pi * avogadro * result.mean * length_unit**3 * 10**(-4) # 10**6 to convert from m3/mol to cm3/mol
    
    print 'Temperature: %d; SVC: %.5f' % (T, SVC)
    print 'Time needed: {0}'.format(time() - start)

    return SVC

def save_data(file_path, temperatures, svcs):
    with open(file_path, mode = 'a') as out:
        for temperature, svc in zip(temperatures, svcs):
            out.write(str(temperature) + ' ' + str(svc) + '\n')

def calculate_svcs():
    temperatures = [100 + i for i in range(400)]

    fig = plt.figure()

    plt.rc('text', usetex = True)

    start_values = range(15)
    end_values = range(1, 16)

    for start, end in zip(start_values, end_values):
        print('Start: {0}; end: {1}'.format(start, end))
        svcs = [cycle(start, end, temperature) for temperature in temperatures]
        save_data(str(int(start)) + str(int(end)) + 'HUTsvc.txt', temperatures, svcs)
    
        l, = plt.plot(temperatures, svcs, linestyle = 'solid')

    plt.grid()
    plt.show()

calculate_svcs()

def read_data(filename):
    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    temperatures, svcs = [], []

    for line in lines:
        if len(line) > 1:
            data = line.split()
            temperatures.append(float(data[0]))
            svcs.append(float(data[1]))

    return svcs, temperatures

def diff(l1, l2):
    res = []
    for e1, e2 in zip(l1, l2):
        res.append(e1 - e2)
    return res

def plot_diff():

    psvc = []

    for i in range(1, 11):
        filename = str(float(i)) + 'svc.txt'
        _, temperatures = read_data(filename)
        psvc.append(_)

    fig = plt.figure()
    plt.rc('text', usetex = True)

    col = 0.1
    for p1, p2, n in zip(psvc[1:], psvc[2:], range(1, 11)):
        col += 0.05
        print('given color: {0}'.format(col))

        d = diff(p1, p2)
        l, = plt.plot(temperatures, d, color = str(col), linestyle = 'solid')

        print('First element: {0}'.format(d[0]))

    plt.grid()
    plt.show()




