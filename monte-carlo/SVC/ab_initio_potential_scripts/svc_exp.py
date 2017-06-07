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
    _integrand = partial(integrand, Temperature = T)
    
    start = time()
    integ = vegas.Integrator([[lb, hb], [0., np.pi]])
    result = integ(_integrand, nitn = 10, neval = 5 * 10**4)
    print 'result = %s Q = %.2f' % (result, result.Q)

    SVC = np.pi * avogadro * result.mean * length_unit**3 * 10**(-4) # 10**6 to convert from m3/mol to cm3/mol
    
    print 'Temperature: %d; SVC: %.5f' % (T, SVC)
    print 'Time needed: {0}'.format(time() - start)

    return SVC

def save_data(file_path, temperatures, svcs):
    with open(file_path, mode = 'w') as out:
        for temperature, svc in zip(temperatures, svcs):
            out.write(str(temperature) + ' ' + str(svc) + '\n')

def calculate_svcs():
    temperatures = [100 + i for i in range(400)]

    start_values = range(5, 6, 1)
    end_values = range(6, 7, 1)

    for start, end in zip(start_values, end_values):
        print('Start: {0}; end: {1}'.format(start, end))
        svcs = [cycle(start, end, temperature) for temperature in temperatures]
        save_data(str(int(start)) + str(int(end)) + 'svc.txt', temperatures, svcs)
    

#calculate_svcs()

def tr(arr):
    return [np.log(el) for el in arr]

def tr2(arr):
    return [np.log(-el) for el in arr]

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

def plot_psvc():
    
    fig = plt.figure()
    plt.rc('text', usetex = True)
    lw = 1.75

    ax = plt.subplot(1, 1, 1)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    plt.tick_params(axis="both", which="both", bottom="off", left="off")

    for start, end in zip(range(5, 6), range(6, 7)):
        filename = str(int(start)) + str(int(end)) + 'svc.txt'
        psvc_ab, temperatures = read_data(filename)

        l1, = ax.plot(temperatures, psvc_ab, color = '0.5', linestyle = 'solid', linewidth = lw)

        filename = str(int(start)) + str(int(end)) + 'HUTsvc.txt'
        psvc_hut, temperatures = read_data(filename)

        l2, = ax.plot(temperatures, psvc_hut, color = 'red', linestyle = 'solid', linewidth = lw)
        
        y_pos = (psvc_ab[-1] + psvc_hut[-1]) / 2
        text = "[" + str(start) + ", " + str(end) + "]"

        if start == 7:
            y_pos += 0.05
        if start == 9:
            y_pos -= 0.05

        plt.text(temperatures[-1] + 0.03, y_pos - 0.05, text, fontsize=12)

    ax.set_xlabel(r'ln \textbf{T}')
    ax.set_ylabel(r'ln \textbf{SVC} (partial)')

    fig.legend((l1, l2), ('Ab-initio', 'Hutson'), 'lower center', ncol = 2, fancybox = True, shadow = True, prop = {'size': 'large'})

    ax.grid(linestyle = ':', alpha = 0.7)
    plt.show()

plot_psvc()


