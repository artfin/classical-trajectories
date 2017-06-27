from __future__ import print_function

import matplotlib.pyplot as plt

def read_data(filename):

    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    temperatures = []
    moments = []

    for line in lines:
        
        data = line.split()

        if len(data) > 1:
            temperatures.append(float(data[0]))
            moments.append(float(data[1]))

    return temperatures, moments

temperatures_t, moments_t = read_data(filename = 'theory.txt')
temperatures_e, moments_e = read_data(filename = 'experiment.txt')

temperatures_e1 = [temperatures_e[0], temperatures_e[2]]
moments_e1 = [moments_e[0], moments_e[2]]

temperatures_e2 = [temperatures_e[1], temperatures_e[3]]
moments_e2 = [moments_e[1], moments_e[3]]

fig = plt.figure()

plt.rc('text', usetex = True)

_size = 45
lw = 1.75
plt.title(r'\Large \textbf{Second spectral moment}')

plt.xlabel(r'\large \textbf{T}, (K)')
plt.ylabel(r'\large \textbf{M}$_0$, (cm$^{-3}$ $\cdot$ amagat$^{-1}$)')

l, = plt.plot(temperatures_t, moments_t, color = '0.6', linewidth = lw)
e1 = plt.scatter(temperatures_e1, moments_e1, color = 'k', marker = '^', s = _size)
e2 = plt.scatter(temperatures_e2, moments_e2, color = 'k', marker = 's', s = _size)

fig.legend((l, e1, e2), ('Calc.', r'I. R. Dagg \textit{et al.}, \textit{Can. J. Phys.} \textbf{64}, 1485, (1986)', r'M. V. Tonkov In: \textit{Collision- and Interaction-Induced Spectroscopy}. Kluwer, AP, 1995'), 'lower center',  ncol = 3, fancybox = True, shadow = True, prop = {'size': 'large'}) 


plt.grid(linestyle = ':', alpha = 0.7)
plt.show()


