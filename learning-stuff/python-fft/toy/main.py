import matplotlib.pyplot as plt
import numpy as np

from pprint import pprint
from itertools import product
import random

import scipy.fftpack as sc

sample_freq = 6 # Hz

T = 300.0 # s
length = T * sample_freq * 10  
N = int(T * sample_freq)

def f(x):
    return np.sin(2 * np.pi * x) + 1.5 * np.sin(2 * np.pi * 3 * x) + 1.1 * np.sin(2 * np.pi * 7 * x) 

random_arr = np.random.random_sample(int(length))

x = np.linspace(0, T, length) 
y = f(x) # + random_arr

xs = x[::sample_freq]
ys = f(xs) # + random_arr[::sample_freq]

fig = plt.figure()
plt.rc('text', usetex = True)

plt.xlabel(r'\large \textbf{Time (s)}')
plt.ylabel(r'\large \textbf{Amplitude}')

lw = 1.5
plt.plot(x, y, color = '0.3', linestyle = 'solid', linewidth = lw)
plt.scatter(xs, ys, color = '0.5', marker='o', s = 50)
plt.grid(linestyle = ':', alpha = 0.7)
plt.show()

fs = np.arange(0, N//2 / T, 1/T)
print(len(fs))
print(length / T)
print(1 /T)

print('*'*30)

print(len(sc.fft(ys)))
ps = 2 * sc.fft(ys)[:N//2] / N 
print(len(ps))
plt.scatter(fs, abs(ps), color = '0.5', marker = 'o', s = 50)

#plt.plot([8.0, 8.0], [0.0, 25], linestyle = 'dotted', color = 'blue')
#plt.text(8.2, 20.0, r'\large \textbf{Nyquist limit}')

#plt.xlabel(r'\large \textbf{Frequency (Hz)}')
#plt.ylabel(r'\large \textbf{Amplitude}')

plt.grid(linestyle = ':', alpha = 0.7)
plt.show()



