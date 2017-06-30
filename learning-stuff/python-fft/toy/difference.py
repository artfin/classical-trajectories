from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

import scipy.fftpack as sc

def computeDFT(inreal):
    
    outreal = np.array([])
    outimag = np.array([])

    n = len(inreal)
    for k in xrange(n):

        sumreal = sumimag = 0
        for t in xrange(n):
            angle = - t * k / n

            sumreal += inreal[t] * np.cos(angle)
            sumimag -= inreal[t] * np.sin(angle)

        outreal = np.append(outreal, sumreal)
        outimag = np.append(outimag, sumimag)

    return outreal, outimag

N = 300
x = np.linspace(0, 10, N)
inreal = np.sin(5 * x)

outreal1, outimag1 = computeDFT(inreal)
pow1 = np.sqrt(outreal1**2 + outimag1**2)

ps = abs(sc.fft(inreal))

lw = 1.5
plt.scatter(range(len(pow1[:N//2])), pow1[:N//2], color = 'k')
plt.scatter(range(len(ps[:N//2])), ps[:N//2], color = 'r')
plt.show()


