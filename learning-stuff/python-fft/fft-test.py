import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft

def fourier_spectrum(X, sample_spacing, min_period):
    ps = np.abs(fft(X))**2
    freqs = np.fft.fftfreq(X.size, sample_spacing)
    idx = np.argsort(freqs)

    plt.plot(freqs[idx], ps[idx])
    plt.show()

 #number of samplepoints
N = 200 
 #sample spacing
T = 1.0

x = np.linspace(0.0, N * T, N)
y = np.sin(2.0 * np.pi * x) 
yf = fft(y)
xf = np.linspace(0.0, 1.0 / (2.0 * T), N/2)

#fig, ax = plt.subplots()
#ax.plot(xf, 2.0/N * np.abs(yf[:N//2]))
#plt.show()

fourier_spectrum(y, x, 0.01)
