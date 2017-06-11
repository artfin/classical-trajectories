import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft

# number of samplepoints
N = 4 
# sample spacing
T = 1.0

x = np.linspace(0.0, N * T, N)
y = np.sin(50.0 * 2.0 * np.pi * x) + 0.5 * np.sin(80.0 * 2.0 * np.pi * x)
yf = fft(y)
xf = np.linspace(0.0, 1.0 / (2.0 * T), N/2)

print(x)
print(y)
print('-'*30)
print(xf)
print(yf)

#fig, ax = plt.subplots()
#ax.plot(xf, 2.0/N * np.abs(yf[:N//2]))
#plt.show()

