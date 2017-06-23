import matplotlib.pyplot as plt
import numpy as np

signal_freq = 1.0 # Hz
sample_freq = 9.0 # Hz

T = 1.0 # s

x = np.linspace(0, T, 100)
y = np.sin(signal_freq * x * 2 * np.pi)

xs = np.linspace(0, T, T * sample_freq)
ys = np.sin(xs * sample_freq * 2 * np.pi)

print(ys)

fig = plt.figure()
plt.rc('text', usetex = True)

plt.xlabel(r'\large \textbf{Frequency (Hz)}')
plt.ylabel(r'\large \textbf{Amplitude}')

#lw = 1.5
#plt.plot(x, y, color = '0.3', linestyle = 'solid', linewidth = lw)
#plt.scatter(xs, ys, color = '0.5', marker='o', s = 50)
#plt.grid(linestyle = ':', alpha = 0.7)
#plt.show()

xs = np.array([0.0, 1.0, 2.0, 3.0])
ps = np.array([0.0, 1.0, 0.0, 0.0])

plt.scatter(xs, ps, color = '0.5', marker = 'o', s = 50)

#plt.plot([4.0, 4.0], [0.0, 2.5], linestyle = 'dotted', color = 'blue')
#plt.text(4.2, 2.0, r'\large \textbf{Nyquist limit}')

plt.grid(linestyle = ':', alpha = 0.7)
plt.show()



