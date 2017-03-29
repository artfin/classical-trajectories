from statistics import mean
import matplotlib.pyplot as plt

with open('simple_diatomics.dat', 'r') as inputfile:
    lines = inputfile.readlines()

temperatures = []
simple_diatomics = []

for line in lines:
    data = line.split()
    temperatures.append(int(data[0]))
    simple_diatomics.append(float(data[1]))

with open('phase_diatomics.dat', mode = 'r') as inputfile:
    lines = inputfile.readlines()

phase_diatomics = []

for line in lines:
    phase_diatomics.append(float(line.split()[1]))


coeffs = [simple / phase for simple, phase in zip(simple_diatomics, phase_diatomics)]

coeffs = [100 * coeff / mean(coeffs) for coeff in coeffs]

plt.scatter(temperatures, coeffs)
#plt.show()
plt.savefig('deviation.png')


