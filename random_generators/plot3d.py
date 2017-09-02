import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

with open('rvectors', 'r') as inputfile:
    lines = inputfile.readlines()

x_comps = []
y_comps = []
z_comps = []

for line in lines:
    ns = line.split()
    x_comps.append(float(ns[0]))
    y_comps.append(float(ns[1]))
    z_comps.append(float(ns[2]))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x_comps, y_comps, z_comps, color = 'k', s=20)
plt.show()
