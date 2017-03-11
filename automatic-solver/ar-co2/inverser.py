import numpy as np

jx = 0.1
jy = -0.0546302489844
jz = 0.1

j = np.sqrt(jx**2 + jy**2 + jz**2)
theta = np.arccos(jz / j)
varphi = np.arccos(jx / j / np.sin(theta))

print 'j: {0}; theta: {1}; varphi: {2}'.format(j, theta, varphi)

jx = j * np.cos(varphi) * np.sin(theta)
jy = j * np.sin(varphi) * np.sin(theta)
jz = j * np.cos(theta)

print 'jx: {0}; jy: {1}; jz: {2}'.format(jx, jy, jz)