from __future__ import print_function
import numpy as np

normalize = lambda v: v / np.linalg.norm(v)

A = normalize(np.array([0.043477, 0.036412, 0.998391], dtype = 'float'))
B = normalize(np.array([0.60958, 0.73540, 0.29597], dtype = 'float'))

product = np.dot(A, B)
cross = np.cross(A, B)
print('cross product: {0}; norm: {1}'.format(cross, np.linalg.norm(cross)))
print('scalar product: {0}'.format(product))

G = np.array([[product, -np.linalg.norm(cross), 0.0], [np.linalg.norm(cross), product, 0.0], [0.0, 0.0, 1.0]], dtype = 'float')
print('G:\n {0}'.format(G))

F = np.array([A, normalize(B - product * A), cross], dtype = 'float')
print('F:\n {0}'.format(F))

U = F * G * np.linalg.inv(F)
print('U:\n {0}'.format(U))

test_B = U.dot(A)
print('test_B: {0}'.format(test_B))
print('B: {0}'.format(B))

