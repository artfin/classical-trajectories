from __future__ import print_function

import numpy as np

def initialize_array(length):
    arr = [1.0]
    for s in range(length - 1):
        arr.append(0.0)
    return np.array(arr)
   
def move_array(arr, step):
    temp = arr[:(len(arr)-step)]
    return prezero(temp, len(arr))

def prezero(arr, to_length):
    nzeros = to_length - len(arr)
    
    if nzeros == 0:
        return arr
    else:
        new_arr = np.array([0] * nzeros)
        return np.append(new_arr, arr)
    
def round_levels(arr, grain):
    return np.array([int(round(_ / grain)) for _ in arr])

def step(T, AT, rfreqs, length):

    for freq in rfreqs:
        temp = move_array(AT, freq)
        print(temp)

        T += temp
        print(T)
        print('*'*30)

    return T

maximum = 1000.0 # cm^-1
grain = 100.0 # cm^-1

LENGTH = int(maximum / grain) + 1
print("LENGTH: {0}".format(LENGTH))

T = initialize_array(LENGTH)
AT = initialize_array(LENGTH)

freqs1 = np.array([205.0, 206.0, 420.0, 651.0])
rfreqs1 = round_levels(freqs1, grain)
print('RFREQS: {0}'.format(rfreqs1))

freqs2 = np.array([1., 20., 175., 510.])
rfreqs2 = round_levels(freqs2, grain)
print('RFREQS: {0}'.format(rfreqs2))

T = step(T, AT, rfreqs1, LENGTH)
AT = np.copy(T)

print('-'*50)
T = step(T, AT, rfreqs2, LENGTH)
