import numpy as np
from itertools import product

HE_MASS = 4.00260325413;
AR_MASS = 39.9623831237; 
PROTON_TO_ELECTRON_RATIO = 1836.15267389; 

MU = HE_MASS * AR_MASS / ( HE_MASS + AR_MASS ) * PROTON_TO_ELECTRON_RATIO; 

# atomic time unit
ATU = 2.418884326509 * 10**(-17) # s
ALU = 0.52917721067 * 10**(-10) # m

AVU = ALU / ATU
#print("AVU: {0}".format(AVU))

V0_MAX = 4600 / AVU # atomic velocity units
#print('V0_MAX: {0}'.format(V0_MAX))
V0_DELTA = 10 / AVU 
print('V0_STEP: {0}'.format(V0_DELTA))
V0_POINTS = V0_MAX / V0_DELTA + 1

B_MAX = 6.0 * 10**(-10) / ALU 
B_DELTA = 0.25 * 10**(-10) / ALU
print('B_STEP: {0}'.format(B_DELTA))
B_POINTS = B_MAX / B_DELTA + 1

V0 = np.linspace( V0_DELTA, V0_MAX, V0_POINTS )
#print(V0)
B = np.linspace( 0.0, B_MAX, B_POINTS )

R = 40.0

counter = 1 
for b, v0 in product( B, V0 ):
    print("{0} {1} {2} {3}".format( counter, R, b, v0 ))
    counter += 1
