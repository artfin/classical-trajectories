import numpy as np
from itertools import product

HE_MASS = 4.00260325413;
AR_MASS = 39.9623831237; 
PROTON_TO_ELECTRON_RATIO = 1836.15267389; 

MU = HE_MASS * AR_MASS / ( HE_MASS + AR_MASS ) * PROTON_TO_ELECTRON_RATIO; 

R = 40.0

B_MIN = 0.0
B_MAX = 12.0

V_MIN = 1.0 / MU
V_MAX = 50.0 / MU

B_STEP = 0.10
V_STEP = 5.000 / MU

b = np.linspace( B_MIN, B_MAX, (B_MAX - B_MIN) / B_STEP + 1 );
#print(b)
v = np.linspace( V_MIN, V_MAX, (V_MAX - V_MIN) / V_STEP + 1 );
#print(v)

counter = 1 
for _b, _v in product( b, v ):
    print("{0} {1} {2} {3}".format(counter, R, _b, - _v))
    counter += 1
