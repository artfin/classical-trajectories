import numpy as np
from itertools import product
import sys
from math import ceil

HE_MASS = 4.00260325413;
AR_MASS = 39.9623831237; 
PROTON_TO_ELECTRON_RATIO = 1836.15267389; 

MU = HE_MASS * AR_MASS / ( HE_MASS + AR_MASS ) * PROTON_TO_ELECTRON_RATIO; 

# atomic time unit
ATU = 2.418884326509 * 10**(-17) # s
ALU = 0.52917721067 * 10**(-10) # m

AVU = ALU / ATU
#print("AVU: {0}".format(AVU))

V0_MAX = 5655 / AVU # atomic velocity units
#print('V0_MAX: {0}'.format(V0_MAX))
V0_DELTA = 5 / AVU 
#print('V0_STEP: {0}'.format(V0_DELTA))
V0_POINTS = V0_MAX / V0_DELTA + 1

if ( int(ceil(V0_POINTS)) % 2 != 0 ):
    print("V0_POINTS = {0} should be even!".format(V0_POINTS))
    sys.exit()

B_MAX = 6.25 * 10**(-10) / ALU 
B_DELTA = 0.25 * 10**(-10) / ALU
#print('B_STEP: {0}'.format(B_DELTA))
B_POINTS = B_MAX / B_DELTA + 1

if ( int(ceil(B_POINTS)) % 2 != 0 ):
    print("B_POINTS = {0} should be even!".format(B_POINTS))
    sys.exit()

V0 = np.linspace( V0_DELTA, V0_MAX, V0_POINTS )
#print(V0)
B = np.linspace( B_DELTA, B_MAX, B_POINTS )

R = 40.0

def b_weight_simpson( b_index ):  
    if ( b_index == 0 or b_index == (len(B) - 1) ):
        b_weight = 1.0
    if ( b_index % 2 == 1 and b_index != (len(B) - 1) ):
        b_weight = 4.0
    if ( b_index % 2 == 0 and b_index != 0 ):
        b_weight = 2.0

    return b_weight

def v0_weight_simpson( v0_index ):
    if ( v0_index == 0 or v0_index == (len(V0) - 1) ):
        v0_weight = 1.0
    if ( v0_index % 2 == 1 and v0_index != (len(V0) - 1) ):
        v0_weight = 4.0
    if ( v0_index % 2 == 0 and v0_index != 0 ):
        v0_weight = 2.0
    
    return v0_weight

def b_weight_trapezoid( b_index ):
    if ( b_index == 0 or b_index == (len(B) - 1) ):
        b_weight = 1.0
    else:
        b_weight = 2.0
    
    return b_weight

def v0_weight_trapezoid( v0_index ):
    if ( v0_index == 0 or v0_index == (len(V0) - 1) ):
        v0_weight = 1.0
    else:
        v0_weight = 2.0
    
    return v0_weight

if ( sys.argv[1] == "simpson" ):
    integration_type = "simpson"

if ( sys.argv[1] == "trapezoid" ):
    integration_type = "trapezoid"

print('{0}'.format(integration_type))
print("{0} {1}".format(B_DELTA, V0_DELTA))

counter = 1 
for b, v0 in product( B, V0 ):
    b_index = np.where( B == b )[0][0]
    v0_index = np.where( V0 == v0 )[0][0]

    if ( integration_type == "simpson" ):
        b_weight = b_weight_simpson( b_index )
        v0_weight = v0_weight_simpson( v0_index ) 

    if ( integration_type == "trapezoid" ):
        b_weight = b_weight_trapezoid( b_index )
        v0_weight = v0_weight_trapezoid( v0_index )

    print("{0} {1} {2} {3} {4} {5}".format( counter, R, b, b_weight, v0, v0_weight ))
    counter += 1
