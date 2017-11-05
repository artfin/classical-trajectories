import numpy as np
from itertools import product

HE_MASS = 4.00260325413;
AR_MASS = 39.9623831237; 
PROTON_TO_ELECTRON_RATIO = 1836.15267389; 

MU = HE_MASS * AR_MASS / ( HE_MASS + AR_MASS ) * PROTON_TO_ELECTRON_RATIO; 

R = 40.0

PR_MIN = -10.0
PR_MAX = -0.05
PR_STEP = 0.05

PT_MIN = -1.0
PT_MAX = 0.0
PT_STEP = 0.05

PR = np.linspace( PR_MIN, PR_MAX, (PR_MAX - PR_MIN) / PR_STEP + 1 );
PT = np.linspace( PT_MIN, PT_MAX, (PT_MAX - PT_MIN) / PT_STEP + 1 );

theta = 0.0

counter = 1 
for pr, pt in product( PR, PT ):
    print("{0} {1} {2} {3} {4}".format(counter, R, pr, theta, pt))
    counter += 1
