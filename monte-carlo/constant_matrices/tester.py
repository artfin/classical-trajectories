import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/N2N2')

from n2n2 import potinit, potn2n2

potinit()

value = potn2n2(10.0, 0.5, 0.5, 0.5)
print "value: {0}".format(value)
