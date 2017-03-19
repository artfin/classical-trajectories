import vegas
from time import time

#compile cfcn, if needed, at import 
import pyximport
pyximport.install(inplace = True)

import cfcn

def main():
    integ = vegas.Integrator(4 * [[0, 1]])
    
    nevals = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
    nitn = 10

    for neval in nevals:

        start = time()
        result = integ(cfcn.f, neval = neval, nitn = nitn)
        print 'neval: {0}; nitn: {1}; result: {2}'.format(neval, nitn, result.mean)
        print 'time needed: {0}s\n'.format(time() - start)
    
if __name__ == '__main__':
    main()
