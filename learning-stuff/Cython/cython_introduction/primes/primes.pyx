# primes.pyx

def primes(int kmax):
    # passing kmax as 'int' type
    # if it's impossible TypeError will be if raised

    cdef int n, k, i # defining local C variables
    cdef int p[10000]
    result = []
    # since the variable hasn't been given a type, it is assumed to hold a Python object

    # the loop will test candidate numbers for primeness until the required number of primes has been found
    # most part of the loop will be translated entirely into C code, and thus runs very fast
    k = 0
    n = 2
    while k < kmax:
        i = 0
        while i < k and n % p[i] != 0:
            i = i + 1
        if i == k:
            p[k] = n
            k = k + 1
            result.append(n)
        n = n + 1
    return result
