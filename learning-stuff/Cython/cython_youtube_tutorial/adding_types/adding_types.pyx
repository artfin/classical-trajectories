def var2(list v):
    cdef double m, x
    m = float(sum(v)) / len(v)
    return sum([(x-m)**2 for x in v])/len(v)
