def var3(list v):
    cdef double m, x, s
    cdef int n

    n = len(v)
    
    s = 0
    for x in v:
        s += x
    m = s/n
    
    s = 0
    for x in v:
        s += (x - m) * (x - m)
    return s / n
