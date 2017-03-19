cdef class DoubleList:
    cdef double* v
    cdef int n

    def __init__(self, v):
        self.n = len(v)
        self.v = <double*> sage_malloc(sizeof(double) * self.n)
        cdef int i

        for i in range(self.n):
            self.v[i] = v[i]
    
    def __del__(self):
        sage_free(self.v)

    def variance(self):
        cdef double m, x, s
        cdef int i, n

        n = self.n
        s = 0

        for i in range(n):
            s += self.v[i]
        m = s/n

        s = 0
        for i in range(n):
            x = self.v[i] - m
            s += x * x
        return s/n

