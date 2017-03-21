cdef extern from "string.h":
    # describe the interface for the functions used
    int strlen(char*c )

def get_len(char* message):
    # strlen can now be used from cython code
    return strlen(message)
