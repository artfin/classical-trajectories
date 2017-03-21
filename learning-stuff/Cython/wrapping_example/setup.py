from distutils.core import setup
from Cython.Build import cythonize

setup(name = "get_len",
        ext_modules = cythonize("len_extern.pyx"))
