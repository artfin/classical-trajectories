from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension(name = "derivatives_cython",      
                sources = ["derivatives_cython.pyx", "derivatives.c"])                  

setup(ext_modules=cythonize(ext))
