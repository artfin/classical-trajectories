from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension(name = "derivatives_wrapper",      
                sources = ["derivatives_wrapper.pyx", "ar_co2_pes_deriv.c"])                  

setup(ext_modules=cythonize(ext))
