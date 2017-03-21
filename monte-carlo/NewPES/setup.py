from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension(name = "new_pes_cython", 
                sources = ["new_pes_cython.pyx", "ar_co2_pes_deriv.cpp"])

setup(ext_modules=cythonize(ext))
