from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension(name = "hamiltonian_wrapper",
        sources = ["hamiltonian_wrapper.pyx", "hamiltonian.cpp"])

setup(ext_modules = cythonize(ext))

