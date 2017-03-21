from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension(name = "potential_wrapper",
                sources = ["potential_wrapper.pyx", "ab_initio_potential.c"])


setup(ext_modules = cythonize(ext))
