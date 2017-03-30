#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

setup(
        cmdclass = {'build_ext': build_ext},
        ext_modules = [Extension
            ("ham",
        sources = ["wrapper.pyx", "hamiltonian.cpp"],
        include_dirs=[numpy.get_include(), "/usr/local/include/eigen3/"],
        language="c++",
        extra_compiles_args=['-O3', '-fPIC', '-march=native', '-mtune=native', '-funroll-loops'],
        )]
)
