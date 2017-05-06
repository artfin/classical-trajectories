#!/usr/bin/env python
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
        cmdclass = {'build_ext': build_ext},
        ext_modules = [Extension
            ("ham",
        sources = ["wrapper.pyx", "hamiltonian.cpp"],
        include_dirs=["/usr/local/include/eigen3/"],
        language="c++",
        extra_compile_args=['-std=c++11', '-O3', '-fPIC', '-march=native', '-mtune=native', '-funroll-loops'],
        )]
)

