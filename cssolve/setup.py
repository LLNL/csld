#!/usr/bin/env python3

import glob
import os
import subprocess

#from ez_setup import use_setuptools
#use_setuptools()
from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy


subprocess.call(['make'])
setup(
     name="cssolve",
     packages=find_packages("cssolve"),
     version="0.7",
     install_requires=["numpy>=1.9", "scipy>=0.13", "matplotlib>=1.4"],
     scripts="csfit.py",
     package_data={'cssolve': ['bregman*.so']},
#     license="MIT",
     description="CS solver for CSLD package"
)

#from numpy.distutils.core import setup, Extension
#setup(ext_modules= [Extension('bregman', 
#    sources=['csld/cs_fitting/bregman.f90'],
#    extra_f90_compile_args=["-cpp", "-heap-arrays"])])

