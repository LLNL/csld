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

#try:
#    from numpy.distutils.misc_util import get_numpy_include_dirs
#except ImportError:
#    print("numpy.distutils.misc_util cannot be imported. Attempting to "
#          "install...")
#    subprocess.call(["easy_install", "numpy"])
#    from numpy.distutils.misc_util import get_numpy_include_dirs

subprocess.call(['make'])

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths
data_files = package_files('examples')

setup(
     name="csld",
     author="Fei Zhou",
     author_email="fei.fzhou@gmail.com",
     url="https://to-be-determined",
     platforms=['any'],
     packages=find_packages(),
     version="1.0",
     install_requires=["numpy>=1.9", "scipy>=0.13", "matplotlib>=1.4", "spglib>=1.9"],
     package_data={"csld.util": ["*.json"],
#       'csld': ['../bregman*.so', '../f_phonon*.so', '../bcs_driver*.so','../f_util*.so','../Makefile', '../css*/*.f90','../compile/f_util/*.f90', '../csld/*/*.f90', '../csld/*/*.f']},
       'csld': ['../bregman*.so', '../f_phonon*.so', '../f_util*.so','../Makefile', '../css*/*.f90','../compile/f_util/*.f90', '../csld/*/*.f90', '../csld/*/*.f']},
     data_files=data_files,
     license="MIT",
     description="CSLD",
     long_description="Compressive sensing lattice dynamics",
     cmdclass={'build_ext': build_ext},
     ext_modules=[Extension("_c_util", sources=glob.glob('compile/c_util/*.pyx')+glob.glob('compile/c_util/[a-zA-Z]*.cpp'),
                            extra_compile_args=['-DNDEBUG','-O3','-g0'], extra_link_args=[],
                            include_dirs=[numpy.get_include(), 'compile/c_util'], language='c++', libraries = ['stdc++'])],
     scripts=glob.glob("scripts/*")
)

# from numpy.distutils.core import setup, Extension
# print(dir(Extension))
# setup(ext_modules= [
#     Extension('bregman', sources=['cssolve/bregman.f90'],
#               extra_link_args=["-llapack"],
#               #f2py_options=['--link-lapack_opt'],
#               extra_f90_compile_args=["-cpp", "-heap-arrays"]),
#     Extension('f_phonon', sources=['csld/phonon/f_phonon.f90'],
#               #f2py_options=['--link-lapack_opt'],
#               extra_f90_compile_args=["-cpp", "-heap-arrays"]),
#     Extension('bcs_driver', sources=['cssolve/bcs_driver.f90',
#             'cssolve/num_types.f90', 'cssolve/matrix_sets.f90', 'cssolve/laplace.f90', 'cssolve/bcs.f90'],
#               extra_f90_compile_args=["-cpp", "-heap-arrays"])
#     ])

