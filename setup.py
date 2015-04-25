from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize(
           "poissondisk.pyx",           # our Cython source
           language="c++",             # generate C++ code
      ))
