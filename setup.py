from distutils.core import setup
from Cython.Build import cythonize


cobj = cythonize("poissondisk.pyx", language="c++")

setup(ext_modules=cobj)
