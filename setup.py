from setuptools import setup
from Cython.Build import cythonize


cobj = cythonize("poissondisk.pyx", language="c++")

setup(name="poissondisk",
      version="1.0",
      ext_modules=cobj)
