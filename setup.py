# Run this script: python setup.py build_ext --inplace

from distutils.core import setup
from Cython.Build import cythonize

setup(name='read_file',
      ext_modules=cythonize("file_io.pyx"))
