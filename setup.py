from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
from os import path, environ


import numpy as np
npgi = np.get_include()

here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md'), 'r') as f:
    long_descr = f.read()

include_dirs = [npgi]
# maybe in a conda env
if 'CONDA_PREFIX' in environ:
    include_dirs.append(path.join(environ['CONDA_PREFIX'], 'include'))

pyeht = Extension(
    name="yaehmop",
    sources=["pyyaehmop/pyeht.pyx"],
    libraries=["yaehmop_eht", "lapack", "blas"],
    #library_dirs=["lib"],
    include_dirs=include_dirs,
    extra_compile_args=['-ffast-math', '-O3'],
)
setup(
    name="pyyaehmop",
    long_description=long_descr,
    long_description_content_type='text/markdown',
    author='Richard J Gowers',
    ext_modules=cythonize([pyeht]),
    packages=find_packages(),
)
