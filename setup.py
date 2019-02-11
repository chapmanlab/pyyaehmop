from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

#import numpy as np
#npgi = np.get_include()

examples_extension = Extension(
    name="yaehmop",
    sources=["pyyaehmop/pyeht.pyx"],
    libraries=["yaehmop_eht", "lapack", "blas"],
    #library_dirs=["lib"],
    #include_dirs=["yaehmop"],
    extra_compile_args=['-ffast-math', '-O3'],
)
setup(
    name="pyeht",
    ext_modules=cythonize([examples_extension])
)
