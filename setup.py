from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

#import numpy as np
#npgi = np.get_include()

examples_extension = Extension(
    name="pyeht",
    sources=["pyeht.pyx"],
    libraries=["yaehmop_eht"],
    library_dirs=["tightbind"],
    include_dirs=["tightbind"],
    extra_compile_args=['-ffast-math', '-O3'],
)
setup(
    name="pyeht",
    ext_modules=cythonize([examples_extension])
)
