from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
from os import path, environ

cython_linetrace = bool(environ.get('CYTHON_TRACE_NOGIL', False))


here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md'), 'r') as f:
    long_descr = f.read()

include_dirs = []
# maybe in a conda env, so add this to path
if 'CONDA_PREFIX' in environ:
    include_dirs.append(path.join(environ['CONDA_PREFIX'], 'include'))

compile_args = ['-ffast-math', '-O3']
if cython_linetrace:
    compile_args.append("-DCYTHON_TRACE_NOGIL")


pyeht = Extension(
    name="pyyaehmop._pyeht",
    sources=["pyyaehmop/pyeht.pyx"],
    libraries=["yaehmop_eht", "lapack", "blas"],
    include_dirs=include_dirs,
    extra_compile_args=compile_args,
)
setup(
    name="pyyaehmop",
    long_description=long_descr,
    long_description_content_type='text/markdown',
    author='Richard J Gowers',
    ext_modules=cythonize([pyeht],
                          compiler_directives={
                              'linetrace': cython_linetrace,
                          },
    ),
    packages=find_packages(),
)
