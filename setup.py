
import sys
import os
import shutil

from distutils.core import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

import numpy

# clean previous build
for root, dirs, files in os.walk("./src/main/python/", topdown=False):
    for name in dirs:
        if (name == "build"):
            shutil.rmtree(name)

include_dirs = [
    numpy.get_include(),
    "./include",
    os.environ['SSHT']+"/include/c",
    os.environ['SO3']+"/include/c"
    ]

extra_link_args=[
    "-L./lib",
    "-L"+os.environ['FFTW']+"/lib",
    "-L"+os.environ['SSHT']+"/lib/c",
    "-L"+os.environ['SO3']+"/lib/c"
    ]

setup(
    name = "pys2let",
    version = "2.0",
    prefix='.',
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize([Extension(
        "src/main/python/pys2let",
        package_dir=['src'],
        sources=["src/main/python/pys2let.pyx"],
        include_dirs=include_dirs,
        libraries=["s2let", "so3", "ssht", "fftw3"],
        extra_link_args=extra_link_args,
        extra_compile_args=[]
    )])
)


