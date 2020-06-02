
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
    "/usr/local/opt/src_ssht"+"/include/c",
    "/usr/local/opt/src_so3"+"/include/c"
    ]

extra_link_args=[
    "-L./build",
    "-L"+"/usr/local/lib",
    "-L"+"/usr/local/opt/src_ssht/build/src/c",
    "-L"+"/usr/local/opt/src_so3/build"
    ]

setup(
    name = "pys2let",
    version = "2.0",
    prefix='.',
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize([Extension(
        "pys2let",
        package_dir=['src'],
        sources=["src/main/python/pys2let.pyx"],
        include_dirs=include_dirs,
        libraries=["s2let", "so3", "ssht", "fftw3"],
        extra_link_args=extra_link_args,
        extra_compile_args=[]
    )])
)


