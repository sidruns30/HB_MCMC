from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy
from distutils.sysconfig import get_python_lib
import sys


setup(
    ext_modules=cythonize(["pyHB.pyx","likelihood2.c"],
                          #include_dirs = [numpy.get_include()],
                          language="c",
                          compiler_directives={'language_level' : "3"}
    )
)

