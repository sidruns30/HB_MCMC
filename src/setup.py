from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy
from distutils.sysconfig import get_python_lib
import sys


setup(name="pyHB",version="0.0.3",
    ext_modules=cythonize(Extension("pyHB",
                                    ["pyHB.pyx"],
                                    #["pyHB.pyx","likelihood3.c"],
                                    #include_dirs = [numpy.get_include()],
                                    #language="c",
                                    compiler_directives={'language_level' : "3"}
                                    )
    )
)

