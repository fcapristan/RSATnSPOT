from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

setup(
      cmdclass = {'build_ext': build_ext},
      ext_modules = [Extension("serialKernel",
                               sources=["kernelQuad.pyx", "serialKernel.cpp","kernelMath.cpp","sortandmerge.cpp","quadtreeCME.cpp","objectQuad.cpp"],
                               include_dirs=[numpy.get_include()],
                               language="c++",)])
