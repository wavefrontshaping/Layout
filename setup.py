try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

import numpy, sys, os


cmdclass = {}
ext_modules = []

if use_cython:
    ext_modules += [
        Extension("Layout.core",  ["Layout/core.pyx"],include_dirs=[numpy.get_include()]),
        Extension("Layout.layouts",  ["Layout/layouts.pyx"])
    ]
    cmdclass.update({'build_ext': build_ext})
else:
    ext_modules += [
        Extension("Layout.core",  ["Layout/core.c"],include_dirs=[numpy.get_include()]),
        Extension("Layout.layouts",  ["Layout/layouts.c"])
    ]



setup(name='Layout',
    version='0.2',
    setup_requires=[
          'setuptools>=18.0',
          'cython',
      ],
    install_requires=[
          'numpy',
          'matplotlib',
      ],
    packages=['Layout'],
    py_modules=['Layout.logger'],
    cmdclass=cmdclass,
    ext_modules=ext_modules
)
