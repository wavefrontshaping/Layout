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
        Extension("SLMlayout.core",  ["SLMlayout/core.pyx"],include_dirs=[numpy.get_include()]),
    ]
    cmdclass.update({'build_ext': build_ext})
else:
    ext_modules += [
        Extension("SLMlayout.core",  ["SLMlayout/core.c"],include_dirs=[numpy.get_include()]),
    ]

setup(name='SLMlayout',
    version='0.2',
    setup_requires=[
          'setuptools>=18.0',
          'cython',
      ],
    install_requires=[
          'numpy',
          'matplotlib',
          'sphinx',
          'sphinx_rtd_theme',
      ],
    packages=['SLMlayout', 'SLMlayout.layouts'],
    cmdclass=cmdclass,
    ext_modules=ext_modules
)
