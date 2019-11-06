try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

from Cython.Build import cythonize
import numpy, sys, os





ext_modules = [
    Extension("Layout.core",  ["Layout/core.pyx"],include_dirs=[numpy.get_include()]),
    Extension("Layout.layouts",  ["Layout/layouts.pyx"])
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
    py_modules = ['Layout.logger'],
      #find_pyx(path='Layout/')
    ext_modules=cythonize(ext_modules, language_level=sys.version_info[0])
    # [
    #     cythonize(Extension("Layout.core", ["Layout/core.pyx"],
    #     include_dirs=[numpy.get_include()])
	# #   ,library_dirs=['/usr/lib/python2.7']
	# 	),
   	# ]
)
