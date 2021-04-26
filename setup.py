try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

setup(name='SLMlayout',
    version=SLMlayout.__version__,
    setup_requires=[
          'setuptools',
      ],
    install_requires=[
          'numpy',
          'matplotlib',
          'sphinx',
          'sphinx_rtd_theme',
          'numba'
      ],
    packages=['SLMlayout', 
              'SLMlayout.layouts'],
)
