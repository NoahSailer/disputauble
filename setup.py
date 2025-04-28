from setuptools import setup
from disputauble import __author__, __version__
import os

file_dir = os.path.abspath(os.path.dirname(__file__))
os.chdir(file_dir)

setup(name="disputauble",
      version=__version__,
      author=__author__,
      license='BSD 2-Clause',
      description='CAMB wrapper for Cobaya to implement effective neutrino masses (i.e. allowing negative masses) via a simple extrapolation.',
      zip_safe=False,  # set to false if you want to easily access bundled package data files
      packages=['disputauble'],
      package_data={'disputauble': ['*.yaml', '*.bibtex', 'yamls/*', 'yamls/**/*', 'defaults/*', 'defaults/**/*']},
      )