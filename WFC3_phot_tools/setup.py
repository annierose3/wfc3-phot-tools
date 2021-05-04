#!/usr/bin/env python

from setuptools import find_packages
from setuptools import setup

setup(name = 'WFC3_phot_tools',
      description = 'Python photometry tools for WFC3 calibration',
      author = 'Space Telescope Science Institute',
      url = 'https://github.com/shannnnnyyy/WFC3_phot_tools',
      packages = find_packages(),
      install_requires = ['astropy', 'matplotlib', 'numpy', 'photutils'])