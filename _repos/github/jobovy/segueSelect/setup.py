from setuptools import setup #, Extension
import os, os.path
import re

longDescription= ""


setup(name='segueSelect',
      version='1.',
      description='SDSS/SEGUE selection function for G and K stars',
      author='Jo Bovy',
      author_email='bovy@ias.edu',
      license='New BSD',
      long_description=longDescription,
      url='https://github.com/jobovy/segueSelect',
      package_dir = {'segueSelect/': ''},
      packages=['segueSelect'],
      install_requires=['numpy','scipy','pyfits']
      )
