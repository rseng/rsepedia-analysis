#!/usr/bin/env python
from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='pyTANSPEC',
      version='0.0.1',
      description='pyTANSPEC package: Python Tool for extracting 1D TANSPEC XD-spectra from 2D image',
      long_description = readme(),
      classifiers=[
          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
          'Programming Language :: Python :: 3.7+',
          'Topic :: Scientific/Engineering :: Astronomy',
      ],
      keywords='Spectrum Extraction For TANSPEC XD-mode',
      url='https://github.com/astrosupriyo/pyTANSPEC',
      author='Supriyo Ghosh',
      author_email='sbuphy2010@gmail.com',
      license='GPLv3+',
      packages=['pyTANSPEC','pyTANSPEC.libs'],
      entry_points = {
          'console_scripts': ['pyxdspec=pyTANSPEC.reduce_specXD:main'],
      },
      install_requires=[
          'numpy',
          'astropy',
          'scikit-image',
          'matplotlib',
          'scipy',
          'ccdproc',
          'SpectrumExtractor @ git+https://github.com/indiajoe/SpectrumExtractor.git@master',
          'WavelengthCalibrationTool @ git+https://github.com/indiajoe/WavelengthCalibrationTool.git@master'
      ],
      include_package_data=True,
      zip_safe=False)
