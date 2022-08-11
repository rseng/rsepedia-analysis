#!/usr/bin/env python

from numpy.distutils.core import setup, Extension

wrapper = Extension('fortran_routines',
                    sources=['ftnroutines.f90'],
                    extra_f90_compile_args=["-std=f2003"])
setup(name='calibrate',
      version='0.1',
      description='DaCapo calibration for CORE',
      author='Maurizio Tomasi',
      author_email='maurizio.tomasi@unimi.it',
      py_modules=['calibrate', 'index'],
      ext_modules=[wrapper]
)
