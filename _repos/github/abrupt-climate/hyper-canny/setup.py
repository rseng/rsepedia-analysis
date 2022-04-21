#!/usr/bin/env python

from distutils.core import setup
# from distutils.command.build_clib import build_clib
from distutils.extension import Extension
# from setuptools import setup

try:
    from Cython.Build import cythonize
except ImportError:
    has_cython = False
else:
    has_cython = True

from os import path
import os
from codecs import open
from glob import glob
import numpy
from pathlib import Path


here = Path(__file__).parent
src_dir = here.joinpath('src')

c_source_files = [str(p) for p in src_dir.glob("base/*.cc")] \
               + [str(p) for p in src_dir.glob("module/*.cc")]

if os.name == 'nt':
    cflags = ['/O2', '/std:c++17']
else:
    cflags = ['-O2', '-std=c++17']

if has_cython:
    ext_modules = cythonize([Extension(
            "hyper_canny.chc",
            sources=c_source_files + ["hyper_canny/chc.pyx"],
            include_dirs=[numpy.get_include(), './src', './include'],
            extra_compile_args=cflags,
            language="c++")])
else:
    ext_modules = [Extension(
            "hyper_canny.chc",
            sources=c_source_files + ["hyper_canny/chc.cpp"],
            include_dirs=[numpy.get_include(), './src', './include'],
            extra_compile_args=cflags,
            language="c++")]


# Get the long description from the README file

## this would be needed for python2...:
here2=str(here)
with open(here2 + '/README.md', encoding='utf-8') as f:    # this version works for Sebastian
#with open(here / 'README.md', encoding='utf-8') as f:       # version from Johan (fails for Sebastian)
    long_description = f.read()

setup(
    name='HyperCanny',
    version='0.0.1',
    description='High-dimensional (>2) edge finding.',
    author='Johan Hidding',
    author_email='j.hidding@esciencecenter.nl',
    url='https://github.com/abrupt-climate/hyper-canny',
    packages=['hyper_canny'],

    classifiers=[
        'License :: OSI Approved :: Apache Software License',
        'Intended Audience :: Science/Research',
        'Environment :: Console',
        'Development Status :: 4 - Beta',
        'Programming Language :: C++',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering'],

    install_requires=[
            'numpy'
        ],
    extras_require={
        'develop': [
            'pytest', 'pytest-cov', 'pep8', 'pyflakes', 'cython'
        ],
    },
    ext_modules=ext_modules
)
