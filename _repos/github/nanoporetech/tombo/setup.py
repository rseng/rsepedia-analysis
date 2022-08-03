from __future__ import unicode_literals

import os
import sys

import re

from setuptools import setup, Extension
from distutils.version import LooseVersion
from setuptools.command.build_ext import build_ext as _build_ext


# Get the version number from _version.py, and exe_path
verstrline = open(os.path.join('tombo', '_version.py'), 'r').readlines()[-1]
vsre = r"^TOMBO_VERSION = ['\"]([^'\"]*)['\"]"
mo = re.search(vsre, verstrline)
if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError('Unable to find version string in "tombo/_version.py".')

def readme():
    with open('README.rst') as f:
        return f.read()

try:
    import numpy as np
    if LooseVersion(np.__version__) >= LooseVersion("1.20.0"):
        sys.stderr.write(
            '*' * 60 + '\nINSTALLATION ERROR:\n'
            '\tTombo is only compatible with numpy < 1.20. Please downgrade ' +
            'numpy and reinstall tombo `pip install numpy<1.20`\n' +
            '*' * 60 + '\n')
        sys.exit()
    include_dirs = [np.get_include()]
except:
    sys.stderr.write(
        '*' * 60 + '\nINSTALLATION ERROR:\n'
        '\tNeed to install numpy before tombo installation.\n' +
        '\tThis is required in order to get maximum efficincy from ' +
        'cython code optimizations.\nTo install run:\n$ pip install numpy\n' +
        '*' * 60 + '\n')
    sys.exit()

extras_require = ['pyfaidx']
if sys.version_info[0] < 3:
    extras_require.append('rpy2<=2.8.6')
else:
    extras_require.append('rpy2<=2.9.5')

ext_modules = [
    Extension(str("tombo._c_dynamic_programming"),
              [str("tombo/_c_dynamic_programming.pyx")],
              include_dirs=include_dirs,
              language="c++"),
    Extension(str("tombo._c_helper"),
              [str("tombo/_c_helper.pyx")],
              include_dirs=include_dirs,
              language="c++")
]

for e in ext_modules:
    e.cython_directives = {"embedsignature": True}

setup(
    name = "ont-tombo",
    version = __version__,
    packages = ["tombo"],
    install_requires = ['h5py < 3', 'numpy < 1.20', 'scipy', 'cython',
                        'setuptools >= 18.0', 'mappy >= 2.10', 'future', 'tqdm'],
    extras_require={'full':extras_require},

    author = "Marcus Stoiber",
    author_email = "marcus.stoiber@nanoporetech.com",
    description='Analysis of raw nanopore sequencing data.',
    long_description = readme(),
    license = "mpl-2.0",
    keywords = "nanopore high-throughput sequencing correction genome",
    url = "https://github.com/nanoporetech/tombo",

    entry_points={
        'console_scripts': [
            'tombo = tombo.__main__:main'
        ]
    },
    include_package_data=True,
    ext_modules=ext_modules,
    test_suite='nose2.collector.collector',
    tests_require=['nose2'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
)
