# -*- coding: utf8 -*-
#
# This file were created by Python Boilerplate. Use Python Boilerplate to start
# simple, usable and best-practices compliant Python projects.
#
# Learn more about it at: http://github.com/fabiommendes/python-boilerplate/
#

import os

from setuptools import setup, find_packages

# Meta information
version = open('VERSION').read().strip()
dirname = os.path.dirname(__file__)

# Save version and author to __meta__.py
path = os.path.join(dirname,'tsp',  '__meta__.py')
data = f'''# Automatically created. Please do not edit.
__version__ = '{version}'
__author__ = 'Nick Brown'
'''
with open(path, 'wb') as F:
    F.write(data.encode())

setup(
    # Basic info
    name='tsp',
    version=version,
    author='Nick Brown',
    author_email='nick.brown@carleton.ca',
    url='https://gitlab.com/permafrostnet/teaspoon',
    description='Making permafrost data effortless',
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Software Development :: Libraries',
    ],
    packages = ['tsp', 'tsp/dataloggers', 'tsp/plots'],
    package_data={'tsp': ['dataloggers/test_files/*', 'data/*']},
    
    install_requires=[
        'pandas',
        'numpy',
        'regex'
    ],
    extras_require={
        'dev': [
            'manuel',
            'pytest',
            'pytest-cov',
            'coverage',
            'mock',
        ],
        'nc':[
            'netCDF4',
            'pfit==0.2.1'
        ],
        'plotting':[
            'matplotlib',
            'scipy']
    })

