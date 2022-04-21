#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

# To update the package version number, edit fdpclient/__version__.py
version = {}
with open(os.path.join(here, 'fdpclient', '__version__.py')) as f:
    exec(f.read(), version)

with open('README.md') as readme_file:
    readme = readme_file.read()

setup(
    name='fairdatapoint-client',
    version=version['__version__'],
    description="FAIR Data Point API client",
    long_description=readme + '\n\n',
    long_description_content_type='text/markdown',
    author="Cunliang Geng",
    author_email='c.geng@esciencecenter.nl',
    url='https://github.com/fair-data/fairdatapoint-client',
    packages=find_packages(),
    include_package_data=True,  # check MANIFEST.in
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='fairdatapoint-client',
    classifiers=[ #check details https://pypi.org/classifiers/
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    install_requires=[
        'rdflib',
        'rdflib-jsonld',
        'requests',
    ],
    extras_require={
        'tests': [
            'pytest',
            'pytest-datadir-ng',
            'pytest-cov',
            'coveralls',
            'requests-mock'
        ],
        'docs': ['Sphinx', 'ipython']
    }
)
