#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))

# To update the package version number, edit CITATION.cff
with open('CITATION.cff', 'r') as cff:
    for line in cff:
        if 'version:' in line:
            version = line.replace('version:', '').strip().strip('"')

with open('README.md') as readme_file:
    readme = readme_file.read()

setup(
    name='fhirtocapacity',
    version=version,
    description="Conversion of FHIR data and upload to CAPACITY",
    long_description=readme + '\n\n',
    author="Djura Smits",
    author_email='d.smits@esciencecenter.nl',
    url='https://github.com/FAIR-data-for-CAPACITY /fhirtocapacity',
    packages=[
        'fhirtocapacity',
    ],
    include_package_data=True,
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='FHIR-to-CAPACITY',
    scripts=['scripts/fhir_to_capacity', 'scripts/fill_server'],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    test_suite='tests',
    install_requires=[
        'fhirclient~=3.2.0',
        'PyCap~=1.1.3',
        'clize~=4.1.1',
        'Requests~=2.25.1',
        'python-dateutil~=2.8.1'
    ],
    setup_requires=[
        # dependency for `python setup.py test`
        'pytest-runner',
    ],
    tests_require=[
        'pytest',
        'pytest-cov',
        'pycodestyle',
    ],
    extras_require={
        'dev': ['prospector[with_pyroma]', 'yapf', 'isort'],
    },
    data_files=[('citation/fhirtocapacity', ['CITATION.cff'])]
)
