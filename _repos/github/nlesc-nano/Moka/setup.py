#!/usr/bin/env python
import os

from setuptools import setup

HERE = os.path.abspath(os.path.dirname(__file__))

version = {}
with open(os.path.join(HERE, 'ceibacli', '__version__.py')) as f:
    exec(f.read(), version)


with open('README.rst') as readme_file:
    readme = readme_file.read()

setup(
    name='ceibacli',
    version=version['__version__'],
    description="command line interface to compute and query molecular properties from a database",
    long_description=readme + '\n\n',
    author="Felipe Zapata",
    author_email='f.zapata@esciencecenter.nl',
    url='https://github.com/nlesc-nano/ceiba-cli',
    packages=[
        'ceibacli',
    ],
    include_package_data=True,
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords=['Simulations', 'graphql'],
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    entry_points={
        'console_scripts': [
            'ceibacli=ceibacli.cli:main',
            'mock_workflow=ceibacli.mock_workflow:main'
        ]
    },
    install_requires=[
        'numpy', 'pandas', 'python-swiftclient', 'pyyaml>=5.1.1', 'requests',
        'schema', 'typing-extensions'],
    extras_require={
        'test': ['coverage', 'mypy', 'pycodestyle', 'pytest>=3.9', 'pytest-cov',
                 'pytest-mock'],
        'docs': ['sphinx', 'sphinx_rtd_theme']
    },
    data_files=[('citation/ceibacli', ['CITATION.cff'])]
)
