#!/usr/bin/env python

import importlib
import os

from setuptools import find_packages, setup

HERE = os.path.abspath(os.path.dirname(__file__))

version = {}
with open(os.path.join(HERE, 'flamingo', '__version__.py')) as f:
    exec(f.read(), version)

with open('README.rst') as readme_file:
    README = readme_file.read()

try:
    importlib.import_module("rdkit")
except ModuleNotFoundError:
    exc = ModuleNotFoundError(
        "'flamingo' requires the 'rdkit' package: https://anaconda.org/conda-forge/rdkit"
    )
    exc.__cause__ = None
    raise exc


setup(
    name='nlesc-flamingo',
    version=version['__version__'],
    description="Compute and filter molecular properties",
    long_description=README + '\n\n',
    long_description_content_type='text/x-rst',
    author="Felipe Zapata",
    author_email='f.zapata@esciencecenter.nl',
    url='https://github.com/https://github.com/nlesc-nano/flamingo',
    packages=find_packages(),
    include_package_data=True,
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='flamingo',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Typing :: Typed',
    ],
    python_requires='>=3.7',
    install_requires=[
        'nlesc-CAT>=0.10.0',
        'nano-CAT>=0.7.0',
        'data-CAT>=0.7.0',
        'plams>=1.5.1',
        'more_itertools',
        'numpy',
        'pandas',
        'pyyaml>=5.1.1',
        'schema!=0.7.5',
        'typing_extensions',
        'h5py',
    ],
    entry_points={
        'console_scripts': [
            'smiles_screener=flamingo.screen:main'
        ]
    },
    package_data={
        'flamingo': [
            'data/scscore/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz',
            'data/scscore/full_reaxys_model_2048bool/model.ckpt-10654.as_numpy.json.gz',
            'py.typed',
        ]
    },
    data_files=[('citation/flamingo', ['CITATION.cff'])],
    extras_require={
        'test': ['coverage', 'mypy', 'pycodestyle', 'pytest>=3.9', 'pytest-cov',
                 'pytest-mock'],
        'doc': ['sphinx', 'sphinx-autodoc-typehints', 'sphinx_rtd_theme',
                'nbsphinx']
    }
)
