#!/usr/bin/env python3

from setuptools import setup

with open('README.rst', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='camap',
    version='0.1.0',
    author='Tariq Daouda',
    author_email='tariq.daouda@gmail.com',
    url='https://epitopes.world',
    description='Codon Arrangement MAP Predictor, predicting MHC-I Associated Peptides presentation from mRNA',
    long_description=long_description,
    license='MIT',
    package_dir={'camap': 'camap'},
    packages=['camap'],
    #package_data={'camap': ['training_datasets/*']},
    entry_points={'console_scripts': ['camap = camap.__main__:main']},
    install_requires=[
        'numpy==1.18.1',
        'torch==1.3.1',
        'tqdm~=4.43.0',
    ],
    python_requires='>=3.5',
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
)
