from setuptools import setup, find_packages
import subprocess

# Run makefile
subprocess.call(['make'])

setup(
    name='EMVC-2',
    version='1.0',
    packages=find_packages(),
    install_requires=[
        'cython>=0.29.17',
        'numpy>=1.16.6,<=1.20.3',
        'argparse>=1.1',
        'scipy>=1.1.0,<1.5.4',
        'tqdm>=4.46.0',
        'scikit-learn>=0.22.2,<=0.24.2',
    ],
)