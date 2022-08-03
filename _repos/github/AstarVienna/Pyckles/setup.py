"""
Pyckles: A wrapper for common astronomical spectral libraries
=======================================================================

How to compile and put these on pip::

    $ python setup.py sdist bdist_wheel
    $ twine upload dist/*

"""

from setuptools import setup, find_packages

# Version number
with open('pyckles/version.py') as f:
    __version__ = f.readline().split("'")[1]

with open('README.md') as f:
    __readme__ = f.read()

setup(
    name='pyckles',
    version=__version__,
    description="Simple interface to the Pickles 1998 stellar spectra catalogue",
    long_description=__readme__,
    long_description_content_type='text/markdown',
    author='Kieran Leschinski',
    author_email='kieran.leschinski@univie.ac.at',
    url='https://github.com/astronomyk/Pyckles',
    license="GNU General Public License",
    include_package_data=True,
    packages=find_packages(exclude=('tests', 'data')),
    install_requires=['numpy>=1.16', 'astropy', 'matplotlib', 'synphot']
    )
