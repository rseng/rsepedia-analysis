"""
setup.py for the package
"""
from setuptools import setup, find_packages
from codecs import open
from os import path
here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description=f.read()

exec(open("acda/_version.py").read())

setup(
    name='acda',
    packages=find_packages(),
    version=__version__,
    description='Implementation of drug synergy prediction algorithms',
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    author='S. Domanskyi, A. Srivastava',
    author_email='sergii.domanskyi@jax.org',
    license='MIT License',
    url='https://github.com/TheJacksonLaboratory/drug-synergy',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: Unix',
        'Topic :: Education',
        'Topic :: Utilities',
        ],
    python_requires='>=3.8',
    install_requires=[
        'numpy>=1.19.1',
        'pandas>=1.0.1',
        'scipy>=1.4.1',
        'matplotlib>=3.1.3',
        'scikit-learn>=0.22.1'],
    zip_safe=False
)