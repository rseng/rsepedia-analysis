# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

exec(open('nanofilt/version.py').read())

setup(
    name='NanoFilt',
    version=__version__,
    description='Filtering and trimming of Oxford Nanopore Sequencing data',
    long_description=open(path.join(here, "README.md")).read(),
    long_description_content_type="text/markdown",
    url='https://github.com/wdecoster/nanofilt',
    author='Wouter De Coster',
    author_email='decosterwouter@gmail.com',
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    keywords='nanopore sequencing processing trimming filtering',
    packages=find_packages() + ['scripts'],
    python_requires='>=3',
    install_requires=['biopython', 'pandas>=0.22.0'],
    package_data={'nanofilt': []},
    package_dir={'nanofilt': 'nanofilt'},
    include_package_data=True,
    entry_points={
        'console_scripts': ['NanoFilt=nanofilt.NanoFilt:main', ]}
)
