#!/usr/bin/python
"""tangle setup script."""

from setuptools import setup, find_packages

# Package Version
from tangle import __version__ as version

setup(
    name='tangle',
    version=version,

    description=("MBS-PBS 10% dataset utilities"),
    long_description=open('README.md').read(),
    author='Samuele Fiorini, Farshid Hajati, Annalisa Barla, Federico Girosi',
    author_email='samuele.fiorini@dibris.unige.it',
    maintainer='Samuele Fiorini',
    maintainer_email='samuele.fiorini@dibris.unige.it',
    url='https://github.com/samuelefiorini/mbspbs10pc',
    download_url='https://github.com/samuelefiorini/tangle/tarball/' + version,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Programming Language :: Python',
        'License :: OSI Approved :: BSD License',
        'Topic :: Software Development',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS'
    ],
    license='FreeBSD',

    packages=find_packages(),
    install_requires=['pandas (>=0.22.0)',
                      'setuptools (>=39.0.1)',
                      'numpy (>=1.14.2)',
                      'matplotlib (>=2.1.1)',
                      'requests (>=2.18.4)',
                      'Keras (>=2.1.5)',
                      'tensorflow (>=1.6.0)',
                      'tqdm (>=4.19.5)',
                      'beautifulsoup4 (>=4.6.0)',
                      'joblib (>=0.11)',
                      'scikit_learn (>=0.19.1)'
                      ],
    scripts=['scripts/assign_labels.py', 'scripts/cross_validate.py',
             'scripts/extract_sequences.py',
             'scripts/get_population_of_interest.py',
             'scripts/matching_step1.py',
             'scripts/matching_step2.R',
             'scripts/single_train.py'],
)
