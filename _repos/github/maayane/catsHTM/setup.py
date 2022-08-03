from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
#with open(path.join(here, 'README.md'), encoding='utf-8') as f:
#    long_description = f.read()

setup(
    name='catsHTM',
    version='0.2.4',
    description='fast access and cross-matching of large astronomical catalogs',
    #long_description=long_description,
    #long_description_content_type='text/markdown',
    url='https://github.com/maayane/catsHTM',  # Optional
    author='Maayane T. Soumagnac, Eran O. Ofek',
    author_email='maayane.soumagnac@weizmann.ac.il',  # Optional
    classifiers=[ 
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 2.7',
	'Programming Language :: Python :: 3.5',
        'Operating System :: Unix',
        'Operating System :: MacOS',
    ],

    keywords='astronomy catalogs cone-search cross-matching',  # Optional

    packages=["catsHTM"],
    install_requires=['h5py','scipy','tqdm','hdf5storage'],  # Optional
    python_requires='>=2.7.10',

    project_urls={ 
	'Preliminary documentation': 'https://webhome.weizmann.ac.il/home/eofek/matlab/doc/catsHTM.html',	
        'Bug Reports': 'https://github.com/maayane/catsHTM/issues',
        'Matlab Version': 'https://webhome.weizmann.ac.il/home/eofek/matlab/doc/install.html',
        'Credit Page': 'https://webhome.weizmann.ac.il/home/eofek/matlab/doc/catsHTMcredit.html',
        'Source': 'https://github.com/maayane/catsHTM',
    },
)

