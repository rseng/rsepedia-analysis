import os
from setuptools import setup
import sys


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def conf_path(name):
    if sys.prefix == '/usr':
        conf_path = os.path.join('/etc', name)
    else:
        conf_path = os.path.join(sys.prefix, 'etc', name)
    return conf_path


setup(
    name="pyfrbcatdb",
    version="2.0.0",
    author="Ronald van Haren, Oscar Martinez-Rubi",
    author_email="r.vanharen@esciencecenter.nl",
    description=("A package for manipulating the frbcatdb and its linking with the VOEvent backbone."),
    license="Apache 2.0",
    keywords="VOEvent, FRBCAT",
    url="https://github.com/TRASAL/frbcatdb",
    packages=['pyfrbcatdb'],
    package_data={'pyfrbcatdb': ['mapping.json', 'zenodo.json']},
    data_files=[(os.path.join(conf_path('pyfrbcatdb')),
                ['pyfrbcatdb/dbase.config'])],
    scripts=['pyfrbcatdb/scripts/decode_VOEvent',
             'pyfrbcatdb/scripts/create_VOEvent',
             'pyfrbcatdb/scripts/frbcatdb-image'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: Apache Software License",
    ],
    install_requires=['voevent-parse', 'python-dateutil',
                      'psycopg2', 'configargparse',
                      'PyYAML', 'astropy', 'requests'],
    setup_requires=['sphinx', 'sphinx-autobuild'],
)
