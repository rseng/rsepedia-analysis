from setuptools import setup, find_packages

import os
import sys
import io
import re


def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    os.system("python setup.py bdist_wheel upload")
    sys.exit()

setup(
    name="mrpy",
    version=find_version("mrpy", "__init__.py"),
    packages=find_packages(),#['mrpy','mrpy.fitting','mrpy.base','mrpy.extra','mrpy.fitting.stan'],
    install_requires=["numpy>=1.6.2",
                      "scipy>=0.12.0",
                      "mpmath",
                      "cached_property"],
    # package_data= {"mrpy.fitting.stan":['mrpy/fitting/stan/stan_models/*']},
    #scripts=[],
    author="Steven Murray",
    author_email="steven.murray@curtin.edu.au",
    description="An efficient alternative halo mass function distribution",
    long_description=read('README.rst'),
    license="MIT",
    keywords="halo mass function bayesian gamma",
    url="https://github.com/steven-murray/mrpy",
    classifiers=["Development Status :: 4 - Beta",
                 "Programming Language :: Python :: 2",
                 "Programming Language :: Python :: 2.7",
                 "License :: OSI Approved :: MIT License",
                 "Intended Audience :: Science/Research",
                 "Natural Language :: English",
                 "Topic :: Scientific/Engineering :: Mathematics",
                 "Topic :: Scientific/Engineering :: Physics"]
    # could also include long_description, download_url, classifiers, etc.
)
