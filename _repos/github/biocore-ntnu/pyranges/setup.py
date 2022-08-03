import os
from distutils.core import setup
from setuptools import find_packages

__version__ = open("pyranges/version.py").readline().split(" = ")[1].replace(
    '"', '').strip()

install_requires = [
    "cython", "pandas", "ncls>=0.0.63", "tabulate", "sorted_nearest>=0.0.33", "pyrle",
    "natsort"] #,

# optional_requires = ["bamread", "pybigwig", "ray"]

if os.getenv("TRAVIS"):
    install_requires.append("bamread pybigwig".split())

setup(
    name="pyranges",
    packages=find_packages(),
    package_data={
        'pyranges': [
            'example_data/*.bed', 'example_data/*.gtf', 'example_data/*.bam',
            'example_data/*.bam.bai', 'example_data/*.gz', 'example_data/*.bw'
        ]
    },
    # data_files = ["LICENSE"],
    include_dirs=["."],
    version=__version__,
    description="GenomicRanges for Python.",
    author="Endre Bakken Stovner",
    author_email="endrebak85@gmail.com",
    url="http://github.com/endrebak/pyranges",
    keywords=["Bioinformatics"],
    license="MIT",
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta", "Environment :: Other Environment",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        'License :: OSI Approved :: MIT License',
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "Topic :: Scientific/Engineering"
    ],
    long_description=("Performant Pythonic GenomicRanges."))
