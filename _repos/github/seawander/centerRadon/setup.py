#!/usr/bin/env python
from setuptools import setup


setup(
    name="radonCenter",
    version="1.0",
    author="Bin Ren",
    author_email="bin.ren@jhu.edu",
    url="https://github.com/seawander/centerRadon",
    py_modules=["radonCenter"],
    description="Center corongraphic images with a Radon transform based routine",
    classifiers=[
        "License :: OSI Approved :: BSD License",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
    install_requires=['numpy', 'scipy']
)
