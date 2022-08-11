# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 10:48:20 2021

@author: irawo
"""

import setuptools
#import py2exe

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="INSANE",
    version="0.0.1",
    author="Ira_Wolfson",
    author_email="irawolfsonprof@gmail.com",
    description="INflationary potential Simulator and ANalysis Engine",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/beastraban/INSANE",
    project_urls={
        "Bug Tracker": "https://github.com/beastraban/INSANE/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    entry_points={'console_scripts':['INSANE=INSANE.Usage_Example:main']},
    # console=['Usage_Example.py'],
    # options={
    #           "py2exe":{
    #                 "optimize": 2,
    #                 "includes": ["MsSolver16.py", "BackgroundGeometry.py"], # List of all the modules you want to import
    #                 "packages": ["INSANE"] # List of the package you want to make sure that will be imported
    #                 }
    #           }
    )