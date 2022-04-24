#!/usr/bin/env python
from setuptools import setup
import shutil

# compilation paper: pandoc paper.md --filter pandoc-citeproc -o paper.pdf
# documentation generation (dans predihood/predihood): pdoc3 --html -o testdoc .
# python3 -m setup bdist_wheel sdist
# python3 -m pip install -e . -r requirements.txt

print("Starting to install predihood and its dependencies.")
print("Note that this installation may take time due to the size of the mongiris dependency (700 MB).")


def delete_dir(directories):  # delete directories (ignore errors such as read-only files)
    for directory in directories:
        shutil.rmtree(directory, ignore_errors=True)


# delete build/config directories because egg-info only updates config files for new versions
delete_dir(["./predihood.egg-info/", "./build/", "./dist/"])

readme = open("README.md").read()

setup(long_description=readme, install_requires=["scikit-learn", "matplotlib", "numpy", "pandas", "requests", "Flask", "setuptools", "area", "seaborn", "StringDist"])
