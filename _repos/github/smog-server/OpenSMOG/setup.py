from setuptools import setup, find_packages
from os import path

this_dir = path.abspath(path.dirname(__file__))
with open(path.join(this_dir, "README.rst")) as f:
    long_description = f.read()

__version__ = "1.1.0"
for line in open(path.join("OpenSMOG", "__init__.py")):
    if line.startswith("__version__"):
        exec(line.strip())

setup(
    name="OpenSMOG",
    version=__version__,
    description="Structure-based Models for Biomolecules using OpenMM",
    url="https://github.com/junioreif/OpenSMOG",
    author=["Antonio Bento de Oliveira Junior","Vinicius de Godoi Contessoto"],
    author_email="antonio.oliveira@rice.edu,contessoto@rice.edu",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        'Programming Language :: Python :: 3.7',
        "Programming Language :: Python :: 3.8",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Natural Language :: English",
    ],
    include_package_data=True,
    packages=find_packages(),
    install_requires=['numpy', 'lxml'],
    entry_points={"console_scripts": ["CLINAME=OpenSMOG._cli:main"]},
    zip_safe=True,
    long_description=long_description,
)
