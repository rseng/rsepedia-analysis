import glob
import os

from setuptools import find_packages, setup

# Read README for PyPI, fallback if it fails.
desc = "Tools for creating and manipulating shapes."
try:
    readme_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "README.rst")
    with open(readme_file) as f:
        readme = f.read()
except ImportError:
    readme = desc

version = "0.6.1"


################################################
# Set up for the various optional dependencies
# that may be installed for additional features.
################################################

test_deps = [
    "hypothesis[numpy]",
    "matplotlib",
    "plato-draw",
    "pytest",
    "pytest-cov",
]

bounding_deps = [
    "miniball",
]

extras = {
    "test": test_deps + bounding_deps,
    "bounding_sphere": bounding_deps,
}

# Acquire package data files.
DATA = [fn.replace("coxeter/", "") for fn in glob.glob("coxeter/families/data/*.json")]

setup(
    name="coxeter",
    version=version,
    description=desc,
    long_description=readme,
    long_description_content_type="text/x-rst",
    url="https://github.com/glotzerlab/coxeter",
    author="Vyas Ramasubramani",
    author_email="vramasub@umich.edu",
    packages=find_packages(),
    package_data={"coxeter": DATA},
    install_requires=["numpy>=1.19.0", "rowan>=1.2.0", "scipy>=1.0.0"],
    tests_require=test_deps,
    extras_require=extras,
    zip_safe=False,
)
