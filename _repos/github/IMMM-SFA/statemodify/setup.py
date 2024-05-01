"""Setup script for the statemodify package."""

import re

from setuptools import find_packages, setup


def readme():
    """Return the contents of the project README file."""
    with open("README.md") as f:
        return f.read()


version = re.search(
    r"__version__ = ['\"]([^'\"]*)['\"]", open("statemodify/__init__.py").read(), re.M
).group(1)

setup(
    name="statemodify",
    version=version,
    packages=find_packages(),
    url="https://github.com/IMMM-SFA/statemodify",
    license="BSD-2-Clause",
    author="Rohini S. Gupta, Chris R. Vernon",
    author_email="",
    description="A package to modify StateMod's input and output files for exploratory modeling",
    long_description=readme(),
    long_description_content_type="text/markdown",
    python_requires=">=3.8",
    include_package_data=True,
    install_requires=[
        "joblib>=1.1.0",
        "numpy>=1.22.3,<2",
        "pandas>=1.4.2",
        "joblib>=1.1.0",
        "SALib>1.4.5",
        "scipy>=1.8.0",
        "matplotlib>=3.5.1",
        "seaborn>=0.12.2",
        "hmmlearn>=0.2.8",
        "statsmodels>=0.13.5",
        "pyarrow>=10.0.1",
        "SALib>1.4.5",
        "tqdm>=4.64.1",
        "pyyaml>=6.0.0",
        "pygame>=2.3.0",
        "pre-commit",
    ],
)
