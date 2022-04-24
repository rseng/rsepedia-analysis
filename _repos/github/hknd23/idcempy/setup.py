from os import path

from setuptools import find_packages, setup

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="idcempy",
    packages=find_packages(where="idcempy"),
    version="0.1.1",
    license='GPLv3',
    description="Inflated Discrete Choice Models",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Nguyen Huynh, Sergio Bejar, Vineeta Yadav, Bumba Mukherjee",
    author_email="nguyenhuynh831@gmail.com",
    url="https://github.com/hknd23/idcempy",
    keywords=["Inflated", "Mixture", "Ordered Probit", "Multinomial Logit"],
    install_requires=["scipy>=1.4.1", "numpy", "pandas"],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Intended Audience :: Science/Research",
    ],
)
