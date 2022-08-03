from setuptools import setup
from codecs import open
from os import path
import re


here = path.abspath(path.dirname(__file__))
with open(path.join(here, "sslf/__init__.py")) as f:
    contents = f.read()
    version_number = re.search(r"__version__ = \"(\S+)\"", contents).group(1)

with open(path.join(here, "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

setup(name="sslf",
      version=version_number,
      description="A simple spectral line finder",
      long_description=long_description,
      url="https://gitlab.com/chjordan/sslf",
      author="Christopher H. Jordan",
      author_email="christopherjordan87@gmail.com",
      license="MIT",
      keywords="signal processing",
      packages=["sslf"],
      install_requires=["numpy>=1.8.0", "scipy", "future"],
      )
