#!/usr/bin/env python3

# ------------------------------------------------------------------------------
# Copyright (C) 2012-2017 Guillaume Sagnol
# Copyright (C) 2018-2021 Maximilian Stahlberg
#
# This file is part of PICOS Release Scripts.
#
# PICOS Release Scripts are free software: you can redistribute them and/or
# modify them under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# PICOS Release Scripts are distributed in the hope that they will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.
# ------------------------------------------------------------------------------

"""Packaging and installation script for PICOS."""

from pathlib import Path

from setuptools import find_packages, setup
from setuptools.command.sdist import sdist

try:
    from version import get_version
except ModuleNotFoundError:  # PyPI strips version.py from the source package.
    with (Path(__file__).parent / Path("picos", ".version")).open() as file:
        VERSION = file.read().strip()
else:
    VERSION = get_version()


assert len(VERSION.split(".")) in (2, 3)


class sdist_with_hardcoded_version(sdist):
    """Variant of sdist to support the distribution of PICOS via pip.

    This changes the content of picos/.version in the release tree from
    MAJOR.MINOR to MAJOR.MINOR.PATCH so that the PATCH bit is available at
    installation time even though version.py has no access to git metadata.

    See https://gitlab.com/picos-api/picos/-/issues/122.
    """

    def make_release_tree(self, base_dir, files):  # noqa
        sdist.make_release_tree(self, base_dir, files)

        with Path(base_dir, "picos", ".version").open("w") as file:
            file.write(VERSION + "\n")


setup(
    name="PICOS",
    version=VERSION,
    description="A Python interface to conic optimization solvers.",
    long_description=open("README.rst", encoding="utf8").read(),
    long_description_content_type="text/x-rst",
    author="G. Sagnol, M. Stahlberg",
    author_email="incoming+picos-api/picos@incoming.gitlab.com",
    classifiers=[
        "Topic :: Scientific/Engineering :: Mathematics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
    ],
    keywords=[
        "conic optimization",
        "convex optimization",
        "robust optimization",
        "semidefinite programming",
    ],
    project_urls={
        "Source": "https://gitlab.com/picos-api/picos",
        "Documentation": "https://picos-api.gitlab.io/picos/",
    },
    packages=find_packages(include=("picos", "picos.*")),
    install_requires=["cvxopt", "numpy"],
    python_requires=">=3.4",
    package_data={"picos": [".version"]},
    cmdclass={"sdist": sdist_with_hardcoded_version},
)
