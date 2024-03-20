"""
Copyright(C) 2023 by
Trey V. Wenger; tvwenger@gmail.com

GNU General Public License v3 (GNU GPLv3)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from setuptools import setup

setup(
    name="kinematic_scaleheight",
    version="2.2",
    description="Kinematic estimates of vertical cloud distribution",
    author="Trey V. Wenger",
    author_email="tvwenger@gmail.com",
    packages=["kinematic_scaleheight"],
    install_requires=["numpy", "scipy", "matplotlib", "pymc==5.8.2", "corner"],
    package_data={"kinematic_scaleheight": ["data/reid19_corr.dat"]},
)
