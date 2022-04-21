# Copyright 2020-2021 Jason Rumengan
# Copyright 2020-2021 Data61/CSIRO
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ------------------------------------------------------------------------

try:
    from skbuild import setup
except ImportError:
    print('pip must be updated to the latest version for the installer to work.\nRun "pip3 install --user --upgrade pip" to do this.')
    raise
try:
    from psutil import cpu_count
    psutil_found = True
except ImportError:
    psutil_found = False
import os

# Versioning information
# Only bump for stable releases or API breaks
major = ("PYARMA_VERSION_MAJOR", "0")

# Bump for functionality additions that DO NOT break past APIs
minor = ("PYARMA_VERSION_MINOR", "500")

# Bump for bugfixes that DO NOT break past APIs
patch = ("PYARMA_VERSION_PATCH", "3")

name = ("PYARMA_VERSION_NAME", "Violet")

# Check for cores and use the physical cores
if psutil_found:
    os.environ["CMAKE_BUILD_PARALLEL_LEVEL"] = str(cpu_count(logical=False))

# Arguments to pass to CMake
cmake_major = "-D" + major[0] + "=" + major[1]
cmake_minor = "-D" + minor[0] + "=" + minor[1]
cmake_patch = "-D" + patch[0] + "=" + patch[1]
cmake_name = "-D" + name[0] + "=" + name[1]
cmake_args = [cmake_major, cmake_minor, cmake_patch, cmake_name]

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="pyarma",
    version=major[1]+"."+minor[1]+"."+patch[1],
    author="Jason Rumengan",
    author_email="jason.rumengan@connect.qut.edu.au",
    description="Linear algebra library for Python, with emphasis on ease of use",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://pyarma.sourceforge.io",
    project_urls={
        'Documentation': 'https://pyarma.sourceforge.io/docs.html',
        'Source': 'https://gitlab.com/jason-rumengan/pyarma',
        'Tracker': 'https://gitlab.com/jason-rumengan/pyarma/issues'
    },
    packages=["pyarma"],
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Programming Language :: C++",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6.0',
    keywords='linear algebra scientific computing pyarma pyarmadillo armadillo arma c++ pybind11 library',
    cmake_args=cmake_args,
    cmake_install_dir="src/pyarma",
    setup_requires=["setuptools", "wheel",
                      "scikit-build", "cmake", "ninja", "psutil"]
)
