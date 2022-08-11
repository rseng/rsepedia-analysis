"""Script for installing hera_opm."""
import os
from glob import glob

from setuptools import setup


def package_files(package_dir, subdirectory):
    """Walk the input package_dir/subdirectory to generate a package_data list.

    Parameters
    ----------
    package_dir : str
        Path to the package directory.
    subdirectory : str
        The subdirectory to investigate.

    Returns
    -------
    paths : list of str
        A list containing all of the relevant package_data.
    """
    paths = []
    directory = os.path.join(package_dir, subdirectory)
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            path = path.replace(package_dir + "/", "")
            paths.append(os.path.join(path, filename))
    return paths


data_files = package_files("hera_opm", "data")

setup_args = {
    "name": "hera_opm",
    "author": "HERA Team",
    "url": "https://github.com/HERA-Team/hera_opm",
    "license": "BSD",
    "description": "offline-processing and pipeline managment for HERA data analysis",
    "package_dir": {"hera_opm": "hera_opm"},
    "packages": ["hera_opm"],
    "include_package_data": True,
    "scripts": glob("scripts/*.py") + glob("scripts/*.sh"),
    "use_scm_version": True,
    "package_data": {"hera_opm": data_files},
    "install_requires": ["toml>=0.9.4"],
    "zip_safe": False,
}

if __name__ == "__main__":
    setup(**setup_args)
