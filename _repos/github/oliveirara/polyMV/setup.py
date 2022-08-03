from setuptools import setup
from setuptools import find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="polyMV",
    version="0.1",
    author="Renan Alves de Oliveira",
    author_email="oliveirara@uel.br",
    description="Obtain Multipole Vectors and Frechet Vectors.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/oliveirara/polyMV",
    project_urls={
        "MPSolve": "https://github.com/robol/MPSolve",
        "Planck Legacy Archive": "https://pla.esac.esa.int/",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=["gmpy2", "healpy", "iteration_utilities", "numpy", "numba",],
    py_modules=["polymv", "mpsolve"],
    package_dir={"": "src"},
    keywords="multipole-vectors frechet-vectors cosmology cosmic-microwave-background",
    license="GPLv3"
)
