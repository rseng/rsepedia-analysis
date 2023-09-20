import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pythia-uq",
    version="4.0.3",
    author="Nando Hegemann",
    author_email="nando.hegemann@ptb.de",
    description=(
        "Toolbox for non-intrusive functional approximation of data "
        + "via (sparse) general polynomial chaos."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/pythia-uq/pythia",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.5.0",
        "psutil>=5.0",
        "sphinx-autodoc-typehints>=1.18.1",
    ],
)
