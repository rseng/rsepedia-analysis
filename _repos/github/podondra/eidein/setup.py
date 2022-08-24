import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="eidein",
    version="0.0.0",
    author="Ondřej Podsztavek",
    author_email="ondrej.podsztavek@gmail.com",
    description="Interactive tool for dimensionality reduction of data and their visualisation.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/podondra/eidein",
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: 1 - Planning",
        "Environment :: Web Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Visualization"],
    python_requires=">=3.7",
    install_requires=[
        "ipympl>=0.7",
        "ipywidgets>=7.5",
        "matplotlib>=3.3",
        "scikit-learn>=0.23",
        "umap-learn>=0.4"])
