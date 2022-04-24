import setuptools

with open("README.md") as fh:
    long_description = fh.read()

setuptools.setup(
    name="circuitgraph",
    version="0.2.0",
    author="Ruben Purdy, Joseph Sweeney",
    author_email="rpurdy@andrew.cmu.edu, joesweeney@cmu.edu",
    description="Tools for working with boolean circuits as graphs.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/circuitgraph/circuitgraph",
    project_urls={
        "Documentation": "https://circuitgraph.github.io/circuitgraph/",
        "Source": "https://github.com/circuitgraph/circuitgraph",
    },
    include_package_data=True,
    packages=["circuitgraph", "circuitgraph.parsing"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=["lark", "networkx",],
)
