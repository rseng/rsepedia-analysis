import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
    name='casbert',
    version='0.0.1',
    author="Yuda Munarko",
    author_email="yuda.munarko@gmail.com",
    description="A tool for cellml, variables, and sedml search. This tool using BERT to represent objects.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    install_requires=[
        'requests',
        'lxml',
        'ipython',
        'torch',
        'sentence-transformers',
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
    package_data={'': ['resources/*', 'sedmlImages/*']},
)
