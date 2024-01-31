import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "bsavi",
    version = "0.4.0",
    author = "James Wen",
    author_email = "jswen@usc.edu",
    description = ("An interactive visualizer to help explore high-dimensional likelihoods and their observables."),
    license = "MIT",
    keywords = "interactive visualizer cosmology",
    url = "http://packages.python.org/bsavi",
    packages=setuptools.find_packages(where='src'),
    # packages=['bsavi', 'bsavi.cosmo', 'bsavi.loaders'],
    package_dir={'': 'src'},
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires='>=3.8, <3.12',
    install_requires=["holoviews>=1.15.4",
                      "bokeh==2.4.3", # locking bokeh here until v3 works with latex
                      "panel==0.14.4", # same here
                      "spatialpandas==0.4.8",
                      "dask<=2023.5.0", # python 3.8 compatibility
                      "param==1.13.0",
                      "numpy>=1.21, <=1.24",
                      "matplotlib==3.7.1"]
    )

