# Continuous integration
test 1 2 3 
# Known issues

ROOT is dynamically linked against glibc. If you experience errors like the following:

``` root: /lib64/libc.so.6: version `GLIBC_2.{some_old_version}' not found 
(required by /anaconda/envs/testenv/bin/../lib/libstdc++.so.6) ```

or

``` ImportError: anaconda/envs/testenv/lib/libPyROOT.so: ELF file OS ABI invalid ```

this means that you deployed ROOT on a machine with a very old glibc version, and you need to upgrade your distro. 
# Rebuilding Conda binaries

The binaries were built on a Scientific Linux 6.7 (Carbon). The default GCC (4.4.7) is too old for ROOT6, so we used the [Developer Toolset (v2) from CERN](http://linux.web.cern.ch/linux/devtoolset). It provides everything you need to rebuild from the recipes: gcc/binutils/git/etc..
If you really want to rebuild or update the Anaconda binaries (typically not needed as we try keep them updated), you can:
 - use a VM: pickup a VirtualBox image from [here](https://virtualboximages.com/VirtualBox+Scientific+Linux+Images), or create your own based on the ISO images from [here](https://www.scientificlinux.org/downloads/). 
Check [this guide](http://perso.crans.org/~raffo/cern-scientific-linux.php) for instructions on creating a VM. 

 - Much easier: [Docker image](https://hub.docker.com/r/nlesc/slc6-devtoolset-anaconda/) ready with the Developer Toolset (v2), cmake, Anaconda installed and this git repository with the recipes pulled.

To be able to build conda packages, you need to have the conda-build package installed:

```
conda install conda-build
```
From the top directory of this cloned github repository of recipes, building a package is as easy as:
```
conda build ./directory-of-package-recipe
```

**For packages which depend on Python (e.g., ROOT)**, you need to build a separate binary for each Python version you want to provide. Simply pass the CONDA_PY variable to the same command:

```
CONDA_PY={python-version} conda build ./directory-of-package-recipe
```
Remember, you need to build separately each such binary.

Working ROOT has been tested on: Ubuntu 11.10, 12.04, 14.04, 15.04, SLC-6.7, SLC-7. Please try it out and let us know if you experience problems. 

**For packages which depend on ROOT (e.g., root-numpy)**, you need to, in addition, pass the ROOT version against which they should be build. For example:

```
CONDA_PY=3.4 ROOT_VERSION=5.34 conda build ./root-numpy
```
# Installing root_numpy and rootpy via conda

We also provide conda distributions of [root-numpy](https://github.com/rootpy/root_numpy) (the interface between [ROOT](https://root.cern.ch/) and [NumPy](http://www.numpy.org/)) and [rootpy](https://github.com/rootpy/rootpy). **When installing root-numpy, ROOT's latest version will be picked up as a dependency:**

```
$ conda install root-numpy
```
as can be seen from the conda installation plan, the currently-latest (6.04) ROOT version will be picked up:

```
The following NEW packages will be INSTALLED:
    ...
    numexpr:       2.4.6-np110py27_1    defaults                                 
    numpy:         1.10.4-py27_1        defaults                                
    readline:      6.2.5-15             https://conda.binstar.org/NLeSC/linux-64/
    root:          6.04-py2.7_gcc4.8.2  https://conda.binstar.org/NLeSC/linux-64/
    root-numpy:    4.4.0-root6.04_py2.7 https://conda.binstar.org/NLeSC/linux-64/    
    ...           

Proceed ([y]/n)? 

```
If you rather want to have another ROOT version picked up as a dependency, specify that version explicitly:

```
$ conda install root-numpy root=5
```
Now conda proposes:
```
The following NEW packages will be INSTALLED:
    ...
    root:          5.34.32-py2.7_gcc4.8.2  https://conda.binstar.org/NLeSC/linux-64/
    root-numpy:    4.4.0-root5.34.32_py2.7 https://conda.binstar.org/NLeSC/linux-64/
    ...             
    
  ```
 **When installing rootpy, both root-numpy and ROOT will be picked up as dependencies automatically.** The same above holds for fixing your ROOT or Python version.
# Updating the recipes
Building a conda package with conda build involves creating a conda recipe. The recipe is a flat directory holding metadata and the scripts needed to build the package.

The files in a conda recipe are:
* **meta.yaml** (metadata file)
* **build.sh** (Unix build script which is executed using bash)
* **bld.bat** (Windows build script which is executed using cmd)

In our conda github repository, each folder contains the recipe for a particular software package. For example:

![](gitbook2.png)
There could also be other optional files, like patches, scripts for testing, and to optionally execute after installation etc. See more details on the official [conda guidelines page](http://conda.pydata.org/).

```meta.yaml``` (more details [here](http://conda.pydata.org/docs/building/meta-yaml.html)) contains all the metadata: the package name, version, source code repository, build number, dependencies, etc.

```build.sh``` is a script used for building the binary, as specified in the ```meta.yaml``` file. 


---


The details of writing a recipe are beyond the scope of these guidelines, here we only point out some rules of thumb for updating the recipes, especially those for building ROOT.

* **Always increase the build number** in ```meta.yaml``` whenever you change anything in the recipe and you plan on re-building and publishing a new conda package. If you do not update the build number, and publish your new package on the Anaconda Cloud, a user will not pickup your fresh binary. Conda will compare the build number of the user's installation with the one on the Anaconda Cloud, and conclude that there is no new release of your package. 
* **Keep the package name in ```meta.yaml``` consistent** (aka don't change it). This is the name that will be used with ```conda install <package-name>``` so it should typically not contain version information or anything else.
* The ```build.sh``` contains instructions like you would typically build a package/binary outside of the conda environment. For example, for ROOT, it contains:


```
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX \
-Dbuiltin_llvm=ON \
-Dbuiltin-lzma=ON \
-Dbuiltin_zlib=ON \
-Dbuiltin_freetype=ON \
-Dcxx11=ON \
-Drpath=ON \
-Droofit=ON \
-Dopengl=OFF \
-Dgviz=OFF \
|| return 1;
#-DCMAKE_C_COMPILER=$PREFIX/bin/gcc \
#-DCMAKE_CXX_COMPILER=$PREFIX/bin/c++ \
#-Dbuiltin_pcre=ON \

#-Dbuiltin_gsl=ON \

make -j2 || return 1;
make install || return 1;

return 0;
```
**If you want to enable or disable some module of ROOT**, this is where you should do it. Python packages have a simpler ```build.sh``` content, could amount to even one line:

```
 $PYTHON setup.py install || return 1;
```
* Keep the conda package version number intact with the version of the software you are building.
 
Below are the most important parts of the ROOT recipe, what they mean, and how to change them:

![](image4198.png)








# Non-standard channels

Binary packages can have different labels, the default one being `main`. Sometimes binaries are labeled with `dev` (e.g., newer versions which are pending certain tests), which makes them invisible from your anaconda client, unless you explicitly add the corresponding channel to that label, to your configuration:

```
$ conda config --add channels https://conda.anaconda.org/nlesc/label/dev
```
**Beware**: means that when you update or install binaries, you may pickup a newer version of a package, which has not been tested and put in "production" yet. **A safer way** to grab a binary from such a non-standard channel is to rather use directly:
```
conda install -c https://conda.anaconda.org/nlesc/label/dev <package>
```
Upon request, we build binaries with updated versions of ROOT and Python and temporarily put them in the `dev` channel, until properly tested. 
Fo example, ROOT 6.06 has new support for jupyter notebooks, or use JSROOT to implement ROOT graphics for web browsers. Currently in the `dev` channel we provide:

|   | ROOT | ROOT | ROOT |
| ---| ------------- |-------------| -------------:|
| **Python**| 2.7 / 6.05.02 | 3.4 / 6.05.02 | 3.5 / 6.05.02| 

The server can be tested with:
```
serv = new THttpServer("http:8080");
```


[![Build Status](https://api.travis-ci.org/NLeSC/root-conda-recipes.svg)](https://travis-ci.org/NLeSC/root-conda-recipes/) [![DOI](https://zenodo.org/badge/20885/NLeSC/root-conda-recipes.svg)](https://zenodo.org/badge/latestdoi/20885/NLeSC/root-conda-recipes)

Deprecated
==========
This package is no longer being actively maintained. Please switch to the conda-forge package maintained and supported by ROOT users and ROOT core developers:
```sh
conda install root -c conda-forge
```

Below the original contents of this readme.


About
=============
This repository contains Conda recipes for building CERN [ROOT](https://root.cern.ch/) binaries and its dependencies. For the needs of the [XENON Dark Matter project](http://xenon.astro.columbia.edu/), the goal is to provide a "pythonic" interface to the ROOT I/O format, primarily for loading and saving Pandas dataframes in the ROOT format. For this purpose, there are also recipes for building conda binaries of [root-numpy](https://github.com/rootpy/root_numpy) and [rootpy](https://github.com/rootpy/rootpy), the community-driven initiative to provide a more pythonic interface with ROOT on top of the existing PyROOT bindings.

The most most important thing for making things work out of the box is the ABI (binary) compatibility between different gcc(libstdc++)/glibc library versions, on various linux distributions. Typically ROOT would even complain when the GCC headers are not of the same version as the one used for building it, so *shipping the full GCC and glibc as a run dependency* of ROOT, seemed like the best solution.

Combine this with the fact that ROOT 6 requires GCC>=4.8, while we want things to work on older platforms with no "sudo" required, **we decided to fix our GCC distribution to (a relatively recent one) 4.8.2, built against a rather old glibc version 2.12**, making it as cross platform as possible. 

Working ROOT has been tested on: Ubuntu 11.10, 12.04, 14.04, 15.04, SLC-6.7, SLC-7. Please try it out and let us know if you experience problems. 

For details on how to maintain and update recipes, please refer to [this manual](https://www.gitbook.com/book/nlesc/cern-root-conda-recipes/details).

Thanks upfront for any feedback!

Daniela

# Deprecated

This package is no longer being actively maintained. Please switch to [the conda-forge package maintained and supported by ROOT users and ROOT core developers](https://github.com/conda-forge/root-feedstock/):

```sh
conda install root -c conda-forge
```

Below the original contents of this page.

# Installing ROOT via conda

Currently the following ROOT binaries with Python support are provided for the following versions in the `main` channel: 

|   | ROOT | ROOT |
| ---| ------------- |:-------------:| 
| **Python**| 2.7 / 5.34.32 | 3.4 / 5.34.32 |
| **Python** | 2.7 / 6.04  | 3.4 / 6.04 |


To install ROOT in your conda environment, decide upon the ROOT and Python version you plan to use. **We discourage** installing everything in your default conda (*root*) environment, and rather creating a separate one. For example:

```
$ conda create --name=testenv root=6 python=3
$ source activate testenv
```
This will install ROOT6 and Python3 and all dependencies to make things work.

Test if ROOT works like it should:

```
$ root -b -q
   ------------------------------------------------------------
  | Welcome to ROOT 6.05/02                http://root.cern.ch |
  |                               (c) 1995-2014, The ROOT Team |
  | Built for linuxx8664gcc                                    |
  | From tag v6-05-02, 14 September 2015                       |
  | Try '.help', '.demo', '.license', '.credits', '.quit'/'.q' |
   ------------------------------------------------------------

root [0] 
```
```
$ python
Python 3.4.3 |Continuum Analytics, Inc.| (default, Oct 19 2015, 21:52:17) 
Type "help", "copyright", "credits" or "license" for more information.
Anaconda is brought to you by Continuum Analytics.
>>> import ROOT
>>> f = ROOT.TFile.Open("test.root","recreate")
```

If you have other channels in your conda configuration (besides the defaults one), make sure that the following packages are picked up from the right one (NLeSC) when you create the new environment.

You may need to set:
``` 
$ conda config --set show_channel_urls yes 
```
in order to see the channels from which the binaries are installed. 

```
....
    fftw:       3.3.4-0                https://conda.anaconda.org/NLeSC/linux-64/
    gmp:        5.1.2-2                https://conda.anaconda.org/NLeSC/linux-64/
    graphviz:   2.38.0-3               https://conda.anaconda.org/NLeSC/linux-64/
    gsl:        1.16-1                 https://conda.anaconda.org/NLeSC/linux-64/
    isl:        0.12.2-0               https://conda.anaconda.org/NLeSC/linux-64/
    mpc:        1.0.1-0                https://conda.anaconda.org/NLeSC/linux-64/
    mpfr:       3.1.2-0                https://conda.anaconda.org/NLeSC/linux-64/
    ncurses:    5.9-5                  https://conda.anaconda.org/NLeSC/linux-64/
    pcre:       8.35-0                 https://conda.anaconda.org/NLeSC/linux-64/
    gcc:        4.8.2-20               https://conda.anaconda.org/NLeSC/linux-64/
    readline:   6.2.5-11               https://conda.anaconda.org/NLeSC/linux-64/
....
```
# Summary

* [Introduction](README.md)
* [Using the conda binary packages](using_the_conda_binary_packages.md)
   * [Installing ROOT via conda](installing_root_via_conda.md)
   * [Non-standard channels](non-standard_channels.md)
   * [Installing root-numpy and rootpy via conda](installing_rootnumpy_and_rootpy_via_conda.md)
* [Updating the recipes](updating_the_recipes.md)
* [Rebuilding Conda binaries](rebuilding_conda_binaries.md)
* [Continuous integration](continuous_integration.md)
   * Travis
   * Jenkins
* [Known issues](known_issues.md)

# Using the conda binary packages

To use the conda binary packages from the NLeSC AnacondaCloud repository, you need to add the appropriate NLeSC ```main``` channel.  
```
$ conda config --add channels https://conda.anaconda.org/NLeSC
```

Currently the following ROOT binaries with Python support are provided for the following versions in the `main` channel: 

|   | ROOT | ROOT |
| ---| ------------- |:-------------:| 
| **Python**| 2.7 / 5.34.32 | 3.4 / 5.34.32 |
| **Python** | 2.7 / 6.04  | 3.4 / 6.04 |


To install ROOT in your conda environment, decide upon the ROOT and Python version you plan to use. **We discourage** installing everything in your default conda (*root*) environment, and rather creating a separate one. For example:

```
$ conda create --name=testenv root=6 python=3
$ source activate testenv
```
This will install ROOT6 and Python3 and all dependencies to make things work.

Test if ROOT works like it should:

```
$ root -b -q
   ------------------------------------------------------------
  | Welcome to ROOT 6.05/02                http://root.cern.ch |
  |                               (c) 1995-2014, The ROOT Team |
  | Built for linuxx8664gcc                                    |
  | From tag v6-05-02, 14 September 2015                       |
  | Try '.help', '.demo', '.license', '.credits', '.quit'/'.q' |
   ------------------------------------------------------------

root [0] 
```
```
$ python
Python 3.4.3 |Continuum Analytics, Inc.| (default, Oct 19 2015, 21:52:17) 
Type "help", "copyright", "credits" or "license" for more information.
Anaconda is brought to you by Continuum Analytics.
>>> import ROOT
>>> f = ROOT.TFile.Open("test.root","recreate")
```

If you have other channels in your conda configuration (besides the defaults one), make sure that the following packages are picked up from the right one (NLeSC) when you create the new environment.

You may need to set:
``` 
$ conda config --set show_channel_urls yes 
```
in order to see the channels from which the binaries are installed. 

```
....
    fftw:       3.3.4-0                https://conda.anaconda.org/NLeSC/linux-64/
    gmp:        5.1.2-2                https://conda.anaconda.org/NLeSC/linux-64/
    graphviz:   2.38.0-3               https://conda.anaconda.org/NLeSC/linux-64/
    gsl:        1.16-1                 https://conda.anaconda.org/NLeSC/linux-64/
    isl:        0.12.2-0               https://conda.anaconda.org/NLeSC/linux-64/
    mpc:        1.0.1-0                https://conda.anaconda.org/NLeSC/linux-64/
    mpfr:       3.1.2-0                https://conda.anaconda.org/NLeSC/linux-64/
    ncurses:    5.9-5                  https://conda.anaconda.org/NLeSC/linux-64/
    pcre:       8.35-0                 https://conda.anaconda.org/NLeSC/linux-64/
    gcc:        4.8.2-20               https://conda.anaconda.org/NLeSC/linux-64/
    readline:   6.2.5-11               https://conda.anaconda.org/NLeSC/linux-64/
....
```

Binary packages can have different labels, the default one being `main`. Sometimes binaries are labeled with `dev` (e.g., newer versions which are pending certain tests), which makes them invisible from your anaconda client, unless you explicitly add the corresponding channel to that label, to your configuration:

```
$ conda config --add channels https://conda.anaconda.org/nlesc/label/dev
```
**Beware**: means that when you update or install binaries, you may pickup a newer version of a package, which has not been tested and put in "production" yet. **A safer way** to grab a binary from such a non-standard channel is to rather use directly:
```
conda install -c https://conda.anaconda.org/nlesc/label/dev <package>
```
Upon request, we build binaries with updated versions of ROOT and Python and temporarily put them in the `dev` channel, until properly tested. 
Fo example, ROOT 6.06 has new support for jupyter notebooks, or use JSROOT to implement ROOT graphics for web browsers. Currently in the `dev` channel we provide:

|   | ROOT | ROOT | ROOT |
| ---| ------------- |-------------| -------------:|
| **Python**| 2.7 / 6.05.02 | 3.4 / 6.05.02 | 3.5 / 6.05.02| 

The server can be tested with:
```
serv = new THttpServer("http:8080");
```

*Please update your environment regularly, for new and more stable package releases*:

```
$ conda update --all 
$ conda update --yes -q conda
```To build ROOT (same holds for root_numpy):

```$ ROOT_VERSION=[5 or 6] CONDA_PY=[2.7 or 3.4] conda build root/ ```

ROOT 6 requires gcc>=4.8, and the recipes are build with cmake(>=3.2)



 

To build: (from the base conda-recipes directory)
---------------
```
CONDA_PY=x ROOT_VERSION=y conda build root-numpy
```

where ```x``` is the version of Python you want to build against (2.7 and 3.4 for now), and ```y``` is the ROOT version (root-numpy is built against
a specific ROOT version)

For example:
```
CONDA_PY=3.4 ROOT_VERSION=5.34 conda build root-numpy
```
To build: (from the base conda-recipes directory)
---------------
```
CONDA_PY=x conda build root5
```

where ```x``` is the version of Python you want to build against (2.7 and 3.4 for now)

To build: (from the base conda-recipes directory)
---------------
```
CONDA_PY=x conda build root6
```

where ```x``` is the version of Python you want to build against (2.7 and 3.4 for now)
