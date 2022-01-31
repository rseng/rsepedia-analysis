![Pyntacle logo](http://pyntacle.css-mendel.it/images/title_joined.png)

A Python package for network analysis based on non canonical
metrics and HPC-Computing

- **Compatibility**: Python 3.7
- **Contributions**: bioinformatics@css-mendel.it
- **Website**: http://pyntacle.css-mendel.it
- **Pypi**: https://pypi.org/project/pyntacle/
- **Conda**: https://anaconda.org/bfxcss/pyntacle [![Anaconda-Server Badge](https://anaconda.org/bfxcss/pyntacle/badges/platforms.svg)](https://anaconda.org/bfxcss/pyntacle) [![Anaconda-Server Badge](https://anaconda.org/bfxcss/pyntacle/badges/downloads.svg)](https://anaconda.org/bfxcss/pyntacle)
- **Docker Hub**: https://hub.docker.com/r/mazzalab/pyntacle
- **Bug report**: https://github.com/mazzalab/pyntacle/issues


## Installing using Pypi
[![PyPI version](https://badge.fury.io/py/pyntacle.svg)](https://badge.fury.io/py/pyntacle)

### [optional] Create and activate a virtualenv
#### Linux and MacOS X
```bash
python -m venv pyntacle_env
source pyntacle_env/bin/activate
```
#### Windows
```bash
python -m venv pyntacle_env
.\pyntacle_env\Scripts\activate
```
### Installation
```bash
pip install pyntacle
```


## Installing using Anaconda or Miniconda [![Anaconda-Server Badge](https://anaconda.org/bfxcss/pyntacle/badges/installer/conda.svg)](https://conda.anaconda.org/bfxcss)
![Anaconda-Server Badge](https://anaconda.org/bfxcss/pyntacle/badges/version.svg) ![Anaconda-Server Badge](https://anaconda.org/bfxcss/pyntacle/badges/latest_release_date.svg)

There are several advantages in using Anaconda to
install not only Pyntacle, but also Python and other packages: it is
cross platform (Linux, MacOS X, Windows), you do not require
administrative rights to install it (it goes in the user home
directory), it allows you to work in virtual environments, which can be
used as safe sandbox-like sub-systems that can be created, used,
exported or deleted at your will.

You can choose between the full [Anaconda](http://docs.continuum.io/anaconda/) and its lite version,
[Miniconda](http://conda.pydata.org/miniconda.html). The difference between the two is that Anaconda comes with
hundreds of packages and can be a bit heavier to install, while
Miniconda allows you to create a minimal, self-contained Python
installation, and then use the [Conda](https://conda.io/docs/) command to install additional
packages of your choice.

In any case, Conda is the package manager that the Anaconda and
Miniconda distributions are built upon. It is both cross-platform and
language agnostic (it can play a similar role to a pip and virtualenv
combination), and you need to set it up by running either the [Anaconda
installer](https://www.anaconda.com/download/) or the
[Miniconda installer](https://conda.io/miniconda.html), choosing the
Python 3.7 version.

The next step is to create a new Conda environment (if you are familiar
with virtual environments, this is analogous to a virtualenv).

#### Linux and MacOS X

Run the following commands from a terminal window:

```bash
conda create -n name_of_my_env python=3.7
```

This will create a minimal environment with only Python 3.7 installed
in it. To put your self inside this environment run:

```bash
source activate name_of_my_env
```

And finally, install the latest version of Pyntacle:

```bash
conda install -y -c bfxcss -c conda-forge pyntacle
```

#### Windows
<aside class="warning">
<b>Warning</b>: Windows users could experience some issues when installing Conda or Miniconda in folders that
whose name contains whitespaces (e.g. "%userprofile%\John Doe\Miniconda"). This is a known bug,
as reported <a href="https://github.com/ContinuumIO/anaconda-issues/issues/1029" target="_blank">here</a> and
<a href="https://groups.google.com/a/continuum.io/forum/#!topic/anaconda/zTQQ0NqqIvk" target="_blank">here</a>. If this occurs,
we recommend to create a new directory with no whitespaces (e.g. "%userprofile%\John<b style="color:red">_</b>Doe\) and install Conda in there.
</aside>


Open a windows prompt or (even better) an
[Anaconda prompt](https://chrisconlan.com/wp-content/uploads/2017/05/anaconda_prompt.png)
, and type:

```bash
conda create -y -n name_of_my_env python=3.7
```

Then, activate the newly created environment:

```bash
conda activate name_of_my_env
```

Finally, install the latest version of Pyntacle:

```bash
conda install -y -c bfxcss -c conda-forge pyntacle
```

### CUDA support (experimental)

Independently of the OS in use, if you need CUDA support, you must
also install the CUDA toolkit by downloading and installing the Toolkit from the
[_NVIDIA website_](https://developer.nvidia.com/cuda-toolkit).

**NOTE** GPU-base processing is an **experimental** feature in the current version (1.3), and is not covered by the command-line interface. This is because of weird behaviors of Numba with some hardware configurations that we were not be able to describe and circumvent so far. Although currently accessible by APIs, the GPU feature will be stable in the release 2.0, when Pyntacle will have covered the possibility to manage huge matrices for which replacing fine-grained parallelism with GPU computing would make sense.



## Release history

Changelog for current and past releases:

### 1.3.2:
Bug fixes:
- \#55 ImportAttributes.import_node_attributes bad attribute

### 1.3.1:
Bug fixes:
- \#47 -nprocs removed in keyplayer kp-info command line
- \#48 empty result with kp-finder
- \#49 -seed argument removed in pyntacle generator
- \#50 --plot-format option removed from 1.2 onward
- \#51 seed argument removed in group degree API
- \#52 bad handling of missing output file names
- \#53 bad handling of empty set due to graph intersection

### 1.3:
Major updates:
- [algorithms] Implementation of the new Stochastic Gradient Descent (SGD) search algorithm
- [tests] Tests for SGD included
- [environment] Upgraded base Python version to 3.7
- [environment] Install igraph ver. 0.8.2 (conda-forge)

Minor updates:
- removed dependency to Cairo and the old plotter

### 1.2:
Major updates:
- [command-line] The algorithm that decides the computing configuration to be used to analyze a give graph was updated to exclude the possibility to run multi-process and multi-threaded at the same time. This is still possible by accessing directly to the APIs.
- [command-line] Renamed option from -T/--threads to -O/--nprocs to avoid clashes with other synonymous options
- [API] Removed all decorator methods that over-checked the sanity of the arguments of methods. These resulted to improve.
- [PyntacleInk] bug #28 "initial value" and "value" are swapped, solved
- [Tests] bug #25 "gr-finder bruteforce test fails", solved

Minor:
- [command-line] bug #23 "the command line option --type m-reach in kp-finder produces no output", solved
- [API] removed the *max_distance* argument from all methods
- [API] removed the seed from each methods. Postponed to later versions the implementation of clever manner of controlling randomness of number generators
- the default number of forked processes is now 1 and not equals to the total number of available processors -1
- removed *shortest_path_modifications.py* file

### 1.1:
New Graph Plotting tool: PyntacleInk
- PyntacleInk is a web-based, javascript-based visualizer based on the [sigmajs](http://sigmajs.org/) library. It is 
designed to integrate with Pyntacle command-line library and replaces the igraph-based plotter. 
- Graphical plots will be now be produced in a html file within the pyntacle results directory, containing detailed 
summarization of each Pyntacle run performed on the input file. This file is updated at each run, making graphical
expolration of results more intuitive.

### 1.0:

Major update of Pyntacle, including:
- New major feature: Group centralities. A redesign of single-node centralities that accounts for the centrality of groups have been added to Pyntacle and a new command, groupcentrality has now been added to the Pyntacle command line. Its behavior is similar to the keyplayer command. Users can compute group centrality indices for predetermined sets of nodes or perform group centrality-based searches for sets of nodes that optimize predetermined group centrality scores
- Octopus redesign
- All graphs imported from either the command line or through the `io_stream` methods now have the same id-->node name pairing, with the exceptions of pickled igraphs (binary)
- Node isolates are now removed from graphs when imported through any Pyntacle import istance
- GPU-based computation of the shortest paths using the Floyd-Warshall algorithm is now an experimental feature and is disabvled in the Pyntacle command line. Users can choose to override this behavior in the Pyntacle library by using the correct Cmode enumerator
- Added exceptions and specific behaviors to unusual group-based search calculations that caused exceptions before
- Minor bugfixes 

### 0.2:

- Added the --input-separator option to all commands.
- Added the --repeat option to pyntacle generate, so that the user can decide how many random graphs need to be created in one run.
- Bugfix in the edgelist importer, when a header is present.
- Bugfix for edgelist when node names are numbers, and now whitelines are skipped.
- Communities gracefully exits when no get_modules are found or all get_modules are filtered out by the user's custom filters.
- Major editing of the main inline help to match the documentation on the website.
- Added warnings in documentation for Windows users that have whitespaces in the Conda installation folder.
- Minor bugfixes

### 0.1.3:

- Bugfixes

### 0.1.2:

- Bugfixes

### 0.1.1:

-  The first release of Pyntacle.



## License

This work is licensed under a [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
pyntacle.algorithms.local\_topology module
==========================================

.. automodule:: pyntacle.algorithms.local_topology
    :members:
    :undoc-members:
    :show-inheritance:
pyntacle.tools.edgelist\_utils module
=====================================

.. automodule:: pyntacle.tools.edgelist_utils
    :members:
    :undoc-members:
    :show-inheritance:
pyntacle.io\_stream.exporter module
===================================

.. automodule:: pyntacle.io_stream.exporter
    :members:
    :undoc-members:
    :show-inheritance:
pyntacle.graph\_operations.set\_operations module
=================================================

.. automodule:: pyntacle.graph_operations.set_operations
    :members: GraphOperations
    :undoc-members:
    :show-inheritance:
pyntacle.algorithms.keyplayer module
====================================

.. automodule:: pyntacle.algorithms.keyplayer
    :members:
    :undoc-members:
    :show-inheritance:
pyntacle.algorithms.bruteforce\_search module
=============================================

.. automodule:: pyntacle.algorithms.bruteforce_search
    :members: BruteforceSearch
    :undoc-members:
    :show-inheritance:
pyntacle.tools package
======================

.. automodule:: pyntacle.tools
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

.. toctree::

   pyntacle.tools.add_attributes
   pyntacle.tools.adjmatrix_utils
   pyntacle.tools.edgelist_utils
   pyntacle.tools.enums
   pyntacle.tools.graph_utils
   pyntacle.tools.octopus

pyntacle.tools.add\_attributes module
=====================================

.. automodule:: pyntacle.tools.add_attributes
    :members:
    :undoc-members:
    :show-inheritance:
pyntacle.graph\_operations package
==================================

.. automodule:: pyntacle.graph_operations
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

.. toctree::

   pyntacle.graph_operations.communities
   pyntacle.graph_operations.set_operations

pyntacle.tools.graph\_utils module
==================================

.. automodule:: pyntacle.tools.graph_utils
    :members:
    :undoc-members:
    :show-inheritance:
pyntacle.graph\_operations.get_modules\_finder module
=================================================

.. automodule:: pyntacle.graph_operations.communities
    :members: CommunityFinder
    :undoc-members:
    :show-inheritance:
pyntacle.tools.enums module
===========================

.. automodule:: pyntacle.tools.enums
    :members: CmodeEnum, KpnegEnum, KpposEnum, GroupCentralityEnum, GroupDistanceEnum, GraphOperationEnum
    :undoc-members:
    :show-inheritance:
pyntacle package
================

.. automodule:: pyntacle
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
-----------

.. toctree::

    pyntacle.algorithms
    pyntacle.graph_operations
    pyntacle.io_stream
    pyntacle.tools

pyntacle.io\_stream package
===========================

.. automodule:: pyntacle.io_stream
    :members:
    :undoc-members:
    :show-inheritance:

Submodules
----------

.. toctree::

   pyntacle.io_stream.converter
   pyntacle.io_stream.export_attributes
   pyntacle.io_stream.exporter
   pyntacle.io_stream.generator
   pyntacle.io_stream.import_attributes
   pyntacle.io_stream.importer

pyntacle.io\_stream.export\_attributes module
=============================================

.. automodule:: pyntacle.io_stream.export_attributes
    :members:
    :undoc-members:
    :show-inheritance:
pyntacle.algorithms.scalefree\_inference module
===============================================

.. automodule:: pyntacle.algorithms.scalefree_inference
    :members:
    :undoc-members:
    :show-inheritance:
pyntacle.algorithms package
===========================

.. automodule:: pyntacle.algorithms
    :members:
    :undoc-members:
    :show-inheritance: False
	:exclude-members:
Submodules
----------

.. toctree::

   pyntacle.algorithms.bruteforce_search
   pyntacle.algorithms.global_topology
   pyntacle.algorithms.greedy_optimization
   pyntacle.algorithms.keyplayer
   pyntacle.algorithms.local_topology
   pyntacle.algorithms.scalefree_inference
   pyntacle.algorithms.shortest_path
   pyntacle.algorithms.sparseness

pyntacle.algorithms.shortest\_path module
=========================================

.. automodule:: pyntacle.algorithms.shortest_path
    :members:
    :undoc-members:
    :show-inheritance:
	:exclude-members: ShortestPath.subtract_count_dist_matrix
pyntacle.algorithms.greedy\_optimization module
===============================================

.. automodule:: pyntacle.algorithms.greedy_optimization
    :members:
    :undoc-members:
    :show-inheritance:
pyntacle.io\_stream.generator module
====================================

.. automodule:: pyntacle.io_stream.generator
    :members:
    :undoc-members:
    :show-inheritance:
############
Installation
############

The easiest way for the majority of users to install pandas is to install it as part of the Anaconda distribution, a cross-platform distribution for data analysis and scientific computing. This is the recommended installation method for most users.

Instructions for installing from source on various Linux distributions and MacOs are also provided.

**********************
Python version support
**********************
Officially, Python **>= 3.5**.

*******************
Installing Pyntacle
*******************

===============================================
Installing Pyntacle using Anaconda or Miniconda
===============================================
Installing Pyntacle and all its dependencies can be challenging for inexperienced users.
There are several advantages in using Anaconda to install not only Pyntacle, but also Python and other packages: it is cross platform (Linux, MacOS X, Windows),
you do not require administrative rights to install it (it goes in the user's home directory),
it allows you to work in *virtual environments*, which can be used as safe sandbox-like sub-systems that can be created, used, exported or deleted at your will.

You can choose between the full `Anaconda <http://docs.continuum.io/anaconda/>`_ and its minified version, `Miniconda <http://conda.pydata.org/miniconda.html>`_.
The difference between the two is that Anaconda comes with hundreds of packages and can be a bit heavier to install,
while Miniconda allows you to create a minimal, self-contained Python installation, and then use the `Conda <https://conda.io/docs/>`_ command to install additional packages of your choice.

In any case, `Conda <https://conda.io/docs/>`_ is the package manager that the Anaconda and Miniconda distributions are built upon.
It is both cross-platform and language agnostic (it can play a similar role to a pip and virtualenv combination), and you need to set it up by running either the `Anaconda installer <https://www.anaconda.com/download/>`_
or the `Miniconda installer <https://conda.io/miniconda.html>`_, choosing the **Python 3.6** version.

The next step is to create a conda environment specifically for Pyntacle (if you are familiar with virtual environments, these are analogous to a virtualenv but they also allow you to specify precisely which Python version to install).
To do so, you need to download our pre-compiled **environment file**, that contain all the dependencies and their versions required to run the tool:

	* for Unix based operative system (Mac and Linux), `you can use this environment  <http://pyntacle.css-mendel.it/resources/envs/unix/pyntacle_latest.yml>`_
	* for Windows, you will need `this one <http://pyntacle.css-mendel.it/resources/envs/win/pyntacle_latest.yml>`_ 

.. note:: To date (November 2018) the environments are identical, but future releases of Pyntacle may include different dependencies for each operative system

------------------
Linux and Mac OS X
------------------

Run the following commands from a terminal window:
::

   conda env create -n pyntacle_env --file pyntacle_latest.yml

This will create a Pyntacle-ready environment running Python v.3.6. You can enter in this environment as follow:

::

  source activate name_of_my_env

And run Pyntacle freely.

-------
Windows
-------

Open a cmd terminal window or - better - an `Anaconda Prompt <https://chrisconlan.com/wp-content/uploads/2017/05/anaconda_prompt.png>`_, and type:

::

  conda env create -n pyntacle_env --file pyntacle_latest.yml

Then, we activate the newly created environment:

::

  activate name_of_my_env

.. warning:: Windows users could experience some issues when installing Conda or Miniconda in folders that have a whitespace in their name (e.g. ``C:\John Doe\Miniconda``). This is a known bug, as reported `here <https://github.com/ContinuumIO/anaconda-issues/issues/1029>`_ and `here <https://groups.google.com/a/continuum.io/forum/#!topic/anaconda/zTQQ0NqqIvk>`_. If this happens, a workaround could be to create a new user without whitespaces (e.g. ``C:\John_doe\``).
 
===============================
Installing Pyntacle from source
===============================

Installing from source is advised for advanced users only. The following instructions were written for Mac OS v.11+ and a few major Linux distros. System requirements can vary for other distros/versions.

The source code can be downloaded from our GitHub `releases <https://github.com/mazzalab/pyntacle/releases>`_ page as a .tar.gz file. Before trying to install Pyntacle, there are system requirements that need to be satisfied on each platform.

--------------
Debian, Ubuntu
--------------

As a user with admin rights, run:

::

 apt-get install -y build-essential linux-headers-$(uname -r) libgl1-mesa-glx libigraph0v5 libigraph0-dev libffi-dev libjpeg-dev libgif-dev libblas-dev liblapack-dev git python3-pip python3-tk

.. note:: **Ubuntu/Debian version <= 16.04**:
   For Ubuntu/Debian 16.04 and older, you also have to install two dependencies from the PyPi repository, by running:

   ::

    apt-get install -y build-essential linux-headers-$(uname -r) libgl1-mesa-glx libigraph0v5 libigraph0-dev libffi-dev libjpeg-dev libgif-dev libblas-dev liblapack-dev git python3-pip python3-tk


Finally, extract the Pyntacle `source tar.gz file <https://github.com/mazzalab/pyntacle/releases/latest>`_ navigate into it and run as an administrator (or add ``--user`` if you do not have admin rights and prefer to install the Pyntacle binary in ``~/.local/bin``):

::

  python3 setup.py install


--------
CentOS 7
--------

As an admin, you need to run:

::

  yum groupinstall -y development kernel-headers-`uname -r` kernel-devel-`uname -r` gcc gcc-c++ yum-utils; yum install -y https://centos7.iuscommunity.org/ius-release.rpm; yum install -y wget python36u-devel.x86_64 igraph-devel.x86_64 atlas-devel.x86_64 libffi-devel.x86_64 python36u-pip python36u-tkinter.x86_64

Finally, extract the Pyntacle `source tar.gz file <https://github.com/mazzalab/pyntacle/releases/latest>`_ navigate into it and run as an administrator (or add ``--user`` if you do not have admin rights and prefer to install the Pyntacle binary in ``~/.local/bin``):


::

  python3.6 setup.py install

--------
Mac OS X
--------

In order to compile from source, you need some of the tools that are conveniently packed in `XCode <https://itunes.apple.com/us/app/xcode/id497799835?mt=12>`_, which has to be downloaded and installed from the Mac App Store.
Once you have XCode - and you have opened at least once -, you will need to install the XCode Command Line Tools, by opening a terminal, typing:

::

  xcode-select --install

and following the prompt on screen.

Additionally, you need other dependencies to compile Pyntacle. You can easily fetch them using the package manager `Mac Ports <https://www.macports.org/install.php>`_.

Once Mac Ports is installed, getting the dependencies is easy:

::

  port install py35-setuptools py35-pandas py35-seaborn py35-colorama py35-xlsxwriter py35-igraph

Note: unfortunately, at the time of writing this guide, Mac Ports does not provide a python3.6 version of the library 'xlsxwriter'; therefore, everything must be downgraded to Python 3.5. This does not affect the performance or the results.

Finally, extract the Pyntacle `source tar.gz file <https://github.com/mazzalab/pyntacle/releases/latest>`_ navigate into it and run as an administrator:

::

  python3.5 setup.py install
  ln -s /opt/local/Library/Frameworks/Python.framework/Versions/3.5/bin/Pyntacle /opt/local/bin


============
CUDA support
============

Independently of the OS in use, if you need CUDA support, you should also install the CUDA toolkit by downloading and installing the Toolkit from the `NVIDIA website <https://developer.nvidia.com/cuda-toolkit>`_.
pyntacle.io\_stream.converter module
====================================

.. automodule:: pyntacle.io_stream.converter
    :members:
    :undoc-members:
    :show-inheritance:
pyntacle.tools.adjmatrix\_utils module
======================================

.. automodule:: pyntacle.tools.adjmatrix_utils
    :members:
    :undoc-members:
    :show-inheritance:
pyntacle.tools.octopus module
=============================

.. automodule:: pyntacle.tools.octopus
    :members: Octopus
    :undoc-members:
    :show-inheritance:
pyntacle.io\_stream.import\_attributes module
=============================================

.. automodule:: pyntacle.io_stream.import_attributes
    :members:
    :undoc-members:
    :show-inheritance:
pyntacle.algorithms.global\_topology module
===========================================

.. automodule:: pyntacle.algorithms.global_topology
    :members:
    :undoc-members:
    :show-inheritance:
pyntacle.io\_stream.importer module
===================================

.. automodule:: pyntacle.io_stream.importer
    :members: PyntacleImporter
    :undoc-members:
    :show-inheritance:
.. Pyntacle documentation master file, created by
   sphinx-quickstart on Wed Nov 21 14:57:02 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Pyntacle's documentation!
====================================

.. toctree::
   :maxdepth: 6
   :caption: Contents:

   installation
   pyntacle


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
pyntacle.algorithms.sparseness module
=====================================

.. automodule:: pyntacle.algorithms.sparseness
    :members:
    :undoc-members:
    :show-inheritance:

.. |br| raw:: html

   <br />
