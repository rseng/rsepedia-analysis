<img src="docs/src/images/Magritte_logo_plain.svg" alt="logo" width="350"/>

[![Build Status](https://travis-ci.com/Magritte-code/Magritte.svg?branch=master)](https://travis-ci.com/Magritte-code/Magritte)
[![Documentation Status](https://readthedocs.org/projects/magritte/badge/?version=stable)](https://magritte.readthedocs.io/en/stable/?badge=stable)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03905/status.svg)](https://doi.org/10.21105/joss.03905)
---

Welcome to the Magritte repository! Magritte is a open-source software library for
simulating radiation transport, developed at
[University College London](https://www.ucl.ac.uk/) (UCL, UK) and
[KU Leuven](https://www.kuleuven.be/english/) (Belgium).

Magritte is currently mainly used for post-processing hydrodynamical simulations of
astrophysical models by creating synthetic observations, but the techniques could
in principle also be applied more general (e.g. plasma and nuclear physics).  

It can either be used as a Python package or as a C++ library.
Magritte uses a deterministic ray-tracer with a formal solver that currently focusses on
line radiative transfer (see
[De Ceuster et al. 2019](https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1812D/abstract)
for more details). The plot below shows synthetic observations made with
Magritte at three different inclinations for a hydro simulation of the stellar wind around
an asymptotic giant branch (AGB) star as it is perturbed by a companion.

<img src="docs/src/_static/movie.gif" alt="movie"/>

The hydro model was created by Jan Bolte using [MPI-AMRVAC](http://amrvac.org/). See the
[examples](https://magritte.readthedocs.io/en/latest/1_examples/index.html) in the
documentation to learn how these synthetic observations were created with Magritte.


## Documentation
Please find our online documentation [here](https://magritte.readthedocs.io).
In particular, see our
[quickstart documentation](https://magritte.readthedocs.io/en/latest/0_getting_started/0_quickstart.html)
for a quick intro.


## Papers about Magritte
The following list of papers might provide further insights in the inner workings of
Magritte:
* _Magritte II: Adaptive ray-tracing, mesh construction and reduction_
([arXiv](https://arxiv.org/abs/2011.14998), [MNRAS](https://doi.org/10.1093/mnras/staa3199));
* _Magritte I: Non-LTE atomic and molecular line modelling_
([arXiv](https://arxiv.org/abs/1912.08445),
[MNRAS](https://doi.org/10.1093/mnras/stz3557));
* _3D Line Radiative Transfer & Synthetic Observations with
Magritte_ ([JOSS](https://doi.org/10.21105/joss.03905)).

Please note that some features presented in these papers might not yet be implemented
and documented in the latest release of Magritte.


## Issues & Contact
Please report any issues with Magritte or its documentation
[here](https://github.com/UCL/Magritte/issues).
If you need any further help, please contact
[Frederik De Ceuster](https://freddeceuster.github.io/).


## Cite
Please contact the authors of the papers referenced above if you want to use
Magritte in your research. We are currently working on documentation and
examples to facilitate its independent use. Until then, please
[contact us](https://freddeceuster.github.io/).


## Developers & Contributors
**Developers**
* Frederik De Ceuster
* Thomas Ceulemans
* Atulit Srivastava

**Scientific & Technical advisors**
* Ward Homan
* Jan Bolte
* Jeremy Yates
* Leen Decin
* Peter Boyle
* James Hetherington

**Contributors**
* Silke Maes
* Jolien Malfait


## Acknowledgements
FDC was for most of the development supported by the EPSRC iCASE studentship programme, Intel Corporation and Cray Inc.
FDC, JB, WH, and LD acknowledge support from the ERC consolidator grant 646758 AEROSOL.
TC is a PhD fellow of the Research Foundation - Flanders (FWO).# Contributing to Magritte

Thank you for considering to contribute to Magritte! We highly appreciate any efforts to improve our software.
In case you have any suggestions, questions, possible bugs, or feature requests, please, open an issue [here](https://github.com/Magritte-code/Magritte/issues).


## Contributing examples and use cases

If you have used Magritte in a way that is different from any of [our examples](https://magritte.readthedocs.io/en/stable/1_examples/index.html), we highly encourage you to contribute to the [examples](https://magritte.readthedocs.io/en/stable/1_examples/index.html) in the documentation.
You can do this by preparing your example as a jupyter notebook and creating a pull request (see below). If you are unsure about your example or how to add it, please contact us and we can help you.


## Contributing new functionality

We have a clear development strategy for Magritte, which is closely related to our reseach. Several features are being developed for several projects, all of which might not yet be publicly visible. Therefore, to prevent development overlap, we advise to contact us (e.g. by creating an [issue](https://github.com/Magritte-code/Magritte/issues)), before starting to work on Magritte yourself. We are happy to share our preliminary features if they could help you, or welcome you to our develoment team to colaboratively work on our projects.


## Contributing bug fixes

For bug fixes or smaller contributions, we gladly accept pull requests (see below).


## Creating pull requests

The advised procedure for creating pull requests for Magritte goes as follows:

- Fork the Magritte repository and create a new branch for your work;
- Add your changes and commit them;
- Run the benchmarks and/or the relevant examples;
- Open a pull request for your branch on Github.

Also, have a look at [this guide](https://docs.github.com/en/get-started/quickstart/contributing-to-projects) on forking in Github.


_Thank you for contributing to Magritte!_# magritte

This directory contains the Magritte python package.# Dependencies

Some installers for dependencies.
# Data

This folder contains line data files form the (LAMDA database)[https://home.strw.leidenuniv.nl/~moldata/].
These are only for testing.
Please always download the most recent data from the LAMDA website for your research and properly cite its origin.
# Models
A directory to store test models.
# Magritte Documentation

Magritte's documentation is build using Doxygen, Breathe and Sphinx.
The inline documentation in the source code is parsed by Doxygen and turned into
xml. The xml is then used by sphinx to create html documentation which is hosted
on readthedocs.
Magritte documentation
######################

Welcome to the Magritte documentation! Magritte is an open-source software
library for simulating radiation transport, developed at `University College
London <https://www.ucl.ac.uk/>`_ (UCL, UK) and `KU Leuven
<https://www.kuleuven.be/english/>`_ (Belgium).

Magritte is currently mainly used for post-processing hydrodynamical simulations by
creating synthetic observations, but the techniques could also be applied more general.
It can either be used as a Python package or as a C++ library.
Magritte uses a deterministic ray-tracer with a formal solver that currently focusses on
line radiative transfer (see
`De Ceuster et al. 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1812D/abstract>`_
for more details). The plot below shows synthetic observations made with
Magritte at three different inclinations for a hydro simulation of the stellar wind around
an asymptotic giant branch (AGB) star as it is perturbed by a companion.

.. raw:: html

    <video width="100%" controls playsinline autoplay muted loop>
      <source src="_static/movie.webm">
      Your browser does not support the video tag.
    </video>

The hydro model was created by Jan Bolte using `MPI-AMRVAC <http://amrvac.org/>`_. See
the :ref:`examples <link-examples>` to learn how these synthetic observations were created
with Magritte.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   0_getting_started/index
   1_examples/index
   2_benchmarks/index
   3_python_api_documentation/index
..   4_cpp_api_documentation/index


Papers about Magritte
*********************

The following list of papers might provide further insights in the inner workings of
Magritte:

* Magritte: **Adaptive ray-tracing, mesh construction and reduction**,
  *F. De Ceuster, J. Bolte, W. Homan, S. Maes, J. Malfait, L. Decin, J. Yates, P. Boyle, J. Hetherington*, 2020
  (`arXiv <https://arxiv.org/abs/2011.14998>`_,
  `MNRAS <https://doi.org/10.1093/mnras/staa3199>`_);


* Magritte: **Non-LTE atomic and molecular line modelling**,
  *F. De Ceuster, W. Homan, J. Yates, L. Decin, P. Boyle, J. Hetherington*, 2019
  (`arXiv <https://arxiv.org/abs/1912.08445>`_,
  `MNRAS <https://doi.org/10.1093/mnras/stz3557>`_);
  
* **3D Line Radiative Transfer & Synthetic Observations with Magritte**,
  *F. De Ceuster, T. Ceulemans, A. Srivastava, W. Homan, J. Bolte, J. Yates, L. Decin, P. Boyle, J., Hetherington*
  (`JOSS <https://doi.org/10.21105/joss.03905>`_).

Please note that some features presented in these papers might not yet be implemented
and documented in the latest release of Magritte.


Issues & Contact
****************

Please report any `issues <https://github.com/Magritte-code/Magritte/issues>`_ with
Magritte or its documentation `here <https://github.com/Magritte-code/Magritte/issues>`_.
If you need any further help, please contact `Frederik De Ceuster
<https://www.kuleuven.be/wieiswie/en/person/00101884>`_.


Developers & Contributors
*************************

**Developers**

* Frederik De Ceuster
* Thomas Ceulemans
* Atulit Srivastava

**Scientific & Technical advisors**

* Ward Homan
* Jan Bolte
* Jeremy Yates
* Leen Decin
* Peter Boyle
* James Hetherington

**Contributors**

* Silke Maes
* Jolien Malfait


Acknowledgements
****************

FDC is supported by the EPSRC iCASE studentship programme, Intel Corporation and Cray Inc.
FDC, JB, WH, and LD acknowledge support from the ERC consolidator grant 646758 AEROSOL.
TC is a PhD fellow of the Research Foundation - Flanders (FWO)... _link-examples:


Examples
########

The following examples demonstrate how to use Magritte as a Python library.

.. Please refer to the :ref:`C++ API documentation <link-cpp_api_documentation>`
.. to see how Magritte can be used as a C++ library.


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   0_creating_models/index
   1_post-processing/index
   2_synthetic_observations/index
.. _link-synthetic_observations:


Synthetic observations
######################

The following examples show how to create synthetic observations with Magritte.

.. Note::
    In these examples we assume to have a Magritte model available. For more information
    on how to create a Magritte model, please refer to the :ref:`creating models
    <link-creating_models>` or :ref:`post-processing examples <link-post-processing>`.

.. toctree::
   :maxdepth: 1
   :caption: Contents
   
   0_image_analytic_disk.ipynb
   1_image_analytic_spiral.ipynb
   2_image_1D_AMRVAC.ipynb
   3_image_3D_AMRVAC.ipynb
   4_image_3D_AMRVAC_red.ipynb
   5_image_3D_Phantom.ipynb
   6_image_3D_Phantom_red.ipynb.. _link-creating_models:


Creating models
###############

The following examples show how to create Magritte models from analytic or tabulated data.

.. Note::
    For specific examples on how to generate Magritte models from snapshots of AMRVAC
    or Phantom hydrodynamics simulations, please refer to the :ref:`post-processing
    examples <link-post-processing>`.


.. toctree::
   :maxdepth: 1
   :caption: Contents:
   
   0_create_analytic_disk.ipynb
   1_create_analytic_spiral.ipynb

.. _link-post-processing:


Post-processing hydro snapshots
###############################

We show how Magritte can be used to post-process snapshots of hydrodynamics simulations.
We consider snapshots of both `MPI-AMRVAC <http://amrvac.org/>`_, a solver based on
adaptive mesh refinement (AMR, `Xia et al. 2018
<https://ui.adsabs.harvard.edu/abs/2018ApJS..234...30X/abstract>`_), and `Phantom
<https://phantomsph.bitbucket.io/>`_ which uses smoothed-particle hydrodynamics (SPH,
`Price et al. 2018 <https://ui.adsabs.harvard.edu/abs/2018PASA...35...31P/abstract>`_).

.. Note::
    For general examples on how to generate Magritte models from analytic or tabulated
    data, please refer to the :ref:`creating model examples <link-creating_models>`.

First we show how to build Magritte models for these snapshots and then we demonstrate
how to reduce them to smaller approximate models for fast radiative transfer
(`De Ceuster et al. 2020
<https://ui.adsabs.harvard.edu/abs/2020MNRAS.499.5194D/abstract>`_). 


.. toctree::
   :maxdepth: 1
   :caption: AMRVAC snapshots:
   
   0_create_AMRVAC_1D.ipynb
   1_create_AMRVAC_3D.ipynb
   2_reduce_AMRVAC_3D.ipynb

.. toctree::
   :maxdepth: 1
   :caption: Phantom snapshots:  

   3_create_Phantom_3D.ipynb
   4_reduce_Phantom_3D.ipynb

.. _link-cpp_api_documentation:

C++ API documentation
#####################

All source code for Magritte is available on `GitHub
<https://github.com/Magritte-code/Magritte>`_.

.. Warning::
    These docs are far from complete.

.. doxygenstruct:: Geometry
   :members:
magritte.setup module
=====================

.. automodule:: magritte.setup
   :members:
   :undoc-members:
   :show-inheritance:
magritte.mesher module
======================

.. automodule:: magritte.mesher
   :members:
   :exclude-members: get_pfs, get_nbs, norm, get_Gs, get_Ls, get_weights
   :special-members: __init__
   :show-inheritance:

Python API documentation
########################

.. toctree::
   :maxdepth: 3

   magrittemagritte package
================

.. automodule:: magritte
   :members:
   :undoc-members:
   :show-inheritance:

Submodules
----------

.. toctree::
   :maxdepth: 0

   magritte.core
   magritte.mesher
   magritte.plot
   magritte.setup
   magritte.tools
magritte.plot module
====================

.. automodule:: magritte.plot
   :members:
   :exclude-members: matplotlib_to_plotly
   :show-inheritance:
magritte.core module
====================

.. automodule:: magritte.core
   :members:
   :exclude-members: BoundaryCondition, Solver, Lambda, MReal, TReal, VReal, VSize, VVector3D, vCollisionPartner, vLineProducingSpecies
   :show-inheritance:
..   :undoc-members:

magritte.tools module
=====================

.. automodule:: magritte.tools
   :members:
   :undoc-members:
   :show-inheritance:
.. _link-quickstart:

Quickstart
##########

.. Warning::
    Please ensure that the :ref:`prerequisites <link-prerequisites>` are in place before continuing with the installation.

.. Warning::
    This is the (too) quick intallation guide. For details, see our
    :ref:`comprehensive installation <link-installation>` guide.


Clone
*****

Magritte has to be compiled from its source code, which can be cloned using:

.. code-block:: shell

    git clone --recursive https://github.com/Magritte-code/Magritte.git

from our `GitHub <https://github.com/Magritte-code/Magritte>`_ repository.


Setup
*****

The Magritte `conda <https://www.anaconda.com/products/individual>`_ environment
can be created from an environment file with:

.. code-block:: shell

    conda env create -f dependencies/conda_env.yml

which installs all required packages in the :literal:`magritte` conda
environment, and can be activated with:

.. code-block:: shell

    conda activate magritte

This environment has to be active whenever Magritte is compiled or used!

.. Note::
    If you didn't configure your shell with the Anaconda installation, the above command won't work, but the same result can be achieved with:

    .. code-block:: shell

        source activate magritte


Build
*****

Magritte can easily be build in its default configuration with:

.. code-block:: shell

    bash build.sh

This will create the necessary files for the Magritte python package.


Run
***

If all the above worked, go download and experiment with some of our :ref:`examples
<link-examples>`!
.. _link-prerequisites:

Prerequisites
#############

Throughout this documentation, we assume a Unix-based operating system (e.g. Linux or MacOS).
For Windows users, we recommend to use the Windows Subsystem Linux (WSL; see e.g. `here <https://www.windowscentral.com/install-windows-subsystem-linux-windows-10>`_).

Below is a list of the software required to be able to compile and run Magritte.

* `GCC <https://gcc.gnu.org/>`_ (version :literal:`5.0.0` or later): to compile the C++ part of Magritte.
* `CMake <https://cmake.org/>`_ (version :literal:`3.18.0` or later): for building Magritte, organising compilation and linking;
* `Anaconda <https://www.anaconda.com/blog/individual-edition-2020-11>`_ (optional): for managing the required Python packages;
* `Gmsh <https://gmsh.info/>`_ (optional): for custructing geometric meshes for the models;

Although Anaconda is optional and you could, in priciple, use any other way of installing the required Python dependencies, we recommend and assume the use of Anaconda.
The list of required Python packages can be found in the `conda environment file <https://github.com/Magritte-code/Magritte/blob/stable/dependencies/conda_env.yml>`_.
Also Gmsh is optional, however, it will be used in most examples in this documentation to create models.

Once these prerequisites are in place, Magritte can be installed following the :ref:`quickstart <link-quickstart>`.

.. Warning::
    On MacOS, gcc by default points to Apple Clang (which is NOT the GNU compiler). Apple Clang won't work, due to an issue with OpenMP.
    Hence, GNU gcc has to be installed separately, e.g. using `Homebrew <https://brew.sh/>`_ or `Macports <https://www.macports.org/>`_.


.. Hint::
    Although we recommend to first try to install the required software in your preferred usual way, the folder `Magritte/dependencies/ <https://github.com/Magritte-code/Magritte/tree/stable/dependencies>`_ contains a set of bash scripts: :literal:`get_conda.sh`, :literal:`get_cmake.sh`,
    and :literal:`get_gmsh.sh`, to quickly obtain the appropriate versions of miniconda3 (for Anaconda), CMake, and
    Gmsh respectively. However, these scripts don't properly install the software, but rather download an executable binary, which is useful, e.g. if you don't have root access.
    
    To use these, first clone the Magritte `GitHub <https://github.com/Magritte-code/Magritte>`_ repository
    
    .. code-block:: shell

        git clone --recursive https://github.com/Magritte-code/Magritte.git
    
    enter the newly created Magritte directory, and from there run the relevant script.
    
    .. code-block:: shell
    
        bash dependencies/get_conda.sh
        
    .. code-block:: shell
    
        bash dependencies/get_cmake.sh
        
    .. code-block:: shell
    
        bash dependencies/get_gmsh.sh

    ⚠️  Don't forget to include the relevant paths to your shell's :literal:`$PATH` variable!
    (For convenience, the scripts will print the relevant path.)

.. Note::
    We are working on a set of pre-compiled versions of Magritte that would greatly simply the installation process.
    `Let us know <https://github.com/Magritte-code/Magritte/issues>`_, if you have any suggestions or want to help!
.. _link-installation:

Installation
############

.. note::

    This is the comprehensive installation guide. For a quick intro, see our
    :ref:`quickstart <link-quickstart>` guide.


Download
********

Magritte has to be compiled from its source code, which can be cloned using:

.. code-block:: shell

    git clone --recursive https://github.com/Magritte-code/Magritte.git

from our `GitHub <https://github.com/Magritte-code/Magritte>`_ repository. This
creates the directory, :literal:`Magritte`, which will be refered to as the
Magritte root directory. Ensure to include the :literal:`--recursive` to also
clone the required submodules.


Dependencies
************

Magritte has several dependencies, some of which are optional.


**Submodules**

* `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_, version :literal:`3.3.7` or later, for some of the linear algebra;

* `pybind11 <https://github.com/pybind/pybind11>`_, version :literal:`2.2.4` or later, for binding C++ to python;

* `googletest <https://github.com/google/googletest>`_ version :literal:`1.10.0` or later, for some of the tests;

* `Paracabs <https://github.com/Magritte-code/Paracabs>`_, our custom parallelization and acceleration abstractions.

**Required**

* `GCC <https://gcc.gnu.org/>`_, version :literal:`5.0.0` or later, to compile the C++ part of Magritte;
* `CMake <https://cmake.org/>`_, version :literal:`3.18.0` or later, for building the library, organising compilation and linking;


**Optional**

* `Anaconda <https://www.anaconda.com/blog/individual-edition-2020-11>`_, for managing the required Python packages;
* `Gmsh <https://gmsh.info/>`_, version :literal:`4.6.0` or later, for meshing model geometries.


Please note that :literal:`Paracabs` might have further dependencies depending
on which paralellization and acceleration libraries are used. See
:ref:`advanced compilation <link-advanced_compilation>` for further details.


Python packages & Environment
*****************************

To make Magritte useful, some python packages have to be installed which
is best done using a `conda <https://www.anaconda.com/products/individual>`_ environment.

The following python packages are used in Magritte, mainly for io and to create
the model files.

* :mod:`numpy`, to bind the Magritte data;
* :mod:`h5py`, to read and write HDF5 data files;
* :mod:`scipy`, for interpolation and spatial functions such as nearest neighbour calculations;
* :mod:`healpy`, to sample directions from a discretized unit sphere;
* :mod:`astropy`, for unit conversions and physical constants;
* :mod:`meshio`, for reading and writing several types of mesh data structures;
* :mod:`vtk`, for reading and writing vtk files;
* :mod:`pyyaml`, for reading and writing yaml files;
* :mod:`mpi4py`, for MPI (Message Passing Interface) functionality in Python;
* :mod:`tqdm`, for progress bars;
* :mod:`numba`, for just-in-time compilation of some Python functions;
* :mod:`palettable`, for nice colourmaps;
* :mod:`matplotlib`, for basic plotting;
* :mod:`plotly`, for advanced plotting;
* :mod:`nodejs`, for interactivity in some advanced plots;
* :mod:`ipywidgets`, for interactive plotting;
* :mod:`jupyterlab`, for convenient use of the jupyter notebooks;

All of these packages can also be found in the `conda environment file <https://github.com/Magritte-code/Magritte/blob/stable/dependencies/conda_env.yml>`_.

.. hint::

    The simplest way to setup the required python packages is using the
    `anaconda <https://www.anaconda.com/products/individual>`_ package manager.
    The Magritte conda environment can be created from the environment
    file :literal:`conda_env.yml` located in the :literal:`dependencies` directory, with

    .. code-block:: shell

        conda env create -f conda_env.yml

    This will download and install all required python packages in a newly created
    :literal:`magritte` conda environment. The environment can be activated with

    .. code-block:: shell

        conda activate magritte

    Please ensure that this environment is active whenever Magritte is compiled or used.

.. warning::

    Magritte uses plotly for some interactive plots. Plotly requires additional
    extensions to be able to render plots in a jupyter notebook or in jupyter lab. Please
    consult their `intstallation notes <https://plotly.com/python/getting-started/>`_ to get
    plotly working with jupyter.


Compilation
***********

Once all dependencies are in place, Magritte can be compiled.

.. hint::

    There is a shortcut script to build Magritte in the default configuration.
    From within the Magritte root directory, run:

    .. code-block:: shell

        bash build.sh

    This will create a :literal:`bin` directory in the Magritte root directory
    containing the library binary files and the executables for the tests. It will
    also create a shared object file :literal:`core.so` in the magritte python package,
    located in the :literal:`magritte` directory.

See :ref:`advanced compilation <link-advanced_compilation>` for further options.



.. _link-advanced_compilation:



Advanced compilation
********************

Compilers
=========

Corrently only the GNU gcc compiler is fully supported.
We are currently further investigating Clang and Intel compiler (:literal:`icc`) support.


Vectorisation
=============

🧑‍💻 Comming soon!


GPU acceleration
================

🧑‍💻 Comming soon!


Getting started
###############

First check if all :ref:`prerequisites <link-prerequisites>` are in place before continuing to install Magritte in the :ref:`quickstart <link-prerequisites>`.
More background on the installation process of Magritte can be found in the :ref:`installation <link-installation>` notes.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   0_prerequisites
   1_quickstart
   2_installation
Numeric benchmarks
##################


Van Zadelhoff et al. (2002)
***************************

The `van Zadelhoff et al. (2002) <https://ui.adsabs.harvard.edu/abs/2002A%26A...395..373V>`_
benchmark comprises two main benchmark problems for non-LTE line radiative transfer.
The results for Magritte on this benchmark and others are discussed and compared with
`RATRAN <https://personal.sron.nl/~vdtak/ratran/frames.html>`_
and `LIME <https://github.com/lime-rt/lime>`_ in
`De Ceuster et al. (2019) <https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1812D>`_.
The scripts to setup and run these benchmarks are provided with Magritte and can be found `here <https://github.com/Magritte-code/Magritte/tree/stable/tests/benchmarks/numeric>`_.


Problem 1
=========

The first problem considers a fictitious two-level species in a spherically symmetric
cloud, without velocity field, with a constant temperature distribution, and a
quadratucally decaying density distribution. A line data file for the fictitious
two-level species (in `LAMDA <https://home.strw.leidenuniv.nl/~moldata/>`_ format)
is provided with Magritte and can be found
`here <https://github.com/Magritte-code/Magritte/blob/master/tests/data/test.txt>`_.


Problem 2
=========

The second problem has a more realistic setup and considers the lines of HCO+ in a snapshot of an inside-out collapse model.
The model consists of 50 logarithmically spaced grid points. In each grid point the radial velocity, gas temperature, micro-turbulence, and both HCO+ and H2 abundances are given.
This data was provided on the benchmark  website (which is offline now), but can still be found `here <https://github.com/Magritte-code/Magritte/blob/stable/tests/benchmarks/numeric/vanZadelhoff_2a.in>`_ and `here <https://github.com/Magritte-code/Magritte/blob/stable/tests/benchmarks/numeric/vanZadelhoff_2b.in>`_.
The line data file for HCO+ comes from the LAMDA data base and can be found `here <https://github.com/Magritte-code/Magritte/blob/stable/tests/data/hco%2B.txt>`_Benchmarks
##########

To demonstrate the validity of Magritte we present some banchmark results.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   0_analytic_benchmarks
   1_numeric_benchmarks
   2_performance_scaling
Analytic benchmarks
###################

Magritte was tested against a set of analytic benchmarks, i.e. models for which an semi-analytic solution is known.
These models and the performance of Magritte on them are described in `De Ceuster et al. (2019) <https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.1812D>`_.
The scripts to setup and run these benchmarks are provided with Magritte and can be found `here <https://github.com/Magritte-code/Magritte/tree/stable/tests/benchmarks/analytic>`_.Performance scaling
###################

The computation time in Magritte is dominated by the computation of the radiation field.
This is computed by solving the radiative transfer equation for each ray through each point
and in each frequency bin. Schematically this computation can be written as:

.. code-block:: C++

    for(direction : directions)
    {
        for(point : points)
        {
            Ray_pair ray_pair = trace (direction, point);

            for(frequency : frequencies)
            {
                solve_radiative_transfer (ray_pair, frequency);
            }
        }
    }


The parallelisation in Magritte is an hybrid between MPI and OpenMP.
The outer loop over the directions is parallelised with MPI,
while the middle loop over the points is parallelised with OpenMP.
The inner loop over the frequencies is ideally suited for vectorisation.
The MPI and OpenMP parall for loops are implemented using the
:literal:`distributed_for` and :literal:`threaded_for` loop respectively
in `Paracabs <https://github.com/Magritte-code/Paracabs>`_.

Multi-threading (OpenMP)
************************

.. image:: strong_scaling_omp.png
   :width: 600

The figure above shows a strong scaling plot for the multi-threading using OpenMP in Magritte.


Message passing (MPI)
*********************

.. image:: strong_scaling_mpi.png
   :width: 600
   
The figure above shows a strong scaling plot for the message passing using MPI in Magritte.
