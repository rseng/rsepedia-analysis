
Quasi-Geostrophic Spectral model (qgs)
======================================

[![DOI](https://zenodo.org/badge/246609584.svg)](https://zenodo.org/badge/latestdoi/246609584)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02597/status.svg)](https://doi.org/10.21105/joss.02597)
[![Documentation Status](https://readthedocs.org/projects/qgs/badge/?version=latest)](https://qgs.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

General Information
-------------------

qgs is a Python implementation of an atmospheric model for midlatitudes.  It models the dynamics of
a 2-layer [quasi-geostrophic](https://en.wikipedia.org/wiki/Quasi-geostrophic_equations) channel
atmosphere on a [beta-plane](https://en.wikipedia.org/wiki/Beta_plane), coupled to a simple land or
[shallow-water](https://en.wikipedia.org/wiki/Shallow_water_equations) ocean component. 

![](./misc/figs/readme.gif)

About
-----

(c) 2020-2021 qgs Developers and Contributors

Part of the code originates from the Python [MAOOAM](https://github.com/Climdyn/MAOOAM) implementation by Maxime Tondeur and Jonathan Demaeyer.

See [LICENSE.txt](./LICENSE.txt) for license information.

**Please cite the code description article if you use (a part of) this software for a publication:**

* Demaeyer J., De Cruz, L. and Vannitsem, S. , (2020). qgs: A flexible Python framework of reduced-order multiscale climate models. 
  *Journal of Open Source Software*, **5**(56), 2597,   [https://doi.org/10.21105/joss.02597](https://doi.org/10.21105/joss.02597).

Please consult the qgs [code repository](http://www.github.com/Climdyn/qgs) for updates.


Installation
------------

The easiest way to run qgs is to use an appropriate environment created through [Anaconda](https://www.anaconda.com/).

First install Anaconda and clone the repository:

    git clone https://github.com/Climdyn/qgs.git

Then install and activate the Python3 Anaconda environment:

    conda env create -f environment.yml
    conda activate qgs

You can then perform a test by running the script

    python qgs_rp.py
    
to see if everything runs smoothly (this should take less than a minute).

#### Note for Windows and MacOS users

Presently, qgs is compatible with Windows and MacOS but users wanting to use qgs inside their Python scripts must guard the main script with a

```
if __name__ == "__main__":
```

clause and add the following lines below

```
  from multiprocessing import freeze_support
  freeze_support()
```

About this usage, see for example the main scripts `qgs_rp.py` and `qgs_maooam.py` in the root folder.
Note that the Jupyter notebooks are not concerned by this recommendation and work perfectly well on both operating systems.

> **Why?** These lines are required to make the multiprocessing library works with these operating systems. See [here](https://docs.python.org/3.8/library/multiprocessing.html) for more details, 
> and in particular [this section](https://docs.python.org/3.8/library/multiprocessing.html#the-spawn-and-forkserver-start-methods).


#### Activating DifferentialEquations.jl optional support

In addition to the qgs builtin Runge-Kutta integrator, the qgs model can alternatively be integrated with a package called [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) written in [Julia](https://julialang.org/), and available through the
[diffeqpy](https://github.com/SciML/diffeqpy) python package.
The diffeqpy package first installation step is done by Anaconda in the qgs environment but then you must [install Julia](https://julialang.org/downloads/) and follow the final manual installation instruction found in the [diffeqpy README](https://github.com/SciML/diffeqpy).

These can be summed up as opening a terminal and doing:
```
conda activate qgs
python
```
and then inside the Python command line interface do:

```
>>> import diffeqpy
>>> diffeqpy.install()
```
which will then finalize the installation. An example of a notebook using this package is available in the documentation and on [readthedocs](https://qgs.readthedocs.io/en/latest/files/examples/diffeq.html).

Documentation
-------------

To build the documentation, please run (with the conda environment activated):

    cd documentation
    make html

You may need to install [make](https://www.gnu.org/software/make/) if it is not already present on your system.
Once built, the documentation is available [here](./documentation/build/html/index.html).

The documentation is also available online on read the docs: [https://qgs.readthedocs.io/](https://qgs.readthedocs.io/)

Usage
-----

qgs can be used by editing and running the script `qgs_rp.py` and `qgs_maooam.py` found in the main folder.

For more advanced usages, please read the [User Guides](https://qgs.readthedocs.io/en/latest/files/user_guide.html).

Examples
--------

Another nice way to run the model is through the use of Jupyter notebooks. 
Simple examples can be found in the [notebooks folder](./notebooks).
For instance, running 

    conda activate qgs
    cd notebooks
    jupyter-notebook
    
will lead you to your favorite browser where you can load and run the examples.

Dependencies
------------

qgs needs mainly:

   * [Numpy](https://numpy.org/) for numeric support
   * [sparse](https://sparse.pydata.org/) for sparse multidimensional arrays support
   * [Numba](https://numba.pydata.org/) for code acceleration
   * [Sympy](https://www.sympy.org/) for symbolic manipulation of inner products
   
Check the yaml file [environment.yml](./environment.yml) for the dependencies.

Forthcoming developments
------------------------

* Scientific development (short-to-mid-term developments)
    + Dynamical equilibrium temperature equations
    + Non-autonomous equation (seasonality, etc...)
    + Quartic T⁴ temperature tendencies
* Technical mid-term developments
    + Dimensionally robust Parameter class operation
* Long-term development track
    + Active advection
    + True quasi-geostrophic ocean when using ocean model version
    + Salinity in the ocean 
    + Symbolic PDE equation specification
    + Numerical basis of function
  
Contributing to qgs
-------------------

If you want to contribute actively to the roadmap detailed above, please contact the main authors.

In addition, if you have made changes that you think will be useful to others, please feel free to suggest these as a pull request on the [qgs Github repository](https://github.com/Climdyn/qgs).

More information and guidance about how to do a pull request for qgs can be found in the documentation [here](https://qgs.readthedocs.io/en/latest/files/general_information.html#contributing-to-qgs).

Other atmospheric models in Python
----------------------------------

Non-exhaustive list:

* [Q-GCM](http://q-gcm.org/): A mid-latitude grid based ocean-atmosphere model like MAOOAM. Code in Fortran,
                                interface is in Python.
* [pyqg](https://github.com/pyqg/pyqg): A pseudo-spectral python solver for quasi-geostrophic systems.
* [Isca](https://execlim.github.io/IscaWebsite/index.html): Research GCM written in Fortran and largely
            configured with Python scripts, with internal coding changes required for non-standard cases.
            
            
---
title: 'qgs: A flexible Python framework of reduced-order multiscale climate models'
tags:
  - Python
  - Numba
  - Idealized atmospheric model
  - Coupled model
  - Mid-latitude climate variability
authors:
  - name: Jonathan Demaeyer
    orcid: 0000-0002-5098-404X 
    affiliation: 1
  - name: Lesley De Cruz
    orcid: 0000-0003-4458-8953
    affiliation: 1
  - name: Stéphane Vannitsem 
    orcid: 0000-0002-1734-1042
    affiliation: 1
affiliations:
 - name: Royal Meteorological Institute of Belgium, Avenue Circulaire, 3, 1180 Brussels, Belgium
   index: 1

date: 7 August 2020
bibliography: joss.bib

---

# Summary

`qgs` is a a Python implementation of a set of idealized reduced-order models representing atmospheric mid-latitude variability. 
It consists of a two-layer *q*uasi-*g*eostrophic *s*pectral (`qgs`) model of the atmosphere on a beta-plane, coupled either to a simple land surface or to a shallow-water ocean. 
The model's dynamical fields include the atmospheric and oceanic streamfunction and temperature fields, and the land temperature field.

* In the case where it is coupled to an ocean, it reproduces the Modular Arbitrary-Order Ocean-Atmosphere Model (MAOOAM), described in @DDV2016. In @VDDG2015, a 36-variable configuration of this model was shown to reproduce a 
low-frequency variability (LFV) typical of the coupled ocean-atmosphere system. 
This coupling consists in both mechanical and heat exchange interactions between the two components. 
The model has already been used in different contexts, in particular for data assimilation [@PBBCDSY2019; @TCVB2020], 
and predictability studies [@VSD2019; @VD2020]
* In the case of a land surface coupling, it emulates the model proposed in @RP1982 and @CT1987 with a simple thermal relaxation toward a climatological temperature and a mechanical coupling due to the 
friction between the land and the atmosphere. It can also emulate the model proposed in @LHHBD2018, with mechanical coupling and heat exchange. In addition, the number of dynamical spectral modes can be configured by the user, as is the case for the MAOOAM model.

In the `qgs` framework, the partial differential equations (PDEs) that govern the time evolution of its fields are projected on a basis of functions defined on its 
spatial domain. 
This kind of decomposition transforms the PDEs into a set of ordinary differential equations (ODEs) which can then be solved with the usual integration techniques.
Presently in `qgs`, the functions of the basis are chosen amongst the orthogonal Fourier modes compatible with the boundary conditions of each subcomponent of the system, namely the atmosphere, and the ocean or the land surface coupled to it. 
A future development is planned that will enable the user to specify the basis of functions for each component, depending on the required boundary conditions.

The model implementation consists of submodules to set up the model's parameters and to compute the tensor that defines the coefficients of the system of ODEs^[More details about the implementation can be found in @DDV2016 and in the *Code Description* section of the included documentation.].
This tensor is used by the code to compute the tendencies function and its Jacobian matrix. These functions can then be fed to the `qgs` built-in Runge-Kutta integrator or 
to another integrator implemented by the user. As an example, the usage of the Julia `DifferentialEquations.jl` [@RN2017] integration package through the Python `diffeqpy` [@diffeqpy] package is provided.
The tangent linear and adjoint models [@K2003] are also available and allow one to easily conduct data assimilation and linear sensitivity analysis experiments.

The model implementation uses NumPy [@vCV2011; @O2006] and SciPy [@scipy] for arrays and computations support, as well as Numba [@numba] and sparse [@sparse] to considerably accelerate the tensor products computation used to compute the tendencies.


# Statement of need

In atmospheric and climate sciences, research and development is often first conducted with a simple idealized system like the Lorenz-$N$ models ($N \in \{63, 84, 96\}$) [@L63; @L84; @L96] which are toy models of atmospheric variability. 
The first two models are heavily truncated systems (3-variable) describing the very large synoptic-scale dynamics of the single-component atmosphere, that neglect the interaction with other components of the climate system and with smaller scales.
The third one is based on reasonable heuristic assumptions on the spatial dynamics along a latitude, which may however lead to unrealistic statistical features. 

Reduced-order spectral quasi-geostrophic models of the atmosphere with a large number of modes offer better representations of the dry atmospheric dynamics [@OB1989]. 
The dynamics thus obtained allow one to identify typical features of the atmospheric circulation, such as blocked and zonal circulation regimes, and low-frequency variability.
However, these models are less often considered in literature, despite their demonstration of more realistic behavior.

`qgs` aims to popularize these systems by providing a fast and easy-to-use Python framework for researchers and teachers to integrate this kind of model. 
For an efficient handling of the model by users, its documentation is conceived such that its equations and parameters are explained and linked to the code.
In the future, its development will be done in a modular fashion which enables the connection of the atmosphere to various other subsystems and the use of built-in and external toolboxes.

The choice to use Python was specifically made to facilitate its use in Jupyter Notebooks [@jupyter] and with the multiple recent machine learning libraries that are available in this
language.

# State of the field

Other software might interest the reader in need of an easy-to-use idealized atmospheric model.

* MAOOAM: The Modular Arbitrary-Order Ocean-Atmosphere Model, a coupled ocean-atmosphere model included in `qgs` [@MAOOAM]. 
          Code available in Lua, Fortran and Python.
* q-gcm: A mid-latitude grid-based quasi-geostrophic ocean-atmosphere model with two oceanic layers. Code in Fortran, interface in Python [@qgcm].
* pyqg: A pseudo-spectral Python solver for quasi-geostrophic systems [@pyqg]. Allow one to create and solve multiple-layers quasi-geostrophic systems.
* Isca: A research General Circulation Model (GCM) to simulate global dynamics. Written in Fortran and largely
        configurable with Python scripts, with internal coding changes required for non-standard cases [@Vetal2018; @Isca].

`qgs` distinguishes itself from these other models by the combination of a simplified and configurable geometry, a spectral discretization, an efficient numerical implementation of the ODE system by a sparse tensor multiplication, 
and the availability of the tangent linear and adjoint models. As such it is very suitable to quickly simulate very long time periods while capturing key aspects of the climate dynamics at mid-latitudes.

The mechanically coupled atmosphere-land configuration of `qgs` was used to test new ideas using response theory to adapt statistical postprocessing schemes to model changes [@DV2020].
The MAOOAM model configuration of `qgs` was recently considered to perform strongly-coupled data assimilation experiments in the ocean-atmosphere system [@CBDGRV2020].

# Performance

The performance of the `qgs` MAOOAM implementation has been benchmarked against the Lua and Fortran implementations of this model [@MAOOAM].
This comparison was done on a recent Intel CPU with 12 cores, with two different model resolutions: one used in @VDDG2015 and one truncated at the wavenumber 6 for both the oceanic and atmospheric components.
The former leads to a 36-dimensional system of ODEs while the latter is higher-dimensional, using 228 variables.
 
In both cases, all the different code implementations have been initialized with the same initial data and parameters, except for the length of the trajectory being computed. 
The low-dimensional system was integrated for 10$^7$ timeunits (roughly $\sim$ 1850 years) while the higher-dimensional one was integrated over 10$^6$ timeunits ($\sim$ 185 years).
In the case of the Fortran implementation, two different compilers (GNU Gfortran and Intel Ifort) with two different levels of optimization (O2 and O3) have been tested, but no significant differences between these compilers and options were found.
In addition, two different built-in integration modules of `qgs` have been considered: a non-parallel integrator located in the module `integrators.integrate` and a parallel one located in the module `integrators.integrator`.
The latter can integrate multiple trajectories simultaneously, but for the purpose of the benchmark, only one trajectory was computed with it, the other implementations being non-parallel.

The results of this benchmark are depicted on \autoref{fig:benchmark} and show that `qgs`, while not the fastest implementation of MAOOAM available, is a fair competitor. 
The time difference is in general not greater than a factor 5 and tends to be less for high-dimensional model configurations, with an integration time 
roughly the same as the Lua implementation. We note that there is also a significant difference between the parallel and non-parallel implementation of `qgs`, but this difference also seems to vanish for higher-resolution model configurations.
In any case, the parallel integrator of `qgs` can straightforwardly integrate multiple trajectories simultaneously and therefore has an advantage over the non-parallel one (provided that multiple CPU cores are available).
A final remark is that the initial Python version of MAOOAM [found in @MAOOAM] takes 283 minutes to integrate the low-resolution model configuration (not shown).

![Computational times in seconds of different MAOOAM implementations: (a) time to compute a 10$^7$ timeunits trajectory with a low-order model configuration (36 variables). (b) time to compute a 10$^6$ timeunits trajectory with a higher-order model configuration (228 variables). \label{fig:benchmark}](timing_results.pdf)
 
In conclusion, `qgs` is a sufficiently fast Python implementation as compared to the other implementations of the MAOOAM model. In addition, it has the benefit of being more flexible, extensible, and easier to use in the general Python scientific ecosystem.

 
# Acknowledgements

This research has been partly supported by EUMETNET (Postprocessing module of the NWP Cooperation Programme).

# References
