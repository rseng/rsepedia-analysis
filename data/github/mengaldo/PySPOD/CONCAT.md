<p align="center">
  <a href="http://mengaldo.github.io/PySPOD/" target="_blank" >
    <img alt="Python Spectral Proper Orthogonal Decomposition" src="readme/PySPOD_logo2.png" width="200" />
  </a>
</p>

<p align="center">
  <a href="https://doi.org/10.21105/joss.02862" target="_blank">
    <img alt="JOSS Paper" src="https://joss.theoj.org/papers/10.21105/joss.02862/status.svg">
  </a>

  <a href="https://github.com/mengaldo/PySPOD/LICENSE" target="_blank">
    <img alt="Software License" src="https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square">
  </a>

  <a href="https://badge.fury.io/py/pyspod">
    <img src="https://badge.fury.io/py/pyspod.svg" alt="PyPI version" height="18">
  </a>

  <a href="https://travis-ci.com/mengaldo/PySPOD" target="_blank">
    <img alt="Build Status" src="https://travis-ci.com/mengaldo/PySPOD.svg?token=sY467gr18pmboZ16AN1x&branch=main">	  
  </a>

  <a href="https://coveralls.io/github/mengaldo/PySPOD?branch=main" target="_blank">
    <img src="https://coveralls.io/repos/github/mengaldo/PySPOD/badge.svg?branch=main" alt="Coverage Status" />
  </a>

  <a href="https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=mengaldo/PySPOD&amp;utm_campaign=Badge_Grade">
    <img src="https://app.codacy.com/project/badge/Grade/7ac24e711aea47df806ad52ab067e3a6"/>
  </a>
</p>

**PySPOD**: Python Spectral Proper Orthogonal Decomposition

## Table of contents

  * [Description](#description)
  * [Installation and dependencies](#installation-and-dependencies)
    * [Installing via PIP](#installing-via-pip)
    * [Installing from source](#installing-from-source)
  * [Documentation](#documentation)
  * [Testing](#testing)
  * [References](#references)
  * [Recent works with PySPOD](#recent-works-with-pyspod)
  * [Authors and contributors](#authors-and-contributors)
  * [License](#license)

## Description
**PySPOD** is a Python package that implements the so-called **Spectral Proper Orthgonal Decomposition** whose name was first conied by [(Picard and Delville 2000)](#picard-and-delville-2000), and goes back to the original work by [(Lumley 1970)](#lumley-1970). The implementation proposed here follows the original contributions by [(Towne et al. 2018)](#towne-et-al-2018), [(Schmidt and Towne 2019)](#schmidt-and-towne-2019).

**Spectral Proper Orthgonal Decomposition (SPOD)** has been extensively used in the past few years to identify spatio-temporal coherent patterns in a variety of datasets, mainly in the fluidmechanics and climate communities. In fluidmechanics it was applied to jets [(Schmidt et al. 2017)](#schmidt-et-al-2017), wakes [(Araya et al. 2017)](#araya-et-al-2017), and boundary layers [(Tutkun and George 2017)](#tutkun-and-george-2017), among others, while in weather and climate it was applied to ECMWF reanalysis datasets under the name Spectral Empirical Orthogonal Function, or SEOF, [(Schmidt et al. 2019)](#schmidt-et-al-2019).

The SPOD approach targets statistically stationary problems and involves the decomposition of the cross-spectral density tensor. This means that the SPOD leads to a set of spatial modes that oscillate in time at a single frequency and that optimally capture the variance of an ensemble of stochastic data [(Towne et al. 2018)](#towne-et-al-2018). Therefore, given a dataset that is statistically stationary, one is able to capture the optimal spatio-temporal coherent structures that explain the variance in the dataset. 

This can help identifying relations to multiple variables or understanding the reduced order behavior of a given phenomenon of interest and represent a powerful tool for the data-driven analysis of nonlinear dynamical systems. The SPOD approach shares some relationships with the dynamic mode decomposition (DMD), and the resolvent analysis,  [(Towne et al. 2018)](#Towne-et-al-2018), that are also widely used approaches for the data-driven analysis of nonlinear systems. SPOD can be used for both experimental and simulation data, and a general description of its key parameters can be found in [(Schmidt and Colonius 2020)](#schmidt-and-colonius-2020).  

In this package we implement three version of SPOD 

  - SPOD_low_storage: that is intended for large RAM machines or small datasets
  - SPOD_low_ram: that is intended for small RAM machines or large datasets, and 
  - SPOD_streaming: that is the algorithm presented in [(Schmidt and Towne 2019)](schmidt-and-towne-2019).

To see how to use the **PySPOD** package and its user-friendly interface, you can look at the [**Tutorials**](tutorials/README.md). 

## Installation and dependencies
**PySPOD** requires the following Python packages: 
`numpy`, `scipy`, `matplotlib`, `xarray`, `netcdf4`, `h5py`, `psutil`, `tdqm`, `future`, `ffmpeg`, `sphinx` (for the documentation). 
Some of the *Climate tutorials*, additionally need `ecmwf_api_client` and `cdsapi`. 

The code is developed and tested for Python 3 only. 
It can be installed using `pip` or directly from the source code.

<!-- NOTE: 
  - to properly install netcdf4, you might need to have a local installation of `hdf5`. 
  - to be able to use the ffmpeg functionalities of the library (generating video of your data), you need a local installation of `ffmpeg` libraries.
  -->
	
### Installing via PIP
Mac and Linux users can install pre-built binary packages using pip.
To install the package just type: 
```bash
    > pip install pyspod 
```
To uninstall the package:
```bash
    > pip uninstall pyspod
```

### Installing from source
The official distribution is on GitHub, and you can clone the repository using
```bash
> git clone https://github.com/mengaldo/PySPOD
```

To install the package just type:
```bash
> python setup.py install
```

To uninstall the package you have to rerun the installation and record the installed files in order to remove them:

```bash
> python setup.py install --record installed_files.txt
> cat installed_files.txt | xargs rm -rf
```

## Get started with a simple analysis
**PySPOD** comes with an extensive suite of [**Tutorials**](tutorials/README.md). 
You can browse the [**Tutorials**](tutorials/README.md) to explore the capabilities 
and various functionalities of the library. However, if you want to get started 
quickly, after you installed the library you can simply copy the following script 
into a file `your_script.py`, and run it with Python 3 (e.g. from a terminal window, 
the run command would look like `> python3 your_script.py`).

```python
import os
import xarray as xr
import numpy  as np

# Import library specific modules
from pyspod.spod_low_storage import SPOD_low_storage
from pyspod.spod_low_ram     import SPOD_low_ram
from pyspod.spod_streaming   import SPOD_streaming
import pyspod.utils_weights as utils_weights


# Let's create some 2D syntetic data

# -- define spatial and time coordinates
x1 = np.linspace(0,10,100) 
x2 = np.linspace(0, 5, 50) 
xx1, xx2 = np.meshgrid(x1, x2)
t = np.linspace(0, 200, 1000)

# -- define 2D syntetic data
s_component = np.sin(xx1 * xx2) + np.cos(xx1)**2 + np.sin(0.1*xx2)
t_component = np.sin(0.1 * t)**2 + np.cos(t) * np.sin(0.5*t)
p = np.empty((t_component.shape[0],)+s_component.shape)
for i, t_c in enumerate(t_component):
    p[i] = s_component * t_c


# Let's define the required parameters into a dictionary
params = dict()

# -- required parameters
params['time_step'   ] = 1              # data time-sampling
params['n_snapshots' ] = t.shape[0]     # number of time snapshots (we consider all data)
params['n_space_dims'] = 2              # number of spatial dimensions 
params['n_variables' ] = 1 		# number of variables
params['n_DFT'       ] = 100          	# length of FFT blocks (100 time-snapshots)

# -- optional parameters
params['overlap'          ] = 0           # dimension block overlap region
params['mean_type'        ] = 'blockwise' # type of mean to subtract to the data
params['normalize_weights'] = False       # normalization of weights by data variance
params['normalize_data'   ] = False       # normalize data by data variance
params['n_modes_save'     ] = 3           # modes to be saved
params['conf_level'       ] = 0.95        # calculate confidence level
params['reuse_blocks'     ] = True        # whether to reuse blocks if present
params['savefft'          ] = False       # save FFT blocks to reuse them in the future (saves time)
params['savedir'          ] = os.path.join('results', 'simple_test') # folder where to save results


# Initialize libraries for the low_storage algorithm
spod = SPOD_low_storage(p, params=params, data_handler=False, variables=['p'])

# and run the analysis
spod.fit()


# Let's plot the data
spod.plot_2D_data(time_idx=[1,2])
spod.plot_data_tracers(coords_list=[(5,2.5)], time_limits=[0,t.shape[0]])


# Show results
T_approx = 10 # approximate period = 10 time units
freq = spod.freq
freq_found, freq_idx = spod.find_nearest_freq(freq_required=1/T_approx, freq=freq)
modes_at_freq = spod.get_modes_at_freq(freq_idx=freq_idx)
spod.plot_eigs()
spod.plot_eigs_vs_period(freq=freq, xticks=[1, 7, 30, 365, 1825])
spod.plot_2D_modes_at_frequency(
	freq_required=freq_found, freq=freq, x1=x2, x2=x1, modes_idx=[0,1], vars_idx=[0])
```
You can change `SPOD_low_storage` to `SPOD_low_ram` and `SPOD_streaming`, 
to run the other two SPOD algorithms available.

## Documentation
**PySPOD** uses [Sphinx](http://www.sphinx-doc.org/en/stable/) for code documentation. 
You can view the documentation online [here](http://mengaldo.github.io/PySPOD/). 
If you want to build the documentation locally on your computer, you can do so 
by:

```bash
> cd docs
> make html
```

This will generate a `docs/build/html` folder, where you can find an `index.html` file. 
Open it with your browser and explore the documentation locally.

## Testing
Regression tests are deployed using Travis CI, that is a continuous intergration framework. 
You can check out the current status of **PySPOD** [here](https://travis-ci.org/mengaldo/PySPOD).

IF you want to run tests locally, you can do so by:

```bash
> cd tests/
> pytest -v
```

## References

#### (Lumley 1970) 
*Stochastic Tools in Turbulence.* [[DOI](https://www.elsevier.com/books/stochastic-tools-in-turbulence/lumey/978-0-12-395772-6?aaref=https%3A%2F%2Fwww.google.com)]

#### (Picard and Delville 2000) 

*Pressure velocity coupling in a subsonic round jet.*
[[DOI](https://www.sciencedirect.com/science/article/abs/pii/S0142727X00000217)]

#### (Tutkun and George 2017)

*Lumley decomposition of turbulent boundary layer at high Reynolds numbers.*
[[DOI](https://aip.scitation.org/doi/10.1063/1.4974746)]

#### (Schmidt et al 2017) 

*Wavepackets and trapped acoustic modes in a turbulent jet: coherent structure eduction and global stability.*
[[DOI](https://doi.org/10.1017/jfm.2017.407)]

#### (Araya et al 2017)

*Transition to bluff-body dynamics in the wake of vertical-axis wind turbines.*
[[DOI]( https://doi.org/10.1017/jfm.2016.862)]

#### (Taira et al 2017) 

*Modal analysis of fluid flows: An overview.*
[[DOI](https://doi.org/10.2514/1.J056060)]

#### (Towne et al 2018)

*Spectral proper orthogonal decomposition and its relationship to dynamic mode decomposition and resolvent analysis.*
[[DOI]( https://doi.org/10.1017/jfm.2018.283)]

#### (Schmidt and Towne 2019)

*An efficient streaming algorithm for spectral proper orthogonal decomposition.*
[[DOI](https://doi.org/10.1016/j.cpc.2018.11.009)]

#### (Schmidt et al 2019)

*Spectral empirical orthogonal function analysis of weather and climate data.*
[[DOI](https://doi.org/10.1175/MWR-D-18-0337.1)]

#### (Schmidt and Colonius 2020)

*Guide to spectral proper orthogonal decomposition.*
[[DOI](https://doi.org/10.2514/1.J058809)]

## Recent works with **PySPOD**

Please, [contact me](mailto:gianmarco.mengaldo@gmail.com) if you used PySPOD for a publication and you want it to be advertised here.

- A. Lario, R. Maulik, G. Rozza, G. Mengaldo, [Neural-Network learning of SPOD latent space](https://arxiv.org/abs/2110.09218)

## Authors and contributors

**PySPOD** is currently developed and mantained by

  * [G. Mengaldo](mailto:mpegim@nus.edu.sg), National University of Singapore (Singapore).

Current active contributors include:

  * [R. Maulik](https://romit-maulik.github.io), Argonne National Laboratory (US).
  * [A. Lario](https://www.math.sissa.it/users/andrea-lario), SISSA (Italy)
  
## How to contribute

Contributions improving code and documentation, as well as suggestions about new features are more than welcome!

The guidelines to contribute are as follows: 
1. open a new issue describing the bug you intend to fix or the feature you want to add.
2. fork the project and open your own branch related to the issue you just opened, and call the branch `fix/name-of-the-issue` if it is a bug fix, or `feature/name-of-the-issue` if you are adding a feature.
3. ensure to use 4 spaces for formatting the code.
4. if you add a feature, it should be accompanied by relevant tests to ensure it functions correctly, while the code continue to be developed.
5. commit your changes with a self-explanatory commit message. 
6. push your commits and submit a pull request. Please, remember to rebase properly in order to maintain a clean, linear git history.

[Contact me](mailto:mpegim@nus.edu.sg) by email for further information or questions about **PySPOD** or ways on how to contribute. 


## License

See the [LICENSE](LICENSE.rst) file for license rights and limitations (MIT).
---
title: 'PySPOD: A Python package for Spectral Proper Orthogonal Decomposition (SPOD)'
tags:
  - Python
  - dynamical systems
  - nonlinear dynamics
  - data-driven dynamics
  - data mining
authors:
  - name: Gianmarco Mengaldo
    orcid: 0000-0002-0157-5477
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Romit Maulik
    affiliation: "2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Department of Mechanical Engineering, National University of Singapore (SG)
   index: 1
 - name: Argonne Leadership Computing Facility, Argonne National Laboratory (USA) 
   index: 2
date: 6 November 2020
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
<!-- aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal. -->
---

# Summary

Large unstructured datasets may contain complex coherent patterns that 
evolve in time and space, and that the human eye cannot grasp. These 
patterns are frequently essential to unlock our understanding of complex 
systems that can arise in nature, such as the evolution of the atmosphere 
in the short (weather prediction) and long term (climate prediction), 
the behavior of turbulent flows, and the dynamics of plate tectonics, 
among several others. Identifying these coherent structures can 
prove crucial to facilitate the construction of modeling tools that can 
help anticipate scenarios that would not otherwise be predictable.

Within this context, dynamical systems theory, complemented with recent 
developments in machine learning and data mining tools, is achieving tremendous 
advances in our ability to acquire actionable information from complex 
data. Singular-value decomposition based techniques, in particular, are 
a promising area that is gaining popularity, due to their links to reduced 
order modeling and dynamical systems. Also, these techniques can be 
used in the context of machine learning as additional inputs to the learning 
architecture, thereby augmenting the dataset and possibly helping in 
the interpretability of the results. 

Several variants of singular-value decomposition (SVD) based techniques 
have been proposed in the literature, this library provides efficient 
implementations of the spectral proper orthogonal decomposition 
(SPOD) [@lumley1970; @towne2017]. SPOD is also referred to as spectral 
empirical orthogonal function (SEOF) in the weather and climate community 
[@schmidt2019a]. SPOD differs from other SVD-based techniques as it is 
derived from a standard (space-time) POD problem for stationary data and 
leads to modes that are (i) time harmonic and oscillate at a single frequency, 
(ii) are coherent in both time and space, (iii) optimally represent the space-time
statistical variability of the underlying stationary random processes, and 
(iv) are both spatially and space-time orthogonal [@schmidt2020]. 
We note that the `PySPOD` implements the Python counterpart of the Matlab 
code [@schmidt-code], with the addition of the streaming algorithm outlined 
by @schmidt2019b. We also acknowledge that there exist other two Python 
packages implementing SPOD. The first, by @spod-code-jburrows, is also a 
Python counterpart of the Matlab code of @schmidt-code. However, our 
implementation provides extensive post-processing capabilities, testing, 
and tutorial. It also adds the streaming version [@schmidt2019b], that 
is not present in @spod-code-jburrows. Similar differences exist between 
`PySPOD` and the Python package presented in @spod-code-loiseau.

# Capabilities 

`PySPOD` is a modular Python package that implements three different 
variants of SPOD, (i) a low storage [@towne2017; @schmidt2019a], 
(ii) a low RAM [@towne2017; @schmidt2019a], and (iii) a streaming version 
[@schmidt2019b]. The three versions differ in terms of I/O and RAM requirements. 
The low storage version allows faster computations, and it is intended for small 
datasets, or high RAM machines. The low RAM version can handle 
large datasets, but it is typically slower than the low storage counterpart. 
The streaming version is a streaming implementation of SPOD.
The API to the library offers a flexible and user-friendly experience, and 
the library can be complemented with additional SPOD algorithms in an easy-to-implement
way. The structure of the library and the use of Python enable efficient 
interfacing with low level and highly optimized libraries (written in C 
or Fortran) for the calculation of e.g. the fast Fourier transform, eigenvalue 
decomposition, and other linear algebra operations. Users can also take advantage 
of the ready-to-use postprocessing tools offered, and they can easily extend 
the postprocessing functionalities to suit their own needs. 

`PySPOD` is designed to be used in different fields of engineering and applied 
science, including weather and climate, fluid mechanics, seismology, among others.
It can be used as a production code, for the analysis of large datasets, as well 
as for experimenting on smaller problems. Users can be students and experts alike.
For an overview of the guidelines one should follow when using SPOD, the reader 
can refer to @schmidt2020.

In \autoref{fig:MEI}, we show the application of this package to identify 
the Multivariate ENSO Index from ECMWF reanalysis datasets (E20C in particular), 
where we used monthly averages of (i) mean sea level pressure (MSL), (ii) Zonal 
component of the surface wind (U10), (iii) Meridional component of the surface 
wind (V10), (iv) Sea surface temperature (SST),(v) 2-meter temperature (T2M), 
and (vi) Total cloud cover (TCC). \autoref{fig:MEI} shows the leading 
modes of the meridional component of the surface wind (left), and of the mean 
seal-level pressure (right). It is possible to appreciate a possible coupling 
between ENSO and the vortices over West Antarctica (that in turn could affect 
the height of the ice shelves [@paolo2018]). For more detail regarding this 
simulation, the interested reader can refer to @schmidt2019a.

![Identification of the Multivariate ENSO Index (MEI) from ECMWF reanalysis data.\label{fig:MEI}](../readme/MEI.png)



# Acknowledgements
G. Mengaldo wants to thank Oliver T. Schmidt for fruitful discussions. 
We also thank the reviewers who helped substantially improve the software package. 



# References
## Tutorials

We provide several tutorials that cover the main features of the PySPOD library. 
These are organized in the form of `jupyter-notebooks`, along with their plain 
`python` implementation.

In particular, we divided the tutorials in such a way that they cover different 
functionalities of the library and practical application areas.

### Basic

#### [Tutorial: Basic 1](basic/methods_comparison/methods_comparison.ipynb)

In this tutorial we give an introduction to the main functionalities 
of the package, by definining a 2D dataset and analyzing it via the 
three SPOD algorithms implemented. In this case, we load the entire 
data in RAM and pass it to the constructor of the SPOD class.

#### [Tutorial: Basic 2](basic/methods_comparison_file/methods_comparison_file.ipynb)

In this tutorial we give an introduction to the main functionalities 
of the package, by definining a 2D dataset and analyzing it via the 
three SPOD algorithms implemented. In contrast to [Tutorial: Basic 1](#tutorial-basic-1), 
in this tutorial we highlight how one can define a function to read 
data and pass it to the constructor of the SPOD class, thereby allowing 
for a reduced use of RAM (for large datasets).

### Climate 

#### [Tutorial: 2D Multivariate ENSO Index](climate/ERA20C_MEI_2D/ERA20C_MEI_2D.ipynb)

This tutorial shows how to download data from an ECMWF reanalysis dataset (ERA20C), 
and use **PySPOD** to identify spatio-temporal coherent structured in multivariate 
2D data. In particular, we seek to identify the multivariate ENSO index (MEI). 
The data is composed by the following monthly-averaged variables: mean sea level 
pressure (MSL), zonal component of the surface wind (U10), meridional component 
of the surface wind (V10), sea surface temperature (SST), 2-meter temperature 
(T2M), and, total cloud cover (TCC), on a 2D longitude-latitude grid.  

#### [Tutorial: 3D Quasi-Bienniel Oscillation](climate/ERA20C_QBO_3D/ERA20C_QBO_3D.ipynb)

This tutorial shows how to download data from an ECMWF reanalysis dataset (ERA20C), 
and use **PySPOD** to identify spatio-temporal coherent structured in univariate 
3D data. In particular, we seek to identify the Quasi-Bienniel Oscillation (QBO). 
The data is composed by the monthly-averages of the zonal-mean zonal winds 
on a 3D longitude, latitude, pressure-levels grid.

#### [Tutorial: 2D ERA5 Mean-Sea Level Pressure](climate/ERA5_MSLP_2D/ERA5_MSLP_2D.ipynb)

This tutorial shows how to download data from an ECMWF reanalysis dataset (ERA5), 
and use **PySPOD** to identify spatio-temporal coherent structured in univariate 
2D data. In particular, we seek to identify spatio-temporal coherent structure in 
high-resolution mean-sea level pressure data from the ERA5 dataset.

#### [Tutorial: 2D NAM Relative Humidity](climate/NAM_2D/NAM_2D.ipynb)

This tutorial explores the NAM dataset provided by NOAA, and in particular, the daily 
relative humidity reanalysis data for a period of ten years (2008-10-28) to (2018-09-20). 
While we use the first few years worth of data for a quick assessment, the readers are 
encouraged to increase the number of snapshots.

### Fluidmechanics 

#### [Tutorial: 2D Jet](fluidmechanics/jet_2D/jet_2D.ipynb)

This tutorial shows a simple 2D application to a turbulent jet, where the variable 
studied is pressure.

### Earthquakes 

#### [Tutorial: 2D Slip Potency](earthquakes/slip_potency_2D/slip_potency_2D.ipynb)

This tutorial shows a simple 2D application seismic data, where the variable studied 
is the slip potency.
