# SNEWPY: Supernova Neutrino Early Warning Models for Python

[![DOI](https://zenodo.org/badge/221705586.svg)](https://zenodo.org/badge/latestdoi/221705586)
[![PyPI](https://img.shields.io/pypi/v/snewpy)](https://pypi.org/project/snewpy/)
![tests](https://github.com/SNEWS2/snewpy/actions/workflows/tests.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/snewpy/badge/?version=latest)](https://snewpy.readthedocs.io/en/latest/?badge=latest)

SNEWPY is a Python package for working with supernova neutrinos. It offers …

* … a simple and unified interface to hundreds of supernova simulations.
* … a large library of flavor transformations that relate neutrino fluxes produced in the supernova to those reaching a detector on Earth.
* … and a Python interface to SNOwGLoBES which lets you estimate and plot event rates in many different neutrino detectors.


## Installation

### For Users
Run `pip install snewpy` to install SNEWPY.

SNEWPY includes a large number of supernova models from different simulation groups. Since these models have a size of several 100 MB, they are not included in the initial install. Instead, after installing, run the following command to download models you want to use:

`python -c 'import snewpy; snewpy.get_models()'`

By default, they will be downloaded to a subdirectory named `SNEWPY-models/<model_name>/` in the current directory.

### For Developers

**Your contributions to SNEWPY are welcome!** For minor changes, simply submit a pull request. If you plan larger changes, it’s probably a good idea to open an issue first to coordinate our work.

To contribute, first clone the repository (`git clone https://github.com/SNEWS2/snewpy.git`), then make changes and install your modified version locally using `pip install .` from the base directory of the repository.
Once you’re happy with your changes, please submit a pull request.
Unit tests will run automatically for every pull request or you can run them locally using `python -m unittest python/snewpy/test/test_*.py`.

### Dependencies 

Some functionality of SNEWPY requires that [SNOwGLoBES](https://github.com/SNOwGLoBES/snowglobes) and its dependency, [GLoBES](https://www.mpi-hd.mpg.de/personalhomes/globes/) are installed.

<details>
<summary>Expand this to view sample instructions for installing SNOwGLoBES and GLoBES</summary>

This is a walkthrough to install GLoBES and SNOwGLoBES locally in the users home directory
(i.e. `~/opt/`). It uses `bash` syntax.

```bash
	cd ~
	mkdir opt
	cd opt
	wget https://www.mpi-hd.mpg.de/personalhomes/globes/download/globes-3.2.17.tar.gz
	tar -zxf globes-3.2.17.tar.gz
	cd globes-3.2.17/
	./configure --prefix=~/opt/globes-3.2.17-install --disable-binary
	make
	make install
	cd ~/opt/globes-3.2.17-install
	export GLB_DIR=${PWD}
	cd ..

	git clone https://github.com/SNOwGLoBES/snowglobes.git
	cd snowglobes
	export SNOWGLOBES=${PWD}
	cd src
	make
	make install
```
</details> 


## Usage and Documentation
Example scripts which show how SNEWPY can be used are available in the
`python/snewpy/scripts/` subfolder as well as notebooks in `doc/nb/`.
Most downloadable models also include a Jupyter notebook with simple usage examples.

A paper describing SNEWPY and the underlying physics is available at [arXiv:2109.08188](https://arxiv.org/abs/2109.08188).

For more, see the [full documentation on Read the Docs](https://snewpy.rtfd.io/).
Data from I. Tamborra, G. Raffelt, F. Hanke, H.-T. Janka, and B. Mueller *Neutrino emission characteristics and detection opportunities based on three-dimensional supernova simulations*, [Phys. Rev. D 90:045032, 2014](https://arxiv.org/abs/1406.0006).

Files taken from the [Garching Supernova archive](http://wwwmpa.mpa-garching.mpg.de/ccsnarchive/data/Tamborra2014/) for 20Msun and 27Msun, with permission to include them into snewpy.

NOTE: same file format as for Bollig_2016 models.
Model from L. Walk, I. Tamborra, H.-T. Janka, A. Summa, D. Kresse, [Neutrino emission characteristics of black hole formation in three-dimensional simulations of stellar collapse](https://arxiv.org/abs/1910.12971), Phys. Rev. D 101:123013, 2019.

Data taken from the [Garching Supernova archive](https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/data/Walk2019/) for 40Msun progenitor, with permission to use them within snewpy.
Models from Takami Kuroda using a full 3D magnetorotational CCSN simulation.

Included with SNEWPY with permission of the authors.

The citation to be used is:
* *Impact of a Magnetic Field on Neutrino-Matter Interactions in Core-collapse Supernovae* by Takami Koruda, [ApJ 906 (2021) 128](https://iopscience.iop.org/article/10.3847/1538-4357/abce61), [arXiv:2009.07733](https://arxiv.org/abs/2009.07733). 
CCSN neutrino models from Nakazato et al., 2013, 2015 and 2021. Data are available publicly
on [their website](http://asphwww.ph.noda.tus.ac.jp/snn/) and were reprocessed for SNEWPY.

The citation for use of the database is: *Supernova Neutrino Light Curves and
Spectra for Various Progenitor Stars: From Core Collapse to Proto-neutron Star
Cooling*, K. Nakazato, K. Sumiyoshi, H. Suzuki, T. Totani, H. Umeda, and S.
Yamada, [Astrophys. J. Supp. 205 (2013)
2](http://dx.doi.org/10.1088/0067-0049/205/1/2), [arXiv:1210.6841](http://arxiv.org/abs/1210.6841).

If the BH model with LS220 EOS is used, the citation is: *Spectrum of the
Supernova Relic Neutrino Background and Metallicity Evolution of Galaxies*, 
K. Nakazato, E. Mochida, Y. Niino, and H. Suzuki, 
[Astrophys. J. 804 (2015) 75](http://dx.doi.org/10.1088/0004-637X/804/1/75), [arXiv:1503.01236](http://arxiv.org/abs/1503.01236).

For the BH model with Togashi EOS, the citation is: *Numerical Study of Stellar Core
Collapse and Neutrino Emission Using the Nuclear Equation of State Obtained by the Variational Method*,
K. Nakazato, K. Sumiyoshi, and H. Togashi, 
[Publ. Astron. Soc. Jpn. 73 (2021) 639-651](https://doi.org/10.1093/pasj/psab026), [arXiv:2103.14386](http://arxiv.org/abs/2103.14386).The files in this folder are for Pair Instability Supernova models. 

The models are described in *Neutrino signal from pair-instability supernovae*, Warren P. Wright, Matthew S. Gilmer, Carla Fröhlich, James P. Kneller, [Phys. Rev. D96 (2017) 103008](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.96.103008), [arXiv:1706.08410](https://arxiv.org/abs/1706.08410).


There are two different progenitors, P150 and P250, and the neutrino spectra are computed using two different Equations of State. These EOS's radically change the weak neutrino emission but not the thermal. 

There is only one line of sight into each model. Flavor transformation due to MSW effects and decoherence are included. Self-interaction effects do not occur in PISN models - this was checked. Earth matter effects are not included. Both mass orderings are considered plus the no-oscillation case. 

The format of each file is the SNOwBLoBES format.
Data from Mirizzi et al., in particular the models that go into Figure 17. Models (s11.2c and s27.0c) taken from the Garching Supernova archive (https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/data/Bollig2016/) with permission for use in SNEWS2.0.

Reference: Mirizzi et al. Rivista del Nuovo Cimento Vol 39 N. 1-2 (2016)  
- doi:10.1393/ncr/i2016-10120-8
- arXiv:1508.00785
Data from O'Connor 2015, black hole forming simulations of a 40 solar mass progenitor from Woosley and Heger 2007 and the LS220 EOS.  

Reference: O'Connor ApJS 219 24 2015
- doi:10.1088/0067-0049/219/2/24
- arXiv:1411.7058
CCSN neutrino estimates from the 3D Fornax model of Vartanyan et al., taken from their [3D simulations index](https://www.astro.princeton.edu/~burrows/nu-emissions.3d/).

The citations to be used are:
* *A Successful 3D Core-Collapse Supernova Explosion Model* David Vartanyan, Adam Burrows, David Radice, Aaron Skinner, and Joshua Dolence, [MNRAS 482 (2019) 351](https://doi.org/10.1093/mnras/sty2585), [arXiv:1809.05106](https://arxiv.org/abs/1809.05106).
* *Three-Dimensional Supernova Explosion Simulations of 9-, 10-, 11-, 12-, and 13-M⊙ Stars*, Adam Burrows, David Radice, and David Vartanyan, [MNRAS 485 (2019) 3153](https://doi.org/10.1093/mnras/stz543), [arXiv:1902.00547](https://arxiv.org/abs/1902.00547).

The response of several neutrino detectors was also calculated by this group and can be found in *Neutrino Signals of Core-Collapse Supernovae in Underground Detectors* by Shaquann Seadrow, Adam Burrows, David Vartanyan, David Radice, and M. Aaron Skinner, [MNRAS 480 (2018) 4710](https://doi.org/10.1093/mnras/sty2164), [arXiv:1804.00689](https://arxiv.org/abs/1804.00689).
CCSN neutrino estimates from the 2D Fornax simulation extended to late times (4.0-4.5 seconds after core bounce) by Burrows and Vartanyan. See their [2D simulations index](https://www.astro.princeton.edu/~burrows/nu-emissions.2d/).

The citations to be used are:
* *Core-Collapse Supernova Explosion Theory* by Adam Burrows and David Vartanyan, [Nature 589 (2021) 29](https://www.nature.com/articles/s41586-020-03059-w), [arXiv:2009.14157](https://arxiv.org/abs/2009.14157).

Data from L. Walk, I. Tamborra, H.-T. Janka, A. Summa, [Identifying rotation in SASI-dominated core-collapse supernovae with a neutrino gyroscope](https://arxiv.org/abs/1807.02366), Phys. Rev. D 98:123001, 2018.

Data taken from the [Garching Supernova archive](https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/data/Walk2018/) for 15Msun progenitor, with permission to use them within snewpy.
Data from "Constraining properties of the next nearby core-collapse supernova with multi-messenger signals: multi-messenger signals" by Warren, MacKenzie; Couch, Sean; O'Connor, Evan; Morozova, Viktoriya.

1D FLASH simulations with STIR, for alpha_lambda = 1.23, 1.25, and 1.27.  Run with SFHo EOS, M1 with 12 energy groups.

For more information on these simulations, see Warren, Couch, O'Connor, & Morozova (arXiv:1912.03328) and Couch, Warren, & O'Connor (2020).

Includes the multi-messenger data from the STIR simulations. The filename indicates the turbulent mixing parameter a and progenitor mass m of the simulation.  Columns are time [s], shock radius [cm], explosion energy [ergs], electron neutrino mean energy [MeV], electron neutrino rms energy [MeV], electron neutrino luminosity [10^51 ergs/s], electron antineutrino mean energy [MeV], electron antineutrino rms energy [MeV], electron antineutrino luminosity [10^51 ergs/s], x neutrino mean energy [MeV], x neutrino rms energy [MeV], x neutrino luminosity [10^51 ergs/s], gravitational wave frequency from eigenmode analysis of the protoneutron star structure [Hz].  Note that the x neutrino luminosity is for one neutrino flavor - to get the total mu/tau neutrino and antineutrino luminosities requires multiplying this number by 4.
CCSN models from O'Connor and Ott, 2013 using 32 stellar models of metallicity with ZAMS masses from 12 to 120 solar masses, assuming two different equations of state. The models cover only the pre-explosion phase since the simulations did not explode.

Data downloaded from [sntheory.org](https://sntheory.org/M1prog) include:

* SNOwGLoBES fluence data for 5 ms bins.
* Luminosities and average energies in 1 ms bins.
* Spectra in 1 ms bins.

See *The Progenitor Dependence of the Preexplosion Neutrino Emission in Core-Collapse Supernovae* E. O'Connor and E. Ott, [Astrophys. J. 762 (2013) 126](https://iopscience.iop.org/article/10.1088/0004-637X/762/2/126), [arXiv:1207.1100](https://arxiv.org/abs/1207.1100).
The files in this folder are for two Type Ia models. 

The DDT model is described in *Neutrinos from type Ia supernovae: The deflagration-to-detonation transition scenario*, Warren P. Wright, Gautam Nagaraj, James P. Kneller, Kate Scholberg, Ivo R.  Seitenzahl, [Phys. Rev. D94 (2016) 025026](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.025026), [arXiv:1605.01408](https://arxiv.org/abs/1605.01408).  

There are 30 snapshots in time and for each snapshot there are 8 lines of sight through the explosion. Flavor transformation due to MSW effects and decoherence are included plus, the no-oscillations case. Self-interaction effects do not occur in Type Ia models. Earth matter effects are not included. Both mass orderings are considered. The format of each data file is the SNOwGLoBES format. 

The GCD model is described in *Neutrinos from type Ia supernovae: The gravitationally confined detonation scenario* Warren P. Wright, James P. Kneller, Sebastian T. Ohlmann, Friedrich K. Röpke, Kate Scholberg, Ivo R. Seitenzahl, [Phys. Rev. D95 (2017) 043006](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.95.043006), [arXiv:1609.07403](https://arxiv.org/abs/1609.07403).  

There are 64 snapshots in time and for each snapshot there are 10 lines of sight through the explosion. Flavor transformation due to MSW effects and decoherence are included, plus the no-oscillations case. Self-interaction effects do not occur in Type Ia models. Earth matter effects are not included. Both mass orderings are considered. The format of each data file is the SNOwGLoBES format. Flavor transformation due to MSW effects and decoherence are included. Self-interaction effects do not occur in Type Ia models. Earth matter effects are not included.
Data from Zha, O'Connor, Schneider 2021, failing core-collapse supernova simulation with a hybrid EOS including a hadron-quark phase transition.  

Reference: Zha, O'Connor, Schneider 2021 ApJ in press
- arXiv:2103.02268
CCSN neutrino models from the MPA Garching CCSN archive based on the paper by Sukhbold et al., 2015. The archive is available on
[their website](https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/data/SEWBJ_2015/index.html), but the data are private and available only upon request.
Note these are the results using the PROMETHEUS-VERTEX code https://ui.adsabs.harvard.edu/abs/2002A%26A...396..361R/abstract. 
The four models are also described in Appendix C of this paper https://arxiv.org/abs/2010.04728

The citation is: *Core-Collapse Supernovae from 9 to 120 Solar Masses Based on Neutrino-powered Explosions*, Tuguldur Sukhbold, T. Ertl, S. E. Woosley, Justin M. Brown, H.-T. Janka, [Astrophys. J. 821 (2016)
38](http://dx.doi.org/10.3847/0004-637X/821/1/38), [arXiv:1510.04643](http://arxiv.org/abs/1510.04643).

---
title: 'SNEWPY: A Data Pipeline from Supernova Simulations to Neutrino Signals'
tags:
  - Python
  - astronomy
  - supernova
  - neutrinos
authors:
  - name: Amanda L. Baxter
    affiliation: 1
  - name: Segev BenZvi
    orcid: 0000-0001-5537-4710
    affiliation: 2
  - name: Joahan Castaneda Jaimes
    affiliation: 3
  - name: Alexis Coleiro
    affiliation: 4
  - name: Marta Colomer Molla
    orcid: 0000-0003-1801-8121
    affiliation: 5
  - name: Damien Dornic
    affiliation: 6
  - name: Spencer Griswold
    orcid: 0000-0002-7321-7513
    affiliation: 2
  - name: Tomer Goldhagen
    affiliation: 7
  - name: Anne Graf
    affiliation: 8
  - name: Alec Habig
    affiliation: 9
    orcid: 0000-0002-1018-9383
  - name: Remington Hill
    affiliation: 10
  - name: Shunsaku Horiuchi
    orcid: 0000-0001-6142-6556
    affiliation: 11
  - name: James P. Kneller^[Corresponding author]
    orcid: 0000-0002-3502-3830
    affiliation: 8
  - name: Mathieu Lamoureux
    orcid: 0000-0002-8860-5826 
    affiliation: 12
  - name: Rafael F. Lang
    affiliation: 1
    orcid: 0000-0001-7594-2746
  - name: Massimiliano Lincetto
    orcid: 0000-0002-1460-3369
    affiliation: 13
  - name: Jost Migenda
    orcid: 0000-0002-5350-8049
    affiliation: 14
  - name: McKenzie Myers
    orcid: 0000-0002-2901-9173
    affiliation: 8
  - name: Evan O'Connor
    affiliation: 15
  - name: Andrew Renshaw
    affiliation: 16
    orcid: 0000-0003-2913-8057
  - name: Kate Scholberg
    orcid: 0000-0002-7007-2021
    affiliation: 17
  - name: Andrey Sheshukov
    affiliation: 18
    orcid: 0000-0001-5128-9279
  - name: Jeff Tseng
    affiliation: 19
    orcid: 0000-0003-1731-5853
  - name: Christopher Tunnell
    orcid: 0000-0001-8158-7795
    affiliation: 20
  - name: Navya Uberoi
    affiliation: 2
  - name: Arkin Worlikar
    affiliation: 21
affiliations:
  - name: Purdue University, West Lafayette, IN, USA
    index: 1
  - name: University of Rochester, Rochester, NY, USA
    index: 2
  - name: California Institute of Technology, Pasadena, CA, USA
    index: 3
  - name: Université de Paris, CNRS, AstroParticule et Cosmologie, Paris, France
    index: 4
  - name: Université Libre de Bruxelles, Brussels, Belgium
    index: 5
  - name: Aix Marseille Univ, CNRS/IN2P3, CPPM, Marseille, France
    index: 6
  - name: University of North Carolina - Chapel Hill, Chapel Hill, NC, USA
    index: 7
  - name: NC State University, Raleigh, NC, USA
    index: 8
  - name: University of Minnesota Duluth, Duluth, MN, USA
    index: 9
  - name: Laurentian University, Sudbury, ON, Canada
    index: 10
  - name: Virginia Tech, Blacksburg, VA, USA
    index: 11
  - name: INFN Sezione di Padova, Padova, Italy
    index: 12
  - name: Ruhr-Universität Bochum, Bochum, Germany
    index: 13
  - name: King’s College London, London, UK
    index: 14
  - name: Stockholm University, Stockholm, Sweden
    index: 15
  - name: University of Houston, Houston, TX, USA
    index: 16
  - name: Duke University, Durham, NC, USA
    index: 17
  - name: Joint Institute for Nuclear Research, Dubna, Russia
    index: 18
  - name: Oxford University, Oxford, UK
    index: 19
  - name: Rice University, Houston, TX, USA
    index: 20
  - name: Georgia Institute of Technology, Atlanta, GA, USA
    index: 21
date: 1 November 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/1538-4357/ac350f # <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal # <- The name of the AAS journal.
---


# Summary

Current neutrino detectors will observe hundreds to thousands of neutrinos
from a Galactic supernova, and future detectors will increase this yield by
an order of magnitude or more. With such neutrino data sets, the next
Galactic supernova will significantly increase our understanding of the
explosions of massive stars, nuclear physics under extreme conditions, and
the fundamental properties of neutrinos. However, there is a gulf
between supernova simulations and the corresponding signals in detectors,
making comparisons between theory and observation, as well as between
different detectors, very difficult. SNEWPY offers a unified interface for
hundreds of supernova simulations, a large library of flux transformations on
the way towards the detector, and an interface to SNOwGLoBES [@SNOwGLoBES],
allowing users to easily calculate and compare expected event rates from many supernova
models in many different neutrino detectors.

# Statement of need

SNEWPY is an open-source software package which bridges the gap between
simulations of supernova neutrinos and the corresponding signals (neutrino
events) one would expect from neutrino detectors here on Earth. The package,
written in Python, is built upon NumPy [@harris2020array] and SciPy
[@Virtanen:2019joe], and makes use of Astropy [@Astropy:2013muo;
@Price-Whelan:2018hus] for model I/O and unit conversions.

![Flowchart showing the complete SNEWPY pipeline. SNEWPY supports a wide variety of input formats and can output results as plots or as a Python dictionary for further analysis.\label{fig:flowchart}](snewpy-flowchart.pdf)

SNEWPY consists of three main modules that together form a complete
simulation pipeline (see \autoref{fig:flowchart}).
The first module, `snewpy.models`, interfaces with supernova simulation data
sets in different formats to extract the neutrino emission produced in the
supernova as a function of time, energy, angle, and neutrino flavor.
The `snewpy.flavor_transformation` module then
convolves the neutrino spectra with a prescription for neutrino flavor
transformation in the mantle of the star and during propagation to Earth.
The third module, `snewpy.snowglobes`, interfaces with SNOwGLoBES itself:
First, it can generate either a time series of neutrino spectra at Earth—the
“neutrinocurve”—or the spectral fluence. The module is then able to
run the generated data files through SNOwGLoBES, which computes the expected
event rates in different neutrino detector models, before collating the output
from SNOwGLoBES into a signal data file per detector per interaction channel.

Instead of using it as a complete simulation pipeline, SNEWPY can also be
integrated into other software thanks to its modular design.
For example, the supernova event generator sntools [@Migenda2021] recently
incorporated SNEWPY as a dependency to provide access to a broad range of
supernova models and flavor transformations.

In addition to the source code, SNEWPY comes with data from several hundred
simulations kindly provided by various modeling groups, a script for
generating a spectral fluence from an analytic prescription, and several
Jupyter notebooks illustrating its capabilities. While SNEWPY has been
developed explicitly for the SuperNova Early Warning System, SNEWS 2.0
[@SNEWS:2020tbu], its object-oriented design makes the addition of new
supernova models and flavor transformations straightforward. We expect that it
will prove broadly useful to modelers and theorists interested in what
neutrino detectors will observe from a supernova simulation, as well as
experimentalists wishing to evaluate the sensitivity of their detector to
supernova neutrinos. 

# Acknowledgements

This work is supported by the National Science Foundation “Windows on the
Universe: the Era of Multi-Messenger Astrophysics” Program: “WoU-MMA:
Collaborative Research: A Next-Generation SuperNova Early Warning System for
Multimessenger Astronomy” through Grant Nos. 1914448, 1914409, 1914447,
1914418, 1914410, 1914416, and 1914426.
This work is also supported at NC State by U.S. Department of Energy grant
DE-FG02-02ER41216, at Stockholm University by the Swedish Research Council
(Project No. 2020-00452), and at King’s College London by STFC.

# References
# SNEWPY Usage Examples

The Jupyter notebooks in this directory contain different examples for how to use SNEWPY. Please see also the example notebooks for each supernova model, which are included in the same directory as the model files when downloading them via `python -c 'import snewpy; snewpy.get_models()'`.

## AnalyticFluence

This notebook demonstrates how to use the `Analytic3Species` class from `snewpy.models` to create an analytic supernova model by specifying the luminosity, mean energy and mean squared energy for three neutrino flavors.

## FlavorTransformation

This notebook demonstrates the flavor transformations available in `snewpy.flavor_transformation`. It was used to produce many of the figures in the SNEWPY ApJ paper.

## SNOwGLoBES_models

This notebook demonstrates how to use the `SNOwGLoBES` class in `snewpy.models`, which can be used with the `Type_Ia` and `PISN` model files that are available for download through SNEWPY.

## SNOwGLoBES_usage

This notebook demonstrates how to use SNEWPY’s `snewpy.snowglobes` module to interact with SNOwGLoBES.
