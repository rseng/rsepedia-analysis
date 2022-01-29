[![Build Status](https://github.com/temken/obscura/workflows/Build%20Status/badge.svg)](https://github.com/temken/obscura/actions)
[![codecov](https://codecov.io/gh/temken/obscura/branch/master/graph/badge.svg)](https://codecov.io/gh/temken/obscura)
[![Documentation Status](https://readthedocs.org/projects/obscura/badge/?version=latest)](https://obscura.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

# obscura - Direct detection of dark matter with nucleus and electron recoil experiments

[![status](https://joss.theoj.org/papers/fd8076268036956d3bf08193c4fc2db9/status.svg)](https://joss.theoj.org/papers/fd8076268036956d3bf08193c4fc2db9)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5665890.svg)](https://doi.org/10.5281/zenodo.5665890)

A modular C++ tool and library for dark matter direct detection computations for both nuclear and electron recoil experiments.

<img src="paper/FlowChart.png" width="500">

You can find more detailed documentation of *obscura* [here](https://obscura.readthedocs.io/en/latest/index.html). The documentation contains e.g. a [guide to get started](https://obscura.readthedocs.io/en/latest/01_Getting_Started.html) and a [list of all included experiments](https://obscura.readthedocs.io/en/latest/08_Experiments.html).

## CITATION

If you decide to use this code, or if you want to add a reference to it, please cite the latest archived version,

> Emken, T., 2021, obscura - A C++ library for dark matter detection computations [Code, v1.0.0] [[DOI:10.5281/zenodo.5665890]](https://zenodo.org/record/5665890).

<details><summary>Bibtex entry</summary>
<p>

```
@software{obscura,
  author = {Emken, Timon},
  title = {{obscura - A C++ library for dark matter detection computations [Code, v1.0.0]}},
  year         = {2021},
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {DOI:10.5281/zenodo.5665890},
  url          = {https://doi.org/10.5281/zenodo.5665890},
  howpublished={The code can be found under \url{https://github.com/temken/obscura}.}
}
```

</p>
</details>

## VERSION HISTORY

- 10.11.2021: Release of version 1.0.0
- 23.02.2021: Release of version 0.1.0

## AUTHORS & CONTACT

The author of *obscura* is [Timon Emken](https://timonemken.com/).

For questions, support, bug reports, or other suggestions, please contact [timon.emken@fysik.su.se](mailto:timon.emken@fysik.su.se) or open an [issue](https://github.com/temken/obscura/issues).


## LICENSE

This project is licensed under the MIT License - see the LICENSE file.
---
title: 'obscura: A modular C++ tool and library for the direct detection of (sub-GeV) dark matter via nuclear and electron recoils'
tags: 
  - c++
  - astroparticle physics
  - dark matter
  - direct detection
authors: 
  - name: Timon Emken
    orcid: 0000-0002-4251-2229
    affiliation: 1
affiliations: 
  - name: The Oskar Klein Centre, Department of Physics, Stockholm University, AlbaNova, SE-10691 Stockholm, Sweden 
    index: 1
date: 15 July 2021
bibliography: paper.bib
---

# Summary

The observation of a large number of gravitational anomalies on astrophysical and cosmological scales have convinced us that the majority of matter in the Universe is invisible [@Bertone:2004pz;@Bertone:2018krk].
This *dark matter* (DM) must be fundamentally different from the known matter we can describe using the Standard Model of Particle Physics (SM).
Its only established property is that it interacts gravitationally and indeed dominates the gravitational potential of galaxies and galaxy clusters.
One of the leading hypothesis is that DM is made up of one or more new particles and that galaxies such as our Milky Way are embedded in gigantic haloes of these as of yet undetected particles.
Our planet would at any moment be penetrated by a stream of these particles without much of an effect.
If these dark particles interact with nuclei and/or electrons via some new force besides gravity, they would on occasion collide with a terrestrial particle.
*Direct detection experiments* search for these kind of interactions and aim to observe DM events within a detector caused by an interaction with target nuclei [@Goodman:1984dc;@Drukier:1986tm;@Wasserman:1986hh] or electrons [@Kopp:2009et;@Essig:2011nj].
These experiments are typically placed deep underground to shield them from possible backgrounds, e.g., due to cosmic rays.

In order to interpret the outcome of direct detection experiments, we need to make predictions for the expected events caused by the incoming DM particles.
In all cases, this requires making a number of assumptions about the possible particle attributes of DM (e.g., mass and interaction strength) and the properties of the galactic DM halo (e.g. the local DM density and their energy distribution) [@Lewin:1995rx;@DelNobile:2021icc].

`obscura` is a tool to make quantitative predictions for direct DM searches, analyse experimental data, and derive, e.g., exclusion limits as seen in \autoref{fig:constraints}.
`obscura` can e.g. be used to compute the expected event rates in terrestrial detectors looking for rare interactions between the DM and nuclei or electrons.
There are many different experimental techniques and targets proposed and applied for direct detection experiments [@Griffin:2019mvc].
Additionally, due to our ignorance about the particle physics of DM there exists a plethora of viable assumptions and models.
The vast variety of viable assumptions is reflected by the modular, polymorphic structure of all modules of the `obscura` library which allows to easily extend `obscura`'s functionality to the users' new idea on the fundamental nature of DM particles, or on a new detection technology.
For example, the library can handle any kind of DM particles of any mass, provided that the scattering is well-described by non-relativistic dynamics, and that the differential (nucleus and/or electron) scattering cross sections depend only on the momentum transfer, the relative speed between DM and target, and at most one additional dynamic parameter such as the center-of-mass energy or the local temperature of the target.
Furthermore, a generic structure also allows applications of (a subset of) the `obscura` classes in a variety of DM research projects even beyond the context of direct detection, e.g., to compute DM capture rates in the Sun [@Emken:2021lgc].

For more details on `obscura` and its implementation in C++, we refer to the [documentation](https://obscura.readthedocs.io)[^1].

[^1]: The latest version of the documentation can be found under [https://obscura.readthedocs.io](https://obscura.readthedocs.io).

![Excluded regions (90% confidence level) of the DM parameter space given by the $(m_\mathrm{DM},\sigma_i)$ plane, where $m_\mathrm{DM}$ is the assumed DM mass and $\sigma_i$ is the interaction cross section with target $i$. For comparison, the dashed lines denote the official results published by the experimental collaborations. Some of the `obscura` results are conservative due to a simplified analysis.\label{fig:constraints}](obscura_DD_Constraints.png){ width=85% } 



# The modular structure of direct detection computations

Making predictions and performing analyses for direct detection experiments involves methods and results from statistics, astrophysics, particle physics, nuclear and atomic physics, and condensed matter physics.
For each of these fields, we need to make choices and assumptions which will affect our interpretation of DM searches.

As an example, let us look at the energy spectrum of DM induced ionization events as derived by @Essig:2015cda.
\begin{equation}
 \frac{\mathrm{d} R_\mathrm{ion}}{\mathrm{d} E_e} = N_T \frac{\rho_\chi}{m_\mathrm{DM}}\sum_{n,\ell} \int \mathrm{d}q^2\int \mathrm{d}v\; v f_\chi(v) \frac{1}{4E_e}\frac{\mathrm{d}\sigma_e}{\mathrm{d}q^2} \left|f_\mathrm{ion}^{n\ell}(q,E_e)\right|^2\, .
\end{equation}
The DM mass $m_\mathrm{DM}$ and the differential DM-electron scattering cross section $\frac{\mathrm{d}\sigma_e}{\mathrm{d}q^2}$ are defined by the assumed particle physics of the hypothetical DM particle the experiment is probing.
The velocity distribution $f_\chi(v)$ and the local DM energy density $\rho_\chi$ are important inputs from astrophysics and cosmology.
Lastly, the ionization form factor $f_\mathrm{ion}^{n\ell}(q,E_e)$ encapsulates the atomic physics of the electronic bound states and describes the probability of an electron with quantum numbers $(n\ell)$ to get ionized by an incoming DM particle.
As we can see, the evaluation of this expression for the electron recoil spectrum is highly modular combining inputs from various fields of research.
This modularity should be reflected in the structure of corresponding research software.

It is our ambition for the `obscura` code that the basic functionality does not rely on specific choices and that no particular assumption is hard-coded.
Instead the basic code's setup is polymorphic and written in terms of generic base classes widely agnostic to specific assumptions.
Of course, a number of standard ideas and models are implemented as derived classes, which also illustrate the usage of the base classes.
The classes are described in more detail in the [documentation](https://obscura.readthedocs.io).

![The class structure of `obscura`.\label{fig:flowchart}](FlowChart.png){ width=70% } 

External research software can use `obscura` by implementing its classes following the dependencies indicated by the flow chart of \autoref{fig:flowchart} and computing standard quantities in the context of direct detection of dark matter.
In addition, these classes are meant to be general-purpose and can be applied in other contexts depending on the research project's main objectives.
It is also possible to exploit the polymorphic structure and extend its functionality by creating new derived classes based on the users' own ideas.
As a final benefit of the polymorphic structure, any research software that is formulated entirely in terms of the abstract base classes can later on be used with any derived classes and allows analyses and research for a broad range of alternative assumptions without changing the core of the scientific code. 
As an example, the `DaMaSCUS-SUN` code uses `obscura` in the context of Monte Carlo simulations [@Emken:2021lgc;@Emken2021].

# Statement of need

For the interpretation of past and future direct searches for DM particles, it is important to be able to provide accurate predictions for event rates and spectra under a variety of possible and viable assumptions in a computationally efficient way.
While there exists a few tools to compute DM induced nuclear recoil spectra, such as DDCalc [@GAMBITDarkMatterWorkgroup:2017fax;@GAMBIT:2018eea] or WimPyDD [@Jeong:2021bpl], `obscura` is not limited to nuclear targets.
Instead its main focus lies on sub-GeV DM searches probing electron recoils which typically requires methods from atomic and condensed matter physics [@Essig:2015cda;@Catena:2019gfa;@Catena:2021qsr].
In the context of sub-GeV DM searches, new ideas such as target materials or detection techniques are being proposed regularly, and the theoretical modelling of these are getting improved continuously [@Griffin:2021znd].
At the same time, currently running experiments continue to publish their results and analyses, setting increasingly strict bounds on the DM parameter space.
In such a dynamic field, `obscura` can be an invaluable tool due to its high level of adaptability and facilitate and accelerate the development of new, reliable research software for the preparation of a DM discovery in the hopefully near future.


# Acknowledgements
The author thanks Radovan Bast for  valuable  discussions and support regarding research software engineering.
The author was supported by the Knut & Alice Wallenberg Foundation (PI, Jan Conrad).

# References# References

## Experimental data

### CRESST-II

> Angloher, G. *et al*, **Results on light dark matter particles with a low-threshold CRESST-II detector**, [Eur.Phys.J. C76 (2016) no.1, 25](https://doi.org/10.1140/epjc/s10052-016-3877-3) , [[arXiv:1509.01515]](https://arxiv.org/abs/1509.01515)

The files in */CRESST-II/* are the arXiv ancillary files of

> Angloher, G. *et al*, **Description of CRESST-II data**, [[arXiv:1701.08157]](https://arxiv.org/abs/1701.08157)

### CRESST-III

> Abdelhameed, A.H. *et al*, **First results from the CRESST-III low-mass dark matter program**, [Phys.Rev.D 100 (2019) 10, 102002](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.100.102002) , [[arXiv:1904.00498]](https://arxiv.org/abs/1904.00498)

The files in */CRESST-III/* are the arXiv ancillary files of

> Abdelhameed, A.H. *et al*, **Description of CRESST-III data**, [[arXiv:1905.07335]](https://arxiv.org/abs/1905.07335)

### CRESST-surface

The files in */CRESST-surface/* are the arXiv ancillary files of

> Angloher, G. *et al*, **Results on MeV-scale dark matter from a gram-scale cryogenic calorimeter operated above ground**, [Eur.Phys.J.C 77 (2017) 9, 637](https://link.springer.com/article/10.1140%2Fepjc%2Fs10052-017-5223-9) , [[arXiv:1707.06749]](https://arxiv.org/abs/1707.06749)

### XENON10-S2

> Trigger efficiencies from [arXiv:1206.2644](https://arxiv.org/abs/1206.2644)

### XENON100-S2

> Acceptance and trigger efficiencies from [arXiv:1605.06262](https://arxiv.org/abs/1605.06262)

### XENON1T-S2
> Efficiency from Fig. 2 of [arXiv:1907.11485](https://arxiv.org/abs/1907.11485) (see also App. E of [arXiv:1912.08204](https://arxiv.org/abs/1912.08204))

## Crystal form factors
The crystal form factors in the folder */Semiconductors/* were taken from [this website](http://ddldm.physics.sunysb.edu/ddlDM/) and were tabulated with [QEdark](https://github.com/adrian-soto/QEdark_repo). For details, we refer to the corresponding publication:

>Essig, R. *et al*, **Direct Detection of sub-GeV Dark Matter with Semiconductor Targets**, [JHEP 05 (2016) 046](https://link.springer.com/article/10.1007/JHEP05(2016)046) , [[arXiv:1509.01598]](https://arxiv.org/abs/1509.01598).

## Nuclear Data
Sources for spin expectation values
		
- F-19, Na-23, Al-27, Si-29,Ge-73, I-127, Xe-12, Xe-131	[[arXiv:1304.7684]](https://arxiv.org/abs/1304.7684)
- Remaining from [[arXiv:hep-ph/0406218]](https://arxiv.org/abs/hep-ph/0406218) and references therein.# Contributing to obscura
Your feedback and contribution to obscura are very welcome. This could be:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features
- Becoming a maintainer

## We Use [Github Flow](https://guides.github.com/introduction/flow/index.html)
We use github to host code, to track issues and feature requests, as well as accept pull requests. Pull requests are the best way to propose changes to the codebase (we use [Github Flow](https://guides.github.com/introduction/flow/index.html)). We actively welcome your pull requests:

1. Fork the repo and create your branch from `master`.
2. If you've added code that should be tested, add tests using [googletest](https://github.com/google/googletest).
3. If you've changed existing code, please check if the documentation requires an update.
4. Ensure all tests pass.
5. Keep the coding style consistent, ideally using `clang-format` (see below).
6. Issue that pull request!

## Any contributions you make will be under the MIT Software License
In short, when you submit code changes, your submissions are understood to be under the same [MIT License](http://choosealicense.com/licenses/mit/) that covers the project. Feel free to contact the maintainers if that's a concern.

## Report bugs using Github's [issues](https://github.com/temken/obscura/issues)
We use GitHub issues to track public bugs. Report a bug by opening an [issue](https://github.com/temken/obscura/issues).

## Write bug reports with detail, background, and sample code
Your bug report should contain:

- A quick summary and/or background
- Steps to reproduce
  - Be specific!
  - Give sample code if you can.
- What you expected would happen
- What actually happens
- Notes (possibly including why you think this might be happening, or stuff you tried that didn't work)


## Use a Consistent Coding Style
We want the coding style to remain consistent throughout the library. For this purpose, we use [clang-format](https://clang.llvm.org/docs/ClangFormat.html). With [this style configuration file](https://github.com/temken/clang-format-style/blob/master/.clang-format) you can find the specific format choices currently used. Most importantly, tabs are used for indentation.

## References
This document was adapted from this open-source contribution [template](https://gist.github.com/briandk/3d2e8b3ec8daf5a27a62).
