# Contributing Guidelines

Contributions to Riroriro are welcome. If you have an issue, a code contribution or a documentation contribution, thanks for helping to improve Riroriro!

## Issues

When creating an issue report, please be as specific as possible when describing how to reproduce an issue, and include both the intended/expected result and what you are actually getting.

## Code contributions

If you are performing a bugfix or adding a feature, please fork the repository, create a dedicated branch for your contribution and create a pull request when you are done.

## Documentation contributions

Riroriro has documentation written in reStructuredText, which can be found in the docs/ directory. If you have an improvement for the documentation, please put it on a dedicated branch as with code contributions. If you are adding a new feature in code contributions, you are encouraged to add documentation about it in the docs/ directory.# riroriro

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4588070.svg)](https://doi.org/10.5281/zenodo.4588070)

**Riroriro** is a set of Python modules containing functions to simulate the gravitational waveforms of mergers of black holes and/or neutron stars, and calculate several properties of these mergers and waveforms, specifically relating to their observability by gravitational wave detectors. Riroriro combines areas covered by previous gravitational wave models (such as gravitational wave simulation, SNR calculation, horizon distance calculation) into a single package with broader scope and versatility in Python, a programming language that is ubiquitous in astronomy. Aside from being a research tool, Riroriro is also designed to be easy to use and modify, and it can also be used as an educational tool for students learning about gravitational waves.

The modules “inspiralfuns”, “mergerfirstfuns”, “matchingfuns”, “mergersecondfuns” and “gwexporter”, in that order, can be used to simulate the strain amplitude and frequency of a merger gravitational waveform. The module “snrcalculatorfuns” can compare such a simulated waveform to a detector noise spectrum to calculate a signal-to-noise ratio (SNR) for that signal for that detector. The module “horizondistfuns” calculates the horizon distance of a merger given its waveform, and the module “detectabilityfuns” evaluates the detectability of a merger given its SNR.

Riroriro is installable via pip:

    pip install riroriro

More information on the pip installation can be found here: https://pypi.org/project/riroriro/

Tutorials for Riroriro can be found here: https://github.com/wvanzeist/riroriro_tutorials

Full documentation of each of the functions of Riroriro can be found here: https://wvanzeist.github.io/

Riroriro is one of several Python packages associated with **BPASS** (Binary Population And Spectral Synthesis), a suite of programs that simulates the evolution of a population of binary and single-star systems from a wide range of initial conditions. Each of these associated packages are named after native animals of New Zealand. The riroriro (*Gerygone igata*, also known as the grey warbler) is a small bird that can be recognised by its distinctive melodious call but is rarely seen, similarly to how black hole binary mergers are detected by their gravitational wave signals rather than visually.

The central website of BPASS, which also contains links to related programs, can be found here: https://bpass.auckland.ac.nz

## Paper

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02968/status.svg)](https://doi.org/10.21105/joss.02968)

A paper describing Riroriro has been published in the Journal of Open Source Software. If you use Riroriro in your work, please cite this paper! https://doi.org/10.21105/joss.02968

    @ARTICLE{2021JOSS....6.2968V,
           author = {{van Zeist}, Wouter G.~J. and {Stevance}, H{\'e}lo{\"i}se F. and {Eldridge}, J.~J.},
            title = "{Riroriro: Simulating gravitational waves and evaluating their detectability in Python}",
          journal = {The Journal of Open Source Software},
         keywords = {Python, neutron stars, astronomy, gravitational waves, black holes, General Relativity and Quantum Cosmology, Astrophysics - High Energy Astrophysical Phenomena},
             year = 2021,
            month = mar,
           volume = {6},
           number = {59},
              eid = {2968},
            pages = {2968},
              doi = {10.21105/joss.02968},
    archivePrefix = {arXiv},
           eprint = {2103.06943},
     primaryClass = {gr-qc},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2021JOSS....6.2968V},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
---
title: 'Riroriro: Simulating gravitational waves and evaluating their detectability in Python'
tags:
  - Python
  - astronomy
  - gravitational waves
  - black holes
  - neutron stars
authors:
  - name: Wouter G. J. van Zeist
    affiliation: 1
  - name: Héloïse F. Stevance
    affiliation: 1
  - name: J. J. Eldridge
    affiliation: 1
affiliations:
  - name: Department of Physics, University of Auckland, New Zealand
    index: 1
date: 2021
bibliography: paper.bib
---

# Summary

`Riroriro` is a Python package to simulate the gravitational waveforms of binary mergers of black holes and/or neutron stars, and calculate several properties of these mergers and waveforms, specifically relating to their observability by gravitational wave detectors.

The gravitational waveform simulation of `Riroriro` is based upon the methods of @buskirk2019, a paper which describes a computational implementation of an earlier theoretical gravitational waveform model by @huerta2017, using post-Newtonian expansions and an approximation called the implicit rotating source to simplify the Einstein field equations and simulate gravitational waves. `Riroriro`'s calculation of signal-to-noise ratios (SNR) of gravitational wave events is based on the methods of @barrett2018, with the simpler gravitational wave model `Findchirp` [@findchirp] being used for comparison and calibration in these calculations.

# Statement of Need

The field of gravitational wave astronomy has been especially active since the first observation of gravitational waves was announced in 2016 [@gw150914discovery]. Observations of gravitational waves from binary mergers can provide unique information about their progenitors and stellar populations, especially when combined with electromagnetic observations in the field called multi-messenger astronomy. A major factor in the successful detection and analysis of gravitational wave signals is the creation of simulations of such signals which observed data can be compared to. Because of this, multiple gravitational wave models have been created over the years. In particular, the gravitational wave observatories LIGO and Virgo have created their own models to use as templates in gravitational wave searches, with the main software for this being `LALSuite` [@lalsuite]. Various research groups have created waveform models, with some examples of recent sophisticated waveform models being `IMRPhenomXPHM` [@imrphenom] and `SEOBNRv4PHM` [@seobnr], which are both also available in `LALSuite`.

We have not tested if the waveform model of `Riroriro` that is based on @huerta2017 and @buskirk2019 is accurate enough to use for parameter estimation of detected gravitational wave transients as this was not within the scope of our project; it is likely that accurate parameter estimation requires careful modelling of the ringdown phase, especially for the most massive mergers. However, we use a level of accuracy adequate for our aim of modelling the detectability of gravitational wave transients predicted by stellar population syntheses. Furthermore, the code of `Riroriro` is structured and commented in such a way that each step in the process of the simulation is individually identifiable and modifiable by users. Users could also substitute in functions from other algorithms or even use the detectability modules on waveforms from other sources, as long as the user puts these in Riroriro’s format.

`Riroriro` combines areas covered by previous models (such as gravitational wave simulation, SNR calculation, horizon distance calculation) into a single package with broader scope and versatility in Python, a programming language that is ubiquitous in astronomy. Aside from being a research tool, `Riroriro` is also designed to be easy to use and modify, and it can also be used as an educational tool for students learning about gravitational waves.

# Features

Features of `Riroriro` include:

- Simulating the gravitational waveform signal from a binary merger of two black holes, two neutron stars or a black hole and a neutron star and outputting the data of this signal in terms of frequency and strain amplitude.
- Using a gravitational wave output and given a detector noise spectrum (such spectra are made publicly available by LIGO), calculating the signal-to-noise ratio (SNR) of the signal at a given distance assuming optimal alignment.
- Calculating the horizon distance (maximum distance at which an event could be observed) for a gravitational wave model and a given detector.
- Given the optimal-alignment SNR of an event, evaluating its detectability, the probability that the event would be detected with a SNR above the commonly used threshold of 8, if the alignment would be arbitrary. These results could then be combined with population synthesis calculations to estimate how many of the predicted mergers would be detected.

In addition, to help users get started with `Riroriro`, we have created Jupyter Notebook tutorials ([see here](https://github.com/wvanzeist/riroriro_tutorials)).

# Research

`Riroriro` has been used for research in conjunction with `BPASS`, a suite of computer programs that simulates the evolution of a population of binary and single-star systems from a wide range of initial conditions and predicts their electromagnetic spectral emission [@bpass1; @bpass2]. There is also a Python interface for `BPASS` called `Hoki` [@hoki]. This research took rates of formation of merging systems from `BPASS` and then evaluated the detectability of the gravitational wave signals from those systems using `Riroriro`. This was done to obtain predictions of the rates at which gravitational waves of different types would be expected to be observed, which can then be directly compared to those events found by the LIGO/Virgo gravitational wave observatories [@massdistribution].

# Acknowledgments

HFS and JJE acknowledge support from the University of Auckland and also the Royal Society of New Zealand Te Apārangi under the Marsden Fund​.

# References
