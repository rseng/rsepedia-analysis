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
=====================
3. The target classes
=====================

The basic hope of direct detection experiments is that the DM particles from the galactic halo occasionally collide with ordinary particles, see e.g. [Nobile2021]_.
The original target of direct DM searches were nuclear recoils.
Later on also electron targets gained more and more attention in the context of sub-GeV dark matter [Essig2012]_.

In *obscura* each target type is represented by a class, that can e.g. be passed to the cross section functions of the ``DM_Particle`` class, see :ref:`Section_DM_Particle`.

---------------
Nuclear targets
---------------

For nuclear recoil experiments, we define two target classes, ``Isotope`` and ``Nucleus`` that are declared in `/include/obscura/Target_Nucleus.hpp <https://github.com/temken/obscura/blob/master/include/obscura/Target_Nucleus.hpp>`_.

^^^^^^^^^^^^^^^^^^^^^
The ``Isotope`` class
^^^^^^^^^^^^^^^^^^^^^

A nuclear isotope is characterized by the number Z of protons and A of nucleons (protons *and* neutrons), its mass, spin, and average spin contribution for protons and neutrons as required e.g. in the context of spin-dependent nuclear interactions.

^^^^^^^^^^^^^^^^^^^^^
The ``Nucleus`` class
^^^^^^^^^^^^^^^^^^^^^

In addition to ``Isotope``, we also define a ``Nucleus`` class which mainly consists of a number of isotopes with given relative abundances.


"""""""""""""""""""""""""""""""
Construction of nuclear targets
"""""""""""""""""""""""""""""""

There are different ways to construct instances of ``Isotope`` and ``Nucleus``.

**Example:** Assume we are interested in oxygen as a target, either the isotope O-16 or the element of various isotopes.

.. code-block:: c++

    #include "obscura/Target_Nucleus.hpp"

    // ...

    // We can define O-16 via the constructor.
    Isotope oxygen_16(8,16);

This instance of an oxygen isotope however has no knowledge of e.g. its spin or relative abundance in nature.

For this purpose, *obscura* contains a nuclear data set, see `/data/Nuclear_Data.txt <https://github.com/temken/obscura/blob/master/data/Nuclear_Data.txt>`_ ([Bednyakov2005]_ [Klos2013]_), which can be accessed through the following function defined in *Target_Nucleus.hpp*.

.. code-block:: c++

   extern Isotope Get_Isotope(unsigned int Z, unsigned int A);
   extern Nucleus Get_Nucleus(unsigned int Z);
   extern Nucleus Get_Nucleus(std::string name);

Using these functions, we can construct isotopes and nuclei simply as

.. code-block:: c++

    #include "obscura/Target_Nucleus.hpp"

    // ...

    Isotope oxygen_16 = Get_Isotope(8,16);
    Nucleus oxygen = Get_Nucleus(8);
    Nucleus oxygen_alternative = Get_Nucleus("O");


The last two lines construct an instance of the ``Nucleus`` class containing all isotopes of oxygen including their relative abundance, spin, and average spin contribution of protons and neutrons.

-------------------------
Electron targets in atoms
-------------------------

For sub-GeV DM searches, an important target are electrons bound in atoms [Essig2012]_.
To take into account the fact that electrons are bound states, we need to evaluate the *ionization form factor* or *atomic response function* for each electronic orbital [Catena2019]_.

The target classes for atomic electrons are declared in `/include/obscura/Target_Atom.hpp <https://github.com/temken/obscura/blob/master/include/obscura/Target_Atom.hpp>`_.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The ``Atomic_Electron`` class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first target class in this context is ``Atomic_Electron``.

By constructing an instance of this class, the tabulated ionization form factor is imported from `/data/Form_Factors_Ionization/ <https://github.com/temken/obscura/tree/master/data/Form_Factors_Ionization>`_.


^^^^^^^^^^^^^^^^^^
The ``Atom`` class
^^^^^^^^^^^^^^^^^^

Having target classes for nuclei and bound electrons, we can combine them into a single atomic target, consisting of a nucleus and a number of bound electrons.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Included ionization form factors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
At this point, *obscura* comes with the ionization form factors of

* Xenon (5p, 5s, 4d, 4p, 4s)
* Argon (3p, 3s, 2p, 2s, 1s)

The tables can be found under `/data/Form_Factors_Ionization/ <https://github.com/temken/obscura/tree/master/data/Form_Factors_Ionization>`_.
They have been tabulated using the `DarkARC <https://github.com/temken/DarkARC>`_ code as described in detail in [Catena2019]_.

The easiest way to access the ionization form factors is by constructing an instance of ``Atom``, as seen in this example.

.. code-block:: c++

    #include "libphysica/Natural_Units.hpp"

    #include "obscura/Target_Atom.hpp"

    using namespace libphysica::natural_units;
    // ... 

    Atom xenon("Xe");
    Atom argon("Ar");

    // For example, to access the ionization form factor of xenon's 5s (quantum numbers n=5, l=0) orbital for a given momentum transfer q and energy E_e:
    int n = 5; int l = 0;
    double q = 0.5 * keV;
    double E_e = 10.0 * eV;

    std::cout << xenon.Electron(n, l).Ionization_Form_Factor(q, E_e) << std::endl;


----------------------------
Electron targets in crystals
----------------------------

One of the most important targets for sub-GeV DM detectors are crystals, such as e.g. semiconductors [Essig2016]_.
The electronic properties of the target material is encapsulated in the crystal form factor which is tabulated and can be found in `/data/Semiconductors/ <https://github.com/temken/obscura/tree/master/data/Form_Factors_Ionization>`_.
The included crystals are

* Silicon semiconductors
* Germanium semiconductors

The tables have been generated using `QEdark <http://ddldm.physics.sunysb.edu/ddlDM/>`_, a module of `Quantum ESPRESSO <https://www.quantum-espresso.org/>`_.

Also for crystals, *obscura* contains a target class ``Crystal`` declared in `/include/obscura/Target_Crystal.hpp <https://github.com/temken/obscura/blob/master/include/obscura/Target_Crystal.hpp>`_.

The crystal form factor, similarly to the ionization form factors, are imported by the class constructor. Here is an example of how to access the crystal form factor.

.. code-block:: c++

  #include "libphysica/Natural_Units.hpp"

  #include "obscura/Target_Crystal.hpp"

  using namespace libphysica::natural_units;

  // ...

  Crystal silicon("Si");
  Crystal germanium("Ge");

  double q = 0.5 * keV;
  double E_e = 10.0 * eV;

  std::cout << silicon.Crystal_Form_Factor(q, E_e) << std::endl;
  std::cout << germanium.Crystal_Form_Factor(q, E_e) << std::endl;==============================
6. The ``DM_Detector`` classes
==============================

The details of a direct detection experiment are summarized in the ``DM_Detector`` class declared in `/include/obscura/Direct_Detection.hpp <https://github.com/temken/obscura/blob/master/include/obscura/Direct_Detection.hpp>`_.
In particular, it responsible for:

1. The statistical methods to compute likelihoods and exclusion limits. Since these are independent of the type of experiment, this functionality is part of the base class ``DM_Detector``. As of now, *obscura* implements the following statistical analyses.
   1. Poisson statistics
   2. Binned Poisson statistics
   3. Maximum gap following [Yellin2002]_. 
2. The detector details, such as detection efficiencies, energy resolution, target particles, etc. These can be very specific and are implemented in classes derived from ``DM_Detector``, e.g. ``DM_Detector_Nucleus``.


We provide a number of examples of how to construct different instances of derived classes of ``DM_Detector``.

--------------------------
Nuclear recoil experiments
--------------------------

For experiments looking for DM induced nuclear recoils, *obscura* contains the ``DM_Detector_Nucleus`` class that is declared in `/include/obscura/Direct_Detection_Nucleus.hpp <https://github.com/temken/obscura/blob/master/include/obscura/Direct_Detection_Nucleus.hpp>`_.

For example, assume we have a nuclear recoil experiment with :math:`\mathrm{CaWO}_4` crystals, an energy threshold of 500 eV, and an exposure of 100 kg days.
This information suffices to define a toy experiment.

.. code-block:: c++

   #include "libphysica/Natural_Units.hpp"

   #include "obscura/Target_Nucleus.hpp"
   #include "obscura/DM_Detector_Nucleus.hpp"

   using namespace libphysica::natural_units;

   // ...

   double exposure = 100.0 * kg * day;
   std::vector<Nucleus> nuclear_targets = {obscura::Get_Nucleus(8), obscura::Get_Nucleus(20), obscura::Get_Nucleus(74)};
   std::vector<double> target_ratios = {4, 1, 1};
   double energy_threshold = 500 * eV;
   obscura::DM_Detector_Nucleus detector("Nuclear recoil experiment", exposure, nuclear_targets, target_ratios);

---------------------------
Electron recoil experiments
---------------------------

For electron recoil experiments with atomic targets, we have to use the ``DM_Detector_Ionization_ER`` class that can be found in `/include/obscura/Direct_Detection_ER.hpp <https://github.com/temken/obscura/blob/master/include/obscura/Direct_Detection_ER.hpp>`_

Here is an example of a xenon target experiment probing DM-electron interactions and DM induced ionizations. As exposure we choose 100 kg days, and we furthermore assume that only events with at least 4 ionized electrons can be detectred.

.. code-block:: c++

   #include "libphysica/Natural_Units.hpp"

   #include "obscura/DM_Detector_ER.hpp"

   using namespace libphysica::natural_units;

   // ...

   double exposure = 100.0 * kg * day;
   obscura::DM_Detector_Ionization_ER xenon_experiment("Electron recoil experiment", exposure, "Xe");
   argon_experiment.Use_Electron_Threshold(4);

Alternatively, many experiments looking for sub-GeV DM use semiconductor crystals as targets. In this case, there is another derived class, ``DM_Detector_Crystal``.

Again we construct an example toy experiment. This time, we assume a silicon crystal target, choose an exposure of 10 g year, and assume that only events with at least 2 electron-hole pairs can trigger the detector.

.. code-block:: c++

   #include "libphysica/Natural_Units.hpp"

   #include "obscura/DM_Detector_Crystal.hpp"

   using namespace libphysica::natural_units;

   // ...

   double exposure = 10.0 * gram * year;
   obscura::DM_Detector_Crystal silicon_experiment("Crystal target experiment", exposure, "Si");
   silicon_experiment.Use_Q_Threshold(2);

.. _Section_DM_Particle:

==============================
4. The ``DM_Particle`` classes
==============================

The ``DM_Particle`` class and its derived classes are responsible for the particle physics aspects of direct detection.
In particular, an instance of ``DM_Particle`` entails the particle properties of a DM candidate particle, such as its mass, spin, and its differential and total interaction cross sections with nuclei or electrons.
The base class's functions provide an interface that is sufficient for the calculation of e.g. event rates at direct DM search experiments.
Furthermore, it contains a number of functions regarding scattering angles, their distributions and sampling, which might not be relevant for direct detection, but can be used in the context of e.g. MC simulations.

--------------------------
The interface / base class
--------------------------

The abstract base class is defined in `/include/obscura/DM_Particle.hpp <https://github.com/temken/obscura/blob/master/include/obscura/DM_Particle.hpp>`_ and all the member functions and parameters can be seen there.

The most important (virtual) functions for direct detection specific calculations are the differential cross sections.

.. code-block:: c++

   //Differential cross sections for nuclear targets
   virtual double dSigma_dq2_Nucleus(double q, const Isotope& target, double vDM, double param = -1.0) const { return 0.0; };
   double dSigma_dER_Nucleus(double ER, const Isotope& target, double vDM, double param = -1.0) const;
   double d2Sigma_dER_dEe_Migdal(double ER, double Ee, double vDM, const Isotope& isotope, Atomic_Electron& shell) const;

   // Differential cross section for electron targets
   virtual double dSigma_dq2_Electron(double q, double vDM, double param = -1.0) const { return 0.0; };
   virtual double d2Sigma_dq2_dEe_Ionization(double q, double Ee, double vDM, Atomic_Electron& shell) const { return 0.0; };
   virtual double d2Sigma_dq2_dEe_Crystal(double q, double Ee, double vDM, Crystal& crystal) const { return 0.0; };

We point out that here we have to pass instances of the target classes discussed in the previous section (i.e. nuclear isotopes, atomic electrons, and electrons in crystals).
Also included is a simple implementation of Migdal scatterings with atomic targets based on [Essig2020]_.

The most standard DM candidate considered in the direct detection literature is a WIMP with SI or SD interactions. *obscura* contains derived classes for each of these scenarios, which are declared in `/include/obscura/DM_Particle_Standard.hpp <https://github.com/temken/obscura/blob/master/include/obscura/DM_Particle_Standard.hpp>`_.

----------------------------------
Spin-Independent (SI) interactions
----------------------------------

The differential cross section for SI nuclear interactions is given by

.. math::
	\frac{\mathrm{d}\sigma^{\rm SI}_{N}}{\mathrm{d} E_R} =\frac{m_N}{2\pi v_\chi^2}\left[f_p Z+f_n(A-Z) \right]^2 \left|F^{\rm SI}_N\left(E_R\right)\right|^2\, .

For details, we refere to e.g. chapter 3.4 of [Emken2019]_.

The class ``DM_Particle_SI`` is derived from ``DM_Particle`` and evaluates the cross sections of SI interactions with nuclei (and electrons).

The following example demonstrates how to

* construct an instance of ``DM_Particle_SI`` that describes a DM particle of 10 GeV mass.
* set the SI proton cross section to :math:`\sigma_p=10^{-40}\mathrm{cm}^2`, the electron cross section of :math:`\sigma_e=10^{-36}\mathrm{cm}^2`.
* to evaluate the differential and total scattering cross section with argon nuclei.

.. code-block:: c++

  #include "libphysica/Natural_Units.hpp"

  #include "obscura/DM_Particle_Standard.hpp"

  using namespace libphysica::natural_units;

  // ...

  // Declare the DM particle
  obscura::DM_Particle_SI dm(10.0 * GeV);
  dm.Set_Sigma_Proton(1.0e-40 * cm * cm);
  dm.Set_Sigma_Electron(1.0e-36 * cm * cm);

  // Define the target
  obscura::Isotope argon = obscura::Get_Isotope(18, 40);

  // Evaluate cross sections
  double E_R = 1.0 * keV;
  double v_DM = 300.0 * km/sec;
  double diff_cross_section = dm.dSigma_dER_Nucleus(E_R, argon, v_DM);
  double tot_cross_section = dm.Sigma_Total_Nucleus(argon, v_DM);

  // Convert to other units
  std::cout <<In_Units(diff_cross_section, cm * cm / keV)<<std::endl;
  std::cout <<In_Units(tot_cross_section, cm * cm)<<std::endl;

--------------------------------
Spin-Dependent (SD) interactions
--------------------------------

The differential cross section for SD nuclear interactions is given by

.. math::
	\frac{\mathrm{d} \sigma_N^{\rm SD}}{\mathrm{d} E_R} = \frac{2m_N}{\pi v_\chi^2}\frac{J+1}{J}\left(f_p \langle S_p\rangle +f_n \langle S_N\rangle\right)^2 \left.F_N^{\rm SD}(E_R)^2\right|

Similarly to ``DM_Particle_SI``, we also define a ``DM_Particle_SD`` class, which evaluates this cross section for nuclear targets with spin :math:`S\neq0`.===============
Release history
===============

* 10.11.2021: Release of `version 1.0.0 <https://github.com/temken/obscura/releases/tag/v1.0.0>`_
* 23.02.2021: Release of `version 0.1.0 <https://github.com/temken/obscura/releases/tag/v0.1.0>`_==================
References
==================

.. .. [ref] author, *title*, `journal <>`_, `[arXiv:xxxx] <https://arxiv.org/abs/xxxx>`_.
.. [Catena2019] R. Catena et al., *Atomic responses to general dark matter-electron interactions*, `Phys.Rev.Res. 2 (2020) 3, 033195 <https://doi.org/10.1103/PhysRevResearch.2.033195>`_, `[arXiv:1912.08204] <https://arxiv.org/abs/1912.08204>`_.
.. [Bednyakov2005] V.A. Bednyakov, *Nuclear spin structure in dark matter search: The Zero momentum transfer limit*, Phys.Part.Nucl. 36 (2005) 131-152, `[arXiv:0406218] <https://arxiv.org/abs/0406218>`_.
.. [Emken2019] T. Emken, *Dark Matter in the Earth and the Sun - Simulating Underground Scatterings for the Direct Detection of Low-Mass Dark Matter*, PhD thesis 2019, `[arXiv:1906.07541] <https://arxiv.org/abs/1906.07541>`_.
.. [Essig2012] R. Essig et al. , *Direct Detection of Sub-GeV Dark Matter*, `Phys.Rev.D 85 (2012) 076007 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.85.076007>`_ , `[arXiv:1108.5383] <https://arxiv.org/abs/1108.5383>`_.
.. [Essig2016] R. Essig et al. , *Direct Detection of sub-GeV Dark Matter with Semiconductor Targets*, `JHEP 05 (2016) 046 <https://doi.org/10.1007/JHEP05(2016)046>`_ , `[arXiv:1509.01598] <https://arxiv.org/abs/1509.01598>`_.
.. [Essig2020] R. Essig et al. , *Relation between the Migdal Effect and Dark Matter-Electron Scattering in Isolated Atoms and Semiconductors*, `Phys.Rev.Lett. 124 (2020) 2, 021801 <https://doi.org/10.1103/PhysRevLett.124.021801>`_ , `[arXiv:1908.10881] <https://arxiv.org/abs/1908.10881>`_.
.. [Evans2019] N.W. Evans et al., *Refinement of the standard halo model for dark matter searches in light of the Gaia Sausage*, `Phys.Rev.D 99 (2019) 2, 023012 <https://doi.org/10.1103/PhysRevD.99.023012>`_, `[arXiv:1810.11468] <https://arxiv.org/abs/1810.11468>`_.
.. [Klos2013] P. Klos et al., *Large-scale nuclear structure calculations for spin-dependent WIMP scattering with chiral effective field theory currents*, `Phys.Rev.D 88 (2013) 8, 083516 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.88.083516>`_, `[arXiv:1304.7684] <https://arxiv.org/abs/1304.7684>`_.
.. [Nobile2021] E. Del Nobile, *Appendiciario -- A hands-on manual on the theory of direct Dark Matter detection*, `[arXiv:2104.12785] <https://arxiv.org/abs/2104.12785>`_.
.. [Yellin2002] S. Yellin, *Finding an upper limit in the presence of unknown background*, `Phys.Rev.D 66 (2002) 032005 <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.66.032005>`_, `[arXiv:0203002] <https://arxiv.org/abs/0203002>`_.
=================
Contact & Support
=================

The author of *obscura* is `Timon Emken <https://timonemken.com/>`_.

For questions, support, bug reports, or other suggestions, please contact timon.emken@fysik.su.se or open an issue on `github <https://github.com/temken/obscura/issues>`_.==================
1. Getting started
==================

------------
Installation
------------

Before building *obscura*, there are a few libraries that need to be installed.

^^^^^^^^^^^^
Dependencies
^^^^^^^^^^^^

""""""""""""""""""""""""""""""""""""
1. `boost <https://www.boost.org/>`_
""""""""""""""""""""""""""""""""""""

To install *boost* on a Mac, we can use `homebrew <https://brew.sh/>`_ ::

	brew install boost

On Linux machines, run::

   sudo apt-get update && sudo apt-get install -yq libboost-all-dev


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
2. `libconfig <https://hyperrealm.github.io/libconfig/>`_
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To install *boost* on a Mac, we can use `homebrew <https://brew.sh/>`_ ::

	brew install libconfig

On Linux machines, you can build `libconfig` via::

	wget https://hyperrealm.github.io/libconfig/dist/libconfig-1.7.2.tar.gz
	tar -xvzf libconfig-1.7.2.tar.gz
	pushd libconfig-1.7.2
	./configure
	make
	sudo make install
	popd

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
3. `libphysica <https://github.com/temken/libphysica>`_
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

*libphysica* does not need to be installed. It will be downloaded and compiled during the CMake build.

^^^^^^^^^^^^^^^^
Download & Build
^^^^^^^^^^^^^^^^

The `obscura` source code can be downloaded by cloning this `git repository <https://github.com/temken/obscura>`_: ::

   git clone https://github.com/temken/obscura.git
   cd obscura

The code is compiled and the executable and library is built by `CMake <https://cmake.org/>`_. To build run the following commands from the repository's root folder.::

	cmake -E make_directory build
	cd build
	cmake -DCMAKE_BUILD_TYPE=Release -DCODE_COVERAGE=OFF ..
	cmake --build . --config Release
	cmake --install .

If everything worked well, the executable and library file are created as::

	bin/obscura
	lib/libobscura.a


-------------------------
Using *obscura* as a tool
-------------------------

*Obscura* can be used as a tool and builds an executable which can be run from */bin/* via::

./obscura config.cfg

As can be seen in the `/src/main.cpp <https://github.com/temken/obscura/blob/master/src/main.cpp>`_ file, this script computes direct detection limits and saves them in the */results/* folder.
The specifications of the exclusion limits (DM physics and halo model, statistics, experiment, mass range,...) are defined in a configuration file, in this case *config.cfg*.
For the handling of configuration files, *obscura* relies on `libconfig <https://hyperrealm.github.io/libconfig/>`_. 

^^^^^^^^^^^^^^^^^^^^^^
The configuration file
^^^^^^^^^^^^^^^^^^^^^^

The configuration file contains all input parameters necessary to define the various *obscura* models.

.. warning::

	The import of these parameters via libconfig is very case-sensitive. A float parameter has to be set to e.g. *1.0*, and **not** just *1*.

.. raw:: html

	<details>
	<summary><a>The full configuration file</a></summary>
 
.. code-block:: c++

   //obscura - Configuration File

   //ID
   	ID		=	"test";

   //Dark matter particle
   	DM_mass		  	=	0.1;		// in GeV
   	DM_spin		  	=	0.5;
   	DM_fraction		=	1.0;		// the DM particle's fractional abundance (set to 1.0 for 100%)
   	DM_light		=	false;		// Options: true or false. low mass mode

   	DM_interaction		=	"SI";		// Options: "SI" or "SD"

   	DM_isospin_conserved		=	true; 		// only relevant for SI and SD
   	DM_relative_couplings		=	(1.0, 0.0); //relation between proton (left) and neutron (right) couplings.
   												//only relevant if 'DM_isospin_conserved' is false.
   	DM_cross_section_nucleon	=	1.0e-36;	//in cm^2
   	DM_cross_section_electron	=	1.0e-36;	//in cm^2 (only relevant for SI and SD)
   	DM_form_factor		=	"Contact";	// Options: "Contact", "Electric-Dipole", "Long-Range", "General"
   												//(only relevant for SI)
   	DM_mediator_mass	=	0.0;		// in MeV (only relevant if 'DM_form_factor' is "General")

   //Dark matter distribution
   	DM_distribution 	=	"SHM";		//Options: "SHM", "SHM++", "File"
   	DM_local_density	=	0.4;		//in GeV / cm^3
   	
   	//Options for "SHM" and "SHM++"
   		SHM_v0		=	220.0;				//in km/sec
   		SHM_vObserver	=	(0.0, 232.0, 0.0);	//in km/sec
   		SHM_vEscape	=	544.0;				//in km/sec
   	//Options for "SHM++"
   		SHMpp_eta	=	0.2;
   		SHMpp_beta	=	0.9;
   	//Options for "File" (The file has to be a 2-column table of format v[km/sec] :: f(v) [sec/km])
   		file_path  = "DM_Speed_PDF.txt";

   //Dark matter detection experiment
   	DD_experiment	=	"Electron recoil";	//Options for nuclear recoils: "Nuclear recoil", "DAMIC_N_2011", "XENON1T_N_2017", "CRESST-II","CRESST-III", "CRESST-surface"
							//Options for electron recoils: "Semiconductor","protoSENSEI@MINOS","protoSENSEI@surface", "SENSEI@MINOS", "CDMS-HVeV_2018", "CDMS-HVeV_2020", "Electron recoil", "XENON10_S2", "XENON100_S2", "XENON1T_S2", "DarkSide-50_S2"

   	//Options for user-defined experiments ("Nuclear recoil", "Electron recoil", and "Semiconductor")
	  //General
	  DD_exposure 		=	1.0;	//in kg years
	  DD_efficiency 		=	1.0;	//flat efficiency
	  DD_observed_events 	=	0;		//observed signal events
	  DD_expected_background 	=	0.0;	//expected background events

	  //Specific options for "Nuclear recoil"
	  DD_targets_nuclear	=	(
	  				(4.0, 8),
	  				(1.0, 20),
	  				(1.0, 74)
	  			);				// Nuclear targets defined by atom ratio/abundances and Z
	  DD_threshold_nuclear    =	4.0;    //in keV
	  DD_Emax_nuclear         =	40.0;	//in keV
	  DD_energy_resolution    =	0.0;    //in keV

	  //Specific options for "Electron recoil" and "Semiconductor:
	  DD_target_electron	=	"Xe";	//Options for "Electron recoil": 	"Xe", "Ar"
	  								//Options for "Semiconductor":	"Si", "Ge"
	  DD_threshold_electron	=	4;		//In number of electrons or electron hole pairs.

   //Computation of exclusion limits
   	constraints_certainty	=	0.95;	//Certainty level
   	constraints_mass_min	=	0.02;	//in GeV										
   	constraints_mass_max	=	1.0;	//in GeV
   	constraints_masses	=	10;										
 
.. raw:: html

	</details>

----------------------------
Using *obscura* as a library
----------------------------

If we want to use *obscura* functions in an external code, we can do so and import it as a library.
We recommend to do this inside your CMake build, where *obscura* can be downloaded, built, included, and linked automatically during the build of your code.


As an instructional example `this repository <https://github.com/temken/template_cpp_cmake_obscura>`_ contains a C++ project template built with CMake that imports and uses the *obscura* library.====================================
7. Examples: Putting it all together
====================================


-------------------------------------------------------------
Computing the recoil spectrum of SI & SD nuclear interactions
-------------------------------------------------------------

As an example for nuclear recoils that also illustrates nicely the modular structure of *obscura*, we compute the nuclear recoil spectrum :math:`\frac{ \mathrm{d}R}{\mathrm{d}E_R}` for a 10 GeV DM particle interacting with xenon nuclei via spin-independent and spin-dependent interactions.

For the definition and details of the nuclear recoil spectrum, see e.g. chapter 3.5 of [Emken2019]_.

1. First we define the DM particle objects that describe SI and SD interactions

.. code-block:: c++

   // 1. DM particle (SI and Sd)
   obscura::DM_Particle_SI dm_SI(10.0 * GeV);
   dm_SI.Set_Sigma_Proton(1.0e-40 * cm * cm);
   dm_SI.Print_Summary();

   obscura::DM_Particle_SD dm_SD(10.0 * GeV);
   dm_SD.Set_Sigma_Proton(1.0e-40 * cm * cm);
   dm_SD.Print_Summary();    

The ``Print_Summary()`` function is a member of many of the classes and provides a terminal output that summarizes the object.

2. For the DM distribution we use the standard halo model with default parameters.

.. code-block:: c++

   // 2. DM distribution
   obscura::Standard_Halo_Model shm;
   shm.Print_Summary();

3. As target nuclei, we choose xenon and import the nuclear data.

.. code-block:: c++

   // 3. Direct detection targets
   obscura::Nucleus xenon = obscura::Get_Nucleus("Xe");
   xenon.Print_Summary();

4. With these three objects, we can compute the differential nuclear recoil spectrum for a given recoil energy :math:`E_R`.

.. code-block:: c++

   double E_R		= 1.0 * keV;
   double dRdER_SI = obscura::dRdER_Nucleus(E_R, dm_SI, shm, xenon);
   double dRdER_SD = obscura::dRdER_Nucleus(E_R, dm_SD, shm, xenon);

5. The results are given in natural units in powers of GeV. To convert it to another unit, we can use the unit functionality of the *libphysica* library.

.. code-block:: c++

   std::cout << "SI-interactions: \tdR/dER (1 keV) = " << In_Units(dRdER_SI, 1.0 / kg / year / keV) << " events / kg / year / keV" << std::endl;
   std::cout << "SD-interactions: \tdR/dER (1 keV) = " << In_Units(dRdER_SD, 1.0 / kg / year / keV) << " events / kg / year / keV" << std::endl;

.. raw:: html

	<details>
	<summary><a>The full main.cpp</a></summary>
 
.. code-block:: c++

   #include <iostream>

   #include "libphysica/Natural_Units.hpp"

   #include "obscura/DM_Halo_Models.hpp"
   #include "obscura/DM_Particle_Standard.hpp"
   #include "obscura/Direct_Detection_Nucleus.hpp"
   #include "obscura/Target_Nucleus.hpp"

   using namespace libphysica::natural_units;

   int main()
   {

   	// 1. DM particle (SI and Sd)
   	obscura::DM_Particle_SI dm_SI(10.0 * GeV);
   	dm_SI.Set_Sigma_Proton(1.0e-40 * cm * cm);
   	dm_SI.Print_Summary();

   	obscura::DM_Particle_SD dm_SD(10.0 * GeV);
   	dm_SD.Set_Sigma_Proton(1.0e-40 * cm * cm);
   	dm_SD.Print_Summary();

   	// 2. DM distribution
   	obscura::Standard_Halo_Model shm;
   	shm.Print_Summary();

   	// 3. Direct detection targets
   	obscura::Nucleus xenon = obscura::Get_Nucleus("Xe");
   	xenon.Print_Summary();

   	// 4. Evalute the nuclear recoil spectrum
   	double E_R		= 1.0 * keV;
   	double dRdER_SI = obscura::dRdER_Nucleus(E_R, dm_SI, shm, xenon);
   	double dRdER_SD = obscura::dRdER_Nucleus(E_R, dm_SD, shm, xenon);

   	std::cout << "SI-interactions: \tdR/dER (1 keV) = " << In_Units(dRdER_SI, 1.0 / kg / year / keV) << " events / kg / year / keV" << std::endl;
   	std::cout << "SD-interactions: \tdR/dER (1 keV) = " << In_Units(dRdER_SD, 1.0 / kg / year / keV) << " events / kg / year / keV" << std::endl;

   	return 0;
   }

.. raw:: html

	</details>

.. raw:: html

	<details>
	<summary><a>The terminal output</a></summary>

.. code-block::

   ----------------------------------------
   DM particle summary:
           Mass:                   10 GeV
           Spin:                   0.5
           Low mass:               [ ]

           Interaction:            Spin-Independent (SI)

           Coupling ratio fixed:   [x]
           Isospin conservation:   [x]
           Coupling ratio:         fn/fp = 1

           Sigma_P[cm^2]:          1e-40
           Sigma_N[cm^2]:          1e-40
           Sigma_E[cm^2]:          1e-40

           Interaction type:       Contact
   ----------------------------------------

   ----------------------------------------
   DM particle summary:
           Mass:                   10 GeV
           Spin:                   0.5
           Low mass:               [ ]
   Interaction:            Spin-Dependent (SD)

           Coupling ratio fixed:   [x]
           Isospin conservation:   [x]
           Coupling ratio:         fn/fp = 1

           Sigma_P[cm^2]:          1e-40
           Sigma_N[cm^2]:          1e-40
           Sigma_E[cm^2]:          1e-40

   ----------------------------------------

   Dark matter distribution - Summary
           Standard halo model (SHM)

           Local DM density[GeV/cm^3]:     0.4
           Speed domain [km/sec]:          [0,777]
           Average DM velocity [km/sec]:   (-11.1 , -232 , -7.3)
           Average DM speed [km/sec]:      330

           Speed dispersion v_0[km/sec]:   220
           Gal. escape velocity [km/sec]:  544
           Observer's velocity [km/sec]:   (11.1 , 232 , 7.3)
           Observer's speed [km/sec]:      233


   Xe
   Isotope Z       A       Abund.[%]       Spin    <sp>    <sn>
   ------------------------------------------------------------
   Xe-124  54      124     0.095           0       0       0
   Xe-126  54      126     0.089           0       0       0
   Xe-128  54      128     1.91            0       0       0
   Xe-129  54      129     26.4            0.5     0.01    0.329
   Xe-130  54      130     4.07            0       0       0
   Xe-131  54      131     21.2            1.5     -0.009  -0.272
   Xe-132  54      132     26.9            0       0       0
   Xe-134  54      134     10.4            0       0       0
   Xe-136  54      136     8.86            0       0       0
   Total:          131     99.999

   SI-interactions:        dR/dER (1 keV) = 13621.8 events / kg / year / keV
   SD-interactions:        dR/dER (1 keV) = 0.132525 events / kg / year / keV

.. raw:: html

   </details>

--------------------------------------------------------------------------
Exclusion limits for a sub-GeV DM particle via electron recoil experiments
--------------------------------------------------------------------------

As a second example for an application of *obscura*, we will compute the 95% confidence level exclusion limit on the DM-electron cross section for a sub-GeV DM particle.

We assume a DM mass of 100 MeV, and two different direct detection experiments.

1. An argon based experiment with an exposure of 100 kg years and an observational threshold of at least 4 ionized electrons.
2. A semiconductor experiment with Si crystal targets, an exposure of 10 gram years, and an observational threshold of minimum 2 electron-hole pairs.

Let us set up the different objects to obtain the limits.

1. First we define the DM particle object with 100 MeV mass.

.. code-block:: c++

   // 1. DM particle
   obscura::DM_Particle_SI dm(100.0 * MeV);
   dm.Print_Summary();   

2. For the DM distribution we again use the standard halo model with default parameters.

.. code-block:: c++

   // 2. DM distribution
   obscura::Standard_Halo_Model shm;
   shm.Print_Summary();

3. For the first experiment, we create an instance of the ``DM_Detector_Ionization_ER`` class and specify the desired detector properties of the toy experiment.

.. code-block:: c++

   // 3. Argon target experiment
   obscura::DM_Detector_Ionization_ER argon_experiment("Argon toy experiment", 100.0 * kg * year, "Ar");
   argon_experiment.Use_Electron_Threshold(4);
   argon_experiment.Print_Summary();

4. The same for the semiconductor experiment:

.. code-block:: c++

   // 4. Si target experiment
   obscura::DM_Detector_Crystal silicon_experiment("Silicon toy experiment", 10.0 * gram * year, "Si");
   silicon_experiment.Use_Q_Threshold(2);
   silicon_experiment.Print_Summary();

4. With these three objects, we can compute the limit on the DM-electron cross section.

.. code-block:: c++

   // 5. Compute the 95% CL exclusion limits for m = 100.0 MeV
   double limit_Ar = argon_experiment.Upper_Limit(dm, shm, 0.95);
   double limit_Si = silicon_experiment.Upper_Limit(dm, shm, 0.95);

5. As in the previous example, the results are given in natural units in powers of GeV. We convert it to :math:`\mathrm{cm}^2`, and print the result on the terminal.

.. code-block:: c++

   std::cout << "Argon experiment: \tsigma_e < " << In_Units(limit_Ar, cm * cm) << " cm^2 (95%CL)" << std::endl;
   std::cout << "Silicon experiment: \tsigma_e < " << In_Units(limit_Si, cm * cm) << " cm^2 (95%CL)" << std::endl;

.. raw:: html

	<details>
	<summary><a>The full main.cpp</a></summary>
 
.. code-block:: c++

   #include <iostream>

   #include "libphysica/Natural_Units.hpp"

   #include "obscura/DM_Halo_Models.hpp"
   #include "obscura/DM_Particle_Standard.hpp"
   #include "obscura/Direct_Detection_Crystal.hpp"
   #include "obscura/Direct_Detection_ER.hpp"
   #include "obscura/Target_Atom.hpp"
   #include "obscura/Target_Crystal.hpp"

   using namespace libphysica::natural_units;

   int main()
   {

   	// 1. DM particle
   	obscura::DM_Particle_SI dm(100.0 * MeV);
   	dm.Print_Summary();

   	// 2. DM distribution
   	obscura::Standard_Halo_Model shm;
   	shm.Print_Summary();

   	// 3. Argon target experiment
   	obscura::DM_Detector_Ionization_ER argon_experiment("Argon toy experiment", 100.0 * kg * year, "Ar");
   	argon_experiment.Use_Electron_Threshold(4);
   	argon_experiment.Print_Summary();

   	// 4. Si target experiment
   	obscura::DM_Detector_Crystal silicon_experiment("Silicon toy experiment", 10.0 * gram * year, "Si");
   	silicon_experiment.Use_Q_Threshold(2);
   	silicon_experiment.Print_Summary();

   	// 5. Compute the 95% CL exclusion limits for m = 100.0 MeV
   	double limit_Ar = argon_experiment.Upper_Limit(dm, shm, 0.95);
   	double limit_Si = silicon_experiment.Upper_Limit(dm, shm, 0.95);

   	std::cout << "Argon experiment: \tsigma_e < " << In_Units(limit_Ar, cm * cm) << " cm^2 (95%CL)" << std::endl;
   	std::cout << "Silicon experiment: \tsigma_e < " << In_Units(limit_Si, cm * cm) << " cm^2 (95%CL)" << std::endl;

   	return 0;
   }

.. raw:: html

	</details>

.. raw:: html

	<details>
	<summary><a>The terminal output</a></summary>
 
.. code-block::
  
   ----------------------------------------
   DM particle summary:
           Mass:                   100 MeV
           Spin:                   0.5
           Low mass:               [ ]

           Interaction:            Spin-Independent (SI)

           Coupling ratio fixed:   [x]
           Isospin conservation:   [x]
           Coupling ratio:         fn/fp = 1

           Sigma_P[cm^2]:          1e-40
           Sigma_N[cm^2]:          1e-40
           Sigma_E[cm^2]:          1e-40

           Interaction type:       Contact
   ----------------------------------------
   Dark matter distribution - Summary
           Standard halo model (SHM)

           Local DM density[GeV/cm^3]:     0.4
           Speed domain [km/sec]:          [0,777]
           Average DM velocity [km/sec]:   (-11.1 , -232 , -7.3)
           Average DM speed [km/sec]:      330

           Speed dispersion v_0[km/sec]:   220
           Gal. escape velocity [km/sec]:  544
           Observer's velocity [km/sec]:   (11.1 , 232 , 7.3)
           Observer's speed [km/sec]:      233


   ----------------------------------------
   Experiment summary:     Argon toy experiment
           Target particles:       Electrons
           Exposure [kg year]:     100
           Flat efficiency [%]:    100
           Observed events:        0
           Expected background:    0
           Statistical analysis:   Poisson


           Electron recoil experiment (ionization).
           Target(s):
                           Ar      (100%)
           Electron bins:          [ ]
           PE (S2) bins:           [ ]
                   Ne threshold:   4
                   Ne max:         15
   ----------------------------------------


   ----------------------------------------
   Experiment summary:     Silicon toy experiment
           Target particles:       Electrons
           Exposure [kg year]:     0.01
           Flat efficiency [%]:    100
           Observed events:        0
           Expected background:    0
           Statistical analysis:   Poisson


           Electron recoil experiment (semiconductor).
           Target:                 Si semiconductor
           eh pair threshold:      2
   ----------------------------------------

   Argon experiment:       sigma_e < 1.67038e-41 cm^2 (95%CL)
   Silicon experiment:     sigma_e < 1.1756e-39 cm^2 (95%CL)

.. raw:: html

	</details>.. image:: https://github.com/temken/obscura/actions/workflows/main.yml/badge.svg?branch=master
   :target: https://github.com/temken/obscura/actions/workflows/main.yml
   :alt: Build Status
.. image:: https://codecov.io/gh/temken/obscura/branch/master/graph/badge.svg?token=1Pe1QMcngr
   :target: https://codecov.io/gh/temken/obscura
   :alt: Code Coverage 
.. image:: https://readthedocs.org/projects/obscura/badge/?version=latest
   :target: https://obscura.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License

========================================================================================
*obscura* - Direct detection of dark matter with nucleus and electron recoil experiments
========================================================================================

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5665890.svg
   :target: https://doi.org/10.5281/zenodo.5665890
   :alt: DOI
.. image:: https://joss.theoj.org/papers/fd8076268036956d3bf08193c4fc2db9/status.svg
   :target: https://joss.theoj.org/papers/fd8076268036956d3bf08193c4fc2db9
   :alt: JOSS paper

A modular C++ tool and library for dark matter direct detection computations for both nuclear and electron recoil experiments.

The purpose of this documentation or manual is to provide insight into the polymorphic class structure of *obscura* and how it can be applied in different contexts.
It should also serve as a guide and describe the usage of *obscura* via code examples.

The documentation does not contain a review of the physics implemented in the library.
For more physics details, we refer to e.g. chapter 3 of [Emken2019]_ or [Nobile2021]_.

If you want to contribute to `obscura`, please check out the `contribution guidelines <https://github.com/temken/obscura/blob/master/docs/CONTRIBUTING.md>`_.

.. image:: https://raw.githubusercontent.com/temken/obscura/master/paper/FlowChart.png
   :width: 500

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   01_Getting_Started
   02_Main_Modules
   03_Targets
   04_DM_Particle
   05_DM_Distribution
   06_DM_Detector
   07_Examples
   08_Experiments
   09_Citations
   10_Release_History
   11_License
   12_Contact
   References

For the interpretation of past and future direct searches for DM particles, it is important to be able to provide accurate predictions for event rates and spectra under a variety of possible and viable assumptions in a computationally efficient way.
While there exists a few tools to compute DM induced nuclear recoil spectra, such as `DDCalc <https://ddcalc.hepforge.org/>`_ or `WimPyDD <https://wimpydd.hepforge.org/>`_, `obscura` is not limited to nuclear targets.
Instead its main focus lies on sub-GeV DM searches probing electron recoils which typically requires methods from atomic and condensed matter physics, see e.g. [Essig2012]_ or [Catena2019]_.
In the context of sub-GeV DM searches, new ideas such as target materials or detection techniques are being proposed regularly, and the theoretical modelling of these are getting improved continuosly.
At the same time, currently running experiments continue to publish their results and analyses, setting increasingly strict bounds on the DM parameter space.
In such a dynamic field, `obscura` can be an invaluable tool due to its high level of adaptability and facilitate and accelerate the development of new, reliable research software for the preparation of a DM discovery in the hopefully near future.================
Citing *obscura*
================

-----------
How to cite
-----------

If you decide to use this code, or if you want to add a reference to it, please cite the latest archived version,

    Emken, T., 2021, obscura - A C++ library for dark matter detection computations [Code, v1.0.0] [DOI:10.5281/zenodo.5665890]

.. raw:: html

	<details>
	<summary><a>Bibtex entry.</a></summary>
 
.. code-block::

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

.. raw:: html

	</details>


.. as well as the corresponding publication.::

..     Emken, T., 2021, obscura - A C++ library for dark matter detection computations


.. .. raw:: html

.. 	<details>
.. 	<summary><a>Bibtex entry.</a></summary>
 
.. .. code-block::


.. .. raw:: html

.. 	</details>

----------------------------------------------
Research and research software using *obscura*
----------------------------------------------

The library *obscura* has been applied to obtain the scientific results of the following papers

#. **Solar reflection of light dark matter with heavy mediators**
  
  Timon Emken

  .. image:: https://img.shields.io/badge/arXiv-2102.12483-B31B1B.svg
      :target: https://arxiv.org/abs/2102.12483
      :alt: [arXiv:2102.12483]


Here is a list of research software using *obscura*:

#. Emken, T., 2021, `Dark Matter Simulation Code for Underground Scatterings - Sun Edition (DaMaSCUS-SUN) <https://github.com/temken/DaMaSCUS-SUN>`_ Astrophysics Source Code Library, record `[ascl:2102.018] <https://ascl.net/2102.018>`_, `[DOI:10.5281/zenodo.4559874] <https://zenodo.org/record/4559874>`_

========================================
2. The modular structure of *obscura*
========================================

The computation of e.g. the electron recoil spectrum probed in direct detection experiments combines inputs from various fields of physics.
We need to specify the assumed *particle physics* of the DM particle.
The properties of the DM halo of the Milky way is an important *astrophysics* input.
For the description of the target particles, and how they react to a kick from an incoming DM particle, we need to include knowledge of *atomic*, *nuclear*, and *condensed matter physics*.
In order to make predictions, we furthermore need to define the *detection* experiments specifications.
Finally, the result of such an experiment needs to be interpreted using *statistics*.

.. image:: https://raw.githubusercontent.com/temken/obscura/master/paper/FlowChart.png
   :width: 500

This high level of modularity in this type of calculation needs to be reflected in the code's polymorphic structure.
The goal of *obscura* is to provide for each of the different inputs one generic interface or abstract base class, that comprises the general required functionalities, without specifying the detailed implementations further.
These depend on a multitude of assumptions which can change in different projects, for different users, etc.

If the base classes are defined properly, it is also possible and straight-forward to 

#. extend *obscura* by implementing further derived classes overriding the virtual functions of the base class.
#. design research software that is agnostic to the detailed implementation and thereby very generally applicable to a variety of scenarios. As long as our scientific functions are formulated in terms of these base functions, they will be able to handle any new implementation that comes in the form of derived classes.

The three most important abstract base classes of *obscura* are

#. ``DM_Particle``
#. ``DM_Distribution``
#. ``DM_Detector``

We will discuss the interface each of these classes provide in more detail.
But first we take a look at the detection targets in direct DM search experiments, namely nuclei, bound electrons in atoms, and bound electrons in crystals.=================================
8. Included experimental analyses
=================================

.. image:: https://raw.githubusercontent.com/temken/obscura/master/paper/obscura_DD_Constraints.png
    :width: 500

The module `Experiments.hpp <https://github.com/temken/obscura/blob/master/include/obscura/Experiments.hpp>`_ contains a series of functions that build a number of experimental analysis as instances of the ``DM_Detector`` class and its derivatives.


For example, for an analysis based on the CRESST-II experiment, we can construct the class instance via

.. code-block:: c++

    #include "obscura/Experiments.hpp"

    // ...
    
    DM_Detector_Nucleus detector = CRESST_II();

    // ...


The following nuclear and electron recoil direct detection experiments are implemented in *obscura*.

--------------------------
Nuclear recoil experiments
--------------------------

CRESST-II
^^^^^^^^^

* **Results on light dark matter particles with a low-threshold CRESST-II detector**

  CRESST Collaboration (G. Angloher et al.)

  .. image:: https://img.shields.io/badge/Eur.Phys.J.-C76(2016)no.1,25-255773.svg
      :target: https://link.springer.com/article/10.1140/epjc/s10052-016-3877-3
      :alt: Eur.Phys.J. C76 (2016) no.1, 25
  .. image:: https://img.shields.io/badge/arXiv-1509.01515-B31B1B.svg
      :target: https://arxiv.org/abs/1509.01515
      :alt: [arXiv:1509.01515]

* **Description of CRESST-II data**

  CRESST Collaboration (G. Angloher et al.)

  .. image:: https://img.shields.io/badge/arXiv-1701.08157-B31B1B.svg
      :target: https://arxiv.org/abs/1701.08157
      :alt: [arXiv:1701.08157]

CRESST-III
^^^^^^^^^^

* **First results on low-mass dark matter from the CRESST-III experiment**
  
  CRESST Collaboration (F. Petricca et al.)

  .. image:: https://img.shields.io/badge/J.Phys.Conf.Ser.-1342(2020)no.1,012076-255773.svg
      :target: https://iopscience.iop.org/article/10.1088/1742-6596/1342/1/012076
      :alt: J.Phys.Conf.Ser. 1342 (2020) no.1, 012076
  .. image:: https://img.shields.io/badge/arXiv-1711.07692-B31B1B.svg
      :target: https://arxiv.org/abs/1711.07692
      :alt: [arXiv:1711.07692]


* **Description of CRESST-III data**
  
  CRESST Collaboration (A.H. Abdelhameed et al.) 

  .. image:: https://img.shields.io/badge/arXiv-1905.07335-B31B1B.svg
      :target: https://arxiv.org/abs/1905.07335
      :alt: [arXiv:1905.07335]


CRESST-surface
^^^^^^^^^^^^^^

* **Results on MeV-scale dark matter from a gram-scale cryogenic calorimeter operated above ground**
  
  CRESST Collaboration (G. Angloher et al.)  

  .. image:: https://img.shields.io/badge/Eur.Phys.J.-C77(2017)no.9,637-255773.svg
      :target: https://link.springer.com/article/10.1140%2Fepjc%2Fs10052-017-5223-9
      :alt: Eur.Phys.J. C77 (2017) no.9, 637
  .. image:: https://img.shields.io/badge/arXiv-1707.06749-B31B1B.svg
      :target: https://arxiv.org/abs/1707.06749
      :alt: [arXiv:1707.06749]


DAMIC_N_2012
^^^^^^^^^^^^

* **Direct Search for Low Mass Dark Matter Particles with CCDs**
  
  DAMIC Collaboration (J. Barreto et al.)  

  .. image:: https://img.shields.io/badge/Phys.Lett.B-711(2012)264-255773.svg
      :target: https://www.sciencedirect.com/science/article/pii/S0370269312003887?via%3Dihub
      :alt: Phys.Lett. B711 (2012) 264
  .. image:: https://img.shields.io/badge/arXiv-1105.5191-B31B1B.svg
      :target: https://arxiv.org/abs/1105.5191
      :alt: [arXiv:1105.5191]


XENON1T_N_2017
^^^^^^^^^^^^^^

* **First Dark Matter Search Results from the XENON1T Experiment**
  
  XENON Collaboration (E. Aprile et al.) 

  .. image:: https://img.shields.io/badge/Phys.Rev.Lett.-119(2017)no.18,181301-255773.svg
      :target: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.181301
      :alt: Phys.Rev.Lett. 119 (2017) no.18, 181301
  .. image:: https://img.shields.io/badge/arXiv-1705.06655-B31B1B.svg
      :target: https://arxiv.org/abs/1705.06655
      :alt: [arXiv:1705.06655]



---------------------------
Electron recoil experiments
---------------------------

CDMS-HVeV_2018
^^^^^^^^^^^^^^

* **First Dark Matter Constraints from a SuperCDMS Single-Charge Sensitive Detector**
  
  SuperCDMS Collaboration (R. Agnese et al.)

  .. image:: https://img.shields.io/badge/Phys.Rev.Lett.-121(2018)no.5,051301-255773.svg
      :target: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.051301
      :alt: Phys.Rev.Lett. 121 (2018) no.5, 051301
  .. image:: https://img.shields.io/badge/arXiv-1804.10697-B31B1B.svg
      :target: https://arxiv.org/abs/1804.10697
      :alt: [arXiv:1804.10697]


CDMS-HVeV_2010
^^^^^^^^^^^^^^

* **Constraints on low-mass, relic dark matter candidates from a surface-operated SuperCDMS single-charge sensitive detector**
  
  SuperCDMS Collaboration (D.W. Amaral et al.) 

  .. image:: https://img.shields.io/badge/Phys.Rev.Lett.-121(2018)no.5,051301-255773.svg
      :target: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.051301
      :alt: Phys.Rev.Lett. 121 (2018) no.5, 051301
  .. image:: https://img.shields.io/badge/arXiv-2005.14067-B31B1B.svg
      :target: https://arxiv.org/abs/2005.14067
      :alt: [arXiv:2005.14067]


DarkSide-50_S2
^^^^^^^^^^^^^^

* **Constraints on Sub-GeV Dark-Matter–Electron Scattering from the DarkSide-50 Experiment**
  
  DarkSide Collaboration (P. Agnes et al.) 

  .. image:: https://img.shields.io/badge/Phys.Rev.Lett.-121(2018)no.11,111303-255773.svg
      :target: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.111303
      :alt: Phys.Rev.Lett. 121 (2018) no.11, 111303
  .. image:: https://img.shields.io/badge/arXiv-1802.06998-B31B1B.svg
      :target: https://arxiv.org/abs/1802.06998
      :alt: [arXiv:1802.06998]


**protoSENSEI@surface**
^^^^^^^^^^^^^^^^^^^^^^^

* **SENSEI: First Direct-Detection Constraints on sub-GeV Dark Matter from a Surface Run**
  
  SENSEI Collaboration (Michael Crisler et al.)

  .. image:: https://img.shields.io/badge/Phys.Rev.Lett.-121(2018)no.6-255773.svg
      :target: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.121.061803
      :alt: Phys.Rev.Lett. 121 (2018) no.6, 061803
  .. image:: https://img.shields.io/badge/arXiv-1804.00088-B31B1B.svg
      :target: https://arxiv.org/abs/1804.00088
      :alt: [arXiv:1804.00088]


**protoSENSEI@MINOS**
^^^^^^^^^^^^^^^^^^^^^

* **SENSEI: Direct-Detection Constraints on Sub-GeV Dark Matter from a Shallow Underground Run Using a Prototype Skipper-CCD**
  
  SENSEI Collaboration (Orr Abramoff et al.) 

  .. image:: https://img.shields.io/badge/Phys.Rev.Lett.-122(2019)no.16,161801-255773.svg
      :target: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.161801
      :alt: Phys.Rev.Lett. 122 (2019) no.16, 161801
  .. image:: https://img.shields.io/badge/arXiv-1901.10478-B31B1B.svg
      :target: https://arxiv.org/abs/1901.10478
      :alt: [arXiv:1901.10478]


**SENSEI@MINOS**
^^^^^^^^^^^^^^^^

* **SENSEI: Direct-Detection Results on sub-GeV Dark Matter from a New Skipper-CCD**
  
  SENSEI Collaboration (Liron Barak et al.)

  .. image:: https://img.shields.io/badge/Phys.Rev.Lett.-125(2020)no.17,171802-255773.svg
      :target: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.171802
      :alt: Phys.Rev.Lett. 125 (2020) 17, 171802
  .. image:: https://img.shields.io/badge/arXiv-2004.11378-B31B1B.svg
      :target: https://arxiv.org/abs/2004.11378
      :alt: [arXiv:2004.11378]


XENON10_S2
^^^^^^^^^^

* **A search for light dark matter in XENON10 data**
  
  XENON10 Collaboration (J. Angle et al.)

  .. image:: https://img.shields.io/badge/Phys.Rev.Lett.-107(2011)051301-255773.svg
      :target: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.107.051301
      :alt: Phys.Rev.Lett. 107 (2011) 051301
  .. image:: https://img.shields.io/badge/arXiv-1104.3088-B31B1B.svg
      :target: https://arxiv.org/abs/1104.3088
      :alt: [arXiv:1104.3088]


* **First Direct Detection Limits on sub-GeV Dark Matter from XENON10**
  
  Rouven Essig, Aaron Manalaysay, Jeremy Mardon, Peter Sorensen, Tomer Volansky

  .. image:: https://img.shields.io/badge/Phys.Rev.Lett.-109(2012)021301-255773.svg
      :target: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.109.021301
      :alt: Phys.Rev.Lett. 109 (2012) 021301
  .. image:: https://img.shields.io/badge/arXiv-1206.2644-B31B1B.svg
      :target: https://arxiv.org/abs/1206.2644
      :alt: [arXiv:1206.2644]


* **New Constraints and Prospects for sub-GeV Dark Matter Scattering off Electrons in Xenon**
  
  Rouven Essig, Tomer Volansky, Tien-Tien Yu 

  .. image:: https://img.shields.io/badge/Phys.Rev.D-96(2017)no.4-255773.svg
      :target: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.96.043017
      :alt: Phys.Rev. D96 (2017) no.4, 043017
  .. image:: https://img.shields.io/badge/arXiv-1703.00910-B31B1B.svg
      :target: https://arxiv.org/abs/1703.00910
      :alt: [arXiv:1703.00910]


XENON100_S2
^^^^^^^^^^^

* **Low-mass dark matter search using ionization signals in XENON100**
  
  XENON Collaboration (E. Aprile et al.)  

  .. image:: https://img.shields.io/badge/Phys.Rev.D-94(2016)no.9-255773.svg
      :target: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.092001
      :alt: Phys.Rev. D94 (2016) no.9, 092001
  .. image:: https://img.shields.io/badge/arXiv-1605.06262-B31B1B.svg
      :target: https://arxiv.org/abs/1605.06262
      :alt: [arXiv:1605.06262]


* **New Constraints and Prospects for sub-GeV Dark Matter Scattering off Electrons in Xenon**
  
  Rouven Essig, Tomer Volansky, Tien-Tien Yu  

  .. image:: https://img.shields.io/badge/Phys.Rev.D-96(2017)no.4-255773.svg
      :target: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.96.043017
      :alt: Phys.Rev. D96 (2017) no.4, 043017
  .. image:: https://img.shields.io/badge/arXiv-1703.00910-B31B1B.svg
      :target: https://arxiv.org/abs/1703.00910
      :alt: [arXiv:1703.00910]


XENON1T_S2
^^^^^^^^^^

* **Light Dark Matter Search with Ionization Signals in XENON1T**
  
  XENON Collaboration (E. Aprile et al.)  

  .. image:: https://img.shields.io/badge/Phys.Rev.Lett.-123(2019)no.25,251801-255773.svg
      :target: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.251801
      :alt: Phys.Rev.Lett. 123 (2019) no.25, 251801
  .. image:: https://img.shields.io/badge/arXiv-1907.11485-B31B1B.svg
      :target: https://arxiv.org/abs/1907.11485
      :alt: [arXiv:1907.11485]

==================================
5. The ``DM_Distribution`` classes
==================================

In order to make predictions for direct detection experiments, the statistical properties of the incoming DM flux need to be specified. In particular, we need to know how many DM particles pass through the detector and with what energy.
In other words, we need to know the DM particle flux, or alternatively the local DM density and the energy distribution.



--------------------------
The interface / base class
--------------------------


The class ``DM_Distribution``, that is declared in `/include/obscura/DM_Distribution.hpp <https://github.com/temken/obscura/blob/master/include/obscura/DM_Distribution.hpp>`_, is an abstract base class or interface that defines all the functions we require to characterize a distribution and flux of DM particles.

Most importantly, the class provides interfaces to probability density functions (PDFs) for the DM particles' velocity or speed, their local energy density, differential particle flux, etc.

-----------------------------
The standard halo model (SHM)
-----------------------------

The conventional assumptions on the halo DM particles' properties is the Standard Halo Model (SHM).
The SHM describes the galactic DM by a truncated Maxwell-Boltzmann distribution. It is characterized by the following 4 parameters (for details see e.g. chapter 3.2 of [Emken2019]_

.. math::
	\rho_\chi, v_0, v_\mathrm{esc}, \mathbf{v}_\mathrm{obs}

In `/include/obscura/DM_Halo_Models.hpp <https://github.com/temken/obscura/blob/master/include/obscura/DM_Halo_Models.hpp>`_ we define the ``Standard_Halo_Model`` class which is an implemenation of this model.
It is a derived class of ``DM_Distribution``.

We can construct the SHM model by the default constructor, which assumes default values for the 4 parameters.

.. code-block:: c++

	#include "obscura/DM_Halo_Models.hpp"

	// ...

	obscura::Standard_Halo_Model shm;

Or we define the parameters explicitly.

.. code-block:: c++

	#include "libphysica/Natural_Units.hpp"

	#include "obscura/DM_Halo_Models.hpp"

	using namespace libphysica::natural_units;

	// ...
	double rho = 0.4 * GeV / cm / cm / cm;
	double v_0 = 230.0 * km / sec;
	double v_esc = 600 * km / sec;
	double v_obs = 232.0 * km / sec;
	obscura::Standard_Halo_Model shm(rho, v_0, v_obs, v_esc);

---------
The SHM++
---------

As a second example for a DM halo model, *obscura* also implements the SHM++ as proposed in [Evans2019]_.

Since it extends the SHM, the corresponding class ``SHM_Plus_Plus`` is a derived class of ``Standard_Halo_Model`` which is in turn derived from ``DM_Distribution``.
The class is also declared in `/include/obscura/DM_Halo_Models.hpp <https://github.com/temken/obscura/blob/master/include/obscura/DM_Halo_Models.hpp>`_.

This halo model can be constructed and used essentially identically to the SHM.

-------------------------
Imported DM distributions
-------------------------

It is also possible to import a DM distribution from a file.
This is the purpose of the ``Imported_DM_Distribution`` class, another derived class of ``DM_Distribution`` which can be found in `/include/obscura/DM_Distribution.hpp <https://github.com/temken/obscura/blob/master/include/obscura/DM_Distribution.hpp>`_.

As input file, we need a two-column table of the DM speed PDF using the format (v[km/sec] :: f(v) [sec/km]).
Additionally we need to specify the local DM density.

Here is an example of using this class assuming a tabulated speed pdf given in the file *DM_Speed_PDF.txt*.


.. code-block:: c++

	#include "libphysica/Natural_Units.hpp"

	#include "obscura/DM_Distribution.hpp"

	using namespace libphysica::natural_units;

	// ...
	double rho = 0.4 * GeV / cm / cm / cm;
	obscura::Imported_DM_Distribution dm_distribution(rho, "DM_Speed_PDF.txt");