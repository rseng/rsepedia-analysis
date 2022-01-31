Contributing to lenstronomy
===========================

GitHub Workflow
---------------

### Fork and Clone the lenstronomy Repository
**You should only need to do this step once**

First *fork* the lenstronomy repository. A fork is your own remote copy of the repository on GitHub. To create a fork:

  1. Go to the [lenstronomy GitHub Repository](https://github.com/sibirrer/lenstronomy)
  2. Click the **Fork** button (in the top-right-hand corner)
  3. Choose where to create the fork, typically your personal GitHub account

Next *clone* your fork. Cloning creates a local copy of the repository on your computer to work with. To clone your fork:

  ```bash
  git clone https://github.com/<your-account>/lenstronomy.git
  ```

Finally add the `lenstronomy` repository as a *remote*. This will allow you to fetch changes made to the codebase. To add the `skypyproject` remote:

  ```bash
  cd lenstronomy
  git remote add lenstronomyproject https://github.com/sibirrer/lenstronomy.git
  ```

### Install your local lenstronomy version

To enable that your new code gets accessible by python also outside of the development environment, 
make sure all previous versions of lenstronomy are uninstalled and then install your version of lenstronomy (aka add the software to the python path)

  ```bash
  cd lenstronomy
  python setup.py develop --user
  ```

Alternatively, create virtual environments for the development (recommended for advanced usage with multiple branches).



### Create a branch for your new feature

Create a *branch* off the `lenstronomyproject` main branch. Working on unique branches for each new feature simplifies the development, review and merge processes by maintining logical separation. To create a feature branch:

  ```bash
  git fetch lenstronomyproject
  git checkout -b <your-branch-name> lenstronomyproject/main
  ```

### Hack away!

Write the new code you would like to contribute and *commit* it to the feature branch on your local repository. Ideally commit small units of work often with clear and descriptive commit messages describing the changes you made. To commit changes to a file:

  ```bash
  git add file_containing_your_contribution
  git commit -m 'Your clear and descriptive commit message'
  ```

*Push* the contributions in your feature branch to your remote fork on GitHub:

  ```bash
  git push origin <your-branch-name>
  ```


**Note:** The first time you *push* a feature branch you will probably need to use `--set-upstream origin` to link to your remote fork:

  
  ```bash
  git push --set-upstream origin <your-branch-name>
  ```

### Open a Pull Request

When you feel that work on your new feature is complete, you should create a *Pull Request*. This will propose your work to be merged into the main lenstronomy repository.

  1. Go to [lenstronomy Pull Requests](https://github.com/sibirrer/lenstronomy/pulls)
  2. Click the green **New pull request** button
  3. Click **compare across forks**
  4. Confirm that the base fork is `sibirrer/lenstronomy` and the base branch is `main`
  5. Confirm the head fork is `<your-account>/lenstronomy` and the compare branch is `<your-branch-name>`
  6. Give your pull request a title and fill out the the template for the description
  7. Click the green **Create pull request** button

### Updating your branch

As you work on your feature, new commits might be made to the `sibirrer/lenstronomy` main branch. You will need to update your branch with these new commits before your pull request can be accepted. You can achieve this in a few different ways:

  - If your pull request has no conflicts, click **Update branch**
  - If your pull request has conflicts, click **Resolve conflicts**, manually resolve the conflicts and click **Mark as resolved**
  - *merge* the `lenstronomyproject` main branch from the command line:
    ```bash
    git fetch lenstronomyproject
    git merge lenstronomyproject/main
    ```
  - *rebase* your feature branch onto the `lenstronomy` main branch from the command line:
    ```bash
    git fetch lenstronomyproject
    git rebase lenstronomyproject/main
    ```

**Warning**: It is bad practice to *rebase* commits that have already been pushed to a remote such as your fork. Rebasing creates new copies of your commits that can cause the local and remote branches to diverge. `git push --force` will **overwrite** the remote branch with your newly rebased local branch. This is strongly discouraged, particularly when working on a shared branch where you could erase a collaborators commits.

For more information about resolving conflicts see the GitHub guides:
  - [Resolving a merge conflict on GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/resolving-a-merge-conflict-on-github)
  - [Resolving a merge conflict using the command line](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/resolving-a-merge-conflict-using-the-command-line)
  - [About Git rebase](https://help.github.com/en/github/using-git/about-git-rebase)

### More Information

More information regarding the usage of GitHub can be found in the [GitHub Guides](https://guides.github.com/).

Coding Guidelines
-----------------

Before your pull request can be merged into the codebase, it will be reviewed by one of the lenstronomy developers and required to pass a number of automated checks. Below are a minimum set of guidelines for developers to follow:

### General Guidelines

- lenstronomy is compatible with Python>=3.5 (see [setup.cfg](setup.cfg)). lenstronomy *does not* support backwards compatibility with Python 2.x; `six`, `__future__` and `2to3` should not be used.
- All contributions should follow the [PEP8 Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/). We recommend using [flake8](https://flake8.pycqa.org/) to check your code for PEP8 compliance.
- Importing lenstronomy should only depend on having [NumPy](https://www.numpy.org), [SciPy](https://www.scipy.org/) and [Astropy](https://www.astropy.org/) installed.
- Code is grouped into submodules e.g. [LensModel](lenstronomy/LensModel), [LightModel](lenstronomy/LightModel) and [ImSim](lenstronomy/ImSim). There is also a [Util](lenstronomy/Util) submodule for general utility functions.
- For more information see the [Astropy Coding Guidelines](http://docs.astropy.org/en/latest/development/codeguide.html)

### Unit Tests

Pull requests will require existing unit tests to pass before they can be merged. Additionally, new unit tests should be written for all new public methods and functions. Unit tests for each submodule are contained in subdirectories called `tests` and you can run them locally using `python setup.py test`. For more information see the [Astropy Testing Guidelines](https://docs.astropy.org/en/stable/development/testguide.html).

### Docstrings

All public classes, methods and functions require docstrings. You can build documentation locally by installing [sphinx-astropy](https://github.com/astropy/sphinx-astropy) and calling `python setup.py build_docs`. Docstrings should include the following sections:

  - Description
  - Parameters
  - Notes
  - Examples
  - References

For more information see the Astropy guide to [Writing Documentation](https://docs.astropy.org/en/stable/development/docguide.html).

This page is inspired by the Contributions guidelines of the [Skypy project](https://github.com/skypyproject/skypy/blob/main/CONTRIBUTING.md).
---
title: 'lenstronomy II: A gravitational lensing software ecosystem'
tags:
  - Python
  - astronomy
  - gravitational lensing
  - image simulations
  - dynamics
authors:
  - name: Simon Birrer
    orcid: 0000-0003-3195-5507
    affiliation: "1, 2"
  - name: Anowar J. Shajib
    orcid: 0000-0002-5558-888X
    affiliation: "3, 4"
  - name: Daniel Gilman
    orcid: 0000-0002-5116-7287
    affiliation: 5
  - name: Aymeric Galan
    orcid: 0000-0003-2547-9815
    affiliation: 6
  - name: Jelle Aalbers
    orcid: 0000-0003-0030-0030
    affiliation: "1, 2"   
  - name: Martin Millon
    orcid: 0000-0001-7051-497X
    affiliation: 6  
  - name: Robert Morgan
    orcid: 0000-0002-7016-5471
    affiliation: "7, 8"  
  - name: Giulia Pagano
    orcid: 0000-0002-3636-0767
    affiliation: 9
  - name: Ji Won Park
    orcid: 0000-0002-0692-1092
    affiliation: "1, 2"
  - name: Luca Teodori
    affiliation: 10
  - name: Nicolas Tessore
    orcid: 0000-0002-9696-7931
    affiliation: 11
  - name: Madison Ueland
    affiliation: 1
  - name: Lyne Van de Vyvere
    orcid: 0000-0002-0585-4203
    affiliation: 12
  - name: Sebastian Wagner-Carena
    orcid: 0000-0001-5039-1685
    affiliation: "1, 2"
  - name: Ewoud Wempe
    orcid: 0000-0001-8232-4188
    affiliation: 13
  - name: Lilan Yang
    orcid: 0000-0002-8434-880X
    affiliation: 14
  - name: Xuheng Ding
    orcid: 0000-0001-8917-2148
    affiliation: 15
  - name: Thomas Schmidt
    orcid: 0000-0002-2772-8160
    affiliation: 4
  - name: Dominique Sluse 
    orcid: 0000-0001-6116-2095
    affiliation: 12
  - name: Ming Zhang
    affiliation: 16
  - name: Adam Amara
    orcid: 0000-0003-3481-3491
    affiliation: 17

 
 
affiliations:
 - name: Kavli Institute for Particle Astrophysics and Cosmology and Department of Physics, Stanford University, Stanford, CA 94305, USA
   index: 1
 - name: SLAC National Accelerator Laboratory, Menlo Park, CA, 94025, USA
   index: 2
 - name: Department of Astronomy & Astrophysics, University of Chicago, Chicago, IL 60637, USA
   index: 3
 - name: Department of Physics and Astronomy, University of California, Los Angeles, CA 90095, USA
   index: 4 
 - name: Department of Astronomy and Astrophysics, University of Toronto, 50 St. George Street, Toronto, ON, M5S 3H4, Canada
   index: 5 
 - name: Institute of Physics, Laboratory of Astrophysics, Ecole Polytechnique Fédérale de Lausanne (EPFL), Switzerland
   index: 6
 - name: Physics Department, University of Wisconsin-Madison, 1150 University Avenue Madison, WI  53706, USA
   index: 7
 - name: Legacy Survey of Space and Time Corporation Data Science Fellowship Program, USA
   index: 8
 - name: Independent Researcher
   index: 9
 - name: Weizmann Institute, 234 Herzl Street, Rehovot, 7610001 Israel
   index: 10
 - name: Department of Physics and Astronomy, University College London, Gower Street, London, WC1E 6BT, UK
   index: 11
 - name: STAR Institute, Université de Liège, Quartier Agora - Allée du six Août, 19c, B-4000 Liège, Belgium
   index: 12
 - name: Kapteyn Astronomical Institute, University of Groningen, PO Box 800, 9700 AV Groningen, the Netherlands
   index: 13
 - name: School of Physics and Technology, Wuhan University, Wuhan 430072, China
   index: 14
 - name: Kavli IPMU (WPI), UTIAS, The University of Tokyo, Kashiwa, Chiba 277-8583, Japan
   index: 15
 - name: Xinjiang Astronomical Observatory, Chinese Academy of Sciences, 150 Science 1-Street, Urumqi 831001, China
   index: 16
 - name: Institute of Cosmology and Gravitation, University of Portsmouth, Portsmouth PO1 3FX, UK
   index: 17

   
   
date: 28 April 2021
codeRepository: https://github.com/sibirrer/lenstronomy
license: MIT
bibliography: paper.bib
---

# Summary

`lenstronomy` is an Astropy-affiliated [@astropy:2013; @astropy:2018] Python package for gravitational lensing simulations and analyses.
`lenstronomy` was introduced by @lenstronomy1 and is based on the linear basis set approach by @Birrer:2015.
The user and developer base of `lenstronomy` has substantially grown since then, and the software has become an integral part of a wide range of recent analyses, such as measuring the Hubble constant with time-delay strong lensing or constraining the nature of dark matter from resolved and unresolved small scale lensing distortion statistics. 
The modular design has allowed the community to incorporate innovative new methods, as well as to develop enhanced software and wrappers with more specific aims on top of the `lenstronomy` API.
Through community engagement and involvement, `lenstronomy` has become a foundation of an ecosystem of affiliated packages extending the original scope of the software and proving its robustness and applicability
at the forefront of the strong gravitational lensing community in an open source and reproducible manner.


![Illustration of the strong gravitational lensing phenomenology and the capabilities of lenstronomy in performing realistic simulations as well as reconstructing lensing and source properties from a given data set. Top row from left to right along the green arrow:
A galaxy is lensed around a foreground massive object, becomes highly distorted, and has components appearing multiple times. Observations of this phenomena are limited in resolution (convolution), depending on the detector (pixelation), and are subject to noise.
Bottom row from right to left along the red arrow: The inverse problem is solved with a linear basis set in the source morphology maximizing the likelihood of the model given the data.
\label{fig:example}](paper_fig.png)

# Background

Gravitational lensing displaces the observed positions and distorts the shapes of apparent objects on the sky due to intervening inhomogeneous matter along the line of sight. Strong gravitational lensing describes the regime where the background source, such as a galaxy or quasar, is lensed by a massive foreground object, such as another galaxy or cluster of galaxies, to produce multiple images of the source in a highly distorted manner. 
The top row of \autoref{fig:example} illustrates such a process from the intrinsic galaxy to the data product at hand, including the lensing distortions, effects of the instrument, observational conditions, and noise.

Analyses of strong gravitational lensing have provided a wealth of key insights into cosmology and astrophysics.
For example, relative time delays of multiply imaged variable sources provided precision measurements on the expansion rate of the Universe [@Wong:2020; @Shajib:2020strides; @Birrer:2020tdcosmoiv]. Small scale distortions in the lensing signal of resolved sources [@Vegetti:2012; @Hezaveh:2016; @Birrer:2017]
and unresolved flux ratios [@Gilman:2020; @Hsueh:2020] constrain the nature of dark matter. Combined strong lensing and kinematic observables constrain the formation and evolution of galaxies [@Sonnenfeld:2015; @Shajib:2021slacs], and the lensing magnification effect provides an otherwise inaccessible angle on the early Universe [@Zheng:2012; @Cava:2018].




# Statement of need

Strong lensing studies have significantly enhanced, and sometimes challenged, our current fundamental understanding of the Universe.
In the near future, with the onset of the next-generation ground and space-based wide and deep astronomical imaging [Rubin, Roman, Euclid observatories; @LSST; @Roman; @Euclid] and interferometric [SKA; @SKA] surveys, the number of discovered lenses of different types will be growing by more than an order of magnitude [@Collett:2015; @OM10].
Such large samples can provide unprecedented statistical precision to stress-test our current understanding and exploit discovery potential.
It is key that these demanding studies, at present and in the future, are conducted by reliable software and supported by reproducible and open-source analysis products to provide the most compelling and transparent evidence required to further our physical understanding.

The primary design goal of `lenstronomy` is to facilitate scientific investigations into the outstanding and most pressing questions in the cosmology and astrophysics community.
`lenstronomy` has been applied throughout its development to the most demanding modeling and inference problems in strong lensing and the software has evolved around the requirements of the scientific applications to facilitate robust analyses. The modular API of the original design of lenstronomy [@lenstronomy1] has accommodated the addition of new features. Code review processes in the development phase have led to additional benefits for the user community at large beyond the specific needs of the developer.

`lenstronomy` provides reliable and well-tested specific functionalities, as well as top-level interfaces, which allow for adaptive and innovative usage in control by the scientific investigator.
Guidance for the user community is provided on multiple levels. First, source code is well documented and provided through [readthedocs.org](http://lenstronomy.readthedocs.org). Second, a set of `jupyter` notebooks are provided in an [extension repository](https://github.com/sibirrer/lenstronomy_extensions). These notebooks demonstrate simplified example use cases, each notebook individually highlighting different specific functionalities of `lenstronomy`, including a [starting guide notebook](https://github.com/sibirrer/lenstronomy_extensions/blob/v1.8.1/lenstronomy_extensions/Notebooks/starting_guide.ipynb)  to introduce the modular design structure of  `lenstronomy`. Third, end-to-end analysis pipelines of some of the published work are publicly available, providing ‘real-life’ examples at advanced levels.


# Track-record of applications

`lenstronomy` has been applied in and contributed to more than 30 peer reviewed publications since its first public release in 2018.
In particular, `lenstronomy` has been used to provide state-of-the-art measurements on real data sets, such as: 
(i) Hubble constant measurements from three quadruly lensed quasars with Hubble Space Telescope (HST) imaging [@Birrer:2016; @Birrer:2019; @Shajib:2020strides], 
dynamical modeling in the hierarchical analysis by @Birrer:2020tdcosmoiv, and modeling of lensed supernovae [@Moertsell:2020]; 
(ii) inference of small scale dark matter properties from detailed studies of both, resolved imaging [@Birrer:2017], and unresolved flux ratio statistics [@Gilman:2020]; 
(iii) decomposition of quasar and host galaxy light in both, lensed and unlensed cases [@Ding:2020; @Bennert:2021]; 
(iv) morphological studies of high-redshift sources in the cluster environment [@Yang:2020; @Yang:2021]; 
(v) internal structure of galaxies [@Shajib:2021slacs; @Shajib:2021AO]; 
(vi) measurements of the weak lensing effect imprinted in Einstein rings [@Birrer:2017cosmos; @Kuhn:2021].
Among the studies, some of them have applied a pipeline to uniformly analyse dozens of lenses of different types [@Shajib:2019; @Shajib:2021slacs; @Shajib:2021AO], 
a milestone in moving towards utilizing thousands of lenses in the near future.

Beyond analyzing data, many theoretical studies have been conducted using `lenstronomy` to investigate statistical robustness in present and anticipated future analyses [@BirrerTreu:2019; @Millon:2020; @vdVyvere:2020; @Li:2021; @Ding:2021transient], 
as well as to provide forecasts for anticipated future constraints for different science cases [@Gilman:2019; @Sengul:2020; @BirrerTreu:2021].
Particularly, three separate teams participated in the blind time-delay lens modeling challenge [@Ding:2021tdlmc] using `lenstronomy`.

`lenstronomy` has seen a substantial development and incorporation of innovations and numerical recipes [@Tessore:2015; @Shajib:2019unified; @Joseph:2019; @Galan:2021; @Birrer:2021arcs], 
and has found applications beyond its original aim due to the robust and high-standard design requirements.


# Ecosystem of affiliated packages

`lenstronomy` has allowed the community to develop third-party analysis products and software products utilizing its core functionalities to provide more targeted and integrated software solutions for a wide range of scientific analyses. 
These open-source [affiliated packages](https://github.com/sibirrer/lenstronomy/blob/1.8.1/AFFILIATEDPACKAGES.rst) effectively create an ecosystem enhancing the capability of `lenstronomy`. 
They provide specified and tested solution for specific scientific investigations, such as plug-ins and direct implementation for innovative source reconstruction algorithms [[SLITronomy](https://github.com/aymgal/SLITronomy); @Joseph:2019; @Galan:2021], 
gravitational wave lensing computations [[lensingGW](https://gitlab.com/gpagano/lensinggw); @Pagano:2020], 
automated pipelines for gravitational lensing reconstruction [[dolphin](https://github.com/ajshajib/dolphin); @Shajib:2021slacs], 
cluster source reconstruction and local perturbative lens modeling [[lenstruction](https://github.com/ylilan/lenstruction); @Yang:2020], 
enhancement in large-scale structure imaging survey simulations [[DESC SLSprinkler](https://github.com/LSSTDESC/SLSprinkler); @LSSTDESC:2021], 
rendering of sub-halos and line-of-sight halos [[pyHalo](https://github.com/dangilman/pyHalo); @Gilman:2020], 
galaxy morphology analysis [[galight](https://github.com/dartoon/galight); @Ding:2020],
and hierarchical analyses to measure the Hubble constant [[hierArc](https://github.com/sibirrer/hierarc); @Birrer:2020tdcosmoiv].
With the rise in popularity and the promises in dealing with ever complex data problems with fast deep-learning methods, 
dedicated tools for simulating large datasets for applying such methods to strong gravitational lensing [[deeplenstronomy](https://github.com/deepskies/deeplenstronomy); @Morgan:2021], [[baobab](https://github.com/jiwoncpark/baobab); @Park:2021], 
as well as end-to-end Bayesian Neural Network training and validation packages for Hubble constant measurements [[h0rton](https://github.com/jiwoncpark/h0rton); @Park:2021], 
and for a hierarchical analysis of galaxy-galaxy lenses [[ovejero](https://github.com/swagnercarena/ovejero); @Wagner-Carena:2021]
 have been developed.
The affiliated packages make best use of the `lenstronomy` modules without duplicating source code and make it possible to combine aspects of multiple affiliated packages in one single analysis.


# Related open source software

- [`lenstronomy`](https://github.com/sibirrer/lenstronomy) [@Birrer:2015; @lenstronomy1]
- [`PyAutoLens`](https://github.com/Jammy2211/PyAutoLens) [@Nightingale:2018; @Nightingale:2021]
- [`gravlens`](http://www.physics.rutgers.edu/~keeton/gravlens/) [@Keeton:2011]
- [`glafic`](https://www.slac.stanford.edu/~oguri/glafic/) [@Oguri:2010]
- [`visilens`](https://github.com/jspilker/visilens) [@spilker16a]
- [`PixeLens`](https://www.physik.uzh.ch/~psaha/lens/pixelens.php) [@PixeLens]
- [`GRALE`](https://github.com/j0r1/GRALE2) [@GRALE]
- [`lenstool`](http://projets.lam.fr/projects/lenstool/wiki) [@Jullo:2009]



# Acknowledgements

Support for this work was provided by the National Science Foundation through NSF AST-1716527. 
AJS was supported by NASA through the STScI grant HST-GO-15320 and by a Dissertation Year Fellowship from the UCLA Graduate Division. 
This research was supported by the U.S. Department of Energy (DOE) Office of Science Distinguished Scientist Fellow Program.
DG is supported by NASA HST-GO-15177.
AG, MM LvdV, DS are supported by COSMICLENS: ERC grant agreement No 787886.
LT is supported by International Helmholtz-Weizmann Research School for Multimessenger Astronomy.
MU is supported by KIPAC and the Stanford Summer Research Program.
XD is supported by NASA HST-GO-15115.
TS is supported by NASA grant HST-GO-15320 and HST-GO-15652.
MZ is supported by the National Science Foundation of China.
AA is supported by a Royal Society Wolfson Fellowship.
We are grateful to the user community for valuable feedback and encouragement in continuing the development.


# References
============
Mailing list
============

You can join the **lenstronomy** mailing list by signing up on the
`google groups page <https://groups.google.com/forum/#!forum/lenstronomy>`_.

The email list is meant to provide a communication platform between users and developers. You can ask questions,
and suggest new features. New releases will be announced via this mailing list.

If you encounter errors or problems with **lenstronomy**, please let us know!===================
Affiliated packages
===================

Here is an (incomplete) list of packages and wrappers that are using lenstronomy in various ways for specific scientific
applications:

- `baobab <https://github.com/jiwoncpark/baobab>`_: Training data generator for hierarchically modeling of strong lenses with Bayesian neural networks.
- `dolphin <https://github.com/ajshajib/dolphin>`_: Automated pipeline for lens modeling based on lenstronomy.
- `hierArc <https://github.com/sibirrer/hierarc>`_: Hierarchical Bayesian time-delay cosmography to infer the Hubble constant and galaxy density profiles in conjunction with lenstronomy.
- `lenstruction <https://github.com/ylilan/lenstruction>`_: Versatile tool for cluster source reconstruction and local perturbative lens modeling.
- `SLITronomy <https://github.com/aymgal/SLITronomy>`_: Updated and improved version of the Sparse Lens Inversion Technique (SLIT), developed within the framework of lenstronomy.
- `LSSTDESC SLSprinkler <https://github.com/LSSTDESC/SLSprinkler>`_: The DESC SL (Strong Lensing) Sprinkler adds strongly lensed AGN and SNe to simulated catalogs and generates postage stamps for these systems.
- `lensingGW <https://gitlab.com/gpagano/lensinggw>`_: A Python package designed to handle both strong and microlensing of compact binaries and the related gravitational-wave signals.
- `ovejero <https://github.com/swagnercarena/ovejero>`_: Conducts hierarchical inference of strongly-lensed systems with Bayesian neural networks.
- `h0rton <https://github.com/jiwoncpark/h0rton>`_: H0 inferences with Bayesian neural network lens modeling.
- `deeplenstronomy <https://github.com/deepskies/deeplenstronomy>`_: Tool for simulating large datasets for applying deep learning to strong gravitational lensing.
- `pyHalo <https://github.com/dangilman/pyHalo>`_: Tool for rendering full substructure mass distributions for gravitational lensing simulations.
- `GaLight <https://github.com/dartoon/galight>`_: Tool to perform two-dimensional model fitting of optical and near-infrared images to characterize surface brightness distributions.



These packages come with their own documentation and examples - so check them out!



Guidelines for affiliated packages
----------------------------------
If you have a package/wrapper/analysis pipeline that is open source and you would like to have it advertised here, please let the developers know!
Before you write your own wrapper and scripts in executing lenstronomy for your purpose check out the list
of existing add-on packages. Affiliated packages should not duplicate the core routines of lenstronomy and whenever possible make use of the lenstronomy modules.
The packages should be maintained to keep up with the development of lenstronomy. Please also make sure the citation guidelines are presented.=======
Credits
=======

Development Lead
----------------

* Simon Birrer <sibirrer@gmail.com> `sibirrer <https://github.com/sibirrer/>`_


Key contributors
----------------
* Anowar Shajib `ajshajib <https://github.com/ajshajib/>`_
* Aymeric Galan `aymgal <https://github.com/aymgal/>`_
* Daniel Gilman `dangilman <https://github.com/dangilman/>`_


Contributors (alphabetic)
-------------------------

* Jelle Aalbers `JelleAalbers <https://github.com/JelleAalbers>`_
* Joel Akeret `jakeret <https://github.com/jakeret/>`_
* Adam Amara `aamara <https://github.com/aamara/>`_
* Xuheng Ding `dartoon <https://github.com/dartoon/>`_
* Sydney Erickson `smericks <https://github.com/smericks/>`_
* Kevin Fusshoeller
* Matthew R. Gomer `mattgomer <https://github.com/mattgomer>`_
* Felix A. Kuhn
* Felix Mayor
* Martin Millon `martin-millon <https://github.com/martin-millon/>`_
* Robert Morgan `rmorgan10 <https://github.com/rmorgan10/>`_
* Anna Nierenberg `amn3142 <https://github.com/amn3142/>`_
* Brian Nord `bnord <https://github.com/bnord/>`_
* Giulia Pagano
* Ji Won Park `jiwoncpark <https://github.com/jiwoncpark/>`_
* Andreas Filipp `andreasfilipp <https://github.com/andreasfilipp/>`_
* Thomas Schmidt `Thomas-01 <https://github.com/Thomas-01/>`_
* Dominique Sluse
* Luca Teodori `lucateo <https://github.com/lucateo/>`_
* Nicolas Tessore `ntessore <https://github.com/ntessore/>`_
* Madison Ueland `mueland <https://github.com/mueland/>`_
* Lyne Van de Vyvere `LyneVdV <https://github.com/LyneVdV/>`_
* Sebastian Wagner-Carena `swagnercarena <https://github.com/swagnercarena>`_
* Cyril Welschen
* Ewoud Wempe `ewoudwempe <https://github.com/ewoudwempe/>`_
* Lilan Yang `ylilan <https://github.com/ylilan/>`_
* Zhiyuan Ma `Jerry-Ma <https://github.com/Jerry-Ma/>`_
.. :changelog:

History
-------

0.0.1 (2018-01-09)
++++++++++++++++++

* First release on PyPI.

0.0.2 (2018-01-16)
++++++++++++++++++

* Improved testing and stability

0.0.6 (2018-01-29)
++++++++++++++++++

* Added feature to align coordinate system of different images

0.1.0 (2018-02-25)
++++++++++++++++++

* Major design update

0.1.1 (2018-03-05)
++++++++++++++++++

* minor update to facilitate options without lensing

0.2.0 (2018-03-10)
++++++++++++++++++

* ellipticity parameter handling changed
* time-delay distance sampling included
* parameter handling for sampling more flexible
* removed redundancies in the light and mass profiles

0.2.1 (2018-03-19)
++++++++++++++++++

* updated documentation
* improved sub-sampling of the PSF

0.2.2 (2018-03-25)
++++++++++++++++++

* improved parameter handling
* minor bugs with parameter handling fixed

0.2.8 (2018-05-31)
++++++++++++++++++

* improved GalKin module
* minor improvements in PSF reconstruction
* mass-to-light ratio parameterization

0.3.1 (2018-07-21)
++++++++++++++++++

* subgrid psf sampling for inner parts of psf exclusively
* minor stability improvements
* cleaner likelihood definition
* additional Chameleon lens and light profiles

0.3.3 (2018-08-21)
++++++++++++++++++
* minor updates, better documentation and handling of parameters

0.4.1-3 (2018-11-27)
++++++++++++++++++++
* various multi-band modelling frameworks added
* lens models added
* Improved fitting sequence, solver and psf iteration

0.5.0 (2019-1-30)
+++++++++++++++++
* Workflow module redesign
* improved parameter handling
* improved PSF subsampling module
* relative astrometric precision of point sources implemented

0.6.0 (2019-2-26)
+++++++++++++++++
* Simulation API module for mock generations
* Multi-source plane modelling

0.7.0 (2019-4-13)
+++++++++++++++++
* New design of Likelihood module
* Updated design of FittingSequence
* Exponential Shapelets implemented

0.8.0 (2019-5-23)
+++++++++++++++++
* New design of Numerics module
* New design of PSF and Data module
* New design of multi-band fitting module

0.8.1 (2019-5-23)
+++++++++++++++++
* PSF numerics improved and redundancies removed.

0.8.2 (2019-5-27)
+++++++++++++++++
* psf_construction simplified
* parameter handling for catalogue modelling improved

0.9.0 (2019-7-06)
+++++++++++++++++
* faster fft convolutions
* re-design of multi-plane lensing module
* re-design of plotting module
* nested samplers implemented
* Workflow module with added features

0.9.1 (2019-7-21)
+++++++++++++++++
* non-linear solver for 4 point sources updated
* new lens models added
* updated Workflow module
* implemented differential extinction

0.9.2 (2019-8-29)
+++++++++++++++++
* non-linear solver for 4 point sources updated
* Moffat PSF for GalKin in place
* Likelihood module for point sources and catalogue data improved
* Design improvements in the LensModel module
* minor stability updates

0.9.3 (2019-9-25)
+++++++++++++++++
* improvements in SimulationAPI design
* improvements in astrometric uncertainty handling of parameters
* local arc lens model description and differentials


1.0.0 (2019-9-25)
+++++++++++++++++
* marking version as 5 - Stable/production mode

1.0.1 (2019-10-01)
++++++++++++++++++
* compatible with emcee 3.0.0
* removed CosmoHammer MCMC sampling

1.1.0 (2019-11-5)
+++++++++++++++++
* plotting routines split in different files
* curved arc parameterization and eigenvector differentials
* numerical differentials as part of the LensModel core class


1.2.0 (2019-11-17)
++++++++++++++++++
* Analysis module re-designed
* GalKin module partially re-designed
* Added cosmography module
* parameterization of cartesian shear coefficients changed


1.2.4 (2020-01-02)
++++++++++++++++++
* First implementation of a LightCone module for numerical ray-tracing
* Improved cosmology sampling from time-delay cosmography measurements
* TNFW profile lensing potential implemented


1.3.0 (2020-01-10)
++++++++++++++++++
* image position likelihood description improved


1.4.0 (2020-03-26)
++++++++++++++++++
* Major re-design of GalKin module, added new anisotropy modeling and IFU aperture type
* Updated design of the Analysis.kinematicsAPI sub-module
* Convention and redundancy in the Cosmo module changed
* NIE, SIE and SPEMD model consistent with their ellipticity and Einstein radius definition
* added cored-Sersic profile
* dependency for PSO to CosmoHammer removed
* MPI and multi-threading for PSO and MCMC improved and compatible with python3


1.5.0 (2020-04-05)
++++++++++++++++++
* Re-naming SPEMD to PEMD, SPEMD_SMOOTH to SPEMD
* adaptive numerics improvement
* multi-processing improvements


1.5.1 (2020-06-20)
++++++++++++++++++
* bug fix in Hession of POINT_SOURCE model
* EPL model from Tessore et al. 2015 implemented
* multi-observation mode for kinematics calculation


1.6.0 (2020-09-07)
++++++++++++++++++
* SLITronomy integration
* observation configuration templates and examples
* lens equation solver arguments in single sub-kwargs
* adapted imports to latest scipy release
* iterative PSF reconstruction improved
* multipole lens model


1.7.0 (2020-12-16)
++++++++++++++++++
* cosmo.NFWParam mass definition changed
* QuadOptimizer re-factored
* interpol light model support for non-square grid
* add functionality to psf error map
* fix in multiband reconstruction
* observational config for ZTF
* short-hand class imports


1.8.0 (2020-03-21)
++++++++++++++++++
* EPL numba version
* numba configuration variables can be set globally with configuration file
* Series of curved arc models available
* single plane hessian return all for differentials
* elliptical density slice lens model
* vectorized lens and light interpolation models
* updated installation description
* fast caustic calculation replacing matplotlib with skitlearn
* multi-patch illustration class and plotting routines
* updated PSF iteration procedure with more settings

1.8.1 (2020-04-19)
++++++++++++++++++
* illustration plots for curved arcs updated
* documentation of elliptical lens models updated


1.8.2 (2020-06-08)
++++++++++++++++++
* JOSS paper added
* improved testing documentation and tox compatibility
* TNFW_ELLIPSE lens model implemented
* ULDM lens model implemented


1.9.0 (2020-07-15)
++++++++++++++++++
* re-defined half light radius in Sersic profile
* re-named parameter in 'CONVERGENCE' profile
* improved numerics in Galkin
* configuration import design changed


1.9.1 (2020-08-27)
++++++++++++++++++
* re-defined amplitude normalization in NIE and CHAMELEON light profiles
* bug fix in sky brightness errors (SimulationAPI)


1.9.2 (2020-12-12)
++++++++++++++++++
* support for astropy v5
* new PSF iteration procedure implemented
* updated caustic plotting feature
* magnification perturbations in point source amplitudes
* analytic point source solver for SIE+shear

1.9.3 (2020-12-22)
++++++++++++++++++
* changed syntax to be compatible with python3 version <3.9===============================
Published work with lenstronomy
===============================

In this section you can find the concept papers **lenstronomy** is based on the list of science publications that made
use of **lenstronomy**. Please let the developers know when you publish a paper that made use of **lenstronomy**.
We are happy to include your publication in this list.



Core lenstronomy methodology and software publications
------------------------------------------------------

* lenstronomy: Multi-purpose gravitational lens modelling software package; `Birrer & Amara 2018 <https://ui.adsabs.harvard.edu/abs/2018PDU....22..189B>`_
    *This is the lenstronomy software paper. Please cite this paper whenever you make use of lenstronomy. The paper gives a design overview and highlights some use cases.*

* lenstronomy II: A gravitational lensing software ecosystem; `Birrer et al. 2021 <https://joss.theoj.org/papers/10.21105/joss.03283>`_
    *JOSS software publication. Please cite this paper whenever you make use of lenstronomy.*

* Gravitational Lens Modeling with Basis Sets; `Birrer et al. 2015 <http://adsabs.harvard.edu/abs/2015ApJ...813..102B>`_
    *This is the method paper lenstronomy is primary based on. Please cite this paper whenever you publish results with lenstronomy by using Shapelet basis sets and/or the PSO and MCMC chain.*


Related software publications
-----------------------------

* A versatile tool for cluster lensing source reconstruction. I. methodology and illustration on sources in the Hubble Frontier Field Cluster MACS J0717.5+3745; `Yang et al. 2020a <https://arxiv.org/abs/2001.07719>`_
    *reconstructing the intrinsic size-mass relation of strongly lensed sources in clusters*

* SLITronomy: towards a fully wavelet-based strong lensing inversion technique; `Galan et al. 2020 <https://arxiv.org/abs/2012.02802>`_
    *This is the method paper presenting SLITromomy, an improved version of the SLIT algorithm fully implemented and compatible with lenstronomy.*

* deeplenstronomy: A dataset simulation package for strong gravitational lensing; `Morgan et al. 2021a <https://arxiv.org/abs/2102.02830>`_
    *Software to simulating large datasets for applying deep learning to strong gravitational lensing.*

* Galaxy shapes of Light (GaLight): a 2D modeling of galaxy images; `Ding et al. 2021b <https://arxiv.org/abs/2111.08721>`_
    *Tool to perform two-dimensional model fitting of optical and near-infrared images to characterize surface brightness distributions.*




Measuring the Hubble constant
-----------------------------

* The mass-sheet degeneracy and time-delay cosmography: analysis of the strong lens RXJ1131-1231; `Birrer et al. 2016 <http://adsabs.harvard.edu/abs/2016JCAP...08..020B>`_
    *This paper performs a cosmographic analysis and applies the Shapelet basis set scaling to marginalize over a major lensing degeneracy.*

* H0LiCOW - IX. Cosmographic analysis of the doubly imaged quasar SDSS 1206+4332 and a new measurement of the Hubble constant; `Birrer et al. 2019 <https://ui.adsabs.harvard.edu/#abs/2018arXiv180901274B/abstract>`_
    *This paper performs a cosmographic analysis with power-law and composite models and covers a range in complexity in the source reconstruction*

* Astrometric requirements for strong lensing time-delay cosmography; `Birrer & Treu 2019 <https://ui.adsabs.harvard.edu/abs/2019MNRAS.tmp.2172B>`_
    *Derives requirements on how well the image positions of time-variable sources has to be known to perform a time-delay cosmographic measurement*

* H0LiCOW XIII. A 2.4% measurement of  H0 from lensed quasars: 5.3σ tension between early and late-Universe probes; `Wong et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019arXiv190704869W>`_
    *Joint analysis of the six H0LiCOW lenses including the lenstronomy analysis of J1206*

* STRIDES: A 3.9 per cent measurement of the Hubble constant from the strongly lensed system DES J0408-5354; `Shajib et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019arXiv191006306S/abstract>`_
    *most precise single lensing constraint on the Hubble constant. This analysis includes two source planes and three lensing planes*

* TDCOSMO. I. An exploration of systematic uncertainties in the inference of H0 from time-delay cosmography `Millon et al. 2020 <https://ui.adsabs.harvard.edu/abs/2019arXiv191208027M/abstract>`_
    *mock lenses to test accuracy on the recovered H0 value*

* Lens modelling of the strongly lensed Type Ia supernova iPTF16geu `Moertsell et al. 2020 <https://ui.adsabs.harvard.edu/abs/2019arXiv190706609M/abstract>`_
    *Modeling of a lensed supernova to measure the Hubble constant*

* The impact of line-of-sight structures on measuring H0 with strong lensing time-delays `Li, Becker and Dye 2020 <https://arxiv.org/abs/2006.08540v1>`_
    *Point source position and time-delay modeling of quads*

* TDCOSMO III: Dark matter substructure meets dark energy -- the effects of (sub)halos on strong-lensing measurements of H0 `Gilman, Birrer and Treu 2020 <https://ui.adsabs.harvard.edu/abs/2020arXiv200701308G/abstract>`_
    *Full line-of-sight halo rendering and time-delay analysis on mock images*

* TDCOSMO IV: Hierarchical time-delay cosmography -- joint inference of the Hubble constant and galaxy density profiles `Birrer et al. 2020 <https://arxiv.org/abs/2007.02941>`_
    *lenstronomy.Galkin for kinematics calculation that folds in the hierarchical analysis*

* TDCOSMO V: strategies for precise and accurate measurements of the Hubble constant with strong lensing `Birrer & Treu 2020 <https://ui.adsabs.harvard.edu/abs/2020arXiv200806157B/abstract>`_
    *lenstronomy.Galkin for kinematics calculation that folds in the hierarchical analysis for a forecast for future Hubble constant constraints*

* Large-Scale Gravitational Lens Modeling with Bayesian Neural Networks for Accurate and Precise Inference of the Hubble Constant `Park et al. 2020 <https://arxiv.org/abs/2012.00042>`_
    *BBN lens model inference using lenstronomy through `baobab <https://github.com/jiwoncpark/baobab>`_ for training set generation.*

* Improved time-delay lens modelling and H0 inference with transient sources `Ding et al. 2021a <https://arxiv.org/abs/2103.08609>`_
    *Simulations and models with and without lensed point sources to perform a time-delay cosmography analysis.*

* Gravitational lensing H0 tension from ultralight axion galactic cores `Blum & Teodori 2021 <https://arxiv.org/abs/2105.10873>`_
    *Investigating the detectability of a cored component with mock imaging modeling and comparison of kinematic modeling.*

* The Hubble constant from strongly lensed supernovae with standardizable magnifications `Birrer, Dhawan, Shajib 2021 <https://arxiv.org/abs/2107.12385>`_
    *Methodology and forecast to use standardizable magnifications to break the mass-sheet degeneracy and hierarchically measure H0.*

* AI-driven spatio-temporal engine for finding gravitationally lensed supernovae `Ramanah et al. 2021 <https://arxiv.org/abs/2107.12399>`_
    *Simulated images with time series of lensed supernovae.*

* Systematic errors induced by the elliptical power-law model in galaxy-galaxy strong lens modeling `Cao et al. 2021 <https://arxiv.org/abs/2110.14554>`_
    *Computing lensing quantities from mass maps.*

* TDCOSMO. VII. Boxyness/discyness in lensing galaxies : Detectability and impact on H0 `Van de Vyvere et al. 2021 <https://arxiv.org/abs/2112.03932>`_
    *Assessment of boxy and discy lens model on the inference of H0.*



Dark Matter substructure
------------------------

* Lensing substructure quantification in RXJ1131-1231: a 2 keV lower bound on dark matter thermal relic mass; `Birrer et al. 2017b <http://adsabs.harvard.edu/abs/2017JCAP...05..037B>`_
    *This paper quantifies the substructure content of a lens by a sub-clump scanning procedure and the application of Approximate Bayesian Computing.*

* Probing the nature of dark matter by forward modelling flux ratios in strong gravitational lenses; `Gilman et al. 2018 <https://ui.adsabs.harvard.edu/abs/2018MNRAS.481..819G>`_
    * *

* Probing dark matter structure down to 10**7 solar masses: flux ratio statistics in gravitational lenses with line-of-sight haloes; `Gilman et al. 2019a <https://ui.adsabs.harvard.edu/abs/2019MNRAS.487.5721G>`_
    * *

* Double dark matter vision: twice the number of compact-source lenses with narrow-line lensing and the WFC3 Grism; `Nierenberg et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019arXiv190806344N>`_
    * *

* Warm dark matter chills out: constraints on the halo mass function and the free-streaming length of dark matter with 8 quadruple-image strong gravitational lenses; `Gilman et al. 2019b <https://ui.adsabs.harvard.edu/abs/2019arXiv190806983G>`_
    * *

* Constraints on the mass-concentration relation of cold dark matter halos with 11 strong gravitational lenses; `Gilman et al. 2019c <https://ui.adsabs.harvard.edu/abs/2019arXiv190902573G>`_
    * *

* Circumventing Lens Modeling to Detect Dark Matter Substructure in Strong Lens Images with Convolutional Neural Networks; `Diaz Rivero & Dvorkin <https://ui.adsabs.harvard.edu/abs/2019arXiv191000015D>`_
    * *

* Dark Matter Subhalos, Strong Lensing and Machine Learning; `Varma, Fairbairn, Figueroa <https://arxiv.org/abs/2005.05353>`_
    * *

* Quantifying the Line-of-Sight Halo Contribution to the Dark Matter Convergence Power Spectrum from Strong Gravitational Lenses; `Sengul et al. 2020 <https://arxiv.org/abs/2006.07383>`_
    * *

* Detecting Subhalos in Strong Gravitational Lens Images with Image Segmentation; `Ostdiek et al. 2020a <https://arxiv.org/abs/2009.06663>`_
    * *

* Extracting the Subhalo Mass Function from Strong Lens Images with Image Segmentation; `Ostdiek et al. 2020b <https://arxiv.org/abs/2009.06639>`_
    * *

* Strong lensing signatures of self-interacting dark matter in low-mass halos; `Gilman et al. 2021 <https://arxiv.org/abs/2105.05259>`_
    * *

* Substructure Detection Reanalyzed: Dark Perturber shown to be a Line-of-Sight Halo; `Sengul et al. 2021 <https://arxiv.org/abs/2112.00749>`_
    *modeling a line-of-sight mini-halo*

* The primordial matter power spectrum on sub-galactic scales; `Gilman et al. 2021 <https://arxiv.org/abs/2112.03293>`_
    *rendering sub- and line-of-sight halos*




Lens searches
-------------

* Strong lens systems search in the Dark Energy Survey using Convolutional Neural Networks; `Rojas et al. 2021 <https://arxiv.org/abs/2109.00014>`_
    *simulating training sets for lens searches*

* On machine learning search for gravitational lenses; `Khachatryan 2021 <https://arxiv.org/abs/2104.01014>`_
    *simulating training sets for lens searches*

* DeepZipper: A Novel Deep Learning Architecture for Lensed Supernovae Identification; `Morgan et al. 2021b <hhttps://arxiv.org/abs/2112.01541>`_
    *simulating deeplenstronomy to simulate lensed supernovae data sets*


Galaxy formation and evolution
------------------------------

* Massive elliptical galaxies at z∼0.2 are well described by stars and a Navarro-Frenk-White dark matter halo; `Shajib et al. 2020a <https://arxiv.org/abs/2008.11724>`_
    *Automatized modeling of 23 SLACS lenses with dolphin, a lenstronomy wrapper*

* High-resolution imaging follow-up of doubly imaged quasars; `Shajib et al. 2020b <https://arxiv.org/abs/2011.01971>`_
    *Modeling of doubly lensed quasars from Keck Adaptive Optics data*

* The evolution of the size-mass relation at z=1-3 derived from the complete Hubble Frontier Fields data set; `Yang et al. 2020b <https://arxiv.org/abs/2011.10059>`_
    *reconstructing the intrinsic size-mass relation of strongly lensed sources in clusters*

* PS J1721+8842: A gravitationally lensed dual AGN system at redshift 2.37 with two radio components; `Mangat et al. 2021 <https://arxiv.org/abs/2109.03253>`_
    *Imaging modeling of a dual lensed AGN with point sources and extended surface brightness*

* RELICS: Small Lensed z≥5.5 Galaxies Selected as Potential Lyman Continuum Leakers; `Neufeld et al. 2021 <https://arxiv.org/abs/2111.14882>`_
    *size measurements of high-z lensed galaxies*

* The size-luminosity relation of lensed galaxies at z=6−9 in the Hubble Frontier Fields; `Yang et al. 2022 <https://arxiv.org/abs/2201.08858>`_
    *size measurements of high-z lensed galaxies*


Automatized Lens Modeling
-------------------------

* Is every strong lens model unhappy in its own way? Uniform modelling of a sample of 12 quadruply+ imaged quasars; `Shajib et al. 2018 <https://ui.adsabs.harvard.edu/abs/2019MNRAS.483.5649S>`_
    *This work presents a uniform modelling framework to model 13 quadruply lensed quasars in three HST bands.*

* Hierarchical Inference With Bayesian Neural Networks: An Application to Strong Gravitational Lensing; `Wagner-Carena et al. 2020 <https://arxiv.org/abs/2010.13787>`_
    *This work conducts hierarchical inference of strongly-lensed systems with Bayesian neural networks.*

* A search for galaxy-scale strong gravitational lenses in the Ultraviolet Near Infrared Optical Northern Survey (UNIONS); `Savary et al. 2021 <https://arxiv.org/abs/2110.11972>`_
    *Automated modeling of best candidates of ground based data.*





Quasar-host galaxy decomposition
--------------------------------


* The mass relations between supermassive black holes and their host galaxies at 1<z<2 with HST-WFC3; `Ding et al. 2019 <https://arxiv.org/abs/1910.11875>`_
    *Quasar host galaxy decomposition at high redshift on HST imaging and marginalization over PSF uncertainties.*

* Testing the Evolution of the Correlations between Supermassive Black Holes and their Host Galaxies using Eight Strongly Lensed Quasars; `Ding et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020arXiv200513550D/abstract>`_
    *Quasar host galaxy decomposition with lensed quasars.*

* A local baseline of the black hole mass scaling relations for active galaxies. IV. Correlations between MBH and host galaxy σ, stellar mass, and luminosity; `Bennert et al. 2021 <https://arxiv.org/abs/2101.10355>`_
    *Detailed measurement of galaxy morphology, decomposing in spheroid, disk and bar, and central AGN*

* The Sizes of Quasar Host Galaxies with the Hyper Suprime-Cam Subaru Strategic Program; `Li et al. 2021a <https://arxiv.org/abs/2105.06568>`_
    *Quasar-host decomposition of 5000 SDSS quasars*

* The eROSITA Final Equatorial-Depth Survey (eFEDS): A multiwavelength view of WISE mid-infrared galaxies/active galactic nuclei; `Toba et al. 2021 <https://arxiv.org/abs/2106.14527>`_
    *Quasar-host decomposition of HSC imaging*

* Synchronized Co-evolution between Supermassive Black Holes and Galaxies Over the Last Seven Billion Years as Revealed by the Hyper Suprime-Cam; `Li et al. 2021b <https://arxiv.org/abs/2109.02751>`_
    *Quasar-host decomposition of SDSS quasars with HSC data*






Lensing of Gravitational Waves
------------------------------
* lensingGW: a Python package for lensing of gravitational waves; `Pagano et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020arXiv200612879P/abstract>`_
    *A Python package designed to handle both strong and microlensing of compact binaries and the related gravitational-wave signals.*

* Localizing merging black holes with sub-arcsecond precision using gravitational-wave lensing; `Hannuksela et al. 2020 <https://arxiv.org/abs/2004.13811v3>`_
    *solving the lens equation with lenstronomy using lensingGW*

* Lensing magnification: gravitational wave from coalescing stellar-mass binary black holes; `Shan & Hu 2020 <https://arxiv.org/abs/2012.08381>`_
    *lensing magnificatoin calculations*

* Identifying Type-II Strongly-Lensed Gravitational-Wave Images in Third-Generation Gravitational-Wave Detectors; `Y. Wang et al. 2021 <https://arxiv.org/abs/2101.08264>`_
    *solving the lens equation*

* Beyond the detector horizon: Forecasting gravitational-wave strong lensing; `Renske et al. 2021 <https://arxiv.org/abs/2106.06303>`_
    *computing image positions, time delays and magnifications for gravitational wave forecasting*




Theory papers
-------------

* Line-of-sight effects in strong lensing: putting theory into practice; `Birrer et al. 2017a <http://adsabs.harvard.edu/abs/2017JCAP...04..049B>`_
    *This paper formulates an effective parameterization of line-of-sight structure for strong gravitational lens modelling and applies this technique to an Einstein ring in the COSMOS field*

* Cosmic Shear with Einstein Rings; `Birrer et al. 2018a <http://adsabs.harvard.edu/abs/2018ApJ...852L..14B>`_
    *Forecast paper to measure cosmic shear with Einstein ring lenses. The forecast is made based on lenstronomy simulations.*

* Unified lensing and kinematic analysis for any elliptical mass profile; `Shajib 2019 <https://ui.adsabs.harvard.edu/abs/2019MNRAS.488.1387S>`_
    *Provides a methodology to generalize the multi-Gaussian expansion to general elliptical mass and light profiles*

* Gravitational lensing formalism in a curved arc basis: A continuous description of observables and degeneracies from the weak to the strong lensing regime; `Birrer 2021 <https://arxiv.org/abs/2104.09522>`_
    *Lensing formalism with curved arc distortion formalism. Link to code repository `here <https://github.com/sibirrer/curved_arcs>`_.*





Simulation products
-------------------

* The LSST DESC DC2 Simulated Sky Survey; `LSST Dark Energy Science Collaboration et al. 2020 <https://arxiv.org/abs/2010.05926v1>`_
    *Strong lensing simulations produced by SLSprinkler utilizing lenstronomy functionalities*

* The impact of mass map truncation on strong lensing simulations; `Van de Vyvere et al. 2020 <https://arxiv.org/abs/2010.13650>`_
    *Uses numerical integration to compute lensing quantities from projected mass maps from simulations.*



Large scale structure
---------------------

* Combining strong and weak lensingestimates in the Cosmos field; `Kuhn et al. 2020 <https://arxiv.org/abs/2010.08680>`_
    *inferring cosmic shear with three strong lenses in the COSMOS field*




Others
------

* Predicting future astronomical events using deep learning; `Singh et al. <https://arxiv.org/abs/2012.15476>`_
    *simulating strongly lensed galaxy merger pairs in time sequence*
====================================================
lenstronomy - gravitational lensing software package
====================================================


.. image:: https://badge.fury.io/py/lenstronomy.png
    :target: https://badge.fury.io/py/lenstronomy

.. image:: https://github.com/sibirrer/lenstronomy/workflows/Tests/badge.svg
    :target: https://github.com/sibirrer/lenstronomy/actions

.. image:: https://readthedocs.org/projects/lenstronomy/badge/?version=latest
        :target: http://lenstronomy.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://coveralls.io/repos/github/sibirrer/lenstronomy/badge.svg?branch=main
        :target: https://coveralls.io/github/sibirrer/lenstronomy?branch=main

.. image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat
    :target: https://github.com/sibirrer/lenstronomy/blob/main/LICENSE

.. image:: https://img.shields.io/badge/arXiv-1803.09746%20-yellowgreen.svg
    :target: https://arxiv.org/abs/1803.09746

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
        :target: http://www.astropy.org
        :alt: Powered by Astropy Badge

.. image:: https://joss.theoj.org/papers/6a562375312c9a9e4466912a16f27589/status.svg
    :target: https://joss.theoj.org/papers/6a562375312c9a9e4466912a16f27589

.. image:: https://raw.githubusercontent.com/sibirrer/lenstronomy/main/docs/figures/readme_fig.png
    :target: https://raw.githubusercontent.com/sibirrer/lenstronomy/main/docs/figures/readme_fig.png


``lenstronomy`` is a multi-purpose package to model strong gravitational lenses. The software package is presented in
`Birrer & Amara 2018 <https://arxiv.org/abs/1803.09746v1>`_ and `Birrer et al. 2021 <https://joss.theoj.org/papers/10.21105/joss.03283>`_ , and is based on `Birrer et al 2015 <http://adsabs.harvard.edu/abs/2015ApJ...813..102B>`_.
``lenstronomy`` finds application for time-delay cosmography and measuring
the expansion rate of the Universe, for quantifying lensing substructure to infer dark matter properties, morphological quantification of galaxies,
quasar-host galaxy decomposition and much more.
A (incomplete) list of publications making use of lenstronomy can be found `at this link <https://github.com/sibirrer/lenstronomy/blob/main/PUBLISHED.rst>`_.


The development is coordinated on `GitHub <https://github.com/sibirrer/lenstronomy>`_ and contributions are welcome.
The documentation of ``lenstronomy`` is available at `readthedocs.org <http://lenstronomy.readthedocs.org/>`_ and
the package is distributed over `PyPI <https://pypi.python.org/pypi/lenstronomy>`_.
``lenstronomy`` is an `affiliated package <https://www.astropy.org/affiliated/>`_ of `astropy <https://www.astropy.org/>`_.



Installation
------------

.. code-block:: bash

    $ pip install lenstronomy --user


Specific instructions for settings and installation requirements for special cases that can provide speed-ups,
we refer to the `documentation <https://lenstronomy.readthedocs.io/en/latest/installation.html>`_ page.



Getting started
---------------

The `starting guide jupyter notebook <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/starting_guide.ipynb>`_
leads through the main modules and design features of ``lenstronomy``. The modular design of ``lenstronomy`` allows the
user to directly access a lot of tools and each module can also be used as stand-alone packages.

If you are new to gravitational lensing, check out the `mini lecture series <https://github.com/sibirrer/strong_lensing_lectures>`_ giving an introduction to gravitational lensing
with interactive Jupyter notebooks in the cloud.



Example notebooks
-----------------

We have made an extension module available at `https://github.com/sibirrer/lenstronomy_extensions <https://github.com/sibirrer/lenstronomy_extensions>`_.
You can find simple examle notebooks for various cases. The latest versions of the notebooks should be compatible with the recent pip version of lenstronomy.

* `Units, coordinate system and parameter definitions in lenstronomy <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/units_coordinates_parameters.ipynb>`_
* `FITS handling and extracting needed information from the data prior to modeling <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/fits_handling.ipynb>`_
* `Modeling a simple Einstein ring <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/simple_ring.ipynb>`_
* `Quadrupoly lensed quasar modelling <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/quad_model.ipynb>`_
* `Double lensed quasar modelling <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/double_model.ipynb>`_
* `Time-delay cosmography <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/time-delay%20cosmography.ipynb>`_
* `Source reconstruction and deconvolution with Shapelets <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/shapelet_source_modelling.ipynb>`_
* `Solving the lens equation <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/lens_equation.ipynb>`_
* `Multi-band fitting <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/multi_band_fitting.ipynb>`_
* `Measuring cosmic shear with Einstein rings <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/EinsteinRingShear_simulations.ipynb>`_
* `Fitting of galaxy light profiles, like e.g. GALFIT <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/galfitting.ipynb>`_
* `Quasar-host galaxy decomposition <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/quasar-host%20decomposition.ipynb>`_
* `Hiding and seeking a single subclump <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/substructure_challenge_simple.ipynb>`_
* `Mock generation of realistic images with substructure in the lens <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/substructure_challenge_mock_production.ipynb>`_
* `Mock simulation API with multi color models <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/simulation_api.ipynb>`_
* `Catalogue data modeling of image positions, flux ratios and time delays <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/catalogue%20modelling.ipynb>`_
* `Example of numerical ray-tracing and convolution options <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/lenstronomy_numerics.ipynb>`_
* `Simulated lenses with populations generated by SkyPy <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/skypy_lenstronomy.ipynb>`_



Affiliated packages
-------------------
Multiple affiliated packages that make use of lenstronomy can be found `here <https://lenstronomy.readthedocs.io/en/latest/affiliatedpackages.html>`_
(not complete) and further packages are under development by the community.


Mailing list and Slack channel
------------------------------

You can join the ``lenstronomy`` mailing list by signing up on the
`google groups page <https://groups.google.com/forum/#!forum/lenstronomy>`_.


The email list is meant to provide a communication platform between users and developers. You can ask questions,
and suggest new features. New releases will be announced via this mailing list.

We also have a `Slack channel <https://lenstronomers.slack.com>`_ for the community.
Please send me an `email <sibirrer@gmail.com>`_ such that I can add you to the channel.


If you encounter errors or problems with ``lenstronomy``, please let us know!



Contribution
------------
Check out the `contributing page <https://lenstronomy.readthedocs.io/en/latest/contributing.html>`_
and become an author of ``lenstronomy``! A big shutout to the current `list of contributors and developers <https://lenstronomy.readthedocs.io/en/latest/authors.html>`_!




Shapelet reconstruction demonstration movies
--------------------------------------------

We provide some examples where a real galaxy has been lensed and then been reconstructed by a shapelet basis set.

* `HST quality data with perfect knowledge of the lens model <http://www.astro.ucla.edu/~sibirrer/video/true_reconstruct.mp4>`_
* `HST quality with a clump hidden in the data <http://www.astro.ucla.edu/~sibirrer/video/clump_reconstruct.mp4>`_
* `Extremely large telescope quality data with a clump hidden in the data <http://www.astro.ucla.edu/~sibirrer/video/TMT_high_res_clump_reconstruct.mp4>`_



Attribution
-----------
The design concept of ``lenstronomy`` is reported by `Birrer & Amara 2018 <https://arxiv.org/abs/1803.09746v1>`_.
The current JOSS software publication is presented by `Birrer et al. 2021 <https://joss.theoj.org/papers/10.21105/joss.03283>`_.
Please cite these two publications when you use lenstronomy in a publication and link to `https://github.com/sibirrer/lenstronomy <https://github.com/sibirrer/lenstronomy>`_.
Please also cite `Birrer et al 2015 <http://adsabs.harvard.edu/abs/2015ApJ...813..102B>`_
when you make use of the ``lenstronomy`` work-flow or the Shapelet source reconstruction and make sure to cite also
the relevant work that was implemented in ``lenstronomy``, as described in the release paper and the documentation.
Don't hesitate to reach out to the developers if you have questions!
lenstronomy.Util package
========================

Submodules
----------

lenstronomy.Util.analysis\_util module
--------------------------------------

.. automodule:: lenstronomy.Util.analysis_util
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Util.class\_creator module
--------------------------------------

.. automodule:: lenstronomy.Util.class_creator
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Util.constants module
---------------------------------

.. automodule:: lenstronomy.Util.constants
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Util.correlation module
-----------------------------------

.. automodule:: lenstronomy.Util.correlation
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Util.data\_util module
----------------------------------

.. automodule:: lenstronomy.Util.data_util
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Util.derivative\_util module
----------------------------------------

.. automodule:: lenstronomy.Util.derivative_util
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Util.image\_util module
-----------------------------------

.. automodule:: lenstronomy.Util.image_util
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Util.kernel\_util module
------------------------------------

.. automodule:: lenstronomy.Util.kernel_util
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Util.mask\_util module
----------------------------------

.. automodule:: lenstronomy.Util.mask_util
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Util.multi\_gauss\_expansion module
-----------------------------------------------

.. automodule:: lenstronomy.Util.multi_gauss_expansion
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Util.numba\_util module
-----------------------------------

.. automodule:: lenstronomy.Util.numba_util
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Util.param\_util module
-----------------------------------

.. automodule:: lenstronomy.Util.param_util
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Util.prob\_density module
-------------------------------------

.. automodule:: lenstronomy.Util.prob_density
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Util.sampling\_util module
--------------------------------------

.. automodule:: lenstronomy.Util.sampling_util
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Util.simulation\_util module
----------------------------------------

.. automodule:: lenstronomy.Util.simulation_util
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Util.util module
----------------------------

.. automodule:: lenstronomy.Util.util
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.Util
    :members:
    :undoc-members:
    :show-inheritance:
===========================
Contributing to lenstronomy
===========================

.. include:: ../CONTRIBUTING.mdlenstronomy.Conf package
========================

Submodules
----------

lenstronomy.Conf.config\_loader module
--------------------------------------

.. automodule:: lenstronomy.Conf.config_loader
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.Conf
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy.ImSim package
=========================

Subpackages
-----------

.. toctree::

    lenstronomy.ImSim.MultiBand
    lenstronomy.ImSim.Numerics

Submodules
----------

lenstronomy.ImSim.de\_lens module
---------------------------------

.. automodule:: lenstronomy.ImSim.de_lens
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.ImSim.image2source\_mapping module
----------------------------------------------

.. automodule:: lenstronomy.ImSim.image2source_mapping
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.ImSim.image\_linear\_solve module
---------------------------------------------

.. automodule:: lenstronomy.ImSim.image_linear_solve
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.ImSim.image\_model module
-------------------------------------

.. automodule:: lenstronomy.ImSim.image_model
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.ImSim
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy.LensModel.Util package
==================================

Submodules
----------

lenstronomy.LensModel.Util.epl\_util module
-------------------------------------------

.. automodule:: lenstronomy.LensModel.Util.epl_util
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: lenstronomy.LensModel.Util
   :members:
   :undoc-members:
   :show-inheritance:
lenstronomy.Sampling package
============================

Subpackages
-----------

.. toctree::

    lenstronomy.Sampling.Likelihoods
    lenstronomy.Sampling.Pool
    lenstronomy.Sampling.Samplers

Submodules
----------

lenstronomy.Sampling.likelihood module
--------------------------------------

.. automodule:: lenstronomy.Sampling.likelihood
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Sampling.parameters module
--------------------------------------

.. automodule:: lenstronomy.Sampling.parameters
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Sampling.sampler module
-----------------------------------

.. automodule:: lenstronomy.Sampling.sampler
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Sampling.special\_param module
------------------------------------------

.. automodule:: lenstronomy.Sampling.special_param
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.Sampling
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy\.LensModel\.Solver package
======================================

Submodules
----------

lenstronomy\.LensModel\.Solver\.lens\_equation\_solver module
-------------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Solver.lens_equation_solver
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy\.LensModel\.Solver\.solver module
---------------------------------------------

.. automodule:: lenstronomy.LensModel.Solver.solver
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy\.LensModel\.Solver\.solver2point module
---------------------------------------------------

.. automodule:: lenstronomy.LensModel.Solver.solver2point
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy\.LensModel\.Solver\.solver4point module
---------------------------------------------------

.. automodule:: lenstronomy.LensModel.Solver.solver4point
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.LensModel.Solver
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy.GalKin package
==========================

Submodules
----------

lenstronomy.GalKin.analytic\_kinematics module
----------------------------------------------

.. automodule:: lenstronomy.GalKin.analytic_kinematics
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.GalKin.anisotropy module
------------------------------------

.. automodule:: lenstronomy.GalKin.anisotropy
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.GalKin.aperture module
----------------------------------

.. automodule:: lenstronomy.GalKin.aperture
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.GalKin.aperture\_types module
-----------------------------------------

.. automodule:: lenstronomy.GalKin.aperture_types
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.GalKin.cosmo module
-------------------------------

.. automodule:: lenstronomy.GalKin.cosmo
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.GalKin.galkin module
--------------------------------

.. automodule:: lenstronomy.GalKin.galkin
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.GalKin.galkin\_model module
---------------------------------------

.. automodule:: lenstronomy.GalKin.galkin_model
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.GalKin.light\_profile module
----------------------------------------

.. automodule:: lenstronomy.GalKin.light_profile
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.GalKin.numeric\_kinematics module
---------------------------------------------

.. automodule:: lenstronomy.GalKin.numeric_kinematics
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.GalKin.observation module
-------------------------------------

.. automodule:: lenstronomy.GalKin.observation
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.GalKin.psf module
-----------------------------

.. automodule:: lenstronomy.GalKin.psf
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.GalKin.velocity\_util module
----------------------------------------

.. automodule:: lenstronomy.GalKin.velocity_util
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.GalKin
    :members:
    :undoc-members:
    :show-inheritance:
.. include:: ../HISTORY.rstlenstronomy\.Data package
=========================

Submodules
----------

lenstronomy\.Data\.coord\_transforms module
-------------------------------------------

.. automodule:: lenstronomy.Data.coord_transforms
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy\.Data\.imaging\_data module
---------------------------------------

.. automodule:: lenstronomy.Data.imaging_data
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy\.Data\.psf module
-----------------------------

.. automodule:: lenstronomy.Data.psf
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.Data
    :members:
    :undoc-members:
    :show-inheritance:
========
Usage
========

To use lenstronomy in a project::

	import lenstronomy

Getting started
---------------

The `starting guide jupyter notebook <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/starting_guide.ipynb>`_
leads through the main modules and design features of **lenstronomy**. The modular design of **lenstronomy** allows the
user to directly access a lot of tools and each module can also be used as stand-alone packages.


Example notebooks
-----------------

We have made an extension module available at `https://github.com/sibirrer/lenstronomy_extensions <https://github.com/sibirrer/lenstronomy_extensions>`_.
You can find simple examle notebooks for various cases. The latest versions of the notebooks should be compatible with the recent pip version of lenstronomy.

* `Units, coordinate system and parameter definitions in lenstronomy <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/units_coordinates_parameters.ipynb>`_
* `FITS handling and extracting needed information from the data prior to modeling <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/fits_handling.ipynb>`_
* `Modeling a simple Einstein ring <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/simple_ring.ipynb>`_
* `Quadrupoly lensed quasar modelling <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/quad_model.ipynb>`_
* `Double lensed quasar modelling <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/double_model.ipynb>`_
* `Time-delay cosmography <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/time-delay%20cosmography.ipynb>`_
* `Source reconstruction and deconvolution with Shapelets <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/shapelet_source_modelling.ipynb>`_
* `Solving the lens equation <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/lens_equation.ipynb>`_
* `Multi-band fitting <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/multi_band_fitting.ipynb>`_
* `Measuring cosmic shear with Einstein rings <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/EinsteinRingShear_simulations.ipynb>`_
* `Fitting of galaxy light profiles, like e.g. GALFIT <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/galfitting.ipynb>`_
* `Quasar-host galaxy decomposition <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/quasar-host%20decomposition.ipynb>`_
* `Hiding and seeking a single subclump <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/substructure_challenge_simple.ipynb>`_
* `Mock generation of realistic images with substructure in the lens <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/substructure_challenge_mock_production.ipynb>`_
* `Mock simulation API with multi color models <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/simulation_api.ipynb>`_
* `Catalogue data modeling of image positions, flux ratios and time delays <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/catalogue%20modelling.ipynb>`_
* `Example of numerical ray-tracing and convolution options <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/lenstronomy_numerics.ipynb>`_
* `Simulated lenses with populations generated by SkyPy <https://github.com/sibirrer/lenstronomy_extensions/blob/main/lenstronomy_extensions/Notebooks/skypy_lenstronomy.ipynb>`_

.. include:: ../AFFILIATEDPACKAGES.rstlenstronomy.LensModel.MultiPlane package
========================================

Submodules
----------

lenstronomy.LensModel.MultiPlane.multi\_plane module
----------------------------------------------------

.. automodule:: lenstronomy.LensModel.MultiPlane.multi_plane
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.MultiPlane.multi\_plane\_base module
----------------------------------------------------------

.. automodule:: lenstronomy.LensModel.MultiPlane.multi_plane_base
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.LensModel.MultiPlane
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy.LensModel.Profiles package
======================================

Submodules
----------

lenstronomy.LensModel.Profiles.arc\_perturbations module
--------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.arc_perturbations
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.base\_profile module
---------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.base_profile
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.chameleon module
-----------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.chameleon
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.cnfw module
------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.cnfw
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.cnfw\_ellipse module
---------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.cnfw_ellipse
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.const\_mag module
------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.const_mag
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.constant\_shift module
-----------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.constant_shift
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.convergence module
-------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.convergence
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.coreBurkert module
-------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.coreBurkert
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.cored\_density module
----------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.cored_density
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.cored\_density\_2 module
-------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.cored_density_2
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.cored\_density\_exp module
---------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.cored_density_exp
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.cored\_density\_mst module
---------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.cored_density_mst
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.curved\_arc\_const module
--------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.curved_arc_const
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.curved\_arc\_sis\_mst module
-----------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.curved_arc_sis_mst
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.curved\_arc\_spp module
------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.curved_arc_spp
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.curved\_arc\_spt module
------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.curved_arc_spt
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.curved\_arc\_tan\_diff module
------------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.curved_arc_tan_diff
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.dipole module
--------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.dipole
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.elliptical\_density\_slice module
----------------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.elliptical_density_slice
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.epl module
-----------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.epl
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.epl\_numba module
------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.epl_numba
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.flexion module
---------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.flexion
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.flexionfg module
-----------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.flexionfg
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.gauss\_decomposition module
----------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.gauss_decomposition
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.gaussian\_ellipse\_kappa module
--------------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.gaussian_ellipse_kappa
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.gaussian\_ellipse\_potential module
------------------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.gaussian_ellipse_potential
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.gaussian\_kappa module
-----------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.gaussian_kappa
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.gaussian\_potential module
---------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.gaussian_potential
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.hernquist module
-----------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.hernquist
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.hernquist\_ellipse module
--------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.hernquist_ellipse
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.hessian module
---------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.hessian
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.interpol module
----------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.interpol
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.multi\_gaussian\_kappa module
------------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.multi_gaussian_kappa
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.multipole module
-----------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.multipole
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.nfw module
-----------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.nfw
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.nfw\_ellipse module
--------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.nfw_ellipse
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.nfw\_mass\_concentration module
--------------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.nfw_mass_concentration
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.nfw\_vir\_trunc module
-----------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.nfw_vir_trunc
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.nie module
-----------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.nie
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.nie\_potential module
----------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.nie_potential
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.numerical\_deflections module
------------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.numerical_deflections
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.p\_jaffe module
----------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.p_jaffe
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.p\_jaffe\_ellipse module
-------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.p_jaffe_ellipse
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.pemd module
------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.pemd
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.point\_mass module
-------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.point_mass
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.sersic module
--------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.sersic
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.sersic\_ellipse\_kappa module
------------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.sersic_ellipse_kappa
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.sersic\_ellipse\_potential module
----------------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.sersic_ellipse_potential
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.sersic\_utils module
---------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.sersic_utils
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.shapelet\_pot\_cartesian module
--------------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.shapelet_pot_cartesian
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.shapelet\_pot\_polar module
----------------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.shapelet_pot_polar
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.shear module
-------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.shear
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.sie module
-----------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.sie
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.sis module
-----------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.sis
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.sis\_truncate module
---------------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.sis_truncate
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.spemd module
-------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.spemd
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.spep module
------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.spep
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.splcore module
---------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.splcore
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.spp module
-----------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.spp
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.tnfw module
------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.tnfw
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.Profiles.uldm module
------------------------------------------

.. automodule:: lenstronomy.LensModel.Profiles.uldm
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.LensModel.Profiles
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy.LensModel package
=============================

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   lenstronomy.LensModel.LightConeSim
   lenstronomy.LensModel.MultiPlane
   lenstronomy.LensModel.Profiles
   lenstronomy.LensModel.QuadOptimizer
   lenstronomy.LensModel.Solver
   lenstronomy.LensModel.Util

Submodules
----------

lenstronomy.LensModel.convergence\_integrals module
---------------------------------------------------

.. automodule:: lenstronomy.LensModel.convergence_integrals
   :members:
   :undoc-members:
   :show-inheritance:

lenstronomy.LensModel.lens\_model module
----------------------------------------

.. automodule:: lenstronomy.LensModel.lens_model
   :members:
   :undoc-members:
   :show-inheritance:

lenstronomy.LensModel.lens\_model\_extensions module
----------------------------------------------------

.. automodule:: lenstronomy.LensModel.lens_model_extensions
   :members:
   :undoc-members:
   :show-inheritance:

lenstronomy.LensModel.lens\_param module
----------------------------------------

.. automodule:: lenstronomy.LensModel.lens_param
   :members:
   :undoc-members:
   :show-inheritance:

lenstronomy.LensModel.profile\_integrals module
-----------------------------------------------

.. automodule:: lenstronomy.LensModel.profile_integrals
   :members:
   :undoc-members:
   :show-inheritance:

lenstronomy.LensModel.profile\_list\_base module
------------------------------------------------

.. automodule:: lenstronomy.LensModel.profile_list_base
   :members:
   :undoc-members:
   :show-inheritance:

lenstronomy.LensModel.single\_plane module
------------------------------------------

.. automodule:: lenstronomy.LensModel.single_plane
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: lenstronomy.LensModel
   :members:
   :undoc-members:
   :show-inheritance:
lenstronomy.PointSource package
===============================

Subpackages
-----------

.. toctree::

    lenstronomy.PointSource.Types

Submodules
----------

lenstronomy.PointSource.point\_source module
--------------------------------------------

.. automodule:: lenstronomy.PointSource.point_source
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.PointSource.point\_source\_cached module
----------------------------------------------------

.. automodule:: lenstronomy.PointSource.point_source_cached
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.PointSource.point\_source\_param module
---------------------------------------------------

.. automodule:: lenstronomy.PointSource.point_source_param
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.PointSource
    :members:
    :undoc-members:
    :show-inheritance:
.. include:: ../AUTHORS.rstlenstronomy.Analysis package
============================

Submodules
----------

lenstronomy.Analysis.kinematics\_api module
-------------------------------------------

.. automodule:: lenstronomy.Analysis.kinematics_api
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Analysis.lens\_profile module
-----------------------------------------

.. automodule:: lenstronomy.Analysis.lens_profile
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Analysis.light2mass module
--------------------------------------

.. automodule:: lenstronomy.Analysis.light2mass
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Analysis.light\_profile module
------------------------------------------

.. automodule:: lenstronomy.Analysis.light_profile
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Analysis.td\_cosmography module
-------------------------------------------

.. automodule:: lenstronomy.Analysis.td_cosmography
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.Analysis
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy.LensModel.LightConeSim package
==========================================

Submodules
----------

lenstronomy.LensModel.LightConeSim.light\_cone module
-----------------------------------------------------

.. automodule:: lenstronomy.LensModel.LightConeSim.light_cone
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.LensModel.LightConeSim
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy.LightModel.Profiles package
=======================================

Submodules
----------

lenstronomy.LightModel.Profiles.chameleon module
------------------------------------------------

.. automodule:: lenstronomy.LightModel.Profiles.chameleon
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LightModel.Profiles.ellipsoid module
------------------------------------------------

.. automodule:: lenstronomy.LightModel.Profiles.ellipsoid
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LightModel.Profiles.gaussian module
-----------------------------------------------

.. automodule:: lenstronomy.LightModel.Profiles.gaussian
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LightModel.Profiles.hernquist module
------------------------------------------------

.. automodule:: lenstronomy.LightModel.Profiles.hernquist
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LightModel.Profiles.interpolation module
----------------------------------------------------

.. automodule:: lenstronomy.LightModel.Profiles.interpolation
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LightModel.Profiles.moffat module
---------------------------------------------

.. automodule:: lenstronomy.LightModel.Profiles.moffat
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LightModel.Profiles.nie module
------------------------------------------

.. automodule:: lenstronomy.LightModel.Profiles.nie
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LightModel.Profiles.p\_jaffe module
-----------------------------------------------

.. automodule:: lenstronomy.LightModel.Profiles.p_jaffe
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LightModel.Profiles.power\_law module
-------------------------------------------------

.. automodule:: lenstronomy.LightModel.Profiles.power_law
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LightModel.Profiles.sersic module
---------------------------------------------

.. automodule:: lenstronomy.LightModel.Profiles.sersic
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LightModel.Profiles.shapelets module
------------------------------------------------

.. automodule:: lenstronomy.LightModel.Profiles.shapelets
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LightModel.Profiles.shapelets\_polar module
-------------------------------------------------------

.. automodule:: lenstronomy.LightModel.Profiles.shapelets_polar
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LightModel.Profiles.uniform module
----------------------------------------------

.. automodule:: lenstronomy.LightModel.Profiles.uniform
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.LightModel.Profiles
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy.Cosmo package
=========================

Submodules
----------

lenstronomy.Cosmo.background module
-----------------------------------

.. automodule:: lenstronomy.Cosmo.background
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Cosmo.cosmo\_solver module
--------------------------------------

.. automodule:: lenstronomy.Cosmo.cosmo_solver
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Cosmo.kde\_likelihood module
----------------------------------------

.. automodule:: lenstronomy.Cosmo.kde_likelihood
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Cosmo.lcdm module
-----------------------------

.. automodule:: lenstronomy.Cosmo.lcdm
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Cosmo.lens\_cosmo module
------------------------------------

.. automodule:: lenstronomy.Cosmo.lens_cosmo
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Cosmo.nfw\_param module
-----------------------------------

.. automodule:: lenstronomy.Cosmo.nfw_param
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.Cosmo
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy.Sampling.Samplers package
=====================================

Submodules
----------

lenstronomy.Sampling.Samplers.base\_nested\_sampler module
----------------------------------------------------------

.. automodule:: lenstronomy.Sampling.Samplers.base_nested_sampler
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Sampling.Samplers.dynesty\_sampler module
-----------------------------------------------------

.. automodule:: lenstronomy.Sampling.Samplers.dynesty_sampler
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Sampling.Samplers.multinest\_sampler module
-------------------------------------------------------

.. automodule:: lenstronomy.Sampling.Samplers.multinest_sampler
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Sampling.Samplers.polychord\_sampler module
-------------------------------------------------------

.. automodule:: lenstronomy.Sampling.Samplers.polychord_sampler
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.Sampling.Samplers
    :members:
    :undoc-members:
    :show-inheritance:
============
Installation
============

At the command line with pip::

    $ pip install lenstronomy

Or, if you have virtualenvwrapper installed::

    $ mkvirtualenv lenstronomy
    $ pip install lenstronomy

You can also clone the github repository for development purposes.


Requirements
------------

Make sure the standard python libraries as specified in the `requirements <https://github.com/sibirrer/lenstronomy/blob/main/requirements.txt>`_.
The standard usage does not require all libraries to be installed, in particular the different posterior samplers are only required when being used.

In the following, a few specific cases are mentioned that may require special attention in the installation and settings, in particular when it comes
to MPI and HPC applications.


MPI
---
MPI support is provided for several sampling techniques for parallel computing. A specific version of the library schwimmbad is required
for the correct support of the moving of the likelihood elements from one processor to another with MPI. Pay attention ot the
`requirements <https://github.com/sibirrer/lenstronomy/blob/main/requirements.txt>`_.


NUMBA
-----
Just-in-time (jit) compilation with numba can provide significant speed-up for certain calculations.
There are specific settings for the settings provided as per default, but these may need to be adopted when running on a HPC cluster.
You can define your own configuration file in your $XDG_CONFIG_HOME/lenstronomy/config.yaml file. E.g. (check your system for the path)::

    $ ~/.conf/lenstronomy/config.yaml

following the format of the default configuration which is `here <https://github.com/sibirrer/lenstronomy/blob/main/lenstronomy/Conf/conf_default.yaml>`_.


FASTELL
-------
The fastell4py package, originally from Barkana (fastell), is required to run the PEMD (power-law elliptical mass distribution) lens model
and can be cloned from: `https://github.com/sibirrer/fastell4py <https://github.com/sibirrer/fastell4py>`_ (needs a fortran compiler).
We recommend using the EPL model as it is a pure python version of the same profile.

.. code-block:: bash

    $ sudo apt-get install gfortran
    $ git clone https://github.com/sibirrer/fastell4py.git <desired location>
    $ cd <desired location>
    $ python setup.py install --user


Check installation by running tests
-----------------------------------

You can check your installation with pytest::

    $ cd <lenstronomy_repo>
    $ py.test

Or you can run a partial test with::

    $ cd <lenstronomy_repo>
    $ py.test/test/test_LensModel/

You can also run the tests with tox in a virtual environment with::

    $ cd <lenstronomy_repo>
    $ tox

Note: tox might have trouble with the PyMultiNest installation and the cmake part of it.
lenstronomy package
===================

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   lenstronomy.Analysis
   lenstronomy.Conf
   lenstronomy.Cosmo
   lenstronomy.Data
   lenstronomy.GalKin
   lenstronomy.ImSim
   lenstronomy.LensModel
   lenstronomy.LightModel
   lenstronomy.Plots
   lenstronomy.PointSource
   lenstronomy.Sampling
   lenstronomy.SimulationAPI
   lenstronomy.Util
   lenstronomy.Workflow

Module contents
---------------

.. automodule:: lenstronomy
   :members:
   :undoc-members:
   :show-inheritance:
lenstronomy.Sampling.Likelihoods package
========================================

Submodules
----------

lenstronomy.Sampling.Likelihoods.image\_likelihood module
---------------------------------------------------------

.. automodule:: lenstronomy.Sampling.Likelihoods.image_likelihood
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Sampling.Likelihoods.position\_likelihood module
------------------------------------------------------------

.. automodule:: lenstronomy.Sampling.Likelihoods.position_likelihood
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Sampling.Likelihoods.prior\_likelihood module
---------------------------------------------------------

.. automodule:: lenstronomy.Sampling.Likelihoods.prior_likelihood
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Sampling.Likelihoods.time\_delay\_likelihood module
---------------------------------------------------------------

.. automodule:: lenstronomy.Sampling.Likelihoods.time_delay_likelihood
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.Sampling.Likelihoods
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy.ImSim.MultiBand package
===================================

Submodules
----------

lenstronomy.ImSim.MultiBand.joint\_linear module
------------------------------------------------

.. automodule:: lenstronomy.ImSim.MultiBand.joint_linear
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.ImSim.MultiBand.multi\_data\_base module
----------------------------------------------------

.. automodule:: lenstronomy.ImSim.MultiBand.multi_data_base
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.ImSim.MultiBand.multi\_linear module
------------------------------------------------

.. automodule:: lenstronomy.ImSim.MultiBand.multi_linear
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.ImSim.MultiBand.single\_band\_multi\_model module
-------------------------------------------------------------

.. automodule:: lenstronomy.ImSim.MultiBand.single_band_multi_model
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.ImSim.MultiBand
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy.Workflow package
============================

Submodules
----------

lenstronomy.Workflow.alignment\_matching module
-----------------------------------------------

.. automodule:: lenstronomy.Workflow.alignment_matching
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Workflow.fitting\_sequence module
---------------------------------------------

.. automodule:: lenstronomy.Workflow.fitting_sequence
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Workflow.psf\_fitting module
----------------------------------------

.. automodule:: lenstronomy.Workflow.psf_fitting
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Workflow.update\_manager module
-------------------------------------------

.. automodule:: lenstronomy.Workflow.update_manager
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.Workflow
    :members:
    :undoc-members:
    :show-inheritance:
.. include:: ../PUBLISHED.rstlenstronomy
===========

.. toctree::
   :maxdepth: 4

   lenstronomy
lenstronomy.LensModel.QuadOptimizer package
===========================================

Submodules
----------

lenstronomy.LensModel.QuadOptimizer.multi\_plane\_fast module
-------------------------------------------------------------

.. automodule:: lenstronomy.LensModel.QuadOptimizer.multi_plane_fast
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.QuadOptimizer.optimizer module
----------------------------------------------------

.. automodule:: lenstronomy.LensModel.QuadOptimizer.optimizer
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.LensModel.QuadOptimizer.param\_manager module
---------------------------------------------------------

.. automodule:: lenstronomy.LensModel.QuadOptimizer.param_manager
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.LensModel.QuadOptimizer
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy.SimulationAPI.ObservationConfig package
===================================================

Submodules
----------

lenstronomy.SimulationAPI.ObservationConfig.DES module
------------------------------------------------------

.. automodule:: lenstronomy.SimulationAPI.ObservationConfig.DES
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.SimulationAPI.ObservationConfig.Euclid module
---------------------------------------------------------

.. automodule:: lenstronomy.SimulationAPI.ObservationConfig.Euclid
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.SimulationAPI.ObservationConfig.HST module
------------------------------------------------------

.. automodule:: lenstronomy.SimulationAPI.ObservationConfig.HST
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.SimulationAPI.ObservationConfig.LSST module
-------------------------------------------------------

.. automodule:: lenstronomy.SimulationAPI.ObservationConfig.LSST
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.SimulationAPI.ObservationConfig
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy\.LightModel package
===============================

Subpackages
-----------

.. toctree::

    lenstronomy.LightModel.Profiles

Submodules
----------

lenstronomy\.LightModel\.light\_model module
--------------------------------------------

.. automodule:: lenstronomy.LightModel.light_model
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy\.LightModel\.light\_param module
--------------------------------------------

.. automodule:: lenstronomy.LightModel.light_param
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.LightModel
    :members:
    :undoc-members:
    :show-inheritance:
.. complexity documentation master file, created by
   sphinx-quickstart on Tue Jul  9 22:26:36 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. include:: ../README.rst

Contents:
=========

.. toctree::
   :maxdepth: 2

   installation
   usage
   lenstronomy
   contributing
   mailinglist
   authors
   published
   affiliatedpackages
   history
   modules

Feedback
========

If you have any suggestions or questions about **lenstronomy** feel free to email me
at sibirrer@gmail.com.

If you encounter any errors or problems with **lenstronomy**, please let me know!lenstronomy.Plots package
=========================

Submodules
----------

lenstronomy.Plots.chain\_plot module
------------------------------------

.. automodule:: lenstronomy.Plots.chain_plot
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Plots.lens\_plot module
-----------------------------------

.. automodule:: lenstronomy.Plots.lens_plot
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Plots.model\_band\_plot module
------------------------------------------

.. automodule:: lenstronomy.Plots.model_band_plot
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Plots.model\_plot module
------------------------------------

.. automodule:: lenstronomy.Plots.model_plot
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Plots.plot\_util module
-----------------------------------

.. automodule:: lenstronomy.Plots.plot_util
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.Plots
    :members:
    :undoc-members:
    :show-inheritance:
.. include:: ../MAILINGLIST.rstlenstronomy.ImSim.Numerics package
==================================

Submodules
----------

lenstronomy.ImSim.Numerics.adaptive\_numerics module
----------------------------------------------------

.. automodule:: lenstronomy.ImSim.Numerics.adaptive_numerics
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.ImSim.Numerics.convolution module
---------------------------------------------

.. automodule:: lenstronomy.ImSim.Numerics.convolution
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.ImSim.Numerics.grid module
--------------------------------------

.. automodule:: lenstronomy.ImSim.Numerics.grid
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.ImSim.Numerics.numba\_convolution module
----------------------------------------------------

.. automodule:: lenstronomy.ImSim.Numerics.numba_convolution
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.ImSim.Numerics.numerics module
------------------------------------------

.. automodule:: lenstronomy.ImSim.Numerics.numerics
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.ImSim.Numerics.partial\_image module
------------------------------------------------

.. automodule:: lenstronomy.ImSim.Numerics.partial_image
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.ImSim.Numerics.point\_source\_rendering module
----------------------------------------------------------

.. automodule:: lenstronomy.ImSim.Numerics.point_source_rendering
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.ImSim.Numerics
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy.Sampling.Pool package
=================================

Submodules
----------

lenstronomy.Sampling.Pool.multiprocessing module
------------------------------------------------

.. automodule:: lenstronomy.Sampling.Pool.multiprocessing
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.Sampling.Pool.pool module
-------------------------------------

.. automodule:: lenstronomy.Sampling.Pool.pool
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.Sampling.Pool
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy.PointSource.Types package
=====================================

Submodules
----------

lenstronomy.PointSource.Types.base\_ps module
---------------------------------------------

.. automodule:: lenstronomy.PointSource.Types.base_ps
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.PointSource.Types.lensed\_position module
-----------------------------------------------------

.. automodule:: lenstronomy.PointSource.Types.lensed_position
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.PointSource.Types.source\_position module
-----------------------------------------------------

.. automodule:: lenstronomy.PointSource.Types.source_position
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.PointSource.Types.unlensed module
---------------------------------------------

.. automodule:: lenstronomy.PointSource.Types.unlensed
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.PointSource.Types
    :members:
    :undoc-members:
    :show-inheritance:
lenstronomy.SimulationAPI package
=================================

Subpackages
-----------

.. toctree::

    lenstronomy.SimulationAPI.ObservationConfig

Submodules
----------

lenstronomy.SimulationAPI.data\_api module
------------------------------------------

.. automodule:: lenstronomy.SimulationAPI.data_api
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.SimulationAPI.model\_api module
-------------------------------------------

.. automodule:: lenstronomy.SimulationAPI.model_api
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.SimulationAPI.observation\_api module
-------------------------------------------------

.. automodule:: lenstronomy.SimulationAPI.observation_api
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.SimulationAPI.observation\_constructor module
---------------------------------------------------------

.. automodule:: lenstronomy.SimulationAPI.observation_constructor
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.SimulationAPI.point\_source\_variability module
-----------------------------------------------------------

.. automodule:: lenstronomy.SimulationAPI.point_source_variability
    :members:
    :undoc-members:
    :show-inheritance:

lenstronomy.SimulationAPI.sim\_api module
-----------------------------------------

.. automodule:: lenstronomy.SimulationAPI.sim_api
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: lenstronomy.SimulationAPI
    :members:
    :undoc-members:
    :show-inheritance:
