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
