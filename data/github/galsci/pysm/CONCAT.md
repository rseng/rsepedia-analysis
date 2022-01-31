---
title: 'The Python Sky Model 3 software'
tags:
  - cosmology
  - astronomy
  - python
authors:
 - name: Andrea Zonca
   orcid: 0000-0001-6841-1058
   affiliation: "1"
 - name: Ben Thorne
   orcid: 0000-0002-0457-0153
   affiliation: "2"
 - name: Nicoletta Krachmalnicoff
   affiliation: "3,4,5"
 - name: Julian Borrill
   affiliation: "6,7"
affiliations:
 - name: San Diego Supercomputer Center, University of California San Diego, San Diego, USA
   index: 1
 - name: Department of Physics, University of California Davis, One Shields Avenue, Davis, CA 95616, USA
   index: 2
 - name: SISSA, Via Bonomea 265, 34136 Trieste, Italy
   index: 3
 - name: INFN, Via Valerio 2, 34127 Trieste, Italy
   index: 4
 - name: IFPU, Via Beirut 2, 34014 Trieste, Italy
   index: 5
 - name: Computational Cosmology Center, Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA
   index: 6
 - name: Space Sciences Laboratory at University of California, 7 Gauss Way, Berkeley, CA 94720
   index: 7
date: 22 July 2021
bibliography: paper.bib
---

# Statement of Need

The Cosmic Microwave Background (CMB) radiation, emitted just 370 thousand years after the Big Bang, is a pristine probe of the Early Universe. After being emitted at high temperatures, the CMB was redshifted by the subsequent 13.8 billion years of cosmic expansion, such that it is brightest at microwave frequencies today.
However, our own Milky Way galaxy also emits in the microwave portion of the spectrum, obscuring our view of the CMB.  Examples of this emission are thermal radiation by interstellar dust grains and synchrotron emission by relativistic electrons spiraling in magnetic fields.
Cosmologists need to create synthetic maps of the CMB and of the galactic emission based on available data and on physical models that extrapolate observations to different frequencies. The resulting maps are useful to test data reduction algorithms, to understand residual systematics, to forecast maps produced by future instruments, to run Monte Carlo analysis for noise estimation, and more.

# Summary

The Python Sky Model (PySM) is a Python package used by Cosmic Microwave Background (CMB) experiments to simulate maps, in HEALPix [@gorski05; @healpy09] pixelization, of the various diffuse astrophysical components of Galactic emission relevant at CMB frequencies (i.e., dust, synchrotron, free-free and Anomalous Microwave Emission), as well as the CMB itself. These maps may be integrated over a given instrument bandpass and smoothed with a given instrument beam.
The template emission maps used by PySM are based on Planck [@planck18] and WMAP [@wmap13] data and are noise-dominated at small scales. Therefore, PySM simulation templates are smoothed to retain the large-scale information, and then supplemented with modulated Gaussian realizations at smaller scales. This strategy allows one to simulate data at higher resolution than the input maps.

PySM 2 [@pysm17], released in 2016, has become the de-facto standard for simulating Galactic emission; it is used, for example, by CMB-S4, Simons Observatory, LiteBird, PICO, CLASS, POLARBEAR, and other CMB experiments, as shown by the [80+ citations of the PySM 2 publication](https://scholar.google.com/scholar?start=0&hl=en&as_sdt=2005&sciodt=0,5&cites=16628417670342266167&scipsc=).
As the resolution of upcoming experiments increases, the PySM 2 software has started to show some limitations:

* Emission templates are provided at 7.9 arcminutes resolution (HEALPix $N_{side}=512$), while the next generation of CMB experiments will require sub-arcminute resolution.
* The software is implemented in pure `numpy`, meaning that it has significant memory overhead and is not multi-threaded, precluding simply replacing the current templates with higher-resolution versions.
* Emission templates are included in the PySM 2 Python package, which is still practical when each of the roughly 40 input maps is ~10 Megabytes, but will not be if they are over 1 Gigabyte.

The solution to these issues was to reimplement PySM from scratch focusing of these features:

* Reimplement all the models with the `numba` [@numba] Just-In-Time compiler for Python to reduce memory overhead and optimize performance: the whole integration loop of a template map over the frequency response of an instrument is performed in a single pass in automatically compiled and multi-threaded Python code.
* Use MPI through `mpi4py` to coordinate execution of PySM 3 across multiple nodes, this allows supporting template maps at a resolution up to 0.4 arcminutes (HEALPix $N_{side}=8192$).
* Rely on `libsharp` [@libsharp], a distributed implementation of spherical harmonic transforms, to smooth the maps with the instrument beam when maps are distributed over multiple nodes with MPI.
* Employ the data utilities infrastructure provided by `astropy` [@astropy2013; @astropy2018] to download the input templates and cache them when requested.

At this stage we strive to maintain full compatibility with PySM 2, therefore we implement the exact same astrophysical emission models with the same naming scheme. In the extensive test suite we compare the output of each PySM 3 model with the results obtained by PySM 2.

# Performance

As an example of the performance improvements achieved with PySM 3 over PySM 2, we run the following configuration:

* An instrument with 3 channels, with different beams, and a top-hat bandpass defined numerically at 10 frequency samples.
* A sky model with the simplest models of dust, synchrotron, free-free and AME [`a1,d1,s1,f1` in PySM terms].
* Execute on a 12-core Intel processor with 12 GB of RAM.

The following tables shows the walltime and peak memory usage of this simulation executed at the native PySM 2 resolution of $N_{side}=512$ and at two higher resolutions:

| Output $N_{side}$ | PySM 3        | PySM 2        |
|-------------------|---------------|---------------|
| 512               | 1m 0.7 GB     | 1m40s 1.45 GB |
| 1024              | 3m30s 2.3 GB  | 7m20s 5.5 GB  |
| 2048              | 16m10s 8.5 GB | Out of memory |

The models at $N_{side}=512$ have been tested to be equal given a relative tolerance of `1e-5`.

At the moment it is not very useful to run at resolutions higher than $N_{side}=512$ because there is no actual template signal at smaller scales. However, this demonstrates the performance improvements that will make working with higher resolution templates possible.

# Future work

PySM 3 opens the way to implement a new category of models at much higher resolution. However, instead of just upgrading the current models to smaller scales, we want to also update them with the latest knowledge of Galactic emission and gather feedback from each of the numerous CMB experiments. For this reason we are collaborating with the Panexperiment Galactic Science group to lead the development of the new class of models to be included in PySM 3.

# How to cite

If you are using PySM 3 for your work, please cite this paper for the software itself; for the actual emission modeling please also cite the original PySM 2 paper [@pysm17]. There will be a future paper on the generation of new PySM 3 astrophysical models.

# Acknowledgments

* This work was supported in part by NASA grant `80NSSC18K1487`.
* The software was tested, in part, on facilities run by the Scientific Computing Core of the Flatiron Institute.
* This research used resources of the National Energy Research Scientific Computing Center (NERSC), a U.S. Department of Energy Office of Science User Facility located at Lawrence Berkeley National Laboratory, operated under Contract No. `DE-AC02-05CH11231`.

# References
3.3.3 (unreleased)
==================

3.3.2 (2021-10-29)
==================

- `Improvements to install documentation <https://github.com/galsci/pysm/pull/93>`_
- Moved Github repository from `healpy/pysm` to `galsci/pysm`, under the Panexperiment Galactic Science group organization
- Changes in this release are related to the `JOSS Review <https://github.com/openjournals/joss-reviews/issues/3783>`_
- Turned `UserWarning` into `logging`, `Pull Request 88 <https://github.com/galsci/pysm/pull/88>`_

3.3.1 (2021-06-30)
==================

- Packaging: removed Poetry, using new Astropy package template and Github actions, `Pull Request 76 <https://github.com/galsci/pysm/pull/76>`_ and `77 <https://github.com/galsci/pysm/pull/77>`_.
- Docs: Reproduce PySM 2 template preprocessing for Synchrotron, `Pull Request 71 <https://github.com/galsci/pysm/pull/71>`_
- Docs: Reproduce PySM 2 template preprocessing for Dust, `Pull Request 66 <https://github.com/galsci/pysm/pull/66>`_

3.3.0 (2020-09-12)
==================

- Avoid an imcompatibility issue with ``numba``, see `Pull Request 63 <https://github.com/galsci/pysm/pull/63>`_
- Fix a severe bug in unit conversion with bandpass integration, which can give an overall scale error of a few percent at high frequency for all components, see `Issue 59 <https://github.com/galsci/pysm/issues/59>`_, also imported all bandpass integration tests from PySM 2 and added a comparison with the `tod2flux` tool by @keskitalo
- Removed support for `has_polarization` in interpolator, always return IQU map

3.2.2 (2020-06-23)
==================

- Fix packaging issue `importlib-resources` for python 3.6 was missing

3.2.1 (2020-06-05)
==================

- Renamed the package to `pysm3`, therefore now need to `import pysm3`
- Using `poetry` to build package and manage dependencies `PR 56 <https://github.com/galsci/pysm/pull/56>`_

3.2.0 (2020-04-15)
==================

First version with all models available in PySM 2

- Implemented HD2017 `d7` dust model `PR 37 <https://github.com/galsci/pysm/pull/37>`_
- Implemented HD2017 `d5` and `d8` dust models `PR 51 <https://github.com/galsci/pysm/pull/51>`_
- Improved documentation about Sky
- Implement local data folder `PR 53 <https://github.com/galsci/pysm/pull/53>`_

3.1.2 (2020-03-27)
==================

HD2017 `d7` dust model still being implemented

- Updated build/test setup to latest Astropy template `PR 47 <https://github.com/galsci/pysm/pull/47>`_
- Bugfix: `d6` model `PR 43 <https://github.com/galsci/pysm/pull/43>`_
- Bugfix: units other than GHz `PR 45 <https://github.com/galsci/pysm/pull/45>`_

3.1.0 (2019-12-11)
==================

- All emissions implemented except HD2017 `d7` dust

3.0.0 (2019-09-23)
==================

- Development release
|CI Tests| |Documentation Status| |PyPI| |Conda| |Astropy| |JOSS|

PySM 3
======

PySM generates full-sky simulations of Galactic emissions in intensity
and polarization relevant to CMB experiments. It is a large refactor of
`PySM 2 <https://github.com/bthorne93/PySM_public>`__ focused on
reducing memory usage, improving performance and run in parallel with
MPI.

See the documentation at https://pysm3.readthedocs.io

See changes in ``CHANGES.rst`` in the repository.

Related scientific papers
-------------------------

See `CITATION <https://github.com/galsci/pysm/blob/main/CITATION>`_

* `The Python Sky Model 3 software (Zonca et al, 2021) <https://arxiv.org/abs/2108.01444>`_
* `The Python Sky Model: software for simulating the Galactic microwave sky (Thorne et al, 2017) <https://arxiv.org/abs/1608.02841>`_

Install
-------

See the `documentation <https://pysm3.readthedocs.io/en/latest/#installation>`_

* Install with ``pip install .`` or with ``pip install .[test]`` to also install the requirements for running tests
* Optionally, if you have an MPI environment available and you would like to test the MPI capabilities of PySM, install ``mpi4py`` and ``libsharp``, check the documentation link above for more details.
* Check code style with ``tox -e codestyle``
* Test with ``pytest`` or ``tox -e test``
* Building docs requires ``pandoc``, not the python package, the actual ``pandoc`` command line tool, install it with conda or your package manager
* Build docs locally with ``tox -e build_docs``

Support
-------

For any question or issue with the software `open an issue <https://github.com/galsci/pysm/issues/>`_.

Release
-------

* Tag the new version with git
* ``pip install build --upgrade``
* ``python -m build --sdist --wheel .``
* ``twine upload dist/*``

.. |CI Tests| image:: https://github.com/galsci/pysm/actions/workflows/ci_tests.yml/badge.svg
   :target: https://github.com/galsci/pysm/actions/workflows/ci_tests.yml
.. |Documentation Status| image:: https://readthedocs.org/projects/pysm3/badge/?version=latest
   :target: https://pysm3.readthedocs.io/en/latest/?badge=latest
.. |PyPI| image:: https://img.shields.io/pypi/v/pysm3
   :target: https://pypi.org/project/pysm3/
.. |Conda| image:: https://img.shields.io/conda/vn/conda-forge/pysm3
   :target: https://anaconda.org/conda-forge/pysm3
.. |Astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
   :target: http://www.astropy.org/
.. |JOSS| image:: https://joss.theoj.org/papers/8f2d6c3bbf6cbeffbb403a1207fa8de7/status.svg
   :target: https://joss.theoj.org/papers/8f2d6c3bbf6cbeffbb403a1207fa8de7
Data directory
==============

This directory contains data files included with the package source
code distribution. Note that this is intended only for relatively small files
- large files should be externally hosted and downloaded as needed.
Contributing
============

When contributing to this repository a large modification of the
software, please first discuss it via issue.

Please note we have a code of conduct, please follow it in all your
interactions with the project.

Pull Request Process
--------------------

1. Describe the intent of the PR in the body
2. Make sure the feature you implement has significant automated test
   coverage
3. Add a line to CHANGES.rst with a description of the modification with
   link to the PR
4. Wait at least 2 weeks for a maintainer to review the PR before
   soliciting

Code of Conduct
---------------

Our Pledge
~~~~~~~~~~

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our
project and our community a harassment-free experience for everyone,
regardless of age, body size, disability, ethnicity, gender identity and
expression, level of experience, nationality, personal appearance, race,
religion, or sexual identity and orientation.

Our Standards
~~~~~~~~~~~~~

Examples of behavior that contributes to creating a positive environment
include:

-  Using welcoming and inclusive language
-  Being respectful of differing viewpoints and experiences
-  Gracefully accepting constructive criticism
-  Focusing on what is best for the community
-  Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

-  The use of sexualized language or imagery and unwelcome sexual
   attention or advances
-  Trolling, insulting/derogatory comments, and personal or political
   attacks
-  Public or private harassment
-  Publishing others’ private information, such as a physical or
   electronic address, without explicit permission
-  Other conduct which could reasonably be considered inappropriate in a
   professional setting

Our Responsibilities
~~~~~~~~~~~~~~~~~~~~

Project maintainers are responsible for clarifying the standards of
acceptable behavior and are expected to take appropriate and fair
corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit,
or reject comments, commits, code, wiki edits, issues, and other
contributions that are not aligned to this Code of Conduct, or to ban
temporarily or permanently any contributor for other behaviors that they
deem inappropriate, threatening, offensive, or harmful.

Scope
~~~~~

This Code of Conduct applies both within project spaces and in public
spaces when an individual is representing the project or its community.
Examples of representing a project or community include using an
official project e-mail address, posting via an official social media
account, or acting as an appointed representative at an online or
offline event. Representation of a project may be further defined and
clarified by project maintainers.

Enforcement
~~~~~~~~~~~

Instances of abusive, harassing, or otherwise unacceptable behavior may
be reported by contacting the project team. All complaints will be
reviewed and investigated and will result in a response that is deemed
necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an
incident. Further details of specific enforcement policies may be posted
separately.

Project maintainers who do not follow or enforce the Code of Conduct in
good faith may face temporary or permanent repercussions as determined
by other members of the project’s leadership.

Attribution
~~~~~~~~~~~

This Code of Conduct is adapted from the `Contributor
Covenant <http://contributor-covenant.org>`__, version 1.4, available at
`http://contributor-covenant.org/version/1/4 <http://contributor-covenant.org/version/1/4/>`__
.. _models:

Summary of Models
*****************

For all details of the models also check the `presets.cfg file <https://github.com/galsci/pysm/blob/main/pysm3/data/presets.cfg>`_.

Input templates
===============

The PySM input templates are not stored in the Python package, they are downloaded on the fly when requested using the machinery from the `astropy.io.data` package.
The templates are stored in the `pysm-data <https://github.com/galsci/pysm-data>`_ repository on Github and published `via the cmb organization at NERSC <https://portal.nersc.gov/project/cmb/pysm-data/>`_.

PySM has an hardcoded path to access the files at NERSC, therefore when running on NERSC supercomputers or on Jupyter@NERSC, files are accessed directly.

When PySM is executed in another environment, the necessary templates are downloaded and cached in the home folder by astropy, see `the astropy.utils.data documentation <https://docs.astropy.org/en/stable/utils/data.html>`_ for the available configuration options.

Another option is to download the whole data repository, currently ~700MB (it will be way larger in the future), and define an environment variable::

    git clone https://github.com/galsci/pysm-data /data/pysm-data
    export PYSM_LOCAL_DATA=/data/pysm-data

Reproduce PySM 2 template preprocessing
=======================================

The PySM 2 paper describes how input data (e.g. component separation maps from Planck) have been processed, for example
most templates have been smoothed to remove high $\ell$ noise and then added small scale fluctuations.
Here we try to reproduce the process to clarify it:

* `Dust polarization templates from Planck Commander component separation outputs <preprocess-templates/reproduce_pysm2_dust_pol.html>`_, used in all dust models from `d0` to `d8` except `d4`
* `Synchrotron polarization templates from WMAP low frequency maps <preprocess-templates/reproduce_pysm2_sync_pol.html>`_, used in all synchrotron models


Dust
====

- **d1**: Thermal dust is modelled as a single-component modified black body (mbb). We use dust templates for emission at 545 GHz in intensity and 353 GHz in polarisation from the Planck-2015 analysis, and scale these to different frequencies with a mbb spectrum using the spatially varying temperature and spectral index obtained from the Planck data using the Commander code (Planck Collaboration 2015, arXiv:1502.01588). Note that it therefore assumes the same spectral index for polarization as for intensity. The input intensity template at 545 GHz is simply the available 2048 product degraded to nside 512. The polarization templates have been smoothed with a Gaussian kernel of FWHM 2.6 degrees, and had small scales added via the procedure described in the accompanying paper.

- **d0**: Simplified version of the **d1** model with a fixed spectral index of 1.54 and a black body temperature of 20 K.

- **d2** (**d3**): emissivity that varies spatially on degree scales, drawn from a Gaussian with beta=1.59 \pm 0.2 (0.3). A Gaussian variation is not physically motivated, but amount of variation consistent with Planck.

- **d4**: a generalization of model 1 to multiple dust populations. It has been found that a two component model is still a good fit to the Planck data. This option uses the two component model from Finkbeiner, D. P., Davis, M., & Schlegel, D. J. 1999, Astrophysical Journal, 524, 867.

- **d5**: implementation of the dust model described in Hensley and Draine 2017.
  
- **d6**: implementation of the frequency decorrelation of dust, modelling the impact of averaging over spatially varying dust spectral indices both unresolved and along the line of sight. We take an analytic frequency covariance (Vansyngel 2016 arXiv:1611.02577) to calculate the resulting frequency dependence. The user specifies a single parameter, the correlation length. The smaller the correlation length, the larger the decorrelation. This parameter is constant across the sky.

- **d7**: modification of `d5` with iron inclusions in the grain composition.

- **d8**: simplified version of `d7` where the interstellar radiation field (ISRF) strength, instead of being a random realization, is fixed at 0.2.  This corresponds reasonably well to a Modifield Black Body model with temperature of 20K and an index of 1.54.

- **d12**: 3D model of polarized dust emission with 6 layers, based on the paper `"A 3-D model of polarised dust emission in the Milky Way" <https://arxiv.org/abs/1706.04162>`_, named MKD based on the names of the authors. Each layer has different templates, spectral index and dust temperature. All maps were generated at N_side 2048 with the Planck Sky Model (PSM) by Jacques Delabrouille.

Synchrotron
===========

- **s1**: A power law scaling is used for the synchrotron emission, with a spatially varying spectral index. The emission templates are the Haslam 408 MHz, 57' resolution data reprocessed by Remazeilles et al 2015 MNRAS 451, 4311, and the WMAP 9-year 23 GHz Q/U maps (Bennett, C.L., et.al., 2013, ApJS, 208, 20B). The polarization maps have been smoothed with a Gaussian kernel of FWHM 5 degrees and had small scales added. The intensity template has had small scales added straight to the template. The details of the small scale procedure is outlined in the accompanying paper. The spectral index map was derived using a combination of the Haslam 408 MHz data and WMAP 23 GHz 7-year data (Miville-Deschenes, M.-A. et al., 2008, A&A, 490, 1093). The same scaling is used for intensity and polarization. This is the same prescription as used in the Planck Sky Model's v1.7.8 'power law' option (Delabrouille et al. A&A 553, A96, 2013), but with the Haslam map updated to the Remazeilles version. A 'curved power law' model is also supported with a single isotropic curvature index. The amplitude of this curvature is taken from Kogut, A. 2012, ApJ, 753, 110.

- **s2**: synchrotron index steepens off the Galactic plane, from -3.0 in the plane to -3.3 off the plane. Consistent with WMAP.

- **s3**: a power law with a curved index. The model uses the same index map as the nominal model, plus a curvature term. We use the best-fit curvature amplitude of -0.052 found in Kogut, A. 2012, ApJ, 753, 110, pivoted at 23 GHz.


AME
===

- **a1**: We model the AME as a sum of two spinning dust populations based on the Commander code (Planck Collaboration 2015, arXiv:1502.01588). A component is defined by a degree-scale emission template at a reference frequency and a peak frequency of the emission law. Both populations have a spatially varying emission template, one population has a spatially varying peak frequency, and the other population has a spatially constant peak frequency. The emission law is generated using the SpDust2 code (Ali-Haimoud 2008). The nominal model is unpolarized. We add small scales to the emission maps, the method is outlined in the accompanying paper.
  
- **a2**: AME has 2% polarization fraction. Polarized maps simulated with thermal dust angles and nominal AME intensity scaled globally by polarization fraction. Within WMAP/Planck bounds.


Free-free
=========

- **f1**: We model the free-free emission using the analytic model assumed in the Commander fit to the Planck 2015 data (Draine 2011 'Physics of the Interstellar and Intergalactic Medium') to produce a degree-scale map of free-free emission at 30 GHz. We add small scales to this using a procedure outlined in the accompanying paper. This map is then scaled in frequency by applying a spatially constant power law index of -2.14.

CMB
===

- **c1**: A lensed CMB realisation is computed using Taylens, a code to compute a lensed CMB realisation using nearest-neighbour Taylor interpolation (`taylens <https://github.com/amaurea/taylens>`_; Naess, S. K. and Louis, T. JCAP 09 001, 2013, astro-ph/1307.0719). This code takes, as an input, a set of unlensed Cl's generated using `CAMB <http://www.camb.info/>`_. The params.ini is in the Ancillary directory. There is a pre-computed CMB map provided at Nside 512.

PySM Documentation
==================

PySM 3 generates full-sky simulations of Galactic foregrounds in intensity and polarization relevant for CMB experiments. The components simulated are: thermal dust, synchrotron, AME, free-free, and CMB at a given HEALPix $N_{side}$, with an option to integrate over a bandpass and to smooth with a given beam.

There is scope for a few options for the model for each component, attempting to be consistent with current data.

Currently much of the available data is limited in resolution at degree-scale. We therefore make efforts to provide reasonable small-scale simulations to extend the data to higher multipoles. The details of the procedures developed can be found in the accompanying paper.

If you are using this code please cite the PySM 2 paper `Thorne at al <https://arxiv.org/abs/1608.02841>`_

This code is based on the large-scale Galactic part of `Planck Sky Model <http://www.apc.univ-paris7.fr/~delabrou/PSM/psm.html>`_ (`Delabrouille 2012 <https://arxiv.org/abs/1207.3675>`_) code and uses some of its inputs.

Models
------

Each model is identified with a letter and a number, the letter indicates the kind of emission and the number the type of model, generally in order of complexity starting at 1. For example for dust we start with `d1` based on Planck commander results, `d4` has 2 dust populations and `d6` implements a model of dust frequency decorrelation.

For example free-free:

    **f1**: We model the free-free emission using the analytic model assumed in the Commander fit to the Planck 2015 data (Draine 2011 Physics of the Interstellar and Intergalactic Medium) to produce a degree-scale map of free-free emission at 30 GHz. We add small scales to this using a procedure outlined in the accompanying paper. This map is then scaled in frequency by applying a spatially constant power law index of -2.14.

* See the PySM 2 paper `Thorne at al <https://arxiv.org/abs/1608.02841>`_
* See the documentation about :ref:`models`

.. toctree::
  :maxdepth: 2

  models

Dependencies
============

PySM is written in pure Python, multi-threading and Just-In-Time compilation are provided by `numba`.
The required packages are:

* `healpy`
* `numba`
* `toml`
* `astropy`
* `importlib_metadata` just for Python 3.7

How many threads are used by `numba` can be controlled by setting::

    export OMP_NUM_THREADS=2
    export NUMBA_NUM_THREADS=2

For debugging purposes, you can also completely disable the `numba` JIT compilation and have the code just use plain
Python, which is slower, but easier to debug::

    export NUMBA_DISABLE_JIT=1

Run PySM with MPI support
-------------------------

PySM 3 is capable of running across multiple nodes using a MPI communicator.

See the details in the `MPI section of the tutorial <mpi.ipynb>`_.

In order to run in parallel with MPI, it also needs a functioning MPI environment and:

* `mpi4py`

MPI-Distributed smoothing (optional) requires `libsharp`, it is easiest to install the conda package::

    conda install -c conda-forge libsharp=*=*openmpi*

It also has a `mpich` version::

    conda install -c conda-forge libsharp=*=*mpich*


Installation
============

The easiest way to install the last release is to use `conda`::

    conda install -c conda-forge pysm3

or `pip`::

    pip install pysm3

See the `changelog on Github <https://github.com/galsci/pysm/blob/main/CHANGES.rst>`_ for details about what is included in each release.

Install at NERSC
----------------

Optionally replace with a newer anaconda environment::

    module load python/3.7-anaconda-2019.10
    conda create -c conda-forge -n pysm3 pysm3 python=3.7 ipython
    conda activate pysm3
    module unload python

Development install
-------------------

The development version is available in the `main` branch of the `GitHub repository <https://github.com/galsci/pysm>`_,
you can clone and install it with::

    pip install .

Create a development installation with::

    pip install -e .

Install the requirements for testing with::

    pip install -e .[test]

Execute the unit tests with::

    pytest

Configure verbosity
-------------------

PySM uses the `logging` module to configure its verbosity,
by default it will only print warnings and errors, to configure logging
you can access the "pysm3" logger with::

    import logging
    log = logging.getLogger("pysm3")

configure the logging level::

    log.setLevel(logging.DEBUG)

redirect the logs to the console::

    handler = logging.StreamHandler()
    log.addHandler(handler)

or customize their format::

    log_format="%(name)s - %(levelname)s - %(message)s"
    formatter = logging.Formatter(log_format)
    handler.setFormatter(formatter)

For more details see the `Python documentation <https://docs.python.org/3/library/logging.html>`_.

Tutorials
=========

All the tutorials are Jupyter Notebooks and can be accessed `from the repository <https://github.com/galsci/pysm/tree/main/docs>`_:

.. toctree::
  :maxdepth: 2

  basic_use
  model_data
  bandpass_integration
  smoothing_coord_rotation
  customize_components
  mpi

Contributing
============

.. toctree::
  :maxdepth: 2

  contributing

Reference/API
=============

.. automodapi:: pysm3
