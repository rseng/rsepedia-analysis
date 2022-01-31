---
title: '21cmFAST v3: A Python-integrated C code for generating 3D realizations of the cosmic 21cm signal.'
tags:
  - Python
  - astronomy
  - cosmology
  - simulation
authors:
  - name: Steven G. Murray
    orcid: 0000-0003-3059-3823
    affiliation: 1
  - name: Bradley Greig
    orcid: 0000-0002-4085-2094
    affiliation: 2, 3
  - name: Andrei Mesinger
    orcid: 0000-0003-3374-1772
    affiliation: 4
  - name: Julian B. Muñoz
    orcid: 0000-0002-8984-0465
    affiliation: 5
  - name: Yuxiang Qin
    orcid: 0000-0002-4314-1810
    affiliation: 4
  - name: Jaehong Park
    orcid: 0000-0003-3095-6137
    affiliation: 4, 7
  - name: Catherine A. Watkinson
    orcid: 0000-0003-1443-3483
    affiliation: 6

affiliations:
 - name: School of Earth and Space Exploration, Arizona State University, Phoenix, USA
   index: 1
 - name: ARC Centre of Excellence for All-Sky Astrophysics in 3 Dimensions (ASTRO 3D)
   index: 2
 - name: School of Physics, University of Melbourne, Parkville, VIC 3010, Australia
   index: 3
 - name: Scuola Normale Superiore, Piazza dei Cavalieri 7, 56126 Pisa, Italy
   index: 4
 - name: Department of Physics, Harvard University, 17 Oxford St., Cambridge, MA, 02138, USA
   index: 5
 - name: School of Physics and Astronomy, Queen Mary University of London, G O Jones Building, 327 Mile End Road, London, E1 4NS, UK
   index: 6
 - name: School of Physics, Korea Institute for Advanced Study, 85 Hoegiro, Dongdaemun-gu, Seoul, 02455, Republic of Korea
   index: 7

date: 24 Sep 2020
bibliography: paper.bib
---

# Summary

The field of 21-cm cosmology -- in which the hyperfine spectral line of neutral hydrogen
(appearing at the rest-frame wavelength of 21 cm) is mapped over large swathes of the
Universe's history -- has developed radically over the last decade.
The promise of the field is to revolutionize our
knowledge of the first stars, galaxies, and black holes through the timing and patterns
they imprint on the cosmic 21-cm signal.
In order to interpret the eventual observational data, a range of physical models have
been developed -- from simple analytic models of the global history of hydrogen reionization,
through to fully hydrodynamical simulations of the 3D evolution of the brightness
temperature of the spectral line.
Between these extremes lies an especially versatile middle-ground: fast semi-numerical
models that approximate the full 3D evolution of the relevant fields: density, velocity,
temperature, ionization, and radiation (Lyman-alpha, neutral hydrogen 21-cm, etc.).
These have the advantage of being comparable to the full first-principles
hydrodynamic simulations, but significantly quicker to run; so much so that they can
be used to produce thousands of realizations on scales comparable to those observable
by upcoming low-frequency radio telescopes, in order to explore the very wide
parameter space that still remains consistent with the data.


Amongst practitioners in the field of 21-cm cosmology, the `21cmFAST` program has become
the *de facto* standard for such semi-numerical simulators.
`21cmFAST` [@mesinger2007; @mesinger2010] is a high-performance C code that uses the
excursion set formalism [@furlanetto2004] to
identify regions of ionized hydrogen atop a cosmological density field evolved using
first- or second-order Lagrangian perturbation theory [@zeldovich1970; @scoccimarro2002],
tracking the thermal and ionization state of the intergalactic medium, and computing
X-ray, soft UV and ionizing UV cosmic radiation fields based on parametrized galaxy models.
For example, the following figure contains slices of lightcones (3D fields in which one
axis corresponds to both spatial *and* temporal evolution) for the various
component fields produced by `21cmFAST`.

![Sample of Component Fields output by 21cmFAST. Cosmic evolution occurs from bottom to top. From left to right, quantities shown are: (i) dark matter overdensity field; (ii) Lyman-alpha flux; (iii) Lyman-Werner flux; (iv) X-ray heating rate; (v) locally-averaged UVB; (vi) critical halo mass for star formation in Atomically Cooled Galaxies; (vii) critical halo mass for star formation in Molecularly Cooled Galaxies; (viii) cumulative number of recombinations per baryon; (ix) neutral hydrogen fraction; and (x) the 21cm brightness temperature. A high-resolution version of this figure is available at http://homepage.sns.it/mesinger/Media/lightcones_minihalo.png](yuxiangs-plot-small.png){height=450px}

However, `21cmFAST` is a highly specialized code, and its implementation has been
quite specific and relatively inflexible.
This inflexibility makes it difficult to modify the behaviour of the code without detailed
knowledge of the full system, or disrupting its workings.
This lack of modularity within the code has led to
widespread code "branching" as researchers hack new physical features of interest
into the C code; the lack of a streamlined API has led derivative codes which run
multiple realizations of `21cmFAST` simulations [such as the Monte Carlo simulator,
`21CMMC`, @greig2015] to re-write large portions of the code in order to serve their purpose.
It is thus of critical importance, as the field moves forward in its understanding -- and
the range and scale of physical models of interest continues to increase -- to
reformulate the `21cmFAST` code in order to provide a fast, modular, well-documented,
well-tested, stable simulator for the community.

# Features of 21cmFAST v3

This paper presents `21cmFAST` v3+, which is formulated to follow these essential
guiding principles.
While keeping the same core functionality of previous versions of `21cmFAST`, it has
been fully integrated into a Python package, with a simple and intuitive interface, and
a great deal more flexibility.
At a higher level, in order to maintain best practices, a community of users and
developers has coalesced into a formal collaboration which maintains the project via a
Github organization.
This allows the code to be consistently monitored for quality, maintaining high test
coverage, stylistic integrity, dependable release strategies and versioning,
and peer code review.
It also provides a single point-of-reference for the community to obtain the code,
report bugs and request new features (or get involved in development).

A significant part of the work of moving to a Python interface has been the
development of a robust series of underlying Python structures which handle the passing
of data between Python and C via the `CFFI` library.
This foundational work provides a platform for future versions to extend the scientific
capabilities of the underlying simulation code.
The primary *new* usability features of `21cmFAST` v3+ are:

* Convenient (Python) data objects which simplify access to and processing of the various
  fields that form the brightness temperature.
* Enhancement of modularity: the underlying C functions for each step of the simulation
  have been de-coupled, so that arbitrary functionality can be injected into the process.
* Conversion of most global parameters to local structs to enable this modularity, and
  also to obviate the requirement to re-compile in order to change parameters.
* Simple `pip`-based installation.
* Robust on-disk caching/writing of data, both for efficiency and simplified reading of
  previously processed data (using HDF5).
* Simple high-level API to generate either coeval cubes (purely spatial 3D fields defined
  at a particular time) or full lightcone data (i.e. those coeval cubes interpolated over
  cosmic time, mimicking actual observations).
* Improved exception handling and debugging.
* Convenient plotting routines.
* Simple configuration management, and also more intuitive management for the
  remaining C global variables.
* Comprehensive API documentation and tutorials.
* Comprehensive test suite (and continuous integration).
* Strict semantic versioning^[https://semver.org].

While in v3 we have focused on the establishment of a stable and extendable infrastructure,
we have also incorporated several new scientific features, appearing in separate papers:

* Generate transfer functions using the `CLASS` Boltzmann code [@Lesgourgues2011].
* Simulate the effects of relative velocities between dark matter and Baryons [@munoz2019a; @munoz2019b].
* Correction for non-conservation of ionizing photons (Park, Greig et al., *in prep*).
* Include molecularly cooled galaxies with distinct properties [@qin2020]
* Calculate rest-frame UV luminosity functions based on parametrized galaxy models.

`21cmFAST` is still in very active development.
Amongst further usability and performance improvements,
future versions will see several new physical models implemented,
including milli-charged dark matter models [@Munoz2018] and forward-modelled CMB
auxiliary data [@qin2020a].

In addition, `21cmFAST` will be incorporated into large-scale inference codes, such as
`21CMMC`, and is being used to create large data-sets for inference via machine learning.
We hope that with this new framework, `21cmFAST` will remain an important component
 of 21-cm cosmology for years to come.

# Examples

`21cmFAST` supports installation using `conda`, which means installation is as simple
as typing `conda install -c conda-forge 21cmFAST`. The following example can then
be run in a Python interpreter.

In-depth examples can be found in the official documentation.
As an example of the simplicity with which a full lightcone may be produced with
the new `21cmFAST` v3, the following may be run in a Python interpreter (or Jupyter
notebook):

```python
import py21cmfast as p21c

lightcone = p21c.run_lightcone(
    redshift=6.0,              # Minimum redshift of lightcone
    max_redshift=30.0,
    user_params={
        "HII_DIM": 150,        # N cells along side in output cube
        "DIM": 400,            # Original high-res cell number
        "BOX_LEN": 300,        # Size of the simulation in Mpc
    },
    flag_options={
        "USE_TS_FLUCT": True,  # Don't assume saturated spin temp
        "INHOMO_RECO": True,   # Use inhomogeneous recombinations
    },
    lightcone_quantities=(     # Components to store as lightcones
        "brightness_temp",
        "xH_box",
        "density"
    ),
    global_quantities=(        # Components to store as mean
        "xH_box",              # values per redshift
        "brightness_temp"
    ),
)

# Save to a unique filename hashing all input parameters
lightcone.save()

# Make a lightcone sliceplot
p21c.plotting.lightcone_sliceplot(lightcone, "brightness_temp")
```

![Brightness temperature lightcone produced by the example code in this paper.](lightcone.pdf){height=300px}

```python
# Plot a global quantity
p21c.plotting.plot_global_history(lightcone, "xH")
```

![Globally volume-averaged hydrogen neutral fraction produced by the example code in this paper.](xH_history.pdf){height=300px}

# Performance

Despite being a Python code, `21cmFAST` v3 does not diminish the performance of previous
pure-C versions. It utilises `CFFI` to provide the interface to the C-code through
Python, which is managed by some custom Python classes that oversee the construction and
memory allocation of each C `struct`.

OpenMP parallelization is enabled within the C-code, providing excellent speed-up for
large simulations when performed on high-performance machines.

A simple performance comparison between v3 and v2.1 (the last pure-C version), running
a light-cone simulation over a redshift range between 35 and 5 (92 snapshots) with spin
temperature fluctuations (`USE_TS_FLUCT`), inhomogeneous recombinations
(`INHOMO_RECO`), FFTW Wisdoms (`USE_FFTW_WISDOM`) and interpolation tables
(`USE_INTERPOLATION_TABLES`),
with a resolution of `HII_DIM=250` cells, and `DIM=1000` cells for the initial conditions,
on an Intel(R) Xeon(R) CPU (E5-4657L v2 @ 2.40GHz) with 16 shared-memory cores, reveals
that a clock time of 7.63(12.63) hours and a maximum RAM of 224(105) gigabytes are needed
for v3(v2.1).

Note that while a full light-cone simulation can be expensive to perform,
it only takes 2-3min to calculate a Coeval box (excluding the initial conditions).
For instance, the aforementioned timing for v3 includes 80 minutes to generate the
initial condition, which also dominates the maximum RAM required, with an additional
~4 minutes per snapshot to calculate all required fields of perturbation, ionization,
spin temperature and brightness temperature.

To guide the user, we list some performance benchmarks for variations on this simulation,
run with `21cmFAST` v3.0.2. Note that these benchmarks are subject to change as new
minor versions are delivered; in particular, operational modes that reduce maximum
memory consumption are planned for the near future.

| Variation                                       | Time (hr) | Memory (GB) |
| ----------------------------------------------- | --------- | ----------- |
| Reference                                       | 7.63      | 224         |
| Single Core                                     | 14.77     | 224         |
| 4 Shared-memory Cores                           | 7.42      | 224         |
| 64 Shared-memory Cores                          | 9.60      | 224         |
| Higher Resolution <br />(HII_DIM=500, DIM=2000) | 68.37     | 1790        |
| Lower Resolution <br />(HII_DIM=125, DIM=500)   | 0.68      | 28          |
| No Spin Temperature                             | 4.50      | 224         |
| Use Mini-Halos                                  | 11.57     | 233         |
| No FFTW Wisdoms                                 | 7.33      | 224         |

At this time, the `21cmFAST` team suggests using 4 or fewer shared-memory cores.
However, it is worth noting that as performance does vary on different machines,
users are recommended to calculate their own scalability.

# Acknowledgements

This work was supported in part by the European Research Council
(ERC) under the European Union’s Horizon 2020 research
and innovation programme (AIDA – #638809). The results
presented here reflect the authors’ views; the ERC is not
responsible for their use. JBM was partially supported by NSF grant AST-1813694.
Parts of this research were supported by the European Research Council under ERC grant
number 638743-FIRSTDAWN. Parts of this research were supported by the Australian Research
Council Centre of Excellence for
All Sky Astrophysics in 3 Dimensions (ASTRO 3D), through project number CE170100013.
JP was supported in part by a KIAS individual Grant (PG078701) at Korea Institute for
Advanced Study.

# References
============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

Bug reports/Feature Requests/Feedback/Questions
===============================================
It is incredibly helpful to us when users report bugs, unexpected behaviour, or request
features. You can do the following:

* `Report a bug <https://github.com/21cmFAST/21cmFAST/issues/new?template=bug_report.md>`_
* `Request a Feature <https://github.com/21cmFAST/21cmFAST/issues/new?template=feature_request.md>`_
* `Ask a Question <https://github.com/21cmFAST/21cmFAST/issues/new?template=question.md>`_

When doing any of these, please try to be as succinct, but detailed, as possible, and use
a "Minimum Working Example" whenever applicable.

Documentation improvements
==========================

``21cmFAST`` could always use more documentation, whether as part of the
official ``21cmFAST`` docs, in docstrings, or even on the web in blog posts,
articles, and such. If you do the latter, take the time to let us know about it!

High-Level Steps for Development
================================

This is an abbreviated guide to getting started with development of ``21cmFAST``,
focusing on the discrete high-level steps to take. See our
`notes for developers <https://21cmfast.readthedocs.org/en/latest/notes_for_developers>`_
for more details about how to get around the ``21cmFAST`` codebase and other
technical details.

There are two avenues for you to develop ``21cmFAST``. If you plan on making significant
changes, and working with ``21cmFAST`` for a long period of time, please consider
becoming a member of the 21cmFAST GitHub organisation (by emailing any of the owners
or admins). You may develop as a member or as a non-member.

The difference between members and non-members only applies to the first step
of the development process.

Note that it is highly recommended to work in an isolated python environment with
all requirements installed from ``environment_dev.txt``. This will also ensure that
pre-commit hooks will run that enforce the ``black`` coding style. If you do not
install these requirements, you must manually run ``black`` before committing your changes,
otherwise your changes will likely fail continuous integration.

As a *member*:

1. Clone the repo::

    git clone git@github.com:21cmFAST/21cmFAST.git

As a *non-member*:

1. First fork ``21cmFAST <https://github.com/21cmFAST/21cmFAST>``_
   (look for the "Fork" button), then clone the fork locally::

    git clone git@github.com:your_name_here/21cmFAST.git

The following steps are the same for both *members* and *non-members*:

2. Install a fresh new isolated environment::

       conda create -n 21cmfast python=3
       conda activate 21cmfast

3. Install the *development* requirements for the project::

    conda env update -f environment_dev.yml

4. Install 21cmFAST. See `Installation <./installation.html>`_ for more details.::

    pip install -e .

4. Install pre-commit hooks::

    pre-commit install

5. Create a branch for local development::

    git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally. **Note: as a member, you _must_ do step 5. If you
   make changes on master, you will _not_ be able to push them**.

6. When you're done making changes, run ``pytest`` to check that your changes didn't
   break things. You can run a single test or subset of tests as well (see pytest docs)::

    pytest

7. Commit your changes and push your branch to GitHub::

    git add .
    git commit -m "Your detailed description of your changes."
    git push origin name-of-your-bugfix-or-feature

   Note that if the commit step fails due to a pre-commit hook, *most likely* the act
   of running the hook itself has already fixed the error. Try doing the ``add`` and
   ``commit`` again (up, up, enter). If it's still complaining, manually fix the errors
   and do the same again.

8. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

If you need some code review or feedback while you're developing the code just make the
pull request. You can mark the PR as a draft until you are happy for it to be merged.
Changelog
=========

dev-version
-----------

Fixed
-----

* docs not compiling on RTD due to missing ``scipy.integrate`` mock module
* Updated matplotlib removed support for setting vmin/vmax and norm. Now passes vmin/vmax
  to the norm() constructor.

v3.1.3 [27 Oct 2021]
----------------------

* Fixed ``FAST_FCOLL_TABLES`` so it only affects MCGs and not ACGs. Added tests of this
  flag for high and low z separately.

v3.1.2 [14 Jul 2021]
----------------------

Internals
~~~~~~~~~
* ``MINIMIZE_MEMORY`` flag significantly reduces memory without affecting performance much,
  by changing the way some arrays are allocated and accessed in C. (#224)

Change
~~~~~~

* Updated ``USE_INTERPOLATION_TABLES`` to be default True. This makes much more sense as
  a default value. Until v4, a warning will be raised if it is not set explicitly.


v3.1.1 [13 Jun 2021]
----------------------

Fixed
~~~~~

* Bug in deployment to PyPI.

v3.1.0 [13 Jun 2021]
----------------------

Added
~~~~~
* Ability to access all evolutionary Coeval components, either from the end Coeval
  class, or the Lightcone.
* Ability to gather all evolutionary antecedents from a Coeval/Lightcone into the one
  file.
* ``FAST_FCOLL_TABLES`` in ``UserParams`` which improves speeds quite significantly for
  ~<10% accuracy decrease.
* Fast and low-memory generation of relative-velocity (vcb) initial conditions. Eliminated hi-res vcb boxes, as they are never needed.
* Also output the mean free path (i.e. MFP_box in IonizedBox).
* Added the effect of DM-baryon relative velocities on PopIII-forming minihaloes. This now provides the correct background evolution jointly with LW feedback. It gives rise to velocity-induced acoustic oscillations (VAOs) from the relative-velocity fluctuations. We also follow a more flexible parametrization for LW feedback in minihaloes, following new simulation results, and add a new index ALPHA_STAR_MINI for minihaloes, now independent of regular ACGs.
* New ``hooks`` keyword to high-level functions, that are run on the completion of each computational step, and can
  be used to more generically write parts of the data to file.
* Ability to pass a function to ``write=`` to write more specific aspects of the data (internally, this will be put into the ``hooks`` dictionary).
* ``run_lightcone`` and ``run_coeval`` use significantly less memory by offloading initial conditions and perturb_field instances to disk if possible.

Fixed
~~~~~
* Bug in 2LPT when ``USE_RELATIVE_VELOCITIES=True`` [Issue #191, PR #192]
* Error raised when redshifts are not in ascending order [Issue #176, PR #177]
* Errors when ``USE_FFTW_WISDOM`` is used on some systems [Issue #174, PR #199]
* Bug in ComputeIonizedBox causing negative recombination rate and ring structure in ``Gamma12_box`` [Issue #194, PR #210]
* Error in determining the wisdom file name [Issue #209, PR#210]
* Bug in which cached C-based memory would be read in and free'd twice.

Internals
~~~~~~~~~

* Added ``dft.c``, which makes doing all the cubic FFTs a lot easier and more consistent. [PR #199]
* More generic way of keeping track of arrays to be passed between C and Python, and their shape in Python, using ``_get_box_structures``.
  This also means that the various boxes can be queried before they are initialized and computed.
* More stringent integration tests that test each array, not just the final brightness temperature.
* Ability to plot the integration test data to more easily identify where things have gone wrong (use ``--plots`` in the ``pytest`` invocation).
* Nicer CLI interface for ``produce_integration_test_data.py``. New options to ``clean`` the ``test_data/`` directory,
  and also test data is saved by user-defined key rather than massive string of variables.
* Nicer debug statements before calls to C, for easily comparing between versions.
* Much nicer methods of keeping track of array state (in memory, on disk, c-controlled, etc.)
* Ability to free C-based pointers in a more granular way.

v3.0.3
------

Added
~~~~~
* ``coeval_callback`` and ``coeval_callback_redshifts`` flags to the ``run_lightcone``.
  Gives the ability to run arbitrary code on ``Coeval`` boxes.
* JOSS paper!
* ``get_fields`` classmethod on all output classes, so that one can easily figure out
  what fields are computed (and available) for that class.

Fixed
~~~~~
* Only raise error on non-available ``external_table_path`` when actually going to use it.

v3.0.2
------

Fixed
-----
* Added prototype functions to enable compilation for some standard compilers on MacOS.

v3.0.1
------
Modifications to the internal code structure of 21cmFAST

Added
~~~~~
* Refactor FFTW wisdom creation to be a python callable function


v3.0.0
------
Complete overhaul of 21cmFAST, including a robust python-wrapper and interface,
caching mechanisms, and public repository with continuous integration. Changes
and equations for minihalo features in this version are found in
https://arxiv.org/abs/2003.04442

All functionality of the original 21cmFAST v2 C-code has been implemented in this
version, including ``USE_HALO_FIELD`` and performing full integration instead of using
the interpolation tables (which are faster).

Added
~~~~~
* Updated the radiation source model: (i) all radiation fields including X-rays, UV
  ionizing, Lyman Werner and Lyman alpha are considered from two seperated population
  namely atomic-cooling (ACGs) and minihalo-hosted molecular-cooling galaxies (MCGs);
  (ii) the turn-over masses of ACGs and MCGs are estimated with cooling efficiency and
  feedback from reionization and lyman werner suppression (Qin et al. 2020). This can
  be switched on using new ``flag_options`` ``USE_MINI_HALOS``.
* Updated kinetic temperature of the IGM with fully ionized cells following equation 6
  of McQuinn (2015) and partially ionized cells having the volume-weightied temperature
  between the ionized (volume: 1-xHI; temperature T_RE ) and neutral components (volume:
  xHI; temperature: temperature of HI). This is stored in IonizedBox as
  temp_kinetic_all_gas. Note that Tk in TsBox remains to be the kinetic temperature of HI.
* Tests: many unit tests, and also some regression tests.
* CLI: run 21cmFAST boxes from the command line, query the cache database, and produce
  plots for standard comparison runs.
* Documentation: Jupyter notebook demos and tutorials, FAQs, installation instructions.
* Plotting routines: a number of general plotting routines designed to plot coeval
  and lightcone slices.
* New power spectrum option (``POWER_SPECTRUM=5``) that uses a CLASS-based transfer
  function. WARNING: If POWER_SPECTRUM==5 the cosmo parameters cannot be altered, they
  are set to the Planck2018 best-fit values for now (until CLASS is added):
  (omegab=0.02237, omegac= 0.120, hubble=0.6736 (the rest are irrelevant for the
  transfer functions, but in case:  A_s=2.100e-9, n_s=0.9649, z_reio = 11.357)
* New ``user_params`` option ``USE_RELATIVE_VELOCITIES``, which produces initial relative
  velocity cubes (option implemented, but not the actual computation yet).
* Configuration management.
* global params now has a context manager for changing parameters temporarily.
* Vastly improved error handling: exceptions can be caught in C code and propagated to
  Python to inform the user of what's going wrong.
* Ability to write high-level data (``Coeval`` and ``Lightcone`` objects) directly to
  file in a simple portable format.

Changed
~~~~~~~
* ``POWER_SPECTRUM`` option moved from ``global_params`` to ``user_params``.
* Default cosmology updated to Planck18.

v2.0.0
------
All changes and equations for this version are found in https://arxiv.org/abs/1809.08995.

Changed
~~~~~~~

* Updated the ionizing source model: (i) the star formation rates and ionizing escape
  fraction are scaled with the masses of dark matter halos and (ii) the abundance of
  active star forming galaxies is exponentially suppressed below the turn-over halo
  mass, M_{turn}, according to a duty cycle of exp(−M_{turn}/M_{h}), where M_{h} is a
  halo mass.
* Removed the mean free path parameter, R_{mfp}. Instead, directly computes
  inhomogeneous, sub-grid recombinations in the intergalactic medium following the
  approach of Sobacchi & Mesinger (2014)




v1.2.0
------
Added
~~~~~
* Support for a halo mass dependent ionizing efficiency: zeta = zeta_0 (M/Mmin)^alpha,
  where zeta_0 corresponds to  HII_EFF_FACTOR, Mmin --> ION_M_MIN,
  alpha --> EFF_FACTOR_PL_INDEX in ANAL_PARAMS.H


v1.12.0
-------
Added
~~~~~
- Code 'redshift_interpolate_boxes.c' to interpolate between comoving cubes,
  creating comoving light cone boxes.
- Enabled openMP threading  for SMP machines.  You can specify the number of threads
  (for best performace, do not exceed the number of processors) in INIT_PARAMS.H. You do
  not need to have an SMP machine to run the code. NOTE: YOU SHOULD RE-INSTALL FFTW to
  use openMP (see INSTALL file)
- Included a threaded driver file 'drive_zscroll_reion_param.c' set-up to perform
  astrophysical parameter studies of reionization
- Included explicit support for WDM cosmologies; see COSMOLOGY.H.  The prescription is
  similar to that discussed in Barkana+2001; Mesinger+2005, madifying the (i) transfer
  function (according to the Bode+2001 formula; and (ii) including the effective
  pressure term of WDM using a Jeans mass analogy.  (ii) is approximated with a sharp
  cuttoff in the EPS barrier, using 60* M_J found in Barkana+2001 (the 60 is an
  adjustment factor found by fitting to the WDM collapsed fraction).
- A Gaussian filtering step of the PT fields to perturb_field.c, in addition to the
  implicit boxcar smoothing.  This avoids having"empty" density cells, i.e. \delta=-1,
  with some small loss in resolution.  Although for most uses \delta=-1 is ok, some Lya
  forest statistics do not like it.
- Added treatment of the risidual electron fraction from X-ray heating when computing
  the ionization field.  Relatedly, modified Ts.c to output all intermediate evolution
  boxes, Tk and x_e.
- Added a missing factor of Omega_b in Ts.c corresponding to eq. 18 in MFC11.  Users who
  used a previous version should note that their results just effecively correspond to a
  higher effective X-ray efficiency, scaled by 1/Omega_baryon.
- Normalization optimization to Ts.c, increasing performace on arge resolution boxes


Fixed
~~~~~
- GSL interpolation error in kappa_elec_pH for GSL versions > 1.15
- Typo in macro definition, which impacted the Lya background calculation in v1.11 (not applicable to earlier releases)
- Outdated filename sytax when calling gen_size_distr in drive_xHIscroll
- Redshift scrolling so that drive_logZscroll_Ts.c and Ts.c are in sync.

Changed
~~~~~~~
- Output format to avoid FFT padding for all boxes
- Filename conventions to be more explicit.
- Small changes to organization and structure


v1.1.0
------
Added
~~~~~
- Wrapper functions mod_fwrite() and mod_fread() in Cosmo_c_progs/misc.c, which
  should fix problems with the library fwrite() and fread() for large files (>4GB) on
  certain operating systems.
- Included print_power_spectrum_ICs.c program which reads in high resolution initial
  conditions and prints out an ASCII file with the associated power spectrum.
- Parameter in Ts.c for the maximum allowed kinetic temperature, which increases
  stability of the code when the redshift step size and the X-ray efficiencies are large.

Fixed
~~~~~
- Oversight adding support for a Gaussian filter for the lower resolution field.
============
Installation
============

The easiest way to install ``21cmFAST`` is to use ``conda``. Simply use
``conda install -c conda-forge 21cmFAST``. With this method, all dependencies are taken
care of, and it should work on either Linux or MacOS. If for some reason this is not
possible for you, read on.

Dependencies
------------
We try to have as many of the dependencies automatically installed as possible.
However, since ``21cmFAST`` relies on some C libraries, this is not always possible.

The C libraries required are:

* ``gsl``
* ``fftw`` (compiled with floating-point enabled, and ``--enable-shared``)
* ``openmp``
* A C-compiler with compatibility with the ``-fopenmp`` flag. **Note:** it seems that on
  OSX, if using ``gcc``, you will need ``v4.9.4+``.

As it turns out, though these are fairly common libraries, getting them installed in a
way that ``21cmFAST`` understands on various operating systems can be slightly non-trivial.

HPC
~~~
These libraries will often be available on a HPC environment by using the
``module load gsl`` and similar commands. Note that just because they are loaded
doesn't mean that ``21cmFAST`` will be able to find them. You may have to point to the
relevant ``lib/`` and ``include/`` folders for both ``gsl`` and ``fftw`` (these should
be available using ``module show gsl`` etc.)

Note also that while ``fftw`` may be available to load, it may not have the correct
compilation options (i.e. float-enabled and multiprocessing-enabled). In this case,
see below.

Linux
~~~~~
Most linux distros come with packages for the requirements, and also ``gcc`` by default,
which supports ``-fopenmp``. As long as these packages install into the standard location,
a standard installation of ``21cmFAST`` will be automatically possible (see below).
If they are installed to a place not on the ``LD_LIBRARY``/``INCLUDE`` paths, then you
must use the compilation options (see below) to specify where they are.

.. note:: there exists the option of installing ``gsl``, ``fftw`` and ``gcc`` using ``conda``.
          This is discussed below in the context of MacOSX, where it is often the
          easiest way to get the dependencies, but it is equally applicable to linux.

MacOSX
~~~~~~
On MacOSX, obtaining ``gsl`` and ``fftw`` is typically more difficult, and in addition,
the newer native ``clang`` does not offer ``-fopenmp`` support.

For ``conda`` users (which we recommend using), the easiest way to get ``gsl`` and ``fftw``
is by doing ``conda install -c conda-forge gsl fftw`` in your environment.

.. note:: if you use ``conda`` to install ``gsl`` and ``fftw``, then you will need to point at
          their location when installing `21cmFAST` (see compiler options below for details).
          In this case, the installation command should simply be *prepended* with::

              LIB=/path/to/conda/env/lib INC=/path/to/conda/env/include

To get ``gcc``, either use ``homebrew``, or again, ``conda``: ``conda install -c anaconda gcc``.
If you get the ``conda`` version, you still need to install the headers::

    xcode-select --install

On older versions then you need to do::

    open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_<input version>.pkg

.. note:: some versions of MacOS will also require you to point to the correct gcc
          compiler using the ``CC`` environment variable. Overall, the point is to NOT
          use ``clang``. If ``gcc --version`` shows that it is actually GCC, then you
          can set ``CC=gcc``. If you use homebrew to install ``gcc``, it is likely that
          you'll have to set ``CC=gcc-11``.

For newer versions, you may need to prepend the following command to your ``pip install`` command
when installing ``21cmFAST`` (see later instructions)::

    CFLAGS="-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX<input version>.sdk"

See `<faqs/installation_faq>`_ for more detailed questions on installation.
If you are on MacOSX and are having trouble with installation (or would like to share
a successful installation strategy!) please see the
`open issue <https://github.com/21cmfast/21cmFAST/issues/84>`_.

With the dependencies installed, follow the instructions below,
depending on whether you are a user or a developer.

For Users
---------

.. note:: ``conda`` users may want to pre-install the following packages before running
          the below installation commands::

            conda install numpy scipy click pyyaml cffi astropy h5py


Then, at the command line::

    pip install git+git://github.com/21cmFAST/21cmFAST.git

If developing, from the top-level directory do::

    pip install -e .

Note the compile options discussed below!

For Developers
--------------
If you are developing ``21cmFAST``, we highly recommend using ``conda`` to manage your
environment, and setting up an isolated environment. If this is the case, setting up
a full environment (with all testing and documentation dependencies) should be as easy
as (from top-level dir)::

    conda env create -f environment_dev.yml

Otherwise, if you are using ``pip``::

    pip install -e .[dev]

The ``[dev]`` "extra" here installs all development dependencies. You can instead use
``[tests]`` if you only want dependencies for testing, or ``[docs]`` to be able to
compile the documentation.

Compile Options
---------------
Various options exist to manage compilation via environment variables. Basically,
any variable with "INC" in its name will add to the includes directories, while
any variable with "lib" in its name will add to the directories searched for
libraries. To change the C compiler, use ``CC``. Finally, if you want to compile
the C-library in dev mode (so you can do stuff like valgrid and gdb with it),
install with DEBUG=True. So for example::

    CC=/usr/bin/gcc DEBUG=True GSL_LIB=/opt/local/lib FFTW_INC=/usr/local/include pip install -e .

.. note:: For MacOS a typical installation command will look like
          ``CC=gcc CFLAGS="-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX<input version>.sdk" pip install .``
          (using either ``gcc`` or ``gcc-11`` depending on how you installed gcc), with
          other compile options possible as well.

In addition, the ``BOXDIR`` variable specifies the *default* directory that any
data produced by 21cmFAST will be cached. This value can be updated at any time by
changing it in the ``$CFGDIR/config.yml`` file, and can be overwritten on a
per-call basis.

While the ``-e`` option will keep your library up-to-date with any (Python)
changes, this will *not* work when changing the C extension. If the C code
changes, you need to manually run ``rm -rf build/*`` then re-install as above.

Logging in C-Code
~~~~~~~~~~~~~~~~~
By default, the C-code will only print to stderr when it encounters warnings or
critical errors. However, there exist several levels of logging output that can be
switched on, but only at compilation time. To enable these, use the following::

    LOG_LEVEL=<log_level> pip install -e .

The ``<log_level>`` can be any non-negative integer, or one of the following
(case-insensitive) identifiers::

    NONE, ERROR, WARNING, INFO, DEBUG, SUPER_DEBUG, ULTRA_DEBUG

If an integer is passed, it corresponds to the above levels in order (starting
from zero). Be careful if the level is set to 0 (or NONE), as useful error
and warning messages will not be printed. By default, the log level is 2 (or
WARNING), unless the DEBUG=1 environment variable is set, in which case the
default is 4 (or DEBUG). Using very high levels (eg. ULTRA_DEBUG) can print out
*a lot* of information and make the run time much longer, but may be useful
in some specific cases.
=======
Authors
=======

* Brad Greig - github.com/BradGreig
* Andrei Mesinger - github.com/andreimesinger
* Steven Murray - github.com/steven-murray


Contributors
============
========
21cmFAST
========

.. start-badges
.. image:: https://travis-ci.org/21cmFAST/21cmFAST.svg
    :target: https://travis-ci.org/21cmFAST/21cmFAST
.. image:: https://coveralls.io/repos/github/21cmFAST/21cmFAST/badge.svg
    :target: https://coveralls.io/github/21cmFAST/21cmFAST
.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/ambv/black
.. image:: https://readthedocs.org/projects/21cmfast/badge/?version=latest
    :target: https://21cmfast.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
.. image:: https://img.shields.io/conda/dn/conda-forge/21cmFAST
    :target: https://github.com/conda-forge/21cmfast-feedstock
    :alt: Conda
.. image:: https://joss.theoj.org/papers/10.21105/joss.02582/status.svg
   :target: https://doi.org/10.21105/joss.02582
.. end-badges

**A semi-numerical cosmological simulation code for the radio 21-cm signal.**

.. image:: joss-paper/yuxiangs-plot-small.png
    :target: http://homepage.sns.it/mesinger/Media/lightcones_minihalo.png


This is the official repository for ``21cmFAST``: a semi-numerical code that is able to
produce 3D cosmological realisations of many physical fields in the early Universe.
It is super-fast, combining the excursion set formalism with perturbation theory to
efficiently generate density, velocity, halo, ionization, spin temperature, 21-cm, and
even ionizing flux fields (see the above lightcones!).
It has been tested extensively against numerical simulations, with excellent agreement
at the relevant scales.

``21cmFAST`` has been widely used, for example, by the Murchison Widefield Array (MWA),
LOw-Frequency ARray (LOFAR) and Hydrogen Epoch of Reionization Array (HERA), to model the
large-scale cosmological 21-cm signal. In particular, the speed of ``21cmFAST`` is important
to produce simulations that are large enough (several Gpc across) to represent modern
low-frequency observations.

As of ``v3.0.0``, ``21cmFAST`` is conveniently wrapped in Python to enable more dynamic code.


New Features in 3.0.0+
======================

* Robust on-disk caching/writing both for efficiency and simplified reading of
  previously processed data (using HDF5).
* Convenient data objects which simplify access to and processing of the various density
  and ionization fields.
* De-coupled functions mean that arbitrary functionality can be injected into the process.
* Improved exception handling and debugging
* Comprehensive documentation
* Comprehensive test suite.
* Strict `semantic versioning <https://semver.org>`_.

Installation
============
We support Linux and MacOS (please let us know if you are successful in installing on
Windows!). On these systems, the simplest way to get ``21cmFAST`` is by using
`conda <https://www.anaconda.com/>`_::

    conda install -c conda-forge 21cmFAST

``21cmFAST`` is also available on PyPI, so that ``pip install 21cmFAST`` also works. However,
it depends on some external (non-python) libraries that may not be present, and so this
method is discouraged unless absolutely necessary. If using ``pip`` to install ``21cmFAST``
(especially on MacOS), we thoroughly recommend reading the detailed
`installation instructions <https://21cmfast.readthedocs.io/en/latest/installation.html>`_.

Basic Usage
===========
``21cmFAST`` can be run both interactively and from the command line (CLI).

Interactive
-----------
The most basic example of running a (very small) coeval simulation at a given redshift,
and plotting an image of a slice through it::

    >>> import py21cmfast as p21c
    >>> coeval = p21c.run_coeval(
    >>>     redshift=8.0,
    >>>     user_params={'HII_DIM': 50, "USE_INTERPOLATION_TABLES": False}
    >>> )
    >>> p21c.plotting.coeval_sliceplot(coeval, kind='brightness_temp')

The coeval object here has much more than just the ``brightness_temp`` field in it. You
can plot the ``density`` field, ``velocity`` field or a number of other fields.
To simulate a full lightcone::

    >>> lc = p21c.run_lightcone(
    >>>     redshift=8.0,
    >>>     max_redshift=15.0,
    >>>     init_box = coeval.init_struct,
    >>> )
    >>> p21c.plotting.lightcone_sliceplot(lc)

Here, we used the already-computed initial density field from ``coeval``, which sets
the size and parameters of the run, but also means we don't have to compute that
(relatively expensive step again). Explore the full range of functionality in the
`API Docs <https://21cmfast.readthedocs.io/en/latest/reference/py21cmfast.html>`_,
or read more `in-depth tutorials <https://21cmfast.readthedocs.io/en/latest/tutorials.html>`_
for further guidance.

CLI
---
The CLI can be used to generate boxes on-disk directly from a configuration file or
command-line parameters. You can run specific steps of the simulation independently,
or an entire simulation at once. For example, to run just the initial density field,
you can do::

    $ 21cmfast init --HII_DIM=100

The (quite small) simulation box produced is automatically saved into the cache
(by default, at ``~/21cmFAST-cache``).
You can list all the files in your cache (and the parameters used in each of the simulations)
with::

    $ 21cmfast query

To run an entire coeval cube, use the following as an example::

    $ 21cmfast coeval 8.0 --out=output/coeval.h5 --HII_DIM=100

In this case all the intermediate steps are cached in the standard cache directory, and
the final ``Coeval`` box is saved to ``output/coeval.h5``. If no ``--out`` is specified,
the coeval box itself is not written, but don't worry -- all of its parts are cached, and
so it can be rebuilt extremely quickly. Every input parameter to any of the
`input classes <https://21cmfast.readthedocs.io/en/latest/reference/_autosummary/py21cmfast.inputs.html>`_
(there are a lot of parameters) can be specified at the end of the call with prefixes of
``--`` (like ``HII_DIM`` here). Alternatively, you can point to a config YAML file, eg.::

    $ 21cmfast lightcone 8.0 --max-z=15.0 --out=. --config=~/.21cmfast/runconfig_example.yml

There is an example configuration file `here <user_data/runconfig_example.yml>`_ that you
can build from. All input parameters are
`documented here <https://21cmfast.readthedocs.io/en/latest/reference/_autosummary/py21cmfast.inputs.html>`_.

Documentation
=============
Full documentation (with examples, installation instructions and full API reference)
found at https://21cmfast.readthedocs.org.

Acknowledging
=============
If you use ``21cmFAST v3+`` in your research please cite both of:

    Murray et al., (2020). 21cmFAST v3: A Python-integrated C code for generating 3D
    realizations of the cosmic 21cm signal. Journal of Open Source Software, 5(54),
    2582, https://doi.org/10.21105/joss.02582

    Andrei Mesinger, Steven Furlanetto and Renyue Cen, "21CMFAST: a fast, seminumerical
    simulation of the high-redshift 21-cm signal", Monthly Notices of the Royal
    Astronomical Society, Volume 411, Issue 2, pp. 955-972 (2011),
    https://ui.adsabs.harvard.edu/link_gateway/2011MNRAS.411..955M/doi:10.1111/j.1365-2966.2010.17731.x

In addition, the following papers introduce various features into ``21cmFAST``. If you use
these features, please cite the relevant papers.

Mini-halos:

    Muñoz, J.B., Qin, Y., Mesinger, A., Murray, S., Greig, B., and Mason, C.,
    "The Impact of the First Galaxies on Cosmic Dawn and Reionization"
    https://arxiv.org/abs/2110.13919
    (for DM-baryon relative velocities)

    Qin, Y., Mesinger, A., Park, J., Greig, B., and Muñoz, J. B.,
    “A tale of two sites - I. Inferring the properties of minihalo-hosted galaxies from
    current observations”, Monthly Notices of the Royal Astronomical Society, vol. 495,
    no. 1, pp. 123–140, 2020. https://doi.org/10.1093/mnras/staa1131.
    (for Lyman-Werner and first implementation)

Mass-dependent ionizing efficiency:

    Park, J., Mesinger, A., Greig, B., and Gillet, N.,
    “Inferring the astrophysics of reionization and cosmic dawn from galaxy luminosity
    functions and the 21-cm signal”, Monthly Notices of the Royal Astronomical Society,
    vol. 484, no. 1, pp. 933–949, 2019. https://doi.org/10.1093/mnras/stz032.
.. include:: ../CONTRIBUTING.rst
.. include:: notes_for_developers.rst
.. include:: ../AUTHORS.rst
Developer Documentation
=======================

If you are new to developing ``21cmFAST``, please read the :ref:`contributing:Contributing`
section *first*, which outlines the general concepts for contributing to development,
and provides a step-by-step walkthrough for getting setup.
This page lists some more detailed notes which may be helpful through the
development process.

Compiling for debugging
-----------------------
When developing, it is usually a good idea to compile the underlying C code in ``DEBUG``
mode. This may allow extra print statements in the C, but also will allow running the C
under ``valgrind`` or ``gdb``. To do this::

    $ DEBUG=True pip install -e .

See :ref:`installation:Installation` for more installation options.

Developing the C Code
---------------------
In this section we outline how one might go about modifying/extending the C code and
ensuring that the extension is compatible with the wrapper provided here. It is
critical that you run all tests after modifying _anything_ (and see the section
below about running with valgrind). When changing C code, before
testing, ensure that the new C code is compiled into your environment by running::

    $ rm -rf build
    $ pip install .

Note that using a developer install (`-e`) is not recommended as it stores compiled
objects in the working directory which don't get updated as you change code, and can
cause problems later.

There are two main purposes you may want to write some C code:

1. An external plugin/extension which uses the output data from 21cmFAST.
2. Modifying the internal C code of 21cmFAST.

21cmFAST currently provides no support for external plugins/extensions. It is entirely
possible to write your own C code to do whatever you want with the output data, but we
don't provide any wrapping structure for you to do this, you will need to write your
own. Internally, 21cmFAST uses the ``cffi`` library to aid the wrapping of the C code into
Python. You don't need to do the same, though we recommend it. If your desired
"extension" is something that needs to operate in-between steps of 21cmFAST, we also
provide no support for this, but it is possible, so long as the next step in the
chain maintains its API. You would be required to re-write the low-level wrapping
function _preceding_ your inserted step as well. For instance, if you had written a
self-contained piece of code that modified the initial conditions box, adding some
physical effect which is not already covered, then you would need to write a low-level
wrapper _and_ re-write the ``initial_conditions`` function to modify the box before
returning it. We provide no easy "plugin" system for doing this currently. If your
external code is meant to be inserted _within_ a basic step of 21cmFAST, this is
currently not possible. You will instead have to modify the source code itself.

Modifying the C-code of 21cmFAST should be relatively simple. If your changes are
entirely internal to a given function, then nothing extra needs to be done. A little
more work has to be done if the modifications add/remove input parameters or the output
structure. If any of the input structures are modified (i.e. an extra parameter
added to it), then the corresponding class in ``py21cmfast.wrapper`` must be modified,
usually simply to add the new parameter to the ``_defaults_`` dict with a default value.
For instance, if a new variable ``some_param`` was added to the ``user_params`` struct
in the ``ComputeInitialConditions`` C function, then the ``UserParams`` class in
the wrapper would be modified, adding ``some_param=<default_value>`` to its ``_default_``
dict. If the default value of the parameter is dependent on another parameter, its
default value in this dict can be set to ``None``, and you can give it a dynamic
definition as a Python ``@property``. For example, the ``DIM`` parameter of
``UserParams`` is defined as::

    @property
    def DIM(self):
        if self._some_param is None:
            return self._DIM or 4 * self.HII_DIM

Note the underscore in ``_DIM`` here: by default, if a dynamic property is defined for
a given parameter, the ``_default_`` value is saved with a prefixed underscore. Here we
return either the explicitly set ``DIM``, or 4 by the ``HII_DIM``. In addition, if the
new parameter is not settable -- if it is completely determined by other parameters --
then don't put it in ``_defaults_`` at all, and just give it a dynamic definition.

If you modify an output struct, which usually house a number of array quantities
(often float pointers, but not necessarily), then you'll again need to modify the
corresponding class in the wrapper. In particular, you'll need to add an entry for that
particular array in the ``_init_arrays`` method for the class. The entry consists of
initialising that array (usually to zeros, but not necessarily), and setting its proper
dtype. All arrays should be single-pointers, even for multi-dimensional data. The latter
can be handled by initalising the array as a 1D numpy array, but then setting its shape
attribute (after creation) to the appropriate n-dimensional shape (see the
``_init_arrays`` method for the ``InitialConditions`` class for examples of this).

Modifying the ``global_params`` struct should be relatively straightforward, and no
changes in the Python are necessary. However, you may want to consider adding the new
parameter to relevant ``_filter_params`` lists for the output struct wrapping classes in
the wrapper. These lists control which global parameters affect which output structs,
and merely provide for more accurate caching mechanisms.

C Function Standards
~~~~~~~~~~~~~~~~~~~~
The C-level functions are split into two groups -- low-level "private" functions, and
higher-level "public" or "API" functions. All API-level functions are callable from
python (but may also be called from other C functions). All API-level functions are
currently prototyped in ``21cmFAST.h``.

To enable consistency of error-checking in Python (and a reasonable standard for any
kind of code), we enforce that any API-level function must return an integer status.
Any "return" objects must be modified in-place (i.e. passed as pointers). This enables
Python to control the memory access of these variables, and also to receive proper
error statuses (see below for how we do exception handling). We also adhere to the
convention that "output" variables should be passed to the function as its last
argument(s). In the case that _only_ the last argument is meant to be "output", there
exists a simple wrapper ``_call_c_simple`` in ``wrapper.py`` that will neatly handle the
calling of the function in an intuitive pythonic way.

Running with Valgrind
~~~~~~~~~~~~~~~~~~~~~
If any changes to the C code are made, it is ideal to run tests under valgrind, and
check for memory leaks. To do this, install ``valgrind`` (we have tested v3.14+),
which is probably available via your package manager. We provide a
suppression file for ``valgrind`` in the ``devel/`` directory of the main repository.

It is ideal if you install a development-version of python especially for running these
tests. To do this, download the version of python you want and then configure/install with::

    $ ./configure --prefix=<your-home>/<directory> --without-pymalloc --with-pydebug --with-valgrind
    $ make; make install

Construct a ``virtualenv`` on top of this installation, and create your environment,
and install all requirements.

If you do not wish to run with a modified version of python, you may continue with your
usual version, but may get some extra cruft in the output. If running with Python
version > 3.6, consider running with environment variable ``PYTHONMALLOC=malloc``
(see https://stackoverflow.com/questions/20112989/how-to-use-valgrind-with-python ).

The general pattern for using valgrind with python is::

    $ valgrind --tool=memcheck --track-origins=yes --leak-check=full --suppressions=devel/valgrind-suppress-all-but-c.supp <python script>

One useful command is to run valgrind over the test suite (from the top-level repo
directory)::

    $ valgrind --tool=memcheck --track-origins=yes --leak-check=full --suppressions=devel/valgrind-suppress-all-but-c.supp pytest

While we will attempt to keep the suppression file updated to the best of our knowledge
so that only relevant leaks and errors are reported, you will likely have to do a bit of
digging to find the relevant parts.

Valgrind will likely run very slowly, and sometimes  you will know already which exact
tests are those which may have problems, or are relevant to your particular changes.
To run these::

    $ PYTHONMALLOC=malloc valgrind --tool=memcheck --track-origins=yes --leak-check=full --suppressions=devel/valgrind-suppress-all-but-c.supp pytest -v tests/<test_file>::<test_func> > valgrind.out 2>&1

Note that we also routed the stderr output to a file, which is useful because it can be
quite voluminous. There is a python script, ``devel/filter_valgrind.py`` which can be run
over the output (`valgrind.out` in the above command) to filter it down to only have
stuff from 21cmfast in it.

Producing Integration Test Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
There are bunch of so-called "integration tests", which rely on previously-produced
data. To produce this data, run ``python tests/produce_integration_test_data.py``.

Furthermore, this data should only be produced with good reason -- the idea is to keep
it static while the code changes, to have something steady to compare to. If a particular
PR fixes a bug which affects a certain tests' data, then that data should be re-run, in
the context of the PR, so it can be explained.

Logging in C
~~~~~~~~~~~~
The C code has a header file ``logging.h``. The C code should *never* contain bare
print-statements -- everything should be formally logged, so that the different levels
can be printed to screen correctly. The levels are defined in ``logging.h``, and include
levels such as ``INFO``, ``WARNING`` and ``DEBUG``. Each level has a corresponding macro
that starts with ``LOG_``. Thus to log run-time information to stdout, you would use
``LOG_INFO("message");``. Note that the message does not require a final newline character.

Exception handling in C
~~~~~~~~~~~~~~~~~~~~~~~
There are various places that things can go wrong in the C code, and they need to be
handled gracefully so that Python knows what to do with it (rather than just quitting!).
We use the simple ``cexcept.h`` header file from http://www.nicemice.net/cexcept/ to
enable a simple form of exception handling. That file itself should **not be edited**.
There is another header -- ``exceptions.h`` -- that defines how we use exceptions
throughout ``21cmFAST``. Any time an error arises that can be understood, the developer
should add a ``Throw <ErrorKind>;`` line. The ``ErrorKind`` can be any of the kinds
defined in ``exceptions.h`` (eg. ``GSLError`` or ``ValueError``). These are just integers.

Any C function that has a header in ``21cmFAST.h`` -- i.e. any function that is callable
directly from Python -- *must* be globally wrapped in a ``Try {} Catch(error_code) {}`` block. See
``GenerateICs.c`` for an example. Most of the code should be in the ``Try`` block.
Anything that does a ``Throw`` at any level of the call stack within that ``Try`` will
trigger a jump to the ``Catch``. The ``error_code`` is the integer that was thrown.
Typically, one will perhaps want to do some cleanup here, and then finally *return* the
error code.

Python knows about the exit codes it can expect to receive, and will raise Python
exceptions accordingly. From the python side, two main kinds of exceptions could be
raised, depending on the error code returned from C. The lesser exception is called a
``ParameterError``, and is supposed to indicate an error that happened merely because
the parameters that were input to the calculation were just too extreme to handle.
In the case of something like an automatic Monte Carlo algorithm that's iterating over
random parameters, one would *usually* want to just keep going at this point, because
perhaps it just wandered too far in parameter space.
The other kind of error is a ``FatalCError``, and this is where things went truly wrong,
and probably will do for any combination of parameters.

If you add a kind of Exception in the C code (to ``exceptions.h``), then be sure to add
a handler for it in the ``_process_exitcode`` function in ``wrapper.py``.


Maintaining Array State
~~~~~~~~~~~~~~~~~~~~~~~
Part of the challenge of maintaining a nice wrapper around the fast C-code is keeping
track of initialized memory, and ensuring that the C structures that require that memory
are pointing to the right place. Most of the arrays that are computed in ``21cmFAST``
are initialized *in Python* (using Numpy), then a pointer to their memory is given to
the C wrapper object.

To make matters more complicated, since some of the arrays are really big, it is sometimes
necessary to write them to disk to relieve memory pressure, and load them back in as required.
That means that any time, a given array in a C-based class may have one of several different "states":

1. Completely Uninitialized
1. Allocated an initialized in memory
1. Computed (i.e. filled with the values defining that array after computation in C)
1. Stored on disk
1. Stored *and* in memory.

It's important to keep track of these states, because when passing the struct to the ``compute()``
function of another struct (as input), we go and check if the array exists in memory, and
initialize it. Of course, we shouldn't initialize it with zeros if in fact it has been computed already
and is sitting on disk ready to be loaded. Thus, the ``OutputStruct`` tries to keep track of these
states for every array in the structure, using the ``_array_state`` dictionary. Every write/read/compute/purge
operation self-consistently modifies the status of the array.

However, one needs to be careful -- you *can* modify the actual state without modifying the ``_array_state``
(eg. simply by doing a ``del object.array``). In the future, we may be able to protect this to some extent,
but for now we rely on the good intent of the user.

Purging/Loading C-arrays to/from Disk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
As of v3.1.0, there are more options for granular I/O, allowing large arrays to be purged from memory
when they are unnecessary for further computation. As a developer, you should be aware of the ``_get_required_input_arrays``
method on all ``OutputStruct`` subclasses. This is available to tell the given class what arrays need to
be available at compute time in any of the input structs. For example, if doing ``PERTURB_ON_HIGH_RES``,
the ``PerturbedField`` requires the hi-res density fields in ``InitialConditions``. This gives indications
as to what boxes can be purged to disk (all the low-res boxes in the ICs, for example).
Currently, this is only used to *check* that all boxes are available at compute time, and is not used
to actually automatically purge anything. Note however that ``InitialConditions`` does have two
custom methods that will purge unnecessary arrays before computing perturb fields or ionization fields.

.. note:: If you add a new quantity to a struct, and it is required input for other structs, you need
          to add it to the relevant ``_get_required_input_arrays`` methods.

Further note that as of v3.1.0, partial structs can be written and read from disk (so you can specify
``keys=['hires_density']`` in the ``.read()`` method to just read the hi-res density field into the object.



Branching and Releasing
-----------------------
The aim is to make 21cmFAST's releases as useful, comprehendible, and automatic
as possible. This section lays out explicitly how this works (mostly for the benefit of
the admin(s)).

Versioning
~~~~~~~~~~
The first thing to mention is that we use strict `semantic versioning <https://semver.org>`_
(since v2.0). Thus the versions are ``MAJOR.MINOR.PATCH``, with ``MAJOR`` including
API-breaking changes, ``MINOR`` including new features, and ``PATCH`` fixing bugs or
documentation etc. If you depend on hmf, you can set your dependency as
``21cmFAST >= X.Y < X+1`` and not worry that we'll break your code with an update.

To mechanically handle versioning within the package, we use two methods that we make
to work together automatically. The "true" version of the package is set with
`setuptools-scm <https://pypi.org/project/setuptools-scm/>`_. This stores the version
in the git tag. There are many benefits to this -- one is that the version is unique
for every single change in the code, with commits on top of a release changing the
version. This means that versions accessed via ``py21cmfast.__version__`` are unique and track
the exact code in the package (useful for reproducing results). To get the current
version from command line, simply do ``python setup.py --version`` in the top-level
directory.

To actually bump the version, we use ``bump2version``. The reason for this is that the
CHANGELOG requires manual intervention -- we need to change the "dev-version" section
at the top of the file to the current version. Since this has to be manual, it requires
a specific commit to make it happen, which thus requires a PR (since commits can't be
pushed to master). To get all this to happen as smoothly as possible, we have a little
bash script ``bump`` that should be used to bump the version, which wraps ``bump2version``.
What it does is:

1. Runs ``bump2version`` and updates the ``major``, ``minor`` or ``patch`` part (passed like
   ``./bump minor``) in the VERSION file.
2. Updates the changelog with the new version heading (with the date),
   and adds a new ``dev-version`` heading above that.
3. Makes a commit with the changes.

.. note:: Using the ``bump`` script is currently necessary, but future versions of
   ``bump2version`` may be able to do this automatically, see
   https://github.com/c4urself/bump2version/issues/133.

The VERSION file might seem a bit redundant, and it is NOT recognized as the "official"
version (that is given by the git tag). Notice we didn't make a git tag in the above
script. That's because the tag should be made directly on the merge commit into master.
We do this using a Github Action (``tag-release.yaml``) which runs on every push to master,
reads the VERSION file, and makes a tag based on that version.


Branching
~~~~~~~~~
For branching, we use a very similar model to `git-flow <https://nvie.com/posts/a-successful-git-branching-model/>`_.
That is, we have a ``master`` branch which acts as the current truth against which to develop,
and ``production`` essentially as a deployment branch.
I.e., the ``master`` branch is where all features are merged (and some
non-urgent bugfixes). ``production`` is always production-ready, and corresponds
to a particular version on PyPI. Features should be branched from ``master``,
and merged back to ``production``. Hotfixes can be branched directly from ``production``,
and merged back there directly, *as well as* back into ``master``.
*Breaking changes* must only be merged to ``master`` when it has been decided that the next
version will be a major version. We do not do any long-term support of releases
(so can't make hotfixes to ``v2.x`` when the latest version is ``2.(x+1)``, or make a
new minor version in 2.x when the latest version is 3.x). We have set the default
branch to ``dev`` so that by default, branches are merged there. This is deemed best
for other developers (not maintainers/admins) to get involved, so the default thing is
usually right.

.. note:: Why not a more simple workflow like Github flow? The simple answer is it just
          doesn't really make sense for a library with semantic versioning. You get into
          trouble straight away if you want to merge a feature but don't want to update
          the version number yet (you want to merge multiple features into a nice release).
          In practice, this happens quite a lot.

.. note:: OK then, why not just use ``production`` to accrue features and fixes until such
          time we're ready to release? The problem here is that if you've merged a few
          features into master, but then realize a patch fix is required, there's no
          easy way to release that patch without releasing all the merged features, thus
          updating the minor version of the code (which may not be desirable). You could
          then just keep all features in their own branches until you're ready to release,
          but this is super annoying, and doesn't give you the chance to see how they
          interact.


Releases
~~~~~~~~
To make a **patch** release, follow these steps:

1. Branch off of ``production``.
2. Write the fix.
3. Write a test that would have broken without the fix.
4. Update the changelog with your changes, under the ``**Bugfixes**`` heading.
5. Commit, push, and create a PR.
6. Locally, run ``./bump patch``.
7. Push.
8. Get a PR review and ensure CI passes.
9. Merge the PR

Note that in the background, Github Actions *should* take care of then tagging ``production``
with the new version, deploying that to PyPI, creating a new PR from master back into
``master``, and accepting that PR. If it fails for one of these steps, they can all be done
manually.

Note that you don't have to merge fixes in this way. You can instead just branch off
``master``, but then the fix won't be included until the next ``minor`` version.
This is easier (the admins do the adminy work) and useful for non-urgent fixes.

Any other fix/feature should be branched from ``master``. Every PR that does anything
noteworthy should have an accompanying edit to the changelog. However, you do not have
to update the version in the changelog -- that is left up to the admin(s). To make a
minor release, they should:

1. Locally, ``git checkout release``
2. ``git merge master``
3. No new features should be merged into ``master`` after that branching occurs.
4. Run ``./bump minor``
5. Make sure everything looks right.
6. ``git push``
7. Ensure all tests pass and get a CI review.
8. Merge into ``production``

The above also works for ``MAJOR`` versions, however getting them *in* to ``master`` is a little
different, in that they should wait for merging until we're sure that the next version
will be a major version.
Tutorials and FAQs
==================

The following introductory tutorials will help you get started with ``21cmFAST``:

.. toctree::
   :maxdepth: 2

   tutorials/coeval_cubes
   tutorials/lightcones
   tutorials/mini-halos
   tutorials/gather_data
   tutorials/relative_velocities

If you've covered the tutorials and still have questions about "how to do stuff" in
``21cmFAST``, consult the FAQs:

.. toctree::
   :maxdepth: 2

   faqs/installation_faq
   faqs/misc
.. include:: ../README.rst
.. include:: ../INSTALLATION.rst
======================================
Design Philosophy and Features for v3+
======================================

Here we describe in broad terms the design philosophy of the *new* ``21cmFAST``,
and some of its new features.
This is useful to get an initial bearing of how to go about using ``21cmFAST``, though
most likely the :doc:`tutorials <tutorials>` will be better for that.
It is also useful for those who have used the "old" ``21cmFAST`` (versions 2.1 and less)
and want to know why they should use this new version (and how to convert).
In doing so, we'll go over some of the key features of ``21cmFAST`` v3+.
To get a more in-depth view of all the options and features available, look at the
very thorough :doc:`API Reference <reference/index>`.


Design Philosophy
=================
The goal of v3 of ``21cmFAST`` is to provide the same computational efficiency and
scientific features of the previous generations, but packaged in a form that adopts the
best modern programming standards, including:

* simple installation
* comprehensive documentation
* comprehensive test suite
* more modular code
* standardised code formatting
* truly open-source and collaborative design (via Github)

Partly to enable these standards, and partly due to the many *extra* benefits it brings,
v3 also has the major distinction of being wrapped entirely in Python. The *extra*
benefits brought by this include:

* a native python library interface (eg. get your output box directly as a ``numpy`` array).
* better file-writing, into the HDF5 format, which saves metadata along with the box data.
* a caching system so that the same data never has to be calculated twice.
* reproducibility: know which exact version of ``21cmFAST``, with what parameters, produced a given dataset.
* significantly improved warnings and error checking/propagation.
* simplicity for adding new additional effects, and inserting them in the calculation pipeline.

We hope that additional features and benefits will be found by the growing community
of ``21cmFAST`` developers and users.

How it Works
============
v3 is *not* a complete rewrite of ``21cmFAST``. Most of the C-code of previous versions
is kept, though it has been modularised and modified in many places. The fundamental
routines are the same (barring bugfixes!).

The major programs of the original version (``init``, ``perturb``, ``ionize`` etc.) have
been converted into modular *functions* in C. Furthermore, most of the global parameters
(and, more often than not, global ``#define`` options) have been modularised and converted
into a series of input "parameter" ``structs``. These get passed into the functions.
Furthermore, each C function, instead of writing a bunch of files, returns an output
``struct`` containing all the stuff it computed.

Each of these functions and structs are wrapped in Python using the ``cffi`` package.
CFFI compiles the C code once upon *installation*. Due to the fact that parameters are
now passed around to the different functions, rather than being global defines, we no
longer need to re-compile every time an option is changed. Python itself can handle
changing the parameters, and can use the outputs in whatever way the user desires.

To maintain continuity with previous versions, a CLI interface is provided (see below)
that acts in a similar fashion to previous versions.

When ``21cmFAST`` is installed, it automatically creates a configuration directory in
the user's home: ``~/.21cmfast``. This houses a number of important configuration
options; usually default values of parameters. At this stage, the location of this
directory is not itself configurable. The config directory contains example
configuration files for the CLI interface (see below), which can also be copied anywhere
on disk and modified. Importantly, the ``config.yml`` file in this directory specifies
some of the more important options for how ``21cmFAST`` behaves by default.
One such option is ``boxdir``, which specifies the directory in which ``21cmFAST`` will
cache results (see below for details). Finally, the config directory houses several data
tables which are used to accelerate several calculations. In principle
these files are over-writeable, but they should only be touched if one knows very well
what they are doing.

Finally, ``21cmFAST`` contains a more robust cataloguing/caching method. Instead of
saving data with a selection of the dependent parameters written into the filename --
a method which is prone to error if a parameter which is not part of that selection is
modified -- ``21cmFAST`` writes all data into a configurable central directory with a hash
filename unique to *all* parameters upon which the data depends. Each kind of dataset has
attached methods which efficiently search this central directory for matching data to be
read when necessary.
Several arguments are available for all library functions which produce such datasets
that control this output. In this way, the data that is being retrieved is always
reliably produced with the desired parameters, and users need not concern themselves
with how and where the data is saved -- it can be retrieved merely by creating an empty
object with the desired parameters and calling ``.read()``, or even better, by calling
the function to *produce* the given dataset, which will by default just read it in if
available.

CLI
===
The CLI interface always starts with the command ``21cmfast``, and has a number of
subcommands. To list the available subcommands, use::

    $ 21cmfast --help

To get help on any subcommand, simply use::

    $ 21cmfast <subcommand> --help

Any subcommand which runs some aspect of ``21cmFAST`` will have a ``--config`` option,
which specifies a configuration file (by default this is
``~/.21cmfast/runconfig_example.yml``). This config file specifies the parameters of the
run. Furthermore, any particular parameter that can be specified in the config file can
be alternatively specified on the command line by appending the command with the
parameter name, eg.::

    $ 21cmfast init --config=my_config.yml --HII_DIM=40 hlittle=0.7 --DIM 100 SIGMA_8 0.9

The above command shows the numerous ways in which these parameters can be specified
(with or without leading dashes, and with or without "=").

The CLI interface, while simple to use, does have the limitation that more complex
arguments than can be passed to the library functions are disallowed. For example,
one cannot pass a previously calculated initial conditions box to the ``perturb``
command. However, if such a box has been calculated with the default option to write it
to disk, then it will automatically be found and used in such a situation, i.e. the
following will not re-calculate the init box::

    $ 21cmfast init
    $ 21cmfast perturb redshift=8.0

This means that almost all of the functionality provided in the library is accessible
via the CLI.
.. include:: readme.rst

Contents
========
.. toctree::
   :maxdepth: 2

   installation
   design
   tutorials
   reference/index
   contributing
   authors
   changelog

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
.. include:: ../CHANGELOG.rst
Miscellaneous FAQs
==================

My run seg-faulted, what should I do?
-------------------------------------
Since ``21cmFAST`` is written in C, there is the off-chance that something
catastrophic will happen, causing a segfault. Typically, if this happens, Python will
not print a traceback where the error occurred, and finding the source of such errors
can be difficult. However, one has the option of using the standard library
`faulthandler <https://docs.python.org/3/library/faulthandler.html>`_. Specifying
``-X faulthandler`` when invoking Python will cause a minimal traceback to be printed
to ``stderr`` if a segfault occurs.

Configuring 21cmFAST
--------------------
``21cmFAST`` has a configuration file located at ``~/.21cmfast/config.yml``. This file
specifies some options to use in a rather global sense, especially to do with I/O.
You can directly edit this file to change how ``21cmFAST`` behaves for you across
sessions.
For any particular function call, any of the options may be overwritten by supplying
arguments to the function itself.
To set the configuration for a particular session, you can also set the global ``config``
instance, for example::

    >>> import py21cmfast as p21
    >>> p21.config['regenerate'] = True
    >>> p21.run_lightcone(...)

All functions that use the ``regenerate`` keyword will now use the value you've set in the
config. Sometimes, you may want to be a little more careful -- perhaps you want to change
the configuration for a set of calls, but have it change back to the defaults after that.
We provide a context manager to do this::

    >>> with p21.config.use(regenerate=True):
    >>>     p21.run_lightcone()
    >>>     print(p21.config['regenerate'])  # prints "True"
    >>> print(p21.config['regenerate'])  # prints "False"

To make the current configuration permanent, simply use the ``write`` method::

    >>> p21.config['direc'] = 'my_own_cache'
    >>> p21.config.write()

Global Parameters
-----------------
There are a bunch of "global" parameters that are used throughout the C code. These are
parameters that are deemed to be constant enough to not expose them through the
regularly-used input structs, but nevertheless may necessitate modification from
time-to-time. These are accessed through the ``global_params`` object::

    >>> from py21cmfast import global_params

Help on the attributes can be obtained via ``help(global_params)`` or
`in the docs <../reference/_autosummary/py21cmfast.inputs.html>`_. Setting the
attributes (which affects them everywhere throughout the code) is as simple as, eg::

    >>> global_params.Z_HEAT_MAX = 30.0

If you wish to use a certain parameter for a fixed portion of your code (eg. for a single
run), it is encouraged to use the context manager, eg.::

    >>> with global_params.use(Z_HEAT_MAX=10):
    >>>    run_lightcone(...)
Installation FAQ
================

Errors with "recompile with -fPIC" for FFTW
-------------------------------------------
Make sure you have installed FFTW with ``--enable-shared``. On Ubuntu, you will
also need to have ``libfftw3-dev`` installed.
API Reference
=============

.. toctree::
    :glob:
    :maxdepth: 3

    py21cmfast*
py21cmfast
==========

.. testsetup::

    from py21cmfast import *

.. autosummary::
    :toctree: _autosummary
    :template: module.rst

    py21cmfast.inputs
    py21cmfast.outputs
    py21cmfast.wrapper
    py21cmfast.plotting
    py21cmfast.cache_tools
{{ fullname | escape | underline }}

.. automodule:: {{ fullname }}

   {% block functions %}
   {% if functions %}
   .. rubric:: Functions

   .. autosummary::
      :toctree: {{ objname }}
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: Classes

   .. autosummary::
      :toctree: {{ objname }}
      :template: class.rst
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   .. rubric:: Exceptions

   .. autosummary::
   {% for item in exceptions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
{{ fullname | escape | underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}

   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
      :toctree: {{ objname }}
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
      :toctree: {{ objname }}
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
