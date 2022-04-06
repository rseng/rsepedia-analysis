<p align="center">
    <img width="500", src="https://raw.githubusercontent.com/TeamLEGWORK/LEGWORK/main/docs/images/legwork.png">
</p>

<h2 align="center">
    The <b>L</b>ISA <b>E</b>volution and <b>G</b>ravitational <b>W</b>ave <b>OR</b>bit <b>K</b>it
    <br>
    <a href="https://github.com/TeamLEGWORK/LEGWORK-paper">
        <img src="https://img.shields.io/badge/release paper-repo-blue.svg?style=flat&logo=GitHub" alt="Read the article"/>
    </a>
    <a href="https://codecov.io/gh/TeamLEGWORK/LEGWORK">
        <img src="https://codecov.io/gh/TeamLEGWORK/LEGWORK/branch/main/graph/badge.svg?token=FUG4RFYCWX"/>
    </a>
    <a href='https://legwork.readthedocs.io/en/latest/?badge=latest'>
        <img src='https://readthedocs.org/projects/legwork/badge/?version=latest' alt='Documentation Status' />
    </a>
    <a href="https://ascl.net/2111.007">
        <img src="https://img.shields.io/badge/ascl-2111.007-blue.svg?colorB=262255" alt="ascl:2111.007" />
    </a>
    <a href="mailto:tomjwagg@gmail.com?cc=kbreivik@flatironinstitute.org">
        <img src="https://img.shields.io/badge/contact-authors-blueviolet.svg?style=flat" alt="Email the authors"/>
    </a>
</h2>

<p align="center">
    A python package that does the <code>LEGWORK</code> for you by evolving binaries,
    computing signal-to-noise ratios for binary systems potentially observable with LISA
    and visualising the results.
</p>

### Installation
Put simply? `pip install legwork`! But we recommend creating a conda environment first to ensure everything goes smoothly! Check out the installation instructions [here](https://legwork.readthedocs.io/en/latest/install.html) to learn exactly how to install LEGWORK

LEGWORK has a couple of dependencies: `numpy`, `astropy`, `numba`, `scipy`, `schwimmbad`, `matplotlib` and `seaborn` (see [requirements.txt](requirements.txt) for the exact version requirements). These will be installed automatically if you follow the installation instructions!

### Documentation
All documentation related to LEGWORK can be found [at this link](https://legwork.readthedocs.io/en/latest/)

### Other quick links
- [Quickstart](https://legwork.readthedocs.io/en/latest/notebooks/Quickstart.html) - New to LEGWORK? Try out our quickstart tutorial!
- [Tutorials](https://legwork.readthedocs.io/en/latest/tutorials.html) - Learn more about what you can do with LEGWORK with our tutorials!
- [Citing LEGWORK](https://legwork.readthedocs.io/en/latest/cite.html) - If you're using LEGWORK for a scientific publication please follow the link for citation intstructions
- [Demos](https://legwork.readthedocs.io/en/latest/demos.html) - Want to see what LEGWORK is capable of? Check out our demos!
- [API reference](https://legwork.readthedocs.io/en/latest/modules.html) - Wondering how you should use a particular function? Go take a look at our full API reference!
- [Feature requests](https://github.com/TeamLEGWORK/LEGWORK/issues/new) - Do you have an idea for adding something to LEGWORK? Create an issue [here](https://github.com/TeamLEGWORK/LEGWORK/issues/new) and let us know! Or, even better, make the change yourself and create a [pull request](https://github.com/TeamLEGWORK/LEGWORK/pulls)!
- [Bug reporting](https://github.com/TeamLEGWORK/LEGWORK/issues/new) - If you see a bug we would love to know about it! Please create an issue [here](https://github.com/TeamLEGWORK/LEGWORK/issues/new)!
- [Release paper](https://arxiv.org/abs/2111.08717) - The LEGWORK release paper is now on the ArXiv and you can also view it directly [in GitHub](https://github.com/TeamLEGWORK/LEGWORK-paper) if you prefer!
# LEGWORK Changelog
This log keeps track of the changes implemented in each version of LEGWORK.

## 0.0.1
- Initial release

## 0.0.2
- Update version number to work with pip

## 0.0.3
*TW, 12/04/21*
- Allow the computation of merger times with ``get_merger_times``
    - After computing, times are used automatically
    in subsequent SNR calculations to avoid doubling up the computation
- Add ``evolve_sources`` function that evolves sources through time, updates merger times and if necessary marks them as merged
    - Merged sources are ignored in other computations (e.g. strain and SNR)
- Add ``ret_snr2_by_harmonic`` to eccentric snr functions to allow the user to get the SNR at each harmonic separately instead of the total
- Add minor fixes to snr for evolving sources for when sources are closes to merging

## 0.0.4
*TW, 19/05/21*
- Change visualisation module to be more flexible with **kwargs (allow any for dist plot and add linewidth to sensitivity curve function)
- Change Source.get_snr() to allow re-interpolation of the sensitivity curve for convenience (and fix the warning so it works properly)

## 0.0.5
*TW, 25/09/21*
- Avoid LSODA warnings by preventing integration from getting near the singularity at the merger
- Allow user to select how long before a merger to stop integration

## 0.0.6
*TW, 26/10/21*
- Avoid plotting merged sources in any of the automatic routines
- Allow source class evolution code to handle sources close to their merger
- Ensure SNR calculation works if some sources have merged and produces no warnings
- Change default behaviour of Source class with interpolate_g - no longer always interpolate, only when the collection of sources is fairly large or it contains eccentric sources
- Add a warning for if all timesteps are too close to the merger (based on `t_before`) and hence evolution can't happen

## 0.1.0
*TW, KB 31/10/21*
Major version change as we've added a significant enhancement with the new non-average SNR calculations.

- Change `snr` module to allow the calculation of non-averaged SNR using exact inclination, sky position and polarisations
- Let users specific inclination, sky position and polarisation in `Source` instantiation
- Add `VerificationBinaries` class to `Source` module for convenient access to LISA verification binary data from Kupfer+18
- Change max line length in code from 80 to 110 to increase readability

## 0.1.1
*TW, 5/11/21*
Small changes to visualisation code and updates to tutorials/demos with the new code
- allow user to specify weights in `Source` and visualisation functions
- allow user to customise sensitivity curve in any use of `plot_sources_on_sc`
- `legwork.__version__` now prints the version number
- Add new demo about verification binaries and other miscellaneous docs fixes

## 0.1.2
*TW 12/11/21*
- Change default values for `small_e_tol` and `large_e_tol` in `get_t_merge_ecc`
- allow users to specify custom confusion noise in sensitivity curves

## 0.1.3
*TW 16/11/21*
- Move `determine_stationarity` from `utils` to `evol` to avoid cyclical imports

## 0.1.4
*TW 17/11/21*
- Add `custom_psd` to the `sc_params` in `Source`
- Change the default parameters of `Source.get_snr()` to use `sc_params`

## 0.1.5
*KB 18/11/21*
- Make confusion noise follow shape of supplied frequency even in case of no confusion

## 0.1.6
*TW 18/11/21*
- Add a factor of 10/3 to the TianQin sensitivity curve to make it consistent with LISA (thanks to Yi-Ming Hu for pointing this out!)

## 0.1.7
*TW 05/01/21*
- Link GitHub releases to Zenodo

## 0.2.0
*TW 05/01/21*
- A couple of changes to how confusion noise is handled
    - End user can now access confusion noise functions directly through `get_confusion_noise`
    - Added confusion noise fit from Huang+20 to be used with TianQin
    - Added confusion noise fit from Thiele+21 which is based on a WDWD population with a metallicity dependent binary fraction
- TianQin psd function now include the confusion noise
- Change defaults used in `source` and `psd`
    - Often defaults for arm length, observation time and confusion noise were previously LISA related, LEGWORK now automatically works out the defaults based on what instrument is chosen
- Bug fix: in `visualisation` avoid mixing floats with Quantities when filling in a sensitivity curve

## 0.2.1
*TW 11/01/22*
- [Issue [#64](https://github.com/TeamLEGWORK/LEGWORK/issues/64)] Remove "auto" option from `interpolate_g` in favour of interpolating by default warning the user if they don't have many samples
- [Issue [#76](https://github.com/TeamLEGWORK/LEGWORK/issues/76)] Make it so that `strain.h_0_n` returns as unitless (same as `source.Source.get_h_0_n`). Same for `snr` functions.
- Fixed an issue introduced in 0.2.0 where automated observation times didn't work in `visualisation.plot_sources_on_sc_circ_stat`
- [Issue [#78](https://github.com/TeamLEGWORK/LEGWORK/issues/78)] Add a warning for when people are evolving past the merger with `avoid_merger=True`

*KB 13/01/22*
- [Issues [#79](https://github.com/TeamLEGWORK/LEGWORK/issues/79), [#80](https://github.com/TeamLEGWORK/LEGWORK/issues/80)] Add discussion of limitiations and scope of legwork snr calculations

## 0.2.2
*KB 18/01/22*
- [Issues [#68](https://github.com/TeamLEGWORK/LEGWORK/issues/68), [#84](https://github.com/TeamLEGWORK/LEGWORK/issues/84)] Removes upper bound on version limits for numpy and numba

## 0.2.3
*TW 23/01/22*
- [Issue [#75](https://github.com/TeamLEGWORK/LEGWORK/issues/75)] Fix mixing quantities with floats when plotting
- [Issue [#86](https://github.com/TeamLEGWORK/LEGWORK/issues/86)] Clarify how notebooks should be run and update the installation instructions

## 0.2.4
*TW 27/01/22*
- Make dependencies in setup.cfg match requirements.txt!

## 0.2.5
*TW 27/01/22*
- [Issue [#89](https://github.com/TeamLEGWORK/LEGWORK/issues/89)]
    - Created environment.yml for the package
    - Updated installation instructions to match
- [Issue [#90](https://github.com/TeamLEGWORK/LEGWORK/issues/90)] Added all psd functions to __all__

## 0.3.0
*TW 01/02/22*
- New major version of LEGWORK after several updates during the JOSS review (see 0.2.0-0.2.5)
---
name: Feature request
about: Got an idea for adding something to LEGWORK? Let us know!
title: ''
labels: enhancement
assignees: ''

---

**Describe the new feature you'd like**
A clear and concise description of what you want to happen. E.g. "I love how LEGWORK can compute SNRs, but I think a great improvement would be for it to also fold my laundry."

**Any suggestions for how to implement this**
Is there a particular way you would implement this? How could this integrate into the existing code?

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]
---
name: Bug report
about: Help us squash those bugs in LEGWORK!
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**Example for reproducing**
Please include a minimum working example that reproduces the error, e.g.
```
import legwork as lw
# do stuff with legwork
# something breaks
```

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. macOS]
 - Version: [e.g. v0.1.6]
 - Run from: [e.g. terminal or Jupyter notebook]
---
name: Documentation problem
about: A part of the documentation is wrong/misleading/unclear
title: ''
labels: documentation
assignees: ''

---

**Which documentation page is causing an issue?**
Please include the URL of the page. E.g. [The quickstart tutorial](https://legwork.readthedocs.io/en/latest/notebooks/Quickstart.html)

**What is the problem?**
Please quote the section and explain either what is misleading/unclear or point out the mistake

**Suggested solution (optional)**
Suggest how this problem would be solved
Scope and Limitations
=====================
LEGWORK is designed to provide quick estimates of the signal-to-noise ratio
of stellar-mass sources of mHz gravitational waves. These calculations can be
especially helpful when determining which sources in a large population containing
millions of binaries will be potentially detectable by LISA. If you are looking for
estimates of non-stellar-origin sources, check out the
`Gravitational Wave Universe Toolbox <http://gw-universe.org/>`_ or
`gwplotter <http://gwplotter.com/>`_. If you are looking for more advanced LISA simulator
tools, check out `lisacattools <https://github.com/tlittenberg/lisacattools>`_ or
`ldasoft <https://github.com/tlittenberg/ldasoft>`_.

The calculations done by LEGWORK apply the lowest-order post-Newtonian
description of gravitational wave emission and carefully follow the derivation of
`Flanagan and Hughes 1998a <https://ui.adsabs.harvard.edu/abs/1998PhRvD..57.4535F/abstract>`_.
This means that any higher order effects like spin-orbit
coupling or radiation reaction are not accounted for and thus if
a source is expected to depend on these higher order
effects the SNRs provided by LEGWORK may not be representative of the true SNR.
LEGWORK's SNRs are appropriate at mHz frequencies for sources with masses less than
a few tens of solar masses.

LEGWORK provides SNRs for two assumption cases: one where the detector orbit,
sky position of the source, and inclination of the source are all averaged and one where
the detector orbit is averaged while the sky position and inclination of of the source
are provided by the user following
`Cornish and Larson 2003 <https://ui.adsabs.harvard.edu/abs/2003PhRvD..67j3001C/abstract>`_.
The calculation which takes the sky position and inclination of the source into account
also accounts for frequency spreading due to doppler modulation from the detector's orbit.
This means that this calculation is only valid for circular sources and will *always* return
SNRs that are lower than the fully averaged SNR calculation because of the effects of the
frequency spreading.
 
LEGWORK modules and API reference
=================================

.. toctree::
   :maxdepth: 4

``LEGWORK`` is composed of 7 different modules, each with a different focus. The diagram below illustrates
how each module is connected to the others as well as listing the general purpose of each module. In particular,
note that the ``source`` module provides a simple interface to the functions in all other modules!

.. image:: https://github.com/TeamLEGWORK/LEGWORK-paper/raw/main/src/static/package_overview.png
   :width: 600
   :alt: Package structure graph
   :align: center

The rest of this page contains the API reference for each individual function in every module, feel free to
use the table of contents on the left to easily navigate to the function you need.

.. automodapi:: legwork.evol

.. tip::
    Feeling a bit spun around by all this binary evolution? Check out our
    tutorial on using functions in the ``evol`` module `here! <notebooks/Evolution.ipynb>`__

.. automodapi:: legwork.psd

.. automodapi:: legwork.snr

.. automodapi:: legwork.source

.. tip::
    Unable to find the source of your issues? Never fear! Check out our
    tutorial on using the Source class `here! <notebooks/Source.ipynb>`__

.. automodapi:: legwork.strain

.. tip::
    Feeling a little strained trying to parse these docs? Check out our
    tutorial on using functions in the ``strain`` module `here! <notebooks/Strains.ipynb>`__

.. automodapi:: legwork.utils

.. automodapi:: legwork.visualisation

.. tip::
    Not quite sure how things are working vis-à-vis visualisation? Check out our
    tutorial on using functions in the ``visualisation`` module `here! <notebooks/Visualisation.ipynb>`__Demos
=====

This section contains a series of demos that highlight ``LEGWORK``'s capabilities
with some example use cases. We'll continue to update this page with new demos
and if you've got an idea for a great use case for ``LEGWORK`` please `open a
pull request <https://www.github.com/TeamLEGWORK/LEGWORK/pulls>`_ and show off your new demo!

.. raw:: html

    <div class="toms-nav-container" style="margin-bottom:40px;">
        <div class="box" data-href="demos/BasicSNRCalculation.html">Calculate SNR of sources</div>
        <div class="box" data-href="demos/TheRoleofEccentricity.html">The Role of Eccentricity</div>
        <div class="box" data-href="demos/CompareSensitivityCurves.html">Compare Sensitivity Curves</div>
        <div class="box" data-href="demos/HorizonDistance.html">Horizon Distance</div>
    </div>

    <div class="toms-nav-container" style="margin-bottom:-50px;">
        <div class="box" data-href="demos/SNROverTime.html">Investigate SNR of source over time</div>
        <div class="box" data-href="demos/VerificationBinaries.html">LISA Verification Binaries</div>
    </div>

.. toctree::
    :titlesonly:
    :maxdepth: 1
    :hidden:

    Calculate SNR of sources <demos/BasicSNRCalculation.ipynb>
    The Role of Eccentricity <demos/TheRoleofEccentricity.ipynb>
    Compare Sensitivity Curves <demos/CompareSensitivityCurves.ipynb>
    Horizon Distance <demos/HorizonDistance.ipynb>
    Investigate SNR of source over time <demos/SNROverTime.ipynb>
    LISA Verification Binaries <demos/VerificationBinaries.ipynb>
.. image:: images/legwork.png
   :width: 400
   :alt: logo for legwork
   :align: center

|

**LEGWORK** (**L**\ ISA **E**\ volution and **G**\ ravitational **W**\ ave **OR**\ bit **K**\ it) is a python package designed to calculate signal-to-noise ratios for GWs emitted from inspiraling binary systems that are potentially observable by LISA.

LEGWORK also contains a plotting module so you can show off those detectable GW sources of yours and treats any and all kinds of stellar-origin inspiralling binaries. 
Circular and stationary? No problem!
Eccentric and chirping? We've got you covered!

Want to determine if *your* favorite source is detectable by LISA? Let LEGWORK do the legwork and keep the LISA literature free from spurious factors of 2!

.. raw:: html

    <div class="toms-nav-container" style="margin-bottom:50px;">
        <div class="box" data-href="install.html">Install LEGWORK</div>
        <div class="box" data-href="tutorials.html">Learn from tutorials</div>
        <div class="box" data-href="modules.html">Explore the modules</div>
        <div class="box" data-href="notebooks/Derivations.html">Delve into the derivations</div>
    </div>

.. toctree::
    :maxdepth: 1
    :hidden:

    Installation <install>
    Citing LEGWORK <cite>
    Tutorials <tutorials>
    Demos <demos>
    Modules <modules>
    Derivations <notebooks/Derivations.ipynb>
    Limitations and Scope <limitations>
    GitHub <https://github.com/TeamLEGWORK/LEGWORK>
    Submit an issue <https://github.com/TeamLEGWORK/LEGWORK/issues/new>
Citing LEGWORK
==============

If you use LEGWORK in a scientific publication we ask that you please cite it. Below we include a ready-made
BibTeX entry:

.. code-block:: bash

    @ARTICLE{Wagg+2021,
        author = {{Wagg}, Tom and {Breivik}, Katelyn and {de Mink}, Selma E.},
        title = "{LEGWORK: A python package for computing the evolution and detectability of stellar-origin gravitational-wave sources with space-based detectors}",
        journal = {arXiv e-prints},
        keywords = {Astrophysics - High Energy Astrophysical Phenomena, General Relativity and Quantum Cosmology},
        year = 2021,
        month = nov,
        eid = {arXiv:2111.08717},
        pages = {arXiv:2111.08717},
        archivePrefix = {arXiv},
        eprint = {2111.08717},
        primaryClass = {astro-ph.HE},
        adsurl = {https://ui.adsabs.harvard.edu/abs/2021arXiv211108717W},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

More citation methods are `available on ADS <https://ui.adsabs.harvard.edu/abs/2021arXiv211108717W/exportcitation>`__ if you'd rather use a different format!Tutorials
=========

This section contains a series of tutorials on how to use and get the most out
of ``LEGWORK``. We are working to add more tutorials over time and are very open
to feedback!

The fastest way to get started is our quickstart tutorial:

.. raw:: html

    <div style="max-width: 400px; margin: 30px auto;">
        <div class="toms-nav-box" data-href="notebooks/Quickstart.html">Quickstart</div>
    </div>

For more in-depth details and examples, check out our other tutorials:

.. raw:: html

    <div class="toms-nav-container" style="margin-bottom:30px;">
        <div class="box" data-href="notebooks/Source.html">Using the Source Class</div>
        <div class="box" data-href="notebooks/Strains.html">Computing Gravitational Wave Strains</div>
        <div class="box" data-href="notebooks/Evolution.html">Binary Evolution</div>
        <div class="box" data-href="notebooks/Visualisation.html">Visualisation</div>
    </div>

    <div style="max-width: 300px; margin: 40px auto 50px auto;">
        <div class="toms-nav-box" data-href="notebooks/Units.html">Units in LEGWORK</div>
    </div>

.. toctree::
    :titlesonly:
    :maxdepth: 1
    :hidden:

    Quickstart <notebooks/Quickstart.ipynb>
    Using the Source Class <notebooks/Source.ipynb>
    Computing Gravitational Wave Strains <notebooks/Strains.ipynb>
    Binary Evolution <notebooks/Evolution.ipynb>
    Visualisation <notebooks/Visualisation.ipynb>
    Understanding units in LEGWORK <notebooks/Units.ipynb>
Installation
============

.. tabs::

    .. tab:: Stable (with conda)

        This is our recommend installation method! Follow the steps below to start using ``LEGWORK``!

        #. :download:`Download the environment.yml file from our repository <https://raw.githubusercontent.com/TeamLEGWORK/LEGWORK/main/environment.yml>`
        #. Create a new conda environment using this file

            .. code-block:: bash

                conda env create -f path/to/environment.yml

        #. Activate the environment by running

            .. code-block:: bash

                conda activate legwork

        and you should be all set! Check out our `quickstart tutorial <notebooks/Quickstart.ipynb>`__ to learn some LEGWORK basics.
        Note that if you also want to work with the notebooks in the tutorials and/or demos you'll also need to install jupyter/ipython in this environment!

    .. tab:: Stable (without conda)

        We don't recommend installing ``LEGWORK`` without a conda environment but if you prefer to do it this
        way then all you need to do is run

        .. code-block:: bash

            pip install legwork

        and you should be all set! Check out our `quickstart tutorial <notebooks/Quickstart.ipynb>`__ to learn some LEGWORK basics.
        Note that if you also want to work with the notebooks in the tutorials and/or demos you'll also need to install jupyter/ipython in this environment!

    .. tab:: Development (from GitHub)
        
        .. warning::

            We don't guarantee that there won't be mistakes or bugs in the development version, use at your own risk!

        The latest development version is available directly from our `GitHub Repo
        <https://github.com/TeamLEGWORK/LEGWORK>`_. To start, clone the repository onto your machine:

        .. code-block:: bash
        
            git clone https://github.com/TeamLEGWORK/LEGWORK
            cd LEGWORK

        Next, we recommend that you create a Conda environment for working with LEGWORK.
        You can do this by running

        .. code-block:: bash

            conda create --name legwork "python>=3.7" pip "numba>=0.50" "numpy>=1.16" "astropy>=4.0" "scipy>=1.5.0" "matplotlib>=3.3.2" "seaborn>=0.11.1" "schwimmbad>=0.3.2" -c conda-forge -c defaults

        And then activate the environment by running

        .. code-block:: bash

            conda activate legwork

        At this point, all that's left to do is install LEGWORK!

        .. code-block:: bash

            pip install .

        and you should be all set! Check out our `quickstart tutorial <notebooks/Quickstart.ipynb>`__ to learn some LEGWORK basics.
        Note that if you also want to work with the notebooks in the tutorials and/or demos you'll also need to install jupyter/ipython in this environment!