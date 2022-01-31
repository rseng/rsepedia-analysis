# Changelog

Changes and updates to EMD are tracked by version on this page.  The format of
this changelog is (mostly) based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this the EMD package uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Changes should be categorised under the following types:

- **Added** for new features.
- **Changed** for changes in existing functionality.
- **Deprecated** for soon-to-be removed features.
- **Removed** for now removed features.
- **Fixed** for any bug fixes.
- **Security** in case of vulnerabilities.

Where appropriate, links to specific Issues & Merge Requests on [our gitlab page](https://gitlab.com/emd-dev/emd).


## Development Version
Work in progress...

    git clone https://gitlab.com/emd-dev/emd.git

---

# Stable Versions

## 0.5.2

    pip install emd==0.5.2
Released 2022-01-13

### Added

- New tutorial on harmonic structures and instantaneous frequency [!74](https://gitlab.com/emd-dev/emd/-/merge_requests/74)

---

## 0.5.1

    pip install emd==0.5.1
Released 2021-12-17

### Notes

This release fixes a bug with the initial v0.5.0 release.

---

## 0.5.0

    pip install emd==0.5.0
Released 2021-12-17

### Notes
This release contains several breaking changes in the emd.spectra submodule
particularly in the emd.spectra.hilberthuang and emd.spectra.holospectrum
functions. Please see the relevant docstrings for a review of new API.

### Added
- Implementation of the Iterated Mask Sift from @marcofabus [!69](https://gitlab.com/emd-dev/emd/-/merge_requests/69)
- Support for Python 3.9 - requirements updated and test build added [!62](https://gitlab.com/emd-dev/emd/-/merge_requests/62)
- Dependancy on Sparse for spectrum computation (https://sparse.pydata.org/en/stable/) [!70](https://gitlab.com/emd-dev/emd/-/merge_requests/70)
- New citations page on website (and assorted website/reference fixups) [!60](https://gitlab.com/emd-dev/emd/-/merge_requests/60)

### Changed
- BREAKING: emd.spectra.hilberthuang and emd.spectra.holospectrum API changed, also now return frequency bin vectors. [!70](https://gitlab.com/emd-dev/emd/-/merge_requests/70)
- Major refector of spectrum code [!70](https://gitlab.com/emd-dev/emd/-/merge_requests/70)
- Major refactor of extrema padding underlying the sift [!71](https://gitlab.com/emd-dev/emd/-/merge_requests/71)

### Fixed
- Average over only the ensembles with the modal number of IMFs in ensemble_sift, previously there could be an error due to difference in nimfs between ensembles [!63](https://gitlab.com/emd-dev/emd/-/merge_requests/63)
- Fix pyyaml dependency for google-colab import [!64](https://gitlab.com/emd-dev/emd/-/merge_requests/64)
- Fix case where mask sift could drop some options [!66](https://gitlab.com/emd-dev/emd/-/merge_requests/66)
- Fix bug in energy ratio sift stopping method [95295a8c](https://gitlab.com/emd-dev/emd/-/merge_requests/71/diffs?commit_id=95295a8ca992a0df13597c39924af61eca7130bc)
- Fix bug in zero-crossing count mask frequency selection method [fe2605c2](https://gitlab.com/emd-dev/emd/-/merge_requests/71/diffs?commit_id=fe2605c2ca0730750912a1d5af1d2c52a27b142a)

### Removed
- emd.spectra.frequency_stats removed - replaced by emd.spectra.frequency_transform
- emd.spectra.hilberthuang_1d merged into emd.spectra.hilberthuang

## 0.4.0

    pip install emd==0.4.0
Released 2021-03-30

### Notes
Many changes in this revision come from the review process at [JOSS](https://github.com/openjournals/joss-reviews/issues/2977)

### Added
- New tutorials
  - Cross-frequency coupling [!52](https://gitlab.com/emd-dev/emd/-/merge_requests/52)
  - Why use EMD? [90a5b5e2](https://gitlab.com/emd-dev/emd/-/commit/90a5b5e2e4ffdd7634cf63e30836843c920fcaa3)
  - Tutorial on code speed [35ae8c82](https://gitlab.com/emd-dev/emd/-/commit/35ae8c82ab72b9c36d641eb4bef4d1ae7c53b0a5)
- Second layer mask sift function [65a05dd2](https://gitlab.com/emd-dev/emd/-/commit/65a05dd2cf1610508d13a68f6753094f07d67e48)
- Add html printing functionality for SiftConfig [5c57781e](https://gitlab.com/emd-dev/emd/-/commit/5c57781e2e8b92b8d2a7e00ceec8bde064bc412b)
- Update contribution and installation details on website - add accordions for better readability [!53](https://gitlab.com/emd-dev/emd/-/merge_requests/53)
- Add new plotting functionality for HHT and Holospectra [!53](https://gitlab.com/emd-dev/emd/-/merge_requests/53)
- Show warning when max_imfs is very high compared to length of time-series [4cd15291](https://gitlab.com/emd-dev/emd/-/commit/4cd15291c25e082cbb9ffb56a2c3812b6b3d391e)

### Changed
- Major refactor in handling of cycles analysis [!56](https://gitlab.com/emd-dev/emd/-/merge_requests/56)
  - Introduce Cycles class
  - Introduce \_cycles\_support module
- Renamed 'References' webpage to 'API' [d8fe93b5](https://gitlab.com/emd-dev/emd/-/commit/d8fe93b520c19ce45f3f5a73294074a4b1d75ce5)

### Fixed
- Widespread fixing of typos and mistakes in documentation & website [!52](https://gitlab.com/emd-dev/emd/-/merge_requests/52)
- Make docstrings pydocstyle compliant and add pydocstyle conventions [!53](https://gitlab.com/emd-dev/emd/-/merge_requests/53)
- Large number of pylint recommended fixes [271d7937](https://gitlab.com/emd-dev/emd/-/commit/271d793731fad64902f16323493ee06893002286)
- Indexing typo fixed in bin_by_phase [c5679432](https://gitlab.com/emd-dev/emd/-/commit/c5679432cfcd011965547144aaa936eee1405f62)
- Improve label alignments in plot_imfs [!54](https://gitlab.com/emd-dev/emd/-/merge_requests/56)

---

## 0.3.3

    pip install emd==0.3.3
Released 2021-02-04

### Added
- New function for computing summary stats from chains of cycles (from marcoFabus) [!46](https://gitlab.com/emd-dev/emd/-/merge_requests/46)

### Changed
- Major updates to tutorials [!40](https://gitlab.com/emd-dev/emd/-/merge_requests/40)
  - Binder notebooks added
  - New sifting tutorials added

### Fixed
- Replaced missing dependencies in setup.py [!42](https://gitlab.com/emd-dev/emd/-/merge_requests/42)

---

## 0.3.2

    pip install emd==0.3.2
Released 2020-11-29

### Added
- Add input array shape ensurance functions and start to use in sift & cycle submodules  [!26](https://gitlab.com/emd-dev/emd/-/merge_requests/26)
- Add more stopping criteria to sift module [!27](https://gitlab.com/emd-dev/emd/-/merge_requests/26)
  - Rilling et al and fixed iterations IMF stopping criteria
  - Energy threshold sift stopping criterion


### Changed
- Refactor some options extrema detection functions [!29](https://gitlab.com/emd-dev/emd/-/merge_requests/29)
- Sift throws an error if an IMF doesn't converge after a specified maximum number of iterations.
- Refactor mask generation in mask sift. Now specifies N masks of different phases and has options for parallel processing.
- SiftConfig yaml file also stores which type of sift the config is for [!35](https://gitlab.com/emd-dev/emd/-/merge_requests/35)
- 18% increase in testing coverage (to 75% total) [!30](https://gitlab.com/emd-dev/emd/-/merge_requests/30)

### Deprecated
- emd.spectra.frequency_stats renamed to emd.spectra.frequency_transform. Original func kept for now.

---

## 0.3.1

    pip install emd==0.3.1
Released 2020-09-06

### Added
- This changelog [!18](https://gitlab.com/emd-dev/emd/-/merge_requests/18)
- support.py submodule with some helper functions for checking installs and running tests [!20](https://gitlab.com/emd-dev/emd/-/merge_requests/20)
- envs subdir containing anaconda install environment config files [!21](https://gitlab.com/emd-dev/emd/-/merge_requests/21)
- Options for reading and writing sift configs to yaml format [!24](https://gitlab.com/emd-dev/emd/-/merge_requests/24)
- major update to webpage [!12](https://gitlab.com/emd-dev/emd/-/merge_requests/24)
  - Reformat page to bootstrap
  - Add & update the tutorials
  - New landing page

### Fixed
- Input array dimensions in phase_align clarified and fixed up [ef28b36c](https://gitlab.com/emd-dev/emd/-/commit/ef28b36cac8be7224280fd7ba02d25b3f084ab30)
- Extrema opts were dropped in get_next_imf [!23](https://gitlab.com/emd-dev/emd/-/merge_requests/23)

### Changed
- get_control_points internal refector [af153ed6](https://gitlab.com/emd-dev/emd/-/commit/af153ed606601f3963c125329c86710e47c06b45)

---

## 0.3.0

    pip install emd==0.3.0
Released on 2020-07-22

### Added
- get_cycle_stat refectored to allow general numpy and user-specified metrics to be computed
- Logger coverage increased, particularly in cycle.py
  - Logger exit message added

### Changed
- Major SiftConfig refactor - API & syntax now much cleaner

---

## 0.2.0

    pip install emd==0.2.0
Released 2020-06-05

### Added
- Tutorials on the sift, hilbert-huang and holospectrum analyses.
- Parabolic extrema interpolation
- Average envelope scaling in sift
- Testing expanded to include python 3.5, 3.6, 3.7 & 3.8


### Changed
- API in sift functions updated for compatabillity with new SiftConfig
  - Expose options for extrema padding to top level sift function
  - Sift relevant util functions moved into sift.py submodule
- Masked sift functions merged into single function
- get_cycle_chain refactor to cleaner internal syntax

---

## 0.1.0

    pip install emd==0.1.0
Released 2019-12-10

### Added
- Everything
A python package for Empirical Mode Decomposition and related spectral analyses.

Please note that this project is in active development for the moment - the API may change relatively quickly between releases!

# Installation

You can install the latest stable release from the PyPI repository

```
pip install emd
```

or clone and install the source code.

```
git clone https://gitlab.com/emd-dev/emd.git
cd emd
pip install .
```

Requirements are specified in requirements.txt. Main functionality only depends
on numpy and scipy for computation and matplotlib for visualisation.

# Quick Start

Full documentation can be found at https://emd.readthedocs.org and development/issue tracking at gitlab.com/emd-dev/emd

Import emd

```python
import emd
```

Define a simulated waveform containing a non-linear wave at 5Hz and a sinusoid at 1Hz.

```python
import emd
sample_rate = 1000
seconds = 10
num_samples = sample_rate*seconds

import numpy as np
time_vect = np.linspace(0, seconds, num_samples)

freq = 5
nonlinearity_deg = .25  # change extent of deformation from sinusoidal shape [-1 to 1]
nonlinearity_phi = -np.pi/4  # change left-right skew of deformation [-pi to pi]
x = emd.utils.abreu2010(freq, nonlinearity_deg, nonlinearity_phi, sample_rate, seconds)
x += np.cos(2*np.pi*1*time_vect)
```

Estimate IMFs

```python
imf = emd.sift.sift(x)
```

Compute instantaneous frequency, phase and amplitude using the Normalised Hilbert Transform Method.

```python
IP, IF, IA = emd.spectra.frequency_transform(imf, sample_rate, 'hilbert')
```
Compute Hilbert-Huang spectrum

```python
freq_edges, freq_bins = emd.spectra.define_hist_bins(0, 10, 100)
hht = emd.spectra.hilberthuang(IF, IA, freq_edges)
```
```
Make a summary plot

```python
import matplotlib.pyplot as plt
plt.figure(figsize=(16, 8))
plt.subplot(211, frameon=False)
plt.plot(time_vect, x, 'k')
plt.plot(time_vect, imf[:, 0]-4, 'r')
plt.plot(time_vect, imf[:, 1]-8, 'g')
plt.plot(time_vect, imf[:, 2]-12, 'b')
plt.xlim(time_vect[0], time_vect[-1])
plt.grid(True)
plt.subplot(212)
plt.pcolormesh(time_vect, freq_bins, hht, cmap='ocean_r')
plt.ylabel('Frequency (Hz)')
plt.xlabel('Time (secs)')
plt.grid(True)
plt.show()
```

# EMD Install Environments

This folder contains some environment config files to help install EMD. These
files can be copied from or directly downloaded and used to install EMD.

### Anaconda envs
- emd_conda_env.yml: Installs the latest EMD version from PyPI (stable version)
- emd-dev_conda_env.yml: Installs the latest EMD version from gitlab (unstable-development version)
- emd_tutorial_conda_env.yml: Installs EMD from PyPI (stable) alongside ipython, spyder and jupyter.
---
title: "EMD: Empirical Mode Decomposition and Hilbert-Huang Spectral Analyses in Python"
tags:
  - Python
  - Time-series
  - Non-linear
  - Dynamics
authors:
  - name: Andrew J. Quinn
    affiliation: 1
  - name: Vitor Lopes-dos-Santos
    affiliation: 2
  - name: David Dupret
    affiliation: 2
  - name: Anna Christina Nobre
    affiliation: "1,3"
  - name: Mark W. Woolrich
    affiliation: 1
affiliations:
  - name: Oxford Centre for Human Brain Activity, Wellcome Centre for Integrative Neuroimaging, Department of Psychiatry, University of Oxford, Oxford, UK
    index: 1
  - name: Medical Research Council Brain Network Dynamics Unit, Nuffield Department of Clinical Neurosciences, University of Oxford, Oxford, OX1 3TH, United Kingdom
    index: 2
  - name: Department of Experimental Psychology, University of Oxford, Oxford, OX2 6GG, UK
    index: 3

date: 30 November 2020
bibliography: paper.bib
---

# Summary

The Empirical Mode Decomposition ([`EMD`](https://emd.readthedocs.io/en/latest/))
package contains Python (>=3.5) functions for analysis of non-linear and
non-stationary oscillatory time series. `EMD` implements a family of sifting
algorithms, instantaneous frequency transformations, power spectrum
construction and single-cycle feature analysis. These implementations are
supported by online documentation containing a range of practical tutorials.

# Statement of Need

Many oscillatory signals contain non-linear or non-sinusoidal features that
change dynamically over time. These complex and dynamic features are often of
analytic interest but can be challenging to isolate and quantify. The Empirical
Mode Decomposition offers a potential solution defined by the *sift algorithm*, a
data-adaptive decomposition that separates a signal into a set of Intrinsic
Mode Functions (IMFs) that permit physically interpretable Hilbert transforms
[@Huang1998] and subsequent analysis of instantaneous frequency. Crucially, the
sift is able to efficiently isolate and describe non-linear and non-stationary
signal features as it works on adaptive, local data segments without
prescribing that features remain consistent across the entire signal.

# Package Features

Empirical Mode Decomposition is defined by the 'sift' algorithm (@Huang1998).
This is a time-domain process which looks to isolate the fastest dynamics in a
time-series by iteratively sifting out slower dynamics.  Any slow dynamics are
removed by subtracting the average of the signal's upper and lower amplitude
envelope until that average is sufficiently close to zero. This isolated signal
component is known as an Intrinsic Mode Function (IMF); it is subtracted from
the original signal and the sifting process repeated to identify the next IMF,
which will contain slower dynamics. This process is repeated until only a trend
remains in the signal.

The sift algorithm is implemented in the `emd.sift` module, including the
classic sift (`emd.sift.sift`; @Huang1998), the Ensemble EMD
(`emd.sift.ensemble_sift`; @Wu2009), Masked EMD (`emd.sift.mask_sift`;
@Deering2005) and the second-level sift (`emd.sift.sift_second_layer`;
@Huang2016). The ensemble and masked sift variants can be optionally
accelerated by parallel processing (though this is not possible in all variants
of the sift algorithm). An example set of Intrinsic Mode Functions isolated by
a Masked-Sift is shown in Figure 1. The sift functions rest upon a range of
lower-level utility functions, which can be customised and used directly if
needed. All levels of the sift computation are customisable from the top-level
sift functions. Users can configure these sift options using a dictionary-like
`emd.sift.SiftConfig` object. This config can then be passed directly to the
sift functions or saved in `YAML` format for later use or sharing.

Each IMF can be analysed in terms of its instantaneous frequency
characteristics at the full temporal resolution of the dataset [@Huang2009].
The Hilbert transform is used to construct an energy-frequency or
energy-frequency-time spectrum known as the Hilbert-Huang Transform (HHT). A
second level decomposition of the amplitude modulations of each IMF extends the
HHT to the Holospectrum, describing signal energy across carrier frequency,
amplitude modulation frequency and time [@Huang2016]. The frequency transforms are
implemented in the `emd.spectra` submodule. `emd.spectra.frequency_stats`
implements a set of methods for computing instantaneous frequency, phase and
amplitude from a set of IMFs. These can be used as inputs to the
`emd.spectra.hilberthuang` or `emd.spectra.holospectrum` to obtain energy
distributions across time and frequency (see examples in Figures 3 and 4). The
Hilbert-Huang and Holospectrum computations can be very large, so these
functions use an efficient sparse array implementation.

The `EMD` toolbox provides a range of functions for the detection of oscillatory
cycles from the IMFs of a signal. Once identified, each cycle can be
characterised by a range of features, including its amplitude, frequency and
waveform shape. Tools are provided for detecting continuous chains of
oscillatory cycles and for matching similar cycles across datasets. The cycle
analysis functions are implemented in `emd.cycle`.

A range of utility and support features are included in the `EMD` toolbox.
Firstly, a customisable logger (implemented in `emd.logger`) is threaded
throughout the toolbox to provide progress output about ongoing computations,
warnings and errors. The logger output may be augmented by the user and any
output can be directed to a specified log file in addition to the console.
Secondly, `EMD` is supported by a range of tests, implemented in the `py.test`
framework. These include both routine usage tests and tests ensuring that the
behaviour of the sift routines meet a set of pre-specified requirements.
Finally, `emd.support` contains a set of functions for running tests and
checking which versions of `EMD` are currently installed and whether updates
are available on [PyPI](https://pypi.org/project/emd/).

# Target Audience

Since its initial publication in 1998, the EMD approach has had a wide impact
across science and engineering, finding applications in turbulence, fluid
dynamics, geology, biophysics and neuroscience amongst many others. The `EMD`
toolbox will be of interest to scientists, engineers and applied mathematicians
looking to characterise signals with rich dynamics with a high temporal and
spectral resolution.

# State of the field

The popularity of the EMD algorithm has led to several
implementations which offer overlapping functionality. Here we include an
incomplete list of these toolboxes providing sift, ensemble sift and HHT
implementations. In Python there are two substantial EMD implementations
available on the PyPI server: [PyEMD](https://github.com/laszukdawid/PyEMD) [@pyemd]
and [PyHHT](https://pyhht.readthedocs.io/en/latest/) [@pyHHT]. Each of these packages
implements a family of sifting routines and frequency transforms. Another
implementation of EMD, in Matlab and C, is available from [Patrick
Flandrin](http://perso.ens-lyon.fr/patrick.flandrin/emd.html) [@flandrin]. This provides a
wide range of sift functions, but limited frequency transform or spectrum
computations. Finally, the basic EMD algorithm and HHT is implemented in the
[MatLab signal processing
toolbox](https://uk.mathworks.com/help/signal/ref/emd.html) [@matlabsignal].

The `EMD` toolbox covers much of the functionality in these packages within a
single computational framework. Beyond these methods, we add fully-featured
implementations of masked sift and second-level sift routines, as well as the
first Python implementation of higher-level Holospectrum analyses. Finally, we
offer a suite of tools designed for analysis of single-cycles of an Intrinsic
Mode Function.

# Installation & Contribution

The `EMD` package is implemented in Python (>=3.5) and is freely available
under a GPL-3 license. The stable version of the package can be installed from
from PyPI.org using ```pip install emd```. Users and developers can also
install from source from [gitlab](https://gitlab.com/emd-dev/emd). Our
[documentation](https://emd.readthedocs.io) provides detailed instructions on
[installation](https://emd.readthedocs.io/en/latest/install.html) and a range
of practical
[tutorials](https://emd.readthedocs.io/en/latest/emd_tutorials/index.html).
Finally, users wishing to submit bug reports or merge-requests are able to do
so on our gitlab page following our [contribution
guidelines](https://emd.readthedocs.io/en/latest/contribute.html).

![A simulated signal with an oscillatory component (black line - top panel) with a set of intrinsic mode functions estimated using a mask sift EMD (coloured lines - lower panels).](figures/emd_joss_example1_sift.png)

![A segment of a simulated signal with its instantaneous amplitude and a time-series containing the maxiumum amplitude of each successive cycle..](figures/emd_joss_example2_amp.png)

![Top panel: An Instrinsic Mode function from a simulated signal (black line) and an amplitude threshold (dotted line). Bottom Panel: 2D Hilbert-Huang Transform. Darker colours indicate greater power and the black lines indicate cycle average instantaneous frequency of large amplitude cycles.](figures/emd_joss_example3_hht.png)

![Top panel: A segment of a simulated containing two nested oscillations and white noise. One 5Hz oscillation with 0.5Hz amplitude modulation and a 37Hz signal whose amplitude is modulated by the lower-frequency 5Hz oscillation. Bottom left: The 1D Hilbert-Huang transform of this signal. Bottom Center: The 2D Hilbert-Huang transform. Bottom Right: The Holospectrum.](figures/emd_joss_example4_holo.png)

\pagebreak

# Acknowledgements

We thank Norden Huang, Chi-Hung Juan, Jia-Rong Yeh and Wei-Kuang
Liang for enjoyable and fruitful discussions on EMD theory and applications in
recent years. We also thank Jasper Hajonides van der Meulen and
Irene Echeverria-Altuna for their time, patience and feedback on early versions
of this toolbox.

This project was supported by the Medical Research Council (RG94383/RG89702)
and by the NIHR Oxford Health Biomedical Research Centre. The Wellcome Centre
for Integrative Neuroimaging is supported by core funding from the Wellcome
Trust (203139/Z/16/Z). V.L.d.S. and D.D. are supported by the Medical Research
Council UK (Programmes MC_UU_12024/3 and MC_UU_00003/4 to D.D.) ACN is
supported by the Wellcome Trust (104571/Z/14/Z) and James S. McDonnell
foundation (220020448). MWW is supported by the Wellcome Trust (106183/Z/14/Z;
215573/Z/19/Z). ACN and MWW are further supported by an EU European Training
Network grant (euSSN; 860563).

# References
Contribute to EMD
=================

Thank you for your interest in contributing to EMD! Development of EMD takes place on `our gitlab page <https://gitlab.com/emd-dev/emd>`_. You can contribute to the developement of the EMD toolbox through gitlab by identifying and raising issues or by submitting new code through merge requests. Note that these will require an active account on `gitlab.com <https://www.gitlab.com>`_.

Both issues and merge requests are very welcome! We will try to resolve all issues but please bear in mind that this is a small open-source project with limited developement time.

Issues
------

Our `issue tracker <https://gitlab.com/emd-dev/emd/-/issues>`_ is the place to submit tickets about things that could or should change in the EMD toolbox. These could be bugs or about any problems you find, or requests for new functionality.

When submitting information about a bug, please include the following information so that the issue can be easily understood and reproduced.

- Expected Behaviour, what should happen when the code runs?
- Actual Behaviour, what is actually happening?
- Steps to Reproduce the Problem, what would another person need to do to see the issue?
- System Specifications, what operating system are you using? which versions of python and emd are you running?

Once the ticket has been submitted, we can take a look at the issue and may use the comments section at the bottom of the issue page to ask for more information and discuss a solution. Once a solution has been found, it will be submitted in a merge request linked to the issue in question.

Merge Requests
--------------

It it not possible to commit directly to the master branch of EMD, so all updates must first appear as `merge requests <https://gitlab.com/emd-dev/emd/-/merge_requests>`_. Each merge requests corresponds to a git branch containing a set of changes that may be added into EMD.

The merge request page provides a lot of useful information about the branch.

- Is the branch up to date with master? if not a rebase or merge may be required.
- Are the tests passing? We cannot merge code which breaks our tests! We test for a range of features including code behaviour, flake8 style compliance and spelling (using ``codespell``, this is particularly important for docstrings and tutorials)

At the start, all merge requests are marked as Work In Progress (WIP), this means that the code is not ready for merge yet. This tests can be run and process tracked as commits are added to the branch. Please feel free to use the comments section to discuss any changes and ask for feedback.

See below for detailed descriptions of the contribution process in a couple of different cases.

Create a merge-request....
**************************

.. container:: toggle body

    .. container:: header body

        .. raw:: html

            <h3 class='installbar'>... from an issue ticket</h3>

    .. container:: installbody body

        This section outlines how to contribute to by the issue-tracker on gitlab.com/emd-dev/emd. This is the best method to use if the changes will be made with contributions from several people.

        1. First, create an issue in the EMD `issue tracker <https://gitlab.com/emd-dev/emd/-/issues>`_. The issue should clearly introduce the potential changes that will be made.
        2. The issue will be read by an emd-dev developer who may ask for more information or suggest another solution in the issue discussion.
        3. If everyone agrees that a change is needed -  the developer will create a new branch and merge request which are specifically linked to the issue.
        4. The branch will be publicly accessible, you can install EMD with this branch using the methods in the `developer section of the install page <file:///Users/andrew/src/emd/doc/build/html/install.html#development-gitlab-com-version>`_.
        5. You and the developer can the work on the branch until the changes are finalised. Feel free to use the discussion section of the merge request.
        6. Once the work is complete, run some final checks on your local branch.
            - Ensure the tests are passing
            - Ensure that the code to be committed is flake8 compliant
            - Briefly describe the changes in the appropriate section of the changelog.
        7. The developer can then approve the merge request and the changes will be accepted into the main EMD branch.


.. container:: toggle body

    .. container:: header body

        .. raw:: html

            <h3 class='installbar'>... from a fork of EMD</h3>

    .. container:: installbody body

        This section outlines how to contribute to EMD from your own fork of the repository. This might be the simplest method if you would like to configure the gitlab.com environment and/or keep any changes private during development.

        1. First, create a fork of EMD from the `gitlab repository <https://gitlab.com/emd-dev/emd>`_.
        2. Create a branch for your changes in the fork of EMD. Any contributions must come from a branch - don't merge the branch into master in the forked repository.
        3. Ensure that runners are enabled so that tests can run on gitlab.com. "Settings -> CI/CD -> Public Pipelines" should be ticked in the gitlab.com settings.
        4. Complete your work on the branch.
        5. Run some final checks on your local branch.
            - Ensure that your branch is up to date with the main branch on emd-dev - this may require updating your fork.
            - Ensure the tests are passing
            - Ensure that the code to be committed is flake8 compliant
            - Briefly describe the changes in the appropriate section of the changelog.
        6. Submit the merge request from your fork of EMD. On your fork of gitlab.com go to "Repository -> Branches" and click 'Merge request" next to corresponding branch.
        7. The request will intially be marked as a Work In Progress (WIP). We will review the changes and potentially request some final changes or tweaks in the discussion on the Merge Request page.
        8. Once the developers are happy that the changes are ready, WIP status will be updated and the branch merged into the main EMD branch.


Quick Start
===========

EMD can be install from `PyPI <https://pypi.org/project/emd/>`_ using pip::

    pip install emd

and used to decompose and describe non-linear timeseries.::

    # Imports
    import emd
    import numpy as np
    import matplotlib.pyplot as plt

    # Definitions
    sample_rate = 1000
    seconds = 3
    time_vect = np.linspace(0,seconds,seconds*sample_rate)

    # A non-linear oscillation
    x = emd.utils.abreu2010( 5, .25, -np.pi/4, sample_rate, seconds )
    # ...plus a linear oscillation
    x += np.cos( 2*np.pi*1*time_vect )

    # Sift
    imf = emd.sift.sift( x )

    # Visualise Intrinsic Mode Functions
    emd.plotting.plot_imfs( imf, scale_y=True, cmap=True )

    # Compute instantaneous spectral stats
    IP,IF,IA = emd.spectra.frequency_transform( imf, sample_rate ,'nht' )

    # Compute Hilbert-Huang transform
    edges,centres = emd.spectra.define_hist_bins(0,10,32)
    hht = emd.spectra.hilberthuang( IF, IA, edges )

    # Visualise time-frequency spectrum
    plt.figure()
    plt.pcolormesh( time_vect, centres, hht, cmap='hot_r')
    plt.colorbar()
    plt.xlabel('Time (seconds)')
    plt.ylabel('Instantaneous  Frequency (Hz)')
.. emd documentation master file, created by
   sphinx-quickstart on Sun Jan 27 23:11:40 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

EMD: Empirical Mode Decomposition
=================================

Python tools for the extraction and analysis of non-linear and non-stationary oscillatory signals.

.. title image, description
.. raw:: html

    <div class="nopad">
      <img src="_static/emd_example.png" style="max-width: 768px;  alt="EMD">
      </div>

Features
========

* Sift algorithms including the ensemble sift, complete ensemble sift and mask sift
* Instantaneous phase, frequency and amplitude computation
* Cycle detection and analysis
* Hilbert-Huang spectrum estimation (1d frequency spectrum or 2d time-frequency spectrum)
* Second layer sift to quantify structure in amplitude modulations
* Holospectrum estimation (3d instantaneous frequency x amplitude modulation frequency x time spectrum)

.. toctree::
   :maxdepth: 2

   Install<install>
   Tutorials<emd_tutorials/index>
   API Reference<api>
   Contribute<contribute>
   Changelog<changelog>
   Cite<cite>
Citing the EMD Package
=================================
|

If this package has significantly contributed to your work, please include the following citation:

.. title image, description
.. raw:: html

    {Quinn2021_joss}

.. highlight:: none
::

    @article{Quinn2021_joss,
      doi = {10.21105/joss.02977},
      url = {https://doi.org/10.21105/joss.02977},
      year = {2021},
      publisher = {The Open Journal},
      volume = {6},
      number = {59},
      pages = {2977},
      author = {Quinn, Andrew J. and Lopes-dos-Santos, Vitor and Dupret, David and Nobre, Anna C. and Woolrich, Mark W.},
      title = {EMD: Empirical Mode Decomposition and Hilbert-Huang Spectral Analyses in Python},
      journal = {Journal of Open Source Software}
    }

|

The Empirical Mode Decomposition and Hilbert-Huang Spectrum were initially developed by `Norden Huang <https://en.wikipedia.org/wiki/Norden_E._Huang>`_ and colleagues in 1998. This paper is the main reference for the motivation, theory and initial practice of EMD.

.. title image, description
.. raw:: html

    {Huang1998}

.. highlight:: none
::

    @article{Huang1998,
      doi = {10.1098/rspa.1998.0193},
      url = {https://doi.org/10.1098/rspa.1998.0193},
      year = {1998},
      month = mar,
      publisher = {The Royal Society},
      volume = {454},
      number = {1971},
      pages = {903--995},
      author = {Norden E. Huang and Zheng Shen and Steven R. Long and Manli C. Wu and Hsing H. Shih and Quanan Zheng and Nai-Chyuan Yen and Chi Chao Tung and Henry H. Liu},
      title = {The empirical mode decomposition and the Hilbert spectrum for nonlinear and non-stationary time series analysis},
      journal = {Proceedings of the Royal Society of London. Series A: Mathematical,  Physical and Engineering Sciences}
    }


|

Reference Library
-----------------

Depending on the analysis in question, you may also consider citing one or more
of the following papers. The `function docstrings <api.html>`_ in the `EMD`
toolbox contains many citations and references for each method. Here we include
a range of the most relevant papers, grouped by topic.

|

EMD, sifting and the Hilbert-Huang Transform
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

    {Huang1998}

Further EMD Theory
^^^^^^^^^^^^^^^^^^

.. raw:: html

    {Flandrin2004_empirical}

    {Rilling2008_one}


Ensemble Sift
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

    {Wu2009_ensemble}

    {Torres2011_complete}


Masked Sift
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

    {Deering2005_masking}

    {Tsai2016_investigating}

    {Fabus2021_automatic}

Instantaneous Frequency
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

    {Huang2009_instantaneous}


Cycle Analysis and Waveform Shape
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

    {Quinn2021_within}

Holospectrum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

    {Huang2016_holo}
Installing EMD
=================================

There are several ways to install the EMD toolbox. The best one to use depends
on how you want to use the code.


Stable PyPI version
*******************

The `stable version of the code <https://pypi.org/project/emd/>`_ is hosted on `PyPI <https://pypi.org>`_ and will be updated relatively slowly. Any updates to PyPI will (hopefully) only contain working changes that have been running without problems on the development versions of the code for a while.

.. container:: toggle body

    .. container:: header body

        .. raw:: html

            <h3 class='installbar'>install from pip</h3>

    .. container:: installbody body

        EMD can be install from `PyPI <https://pypi.org/project/emd/>`_ using pip::

            pip install emd

        pip will install the latest version of EMD from PyPI alongside any missing dependencies into the current python environment. You can install a specific version by specifying the version number::

            pip install emd==0.5.2


.. container:: toggle body

    .. container:: header body

        .. raw:: html

            <h3 class='installbar'>install in conda environment</h3>

    .. container:: installbody body

        If you want to create a conda environment containing EMD, you can use the following yaml config::

            name: emd
            channels:
            dependencies:
               - pip
               - pip:
                 - emd

        This can be adapted to specify a particular release of EMD by adding the version number to the emd line::

            name: emd
            channels:
            dependencies:
               - pip
               - pip:
                 - emd==0.5.2

        This environment can be customised to include any other packages that you might be working with. The last two lines can also be added to an existing conda environment configuration file to include emd in that env.

        This env can be downloaded `HERE (emd_conda_env.yml) <https://gitlab.com/emd-dev/emd/-/blob/master/envs/emd_conda_env.yml>`_. You can download the config and install the enviromnent by changing directory to the install location and calling these commands::

            curl https://gitlab.com/emd-dev/emd/-/raw/master/envs/emd_conda_env.yml > emd_conda_env.yml
            conda env create -f emd_conda_env.yml

        this will automatically install the required dependancies alongside EMD. The environment can then be activated by calling::

            source activate emd



Development gitlab.com version
******************************

You can also install the `latest development version of EMD
<https://gitlab.com/emd-dev/emd>`_ from gitlab.com using a conda environment.
This version is less stable and likely to change quickly during active
development - however you will get access to new bug-fixes, features and bugs
more quickly.


.. container:: toggle body

    .. container:: header body

        .. raw:: html

            <h3 class='installbar'>install in conda environment</h3>

    .. container:: installbody body

        A conda environment config file can be specified pointing at the development version of EMD on gitlab::

            name: emd
            channels:
            dependencies:
               - pip
               - pip:
                 - git+https://gitlab.com/emd-dev/emd.git

        The env can be downloaded `HERE (emd-dev_conda_env.yml) <https://gitlab.com/emd-dev/emd/-/blob/master/envs/emd-dev_conda_env.yml>`_. You can download the config and install the enviromnent by changing directory to the install location and calling these commands::

            curl https://gitlab.com/emd-dev/emd/-/raw/master/envs/emd-dev_conda_env.yml > emd-dev_conda_env.yml
            conda env create -f emd-dev_conda_env.yml

        this will automatically install the required dependancies alongside EMD. The environment can then be activated by calling::

            source activate emd-dev


.. container:: toggle body

    .. container:: header body

        .. raw:: html

            <h3 class='installbar'>install development branch in conda environment</h3>

    .. container:: installbody body

        A conda environment config file can be specified pointing at the development version of EMD on gitlab. A specific branch can be indicated by adding the branch name after an @ sign in the line specifying the git repo. Here is an example which installs a branch called 'new_feature'::

            name: emd
            channels:
            dependencies:
               - pip
               - pip:
                 - git+https://gitlab.com/emd-dev/emd.git@new_feature

        We provide `an example env here (emd-dev_conda_env.yml) <https://gitlab.com/emd-dev/emd/-/blob/master/envs/emd-dev_conda_env.yml>`_. You can download the config and add the branch name to the right line. Finally, you can install the enviromnent by changing directory to the install location and calling these commands::

            curl https://gitlab.com/emd-dev/emd/-/raw/master/envs/emd-dev_conda_env.yml > emd-dev_conda_env.yml
            conda env create -f emd-dev_conda_env.yml

        this will automatically install the required dependancies alongside EMD. The environment can then be activated by calling::

            source activate emd-dev

.. container:: toggle body

    .. container:: header body

        .. raw:: html

            <h3 class='installbar'>install from source code</h3>

    .. container:: installbody body

        If you plan to actively contribute to EMD, you will need to install EMD directly from source using git. From the terminal, change into the directory you want to install emd into and run the following command::

            cd /home/andrew/src
            git clone https://gitlab.com/emd-dev/emd.git
            cd emd
            python setup.py install

        You will then be able to use git as normal to switch between development branches of EMD and contribute your own.
Citing the EMD Package
=================================
|

If this package has significantly contributed to your work, please include the following citation:

.. title image, description
.. raw:: html

    
    <div class="container" style="margin-bottom:10px; padding-left: 35px; text-indent: -40px">
      Andrew J. Quinn, Vitor Lopes-dos-Santos, David Dupret, Anna Christina Nobre & Mark W. Woolrich (2021)<br>
      <strong><font size="3px">EMD: Empirical Mode Decomposition and Hilbert-Huang Spectral Analyses in Python</font></strong><br>
      Journal of Open Source Software <a href=https://www.doi.org/10.21105/joss.02977>10.21105/joss.02977</a>
    </div>


.. highlight:: none
::

    @article{Quinn2021_joss,
      doi = {10.21105/joss.02977},
      url = {https://doi.org/10.21105/joss.02977},
      year = {2021},
      publisher = {The Open Journal},
      volume = {6},
      number = {59},
      pages = {2977},
      author = {Quinn, Andrew J. and Lopes-dos-Santos, Vitor and Dupret, David and Nobre, Anna C. and Woolrich, Mark W.},
      title = {EMD: Empirical Mode Decomposition and Hilbert-Huang Spectral Analyses in Python},
      journal = {Journal of Open Source Software}
    }

|

The Empirical Mode Decomposition and Hilbert-Huang Spectrum were initially developed by `Norden Huang <https://en.wikipedia.org/wiki/Norden_E._Huang>`_ and colleagues in 1998. This paper is the main reference for the motivation, theory and initial practice of EMD.

.. title image, description
.. raw:: html

    
    <div class="container" style="margin-bottom:10px; padding-left: 35px; text-indent: -40px">
      Norden E. Huang, Zheng Shen, Steven R. Long, Manli C. Wu, Hsing H. Shih, Quanan Zheng, Nai-Chyuan Yen, Chi Chao Tung & Henry H. Liu (1998)<br>
      <strong><font size="3px">The empirical mode decomposition and the Hilbert spectrum for nonlinear and non-stationary time series analysis</font></strong><br>
      Proceedings of the Royal Society of London. Series A: Mathematical,  Physical and Engineering Sciences <a href=https://www.doi.org/10.1098/rspa.1998.0193>10.1098/rspa.1998.0193</a>
    </div>


.. highlight:: none
::

    @article{Huang1998,
      doi = {10.1098/rspa.1998.0193},
      url = {https://doi.org/10.1098/rspa.1998.0193},
      year = {1998},
      month = mar,
      publisher = {The Royal Society},
      volume = {454},
      number = {1971},
      pages = {903--995},
      author = {Norden E. Huang and Zheng Shen and Steven R. Long and Manli C. Wu and Hsing H. Shih and Quanan Zheng and Nai-Chyuan Yen and Chi Chao Tung and Henry H. Liu},
      title = {The empirical mode decomposition and the Hilbert spectrum for nonlinear and non-stationary time series analysis},
      journal = {Proceedings of the Royal Society of London. Series A: Mathematical,  Physical and Engineering Sciences}
    }


|

Reference Library
-----------------

Depending on the analysis in question, you may also consider citing one or more
of the following papers. The `function docstrings <api.html>`_ in the `EMD`
toolbox contains many citations and references for each method. Here we include
a range of the most relevant papers, grouped by topic.

|

EMD, sifting and the Hilbert-Huang Transform
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

    
    <div class="container" style="margin-bottom:10px; padding-left: 35px; text-indent: -40px">
      Norden E. Huang, Zheng Shen, Steven R. Long, Manli C. Wu, Hsing H. Shih, Quanan Zheng, Nai-Chyuan Yen, Chi Chao Tung & Henry H. Liu (1998)<br>
      <strong><font size="3px">The empirical mode decomposition and the Hilbert spectrum for nonlinear and non-stationary time series analysis</font></strong><br>
      Proceedings of the Royal Society of London. Series A: Mathematical,  Physical and Engineering Sciences <a href=https://www.doi.org/10.1098/rspa.1998.0193>10.1098/rspa.1998.0193</a>
    </div>


Further EMD Theory
^^^^^^^^^^^^^^^^^^

.. raw:: html

    
    <div class="container" style="margin-bottom:10px; padding-left: 35px; text-indent: -40px">
      P. Flandrin, G. Rilling & P. Goncalves (2004)<br>
      <strong><font size="3px">Empirical Mode Decomposition as a Filter Bank</font></strong><br>
      IEEE} Signal Processing Letters <a href=https://www.doi.org/10.1109/lsp.2003.821662>10.1109/lsp.2003.821662</a>
    </div>


    
    <div class="container" style="margin-bottom:10px; padding-left: 35px; text-indent: -40px">
      G. Rilling & P. Flandrin (2008)<br>
      <strong><font size="3px">One or Two Frequencies? The Empirical Mode Decomposition Answers</font></strong><br>
      IEEE} Transactions on Signal Processing <a href=https://www.doi.org/10.1109/tsp.2007.906771>10.1109/tsp.2007.906771</a>
    </div>



Ensemble Sift
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

    
    <div class="container" style="margin-bottom:10px; padding-left: 35px; text-indent: -40px">
      Zhaohua Wu & Norden E. Huang (2009)<br>
      <strong><font size="3px">Ensemble Empirical Mode Decomposition: A Noise-Assisted Data Analysis Method</font></strong><br>
      Advances in Adaptive Data Analysis <a href=https://www.doi.org/10.1142/s1793536909000047>10.1142/s1793536909000047</a>
    </div>


    
    <div class="container" style="margin-bottom:10px; padding-left: 35px; text-indent: -40px">
      Maria E. Torres, Marcelo A. Colominas, Gaston Schlotthauer & Patrick Flandrin (2011)<br>
      <strong><font size="3px">A complete ensemble empirical mode decomposition with adaptive noise</font></strong><br>
      2011 IEEE International Conference on Acoustics,  Speech and Signal Processing (ICASSP) <a href=https://www.doi.org/10.1109/icassp.2011.5947265>10.1109/icassp.2011.5947265</a>
    </div>



Masked Sift
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

    
    <div class="container" style="margin-bottom:10px; padding-left: 35px; text-indent: -40px">
      Ryan Deering &  James F. Kaiser (2005)<br>
      <strong><font size="3px">The Use of a Masking Signal to Improve Empirical Mode Decomposition</font></strong><br>
      2005 IEEE International Conference on Acoustics,  Speech and Signal Processing (ICASSP) <a href=https://www.doi.org/10.1109/icassp.2005.1416051>10.1109/icassp.2005.1416051</a>
    </div>


    
    <div class="container" style="margin-bottom:10px; padding-left: 35px; text-indent: -40px">
      Feng-Fang Tsai, Shou-Zen Fan, Yi-Shiuan Lin, Norden E. Huang & Jia-Rong Yeh (2016)<br>
      <strong><font size="3px">Investigating Power Density and the Degree of Nonlinearity in Intrinsic Components of Anesthesia EEG by the Hilbert-Huang Transform: An Example Using Ketamine and Alfentanil</font></strong><br>
      PLOS-ONE <a href=https://www.doi.org/10.1371/journal.pone.0168108>10.1371/journal.pone.0168108</a>
    </div>


    
    <div class="container" style="margin-bottom:10px; padding-left: 35px; text-indent: -40px">
      Marco S. Fabus, Andrew J. Quinn, Catherine E. Warnaby & Mark W. Woolrich (2021)<br>
      <strong><font size="3px">Automatic decomposition of electrophysiological data into distinct nonsinusoidal oscillatory modes</font></strong><br>
      Journal of Neurophysiology <a href=https://www.doi.org/10.1152/jn.00315.2021>10.1152/jn.00315.2021</a>
    </div>


Instantaneous Frequency
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

    
    <div class="container" style="margin-bottom:10px; padding-left: 35px; text-indent: -40px">
      Norden E. Huang, Zhaohua Wu, Steven R. Long, Kenneth C. Arnold, Xianyao Chen & Karin Blank (2009)<br>
      <strong><font size="3px">On Instantaneous Frequency</font></strong><br>
      Advances in Adaptive Data Analysis <a href=https://www.doi.org/10.1142/s1793536909000096>10.1142/s1793536909000096</a>
    </div>



Cycle Analysis and Waveform Shape
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

    
    <div class="container" style="margin-bottom:10px; padding-left: 35px; text-indent: -40px">
      Quinn, Andrew J., Lopes-dos-Santos, Vítor, Huang, Norden, Liang, Wei-Kuang, Juan, Chi-Hung, Yeh, Jia-Rong, Nobre, Anna C., Dupret, David, Woolrich & Mark W. (2021)<br>
      <strong><font size="3px">Within-cycle instantaneous frequency profiles report oscillatory waveform dynamics</font></strong><br>
      Journal of Neurophysiology <a href=https://www.doi.org/10.1152/jn.00201.2021>10.1152/jn.00201.2021</a>
    </div>


Holospectrum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. raw:: html

    
    <div class="container" style="margin-bottom:10px; padding-left: 35px; text-indent: -40px">
      Norden E. Huang, Kun Hu, Albert C. C. Yang, Hsing-Chih Chang, Deng Jia, Wei-Kuang Liang, Jia Rong Yeh, Chu-Lan Kao, Chi-Hung Juan, Chung Kang Peng, Johanna H. Meijer, Yung-Hung Wang, Steven R. Long & Zhauhua Wu (2016)<br>
      <strong><font size="3px">On Holo-Hilbert spectral analysis: a full informational spectral representation for nonlinear and non-stationary data</font></strong><br>
      Philosophical Transactions of the Royal Society A: Mathematical,  Physical and Engineering Sciences <a href=https://www.doi.org/10.1098/rsta.2015.0206>10.1098/rsta.2015.0206</a>
    </div>

API
================


Sift Functions
*********************

Primary user-level functions for running the sift.

.. autosummary::
     :toctree: stubs

     emd.sift.sift
     emd.sift.ensemble_sift
     emd.sift.complete_ensemble_sift
     emd.sift.mask_sift
     emd.sift.iterated_mask_sift
     emd.sift.sift_second_layer
     emd.sift.mask_sift_second_layer


Sift Utilities
*********************

Low-level utility functions used by the sift routines.

.. autosummary::
     :toctree: stubs

     emd.sift.get_config
     emd.sift.get_next_imf
     emd.sift.get_next_imf_mask
     emd.sift.interp_envelope
     emd.sift.get_padded_extrema
     emd.sift.fixed_stop
     emd.sift.sd_stop
     emd.sift.rilling_stop
     emd.sift.energy_stop
     emd.sift.is_imf

Frequency Functions
*********************

Computing frequency transforms from narrow band oscillations (IMFs).

.. autosummary::
     :toctree: stubs

     emd.spectra.frequency_transform
     emd.spectra.phase_from_complex_signal
     emd.spectra.freq_from_phase
     emd.spectra.phase_from_freq

Spectrum Functions
*********************

Compute Hilbert-Huang and Holospectra from instantaneous frequency data.

.. autosummary::
     :toctree: stubs

     emd.spectra.hilberthuang
     emd.spectra.holospectrum

Spectrum Utilities
*********************

Low-level helper functions for spectrum computations.

.. autosummary::
     :toctree: stubs

     emd.spectra.define_hist_bins
     emd.spectra.define_hist_bins_from_data

Cycle Analysis
*********************

Identify and analyse single cycles of an oscillation.

.. autosummary::
     :toctree: stubs

     emd.cycles.Cycles
     emd.cycles.get_cycle_vector
     emd.cycles.get_cycle_stat
     emd.cycles.get_control_points
     emd.cycles.phase_align
     emd.cycles.normalised_waveform
     emd.cycles.bin_by_phase
     emd.cycles.mean_vector
     emd.cycles.kdt_match

Package Utilities
*********************

Routines related to python, logging and installation.

.. autosummary::
     :toctree: stubs

     emd.support.get_install_dir
     emd.support.get_installed_version
     emd.logger.set_up
