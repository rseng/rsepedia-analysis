# Elephant - Electrophysiology Analysis Toolkit

[![Build Status](https://travis-ci.org/NeuralEnsemble/elephant.svg?branch=master)](https://travis-ci.org/NeuralEnsemble/elephant)
[![Coverage Status](https://coveralls.io/repos/github/NeuralEnsemble/elephant/badge.svg?branch=master)](https://coveralls.io/github/NeuralEnsemble/elephant?branch=master)
[![Documentation Status](https://readthedocs.org/projects/elephant/badge/?version=latest)](https://elephant.readthedocs.io/en/latest/?badge=latest)
[![![PyPi]](https://img.shields.io/pypi/v/elephant)](https://pypi.org/project/elephant/)
[![Statistics](https://img.shields.io/pypi/dm/elephant)](https://seladb.github.io/StarTrack-js/#/preload?r=neuralensemble,elephant)
[![Gitter](https://badges.gitter.im/python-elephant/community.svg)](https://gitter.im/python-elephant/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

*Elephant* package analyses all sorts of neurophysiological data:
spike trains, LFP, analog signals. The input-output data format is either
[Neo](https://github.com/NeuralEnsemble/python-neo), Quantity or Numpy array.


### More Information

* Documentation: https://elephant.readthedocs.io/en/latest/
* Mailing list: https://groups.google.com/group/neuralensemble


#### Visualization of Elephant analysis objects

Viziphant package (https://github.com/INM-6/viziphant) is developed by Elephant
team for plotting and visualization of the output of Elephant functions in a
few lines of code.


#### License
 
Modified BSD License, see [LICENSE.txt](LICENSE.txt) for details.


#### Copyright

:copyright: 2014-2021 by the [Elephant team](doc/authors.rst).


#### Acknowledgments

See [acknowledgments](doc/acknowledgments.rst).


#### Citation

See [citations](doc/citation.rst).

Here, are CSD methods for different electrode configurations.

Keywords: Local field potentials; Current-source density; CSD;
Multielectrode; Laminar electrode; Barrel cortex

1D - laminar probe like electrodes. 
2D - Microelectrode Array like
3D - UtahArray or multiple laminar probes.

The following methods have been implemented here, for, 

1D - StandardCSD, DeltaiCSD, SplineiCSD, StepiCSD, KCSD1D
2D - KCSD2D, MoIKCSD (Saline layer on top of slice)
3D - KCSD3D

Each of these methods listed have some advantages - except StandardCSD which is
not recommended. The KCSD methods can handle broken or irregular electrode
configurations electrode

iCSD
----
Python-implementation of the inverse current source density (iCSD) methods from
http://software.incf.org/software/csdplotter

The Python iCSD toolbox lives on GitHub as well:
https://github.com/espenhgn/iCSD

The methods were originally developed by Klas H. Pettersen, as described in:
Klas H. Pettersen, Anna Devor, Istvan Ulbert, Anders M. Dale, Gaute
T. Einevoll, Current-source density estimation based on inversion of
electrostatic forward solution: Effects of finite extent of neuronal activity
and conductivity discontinuities, Journal of Neuroscience Methods, Volume 154,
Issues 1Ð2, 30 June 2006, Pages 116-133, ISSN 0165-0270,
http://dx.doi.org/10.1016/j.jneumeth.2005.12.005.
(http://www.sciencedirect.com/science/article/pii/S0165027005004541)

To see an example of usage of the methods, see
[demo_icsd.py](https://github.com/espenhgn/iCSD/blob/master/demo_icsd.py)

KCSD 
---- 
This is 1.0 version of kCSD inverse method proposed in

J. Potworowski, W. Jakuczun, S. Łęski, D. K. Wójcik
"Kernel Current Source Density Method"
Neural Computation 24 (2012), 541–575

Some key advantages for KCSD methods are
-- irregular grid of electrodes - accepts arbitrary electrode placement.
-- crossvalidation to ensure no over fitting
-- CSD is not limited to electrode positions - it can obtained at any location

For citation purposes, 
If you use this software in published research please cite the following work
- kCSD1D - [1, 2]
- kCSD2D - [1, 3]
- kCSD3D - [1, 4]
- MoIkCSD - [1, 3, 5]

[1] Potworowski, J., Jakuczun, W., Łęski, S. & Wójcik, D. (2012) 'Kernel
current source density method.' Neural Comput 24(2), 541-575.

[2] Pettersen, K. H., Devor, A., Ulbert, I., Dale, A. M. & Einevoll,
G. T. (2006) 'Current-source density estimation based on inversion of
electrostatic forward solution: effects of finite extent of neuronal activity
and conductivity discontinuities.' J Neurosci Methods 154(1-2), 116-133.

[3] Łęski, S., Pettersen, K. H., Tunstall, B., Einevoll, G. T., Gigg, J. &
Wójcik, D. K. (2011) 'Inverse Current Source Density method in two dimensions:
Inferring neural activation from multielectrode recordings.' Neuroinformatics
9(4), 401-425.

[4] Łęski, S., Wójcik, D. K., Tereszczuk, J., Świejkowski, D. A., Kublik, E. &
Wróbel, A. (2007) 'Inverse current-source density method in 3D: reconstruction
fidelity, boundary effects, and influence of distant sources.' Neuroinformatics
5(4), 207-222.

[5] Ness, T. V., Chintaluri, C., Potworowski, J., Łeski, S., Głabska, H.,
Wójcik, D. K. & Einevoll, G. T. (2015) 'Modelling and Analysis of Electrical
Potentials Recorded in Microelectrode Arrays (MEAs).' Neuroinformatics 13(4),
403-426.

For your research interests of Kernel methods of CSD please see,
https://github.com/Neuroinflab/kCSD-python 

Contact: Prof. Daniel K. Wojcik

Here (https://github.com/Neuroinflab/kCSD-python/tree/master/tests), are
scripts to compare different KCSD methods with different CSD sources. You can
play around with the different parameters of the methods.

The implentation is based on the Matlab version at INCF
(http://software.incf.org/software/kcsd), which is now out-dated. A python
version based on this was developed by Grzegorz Parka
(https://github.com/INCF/pykCSD), which is also not supported at this
point. This current version of KCSD methods in elephant is a mirror of
https://github.com/Neuroinflab/kCSD-python/commit/8e2ae26b00da7b96884f2192ec9ea612b195ec30
---
name: Bug report
about: Report python-elephant bugs/errors here
title: "[Bug] "
labels: ''
assignees: ''

---

**Describe the bug**
<!-- A clear and concise description of what the bug is. -->

**To Reproduce**
1. 
2.

<!-- If you have a code sample, error messages, or stack traces, please provide it here as well. -->


**Expected behavior**
<!-- A clear and concise description of what you expected to happen. -->

**Environment**
 - OS (e.g., Linux):
 - How you installed elephant (`conda`, `pip`, source):
 - Python version:
 - `neo` python package version: 
 - `elephant` python package version:
 - (Any additional python package you want to include here)
---
name: Questions/Help/Support
about: Do you have a question or need help to start running elephant?
title: ''
labels: ''
assignees: ''

---

**Questions and Help**

_Please note that this issue tracker is not a help form and this issue may be closed and redirected to a more appropriate communication channel. Before asking a new question, please check our [NeuralEnsembe forum](https://groups.google.com/forum/#!forum/neuralensemble), where you can discuss any NeuralEnsemble resource, including Elephant._
---
name: Feature request
about: Suggest an idea for this project
title: "[Feature] "
labels: ''
assignees: ''

---

**Feature**
<!-- A clear and concise description of the feature proposal. -->

**Motivation**
<!-- Please outline the motivation for the proposal. Is your feature request related to a problem? E.g., I'm always frustrated when [...]. If this problem is related to another GitHub issue, please link here too. -->

**Alternatives**
<!-- A clear and concise description of any alternative solutions or features you've considered, if any. -->

**Additional context**
<!-- Add any other context about the feature request here. This may be screenshots, online resources, potential experts to contact, scientific publications, demos, etc.-->
=============
Release Notes
=============

Elephant 0.10.0 release notes
=============================

Documentation
-------------
The documentation is revised and restructured by categories (https://github.com/NeuralEnsemble/elephant/pull/386) to simplify navigation on readthedocs and improve user experience. All citations used in Elephant are stored in a single [BibTex file](https://github.com/NeuralEnsemble/elephant/blob/master/doc/bib/elephant.bib).

Optimizations
-------------

CUDA and OpenCL support
***********************
[Analysis of Sequences of Synchronous EvenTs](https://elephant.readthedocs.io/en/latest/reference/asset.html) has become the first module in Elephant that supports CUDA and OpenCL (https://github.com/NeuralEnsemble/elephant/pull/351, https://github.com/NeuralEnsemble/elephant/pull/404, https://github.com/NeuralEnsemble/elephant/pull/399). Whether you have an Nvidia GPU or just run the analysis on a laptop with a built-in Intel graphics card, the speed-up is **X100** and **X1000** compared to a single CPU core. The computations are optimized to a degree that you can analyse and look for spike patterns in real data in several minutes of compute time on a laptop. The installation instructions are described in the [install](https://elephant.readthedocs.io/en/latest/install.html) section.

Other optimizations
*******************
* Surrogates: sped up bin shuffling (https://github.com/NeuralEnsemble/elephant/pull/400) and reimplemented the continuous time version (https://github.com/NeuralEnsemble/elephant/pull/397)
* Improved memory efficiency of creating a BinnedSpikeTrain (https://github.com/NeuralEnsemble/elephant/pull/395)

New functionality and features
------------------------------
* Synchrofact detection (https://github.com/NeuralEnsemble/elephant/pull/322) is a method to detect highly synchronous spikes (at the level of sampling rate precision with an option to extend this to jittered synchrony) and annotate or optionally remove them.
* Added `phase_locking_value`, `mean_phase_vector`, and `phase_difference` functions (https://github.com/NeuralEnsemble/elephant/pull/385/files)
* BinnedSpikeTrain:
  - added `to_spike_trains` and `time_slice` functions (https://github.com/NeuralEnsemble/elephant/pull/390). Now you can slice a binned spike train as `bst[:, i:j]` or `bst.time_slice(t_start, t_stop)`. Also, with `to_spike_trains` function, you can generate a realization of spike trains that maps to the same BinnedSpikeTrain object when binned.
  - optional CSC format (https://github.com/NeuralEnsemble/elephant/pull/402)
  - the `copy` parameter (False by default) in the `binarize` function makes a *shallow* copy, if set to True, of the output BinnedSpikeTrain object (https://github.com/NeuralEnsemble/elephant/pull/402)
* Granger causality tutorial notebook (https://github.com/NeuralEnsemble/elephant/pull/393)
* Unitary Event Analysis support multiple pattern hashes (https://github.com/NeuralEnsemble/elephant/pull/387)

Bug fixes
---------
* Account for unidirectional spiketrain->segment links in synchrofact deletion (https://github.com/NeuralEnsemble/elephant/pull/398)
* Joint-ISI dithering: fixed a bug regarding first ISI bin (https://github.com/NeuralEnsemble/elephant/pull/396)
* Fix LvR values from being off when units are in seconds (https://github.com/NeuralEnsemble/elephant/pull/389)


Elephant 0.9.0 release notes
****************************

This release is titled to accompany the [2nd Elephant User Workshop](https://www.humanbrainproject.eu/en/education/participatecollaborate/infrastructure-events-trainings/2nd-elephant-user-workshop/)

Viziphant
---------
Meet Viziphant, the visualization of Elephant analysis methods, at https://viziphant.readthedocs.io/en/latest/. This package provides support to easily plot and visualize the output of Elephant functions in a few lines of code.

Provenance tracking
-------------------
Provenance is becoming a separate direction in Elephant. Many things are still to come, and we started with annotating `time_histogram`, `instantaneous_rate` and `cross_correlation_histogram` outputs to carry the information about the parameters these functions used. This allowed Viziphant, the visualization of Elephant analyses, to look for the `.annotations` dictionary of the output of these function to "understand" how the object has been generated and label the plot axes accordingly.

New functionality and features
------------------------------
* Time-domain pairwise and conditional pairwise Granger causality measures (https://github.com/NeuralEnsemble/elephant/pull/332, https://github.com/NeuralEnsemble/elephant/pull/359)
* Spike contrast function that measures the synchrony of spike trains (https://github.com/NeuralEnsemble/elephant/pull/354; thanks to @Broxy7 for bringing this in Elephant).
* Revised local variability LvR (https://github.com/NeuralEnsemble/elephant/pull/346) as an alternative to the LV measure.
* Three surrogate methods: Trial-shifting, Bin Shuffling, ISI dithering (https://github.com/NeuralEnsemble/elephant/pull/343).
* Added a new function to generate spike trains: `inhomogeneous_gamma_process` (https://github.com/NeuralEnsemble/elephant/pull/339).
* The output of `instantaneous_rate` function is now a 2D matrix of shape `(time, len(spiketrains))` (https://github.com/NeuralEnsemble/elephant/issues/363). Not only can the users assess the averaged instantaneous rate (`rates.mean(axis=1)`) but also explore how much the instantaneous rate deviates from trial to trial (`rates.std(axis=1)`) (originally asked in https://github.com/NeuralEnsemble/elephant/issues/363).

Python 3 only
-------------
* Python 2.7 and 3.5 support is dropped. You can still however enjoy the features of Elephant v0.9.0 with Python 2.7 or 3.5 by installing Elephant from [this](https://github.com/NeuralEnsemble/elephant/tree/295c6bd7fea196cf9665a78649fafedab5840cfa) commit `pip install git+https://github.com/NeuralEnsemble/elephant@295c6bd7fea196cf9665a78649fafedab5840cfa#egg=elephant[extras]`
* Added Python 3.9 support.

Optimization
------------
* You have been asking for direct numpy support for years. Added `_t_start`, `_t_stop`, and `_bin_size` attributes of BinnedSpikeTrain are guaranteed to be of the same units and hence are unitless (https://github.com/NeuralEnsemble/elephant/pull/378). It doesn't mean though that you need to care about units on your own: `t_start`, `t_stop`, and `bin_size` properties are still quantities with units. The `.rescale()` method of a BinnedSpikeTrain rescales the internal units to new ones in-place. The following Elephant functions are optimized with unitless BinnedSpikeTrain:
  - cross_correlation_histogram
  - bin_shuffling (one of the surrogate methods)
  - spike_train_timescale
* X4 faster binning and overall BinnedSpikeTrain object creation (https://github.com/NeuralEnsemble/elephant/pull/368).
* `instantaneous_rate` function is vectorized to work with a list of spike train trials rather than computing them in a loop (previously, `for spiketrain in spiketrains; do compute instantaneous_rate(spiketrain); done`), which brought X25 speedup (https://github.com/NeuralEnsemble/elephant/pull/362; thanks to @gyyang for the idea and original implementation).
* Memory-efficient `zscore` function (https://github.com/NeuralEnsemble/elephant/pull/372).
* Don't sort the input array in ISI function (https://github.com/NeuralEnsemble/elephant/pull/371), which reduces function algorithmic time complexity from `O(N logN)` to linear `O(N)`. Now, when the input time array is not sorted, a warning is shown.
* Vectorized Current Source Density `generate_lfp` function (https://github.com/NeuralEnsemble/elephant/pull/358).

Breaking changes
----------------
* mpi4py package is removed from the extra requirements to allow `pip install elephant[extras]` on machines without MPI installed system-wide. Refer to [MPI support](https://elephant.readthedocs.io/en/latest/install.html#mpi-support) installation page in elephant.
* BinnedSpikeTrain (https://github.com/NeuralEnsemble/elephant/pull/368, https://github.com/NeuralEnsemble/elephant/pull/377):
  - previously, when t_start/stop, if set manually, was outside of the shared time interval, only the shared [t_start_shared=max(t_start), t_stop_shared=min(t_stop)] interval was implicitly considered without any warnings. Now an error is thrown with a description on how to fix it.
  - removed `lst_input`, `input_spiketrains`, `matrix_columns`, `matrix_rows` (in favor of the new attribute - `shape`), `tolerance`, `is_spiketrain`, `is_binned` attributes from BinnedSpikeTrain class. Part of them are confusing (e.g., `is_binned` was just the opposite of `is_spiketrain`, but one can erroneously think that it's data is clipped to 0 and 1), and part of them - `lst_input`, `input_spiketrains` input data - should not have been saved as attributes of an object in the first place because the input spike trains are not used after the sparse matrix is created.
  - now the users can directly access `.sparse_matrix` attribute of BinnedSpikeTrain to do efficient (yet unsafe in general) operations. For this reason, `to_sparse_array()` function, which does not make a copy, as one could think of, is deprecated.
* `instantaneous_rate` function (https://github.com/NeuralEnsemble/elephant/pull/362):
  - in case of multiple input spike trains, the output of the instantaneous rate function is (always) a 2D matrix of shape `(time, len(spiketrains))` instead of a pseudo 1D array (previous behavior) of shape `(time, 1)` that contained the instantaneous rate summed across input spike trains;
  - in case of multiple input spike trains, the user needs to manually provide the input kernel instead of `auto`, which is set by default, for the reason that it's currently not clear how to estimate the common kernel for a set of spike trains. If you have an idea how to do this, we`d appreciate if you let us know by [getting in touch with us](https://elephant.readthedocs.io/en/latest/get_in_touch.html).

Other changes
-------------
* `waveform_snr` function now directly takes a 2D or 3D waveforms matrix rather than a spike train (deprecated behavior).
* Added a warning in fanofactor function when the input spiketrains vary in their durations (https://github.com/NeuralEnsemble/elephant/pull/341).
* SPADE: New way to count patterns for multiple testing (https://github.com/NeuralEnsemble/elephant/pull/347)
* GPFA renamed 'xsm' -> 'latent_variable' and 'xorth' -> 'latent_variable_orth'

Bug fixes
---------
* Instantaneous rate arrays were not centered at the origin for spike trains that are symmetric at t=0 with `center_kernel=True` option (https://github.com/NeuralEnsemble/elephant/pull/362).
* The number of discarded spikes that fall into the last bin of a BinnedSpikeTrain object was incorrectly calculated (https://github.com/NeuralEnsemble/elephant/pull/368).
* Fixed index selection in `spike_triggered_phase` (https://github.com/NeuralEnsemble/elephant/pull/382)
* Fixed surrogates bugs:
  - `joint-ISI` and `shuffle ISI` output spike trains were not sorted in time (https://github.com/NeuralEnsemble/elephant/pull/364);
  - surrogates get arbitrary sampling_rate (https://github.com/NeuralEnsemble/elephant/pull/353), which relates to the provenance tracking issue;



Elephant 0.8.0 release notes
****************************

New features
------------
* The `parallel` module is a new experimental module (https://github.com/NeuralEnsemble/elephant/pull/307) to run python functions concurrently. Supports native (pythonic) ProcessPollExecutor and MPI. Not limited to Elephant functional.
* Added an optional `refractory_period` argument, set to None by default, to `dither_spikes` function (https://github.com/NeuralEnsemble/elephant/pull/297).
* Added `cdf` and `icdf` functions in Kernel class to correctly estimate the median index, needed for `instantaneous_rate` function in statistics.py (https://github.com/NeuralEnsemble/elephant/pull/313).
* Added an optional `center_kernel` argument, set to True by default (to behave as in Elephant <0.8.0 versions) to `instantaneous_rate` function in statistics.py (https://github.com/NeuralEnsemble/elephant/pull/313).

New tutorials
-------------
* Analysis of Sequences of Synchronous EvenTs (ASSET) tutorial: https://elephant.readthedocs.io/en/latest/tutorials/asset.html
* Parallel module tutorial: https://elephant.readthedocs.io/en/latest/tutorials/parallel.html

Optimization
------------
* Optimized ASSET runtime by a factor of 10 and more (https://github.com/NeuralEnsemble/elephant/pull/259, https://github.com/NeuralEnsemble/elephant/pull/333).

Python 2.7 and 3.5 deprecation
------------------------------
Python 2.7 and 3.5 are deprecated and will not be maintained by the end of 2020. Switch to Python 3.6+.

Breaking changes
----------------
* Naming convention changes (`binsize` -> `bin_size`, etc.) in almost all Elephant functions (https://github.com/NeuralEnsemble/elephant/pull/316).

Elephant 0.7.0 release notes
****************************

Breaking changes
----------------
* [gpfa] GPFA dimensionality reduction method is rewritten in easy-to-use scikit-learn class style format (https://github.com/NeuralEnsemble/elephant/pull/287):

.. code-block:: python

    gpfa = GPFA(bin_size=20*pq.ms, x_dim=8)
    results = gpfa.fit_transform(spiketrains, returned_data=['xorth', 'xsm'])

New tutorials
-------------
* GPFA dimensionality reduction method: https://elephant.readthedocs.io/en/latest/tutorials/gpfa.html
* Unitary Event Analysis of coordinated spiking activity: https://elephant.readthedocs.io/en/latest/tutorials/unitary_event_analysis.html
* (Introductory) statistics module: https://elephant.readthedocs.io/en/latest/tutorials/statistics.html

Deprecations
------------
* **Python 2.7 support will be dropped on Dec 31, 2020.** Please switch to Python 3.6, 3.7, or 3.8.
* [spike train generation] `homogeneous_poisson_process_with_refr_period()`, introduced in v0.6.4, is deprecated and will be deleted in v0.8.0. Use `homogeneous_poisson_process(refractory_period=...)` instead.
* [pandas bridge] pandas\_bridge module is deprecated and will be deleted in v0.8.0.

New features
------------
* New documentation style, guidelines, tutorials, and more (https://github.com/NeuralEnsemble/elephant/pull/294).
* Python 3.8 support (https://github.com/NeuralEnsemble/elephant/pull/282).
* [spike train generation] Added `refractory_period` flag in `homogeneous_poisson_process()` (https://github.com/NeuralEnsemble/elephant/pull/292) and `inhomogeneous_poisson_process()` (https://github.com/NeuralEnsemble/elephant/pull/295) functions. The default is `refractory_period=None`, meaning no refractoriness.
* [spike train correlation] `cross_correlation_histogram()` supports different t_start and t_stop of input spiketrains.
* [waveform features] `waveform_width()` function extracts the width (trough-to-peak TTP) of a waveform (https://github.com/NeuralEnsemble/elephant/pull/279).
* [signal processing] Added `scaleopt` flag in `pairwise_cross_correlation()` to mimic the behavior of Matlab's `xcorr()` function (https://github.com/NeuralEnsemble/elephant/pull/277). The default is `scaleopt=unbiased` to be consistent with the previous versions of Elephant.
* [spike train surrogates] Joint-ISI dithering method via `JointISI` class (https://github.com/NeuralEnsemble/elephant/pull/275).

Bug fixes
---------
* [spike train correlation] Fix CCH Border Correction (https://github.com/NeuralEnsemble/elephant/pull/298). Now, the border correction in `cross_correlation_histogram()` correctly reflects the number of bins used for the calculation at each lag. The correction factor is now unity at full overlap.
* [phase analysis] `spike_triggered_phase()` incorrect behavior when the spike train and the analog signal had different time units (https://github.com/NeuralEnsemble/elephant/pull/270).

Performance
-----------
* [spade] SPADE x7 speedup (https://github.com/NeuralEnsemble/elephant/pull/280, https://github.com/NeuralEnsemble/elephant/pull/285, https://github.com/NeuralEnsemble/elephant/pull/286). Moreover, SPADE is now able to handle all surrogate types that are available in Elephant, as well as more types of statistical corrections.
* [conversion] Fast & memory-efficient `covariance()` and Pearson `corrcoef()` (https://github.com/NeuralEnsemble/elephant/pull/274). Added flag `fast=True` by default in both functions.
* [conversion] Use fast fftconvolve instead of np.correlate in `cross_correlation_histogram()` (https://github.com/NeuralEnsemble/elephant/pull/273).


Elephant 0.6.4 release notes
****************************

This release has been made for the "1st Elephant User Workshop" (https://www.humanbrainproject.eu/en/education/participatecollaborate/infrastructure-events-trainings/1st-elephant-user-workshop-accelerate-structured-and-reproducibl).


Main features
-------------
* neo v0.8.0 compatible


New modules
-----------
* GPFA - Gaussian-process factor analysis - dimensionality reduction method for neural trajectory visualization (https://github.com/NeuralEnsemble/elephant/pull/233). _Note: the API could change in the future._


Bug fixes
---------
* [signal processing] Keep `array_annotations` in the output of signal processing functions (https://github.com/NeuralEnsemble/elephant/pull/258).
* [SPADE] Fixed the calculation of the duration of a pattern in the output (https://github.com/NeuralEnsemble/elephant/pull/254).
* [statistics] Fixed automatic kernel selection yields incorrect values (https://github.com/NeuralEnsemble/elephant/pull/246).


Improvements
------------
* Vectorized `spike_time_tiling_coefficient()` function - got rid of a double for-loop (https://github.com/NeuralEnsemble/elephant/pull/244)
* Reduced the number of warnings during the tests (https://github.com/NeuralEnsemble/elephant/pull/238).
* Removed unused debug code in `spade/fast_fca.py` (https://github.com/NeuralEnsemble/elephant/pull/249).
* Improved doc string of `covariance()` and `corrcoef()` (https://github.com/NeuralEnsemble/elephant/pull/260).



Elephant 0.6.3 release notes
****************************
July 22nd 2019

The release v0.6.3 is mostly about improving maintenance.

New functions
-------------
* `waveform_features` module
    * Waveform signal-to-noise ratio (https://github.com/NeuralEnsemble/elephant/pull/219).
* Added support for Butterworth `sosfiltfilt` - numerically stable (in particular, higher order) filtering (https://github.com/NeuralEnsemble/elephant/pull/234).

Bug fixes
---------
* Fixed neo version typo in requirements file (https://github.com/NeuralEnsemble/elephant/pull/218)
* Fixed broken docs (https://github.com/NeuralEnsemble/elephant/pull/230, https://github.com/NeuralEnsemble/elephant/pull/232)
* Fixed issue with 32-bit arch (https://github.com/NeuralEnsemble/elephant/pull/229)

Other changes
-------------
* Added issue templates (https://github.com/NeuralEnsemble/elephant/pull/226)
* Single VERSION file (https://github.com/NeuralEnsemble/elephant/pull/231)

Elephant 0.6.2 release notes
****************************
April 23rd 2019

New functions
-------------
* `signal_processing` module
    * New functions to calculate the area under a time series and the derivative of a time series.

Other changes
-------------
* Added support to initialize binned spike train representations with a matrix
* Multiple bug fixes


Elephant 0.6.1 release notes
****************************
April 1st 2019

New functions
-------------
* `signal_processing` module
    * New function to calculate the cross-correlation function for analog signals.
* `spade` module
    * Spatio-temporal spike pattern detection now includes the option to assess significance also based on time-lags of patterns, in addition to patterns size and frequency (referred to as 3D pattern spectrum).

Other changes
-------------
* This release fixes a number of compatibility issues in relation to API breaking changes in the Neo library.
* Fixed error in STTC calculation (spike time tiling coefficient)
* Minor bug fixes


Elephant 0.6.0 release notes
****************************
October 12th 2018

New functions
-------------
* `cell_assembly_detection` module
    * New function to detect higher-order correlation structures such as patterns in parallel spike trains based on Russo et al, 2017.
*  **wavelet_transform()** function in `signal_prosessing.py` module
    * Function for computing wavelet transform of a given time series based on Le van Quyen et al. (2001)

Other changes
-------------
* Switched to multiple `requirements.txt` files which are directly read into the `setup.py`
* `instantaneous_rate()` accepts now list of spiketrains
* Minor bug fixes


Elephant 0.5.0 release notes
****************************
April 4nd 2018

New functions
-------------
* `change_point_detection` module:
    * New function to detect changes in the firing rate
* `spike_train_correlation` module:
    * New function to calculate the spike time tiling coefficient
* `phase_analysis` module:
    * New function to extract spike-triggered phases of an AnalogSignal
* `unitary_event_analysis` module:
    * Added new unit test to the UE function to verify the method based on data of a recent [Re]Science publication

Other changes
-------------
* Minor bug fixes


Elephant 0.4.3 release notes
****************************
March 2nd 2018

Other changes
-------------
* Bug fixes in `spade` module:
    * Fixed an incompatibility with the latest version of an external library


Elephant 0.4.2 release notes
****************************
March 1st 2018

New functions
-------------
* `spike_train_generation` module:
    * **inhomogeneous_poisson()** function
* Modules for Spatio Temporal Pattern Detection (SPADE) `spade_src`:
    * Module SPADE: `spade.py`
* Module `statistics.py`:
    * Added CV2 (coefficient of variation for non-stationary time series)
* Module `spike_train_correlation.py`:
    * Added normalization in **cross-correlation histogram()** (CCH)

Other changes
-------------
* Adapted the `setup.py` to automatically install the spade modules including the compiled `C` files `fim.so`
* Included testing environment for MPI in `travis.yml`
* Changed function arguments  in `current_source_density.py` to `neo.AnalogSignal` instead list of `neo.AnalogSignal` objects
* Fixes to travis and setup configuration files
* Fixed bug in ISI function `isi()`, `statistics.py` module
* Fixed bug in `dither_spikes()`, `spike_train_surrogates.py`
* Minor bug fixes


Elephant 0.4.1 release notes
****************************
March 23rd 2017

Other changes
-------------
* Fix in `setup.py` to correctly import the current source density module


Elephant 0.4.0 release notes
****************************
March 22nd 2017

New functions
-------------
* `spike_train_generation` module:
    * peak detection: **peak_detection()**
* Modules for Current Source Density: `current_source_density_src`
    * Module Current Source Density: `KCSD.py`
    * Module for Inverse Current Source Density: `icsd.py`

API changes
-----------
* Interoperability between Neo 0.5.0 and Elephant
    * Elephant has adapted its functions to the changes in Neo 0.5.0,
      most of the functionality behaves as before
    * See Neo documentation for recent changes: http://neo.readthedocs.io/en/latest/whatisnew.html

Other changes
-------------
* Fixes to travis and setup configuration files.
* Minor bug fixes.
* Added module `six` for Python 2.7 backwards compatibility


Elephant 0.3.0 release notes
****************************
April 12st 2016

New functions
-------------
* `spike_train_correlation` module:
    * cross correlation histogram: **cross_correlation_histogram()**
* `spike_train_generation` module:
    * single interaction process (SIP): **single_interaction_process()**
    * compound Poisson process (CPP): **compound_poisson_process()**
* `signal_processing` module:
    * analytic signal: **hilbert()**
* `sta` module:
    * spike field coherence: **spike_field_coherence()**
* Module to represent kernels: `kernels` module
* Spike train metrics / dissimilarity / synchrony measures: `spike_train_dissimilarity` module
* Unitary Event (UE) analysis: `unitary_event_analysis` module
* Analysis of Sequences of Synchronous EvenTs (ASSET): `asset` module

API changes
-----------
* Function **instantaneous_rate()** now uses kernels as objects defined in the `kernels` module. The previous implementation of the function using the `make_kernel()` function is deprecated, but still temporarily available as `oldfct_instantaneous_rate()`.

Other changes
-------------
* Fixes to travis and readthedocs configuration files.


Elephant 0.2.1 release notes
****************************
February 18th 2016

Other changes
-------------
Minor bug fixes.


Elephant 0.2.0 release notes
****************************
September 22nd 2015

New functions
-------------
* Added covariance function **covariance()** in the `spike_train_correlation` module
* Added complexity pdf **complexity_pdf()** in the `statistics` module
* Added spike train extraction from analog signals via threshold detection the in `spike_train_generation` module
* Added **coherence()** function for analog signals in the `spectral` module
* Added **Cumulant Based Inference for higher-order of Correlation (CuBIC)** in the `cubic` module for correlation analysis of parallel recorded spike trains

API changes
-----------
* **Optimized kernel bandwidth** in `rate_estimation` function: Calculates the optimized kernel width when the paramter kernel width is specified as `auto`

Other changes
-------------
* **Optimized creation of sparse matrices**: The creation speed of the sparse matrix inside the `BinnedSpikeTrain` class is optimized
* Added **Izhikevich neuron simulator** in the `make_spike_extraction_test_data` module
* Minor improvements to the test and continous integration infrastructure
***************
Acknowledgments
***************

This open source software code was developed in part or in whole in the Human Brain Project, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under Specific Grant Agreements No. 720270, No. 785907 and No. 945539 (Human Brain Project SGA1, SGA2 and SGA3).
.. _contribute:

========================
Contributing to Elephant
========================

You are here to help with Elephant? Awesome, feel welcome to get in touch with
us by asking questions, proposing features and improvements to Elephant.

For guidelines on code documentation, please see the :ref:`documentation_guide`
below.


.. note::

    We highly recommend to get in touch with us *before* starting to implement a
    new feature in Elephant. This way, we can point out synergies with complementary
    efforts and help in designing your implementation such that its integration
    into Elephant will be an easy process.


.. _get_in_touch:

************
Get in touch
************

Using the mailing list
----------------------

General discussion of Elephant development takes place in the
`NeuralEnsemble Google group <http://groups.google.com/group/neuralensemble>`_.

Please raise concrete issues, bugs, and feature requests on the `issue tracker`_.


Using the issue tracker
-----------------------

If you find a bug in Elephant, please create a new ticket on the
`issue tracker`_. Choose one of the available templates - "Bug report",
"Feature request", or "Questions".
Choose an issue title that is as specific as possible to the problem you've found, and
in the description give as much information as you think is necessary to
recreate the problem. The best way to do this is to create the shortest possible
Python script that demonstrates the problem, and attach the file to the ticket.

If you have an idea for an improvement to Elephant, create a ticket with type
"enhancement". If you already have an implementation of the idea, open a
`pull request <https://github.com/NeuralEnsemble/elephant/pulls>`_.

.. _set_up_an_environment:

************************************
Setting up a development environment
************************************

In order to contribute to the Elephant code, you must set up a Python environment on
your local machine.

1. Fork `Elephant <https://github.com/NeuralEnsemble/elephant>`_ as described
   in `Fork a repo <https://help.github.com/en/github/getting-started-with-github/fork-a-repo>`_.
   Download Elephant source code from your forked repo::

    $ git clone git://github.com/<your-github-profile>/elephant.git
    $ cd elephant

2. Set up the virtual environment (either via pip or conda):

.. tabs::

    .. tab:: conda (recommended)

        .. code-block:: sh

            conda env create -f requirements/environment.yml
            conda activate elephant
            pip install -r requirements/requirements-tests.txt

    .. tab:: pip

        .. code-block:: sh

            python -m venv elephant-env
            source elephant-env/bin/activate

            pip install -r requirements/requirements.txt
            pip install -r requirements/requirements-extras.txt  # optional
            pip install -r requirements/requirements-tests.txt


3. Before you make any changes, run the test suite to make sure all the tests
   pass on your system::

    $ nosetests .

   You can specify a particular module to test, for example
   ``test_statistics.py``::

    $ nosetests elephant/test/test_statistics.py

   At the end, if you see "OK", then all the tests passed (or were skipped
   because certain dependencies are not installed), otherwise it will report
   on tests that failed or produced errors.


************************
Adding new functionality
************************

After you've set up a development environment, implement a functionality you
want to add in Elephant. This includes, in particular:

   * fixing a bug;
   * improving the documentation;
   * adding a new functionality.

Writing code
------------

Imagine that you want to add a novel
statistical measure of a list of spike trains that returns an analog signal in the existing module ``elephant/statistics.py``.
Let's call it ``statistics_x``.

Elephant relies on Neo structures that use Quantities extensively, allowing
the users to conveniently specify physical units (e.g., seconds and milliseconds)
of input objects. Typically, to avoid computationally expensive quantities
rescaling operation on large input arrays, we check that the main data objects
- spike trains or analog signals - share the same units and rescale additional
parameters (t_start, t_stop, bin_size, etc.) to the units of input objects.

.. code-block:: python

    import neo
    import quantities as pq

    from elephant.utils import check_same_units


    def statistics_x(spiketrains, t_start=None, t_stop=None):
        """
        Compute the X statistics of spike trains.

        Parameters
        ----------
        spiketrains : list of neo.SpikeTrain
            Input spike trains.
        t_start, t_stop : pq.Quantity or None
            Start and stop times to compute the statistics over the specified
            interval. If None, extracted from the input spike trains.

        Returns
        -------
        signal : neo.AnalogSignal
            The X statistics of input spike trains.
            (More description follows.)

        """
        check_same_units(spiketrains, object_type=neo.SpikeTrain)

        # alternatively, if spiketrains are required to be aligned in time,
        # when t_start and t_stop are not specified, use 'check_neo_consistency'
        # check_neo_consistency(spiketrains, object_type=neo.SpikeTrain, t_start=t_start, t_stop=t_stop)

        # convert everything to spiketrain units and strip off the units
        if t_start is None:
            t_start = spiketrains[0].t_start
        if t_stop is None:
            t_stop = spiketrains[0].t_stop
        units = spiketrains[0].units
        t_start = t_start.rescale(units).item()
        t_stop = t_stop.rescale(units).item()
        spiketrains = [spiketrain.magnitude for spiketrain in spiketrains]

        # do the analysis here on unit-less spike train arrays
        x = ...

        signal = neo.AnalogSignal(x,
                                  units=...,
                                  t_start=t_start,
                                  sampling_rate=...,
                                  name="X statistics of spiketrains",
                                  ...)
        return signal


Testing code
------------

Write at least one test in ``elephant/test/test_module_name.py`` file that
covers the functionality.

For example, to check the correctness of the implemented ``statistics_x``
function, we add unittest code in ``elephant/test/test_statistics.py``,
something like

.. code-block:: python

    import unittest

    import neo
    import quantities as pq
    from numpy.testing import assert_array_almost_equal

    from elephant.statistics import statistics_x


    class StatisticsXTestCase(unittest.TestCase):
        def test_statistics_x_correctness(self):
            spiketrain1 = neo.SpikeTrain([0.3, 4.5, 7.8], t_stop=10, units='s')
            spiketrain2 = neo.SpikeTrain([2.4, 5.6], t_stop=10, units='s')
            result = statistics_x([spiketrain1, spiketrain2])
            self.assertIsInstance(result, neo.AnalogSignal)
            self.assertEqual(result.t_start, 0 * pq.s)
            expected_magnitude = [0, 1, 2]
            assert_array_almost_equal(result.magnitude, expected_magnitude)
            ...  # more checking


Pushing the changes and creating a pull request
-----------------------------------------------

Now you're ready to share the code publicly.

1.  Commit your changes:

    .. code-block:: sh

        git add .
        git commit -m "informative commit message"
        git push

    If this is your first commitment to Elephant, please add your name and
    affiliation/employer in :file:`doc/authors.rst`

2.  Open a `pull request <https://github.com/NeuralEnsemble/elephant/pulls>`_.
    Then we will guide you through the process of merging your code into Elephant.

That's all! We're happy to assist you throughout the process of contributing.

If you experience any problems during one of the steps below, please contact us
and we'll help you.


.. _documentation_guide:

*******************
Documentation Guide
*******************


Writing documentation
---------------------

Each module (python source file) should start with a short description of the
listed functionality. Class and function docstrings should conform to the
`NumPy docstring standard <https://numpydoc.readthedocs.io/en/latest/format.html>`_.

.. note:: Please refer to our :doc:`style_guide`.


Building documentation
----------------------

The documentation in :file:`doc/` folder is written in `reStructuredText
<http://docutils.sourceforge.net/rst.html>`_, using the
`Sphinx <http://sphinx-doc.org/>`_ documentation system. To build the
documentation locally on a Linux system, follow these steps:

1. Install requirements-docs.txt and requirements-tutorials.txt the same way
   it's explained in :ref:`set_up_an_environment` step 3:

   .. code-block:: sh

        pip install -r requirements/requirements-docs.txt
        pip install -r requirements/requirements-tutorials.txt

2. Build the documentation:


   .. code-block:: sh

        cd doc
        export PYTHONPATH=${PYTHONPATH}.:../..
        make html

   ``PYTHONPATH`` environmental variable is set in order to find Elephant
   package while executing jupyter notebooks that are part of the documentation.
   You may also need to install LaTeX support:

   .. code-block:: sh

        sudo apt-get install texlive-full

3. Open :file:`_build/html/index.html` in your browser.

4. (Optional) To check that all URLs in the documentation are correct, run:

   .. code-block:: sh

        make linkcheck


Citations
---------

The citations are in BibTeX format, stored in `doc/bib/elephant.bib
<https://github.com/NeuralEnsemble/elephant/blob/master/doc/bib/elephant.bib>`_.

To cite Elephant, refer to :doc:`citation`.

Each module in ``doc/reference`` folder ends with the reference section:

.. code-block:: rst

    References
    ----------

    .. bibliography:: ../bib/elephant.bib
       :labelprefix: <module name shortcut>
       :keyprefix: <module name>-
       :style: unsrt

where ``<module name>`` is (by convention) the Python source file name, and
``<module name shortcut>`` is what will be displayed to the users.

For example, ``:cite:'spade-Torre2013_132'`` will be rendered as ``sp1`` in
the built HTML documentation, if ``<module name shortcut>`` is set to ``sp``
and ``<module name>`` - to ``spade``.

.. _Issue tracker: https://github.com/NeuralEnsemble/elephant/issues
************************
Authors and contributors
************************

The following people have contributed code and/or ideas to the current version
of Elephant. The institutional affiliations are those at the time of the
contribution, and may not be the current affiliation of a contributor.

* Alper Yegenoglu [1]
* Andrew Davison [2]
* Björn Müller [1]
* Detlef Holstein [2]
* Eilif Muller [3, 4]
* Emiliano Torre [1]
* Espen Hagen [1]
* Jeffrey Gill [11]
* Jan Gosmann [6, 8]
* Julia Sprenger [1]
* Junji Ito [1]
* Michael Denker [1]
* Paul Chorley [1]
* Pierre Yger [2]
* Pietro Quaglio [1]
* Richard Meyes [1]
* Vahid Rostami [1]
* Subhasis Ray [5]
* Robert Pröpper [6]
* Richard C Gerkin [7]
* Bartosz Telenczuk [2]
* Chaitanya Chintaluri [9]
* Michał Czerwiński [9]
* Michael von Papen [1]
* Robin Gutzen [1]
* Felipe Méndez [10]
* Simon Essink [1]
* Alessandra Stella [1]
* Peter Bouss [1]
* Alexander van Meegen [1]
* Aitor Morales-Gregorio [1]
* Cristiano Köhler [1]
* Paulina Dąbrowska [1]
* Jan Lewen [1]
* Alexander Kleinjohann [1]
* Danylo Ulianych [1]
* Anno Kurth [1]
* Regimantas Jurkus [1]
* Philipp Steigerwald [12]
* Manuel Ciba [12]
* Maximilian Kramer [1]
* Florian Porrmann [13]
* Sarah Pilz [13]

1. Institute of Neuroscience and Medicine (INM-6), Computational and Systems Neuroscience & Institute for Advanced Simulation (IAS-6), Theoretical Neuroscience, Jülich Research Centre and JARA, Jülich, Germany
2. Unité de Neurosciences, Information et Complexité, CNRS UPR 3293, Gif-sur-Yvette, France
3. Electronic Visions Group, Kirchhoff-Institute for Physics, University of Heidelberg, Germany
4. Brain-Mind Institute, Ecole Polytechnique Fédérale de Lausanne, Switzerland
5. NIH–NICHD, Laboratory of Cellular and Synaptic Physiology, Bethesda, Maryland 20892, USA
6. Neural Information Processing Group, Institute of Software Engineering and Theoretical Computer Science, Technische Universität Berlin, Germany
7. Arizona State University School of Life Sciences, USA
8. Computational Neuroscience Research Group (CNRG), Waterloo Centre for Theoretical Neuroscience, Waterloo, Canada
9. Nencki Institute of Experimental Biology, Warsaw, Poland
10. Instituto de Neurobiología, Universidad Nacional Autónoma de México, Mexico City, Mexico
11. Case Western Reserve University (CWRU), Cleveland, OH, USA
12. BioMEMS Lab, TH Aschaffenburg University of applied sciences, Germany
13. Cognitronics and Sensor Systems, CITEC, Bielefeld University, Bielefeld, Germany

If we've somehow missed you off the list we're very sorry - please let us know.
***************
Citing Elephant
***************

To refer to the Elephant software package in publications, please use:
**Elephant (RRID:SCR_003833)**.

    
To cite Elephant, please use:

**Denker M, Yegenoglu A, Grün S (2018) Collaborative HPC-enabled workflows on
the HBP Collaboratory using the Elephant framework. Neuroinformatics 2018, P19.
doi: 10.12751/incf.ni2018.0019**

A BibTeX entry for LaTeX users is:

.. code-block:: none

    @conference{elephant18,
        author = {Denker, M. and Yegenoglu, A. and Grün, S.},
        booktitle = {Neuroinformatics 2018},
        title = {{C}ollaborative {HPC}-enabled workflows on the {HBP} {C}ollaboratory using the {E}lephant framework},
        pages = {P19},
        year = {2018},
        doi = {10.12751/incf.ni2018.0019},
        url = {https://abstracts.g-node.org/conference/NI2018/abstracts#/uuid/023bec4e-0c35-4563-81ce-2c6fac282abd},
    }


Further publications directly related to Elephant development 
:cite:`citations-Rostami17_3,citations-Stella19_104022` (see a list of full
`BibTex references <https://github.com/NeuralEnsemble/elephant/blob/master/doc/bib/elephant.bib>`_
used in Elephant documentation).

.. bibliography:: bib/elephant.bib
   :labelprefix: ref
   :keyprefix: citations-
   :style: unsrt
*********
Tutorials
*********

These tutorials provide narrative explanations, sample code, and expected
output for some of the neurophysiological analyses in Elephant. You can browse
the tutorials or launch them in mybinder to change and interact with the code.


Introductory
------------

* Statistics of spike trains.

  Covers ``statistics`` module with an introduction to
  `Neo <https://neo.readthedocs.io/en/stable/>`_ input-output data types like
  `SpikeTrain`, `AnalogSignal`, etc.

  :doc:`View the notebook <../tutorials/statistics>` or run interactively:

  .. image:: https://mybinder.org/badge.svg
     :target: https://mybinder.org/v2/gh/NeuralEnsemble/elephant/master?filepath=doc/tutorials/statistics.ipynb


Advanced
--------

* Unitary Event Analysis.

  The analysis detects coordinated spiking activity that occurs significantly
  more often than predicted by the firing rates of neurons alone. It's superior
  to simple statistics.

  :doc:`View the notebook <../tutorials/unitary_event_analysis>` or run
  interactively:

  .. image:: https://mybinder.org/badge.svg
     :target: https://mybinder.org/v2/gh/NeuralEnsemble/elephant/master?filepath=doc/tutorials/unitary_event_analysis.ipynb

* Gaussian Process Factor Analysis (GPFA).

  GPFA is a dimensionality reduction method for neural trajectory visualization
  of parallel spike trains.

  :doc:`View the notebook <../tutorials/gpfa>` or run interactively:

  .. image:: https://mybinder.org/badge.svg
     :target: https://mybinder.org/v2/gh/NeuralEnsemble/elephant/master?filepath=doc/tutorials/gpfa.ipynb

* Spike Pattern Detection and Evaluation (SPADE)

  :doc:`View the notebook <../tutorials/spade>` or run interactively:

  .. image:: https://mybinder.org/badge.svg
     :target: https://mybinder.org/v2/gh/NeuralEnsemble/elephant/master?filepath=doc/tutorials/spade.ipynb

* Analysis of Sequences of Synchronous EvenTs (ASSET)

  :doc:`View the notebook <../tutorials/asset>` or run interactively:

  .. image:: https://mybinder.org/badge.svg
     :target: https://mybinder.org/v2/gh/NeuralEnsemble/elephant/master?filepath=doc/tutorials/asset.ipynb

* Granger causality

  :doc:`View the notebook <../tutorials/granger_causality>` or run interactively:

  .. image:: https://mybinder.org/badge.svg
     :target: https://mybinder.org/v2/gh/NeuralEnsemble/elephant/master?filepath=doc/tutorials/granger_causality.ipynb


Additional
----------

* Parallel

  ``elephant.parallel`` module provides a simple interface to parallelize
  multiple calls to any user-specified function.

  :doc:`View the notebook <../tutorials/parallel>` or run interactively:

  .. image:: https://mybinder.org/badge.svg
     :target: https://mybinder.org/v2/gh/NeuralEnsemble/elephant/master?filepath=doc/tutorials/parallel.ipynb

..
    Index the notebooks in a hidden toctree to avoid sphinx warnings.

.. toctree::
    :hidden:

    tutorials/asset.ipynb
    tutorials/gpfa.ipynb
    tutorials/parallel.ipynb
    tutorials/statistics.ipynb
    tutorials/unitary_event_analysis.ipynb
    tutorials/granger_causality.ipynb
:orphan:

********************************
Coding Style Guide with Examples
********************************

This guide follows mostly:
https://github.com/numpy/numpy/blob/master/doc/example.py.

In the Python code blocks below, some remarks are included using JavaScript
style notation <!-- comment -->. They provide additional information regarding
the code that is being presented as example.

Module docstring
----------------

The module should start with its own docstring in the first line.
The example below illustrates show it should be implemented.

.. code-block:: python

    """
    Here comes the module description

    Detailed introduction and tutorial of method goes here.

    You can even link tutorials here, if appropriate.

    If you have references you can insert them here, otherwise in the corresponding
    functions, methods or classes.

    For a working example see
    :py:mod:`unitary_event_analysis.py<elephant.unitary_event_analysis>`
    """


Function docstring. Naming conventions
--------------------------------------

The function below illustrates how arguments and functions should be named
throughout Elephant.

.. code-block:: python

    def pair_of_signals_example(spiketrain_i, spiketrain_j):
        # Add '_i' and '_j' suffixes to a pair of signals, spiketrains or any
        # other variables that come in pairs.

    def perfect_naming_of_parameters(spiketrains, spiketrain, reference_spiketrain,
                         target_spiketrain, signal, signals, max_iterations,
                         min_threshold, n_bins, n_surrogates, bin_size, max_size,
                         time_limits, time_range, t_start, t_stop, period, order,
                         error, capacity, source_matrix, cov_matrix,
                         selection_method='aic'
                         ):
        r"""
        Full example of the docstring and naming conventions.

        Function names should be in lowercase, with words written in full, and
        separated by underscores. Exceptions are for common abbreviations, such
        as "psd" or "isi". But words such as "coherence" must be written in full,
        and not truncated (e.g. "cohere").

        If the truncation or abbreviation not in conformity to this naming
        convention was adopted to maintain similarity to a function used
        extensively in another language or package, mention this in the "Notes"
        section, like the comment below:
        <!--
        Notes
        -----
        This function is similar to `welch_cohere` function in MATLAB.
        -->

        The rationale for the naming of each parameter in this example will be
        explained in the relevant "Parameters" section. Class parameters and
        attributes also follow the same naming convention.

        Parameters
        ----------
        <!-- As a general rule, each word is written in full lowercase, separated
        by underscores. Special cases apply according to the examples below -->
        spiketrains : neo.SpikeTrain or list of neo.SpikeTrain
            Within Elephant, this is how to name an input parameter that contains
            at least one spike train. The parameter name is in plural (i.e.,
            `spiketrains`). The function will deal with converting a single
            `neo.SpikeTrain` to a list of `neo.SpikeTrain` if needed.
            Note that although these are two words, they are NOT separated by
            underscore because Neo does not use underscore, and Elephant must keep
            compatibility. Do not use names such as `sts`, `spks`, or
            `spike_trains`.
        spiketrain: neo.SpikeTrain
            If the function EXPLICITLY requires only a single spike train, then
            the parameter should be named in singular (i.e., `spiketrain`). Do
            not use names such as `st`, `spk`, or `spike_train`.
        reference_spiketrain : neo.SpikeTrain
            If a function uses more than one parameter with single spike trains,
            then each parameter name begins with a meaningful name,
            followed by "_spiketrain" in singular form.
        target_spiketrain: neo.SpikeTrain
            Second parameter that is a single spike train. Note that the difference
            from `reference_spiketrain` is indicated by a meaningful name at the
            beginning.
        signal : neo.AnalogSignal
            If a single `neo.AnalogSignal` object is passed to the function, even if
            it contains several signals (arrays).
        signals : list of neo.AnalogSignal
            If the parameter is a container that has at least one `neo.AnalogSignal`
            object. The name of the parameter is `signals` (plural).
        max_iterations : int
            Parameters that represent a maximum value should start with "max_"
            prefix, followed by the description as a full word. Therefore, do not
            use names such as `max_iter` or `maxiter`.
        min_threshold : float
            Same case as for maximum. Parameters that represent a minimum value
            should start with "min_" prefix, followed by the description as a full
            word. Therefore, do not use names such as `min_thr` or `minthr`.
        n_bins : int
            Parameters that represent a number should start with the prefix "n_".
            Do not use `numbins`, `bin_number`, or `num_bins`. The prefix should
            be followed by a meaningful word in full.
        n_surrogates : int
            The description should always be meaningful an without abbreviations.
            Therefore, do not use terms as `n` or `n_surr`, that are not
            immediately understood.
        bin_size : pq.Quantity or int
            Separate the words by underscore. Do not use `bin_size`. Old functions
            which use `binsize` are deprecated.
        max_size : float
            Another example showing that words should be separated by underscores.
            This intersects with the naming convention for a maximum value.
        time_limits: list or tuple
            For parameters that define minimum and maximum values as a list or
            tuple (e.g., [-2, 2]), the parameter must start with a meaningful
            word followed by the suffix "_limits". Preferentially, one should use
            two separated parameters (e.g., `max_time` and `min_time` following
            the convention for maximum and minimum already mentioned). But should
            the function require the definition of limits in this form, use the
            name `_limits` and not `_range` (see next parameter).
        time_range: list
            For parameters that behave like a Python range (e.g. [1, 2, 3, 4])), in
            the sense that it is a sequence, not only the lower and upper limits
            as in the example above, the parameter should start with a meaningful
            name followed by the suffix "_range".
        t_start : pq.Quantity
            Standard name within Elephant for defining starting times.
        t_stop : pq.Quantity
            Standard name within Elephant for defining stopping times.
        period : pq.Quantity
            Oscillation period.
            Always use informative names. In this case, one could name the
            parameter as simply as `T`, since this is standard for referring to
            periods. If the function is implementing computations based on a paper
            that has a formula with a variable "T", acknowledge this after
            describing the formula in the docstring. Therefore, write a sentence
            like "`period` refers to :math:`T`"
            If the Elephant function uses an external function (such as from
            `scipy`), and such function has an argument named `T`, also
            acknowledge this in the docstring. Therefore, write a sentence like
            "`period` is forwarded as argument `T` of `scipy.uses_T` function".
            If the external function already has an informative parameter name
            (such as `period`), the same parameter name can be used in the Elephant
            function if forwarded.
            If several input parameters are forwarded or are members
            of a formula, the docstring can present them together as a list.
            But always use informative names, not single letter names if this is
            how they are described in the paper or implemented in another function.
        order : int
            Order of the Butterworth filter.
            This is an example of how the `N` parameter of `scipy.signal.butter`
            function could be provided by the user of the Elephant function.
            The docstring would present a text similar to
            "`order` is passed as the `N` argument for `scipy.signal.butter` function".
            Also, in the code implementation, use keyword arguments to make this
            explicit (see the implementation of the function below)
        error : float
            In the case the function has an input parameter that corresponds to a
            greek letter in a formula (in a paper, for instance) always use the
            meaning of the greek letter. Therefore, should :math:`\epsilon` refer
            to the error in the formula, the parameter should be named `error`. As
            already mentioned, this is acknowledged in the docstring after the
            description of the formula.
        capacity : float
            Capacity value.
            When using parameters based on a paper (which, e.g., derives some
            formula), and the parameter's name in this paper is a single letter
            (such as `C` for capacity), always use the meaning
            of the letter. Therefore, the parameter should be named `capacity`,
            not `C`. Acknowledge this in the docstring as already mentioned.
        source_matrix: np.ndarray
            Parameters that are matrices should end with the suffix "_matrix", and
            start with a meaningful name.
        cov_matrix: np.ndarray
            A few exceptions allow the use of abbreviations instead of full words
            in the name of the parameter. These are:
            * "cov" for "covariance" (e.g., `cov_matrix`)
            * "lfp" for "local_field_potential" (e.g. `lfp_signal`)
            * "corr" for "correlation" (e.g. `corr_matrix`).
            THESE EXCEPTIONS ARE NOT ACCEPTED FOR FUNCTION NAMES. Therefore, a
            parameter would be named `cov_matrix`, but the function would be named
            `calculate_covariance_matrix`. If the function name becomes very long,
            then an alias may be created and described appropriately in the "Notes"
            section, as mentioned above. For aliases, see example below.
        selection_method : {'aic', 'bic'}
            Metric for selecting the autoregressive model.
            If 'aic', uses the Akaike Information Criterion (AIC).
            If 'bic', uses Bayesian Information Criterion (BIC).
            Default: 'bic', because it is more reliable than AIC due to the
            mathematical properties (see Notes [3]).
            <!-- Note that the default value that comes in the last line is
            followed by comma and a brief reasoning for defining the default
            `selection_method`). -->

        <!-- Other remarks:
        1. Do not use general parameter names, such as `data` or `matrix`.
        2. Do not use general result names, such as `result` or `output`.
        3. Avoid capitalization (such as the examples mentioned for parameters
           such as `T` for period, or `C` for capacity or a correlation matrix.
        -->

        Returns
        -------
            frequency : float
                The frequency of the signal.
            filtered_signal : np.ndarray
                Signal filtered using Butterworth filter.

        Notes
        -----
        1. Frequency is defined as:

        .. math::

            f = \frac{1}{T}

           `period` corresponds to :math:`T`

        2. `order` is passed as the `N` parameter when calling
           `scipy.signal.butter`.
        3. According to [1]_, BIC should be used instead of AIC for this
           computation. The brief rationale is .......

        References
        ----------
        .. [1] Author, "Why BIC is better than AIC for AR model", Statistics,
               vol. 1, pp. 1-15, 1996.

        """
        # We use Butterworth filter from scipy to perform some calculation.
        # Note that parameter `N` is passed using keys, taking the value of the
        # `order` input parameter
        filtered_signal = scipy.signal.butter(N=order, ...)

        # Here we calculate a return value using a function variable. Note that
        # this variable is named in the "Returns" section
        frequency = 1 / period
        return frequency, filtered_signal



Class docstring
---------------

Class docstrings follow function docstring format. Here is an example.

.. code-block:: python

    class MyClass(object):  # Classes use CamelCase notation
        """
        One line description of class.

        Long description of class, may span several lines. Possible sections are
        the same as for a function doc, with additional "Attributes" and "Methods"
        after "Parameters" (cf. numpydoc guide). Do not put a blank line after
        section headers, do put a blank line at the end of a long docstring.

        When explaining the algorithm, you can use mathematical notation, e.g.:

        .. math::

            E = m c^2

        To insert an equation use `.. math::` and surround the whole expression in
        blank lines. To use math notation in-line, write :math:`E` corresponds to
        energy and :math:`m` to mass. Embed expressions after `:math:` in
        backticks, e.g. :math:`x^2 + y^2 = z^2`.

        To refer to a paper in which formula is described, use the expression
        "see [1]_" - it will become an interactive link on the readthedocs website.
        The underscore after the closing bracket is mandatory for the link to
        work.

        To refer to a note in the "Notes" section, simply write "see Notes [1]".

        Variable, module, function, and class names should be written
        between single back-ticks (`kernels.AlphaKernel`), NOT *bold*.

        For common modules such as Numpy and Quantities, use the notation
        according to the import statement. For example:
        "this function uses `np.diff`", not "uses `numpy.diff`".

        Prefixes for common packages in Elephant are the following:

        1. Neo objects = neo (e.g. `neo.SpikeTrain`)
        2. Numpy = np (e.g. `np.ndarray`)
        3. Quantities = pq (e.g. `pq.Quantity`)

        For other objects, list the full path to the object (e.g., for the
        BinnedSpikeTrain, this would be `elephant.conversion.BinnedSpikeTrain`)

        For None and NaNs, do not use backticks. NaN is referred as np.NaN (i.e.,
        with the Numpy prefix "np").

        Use backticks also when referring to arguments of a function (e.g., `x` or
        `y`), and :attr:`attribute_name` when referring to attributes of a class
        object in docstrings of this class.

        To refer to attributes of other objects, write
        `other_object.relevant_attribute` (e.g. `neo.SpikeTrain.t_stop`).

        When mentioning a function from other module, type `other_module.function`
        (without parentheses after the function name; e.g., `scipy.signal.butter`).

        If you refer values to True/False/None, do not use backticks, unless an
        emphasis is needed. In this case, write `True` and not bold, like **True**.

        Parameters
        ----------
        <!-- List the arguments of the constructor (__init__) here!
        Arguments must come in the same order as in the constructor or function -->
        parameter : int or float
            Description of parameter `parameter`. Enclose variables in single
            backticks. The colon must be preceded by a space.
        no_type_parameter
            Colon omitted if the type is absent.
        x : float
            The X coordinate.
        y : float
            The Y coordinate.
            Default: 1.0.  <!-- not "Default is 1.0." (it is just a convention) -->
        z : float or int or pq.Quantity
            This is Z coordinate.
            If it can take multiple types, separate them by "or", do not use commas
            (numpy style).
            If different actions will happen depending on the type of `z`, explain
            it briefly here, not in the main text of the function/class docstring.
        s : {'valid', 'full', 'other'}
            This is the way to describe a list of possible argument values, if the
            list is discrete and predefined (typically concerns strings).
            If 'valid', the function performs some action.
            If 'full', the function performs another action.
            If 'other', the function will ignore the value defined in `z`.
            Default: 'valid'.
        spiketrains : neo.SpikeTrain or list of neo.SpikeTrain or np.ndarray
            When the parameter can be a container (such as list or tuple), you can
            specify the type of elements using "of". But use the Python type name
            (do not add "s" to make it plural; e.g., do not write
            "list of neo.SpikeTrains" or "list of neo.SpikeTrain objects").
        counts_matrix : (N, M) np.ndarray
            This is the way to indicate dimensionality of the required array
            (i.e.,if the function only works with 2D-arrays). `N` corresponds to
            the number of rows and `M` to the number of columns. Refer to the same
            `N` and `M` to describe the dimensions of the returned values when
            they are determined by the dimensions of the parameter.
        is_true : bool
            True, if 1.
            False, if 0.
            Default: True.
        other_parameter : int
            Some value.
            If value is None and the function takes some specific action (e.g.,
            calculate some value based on the other inputs), describe here.
            Default: None.

        Attributes
        ----------
        <!-- Here list the attributes of class object which are not simply copies
        of the constructor parameters. Property decorators (@property) are also
        considered attributes -->
        a : list
            This is calculated based on `x` and `y`.
        b : int
            This is calculated on the way, during some operations.

        Methods
        -------
        <!--  Here list the most important/useful class methods (not all the
        methods) -->

        Returns
        -------
        <!-- This section is rarely used in class docstrings, but often in
        function docs. Follow the general recommendation of numpydoc.
        If there is more than one returned value, use variable names for the
        returned value, like `error_matrix` below. -->
        error_matrix : np.ndarray
            A matrix is stored in a variable called `error_matrix`, containing
            errors estimated from some calculations. The function "return"
            statement then returns the variable (e.g. "return error_matrix").
            Format is the same as for any parameter in section "Parameters".
            Use meaningful names, not general names such as `output` or `result`.
        list
            The returned object is created on the fly and is never assigned to
            a variable (e.g. "return [1, 2, 3]"). Simply name the type and
            describe the content. This should be used only if the function returns
            a single value.
        dict
            key_1 : type
                Description of key_1, formatted the same as in "Parameters".
            key_2 : type
                Description of key_2
        particular_matrix : (N, N, M) np.ndarray
            The dimensionality of this array depends on the dimensionality of
            `counts_matrix` input parameter. Note that `N` and `M` are used since
            these were the names of the dimensions of `counts_matrix` in the
            "Parameters" section.
        list_variable : list of np.ndarray
            Returns a list of numpy arrays.
        signal : int
            Description of `signal`.

        Raises
        ------
        <!-- List the errors explicitly raised by the constructor (raise
        statements), even if they are in fact raised by other Elephant functions
        called inside the constructor. Enumerate them in alphabetical order. -->
        TypeError
            If `x` is an `int` or None.
            If `y` is not a `float`.
        ValueError
            If this and that happens.

        Warns
        -----
        <!-- Here apply the same rules as for "Raises". -->
        UserWarning
            If something may be wrong but does not prevent execution of the code.
            The default warning type is UserWarning.

        Warning
        -------
        <!-- Here write a message to the users to warn them about something
        important.
        Do not enumerate Warnings in this section! -->

        See Also
        --------
        <!-- Here refer to relevant functions (also from other modules). Follow
        numpydoc recommendations.
        If the function name is not self-explanatory, you can add a brief
        explanation using a colon separated by space.
        This items will be placed as links to the documentation of the function
        referred.
        -->
        statistics.isi
        scipy.signal.butter : Butterworth filter

        Notes
        -----
        <!-- Here you can add some additional explanations etc. If you have several
        short notes (at least two), use a list -->
        1. First remark.
        2. Second much longer remark, which will span several lines. To refer to a
           note in other parts of the docstring, use a phrase like "See Notes [2]".
           To make sure that the list displays correctly, keep the indentation to
           match the first word after the point (as in this text).
        3. If you want to explain why the default value of an argument is
           something particular, you can give a more elaborate explanation here.
        4. If the function has an alias (see the last function in this file), the
           information about it should be in this section in the form:
           Alias: bla.
           Aliases should be avoided.
        5. Information about validation should be here, and insert bibliographic
           citation in the "References". Also specify in parenthesis the unit test
           that implements the validation. Example:
           "This function reproduces the paper Riehle et al., 1997 [2]_.
           (`UETestCase.test_Riehle_et_al_97_UE`)."
        6. Do not create new section names, because they will not be displayed.
           Place the relevant information here instead.
        7. This is an optional section that provides additional information about
           the code, possibly including a discussion of the algorithm. This
           section may include mathematical equations, written in LaTeX format.
           Inline: :math:`x^2`. An equation:

           .. math::

           x(n) * y(n) \Leftrightarrow X(e^{j\omega } )Y(e^{j\omega } )

        8. Python may complain about backslashes in math notation in docstrings.
           To prevent the complains, precede the whole docstring with "r" (raw
           string).
        9. Images are allowed, but should not be central to the explanation;
           users viewing the docstring as text must be able to comprehend its
           meaning without resorting to an image viewer. These additional
           illustrations are included using:

            .. image:: filename

        References
        ----------
        .. [1] Smith J., "Very catchy title," Elephant 1.0.0, 2020. The ".." in
               front makes the ref referencable in other parts of the docstring.
               The indentation should match the level of the first word AFTER the
               number (in this case "Smith").

        Examples
        --------
        <!-- If applicable, provide some brief description of the example, then
        leave a blank line.
        If the second example uses an import that was already used in the first
        example, do not write the import again.
        Examples should be very brief, and should avoid plotting. If plotting
        is really needed, use simple matplotlib plots, that take only few lines.
        More complex examples, that require lots of plotting routines (e.g.,
        similar to Jupyter notebooks), should be placed as tutorials, with links
        in the docstring. Examples should not load any data, but only use easy
        generated data.
        Finally, avoid using abbreviations in examples, such as
        "import elephant.conversion as conv" -->

        >>> import neo
        >>> import numpy as np
        >>> import quantities as pq
        ...
        ... # This is a way to make a blank line within the example code.
        >>> st = neo.SpikeTrain([0, 1, 2, 3] * pq.ms, t_start=0 * pq.ms,
        ...                     t_stop=10 * pq.ms, sampling_rate=1 * pq.Hz)
        ... # Use "..." also as a continuation line.
        >>> print(st)
        SpikeTrain

        Here provide a brief description of a second example. Separate examples
        with a blank line even if you do not add any description.

        >>> import what_you_need
        ...
        >>> st2 = neo.SpikeTrain([5, 6, 7, 8] * pq.ms, t_start=0 * pq.ms,
        ...                      t_stop=10 * pq.ms, sampling_rate=1 * pq.Hz)
        >>> sth = what_you_need.function(st2)
        >>> sth_else = what_you_need.interesting_function(sth)

        """

        def __init__(self, parameter):
            """
            Constructor
            (actual documentation is in class documentation, see above!)
            """
            self.parameter = parameter
            self.function_a()  # creates new attribute of self 'a'

        def function_a(self, parameter, no_type_parameter, spiketrains,
                       is_true=True, string_parameter='C', other_parameter=None):
            """
            One-line short description of the function.

            Long description of the function. Details of what the function is doing
            and how it is doing it. Used to clarify functionality, not to discuss
            implementation detail or background theory, which should rather be
            explored in the "Notes" section below. You may refer to the parameters
            and the function name, but detailed parameter descriptions still
            belong in the "Parameters" section.

            Parameters
            ----------
            <!-- See class docstring above -->

            Returns
            -------
            <!-- See class docstring above -->

            Raises
            ------
            <!-- See class docstring above.
            List only exceptions explicitly raised by the function -->

            Warns
            -----
            <!-- See class docstring above. -->

            See Also
            --------
            <!-- See class docstring above  -->

            Notes
            -----
            <!-- See class docstring above -->

            References
            ----------
            <!-- See class docstring above -->

            Examples
            --------
            <!-- See class docstring above -->

            """

            # Variables use underscore notation
            dummy_variable = 1
            a = 56  # This mini comment uses two spaces after the code!

            # Textual strings use double quotes
            error = "An error occurred. Please fix it!"
            # Textual strings are usually meant to be printed, returned etc.

            # Non-textual strings use single quotes
            default_character = 'a'
            # Non textual strings are single characters, dictionary keys and other
            # strings not meant to be returned or printed.

            # Normal comments are proceeded by a single space, and begin with a
            # capital letter
            dummy_variable += 1

            # Longer comments can have several sentences. These should end with a
            # period. Just as in this example.
            dummy_variable += 1

        # Class functions need only 1 blank line.
        # This function is deprecated. Add a warning!
        def function_b(self, **kwargs):
            """
            This is a function that does b.

            .. deprecated:: 0.4
              `function_b` will be removed in elephant 1.0, it is replaced by
              `function_c` because the latter works also with Numpy Ver. 1.6.

            Parameters
            ----------
            kwargs : dict
                kwarg1 : type
                    Same style as docstring of class `MyClass`.
                kwarg2 : type
                    Same style as docstring of class `MyClass`.

            """
            pass
============================
Function Reference by Module
============================


****************************************************
Local field potentials (LFPs) and population signals
****************************************************


.. toctree::
    :maxdepth: 1

    reference/signal_processing
    reference/spectral
    reference/causality
    reference/current_source_density


************
Spike trains
************

.. toctree::
    :maxdepth: 2

    reference/statistics

.. toctree::
    :maxdepth: 2

    reference/_spike_train_processing
    reference/_spike_train_patterns

.. toctree::
    :maxdepth: 1

    reference/change_point_detection
    reference/gpfa
.. toctree::
    :maxdepth: 1

    reference/spike_train_surrogates

.. toctree::
    :maxdepth: 2

    reference/spike_train_generation
    

********************************
LFPs and spike trains (combined)
********************************

.. toctree::
    :maxdepth: 1

    reference/sta
    reference/phase_analysis

*******
Kernels
*******

.. toctree::
    :maxdepth: 2

    reference/kernels


*********
Waveforms
*********

.. toctree::
    :maxdepth: 1

    reference/waveform_features

********************************
Alternative data representations
********************************

.. toctree::
    :maxdepth: 1

    reference/conversion

*************
Miscellaneous
*************

.. toctree::
    :maxdepth: 1

    reference/neo_tools
    reference/utils
    reference/pandas_bridge
    reference/parallel
=================
Maintainers guide
=================

This guide is for Elephant maintainers only.


Python 3
--------

Backward compatibility is achieved by putting a few future imports at the
beginning of each source file:

.. code-block:: python

    from __future__ import division, print_function, unicode_literals

All code should conform as much as possible to
`PEP 8 <http://www.python.org/dev/peps/pep-0008/>`_.

Each source file should have a copyright and a license note.


Structure of the doc/ folder
----------------------------

Documentation in Elephant is written exclusively using the ``sphinx`` package
and resides in the ``doc`` folder, in addition to the docstrings contained of
the ``elephant`` modules. In the following, we outline the main components of
the Elephant documentation.


Top-level documentation
~~~~~~~~~~~~~~~~~~~~~~~

General information about the Elephant package and a gentle introduction are
contained in various ``.rst`` files in the top-level directory of the Elephant
package. Here, :file:`index.rst` is the central starting point, and the hierarchical
document structure is specified using the ``toctree`` directives. In particular,
these files contain a general introduction and tutorial on Elephant, the
release notes of Elephant versions, and this development guide.


Module and function reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All modules in Elephant are semi-automatically documented. To this end, for
each module ``x`` a file ``doc/reference/x.rst`` exists with the following
contents:

.. code:: rst

    ============================
    `x` - Short descriptive name
    ============================

    .. automodule:: elephant.x

This instructs Sphinx to add the module documentation in the module docstring
into the file.

The module docstring of ``elephant/x.py`` is also standardized in its structure:

.. code:: rst

    .. include:: x-overview.rst

    .. current_module elephant.x

    Overview of Functions
    ---------------------

    <<Heading 1 to group functions (optional)>>
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    .. autosummary::
        :toctree: toctree/x/

        function1
        function2

    <<Heading 2 to group functions (optional)>>
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    .. autosummary::
        :toctree: toctree/x/

        function3
        function4


Each module documentation starts with a short, understandable introduction to
the functionality of the module, the "Overview". This text is written in a
separate file residing in `doc/reference/x-overview.rst`, and is included on
the first line. This keeps the docstring in the code short.

Next, we specify the current module as `x`, in order to avoid confusion if
a module uses submodules.

In the following, all functions of the module are listed in an order that is
intuitive for users. Where it makes sense, these functions can be thematically
grouped using third-level headings. For small modules, no such headings are
needed.



Making a release
----------------

1. Increment the Elephant package version in :file:`elephant/VERSION`.

2. Add a section in :file:`doc/release_notes.rst`, describing in short the
   changes made from the previous release.

3. Check that the copyright statement (in :file:`LICENSE.txt`,
   :file:`README.md`, and :file:`doc/conf.py`) is correct.

4. If there is a new module do not forget to add the module name to the
   :file:`doc/modules.rst` and make a file with a short description in
   :file:`doc/reference/<modulename>.rst`.

5. Push the commit with release notes and version updated to github.

6. Remove :file:`elephant/spade_src/fim.so`. Otherwise, it'll be included in
   the built package (it should be downloaded at pip install).

7. Build a source package and upload it to PyPi.

   Build a source package (see `Packaging Python Projects
   <https://packaging.python.org/tutorials/packaging-projects/#generating-distribution-archives>`_)::

    $ pip install --user --upgrade twine
    $ python setup.py sdist

   To upload the package to `PyPI <http://pypi.python.org>`_
   (if you have the necessary permissions)::

    $ python -m twine upload dist/elephant-X.Y.Z.tar.gz

8. Finally, make a release on GitHub UI page and copy-paste the release notes.
   Then tag the release in the Git repository and push it::

    $ git tag <version>
    $ git push --tags upstream

   Here, version should be of the form ``vX.Y.Z``.
=============================================
Elephant - Electrophysiology Analysis Toolkit
=============================================

*Elephant* (Electrophysiology Analysis Toolkit) is an emerging open-source,
community centered library for the analysis of electrophysiological data in
the Python programming language.

The focus of Elephant is on generic analysis functions for spike train data and
time series recordings from electrodes, such as the local field potentials
(LFP) or intracellular voltages. In addition to providing a common platform for
analysis codes from different laboratories, the Elephant project aims to
provide a consistent and homogeneous analysis framework that is built on a
modular foundation. Elephant is the direct successor to Neurotools_ and
maintains ties to complementary projects such as ephyviewer_ and
neurotic_ for raw data visualization.

The input-output data format is either Neo_, Quantity_ or Numpy_ array.
Quantity is a Numpy-wrapper package for handling physical quantities like
seconds, milliseconds, Hz, volts, etc. Quantity is used in both Neo and
Elephant.


**Visualization of Elephant analysis objects**

`Viziphant <https://viziphant.readthedocs.io/en/latest/>`_ package is developed
by Elephant team and provides a high-level API to easily generate plots and
interactive visualizations of neuroscientific data and analysis results.
The API uses and extends the same structure as in Elephant to ensure intuitive
usage for scientists that are used to Elephant.


*****************
Table of Contents
*****************

* :doc:`install`
* :doc:`tutorials`
* :doc:`modules`
* :doc:`contribute`
* :doc:`release_notes`
* :doc:`acknowledgments`
* :doc:`authors`
* :doc:`citation`


.. toctree::
    :maxdepth: 2
    :hidden:

    install
    tutorials
    modules
    contribute
    release_notes
    acknowledgments
    authors
    citation


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`


.. _Neurotools:  http://neuralensemble.org/NeuroTools/
.. _ephyviewer:  https://ephyviewer.readthedocs.io/en/latest/
.. _neurotic:  https://neurotic.readthedocs.io/en/latest/
.. _Neo: http://neuralensemble.org/neo/
.. _Numpy: http://www.numpy.org/
.. _Quantity: https://python-quantities.readthedocs.io/en/latest/


.. |date| date::
.. |time| date:: %H:%M
.. _install:

============
Installation
============

The easiest way to install Elephant is by creating a conda environment, followed by ``pip install elephant``.
Below is the explanation of how to proceed with these two steps.


.. _prerequisites:

*************
Prerequisites
*************

Elephant requires `Python <http://python.org/>`_ 3.6, 3.7, 3.8, or 3.9.

.. tabs::


    .. tab:: (recommended) Conda (Linux/MacOS/Windows)

        1. Create your conda environment (e.g., `elephant`):

           .. code-block:: sh

              conda create --name elephant python=3.7 numpy scipy tqdm

        2. Activate your environment:

           .. code-block:: sh

              conda activate elephant


    .. tab:: Debian/Ubuntu

        Open a terminal and run:

        .. code-block:: sh

           sudo apt-get install python-pip python-numpy python-scipy python-pip python-six python-tqdm


************
Installation
************

.. tabs::


    .. tab:: Stable release version

        The easiest way to install Elephant is via `pip <http://pypi.python.org/pypi/pip>`_:

           .. code-block:: sh

              pip install elephant

        If you want to use advanced features of Elephant, install the package
        with extras:

           .. code-block:: sh

              pip install elephant[extras]


        To upgrade to a newer release use the ``--upgrade`` flag:

           .. code-block:: sh

              pip install --upgrade elephant

        If you do not have permission to install software systemwide, you can
        install into your user directory using the ``--user`` flag:

           .. code-block:: sh

              pip install --user elephant


    .. tab:: Development version

        If you have `Git <https://git-scm.com/>`_ installed on your system,
        it is also possible to install the development version of Elephant.

        1. Before installing the development version, you may need to uninstall
           the previously installed version of Elephant:

           .. code-block:: sh

              pip uninstall elephant

        2. Clone the repository and install the local version:

           .. code-block:: sh

              git clone git://github.com/NeuralEnsemble/elephant.git
              cd elephant

        .. tabs::

            .. tab:: Minimal setup

                .. code-block:: sh

                    pip install -e .


            .. tab:: conda (with extras)

                .. code-block:: sh

                    conda remove -n elephant --all  # remove the previous environment
                    conda env create -f requirements/environment.yml
                    conda activate elephant
                    pip install -e .

***********
MPI support
***********

Some Elephant modules (ASSET, SPADE, etc.) are parallelized to run with MPI.
In order to make use of MPI parallelization, you need to install ``mpi4py``
package:

.. tabs::

    .. tab:: conda (easiest)

        .. code-block:: sh

            conda install -c conda-forge mpi4py

    .. tab:: pip (Linux)

        .. code-block:: sh

            sudo apt install -y libopenmpi-dev openmpi-bin
            pip install mpi4py

To run a python script that supports MPI parallelization, run in a terminal:

.. code-block:: sh

    mpiexec -n numprocs python -m mpi4py pyfile [arg] ...

For more information, refer to `mpi4py
<https://mpi4py.readthedocs.io/en/stable/mpi4py.run.html>`_ documentation.


***********************
CUDA and OpenCL support
***********************

:ref:`asset` module supports CUDA and OpenCL. These are experimental features.
You can have one, both, or none installed in your system.

.. tabs::

    .. tab:: CUDA

        To leverage CUDA acceleration on an NVIDIA GPU card, `CUDA toolkit
        <https://developer.nvidia.com/cuda-downloads>`_ must installed on
        your system. Then run the following command in a terminal:

        .. code-block:: sh

            pip install pycuda

        In case you experience issues installing PyCUDA, `this guide
        <https://medium.com/leadkaro/setting-up-pycuda-on-ubuntu-18-04-for-
        gpu-programming-with-python-830e03fc4b81>`_ offers a step-by-step
        installation manual.

        If PyCUDA is detected and installed, CUDA backend is used by default in
        Elephant ASSET module. To turn off CUDA support, set ``ELEPHANT_USE_CUDA``
        environment flag to ``0``.


    .. tab:: OpenCL

        If you have a laptop with a built-in Intel Graphics Card, you can still
        leverage significant performance optimization with OpenCL backend.
        The simplest way to install PyOpenCL is to run a conda command:

        .. code-block:: sh

            conda install -c conda-forge pyopencl intel-compute-runtime

        However, if you have root (sudo) privileges, it's recommended to install
        up-to-date `Intel Graphics Compute Runtime
        <https://github.com/intel/compute-runtime/releases>`_ system-wide and then
        install PyOpenCL as follows:

        .. code-block:: sh

            conda install -c conda-forge pyopencl ocl-icd-system

        Set ``ELEPHANT_USE_OPENCL`` environment flag to ``0`` to turn off
        PyOpenCL support.

        .. note::

            Make sure you've disabled GPU Hangcheck as described in the
            `Intel GPU developers documentation <https://software.intel.com/
            content/www/us/en/develop/documentation/get-started-with-intel-
            oneapi-base-linux/top/before-you-begin.html>`_. Do it with caution -
            using your graphics card to perform computations may make the system
            unresponsive until the compute program terminates.


************
Dependencies
************

Elephant relies on two special packages, installed by default:

    * `quantities <http://pypi.python.org/pypi/quantities>`_ - support for physical quantities with units (mV, ms, etc.)
    * `neo <http://pypi.python.org/pypi/neo>`_ - electrophysiology data manipulations
{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :members:
   :inherited-members:
   :special-members: __call__, __getitem__, __setitem__

   {% block methods %}

   {% if methods %}
   .. rubric:: {{ _('Methods') }}

   .. autosummary::
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Attributes') }}

   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}====================================
Correlative measures on spike trains
====================================

***********************
Spike train correlation
***********************

.. automodule:: elephant.spike_train_correlation


*************************
Spike train dissimilarity
*************************

.. automodule:: elephant.spike_train_dissimilarity


*********************
Spike train synchrony
*********************

.. automodule:: elephant.spike_train_synchrony


References
----------

.. bibliography:: ../bib/elephant.bib
   :labelprefix: cor
   :keyprefix: correlation-
   :style: unsrt

.. bibliography:: ../bib/elephant.bib
   :labelprefix: ds
   :keyprefix: dissimilarity-
   :style: unsrt

.. bibliography:: ../bib/elephant.bib
   :labelprefix: syn
   :keyprefix: synchrony-
   :style: unsrt
======================
Spike train generation
======================


.. automodule:: elephant.spike_train_generation
======================
Spike train surrogates
======================


.. automodule:: elephant.spike_train_surrogates


References
----------

.. bibliography:: ../bib/elephant.bib
   :labelprefix: sr
   :keyprefix: surrogates-
   :style: unsrt
=================
Signal processing
=================


.. automodule:: elephant.signal_processing


References
----------

.. bibliography:: ../bib/elephant.bib
   :labelprefix: sig
   :keyprefix: signal-
   :style: unsrt
.. _asset:

===================================================
Analysis of Sequences of Synchronous EvenTs (ASSET)
===================================================

.. automodule:: elephant.asset.asset


References
----------

.. bibliography:: ../bib/elephant.bib
   :labelprefix: asset
   :keyprefix: asset-
   :style: unsrt
=====================
Neo objects utilities
=====================

.. automodule:: elephant.neo_tools
=====================================
Detection of non-stationary processes
=====================================

.. automodule:: elephant.change_point_detection


References
----------

.. bibliography:: ../bib/elephant.bib
   :labelprefix: cpd
   :keyprefix: cpd-
   :style: unsrt
==================
Causality measures
==================

.. automodule:: elephant.causality.granger


References
----------

.. bibliography:: ../bib/elephant.bib
   :labelprefix: gr
   :keyprefix: granger-
   :style: unsrt
=================
Utility functions
=================

.. automodule:: elephant.utils
=================
Spectral analysis
=================

.. automodule:: elephant.spectral
=============================
BinnedSpikeTrain (conversion)
=============================

.. automodule:: elephant.conversion
=======
Kernels
=======

.. automodule:: elephant.kernels
==============================================
Spike Pattern Detection and Evaluation (SPADE)
==============================================

.. automodule:: elephant.spade


References
----------

.. bibliography:: ../bib/elephant.bib
   :labelprefix: spade
   :keyprefix: spade-
   :style: unsrt
=======================================
Gaussian-Process Factor Analysis (GPFA)
=======================================

.. automodule:: elephant.gpfa.gpfa


References
----------

.. bibliography:: ../bib/elephant.bib
   :labelprefix: gpfa
   :keyprefix: gpfa-
   :style: unsrt
==========================
Statistics of spike trains
==========================

.. automodule:: elephant.statistics
=================
Waveform features
=================

.. automodule:: elephant.waveform_features


References
----------

.. bibliography:: ../bib/elephant.bib
   :labelprefix: wf
   :keyprefix: waveforms-
   :style: unsrt
=========================
Spike-triggered LFP phase
=========================

.. automodule:: elephant.phase_analysis
============================================================
Cumulant Based Inference of higher-order Correlation (CuBIC)
============================================================

.. automodule:: elephant.cubic


References
----------

.. bibliography:: ../bib/elephant.bib
   :labelprefix: cubic
   :keyprefix: cubic-
   :style: unsrt
=============================
Cell assembly detection (CAD)
=============================

.. automodule:: elephant.cell_assembly_detection


References
----------

.. bibliography:: ../bib/elephant.bib
   :labelprefix: cad
   :keyprefix: cad-
   :style: unsrt
=======================
Spike-triggered average
=======================

.. automodule:: elephant.sta
===========================
Unitary Event Analysis (UE)
===========================

.. automodule:: elephant.unitary_event_analysis

Author Contributions
--------------------

- Vahid Rostami (VR)
- Sonja Gruen (SG)
- Markus Diesmann (MD)

VR implemented the method, SG and MD provided guidance.


References
----------

.. bibliography:: ../bib/elephant.bib
   :labelprefix: ue
   :keyprefix: unitary_event_analysis-
   :style: unsrt
========
Parallel
========

.. automodule:: elephant.parallel
============================
Bridge to the pandas library
============================

.. automodule:: elephant.pandas_bridge
=======================================
Detection of synchronous spike patterns
=======================================

.. toctree::
    :maxdepth: 1

    cell_assembly_detection
    unitary_event_analysis
    asset
    spade
    cubic
===============================
Current source density analysis
===============================

.. automodule:: elephant.current_source_density

Keywords
--------
LFP, CSD, multielectrode, Laminar electrode, Barrel cortex.


Citation Policy
---------------
See `current_source_density_src/README.md
<https://github.com/NeuralEnsemble/elephant/blob/master/elephant/
current_source_density_src/README.md>`_


Author Contributions
--------------------
- Chaitanya Chintaluri (CC)
- Espen Hagen (EH)
- Michał Czerwinski (MC)

EH implemented the iCSD methods and StandardCSD.
CC implemented the kCSD methods, kCSD1D (MC and CC).
CC and EH developed the interface to Elephant.
