# EMG Signal Processing Pipeline

**pyemgpipeline** is an electromyography (EMG) signal processing pipeline package.

This package implements internationally accepted EMG processing conventions
and provides a high-level interface for ensuring user adherence to those conventions,
in terms of (1) processing parameter values, (2) processing steps, and (3) processing step order.

The processing steps included in the package are
DC offset removal, bandpass filtering, full wave rectification, linear envelope,
end frame cutting, amplitude normalization, and segmentation.

## Scope

This package defines the processing pipeline for both surface EMG and
intramuscular EMG but not for high density EMG.
The EMG recording requires that the minimum sample rate be at least twice the
highest cutoff frequency of the bandpass filter based on the Nyquist theorem.

## Overview

In **pyemgpipeline**, class `DataProcessingManager` in module `wrappers` is
designed as the main wrapper for high-level, guided processing,
and users are encouraged to use it to adhere to accepted EMG processing conventions.
The other classes, methods, and functions are considered as lower level processing
options.

The package is organized in modules `processors`, `wrappers`, and `plots`.

Module `processors` includes the base class `BaseProcessor` of all signal
processors and seven classes for different processing steps:
`DCOffsetRemover`, `BandpassFilter`, `FullWaveRectifier`, `LinearEnvelope`,
`EndFrameCutter`, `AmplitudeNormalizer`, and `Segmenter`.

Module `wrappers` includes three wrapper classes to facilitate the signal
processing by integrating data and individual processors.
Class `EMGMeasurement` works for data of a single trial,
class `EMGMeasurementCollection` works for data of multiple trials,
and class `DataProcessingManager` is the high-level, guided processing wrapper
with EMG processing conventions.

Module `plots` includes
the function `plot_emg` to plot EMG signals on `matplotlib` figures
and the class `EMGPlotParams` to manage the plot-related parameters.

## Documentation

The [documentation](https://aalhossary.github.io/pyemgpipeline/)
describes how to use this package, including
package installation, quick start, examples explaining the breadth of the package’s functionality,
and API reference.

## Community Guidelines

For contribution, please clone the repository, make changes, and create a pull request.

For reporting any issues, please use
[github issues](https://github.com/aalhossary/pyemgpipeline/issues).

For support, please contact the authors via their emails
or [github issues](https://github.com/aalhossary/pyemgpipeline/issues).

## Citation

If you use this package in your project, please cite this work.


---
title: 'pyemgpipeline: A Python package for electromyography processing'
tags:
  - Python
  - electromyography
  - EMG processing
authors:
  - name: Tsung-Lin Wu
    orcid: 0000-0002-1823-6818
    affiliation: 1
  - name: Amr A. Alhossary
    orcid: 0000-0002-4470-5817
    affiliation: 2
  - name: Todd C. Pataky
    orcid: 0000-0002-8292-7189
    affiliation: 3
  - name: Wei Tech Ang
    orcid: 0000-0002-5778-7719
    affiliation: "1, 2"
  - name: Cyril J. Donnelly^[corresponding author]
    orcid: 0000-0003-2443-5212
    affiliation: 2
affiliations:
  - name: Nanyang Technological University, School of Mechanical & Aerospace Engineering
    index: 1
  - name: Nanyang Technological University, Rehabilitation Research Institute of Singapore
    index: 2
  - name: Kyoto University, Department of Human Health Sciences
    index: 3
date: 17 January 2022
bibliography: paper.bib
---

# Summary

We have developed an electromyography (EMG) signal processing pipeline package called `pyemgpipeline`, which is suitable for both surface EMG and intramuscular EMG processing.  `pyemgpipeline` implements internationally accepted EMG processing conventions and provides a high-level interface for ensuring user adherence to those conventions, in terms of (1) processing parameter values, (2) processing steps, and (3) processing step order.  The international standards are from surface EMG for non-invasive assessment of muscles (SENIAM) [@Stegeman:1999].  The seven processing steps included in the package are DC offset removal, bandpass filtering, full wave rectification, linear envelope, end frame cutting, amplitude normalization, and segmentation.  In sport tasks particularly, it has been observed that amplitudes greater that 100% maximum voluntary contraction (MVC) can be observed [@Devaprakash:2016].  Therefore, we will be using amplitude normalization recommendations from @Devaprakash:2016 for amplitude normalization – this method guarantees a muscle will not exceed 100% MVC as all EMG trials from the experiment are used to identify maximal muscle activation.  What is not included in the package is time event detection and time normalization as different laboratories are interested in different phases of gait and prefer to use either linear or nonlinear time normalization techniques.

As stated by @Stegeman:1999 although EMG is easy to use it is also easy to abuse.  Researchers have thus tried to standardize the use of EMG regarding implementation, analysis, and reporting [@Stegeman:1999; @Merletti:1999; @JEK:2018]. These fine works however have arguably been unsuccessfully adopted by many researchers to date.  Based on the published standards and decades of practical experience, we summarize EMG signal processing conventions into the seven steps mentioned above and implement them within an easy-to-use package so that researchers and clinicians with all levels of experience process their EMG signals using correct conventions.  In addition to the processing steps and their order, the choice of processing/filtering parameters can be changed to match the needs of the research or clinician.  For example, within the bandpass filtering step, proper values of the low and high cutoff frequencies of the Butterworth bandpass filter depend on different use cases.  We carefully set default values in the package, and we also include suggested values under different scenarios in [related API documentation](https://aalhossary.github.io/pyemgpipeline/api/pyemgpipeline.processors.html#bandpassfilter) [@JEK:2018; @Merletti:1999; @Stegeman:1999; @Drake:2006].

After the aforementioned seven processing steps, sometimes it is necessary to perform the time normalization step before further analysis.  Time normalization can be executed by an external package `mwarp1d` [@Pataky:2019].

Within this package, class [DataProcessingManager](https://aalhossary.github.io/pyemgpipeline/api/pyemgpipeline.wrappers.html#dataprocessingmanager) is designed for high-level, guided processing, and users are encouraged to use it to adhere to accepted EMG processing conventions.

The depiction of `pyemgpipeline` data processing flow is shown in \autoref{fig:example}, which includes the original signal and the processed signals of all processing steps.  Please refer to the [full documentation](https://aalhossary.github.io/pyemgpipeline/) and the [source code](https://github.com/aalhossary/pyemgpipeline) for detailed information.

# Statement of need

**Research purpose**: `pyemgpipeline` aims to provide software tools for electromyography (EMG) data processing, while ensuring adherence to internationally accepted EMG processing conventions.

**Problem solved**: `pyemgpipeline` implements internationally accepted EMG processing conventions and provides a high-level interface for ensuring user adherence to those conventions. It facilitates the convenience and correctness of processing EMG data. To our best knowledge, no other package provides tools of EMG processing pipeline to ensure users adhere to accepted conventions in terms of (1) processing parameter values, (2) processing steps, (3) processing step order, and (4) amplitude normalization.

**Target audience**: The target audience is anyone working with surface or intramuscular EMG data such as gait biomechanics, sports science, rehabilitation, and robotics.

# Figures

![Depiction of `pyemgpipeline` data processing flow. (x-axis: seconds, y-axis: amplitude.)
\label{fig:example}](example_fig_in_paper.png){ height=80% }

# References
Examples
========================================

.. toctree::
   :maxdepth: 1

   notebooks/ex0_input_data_description.ipynb
   notebooks/ex1_EMGMeasurement.ipynb
   notebooks/ex2_EMGMeasurementCollection.ipynb
   notebooks/ex3_DataProcessingManager.ipynb
   notebooks/ex4_DataProcessingManager.ipynb
pyemgpipeline Documentation
========================================

**pyemgpipeline** is an electromyography (EMG) signal processing pipeline package.
It implements internationally accepted EMG processing conventions
and provides a high-level interface for ensuring user adherence to those conventions,
in terms of (1) processing parameter values, (2) processing steps, and (3) processing step order.
It is suitable for both surface EMG and intramuscular EMG.

Within this package, class :ref:`DataProcessingManager` is
designed for high-level, guided processing,
and users are encouraged to use it to adhere to accepted EMG processing conventions.

The source code is available in the
`GitHub repository <https://github.com/aalhossary/pyemgpipeline>`_.

Contents
----------------------------------------

.. toctree::
   :maxdepth: 2

   installation.rst
   quickstart.rst
   examples.rst
   api.rst
Installation
========================================

**pyemgpipeline** can be installed from the PyPI repository:

.. code::

   pip install pyemgpipeline

The dependencies of this package include
`numpy <https://numpy.org/>`_,
`scipy <https://scipy.org/>`_,
and `matplotlib <https://matplotlib.org/>`_.

After installation, the package can be launched by:

.. code:: python

   >>> import pyemgpipeline as pep
API
========================================

.. toctree::
   :maxdepth: 2

   api/pyemgpipeline.plots.rst
   api/pyemgpipeline.processors.rst
   api/pyemgpipeline.wrappers.rst
Quick Start
========================================

Input data format
----------------------------------------

Signal data of each trial should be stored as a 2d ndarray with shape
*(n_samples, n_channels)*,
where each column represents data of one channel.
If only one channel is presented, it can also be stored as
a 1d ndarray with shape *(n_samples,)*.

Timestamp data of each trial should be stored as a 1d ndarray with shape
*(n_samples,)*.
If timestamp data is not provided, it will be generated automatically,
starting from 0 and in increments of *1/hz*,
where *hz* is the sample rate.

Data of multiple trials should be organized in a *list*.

Demo code
----------------------------------------

The demo code simply creates example signal data with shape *(20,)*,
assumes sample rate 1000 Hz, applies seven processing steps to the data
using class :ref:`EMGMeasurement`, and plots the final result.

.. plot::
   :align: center
   :include-source:

   >>> import numpy as np
   >>> import pyemgpipeline as pep
   >>> data = np.array([20.3, 41.0, 53.9, 63.3, 39.5, 24.9, 26.1, 24.0, 44.1, 42.0,
   ...                  37.4, 24.6, -21.8, -56.3, -48.1, -45.0, -29.1, -9.6, 5.3, 1.4])
   >>> hz = 1000
   >>> m = pep.wrappers.EMGMeasurement(data, hz=hz)
   >>> m.apply_dc_offset_remover()
   >>> m.apply_bandpass_filter()
   >>> m.apply_full_wave_rectifier()
   >>> m.apply_linear_envelope()
   >>> m.apply_end_frame_cutter(n_end_frames=2)
   >>> m.apply_amplitude_normalizer(max_amplitude=8.5)
   >>> m.apply_segmenter(beg_ts=0, end_ts=0.015)
   >>> m.plot()
pyemgpipeline.plots
========================================

.. py:currentmodule:: pyemgpipeline.plots

Summary
----------------------------------------

Class
++++++++++++++++++++

.. autosummary::
   :nosignatures:

   EMGPlotParams

Functions
++++++++++++++++++++

.. autosummary::
   :nosignatures:

   plot_emg
   plot_emg_overlapping_trials

EMGPlotParams
----------------------------------------

.. autoclass:: EMGPlotParams
   :members:
   :undoc-members:

plot_emg
----------------------------------------

.. autofunction:: plot_emg

plot_emg_overlapping_trials
----------------------------------------

.. autofunction:: plot_emg_overlapping_trials
pyemgpipeline.wrappers
========================================

.. py:currentmodule:: pyemgpipeline.wrappers

Summary
----------------------------------------

Classes
++++++++++++++++++++

.. autosummary::
   :nosignatures:

   DataProcessingManager
   EMGMeasurement
   EMGMeasurementCollection

.. _DataProcessingManager:

DataProcessingManager
----------------------------------------

.. autoclass:: DataProcessingManager
   :members:
   :undoc-members:

Example
++++++++++++++++++++

.. plot:: api/ex_dpm.py
   :align: center
   :include-source:

.. _EMGMeasurement:

EMGMeasurement
----------------------------------------

.. autoclass:: EMGMeasurement
   :members:
   :undoc-members:

Example
++++++++++++++++++++

.. plot:: api/ex_em.py
   :scale: 80
   :align: center
   :include-source:

.. _EMGMeasurementCollection:

EMGMeasurementCollection
----------------------------------------

.. autoclass:: EMGMeasurementCollection
   :members:
   :undoc-members:
   :special-members: __getitem__

Example
++++++++++++++++++++

.. plot:: api/ex_emc.py
   :scale: 80
   :align: center
   :include-source:
﻿pyemgpipeline.processors
========================================

.. py:currentmodule:: pyemgpipeline.processors

Summary
----------------------------------------

Base Class
++++++++++++++++++++

.. autosummary::
   :nosignatures:

   BaseProcessor

Classes
++++++++++++++++++++

.. autosummary::
   :nosignatures:

   AmplitudeNormalizer
   BandpassFilter
   DCOffsetRemover
   EndFrameCutter
   FullWaveRectifier
   LinearEnvelope
   Segmenter

AmplitudeNormalizer
----------------------------------------

.. autoclass:: AmplitudeNormalizer
   :members:
   :undoc-members:
   :show-inheritance:

BandpassFilter
----------------------------------------

.. autoclass:: BandpassFilter
   :members:
   :undoc-members:
   :show-inheritance:

BaseProcessor
----------------------------------------

.. autoclass:: BaseProcessor
   :members:
   :undoc-members:

DCOffsetRemover
----------------------------------------

.. autoclass:: DCOffsetRemover
   :members:
   :undoc-members:
   :show-inheritance:

EndFrameCutter
----------------------------------------

.. autoclass:: EndFrameCutter
   :members:
   :undoc-members:
   :show-inheritance:

FullWaveRectifier
----------------------------------------

.. autoclass:: FullWaveRectifier
   :members:
   :undoc-members:
   :show-inheritance:

LinearEnvelope
----------------------------------------

.. autoclass:: LinearEnvelope
   :members:
   :undoc-members:
   :show-inheritance:

Segmenter
----------------------------------------

.. autoclass:: Segmenter
   :members:
   :undoc-members:
   :show-inheritance:
