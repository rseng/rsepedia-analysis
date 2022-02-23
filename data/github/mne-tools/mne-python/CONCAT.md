# Security Policy

## Supported Versions

New minor versions of MNE-Python are typically released twice per year.
Only the most current stable release is officially supported.
The unreleased, unstable "dev version" is also supported, though users
should beware that the API of the dev version is subject to change
without a proper 6-month deprecation cycle.

| Version | Supported                |
| ------- | ------------------------ |
| 0.25.x  | :heavy_check_mark: (dev) |
| 0.24.x  | :heavy_check_mark:       |
| < 0.24  | :x:                      |

## Reporting a Vulnerability

MNE-Python is software for analysis and visualization of brain activity
recorded with a variety of devices/modalities (EEG, MEG, ECoG, fNIRS, etc).
It is not expected that using MNE-Python will lead to security
vulnerabilities under normal use cases (i.e., running without administrator
privileges). However, if you think you have found a security vulnerability
in MNE-Python, **please do not report it as a GitHub issue**, in order to 
keep the vulnerability confidential. Instead, please report it to
mne-core-dev-team@groups.io and include a description and proof-of-concept
that is [short and self-contained](http://www.sscce.org/).

Generally you will receive a response within one week. MNE-Python does not
award bounties for security vulnerabilities.
Contributing to MNE-Python
==========================

MNE-Python is maintained by a community of scientists and research labs. The project accepts contributions in the form of bug reports, fixes, feature additions, and documentation improvements (including typo corrections). The best way to start contributing is by [opening an issue](https://github.com/mne-tools/mne-python/issues/new/choose) on our GitHub page to discuss ideas for changes or enhancements, or to tell us about behavior that you think might be a bug. For *general troubleshooting* or *usage questions*, please consider posting your questions on our [MNE Forum](https://mne.discourse.group).

Users and contributors to MNE-Python are expected to follow our [code of conduct](https://github.com/mne-tools/.github/blob/main/CODE_OF_CONDUCT.md).

The [contributing guide](https://mne.tools/dev/install/contributing.html) has details on the preferred contribution workflow
and the recommended system configuration for a smooth contribution/development experience.
# MNE-Python logos

The following logos used for MNE-Python are made available under a
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
license, see the [LICENSE](./LICENSE) file in this directory:

- The main MNE-Python logo: [mne_logo.svg](../doc/_static/mne_logo.svg)
- The hexagonal version of the MNE-Python logo: [logo_hex.svg](./logo_hex.svg)
- The small form of the MNE-Python logo: [mne_logo_small.svg](../doc/_static/mne_logo_small.svg)
- The MNE-Python "helmet": [mne_helmet.png](../doc/_static/mne_helmet.png)
- The MNE-Python "helmet" thumbnail: [favicon.ico](../doc/_static/favicon.ico)

See also the code to create some of these logos:

- [generate_mne_logos.py](./generate_mne_logos.py)
Thanks for contributing a pull request! Please make sure you have read the
[contribution guidelines](https://mne.tools/dev/install/contributing.html)
before submitting.

Please be aware that we are a loose team of volunteers so patience is
necessary. Assistance handling other issues is very welcome. We value
all user contributions, no matter how minor they are. If we are slow to
review, either the pull request needs some benchmarking, tinkering,
convincing, etc. or more likely the reviewers are simply busy. In either
case, we ask for your understanding during the review process.

Again, thanks for contributing!

#### Reference issue
Example: Fixes #1234.


#### What does this implement/fix?
Explain your changes.


#### Additional information
Any additional information you think is important.
# Contributing to MNE-Python

First off, thanks for taking the time to contribute!

The following is a quick summary to a set of guidelines for contributing to [MNE-Python](https://github.com/mne-tools/mne-python) on GitHub. These are mostly guidelines, not rules. Use your best judgment, and feel free to propose changes to this document in a pull request.

## Code of Conduct

This project and everyone participating in it is governed by the [MNE-Python's Code of Conduct](https://github.com/mne-tools/.github/blob/main/CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to [mne-conduct@googlegroups.com](mailto:mne-conduct@googlegroups.com).

## How to contribute

Before contributing make sure you are familiar with [our contributing guide](https://mne.tools/dev/install/contributing.html).
---
name: Bug report
about: Tell us about broken, incorrect, or inconsistent behavior.
title: ''
labels: BUG
assignees: ''
---

**READ THIS FIRST:** If you are having trouble getting MNE-Python to work with
your own data, you should ask for help on the
[MNE Forum](https://mne.discourse.group).

Our GitHub issue tracker is only used to report bugs and suggest improvements
to MNE-Python. For any other questions, please use the forum.
Usage questions that are posted as GitHub issues are usually closed without
being answered. See
[the FAQ entry on filing bug reports](https://mne.tools/dev/overview/faq.html#i-think-i-found-a-bug-what-do-i-do)
for more guidance. If you're pretty sure your problem is a bug in MNE-Python,
please **delete this section** and fill in the headings below, replacing the
placeholder text with your own responses. Surround any code samples with triple
backticks above and below the code block (see
[the GitHub guide to markdown](https://guides.github.com/features/mastering-markdown/#GitHub-flavored-markdown)
for help with issue formatting). Alternatively, you can put your MWE in a
[public gist](https://gist.github.com) and link to it in this issue.


#### Describe the bug
*Replace this text with a description of the bug.*


#### Steps to reproduce
*Replace this text with a code snippet or minimal working example [MWE] to
replicate your problem, using one of the [built-in datasets], preferably the
one called [sample]. If you can't replicate on a built-in dataset, provide also
a link to a small, anonymized portion of your data that does yield the error.*

[MWE]: https://en.wikipedia.org/wiki/Minimal_Working_Example
[built-in datasets]: https://mne.tools/dev/overview/datasets_index.html
[sample]: https://mne.tools/dev/overview/datasets_index.html#sample


#### Expected results
*Replace this text with a description of what you expected to happen.*


#### Actual results
*Replace this text with the actual output, traceback, screenshot, or other
description of the results.*


#### Additional information
*Replace this text with the output of `mne.sys_info()`.*
---
name: Feature request
about: Suggest additions or improvements to MNE-Python.
title: ''
labels: ENH
assignees: ''
---

#### Describe the new feature or enhancement
*Please provide a clear and concise description of what you want to add or
change.*


#### Describe your proposed implementation
*Describe how you think the feature or improvement should be implemented (e.g., 
as a new method on an existing class? as new capability added to an existing
method?) If you're not sure, please delete this section and the next section.*


#### Describe possible alternatives
*If you've suggested an implementation above, list here any alternative
implementations you can think of, and brief comments explaining why the chosen
implementation is better.*


#### Additional comments
*Add any other context or screenshots about the feature request here.*
---
name: Documentation
about: Request a new tutorial, an improvement to an existing tutorial, or a new term in our glossary.
title: ''
labels: DOC
assignees: ''
---

*Please make the issue title a succinct description of your proposed change.
For example: "add 'foo' to the glossary" or "add my favorite visualization to
the ERP tutorial".*

#### Proposed documentation ehancement
*Describe your proposed enhancement in detail. If you are requesting a new
glossary term, please include a proposed definition.*
.. -*- mode: rst -*-

|GH-Linux|_ |GH-macOS|_ |Azure|_ |Circle|_ |Codecov|_ |PyPI|_ |conda-forge|_ |Zenodo|_

|MNE|_

.. |GH-Linux| image:: https://github.com/mne-tools/mne-python/workflows/linux%20/%20conda/badge.svg?branch=main
.. _GH-Linux: https://github.com/mne-tools/mne-python/actions?query=branch:main+event:push

.. |GH-macOS| image:: https://github.com/mne-tools/mne-python/workflows/macos%20/%20conda/badge.svg?branch=main
.. _GH-macOS: https://github.com/mne-tools/mne-python/actions?query=branch:main+event:push

.. |Azure| image:: https://dev.azure.com/mne-tools/mne-python/_apis/build/status/mne-tools.mne-python?branchName=main
.. _Azure: https://dev.azure.com/mne-tools/mne-python/_build/latest?definitionId=1&branchName=main

.. |Circle| image:: https://circleci.com/gh/mne-tools/mne-python.svg?style=shield
.. _Circle: https://circleci.com/gh/mne-tools/mne-python

.. |Codecov| image:: https://codecov.io/gh/mne-tools/mne-python/branch/main/graph/badge.svg
.. _Codecov: https://codecov.io/gh/mne-tools/mne-python

.. |PyPI| image:: https://img.shields.io/pypi/dm/mne.svg?label=PyPI%20downloads
.. _PyPI: https://pypi.org/project/mne/

.. |conda-forge| image:: https://img.shields.io/conda/dn/conda-forge/mne.svg?label=Conda%20downloads
.. _conda-forge: https://anaconda.org/conda-forge/mne

.. |Zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.592483.svg
.. _Zenodo: https://doi.org/10.5281/zenodo.592483

.. |MNE| image:: https://mne.tools/stable/_static/mne_logo.svg
.. _MNE: https://mne.tools/dev/

MNE-Python
==========

`MNE-Python software`_ is an open-source Python package for exploring,
visualizing, and analyzing human neurophysiological data such as MEG, EEG, sEEG,
ECoG, and more. It includes modules for data input/output, preprocessing,
visualization, source estimation, time-frequency analysis, connectivity analysis,
machine learning, and statistics.


Documentation
^^^^^^^^^^^^^

`MNE documentation`_ for MNE-Python is available online.


Installing MNE-Python
^^^^^^^^^^^^^^^^^^^^^

To install the latest stable version of MNE-Python, you can use pip_ in a terminal:

.. code-block:: bash

    pip install -U mne

- MNE-Python 0.17 was the last release to support Python 2.7
- MNE-Python 0.18 requires Python 3.5 or higher
- MNE-Python 0.21 requires Python 3.6 or higher
- MNE-Python 0.24 requires Python 3.7 or higher

For more complete instructions and more advanced installation methods (e.g. for
the latest development version), see the `installation guide`_.


Get the latest code
^^^^^^^^^^^^^^^^^^^

To install the latest version of the code using pip_ open a terminal and type:

.. code-block:: bash

    pip install -U https://github.com/mne-tools/mne-python/archive/main.zip

To get the latest code using `git <https://git-scm.com/>`__, open a terminal and type:

.. code-block:: bash

    git clone https://github.com/mne-tools/mne-python.git

Alternatively, you can also download a
`zip file of the latest development version <https://github.com/mne-tools/mne-python/archive/main.zip>`__.


Dependencies
^^^^^^^^^^^^

The minimum required dependencies to run MNE-Python are:

- Python >= 3.7
- NumPy >= 1.18.1
- SciPy >= 1.4.1
- Matplotlib >= 3.1.0
- pooch >= 1.5
- tqdm
- Jinja2
- decorator

For full functionality, some functions require:

- Scikit-learn >= 0.22.0
- Numba >= 0.48.0
- NiBabel >= 2.5.0
- Pandas >= 1.0.0
- Picard >= 0.3
- CuPy >= 7.1.1 (for NVIDIA CUDA acceleration)
- DIPY >= 1.1.0
- Imageio >= 2.6.1
- PyVista >= 0.32
- pyvistaqt >= 0.4
- mffpy >= 0.5.7
- h5py
- h5io
- pymatreader

Contributing to MNE-Python
^^^^^^^^^^^^^^^^^^^^^^^^^^

Please see the documentation on the MNE-Python homepage:

https://mne.tools/dev/install/contributing.html


Forum
^^^^^^

https://mne.discourse.group


Licensing
^^^^^^^^^

MNE-Python is **BSD-licenced** (3 clause):

    This software is OSI Certified Open Source Software.
    OSI Certified is a certification mark of the Open Source Initiative.

    Copyright (c) 2011-2022, authors of MNE-Python.
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    * Neither the names of MNE-Python authors nor the names of any
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.

    **This software is provided by the copyright holders and contributors
    "as is" and any express or implied warranties, including, but not
    limited to, the implied warranties of merchantability and fitness for
    a particular purpose are disclaimed. In no event shall the copyright
    owner or contributors be liable for any direct, indirect, incidental,
    special, exemplary, or consequential damages (including, but not
    limited to, procurement of substitute goods or services; loss of use,
    data, or profits; or business interruption) however caused and on any
    theory of liability, whether in contract, strict liability, or tort
    (including negligence or otherwise) arising in any way out of the use
    of this software, even if advised of the possibility of such
    damage.**


.. _MNE-Python software: https://mne.tools/dev/
.. _MNE documentation: https://mne.tools/dev/overview/index.html
.. _installation guide: https://mne.tools/dev/install/index.html
.. _pip: https://pip.pypa.io/en/stable/

.. _api_decoding:

Decoding
========

:py:mod:`mne.decoding`:

.. automodule:: mne.decoding
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/

   CSP
   EMS
   FilterEstimator
   LinearModel
   PSDEstimator
   Scaler
   TemporalFilter
   TimeFrequency
   UnsupervisedSpatialFilter
   Vectorizer
   ReceptiveField
   TimeDelayingRidge
   SlidingEstimator
   GeneralizingEstimator
   SPoC
   SSD

Functions that assist with decoding and model fitting:

.. autosummary::
   :toctree: generated/

   compute_ems
   cross_val_multiscore
   get_coef
.. _api_reference:

====================
Python API Reference
====================

This is the reference for classes (``CamelCase`` names) and functions
(``underscore_case`` names) of MNE-Python, grouped thematically by analysis
stage. Functions and classes that are not
below a module heading are found in the ``mne`` namespace.

MNE-Python also provides multiple command-line scripts that can be called
directly from a terminal, see :ref:`python_commands`.

.. container:: d-none

   :py:mod:`mne`:

   .. automodule:: mne
      :no-members:
      :no-inherited-members:

.. toctree::
    :maxdepth: 2

    most_used_classes
    reading_raw_data
    file_io
    creating_from_arrays
    export
    datasets
    visualization
    preprocessing
    events
    sensor_space
    covariance
    mri
    forward
    inverse
    source_space
    time_frequency
    connectivity
    statistics
    simulation
    decoding
    realtime
    report
    logging
Reading raw data
================

:py:mod:`mne.io`:

.. currentmodule:: mne.io

.. automodule:: mne.io
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/

   anonymize_info
   read_raw
   read_raw_artemis123
   read_raw_bti
   read_raw_cnt
   read_raw_ctf
   read_raw_curry
   read_raw_edf
   read_raw_bdf
   read_raw_gdf
   read_raw_kit
   read_raw_nedf
   read_raw_nicolet
   read_raw_hitachi
   read_raw_nirx
   read_raw_snirf
   read_raw_eeglab
   read_raw_brainvision
   read_raw_egi
   read_raw_fif
   read_raw_eximia
   read_raw_fieldtrip
   read_raw_boxy
   read_raw_persyst
   read_raw_nihon

Base class:

.. autosummary::
   :toctree: generated

   BaseRaw

:py:mod:`mne.io.kit`:

.. currentmodule:: mne.io.kit

.. automodule:: mne.io.kit
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/

   read_mrk
Most-used classes
=================

.. currentmodule:: mne

.. autosummary::
   :toctree: generated/

   io.Raw
   Epochs
   Evoked
   Info
.. _whats_new:

What's new
==========

.. currentmodule:: mne

.. include:: changes/latest.inc
.. include:: changes/0.24.inc
.. include:: changes/0.23.inc
.. include:: changes/0.22.inc
.. include:: changes/0.21.inc
.. include:: changes/0.20.inc
.. include:: changes/0.19.inc
.. include:: changes/0.18.inc
.. include:: changes/0.17.inc
.. include:: changes/0.16.inc
.. include:: changes/0.15.inc
.. include:: changes/0.14.inc
.. include:: changes/0.13.inc
.. include:: changes/0.12.inc
.. include:: changes/0.11.inc
.. include:: changes/0.10.inc
.. include:: changes/0.9.inc
.. include:: changes/0.8.inc
.. include:: changes/0.7.inc
.. include:: changes/0.6.inc
.. include:: changes/0.5.inc
.. include:: changes/0.4.inc
.. include:: changes/0.3.inc
.. include:: changes/0.2.inc
.. include:: changes/0.1.inc
.. include:: changes/names.inc
.. include:: links.inc

Sensor Space Data
=================

.. currentmodule:: mne

.. autosummary::
   :toctree: generated/

   combine_evoked
   concatenate_raws
   equalize_channels
   grand_average
   pick_channels
   pick_channels_cov
   pick_channels_forward
   pick_channels_regexp
   pick_types
   pick_types_forward
   pick_info
   read_epochs
   read_reject_parameters
   read_vectorview_selection
   rename_channels

:py:mod:`mne.baseline`:

.. automodule:: mne.baseline
   :no-members:
   :no-inherited-members:

.. currentmodule:: mne.baseline

.. autosummary::
   :toctree: generated/

   rescale

MNE-Report
==========

:py:mod:`mne`:

.. currentmodule:: mne

.. autosummary::
   :toctree: generated/

   Report
   open_report

MRI Processing
==============

.. currentmodule:: mne

Step by step instructions for using :func:`gui.coregistration`:

 - `Coregistration for subjects with structural MRI
   <https://www.slideshare.net/mne-python/mnepython-coregistration>`_
 - `Scaling a template MRI for subjects for which no MRI is available
   <https://www.slideshare.net/mne-python/mnepython-scale-mri>`_

.. autosummary::
   :toctree: generated/

   coreg.get_mni_fiducials
   coreg.estimate_head_mri_t
   io.read_fiducials
   io.write_fiducials
   get_montage_volume_labels
   gui.coregistration
   gui.locate_ieeg
   create_default_subject
   head_to_mni
   head_to_mri
   read_freesurfer_lut
   read_lta
   read_talxfm
   scale_mri
   scale_bem
   scale_labels
   scale_source_space
   transforms.apply_volume_registration
   transforms.compute_volume_registration
   vertex_to_mni
   warp_montage_volume
   coreg.Coregistration
Glossary
========

.. currentmodule:: mne

The Glossary provides short definitions of vocabulary specific to MNE-Python and
general neuroimaging concepts. If you think a term is missing, please consider
`creating a new issue`_ or `opening a pull request`_ to add it.

.. glossary::
    :sorted:


    annotations
        An annotation is defined by an onset, a duration, and a textual
        description. It can contain information about the experiment, but
        also details on signals marked by a human such as bad data segments,
        sleep stages, sleep events (spindles, K-complex), and so on.
        An :class:`Annotations` object is a container for multiple annotations,
        which is available as the ``annotations`` attribute of :class:`~io.Raw`
        objects. See :class:`Annotations` for the class definition and
        :ref:`tut-events-vs-annotations` for a short tutorial.
        See also :term:`events`.

    array-like
        Something that acts like – or can be converted to – a
        :class:`NumPy array <numpy.ndarray>`.
        This includes (but is not limited to)
        :class:`arrays <numpy.ndarray>`, `lists <list>`, and
        `tuples <tuple>`.

    beamformer
        A beamformer is a popular source estimation approach that uses a set of
        spatial filters (beamformer weights) to compute time courses of sources
        at predefined locations. See :class:`beamformer.Beamformer` for the class
        definition. See also :term:`LCMV`.

    BEM
    boundary element model
    boundary element method
        BEM is the acronym for boundary element method or boundary element
        model. Both are related to the definion of the conductor model in the
        forward model computation. The boundary element model consists of surfaces
        such as the inner skull, outer skull, and outer skin (scalp) that define
        compartments of tissues of the head. You can compute the BEM surfaces with
        :func:`bem.make_watershed_bem` or :func:`bem.make_flash_bem`.
        See :ref:`tut-forward` for a usage demo.

    channels
        Channels refer to MEG sensors, EEG electrodes or other sensors such as
        EOG, ECG, sEEG, ECoG, etc. Channels usually have
        a type (such as gradiometer), and a unit (such as T/m) used e.g. for
        plotting. See also :term:`data channels`.

    data channels
        Many functions in MNE-Python operate on "data channels" by default. These
        are channels that contain electrophysiological data from the brain,
        as opposed to other channel types such as EOG, ECG, stimulus/trigger,
        or acquisition system status data. The set of channels considered
        "data channels" in MNE contains the following types (together with scale
        factors for plotting):

        .. mne:: data channels list

    DICS
    dynamic imaging of coherent sources
        Dynamic Imaging of Coherent Sources is a method for computing source
        power in different frequency bands. See :ref:`ex-inverse-source-power`
        and :func:`beamformer.make_dics` for more details.

    digitization
        Digitization is a procedure of recording the head shape and locations of
        fiducial coils (or :term:`HPI`) and/or EEG electrodes on the head. They
        are represented as a set of points in 3D space.
        See :ref:`reading-dig-montages` and :ref:`dig-formats`.

    dipole
    ECD
    equivalent current dipole
        An equivalent current dipole (ECD) is an approximate representation of
        post-synaptic activity in a small cortical region. The intracellular
        currents that give rise to measurable EEG/MEG signals are thought to
        originate in populations of cortical pyramidal neurons aligned
        perpendicularly to the cortical surface. Because the length of such
        current sources is very small relative to the distance between the
        cortex and the EEG/MEG sensors, the fields measured by these techniques
        are well approximated by (i.e., equivalent to) fields generated by
        idealized point sources (dipoles) located on the cortical surface.

    dSPM
    dynamic statistical parametric mapping
        Dynamic statistical parametric mapping (dSPM) gives a noise-normalized
        minimum-norm estimate at a given source location. It is calculated by
        dividing the activity estimate at each source location by the baseline
        standard deviation of the noise.

    eLORETA
    sLORETA
        eLORETA and sLORETA (exact and standardized low resolution brain
        electromagnetic tomography) are linear source estimation techniques
        like :term:`dSPM` and :term:`MNE`. sLORETA outputs
        standardized values (like dSPM), while eLORETA generates normalized
        current estimates. See :func:`minimum_norm.apply_inverse`,
        :ref:`tut-inverse-methods`, and :ref:`example-sLORETA`.

    epochs
        Epochs (sometimes called "trials" in other software packages) are
        equal-length segments of data extracted from continuous data. Usually,
        epochs are extracted around stimulus events or responses,
        though sometimes sequential or overlapping epochs are used (e.g.,
        for analysis of resting-state activity). See :class:`Epochs` for the
        class definition and :ref:`tut-epochs-class` for a narrative overview.

    events
        Events correspond to specific time points in raw data, such as triggers,
        experimental condition events, etc. MNE-Python represents events with
        integers stored in NumPy arrays of shape ``(n_events, 3)``. The first
        column contains the event onset (in samples) with :term:`first_samp`
        included. The last column contains the event code. The second
        column contains the signal value of the immediately preceding sample,
        and reflects the fact that event arrays sometimes originate from
        analog voltage channels ("trigger channels" or "stim channels"). In
        most cases, the second column is all zeros and can be ignored.
        Event arrays can be created with :func:`mne.make_fixed_length_events`,
        :func:`mne.read_events`, and :func:`mne.find_events`.
        See :ref:`tut-events-vs-annotations` for a short tutorial.
        See also :term:`events`.

    evoked
        Evoked data are obtained by averaging epochs. Typically, an evoked object
        is constructed for each subject and each condition, but it can also be
        obtained by averaging a list of evoked objects over different subjects.
        See :class:`EvokedArray` for the class definition and
        :ref:`tut-evoked-class` for a narrative overview.

    fiducial
    fiducial point
    anatomical landmark
        Fiducials are objects placed in the field of view of an imaging system
        to act as known spatial references that are easy to localize.
        In neuroimaging, fiducials are often placed on anatomical landmarks
        such as the nasion (NAS) or left/right preauricular points (LPA and
        RPA).

        These known reference locations are used to define a coordinate system
        for localizing sensors (hence NAS, LPA and RPA are often
        called "cardinal points" because they define the cardinal directions of
        the head coordinate system). The cardinal points are also useful when
        co-registering measurements in different coordinate systems (such as
        aligning EEG sensor locations to an MRI of the head).

        Due to the common neuroimaging practice of placing fiducial objects on
        anatomical landmarks, the terms "fiducial", "anatomical landmark", and
        "cardinal point" are often (erroneously) used interchangeably.

    first_samp
        The :attr:`~io.Raw.first_samp` attribute of :class:`~io.Raw`
        objects is an integer representing the number of time samples that
        passed between the onset of the hardware acquisition system and the
        time when data recording started. This approach to sample
        numbering is a peculiarity of VectorView MEG systems, but for
        consistency it is present in all :class:`~io.Raw` objects
        regardless of the source of the data. In other words,
        :attr:`~io.Raw.first_samp` will be ``0`` in :class:`~io.Raw`
        objects loaded from non-VectorView data files. See also
        :term:`last_samp`.

    forward
    forward solution
        The forward solution is a linear operator capturing the
        relationship between each dipole location in the :term:`source space`
        and the corresponding field distribution measured by the sensors
        (the "lead field matrix"). Calculating a forward solution requires a
        conductivity model of the head, which encapsulates the geometries and
        electrical conductivities of the different tissue compartments (see
        :term:`boundary element model` and :class:`bem.ConductorModel`).

    GFP
    global field power
        Global Field Power (GFP) is a measure of the (non-)uniformity
        of the electromagnetic field at the sensors. It is typically calculated
        as the standard deviation of the sensor values at each time point. Thus,
        it is a one-dimensional time series capturing the spatial variability
        of the signal across sensor locations.

    HED
    hierarchical event descriptors
        Hierarchical event descriptors (HED) are tags that use
        keywords separated by slashes (/) to describe different types of
        experimental events (for example, ``stimulus/circle/red/left`` and
        ``stimulus/circle/blue/left``). These tags can be used to group
        experimental events and select event types for analysis.

    HPI
    cHPI
    head position indicator
        Head position indicators (HPI, sometimes cHPI for
        *continuous* head position indicators) are small coils attached to a
        subject's head during MEG acquisition. Each coil emits a sinusoidal
        signal of a different frequency, which is picked up by the MEG sensors
        and can be used to infer the head position. With cHPI, the sinusoidal
        signals are typically set at frequencies above any neural signal of
        interest, and thus can be removed after head position correction via
        low-pass filtering. See :ref:`tut-head-pos`.

    info
    measurement info
        A "measurement info" (or short "info") object is a collection of metadata
        related to :class:`~io.Raw`, :class:`Epochs`, or :class:`Evoked`
        objects. It contains channel locations and types, sampling frequency,
        preprocessing history such as filters, etc.
        See :ref:`tut-info-class` for a narrative overview.

    inverse
    inverse operator
        The inverse operator is an :math:`M \times N` matrix (:math:`M` source
        locations by :math:`N` sensors) that, when applied to the sensor
        signals, yields estimates of the brain activity that gave rise to the
        observed sensor signals. Inverse operators are available for the linear
        inverse methods :term:`MNE`, :term:`dSPM`, :term:`sLORETA`, and
        :term:`eLORETA`. See :func:`minimum_norm.apply_inverse`.

    label
        A :class:`Label` refers to a defined region in the cortex, often called
        a region of interest (ROI) in the literature. Labels can be defined
        anatomically (based on the physical structure of the cortex) or functionally
        (based on cortical responses to specific stimuli). See also :term:`ROI`.

    last_samp
        The :attr:`~io.Raw.last_samp` attribute of :class:`~io.Raw`
        objects is an integer representing the number of time samples that
        passed between the start and end of data recording. This approach to sample
        numbering is a peculiarity of VectorView MEG systems, but for
        consistency it is present in all :class:`~io.Raw` objects
        regardless of the source of the data. See also :term:`first_samp`.

    layout
        A :class:`~channels.Layout` gives sensor positions in two
        dimensions (defined by ``x``, ``y``, ``width``, and ``height`` values for
        each sensor). It is primarily used for illustrative purposes (i.e., making
        diagrams of approximate sensor positions in cartoons of the head,
        so-called topographies or topomaps). See also :term:`montage`.

    LCMV
    LCMV beamformer
        Linearly constrained minimum variance beamformer attempt to
        estimate activity for a given source while suppressing cross-talk from
        other regions (:func:`beamformer.make_lcmv`). See also
        :term:`beamformer`.

    FreeSurfer LUT
    LUT
        A FreeSurfer lookup table (LUT) provides a mapping between a given
        volumetric atlas or surface label name, its integer value
        (e.g., in ``aparc+aseg.mgz``), and its standard color (see the
        `FreeSurfer wiki <https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT>`__
        for more information). Custom LUTs can be also be created from different
        surface parcellations, see for example `this comment about HCPMMP
        <https://github.com/mne-tools/mne-python/pull/7639#issuecomment-625907891>`__.

    maximum intensity projection
        A method to display pixel-wise activity within some volume by
        finding the maximum value along a vector from the viewer to the pixel
        (i.e., along the vector pependicular to the view plane).

    MNE
    minimum-norm estimate
    minimum-norm estimation
        Minimum-norm estimation (MNE) can be used to generate a distributed
        map of activation on a :term:`source space` (usually on a cortical surface).
        MNE uses a linear :term:`inverse operator` to project sensor measurements
        into the source space. The :term:`inverse operator` is computed from the
        :term:`forward solution` for a subject and an estimate of the
        :term:`noise covariance` of sensor measurements.

    montage
        EEG channel names and relative positions of sensors on the scalp.
        While layouts are 2D locations, montages are 3D locations. A montage
        can also contain locations for HPI points, fiducial points, or
        extra head shape points.
        See :class:`~channels.DigMontage` for the class definition. See also
        :term:`layout`.

    morphing
        Morphing refers to the operation of transferring source estimates from
        one anatomy to another. It is known as realignment in the fMRI
        literature. This operation is necessary for group studies to get the
        data into a common space for statistical analysis.
        See :ref:`ch_morph` for more details.

    noise covariance
        The noise covariance is a matrix that contains the covariance between data
        channels. It is a square matrix with shape ``n_channels`` :math:`\times`
        ``n_channels``. It is especially useful when working with multiple sensor
        types (e.g. EEG and MEG). In practice, the matrix is estimated from baseline
        periods or empty room measurements, and it also provides a noise model
        that can be used for subsequent analysis (like source imaging).

    path-like
        Something that acts like a path in a file system. This can be a `str`
        or a `pathlib.Path`.

    pick
        An integer that is the index of a channel in the :term:`measurement info`.
        It allows to obtain the information on a channel in the list of channels
        available in ``info['chs']``.

    projector
    SSP
        A projector, also referred to as Signal Space
        Projection (SSP), defines a linear operation applied spatially to EEG
        or MEG data. A matrix multiplication of an SSP projector with the data
        will reduce the rank of the data by projecting it to a
        lower-dimensional subspace. Such projections are typically applied to
        both the data and the forward operator when performing
        source localization. Note that EEG average referencing can be done
        using such a projection operator. Projectors are stored alongside data
        in the :term:`measurement info` in the field ``info['projs']``.

    raw
        `~io.Raw` objects hold continuous data (preprocessed or not), typically
        obtained from reading recordings stored in a file.
        See :class:`~io.RawArray` for the class definition and :ref:`tut-raw-class`
        for a narrative overview.

    RAS
        Right-Anterior-Superior, denoting the standard way to define coordinate
        frames in MNE-Python:

        R
            +X is right, -X is left
        A
            +Y is anterior (front), -Y is posterior (rear)
        S
            +Z is superior (top), -Z is inferior (bottom)

    ROI
    region of interest
        A spatial region where an experimental effect is expected to manifest.
        This can be a collection of sensors or, when performing inverse imaging,
        a set of vertices on the cortical surface or within the cortical volume.
        See also :term:`label`.

    selection
        A selection is a set of picked channels (for example, all sensors
        falling within a :term:`region of interest`).

    STC
    source estimate
    source time course
        Source estimates, commonly referred to as STC (Source Time Courses),
        are obtained from source localization methods such as :term:`dSPM`,
        :term:`sLORETA`, :term:`LCMV`, or MxNE.
        STCs contain the amplitudes of the neural sources over time.
        In MNE-Python, :class:`SourceEstimate` objects only store the
        amplitudes of activation but not the locations of the sources. The
        locations are stored separately in the :class:`SourceSpaces` object
        that was used to compute the forward operator.
        See :class:`SourceEstimate`, :class:`VolSourceEstimate`,
        :class:`VectorSourceEstimate`, and :class:`MixedSourceEstimate`.

    source space
        A source space specifies where in the brain source amplitudes are
        estimated. It corresponds to locations of a set of
        candidate :term:`equivalent current dipoles<ECD>`. MNE-Python mostly works
        with source spaces defined on the cortical surfaces estimated
        by FreeSurfer from a T1-weighted MRI image. See :ref:`tut-forward`
        to read about how to compute a forward operator in a source space.
        See :class:`SourceSpaces` for the class definition.

    stim channel
    trigger channel
        A stim channel or trigger channel is a channel that encodes
        events during the recording. It is typically a channel that is always
        zero and takes positive values when something happens (such as the
        onset of a stimulus or a subject response). Stim channels are often
        prefixed with ``STI`` to distinguish them from other channel types. See
        :ref:`stim-channel-defined` for more details.

    tfr
        A time-frequency representation (TFR) is often a spectrogram (STFT) or
        scaleogram (wavelet) showing the frequency content as a function of
        time.

    trans
        A coordinate frame affine transformation, usually between the Neuromag head
        coordinate frame and the MRI Surface RAS coordinate frame used by Freesurfer.

    whitening
        A linear operation that transforms data with a known covariance
        structure into "whitened data", which has a covariance structure equal to
        the identity matrix. In other words, whitening creates virtual channels that
        are uncorrelated and have unit variance. This is also known as a
        sphering transformation.

        The term "whitening" comes from the fact that light with a flat
        frequency spectrum in the visible range is white, whereas
        non-uniform frequency spectra lead to perception of different colors
        (e.g., "pink noise" has a ``1/f`` characteristic, which for visible
        light would appear pink).

.. LINKS

.. _`creating a new issue`:
   https://github.com/mne-tools/mne-python/issues/new?template=glossary.md
.. _`opening a pull request`:
   https://github.com/mne-tools/mne-python/pull/new/main

Visualization
=============

.. currentmodule:: mne.viz

:py:mod:`mne.viz`:

.. automodule:: mne.viz
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/

   Brain
   ClickableImage
   Figure3D
   add_background_image
   centers_to_edges
   compare_fiff
   circular_layout
   iter_topography
   mne_analyze_colormap
   plot_bem
   plot_brain_colorbar
   plot_chpi_snr
   plot_cov
   plot_channel_labels_circle
   plot_csd
   plot_dipole_amplitudes
   plot_dipole_locations
   plot_drop_log
   plot_epochs
   plot_epochs_psd_topomap
   plot_events
   plot_evoked
   plot_evoked_image
   plot_evoked_topo
   plot_evoked_topomap
   plot_evoked_joint
   plot_evoked_field
   plot_evoked_white
   plot_filter
   plot_head_positions
   plot_ideal_filter
   plot_compare_evokeds
   plot_ica_sources
   plot_ica_components
   plot_ica_properties
   plot_ica_scores
   plot_ica_overlay
   plot_epochs_image
   plot_layout
   plot_montage
   plot_projs_topomap
   plot_raw
   plot_raw_psd
   plot_sensors
   plot_snr_estimate
   plot_source_estimates
   link_brains
   plot_volume_source_estimates
   plot_vector_source_estimates
   plot_sparse_source_estimates
   plot_tfr_topomap
   plot_topo_image_epochs
   plot_topomap
   plot_alignment
   snapshot_brain_montage
   plot_arrowmap
   set_3d_backend
   get_3d_backend
   use_3d_backend
   set_3d_options
   set_3d_view
   set_3d_title
   create_3d_figure
   close_3d_figure
   close_all_3d_figures
   get_brain_class
   set_browser_backend
   get_browser_backend
   use_browser_backend

Inverse Solutions
=================

:py:mod:`mne.minimum_norm`:

.. automodule:: mne.minimum_norm
   :no-members:
   :no-inherited-members:

.. currentmodule:: mne.minimum_norm

.. autosummary::
   :toctree: generated/

   InverseOperator
   apply_inverse
   apply_inverse_cov
   apply_inverse_epochs
   apply_inverse_raw
   compute_source_psd
   compute_source_psd_epochs
   compute_rank_inverse
   estimate_snr
   make_inverse_operator
   prepare_inverse_operator
   read_inverse_operator
   source_band_induced_power
   source_induced_power
   write_inverse_operator
   make_inverse_resolution_matrix
   resolution_metrics
   get_cross_talk
   get_point_spread

:py:mod:`mne.inverse_sparse`:

.. automodule:: mne.inverse_sparse
   :no-members:
   :no-inherited-members:

.. currentmodule:: mne.inverse_sparse

.. autosummary::
   :toctree: generated/

   mixed_norm
   tf_mixed_norm
   gamma_map
   make_stc_from_dipoles

:py:mod:`mne.beamformer`:

.. automodule:: mne.beamformer
   :no-members:
   :no-inherited-members:

.. currentmodule:: mne.beamformer

.. autosummary::
   :toctree: generated/

   Beamformer
   read_beamformer
   make_lcmv
   apply_lcmv
   apply_lcmv_epochs
   apply_lcmv_raw
   apply_lcmv_cov
   make_dics
   apply_dics
   apply_dics_csd
   apply_dics_epochs
   rap_music
   make_lcmv_resolution_matrix

.. currentmodule:: mne

.. autosummary::
   :toctree: generated/

   Dipole
   DipoleFixed
   fit_dipole

:py:mod:`mne.dipole`:

.. automodule:: mne.dipole
   :no-members:
   :no-inherited-members:

.. currentmodule:: mne.dipole

.. autosummary::
   :toctree: generated/

   get_phantom_dipoles
.. title:: MNE

.. The page title must be in rST for it to show in next/prev page buttons.
   Therefore we add a special style rule to only this page that hides h1 tags

.. raw:: html

    <style type="text/css">h1 {display:none;}</style>

MNE-Python Homepage
===================

.. LOGO

.. image:: _static/mne_logo.svg
   :alt: MNE-Python
   :class: logo
   :align: center


.. rst-class:: h4 text-center font-weight-light my-4

   Open-source Python package for exploring, visualizing, and analyzing
   human neurophysiological data: MEG, EEG, sEEG, ECoG, NIRS, and more.

.. frontpage gallery is added by a conditional in _templates/layout.html

.. toctree::
   :hidden:

   Install<install/index>
   Documentation<overview/index>
   API Reference<python_reference>
   Get help<overview/get_help>
   Development<overview/development>

Events
======

.. currentmodule:: mne

.. autosummary::
   :toctree: generated/

   Annotations
   AcqParserFIF
   concatenate_events
   find_events
   find_stim_steps
   make_fixed_length_events
   make_fixed_length_epochs
   merge_events
   parse_config
   pick_events
   read_annotations
   read_events
   write_events
   concatenate_epochs
   events_from_annotations
   annotations_from_events

:py:mod:`mne.event`:

.. automodule:: mne.event
   :no-members:
   :no-inherited-members:

.. currentmodule:: mne.event

.. autosummary::
   :toctree: generated/

   define_target_events
   match_event_names
   shift_time_events

:py:mod:`mne.epochs`:

.. automodule:: mne.epochs
   :no-members:
   :no-inherited-members:

.. currentmodule:: mne.epochs

.. autosummary::
   :toctree: generated/

   add_channels_epochs
   average_movements
   combine_event_ids
   equalize_epoch_counts
   make_metadata
Source Space Data
=================

.. currentmodule:: mne

.. autosummary::
   :toctree: generated/

   BiHemiLabel
   Label
   MixedSourceEstimate
   MixedVectorSourceEstimate
   SourceEstimate
   VectorSourceEstimate
   VolSourceEstimate
   VolVectorSourceEstimate
   SourceMorph
   compute_source_morph
   extract_label_time_course
   grade_to_tris
   grade_to_vertices
   label.select_sources
   grow_labels
   label_sign_flip
   labels_to_stc
   morph_labels
   random_parcellation
   read_labels_from_annot
   read_dipole
   read_label
   read_source_estimate
   read_source_morph
   split_label
   stc_to_label
   stc_near_sensors
   transform_surface_to
   write_labels_to_annot
   write_label
   source_space.compute_distance_to_sensors

Realtime
========

Realtime functionality has moved to the standalone module :mod:`mne_realtime`.
File I/O
========

.. currentmodule:: mne

.. autosummary::
   :toctree: generated

   channel_type
   channel_indices_by_type
   get_head_surf
   get_meg_helmet_surf
   get_volume_labels_from_aseg
   get_volume_labels_from_src
   parse_config
   read_labels_from_annot
   read_bem_solution
   read_bem_surfaces
   read_cov
   read_dipole
   read_epochs
   read_epochs_kit
   read_epochs_eeglab
   read_epochs_fieldtrip
   read_events
   read_evokeds
   read_evoked_fieldtrip
   read_evokeds_mff
   read_freesurfer_lut
   read_forward_solution
   read_label
   read_morph_map
   read_proj
   read_reject_parameters
   read_source_estimate
   read_source_spaces
   read_surface
   read_trans
   read_tri
   write_labels_to_annot
   write_bem_solution
   write_bem_surfaces
   write_head_bem
   write_cov
   write_events
   write_evokeds
   write_forward_solution
   write_label
   write_proj
   write_source_spaces
   write_surface
   write_trans
   what
   io.read_info
   io.show_fiff

Base class:

.. autosummary::
   :toctree: generated

   BaseEpochs

Forward Modeling
================

.. currentmodule:: mne

.. autosummary::
   :toctree: generated/

   Forward
   SourceSpaces
   add_source_space_distances
   apply_forward
   apply_forward_raw
   average_forward_solutions
   convert_forward_solution
   decimate_surface
   dig_mri_distances
   forward.compute_depth_prior
   forward.compute_orient_prior
   forward.restrict_forward_to_label
   forward.restrict_forward_to_stc
   make_bem_model
   make_bem_solution
   make_forward_dipole
   make_forward_solution
   make_field_map
   make_sphere_model
   morph_source_spaces
   read_bem_surfaces
   read_forward_solution
   read_trans
   read_source_spaces
   read_surface
   sensitivity_map
   setup_source_space
   setup_volume_source_space
   surface.complete_surface_info
   surface.read_curvature
   use_coil_def
   write_bem_surfaces
   write_trans

:py:mod:`mne.bem`:

.. automodule:: mne.bem
   :no-members:
   :no-inherited-members:

.. currentmodule:: mne.bem

.. autosummary::
   :toctree: generated/

   ConductorModel
   fit_sphere_to_headshape
   get_fitting_dig
   make_watershed_bem
   make_flash_bem
   make_scalp_surfaces
   convert_flash_mris
:orphan:

.. _general_bibliography:

General bibliography
====================

The references below are arranged alphabetically by first author.

.. bibliography:: ./references.bib
   :all:
   :list: enumerated

Covariance computation
======================

.. currentmodule:: mne

.. autosummary::
   :toctree: generated/

   Covariance
   compute_covariance
   compute_raw_covariance
   cov.compute_whitener
   cov.prepare_noise_cov
   cov.regularize
   compute_rank
   make_ad_hoc_cov
   read_cov
   write_cov

Preprocessing
=============

Projections:

.. currentmodule:: mne

.. autosummary::
   :toctree: generated/

   Projection
   compute_proj_epochs
   compute_proj_evoked
   compute_proj_raw
   read_proj
   write_proj

:py:mod:`mne.channels`:

.. currentmodule:: mne.channels

.. automodule:: mne.channels
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/

   Layout
   DigMontage
   compute_native_head_t
   fix_mag_coil_types
   read_polhemus_fastscan
   get_builtin_montages
   make_dig_montage
   read_dig_polhemus_isotrak
   read_dig_captrak
   read_dig_dat
   read_dig_egi
   read_dig_fif
   read_dig_hpts
   read_dig_localite
   make_standard_montage
   read_custom_montage
   compute_dev_head_t
   read_layout
   find_layout
   make_eeg_layout
   make_grid_layout
   find_ch_adjacency
   read_ch_adjacency
   equalize_channels
   rename_channels
   generate_2d_layout
   make_1020_channel_selections
   combine_channels

:py:mod:`mne.preprocessing`:

.. currentmodule:: mne.preprocessing

.. automodule:: mne.preprocessing
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/

   ICA
   Xdawn
   annotate_amplitude
   annotate_break
   annotate_movement
   annotate_muscle_zscore
   annotate_nan
   compute_average_dev_head_t
   compute_current_source_density
   compute_fine_calibration
   compute_maxwell_basis
   compute_proj_ecg
   compute_proj_eog
   cortical_signal_suppression
   create_ecg_epochs
   create_eog_epochs
   find_bad_channels_maxwell
   find_ecg_events
   find_eog_events
   fix_stim_artifact
   ica_find_ecg_events
   ica_find_eog_events
   infomax
   equalize_bads
   maxwell_filter
   oversampled_temporal_projection
   peak_finder
   read_ica
   realign_raw
   regress_artifact
   corrmap
   read_ica_eeglab
   read_fine_calibration
   write_fine_calibration

:py:mod:`mne.preprocessing.nirs`:

.. currentmodule:: mne.preprocessing.nirs

.. automodule:: mne.preprocessing.nirs
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/

   optical_density
   beer_lambert_law
   source_detector_distances
   short_channels
   scalp_coupling_index
   temporal_derivative_distribution_repair

:py:mod:`mne.preprocessing.ieeg`:

.. currentmodule:: mne.preprocessing.ieeg

.. automodule:: mne.preprocessing.ieeg
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/

   project_sensors_onto_brain

EEG referencing:

.. currentmodule:: mne

.. autosummary::
   :toctree: generated/

   add_reference_channels
   set_bipolar_reference
   set_eeg_reference

:py:mod:`mne.filter`:

.. currentmodule:: mne.filter

.. automodule:: mne.filter
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/

   construct_iir_filter
   create_filter
   estimate_ringing_samples
   filter_data
   notch_filter
   resample

:py:mod:`mne.chpi`

.. currentmodule:: mne.chpi

.. automodule:: mne.chpi
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/

   compute_chpi_amplitudes
   compute_chpi_snr
   compute_chpi_locs
   compute_head_pos
   extract_chpi_locs_ctf
   extract_chpi_locs_kit
   filter_chpi
   get_chpi_info
   head_pos_to_trans_rot_t
   read_head_pos
   write_head_pos

:py:mod:`mne.transforms`

.. currentmodule:: mne.transforms

.. automodule:: mne.transforms
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/

   Transform
   quat_to_rot
   rot_to_quat
   read_ras_mni_t
:orphan:

.. _mne_cpp:

======================
MNE with CPP
======================

MNE-CPP is a cross-platform application and algorithm C++ framework
for MEG/EEG data acquisition, analysis and visualization. It
provides a modular structure with sub-libraries. The MNE-CPP API
can be integrated into other stand-alone projects to, e.g. provide
full I/O support for the FIF-file format or files generated by the
MNE and Freesurfer suite. MNE-CPP’s 3D visualization library is based
on the Qt3D module, which provides tools for online data displaying
with OpenGL.

MNE-CPP ships with built-in stand-alone applications, some of which
are closely connected to well-known MNE-C applications. MNE Browse can
be used to inspect and process pre-recorded data. Among others, dipole
fitting and the computation of forward solutions have been ported from
the MNE-C library, including the same command line interfaces. With MNE
Scan the MNE-CPP project provides an application for acquiring and
processing MEG/EEG data in real-time. Supported MEG devices include
the Elekta Neuromag VectorView and BabyMEG system. Several EEG amplifiers
(TMSI Refa, BrainAmp, ANT eegosports, gUSBamp) are supported as well.

For further information please visit the MNE-CPP project pages:

  * `Project Page <https://www.mne-cpp.org/>`_
  * `GitHub Sources <https://github.com/mne-tools/mne-cpp/>`_

.. raw:: html

    <div><script type="text/javascript" src="http://www.openhub.net/p/687714/widgets/project_basic_stats.js"></script></div>
.. _cited:

Papers citing MNE-Python
========================

Estimates provided by Google Scholar as of 02 November 2021:

- `MNE (1100) <https://scholar.google.com/scholar?cites=12188330066413208874&as_ylo=2014>`_
- `MNE-Python (1060) <https://scholar.google.com/scholar?cites=1521584321377182930&as_ylo=2013>`_

Time-Frequency
==============

:py:mod:`mne.time_frequency`:

.. automodule:: mne.time_frequency
   :no-members:
   :no-inherited-members:

.. currentmodule:: mne.time_frequency

.. autosummary::
   :toctree: generated/

   AverageTFR
   EpochsTFR
   CrossSpectralDensity

Functions that operate on mne-python objects:

.. autosummary::
   :toctree: generated/

   csd_fourier
   csd_multitaper
   csd_morlet
   pick_channels_csd
   read_csd
   fit_iir_model_raw
   psd_welch
   psd_multitaper
   tfr_morlet
   tfr_multitaper
   tfr_stockwell
   read_tfrs
   write_tfrs

Functions that operate on ``np.ndarray`` objects:

.. autosummary::
   :toctree: generated/

   csd_array_fourier
   csd_array_multitaper
   csd_array_morlet
   dpss_windows
   morlet
   stft
   istft
   stftfreq
   psd_array_multitaper
   psd_array_welch
   tfr_array_morlet
   tfr_array_multitaper
   tfr_array_stockwell


:py:mod:`mne.time_frequency.tfr`:

.. automodule:: mne.time_frequency.tfr
   :no-members:
   :no-inherited-members:

.. currentmodule:: mne.time_frequency.tfr

.. autosummary::
   :toctree: generated/

   cwt
   morlet

Datasets
========

.. currentmodule:: mne.datasets

:py:mod:`mne.datasets`:

.. automodule:: mne.datasets
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/

   fetch_dataset
   has_dataset
   brainstorm.bst_auditory.data_path
   brainstorm.bst_resting.data_path
   brainstorm.bst_raw.data_path
   eegbci.load_data
   eegbci.standardize
   fetch_aparc_sub_parcellation
   fetch_fsaverage
   fetch_hcp_mmp_parcellation
   fetch_infant_template
   fetch_phantom
   fnirs_motor.data_path
   hf_sef.data_path
   kiloword.data_path
   limo.load_data
   misc.data_path
   mtrf.data_path
   multimodal.data_path
   opm.data_path
   sleep_physionet.age.fetch_data
   sleep_physionet.temazepam.fetch_data
   sample.data_path
   somato.data_path
   spm_face.data_path
   visual_92_categories.data_path
   phantom_4dbti.data_path
   refmeg_noise.data_path
   ssvep.data_path
   erp_core.data_path
   epilepsy_ecog.data_path:orphan:

.. _inside_martinos:

Martinos Center setup
---------------------

For people within the MGH/MIT/HMS Martinos Center, MNE is available on the network.

In a terminal do:

.. code-block:: console

    $ setenv PATH /usr/pubsw/packages/python/anaconda/bin:${PATH}

If you use Bash replace the previous instruction with:

.. code-block:: console

    $ export PATH=/usr/pubsw/packages/python/anaconda/bin:${PATH}

Then start the python interpreter with:

.. code-block:: console

    $ ipython

Then type::

    >>> import mne

If you get a new prompt with no error messages, you should be good to go.

We encourage all Martinos center Python users to subscribe to the
`Martinos Python mailing list`_.

.. _Martinos Python mailing list: https://mail.nmr.mgh.harvard.edu/mailman/listinfo/martinos-python
:orphan:

.. _funding:

Funding and other support
=========================

Development of MNE-Python has been supported by:

.. rst-class:: list-unstyled funders

- |nih| **National Institutes of Health:**
  `R01-EB009048 <https://reporter.nih.gov/project-details/9053482>`_,
  `R01-EB006385 <https://reporter.nih.gov/project-details/8105475>`_,
  `R01-HD040712 <https://reporter.nih.gov/project-details/8511739>`_,
  `R01-NS044319 <https://reporter.nih.gov/project-details/6924553>`_,
  `R01-NS037462 <https://reporter.nih.gov/project-details/9083237>`_,
  `R01-NS104585 <https://reporter.nih.gov/project-details/10175064>`_,
  `P41-EB015896 <https://reporter.nih.gov/project-details/9518908>`_,
  `P41-RR014075 <https://reporter.nih.gov/project-details/8098820>`_
- |nsf| **US National Science Foundation:**
  `0958669 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=0958669>`_,
  `1042134 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=1042134>`_
- |erc| **European Research Council:**
  `YStG-263584 <https://erc.easme-web.eu/?p=263584>`_,
  `YStG-676943 <https://erc.easme-web.eu/?p=676943>`_
- |doe| **US Department of Energy:** DE-FG02-99ER62764 (MIND)
- |anr| **Agence Nationale de la Recherche:**
  `14-NEUC-0002-01 <https://anr.fr/Project-ANR-14-NEUC-0002>`_,
  **IDEX** Paris-Saclay
  `11-IDEX-0003-02 <https://anr.fr/ProjetIA-11-IDEX-0003>`_
- |cds| **Paris-Saclay Center for Data Science:**
  `PARIS-SACLAY <http://www.datascience-paris-saclay.fr>`_
- |goo| **Google:**
  Summer of code (×7 years)
- |ama| **Amazon:**
  AWS Research Grants
- |czi| **Chan Zuckerberg Initiative:**
  `EOSS2`_,
  `EOSS4`_


.. _supporting-institutions:

Institutional partners
----------------------

Additionally, many universities or research institutions have supported their
employees’ contributions to MNE-Python as part of normal work duties. These
institutions include:

.. include:: _includes/institutional-partners.rst
   :start-after: institutional-partners-begin-content


.. |nih| image:: _static/funding/nih.png
.. |nsf| image:: _static/funding/nsf.png
.. |erc| image:: _static/funding/erc.svg
.. |doe| image:: _static/funding/doe.svg
.. |anr| image:: _static/funding/anr.svg
.. |cds| image:: _static/funding/cds.png
.. |goo| image:: _static/funding/google.svg
.. |ama| image:: _static/funding/amazon.svg
.. |czi| image:: _static/funding/czi.svg

.. include:: links.inc
.. _api_reference_statistics:

Statistics
==========

:py:mod:`mne.stats`:

.. automodule:: mne.stats
   :no-members:
   :no-inherited-members:

.. currentmodule:: mne.stats

Parametric statistics (see :mod:`scipy.stats` and :mod:`statsmodels` for more
options):

.. autosummary::
   :toctree: generated/

   ttest_1samp_no_p
   ttest_ind_no_p
   f_oneway
   f_mway_rm
   f_threshold_mway_rm
   linear_regression
   linear_regression_raw

Mass-univariate multiple comparison correction:

.. autosummary::
   :toctree: generated/

   bonferroni_correction
   fdr_correction

Non-parametric (clustering) resampling methods:

.. autosummary::
   :toctree: generated/

   combine_adjacency
   permutation_cluster_test
   permutation_cluster_1samp_test
   permutation_t_test
   spatio_temporal_cluster_test
   spatio_temporal_cluster_1samp_test
   summarize_clusters_stc
   bootstrap_confidence_interval

Compute ``adjacency`` matrices for cluster-level statistics:

.. currentmodule:: mne

.. autosummary::
   :toctree: generated/

   channels.find_ch_adjacency
   channels.read_ch_adjacency
   spatial_dist_adjacency
   spatial_src_adjacency
   spatial_tris_adjacency
   spatial_inter_hemi_adjacency
   spatio_temporal_src_adjacency
   spatio_temporal_tris_adjacency
   spatio_temporal_dist_adjacency

Logging and Configuration
=========================

.. currentmodule:: mne

.. autosummary::
   :toctree: generated/

   get_config_path
   get_config
   open_docs
   set_log_level
   set_log_file
   set_config
   set_cache_dir
   sys_info
   use_log_level
   verbose

:py:mod:`mne.utils`:

.. currentmodule:: mne.utils

.. automodule:: mne.utils
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/

   deprecated
   warn

:py:mod:`mne.cuda`:

.. currentmodule:: mne.cuda

.. automodule:: mne.cuda
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/

   get_cuda_memory
   init_cuda
   set_cuda_device

Simulation
==========

:py:mod:`mne.simulation`:

.. automodule:: mne.simulation
   :no-members:
   :no-inherited-members:

.. currentmodule:: mne.simulation

.. autosummary::
   :toctree: generated/

   add_chpi
   add_ecg
   add_eog
   add_noise
   simulate_evoked
   simulate_raw
   simulate_stc
   simulate_sparse_stc
   select_source_in_label
   SourceSimulator

Creating data objects from arrays
=================================

.. currentmodule:: mne

.. autosummary::
   :toctree: generated/

   EvokedArray
   EpochsArray
   io.RawArray
   create_info

Connectivity Estimation
=======================

As of 0.24, connectivity functionality has been moved to the separate package
:mod:`mne-connectivity:mne_connectivity`.

Exporting
================

:py:mod:`mne.export`:

.. automodule:: mne.export
   :no-members:
   :no-inherited-members:

.. currentmodule:: mne.export

.. autosummary::
   :toctree: generated/

   export_epochs
   export_evokeds
   export_evokeds_mff
   export_raw
Installing FreeSurfer
=====================

`FreeSurfer <fs-wiki_>`_ is software for analysis and visualization of MRI data.
In the MNE ecosystem, freesurfer is used to convert structural MRI scans into
models of the scalp, inner/outer skull, and cortical surfaces, which are used
to

1. model how changes in the electrical and magnetic field caused by neural
   activity propagate to the sensor locations (part of computing the "forward
   solution"), and

2. constrain the estimates of where brain activity may have occurred (in the
   "inverse imaging" step of source localization).

System requirements, setup instructions, and test scripts are provided on the
`FreeSurfer download page`_. Note that if you don't already have it, you will
need to install ``tcsh`` for FreeSurfer to work; ``tcsh`` is usually
pre-installed with macOS, and is available in the package repositories for
Linux-based systems (e.g., ``sudo apt install tcsh`` on Ubuntu-like systems).

.. LINKS

.. _fs-wiki: https://surfer.nmr.mgh.harvard.edu/fswiki/
.. _`FreeSurfer download page`: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall
.. _contributing:

Contributing guide
==================

.. include:: ../links.inc
.. highlight:: console

Thanks for taking the time to contribute! MNE-Python is an open-source project
sustained mostly by volunteer effort. We welcome contributions from anyone as
long as they abide by our `Code of Conduct`_.

There are lots of ways to contribute, such as:

.. rst-class:: icon-bullets

- |bug| Use the software, and when you find bugs, tell us about them! We can
  only fix the bugs we know about.
- |discourse| Answer questions on `our user forum`_.
- |comment| Tell us about parts of the documentation that you find confusing or
  unclear.
- |hand-sparkles| Tell us about things you wish MNE-Python could do, or things
  it can do but you wish they were easier.
- |universal-access| Improve the accessibility of our website.
- |fix-bug| Fix bugs.
- |remove-format| Fix mistakes in our function documentation strings.
- |magic| Implement new features.
- |pencil-alt| Improve existing tutorials or write new ones.
- |python| Contribute to one of the many Python packages that MNE-Python
  depends on.

To *report* bugs, *request* new features, or *ask about* confusing
documentation, it's usually best to open a new issue on `our user forum`_
first; you'll probably get help fastest that way, and it helps keep our GitHub
issue tracker focused on things that we *know* will require changes to our
software (as opposed to problems that can be fixed in the user's code). We may
ultimately ask you to open an issue on GitHub too, but starting on the forum
helps us keep things organized. For fastest results, be sure to include
information about your operating system and MNE-Python version, and (if
applicable) include a reproducible code sample that is as short as possible and
ideally uses one of :ref:`our example datasets <datasets>`.

If you want to *fix* bugs, *add* new features, or *improve* our
docstrings/tutorials/website, those kinds of contributions are made through
`our GitHub repository <MNE-Python GitHub_>`_. The rest of this page explains
how to set up your workflow to make contributing via GitHub as easy as
possible.

.. collapse:: |rocket| Want an example to work through?
   :class: success

   Feel free to just read through the rest of the page, but if you find it
   easier to "learn by doing", take a look at our
   `GitHub issues marked "easy"`_, pick one that looks interesting, and work
   through it while reading this guide!

.. _`GitHub issues marked "easy"`: https://github.com/mne-tools/mne-python/issues?q=is%3Aissue+is%3Aopen+label%3AEASY


Overview of contribution process
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: Reminder: all contributors are expected to follow our
          `code of conduct`_.

Changes to MNE-Python are typically made by `forking`_ the MNE-Python
repository, making changes to your fork (usually by `cloning`_ it to your
personal computer, making the changes locally, and then `pushing`_ the local
changes up to your fork on GitHub), and finally creating a `pull request`_ to incorporate
your changes back into the shared "upstream" version of the codebase.

In general you'll be working with three different copies of the MNE-Python
codebase: the official remote copy at https://github.com/mne-tools/mne-python
(usually called ``upstream``), your remote `fork`_ of the upstream repository
(similar URL, but with your username in place of ``mne-tools``, and usually
called ``origin``), and the local copy of the codebase on your computer. The
typical contribution process is to:

1. synchronize your local copy with ``upstream``

2. make changes to your local copy

3. `push`_ your changes to ``origin`` (your remote fork of the upstream)

4. submit a `pull request`_ from your fork into ``upstream``

The sections :ref:`basic-git` and :ref:`github-workflow` (below) describe this
process in more detail.


Setting up your local development environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Configuring git
~~~~~~~~~~~~~~~

.. sidebar:: Git GUI alternative

    `GitHub desktop`_ is a GUI alternative to command line git that some users
    appreciate; it is available for |windows| Windows and |apple| MacOS.

To get set up for contributing, make sure you have git installed on your local
computer:

- On Linux, the command ``sudo apt install git`` is usually sufficient; see the
  `official Linux instructions`_ for more options.

- On MacOS, download `the .dmg installer`_; Atlassian also provides `more
  detailed instructions and alternatives`_ such as using MacPorts or Homebrew.

- On Windows, download and install `git for Windows`_. With Git BASH it provides its own shell that
  includes many Linux-equivalent command line programs that are useful for development.

  *Windows 10 also offers the* `Windows subsystem for Linux`_ *that offers similar
  functionality to git BASH, but has not been widely tested by MNE-Python
  developers yet and may still pose problems with graphical output (e.g. building the documentation)*


Once git is installed, the only absolutely necessary configuration step is
identifying yourself and your contact info::

   $ git config --global user.name "Your Name"
   $ git config --global user.email you@yourdomain.example.com

Make sure that the same email address is associated with your GitHub account
and with your local git configuration. It is possible to associate multiple
emails with a GitHub account, so if you initially set them up with different
emails, you can add the local email to the GitHub account.

Sooner or later, git is going to ask you what text editor you want it to use
when writing commit messages, so you might as well configure that now too::

   $ git config --global core.editor emacs    # or vim, or nano, or subl, or...

There are many other ways to customize git's behavior; see `configuring git`_
for more information.


GNU Make
~~~~~~~~

We use `GNU Make`_ to organize commands or short scripts that are often needed
in development. These are stored in files with the name :file:`Makefile`.
MNE-Python has two Makefiles, one in the package's root directory (containing
mainly testing commands) and one in :file:`doc/` (containing recipes for
building our documentation pages in different ways).

To check if make is already installed type ::

   $ make

into a terminal and you should see ::

   make: *** No targets specified and no makefile found.  Stop.

If you don't see this or something similar:

.. sidebar:: If you get:

   *bash: conda: command not found*

   you need to add

   - :file:`{path_to_Anaconda}`
   - :file:`{path_to_Anaconda}\\Scripts`

   to Windows-PATH.

- For Linux/MacOS, get `GNU Make`_
- For Windows, you can install make for git BASH (which comes with `git for Windows`_):

  1. Download :file:`make-{newest.version}-without-guile-w32-bin.zip` from `ezwinports`_
  2. Extract zip-folder
  3. Copy the contents into :file:`{path_to_git}\\mingw64\\` (e.g. by merging the
     folders with the equivalent ones already inside)
  4. For the first time using git BASH, you need to run once (to be able to
     activate your mnedev-environment): ::

      $ conda init bash


Forking the MNE-Python repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you have git installed and configured, and before creating your local copy
of the codebase, go to the `MNE-Python GitHub`_ page and create a `fork`_ into
your GitHub user account.

.. image:: https://help.github.com/assets/images/help/repository/fork_button.jpg

This will create a copy of the MNE-Python codebase inside your GitHub user
account (this is called "your fork"). Changes you make to MNE-Python will
eventually get "pushed" to your fork, and will be incorporated into the
official version of MNE-Python (often called the "upstream version") through a
"pull request". This process will be described in detail below; a summary
of how that structure is set up is given here:

.. graphviz:: ../_static/diagrams/git_setup.dot
   :alt: Diagram of recommended git setup
   :align: left


Creating the virtual environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. sidebar:: Supported Python environments

    We strongly recommend the `Anaconda`_ or `Miniconda`_ environment managers
    for Python. Other setups are possible but are not officially supported by
    the MNE-Python development team; see discussion :ref:`here
    <other-py-distros>`. These instructions use  ``conda`` where possible;
    experts may replace those lines with some combination of ``git`` and
    ``pip``.

These instructions will set up a Python environment that is separated from your
system-level Python and any other managed Python environments on your computer.
This lets you switch between different versions of Python (MNE-Python requires
version 3.7 or higher) and also switch between the stable and development
versions of MNE-Python (so you can, for example, use the same computer to
analyze your data with the stable release, and also work with the latest
development version to fix bugs or add new features). Even if you've already
followed the :ref:`installation instructions <install-python>` for the stable
version of MNE-Python, you should now repeat that process to create a new,
separate environment for MNE-Python development (here we'll give it the name
``mnedev``)::

    $ curl --remote-name https://raw.githubusercontent.com/mne-tools/mne-python/main/environment.yml
    $ conda env create --file environment.yml --name mnedev
    $ conda activate mnedev

Now you'll have *two* MNE-Python environments: ``mne`` (or whatever custom
name you used when installing the stable version of MNE-Python) and ``mnedev``
that we just created. At this point ``mnedev`` also has the stable version of
MNE-Python (that's what the :file:`environment.yml` file installs), but we're
about to remove the stable version from ``mnedev`` and replace it with the
development version. To do that, we'll `clone`_ the MNE-Python repository from
your remote fork, and also connect the local copy to the ``upstream`` version
of the codebase, so you can stay up-to-date with changes from other
contributors. First, edit these two variables for your situation::

    $ GITHUB_USERNAME="insert_your_actual_GitHub_username_here"
    $ # pick where to put your local copy of MNE-Python development version:
    $ INSTALL_LOCATION="/opt"

.. note::
   On Windows, add ``set`` before the variable names (``set GITHUB_USERNAME=...``, etc.).

Then make a local clone of your remote fork (``origin``)::

    $ cd $INSTALL_LOCATION
    $ git clone https://github.com/$GITHUB_USERNAME/mne-python.git

Finally, set up a link between your local clone and the official repository
(``upstream``)::

    $ cd mne-python
    $ git remote add upstream https://github.com/mne-tools/mne-python.git
    $ git fetch --all

Now we'll remove the *stable* version of MNE-Python and replace it with the
*development* version (the clone we just created with git). Make sure you're in
the correct environment first (``conda activate mnedev``), and then do::

    $ cd $INSTALL_LOCATION/mne-python    # make sure we're in the right folder
    $ pip uninstall -y mne
    $ pip install -e .

The command ``pip install -e .`` installs a python module into the current
environment by creating a link to the source code directory (instead of copying
the code to pip's :file:`site_packages` directory, which is what normally
happens). This means that any edits you make to the MNE-Python source code will
be reflected the next time you open a Python interpreter and ``import mne``
(the ``-e`` flag of ``pip`` stands for an "editable" installation).

Finally, we'll add a few dependencies that are not needed for running
MNE-Python, but are needed for locally running our test suite::

    $ pip install -r requirements_testing.txt

And for building our documentation::

    $ pip install -r requirements_doc.txt
    $ conda install graphviz

.. note::
   On Windows, if you installed graphviz using the conda command above but still get an error like this::

      WARNING: dot command 'dot' cannot be run (needed for graphviz output), check the graphviz_dot setting

   try adding the graphviz folder to path::

      $ PATH=$CONDA_PREFIX\\Library\\bin\\graphviz:$PATH

To build documentation, you will also require `optipng`_:

- On Linux, use the command ``sudo apt install optipng``.

- On MacOS, optipng can be installed using Homebrew.

- On Windows, unzip :file:`optipng.exe` from the `optipng for Windows`_ archive
  into the :file:`doc/` folder. This step is optional for Windows users.

You can also choose to install some optional linters for reStructuredText::

    $ conda install -c conda-forge sphinx-autobuild doc8


.. _basic-git:

Basic git commands
~~~~~~~~~~~~~~~~~~

Learning to work with git can take a long time, because it is a complex and
powerful tool for managing versions of files across multiple users, each of
whom have multiple copies of the codebase. We've already seen in the setup
commands above a few of the basic git commands useful to an MNE-Python
developer:

- :samp:`git clone {<URL_OF_REMOTE_REPO>}` (make a local copy of a repository)

- :samp:`git remote add {<NICKNAME_OF_REMOTE>} {<URL_OF_REMOTE_REPO>}` (connect
  a local copy to an additional remote)

- ``git fetch --all`` (get the current state of connected remote repos)

Other commands that you will undoubtedly need relate to `branches`_. Branches
represent multiple copies of the codebase *within a local clone or remote
repo*. Branches are typically used to experiment with new features while still
keeping a clean, working copy of the original codebase that you can switch back
to at any time. The default branch of any repo is called ``main``, and
it is recommended that you reserve the ``main`` branch to be that clean copy
of the working ``upstream`` codebase. Therefore, if you want to add a new
feature, you should first synchronize your local ``main`` branch with the
``upstream`` repository, then create a new branch based off of ``main`` and
`check it out`_ so that any changes you make will exist on that new branch
(instead of on ``main``)::

    $ git checkout main            # switch to local main branch
    $ git fetch upstream             # get the current state of the remote upstream repo
    $ git merge upstream/main      # synchronize local main branch with remote upstream main branch
    $ git checkout -b new-feature-x  # create local branch "new-feature-x" and check it out

.. sidebar:: Alternative

    You can save some typing by using ``git pull upstream/main`` to replace
    the ``fetch`` and ``merge`` lines above.

Now that you're on a new branch, you can fix a bug or add a new feature, add a
test, update the documentation, etc. When you're done, it's time to organize
your changes into a series of `commits`_. Commits are like snapshots of the
repository — actually, more like a description of what has to change to get
from the most recent snapshot to the current snapshot.

Git knows that people often work on multiple changes in multiple files all at
once, but that ultimately they should separate those changes into sets of
related changes that are grouped together based on common goals (so that it's
easier for their colleagues to understand and review the changes). For example,
you might want to group all the code changes together in one commit, put new
unit tests in another commit, and changes to the documentation in a third
commit.  Git makes this possible with something called the `stage`_ (or
*staging area*). After you've made some changes to the codebase, you'll have
what git calls "unstaged changes", which will show up with the `status`_
command::

    $ git status    # see what state the local copy of the codebase is in

Those unstaged changes can be `added`_ to the stage one by one, by either
adding a whole file's worth of changes, or by adding only certain lines
interactively::

    $ git add mne/some_file.py      # add all the changes you made to this file
    $ git add mne/some_new_file.py  # add a completely new file in its entirety
    $ # enter interactive staging mode, to add only portions of a file:
    $ git add -p mne/viz/some_other_file.py

Once you've collected all the related changes together on the stage, the ``git
status`` command will now refer to them as "changes staged for commit". You can
commit them to the current branch with the `commit`_ command. If you just type
``git commit`` by itself, git will open the text editor you configured it to
use so that you can write a *commit message* — a short description of the
changes you've grouped together in this commit. You can bypass the text editor
by passing a commit message on the command line with the ``-m`` flag. For
example, if your first commit adds a new feature, your commit message might be::

    $ git commit -m 'ENH: adds feature X to the Epochs class'

Once you've made the commit, the stage is now empty, and you can repeat the
cycle, adding the unit tests and documentation changes::

    $ git add mne/tests/some_testing_file.py
    $ git commit -m 'add test of new feature X of the Epochs class'
    $ git add -p mne/some_file.py mne/viz/some_other_file.py
    $ git commit -m 'DOC: update Epochs and BaseEpochs docstrings'
    $ git add tutorials/new_tutorial_file.py
    $ git commit -m 'DOC: adds new tutorial about feature X'

When you're done, it's time to run the test suite to make sure your changes
haven't broken any existing functionality, and to make sure your new test
covers the lines of code you've added (see :ref:`run-tests` and
:ref:`build-docs`, below). Once everything looks good, it's time to push your
changes to your fork::

    $ # push local changes to remote branch origin/new-feature-x
    $ # (this will create the remote branch if it doesn't already exist)
    $ git push origin new-feature-x

Finally, go to the `MNE-Python GitHub`_ page, click on the pull requests tab,
click the "new pull request" button, and choose "compare across forks" to
select your new branch (``new-feature-x``) as the "head repository".  See the
GitHub help page on `creating a PR from a fork`_ for more information about
opening pull requests.

If any of the tests failed before you pushed your changes, try to fix them,
then add and commit the changes that fixed the tests, and push to your fork. If
you're stuck and can't figure out how to fix the tests, go ahead and push your
commits to your fork anyway and open a pull request (as described above), then
in the pull request you should describe how the tests are failing and ask for
advice about how to fix them.

To learn more about git, check out the `GitHub help`_ website, the `GitHub
Learning Lab`_ tutorial series, and the `pro git book`_.


.. _github-ssh:

Connecting to GitHub with SSH (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One easy way to speed up development is to reduce the number of times you have
to type your password. SSH (secure shell) allows authentication with pre-shared
key pairs. The private half of your key pair is kept secret on your computer,
while the public half of your key pair is added to your GitHub account; when
you connect to GitHub from your computer, the local git client checks the
remote (public) key against your local (private) key, and grants access your
account only if the keys fit. GitHub has `several help pages`_ that guide you
through the process.

Once you have set up GitHub to use SSH authentication, you should change the
addresses of your MNE-Python GitHub remotes, from ``https://`` addresses to
``git@`` addresses, so that git knows to connect via SSH instead of HTTPS. For
example::

    $ git remote -v  # show existing remote addresses
    $ git remote set-url origin git@github.com:$GITHUB_USERNAME/mne-python.git
    $ git remote set-url upstream git@github.com:mne-tools/mne-python.git


MNE-Python coding conventions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

General requirements
~~~~~~~~~~~~~~~~~~~~

All new functionality must have test coverage
---------------------------------------------

For example, a new `mne.Evoked` method in :file:`mne/evoked.py` should
have a corresponding test in :file:`mne/tests/test_evoked.py`.


All new functionality must be documented
----------------------------------------

This includes thorough docstring descriptions for all public API changes, as
well as how-to examples or longer tutorials for major contributions. Docstrings
for private functions may be more sparse, but should usually not be omitted.


Avoid API changes when possible
-------------------------------

Changes to the public API (e.g., class/function/method names and signatures)
should not be made lightly, as they can break existing user scripts. Changes to
the API require a deprecation cycle (with warnings) so that users have time to
adapt their code before API changes become default behavior. See :ref:`the
deprecation section <deprecating>` and `mne.utils.deprecated` for
instructions. Bug fixes (when something isn't doing what it says it will do) do
not require a deprecation cycle.

Note that any new API elements should be added to the main reference;
classes, functions, methods, and attributes cannot be cross-referenced unless
they are included in the :ref:`api_reference`
(:file:`doc/python_reference.rst`).


.. _deprecating:

Deprecate with a decorator or a warning
---------------------------------------

MNE-Python has a :func:`~mne.utils.deprecated` decorator for classes and
functions that will be removed in a future version:

.. code-block:: python

    from mne.utils import deprecated

    @deprecated('my_function is deprecated and will be removed in 0.XX; please '
                'use my_new_function instead.')
    def my_function():
       return 'foo'

If you need to deprecate a parameter, use :func:`mne.utils.warn`. For example,
to rename a parameter from ``old_param`` to ``new_param`` you can do something
like this:

.. code-block:: python

    from mne.utils import warn

    def my_other_function(new_param=None, old_param=None):
        if old_param is not None:
            depr_message = ('old_param is deprecated and will be replaced by '
                            'new_param in 0.XX.')
            if new_param is None:
                new_param = old_param
                warn(depr_message, DeprecationWarning)
            else:
                warn(depr_message + ' Since you passed values for both '
                     'old_param and new_param, old_param will be ignored.',
                     DeprecationWarning)
        # Do whatever you have to do with new_param
        return 'foo'

When deprecating, you should also add corresponding test(s) to the relevant
test file(s), to make sure that the warning(s) are being issued in the
conditions you expect:

.. code-block:: python

    # test deprecation warning for function
    with pytest.warns(DeprecationWarning, match='my_function is deprecated'):
        my_function()

    # test deprecation warning for parameter
    with pytest.warns(DeprecationWarning, match='values for both old_param'):
        my_other_function(new_param=1, old_param=2)
    with pytest.warns(DeprecationWarning, match='old_param is deprecated and'):
        my_other_function(old_param=2)

You should also search the codebase for any cases where the deprecated function
or parameter are being used internally, and update them immediately (don't wait
to the *end* of the deprecation cycle to do this). Later, at the end of the
deprecation period when the stated release is being prepared:

- delete the deprecated functions
- remove the deprecated parameters (along with the conditional branches of
  ``my_other_function`` that handle the presence of ``old_param``)
- remove the deprecation tests
- double-check for any other tests that relied on the deprecated test or
  parameter, and (if found) update them to use the new function / parameter.


Describe your changes in the changelog
--------------------------------------

Include in your changeset a brief description of the change in the
:ref:`changelog <whats_new>` (:file:`doc/changes/latest.inc`; this can be
skipped for very minor changes like correcting typos in the documentation).

There are different sections of the changelog for each release, and separate
**subsections for bugfixes, new features, and changes to the public API.**
Please be sure to add your entry to the appropriate subsection.

The styling and positioning of the entry depends on whether you are a
first-time contributor or have been mentioned in the changelog before.

First-time contributors
"""""""""""""""""""""""

Welcome to MNE-Python! We're very happy to have you here. 🤗 And to ensure you
get proper credit for your work, please add a changelog entry with the
following pattern **at the top** of the respective subsection (bugfix,
new feature etc.):

.. code-block:: rst


  Bug
  ---

  .. |Your Name| replace:: **Your Name**

  - Short description of the changes (:gh:`0000` **by new contributor** |Your Name|_)

  - ...

where ``0000`` must be replaced with the respective GitHub pull request (PR)
number.

It is usually best to wait to add a line to the changelog until your PR is
finalized, to avoid merge conflicts (since the changelog is updated with
almost every PR).

Lastly, make sure that your name is included in the list of authors in
:file:`doc/changes/names.inc`, otherwise the documentation build will fail.
To add an author name, append a line with the following pattern (note
how the syntax is different from that used in the changelog):

.. code-block:: rst

  .. _Your Name: https://www.your-website.com/

Many contributors opt to link to their GitHub profile that way. Have a look
at the existing entries in the file to get some inspiration.

Recurring contributors
""""""""""""""""""""""

The changelog entry should follow the following patterns:

.. code-block:: rst

    - Short description of the changes from one contributor (:gh:`0000` by `Contributor Name`_)
    - Short description of the changes from several contributors (:gh:`0000` by `Contributor Name`_, `Second Contributor`_, and `Third Contributor`_)

where ``0000`` must be replaced with the respective GitHub pull request (PR)
number. Mind the Oxford comma in the case of multiple contributors.

Sometimes, changes that shall appear as a single changelog entry are spread out
across multiple PRs. In this case, name all relevant PRs, separated by
commas:

.. code-block:: rst

    - Short description of the changes from one contributor in multiple PRs (:gh:`0000`, :gh:`1111` by `Contributor Name`_)
    - Short description of the changes from several contributors in multiple PRs (:gh:`0000`, :gh:`1111` by `Contributor Name`_, `Second Contributor`_, and `Third Contributor`_)

Test locally before opening pull requests (PRs)
-----------------------------------------------

MNE-Python uses `continuous integration`_ (CI) to ensure code quality and
test across multiple installation targets. However, the CIs are often slower
than testing locally, especially when other contributors also have open PRs
(which is basically always the case). Therefore, do not rely on the CIs to
catch bugs and style errors for you; :ref:`run the tests locally <run-tests>`
instead before opening a new PR and before each time you push additional
changes to an already-open PR.


Make tests fast and thorough
----------------------------

Whenever possible, use the testing dataset rather than one of the sample
datasets when writing tests; it includes small versions of most MNE-Python
objects (e.g., `~mne.io.Raw` objects with short durations and few
channels). You can also check which lines are missed by the tests, then modify
existing tests (or write new ones) to target the missed lines. Here's an
example that reports which lines within ``mne.viz`` are missed when running
:file:`test_evoked.py` and :file:`test_topo.py`::

    $ pytest --cov=mne.viz --cov-report=term-missing mne/viz/tests/test_evoked.py mne/viz/tests/test_topo.py

You can also use ``pytest --durations=5`` to ensure new or modified tests will
not slow down the test suite too much.


Code style
~~~~~~~~~~

Adhere to standard Python style guidelines
------------------------------------------

All contributions to MNE-Python are checked against style guidelines described
in `PEP 8`_. We also check for common coding errors (such as variables that are
defined but never used). We allow very few exceptions to these guidelines, and
use tools such as pep8_, pyflakes_, and flake8_ to check code style
automatically. From the :file:`mne-python` root directory, you can check for
style violations by running::

    $ make flake

in the shell. Several text editors or IDEs also have Python style checking,
which can highlight style errors while you code (and train you to make those
errors less frequently). This functionality is built-in to the Spyder_ IDE, but
most editors have plug-ins that provide similar functionality. Search for
:samp:`python linter <name of your favorite editor>` to learn more.


Use consistent variable naming
------------------------------

Classes should be named using ``CamelCase``. Functions and instances/variables
should use ``snake_case`` (``n_samples`` rather than ``nsamples``). Avoid
single-character variable names, unless inside a :term:`comprehension <list
comprehension>` or :ref:`generator <tut-generators>`.


We (mostly) follow NumPy style for docstrings
---------------------------------------------

In most cases you can look at existing MNE-Python docstrings to figure out how
yours should be formatted. If you can't find a relevant example, consult the
`Numpy docstring style guidelines`_ for examples of more complicated formatting
such as embedding example code, citing references, or including rendered
mathematics.  Note that we diverge from the NumPy docstring standard in a few
ways:

1. We use a module called ``sphinxcontrib-bibtex`` to render citations. Search
   our source code (``git grep footcite`` and ``git grep footbibliography``) to
   see examples of how to add in-text citations and formatted references to
   your docstrings, examples, or tutorials. The structured bibliographic data
   lives in :file:`doc/references.bib`; please follow the existing key scheme
   when adding new references (e.g., ``Singleauthor2019``,
   ``AuthoroneAuthortwo2020``, ``FirstauthorEtAl2021a``,
   ``FirstauthorEtAl2021b``).
2. We don't explicitly say "optional" for optional keyword parameters (because
   it's clear from the function or method signature which parameters have
   default values).
3. For parameters that may take multiple types, we use pipe characters instead
   of the word "or", like this: ``param_name : str | None``.
4. We don't include a ``Raises`` or ``Warns`` section describing
   errors/warnings that might occur.


Private function/method docstrings may be brief for simple functions/methods,
but complete docstrings are appropriate when private functions/methods are
relatively complex. To run some basic tests on documentation, you can use::

    $ pytest mne/tests/test_docstring_parameters.py
    $ make docstyle


Cross-reference everywhere
--------------------------

Both the docstrings and dedicated documentation pages (tutorials, how-to
examples, discussions, and glossary) should include cross-references to any
mentioned module, class, function, method, attribute, or documentation page.
There are sphinx roles for all of these (``:mod:``, ``:class:``,
``:func:``, ``:meth:``, ``:attr:``, ``:doc:``) as well as a generic
cross-reference directive (``:ref:``) for linking to specific sections of a
documentation page.

.. warning::

    Some API elements have multiple exposure points (for example,
    ``mne.set_config`` and ``mne.utils.set_config``). For cross-references to
    work, they must match an entry in :file:`doc/python_reference.rst` (thus
    ``:func:`mne.set_config``` will work but ``:func:`mne.utils.set_config```
    will not).

MNE-Python also uses Intersphinx_, so you can (and should)
cross-reference to Python built-in classes and functions as well as API
elements in :mod:`NumPy <numpy>`, :mod:`SciPy <scipy>`, etc. See the Sphinx
configuration file (:file:`doc/conf.py`) for the list of Intersphinx projects
we link to. Their inventories can be examined using a tool like `sphobjinv`_ or
dumped to file with commands like::

    $ python -m sphinx.ext.intersphinx https://docs.python.org/3/objects.inv > python.txt

Note that anything surrounded by single backticks that is *not* preceded by one
of the API roles (``:class:``, ``:func:``, etc) will be assumed to be
in the MNE-Python namespace. This can save some typing especially in
tutorials; instead of ``see :func:`mne.io.Raw.plot_psd` for details`` you can
instead type ``see `mne.io.Raw.plot_psd` for details``.


Other style guidance
--------------------

- Use single quotes whenever possible.

- Prefer :ref:`generators <tut-generators>` or
  :term:`comprehensions <list comprehension>` over :func:`filter`, :func:`map`
  and other functional idioms.

- Use explicit functional constructors for builtin containers to improve
  readability (e.g., :ref:`list() <func-list>`, :ref:`dict() <func-dict>`,
  :ref:`set() <func-set>`).

- Avoid nested functions or class methods if possible — use private functions
  instead.

- Avoid ``*args`` and ``**kwargs`` in function/method signatures.


Code organization
~~~~~~~~~~~~~~~~~

Importing
---------

Import modules in this order, preferably alphabetized within each subsection:

1. Python built-in (``copy``, ``functools``, ``os``, etc.)
2. NumPy (``numpy as np``) and, in test files, pytest (``pytest``)
3. MNE-Python imports (e.g., ``from .pick import pick_types``)

When importing from other parts of MNE-Python, use relative imports in the main
codebase and absolute imports in tests, tutorials, and how-to examples. Imports
for ``matplotlib``, ``scipy``, and optional modules (``sklearn``, ``pandas``,
etc.) should be nested (i.e., within a function or method, not at the top of a
file). This helps reduce import time and limit hard requirements for using MNE.


Return types
------------

Methods should modify inplace and return ``self``, functions should return
copies (where applicable). Docstrings should always give an informative name
for the return value, even if the function or method's return value is never
stored under that name in the code.


Visualization
-------------

Visualization capabilities should be made available in both function and method
forms. Add public visualization functions to the :mod:`mne.viz` submodule, and
call those functions from the corresponding object methods. For example, the
method :meth:`mne.Epochs.plot` internally calls the function
:func:`mne.viz.plot_epochs`.

All visualization functions must accept a boolean ``show`` parameter and
typically return a :class:`matplotlib.figure.Figure` (or a list of
:class:`~matplotlib.figure.Figure` objects). 3D visualization functions return
a :class:`mne.viz.Figure3D`, :class:`mne.viz.Brain`, or other return type
as appropriate.

Visualization functions should default to the colormap ``RdBu_r`` for signed
data with a meaningful middle (zero-point) and ``Reds`` otherwise. This applies
to both visualization functions and tutorials/examples.


.. _run_tests:

Running the test suite
~~~~~~~~~~~~~~~~~~~~~~

.. sidebar:: pytest flags

    The ``-x`` flag exits the pytest run when any test fails; this can speed
    up debugging when running all tests in a file or module.

    The ``--pdb`` flag will automatically start the python debugger upon test
    failure.

The full test suite can be run by calling ``make test`` from the
``mne-python`` root folder. Testing the entire module can be quite
slow, however, so to run individual tests while working on a new feature, you
can run the following line::

    $ pytest mne/tests/test_evoked.py::test_io_evoked --verbose

Or alternatively::

    $ pytest mne/tests/test_evoked.py -k test_io_evoked --verbose

Make sure you have the testing dataset, which you can get by running this in
a Python interpreter:

.. code-block:: python

    >>> mne.datasets.testing.data_path(verbose=True)  # doctest: +SKIP


.. _build-docs:

Building the documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Our documentation (including docstrings in code files) is in
reStructuredText_ format and is built using Sphinx_ and `Sphinx-Gallery`_.
The easiest way to ensure that your contributions to the documentation are
properly formatted is to follow the style guidelines on this page, imitate
existing documentation examples, refer to the Sphinx and Sphinx-Gallery
reference materials when unsure how to format your contributions, and build the
docs locally to confirm that everything looks correct before submitting the
changes in a pull request.

You can build the documentation locally using `GNU Make`_ with
:file:`doc/Makefile`. From within the :file:`doc` directory, you can test
formatting and linking by running::

    $ make html_dev-noplot

This will build the documentation *except* it will format (but not execute) the
tutorial and example files. If you have created or modified an example or
tutorial, you should instead run
:samp:`PATTERN={<REGEX_TO_SELECT_MY_TUTORIAL>} make html_dev-pattern` to render
all the documentation and additionally execute just your example or tutorial
(so you can make sure it runs successfully and generates the output / figures
you expect).

.. note::
   On Windows, to use the pattern approach, use the following two lines:

   .. code-block:: python

      set PATTERN={<REGEX_TO_SELECT_MY_TUTORIAL>}
      make html_dev-pattern

After either of these commands completes, ``make show`` will open the
locally-rendered documentation site in your browser. Additional ``make``
recipes are available; run ``make help`` from the :file:`doc` directory or
consult the `Sphinx-Gallery`_ documentation for additional details.


Modifying command-line tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MNE-Python provides support for a limited set of :ref:`python_commands`.
These are typically used with a call like::

    $ mne browse_raw ~/mne_data/MNE-sample-data/MEG/sample/sample_audvis_raw.fif

These are generally available for convenience, and can be useful for quick
debugging (in this case, for `mne.io.Raw.plot`).

If a given command-line function fails, they can also be executed as part of
the ``mne`` module with ``python -m``. For example::

    $ python -i -m mne browse_raw ...

Because this was launched with ``python -i``, once the script completes
it will drop to a Python terminal. This is useful when there are errors,
because then you can drop into a :func:`post-mortem debugger <python:pdb.pm>`:

.. code-block:: python

    >>> import pdb; pdb.pm()  # doctest:+SKIP


.. _`github-workflow`:

GitHub workflow
~~~~~~~~~~~~~~~

Nearly everyone in the community of MNE-Python contributors and maintainers is
a working scientist, engineer, or student who contributes to MNE-Python in
their spare time. For that reason, a set of best practices have been adopted to
streamline the collaboration and review process. Most of these practices are
common to many open-source software projects, so learning to follow them while
working on MNE-Python will bear fruit when you contribute to other projects
down the road. Here are the guidelines:

- Search the `MNE-Python issues page`_ (both open and closed issues) in case
  someone else has already started work on the same bugfix or feature. If you
  don't find anything, `open a new issue`_ to discuss changes with maintainers
  before starting work on your proposed changes.

- Implement only one new feature or bugfix per pull request (PR). Occasionally
  it may make sense to fix a few related bugs at once, but this makes PRs
  harder to review and test, so check with MNE-Python maintainers first before
  doing this. Avoid purely cosmetic changes to the code; they make PRs harder
  to review.

- It is usually better to make PRs *from* branches other than your main
  branch, so that you can use your main branch to easily get back to a
  working state of the code if needed (e.g., if you're working on multiple
  changes at once, or need to pull in recent changes from someone else to get
  your new feature to work properly).

- In most cases you should make PRs *into* the upstream's main branch, unless
  you are specifically asked by a maintainer to PR into another branch (e.g.,
  for backports or maintenance bugfixes to the current stable version).

- Don't forget to include in your PR a brief description of the change in the
  :ref:`changelog <whats_new>` (:file:`doc/whats_new.rst`).

- Our community uses the following commit tags and conventions:

  - Work-in-progress PRs should be created as `draft PRs`_ and the PR title
    should begin with ``WIP``.

  - When you believe a PR is ready to be reviewed and merged, `convert it
    from a draft PR to a normal PR`_, change its title to begin with ``MRG``,
    and add a comment to the PR asking for reviews (changing the title does not
    automatically notify maintainers).

  - PRs that only affect documentation should additionally be labelled
    ``DOC``, bugfixes should be labelled ``FIX``, and new features should be
    labelled ``ENH`` (for "enhancement"). ``STY`` is used for style changes
    (i.e., improving docstring consistency or formatting without changing its
    content).

  - the following commit tags are used to interact with our
    `continuous integration`_ (CI) providers. Use them judiciously; *do not
    skip tests simply because they are failing*:

    - ``[skip circle]`` Skip `CircleCI`_, which tests successful building of
      our documentation.

    - ``[skip actions]`` Skip our `GitHub Actions`_, which test installation
      and execution on Linux and macOS systems.

    - ``[skip azp]`` Skip `azure`_ which tests installation and execution on
      Windows systems.

    - ``[ci skip]`` is an alias for ``[skip actions][skip azp][skip circle]``.
      Notice that ``[skip ci]`` is not a valid tag.

    - ``[circle full]`` triggers a "full" documentation build, i.e., all code
      in tutorials and how-to examples will be *executed* (instead of just
      nicely formatted) and the resulting output and figures will be rendered
      as part of the tutorial/example.

`This sample pull request`_ exemplifies many of the conventions listed above:
it addresses only one problem; it started with an issue to discuss the problem
and some possible solutions; it is a PR from the user's non-main branch into
the upstream main branch; it separates different kinds of changes into
separate commits and uses labels like ``DOC``, ``FIX``, and ``STY`` to make it
easier for maintainers to review the changeset; etc. If you are new to GitHub
it can serve as a useful example of what to expect from the PR review process.


.. MNE

.. _MNE-Python GitHub: https://github.com/mne-tools/mne-python
.. _MNE-Python issues page: https://github.com/mne-tools/mne-python/issues
.. _open a new issue: https://github.com/mne-tools/mne-python/issues/new/choose
.. _This sample pull request: https://github.com/mne-tools/mne-python/pull/6230
.. _our user forum: https://mne.discourse.group

.. git installation

.. _the .dmg installer: https://git-scm.com/download/mac
.. _git for Windows: https://gitforwindows.org/
.. _official Linux instructions: https://git-scm.com/download/linux
.. _more detailed instructions and alternatives: https://www.atlassian.com/git/tutorials/install-git
.. _Windows subsystem for Linux: https://docs.microsoft.com/en-us/windows/wsl/about
.. _GitHub desktop: https://desktop.github.com/
.. _GNU Make: https://www.gnu.org/software/make/
.. _ezwinports: https://sourceforge.net/projects/ezwinports/files/

.. github help pages

.. _GitHub Help: https://help.github.com
.. _GitHub learning lab: https://lab.github.com/
.. _fork: https://help.github.com/en/articles/fork-a-repo
.. _clone: https://help.github.com/en/articles/cloning-a-repository
.. _push: https://help.github.com/en/articles/pushing-to-a-remote
.. _forking: https://help.github.com/en/articles/fork-a-repo
.. _cloning: https://help.github.com/en/articles/cloning-a-repository
.. _pushing: https://help.github.com/en/articles/pushing-to-a-remote
.. _branches: https://help.github.com/en/articles/about-branches
.. _several help pages: https://help.github.com/en/articles/connecting-to-github-with-ssh
.. _draft PRs: https://help.github.com/en/articles/about-pull-requests#draft-pull-requests
.. _convert it from a draft PR to a normal PR: https://help.github.com/en/articles/changing-the-stage-of-a-pull-request
.. _pull request: https://help.github.com/en/articles/creating-a-pull-request-from-a-fork
.. _creating a PR from a fork: https://help.github.com/en/articles/creating-a-pull-request-from-a-fork

.. git docs

.. _check it out: https://git-scm.com/docs/git-checkout
.. _added: https://git-scm.com/docs/git-add
.. _commits: https://git-scm.com/docs/git-commit
.. _commit: https://git-scm.com/docs/git-commit
.. _status: https://git-scm.com/docs/git-status

.. git book

.. _stage: https://git-scm.com/book/en/v2/Git-Tools-Interactive-Staging
.. _configuring git: https://www.git-scm.com/book/en/v2/Customizing-Git-Git-Configuration

.. sphinx

.. _sphinx-gallery: https://sphinx-gallery.github.io
.. _reStructuredText: https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html
.. _intersphinx: https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html
.. _sphobjinv: https://sphobjinv.readthedocs.io/en/latest/

.. linting

.. _NumPy docstring style guidelines: https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt
.. _PEP 8: https://www.python.org/dev/peps/pep-0008/
.. _pyflakes: https://pypi.org/project/pyflakes
.. _Flake8: http://flake8.pycqa.org/

.. misc

.. _miniconda: https://conda.io/en/latest/miniconda.html
.. _Spyder: https://www.spyder-ide.org/
.. _continuous integration: https://en.wikipedia.org/wiki/Continuous_integration
.. _matplotlib: https://matplotlib.org/
.. _github actions: https://docs.github.com/en/free-pro-team@latest/actions/learn-github-actions
.. _azure: https://dev.azure.com/mne-tools/mne-python/_build/latest?definitionId=1&branchName=main
.. _CircleCI: https://circleci.com/gh/mne-tools/mne-python

.. optipng

.. _optipng: http://optipng.sourceforge.net/
.. _optipng for Windows: http://prdownloads.sourceforge.net/optipng/optipng-0.7.7-win32.zip?download
Updating MNE-Python
===================

If you want to update MNE-Python to a newer version, there are a few different
options, depending on how you originally installed it.

.. warning::

    Before performing package upgrade operations, check to make sure that the
    environment you wish to modify has been activated (and if not, call
    ``conda activate name_of_environment`` first).


Upgrading MNE-Python only
^^^^^^^^^^^^^^^^^^^^^^^^^

If you wish to update MNE-Python only and leave other packages in their current
state, you can usually safely do this with ``pip``, even if you originally
installed via conda. With the ``mne`` environment active, do:

.. code-block:: console

    $ pip install -U mne


Upgrading all packages
^^^^^^^^^^^^^^^^^^^^^^

Generally speaking, if you want to upgrade *your whole software stack*
including all the dependencies, the best approach is to re-create it as a new
virtual environment, because neither conda nor pip are fool-proof at making
sure all packages remain compatible with one another during upgrades.

Here we'll demonstrate renaming the old environment first, as a safety measure.
We'll assume that the existing environment is called ``mne`` and you want to
rename the old one so that the new, upgraded environment can be called ``mne``
instead. Unfortunately ``conda`` doesn't have a "rename" command so we'll first
clone the old one with a new name (``old_mne``), then delete the original, then
create the new, updated environment re-using the original name. In the first
step we'll also use conda in ``--offline`` mode so that it uses cached
copies of all the packages instead of re-downloading them.

.. code-block:: console

    $ conda create --name old_mne --clone mne --offline  # copy with new name,
    $ conda env remove --name mne --all                  # remove original,
    $ conda create --name mne --channel conda-forge mne  # replace with updated

.. note::

    If you installed extra packages into your old ``mne`` environment,
    you'll need to repeat that process after re-creating the updated
    environment. Comparing the output of ``conda list --name old_mne`` versus
    ``conda list --name mne`` will show you what is missing from the new
    environment. On Linux, you can automate that comparison like this:

    .. code-block:: console

        $ diff <(conda list -n mne | cut -d " " -f 1 | sort) <(conda list -n old_mne | cut -d " " -f 1 | sort) | grep "^>" | cut -d " " -f 2


.. _installing_main:

Upgrading to the development version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning:: In between releases, function and class APIs can change without
    warning.

Sometimes, new features or bugfixes become available that are important to your
research and you just can't wait for the next official release of MNE-Python to
start taking advantage of them. In such cases, you can use ``pip`` to install
the *development version* of MNE-Python:

.. code-block:: console

    $ pip install -U --no-deps https://github.com/mne-tools/mne-python/archive/main.zip
:orphan:

.. include:: ../links.inc

.. _install_mne_c:

Installing MNE-C
================

System requirements
^^^^^^^^^^^^^^^^^^^

MNE-C runs on macOS (version 10.5 "Leopard" or later) and Linux (kernel 2.6.9
or later). Both 32- and 64-bit operating systems are supported; a PowerPC
version for macOS can be provided upon request. At least 2 GB of memory is
required, 4 GB or more is recommended. The software requires at least 80 MB of
disk space. MATLAB is an optional dependency; the free `MATLAB runtime`_ is
sufficient. If MATLAB is not present, the utilities ``mne_convert_mne_data``,
``mne_epochs2mat``, ``mne_raw2mat``, and ``mne_simu`` will not work.

For boundary-element model (BEM) mesh generation, and for accessing the ``tkmedit``
program from ``mne_analyze``, MNE-C needs access to a
working installation of :doc:`FreeSurfer <freesurfer>`, including the
environment variables ``FREESURFER_HOME``, ``SUBJECTS_DIR``, and ``SUBJECT``.

.. admonition:: |apple| macOS
  :class: note

  For installation on macOS, you also need:

  - the `XCode developer tools`_.
  - an X Window System such as XQuartz_. Version 2.7.9 of XQuartz should work
    out of the box; the most current version (2.7.11, as of May 2019) may
    require these additional steps to work:

    .. code-block:: console

        $ cd /opt/X11/lib
        $ sudo cp libXt.6.dylib libXt.6.dylib.bak
        $ cd flat_namespace/
        $ sudo cp libXt.6.dylib ../.

  - the netpbm_ library. The recommended way to get netpbm is to install
    Homebrew_, and run ``brew install netpbm`` in the Terminal app.
    Alternatively, if you prefer to use MacPorts_, you can run
    ``sudo port install netpbm`` in the Terminal app.


Downloading and Installing MNE-C
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MNE-C is distributed as either a compressed tar archive (.tar.gz) or a macOS
disk image (.dmg). The `MNE-C download page`_ requires registration with a
valid email address.  The current stable version is 2.7.3; "nightly" builds of
the development version are also available on the download page.

To install from the compressed tar archive, change directory to the desired
install location, and unpack the software using ``tar``:

.. code-block:: console

    $ cd <path_to_desired_install_location>
    $ tar zxvf <path_to_archive_file>

To install from the macOS disk image, double-click the downloaded .dmg file. In
the window that opens, double-click the installer package file (.pkg) to launch
the installer, and follow its instructions. In newer versions of macOS, if you
see an error that the app is from an untrusted developer, you can override this
warning by opening it anyway from the Security & Privacy pane within the
computer's System Preferences.

.. _user_environment:

Configuring MNE-C
^^^^^^^^^^^^^^^^^

MNE-C requires two environment variables to be defined manually:

- ``MNE_ROOT`` should give the path to the folder where MNE-C is installed
- ``MATLAB_ROOT`` should give the path to your MATLAB binary (e.g.,
  ``/opt/MATLAB/R2018b`` or similar).  If you do not have MATLAB or the MATLAB
  runtime, leave ``MATLAB_ROOT`` undefined.

Other environment variables are defined by setup scripts provided with MNE-C.
You may either run the setup script each time you use MNE-C, or (recommended)
configure your shell to run it automatically each time you open a terminal. For
bash compatible shells, e.g., sh/bash/zsh, the script to source is
``$MNE_ROOT/bin/mne_setup_sh``.  For C shells, e.g., csh/tcsh, the script to
source is ``$MNE_ROOT/bin/mne_setup``.  If you don't know what shell you are
using, you can run the following command to find out:

.. code-block:: console

    $ echo $SHELL

To configure MNE-C automatically for ``bash`` or ``sh`` shells, add this to
your ``.bashrc``:

.. code-block:: sh

    export MNE_ROOT=<path_to_MNE>
    export MATLAB_ROOT=<path_to_MATLAB>
    source $MNE_ROOT/bin/mne_setup_sh

where ``<path_to_MNE>`` and ``<path_to_MATLAB>`` are replaced by the absolute
paths to MNE-C and MATLAB, respectively. If you don't have MATLAB, you should
still include the ``export MATLAB_ROOT=`` statement, but leave
``<path_to_MATLAB>`` blank.

To configure MNE-C automatically for ``zsh``, use the built-in ``emulate``
command in your ``.zshrc`` file:

.. code-block:: sh

    export MNE_ROOT=<path_to_MNE>
    export MATLAB_ROOT=<path_to_MATLAB>
    emulate sh -c 'source $MNE_ROOT/bin/mne_setup_sh'

To configure MNE-C automatically for ``csh`` or ``tcsh`` shells, the
corresponding commands in the ``.cshrc`` / ``.tcshrc`` file are:

.. code-block:: tcsh

    setenv MNE_ROOT <path_to_MNE>
    setenv MATLAB_ROOT <path_to_MATLAB>
    source $MNE_ROOT/bin/mne_setup

If you have done this correctly, the command ``ls $MNE_ROOT/bin/mne_setup_sh``
should succeed when run in a new terminal.

Testing MNE-C installation
^^^^^^^^^^^^^^^^^^^^^^^^^^

An easy way to verify whether your installation of MNE-C is working is to test
the OpenGL graphics performance:

.. code-block:: console

    $ $MNE_ROOT/bin/mne_opengl_test

This will render an inflated brain surface repeatedly, rotating it by 5 degrees
around the z-axis between redraws. The time spent for each full revolution is
printed to the terminal window where ``mne_opengl_test`` was invoked.  Switch
focus to that terminal window and use the interrupt key (usually control-c) to
halt the test.

The best graphics performance occurs when MNE-C renders to a local display on a
computer with hardware acceleration enabled. The ``mne_analyze`` GUI has a menu
item "On GLX..." in the Help menu; if the GLX dialog says "Direct rendering
context" then hardware acceleration is in use. If you are rendering to a local
display and see "Nondirect rendering context", it is recommended that you
enable hardware acceleration (consult a search engine or your local IT support
staff for assistance). If you are rendering to a remote display or using a VNC
connection, "Nondirect rendering context" is normal.

On the fastest graphics cards, the time per revolution in the
``mne_opengl_test`` is well below 1 second. If your time per revolution is
longer than 10 seconds, either the graphics hardware acceleration is not in
effect or you need a faster graphics adapter.

Troubleshooting MNE-C installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If MNE-C can't find ``libxp.so.6``, download libxp6 from debian_ or similar and
install it:

.. code-block:: console

    $ sudo dpkg -i libxp6_1.0.2-1ubuntu1_amd64.deb

If MNE-C can't find ``libgfortran.so.1``, you can probably safely link that
filename to the current version of libfortran that came with your system. On
a typical 64-bit Ubuntu-like system this would be accomplished by:

.. code-block:: console

    $ cd /usr/lib/x86_64-linux-gnu
    $ sudo ln -s libgfortran.so.1 $(find . -maxdepth 1 -type f -name libgfortran.so*)

If you encounter other errors installing MNE-C, please post a message to the
`MNE Forum`_.

.. links

.. _MNE-C download page: http://www.nmr.mgh.harvard.edu/martinos/userInfo/data/MNE_register/index.php
.. _MATLAB runtime: https://www.mathworks.com/products/compiler/matlab-runtime.html
.. _netpbm: http://netpbm.sourceforge.net/
.. _MacPorts: https://www.macports.org/
.. _Homebrew: https://brew.sh/
.. _XCode developer tools: https://developer.apple.com/xcode/
.. _xquartz: https://www.xquartz.org/
.. _debian: https://packages.debian.org/jessie/amd64/libxp6/download
.. include:: ../links.inc

.. _quick-start:

Quick start
===========

MNE-Python requires Python version |min_python_version| or higher. If you've
never worked with Python before, skip ahead to the last paragraph of this page.
For users already familiar with Python:

- If you only need MNE-Python's computational functions, only hard dependencies
  will be included when running:

  .. code-block:: console

      $ pip install mne

- If you plan to use MNE-Python's functions that use HDF5-based I/O (e.g.,
  :func:`mne.io.read_raw_eeglab`, :meth:`mne.SourceMorph.save`, etc.),
  you should run:

  .. code-block:: console

      $ pip install mne[hdf5]

  This will pull in additional dependencies pymatreader_, h5io_, and h5py_.

- If you need MNE-Python's 3D rendering capabilities (e.g., plotting estimated
  source activity on a cortical surface) it is a good idea to install
  MNE-Python into its own virtual environment. To do this with
  `conda <anaconda_>`_:

  .. code-block:: console

      $ conda create --name=mne --channel=conda-forge mne
      $ #                   ↑↑↑                       ↑↑↑
      $ #             environment name            package name

  This will create a new ``conda`` environment called ``mne`` and install all
  dependencies into it. If you need to convert structural MRI scans into models
  of the scalp, inner/outer skull, and cortical surfaces you also need
  :doc:`FreeSurfer <freesurfer>`.

For users unfamiliar with Python, the :ref:`standard_instructions` page has
detailed instructions for different
operating systems, and there are instructions for :ref:`install-python`
if you don't already have it. The :ref:`advanced_setup` page has additional
tips and tricks for special situations (servers, notebooks, CUDA, installing
the development version, etc). The :ref:`contributing` has additional
installation instructions for (future) contributors to MNE-Python (e.g, extra
dependencies for running our tests and building our docs).

.. toctree::
    :hidden:

    pre_install
    install_python
    mne_python
    updating
    freesurfer
    advanced
.. include:: ../links.inc

.. _advanced_setup:

Advanced setup
==============

Using with IPython / Jupyter notebooks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using MNE-Python within IPython or a Jupyter notebook, we strongly
recommend using the Qt matplotlib backend for fast and correct rendering. On
Linux, for example, Qt is the only matplotlib backend for which 3D rendering
will work correctly. On macOS, certain matplotlib functions might not work as
expected on backends other than Qt. Enabling Qt can be accomplished when
starting IPython from a terminal:

.. code-block:: console

    $ ipython --matplotlib=qt

or in a Jupyter Notebook, you can use the "magic" command:

.. code-block:: ipython

    In [1]: %matplotlib qt

This will create separate pop-up windows for each figure, and has the advantage
that the 3D plots will retain rich interactivity (so, for example, you can
click-and-drag to rotate cortical surface activation maps).

If you are creating a static notebook or simply prefer Jupyter's inline plot
display, MNE-Python will work with the standard "inline" magic:

.. code-block:: ipython

    In [1]: %matplotlib inline

but some functionality will be lost. For example, PyVista scenes will still
pop-up a separate window, but only one window at a time is possible, and
interactivity within the scene is limited in non-blocking plot calls.

.. admonition:: |windows| Windows
  :class: note

  If you are using MNE-Python on Windows through IPython or Jupyter, you might
  also have to use the IPython magic command ``%gui qt`` (see `here
  <https://github.com/ipython/ipython/issues/10384>`_). For example:

  .. code-block:: ipython

     In [2]: %gui qt

If you installed the ``nb_conda_kernels`` package into your ``base``
environment (as recommended), you should be able to launch ``mne``-capable
notebooks from within the Anaconda Navigator GUI without having to explicitly
switch to the ``mne`` environment first; look for ``Python [conda env:mne]``
when choosing which notebook kernel to use. Otherwise, be sure to activate the
``mne`` environment before launching the notebook.

If you use another Python setup and you encounter some difficulties please
report them on the `MNE Forum`_ or on the `GitHub issues page`_ to get
assistance.

It is also possible to interact with the 3D plots without installing Qt by using
the notebook 3d backend:

.. code-block:: ipython

   In [1]: import mne
   In [2]: mne.viz.set_3d_backend("notebook")


The notebook 3d backend requires PyVista to be installed along with other packages,
please follow :doc:`mne_python`.


Using the development version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See :ref:`installing_main` for how to do a one-time update to the latest
development version of MNE-Python. If you plan to contribute to MNE-Python, or
just prefer to use git rather than pip to make frequent updates, there are
instructions for installing from a ``git clone`` in the :ref:`contributing`.


.. _other-py-distros:

Other Python distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^

While the `Anaconda`_ Python distribution provides many conveniences, other
distributions of Python should also work with MNE-Python.  In particular,
`Miniconda`_ is a lightweight alternative to Anaconda that is fully compatible;
like Anaconda, Miniconda includes the ``conda`` command line tool for
installing new packages and managing environments; unlike Anaconda, Miniconda
starts off with a minimal set of around 30 packages instead of Anaconda's
hundreds. See the `installation instructions for Miniconda`_ for more info.
A similar alternative is `MiniForge`_, which uses the ``conda-forge`` channel
as the default source for package installation (saving you the trouble of
typing ``--channel=conda-forge`` with each ``conda install`` command).

.. warning::

    If you have the ``PYTHONPATH`` or ``PYTHONHOME`` environment variables set,
    you may run into difficulty using Anaconda. See the
    `Anaconda troubleshooting guide`_ for more information. Note that it is
    easy to switch between ``conda``-managed Python installations and the
    system Python installation using the ``conda activate`` and ``conda
    deactivate`` commands, so you may find that after adopting Anaconda it is
    possible (indeed, preferable) to leave ``PYTHONPATH`` and ``PYTHONHOME``
    permanently unset.


It is also possible to use a system-level installation of Python (version
|min_python_version| or higher) and use ``pip`` to install MNE-Python and its
dependencies, using the provided `requirements file`_:

.. code-block:: console

    curl --remote-name https://raw.githubusercontent.com/mne-tools/mne-python/main/requirements.txt
    pip install --user -r requirements.txt

Other configurations will probably also work, but we may be unable to offer
support if you encounter difficulties related to your particular Python
installation choices.

.. _CUDA:

GPU acceleration with CUDA
^^^^^^^^^^^^^^^^^^^^^^^^^^

MNE-Python can utilize `NVIDIA CUDA GPU processing`_ to speed up some
operations (e.g. FIR filtering) by roughly an order of magnitude. To use CUDA,
first  ensure that you are running the `NVIDIA proprietary drivers`_ on your
operating system, and then do:

.. code-block:: console

    $ conda install cupy
    $ MNE_USE_CUDA=true python -c "import mne; mne.cuda.init_cuda(verbose=True)"
    Enabling CUDA with 1.55 GB available memory

If you receive a message reporting the GPU's available memory, CuPy_
is working properly. To permanently enable CUDA in MNE, you can do::

    >>> mne.utils.set_config('MNE_USE_CUDA', 'true')  # doctest: +SKIP

You can then test MNE CUDA support by running the associated test:

.. code-block:: console

    $ pytest mne/tests/test_filter.py -k cuda

If the tests pass, then CUDA should work in MNE. You can use CUDA in methods
that state that they allow passing ``n_jobs='cuda'``, such as
:meth:`mne.io.Raw.filter` and :meth:`mne.io.Raw.resample`,
and they should run faster than the CPU-based multithreading such as
``n_jobs=8``.

Off-screen rendering with MESA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On remote Linux systems, it might be possible to use MESA software rendering
(such as ``llvmpipe`` or ``swr``) for 3D visualization (with some tweaks).
For example, on CentOS 7.5 you might be able to use an environment variable
to force MESA to use modern OpenGL by using this before executing
``spyder`` or ``python``:

.. code-block:: console

    $ export MESA_GL_VERSION_OVERRIDE=3.3

Also, it's possible that different software rending backends might perform
better than others, such as using the ``llvmpipe`` backend rather than ``swr``.

MESA also can have trouble with full-screen antialiasing, which you can
disable with:

.. code-block:: console

    $ export MNE_3D_OPTION_ANTIALIAS=false

or by doing
:func:`mne.viz.set_3d_options(antialias=False) <mne.viz.set_3d_options>` within
a given Python session.

Another issue that may come up is that the MESA software itself may be out of date
in certain operating systems, for example CentOS. This may lead to incomplete
rendering of some 3D plots. A solution is described in this `Github comment <https://github.com/mne-tools/mne-python/issues/7977#issuecomment-729921035>`_.
It boils down to building a newer version (e.g., 18.3.6)
locally following a variant of `these instructions <https://xorg-team.pages.debian.net/xorg/howto/build-mesa.html#_preparing_mesa_sources>`_.
If you have CentOS 7 or newer, you can also try some `prebuilt binaries <https://osf.io/sp9qg/download>`_ we made.
After downloading the files, untar them and add them to the appropriate library paths
using the following commands:

.. code-block:: console

    $ tar xzvf mesa_18.3.6_centos_lib.tgz
    $ export LIBGL_DRIVERS_PATH="${PWD}/lib"
    $ export LD_LIBRARY_PATH="${PWD}/lib"

To check that everything went well, type the following:

.. code-block:: console

    $ glxinfo | grep "OpenGL core profile version"

which should give::

    OpenGL core profile version string: 3.3 (Core Profile) Mesa 18.3.6

Another way to check is to type:

.. code-block:: console

    $ mne sys_info

and it should show the right version of MESA::

    ...
    pyvista:       0.27.4 {pyvistaqt=0.2.0, OpenGL 3.3 (Core Profile) Mesa 18.3.6 via llvmpipe (LLVM 3.4, 256 bits)}
    ...

.. _troubleshoot_3d:

Troubleshooting 3D plots
^^^^^^^^^^^^^^^^^^^^^^^^

3D plotting trouble after upgrade on macOS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When upgrading MNE-Python from version 0.19 or lower, some macOS users may end up with
conflicting versions of some of the 3D plotting dependencies. If you plot using
the pyvista 3D backend and find that you can click-drag to rotate the brain,
but cannot adjust any of the settings sliders, it is likely that your versions
of VTK and/or QT are incompatible. This series of commands should fix it:

.. code-block:: console

    $ conda uninstall vtk
    $ pip uninstall -y pyvista
    $ conda install vtk
    $ pip install --no-cache pyvista

If you installed VTK using ``pip`` rather than ``conda``, substitute the first
line for ``pip uninstall -y vtk``.
.. include:: ../links.inc

Overview of the MNE tools suite
===============================

MNE-Python is an open-source Python module for processing, analysis, and
visualization of functional neuroimaging data (EEG, MEG, sEEG, ECoG, and
fNIRS). There are several related or interoperable software packages that you
may also want to install, depending on your analysis needs.

Related software
^^^^^^^^^^^^^^^^

- MNE-C was the initial stage of this project,
  providing a set of interrelated command-line and GUI programs focused on
  computing cortically constrained Minimum Norm Estimates from MEG and EEG
  data. These tools were written in C by Matti Hämäläinen, and are
  documented `here <MNE-C manual_>`_. See :ref:`install_mne_c` for installation
  instructions.

- MNE-Python reimplements the functionality of MNE-C, extends considerably the
  analysis and visualization capabilities, and adds support for additional data
  types like functional near-infrared spectroscopy (fNIRS). MNE-Python is
  collaboratively developed and has more than 200 contributors.

- :ref:`MNE MATLAB <mne_matlab>` provides a MATLAB interface to the .fif file
  format and other MNE data structures, and provides example MATLAB
  implementations of some of the core analysis functionality of MNE-C. It is
  distributed alongside MNE-C, and can also be downloaded from the `MNE-MATLAB
  git repository`_.

- :ref:`MNE-CPP <mne_cpp>` provides core MNE functionality implemented in
  C++ and is primarily intended for embedded and real-time applications.

There is also a growing ecosystem of other Python packages that work alongside
MNE-Python, including packages for:

.. sidebar:: Something missing?

    If you know of a package that is related but not listed here, feel free to
    :ref:`make a pull request <contributing>` to add it to this list.

- a graphical user interface for MNE-Python (`MNELAB`_)
- easily importing MEG data from the Human Connectome Project for
  use with MNE-Python (`MNE-HCP`_)
- managing MNE projects so that they comply with the `Brain
  Imaging Data Structure`_ specification (`MNE-BIDS`_)
- automatic bad channel detection and interpolation (`autoreject`_)
- convolutional sparse dictionary learning and waveform shape estimation
  (`alphaCSC`_)
- independent component analysis (ICA) with good performance on real data
  (`PICARD`_)
- phase-amplitude coupling (`pactools`_)
- representational similarity analysis (`rsa`_)
- microstate analysis (`microstate`_)
- connectivity analysis using dynamic imaging of coherent sources (DICS)
  (`conpy`_)
- general-purpose statistical analysis of M/EEG data (`eelbrain`_)
- post-hoc modification of linear models (`posthoc`_)
- a python implementation of the Preprocessing Pipeline (PREP) for EEG data
  (`pyprep`_)
- automatic multi-dipole localization and uncertainty quantification with
  the Bayesian algorithm SESAME (`sesameeg`_)
- GLM and group level analysis of near-infrared spectroscopy data (`mne-nirs`_)
- M/EEG inverse solutions using artificial neural networks (`ESINet`_)
- All-Resolutions Inference (ARI) for statistically valid circular inference
  and effect localization (`MNE-ARI`_)


What should I install?
^^^^^^^^^^^^^^^^^^^^^^

If you intend only to perform ERP, ERF, or other sensor-level analyses,
:doc:`MNE-Python <mne_python>` is all you need. If you prefer to work with
shell scripts and the Unix command line, or prefer MATLAB over Python, probably
all you need is :doc:`MNE-C <mne_c>` — the MNE MATLAB toolbox is distributed
with it — although note that the C tools and the MATLAB toolbox are less
actively developed than the MNE-Python module, and hence are considerably less
feature-complete.

If you want to transform sensor recordings into estimates of localized brain
activity, you will need MNE-Python, plus :doc:`FreeSurfer <freesurfer>` to
convert structural MRI scans into models of the scalp, inner/outer skull, and
cortical surfaces (specifically, for command-line functions
:ref:`mne flash_bem`, :ref:`mne watershed_bem`, and
:ref:`mne make_scalp_surfaces`).


Getting help
^^^^^^^^^^^^

Help with installation is available through the `MNE Forum`_. See the
:ref:`help` page for more information.


.. LINKS:

.. _MNELAB: https://github.com/cbrnr/mnelab
.. _autoreject: https://autoreject.github.io/
.. _alphaCSC: https://alphacsc.github.io/
.. _picard: https://pierreablin.github.io/picard/
.. _pactools: https://pactools.github.io/
.. _rsa: https://github.com/wmvanvliet/rsa
.. _microstate: https://github.com/wmvanvliet/mne_microstates
.. _conpy: https://aaltoimaginglanguage.github.io/conpy/
.. _eelbrain: https://eelbrain.readthedocs.io/en/stable/index.html
.. _posthoc: https://users.aalto.fi/~vanvlm1/posthoc/python/
.. _pyprep: https://github.com/sappelhoff/pyprep
.. _sesameeg: https://pybees.github.io/sesameeg
.. _mne-nirs: https://github.com/mne-tools/mne-nirs
.. _ESINet: https://github.com/LukeTheHecker/ESINet
.. _MNE-ARI: https://github.com/john-veillette/mne_ari
.. include:: ../links.inc

.. _install-python:

Installing Python
=================

MNE-Python requires Python and several Python packages. MNE-Python
version |version| requires Python version |min_python_version| or higher. We
recommend the `Anaconda`_ distribution of Python, which comes with more than
250 scientific packages pre-bundled and includes the ``conda`` command line
tool for installing new packages and managing different package sets
("environments") for different projects.

To get started, follow the `installation instructions for Anaconda`_.
When you are done, if you type the following commands in a command shell,
you should see outputs similar to the following (assuming you installed
conda to ``/home/user/anaconda3``):

.. collapse:: |linux| Linux

    .. code-block:: console

        $ conda --version && python --version
        conda 4.9.2
        Python 3.7.7 :: Anaconda, Inc.
        $ which python
        /home/user/anaconda3/bin/python
        $ which pip
        /home/user/anaconda3/bin/pip


.. collapse:: |apple| macOS

    .. code-block:: console

        $ conda --version && python --version
        conda 4.9.2
        Python 3.7.7
        $ which python
        /Users/user/opt/anaconda3/bin/python
        $ which pip
        /Users/user/opt/anaconda3/bin/pip


.. collapse:: |windows| Windows

    Most of our instructions start with ``$``, which indicates
    that the commands are designed to be run from a ``bash`` command shell.

    Windows command prompts do not expose the same command-line tools as
    ``bash`` shells, so commands like ``which`` will not work. You can test
    your installation in Windows ``cmd.exe`` shells with ``where`` instead:

    .. code-block:: doscon

        > where python
        C:\Users\user\anaconda3\python.exe
        > where pip
        C:\Users\user\anaconda3\Scripts\pip.exe

.. raw:: html

    <div width="100%" height="0 px" style="margin: 0 0 15px;"></div>

.. javascript below adapted from nilearn

.. raw:: html

     <script type="text/javascript">
     var OSName="linux-linux";
     if (navigator.userAgent.indexOf("Win")!=-1) OSName="windows-windows";
     if (navigator.userAgent.indexOf("Mac")!=-1) OSName="apple-macos";
     $(document).ready(function(){
         var element = document.getElementById("collapse_" + OSName);
         element.className += " show";
         element.setAttribute("aria-expanded", "true");
     });
     </script>


.. collapse:: |hand-paper| If you get an error or these look incorrect...
    :class: danger

    .. rubric:: If you see something like:

    ::

        conda: command not found

    It means that your ``PATH`` variable (what the system uses to find
    programs) is not set properly. In a correct installation, doing::

        $ echo $PATH
        ...:/home/user/anaconda3/bin:...

    Will show the Anaconda binary path (above) somewhere in the output
    (probably at or near the beginning), but the ``command not found`` error
    suggests that it is missing.

    On Linux or macOS, the installer should have put something
    like the following in your ``~/.bashrc`` or ``~/.bash_profile`` (or your
    ``.zprofile`` if you're using macOS Catalina or later, where the default
    shell is ``zsh``):

    .. code-block:: bash

        # >>> conda initialize >>>
        # !! Contents within this block are managed by 'conda init' !!
        __conda_setup= ...
        ...
        # <<< conda initialize <<<

    If this is missing, it is possible that you are not on the same shell that
    was used during the installation. You can verify which shell you are on by
    using the command::

        $ echo $SHELL

    If you do not find this line in the configuration file for the shell you
    are using (bash, zsh, tcsh, etc.), try running::

        conda init

    in your command shell. If your shell is not ``cmd.exe`` (Windows) or
    ``bash`` (Linux, macOS) you will need to pass the name of the shell to the
    ``conda init`` command. See ``conda init --help`` for more info and
    supported shells.

    You can also consult the Anaconda documentation and search for
    Anaconda install tips (`Stack Overflow`_ results are often helpful)
    to fix these or other problems when ``conda`` does not work.
.. include:: ../links.inc

.. _standard_instructions:

Installing MNE-Python
=====================

.. highlight:: console

Once you have Python/Anaconda installed, you have a few choices for how to
install MNE-Python.

2D plotting and sensor-level analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you only need 2D plotting capabilities with MNE-Python (i.e., most EEG/ERP
or other sensor-level analyses), you can install all you need by running
``pip install mne`` in a terminal window (on Windows, use the "Anaconda Prompt"
from the Start menu, or the "CMD.exe prompt" from within the Anaconda Navigator
GUI). This will install MNE-Python into the "base" conda environment, which
should be active by default and should already have the necessary dependencies
(``numpy``, ``scipy``, and ``matplotlib``). If you want to make use of
MNE-Python's dataset downloading functions, run ``pip install mne[data]``
instead.

A second option is to install MNE-Python into its own virtual environment
(instead of installing into conda's "base" environment). This can be done via::

    $ conda create --name=new_environment_name python=3
    $ conda activate new_environment_name
    $ pip install mne

This approach is a good choice if you want to keep a separate virtual
environment for each project. This helps with reproducibility, since each
project-specific environment will have a record of which versions of the
various software packages are installed in it (accessible with ``conda list``).

3D plotting and source analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you need MNE-Python's 3D rendering capabilities (e.g., plotting estimated
source activity on a cortical surface) it is best to install
MNE-Python into its own virtual environment, so that the extra dependencies
needed for 3D plotting stay in sync (i.e., they only get updated to versions
that are compatible with MNE-Python). See the detailed instructions below for
your operating system.

.. collapse:: |linux| Linux

   Install MNE-Python from conda-forge::

       $ conda create --name=mne --channel=conda-forge mne

.. collapse:: |apple| macOS

    Install MNE-Python into a new environment (here called ``mne``, but you can
    name the environment whatever you want)::

        $ conda create --name=mne --channel=conda-forge mne

    If you like using Jupyter notebooks, you should also update the "base"
    conda environment to include the ``nb_conda_kernels`` package; this will
    make it easier to use MNE-Python in Jupyter Notebooks launched from the
    Anaconda GUI::

        $ conda install --name=base nb_conda_kernels


.. collapse:: |windows| Windows

    Open an Anaconda command prompt, and run:

    .. code-block:: doscon

        > conda create --name=mne --channel=conda-forge mne

    If you like using Jupyter notebooks, you should also update the "base"
    conda environment to include the ``nb_conda_kernels`` package; this will
    make it easier to use MNE-Python in Jupyter Notebooks launched from the
    Anaconda GUI:

    .. code-block:: doscon

        > conda install --name base nb_conda_kernels

.. raw:: html

   <div width="100%" height="0 px" style="margin: 0 0 15px;"></div>

.. javascript below adapted from nilearn

.. raw:: html

    <script type="text/javascript">
    var OSName="linux-linux";
    if (navigator.userAgent.indexOf("Win")!=-1) OSName="windows-windows";
    if (navigator.userAgent.indexOf("Mac")!=-1) OSName="apple-macos";
    $(document).ready(function(){
        var element = document.getElementById("collapse_" + OSName);
        element.className += " show";
        element.setAttribute("aria-expanded", "true");
    });
    </script>

Installing to a headless server
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. collapse:: |server| If you are installing on a headless server...
    :class: danger

    With `pyvista`_:
    Download the `server environment file`_ and use it to create the conda
    environment::

        $ curl --remote-name https://raw.githubusercontent.com/mne-tools/mne-python/main/server_environment.yml
        $ conda create --name=mne --file=server_environment.yml

Testing your installation
^^^^^^^^^^^^^^^^^^^^^^^^^

To make sure MNE-Python installed itself and its dependencies correctly,
type the following command in a terminal::

    $ python -c "import mne; mne.sys_info()"

This should display some system information along with the versions of
MNE-Python and its dependencies. Typical output looks like this::

    Platform:      Linux-5.0.0-1031-gcp-x86_64-with-glibc2.2.5
    Python:        3.8.1 (default, Dec 20 2019, 10:06:11)  [GCC 7.4.0]
    Executable:    /home/travis/virtualenv/python3.8.1/bin/python
    CPU:           x86_64: 2 cores
    Memory:        7.8 GB

    mne:           0.21.dev0
    numpy:         1.19.0.dev0+8dfaa4a {blas=openblas, lapack=openblas}
    scipy:         1.5.0.dev0+f614064
    matplotlib:    3.2.1 {backend=Qt5Agg}

    sklearn:       0.22.2.post1
    numba:         0.49.0
    nibabel:       3.1.0
    cupy:          Not found
    pandas:        1.0.3
    dipy:          1.1.1
    pyvista:       0.25.2 {pyvistaqt=0.1.0}
    vtk:           9.0.0
    PyQt5:         5.14.1


.. collapse:: |hand-paper| If you get an error...
    :class: danger

    .. rubric:: If you see an error like:

    ::

        Traceback (most recent call last):
          File "<string>", line 1, in <module>
        ModuleNotFoundError: No module named 'mne'

    This suggests that your environment containing MNE-Python is not active.
    If you followed the setup for 3D plotting/source analysis (i.e., you
    installed to a new ``mne`` environment instead of the ``base`` environment)
    try running ``conda activate mne`` first, and try again. If this works,
    you might want to set your terminal to automatically activate the
    ``mne`` environment each time you open a terminal::

        $ echo conda activate mne >> ~/.bashrc    # for bash shells
        $ echo conda activate mne >> ~/.zprofile  # for zsh shells

If something else went wrong during installation and you can't figure it out,
check out the :doc:`advanced` page to see if your problem is discussed there.
If not, the `MNE Forum`_ is a good resources for troubleshooting installation
problems.


Python IDEs
^^^^^^^^^^^

Most users find it convenient to write and run their code in an `Integrated
Development Environment`_ (IDE). Some popular choices for scientific
Python development are:

- `Spyder`_ is a free and open-source IDE developed by and for scientists who
  use Python. It is included by default in the ``base`` environment when you
  install Anaconda, and can be started from a terminal with the command
  ``spyder`` (or on Windows or macOS, launched from the Anaconda Navigator GUI).
  It can also be installed with `dedicated installers <https://www.spyder-ide.org/#section-download>`_.
  To avoid dependency conflicts with Spyder, you should install ``mne`` in a
  separate environment, like explained in the earlier sections. Then, set
  Spyder to use the ``mne`` environment as its default interpreter by opening
  Spyder and navigating to
  :samp:`Tools > Preferences > Python Interpreter > Use the following interpreter`.
  There, paste the output of the following terminal commands::

      $ conda activate mne
      $ python -c "import sys; print(sys.executable)"

  It should be something like ``C:\Users\user\anaconda3\envs\mne\python.exe``
  (Windows) or ``/Users/user/opt/anaconda3/envs/mne/bin/python`` (macOS).

  If the Spyder console can not start because ``spyder-kernels`` is missing,
  install the required version in the ``mne`` environment with the following
  commands in the terminal::

      $ conda activate mne
      $ conda install spyder-kernels=HERE_EXACT_VERSION -c conda-forge

  Refer to the `spyder documentation <https://docs.spyder-ide.org/current/troubleshooting/common-illnesses.html#spyder-kernels-not-installed-incompatible>`_
  for more information about ``spyder-kernels`` and the version matching.

  If the Spyder graphic backend is not set to ``inline`` but to e.g. ``Qt5``,
  ``pyqt`` must be installed in the ``mne`` environment.

- `Visual Studio Code`_ (often shortened to "VS Code" or "vscode") is a
  development-focused text editor that supports many programming languages in
  addition to Python, includes an integrated terminal console, and has a rich
  ecosystem of packages to extend its capabilities. Installing
  `Microsoft's Python Extension
  <https://marketplace.visualstudio.com/items?itemName=ms-python.python>`__ is
  enough to get most Python users up and running. VS Code is free and
  open-source.
- `Atom`_ is a text editor similar to vscode, with a package ecosystem that
  includes a `Python IDE package <https://atom.io/packages/ide-python>`__ as
  well as `several <https://atom.io/packages/atom-terminal>`__
  `packages <https://atom.io/packages/atom-terminal-panel>`__
  `for <https://atom.io/packages/terminal-plus>`__
  `integrated <https://atom.io/packages/platformio-ide-terminal>`__
  `terminals <https://atom.io/packages/term3>`__. Atom is free and open-source.
- `SublimeText`_ is a general-purpose text editor that is fast and lightweight,
  and also has a rich package ecosystem. There is a package called `Terminus`_
  that provides an integrated terminal console, and a (confusingly named)
  package called "anaconda"
  (`found here <https://packagecontrol.io/packages/Anaconda>`__) that provides
  many Python-specific features. SublimeText is free (closed-source shareware).
- `PyCharm`_ is an IDE specifically for Python development that provides an
  all-in-one installation (no extension packages needed). PyCharm comes in a
  free "community" edition and a paid "professional" edition, and is
  closed-source.

.. highlight:: python

.. LINKS

.. _environment file: https://raw.githubusercontent.com/mne-tools/mne-python/main/environment.yml
.. _server environment file: https://raw.githubusercontent.com/mne-tools/mne-python/main/server_environment.yml
.. _`pyvista`: https://docs.pyvista.org/
.. _`X server`: https://en.wikipedia.org/wiki/X_Window_System
.. _`xvfb`: https://en.wikipedia.org/wiki/Xvfb
.. _`integrated development environment`: https://en.wikipedia.org/wiki/Integrated_development_environment
.. _`spyder`: https://www.spyder-ide.org/
.. _`visual studio code`: https://code.visualstudio.com/
.. _`sublimetext`: https://www.sublimetext.com/
.. _`terminus`: https://packagecontrol.io/packages/Terminus
.. _`pycharm`: https://www.jetbrains.com/pycharm/
.. _`atom`: https://atom.io/
.. _governance:

==================
Project Governance
==================

The purpose of this document is to formalize the governance process
used by the MNE-Python project in both ordinary and extraordinary
situations, and to clarify how decisions are made and how the various
elements of our community interact, including the relationship between
open source collaborative development and work that may be funded by
for-profit or non-profit entities.


The Project
===========

The MNE-Python Project (The Project) is an open source software project. The
goal of The Project is to develop open source software for analysis of
neuroscience data in Python. The Project is released under the BSD (or similar)
open source license, developed openly and is hosted publicly under the
``mne-tools`` GitHub organization.

The Project is developed by a team of distributed developers, called
Contributors. Contributors are individuals who have contributed code,
documentation, designs, or other work to the Project. Anyone can be a
Contributor. Contributors can be affiliated with any legal entity or
none. Contributors participate in the project by submitting, reviewing,
and discussing GitHub Pull Requests and Issues and participating in open
and public Project discussions on GitHub, Discourse, and other
channels. The foundation of Project participation is openness and
transparency.

The Project Community consists of all Contributors and Users of the
Project. Contributors work on behalf of and are responsible to the
larger Project Community and we strive to keep the barrier between
Contributors and Users as low as possible.

The Project is not a legal entity, nor does it currently have any formal
relationships with legal entities.


Governance model
================

This section describes the governance and leadership model of The
Project.

The foundations of Project governance are:

-  openness and transparency
-  active contribution
-  institutional neutrality


Traditionally, Project leadership was provided by a subset of Contributors,
informally called Core Developers, whose active and consistent contributions
were rewarded by granting them “commit rights” to the Project GitHub
repositories. In general, all Project decisions are made through consensus among
the Core Developers with input from the Community.

While this approach has served us well, as the Project grows we see a need for
a more formal governance model. The MNE-Python Core Developers expressed a
preference for a leadership model which includes a BDFL (Benevolent Dictator
for Life). Therefore, moving forward The Project leadership will consist of a
BDFL and Steering Council.

BDFL
----

The Project will have a BDFL (Benevolent Dictator for Life), who is currently
Alexandre Gramfort. As Dictator, the BDFL has the authority to make all final
decisions for The Project. As Benevolent, the BDFL, in practice, chooses to
defer that authority to the consensus of the community discussion channels and
the Steering Council (see below). It is expected, and in the past has been the
case, that the BDFL will only rarely assert their final authority. Because
rarely used, we refer to BDFL’s final authority as a “special” or “overriding”
vote. When it does occur, the BDFL override typically happens in situations
where there is a deadlock in the Steering Council or if the Steering Council
asks the BDFL to make a decision on a specific matter. To ensure the
benevolence of the BDFL, The Project encourages others to fork the project if
they disagree with the overall direction the BDFL is taking. The BDFL may
delegate their authority on a particular decision or set of decisions to
any other Council member at their discretion.

The BDFL can appoint their successor, but it is expected that the Steering
Council would be consulted on this decision. If the BDFL is unable to appoint a
successor, the Steering Council will make this decision — preferably by
consensus, but if needed, by a majority vote.

Note that the BDFL can step down at any time, and acting in good faith, will
also listen to serious calls to do so. Also note that the BDFL is more a role
for fallback decision making rather than that of a director/CEO.

Steering Council
----------------

The Project will have a Steering Council that consists of Project Contributors
who have produced contributions that are substantial in quality and quantity,
and sustained over at least one year. The overall role of the Council is to
ensure, through working with the BDFL and taking input from the Community, the
long-term well-being of the project, both technically and as a community.

During the everyday project activities, Council Members participate in
discussions, code review, and other project activities as peers with all other
Contributors and the Community. In these everyday activities, Council Members
do not have any special power or privilege through their membership on the
Council. However, it is expected that because of the quality and quantity of
their contributions and their expert knowledge of the Project Software and
Services, Council Members will provide useful guidance, both technical and
in terms of project direction, to potentially less experienced contributors.

The Steering Council and its Members play a special role in certain situations.
In particular, the Council may:

- Make decisions about the overall scope, vision, and direction of the project.
- Make decisions about strategic collaborations with other organizations or
  individuals.
- Make decisions about specific technical issues, features, bugs, and pull
  requests. They are the primary mechanism of guiding the code review process
  and merging pull requests.
- Make decisions about the Services that are run by The Project and manage
  those Services for the benefit of the Project and Community.
- Make decisions when regular community discussion does not produce consensus
  on an issue in a reasonable time frame.
- Update policy documents, such as this one.

Council membership
~~~~~~~~~~~~~~~~~~

To become eligible for being a Steering Council Member, an individual must be a
Project Contributor who has produced contributions that are substantial in
quality and quantity, and sustained over at least one year. Potential Council
Members are nominated by existing Council members and voted upon by the
existing Council after asking if the potential Member is interested and willing
to serve in that capacity. The Council will be initially formed from the set of
existing Core Developers who, as of May 2021, have been significantly
active over the last two years.

When considering potential Members, the Council will look at candidates with a
comprehensive view of their contributions. This will include, but is not limited
to, code, code review, infrastructure work, Discourse participation,
community help/building, education and outreach, design work, etc. We are
deliberately not setting arbitrary quantitative metrics (like “100 commits in
this repo”) to avoid encouraging behavior that plays to the metrics rather than
The Project’s overall well-being. We want to encourage a diverse array of
backgrounds, viewpoints, and talents in our team, which is why we explicitly do
not define code as the sole metric on which council membership will be
evaluated.

If a Council Member becomes inactive in the project for a period of one year,
they will be considered for removal from the Council. Before removal, inactive
Member will be approached to see if they plan on returning to active
participation. If not, they will be removed immediately upon a Council
vote. If they plan on returning to active participation soon, they will be
given a grace period of one year. If they don’t return to active participation
within that time period they will be removed by vote of the Council without
further grace period. All former Council Members can be considered for
membership again at any time in the future, like any other Project Contributor.
Retired Council Members will be listed on the project website, acknowledging
the period during which they were active in the Council.

The Council reserves the right to eject current Members, other than the BDFL,
if they are deemed to be actively harmful to the project’s well-being, and
attempts at communication and conflict resolution have failed.

A list of current Steering Council Members is maintained at the
page :ref:`governance-people`.

Conflict of interest
~~~~~~~~~~~~~~~~~~~~

It is expected that the BDFL and Council Members will be employed at a wide
range of companies, universities, and non-profit organizations. Because of this,
it is possible that Members will have a conflict of interest. Such conflicts of
interest include, but are not limited to:

- Financial interest, such as investments, employment or contracting work,
  outside of The Project that may influence their work on The Project.
- Access to proprietary information of their employer that could potentially
  leak into their work with the Project.

All members of the Council, BDFL included, shall disclose to the rest of the
Council any conflict of interest they may have. Members with a conflict of
interest in a particular issue may participate in Council discussions on that
issue, but must recuse themselves from voting on the issue. If the BDFL has
recused themself for a particular decision, the Council will appoint a
substitute BDFL for that decision.

Private communications of the Council
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unless specifically required, all Council discussions and activities will be
public and done in collaboration and discussion with the Project Contributors
and Community. The Council will have a private communication channel that will be used
sparingly and only when a specific matter requires privacy. When private
communications and decisions are needed, the Council will do its best to
summarize those to the Community after removing personal/private/sensitive
information that should not be posted to the public internet.

Council decision making
~~~~~~~~~~~~~~~~~~~~~~~

If it becomes necessary for the Steering Council to produce a formal
decision, then they will use a form of the `Apache Foundation voting
process <https://www.apache.org/foundation/voting.html>`_. This is a
formalized version of consensus, in which +1 votes indicate agreement,
-1 votes are vetoes (and must be accompanied with a rationale),
and one can also vote fractionally (e.g. -0.5, +0.5) if one
wishes to express an opinion without registering a full veto. These
numeric votes are also often used informally as a way of getting a
general sense of people's feelings on some issue, and should not
normally be taken as formal votes. A formal vote only occurs if
explicitly declared, and if this does occur, then the vote should be held
open for long enough to give all interested Council Members a chance to
respond — at least one week.

In practice, we anticipate that for most Steering Council decisions
(e.g., voting in new members) a more informal process will suffice.


Institutional Partners and funding
==================================

The Steering Council is the primary leadership for the project. No
outside institution, individual, or legal entity has the ability to own,
control, usurp, or influence the project other than by participating in
the Project as Contributors and Council Members. However, because
institutions can be an important funding mechanism for the project, it
is important to formally acknowledge institutional participation in the
project. These are Institutional Partners.

An Institutional Contributor is any individual Project Contributor who
contributes to the project as part of their official duties at an
Institutional Partner. Likewise, an Institutional Council Member is any
Project Steering Council Member who contributes to the project as part
of their official duties at an Institutional Partner.

With these definitions, an Institutional Partner is any recognized legal
entity in any country that employs at least 1 Institutional Contributor or
Institutional Council Member. Institutional Partners can be for-profit or
non-profit entities.

Institutions become eligible to become an Institutional Partner by
employing individuals who actively contribute to The Project as part of
their official duties. To state this another way, the only way for a
Partner to influence the project is by actively contributing to the open
development of the project, in equal terms to any other member of the
community of Contributors and Council Members. Merely using Project
Software in institutional context does not allow an entity to become an
Institutional Partner. Financial gifts do not enable an entity to become
an Institutional Partner. Once an institution becomes eligible for
Institutional Partnership, the Steering Council must nominate and
approve the Partnership.

If, at some point, an existing Institutional Partner stops having any
contributing employees, then a one year grace period commences. If, at
the end of this one-year period, they continue not to have any
contributing employees, then their Institutional Partnership will
lapse, and resuming it will require going through the normal process
for new Partnerships.

An Institutional Partner is free to pursue funding for their work on The
Project through any legal means. This could involve a non-profit
organization raising money from private foundations and donors or a
for-profit company building proprietary products and services that
leverage Project Software and Services. Funding acquired by
Institutional Partners to work on The Project is called Institutional
Funding. However, no funding obtained by an Institutional Partner can
override the Steering Council. If a Partner has funding to do MNE-Python work
and the Council decides to not pursue that work as a project, the
Partner is free to pursue it on their own. However, in this situation,
that part of the Partner’s work will not be under the MNE-Python umbrella and
cannot use the Project trademarks in any way that suggests a formal
relationship.

Institutional Partner benefits are:

- optional acknowledgement on the MNE-Python website and in talks
- ability to acknowledge their own funding sources on the MNE-Python
  website and in talks
- ability to influence the project through the participation of their
  Council Member
- invitation of the Council Members to MNE-Python Developer Meetings

A list of current Institutional Partners is maintained at the page
:ref:`supporting-institutions`.


Document history
================

https://github.com/mne-tools/mne-python/commits/main/doc/overview/governance.rst


Acknowledgements
================

Substantial portions of this document were adapted from the
`SciPy project's governance document
<https://github.com/scipy/scipy/blob/main/doc/source/dev/governance.rst>`_,
which in turn was adapted from
`Jupyter/IPython project's governance document
<https://github.com/jupyter/governance/blob/master/governance.md>`_ and
`NumPy's governance document
<https://github.com/numpy/numpy/blob/master/doc/source/dev/governance/governance.rst>`_.

License
=======

To the extent possible under law, the authors have waived all
copyright and related or neighboring rights to the MNE-Python project
governance document, as per the `CC-0 public domain dedication / license
<https://creativecommons.org/publicdomain/zero/1.0/>`_.
.. include:: ../links.inc

.. _learn-python:

Getting started with Python
===========================

`Python`_ is a modern general-purpose object-oriented high-level programming
language. There are many general introductions to Python online; here are a
few:

- The official `Python tutorial <https://docs.python.org/3/tutorial/index.html>`__
- W3Schools `Python tutorial <https://www.w3schools.com/python/>`__
- Software Carpentry's `Python lesson <http://swcarpentry.github.io/python-novice-inflammation/>`_

Additionally, here are a couple tutorials focused on scientific programming in
Python:

- the `SciPy Lecture Notes <http://scipy-lectures.org/>`_
- `NumPy for MATLAB users <https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html>`_

There are also many video tutorials online, including `videos from the annual
SciPy conferences
<https://www.youtube.com/user/EnthoughtMedia/playlists?shelf_id=1&sort=dd&view=50>`_.
One of those is a `Python introduction for complete beginners
<https://www.youtube.com/watch?v=Xmxy2NU9LOI>`_, but there are many more
lectures on advanced topics available as well.
.. _implementation:

Algorithms and other implementation details
===========================================

This page describes some of the technical details of MNE-Python implementation.

.. _units:

Internal representation (units)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../_includes/units.rst
   :start-after: units-begin-content


.. _precision:

Floating-point precision
^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../_includes/precision.rst
   :start-after: precision-begin-content


.. _channel-types:

Supported channel types
^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../_includes/channel_types.rst
   :start-after: channel-types-begin-content


.. _data-formats:

Supported data formats
^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../_includes/data_formats.rst
   :start-after: data-formats-begin-content


.. _dig-formats:

Supported formats for digitized 3D locations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../_includes/dig_formats.rst
   :start-after: dig-formats-begin-content


.. _memory:

Memory-efficient I/O
^^^^^^^^^^^^^^^^^^^^

.. include:: ../_includes/memory.rst
   :start-after: memory-begin-content


.. _channel-interpolation:

Bad channel repair via interpolation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../_includes/channel_interpolation.rst
   :start-after: channel-interpolation-begin-content
   :end-before: channel-interpolation-end-content


.. _maxwell:

Maxwell filtering
^^^^^^^^^^^^^^^^^

MNE-Python's implementation of Maxwell filtering is described in the
:ref:`tut-artifact-sss` tutorial.


.. _ssp-method:

Signal-Space Projection (SSP)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../_includes/ssp.rst
   :start-after: ssp-begin-content


.. _bem-model:

The Boundary Element Model (BEM)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../_includes/bem_model.rst
   :start-after: bem-begin-content


.. _ch_forward:

The forward solution
^^^^^^^^^^^^^^^^^^^^

.. include:: ../_includes/forward.rst
   :start-after: forward-begin-content
   :end-before: forward-end-content

.. _minimum_norm_estimates:

The minimum-norm current estimates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../_includes/inverse.rst
   :start-after: inverse-begin-content
   :end-before: inverse-end-content


.. _ch_morph:

Morphing and averaging source estimates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: ../_includes/morph.rst
   :start-after: morph-begin-content


References
^^^^^^^^^^
.. footbibliography::
.. _cite:

How to cite MNE-Python
======================

Citing the software
-------------------

To cite specific version numbers of the software, you can use the DOIs provided
by `Zenodo <https://doi.org/10.5281/zenodo.592483>`_. Additionally, we ask that
when citing the MNE-Python package, you cite the canonical journal article
reference :footcite:`GramfortEtAl2013a`:

.. footbibliography::

.. collapse:: |quote-left| BibTeX for MNE-Python
    :class: info

    .. include:: ../references.bib
        :code: bibtex
        :start-after: % MNE-Python reference
        :end-before: % everything else


Citing the inverse imaging algorithms
-------------------------------------

To cite MNE-C or the inverse imaging implementations provided by the MNE
software, please use :footcite:`GramfortEtAl2014`:

.. footbibliography::


.. collapse:: |quote-left| BibTeX for inverse algorithms / MNE-C
    :class: info

    .. include:: ../references.bib
        :code: bibtex
        :start-after: % MNE-C reference
        :end-before: % MNE-Python reference


Citing other algorithms
-----------------------

Depending on your research topic, it may also be appropriate to cite related
method papers, some of which are listed in the documentation strings of the
relevant functions or methods. All references cited in the MNE-Python codebase
and documentation are collected in the :ref:`general_bibliography`.
.. include:: ../links.inc

.. _documentation_overview:

Documentation overview
======================

.. note::

   If you haven't already installed Python and MNE-Python, here are the
   installation instructions for :ref:`Python <install-python>` and
   :ref:`MNE-Python <standard_instructions>`, and some
   resources for :doc:`learn_python`.


The documentation for MNE-Python is divided into four main sections:

1. The :doc:`../auto_tutorials/index` provide narrative explanations, sample
   code, and expected output for the most common MNE-Python analysis tasks. The
   emphasis is on thorough explanations that get new users up to speed quickly,
   at the expense of covering only a limited number of topics.

2. The :doc:`How-to Examples <../auto_examples/index>` provides working code
   samples demonstrating various analysis and visualization techniques. These
   examples often lack the narrative explanations seen in the tutorials, but
   can be a useful way to discover new analysis or plotting ideas, or to see
   how a particular technique you've read about can be applied using
   MNE-Python.

3. The :doc:`../glossary` provides short definitions of MNE-Python-specific
   vocabulary and general neuroimaging concepts. The glossary is often a good
   place to look if you don't understand a term or acronym used somewhere else
   in the documentation.

4. The :doc:`API reference <../python_reference>` provides documentation for
   the classes, functions and methods in the MNE-Python codebase. This is the
   same information that is rendered when running
   :samp:`help(mne.{<function_name>})` in an interactive Python session, or
   when typing :samp:`mne.{<function_name>}?` in an IPython session or Jupyter
   notebook.

The rest of the MNE-Python documentation pages (parts outside of the four
categories above) are shown in the navigation menu, including the
:ref:`list of example datasets<datasets>`,
:ref:`implementation details<implementation>`, and more.
Documentation for the related C and MATLAB tools are available here:

- :ref:`MNE-MATLAB <mne_matlab>` (HTML)
- `MNE-C <MNE-C manual_>`_ (PDF)

.. toctree::
   :hidden:

   Tutorials<../auto_tutorials/index>
   Examples<../auto_examples/index>
   ../glossary
   Implementation details<implementation>
   design_philosophy
   Example datasets<datasets_index>
   Command-line tools<../generated/commands>
   migrating
   cookbook
   cite
   ../cited
.. include:: ../links.inc

Roadmap
=======

This page describes some of the major medium- to long-term goals for
MNE-Python. These are goals that require substantial effort and/or
API design considerations. Some of these may be suitable for Google Summer of
Code projects, while others require more extensive work.

.. contents:: Page contents
   :local:

Open
----

.. _time-frequency-viz:

Time-frequency visualization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We should implement a viewer for interactive visualization of volumetric
source-time-frequency (5-D) maps on MRI slices (orthogonal 2D viewer).
`NutmegTrip <https://github.com/fieldtrip/fieldtrip/tree/master/contrib/nutmegtrip>`__
(written by Sarang Dalal) provides similar functionality in Matlab in
conjunction with FieldTrip. Example of NutmegTrip's source-time-frequency mode
in action (click for link to YouTube):

.. image:: https://i.ytimg.com/vi/xKdjZZphdNc/maxresdefault.jpg
   :target: https://www.youtube.com/watch?v=xKdjZZphdNc
   :width: 50%

Clustering statistics API
^^^^^^^^^^^^^^^^^^^^^^^^^
The current clustering statistics code has limited functionality. It should be
re-worked to create a new ``cluster_based_statistic`` or similar function.
In particular, the new API should:

1. Support mixed within- and between-subjects designs, different statistical
   functions, etc. This should be done via a ``design`` argument that mirrors
   :func:`patsy.dmatrices` or similar community standard (e.g., this is what
   is used by :class:`statsmodels.regression.linear_model.OLS`).
2. Have clear tutorials showing how different contrasts can be done (toy data).
3. Have clear tutorials showing some common analyses on real data (time-freq,
   sensor space, source space, etc.)
4. Not introduce any significant speed penalty (e.g., < 10% slower) compared
   to the existing, more specialized/limited functions.

More details are in :gh:`4859`.

Access to open EEG/MEG databases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We should improve the access to open EEG/MEG databases via the
:mod:`mne.datasets` module, in other words improve our dataset fetchers.
We have physionet, but much more. Having a consistent API to access multiple
data sources would be great. See :gh:`2852` and :gh:`3585` for some ideas,
as well as:

- `OpenNEURO <https://openneuro.org>`__
    "A free and open platform for sharing MRI, MEG, EEG, iEEG, and ECoG data."
    See for example :gh:`6687`.
- `Human Connectome Project Datasets <http://www.humanconnectome.org/data>`__
    Over a 3-year span (2012-2015), the Human Connectome Project (HCP) scanned
    1,200 healthy adult subjects. The available data includes MR structural
    scans, behavioral data and (on a subset of the data) resting state and/or
    task MEG data.
- `MMN dataset <http://www.fil.ion.ucl.ac.uk/spm/data/eeg_mmn>`__
    Used for tutorial/publications applying DCM for ERP analysis using SPM.
- Kymata datasets
    Current and archived EMEG measurement data, used to test hypotheses in the
    Kymata atlas. The participants are healthy human adults listening to the
    radio and/or watching films, and the data is comprised of (averaged) EEG
    and MEG sensor data and source current reconstructions.
- `BNCI Horizon <http://bnci-horizon-2020.eu/database/data-sets>`__
    BCI datasets.

Integrate OpenMEEG via improved Python bindings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`OpenMEEG <http://openmeeg.github.io>`__ is a state-of-the art solver for
forward modeling in the field of brain imaging with MEG/EEG. It solves
numerically partial differential equations (PDE). It is written in C++ with
Python bindings written in `SWIG <https://github.com/openmeeg/openmeeg>`__.
The ambition of the project is to integrate OpenMEEG into MNE offering to MNE
the ability to solve more forward problems (cortical mapping, intracranial
recordings, etc.). Some software tasks that shall be completed:

- Cleanup Python bindings (remove useless functions, check memory managements,
  etc.)
- Write example scripts for OpenMEEG that automatically generate web pages as
  for `MNE <http://martinos.org/mne/stable/auto_examples/index.html>`__
- Understand how MNE encodes info about sensors (location, orientation,
  integration points etc.) and allow OpenMEEG to be used.
- Help package OpenMEEG for Debian/Ubuntu
- Help manage `the continuous integration system
  <https://ci.inria.fr/>`__


In progress
-----------

Diversity, Equity, and Inclusion (DEI)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MNE-Python is committed to recruiting and retaining a diverse pool of
contributors, see :gh:`8221`.

First-class OPM support
^^^^^^^^^^^^^^^^^^^^^^^
MNE-Python has support for reading some OPM data formats such as FIF, but
support is still rudimentary. Support should be added for other manufacturers,
and standard (and/or novel) preprocessing routines should be added to deal with
coregistration adjustment, forward modeling, and OPM-specific artifacts.

Deep source modeling
^^^^^^^^^^^^^^^^^^^^
Existing source modeling and inverse routines are not explicitly designed to
deal with deep sources. Advanced algorithms exist from MGH for enhancing
deep source localization, and these should be implemented and vetted in
MNE-Python.

Better sEEG/ECoG/DBS support
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Some support already exists for iEEG electrodes in MNE-Python thanks in part
to standard abstractions. However, iEEG-specific pipeline steps (e.g.,
electrode localization) and visualizations (e.g., per-shaft topo plots,
:ref:`time-frequency-viz`) are missing. MNE-Python should work with members of
the ECoG/sEEG community to work with or build in existing tools, and extend
native functionality for depth electrodes.

Time-frequency classes
^^^^^^^^^^^^^^^^^^^^^^
Our current codebase implements classes related to :term:`TFRs <tfr>` that
remain incomplete. We should implement new classes from the ground up
that can hold frequency data (``Spectrum``), cross-spectral data
(``CrossSpectrum``), multitaper estimates (``MultitaperSpectrum``), and
time-varying estimates (``Spectrogram``). These should work for
continuous, epoched, and averaged sensor data, as well as source-space brain
data.

See related issues :gh:`6290`, :gh:`7671`, :gh:`8026`, :gh:`8724`, :gh:`9045`,
and PRs :gh:`6609`, :gh:`6629`, :gh:`6672`, :gh:`6673`, :gh:`8397`, and
:gh:`8892`.

Pediatric and clinical MEG pipelines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MNE-Python is in the process of providing automated analysis of BIDS-compliant
datasets, see `MNE-BIDS-Pipeline`_. By incorporating functionality from the
`mnefun <https://labsn.github.io/mnefun/overview.html>`__ pipeline,
which has been used extensively for pediatric data analysis at `I-LABS`_,
better support for pediatric and clinical data processing can be achieved.
Multiple processing steps (e.g., eSSS), sanity checks (e.g., cHPI quality),
and reporting (e.g., SSP joint plots, SNR plots) will be implemented.

Statistics efficiency
^^^^^^^^^^^^^^^^^^^^^
A key technique in functional neuroimaging analysis is clustering brain
activity in adjacent regions prior to statistical analysis. An important
clustering algorithm — threshold-free cluster enhancement (TFCE) — currently
relies on computationally expensive permutations for hypothesis testing.
A faster, probabilistic version of TFCE (pTFCE) is available, and we are in the
process of implementing this new algorithm.

3D visualization
^^^^^^^^^^^^^^^^
Historically we have used Mayavi for 3D visualization, but have faced
limitations and challenges with it. We should work to use some other backend
(e.g., PyVista) to get major improvements, such as:

1. *Proper notebook support (through ipyvtklink)* (complete)
2. *Better interactivity with surface plots* (complete)
3. Time-frequency plotting (complementary to volume-based
   :ref:`time-frequency-viz`)
4. Integration of multiple functions as done in ``mne_analyze``, e.g.,
   simultaneous source estimate viewing, field map
   viewing, head surface display, etc. These are all currently available in
   separate functions, but we should be able to combine them in a single plot
   as well.

The meta-issue for tracking to-do lists for surface plotting is :gh:`7162`.

.. _documentation-updates:

Documentation updates
^^^^^^^^^^^^^^^^^^^^^
Our documentation has many minor issues, which can be found under the tag
:gh:`labels/DOC`.


Completed
---------

Distributed computing support
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`MNE-BIDS-Pipeline`_ has been enhanced with support for cloud computing
via `Dask`_ and :doc:`joblib <joblib:auto_examples/parallel/distributed_backend_simple>`.
After configuring Dask to use local or remote distributed computing resources,
MNE-BIDS-Pipeline can readily make use of remote workers to parallelize
processing across subjects.

2D visualization
^^^^^^^^^^^^^^^^
`This goal <https://mne.tools/0.22/overview/roadmap.html#2d-visualization>`__
was completed under CZI `EOSS2`_. Some additional enhancements that could also
be implemented are listed in :gh:`7751`.

Tutorial / example overhaul
^^^^^^^^^^^^^^^^^^^^^^^^^^^
`This goal <https://mne.tools/0.22/overview/roadmap.html#tutorial-example-overhaul>`__
was completed under CZI `EOSS2`_. Ongoing documentation needs are listed in
:ref:`documentation-updates`.

Cluster computing images
^^^^^^^^^^^^^^^^^^^^^^^^
As part of `this goal <https://mne.tools/0.22/overview/roadmap.html#cluster-computing>`__,
we created docker images suitable for cloud computing via `MNE-Docker`_.

.. _I-LABS: http://ilabs.washington.edu/
.. _cookbook:

==========================
The typical M/EEG workflow
==========================

Overview
========

This section describes a typical MEG/EEG workflow, eventually up to source
reconstruction. The workflow is summarized in :ref:`flow_diagram`.
References below refer to Python functions and objects.

.. _flow_diagram:

.. figure:: images/flow_diagram.svg
    :alt: MNE Workflow Flowchart
    :align: center

    **Workflow of the MNE software**


Preprocessing
=============
The following MEG and EEG data preprocessing steps are recommended:

- Bad channels in the MEG and EEG data must be identified, see :ref:`marking_bad_channels`.

- The data has to be filtered to the desired passband.

- Artifacts should be suppressed (e.g., using ICA or SSP).

.. _marking_bad_channels:

Marking bad channels
--------------------

Sometimes some MEG or EEG channels are not functioning properly
for various reasons. These channels should be excluded from
analysis by marking them bad as::

    >>> raw.info['bads'] = ['MEG2443']  # doctest: +SKIP

Especially if a channel does not show
a signal at all (flat) it is important to exclude it from the
analysis, since its noise estimate will be unrealistically low and
thus the current estimate calculations will give a strong weight
to the zero signal on the flat channels and will essentially vanish.
It is also important to exclude noisy channels because they can
possibly affect others when signal-space projections or EEG average electrode
reference is employed. Noisy bad channels can also adversely affect
averaging and noise-covariance matrix estimation by causing
unnecessary rejections of epochs.

Recommended ways to identify bad channels are:

- Observe the quality of data during data
  acquisition and make notes of observed malfunctioning channels to
  your measurement protocol sheet.

- View the on-line averages and check the condition of the channels.

- Compute preliminary off-line averages with artifact rejection,
  SSP/ICA, and EEG average electrode reference computation
  off and check the condition of the channels.

- View raw data with :func:`mne.io.Raw.plot` without SSP/ICA
  enabled and identify bad channels.

.. note:: It is strongly recommended that bad channels are identified and
          marked in the original raw data files. If present in the raw data
          files, the bad channel selections will be automatically transferred
          to averaged files, noise-covariance matrices, forward solution
          files, and inverse operator decompositions.

Artifact suppression
--------------------

SSP
###

The Signal-Space Projection (SSP) is one approach to rejection
of external disturbances in software. Unlike many other
noise-cancellation approaches, SSP does
not require additional reference sensors to record the disturbance
fields. Instead, SSP relies on the fact that the magnetic field
distributions generated by the sources in the brain have spatial
distributions sufficiently different from those generated by external
noise sources. Furthermore, it is implicitly assumed that the linear
space spanned by the significant external noise patterns has a low
dimension.

SSP-based rejection is often done using the
:func:`mne.preprocessing.compute_proj_ecg` and
:func:`mne.preprocessing.compute_proj_eog` methods, see
:ref:`tut-projectors-background` and :ref:`tut-artifact-ssp` for more
information.

ICA
###

Many M/EEG signals including biological artifacts reflect non-Gaussian
processes. Therefore PCA-based artifact rejection will likely perform worse at
separating the signal from noise sources.

ICA-based artifact rejection is done using the :class:`mne.preprocessing.ICA`
class, see the :ref:`ica` section for more information.


Epoching and evoked data
========================

Epoching of raw data is done using events, which define a ``t=0`` for your
data chunks. Event times stamped to the acquisition software can be extracted
using :func:`mne.find_events`::

    >>> events = mne.find_events(raw)  # doctest: +SKIP

The ``events`` array can then be modified, extended, or changed if necessary.
If the original trigger codes and trigger times are correct for the analysis
of interest, :class:`mne.Epochs` for the first event type (``1``) can be
constructed using::

    >>> reject = dict(grad=4000e-13, mag=4e-12, eog=150e-6)  # doctest: +SKIP
    >>> epochs = mne.Epochs(raw, events, event_id=1, tmin=-0.2, tmax=0.5,  # doctest: +SKIP
    >>>                     proj=True, picks=picks, baseline=(None, 0),  # doctest: +SKIP
    >>>                     preload=True, reject=reject)  # doctest: +SKIP

.. note:: The rejection thresholds (set with argument ``reject``) are defined
          in T / m for gradiometers, T for magnetometers and V for EEG and EOG
          channels.


Rejection using annotations
---------------------------

The reject keyword of :class:`mne.Epochs` is used for rejecting bad epochs
based on peak-to-peak thresholds. Bad segments of data can also be rejected
by marking segments of raw data with annotations. See
:ref:`tut-reject-data-spans` and :class:`mne.Annotations` for more .

Once the :class:`mne.Epochs` are constructed, they can be averaged to obtain
:class:`mne.Evoked` data as::

    >>> evoked = epochs.average()  # doctest: +SKIP


Source localization
===================

MNE makes extensive use of the FreeSurfer file structure for analysis.
Before starting data analysis, we recommend setting up the environment
variable ``SUBJECTS_DIR`` (or set it permanently using :func:`mne.set_config`)
to select the directory under which the anatomical MRI data are stored.
This makes it so that the ``subjects_dir`` argument does not need to
be passed to many functions.

Anatomical information
----------------------

.. _CHDBBCEJ:

Cortical surface reconstruction with FreeSurfer
###############################################

The first processing stage is the creation of various surface
reconstructions with FreeSurfer. The recommended FreeSurfer workflow
is summarized on the `FreeSurfer wiki pages <https://surfer.nmr.mgh.harvard.edu/fswiki/RecommendedReconstruction>`_. See
also this information :ref:`tut-freesurfer-reconstruction`.

.. _setting_up_source_space:

Setting up the source space
###########################

This stage consists of the following:

- Creating a suitable decimated dipole grid on the white matter surface.

- Creating the source space file in fif format.

This is accomplished with using :func:`mne.setup_source_space` and
:func:`mne.write_source_spaces`. These assume that the anatomical MRI processing
has been completed as described in :ref:`CHDBBCEJ`.

.. _BABGCDHA:

.. table:: Recommended subdivisions of an icosahedron and an octahedron for
           the creation of source spaces. The approximate source spacing and
           corresponding surface area have been calculated assuming a
           1000-cm2 surface area per hemisphere.

    ===========  ======================  ===================  =============================
    ``spacing``  Sources per hemisphere  Source spacing / mm  Surface area per source / mm2
    ===========  ======================  ===================  =============================
    ``'oct5'``   1026                    9.9                  97
    ``'ico4'``   2562                    6.2                  39
    ``'oct6'``   4098                    4.9                  24
    ``'ico5'``   10242                   3.1                  9.8
    ===========  ======================  ===================  =============================

For example, to create the reconstruction geometry for ``subject='sample'``
with a ~5-mm spacing between the grid points, say::

    >>> src = setup_source_space('sample', spacing='oct6')  # doctest: +SKIP
    >>> write_source_spaces('sample-oct6-src.fif', src)  # doctest: +SKIP

This creates the source spaces and writes them to disk.

:ref:`plot_forward_source_space` illustrates how the source space is used to
compute the forward model.

.. _CHDBJCIA:

Creating the BEM model meshes
#############################

Calculation of the forward solution using the boundary-element
model (BEM) requires that the surfaces separating regions of different
electrical conductivities are tessellated with suitable surface
elements. Our BEM software employs triangular tessellations. Therefore,
prerequisites for BEM calculations are the segmentation of the MRI
data and the triangulation of the relevant surfaces.

For MEG computations, a reasonably accurate solution can
be obtained by using a single-compartment BEM assuming the shape
of the intracranial volume. For EEG, the standard model contains
the intracranial space, the skull, and the scalp.

At present, no bulletproof method exists for creating the
triangulations. Feasible approaches are described in :ref:`bem-model`.

.. _BABDBBFC:

Setting up the head surface triangulation files
###############################################

The segmentation algorithms described in :ref:`bem-model` produce
either FreeSurfer surfaces or triangulation
data in text. Before proceeding to the creation of the boundary
element model, standard files for FreeSurfer surfaces must be present:

1. **inner_skull.surf** contains the inner skull triangulation.

2. **outer_skull.surf** contains the outer skull triangulation.

3. **outer_skin.surf** contains the head surface triangulation.

.. _CIHDBFEG:

Setting up the boundary-element model
#####################################

This stage sets up the subject-dependent data for computing
the forward solutions:"

    >>> model = make_bem_model('sample')  # doctest: +SKIP
    >>> write_bem_surfaces('sample-5120-5120-5120-bem.fif', model)  # doctest: +SKIP

Where ``surfaces`` is a list of BEM surfaces that have each been read using
:func:`mne.read_surface`. This step also checks that the input surfaces
are complete and that they are topologically correct, *i.e.*,
that the surfaces do not intersect and that the surfaces are correctly
ordered (outer skull surface inside the scalp and inner skull surface
inside the outer skull).

This step assigns the conductivity values to the BEM compartments.
For the scalp and the brain compartments, the default is 0.3 S/m.
The default skull conductivity is 50 times smaller, *i.e.*,
0.006 S/m. Recent publications report a range of skull conductivity ratios
ranging from 1:15 :footcite:`OostendorpEtAl2000` to 1:25 - 1:50
:footcite:`GoncalvesEtAl2003,LewEtAl2009`. The MNE default ratio 1:50 is based
on the typical values reported in :footcite:`GoncalvesEtAl2003`, since their
approach is based on comparison of SEF/SEP measurements in a BEM model.
The variability across publications may depend on individual variations
but, more importantly, on the precision of the skull compartment
segmentation.

.. note:: To produce single layer BEM models (--homog flag in the C command
          line tools) pass a list with one single conductivity value,
          e.g. ``conductivities=[0.3]``.

Using this model, the BEM solution can be computed using
:func:`mne.make_bem_solution` as::

    >>> bem_sol = make_bem_solution(model)  # doctest: +SKIP
    >>> write_bem_solution('sample-5120-5120-5120-bem-sol.fif', bem_sol)  # doctest: +SKIP

After the BEM is set up it is advisable to check that the
BEM model meshes are correctly positioned using *e.g.*
:func:`mne.viz.plot_alignment` or :class:`mne.Report`.

.. note:: Up to this point all processing stages depend on the
          anatomical (geometrical) information only and thus remain
          identical across different MEG studies.

.. note:: If you use custom head models you might need to set the ``ico=None``
          parameter to ``None`` and skip subsampling of the surface.


.. _CHDBEHDC:

Aligning coordinate frames
--------------------------

The calculation of the forward solution requires knowledge
of the relative location and orientation of the MEG/EEG and MRI
coordinate systems (see :ref:`head_device_coords`). The head coordinate
frame is defined by identifying the fiducial landmark locations,
making the origin and orientation of the head coordinate system
slightly user dependent. As a result, it is safest to reestablish
the definition of the coordinate transformation computation
for each experimental session, *i.e.*, each time when new head
digitization data are employed.

The corregistration is stored in ``-trans.fif`` file. If is present,
you can follow :ref:`tut-source-alignment` to validate its correctness.
If the ``-trans.fif`` is not present or the alignment is not correct
you need to use :func:`mne.gui.coregistration` (or its convenient command line
equivalent :ref:`mne coreg`) to generate it.

.. XXX: It would be good to link to the ``-trans.fif`` file description

.. warning:: This step is important. If the alignment of the
             coordinate frames is inaccurate all subsequent processing
             steps suffer from the error. Therefore, this step should be
             performed by the person in charge of the study or by a trained
             technician. Written or photographic documentation of the alignment
             points employed during the MEG/EEG acquisition can also be
             helpful.

.. _computing_the_forward_solution:

Computing the forward solution
------------------------------

After the MRI-MEG/EEG alignment has been set, the forward
solution, *i.e.*, the magnetic fields and electric
potentials at the measurement sensors and electrodes due to dipole
sources located on the cortex, can be calculated with help of
:func:`mne.make_forward_solution` as::

    >>> fwd = make_forward_solution(raw.info, fname_trans, src, bem_sol)  # doctest: +SKIP

Computing the noise-covariance matrix
-------------------------------------

The MNE software employs an estimate of the noise-covariance
matrix to weight the channels correctly in the calculations. The
noise-covariance matrix provides information about field and potential
patterns representing uninteresting noise sources of either human
or environmental origin.

The noise covariance matrix can be calculated in several
ways:

- Employ the individual epochs during
  off-line averaging to calculate the full noise covariance matrix.
  This is the recommended approach for evoked responses, *e.g.* using
  :func:`mne.compute_covariance`::

      >>> cov = mne.compute_covariance(epochs, method='auto')  # doctest: +SKIP

- Employ empty room data (collected without the subject) to
  calculate the full noise covariance matrix. This is recommended
  for analyzing ongoing spontaneous activity. This can be done using
  :func:`mne.compute_raw_covariance` as::

      >>> cov = mne.compute_raw_covariance(raw_erm)  # doctest: +SKIP

- Employ a section of continuous raw data collected in the presence
  of the subject to calculate the full noise covariance matrix. This
  is the recommended approach for analyzing epileptic activity. The
  data used for this purpose should be free of technical artifacts
  and epileptic activity of interest. The length of the data segment
  employed should be at least 20 seconds. One can also use a long
  (``*> 200 s``) segment of data with epileptic spikes present provided
  that the spikes occur infrequently and that the segment is apparently
  stationary with respect to background brain activity. This can also
  use :func:`mne.compute_raw_covariance`.

.. _CIHCFJEI:

Calculating the inverse operator
--------------------------------

The MNE software doesn't calculate the inverse operator
explicitly but rather computes an SVD of a matrix composed of the
noise-covariance matrix, the result of the forward calculation,
and the source covariance matrix. This approach has the benefit
that the regularization parameter ('SNR') can
be adjusted easily when the final source estimates or dSPMs are
computed. For mathematical details of this approach,
please consult :ref:`minimum_norm_estimates`.

This computation stage can be done by using
:func:`mne.minimum_norm.make_inverse_operator` as::

    >>> inv = mne.minimum_norm.make_inverse_operator(raw.info, fwd, cov, loose=0.2)  # doctest: +SKIP

Creating source estimates
-------------------------

Once all the preprocessing steps described above have been
completed, the inverse operator computed can be applied to the MEG
and EEG data as::

    >>> stc = mne.minimum_norm.apply_inverse(evoked, inv, lambda2=1. / 9.)  # doctest: +SKIP

And the results can be viewed as::

    >>> stc.plot()  # doctest: +SKIP

Group analyses
--------------

Group analysis is facilitated by morphing source estimates, which can be
done *e.g.*, to ``subject='fsaverage'`` as::

    >>> morph = mne.compute_source_morph(stc, subject_from='sample', subject_to='fsaverage')  # doctest: +SKIP
    >>> stc_fsaverage = morph.apply(stc)  # doctest: +SKIP

See :ref:`ch_morph` for more information.


References
==========

.. footbibliography::
.. _datasets:

Datasets Overview
#################

.. sidebar:: Contributing datasets to MNE-Python

    Do not hesitate to contact MNE-Python developers on the
    `MNE Forum <https://mne.discourse.group>`_ to discuss the possibility of
    adding more publicly available datasets.

All the dataset fetchers are available in :mod:`mne.datasets`. To download any of the datasets,
use the ``data_path`` (fetches full dataset) or the ``load_data`` (fetches dataset partially) functions.

All fetchers will check the default download location first to see if the dataset
is already on your computer, and only download it if necessary. The default
download location is also configurable; see the documentation of any of the
``data_path`` functions for more information.

.. _sample-dataset:

Sample
======
:func:`mne.datasets.sample.data_path`

These data were acquired with the Neuromag
Vectorview system at MGH/HMS/MIT Athinoula A. Martinos Center Biomedical
Imaging. EEG data from a 60-channel electrode cap was acquired simultaneously with
the MEG. The original MRI data set was acquired with a Siemens 1.5 T
Sonata scanner using an MPRAGE sequence.

.. note:: These data are provided solely for the purpose of getting familiar
          with the MNE software. The data should not be used to evaluate the
          performance of the MEG or MRI system employed.

In this experiment, checkerboard patterns were presented to the subject
into the left and right visual field, interspersed by tones to the
left or right ear. The interval between the stimuli was 750 ms. Occasionally
a smiley face was presented at the center of the visual field.
The subject was asked to press a key with the right index finger
as soon as possible after the appearance of the face.

.. table:: Trigger codes for the sample data set.

    =========  =====  ==========================================
    Name              Contents
    =========  =====  ==========================================
    LA         1      Response to left-ear auditory stimulus
    RA         2      Response to right-ear auditory stimulus
    LV         3      Response to left visual field stimulus
    RV         4      Response to right visual field stimulus
    smiley     5      Response to the smiley face
    button     32     Response triggered by the button press
    =========  =====  ==========================================

Contents of the data set
^^^^^^^^^^^^^^^^^^^^^^^^

The sample data set contains two main directories: ``MEG/sample`` (the MEG/EEG
data) and ``subjects/sample`` (the MRI reconstructions).
In addition to subject ``sample``, the MRI surface reconstructions from another
subject, morph, are provided to demonstrate morphing capabilities.

.. table:: Contents of the MEG/sample directory.

    ========================  =====================================================================
    File                      Contents
    ========================  =====================================================================
    sample/audvis_raw.fif     The raw MEG/EEG data
    audvis.ave                A template script for off-line averaging
    auvis.cov                 A template script for the computation of a noise-covariance matrix
    ========================  =====================================================================

.. table:: Overview of the contents of the subjects/sample directory.

    =======================  ======================================================================
    File / directory         Contents
    =======================  ======================================================================
    bem                      Directory for the forward modelling data
    bem/watershed            BEM surface segmentation data computed with the watershed algorithm
    bem/inner_skull.surf     Inner skull surface for BEM
    bem/outer_skull.surf     Outer skull surface for BEM
    bem/outer_skin.surf      Skin surface for BEM
    sample-head.fif          Skin surface in fif format for mne_analyze visualizations
    surf                     Surface reconstructions
    mri/T1                   The T1-weighted MRI data employed in visualizations
    =======================  ======================================================================

The following preprocessing steps have been already accomplished
in the sample data set:

- The MRI surface reconstructions have
  been computed using the FreeSurfer software.

- The BEM surfaces have been created with the watershed algorithm,
  see :ref:`bem_watershed_algorithm`.

The **sample** dataset is distributed with :ref:`fsaverage` for convenience.

Brainstorm
==========
Dataset fetchers for three Brainstorm tutorials are available. Users must agree to the
license terms of these datasets before downloading them. These files are recorded in a CTF 275 system
and are provided in native CTF format (.ds files).

Auditory
^^^^^^^^
:func:`mne.datasets.brainstorm.bst_raw.data_path`.

Details about the data can be found at the Brainstorm `auditory dataset tutorial`_.

.. topic:: Examples

    * :ref:`tut-brainstorm-auditory`: Partially replicates the original Brainstorm tutorial.

Resting state
^^^^^^^^^^^^^
:func:`mne.datasets.brainstorm.bst_resting.data_path`

Details can be found at the Brainstorm `resting state dataset tutorial`_.

.. topic:: Examples

    * :ref:`mne-connectivity:ex-envelope-correlation`

Median nerve
^^^^^^^^^^^^
:func:`mne.datasets.brainstorm.bst_raw.data_path`

Details can be found at the Brainstorm `median nerve dataset tutorial`_.

.. topic:: Examples

    * :ref:`ex-brainstorm-raw`

SPM faces
=========
:func:`mne.datasets.spm_face.data_path`

The `SPM faces dataset`_ contains EEG, MEG and fMRI recordings on face perception.

.. topic:: Examples

    * :ref:`ex-spm-faces` Full pipeline including artifact removal, epochs averaging, forward model computation and source reconstruction using dSPM on the contrast: "faces - scrambled".

EEGBCI motor imagery
====================
:func:`mne.datasets.eegbci.load_data`

The EEGBCI dataset is documented in :footcite:`SchalkEtAl2004`. The data set is
available at PhysioNet :footcite:`GoldbergerEtAl2000`. The dataset contains
64-channel EEG recordings from 109 subjects and 14 runs on each subject in EDF+
format. The recordings were made using the BCI2000 system. To load a subject,
do::

    from mne.io import concatenate_raws, read_raw_edf
    from mne.datasets import eegbci
    raw_fnames = eegbci.load_data(subject, runs)
    raws = [read_raw_edf(f, preload=True) for f in raw_fnames]
    raw = concatenate_raws(raws)

.. topic:: Examples

    * :ref:`ex-decoding-csp-eeg`

.. _somato-dataset:

Somatosensory
=============
:func:`mne.datasets.somato.data_path`

This dataset contains somatosensory data with event-related synchronizations
(ERS) and desynchronizations (ERD).

.. topic:: Examples

    * :ref:`tut-sensors-time-freq`
    * :ref:`ex-inverse-source-power`
    * :ref:`ex-time-freq-global-field-power`

Multimodal
==========
:func:`mne.datasets.multimodal.data_path`

This dataset contains a single subject recorded at Otaniemi (Aalto University)
with auditory, visual, and somatosensory stimuli.

.. topic:: Examples

    * :ref:`ex-io-ave-fiff`

.. _fnirs-motor-dataset:

fNIRS motor
===========
:func:`mne.datasets.fnirs_motor.data_path`

This dataset contains a single subject recorded at Macquarie University.
It has optodes placed over the motor cortex. There are three conditions:

- tapping the left thumb to fingers
- tapping the right thumb to fingers
- a control where nothing happens

The tapping lasts 5 seconds, and there are 30 trials of each condition.

.. topic:: Examples

    * :ref:`tut-fnirs-processing`

High frequency SEF
==================
:func:`mne.datasets.hf_sef.data_path()`

This dataset contains somatosensory evoked fields (median nerve stimulation)
with thousands of epochs. It was recorded with an Elekta TRIUX MEG device at
a sampling frequency of 3 kHz. The dataset is suitable for investigating
high-frequency somatosensory responses. Data from two subjects are included
with MRI images in DICOM format and FreeSurfer reconstructions.

.. topic:: Examples

    * :ref:`high-frequency SEF responses <ex-hf-sef-data>`.

Visual 92 object categories
===========================
:func:`mne.datasets.visual_92_categories.data_path`.

This dataset is recorded using a 306-channel Neuromag vectorview system.

Experiment consisted in the visual presentation of 92 images of human, animal
and inanimate objects either natural or artificial :footcite:`CichyEtAl2014`.
Given the high number of conditions this dataset is well adapted to an approach
based on Representational Similarity Analysis (RSA).

.. topic:: Examples

    * :ref:`Representational Similarity Analysis (RSA) <ex-rsa-noplot>`: Partially replicates the results from :footcite:`CichyEtAl2014`.


mTRF Dataset
============
:func:`mne.datasets.mtrf.data_path`.

This dataset contains 128 channel EEG as well as natural speech stimulus features,
which is also available `here <https://sourceforge.net/projects/aespa/files/>`_.

The experiment consisted of subjects listening to natural speech.
The dataset contains several feature representations of the speech stimulus,
suitable for using to fit continuous regression models of neural activity.
More details and a description of the package can be found in
:footcite:`CrosseEtAl2016`.

.. topic:: Examples

    * :ref:`Receptive Field Estimation and Prediction <ex-receptive-field-mtrf>`: Partially replicates the results from :footcite:`CrosseEtAl2016`.


.. _kiloword-dataset:

Kiloword dataset
================
:func:`mne.datasets.kiloword.data_path`.

This dataset consists of averaged EEG data from 75 subjects performing a
lexical decision task on 960 English words :footcite:`DufauEtAl2015`. The words
are richly annotated, and can be used for e.g. multiple regression estimation
of EEG correlates of printed word processing.


4D Neuroimaging / BTi dataset
=============================
:func:`mne.datasets.phantom_4dbti.data_path`.

This dataset was obtained with a phantom on a 4D Neuroimaging / BTi system at
the MEG center in La Timone hospital in Marseille.

.. topic:: Examples

    * :ref:`tut-phantom-4Dbti`

OPM
===
:func:`mne.datasets.opm.data_path`

OPM data acquired using an Elekta DACQ, simply piping the data into Elekta
magnetometer channels. The FIF files thus appear to come from a TRIUX system
that is only acquiring a small number of magnetometer channels instead of the
whole array.

The OPM ``coil_type`` is custom, requiring a custom ``coil_def.dat``.
The new ``coil_type`` is 9999.

OPM co-registration differs a bit from the typical SQUID-MEG workflow.
No ``-trans.fif`` file is needed for the OPMs, the FIF files include proper
sensor locations in MRI coordinates and no digitization of RPA/LPA/Nasion.
Thus the MEG<->Head coordinate transform is taken to be an identity matrix
(i.e., everything is in MRI coordinates), even though this mis-identifies
the head coordinate frame (which is defined by the relationship of the
LPA, RPA, and Nasion).

Triggers include:

* Median nerve stimulation: trigger value 257.
* Magnetic trigger (in OPM measurement only): trigger value 260.
  1 second before the median nerve stimulation, a magnetic trigger is piped into the MSR.
  This was to be able to check the synchronization between OPMs retrospectively, as each
  sensor runs on an independent clock. Synchronization turned out to be satisfactory.

.. topic:: Examples

    * :ref:`ex-opm-somatosensory`
    * :ref:`ex-opm-resting-state`

The Sleep PolySomnoGraphic Database
===================================
:func:`mne.datasets.sleep_physionet.age.fetch_data`
:func:`mne.datasets.sleep_physionet.temazepam.fetch_data`

The sleep PhysioNet database contains 197 whole-night PolySomnoGraphic sleep
recordings, containing EEG, EOG, chin EMG, and event markers. Some records also
contain respiration and body temperature. Corresponding hypnograms (sleep
patterns) were manually scored by well-trained technicians according to the
Rechtschaffen and Kales manual, and are also available. If you use these
data please cite :footcite:`KempEtAl2000` and :footcite:`GoldbergerEtAl2000`.

.. topic:: Examples

    * :ref:`tut-sleep-stage-classif`

Reference channel noise MEG data set
====================================
:func:`mne.datasets.refmeg_noise.data_path`.

This dataset was obtained with a 4D Neuroimaging / BTi system at
the University Clinic - Erlangen, Germany. There are powerful bursts of
external magnetic noise throughout the recording, which make it a good
example for automatic noise removal techniques.

.. topic:: Examples

    * :ref:`ex-megnoise_processing`

Miscellaneous Datasets
======================
These datasets are used for specific purposes in the documentation and in
general are not useful for separate analyses.

.. _fsaverage:

fsaverage
^^^^^^^^^
:func:`mne.datasets.fetch_fsaverage`

For convenience, we provide a function to separately download and extract the
(or update an existing) fsaverage subject.

.. topic:: Examples

    :ref:`tut-eeg-fsaverage-source-modeling`

Infant template MRIs
^^^^^^^^^^^^^^^^^^^^
:func:`mne.datasets.fetch_infant_template`

This function will download an infant template MRI from
:footcite:`OReillyEtAl2021` along with MNE-specific files.

ECoG Dataset
^^^^^^^^^^^^
:func:`mne.datasets.misc.data_path`. Data exists at ``/ecog/``.

This dataset contains a sample electrocorticography (ECoG) dataset. It includes
two grids of electrodes and ten shaft electrodes with simulated motor data (actual data
pending availability).

.. topic:: Examples

    * :ref:`ex-electrode-pos-2d`: Demonstrates how to project a 3D electrode location onto a 2D image, a common procedure in ECoG analyses.
    * :ref:`tut-ieeg-localize`: Demonstrates how to use a graphical user interface to locate electrode contacts as well as warp them to a common atlas.

sEEG Dataset
^^^^^^^^^^^^
:func:`mne.datasets.misc.data_path`. Data exists at ``/seeg/``.

This dataset contains a sample stereoelectroencephalography (sEEG) dataset.
It includes 21 shaft electrodes during a two-choice movement task on a keyboard.

.. topic:: Examples

    * :ref:`tut-ieeg-localize`: Demonstrates how to use a graphical user interface to locate electrode contacts as well as warp them to a common atlas.
    * :ref:`tut-working-with-seeg`: Demonstrates ways to plot sEEG anatomy and results.

.. _limo-dataset:

LIMO Dataset
^^^^^^^^^^^^
:func:`mne.datasets.limo.load_data`.

In the original LIMO experiment (see :footcite:`RousseletEtAl2010`), participants
performed a
two-alternative forced choice task, discriminating between two face stimuli.
Subjects discriminated the same two faces during the whole experiment.
The critical manipulation consisted of the level of noise added to the
face-stimuli during the task, making the faces more or less discernible to the
observer.

The presented faces varied across a noise-signal (or phase-coherence) continuum
spanning from 0 to 100% in increasing steps of 10%. In other words, faces with
high phase-coherence (e.g., 90%) were easy to identify, while faces with low
phase-coherence (e.g., 10%) were hard to identify and by extension hard to
discriminate.

.. topic:: Examples

    * :ref:`Single trial linear regression analysis with the LIMO dataset
      <ex-limo-data>`: Explores data from a single subject of the LIMO dataset
      and demonstrates how to fit a single trial linear regression using the
      information contained in the metadata of the individual datasets.

.. _erp-core-dataset:

ERP CORE Dataset
^^^^^^^^^^^^^^^^
:func:`mne.datasets.erp_core.data_path`

The original `ERP CORE dataset`_ :footcite:`Kappenman2021` contains data from
40 participants who completed 6 EEG experiments, carefully crafted to evoke
7 well-known event-related potential (ERP) components.

Currently, the MNE-Python ERP CORE dataset only provides data from one
participant (subject ``001``) of the Flankers paradigm, which elicits the
lateralized readiness potential (LRP) and error-related negativity (ERN). The
data provided is **not** the original data from the ERP CORE dataset, but
rather a slightly modified version, designed to demonstrate the Epochs metadata
functionality. For example, we already set the references and montage
correctly, and stored events as Annotations. Data is provided in ``FIFF``
format.

.. topic:: Examples

    * :ref:`tut-autogenerate-metadata`: Learn how to auto-generate
      `~mne.Epochs` metadata, and visualize the error-related negativity (ERN)
      ERP component.

.. _ssvep-dataset:

SSVEP
=====
:func:`mne.datasets.ssvep.data_path`

This is a simple example dataset with frequency tagged visual stimulation:
N=2 participants observed checkerboards patterns inverting with a constant
frequency of either 12.0 Hz of 15.0 Hz. 10 trials of 20.0 s length each.
32 channels wet EEG was recorded.

Data format: BrainVision .eeg/.vhdr/.vmrk files organized according to BIDS
standard.

.. topic:: Examples

    * :ref:`tut-ssvep`

References
==========

.. footbibliography::


.. LINKS

.. _auditory dataset tutorial: https://neuroimage.usc.edu/brainstorm/DatasetAuditory
.. _resting state dataset tutorial: https://neuroimage.usc.edu/brainstorm/DatasetResting
.. _median nerve dataset tutorial: https://neuroimage.usc.edu/brainstorm/DatasetMedianNerveCtf
.. _SPM faces dataset: https://www.fil.ion.ucl.ac.uk/spm/data/mmfaces/
.. _ERP-CORE dataset: https://erpinfo.org/erp-core
:orphan:

.. include:: ../links.inc

.. _mne_matlab:

========================
MNE-MATLAB documentation
========================

.. note:: The MNE MATLAB Toolbox is compatible with Matlab versions 7.0 or later.

Overview
########

The MNE software contains a collection Matlab ``.m``-files to
facilitate interfacing with binary file formats of the MNE software.
The toolbox is located at ``$MNE_ROOT/share/matlab`` . The
names of the MNE Matlab toolbox functions begin either with ``mne_`` or
with ``fiff_`` . When you source the ``mne_setup`` script
as described in :ref:`user_environment`, one of the following actions
takes place:

- If you do not have the Matlab startup.m
  file, it will be created and lines allowing access to the MNE Matlab
  toolbox are added.

- If you have startup.m and it does not have the standard MNE
  Matlab toolbox setup lines, you will be instructed to add them manually.

- If you have startup.m and the standard MNE Matlab toolbox
  setup lines are there, nothing happens.

A summary of the available routines is provided in the `MNE-C manual`_. The
toolbox also contains a set of examples which may be useful starting points
for your own development. The names of these functions start with ``mne_ex``.

.. note::

   The MATLAB function ``fiff_setup_read_raw`` has a significant change. The
   sample numbers now take into account possible initial skip in the file,
   *i.e.*, the time between the start of the data acquisition and the start of
   saving the data to disk. The ``first_samp`` member of the returned structure
   indicates the initial skip in samples. If you want your own routines, which
   assume that initial skip has been removed, perform identically with the
   previous version, subtract ``first_samp`` from the sample numbers you
   specify to ``fiff_read_raw_segment``. Furthermore, ``fiff_setup_read_raw``
   has an optional argument to allow reading of unprocessed MaxShield data
   acquired with the Elekta MEG systems.

.. tabularcolumns:: |p{0.3\linewidth}|p{0.6\linewidth}|
.. _BGBCGHAG:
.. table:: High-level reading routines.

    +--------------------------------+--------------------------------------------------------------+
    | Function                       | Purpose                                                      |
    +================================+==============================================================+
    | fiff_find_evoked               | Find all evoked data sets from a file.                       |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_read_bad_channels         | Read the bad channel list.                                   |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_read_ctf_comp             | Read CTF software gradient compensation data.                |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_read_evoked               | Read evoked-response data.                                   |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_read_evoked_all           | Read all evoked-response data from a file.                   |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_read_meas_info            | Read measurement information.                                |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_read_mri                  | Read an MRI description file.                                |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_read_proj                 | Read signal-space projection data.                           |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_read_raw_segment          | Read a segment of raw data with time limits are specified    |
    |                                | in samples.                                                  |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_read_raw_segment_times    | Read a segment of raw data with time limits specified        |
    |                                | in seconds.                                                  |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_setup_read_raw            | Set up data structures before using fiff_read_raw_segment    |
    |                                | or fiff_read_raw_segment_times.                              |
    +--------------------------------+--------------------------------------------------------------+


.. tabularcolumns:: |p{0.3\linewidth}|p{0.6\linewidth}|
.. table:: Channel selection utilities.

    +--------------------------------+--------------------------------------------------------------+
    | Function                       | Purpose                                                      |
    +================================+==============================================================+
    | fiff_pick_channels             | Create a selector to pick desired channels from data         |
    |                                | according to include and exclude lists.                      |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_pick_channels_evoked      | Pick desired channels from evoked-response data according    |
    |                                | to include and exclude lists.                                |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_pick_info                 | Modify measurement info to include only selected channels.   |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_pick_types                | Create a selector to pick desired channels from data         |
    |                                | according to channel types (MEG, EEG, STIM) in combination   |
    |                                | with include and exclude lists.                              |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_pick_types_evoked         | Pick desired channels from evoked-response data according    |
    |                                | to channel types (MEG, EEG, STIM) in combination with        |
    |                                | include and exclude lists.                                   |
    +--------------------------------+--------------------------------------------------------------+


.. tabularcolumns:: |p{0.3\linewidth}|p{0.6\linewidth}|
.. table:: Coordinate transformation utilities.

    +--------------------------------+--------------------------------------------------------------+
    | Function                       | Purpose                                                      |
    +================================+==============================================================+
    | fiff_invert_transform          | Invert a coordinate transformation structure.                |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_reset_ch_pos              | Reset channel position transformation to the default values  |
    |                                | present in the file.                                         |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_transform_eeg_chs         | Transform electrode positions to another coordinate frame.   |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_transform_meg_chs         | Apply a coordinate transformation to the sensor location     |
    |                                | data to bring the integration points to another coordinate   |
    |                                | frame.                                                       |
    +--------------------------------+--------------------------------------------------------------+


.. tabularcolumns:: |p{0.3\linewidth}|p{0.6\linewidth}|
.. table:: Basic reading routines.

    +--------------------------------+--------------------------------------------------------------+
    | Function                       | Purpose                                                      |
    +================================+==============================================================+
    | fiff_define_constants          | Define a structure which contains the constant relevant      |
    |                                | to fif files.                                                |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_dir_tree_find             | Find nodes of a given type in a directory tree structure.    |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_list_dir_tree             | List a directory tree structure.                             |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_make_dir_tree             | Create a directory tree structure.                           |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_open                      | Open a fif file and create the directory tree structure.     |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_read_named_matrix         | Read a named matrix from a fif file.                         |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_read_tag                  | Read one tag from a fif file.                                |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_read_tag_info             | Read the info of one tag from a fif file.                    |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_split_name_list           | Split a colon-separated list of names into a cell array      |
    |                                | of strings.                                                  |
    +--------------------------------+--------------------------------------------------------------+


.. tabularcolumns:: |p{0.3\linewidth}|p{0.6\linewidth}|
.. table:: Writing routines.

    +--------------------------------+--------------------------------------------------------------+
    | Function                       | Purpose                                                      |
    +================================+==============================================================+
    | fiff_end_block                 | Write a FIFF_END_BLOCK tag.                                  |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_end_file                  | Write the standard closing.                                  |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_start_block               | Write a FIFF_START_BLOCK tag.                                |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_start_file                | Write the appropriate beginning of a file.                   |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_ch_info             | Write a channel information structure.                       |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_coord_trans         | Write a coordinate transformation structure.                 |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_ctf_comp            | Write CTF compensation data.                                 |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_dig_point           | Write one digitizer data point.                              |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_complex             | Write single-precision complex numbers.                      |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_complex_matrix      | Write a single-precision complex matrix.                     |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_double              | Write double-precision floats.                               |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_double_complex      | Write double-precision complex numbers.                      |
    +--------------------------------+--------------------------------------------------------------+
    |fiff_write_double_complex_matrix| Write a double-precision complex matrix.                     |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_double_matrix       | Write a double-precision matrix.                             |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_evoked              | Write an evoked-reponse data file.                           |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_float               | Write single-precision floats.                               |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_float_matrix        | Write a single-precision matrix.                             |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_id                  | Write an id tag.                                             |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_int                 | Write 32-bit integers.                                       |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_int_matrix          | Write a matrix of 32-bit integers.                           |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_name_list           | Write a name list.                                           |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_named_matrix        | Write a named matrix.                                        |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_proj                | Write SSP data.                                              |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_short               | Write 16-bit integers.                                       |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_string              | Write a string.                                              |
    +--------------------------------+--------------------------------------------------------------+


.. tabularcolumns:: |p{0.3\linewidth}|p{0.6\linewidth}|
.. table:: High-level data writing routines.

    +--------------------------------+--------------------------------------------------------------+
    | Function                       | Purpose                                                      |
    +================================+==============================================================+
    | fiff_write_evoked              | Write an evoked-response data file.                          |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_finish_writing_raw        | Write the closing tags to a raw data file.                   |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_start_writing_raw         | Start writing raw data file, *i.e.*, write the measurement   |
    |                                | information.                                                 |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_dig_file            | Write a fif file containing digitization data.               |
    +--------------------------------+--------------------------------------------------------------+
    | fiff_write_raw_buffer          | Write one raw data buffer. This is used after a call to      |
    |                                | fiff_start_writing_raw.                                      |
    +--------------------------------+--------------------------------------------------------------+


.. tabularcolumns:: |p{0.3\linewidth}|p{0.6\linewidth}|
.. table:: Coil definition utilities.

    +--------------------------------+--------------------------------------------------------------+
    | Function                       | Purpose                                                      |
    +================================+==============================================================+
    | mne_add_coil_defs              | Add coil definitions to an array of channel information      |
    |                                | structures.                                                  |
    +--------------------------------+--------------------------------------------------------------+
    | mne_load_coil_def              | Load a coil definition file.                                 |
    +--------------------------------+--------------------------------------------------------------+

.. tabularcolumns:: |p{0.3\linewidth}|p{0.6\linewidth}|
.. table:: Routines for software gradient compensation and signal-space projection.

    +--------------------------------+--------------------------------------------------------------+
    | Function                       | Purpose                                                      |
    +================================+==============================================================+
    | mne_compensate_to              | Apply or remove CTF software gradient compensation from      |
    |                                | evoked-response data.                                        |
    +--------------------------------+--------------------------------------------------------------+
    | mne_get_current_comp           | Get the state of software gradient compensation from         |
    |                                | measurement info.                                            |
    +--------------------------------+--------------------------------------------------------------+
    | mne_make_compensator           | Make a compensation matrix which switches the status of      |
    |                                | CTF software gradient compensation from one state to another.|
    +--------------------------------+--------------------------------------------------------------+
    | mne_make_projector_info        | Create a signal-space projection operator with the           |
    |                                | projection item definitions and cell arrays of channel names |
    |                                | and bad channel names as input.                              |
    +--------------------------------+--------------------------------------------------------------+
    | mne_make_projector_info        | Like mne_make_projector but uses the measurement info        |
    |                                | structure as input.                                          |
    +--------------------------------+--------------------------------------------------------------+
    | mne_set_current_comp           | Change the information about the compensation status in      |
    |                                | measurement info.                                            |
    +--------------------------------+--------------------------------------------------------------+


.. tabularcolumns:: |p{0.3\linewidth}|p{0.6\linewidth}|
.. table:: High-level routines for reading MNE data files.

    +--------------------------------+--------------------------------------------------------------+
    | Function                       | Purpose                                                      |
    +================================+==============================================================+
    | mne_pick_channels_cov          | Pick desired channels from a sensor covariance matrix.       |
    +--------------------------------+--------------------------------------------------------------+
    | mne_pick_channels_forward      | Pick desired channels (rows) from a forward solution.        |
    +--------------------------------+--------------------------------------------------------------+
    | mne_read_bem_surfaces          | Read triangular tessellations of surfaces for                |
    |                                | boundary-element models.                                     |
    +--------------------------------+--------------------------------------------------------------+
    | mne_read_cov                   | Read a covariance matrix.                                    |
    +--------------------------------+--------------------------------------------------------------+
    | mne_read_epoch                 | Read an epoch of data from the output file of mne_epochs2mat.|
    +--------------------------------+--------------------------------------------------------------+
    | mne_read_events                | Read an event list from a fif file produced by               |
    |                                | mne_browse_raw or mne_process_raw.                           |
    +--------------------------------+--------------------------------------------------------------+
    | mne_read_forward_solution      | Read a forward solution from a fif file.                     |
    +--------------------------------+--------------------------------------------------------------+
    | mne_read_inverse_operator      | Read an inverse operator from a fif file.                    |
    +--------------------------------+--------------------------------------------------------------+
    | mne_read_morph_map             | Read an morphing map produced with mne_make_morph_maps.      |
    +--------------------------------+--------------------------------------------------------------+
    | mne_read_noise_cov             | Read a noise-covariance matrix from a fif file.              |
    +--------------------------------+--------------------------------------------------------------+
    | mne_read_source_spaces         | Read source space information from a fif file.               |
    +--------------------------------+--------------------------------------------------------------+


.. tabularcolumns:: |p{0.3\linewidth}|p{0.6\linewidth}|
.. table:: High-level routines for writing MNE data files.

    +--------------------------------+--------------------------------------------------------------+
    | Function                       | Purpose                                                      |
    +================================+==============================================================+
    | mne_write_cov                  | Write a covariance matrix to an open file.                   |
    +--------------------------------+--------------------------------------------------------------+
    | mne_write_cov_file             | Write a complete file containing just a covariance matrix.   |
    +--------------------------------+--------------------------------------------------------------+
    | mne_write_events               | Write a fif format event file compatible with mne_browse_raw |
    |                                | and mne_process_raw.                                         |
    +--------------------------------+--------------------------------------------------------------+
    | mne_write_inverse_sol_stc      | Write stc files containing an inverse solution or other      |
    |                                | dynamic data on the cortical surface.                        |
    +--------------------------------+--------------------------------------------------------------+
    | mne_write_inverse_sol_w        | Write w files containing an inverse solution or other static |
    |                                | data on the cortical surface.                                |
    +--------------------------------+--------------------------------------------------------------+


.. tabularcolumns:: |p{0.3\linewidth}|p{0.6\linewidth}|
.. _BABBDDAI:
.. table:: Routines related to stc, w, and label files.

    +--------------------------------+--------------------------------------------------------------+
    | Function                       | Purpose                                                      |
    +================================+==============================================================+
    | mne_read_stc_file              | Read data from one stc file. The vertex numbering in the     |
    |                                | returned structure will start from 0.                        |
    +--------------------------------+--------------------------------------------------------------+
    | mne_read_stc_file1             | Read data from one stc file. The vertex numbering in the     |
    |                                | returned structure will start from 1.                        |
    +--------------------------------+--------------------------------------------------------------+
    | mne_read_w_file                | Read data from one w file. The vertex numbering in the       |
    |                                | returned structure will start from 0.                        |
    +--------------------------------+--------------------------------------------------------------+
    | mne_read_w_file1               | Read data from one w file. The vertex numbering in the       |
    |                                | returned structure will start from 1.                        |
    +--------------------------------+--------------------------------------------------------------+
    | mne_write_stc_file             | Write a new stc file. It is assumed the the vertex numbering |
    |                                | in the input data structure containing the stc information   |
    |                                | starts from 0.                                               |
    +--------------------------------+--------------------------------------------------------------+
    | mne_write_stc_file1            | Write a new stc file. It is assumed the the vertex numbering |
    |                                | in the input data structure containing the stc information   |
    |                                | starts from 1.                                               |
    +--------------------------------+--------------------------------------------------------------+
    | mne_write_w_file               | Write a new w file. It is assumed the the vertex numbering   |
    |                                | in the input data structure containing the w file            |
    |                                | information starts from 0.                                   |
    +--------------------------------+--------------------------------------------------------------+
    | mne_write_w_file1              | Write a new w file. It is assumed the the vertex numbering   |
    |                                | in the input data structure containing the w file            |
    |                                | information starts from 1.                                   |
    +--------------------------------+--------------------------------------------------------------+
    | mne_read_label_file            | Read a label file (ROI).                                     |
    +--------------------------------+--------------------------------------------------------------+
    | mne_write_label_file           | Write a label file (ROI).                                    |
    +--------------------------------+--------------------------------------------------------------+
    | mne_label_time_courses         | Extract time courses corresponding to a label from an        |
    |                                | stc file.                                                    |
    +--------------------------------+--------------------------------------------------------------+


.. tabularcolumns:: |p{0.3\linewidth}|p{0.6\linewidth}|
.. table:: Routines for reading FreeSurfer surfaces.

    +--------------------------------+--------------------------------------------------------------+
    | Function                       | Purpose                                                      |
    +================================+==============================================================+
    | mne_read_curvature             | Read a curvature file.                                       |
    +--------------------------------+--------------------------------------------------------------+
    | mne_read_surface               | Read one surface, return the vertex locations and            |
    |                                | triangulation info.                                          |
    +--------------------------------+--------------------------------------------------------------+
    | mne_read_surfaces              | Read surfaces corresponding to one or both hemispheres.      |
    |                                | Optionally read curvature information and add derived        |
    |                                | surface data.                                                |
    +--------------------------------+--------------------------------------------------------------+
    | mne_reduce_surface             | Reduce the number of triangles on a surface using the        |
    |                                | reducepatch Matlab function.                                 |
    +--------------------------------+--------------------------------------------------------------+
    | mne_write_surface              | Write a FreeSurfer surface file.                             |
    +--------------------------------+--------------------------------------------------------------+


.. tabularcolumns:: |p{0.3\linewidth}|p{0.6\linewidth}|
.. _BGBEGFBD:
.. table:: Utility functions.

    +--------------------------------+--------------------------------------------------------------+
    | Function                       | Purpose                                                      |
    +================================+==============================================================+
    | mne_block_diag                 | Create a sparse block-diagonal matrix out of a vector.       |
    +--------------------------------+--------------------------------------------------------------+
    | mne_combine_xyz                | Calculate the square sum of the three Cartesian components   |
    |                                | of several vectors listed in one row or column vector.       |
    +--------------------------------+--------------------------------------------------------------+
    | mne_file_name                  | Compose a file name relative to $MNE_ROOT.                   |
    +--------------------------------+--------------------------------------------------------------+
    | mne_find_channel               | Find a channel by name from measurement info.                |
    +--------------------------------+--------------------------------------------------------------+
    | mne_find_source_space_hemi     | Determine whether a given source space belongs to the left   |
    |                                | or right hemisphere.                                         |
    +--------------------------------+--------------------------------------------------------------+
    | mne_fread3                     | Read a three-byte integer.                                   |
    +--------------------------------+--------------------------------------------------------------+
    | mne_fwrite3                    | Write a three-byte integer.                                  |
    +--------------------------------+--------------------------------------------------------------+
    | mne_make_combined_event_file   | Combine data from several trigger channels into one event    |
    |                                | file.                                                        |
    +--------------------------------+--------------------------------------------------------------+
    | mne_omit_first_line            | Omit first line from a multi-line message. This routine is   |
    |                                | useful for formatting error messages.                        |
    +--------------------------------+--------------------------------------------------------------+
    | mne_prepare_inverse_operator   | Prepare inverse operator data for calculating L2             |
    |                                | minimum-norm solutions and dSPM.                             |
    +--------------------------------+--------------------------------------------------------------+
    | mne_setup_toolbox              | Set up the MNE Matlab toolbox.                               |
    +--------------------------------+--------------------------------------------------------------+
    | mne_transform_coordinates      | Transform locations between different coordinate systems.    |
    |                                | This function uses the output file from                      |
    |                                | ``mne_collect_transforms``.                                  |
    +--------------------------------+--------------------------------------------------------------+
    | mne_transpose_named_matrix     | Create a transpose of a named matrix.                        |
    +--------------------------------+--------------------------------------------------------------+
    | mne_transform_source_space_to  | Transform source space data to another coordinate frame.     |
    +--------------------------------+--------------------------------------------------------------+


.. tabularcolumns:: |p{0.3\linewidth}|p{0.6\linewidth}|
.. _BGBEFADJ:
.. table:: Examples demonstrating the use of the toolbox.

    +--------------------------------+--------------------------------------------------------------+
    | Function                       | Purpose                                                      |
    +================================+==============================================================+
    | mne_ex_average_epochs          | Example of averaging epoch data produced by mne_epochs2mat.  |
    +--------------------------------+--------------------------------------------------------------+
    | mne_ex_cancel_noise            | Example of noise cancellation procedures.                    |
    +--------------------------------+--------------------------------------------------------------+
    | mne_ex_compute_inverse         | Example of computing a L2 minimum-norm estimate or a dSPM    |
    |                                | solution.                                                    |
    +--------------------------------+--------------------------------------------------------------+
    | mne_ex_data_sets               | Example of listing evoked-response data sets.                |
    +--------------------------------+--------------------------------------------------------------+
    | mne_ex_evoked_grad_amp         | Compute tangential gradient amplitudes from planar           |
    |                                | gradiometer data.                                            |
    +--------------------------------+--------------------------------------------------------------+
    | mne_ex_read_epochs             | Read epoch data from a raw data file.                        |
    +--------------------------------+--------------------------------------------------------------+
    | mne_ex_read_evoked             | Example of reading evoked-response data.                     |
    +--------------------------------+--------------------------------------------------------------+
    | mne_ex_read_raw                | Example of reading raw data.                                 |
    +--------------------------------+--------------------------------------------------------------+
    | mne_ex_read_write_raw          | Example of processing raw data (read and write).             |
    +--------------------------------+--------------------------------------------------------------+

.. note:: In order for the inverse operator calculation to work correctly with data processed with the Elekta-Neuromag Maxfilter (TM) software, the so-called *processing history* block must be included in data files. Previous versions of the MNE Matlab functions did not copy processing history to files saved. As of March 30, 2009, the Matlab toolbox routines fiff_start_writing_raw and fiff_write_evoked have been enhanced to include these data to the output file as appropriate. If you have older raw data files created in Matlab from input which has been processed Maxfilter, it is necessary to copy the *processing history* block from the original to modified raw data file using the ``mne_copy_processing_history`` utility. The raw data processing programs mne_browse_raw and mne_process_raw have handled copying of the processing history since revision 2.5 of the MNE software.

Some data structures
####################

The MNE Matlab toolbox relies heavily on structures to organize
the data. This section gives detailed information about fields in
the essential data structures employed in the MNE Matlab toolbox.
In the structure definitions, data types referring to other MNE
Matlab toolbox structures are shown in italics. In addition, :ref:`matlab_fif_constants`
lists the values of various FIFF constants defined by fiff_define_constants.m .
The documented structures are:

**tag**

    Contains one tag from the fif file, see :ref:`BGBGIIGD`.

**taginfo**

    Contains the information about one tag, see :ref:`BGBBJBJJ`.

**directory**

    Contains the tag directory as a tree structure, see :ref:`BGBEDHBG`.

**id**

    A fif ID, see :ref:`BGBDAHHJ`.

**named matrix**

    Contains a matrix with names for rows and/or columns, see :ref:`BGBBEDID`.
    A named matrix is used to store, *e.g.*, SSP vectors and forward solutions.

**trans**

    A 4 x 4 coordinate-transformation matrix operating on augmented column
    vectors. Indication of the coordinate frames to which this transformation
    relates is included, see :ref:`BGBDHBIF`.

**dig**

    A Polhemus digitizer data point, see :ref:`BGBHDEDG`.

**coildef**

    The coil definition structure useful for forward calculations and array
    visualization, see :ref:`BGBGBEBH`. For more detailed information on
    coil definitions, see :ref:`coil_geometry_information`.

**ch**

    Channel information structure, see :ref:`BGBIABGD`.

**proj**

    Signal-space projection data, see :ref:`BGBCJHJB`.

**comp**

    Software gradiometer compensation data, see :ref:`BGBJDIFD`.

**measurement info**

    Translation of the FIFFB_MEAS_INFO entity, see :ref:`BGBFHDIJ` and
    :class:`mne.Info`. This data structure is returned by fiff_read_meas_info,
    will not be as complete as :class:`mne.Info`.

**surf**

    Used to represent triangulated surfaces and cortical source spaces, see :ref:`BGBEFJCB`.

**cov**

    Used for storing covariance matrices, see :ref:`BGBJJIED`.

**fwd**

    Forward solution data returned by mne_read_forward_solution ,
    see :ref:`BGBFJIBJ`.

**inv**

    Inverse operator decomposition data returned by mne_read_inverse_operator.
    For more information on inverse operator
    decomposition, see :ref:`minimum_norm_estimates`. For an example on how to
    compute inverse solution using this data, see the sample routine mne_ex_compute_inverse .

.. note:: The MNE Matlab toolbox tries it best to employ vertex numbering starting from 1 as opposed to 0 as recorded in the data files. There are, however, two exceptions where explicit attention to the vertex numbering convention is needed. First, the standard stc and w file reading and writing routines return and    assume zero-based vertex numbering. There are now versions with names ending with '1', which return and assume one-based vertex numbering, see :ref:`BABBDDAI`. Second, the logno field of the channel information in the data files produced by mne_compute_raw_inverse is the zero-based number of the vertex whose source space signal is contained on this channel.


.. tabularcolumns:: |p{0.38\linewidth}|p{0.06\linewidth}|p{0.46\linewidth}|
.. _matlab_fif_constants:
.. table:: FIFF constants.

    +-------------------------------+-------+----------------------------------------------------------+
    | Name                          | Value | Purpose                                                  |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MEG_CH                  | 1     | This is a MEG channel.                                   |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_REF_MEG_CH              | 301   | This a reference MEG channel, located far away from the  |
    |                               |       | head.                                                    |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_EEF_CH                  | 2     | This is an EEG channel.                                  |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MCG_CH                  | 201   | This a MCG channel.                                      |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_STIM_CH                 | 3     | This is a digital trigger channel.                       |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_EOG_CH                  | 202   | This is an EOG channel.                                  |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_EMG_CH                  | 302   | This is an EMG channel.                                  |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_ECG_CH                  | 402   | This is an ECG channel.                                  |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MISC_CH                 | 502   | This is a miscellaneous analog channel.                  |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_RESP_CH                 | 602   | This channel contains respiration monitor output.        |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_COORD_UNKNOWN           | 0     | Unknown coordinate frame.                                |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_COORD_DEVICE            | 1     | The MEG device coordinate frame.                         |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_COORD_ISOTRAK           | 2     | The Polhemus digitizer coordinate frame (does not appear |
    |                               |       | in data files).                                          |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_COORD_HPI               | 3     | HPI coil coordinate frame (does not appear in data       |
    |                               |       | files).                                                  |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_COORD_HEAD              | 4     | The MEG head coordinate frame (Neuromag convention).     |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_COORD_MRI               | 5     | The MRI coordinate frame.                                |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_COORD_MRI_SLICE         | 6     | The coordinate frame of a single MRI slice.              |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_COORD_MRI_DISPLAY       | 7     | The preferred coordinate frame for displaying the MRIs   |
    |                               |       | (used by MRIlab).                                        |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_COORD_DICOM_DEVICE      | 8     | The DICOM coordinate frame (does not appear in files).   |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_COORD_IMAGING_DEVICE    | 9     | A generic imaging device coordinate frame (does not      |
    |                               |       | appear in files).                                        |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_COORD_TUFTS_EEG     | 300   | The Tufts EEG data coordinate frame.                     |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_COORD_CTF_DEVICE    | 1001  | The CTF device coordinate frame (does not appear in      |
    |                               |       | files).                                                  |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_COORD_CTF_HEAD      | 1004  | The CTF/4D head coordinate frame.                        |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_ASPECT_AVERAGE          | 100   | Data aspect: average.                                    |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_ASPECT_STD_ERR          | 101   | Data aspect: standard error of mean.                     |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_ASPECT_SINGLE           | 102   | Single epoch.                                            |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_ASPECT_SUBAVERAGE       | 103   | One subaverage.                                          |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_ASPECT_ALTAVERAGE       | 104   | One alternating (plus-minus) subaverage.                 |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_ASPECT_SAMPLE           | 105   | A sample cut from raw data.                              |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_ASPECT_POWER_DENSITY    | 106   | Power density spectrum.                                  |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_ASPECT_DIPOLE_WAVE      | 200   | The time course of an equivalent current dipole.         |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_BEM_SURF_ID_UNKNOWN     | -1    | Unknown BEM surface.                                     |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_BEM_SURF_ID_BRAIN       | 1     | The inner skull surface                                  |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_BEM_SURF_ID_SKULL       | 3     | The outer skull surface                                  |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_BEM_SURF_ID_HEAD        | 4     | The scalp surface                                        |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_SURF_LEFT_HEMI      | 101   | Left hemisphere cortical surface                         |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_SURF_RIGHT_HEMI     | 102   | Right hemisphere cortical surface                        |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_POINT_CARDINAL          | 1     | Digitization point which is a cardinal landmark a.k.a.   |
    |                               |       | fiducial point                                           |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_POINT_HPI               | 2     | Digitized HPI coil location                              |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_POINT_EEG               | 3     | Digitized EEG electrode location                         |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_POINT_ECG               | 3     | Digitized ECG electrode location                         |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_POINT_EXTRA             | 4     | Additional head surface point                            |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_POINT_LPA               | 1     | Identifier for left auricular landmark                   |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_POINT_NASION            | 2     | Identifier for nasion                                    |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_POINT_RPA               | 3     | Identifier for right auricular landmark                  |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_FIXED_ORI           | 1     | Fixed orientation constraint used in the computation of  |
    |                               |       | a forward solution.                                      |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_FREE_ORI            | 2     | No orientation constraint used in the computation of     |
    |                               |       | a forward solution                                       |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_MEG                 | 1     | Indicates an inverse operator based on MEG only          |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_EEG                 | 2     | Indicates an inverse operator based on EEG only.         |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_MEG_EEG             | 3     | Indicates an inverse operator based on both MEG and EEG. |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_UNKNOWN_COV         | 0     | An unknown covariance matrix                             |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_NOISE_COV           | 1     | Indicates a noise covariance matrix.                     |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_SENSOR_COV          | 1     | Synonym for FIFFV_MNE_NOISE_COV                          |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_SOURCE_COV          | 2     | Indicates a source covariance matrix                     |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_FMRI_PRIOR_COV      | 3     | Indicates a covariance matrix associated with fMRI priors|
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_SIGNAL_COV          | 4     | Indicates the data (signal + noise) covariance matrix    |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_DEPTH_PRIOR_COV     | 5     | Indicates the depth prior (depth weighting) covariance   |
    |                               |       | matrix                                                   |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_MNE_ORIENT_PRIOR_COV    | 6     | Indicates the orientation (loose orientation constrain)  |
    |                               |       | prior covariance matrix                                  |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_PROJ_ITEM_NONE          | 0     | The nature of this projection item is unknown            |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_PROJ_ITEM_FIELD         | 1     | This is projection item is a generic field pattern or    |
    |                               |       | field patterns.                                          |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_PROJ_ITEM_DIP_FIX       | 2     | This projection item is the field of one dipole          |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_PROJ_ITEM_DIP_ROT       | 3     | This projection item corresponds to the fields of three  |
    |                               |       | or two orthogonal dipoles at some location.              |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_PROJ_ITEM_HOMOG_GRAD    | 4     | This projection item contains the homogeneous gradient   |
    |                               |       | fields as seen by the sensor array.                      |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_PROJ_ITEM_HOMOG_FIELD   | 5     | This projection item contains the three homogeneous field|
    |                               |       | components as seen by the sensor array.                  |
    +-------------------------------+-------+----------------------------------------------------------+
    | FIFFV_PROJ_ITEM_EEG_AVREF     | 10    | This projection item corresponds to the average EEG      |
    |                               |       | reference.                                               |
    +-------------------------------+-------+----------------------------------------------------------+

.. _BGBGIIGD:

.. table:: The tag structure.

    =======  ===========  ============================================
    Field    Data type    Description
    =======  ===========  ============================================
    kind     int32        The kind of the data item.
    type     uint32       The data type used to represent the data.
    size     int32        Size of the data in bytes.
    next     int32        Byte offset of the next tag in the file.
    data     various      The data itself.
    =======  ===========  ============================================

.. _BGBBJBJJ:

.. table:: The taginfo structure.

    =======  ===========  ============================================
    Field    Data type    Description
    =======  ===========  ============================================
    kind     double       The kind of the data item.
    type     double       The data type used to represent the data.
    size     double       Size of the data in bytes.
    pos      double       Byte offset to this tag in the file.
    =======  ===========  ============================================

.. _BGBEDHBG:

.. table:: The directory structure.

    ============  ============  ================================================================
    Field         Data type     Description
    ============  ============  ================================================================
    block         double        The block id of this directory node.
    id            id            The unique identifier of this node.
    parent_id     id            The unique identifier of the node this node was derived from.
    nent          double        Number of entries in this node.
    nchild        double        Number of children to this node.
    dir           taginfo       Information about tags in this node.
    children      directory     The children of this node.
    ============  ============  ================================================================

.. _BGBDAHHJ:

.. table:: The id structure.

    ==========  ===========  ============================================================
    Field       Data type    Description
    ==========  ===========  ============================================================
    version     int32        The fif file version (major  < < 16 | minor).
    machid      int32(2)     Unique identifier of the computer this id was created on.
    secs        int32        Time since January 1, 1970 (seconds).
    usecs       int32        Time since January 1, 1970 (microseconds past secs ).
    ==========  ===========  ============================================================

.. _BGBBEDID:

.. table:: The named matrix structure.

    ============  ===========  ======================================================================
    Field         Data type    Description
    ============  ===========  ======================================================================
    nrow          int32        Number of rows.
    ncol          int32        Number of columns.
    row_names     cell(*)      The names of associated with the rows. This member may be empty.
    col_names     cell(*)      The names of associated with the columns. This member may be empty.
    data          various      The matrix data, usually of type single or double.
    ============  ===========  ======================================================================


.. tabularcolumns:: |p{0.2\linewidth}|p{0.2\linewidth}|p{0.55\linewidth}|
.. _BGBDHBIF:
.. table:: The trans structure.

    +---------------------------+-----------+----------------------------------------------------------+
    | Field                     | Data Type | Description                                              |
    +===========================+===========+==========================================================+
    | from                      | int32     | The source coordinate frame, see                         |
    |                           |           | :ref:`matlab_fif_constants`. Look                        |
    |                           |           | for entries starting with FIFFV_COORD or FIFFV_MNE_COORD.|
    +---------------------------+-----------+----------------------------------------------------------+
    | to                        | int32     | The destination coordinate frame.                        |
    +---------------------------+-----------+----------------------------------------------------------+
    | trans                     |double(4,4)| The 4-by-4 coordinate transformation matrix. This        |
    |                           |           | operates from augmented position column vectors given in |
    |                           |           | *from* coordinates to give results in *to* coordinates.  |
    +---------------------------+-----------+----------------------------------------------------------+


.. tabularcolumns:: |p{0.2\linewidth}|p{0.2\linewidth}|p{0.55\linewidth}|
.. _BGBHDEDG:
.. table:: The dig structure.

    +---------------------------+-----------+----------------------------------------------------------+
    | Field                     | Data Type | Description                                              |
    +===========================+===========+==========================================================+
    | kind                      | int32     | The type of digitizing point. Possible values are listed |
    |                           |           | in :ref:`matlab_fif_constants`. Look for entries         |
    |                           |           | starting with FIFF_POINT.                                |
    +---------------------------+-----------+----------------------------------------------------------+
    | ident                     | int32     | Identifier for this point.                               |
    +---------------------------+-----------+----------------------------------------------------------+
    | r                         | single(3) | The location of this point.                              |
    +---------------------------+-----------+----------------------------------------------------------+


.. tabularcolumns:: |p{0.2\linewidth}|p{0.2\linewidth}|p{0.55\linewidth}|
.. _BGBGBEBH:
.. table:: The coildef structure. For more detailed information, see :ref:`coil_geometry_information`.

    +-------------------+-------------------+----------------------------------------------------------+
    | Field             | Data Type         | Description                                              |
    +===================+===================+==========================================================+
    | class             | double            | The coil (or electrode) class.                           |
    +-------------------+-------------------+----------------------------------------------------------+
    | id                | double            | The coil (or electrode) id.                              |
    +-------------------+-------------------+----------------------------------------------------------+
    | accuracy          | double            | Representation accuracy.                                 |
    +-------------------+-------------------+----------------------------------------------------------+
    | num_points        | double            | Number of integration points.                            |
    +-------------------+-------------------+----------------------------------------------------------+
    | size              | double            | Coil size.                                               |
    +-------------------+-------------------+----------------------------------------------------------+
    | baseline          | double            | Coil baseline.                                           |
    +-------------------+-------------------+----------------------------------------------------------+
    | description       | char(*)           | Coil description.                                        |
    +-------------------+-------------------+----------------------------------------------------------+
    | coildefs          | double            | Each row contains the integration point weight, followed |
    |                   | (num_points,7)    | by location [m] and normal.                              |
    +-------------------+-------------------+----------------------------------------------------------+
    | FV                | struct            | Contains the faces and vertices which can be used to     |
    |                   |                   | draw the coil for visualization.                         |
    +-------------------+-------------------+----------------------------------------------------------+


.. tabularcolumns:: |p{0.2\linewidth}|p{0.2\linewidth}|p{0.55\linewidth}|
.. _BGBIABGD:
.. table:: The ch structure.

    +---------------------------+-----------+----------------------------------------------------------+
    | Field                     | Data Type | Description                                              |
    +===========================+===========+==========================================================+
    | scanno                    | int32     | Scanning order number, starting from 1.                  |
    +---------------------------+-----------+----------------------------------------------------------+
    | logno                     | int32     | Logical channel number, conventions in the usage of this |
    |                           |           | number vary.                                             |
    +---------------------------+-----------+----------------------------------------------------------+
    | kind                      | int32     | The channel type (FIFFV_MEG_CH, FIFF_EEG_CH, etc., see   |
    |                           |           | :ref:`matlab_fif_constants` ).                           |
    +---------------------------+-----------+----------------------------------------------------------+
    | range                     | double    | The hardware-oriented part of the calibration factor.    |
    |                           |           | This should be only applied to the continuous raw data.  |
    +---------------------------+-----------+----------------------------------------------------------+
    | cal                       | double    | The calibration factor to bring the channels to physical |
    |                           |           | units.                                                   |
    +---------------------------+-----------+----------------------------------------------------------+
    | loc                       | double(12)| The channel location. The first three numbers indicate   |
    |                           |           | the location [m], followed by the three unit vectors of  |
    |                           |           | the channel-specific coordinate frame. These data contain|
    |                           |           | the values saved in the fif file and should not be       |
    |                           |           | changed. The values are specified in device coordinates  |
    |                           |           | for MEG and in head coordinates for EEG channels,        |
    |                           |           | respectively.                                            |
    +---------------------------+-----------+----------------------------------------------------------+
    | coil_trans                |double(4,4)| Initially, transformation from the channel coordinates   |
    |                           |           | to device coordinates. This transformation is updated by |
    |                           |           | calls to fiff_transform_meg_chs and                      |
    |                           |           | fiff_transform_eeg_chs.                                  |
    +---------------------------+-----------+----------------------------------------------------------+
    | eeg_loc                   | double(6) | The location of the EEG electrode in coord_frame         |
    |                           |           | coordinates. The first three values contain the location |
    |                           |           | of the electrode [m]. If six values are present, the     |
    |                           |           | remaining ones indicate the location of the reference    |
    |                           |           | electrode for this channel.                              |
    +---------------------------+-----------+----------------------------------------------------------+
    | coord_frame               | int32     | Initially, the coordinate frame is FIFFV_COORD_DEVICE    |
    |                           |           | for MEG channels and FIFFV_COORD_HEAD for EEG channels.  |
    +---------------------------+-----------+----------------------------------------------------------+
    | unit                      | int32     | Unit of measurement. Relevant values are: 201 = T/m,     |
    |                           |           | 112 = T, 107 = V, and 202 = Am.                          |
    +---------------------------+-----------+----------------------------------------------------------+
    | unit_mul                  | int32     | The data are given in unit s multiplied by 10unit_mul.   |
    |                           |           | Presently, unit_mul is always zero.                      |
    +---------------------------+-----------+----------------------------------------------------------+
    | ch_name                   | char(*)   | Name of the channel.                                     |
    +---------------------------+-----------+----------------------------------------------------------+
    | coil_def                  | coildef   | The coil definition structure. This is present only if   |
    |                           |           | mne_add_coil_defs has been successfully called.          |
    +---------------------------+-----------+----------------------------------------------------------+


.. tabularcolumns:: |p{0.2\linewidth}|p{0.2\linewidth}|p{0.55\linewidth}|
.. _BGBCJHJB:
.. table:: The proj structure.

    +---------------------------+-----------+----------------------------------------------------------+
    | Field                     | Data Type | Description                                              |
    +===========================+===========+==========================================================+
    | kind                      | int32     | The type of the projection item. Possible values are     |
    |                           |           | listed in :ref:`matlab_fif_constants`. Look for entries  |
    |                           |           | starting with FIFFV_PROJ_ITEM or FIFFV_MNE_PROJ_ITEM.    |
    +---------------------------+-----------+----------------------------------------------------------+
    | active                    | int32     | Is this item active, i.e., applied or about to be        |
    |                           |           | applied to the data.                                     |
    +---------------------------+-----------+----------------------------------------------------------+
    | data                      | named     | The projection vectors. The column names indicate the    |
    |                           | matrix    | names of the channels associated to the elements of the  |
    |                           |           | vectors.                                                 |
    +---------------------------+-----------+----------------------------------------------------------+



.. tabularcolumns:: |p{0.2\linewidth}|p{0.2\linewidth}|p{0.55\linewidth}|
.. _BGBJDIFD:
.. table:: The comp structure.

    +---------------------------+-----------+----------------------------------------------------------+
    | Field                     | Data Type | Description                                              |
    +===========================+===========+==========================================================+
    | ctfkind                   | int32     | The kind of the compensation as stored in file.          |
    +---------------------------+-----------+----------------------------------------------------------+
    | kind                      | int32     | ctfkind mapped into small integer numbers.               |
    +---------------------------+-----------+----------------------------------------------------------+
    | save_calibrated           | logical   | Were the compensation data saved in calibrated form. If  |
    |                           |           | this field is false, the matrix will be decalibrated     |
    |                           |           | using the fields row_cals and col_cals when the          |
    |                           |           | compensation data are saved by the toolbox.              |
    +---------------------------+-----------+----------------------------------------------------------+
    | row_cals                  | double(*) | Calibration factors applied to the rows of the           |
    |                           |           | compensation data matrix when the data were read.        |
    +---------------------------+-----------+----------------------------------------------------------+
    | col_cals                  | double(*) | Calibration factors applied to the columns of the        |
    |                           |           | compensation data matrix when the data were read.        |
    +---------------------------+-----------+----------------------------------------------------------+
    | data                      | named     | The compensation data matrix. The row_names list the     |
    |                           | matrix    | names of the channels to which this compensation applies |
    |                           |           | and the col_names the compensation channels.             |
    +---------------------------+-----------+----------------------------------------------------------+


.. tabularcolumns:: |p{0.2\linewidth}|p{0.2\linewidth}|p{0.55\linewidth}|
.. _BGBFHDIJ:
.. table:: The meas info structure.

    +---------------------------+-----------+----------------------------------------------------------+
    | Field                     | Data Type | Description                                              |
    +===========================+===========+==========================================================+
    | file_id                   | id        | The fif ID of the measurement file.                      |
    +---------------------------+-----------+----------------------------------------------------------+
    | meas_id                   | id        | The ID assigned to this measurement by the acquisition   |
    |                           |           | system or during file conversion.                        |
    +---------------------------+-----------+----------------------------------------------------------+
    | nchan                     | int32     | Number of channels.                                      |
    +---------------------------+-----------+----------------------------------------------------------+
    | sfreq                     | double    | Sampling frequency.                                      |
    +---------------------------+-----------+----------------------------------------------------------+
    | highpass                  | double    | Highpass corner frequency [Hz]. Zero indicates a DC      |
    |                           |           | recording.                                               |
    +---------------------------+-----------+----------------------------------------------------------+
    | lowpass                   | double    | Lowpass corner frequency [Hz].                           |
    +---------------------------+-----------+----------------------------------------------------------+
    | chs                       | ch(nchan) | An array of channel information structures.              |
    +---------------------------+-----------+----------------------------------------------------------+
    | ch_names                  |cell(nchan)| Cell array of channel names.                             |
    +---------------------------+-----------+----------------------------------------------------------+
    | dev_head_t                | trans     | The device to head transformation.                       |
    +---------------------------+-----------+----------------------------------------------------------+
    | ctf_head_t                | trans     | The transformation from 4D/CTF head coordinates to       |
    |                           |           | Neuromag head coordinates. This is only present in       |
    |                           |           | 4D/CTF data.                                             |
    +---------------------------+-----------+----------------------------------------------------------+
    | dev_ctf_t                 | trans     | The transformation from device coordinates to 4D/CTF     |
    |                           |           | head coordinates. This is only present in 4D/CTF data.   |
    +---------------------------+-----------+----------------------------------------------------------+
    | dig                       | dig(*)    | The Polhemus digitization data in head coordinates.      |
    +---------------------------+-----------+----------------------------------------------------------+
    | bads                      | cell(*)   | Bad channel list.                                        |
    +---------------------------+-----------+----------------------------------------------------------+
    | projs                     | proj(*)   | SSP operator data.                                       |
    +---------------------------+-----------+----------------------------------------------------------+
    | comps                     | comp(*)   | Software gradient compensation data.                     |
    +---------------------------+-----------+----------------------------------------------------------+


.. tabularcolumns:: |p{0.2\linewidth}|p{0.2\linewidth}|p{0.55\linewidth}|
.. _BGBEFJCB:

.. table:: The surf structure.

    +---------------------------+-----------+----------------------------------------------------------+
    | Field                     | Data Type | Description                                              |
    +===========================+===========+==========================================================+
    | id                        | int32     | The surface ID.                                          |
    +---------------------------+-----------+----------------------------------------------------------+
    | sigma                     | double    | The electrical conductivity of the compartment bounded by|
    |                           |           | this surface. This field is present in BEM surfaces only.|
    +---------------------------+-----------+----------------------------------------------------------+
    | np                        | int32     | Number of vertices on the surface.                       |
    +---------------------------+-----------+----------------------------------------------------------+
    | ntri                      | int32     | Number of triangles on the surface.                      |
    +---------------------------+-----------+----------------------------------------------------------+
    | coord_frame               | int32     | Coordinate frame in which the locations and orientations |
    |                           |           | are expressed.                                           |
    +---------------------------+-----------+----------------------------------------------------------+
    | rr                        | double    | The vertex locations.                                    |
    |                           | (np,3)    |                                                          |
    +---------------------------+-----------+----------------------------------------------------------+
    | nn                        | double    | The vertex normals. If derived surface data was not      |
    |                           | (np,3)    | requested, this is empty.                                |
    +---------------------------+-----------+----------------------------------------------------------+
    | tris                      | int32     | Vertex numbers of the triangles in counterclockwise      |
    |                           | (ntri,3)  | order as seen from the outside.                          |
    +---------------------------+-----------+----------------------------------------------------------+
    | nuse                      | int32     | Number of active vertices, *i.e.*, vertices included in  |
    |                           |           | a decimated source space.                                |
    +---------------------------+-----------+----------------------------------------------------------+
    | inuse                     | int32(np) | Which vertices are in use.                               |
    +---------------------------+-----------+----------------------------------------------------------+
    | vertno                    |int32(nuse)| Indices of the vertices in use.                          |
    +---------------------------+-----------+----------------------------------------------------------+
    | curv                      | double(np)| Curvature values at the vertices. If curvature           |
    |                           |           | information was not requested, this field is empty or    |
    |                           |           | absent.                                                  |
    +---------------------------+-----------+----------------------------------------------------------+
    | tri_area                  | double    | The triangle areas in m2.If derived surface data was not |
    |                           | (ntri)    | requested, this field will be missing.                   |
    +---------------------------+-----------+----------------------------------------------------------+
    | tri_cent                  | double    | The triangle centroids. If derived surface data was not  |
    |                           | (ntri,3)  | requested, this field will be missing.                   |
    +---------------------------+-----------+----------------------------------------------------------+
    | tri_nn                    | double    | The triangle normals. If derived surface data was not    |
    |                           | (ntri,3)  | requested, this field will be missing.                   |
    +---------------------------+-----------+----------------------------------------------------------+
    | nuse_tri                  | int32     | Number of triangles in use. This is present only if the  |
    |                           |           | surface corresponds to a source space created with the   |
    |                           |           | ``--ico`` option.                                        |
    +---------------------------+-----------+----------------------------------------------------------+
    | use_tris                  | int32     | The vertices of the triangles in use in the complete     |
    |                           | (nuse_tri)| triangulation. This is present only if the surface       |
    |                           |           | corresponds to a source space created with the           |
    |                           |           | ``--ico`` option.                                        |
    +---------------------------+-----------+----------------------------------------------------------+
    | nearest                   | int32(np) | This field is present only if patch information has been |
    |                           |           | computed for a source space. For each vertex in the      |
    |                           |           | triangulation, these values indicate the nearest active  |
    |                           |           | source space vertex.                                     |
    +---------------------------+-----------+----------------------------------------------------------+
    | nearest_dist              | double(np)| This field is present only if patch information has been |
    |                           |           | computed for a source space. For each vertex in the      |
    |                           |           | triangulation, these values indicate the distance to the |
    |                           |           | nearest active source space vertex.                      |
    +---------------------------+-----------+----------------------------------------------------------+
    | dist                      | double    | Distances between vertices on this surface given as a    |
    |                           | (np,np)   | sparse matrix. A zero off-diagonal entry in this matrix  |
    |                           |           | indicates that the corresponding distance has not been   |
    |                           |           | calculated.                                              |
    +---------------------------+-----------+----------------------------------------------------------+
    | dist_limit                | double    | The value given to mne_add_patch_info with the ``--dist``|
    |                           |           | option. This value is presently                          |
    |                           |           | always negative, indicating that only distances between  |
    |                           |           | active source space vertices, as indicated by the vertno |
    |                           |           | field of this structure, have been calculated.           |
    +---------------------------+-----------+----------------------------------------------------------+


.. tabularcolumns:: |p{0.2\linewidth}|p{0.2\linewidth}|p{0.55\linewidth}|
.. _BGBJJIED:

.. table:: The cov structure.

    +---------------------------+-----------+----------------------------------------------------------+
    | Field                     | Data Type | Description                                              |
    +===========================+===========+==========================================================+
    | kind                      | double    | What kind of a covariance matrix (1 = noise covariance,  |
    |                           |           | 2 = source covariance).                                  |
    +---------------------------+-----------+----------------------------------------------------------+
    | diag                      | double    | Is this a diagonal matrix.                               |
    +---------------------------+-----------+----------------------------------------------------------+
    | dim                       | int32     | Dimension of the covariance matrix.                      |
    +---------------------------+-----------+----------------------------------------------------------+
    | names                     | cell(*)   | Names of the channels associated with the entries        |
    |                           |           | (may be empty).                                          |
    +---------------------------+-----------+----------------------------------------------------------+
    | data                      | double    | The covariance matrix. This a double(dim) vector for a   |
    |                           | (dim,dim) | diagonal covariance matrix.                              |
    +---------------------------+-----------+----------------------------------------------------------+
    | projs                     | proj(*)   | The SSP vectors applied to these data.                   |
    +---------------------------+-----------+----------------------------------------------------------+
    | bads                      | cell(*)   | Bad channel names.                                       |
    +---------------------------+-----------+----------------------------------------------------------+
    | nfree                     | int32     | Number of data points used to compute this matrix.       |
    +---------------------------+-----------+----------------------------------------------------------+
    | eig                       |double(dim)| The eigenvalues of the covariance matrix. This field may |
    |                           |           | be empty for a diagonal covariance matrix.               |
    +---------------------------+-----------+----------------------------------------------------------+
    | eigvec                    | double    | The eigenvectors of the covariance matrix.               |
    |                           | (dim,dim) |                                                          |
    +---------------------------+-----------+----------------------------------------------------------+


.. tabularcolumns:: |p{0.2\linewidth}|p{0.2\linewidth}|p{0.55\linewidth}|
.. _BGBFJIBJ:

.. table:: The fwd structure.

    +-------------------------+-------------+----------------------------------------------------------+
    | Field                   | Data Type   | Description                                              |
    +=========================+=============+==========================================================+
    | source_ori              | int32       | Has the solution been computed for the current component |
    |                         |             | normal to the cortex only (1) or all three source        |
    |                         |             | orientations (2).                                        |
    +-------------------------+-------------+----------------------------------------------------------+
    | coord_frame             | int32       | Coordinate frame in which the locations and orientations |
    |                         |             | are expressed.                                           |
    +-------------------------+-------------+----------------------------------------------------------+
    | nsource                 | int32       | Total number of source space points.                     |
    +-------------------------+-------------+----------------------------------------------------------+
    | nchan                   | int32       | Number of channels.                                      |
    +-------------------------+-------------+----------------------------------------------------------+
    | sol                     | named       | The forward solution matrix.                             |
    |                         | matrix      |                                                          |
    +-------------------------+-------------+----------------------------------------------------------+
    | sol_grad                | named       | The derivatives of the forward solution with respect to  |
    |                         | matrix      | the dipole location coordinates.                         |
    |                         |             | This field is present only if the forward solution was   |
    |                         |             | computed with the ``--grad`` option in MNE-C.            |
    +-------------------------+-------------+----------------------------------------------------------+
    | mri_head_t              | trans       | Transformation from the MRI coordinate frame to the      |
    |                         |             | (Neuromag) head coordinate frame.                        |
    +-------------------------+-------------+----------------------------------------------------------+
    | src                     | surf(:)     | The description of the source spaces.                    |
    +-------------------------+-------------+----------------------------------------------------------+
    | source_rr               | double      | The source locations.                                    |
    |                         | (nsource,3) |                                                          |
    +-------------------------+-------------+----------------------------------------------------------+
    | source_nn               | double(:,3) | The source orientations. Number of rows is either        |
    |                         |             | nsource (fixed source orientations) or 3*nsource         |
    |                         |             | (all source orientations).                               |
    +-------------------------+-------------+----------------------------------------------------------+


.. tabularcolumns:: |p{0.2\linewidth}|p{0.2\linewidth}|p{0.55\linewidth}|
.. _BGBIEIJE:

.. table:: The inv structure. Note: The fields proj, whitener, reginv, and noisenorm are filled in by the routine mne_prepare_inverse_operator.

    +---------------------+-------------+----------------------------------------------------------+
    | Field               | Data Type   | Description                                              |
    +=====================+=============+==========================================================+
    | methods             | int32       | Has the solution been computed using MEG data (1), EEG   |
    |                     |             | data (2), or both (3).                                   |
    +---------------------+-------------+----------------------------------------------------------+
    | source_ori          | int32       | Has the solution been computed for the current component |
    |                     |             | normal to the cortex only (1) or all three source        |
    |                     |             | orientations (2).                                        |
    +---------------------+-------------+----------------------------------------------------------+
    | nsource             | int32       | Total number of source space points.                     |
    +---------------------+-------------+----------------------------------------------------------+
    | nchan               | int32       | Number of channels.                                      |
    +---------------------+-------------+----------------------------------------------------------+
    | coord_frame         | int32       | Coordinate frame in which the locations and orientations |
    |                     |             | are expressed.                                           |
    +---------------------+-------------+----------------------------------------------------------+
    | source_nn           | double(:,3) | The source orientations. Number of rows is either        |
    |                     |             | nsource (fixed source orientations) or 3*nsource (all    |
    |                     |             | source orientations).                                    |
    +---------------------+-------------+----------------------------------------------------------+
    | sing                | double      | The singular values, *i.e.*, the diagonal values of      |
    |                     | (nchan)     | :math:`\Lambda`, see :ref:`mne_solution`.                |
    +---------------------+-------------+----------------------------------------------------------+
    | eigen_leads         | double      | The matrix :math:`V`, see :ref:`mne_solution`.           |
    |                     | (:,nchan)   |                                                          |
    +---------------------+-------------+----------------------------------------------------------+
    | eigen_fields        | double      | The matrix :math:`U^\top`, see                           |
    |                     | (nchan,     | :ref:`mne_solution`.                                     |
    |                     | nchan)      |                                                          |
    +---------------------+-------------+----------------------------------------------------------+
    | noise_cov           | cov         | The noise covariance matrix :math:`C`.                   |
    +---------------------+-------------+----------------------------------------------------------+
    | source_cov          | cov         | The source covariance matrix :math:`R`.                  |
    +---------------------+-------------+----------------------------------------------------------+
    | src                 | surf(:)     | The description of the source spaces.                    |
    +---------------------+-------------+----------------------------------------------------------+
    | mri_head_t          | trans       | Transformation from the MRI coordinate frame to the      |
    |                     |             | (Neuromag) head coordinate frame.                        |
    +---------------------+-------------+----------------------------------------------------------+
    | nave                | double      | The number of averages.                                  |
    +---------------------+-------------+----------------------------------------------------------+
    | projs               | proj(:)     | The SSP vectors which were active when the decomposition |
    |                     |             | was computed.                                            |
    +---------------------+-------------+----------------------------------------------------------+
    | proj                | double      | The projection operator computed using projs.            |
    |                     | (nchan)     |                                                          |
    +---------------------+-------------+----------------------------------------------------------+
    | whitener            |             | A sparse matrix containing the noise normalization       |
    |                     |             | factors. Dimension is either nsource (fixed source       |
    |                     |             | orientations) or 3*nsource (all source orientations).    |
    +---------------------+-------------+----------------------------------------------------------+
    | reginv              | double      | The diagonal matrix :math:`\Gamma`, see                  |
    |                     | (nchan)     | :ref:`mne_solution`.                                     |
    +---------------------+-------------+----------------------------------------------------------+
    | noisenorm           | double(:)   | A sparse matrix containing the noise normalization       |
    |                     |             | factors. Dimension is either nsource (fixed source       |
    |                     |             | orientations) or 3*nsource (all source orientations).    |
    +---------------------+-------------+----------------------------------------------------------+


On-line documentation for individual routines
#############################################

Each of the routines listed in Tables :ref:`BGBCGHAG` - :ref:`BGBEFADJ` has on-line documentation accessible by saying ``help`` <*routine name*> in Matlab.
.. include:: ../links.inc

.. _help:

Getting help
^^^^^^^^^^^^

There are several places to obtain help with MNE software tools.

- The `MNE Forum`_ is a good placed to go for both troubleshooting and general
  questions.
- The :ref:`faq` page has some troubleshooting tips, and is a good source of
  general information. There are also some troubleshooting tips built into
  the :ref:`Python <install-python>` and
  :ref:`MNE-Python <standard_instructions>` installation pages (look for the
  |hand-paper| symbols), and some tips related to 3D plotting problems on
  :ref:`the advanced setup page <troubleshoot_3d>`.
- If you want to request new features or if you're confident that you have
  found a bug, please create a new issue on the `GitHub issues page`_.
  When reporting bugs, please try to replicate the bug with the MNE-Python
  :ref:`sample data <sample-dataset>`, and make every effort to simplify your
  example script to only the elements necessary to replicate the bug.


.. toctree::
   :hidden:

   learn_python
   faq
:orphan:

.. include:: ../changes/names.inc

.. _governance-people:

Current steering council and institutional partners
===================================================

Benevolent Dictator for Life
----------------------------

Alexandre Gramfort is the Benevolent Dictator for Life (BDFL)


Steering Council
----------------

* `Adam Li`_
* `Alex Gramfort`_
* `Alex Rockhill`_
* `Britta Westner`_
* `Clemens Brunner`_
* `Daniel McCloy`_
* `Denis Engemann`_
* `Eric Larson`_
* `Guillaume Favelier`_
* `Luke Bloy`_
* `Mainak Jas`_
* `Marijn van Vliet`_
* `Mikołaj Magnuski`_
* `Richard Höchenberger`_
* `Robert Luke`_
* `Stefan Appelhoff`_

Institutional Partners
----------------------

.. include:: ../_includes/institutional-partners.rst
   :start-after: institutional-partners-begin-content


Document history
----------------

https://github.com/mne-tools/mne-python/commits/main/doc/overview/people.rst
.. include:: ../links.inc

.. _faq:

================================
Frequently Asked Questions (FAQ)
================================

.. highlight:: python

General MNE-Python issues
=========================


Help! I can't get Python and MNE-Python working!
------------------------------------------------

Check out our installation instructions for :ref:`Python <install-python>` and
:ref:`MNE-Python <standard_instructions>`.


I still can't get it to work!
-----------------------------

See :ref:`help`.


I can't get PyVista/3D plotting to work under Windows
-----------------------------------------------------

If PyVista plotting in Jupyter Notebooks doesn't work well, using the IPython
magic ``%gui qt`` should `help
<https://github.com/ipython/ipython/issues/10384>`_.

.. code-block:: ipython

   %gui qt

Python runs on macOS extremely slow even on simple commands!
------------------------------------------------------------

Python uses some backends that interfere with the macOS energy saver when
using an IDE such as Spyder or PyCharm. To test it, import ``time`` and run::

    start = time.time(); time.sleep(0.0005); print(time.time() - start)

If it takes several seconds you can either:

- Install the module ``appnope`` and run in your script::

      import appnope
      appnope.nope()

- Change the configuration defaults by running in your terminal:

  .. code-block:: console

      $ defaults write org.python.python NSAppSleepDisabled -bool YES


How do I cite MNE?
------------------

See :ref:`cite`.


I'm not sure how to do *X* analysis step with my *Y* data...
------------------------------------------------------------

Knowing "the right thing" to do with EEG and MEG data is challenging. We use
the `MNE Forum`_ to discuss analysis strategies for different kinds of
data. It's worth searching the archives to see if there have been relevant
discussions in the past, but don't hesitate to ask a new question if the answer
isn't out there already.


I think I found a bug, what do I do?
------------------------------------

When you encounter an error message or unexpected results, it can be hard to
tell whether it happened because of a bug in MNE-Python, a mistake in user
code, a corrupted data file, or irregularities in the data itself. Your first
step when asking for help should be the `MNE Forum`_, not GitHub. This bears
repeating: *the GitHub issue tracker is not for usage help* — it is for
software bugs, feature requests, and improvements to documentation. If you
open an issue that contains only a usage question, we will close the issue and
direct you to the forum. If you're pretty sure the problem you've encountered
is a software bug (not bad data or user error):

- Make sure you're using `the most current version`_. You can check it locally
  at a shell prompt with:

  .. code-block:: console

      $ mne sys_info

  which will also give you version info about important MNE-Python
  dependencies.

- If you're already on the most current version, if possible try using
  :ref:`the latest development version <installing_main>`, as the bug may
  have been fixed already since the latest release. If you can't try the latest
  development version, search the GitHub issues page to see if the problem has
  already been reported and/or fixed.

- Try to replicate the problem with one of the :ref:`MNE sample datasets
  <datasets>`. If you can't replicate it with a built-in dataset, provide a
  link to a small, anonymized portion of your data that does yield the error.

If the problem persists, `open a new issue
<https://github.com/mne-tools/mne-python/issues/new?template=bug_report.md>`__
and include the *smallest possible* code sample that replicates the error
you're seeing. Paste the code sample into the issue, with a line containing
three backticks (\`\`\`) above and below the lines of code. This
`minimal working example`_ should be self-contained, which means that
MNE-Python contributors should be able to copy and paste the provided snippet
and replicate the bug on their own computers.

Why is it dangerous to "pickle" my MNE-Python objects and data for later use?
-----------------------------------------------------------------------------

`Pickling <https://docs.python.org/3/library/pickle.html>`_ data and MNE-Python
objects for later use can be tempting due to its simplicity and generality, but
it is usually not the best option. Pickling is not designed for stable
persistence, and it is likely that you will not be able to read your data in
the not-too-distant future. For details, see:

- http://www.benfrederickson.com/dont-pickle-your-data/
- https://stackoverflow.com/questions/21752259/python-why-pickle

MNE-Python is designed to provide its own file saving formats (often based on
the FIF standard) for its objects usually via a ``save`` method or ``write_*``
method, e.g. :func:`mne.io.Raw.save`, :func:`mne.Epochs.save`,
:func:`mne.write_evokeds`, :func:`mne.SourceEstimate.save`. If you have some
data that you want to save but can't figure out how, post to the `MNE Forum`_
or to the `GitHub issues page`_.

If you want to write your own data to disk (e.g., subject behavioral scores),
we strongly recommend using h5io_, which is based on the `HDF5 format
<https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_ and h5py_, to save
data in a fast, future-compatible, standard format.


I downloaded a dataset once, but MNE-Python is asking to download it again. Why?
--------------------------------------------------------------------------------

The default location for the MNE-sample data is ``~/mne_data``. If you
downloaded data and an example asks you whether to download it again, make sure
the data reside in the examples directory and that you run the script from its
current directory:

.. code-block:: console

  $ cd examples/preprocessing

Then in Python you can do::

  In [1]: %run plot_find_ecg_artifacts.py

See :ref:`datasets` for a list of all available datasets and some advanced
configuration options, e.g. to specify a custom location for storing the
datasets.


.. _faq_cpu:

A function uses multiple CPU cores even though I didn't tell it to. Why?
------------------------------------------------------------------------

Ordinarily in MNE-python the ``parallel`` module is used to deploy multiple
cores via the ``n_jobs`` variable. However, functions like
:func:`mne.preprocessing.maxwell_filter` that use :mod:`scipy.linalg` do not
have an ``n_jobs`` flag but may still use multiple cores. This is because
:mod:`scipy.linalg` is built with linear algebra libraries that natively
support multithreading:

- `OpenBLAS <http://www.openblas.net/>`_
- `Intel Math Kernel Library (MKL) <https://software.intel.com/en-us/mkl>`_,
  which uses `OpenMP <https://www.openmp.org/>`_

To control how many cores are used for linear-algebra-heavy functions like
:func:`mne.preprocessing.maxwell_filter`, you can set the ``OMP_NUM_THREADS``
or ``OPENBLAS_NUM_THREADS`` environment variable to the desired number of cores
for MKL or OpenBLAS, respectively. This can be done before running Python, or
inside Python you can achieve the same effect by, e.g.::

    >>> import os
    >>> num_cpu = '4' # Set as a string
    >>> os.environ['OMP_NUM_THREADS'] = num_cpu

This must be done *before* running linear algebra functions; subsequent
changes in the same Python session will have no effect.


I have a mystery FIF file, how do I read it?
--------------------------------------------

The :func:`mne.what` function can be called on any :file:`.fif` file to
identify the kind of data contained in the file. This will help you determine
whether to use :func:`mne.read_cov`, :func:`mne.read_epochs`,
:func:`mne.read_evokeds`, etc. There is also a corresponding command line tool
:ref:`mne what`:

.. code-block:: console

    $ mne what sample_audvis_eog-eve.fif
    events


.. _resampling-and-decimating:

Resampling and decimating data
==============================

What are all these options for resampling, decimating, and binning data?
------------------------------------------------------------------------

There are many functions in MNE-Python for changing the effective sampling rate
of data. We'll discuss some major ones here, with some of their implications:

- :func:`mne.io.Raw.resample` is used to resample (typically downsample) raw
  data. Resampling is the two-step process of applying a low-pass FIR filter
  and subselecting samples from the data.

  Using this function to resample data before forming :class:`mne.Epochs`
  for final analysis is generally discouraged because doing so effectively
  loses precision of (and jitters) the event timings, see
  `this gist <https://gist.github.com/larsoner/01642cb3789992fbca59>`_ as
  a demonstration. However, resampling raw data can be useful for
  (at least):

    - Computing projectors in low- or band-passed data
    - Exploring data

- :func:`mne.preprocessing.ICA.fit` decimates data without low-passing,
  but is only used for fitting a statistical model to the data.

- :func:`mne.Epochs.decimate`, which does the same thing as the
  ``decim`` parameter in the :class:`mne.Epochs` constructor, sub-selects every
  :math:`N^{th}` sample before and after each event. This should only be
  used when the raw data have been sufficiently low-passed e.g. by
  :func:`mne.io.Raw.filter` to avoid aliasing artifacts.

- :func:`mne.Epochs.resample`, :func:`mne.Evoked.resample`, and
  :func:`mne.SourceEstimate.resample` all resample data.
  This process avoids potential aliasing artifacts because the
  resampling process applies a low-pass filter. However, this filtering
  introduces edge artifacts. Edge artifacts also exist when using
  :func:`mne.io.Raw.resample`, but there the edge artifacts are constrained
  to two times: the start and end of the recording. With these three methods,
  edge artifacts are introduced to the start and end of every epoch
  of data (or the start and end of the :class:`mne.Evoked` or
  :class:`mne.SourceEstimate` data), which often has a more pronounced
  effect on the data.

- :func:`mne.SourceEstimate.bin` can be used to decimate, with or without
  "binning" (averaging across data points). This is equivalent to applying
  a moving-average (boxcar) filter to the data and decimating. A boxcar in
  time is a `sinc <https://en.wikipedia.org/wiki/Sinc_function>`_ in
  frequency, so this acts as a simplistic, non-ideal low-pass filter;
  this will reduce but not eliminate aliasing if data were not sufficiently
  low-passed. In the case where the "filter" or bin-width is a single sample
  (i.e., an impulse) this operation simplifies to decimation without filtering.


Resampling raw data is taking forever! What do I do?
----------------------------------------------------

:func:`mne.io.Raw.resample` has a parameter ``npad=='auto'``. This is the
default, but if you've changed it you could try changing it back to ``'auto'``,
it might help.

If you have an NVIDIA GPU you could also try using :ref:`CUDA`, which can
sometimes speed up filtering and resampling operations by an order of
magnitude.


Forward and Inverse Solution
============================


How should I regularize the covariance matrix?
----------------------------------------------

The estimated covariance can be numerically unstable and tends to induce
correlations between estimated source amplitudes and the number of samples
available. It is thus suggested to regularize the noise covariance
matrix (see :ref:`cov_regularization_math`), especially if only few samples
are available. Unfortunately it is not easy to tell the effective number of
samples, hence, to choose the appropriate regularization. In MNE-Python,
regularization is done using advanced regularization methods described in
:footcite:`EngemannGramfort2015`. For this the 'auto' option can be used. With
this option cross-validation will be used to learn the optimal regularization::

    >>> import mne
    >>> epochs = mne.read_epochs(epochs_path) # doctest: +SKIP
    >>> cov = mne.compute_covariance(epochs, tmax=0., method='auto') # doctest: +SKIP

This procedure evaluates the noise covariance quantitatively by how well it
whitens the data using the negative log-likelihood of unseen data. The final
result can also be visually inspected. Under the assumption that the baseline
does not contain a systematic signal (time-locked to the event of interest),
the whitened baseline signal should be follow a multivariate Gaussian
distribution, i.e., whitened baseline signals should be between -1.96 and 1.96
at a given time sample. Based on the same reasoning, the expected value for the
:term:`global field power` (GFP) is 1 (calculation of the :term:`GFP`
should take into account the true degrees of freedom, e.g. ``ddof=3`` with 2
active SSP vectors)::

    >>> evoked = epochs.average() # doctest: +SKIP
    >>> evoked.plot_white(cov) # doctest: +SKIP

This plot displays both, the whitened evoked signals for each channels and the
whitened :term:`GFP`. The numbers in the :term:`GFP` panel represent the
estimated rank of the data, which amounts to the effective degrees of freedom
by which the squared sum across sensors is divided when computing the whitened
:term:`GFP`. The whitened :term:`GFP` also helps detecting spurious late evoked
components which can be the consequence of over- or under-regularization.

Note that if data have been processed using signal space separation (SSS)
:footcite:`TauluEtAl2005`, gradiometers and magnetometers will be displayed
jointly because both are reconstructed from the same SSS basis vectors with the
same numerical rank. This also implies that both sensor types are not any
longer linearly independent.

These methods for evaluation can be used to assess model violations. Additional
introductory materials can be found `here
<https://speakerdeck.com/dengemann/eeg-sensor-covariance-using-cross-validation>`_.

For expert use cases or debugging the alternative estimators can also be
compared::

    >>> covs = mne.compute_covariance(epochs, tmax=0., method='auto', return_estimators=True) # doctest: +SKIP
    >>> evoked = epochs.average() # doctest: +SKIP
    >>> evoked.plot_white(covs) # doctest: +SKIP

This will plot the whitened evoked for the optimal estimator and display the
:term:`GFP` for all estimators as separate lines in the related panel.


.. _faq_watershed_bem_meshes:

My watershed BEM meshes look incorrect
--------------------------------------

After using :ref:`mne watershed_bem` or :func:`mne.bem.make_watershed_bem`
you might find that the BEM meshes for the brain, inner skull, outer skull,
and/or scalp surfaces do not look correct in :func:`mne.viz.plot_alignment`
and :func:`mne.viz.plot_bem`.

MNE relies on FreeSurfer's mri_watershed_ to compute the BEM meshes.
Freesurfer's watershed bem strategy is to:

1. Compute the outer skin (scalp) surface
2. Shrink outer skin inward make the "outer skull"
3. Compute brain surface
4. Expand brain surface outward to make the "inner skull"

A common problem is to see:

    the surface inner skull is not completely inside surface outer skull

When looking at the meshes, the inner skull surface (expanded brain surface)
will have defects, and these defects will protrude into the outer skull surface
(shrunken scalp surface). In these cases, you can try (in rough ascending
order of difficulty):

.. highlight:: console

1. Changing the ``--preflood`` / ``-p`` parameter in
   :ref:`mne watershed_bem`.
2. Changing the ``--atlas`` and ``--gcaatlas`` options of
   :ref:`mne watershed_bem`.
3. Manually editing the meshes (see :ref:`this tutorial <tut-fix-meshes>`).
4. Manually running mri_watershed_ with various FreeSurfer flags (e.g.,
   ``-less`` to fix the output).
5. Going farther back in your Freesurfer pipeline to fix the problem.
   In particular, ``mri/brainmask.mgz`` could be incorrectly generated by the
   autorecon1_ step and contain some dura and/or skull within the brain mask.
   You can check by using freeview_ or some other MRI-viewing tool.

   - Consult the Freesurfer docs on `fixing errors
     <https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/TroubleshootingDataV6.0#Fixingerrors>`__.
   - Try tweaking the mri_normalize_ parameters `via xopts
     <https://www.mail-archive.com/freesurfer@nmr.mgh.harvard.edu/msg20991.html>`__,
     e.g.::

         $ mri_normalize -mprage -b 20 -n 5

   - Try `manually setting the control points and/or using -gentle
     <https://www.mail-archive.com/freesurfer@nmr.mgh.harvard.edu/msg11658.html>`__.
   - Examine the talairach transformation to see if it's not quite right,
     and if it's not, `adjust it manually
     <https://surfer.nmr.mgh.harvard.edu/fswiki/Edits>`__.
   - Search the `FreeSurfer listserv`_ for other ideas

   It can be helpful to run ``recon_all -autorecon1 -xopts xopts.txt`` in a
   clean directory first to see if this fixes everything, and, if not, then
   resorting to manual control point setting and/or talairach adjustment.
   Once everything looks good at the end of ``-autorecon1``, you can then run
   :ref:`mne watershed_bem` to see if the output is good. Once it is
   (and once brainmask.mgz is correct), you can then proceed with
   ``recon_all -autorecon2`` and ``recon_all -autorecon3`` to effectively
   complete all ``recon_all`` steps.

.. highlight:: python


References
----------

.. footbibliography::

.. LINKS

.. _`the most current version`: https://github.com/mne-tools/mne-python/releases/latest
.. _`minimal working example`: https://en.wikipedia.org/wiki/Minimal_Working_Example
.. _mri_watershed: https://freesurfer.net/fswiki/mri_watershed
.. _mri_normalize: https://surfer.nmr.mgh.harvard.edu/fswiki/mri_normalize
.. _freeview: https://surfer.nmr.mgh.harvard.edu/fswiki/FreeviewGuide/FreeviewIntroduction
.. _`FreeSurfer listserv`: https://www.mail-archive.com/freesurfer@nmr.mgh.harvard.edu/
.. _autorecon1: https://surfer.nmr.mgh.harvard.edu/fswiki/ReconAllDevTable
.. include:: ../links.inc

.. _design_philosophy:

Design philosophy
=================

Interactive versus scripted analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MNE-Python has some great interactive plotting abilities that can help you
explore your data, and there are a few GUI-like interactive plotting commands
(like browsing through the raw data and clicking to mark bad channels, or
click-and-dragging to annotate bad temporal spans). But in general it is not
possible to use MNE-Python to mouse-click your way to a finished, publishable
analysis. MNE-Python works best when you assemble your analysis pipeline into
one or more Python scripts. On the plus side, your scripts act as a record of
everything you did in your analysis, making it easy to tweak your analysis
later and/or share it with others (including your future self).


Integration with the scientific python stack
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MNE-Python also integrates well with other standard scientific Python
libraries. For example, MNE-Python objects underlyingly store their data in
NumPy arrays, making it easy to apply custom algorithms or pass your data into
one of `scikit-learn's <scikit-learn_>`_ machine learning pipelines.
MNE-Python's 2-D plotting functions also return `matplotlib`_
:class:`~matplotlib.figure.Figure` objects, and the 3D plotting functions
return :class:`mne.viz.Figure3D` classes with a ``.plotter`` attribute
pointing to :class:`pyvista.Plotter` instances,
so you can customize your MNE-Python plots using any
of matplotlib or PyVista's plotting commands. The intent is that MNE-Python
will get most neuroscientists 90% of the way to their desired analysis goal,
and other packages can get them over the finish line.


Submodule-based organization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A useful-to-know organizing principle is that MNE-Python objects and functions
are separated into submodules. This can help you discover related functions if
you're using an editor that supports tab-completion. For example, you can type
:samp:`mne.preprocessing.{<TAB>}` to see all the functions in the preprocessing
submodule; similarly for visualization functions (:mod:`mne.viz`), functions
for reading and writing data (:mod:`mne.io`), statistics (:mod:`mne.stats`),
etc.  This also helps save keystrokes — instead of::

    import mne
    mne.preprocessing.eog.peak_finder(...)
    mne.preprocessing.eog.find_eog_events(...)
    mne.preprocessing.eog.create_eog_epochs(...)

you can import submodules directly, and use just the submodule name to access
its functions::

    from mne.preprocessing import eog
    eog.peak_finder(...)
    eog.find_eog_events(...)
    eog.create_eog_epochs(...)


(Mostly) unified API
^^^^^^^^^^^^^^^^^^^^

Whenever possible, we've tried to provide a unified API for the different data
classes. For example, the :class:`~mne.io.Raw`, :class:`~mne.Epochs`,
:class:`~mne.Evoked`, and :class:`~mne.SourceEstimate` classes all have a
``plot()`` method that can typically be called with no parameters specified and
still yield an informative plot of the data. Similarly, they all have the
methods ``copy()``, ``crop()``, ``resample()`` and ``save()`` with similar or
identical method signatures. The sensor-level classes also all have an ``info``
attribute containing an :class:`~mne.Info` object, which keeps track of channel
names and types, applied filters, projectors, etc. See :ref:`tut-info-class`
for more info.


.. _sect-meth-chain:

In-place operation
^^^^^^^^^^^^^^^^^^

Because neuroimaging datasets can be quite large, MNE-Python tries very hard to
avoid making unnecessary copies of your data behind-the-scenes. To further
improve memory efficiency, many object methods operate in-place (and silently
return their object to allow `method chaining`_). In-place operation may lead
you to frequent use of the ``copy()`` method during interactive, exploratory
analysis — so you can try out different preprocessing approaches or parameter
settings without having to re-load the data each time — but it can also be a
big memory-saver when applying a finished script to dozens of subjects' worth
of data.



.. LINKS

.. _`method chaining`: https://en.wikipedia.org/wiki/Method_chaining
.. _migrating:

Migrating from other analysis software
======================================

Here we offer some tips on how to migrate from other analysis software.

EEGLAB
^^^^^^

To read in data exported from EEGLAB, MNE-Python includes an :file:`.edf`
reader :func:`mne.io.read_raw_edf` and a ``set`` file reader. To read in
``set`` files containing ``raw`` data, use :func:`mne.io.read_raw_eeglab` and
to read in ``set`` files containing ``epochs`` data, use
:func:`mne.read_epochs_eeglab`.

This table summarizes the equivalent EEGLAB and MNE-Python code for some of the
most common analysis tasks. For the sake of clarity, the table below assumes
the following variables exist: the file name ``fname``, time interval of the
epochs ``tmin`` and ``tmax``, and the experimental conditions ``cond1`` and
``cond2``. The variables ``l_freq`` and ``h_freq`` are the frequencies (in Hz)
below which and above which to filter out data.

.. cssclass:: table-bordered
.. rst-class:: midvalign

+---------------------+----------------------------------------------------------+--------------------------------------------------------------------------------------------------+
| Processing step     | EEGLAB function                                          | MNE-Python                                                                                       |
+=====================+==========================================================+==================================================================================================+
| Get started         | | ``addpath(...);``                                      | | :mod:`import mne <mne>`                                                                        |
|                     | | ``eeglab;``                                            | | :mod:`from mne import io, <mne.io>` :class:`~mne.Epochs`                                       |
|                     | |                                                        | | :mod:`from mne.preprocessing <mne.preprocessing>` :class:`import ICA <mne.preprocessing.ICA>`  |
+---------------------+----------------------------------------------------------+--------------------------------------------------------------------------------------------------+
| Import data         | | ``EEG = pop_fileio(fname);``                           | | :func:`raw = io.read_raw_fif(fname) <mne.io.read_raw_fif>`                                     |
|                     | |                                                        | | :func:`raw = io.read_raw_edf(fname) <mne.io.read_raw_edf>`                                     |
|                     | |                                                        | | :func:`raw = io.read_raw_eeglab(fname) <mne.io.read_raw_eeglab>` ``(set file)``                |
|                     | |                                                        | |                                                                                                |
+---------------------+----------------------------------------------------------+--------------------------------------------------------------------------------------------------+
| Filter data         | | ``EEG = pop_eegfiltnew(EEG, l_freq, h_freq);``         | | :func:`raw.filter(l_freq, h_freq) <mne.io.Raw.filter>`                                         |
+---------------------+----------------------------------------------------------+--------------------------------------------------------------------------------------------------+
| Common Average      | | ``EEG= pop_averef;``                                   | | :func:`raw.set_eeg_reference("average") <mne.io.Raw.set_eeg_reference>`                        |
| referencing         | |                                                        | |                                                                                                |
+---------------------+----------------------------------------------------------+--------------------------------------------------------------------------------------------------+
| Remove channels     | | ``pop_select.m``                                       | | :func:`raw.drop_channels() <mne.io.Raw.drop_channels>`                                         |
|                     | |                                                        | |                                                                                                |
+---------------------+----------------------------------------------------------+--------------------------------------------------------------------------------------------------+
| Run ICA             | | ``EEG = pop_runica(EEG, 'pca', n);``                   | | :func:`ica.fit(raw) <mne.preprocessing.ICA.fit>`                                               |
|                     | |                                                        | |                                                                                                |
|                     | | ``EEG = pop_binica(EEG, 'pca', n);``                   | | :func:`mne.preprocessing.infomax`                                                              |
+---------------------+----------------------------------------------------------+--------------------------------------------------------------------------------------------------+
| Plot ICA properties | | ``pop_compprop( EEG, comp_num, winhandle);``           | | :func:`ica.plot_properties(raw, picks) <mne.preprocessing.ICA.plot_properties>`                |
+---------------------+----------------------------------------------------------+--------------------------------------------------------------------------------------------------+
| Plot ICA components | | ``compheads()``                                        | | :func:`ica.plot_components(raw, picks) <mne.preprocessing.ICA.plot_components>`                |
+---------------------+----------------------------------------------------------+--------------------------------------------------------------------------------------------------+
| Exclude components  | | ``pop_selectcomps()``                                  | | ``ica.exclude = list_of_components_to_exclude``                                                |
+---------------------+----------------------------------------------------------+--------------------------------------------------------------------------------------------------+
| Epoch data          | | ``event_id = {'cond1', 'cond2'};``                     | | :func:`events = mne.find_events(raw) <mne.find_events>`                                        |
|                     | | ``Epochs = pop_epochs(EEG, event_id, [tmin, tmax]);``  | | :class:`event_id = dict(cond1=32, cond2=64) <dict>`                                            |
|                     | |                                                        | | :class:`epochs = Epochs(raw, events, event_id, tmin, tmax) <mne.Epochs>`                       |
+---------------------+----------------------------------------------------------+--------------------------------------------------------------------------------------------------+
| Selecting epochs    | | ``Epochs = pop_epochs(EEG_epochs, {cond2});``          | | :class:`epochs[cond2] <mne.Epochs>`                                                            |
+---------------------+----------------------------------------------------------+--------------------------------------------------------------------------------------------------+
| ERP butterfly plot  | | ``pop_timtopo(EEG_epochs, ...);``                      | | :meth:`evoked = epochs[cond2].average() <mne.Epochs.average>`                                  |
|                     | |                                                        | | :func:`evoked.plot() <mne.Evoked.plot>`                                                        |
|                     | |                                                        | | :func:`evoked.plot_joint() <mne.Evoked.plot_joint>`                                            |
+---------------------+----------------------------------------------------------+--------------------------------------------------------------------------------------------------+
| Contrast ERPs       | | ``pop_compareerps(EEG_epochs1, EEG_epochs2);``         | | :func:`mne.combine_evoked([evoked1, -evoked2], weights='equal').plot() <mne.combine_evoked>`   |
|                     | |                                                        | | :func:`mne.viz.plot_compare_evokeds([evoked1, evoked2]) <mne.viz.plot_compare_evokeds>`        |
+---------------------+----------------------------------------------------------+--------------------------------------------------------------------------------------------------+
| Save data           | | ``EEG = pop_saveset(EEG, fname);``                     | | :func:`raw.save(fname) <mne.io.Raw.save>`                                                      |
|                     | |                                                        | | :func:`epochs.save(fname) <mne.Epochs.save>`                                                   |
|                     | |                                                        | | :func:`evoked.save(fname) <mne.Evoked.save>`                                                   |
+---------------------+----------------------------------------------------------+--------------------------------------------------------------------------------------------------+

Potential pitfalls
~~~~~~~~~~~~~~~~~~

- Many of the MNE-Python objects have methods that operate in-place to save
  memory (i.e., the data in the :class:`~mne.io.Raw` object is changed when you
  call :meth:`raw.filter(lfreq, hfreq) <mne.io.Raw.filter>`). If you do not
  want this, it is always possible to first call the object's
  :meth:`~mne.io.Raw.copy` method (e.g., ``filtered_raw =
  raw.copy().filter(lfreq, hfreq)``). In addition, some MNE-Python functions
  have a boolean ``copy`` parameter that achieves the same purpose.

- The concept of channel types is critical in MNE because it supports analysis
  of multimodal data (e.g., EEG, MEG, EOG, Stim channel, etc) whereas most
  EEGLAB functions assume all channels are of the same type (EEG). To restrict
  channels to a single type, see :func:`mne.pick_types`, :meth:`raw.pick_types
  <mne.io.Raw.pick_types>`, :meth:`epochs.pick_types <mne.Epochs.pick_types>`,
  :meth:`evoked.pick_types <mne.Evoked.pick_types>`, etc.
MNE-Python Development
======================

.. NOTE: this first section (up until "overview of contribution process") is
   basically a copy/paste of CONTRIBUTING.md from the repository root, with one
   sentence deleted to avoid self-referential linking. Changes made here should
   be mirrored there, and vice-versa.

MNE-Python is maintained by a community of scientists and research labs. The
project accepts contributions in the form of bug reports, fixes, feature
additions, and documentation improvements (including typo corrections). The
best way to start contributing is by `opening an issue`_ on our GitHub page to
discuss ideas for changes or enhancements, or to tell us about behavior that
you think might be a bug. For *general troubleshooting* or *usage questions*,
please consider posting your questions on our `MNE Forum`_.

Users and contributors to MNE-Python are expected to follow our
`code of conduct`_.

The `contributing guide`_ has details on the preferred contribution workflow
and the recommended system configuration for a smooth contribution/development
experience.

.. _`opening an issue`: https://github.com/mne-tools/mne-python/issues/new/choose
.. _`MNE Forum`: https://mne.discourse.group
.. _`code of conduct`: https://github.com/mne-tools/.github/blob/main/CODE_OF_CONDUCT.md
.. _`contributing guide`: https://mne.tools/dev/install/contributing.html

.. toctree::
   :hidden:

   ../install/contributing
   ../whats_new
   roadmap
   governance
:orphan:

Supported channel types
=======================

.. NOTE: part of this file is included in doc/overview/implementation.rst.
   Changes here are reflected there. If you want to link to this content, link
   to :ref:`channel-types` to link to that section of the implementation.rst
   page. The next line is a target for :start-after: so we can omit the title
   from the include:
   channel-types-begin-content

Channel types are represented in MNE-Python with shortened or abbreviated
names. This page lists all supported channel types, their abbreviated names,
and the measurement unit used to represent data of that type. Where channel
types occur in two or more sub-types, the sub-type abbreviations are given in
parentheses. More information about measurement units is given in the
:ref:`units` section.

.. NOTE: To include only the table, here's a different target for :start-after:
   channel-types-begin-table

.. cssclass:: table-bordered
.. rst-class:: midvalign

=============  ========================================= =================
Channel type    Description                              Measurement unit
=============  ========================================= =================
eeg            scalp electroencephalography (EEG)        Volts

meg (mag)      Magnetoencephalography (magnetometers)    Teslas

meg (grad)     Magnetoencephalography (gradiometers)     Teslas/meter

ecg            Electrocardiography (ECG)                 Volts

seeg           Stereotactic EEG channels                 Volts

dbs            Deep brain stimulation (DBS)              Volts

ecog           Electrocorticography (ECoG)               Volts

fnirs (hbo)    Functional near-infrared spectroscopy     Moles/liter
               (oxyhemoglobin)

fnirs (hbr)    Functional near-infrared spectroscopy     Moles/liter
               (deoxyhemoglobin)

emg            Electromyography (EMG)                    Volts

bio            Miscellaneous biological channels (e.g.,  Arbitrary units
               skin conductance)

stim           stimulus (a.k.a. trigger) channels        Arbitrary units

resp           response-trigger channel                  Arbitrary units

chpi           continuous head position indicator        Teslas
               (HPI) coil channels

exci           Flux excitation channel

ias            Internal Active Shielding data
               (Triux systems only?)

syst           System status channel information
               (Triux systems only)
=============  ========================================= =================
:orphan:

Internal representation (units)
===============================

.. NOTE: part of this file is included in doc/manual/io.rst and
   doc/overview/implementation.rst. Changes here are reflected there. If you
   want to link to this content, link to :ref:`manual-units` for the manual or
   :ref:`units` for the implementation page. The next line is a target for
   :start-after: so we can omit what's above:
   units-begin-content

Irrespective of the units used in your manufacturer's format, when importing
data, MNE-Python will always convert measurements to the same standard units.
Thus the in-memory representation of data are always in:

- Volts (eeg, eog, seeg, emg, ecg, bio, ecog, dbs)
- Teslas (magnetometers)
- Teslas/meter (gradiometers)
- Amperes*meter (dipole fits, minimum-norm estimates, etc.)
- Moles/liter ("molar"; fNIRS data: oxyhemoglobin (hbo), deoxyhemoglobin (hbr))
- Arbitrary units (various derived unitless quantities)

.. NOTE: this is a target for :end-before: units-end-of-list

Note, however, that most MNE-Python plotting functions will scale the data when
plotted to yield nice-looking axis annotations in a sensible range; for
example, :meth:`mne.io.Raw.plot_psd` will convert teslas to femtoteslas (fT)
and volts to microvolts (µV) when plotting MEG and EEG data.

The units used in internal data representation are particularly important to
remember when extracting data from MNE-Python objects and manipulating it
outside MNE-Python (e.g., when using methods like :meth:`~mne.io.Raw.get_data`
or :meth:`~mne.Epochs.to_data_frame` to convert data to :class:`NumPy arrays
<numpy.ndarray>` or :class:`Pandas DataFrames <pandas.DataFrame>` for analysis
or plotting with other Python modules).
:orphan:

Floating-point precision
========================

.. NOTE: part of this file is included in doc/manual/io.rst and
   doc/overview/implementation.rst. Changes here are reflected there. If you
   want to link to this content, link to :ref:`manual-precision` for the manual
   or :ref:`precision` for the implementation page. The next line is a target
   for :start-after: so we can omit the title above:
   precision-begin-content

MNE-Python performs all computation in memory using the double-precision 64-bit
floating point format. This means that the data is typecast into float64 format
as soon as it is read into memory. The reason for this is that operations such
as filtering and preprocessing are more accurate when using the 64-bit format.
However, for backward compatibility, MNE-Python writes :file:`.fif` files in a
32-bit format by default. This reduces file size when saving data to disk, but
beware that *saving intermediate results to disk and re-loading them from disk
later may lead to loss in precision*. If you would like to ensure 64-bit
precision, there are two possibilities:

- Chain the operations in memory and avoid saving intermediate results.

- Save intermediate results but change the :class:`~numpy.dtype` used for
  saving, by using the ``fmt`` parameter of :meth:`mne.io.Raw.save` (or
  :meth:`mne.Epochs.save`, etc). However, note that this may render the
  :file:`.fif` files unreadable in software packages other than MNE-Python.
:orphan:

Bad channel repair via interpolation
====================================

Spherical spline interpolation (EEG)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. NOTE: part of this file is included in doc/overview/implementation.rst.
   Changes here are reflected there. If you want to link to this content, link
   to :ref:`channel-interpolation` to link to that section of the
   implementation.rst page. The next line is a target for :start-after: so we
   can omit the title from the include:
   channel-interpolation-begin-content

In short, data repair using spherical spline interpolation :footcite:`PerrinEtAl1989` consists of the following steps:

* Project the good and bad electrodes onto a unit sphere
* Compute a mapping matrix that maps :math:`N` good channels to :math:`M` bad channels
* Use this mapping matrix to compute interpolated data in the bad channels

Spherical splines assume that the potential :math:`V(\boldsymbol{r_i})` at any point :math:`\boldsymbol{r_i}` on the surface of the sphere can be represented by:

.. math:: V(\boldsymbol{r_i}) = c_0 + \sum_{j=1}^{N}c_{i}g_{m}(cos(\boldsymbol{r_i}, \boldsymbol{r_{j}}))
   :label: model

where the :math:`C = (c_{1}, ..., c_{N})^{T}` are constants which must be estimated. The function :math:`g_{m}(\cdot)` of order :math:`m` is given by:

.. math:: g_{m}(x) = \frac{1}{4 \pi}\sum_{n=1}^{\infty} \frac{2n + 1}{(n(n + 1))^m}P_{n}(x)
   :label: legendre

where :math:`P_{n}(x)` are `Legendre polynomials`_ of order :math:`n`.

.. _Legendre polynomials: https://en.wikipedia.org/wiki/Legendre_polynomials

To estimate the constants :math:`C`, we must solve the following two equations simultaneously:

.. math:: G_{ss}C + T_{s}c_0 = X
   :label: matrix_form

.. math:: {T_s}^{T}C = 0
   :label: constraint

where :math:`G_{ss} \in R^{N \times N}` is a matrix whose entries are :math:`G_{ss}[i, j] = g_{m}(cos(\boldsymbol{r_i}, \boldsymbol{r_j}))` and :math:`X \in R^{N \times 1}` are the potentials :math:`V(\boldsymbol{r_i})` measured at the good channels. :math:`T_{s} = (1, 1, ..., 1)^\top` is a column vector of dimension :math:`N`. Equation :eq:`matrix_form` is the matrix formulation of Equation :eq:`model` and equation :eq:`constraint` is like applying an average reference to the data. From equation :eq:`matrix_form` and :eq:`constraint`, we get:

.. math:: \begin{bmatrix} c_0 \\ C \end{bmatrix} = {\begin{bmatrix} {T_s}^{T} && 0 \\ T_s && G_{ss} \end{bmatrix}}^{-1} \begin{bmatrix} 0 \\ X \end{bmatrix} = C_{i}X
   :label: estimate_constant

:math:`C_{i}` is the same as matrix :math:`{\begin{bmatrix} {T_s}^{T} && 0 \\ T_s && G_{ss} \end{bmatrix}}^{-1}` but with its first column deleted, therefore giving a matrix of dimension :math:`(N + 1) \times N`.

Now, to estimate the potentials :math:`\hat{X} \in R^{M \times 1}` at the bad channels, we have to do:

.. math:: \hat{X} = G_{ds}C + T_{d}c_0
   :label: estimate_data

where :math:`G_{ds} \in R^{M \times N}` computes :math:`g_{m}(\boldsymbol{r_i}, \boldsymbol{r_j})` between the bad and good channels. :math:`T_{d} = (1, 1, ..., 1)^\top` is a column vector of dimension :math:`M`. Plugging in equation :eq:`estimate_constant` in :eq:`estimate_data`, we get

.. math:: \hat{X} = \begin{bmatrix} T_d && G_{ds} \end{bmatrix} \begin{bmatrix} c_0 \\ C \end{bmatrix} = \underbrace{\begin{bmatrix} T_d && G_{ds} \end{bmatrix} C_{i}}_\text{mapping matrix}X


To interpolate bad channels, one can simply do:

	>>> evoked.interpolate_bads(reset_bads=False)  # doctest: +SKIP

and the bad channel will be fixed.

.. target for :end-before: channel-interpolation-end-content

.. topic:: Examples:

	* :ref:`ex-interpolate-bad-channels`
:orphan:
.. _dig-formats:

Supported formats for digitized 3D locations
============================================

.. NOTE: If you want to link to this content, link to :ref:`dig-formats`
   for the implementation page. The next line is
   a target for :start-after: so we can omit the title above:
   dig-formats-begin-content

MNE-Python can load 3D point locations obtained by digitization systems.
Such files allow to obtain a :class:`montage <mne.channels.DigMontage>`
that can then be added to :class:`~mne.io.Raw` objects with the
:meth:`~mne.io.Raw.set_montage`. See the documentation for each reader
function for more info on reading specific file types.

.. NOTE: To include only the table, here's a different target for :start-after:
   dig-formats-begin-table

.. cssclass:: table-bordered
.. rst-class:: midvalign

=================  ================  ==============================================
Vendor             Extension(s)      MNE-Python function
=================  ================  ==============================================
Neuromag           .fif              :func:`mne.channels.read_dig_fif`

Polhemus ISOTRAK   .hsp, .elp, .eeg  :func:`mne.channels.read_dig_polhemus_isotrak`

EGI                .xml              :func:`mne.channels.read_dig_egi`

MNE-C              .hpts             :func:`mne.channels.read_dig_hpts`

Brain Products     .bvct             :func:`mne.channels.read_dig_captrak`

Compumedics        .dat              :func:`mne.channels.read_dig_dat`
=================  ================  ==============================================

To load Polhemus FastSCAN files you can use
:func:`montage <mne.channels.read_polhemus_fastscan>`.

It is also possible to make a :class:`montage <mne.channels.DigMontage>`
from arrays with :func:`mne.channels.make_dig_montage`.
:orphan:

Memory-efficient I/O
====================

.. NOTE: part of this file is included in doc/overview/implementation.rst.
   Changes here are reflected there. If you want to link to this content, link
   to :ref:`memory` to link to that section of the implementation.rst
   page. The next line is a target for :start-after: so we can omit the title
   from the include:
   memory-begin-content


Preloading continuous (raw) data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MNE-Python can read data on-demand using the ``preload`` option provided in
raw reading functions. For example::

    from mne import io
    from mne.datasets import sample
    data_path = sample.data_path()
    raw_fname = data_path / 'MEG' / 'sample' / 'sample_audvis_filt-0-40_raw.fif'
    raw = io.read_raw_fif(raw_fname, preload=False)

.. note:: Filtering, resampling and dropping or selecting channels does not
          work with ``preload=False``.


Preloading epoched data
~~~~~~~~~~~~~~~~~~~~~~~

Similarly, epochs can also be be read from disk on-demand. For example::

    import mne
    events = mne.find_events(raw)
    event_id, tmin, tmax = 1, -0.2, 0.5
    picks = mne.pick_types(raw.info, meg=True, eeg=True, stim=False, eog=True)
    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
                        baseline=(None, 0), reject=dict(eeg=80e-6, eog=150e-6),
                        preload=False)

When ``preload=False``, the epochs data is loaded from the disk on-demand. Note
that ``preload=False`` for epochs will work even if the ``raw`` object has been
loaded with ``preload=True``. Preloading is also supported for
:func:`mne.read_epochs`.

.. warning:: This comes with a caveat. When ``preload=False``, data rejection
             based on peak-to-peak thresholds is executed when the data is
             loaded from disk, *not* when the ``Epochs`` object is created.

To explicitly reject artifacts with ``preload=False``, use the function :func:`mne.Epochs.drop_bad`.


Loading data explicitly
~~~~~~~~~~~~~~~~~~~~~~~

To load the data if ``preload=False`` was initially selected, use the functions :func:`mne.io.Raw.load_data` and :func:`mne.Epochs.load_data`.


Accessing data as NumPy arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you just want your raw data as a :class:`Numpy array <numpy.ndarray>` to
work with it in a different framework you can use slicing syntax::

    first_channel_data, times = raw[0, :]
    channels_3_and_4, times = raw[3:5, :]
.. _ch_mne:

The minimum-norm current estimates
==================================

.. NOTE: part of this file is included in doc/overview/implementation.rst.
   Changes here are reflected there. If you want to link to this content, link
   to :ref:`ch_mne` to link to that section of the implementation.rst page.
   The next line is a target for :start-after: so we can omit the title from
   the include:
   inverse-begin-content


This section describes the mathematical details of the calculation of
minimum-norm estimates. In Bayesian sense, the ensuing current distribution is
the maximum a posteriori (MAP) estimate under the following assumptions:

- The viable locations of the currents are constrained to the cortex.
  Optionally, the current orientations can be fixed to be normal to the
  cortical mantle.

- The amplitudes of the currents have a Gaussian prior distribution with a
  known source covariance matrix.

- The measured data contain additive noise with a Gaussian distribution with a
  known covariance matrix. The noise is not correlated over time.

Computing the inverse operator is accomplished using
:func:`mne.minimum_norm.make_inverse_operator` and
:func:`mne.minimum_norm.apply_inverse`. The use of these functions is presented
in the tutorial :ref:`tut-inverse-methods`.

The linear inverse operator
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The measured data in the source estimation procedure consists of MEG and EEG
data, recorded on a total of N channels. The task is to estimate a total of
:math:`Q`
strengths of sources located on the cortical mantle. If the number of source
locations is :math:`P`, :math:`Q = P` for fixed-orientation sources and
:math:`Q = 3P` if the source
orientations are unconstrained. The regularized linear inverse operator
following from regularized maximal likelihood of the above probabilistic model
is given by the :math:`Q \times N` matrix

.. math::    M = R' G^\top (G R' G^\top + C)^{-1}\ ,

where :math:`G` is the gain matrix relating the source strengths to the measured
MEG/EEG data, :math:`C` is the data noise-covariance matrix and :math:`R'` is
the source covariance matrix. The dimensions of these matrices are :math:`N
\times Q`, :math:`N \times N`, and :math:`Q \times Q`, respectively. The
:math:`Q \times 1` source-strength vector is obtained by multiplying the
:math:`Q \times 1` data vector by :math:`Q`.

The expected value of the current amplitudes at time *t* is then given by
:math:`\hat{j}(t) = Mx(t)`, where :math:`x(t)` is a vector containing the
measured MEG and EEG data values at time *t*.

For computational convenience, the linear inverse operator is
not computed explicitly. See :ref:`mne_solution` for mathematical
details, and :ref:`CIHCFJEI` for a detailed example.

.. _mne_regularization:

Regularization
~~~~~~~~~~~~~~

The a priori variance of the currents is, in practice, unknown. We can express
this by writing :math:`R' = R/ \lambda^2 = R \lambda^{-2}`, which yields the
inverse operator

.. math::
   :label: inv_m

    M &= R' G^\top (G R' G^\top + C)^{-1} \\
      &= R \lambda^{-2} G^\top (G R \lambda^{-2} G^\top + C)^{-1} \\
      &= R \lambda^{-2} G^\top \lambda^2 (G R G^\top + \lambda^2 C)^{-1} \\
      &= R G^\top (G R G^\top + \lambda^2 C)^{-1}\ ,

where the unknown current amplitude is now interpreted in terms of the
regularization parameter :math:`\lambda^2`. Larger :math:`\lambda^2` values
correspond to spatially smoother and weaker current amplitudes, whereas smaller
:math:`\lambda^2` values lead to the opposite.

We can arrive at the regularized linear inverse operator also by minimizing a
cost function :math:`S` with respect to the estimated current :math:`\hat{j}`
(given the measurement vector :math:`x` at any given time :math:`t`) as

.. math::

    \min_\hat{j} \Bigl\{ S \Bigr\} &= \min_\hat{j} \Bigl\{ \tilde{e}^\top \tilde{e} + \lambda^2 \hat{j}^\top R^{-1} \hat{j} \Bigr\} \\
                                   &= \min_\hat{j} \Bigl\{ (x - G\hat{j})^\top C^{-1} (x - G\hat{j}) + \lambda^2 \hat{j}^\top R^{-1} \hat{j} \Bigr\} \,

where the first term consists of the difference between the whitened measured
data (see :ref:`whitening_and_scaling`) and those predicted by the model while the
second term is a weighted-norm of the current estimate. It is seen that, with
increasing :math:`\lambda^2`, the source term receive more weight and larger
discrepancy between the measured and predicted data is tolerable.

.. _whitening_and_scaling:

Whitening and scaling
~~~~~~~~~~~~~~~~~~~~~

The MNE software employs data whitening so that a 'whitened' inverse operator
assumes the form

.. math::    \tilde{M} = M C^{^1/_2} = R \tilde{G}^\top (\tilde{G} R \tilde{G}^\top + \lambda^2 I)^{-1}\ ,
   :label: inv_m_tilde

where

.. math:: \tilde{G} = C^{-^1/_2}G
   :label: inv_g_tilde

is the spatially whitened gain matrix. We arrive at the whitened inverse
operator equation :eq:`inv_m_tilde` by making the substitution for
:math:`G` from :eq:`inv_g_tilde` in :eq:`inv_m` as

.. math::

    \tilde{M} = M C^{^1/_2} &= R G^\top (G R G^\top + \lambda^2 C)^{-1} C^{^1/_2} \\
                             &= R \tilde{G}^\top C^{^1/_2} (C^{^1/_2} \tilde{G} R \tilde{G}^\top C^{^1/_2} + \lambda^2 C)^{-1} C^{^1/_2} \\
                             &= R \tilde{G}^\top C^{^1/_2} (C^{^1/_2} (\tilde{G} R \tilde{G}^\top + \lambda^2 I) C^{^1/_2})^{-1} C^{^1/_2} \\
                             &= R \tilde{G}^\top C^{^1/_2} C^{-^1/_2} (\tilde{G} R \tilde{G}^\top + \lambda^2 I)^{-1} C^{-^1/_2} C^{^1/_2} \\
                             &= R \tilde{G}^\top (\tilde{G} R \tilde{G}^\top + \lambda^2 I)^{-1}\ .

The expected current values are

.. math::
   :label: inv_j_hat_t

    \hat{j}(t) &= Mx(t) \\
               &= M C^{^1/_2} C^{-^1/_2} x(t) \\
               &= \tilde{M} \tilde{x}(t)

knowing :eq:`inv_m_tilde` and taking

.. math::
   :label: inv_tilde_x_t

    \tilde{x}(t) = C^{-^1/_2}x(t)

as the whitened measurement vector at time *t*. The spatial
whitening operator :math:`C^{-^1/_2}` is obtained with the help of the
eigenvalue decomposition
:math:`C = U_C \Lambda_C^2 U_C^\top` as :math:`C^{-^1/_2} = \Lambda_C^{-1} U_C^\top`.
In the MNE software the noise-covariance matrix is stored as the one applying
to raw data. To reflect the decrease of noise due to averaging, this matrix,
:math:`C_0`, is scaled by the number of averages, :math:`L`, *i.e.*, :math:`C =
C_0 / L`.

As shown above, regularization of the inverse solution is equivalent to a
change in the variance of the current amplitudes in the Bayesian *a priori*
distribution.

A convenient choice for the source-covariance matrix :math:`R` is such that
:math:`\text{trace}(\tilde{G} R \tilde{G}^\top) / \text{trace}(I) = 1`. With this
choice we can approximate :math:`\lambda^2 \sim 1/\rm{SNR}^2`, where SNR is the
(amplitude) signal-to-noise ratio of the whitened data.

.. note::
   The definition of the signal to noise-ratio/ :math:`\lambda^2` relationship
   given above works nicely for the whitened forward solution. In the
   un-whitened case scaling with the trace ratio :math:`\text{trace}(GRG^\top) /
   \text{trace}(C)` does not make sense, since the diagonal elements summed
   have, in general, different units of measure. For example, the MEG data are
   expressed in T or T/m whereas the unit of EEG is Volts.

See :ref:`tut-compute-covariance` for example of noise covariance computation
and whitening.

.. _cov_regularization_math:

Regularization of the noise-covariance matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since finite amount of data is usually available to compute an estimate of the
noise-covariance matrix :math:`C`, the smallest eigenvalues of its estimate are
usually inaccurate and smaller than the true eigenvalues. Depending on the
seriousness of this problem, the following quantities can be affected:

- The model data predicted by the current estimate,

- Estimates of signal-to-noise ratios, which lead to estimates of the required
  regularization, see :ref:`mne_regularization`,

- The estimated current values, and

- The noise-normalized estimates, see :ref:`noise_normalization`.

Fortunately, the latter two are least likely to be affected due to
regularization of the estimates. However, in some cases especially the EEG part
of the noise-covariance matrix estimate can be deficient, *i.e.*, it may
possess very small eigenvalues and thus regularization of the noise-covariance
matrix is advisable.

Historically, the MNE software accomplishes the regularization by replacing a
noise-covariance matrix estimate :math:`C` with

.. math::    C' = C + \sum_k {\varepsilon_k \bar{\sigma_k}^2 I^{(k)}}\ ,

where the index :math:`k` goes across the different channel groups (MEG planar
gradiometers, MEG axial gradiometers and magnetometers, and EEG),
:math:`\varepsilon_k` are the corresponding regularization factors,
:math:`\bar{\sigma_k}` are the average variances across the channel groups, and
:math:`I^{(k)}` are diagonal matrices containing ones at the positions
corresponding to the channels contained in each channel group.

See :ref:`plot_compute_covariance_howto` for details on computing and
regularizing the channel covariance matrix.

.. _mne_solution:

Computation of the solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most straightforward approach to calculate the MNE is to employ the
expression of the original or whitened inverse operator directly. However, for
computational convenience we prefer to take another route, which employs the
singular-value decomposition (SVD) of the matrix

.. math::
   :label: inv_a

    A &= \tilde{G} R^{^1/_2} \\
      &= U \Lambda V^\top

where the superscript :math:`^1/_2` indicates a square root of :math:`R`. For a
diagonal matrix, one simply takes the square root of :math:`R` while in the
more general case one can use the Cholesky factorization :math:`R = R_C R_C^\top`
and thus :math:`R^{^1/_2} = R_C`.

Combining the SVD from :eq:`inv_a` with the inverse equation :eq:`inv_m` it is
easy to show that

.. math::
   :label: inv_m_tilde_svd

    \tilde{M} &= R \tilde{G}^\top (\tilde{G} R \tilde{G}^\top + \lambda^2 I)^{-1} \\
              &= R^{^1/_2} A^\top (A A^\top + \lambda^2 I)^{-1} \\
              &= R^{^1/_2} V \Lambda U^\top (U \Lambda V^\top V \Lambda U^\top + \lambda^2 I)^{-1} \\
              &= R^{^1/_2} V \Lambda U^\top (U (\Lambda^2 + \lambda^2 I) U^\top)^{-1} \\
              &= R^{^1/_2} V \Lambda U^\top U (\Lambda^2 + \lambda^2 I)^{-1} U^\top \\
              &= R^{^1/_2} V \Lambda (\Lambda^2 + \lambda^2 I)^{-1} U^\top \\
              &= R^{^1/_2} V \Gamma U^\top

where the elements of the diagonal matrix :math:`\Gamma` are simply

.. `reginv` in our code:

.. math::
   :label: inv_gamma_k

    \gamma_k = \frac{\lambda_k}{\lambda_k^2 + \lambda^2}\ .

From our expected current equation :eq:`inv_j_hat_t` and our whitened
measurement equation :eq:`inv_tilde_x_t`, if we take

.. math::
   :label: inv_w_t

    w(t) &= U^\top \tilde{x}(t) \\
         &= U^\top C^{-^1/_2} x(t)\ ,

we can see that the expression for the expected current is just

.. math::
   :label: inv_j_hat_t_svd

    \hat{j}(t) &= R^{^1/_2} V \Gamma w(t) \\
               &= \sum_k {\bar{v_k} \gamma_k w_k(t)}\ ,

where :math:`\bar{v_k} = R^{^1/_2} v_k`, with :math:`v_k` being the
:math:`k` th column of :math:`V`. It is thus seen that the current estimate is
a weighted sum of the "weighted" eigenleads :math:`v_k`.

It is easy to see that :math:`w(t) \propto \sqrt{L}`. To maintain the relation
:math:`(\tilde{G} R \tilde{G}^\top) / \text{trace}(I) = 1` when :math:`L` changes
we must have :math:`R \propto 1/L`. With this approach, :math:`\lambda_k` is
independent of  :math:`L` and, for fixed :math:`\lambda`, we see directly that
:math:`j(t)` is independent of :math:`L`.

The minimum-norm estimate is computed using this procedure in
:func:`mne.minimum_norm.make_inverse_operator`, and its usage is illustrated
in :ref:`CIHCFJEI`.


.. _noise_normalization:

Noise normalization
~~~~~~~~~~~~~~~~~~~

Noise normalization serves three purposes:

- It converts the expected current value into a dimensionless statistical test
  variable. Thus the resulting time and location dependent values are often
  referred to as dynamic statistical parameter maps (dSPM).

- It reduces the location bias of the estimates. In particular, the tendency of
  the MNE to prefer superficial currents is eliminated.

- The width of the point-spread function becomes less dependent on the source
  location on the cortical mantle. The point-spread is defined as the MNE
  resulting from the signals coming from a point current source (a current
  dipole) located at a certain point on the cortex.

In practice, noise normalization is implemented as a division by the square
root of the estimated variance of each voxel. In computing these noise
normalization factors, it's convenient to reuse our "weighted eigenleads"
definition from equation :eq:`inv_j_hat_t` in matrix form as

.. math::
   :label: inv_eigenleads_weighted

    \bar{V} = R^{^1/_2} V\ .

dSPM
----

Noise-normalized linear estimates introduced by Dale et al.
:footcite:`DaleEtAl1999` require division of the expected current amplitude by
its variance. In practice, this requires the computation of the diagonal
elements of the following matrix, using SVD equation :eq:`inv_m_tilde` and
:eq:`inv_eigenleads_weighted`:

.. math::

    M C M^\top &= M C^{^1/_2} C^{^1/_2} M^\top \\
            &= \tilde{M} \tilde{M}^\top \\
            &= R^{^1/_2} V \Gamma U^\top U \Gamma V^\top R^{^1/_2} \\
            &= \bar{V} \Gamma^2 \bar{V}^\top\ .

Because we only care about the diagonal entries here, we can find the
variances for each source as

.. math::

    \sigma_k^2 = \gamma_k^2

Under the conditions expressed at the end of :ref:`mne_solution`, it
follows that the *t*-statistic values associated with fixed-orientation
sources) are thus proportional to :math:`\sqrt{L}` while the *F*-statistic
employed with free-orientation sources is proportional to :math:`L`,
correspondingly.

.. note::
   The MNE software usually computes the *square roots* of the F-statistic to
   be displayed on the inflated cortical surfaces. These are also proportional
   to :math:`\sqrt{L}`.

sLORETA
-------
sLORETA :footcite:`Pascual-Marqui2002` estimates the current variances as the
diagonal entries of the
resolution matrix, which is the product of the inverse and forward operators.
In other words, the diagonal entries of (using :eq:`inv_m_tilde_svd`,
:eq:`inv_g_tilde`, and :eq:`inv_a`)

.. math::

    M G &= M C^{^1/_2} C^{-^1/_2} G \\
        &= \tilde{M} \tilde{G} \\
        &= R^{^1/_2} V \Gamma U^\top \tilde{G} R^{^1/_2} R^{-^1/_2} \\
        &= R^{^1/_2} V \Gamma U^\top U \Lambda V^\top R^{-^1/_2} \\
        &= R^{^1/_2} V \Gamma U^\top U \Lambda V^\top R^{^1/_2} R^{-1} \\
        &= \bar{V} \Gamma U^\top U \Lambda \bar{V}^\top R^{-1} \\
        &= \bar{V} \Gamma \Lambda \bar{V}^\top R^{-1}\ .

Because :math:`R` is diagonal and we only care about the diagonal entries,
we can find our variance estimates as

.. math::

    \sigma_k^2 &= \gamma_k \lambda_k R_{k,k}^{-1} \\
               &= \left(\frac{\lambda_k}{(\lambda_k^2 + \lambda^2)}\right) \left(\frac{\lambda_k}{1}\right) \left(\frac{1}{\lambda^2}\right) \\
               &= \frac{\lambda_k^2}{(\lambda_k^2 + \lambda^2) \lambda^2} \\
               &= \left(\frac{\lambda_k^2}{(\lambda_k^2 + \lambda^2)^2}\right) \left(\frac{\lambda^2 + \lambda_k^2}{\lambda^2}\right) \\
               &= \left(\frac{\lambda_k}{\lambda_k^2 + \lambda^2}\right)^2 \left(1 + \frac{\lambda_k^2}{\lambda^2}\right) \\
               &= \gamma_k^2 \left(1 + \frac{\lambda_k^2}{\lambda^2}\right)\ .

eLORETA
~~~~~~~
While dSPM and sLORETA solve for noise normalization weights
:math:`\sigma^2_k` that are applied to standard minimum-norm estimates
:math:`\hat{j}(t)`, eLORETA :footcite:`Pascual-Marqui2011` instead solves for
a source covariance
matrix :math:`R` that achieves zero localization bias. For fixed-orientation
solutions the resulting matrix :math:`R` will be a diagonal matrix, and for
free-orientation solutions it will be a block-diagonal matrix with
:math:`3 \times 3` blocks.

.. In https://royalsocietypublishing.org/doi/full/10.1098/rsta.2011.0081
.. eq. 2.10 (classical min norm), their values map onto our values as:
..
.. - α=λ²
.. - W=R⁻¹ (pos semidef weight matrix)
.. - K=G
.. - ϕ=x
.. - C=H
..

In :footcite:`Pascual-Marqui2011` eq. 2.13 states that the following system
of equations can be used to find the weights, :math:`\forall i \in {1, ..., P}`
(note that here we represent the equations from that paper using our notation):

.. math:: r_i = \left[ G_i^\top \left( GRG^\top + \lambda^2C \right)^{-1} G_i \right] ^{-^1/_2}

And an iterative algorithm can be used to find the values for the weights
:math:`r_i` that satisfy these equations as:

1. Initialize identity weights.
2. Compute :math:`N= \left( GRG^\top + \lambda^2C \right)^{-1}`.
3. Holding :math:`N` fixed, compute new weights :math:`r_i = \left[ G_i^\top N G_i \right]^{-^1/_2}`.
4. Using new weights, go to step (2) until convergence.

In particular, for step (2) we can use our substitution from :eq:`inv_g_tilde`
as:

.. math::

    N &= (G R G^\top + \lambda^2 C)^{-1} \\
      &= (C^{^1/_2} \tilde{G} R \tilde{G}^\top C^{^1/_2} + \lambda^2 C)^{-1} \\
      &= (C^{^1/_2} (\tilde{G} R \tilde{G}^\top + \lambda^2 I) C^{^1/_2})^{-1} \\
      &= C^{-^1/_2} (\tilde{G} R \tilde{G}^\top + \lambda^2 I)^{-1} C^{-^1/_2} \\
      &= C^{-^1/_2} (\tilde{G} R \tilde{G}^\top + \lambda^2 I)^{-1} C^{-^1/_2}\ .

Then defining :math:`\tilde{N}` as the whitened version of :math:`N`, i.e.,
the regularized pseudoinverse of :math:`\tilde{G}R\tilde{G}^\top`, we can
compute :math:`N` as:

.. math::

    N &= C^{-^1/_2} (U_{\tilde{G}R\tilde{G}^\top} \Lambda_{\tilde{G}R\tilde{G}^\top} V_{\tilde{G}R\tilde{G}^\top}^\top + \lambda^2 I)^{-1} C^{-^1/_2} \\
      &= C^{-^1/_2} (U_{\tilde{G}R\tilde{G}^\top} (\Lambda_{\tilde{G}R\tilde{G}^\top} + \lambda^2 I) V_{\tilde{G}R\tilde{G}^\top}^\top)^{-1} C^{-^1/_2} \\
      &= C^{-^1/_2} V_{\tilde{G}R\tilde{G}^\top} (\Lambda_{\tilde{G}R\tilde{G}^\top} + \lambda^2 I)^{-1} U_{\tilde{G}R\tilde{G}^\top}^\top C^{-^1/_2} \\
      &= C^{-^1/_2} \tilde{N} C^{-^1/_2}\ .

In step (3) we left and right multiply with subsets of :math:`G`, but making
the substitution :eq:`inv_g_tilde` we see that we equivalently compute:

.. math::

    r_i &= \left[ G_i^\top N G_i \right]^{-^1/_2} \\
        &= \left[ (C^{^1/_2} \tilde{G}_i)^\top N C^{^1/_2} \tilde{G}_i \right]^{-^1/_2} \\
        &= \left[ \tilde{G}_i^\top C^{^1/_2} N C^{^1/_2} \tilde{G}_i \right]^{-^1/_2} \\
        &= \left[ \tilde{G}_i^\top C^{^1/_2} C^{-^1/_2} \tilde{N} C^{-^1/_2} C^{^1/_2} \tilde{G}_i \right]^{-^1/_2} \\
        &= \left[ \tilde{G}_i^\top \tilde{N} \tilde{G}_i \right]^{-^1/_2}\ .

For convenience, we thus never need to compute :math:`N` itself but can instead
compute the whitened version :math:`\tilde{N}`.

Predicted data
~~~~~~~~~~~~~~

Under noiseless conditions the SNR is infinite and thus leads to
:math:`\lambda^2 = 0` and the minimum-norm estimate explains the measured data
perfectly. Under realistic conditions, however, :math:`\lambda^2 > 0` and there
is a misfit between measured data and those predicted by the MNE. Comparison of
the predicted data, here denoted by :math:`x(t)`, and measured one can give
valuable insight on the correctness of the regularization applied.

In the SVD approach we easily find

.. math::    \hat{x}(t) = G \hat{j}(t) = C^{^1/_2} U \Pi w(t)\ ,

where the diagonal matrix :math:`\Pi` has elements :math:`\pi_k = \lambda_k
\gamma_k` The predicted data is thus expressed as the weighted sum of the
'recolored eigenfields' in :math:`C^{^1/_2} U`.

Cortical patch statistics
~~~~~~~~~~~~~~~~~~~~~~~~~

If the ``add_dists=True`` option was used in source space creation,
the source space file will contain
Cortical Patch Statistics (CPS) for each vertex of the cortical surface. The
CPS provide information about the source space point closest to it as well as
the distance from the vertex to this source space point. The vertices for which
a given source space point is the nearest one define the cortical patch
associated with with the source space point. Once these data are available, it
is straightforward to compute the following cortical patch statistics for each
source location :math:`d`:

- The average over the normals of at the vertices in a patch,
  :math:`\bar{n_d}`,

- The areas of the patches, :math:`A_d`, and

- The average deviation of the vertex normals in a patch from their average,
  :math:`\sigma_d`, given in degrees.

``use_cps`` parameter in :func:`mne.convert_forward_solution`, and
:func:`mne.minimum_norm.make_inverse_operator` controls whether to use
cortical patch statistics (CPS) to define normal orientations or not (see
:ref:`CHDBBCEJ`).

.. _inverse_orientation_constraints:

Orientation constraints
~~~~~~~~~~~~~~~~~~~~~~~

The principal sources of MEG and EEG signals are generally believed to be
postsynaptic currents in the cortical pyramidal neurons. Since the net primary
current associated with these microscopic events is oriented normal to the
cortical mantle, it is reasonable to use the cortical normal orientation as a
constraint in source estimation. In addition to allowing completely free source
orientations, the MNE software implements three orientation constraints based
of the surface normal data:

- Source orientation can be rigidly fixed to the surface normal direction by
  specifying ``fixed=True`` in :func:`mne.minimum_norm.make_inverse_operator`.
  If cortical patch statistics are available the average
  normal over each patch, :math:`\bar{n_d}`, are used to define the source
  orientation. Otherwise, the vertex normal at the source space location is
  employed.

- A *location independent or fixed loose orientation constraint* (fLOC) can be
  employed by specifying ``fixed=False`` and ``loose=1.0`` when
  calling :func:`mne.minimum_norm.make_inverse_operator` (see
  :ref:`plot_dipole_orientations_fLOC_orientations`).
  In this approach, a source coordinate
  system based on the local surface orientation at the source location is
  employed. By default, the three columns of the gain matrix G, associated with
  a given source location, are the fields of unit dipoles pointing to the
  directions of the :math:`x`, :math:`y`, and :math:`z` axis of the coordinate
  system employed in the forward calculation (usually the :ref:`MEG head
  coordinate frame <head_device_coords>`). For LOC the orientation is changed so
  that the first two source components lie in the plane normal to the surface
  normal at the source location and the third component is aligned with it.
  Thereafter, the variance of the source components tangential to the cortical
  surface are reduced by a factor defined by the ``--loose`` option.

- A *variable loose orientation constraint* (vLOC) can be employed by
  specifying ``fixed=False`` and ``loose`` parameters when calling
  :func:`mne.minimum_norm.make_inverse_operator` (see
  :ref:`plot_dipole_orientations_vLOC_orientations`). This
  is similar to *fLOC* except that the value given with the ``loose``
  parameter will be multiplied by :math:`\sigma_d`, defined above.

Depth weighting
~~~~~~~~~~~~~~~

The minimum-norm estimates have a bias towards superficial currents. This
tendency can be alleviated by adjusting the source covariance matrix :math:`R`
to favor deeper source locations. In the depth weighting scheme employed in MNE
analyze, the elements of :math:`R` corresponding to the :math:`p` th source
location are be scaled by a factor

.. math::    f_p = (g_{1p}^\top g_{1p} + g_{2p}^\top g_{2p} + g_{3p}^\top g_{3p})^{-\gamma}\ ,

where :math:`g_{1p}`, :math:`g_{2p}`, and :math:`g_{3p}` are the three columns
of :math:`G` corresponding to source location :math:`p` and :math:`\gamma` is
the order of the depth weighting, which is specified via the ``depth`` option
in :func:`mne.minimum_norm.make_inverse_operator`.

Effective number of averages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is often the case that the epoch to be analyzed is a linear combination over
conditions rather than one of the original averages computed. As stated above,
the noise-covariance matrix computed is originally one corresponding to raw
data. Therefore, it has to be scaled correctly to correspond to the actual or
effective number of epochs in the condition to be analyzed. In general, we have

.. math::    C = C_0 / L_{eff}

where :math:`L_{eff}` is the effective number of averages. To calculate
:math:`L_{eff}` for an arbitrary linear combination of conditions

.. math::    y(t) = \sum_{i = 1}^n {w_i x_i(t)}

we make use of the the fact that the noise-covariance matrix

.. math::    C_y = \sum_{i = 1}^n {w_i^2 C_{x_i}} = C_0 \sum_{i = 1}^n {w_i^2 / L_i}

which leads to

.. math::    1 / L_{eff} = \sum_{i = 1}^n {w_i^2 / L_i}

An important special case  of the above is a weighted average, where

.. math::    w_i = L_i / \sum_{i = 1}^n {L_i}

and, therefore

.. math::    L_{eff} = \sum_{i = 1}^n {L_i}

Instead of a weighted average, one often computes a weighted sum, a simplest
case being a difference or sum of two categories. For a difference :math:`w_1 =
1` and :math:`w_2 = -1` and thus

.. math::    1 / L_{eff} = 1 / L_1 + 1 / L_2

or

.. math::    L_{eff} = \frac{L_1 L_2}{L_1 + L_2}

Interestingly, the same holds for a sum, where :math:`w_1 = w_2 = 1`.
Generalizing, for any combination of sums and differences, where :math:`w_i =
1` or :math:`w_i = -1`, :math:`i = 1 \dotso n`, we have

.. math::    1 / L_{eff} = \sum_{i = 1}^n {1/{L_i}}

.. target for :end-before: inverse-end-content
:orphan:

Morphing and averaging source estimates
=======================================

The spherical morphing of BEM surfaces accomplished by FreeSurfer can be
employed to bring data from different subjects into a common anatomical frame.
This page describes utilities which make use of the spherical :term:`morphing`
procedure. :func:`mne.morph_labels` morphs label files between subjects
allowing the definition of labels in a one brain and transforming them to
anatomically analogous labels in another. :meth:`mne.SourceMorph.apply` offers
the capability to transform all subject data to the same space and,
e.g., compute averages of data across subjects.

.. NOTE: part of this file is included in doc/overview/implementation.rst.
   Changes here are reflected there. If you want to link to this content, link
   to :ref:`ch_morph` to link to that section of the implementation.rst page.
   The next line is a target for :start-after: so we can omit the title from
   the include:
   morph-begin-content


Why morphing?
~~~~~~~~~~~~~

.. sidebar:: Morphing examples in MNE-Python

   Examples of morphing in MNE-Python include :ref:`this tutorial
   <tut-mne-fixed-free>` on surface source estimation or these examples on
   :ref:`surface <ex-morph-surface>` and :ref:`volumetric <ex-morph-volume>`
   source estimation.

Modern neuroimaging techniques, such as source reconstruction or fMRI analyses,
make use of advanced mathematical models and hardware to map brain activity
patterns into a subject-specific anatomical brain space. This enables the study
of spatio-temporal brain activity. The representation of spatio-temporal brain
data is often mapped onto the anatomical brain structure to relate functional
and anatomical maps. Thereby activity patterns are overlaid with anatomical
locations that supposedly produced the activity. Anatomical MR images are often
used as such or are transformed into an inflated surface representations to
serve as  "canvas" for the visualization.

In order to compute group-level statistics, data representations across
subjects must be morphed to a common frame, such that anatomically and
functional similar structures are represented at the same spatial location for
*all subjects equally*. Since brains vary, :term:`morphing` comes into play to
tell us how the data produced by subject A would be represented on the brain of
subject B (and vice-versa).


The morphing maps
~~~~~~~~~~~~~~~~~

The MNE software accomplishes morphing with help of morphing maps.
The morphing is performed with help of the registered
spherical surfaces (``lh.sphere.reg`` and ``rh.sphere.reg`` ) which must be
produced in FreeSurfer. A morphing map is a linear mapping from cortical
surface values in subject A (:math:`x^{(A)}`) to those in another subject B
(:math:`x^{(B)}`)

.. math::    x^{(B)} = M^{(AB)} x^{(A)}\ ,

where :math:`M^{(AB)}` is a sparse matrix with at most three nonzero elements
on each row. These elements are determined as follows. First, using the aligned
spherical surfaces, for each vertex :math:`x_j^{(B)}`, find the triangle
:math:`T_j^{(A)}` on the spherical surface of subject A which contains the
location :math:`x_j^{(B)}`. Next, find the numbers of the vertices of this
triangle and set the corresponding elements on the *j* th row of
:math:`M^{(AB)}` so that :math:`x_j^{(B)}` will be a linear interpolation
between the triangle vertex values reflecting the location :math:`x_j^{(B)}`
within the triangle :math:`T_j^{(A)}`.

It follows from the above definition that in general

.. math::    M^{(AB)} \neq (M^{(BA)})^{-1}\ ,

*i.e.*,

.. math::    x_{(A)} \neq M^{(BA)} M^{(AB)} x^{(A)}\ ,

even if

.. math::    x^{(A)} \approx M^{(BA)} M^{(AB)} x^{(A)}\ ,

*i.e.*, the mapping is *almost* a bijection.


About smoothing
~~~~~~~~~~~~~~~

The current estimates are normally defined only in a decimated grid which is a
sparse subset of the vertices in the triangular tessellation of the cortical
surface. Therefore, any sparse set of values is distributed to neighboring
vertices to make the visualized results easily understandable. This procedure
has been traditionally called smoothing but a more appropriate name might be
smudging or blurring in accordance with similar operations in image processing
programs.

In MNE software terms, smoothing of the vertex data is an iterative procedure,
which produces a blurred image :math:`x^{(N)}` from the original sparse image
:math:`x^{(0)}` by applying in each iteration step a sparse blurring matrix:

.. math::    x^{(p)} = S^{(p)} x^{(p - 1)}\ .

On each row :math:`j` of the matrix :math:`S^{(p)}` there are :math:`N_j^{(p -
1)}` nonzero entries whose values equal :math:`1/N_j^{(p - 1)}`. Here
:math:`N_j^{(p - 1)}` is the number of immediate neighbors of vertex :math:`j`
which had non-zero values at iteration step :math:`p - 1`. Matrix
:math:`S^{(p)}` thus assigns the average of the non-zero neighbors as the new
value for vertex :math:`j`. One important feature of this procedure is that it
tends to preserve the amplitudes while blurring the surface image.

Once the indices non-zero vertices in :math:`x^{(0)}` and the topology of the
triangulation are fixed the matrices :math:`S^{(p)}` are fixed and independent
of the data. Therefore, it would be in principle possible to construct a
composite blurring matrix

.. math::    S^{(N)} = \prod_{p = 1}^N {S^{(p)}}\ .

However, it turns out to be computationally more effective to do blurring with
an iteration. The above formula for :math:`S^{(N)}` also shows that the
smudging (smoothing) operation is linear.
:orphan:

Institutional partners
----------------------

.. NOTE: this file is included in doc/funding.rst and doc/overview/people.rst.
   Changes here are reflected there. If you want to link to this content, link
   to :ref:`supporting-institutions` to link to that section of the funding.rst
   page. The next line is a target for :start-after: so we can omit the title
   from the include:
   institutional-partners-begin-content

Current partners
~~~~~~~~~~~~~~~~

- `Aalto-yliopiston perustieteiden korkeakoulu <https://sci.aalto.fi/>`_
- `Aarhus Universitet <https://www.au.dk/>`_
- `Athinoula A. Martinos Center for Biomedical Imaging <https://martinos.org/>`_
- `Children’s Hospital of Philadelphia Research Institute <https://imaging.research.chop.edu/>`_
- `Harvard Medical School <https://hms.harvard.edu/>`_
- `Institut national de recherche en informatique et en automatique <https://www.inria.fr/>`_
- `Karl-Franzens-Universität Graz <https://www.uni-graz.at/>`_
- `Macquarie University <https://www.mq.edu.au/>`_
- `Massachusetts General Hospital <https://www.massgeneral.org/>`_
- `Max-Planck-Institut für Bildungsforschung <https://www.mpib-berlin.mpg.de/>`_
- `SWPS Uniwersytet Humanistycznospołeczny <https://www.swps.pl/>`_
- `University of Washington <https://www.washington.edu/>`_

Former partners
~~~~~~~~~~~~~~~

- `Berkeley Institute for Data Science <https://bids.berkeley.edu/>`_
- `Boston University <https://www.bu.edu/>`_
- `Commissariat à l’énergie atomique et aux énergies alternatives <http://www.cea.fr/>`_
- `Forschungszentrum Jülich <https://www.fz-juelich.de/>`_
- `Institut du Cerveau et de la Moelle épinière <https://icm-institute.org/>`_
- `Institut national de la santé et de la recherche médicale <https://www.inserm.fr/>`_
- `Massachusetts Institute of Technology <https://web.mit.edu/>`_
- `New York University <https://www.nyu.edu/>`_
- `Technische Universität Ilmenau <https://www.tu-ilmenau.de/>`_
- `Télécom ParisTech <https://www.telecom-paris.fr/>`_
:orphan:

The forward solution
====================

This page covers the definitions of different coordinate systems employed in
MNE software and FreeSurfer, the details of the computation of the forward
solutions, and the associated low-level utilities.

.. NOTE: part of this file is included in doc/overview/implementation.rst.
   Changes here are reflected there. If you want to link to this content, link
   to :ref:`ch_forward` to link to that section of the implementation.rst page.
   The next line is a target for :start-after: so we can omit the title from
   the include:
   forward-begin-content


.. _coordinate_systems:

MEG/EEG and MRI coordinate systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. sidebar:: Coordinate systems in MNE-Python

   In some MNE-Python objects (e.g., :class:`~mne.Forward`,
   :class:`~mne.SourceSpaces`, etc), information about the coordinate frame is
   encoded as a constant integer value. The meaning of those integers is
   determined `in the source code
   <https://github.com/mne-tools/mne-python/blob/main/mne/io/constants.py#L186-L197>`__.

The coordinate systems used in MNE software (and FreeSurfer) and their
relationships are depicted in :ref:`coordinate_system_figure`. Except for the
*sensor coordinates*, all of the coordinate systems are Cartesian and have the
"RAS" (Right-Anterior-Superior) orientation, *i.e.*, the :math:`x` axis points
to the right, the :math:`y` axis to the front, and the :math:`z` axis up.

.. _coordinate_system_figure:

.. figure:: ../_static/CoordinateSystems.png
    :alt: MEG/EEG and MRI coordinate systems

    MEG/EEG and MRI coordinate systems

    The coordinate transforms present in the fif files in MNE and the
    FreeSurfer files as well as those set to fixed values are indicated with
    :math:`T_x`, where :math:`x` identifies the transformation.

The coordinate systems related to MEG/EEG data are:

**Head coordinates**

    This is a coordinate system defined with help of the fiducial landmarks
    (nasion and the two auricular points). In fif files, EEG electrode
    locations are given in this coordinate system. In addition, the head
    digitization data acquired in the beginning of an MEG, MEG/EEG, or EEG
    acquisition are expressed in head coordinates. For details, see
    :ref:`coordinate_systems`.

**Device coordinates**

    This is a coordinate system tied to the MEG device. The relationship of the
    Device and Head coordinates is determined during an MEG measurement by
    feeding current to three to five head-position indicator (HPI) coils and by
    determining their locations with respect to the MEG sensor array from the
    magnetic fields they generate.

**Sensor coordinates**

    Each MEG sensor has a local coordinate system defining the orientation and
    location of the sensor. With help of this coordinate system, the numerical
    integration data needed for the computation of the magnetic field can be
    expressed conveniently as discussed in :ref:`coil_geometry_information`.
    The channel information data in the fif files contain the information to
    specify the coordinate transformation between the coordinates of each
    sensor and the MEG device coordinates.

The coordinate systems related to MRI data are:

**Surface RAS coordinates**

    The FreeSurfer surface data are expressed in this coordinate system. The
    origin of this coordinate system is at the center of the conformed
    FreeSurfer MRI volumes (usually 256 x 256 x 256 isotropic 1-mm3  voxels)
    and the axes are oriented along the axes of this volume. The BEM surface
    and the locations of the sources in the source space are usually expressed
    in this coordinate system in the fif files. In this manual, the *Surface
    RAS coordinates* are usually referred to as *MRI coordinates* unless there
    is need to specifically discuss the different MRI-related coordinate
    systems.

**RAS coordinates**

    This coordinate system has axes identical to the Surface RAS coordinates
    but the location of the origin is different and defined by the original MRI
    data, i.e. , the origin is in a scanner-dependent location. There is hardly
    any need to refer to this coordinate system explicitly in the analysis with
    the MNE software. However, since the Talairach coordinates, discussed
    below, are defined with respect to *RAS coordinates* rather than the
    *Surface RAS coordinates*, the RAS coordinate system is implicitly involved
    in the transformation between Surface RAS coordinates and the two
    *Talairach* coordinate systems.

**MNI Talairach coordinates**

    The definition of this coordinate system is discussed, e.g., in
    https://imaging.mrc-cbu.cam.ac.uk/imaging/MniTalairach. This transformation
    is determined during the FreeSurfer reconstruction process. These
    coordinates are in MNI305 space.

**FreeSurfer Talairach coordinates**

    The problem with the MNI Talairach coordinates is that the linear MNI
    Talairach transform does not match the brains completely to the Talairach
    brain. This is probably because the Talairach atlas brain is a rather odd
    shape, and as a result, it is difficult to match a standard brain to the
    atlas brain using an affine transform. As a result, the MNI brains are
    slightly larger (in particular higher, deeper and longer) than the
    Talairach brain. The differences are larger as you get further from the
    middle of the brain, towards the outside. The FreeSurfer Talairach
    coordinates mitigate this problem by additing a an additional
    transformation, defined separately for negatice and positive MNI Talairach
    :math:`z` coordinates. These two transformations, denoted by :math:`T_-`
    and :math:`T_+` in :ref:`coordinate_system_figure`, are fixed as discussed in
    https://imaging.mrc-cbu.cam.ac.uk/imaging/MniTalairach (*Approach 2*).

The different coordinate systems are related by coordinate transformations
depicted in :ref:`coordinate_system_figure`. The arrows and coordinate
transformation symbols (:math:`T_x`) indicate the transformations actually
present in the FreeSurfer files. Generally,

.. math::    \begin{bmatrix}
		x_2 \\
		y_2 \\
		z_2 \\
		1
	        \end{bmatrix} = T_{12} \begin{bmatrix}
		x_1 \\
		y_1 \\
		z_1 \\
		1
	        \end{bmatrix} = \begin{bmatrix}
		R_{11} & R_{12} & R_{13} & x_0 \\
		R_{21} & R_{22} & R_{23} & y_0 \\
		R_{31} & R_{32} & R_{33} & z_0 \\
		0 & 0 & 0 & 1
	        \end{bmatrix} \begin{bmatrix}
		x_1 \\
		y_1 \\
		z_1 \\
		1
	        \end{bmatrix}\ ,

where :math:`x_k`, :math:`y_k`,and :math:`z_k` are the location coordinates in
two coordinate systems, :math:`T_{12}` is the coordinate transformation from
coordinate system "1" to "2", :math:`x_0`, :math:`y_0`, and :math:`z_0` is the
location of the origin of coordinate system "1" in coordinate system "2", and
:math:`R_{jk}` are the elements of the rotation matrix relating the two
coordinate systems. The coordinate transformations are present in different
files produced by FreeSurfer and MNE.
The fixed transformations :math:`T_-` and :math:`T_+` are:

.. math::    T_{-} = \begin{bmatrix}
		0.99 & 0 & 0 & 0 \\
		0 & 0.9688 & 0.042 & 0 \\
		0 & -0.0485 & 0.839 & 0 \\
		0 & 0 & 0 & 1
	        \end{bmatrix}

and

.. math::    T_{+} = \begin{bmatrix}
		0.99 & 0 & 0 & 0 \\
		0 & 0.9688 & 0.046 & 0 \\
		0 & -0.0485 & 0.9189 & 0 \\
		0 & 0 & 0 & 1
	        \end{bmatrix}

.. note::
   This section does not discuss the transformation between the MRI voxel
   indices and the different MRI coordinates. However, it is important to note
   that in FreeSurfer, MNE, as well as in Neuromag software an integer voxel
   coordinate corresponds to the location of the center of a voxel. Detailed
   information on the FreeSurfer MRI systems can be found at
   https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems.
   The symbols :math:`T_x` are defined in :ref:`coordinate_system_figure`.

.. tabularcolumns:: |p{0.2\linewidth}|p{0.3\linewidth}|p{0.5\linewidth}|
.. table:: Coordinate transformations in FreeSurfer and MNE software packages.

    +------------------------------+-------------------------------+-------------------------------------------------+
    | Transformation               | FreeSurfer                    | MNE                                             |
    +------------------------------+-------------------------------+-------------------------------------------------+
    | :math:`T_1`                  | Not present                   | | Measurement data files                        |
    |                              |                               | | Forward solution files (``*fwd.fif``)         |
    |                              |                               | | Inverse operator files (``*inv.fif``)         |
    +------------------------------+-------------------------------+-------------------------------------------------+
    | :math:`T_{s_1}\dots T_{s_n}` | Not present                   | Channel information in files                    |
    |                              |                               | containing :math:`T_1`.                         |
    +------------------------------+-------------------------------+-------------------------------------------------+
    | :math:`T_2`                  | Not present                   | | MRI description filesSeparate                 |
    |                              |                               | | Separate ``-trans.fif`` files                 |
    |                              |                               | | from :ref:`mne coreg`                         |
    |                              |                               | | Forward solution files                        |
    |                              |                               | | Inverse operator files                        |
    +------------------------------+-------------------------------+-------------------------------------------------+
    | :math:`T_3`                  | ``mri/*mgz`` files            | :class:`nibabel.freesurfer.mghformat.MGHImage`  |
    +------------------------------+-------------------------------+-------------------------------------------------+
    | :math:`T_4`                  | mri/transforms/talairach.xfm  | Internal reading                                |
    +------------------------------+-------------------------------+-------------------------------------------------+
    | :math:`T_-`                  | Hardcoded in software         | Hardcoded in software.                          |
    +------------------------------+-------------------------------+-------------------------------------------------+
    | :math:`T_+`                  | Hardcoded in software         | Hardcoded in software.                          |
    +------------------------------+-------------------------------+-------------------------------------------------+

.. _head_device_coords:

The head and device coordinate systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: ../_static/HeadCS.png
    :alt: Head coordinate system

    The head coordinate system

The MEG/EEG head coordinate system employed in the MNE software is a
right-handed Cartesian coordinate system. The direction of :math:`x` axis is
from left to right, that of :math:`y` axis to the front, and the :math:`z` axis
thus points up.

The :math:`x` axis of the head coordinate system passes through the two
periauricular or preauricular points digitized before acquiring the data with
positive direction to the right. The :math:`y` axis passes through the nasion
and is normal to the :math:`x` axis. The :math:`z` axis points up according to
the right-hand rule and is normal to the :math:`xy` plane.

The origin of the MEG device coordinate system is device dependent. Its origin
is located approximately at the center of a sphere which fits the occipital
section of the MEG helmet best with :math:`x` axis axis going from left to
right and :math:`y` axis pointing front. The :math:`z` axis is, again, normal
to the :math:`xy` plane with positive direction up.

.. note::
   The above definition is identical to that of the Neuromag MEG/EEG (head)
   coordinate system. However, in 4-D Neuroimaging and CTF MEG systems the head
   coordinate frame definition is different. The origin of the coordinate
   system is at the midpoint of the left and right auricular points. The
   :math:`x` axis passes through the nasion and the origin with positive
   direction to the front. The :math:`y` axis is perpendicular to the :math:`x`
   axis on the and lies in the plane defined by the three fiducial landmarks,
   positive direction from right to left. The :math:`z` axis is normal to the
   plane of the landmarks, pointing up. Note that in this convention the
   auricular points are not necessarily located on :math:`y` coordinate axis.
   The file conversion utilities take care of these idiosyncrasies and convert
   all coordinate information to the MNE software head coordinate frame.

Creating a surface-based source space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fif format source space files containing the dipole locations and
orientations are created with :func:`mne.setup_source_space`.

Creating a volumetric or discrete source space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to source spaces confined to a surface, the MNE software provides
some support for three-dimensional source spaces bounded by a surface as well
as source spaces comprised of discrete, arbitrarily located source points. The
:func:`mne.setup_volume_source_space` utility assists in generating such source
spaces.

Creating the BEM meshes
~~~~~~~~~~~~~~~~~~~~~~~

See :ref:`bem-model`.

Topology checks
---------------

The following topology checks are performed during the creation of BEM models:

- The completeness of each surface is confirmed by calculating the total solid
  angle subtended by all triangles from a point inside the triangulation. The
  result should be very close to :math:`4 \pi`. If the result is :math:`-4 \pi`
  instead, it is conceivable that the ordering of the triangle vertices is
  incorrect and the ``--swap`` option should be specified.

- The correct ordering of the surfaces is verified by checking that the
  surfaces are inside each other as expected. This is accomplished by checking
  that the sum solid angles subtended by triangles of a surface :math:`S_k` at
  all vertices of another surface :math:`S_p` which is supposed to be inside it
  equals :math:`4 \pi`. Naturally, this check is applied only if the model has
  more than one surface. Since the surface relations are transitive, it is
  enough to check that the outer skull surface is inside the skin surface and
  that the inner skull surface is inside the outer skull one.

- The extent of each of the triangulated volumes is checked. If the extent is
  smaller than 50mm, an error is reported. This may indicate that the vertex
  coordinates have been specified in meters instead of millimeters.


Computing the BEM geometry data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The utility :func:`mne.make_bem_solution` computes the geometry information for
BEM.

.. _coil_geometry_information:

Coil geometry information
~~~~~~~~~~~~~~~~~~~~~~~~~

This Section explains the presentation of MEG detection coil geometry
information the approximations used for different detection coils in MNE
software. Two pieces of information are needed to characterize the detectors:

- The location and orientation a local coordinate system for each detector.

- A unique identifier, which has an one-to-one correspondence to the
  geometrical description of the coil.

.. note:: MNE ships with several coil geometry configurations. They can be
          found in ``mne/data``. See :ref:`ex-plot-meg-sensors` for a
          comparison between different coil geometries, and
          :ref:`implemented_coil_geometries` for detailed information regarding
          the files describing Neuromag coil geometries.


The sensor coordinate system
----------------------------

The sensor coordinate system is completely characterized by the location of its
origin and the direction cosines of three orthogonal unit vectors pointing to
the directions of the x, y, and z axis. In fact, the unit vectors contain
redundant information because the orientation can be uniquely defined with
three angles. The measurement fif files list these data in MEG device
coordinates. Transformation to the MEG head coordinate frame can be easily
accomplished by applying the device-to-head coordinate transformation matrix
available in the data files provided that the head-position indicator was used.
Optionally, the MNE software forward calculation applies another coordinate
transformation to the head-coordinate data to bring the coil locations and
orientations to the MRI coordinate system.

If :math:`r_0` is a row vector for the origin of the local sensor coordinate
system and :math:`e_x`, :math:`e_y`, and :math:`e_z` are the row vectors for
the three orthogonal unit vectors, all given in device coordinates, a location
of a point :math:`r_C` in sensor coordinates is transformed to device
coordinates (:math:`r_D`) by

.. math::    [r_D 1] = [r_C 1] T_{CD}\ ,

where

.. math::    T = \begin{bmatrix}
		e_x & 0 \\
		e_y & 0 \\
		e_z & 0 \\
		r_{0D} & 1
	        \end{bmatrix}\ .

Calculation of the magnetic field
---------------------------------

The forward calculation in the MNE software computes the signals detected by
each MEG sensor for three orthogonal dipoles at each source space location.
This requires specification of the conductor model, the location and
orientation of the dipoles, and the location and orientation of each MEG sensor
as well as its coil geometry.

The output of each SQUID sensor is a weighted sum of the magnetic fluxes
threading the loops comprising the detection coil. Since the flux threading a
coil loop is an integral of the magnetic field component normal to the coil
plane, the output of the k :sup:`th` MEG channel, :math:`b_k` can be
approximated by:

.. math::    b_k = \sum_{p = 1}^{N_k} {w_{kp} B(r_{kp}) \cdot n_{kp}}

where :math:`r_{kp}` are a set of :math:`N_k` integration points covering the
pickup coil loops of the sensor, :math:`B(r_{kp})` is the magnetic field due to
the current sources calculated at :math:`r_{kp}`, :math:`n_{kp}` are the coil
normal directions at these points, and :math:`w_{kp}` are the weights
associated to the integration points. This formula essentially presents
numerical integration of the magnetic field over the pickup loops of sensor
:math:`k`.

There are three accuracy levels for the numerical integration expressed above.
The *simple* accuracy means the simplest description of the coil. This accuracy
is not used in the MNE forward calculations. The *normal* or *recommended*
accuracy typically uses two integration points for planar gradiometers, one in
each half of the pickup coil and four evenly distributed integration points for
magnetometers. This is the default accuracy used by MNE. If the ``--accurate``
option is specified, the forward calculation typically employs a total of eight
integration points for planar gradiometers and sixteen for magnetometers.
Detailed information about the integration points is given in the next section.


.. _implemented_coil_geometries:

Implemented coil geometries
---------------------------

This section describes the coil geometries currently implemented
in MNE. The coil types fall in two general categories:

- Axial gradiometers and planar gradiometers
  and

- Planar magnetometers.

For axial sensors, the *z* axis of the local coordinate system is parallel to
the field component detected, *i.e.*, normal to the coil plane.For circular
coils, the orientation of the *x* and *y* axes on the plane normal to the z
axis is irrelevant. In the square coils employed in the Vectorview (TM) system
the *x* axis is chosen to be parallel to one of the sides of the magnetometer
coil. For planar sensors, the *z* axis is likewise normal to the coil plane and
the x axis passes through the centerpoints of the two coil loops so that the
detector gives a positive signal when the normal field component increases
along the *x* axis.

:ref:`normal_coil_descriptions` lists the parameters of the *normal* coil
geometry descriptions :ref:`accurate_coil_descriptions` lists the *accurate*
descriptions. For simple accuracy, please consult the coil definition file, see
:ref:`coil_definition_file`. The columns of the tables contain the following
data:

- The number identifying the coil id.
  This number is used in the coil descriptions found in the FIF files.

- Description of the coil.

- Number of integration points used

- The locations of the integration points in sensor coordinates.

- Weights assigned to the field values at the integration points.
  Some formulas are listed instead of the numerical values to demonstrate
  the principle of the calculation. For example, in the normal coil
  descriptions of the planar gradiometers the weights are inverses
  of the baseline of the gradiometer to show that the output is in
  T/m.

.. note:: The coil geometry information is stored in the file
          :file:`mne/data/coil_def.dat`, which is
          automatically created by the MNE-C utility ``mne_list_coil_def``.

.. tabularcolumns:: |p{0.1\linewidth}|p{0.3\linewidth}|p{0.1\linewidth}|p{0.25\linewidth}|p{0.2\linewidth}|
.. _normal_coil_descriptions:
.. table:: Normal coil descriptions.

    +------+-------------------------+----+----------------------------------+----------------------+
    | Id   | Description             | n  | r/mm                             | w                    |
    +======+=========================+====+==================================+======================+
    | 2    | Neuromag-122            | 2  | (+/-8.1, 0, 0) mm                | +/-1 ⁄ 16.2mm        |
    |      | planar gradiometer      |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 2000 | A point magnetometer    | 1  | (0, 0, 0)mm                      | 1                    |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 3012 | Vectorview type 1       | 2  | (+/-8.4, 0, 0.3) mm              | +/-1 ⁄ 16.8mm        |
    |      | planar gradiometer      |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 3013 | Vectorview type 2       | 2  | (+/-8.4, 0, 0.3) mm              | +/-1 ⁄ 16.8mm        |
    |      | planar gradiometer      |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 3022 | Vectorview type 1       | 4  | (+/-6.45, +/-6.45, 0.3)mm        | 1/4                  |
    |      | magnetometer            |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 3023 | Vectorview type 2       | 4  | (+/-6.45, +/-6.45, 0.3)mm        | 1/4                  |
    |      | magnetometer            |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 3024 | Vectorview type 3       | 4  | (+/-5.25, +/-5.25, 0.3)mm        | 1/4                  |
    |      | magnetometer            |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 2000 | An ideal point          | 1  | (0.0, 0.0, 0.0)mm                | 1                    |
    |      | magnetometer            |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 4001 | Magnes WH               | 4  | (+/-5.75, +/-5.75, 0.0)mm        | 1/4                  |
    |      | magnetometer            |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 4002 | Magnes WH 3600          | 8  | (+/-4.5, +/-4.5, 0.0)mm          | 1/4                  |
    |      | axial gradiometer       |    | (+/-4.5, +/-4.5, 50.0)mm         | -1/4                 |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 4003 | Magnes reference        | 4  | (+/-7.5, +/-7.5, 0.0)mm          | 1/4                  |
    |      | magnetometer            |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 4004 | Magnes reference        | 8  | (+/-20, +/-20, 0.0)mm            | 1/4                  |
    |      | gradiometer measuring   |    | (+/-20, +/-20, 135)mm            | -1/4                 |
    |      | diagonal gradients      |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 4005 | Magnes reference        | 8  | (87.5, +/-20, 0.0)mm             | 1/4                  |
    |      | gradiometer measuring   |    | (47.5, +/-20, 0.0)mm             | -1/4                 |
    |      | off-diagonal gradients  |    | (-87.5, +/-20, 0.0)mm            | 1/4                  |
    |      |                         |    | (-47.5, +/-20, 0.0)mm            | -1/4                 |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 5001 | CTF 275 axial           | 8  | (+/-4.5, +/-4.5, 0.0)mm          | 1/4                  |
    |      | gradiometer             |    | (+/-4.5, +/-4.5, 50.0)mm         | -1/4                 |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 5002 | CTF reference           | 4  | (+/-4, +/-4, 0.0)mm              | 1/4                  |
    |      | magnetometer            |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 5003 | CTF reference           | 8  | (+/-8.6, +/-8.6, 0.0)mm          | 1/4                  |
    |      | gradiometer measuring   |    | (+/-8.6, +/-8.6, 78.6)mm         | -1/4                 |
    |      | diagonal gradients      |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+

.. note:: If a plus-minus sign occurs in several coordinates, all possible
          combinations have to be included.

.. tabularcolumns:: |p{0.1\linewidth}|p{0.3\linewidth}|p{0.05\linewidth}|p{0.25\linewidth}|p{0.15\linewidth}|
.. _accurate_coil_descriptions:
.. table:: Accurate coil descriptions

    +------+-------------------------+----+----------------------------------+----------------------+
    | Id   | Description             | n  | r/mm                             | w                    |
    +======+=========================+====+==================================+======================+
    | 2    | Neuromag-122 planar     | 8  | +/-(8.1, 0, 0) mm                | +/-1 ⁄ 16.2mm        |
    |      | gradiometer             |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 2000 | A point magnetometer    | 1  | (0, 0, 0) mm                     | 1                    |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 3012 | Vectorview type 1       | 2  | (+/-8.4, 0, 0.3) mm              | +/-1 ⁄ 16.8mm        |
    |      | planar gradiometer      |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 3013 | Vectorview type 2       | 2  | (+/-8.4, 0, 0.3) mm              | +/-1 ⁄ 16.8mm        |
    |      | planar gradiometer      |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 3022 | Vectorview type 1       | 4  | (+/-6.45, +/-6.45, 0.3)mm        | 1/4                  |
    |      | magnetometer            |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 3023 | Vectorview type 2       | 4  | (+/-6.45, +/-6.45, 0.3)mm        | 1/4                  |
    |      | magnetometer            |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 3024 | Vectorview type 3       | 4  | (+/-5.25, +/-5.25, 0.3)mm        | 1/4                  |
    |      | magnetometer            |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 4001 | Magnes WH magnetometer  | 4  | (+/-5.75, +/-5.75, 0.0)mm        | 1/4                  |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 4002 | Magnes WH 3600          | 4  | (+/-4.5, +/-4.5, 0.0)mm          | 1/4                  |
    |      | axial gradiometer       |    | (+/-4.5, +/-4.5, 0.0)mm          | -1/4                 |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 4004 | Magnes reference        | 8  | (+/-20, +/-20, 0.0)mm            | 1/4                  |
    |      | gradiometer measuring   |    | (+/-20, +/-20, 135)mm            | -1/4                 |
    |      | diagonal gradients      |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 4005 | Magnes reference        | 8  | (87.5, +/-20, 0.0)mm             | 1/4                  |
    |      | gradiometer measuring   |    | (47.5, +/-20, 0.0)mm             | -1/4                 |
    |      | off-diagonal gradients  |    | (-87.5, +/-20, 0.0)mm            | 1/4                  |
    |      |                         |    | (-47.5, +/-20, 0.0)mm            | -1/4                 |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 5001 | CTF 275 axial           | 8  | (+/-4.5, +/-4.5, 0.0)mm          | 1/4                  |
    |      | gradiometer             |    | (+/-4.5, +/-4.5, 50.0)mm         | -1/4                 |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 5002 | CTF reference           | 4  | (+/-4, +/-4, 0.0)mm              | 1/4                  |
    |      | magnetometer            |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 5003 | CTF 275 reference       | 8  | (+/-8.6, +/-8.6, 0.0)mm          | 1/4                  |
    |      | gradiometer measuring   |    | (+/-8.6, +/-8.6, 78.6)mm         | -1/4                 |
    |      | diagonal gradients      |    |                                  |                      |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 5004 | CTF 275 reference       | 8  | (47.8, +/-8.5, 0.0)mm            | 1/4                  |
    |      | gradiometer measuring   |    | (30.8, +/-8.5, 0.0)mm            | -1/4                 |
    |      | off-diagonal gradients  |    | (-47.8, +/-8.5, 0.0)mm           | 1/4                  |
    |      |                         |    | (-30.8, +/-8.5, 0.0)mm           | -1/4                 |
    +------+-------------------------+----+----------------------------------+----------------------+
    | 6001 | MIT KIT system axial    | 8  | (+/-3.875, +/-3.875, 0.0)mm      | 1/4                  |
    |      | gradiometer             |    | (+/-3.875, +/-3.875, 0.0)mm      | -1/4                 |
    +------+-------------------------+----+----------------------------------+----------------------+


.. _coil_definition_file:

The coil definition file
------------------------

The coil geometry information is stored in the text file
:file:`{$MNE_ROOT}/share/mne/coil_def.dat`. In this file, any lines starting
with the pound sign (#) are comments. A coil definition starts with a
description line containing the following fields:

- :samp:`{<class>}`: A number indicating class of this coil.

- :samp:`{<id>}`: Coil ID value. This value is listed in the first column of
  Tables :ref:`normal_coil_descriptions` and :ref:`accurate_coil_descriptions`.

- :samp:`{<accuracy>}`: The coil representation accuracy. Possible values and
  their meanings are listed in :ref:`coil_accuracies`.

- :samp:`{<np>}`: Number of integration points in this representation.

- :samp:`{<size/m>}`: The size of the coil. For circular coils this is the
  diameter of the coil and for square ones the side length of the square. This
  information is mainly included to facilitate drawing of the coil geometry. It
  should not be employed to infer a coil approximation for the forward
  calculations.

- :samp:`{<baseline/m>}`: The baseline of a this kind of a coil. This will be
  zero for magnetometer coils. This information is mainly included to
  facilitate drawing of the coil geometry. It should not be employed to infer
  a coil approximation for the forward calculations.

- :samp:`{<description>}`: Short description of this kind of a coil. If the
  description contains several words, it is enclosed in quotes.


.. tabularcolumns:: |p{0.1\linewidth}|p{0.5\linewidth}|
.. _coil_accuracies:
.. table:: Coil representation accuracies.

    =======  ====================================================================================
    Value    Meaning
    =======  ====================================================================================
    1        The simplest representation available
    2        The standard or *normal* representation (see :ref:`normal_coil_descriptions`)
    3        The most *accurate* representation available (see :ref:`accurate_coil_descriptions`)
    =======  ====================================================================================

Each coil description line is followed by one or more integration point lines,
consisting of seven numbers:

- :samp:`{<weight>}`: Gives the weight for this integration point (last column
  in Tables :ref:`normal_coil_descriptions` and
  :ref:`accurate_coil_descriptions`).

- :samp:`{<x/m>} {<y/m>} {<z/m>}`: Indicates the location of the integration
  point (fourth column in Tables :ref:`normal_coil_descriptions` and
  :ref:`accurate_coil_descriptions`).

- :samp:`{<nx>} {<ny>} {<nz>}`: Components of a unit vector indicating the
  field component to be selected. Note that listing a separate unit vector for
  each integration points allows the implementation of curved coils and coils
  with the gradiometer loops tilted with respect to each other.


Computing the forward solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Purpose
-------

Examples on how to compute the forward solution in MNE-Python using
:func:`mne.make_forward_solution` can be found
:ref:`plot_forward_compute_forward_solution` and
:ref:`computing_the_forward_solution`.

Implementation of software gradient compensation
------------------------------------------------

Accounting for noise cancellation in MNE-Python is accomplished in
:meth:`mne.io.Raw.apply_gradient_compensation`. See
:ref:`plot_brainstorm_phantom_ctf` for an example.

CTF and 4D Neuroimaging data may have been subjected to noise cancellation
employing the data from the reference sensor array. Even though these sensor
are rather far away from the brain sources, :func:`mne.make_forward_solution`
takes them into account in the computations. If the data file has software
gradient compensation activated, it computes the field of at the reference
sensors in addition to the main MEG sensor array and computes a compensated
forward solution.

The EEG sphere model definition file
------------------------------------

In MNE-Python, different sphere models can be specified through
:func:`mne.make_sphere_model`. The default model has the following structure:

.. tabularcolumns:: |p{0.1\linewidth}|p{0.25\linewidth}|p{0.2\linewidth}|
.. table:: Structure of the default EEG model

    ========  =======================  =======================
    Layer     Relative outer radius    :math:`\sigma` (S/m)
    ========  =======================  =======================
    Head      1.0                      0.33
    Skull     0.97                     0.04
    CSF       0.92                     1.0
    Brain     0.90                     0.33
    ========  =======================  =======================

Although it is not BEM model per se the ``sphere`` structure describes the head
geometry so it can be passed as ``bem`` parameter in MNE-Python functions such
as :func:`mne.fit_dipole`, :func:`mne.viz.plot_alignment` or
:func:`mne.make_forward_solution`.

.. _eeg_sphere_model:

EEG forward solution in the sphere model
----------------------------------------

.. sidebar:: Sphere-model examples in MNE-Python

   For examples of using the sphere model when computing the forward model
   (using :func:`mne.make_forward_solution`), see :ref:`Brainstorm CTF phantom
   dataset tutorial <plt_brainstorm_phantom_ctf_eeg_sphere_geometry>`,
   :ref:`Brainstorm Elekta phantom dataset tutorial
   <plt_brainstorm_phantom_elekta_eeg_sphere_geometry>`, and
   :ref:`tut-source-alignment-without-mri`.

When the sphere model is employed, the computation of the EEG solution can be
substantially accelerated by using approximation methods described by Mosher
:footcite:`MosherEtAl1999`, Zhang :footcite:`Zhang1995`, and Berg
:footcite:`BergScherg1994`.
:func:`mne.make_forward_solution` approximates the solution with three dipoles
in a homogeneous sphere whose locations and amplitudes are determined by
minimizing the cost function:

.. math::
   S(r_1,\dotsc,r_m\ ,\ \mu_1,\dotsc,\mu_m) = \int_{scalp} {(V_{true} - V_{approx})}\,dS

where :math:`r_1,\dotsc,r_m` and :math:`\mu_1,\dotsc,\mu_m` are the locations
and amplitudes of the approximating dipoles and :math:`V_{true}` and
:math:`V_{approx}` are the potential distributions given by the true and
approximative formulas, respectively. It can be shown that this integral can be
expressed in closed form using an expansion of the potentials in spherical
harmonics. The formula is evaluated for the most superficial dipoles, *i.e.*,
those lying just inside the inner skull surface.

Averaging forward solutions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

One possibility to make a grand average over several runs of a experiment is to
average the data across runs and average the forward solutions accordingly. For
this purpose, :func:`mne.average_forward_solutions` computes a weighted average
of several forward solutions. The function averages both MEG and EEG forward
solutions. Usually the EEG forward solution is identical across runs because
the electrode locations do not change.

.. target for :end-before: forward-end-content
:orphan:

Creating the BEM meshes
=======================

.. NOTE: part of this file is included in doc/overview/implementation.rst.
   Changes here are reflected there. If you want to link to this content, link
   to :ref:`bem-model` to link to that section of the implementation.rst page.
   The next line is a target for :start-after: so we can omit the title from
   the include:
   bem-begin-content

.. _bem_watershed_algorithm:

Using the watershed algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The watershed algorithm [Segonne *et al.*,
2004] is part of the FreeSurfer software.
The name of the program is mri_watershed .
Its use in the MNE environment is facilitated by the script
:ref:`mne watershed_bem`.

After ``mne watershed_bem`` has completed, the following files appear in the
subject's :file:`bem/watershed` directory:

- :file:`{<subject>}_brain_surface` contains the brain surface triangulation.

- :file:`{<subject>}_inner_skull_surface` contains the inner skull
  triangulation.

- :file:`{<subject>}_outer_skull_surface` contains the outer skull
  triangulation.

- :file:`{<subject>}_outer_skin_surface` contains the scalp triangulation.

All of these surfaces are in the FreeSurfer format. In addition, there will be
a file called :file:`bem/watershed/ws.mgz` which contains the brain MRI
volume. Furthermore, ``mne watershed_bem`` script converts the scalp surface to
fif format and saves the result to :file:`bem/{<subject>}-head.fif`.


Using FLASH images
~~~~~~~~~~~~~~~~~~

This method depends on the availablily of MRI data acquired with a multi-echo
FLASH sequence at two flip angles (5 and 30 degrees). These data can be
acquired separately from the MPRAGE data employed in FreeSurfer cortical
reconstructions but it is strongly recommended that they are collected at the
same time with the MPRAGEs or at least with the same scanner. For easy
co-registration, the images should have FOV, matrix, slice thickness, gap, and
slice orientation as the MPRAGE data. For information on suitable pulse
sequences, see :footcite:`FischlEtAl2004`.

Creation of the BEM meshes using this method involves the following steps:

- Creating a synthetic 5-degree flip angle FLASH volume, register
  it with the MPRAGE data, and run the segmentation and meshing program.
  This step is accomplished by running the script :ref:`mne flash_bem`.

- Inspecting the meshes with tkmedit, see :ref:`inspecting-meshes`.

.. note:: Different methods can be employed for the creation of the
          individual surfaces. For example, it may turn out that the
          watershed algorithm produces are better quality skin surface than
          the segmentation approach based on the FLASH images. If this is
          the case, ``outer_skin.surf`` can set to point to the corresponding
          watershed output file while the other surfaces can be picked from
          the FLASH segmentation data.


Organizing MRI data into directories
------------------------------------

Since all images comprising the multi-echo FLASH data are contained in a single
series, it is necessary to organize the images according to the echoes before
proceeding to the BEM surface reconstruction. This can be accomplished by using
`dcm2niix <https://www.nitrc.org/plugins/mwiki/index.php/dcm2nii:MainPage>`__
or the MNE-C tool ``mne_organize_dicom`` if necessary, then use
:func:`mne.bem.convert_flash_mris`.

Creating the surface tessellations
----------------------------------

The BEM surface segmentation and tessellation is automated with the script
:ref:`mne flash_bem`. It assumes that a FreeSurfer reconstruction for this
subject is already in place.

Before running :ref:`mne flash_bem` do the following:

- Create symbolic links from the directories containing the 5-degree and
  30-degree flip angle FLASH series to ``flash05`` and ``flash30``,
  respectively:

  - :samp:`ln -s {<FLASH 5 series dir>} flash05`

  - :samp:`ln -s {<FLASH 30 series dir>} flash30`

- Some partition formats (e.g. FAT32) do not support symbolic links. In this
  case, copy the file to the appropriate series:

  - :samp:`cp {<FLASH 5 series dir>} flash05`

  - :samp:`cp {<FLASH 30 series dir>} flash30`

- Set the ``SUBJECTS_DIR`` and ``SUBJECT`` environment variables or pass
  the ``--subjects-dir`` and ``--subject`` options to ``mne flash_bem``

.. note:: If ``mne flash_bem`` is run with the ``--noflash30`` option, the
   :file:`flash30` directory is not needed, *i.e.*, only the 5-degree flip
   angle flash data are employed.

It may take a while for ``mne flash_bem`` to complete. It uses the FreeSurfer
directory structure under ``$SUBJECTS_DIR/$SUBJECT``. The script encapsulates
the following processing steps:

- It creates an mgz file corresponding to each of the eight echoes in each of
  the FLASH directories in ``mri/flash``. The files will be called
  :file:`mef {<flip-angle>}_{<echo-number>}.mgz`.

- If the ``unwarp=True`` option is specified, run grad_unwarp and produce
  files :file:`mef {<flip-angle>}_{<echo-number>}u.mgz`. These files will be
  then used in the following steps.

- It creates parameter maps in :file:`mri/flash/parameter_maps` using
  ``mri_ms_fitparms``.

- It creates a synthetic 5-degree flip angle volume in
  :file:`mri/flash/parameter_maps/flash5.mgz` using ``mri_synthesize``.

- Using ``fsl_rigid_register``, it creates a registered 5-degree flip angle
  volume ``mri/flash/parameter_maps/flash5_reg.mgz`` by registering
  :file:`mri/flash/parameter_maps/flash5.mgz` to the *T1* volume under ``mri``.

- Using ``mri_convert``, it converts the flash5_reg volume to COR format under
  ``mri/flash5``. If necessary, the T1 and brain volumes are also converted
  into the COR format.

- It runs ``mri_make_bem_surfaces`` to create the BEM surface tessellations.

- It creates the directory :file:`bem/flash`, moves the tri-format
  tringulations there and creates the corresponding FreeSurfer surface files
  in the same directory.

- The COR format volumes created by ``mne flash_bem`` are removed.

If the ``--noflash30`` option is specified to ``mne flash_bem``,
steps 3 and 4 in the above are replaced by averaging over the different
echo times in 5-degree flip angle data.

.. _inspecting-meshes:

Inspecting the meshes
---------------------

It is advisable to check the validity of the BEM meshes before
using them. This can be done with:

- the ``--view`` option of :ref:`mne flash_bem`
- calling :func:`mne.viz.plot_bem` directly
- Using FreeSurfer tools ``tkmedit`` or ``freeview``
:orphan:

Supported data formats
======================

.. NOTE: part of this file is included in doc/overview/implementation.rst.
   Changes here are reflected there. If you want to link to this content,
   link to :ref:`data-formats`. The next line is
   a target for :start-after: so we can omit the title above:
   data-formats-begin-content

When MNE-Python loads sensor data, the data are stored in a Python object of
type :class:`mne.io.Raw`. Specialized loading functions are provided for the
raw data file formats from a variety of equipment manufacturers. All raw data
input/output functions in MNE-Python are found in :mod:`mne.io` and start
with :samp:`read_raw_{*}`; see the documentation for each reader function for
more info on reading specific file types.

As seen in the table below, there are also a few formats defined by other
neuroimaging analysis software packages that are supported (EEGLAB,
FieldTrip). Like the equipment-specific loading functions, these will also
return an object of class :class:`~mne.io.Raw`; additional functions are
available for reading data that has already been epoched or averaged (see
table).

.. NOTE: To include only the table, here's a different target for :start-after:
   data-formats-begin-table

.. cssclass:: table-bordered
.. rst-class:: midvalign

============  ============================================  =========  ===================================
Data type     File format                                   Extension  MNE-Python function
============  ============================================  =========  ===================================
MEG           :ref:`Artemis123 <import-artemis>`            .bin       :func:`mne.io.read_raw_artemis123`

MEG           :ref:`4-D Neuroimaging / BTi <import-bti>`    <dir>      :func:`mne.io.read_raw_bti`

MEG           :ref:`CTF <import-ctf>`                       <dir>      :func:`mne.io.read_raw_ctf`

MEG and EEG   :ref:`Elekta Neuromag <import-neuromag>`      .fif       :func:`mne.io.read_raw_fif`

MEG           :ref:`KIT <import-kit>`                       .sqd       :func:`mne.io.read_raw_kit`,
                                                                       :func:`mne.read_epochs_kit`

MEG and EEG   :ref:`FieldTrip <import-fieldtrip>`           .mat       :func:`mne.io.read_raw_fieldtrip`,
                                                                       :func:`mne.read_epochs_fieldtrip`,
                                                                       :func:`mne.read_evoked_fieldtrip`

EEG           :ref:`Brainvision <import-bv>`                .vhdr      :func:`mne.io.read_raw_brainvision`

EEG           :ref:`Biosemi data format <import-biosemi>`   .bdf       :func:`mne.io.read_raw_bdf`

EEG           :ref:`Neuroscan CNT <import-cnt>`             .cnt       :func:`mne.io.read_raw_cnt`

EEG           :ref:`European data format <import-edf>`      .edf       :func:`mne.io.read_raw_edf`

EEG           :ref:`EEGLAB <import-set>`                    .set       :func:`mne.io.read_raw_eeglab`,
                                                                       :func:`mne.read_epochs_eeglab`

EEG           :ref:`EGI simple binary <import-egi>`         .egi       :func:`mne.io.read_raw_egi`

EEG           :ref:`EGI MFF format <import-mff>`            .mff       :func:`mne.io.read_raw_egi`

EEG           :ref:`eXimia <import-nxe>`                    .nxe       :func:`mne.io.read_raw_eximia`

EEG           :ref:`General data format <import-gdf>`       .gdf       :func:`mne.io.read_raw_gdf`

EEG           :ref:`Nicolet <import-nicolet>`               .data      :func:`mne.io.read_raw_nicolet`

EEG           :ref:`Persyst <import-persyst>`               .lay       :func:`mne.io.read_raw_persyst`

NIRS          :ref:`NIRx <import-nirx>`                     directory  :func:`mne.io.read_raw_nirx`

NIRS          :ref:`BOXY <import-boxy>`                     directory  :func:`mne.io.read_raw_boxy`
============  ============================================  =========  ===================================

More details are provided in the tutorials in the :ref:`tut-data-formats`
section.
:orphan:

The Signal-Space Projection (SSP) method
========================================

.. NOTE: part of this file is included in doc/overview/implementation.rst.
   Changes here are reflected there. If you want to link to this content, link
   to :ref:`ssp-method` to link to that section of the implementation.rst
   page. The next line is a target for :start-after: so we can omit the title
   from the include:
   ssp-begin-content

The Signal-Space Projection (SSP) is one approach to rejection of external
disturbances in software. The section presents some relevant details of this
method. For practical examples of how to use SSP for artifact rejection, see
:ref:`tut-artifact-ssp`.

General concepts
~~~~~~~~~~~~~~~~

Unlike many other noise-cancellation approaches, SSP does not require
additional reference sensors to record the disturbance fields. Instead, SSP
relies on the fact that the magnetic field distributions generated by the
sources in the brain have spatial distributions sufficiently different from
those generated by external noise sources. Furthermore, it is implicitly
assumed that the linear space spanned by the significant external noise patterns
has a low dimension.

Without loss of generality we can always decompose any :math:`n`-channel
measurement :math:`b(t)` into its signal and noise components as

.. math::    b(t) = b_s(t) + b_n(t)
   :label: additive_model

Further, if we know that :math:`b_n(t)` is well characterized by a few field
patterns :math:`b_1 \dotso b_m`, we can express the disturbance as

.. math::    b_n(t) = Uc_n(t) + e(t)\ ,
   :label: pca

where the columns of :math:`U` constitute an orthonormal basis for :math:`b_1
\dotso b_m`, :math:`c_n(t)` is an :math:`m`-component column vector, and the
error term :math:`e(t)` is small and does not exhibit any consistent spatial
distributions over time, *i.e.*, :math:`C_e = E \{e e^\top\} = I`. Subsequently,
we will call the column space of :math:`U` the noise subspace. The basic idea
of SSP is that we can actually find a small basis set :math:`b_1 \dotso b_m`
such that the conditions described above are satisfied. We can now construct
the orthogonal complement operator

.. math::    P_{\perp} = I - UU^\top
   :label: projector

and apply it to :math:`b(t)` in Equation :eq:`additive_model` yielding

.. math::    b_{s}(t) \approx P_{\perp}b(t)\ ,
   :label: result

since :math:`P_{\perp}b_n(t) = P_{\perp}(Uc_n(t) + e(t)) \approx 0` and
:math:`P_{\perp}b_{s}(t) \approx b_{s}(t)`. The projection operator
:math:`P_{\perp}` is called the **signal-space projection operator** and
generally provides considerable rejection of noise, suppressing external
disturbances by a factor of 10 or more. The effectiveness of SSP depends on two
factors:

- The basis set :math:`b_1 \dotso b_m` should be able to characterize the
  disturbance field patterns completely and

- The angles between the noise subspace space spanned by :math:`b_1 \dotso b_m`
  and the signal vectors :math:`b_s(t)` should be as close to :math:`\pi / 2`
  as possible.

If the first requirement is not satisfied, some noise will leak through because
:math:`P_{\perp}b_n(t) \neq 0`. If the any of the brain signal vectors
:math:`b_s(t)` is close to the noise subspace not only the noise but also the
signal will be attenuated by the application of :math:`P_{\perp}` and,
consequently, there might by little gain in signal-to-noise ratio.

Since the signal-space projection modifies the signal vectors originating in
the brain, it is necessary to apply the projection to the forward solution in
the course of inverse computations.

For more information on SSP, please consult the references listed in
:footcite:`TescheEtAl1995,UusitaloIlmoniemi1997`.

Estimation of the noise subspace
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As described above, application of SSP requires the estimation of the signal
vectors :math:`b_1 \dotso b_m` constituting the noise subspace. The most common
approach, also implemented in :func:`mne.compute_proj_raw`
is to compute a covariance matrix
of empty room data, compute its eigenvalue decomposition, and employ the
eigenvectors corresponding to the highest eigenvalues as basis for the noise
subspace. It is also customary to use a separate set of vectors for
magnetometers and gradiometers in the Vectorview system.

EEG average electrode reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The EEG average reference is the mean signal over all the sensors. It is
typical in EEG analysis to subtract the average reference from all the sensor
signals :math:`b^{1}(t), ..., b^{n}(t)`. That is:

.. math::	{b}^{j}_{s}(t) = b^{j}(t) - \frac{1}{n}\sum_{k}{b^k(t)}
   :label: eeg_proj

where the noise term :math:`b_{n}^{j}(t)` is given by

.. math:: 	b_{n}^{j}(t) = \frac{1}{n}\sum_{k}{b^k(t)}
   :label: noise_term

Thus, the projector vector :math:`P_{\perp}` will be given by
:math:`P_{\perp}=\frac{1}{n}[1, 1, ..., 1]`

.. warning::
   When applying SSP, the signal of interest can also be sometimes removed.
   Therefore, it's always a good idea to check how much the effect of interest
   is reduced by applying SSP. SSP might remove *both* the artifact and signal
   of interest.
{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autofunction:: {{ objname }}

.. _sphx_glr_backreferences_{{ fullname }}:

.. minigallery:: {{ fullname }}
    :add-heading:
{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :special-members: __contains__,__getitem__,__iter__,__len__,__add__,__sub__,__mul__,__div__,__neg__
   :members:

.. _sphx_glr_backreferences_{{ fullname }}:

.. minigallery:: {{ fullname }}
    :add-heading:
.. -*- mode: rst -*-


Documentation
=============

The icons are used in ``mne/viz/_brain/_brain.py`` for the toolbar.
It is necessary to compile those icons into a resource file for proper use by
the application.

The resource configuration file ``mne/icons/mne.qrc`` describes the location of
the resources in the filesystem and also defines aliases for their use in the code.

To automatically generate the resource file in ``mne/icons``:

.. code-block:: bash

    pyrcc5 -o mne/icons/resources.py mne/icons/mne.qrc

These Material design icons are provided by Google under the `Apache 2.0`_ license.


.. _Apache 2.0: https://github.com/google/material-design-icons/blob/master/LICENSE
