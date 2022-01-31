[![DOI](https://joss.theoj.org/papers/10.21105/joss.02404/status.svg)](https://doi.org/10.21105/joss.02404)

# seg1d

Automated segmentation of one-dimensional (1D) data


## Overview

seg1d is an open-source Python package for the automated segmentation of one-dimensional data using one or more reference segments. The segmentation process allows users to apply various methods and parameters for the process through weighted features (i.e., additional data attributed to the same set) in a rolling correlation size-varying window of a reference(s) to a target. Correlations can be averaged across the references and a peak detection algorithm finds the most prominent segments. Non-overlapping segments are identified and a clustering algorithm groups the most similar segments within the target. The package was developed for movement sciences but should be useful to anyone interested in extracting correlated subsequences from a dataset. 

![seg1d](https://raw.githubusercontent.com/cadop/seg1d/master/docs/build/plot_directive/api_basic-1.png)

Example of the segmentation algorithm using a portion of a sine wave. The initial reference segment is found, along with any additional segments. 

### Documentation

Full documentation available online: https://cadop.github.io/seg1d/


### Alternatives

There are existing libraries that provide clustering and similarity measures for querying subseries data. 

  * [The UCR Suite](https://www.cs.ucr.edu/~eamonn/UCRsuite.html) (Ultrafast subsequence search under both Dynamic Time Warping (DTW) and Euclidean Distance (ED))

  * [tslearn](https://github.com/rtavenar/tslearn)(Machine learning tools for the analysis of time series)
  
  * [seglearn](https://dmbee.github.io/seglearn/)(Provides a flexible approach to multivariate time series and related contextual (meta) data for classification, regression, and forecasting problems)

The advantage of seg1d is a simpler API which makes getting started and using the code quicker for non-experts when interested in purely extracting timestamps (segments) from a dataset. The API was built with motion capture data in mind, making the addition of features (e.g., marker trajectories) and sets of features (e.g., multiple subjects) easy. The output of seg1d is also geared towards users that are interested in comparing the segments found, such as returning masked arrays of segments. 


## Quickstart 


### Minimum Dependencies

Currently tested on ``Python 3.8`` on Ubuntu 18.04 and Windows 10. (Should work on ``Python 3.6`` and above)

Required Packages:

``numpy>=1.15``, ``scipy>=1.0.0``, ``sklearn>=0.2``, ``numba>=0.40``

For documentation:

``sphinx>=2``, ``numpydoc>=0.9.2``

For examples:

``matplotlib>=3.2.0``

### Installation

```pip install seg1d```


### Example usage

The documentation contains examples using data of both generated data (e.g., sine wave) and real-world examples (i.e., motion capture data). 

To quickly get started, try importing the seg1d module and using the provided sample data. 

```python

import seg1d 

#retrieve the sample reference, target, and weight data
r,t,w = seg1d.sampleData()
# define scaling percentage and rolling step size
minW, maxW, step  = 70, 150, 1 
#call the segmentation algorithm
a = seg1d.segment_data(r,t,w,minW,maxW,step)
print(a)
# Should output an array equal to:
# array([[207.       , 240.       ,   0.9124224],
#        [342.       , 381.       ,   0.8801901],
#        [ 72.       , 112.       ,   0.8776795]])
# where each array is of form [start index, end index, correlation]

```

For more examples, please refer to the full documention. 

## Project Info

To cite: 
Schwartz et al., (2020). seg1d: A Python package for Automated segmentation of one-dimensional (1D) data. Journal of Open Source Software, 5(52), 2404, https://doi.org/10.21105/joss.02404

Bibtex
```
@article{Schwartz2020,
  doi = {10.21105/joss.02404},
  url = {https://doi.org/10.21105/joss.02404},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {52},
  pages = {2404},
  author = {Mathew Schwartz and Todd C. Pataky and Cyril J. Donnelly},
  title = {seg1d: A Python package for Automated segmentation of one-dimensional (1D) data},
  journal = {Journal of Open Source Software}
}
```

This project was used for the following paper -  [Paper Link](https://commons.nmu.edu/isbs/vol38/iss1/231/) : 

Schwartz, Mathew; Pataky, Todd; Sui Geok Karen, CHUA; Wei Tech, ANG; and Donnelly, Cyril (2020) "AUTOMATED MULTI-FEATURE SEGMENTATION OF TREADMILL RUNNING," ISBS Proceedings Archive: Vol. 38 : Iss. 1 , Article 231.

### Community

Issues and feature requests should be submitted on [github](https://github.com/cadop/seg1d/issues). 

Please check the Community Guidelines for further details on contributing. 


### Credits

Code written and developed by Mathew Schwartz (New Jersey Institute of Technology).

Data used in sample provided by Precision Rehab, Rehabilitation Research Institute of Singapore.

Project was funded in part by the Agency for Science, Technology and Research (A\*STAR), Nanyang Technological University (NTU) and the National Health Group (NHG) (RRG3: 2019/19002).

### License

Please refer to the full [LICENSE](https://github.com/cadop/seg1d/blob/master/LICENSE.txt).
---
title: 'seg1d: A Python package for Automated segmentation of one-dimensional (1D) data'
tags:
  - Python
  - biomechanics
  - data series
  - movement
  - motion capture
authors:
  - name: Mathew Schwartz
    orcid: 0000-0003-3662-7203
    affiliation: 1
  - name: Todd C. Pataky
    orcid: 0000-0002-8292-7189
    affiliation: 2
  - name: Cyril J. Donnelly
    affiliation: 3
affiliations:
 - name: New Jersey Institute of Technology
   index: 1
 - name: Kyoto University, Department of Human Health Sciences
   index: 2
 - name: Nanyang Technological University, Rehabilitation Research Institute of Singapore
   index: 3
date: 25 March 2020
bibliography: paper.bib

---

# Summary

The term segmentation refers to the division of a data series into segments. In this paper, segments refer to a contiguous time series which is a subset of another time series [@bouchard2006automated]. Identifying these segments requires a reference series to be found in the larger series, which can be an exact copy or an approximation. When an approximation of the desired segment is used, the metric and method for identifying similarity of all subsets are essential to achieving the desired result. Additionally, multiple references may be desired to further describe the aspects of segment similarity, some of which may be more important than others. 

Segmenting human motion is a topic studied in various fields such as robotics (e.g., humanoids), biomechanics, and computer graphics (e.g., gaming and animation). In the human movement sciences, segmentation is often performed as the labeling of specific events, such as the components of a gait cycle, for analysis and interpretation of the data within these labeled bounds. Subjectively defining what a movement or a phase of a movement is can be particularly difficult due to variations in what one may define as a single movement. As such, the points at which the movement or phases of a movement starts and ends can be ambiguous. The averaging of multiple features (e.g., marker trajectories, joint angles, or other information derived from the data) of a movement or even multiple movements (e.g., multiple marker trajectories from multiple observation sets) allows for tolerance to some individual features failing to provide an expected characteristic (e.g., a signal above or below a threshold value from a force plate) that may normally be relied upon for identifying an event. Defining weights of individual features, with either algorithmic approaches or through expert knowledge, further facilitates segmentation of similar movements. 

With feature-rich data such as multiple marker trajectories from motion capture, the reduction of meaningful features is important for segmentation performance [@bouchard2015segmenting]. The use of marker trajectory and ground reaction force, without computing kinematics, has been shown to be sufficient in movement segmentation tasks [@lin2016comparison]. Unique feature creation, extraction, and storage can be used to additionally index databases for fast movement-based time-series data retrieval [@kapadia2013efficient]. Subseries searching of databases has often been performed with similarity metrics, each with their own individual downfalls, e.g., Longest Common Subseries (LCSS), Euclidean Distance (ED), and Dynamic Time Warping (DTW) [@vlachos2003indexing]. Discrete Fourier Transformation (DFT) has also been used as a method for improving the efficiency of windowed correlations [@zhu2002statstream]. Data reduction techniques for optimization can also be used in windowed correlations of generic time-series data [@cole2005fast].

Peak detection has been used to identify gait events through the inference of cyclic motion and reducing reliance on physical meanings of the signal by searching for the assumed cyclic pattern rather than a given threshold value or calculated joint angle [@jiang2017robust]. Peak detection from cross-correlated data has been used for gait event detection in accelerometry data [@yoneyama2013accelerometry]. Alternative methods for segmenting data have used physical devices to act as switches on foot contact [@agostini2013segmentation]. Most recently, deep neural networks have been used to predict foot contact but require a large amount (i.e., thousands of trials) of training data [@kidzinski2019automatic]. Using DTW, [@sarsfield2019segmentation] was able to identify movements in realtime with a single reference segment. Allowing both single and multiple reference segments, as well as multiple features and optional weights, the segmentation of various data is more practical in a single package.  

[//]: # Segmentation of movement data, when automated, is often built with specific features such as the velocity of a marker or kinematic data in software such as Visual3D. 


# Statement of Need

Subsequence identification and similarity between a reference(s) and target data items is a commonly desired task done both manually and automatically in a variety of fields. The ability to further automate and create reliable, consistent results, is of importance for many data processing related tasks. For example, in typical motion capture sessions of walking gait in a lab, embedded force plates provide high fidelity measurements for foot-strike and toe-off. However, the cycles before and after this event are often discarded. By using a collection of features from a known segment (i.e., the cycle over the force-plate), similar sequences within a trial can be found and used for the study. Furthermore, some movements are not so well defined with these external sensing tools, and rather a template movement selected by a human is the most reasonable way to identify the sequence of data describing a particular motion. 

seg1d is an open-source Python package for the automated segmentation and extraction of time series data using one or more reference sequences. The segmentation process allows users to apply various methods and parameters for the process through weighted reference features in a rolling correlation size-varying window of any scale below the length of the targeted data. Correlations can be averaged across the references, and a peak detection algorithm finds individual segments. Non-overlapping segments are identified, and a clustering algorithm groups the most similar subsequence movements within the target. The package was developed for movement sciences but can be useful to anyone interested in extracting correlated subsequences from a dataset. 

![Sample segments in a timeseries from a reference \label{fig:example}](api_basic-1.png)


# Acknowledgements

This work was supported in part by the Agency for Science, Technology and Research (A\*STAR), Nanyang Technological University (NTU) and the National Health Group (NHG) (RRG3: 2019/19002).

Todd Pataky is supported through the Kiban B Grant 17H02151 (Japan Society for the Promotion of Science).

# References
Features Class
==============

.. currentmodule:: seg1d.processing.Features


.. autosummary::
	:toctree: generated/

	match_len
    center
    shared
    gen_weights
    meaningful
===============
Getting Started
===============

In this section, a few simple examples are given to ensure the installation is working
and you are able to segment data.  These examples call a method provided that uses
default parameters for the algorithms within the segmentation process but allows
user input for the data related parameters such as scaling and step size for the 
matching process. 


Sample Data
===========

A basic example of checking the installation using the included sample data.

.. automodule:: seg1d.examples.ex_simple
	:members:
	

Sine Wave
=========

.. automodule:: seg1d.examples.ex_sine
	:members:
	

Gauss
=====

In this example a Gaussian pulse is used to show segmentation on the varying shape of different amplitude. 
As the center arc is given as reference, the multiple extending arcs are found as well. Through 
the output of the segments, the correlation values can be seen to decrease, although still
clustered within the group. 

.. automodule:: seg1d.examples.ex_gauss
	:members:==============
Code Reference
==============


Simple Interface
================

This method is shown in the getting started section of the documentation. It provides an 
easy-to-use interface for segmenting data, although it has limited functionality in 
fine-tuning parameters of the underlying algorithms. 

.. toctree::
    code_simple.rst


Class
=====

The Class reference should be used for any advanced usage or when needing 
convenience functions for adding data and tuning parameters of the segmentation process. 


.. toctree::
   :maxdepth: 2

   segmenter.rst
   segmenter_meth.rst   
   processing.rst


Methods
=======

The Methods reference explains the individual routines that are used in the segmentation process. 
Although they are accessible directly, it is likely easier to (and recommended) to use the Class. 


.. toctree::
   :maxdepth: 2

   routines.rst



.. automodule:: seg1d

.. automodule:: seg1d.segment

.. automodule:: seg1d.algorithm

.. automodule:: seg1d.optimized_funcs

.. automodule:: seg1d.examplesSegmenter Methods
=================

These methods are used to handle data and run the segmentation process. 

.. currentmodule:: seg1d.Segmenter

.. autosummary::
   :toctree: generated/

    set_target
    add_reference
    clear_reference
    segment


.. rubric:: See Also

:class:`seg1d.Segmenter`Simple Interface
================

.. currentmodule:: seg1d

.. autofunction:: segment_data
    :noindex:
Algorithm Methods
=================

.. currentmodule:: seg1d.algorithm

.. rubric:: General Functions

.. autosummary::
   :toctree: generated/

   rolling_corr
   combine_corr
   uniques
   get_peaks
   cluster
   resample


.. rubric:: Optimized Functions

.. currentmodule:: seg1d.optimized_funcs

.. autosummary::
   :toctree: generated/

   rcor
   vcor
Feature Processing
==================

If there is a series of references that have different lengths, or the weighting is unknown, 
this example can help guide users to process the data. 

.. automodule:: seg1d.examples.ex_feature_processing
	:members:Feature Inclusion
=================

There may be a case where an original dataset has multiple features, but only a subset 
of these features are wanted to be included in the segmentation process.  

.. automodule:: seg1d.examples.ex_segmenter_features
	:members:Segmenter Class
===============

.. currentmodule:: seg1d

.. autosummary::
	:toctree: generated/
	:nosignatures:

		Segmenter


.. autosummary::
	:toctree: generated/
		Segmenter.clusters
		Segmenter.corrs
		Segmenter.groups
		Segmenter.combined
		Segmenter.peaks
		Segmenter.t_masked
		Segmenter.t_segments====================
Community Guidelines
====================

Issues
------

Issues and feature requests should be submitted on `github <https://github.com/cadop/seg1d/issues>`_ . 


Contributing
------------

Please follow the "fork-and-pull" Git workflow. 
However, it is suggested to create an issue first to confirm the contributions align with the future of the module.

1. **Fork** the repo on GitHub
2. **Clone** the project to your own machine
3. **Commit** changes to your own branch
4. **Push** your work back up to your fork
5. Submit a **Pull request** so that we can review your changes

Documentation is on the main branch of the repository (rather than gh-pages) and should be built as-is. 
A redirect ``index.html`` file handles moving the github pages to the build directory and a ``.nojekyll`` file 
preserves folder types. 

Copyright and Licensing
-----------------------

Please refer to the full `LICENSE <https://github.com/cadop/seg1d/blob/master/LICENSE.txt>`_ text.====================================
seg1d: subsequence segmentation
====================================

Contents
========


.. toctree::
   :maxdepth: 2
   
   install.rst
   start.rst
   api.rst
   code.rst
   community.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
Basic Use
=========

.. automodule:: seg1d.examples.ex_segmenter_sine
	:members:Parameter Tuning
================

In some cases, the default values used in the segmentation process
do not result in the desired results. In this case, the various parameters
that are invovled can be manually set by the user. These parameters are 
all defined and available through the ``seg1d.Segmenter()`` class. 


.. automodule:: seg1d.examples.ex_sine_noise
	:members:============
Installation
============

Check dependencies
------------------
Currently tested on ``Python 3.8+``.

For the package: ``numpy>=1.15``, ``scipy>=1.0.0``, ``sklearn>=0.2``, ``numba>=0.40``

For building documentation: ``sphinx``, ``scipy`` html_theme, ``numpydoc``, ``matplotlib``

Get seg1d
-------------

``pip install seg1d``

Building Documentation
----------------------

Documentation is built using sphinx with the scipy html_theme in the `docs` folder. 

To get the requirements for building documentation you can run:

``pip install sphinx numpydoc matplotlib``

If you have not installed ``seg1d`` you can then run the following to get the remaining requirements

``pip install numpy scipy sklearn numba``

Finally, the following two lines will build the documentation. 

``make clean``

``make html``


Testing
--------

Docstrings examples should be compliant with ``doctest``. 
Navigate to the `docs` folder and run:

``make doctest``

A summary report will be generated and should not have any failed conditions. 
Note, this should be done after building the documentation with ``make html``.

End-users can test the installation with a similar method provided through the
examples. 

After installation, start a python interactive console and import:

``from seg1d.examples import test``

Then run the tests with:

``test.run()``

A series of tests will be performed to check for intended output. 
While ``matplotlib`` is not a required dependency, if it is not installed, some
of the tests will fail. This will be obvious as the test will show in the 
traceback that it failed due to ``No module named 'matplotlib'``. 

Due to variations in OS, many of the doctests use ``np.around`` to ensure matching. 
However, sorting of correlation values is dependent on a few platform variables, 
so tests that output multiple array values that are identical may show as failed. 
The documentation is always tested on Windows 10, Python 3.8 for ensuring compliance. 
============
API Examples
============

seg1d can be used through both the routines and through the ``Segmenter`` class,
which provides various methods for constructing a segmenter by setting parameters and data. 

This section contains the more advanced usage of the package by creating an instance
of the class. The usage is explained through different sample datasets, both of real
and generated data. Each example emphasizes a different aspect of the segmentation process. 


.. rubric:: API Examples 

.. toctree::
	api_basic.rst
	api_ecg.rst
	api_feat.rst
	api_tune.rst
	api_processing.rst


.. rubric:: See Also

:class:`seg1d.Segmenter`

ECG 
===

.. automodule:: seg1d.examples.ex_ecg
	:members:
seg1d.algorithm.combine\_corr
=============================

.. currentmodule:: seg1d.algorithm

.. autofunction:: combine_corrseg1d.processing.Features.gen\_weights
======================================

.. currentmodule:: seg1d.processing

.. automethod:: Features.gen_weightsseg1d.Segmenter.set\_target
===========================

.. currentmodule:: seg1d

.. automethod:: Segmenter.set_targetseg1d.Segmenter.clusters
========================

.. currentmodule:: seg1d

.. autoproperty:: Segmenter.clustersseg1d.Segmenter.t\_segments
===========================

.. currentmodule:: seg1d

.. autoproperty:: Segmenter.t_segmentsseg1d.algorithm.resample
========================

.. currentmodule:: seg1d.algorithm

.. autofunction:: resampleseg1d.processing.Features.center
================================

.. currentmodule:: seg1d.processing

.. automethod:: Features.centerseg1d.Segmenter.clear\_reference
================================

.. currentmodule:: seg1d

.. automethod:: Segmenter.clear_referenceseg1d.processing.Features.meaningful
====================================

.. currentmodule:: seg1d.processing

.. automethod:: Features.meaningfulseg1d.Segmenter.peaks
=====================

.. currentmodule:: seg1d

.. autoproperty:: Segmenter.peaksseg1d.Segmenter.combined
========================

.. currentmodule:: seg1d

.. autoproperty:: Segmenter.combinedseg1d.Segmenter.t\_masked
=========================

.. currentmodule:: seg1d

.. autoproperty:: Segmenter.t_maskedseg1d.Segmenter
===============

.. currentmodule:: seg1d

.. autoclass:: Segmenter

   
   .. automethod:: __init__

   
   .. rubric:: Methods

   .. autosummary::
   
      ~Segmenter.__init__
      ~Segmenter.add_reference
      ~Segmenter.clear_reference
      ~Segmenter.segment
      ~Segmenter.set_target
   
   

   
   
   .. rubric:: Attributes

   .. autosummary::
   
      ~Segmenter.clusters
      ~Segmenter.combined
      ~Segmenter.corrs
      ~Segmenter.groups
      ~Segmenter.peaks
      ~Segmenter.t_masked
      ~Segmenter.t_segments
   
   seg1d.Segmenter.segment
=======================

.. currentmodule:: seg1d

.. automethod:: Segmenter.segmentseg1d.Segmenter.add\_reference
==============================

.. currentmodule:: seg1d

.. automethod:: Segmenter.add_referenceseg1d.optimized\_funcs.vcor
===========================

.. currentmodule:: seg1d.optimized_funcs

.. autofunction:: vcorseg1d.processing.Features.shared
================================

.. currentmodule:: seg1d.processing

.. automethod:: Features.sharedseg1d.processing.Features.match\_len
====================================

.. currentmodule:: seg1d.processing

.. automethod:: Features.match_lenseg1d.optimized\_funcs.rcor
===========================

.. currentmodule:: seg1d.optimized_funcs

.. autofunction:: rcorseg1d.algorithm.get\_peaks
==========================

.. currentmodule:: seg1d.algorithm

.. autofunction:: get_peaksseg1d.algorithm.uniques
=======================

.. currentmodule:: seg1d.algorithm

.. autofunction:: uniquesseg1d.algorithm.cluster
=======================

.. currentmodule:: seg1d.algorithm

.. autofunction:: clusterseg1d.Segmenter.groups
======================

.. currentmodule:: seg1d

.. autoproperty:: Segmenter.groupsseg1d.Segmenter.corrs
=====================

.. currentmodule:: seg1d

.. autoproperty:: Segmenter.corrsseg1d.algorithm.rolling\_corr
=============================

.. currentmodule:: seg1d.algorithm

.. autofunction:: rolling_corr