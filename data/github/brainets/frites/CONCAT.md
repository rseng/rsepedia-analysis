.. image:: https://github.com/brainets/frites/actions/workflows/test_doc.yml/badge.svg
    :target: https://github.com/brainets/frites/actions/workflows/test_doc.yml

.. image:: https://github.com/brainets/frites/actions/workflows/flake.yml/badge.svg
    :target: https://github.com/brainets/frites/actions/workflows/flake.yml

.. image:: https://travis-ci.org/brainets/frites.svg?branch=master
    :target: https://travis-ci.org/brainets/frites

.. image:: https://codecov.io/gh/brainets/frites/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/brainets/frites

.. image:: https://badge.fury.io/py/frites.svg
    :target: https://badge.fury.io/py/frites

.. image:: https://pepy.tech/badge/frites
    :target: https://pepy.tech/project/frites

.. figure::  https://github.com/brainets/frites/blob/master/docs/source/_static/logo_desc.png
    :align:  center

======
Frites
======

Description
-----------

`Frites <https://brainets.github.io/frites/>`_ is a Python toolbox for assessing information-theorical measures on human and animal neurophysiological data (M/EEG, Intracranial). The aim of Frites is to extract task-related cognitive brain networks (i.e modulated by the task). The toolbox also includes directed and undirected connectivity metrics such as group-level statistics.

.. figure::  https://github.com/brainets/frites/blob/master/docs/source/_static/network_framework.png
    :align:  center

Documentation
-------------

Frites documentation is available online at https://brainets.github.io/frites/

Installation
------------

Run the following command into your terminal to get the latest stable version :

.. code-block:: shell

    pip install -U frites


You can also install the latest version of the software directly from Github :

.. code-block:: shell

    pip install git+https://github.com/brainets/frites.git


For developers, you can install it in develop mode with the following commands :

.. code-block:: shell

    git clone https://github.com/brainets/frites.git
    cd frites
    python setup.py develop
    # or : pip install -e .

Dependencies
++++++++++++

The main dependencies of Frites are :

* `Numpy <https://numpy.org/>`_
* `Scipy <https://www.scipy.org/>`_
* `MNE Python <https://mne.tools/stable/index.html>`_
* `Xarray <http://xarray.pydata.org/en/stable/>`_
* `Joblib <https://joblib.readthedocs.io/en/latest/>`_

In addition to the main dependencies, here's the list of additional packages that you might need :

* `Numba <http://numba.pydata.org/>`_ : speed up the computations of some functions
* `Dcor <https://dcor.readthedocs.io/en/latest/>`_ for fast implementation of distance correlation
* `Matplotlib <https://matplotlib.org/>`_, `Seaborn <https://seaborn.pydata.org/>`_ and `Networkx <https://networkx.github.io/>`_ for plotting the examples
* Some example are using `scikit learn <https://scikit-learn.org/stable/index.html>`_ estimators
Frites
======

.. image:: https://github.com/brainets/frites/actions/workflows/test_doc.yml/badge.svg
    :target: https://github.com/brainets/frites/actions/workflows/test_doc.yml

.. image:: https://github.com/brainets/frites/actions/workflows/flake.yml/badge.svg
    :target: https://github.com/brainets/frites/actions/workflows/flake.yml

.. image:: https://travis-ci.org/brainets/frites.svg?branch=master
    :target: https://travis-ci.org/brainets/frites

.. image:: https://codecov.io/gh/brainets/frites/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/brainets/frites

.. image:: https://badge.fury.io/py/frites.svg
    :target: https://badge.fury.io/py/frites

.. image:: https://pepy.tech/badge/frites
    :target: https://pepy.tech/project/frites



.. figure::  _static/logo_desc.png
    :align:  center


Description
+++++++++++

**Frites** is a Python toolbox for assessing information-based measures on human and animal neurophysiological data (M/EEG, Intracranial). The toolbox also includes directed and undirected connectivity metrics such as group-level statistics on measures of information (information-theory, machine-learning and measures of distance).


Highlights
++++++++++

.. panels::
    :container: container-lg pb-1
    :header: text-center
    :card: shadow

    Measures of information
    ^^^^^^^^^^^^^^^^^^^^^^^

    Extract cognitive brain networks by linking brain data to an external task
    variable by means of powerful, potentially multivariate, measures of
    information from fields like `information theory <https://brainets.github.io/frites/api/api_estimators.html#information-theoretic-estimators>`_, machine-learning,
    `measures of distances <https://brainets.github.io/frites/api/api_estimators.html#distance-estimators>`_
    or by defining your own `custom estimator <https://brainets.github.io/frites/api/generated/frites.estimator.CustomEstimator.html#frites.estimator.CustomEstimator>`_.

    +++

    .. link-button:: https://brainets.github.io/frites/api/api_estimators.html
        :text: List of estimators
        :classes: btn-outline-primary btn-block

    ---

    Group-level statistics
    ^^^^^^^^^^^^^^^^^^^^^^

    Combine measures of information with group-level statistics to extract
    robust effects across a population of subjects or sessions. Use either fully
    automatized non-parametric permutation-based `workflows <https://brainets.github.io/frites/api/generated/frites.workflow.WfMi.html>`_ or `manual ones <https://brainets.github.io/frites/api/generated/frites.workflow.WfStats.html>`_

    +++

    .. link-button:: https://brainets.github.io/frites/api/api_workflow.html
        :text: List of workflows
        :classes: btn-outline-primary btn-block

    ---

    Measures of connectivity
    ^^^^^^^^^^^^^^^^^^^^^^^^

    Estimate whole-brain pairwise undirected and directed connectivity.

    +++

    .. link-button:: https://brainets.github.io/frites/api/api_connectivity.html
        :text: Connectivity metrics
        :classes: btn-outline-primary btn-block stretched-link

    ---

    Link with external toolboxes
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    Frites supports inputs from standard libraries like `Numpy <https://numpy.org/>`_,
    `MNE Python <https://mne.tools/stable/index.html>`_ or more recent ones like
    labelled `Xarray <http://xarray.pydata.org/en/stable/>`_ objects.

    +++

    .. link-button:: https://brainets.github.io/frites/auto_examples/index.html#xarray
        :text: Xarray quick tour
        :classes: btn-outline-primary btn-block

.. toctree::
   :maxdepth: 2
   :hidden:

   Overview <overview/index>
   Installation <install>
   API Reference <api/index>
   Examples <auto_examples/index>
Installation
------------

.. contents::
   :local:
   :depth: 2

Dependencies
++++++++++++

The main dependencies of Frites are :

* `Numpy <https://numpy.org/>`_
* `Scipy <https://www.scipy.org/>`_
* `MNE <https://mne.tools/stable/index.html>`_
* `Xarray <http://xarray.pydata.org/en/stable/>`_
* `Joblib <https://joblib.readthedocs.io/en/latest/>`_

In addition to the main dependencies, here's the list of additional packages that you might need :

* `Numba <http://numba.pydata.org/>`_ : speed computations of some functions
* `Matplotlib <https://matplotlib.org/>`_, `Seaborn <https://seaborn.pydata.org/>`_ and `Networkx <https://networkx.github.io/>`_ for plotting the examples


Installation from pip
+++++++++++++++++++++

Frites can be installed (and/or updated) via pip with the following command :

.. code-block:: shell

    pip install -U frites

If you choose this installation method, you'll get stable Frites' release. If you want the latest features and patches, you can directly install from Github (see bellow)


Installation from Github
++++++++++++++++++++++++

Run the following line in your terminal to install the latest version of Frites hosted on Github :

.. code-block:: shell

    pip install git+https://github.com/brainets/frites.git


Developer installation
++++++++++++++++++++++

For developers, you can install frites in develop mode with the following commands :

.. code-block:: shell

    git clone https://github.com/brainets/frites.git
    cd frites
    python setup.py develop
    # or : pip install -e .{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autofunction:: {{ objname }}

.. _sphx_glr_backreferences_{{ fullname }}:

.. minigallery:: {{ fullname }}
    :add-heading:
{{ fullname | escape | underline}}


.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :special-members: __contains__,__getitem__,__iter__,__len__,__add__,__sub__,__mul__,__div__,__neg__,__hash__
   :members:

.. _sphx_glr_backreferences_{{ fullname }}:

.. minigallery:: {{ fullname }}
    :add-heading:
.. currentmodule:: frites

What's new
==========

v0.4.2
------

New Features
++++++++++++
* New function :func:`frites.simulations.sim_ground_truth` for simulating spatio-temporal ground-truths (:commit:`ba44a424`)
* New function :func:`frites.conn.conn_spec` for computing the single-trial spectral connectivity (:commit:`8151486`) - :author:`ViniciusLima94`
* New method :class:`frites.workflow.WfMi.confidence_interval` method to estimate the confidence interval (:commit:`ad0391987`, :commit:`8189622b`, :commit:`fc584756`, :commit:`fc584756`)
* New function :func:`frites.conn.conn_net` for computing the net connectivity (:commit:`c86b19f0`)

v0.4.1
------

New Features
++++++++++++
* New :class:`frites.estimator.CustomEstimator` for defining custom estimators (:commit:`e473c713`, :commit:`5584654c`)
* New function :func:`frites.conn.conn_fcd_corr` for computing the temporal correlation across networks (:commit:`2001f0c0`)
* New function :func:`frites.utils.acf` for computing the auto-correlation (:commit:`48ef0a03`)
* New function :func:`frites.conn.conn_ccf` for computing the cross-correlation (:commit:`43fceb00`)

Bug fixes
+++++++++
* Fix attribute conversion in connectivity functions (:commit:`b990c76`)

v0.4.0
------

New Features
++++++++++++
* New estimators (:class:`frites.estimator.CorrEstimator`, :class:`frites.estimator.DcorrEstimator`) for continuous / continuous relationships (:commit:`73ed8bbb`, :commit:`bc370a93`, :commit:`cf7a3456f`)
* :func:`frites.conn.conn_dfc` supports passing other estimators (:commit:`a864a7b05b`)
* :func:`frites.utils.time_to_sample` and :func:`frites.utils.get_closest_sample` conversion functions (:commit:`7c44478e`)
* :func:`frites.conn.conn_ravel_directed` reshaping function (:commit:`f9b9d272`)
* New :class:`frites.workflow.WfMi.copy` for internal workflow copy (:commit:`0c2228c7`, :commit:`860f3d45`)
* New :class:`frites.workflow.WfMiCombine` and example class for combining workflows (:commit:`62072ee52`)
* New :class:`frites.estimator.ResamplingEstimator` trial-resampling estimator (:commit:`13f6271e`)

Bug fixes
+++++++++
* Fix :class:`frites.workflow.WfMi.get_params` when returning a single output (:commit:`3bde82e6`)
* Improve attributes conversion for saving netcdf files (bool and dict) (:commit:`8e7dddb1`, :commit:`c6f7a4db`)
* Fix Numpy `np.float` and `np.int` warnings related (:commit:`896a198a`, :commit:`7f2a1caef`, :commit:`0d1a1223`)

v0.3.9
------

New Features
++++++++++++
* :func:`frites.conn.conn_dfc` supports multivariate data + improve computing efficiency (:commit:`1aed842`, :commit:`c4ac490`)
* Reshaping connectivity arrays support elements on the diagonal + internal drop of duplicated elements (:commit:`daac241f`)
* :func:`frites.conn.conn_dfc` supports better channel aggregation (:commit:`a66faa77`)

Internal Changes
++++++++++++++++
* Connectivity metric now use the :class:`frites.dataset.SubjectEphy` for internal conversion of the input data
* :class:`frites.workflow.WfMi.get_params` returns single-subject MI and permutations with dimension name 'subject' (instead of subjects) (:commit:`85884f3a`)
* All connectivity metrics now use :func:`frites.conn.conn_io` to convert inputs into a similar format
* Improve CI

Bug fixes
+++++++++
* Fix :class:`frites.dataset.SubjectEphy` when the data contains a single time point (:commit:`a33d4437`)
* Fix attributes of :func:`frites.conn.conn_covgc` (:commit:`c120626`)
* Fix :class:`frites.dataset.DatasetEphy` representation without data copy + html representation (:commit:`b3ae7b8ea`, :issue:`16`)
* Fix passing `tail` input to the :class:`frites.workflow.WfMi` (:commit:`6df86d1e`)


v0.3.8
------

New Features
++++++++++++
* new :class:`frites.io.Attributes` class for managing and printing datasets' and workflow's attributes (:commit:`be046b1`)
* new :class:`frites.dataset.SubjectEphy` single-subject container (:commit:`ac22cf4`)
* new estimators of mutual-information, :class:`frites.estimator.GCMIEstimator` (:commit:`901b3cbf`, :commit:`65d1e08`, :commit:`0015bf58`, :commit:`beed6a09`), :class:`frites.estimator.BinMIEstimator` (:commit:`beed6a09`)
* new kernel smoothing function :func:`frites.utils.kernel_smoothing`

Internal Changes
++++++++++++++++
* Removed files (:commit:`cdff9b4`, :commit:`9e96f8e`, :commit:`14961aa0`)
* :class:`frites.dataset.DatasetEphy` don't perform internal data copy when getting the data in a specific ROI (:commit:`2da73ef`)
* Compatibility of MI estimators with workflows (:commit:`7dc76ee9`, :commit:`e7a9c23f`)
* Improve the way to manage pairs of brain regions (:commit:`8b955a16`, :commit:`bfdf2dba`, :commit:`57c1e4ba`, :commit:`b1ff8c3d`)

Breaking changes
++++++++++++++++
* :class:`frites.dataset.SubjectEphy` and :class:`frites.dataset.DatasetEphy` to specify whether channels should be aggregated (default `agg_ch=True`) or not (`agg=False`) when computing MI. The `agg_ch` replace `sub_roi` (:commit:`18d4e24`)
* The workflow `WfComod` has been renamed :class:`frites.workflow.WfConnComod` (:commit:`b7b58248`)

Bug fixes
+++++++++
* Bug fixing according to the new version of :class:`frites.dataset.DatasetEphy` (:commit:`1a15e05`, :commit:`7b83a3d`, :commit:`abd1b281`, :commit:`70bfefb`, :commit:`5879950`, :commit:`66acdf2`, :commit:`4309be9c5`, :commit:`6dc2fbf8`)

v0.3.6
------

New Features
++++++++++++
* :class:`frites.dataset.DatasetEphy` support multi-level anatomical informations (:commit:`3a9ce540`)

Internal Changes
++++++++++++++++
* Connectivity functions have a better support of Xarray inputs (:commit:`60cc16`, :commit:`b72d1519`, :commit:`3d24c98`, :commit:`65bf08`)
* Replace every string comparison 'is' with '==' to test the content (:commit:`1337aa6e`)

v0.3.5
------

New Features
++++++++++++
* New function for reshaping undirected connectivity arrays (like DFC) :func:`frites.conn.conn_reshape_undirected` (:commit:`ffcae34`, :commit:`56515fe`)
* New function :func:`frites.utils.savgol_filter` that works on DataArray (:commit:`3e0e256`)
* New function for reshaping directed connectivity arrays (like COVGC) :func:`frites.conn.conn_reshape_directed` (:commit:`8c2bb63`)
* New method :class:`frites.workflow.WfMi.get_params` in order to get the internal arrays formatted as DataArray (:commit:`03dd2f3`)
* Integration of MNE's progress bar (:commit:`74dc66`, :commit:`8ec636d`, :commit:`2bb7e75`)
* Possibility to cache computations in the parallel function

Internal Changes
++++++++++++++++
* DataArray now contain a name such as a type to make it clear what is it (:commit:`2d6a61`)
* sigma parameter when performing the t-test can be changed though the CONFIG file (:commit:`5d7ba9f`)

Bug fixes
+++++++++
* Fix ttested attribute when saving (:commit:`3d029be`)
* Fix computing the sigma across all ROI since it uses the maximum over every axes (:commit:`ffaca1e`)
* Fix high RAM consumption when computing the `pop_mean_surr` (:commit:`fe6e4d1`)


.. raw:: html

    <hr>

v0.3.4
------

Bug fixes
+++++++++
* Fix :class:`frites.workflow.WfMi.conjunction_analysis` for seeg data (:commit:`6ca64c5`)

Breaking changes
++++++++++++++++
* Xarray is now the default output type as it supports the definition of multi-dimensional containers with label to each coordinates (:commit:`8174e4`) + illustrating examples of how to use xarray (:commit:`f5d28e`)
* Every connectivity measures have been moved to `frites.conn` (:commit:`a33d536`)
* Deep support for correcting for multiple-comparisons using multi-dimensional data (:commit:`840327`, :commit:`e0e1a5b`) + support for :class:`frites.workflow.WfStatsEphy` (:commit:`563910d`) + support for :class:`frites.dataset.DatasetEphy` (:commit:`10a3697`) + support for :class:`frites.workflow.WfMi` (:commit:`9c6165`)

New Features
++++++++++++
* New class :class:`frites.simulations.StimSpecAR` for generating Auto-Regressive Models (:commit:`4688539`, :commit:`5fd1199`)
* Conditional Covariance based Granger Causality + example (:commit:`a148310`)

Internal Changes
++++++++++++++++
* Improve testings (:commit:`6644cd`)

.. raw:: html

    <hr>

v0.3.3
------

Internal Changes
++++++++++++++++
* :class:`frites.workflow.WfFit` and :class:`frites.workflow.WfConn` are now using :class:`frites.dataset.DatasetEphy.get_connectivity_pairs` (:commit:`18e1`)
* Improve warning messages + assert error for negative determinant for :class:`frites.core.covgc`

New Features
++++++++++++
* New method :class:`frites.dataset.DatasetEphy.get_connectivity_pairs` in order to get possible connectivity pairs (:commit:`4e77`)
* New function :func:`frites.utils.define_windows` and :func:`frites.utils.plot_windows` in order to generate and plot slicing windows + tests + example (:commit:`6fcf`)
* New function for computing the DFC :class:`frites.core.dfc_gc` (:commit:`8f6f`)
* When using DataArray with :class:`frites.core.dfc_gc` and :class:`frites.core.covgc`, temporal attributes are added (:commit:`5df6`, :commit:`6266`)
* New function for computing the covgc :class:`frites.core.covgc` (:commit:`ea26`)
* Step parameter for :class:`frites.core.covgc` (:commit:`874af`)
* :class:`frites.core.covgc` can no be computed using Gaussian-Copula (:commit:`aea6a8b`)
* Add :class:`frites.workflow.WfMi.conjunction_analysis` for performing conjunction analysis + example (:commit:`7dd0bea`)

Bug fixes
+++++++++
* Fix when data contains a single time point (:commit:`c815e40`)
* Fix mi model and mixture (1d only) (:commit:`2b119d4`)

.. raw:: html

    <hr>

v0.3.2
------

Breaking changes
++++++++++++++++
* Avoid duplicates dataset construction when using MNE / xarray (:commit:`fef9`, :commit:`9dde`)
* :class:`frites.dataset.DatasetEphy` supports None for the y input (:commit:`a53a`)

Internal Changes
++++++++++++++++
* Dtypes of `y` and `z` inputs are systematically check in :class:`frites.dataset.DatasetEphy` in order to define which MI can then be computed (:commit:`7cc1`)

New Features
++++++++++++
* :class:`frites.dataset.DatasetEphy` supports Xarray inputs + selection though coordinates (:commit:`7418`)
* New workflow for computing pairwise connectivity :class:`frites.workflow.WfConn` (:commit:`65ae`)

Documentation
+++++++++++++
* Adding new examples for creating datasets (:commit:`4df9`)


.. raw:: html

    <hr>

v0.3.1
------

Breaking changes
++++++++++++++++
* change :class:`frites.workflow.WfFit` input `directed` for `net` (:commit:`3ec0`, :issue:`1`)
* The GCRN is automatically defined (per subject when RFX / across subjects when FFX) (:commit:`b8a9`)
* Remove the `level` input parameter + only `mcp` is used + only maxstat when performing cluster based (:commit:`d8fb`)

New Features
++++++++++++
* :class:`frites.dataset.DatasetEphy` support spatio-temporal slicing (:commit:`60d1`), resampling (:commit:`3785`) and savitzki-golay filter (:commit:`9707`)
* Support setting `random_state` in :class:`frites.workflow.WfMi` and :class:`frites.workflow.WfFit` (:commit:`0688`)
* DataArray outputs contains attributes that reflect the configuration of the workflow (:commit:`18181`)

Bug fixes
+++++++++
* Fix multi-sites concatenation (:commit:`3bcc`) in :class:`frites.dataset.DatasetEphy`
* Fix p-values to zeros in :class:`frites.workflow.WfFit` (:commit:`2062`)
* Fix FIT outputs for 3D arrays and DataArray (:commit:`18181`)

Internal Changes
++++++++++++++++
* Remap multiple conditions when integers (:commit:`6092`) in :class:`frites.dataset.DatasetEphy`
* Workflows now have an internal configuration (:commit:`18181`)

Documentation
+++++++++++++
* Reformat examples gallery (:commit:`4d4c`)
References
----------

.. bibliography:: ../refs.bib
    :style: plain
Methods
-------

In this section we are going to cover the basis of the methods used in Frites. In particular two main topics :

1. The use of information theoretical measures for brain data
2. The implemented statistical approaches

Note that this methodological overview has been written for scientists who are not necessarily familiar with information theory and statistics. Therefore, it is not going to contain the mathematics behind the implemented methods, but instead, the references to the papers in question.

Data analysis within the information theoretical framework
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The first question that could be asked is *Why using the information theory to analyze brain data?*

There are several advantages of using information theoretic measures. Here's a non exhaustive list of the main advantages :

1. **No assumption about the distribution of the data :** standard IT measures (based on binning) are model-free in the sense that there's no assumptions about how the data are distributed
2. **Information theory as a unified framework of analysis :** most of the statistical measures we used to analyze data (e.g correlation, machine-learning etc.) have a equivalent using IT metrics. This means that in theory it is possible to analyze the data from A-Z within the IT framework
3. **Support for any data type :** by construction, IT measures support both continuous and discrete variables such as univariate and multivariate.

On the other hand, IT measures also present some disadvantages :

1. **Standard approaches usually require large data :** indeed, methods based on binning need to evaluate the probability distribution (e.g the probability that your data take values between [-1, 1]). And the larger the data, the better the estimation. In addition, binning based IT measures suffer from the curse of dimensionality.
2. **Strictly positive measures :** IT measures returned the amount of information that is shared between two variables (in bits). And this quantity is superior or equal to zero. As an example, two variables might either be correlated or anti-correlated however, in both cases the amount of information is going to be strictly positive.
3. **Computationally expansive :** depending on the implementation, the evaluation of the mutual information can be relatively slow

.. figure::  ../_static/gcmi_corr.png
    :align:  center

    Correlation (orange) vs. strictly positive Mutual information (purple). Extracted from **Ince et al. 2017** :cite:`ince2017statistical`

Gaussian Copula Mutual Information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Frites we use a special type of mutual information that has been carried recently to neuroscience (see Ince et al 2017 :cite:`ince2017statistical`) : the Gaussian Copula Mutual Information (GCMI).

The GCMI provide an estimation of the mutual information **without the necessity of binning the data**. In addition to this point, the GCMI is really fast to compute which makes this method a good candidate when it comes to performing numerous permutations for the statistics. However, the GCMI make the assumption that there's a **monotonic relation between the two variables** (i.e if one variable is increasing, the other should also increase or decrease but not both).


The different types of mutual information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this section we are going to make the difference between two types of variables :

* **Discrete variables :** discrete variables are composed with **integers** where each number can for example refer to an experimental condition (e.g 0={all the trials that belong to the "go" stimulus} and 1={all the trials that belong to the "nogo" stimulus})
* **Continuous variables :** continuous variables are made of **floating points**. Therefore, the brain data are considered as a continuous variable.

With this in mind, the figure bellow summarizes the different types of mutual information :

.. figure::  ../_static/cc_cd_ccd.png
    :align:  center

    **Top left** : continuous / discrete case (``I(C; D); mi_type='cd'``) for contrasting two experimental conditions (e.g two visual stimulus); **Top right** continuous / continuous case (``I(C; C); mi_type='cc'``) for correlating the brain data with, for example, a model of the behavior (e.g Prediction Error / Learning); **Bottom** continuous / continuous | discrete case (``I(C; C | D); mi_type='ccd'``) for correlating the brain data with a continuous while removing the influence of a discrete variable (e.g experimental conditions)


Equivalence with other statistical measures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Find bellow the table showing the equivalence between information theoretic quantities and other statistical measures. Extracted from **Ince et al. 2017** :cite:`ince2017statistical`

+-------------------------------------------------------+----------------------------------------------------------------------+
| Information theoretic quantity                        | Other statistical approaches                                         |
+=======================================================+======================================================================+
| MI (discrete; discrete)                               | - Chi-square test of independence                                    |
|                                                       | - Fishers exact test                                                 |
+-------------------------------------------------------+----------------------------------------------------------------------+
| MI (univariate continuous; discrete)                  | - 2 classes: t test, KS test, Mann–Whitney U test                    |
|                                                       | - ANOVA                                                              |
+-------------------------------------------------------+----------------------------------------------------------------------+
| MI (multivariate continuous; discrete)                | - 2 classes: Hoteling T2 test                                        |
|                                                       | - Decoding (cross-validated classifier)                              |
+-------------------------------------------------------+----------------------------------------------------------------------+
| MI (univariate continuous; univariate continuous)     | - Pearson correlation                                                |
|                                                       | - Spearman rank correlation                                          |
|                                                       | - Kendall rank correlation                                           |
+-------------------------------------------------------+----------------------------------------------------------------------+
| MI (multivariate continuous; univariate continuous)   | - Generalized Linear Model framework                                 |
|                                                       | - Decoding (cross-validated regression)                              |
+-------------------------------------------------------+----------------------------------------------------------------------+
| MI (multivariate continuous; multivariate continuous) | - Canonical correlation analysis                                     |
|                                                       | - Distance correlation                                               |
+-------------------------------------------------------+----------------------------------------------------------------------+
| Conditional mutual information                        | - Partial correlation (continuous variables and linear effects only) |
+-------------------------------------------------------+----------------------------------------------------------------------+
| Directed information (transfer entropy)               | - Granger Causality                                                  |
+-------------------------------------------------------+----------------------------------------------------------------------+
| Directed feature information                          | - Dynamic Causal Modeling                                            |
|                                                       | - Psychophysiological interactions                                   |
+-------------------------------------------------------+----------------------------------------------------------------------+
| Interaction information                               | - Representational similarity analysis (redundancy only)             |
|                                                       | - Cross-classification decoding (redundancy only)                    |
|                                                       | - Mediation analysis                                                 |
+-------------------------------------------------------+----------------------------------------------------------------------+


References
~~~~~~~~~~

* Ince et al. 2017 :cite:`ince2017statistical`
* Timme and Lapish 2018 :cite:`timme_tutorial_2018`


Statistical analyses
++++++++++++++++++++

In addition to the evaluation of the amount of information shared between the data and a feature (stimulus / behavior), Frites also contains a statistical pipeline to evaluate whether an effect can be considered as significant.

Subject and group-level statistical inferences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From a statistical standpoint, Frites gives the possibility to draw inferences either at the **subject-level** either at the **group-level** :

* **Subject-level :** the measures of informations and the p-values are computed per-subject. At the end, you can get the number (or proportion) of subjects that have the effect. This can be particularly useful when working with sEEG data where the number of subjects per brain regions might be different.
* **Group-level :** draw conclusions about a population of subjects. This group-level section is subdivided into two group-level strategies :

    * **Fixed Effect (FFX) :** the mutual information is estimated **across subjects**. By concatenating the subjects, you make the hypothesis that the effect is quite reproducible and stable across your population of subjects. One advantage of the FFX is that it usually requires a lower number of subjects and provides a good sensibility as soon as the hypothesis of stable effect across subjects is verified. On the other hand, two disadvantages are that, first you don't take into account the inter-subject variability and second, the inferences you are allow to make only concern **your population** of subjects. It can't generalize to new subjects
    * **Random Effect (RFX) :** the mutual information is estimated **per subject** and then a model of how the effect is distributed across the subjects is made. To be more precise, we make the assumption that the effect is normally distributed across the subjects. Advantages and disadvantages are basically the opposite of the FFX. The RFX takes into consideration the inter-subject variability which make it more suitable to detect the effect if this effect is slightly different from subject to subject. Building a model based is like considering that our population is in fact a sub-group of a broader population which, in other terms, means that if new subjects are included in the analysis, the estimations of mutual information should still be included in the model. However, the main disadvantage of the RFX is that it requires more subjects to build a reliable model.

.. figure::  ../_static/group_level.png
    :align:  center

    Subject and group-level statistics on information-based measures


Correction for multiple comparisons
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Internally, we used permutations in order to estimate the distribution of mutual-information values that could be obtained by chance and also, to correct the p-values for multiple comparisons. There's three ways to correct the p-values, all using `MNE Python <https://mne.tools/stable/index.html>`_ : 1) cluster-based statistics (default), 2) maximum statistics and 3) FDR. About the cluster-based, since some electrophysiological recordings don't have a spatial contiguity (like sEEG), the clusters are detected across time-points or across time-frequency maps. Note that you can also use the Threshold Free Cluster Enhancement (TFCE) still through MNE.

References
~~~~~~~~~~

* Friston et al. 1999 :cite:`friston1999many`, 1996 :cite:`friston_detecting_1996`
* Wilcox and Rousselet 2017, :cite:`wilcox_guide_2017`
* Nicohls and Holmes 2001 :cite:`nichols_nonparametric_2002`
* Maris and Oostenveld 2007 :cite:`maris2007nonparametric`
* Smith and Nichols, 2009 :cite:`smith2009threshold`
* Cao and Zhang 2014 :cite:`cao_multiple_2014`
* Combrisson et al. 2015 :cite:`combrisson_exceeding_2015`
* Giordano et al. 2017 :cite:`giordano2017contributions`
Goals of Frites
---------------

Frites is an open-source Python software started in 2019 to **analyze neurophysiological data using information theoretical (IT) measures**. Overall, the goals of Frites could be summarized as follow :

1. **Extract task-related cognitive brain networks :** i.e brain regions and connectivity between brain regions that are modulated according to the task
2. **Group-level statistical inferences :** being able to draw conclusions on a population of subjects while taking into account the the richness of data at the single subject level


Task-related and feature-specificity
++++++++++++++++++++++++++++++++++++

The stimulus of the task (e.g visual stimulus like color circles) is first going to be processed by the brain, occurring some changes in the brain activity that are related to the stimulus. This is what we called **stimulus-specific** brain activity. Once processed, this stimulus is going to introduce some changes in the behavior of the subject (e.g decision making). Following the same reasoning as above we called **behavior-specific** the modulations in the activity of the brain that can explain a change in the behavior of the subject.

.. figure::  ../_static/feat_spec.png
    :align:  center

Now, with Frites we try to isolate the brain activity that is either related to the stimulus or to the behavior. In a more general way, this brain activity is called **task-related** or **feature-specific** where *feature* can either refer to the *stimulus* or to the *behavior*.

IT framework to extract task-related cognitive brain networks
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

A brain network is composed of two things :

1. **Nodes :** the nodes of a network are the brain regions (e.g Insula / vmPFC / dlPFC form a three nodes network)
2. **Edges :** the connectivity links between the nodes (e.g comodulations or information flow)

Frites contains a collection of statistical methods to test whether :

1. The **local neural activity** (i.e the activity of a particular node of the network) is modulated by the task
2. The **connectivity** (which can be undirected or directed) between brain regions is also modulated by the task

.. figure::  ../_static/network_framework.png
    :align:  center

    Framework using information theoretical measures on neurophysiological data to extract task-related local neural activity and connectivityCommunity
=========

Events
------

Brainhack Marseille 2020
++++++++++++++++++++++++

Frites was one of the proposed projects during the `Brainhack Marseille 2020 <https://brainhack-marseille.github.io/#project1>`_. During the event we worked on porting the core functions from a CPU implementation to a GPU implementation.

.. figure::  ../_static/bhk_mrs.png
    :align:  center

    Credit for the logo : `Julia Sprenger <https://github.com/JuliaSprenger>`_

Software presentations
----------------------

* Award for the best poster during the Human Brain Project summit 2021 (`link for downloading the poster <https://amubox.univ-amu.fr/s/xf2AjeDXMCiWZEz>`_)
* `Live MEEG 2020 <https://www.crowdcast.io/e/live-meeg-2020/7>`_
* Neuromatch 3.0


Authors and contributors
------------------------

Thanks goes to these wonderful people (`emoji key`_):

.. _emoji key: https://allcontributors.org/docs/en/emoji-key

.. raw:: html

    <!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
    <!-- prettier-ignore-start -->
    <!-- markdownlint-disable -->
    <table>
      <tr>
        <td align="center"><a href="https://github.com/EtienneCmb"><img src="https://avatars3.githubusercontent.com/u/15892073?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Etienne Combrisson</b></sub></a><br /><a href="https://github.com/brainets/frites/commits?author=EtienneCmb" title="Code">💻</a> <a href="#design-EtienneCmb" title="Design">🎨</a> <a href="#example-EtienneCmb" title="Examples">💡</a> <a href="#maintenance-EtienneCmb" title="Maintenance">🚧</a> <a href="#mentoring-EtienneCmb" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-EtienneCmb" title="Project Management">📆</a></td>
        <td align="center"><a href="http://andrea-brovelli.net/"><img src="https://avatars0.githubusercontent.com/u/19585963?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Andrea Brovelli</b></sub></a><br /><a href="https://github.com/brainets/frites/commits?author=brovelli" title="Code">💻</a> <a href="#ideas-brovelli" title="Ideas, Planning, & Feedback">🤔</a> <a href="#mentoring-brovelli" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-brovelli" title="Project Management">📆</a></td>
        <td align="center"><a href="https://github.com/StanSStanman"><img src="https://avatars1.githubusercontent.com/u/26648765?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Ruggero Basanisi</b></sub></a><br /><a href="https://github.com/brainets/frites/commits?author=StanSStanman" title="Code">💻</a> <a href="#design-StanSStanman" title="Design">🎨</a></td>
        <td align="center"><a href="https://github.com/ViniciusLima94"><img src="https://avatars3.githubusercontent.com/u/17538901?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Vinícius Lima</b></sub></a><br /><a href="https://github.com/brainets/frites/commits?author=ViniciusLima94" title="Code">💻</a></td>
        <td align="center"><a href="https://github.com/tprzybyl"><img src="https://avatars1.githubusercontent.com/u/58084045?v=4?s=100" width="100px;" alt=""/><br /><sub><b>tprzybyl</b></sub></a><br /><a href="https://github.com/brainets/frites/commits?author=tprzybyl" title="Code">💻</a></td>
        <td align="center"><a href="https://github.com/micheleallegra"><img src="https://avatars0.githubusercontent.com/u/23451833?v=4?s=100" width="100px;" alt=""/><br /><sub><b>micheleallegra</b></sub></a><br /><a href="https://github.com/brainets/frites/commits?author=micheleallegra" title="Code">💻</a> <a href="#ideas-micheleallegra" title="Ideas, Planning, & Feedback">🤔</a></td>
        <td align="center"><a href="http://www.robinince.net/"><img src="https://avatars0.githubusercontent.com/u/63155?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Robin Ince</b></sub></a><br /><a href="https://github.com/brainets/frites/commits?author=robince" title="Code">💻</a> <a href="#ideas-robince" title="Ideas, Planning, & Feedback">🤔</a></td>
      </tr>
      <tr>
        <td align="center"><a href="https://github.com/samuelgarcia"><img src="https://avatars1.githubusercontent.com/u/815627?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Garcia Samuel</b></sub></a><br /><a href="#ideas-samuelgarcia" title="Ideas, Planning, & Feedback">🤔</a></td>
        <td align="center"><a href="https://github.com/brungio"><img src="https://avatars0.githubusercontent.com/u/33055790?v=4?s=100" width="100px;" alt=""/><br /><sub><b>brungio</b></sub></a><br /><a href="https://github.com/brainets/frites/commits?author=brungio" title="Code">💻</a> <a href="#ideas-brungio" title="Ideas, Planning, & Feedback">🤔</a> <a href="#mentoring-brungio" title="Mentoring">🧑‍🏫</a> <a href="#projectManagement-brungio" title="Project Management">📆</a></td>
      </tr>
    </table>

    <!-- markdownlint-restore -->
    <!-- prettier-ignore-end -->

    <!-- ALL-CONTRIBUTORS-LIST:END -->
Get help
========

If you have questions or if you want to start discussing about new features and ideas, use the `Discussion <https://github.com/brainets/frites/discussions>`_ tab directly from the Github repository of Frites.Cite Frites
===========

If you use Frites in your data analysis, please use the following citation in your publications :

.. code-block:: latex

    @article{combrisson2021group,
      title={Group-level inference of information-based measures for the analyses of cognitive brain networks from neurophysiological data},
      author={Combrisson, Etienne and Allegra, Michele and Basanisi, Ruggero and Ince, Robin AA and Giordano, Bruno and Bastin, Julien and Brovelli, Andrea},
      journal={bioRxiv},
      year={2021},
      publisher={Cold Spring Harbor Laboratory}
    }Start analyzing your data with Frites
-------------------------------------

In this section we are going to cover the basics knowledges if you want to start analyzing your data with Frites.

.. note::

    The gallery of examples contains a `tutorial <https://brainets.github.io/frites/auto_examples/index.html#tutorials>`_ section that explain and illustrate the main usages and features of Frites, in a notebook fashion.

Package organization
++++++++++++++++++++

As other Python packages, Frites contains subject-specific submodules. Here's a short description of the main submodules for new users :

* `frites.dataset <https://brainets.github.io/frites/api.html#module-frites.dataset>`_ : container for the electrophysiological data coming from multiple subjects (see also `those examples <https://brainets.github.io/frites/auto_examples/index.html#multi-subjects-dataset>`_ that explain how to define a container depending on your data type)
* `frites.workflow <https://brainets.github.io/frites/api.html#module-frites.workflow>`_ : the workflows perform a series of analyzes (usually a first analysis to extract information using IT measures follow by a second statistical analysis)
* `frites.conn <https://brainets.github.io/frites/api.html#module-frites.conn>`_ : directed and undirected connectivity metrics that can either be computed across trials or at the single trial level

In addition to the main modules above, the  `gallery of examples <https://brainets.github.io/frites/auto_examples/index.html>`_ illustrate the main functionalities of Frites using simulated data. Those functions can be imported from `frites.simulations <https://brainets.github.io/frites/api.html#module-frites.simulations>`_ and can be used to simulate local neural activity modulated by the task such as stimulus-specific brain networks using `autoregressive models <https://brainets.github.io/frites/api.html#autoregressive-model>`_.

Finally, for the developers here the `frites.core <https://brainets.github.io/frites/api.html#module-frites.core>`_ module include the very low level functions to estimate the mutual information (vector and tensor based implementations).


Main Frites' workflows
++++++++++++++++++++++

Frites contains two centrals workflows :

1. `WfMi <https://brainets.github.io/frites/api/generated/frites.workflow.WfMi.html#frites.workflow.WfMi>`_ : the main workflow of mutual-information that is used to **extract feature-specific brain networks**
2. `WfStats <https://brainets.github.io/frites/api/generated/frites.workflow.WfStats.html#frites.workflow.WfStats>`_ : the workflow of statistics we used to perform group-level inferences

Actually, under the hood, the **WfMi** is also using the workflow of statistics so that it return the estimation of mutual-information between a feature and the data but also the corresponding p-values.

.. warning::

    This workflow of mutual-information can be used for many things, including :

        * Provide a meaningful measure of effect-size between the data and a feature
        * To find some significant differences between multiple experimental conditions (``mi_type='cd``)
        * To correlate the brain data with a behavioral model (``mi_type='cc``) and possibly, independently of some experimental conditions (``mi_type='ccd``)
        * Get the corrected p-values to see whether you have a significant amount of information between your data and a feature, either at the group-level or at the subject level
        * Working on sparse data like the sEEG / iEEG where the number of subjects per brain region varies

    In a nutshell, the **WfMi** is the core workflow to find task-related activity. In addition, this workflow can either be used to quantify the amount of information between local neural activity and a feature **(i.e for each brain region)** but also on **connectivity links (i.e on pairs of brain regions)** to see whether there's some differences in term of connectivity between experimental conditions or ultimately, if the connectivity correlates with a behavioral model.

    To see this workflow in action, checkout `those examples <https://brainets.github.io/frites/auto_examples/index.html#mutual-information>`_


Deep integration with Xarray
++++++++++++++++++++++++++++

`Xarray <http://xarray.pydata.org/en/stable/>`_ is a recent python package to handle of multi-dimensional arrays using labels. For those who are familiar with `Pandas <https://pandas.pydata.org/>`_, you can see Xarray as a generalization of Pandas for multi-dimensional arrays.

Xarray provide a container for the data called `DataArray <http://xarray.pydata.org/en/stable/generated/xarray.DataArray.html#xarray.DataArray>`_. This structure comes with two important inputs : 1) `dims` which describe the name of each dimension of the array and 2) `coords` the values taken by each dimension. For example you can define a DataArray with the dimensions ``('roi', 'times')`` where ``roi = ['insula', 'vmPFC', 'dlPFC']`` and ``times = np.linspace(-3, 3, 1000)``. After that, the manipulation of the data happen using the values of the coordinates. Bellow, a minimal slicing example :

.. code-block:: python

    """
    `da` is a xarray.DataArray. With the code line below, we select the data
    coming from the two brain regions Insula and vmPFC. Then we also select
    every time points comprised between [-1.5, 1.5] seconds
    """
    da.sel(roi=['insula', 'vmPFC'], times=slice(-1.5, 1.5))


The example above only show how to slice the data but actually Xarray contains most of the operations using the same label-based syntax.

.. note::

    Frites make an extensive use of Xarray as most of the outputs returned are DataArrays. Since it's a relatively recent package, we wrote `two mini tutorials <https://brainets.github.io/frites/auto_examples/index.html#xarray>`_ to start working with it.
Overview
========

This section introduces the basis for starting analyzing your data with Frites.

.. toctree::
    :maxdepth: 2

    ovw_goals
    ovw_methods
    ovw_frites
    ovw_help
    ovw_cite
    ovw_community
    ovw_refs
    ovw_whatsnewLow-level core functions
------------------------

:py:mod:`frites.core`:

.. currentmodule:: frites.core

.. automodule:: frites.core
   :no-members:
   :no-inherited-members:


Gaussian-Copula (1d)
++++++++++++++++++++

Gaussian-Copula mutual-information supporting univariate / multivariate 2D inputs

.. autosummary::
   :toctree: generated/

   copnorm_1d
   copnorm_cat_1d
   ent_1d_g
   mi_1d_gg
   mi_model_1d_gd
   mi_mixture_1d_gd
   cmi_1d_ggg
   gcmi_1d_cc
   gcmi_model_1d_cd
   gcmi_mixture_1d_cd
   gccmi_1d_ccc
   gccmi_1d_ccd

Gaussian-copula (Nd)
++++++++++++++++++++

Gaussian-Copula mutual-information supporting univariate / multivariate multi-dimensional inputs

.. autosummary::
   :toctree: generated/

   copnorm_nd
   copnorm_cat_nd
   mi_nd_gg
   mi_model_nd_gd
   cmi_nd_ggg
   gcmi_nd_cc
   gcmi_model_nd_cd
   gccmi_nd_ccnd
   gccmi_model_nd_cdnd
   gccmi_nd_cccDataset
-------

:py:mod:`frites.dataset`:

.. currentmodule:: frites.dataset

.. automodule:: frites.dataset
   :no-members:
   :no-inherited-members:

.. autosummary::
   :toctree: generated/

   SubjectEphy
   DatasetEphy
Connectivity
------------

:py:mod:`frites.conn`:

.. currentmodule:: frites.conn

.. automodule:: frites.conn
   :no-members:
   :no-inherited-members:

Connectivity metrics
++++++++++++++++++++

.. autosummary::
   :toctree: generated/

   conn_dfc
   conn_covgc
   conn_ccf
   conn_spec
   conn_transfer_entropy

Utility functions
+++++++++++++++++

Reshaping connectivity outputs
==============================

.. autosummary::
   :toctree: generated/

   conn_reshape_undirected
   conn_reshape_directed
   conn_ravel_directed
   conn_net
   conn_get_pairs

Metrics to apply on FC
======================

.. autosummary::
   :toctree: generated/

   conn_fcd_corr

Define sliding windows
======================

.. autosummary::
   :toctree: generated/

   define_windows
   plot_windowsConfiguration
-------------

:py:mod:`frites.config`:

.. currentmodule:: frites

.. autosummary::
   :toctree: generated/

   get_config
   set_configUtility functions
-----------------

:py:mod:`frites.utils`:

.. currentmodule:: frites.utils

.. automodule:: frites.utils
   :no-members:
   :no-inherited-members:

Time-series processing
++++++++++++++++++++++

.. autosummary::
   :toctree: generated/

   acf

Data smoothing
++++++++++++++

.. autosummary::
   :toctree: generated/

   savgol_filter
   kernel_smoothing

Data selection
++++++++++++++

.. autosummary::
   :toctree: generated/

   time_to_sample
   get_closest_sampleWorkflow
--------

:py:mod:`frites.workflow`:

.. currentmodule:: frites.workflow

.. automodule:: frites.workflow
   :no-members:
   :no-inherited-members:

Task-related workflows
++++++++++++++++++++++

.. autosummary::
   :toctree: generated/

   WfMi
   WfMiCombine

Connectivity workflows
++++++++++++++++++++++

.. autosummary::
   :toctree: generated/

   WfConnComod


Statistical workflows
+++++++++++++++++++++

.. autosummary::
   :toctree: generated/

   WfStatsStatistics
----------

:py:mod:`frites.stats`:

.. currentmodule:: frites.stats

.. automodule:: frites.stats
   :no-members:
   :no-inherited-members:

Non-parametric statistics
+++++++++++++++++++++++++

.. autosummary::
   :toctree: generated/

   permute_mi_vector
   permute_mi_trials
   bootstrap_partitions

Random-effect (rfx)
+++++++++++++++++++

.. autosummary::
   :toctree: generated/

   ttest_1samp
   rfx_ttest

Correction for multiple comparisons
+++++++++++++++++++++++++++++++++++

.. autosummary::
   :toctree: generated/

   cluster_correction_mcp
   testwise_correction_mcp
   cluster_thresholdData simulations
----------------

:py:mod:`frites.simulations`:

.. currentmodule:: frites.simulations

.. automodule:: frites.simulations
   :no-members:
   :no-inherited-members:

Stimulus-specific autoregressive model
++++++++++++++++++++++++++++++++++++++

.. autosummary::
   :toctree: generated/

   StimSpecAR


Single and multi-subjects gaussian-based simulations
++++++++++++++++++++++++++++++++++++++++++++++++++++

.. autosummary::
   :toctree: generated/

   sim_ground_truth
   sim_local_cc_ss
   sim_local_cc_ms
   sim_local_cd_ss
   sim_local_cd_ms
   sim_local_ccd_ms
   sim_local_ccd_ssInformation-based estimators
----------------------------

:py:mod:`frites.estimator`:

.. currentmodule:: frites.estimator

.. automodule:: frites.estimator
   :no-members:
   :no-inherited-members:

Information-theoretic estimators
++++++++++++++++++++++++++++++++

.. autosummary::
   :toctree: generated/

   GCMIEstimator
   BinMIEstimator

Distance estimators
+++++++++++++++++++

.. autosummary::
   :toctree: generated/

   DcorrEstimator
   CorrEstimator

Estimator utility functions
+++++++++++++++++++++++++++

.. autosummary::
   :toctree: generated/

   CustomEstimator
   ResamplingEstimatorAPI
===

This section contains the list modules, classes and functions of Frites.

Overview
++++++++

.. panels::
    :container: container-lg pb-1
    :header: text-center
    :card: shadow

    Dataset
    ^^^^^^^

    The :mod:`frites.dataset` module contains data containers (a similar idea
    to MNE-Python with objects like :mod:`mne.Raw` or :mod:`mne.Epochs`) for
    single and multi subjects

    ---

    Workflows
    ^^^^^^^^^

    The :mod:`frites.workflow` module contains automated pipelines that are
    usually composed of two steps : `(i)` a first step consisting in estimating
    a measure of effect size (e.g. modulation of neural activity according to
    task conditions) followed by `(ii)` a statistical layer

    ---
    :column: col-lg-12 p-2

    Connectivity
    ^^^^^^^^^^^^

    The :mod:`frites.conn` module contains functions to estimate undirected and
    directed functional connectivity, potentially dynamic, at the single trial
    level

    ---
    :column: col-lg-12 p-2

    Estimators
    ^^^^^^^^^^

    The :mod:`frites.estimator` module contains information-based estimators
    for linking brain data to an external variable (e.g. stimulus type,
    behavioral models etc.). This includes metrics from the information-theory,
    machine-learning or measures of distances

    ---

    Simulated data
    ^^^^^^^^^^^^^^

    The :mod:`frites.simulations` module contains functions and classes to
    generate simulated data (single or multi-subjects)

    ---

    Utils
    ^^^^^

    The :mod:`frites.utils` module contains utility functions (e.g. data
    smoothing)


List of classes and functions
+++++++++++++++++++++++++++++

High-level API for users
^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 2

    api_dataset
    api_workflow
    api_connectivity
    api_estimators
    api_simulations
    api_utils

Low-level API for developers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 2

    api_statistics
    api_core
    api_conf
﻿frites.utils.kernel\_smoothing
==============================

.. currentmodule:: frites.utils

.. autofunction:: kernel_smoothing

.. _sphx_glr_backreferences_frites.utils.kernel_smoothing:

.. minigallery:: frites.utils.kernel_smoothing
    :add-heading:﻿frites.utils.time\_to\_sample
=============================

.. currentmodule:: frites.utils

.. autofunction:: time_to_sample

.. _sphx_glr_backreferences_frites.utils.time_to_sample:

.. minigallery:: frites.utils.time_to_sample
    :add-heading:﻿frites.simulations.sim\_local\_cd\_ss
=====================================

.. currentmodule:: frites.simulations

.. autofunction:: sim_local_cd_ss

.. _sphx_glr_backreferences_frites.simulations.sim_local_cd_ss:

.. minigallery:: frites.simulations.sim_local_cd_ss
    :add-heading:﻿frites.core.gccmi\_model\_nd\_cdnd
==================================

.. currentmodule:: frites.core

.. autofunction:: gccmi_model_nd_cdnd

.. _sphx_glr_backreferences_frites.core.gccmi_model_nd_cdnd:

.. minigallery:: frites.core.gccmi_model_nd_cdnd
    :add-heading:﻿frites.core.gccmi\_1d\_ccc
==========================

.. currentmodule:: frites.core

.. autofunction:: gccmi_1d_ccc

.. _sphx_glr_backreferences_frites.core.gccmi_1d_ccc:

.. minigallery:: frites.core.gccmi_1d_ccc
    :add-heading:﻿frites.core.cmi\_nd\_ggg
========================

.. currentmodule:: frites.core

.. autofunction:: cmi_nd_ggg

.. _sphx_glr_backreferences_frites.core.cmi_nd_ggg:

.. minigallery:: frites.core.cmi_nd_ggg
    :add-heading:﻿frites.simulations.sim\_local\_ccd\_ss
======================================

.. currentmodule:: frites.simulations

.. autofunction:: sim_local_ccd_ss

.. _sphx_glr_backreferences_frites.simulations.sim_local_ccd_ss:

.. minigallery:: frites.simulations.sim_local_ccd_ss
    :add-heading:﻿frites.conn.plot\_windows
=========================

.. currentmodule:: frites.conn

.. autofunction:: plot_windows

.. _sphx_glr_backreferences_frites.conn.plot_windows:

.. minigallery:: frites.conn.plot_windows
    :add-heading:﻿frites.core.gcmi\_model\_nd\_cd
===============================

.. currentmodule:: frites.core

.. autofunction:: gcmi_model_nd_cd

.. _sphx_glr_backreferences_frites.core.gcmi_model_nd_cd:

.. minigallery:: frites.core.gcmi_model_nd_cd
    :add-heading:﻿frites.conn.conn\_dfc
=====================

.. currentmodule:: frites.conn

.. autofunction:: conn_dfc

.. _sphx_glr_backreferences_frites.conn.conn_dfc:

.. minigallery:: frites.conn.conn_dfc
    :add-heading:﻿frites.core.gcmi\_nd\_cc
========================

.. currentmodule:: frites.core

.. autofunction:: gcmi_nd_cc

.. _sphx_glr_backreferences_frites.core.gcmi_nd_cc:

.. minigallery:: frites.core.gcmi_nd_cc
    :add-heading:﻿frites.estimator.CorrEstimator
==============================


.. currentmodule:: frites.estimator

.. autoclass:: CorrEstimator
   :special-members: __contains__,__getitem__,__iter__,__len__,__add__,__sub__,__mul__,__div__,__neg__,__hash__
   :members:

.. _sphx_glr_backreferences_frites.estimator.CorrEstimator:

.. minigallery:: frites.estimator.CorrEstimator
    :add-heading:﻿frites.set\_config
==================

.. currentmodule:: frites

.. autofunction:: set_config

.. _sphx_glr_backreferences_frites.set_config:

.. minigallery:: frites.set_config
    :add-heading:﻿frites.dataset.SubjectEphy
==========================


.. currentmodule:: frites.dataset

.. autoclass:: SubjectEphy
   :special-members: __contains__,__getitem__,__iter__,__len__,__add__,__sub__,__mul__,__div__,__neg__,__hash__
   :members:

.. _sphx_glr_backreferences_frites.dataset.SubjectEphy:

.. minigallery:: frites.dataset.SubjectEphy
    :add-heading:﻿frites.workflow.WfConnComod
===========================


.. currentmodule:: frites.workflow

.. autoclass:: WfConnComod
   :special-members: __contains__,__getitem__,__iter__,__len__,__add__,__sub__,__mul__,__div__,__neg__,__hash__
   :members:

.. _sphx_glr_backreferences_frites.workflow.WfConnComod:

.. minigallery:: frites.workflow.WfConnComod
    :add-heading:﻿frites.core.mi\_model\_nd\_gd
=============================

.. currentmodule:: frites.core

.. autofunction:: mi_model_nd_gd

.. _sphx_glr_backreferences_frites.core.mi_model_nd_gd:

.. minigallery:: frites.core.mi_model_nd_gd
    :add-heading:﻿frites.core.copnorm\_cat\_nd
============================

.. currentmodule:: frites.core

.. autofunction:: copnorm_cat_nd

.. _sphx_glr_backreferences_frites.core.copnorm_cat_nd:

.. minigallery:: frites.core.copnorm_cat_nd
    :add-heading:﻿frites.simulations.sim\_local\_cc\_ss
=====================================

.. currentmodule:: frites.simulations

.. autofunction:: sim_local_cc_ss

.. _sphx_glr_backreferences_frites.simulations.sim_local_cc_ss:

.. minigallery:: frites.simulations.sim_local_cc_ss
    :add-heading:﻿frites.utils.savgol\_filter
===========================

.. currentmodule:: frites.utils

.. autofunction:: savgol_filter

.. _sphx_glr_backreferences_frites.utils.savgol_filter:

.. minigallery:: frites.utils.savgol_filter
    :add-heading:﻿frites.core.copnorm\_cat\_1d
============================

.. currentmodule:: frites.core

.. autofunction:: copnorm_cat_1d

.. _sphx_glr_backreferences_frites.core.copnorm_cat_1d:

.. minigallery:: frites.core.copnorm_cat_1d
    :add-heading:﻿frites.stats.rfx\_ttest
=======================

.. currentmodule:: frites.stats

.. autofunction:: rfx_ttest

.. _sphx_glr_backreferences_frites.stats.rfx_ttest:

.. minigallery:: frites.stats.rfx_ttest
    :add-heading:﻿frites.core.gcmi\_model\_1d\_cd
===============================

.. currentmodule:: frites.core

.. autofunction:: gcmi_model_1d_cd

.. _sphx_glr_backreferences_frites.core.gcmi_model_1d_cd:

.. minigallery:: frites.core.gcmi_model_1d_cd
    :add-heading:﻿frites.conn.conn\_covgc
=======================

.. currentmodule:: frites.conn

.. autofunction:: conn_covgc

.. _sphx_glr_backreferences_frites.conn.conn_covgc:

.. minigallery:: frites.conn.conn_covgc
    :add-heading:﻿frites.conn.conn\_reshape\_undirected
=====================================

.. currentmodule:: frites.conn

.. autofunction:: conn_reshape_undirected

.. _sphx_glr_backreferences_frites.conn.conn_reshape_undirected:

.. minigallery:: frites.conn.conn_reshape_undirected
    :add-heading:﻿frites.core.mi\_mixture\_1d\_gd
===============================

.. currentmodule:: frites.core

.. autofunction:: mi_mixture_1d_gd

.. _sphx_glr_backreferences_frites.core.mi_mixture_1d_gd:

.. minigallery:: frites.core.mi_mixture_1d_gd
    :add-heading:﻿frites.estimator.CustomEstimator
================================


.. currentmodule:: frites.estimator

.. autoclass:: CustomEstimator
   :special-members: __contains__,__getitem__,__iter__,__len__,__add__,__sub__,__mul__,__div__,__neg__,__hash__
   :members:

.. _sphx_glr_backreferences_frites.estimator.CustomEstimator:

.. minigallery:: frites.estimator.CustomEstimator
    :add-heading:﻿frites.core.ent\_1d\_g
======================

.. currentmodule:: frites.core

.. autofunction:: ent_1d_g

.. _sphx_glr_backreferences_frites.core.ent_1d_g:

.. minigallery:: frites.core.ent_1d_g
    :add-heading:﻿frites.stats.ttest\_1samp
=========================

.. currentmodule:: frites.stats

.. autofunction:: ttest_1samp

.. _sphx_glr_backreferences_frites.stats.ttest_1samp:

.. minigallery:: frites.stats.ttest_1samp
    :add-heading:﻿frites.conn.conn\_get\_pairs
============================

.. currentmodule:: frites.conn

.. autofunction:: conn_get_pairs

.. _sphx_glr_backreferences_frites.conn.conn_get_pairs:

.. minigallery:: frites.conn.conn_get_pairs
    :add-heading:﻿frites.utils.get\_closest\_sample
=================================

.. currentmodule:: frites.utils

.. autofunction:: get_closest_sample

.. _sphx_glr_backreferences_frites.utils.get_closest_sample:

.. minigallery:: frites.utils.get_closest_sample
    :add-heading:﻿frites.simulations.sim\_local\_ccd\_ms
======================================

.. currentmodule:: frites.simulations

.. autofunction:: sim_local_ccd_ms

.. _sphx_glr_backreferences_frites.simulations.sim_local_ccd_ms:

.. minigallery:: frites.simulations.sim_local_ccd_ms
    :add-heading:﻿frites.stats.cluster\_threshold
===============================

.. currentmodule:: frites.stats

.. autofunction:: cluster_threshold

.. _sphx_glr_backreferences_frites.stats.cluster_threshold:

.. minigallery:: frites.stats.cluster_threshold
    :add-heading:﻿frites.core.copnorm\_1d
=======================

.. currentmodule:: frites.core

.. autofunction:: copnorm_1d

.. _sphx_glr_backreferences_frites.core.copnorm_1d:

.. minigallery:: frites.core.copnorm_1d
    :add-heading:﻿frites.core.cmi\_1d\_ggg
========================

.. currentmodule:: frites.core

.. autofunction:: cmi_1d_ggg

.. _sphx_glr_backreferences_frites.core.cmi_1d_ggg:

.. minigallery:: frites.core.cmi_1d_ggg
    :add-heading:﻿frites.core.gcmi\_1d\_cc
========================

.. currentmodule:: frites.core

.. autofunction:: gcmi_1d_cc

.. _sphx_glr_backreferences_frites.core.gcmi_1d_cc:

.. minigallery:: frites.core.gcmi_1d_cc
    :add-heading:﻿frites.core.mi\_1d\_gg
======================

.. currentmodule:: frites.core

.. autofunction:: mi_1d_gg

.. _sphx_glr_backreferences_frites.core.mi_1d_gg:

.. minigallery:: frites.core.mi_1d_gg
    :add-heading:﻿frites.conn.conn\_ravel\_directed
=================================

.. currentmodule:: frites.conn

.. autofunction:: conn_ravel_directed

.. _sphx_glr_backreferences_frites.conn.conn_ravel_directed:

.. minigallery:: frites.conn.conn_ravel_directed
    :add-heading:﻿frites.core.mi\_model\_1d\_gd
=============================

.. currentmodule:: frites.core

.. autofunction:: mi_model_1d_gd

.. _sphx_glr_backreferences_frites.core.mi_model_1d_gd:

.. minigallery:: frites.core.mi_model_1d_gd
    :add-heading:﻿frites.simulations.sim\_local\_cc\_ms
=====================================

.. currentmodule:: frites.simulations

.. autofunction:: sim_local_cc_ms

.. _sphx_glr_backreferences_frites.simulations.sim_local_cc_ms:

.. minigallery:: frites.simulations.sim_local_cc_ms
    :add-heading:﻿frites.estimator.BinMIEstimator
===============================


.. currentmodule:: frites.estimator

.. autoclass:: BinMIEstimator
   :special-members: __contains__,__getitem__,__iter__,__len__,__add__,__sub__,__mul__,__div__,__neg__,__hash__
   :members:

.. _sphx_glr_backreferences_frites.estimator.BinMIEstimator:

.. minigallery:: frites.estimator.BinMIEstimator
    :add-heading:﻿frites.core.gcmi\_mixture\_1d\_cd
=================================

.. currentmodule:: frites.core

.. autofunction:: gcmi_mixture_1d_cd

.. _sphx_glr_backreferences_frites.core.gcmi_mixture_1d_cd:

.. minigallery:: frites.core.gcmi_mixture_1d_cd
    :add-heading:﻿frites.conn.conn\_transfer\_entropy
===================================

.. currentmodule:: frites.conn

.. autofunction:: conn_transfer_entropy

.. _sphx_glr_backreferences_frites.conn.conn_transfer_entropy:

.. minigallery:: frites.conn.conn_transfer_entropy
    :add-heading:﻿frites.estimator.GCMIEstimator
==============================


.. currentmodule:: frites.estimator

.. autoclass:: GCMIEstimator
   :special-members: __contains__,__getitem__,__iter__,__len__,__add__,__sub__,__mul__,__div__,__neg__,__hash__
   :members:

.. _sphx_glr_backreferences_frites.estimator.GCMIEstimator:

.. minigallery:: frites.estimator.GCMIEstimator
    :add-heading:﻿frites.workflow.WfMiCombine
===========================


.. currentmodule:: frites.workflow

.. autoclass:: WfMiCombine
   :special-members: __contains__,__getitem__,__iter__,__len__,__add__,__sub__,__mul__,__div__,__neg__,__hash__
   :members:

.. _sphx_glr_backreferences_frites.workflow.WfMiCombine:

.. minigallery:: frites.workflow.WfMiCombine
    :add-heading:﻿frites.conn.conn\_reshape\_directed
===================================

.. currentmodule:: frites.conn

.. autofunction:: conn_reshape_directed

.. _sphx_glr_backreferences_frites.conn.conn_reshape_directed:

.. minigallery:: frites.conn.conn_reshape_directed
    :add-heading:﻿frites.simulations.StimSpecAR
=============================


.. currentmodule:: frites.simulations

.. autoclass:: StimSpecAR
   :special-members: __contains__,__getitem__,__iter__,__len__,__add__,__sub__,__mul__,__div__,__neg__,__hash__
   :members:

.. _sphx_glr_backreferences_frites.simulations.StimSpecAR:

.. minigallery:: frites.simulations.StimSpecAR
    :add-heading:﻿frites.conn.define\_windows
===========================

.. currentmodule:: frites.conn

.. autofunction:: define_windows

.. _sphx_glr_backreferences_frites.conn.define_windows:

.. minigallery:: frites.conn.define_windows
    :add-heading:﻿frites.workflow.WfStats
=======================


.. currentmodule:: frites.workflow

.. autoclass:: WfStats
   :special-members: __contains__,__getitem__,__iter__,__len__,__add__,__sub__,__mul__,__div__,__neg__,__hash__
   :members:

.. _sphx_glr_backreferences_frites.workflow.WfStats:

.. minigallery:: frites.workflow.WfStats
    :add-heading:﻿frites.stats.testwise\_correction\_mcp
======================================

.. currentmodule:: frites.stats

.. autofunction:: testwise_correction_mcp

.. _sphx_glr_backreferences_frites.stats.testwise_correction_mcp:

.. minigallery:: frites.stats.testwise_correction_mcp
    :add-heading:﻿frites.dataset.DatasetEphy
==========================


.. currentmodule:: frites.dataset

.. autoclass:: DatasetEphy
   :special-members: __contains__,__getitem__,__iter__,__len__,__add__,__sub__,__mul__,__div__,__neg__,__hash__
   :members:

.. _sphx_glr_backreferences_frites.dataset.DatasetEphy:

.. minigallery:: frites.dataset.DatasetEphy
    :add-heading:﻿frites.get\_config
==================

.. currentmodule:: frites

.. autofunction:: get_config

.. _sphx_glr_backreferences_frites.get_config:

.. minigallery:: frites.get_config
    :add-heading:﻿frites.estimator.ResamplingEstimator
====================================


.. currentmodule:: frites.estimator

.. autoclass:: ResamplingEstimator
   :special-members: __contains__,__getitem__,__iter__,__len__,__add__,__sub__,__mul__,__div__,__neg__,__hash__
   :members:

.. _sphx_glr_backreferences_frites.estimator.ResamplingEstimator:

.. minigallery:: frites.estimator.ResamplingEstimator
    :add-heading:﻿frites.core.mi\_nd\_gg
======================

.. currentmodule:: frites.core

.. autofunction:: mi_nd_gg

.. _sphx_glr_backreferences_frites.core.mi_nd_gg:

.. minigallery:: frites.core.mi_nd_gg
    :add-heading:﻿frites.core.copnorm\_nd
=======================

.. currentmodule:: frites.core

.. autofunction:: copnorm_nd

.. _sphx_glr_backreferences_frites.core.copnorm_nd:

.. minigallery:: frites.core.copnorm_nd
    :add-heading:﻿frites.core.gccmi\_nd\_ccnd
===========================

.. currentmodule:: frites.core

.. autofunction:: gccmi_nd_ccnd

.. _sphx_glr_backreferences_frites.core.gccmi_nd_ccnd:

.. minigallery:: frites.core.gccmi_nd_ccnd
    :add-heading:﻿frites.workflow.WfMi
====================


.. currentmodule:: frites.workflow

.. autoclass:: WfMi
   :special-members: __contains__,__getitem__,__iter__,__len__,__add__,__sub__,__mul__,__div__,__neg__,__hash__
   :members:

.. _sphx_glr_backreferences_frites.workflow.WfMi:

.. minigallery:: frites.workflow.WfMi
    :add-heading:﻿frites.estimator.DcorrEstimator
===============================


.. currentmodule:: frites.estimator

.. autoclass:: DcorrEstimator
   :special-members: __contains__,__getitem__,__iter__,__len__,__add__,__sub__,__mul__,__div__,__neg__,__hash__
   :members:

.. _sphx_glr_backreferences_frites.estimator.DcorrEstimator:

.. minigallery:: frites.estimator.DcorrEstimator
    :add-heading:﻿frites.core.gccmi\_1d\_ccd
==========================

.. currentmodule:: frites.core

.. autofunction:: gccmi_1d_ccd

.. _sphx_glr_backreferences_frites.core.gccmi_1d_ccd:

.. minigallery:: frites.core.gccmi_1d_ccd
    :add-heading:﻿frites.simulations.sim\_local\_cd\_ms
=====================================

.. currentmodule:: frites.simulations

.. autofunction:: sim_local_cd_ms

.. _sphx_glr_backreferences_frites.simulations.sim_local_cd_ms:

.. minigallery:: frites.simulations.sim_local_cd_ms
    :add-heading:﻿frites.stats.cluster\_correction\_mcp
=====================================

.. currentmodule:: frites.stats

.. autofunction:: cluster_correction_mcp

.. _sphx_glr_backreferences_frites.stats.cluster_correction_mcp:

.. minigallery:: frites.stats.cluster_correction_mcp
    :add-heading:﻿frites.core.gccmi\_nd\_ccc
==========================

.. currentmodule:: frites.core

.. autofunction:: gccmi_nd_ccc

.. _sphx_glr_backreferences_frites.core.gccmi_nd_ccc:

.. minigallery:: frites.core.gccmi_nd_ccc
    :add-heading: