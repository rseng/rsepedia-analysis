############################
Contributing guidelines
############################

We welcome any kind of contribution to our software, from simple comment or question to a full fledged `pull request <https://help.github.com/articles/about-pull-requests/>`_. Please read and follow our `Code of Conduct <CODE_OF_CONDUCT.rst>`_.

A contribution can be one of the following cases:

#. you have a question;
#. you think you may have found a bug (including unexpected behavior);
#. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

You have a question
*******************

#. use the search functionality `here <https://github.com/phenology/cgc/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/phenology/cgc/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the `SHA hashcode <https://help.github.com/articles/autolinked-references-and-urls/#commit-shas>`_ of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
#. apply relevant labels to the newly created issue.

You want to make some kind of change to the code base
*****************************************************

#. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
#. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
#. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions `here <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`__ and `here <https://help.github.com/articles/syncing-a-fork/>`__);
#. make sure the existing tests still work by running ``python setup.py test``;
#. add your own tests (if necessary);
#. update or expand the documentation;
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the CGC repository on GitHub;
#. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
###########
Change Log
###########

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.

[Unreleased]
************

[0.6.1] - 2021-12-17
********************

Fixed
-----
* Fixing README - to be used as long_description on PyPI

[0.6.0] - 2021-12-17
********************

Added
-----
* k-means refinement also return refined-cluster labels

Fixed
-----
* Fixed bug in calculate_cluster_features, affecting kmeans and the calculation of the tri-cluster averages for particular ordering of the dimensions
* Number of converged runs in tri-cluster is updated

Changed
-------
* Numerical parameter epsilon is removed, which should lead to some improvement in the algorithm when empty clusters are present
* The refined cluster averages are not computed anymore over co-/tri-cluster averages but over all corresponding elements
* Dropped non-Numba powered low-mem version of co-clustering

[0.5.0] - 2021-09-23
********************

Added
-----
* k-means implementation for tri-clustering
* utility functions to calculate cluster-based averages for tri-clustering

Changed
-------
* Best k value in k-means is now selected automatically using the Silhouette score

[0.4.0] - 2021-07-29
********************

Added
-----
* utility function to estimate memory peak for numpy-based coclustering
* utility function to calculate cluster-based averages
* added Dask-based tri-clustering implementation


Fixed
-----
* k-means setup is more robust with respect to setting the range of k values and the threshold on the variance
* calculation of k-means statistics is faster


Changed
-------
* new version of tri-clustering algorithm implemented, old version moved to legacy folder


[0.3.0] - 2021-04-30
********************

Fixed
-----

* Reduced memory footprint of low-memory Dask-based implementation
* Fixed error handling in high-performance Dask implementation


Changed
-------

* Dropped tests on Python 3.6, added tests for Python 3.9 (following Dask)


[0.2.1] - 2020-09-18
********************

Fixed
-----

* Solve dependency issue: fail to install requirements with `pip`


[0.2.0] - 2020-09-17
********************

Added
-----

* Low-memory version for numpy-based coclustering, significantly reducing the memory footprint of the code
* Numba-accelerated version of the low-memory version of the numpy-based co-clustering
* Results objects include input_parameters dictionary and other metadata

Fixed
-----

* Solve issue in increasingly large Dask graph for increasing iterations

Changed
-------

* Main calculator classes stores results in dedicated object

[0.1.1] - 2020-08-27
********************

Added
-----

* Cluster results of co-/tri-clustring are now serialized to a file

Fixed
-----

* Improved output
* Bug fix in selecting minimum error run in co- and tri-clustering

Changed
-------

* K-means now loop over multiple k-values

[0.1.0] - 2020-08-11
********************

Added
-----

* First version of the CGC package, including minimal docs and tests
###############################################################################
Contributor Covenant Code of Conduct
###############################################################################

Our Pledge
**********

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

Our Standards
*************

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
  advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

Our Responsibilities
********************

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

Scope
*****

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

Enforcement
***********

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at team-atlas@esciencecenter.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

Attribution
***********

This Code of Conduct is adapted from the `Contributor Covenant <https://www.contributor-covenant.org>`_, version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html
.. list-table::
   :widths: 25 25
   :header-rows: 1

   * - `fair-software.nl <https://fair-software.nl>`_ recommendations
     - Badges
   * - \1. Code repository
     - |GitHub Badge|
   * - \2. License
     - |License Badge|
   * - \3. Community Registry
     - |PyPI Badge|
   * - \4. Enable Citation
     - |Zenodo Badge|
   * - \5. Checklist
     - |CII Best Practices Badge|
   * - **Other best practices**
     -
   * - Continuous integration
     - |Python Build| |Python Publish|
   * - Documentation
     - |Documentation Status|

.. |GitHub Badge| image:: https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue
   :target: https://github.com/phenology/cgc
   :alt: GitHub Badge

.. |License Badge| image:: https://img.shields.io/github/license/phenology/cgc
   :target: https://github.com/phenology/cgc
   :alt: License Badge

.. |PyPI Badge| image:: https://img.shields.io/pypi/v/clustering-geodata-cubes.svg?colorB=blue
   :target: https://pypi.python.org/project/clustering-geodata-cubes/
   :alt: PyPI Badge

.. |Zenodo Badge| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3979172.svg
   :target: https://doi.org/10.5281/zenodo.3979172
   :alt: Zenodo Badge

.. |CII Best Practices Badge| image:: https://bestpractices.coreinfrastructure.org/projects/4167/badge
   :target: https://bestpractices.coreinfrastructure.org/projects/4167
   :alt: CII Best Practices Badge

.. |Python Build| image:: https://github.com/phenology/cgc/workflows/Build/badge.svg
   :target: https://github.com/phenology/cgc/actions?query=workflow%3A%22Build%22
   :alt: Python Build

.. |Python Publish| image:: https://github.com/phenology/cgc/workflows/Publish/badge.svg
   :target: https://github.com/phenology/cgc/actions?query=workflow%3A%22Publish%22
   :alt: Python Publish

.. |Documentation Status| image:: https://readthedocs.org/projects/cgc/badge/?version=latest
   :target: https://cgc.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

CGC: Clustering Geo-Data Cubes
==============================

The Clustering Geo-Data Cubes (CGC) package focuses on the needs of geospatial data scientists who require tools to make sense of multi-dimensional data cubes. It provides the functionality to perform **co-cluster** and **tri-cluster** analyses on both local and distributed systems.

Installation
------------

To install CGC, do:

.. code-block:: console

  pip install clustering-geodata-cubes

Alternatively, you can clone this repository and install it using `pip`:

.. code-block:: console

  git clone https://github.com/phenology/cgc.git
  cd cgc
  pip install .


Run tests (including coverage) with:

.. code-block:: console

  python setup.py test

Documentation
-------------

The project's full API documentation can be found `online <https://cgc.readthedocs.io/en/latest/>`_. Including:

- `Co-clustering <https://cgc.readthedocs.io/en/latest/coclustering.html>`_
- `Tri-clustering <https://cgc.readthedocs.io/en/latest/triclustering.html>`_
- `K-means refinement <https://cgc.readthedocs.io/en/latest/kmeans.html>`_
- `Utility Functions <https://cgc.readthedocs.io/en/latest/utils.html>`_

Examples of CGC applications on real geo-spatial data:

- `Co-clustering application <https://cgc-tutorial.readthedocs.io/en/latest/notebooks/coclustering.html>`_
- `Tri-clustering application <https://cgc-tutorial.readthedocs.io/en/latest/notebooks/triclustering.html>`_

Tutorial
--------

The tutorial of CGC can be found  `here <https://cgc-tutorial.readthedocs.io/en/latest/index.html>`_.


Contributing
------------

If you want to contribute to the development of cgc, have a look at the `contribution guidelines`_.

.. _contribution guidelines: https://github.com/phenology/cgc/tree/master/CONTRIBUTING.rst

License
-------

Copyright (c) 2020-2021,

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Credits
-------

The code has been developed as a collaborative effort between the `ITC, University of Twente`_ and
`the Netherlands eScience Center`_ within the generalization of the project
`High spatial resolution phenological modelling at continental scales`_.

.. _ITC, University of Twente: https://www.itc.nl
.. _High spatial resolution phenological modelling at continental scales: https://research-software.nl/projects/1334
.. _the Netherlands eScience Center: https://www.esciencecenter.nl

This package was created with `Cookiecutter <https://github.com/audreyr/cookiecutter>`_ and the
`NLeSC/python-template <https://github.com/NLeSC/python-template>`_.
Utility Functions
=================

API
---

.. currentmodule:: cgc.utils

.. autofunction:: mem_estimate_coclustering_numpy

.. autofunction:: calculate_cocluster_averages

.. autofunction:: calculate_tricluster_averages

.. autofunction::  calculate_cluster_featureTri-clustering
==============

Introduction
------------

The ``triclustering`` module provides a generalization of the co-clustering algorithm to three-dimensional arrays (see
Ref. [#]_). For geospatial data, tri-clustering analyses allow extending the search for similarity patterns in
data cubes, thus accounting for an extra dimension (the 'band' dimension) in addition to space and time.

Setup the Analysis
------------------

The tri-clustering analysis of a three-dimensional array ``Z``:

.. code-block:: python

    import numpy as np

    Z = np.array([[[1., 1., 2., 4.],
                   [1., 1., 2., 4.]],
                  [[5., 5., 8., 8.],
                   [5., 5., 8., 8.]],
                  [[6., 7., 8., 9.],
                   [6., 7., 9., 8.]]])

is setup by creating an instance of ``Triclustering``:

.. code-block:: python

    from cgc.triclustering import Triclustering
    
    tc = Triclustering(
        Z,  # data array (3D)
        nclusters_row=4,  # number of row clusters
        nclusters_col=3,  # number of column clusters
        nclusters_bnd=2,  # number of band clusters
        max_iterations=100,  # maximum number of iterations
        conv_threshold=1.e-5,  # error convergence threshold 
        nruns=10,  # number of differently-initialized runs
        output_filename='results.json'  # JSON file where to write output
    )

The input arguments of ``Triclustering`` are identical to the ``Coclustering`` ones (see :doc:`coclustering`) -
``nclusters_bnd`` is the only additional argument, which sets the maximum number of clusters along the 'band' dimension. 
Note that a lower number of clusters can be identified by the algorithm (some of the clusters may remain empty).

.. NOTE::
    The first axis of ``Z`` is assumed to represent the 'band' dimension.

Tri-clustering Implementations
------------------------------

Local (Numpy-based)
*******************

As for the co-clustering algorithm (see :doc:`coclustering`), multiple runs of the tri-clustering algorithm can be
efficiently computed in parallel using threads. In order to run the tri-clustering analysis using 4 threads:

.. code-block:: python

    results = tc.run_with_threads(nthreads=4)

Distributed (Dask-based)
************************

Also for the tri-clustering, analysis on distributed systems can be carried out using Dask (see also
:doc:`coclustering`). Once the connection to a Dask cluster is setup:

.. code-block:: python

    from dask.distributed import Client

    client = Client('tcp://daskscheduler:8786')  # connect to the Dask scheduler


the tri-clustering analysis is carried out as:

.. code-block:: python

    results = tc.run_with_dask(client)

Results
-------

The ``TriclusteringResults`` object returned by ``Triclustering.run_with_threads`` and ``Triclustering.run_with_dask``
contains the final row, column, and band cluster assignments (``results.row_clusters``, ``results.col_clusters``, and
``results.bnd_clusters``, respectively) as well as the approximation error of the tri-clustering (``results.error``).
Few other metadata are also present, including the input parameters employed to setup the analysis
(``results.input_parameters``).


API
---

.. currentmodule:: cgc.triclustering

.. autoclass:: Triclustering
    :members:
    :undoc-members:

.. autoclass:: TriclusteringResults

References
----------

.. [#] Xiaojing Wu, Raul Zurita-Milla, Emma Izquierdo Verdiguier, Menno-Jan Kraak, Triclustering Georeferenced Time
 Series for Analyzing Patterns of Intra-Annual Variability in Temperature, Annals of the American Association of
 Geographers 108, 71 (2018)

.. include:: ../README.rst.. CGC documentation master file, created by
   sphinx-quickstart on Thu Jun 21 11:07:11 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CGC's documentation!
===============================

.. toctree::
   :maxdepth: 1
   :caption: Contents

   readme
   coclustering
   triclustering
   kmeans
   utils
   Tutorials <https://cgc-tutorial.readthedocs.io>
Co-clustering
=============

Introduction
------------

The ``coclustering`` module provides the functionality to perform the co-clustering analysis of a positive data matrix
with real-valued elements. The code implements the Bregman block-average co-clustering (BBAC) algorithm from Ref. [#]_
and it was inspired by the Matlab `code`_ by Srujana Merugu and Arindam Banerjee.

.. _code: http://www.ideal.ece.utexas.edu/software.html

The code was designed for geospatial applications (see the `Tutorial`_ section for some examples), where the array
dimensions typically correspond to space and time.

.. _Tutorial: https://cgc-tutorial.readthedocs.io

Setup the Analysis
------------------

For an array ``Z``:

.. code-block:: python

    import numpy as np

    Z = np.array([[1., 1., 2., 4.],
                  [1., 1., 2., 4.],
                  [3., 3., 3., 5.]])

the co-clustering analysis is setup by initializing a ``Coclustering`` object:

.. code-block:: python

    from cgc.coclustering import Coclustering
    
    cc = Coclustering(
        Z,  # data matrix
        nclusters_row=2, # number of row clusters
        nclusters_col=3,  # number of column clusters
        max_iterations=100,  # maximum number of iterations
        conv_threshold=1.e-5,  # error convergence threshold 
        nruns=10,  # number of differently-initialized runs
        output_filename='results.json'  # JSON file where to write output
    )

Here, we have set the maximum number of row and column clusters to 2 and 3, respectively. However, a lower number of
clusters can be identified by the algorithm (some of the clusters may remain empty). The algorithm entails an iterative
procedure that is considered converged when the error of two consecutive iterations differs by less than a threshold
(the default value is 1.e-5).

Multiple runs should be performed in order to limit the influence of the choice of initial cluster assignment on the
result. A numerical parameter guarantees that no zero-valued arguments are encountered in the logarithm that appears in
the I-divergence expression, which is employed as an objective function. Results are (optionally) written to a JSON file.

Co-clustering Implementations
-----------------------------

Local (Numpy-based)
*******************

The first one, based on `Numpy`_, is suitable to run the algorithm on a single machine. To make efficient use of
architectures with multi-core CPUs, the various differently-initialized co-clustering runs can be executed as multiple
threads. They are, in fact, embarrassingly parallel tasks that require no communication between each other. The
co-clustering analysis is run using e.g. 4 threads as:

.. code-block:: python

    results = cc.run_with_threads(nthreads=4)

This first implementation makes use of (fast) matrix multiplications to calculate cluster-based properties, such as
averages and distances. However, if ``Z``'s dimensions are large, large auxiliary matrices need to be stored into
memory, so that the memory requirement of this implementation quickly becomes a bottleneck.

.. _Numpy: https://numpy.org

Local (Numpy-based), low-memory footprint
*****************************************

A second Numpy-based implementation makes use of an algorithm with a much lower memory footprint, and can be selected
with the optional flag ``low_memory``:

.. code-block:: python

    results = cc.run_with_threads(nthreads=4, low_memory=True)

The reduced memory requirement comes at a certain cost of performance. Fortunately, we also applied `Numba`_'s just-in-time
compilation feature in the ``low_memory`` option. Thanks to this feature, the performance cost is significantly reduced.

.. _Numba: https://numba.pydata.org

Distributed (Dask-based)
************************

An alternative implementation makes use of `Dask`_ and is thus suitable to run the co-clustering algorithm on
distributed systems (e.g. on a cluster of compute nodes). Dask arrays are employed to process the data in chunks, which
are distributed across the cluster. This approach is thus suitable to tackle large matrices that do not fit the memory
of a single node.

If a Dask cluster is already running, we can connect to it and run the co-clustering analysis in the following way:

.. code-block:: python

    from dask.distributed import Client

    client = Client('tcp://daskscheduler:8786')  # connect to the Dask scheduler
    results = cc.run_with_dask(client)
    
.. _Dask: https://dask.org

Dask clusters can be run on different types of distributed systems: clusters of nodes connected by SSH, HPC systems,
Kubernetes clusters on cloud services. A local Dask cluster (``LocalCluster``) allows one to make use of the same
framework but using the local (multi-core) CPU(s).

In a second Dask-based implementation, the various co-clustering runs are submitted to the Dask scheduler, which
distributes them across the cluster. This implementation, which is activated by setting ``low_memory=False``, is
experimental and it typically leads to very large memory usages.

Results
-------

The ``CoclusteringResults`` object returned by ``Coclustering.run_with_threads`` and ``Coclustering.run_with_dask``
contains the final row and column cluster assignments (``results.row_clusters`` and ``results.col_clusters``,
respectively) as well as the approximation error of the co-clustering (``results.error``). Few other metadata are also
present, including the input parameters employed to setup the analysis (``results.input_parameters``).

API
---

.. currentmodule:: cgc.coclustering

.. autoclass:: Coclustering
    :members:
    :undoc-members:

.. autoclass:: CoclusteringResults

References
----------

.. [#] Arindam Banerjee, Inderjit Dhillon, Joydeep Ghosh, Srujana Merugu, Dharmendra S. Modha, A Generalized Maximum Entropy Approach to Bregman Co-clustering and Matrix Approximation, Journal of Machine Learning Research 8, 1919 (2007)
K-means refinement
==================

Introduction
------------

The `Kmeans` module is an implementation of the `k-means clustering`_ to refine the results of a co-clustering or
tri-clustering calculation. This k-mean refinement allows identifying similarity patterns between co- or tri-clusters.
The following pre-defined features, computed over all elements belonging to the same co- or tri-cluster, are employed
for the k-means clustering:

#. Mean value;
#. Standard deviation;
#. Minimum value;
#. Maximum value;
#. 5th percentile;
#. 95th percentile;

The implementation, which is based on the `scikit-learn`_ package, tests a range of k values and select the optimal one
based on the `Silhouette coefficient`_.

.. _scikit-learn: https://scikit-learn.org/stable/index.html
.. _Silhouette coefficient: https://en.wikipedia.org/wiki/Silhouette_(clustering)

Running the refinement
----------------------

The k-means refinement should be based on existing co- or tri-clustering results:

.. code-block:: python

    import numpy as np

    Z = np.array([[1., 1., 2., 4.],
                  [1., 1., 2., 4.],
                  [3., 3., 3., 5.]])
    row_clusters = np.array([0, 0, 1, 2])  # 3 clusters
    col_cluster = np.array([0, 0, 1])  # 2 clusters

One can then setup ``Kmeans`` in the following way:

.. code-block:: python

    from cgc.kmeans import Kmeans

    km = Kmeans(
        Z,
        clusters=(row_clusters, col_cluster),
        nclusters=(3, 2)
        k_range=range(2, 5),
        kmean_max_iter=100,
        output_filename='results.json' # JSON file where to write output
    )

Here ``k_range`` is the range of ``k`` values to investigate. If not provided, a sensible range will be setup (from 2 to
a fraction of the number of co- or tri-clusters - the optional `max_k_ratio` argument allows for additional control, see
:ref:`API<API>`). ``kmean_max_iter`` is the maximum number of iterations employed for the k-means clustering.

The ``compute`` function is then called to run the k-means refinement:

.. code-block:: python

    results = km.compute()

Results
-------

The optimal ``k`` value and the refined cluster averages computed over all elements assigned to the co- and tri-clusters
are stored in the ``KmeansResults`` object:

.. code-block:: python

    results.k_value
    results.cluster_averages


.. _k-means clustering: https://en.wikipedia.org/wiki/K-means_clustering

.. _API:

API
---

.. currentmodule:: cgc.kmeans

.. autoclass:: Kmeans
    :members:
    :undoc-members:

.. autoclass:: KmeansResults
