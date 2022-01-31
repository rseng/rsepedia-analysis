# Welcome to pypillometry!

[![Build Status](https://travis-ci.com/ihrke/pypillometry.svg?branch=master)](https://travis-ci.com/ihrke/pypillometry)
[![GitHub issues](https://img.shields.io/github/issues/ihrke/pypillometry)](https://github.com/ihrke/pypillometry/issues)
[![status](https://joss.theoj.org/papers/3b06f4f3d5b703fd99c7e622b7edebe4/status.svg)](https://joss.theoj.org/papers/3b06f4f3d5b703fd99c7e622b7edebe4)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3925528.svg)](https://doi.org/10.5281/zenodo.3925528)


[![GitHub release](https://img.shields.io/github/release/ihrke/pypillometry.svg)](https://gitHub.com/ihrke/pypillometry/releases/)
[![PyPI version](https://badge.fury.io/py/pypillometry.svg)](https://badge.fury.io/py/pypillometry)


![](https://raw.githubusercontent.com/ihrke/pypillometry/master/logo/pypillometry_logo_200x200.png?token=AAIWMEINEM6MUOAPT2NV4I252K5QW)


This package implements functions for the analysis of pupillometric data. Features include preprocessing, blink handling, event-related pupil-dilation, plotting and signal modeling.

- Github-repository: <https://github.com/ihrke/pypillometry>

- Homepage/Documentation: <http://ihrke.github.com/pypillometry>


---

The project is maintained by Matthias Mittner <matthias.mittner@uit.no>

To report a bug or recommend a new feature, please use [github's issue tracker](https://github.com/ihrke/pypillometry/issues).

For more details about how to contribute, please see the [contribution-section](https://ihrke.github.io/pypillometry/html/docs/dev.html) in the documentation.---
title: 'pypillometry: A Python package for pupillometric analyses'
tags:
  - Python
  - pupil
  - pupillometry
  - psychology
  - neuroscience
  - LC-NE
authors:
  - name: Matthias Mittner
    orcid: 0000-0003-0205-7353
    affiliation: 1
affiliations:
 - name: Institute for Psychology, UiT The Arctic University of Norway, Norway
   index: 1
date: 19 May 2020
bibliography: paper.bib

---

# Summary

The size of the human pupil is controlled by pairs of constrictor and dilator muscles that allow its opening (dilation) and closing (constriction) in response to varying lighting conditions [@mathot2018pupillometry]. Importantly, it has long been known that the pupil also reacts to psychological important stimuli [@hess1960pupil] and has been a firmly established tool for studying "mental effort" in the research kit of psychologists for many decades [@laeng2012pupillometry]. More recently, pupil-size has been linked to the norepinephrinergic (NE) system originating from area *locus coeruleus* (LC) in the brainstem [@aston2005integrative], a link that has been substantiated experimentally by direct recordings in the brainstem of monkeys [@joshi2016relationships]. This finding of a correlation between NE activity in the brainstem and pupil-dilation has opened the way for researchers investigating the relationship between the LC-NE system and many cognitive functions, such as cognitive control [@gilzenrat2010pupil] and mind wandering [@mittner2016neural]. Advancing this emerging field requires the decomposition of the pupillometric signal into tonic (baseline) and phasic (response) components that relate to different processing regimes of the LC-NE system. 

The Python package `pypillometry` is a comprehensive library implementing preprocessing, plotting and advanced analysis tools in a coherent and extensible, object-oriented framework. It is oriented towards researchers in psychology and neuroscience that wish to analyze data from pupillometric experiments. `pypillometry` implements an intuitive, pipeline-based processing strategy where an analysis pipeline can be dynamically chained together. All operations and parameters applied to a dataset are stored in its history. This allows (1) a transparent and comprehensive logging of the operations applied for an analysis which is valuable for reproducible analyses, (2) the ability to "roll-back" any changes made to any point in the history of the dataset and (3) to easily generalize a processing pipeline to multiple datasets. The package contains pre-processing facilities implementing algorithms for blink detection and interpolation, filtering and resampling. All parameters are clearly documented, accessible and set to sensible default-parameters. A focus of the package is to provide extensive visualization features in order to facilitate dynamic exploration of the data as well as iterative adjustment of the pre-processing parameters. As the time-series of pupillometric data can be quite long, this requires separation into several plots or dynamically adjustable plot-axes. Both strategies are implemented in this package by allowing interactive plots if run from a Jupyter-Notebook [@kluyver2016jupyter] or storing a multi-page PDF document, allowing both interactive and scripted use. The `pypillometry` package also implements functions for event-related pupil-dilation (ERPD) analyses both at the individual and the group-level. Finally, the package implements novel algorithms for decomposing the pupillometric signal into tonic and phasic components. This approach allows to simultaneously quantify dynamic changes of both baseline and response-strength that can be related to the tonic and phasic processing regimes of the LC-NE system.

`Pypillometry` was already used for the analyses of several pupillometric datasets in our department. Several software packages with similar goals are available in R [e.g., @geller_gazer_2020; @forbes2020]. However, to date, no comprehensive Python-based solution besides `pypillometry` exists. None of the other packages provides facilities to estimate tonic and phasic components of the pupillometric signal. 

# Acknowledgements

I would like to thank the members of my research group for stimulating discussions and critical advice. In particular, I would like to thank Isabel Viola Kreis, Josephine Maria Groot, Gábor Csifcsák and Carole Schmidt for their input. I would also like to thank Bruno Laeng for inspiring discussions.

# References.. mdinclude:: ./README.md


Documentation
=============

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   docs/installation
   docs/usage
   docs/notebooks
   docs/api
   docs/dev
   docs/future

.. include:: /docs/notebooks.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
.. currentmodule:: pypillometry.pupildata

Overview
========
   
The package is divided into separate sub-modules that solve specific subtasks (pre-processing, baseline estimation etc). The functions in these modules operate on plain :py:mod:`numpy.array`'s or standard Python-structures. While these functions can be used directly, an easier approach is to use the class :py:class:`~pypillometry.pupildata.PupilData` which wraps these functions in a convenient way. Each object of class :py:class:`~pypillometry.pupildata.PupilData` represents one dataset of pupillometric data, including time, signal and external events (trial onsets, stimulus presentations, responses, etc). By calling the member functions of such an object, the corresponding function from one of the sub-modules is called using the appropriate arguments.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage

Reading/writing pupillometric data
----------------------------------

So far, reading in data is not part of the :py:mod:`pypillometry`-package. This is because the format of the eyetracker-files will vary depending on the setup of the eyetracker (there are many ways to represent event-triggers) and the actual model of the eyetracker. Python provides excellent functionality to parse text-based datafiles and we therefore give guidance how to use these tools rather than trying to implement that functionality in our package.

There are many ways in which pupillometric data can be read into Python. For example, Eyelink's `ASC-format <http://download.sr-support.com/dispdoc/page25.html>`_ generated by the EDF2ASC conversion tool outputs space-separated data that can be easily loaded using the `I/O functionality of the pandas package <https://pandas.pydata.org/pandas-docs/stable/reference/io.html>`_ . 

Data is input into :mod:`pypillometry` using the :func:`constructor <pypillometry.pupildata.PupilData.__init__>` of the :class:`~pypillometry.pupildata.PupilData` object.

However, data that has been converted into :class:`PupilData`-objects can be easily saved and restored (using :mod:`shelve`).

.. code-block:: python

    d=PupilData(...) # create PupilData object after manually loading data
    # save the dataset into a shelve-file
    d.write_file("dataset.pd")
    
    # load it from the file
    d2=PupilData.from_file("dataset.pd")
    


:ref:`An example for importing data from Eyelink EDF-files </docs/importdata.ipynb>`


Pipeline-based processing
-------------------------

:py:mod:`pypillometry` implements a pipeline-like approach where each operation executed on a :class:`~pypillometry.pupildata.PupilData`-object returns a copy of the (modified) object. This enables the "chaining" of commands as follows:

.. code-block:: python

    d=PupilData.from_file("data/test.pd").blinks_detect().blinks_merge().lowpass_filter(3).downsample(50)


This command loads a data-file (`test.pd`), applies a 3Hz low-pass filter to it, downsamples the signal to 50 Hz, detects blinks in the signal and merges short, successive blinks together. The final result of this processing-pipeline is stored in object `d`. This object stores also the complete history of the operations applied to the dataset and allows to transfer it to a new dataset.

See the following page more on this: :ref:`Pipeline-based processing in pypillometry </docs/pipes.ipynb>`

Pre-processing data
-------------------

Assuming you have generated a :class:`~pypillometry.pupildata.PupilData` object, a range of pre-processing functions are available. The main pre-processing issues with pupillometric data are:

- artifacts and missing data due to blinks (these can usually be corrected/interpolated)
- missing data/artifacts from other sources (e.g., looking away, eyetracker losing pupil for other reasons)
- smoothing/downsampling to get rid of high-freq low-amp noise



Handling Blinks
^^^^^^^^^^^^^^^

Pupillometric data usually contain blinks which show up as missing data in the signal where the eyetracker is unable to record the size of the pupil.
A range of functions are available for detecting and interpolating blinks. 

More details and an example can be found in the notebook: :ref:`An example for how to handle blinks </docs/blinks.ipynb>`

A fully worked-out example of a real study can be found in this notebook: :ref:`Preprocessing of a full dataset with multiple subjects </docs/preproc_example_pavlov.ipynb>`

The following is a list of functions for that purpose. Note that the functions take multiple arguments that control the algorithms behaviour. It is often crucial to adjust the parameters on an individual level since the artifacts tend to be quite dissimilar between subjects (but usually stable within-subject). All arguments are documented in the :ref:`API-docs </docs/api.rst>`.

.. autosummary::

    PupilData.blinks_detect
    PupilData.blinks_interpolate    
    PupilData.blinks_interp_mahot
    PupilData.blinks_merge
    PupilData.blinks_plot

Smoothing/low-pass filtering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In most cases, pupillometric data should be low-pass filtered (e.g., using a cutoff of 4 Hz `Jackson & Sirois, 2009 <https://doi.org/10.1111/j.1467-7687.2008.00805.x>`_) or smoothed in other ways (e.g., with a running-window). 

Tge following is a list of functions for smoothing:

.. autosummary::

    PupilData.lowpass_filter
    PupilData.smooth_window
    PupilData.downsample

Changing/Slicing data
^^^^^^^^^^^^^^^^^^^^^

Often, pupillometric data needs to be trimmed, e.g., to remove pre-experiment recordings or to remove unusable parts of the data (:func:`PupilData.sub_slice`). The timing should usually be realigned to the start of the experiment (:func:`PupilData.reset_time`). Furthermore, a scaling (e.g., Z-transform) of the pupil-data can be useful for comparing multiple subjects (:func:`PupilData.scale`).

The following is a list of available functions for these purposes:

.. autosummary::

    PupilData.sub_slice
    PupilData.copy
    PupilData.scale
    PupilData.unscale
    PupilData.reset_time

Plotting/Summarizing Data
-------------------------


Plotting
^^^^^^^^

It is crucial to validate preprocessing steps by visually inspecting the results using plots. Therefore, :mod:`pypillometry` implements several plotting facilities that encourage active exploration of the dataset. 

Please see the tutorial :ref:`Plotting of pupillometric data </docs/plotting.ipynb>` for more details.


.. autosummary::

    PupilData.plot
    PupilData.plot_segments
    PupilData.blinks_plot
    plotpd
    plotpd_ia
    PupilData.get_erpd



Inspecting/Summarizing
^^^^^^^^^^^^^^^^^^^^^^

The package also provides several functions for summarizing datasets. Simply `print()`ing a :class:`PupilData` object gives a readable summary of the main properties of the dataset and also prints the complete history of the results. By calling :func:`PupilData.summary`, summary data can be arranged and summarized in tabular form. 

See the notebook :ref:`Summarizing pupillometric data </docs/summary.ipynb>` for more details.

.. autosummary::

    PupilData.summary
    PupilData.stat_per_event 
    PupilData.get_erpd
    
Event-Related Pupil Dilation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The average pupillometric signal, timelocked to repeating events during the experiment, is referred to as "Event-related pupil dilation" or ERPD. In :mod:`pypillometry`, the functionality for this is implemented in the class :class:`pypillometry.erpd.ERPD`.

Running :func:`PupilData.get_erpd` returns an Object of class :class:`ERPDSingleSubject`. This object has functions for plotting and summarising the event-related pupillary dilation. 

:ref:`Here is an example notebook for for how to work with ERPDs </docs/erpds.ipynb>`.


.. autosummary::

    PupilData.get_erpd
    
.. currentmodule:: pypillometry.erpd

.. autosummary::

    group_erpd
    ERPD.summary
    ERPD.plot

    

Modeling the pupillometric signal
---------------------------------

For some applications, it is interesting to model the full pupillometric signal as consisting of a (tonic) baseline and a (phasic) response component. 
The package implements novel algorithms developed in our lab and documentation will become available here.

More details are availabel in this notebook: :ref:`Modeling the pupillometric signal </docs/modeling.ipynb>`.

.. currentmodule:: pypillometry.pupildata

.. autosummary::

    PupilData.estimate_baseline
    PupilData.estimate_response    
    PupilData.stat_per_event
    
Artificial Data
---------------

For validation and testing purposes, it can be useful to generate artificial datasets. The package implements a :class:`FakePupilData` as inheriting from regular :class:`PupilData` and therefore shares all its functionality. In addition to that, :class:`FakePupilData` stores "ground-truth" data and parameters that was used while creating the artificial data.

The function :func:`create_fake_pupildata` allows to quickly create datasets generating according to a provided experimental structure (see the functions documentation for an overview over the many available options).

.. autosummary::
    FakePupilData
    create_fake_pupildata


Tutorials/Example Notebooks
---------------------------

These examples are available as executable Jupyter Notebooks from the `Github repository <https://github.com/ihrke/pypillometry/tree/master/docs>`_ (the files ending in `.ipynb`). The link to the notebook files is given at the end of each page.

Each Notebook comes with a link that looks like this:

.. image:: https://mybinder.org/badge_logo.svg

When clicking on this link, the notebook will open in an interactive session in your browser without the need to install :py:mod:`pypillometry` at all.

.. toctree::
   :maxdepth: 2
   :caption: List of Notebooks:

   /docs/importdata
   /docs/pipes
   /docs/blinks
   /docs/plotting
   /docs/summary
   /docs/erpds
   /docs/modeling
   /docs/preproc_example_pavlovInstallation
============

Installing :mod:`pypillometry` and its dependencies is automated and can be done by running the following lines (on Mac OS X or Linux). 

.. code-block:: bash

    $ git clone https://github.com/ihrke/pypillometry.git
    $ cd pypillometry
    $ pip install -r requirements.txt
    $ python setup.py install

:mod:`pypillometry` is on `PyPI <https://pypi.org/>`_ and released versions can be installed with `pip` (this will also install the dependencies automatically):

.. code-block:: bash

    $ pip install pypillometry

(`link to the PyPI project page <https://pypi.org/project/pypillometry/>`_).

It is also possible to install the developer's version directly from github using `pip`

.. code-block:: bash

    $ pip install git+https://github.com/ihrke/pypillometry.git


Requirements
------------

:mod:`pypillometry` requires Python3 and a range of standard numerical computing packages (all of which listed in the file `requirements.txt`)

- :mod:`numpy`, :mod:`scipy` and :mod:`matplotlib`
- :mod:`pystan` 

It is useful to access :mod:`pypillometry` through Jupyter or Jupyter Notebook, so installing those packages is also useful but not necessary.

All requirements can be installed by running `pip install -r requirements.txt`.

Virtual environments
--------------------

It can sometimes be useful to install a new package in a new virtual environment using either `Python's virtual environments <https://docs.python.org/3/tutorial/venv.html>`_ or `conda <https://docs.conda.io/en/latest/>`_. 


.. code-block:: bash

    $ conda create -n pypil python=3
    $ conda activate pypil
    $ conda install anaconda 

The ``anaconda`` package contains all the requirements except :mod:`pystan` which can be installed from `conda-forge <https://anaconda.org/conda-forge/pystan>`_

.. code-block:: bash

    $ conda install -c conda-forge pystan


Pystan 
------

Note that the installation of :mod:`pystan` may cause trouble on Windows-systems (you may need to install a compiler). Please follow the instructions on `the Pystan-webpage <https://pystan.readthedocs.io/en/latest/getting_started.html>`_ should you encounter any trouble.


Notes/Potential Problems
-------------------------

Under Linux, I encountered a problem where :mod:`pystan` crashed the `Jupyter kernel <https://jupyter.org/>`_.
To circumvent this issue, I needed to install :mod:`pystan` using

.. code-block:: bash

    $ pip install pystan
    $ conda install gcc_linux-64
    $ conda install gxx_linux-64

otherwise, there were random crashes of the jupyter kernel for some reason. 

On Mac OS X, I had some trouble getting the compiler to work with PyStan. See `this issue <https://github.com/stan-dev/pystan/issues/622#issuecomment-518825883>`_ for a solution that worked for me.

To enable interactive plotting widgets in jupyter notebook and jupyter lab, widgets need to be enabled in the notebook.

.. code-block:: bash

    $ conda install ipywidgets nodejs
    $ jupyter nbextension enable --py widgetsnbextension
    $ jupyter labextension install @jupyter-widgets/jupyterlab-manager

Future Plans/Todo
=================

- implement functionality for binocular data
    - substitute data missing in only one eye by using data from the other
    - implement averaging to reduce to monocoular data
- extend fake-data
    - fake blinks
    - other fake artifacts
- implement pupil-size correction due to eye-movements
- `Requested enhancements at GitHub <https://github.com/ihrke/pypillometry/labels/enhancement>`_.. currentmodule:: pypillometry.pupildata

Developers
==========

.. image:: https://travis-ci.com/ihrke/pypillometry.svg?branch=master
    :target: https://travis-ci.com/ihrke/pypillometry

This package is developed and maintained by the `Cognitive Neuroscience Research Group <http://uit.no/research/cognitive-neuroscience>`_ at the `University of Tromsø <http://uit.no>`_. We encourage everyone to become a member of the team and to contribute to the package's development, either by `reporting bugs  <https://github.com/ihrke/pypillometry/issues>`_, `providing enhancements <https://github.com/ihrke/pypillometry/pulls>`_ or otherwise.

How to contribute
-----------------

Development of this package happens on GitHub: https://github.com/ihrke/pypillometry

In order to contribute, please use

- `github's issue tracker <https://github.com/ihrke/pypillometry/issues>`_ for reporting bugs or proposing enhancements
- `github's pull-request (PR) feature <https://github.com/ihrke/pypillometry/pulls>`_ for contributing to the package and/or documentation

If you are unfamiliar with Github's system, you can read up about `collaborating with issues and pull-requests on Github <https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests>`_.


Coding style/package architecture
---------------------------------

We chose to put most functionality directly below :mod:`pypillometry`'s main class :class:`PupilData`. This is mostly for convenience,
as users can easily find and chain commands together (see the following page more on this: :ref:`Pipeline-based processing in pypillometry </docs/pipes.ipynb>`).

The implementation of the functionality, however, is coded in submodules using "stateless" functions that take :mod:`numpy`-arrays as arguments. So, in order to create a new function that changes a :class:`PupilData`-object, you will first want to implement a function working purely on :mod:`numpy`-arrays and then create a thin wrapper function below :class:`PupilData` that calls this low-level function with the correct arguments.

Unit-testing
------------

Every newly developed function should be provided with `Unit-testing functions <https://en.wikipedia.org/wiki/Unit_testing>`_. As a minimum, the function should be called with some reasonable arguments to ensure that the function does not crash. Better yet, the results of the functions should be tested against the desired output within these test functions. All unit-tests are located in `/pypillometry/tests/` and are run automatically for every commit to Github using `Travis CI <https://travis-ci.com/ihrke/pypillometry>`_. 

Building the packages documentation
-----------------------------------

This is only necessary when extending `pypillometry`'s documentation (i.e., updating the website at https://ihrke.github.io/pypillometry).

The package uses the `Sphinx documentation generator <https://www.sphinx-doc.org/>`_ to create this website. All source-files are located under `/docs/` (both `.rst` and `.ipynb` files are being used). In addition, the API-documentation is placed within each function or classes' docstring.

To compile the documentation, sphinx must be installed. In addition, a couple of non-standard sphinx extensions are being used:

.. code-block:: python

    ['sphinx_math_dollar',  "m2r",'sphinx_autodoc_typehints', 'nbsphinx']


These can be installed by using the following commands:


.. code-block:: bash

    # jupyter notebooks
    conda install -c conda-forge nbsphinx 
    # math
    conda install -c conda-forge sphinx-math-dollar
    # markdown 2 rst (for README)
    pip install m2r
    # for the typehints
    conda install -c conda-forge sphinx-autodoc-typehints

Finally, the documentation can be created by running

.. code-block:: bash

    make html
    
in the packages' root-directory. 

.. warning::

    don't run `make clean`, it will delete the source-files from `docs`


Creating Releases
-----------------


A release should be created whenever crucial bugs have been fixed or new functionality has been added.
:mod:`pypillometry` goes with version numbers `x.y.z`. Increment `z` for bug-fixes, `y` for new features and
`x` for releases that break backward-compatibility.

.. note::

    the process of uploading to PyPI has been automatized using Github actions; when a
    release is created on `Github-releases <https://github.com/ihrke/pypillometry/releases>`_, the file 
    `VERSION` is updated with the most recent tag
    (must be `PEP-440 <https://www.python.org/dev/peps/pep-0440/#version-scheme>`_-compliant)
    and a `PyPI <https://pypi.org/>`_ release is being issued

   

API documentation
==================

.. automodule:: pypillometry.pupildata
    :members:
    :special-members:

.. automodule:: pypillometry.erpd
    :members:

.. automodule:: pypillometry.io
    :members:

.. automodule:: pypillometry.preproc
    :members:

.. automodule:: pypillometry.fakedata
    :members:

.. automodule:: pypillometry.baseline
    :members:

.. automodule:: pypillometry.pupil
    :members:

.. automodule:: pypillometry.convenience
    :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`