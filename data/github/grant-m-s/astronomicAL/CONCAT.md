---
title: 'AstronomicAL: an interactive dashboard for visualisation, integration and classification of data with Active Learning'
tags:
 - Python
 - astronomy
 - machine learning
 - active learning
 - labelling
 - classification
 - software
authors:
 - name: Grant Stevens
   orcid: 0000-0002-8885-4443
   affiliation: 1
 - name: Sotiria Fotopoulou
   orcid: 0000-0002-9686-254X
   affiliation: 2
 - name: Malcolm N. Bremer
   affiliation: 2
 - name: Oliver Ray
   affiliation: 1

affiliations:
- name: Department of Computer Science, Merchant Venturers Building, University of Bristol, Woodland Road, Bristol, BS8 1UB
  index: 1
- name: School of Physics, HH Wills Physics Laboratory, University of Bristol, Tyndall Avenue, Bristol, BS8 1TL
  index: 2
date: 5 August 2021
bibliography: paper/paper.bib

---

# Summary

AstronomicAL is a human-in-the-loop interactive labelling and training dashboard that allows users to create reliable datasets and robust classifiers using active learning. This technique prioritises data that offer high information gain, leading to improved performance using substantially less data. The system allows users to visualise and integrate data from different sources and deal with incorrect or missing labels and imbalanced class sizes. AstronomicAL enables experts to visualise domain-specific plots and key information relating both to broader context and details of a point of interest drawn from a variety of data sources, ensuring reliable labels. In addition, AstronomicAL provides functionality to explore all aspects of the training process, including custom models and query strategies. This makes the software a tool for experimenting with both domain-specific classifications and more general-purpose machine learning strategies. We illustrate using the system with an astronomical dataset due to the field’s immediate need; however, AstronomicAL has been designed for datasets from any discipline. Finally, by exporting a simple configuration file, entire layouts, models, and assigned labels can be shared with the community. This allows for complete transparency and ensures that the process of reproducing results is effortless.

![AstronomicAL workflow from the perspective of an astronomy user. Each part of the workflow improves the reliability of labels, leading to enhanced performance of any classifiers produced. AstronomicAL allows the user to tailor this workflow specifically to their domain with its modular and extensible design.  \label{fig:workflow}](paper/astronomicAL_diagram_portrait.png){ width=86% }

# Statement of Need

Active learning has proven to be an effective method in machine learning over the last few decades [@settles2008analysis; @baldridge2009well; @gal2017deep]. In more recent years, there has been a growing adoption in all domains of scientific research [@walmsley2020galaxy; @eisenstein2020active]. This is evident by the increasing popularity of machine learning frameworks offering an active learning element^[https://github.com/topics/active-learning]. However, the purpose of these frameworks is to allow users to integrate active learning into their machine learning workflows. Any additional functionality such as data exploration and interactive labelling is left to the user to implement. Due to each domain's varying requirements, researchers will often need multiple specialised but disjoint tools that hinder effective and efficient interactive labelling. When data from multiple sources is required, presenting all the information in one system is essential. Without a tool that combines these key features, a significant barrier exists to the adoption and effectiveness of active learning in research. AstronomicAL's goal is to remove this barrier to allow more researchers to access the benefits of active learning.

Active learning [@settles2012active] removes the requirement for large amounts of labelled training data whilst still producing high accuracy models. This is extremely important as with ever-growing datasets; it is becoming impossible to manually inspect and verify ground truth used to train machine learning systems. The reliability of the training data limits the performance of any supervised learning model, so consistent classifications become more problematic as data sizes increase. The problem is exacerbated when a dataset does not contain any labelled data, preventing supervised learning techniques entirely. AstronomicAL has been developed to tackle these issues head-on and provide a solution for any large scientific dataset.

It is common for active learning to query areas of high uncertainty; these are often in the boundaries between classes where the expert’s knowledge is required. To facilitate this human-in-the-loop process, AstronomicAL provides users with the functionality to fully explore each data point chosen. This allows them to inject their domain expertise directly into the training process, ensuring that assigned labels are both accurate and reliable.


AstronomicAL has been extensively validated on astronomy datasets. These are highly representative of the issues that we anticipate will be found in other domains for which the tool is designed to be easily customisable. Such issues include the volume of data (millions of sources per survey), vastly imbalanced classes and ambiguous class definitions leading to inconsistent labelling. AstronomicAL has been developed to be sufficiently general for any tabular data and can be customised for any domain. For example, we provide the functionality for data fusion of catalogued data and online cutout services for astronomical datasets.

Using its modular and extensible design, researchers can quickly adapt AstronomicAL for their research to allow for domain-specific plots, novel query strategies, and improved models. Furthermore, there is no requirement to be well-versed in the underlying libraries that the software uses. This is due to large parts of the complexity being abstracted whilst allowing more experienced users to access full customisability.

As the software runs entirely locally on the user’s system, AstronomicAL provides a private space to experiment whilst providing a public mechanism to share results. By sharing only the configuration file, users remain in charge of distributing their potentially sensitive data, enabling collaboration whilst respecting privacy. The full workflow for AstronomicAL is shown in \autoref{fig:workflow}.

![The dashboard is set up in 3 main areas. A) Active learning and labelling panel which controls the machine learning functionality. Users can train separate models for each class and see how each performs as additional points are added to the model. B) User-defined information shown for the currently active data point, such as specific columns, local images and online cutout services. C) Customisable domain-specific plots that can render millions of data points. \label{fig:full_layout}](paper/full_layout_browser_box.png)

# Active Learning and Classification

Active learning utilises methods for predicting which part of the search space of unseen examples would be most likely to improve classification accuracy. By querying these key examples and presenting them directly to the user for inspection and labelling, active learning removes the reliance on non-verified data. In the optimal case, it provides users with simple classifications that have a high impact on model performance, leading to a dramatic reduction in the number of data points required for training a classifier.

By making use of the modAL [@modAL2018] active learning framework, AstronomicAL allows users to take full advantage of active learning techniques, leading to improved classifiers. Scikit-Learn [@scikit-learn] models are fully supported in AstronomicAL, and the software will be extended to support Pytorch [@pytorch] and Tensorflow [@tensorflow2015-whitepaper] in the future.

Astronomy datasets are frequently the result of the fusion of multiple heterogeneous data sources. Often in machine learning, only information that is available across all data points is usable during training. However, in astronomy, where *any* available information could be critical in determining an object’s true class, it must be used whenever this information is available. With the flexibility to choose which information is essential for classification, AstronomicAL ensures that users are presented with this critical information whenever available - directly leading to improved confidence and justifiability for any assigned labels. Although this extra information remains unseen by the model, by using AstronomicAL, we inject this information into the training process by labelling the data based on all the available information rather than just the features the model sees.

In \autoref{fig:full_layout} we show an example of a classification task using an extract of the astronomical dataset used in @Fotopoulou2018CPz. \autoref{fig:full_layout} (A) shows the benefits of active learning where with only 30 data points given as training data to the model, we can classify Stars with 97.5% accuracy. \autoref{fig:full_layout} (B) shows the extra information chosen by the user to improve confidence in the assigned labels. This information can be retrieved from web services or local files. In this example, we are showing the SDSS^[http://skyserver.sdss.org/dr16/en/help/docs/api.aspx#imgcutout] [@Ahumada2020DR16] and FIRST^[https://third.ucllnl.org/cgi-bin/firstcutout] [@Becker1995Radio] cutout services. \autoref{fig:full_layout} (C) shows how domain-specific plots can be created and visualised to gain extra insight into the data, aiding the researcher in their ability and confidence in labelling data points. \autoref{fig:active_learning} shows the model performance plots available to the user for each classifier being trained. These include training and validation correctness plots, tracking of performance scores as new data points are added, and a visualisation of the model's confidence of each data point.

![AstronomicAL allows users to view all key information during training. Top Left: Training set showing points that have been trained on, current queried point and correct and incorrect predictions. Top Right: Visualisation of the model’s confidence of each data point. Bottom Left: Validation set showing correct and incorrect predictions of the model. Bottom Right: The performance scores of the model as new points are added to the training set. \label{fig:active_learning}](paper/active_learning_1.png){ width=65% }

# Interactivity

Using the Panel dashboard library [@philipp_rudiger_panel_2021_4573728], paired with Holoviews [@philipp_rudiger_holoviews_2021_4581995], Bokeh [@Bokeh2018] and Datashader [@james_a_bednar_datashader_2020_3987379] visualisation libraries, AstronomicAL provides users with interactive plots, enabling zooming and panning of millions of data points in real-time whilst also giving the ability to rearrange and resize plots dynamically to optimise screen layout.

Plots are also interconnected, with any updates projected to all panels simultaneously, allowing analysis of how the properties of a particular point fit in the context of the whole dataset. This enables researchers to gain insight and identify trends - crucial for reliable labelling.

![AstronomicAL enables the creation of curated test sets. The sample is selected according to user-defined criteria and verified visually. This test set can then be used during the active learning training process to ensure that the model will be sufficiently generalisable to future data. \label{fig:labelling_mode}](paper/labelling_mode.png){ width=85% }

The interactive labelling functionality is not limited to only the training stage; as shown in \autoref{fig:labelling_mode}, users can curate a labelled test set to sufficiently demonstrate the validity and generalisability of their model. Furthermore, to facilitate and encourage only assigning labels that the user trusts, AstronomicAL allows users to mark any example as *unsure*. Such examples are excluded from the training set, ensuring that all training data are of high quality.

In summary, all design decisions have been in response to user feedback and were explicitly tailored to improve the visualisation of data, irrespective of the data’s specific nature, to ensure that AstronomicAL is a tool for any discipline.

# Acknowledgements

Grant Stevens acknowledges financial support from the UKRI for a Centre for Doctoral Training studentship in Interactive Artificial Intelligence at the University of Bristol. (EP/S022937/1)

# References
[![Build Status](https://travis-ci.com/grant-m-s/astronomicAL.svg?token=upRGxrMseZqj7kT3bSGx&branch=master)](https://travis-ci.com/grant-m-s/astronomicAL) [![codecov](https://codecov.io/gh/grant-m-s/astronomicAL/branch/master/graph/badge.svg?token=TCO9J2AD1Z)](https://codecov.io/gh/grant-m-s/astronomicAL) [![Documentation Status](https://readthedocs.org/projects/astronomical/badge/?version=latest)](https://astronomical.readthedocs.io/en/latest/?badge=latest)

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03635/status.svg)](https://doi.org/10.21105/joss.03635)

# AstronomicAL

## An interactive dashboard for visualisation, integration and classification of data using Active Learning.

AstronomicAL is a human-in-the-loop interactive labelling and training dashboard that allows users to create reliable datasets and robust classifiers using active learning. The system enables users to visualise and integrate data from different sources and deal with incorrect or missing labels and imbalanced class sizes by using active learning to help the user focus on correcting the labels of a few key examples. Combining the use of the [Panel](https://panel.holoviz.org/), [Bokeh](https://docs.bokeh.org/en/latest/index.html), [modAL](https://github.com/modAL-python/modAL) and [SciKit Learn](https://scikit-learn.org/stable/) packages, AstronomicAL enables researchers to take full advantage of the benefits of active learning: high accuracy models using just a fraction of the total data, without the requirement of being well versed in underlying libraries.

![Load Configuration](docs/source/images/AstronomicAL_demo.gif)

### Statement of Need

Active learning [(Settles, 2012)](https://www.morganclaypool.com/doi/abs/10.2200/S00429ED1V01Y201207AIM018) removes the requirement for large amounts of labelled training data whilst still producing high accuracy models. This is extremely important as with ever-growing datasets; it is becoming impossible to manually inspect and verify ground truth used to train machine learning systems. The reliability of the training data limits the performance of any supervised learning model, so consistent classifications become more problematic as data sizes increase. The problem is exacerbated when a dataset does not contain any labelled data, preventing supervised learning techniques entirely. AstronomicAL has been developed to tackle these issues head-on and provide a solution for any large scientific dataset.

It is common for active learning to query areas of high uncertainty; these are often in the boundaries between classes where the expert's knowledge is required. To facilitate this human-in-the-loop process, AstronomicAL provides users with the functionality to fully explore each data point chosen. This allows them to inject their domain expertise directly into the training process, ensuring that asigned labels are both accurate and reliable.

AstronomicAL has been extensively validated on astronomy datasets. These are highly representative of the issues that we anticipate will be found in other domains for which the tool is designed to be easily customisable. Such issues include the volume of data (millions of sources per survey), vastly imbalanced classes and ambiguous class definitions leading to inconsistent labelling. AstronomicAL has been developed to be sufficiently general for any tabular data and can be customised for any domain. For example, we provide the functionality for data fusion of catalogued data and online cutout services for astronomical datasets.

Using its modular and extensible design, researchers can quickly adapt AstronomicAL for their research to allow for domain-specific plots, novel query strategies, and improved models. Furthermore, there is no requirement to be well-versed in the underlying libraries that the software uses. This is due to large parts of the complexity being abstracted whilst allowing more experienced users to access full customisability.

As the software runs entirely locally on the user's system, AstronomicAL provides a private space to experiment whilst providing a public mechanism to share results. By sharing only the configuration file, users remain in charge of distributing their potentially sensitive data, enabling collaboration whilst respecting privacy.

### Documentation

The documentation for AstronomicAL can be found [here](https://astronomical.readthedocs.io).

## Installation

To install AstronomicAL and its dependencies, the user can clone the repository and from within the repo folder run `pip install -r requirements.txt`. . It is recommended that the user creates a virtual environment using tools such as [Virtualenv](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/#installing-virtualenv) or [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html), to prevent any conflicting package versions.

```
    git clone https://github.com/grant-m-s/AstronomicAL.git
    cd AstronomicAL
    conda config --add channels conda-forge
    conda create --name astronomical --file requirements.txt
    conda activate astronomical
```

### Quickstart Instructions

To begin using the software, run `bokeh serve astronomicAL --show` and your browser should automatically open to [localhost:5006/astronomicAL](localhost:5006/astronomicAL>`)

AstronomicAL provides both an example dataset and an example configuration file to allow you to jump right into the software and give it a test run.

![Load Configuration](docs/source/images/Load_config_AL.gif)

To begin training you simply have to select **Load Custom Configuration** checkbox and select your config file. Here we have chosen to use the `example_config.json` file.

The **Load Config Select** option allows use to choose the extent to which to reload the configuration.

## Contributing to AstronomicAL

### Reporting Bugs

If you encounter a bug, you can directly report it in the [issues section](https://github.com/grant-m-s/AstronomicAL/issues).

Please describe how to reproduce the bug and include as much information as possible that can be helpful for fixing it.

**Are you able to fix a bug?**

You can open a new pull request or include your suggested fix in the issue.

### Submission of extensions

**Have you created an extension that you want to share with the community?**

Create a pull request describing your extension and how it can improve research for others.

### Support and Feedback

We would love to hear your thoughts on AstronomicAL.

Are there any features that would improve the effectiveness and usability of AstronomicAL? Let us know!

Any feedback can be submitted as an [issue](https://github.com/grant-m-s/AstronomicAL/issues).

## Referencing the Package

Please remember to cite our software and user guide whenever relevant.

See the [Citing page](https://astronomical.readthedocs.io/en/latest/content/other/citing.html) in the documentation for instructions about referencing and citing the astronomicAL software.
Welcome to AstronomicAL's documentation!
========================================
An interactive dashboard for visualisation, integration and classification of data using Active Learning.
--------------------------------------------------------------------

.. image:: https://travis-ci.com/grant-m-s/astronomicAL.svg?token=upRGxrMseZqj7kT3bSGx&branch=master
    :target: https://travis-ci.com/grant-m-s/astronomicAL


.. image:: https://codecov.io/gh/grant-m-s/astronomicAL/branch/master/graph/badge.svg?token=TCO9J2AD1Z
    :target: https://codecov.io/gh/grant-m-s/astronomicAL

.. image:: https://readthedocs.org/projects/astronomical/badge/?version=latest
    :target: https://astronomical.readthedocs.io

.. image:: https://joss.theoj.org/papers/10.21105/joss.03635/status.svg
   :target: https://doi.org/10.21105/joss.03635

AstronomicAL is a human-in-the-loop interactive labelling and training dashboard that allows users to create reliable datasets and robust classifiers using active learning. The system enables users to visualise and integrate data from different sources and deal with incorrect or missing labels and imbalanced class sizes by using active learning to help the user focus on correcting the labels of a few key examples. Combining the use of the Panel_, Bokeh_, modAL_ and `SciKit Learn`_ packages, AstronomicAL enables researchers to take full advantage of the benefits of active learning: high accuracy models using just a fraction of the total data, without the requirement of being well versed in underlying libraries.

.. _Panel: https://panel.holoviz.org/
.. _Bokeh: https://docs.bokeh.org/en/latest/index.html
.. _modAL: https://github.com/modAL-python/modAL
.. _`SciKit Learn`: https://scikit-learn.org/stable/

.. figure:: images/AstronomicAL_demo.gif

Statement of Need
*****************
With ever-growing datasets, it is becoming impossible to manually inspect and verify ground truth used to train machine learning systems. The reliability of the training data limits the performance of any supervised learning model, so consistent classifications become more problematic as data sizes increase. The problem is exacerbated when a dataset does not contain any labelled data, preventing supervised learning techniques entirely. Active learning `(Settles, 2012)`_ addresses these issues by removing the requirement for large amounts of labelled training data whilst still producing high accuracy models.

Although initially designed for astronomers, by providing the functionality for data fusion of catalogued data and online cutout services, AstronomicAL has been developed to be sufficiently general for any tabular data. Large datasets, unreliable labels and vastly imbalanced classes make astronomy data the ideal vehicle to develop this software. Each of these issues is an examplar of more generalised problems that active learning could solve in any dataset.

Using its modular and extensible design, researchers can quickly adapt AstronomicAL for their research to allow for domain-specific plots, novel query strategies, and improved models. Further, there is no requirement to be well-versed in the underlying libraries that the software uses as large parts of the complexity are abstracted whilst still allowing more experienced users to access full customisability.

As the software runs entirely locally on the user’s system, AstronomicAL provides a private space to experiment whilst providing a public mechanism to share results. By sharing only the configuration file, users remain in charge of distributing their potentially sensitive data, enabling collaboration whilst respecting privacy.

.. _`(Settles, 2012)`: https://www.morganclaypool.com/doi/abs/10.2200/S00429ED1V01Y201207AIM018

Installation
------------------
To install AstronomicAL and its dependencies, the user can clone the repository and from within the repo folder run :code:`pip install -r requirements.txt`. It is recommended that the user creates a virtual environment using tools such as Virtualenv_ or Conda_, to prevent any conflicting package versions.

.. code-block:: bash

    git clone https://github.com/grant-m-s/AstronomicAL.git
    cd AstronomicAL
    conda config --add channels conda-forge
    conda create --name astronomical --file requirements.txt
    conda activate astronomical

.. _Virtualenv: https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/#installing-virtualenv
.. _Conda: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

Quickstart Instructions
-----------------------
To begin using the software, run :code:`bokeh serve astronomicAL --show`, and your browser should automatically open to `localhost:5006/astronomicAL
<localhost:5006/astronomicAL>`_.

AstronomicAL provides an example dataset and an example configuration file to allow you to jump right into the software and give it a test run.

.. figure:: images/Load_config_AL.gif

    AstronomicAL makes it easy to start training your classifier or reload a previous checkpoint.

To begin training, you simply have to select **Load Custom Configuration** checkbox and select your config file. Here we have chosen to use the :code:`example_config.json` file.

The **Load Config Select** option allows users to choose the extent to which to reload the configuration.

.. raw:: html

   <hr>


Contributing to AstronomicAL
-------------------------

Reporting Bugs
*****************

If you encounter a bug, you can directly report it in the `issues section <https://github.com/grant-m-s/AstronomicAL/issues>`_.

Please describe how to reproduce the bug and include as much information as possible that can be helpful for fixing it.

**Are you able to fix a bug?**

You can open a new pull request or include your suggested fix in the issue.

Submission of extensions
*****************

**Have you created an extension that you want to share with the community?**

Create a pull request describing your extension and how it can improve research for others.

Support and Feedback
*****************

We would love to hear your thoughts on AstronomicAL.

Are there any features that would improve the effectiveness and usability of AstronomicAL? Let us know!

Any feedback can be submitted as an `issue <https://github.com/grant-m-s/AstronomicAL/issues>`_.

.. raw:: html

   <hr>

Referencing the Package
-------------------------

Please remember to cite our software and user guide whenever relevant.

See the :ref:`citing <citing>` page for instructions about referencing and citing the AstronomicAL software.


.. raw:: html

   <hr>

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. toctree::
    :glob:
    :maxdepth: 1
    :caption: API reference

    content/apireference/active_learning.rst
    content/apireference/dashboard.rst
    .. content/apireference/extensions.rst
    content/apireference/settings.rst
    .. content/apireference/utils.rst

.. toctree::
    :maxdepth: 1
    :caption: Tutorials

    content/tutorials/preparing_dataset.rst
    content/tutorials/settings.rst
    content/tutorials/active_learning.rst
    content/tutorials/reload_config.rst
    content/tutorials/labelling_test_set.rst
    content/tutorials/plots.rst
    content/tutorials/feature_generation.rst
    content/tutorials/using_model.rst

.. toctree::
    :maxdepth: 1
    :caption: Other

    content/other/contributors.rst
    content/other/citing.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. image:: images/CDT-UOB-logo.png
   :width: 38%
   :target: http://www.bristol.ac.uk/cdt/interactive-ai/


.. image:: images/EPSRC+logo.png
   :width: 56%
   :target: https://gtr.ukri.org/projects?ref=studentship-2466020
Contributors
========================


Developers
------------------------

Core Development
**********************

* `Grant Stevens <https://research-information.bris.ac.uk/en/persons/grant-stevens>`_

Technical Advisors
**********************

* `Sotiria Fotopoulou <https://research-information.bris.ac.uk/en/persons/sotiria-fotopoulou>`_
* `Malcolm Bremer <https://research-information.bris.ac.uk/en/persons/malcolm-n-bremer>`_
* `Oliver Ray <https://research-information.bris.ac.uk/en/persons/oliver-ray>`_

.. _funding:

Funding
-------------------------

Grant Stevens acknowledges financial support from the `UKRI`_ for an `Interactive AI Centre for Doctoral Training`_ studentship at the
`University of Bristol`_.

.. image:: ../../images/CDT-UOB-logo.png
   :width: 38%
   :target: http://www.bristol.ac.uk/cdt/interactive-ai/

.. image:: ../../images/EPSRC+logo.png
   :width: 56%
   :target: https://gtr.ukri.org/projects?ref=studentship-2466020

.. raw:: html

 <hr>

Referencing the Package
-------------------------

Please see the :ref:`citing <citing>` page for instructions about referencing and citing
the AstronomicAL software.

.. _`University of Bristol`: https://www.bristol.ac.uk/
.. _`UKRI`: https://gtr.ukri.org/projects?ref=studentship-2466020
.. _`Interactive AI Centre for Doctoral Training`: http://www.bristol.ac.uk/cdt/interactive-ai/
.. title:: Citing AstronomicAL
.. _citing:


Please remember to cite our software and user guide whenever relevant.

Citing the Software
---------------------

.. image:: https://joss.theoj.org/papers/10.21105/joss.03635/status.svg
   :target: https://doi.org/10.21105/joss.03635

To reference the software and its functionality, consider citing the paper
published in JOSS_.

.. code :: latex

  @article{Stevens2021,
    doi = {10.21105/joss.03635},
    url = {https://doi.org/10.21105/joss.03635},
    year = {2021},
    publisher = {The Open Journal},
    volume = {6},
    number = {65},
    pages = {3635},
    author = {Grant Stevens and Sotiria Fotopoulou and Malcolm N. Bremer and Oliver Ray},
    title = {AstronomicAL: an interactive dashboard for visualisation, integration and classification of data with Active Learning},
    journal = {Journal of Open Source Software}
  }

.. _JOSS: https://joss.theoj.org/papers/10.21105/joss.03635


The paper is also available through arXiv_.

.. _arXiv: https://arxiv.org/abs/2109.05207
astronomicAL.utils
======================================

.. automodule:: astronomicAL.utils.load_config
   :members: verify_import_config, update_config_settings, create_layout_from_file, create_default_layout

.. automodule:: astronomicAL.utils.optimise
   :members: optimise_floats, optimise_ints, optimise_objects, optimise

.. automodule:: astronomicAL.utils.save_config
   :members: save_config_file
astronomicAL.extensions
======================================

.. autoclass:: astronomicAL.extensions.extension_plots.CustomPlot
   :members: create_settings, render, plot

.. automodule:: astronomicAL.extensions.extension_plots
  :members: get_plot_dict, create_plot, bpt_plot, mateos_2012_wedge

.. automodule:: astronomicAL.extensions.feature_generation
   :members: get_oper_dict

.. automodule:: astronomicAL.extensions.models
   :members: get_classifiers

.. automodule:: astronomicAL.extensions.query_strategies
   :members: get_strategy_dict
astronomicAL.active_learning
======================================

.. autoclass:: astronomicAL.active_learning.active_learning.ActiveLearningModel
   :members: assign_global_data, remove_from_pool, save_model, show_queried_point, iterate_AL, query_new_point, split_x_y_ids, exclude_unclassified_labels, train_val_test_split, reconstruct_tailored_sets, scale_data, split_y_ids, create_pool, setup_learners, generate_features, setup_panel, panel
astronomicAL.dashboard
======================================

.. autoclass:: astronomicAL.dashboard.active_learning.ActiveLearningDashboard
   :members: add_active_learning, panel

.. autoclass:: astronomicAL.dashboard.dashboard.Dashboard
   :members: set_contents, panel

.. autoclass:: astronomicAL.dashboard.dashboard.LabellingDashboard
   :members: _update_variable_lists, plot, get_previous_labels, save_label, select_random_point, get_current_index_in_labelled_data, panel

.. autoclass:: astronomicAL.dashboard.menu.MenuDashboard
   :members: panel

.. autoclass:: astronomicAL.dashboard.plot.PlotDashboard
  :members: update_df,update_variable_lists, plot, panel

.. autoclass:: astronomicAL.dashboard.selected_source.SelectedSourceDashboard
  :members: empty_selected, check_required_column, panel

.. autoclass:: astronomicAL.dashboard.settings_dashboard.SettingsDashboard
  :members: create_pipeline, get_settings, panel
astronomicAL.settings
======================================

.. autoclass:: astronomicAL.settings.active_learning.ActiveLearningSettings
   :members: update_data, get_default_variables, get_df, is_complete, panel

.. autoclass:: astronomicAL.settings.data_selection.DataSelection
   :members: get_dataframe_from_fits_file, add_ra_dec_col, get_df, panel

.. autoclass:: astronomicAL.settings.param_assignment.ParameterAssignment
   :members: update_data, _update_labels_cb, update_colours, get_id_column, get_label_column, get_label_colours, get_label_strings, get_settings, is_complete, panel
.. _preparing-data:
Preparing Your Data For Use in AstronomicAL
============================================

AstronomicAL has been designed so that researchers from any discipline can take advantage of the benefits of active learning with any tabular dataset. For this reason, we have placed very few requirements on the pre-processing required for your dataset.

Dataset Columns
-----------------

AstronomicAL only requires two columns to be included in your dataset:
  - **ID Column**: This column must contain unique identifiers for each data point.
  - **Label Column**: This column must contain integer labels for each data point.

.. note::
    We put no naming requirements for these columns (or any of the other columns in your dataset) as we ask you to select these corresponding columns in the settings panel.

*Which columns should I include?*
##########################################

As well as the two columns above, you should include:
    - All columns that you want to input as features to your classifiers.
    - Any columns that you require for axes in domain-specific plots.
    - Any columns that feature key information that would improve your ability to classify data points.

.. note::

	Active learning removes the requirement for large amounts of labelled data; this opens up possibilities for working with unlabelled datasets that are much larger. Although AstronomicAL can render millions of points, such large datasets can be demanding for memory consumption. **It is recommended that any columns not explicitly used should not be included in your dataset**.

Dataset Labels
-----------------
Each data point within your **Label Column** should have an integer assigned that corresponds to a particular class label.

In our :code:`example_dataset.fits` file, we use the following labels:
    - :code:`0` is a Star
    - :code:`1` is a Galaxy
    - :code:`2` is a QSO

AstronomicAL interprets the label :code:`-1` as an **unknown label**, and so you should only assign this value to data points that you do not yet know the correct label for (see below for more details).

For the known labels, it does not matter which value you assign to which class, as long as you're consistent across your dataset.

.. note::
    To interpret these labels easier, we allow you to assign name aliases for these values in the settings. These names will then be shown to you in plots and any labelling buttons.

Unknown Labels
##########################################

Unknown labels are handled differently from the other labels in the following ways:

    - None of the unknown labels are passed to the performance metrics as the models would interpret :code:`-1` as just another label leading to incorrect classifications due to your model never predicting :code:`-1`.
    - As AstronomicAL allows you to relabel any queried point, unknown labels are only valuable when used in the training pool. Therefore all unlabelled points are always assigned to the training set.

*What if I don't have any labels?*
##########################################

One of the main benefits of active learning is that you don't need large amounts of labelled data to produce high accuracy classifiers. Due to AstronomicAL's labelling framework, all that's required is one instance of every label in your classification task for you to start the rest of your labelling in either :ref:`active learning mode <active-learning>` or :ref:`labelling mode <labelling>`.
Using Your Models Outside of AstronomicAL
=============================================

Loading your model
---------------------------

The other tutorials show how you can use AstronomicAL to create accurate and robust classifiers using significantly less training data than normally required. This tutorial will explain the steps that will allow you to use your trained classifiers outside of AstronomicAL.

After every iteration of active learning or whenever you press the :code:`checkpoint` button during training, your model is saved inside :code:`models/`. As the models are created using `Scikit Learn`_ it is a quick and straightforward process to reload the models and use them on any new data.

.. _`SciKit Learn`: https://scikit-learn.org/stable/

To load your model in another python file or jupyter notebook run the following code:

.. code-block:: python
  :linenos:

  from joblib import load
  clf = load('models/your_classifier_model.joblib')

You can then simple run :code:`clf.predict(some_new_data)` to generate new predictions. If you want the probability outputs of the predictions you can run :code:`clf.predict_proba(some_new_data)` instead.

.. note::

	You must use data that contains the same number of features that you trained with your model. This includes any feature combinations that you generated in the settings.

    If you do not use the same amount of features, you will receive the following error:

        **ValueError: X has # features, but Classifier is expecting # features as input.**

Loading a Committee
---------------------------

If you created a committee of classifiers during training, you would need to continue using the committee as a whole to get the same performance. When your committee is saved, rather than just saving as a single model, each individual model is saved and bundled into a directory within :code:`data/`.

By default, the committees in AstronomicAL return the class probabilities averaged across each learner (known as *consensus probabilities*) and then the maximum probability is chosen to assign the prediction. It is recommended to recreate this behaviour outside of AstronomicAL to ensure that the model performs as expected.

.. code-block:: python
    :linenos:

    from joblib import load
    import glob
    import numpy as np

    classifiers = []

    committee_dir = "models/committee_dir"

    classifiers = {} # Dictionary to hold all the individual classifiers in the committee
    for i, filename in enumerate(glob.iglob(f"{committee_dir}/*.joblib")):

        # This is only required if you are using a scalar during training
        if "SCALER" in filename:
            continue

        classifiers[f"{i}"] = load(filename)

    for i, clf in enumerate(classifiers):
        pred_clf = classifiers[clf].predict_proba(some_new_data) # get the probability output for each classifier
        if i == 0:
            predictions_total = pred_clf
        else:
            predictions_total = predictions_total + pred_clf

    predictions_avg = predictions_total / len(classifiers.keys()) # This is the averaged probabilities for each class
    predictions = np.argmax(predictions_avg, axis=1) # This is the final predictions made

Scaling New Data
---------------------------
During the :ref:`settings <settings>` assignment, you had the option to scale your features. If this option was selected, you must scale any new data with the same scaler. Whenever your model is saved, the scaler (if used) is saved along with your updated models. The code below explains how to apply this scaler to your new data.

.. code-block:: python

    scaler = load("models/your_classifier_model-SCALER.joblib")

    new_data_scaled = scaler.transform(some_new_data)
.. _labelling:
Creating a Labelled Test Set
========================================

Overview
----------
Pre-assigned labels in surveys are often unreliable due to ambiguous class boundaries and data-limited information about each source. This is especially true when dealing with rarer classes. AstronomicAL allows you to create your own test set by manually labelling the points whilst viewing all the required information to ensure that the labels are reliable and sufficient for proving your model is accurate and generalisable.

.. figure:: ../../images/labelling_mode_load.gif

Load Configuration
**********************************
It is quick and easy to load a previous configuration file for **Labelling Mode** and uses the same settings setup as in **Active Learning Mode**, allowing you to jump between them when necessary.

Search Space
**********************************
By default, AstronomicAL will choose a point at random within your dataset when selecting a new point to label. However, it is often the case that sources with certain information available are more desirable to have in your test set than those missing this information.

.. figure:: ../../images/labelling_mode_demo.gif

AstronomicAL allows you to restrict your search space to only sources that match your specified criteria. Each addition criteria is matched with an *AND* operator, and so only when a source matches all criteria will it be available.

To remove a chosen criterion, simply select it from the **Criterion to Remove** dropdown and press the **Remove Criterion** button.

Labelling Your Data
**********************************
Pressing the **New** button will select a random source from your selection criteria, updating all the plots you have visible with this updated source. Once you are happy with the label you have selected, and press **Assign Label**. These labels are automatically saved to :code:`data/test_set.json`.

Using Your New Test Set
**********************************
Once you have labelled a sufficient number of sources, you can use your test set by checking the **Should data/test_set.json be the used test set?** box in the Active Learning Settings panel.

If you have already begun training a classifier and want to keep your progress, simply edit the configuration file you are using and set the :code:`"test_set_file"` flag to :code:`true`.

.. figure:: ../../images/use_test_set.png

If you don't set this flag, the code will automatically generate a test set from your data using the labels in your :code:`default label` column.
.. _custom-plots:
Creating Custom Plots
========================================

Overview
----------
To assign accurate labels to sources, the user must be presented with as full of a picture about each source as possible. Therefore the ability to create specialised plots is paramount.

Simple
**********************************
astronomicAL allows for fast plot creation using the :code:`create_plot` function, which abstracts away the potentially complicated plotting code typically required and replaces it with a simple one-line call.

Optimised
**********************************
By default, all custom plots have Datashader_ implemented, allowing millions of data points to be plotted at once whilst remaining responsive.

.. _Datashader: http://holoviews.org/user_guide/Large_Data.html

Collaborative
**********************************
Allowing researchers to share their designs is pivotal for collaborative research. For this reason, we make it seamless to share your designs with others by automatically creating a settings page to allow new users to select which columns from their own datasets correspond with the naming convention of the author.

Custom Plots
--------------------------
Custom plots can be added as a new function to :doc:`astronomicAL.extensions.extension_plots <../apireference/extensions>`.

There are some requirements when declaring a new feature generation function:

1. The new function must have 2 input parameters:
  - :code:`data` - The dataframe containing the entire dataset.
  - :code:`selected` - The currently selected points (Default:None)

2. The function must return the following:
  - :code:`plot` - The final plot to be rendered

3. The created function must be added in the form :code:`CustomPlot(new_plot_function, list_of_columns_used)` as a new value to the :code:`plot_dict` dictionary within the :code:`get_plot_dict` function, with a brief string key identifying the plot.

.. note::
    When using :code:`create_plot` within your custom plot functions, whenever you make use of a column name from your dataset, you must reference it as a key from the :code:`config.settings` dictionary.

    For example: :code:`"x_axis_feature"` should always be written as :code:`config.settings["x_axis_feature"]`.

    If this is not done, you will likely cause a :code:`KeyError` for any external researchers who want to use your plot but have a different dataset with different column names.


Example: Plotting :math:`Y=X^2`
-----------------------------------
In this example, we will show the simple case of creating a custom plot which shows :math:`Y=X^2` plotted along with the data.

.. code-block:: python
  :linenos:

  def x_squared(data, selected=None): # The function must include the parameters df, selected=None

      # This plots all the data on your chosen axes
      plot_data = create_plot(
          data,
          config.settings["x_axis"], # always use dataset column names as keys to config.settings
          config.settings["y_axis"],
      )

      line_x = np.linspace(
          np.min(data[config.settings["x_axis"]]),
          np.max(data[config.settings["x_axis"]]),
          30,
      ) # define the range of values for x

      line_y = np.square(line_x) # y=x^2

      line_data = pd.DataFrame(
          np.array([line_x, line_y]).T,
          columns=["x", "y"]
      ) # create dataframe with the newly created data

      # This plots the X^2 line
      line = create_plot(
          line_data,
          "x","y", # config.settings is not required here as these column names are not in reference to the main dataset
          plot_type="line", # we want a line drawn
          legend=False, # we don't need a legend for this data
          colours=False # Use default colours
      )

      x_squared_plot = plot_data * line # The * symbol combines multiple plots onto the same figure

      return x_squared_plot # The function must return the plot that is going to be rendered

Finally, add the new entry in the :code:`plot_dict` dictionary, **without specifying the parameters of the plotting function**:

.. code-block:: python

  def get_plot_dict():

      plot_dict = {
          "Mateos 2012 Wedge": CustomPlot(
              mateos_2012_wedge, ["Log10(W3_Flux/W2_Flux)", "Log10(W2_Flux/W1_Flux)"]
          ),
          "BPT Plots": CustomPlot(
              bpt_plot,
              [
                  "Log10(NII_6584_FLUX/H_ALPHA_FLUX)",
                  "Log10(SII_6717_FLUX/H_ALPHA_FLUX)",
                  "Log10(OI_6300_FLUX/H_ALPHA_FLUX)",
                  "Log10(OIII_5007_FLUX/H_BETA_FLUX)",
              ],
          ),
          "X^2": CustomPlot(
              x_squared,
              ["x_axis", "y_axis"],
          ),
      }

      return plot_dict

And that is all that is required. The new :code:`x_squared` plot is now available to use in AstronomicAL:

.. image:: ../../images/x_squared_in_plot_list.png

A settings page has automatically been generated, allowing users to select which of their dataset columns correspond to the author's specified column.

.. image:: ../../images/x_squared_settings.png

Once the columns have been chosen, the user is presented with the brand new :code:`x_squared` plot:

.. image:: ../../images/x_squared_example.png

Optional Plot Flags
-------------------

The :code:`create_plot` function allows users to specify the number of flags to ensure that the plot is as informative as possible.

The following pairs of images are arranged so that the *Flag=On* is on the left and *Flag=Off* on the right.

Colours
********************

.. image:: ../../images/fig_flags_colours.png
    :width: 47%
.. image:: ../../images/fig_flags_coloursN.png
    :width: 47%

The :code:`colours` flag will assign the colours the user specified in the opening settings. By choosing :code:`False`, all points remain the default colour set by Datashader.

The default for this value is :code:`True`.

.. raw:: html

   <hr>

Legends
*******************

.. image:: ../../images/fig_flags_legend.png
    :width: 47%
.. image:: ../../images/fig_flags_legendN.png
    :width: 47%

The :code:`legend` flag will include the plotted points in the plot legend. If all plots have this flag set to :code:`False` then no legend will be rendered.

The default for this value is :code:`True`.

.. raw:: html

   <hr>

Legend Positions
*******************

.. image:: ../../images/fig_flags_legend_position_inside.png
   :width: 47%
.. image:: ../../images/fig_flags_legend_position_outside.png
   :width: 47%

The :code:`legend_position` option allows you to position the legend in a more suitable place than the default positioning.

To keep the legend within the plot window you can choose between the following options: :code:`["top_left","top_right","bottom_left","bottom_right"]`.

To position the legend outside of the plot window, you can use one of the following options: :code:`["top","bottom","left","right"]`.

The examples above show :code:`["bottom_right"]` and :code:`["left"]` positions.

.. note::
    If all plots have the :code:`legend` flag set to :code:`False`, then the :code:`legend_position` flag is ignored, and no legend is rendered.

The default for this value is :code:`None`.

.. raw:: html

  <hr>

Smaller Axes Limits
*************************

.. image:: ../../images/fig_flags_smalleraxes.png
    :width: 47%
.. image:: ../../images/fig_flags_smalleraxesN.png
    :width: 47%

The :code:`smaller_axes_limits` flag will reduce the x and y axes limits so that the default ranges are between 4 standard deviations of the mean values. This can be used to reduce the negative impact on viewing from large outliers in the data, as can be seen above. However, all the data remains and values outside this range can still be viewed by interacting with the plot. If the minimum or maximum of an axis is already within four standard deviations of the mean, this will remain the limit for that axis.

.. note::

	If a selected source falls outside the range of the new axes limits, the axes ranges will extend to show the user that selected point so that the user does not miss out on potentially vital information when labelling.

The default for this value is :code:`None`, and so no axes limits are changed.

.. raw:: html

   <hr>

Bounded Axes
**********************

.. image:: ../../images/fig_flags_bounded.png
    :width: 47%
.. image:: ../../images/fig_flags_boundedN.png
    :width: 47%


The :code:`bounds` parameter, much like :code:`smaller_axes_limits`, will reduce the x and y axes limits; however, it does this much more abruptly, and any data points not within the specified bounds will be removed from the plot completely. The bound is specified as follows :code:`[xmin,ymax,xmax,ymin]` using the *[left,top,right,bottom]* style.

In the example above, we have assigned :code:`bounds=[0,1,1,0]` and as you can see below, if you zoom out, there are no points rendered outside this region.

.. image:: ../../images/fig_flags_bounded_1.png
    :width: 47%
    :align: center

This parameter is useful when you have missing data that default to extreme values, allowing you to specify the region representing realistic values.

If a selected source falls outside this region and is not shown on the plot, you can use this to indicate that the data for the chosen axes are not available for that data point.

The default for this value is :code:`None`, and so no axes limits are changed.

Slow Render
*******************

.. image:: ../../images/fig_flags_slow_render.png
    :width: 47%
.. image:: ../../images/fig_flags_legendN.png
    :width: 47%

The :code:`slow_render` flag removes all optimisations applied by the Datashader_ library and renders points using solely Bokeh. These points provide the user with much more customisability regarding glyph shapes and styles (see `Holoviews documentation`_ for more details).

.. _`Holoviews documentation`: http://holoviews.org/user_guide/Plotting_with_Bokeh.html

.. caution::
    Rendering points without Datashader requires substantially more processing power. As such, if you are rendering more than a few tens of thousands of points, you may notice the plots become laggy and unresponsive.

    It is recommended that this is only used when you have only a small sample of points that you want to emphasise in your plot.

    An example of this is when we render selected or queried points.

The default for this value is :code:`False`.

.. raw:: html

   <hr>
Reloading a Previous Configuration
========================================
AstronomicAL makes it easy to load a previous configuration, allowing you to continue training your model from an earlier checkpoint or to verify the results of someone else's classifier. Your configuration is automatically saved at each active learning iteration or by clicking the :code:`Save Current Configuration` button. It keeps track of your entire layout, settings, models and assigned labels and stores them within a small and easily sharable JSON file.

.. figure:: ../../images/Load_config_AL.gif

Below we can see the :code:`config/example_config.json` provided in the repository:

.. code-block:: JSON

    {
        "Author": "",
        "doi": "",
        "dataset_filepath": "data/example_dataset.fits",
        "optimise_data": true,
        "layout": {
            "0": {
                "x": 0,
                "y": 0,
                "w": 6,
                "h": 4,
                "contents": "Active Learning"
            },
            "1": {
                "x": 6,
                "y": 0,
                "w": 6,
                "h": 4,
                "contents": "Selected Source Info"
            },
            "2": {
                "x": 0,
                "y": 4,
                "w": 3,
                "h": 3,
                "contents": "Basic Plot",
                "panel_contents": [
                    "g-j",
                    "y-w1"
                ]
            },
            "3": {
                "x": 3,
                "y": 4,
                "w": 3,
                "h": 3,
                "contents": "SED Plot"
            },
            "4": {
                "x": 6,
                "y": 4,
                "w": 3,
                "h": 3,
                "contents": "BPT Plots"
            },
            "5": {
                "x": 9,
                "y": 4,
                "w": 3,
                "h": 3,
                "contents": "Mateos 2012 Wedge"
            }
        },
        "id_col": "ID",
        "label_col": "labels",
        "default_vars": [
            "g-j",
            "y-w1"
        ],
        "labels": [
            0,
            1,
            2
        ],
        "label_colours": {
            "0": "#ff7f0e",
            "1": "#2ca02c",
            "2": "#d62728",
            "-1": "#1f77b4"
        },
        "labels_to_strings": {
            "0": "Star",
            "1": "Galaxy",
            "2": "QSO",
            "-1": "Unknown"
        },
        "strings_to_labels": {
            "Unknown": -1,
            "Star": 0,
            "Galaxy": 1,
            "QSO": 2
        },
        "extra_info_cols": [
            "Lx",
            "redshift"
        ],
        "extra_image_cols": [
            "png_path_DR16"
        ],
        "labels_to_train": [
            "Star",
            "Galaxy",
            "QSO"
        ],
        "features_for_training": [
            "u",
            "g",
            "r",
            "i",
            "z",
            "y",
            "j",
            "h",
            "k",
            "w1",
            "w2"
        ],
        "exclude_labels": false,
        "exclude_unknown_labels": true,
        "unclassified_labels": [],
        "scale_data": false,
        "feature_generation": [
            [
                "subtract (a-b)",
                2
            ]
        ],
        "test_set_file": false,
        "classifiers": {
            "0": {
                "classifier": [
                    "RForest"
                ],
                "query": [
                    "Uncertainty Sampling"
                ],
                "id": [
                    "VIPERS 123044938",
                    "1032524471874906112",
                    "1031339472848971776",
                    "0185-00395",
                    "659811189332666368",
                    "0211-00169",
                    "2880063850911131648",
                    "VIPERS 123100497",
                    "0198-02511",
                    "3259593556628629504",
                    "4536487763679834112",
                    "3640086982559899648",
                    "601389188310394880",
                    "8391483576654950400",
                    "376114373996865536",
                    "VVDS-J022540.39-042204.3",
                    "326564606740817920",
                    "8390233431867060224",
                    "3723381311881191424",
                    "316540078866327552",
                    "3294552529806845952",
                    "0121-00286",
                    "3722278500913225728",
                    "4536478967586811904",
                    "336789235715565568",
                    "334496204081620992",
                    "1034778745740748800",
                    "1035856267169523712",
                    "VIPERS 124056894",
                    "8145996712583557120"
                ],
                "y": [
                    1,
                    1,
                    0,
                    1,
                    1,
                    1,
                    0,
                    1,
                    1,
                    0,
                    2,
                    0,
                    1,
                    0,
                    1,
                    0,
                    1,
                    0,
                    0,
                    -1,
                    0,
                    0,
                    0,
                    0,
                    -1,
                    1,
                    1,
                    0,
                    0,
                    0
                ]
            }
        }
    }


As you can see, all the chosen parameters from the :ref:`settings <settings>` panel are included in this document, which enables AstronomicAL to recreate the entire dashboard exactly how it was when the configuration was saved. This includes classifiers and the data they have been trained on, ensuring total reproducibility.

Editing The Configuration File
--------------------------------------
The configuration is entirely editable, and some users may find this a more convenient way of amending their settings.

Editing Parameters
##################################
Many changes, such as :code:`labels_to_strings` and :code:`label_colours` will only have cosmetic effects on the dashboard, making it simple to create plots that can be perfect for publications.

Changing parameters such as :code:`scale_data` or :code:`features_for_training` will affect the model specifically, so your model's performance may be drastically different. However, this can be useful to allow for quick prototyping of model parameters.

Editing Layout
##################################
AstronomicAL allows you to rearrange any of the panels to create a dashboard that works best with your workflow. Due to a limit in the React layout we currently use, which enables movable and expandable panels, it does not currently allow for dynamic adding and removing of panels. Therefore if you would like to have extra plots on your dashboard, then editing the :code:`layout` inside your configuration is the best way to do it.

Each panel is represented as follows:

.. code-block:: JSON

    {
        "5": {
            "x": 9,
            "y": 4,
            "w": 3,
            "h": 3,
            "contents": "Mateos 2012 Wedge"
        }

    }

The different parameters mean the following:
    - Each panel will have a string ID
    - :code:`x`: This is the x coordinate of the top left corner of the panel
    - :code:`y`: This is the y coordinate of the top left corner of the panel
    - :code:`w`: This is the width of the panel
    - :code:`h`: This is the height of the panel
    - :code:`contents`: The plot that will be displayed when configuration loaded

To add a new panel to your layout, you will need to assign the panel to an :code:`x` and :code:`y` coordinate that isn't already used. There is a width limit of 12 and so :code:`x + w` <= 12. There is no limit on the :code:`y` coordinate, but you will need to scroll down the page as you increase this value.

The :code:`contents` parameter is the type of plot displayed when you load the configuration file. These can be any of the custom plots defined in :ref:`this tutorial <custom-plots>`. If you do not know in advance what plot you want to view or you do not have the code for a particular domain-specific plot, then you can assign it as :code:`Menu`.
.. _active-learning:
Training a Classifier
=======================================================

Creating an Active Learning-based Classifier
---------------------------------------------

.. image:: ../../images/active_learning_start_panel.png

The Active Learning Dashboard is arranged in two sets of tabs:

  1. A Tab for each label classifier you chose in the previous settings will encapsulate everything you need to know about that particular one-vs-rest classifier.

  2. Inside each classifier tab, you have a set of tabs showing different plots related to how the classifier performs.

In this example, we will only be using the :code:`Star` tab as we are only training the Star classifier, but the same steps will apply for each classifier you choose to train.

.. raw:: html

   <hr>

Choosing your model
****************************************

For each classifier tab, you can assign which classifiers you want to use for the active learning process.

.. image:: ../../images/training_tutorial_AL_1.png

.. raw:: html

   <hr>

Choosing a Query strategy
**************************************

The main idea of Active Learning is that instead of piling as much data as possible onto a model to train on, you can get equal or better performance with substantially less training data if you analytically choose the most informative data points according to a metric. The query strategy *is* that metric.

Each classifier that you select is paired with your chosen query strategy.

In this run, we are using a Random Forest classifier with an Uncertainty Sampling query strategy.

.. image:: ../../images/training_tutorial_AL_4.png

.. raw:: html

   <hr>

Creating a Committee
*****************************
Even though we have only used a single classifier in this example, you are not restricted to only one. You can use any number of classifiers for your model, leading to an ensemble of classifiers known as a committee in Active Learning.

If you choose to create a committee, each classifier will have to retrain at each iteration of Active Learning, increasing the waiting times between queries.

When using a committee, when the model is saved, rather than being a single file for the classifier, it is saved as a folder of classifier files, which would need to continue being used together as an ensemble.

.. note::

	When adding multiple classifiers, you will still add a different query strategy for each; however, these are not used during training. Instead, the query strategy becomes the *vote entropy*, where the most informative point is the one that has the most significant disagreement between classifiers.

.. raw:: html

   <hr>

How Many Initial Points?
***************************
.. image:: ../../images/training_tutorial_AL_2.png
  :align: center

The final decision is to choose how many initial points to start the training process with. These points are chosen randomly, so choosing a high number may negatively impact Active Learning effectiveness and may reduce the maximum possible performance of your model. However, you will likely reduce the time spent on training and labelling.

Choosing a low number has the benefit that at the point you stop training, the majority of the points accumulated in your model's training set will have been selected based on their informativeness. However, to get good performance, you will have to hand-label more points which could become time-consuming.

.. raw:: html

   <hr>

.. image:: ../../images/training_tutorial_AL_5.png

Here is the final setup for our Star classifier.

Let the training begin!
-------------------------------------

Immediately after confirming the classifier settings, your model will begin training on some randomly selected points; how many is determined by the number you set in the previous section.

The First Results
***************************

.. image:: ../../images/active_learning_initial_train.png

Once trained, you will be presented with a colourful plot showing all your training data, with the axes specified in your :code:`default_x_variable` and :code:`default_y_variable` parameters from :ref:`settings <choosing-default-axis-variables>`. All the green points are your model's correct predictions, and red is your incorrect predictions. The blue dots are the five randomly chosen initial points, and the yellow point is the most informative data point based on the chosen classifier and query strategy.

Good Results, Bad Predictions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
At first glance, it seems like the model is performing pretty well: nearly 80% accuracy in both training and validation sets using only 5 points! However, these results may be deceiving.

The split between Stars and non-Stars here is roughly 33:66. Due to this imbalance, the classifier may develop a preference for predicting non-Star. However, one of the benefits of active learning is that by training on smaller amounts of data, it becomes easier to avoid the adverse effects of imbalanced datasets.

If we look at the correct and incorrect predictions (green and red areas), we see that the points from bottom left to top right look almost all correct. However, this is a very dense area, meaning that they may be overwhelmed by correct predictions if there are incorrect predictions (or vice-versa). To check how many incorrect points are actually there, we can hide the correct points by toggling off the :code:`Show Correct` button.

.. image:: ../../images/active_learning_toggle_correct.png

After removing the correct points, it is much easier to see only a couple of incorrect points are in the centre region. It is even more apparent that the problem lies in the two *branches* appearing from the bottom right.

The Labelling Process
-------------------------------------


We will need to add some more labelled data for the model to train on to improve our results across all metrics. However, we know very little about the current queried point and cannot make a confident classification without more information about the source.

So let's get more information about the source.

Exploring each source
***************************

.. image:: ../../images/choose_plot_type.png
    :align: center

Throughout the UI, you will have already noticed several **Choose plot type** panels. This is where the user can see more about each source at one time than would typically be possible.

Selected Source Information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One of the available plots is the :code:`Selected Source Information`, which is more of a mini dashboard than a plot, but it allows us to see critical information about the selected source.

.. image:: ../../images/selected_source_info.png

As you can see, we now get the crucial information required to make a confident classification of the source.

By default, the source's optical and Radio images are pulled from the SDSS_ and FIRST_ cutout services, respectively. These are provided free to the user as all that is required is the source's location (RA and Dec columns). As long as that area of the sky has been sampled, the images will be provided. (If you do not have these columns or are not using an astronomical dataset, these images will not be shown)

.. _SDSS: http://skyserver.sdss.org/dr16/en/help/docs/api.aspx#imgcutout
.. _FIRST: https://third.ucllnl.org/cgi-bin/firstcutout

We also see the two columns we specified in the settings earlier and the ID of the datapoint.

Sometimes, however, even this information may not be enough, and that is where the other plots are extremely useful.

The Basic Plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


The basic plot allows you to choose any :code:`X` and :code:`Y` axes from all the columns in your original file, as well as the additional features you created earlier.

When you first load the basic plot, the axes displayed will be your specified :code:`default_x_variable` and :code:`default_y_variable`, along with the colours you chose at the beginning for each label. As these are the same axes displayed in the Active Learning panel, we can now take a more detailed look at where things are going wrong.

.. image:: ../../images/basic_plot_small_stretched.png
  :width: 48%

.. image:: ../../images/toggle_correct_cropped.png
  :width: 46%

It is now much more apparent why we have the two branches of incorrect values. The branch trailing off the right are majoritively Stars, whereas the centre regions of majoritively Galaxies. The classifier is likely using the labels from the three trained on centre points (which will be labelled as non-Star as this is a one-vs-rest classifier) and labelling the Stars as non-Stars.

The branch on the left, which, as you approach the top half of the plot, are majoritively QSOs, is being classed as Stars. This is likely due to no QSOs being included in the classifier yet, leading it to view Stars as its closest match. Once a point is queried in that area and labelled as a non-Star, many red points will likely turn green.

.. raw:: html

   <hr>

Let's look at some of the other generated features and see if they can separate the data.

.. image:: ../../images/training_tutorial_AL_12.png
.. image:: ../../images/basic_plot_alternative_large.png
  :width: 70%

All plots are rendered using Bokeh_ and optimised using Datashader_, enabling you to plot millions of points at once whilst remaining responsive.

.. _Datashader: http://holoviews.org/user_guide/Large_Data.html
.. _Bokeh: https://docs.bokeh.org/en/latest/index.html

.. image:: ../../images/basic_plot_interactive.gif

.. raw:: html

   <hr>

Once again, we can see clear clustering between the sets of objects; however, the overlapping boundary is potentially still problematic. Taking all the information into consideration, we can assign this point the Galaxy label and continue training.

.. image:: ../../images/assigned_label.png

.. image:: ../../images/classifier_training.png

.. raw:: html

   <hr>

The model has now been retrained with the addition of the extra Galaxy, and now a new *most informative* point has been queried.

.. image:: ../../images/updated_queried_point.png

As you can see, the left branch of incorrect points has been largely reduced, and so our accuracy has increased up to 86% for both training and validation.

Analysing the Performance Scores
-------------------------------------

Performance Metrics
***************************

If you look at the score for this iteration, you can see that although accuracy, precision and F1-score increased, recall dropped by nearly 0.1. Is this something we should worry about?

Let's first begin with the definition for each of the metrics:

.. math::

    Precision = \frac{TP}{TP+FP}

.. math::

    Recall = \frac{TP}{TP+FN}

.. math::

    F1 Score = \frac{2 * Precision * Recall}{Precision+Recall}

Where TP is True Positives, FP is False Positives, and FN is False Negatives.

If we look at the change in precision, an increase of nearly 0.25 shows that the classifier is better at labelling only Stars as Stars. Unfortunately, the drop in recall shows that we are now misclassifying more Stars than before. Pair these with the increase in accuracy, which indicates that we are predicting more points to be non-Stars overall, and due to the higher number of non-Stars, this leads us to predict more points correctly. This is confirmed when we view the confusion matrices where the bottom right (TP) has decreased, whereas the sum of the left-hand columns (Predicted 0) has increased by around 2000 points.

.. image:: ../../images/conf_mat_train_orig.png
  :width: 60%
  :align: center

.. image:: ../../images/conf_mat_train_updated.png
  :width: 60%
  :align: center

.. raw:: html

   <hr>

Checking Training Progress
********************************

Now is a good time to look at the rest of the plots available to us within the Active Learning panel.

Training Set Plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: ../../images/training_set_plot.png

As we have seen already, this plot shows all the data within our training set, plotted according to whether our current model has predicted correctly. We also see which points the model has so far trained on and also the queried point, which would provide the most information to the model if it was labelled.

.. note::

	It is easy to get confused by the difference between the **training set** and the **points the model has trained on**.

  To clarify:

    **training set** = **training pool** + **points the model has trained on**

  Where the **training pool** are all the points the model gets to choose from when querying its next source.

.. raw:: html

   <hr>


Metric Plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: ../../images/metric_plot.png

Arguably the most interesting of the plots to look at is the metric plot, which is a visualisation of the query strategy and the driver for the active learning process. As we are using **Uncertainty Sampling**, this visualises the certainty the model has in its prediction. Green means the model is very confident in its prediction; Red means it's very unsure and can't decide whether the source is a Star or Galaxy.

.. caution::

	It is important to note that it will, at times, look as though the areas of high uncertainty match the areas of incorrect predictions from the model. However, with the query strategies we are using, the Active Learning query process completely ignores which label the model assigns to a source and therefore is not affected by correctness.

  It is easy to misunderstand this as *Active Learning improves your model's accuracy* when all it is doing is reducing the uncertainty of the most uncertain point at that particular iteration. It just so happens that, for many cases, the accuracy and other performance scores increase as a byproduct.

.. raw:: html

   <hr>


Validation Set Plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Next, we have the validation set plot, which is plotted according to its correctness, just like the training set plot.

The plot looks less densely packed because it is only 20% of the original dataset, whereas the training set is 60% of the original dataset.

.. image:: ../../images/val_set_plot.png

.. raw:: html

   <hr>


Score Tracking Plot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: ../../images/scores_plot.png

The score tracking plot does exactly that - tracks scores. This is useful for seeing the overall trend of your models improvement. As is commonly the case, you may start to notice your scores make smaller and smaller gains as more labelled data are added to the model, eventually leading to a nearly flat line over multiple iterations. Although there aren't specific stopping criteria for active learning, having your scores converge in this way with no changes in performance as you add more data might be a good time to stop.


.. raw:: html

   <hr>


It's OK to be Unsure
-------------------------------------

As you query more points, there will inevitably be a time when you are presented with an inconclusive point. This may be caused by specific features giving conflicting results or just that a particular source is missing too much key information for you to assign a confident and justifiable label.

Given that the model is likely to be training on such a small amount of data, it is not worth risking a potential incorrect label that *may* dramatically affect our models' performance.


.. image:: ../../images/assign_unsure.png

Labelling a point as unsure removes this point from the training set and then re-queries the training pool for the following most informative source.

No harm done!

Seeing the Results
-----------------------------

Training a little further (up to 20 points), let's see how our Star classifier has performed.

.. image:: ../../images/after_20_points_score.png
  :width: 47%

.. image:: ../../images/after_20_points_train.png
  :width: 49%

As you can see, the performance overall continues to improve. There are occasional drops, likely due to a queried point being in a part of the search space that has yet to be explored and causing local points to change labels abruptly; however, they bounce back almost immediately.

Saving your model
----------------------------

Now that the model has reached a suitable performance for us to apply it to new and unseen data, it is important that we save it for reusability and portability.

Well, the good news is that after each iteration of active learning, AstronomicAL automatically saves a copy of your model inside the :code:`models/` directory in the form :code:`label-Classifier_QueryStrategy.joblib`. This gets overwritten at each iteration, so it is always the most up-to-date. However, when you require something more permanent, you can use the :code:`Checkpoint` button.

.. image:: ../../images/training_tutorial_AL_28.png

This can be pressed once per iteration and will save your current model in the form :code:`label-Classifier_QueryStrategy-iteration-validationF1score-YYYYMMDD_HH:MM:SS.joblib`
to allow you to choose your best performing or most recent model quickly.

What About The Other Classifiers?
----------------------------------

In this example, we only used the Star classifier; well, what about the Galaxy classifier?

.. image:: ../../images/galaxy_al_panel.png

As you can see, each classifier tab is independent of the others, allowing you to tailor each classifier for each label. The workflow for training multiple classifiers is down to preference. You could focus on a single classifier until you are happy with its performance, then move on to the next. Alternatively, you could assign a label for a source on one classifier, switch tabs and label a source on one of the other classifiers. Each will produce the same results.

.. raw:: html

   <hr>

.. image:: ../../images/currently_selected.png

.. image:: ../../images/currently_not_selected.png

If you lose track of which tab the selected source is from, it is always shown at the bottom of each classifier tab whether the selected point is that classifier's queried point. If not, you can simply press the :code:`Show Queried` button to reselect the current classifier's queried point.
.. _custom_features:


Adding Custom Features to Your Model
====================================================

The features given to the model can often be the deciding factor for how well a model can produce accurate predictions. This is arguably even more so when approaching the problem using a method such as Active Learning, where you may only be using a tiny fraction of your entire dataset.

Custom Features
---------------------------------------------------
Custom generated features can be added as a new function to :doc:`astronomicAL.extensions.feature_generation <../apireference/extensions>`.

There are some requirements when declaring a new feature generation function:

1. The new function must have two input parameters:
  - :code:`df` - The DataFrame containing the entire dataset.
  - :code:`n` - The number of features involved in the operation.

  .. note::
      If your particular operation does not easily scale to more than two features at a time, then you can simply not make use of the :code:`n` input parameter inside your function.
      **However you must still include** :code:`n` **as an input parameter, even if you don't use it.**

2. The function must return the following:
  - :code:`df` - The updated dataframe with the newly generated features.
  - :code:`generated_features` - a list containing the updated column names.

3. The created function must be added as a new value to the :code:`oper` dictionary within the :code:`get_oper_dict` function, with a brief string key identifying the operation.

Within the function, you can generate the combinations of features using:

.. code-block:: python
  :linenos:

  base_features = config.settings["features_for_training"]
  combs = list(combinations(base_features, n))

Creating Colours
********************************************
Given the prevalence of photometry data in astronomical datasets, the most common additional features to create are colours. In AstronomicAL, these are provided with the default `subtract (a-b)` with a combination value of 2.


Example: Max(a, b)
-----------------------------------
In this example, we will show how we would create a new :code:`max` function. Although the produced features from this specific function may not be particularly useful for improving this model's performance, it works well as an example.

.. code-block:: python
  :linenos:

  def max_oper(df, n): # The function must include the parameters df and n

      np.random.seed(0) # set random seed if required for reproducability

      base_features = config.settings["features_for_training"] # get the base features chosen by the user

      combs = list(combinations(base_features, n)) # a list of all the combinations of n base_features

      cols = list(df.columns) # all the columns in the dataset
      generated_features = [] # The list that will keep track of all the new feature names

      for comb in combs: #loop over all combination tuples

          # This loop is to create the feature name string for the dataframe
          new_feature_name = "max(" # start of feature name
          for i in range(n): # loop over each feature in tuple
              new_feature_name = new_feature_name + f"{comb[i]}" # add each feature in operation
              if i != (n - 1):
                  new_feature_name = new_feature_name + "," # seperate features by a comma
              else:
                  new_feature_name = new_feature_name + ")"

          generated_features.append(new_feature_name) # add new feature name which is the form: max(f_1,f_2,...,f_n)


          if new_feature_name not in cols: # if the feature already exists in the data, dont recalculate

              # This loop applies the operation over all the feature in the combination and adds it as the new column in the dataframe
              for i in range(n): # Loop of each individual feature in comb
                  if i == 0:
                      df[new_feature_name] = df[comb[i]] # add the new column and set its value to the starting feature (without this you will get a KeyError)
                  else:
                      df[new_feature_name] = np.maximum(df[new_feature_name], df[comb[i]]) #calculate the running maximum

      return df, generated_features  # The function must return the updated dataframe and the list of generated features

Finally, adding the new entry in the :code:`oper` dictionary, **without specifying the parameters**:

.. code-block:: python

  def get_oper_dict():

      oper = {
          "subtract (a-b)": subtract,
          "add (a+b)": add,
          "multiply (a*b)": multiply,
          "divide (a/b)": divide,
          "max(a,b)": max_oper, # Newly created function
      }

      return oper

And that is all that is required. The new :code:`max_oper` function is now available to use in AstronomicAL:

.. image:: ../../images/create_feature_comb_list_max.png
.. _settings:
Completing the Data Setup
-----------------------------------

.. image:: ../../images/settings_completed.gif

Choosing your mode
***********************

When you first open AstronomicAL, you are given the option to enter Labelling mode or Active Learning mode. This choice determines the main screen shown after assigning settings. Both modes will allow you to assign all your required settings.

Once the settings have been completed, you can save your choices to a configuration file using the **Save Current Configuration** button at the top of the dashboard. You can then restart the software and load the configuration file in the other mode when required.

Loading your dataset
***********************

The first step is to load in your dataset. All files within the :code:`data/` directory will be selectable. If you have files in a different directory, creating a symlink/shortcut to your file inside :code:`data/` will suffice.

To begin with, we are going to load in the :code:`example_dataset.fits` file.

.. image:: ../../images/load_dataset.png
  :width: 50%
  :align: center

.. note::
  The **Optimise for memory** option allows the user to reduce the overall memory consumption of the dataset; however, this will increase the time it takes to load in the file and delay the user from progressing to the next stage. The more columns you have in your dataset, the more significant the overall improvement on memory; however, this also results in the most extended loading times.

Assigning Parameters
**********************

To make AstronomicAL as accessible as possible, there are very few requirements on columns. Where there are requirements, the user can specify the corresponding columns in the following dropdowns.

.. image:: ../../images/id_label_assignment.png
  :align: center

When the user chooses their labels column, the UI will autofill extra options for each unique label.

.. image:: ../../images/assign_params.png
  :align: center

Each label can have a unique colour assigned and a custom string identifier, all of which will be used throughout the UI.

.. image:: ../../images/colours_aliases.png
  :align: center

For this classification task, we have assigned Unknown, Star, Galaxy and QSO labels.

Extra Information Columns
#############################

Finally, we can specify which column values we would like to be displayed when inspecting a queried source. These extra bits of information can be especially useful when you have missing data that can't be trained on accurately but would improve the confidence in a classification if a source has the information available. Although this extra information remains unseen by the model, you can inject this information into the training process by labelling the data based on all the available information rather than just the features the model sees.

.. image:: ../../images/extra_columns.png
  :align: center

Here we want to see :code:`redshift` as the distance of an object can be very informative for its classification, as well as :code:`Lx` (X-ray luminosity) a feature only around 10% of our data has but is often the deciding factor on an object's classification. When a data point has this value available, we must see it.

.. image:: ../../images/extra_image_columns.png
  :align: center

Images can provide invaluable information that is often left out when only relying on trainable features. Each image column should contain either web addresses or local paths to the corresponding image for each data point (when available).

To summarise, our final parameter assignment settings look like the following:

.. image:: ../../images/assign_param_confirmed.png

Active Learning Settings
*************************

At first sight, the Active Learning settings panel looks a little overwhelming, but we will go through each part to ensure you can tune your models precisely as needed.

.. image:: ../../images/active_learning_settings.png

.. raw:: html

   <hr>

Selecting the Classifiers
##################################

The first step is to decide which classifiers we require to train. Currently, AstronomicAL produces a separate one-vs-rest classifier for each label as we have found this often produces more accurate predictions than a single multiclass classifier.

.. image:: ../../images/choose_labels_0.png
  :align: center
.. image:: ../../images/choose_labels_1.png
  :align: center
.. image:: ../../images/choose_labels_2.png
  :align: center

.. raw:: html

   <hr>

We are also given four options, all of which could have significant implications on each models performance.

.. image:: ../../images/al_settings_checkboxes.png


Should Remaining Labels be Removed from Active Learning Datasets?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you choose not to create a classifier for a particular label (by leaving it in the left column above), this option decides whether data with the left-out label should be removed from the training, validation, and test sets.

.. note::

  Selecting this option will only remove the data that gets inputted into the Active Learning process. The full dataset, including plots, will still contain all the data regardless of whether this option has been selected or not.

Should Features be scaled?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In many cases scaling your data can be effective at improving model performance. All data given to the Active Learning models will be scaled according to the training set when this option is selected.

The system handles this scaling; however, **if selected, the user must scale any new data they want predictions from according to the original training data**. For this reason, during the training process, AstronomicAL will save the scaler produced alongside the current model.

Should :code:`data/test_set.json` be used as the test set?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This allows you to use your hand-labelled points from :ref:`labelling mode <labelling>` as a test set. Such a test set ensures that you have reliable ground truth for which to test your model. If your model performs well on this set, you can be confident in its robustness and its generalisability to future data.

When selected, all data points included in this set are assigned as the test set; all remaining data points are then shared amongst the training and test sets.

When this option is not selected, the data is split into training, validation and test sets randomly (stratified sampling is used to ensure each set has an equal share of labels). As discussed in :ref:`Preparing your data <preparing-data>`, all unknown labels are also assigned to the training set.

As we have yet to use labelling mode in this example, we haven't created a curated test set, so this option remains disabled.

.. raw:: html

  <hr>

Selecting Your features
#################################################

.. image:: ../../images/classifier_features_selection.png

Given that your fits files will likely contain many more columns than those you require for training, you must select which columns will become your features in your model.

In this example, our base features will be photometry bands :code:`u`-:code:`W2`.


.. note::

	If you want to train on features that are combinations of each other, for example, when creating colours with photometry data, you don't need to include them in your fits file. The only features you need to include are base features that cannot be created from the combination of any other features.

.. raw:: html

   <hr>

Creating Feature Combinations
#####################################

The next step is to create any of the feature combinations we require. By default, AstronomicAL allows you to add, subtract, multiply and divide any :code:`n` features. For instructions on how to create custom feature combination functions, see :ref:`here<custom_features>`.

.. image:: ../../images/feature_combinations.png

.. caution::

	To find all the combinations of :math:`r` features of out all your baseline features of size :math:`n`, the following equation is used:
  .. math::

    \frac{n!}{r!(n-r)!}

  This quickly results in a huge number of additional features as :math:`r` and :math:`n` get larger. Please bear this in mind when adding these features, as this can increase training times substantially and have a negative impact on the performance of the dashboard.

.. note::

	Even though subtraction and division are not commutative or associative, we thought it was useful to the user to have the option to apply these operations to more than two features, especially when :ref:`creating a custom feature generation function<custom_features>`.


For all the combinations you add, which are displayed to you on the right-hand side of the row, all of the produced features will be available in both the Active Learning data as well as being plottable in the basic plot.

In this run, we have generated the colours from the photometric bands we chose earlier.

.. _choosing-default-axis-variables:
Choosing Default Axis Variables
#######################################

Next, we must assign the default x and y-axis variables. The columns chosen will become our default axes for many of the Active Learning plots, as well as being the opening axes in the :code:`basic plot`.

.. image:: ../../images/default_variables.png
  :align: center

The list of axes available automatically updates to include any of the generated features from your feature combination list. Here we have chosen :code:`g-j` and :code:`y-w1` as our default axes, both of which were generated from the subtract operation we selected above.

.. raw:: html

And this brings us to the end of the settings panel. We are now presented with a close button that will initialise either the labelling panel or active learning panel depending on your chosen mode at the beginning.

.. raw:: html

   <hr>

.. image:: ../../images/training_tutorial_settings_assign_params_20.png
