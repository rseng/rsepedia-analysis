---
title: 'ddop: A python package for data-driven operations management'

tags:
  - Python
  - scikit-learn
  - data-driven operations management
  - newsvendor

authors:
  - name: Andreas Philippi
    orcid: 0000-0002-6508-9128
    affiliation: 1
  - name: Simone Buttler
    orcid: 0000-0003-3986-057X
    affiliation: 1
  - name: Nikolai Stein
    orcid: 0000-0001-9847-3444
    affiliation: 1

affiliations:
 - name: Chair of Logistics and Quantitative Methods, Julius-Maximilians-Universität Würzburg,  
         Sandering 2, Würzburg 97070, Germany
   index: 1

date: 18 June 2021

bibliography: paper.bib

---

# Summary
In today's fast-paced world, companies face considerable uncertainty when making important decisions in operations management, for example, when deciding upon capacity, inventory levels, transportation, and production schedules. However, with the rise of digitization, companies have gained unprecedented access to data related to their particular decision problem, offering the opportunity to reduce the degree of uncertainty. For example, in inventory management the decision maker may have access to historical demand data as well as additional side information, such as social media data, customer behaviour, weather forecasts or calendar data. Driven by the availability of such rich data sources there has recently emerged a stream of literature in operations management research called “data-driven operations management” (DDOM). The focus of DDOM is to combine machine learning and traditional optimization techniques to prescribe cost optimal decisions directly from data. Various models have been developed and shown great performance on the dataset used. However, what is missing is efficient access to open-source code and datasets.
With *ddop*, we provide a Python library that integrates well-established algorithms form the field of data-driven operations management, as well as standard benchmark datasets. Thus, *ddop* helps researchers in two ways:

* Researchers can efficiently apply and compare well-established DDOM models.
* Researchers can test new developed models on benchmark datasets provided in the package. 

The application programming interface (API) of *ddop* is designed to be consistent, easy-to-use, and accessible even for non-experts. With only a few lines of code, one can build and compare various models. In *ddop* all models are offered as objects implementing the estimator interface from scikit-learn [@buitinck2013api]. We thus not only provide a uniform API for our models, but also ensure that they safely interact with scikit-learn pipelines, model evaluation and selection tools.  

The library is distributed under the 3-Clause BSD license, encouraging its use in both academic and commercial settings. The full source code is available at https://github.com/opimwue/ddop. The package can be installed via the Python Package Index using `pip install ddop`.  A detailed documentation providing all information required to work with the API can be found at https://opimwue.github.io/ddop/. 

# Statement of need

With the growing number of publications in the field of data-driven operations management, comparability is becoming increasingly difficult. The reasons for this are twofold: One, most scientists work with proprietary company data which cannot be shared. Two, it is not yet standard that researchers share code used to implement their models. Consequently, results are not directly reproducible and models have to be re-implemented every time a researcher wants to benchmark a new approach. This not only takes a lot of time but can also be a demanding process since such complex models are often challenging to implement. Against this background, there has recently been a call to take inspiration from the machine learning community, where great APIs like scikit-learn [@buitinck2013api], fastai [@howard2020fastai], or Hugging Face [@wolf2019huggingface] have been developed that allow previous developed ML models to be effectively applied on different dataset. Following up on this, *ddop* is the first of its kind to integrate well-established data-driven models for operations management tasks. At the current state, this includes various approaches to solve the data-driven newsvendor problem, such as weighted sample average approximation [@bertsimas2020predictive], empirical risk minimization [@ban2019big], and a deep learning based approach [@oroojlooyjadid2020applying]. In addition, the library provides different real-world datasets that can be used to quickly illustrate the behaviour of the available models or as a benchmark for testing new models. *ddop's* aim is to make data-driven operations management accessible and reproducible. 


# Usage
Since all models in *ddop* implement the estimator interface from *scikit-learn* consisting of a *fit*, *predict*, and *score* method, usage follows the standard procedure of an *scikit-learn* regressor. First, a model is initialized by calling the class constructor from a given set of constant hyper-parameter values, each describing the model or the optimisation problem the estimator tries to solve. Note that for ease of use, all estimators use reasonable default values. It is therefore not necessary to pass any parameter to the constructor. However, it is recommended to tune them for the respective application, since this can often improve decision quality. After the model has been initialized, the *fit* method is used to learn a decision model from the training data (*X_train*, *y_train*). Once the training process is completed, the function returns the fitted model, which can then be used to make decisions for new data (*X_test*) by using the *predict* method. Finally, the score method can be used to access the decision quality of a model. The method takes as input *X_test* as well as the corresponding true values *y_test* and computes the average costs between *y_test* and *predict(X_test)*. Because all estimators follow the same interface, using a different model is as simple as replacing the constructor.

# Future Work
There are several directions that the ddop project aims to focus on in future development. While at the current state there are only algorithms available to solve the newsvendor problem, the goal is to include models to solve other operations management task like multi-period inventory management or capacity management. In addition, we aim to extend the library in terms of available datasets and tutorials.

# References 
.. -*- mode: rst -*-

.. image:: https://travis-ci.com/opimwue/ddop.svg?branch=master
    :target: https://travis-ci.com/opimwue/ddop

.. image:: https://d25lcipzij17d.cloudfront.net/badge.svg?id=py&type=6&v=0.6.5&x2=0
    :target: https://badge.fury.io/py/ddop

.. image:: https://img.shields.io/github/license/andreasphilippi/ddop
    :target: https://github.com/andreasphilippi/ddop/blob/master/LICENSE
    
.. image:: https://www.code-inspector.com/project/22456/status/svg
    :target: https://frontend.code-inspector.com/public/project/22456/ddop/dashboard
    
.. image:: https://joss.theoj.org/papers/0de119f95840b69fcea94309c18058e4/status.svg
    :target: https://joss.theoj.org/papers/0de119f95840b69fcea94309c18058e4   
    

----------------------


Welcome to ddop!
====================

.. image:: /docsrc/logos/logo.png
    :width: 300

``ddop`` is a Python library for data-driven operations management. The goal of ``ddop`` is to provide well-established
data-driven operations management tools within a programming environment that is accessible and easy to use even
for non-experts. At the current state ``ddop`` contains well known data-driven newsvendor models, a set of
performance metrics that can be used for model evaluation and selection, as well as datasets that are useful to
quickly illustrate the behavior of the various algorithms implemented in ``ddop`` or as benchmark for testing new
models. Through its consistent and easy-to-use interface one can run and compare provided models with only a few
lines of code.

------------------------------------------------------------

Installation
------------

ddop is available via PyPI using:

.. code-block:: bash

    pip install ddop

The installation requires the following dependencies:

- numpy==1.18.2
- scipy==1.4.1
- pandas==1.1.4
- statsmodels==0.11.1
- scikit-learn==0.23.0
- tensorflow==2.4.1
- pulp==2.0
- mpmath

Note: The package is actively developed and conflicts with other packages may occur during
installation. To avoid any installation conflicts we therefore recommend to install the
package in an empty environment with the above mentioned dependencies

Quickstart
----------
``ddop`` provides a varity of newsvendor models. The following example
shows how to use one of these models for decision making. It assumes
a very basic knowledge of data-driven operations management practices.

As first step we initialize the model we want to use. In this example
`LinearRegressionNewsvendor <https://opimwue.github.io/ddop/modules/auto_generated/ddop.newsvendor.LinearRegressionNewsvendor.html#ddop.newsvendor.LinearRegressionNewsvendor>`__.

.. code-block:: python

    >>> from ddop.newsvendor import LinearRegressionNewsvendor
    >>> mdl = LinearRegressionNewsvendor(cu=2,co=1)

A model can take a set of parameters, each describing the model or the optimization
problem it tries to solve. Here we set the underage costs ``cu`` to 2 and
the overage costs ``co`` to 1.

As next step we load the `Yaz Dataset <https://opimwue.github.io/ddop/modules/auto_generated/ddop.datasets.load_yaz.html#ddop.datasets.load_yaz>`__ and split it into train and test set.

.. code-block:: python

    >>> from ddop.datasets import load_yaz
    >>> from sklearn.model_selection import train_test_split
    >>> X, y = load_yaz(one_hot_encoding=True, return_X_y=True)
    >>> X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle=False, random_state=0)

After the model is initialized, the ``fit`` method can be used to learn a decision model from the training data ``X_train``, ``y_train``.

.. code-block:: python

    >>> mdl.fit(X_train, y_train)

We can then use the ``predict`` method to make a decision for new data samples.

.. code-block:: python

    >>> mdl.predict(X_test)
    >>> array([[ 8.32..,  7.34.., 16.92.., ..]])

To get a representation of the model's decision quality we can use the ``score`` function, which takes as input
``X_test`` and  ``y_test``. The score function makes a decision for each sample in ``X_test`` and calculates
the negated average costs with respect to the true values ``y_test`` and the overage and underage costs.

.. code-block:: python

    >>> mdl.score(X_test,y_test)
    -7.05..

------------------------------------------------------------

See also
-----------
* Follow the `API reference <https://opimwue.github.io/ddop/api_reference.html>`__ to get an overview of available functionalities and for detailed class and function information.
* To get familiar with ``ddop`` and to learn more about data-driven operations management check out our `Tutorials <https://opimwue.github.io/ddop/tutorial.html>`__.

------------------------------------------------------------
.. _api_reference:

=============
API Reference
=============

This is the class and function reference of ddop.

.. _newsvendor_ref:

:mod:`ddop.newsvendor`: Newsvendor decision making
===================================================

The ``ddop.newsvendor`` module contains different newsvendor approaches for
decision making.

.. automodule:: ddop.newsvendor
    :no-members:
    :no-inherited-members:

.. currentmodule:: ddop.newsvendor

Sample Average Approximation (SAA)
-----------------------------------

.. autosummary::
    :toctree: modules/auto_generated/
    :template: class.rst

    SampleAverageApproximationNewsvendor


Weighted SAA (wSAA)
----------------------

.. autosummary::
    :toctree: modules/auto_generated/
    :template: class.rst

    DecisionTreeWeightedNewsvendor
    RandomForestWeightedNewsvendor
    KNeighborsWeightedNewsvendor
    GaussianWeightedNewsvendor


Empirical Risk Minimization (ERM)
---------------------------------

.. autosummary::
    :toctree: modules/auto_generated/
    :template: class.rst

    LinearRegressionNewsvendor
    DeepLearningNewsvendor


------------------------------------------------------------

.. _metrics_ref:

:mod:`ddop.metrics`: Evaluation metrics
========================================

The ``ddop.metrics`` module includes different performance metrics that can be used for model selection and
evaluation.

.. automodule:: ddop.metrics
    :no-members:
    :no-inherited-members:

.. currentmodule:: ddop.metrics

.. autosummary::
    :toctree: modules/auto_generated/

    total_costs
    average_costs
    prescriptiveness_score

All performance metrics can also be used with scikit-learn model selection tools. However, therefore
a proper scoring object has to be generated by using the make_scorer function.

.. autosummary::
    :toctree: modules/auto_generated/

    make_scorer

Moreover, the module contains a function to calculate the pairwise costs.

.. autosummary::
    :toctree: modules/auto_generated/

    pairwise_costs


------------------------------------------------------------

.. _datasets_ref:

:mod:`ddop.datasets`: Datasets
==============================

``ddop`` comes with a few default datasets that can be loaded using the ``ddop.datasets`` module.


.. automodule:: ddop.datasets
    :no-members:
    :no-inherited-members:

Loaders
-----------

.. currentmodule:: ddop.datasets

.. autosummary::
    :toctree: modules/auto_generated/

    load_yaz
    load_bakery
    load_SID

These datasets are useful to quickly illustrate the behavior of the various algorithms implemented in ddop.Tutorial
==========

The following tutorial is designed especially for developers without prior knowledge of the domain covered by the API.
Within the tutorial you can run through a complete application scenario, from the installation to the comparison of
the various newsvendor models included in ddop. In this context, implementation details, as well as theoretical principles
are discussed.

.. toctree::
    :maxdepth: 1

    The Data-Driven Newsvendor <tutorial_modules/tutorial.ipynb>
.. ddop documentation master file, created by
   sphinx-quickstart on Fri Jul 17 13:31:59 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ddop!
=================

.. image:: ../docsrc/logos/logo.png
    :width: 300

``ddop`` is a Python library for data-driven operations management. The goal of ``ddop`` is to provide well-established
data-driven operations management tools within a programming environment that is accessible and easy to use even
for non-experts. At the current state ``ddop`` contains well known data-driven newsvendor models, a set of
performance metrics that can be used for model evaluation and selection, as well as datasets that are useful to
quickly illustrate the behavior of the various algorithms implemented in ``ddop`` or as benchmark for testing new
models. Through its consistent and easy-to-use interface one can run and compare provided models with only a few
lines of code.

------------------------------------------------------------

Installation
------------

ddop is available via PyPI using:

.. code-block:: bash

    pip install ddop

The installation requires the following dependencies:

- numpy==1.18.2
- scipy==1.4.1
- pandas==1.1.4
- statsmodels==0.11.1
- scikit-learn==0.23.0
- tensorflow==2.4.1
- pulp==2.0
- mpmath

Note: The package is actively developed and conflicts with other packages may occur during
installation. To avoid any installation conflicts we therefore recommend to install the
package in an empty environment with the above mentioned dependencies

Quickstart
----------
``ddop`` provides a varity of newsvendor models. The following example
shows how to use one of these models for decision making. It assumes
a very basic knowledge of data-driven operations management practices.

As first step we initialize the model we want to use. In this example
`LinearRegressionNewsvendor <https://andreasphilippi.github.io/ddop/modules/auto_generated/ddop.newsvendor.LinearRegressionNewsvendor.html#ddop.newsvendor.LinearRegressionNewsvendor>`__.

.. code-block:: python

    >>> from ddop.newsvendor import LinearRegressionNewsvendor
    >>> mdl = LinearRegressionNewsvendor(cu=2,co=1)

A model can take a set of parameters, each describing the model or the optimization
problem it tries to solve. Here we set the underage costs ``cu`` to 2 and
the overage costs ``co`` to 1.

As next step we load the `Yaz Dataset <https://andreasphilippi.github.io/ddop-kit/modules/auto_generated/ddop.datasets.load_yaz.html#ddop.datasets.load_yaz>`__ and split it into train and test set.

.. code-block:: python

    >>> from ddop.datasets import load_yaz
    >>> from sklearn.model_selection import train_test_split
    >>> X, y = load_yaz(one_hot_encoding=True, return_X_y=True)
    >>> X_train, X_test, y_train, y_test = train_test_split(X, y, shuffle=False, random_state=0)

After the model is initialized, the ``fit`` method can be used to learn a decision model from the training data ``X_train``, ``y_train``.

.. code-block:: python

    >>> mdl.fit(X_train, y_train)

We can then use the ``predict`` method to make a decision for new data samples.

.. code-block:: python

    >>> mdl.predict(X_test)
    >>> array([[ 8.32..,  7.34.., 16.92.., ..]])

To get a representation of the model's decision quality we can use the ``score`` function, which takes as input
``X_test`` and  ``y_test``. The score function makes a decision for each sample in ``X_test`` and calculates
the negated average costs with respect to the true values ``y_test`` and the overage and underage costs.

.. code-block:: python

    >>> mdl.score(X_test,y_test)
    -7.05..

------------------------------------------------------------

See also
-----------
Follow the `API reference <https://andreasphilippi.github.io/ddop/api_reference.html>`__ to get an overview of available functionalities and for detailed class and function information.
To get familiar with ``ddop`` and to learn more about data-driven operations management check out our `Tutorials <https://andreasphilippi.github.io/ddop/tutorial.html>`__.

.. toctree::
   :maxdepth: 1

   api_reference
   tutorial

:mod:`{{module}}`.{{objname}}
{{ underline }}==============

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}
    {% if methods %}
   .. rubric:: Methods

   .. autosummary::
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
.. _bakery_dataset:

Bakery dataset
----------------

The bakery dataset contains the demand for a number of products from different stores. Moreover, it stores a number of
demand features. A description of targets and features is given below.


**Dataset Characteristics:**

    :Number of Instances: 127575

    :Number of Targets: 1

    :Number of Features: 13

    :Target Information:
        - 'demand' the corresponding demand observation

    :Feature Information:
        - 'date' the date
        - 'weekday' the day of the week,
        - 'month' the month of the year,
        - 'year' the year,
        - 'is_holiday' whether or not it is a national holiday,
        - 'is_holiday_next2days' whether or not it is a national holiday in the next two days,
        - 'is_schoolholiday' whether or not it is a school holiday,
        - 'store' the store id,
        - 'product' the product id,
        - 'rain' the amount of rain,
        - 'temperature' the average temperature in °C,
        - 'promotion_currentweek' whether or not there is a promotion this week
        - 'promotion_lastweek' whether there was a promotion last week

    Note: By default the date feature is not included when loading the data. You can include it
    by setting the parameter `include_date` to `True`.



.. _yaz_dataset:

YAZ dataset
----------------

This is a real world dataset from Yaz. Yaz is a fast casual restaurant in Stuttgart providing good service
and food at short waiting times. The dataset contains the demand for the main ingredients at YAZ.
Moreover, it stores a number of demand features. These features include information about the day, month, year,
lag demand, weather conditions and more. A full description of targets and features is given below.


**Dataset Characteristics:**

    :Number of Instances: 765

    :Number of Targets: 7

    :Number of Features: 12

    :Target Information:
        - 'calamari' the demand for calamari
        - 'fish' the demand for fish
        - 'shrimp' the demand for shrimps
        - 'chicken' the demand for chicken
        - 'koefte' the demand for koefte
        - 'lamb' the demand for lamb
        - 'steak' the demand for steak

    :Feature Information:
        - 'date' the date,
        - 'weekday' the day of the week,
        - 'month' the month of the year,
        - 'year' the year,
        - 'is_holiday' whether or not it is a national holiday,
        - 'is_closed' whether or not the restaurant is closed,
        - 'weekend' whether or not it is weekend,
        - 'wind' the wind force,
        - 'clouds' the cloudiness degree,
        - 'rain' the amount of rain,
        - 'sunshine' the sunshine hours,
        - 'temperature' the outdoor temperature

    Note: By default the date feature is not included when loading the data. You can include it
    by setting the parameter `include_date` to `True`.





.. _SID_dataset:

SID dataset
-------------

This dataset contains 5 zears of store-item demand data

**Dataset Characteristics:**

    :Number of Instances: 887284

    :Number of Targets: 1

    :Number of Features: 6

    :Target Information:
       - 'demand' the corresponding demand observation

    :Feature Information:
        - 'date' the date
        - 'weekday' the day of the week,
        - 'month' the month of the year,
        - 'year' the year,
        - 'store' the store id,
        - 'item' the item id

    Note: By default the date feature is not included when loading the data. You can include it
    by setting the parameter `include_date` to `True`.





