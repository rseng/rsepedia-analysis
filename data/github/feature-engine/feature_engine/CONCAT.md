# Feature Engine

![PythonVersion](https://img.shields.io/badge/python-3.6%20|3.7%20|%203.8%20|%203.9-success)
[![License https://github.com/feature-engine/feature_engine/blob/master/LICENSE.md](https://img.shields.io/badge/license-BSD-success.svg)](https://github.com/feature-engine/feature_engine/blob/master/LICENSE.md)
[![PyPI version](https://badge.fury.io/py/feature-engine.svg)](https://badge.fury.io/py/feature-engine)
[![Conda https://anaconda.org/conda-forge/feature_engine](https://anaconda.org/conda-forge/feature_engine/badges/installer/conda.svg)](https://anaconda.org/conda-forge/feature_engine)
[![CircleCI https://app.circleci.com/pipelines/github/feature-engine/feature_engine?branch=1.1.X](https://img.shields.io/circleci/build/github/feature-engine/feature_engine)](https://app.circleci.com/pipelines/github/feature-engine/feature_engine?branch=1.1.X)
[![Documentation Status https://feature-engine.readthedocs.io/en/latest/index.html](https://readthedocs.org/projects/feature-engine/badge/?version=latest)](https://feature-engine.readthedocs.io/en/latest/index.html)
[![Join the chat at https://gitter.im/feature_engine/community](https://badges.gitter.im/feature_engine/community.svg)](https://gitter.im/feature_engine/community)
[![Sponsorship https://www.trainindata.com/](https://img.shields.io/badge/Powered%20By-TrainInData-orange.svg)](https://www.trainindata.com/)
[![Downloads](https://pepy.tech/badge/feature-engine)](https://pepy.tech/project/feature-engine)
[![Downloads](https://pepy.tech/badge/feature-engine/month)](https://pepy.tech/project/feature-engine)
[![DOI](https://zenodo.org/badge/163630824.svg)](https://zenodo.org/badge/latestdoi/163630824)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03642/status.svg)](https://doi.org/10.21105/joss.03642)


[<img src="./docs/images/logo/FeatureEngine.png" width="248">](http://feature-engine.readthedocs.io)

Feature-engine is a Python library with multiple transformers to engineer and select features for use in machine learning models. 
Feature-engine's transformers follow Scikit-learn's functionality with fit() and transform() methods to learn the 
transforming parameters from the data and then transform it.


## Feature-engine features in the following resources

* [Feature Engineering for Machine Learning, Online Course](https://www.udemy.com/course/feature-engineering-for-machine-learning/?referralCode=A855148E05283015CF06)

* [Feature Selection for Machine Learning, Online Course](https://www.udemy.com/course/feature-selection-for-machine-learning/?referralCode=186501DF5D93F48C4F71)

* [Deployment of Machine Learning Models, Online Course](https://www.udemy.com/course/deployment-of-machine-learning-models/?referralCode=D4FE5EA129FFD203CFF4)

* [Python Feature Engineering Cookbook](https://packt.link/python)


## Blogs about Feature-engine

* [Feature-engine: A new open-source Python package for feature engineering](https://trainindata.medium.com/feature-engine-a-new-open-source-python-package-for-feature-engineering-29a0ab88ea7c)

* [Practical Code Implementations of Feature Engineering for Machine Learning with Python](https://towardsdatascience.com/practical-code-implementations-of-feature-engineering-for-machine-learning-with-python-f13b953d4bcd)


## En Español

* [Ingeniería de variables para machine learning, Curso Online](https://www.udemy.com/course/ingenieria-de-variables-para-machine-learning/?referralCode=CE398C784F17BD87482C)

* [Ingeniería de variables, MachinLenin, charla online](https://www.youtube.com/watch?v=NhCxOOoFXds)


## Documentation

* [Documentation](http://feature-engine.readthedocs.io)


## Current Feature-engine's transformers include functionality for:

* Missing Data Imputation
* Categorical Encoding
* Discretisation
* Outlier Capping or Removal
* Variable Transformation
* Variable Creation
* Variable Selection
* Datetime Feature Extraction
* Preprocessing
* Scikit-learn Wrappers

### Imputation Methods
* MeanMedianImputer
* RandomSampleImputer
* EndTailImputer
* AddMissingIndicator
* CategoricalImputer
* ArbitraryNumberImputer
* DropMissingData

### Encoding Methods
* OneHotEncoder
* OrdinalEncoder
* CountFrequencyEncoder
* MeanEncoder
* WoEEncoder
* PRatioEncoder
* RareLabelEncoder
* DecisionTreeEncoder

### Discretisation methods
* EqualFrequencyDiscretiser
* EqualWidthDiscretiser
* DecisionTreeDiscretiser
* ArbitraryDiscreriser

### Outlier Handling methods
* Winsorizer
* ArbitraryOutlierCapper
* OutlierTrimmer

### Variable Transformation methods
* LogTransformer
* LogCpTransformer
* ReciprocalTransformer
* PowerTransformer
* BoxCoxTransformer
* YeoJohnsonTransformer

### Variable Creation:
 * MathematicalCombination
 * CombineWithReferenceFeature
 * CyclicalTransformer

### Feature Selection:
 * DropFeatures
 * DropConstantFeatures
 * DropDuplicateFeatures
 * DropCorrelatedFeatures
 * SmartCorrelationSelection
 * ShuffleFeaturesSelector
 * SelectBySingleFeaturePerformance
 * SelectByTargetMeanPerformance
 * RecursiveFeatureElimination
 * RecursiveFeatureAddition
 * DropHighPSIFeatures

### Preprocessing
 * MatchVariables
 
### Wrappers:
 * SklearnTransformerWrapper

## Installation

From PyPI using pip:

```
pip install feature_engine
```

From Anaconda:

```
conda install -c conda-forge feature_engine
```

Or simply clone it:

```
git clone https://github.com/feature-engine/feature_engine.git
```

## Example Usage

```python
>>> import pandas as pd
>>> from feature_engine.encoding import RareLabelEncoder

>>> data = {'var_A': ['A'] * 10 + ['B'] * 10 + ['C'] * 2 + ['D'] * 1}
>>> data = pd.DataFrame(data)
>>> data['var_A'].value_counts()
```

```
Out[1]:
A    10
B    10
C     2
D     1
Name: var_A, dtype: int64
```
    
```python 
>>> rare_encoder = RareLabelEncoder(tol=0.10, n_categories=3)
>>> data_encoded = rare_encoder.fit_transform(data)
>>> data_encoded['var_A'].value_counts()
```

```
Out[2]:
A       10
B       10
Rare     3
Name: var_A, dtype: int64
```

Find more examples in our [Jupyter Notebook Gallery](https://nbviewer.org/github/feature-engine/feature-engine-examples/tree/main/) 
or in the [documentation](http://feature-engine.readthedocs.io).

## Contribute

Details about how to contribute can be found in the [Contribute Page](https://feature-engine.readthedocs.io/en/latest/contribute/index.html)

Briefly:

- Fork the repo
- Clone your fork into your local computer: ``git clone https://github.com/<YOURUSERNAME>/feature_engine.git``
- navigate into the repo folder ``cd feature_engine``
- Install Feature-engine as a developer: ``pip install -e .``
- Optional: Create and activate a virtual environment with any tool of choice
- Install Feature-engine dependencies: ``pip install -r requirements.txt`` and ``pip install -r test_requirements.txt``
- Create a feature branch with a meaningful name for your feature: ``git checkout -b myfeaturebranch``
- Develop your feature, tests and documentation
- Make sure the tests pass
- Make a PR

Thank you!!


### Documentation

Feature-engine documentation is built using [Sphinx](https://www.sphinx-doc.org) and is hosted on [Read the Docs](https://readthedocs.org/).

To build the documentation make sure you have the dependencies installed: from the root directory: ``pip install -r docs/requirements.txt``.

Now you can build the docs using: ``sphinx-build -b html docs build``


## License

BSD 3-Clause

## Donate

[Sponsor the maintainer](https://github.com/sponsors/solegalli) to support her continue expanding 
Feature-engine.---
title: 'Feature-engine: A Python package for feature engineering for machine learning'
tags:
  - python
  - feature engineering
  - feature selection
  - machine learning
  - data science
authors:
  - name: Soledad Galli
    affiliation: 1
affiliations:
 - name: Train in Data
   index: 1
date: 6 August 2021
bibliography: paper.bib
---

# Summary

Feature-engine is an open source Python library with the most exhaustive battery of 
transformations to engineer and select features for use in machine learning. Feature-engine 
supports several techniques to impute missing data, encode categorical variables, transform 
variables mathematically, perform discretization, remove or censor outliers, and combine 
variables into new features. Feature-engine also hosts an array of algorithms for feature 
selection.

The primary goal of Feature-engine is to make commonly used data transformation procedures 
accessible to researchers and data scientists, focusing on creating user-friendly and 
intuitive classes, compatible with existing machine learning libraries, like Scikit-learn 
[@sklearn] and Pandas [@pandas].

Many feature transformation techniques learn parameters from data, like the values for 
imputation or the mappings for encoding. Feature-engine classes learn these parameters 
from the data and store them in their attributes to transform future data. Feature-engine’s 
transformers preserve Scikit-learn’s functionality with the methods fit() and transform() 
to learn parameters from and then transform data. Feature-engine's transformers can be 
incorporated into a Scikit-learn Pipeline to streamline data transformation and facilitate 
model deployment, by allowing the serialization of the entire pipeline in one pickle.

When pre-processing a dataset different feature transformations are applied to different 
variable groups. Feature-engine classes allow the user to select which variables to transform 
within each class, therefore, while taking the entire dataframe as input, only the indicated 
variables are modified. Data pre-processing and feature engineering are commonly done 
together with data exploration. Feature-engine transformers return dataframes as output, 
thus, users can continue to leverage the power of Pandas for data analysis and visualization 
after transforming the data set.

In summary, Feature-engine supports a large variety of commonly used data transformation 
techniques [@data_prep; @boxcox; @yeojohnson; @kdd_2009_competition; 
@beatingkaggle; @micci_mean_encoder], as well as techniques that were developed 
in data science competitions [@niculescu09_kdd], including those for feature selection 
[@miller09_kdd]. Thus, Feature-engine builds upon and extends the capabilities of 
Python's current scientific computing stack and makes accessible transformations that 
are otherwise not easy to find, understand or code, to data scientist and data 
practitioners.



# Statement of need

Data scientists spend an enormous amount of time on data pre-processing and transformation 
ahead of training machine learning models [@domingos]. While some feature engineering 
processes can be domain-specific, a large variety of transformations are commonly applied 
across datasets. For example, data scientists need to impute or remove missing values or 
transform categories into numbers, to train machine learning models using Scikit-learn, 
the main library for machine learning. Yet, depending on the nature of the variable and 
the characteristics of the machine learning model, they may need to use different techniques. 

Feature-engine gathers the most frequently used data pre-processing techniques, as well as 
bespoke techniques developed in data science competitions, in a library, from which users can pick 
and choose the transformation that they need, and use it just like they would use any other 
Scikit-learn class. As a result, users are spared of manually creating a lot of code, which 
is often repetitive, as the same procedures are applied to different datasets. In addition, 
Feature-engine classes are written to production standards, which ensures classes return 
the expected result, and maximizes reproducibility between research and production 
environments through version control.

In the last few years, a number of open source Python libraries that support feature 
engineering techniques have emerged, highlighting the importance of making feature 
engineering and creation accessible and, as much as possible, automated. Among these, 
Featuretools [@kanter2015deep] creates features from temporal and relational datasets, 
tsfresh [@christ_tsfresh] extracts features from time series, Category encoders 
[@category_encoders] supports a comprehensive list of methods to encode categorical 
variables, and Scikit-learn [@sklearn] implements a number of data transformation 
techniques, with the caveat that the transformations are applied to the entire dataset, 
and the output are NumPy arrays. Feature-engine extends the capabilities of the current 
Python’s scientific computing stack by allowing the application of the transformations 
to subsets of variables in the dataset, returning dataframes for data exploration, and 
supporting transformations not currently available in other libraries, like those for 
outlier censoring or removal, besides additional techniques for discretization and 
feature selection that were developed by data scientist working in the industry or data 
science competitions.


# Acknowledgements

I would like to acknowledge all of the contributors and users of Feature-engine, who helped 
with valuable feedback, bug fixes, and additional functionality to further improve the library. 
A special thanks to Christopher Samiullah for continuous support on code quality and 
architecture. A list of  Feature-engine contributors is available at 
https://github.com/feature-engine/feature_engine/graphs/contributors.

# References---
name: Docs
about: What documentation is missing?
title: ''
labels: ''
assignees: ''

---

Please let us know if you think there is information missing, or how else we can improve the documentation from Feature-engine.

If you are referring to an existing page, please paste the url.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: Jupyter notebook examples
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

Please let us know what is missing from existing Jupyter notebook demos, or suggest which new demo you think it would be useful for the community.
.. feature_engine documentation master file, created by
   sphinx-quickstart on Wed Jan 10 14:43:38 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Feature-engine
==============

A Python library for Feature Engineering and Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure::  images/logo/FeatureEngine.png
   :align:   center

   **Feature-engine rocks!**

Feature-engine is a Python library with multiple transformers to engineer and select
features to use in machine learning models. Feature-engine preserves Scikit-learn
functionality with methods `fit()` and `transform()` to learn parameters from and then
transform the data.

Feature-engine includes transformers for:

- Missing data imputation
- Categorical encoding
- Discretisation
- Outlier capping or removal
- Variable transformation
- Variable combination
- Variable selection
- Datetime features
- Preprocessing

Feature-engine allows you to select the variables you want to transform **within** each
transformer. This way, different engineering procedures can be easily applied to
different feature subsets.

Feature-engine transformers can be assembled within the Scikit-learn pipeline,
therefore making it possible to save and deploy one single object (.pkl) with the
entire machine learning pipeline. Check :ref:`**Quick Start** <quick_start>` for an
example.

What is unique about Feature-engine?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following characteristics make Feature-engine unique:

- Feature-engine contains the most exhaustive battery of feature engineering transformations.
- Feature-engine can transform a specific group of variables in the dataframe.
- Feature-engine returns dataframes, hence suitable for data exploration and model deployment.
- Feature-engine is compatible with the Scikit-learn pipeline.
- Feature-engine automatically recognizes numerical and categorical variables.
- Feature-engine alerts you if a transformation is not possible, e.g., if applying logarithm to negative variables or divisions by 0.

If you want to know more about what makes Feature-engine unique, check this
`article <https://trainindata.medium.com/feature-engine-a-new-open-source-python-package-for-feature-engineering-29a0ab88ea7c>`_.


Installation
~~~~~~~~~~~~

Feature-engine is a Python 3 package and works well with 3.6 or later. Earlier versions
have not been tested. The simplest way to install Feature-engine is from PyPI with pip:

.. code-block:: bash

    $ pip install feature-engine

Note, you can also install it with a _ as follows:

.. code-block:: bash

    $ pip install feature_engine

Feature-engine is an active project and routinely publishes new releases. To upgrade
Feature-engine to the latest version, use pip like this:

.. code-block:: bash

    $ pip install -U feature-engine

If you’re using Anaconda, you can install the
`Anaconda Feature-engine package <https://anaconda.org/conda-forge/feature_engine>`_:

.. code-block:: bash

    $ conda install -c conda-forge feature_engine


Feature-engine features in the following resources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `Feature Engineering for Machine Learning <https://www.udemy.com/course/feature-engineering-for-machine-learning/?referralCode=A855148E05283015CF06>`_, Online Course.
- `Feature Selection for Machine Learning <https://www.udemy.com/course/feature-selection-for-machine-learning/?referralCode=186501DF5D93F48C4F71>`_, Online Course.
- `Python Feature Engineering Cookbook <https://packt.link/python>`_.
- `Feature-engine: A new open-source Python package for feature engineering <https://trainindata.medium.com/feature-engine-a-new-open-source-python-package-for-feature-engineering-29a0ab88ea7c/>`_.
- `Practical Code Implementations of Feature Engineering for Machine Learning with Python <https://towardsdatascience.com/practical-code-implementations-of-feature-engineering-for-machine-learning-with-python-f13b953d4bcd>`_.

More learning resources in the :ref:`**Learning Resources** <learning_resources>`.


Feature-engine's Transformers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Feature-engine hosts the following groups of transformers:

Missing Data Imputation: Imputers
---------------------------------

- :doc:`api_doc/imputation/MeanMedianImputer`: replaces missing data in numerical variables by the mean or median
- :doc:`api_doc/imputation/ArbitraryNumberImputer`: replaces missing data in numerical variables by an arbitrary number
- :doc:`api_doc/imputation/EndTailImputer`: replaces missing data in numerical variables by numbers at the distribution tails
- :doc:`api_doc/imputation/CategoricalImputer`: replaces missing data with an arbitrary string or by the most frequent category
- :doc:`api_doc/imputation/RandomSampleImputer`: replaces missing data by random sampling observations from the variable
- :doc:`api_doc/imputation/AddMissingIndicator`: adds a binary missing indicator to flag observations with missing data
- :doc:`api_doc/imputation/DropMissingData`: removes observations (rows) containing missing values from dataframe

Categorical Encoders: Encoders
------------------------------

- :doc:`api_doc/encoding/OneHotEncoder`: performs one hot encoding, optional: of popular categories
- :doc:`api_doc/encoding/CountFrequencyEncoder`: replaces categories by the observation count or percentage
- :doc:`api_doc/encoding/OrdinalEncoder`: replaces categories by numbers arbitrarily or ordered by target
- :doc:`api_doc/encoding/MeanEncoder`: replaces categories by the target mean
- :doc:`api_doc/encoding/WoEEncoder`: replaces categories by the weight of evidence
- :doc:`api_doc/encoding/PRatioEncoder`: replaces categories by a ratio of probabilities
- :doc:`api_doc/encoding/DecisionTreeEncoder`: replaces categories by predictions of a decision tree
- :doc:`api_doc/encoding/RareLabelEncoder`: groups infrequent categories

Variable Discretisation: Discretisers
-------------------------------------

- :doc:`api_doc/discretisation/ArbitraryDiscretiser`: sorts variable into intervals defined by the user
- :doc:`api_doc/discretisation/EqualFrequencyDiscretiser`: sorts variable into equal frequency intervals
- :doc:`api_doc/discretisation/EqualWidthDiscretiser`: sorts variable into equal width intervals
- :doc:`api_doc/discretisation/DecisionTreeDiscretiser`: uses decision trees to create finite variables

Outlier Capping or Removal
--------------------------

-  :doc:`api_doc/outliers/ArbitraryOutlierCapper`: caps maximum and minimum values at user defined values
-  :doc:`api_doc/outliers/Winsorizer`: caps maximum or minimum values using statistical parameters
-  :doc:`api_doc/outliers/OutlierTrimmer`: removes outliers from the dataset

Numerical Transformation: Transformers
--------------------------------------

- :doc:`api_doc/transformation/LogTransformer`: performs logarithmic transformation of numerical variables
- :doc:`api_doc/transformation/LogCpTransformer`: performs logarithmic transformation after adding a constant value
- :doc:`api_doc/transformation/ReciprocalTransformer`: performs reciprocal transformation of numerical variables
- :doc:`api_doc/transformation/PowerTransformer`: performs power transformation of numerical variables
- :doc:`api_doc/transformation/BoxCoxTransformer`: performs Box-Cox transformation of numerical variables
- :doc:`api_doc/transformation/YeoJohnsonTransformer`: performs Yeo-Johnson transformation of numerical variables

Mathematical Combination:
-------------------------

-  :doc:`api_doc/creation/MathematicalCombination`: creates new variables by combining features with mathematical operations
-  :doc:`api_doc/creation/CombineWithReferenceFeature`: combines variables with reference features
-  :doc:`api_doc/creation/CyclicalTransformer`: creates variables using sine and cosine, suitable for cyclical features

Feature Selection:
------------------

- :doc:`api_doc/selection/DropFeatures`: drops an arbitrary subset of variables from a dataframe
- :doc:`api_doc/selection/DropConstantFeatures`: drops constant and quasi-constant variables from a dataframe
- :doc:`api_doc/selection/DropDuplicateFeatures`: drops duplicated variables from a dataframe
- :doc:`api_doc/selection/DropCorrelatedFeatures`: drops correlated variables from a dataframe
- :doc:`api_doc/selection/SmartCorrelatedSelection`: selects best features from correlated groups
- :doc:`api_doc/selection/DropHighPSIFeatures`: selects features based on the Population Stability Index (PSI)
- :doc:`api_doc/selection/SelectByShuffling`: selects features by evaluating model performance after feature shuffling
- :doc:`api_doc/selection/SelectBySingleFeaturePerformance`: selects features based on their performance on univariate estimators
- :doc:`api_doc/selection/SelectByTargetMeanPerformance`: selects features based on target mean encoding performance
- :doc:`api_doc/selection/RecursiveFeatureElimination`: selects features recursively, by evaluating model performance
- :doc:`api_doc/selection/RecursiveFeatureAddition`: selects features recursively, by evaluating model performance

Datetime:
---------
- :doc:`api_doc/datetime/DatetimeFeatures`: extract features from datetime variables

Preprocessing:
--------------

- :doc:`api_doc/preprocessing/MatchVariables`: ensures that columns in test set match those in train set

Scikit-learn Wrapper:
---------------------

-  :doc:`api_doc/wrappers/Wrapper`: applies Scikit-learn transformers to a selected subset of features


Getting Help
~~~~~~~~~~~~

Can't get something to work? Here are places where you can find help.

1. The :ref:`**User Guide** <user_guide>` in the docs.
2. `Stack Overflow <https://stackoverflow.com/search?q=feature_engine>`_. If you ask a question, please mention "feature_engine" in it.
3. If you are enrolled in the `Feature Engineering for Machine Learning course <https://www.udemy.com/course/feature-engineering-for-machine-learning/?referralCode=A855148E05283015CF06>`_ , post a question in a relevant section.
4. If you are enrolled in the `Feature Selection for Machine Learning course <https://www.udemy.com/course/feature-selection-for-machine-learning/?referralCode=186501DF5D93F48C4F71>`_ , post a question in a relevant section.
5. Join our `gitter community <https://gitter.im/feature_engine/community>`_. You an ask questions here as well.
6. Ask a question in the repo by filing an `issue <https://github.com/feature-engine/feature_engine/issues/>`_ (check before if there is already a similar issue created :) ).


Contributing
~~~~~~~~~~~~

Interested in contributing to Feature-engine? That is great news!

Feature-engine is a welcoming and inclusive project and we would be delighted to have you
on board. We follow the
`Python Software Foundation Code of Conduct <http://www.python.org/psf/codeofconduct/>`_.

Regardless of your skill level you can help us. We appreciate bug reports, user testing,
feature requests, bug fixes, addition of tests, product enhancements, and documentation
improvements. We also appreciate blogs about Feature-engine. If you happen to have one,
let us know!

For more details on how to contribute check the contributing page. Click on the
:ref:`**Contribute** <contribute>` guide.

Donate
~~~~~~

`Sponsor us <https://github.com/sponsors/solegalli>`_ to support us continue expanding
Feature-engine.

Open Source
~~~~~~~~~~~

Feature-engine's `license <https://github.com/feature-engine/feature_engine/blob/master/LICENSE.md>`_
is an open source BSD 3-Clause.

Feature-engine is hosted on `GitHub <https://github.com/feature-engine/feature_engine/>`_.
The `issues <https://github.com/feature-engine/feature_engine/issues/>`_ and
`pull requests <https://github.com/feature-engine/feature_engine/pulls>`_ are tracked there.


Table of Contents
~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 2

   quickstart/index
   user_guide/index
   api_doc/index
   resources/index
   contribute/index
   about/index
   whats_new/index
{{index}}
{{summary}}
{{extended_summary}}
{{parameters}}
{{returns}}
{{yields}}
{{other_parameters}}
{{attributes}}
{{raises}}
{{warns}}
{{warnings}}
{{see_also}}
{{notes}}
{{references}}
{{examples}}
{{methods}}{{objname}}
{{ underline }}==============

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}

   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
   {% for item in methods %}
      {% if '__init__' not in item %}
        ~{{ name }}.{{ item }}
      {% endif %}
   {%- endfor %}
   {% endif %}
   {% endblock %}

.. include:: {{module}}.{{objname}}.examples

.. raw:: html

    <div style='clear:both'></div>
.. _governance:

Governance
==========

The purpose of this document is to formalize the governance process used by the
Feature-engine project and clarify how decisions are made and how the community works
together. This is the first version of our governance policy and will be updated as our
community grows and more of us take on different roles.

Roles and Responsibilities
--------------------------

Contributors
~~~~~~~~~~~~

Contributors are community members who contribute in various ways to the project.
Anyone can become a contributor, and contributions can be of various forms, not just
code. To see how you can help check the :ref:`**Contribute page** <contribute>`.


Core Contributors
~~~~~~~~~~~~~~~~~

Core Contributors are community members who are dedicated to the continued development
of the project through ongoing engagement with the community. Core Contributors are
expected to review code contributions, can approve and merge pull requests, can decide
on the fate of pull requests, and can be involved in deciding major changes to the
Feature-engine API. Core Contributors determine who can join as a Core Contributor.


Founder and Leadership
~~~~~~~~~~~~~~~~~~~~~~

Feature-engine was founded by `Soledad Galli <https://www.linkedin.com/in/soledad-galli/>`_
who at the time was solely responsible for the initial prototypes, documentation, and
dissemination of the project. In the tradition of Python, Sole is referred to as the
“benevolent dictator for life” (BDFL) of the project, or simply, “the founder”. From a
governance perspective, the BDFL has a special role in that they provide vision,
thought leadership, and high-level direction for the project’s members and contributors.
The BDFL has the authority to make all final decisions for the Feature-engine Project.
However, in practice, the BDFL, chooses to defer that authority to the consensus of the
community discussion channels and the Core Contributors. The BDFL can also propose and
vote for new Core Contributors.


Join the community
------------------

Feature-engine is currently looking to expand the team of Core Contributors, if you are
interested, please get in touch.

If you want to Contribute to the project in any other way, get in touch using our Github
issues page or through Gitter:

1. `Github issues <https://github.com/feature-engine/feature_engine/issues/>`_.
2. `Gitter community <https://gitter.im/feature_engine/community>`_... raw :: html

    <!-- Generated by generate_authors_table.py -->
    <div class="authors-container">
        <style>
          img.avatar {border-radius: 15px; padding: 10px;}
          .author {text-align: center;}
        </style>
        <div class="author">
            <a href='https://github.com/solegalli'><img src='https://avatars.githubusercontent.com/solegalli?v=4' class='avatar' width="120" height="120" /></a> <br />
            <p>Soledad Galli</p>
        </div>
        <div class="author">
            <a href='https://github.com/christophergs'><img src='https://avatars.githubusercontent.com/christophergs?v=4' width="120" height="120"class='avatar' /></a> <br />
            <p>Chris Samiullah</p>
        </div>
        <div class="author">
            <a href='https://github.com/nicogalli'><img src='https://avatars.githubusercontent.com/nicogalli?v=4' class='avatar'width="120" height="120"/></a> <br />
            <p>Nicolas Galli</p>
        </div>
    </div>
Roadmap
=======

This document provides general directions on what the core contributors would like to
see developed in Feature-engine. As resources are limited, we can't promise when or if
the transformers listed here will be included in the library. We welcome all the help
we can get to support this vision. If you are interested in contributing, please get in
touch.

Purpose
-------

Feature-engine's mission is to simplify and streamline the implementation of end-to-end
feature engineering pipelines. It aims to help users both during the research phase and
while putting a model in production.

Feature-engine makes data engineering easy by allowing the selection of feature subsets
directly within its transformers. It also interlaces well with exploratory data analysis
(EDA) by returning dataframes for easy data exploration.

Feature-engine’s transformers preserve Scikit-learn functionality with the methods fit()
and transform() and can be integrated into a Pipeline to simplify putting the model in
production.

Feature-engine was designed to be used in real settings. Each transformer has a concrete
aim, and is tailored to certain variables and certain data. Transformers raise errors
and warnings to support the user to use a suitable transformation given the data.
These errors help avoid inadvertedly incorporating missing values to the dataframe at
unwanted stages of the development.


Vision
------

At the moment, Feature-engine's functionality is tailored to cross-sectional or tabular
data, mostly numerical or categorical. But we would like to extend its functionality
to work with datetime, text and time series. In the following figure we show how we
would like the overall structure of Feature-engine to look like:

.. figure::  ../images/FeatureEnginePackageStructure.png
   :align:   center

   Feature-engine structure

Current functionality
---------------------

Most of the functionality for cross-sectional data is already included in the package.
We expand and update this arm of the library, based on user feedback and suggestions
and our own research in the field. In grey, the transformers that are not yet included
in the package:

.. figure::  ../images/FeatureEnginePackageStructureCrossSectional.png
   :align:   center

   Transformers for cross-sectional data

The current transformations supported by Feature-engine return features that are easy
to interpret, and the effects of the transformations are clear and easy to understand.
The original aim of Feature-engine was to provide technology that is suitable to create
models that will be used in real settings, and return understandable variables.

Having said this, more and more, users are requesting features to combine or transform
variables in ways that would return features that are not human readable, in an attempt
to improve model performance and perhaps have an edge in data science competitions. We
are currently contemplating the incorporation of this functionality to the package.

Wanted functionality
--------------------

We are interested in adding a module that creates date and time related features from
datetime variables. This module would include transformers to extract all possible date
and time related features, like hr, min, sec, day, year, is_weekend, etc. And it would
also include transformers to capture elapsed time between 2 or more variables.

We would also like to add a module that returns straightforward features from simple
text variables, to capture text complexity, like for example counting the number
of words, unique words, lexical complexity, number of paragraphs and sentences. We would
also consider integrating the Bag of Words and TFiDF from sklearn with a wrapper that
returns a dataframe ready to use to train machine learning models. Below we show more
detail into these new modules.

.. figure::  ../images/FeatureEnginePackageStructureDatetimeText.png
   :align:   center

   New models wanted: datetime and text

In addition, we are evaluating whether including a module to extract features from time
series is possible, within the current design of the package, and if it adds real value
compared to the functionality already existing in pandas and Scipy, and in other well
established open source projects like `tsfresh <https://tsfresh.readthedocs.io/en/latest/>`_
and `featuretools <https://featuretools.alteryx.com/en/stable/index.html>`_.
The transformations we are considering are shown in this image:

.. figure::  ../images/FeatureEnginePackageStructureTimeseries.png
   :align:   center

   Time series module and the transformations envisioned


Goals
-----

Our main goals are:

- Continue maintaining a high-quality, well-documented collection of canonical tools for data processing.
- Expand the documentation with more examples about Feature-engine's functionality.
- Expand the documentation with more detail on how to contribute to the package.
- Expand the library's functionality as per the precedent paragraphs.

For more fine-grained goals and current and lined-up issues please visit the `issues <https://github.com/feature-engine/feature_engine/issues/>`_
section in our repo.

.. -*- mode: rst -*-
.. _about:

About
=====

In this section you will find information about the Feature-engine's origin, main
developers, roadmap and overall vision for the package. You will also find information
about how to cite Feature-engine and our main sponsors.

.. toctree::
   :maxdepth: 1

   about
   governance
   roadmap.. -*- mode: rst -*-

About
=====

History
-------

Data scientists spend a huge amount of time on data pre-processing and transformation.
It would be great (we thought back in the day) to gather the most frequently used data
pre-processing techniques and transformations in a library, from which we could pick
and choose the transformation that we need, and use it just like we would use any other
sklearn class. This was the original vision for Feature-engine.

Feature-engine is an open source Python package originally designed to support the online
course `Feature Engineering for Machine Learning in Udemy <https://www.udemy.com/feature-engineering-for-machine-learning/?couponCode=FEATENGREPO>`_,
but has now gained popularity and supports transformations beyond those taught in the
course. It was launched in 2017, and since then, several releases have appeared and a
growing international community is beginning to lead the development.

Governance
----------

The decision making process and governance structure of Feature-engine is laid out in
the :ref:`**governance document** <governance>`.

Core contributors
-----------------

The following people are currently core contributors to Feature-engine’s development
and maintenance:

.. include:: authors.rst

Contributors
------------

A growing international community is beginning to lead Feature-engine's development.
You can learn more about Feature-engine's Contributors in the
`GitHub contributors page <https://github.com/feature-engine/feature_engine/graphs/contributors>`_.

Citing Feature-engine
---------------------

.. image:: https://zenodo.org/badge/163630824.svg
   :target: https://zenodo.org/badge/latestdoi/163630824

.. image:: https://joss.theoj.org/papers/10.21105/joss.03642/status.svg
   :target: https://joss.theoj.org/papers/10.21105/joss.03642

|

If you use Feature-engine in a scientific publication, you can cite the following paper:
Galli, S., (2021). `Feature-engine: A Python package for feature engineering for machine learning. <https://joss.theoj.org/papers/10.21105/joss.03642>`_
Journal of Open Source Software, 6(65), 3642.

Bibtex entry:

.. code-block:: bibtex

    @article{Galli2021,
    doi = {10.21105/joss.03642},
    url = {https://doi.org/10.21105/joss.03642},
    year = {2021},
    publisher = {The Open Journal},
    volume = {6},
    number = {65},
    pages = {3642},
    author = {Soledad Galli},
    title = {Feature-engine: A Python package for feature engineering for machine learning},
    journal = {Journal of Open Source Software}
    }



You can also find a DOI (digital object identifier) for every version of Feature-engine
on `zenodo.org <https://zenodo.org/badge/latestdoi/163630824>`_; use the BibTeX on this
site to reference specific versions of the software.


Artwork
-------

High quality PNG and SVG logos are available in the `docs/images/ <https://github.com/feature-engine/feature_engine/tree/main/docs/images/logo>`_
source directory of the repository.

.. figure::  ../images/logo/FeatureEngine.png
   :width: 200
   :figclass: align-center
   :align: center


Sponsors
--------

Feature-engine is a community driven project, however institutional and private grants
help to assure its sustainability. The project would like to thank the following
sponsors:

|
|

.. raw:: html

   </div>
   </div>

........

.. raw:: html

   <div class="sk-sponsor-div">
   <div class="sk-sponsor-div-box">

Soledad Galli spends a big part of her time at `Train in Data <https://www.trainindata.com/>`_
maintaining the project.

.. raw:: html

   </div>

   <div class="sk-sponsor-div-box">

.. image:: ../images/sponsors/trainindata.png
   :width: 150pt
   :align: center
   :target:  https://www.trainindata.com/

.. raw:: html

   </div>
   </div>

..........

|
|.. -*- mode: rst -*-
.. _user_guide:

User Guide
==========

In this section you will find additional information about Feature-engine's transformers
and feature engineering transformations in general, as well as additional examples.

.. toctree::
   :maxdepth: 1

   imputation/index
   encoding/index
   discretisation/index
   outliers/index
   transformation/index
   creation/index
   selection/index
   datetime/index
   preprocessing/index
   wrappers/index
.. _math_combination:

.. currentmodule:: feature_engine.creation

MathematicalCombination
=======================

:class:`MathematicalCombination()` applies basic mathematical operations to multiple
features, returning one or more additional features as a result. That is, it sums,
multiplies, takes the average, finds the maximum, minimum or standard deviation of a
group of variables and returns the result into new variables.

For example, if we have the variables:

- **number_payments_first_quarter**,
- **number_payments_second_quarter**,
- **number_payments_third_quarter** and
- **number_payments_fourth_quarter**,

we can use :class:`MathematicalCombination()` to calculate the total number of payments
and mean number of payments as follows:

.. code-block:: python

    transformer = MathematicalCombination(
        variables_to_combine=[
            'number_payments_first_quarter',
            'number_payments_second_quarter',
            'number_payments_third_quarter',
            'number_payments_fourth_quarter'
        ],
        math_operations=[
            'sum',
            'mean'
        ],
        new_variables_name=[
            'total_number_payments',
            'mean_number_payments'
        ]
    )

    Xt = transformer.fit_transform(X)


The transformed dataset, Xt, will contain the additional features
**total_number_payments** and **mean_number_payments**, plus the original set of
variables.

The variable **total_number_payments** is obtained by adding up the features
indicated in `variables_to_combine`, whereas the variable **mean_number_payments** is
the mean of those 4 features.

Below we show another example using the House Prices Dataset (more details about the
dataset :ref:`here <datasets>`). In this example, we sum 2 variables: 'LotFrontage' and
'LotArea' to obtain 'LotTotal'.

.. code:: python

    import pandas as pd
    from sklearn.model_selection import train_test_split

    from feature_engine.creation import MathematicalCombination

    data = pd.read_csv('houseprice.csv').fillna(0)

    X_train, X_test, y_train, y_test = train_test_split(
        data.drop(['Id', 'SalePrice'], axis=1),
        data['SalePrice'],
        test_size=0.3,
        random_state=0
    )

    math_combinator = MathematicalCombination(
        variables_to_combine=['LotFrontage', 'LotArea'],
        math_operations = ['sum'],
        new_variables_names = ['LotTotal']
    )

    math_combinator.fit(X_train, y_train)

    X_train_ = math_combinator.transform(X_train)


In the attribute `combination_dict_` the transformer stores the variable name and the
operation used to obtain that variable. This way, we can easily identify which variable
is the result of which transformation.

.. code:: python

    print(math_combinator.combination_dict_)

.. code:: python

    {'LotTotal': 'sum'}

We can see that the transformed dataset contains the additional variable:

.. code:: python

    print(X_train_.loc[:,['LotFrontage', 'LotArea', 'LotTotal']].head())

.. code:: python

          LotFrontage  LotArea  LotTotal
    64            0.0     9375    9375.0
    682           0.0     2887    2887.0
    960          50.0     7207    7257.0
    1384         60.0     9060    9120.0
    1100         60.0     8400    8460.0

**new_variables_names**

Even though the transfomer allows to combine variables automatically, it was originally
designed to combine variables with domain knowledge. In this case, we normally want to
give meaningful names to the variables. We can do so through the parameter
`new_variables_names`.

`new_variables_names` takes a list of strings, with the new variable names. In this
parameter, you need to enter a name or a list of names for the newly created features
(recommended). You must enter one name for each mathematical transformation indicated
in the `math_operations` parameter. That is, if you want to perform mean and sum of
features, you should enter 2 new variable names. If you perform only mean of features,
enter 1 variable name. Alternatively, if you chose to perform all mathematical
transformations, enter 6 new variable names.

The name of the variables should coincide with the order in which the
mathematical operations are initialised in the transformer. That is, if you set
math_operations = ['mean', 'prod'], the first new variable name will be
assigned to the mean of the variables and the second variable name
to the product of the variables.

More details
^^^^^^^^^^^^

You can find creative ways to use the :class:`MathematicalCombination()` in the
following Jupyter notebooks and Kaggle kernels.

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/creation/MathematicalCombination.ipynb>`_
- `Kaggle kernel - Wine Quality <https://www.kaggle.com/solegalli/create-new-features-with-feature-engine>`_
- `Kaggle kernel - House Price <https://www.kaggle.com/solegalli/feature-engineering-and-model-stacking>`_

.. _combine_with_ref:

.. currentmodule:: feature_engine.creation

CombineWithReferenceFeature
===========================

:class:`CombineWithReferenceFeature()` combines a group of variables with a group of
reference variables utilizing basic mathematical operations (subtraction, division,
addition and multiplication). It returns one or more additional features in the
dataframe as a result of these operations.

In other words, :class:`CombineWithReferenceFeature()` sums, multiplies, subtracts or
divides a group of features (indicated in `variables_to_combine`) to or by a group of
reference variables (indicated in `reference_variables`), and returns the
result as new variables in the dataframe.

For example, if we have the variables:

- **number_payments_first_quarter**,
- **number_payments_second_quarter**,
- **number_payments_third_quarter**,
- **number_payments_fourth_quarter**, and
- **total_payments**,

we can use :class:`CombineWithReferenceFeature()` to determine the percentage of
payments per quarter as follows:

.. code-block:: python

    transformer = CombineWithReferenceFeature(
        variables_to_combine=[
            'number_payments_first_quarter',
            'number_payments_second_quarter',
            'number_payments_third_quarter',
            'number_payments_fourth_quarter',
        ],

        reference_variables=['total_payments'],

        operations=['div'],

        new_variables_name=[
            'perc_payments_first_quarter',
            'perc_payments_second_quarter',
            'perc_payments_third_quarter',
            'perc_payments_fourth_quarter',
        ]
    )

    Xt = transformer.fit_transform(X)

The precedent code block will return a new dataframe, Xt, with 4 new variables, those
indicated in `new_variables_name`, that are calculated as the division of each one of
the variables in `variables_to_combine` and 'total_payments'.

Below we show another example using the House Prices Dataset (more details about the
dataset :ref:`here <datasets>`). In this example, we subtract `LotFrontage` from
`LotArea`.

.. code:: python

    import pandas as pd
    from sklearn.model_selection import train_test_split

    from feature_engine.creation import CombineWithReferenceFeature

    data = pd.read_csv('houseprice.csv').fillna(0)

    X_train, X_test, y_train, y_test = train_test_split(
    data.drop(['Id', 'SalePrice'], axis=1),
    data['SalePrice'],
    test_size=0.3,
    random_state=0
    )

    combinator = CombineWithReferenceFeature(
        variables_to_combine=['LotArea'],
        reference_variables=['LotFrontage'],
        operations = ['sub'],
        new_variables_names = ['LotPartial']
        )

    combinator.fit(X_train, y_train)

    X_train = combinator.transform(X_train)

We can see the newly created variable in the following code blocks:

.. code:: python

    print(X_train[["LotPartial","LotFrontage","LotArea"]].head())

.. code:: python

        LotTotal  LotFrontage  LotArea
    64      9375.0          0.0     9375
    682     2887.0          0.0     2887
    960     7157.0         50.0     7207
    1384    9000.0         60.0     9060
    1100    8340.0         60.0     8400

**new_variables_names**

Even though the transfomer allows to combine variables automatically, it was originally
designed to combine variables with domain knowledge. In this case, we normally want to
give meaningful names to the variables. We can do so through the parameter
`new_variables_names`.

`new_variables_names` takes a list of strings, with the new variable names. In this
parameter, you need to enter as many names as new features are created by the
transformer. The number of new features is the number of operations, times the number
of reference variables, times the number of variables to combine.

Thus, if you want to perform 2 operations, sub and div, combining 4 variables
with 2 reference variables, you should enter 2 X 4 X 2 new variable names.

The name of the variables should coincide with the order in which the operations are
performed by the transformer. The transformer will first carry out 'sub', then 'div',
then 'add' and finally 'mul'.

More details
^^^^^^^^^^^^

You can find creative ways to use the :class:`CombineWithReferenceFeature()` in the
following Jupyter notebooks and Kaggle kernels.

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/creation/CombineWithReferenceFeature.ipynb>`_
- `Kaggle kernel - Wine Quality <https://www.kaggle.com/solegalli/create-new-features-with-feature-engine>`_
- `Kaggle kernel - House Price <https://www.kaggle.com/solegalli/feature-engineering-and-model-stacking>`_
.. _cyclical_features:

.. currentmodule:: feature_engine.creation

CyclicalTransformer
===================

The :class:`CyclicalTransformer()` applies cyclical transformations to numerical
variables. The transformations return 2 new features per variable, according to:

- var_sin = sin(variable * (2. * pi / max_value))
- var_cos = cos(variable * (2. * pi / max_value))

where max_value is the maximum value in the variable, and pi is 3.14...

**Motivation**

There are some features that are cyclic by nature. For example the
hours of a day or the months in a year. In these cases, the higher values of
the variable are closer to the lower values. For example, December (12) is closer
to January (1) than to June (6). By applying a cyclical transformation we capture
this cycle or proximity between values.

**Examples**

In the code example below, we show how to obtain cyclical features from days and months
in a toy dataframe.

We first create a toy dataframe with the variables "days" and "months":

.. code:: python

    import pandas as pd
    from sklearn.model_selection import train_test_split

    from feature_engine.creation import CyclicalTransformer

    df = pd.DataFrame({
        'day': [6, 7, 5, 3, 1, 2, 4],
        'months': [3, 7, 9, 12, 4, 6, 12],
        })

Now we set up the transformer to find the maximum value automatically:

.. code:: python

    cyclical = CyclicalTransformer(variables=None, drop_original=True)

    X = cyclical.fit_transform(df)

The maximum values used for the transformation are stored in the attribute `max_values_`:

.. code:: python

    print(cyclical.max_values_)

.. code:: python

    {'day': 7, 'months': 12}

We can now see the new variables in the dataframe. Note that we set `drop_original=True`,
which means that the transformer will drop the original variables after the transformation.
If we had chosen False, the new variables will be added alongside the original ones.

.. code:: python

    print(X.head())

.. code:: python

          day_sin     day_cos  months_sin  months_cos
    1    -0.78183	  0.62349	      1.0	      0.0
    2         0.0	      1.0	     -0.5	 -0.86603
    3    -0.97493	-0.222521	     -1.0	     -0.0
    4     0.43388	-0.900969	      0.0	      1.0
    5     0.78183	  0.62349	  0.86603	     -0.5
    6     0.97493	-0.222521	      0.0	     -1.0
    7    -0.43388	-0.900969	      0.0	      1.0





.. -*- mode: rst -*-

Feature Creation
================

Feature-engine's creation transformers create and add new features to the dataframe
by either combining or transforming existing features.

.. toctree::
   :maxdepth: 1

   MathematicalCombination
   CombineWithReferenceFeature
   CyclicalTransformer.. _arbitrary_number_imputer:

.. currentmodule:: feature_engine.imputation

ArbitraryNumberImputer
======================

The :class:`ArbitraryNumberImputer()` replaces missing data with an arbitrary numerical
value determined by the user. It works only with numerical variables.

The :class:`ArbitraryNumberImputer()` can find and impute all numerical variables
automatically. Alternatively, you can pass a list of the variables you want to impute
to the `variables` parameter.

You can impute all variables with the same number, in which case you need to define
the variables to impute in the `variables` parameter and the imputation number in
`arbitrary_number` parameter. For example, you can impute varA and varB with 99
like this:

.. code-block:: python

    transformer = ArbitraryNumberImputer(
            variables = ['varA', 'varB'],
            arbitrary_number = 99
            )

    Xt = transformer.fit_transform(X)

You can also impute different variables with different numbers. To do this, you need to
pass a dictionary with the variable names and the numbers to use for their imputation
to the `imputer_dict` parameter. For example, you can impute varA with 1 and varB
with 99 like this:

.. code-block:: python

    transformer = ArbitraryNumberImputer(
            imputer_dict = {'varA' : 1, 'varB': 99]
            )

    Xt = transformer.fit_transform(X)


Below a code example using the House Prices Dataset (more details about the dataset
:ref:`here <datasets>`).

First, let's load the data and separate it into train and test:

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.imputation import ArbitraryNumberImputer

	# Load dataset
	data = pd.read_csv('houseprice.csv')

	# Separate into train and test sets
	X_train, X_test, y_train, y_test = train_test_split(
                                    data.drop(['Id', 'SalePrice'], axis=1),
                                    data['SalePrice'],
                                    test_size=0.3,
                                    random_state=0,
                                    )

Now we set up the :class:`ArbitraryNumberImputer()` to impute 2 variables from the
dataset with the number -999:

.. code:: python

	# set up the imputer
	arbitrary_imputer = ArbitraryNumberImputer(
            arbitrary_number=-999,
            variables=['LotFrontage', 'MasVnrArea'],
            )

	# fit the imputer
	arbitrary_imputer.fit(X_train)

With `fit()`, the transformer does not learn any parameter. It just assigns the imputation
values to each variable, which can be found in the attribute `imputer_dict_`.

With transform, we replace the missing data with the arbitrary values both in train and
test sets:

.. code:: python

	# transform the data
	train_t= arbitrary_imputer.transform(X_train)
	test_t= arbitrary_imputer.transform(X_test)

Note that after the imputation, if the percentage of missing values is relatively big,
the variable distribution will differ from the original one (in red the imputed
variable):

.. code:: python

	fig = plt.figure()
	ax = fig.add_subplot(111)
	X_train['LotFrontage'].plot(kind='kde', ax=ax)
	train_t['LotFrontage'].plot(kind='kde', ax=ax, color='red')
	lines, labels = ax.get_legend_handles_labels()
	ax.legend(lines, labels, loc='best')

.. image:: ../../images/arbitraryvalueimputation.png

More details
^^^^^^^^^^^^

In the following Jupyter notebook you will find more details on the functionality of the
:class:`ArbitraryNumberImputer()`, including how to select numerical variables automatically.
You will also see how to navigate the different attributes of the transformer.

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/imputation/ArbitraryNumberImputer.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _drop_missing_data:

.. currentmodule:: feature_engine.imputation

DropMissingData
===============


The :class:`DropMissingData()` will delete rows containing missing values. It provides
similar functionality to `pandas.drop_na()`. The transformer has however some
advantages over pandas:

- it learns and stores the variables for which the rows with na should be deleted
- it can be used within the Scikit-learn pipeline

It works with numerical and categorical variables. You can pass a list of variables to
impute, or the transformer will select and impute all variables.

The trasformer has the option to learn the variables with missing data in the train set,
and then remove observations with NA only in those variables. Or alternatively remove
observations with NA in all variables. You can change the behaviour using the parameter
`missing_only`.

This means that if you pass a list of variables to impute and set `missing_only=True`,
and some of the variables in your list do not have missing data in the train set,
missing data will not be removed during transform for those particular variables. In
other words, when `missing_only=True`, the transformer "double checks" that the entered
variables have missing data in the train set. If not, it ignores them during
`transform()`.

It is recommended to use `missing_only=True` when not passing a list of variables to
impute.

Below a code example using the House Prices Dataset (more details about the dataset
:ref:`here <datasets>`).

First, let's load the data and separate it into train and test:

.. code:: python

    import numpy as np
    import pandas as pd
    from sklearn.model_selection import train_test_split

    from feature_engine.imputation import DropMissingData

    # Load dataset
    data = pd.read_csv('houseprice.csv')

    # Separate into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
    data.drop(['Id', 'SalePrice'], axis=1),
    data['SalePrice'],
    test_size=0.3,
    random_state=0)

Now, we set up the imputer to remove observations if they have missing data in any of
the variables indicated in the list.

.. code:: python

    # set up the imputer
    missingdata_imputer = DropMissingData(variables=['LotFrontage', 'MasVnrArea'])

    # fit the imputer
    missingdata_imputer.fit(X_train)


Now, we can go ahead and add the missing indicators:

.. code:: python

    # transform the data
    train_t= missingdata_imputer.transform(X_train)
    test_t= missingdata_imputer.transform(X_test)

We can explore the number of observations with NA in the variable `LotFrontage` before
the imputation:

.. code:: python

    # Number of NA before the transformation
    X_train['LotFrontage'].isna().sum()

.. code:: python

    189

And after the imputation we should not have observations with NA:

.. code:: python

    # Number of NA after the transformation:
    train_t['LotFrontage'].isna().sum()

.. code:: python

    0

We can go ahead and compare the shapes of the different dataframes, before and after
the imputation, and we will see that the imputed data has less observations, because
those with NA in any of the 2 variables of interest were removed.

.. code:: python

    # Number of rows before and after transformation
    print(X_train.shape)
    print(train_t.shape)

.. code:: python

    (1022, 79)
    (829, 79)

Drop partially complete rows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default behaviour of :class:`DropMissingData()` will drop rows in NA is present in
any of the variables indicated in the list.

We have the option of dropping rows only if a certain percentage of values is missing
across all variables.

For example, if we set the parameter `threshold=0.5`, a row will be dropped if data is
missing in 50% of the variables. If we set the parameter `threshold=0.01`, a row will
be dropped if data is missing in 1% of the variables. If we set the parameter
`threshold=1`, a row will be dropped if data is missing in all the variables.


More details
^^^^^^^^^^^^
In the following Jupyter notebook you will find more details on the functionality of the
:class:`DropMissingData()`, including how to select numerical variables automatically.

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/imputation/DropMissingData.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _mean_median_imputer:

.. currentmodule:: feature_engine.imputation

MeanMedianImputer
=================

The :class:`MeanMedianImputer()` replaces missing data with the mean or median of the variable.
It works only with numerical variables. You can pass the list of variables you want to impute,
or alternatively, the imputer will automatically select all numerical variables in the
train set.

Note that in symetrical distributions, the mean and the median are very similar. But in
skewed distributions, the median is a better representation of the majority, as the mean
is biased to extreme values. The following image was taken from Wikipedia. The image links
to the use license.

.. figure::  ../../images/1024px-Relationship_between_mean_and_median_under_different_skewness.png
   :align:   center
   :target: https://commons.wikimedia.org/wiki/File:Relationship_between_mean_and_median_under_different_skewness.png


With the `fit()` method, the transformer learns and stores the mean or median values per
variable. Then it uses these values in the `transform()` method to transform the data.

Below a code example using the House Prices Dataset (more details about the dataset
:ref:`here <datasets>`).

First, let's load the data and separate it into train and test:

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.imputation import MeanMedianImputer

	# Load dataset
	data = pd.read_csv('houseprice.csv')

	# Separate into train and test sets
	X_train, X_test, y_train, y_test = train_test_split(
                                            data.drop(['Id', 'SalePrice'], axis=1),
                                            data['SalePrice'],
                                            test_size=0.3,
                                            random_state=0,
                                            )


Now we set up the :class:`MeanMedianImputer()` to impute in this case with the median
and only 2 variables from the dataset.

.. code:: python

	# set up the imputer
	median_imputer = MeanMedianImputer(
                           imputation_method='median',
                           variables=['LotFrontage', 'MasVnrArea']
                           )

	# fit the imputer
	median_imputer.fit(X_train)

With fit, the :class:`MeanMedianImputer()` learned the median values for the indicated
variables and stored it in one of its attributes. We can now go ahead and impute both
the train and the test sets.

.. code:: python

	# transform the data
	train_t= median_imputer.transform(X_train)
	test_t= median_imputer.transform(X_test)

Note that after the imputation, if the percentage of missing values is relatively big,
the variable distribution will differ from the original one (in red the imputed
variable):

.. code:: python

	fig = plt.figure()
	ax = fig.add_subplot(111)
	X_train['LotFrontage'].plot(kind='kde', ax=ax)
	train_t['LotFrontage'].plot(kind='kde', ax=ax, color='red')
	lines, labels = ax.get_legend_handles_labels()
	ax.legend(lines, labels, loc='best')

.. image:: ../../images/medianimputation.png

More details
^^^^^^^^^^^^

In the following Jupyter notebook you will find more details on the functionality of the
:class:`MeanMedianImputer()`, including how to select numerical variables automatically.
You will also see how to navigate the different attributes of the transformer to find the
mean or median values of the variables.

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/imputation/MeanMedianImputer.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _categorical_imputer:

.. currentmodule:: feature_engine.imputation

CategoricalImputer
==================

The :class:`CategoricalImputer()` replaces missing data in categorical variables with an
arbitrary value, like the string 'Missing' or by the most frequent category.

You can indicate which variables to impute passing the variable names in a list, or the
imputer automatically finds and selects all variables of type object and categorical.

Originally, we designed this imputer to work only with categorical variables. From version
1.1.0 we introduced the parameter `ignore_format` to allow the imputer to also impute
numerical variables with this functionality. This is, because in some cases, variables
that are by nature categorical, have numerical values.

Below a code example using the House Prices Dataset (more details about the dataset
:ref:`here <datasets>`).

In this example, we impute 2 variables from the dataset with the string 'Missing', which
is the default functionality of the transformer:

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.imputation import CategoricalImputer

	# Load dataset
	data = pd.read_csv('houseprice.csv')

	# Separate into train and test sets
	X_train, X_test, y_train, y_test = train_test_split(
    	data.drop(['Id', 'SalePrice'], axis=1), data['SalePrice'], test_size=0.3, random_state=0)

	# set up the imputer
	imputer = CategoricalImputer(variables=['Alley', 'MasVnrType'])

	# fit the imputer
	imputer.fit(X_train)

	# transform the data
	train_t= imputer.transform(X_train)
	test_t= imputer.transform(X_test)

	test_t['MasVnrType'].value_counts().plot.bar()

Note in the plot the presence of the category "Missing" which is added after the imputation:

.. image:: ../../images/missingcategoryimputer.png

More details
^^^^^^^^^^^^
In the following Jupyter notebook you will find more details on the functionality of the
:class:`EndTailImputer()`, including how to select numerical variables automatically.
You will also find demos on how to impute using the maximum value or the interquartile
range proximity rule.

Check also this `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/imputation/CategoricalImputer.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _end_tail_imputer:

.. currentmodule:: feature_engine.imputation

EndTailImputer
==============


The :class:`EndTailImputer()` replaces missing data with a value at the end of the distribution.
The value can be determined using the mean plus or minus a number of times the standard
deviation, or using the inter-quartile range proximity rule. The value can also be
determined as a factor of the maximum value.

You decide whether the missing data should be placed at the right or left tail of
the variable distribution.

In a sense, the :class:`EndTailImputer()` **"automates"** the work of the
:class:`ArbitraryNumberImputer()` because it will find automatically "arbitrary values"
far out at the end of the variable distributions.

:class:`EndTailImputer()` works only with numerical variables. You can impute only a
subset of the variables in the data by passing the variable names in a list. Alternatively,
the imputer will automatically select all numerical variables in the train set.


Below a code example using the House Prices Dataset (more details about the dataset
:ref:`here <datasets>`).

First, let's load the data and separate it into train and test:

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.imputation import EndTailImputer

	# Load dataset
	data = pd.read_csv('houseprice.csv')

	# Separate into train and test sets
	X_train, X_test, y_train, y_test = train_test_split(
                                            data.drop(['Id', 'SalePrice'], axis=1),
                                            data['SalePrice'],
                                            test_size=0.3,
                                            random_state=0,
                                            )


Now we set up the :class:`EndTailImputer()` to impute in this case only 2 variables
from the dataset. We instruct the imputer to find the imputation values using the mean
plus 3 times the standard deviation as follows:

.. code:: python

	# set up the imputer
	tail_imputer = EndTailImputer(imputation_method='gaussian',
                                  tail='right',
                                  fold=3,
                                  variables=['LotFrontage', 'MasVnrArea'])
	# fit the imputer
	tail_imputer.fit(X_train)


With fit, the :class:`EndTailImputer()` learned the imputation values for the indicated
variables and stored it in one of its attributes. We can now go ahead and impute both
the train and the test sets.

.. code:: python

	# transform the data
	train_t= tail_imputer.transform(X_train)
	test_t= tail_imputer.transform(X_test)


Note that after the imputation, if the percentage of missing values is relatively big,
the variable distribution will differ from the original one (in red the imputed
variable):

.. code:: python

	fig = plt.figure()
	ax = fig.add_subplot(111)
	X_train['LotFrontage'].plot(kind='kde', ax=ax)
	train_t['LotFrontage'].plot(kind='kde', ax=ax, color='red')
	lines, labels = ax.get_legend_handles_labels()
	ax.legend(lines, labels, loc='best')

.. image:: ../../images/endtailimputer.png

More details
^^^^^^^^^^^^

In the following Jupyter notebook you will find more details on the functionality of the
:class:`CategoricalImputer()`, including how to select numerical variables automatically,
how to impute with the most frequent category, and how to impute with a used defined
string.

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/imputation/EndTailImputer.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.

.. -*- mode: rst -*-

Missing Data Imputation
=======================

Feature-engine's missing data imputers replace missing data by parameters estimated
from data or arbitrary values pre-defined by the user. The following image summarizes
the main imputer's functionality.

.. figure::  ../../images/summary/imputersSummary.png
   :align:   center

|

In this guide, you will find code snippets to quickly be able to apply the imputers
to your datasets, as well as general knowledge and guidance on the imputation
techniques.


Imputers
~~~~~~~~

.. toctree::
   :maxdepth: 1

   MeanMedianImputer
   ArbitraryNumberImputer
   EndTailImputer
   CategoricalImputer
   RandomSampleImputer
   AddMissingIndicator
   DropMissingData.. _add_missing_indicator:

.. currentmodule:: feature_engine.imputation

AddMissingIndicator
===================


The :class:`AddMissingIndicator()` adds a binary variable indicating if observations are
missing (missing indicator). It adds missing indicators to both categorical and numerical
variables.

You can select the variables for which the missing indicators should be created passing
a variable list to the `variables` parameter. Alternatively, the imputer will
automatically select all variables.

The imputer has the option to add missing indicators to all variables or only to those
that have missing data in the train set. You can change the behaviour using the
parameter `missing_only`.

If `missing_only=True`, missing indicators will be added only to those variables with
missing data in the train set. This means that if you passed a variable list to
`variables` and some of those variables did not have missing data, no missing indicators
will be added to them. If it is paramount that all variables in your list get their
missing indicators, make sure to set `missing_only=False`.

It is recommended to use `missing_only=True` when not passing a list of variables to
impute.

Below a code example using the House Prices Dataset (more details about the dataset
:ref:`here <datasets>`).

First, let's load the data and separate it into train and test:

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.imputation import AddMissingIndicator

	# Load dataset
	data = pd.read_csv('houseprice.csv')


	# Separate into train and test sets
	X_train, X_test, y_train, y_test = train_test_split(
    	data.drop(['Id', 'SalePrice'], axis=1), data['SalePrice'], test_size=0.3, random_state=0)


Now we set up the imputer to add missing indicators to the 4 indicated variables:

.. code:: python

	# set up the imputer
	addBinary_imputer = AddMissingIndicator(
        variables=['Alley', 'MasVnrType', 'LotFrontage', 'MasVnrArea'],
        )

	# fit the imputer
	addBinary_imputer.fit(X_train)

Because we left the default value for `missing_only`, the :class:`AddMissingIndicator()`
will check if the variables indicated above have missing data in X_train. If they do,
missing indicators will be added for all 4 variables looking forward. If one of them
had not had missing data in X_train, missing indicators would have been added to the
remaining 3 variables only.

We can know which variables will have missing indicators by looking at the variable list
in the :class:`AddMissingIndicator()`'s attribute `variables_`.

Now, we can go ahead and add the missing indicators:

.. code:: python

	# transform the data
	train_t = addBinary_imputer.transform(X_train)
	test_t = addBinary_imputer.transform(X_test)

	train_t[['Alley_na', 'MasVnrType_na', 'LotFrontage_na', 'MasVnrArea_na']].head()


.. image:: ../../images/missingindicator.png
   :width: 500

Note that after adding missing indicators, we still need to replace NA in the original
variables if we plan to use them to train machine learning models.

Tip
^^^

Missing indicators are commonly used together with random sampling, mean or median
imputation, or frequent category imputation.

More details
^^^^^^^^^^^^

In the following Jupyter notebook you will find more details on the functionality of the
:class:`AddMissingIndicator()`, including how to use the parameter `missing_indicator` and
how to select the variables automatically.

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/imputation/AddMissingIndicator.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.

.. _random_sample_imputer:

.. currentmodule:: feature_engine.imputation

RandomSampleImputer
===================

The :class:`RandomSampleImputer()` replaces missing data with a random sample extracted from the
variable. It works with both numerical and categorical variables. A list of variables
can be indicated, or the imputer will automatically select all variables in the train
set.

**Note**

The random samples used to replace missing values may vary from execution to
execution. This may affect the results of your work. Thus, it is advisable to set a
seed.

Setting the seed
~~~~~~~~~~~~~~~~

There are 2 ways in which the seed can be set in the :class:`RandomSampleImputer()`:

If `seed = 'general'` then the random_state can be either `None` or an integer.
The `random_state` then provides the seed to use in the imputation. All observations will
be imputed in one go with a single seed. This is equivalent to
`pandas.sample(n, random_state=seed)` where `n` is the number of observations with
missing data and `seed` is the number you entered in the `random_state`.

If `seed = 'observation'`, then the random_state should be a variable name
or a list of variable names. The seed will be calculated observation per
observation, either by adding or multiplying the values of the variables
indicated in the `random_state`. Then, a value will be extracted from the train set
using that seed and used to replace the NAN in that particular observation. This is the
equivalent of `pandas.sample(1, random_state=var1+var2)` if the `seeding_method` is
set to `add` or `pandas.sample(1, random_state=var1*var2)` if the `seeding_method`
is set to `multiply`.

For example, if the observation shows variables color: np.nan, height: 152, weight:52,
and we set the imputer as:

.. code:: python

	RandomSampleImputer(random_state=['height', 'weight'],
                                  seed='observation',
                                  seeding_method='add'))

the np.nan in the variable colour will be replaced using pandas sample as follows:

.. code:: python

	observation.sample(1, random_state=int(152+52))

For more details on why this functionality is important refer to the course
`Feature Engineering for Machine Learning <https://www.udemy.com/feature-engineering-for-machine-learning/>`_.

You can also find more details about this imputation in the following
`notebook <https://github.com/solegalli/feature-engineering-for-machine-learning/blob/master/Section-04-Missing-Data-Imputation/04.07-Random-Sample-Imputation.ipynb>`_.

Note, if the variables indicated in the `random_state` list are not numerical
the imputer will return an error. In addition, the variables indicated as seed
should not contain missing values themselves.

Important for GDPR
~~~~~~~~~~~~~~~~~~

This estimator stores a copy of the training set when the `fit()` method is
called. Therefore, the object can become quite heavy. Also, it may not be GDPR
compliant if your training data set contains Personal Information. Please check
if this behaviour is allowed within your organisation.

Below a code example using the House Prices Dataset (more details about the dataset
:ref:`here <datasets>`).

First, let's load the data and separate it into train and test:

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.imputation import RandomSampleImputer

	# Load dataset
	data = pd.read_csv('houseprice.csv')

	# Separate into train and test sets
	X_train, X_test, y_train, y_test = train_test_split(
            data.drop(['Id', 'SalePrice'], axis=1),
            data['SalePrice'],
            test_size=0.3,
            random_state=0
        )

In this example, we sample values at random, observation per observation, using as seed
the value of the variable 'MSSubClass' plus the value of the variable 'YrSold'. Note
that this value might be different for each observation.

The :class:`RandomSampleImputer()` will impute all variables in the data, as we left the
default value of the parameter `variables` to `None`.

.. code:: python

	# set up the imputer
	imputer = RandomSampleImputer(
                random_state=['MSSubClass', 'YrSold'],
                seed='observation',
                seeding_method='add'
            )

	# fit the imputer
	imputer.fit(X_train)

With `fit()` the imputer stored a copy of the X_train. And with transform, it will extract
values at random from this X_train to replace NA in the datasets indicated in the `transform()`
methods.

.. code:: python

	# transform the data
	train_t = imputer.transform(X_train)
	test_t = imputer.transform(X_test)

The beauty of the random sampler is that it preserves the original variable distribution:

.. code:: python

	fig = plt.figure()
	ax = fig.add_subplot(111)
	X_train['LotFrontage'].plot(kind='kde', ax=ax)
	train_t['LotFrontage'].plot(kind='kde', ax=ax, color='red')
	lines, labels = ax.get_legend_handles_labels()
	ax.legend(lines, labels, loc='best')

.. image:: ../../images/randomsampleimputation.png

More details
^^^^^^^^^^^^

In the following Jupyter notebook you will find more details on the functionality of the
:class:`RandomSampleImputer()`, including how to set the different types of seeds.

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/imputation/RandomSampleImputer.ipynb>`_

All Feature-engine notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.

You will also find a lot of information on why the seed matters in this notebook:

- `notebook <https://github.com/solegalli/feature-engineering-for-machine-learning/blob/master/Section-04-Missing-Data-Imputation/04.07-Random-Sample-Imputation.ipynb>`_

And finally, there is also a lot of information about this and other imputation techniques
in this online course:

- `Feature Engineering for Machine Learning <https://www.udemy.com/feature-engineering-for-machine-learning/>`_
.. _winsorizer:

.. currentmodule:: feature_engine.outliers

Winsorizer
==========

The :class:`Winsorizer()` caps maximum and/or minimum values of a variable at automatically
determined values. The minimum and maximum values can be calculated in 1 of 3 different ways:

Gaussian limits:

- right tail: mean + 3* std
- left tail: mean - 3* std

IQR limits:

- right tail: 75th quantile + 3* IQR
- left tail:  25th quantile - 3* IQR

where IQR is the inter-quartile range: 75th quantile - 25th quantile.

percentiles or quantiles:

- right tail: 95th percentile
- left tail:  5th percentile

**Example**

Let's cap some outliers in the Titanic Dataset. First, let's load the data and separate
it into train and test:

.. code:: python

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.model_selection import train_test_split

    from feature_engine.outliers import Winsorizer

    # Load dataset
    def load_titanic():
        data = pd.read_csv('https://www.openml.org/data/get_csv/16826755/phpMYEkMl')
        data = data.replace('?', np.nan)
        data['cabin'] = data['cabin'].astype(str).str[0]
        data['pclass'] = data['pclass'].astype('O')
        data['embarked'].fillna('C', inplace=True)
        data['fare'] = data['fare'].astype('float')
        data['fare'].fillna(data['fare'].median(), inplace=True)
        data['age'] = data['age'].astype('float')
        data['age'].fillna(data['age'].median(), inplace=True)
        return data
	
    data = load_titanic()

    # Separate into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
		data.drop(['survived', 'name', 'ticket'], axis=1),
		data['survived'], test_size=0.3, random_state=0)

Now, we will set the :class:`Winsorizer()` to cap outliers at the right side of the
distribution only (param `tail`). We want the maximum values to be determined using the
mean value of the variable (param `capping_method`) plus 3 times the standard deviation
(param `fold`). And we only want to cap outliers in 2 variables, which we indicate in a
list.

.. code:: python

    # set up the capper
    capper = Winsorizer(capping_method='gaussian', tail='right', fold=3, variables=['age', 'fare'])

    # fit the capper
    capper.fit(X_train)

With `fit()`, the :class:`Winsorizer()` finds the values at which it should cap the variables.
These values are stored in its attribute:

.. code:: python

    capper.right_tail_caps_

.. code:: python

	{'age': 67.49048447470315, 'fare': 174.78162171790441}

We can now go ahead and censor the outliers:

.. code:: python

    # transform the data
    train_t= capper.transform(X_train)
    test_t= capper.transform(X_test)
    
If we evaluate now the maximum of the variables in the transformed datasets, they should
coincide with the values observed in the attribute `right_tail_caps_`:

.. code:: python

    train_t[['fare', 'age']].max()

.. code:: python

    fare    174.781622
    age      67.490484
    dtype: float64

More details
^^^^^^^^^^^^

You can find more details about the :class:`Winsorizer()` functionality in the following
notebook:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/outliers/Winsorizer.ipynb>`_
.. _arbitrary_capper:

.. currentmodule:: feature_engine.outliers

ArbitraryOutlierCapper
======================

The :class:`ArbitraryOutlierCapper()` caps the maximum or minimum values of a variable
at an arbitrary value indicated by the user. The maximum or minimum values should be
entered in a dictionary with the form {feature:capping value}.

Let's look at this in an example. First we load the Titanic dataset, and separate it
into a train and a test set:

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.outliers import ArbitraryOutlierCapper

	# Load dataset
	def load_titanic():
		data = pd.read_csv('https://www.openml.org/data/get_csv/16826755/phpMYEkMl')
		data = data.replace('?', np.nan)
		data['cabin'] = data['cabin'].astype(str).str[0]
		data['pclass'] = data['pclass'].astype('O')
		data['embarked'].fillna('C', inplace=True)
		data['fare'] = data['fare'].astype('float')
		data['fare'].fillna(data['fare'].median(), inplace=True)
		data['age'] = data['age'].astype('float')
		data['age'].fillna(data['age'].median(), inplace=True)
		return data
	
	data = load_titanic()

	# Separate into train and test sets
	X_train, X_test, y_train, y_test = train_test_split(
			data.drop(['survived', 'name', 'ticket'], axis=1),
			data['survived'], test_size=0.3, random_state=0)

Now, we set up the :class:`ArbitraryOutlierCapper()` indicating that we want to cap the
variable 'age' at 50 and the variable 'Fare' at 200. We do not want to cap these variables
on the left side of their distribution.

.. code:: python

	# set up the capper
	capper = ArbitraryOutlierCapper(max_capping_dict={'age': 50, 'fare': 200}, min_capping_dict=None)

	# fit the capper
	capper.fit(X_train)

With `fit()` the transformer does not learn any parameter. It just reassigns the entered
dictionary to the attribute that will be used in the transformation:

.. code:: python

	capper.right_tail_caps_

.. code:: python

	{'age': 50, 'fare': 200}

Now, we can go ahead and cap the variables:

.. code:: python

	# transform the data
	train_t= capper.transform(X_train)
	test_t= capper.transform(X_test)

If we now check the maximum values in the transformed data, they should be those entered
in the dictionary:

.. code:: python

    train_t[['fare', 'age']].max()

.. code:: python

    fare    200
    age      50
    dtype: float64


More details
^^^^^^^^^^^^

You can find more details about the :class:`ArbitraryOutlierCapper()` functionality in the following
notebook:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/outliers/ArbitraryOutlierCapper.ipynb>`_

.. -*- mode: rst -*-

Outlier Handling
================

Feature-engine's outlier cappers cap maximum or minimum values of a variable at an
arbitrary or derived value. The OutlierTrimmer removes outliers from the dataset.

.. toctree::
   :maxdepth: 1

   Winsorizer
   ArbitraryOutlierCapper
   OutlierTrimmer.. _outlier_trimmer:

.. currentmodule:: feature_engine.outliers

OutlierTrimmer
==============

The :class:`OutlierTrimmer()` removes values beyond an automatically generated
minimum and/or maximum values. The minimum and maximum values can be calculated in 1 of
3 ways:

Gaussian limits:

- right tail: mean + 3* std
- left tail: mean - 3* std

IQR limits:

- right tail: 75th quantile + 3* IQR
- left tail:  25th quantile - 3* IQR

where IQR is the inter-quartile range: 75th quantile - 25th quantile.

percentiles or quantiles:

- right tail: 95th percentile
- left tail:  5th percentile

**Example**

Let's remove some outliers in the Titanic Dataset. First, let's load the data and separate
it into train and test:

.. code:: python

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.model_selection import train_test_split

    from feature_engine.outliers import OutlierTrimmer

    # Load dataset
    def load_titanic():
        data = pd.read_csv('https://www.openml.org/data/get_csv/16826755/phpMYEkMl')
        data = data.replace('?', np.nan)
        data['cabin'] = data['cabin'].astype(str).str[0]
        data['pclass'] = data['pclass'].astype('O')
        data['embarked'].fillna('C', inplace=True)
        data['fare'] = data['fare'].astype('float')
        data['fare'].fillna(data['fare'].median(), inplace=True)
        data['age'] = data['age'].astype('float')
        data['age'].fillna(data['age'].median(), inplace=True)
        return data

    data = load_titanic()

    # Separate into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
		data.drop(['survived', 'name', 'ticket'], axis=1),
		data['survived'], test_size=0.3, random_state=0)

Now, we will set the :class:`OutlierTrimmer()` to remove outliers at the right side of the
distribution only (param `tail`). We want the maximum values to be determined using the
75th quantile of the variable (param `capping_method`) plus 1.5 times the IQR
(param `fold`). And we only want to cap outliers in 2 variables, which we indicate in a
list.

.. code:: python

    # set up the capper
    capper = OutlierTrimmer(capping_method='iqr', tail='right', fold=1.5, variables=['age', 'fare'])

    # fit the capper
    capper.fit(X_train)

With `fit()`, the :class:`OutlierTrimmer()` finds the values at which it should cap the variables.
These values are stored in its attribute:

.. code:: python

    capper.right_tail_caps_

.. code:: python

	{'age': 53.0, 'fare': 66.34379999999999}

We can now go ahead and remove the outliers:

.. code:: python

    # transform the data
    train_t= capper.transform(X_train)
    test_t= capper.transform(X_test)

If we evaluate now the maximum of the variables in the transformed datasets, they should
be <= the values observed in the attribute `right_tail_caps_`:

.. code:: python

    train_t[['fare', 'age']].max()

.. code:: python

    fare    65.0
    age     53.0
    dtype: float64

More details
^^^^^^^^^^^^

You can find more details about the :class:`OutlierTrimmer()` functionality in the following
notebook:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/outliers/OutlierTrimmer.ipynb>`_


.. _woe_encoder:

.. currentmodule:: feature_engine.encoding

WoEEncoder
==========

The :class:`WoEEncoder()` replaces categories by the weight of evidence
(WoE). The WoE was used primarily in the financial sector to create credit risk
scorecards.

The weight of evidence is given by:

.. math::

    log( p(X=xj|Y = 1) / p(X=xj|Y=0) )



**The WoE is determined as follows:**

We calculate the percentage positive cases in each category of the total of all
positive cases. For example 20 positive cases in category A out of 100 total
positive cases equals 20 %. Next, we calculate the percentage of negative cases in
each category respect to the total negative cases, for example 5 negative cases in
category A out of a total of 50 negative cases equals 10%. Then we calculate the
WoE by dividing the category percentages of positive cases by the category
percentage of negative cases, and take the logarithm, so for category A in our
example WoE = log(20/10).

**Note**

- If WoE values are negative, negative cases supersede the positive cases.
- If WoE values are positive, positive cases supersede the negative cases.
- And if WoE is 0, then there are equal number of positive and negative examples.

**Encoding into WoE**:

- Creates a monotonic relationship between the encoded variable and the target
- Returns variables in a similar scale

**Note**

This categorical encoding is exclusive for binary classification.

Let's look at an example using the Titanic Dataset.

First, let's load the data and separate it into train and test:


.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.encoding import WoEEncoder, RareLabelEncoder

	# Load dataset
	def load_titanic():
		data = pd.read_csv('https://www.openml.org/data/get_csv/16826755/phpMYEkMl')
		data = data.replace('?', np.nan)
		data['cabin'] = data['cabin'].astype(str).str[0]
		data['pclass'] = data['pclass'].astype('O')
		data['embarked'].fillna('C', inplace=True)
		return data
	
	data = load_titanic()

	# Separate into train and test sets
	X_train, X_test, y_train, y_test = train_test_split(
			data.drop(['survived', 'name', 'ticket'], axis=1),
			data['survived'], test_size=0.3, random_state=0)

Before we encode the variables, I would like to group infrequent categories into one
category, called 'Rare'. For this, I will use the :class:`RareLabelEncoder()` as follows:

.. code:: python

	# set up a rare label encoder
	rare_encoder = RareLabelEncoder(tol=0.03, n_categories=2, variables=['cabin', 'pclass', 'embarked'])

	# fit and transform data
	train_t = rare_encoder.fit_transform(X_train)
	test_t = rare_encoder.transform(X_train)

Now, we set up the :class:`WoEEncoder()` to replace the categories by the weight of the
evidence, only in the 3 indicated variables:

.. code:: python

	# set up a weight of evidence encoder
	woe_encoder = WoEEncoder(variables=['cabin', 'pclass', 'embarked'])

	# fit the encoder
	woe_encoder.fit(train_t, y_train)

With `fit()` the encoder learns the weight of the evidence for each category, which are stored in
its `encoder_dict_` parameter:

.. code:: python

	woe_encoder.encoder_dict_

In the `encoder_dict_` we find the WoE for each one of the categories of the
variables to encode. This way, we can map the original values to the new value.

.. code:: python

    {'cabin': {'B': 1.6299623810120747,
    'C': 0.7217038208351837,
    'D': 1.405081209799324,
    'E': 1.405081209799324,
    'Rare': 0.7387452866900354,
    'n': -0.35752781962490193},
    'pclass': {1: 0.9453018143294478,
    2: 0.21009172435857942,
    3: -0.5841726684724614},
    'embarked': {'C': 0.6999054533737715,
    'Q': -0.05044494288988759,
    'S': -0.20113381737960143}}

Now, we can go ahead and encode the variables:

.. code:: python

	# transform
	train_t = woe_encoder.transform(train_t)
	test_t = woe_encoder.transform(test_t)


**WoE for continuous variables**

In credit scoring, continuous variables are also transformed using the WoE. To do
this, first variables are sorted into a discrete number of bins, and then these
bins are encoded with the WoE as explained here for categorical variables. You can
do this by combining the use of the equal width, equal frequency or arbitrary
discretisers.

More details
^^^^^^^^^^^^

In the following notebooks, you can find more details into the :class:`WoEEncoder()`
functionality and example plots with the encoded variables:

- `WoE in categorical variables <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/encoding/WoEEncoder.ipynb>`_
- `WoE in numerical variables <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/discretisation/EqualFrequencyDiscretiser_plus_WoEEncoder.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _decisiontree_encoder:

.. currentmodule:: feature_engine.encoding

DecisionTreeEncoder
===================

The :class:`DecisionTreeEncoder()` replaces categories in the variable with
the predictions of a decision tree.

The transformer first encodes categorical variables into numerical variables using
:class:`OrdinalEncoder()`. You have the option to have the integers assigned to the
categories as they appear in the variable, or ordered by the mean value of the target
per category. You can regulate this behaviour with the parameter `encoding_method`. As
decision trees are able to pick non-linear relationships, replacing categories by
arbitrary numbers should be enough in practice.

After this, the transformer fits with this numerical variable a decision tree to predict
the target variable. Finally, the original categorical variable is replaced by the
predictions of the decision tree.

The motivation of the :class:`DecisionTreeEncoder()` is to try and create monotonic
relationships between the categorical variables and the target.

Let's look at an example using the Titanic Dataset.

First, let's load the data and separate it into train and test:

.. code:: python

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.model_selection import train_test_split

    from feature_engine.encoding import DecisionTreeEncoder

    # Load dataset
    def load_titanic():
            data = pd.read_csv('https://www.openml.org/data/get_csv/16826755/phpMYEkMl')
            data = data.replace('?', np.nan)
            data['cabin'] = data['cabin'].astype(str).str[0]
            data['pclass'] = data['pclass'].astype('O')
            data['embarked'].fillna('C', inplace=True)
            return data

    data = load_titanic()

    # Separate into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
                    data.drop(['survived', 'name', 'ticket'], axis=1),
                    data['survived'], test_size=0.3, random_state=0)

    X_train[['cabin', 'pclass', 'embarked']].head(10)

We will encode the following categorical variables:

.. code:: python

         cabin pclass embarked
   501      n      2        S
   588      n      2        S
   402      n      2        C
   1193     n      3        Q
   686      n      3        Q
   971      n      3        Q
   117      E      1        C
   540      n      2        S
   294      C      1        C
   261      E      1        S

We set up the encoder to encode the variables above with 3 fold cross-validation, using
a grid search to find the optimal depth of the decision tree (this is the default
behaviour of the :class:`DecisionTreeEncoder()`). In this example, we optimize the
tree using the roc-auc metric.

.. code:: python

    # set up the encoder
    encoder = DecisionTreeEncoder(
         variables=['cabin', 'pclass', 'embarked'],
         regression=False,
         scoring='roc_auc',
         cv=3,
         random_state=0)

    # fit the encoder
    encoder.fit(X_train, y_train)

With `fit()` the :class:`DecisionTreeEncoder()` fits 1 decision tree per variable. Now we can go ahead and
transform the categorical variables into numbers, using the predictions of these trees:

.. code:: python

    # transform the data
    train_t = encoder.transform(X_train)
    test_t = encoder.transform(X_test)

    train_t[['cabin', 'pclass', 'embarked']].head(10)

We can see the encoded variables below:

.. code:: python

             cabin    pclass  embarked
    501   0.304843  0.436170  0.338957
    588   0.304843  0.436170  0.338957
    402   0.304843  0.436170  0.558011
    1193  0.304843  0.259036  0.373494
    686   0.304843  0.259036  0.373494
    971   0.304843  0.259036  0.373494
    117   0.611650  0.617391  0.558011
    540   0.304843  0.436170  0.338957
    294   0.611650  0.617391  0.558011
    261   0.611650  0.617391  0.338957


More details
^^^^^^^^^^^^

In the following notebook, you can find more details into the :class:`DecisionTreeEncoder()`
functionality and example plots with the encoded variables:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/encoding/DecisionTreeEncoder.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _rarelabel_encoder:

.. currentmodule:: feature_engine.encoding

RareLabelEncoder
================

The :class:`RareLabelEncoder()` groups infrequent categories into one new category
called 'Rare' or a different string indicated by the user. We need to specify the
minimum percentage of observations a category should have to be preserved and the
minimum number of unique categories a variable should have to be re-grouped.

**tol**

In the parameter `tol` we indicate the minimum proportion of observations a category
should have, not to be grouped. In other words, categories which frequency, or proportion
of observations is <= `tol` will be grouped into a unique term.

**n_categories**

In the parameter `n_categories` we indicate the minimum cardinality of the categorical
variable in order to group infrequent categories. For example, if `n_categories=5`,
categories will be grouped only in those categorical variables with more than 5 unique
categories. The rest of the variables will be ignored.

This parameter is useful when we have big datasets and do not have time to examine all
categorical variables individually. This way, we ensure that variables with low cardinality
are not reduced any further.

**max_n_categories**

In the parameter `max_n_categories` we indicate the maximum number of unique categories
that we want in the encoded variable. If `max_n_categories=5`, then the most popular 5
categories will remain in the variable after the encoding, all other will be grouped into
a single category.

This parameter is useful if we are going to perform one hot encoding at the back of it,
to control the expansion of the feature space.

**Example**

Let's look at an example using the Titanic Dataset.

First, let's load the data and separate it into train and test:

.. code:: python

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.model_selection import train_test_split

    from feature_engine.encoding import RareLabelEncoder

    def load_titanic():
        data = pd.read_csv(
            'https://www.openml.org/data/get_csv/16826755/phpMYEkMl')
        data = data.replace('?', np.nan)
        data['cabin'] = data['cabin'].astype(str).str[0]
        data['pclass'] = data['pclass'].astype('O')
        data['embarked'].fillna('C', inplace=True)
        return data

    data = load_titanic()

    # Separate into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
        data.drop(['survived', 'name', 'ticket'], axis=1),
        data['survived'], test_size=0.3, random_state=0)

Now, we set up the :class:`RareLabelEncoder()` to group categories shown by less than 3%
of the observations into a new group or category called 'Rare'. We will group the
categories in the indicated variables if they have more than 2 unique categories each.

.. code:: python

    # set up the encoder
    encoder = RareLabelEncoder(tol=0.03, n_categories=2, variables=['cabin', 'pclass', 'embarked'],
                               replace_with='Rare')

    # fit the encoder
    encoder.fit(X_train)

With `fit()`, the :class:`RareLabelEncoder()` finds the categories present in more than
3% of the observations, that is, those that will not be grouped. These categories can
be found in the `encoder_dict_` attribute.

.. code:: python

    encoder.encoder_dict_

In the `encoder_dict_` we find the most frequent categories per variable to encode.
Any category that is not in this dictionary, will be grouped.

.. code:: python

	{'cabin': Index(['n', 'C', 'B', 'E', 'D'], dtype='object'),
	 'pclass': array([2, 3, 1], dtype='int64'),
	 'embarked': array(['S', 'C', 'Q'], dtype=object)}

Now we can go ahead and transform the variables:

.. code:: python

    # transform the data
    train_t = encoder.transform(X_train)
    test_t = encoder.transform(X_test)

We can also specify the maximum number of categories that can be considered frequent
using the `max_n_categories` parameter.

Let's begin by creating a toy dataframe and count the values of observations per
category:

.. code:: python

    from feature_engine.encoding import RareLabelEncoder
    import pandas as pd
    data = {'var_A': ['A'] * 10 + ['B'] * 10 + ['C'] * 2 + ['D'] * 1}
    data = pd.DataFrame(data)
    data['var_A'].value_counts()

.. code:: python

    A    10
    B    10
    C     2
    D     1
    Name: var_A, dtype: int64

In this block of code, we group the categories only for variables with more than 3
unique categories and then we plot the result:

.. code:: python

    rare_encoder = RareLabelEncoder(tol=0.05, n_categories=3)
    rare_encoder.fit_transform(data)['var_A'].value_counts()

.. code:: python

    A       10
    B       10
    C        2
    Rare     1
    Name: var_A, dtype: int64

Now, we retain the 2 most frequent categories of the variable and group the rest into
the 'Rare' group:

.. code:: python

    rare_encoder = RareLabelEncoder(tol=0.05, n_categories=3, max_n_categories=2)
    Xt = rare_encoder.fit_transform(data)
    Xt['var_A'].value_counts()

.. code:: python

    A       10
    B       10
    Rare     3
    Name: var_A, dtype: int64

Tips
^^^^

The :class:`RareLabelEncoder()` can be used to group infrequent categories and like this
control the expansion of the feature space if using one hot encoding.

Some categorical encodings will also return NAN if a category is present in the test
set, but was not seen in the train set. This inconvenient can usually be avoided if we
group rare labels before training the encoders.

Some categorical encoders will also return NAN if there is not enough observations for
a certain category. For example the :class:`WoEEncoder()` and the :class:`PRatioEncoder()`.
This behaviour can be also prevented by grouping infrequent labels before the encoding
with the :class:`RareLabelEncoder()`.


More details
^^^^^^^^^^^^

In the following notebook, you can find more details into the :class:`RareLabelEncoder()`
functionality and example plots with the encoded variables:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/encoding/RareLabelEncoder.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _ordinal_encoder:

.. currentmodule:: feature_engine.encoding

OrdinalEncoder
==============


The :class:`OrdinalEncoder()` replaces the categories by digits, starting from 0 to k-1,
where k is the number of different categories. If you select **"arbitrary"** in the
`encoding_method`, then the encoder will assign numbers as the labels appear in the
variable (first come first served). If you select **"ordered"**, the encoder will assign
numbers following the mean of the target value for that label. So labels for which the
mean of the target is higher will get the number 0, and those where the mean of the
target is smallest will get the number k-1. This way, we create a monotonic relationship
between the encoded variable and the target.

**Arbitrary vs ordered encoding**

**Ordered ordinal encoding**: for the variable colour, if the mean of the target
for blue, red and grey is 0.5, 0.8 and 0.1 respectively, blue is replaced by 1,
red by 2 and grey by 0.

The motivation is to try and create a monotonic relationship between the target and the
encoded categories. This tends to help improve performance of linear models.

**Arbitrary ordinal encoding**: the numbers will be assigned arbitrarily to the
categories, on a first seen first served basis.

Let's look at an example using the Titanic Dataset.

First, let's load the data and separate it into train and test:

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.encoding import OrdinalEncoder

	# Load dataset
	def load_titanic():
		data = pd.read_csv('https://www.openml.org/data/get_csv/16826755/phpMYEkMl')
		data = data.replace('?', np.nan)
		data['cabin'] = data['cabin'].astype(str).str[0]
		data['pclass'] = data['pclass'].astype('O')
		data['embarked'].fillna('C', inplace=True)
		return data
	
	data = load_titanic()

	# Separate into train and test sets
	X_train, X_test, y_train, y_test = train_test_split(
			data.drop(['survived', 'name', 'ticket'], axis=1),
			data['survived'], test_size=0.3, random_state=0)


Now, we set up the :class:`OrdinalEncoder()` to replace the categories by strings based
on the target mean value and only in the 3 indicated variables:

.. code:: python

	# set up the encoder
	encoder = OrdinalEncoder(encoding_method='ordered', variables=['pclass', 'cabin', 'embarked'])

	# fit the encoder
	encoder.fit(X_train, y_train)

With `fit()` the encoder learns the mappings for each category, which are stored in
its `encoder_dict_` parameter:

.. code:: python

	encoder.encoder_dict_

In the `encoder_dict_` we find the integers that will replace each one of the categories
of each variable that we want to encode. This way, we can map the original value of the
variable to the new value.

.. code:: python

	{'pclass': {3: 0, 2: 1, 1: 2},
	 'cabin': {'T': 0,
	  'n': 1,
	  'G': 2,
	  'A': 3,
	  'C': 4,
	  'F': 5,
	  'D': 6,
	  'E': 7,
	  'B': 8},
	 'embarked': {'S': 0, 'Q': 1, 'C': 2}}

We can now go ahead and replace the original strings with the numbers:

.. code:: python

	# transform the data
	train_t= encoder.transform(X_train)
	test_t= encoder.transform(X_test)




More details
^^^^^^^^^^^^
In the following notebook, you can find more details into the :class:`OrdinalEncoder()`'s
functionality and example plots with the encoded variables:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/encoding/OrdinalEncoder.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _mean_encoder:

.. currentmodule:: feature_engine.encoding

MeanEncoder
===========

The :class:`MeanEncoder()` replaces categories with the mean of the target per category.
For example, if we are trying to predict default rate, and our data has the variable city,
with categories, London, Manchester and Bristol, and the default rate per city is 0.1,
0.5, and 0.3, respectively, the encoder will replace London by 0.1, Manchester by 0.5
and Bristol by 0.3.

The motivation is to try and create a monotonic relationship between the target and
the encoded categories. This tends to help improve performance of linear models.

Let's look at an example using the Titanic Dataset.

First, let's load the data and separate it into train and test:

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.encoding import MeanEncoder

	# Load dataset
	def load_titanic():
		data = pd.read_csv('https://www.openml.org/data/get_csv/16826755/phpMYEkMl')
		data = data.replace('?', np.nan)
		data['cabin'] = data['cabin'].astype(str).str[0]
		data['pclass'] = data['pclass'].astype('O')
		data['embarked'].fillna('C', inplace=True)
		return data
	
	data = load_titanic()

	# Separate into train and test sets
	X_train, X_test, y_train, y_test = train_test_split(
			data.drop(['survived', 'name', 'ticket'], axis=1),
			data['survived'], test_size=0.3, random_state=0)

Now, we set up the :class:`MeanEncoder()` to replace the categories only in the 3
indicated variables:

.. code:: python

	# set up the encoder
	encoder = MeanEncoder(variables=['cabin', 'pclass', 'embarked'])

	# fit the encoder
	encoder.fit(X_train, y_train)

With `fit()` the encoder learns the target mean value for each category, which are stored in
its `encoder_dict_` parameter:

.. code:: python

	encoder.encoder_dict_

The `encoder_dict_` contains the mean value of the target per category, per variable.
So we can easily use this dictionary to map the numbers to the original labels.

.. code:: python

	{'cabin': {'A': 0.5294117647058824,
	  'B': 0.7619047619047619,
	  'C': 0.5633802816901409,
	  'D': 0.71875,
	  'E': 0.71875,
	  'F': 0.6666666666666666,
	  'G': 0.5,
	  'T': 0.0,
	  'n': 0.30484330484330485},
	 'pclass': {1: 0.6173913043478261,
	  2: 0.43617021276595747,
	  3: 0.25903614457831325},
	 'embarked': {'C': 0.5580110497237569,
	  'Q': 0.37349397590361444,
	  'S': 0.3389570552147239}}

We can now go ahead and replace the original strings with the numbers:

.. code:: python

	# transform the data
	train_t= encoder.transform(X_train)
	test_t= encoder.transform(X_test)


Handling Cardinality
^^^^^^^^^^^^^^^^^^^^

The :class:`MeanEncoder()` replaces categories with the mean of the target per category.
If the variable has low cardinality, then there is a fair representation of each label
in the dataset, and the mean target value per category can be determined with some certainty.
However, if variables are highly cardinal, with only very few observations for some labels,
then the mean target value for those categories will be unreliable.

To encode highly cardinal variables using target mean encoding, we could either group
infrequent categories first using the :class:`RareLabelEncoder()`. Alternatively, we
may want to choose different encoding methods that use blends of probabilities to try and
better estimate the encoding mappings, like those available in the open-source package
Category encoders through the transformers
`M-estimate <https://contrib.scikit-learn.org/category_encoders/mestimate.html>`_ and
`Target Encoder <https://contrib.scikit-learn.org/category_encoders/targetencoder.html>`_.


More details
^^^^^^^^^^^^

In the following notebook, you can find more details into the :class:`MeanEncoder()`
functionality and example plots with the encoded variables:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/encoding/MeanEncoder.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _count_freq_encoder:

.. currentmodule:: feature_engine.encoding

CountFrequencyEncoder
=====================

The :class:`CountFrequencyEncoder()` replaces categories by either the count or the
percentage of observations per category. For example in the variable colour, if 10
observations are blue, blue will be replaced by 10. Alternatively, if 10% of the
observations are blue, blue will be replaced by 0.1.

Let's look at an example using the Titanic Dataset.

First, let's load the data and separate it into train and test:

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.encoding import CountFrequencyEncoder

	# Load dataset
	def load_titanic():
		data = pd.read_csv('https://www.openml.org/data/get_csv/16826755/phpMYEkMl')
		data = data.replace('?', np.nan)
		data['cabin'] = data['cabin'].astype(str).str[0]
		data['pclass'] = data['pclass'].astype('O')
		data['embarked'].fillna('C', inplace=True)
		return data
	
	data = load_titanic()

	# Separate into train and test sets
	X_train, X_test, y_train, y_test = train_test_split(
			data.drop(['survived', 'name', 'ticket'], axis=1),
			data['survived'], test_size=0.3, random_state=0)


Now, we set up the :class:`CountFrequencyEncoder()` to replace the categories by their
frequencies, only in the 3 indicated variables:

.. code:: python

	# set up the encoder
	encoder = CountFrequencyEncoder(encoding_method='frequency',
				 variables=['cabin', 'pclass', 'embarked'])

	# fit the encoder
	encoder.fit(X_train)

With `fit()` the encoder learns the frequencies of each category, which are stored in
its `encoder_dict_` parameter:

.. code:: python

	encoder.encoder_dict_

In the `encoder_dict_` we find the frequencies for each one of the categories of each
variable that we want to encode. This way, we can map the original value to the new
value.

.. code:: python

	{'cabin': {'n': 0.7663755458515283,
	  'C': 0.07751091703056769,
	  'B': 0.04585152838427948,
	  'E': 0.034934497816593885,
	  'D': 0.034934497816593885,
	  'A': 0.018558951965065504,
	  'F': 0.016375545851528384,
	  'G': 0.004366812227074236,
	  'T': 0.001091703056768559},
	 'pclass': {3: 0.5436681222707423,
	  1: 0.25109170305676853,
	  2: 0.2052401746724891},
	 'embarked': {'S': 0.7117903930131004,
	  'C': 0.19759825327510916,
	  'Q': 0.0906113537117904}}

We can now go ahead and replace the original strings with the numbers:

.. code:: python

	# transform the data
	train_t= encoder.transform(X_train)
	test_t= encoder.transform(X_test)


More details
^^^^^^^^^^^^

In the following notebook, you can find more details into the :class:`CountFrequencyEncoder()`
functionality and example plots with the encoded variables:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/encoding/CountFrequencyEncoder.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _onehot_encoder:

.. currentmodule:: feature_engine.encoding


OneHotEncoder
=============

The :class:`OneHotEncoder()` performs one hot encoding. One hot encoding consists in
replacing the categorical variable by a group of binary variables which take value 0 or
1, to indicate if a certain category is present in an observation. The binary variables
are also known as dummy variables.

For example, from the categorical variable "Gender" with categories "female" and
"male", we can generate the boolean variable "female", which takes 1 if the
observation is female or 0 otherwise. We can also generate the variable "male",
which takes 1 if the observation is "male" and 0 otherwise. By default, the
:class:`OneHotEncoder()` will return both binary variables from "Gender": "female" and
"male".

**Binary variables**

When a categorical variable has only 2 categories, like "Gender" in our previous example, then
the second dummy variable created by one hot encoding can be completely redundant. We
can drop automatically the last dummy variable for those variables that contain only 2
categories by setting the parameter `drop_last_binary=True`. This will ensure that for
every binary variable in the dataset, only 1 dummy is created. This is recommended,
unless we suspect that the variable could, in principle take more than 2 values.

**k vs k-1 dummies**

From a categorical variable with k unique categories, the :class:`OneHotEncoder()` can
create k binary variables, or alternatively k-1 to avoid redundant information. This
behaviour can be specified using the parameter `drop_last`. Only k-1 binary variables
are necessary to encode all of the information in the original variable. However, there
are situations in which we may choose to encode the data into k dummies.

Encode into k-1 if training linear models: Linear models evaluate all features during
fit, thus, with k-1 they have all information about the original categorical variable.

Encode into k if training decision trees or performing feature selection: tree based
models and many feature selection algorithms evaluate variables or groups of variables
separately. Thus, if encoding into k-1, the last category will not be examined. That is,
we lose the information contained in that category.

**Encoding only popular categories**

The encoder can also create binary variables for the n most popular categories, n being
determined by the user. For example, if we encode only the 6 more popular categories, by
setting the parameter `top_categories=6`, the transformer will add binary variables only
for the 6 most frequent categories. The most frequent categories are those with the biggest
number of observations. The remaining categories will not be encoded into dummies. Thus,
if an observation presents a category other than the most frequent ones, it will have a
0 value in each one of the derived dummies. This behaviour is useful when the categorical
variables are highly cardinal, to control the expansion of the feature space.

**Note**

Only when creating binary variables for all categories of the variable (instead of the
most popular ones), we can specify if we want to encode into k or k-1 binary variables,
where k is the number if unique categories. If we encode only the top n most popular
categories, the encoder will create only n binary variables per categorical variable.
Observations that do not show any of these popular categories, will have 0 in all
the binary variables.

Let's look at an example using the Titanic Dataset. First we load the data and divide it
into a train and a test set:

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.encoding import OneHotEncoder

	# Load dataset
	def load_titanic():
		data = pd.read_csv('https://www.openml.org/data/get_csv/16826755/phpMYEkMl')
		data = data.replace('?', np.nan)
		data['cabin'] = data['cabin'].astype(str).str[0]
		data['pclass'] = data['pclass'].astype('O')
		data['embarked'].fillna('C', inplace=True)
		return data
	
	data = load_titanic()

	# Separate into train and test sets
	X_train, X_test, y_train, y_test = train_test_split(
				data.drop(['survived', 'name', 'ticket'], axis=1),
				data['survived'], test_size=0.3, random_state=0)

Now, we set up the encoder to encode only the 2 most frequent categories of each of the
3 indicated categorical variables:

.. code:: python

	# set up the encoder
	encoder = OneHotEncoder(top_categories=2, variables=['pclass', 'cabin', 'embarked'])

	# fit the encoder
	encoder.fit(X_train)

With `fit()` the encoder will learn the most popular categories of the variables, which
are stored in the attribute `encoder_dict_`.

.. code:: python

	encoder.encoder_dict_

.. code:: python

	{'pclass': [3, 1], 'cabin': ['n', 'C'], 'embarked': ['S', 'C']}

The `encoder_dict_` contains the categories that will derive dummy variables for each
categorical variable.

With transform, we go ahead and encode the variables. Note that by default, the
:class:`OneHotEncoder()` will drop the original variables.

.. code:: python

	# transform the data
	train_t= encoder.transform(X_train)
	test_t= encoder.transform(X_test)

If you do not want to drop the original variables, consider using the OneHotEncoder
from Scikit-learn and wrap it with the :ref:`SklearnTransformerWrapper <sklearn_wrapper>`.

**Feature space and duplication**

If the categorical variables are highly cardinal, we may end up with very big datasets
after one hot encoding. In addition, if some of these variables are fairly constant or
fairly similar, we may end up with one hot encoded features that are highly correlated
if not identical.

Consider checking this up and dropping redundant features with the transformers from the
:ref:`selection module <selection_user_guide>`.

More details
^^^^^^^^^^^^

For more details into :class:`OneHotEncoder()`'s functionality visit:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/encoding/OneHotEncoder.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. -*- mode: rst -*-

.. currentmodule:: feature_engine.encoding


Categorical Encoding
====================

Feature-engine's categorical encoders replace variable strings by estimated or
arbitrary numbers. The following image summarizes the main encoder’s functionality.

.. figure::  ../../images/summary/categoricalSummary.png
   :align:   center

   Summary of Feature-engine's encoders main characteristics

Feature-engine's categorical encoders work only with categorical variables by default.
From version 1.1.0, you have the option to set the parameter ignore_format to False,
and make the transformers also accept numerical variables as input.

**Monotonicity**

Most Feature-engine's encoders will return, or attempt to return monotonic relationships
between the encoded variable and the target. A monotonic relationship is one in which
the variable value increases as the values in the other variable increase, or decrease.
See the following illustration as examples:

.. figure::  ../../images/monotonic.png
   :align:   center
   :width: 400

Monotonic relationships tend to help improve the performance of linear models and build
shallower decision trees.

**Regression vs Classification**

Most Feature-engine's encoders are suitable for both regression and classification, with
the exception of the :class:`WoEEncoder()` and the :class:`PRatioEncoder()` which are
designed solely for **binary** classification.

**Multi-class classification**

Finally, some Feature-engine's encoders can handle multi-class targets off-the-shelf for
example the :class:`OneHotEncoder()`, the :class:CountFrequencyEncoder()` and the
:class:`DecisionTreeEncoder()`.

Note that while the :class:`MeanEncoder()` and the :class:`OrdinalEncoder()` will operate
with multi-class targets, but the mean of the classes may not be significant and this will
defeat the purpose of these encoding techniques.

**Encoders**

.. toctree::
   :maxdepth: 1

   OneHotEncoder
   CountFrequencyEncoder
   OrdinalEncoder
   MeanEncoder
   WoEEncoder
   PRatioEncoder
   DecisionTreeEncoder
   RareLabelEncoder


Additional categorical encoding transformations ara available in the open-source package
`Category encoders <https://contrib.scikit-learn.org/category_encoders/>`_.
   
.. _pratio_encoder:

.. currentmodule:: feature_engine.encoding

PRatioEncoder
=============

The :class:`PRatioEncoder()` replaces categories by the ratio of the probability of the
target = 1 and the probability of the target = 0.

The target probability ratio is given by:

.. math::

    p(1) / p(0)

The log of the target probability ratio is:

.. math::

    log( p(1) / p(0) )


For example in the variable colour, if the mean of the target = 1 for blue
is 0.8 and the mean of the target = 0  is 0.2, blue will be replaced by:
0.8 / 0.2 = 4 if ratio is selected, or log(0.8/0.2) = 1.386 if log_ratio
is selected.

**Note**

This categorical encoding is exclusive for binary classification.

Let's look at an example using the Titanic Dataset.

First, let's load the data and separate it into train and test:

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.encoding import PRatioEncoder, RareLabelEncoder

	# Load dataset
	def load_titanic():
		data = pd.read_csv('https://www.openml.org/data/get_csv/16826755/phpMYEkMl')
		data = data.replace('?', np.nan)
		data['cabin'] = data['cabin'].astype(str).str[0]
		data['pclass'] = data['pclass'].astype('O')
		data['embarked'].fillna('C', inplace=True)
		return data

	data = load_titanic()

	# Separate into train and test sets
	X_train, X_test, y_train, y_test = train_test_split(
			data.drop(['survived', 'name', 'ticket'], axis=1),
			data['survived'], test_size=0.3, random_state=0)

Before we encode the variables, I would like to group infrequent categories into one
category, called 'Rare'. For this, I will use the :class:`RareLabelEncoder()` as follows:

.. code:: python

	# set up a rare label encoder
	rare_encoder = RareLabelEncoder(tol=0.03, n_categories=2, variables=['cabin', 'pclass', 'embarked'])

	# fit and transform data
	train_t = rare_encoder.fit_transform(X_train)
	test_t = rare_encoder.transform(X_train)

Now, we set up the :class:`PRatioEncoder()` to replace the categories by the probability
ratio, only in the 3 indicated variables:

.. code:: python

	# set up a weight of evidence encoder
	pratio_encoder = PRatioEncoder(encoding_method='ratio', variables=['cabin', 'pclass', 'embarked'])

	# fit the encoder
	pratio_encoder.fit(train_t, y_train)

With `fit()` the encoder learns the values to replace each category, which are stored in
its `encoder_dict_` parameter:

.. code:: python

	pratio_encoder.encoder_dict_

In the `encoder_dict_` we find the probability ratio for each category in each
variable to encode. This way, we can map the original value to the new value.

.. code:: python

    {'cabin': {'B': 3.1999999999999993,
     'C': 1.2903225806451615
     'D': 2.5555555555555554,
     'E': 2.5555555555555554,
     'Rare': 1.3124999999999998,
     'n': 0.4385245901639344},
     'pclass': {1: 1.6136363636363635,
      2: 0.7735849056603774,
      3: 0.34959349593495936},
      'embarked': {'C': 1.2625000000000002,
      'Q': 0.5961538461538461,
      'S': 0.5127610208816704}}

Now, we can go ahead and encode the variables:

.. code:: python

	# transform
	train_t = pratio_encoder.transform(train_t)
	test_t = pratio_encoder.transform(test_t)


More details
^^^^^^^^^^^^

In the following notebook, you can find more details into the :class:`PRatioEncoder()`
functionality and example plots with the encoded variables:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/encoding/PRatioEncoder.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _sklearn_wrapper:

.. currentmodule:: feature_engine.wrappers

SklearnTransformerWrapper
=========================

The :class:`SklearnTransformerWrapper()` applies Scikit-learn transformers to a selected
group of variables. It works with transformers like the SimpleImputer, OrdinalEncoder,
OneHotEncoder, KBinsDiscretizer, all scalers and also transformers for feature selection.
Other transformers have not been tested, but we think it should work with most of them.

The :class:`SklearnTransformerWrapper()` offers similar functionality to the
`ColumnTransformer <https://scikit-learn.org/stable/modules/generated/sklearn.compose.ColumnTransformer.html>`_
class available in Scikit-learn. They differ in the implementation to select the
variables and the output.

The :class:`SklearnTransformerWrapper()` returns a pandas dataframe with the variables
in the order of the original data. The
`ColumnTransformer <https://scikit-learn.org/stable/modules/generated/sklearn.compose.ColumnTransformer.html>`_
returns a Numpy array, and the order of the variables may not coincide with that of the
original dataset.

In the next code snippet we show how to wrap the SimpleImputer from Scikit-learn to
impute only the selected variables.

.. code:: python

    import pandas as pd
    import numpy as np
    from sklearn.model_selection import train_test_split
    from sklearn.impute import SimpleImputer
    from feature_engine.wrappers import SklearnTransformerWrapper
	
    # Load dataset
    data = pd.read_csv('houseprice.csv')
    
    # Separate into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
    	data.drop(['Id', 'SalePrice'], axis=1),
    	data['SalePrice'], test_size=0.3, random_state=0)
    	
    # set up the wrapper with the SimpleImputer
    imputer = SklearnTransformerWrapper(transformer = SimpleImputer(strategy='mean'),
                                        variables = ['LotFrontage', 'MasVnrArea'])
    
    # fit the wrapper + SimpleImputer                              
    imputer.fit(X_train)
	
    # transform the data
    X_train = imputer.transform(X_train)
    X_test = imputer.transform(X_test)


In the next snippet of code we show how to wrap the StandardScaler from Scikit-learn
to standardize only the selected variables.

.. code:: python

    import pandas as pd
    import numpy as np
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler
    from feature_engine.wrappers import SklearnTransformerWrapper

    # Load dataset
    data = pd.read_csv('houseprice.csv')

    # Separate into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
    	data.drop(['Id', 'SalePrice'], axis=1),
    	data['SalePrice'], test_size=0.3, random_state=0)

    # set up the wrapper with the StandardScaler
    scaler = SklearnTransformerWrapper(transformer = StandardScaler(),
                                        variables = ['LotFrontage', 'MasVnrArea'])

    # fit the wrapper + StandardScaler
    scaler.fit(X_train)

    # transform the data
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)


In the next snippet of code we show how to wrap the SelectKBest from Scikit-learn
to select only a subset of the variables.

.. code:: python

    import pandas as pd
    import numpy as np
    from sklearn.model_selection import train_test_split
    from sklearn.feature_selection import f_regression, SelectKBest
    from feature_engine.wrappers import SklearnTransformerWrapper

    # Load dataset
    data = pd.read_csv('houseprice.csv')

    # Separate into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
    	data.drop(['Id', 'SalePrice'], axis=1),
    	data['SalePrice'], test_size=0.3, random_state=0)

    cols = [var for var in X_train.columns if X_train[var].dtypes !='O']

    # let's apply the standard scaler on the above variables

    selector = SklearnTransformerWrapper(
        transformer = SelectKBest(f_regression, k=5),
        variables = cols)

    selector.fit(X_train.fillna(0), y_train)

    # transform the data
    X_train_t = selector.transform(X_train.fillna(0))
    X_test_t = selector.transform(X_test.fillna(0))


More details
^^^^^^^^^^^^

In the following Jupyter notebooks you can find more details about how to navigate the
parameters of the :class:`SklearnTransformerWrapper()` and also access the parameters
of the Scikit-learn transformer wrapped, as well as the output of the transformations.

- `Wrap sklearn categorical encoder <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/wrappers/Sklearn-wrapper-plus-Categorical-Encoding.ipynb>`_
- `Wrap sklearn KBinsDiscretizer <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/wrappers/Sklearn-wrapper-plus-KBinsDiscretizer.ipynb>`_
- `Wrap sklearn SimpleImputer <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/wrappers/Sklearn-wrapper-plus-SimpleImputer.ipynb>`_
- `Wrap sklearn feature selectors <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/wrappers/Sklearn-wrapper-plus-feature-selection.ipynb>`_
- `Wrap sklearn scalers <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/wrappers/Sklearn-wrapper-plus-scalers.ipynb>`_

The notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. -*- mode: rst -*-

Scikit-learn Wrapper
====================

Feature-engine's Scikit-learn wrappers wrap Scikit-learn transformers allowing their
implementation only on a selected subset of features.

.. toctree::
   :maxdepth: 1

   Wrapper.. _equal_width_discretiser:

.. currentmodule:: feature_engine.discretisation

EqualWidthDiscretiser
=====================

The :class:`EqualWidthDiscretiser()` sorts the variable values into contiguous intervals
of equal size. The size of the interval is calculated as:

( max(X) - min(X) ) / bins

where bins, which is the number of intervals, should be determined by the user. The
interval limits are determined using `pandas.cut()`.

**A note on number of intervals**

Common values are 5 and 10. Note that if the variable is highly skewed or not continuous
smaller intervals maybe required. Otherwise, the transformer will introduce np.nan.

The :class:`EqualWidthDiscretiser()` works only with numerical variables. A list of
variables to discretise can be indicated, or the discretiser will automatically select
all numerical variables in the train set.

**Example**

Let's look at an example using the House Prices Dataset (more details about the
dataset :ref:`here <datasets>`).

Let's load the house prices dataset and  separate it into train and test sets:

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.discretisation import EqualWidthDiscretiser

	# Load dataset
	data = data = pd.read_csv('houseprice.csv')

	# Separate into train and test sets
	X_train, X_test, y_train, y_test =  train_test_split(
		    data.drop(['Id', 'SalePrice'], axis=1),
		    data['SalePrice'], test_size=0.3, random_state=0)


Now we want to discretise the 2 variables indicated below into 10 intervals of equal
width:

.. code:: python

	# set up the discretisation transformer
	disc = EqualWidthDiscretiser(bins=10, variables=['LotArea', 'GrLivArea'])

	# fit the transformer
	disc.fit(X_train)

With `fit()` the transformer learns the boundaries of each interval. Then, we can go
ahead and sort the values into the intervals:

.. code:: python

	# transform the data
	train_t= disc.transform(X_train)
	test_t= disc.transform(X_test)

The `binner_dict_` stores the interval limits identified for each variable.

.. code:: python

	disc.binner_dict_

.. code:: python

	'LotArea': [-inf,
	  22694.5,
	  44089.0,
	  65483.5,
	  86878.0,
	  108272.5,
	  129667.0,
	  151061.5,
	  172456.0,
	  193850.5,
	  inf],
	 'GrLivArea': [-inf,
	  768.2,
	  1202.4,
	  1636.6,
	  2070.8,
	  2505.0,
	  2939.2,
	  3373.4,
	  3807.6,
	  4241.799999999999,
	  inf]}

With equal width discretisation, each bin does not necessarily contain the same number of observations.

.. code:: python

	train_t.groupby('GrLivArea')['GrLivArea'].count().plot.bar()
	plt.ylabel('Number of houses')

We can see below that the intervals contain different number of observations.

.. image:: ../../images/equalwidthdiscretisation.png

|

**Discretisation plus encoding**

If we return the interval values as integers, the discretiser has the option to return
the transformed variable as integer or as object. Why would we want the transformed
variables as object?

Categorical encoders in Feature-engine are designed to work with variables of type
object by default. Thus, if you wish to encode the returned bins further, say to try and
obtain monotonic relationships between the variable and the target, you can do so
seamlessly by setting `return_object` to True. You can find an example of how to use
this functionality `here <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/discretisation/EqualWidthDiscretiser_plus_OrdinalEncoder.ipynb>`_.

More details
^^^^^^^^^^^^

Check also for more details on how to use this transformer:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/discretisation/EqualWidthDiscretiser.ipynb>`_
- `Jupyter notebook - Discretiser plus Ordinal encoding <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/discretisation/EqualWidthDiscretiser_plus_OrdinalEncoder.ipynb>`_
.. _arbitrary_discretiser:

.. currentmodule:: feature_engine.discretisation

ArbitraryDiscretiser
====================

The :class:`ArbitraryDiscretiser()` sorts the variable values into contiguous intervals
which limits are arbitrarily defined by the user. Thus, you must provide a dictionary
with the variable names as keys and the limits of the intervals in a list as values,
when setting up the discretiser.

The :class:`ArbitraryDiscretiser()` works only with numerical variables. The discretiser
will check that the variables entered by the user are present in the train set and cast
as numerical.

**Example**

Let's take a look at how this transformer works. First, let's load a dataset and plot a
histogram of a continuous variable. We use the boston house prices dataset that comes
with Scikit-learn.

.. code:: python

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.datasets import load_boston
    from feature_engine.discretisation import ArbitraryDiscretiser

    boston_dataset = load_boston()
    data = pd.DataFrame(boston_dataset.data, columns=boston_dataset.feature_names)

    data['LSTAT'].hist(bins=20)
    plt.xlabel('LSTAT')
    plt.ylabel('Number of observations')
    plt.title('Histogram of LSTAT')
    plt.show()

.. image:: ../../images/lstat_hist.png

Now, let's discretise the variable into arbitrarily determined intervals. We want the
interval names as integers, so we set `return_boundaries` to False.

.. code:: python

    user_dict = {'LSTAT': [0, 10, 20, 30, np.Inf]}

    transformer = ArbitraryDiscretiser(
        binning_dict=user_dict, return_object=False, return_boundaries=False)

    X = transformer.fit_transform(data)

Now, we can go ahead and plot the variable after the transformation:

.. code:: python

    X['LSTAT'].value_counts().plot.bar()
    plt.xlabel('LSTAT - bins')
    plt.ylabel('Number of observations')
    plt.title('Discretised LSTAT')
    plt.show()

.. image:: ../../images/lstat_disc_arbitrarily.png

Note that in the above figure the intervals are represented by digits.

Alternatively, we can return the interval limits in the discretised variable by
setting `return_boundaries` to True.

.. code:: python

    transformer = ArbitraryDiscretiser(
        binning_dict=user_dict, return_object=False, return_boundaries=True)
    X = transformer.fit_transform(data)

    X['LSTAT'].value_counts().plot.bar(rot=0)
    plt.xlabel('LSTAT - bins')
    plt.ylabel('Number of observations')
    plt.title('Discretised LSTAT')
    plt.show()

.. image:: ../../images/lstat_disc_arbitrarily2.png

**Discretisation plus encoding**

If we return the interval values as integers, the discretiser has the option to return
the transformed variable as integer or as object. Why would we want the transformed
variables as object?

Categorical encoders in Feature-engine are designed to work with variables of type
object by default. Thus, if you wish to encode the returned bins further, say to try and
obtain monotonic relationships between the variable and the target, you can do so
seamlessly by setting `return_object` to True. You can find an example of how to use
this functionality `here <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/discretisation/ArbitraryDiscretiser_plus_MeanEncoder.ipynb>`_.

More details
^^^^^^^^^^^^

Check also:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/discretisation/ArbitraryDiscretiser.ipynb>`_
- `Jupyter notebook - Discretiser plus Mean Encoding <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/discretisation/ArbitraryDiscretiser_plus_MeanEncoder.ipynb>`_
.. _decisiontree_discretiser:

.. currentmodule:: feature_engine.discretisation

DecisionTreeDiscretiser
=======================

The :class:`DecisionTreeDiscretiser()` replaces numerical variables by discrete, i.e.,
finite variables, which values are the predictions of a decision tree. The method is
based on the winning solution of the KDD 2009 competition:

`Niculescu-Mizil, et al. "Winning the KDD Cup Orange Challenge with Ensemble
Selection". JMLR: Workshop and Conference Proceedings 7: 23-34. KDD 2009
<http://proceedings.mlr.press/v7/niculescu09/niculescu09.pdf>`_.

In the original article, each feature in the challenge dataset was re-coded by training
a decision tree of limited depth (2, 3 or 4) using that feature alone, and letting the
tree predict the target. The probabilistic predictions of this decision tree were used
as an additional feature, that was now linearly (or at least monotonically) correlated
with the target.

According to the authors, the addition of these new features had a significant impact
on the performance of linear models.

**Example**

In the following example, we re-code 2 numerical variables using decision trees.

First we load the data and separate it into train and test:

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.discretisation import DecisionTreeDiscretiser

	# Load dataset
	data = data = pd.read_csv('houseprice.csv')

	# Separate into train and test sets
	X_train, X_test, y_train, y_test =  train_test_split(
		    data.drop(['Id', 'SalePrice'], axis=1),
		    data['SalePrice'], test_size=0.3, random_state=0)

Now we set up the discretiser. We will optimise the decision tree's depth using 3 fold
cross-validation.

.. code:: python

	# set up the discretisation transformer
	disc = DecisionTreeDiscretiser(cv=3,
                                  scoring='neg_mean_squared_error',
                                  variables=['LotArea', 'GrLivArea'],
                                  regression=True)

	# fit the transformer
	disc.fit(X_train, y_train)


With `fit()` the transformer fits a decision tree per variable. Then, we can go
ahead replace the variable values by the predictions of the trees:

.. code:: python

	# transform the data
	train_t= disc.transform(X_train)
	test_t= disc.transform(X_test)

The `binner_dict_` stores the details of each decision tree.

.. code:: python

	disc.binner_dict_


.. code:: python

	{'LotArea': GridSearchCV(cv=3, error_score='raise-deprecating',
	              estimator=DecisionTreeRegressor(criterion='mse', max_depth=None,
	                                              max_features=None,
	                                              max_leaf_nodes=None,
	                                              min_impurity_decrease=0.0,
	                                              min_impurity_split=None,
	                                              min_samples_leaf=1,
	                                              min_samples_split=2,
	                                              min_weight_fraction_leaf=0.0,
	                                              presort=False, random_state=None,
	                                              splitter='best'),
	              iid='warn', n_jobs=None, param_grid={'max_depth': [1, 2, 3, 4]},
	              pre_dispatch='2*n_jobs', refit=True, return_train_score=False,
	              scoring='neg_mean_squared_error', verbose=0),
	 'GrLivArea': GridSearchCV(cv=3, error_score='raise-deprecating',
	              estimator=DecisionTreeRegressor(criterion='mse', max_depth=None,
	                                              max_features=None,
	                                              max_leaf_nodes=None,
	                                              min_impurity_decrease=0.0,
	                                              min_impurity_split=None,
	                                              min_samples_leaf=1,
	                                              min_samples_split=2,
	                                              min_weight_fraction_leaf=0.0,
	                                              presort=False, random_state=None,
	                                              splitter='best'),
	              iid='warn', n_jobs=None, param_grid={'max_depth': [1, 2, 3, 4]},
	              pre_dispatch='2*n_jobs', refit=True, return_train_score=False,
	              scoring='neg_mean_squared_error', verbose=0)}


With tree discretisation, each bin, that is, each prediction value, does not necessarily
contain the same number of observations.

.. code:: python

	# with tree discretisation, each bin does not necessarily contain
	# the same number of observations.
	train_t.groupby('GrLivArea')['GrLivArea'].count().plot.bar()
	plt.ylabel('Number of houses')


.. image:: ../../images/treediscretisation.png

**Note**

Our implementation of the :class:`DecisionTreeDiscretiser()` will replace the original
values of the variable by the predictions of the trees. This is not strictly identical
to what the winners of the KDD competition did. They added the predictions of the features
as new variables, while keeping the original ones.

More details
^^^^^^^^^^^^

Check also for more details on how to use this transformer:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/discretisation/DecisionTreeDiscretiser.ipynb>`_
- `tree_pipe in cell 21 of this Kaggle kernel <https://www.kaggle.com/solegalli/feature-engineering-and-model-stacking>`_
.. -*- mode: rst -*-

Variable Discretisation
=======================

Feature-engine's variable discretisation transformers transform continuous numerical
variables into discrete variables. The discrete variables will contain contiguous
intervals in the case of the equal frequency and equal width transformers. The
Decision Tree discretiser will return a discrete variable, in the sense that the
new feature takes a finite number of values.

The following illustration shows the process of discretisation:

.. figure::  ../../images/Discretisation.png
   :align:   center
   :width: 500


With discretisation, sometimes we can obtain a more homogeneous value spread from an
originally skewed variable. But this is not always possible.

**Discretisation plus encoding**

Very often, after we discretise the numerical continuous variables into discrete intervals
we want to proceed their engineering as if they were categorical. This is common practice.
Throughout the user guide, we point to jupyter notebooks that showcase this functionality.

**Discretisers**

.. toctree::
   :maxdepth: 1

   EqualFrequencyDiscretiser
   EqualWidthDiscretiser
   ArbitraryDiscretiser
   DecisionTreeDiscretiser
.. _equal_freq_discretiser:

.. currentmodule:: feature_engine.discretisation

EqualFrequencyDiscretiser
=========================

The :class:`EqualFrequencyDiscretiser()` sorts continuous numerical variables into
contiguous equal frequency intervals, that is, intervals that contain approximately the
same proportion of observations. The limits of the intervals are calculated according
to percentiles or quantiles utilising `pandas.qcut()`. You decide the number of
intervals.

**A note on number of intervals**

Common values are 5 and 10. Note that if the variable is highly skewed or not continuous
smaller intervals maybe required. Otherwise, the transformer will introduce np.nan.

The :class:`EqualFrequencyDiscretiser()` works only with numerical variables. A list of
variables can be indicated, or the discretiser will automatically select all numerical
variables in the train set.

**Example**

Let's look at an example using the House Prices Dataset (more details about the
dataset :ref:`here <datasets>`).

Let's load the house prices dataset and  separate it into train and test sets:

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.discretisation import EqualFrequencyDiscretiser

	# Load dataset
	data = data = pd.read_csv('houseprice.csv')

	# Separate into train and test sets
	X_train, X_test, y_train, y_test =  train_test_split(
		    data.drop(['Id', 'SalePrice'], axis=1),
		    data['SalePrice'], test_size=0.3, random_state=0)

Now we want to discretise the 2 variables indicated below into 10 intervals of equal
number of observations:

.. code:: python

	# set up the discretisation transformer
	disc = EqualFrequencyDiscretiser(q=10, variables=['LotArea', 'GrLivArea'])

	# fit the transformer
	disc.fit(X_train)

With `fit()` the transformer learns the boundaries of each interval. Then, we can go
ahead and sort the values into the intervals:

.. code:: python

	# transform the data
	train_t= disc.transform(X_train)
	test_t= disc.transform(X_test)

The `binner_dict_` stores the interval limits identified for each variable.

.. code:: python

	disc.binner_dict_

.. code:: python

	{'LotArea': [-inf,
	  5007.1,
	  7164.6,
	  8165.700000000001,
	  8882.0,
	  9536.0,
	  10200.0,
	  11046.300000000001,
	  12166.400000000001,
	  14373.9,
	  inf],
	 'GrLivArea': [-inf,
	  912.0,
	  1069.6000000000001,
	  1211.3000000000002,
	  1344.0,
	  1479.0,
	  1603.2000000000003,
	  1716.0,
	  1893.0000000000005,
	  2166.3999999999996,
	  inf]}


With equal frequency discretisation, each bin contains approximately the same number of observations.

.. code:: python

	train_t.groupby('GrLivArea')['GrLivArea'].count().plot.bar()
	plt.ylabel('Number of houses')

We can see below that the intervals contain approximately the same number of
observations.

.. image:: ../../images/equalfrequencydiscretisation.png

|

**Discretisation plus encoding**

If we return the interval values as integers, the discretiser has the option to return
the transformed variable as integer or as object. Why would we want the transformed
variables as object?

Categorical encoders in Feature-engine are designed to work with variables of type
object by default. Thus, if you wish to encode the returned bins further, say to try and
obtain monotonic relationships between the variable and the target, you can do so
seamlessly by setting `return_object` to True. You can find an example of how to use
this functionality `here <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/discretisation/EqualFrequencyDiscretiser_plus_WoEEncoder.ipynb>`_.

More details
^^^^^^^^^^^^

Check also for more details on how to use this transformer:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/discretisation/EqualFrequencyDiscretiser.ipynb>`_
- `Jupyter notebook - Discretiser plus Weight of Evidence encoding <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/discretisation/EqualFrequencyDiscretiser_plus_WoEEncoder.ipynb>`_
.. _drop_correlated:

.. currentmodule:: feature_engine.selection

DropCorrelatedFeatures
======================

The :class:`DropCorrelatedFeatures()` finds and removes correlated variables from a dataframe.
Correlation is calculated with `pandas.corr()`. All correlation methods supported by `pandas.corr()`
can be used in the selection, including Spearman, Kendall, or Spearman. You can also pass a
bespoke correlation function, provided it returns a value between -1 and 1.

Features are removed on first found first removed basis, without any further insight. That is,
the first feature will be retained an all subsequent features that are correlated with this, will
be removed.

The transformer will examine all numerical variables automatically. Note that you could pass a
dataframe with categorical and datetime variables, and these will be ignored automatically.
Alternatively, you can pass a list with the variables you wish to evaluate.

**Example**

Let's create a toy dataframe where 4 of the features are correlated:

.. code:: python

    import pandas as pd
    from sklearn.datasets import make_classification
    from feature_engine.selection import DropCorrelatedFeatures

    # make dataframe with some correlated variables
    def make_data():
        X, y = make_classification(n_samples=1000,
                               n_features=12,
                               n_redundant=4,
                               n_clusters_per_class=1,
                               weights=[0.50],
                               class_sep=2,
                               random_state=1)

        # trasform arrays into pandas df and series
        colnames = ['var_'+str(i) for i in range(12)]
        X = pd.DataFrame(X, columns =colnames)
        return X

    X = make_data()

Now, we set up :class:`DropCorrelatedFeatures()` to find and remove variables which
(absolute) correlation coefficient is bigger than 0.8:

.. code:: python

    tr = DropCorrelatedFeatures(variables=None, method='pearson', threshold=0.8)


With `fit()` the transformer finds the correlated variables and with `transform()` it drops
them from the dataset:

.. code:: python

    Xt = tr.fit_transform(X)

The correlated feature groups are stored in the transformer's attributes:

.. code:: python

    tr.correlated_feature_sets_


.. code:: python

    [{'var_0', 'var_8'}, {'var_4', 'var_6', 'var_7', 'var_9'}]


As well as the features that will be removed from the dataset:

..  code:: python

    tr.features_to_drop_

.. code:: python

    {'var_6', 'var_7', 'var_8', 'var_9'}

If we now go ahead and print the transformed data, we see that the correlated features
have been removed.

.. code:: python

    print(print(Xt.head()))

.. code:: python

              var_0     var_1     var_2     var_3     var_4     var_5    var_10  \
    0  1.471061 -2.376400 -0.247208  1.210290 -3.247521  0.091527  2.070526
    1  1.819196  1.969326 -0.126894  0.034598 -2.910112 -0.186802  1.184820
    2  1.625024  1.499174  0.334123 -2.233844 -3.399345 -0.313881 -0.066448
    3  1.939212  0.075341  1.627132  0.943132 -4.783124 -0.468041  0.713558
    4  1.579307  0.372213  0.338141  0.951526 -3.199285  0.729005  0.398790

         var_11
    0 -1.989335
    1 -1.309524
    2 -0.852703
    3  0.484649
    4 -0.186530


More details
^^^^^^^^^^^^

In this notebook, we show how to use :class:`DropCorrelatedFeatures()` with a different
relation metric:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/selection/Drop-Correlated-Features.ipynb>`_
.. _recursive_addition:

.. currentmodule:: feature_engine.selection

RecursiveFeatureAddition
========================

:class:`RecursiveFeatureAddition` implements recursive feature addition. Recursive
feature addition (RFA) is a forward feature selection process.

This technique begins by building a model on the entire set of variables and computing
an importance score for each variable. Features are ranked by the model’s `coef_` or
`feature_importances_` attributes.

In the next step, it trains a model only using the feature with the highest importance and
stores the model performance.

Then, it adds the second most important, trains a new model and determines a new performance
metric. If the performance increases beyond the threshold, compared to the previous model,
then that feature is important and will be kept. Otherwise, that feature is removed.

It proceeds to evaluate the next most important feature, and so on, until all features
are evaluated.

Note that feature importance is used just to rank features and thus determine the order
in which the features will be added. But whether to retain a feature is determined based
on the increase in the performance of the model after the feature addition.

**Parameters**

Feature-engine's RFA has 2 parameters that need to be determined somewhat arbitrarily by
the user: the first one is the machine learning model which performance will be evaluated. The
second is the threshold in the performance increase that needs to occur, to keep a feature.

RFA is not machine learning model agnostic, this means that the feature selection depends on
the model, and different models may have different subsets of optimal features. Thus, it is
recommended that you use the machine learning model that you finally intend to build.

Regarding the threshold, this parameter needs a bit of hand tuning. Higher thresholds will
of course return fewer features.

**Example**

Let's see how to use this transformer with the diabetes dataset that comes in Scikit-learn.
First, we load the data:

.. code:: python

    import pandas as pd
    from sklearn.datasets import load_diabetes
    from sklearn.linear_model import LinearRegression
    from feature_engine.selection import RecursiveFeatureElimination

    # load dataset
    diabetes_X, diabetes_y = load_diabetes(return_X_y=True)
    X = pd.DataFrame(diabetes_X)
    y = pd.DataFrame(diabetes_y)

Now, we set up :class:`RecursiveFeatureAddition` to select features based on the r2
returned by a Linear Regression model, using 3 fold cross-validation. In this case,
we leave the parameter `threshold` to the default value which is 0.01.

.. code:: python

    # initialize linear regresion estimator
    linear_model = LinearRegression()

    # initialize feature selector
    tr = RecursiveFeatureElimination(estimator=linear_model, scoring="r2", cv=3)

With `fit()` the model finds the most useful features, that is, features that when added
cause an increase in model performance bigger than 0.01. With `transform()`, the transformer
removes the features from the dataset.

.. code:: python

    # fit transformer
    Xt = tr.fit_transform(X, y)

:class:`RecursiveFeatureAddition` stores the performance of the model trained using all
the features in its attribute:

.. code:: python

    # get the initial linear model performance, using all features
    tr.initial_model_performance_

.. code:: python

    0.488702767247119

:class:`RecursiveFeatureAddition`  also stores the change in the performance caused by
adding each feature.

..  code:: python

    # Get the performance drift of each feature
    tr.performance_drifts_

..  code:: python

    {4: 0,
     8: 0.2837159006046677,
     2: 0.1377700238871593,
     5: 0.0023329006089969906,
     3: 0.0187608758643259,
     1: 0.0027994385024313617,
     7: 0.0026951300105543807,
     6: 0.002683967832484757,
     9: 0.0003040126429713075,
     0: -0.007386876030245182}

:class:`RecursiveFeatureAddition` also stores the features that will be dropped based
n the given threshold.

..  code:: python

    # the features to drop
    tr.features_to_drop_

..  code:: python

    [0, 6, 7, 9]

If we now print the transformed data, we see that the features above were removed.

..  code:: python

    print(Xt.head())

..  code:: python

              4         8         2         3
    0 -0.044223  0.019908  0.061696  0.021872
    1 -0.008449 -0.068330 -0.051474 -0.026328
    2 -0.045599  0.002864  0.044451 -0.005671
    3  0.012191  0.022692 -0.011595 -0.036656
    4  0.003935 -0.031991 -0.036385  0.021872

.. _single_feat_performance:

.. currentmodule:: feature_engine.selection

SelectBySingleFeaturePerformance
================================

The :class:`SelectBySingleFeaturePerformance()` selects features based on the performance of
machine learning models trained using individual features. That is, it selects
features based on their individual performance. In short, the selection algorithms works
as follows:

1. Train a machine learning model per feature (using only 1 feature)
2. Determine the performance metric of choice
3. Retain features which performance is above a threshold

If the parameter `threshold` is left to None, it will select features which performance is
above the mean performance of all features.

**Example**

Let's see how to use this transformer with the diabetes dataset that comes in Scikit-learn.
First, we load the data:

.. code:: python

    import pandas as pd
    from sklearn.datasets import load_diabetes
    from sklearn.linear_model import LinearRegression
    from feature_engine.selection import SelectBySingleFeaturePerformance

    # load dataset
    diabetes_X, diabetes_y = load_diabetes(return_X_y=True)
    X = pd.DataFrame(diabetes_X)
    y = pd.DataFrame(diabetes_y)

Now, we start :class:`SelectBySingleFeaturePerformance()` to select features based on the
r2 returned by a Linear regression, using 3 fold cross-validation. We want to select features
which r2 > 0.01.

.. code:: python

    # initialize feature selector
    sel = SelectBySingleFeaturePerformance(
            estimator=LinearRegression(), scoring="r2", cv=3, threshold=0.01)

With `fit()` the transformer fits 1 model per feature, determines the performance and
selects the important features:

.. code:: python

    # fit transformer
    sel.fit(X, y)

The features that will be dropped are stored in an attribute:

.. code:: python

    sel.features_to_drop_

.. code:: python

    [1]

:class:`SelectBySingleFeaturePerformance()` also stores the performace of each one of the
models, in case we want to study those further:

..  code:: python

    sel.feature_performance_

.. code:: python

    {0: 0.029231969375784466,
     1: -0.003738551760264386,
     2: 0.336620809987693,
     3: 0.19219056680145055,
     4: 0.037115559827549806,
     5: 0.017854228256932614,
     6: 0.15153886177526896,
     7: 0.17721609966501747,
     8: 0.3149462084418813,
     9: 0.13876602125792703}

With `transform()` we go ahead and remove the features from the dataset:

.. code:: python

    # drop variables
    Xt = sel.transform(X)


More details
^^^^^^^^^^^^

Check also:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/selection/Select-by-Single-Feature-Performance.ipynb>`_

.. _feature_shuffling:

.. currentmodule:: feature_engine.selection

SelectByShuffling
=================

The :class:`SelectByShuffling()` selects important features if a random permutation of
their values decreases the model performance. If the feature is predictive, a random
shuffle of the values across the rows, should return predictions that are off the truth.
If the feature is not predictive, their values should have a minimal impact on the prediction.

Procedure
---------

The algorithm works as follows:

1. Train a machine learning model using all features
2. Determine a model performance metric of choice
3. Shuffle the order of 1 feature values
4. Use the model trained in 1 to obtain new predictions
5. Determine the performance with the predictions in 4
6. If there is a drop in performance beyond a threshold, keep the feature.
7. Repeat 3-6 until all features are examined.

**Example**

Let's see how to use this transformer with the diabetes dataset that comes in Scikit-learn.
First, we load the data:

.. code:: python

    import pandas as pd
    from sklearn.datasets import load_diabetes
    from sklearn.linear_model import LinearRegression
    from feature_engine.selection import SelectByShuffling

    # load dataset
    diabetes_X, diabetes_y = load_diabetes(return_X_y=True)
    X = pd.DataFrame(diabetes_X)
    y = pd.DataFrame(diabetes_y)

Now, we set up the model for which we want to have the performance drop evaluated:

.. code:: python

    # initialize linear regresion estimator
    linear_model = LinearRegression()

Now, we instantiate :class:`SelectByShuffling()` to select features by shuffling, based on
the r2 of the model from the previous cell, using 3 fold cross-validation. The parameter
`threshold` was left to None, which means that features will be selected if the performance
drop is bigger than the mean drop caused by all features.

.. code:: python

    # initialize feature selector
    tr = SelectByShuffling(estimator=linear_model, scoring="r2", cv=3)


With `fit()` the transformer finds the important variables, that is, those which values
permutations caused a drop in the model performance. With `transform()` it drops them
from the dataset:

.. code:: python

    # fit transformer
    Xt = tr.fit_transform(X, y)

:class:`SelectByShuffling()` stores the performance of the model trained using all the features
in its attribute:

.. code:: python

    tr.initial_model_performance_

.. code:: python

    0.488702767247119

:class:`SelectByShuffling()` also stores the performance change caused by every single
feature after shuffling. In case you are not satisfied with the threshold used, you can get
an idea of where the threshold could be by looking at these values:

..  code:: python

    tr.performance_drifts_

.. code:: python

    {0: -0.02368121940502793,
     1: 0.017909161264480666,
     2: 0.18565460365508413,
     3: 0.07655405817715671,
     4: 0.4327180164470878,
     5: 0.16394693824418372,
     6: -0.012876023845921625,
     7: 0.01048781540981647,
     8: 0.3921465005640224,
     9: -0.01427065640301245}

:class:`SelectByShuffling()` also stores the features that will be dropped based on the
threshold indicated.

.. code:: python

    tr.features_to_drop_

.. code:: python

    [0, 1, 3, 6, 7, 9]

.. _drop_features:

.. currentmodule:: feature_engine.selection

DropFeatures
=============

The :class:`DropFeatures()` drops a list of variables indicated by the user from the original
dataframe. The user can pass a single variable as a string or list of variables to be
dropped.

:class:`DropFeatures()` offers similar functionality to `pandas.dataframe.drop <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.drop.html>`_,
but the difference is that :class:`DropFeatures()` can be integrated into a Scikit-learn
pipeline.


**When is this transformer useful?**

Sometimes, we create new variables combining other variables in the dataset, for
example, we obtain the variable `age` by subtracting `date_of_application` from
`date_of_birth`. After we obtained our new variable, we do not need the date
variables in the dataset any more. Thus, we can add :class:`DropFeatures()` in the Pipeline
to have these removed.

**Example**

Let's see how to use :class:`DropFeatures()` in an example with the Titanic dataset. We
first load the data and separate it into train and test:

.. code:: python

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.model_selection import train_test_split

    from feature_engine.selection import DropFeatures

    # Load dataset
    def load_titanic():
            data = pd.read_csv('https://www.openml.org/data/get_csv/16826755/phpMYEkMl')
            data = data.replace('?', np.nan)
            data['cabin'] = data['cabin'].astype(str).str[0]
            data['pclass'] = data['pclass'].astype('O')
            data['embarked'].fillna('C', inplace=True)
            return data

    # load data as pandas dataframe
    data = load_titanic()

    # Separate into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
                data.drop(['survived', 'name'], axis=1),
                data['survived'], test_size=0.3, random_state=0)

Now, we go ahead and print the dataset column names:

.. code:: python

    # original columns
    X_train.columns

.. code:: python

    Index(['pclass', 'sex', 'age', 'sibsp', 'parch', 'ticket', 'fare', 'cabin',
           'embarked', 'boat', 'body', 'home.dest'],
          dtype='object')

Now, with :class:`DropFeatures()` we can very easily drop a group of variables. Below
we set up the transformer to drop a list of 6 variables:

.. code:: python

    # set up the transformer
    transformer = DropFeatures(
        features_to_drop=['sibsp', 'parch', 'ticket', 'fare', 'body', 'home.dest']
    )

    # fit the transformer
    transformer.fit(X_train)

With `fit()` this transformer does not learn any parameter. We can go ahead and remove
the variables as follows:

.. code:: python

    # transform the data
    train_t = transformer.transform(X_train)

And now, if we print the variable names of the transformed dataset, we see that it has
been reduced:

.. code:: python

    train_t.columns

.. code:: python

    Index(['pclass', 'sex', 'age', 'cabin', 'embarked' 'boat'],
          dtype='object')


More details
^^^^^^^^^^^^

In this Kaggle kernel we feature 3 different end-to-end machine learning pipelines using
:class:`DropFeatures()`:

- `Kaggle Kernel <https://www.kaggle.com/solegalli/feature-engineering-and-model-stacking>`_
6).. _psi_selection:

.. currentmodule:: feature_engine.selection

DropHighPSIFeatures
===================

The :class:`DropHighPSIFeatures()` finds and removes features with changes in their
distribution, i.e. "unstable values", from a pandas dataframe.
The stability of the distribution is computed using the **Population Stability
Index (PSI)** and all features having a PSI value above a given threshold are removed.

Unstable features may introduce an additional bias in a model if the training population
significantly differs from the population in production. Removing features for
which a shift in the distribution is suspected leads to
more robust models and therefore to better performance. In the field of Credit Risk
modelling, eliminating features with high PSI is common practice and usually required by the
Regulator.

Population Stability Index - PSI
--------------------------------

The PSI is a measure of how much a population has changed in time or how different the distributions
are between two different population samples.

To determine the PSI, continuous features are sorted into discrete intervals, the
fraction of observations per interval is then determined, and finally those values
are compared between the 2 groups, or as we call them in Feature-engine, between the
basis and test sets, to obtain the PSI.

In other words, the PSI is computed as follows:

- Define the intervals into which the observations will be sorted.
- Sort the feature values into those intervals.
- Determine the fraction of observations within each interval.
- Compute the PSI.

The PSI is determined as:

.. math::

    PSI = \sum_{i=1}^n (test_i - basis_i) . ln(\frac{test_i}{basis_i})

where `basis` and `test` are the "reference" and "evaluation" datasets, respectively, and `i`
refers to the interval.

In other words, the PSI determines the difference in the proportion of observations in each
interval, between the reference (aka, original) and test datasets.

In the PSI equation, `n` is the total number of intervals.

Important
~~~~~~~~~

When working with the PSI it is worth highlighting the following:

- The PSI is not symmetric; switching the order of the basis and test dataframes in the PSI calculation will lead to different values.
- The number of bins used to define the distributions has an impact on the PSI values.
- The PSI is a suitable metric for numerical features (i.e., either continuous or with high cardinality).
- For categorical or discrete features, the change in distributions is better assessed with Chi-squared.

Threshold
~~~~~~~~~

Different thresholds can be used to assess the magnitude of the distribution shift according
to the PSI value. The most commonly used thresholds are:

- Below 10%, the variable has not experienced a significant shift.
- Above 25%, the variable has experienced a major shift.
- Between those two values, the shift is intermediate.


Procedure
---------

To compute the PSI, the :class:`DropHighPSIFeatures()` splits the input dataset in
two: a basis data set (aka the reference data) and a test set. The basis data set is assumed to contain
the expected or original feature distributions. The test set will be assessed
against the basis data set.

In the next step, the interval boundaries are determined based on the features in the basis
or reference data. These intervals can be determined to be of equal with, or equal number
of observations.

Next, :class:`DropHighPSIFeatures()` sorts each of the variable values into those intervals, both in the
basis and test datasets, and then determines the proportion (percentage) of observations
within each interval.

Finally, the PSI is determined as indicated in the previous paragraph for each feature.
With the PSI value per feature, :class:`DropHighPSIFeatures()` can now select the features that are unstable and
drop them, based on a threshold.


Splitting the data
------------------

:class:`DropHighPSIFeatures()` allows us to determine how much a feature distribution has
changed in time, or how much it differs between 2 groups.

If we want to evaluate the distribution change in time, we can use a datetime variable as splitting
reference and provide a datetime cut-off as split point.

If we want to compare the distribution change between 2 groups, :class:`DropHighPSIFeatures()`
offers 3 different approaches to split the input dataframe:

- Based on proportion of observations.
- Based on proportions of unique observations.
- Using a cut-off value.


Proportion of observations
~~~~~~~~~~~~~~~~~~~~~~~~~~

Splitting by proportion of observations will result in a certain proportion of observations
allocated to either the reference and test datasets. For example, if we set `split_frac=0.75`,
then 75% and 25% of the observations will be put into the reference and test data, respectively.

If we select this method, we can pass a variable in the parameter `split_col` or leave it to None.

Note that the data split is not done at random, but instead guided by the values in the reference
variable indicated in `split_col`. Under the hood, the reference variable indicated in `split_col`
is ordered, and the percentage of observations is determined with NumPy quantile. This means
that the observations with smaller values in `split_col` will land in the reference dataset,
and those with bigger values will go to the test set.

If the rows in your dataset are sorted in time, this could be a good default option to split the
dataframe in 2 and compute the PSI. This will for example be the case if your data set contains
daily (or any other frequency) sales information on a company's products.

Proportions of unique observations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If we split based on proportion of unique observations, it is important that we indicate which
column we want to use as reference in the `split_col` parameter, to make a meaningful split. If
we leave this to None, :class:`DropHighPSIFeatures()` will use the dataframe index as reference.
This makes sense only if the index in the dataframe has meaningful values.

:class:`DropHighPSIFeatures()` will first identify the unique values of the variable in
`split_col`. Then it will put a certain proportion of those values into the reference
dataset and the remaining to the test dataset. The proportion is indicated in the parameter
`split_frac`.

Under the hood, :class:`DropHighPSIFeatures()` will sort the unique values of the reference
variable, and then use NumPy quantiles to determine the fraction that should be allocated to the
reference and test sets. Thus, it is important to consider that the order of the unique values
matters in the split.

This split makes sense when we have for example unique customer identifiers and multiple rows
per customer in the dataset. We want to make sure that all rows belonging to the same customer
are allocated either in the reference or test data, but the same customer cannot be in both
data sets. This way of splitting the data will also ensure that we have a certain percentage,
indicated in `split_frac` of customers in either data set after the split.

Thus, if `split_frac=0.6` and `split_distinct=True`, :class:`DropHighPSIFeatures()` will send
the first 60% of customers to the reference data set, and the remaining 40% to the test set. And it will
ensure that rows belonging to the same customer are just in one of the 2 data sets.

Using a cut-off value
~~~~~~~~~~~~~~~~~~~~~

We have the option to pass a reference variable to use to split the dataframe using `split_col` and
also a cut-off value in the `cut_off` parameter. The cut-off value can be a number, integer or float,
a date or a list of values.

If we pass a datetime column in `split_col` and a datetime value in the `cut_off`, we can split the
data in a temporal manner. Observations collected before the time indicated will be sent to the reference
dataframe, and the remaining to the test set.

If we pass a list of values in the `cut_off` all observations which values are included in the
list will go into the reference data set, and the remaining to the test set. This split is useful
if we have a categorical variable indicating a portfolio from which the observations have been collected.
For example, if we set `split_col='portfolio'` and `cut_off=['port_1', 'port_2']`, all observations
that belong to the first and second portfolio will be sent to the reference data set, and the observations
from other portfolios to the test set.

Finally, if we pass a number to `cut_off`, all observations which value in the variable indicated
in `split_col` is <= cut-off, will be sent to the reference data set, alternatively to the test set.
This can be useful for example when dates are defined as integer (for example 20200411) or when
using an ordinal customer segmentation to split the dataframe (1: retail customers, 2: private
banking customers, 3: SME and 4: Wholesale).

split_col
~~~~~~~~~

To split the data set, we recommend that you indicate which column you want to use as
reference in the `split_col` parameter. If you don't, the split will be done based on the
values of the dataframe index. This might be a good option if the index contains meaningful
values or if splitting just based on `split_frac`.


Examples
--------

The versatility of the class lies in the different options to split the input dataframe
in a reference or basis data set with the "expected" distributions, and a test set which
will be evaluated against the reference.

After splitting the data, :class:`DropHighPSIFeatures()` goes ahead and compares the
feature distributions in both data sets by computing the PSI.

To illustrate how to best use :class:`DropHighPSIFeatures()` depending on your data, we
provide various examples illustrating the different possibilities.

Case 1: split data based on proportions (split_frac)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this case, :class:`DropHighPSIFeatures()` will split the dataset in 2, based on the
indicated proportion. The proportion is indicated in the `split_frac` parameter. You have
the option to select a variable in `split_col` or leave it to None. In the latter, the
dataframe index will be used to split.

Let's first create a toy dataframe containing 5 random variables and 1 variable
with a shift in its distribution (*var_3* in this case).

.. code:: python

    import pandas as pd
    import seaborn as sns

    from sklearn.datasets import make_classification
    from feature_engine.selection import DropHighPSIFeatures

    # Create a dataframe with 500 observations and 6 random variables
    X, y = make_classification(
        n_samples=500,
        n_features=6,
        random_state=0
    )

    colnames = ["var_" + str(i) for i in range(6)]
    X = pd.DataFrame(X, columns=colnames)

    # Add a column with a shift.
    X['var_3'][250:] = X['var_3'][250:] + 1

The default approach in :class:`DropHighPSIFeatures()` is to split the
input dataframe `X` in two equally sized data sets. You can adjust the proportions by changing
the value in the `split_frac` parameter.

For example, let's split the input dataframe into a reference data set containing 60% of
the observations and a test set containing 40% of the observations.

.. code:: python

    # Remove the features with high PSI values using a 60-40 split.

    transformer = DropHighPSIFeatures(split_frac=0.6)
    transformer.fit(X)

The value of `split_frac` tells :class:`DropHighPSIFeatures()` to split X according to a
60% - 40% ratio. The `fit()` method performs the split of the dataframe and the calculation
of the PSI.

Because we created random variables, these features will have low PSI values (i.e., no
distribution change). However, we manually added a distribution shift in the variable *var_3*
and therefore expect the PSI for this particular feature to be above the
0.25 PSI threshold.

The PSI values are accessible through the `psi_values_` attribute:

.. code:: python

    transformer.psi_values_

The analysis of the PSI values below shows that only feature 3 (called `var_3`)
has a PSI above the 0.25 threshold (default value) and will be removed
by the `transform` method.

.. code:: python

    {'var_0': 0.07405459925568803,
    'var_1': 0.09124093185820083,
    'var_2': 0.16985790067687764,
    'var_3': 1.342485289730313,
    'var_4': 0.0743442762545251,
    'var_5': 0.06809060587241555}

From the output, we see that the PSI value for *var_0* is around 7%. This means
that, when comparing the first 300 and the last 200 observations of the dataframe,
there is only a small difference in the distribution of the *var_0* feature.
A similar conclusion applies to *var_1, var_2, var_4* and *var_5*.
Looking at the PSI value for *var_3*, we see that it exceeds by far the 0.25
threshold. We can then conclude the population of this feature has shifted and
it is wise not to include it in the feature set for modelling.

The cut-off value used to split the dataframe is stored in the `cut_off_` attribute:

.. code:: python

    transformer.cut_off_

This yields the following answer

.. code:: python

    299.4

The value of 299.4 means that observations with index from 0 to 299 are used
to define the basis data set. This corresponds to 60% (300 / 500) of the original dataframe
(X).
The value of 299.4 may seem strange because it is not one of the value present in (the
(index of) the dataframe. Intuitively, we would expect the cut_off to be an integer
in the present case. However, the cut_off is computed using quantiles and the quantiles
are computed using extrapolation.

Splitting with proportions will order the index or the reference column first, and then
determine the data that will go into each dataframe. In other words, the order of the index
or the variable indicated in `split_col` matters. Observations with the lowest values will
be sent to the basis dataframe and the ones with the highest values to the test set.

The `features_to_drop_` attribute provides the list with the features to
be dropped when executing the `transform` method.

The command

.. code:: python

    transformer.features_to_drop_

Yields the following result:

.. code:: python

    ['var_3']

That the *var_3* feature is dropped during the procedure is illustrated when
looking at the columns from the `X_transformed` dataframe.

.. code:: python

    X_transformed = transformer.transform(X)

    X_transformed.columns

    Index(['var_0', 'var_1', 'var_2', 'var_4', 'var_5'], dtype='object')

:class:`DropHighPSIFeatures()` also contains a `fit_transform` method that combines
the `fit` and the `transform` methods.

The difference in distribution between a non-shifted and
a shifted distribution is clearly visible when plotting the cumulative density
function.

For the shifted variable:

.. code:: python

    X['above_cut_off'] = X.index > transformer.cut_off_
    sns.ecdfplot(data=X, x='var_3', hue='above_cut_off')

and a non-shifted variable (for example *var_1*)

.. code:: python

    sns.ecdfplot(data=X, x='var_1', hue='above_cut_off')


.. image:: ../../images/PSI_distribution_case1.png


Case 2: split data based on variable (numerical cut_off)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


In the previous example, we wanted to split the input dataframe in 2 datasets, with the
reference dataset containing 60% of the observations. We let :class:`DropHighPSIFeatures()`
find the cut-off to achieve this.

We can instead, provide ourselves the numerical cut-off that determines which observations will
go to the reference or basis data set, and which to the test set. Using the `cut_off` parameter,
we can define the specific threshold for the split.

A real life example for this case is the use of the customer ID or contract ID
to split the dataframe. These IDs are often increasing in value over time which justifies
their use to assess distribution shifts in the features.

Let's create a toy dataframe representing the customers' characteristics of a
company. This dataset contains six random variables (in
real life this are variables like age or postal code), the seniority of the customer
(i.e. the number of months since the start of the relationship between the
customer and the company) and the customer ID (i.e. the number (integer) used
to identify the customer). Generally the customer ID grows
over time which means that early customers have a lower customer ID than late
customers.

From the definition of the variables, we expect the *seniority* to increase with
the customer ID and therefore to have a high PSI value when comparing early and
late customer,

.. code:: python

    import pandas as pd
    from sklearn.datasets import make_classification
    from feature_engine.selection import DropHighPSIFeatures

    X, y = make_classification(
            n_samples=500,
            n_features=6,
            random_state=0
        )

    colnames = ["var_" + str(i) for i in range(6)]
    X = pd.DataFrame(X, columns=colnames)

    # Let's add a variable for the customer ID
    X['customer_id'] = [customer_id for customer_id in range(1, 501)]

    # Add a column with the seniority... that is related to the customer ID
    X['seniority'] = 100 - X['customer_id'] // 10

    transformer = DropHighPSIFeatures(split_col='customer_id', cut_off=250)
    transformer.fit(X)

In this case, :class:`DropHighPSIFeatures()` will allocate in the basis or reference data
set, all observations which values in `customer_id` are <= 250. The test dataframe contains the
remaining observations.

The method `fit()` will determine the PSI values, which are stored in the class:

.. code:: python

    transformer.psi_values_

We see that :class:`DropHighPSIFeatures()` does not provide any PSI value for
the `customer_id` feature, because this variable was used as a reference to split the data.

.. code:: python

    {'var_0': 0.07385590683974477,
    'var_1': 0.061155637727757485,
    'var_2': 0.1736694458621651,
    'var_3': 0.044965387331530465,
    'var_4': 0.0904519893659045,
    'var_5': 0.027545195437270797,
    'seniority': 7.8688986006052035}


.. code:: python

    transformer.features_to_drop_

Gives

.. code:: python

    ['seniority']


Executing the dataframe transformation leads to the exclusion of the *seniority*
feature but not to the exclusion of the *customer_id*.

.. code:: python

    X_transformed = transformer.transform(X)

    X_transformed.columns

    Index(['var_0', 'var_1', 'var_2', 'var_3', 'var_4', 'var_5', 'customer_id'], dtype='object')


Case 3: split data based on time (date as cut_off)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`DropHighPSIFeatures()` can handle different types of `split_col`
variables. The following case illustrates how it works with a date variable. In fact,
we often want to determine if the distribution of a feature changes in time, for example
after a certain event like the start of the Covid-19 pandemic.

This is how to do it. Let's create a toy dataframe with 6 random numerical variables
and two date variables. One will be use to specific the split of the dataframe
while the second one is expected to have a high PSI value.

.. code:: python

    import pandas as pd
    from datetime import date
    from sklearn.datasets import make_classification
    from feature_engine.selection import DropHighPSIFeatures

    X, y = make_classification(
            n_samples=1000,
            n_features=6,
            random_state=0
        )

    colnames = ["var_" + str(i) for i in range(6)]
    X = pd.DataFrame(X, columns=colnames)

    # Add two time variables to the dataframe
    X['time'] = [date(year, 1, 1) for year in range(1000, 2000)]
    X['century'] = X['time'].apply(lambda x: ((x.year - 1) // 100) + 1)

    # Let's shuffle the dataframe and reset the index to remove the correlation
    # between the index and the time variables.

    X = X.sample(frac=1).reset_index(drop=True)

Dropping features with high PSI values comparing two periods of time is done simply
by providing the name of the column with the date and a cut-off date.
In the example below the PSI calculations
will be done comparing the periods up to the French revolution and after.

.. code:: python

    transformer = DropHighPSIFeatures(split_col='time', cut_off=date(1789, 7, 14))
    transformer.fit(X)

**Important**: if the date variable is in pandas or NumPy datetime format, you may need
to pass the cut_off value as `pd.to_datetime(1789-07-14)`.

The PSI values shows the *century* variables in unstable as its value is above
the 0.25 threshold.

.. code:: python

    transformer.psi_values_

    {'var_0': 0.0181623637463045,
    'var_1': 0.10595496570984747,
    'var_2': 0.05425659114295842,
    'var_3': 0.09720689210928271,
    'var_4': 0.07917647542638032,
    'var_5': 0.10122468631060424,
    'century': 8.272395772368412}

The class has correctly identified the feature to be dropped.

.. code:: python

    transformer.features_to_drop_

    ['century']

And the transform method correctly removes the feature.

.. code:: python

    X_transformer = transformer.transform(X)

    X_transformed.columns

    Index(['var_0', 'var_1', 'var_2', 'var_3', 'var_4', 'var_5', 'time'], dtype='object')

The difference in distribution between a non-shifted and
a shifted distribution is clearly visible when plotting the cumulative density
function for each of the group.

We can plot the cumulative distribution of the shifted variable like this:

.. code:: python

    X['above_cut_off'] = X.time > pd.to_datetime(transformer.cut_off_)
    sns.ecdfplot(data=X, x='century', hue='above_cut_off')

and the distribution of a non-shifted variable, for example *var_2*, like this:

.. code:: python

    sns.ecdfplot(data=X, x='var_2', hue='above_cut_off')

And below we can compare both plots:

.. image:: ../../images/PSI_distribution_case3.png


Case 4: split data based on a categorical variable (category or list as cut_off)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`DropHighPSIFeatures()` can also split the original dataframe based on
a categorical variable. The cut-off can then be defined in two ways:

- Using a single string.
- Using a list of values.

In the first case, the column with the categorical variable is sorted alphabetically and
the split is determined by the cut-off. We recommend being very careful when using a single
category as cut-off, because alphabetical sorting in combination with a cut-off does
not always provide obvious results. In other words, for this way of splitting the data to
be meaningful, the alphabetical order of the categories in the reference variable should have
an intrinsic meaning.

A better purpose for splitting the data based on a categorical variable would be to pass a
list with the values of the variable that want in the reference dataframe. A real life
example for this case is the computation of the PSI between different customer segments
like 'Retail', 'SME' or 'Wholesale'. In this case, if we indicate ['Retail'] as
cut-off, observations for Retail will be sent to the basis data set, and those for 'SME'
and 'Wholesale' will be added to the test set.

Split passing a category value
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's show how to set up the transformer in this case. The example data set
contains 6 randoms variables, a categorical variable with the labels of the
different categories and 2 category related features.

.. code:: python

    import pandas as pd
    import seaborn as sns

    from sklearn.datasets import make_classification
    from feature_engine.selection import DropHighPSIFeatures

    X, y = make_classification(
        n_samples=1000,
        n_features=6,
        random_state=0
    )

    colnames = ["var_" + str(i) for i in range(6)]
    X = pd.DataFrame(X, columns=colnames)

    # Add a categorical column
    X['group'] = ["A", "B", "C", "D", "E"] * 200

    # And two category related features
    X['group_means'] = X.group.map({"A": 1, "B": 2, "C": 0, "D": 1.5, "E": 2.5})
    X['shifted_feature'] = X['group_means'] + X['var_2']

We can define a simple cut-off value (for example the letter C). In this case, observations
with values that come before C, alphabetically, will be allocated to the reference data set.

.. code:: python

    transformer = DropHighPSIFeatures(split_col='group', cut_off='C')
    X_transformed = transformer.fit_transform(X)

The PSI values are provided in the `psi_values_` attribute.

.. code:: python

    transformer.psi_values_

    {'var_0': 0.06485778974895254,
    'var_1': 0.03605540598761757,
    'var_2': 0.040632784917352296,
    'var_3': 0.023845405645510645,
    'var_4': 0.028007185972248064,
    'var_5': 0.07009152672971862,
    'group_means': 6.601444547497699,
    'shifted_feature': 0.48428009522119164}

From these values we see that the last 2 features should be removed. We can corroborate
that in the `features_to_drop_` attribute:


.. code:: python

    transformer.features_to_drop_

    ['group_means', 'shifted_feature']

And these columns are removed from the original dataframe by the transform
method that, in the present case, has been applied through the fit_transform
method a couple of block cells above.


.. code:: python

    X_transformed.columns

    Index(['var_0', 'var_1', 'var_2', 'var_3', 'var_4', 'var_5', 'group'], dtype='object')


Split passing a list of categories
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instead of passing a category value, we can instead pass a list of values to the `cut_off`.
Using the same data set let's set up the :class:`DropHighPSIFeatures()` to split
the dataframe according to the list ['A', 'C', 'E'] for the categorical variable *group*.

In this case, the PSI's will be computed by comparing two dataframes: the first one
containing only the values A, C and E for the *group* variable and the second one containing
only the values B and D.

.. code:: python

    trans = DropHighPSIFeatures(split_col='group', cut_off=['A', 'C', 'E'])
    X_no_drift = trans.fit_transform(X)


.. code:: python

    trans.psi_values_

    'var_0': 0.04322345673014104,
    'var_1': 0.03534439253617049,
    'var_2': 0.05220272785661243,
    'var_3': 0.04550964862452317,
    'var_4': 0.04492720670343145,
    'var_5': 0.044886435640028144,
    'group_means': 6.601444547497699,
    'shifted_features': 0.3683642099948127}


Here again, the object will remove the *group_means* and the *shifted_features* columns
from the dataframe.


.. code:: python

    trans.features_to_drop_

    ['group_means', 'shifted_features']

And these columns are removed from the original dataframe by the transform
method that has been applied through the `fit_transform` method.


.. code:: python

    X_transformed.columns

    Index(['var_0', 'var_1', 'var_2', 'var_3', 'var_4', 'var_5', 'group'], dtype='object')


In the following plots, we can compare the distribution of a feature with high PSI and
one with low PSI, in the different categories of the categorical variable.

With this code we plot the cumulative distribution of a feature which distribution is
different among the different categories of the variable:

.. code:: python

    sns.ecdfplot(data=X, x='shifted_feature', hue='group')

With this code we plot the cumulative distribution of a feature which distribution is
the same across the different categories of the categorical variable:

.. code:: python

    sns.ecdfplot(data=X, x='var_0', hue='group')

And below we can compare the plots of both features:

.. image:: ../../images/PSI_distribution_case4.png

Case 5: split data based on unique values (split_distinct)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A variant to the previous example is the use of the `split_distinct` functionality.
In this case, the split is not done based on the number observations from
`split_col` but from the number of distinct values in the reference variable indicated
in `split_col`.

A real life example for this case is when dealing with groups of different sizes
like customers income classes ('1000', '2000', '3000', '4000', ...).
Split_distinct allows to control the numbers of classes in the basis and test
dataframes regardless of the number of observations in each class.

This case is illustrated in the toy data for this case. The data set contains
6 random variable and 1 income variable that is larger for one of the 6 group
defined (the F group).

.. code:: python

    import numpy as np
    import pandas as pd
    import seaborn as sns

    from sklearn.datasets import make_classification
    from feature_engine.selection import DropHighPSIFeatures

    X, y = make_classification(
        n_samples=1000,
        n_features=6,
        random_state=0
    )

    colnames = ["var_" + str(i) for i in range(6)]
    X = pd.DataFrame(X, columns=colnames)

    # Add a categorical column
    X['group'] = ["A", "B", "C", "D", "E"] * 100 + ["F"] * 500

    # And an income variable that is category dependent.
    np.random.seed(0)
    X['income'] = np.random.uniform(1000, 2000, 500).tolist() +
                  np.random.uniform(1250, 2250, 500).tolist()

    # Shuffle the dataframe to make the dataset more real life case.
    X = X.sample(frac=1).reset_index(drop=True)


The `group` column contains 500 observations in the (A, B, C, D, E)
group and 500 in the (F) group.

When we pass `split_distinct=True` when initializing
the `DropHighPSIFeatures` object, the two dataframes used to compute the
PSI will contain the same number of **unique** values in the `group`
column (i.e., one dataframe will contain 300 rows associated to groups A, B and C
while the other will contain 700 rows associated to groups D, E and F).

.. code:: python

    transformer = DropHighPSIFeatures(split_col='group', split_distinct=True)
    transformer.fit(X)

    transformer.psi_values_

This yields the following PSI values:

.. code:: python

    {'var_0': 0.014825303242393804,
    'var_1': 0.03818316821350485,
    'var_2': 0.029635981271458896,
    'var_3': 0.021700399485890084,
    'var_4': 0.061194837255216114,
    'var_5': 0.04119583769297253,
    'income': 0.46191580731264914}

And we can find the feature that will be dropped, income, here:

.. code:: python

    transformer.features_to_drop_

        ['income']


The former feature will be removed from the dataset when calling the `transform()` method.


.. code:: python

    X_transformed = transformer.transform(X)

    X_transformed.columns

    Index(['var_0', 'var_1', 'var_2', 'var_3', 'var_4', 'var_5', 'group'], dtype='object')



The difference in distribution between a non-shifted and
a shifted distribution is clearly visible when plotting the cumulative density
function for each of the group.

For the shifted variable (income):

.. code:: python

    sns.ecdfplot(data=X, x='income', hue='group')

and a non-shifted variable (for example *var_4*)

.. code:: python

    sns.ecdfplot(data=X, x="var_4", hue="group")


.. image:: ../../images/PSI_distribution_case5.png


More details
^^^^^^^^^^^^

In this notebook, we show how to use :class:`DropHighPSIFeatures` on a real dataset and
give more detail about the underlying base and reference sub-dataframes used to
determine the PSI.

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/selection/Drop-High-PSI-Features.ipynb>`_

.. _drop_constant:

.. currentmodule:: feature_engine.selection

DropConstantFeatures
====================


The :class:`DropConstantFeatures()` drops constant and quasi-constant variables from a dataframe.
By default, it drops only constant variables. Constant variables have a single
value. Quasi-constant variables have a single value in most of its observations.

This transformer works with numerical and categorical variables, and it offers a pretty straightforward
way of reducing the feature space. Be mindful though, that depending on the context, quasi-constant
variables could be useful.

**Example**

Let's see how to use :class:`DropConstantFeatures()` in an example with the Titanic dataset. We
first load the data and separate it into train and test:

.. code:: python

    import numpy as np
    import pandas as pd
    from sklearn.model_selection import train_test_split

    from feature_engine.selection import DropConstantFeatures

    # Load dataset
    def load_titanic():
            data = pd.read_csv('https://www.openml.org/data/get_csv/16826755/phpMYEkMl')
            data = data.replace('?', np.nan)
            data['cabin'] = data['cabin'].astype(str).str[0]
            data['pclass'] = data['pclass'].astype('O')
            data['embarked'].fillna('C', inplace=True)
            return data

    # load data as pandas dataframe
    data = load_titanic()

    # Separate into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
                data.drop(['survived', 'name', 'ticket'], axis=1),
                data['survived'], test_size=0.3, random_state=0)

Now, we set up the :class:`DropConstantFeatures()` to remove features that show the same
value in more than 70% of the observations:

.. code:: python

    # set up the transformer
    transformer = DropConstantFeatures(tol=0.7, missing_values='ignore')


With `fit()` the transformer finds the variables to drop:

.. code:: python

    # fit the transformer
    transformer.fit(X_train)

The variables to drop are stored in the attribute `features_to_drop_`:

.. code:: python

    transformer.features_to_drop_

.. code:: python

    ['parch', 'cabin', 'embarked']


We see in the following code snippets that for the variables parch and embarked, more
than 70% of the observations displayed the same value:

.. code:: python

    X_train['embarked'].value_counts() / len(X_train)

.. code:: python

    S    0.711790
    C    0.197598
    Q    0.090611
    Name: embarked, dtype: float64


71% of the passengers embarked in S.

.. code:: python

    X_train['parch'].value_counts() / len(X_train)

.. code:: python

    0    0.771834
    1    0.125546
    2    0.086245
    3    0.005459
    4    0.004367
    5    0.003275
    6    0.002183
    9    0.001092
    Name: parch, dtype: float64

77% of the passengers had 0 parent or child. Because of this, these features were
deemed constant and removed.

With `transform()`, we can go ahead and drop the variables from the data:

.. code:: python

    # transform the data
    train_t = transformer.transform(X_train)

More details
^^^^^^^^^^^^

In this Kaggle kernel we use :class:`DropConstantFeatures()` together with other feature selection algorithms:

- `Kaggle kernel <https://www.kaggle.com/solegalli/feature-selection-with-feature-engine>`_
.. _smart_correlation:

.. currentmodule:: feature_engine.selection

SmartCorrelatedSelection
========================

When we have big datasets, more than 2 features can be correlated. We could have 3, 4 or
more features that are correlated. Thus, which one should be keep and which ones should we
drop?

:class:`SmartCorrelatedSelection` tries to answer this question.

From a group of correlated variables, the :class:`SmartCorrelatedSelection` will retain
the one with:

- the highest variance
- the highest cardinality
- the least missing data
- the most important (based on embedded selection methods)

And drop the rest.

Features with higher diversity of values (higher variance or cardinality), tend to be more
predictive, whereas features with least missing data, tend to be more useful.

Procedure
---------

:class:`SmartCorrelatedSelection` will first find correlated feature groups using any
correlation method supported by `pandas.corr()`, or a user defined function that returns
a value between -1 and 1.

Then, from each group of correlated features, it will try and identify the best candidate
based on the above criteria.

If the criteria is based on feature importance, :class:`SmartCorrelatedSelection` will
train a machine learning model using the correlated feature group, derive the feature importance
from this model, end then keep the feature with the highest important.

:class:`SmartCorrelatedSelection` works with machine learning models that derive coefficients
or feature importance values.

If the criteria is based on variance or cardinality, :class:`SmartCorrelatedSelection` will
determine these attributes for each feature in the group and retain that one with the highest.

If the criteria is based on missing data, :class:`SmartCorrelatedSelection` will determine the
number of NA in each feature from the correlated group and keep the one with less NA.

**Example**

Let's see how to use :class:`SmartCorrelatedSelection` in a toy example. Let's create a
toy dataframe with 4 correlated features:

.. code:: python

    import pandas as pd
    from sklearn.datasets import make_classification
    from feature_engine.selection import SmartCorrelatedSelection

    # make dataframe with some correlated variables
    def make_data():
        X, y = make_classification(n_samples=1000,
                                   n_features=12,
                                   n_redundant=4,
                                   n_clusters_per_class=1,
                                   weights=[0.50],
                                   class_sep=2,
                                   random_state=1)

        # transform arrays into pandas df and series
        colnames = ['var_'+str(i) for i in range(12)]
        X = pd.DataFrame(X, columns=colnames)
        return X

    X = make_data()

Now, we set up :class:`SmartCorrelatedSelection` to find features groups which (absolute)
correlation coefficient is >0.8. From these groups, we want to retain the feature with
highest variance:

.. code:: python

    # set up the selector
    tr = SmartCorrelatedSelection(
        variables=None,
        method="pearson",
        threshold=0.8,
        missing_values="raise",
        selection_method="variance",
        estimator=None,
    )

With `fit()` the transformer finds the correlated variables and selects the one to keep.
With `transform()` it drops them from the dataset:

.. code:: python

    Xt = tr.fit_transform(X)

The correlated feature groups are stored in the transformer's attributes:

.. code:: python

    tr.correlated_feature_sets_

Note that in the second group, 4 features are correlated among themselves.

.. code:: python

    [{'var_0', 'var_8'}, {'var_4', 'var_6', 'var_7', 'var_9'}]

In the following attribute we find the features that will be removed from the dataset:

..  code:: python

    tr.features_to_drop_

.. code:: python

   ['var_0', 'var_4', 'var_6', 'var_9']

If we now go ahead and print the transformed data, we see that the correlated features
have been removed.

.. code:: python

    print(print(Xt.head()))

.. code:: python

          var_1     var_2     var_3     var_5    var_10    var_11     var_8  \
    0 -2.376400 -0.247208  1.210290  0.091527  2.070526 -1.989335  2.070483
    1  1.969326 -0.126894  0.034598 -0.186802  1.184820 -1.309524  2.421477
    2  1.499174  0.334123 -2.233844 -0.313881 -0.066448 -0.852703  2.263546
    3  0.075341  1.627132  0.943132 -0.468041  0.713558  0.484649  2.792500
    4  0.372213  0.338141  0.951526  0.729005  0.398790 -0.186530  2.186741

          var_7
    0 -2.230170
    1 -1.447490
    2 -2.240741
    3 -3.534861
    4 -2.053965


More details
^^^^^^^^^^^^

In this notebook, we show how to use :class:`SmartCorrelatedSelection` with a different
relation metric:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/selection/Smart-Correlation-Selection.ipynb>`_

.. _recursive_elimination:

.. currentmodule:: feature_engine.selection

RecursiveFeatureElimination
============================

:class:`RecursiveFeatureElimination` implements recursive feature elimination. Recursive
feature elimination (RFE) is a backward feature selection process. In Feature-engine's
implementation of RFE, a feature will be kept or removed based on the performance of a
machine learning model without that feature. This differs from Scikit-learn's implementation of
`RFE <https://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.RFE.html>`_
where a feature will be kept or removed based on the feature importance.

This technique begins by building a model on the entire set of variables, then calculates and
stores a model performance metric, and finally computes an importance score for each variable.
Features are ranked by the model’s `coef_` or `feature_importances_` attributes.

In the next step, the least important feature is removed, the model is re-built, and a new performance
metric is determined. If this performance metric is worse than the original one, then,
the feature is kept, (because eliminating the feature clearly caused a drop in model
performance) otherwise, it removed.

The procedure removes now the second to least important feature, trains a new model, determines a
new performance metric, and so on, until it evaluates all the features, from the least
to the most important.

Note that, in Feature-engine's implementation of RFE, the feature importance is used
just to rank features and thus determine the order
in which the features will be eliminated. But whether to retain a feature is determined
based on the decrease in the performance of the model after the feature elimination.

By recursively eliminating features, RFE attempts to eliminate dependencies and
collinearity that may exist in the model.

**Parameters**

Feature-engine's RFE has 2 parameters that need to be determined somewhat arbitrarily by
the user: the first one is the machine learning model which performance will be evaluated. The
second is the threshold in the performance drop that needs to occur, to remove a feature.

RFE is not machine learning model agnostic, this means that the feature selection depends on
the model, and different models may have different subsets of optimal features. Thus, it is
recommended that you use the machine learning model that you finally intend to build.

Regarding the threshold, this parameter needs a bit of hand tuning. Higher thresholds will
of course return fewer features.

**Example**

Let's see how to use this transformer with the diabetes dataset that comes in Scikit-learn.
First, we load the data:


.. code:: python

    import pandas as pd
    from sklearn.datasets import load_diabetes
    from sklearn.linear_model import LinearRegression
    from feature_engine.selection import RecursiveFeatureElimination

    # load dataset
    diabetes_X, diabetes_y = load_diabetes(return_X_y=True)
    X = pd.DataFrame(diabetes_X)
    y = pd.DataFrame(diabetes_y)

Now, we set up :class:`RecursiveFeatureElimination` to select features based on the r2
returned by a Linear Regression model, using 3 fold cross-validation. In this case,
we leave the parameter `threshold` to the default value which is 0.01.

.. code:: python

    # initialize linear regresion estimator
    linear_model = LinearRegression()

    # initialize feature selector
    tr = RecursiveFeatureElimination(estimator=linear_model, scoring="r2", cv=3)

With `fit()` the model finds the most useful features, that is, features that when removed
cause a drop in model performance bigger than 0.01. With `transform()`, the transformer
removes the features from the dataset.

.. code:: python

    # fit transformer
    Xt = tr.fit_transform(X, y)


:class:`RecursiveFeatureElimination` stores the performance of the model trained using all
the features in its attribute:

.. code:: python

    # get the initial linear model performance, using all features
    tr.initial_model_performance_
    
.. code:: python

    0.488702767247119

:class:`RecursiveFeatureElimination`  also stores the change in the performance caused by
removing every feature.

..  code:: python

    # Get the performance drift of each feature
    tr.performance_drifts_
    
..  code:: python

    {0: -0.0032796652347705235,
     9: -0.00028200591588534163,
     6: -0.0006752869546966522,
     7: 0.00013883578730117252,
     1: 0.011956170569096924,
     3: 0.028634492035512438,
     5: 0.012639090879036363,
     2: 0.06630127204137715,
     8: 0.1093736570697495,
     4: 0.024318093565432353}

:class:`RecursiveFeatureElimination` also stores the features that will be dropped based
n the given threshold.

..  code:: python

    # the features to remove
    tr.features_to_drop_

..  code:: python

    [0, 6, 7, 9]

If we now print the transformed data, we see that the features above were removed.

..  code:: python

    print(Xt.head())

..  code:: python

              1         3         5         2         8         4
    0  0.050680  0.021872 -0.034821  0.061696  0.019908 -0.044223
    1 -0.044642 -0.026328 -0.019163 -0.051474 -0.068330 -0.008449
    2  0.050680 -0.005671 -0.034194  0.044451  0.002864 -0.045599
    3 -0.044642 -0.036656  0.024991 -0.011595  0.022692  0.012191
    4 -0.044642  0.021872  0.015596 -0.036385 -0.031991  0.003935
.. _target_mean_selection:

.. currentmodule:: feature_engine.selection


SelectByTargetMeanPerformance
=============================

:class:`SelectByTargetMeanPerformance()` uses the mean value of the target per unique
category, or per interval if the variable is numerical, as proxy for prediction. And with
this prediction, it determines a performance metric against the target.

:class:`SelectByTargetMeanPerformance()` splits the data into two halves. The first half
is used as training set to determine the mappings from category to mean target value or from
interval to mean target value. The second half is the test set, where the categories and
intervals will be mapped to the determined values, these will be considered a prediction,
and assessed against the target to determine a performance metric.

These feature selection idea is very simple; it involves taking the mean of the
responses (target) for each level (category or interval), and so amounts to a least
squares fit on a single categorical variable against a response variable, with the
categories in the continuous variables defined by intervals.

Despite its simplicity, the method has a number of advantages:

- Speed: Computing means and intervals is fast, straightforward and efficient
- Stability with respect to scale: Extreme values for continuous variables do not skew predictions as they would in many models
- Comparability between continuous and categorical variables.
- Accommodation of non-linearities.
- Does not require encoding categorical variables into numbers.

The methods has as well some limitations. First, the selection of the number of intervals
as well as the threshold are arbitrary. And also, rare categories and
very skewed variables will cause over-fitting or the model to raise an error if NAN
are accidentally introduced.

**Example**

Let's see how to use this method to select variables in the Titanic dataset. We choose
this data because it has a mix of numerical and categorical variables.

Let's go ahead and load the data and separate it into train and test:

.. code:: python

    import pandas as pd
    import numpy as np
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import roc_auc_score
    from feature_engine.selection import SelectByTargetMeanPerformance

    # load data
    data = pd.read_csv('../titanic.csv')

    # extract cabin letter
    data['cabin'] = data['cabin'].str[0]

    # replace infrequent cabins by N
    data['cabin'] = np.where(data['cabin'].isin(['T', 'G']), 'N', data['cabin'])

    # cap maximum values
    data['parch'] = np.where(data['parch']>3,3,data['parch'])
    data['sibsp'] = np.where(data['sibsp']>3,3,data['sibsp'])

    # cast variables as object to treat as categorical
    data[['pclass','sibsp','parch']] = data[['pclass','sibsp','parch']].astype('O')

    # separate train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
        data.drop(['survived'], axis=1),
        data['survived'],
        test_size=0.3,
        random_state=0)

Now, we set up :class:`SelectByTargetMeanPerformance()`. We will examine the roc-auc
using 2 fold cross-validation. We will separate numerical variables into
equal frequency intervals. And we will retain those variables where the roc-auc is bigger
than 0.6.

.. code:: python

    # feature engine automates the selection for both categorical and numerical
    # variables
    sel = SelectByTargetMeanPerformance(
        variables=None,
        scoring="roc_auc_score",
        threshold=0.6,
        bins=3,
        strategy="equal_frequency", 
        cv=2,# cross validation
        random_state=1, #seed for reproducibility
    )

With `fit()` the transformer:

- replaces categories by the target mean
- sorts numerical variables into equal frequency bins
- replaces bins by the target mean
- using the target mean encoded variables returns the roc-auc
- selects features which roc-auc >0.6

.. code:: python

    # find important features
    sel.fit(X_train, y_train)

The transformer stores the categorical variables identified in the data:

.. code:: python

    sel.variables_categorical_

.. code:: python

    ['pclass', 'sex', 'sibsp', 'parch', 'cabin', 'embarked']

The transformer also stores the numerical variables:

.. code:: python

    sel.variables_numerical_

.. code:: python

    ['age', 'fare']

:class:`SelectByTargetMeanPerformance()` also stores the roc-auc per feature:

.. code:: python

    sel.feature_performance_

.. code:: python

    {'pclass': 0.6802934787230475,
     'sex': 0.7491365252482871,
     'age': 0.5345141148737766,
     'sibsp': 0.5720480307315783,
     'parch': 0.5243557188989476,
     'fare': 0.6600883312700917,
     'cabin': 0.6379782658154696,
     'embarked': 0.5672382248783936}

And the features that will be dropped from the data:

.. code:: python

    sel.features_to_drop_

.. code:: python

    ['age', 'sibsp', 'parch', 'embarked']

With `transform()` we can go ahead and drop the features:

.. code:: python

    # remove features
    X_train = sel.transform(X_train)
    X_test = sel.transform(X_test)

    X_train.shape, X_test.shape

.. code:: python

    ((914, 4), (392, 4))


More details
^^^^^^^^^^^^

Check also:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/selection/Select-by-Target-Mean-Encoding.ipynb>`_
.. _drop_duplicate:

.. currentmodule:: feature_engine.selection

DropDuplicateFeatures
=====================

The :class:`DropDuplicateFeatures()` finds and removes duplicated variables from a dataframe.
Duplicated features are identical features, regardless of the variable or column name. If
they show the same values for every observation, then they are considered duplicated.

The transformer will automatically evaluate all variables, or alternatively, you can pass a
list with the variables you wish to have examined. And it works with numerical and categorical
features.

**Example**

Let's see how to use :class:`DropDuplicateFeatures()` in an example with the Titanic dataset.
These dataset does not have duplicated features, so we will add a few manually:

.. code:: python

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.model_selection import train_test_split

    from feature_engine.selection import DropDuplicateFeatures

    def load_titanic():
            data = pd.read_csv('https://www.openml.org/data/get_csv/16826755/phpMYEkMl')
            data = data.replace('?', np.nan)
            data['cabin'] = data['cabin'].astype(str).str[0]
            data = data[['pclass', 'survived', 'sex', 'age', 'sibsp', 'parch', 'cabin', 'embarked']]
            data = pd.concat([data, data[['sex', 'age', 'sibsp']]], axis=1)
            data.columns = ['pclass', 'survived', 'sex', 'age', 'sibsp', 'parch', 'cabin', 'embarked',
                            'sex_dup', 'age_dup', 'sibsp_dup']
            return data

    # load data as pandas dataframe
    data = load_titanic()
    data.head()

    # Separate into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
                data.drop(['survived'], axis=1),
                data['survived'], test_size=0.3, random_state=0)

Now, we set up :class:`DropDuplicateFeatures()` to find the duplications:

.. code:: python

    # set up the transformer
    transformer = DropDuplicateFeatures()

With `fit()` the transformer finds the duplicated features, With `transform()` it removes
them:

.. code:: python

    # fit the transformer
    transformer.fit(X_train)

    # transform the data
    train_t = transformer.transform(X_train)

If we examine the variable names of the transformed dataset, we see that the duplicated
features are not present:

.. code:: python

    train_t.columns

.. code:: python

    Index(['pclass', 'sex', 'age', 'sibsp', 'parch', 'cabin', 'embarked'], dtype='object')

The features that are removed are stored in the transformer's attribute:

..  code:: python

    transformer.features_to_drop_

.. code:: python

    {'age_dup', 'sex_dup', 'sibsp_dup'}

And the transformer also stores the groups of duplicated features, which could be useful
if we have groups where more than 2 features are identical.

.. code:: python

    transformer.duplicated_feature_sets_

.. code:: python

    [{'sex', 'sex_dup'}, {'age', 'age_dup'}, {'sibsp', 'sibsp_dup'}]


More details
^^^^^^^^^^^^

In this Kaggle kernel we use :class:`DropDuplicateFeatures()` together with other feature selection algorithms:

- `Kaggle kernel <https://www.kaggle.com/solegalli/feature-selection-with-feature-engine>`_

.. -*- mode: rst -*-
.. _selection_user_guide:

.. currentmodule:: feature_engine.selection


Feature Selection
=================

Feature-engine's feature selection transformers identify features with low predictive
performance and drop them from the dataset. To our knowledge, the feature selection
algorithms supported by Feature-engine are not yet available in other libraries. These
algorithms have been gathered from data science competitions or used in the industry.


Selection Mechanism Overview
----------------------------

Feature-engine's transformers select features based on 2 mechanism. The first mechanism
involves selecting features based on the features intrinsic characteristics like distributions
or their relationship with other features. The second mechanism involves selecting features
based on their impact on the machine learning model performance. In this context, features
are evaluated individually or as part of a feature group by different algorithms.

.. figure::  ../../images/selectionChart.png
   :align:   center

   Selection mechanisms - Overview

For example, in the first pillar, features will be selected based on the diversity of their
values, changes in their distribution or their relation to other features. This way,
features that show the same value in all or almost all the observations will be dropped,
features which distribution changes in time will be dropped, or duplicated or correlated
features will be dropped.

Algorithms that select features based on individual feature performance will select features
by either training a machine learning model using an individual feature, or estimating model
performance with a single feature using a prediction proxy.

Algorithms that select features based on their performance within a group of variables, will
normally train a model with all the features, and then remove or add or shuffle a feature and
re-evaluate the model performance.

These methods are normally geared to improve the overall performance of the final machine learning model
as well as reducing the feature space.


Selectors Characteristics Overview
----------------------------------

Some Feature-engine's selectors work with categorical variables off-the-shelf and/or allow
missing data in the variables. These gives you the opportunity to quickly screen features
before jumping into any feature engineering.

In the following tables we highlight the main Feature-engine selectors characteristics:

Selection based on feature characteristics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

============================================ ======================= ============= ====================================================================================
    Transformer                                Categorical variables   Allows NA	    Description
============================================ ======================= ============= ====================================================================================
:class:`DropFeatures()`                         √	                      √	            Drops arbitrary features determined by user
:class:`DropConstantFeatures()`  	            √	                      √	            Drops constant and quasi-constant features
:class:`DropDuplicateFeatures()`                √	                      √             Drops features that are duplicated
:class:`DropCorrelatedFeatures()`               ×	                      √	            Drops features that are correlated
:class:`SmartCorrelatedSelection()`	            ×	                      √	            From a correlated feature group drops the less useful features
:class:`DropHighPSIFeatures()`	                ×	                      √	            Drops features with high Population Stability Index
============================================ ======================= ============= ====================================================================================

Selection based on model performance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

============================================ ======================= ============= ====================================================================================
    Transformer                                Categorical variables   Allows NA	    Description
============================================ ======================= ============= ====================================================================================
:class:`SelectByShuffling()`	                ×	                      ×	            Selects features if shuffling their values causes a drop in model performance
:class:`SelectBySingleFeaturePerformance()`	    ×	                      ×	            Removes observations with missing data from the dataset
:class:`SelectByTargetMeanPerformance()`        √                         ×             Using the target mean as performance proxy, selects high performing features
:class:`RecursiveFeatureElimination()`          ×                         ×             Removes features recursively by evaluating model performance
:class:`RecursiveFeatureAddition()`             ×                         ×             Adds features recursively by evaluating model performance
============================================ ======================= ============= ====================================================================================

In short, selection procedures that require training a machine learning model from Scikit-learn
require numerical variables without missing data. Selection procedures based on correlation work
only with numerical variables but allow missing data. Methods that determine duplication or
the number of unique values can work with both numerical and categorical variables and support
missing data as well.

The :class:`SelectBySingleFeaturePerformance()` uses the target mean value as proxy for prediction,
replacing categories or variable intervals by these values and then determining a performance metric.
Thus, it is suitable for both categorical and numerical variables. In its current implementation,
it does not support missing data.

:class:`DropHighPSIFeatures()` allows to remove features with changes in their distribution. This is done by
splitting the input dataframe in two parts and comparing the distribution of each feature in the two
parts. The metric used to assess distribution shift is the Population Stability Index (PSI). Removing
unstable features may lead to more robust models. In fields like Credit Risk Modelling, the Regulator
often requires the PSI of the final feature set to be below are given threshold.

Throughout the user guide, you will find more details about each of the feature selection procedures.

Feature Selection Algorithms
----------------------------

Click below to find more details on how to use each one of the transformers.

.. toctree::
   :maxdepth: 1

   DropFeatures
   DropConstantFeatures
   DropDuplicateFeatures
   DropCorrelatedFeatures
   SmartCorrelatedSelection
   DropHighPSIFeatures
   SelectByShuffling
   SelectBySingleFeaturePerformance
   SelectByTargetMeanPerformance
   RecursiveFeatureElimination
   RecursiveFeatureAddition


Additional Resources
--------------------

More details about feature selection can be found in the following resources:

- `Feature Selection Online Course <https://www.udemy.com/course/feature-selection-for-machine-learning/?referralCode=186501DF5D93F48C4F71>`_
- `Feature Selection for Machine Learning: A comprehensive Overview <https://trainindata.medium.com/feature-selection-for-machine-learning-a-comprehensive-overview-bd571db5dd2d>`_

.. _log_cp:

.. currentmodule:: feature_engine.transformation

LogCpTransformer
================

The :class:`LogCpTransformer()` applies the transformation log(x + C), where C is a
positive constant.

You can enter the positive quantity to add to the variable. Alternatively, the transformer
will find the necessary quantity to make all values of the variable positive.

**Example**

Let's load the boston house prices dataset that comes baked into Scikit-learn and
separate it into train and test sets.

.. code:: python

    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn.model_selection import train_test_split
    from sklearn.datasets import load_boston

    from feature_engine import transformation as vt

    # Load dataset
    X, y = load_boston(return_X_y=True)
    X = pd.DataFrame(X)

    # Separate into train and test sets
    X_train, X_test, y_train, y_test =  train_test_split(X, y, test_size=0.3, random_state=0)


Now we want to apply the logarithm to 2 of the variables in the dataset using the
:class:`LogCpTransformer()`. We want the transformer to detect automatically the
quantity "C" that needs to be added to the variable:

.. code:: python

    # set up the variable transformer
    tf = vt.LogCpTransformer(variables = [7, 12], C="auto")

    # fit the transformer
    tf.fit(X_train)

With `fit()` the :class:`LogCpTransformer()` learns the quantity "C" and stores it as
an attribute. We can visualise the learned parameters as follows:

.. code:: python

    # learned constant C
    tf.C_

.. code:: python

    {7: 2.1742, 12: 2.73}

We can now go ahead and transform the variables:

.. code:: python

    # transform the data
    train_t= tf.transform(X_train)
    test_t= tf.transform(X_test)

Then we can plot the original variable distribution:

.. code:: python

    # un-transformed variable
    X_train[12].hist()

.. image:: ../../images/logcpraw.png

And the distribution of the transformed variable:

.. code:: python

    # transformed variable
    train_t[12].hist()

.. image:: ../../images/logcptransform.png

More details
^^^^^^^^^^^^

You can find more details about the :class:`LogCpTransformer()` here:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/transformation/LogCpTransformer.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _yeojohnson:

.. currentmodule:: feature_engine.transformation

YeoJohnsonTransformer
=====================

The :class:`YeoJohnsonTransformer()` applies the Yeo-Johnson transformation to the numerical variables.

The Yeo-Johnson transformation is defined as:

.. image:: ../../images/yeojohnsonformula.png

where Y is the response variable and λ is the transformation parameter.

The Yeo-Johnson transformation implemented by this transformer is that of
`SciPy.stats <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.yeojohnson.html>`_.

**Example**

Let's load the house prices dataset and  separate it into train and test sets (more
details about the dataset :ref:`here <datasets>`).

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine import transformation as vt

	# Load dataset
	data = data = pd.read_csv('houseprice.csv')

	# Separate into train and test sets
	X_train, X_test, y_train, y_test =  train_test_split(
		    data.drop(['Id', 'SalePrice'], axis=1),
		    data['SalePrice'], test_size=0.3, random_state=0)

Now we apply the Yeo-Johnson transformation to the 2 indicated variables:

.. code:: python

	# set up the variable transformer
	tf = vt.YeoJohnsonTransformer(variables = ['LotArea', 'GrLivArea'])

	# fit the transformer
	tf.fit(X_train)

With `fit()`, the :class:`YeoJohnsonTransformer()` learns the optimal lambda for the transformation.
Now we can go ahead and trasnform the data:

.. code:: python

	# transform the data
	train_t= tf.transform(X_train)
	test_t= tf.transform(X_test)

Next, we make a histogram of the original variable distribution:

.. code:: python

	# un-transformed variable
	X_train['LotArea'].hist(bins=50)

.. image:: ../../images/lotarearaw.png

And now, we can explore the distribution of the variable after the transformation:

.. code:: python

	# transformed variable
	train_t['LotArea'].hist(bins=50)


.. image:: ../../images/lotareayeojohnson.png


More details
^^^^^^^^^^^^

You can find more details about the :class:`YeoJohnsonTransformer()` here:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/transformation/YeoJohnsonTransformer.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _power:

.. currentmodule:: feature_engine.transformation

PowerTransformer
================

The :class:`PowerTransformer()` applies power or exponential transformations to numerical
variables.

Let's load the house prices dataset and  separate it into train and test sets (more
details about the dataset :ref:`here <datasets>`).

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine import transformation as vt

	# Load dataset
	data = data = pd.read_csv('houseprice.csv')

	# Separate into train and test sets
	X_train, X_test, y_train, y_test =  train_test_split(
		    data.drop(['Id', 'SalePrice'], axis=1),
		    data['SalePrice'], test_size=0.3, random_state=0)

Now we want to apply the square root to 2 variables in the dataframe:

.. code:: python

	# set up the variable transformer
	tf = vt.PowerTransformer(variables = ['LotArea', 'GrLivArea'], exp=0.5)

	# fit the transformer
	tf.fit(X_train)

The transformer does not learn any parameters. So we can go ahead and transform the
variables:

.. code:: python

	# transform the data
	train_t= tf.transform(X_train)
	test_t= tf.transform(X_test)

Finally, we can plot the original variable distribution:

.. code:: python

	# un-transformed variable
	X_train['LotArea'].hist(bins=50)

.. image:: ../../images/lotarearaw.png

And now the distribution after the transformation:

.. code:: python

	# transformed variable
	train_t['LotArea'].hist(bins=50)


.. image:: ../../images/lotareapower.png

More details
^^^^^^^^^^^^

You can find more details about the :class:`PowerTransformer()` here:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/transformation/PowerTransformer.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _log_transformer:

.. currentmodule:: feature_engine.transformation

LogTransformer
==============

The :class:`LogTransformer()` will apply the logarithm to the indicated variables. Note
that the logarithm can only be applied to positive values. Thus, if the variable contains
0 or negative variables, this transformer will return and error.

**Example**

Let's load the house prices dataset and  separate it into train and test sets (more
details about the dataset :ref:`here <datasets>`).

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine import transformation as vt

	# Load dataset
	data = pd.read_csv('houseprice.csv')

	# Separate into train and test sets
	X_train, X_test, y_train, y_test =  train_test_split(
		    data.drop(['Id', 'SalePrice'], axis=1),
		    data['SalePrice'], test_size=0.3, random_state=0)

Now we want to apply the logarithm to 2 of the variables in the dataset using the
:class:`LogTransformer()`.

.. code:: python

	# set up the variable transformer
	tf = vt.LogTransformer(variables = ['LotArea', 'GrLivArea'])

	# fit the transformer
	tf.fit(X_train)

With `fit()`, this transformer does not learn any parameters. We can go ahead not an
transform the variables.

.. code:: python

	# transform the data
	train_t= tf.transform(X_train)
	test_t= tf.transform(X_test)

Next, we make a histogram of the original variable distribution:

.. code:: python

	# un-transformed variable
	X_train['LotArea'].hist(bins=50)

.. image:: ../../images/lotarearaw.png

And now, we can explore the distribution of the variable after the logarithm transformation:

.. code:: python

	# transformed variable
	train_t['LotArea'].hist(bins=50)


.. image:: ../../images/lotarealog.png

Note that the transformed variable has a more Gaussian looking distribution.

More details
^^^^^^^^^^^^

You can find more details about the :class:`LogTransformer()` here:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/transformation/LogTransformer.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _reciprocal:

.. currentmodule:: feature_engine.transformation

ReciprocalTransformer
=====================

The :class:`ReciprocalTransformer()` applies the reciprocal transformation 1 / x to
numerical variables.

The :class:`ReciprocalTransformer()` only works with numerical variables with non-zero
values. If a variable contains the value 0, the transformer will raise an error.

Let's load the house prices dataset and  separate it into train and test sets (more
details about the dataset :ref:`here <datasets>`).

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine import transformation as vt

	# Load dataset
	data = data = pd.read_csv('houseprice.csv')

	# Separate into train and test sets
	X_train, X_test, y_train, y_test =  train_test_split(
		    data.drop(['Id', 'SalePrice'], axis=1),
		    data['SalePrice'], test_size=0.3, random_state=0)

Now we want to apply the reciprocal transformation to 2 variables in the dataframe:

.. code:: python

	# set up the variable transformer
	tf = vt.ReciprocalTransformer(variables = ['LotArea', 'GrLivArea'])

	# fit the transformer
	tf.fit(X_train)

The transformer does not learn any parameters. So we can go ahead and transform the
variables:

.. code:: python

	# transform the data
	train_t= tf.transform(X_train)
	test_t= tf.transform(X_test)

Finally, we can plot the original variable distribution:

.. code:: python

	# un-transformed variable
	X_train['LotArea'].hist(bins=50)

.. image:: ../../images/lotarearaw.png

And now the distribution after the transformation:

.. code:: python

	# transformed variable
	train_t['LotArea'].hist(bins=50)


.. image:: ../../images/lotareareciprocal.png

More details
^^^^^^^^^^^^

You can find more details about the :class:`ReciprocalTransformer()` here:


- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/transformation/ReciprocalTransformer.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. _box_cox:

.. currentmodule:: feature_engine.transformation

BoxCoxTransformer
=================

The :class:`BoxCoxTransformer()` applies the BoxCox transformation to numerical variables.

The Box-Cox transform is given by:

.. code:: python

   y = (x**lmbda - 1) / lmbda,  for lmbda != 0
   log(x),                      for lmbda = 0

The BoxCox transformation implemented by this transformer is that of
`SciPy.stats <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.boxcox.html>`_.

The BoxCox transformation works only for strictly positive variables (>=0). If the variable
contains 0 or negative values, the :class:`BoxCoxTransformer()` will return an error.

If the variable contains values <=0, you should try using the :class:`YeoJohnsonTransformer()`
instead.

**Example**

Let's load the house prices dataset and  separate it into train and test sets (more
details about the dataset :ref:`here <datasets>`).

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine import transformation as vt

	# Load dataset
	data = data = pd.read_csv('houseprice.csv')

	# Separate into train and test sets
	X_train, X_test, y_train, y_test =  train_test_split(
		    data.drop(['Id', 'SalePrice'], axis=1),
		    data['SalePrice'], test_size=0.3, random_state=0)


Now we apply the BoxCox transformation to the 2 indicated variables:

.. code:: python

	# set up the variable transformer
	tf = vt.BoxCoxTransformer(variables = ['LotArea', 'GrLivArea'])

	# fit the transformer
	tf.fit(X_train)

With `fit()`, the :class:`BoxCoxTransformer()` learns the optimal lambda for the transformation.
Now we can go ahead and trasnform the data:

.. code:: python

	# transform the data
	train_t= tf.transform(X_train)
	test_t= tf.transform(X_test)

Next, we make a histogram of the original variable distribution:

.. code:: python

	# un-transformed variable
	X_train['LotArea'].hist(bins=50)

.. image:: ../../images/lotarearaw.png

And now, we can explore the distribution of the variable after the transformation:

.. code:: python

	# transformed variable
	train_t['GrLivArea'].hist(bins=50)


.. image:: ../../images/lotareaboxcox.png


More details
^^^^^^^^^^^^

You can find more details about the :class:`BoxCoxTransformer()` here:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/transformation/BoxCoxTransformer.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. -*- mode: rst -*-

Variable Transformation
=======================

Feature-engine's variable transformers transform numerical variables with various
mathematical transformations.

Variable transformations are commonly used to spread the values of the original variables
over a wider value range. See the following illustration:

.. figure::  ../../images/Variable_Transformation.png
   :align:   center

**Note**

Note however, that improving the value spread is not always possible and it depends
on the nature of the variable.

**Transformers**

.. toctree::
   :maxdepth: 1

   LogTransformer
   LogCpTransformer
   ReciprocalTransformer
   PowerTransformer
   BoxCoxTransformer
   YeoJohnsonTransformer.. _match_variables:

.. currentmodule:: feature_engine.preprocessing

MatchVariables
==============

:class:`MatchVariables()` ensures that the columns in the test set are identical to those
in the train set.

If the test set contains additional columns, they are dropped. Alternatively, if the
test set lacks columns that were present in the train set, they will be added with a
value determined by the user, for example np.nan. :class:`MatchVariables()` will also
return the variables in the order seen in the train set.

Let's explore this with an example. First we load the Titanic dataset and split it into
a train and a test set:

.. code:: python

    import numpy as np
    import pandas as pd

    from feature_engine.preprocessing import MatchVariables


    # Load dataset
    def load_titanic():
        data = pd.read_csv('https://www.openml.org/data/get_csv/16826755/phpMYEkMl')
        data = data.replace('?', np.nan)
        data['cabin'] = data['cabin'].astype(str).str[0]
        data['pclass'] = data['pclass'].astype('O')
        data['age'] = data['age'].astype('float')
        data['fare'] = data['fare'].astype('float')
        data['embarked'].fillna('C', inplace=True)
        data.drop(
            labels=['name', 'ticket', 'boat', 'body', 'home.dest'],
            axis=1, inplace=True,
        )
        return data

    # load data as pandas dataframe
    data = load_titanic()

    # Split test and train
    train = data.iloc[0:1000, :]
    test = data.iloc[1000:, :]

Now, we set up :class:`MatchVariables()` and fit it to the train set.

.. code:: python

    # set up the transformer
    match_cols = MatchVariables(missing_values="ignore")

    # learn the variables in the train set
    match_cols.fit(train)

:class:`MatchVariables()` stores the variables from the train set in its attribute:

.. code:: python

    # the transformer stores the input variables
    match_cols.input_features_

.. code:: python

    ['pclass',
     'survived',
     'sex',
     'age',
     'sibsp',
     'parch',
     'fare',
     'cabin',
     'embarked']

Now, we drop some columns in the test set.

.. code:: python

    # Let's drop some columns in the test set for the demo
    test_t = test.drop(["sex", "age"], axis=1)

    test_t.head()

.. code:: python

         pclass  survived  sibsp  parch     fare cabin embarked
    1000      3         1      0      0   7.7500     n        Q
    1001      3         1      2      0  23.2500     n        Q
    1002      3         1      2      0  23.2500     n        Q
    1003      3         1      2      0  23.2500     n        Q
    1004      3         1      0      0   7.7875     n        Q

If we transform the dataframe with the dropped columns using :class:`MatchVariables()`,
we see that the new dataframe contains all the variables, and those that were missing
are now back in the data, with np.nan values as default.

.. code:: python

    # the transformer adds the columns back
    test_tt = match_cols.transform(test_t)

    test_tt.head()

.. code:: python

    The following variables are added to the DataFrame: ['sex', 'age']

         pclass  survived  sex  age  sibsp  parch     fare cabin embarked
    1000      3         1  NaN  NaN      0      0   7.7500     n        Q
    1001      3         1  NaN  NaN      2      0  23.2500     n        Q
    1002      3         1  NaN  NaN      2      0  23.2500     n        Q
    1003      3         1  NaN  NaN      2      0  23.2500     n        Q
    1004      3         1  NaN  NaN      0      0   7.7875     n        Q



Note how the missing columns were added back to the transformed test set, with
missing values, in the position (i.e., order) in which they were in the train set.

Similarly, if the test set contained additional columns, those would be removed. To
test that, let's add some extra columns to the test set:

.. code:: python

    # let's add some columns for the demo
    test_t[['var_a', 'var_b']] = 0

    test_t.head()

.. code:: python

         pclass  survived  sibsp  parch     fare cabin embarked  var_a  var_b
    1000      3         1      0      0   7.7500     n        Q      0      0
    1001      3         1      2      0  23.2500     n        Q      0      0
    1002      3         1      2      0  23.2500     n        Q      0      0
    1003      3         1      2      0  23.2500     n        Q      0      0
    1004      3         1      0      0   7.7875     n        Q      0      0


And now, we transform the data with :class:`MatchVariables()`:

.. code:: python

    test_tt = match_cols.transform(test_t)

    test_tt.head()

.. code:: python

    The following variables are added to the DataFrame: ['age', 'sex']
    The following variables are dropped from the DataFrame: ['var_a', 'var_b']

         pclass  survived  sex  age  sibsp  parch     fare cabin embarked
    1000      3         1  NaN  NaN      0      0   7.7500     n        Q
    1001      3         1  NaN  NaN      2      0  23.2500     n        Q
    1002      3         1  NaN  NaN      2      0  23.2500     n        Q
    1003      3         1  NaN  NaN      2      0  23.2500     n        Q
    1004      3         1  NaN  NaN      0      0   7.7875     n        Q


Now, the transformer simultaneously added the missing columns with NA as values and
removed the additional columns from the resulting dataset.

By default, :class:`MatchVariables()` will print out messages indicating which variables
were added or removed. We can switch off the messages through the parameter `verbose`.

When to use the transformer
^^^^^^^^^^^^^^^^^^^^^^^^^^^

These transformer is useful in "predict then optimize type of problems". In such cases,
a machine learning model is trained on a certain dataset, with certain input features.
Then, test sets are "post-processed" according to scenarios that want to be modelled.
For example, "what would have happened if the customer received an email campaign"?
where the variable "receive_campaign" would be turned from 0 -> 1.

While creating these modelling datasets, a lot of meta data e.g., "scenario number",
"time scenario was generated", etc, could be added to the data. Then we need to pass
these data over to the model to obtain the modelled prediction.

:class:`MatchVariables()` provides an easy an elegant way to remove the additional metadeta,
while returning datasets with the input features in the correct order, allowing the
different scenarios to be modelled directly inside a machine learning pipeline.

More details
^^^^^^^^^^^^

You can also find a similar implementation of the example shown in this page in the
following Jupyter notebook:

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/preprocessing/MatchVariables.ipynb>`_

All notebooks can be found in a `dedicated repository <https://github.com/feature-engine/feature-engine-examples>`_.
.. -*- mode: rst -*-

Preprocessing
=============

Feature-engine's preprocessing transformers apply general data pre-processing
and transformation procedures.

.. toctree::
   :maxdepth: 1

   MatchVariables
.. _datetime_features:

.. currentmodule:: feature_engine.datetime

DatetimeFeatures
================

:class:`DatetimeFeatures()` extracts several datetime features from datetime
variables. It works with variables whose original dtype is datetime, and also with
object-like and categorical variables, provided that they can be parsed into datetime
format. It *cannot* extract features from numerical variables.

Oftentimes, datasets contain information related to dates and/or times at which an event
occurred. In pandas dataframes, these datetime variables can be cast as datetime or,
more generically, as object.

Datetime variables, in their raw format, are generally not suitable to train machine
learning models. Yet, an enormous amount of information can be extracted from them.

:class:`DatetimeFeatures()` is able to extract many numerical and binary date and
time features from these datetime variables. Among these features we can find the month
in which an event occurred, the day of the week, or whether that day was a weekend day.

With :class:`DatetimeFeatures()` you can choose which date and time features
to extract from your datetime variables. You can also extract date and time features
from a subset of datetime variables in your data.

Examples
--------
Through the following examples we highlight the functionality and versatility of
:class:`DatetimeFeatures()`.

Extract date features
~~~~~~~~~~~~~~~~~~~~~

In this example, we are going to extract three date features from a
specific variable in the dataframe. In particular, we are interested
in the month, the day of the year, and whether that day was the last
day of its correspondent month.

First, we will create a toy dataframe with 2 date variables.

.. code:: python

    import pandas as pd
    from feature_engine.datetime import DatetimeFeatures

    toy_df = pd.DataFrame({
        "var_date1": ['May-1989', 'Dec-2020', 'Jan-1999', 'Feb-2002'],
        "var_date2": ['06/21/12', '02/10/98', '08/03/10', '10/31/20'],
    })

Now, we will extract the variables month, month-end and the day of the year from the
second datetime variable in our dataset.

.. code:: python

    dtfs = DatetimeFeatures(
        variables="var_date2",
        features_to_extract=["month", "month_end", "day_of_year"]
    )

    df_transf = dtfs.fit_transform(toy_df)

    print(df_transf)

.. code:: python

      var_date1  var_date2_month  var_date2_month_end  var_date2_day_of_year
    0  May-1989                6                    0                    173
    1  Dec-2020                2                    0                     41
    2  Jan-1999                8                    0                    215
    3  Feb-2002               10                    1                    305


With `transform()`, the features extracted from the datetime variable are added to the
dataframe.

By default, :class:`DatetimeFeatures()` drops the variable from which the date and time
features were extracted, in this case, *var_date2*. To keep the variable, we just need
to indicate `drop_original=False` when initializing the transformer.

Extract time features
~~~~~~~~~~~~~~~~~~~~~

In this example, we are going to extract the feature *minute* from the two time
variables in our dataset.

First, let's create a toy dataset with 2 time variables and an object variable.

.. code:: python 

    import pandas as pd
    from feature_engine.datetime import DatetimeFeatures

    toy_df = pd.DataFrame({
        "not_a_dt": ['not', 'a', 'date', 'time'],
        "var_time1": ['12:34:45', '23:01:02', '11:59:21', '08:44:23'],
        "var_time2": ['02:27:26', '10:10:55', '17:30:00', '18:11:18'],
    })

:class:`DatetimeFeatures()` automatically finds all variables that can be parsed to
datetime. So if we want to extract time features from all our datetime variables, we
don't need to specify them.

.. code:: python

    dfts = DatetimeFeatures(features_to_extract=["minute"])

    df_transf = dfts.fit_transform(toy_df)

    print(df_transf)

.. code:: python

      not_a_dt  var_time1_minute  var_time2_minute
    0      not                34                27
    1        a                 1                10
    2     date                59                30
    3     time                44                11


The transformer found two variales in the dataframe that can be cast to datetime and
proceeded to extract the requested feature from them.

We can find the variables detected as datetime by the transformer in one of its
attributes.

.. code:: python

    dfts.variables_

.. code:: python

    ['var_time1', 'var_time2']

Again, the original datetime variables are dropped from the data by default. If we
want to keep them, we just need to indicate `drop_original=False` when initializing
the transformer.

Extract date and time features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, we will combine what we have seen in the previous two examples
and extract a date feature - *year* - and time feature - *hour* -
from two variables that contain both date and time information.

Let's go ahead and create a toy dataset.

.. code:: python

    import pandas as pd
    from feature_engine.datetime import DatetimeFeatures

    toy_df = pd.DataFrame({
        "var_dt1": pd.date_range("2018-01-01", periods=3, freq="H"),
        "var_dt2": ['08/31/00 12:34:45', '12/01/90 23:01:02', '04/25/01 11:59:21'],
        "var_dt3": ['03/02/15 02:27:26', '02/28/97 10:10:55', '11/11/03 17:30:00'],
    })

Now, we set up the :class:`DatetimeFeatures()` to extract features only from 2 of the
variables. In this case, we do not want to drop the datetime variable after extracting
the features.

.. code:: python

    dfts = DatetimeFeatures(
        variables=["var_dt1", "var_dt3"],
        features_to_extract=["year", "hour"],
        drop_original=False,
    )
    df_transf = dfts.fit_transform(toy_df)

    print(df_transf)

.. code:: python

                  var_dt1            var_dt2            var_dt3  var_dt1_year  \
    0 2018-01-01 00:00:00  08/31/00 12:34:45  03/02/15 02:27:26          2018
    1 2018-01-01 01:00:00  12/01/90 23:01:02  02/28/97 10:10:55          2018
    2 2018-01-01 02:00:00  04/25/01 11:59:21  11/11/03 17:30:00          2018

       var_dt1_hour  var_dt3_year  var_dt3_hour
    0             0          2015             2
    1             1          1997            10
    2             2          2003            17

And that is it. The additional features are now added in the dataframe.

Important
~~~~~~~~~

We highly recommend specifying the date and time features that you would like to extract
from your datetime variables. If you have too many time variables, this might not be
possible. In this case, keep in mind that if you extract date features from variables
that have only time, or time features from variables that have only dates, your features
will be meaningless.

Let's explore the outcome with an example. We create a dataset with only time variables.

.. code:: python

    import pandas as pd
    from feature_engine.datetime import DatetimeFeatures

    toy_df = pd.DataFrame({
        "not_a_dt": ['not', 'a', 'date', 'time'],
        "var_time1": ['12:34:45', '23:01:02', '11:59:21', '08:44:23'],
        "var_time2": ['02:27:26', '10:10:55', '17:30:00', '18:11:18'],
    })

And now we mistakenly extract only date features.

.. code:: python

    dfts = DatetimeFeatures(
        features_to_extract=["year", "month", "day_of_week"],
    )
    df_transf = dfts.fit_transform(toy_df)

    print(df_transf)

.. code:: python

      not_a_dt  var_time1_year  var_time1_month  var_time1_day_of_week  var_time2_year \
    0      not            2021               12                      2            2021
    1        a            2021               12                      2            2021
    2     date            2021               12                      2            2021
    3     time            2021               12                      2            2021

       var_time2_month  var_time2_day_of_week
    0               12                      2
    1               12                      2
    2               12                      2
    3               12                      2

The transformer will still create features derived from today's date (the date of
creating the docs).

If instead we have a dataframe with only date variables.

.. code:: python

    import pandas as pd
    from feature_engine.datetime import DatetimeFeatures

    toy_df = pd.DataFrame({
        "var_date1": ['May-1989', 'Dec-2020', 'Jan-1999', 'Feb-2002'],
        "var_date2": ['06/21/12', '02/10/98', '08/03/10', '10/31/20'],
    })

And mistakenly extract the hour and the minute.

.. code:: python

    dfts = DatetimeFeatures(
        features_to_extract=["hour", "minute"],
    )
    df_transf = dfts.fit_transform(toy_df)

    print(df_transf)

.. code:: python

       var_date1_hour  var_date1_minute  var_date2_hour  var_date2_minute
    0               0                 0               0                 0
    1               0                 0               0                 0
    2               0                 0               0                 0
    3               0                 0               0                 0

The new features will contain the value 0.

Automating feature extraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can indicate which features we want to extract from the datetime variables as we did
in the previous examples, passing the feature names in lists. Alternatively,
:class:`DatetimeFeatures()` has default options to extract a group of commonly used
features, or all supported features.

Let's first create a toy dataframe.

.. code:: python

    import pandas as pd
    from feature_engine.datetime import DatetimeFeatures

    toy_df = pd.DataFrame({
        "var_dt1": pd.date_range("2018-01-01", periods=3, freq="H"),
        "var_dt2": ['08/31/00 12:34:45', '12/01/90 23:01:02', '04/25/01 11:59:21'],
        "var_dt3": ['03/02/15 02:27:26', '02/28/97 10:10:55', '11/11/03 17:30:00'],
    })

Now, we will extract the most common date and time features from one of the variables.
To do this, we leave the parameter `features_to_extract` to `None`.

.. code:: python

    dfts = DatetimeFeatures(
        variables=["var_dt1"],
        features_to_extract=None,
        drop_original=False,
    )

    df_transf = dfts.fit_transform(toy_df)

    print(df_transf)

.. code:: python

                  var_dt1            var_dt2            var_dt3  var_dt1_month  \
    0 2018-01-01 00:00:00  08/31/00 12:34:45  03/02/15 02:27:26              1
    1 2018-01-01 01:00:00  12/01/90 23:01:02  02/28/97 10:10:55              1
    2 2018-01-01 02:00:00  04/25/01 11:59:21  11/11/03 17:30:00              1

       var_dt1_year  var_dt1_day_of_week  var_dt1_day_of_month  var_dt1_hour  \
    0          2018                    0                     1             0
    1          2018                    0                                   1
    2          2018                    0                  1                2

        var_dt1_minute    var_dt1_second
    0               0                  0
    1               0                  0
    2               0                  0

Our new dataset contains the original features plus the new variables extracted
from them.

We can find the group of features extracted by the transformer in its attribute.

.. code:: python

    dfts.features_to_extract_

.. code:: python

    ['month',
     'year',
     'day_of_week',
     'day_of_month',
     'hour',
     'minute',
     'second']

We can also extract all supported features automatically.

.. code:: python

    dfts = DatetimeFeatures(
        variables=["var_dt1"],
        features_to_extract='all',
        drop_original=False,
    )

    df_transf = dfts.fit_transform(toy_df)

    print(df_transf)

.. code:: python

                  var_dt1            var_dt2            var_dt3  var_dt1_month  \
    0 2018-01-01 00:00:00  08/31/00 12:34:45  03/02/15 02:27:26              1
    1 2018-01-01 01:00:00  12/01/90 23:01:02  02/28/97 10:10:55              1
    2 2018-01-01 02:00:00  04/25/01 11:59:21  11/11/03 17:30:00              1

       var_dt1_quarter  var_dt1_semester  var_dt1_year  \
    0                1                 1          2018
    1                1                 1          2018
    2                1                 1          2018

       var_dt1_week  var_dt1_day_of_week  ...  var_dt1_month_end  var_dt1_quarter_start  \
    0             1                    0  ...                  0                      1
    1             1                    0  ...                  0                      1
    2             1                    0  ...                  0                      1

       var_dt1_quarter_end  var_dt1_year_start  var_dt1_year_end  \
    0                    0                   1                 0
    1                    0                   1                 0
    2                    0                   1                 0

       var_dt1_leap_year  var_dt1_days_in_month  var_dt1_hour  var_dt1_minute  \
    0                  0                     31             0               0
    1                  0                     31             1               0
    2                  0                     31             2               0

       var_dt1_second
    0               0
    1               0
    2               0

We can find the group of features extracted by the transformer in its attribute.

.. code:: python

    dfts.features_to_extract_

.. code:: python

    ['month',
     'quarter',
     'semester',
     'year',
     'week',
     'day_of_week',
     'day_of_month',
     'day_of_year',
     'weekend',
     'month_start',
     'month_end',
     'quarter_start',
     'quarter_end',
     'year_start',
     'year_end',
     'leap_year',
     'days_in_month',
     'hour',
     'minute',
     'second']

Extract and select features automatically
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If we have a dataframe with a date variables, time variables and date and time variables,
we can extract all features, or the most common features from all the variables, and then
go ahead and remove the irrelevant features with the `DropConstantFeatures()` class.

Let's create a dataframe with a mix of datetime variables.

.. code:: python

    import pandas as pd
    from sklearn.pipeline import Pipeline
    from feature_engine.datetime import DatetimeFeatures
    from feature_engine.selection import DropConstantFeatures

    toy_df = pd.DataFrame({
        "var_date": ['06/21/12', '02/10/98', '08/03/10', '10/31/20'],
        "var_time1": ['12:34:45', '23:01:02', '11:59:21', '08:44:23'],
        "var_dt": ['08/31/00 12:34:45', '12/01/90 23:01:02', '04/25/01 11:59:21', '04/25/01 11:59:21'],
    })

Now, we align in a Scikit-learn pipeline the :class:`DatetimeFeatures` and the
`DropConstantFeatures()`. The :class:`DatetimeFeatures` will create date features
derived from today for the time variable, and time features with the value 0 for the
date only variable. `DropConstantFeatures()` will identify and remove these features
from the dataset.

.. code:: python

    pipe = Pipeline([
        ('datetime', DatetimeFeatures()),
        ('drop_constant', DropConstantFeatures()),
    ])

    pipe.fit(toy_df)

.. code:: python

    Pipeline(steps=[('datetime', DatetimeFeatures()),
                    ('drop_constant', DropConstantFeatures())])

.. code:: python

    df_transf = pipe.transform(toy_df)

    print(df_transf)

.. code:: python

       var_date_month  var_date_year  var_date_day_of_week  var_date_day_of_month  \
    0               6           2012                     3                     21
    1               2           1998                     1                     10
    2               8           2010                     1                      3
    3              10           2020                     5                     31

       var_time1_hour  var_time1_minute  var_time1_second  var_dt_month  \
    0              12                34                45             8
    1              23                 1                 2            12
    2              11                59                21             4
    3               8                44                23             4

       var_dt_year  var_dt_day_of_week  var_dt_day_of_month  var_dt_hour  \
    0         2000                   3                   31           12
    1         1990                   5                    1           23
    2         2001                   2                   25           11
    3         2001                   2                   25           11

       var_dt_minute   var_dt_second
    0             34              45
    1              1               2
    2             59              21
    3             59              21

As you can see, we do not have the constant features in the transformed dataset.

Extract features from time-aware variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Time-aware datetime variables can be particularly cumbersome to work with as far
as the format goes. We will briefly show how :class:`DatetimeFeatures()` deals
with such variables in three different scenarios.

**Case 1**: our dataset contains a time-aware variable in object format,
with potentially different timezones across different observations. 
We pass `utc=True` when initializing the transformer to make sure it
converts all data to UTC timezone.

.. code:: python

    import pandas as pd
    from feature_engine.datetime import DatetimeFeatures

    toy_df = pd.DataFrame({
        "var_tz": ['12:34:45+3', '23:01:02-6', '11:59:21-8', '08:44:23Z']
    })

    dfts = DatetimeFeatures(
        features_to_extract=["hour", "minute"],
        drop_original=False,
        utc=True
    )

    df_transf = dfts.fit_transform(toy_df)

    print(df_transf)

.. code:: python

           var_tz  var_tz_hour  var_tz_minute
    0  12:34:45+3            9             34
    1  23:01:02-6            5              1
    2  11:59:21-8           19             59
    3   08:44:23Z            8             44


**Case 2**: our dataset contains a variable that is cast as a localized
datetime in a particular timezone. However, we decide that we want to get all
the datetime information extracted as if it were in UTC timezone.

.. code:: python

    import pandas as pd
    from feature_engine.datetime import DatetimeFeatures

    var_tz = pd.Series(['08/31/00 12:34:45', '12/01/90 23:01:02', '04/25/01 11:59:21'])
    var_tz = pd.to_datetime(var_tz)
    var_tz = var_tz.dt.tz_localize("US/eastern")
    var_tz

.. code:: python

    0   2000-08-31 12:34:45-04:00
    1   1990-12-01 23:01:02-05:00
    2   2001-04-25 11:59:21-04:00
    dtype: datetime64[ns, US/Eastern]

We need to pass `utc=True` when initializing the transformer to revert back to the UTC
timezone.

.. code:: python

    toy_df = pd.DataFrame({"var_tz": var_tz})

    dfts = DatetimeFeatures(
        features_to_extract=["day_of_month", "hour"],
        drop_original=False,
        utc=True,
    )

    df_transf = dfts.fit_transform(toy_df)

    print(df_transf)

.. code:: python

                         var_tz  var_tz_day_of_month  var_tz_hour
    0 2000-08-31 12:34:45-04:00                   31           16
    1 1990-12-01 23:01:02-05:00                    2            4
    2 2001-04-25 11:59:21-04:00                   25           15


**Case 3**: given a variable like *var_tz* in the example above, we now want
to extract the features keeping the original timezone localization,
therefore we pass `utc=False` or `None`. In this case, we leave it to `None` which
is the default option.

.. code:: python

    dfts = DatetimeFeatures(
        features_to_extract=["day_of_month", "hour"],
        drop_original=False,
        utc=None,
    )

    df_transf = dfts.fit_transform(toy_df)

    print(df_transf)

.. code:: python

                         var_tz  var_tz_day_of_month  var_tz_hour
    0 2000-08-31 12:34:45-04:00                   31           12
    1 1990-12-01 23:01:02-05:00                    1           23
    2 2001-04-25 11:59:21-04:00                   25           11

Note that the hour extracted from the variable differ in this dataframe respect to the
one obtained in **Case 2**.

More details
^^^^^^^^^^^^

You can find additional examples with a real dataset on how to use
:class:`DatetimeFeatures()` in the following Jupyter notebook.

- `Jupyter notebook <https://nbviewer.org/github/feature-engine/feature-engine-examples/blob/main/datetime/DatetimeFeatures.ipynb>`_
.. -*- mode: rst -*-

Datetime Features
=================

Feature-engine’s datetime transformers are able to extract a wide variety of datetime
features from existing datetime or object-like data.

.. toctree::
   :maxdepth: 1

   DatetimeFeaturesVersion 1.1.X
=============

Version 1.1.2
-------------

Deployed: 31th August 2021

Contributors
~~~~~~~~~~~~

    - Soledad Galli

This small release fixes a Bug in how the OneHotEncoder handles binary categorical variables
when the parameter `drop_last_binary` is set to True. It also ensures that the values in the
`OneHotEncoder.encoder_dict_` are lists of categories and not arrays. These bugs were
introduced in v1.1.0.

Bug fix
~~~~~~~
    - **OneHotEncoder**: drop_last_binary now outputs 1 dummy variable per binary variable when set to true

Version 1.1.1
-------------

Deployed: 6th August 2021

Contributors
~~~~~~~~~~~~

    - Miguel Trema Marrufo
    - Nicolas Galli
    - Soledad Galli

In this release, we add a new transformer, expand the functionality of 2 other
transformers and migrate the repo to its own organisation!

Mayor changes
~~~~~~~~~~~~~
    - Feature-engine is now hosted in its `own Github organisation <https://github.com/feature-engine/feature_engine>`_

New transformer
~~~~~~~~~~~~~~~
    - **LogCpTransformer**: applies the logarithm transformation after adding a constant (**Miguel Trema Marrufo**)

Minor changes
~~~~~~~~~~~~~
    - Expands functionality of `DropCorrelatedFeatures` and `SmartCorrelationSelectionFeature` to accept callables as a correlation function (**Miguel Trema Marrufo**)
    - Adds `inverse_transform` to all transformers from the transformation module (**Nicolas Galli**).

Documentation
~~~~~~~~~~~~~
    - Migrates main repo to `Feature-engine's Github organisation <https://github.com/feature-engine/feature_engine>`_
    - Migrates example jupyter notebooks to `separate repo <https://github.com/feature-engine/feature-engine-examples>`_
    - Adds Roadmap


Version 1.1.0
-------------

Deployed: 22st June 2021

Contributors
~~~~~~~~~~~~
    - Hector Patino
    - Andrew Tan
    - Shubhmay Potdar
    - Agustin Firpo
    - Indy Navarro Vidal
    - Ashok Kumar
    - Chris Samiullah
    - Soledad Galli

In this release, we enforce compatibility with Scikit-learn by adding the
`check_estimator <https://scikit-learn.org/stable/developers/develop.html>`_ tests to
**all transformers** in the package.

In order to pass the tests, we needed to modify some of the internal functionality of
Feature-engine transformers and create new attributes. We tried not to break backwards
compatibility as much as possible.

Mayor changes
~~~~~~~~~~~~~
    - Most transformers have now the additional attribute `variables_` containing the variables that will be modified. The former attribute `variables` is retained. `variables_` will almost always be identical to `variables` except when the transformer is initialised with `variables=None`.
    - The parameter `transformer` in the SklearnTransformerWrapper and the parameter `estimator` in the SelectBySingleFeaturePerformance, SelectByShuffling, RecursiveFeatureElimination and RecursiveFeatureAddition need a compulsory entry, and cannot be left blank when initialising the transformers.
    - Categorical encoders support now variables cast as `category` as well as `object` (**Shubhmay Potdar and Soledad Galli**)
    - Categorical encoders have now the parameter `ignore_format` to allow the transformer to work with any variable type, and not just object or categorical.
    - `CategoricalImputer` has now the parameter `ignore_format` to allow the transformer to work with any variable type, and not just object or categorical.
    - All transformers have now the new attribute `n_features_in` with captures the number of features in the dataset used to train the transformer (during fit()).

Minor changes
~~~~~~~~~~~~~
    - Feature selection transformers support now all cross-validation schemes in the `cv` parameter, and not just an integer. That is, you can initialize the transformer with LOOCV, or StratifiedCV for example.
    - The OneHotEncoder includes additional functionality to return just 1 dummy variable for categorical variables that contain only 2 categories. In the new attribute `variables_binary_` you can identify the original binary variables.
    - MathematicalCombinator now supports use of dataframes with null values (**Agustin Firpo**).

New transformer
~~~~~~~~~~~~~~~
    - **CyclicalTransformer**: applies a cyclical transformation to numerical variables (**Hector Patino**)

Code improvement
~~~~~~~~~~~~~~~~
    - Tests from check_estimator added to all transformers
    - Test for compatibility with Python 3.9 added to circleCI (**Chris Samiullah and Soledad Galli**)
    - Automatic black8 and linting added to tox
    - Additional code fixes (**Andrew Tan and Indy Navarro Vidal**).

Documentation
~~~~~~~~~~~~~
    - Additional comparison tables for imputers and encoders.
    - Updates Readme with new badges and resources.
    - Expanded SklearnWrapper demos in Jupyter notebooks.
    - Expanded outlier transformer demos in Jupyter notebooks (**Ashok Kumar**)
    - Expanded Pipeline demos in Jupyter notebooks.

Community
~~~~~~~~~
    - Created Gitter community to support users and foster knowledge exchange


Version 1.0.2
-------------

Deployed: 22th January 2021

Contributors
~~~~~~~~~~~~
    - Nicolas Galli
    - Pradumna Suryawanshi
    - Elamraoui Sohayb
    - Soledad Galli

New transformers
~~~~~~~~~~~~~~~~
    - **CombineWithReferenceFeatures**: applies mathematical operations between a group of variables and reference variables (**by Nicolas Galli**)
    - **DropMissingData**: removes missing observations from a dataset (**Pradumna Suryawanshi**)

Bug Fix
~~~~~~~
    - Fix bugs in SelectByTargetMeanPerformance.
    - Fix documentation and jupyter notebook typos.

Tutorials
~~~~~~~~~

    - **Creation**: updated "how to" examples on how to combine variables into new features (**by Elamraoui Sohayb and Nicolas Galli**)
    - **Kaggle Kernels**: include links to Kaggle kernels


Version 1.0.1
-------------

Deployed: 11th January 2021

Bug Fix
~~~~~~~
    - Fix use of r2 in SelectBySingleFeaturePerformance and SelectByTargetMeanPerformance.
    - Fix documentation not showing properly in readthedocs.


Version 1.0.0
-------------

Deployed: 31st December 2020

Contributors
~~~~~~~~~~~~
    - Ashok Kumar
    - Christopher Samiullah
    - Nicolas Galli
    - Nodar Okroshiashvili
    - Pradumna Suryawanshi
    - Sana Ben Driss
    - Tejash Shah
    - Tung Lee
    - Soledad Galli


In this version, we made a major overhaul of the package, with code quality improvement
throughout the code base, unification of attributes and methods, addition of new
transformers and extended documentation. Read below for more details.

New transformers for Feature Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We included a whole new module with multiple transformers to select features.

    - **DropConstantFeatures**: removes constant and quasi-constant features from a dataframe (**by Tejash Shah**)
    - **DropDuplicateFeatures**: removes duplicated features from a dataset (**by Tejash Shah and Soledad Galli**)
    - **DropCorrelatedFeatures**: removes features that are correlated (**by Nicolas Galli**)
    - **SmartCorrelationSelection**: selects feature from group of correlated features based on certain criteria (**by Soledad Galli**)
    - **ShuffleFeaturesSelector**: selects features by drop in machine learning model performance after feature's values are randomly shuffled (**by Sana Ben Driss**)
    - **SelectBySingleFeaturePerformance**: selects features based on a ML model performance trained on individual features (**by Nicolas Galli**)
    - **SelectByTargetMeanPerformance**: selects features encoding the categories or intervals with the target mean and using that as proxy for performance (**by Tung Lee and Soledad Galli**)
    - **RecursiveFeatureElimination**: selects features recursively, evaluating the drop in ML performance, from the least to the most important feature (**by Sana Ben Driss**)
    - **RecursiveFeatureAddition**: selects features recursively, evaluating the increase in ML performance, from the most to the least important feature (**by Sana Ben Driss**)


Renaming of Modules
~~~~~~~~~~~~~~~~~~~

Feature-engine transformers have been sorted into submodules to smooth the development
of the package and shorten import syntax for users.

    - **Module imputation**: missing data imputers are now imported from ``feature_engine.imputation`` instead of ``feature_engine.missing_data_imputation``.
    - **Module encoding**: categorical variable encoders are now imported from ``feature_engine.encoding`` instead of ``feature_engine_categorical_encoders``.
    - **Module discretisation**: discretisation transformers are now imported from ``feature_engine.discretisation`` instead of ``feature_engine.discretisers``.
    - **Module transformation**: transformers are now imported from ``feature_engine.transformation`` instead of ``feature_engine.variable_transformers``.
    - **Module outliers**: transformers to remove or censor outliers are now imported from ``feature_engine.outliers`` instead of ``feature_engine.outlier_removers``.
    - **Module selection**: new module hosts transformers to select or remove variables from a dataset.
    - **Module creation**: new module hosts transformers that combine variables into new features using mathematical or other operations.

Renaming of Classes
~~~~~~~~~~~~~~~~~~~

We shortened the name of categorical encoders, and also renamed other classes to
simplify import syntax.

    - **Encoders**: the word ``Categorical`` was removed from the classes name. Now, instead of ``MeanCategoricalEncoder``, the class is called ``MeanEncoder``. Instead of ``RareLabelCategoricalEncoder`` it is ``RareLabelEncoder`` and so on. Please check the encoders documentation for more details.
    - **Imputers**: the ``CategoricalVariableImputer`` is now called ``CategoricalImputer``.
    - **Discretisers**: the ``UserInputDiscretiser`` is now called ``ArbitraryDiscretiser``.
    - **Creation**: the ``MathematicalCombinator`` is not called ``MathematicalCombination``.
    - **WoEEncoder and PRatioEncoder**: the ``WoEEncoder`` now applies only encoding with the weight of evidence. To apply encoding by probability ratios, use a different transformer: the ``PRatioEncoder`` (**by Nicolas Galli**).

Renaming of Parameters
~~~~~~~~~~~~~~~~~~~~~~

We renamed a few parameters to unify the nomenclature across the Package.

    - **EndTailImputer**: the parameter ``distribution`` is now called ``imputation_method`` to unify convention among imputers. To impute using the IQR, we now need to pass ``imputation_method="iqr"`` instead of ``imputation_method="skewed"``.
    - **AddMissingIndicator**: the parameter ``missing_only`` now takes the boolean values ``True`` or ``False``.
    - **Winzoriser and OutlierTrimmer**: the parameter ``distribution`` is now called ``capping_method`` to unify names across Feature-engine transformers.


Tutorials
~~~~~~~~~

    - **Imputation**: updated "how to" examples of missing data imputation (**by Pradumna Suryawanshi**)
    - **Encoders**: new and updated "how to" examples of categorical encoding (**by Ashok Kumar**)
    - **Discretisation**: new and updated "how to" examples of discretisation (**by Ashok Kumar**)
    - **Variable transformation**: updated "how to" examples on how to apply mathematical transformations to variables (**by Pradumna Suryawanshi**)


For Contributors and Developers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Code Architecture
^^^^^^^^^^^^^^^^^

    - **Submodules**: transformers have been grouped within relevant submodules and modules.
    - **Individual tests**: testing classes have been subdivided into individual tests
    - **Code Style**: we adopted the use of flake8 for linting and PEP8 style checks, and black for automatic re-styling of code.
    - **Type hint**: we rolled out the use of type hint throughout classes and functions (**by Nodar Okroshiashvili, Soledad Galli and Chris Samiullah**)

Documentation
^^^^^^^^^^^^^

    - Switched fully to numpydoc and away from Napoleon
    - Included more detail about methods, parameters, returns and raises, as per numpydoc docstring style (**by Nodar Okroshiashvili, Soledad Galli**)
    - Linked documentation to github repository
    - Improved layout

Other Changes
~~~~~~~~~~~~~

    - **Updated documentation**: documentation reflects the current use of Feature-engine transformers
    - **Typo fixes**: Thank you to all who contributed to typo fixes (Tim Vink, Github user @piecot).. -*- mode: rst -*-

Version 0.6.X
=============

Version 0.6.1
-------------
Deployed: Friday, September 18, 2020

Contributors: Soledad Galli

Minor Changes:
    - **Updated docs**: updated and expanded Contributing guidelines, added Governance, updated references to Feature-engine online.
    - **Updated Readme**: updated and expanded readme.


Version 0.6.0
-------------
Deployed: Friday, August 14, 2020

Contributors: 
    - Michał Gromiec
    - Surya Krishnamurthy
    - Gleb Levitskiy
    - Karthik Kothareddy
    - Richard Cornelius Suwandi
    - Chris Samiullah
    - Soledad Galli


Major Changes:
    - **New Transformer**: the ``MathematicalCombinator`` allows you combine multiple features into new variables by performing mathematical operations like sum, product, mean, standard deviation, or finding the minimum and maximum values (by Michał Gromiec).
    - **New Transformer**: the ``DropFeatures`` allows you remove specified variables from a dataset (by Karthik Kothareddy).
    - **New Transformer**: the ``DecisionTreeCategoricalEncoder`` encodes categorical variables with a decision tree (by Surya Krishnamurthy).
    - **Bug fix**: the ``SklearnTransformerWrapper`` can now automatically select numerical or numerical and categorical variables depending on the Scikit-learn transformer the user implements (by Michał Gromiec).
    - **Bug fix**: the ``SklearnTransformerWrapper`` can now wrap Scikit-learn's OneHotEncoder and concatenate the binary features back to the original dataframe (by Michał Gromiec).
    - **Added functionality**: the ``ArbitraryNumberImputer`` can now take a dictionary of variable, arbitrary number pairs, to impute different variables with different numbers (by Michał Gromiec).
    - **Added functionality**: the ``CategoricalVariableImputer`` can now replace missing data in categorical variables by a string defined by the user (by Gleb Levitskiy).
    - **Added functionality**: the ``RareLabelEnoder`` now allows the user to determine the maximum number of categories that the variable should have when grouping infrequent values (by Surya Krishnamurthy).


Minor Changes:
    - **Improved docs**: fixed typos and tidy Readme.md (by Richard Cornelius Suwandi)
    - **Improved engineering practices**: added Manifest.in to include md and licenses in tar ball in pypi (by Chris Samiullah)
    - **Improved engineering practices**: updated circleci yaml and created release branch for orchestrated release of new versions with significant changes (by Soledad Galli and Chris Samiullah)
    - **Improved engineering practices**: added test for doc build in circleci yaml (by Soledad Galli and Chris Samiullah)
    - **Transformer fix**: removed parameter return_object from the RareLabelEncoder as it was not working as intended(by Karthik Kothareddy and Soledad Galli)


Version 0.5.0
-------------

* Deployed: Friday, July 10, 2020
* Contributors: Soledad Galli

Major Changes:
    - **Bug fix**: fixed error in weight of evidence formula in the ``WoERatioCategoricalEncoder``. The old formula, that is np.log( p(1) / p(0) ) is preserved, and can be obtained by setting the ``encoding_method`` to 'log_ratio'. If ``encoding_method`` is set to 'woe', now the correct formula will operate.
	- **Added functionality**: most categorical encoders have the option ``inverse_transform``, to obtain the original value of the variable from the transformed dataset.
    - **Added functionality**: the `'Winsorizer``, ``OutlierTrimmer`` and ``ArbitraryOutlierCapper`` have now the option to ignore missing values, and obtain the parameters from the original variable distribution, or raise an error if the dataframe contains na, by setting the parameter ``missing_values`` to ``raise`` or ``ignore``.
    - **New Transformer**: the ``UserInputDiscretiser`` allows users to discretise numerical variables into arbitrarily defined buckets.


Version 0.4.3
-------------

* Deployed: Friday, May 15, 2020
* Contributors: Soledad Galli, Christopher Samiullah

Major Changes:
	- **New Transformer**: the `'SklearnTransformerWrapper`` allows you to use most Scikit-learn transformers just on a subset of features. Works with the SimpleImputer, the OrdinalEncoder and most scalers.

Minor Changes:
    - **Added functionality**: the `'EqualFrequencyDiscretiser`` and ``EqualWidthDiscretiser`` now have the ability to return interval boundaries as well as integers, to identify the bins. To return boundareis set the parameter ``return_boundaries=True``.
    - **Improved docs**: added contibuting section, where you can find information on how to participate in the development of Feature-engine's code base, and more.


Version 0.4.0
-------------
* Deployed: Monday, April 04, 2020
* Contributors: Soledad Galli, Christopher Samiullah

Major Changes:
    - **Deprecated**: the ``FrequentCategoryImputer`` was integrated into the class ``CategoricalVariableImputer``. To perform frequent category imputation now use: ``CategoricalVariableImputer(imputation_method='frequent')``
    - **Renamed**: the ``AddNaNBinaryImputer`` is now called ``AddMissingIndicator``.
    - **New**: the ``OutlierTrimmer`` was introduced into the package and allows you to remove outliers from the dataset

Minor Changes:
    - **Improved**: the ``EndTailImputer`` now has the additional option to place outliers at a factor of the maximum value.
    - **Improved**: the ``FrequentCategoryImputer`` has now the functionality to return numerical variables cast as object, in case you want to operate with them as if they were categorical. Set ``return_object=True``.
    - **Improved**: the ``RareLabelEncoder`` now allows the user to define the name for the label that will replace rare categories.
    - **Improved**: All feature engine transformers (except missing data imputers) check that the data sets do not contain missing values.
    - **Improved**: the ``LogTransformer`` will raise an error if a variable has zero or negative values.
    - **Improved**: the ``ReciprocalTransformer`` now works with variables of type integer.
    - **Improved**: the ``ReciprocalTransformer`` will raise an error if the variable contains the value zero.
    - **Improved**: the ``BoxCoxTransformer`` will raise an error if the variable contains negative values.
    - **Improved**: the ``OutlierCapper`` now finds and removes outliers based of percentiles.
    - **Improved**: Feature-engine is now compatible with latest releases of Pandas and Scikit-learn.


Version 0.3.0
-------------
* Deployed: Monday, August 05, 2019
* Contributors: Soledad Galli.

Major Changes:
    - **New**: the ``RandomSampleImputer`` now has the option to set one seed for batch imputation or set a seed observation per observations based on 1 or more additional numerical variables for that observation, which can be combined with multiplication or addition.
    - **New**: the ``YeoJohnsonTransfomer`` has been included to perform Yeo-Johnson transformation of numerical variables.
    - **Renamed**: the  ``ExponentialTransformer`` is now called ``PowerTransformer``.
    - **Improved**: the ``DecisionTreeDiscretiser`` now allows to provide a grid of parameters to tune the decision trees which is done with a GridSearchCV under the hood.
    - **New**: Extended documentation for all Feature-engine's transformers.
    - **New**:  *Quickstart* guide to jump on straight onto how to use Feature-engine.
    - **New**: *Changelog* to track what is new in Feature-engine.
    - **Updated**: new ``Jupyter notebooks`` with examples on how to use Feature-engine's transformers.

Minor Changes:
    - **Unified**: dictionary attributes in transformers, which contain the transformation mappings, now end with ``_``, for example ``binner_dict_``.
Version 1.2.X
=============

Version 1.2.0
-------------

Deployed: 4th January 2022

Contributors
~~~~~~~~~~~~

    - `Edoardo Argiolas <https://github.com/dodoarg>`_
    - `gverbock <https://github.com/gverbock>`_
    - `Thibault Blanc <https://github.com/thibaultbl>`_
    - `David Cortes <https://github.com/david-cortes>`_
    - `Morgan Sell <https://github.com/Morgan-Sell>`_
    - `Kevin Kurek <https://github.com/kevinkurek>`_
    - `Soledad Galli <https://github.com/solegalli>`_

In this big release, we add 3 new transformers, we expand the functionality of existing
classes, we add information about citing Feature-engine and we expand the documentation
with a new look, extended user guide with examples, and more details on how to
contribute to the project.

Thank you so much to the contributors for making this massive release possible!

Thank you to reviewers `Nicolas Galli <https://github.com/nicogalli>`_ and
`Chris Samiullah <https://github.com/christophergs>`_ for useful advice on
various PRs.

New transformers
~~~~~~~~~~~~~~~~

    - **DatetimeFeatures**: extracts date and time features from datetime variables (`Edoardo Argiolas <https://github.com/dodoarg>`_)
    - **DropHishPSIFeatures**: finds and drops features with high population stability index (`gverbock <https://github.com/gverbock>`_)
    - **Matchvariables**: ensures that the same variables observed in the train set are present in the test set (`Thibault Blanc <https://github.com/thibaultbl>`_)

Enhancements
~~~~~~~~~~~~

    - The **Winsorizer** can now add binary variable indicators to flag outlier values (`David Cortes <https://github.com/david-cortes>`_)
    - The **DropMissingData** now allows to drop rows based on % of missing data (`Kevin Kurek <https://github.com/kevinkurek>`_)
    - **Categorical encoders** can now either raise a warning or an error when encoding categories not seen in the train set (`Morgan Sell <https://github.com/Morgan-Sell>`_)
    - The **ArbitraryDiscretiser** can now either raise a warning or an error when values fall outside the limits entered by the user (`Morgan Sell <https://github.com/Morgan-Sell>`_)
    - **CombineWithReferenceFeature** and **MathematicalCombination** have now the option to drop the original input variables after the feature creation (`Edoardo Argiolas <https://github.com/dodoarg>`_)

Bug fixes
~~~~~~~~~

    - All **Encoders** are now able to exclude datetime variables cast as object or categorical when searching for categorical variables automatically (`Edoardo Argiolas <https://github.com/dodoarg>`_)
    - All transformers will now raise an error when users pass an empty list to the variables parameter (`Edoardo Argiolas <https://github.com/dodoarg>`_)
    - All transformers now check the variable type when user passes a single variable to the variables parameter (`Edoardo Argiolas <https://github.com/dodoarg>`_)


Documentation
~~~~~~~~~~~~~
    - We changed the template to pydata (`Soledad Galli <https://github.com/solegalli>`_)
    - We split the information about transformers into a user guide and an API (`Soledad Galli <https://github.com/solegalli>`_)
    - The API documentation shows how to use the transformers (`Soledad Galli <https://github.com/solegalli>`_)
    - The user guide expands the API docs with plenty of examples and tips on when and how to use the transformers (`Soledad Galli <https://github.com/solegalli>`_)
    - We expanded the contribute section with plenty of details on how to make a contribution and how to check your code is top notch (`Soledad Galli <https://github.com/solegalli>`_)
    - You can now sponsor Feature-engine (`Soledad Galli <https://github.com/solegalli>`_)
    - You can now cite our JOSS article when using Feature-engine (`Soledad Galli <https://github.com/solegalli>`_)
    - We added plenty of examples on how to use the new class DropHighPSIFeatures (`gverbock <https://github.com/gverbock>`_)
    - We included various examples on how to extract date and time features using the new DatetimeFeatures class (`Edoardo Argiolas <https://github.com/dodoarg>`_)
    - We included examples on how to use the new class MatchVariables (`Thibault Blanc <https://github.com/thibaultbl>`_)
    - We added a Jupyter notebook with a demo of the new DatetimeFeatures class (`Edoardo Argiolas <https://github.com/dodoarg>`_)
    - We added a Jupyter notebook with a demo of the new DropHighPSIFeatures class (`Soledad Galli <https://github.com/solegalli>`_)
.. -*- mode: rst -*-

What's new
==========

Find out what's new in each new version release.

.. toctree::
   :maxdepth: 2

   v_120
   v_1
   v_06.. _api:

API
===

Full API documentation for Feature-engine transformers.

.. toctree::
   :maxdepth: 1

   imputation/index
   encoding/index
   discretisation/index
   outliers/index
   transformation/index
   creation/index
   selection/index
   datetime/index
   preprocessing/index
   wrappers/indexMathematicalCombination
=======================

.. autoclass:: feature_engine.creation.MathematicalCombination
    :members:

CombineWithReferenceFeature
===========================

.. autoclass:: feature_engine.creation.CombineWithReferenceFeature
    :members:

CyclicalTransformer
===================

.. autoclass:: feature_engine.creation.CyclicalTransformer
    :members:
.. -*- mode: rst -*-

Feature Creation
================

Feature-engine's creation transformers create and add new features to the dataframe
by either combining or transforming existing features.

.. toctree::
   :maxdepth: 2

   MathematicalCombination
   CombineWithReferenceFeature
   CyclicalTransformer

Transformers in other Libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check also the following transformer from Scikit-learn:

* `PolynomialFeatures <https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.PolynomialFeatures.html#sklearn.preprocessing.PolynomialFeatures>`_
ArbitraryNumberImputer
======================

.. autoclass:: feature_engine.imputation.ArbitraryNumberImputer
    :members:

DropMissingData
===============

.. autoclass:: feature_engine.imputation.DropMissingData
    :members:

MeanMedianImputer
=================

.. autoclass:: feature_engine.imputation.MeanMedianImputer
    :members:

CategoricalImputer
==================

.. autoclass:: feature_engine.imputation.CategoricalImputer
    :members:

EndTailImputer
==============

.. autoclass:: feature_engine.imputation.EndTailImputer
    :members:

.. -*- mode: rst -*-

.. currentmodule:: feature_engine.imputation

Missing Data Imputation
=======================

Feature-engine's missing data imputers replace missing data by parameters estimated
from data or arbitrary values pre-defined by the user.


**Summary of Feature-engine's imputers main characteristics**

================================== ===================== ======================= ====================================================================================
    Transformer                     Numerical variables	  Categorical variables	    Description
================================== ===================== ======================= ====================================================================================
:class:`MeanMedianImputer()`	        √	                 ×	                    Replaces missing values by the mean or median
:class:`ArbitraryNumberImputer()`	    √	                 x	                    Replaces missing values by an arbitrary value
:class:`EndTailImputer()`	            √	                 ×	                    Replaces missing values by a value at the end of the distribution
:class:`CategoricalImputer()`           √	                 √	                    Replaces missing values by the most frequent category or by an arbitrary value
:class:`RandomSampleImputer()`	        √	                 √	                    Replaces missing values by random value extractions from the variable
:class:`AddMissingIndicator()`	        √	                 √	                    Adds a binary variable to flag missing observations
:class:`DropMissingData()`	            √	                 √	                    Removes observations with missing data from the dataset
================================== ===================== ======================= ====================================================================================


The :class:`CategoricalImputer()` performs procedures suitable for categorical variables. From
version 1.1.0 it also accepts numerical variables as input, for those cases were
categorical variables by nature are coded as numeric.


.. toctree::
   :maxdepth: 2
   :hidden:

   MeanMedianImputer
   ArbitraryNumberImputer
   EndTailImputer
   CategoricalImputer
   RandomSampleImputer
   AddMissingIndicator
   DropMissingDataAddMissingIndicator
===================

.. autoclass:: feature_engine.imputation.AddMissingIndicator
    :members:

RandomSampleImputer
===================

.. autoclass:: feature_engine.imputation.RandomSampleImputer
    :members:

Winsorizer
==========

.. autoclass:: feature_engine.outliers.Winsorizer
    :members:

ArbitraryOutlierCapper
======================

.. autoclass:: feature_engine.outliers.ArbitraryOutlierCapper
    :members:
.. -*- mode: rst -*-

.. currentmodule:: feature_engine.outliers

Outlier Handling
================

Feature-engine's outlier transformers cap maximum or minimum values of a variable at an
arbitrary or derived value. The OutlierTrimmer removes outliers from the dataset.

=================================== ==============================================================
 Transformer                          Description
=================================== ==============================================================
:class:`Winsorizer()`                 Caps variables at automatically determined extreme values
:class:`ArbitraryOutlierCapper()`     Caps variables at values determined by the user
:class:`OutlierTrimmer()`             Removes outliers from the dataframe
=================================== ==============================================================

.. toctree::
   :maxdepth: 2
   :hidden:

   Winsorizer
   ArbitraryOutlierCapper
   OutlierTrimmerOutlierTrimmer
==============

.. autoclass:: feature_engine.outliers.OutlierTrimmer
    :members:

WoEEncoder
==========

.. autoclass:: feature_engine.encoding.WoEEncoder
    :members:

DecisionTreeEncoder
===================

.. autoclass:: feature_engine.encoding.DecisionTreeEncoder
    :members:

RareLabelEncoder
================


.. autoclass:: feature_engine.encoding.RareLabelEncoder
    :members:

OrdinalEncoder
==============

.. autoclass:: feature_engine.encoding.OrdinalEncoder
    :members:

MeanEncoder
===========

.. autoclass:: feature_engine.encoding.MeanEncoder
    :members:

CountFrequencyEncoder
=====================

.. autoclass:: feature_engine.encoding.CountFrequencyEncoder
    :members:
OneHotEncoder
=============

.. autoclass:: feature_engine.encoding.OneHotEncoder
    :members:

.. -*- mode: rst -*-

Categorical Encoding
====================

Feature-engine's categorical encoders replace variable strings by estimated or
arbitrary numbers.

.. figure::  ../../images/summary/categoricalSummary.png
   :align:   center

   Summary of Feature-engine's encoders main characteristics

Feature-engine's categorical encoders work only with categorical variables by default.
From version 1.1.0, you have the option to set the parameter `ignore_format` to True,
and make the transformers also accept numerical variables as input.


.. toctree::
   :maxdepth: 2

   OneHotEncoder
   CountFrequencyEncoder
   OrdinalEncoder
   MeanEncoder
   WoEEncoder
   PRatioEncoder
   DecisionTreeEncoder
   RareLabelEncoder

Other categorical encoding libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For additional categorical encoding transformations, visit the open-source package
`Category encoders <https://contrib.scikit-learn.org/category_encoders/>`_.
PRatioEncoder
=============

.. autoclass:: feature_engine.encoding.PRatioEncoder
    :members:

SklearnTransformerWrapper
=========================

.. autoclass:: feature_engine.wrappers.SklearnTransformerWrapper
    :members:

.. -*- mode: rst -*-

.. currentmodule:: feature_engine.wrappers

Scikit-learn Wrapper
====================

Feature-engine's Scikit-learn wrappers wrap Scikit-learn transformers allowing their
implementation only on a selected subset of features.

.. toctree::
   :maxdepth: 2

   Wrapper

Other wrappers
~~~~~~~~~~~~~~

The :class:`SklearnTransformerWrapper()` offers a similar function to the
`ColumnTransformer <https://scikit-learn.org/stable/modules/generated/sklearn.compose.ColumnTransformer.html>`_
class available in Scikit-learn. They differ in the implementation to select the
variables.EqualWidthDiscretiser
=====================

.. autoclass:: feature_engine.discretisation.EqualWidthDiscretiser
    :members:
ArbitraryDiscretiser
====================

.. autoclass:: feature_engine.discretisation.ArbitraryDiscretiser
    :members:

DecisionTreeDiscretiser
=======================

.. autoclass:: feature_engine.discretisation.DecisionTreeDiscretiser
    :members:
.. -*- mode: rst -*-
.. currentmodule:: feature_engine.discretisation

Variable Discretisation
=======================

Feature-engine's discretisation transformers transform continuous variables into
discrete features. This is accomplished, in general, by sorting the variable values
into continuous intervals.

**Summary**

=====================================  ========================================================================
      Transformer                           Functionality
=====================================  ========================================================================
:class:`EqualFrequencyDiscretiser()`     Sorts values into intervals with similar number of observations.
:class:`EqualWidthDiscretiser()`         Sorts values into intervals of equal size.
:class:`ArbitraryDiscretiser()`          Sorts values into intervals predefined by the user.
:class:`DecisionTreeDiscretiser()`       Replaces values by predictions of a decision tree, which are discrete
=====================================  ========================================================================


.. toctree::
   :maxdepth: 2
   :hidden:

   EqualFrequencyDiscretiser
   EqualWidthDiscretiser
   ArbitraryDiscretiser
   DecisionTreeDiscretiser

Additional transformers for discretisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For discretisation using K-means, check Scikit-learn's
`KBinsDiscretizer <https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.KBinsDiscretizer.html>`_.EqualFrequencyDiscretiser
=========================

.. autoclass:: feature_engine.discretisation.EqualFrequencyDiscretiser
    :members:

DropCorrelatedFeatures
======================

.. autoclass:: feature_engine.selection.DropCorrelatedFeatures
    :members:
RecursiveFeatureAddition
========================


.. autoclass:: feature_engine.selection.RecursiveFeatureAddition
    :members:

SelectBySingleFeaturePerformance
================================


.. autoclass:: feature_engine.selection.SelectBySingleFeaturePerformance
    :members:
SelectByShuffling
=================

.. autoclass:: feature_engine.selection.SelectByShuffling
    :members:
DropFeatures
=============

.. autoclass:: feature_engine.selection.DropFeatures
    :members:

DropHighPSIFeatures
===================


.. autoclass:: feature_engine.selection.DropHighPSIFeatures
    :members:DropConstantFeatures
====================

.. autoclass:: feature_engine.selection.DropConstantFeatures
    :members:

SmartCorrelatedSelection
========================


.. autoclass:: feature_engine.selection.SmartCorrelatedSelection
    :members:

RecursiveFeatureElimination
============================


.. autoclass:: feature_engine.selection.RecursiveFeatureElimination
    :members:

SelectByTargetMeanPerformance
=============================


.. autoclass:: feature_engine.selection.SelectByTargetMeanPerformance
    :members:

DropDuplicateFeatures
=====================


.. autoclass:: feature_engine.selection.DropDuplicateFeatures
    :members:

.. -*- mode: rst -*-
.. currentmodule:: feature_engine.selection

Feature Selection
=================

Feature-engine's feature selection transformers are used to drop subsets of variables,
or in other words, to select subsets of variables. Feature-engine hosts selection
algorithms that are, in general, not available in other libraries. These algorithms have
been gathered from data science competitions or used in the industry.

Feature-engine's transformers select features based on 2 strategies. They either select
features by looking at the features intrinsic characteristics, like distributions or their
relationship with other features. Or they select features based on their impact on the
machine learning model performance.

In the following tables you find the algorithms that belong to either category.

Selection based on feature characteristics
------------------------------------------

============================================ ======================= ============= ====================================================================================
    Transformer                                Categorical variables   Allows NA	    Description
============================================ ======================= ============= ====================================================================================
:class:`DropFeatures()`                         √	                      √	            Drops arbitrary features determined by user
:class:`DropConstantFeatures()`  	            √	                      √	            Drops constant and quasi-constant features
:class:`DropDuplicateFeatures()`                √	                      √             Drops features that are duplicated
:class:`DropCorrelatedFeatures()`               ×	                      √	            Drops features that are correlated
:class:`SmartCorrelatedSelection()`	            ×	                      √	            From a correlated feature group drops the less useful features
:class:`DropHighPSIFeatures()`	                ×	                      √	            Drops features with high Population Stability Index
============================================ ======================= ============= ====================================================================================

Selection based on model performance
------------------------------------

============================================ ======================= ============= ====================================================================================
    Transformer                                Categorical variables   Allows NA	    Description
============================================ ======================= ============= ====================================================================================
:class:`SelectByShuffling()`	                ×	                      ×	            Selects features if shuffling their values causes a drop in model performance
:class:`SelectBySingleFeaturePerformance()`	    ×	                      ×	            Removes observations with missing data from the dataset
:class:`SelectByTargetMeanPerformance()`        √                         ×             Using the target mean as performance proxy, selects high performing features
:class:`RecursiveFeatureElimination()`          ×                         ×             Removes features recursively by evaluating model performance
:class:`RecursiveFeatureAddition()`             ×                         ×             Adds features recursively by evaluating model performance
============================================ ======================= ============= ====================================================================================


.. toctree::
   :maxdepth: 2
   :hidden:

   DropFeatures
   DropConstantFeatures
   DropDuplicateFeatures
   DropCorrelatedFeatures
   SmartCorrelatedSelection
   DropHighPSIFeatures
   SelectByShuffling
   SelectBySingleFeaturePerformance
   SelectByTargetMeanPerformance
   RecursiveFeatureElimination
   RecursiveFeatureAddition

Other Feature Selection Libraries
---------------------------------

For additional feature selection algorithms visit the following open-source libraries:

* `Scikit-learn selection <https://scikit-learn.org/stable/modules/feature_selection.html>`_
* `MLXtend selection <http://rasbt.github.io/mlxtend/api_subpackages/mlxtend.feature_selection/>`_

Scikit-learn hosts multiple filter and embedded methods that select features based on
statistical tests or machine learning model derived importance. MLXtend hosts greedy
(wrapper) feature selection methods.
LogCpTransformer
================

.. autoclass:: feature_engine.transformation.LogCpTransformer
    :members:
YeoJohnsonTransformer
=====================

.. autoclass:: feature_engine.transformation.YeoJohnsonTransformer
    :members:

PowerTransformer
================


.. autoclass:: feature_engine.transformation.PowerTransformer
    :members:

LogTransformer
==============


.. autoclass:: feature_engine.transformation.LogTransformer
    :members:

ReciprocalTransformer
=====================


.. autoclass:: feature_engine.transformation.ReciprocalTransformer
    :members:

BoxCoxTransformer
=================

.. autoclass:: feature_engine.transformation.BoxCoxTransformer
    :members:

.. -*- mode: rst -*-

Variable Transformation
=======================

Feature-engine's variable transformers transform numerical variables with various
mathematical transformations.

.. toctree::
   :maxdepth: 2

   LogTransformer
   LogCpTransformer
   ReciprocalTransformer
   PowerTransformer
   BoxCoxTransformer
   YeoJohnsonTransformer


Transformers in other Libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These and additional transformations can be obtained with the following Scikit-learn
classes:

* `FunctionTransformer <https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.FunctionTransformer.html>`_
* `PowerTransformer <https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.PowerTransformer.html>`_

Note that Scikit-klearn classes return Numpy arrays and are applied to the entire dataset.
MatchVariables
==============

.. autoclass:: feature_engine.preprocessing.MatchVariables
    :members:

.. -*- mode: rst -*-

Preprocessing
=============

Feature-engine's preprocessing transformers apply general data pre-processing
and transformation procedures.

.. toctree::
   :maxdepth: 2

   MatchVariables
DatetimeFeatures
================

.. autoclass:: feature_engine.datetime.DatetimeFeatures
    :members:

.. -*- mode: rst -*-

Datetime Features
=================

Feature-engine's datetime transformers are able to extract a wide variety of datetime
features from existing datetime or object-like data.

.. toctree::
   :maxdepth: 2

   DatetimeFeatures

.. _featureengine_blogs:

Blogs, Videos and More
======================

Here you find articles and videos about Feature-engine and feature engineering and
selection in general.

Blogs
-----

- `Feature-engine: A new open-source Python package for feature engineering <https://trainindata.medium.com/feature-engine-a-new-open-source-python-package-for-feature-engineering-29a0ab88ea7c/>`_.
- `Practical Code Implementations of Feature Engineering for Machine Learning with Python <https://towardsdatascience.com/practical-code-implementations-of-feature-engineering-for-machine-learning-with-python-f13b953d4bcd>`_.
- `Streamlining Feature Engineering Pipelines with Feature-engine <https://towardsdatascience.com/streamlining-feature-engineering-pipelines-with-feature-engine-e781d551f470?gi=e0fa6e5c0c1a/>`_.
- `Feature Engineering for Machine Learning: A comprehensive Overview <https://trainindata.medium.com/feature-engineering-for-machine-learning-a-comprehensive-overview-a7ad04c896f8>`_.
- `Feature Selection for Machine Learning: A comprehensive Overview <https://trainindata.medium.com/feature-selection-for-machine-learning-a-comprehensive-overview-bd571db5dd2d>`_.


Videos
------

- `Optimising Feature Engineering Pipelines with Feature-engine <https://www.youtube.com/watch?v=qT-3KUaFYmk/>`_, Pydata Cambridge 2020, from minute 51:43.


Podcasts
--------

- `Build Composable And Reusable Feature Engineering Pipelines with Feature-Engine <https://www.pythonpodcast.com/feature-engine-feature-engineering-pipelines-episode-338/>`_.

En Español
----------

- `Ingeniería de variables para machine learning <https://www.udemy.com/course/ingenieria-de-variables-para-machine-learning/?referralCode=CE398C784F17BD87482C>`_, Curso Online.
- `Ingeniería de variables, MachinLenin <https://www.youtube.com/watch?v=NhCxOOoFXds>`_, charla con video online.

More resources will be added as they appear online. If you know of a good resource, let us know.
Books
=====

You can learn more about how to use Feature-engine and feature engineering in general
in the following books:

.. figure::  ../images/cookbook.png
   :width: 200
   :figclass: align-center
   :align: left
   :target: https://packt.link/python

   Python Feature Engineering CookbookTutorials
=========

How To
------

Check our `jupyter notebooks <https://nbviewer.jupyter.org/github/feature-engine/feature-engine-examples/tree/main/>`_
showcasing the functionality of each Feature-engine transformer.

Kaggle Kernels
--------------

We also prepared Kaggle kernels with demos mixing data exploration, feature engineering,
feature creation, feature selection and hyperparameter optimization of entire pipelines.

- `Feature selection for bank customer satisfaction prediction <https://www.kaggle.com/solegalli/feature-selection-with-feature-engine>`_
- `Feature engineering and selection for house price prediction <https://www.kaggle.com/solegalli/predict-house-price-with-feature-engine>`_
- `Feature creation for wine quality prediction <https://www.kaggle.com/solegalli/create-new-features-with-feature-engine>`_
- `Feature engineering and model stacking for house price modelling <https://www.kaggle.com/solegalli/feature-engineering-and-model-stacking>`_
- `Feature engineering with Feature-engine and Randomized search <https://www.kaggle.com/solegalli/feature-engineering-with-randomized-search>`_
- `Feature engineering with Feature-engine and Grid search <https://www.kaggle.com/solegalli/feature-engineering-pipeline-and-hyperparam-tuning>`_



Video tutorials
---------------

You can find some videos on how to use Feature-engine in the
`Feature-engine playlist <https://www.youtube.com/playlist?list=PL_7uaHXkQmKVlqlvgQJuaWEKjagHbERtp>`_
in Train in Data's YouTube channel. The list is a bit short at the moment, apologies... -*- mode: rst -*-
.. _learning_resources:

Resources
=========

Here you find learning resources to know more about Feature-engine and feature
engineering and selection in general.

We have gathered online courses, books, blogs, videos, podcasts, jupyter notebook and
kaggle kernels, so you can follow the resource with the way of learning that you like
the most.

.. toctree::
   :maxdepth: 1

   courses
   books
   blogs
   tutorialsCourses
=======

You can learn more about how to use Feature-engine and, feature engineering and feature
selection in general in the following online courses:

.. figure::  ../images/feml.png
   :width: 200
   :figclass: align-center
   :align: left
   :target: https://www.udemy.com/course/feature-engineering-for-machine-learning/?referralCode=A855148E05283015CF06

   Feature Engineering for Machine Learning

.. figure::  ../images/fsml.png
   :width: 200
   :figclass: align-center
   :align: center
   :target: https://www.udemy.com/course/feature-selection-for-machine-learning/?referralCode=186501DF5D93F48C4F71

   Feature Selection for Machine Learning

|
|

We also show how to put an entire machine learning pipeline in production, using
Feature-engine and Scikit-learn in this course:


.. figure::  ../images/dmlm.png
   :width: 200
   :figclass: align-center
   :align: left
   :target: https://www.udemy.com/course/deployment-of-machine-learning-models/?referralCode=D4FE5EA129FFD203CFF4

   Deployment of Machine Learning Models

|
|
|
|
|
|
|
|
|
|

Note that Soledad Galli, the main developer of Feature-engine will earn some money if you
choose to enroll in one of the above courses.
.. -*- mode: rst -*-

Other ways to contribute
========================

A common misconception about contributing to open source is that you need to contribute code.
Equally important to the code contributions are contributions to the documentation or to
the gallery of examples on how to use the project, which we already discussed.

However, a good project does no-one any good if people don't know about it. So spreading the word
about Feature-engine is extremely important. You will do the project a big favor if you help
us spread the word about Feature-engine.

Here are some examples of how you could do that:

Spread the word
---------------

1. Star `Feature-engine's Repository <https://github.com/feature-engine/feature_engine>`_.
2. Share Feature-engine or any of our :ref:`blogs and videos <featureengine_blogs>` on social media.
3. Write a blog about Feature-engine.
4. Give a talk about Feature-engine or mention it in one of your talks.
5. If you teach, use Feature-engine in your lectures.
6. Share Feature-engine with your colleagues.

If you write a blog or give a talk that is publicly available, and it features Feature-engine,
let us know so we link it to the project or share it on social.

If you teach Feature-engine, give a shout-out on Twiter or LinkedIn and tag the maintainers.

If you have other ideas on how to spread the word, let us know or update this page
straightaway.

We are also happy to talk about the project in talks and podcasts. If you host a meetup
or a podcast channel, do get in touch.

Donate
------

If Feature-engine is helping your organization save money, time, effort, and frustrations,
and be more productive, consider `sponsoring us <https://github.com/sponsors/solegalli>`_.
This is a great way to show your appreciation. Thank you!

For all of this, get in touch with the maintainer via `LinkedIn <https://www.linkedin.com/in/soledad-galli/>`_.
.. -*- mode: rst -*-

Contribute Jupyter notebooks
============================

We created a collection of Jupyter notebooks that showcase the main functionality of
Feature-engine's transformers. We link these notebooks throughout the main documentation
to offer users more examples and details about transformers and how to use them.

**Note** that the Jupyter notebooks are hosted in a separate
`Github repository <https://github.com/feature-engine/feature-engine-examples>`_.

Here are some guidelines on how to add a new notebook or update an existing one. The
contribution workflow is the same we use for the main source code base.

Jupyter contribution workflow
-----------------------------

1. Fork the `Github repository <https://github.com/feature-engine/feature-engine-examples>`_.
2. Clone your fork into your local computer: `git clone https://github.com/<YOURUSERNAME>/feature-engine-examples.git`.
3. Navigate into the project directory: `cd feature-engine-examples`.
4. If you haven't done so yet, install feature-engine: `pip install feature_engine`.
5. Create a feature branch with a meaningful name: `git checkout -b mynotebookbranch`.
6. Develop your notebook
7. Add the changes to your copy of the fork: `git add .`, `git commit -m "a meaningful commit message"`, `git pull origin mynotebookbranch`.
8. Go to your fork on Github and make a PR to this repo
9. Done

The review process for notebooks is usually much faster than for the main source code base.

Jupyter creation guidelines
---------------------------

If you want to add a new Jupyter notebook, there are a few things to note:

- Make sure that the dataset you use is publicly available and with a clear license that it is free to use
- Do not upload datasets to the repository
- Add instructions on how to obtain and prepare the data for the demo
- Throughout the notebook, add guidelines on what you are going to do next, and what is the conclusion of the output

That's it! Fairly straightforward.

We look forward to your contribution :).. -*- mode: rst -*-
.. _contribute:

Contribute
==========

Feature-engine is an open source project, originally designed to support the online
course `Feature Engineering for Machine Learning in Udemy <https://www.udemy.com/feature-engineering-for-machine-learning/?couponCode=FEATENGREPO>`_,
but has now gained popularity and supports transformations beyond those taught in the
course.

Feature-engine is currently supported by a small community and we will be delighted to
accept contributions, large or small, that you wish to make to the project.
Contributing to open-source is a great way to learn and improve coding skills, and also
a fun thing to do. If you've never contributed to an open source project, we hope to
make it easy for you with the following guidelines.

Read more about `"Why Contribute to Open-Source" <https://opensource.guide/how-to-contribute/#why-contribute-to-open-source>`_.

Ways to contribute
------------------

There are many ways to contribute to Feature-engine:

- Create a new transformer
- Enhance functionality of current transformers
- Fix a bug
- If you find a bug, let us know by creating an `issues <https://github.com/feature-engine/feature_engine/issues/>`_ on Github.
- If you would like additional functionality or a new feature, create an `issue <https://github.com/feature-engine/feature_engine/issues/>`_ on Github.
- Add a Jupyter notebook to our `Jupyter notebooks example gallery <https://github.com/feature-engine/feature_engine/tree/master/examples>`_.
- Improve our documentation, i.e., fix typos, improve grammar, or add more code examples.
- Write a blog, tweet, or share our project with others.
- Use Feature-engine in your lectures if you teach
- `Sponsor us <https://github.com/sponsors/solegalli>`_.

With plenty of ways to get involved, we would be happy for you to support the project.
You only need to abide by the principles of openness, respect, and consideration of
others, as described in the
`Python Software Foundation Code of Conduct <http://www.python.org/psf/codeofconduct/>`_
and you are ready to go!.

Read more about `"Ways to Contribute to Open Source" <https://opensource.guide/how-to-contribute/#what-it-means-to-contribute>`_.

Getting in touch
----------------

We prefer to handle most contributions through the github repository. You can also join
our Gitter community.

1. `Open issues <https://github.com/feature-engine/feature_engine/issues/>`_.
2. `Gitter community <https://gitter.im/feature_engine/community>`_.


Contributing Guide
------------------

.. toctree::
    :maxdepth: 2

    contribute_code
    contribute_docs
    contribute_jup
    contribute_other
    code_of_conduct
.. -*- mode: rst -*-

.. _contribute_docs:

Contribute Docs
===============

If you contribute a new transformer, or enhance the functionality of a current transformer,
most likely, you would have to add or update the documentation as well.

This is Feature-engine's documentation ecosystem:

- Feature-engine documentation is built using `Sphinx <https://www.sphinx-doc.org>`_ and is hosted on `Read the Docs <https://readthedocs.org/>`_.
- We use the `pydata sphinx theme <https://pypi.org/project/pydata-sphinx-theme/>`_.
- We follow `PEP 257 <https://www.python.org/dev/peps/pep-0257/>`_ for doscstring conventions and use `numpydoc docstring style <https://numpydoc.readthedocs.io/en/latest/format.html>`_.
- All documentation files are located within the `docs folder <https://github.com/feature-engine/feature_engine/tree/main/docs>`_ in the repository.

To learn more about Sphinx check the `Sphinx Quickstart documentation <https://www.sphinx-doc.org/en/master/usage/quickstart.html>`_.

Documents organisation
----------------------

Feature-engine has just adopted Scikit-learn's documentation style, were we offer API
documentation, as well as, a User Guide with examples on how to use the different transformers.

The API documentation is built directly from the docstrings from each transformer. If you
are adding a new transformer, you need to reference it in a new rst file placed within
the `api_doc folder <https://github.com/feature-engine/feature_engine/tree/main/docs/api_doc>`_.

If you would like to add additional examples, you need to update the rst files located
in the `user_guide folder <https://github.com/feature-engine/feature_engine/tree/main/docs/user_guide>`_.

Docstrings
----------

The quickest way to get started with writing the transformer docstrings, is too look at the docstrings
of some of the classes we already have in Feature-engine. Then simply copy and paste
those docstrings and edit the bits that you need. If you copy and paste, make sure to delete
irrelevant parameters and methods.

Link a new transformer
----------------------

If you coded a new transformer from scratch, you need to update the following files to make
sure users can find information on how to use the class correctly:

Add the name of your transformer in these files:

- `Readme <https://github.com/feature-engine/feature_engine/blob/main/README.md>`_.
- `main index <https://github.com/feature-engine/feature_engine/blob/main/docs/index.rst>`_.
- `api index <https://github.com/feature-engine/feature_engine/tree/main/docs/api_doc/index.rst>`_.
- `user guide index <https://github.com/feature-engine/feature_engine/tree/main/docs/user_guide/index.rst>`_.

Add an rst file with the name of your transformer in these folders:

- `api_doc folder <https://github.com/feature-engine/feature_engine/tree/main/docs/api_doc>`_.
- `user_guide folder <https://github.com/feature-engine/feature_engine/tree/main/docs/user_guide>`_.

That's it!

Expand the User Guide
---------------------

You can add more examples or more details to our current User Guide examples. First, find
the relevant rst file for the transformer you would like to work with. Feel free to add more
details on the description of the method, expand the code showcasing other parameters or
whatever you see fit.

We normally run the code on jupyter notebooks, and then copy and paste the code and the
output in the rst files.

Build the documentation
-----------------------

To build the documentation, make sure you have properly installed Sphinx and the required
dependencies. If you set up the development environment as we described in the
:ref:`contribute code guide <contribute_code>`, you should have those installed already.

Alternatively, first activate your environment. Then navigate to the root folder of
Feature-engine. And now install the requirements for the documentation::

        $ pip install -r docs/requirements.txt

To build the documentation (and test if it is working properly) run::

    $ sphinx-build -b html docs build

This command tells sphinx that the documentation files are within the "docs" folder, and
the html files should be placed in the "build" folder.

If everything worked fine, you can open the html files located in build using your browser.
Alternatively, you need to troubleshoot through the error messages returned by sphinx.

Good luck and get in touch if stuck!Code of Conduct
===============

Feature-engine is an open source Python project. We follow the
`Python Software Foundation Code of Conduct <http://www.python.org/psf/codeofconduct/>`_.
All interactions among members of the Feature-engine community must meet those
guidelines. This includes (but is not limited to) interactions through the mailing
list, GitHub and StackOverflow.

Everyone is expected to be open, considerate, and respectful of others no matter what
their position is within the project. We show gratitude for any contribution, big or
small. We welcome feedback and participation. We want to make Feature-engine a nice,
welcoming and safe place for you to do your first contribution to open source, and why
not the second, the third and so on :).
.. -*- mode: rst -*-

.. _contribute_code:

Contribute Code
===============

Contributing code to Feature-engine is fun and easy. If you want to make a code contribution,
you can check the `issue tracker <https://github.com/feature-engine/feature_engine/issues/>`_
for already requested and wanted functionality. Alternatively, you can create a new issue
with functionality you would like to see included in Feature-engine and then work it through.


Contributing workflow
---------------------

A typical contributing workflow goes like this:

1. **Suggest** new functionality or **pick up** an issue from the `issue tracker <https://github.com/feature-engine/feature_engine/issues/>`_.
2. **Mention** in the issue that you are "working on it".
3. **Fork** the repository into your GitHub account.
4. **Clone** your fork into your local computer.
5. Set up the **development environment**.
6. **Create** a new branch with the name of your feature
7. **Code** the feature, the tests and update or add the documentation.
8. **Commit** the changes to your fork.
9. Make a **Pull Request (PR)** with your changes from your fork to the main repo.
10. **Test** the code.
11. **Review** the code with us.
12. Make the **changes** and commit them to your fork, using the same branch created in 5.
13. When it is ready, we will **merge** your contribution into Feature-engine's source code base.

To avoid extra work or duplication, it is important that we communicate right from the
beginning, so we can build together a clear understanding of how you would like to get involved
and what is needed to complete the task. This is particularly important for big code additions.

In the rest of the document, we will describe steps 3 to 13 in more detail.

Fork the Repository
-------------------

When you fork the repository, you create a copy of Feature-engine's source code into
your Github account, which you can edit. To fork Feature-engine's repository, click the
**fork** button in the upper right corner of
`Feature-engine's Repository <https://github.com/feature-engine/feature_engine>`_.

.. figure::  ../images/fork.png
   :figclass: align-center
   :align:   center

Clone the Repository
--------------------

Once you forked the repository, clone the fork to your local machine.

1. Clone your fork into your local machine::

    $ git clone https://github.com/<YOURUSERNAME>/feature_engine

2. Set up an ``upstream`` remote from where you can pull the latest code changes occurring in the main Feature-engine repository::

    $ git remote add upstream https://github.com/feature-engine/feature_engine.git

3. Check that the remote was set correctly::

    $ git remote -v

You should see both your fork (origin) and the main repository (upstream) linked to your local copy::

    origin    https://github.com/YOUR_USERNAME/feature_engine.git (fetch)
    origin    https://github.com/YOUR_USERNAMEfeature_engine.git (push)
    upstream  https://github.com/feature-engine/feature_engine.git (fetch)
    upstream  https://github.com/feature-engine/feature_engine.git (push)

Keep in mind that Feature-engine is being actively developed, so you may need to update
your fork regularly. See below how to **Keep your fork up to date**.

Set up the Development Environment
----------------------------------

After creating a local copy of the repo, you need to set up the development environment.
Setting up a development environment will ensure that you have all the libraries
you need for the development, no more and no less. These libraries include
`Feature-engine dependencies <https://github.com/feature-engine/feature_engine/blob/main/requirements.txt>`_,
like Pandas, NumPy and Scikit-learn and also
`software development libraries <https://github.com/feature-engine/feature_engine/blob/main/test_requirements.txt>`_
like pytest, mypy, flake8, isort and black.

It is optional but highly advisable that you create a virtual environment. A virtual environment
is a "separate space", where you can install Feature-engine's dependencies. To create a virtual
environment, use any virtual environment tool of your choice. Some examples include:

    1. `venv <https://docs.python.org/3/library/venv.html>`_
    2. `conda environments <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_

In the previous links, you find details on how to create the environments as well. We also
provide some guidelines below.

venv
~~~~

If you use venv, from the windows cmd or Mac terminal, create and activate the environment
like this::

    python -m venv /path/to/new/virtual/environment

For example I would do::

    python -m venv Documents/Repositories/envs/featureengine

where "featureengine" will be the name of the environment and "Documents/Repositories/envs"
the location where the environment will be created.

Then, to activate the environment, run::

    Documents/Repositories/envs/featureengine/Scripts/activate

**Note for windows users:** you may need to use \\ instead of /.

conda
~~~~~

If you are using anaconda, from your conda prompt, create and activate the environment
like this::

    conda create --name myenv

where "myenv" will be the name of the environment, so you probably want to change that to
something more meaningful.

Then, to activate the environment, run::

    conda activate myenv


Install dependencies
~~~~~~~~~~~~~~~~~~~~

Now, you are ready to install all dependencies, that is, all the Python libraries used by
Feature-engine. First, navigate to your clone of Feature-engine::

        $ cd feature_engine

Now, install Feature_engine in developer mode::

        $ pip install -e .

Don't forget the `.` after the `-e`. This will add Feature-engine to your PYTHONPATH so your code edits
are automatically picked up, and there is no need to re-install the package after each
code change. This will also install Feature'engine's dependencies.
    
Finally, install the additional dependencies for tests and documentation::

        $ pip install -r test_requirements.txt
        $ pip install -r docs/requirements.txt

Make sure that your local main branch is up to date with the remote main branch::

        $ git pull --rebase upstream main

If you just cloned your fork, your local main branch should be up to date. If you cloned
your fork a time ago, probably the main repository had some code changes. To sync your
fork main branch to the main repository, read below the section **Keep your fork up
to date**.

Create a branch
---------------

It is important to create a new branch, different from main, where you will code your
changes. It is advisable, almost never to work on the main branch.

Create a new branch where you will develop your feature::

    $ git checkout -b myfeaturebranch

where "myfeaturebranch" is the name you choose for your branch.

There are 3 things to keep in mind when creating a feature branch:

1. Give the branch a name that identifies the feature you are going to build.
2. Make sure you checked out your branch from the main branch.
3. Make sure your local main branch was updated with the upstream main branch.

Code your feature
-----------------

Now, you are ready to make your code changes. When you develop a new feature, fix a bug, or
make any code contribution, there are a few things to consider:

1. Make regular code commits to your branch, locally.
2. Give clear messages to your commits, indicating which changes were made at each commit (use present tense).
3. Try and push regularly to your fork, so that you don't lose your changes.

Commit
~~~~~~

Make small changes and commit immediately. This way it is easier to track what was changed.
To commit changes do the following::

    $ git add .
    $ git commit -m "my commit message"

and make sure to include an informative but succinct commit message in the present tense,
for example "fixes style in imputation error message".

The previous commands will commit all files that have changes. If you want to commit just 1
or 2 files, you can do so as follows::

    $ git add file1.py file2.py
    $ git commit -m "my commit message"

It is important that you commit only the files relevant to your feature, and not others
that may have been accidentally changed, for example through code styling (more on this in
**Test the Code** below).

After making a few commits, push your changes to your fork::

    $ git push origin myfeaturebranch

This will automatically create a branch in your remote fork called "myfeaturebranch"
containing all your changes.

Make a Pull Request
~~~~~~~~~~~~~~~~~~~

After pushing the first changes, go to your fork in Github. You will see the branch you
just pushed and next to it a button to create a PR (Pull Request). Go ahead and create a PR from your
feature branch to Feature_engine's **main branch**. In the PR message, describe what the overall
aim of the PR is, and if it resolves an issue, link the issue in the message. This will
notify us of your changes.

Don't worry, you can continue making changes and committing more code to the branch. You
basically need to repeat these steps as often as you need::

    $ git add .
    $ git commit -m "my commit message"
    $ git push origin myfeaturebranch

Once you think your code is ready to review, leave a message in the PR saying "please review"
or something similar.

Create Docstrings
~~~~~~~~~~~~~~~~~

If you are coding an entire new class, make sure you follow our :ref:`guidelines to create
the docstrings <contribute_docs>`.

Test the Code
-------------

The code you submit must pass any test you add plus all current tests in the library.
The tests are triggered automatically when you first make a PR, and then any
time you commit new changes to the PR. It is important that the tests pass when you ask
us for review.

We have tests for:

1. Functionality, using pytest
2. Code style, using flake8
3. Typehints, using mypy
4. Documentation, using sphinx.

In the following paragraphs, we will take you through how to test each of the above individually,
and then altogether.

Test functionality
~~~~~~~~~~~~~~~~~~

We use pytest to create and run our tests. If you set up the development environment as
we described previously, you should have pytest installed. Alternatively, run from the
windows cmd or mac terminal::

    $ pip install pytest

You can now run the tests from your command line interface. Make sure you are within the
feature-engine folder. Then run::

    $ pytest

These command will run all the test scripts within the test folder.

Alternatively, you can run a specific script as follows::

    $ pytest tests/test_encoding/test_categorical_encoder.py

So if you just want to run the code you created, you would do::

    $ pytest tests/test_my_new_feature_folder/test_my_new_feature.py

where test_my_new_feature.py is the name of your test script, and it is located in the
test_my_new_feature_folder.

If you are using Pycharm, this is even easier:

1. In your project directory (where you have all the files and scripts), click with the mouse right button on the folder "tests".
2. Select "Run pytest in tests".
3. Done!!

Sweet, isn't it?

With the above procedure you can also "click" on your individual test script and run only
those tests.


Test Code Style
~~~~~~~~~~~~~~~

We follow `PEP8 <https://pep8.org/>`_ and we keep our code lines up to 88 characters.
Before testing the code style, make sure to automatically
fix anything that might not abide by PEP8 with `**black** <https://pypi.org/project/black/>`_
and `**isort** <https://pypi.org/project/isort/>`_.

If you set up the development environment as we described previously, you should have these
libraries installed. Alternatively, run from the windows cmd or mac terminal::

    $ pip install black
    $ pip install isort

Then, you can sort the imports alphabetically by running::

    $ isort my_new_script.py

You can fix code style by running::

    $ black my_new_script.py


**You need to run isort and black on both code files and test files.**

Black and isort may make changes to your file. Don't forget to commit those changes::

    $ git add my_new_script.py
    $ git commit -m "fixes code styling"
    $ git push origin my_feature_branch

Now, you can go ahead and test that your scripts pass the code styling tests. To do so,
execute from the command line::

    $ flake8 my_new_script.py

If the flake8 test pass, you are good to go. Alternatively, you will get an error, indicating
which line of code is not following the coding convention.

Test Typehint
~~~~~~~~~~~~~

We use `Typehint <https://www.python.org/dev/peps/pep-0484/>`_. To test typehinting we use
`**mypy** <http://mypy-lang.org/>`_.

If you set up the development environment as we described previously, you should have
mypy installed. Alternatively, run from the windows cmd or mac terminal::

    $ pip install mypy

now, you test typehint by running::

    $ mypy feature_engine

A few things to notice:

- We use typehint only on the code base and not on the tests.
- You need to run mypy on the entire module and not just your script.

Otherwise, you will most likely get an error.

Test the docs
~~~~~~~~~~~~~

If after running pytest, black and mypy you do not get errors, you are only left with testing
that the documentation builds correctly.

To do this, first make sure you have all the documentation dependencies installed. If you
set up the environment as we described previously, they should be installed. Alternatively,
from the windows cmd or mac terminal, run::

    $ pip install -r docs/requirements.txt

Make sure you are within the feature_engine module when you run the previous command.

Now, you can go ahead and build the documentation::

    $ sphinx-build -b html docs build

This will trigger the building of the docs, which will be stored in html format in the
"build" folder within the repository. You can open those up with your browser. But the
important thing is that you do not get any red warning during the build process.

Using tox
~~~~~~~~~

In Feature-engine, we use tox to run all our tests automatically. If you want to run all
the tests using tox locally:

1. Install tox in your development environment::

    $ pip install tox

2. Make sure you are in the repository folder, alternatively::

    $ cd feature_engine

3. Run the tests in tox::

    $ tox

Just writing `tox`, will trigger automatically the functionality tests, code styling tests,
typehint tests and documentation test. These will test the entire Feature-engine ecosystem
and not just your new scripts, so it will be more time consuming.

If the tests pass, the code is in optimal condition :)

**A few things to note:**

Tox runs our tests in Python versions 3.6, 3.7, 3.8 and 3.9. However, it will only be able to
run the tests in the version you have installed locally. All others will fail. This is OK.
As long as the tests in the Python version you have installed pass, you are good to go.

Tox may modify some local files that are not relevant to your feature. Please **DO NOT** add
those files to your PR.

If you want to know more about tox check this `link <https://tox.readthedocs.io>`_. If
you want to know why we prefer tox, this
`article <https://christophergs.com/python/2020/04/12/python-tox-why-use-it-and-tutorial/>`_
will tell you everything ;)


Review Process
--------------

Once your contribution contains the new code, the tests and the documentation, you can
request a review by mentioning that in a comment in the Pull Request. Likely, there will
be some back and forth until the final submission. We will work together to get the code
in the final shape.

If you feel you would like to clarify something before the first draft is done, or if
you can't get some tests to pass, do not hesitate to mention that in a comment, and we
will try to help.

**We aim to review PRs within a week. If for some reason we can't, we will let you know
through the PR as well.**

Once the submission is reviewed and provided the continuous integration tests have
passed and the code is up to date with Feature-engine's main branch, we will be ready
to "Squash and Merge" your contribution into the ``main`` branch of Feature-engine.
"Squash and Merge" combines all of your commits into a single commit which helps keep
the history of the repository clean and tidy.

Once your contribution has been merged into main, you will be listed as a
Feature-engine contributor :)


Merge Pull Requests
-------------------

Only Core contributors have write access to the repository, can review and merge
pull requests. Some preferences for commit messages when merging in pull requests:

- Make sure to use the “Squash and Merge” option in order to create a Git history that is understandable.
- Keep the title of the commit short and descriptive; be sure it links all related issues.


Releases
--------

After a few features have been added to the main branch by yourself and other
contributors, we will merge main into a release branch, e.g. 1.2.X, to release a new
version of Feature-engine to PyPI and conda-forge.


Keep your Fork up to Date
-------------------------

When you're collaborating using forks, it's important to update your fork to capture
changes that have been made by other collaborators.

If your feature takes a few weeks or months to develop, it may happen that new code
changes are made to Feature_engine's main branch by other contributors. Some of the
files that are changed maybe the same files you are working on. Thus, it is really
important that you pull and rebase the upstream main branch into your feature branch.
To keep your branches up to date:

1. Check out your local main branch::

    $ git checkout main

If your feature branch has uncommitted changes, it will ask you to commit or stage those
first. Refer to the commit guidelines we described above.

2. Pull and rebase the upstream main branch on your local main branch::

    $ git pull --rebase upstream main

Your main should be a copy of the upstream main after this. If was is not, there may appear
some conflicting files. You will need to resolve these conflicts and continue the rebase.

3. Pull the changes to your fork::

    $ git push -f origin main

The previous command will update your fork (remote) so that your fork's main branch is in sync with
Feature-engine's main. Now, you need to rebase main onto your feature branch.

4. Check out your feature branch::

    $ git checkout myfeaturebranch

5. Rebase main onto it::

    $ git rebase main

Again, if conflicts arise, try and resolve them and continue the rebase.

Now you are good to go to continue developing your feature... _datasets:

Datasets
========

The user guide and examples included in Feature-engine's documentation are based on
these 3 datasets:

Titanic dataset
~~~~~~~~~~~~~~~

We use the dataset available in `openML <https://www.openml.org/d/40945>`_ which can be
downloaded from `here <https://www.openml.org/data/get_csv/16826755/phpMYEkMl>`_.

Ames House Prices dataset
~~~~~~~~~~~~~~~~~~~~~~~~~

We use the data set created by Professor Dean De Cock:
* Dean De Cock (2011) Ames, Iowa: Alternative to the Boston Housing
* Data as an End of Semester Regression Project, Journal of Statistics Education, Vol.19, No. 3.

The examples are based on a copy of the dataset available on
`Kaggle <https://www.kaggle.com/c/house-prices-advanced-regression-techniques/data>`_.

The original data and documentation can be found here:

* `Documentation <http://jse.amstat.org/v19n3/decock/DataDocumentation.txt>`_

* `Data <http://jse.amstat.org/v19n3/decock/AmesHousing.xls>`_

Credit Approval dataset
~~~~~~~~~~~~~~~~~~~~~~~

We use the Credit Approval dataset from the UCI Machine Learning Repository:

Dua, D. and Graff, C. (2019). `UCI Machine Learning Repository <http://archive.ics.uci.edu/ml>`_.
Irvine, CA: University of California, School of Information and Computer Science.

To download the dataset visit this
`website <http://archive.ics.uci.edu/ml/machine-learning-databases/credit-screening/>`_
and click on "crx.data" to download the data set.

To prepare the data for the examples:

.. code:: python

    import random
    import pandas as pd
    import numpy as np

    # load data
    data = pd.read_csv('crx.data', header=None)

    # create variable names according to UCI Machine Learning information
    varnames = ['A'+str(s) for s in range(1,17)]
    data.columns = varnames

    # replace ? by np.nan
    data = data.replace('?', np.nan)

    # re-cast some variables to the correct types
    data['A2'] = data['A2'].astype('float')
    data['A14'] = data['A14'].astype('float')

    # encode target to binary
    data['A16'] = data['A16'].map({'+':1, '-':0})

    # save the data
    data.to_csv('creditApprovalUCI.csv', index=False).. -*- mode: rst -*-
.. _quick_start:

Quick Start
===========

If you're new to Feature-engine this guide will get you started. Feature-engine
transformers have the methods `fit()` and `transform()` to learn parameters from the
data and then modify the data. They work just like any Scikit-learn transformer.


Installation
------------

Feature-engine is a Python 3 package and works well with 3.6 or later. Earlier versions
have not been tested. The simplest way to install Feature-engine is from PyPI with pip:

.. code-block:: bash

    $ pip install feature-engine


Note, you can also install it with a _ as follows:

.. code-block:: bash

    $ pip install feature_engine


Note that Feature-engine is an active project and routinely publishes new releases. In
order to upgrade Feature-engine to the latest version, use ``pip`` as follows.

.. code-block:: bash

    $ pip install -U feature-engine

If you’re using Anaconda, you can install the
`Anaconda Feature-engine package <https://anaconda.org/conda-forge/feature_engine>`_:

.. code-block:: bash

    $ conda install -c conda-forge feature_engine

Once installed, you should be able to import Feature-engine without an error, both in
Python and in Jupyter notebooks.


Example Use
-----------
This is an example of how to use Feature-engine's transformers to perform missing data
imputation.

.. code:: python

	import numpy as np
	import pandas as pd
	import matplotlib.pyplot as plt
	from sklearn.model_selection import train_test_split

	from feature_engine.imputation import MeanMedianImputer

	# Load dataset
	data = pd.read_csv('houseprice.csv')

	# Separate into train and test sets
	X_train, X_test, y_train, y_test = train_test_split(
    	    data.drop(['Id', 'SalePrice'], axis=1),
            data['SalePrice'],
            test_size=0.3,
            random_state=0
        )

	# set up the imputer
	median_imputer = MeanMedianImputer(
            imputation_method='median', variables=['LotFrontage', 'MasVnrArea']
            )

	# fit the imputer
	median_imputer.fit(X_train)

	# transform the data
	train_t = median_imputer.transform(X_train)
	test_t = median_imputer.transform(X_test)

	fig = plt.figure()
	ax = fig.add_subplot(111)
	X_train['LotFrontage'].plot(kind='kde', ax=ax)
	train_t['LotFrontage'].plot(kind='kde', ax=ax, color='red')
	lines, labels = ax.get_legend_handles_labels()
	ax.legend(lines, labels, loc='best')

.. image:: ../images/medianimputation.png


Feature-engine with the Scikit-learn's pipeline
-----------------------------------------------

Feature-engine's transformers can be assembled within a Scikit-learn pipeline. This
way, we can store our entire feature engineering pipeline in one single object or
pickle (.pkl). Here is an example of how to do it:

.. code:: python

    from math import sqrt
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    from sklearn.linear_model import Lasso
    from sklearn.metrics import mean_squared_error
    from sklearn.model_selection import train_test_split
    from sklearn.pipeline import Pipeline as pipe
    from sklearn.preprocessing import MinMaxScaler
    
    from feature_engine.encoding import RareLabelEncoder, MeanEncoder
    from feature_engine.discretisation import DecisionTreeDiscretiser
    from feature_engine.imputation import (
        AddMissingIndicator,
        MeanMedianImputer,
        CategoricalImputer,
    )

    # load dataset
    data = pd.read_csv('houseprice.csv')

    # drop some variables
    data.drop(
        labels=['YearBuilt', 'YearRemodAdd', 'GarageYrBlt', 'Id'],
        axis=1,
        inplace=True
    )

    # make a list of categorical variables
    categorical = [var for var in data.columns if data[var].dtype == 'O']

    # make a list of numerical variables
    numerical = [var for var in data.columns if data[var].dtype != 'O']

    # make a list of discrete variables
    discrete = [ var for var in numerical if len(data[var].unique()) < 20]

    # categorical encoders work only with object type variables
    # to treat numerical variables as categorical, we need to re-cast them
    data[discrete]= data[discrete].astype('O')

    # continuous variables
    numerical = [
        var for var in numerical if var not in discrete
        and var not in ['Id', 'SalePrice']
        ]

    # separate into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
                                            data.drop(labels=['SalePrice'], axis=1),
                                            data.SalePrice,
                                            test_size=0.1,
                                            random_state=0
                                            )

    # set up the pipeline
    price_pipe = pipe([
        # add a binary variable to indicate missing information for the 2 variables below
        ('continuous_var_imputer', AddMissingIndicator(variables=['LotFrontage'])),

        # replace NA by the median in the 2 variables below, they are numerical
        ('continuous_var_median_imputer', MeanMedianImputer(
            imputation_method='median', variables=['LotFrontage', 'MasVnrArea']
        )),

        # replace NA by adding the label "Missing" in categorical variables
        ('categorical_imputer', CategoricalImputer(variables=categorical)),

        # disretise continuous variables using trees
        ('numerical_tree_discretiser', DecisionTreeDiscretiser(
            cv=3,
            scoring='neg_mean_squared_error',
            variables=numerical,
            regression=True)),

        # remove rare labels in categorical and discrete variables
        ('rare_label_encoder', RareLabelEncoder(
            tol=0.03, n_categories=1, variables=categorical+discrete
        )),

        # encode categorical and discrete variables using the target mean
        ('categorical_encoder', MeanEncoder(variables=categorical+discrete)),

        # scale features
        ('scaler', MinMaxScaler()),

        # Lasso
        ('lasso', Lasso(random_state=2909, alpha=0.005))

    ])

    # train feature engineering transformers and Lasso
    price_pipe.fit(X_train, np.log(y_train))

    # predict
    pred_train = price_pipe.predict(X_train)
    pred_test = price_pipe.predict(X_test)

    # Evaluate
    print('Lasso Linear Model train mse: {}'.format(
        mean_squared_error(y_train, np.exp(pred_train))))
    print('Lasso Linear Model train rmse: {}'.format(
        sqrt(mean_squared_error(y_train, np.exp(pred_train)))))
    print()
    print('Lasso Linear Model test mse: {}'.format(
        mean_squared_error(y_test, np.exp(pred_test))))
    print('Lasso Linear Model test rmse: {}'.format(
        sqrt(mean_squared_error(y_test, np.exp(pred_test)))))


.. code:: python

    Lasso Linear Model train mse: 949189263.8948538
    Lasso Linear Model train rmse: 30808.9153313591

    Lasso Linear Model test mse: 1344649485.0641894
    Lasso Linear Model train rmse: 36669.46256852136

.. code:: python

    plt.scatter(y_test, np.exp(pred_test))
    plt.xlabel('True Price')
    plt.ylabel('Predicted Price')
    plt.show()

.. image:: ../images/pipelineprediction.png


More examples
~~~~~~~~~~~~~

More examples can be found in:

- :ref:`User Guide <user_guide>`
- :ref:`Learning Resources <learning_resources>`
- `Jupyter notebooks <https://nbviewer.jupyter.org/github/feature-engine/feature-engine-examples/tree/main/>`_

.. toctree::
   :maxdepth: 1
   :hidden:

   datasets