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
