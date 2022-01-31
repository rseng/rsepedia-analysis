---
title: 'AutoFunc: A Python package for automating and verifying functional modeling'
tags:
  - Python
  - engineering design
  - functional modeling
  - design repository
  - data mining
authors:
  - name: Alex Mikes
    affiliation: 1
  - name: Katherine Edmonds
    affiliation: 1
  - name: Robert B. Stone
    affiliation: 1
  - name: Bryony DuPont
    affiliation: 1
affiliations:
 - name: Design Engineering Lab, Oregon State University
   index: 1
date: 2 March 2020
bibliography: paper.bib

---

# Statement of Need

Engineering design is a multi-step process that uses various iterative tools to help improve products. Each component 
in a product performs a corresponding set of subfunctions that contribute to the overall functionality 
of the product. Designers often store this product information, including components and subfunction relationships, in
a database known as a design repository. In addition to storing product information, it also helpful to visualize it
in a graphical representation known as a functional model. Functional modeling is a popular tool in the early design
phases that helps designers ensure the product adheres to the customer requirements while maintaining the 
desired functionality. While significant work has been done to help increase consistency in the structure, syntax, 
and formatting of functional models, they are still highly subjective and time-consuming to create [@Stone2000; @Hirtz2002]. 
Because of the time requirements, inconsistencies, and inaccuracies involved with making them, functional models are 
often omitted from the concept generation process, despite their useful contributions to the early stages of 
engineering design [@Kurfman2003]. 

# Summary

``AutoFunc`` is a Python package that automatically generates the functional representations of components based on data from 
design repositories. The functional representations of components can be connected to form a complete functional model. 
``AutoFunc`` also contains methods to validate and optimize the automation algorithm. A designer can use this software to 
input a list of components in their product, and it will automatically generate the functional representations for those 
components based on the most commonly seen functions and flows from previous products in the design repository. 
The package uses common data-mining techniques for finding information and classifying new observations based on 
that data. ``AutoFunc`` also uses the common methods of cross-validation and the F1 score to find the accuracy at 
different values for the threshold variables.

``AutoFunc`` is intended for use by engineering design researchers, students, and professionals. It has been used in 
several engineering design publications and presentations [@mikes2020;@edmonds2020_2]. Further development is required to 
automate a complete functional model, but this software is a significant step in that direction. Automating functional 
modeling will help standardize the format and syntax, decrease the time required to make them, and increase the 
prevalence and accuracy of functional models in engineering design and design repositories. ``AutoFunc`` has been 
archived to Zenodo with the linked DOI [@alexmikes]


# Acknowledgements

Thanks to Kyle Niemeyer for support with open software and open science.

This material is based upon work supported by the National Science Foundation under Grant No. CMMI-1826469. 
Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and 
do not necessarily reflect the views of the National Science Foundation.

# References
# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

## 0.1.1 - 2019-06-11
### Added
- Final Docstrings and tests for project submission

## 0.0.1 - 2019-05-08
### Added
- This CHANGELOG file 

# AutoFunc
Data Mining for Automated Functional Representations

[![Build Status](https://travis-ci.org/AlexMikes/AutoFunc.svg?branch=master)](https://travis-ci.org/AlexMikes/AutoFunc)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3243689.svg)](https://doi.org/10.5281/zenodo.3243689)

``AutoFunc`` is a Python package that automatically generates the functional representations of components based on data from 
design repositories. ``AutoFunc`` also contains methods to validate and optimize the automation algorithm. A designer can use this software to 
input a list of components in their product, and it will automatically generate the functional representations for those 
components based on the most commonly seen functions and flows from previous products in the design repository. 
The package uses common data-mining techniques for finding information and classifying new observations based on 
that data. ``AutoFunc`` also uses the common methods of cross-validation and the F1 score to find the accuracy at 
different values for the threshold variables.

`AutoFunc` was developed for use with the Design Repository housed at Oregon State University. A rudimentary 
web interface can be found here: http://ftest.mime.oregonstate.edu/repo/browse/

## Installation

`autofunc` has been tested on Linux Python 3.6, 3.7, and 3.8

### Pip

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install autofunc.

The package is not yet on PyPI, so it must be downloaded from here as a .zip file: https://github.com/AlexMikes/AutoFunc

Once downloaded as a .zip file, install with:

```bash
pip install /path/to/file/AutoFunc-master.zip
```

### From repository

To install from this repository:

```bash
git clone https://github.com/AlexMikes/AutoFunc.git
cd autofunc
python setup.py install
```

## Dependencies

This package uses Pandas (Python Data Analysis Library). It can be installed with pip using:

```bash
pip install pandas
```
Many of the examples also use Matplotlib for plotting. While not required to use the AutoFunc modules, it is required to run the examples. It can be installed with:

```bash
pip install -U matplotlib
```

## Usage

Every module has NumPy formatted docstrings to explain the inputs, outputs, and usage of each of them.
Rudimentary API documentation can be found here: https://autofunc.readthedocs.io/

Example files are provided in the examples folder. Autofunc will automate the functional representations of components
as  long as the format of the .csv file has the component in column 1 and the function-flow in column 2

More information on the methods used in these files can found in the various research papers that this software supports, especcially IDETC2020-22346
"OPTIMIZING AN ALGORITHM FOR DATA MINING A DESIGN REPOSITORY TO AUTOMATE FUNCTIONAL MODELING". All of the plots for this paper were created in the ```example_optimize_with_comp_ratio.py``` file.

The following lists the examples included, with their expected functionality and outputs:

1. ```example_cross_validation.py``` uses the k-fold cross validation functionality to find the accuracy of a data mining classifer. This example will print the maximum and average accuracies using this verification method.

1. ```example_find_f1_from_file.py``` finds the F1 score of a single product when the component-function-flow combinations for that product are in a separate .csv file. This example will print the Recall, Precision, and F1 score for that testing product.

1. ```example_find_f1_from_id.py``` finds the F1 score of a single product using that product's ID number from the original dataset. Any number of IDs can be used. This example will print the testing ID(s) used, and the recall, precision, and F1 score for those testing IDs.

1. ```example_find_similarity.py``` will create a similarity matrix for the training dataset. This is the percent of similar components between each product in the dataset. The main diagonal of this matrix consists of ones because every product is 100% similar to itself, but the matrix is not symmetric because each product can contain a different number of components. For example, consider a case where Product 1 has 20 components and Product 2 has 40 components. If they have 10 components in common, the similarity between Product 1 and Product 2 is 10/20 = 50%, but the similarity between Product 2 and Product 1 is 10/40 = 25%. The first product of the pair is known as the ”generating” product, which is the product in the column of this matrix. This example will create a Pandas dataframe of the similarity matrix and write this to a .csv file.

1. ```example_get_func_rep.py``` will create a functional representation of the components in the input file using data mining and a classification threshold. This can be used to automate functional modeling by connecting the functions and flows at the interface of components in a product. This example will write a .csv file with the results of component-function-flow and optional frequency.

1. ```example_optimization.py``` incorporates all of the main modules and optimizes the similarity and classification thresholds. This example will display a lot of plots and print optimum values for thresholds.

1. ```example_optimize_with_comp_ratio.py``` begins with ```example_optimization.py``` and also includes the stratification and optimization of a training set. This example will display a lot of plots and print optimum values for thresholds.

1. ```example_try_best_ids.py``` is a subset of ```example_optimize_with_comp_ratio.py``` which only includes the stratified training set and some
relevant plots. This example will display a plot of the F1 scores vs. Classification threshold of the stratified dataset.


This is the ```example_get_func_rep.py``` file:

```python
from autofunc.get_top_results import get_top_results
from autofunc.counter_pandas import counter_pandas
from autofunc.get_func_rep import get_func_rep
import os.path
import pandas as pd


""" Example showing how to automate functional representation """

# Dataset used for data mining
script_dir = os.path.dirname(__file__)
file_to_learn = os.path.join(script_dir, '../assets/consumer_systems.csv')

include_frequencies = True

train_data = pd.read_csv(file_to_learn)
combos_sorted = counter_pandas(train_data)

# Use a threshold to get the top XX% of confidence values
threshold = 0.5
thresh_results = get_top_results(combos_sorted, threshold)

# Use a known product for verification
input_file = os.path.join(script_dir, '../assets/InputExample.csv')

# Get dictionary of functions and flows for each component based on data mining
results, unmatched = get_func_rep(thresh_results, input_file, include_frequencies)


# Optional write to file - uncomment and rename to write file
write_results_from_dict(results, 'test1.csv')
```


Run from within ```examples``` folder:

```bash
python example_get_func_rep.py
```

And it will generate a file ```test1.csv``` with the results of the automated functional representation of the 
 components in the ```input_file``` based on the data from the ```file_to_learn``` in the ```assets``` folder.
 
## Testing
All tests are automated through [Travis CI](https://travis-ci.org/). Visit [this page](https://travis-ci.org/github/AlexMikes/AutoFunc) to view the results.


## Support
Please submit requests for support or problems with software as issues in the repository.

## Contributing

We welcome contributions to the `autofunc` package in the form of [pull requests](https://github.com/AlexMikes/AutoFunc/pulls) and [issues](https://github.com/AlexMikes/AutoFunc/issues) made in the repository.

If you are having any problems using `autofunc`, please open an issue.
If there is some functionality you would like to see added to `autofunc`, you can also open an issue up to discuss that.

If you have a feature that you would like to propose be integrated into `autofunc`, then you should open a pull request.

## License
[MIT](https://choosealicense.com/licenses/mit/)
autofunc
========

.. toctree::
   :maxdepth: 4

   autofunc
.. autofunc documentation master file, created by
   sphinx-quickstart on Sun Nov  1 07:53:07 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to autofunc's documentation!
====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
autofunc package
================

Submodules
----------

autofunc.counter\_pandas module
-------------------------------

.. automodule:: autofunc.counter_pandas
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.counter\_pandas\_with\_counts module
---------------------------------------------

.. automodule:: autofunc.counter_pandas_with_counts
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.df\_to\_list module
----------------------------

.. automodule:: autofunc.df_to_list
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.find\_similarities module
----------------------------------

.. automodule:: autofunc.find_similarities
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.get\_func\_rep module
------------------------------

.. automodule:: autofunc.get_func_rep
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.get\_precision\_recall module
--------------------------------------

.. automodule:: autofunc.get_precision_recall
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.get\_top\_results module
---------------------------------

.. automodule:: autofunc.get_top_results
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.make\_df module
------------------------

.. automodule:: autofunc.make_df
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.split\_learning\_verification module
---------------------------------------------

.. automodule:: autofunc.split_learning_verification
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.write\_results module
------------------------------

.. automodule:: autofunc.write_results
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: autofunc
   :members:
   :undoc-members:
   :show-inheritance:
autofunc
========

.. toctree::
   :maxdepth: 4

   autofunc
autofunc package
================

Submodules
----------

autofunc.counter\_pandas module
-------------------------------

.. automodule:: autofunc.counter_pandas
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.counter\_pandas\_with\_counts module
---------------------------------------------

.. automodule:: autofunc.counter_pandas_with_counts
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.df\_to\_list module
----------------------------

.. automodule:: autofunc.df_to_list
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.find\_similarities module
----------------------------------

.. automodule:: autofunc.find_similarities
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.get\_func\_rep module
------------------------------

.. automodule:: autofunc.get_func_rep
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.get\_precision\_recall module
--------------------------------------

.. automodule:: autofunc.get_precision_recall
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.get\_top\_results module
---------------------------------

.. automodule:: autofunc.get_top_results
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.make\_df module
------------------------

.. automodule:: autofunc.make_df
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.split\_learning\_verification module
---------------------------------------------

.. automodule:: autofunc.split_learning_verification
   :members:
   :undoc-members:
   :show-inheritance:

autofunc.write\_results module
------------------------------

.. automodule:: autofunc.write_results
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: autofunc
   :members:
   :undoc-members:
   :show-inheritance:
