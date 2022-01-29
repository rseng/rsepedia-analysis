# Change log

## 0.4.0 (03/07/2021)

* added pre-commit
* updated to laspy 2.0.0 (with optional laszip)
* updated bumpversion to double quotes (Python black)
* updated docs build
* tests and release to pypi and github with github actions
* added github issue templates, pull request template and updated contributing docs
* using pydocstyle and black
* added lidar to canopy module
* using keyword-only arguments

## 0.3.1 (11/03/2020)

* drop `num_blocks` argument of `compute_image_descriptor` and `compute_image_descriptor_from_filepath`

## 0.3.0 (02/03/2020)

* set default post-classification refinement parameter `refine_beta` to 50 (instead of 100)
* keyword arguments to `PixelFeaturesBuilder` and `PixelResponseBuilder` can be explicitly provided to the initialization of `ClassifierTrainer`, and are documented there
* raise a `ValueError` when a provided response is not a binary tree/non-tree image

## 0.2.0 (11/12/2019)

* correction (typo) `keep_emtpy_tiles` -> `keep_empty_tiles` in `split_into_tiles`

## 0.1.0 (14/11/2019)

* initial release
[![PyPI version fury.io](https://badge.fury.io/py/detectree.svg)](https://pypi.python.org/pypi/detectree/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/detectree.svg)](https://anaconda.org/conda-forge/detectree)
[![Documentation Status](https://readthedocs.org/projects/detectree/badge/?version=latest)](https://detectree.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://github.com/martibosch/detectree/workflows/tests/badge.svg?branch=main)](https://github.com/martibosch/detectree/actions?query=workflow%3A%22tests%22)
[![codecov](https://codecov.io/gh/martibosch/detectree/branch/main/graph/badge.svg?token=ZTZK2LFR6T)](https://codecov.io/gh/martibosch/detectree)
[![GitHub license](https://img.shields.io/github/license/martibosch/detectree.svg)](https://github.com/martibosch/detectree/blob/master/LICENSE)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02172/status.svg)](https://doi.org/10.21105/joss.02172)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3908338.svg)](https://doi.org/10.5281/zenodo.3908338)

# DetecTree

## Overview

DetecTree is a Pythonic library to classify tree/non-tree pixels from aerial imagery, following the methods of Yang et al. [1]. The target audience is researchers and practitioners in GIS that are interested in two-dimensional aspects of trees, such as their proportional abundance and spatial distribution throughout a region of study. These measurements can be used to assess important aspects of urban planning such as the provision of urban ecosystem services. The approach is of special relevance when LIDAR data is not available or it is too costly in monetary or computational terms.

```python
import detectree as dtr
import matplotlib.pyplot as plt
import rasterio as rio
from rasterio import plot

# select the training tiles from the tiled aerial imagery dataset
ts = dtr.TrainingSelector(img_dir='data/tiles')
split_df = ts.train_test_split(method='cluster-I')

# train a tree/non-tree pixel classfier
clf = dtr.ClassifierTrainer().train_classifier(
    split_df=split_df, response_img_dir='data/response_tiles')

# use the trained classifier to predict the tree/non-tree pixels
test_filepath = split_df[~split_df['train'].sample(1).iloc[0]['img_filepath']
y_pred = dtr.Classifier().classify_img(test_filepath, clf)

# side-by-side plot of the tile and the predicted tree/non-tree pixels
figwidth, figheight = plt.rcParams['figure.figsize']
fig, axes = plt.subplots(1, 2, figsize=(2 * figwidth, figheight))

with rio.open(img_filepath) as src:
    plot.show(src.read(), ax=axes[0])
axes[1].imshow(y_pred)
```

![Example](figures/example.png)

A full example application of DetecTree to predict a tree canopy map for the Aussersihl district in Zurich [is available as a Jupyter notebook](https://github.com/martibosch/detectree-example/blob/master/notebooks/aussersihl-canopy.ipynb). See also [the API reference documentation](https://detectree.readthedocs.io/en/latest/?badge=latest) and the [example repository](https://github.com/martibosch/detectree-example) for more information on the background and some example notebooks.

## Citation

Bosch M. 2020. “DetecTree: Tree detection from aerial imagery in Python”. *Journal of Open Source Software, 5(50), 2172.* [doi.org/10.21105/joss.02172](https://doi.org/10.21105/joss.02172)

Note that DetecTree is based on the methods of Yang et al. [1], therefore it seems fair to reference their work too. An example citation in an academic paper might read as follows:

> The classification of tree pixels has been performed with the Python library DetecTree (Bosch, 2020), which is based on the approach of Yang et al. (2009).

## Installation

### With conda

The easiest way to install `detectree` is with conda as in:

``` bash
conda install -c conda-forge detectree
```

### With pip

You can install `detectree` with pip as in:

``` bash
pip install detectree
```

If you want to be able to read compressed LAZ files, you will need [the Python bindings for `laszip`](https://github.com/tmontaigu/laszip-python). Note that the latter require [`laszip`], which can be installed using conda (which is automatically handled when installing `detectree` with conda as shown above) or downloaded from [laszip.org](https://laszip.org/). Then, detectree and the Python bindings for `laszip` can be installed with pip as in:

``` bash
pip install detectree[laszip]
```

### Development install

To install a development version of detectree, you can first use conda to create an environment with all the dependencies - with the [`environment-dev.yml` file](https://github.com/martibosch/detectree/blob/main/environment-dev.yml) - and activate it as in:

``` bash
conda env create -f environment-dev.yml
conda activate detectree-dev
```

and then clone the repository and use pip to install it in development mode

```bash
git clone git@github.com:martibosch/detectree.git
cd detectree/
pip install -e .
```

This will also install the dependencies required for running tests, linting the code and building the documentation. Additionally, you can activate [pre-commit](https://pre-commit.com/) so that the latter are run as pre-commit hooks as in:

```bash
pre-commit install
```

## See also

* [lausanne-tree-canopy](https://github.com/martibosch/lausanne-tree-canopy): example computational workflow to get the tree canopy of Lausanne with DetecTree
* [A video of a talk about DetecTree](https://www.youtube.com/watch?v=USwF2KyxVjY) in the [Applied Machine Learning Days of EPFL (2020)](https://appliedmldays.org/) and [its respective slides](https://martibosch.github.io/detectree-amld-2020)

## Acknowledgments

* With the support of the École Polytechnique Fédérale de Lausanne (EPFL)

## References

1. Yang, L., Wu, X., Praun, E., & Ma, X. (2009). Tree detection from aerial imagery. In Proceedings of the 17th ACM SIGSPATIAL International Conference on Advances in Geographic Information Systems (pp. 131-137). ACM.
# Contributing

Contributions are always greatly appreciated and credit will always be given.

## Types of contributions

### Report bugs

Report bugs at https://github.com/martibosch/detectree/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

### Fix bugs

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help wanted" is open to whoever wants to implement it.

### Implement features

Look through the GitHub issues for features. Anything tagged with "enhancement" and "help wanted" is open to whoever wants to implement it.

## Pull request guidelines

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put your new functionality into a function with a docstring, and add the feature to the list in README.md.
3. The pull request should work for Python 3.6, 3.7, 3.8 and 3.9. Check https://travis-ci.org/martibosch/detectree/pull_requests and make sure that the tests pass for all supported Python versions.
4. Adhere to the following project standards:
    * `black` code style with max line length of 79
    * `isort` sorted imports
    * `numpy` style docstrings
---
title: 'DetecTree: Tree detection from aerial imagery in Python'
tags:
  - Python
  - tree detection
  - image segmentation  
  - remote sensing images
  - GIS
authors:
 - name: Martí Bosch
   orcid: 0000-0001-8735-9144
   affiliation: 1
affiliations:
 - name: Urban and Regional Planning Community, École Polytechnique Fédérale de Lausanne, Switzerland
   index: 1
date: 4 March 2020
bibliography: paper.bib
---

# Summary

Urban tree canopy datasets are important to a wide variety of social and environmental studies in cities. Surveys to collect such information manually are costly and hard to maintain, which has motivated a growing interest in automated approaches to tree detection in recent years. To that end, LIDAR point clouds are a very appropriate data source, especially given its capability to represent spatial features in three dimensions. However, collecting LIDAR data requires expensive equipment and raw datasets are rarely made openly available. On the other hand, aerial imagery is another major source of data that is able to capture a wide range of objects on Earth, including trees. Although aerial imagery depicts the spatial features in only two dimensions, its main advantage with respect to LIDAR is its greater availability.

The aim of DetecTree is therefore to provide an open source library that performs a binary classification of tree/non-tree pixels from aerial imagery. The target audience is researchers and practitioners in GIS that are interested in two-dimensional aspects of trees, such as their proportional abundance and spatial distribution throughout a region of study. These measurements can be used to assess important aspects of urban planning such as the provision of urban ecosystem services. The approach is of special relevance when LIDAR data is not available or it is too costly in monetary or computational terms.

# Methodology

DetecTree is based on the supervised learning approach of @yang2009tree, which requires an RGB aerial imagery dataset as the only input, and consists of the following steps:

* **Step 0**: split of the dataset into image tiles. Since aerial imagery datasets often already come as a mosaic of image tiles, this step might not be necessary. In any case, DetecTree provides a `split_into_tiles` function that can be used to divide a large image into a mosaic of tiles of a specified dimension.
* **Step 1**: selection of the tiles to be used for training a classifier. As a supervised learning task, the ground-truth maps must be provided for some subset of the dataset. Since this part is likely to involve manual work, it is crucial that the training set has as few tiles as possible. At the same time, to enhance the classifier's ability to detect trees in the diverse scenes of the dataset, the training set should contain as many of the diverse geographic features as possible. Thus, in order to optimize the representativity of the training set, the training tiles are selected according to their GIST descriptor [@oliva2001modeling], *i.e.*, a vector describing the key semantics of the tile's scene. More precisely, *k*-means clustering is applied to the GIST descriptors of all the tiles, with the number of clusters *k* set to the number of tiles of the training set (by default, one percent of the tiles is used). Then, for each cluster, the tile whose GIST descriptor is closest to the cluster's centroid is added to the training set. In DetecTree, this is done by the `train_test_split` method of the `TrainingSelector` class.
* **Step 2**: provision of the ground truth tree/non-tree masks for the training tiles. For each tile of the training set, the ground-truth tree/non-tree masks must be provided to get the pixel-level responses that will be used to train the classifier. To that end, an image editing software such as GIMP or Adobe Photoshop might be used. Alternatively, if LIDAR data for the training tiles is available, it might also be exploited to create the ground truth masks.
* **Step 3**: train a binary pixel-level classifier. For each pixel of the training tiles, a vector of 27 features is computed, where 6, 18 and 3 features capture characteristics of color, texture and entropy respectively. A binary AdaBoost classifier [@freund1995desicion] is then trained by mapping the feature vector of each pixel to its class in the ground truth masks (*i.e.*, tree or non-tree).
* **Step 4**: tree detection in the testing tiles. Given a trained classifier, the `classify_img` and `classify_imgs` methods of the `Classifier` class can respectively be used to classify the tree pixels of a single image tile or of multiple image tiles at scale. For each image tile, the pixel-level classification is refined by means of a graph cuts algorithm [@boykov2004experimental] to avoid sparse pixels classified as trees by enforcing consistency between adjacent tree pixels. An example of an image tile, its pre-refinement pixel-level classification and the final refined result is displayed below:

![Example of an image tile (left), its pre-refinement pixel-level classification (center) and the final refined result (right).](figure.png)

Similar methods of tree classification from aerial imagery include the work of @jain2019efficient, who follow the train/test split method based on GIST descriptors as proposed by @yang2009tree but rely on the Mask R-CNN framework [@he2017mask] instead of the AdaBoost classifier. Another approach by @tianyang2018single employs a cascade neural network over texture and color features which detects single trees in a variety of forest images. Nonetheless, since the former approaches ultimately aim at single tree detection, the accuracy evaluation metrics that they provide are hard to compare with the pixel-level classification accuracy of DetecTree. The experiments performed by @yang2009tree in New York achieve a pixel classification accuracy of 91.7%, whereas the example applications of DetecTree in Zurich and Lausanne achieve accuracies of 85.98% and 91.75% respectively.

The code of DetecTree is organized following an object-oriented approach, and relies on NumPy [@van2011numpy] to represent most data structures and perform operations upon them in a vectorized manner. The Scikit-learn library [@pedregosa2011scikit] is used to implement the AdaBoost pixel-level classifier as well as to perform the *k*-means clustering to select the training tiles. The computation of pixel-level features and GIST descriptors makes use of various features provided by the Scikit-image [@van2014scikit] and SciPy [@virtanen2020scipy] libraries. On the other hand, the classification refinement employs the graph cuts algorithm implementation provided by the library [PyMaxFlow](https://github.com/pmneila/PyMaxflow). Finally, when possible, DetecTree uses the Dask library [@rocklin2015dask] to perform various computations in parallel.

The features of DetecTree are implemented in a manner that enhances the flexibility of the library so that the user can integrate it into complex computational workflows, and also provide custom arguments for the technical aspects. Furthermore, the functionalities of DetecTree can be used through its Python API as well as through its command-line interface (CLI), which is implemented by means of the Click Python package.


# Availability

The source code of DetecTree is fully available at [a GitHub repository](https://github.com/martibosch/detectree). A dedicated Python package has been created and is hosted at the [Python Package Index (PyPI)](https://pypi.org/project/detectree/). The documentation site is hosted at [Read the Docs](https://detectree.readthedocs.io/), and an example repository with Jupyter notebooks of an example application to an openly-available orthophoto of Zurich is provided at a [dedicated GitHub repository](https://github.com/martibosch/detectree-example), which can be executed interactively online by means of the Binder web service [@jupyter2018binder]. An additional example use of DetecTree can be found at a [dedicated GitHub repository](https://github.com/martibosch/lausanne-tree-canopy) with the materials to obtain a tree canopy map for the urban agglomeration of Lausanne from the SWISSIMAGE 2016 orthophoto [@swisstopo2019swissimage].

Unit tests are run within the [Travis CI](https://travis-ci.org/martibosch/detectree) platform every time that new commits are pushed to the GitHub repository. Additionally, test coverage [is reported on Coveralls](https://coveralls.io/github/martibosch/detectree?branch=master).


# Acknowledgments

This research has been supported by the École Polytechnique Fédérale de Lausanne.


# References
Please don't file issues directly, use one of our available templates:

Bugs: https://github.com/martibosch/detectree/issues/new?template=bug.md
Feature Requests: https://github.com/martibosch/detectree/issues/new?template=feature.md

Thank you.
**Read these instructions carefully**

Before you proceed, review the contributing guidelines in the "Pull request guidelines" of the `CONTRIBUTING.md` file.

In this pull request, please include:

* a reference to related issue(s)
* a description of the changes proposed in the pull request
* an example code snippet illustrating usage of the new functionality
* detectree version:
* Python version:
* Operating system:

### Description

Describe what you were trying to get done.
Tell us what happened, what went wrong, and what you expected to happen.

### What I Did

```
Paste the command(s) you ran and the output.
If there was a crash, please include the traceback here.
```
### Description

Describe your feature request
