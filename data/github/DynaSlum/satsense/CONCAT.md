# Changelog

## version 0.9

- Added a large amount of documentation
  - Available at https://satsense.readthedocs.io/en/latest/
  - Includes:
    - Installation instructions
    - Example notebook for feature extraction
    - API Documentation and docstrings

- Bug fixes:
  - Histogram of Gradients
    - fixed the absolute sine difference calculation
  - Fixed the padding around the data when splitting the generator.
  - Fixed the generation of windows

- Development:
  - Added automated versioning
  - Increased code maintainability

## version 0.8
- Initial release
- Features included:
  - Histogram of Gradients
  - Pantex
  - NDVI
    - also available:
    - RgNDVI (Red-green based)
    - RbNDVI (Red-blue based)
    - NDSI (Snow Cover Index)
    - NDWI (Water Cover Index)
    - WVSI (Soil Cover Index)
  - Lacunarity
  - SIFT
  - TextonContributions are very welcome. Please make sure there is a github issue
associated with with every pull request. Creating an issue is also a good
way to propose/discuss new features or get help with using satsense.

# Installation for development

Please follow the installation instructions on
[readthedocs](https://satsense.readthedocs.io/en/latest/installation.html)
to get started.

# Testing

Please add unit tests for the code you are writing (e.g. when fixing a bug, implement
a test that demonstrates the bug is fixed). You can run the unit tests locally
with the command

```python
python setup.py test
```

# Coding style

Please make sure your code is formatted according to
[PEP8](https://www.python.org/dev/peps/pep-0008/) and docstrings are written
according to [PEP257](https://www.python.org/dev/peps/pep-0257/). Publicly visible
functions should have
[numpy style docstrings](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html).

Please autoformat your code with the following commands before making a pull request:

```bash
isort satsense/my_file.py
yapf -i satsense/my_file.py
```

Please use prospector to check that your code meets our standards:

```bash
prospector satsense/my_file.py
```

# Pull requests

Please create a pull request early, to keep other developers informed of what you're doing.
Limit the amount of work in a pull request to fixing a single bug or adding a single new feature.
Make sure the unit tests on Travis pass and review the comments by Codacy (click the Travis/Codacy
buttons below your pull request). Note that Codacy occasionally reports false positives, ask if in
doubt.

# Documentation

All public functions should have
[numpy style docstrings](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html).
You can build the documentation locally by running

```bash
python setup.py build_sphinx
```

Use

```bash
python setup.py build_sphinx -Ea
```

to build everying from scratch. Please check that there are no warnings.

## Converting Notebooks for documentation

If you update the notebooks please update their counterparts in the doc folder by using `jupyter nbconvert`

From the root of the project:
```bash
jupyter nbconvert --to rst notebooks/**/*.ipynb --output-dir=doc/notebooks/
```

# Creating a release

Make sure to update the version number and release date in CITATION.cff.
---
name: User Story
about: Create a scrum user story
title: ''
labels: ''
assignees: ''

---

**As a** <!-- user -->
**I want to** <!-- what -->
**because** <!-- why -->
# Notebooks accompanying the satsense Python package

## FeatureExtraction
Notebooks demonstrating how to extract features from satellite images using
satsense.

## Performance
Notebooks demonstrating how to use the performance metrics with satsense and
the utility functions to convert ground truth data to and from mask files and
shapefiles.

## Classification

Notebooks demonstrating classification using features calculated with Satsense.
Satsense
========

|Build Status| |Codacy Badge| |Maintainability| |Test Coverage|
|Documentation Status| |DOI|

Satsense is an open source Python library for patch based land-use and
land-cover classification, initially developed for a project on deprived
neighborhood detection. However, many of the algorithms made available
through Satsense can be applied in other domains, such as ecology and
climate science.

Satsense is based on readily available open source libraries, such as
opencv for machine learning and the rasterio/gdal and netcdf libraries
for data access. It has a modular design that makes it easy to add your
own hand-crafted feature or use deep learning instead.

Detection of deprived neighborhoods is a land-use classification problem
that is traditionally solved using hand crafted features like HoG,
Lacunarity, NDXI, Pantex, Texton, and SIFT, computed from very high
resolution satellite images. One of the goals of Satsense is to
facilitate assessing the performance of these features on practical
applications. To achieve this Satsense provides an easy to use open
source reference implementation for these and other features, as well as
facilities to distribute feature computation over multiple cpu’s. In the
future the library will also provide easy access to metrics for
assessing algorithm performance.

-  satsense - library for analysing satellite images, performance
   evaluation, etc.
-  notebooks - IPython notebooks for illustrating and testing the usage
   of Satsense

We are using python 3.6/3.7 and jupyter notebook for our code.

Documentation
-------------
Can be found on `readthedocs <https://satsense.readthedocs.io>`__.

Installation
------------

Please see the `installation guide on readthedocs <https://satsense.readthedocs.io/en/latest/installation.html#installation>`__.

Contributing
------------

Contributions are very welcome! Please see
`CONTRIBUTING.md <https://github.com/DynaSlum/satsense/blob/master/CONTRIBUTING.md>`__
for our contribution guidelines.

Citing Satsense
---------------

If you use Satsense for scientific research, please cite it. You can
download citation files from
`research-software.nl <https://www.research-software.nl/software/satsense>`__.

References
----------

The collection of algorithms made available trough this package is
inspired by

    J. Graesser, A. Cheriyadat, R. R. Vatsavai, V. Chandola,
    J. Long and E. Bright, "Image Based Characterization of Formal and
    Informal Neighborhoods in an Urban Landscape", in IEEE Journal of
    Selected Topics in Applied Earth Observations and Remote Sensing,
    vol. 5, no. 4, pp. 1164-1176, Aug. 2012. doi:
    10.1109/JSTARS.2012.2190383

Jordan Graesser himself also maintains `a
library <https://github.com/jgrss/spfeas>`__ with many of these
algorithms.

Test Data
~~~~~~~~~

The test data has been extracted from the Copernicus Sentinel data 2018.

.. |Build Status| image:: https://travis-ci.com/DynaSlum/satsense.svg?branch=master
   :target: https://travis-ci.com/DynaSlum/satsense
.. |Codacy Badge| image:: https://api.codacy.com/project/badge/Grade/458c8543cd304b8387b7b114218dc57c
   :target: https://www.codacy.com/app/DynaSlum/satsense?utm_source=github.com&utm_medium=referral&utm_content=DynaSlum/satsense&utm_campaign=Badge_Grade
.. |Maintainability| image:: https://api.codeclimate.com/v1/badges/ed3655f6056f89f5e107/maintainability
   :target: https://codeclimate.com/github/DynaSlum/satsense/maintainability
.. |Test Coverage| image:: https://api.codeclimate.com/v1/badges/ed3655f6056f89f5e107/test_coverage
   :target: https://codeclimate.com/github/DynaSlum/satsense/test_coverage
.. |Documentation Status| image:: https://readthedocs.org/projects/satsense/badge/?version=latest
   :target: https://satsense.readthedocs.io/en/latest/?badge=latest
.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1463015.svg
   :target: https://doi.org/10.5281/zenodo.1463015
Installation
============

Installing the dependencies
---------------------------
Satsense has a few dependencies that cannot be installed from PyPI:

- the dependencies of the `GDAL <https://pypi.org/project/GDAL/>`_ Python package
- the dependencies of the `netCDF4 <http://unidata.github.io/netcdf4-python/>`_ Python package

**Ubuntu Linux 18.04 and later**

To install the above mentioned dependencies, run

.. code-block:: bash

   sudo apt install libgdal-dev libnetcdf-dev

this probably also works for other Debian-based Linux distributions.

**RPM-based Linux distributions**

To install the above mentioned dependencies, run

.. code-block:: bash

   sudo yum install gdal-devel netcdf-devel


**Conda**

Assuming you have `conda <https://conda.io>`_ installed and have downloaded
the satsense
`environment.yml <https://github.com/DynaSlum/satsense/blob/master/environment.yml>`_
file to the current working directory, you can install
all dependencies by running:

.. code-block:: bash

   conda env create --file environment.yml --name satsense

or you can install just the minimal dependencies by running

.. code-block:: bash

   conda create --name satsense libgdal libnetcdf nb_conda

Make sure to activate the environment after installation:

.. code-block:: bash

   conda activate satsense


Installing Satsense from PyPI
-----------------------------

If you did not use conda to install the dependencies, you may still
want to create and activate a virtual environment for satsense, e.g. using
`venv <https://docs.python.org/3/library/venv.html>`_

.. code-block:: bash

   python3 -m venv ~/venv/satsense
   source ~/venv/satsense/bin/activate

Next, install satsense by running

.. code-block:: bash

   pip install satsense

If you are planning on using the :ref:`notebooks`, you can
install the required extra dependencies with

.. code-block:: bash

   pip install satsense[notebooks]

Installing Satsense from source for development
-----------------------------------------------

Clone the `satsense repository <https://github.com/DynaSlum/satsense>`_,
install the dependencies as described above, go to the directory where
you have checked out satsense and run

.. code-block:: bash

   pip install -e .[dev]

or

.. code-block:: bash

   pip install -e .[dev,notebooks]

if you would also like to use the :ref:`notebooks`.

Please read our
`contribution guidelines <https://github.com/DynaSlum/satsense/blob/master/CONTRIBUTING.md>`_
before starting development.

Known installation issues
-------------------------
If you are experiencing 'NetCDF: HDF errors' after installation with pip,
this may be resolved by using the following command to install

.. code-block:: bash

   pip install --no-binary netcdf4 satsense

see `this rasterio issue <https://github.com/rasterio/rasterio-wheels/issues/12>`_
for more information.
Welcome to satsense's documentation!
====================================

Satsense is a Python library for land use/cover classification using satellite imagery.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   notebooks/index
   api/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Feature Extraction Example
--------------------------

In this example we will extract the Histogram of Gradients (HoG),
Normalized Difference Vegetation Index (NDVI) and the Pantex features
from a test satelite image.

-  The HoG feature captures the distribution of structure orientations.
-  The NDVI feature captures the level of vegetation.
-  The Pantex feature captures the level of built-up structures.

The image will be split into blocks, in this example 20 by 20 pixels,
and each feature is calculated for this block using a certain amount of
context information called a window. A feature can be calculated on
multiple windows to allow for context at different scales.

In this example
~~~~~~~~~~~~~~~

-  First we will define the Features we would like to extract and with
   which window shapes.
-  We will then load the image using the ``Image`` class.
-  Then we will split the image into blocks using the ``FullGenerator``
   Class.
-  Then we will extract the features using the ``extract_features``
   function.

Live iPython Notebook
^^^^^^^^^^^^^^^^^^^^^

If you are reading this example on readthedocs.io a notebook of this
example is available `in the
repository <https://github.com/DynaSlum/satsense/blob/master/notebooks/FeatureExtraction/feature_extraction.ipynb>`__

.. code:: ipython3

    # General imports
    import numpy as np

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    %matplotlib inline

    # Satsense imports
    from satsense import Image
    from satsense.generators import FullGenerator
    from satsense.extract import extract_features
    from satsense.features import NirNDVI, HistogramOfGradients, Pantex

Define the features to calculate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First we define a list of windows for each of the features to use.

Hog and Pantex will be calculated on 2 windows of 25x25 pixels and 23x27
pixels. NDVI will be calculated on one window with 37x37 pixels.

These window shapes are chose arbitrarily to show the capabilities of
satsense, for your own feature extraction you should think and
experiment with these window shapes to give you the best results.

N.B. The NDVI feature here is called NirNDVI because that implementation
uses the near-infrared band of the image, there are several other
implementations of NDVI available in satsense, see `the
documentation <https://satsense.readthedocs.io/en/latest/api/satsense.features.html>`__

.. code:: ipython3

    # Multiple windows
    two_windows = [(25, 25), (23, 37)]
    # Single window
    one_window = [(37, 37),]
    features = [
        HistogramOfGradients(two_windows),
        NirNDVI(one_window),
        Pantex(two_windows),
    ]

Load the image
~~~~~~~~~~~~~~

Here we load the image and normalize it to values between 0 and 1.
Normalization by default is performed per band using the 2nd and 98th
percentiles.

The image class can provide the individual bands, or a number of useful
derivatives such as the RGB image or Grayscale, we call these base
images. More advanced base images are also available, for instance Canny
Edge

.. code:: ipython3

    image = Image('../../test/data/source/section_2_sentinel.tif',
                    'quickbird')
    image.precompute_normalization()

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 8), sharey=True)

    ax1.axis('off')
    ax1.imshow(image['rgb'])
    ax1.set_title('RGB image')

    ax2.axis('off')
    ax2.imshow(image['grayscale'], cmap="gray")
    ax2.set_title('Grayscale image')

    ax3.axis('off')
    ax3.imshow(image['canny_edge'], cmap="gray")
    ax3.set_title('Canny Edge Image')

    plt.show()



.. image:: feature_extraction_files/feature_extraction_5_0.png


Generator
~~~~~~~~~

Next we create a FullGenerator which creates patches of the image in
steps of 20x20 pixels.

In this cell we also show the images, therefore we load the rgb base
image into the generator. This is only needed here so we can show the
blocks using matplotlib. In the next section we will be using the
``extract_features`` function to extract features, which will be loading
the correct base images for you based on the features that will be
calculated.

The patch sizes are determined by the list of window shapes that you
supply the ``load_image`` function. This is normally also provided by
the ``extract_features`` function.

.. code:: ipython3

    generator = FullGenerator(image, (20, 20))
    print("The generator is {} by {}".format(*generator.shape), " blocks")

    # Create a gridspec to show the images
    gs = gridspec.GridSpec(*generator.shape)
    gs.update(wspace=0.05, hspace=0.05)

    # Load a baseimage into the generator.
    # The window is the same as the block size to show the blocks used
    generator.load_image('rgb', ((20, 20),))

    fig = plt.figure(figsize=(8, 8))
    for i, img in enumerate(generator):
        ax = plt.subplot(gs[i])
        ax.imshow(img.filled(0.5))
        ax.axis('off')


.. parsed-literal::

    The generator is 8 by 8  blocks



.. image:: feature_extraction_files/feature_extraction_7_1.png


Calculate all the features and append them to a vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this cell we use the ``extract_features`` function from satsense to
extract all features.

``extract_features`` returns a python generator that we can loop over.
Each invocation of this generator returns the feature vector for one
feature in the order of the features list. The shape of this vector is
(x, y, w, v) where:

    - x is the number of blocks of the generator in the x direction
    - y is the number of blocks of the generator in the y direction
    - w is the number of windows the feature is calculated on
    - v is the length of the feature per window

We use a little numpy reshaping to merge these feature vectors into a
single feature vector of shape (x, y, n) where n is the total length of
all features over all windows. In this example it will be (8, 8, 13)
because:

    - HoG has 5 numbers per window and 2 windows:   10
    - NirNDVI has 1 number per window and 1 window:  1
    - Pantex has 1 number per window and2 windows:   2
    -                                        Total: 13

.. code:: ipython3

    vector = []
    for feature_vector in extract_features(features, generator):
        # The shape returned is (x, y, w, v)
        # Reshape the resulting vector so it is (x, y, w * v)
        # e.g. flattened along the windows and features
        data = feature_vector.vector.reshape(
                    *feature_vector.vector.shape[0:2], -1)
        vector.append(data)
    # dstack reshapes the vector into and (x, y, n)
    # where n is the total length of all features
    featureset = np.dstack(vector)

    print("Feature set has shape:", featureset.shape)


.. parsed-literal::

    Feature set has shape: (8, 8, 13)


Showing the resulting features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below we show the results for the calculated features.

In the result images you can see the edges of the feature vector have
been masked as the windows at the edge of the original image contain
masked values. Furthermore, please keep in mind that the value for the
feature in each block depends on an area around the block.

HoG
^^^

Here is the result of the HoG feature, we display the first value for
each window.

Histogram of Gradients is a feature that first calculates a histogram of
the gradient orientations in the window. Using this histogram 5 values
are calculated. This first value is the 1st heaved central shift moment.
Heaved central shift moments are a measure of spikiness of a histogram.

The other values are: the 2nd heaved central shift moment, the
orientation of the highest and second highest peaks and the sine of the
absolute difference between the highest and second highest peak (this is
1 for right angles).

.. code:: ipython3

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 8))

    ax1.axis('off')
    ax1.imshow(image['rgb'])
    ax1.set_title('Input image')

    ax2.axis('off')
    ax2.imshow(featureset[:, :, 0], cmap="PRGn")
    ax2.set_title('Hog Feature for window {}'.format(two_windows[0]))

    ax3.axis('off')
    ax3.imshow(featureset[:, :, 5], cmap="PRGn")
    ax3.set_title('Hog Feature for window {}'.format(two_windows[1]))

    plt.show()



.. image:: feature_extraction_files/feature_extraction_11_0.png


Normalized Difference Vegetation Index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here we show the result for the NDVI feature The NDVI feature captures
the level of vegetation, purple means very little vegetation, green
means a lot of vegetation.

.. code:: ipython3

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    ax1.axis('off')
    ax1.imshow(image['rgb'])
    ax1.set_title('Input image')

    ax2.axis('off')
    ax2.imshow(featureset[:, :, 10], cmap="PRGn")
    ax2.set_title('NirNDVI for window {}'.format(one_window[0]))

    plt.show()



.. image:: feature_extraction_files/feature_extraction_13_0.png


Pantex
^^^^^^

Here we show the results for the Pantex feature. The Pantex feature
captures the level of built-up structures, purple means very little
built-up while green means very built-up.

.. code:: ipython3

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 8))

    ax1.axis('off')
    ax1.imshow(image['rgb'])
    ax1.set_title('Input image')

    ax2.axis('off')
    ax2.imshow(featureset[:, :, 11], cmap="PRGn")
    ax2.set_title('Pantex for window {}'.format(two_windows[0]))

    ax3.axis('off')
    ax3.imshow(featureset[:, :, 12], cmap="PRGn")
    ax3.set_title('Pantex for window {}'.format(two_windows[1]))

    plt.show()



.. image:: feature_extraction_files/feature_extraction_15_0.png


.. _notebooks:

Demonstration Jupyter notebooks
===============================

There are a number of demonstration `Jupyter Notebooks <http://jupyter.org/>`_
available to help you get started with satsense. They can be found in the
`notebooks folder of our github repository <https://github.com/DynaSlum/satsense/tree/master/notebooks>`_.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   feature_extraction
satsense.bands
==============

 .. automodule:: satsense.bands
    :members:
    :undoc-members:
    :show-inheritance:satsense.util
=================

.. automodule:: satsense.util
    :members:
    :undoc-members:
    :show-inheritance:
satsense.features
=================

.. automodule:: satsense.features
    :members:
    :undoc-members:
    :show-inheritance:
satsense.extract
================

 .. automodule:: satsense.extract
    :members:
    :undoc-members:
    :show-inheritance:satsense.generators
===================

 .. automodule:: satsense.generators
    :members:
    :undoc-members:
    :show-inheritance:satsense package
================

Submodules
----------

.. toctree::

    satsense.image
    satsense.bands
    satsense.features
    satsense.generators
    satsense.extract
    satsense.performance
    satsense.util

Module contents
---------------

.. automodule:: satsense
   :members:
   :undoc-members:
   :show-inheritance:
Python API Reference
====================

.. toctree::

   satsense
satsense.image
==============

 .. automodule:: satsense.image
    :members:
    :undoc-members:
    :show-inheritance:
    :special-members: __getitem__satsense.performance
====================

.. automodule:: satsense.performance
    :members:
    :undoc-members:
    :show-inheritance:
