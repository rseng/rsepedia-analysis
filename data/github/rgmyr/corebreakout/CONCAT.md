---
title: 'CoreBreakout: Subsurface Core Images to Depth-Registered Datasets'
tags:
  - Python
  - image processing
  - deep learning
  - geology
  - geoscience
  - subsurface
authors:
  - name: Ross G. Meyer
    orcid: 0000-0003-2344-554X
    affiliation: 1
  - name: Thomas P. Martin
    orcid: 0000-0002-4171-0004
    affiliation: 1
  - name: Zane R. Jobe
    orcid: 0000-0002-7654-4528
    affiliation: 1
affiliations:
 - name: Department of Geology and Geological Engineering, Colorado School of Mines
   index: 1
date: 9 December 2019
bibliography: paper.bib
---

# Summary

Core samples -- cylindrical rock samples taken from subsurface boreholes -- are commonly used by Earth scientists to study geologic history and processes. Core is usually cut into one-meter segments, slabbed lengthwise to expose a flat surface, and stored in cardboard or wooden boxes which are then photographed to enable remote inspection. Unlike other common sources of borehole data [e.g., well logs @Rider:2011], core is the only data that preserves true geologic scale and heterogeneity.

A geologist will often describe core by visual inspection and hand-draw a graphic log of the vertical changes in grain size and other rock properties [e.g., @Jobe:2017]. This description process is time consuming and subjective, and the resulting data is analog. The digitization and structuring of core image data allows for the development of automated and semi-automated workflows, which can in turn facilitate quantitative analysis of the millions of meters of core stored in public and private repositories around the world.

``corebreakout`` is a Python package that provides two main functionalities: (1) a deep learning workflow for transforming raw images of geological core sample boxes into depth-registered datasets, and (2) a `CoreColumn` data structure for storing and manipulating the depth-registered image data. The former uses the Mask R-CNN algorithm [@He:2017] for instance segmentation, and is built around the open source TensorFlow and Keras implementation released by Matterport, Inc. [@Abdulla:2017].


## Mask R-CNN Workflow

The primary user workflow enabled by ``corebreakout`` is depicted in Figure 1. It is straightforward for geologists to add their own labeled training images using ``LabelMe`` [@Wada:2016; @Russell:2007], configure and train new Mask R-CNN models on the labeled images, and subsequently use the trained models to process their own unlabeled images and compile depth-aligned datasets.

![Primary User Workflow](JOSS_figure_workflow.png)

In the future, we would like to train a more generalized model, but for now we anticipate that most users will have to train their own segmentation models. We have found labeling 25-30 images to be the point of diminishing returns for segmentation accuracy, but this number is likely dependent on the consistency of image layout and core material within a given dataset.

Trained models can be loaded using the `CoreSegmenter` class, which provides methods for processing images using the model and according to user-specified layout parameters.

## `CoreColumn` Data Structure

The other main piece of functionality provided by `corebreakout` is the `CoreColumn` class, which is a container for depth-registered, single-column images of core material, allowing for joint manipulation of images and associated depth arrays. `CoreColumn`s may be sliced, stacked, and iterated over, and they include saving, loading, and plotting functionality. Usage details can be found in the documentation and the provided `CoreColumn` tutorial.   

## General Functionality

``corebreakout`` supports standard vertical and horizontal core image layouts, and provides several methods for measuring and assigning depths to core sample columns, including by labeling arbitrary "measuring stick" objects (e.g., rulers, empty trays). We provide a labeled dataset courtesy of the [British Geological Survey's OpenGeoscience project](https://www.bgs.ac.uk/data/bmd.html), as well as a Mask R-CNN model trained on this dataset for testing and demonstration.

In addition to the core Python package, the source code includes scripts for training models, extracting text meta-data from images with optical character recognition [@Smith:2007], and processing directories of images with saved models.

``corebreakout`` has been used to compile a large image dataset for ongoing work in image-based lithology classification [@Martin:2019]. We plan to release our modeling code as a separate project that uses the `CoreColumn` data structure to combine depth-registered image data, sampled well log data, and interval labels into multi-modal datasets for sequence prediction.


# Acknowledgements

We would like to acknowledge the contribution of open source subsurface core images from the British Geological Survey (https://bgs.ac.uk/), and financial support from Chevron through the Chevron Center of Research Excellence at the Colorado School of Mines (https://core.mines.edu/).


# References
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## 0.4 [Unreleased]

### To-Do

- Refactor `PolygonDataset` to use `labelme`'s `imageData` field rather than require jpegs

## 0.3

### Changed

- Increased flexibility of included scripts
- Converted `docs/*.md` files to `*.rst`
- Minor paper and bibliography corrections/fixes

### Added

- Built `sphinx`-generated html docs, hosted on `readthedocs`
- Implemented `CoreColumn.iter_chunks()`
- Added demo/inspection notebooks:
  - `column_demo.ipynb`
  - `segmenter_demo.ipynb`
  - `inspect_dataset.ipynb`
- Added explanatory Markdown to `select_model.ipynb`
- Added restriction in `setup.py` to require `tensorflow` before install

## 0.2

### Changed

- Modified code format using `black`
- `README.md` installation instructions
- moved paper files to root directory

### Added

- `requirements.txt` (except `imgaug`, which is included in `setup.py`)
- CONTRIBUTING, LICENSE, and CODE_OF_CONDUCT files
- `tests/notebooks`
- tests for `CoreColumn` saving and loading

### Removed

- superfluous notebooks and scripts

## 0.1

Initial private version.
# CoreBreakout

[![status](https://joss.theoj.org/papers/add2021f95268fd4cd2850b105f3d570/status.svg)](https://joss.theoj.org/papers/add2021f95268fd4cd2850b105f3d570)

Requirements, installation, and contribution guidelines can be found below. Our full usage and API documentation can be found at: [corebreakout.readthedocs.io](https://corebreakout.readthedocs.io/en/latest/)

### Overview

`corebreakout` is a Python package built around [matterport/Mask\_RCNN](https://github.com/matterport/Mask_RCNN) for the segmentation and depth-alignment of geological core sample images. It provides utilities and an API to enable the workflow depicted in the figure below, as well as a `CoreColumn` data structure to manage and manipulate the resulting depth-registered image data:

![](JOSS_figure_workflow.png)

We are currently using this package to enable research on [Lithology Prediction of Slabbed Core Photos Using Machine Learning Models](https://figshare.com/articles/Lithology_Prediction_of_Slabbed_Core_Photos_Using_Machine_Learning_Models/8023835/2), and are working on getting a DOI for the project through the [Journal of Open Source Software](https://joss.theoj.org/).

## Getting Started

### Target Platform

This package was developed on Linux (Ubuntu, PopOS), and has also been tested on OS X. It may work on other platforms, but we make no guarantees.

### Requirements

In addition to Python`>=3.6`, the packages listed in [requirements.txt](requirements.txt) are required. Notable exceptions to the list are:

- `1.3<=tensorflow-gpu<=1.14` (or possibly just `tensorflow`)
- `mrcnn` via [submodule: matterport/Mask\_RCNN](https://github.com/matterport/Mask_RCNN/tree/3deaec5d902d16e1daf56b62d5971d428dc920bc)

The TensorFlow requirement is not explicitly listed in `requirements.txt` due to the ambiguity between `tensorflow` and `tensorflow-gpu` in versions `<=1.14`. The latter is almost certainly required for training new models, although it may be possible to perform inference with saved models on CPU, and use of the `CoreColumn` data structure does not require a GPU.

Note that TensorFlow GPU capabilities are implemented with [CUDA](https://developer.nvidia.com/cuda-zone), which requires a [supported NVIDIA GPU](https://developer.nvidia.com/cuda-gpus).

#### Additional (Optional) Requirements

Optionally, `jupyter` is required to run demo and test notebooks, and `pytest` is required to run unit tests. Both of these should be manually installed if you plan to modify or contribute to the package source code.

We also provide a script for extraction of top/base depths from core image text using `pytesseract`. After installing the [Tesseract OCR Engine](https://github.com/tesseract-ocr/tesseract) on your machine, you can install the `pytesseract` package with `conda` or `pip`.

### Download code

```
$ git clone --recurse-submodules https://github.com/rgmyr/corebreakout.git
$ cd corebreakout
```

### Download data (optional)

To make use of the provided dataset and model, or to train new a model starting from the pretrained COCO weights, you will need to download the `assets.zip` folder from the [v0.2 Release](https://github.com/rgmyr/corebreakout/releases/tag/v0.2).

Unzip and place this folder in the root directory of the repository (its contents will be ignored by `git` -- see the `.gitignore`). If you would like to place it elsewhere, you should modify the paths in [corebreakout/defaults.py](https://github.com/rgmyr/corebreakout/blob/master/corebreakout/defaults.py) to point to your preferred location.

The current version of `assets/data` has JSON annotation files which include an `imageData` field representing the associated images as strings. For now you can delete this field and reduce the size of the data with `scripts/prune_imageData.py`:

```
$ python scripts/prune_imageData.py assets/
```

### Installation

We recommend installing `corebreakout` and its dependencies in an isolated environment, and further recommend the use of `conda`. See [Conda: Managing environments](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

---

To create a new `conda` environment called `corebreakout-env` and activate it:

```
$ conda create -n corebreakout-env python=3.6 tensorflow-gpu=1.14
$ conda activate corebreakout-env
```

**Note:** If you want to try a CPU-only installation, then replace `tensorflow-gpu` with `tensorflow`. You may also lower the version number if you are on a machine with `CUDA<10.0` (required for TensorFlow`>=1.13`). See [TensorFlow GPU requirements](https://www.tensorflow.org/install/gpu#software_requirements) for more compatibility details.

---

Then install the rest of the required packages into the environment:

```
$ conda install --file requirements.txt
```

---

Finally, install `mrcnn` and `corebreakout` using `pip`. Develop mode installation (`-e`) is recommended (but not required) for `corebreakout`, since many users will want to change some of the default parameters to suit their own data without having to reinstall afterward:

```
$ pip install ./Mask_RCNN
$ pip install -e .
```

## Usage

Please refer to our [readthedocs page](https://corebreakout.readthedocs.io/en/latest/) for full documentation!

## Development and Community Guidelines

### Submit an Issue

- Navigate to the repository's [issue tab](https://github.com/rgmyr/corebreakout/issues)
- Search for existing related issues
- If necessary, create and submit a new issue

### Contributing

- Please see [`CONTRIBUTING.md`](CONTRIBUTING.md) and the [Code of Conduct](CODE_OF_CONDUCT.md) for how to contribute to the project

### Testing

- Most `corebreakout` functionality not requiring trained model weights can be verified with `pytest`:

```
$ cd <root_directory>
$ pytest .
```

- Model usage via the `CoreSegmenter` class can be verified by running `tests/notebooks/test_inference.ipynb` (requires saved model weights)
- Plotting of `CoreColumn`s can be verified by running `tests/notebooks/test_plotting.ipynb`
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at ross.meyer@utexas.edu. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# Contributing

1. Fork it!
2. Create your feature branch
3. Pass all `pytest` tests and verify `tests/notebooks/` run without errors
4. Update `CHANGELOG.md`
    - Describe changes under `[Unreleased]` section
5. Submit a pull request

These scripts contain utilities for generating and handling Pascal-VOC style XML data labels (create labels with e.g., `LabelImg`). We have used this scheme to assign interval labels to aggregated `CoreColumn` images, for use in downstream ML pipelines (related code will be released concurrently with a paper at a later date.)

This is entirely separate from `corebreakout` functionality, but may be useful to people interested in undertaking similar projects (feel free to raise an issue or email ross.meyer<at>utexas.edu for guidance). We may also integrate some of this functionality into the main package (and have started to do so, with *e.g.*, `CoreColumn.iter_chunks()`).

Basically, we split saved column image arrays such that they are small enough to be saved as JPEGs (which have a dimension limit of `2^16`). We then assign labels to intervals in `LabelImg`(link here), and finally snap + concatenate the XML labels and convert to a row-wise labels array.

### `split_npy_image.py`

Takes `image.npy` and `depth.npy` files from `src`, writes a set of images (under jpeg size limit) to `<dst>/well/`.

### `join_xml_labels.py`

Concatenates XML label files from a `src` directory, writes a row labels file (required for modeling) to `dst`.
Scripts Reference
=================

``train_mrcnn_model.py``
------------------------
Train a new ``mrcnn`` model starting from pretrained COCO weights

optional arguments:
  -h, --help            show this help message and exit
  --steps STEPS         1, 2, or 3. How many of the steps to train: (heads,
                        4+, entire model)
  --model_dir MODEL_DIR
                        Directory in which to create new training
                        subdirectory for checkpoints and tensorboard logs.
  --data_dir DATA_DIR   Directory in which to find ``train`` and ``test``
                        subdirectories containing labeled images.

                        Extract OCR top and base depths from images with `pytesseract`

``get_ocr_depths.py``
---------------------
Extract OCR top and base depths from images with ``pytesseract``

optional arguments:
  -h, --help            show this help message and exit
  --root_dir ROOT_DIR   A common parent directory of all target ``<subdir>`` directories.
  --subdir SUBDIR       A string contained in the name of all target subdirectories.
  --save_name SAVE_NAME
                        Name of depths csv file(s) to be saved in matching subdirs.
  --force               Flag to force overwrite of any existing ``<save_name>.csv`` files.
  --inspect             Flag to inspect images and print OCR output whenever there is an issue.

As an example, you can test the script by running:

.. code:: none

  $ cd scripts
  $ python get_ocr_depths.py --root_dir ../tests/data --subdir two_image_dataset --save_name auto_depths_test

This should save a new file at ``tests/data/two_image_dataset/auto_depths_test.csv``, with contents like:

.. code::

  ,top,bottom
  S00101409.jpeg,2348.0,2350.0
  S00111582.jpeg,7716.0,2220.0

Note that ``7716.0`` is a misread, and should have been ``2218.0``. At least with our BGS images, some manual corrections are usually required, but this provides a template for the ``--depth_csv`` file required to run ``process_directory.py``.


``process_directory.py``
------------------------
Process directory of raw images with Mask R-CNN and save results as a ``CoreColumn``.

The ``path`` given should contain images as jpeg files, and a ``depth_csv`` file in the format:

.. code::

           ,    top,    bottom
  <filename1>, <top1>, <bottom1>
  ...
  <filenameN>, <topN>, <bottomN>

**NOTE**: model ``Config``, ``class_names``, and segmentation ``layout_params`` can only be
changed manually at the top of script, and default to those configured in `defaults.py <https://github.com/rgmyr/corebreakout/blob/master/corebreakout/defaults.py>`_

positional arguments:
  path                 Path to directory of images (and depth information csv) to process.

optional arguments:
 -h, --help            show help message and exit
 --model_dir MODEL_DIR
                       Directory to load ``mrcnn`` model from.
                       Default=``defaults.MODEL_DIR``
 --weights_path WEIGHTS_PATH
                       Path to model weights to load.
                       Default=``defaults.CB_MODEL_PATH``
 --add_tol ADD_TOL     Gap tolerance when adding ``CoreColumn`` objects,
                       default=5.0.
 --add_mode ADD_MODE   ``CoreColumn.add_mode``. One of {'fill', 'collapse'}.
 --depth_csv DEPTH_CSV
                       Name of filename + (top, bottom) csv to read from
                       ``path``, default=``'auto_depths.csv'``
 --save_dir SAVE_DIR   Path to save ``CoreColumn`` to, default=None will save to
                       ``path``
 --save_name SAVE_NAME
                       Name to use for ``CoreColumn.save``, default=None
                       results in ``CoreColumn_<top>_<base>``
 --save_mode SAVE_MODE
                       One of {'pickle', 'numpy'}. Whether to save as single
                       ``pkl`` file or multiple ``npy`` files

Assuming you've downloaded and unzipped the `assets` folder in the default location, you can test the script with default parameters by running:

.. code::

  $ cd scripts
  $ python process_directory.py ../tests/data/two_image_dataset --depth_csv dummy_depths.csv

This should save the aggregated ``CoreColumn`` to ``tests/data/two_image_dataset/CoreColumn_1.00_5.00.pkl``.

``prune_imageData.py``
----------------------
Remove the ``imageData`` field from all JSON files in tree below ``path``:

positional arguments:
  path        Path to parent of all target JSON files.

optional arguments:
  -h, --help  show this help message and exit
.. _model-building:

Building Mask R-CNN Models
==========================

For best results, most users will want to train models on some of
their own data. See :ref:`creating-datasets` for guidelines.

We may release a more general pretrained model in the future, pending demand and the availability of open source datasets. Feel free to contact us if you would like to contribute your data for that purpose.

The rough outline of model construction and training looks like this:

.. code:: python

   import mrcnn.model as modellib
   from corebreakout import defaults, datasets

   model_config = defaults.DefaultConfig()

   train_dataset = datasets.PolygonDataset(...)
   test_dataset = datasets.PolygonDataset(...)

   model = modellib.MaskRCNN(...)

   model.train(train_dataset, test_dataset, ...)

For the finer details see `scripts/train_mrcnn_model.py <https://github.com/rgmyr/corebreakout/blob/master/scripts/train_mrcnn_model.py>`__

``mrcnn`` Model Configuration
-----------------------------

Models are created with a subclass of ``mrcnn.config.Config``. See
``corebreakout.defaults.DefaultConfig`` for our latest model
configuration.

To see all available configuration parameters, see
`mrcnn/config.py <https://github.com/matterport/Mask_RCNN/blob/3deaec5d902d16e1daf56b62d5971d428dc920bc/mrcnn/config.py>`__

The obvious parameters that a user might want to change include:

  - ``NAME``
  - ``RPN_ANCHOR_RATIOS`` : defaults are set up for horizontal columns. Something like ``[1.0, 3.0, 7.0]`` would make more sense for vertical columns.
  - ``STEPS_PER_EPOCH``, ``VALIDATION_STEPS`` : batches per training epoch and validation step. Does not necessarily need to match dataset sizes.
  - ``IMAGES_PER_GPU`` : you can try to increase if you have a large GPU.
  - ``GPU_COUNT`` : you can increase if you multiple GPUs.

You can either modify ``DefaultConfig`` directly, or instantiate your own ``Config`` subclass.

Model Training
--------------

The simplest way to train a model is by running
`scripts/train_mrcnn_model.py <https://github.com/rgmyr/corebreakout/blob/master/scripts/train_mrcnn_model.py>`__.
This script loads pretrained ``COCO`` weights, and executes a three step
training + tuning run.

It also serves as a demonstration of ``Dataset`` collection, and
instantiating and training a ``mrcnn.modellib.MaskRCNN``.

While training, ``mrcnn`` logs ``tensorboard`` files to the specified
``model_dir``. You can view the files by running:

::

   $ tensorboard --logdir <model_dir>

Model Selection
---------------

We recommend viewing the ``tensorboard`` files (and particularly the
``val_loss`` scalar) to select candidate models.

`notebooks/select_model.ipynb <https://github.com/rgmyr/corebreakout/blob/master/notebooks/select_model.ipynb>`__
provides a template for viewing the output of candidate models on the
test dataset.

**Note**: ``mrcnn`` saves Checkpoints each epoch starting at ``0001``,
while ``tensorboard`` logs epochs starting from ``0``. So, if epoch
``X`` looks good on ``tensorboard``, you will want to reference epoch
``X+1`` in your list of candidates to load the corresponding weights.

Using a Model
-------------

Once you have trained and selected a new model, you may want to change
the ``*PATH`` variables in ``corebreakout/defaults.py`` to point to the location of the
new model weights (these paths are what get referenced by default in
``scripts/process_directory.py``, etc.).

Alternatively, you can always pass whatever ``model_dir`` and ``weights_path`` (and
``model_config`` instance) you like when constructing a ``CoreSegmenter``.
.. CoreBreakout documentation master file, created by
   sphinx-quickstart on Fri Apr 17 06:21:33 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CoreBreakout's documentation!
========================================

Project Repository: `rgmyr/corebreakout <https://github.com/rgmyr/corebreakout>`_

``mrcnn`` Repository: `matterport/Mask_RCNN <https://github.com/matterport/Mask_RCNN>`_

``corebreakout`` provides two main functionalities: **(1)** a deep learning workflow for transforming raw images of geological core sample boxes into depth-registered datasets, which is facilitated by the ``CoreSegmenter`` API, and **(2)** a ``CoreColumn`` data structure for storing and manipulating depth-registered image data.

|workflowfigure|

This documentation covers usage of ``corebreakout``. To dig into the finer details of the Mask R-CNN model implementation, please refer to ``mrcnn``'s documentation.

Provided Data
=============

We provide labeled images from the British Geological Survey's North Sea collection, as well as a saved Mask R-CNN model trained on this data, and a pretrained COCO model from which to start new training runs.

To use the data or models, download the ``assets.zip`` folder from the `Releases Page <https://github.com/rgmyr/corebreakout/releases>`_. Unzip it in the root directory of the project, or modify the paths in `defaults.py <https://github.com/rgmyr/corebreakout/blob/master/corebreakout/defaults.py>`_ to point to the location.

JSON annotations in this data currently contain a superfluous field called ``imageData`` which takes up most of the memory for these files. You can delete this field and reduce the file sizes with ``scripts/prune_imageData.py``:

.. code::

  $ python scripts/prune_imageData.py assets/

If you want to use your own data, then label some images for Mask R-CNN training, following the guidelines in :ref:`creating-datasets`. We recommend starting with 20-30 images.

Overview
========

Image Processing Workflow
-------------------------

(1) If you're looking to use your own data, you will probably need to label some of your images for best results. Follow the guidelines in :ref:`creating-datasets`.

(2) Train a Mask R-CNN model using your labeled data. Model configuration, training, and selection are explained in :ref:`model-building`.

(3) Use the trained model to process directories of unlabeled images and save the results as a ``CoreColumn``. We provide `scripts/process_directory.py <https://github.com/rgmyr/corebreakout/blob/master/scripts/process_directory.py>`__ to facilitate this step. It requires saved model weights, and a csv file listing the top and bottom depths for each image in the directory. To make creating these csv's easier, we provide `scripts/get_ocr_depths.py <https://github.com/rgmyr/corebreakout/blob/master/scripts/get_ocr_depths.py>`_.

The ``CoreSegmenter`` Class
---------------------------

`scripts/process_directory.py <https://github.com/rgmyr/corebreakout/blob/master/scripts/process_directory.py>`__ uses the ``corebreakout.CoreSegmenter`` API to handle converting images to ``CoreColumns``, and you may also use this class directly:

.. code:: python

  segmenter = corebreakout.CoreSegmenter(
        model_dir,
        weights_path,
        model_config  = corebreakout.defaults.DefaultConfig,
        class_names   = corebreakout.defaults.CLASSES,
        layout_params = corebreakout.defaults.LAYOUT_PARAMS
  )

  # `img` can be an array or a path to an image
  column = segmenter.segment(img, [top_depth, base_depth], **kwargs)

  # for iterables of images (or paths) and depth range pairs
  column = segmenter.segment_all(imgs, depth_ranges, **kwargs)

``class_names`` should correspond to those in the dataset on which the model was trained, and ``layout_params`` are explained in detail in the :ref:`layout-parameters` documentation.

The ``CoreColumn`` Class
------------------------

This object is a container for depth-registered image data. Columns can be added, sliced, plotted, iterated over in chunks, saved, and loaded (in either single-file ``.pkl`` format or multi-file ``.npy`` format).

For a demonstration of the ``CoreColumn`` API, see: `notebooks/column_demo.py <https://github.com/rgmyr/corebreakout/blob/master/notebooks/column_demo.ipynb>`__

.. toctree::
   :maxdepth: 2
   :caption: Usage Documentation:

   creating_datasets
   model_building
   layout_parameters
   scripts_reference


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. |workflowfigure| image:: images/JOSS_figure_workflow.png
.. _creating-datasets:

Creating ``Datasets``
=============================

Using ``labelme``
-----------------

The recommended way to add a new set of labeled training images is to annotate them using `wkentaro/labelme <https://github.com/wkentaro/labelme>`__. The ``labelme`` GUI allows the user to draw any number of labeled polygons on an image, and saves the labels and coordinates in a JSON annotation file.

You can easily install ``labelme`` with ``pip``, and then open our small test dataset:

.. code::

  $ pip install labelme
  $ labelme tests/data/two_image_dataset

This should open a window that looks like this:

|image0|

New polygons can be drawn by clicking ``Create Polygons`` and then clicking points on the image to add vertices. When you get back around to the starting vertex to close the loop of points and create a polygon, ``labelme`` will ask you to assign it a text label.

Unique text labels are aggregated on the right side, and the individual polygons shown in the image are colored and listed with their assigned labels. They happen to be the same in this example, but there could be more labels in the ``Label List`` than polygons in the ``Polygon Labels`` list for any given image in the directory.

When you finish adding polygons, click the ``Save`` button to export the annotations as a ``.json`` file, which by default will save to the current directory.

The ``File List`` shows images in the currently directory, with the checkbox indicating whether an annotation has been saved for that file yet. You can use the ``Next Image`` and ``Previous Image`` buttons to traverse through the images (and corresponding annotations) in the directory.


Inspect Our Annotations
-----------------------

For a more hands-on exploration of annotation structure and the ``Dataset`` API, see `notebooks/inspect_dataset.py <https://github.com/rgmyr/corebreakout/blob/master/notebooks/inspect_dataset.ipynb>`__


How to Name Things
------------------

As you can see above, our text labels for individual objects have the form ``<class><number>``. The ``class`` refers to a general object type, the ``number`` differentiates between individual `instances` of that class within an image. This isn't the only possible naming scheme, but it's the one you should use if you want your annotations to be compatible with our ``PolygonDataset`` implementation.

In our data, we use two classes, **col** and **tray**, to label two different types of objects:

- **col** is for `columns` of core material -- the individual objects that we want to segment, crop, and stack on top of each other.
- **tray** is for empty trays, which for our data are the most reliable "measuring stick" to use for defining the "top" and "base" positions of columns within an image

The first type of object is required -- the whole purpose of this workflow is to segment core columns, so you will need to at least label those.

The second type of object is optional -- using another object (empty tray, measuring stick, etc.) to define top and base positions of columns is `one` option for the ``endpts`` setting in :ref:`layout-parameters`, but there are others which do not require the existence of any objects other than core columns.

If you're making your own dataset, your class names can be whatever you want them to be. You just need at least some name for core columns, and to keep general consistency between:

- The ``PolygonDataset`` (s) you train your model with, which use the annotations you saved.
- The ``kwargs`` you provide when instantiating a ``CoreSegmenter``. Namely:

  - The ``classes`` argument (it should be the same as for ``PolygonDataset``)
  - The ``col_class`` (and possibly ``endpts``) fields of the ``layout_params``

If you decide to use different class names (or layout parameters), we'd recommend changing the ``CLASSES`` and ``LAYOUT_PARAMS`` in `defaults.py <https://github.com/rgmyr/corebreakout/blob/master/corebreakout/defaults.py>`__ so that you don't have to specify them as often.

Removing ``imageData``
----------------------

``labelme`` saves a field called ``imageData`` which encodes the entire image as a string, consuming quite a bit of unnecessary memory. Our ``PolygonDataset`` class doesn't make use of this field, and we've provided a script to delete it from all JSON files below a root ``path``:

.. code::

  $ python scripts/prune_imageData.py <path>


Summary of Guidelines
---------------------

To be able to use the built-in ``corebreakout.datasets.PolygonDataset`` class with your training data, you want to adhere to these guidelines:

-  Save ``<fname>.json`` annotations in a flat directory with
   corresponding ``<fname>.jpeg`` files (this is ``labelme``\ ’s default
   behavior)
-  You may label any number of classes. You will have to supply a list
   of these classes to the ``PolygonDataset`` and ``CoreSegmenter`` constructors, or modify
   ``defaults.CLASSES``.

   - At a minimum, you will want some class name to represent columns of core. We call ours ``'col'``, but it doesn't matter what you want to call it if you're using your own data.
   - You may also create a class (*e.g.*, ``'tray'``) for any objects that consistently demarcate the top and base positions of core columns better than the columns themselves.
-  Different instances of the same class should begin with the class
   name and be differentiated afterward (*e.g.*, ``col1, col2, col3``)

   -  The corollary is that no class name can be a substring of any
      other class name (*e.g.*, ``col, col_tray`` would not be allowed)
   -  Multiple polygons may belong to a single instance (for example, if there's a large gap in the middle of a column)

-  After annotating images, split into sibling ``'train'`` and
   ``'test'`` subdirectories

**Note:** We've found that the point of diminishing returns happens somewhere in the range of 20-30 training images, which probably corresponds to 30-50 column instances for this dataset. YMMV.

After compiling the annotations, you may wish to modify ``defaults.DATASET_DIR`` to avoid need to explicitly specify the data location.

``corebreakout.datasets.PolygonDataset``
----------------------------------------

This is a subclass of ``mrcnn.utils.Dataset`` for instance segmentation
annotations in the default JSON format of
`wkentaro/labelme <https://github.com/wkentaro/labelme>`__.

Usage
~~~~~

::

   from corebreakout.datasets import PolygonDataset

   data_dir = defaults.DEFAULT_DATA_DIR    # parent of any separate annotation data directories
   subset = 'train'                        # which subdirectory to read from

   dataset = PolygonDataset(classes=defaults.DEFAULT_CLASSES)

   # Collect all of the requied ID + path information
   dataset.collect_annotated_images(data_dir, subset)

   # Set all of the attrs required for use
   dataset.prepare()

   print(dataset)

Two ``Dataset`` objects (train, test) are required in calls to ``model.train(...)``, which is why we split them into separate directories.

Subclassing ``mrcnn.utils.Dataset``
-----------------------------------

If you want to use a different annotation format, you can inherit from
the base ``mrcnn.utils.Dataset`` class.

You will need to write some user-called method to collect file
information:

- *e.g.*, ``collect_annotated_images(data_dir, subset)``: Register ``image_id``, ``path``, and ``ann_path`` for each (image, annotation) file pair in ``<data_dir>/<subset>`` directory.

And then override at least these two methods:

- ``load_mask(image_id)``: Given an ``image_id``, load (and compute, if necessary) the corresponding mask. For an with ``N`` objects (not including the background), the return value from this function should be ``(mask, class_ids)``, where ``mask`` is boolean array of shape ``(H,W,N)`` and ``class_ids`` is a 1D integer array of size ``N`` with one ``class_id`` for each channel in ``mask``.
- ``image_reference(image_id)``: Return the path of an image, a link to it, or some other unique property to help in looking it up or debugging it.

.. |image0| image:: images/labelme1.png
.. _layout-parameters:

Layout Parameters
=================

``corebreakout`` allows for processing core images with different
layouts, and provides several different methods for finding and
cropping the bounding boxes of individual columns within an image.

Layout parameters are used within the ``CoreSegmenter.segment()``
method.

Any and all default parameters can be overridden and updated in the
``CoreSegmenter`` constructor, or in any single call to the
``segment()`` method. In either case, pass any new parameters in a dict
as the ``layout_params`` keyword argument.

The default parameters reflect the characteristics of the BGS dataset on
which the example model is trained:

.. code:: python

   # corebreakout/defaults.py

   LAYOUT_PARAMS = {
     'order' : 't2b',
     'orientation' : 'l2r',
     'col_height' : 1.0,
     'col_class' : 'col',
     'endpts' : 'tray'
   }

**DEVELOPMENT NOTE:** The only allowed values for ``order`` and
``orientation`` are ``'t2b'`` and ``'l2r'``. This covers all
conventional core image layouts that we are aware of, but we would
consider adding ``'b2t'`` and ``'r2l'`` if provided with use-cases. If
you have one, please open an issue (or submit a pull request :-).

``order`` and ``orientation``
-----------------------------

The ``'order'`` parameter specifies the depth order by which to sort the
set of columns detected in an image:

- ``'t2b'`` implies that columns are laid out horizontally, with the uppermost column coming first in order of depth.
- ``'l2r'`` implies that columns are laid out vertically, with the leftmost column coming first in order of depth.

The ``'orientation'`` parameter specifies the depth orientation of
columns. This should be the converse of ``'order'``:

- ``'t2b'`` implies that the top of a (vertical) column is toward the top of the image.
- ``'l2r'`` implies that the top of a (horizontal) column is toward the left side of the image.

Since it is required that ``order`` be one of these options and
``orientation`` be the other, requiring both *is* redundant. However,
making both of them explicit improves code readability, and will make it
easier to add other options should we choose to do so in future
releases.

``col_height``
--------------

The ``'col_height'`` parameter specifies the height in depth units
(usually meters or feet) of individual (and complete) columns.

This value is used in conjunction with the ``depth_range`` positional
argument to find the number of expected columns in an image when calling
``CoreSegmenter.segment(img, depth_range, **kwargs)``.

``col_class``
-------------

The name of the class representing core sample columns in the M-RCNN
model used by ``CoreSegmenter`` instance.

``endpts``
----------

The ``'endpts'`` parameter determines the method for making sure that
before cropping, the ``top`` and ``base`` of partial columns are
extended to locations that are approximately ``'col_height'`` apart.
Different options may work better or worse depending on how clean the
samples are and how consistent the layout is.

Predicted masks tend to be subsets of the ‘true’ masks, so **short
columns are extended**, but **longer columns are NOT shortened**. You
can see this in the example images below, where the computed minimal
endpoint locations are shown as solid yellow lines, and the resulting
column bounding boxes are shown as green dashed lines.

Allowed values of ``'endpts'`` include:

- The name of a class (*e.g.*, ``'tray'``)
    - Results in columns being extended to the ``top`` and ``base`` of the strongest detection of this class
    - Must be found in the ``class_names`` attribute of the ``CoreSegmenter`` instance
    - Typical choices would be empty trays, or the measuring sticks commonly placed next to boxes of core

|image0|

-  One of the keywords ``'auto'`` or ``'auto_all'``

   -  Results in columns being extended to the min/max coordinates of a
      set of detected objects
   -  ``'auto'`` will use only objects of ``'col_class'`` as the
      relevant set (*e.g.*, all ``'col'`` detections – first example
      below)
   -  ``'auto_all'`` will use all objects in the image (*e.g.*, all
      ``'col'`` **and** ``'tray'`` detections – second example below)

|image1|

|image2|

-  A 2-tuple of explicit integer endpoint coordinates (*e.g.*,
   ``(100, 6900)``)

   -  Results in columns being extended to at least these min/max
      coordinates

|image3|

.. |image0| image:: images/endpts_tray.png
.. |image1| image:: images/endpts_auto.png
.. |image2| image:: images/endpts_auto_all.png
.. |image3| image:: images/endpts_explicit.png
corebreakout package
====================

Subpackages
-----------

.. toctree::

   corebreakout.datasets

Submodules
----------

corebreakout.column module
--------------------------

.. automodule:: corebreakout.column
   :members:
   :undoc-members:
   :show-inheritance:

corebreakout.defaults module
----------------------------

.. automodule:: corebreakout.defaults
   :members:
   :undoc-members:
   :show-inheritance:

corebreakout.segmenter module
-----------------------------

.. automodule:: corebreakout.segmenter
   :members:
   :undoc-members:
   :show-inheritance:

corebreakout.utils module
-------------------------

.. automodule:: corebreakout.utils
   :members:
   :undoc-members:
   :show-inheritance:

corebreakout.viz module
-----------------------

.. automodule:: corebreakout.viz
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: corebreakout
   :members:
   :undoc-members:
   :show-inheritance:
corebreakout.datasets package
=============================

Submodules
----------

corebreakout.datasets.polygondataset module
-------------------------------------------

.. automodule:: corebreakout.datasets.polygondataset
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: corebreakout.datasets
   :members:
   :undoc-members:
   :show-inheritance:
corebreakout
============

.. toctree::
   :maxdepth: 4

   corebreakout
