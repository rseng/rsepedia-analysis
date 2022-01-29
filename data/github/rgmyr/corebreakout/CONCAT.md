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
