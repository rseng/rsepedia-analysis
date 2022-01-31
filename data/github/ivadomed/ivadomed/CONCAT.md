---
title: 'ivadomed: A Medical Imaging Deep Learning Toolbox'
tags:
  - Deep Learning
  - Medical Imaging
  - Segmentation
  - Open-source
  - Pipeline
authors:
  - name: Charley Gros
    orcid: 0000-0003-4318-0024
    affiliation: "1, 2"
  - name: Andreanne Lemay
    orcid: 0000-0001-8581-2929
    affiliation: "1, 2"
  - name: Olivier Vincent
    orcid: 0000-0002-5554-8108
    affiliation: "1, 2"
  - name: Lucas Rouhier
    affiliation: 1
  - name: Marie-Helene Bourget
    orcid: 0000-0001-9786-9332
    affiliation: 1
  - name: Anthime Bucquet
    affiliation: 1
  - name: Joseph Paul Cohen
    affiliation: "2, 3"
  - name: Julien Cohen-Adad
    orcid: 0000-0003-3662-9532
    affiliation: "1, 2, 4"
affiliations:
 - name: NeuroPoly Lab, Institute of Biomedical Engineering, Polytechnique Montreal, Montreal, Canada
   index: 1
 - name: Mila, Quebec AI Institute, Montreal, QC, Canada
   index: 2
 - name: AIMI, Stanford University, Stanford, CA, USA
   index: 3
 - name:  Functional Neuroimaging Unit, CRIUGM, Université de Montréal, Montreal, QC, Canada
   index: 4
date: 3 November 2020
bibliography: paper.bib

---

# Summary

`ivadomed` is an open-source Python package for designing, end-to-end training, and evaluating deep learning models applied to medical imaging data. The package includes APIs, command-line tools, documentation, and tutorials. `ivadomed` also includes pre-trained models such as spinal tumor segmentation and vertebral labeling. Original features of `ivadomed` include a data loader that can parse image and subject metadata for custom data splitting or extra information during training and evaluation. Any dataset following the [Brain Imaging Data Structure (BIDS)](https://bids.neuroimaging.io/) convention will be compatible with `ivadomed`. Beyond the traditional deep learning methods, `ivadomed` features cutting-edge architectures, such as FiLM [@perez2017film] and HeMis [@havaei2016hemis], as well as various uncertainty estimation methods (aleatoric and epistemic), and losses adapted to imbalanced classes and non-binary predictions. Example applications of `ivadomed` include MRI object detection, segmentation, and labeling of anatomical and pathological structures. `ivadomed`'s main project page is available at https://ivadomed.org.

# Statement of need

Deep learning applied to medical imaging presents many challenges: datasets are often not publicly-available, ground truth labels are difficult to produce, and needs are usually tailored to particular datasets and tasks [@kim_deep_2019]. There already exists a few deep learning frameworks for medical imaging, but [each has their pros & cons](https://ivadomed.org/en/latest/purpose.html#comparison-with-other-projects). `ivadomed` notably addresses unmet needs in terms of data management, readily-available uncertainty outputs, missing modalities (in case of multi-channel training) and model comparison, to only name few of the original features.

Another challenge of medical imaging is the heterogeneity of the data across hospitals (e.g., contrast, resolution), making it difficult to create models that can generalize well. Recent cutting-edge methods address this problem, such as FiLM [@perez2017film] and HeMis [@havaei2016hemis], however they are usually not implemented within a comprehensive framework. In addition to providing these architectures, `ivadomed` features losses adapted to imbalanced classes and non-binary predictions.

![Overview of `ivadomed` main features.\label{fig:overview}](https://raw.githubusercontent.com/ivadomed/doc-figures/main/index/overview_title.png)

## Loader

**Standard data structure:** In machine learning, lots of time is spent curating the data (renaming, filtering per feature) before entering the training pipeline [@Willemink2020-au]. `ivadomed` features an advanced data loader compatible with a standardized data structure in neuroimaging: the Brain Imaging Data Structure (BIDS) [@bids_2016]. Thus, any dataset following the BIDS convention can readily be used by `ivadomed`. BIDS convention is originally designed around MRI data and accepts NIfTI file formats, but the BIDS community is actively expanding its specifications to other modalities (CT, MEG/EEG, microscopy), which `ivadomed` will then be able to accommodate.

**Access to metadata:** One benefit of the BIDS convention is that each image file is associated with a JSON file containing metadata. `ivadomed`'s loader can parse image metadata (e.g., acquisition parameters, image contrast, resolution) and subject metadata (e.g., pathology, age, sex) for custom data splitting or extra information during training and evaluation. It is possible to modulate specific layers of a convolutional neural network using metadata information to tailor it towards a particular data domain or to enable experiments with architectures such as FiLM [@perez2017film]. Metadata could also be useful to mitigate class imbalance via data balancing techniques.

## Training

**Architectures:** `ivadomed` includes all the necessary components for training segmentation models from start to finish, including data augmentation transformations and transfer learning. Available architectures include: 2D U-Net [@Ronneberger2015unet], 3D U-Net [@isensee2017brain], ResNet [@he2016resnet], DenseNet [@Huang2017densenet], Count-ception [@Cohen2017countception], and HeMIS U-Net. These models can easily be enriched via attention blocks [@oktay2018attention] or FiLM layers (which modulate U-Net features using metadata).

**Losses and class imbalance:** Popular losses are available: Dice coefficient [@milletari2016v], cross-entropy, and L2 norm, including some adapted to medical imaging challenges, such as the adaptive wing loss [@wang_adaptive_2019] for soft labels and the focal Dice loss [@wong20183d] for class imbalance. Useful to alleviate overfitting, mixup [@zhang2017mixup] was modified to handle segmentation tasks. To mitigate class imbalance, `ivadomed` supports cascaded architectures. With a single inference, it is possible to narrow down the region of interest via object detection and then segment a specific structure, as illustrated in \autoref{fig:lemay2020}.

**Getting started:** It can be overwhelming to get started and choose across all the available models, losses, and parameters. `ivadomed`'s repository includes the script `ivadomed_automate_training` to configure and launch multiple trainings across GPUs. In case of interruption during training, all parameters are saved after each epoch so that training can be resumed at any time.

## Evaluation

**Uncertainty:** Aleatoric [@wang_aleatoric_2019] and/or epistemic [@nair_exploring_2018] uncertainties can be computed voxel-wise and/or object-wise [@roy_quicknat_2018]. Multiple metrics are available, including entropy and coefficient of variation.

**Post-processing:** Predictions can be conveniently filtered using popular methods, e.g., fill holes, remove small objects, threshold using uncertainty. It is also possible to compute metrics for specific object sizes (e.g., small vs. large lesions). `ivadomed` has a module to find the optimal threshold value on the soft output prediction, via a grid-search finding applied to evaluation metrics or ROC curve.

**Visualize performance:** Convenient visualization tools are available for model design and optimization: GIF animations across training epochs, visual quality control of data augmentation, training curve plots, integration of the TensorBoard module, and output images with true/false positive labels. See [this example tutorial](https://ivadomed.org/en/latest/tutorials/cascaded_architecture.html#visualize-training-data).

# Usage

Past and ongoing research projects using `ivadomed` are listed [here](https://ivadomed.org/en/latest/use_cases.html). The figure below illustrates a cascaded architecture for segmenting spinal tumors on MRI data [@lemay_fully_2020].

![Fully automatic spinal cord tumor segmentation framework..\label{fig:lemay2020}](https://raw.githubusercontent.com/ivadomed/doc-figures/main/use_cases/lemay_2020.png)

# Acknowledgements

The authors thank Ainsleigh Hill, Giselle Martel, Alexandru Jora, Nick Guenther, Joshua Newton, Christian Perone, Konstantinos Nasiotis, Valentine Louis-Lucas, Benoît Sauty-De-Chalon, Alexandru Foias and Leander Van Eekelen for their useful contributions, and Guillaume Dumas for proof-reading the manuscript. Funded by IVADO, the Canada Research Chair in Quantitative Magnetic Resonance Imaging [950-230815, 950-233166], CIHR [FDN-143263], CFI [32454, 34824], FRQS [28826], NSERC [RGPIN-2019-07244], FRQNT [2020‐RS4‐265502 UNIQUE] and TransMedTech. C.G. has a fellowship from IVADO [EX-2018-4], A.L. has a fellowship from NSERC, FRQNT and UNIQUE, O.V. has a fellowship from NSERC, FRQNT and UNIQUE, M.H.B. has a fellowship from IVADO and FRQNT.

# References
  
![ivadomed Overview](https://raw.githubusercontent.com/ivadomed/doc-figures/main/index/overview_title.png)

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02868/status.svg)](https://doi.org/10.21105/joss.02868)
[![Coverage Status](https://coveralls.io/repos/github/ivadomed/ivadomed/badge.svg?branch=master)](https://coveralls.io/github/ivadomed/ivadomed?branch=master)
[![test status](https://github.com/ivadomed/ivadomed/workflows/Run%20tests/badge.svg)](https://github.com/ivadomed/ivadomed/actions?query=workflow%3A%22Run+tests%22)
[![publish package](https://github.com/ivadomed/ivadomed/workflows/Publish%20Package/badge.svg)](https://github.com/ivadomed/ivadomed/actions?query=workflow%3A%22Publish+Package%22)
[![Documentation Status](https://readthedocs.org/projects/ivado-medical-imaging/badge/?version=latest)](https://ivadomed.org/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)
[![Twitter Follow](https://img.shields.io/twitter/follow/ivadomed.svg?style=social&label=Follow)](https://twitter.com/ivadomed)

`ivadomed` is an integrated framework for medical image analysis with deep learning.

The technical documentation is available [here](https://ivadomed.org).  The more detailed installation instruction is available [there](https://ivadomed.org/installation.html)

## Installation

``ivadomed`` requires Python >= 3.6 and < 3.10 as well as PyTorch == 1.8. We recommend working under a virtual environment, which could be set as follows:

```bash
python -m venv ivadomed_env
source ivadomed_env/bin/activate
```

### Install from release (recommended)

Install ``ivadomed`` and its requirements from `Pypi <https://pypi.org/project/ivadomed/>`__:

```bash
pip install --upgrade pip
pip install ivadomed
```

### Install from source

Bleeding-edge developments builds are available on the project's master branch on Github. Installation procedure is the following:

```bash
git clone https://github.com/neuropoly/ivadomed.git
cd ivadomed
pip install -e .
```

### Install from source

Make sure to install Torch1.8 following commands [here](https://pytorch.org/get-started/previous-versions/#v180) as pip is not able to auto infer GPU/CPU support on your behalf.
Again, the more comprehensive installation instruction is available [there](https://ivadomed.org/installation.html).

## Contributors
<p float="left">
  <img src="https://raw.githubusercontent.com/ivadomed/doc-figures/main/contributors/neuropoly_logo.png" height="80" />
  <img src="https://raw.githubusercontent.com/ivadomed/doc-figures/main/contributors/mila_logo.png" height="80" />
  <img src="https://raw.githubusercontent.com/ivadomed/doc-figures/main/contributors/ivado_logo.png" height="80" />
</p>

This project results from a collaboration between the [NeuroPoly Lab](https://www.neuro.polymtl.ca/)
and the [Mila](https://mila.quebec/en/). The main funding source is [IVADO](https://ivado.ca/en/).

[List of contributors](https://github.com/neuropoly/ivadomed/graphs/contributors)

## Consult our Wiki(https://github.com/ivadomed/ivadomed/wiki) here for more help

# Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leader responsible for enforcement at
jcohen [at] polymtl.ca.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
[https://www.contributor-covenant.org/version/2/0/code_of_conduct.html][v2.0].

Community Impact Guidelines were inspired by 
[Mozilla's code of conduct enforcement ladder][Mozilla CoC].

For answers to common questions about this code of conduct, see the FAQ at
[https://www.contributor-covenant.org/faq][FAQ]. Translations are available 
at [https://www.contributor-covenant.org/translations][translations].

[homepage]: https://www.contributor-covenant.org
[v2.0]: https://www.contributor-covenant.org/version/2/0/code_of_conduct.html
[Mozilla CoC]: https://github.com/mozilla/diversity
[FAQ]: https://www.contributor-covenant.org/faq
[translations]: https://www.contributor-covenant.org/translations
## v2.9.2 (2022-01-18)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.9.1...release)

**FEATURE**

- Implementation of Random Blur Augmentation.  [View pull request](https://github.com/ivadomed/ivadomed/pull/1034)
- Implementation of Random Bias Field Augmentation.  [View pull request](https://github.com/ivadomed/ivadomed/pull/1033)
- Implementation of Random Gamma Contrast Augmentation.  [View pull request](https://github.com/ivadomed/ivadomed/pull/1015)

**DEPENDENCIES**

- Unpin `tensorboard` to avoid conflict with downstream SCT requirements.  [View pull request](https://github.com/ivadomed/ivadomed/pull/1048)

**BUG**

- Rename prediction filenames: add class index and compat. for multi-rater.  [View pull request](https://github.com/ivadomed/ivadomed/pull/1043)
- Fix pixel size keyword in run_segment_command.  [View pull request](https://github.com/ivadomed/ivadomed/pull/1024)
- Replaced flip_axes with the correct bool element at index..  [View pull request](https://github.com/ivadomed/ivadomed/pull/1013)

**DOCUMENTATION**

- Add microscopy tutorial.  [View pull request](https://github.com/ivadomed/ivadomed/pull/1036)
- Removed one child-headings for clarity.  [View pull request](https://github.com/ivadomed/ivadomed/pull/1028)
- Typo fix for URL that is bricking the Colab link.  [View pull request](https://github.com/ivadomed/ivadomed/pull/1021)
- Experimental incorporation of tutorial jupyter notebooks open in Colab path.  [View pull request](https://github.com/ivadomed/ivadomed/pull/998)

## v2.9.1 (2021-12-13)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.9.0...release)

**ENHANCEMENT**

- Add forced indexation of "micr" datatype.  [View pull request](https://github.com/ivadomed/ivadomed/pull/995)
- Apply transforms on 2D patches.  [View pull request](https://github.com/ivadomed/ivadomed/pull/982)

**DOCUMENTATION**

- Update Tutorial 1/2/3 and readme.md to fix minor display issues.  [View pull request](https://github.com/ivadomed/ivadomed/pull/992)
- Update installation instruction to fit recent CUDA11 and torch 1.8+ push.  [View pull request](https://github.com/ivadomed/ivadomed/pull/969)

**REFACTORING**

- Fully Remove HeMIS model, Adaptive and h5py/HDF5.  [View pull request](https://github.com/ivadomed/ivadomed/pull/984)
- Use keywords for the rest of the files.  [View pull request](https://github.com/ivadomed/ivadomed/pull/946)


## v2.9.0 (2021-11-14)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.8.0...release)

**ENHANCEMENT**

- Make ivadomed be compatible with python 3.9 and PyTorch 1.8.  [View pull request](https://github.com/ivadomed/ivadomed/pull/819)

**DEPENDENCIES**

- Pin to CUDA-11.  [View pull request](https://github.com/ivadomed/ivadomed/pull/951)

**BUG FIXES**

- Pin PyParsing version to be compatible with pip 20.  [View pull request](https://github.com/ivadomed/ivadomed/pull/987)
- Fix pytest test_download_data_no_dataset_specified fail bug.  [View pull request](https://github.com/ivadomed/ivadomed/pull/968)
- Fix GeneralizedDiceLoss with `include_background=true` and `batch_size>1` .  [View pull request](https://github.com/ivadomed/ivadomed/pull/962)
- Fix undo_transforms in volume reconstruction.  [View pull request](https://github.com/ivadomed/ivadomed/pull/957)
- Fix undo_transforms in image reconstruction.  [View pull request](https://github.com/ivadomed/ivadomed/pull/956)
- add metadata to create_metadata_dict.  [View pull request](https://github.com/ivadomed/ivadomed/pull/954)
- Update scripts in `dev/prepare_data` to use new SCT config syntax (`.yml`).  [View pull request](https://github.com/ivadomed/ivadomed/pull/949)
- Fix config loading errors.  [View pull request](https://github.com/ivadomed/ivadomed/pull/944)
- Fix dropout_rate key in models.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/937)
- Add additional check for incorrect final_activation value.  [View pull request](https://github.com/ivadomed/ivadomed/pull/933)
- Make ivadomed be compatible with python3.9 and PyTorch 1.8.  [View pull request](https://github.com/ivadomed/ivadomed/pull/819)

**DOCUMENTATION**

- Minor modifications to the documentation for tutorial 3.  [View pull request](https://github.com/ivadomed/ivadomed/pull/988)
- Fix resample axis order in documentation.  [View pull request](https://github.com/ivadomed/ivadomed/pull/978)
- Update help.rst.  [View pull request](https://github.com/ivadomed/ivadomed/pull/967)
- Fixing issues in estimate uncertainty tutorial.  [View pull request](https://github.com/ivadomed/ivadomed/pull/936)
- Fix link to data file in ivadomed instructions.  [View pull request](https://github.com/ivadomed/ivadomed/pull/929)
- Fixes object detection path in cascaded architecture tutorial.  [View pull request](https://github.com/ivadomed/ivadomed/pull/922)
- Make ivadomed be compatible with python3.9 and PyTorch 1.8.  [View pull request](https://github.com/ivadomed/ivadomed/pull/819)

**REFACTORING**

- Fully Remove HeMIS model, Adaptive and h5py/HDF5.  [View pull request](https://github.com/ivadomed/ivadomed/pull/984)
- Fix path_output in automated training.  [View pull request](https://github.com/ivadomed/ivadomed/pull/914)
- Using keywords for ivadomed/scripts folder.  [View pull request](https://github.com/ivadomed/ivadomed/pull/934)
- Keywords refactoring Phase II: loader focus.  [View pull request](https://github.com/ivadomed/ivadomed/pull/909)
- Adopting pathllib for loader/bids_dataframe.  [View pull request](https://github.com/ivadomed/ivadomed/pull/947)
- Adopting pathlib for tests.  [View pull request](https://github.com/ivadomed/ivadomed/pull/901)
- Adopting pathlib training.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/897)
- Adopting pathlib for main.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/892)
- Adopting pathlib for loader/utils.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/879)

**TESTING**

- Fix pytest test_download_data_no_dataset_specified fail bug.  [View pull request](https://github.com/ivadomed/ivadomed/pull/968)

**CI**

- Update Sphinx dependency version and check RTD.org performance.  [View pull request](https://github.com/ivadomed/ivadomed/pull/974)
- Fix pytest problem.  [View pull request](https://github.com/ivadomed/ivadomed/pull/968)
- Update to GitHub Action to use `setup-python@v2`.  [View pull request](https://github.com/ivadomed/ivadomed/pull/959)
- Make ivadomed be compatible with python3.9 and PyTorch 1.8.  [View pull request](https://github.com/ivadomed/ivadomed/pull/819)

## v2.8.0 (2021-08-31)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.7.4...v.2.8.0)

**FEATURE**

- Add image reconstruction from 2D patches.  [View pull request](https://github.com/ivadomed/ivadomed/pull/782)
- Add sha256 for training data.  [View pull request](https://github.com/ivadomed/ivadomed/pull/760)

**CI**

- Exclude testing directory in coveralls.  [View pull request](https://github.com/ivadomed/ivadomed/pull/776)
- Improve current GitHub Action CI with multi OS support.  [View pull request](https://github.com/ivadomed/ivadomed/pull/757)

**BUG**

- Fix training_curve.py output.  [View pull request](https://github.com/ivadomed/ivadomed/pull/923)
- Fix inverted dimensions in microscopy pixelsize.  [View pull request](https://github.com/ivadomed/ivadomed/pull/916)
- Fix segment functions for models without pre-processing transforms.  [View pull request](https://github.com/ivadomed/ivadomed/pull/874)
- Fix microscopy ground-truth range of values.  [View pull request](https://github.com/ivadomed/ivadomed/pull/870)
- `utils.py`: Only raise ArgParseException for non-zero SystemExits.  [View pull request](https://github.com/ivadomed/ivadomed/pull/854)
- Remove `anaconda` from explicit dependencies.  [View pull request](https://github.com/ivadomed/ivadomed/pull/845)
- Fix multiclass evaluation bug.  [View pull request](https://github.com/ivadomed/ivadomed/pull/837)
- Fix last slice missing in testing bug.  [View pull request](https://github.com/ivadomed/ivadomed/pull/835)
- Skip all NumpyToTensor transformation for retrocompatibility.  [View pull request](https://github.com/ivadomed/ivadomed/pull/830)
- Remove all NumpyToTensor configs keys.  [View pull request](https://github.com/ivadomed/ivadomed/pull/826)
- Add missing "-r" flags to installation.rst.  [View pull request](https://github.com/ivadomed/ivadomed/pull/820)
- Call NumpyToTensor last.  [View pull request](https://github.com/ivadomed/ivadomed/pull/818)
- Fix bug in loader for multiple raters.  [View pull request](https://github.com/ivadomed/ivadomed/pull/806)
- Hot patch to address Inference issue #803.  [View pull request](https://github.com/ivadomed/ivadomed/pull/804)
- Add tmp and log file to gitignore.  [View pull request](https://github.com/ivadomed/ivadomed/pull/794)

**INSTALLATION**

- Remove `anaconda` from explicit dependencies.  [View pull request](https://github.com/ivadomed/ivadomed/pull/845)

**DOCUMENTATION**

- Fix neuropoly guidelines link in ivadomed contribution guidelines document.  [View pull request](https://github.com/ivadomed/ivadomed/pull/924)
- Change readme to point to the latest build version.  [View pull request](https://github.com/ivadomed/ivadomed/pull/875)
- Installation instruction steps explicity recommended for MacOS but not Linux.  [View pull request](https://github.com/ivadomed/ivadomed/pull/847)
- Clarified step 2 for pytorch/torchvision.  [View pull request](https://github.com/ivadomed/ivadomed/pull/842)
- Add missing "-r" flags to installation.rst.  [View pull request](https://github.com/ivadomed/ivadomed/pull/820)
- Update one class segmentation tutorial's output and segmentation image.  [View pull request](https://github.com/ivadomed/ivadomed/pull/779)
- Update documentation with the solution to failing test_adaptive.py on MacOS.  [View pull request](https://github.com/ivadomed/ivadomed/pull/771)
- Added link to JOSS paper.  [View pull request](https://github.com/ivadomed/ivadomed/pull/748)

**DEPENDENCIES**

- Remove `anaconda` from explicit dependencies.  [View pull request](https://github.com/ivadomed/ivadomed/pull/845)

**ENHANCEMENT**

- Fix training_curve.py output.  [View pull request](https://github.com/ivadomed/ivadomed/pull/923)
- Fix microscopy ground-truth range of values.  [View pull request](https://github.com/ivadomed/ivadomed/pull/870)
- Fix generate_sha_256 for joblib files.  [View pull request](https://github.com/ivadomed/ivadomed/pull/866)
- Add microscopy config file.  [View pull request](https://github.com/ivadomed/ivadomed/pull/850)
- Add the inference steps for PNG/TIF microscopy data.  [View pull request](https://github.com/ivadomed/ivadomed/pull/834)
- New loader: Load PNG/TIF/JPG microscopy files as Nibabel objects.  [View pull request](https://github.com/ivadomed/ivadomed/pull/813)
- Speed up IvadoMed Import Speed.  [View pull request](https://github.com/ivadomed/ivadomed/pull/793)
- Remove data dependencies from `if` statements in the `Decoder()` forward pass.  [View pull request](https://github.com/ivadomed/ivadomed/pull/752)

**TESTING**

- Unsilence test_rbg.  [View pull request](https://github.com/ivadomed/ivadomed/pull/832)
- Fix test_sampler.  [View pull request](https://github.com/ivadomed/ivadomed/pull/831)
- Fix bug in loader for multiple raters.  [View pull request](https://github.com/ivadomed/ivadomed/pull/806)
- Exclude testing directory in coveralls.  [View pull request](https://github.com/ivadomed/ivadomed/pull/776)
- Update documentation with the solution to failing test_adaptive.py on MacOS.  [View pull request](https://github.com/ivadomed/ivadomed/pull/771)
- Migrate test_segment_volume.py from unit_tests to functional_tests.  [View pull request](https://github.com/ivadomed/ivadomed/pull/767)
- Improve current GitHub Action CI with multi OS support.  [View pull request](https://github.com/ivadomed/ivadomed/pull/757)

**REFACTORING**

- Extract class SliceFilter, BalancedSample and SampleMetaData from loader.util.  [View pull request](https://github.com/ivadomed/ivadomed/pull/928)
- Extracted BidsDataFrame class outside of loader/utils.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/917)
- Initialize the adoption of centralized management of keywords via keywords.py (Phase I: compilation of all keywords).  [View pull request](https://github.com/ivadomed/ivadomed/pull/904)
- Fix empty list default parameter antipattern..  [View pull request](https://github.com/ivadomed/ivadomed/pull/903)
- Pathlib adoption for visualize.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/900)
- Pathlib adoption for utils.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/899)
- Pathlib adoption for uncertainty.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/898)
- Pathlib adoption for testing.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/896)
- Pathlib adoption for postprocessing.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/895)
- Pathlib adoption for models.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/894)
- Pathlib adoption for mixup.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/893)
- Pathlib adoption for inference.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/891)
- Pathlib adoption for evaluation.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/890)
- pathlib config_manager.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/889)
- Pathlib adoption for visualize_transform.  [View pull request](https://github.com/ivadomed/ivadomed/pull/888)
- Pathlib adoption for visualize_and_compare_testing_models.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/887)
- Pathlib adoption for script/training_curve.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/886)
- Pathlib adoption for script/prepare_dataset_vertibral_labeling.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/885)
- Pathlib adoption for extract_small_dataset.  [View pull request](https://github.com/ivadomed/ivadomed/pull/884)
- pathlib for download_data.  [View pull request](https://github.com/ivadomed/ivadomed/pull/883)
- pathlib for script/automate_training.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/881)
- pathlib change for object_detection/utils.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/880)
- Pathlib adoption for loader/segmentation_pair.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/878)
- Pathlib adoption for loader/film.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/877)
- Pathlib adoption for adaptative.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/876)
- Update config_bids.json following changes in microscopy BEP.  [View pull request](https://github.com/ivadomed/ivadomed/pull/838)
- Extracted Loader Classes into separate files.  [View pull request](https://github.com/ivadomed/ivadomed/pull/828)
- Refactor segment_volume to reduce complexity.  [View pull request](https://github.com/ivadomed/ivadomed/pull/791)
- Refactoring: BidsDataset __init__ reduce complexity .  [View pull request](https://github.com/ivadomed/ivadomed/pull/765)
- Refactoring: reduce complexity of BIDStoHDF5 _load_filenames.  [View pull request](https://github.com/ivadomed/ivadomed/pull/737)

## v2.7.4 (2021-03-15)

See `2.7.3`. We had to re-release because the GitHub Action didn't get triggered to push the release
to `PyPI` as it started as a draft. See here for more details:

[GitHub Actions Bug](https://github.community/t/workflow-set-for-on-release-not-triggering-not-showing-up/16286)


## v2.7.3 (2021-03-15)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.7.2...release)

**BUG**

 - Copy nibabel header when creating output prediction.  [View pull request](https://github.com/ivadomed/ivadomed/pull/714)
 - Dynamically write dataset_description.json file to suppress pybids warning.  [View pull request](https://github.com/ivadomed/ivadomed/pull/690)

**DOCUMENTATION**

 - Change archive links to repository links for pre-trained models.  [View pull request](https://github.com/ivadomed/ivadomed/pull/700)

**ENHANCEMENT**

 - New loader: Refactor BidsDataset classes.  [View pull request](https://github.com/ivadomed/ivadomed/pull/691)


## v2.7.2 (2021-02-19)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.7.1...v2.7.2)

**BUG**

 - Multiclass ignored during inference if n_input and n_output are different.  [View pull request](https://github.com/ivadomed/ivadomed/pull/688)
 - Merged participants.tsv file saving bug correction.  [View pull request](https://github.com/ivadomed/ivadomed/pull/684)
 - Make change_keys method from ConfigurationManager compatible with python3.8.  [View pull request](https://github.com/ivadomed/ivadomed/pull/681)

**DOCUMENTATION**

 - Add DOI JOSS.  [View pull request](https://github.com/ivadomed/ivadomed/pull/683)
 - Adding Zenodo DOI.  [View pull request](https://github.com/ivadomed/ivadomed/pull/677)

**ENHANCEMENT**

 - New loader: input from multiple BIDS datasets.  [View pull request](https://github.com/ivadomed/ivadomed/pull/687)
 - Add pre-commit hooks to limit file size to 500KB .  [View pull request](https://github.com/ivadomed/ivadomed/pull/682)
 - Shared weights for the two first FiLM generator layers.  [View pull request](https://github.com/ivadomed/ivadomed/pull/679)
 - Allow for non-dictionary hyperparameters in automate_training.  [View pull request](https://github.com/ivadomed/ivadomed/pull/661)

**FEATURE**

 - Enable the pipeline to run with inputs from multiple BIDS datasets.  [View pull request](https://github.com/ivadomed/ivadomed/pull/588)

## v2.7.1 (2021-02-09)
[View change](https://github.com/ivadomed/ivadomed/pull/676)

## v2.7.0 (2021-02-09)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.6.1...v2.7.0)

**BUG**

 - Fix structure wise uncertainty computation.  [View pull request](https://github.com/ivadomed/ivadomed/pull/664)
 - Fix bugs in plot film params.  [View pull request](https://github.com/ivadomed/ivadomed/pull/646)
 - Change condition to save FiLM parameters .  [View pull request](https://github.com/ivadomed/ivadomed/pull/645)
 - Fix store film params.  [View pull request](https://github.com/ivadomed/ivadomed/pull/642)
 - soft_gt param: only active after Data Augmentation.  [View pull request](https://github.com/ivadomed/ivadomed/pull/624)
 - AnimatedGIf import and documentation.  [View pull request](https://github.com/ivadomed/ivadomed/pull/623)
 - Fix pandas typecast issue in test_split_dataset.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/606)
 - Make sure test_HeMIS runs tests in order.  [View pull request](https://github.com/ivadomed/ivadomed/pull/602)
 - Fix loader/adaptative.py code with reading/writing HDF5 files.  [View pull request](https://github.com/ivadomed/ivadomed/pull/592)
 - Automate_training: fix bug for multiple parameters.  [View pull request](https://github.com/ivadomed/ivadomed/pull/586)
 - Load 2D GT slice as uint if not soft training.  [View pull request](https://github.com/ivadomed/ivadomed/pull/582)

**DOCUMENTATION**

 - Updated affiliations, Added Marie-Helene.  [View pull request](https://github.com/ivadomed/ivadomed/pull/674)
 - Fix missing dilate-gt.png.  [View pull request](https://github.com/ivadomed/ivadomed/pull/653)
 - Reformat configuration_file.rst for docs.  [View pull request](https://github.com/ivadomed/ivadomed/pull/650)
 - Add metavar to parser.  [View pull request](https://github.com/ivadomed/ivadomed/pull/641)
 - Add a README for the Sphinx docs.  [View pull request](https://github.com/ivadomed/ivadomed/pull/626)
 - Add documentation on packaged model format.  [View pull request](https://github.com/ivadomed/ivadomed/pull/625)
 - Add the Twitter badge.  [View pull request](https://github.com/ivadomed/ivadomed/pull/622)
 - Add new custom css rule for table in purpose section (#617).  [View pull request](https://github.com/ivadomed/ivadomed/pull/619)
 - Add DeepReg to comparison table.  [View pull request](https://github.com/ivadomed/ivadomed/pull/618)
 - Update PyTorch Ref.  [View pull request](https://github.com/ivadomed/ivadomed/pull/616)
 - Small clarifications and typos fixes in the Unet tutorial.  [View pull request](https://github.com/ivadomed/ivadomed/pull/610)
 - Added warning on installation to make sure proper Python version is installed.  [View pull request](https://github.com/ivadomed/ivadomed/pull/607)
 - Made the BIDS example more general for the audience.  [View pull request](https://github.com/ivadomed/ivadomed/pull/597)

**ENHANCEMENT**

 - Add new keys config manager.  [View pull request](https://github.com/ivadomed/ivadomed/pull/668)
 - Store FiLM parameters during testing instead of training.  [View pull request](https://github.com/ivadomed/ivadomed/pull/663)
 - Externalize command, log_directory, and bids_path fields from JSON config files to CLI.  [View pull request](https://github.com/ivadomed/ivadomed/pull/652)
 - New loader: BidsDataframe class.  [View pull request](https://github.com/ivadomed/ivadomed/pull/648)
 - version_info.log  added in the log directory.  [View pull request](https://github.com/ivadomed/ivadomed/pull/639)
 - Indicate folder created after running ivadomed_download_data.  [View pull request](https://github.com/ivadomed/ivadomed/pull/609)
 - Add explanation for Windows incompatibility in installation docs.  [View pull request](https://github.com/ivadomed/ivadomed/pull/605)
 - Specify Python version in setup.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/603)
 - Add new filter to SliceFilter class.  [View pull request](https://github.com/ivadomed/ivadomed/pull/594)
 - New loader: Adapt splitting methods.  [View pull request](https://github.com/ivadomed/ivadomed/pull/591)

**TESTING**

 - Add functional test for automate_training run_test flag.  [View pull request](https://github.com/ivadomed/ivadomed/pull/647)
 - Add test template files.  [View pull request](https://github.com/ivadomed/ivadomed/pull/638)
 - Remove the testing_data folder from ivadomed.  [View pull request](https://github.com/ivadomed/ivadomed/pull/631)
 - Bug in Coveralls release 3.0.0.  [View pull request](https://github.com/ivadomed/ivadomed/pull/628)
 - Add tests for create_bids_dataframe function.  [View pull request](https://github.com/ivadomed/ivadomed/pull/584)

**REFACTORING**

 - Reformat configuration_file.rst for docs.  [View pull request](https://github.com/ivadomed/ivadomed/pull/650)
 - New loader: BidsDataframe class.  [View pull request](https://github.com/ivadomed/ivadomed/pull/648)
 - Standardize the gpu ID argument.  [View pull request](https://github.com/ivadomed/ivadomed/pull/644)
 - Unit Test cleanup.  [View pull request](https://github.com/ivadomed/ivadomed/pull/636)
 - Remove test_script and ivado_functional_test files.  [View pull request](https://github.com/ivadomed/ivadomed/pull/634)

## v2.6.1 (2020-12-15)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.6.0...v2.6.1)

**BUG**

 - Fix missing attribute softmax.  [View pull request](https://github.com/ivadomed/ivadomed/pull/547)
 - Split_dataset: consider center_list when per_patient is used.  [View pull request](https://github.com/ivadomed/ivadomed/pull/537)

**DOCUMENTATION**

 - Make usage clearer.  [View pull request](https://github.com/ivadomed/ivadomed/pull/578)
 - Removing support for Python 3.9 (for now).  [View pull request](https://github.com/ivadomed/ivadomed/pull/562)
 - Updating comparison table after review.  [View pull request](https://github.com/ivadomed/ivadomed/pull/560)

**ENHANCEMENT**

 - Remove small for multiclass.  [View pull request](https://github.com/ivadomed/ivadomed/pull/570)
 - Save config file before training.  [View pull request](https://github.com/ivadomed/ivadomed/pull/569)
 - Apply bounding box safety factor in segment volume.  [View pull request](https://github.com/ivadomed/ivadomed/pull/549)
 - Multichannel support for convert_to_onnx script.  [View pull request](https://github.com/ivadomed/ivadomed/pull/544)

**FEATURE**

 - Select subjects for training based on metadata.  [View pull request](https://github.com/ivadomed/ivadomed/pull/534)

## v2.6.0 (2020-11-23)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.5.0...v2.6.0)

**BUG**

 - Make is_2d retrocompatibility.  [View pull request](https://github.com/ivadomed/ivadomed/pull/535)
 - Support multiclass if first class missing.  [View pull request](https://github.com/ivadomed/ivadomed/pull/522)

**DOCUMENTATION**

 - AdapWing 3D: fix comment.  [View pull request](https://github.com/ivadomed/ivadomed/pull/531)
 - paper.md: overview_title.png path.  [View pull request](https://github.com/ivadomed/ivadomed/pull/529)
 - paper.bib: correct typo.  [View pull request](https://github.com/ivadomed/ivadomed/pull/528)
 - Fix DOIs in paper.bib.  [View pull request](https://github.com/ivadomed/ivadomed/pull/527)
 - Redirect to DokuWiki/GitHub from the contributing guidelines.  [View pull request](https://github.com/ivadomed/ivadomed/pull/523)
 - Change path for images.  [View pull request](https://github.com/ivadomed/ivadomed/pull/521)

**ENHANCEMENT**

 - automate_training: add new parameter to change multiple params.  [View pull request](https://github.com/ivadomed/ivadomed/pull/533)
 - Softseg multiclass.  [View pull request](https://github.com/ivadomed/ivadomed/pull/530)
 - Multiclass and multichannel support for segment volume.  [View pull request](https://github.com/ivadomed/ivadomed/pull/524)

**FEATURE**

 - Create sample to balance metadata.  [View pull request](https://github.com/ivadomed/ivadomed/pull/503)

## v2.5.0 (2020-11-10)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.4.0...v2.5.0)

**BUG**

 - paper.md: Fixed broken link.  [View pull request](https://github.com/ivadomed/ivadomed/pull/517)
 - Change default value of config.json.  [View pull request](https://github.com/ivadomed/ivadomed/pull/514)

**DEPENDENCIES**

 - Requirements.txt: force onnxruntime version.  [View pull request](https://github.com/ivadomed/ivadomed/pull/505)
 - set h5py version in requirements.txt.  [View pull request](https://github.com/ivadomed/ivadomed/pull/500)

**DOCUMENTATION**

 - JOSS submission.  [View pull request](https://github.com/ivadomed/ivadomed/pull/502)

**ENHANCEMENT**

 - Some fixes to logging.  [View pull request](https://github.com/ivadomed/ivadomed/pull/509)

**FEATURE**

 - Training without test set.  [View pull request](https://github.com/ivadomed/ivadomed/pull/498)
 - FiLM for 3D Unet.  [View pull request](https://github.com/ivadomed/ivadomed/pull/491)

**REFACTORING**

 - Refactor utils.py.  [View pull request](https://github.com/ivadomed/ivadomed/pull/497)

## v2.4.0 (2020-10-27)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.3.1...v2.4.0)

**BUG**

 - Fix missing version.txt in wheels package.  [View pull request](https://github.com/ivadomed/ivadomed/pull/488)

**DOCUMENTATION**

 - Added reference to arXiv citation.  [View pull request](https://github.com/ivadomed/ivadomed/pull/485)
 - Documenting release workflow.  [View pull request](https://github.com/ivadomed/ivadomed/pull/483)

**ENHANCEMENT**

 - Option to override postprocessing in segment volume.  [View pull request](https://github.com/ivadomed/ivadomed/pull/486)
 - Configuration File Manager.  [View pull request](https://github.com/ivadomed/ivadomed/pull/484)


## v2.3.1 (2020-10-19)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.3.0...v2.3.1)

**BUG**

 - Version format.  [View pull request](https://github.com/ivadomed/ivadomed/pull/481)

## v2.3.0 (2020-10-19)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.2.1...v2.3.0)

**BUG**

 - Adapt all metrics to multiclass predictions.  [View pull request](https://github.com/ivadomed/ivadomed/pull/472)
 - fix run_test gpu assignation.  [View pull request](https://github.com/ivadomed/ivadomed/pull/453)

**DOCUMENTATION**

 - Improving documentation.  [View pull request](https://github.com/ivadomed/ivadomed/pull/477)
 - Tutorial fix.  [View pull request](https://github.com/ivadomed/ivadomed/pull/461)

**ENHANCEMENT**

 - Download data: Add models.  [View pull request](https://github.com/ivadomed/ivadomed/pull/476)
 - Refactoring: Changing print and exit to raise error.  [View pull request](https://github.com/ivadomed/ivadomed/pull/467)
 - Remove "eval" cmd.  [View pull request](https://github.com/ivadomed/ivadomed/pull/465)
 - Custom final activation.  [View pull request](https://github.com/ivadomed/ivadomed/pull/458)
 - Display version.  [View pull request](https://github.com/ivadomed/ivadomed/pull/456)

**FEATURE**

 - Use custom data for film.  [View pull request](https://github.com/ivadomed/ivadomed/pull/460)
 - Uncertainty as post-processing step.  [View pull request](https://github.com/ivadomed/ivadomed/pull/459)

## v2.2.1 (2020-09-22)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.2.0...v2.2.1)

**BUG**

 - Cover image path change on README.  [View pull request](https://github.com/ivadomed/ivadomed/pull/451)

## v2.2.0 (2020-09-22)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.1.0...v2.2.0)

**BUG**

 - Minor fixes prior release.  [View pull request](https://github.com/ivadomed/ivadomed/pull/449)

**DEPENDENCIES**

 - Modify scripts/training_curve.py to avoid tensorflow dependency.  [View pull request](https://github.com/ivadomed/ivadomed/pull/396)

**DOCUMENTATION**

 - Updating documentation.  [View pull request](https://github.com/ivadomed/ivadomed/pull/425)
 - Tutorial on uncertainty estimation.  [View pull request](https://github.com/ivadomed/ivadomed/pull/399)
 - Tutorial cascaded architecture.  [View pull request](https://github.com/ivadomed/ivadomed/pull/389)

**ENHANCEMENT**

 - Retrain model without resetting weights.  [View pull request](https://github.com/ivadomed/ivadomed/pull/447)
 - Normalized ReLU.  [View pull request](https://github.com/ivadomed/ivadomed/pull/384)
 - Create Ivadomed download function.  [View pull request](https://github.com/ivadomed/ivadomed/pull/379)

**FEATURE**

 - Evenly distribute subjects according to metadata.  [View pull request](https://github.com/ivadomed/ivadomed/pull/423)
 - Resume training.  [View pull request](https://github.com/ivadomed/ivadomed/pull/416)
 - Find optimal threshold with ROC analysis.  [View pull request](https://github.com/ivadomed/ivadomed/pull/383)
 - Generate GIF during training.  [View pull request](https://github.com/ivadomed/ivadomed/pull/374)
 - Add classifier model .  [View pull request](https://github.com/ivadomed/ivadomed/pull/278)

**TESTING**

 - Create coverage and improve testing.  [View pull request](https://github.com/ivadomed/ivadomed/pull/385)

## v2.1.0 (2020-07-21)
[View detailed changelog](https://github.com/ivadomed/ivadomed/compare/v2.0.2...v2.1.0)

**BUG**

 - Automate training seed.  [View pull request](https://github.com/ivadomed/ivadomed/pull/366)
 - Automate training bug.  [View pull request](https://github.com/ivadomed/ivadomed/pull/363)
 - Apply preprocessing after filter ROI.  [View pull request](https://github.com/ivadomed/ivadomed/pull/342)
 - Fix bug in automate training.  [View pull request](https://github.com/ivadomed/ivadomed/pull/339)
 - Transformations at test time: minor fixes.  [View pull request](https://github.com/ivadomed/ivadomed/pull/335)

**DOCUMENTATION**

 - Documentation: metric more formal defintion.  [View pull request](https://github.com/ivadomed/ivadomed/pull/357)
 - Fix few documentation issues, add content.  [View pull request](https://github.com/ivadomed/ivadomed/pull/341)
 - Soft training: minor fixes.  [View pull request](https://github.com/ivadomed/ivadomed/pull/334)
 - Tutorial 01: One class segmentation 2D Unet.  [View pull request](https://github.com/ivadomed/ivadomed/pull/309)

**ENHANCEMENT**

 - Split dataset with no test center specified.  [View pull request](https://github.com/ivadomed/ivadomed/pull/370)
 - showing time after training (begin/end/duration).  [View pull request](https://github.com/ivadomed/ivadomed/pull/365)
 - Optimize binarization.  [View pull request](https://github.com/ivadomed/ivadomed/pull/364)
 - Automate training improvement.  [View pull request](https://github.com/ivadomed/ivadomed/pull/362)
 - Simplify code when filtering ROI.  [View pull request](https://github.com/ivadomed/ivadomed/pull/361)
 - Scripts: Add entry points, modify doc display, and started to add github action testing.  [View pull request](https://github.com/ivadomed/ivadomed/pull/328)
# Testing

We have divided our tests into two categories, `functional_tests` and `unit_tests`. In each
folder, you will find a `t_utils` file with some helper functions, and a `t_template` file,
which provides a testing template.

## Running in GitHub

Checkout `ivadomed/.github/workflows/run_tests.yml` to see how tests are run on pull requests.

## Running Locally

1. Download the required dataset(s) using the `ivadomed` command line tools:
```
cd ivadomed  # root of the repo
ivadomed_download_data -d data_testing -o data_testing  # for unit tests
ivadomed_download_data -d data_functional_testing -o data_functional_testing  # for functional tests
```
2. To run all tests:
```
pytest
```
or, to run specific tests:
```
pytest testing/functional_tests/
pytest testing/unit_tests/
pytest testing/functional_tests/test_example.py
```

## Wiki

You can read more about our testing here: https://github.com/ivadomed/ivadomed/wiki/Tests
[Tutorial 1 2D Segmentation UNet ![Open Tutorial 1 on Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ivadomed/ivadomed/blob/master/testing/tutorials/tutorial_1_2d_segmentation_unet.ipynb)

[Tutorial 3 Uncertainty Estimation ![Open Tutorial 3 on Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ivadomed/ivadomed/blob/master/testing/tutorials/tutorial_3_uncertainty_estimation.ipynb)


<!-- Hi, and thank you for submitting a Pull Request! The checklist below is a brief summary of steps found in the NeuroPoly Contributing Guidelines, which can be found here: https://www.neuro.polymtl.ca/software/contributing. 
-->

## Checklist

#### GitHub

- [ ] I've given this PR a concise, self-descriptive, and meaningful title
- [ ] I've linked relevant issues in the PR body
- [ ] I've applied [the relevant labels](https://www.neuro.polymtl.ca/software/contributing#pr_labels) to this PR
- [ ] I've assigned a reviewer

<!-- For the title, please observe the following rules:
	- Provide a concise and self-descriptive title
	- Do not include the applicable issue number in the title, do it in the PR body
	- If the PR is not ready for review, convert it to a draft.
-->

#### PR contents

- [ ] I've consulted [ivadomed's internal developer documentation](https://github.com/ivadomed/ivadomed/wiki) to ensure my contribution is in line with any relevant design decisions
- [ ] I've added [relevant tests](https://github.com/ivadomed/ivadomed/wiki/tests) for my contribution
- [ ] I've updated the [relevant documentation](https://github.com/ivadomed/ivadomed/wiki/documentation) for my changes, including argparse descriptions, docstrings, and ReadTheDocs tutorial pages

## Description
<!-- describe what the PR is about. Explain the approach and possible drawbacks.It's ok to repeat some text from the related issue. -->

## Linked issues
<!-- If the PR fixes any issues, indicate it here with issue-closing keywords: e.g. Resolves #XX, Fixes #XX, Addresses #XX. Note that if you want multiple issues to be autoclosed on PR merge, you must use the issue-closing verb before each relevant issue: e.g. Resolves #1, Resolves #2 -->
---
name: Bug report
about: Report a bug or error to help improve ivadomed

---
<!-- Hi, and thank you for reporting an issue! Please take the time to first consider NeuroPoly's guidelines on issues (titles, descriptions, issue labels) before continuing:
https://www.neuro.polymtl.ca/software/contributing#opening_an_issue
-->


## Issue description

### Current behavior
<!-- Please provide a brief overview of the bug or error you've encountered.
-->


### Expected behavior
<!-- Please provide a description of what you expected to happen instead.
--> 


### Steps to reproduce
<!-- Please provide the exact commands that lead to the issue (and their resulting output), as well as any configuration files or data files you were working with, if possible. 
Note: When providing lengthy terminal output, please consider using GitHub's "details" tag to keep the report tidy. (See the example tag formatting below.)

<details>
<summary> An example summary message </summary>

```
Example terminal output.
```
</details>
-->


## Environment

### System description
<!-- Please describe the device and system you're working on, including the OS type and its version (e.g. Windows, macOS, Linux).
-->


### Installed packages
<!-- Please provide the output of `pip freeze` in the space between the two ``` blocks below.
-->

<details open>
<summary> Output of <code>pip freeze</code> </summary>

```

```

</details>

---
name: Feature request
about: Want improvements or new features in ivadomed?

---

<!-- Hi, and thank you for requesting a feature! Please take the time to first consider NeuroPoly's guidelines on issues (titles, descriptions, issue labels) before continuing:
https://www.neuro.polymtl.ca/software/contributing#opening_an_issue
-->

### Motivation for the feature
<!-- Is your feature request related to a problem? Please provide a clear and concise description of what the problem is. Ex. I'm always frustrated when [...]
-->

### Description of the feature
<!-- A clear and concise description of the feature you would like added, and any implementation details you've thought of.
-->

### Alternatives
<!-- A clear and concise description of any alternative solutions or features you've considered.
-->

### Additional context
<!-- Add any other context (references, screenshots...) about the feature request here.
-->
Development Scripts
===================

These scripts help the development and testing process.


Shortcut
--------

```{bash}
. dev/activate
```

will add the scripts your $PATH so they can be run from anywhere in the project.
(this is a work in progress; add more)
# Data preparation

These scripts prepare the data for training. It takes as input the [Spine Generic Public Database (Multi-Subject)](https://github.com/spine-generic/data-multi-subject) and outputs BIDS-compatible datasets with segmentation labels for each subject. More specifically, for each subject, the segmentation is run in one volume (T1w), then all volumes are registered to the T1w volume so that all volumes are in the same voxel space and the unique segmentation can be used across volumes.

## Dependencies

In its current state, this pipeline uses [SCT development version](https://github.com/neuropoly/spinalcordtoolbox#install-from-github-development). Once the pipeline is finalized, a stable version of SCT will be associated with this pipeline and indicated here. For now, please use the latest development version of SCT.

## How to run

#### Activate environment

See [README](../README.md)
~~~
source PATH_TO_YOUR_VENV/venv-ivadomed/bin/activate
~~~

#### Initial steps, check for folder integrity

- Copy the file `config_template.yml` and rename it as `config.yml`.
- Edit the file `config.yml` and modify the values according to your needs.
- Make sure input files are present:
~~~
sct_run_batch -script check_input_files.sh -config config.yml
~~~

#### Run first processing

Loop across subjects and run full processing:

~~~
sct_run_batch -script prepare_data.sh -config config.yml
~~~

#### Perform QC

##### Spinal cord segmentations

- Open qc/index.html
- Search only for "deepseg" QC entries (use "search" field)
- Take a screenshot of the browser when you spot a problem (wait for the segmentation to appear before taking the screenshot)
- If the data are of **very** bad quality, also take a screenshot (this time, wait for the segmentation to disappear)
- Copy all screenshots under qc_feedback/

##### Registration of MT scans

- Search for "register_multimodal"
- Take a screenshot of the browser when you spot a problem (wait for the segmentation to appear before taking the screenshot)
- If the data are of **very** bad quality, also take a screenshot (this time, wait for the segmentation to disappear)
- Copy all screenshots under qc_feedback/

##### Manually correct the segmentations

Check the following files under e.g. `result/sub-balgrist01/anat/tmp`:

| Image  | Segmentation  |
|:---|:---|
| sub-XX_acq-T1w_MTS_crop_r.nii.gz | sub-XX_acq-T1w_MTS_crop_r_seg.nii.gz|
| sub-XX_T1w_reg.nii.gz | sub-XX_T1w_reg_seg.nii.gz |
| sub-XX_T2w_reg.nii.gz | sub-XX_T2w_reg_seg.nii.gz |
| sub-XX_T2star_mean_reg.nii.gz | sub-XX_T2star_mean_reg_seg.nii.gz |

- Open the segmentation with `fsleyes`
- Manually correct it:
  - If the segmentation is leaking, remove the leak (use CMD+F to switch the overlay on/off)
  - If the segmentation exists in one slice but only consists of a few pixels, because the image quality is bad or because it is no more covering the cord (e.g. brainstem), remove all pixels in the current slice (better to have no segmentation than partial segmentation).
  - If the spinal cord is only partially visible (this can happen in T2star scans due to the registration), zero all pixels in the slice.
- Save with suffix `-manual`.
- Move to a folder named seg_manual/$FILENAME. E.g.: `~/data-multi-subject/derivatives/seg_manual/sub-amu01_acq-T1w_MTS_crop_r_seg-manual.nii.gz`

#### Exclude images

If some images are of unacceptable quality, they could be excluded from the final output dataset. List images to exclude in **config.yml** using the field `exclude-list`. 

#### Re-run processing (using manually-corrected segmentations)

Make sure to place your manually-corrected segmentations in the directory specified by `config.yml`, then re-run:

~~~
sct_run_batch -script prepare_data.sh -config config.yml
~~~

#### Copy files, final QC

Copy final files to anat/, copy json sidecars, move segmentations to derivatives/ and generate another QC:

~~~
sct_run_batch -script final_qc.sh config.yml
~~~

- Open the new QC: qc2/index.html
- Make sure that:
  - the final segmentation properly overlays on each contrast,
  - there is no missing slice (can happen for t2s data),
  - each contrast has sufficient image quality.
- If you spot any problem, take a screenshot of the browser and copy screenshots under qc2_feedback/

#### Clean temporary files

Once QC and manual correction is done, remove tmp/ folder:

~~~
sct_run_batch -script delete_tmp_files.sh -config config.yml
~~~
# Ivadomed Docs

## Installation

Assuming you have already installed the `ivadomed` package, you will also need to install:

```
pip install sphinx_rtd_theme
```

## Build

To create the html pages from the `.rst` files:

```
cd ivadomed/docs
make html
```

Check out the `Makefile` for more information.

## View

Under `docs/build/html`, open `index.html` in your browser to preview
the `Sphinx` documentation.
Contributing to ivadomed
========================

Thank you for your interest in contributing to ivadomed! This project uses the following pages to guide new contributions:

  * The `ivadomed GitHub repository <https://github.com/ivadomed/ivadomed>`_ is where the source code for the project is maintained, and where new contributions are submitted to.
  * The `NeuroPoly Contributing Guidelines <https://intranet.neuro.polymtl.ca/software-development/contributing>`_ provide instructions for development workflows, such as reporting issues or submitting pull requests.
  * The `ivadomed Developer Wiki <https://github.com/ivadomed/ivadomed/wiki>`_ acts as a knowledge base for documenting internal design decisions specific to the ivadomed codebase. It also contains step-by-step walkthroughs for common ivadomed maintainer tasks.Contributing to ivadomed
========================

Thank you for your interest in contributing to ivadomed! This project uses the following pages to guide new contributions:

  * The `ivadomed GitHub repository <https://github.com/ivadomed/ivadomed>`_ is where the source code for the project is maintained, and where new contributions are submitted to.
  * The `NeuroPoly Contributing Guidelines <https://intranet.neuro.polymtl.ca/software-development/contributing>`_ provide instructions for development workflows, such as reporting issues or submitting pull requests.
  * The `ivadomed Developer Wiki <https://github.com/ivadomed/ivadomed/wiki>`_ acts as a knowledge base for documenting internal design decisions specific to the ivadomed codebase. It also contains step-by-step walkthroughs for common ivadomed maintainer tasks.Help
====

If you need help using ``ivadomed``, please don't hesitate to
`post a question on our discussion forum <https://github.com/ivadomed/ivadomed/discussions>`__
and we will be happy to assist you.
.. |yes| raw:: html

   <style> .line {text-align:centers;} </style>
   <p style="color:green" align="center">	&#10004;</p>

.. |no| raw:: html

   <style> .line {text-align:centers;} </style>
   <p style="color:red" align="center">	&#10007;</p>

.. |cent| raw:: html

  <style> .line {text-align:center;} </style>


Purpose
=======

The purpose of the ``ivadomed`` project is to:

* Provide researchers with an open-source framework for training deep learning models for applications in medical imaging;

* Provide ready-to-use :doc:`pretrained_models` trained on multi-center data.

Comparison with other projects
------------------------------

We acknowledge the existence of projects with similar purposes. The table below compares some features across some
of the existing projects. This table was mostly based on the existing documentation for each project. We
understand that the field is rapidly evolving, and that this table might reflect the reality. If you notice
inconsistencies, please let us know by `opening an issue <https://github.com/ivadomed/ivadomed/issues>`_.

..
  If you wish to modify the csv table please modify https://docs.google.com/spreadsheets/d/1_MydnHnlOAuYzJ9QBCvPC9Jq2xUmPWI-XttTfcdtW2Y/edit#gid=0

.. csv-table::
   :file: comparison_other_projects_table.csv

(1): "BIDS" stands for the `Brain Imaging Data Structure <https://bids.neuroimaging.io/>`_, which is a convention initiated by the neuroimaging community to organize datasets (filenames, metadata, etc.). This facilitates the sharing of datasets and minimizes the burden of organizing datasets for training.

(2): **Class**: Classification | **Seg**: Segmentation | **Detect**: Detection | **Gen**: Generation | **Clust**: Clustering | **Reg**: Registration

Contributors
============

.. image:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/contributors/neuropoly_logo.png
  :height: 80
  :alt: Alternative text

.. image:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/contributors/mila_logo.png
  :height: 80
  :alt: Alternative text

This project results from a collaboration between the
`NeuroPoly Lab <https://www.neuro.polymtl.ca>`_ and `Mila <https://mila.quebec/en/>`_.

A list of contributors is available `here <https://github.com/neuropoly/ivadomed/graphs/contributors>`_.

Sponsors
--------
.. image:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/contributors/ivado_logo.png
  :height: 80
  :alt: Alternative text

If you wish to sponsor this project, please consider `donating <https://github.com/sponsors/neuropoly>`_.
Usage
=====

.. _usage:

Command line tools
------------------

New model can be generated using the command-line tool from the
terminal:

::

    ivadomed [command] -c path/to/config_file.json --path-data path/to/bids/data --path-output path/to/output/directory

``[command]`` represents the following choice of flags:

    ``--train``: train a model on a training/validation sub-datasets

    ``--test``: evaluate a trained model on a testing sub-dataset

    ``--segment``: segment a entire dataset using a trained model. Note that you may only specify one command flag at a time.

Note that the command CLI flag is optional and can be specified instead via the configuration file (see :ref:`configuration_file:Configuration File` ).
If not set via CLI, then you MUST specify this field in the configuration file.

``config_file.json`` is a configuration file, which parameters are
described in the :ref:`configuration_file:Configuration File`. This flag is *required*.

``path/to/bids/data`` is the location of the dataset. As discussed in :doc:`Data <data>`, the dataset
should conform to the BIDS standard. Modify the path so it points to the location of the downloaded dataset.

``path/to/output/directory`` is the folder name that will contain the output files (e.g., trained model, predictions, results)

Note that both path CLI flags are optional and can be specified instead via the configuration file.
If not set via CLI, then you MUST specify this field in the configuration file.

Please see section ``TUTORIALS`` to run this command on an example dataset.
Use cases
=========

Use case #1 - Spinal Cord Toolbox:
----------------------------------

`Spinal cord toolbox <http://spinalcordtoolbox.com/>`__ (SCT) is an open-source analysis software package for processing MRI data of the spinal cord `[De Leener et al. 2017] <https://doi.org/10.1016/j.neuroimage.2016.10.009>`__. `ivadomed` is SCT's backbone for the automated segmentation of the spinal cord, gray matter, tumors, and multiple sclerosis lesions, as well as for the labeling of intervertebral discs.

Use case 2 - Creation of anatomical template:
---------------------------------------------

`ivadomed` was used in the generation of several high-resolution anatomical MRI templates `[Calabrese et al. 2018 <https://doi.org/10.1038/s41598-018-24304-3>`__, `Gros et al. 2020] <https://github.com/sct-pipeline/exvivo-template>`__. To make anatomical templates, it is sometimes necessary to segment anatomical regions, such as the spinal cord white matter. When dealing with high resolution data, there may be several thousand 2D slices to segment, stressing the need for a fully-automated and robust solution. In these studies, only a handful of slices were manually-segmented and used to train a specific model. The model then predicted reliably and with high accuracy (Dice score > 90%) the delineation of anatomical structures for the thousands of remaining unlabeled slices.

Use case 3 - Tumor segmentation:
--------------------------------

`ivadomed` also proves to be useful in the context of clinical radiology routine REF, where clinicians need to segment tumors, edema, and cavity to establish prognosis and monitor the outcome. The framework is composed of a cascaded architecture that detects the spinal cord, crops the image around the region of interest, and segments the tumor (Figure herebelow). The resulting model can be applied to new data using only CPUs, which is more convenient in the clinical setting. The advanced features and architectures available in `ivadomed`, such as FiLM, were pivotal in obtaining encouraging results despite the difficulty of the task and the relatively low number of images.

.. figure:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/use_cases/lemay_2020.png
   :alt: Figure tumor segmentation

   Fully automatic spinal cord tumor segmentation framework. Step 1: The spinal cord is localized using a 3D U-Net and the image is cropped around the generated mask. Step 2: The spinal cord tumors are segmented.

   Figure tumor segmentation
Data
====

To facilitate the organization of data, ``ivadomed`` requires the data to be
organized according to the `Brain Imaging Data Structure (BIDS) <http://bids.neuroimaging.io/>`__ convention.
An example of this organization is shown below:

.. image:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/data/1920px-BIDS_Logo.png
    :alt: BIDS_Logo

::

    dataset/
    └── dataset_description.json
    └── participants.tsv
    └── sub-01  <--------------------- Folder enclosing data for subject 1
    └── sub-02
    └── sub-03
        └── anat
            └── sub-03_T1w.nii.gz  <-- MRI image in NIfTI format
            └── sub-03_T1w.json  <---- Metadata including image parameters, MRI vendor, etc.
            └── sub-03_T2w.nii.gz
            └── sub-03_T2w.json
    └── derivatives
        └── labels
            └── sub-03
                └── anat
                    └── sub-03_seg-tumor-manual.nii.gz  <-- Manually-corrected segmentation
                    └── sub-03_seg-tumor-manual.json  <---- Metadata including author who performed the labeling and date

.. note:: ``participants.tsv`` should, at least, include a column ``participant_id``, which is used when loading the dataset.

.. note:: For an exhaustive list of derivatives used in ``ivadomed``, please see our `wiki <https://github.com/ivadomed/ivadomed/wiki/repositories#derivatives>`_
.. _architectures:

Architectures
=============

The following architectures are availabled in ``ivadomed``.

:mod:`ResNet`
^^^^^^^^^^^^^

.. autoclass:: ivadomed.models.ResNet


:mod:`DenseNet`
^^^^^^^^^^^^^^^

.. autoclass:: ivadomed.models.DenseNet


:mod:`Unet`
^^^^^^^^^^^

.. autoclass:: ivadomed.models.Unet


:mod:`FiLMedUnet`
^^^^^^^^^^^^^^^^^

.. autoclass:: ivadomed.models.FiLMedUnet


:mod:`HeMISUnet`
^^^^^^^^^^^^^^^^

.. autoclass:: ivadomed.models.HeMISUnet


:mod:`Modified3DUNet`
^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: ivadomed.models.Modified3DUNet


:mod:`Countception`
^^^^^^^^^^^^^^^^^^^

.. autoclass:: ivadomed.models.Countception
Pre-trained models
==================

For convenience, the following pre-trained models are ready-to-use:

- `t2-tumor <https://github.com/ivadomed/t2_tumor>`_: Cord tumor segmentation model, trained on T2-weighted contrast.
- `t2star_sc <https://github.com/ivadomed/t2star_sc>`_: Spinal cord segmentation model, trained on T2-star contrast.
- `mice_uqueensland_gm <https://github.com/ivadomed/mice_uqueensland_gm>`_: Gray matter segmentation model on mouse MRI. Data from University of Queensland.
- `mice_uqueensland_sc <https://github.com/ivadomed/mice_uqueensland_sc>`_: Cord segmentation model on mouse MRI. Data from University of Queensland.
- `findcord_tumor <https://github.com/ivadomed/findcord_tumor>`_: Cord localisation model, trained on T2-weighted images with tumor.
- `model_find_disc_t1 <https://github.com/ivadomed/model_find_disc_t1>`_: Intervertebral disc detection model trained on T1-weighted images.
- `model_find_disc_t2 <https://github.com/ivadomed/model_find_disc_t2>`_: Intervertebral disc detection model trained on T2-weighted images.

Packaged model format
---------------------
Each folder contains a model (.pt or .onnx) with its corresponding configuration file (.json). The packaged model is
automatically generated during training. The folder containing the packaged model will be saved at the path specified by
the CLI flag ``--path-output`` or "path_output" in your config file. The packaged model, the configuration file, and the model file will
be named by the string specified by the key ``model_name`` in the configuration file.

.. code-block:: xml

   my_model
   ├── my_model.json
   └── my_model.onnx
Configuration File
==================

All parameters used for loading data, training and predicting are contained
within a single JSON configuration file. This section describes how to set up
this configuration file.

For convenience, here is an generic configuration file:
`config\_config.json <https://raw.githubusercontent.com/ivadomed/ivadomed/master/ivadomed/config/config.json>`__.

Below are other, more specific configuration files:

- `config\_classification.json <https://raw.githubusercontent.com/ivadomed/ivadomed/master/ivadomed/config/config_classification.json>`__: Trains a classification model.

- `config\_sctTesting.json <https://raw.githubusercontent.com/ivadomed/ivadomed/master/ivadomed/config/config_sctTesting.json>`__: Trains a 2D segmentation task with the U-Net architecture.

- `config\_spineGeHemis.json <https://raw.githubusercontent.com/ivadomed/ivadomed/master/ivadomed/config/config_spineGeHemis.json>`__: Trains a segmentation task with the HeMIS-UNet architecture.

- `config\_tumorSeg.json <https://raw.githubusercontent.com/ivadomed/ivadomed/master/ivadomed/config/config_tumorSeg.json>`__: Trains a segmentation task with a 3D U-Net architecture.


General Parameters
------------------


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "command",
        "description": "Run the specified command.",
        "type": "string",
        "options": {"train": "train a model on a training/validation sub-datasets",
                    "test": "evaluate a trained model on a testing sub-dataset",
                    "segment": "segment a entire dataset using a trained model"
        }
    }

.. code-block:: JSON

    {
        "command": "train"
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "gpu_ids",
        "description": "List of IDs of one or more GPUs to use.",
        "type": "list * integer"
    }

.. code-block:: JSON

    {
        "gpu_ids": [1,2,3]
    }

.. note::
    Currently only ``ivadomed_automate_training`` supports the use of more than one GPU.


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "log_directory",
        "description": "Folder name that will contain the output files (e.g., trained model,
            predictions, results).",
        "type": "string"
    }




.. code-block:: JSON

    {
        "path_output": "tmp/spineGeneric"
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "model_name",
        "description": "Folder name containing the trained model (ONNX format) and its configuration
            file, located within ``log_directory/``",
        "type": "string"
    }



.. code-block:: sh

    "log_directory/seg_gm_t2star/seg_gm_t2star.onnx"
    "log_directory/seg_gm_t2star/seg_gm_t2star.json"

When possible, the folder name will follow the following convention:
``task_(animal)_region_(contrast)`` with

.. code-block:: sh

   task = {seg, label, find}
   animal = {human, dog, cat, rat, mouse, ...}
   region = {sc, gm, csf, brainstem, ...}
   contrast = {t1, t2, t2star, dwi, ...}


.. code-block:: JSON

   {
       "model_name": "seg_gm_t2star"
   }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "debugging",
        "description": "Extended verbosity and intermediate outputs.",
        "type": "boolean"
    }



.. code-block:: JSON

    {
        "debugging": true
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "log_file",
        "description": "Name of the file to be logged to, located within ``log_directory/``",
        "type": "string"
    }



.. code-block:: JSON

    {
        "log_file": "log"
    }


Loader Parameters
-----------------

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "path_data",
        "description": "Path(s) of the BIDS folder(s).",
        "type": "list or str"
    }


.. code-block:: JSON

    {
        "loader_parameters": {
            "path_data": ["path/to/data_example_spinegeneric", "path/to/other_data_example"]
        }
    }

Alternatively:


.. code-block:: JSON

    {
        "loader_parameters": {
            "path_data": "path/to/data_example_spinegeneric"
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "bids_config",
        "description": "(Optional). Path of the custom BIDS configuration file for
            BIDS non-compliant modalities",
        "type": "string"
    }



.. code-block:: JSON

    {
        "loader_parameters": {
            "bids_config": "ivadomed/config/config_bids.json"
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "subject_selection",
        "description": "Used to specify a custom subject selection from a dataset.",
        "type": "dict",
        "options": {
            "n": {
                "description": "List containing the number subjects of each metadata.",
                "type": "list"
            },
            "metadata": {
                "$$description": [
                    "List of metadata used to select the subjects. Each metadata should be the name\n",
                    "of a column from the participants.tsv file."
                ],
                "type": "list"
            },
            "value": {
                "description": "List of metadata values of the subject to be selected.",
                "type": "list"
            }
        }
    }




.. code-block:: JSON

    {
        "loader_parameters": {
            "subject_selection": {
                "n": [5, 10],
                "metadata": ["disease", "disease"],
                "value": ["healthy", "ms"]
            }
        }
    }

In this example, a subdataset composed of 5 healthy subjects and 10 ms subjects will be selected
for training/testing.


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "target_suffix",
        "description": "Suffix list of the derivative file containing the ground-truth of interest.",
        "type": "list * string"
    }



.. code-block:: JSON

    {
        "loader_parameters": {
            "target_suffix": ["_seg-manual", "_lesion-manual"]
        }
    }

The length of this list controls the number of output channels of the model (i.e.
``out_channel``). If the list has a length greater than 1, then a
multi-class model will be trained. If a list of list(s) is input for a
training, (e.g. [[``"_seg-manual-rater1"``, ``"_seg-manual-rater2"``],
[``"_lesion-manual-rater1"``, ``"_lesion-manual-rater2"``]), then each
sublist is associated with one class but contains the annotations from
different experts: at each training iteration, one of these annotations
will be randomly chosen.


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "extensions",
        "description": "Used to specify a list of file extensions to be selected for
            training/testing. If not specified, then `.nii` and `.nii.gz` will be used by default.",
        "type": "list, string"
    }



.. code-block:: JSON

    {
        "loader_parameters": {
            "extensions": [".png"]
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "contrast_params",
        "type": "dict",
        "options": {
            "train_validation": {
                "type": "list, string",
                "$$description": [
                    "List of image contrasts (e.g. ``T1w``, ``T2w``) loaded for the training and\n",
                    "validation. If ``multichannel`` is ``true``, this list represents the different\n",
                    "channels of the input tensors (i.e. its length equals model's ``in_channel``).\n",
                    "Otherwise, the contrasts are mixed and the model has only one input channel\n",
                    "(i.e. model's ``in_channel=1``)"
                ]
            },
            "test": {
                "type": "list, string",
                "$$description": [
                    "List of image contrasts (e.g. ``T1w``, ``T2w``) loaded in the testing dataset.\n",
                    "Same comment as for ``train_validation`` regarding ``multichannel``."
                ]
            },
            "balance": {
                "type": "dict",
                "$$description": [
                    "Enables to weight the importance of specific channels (or contrasts) in the\n",
                    "dataset: e.g. ``{'T1w': 0.1}`` means that only 10% of the available ``T1w``\n",
                    "images will be included into the training/validation/test set. Please set\n",
                    "``multichannel`` to ``false`` if you are using this parameter."
                ]
            }
        }
    }



.. code-block:: JSON

    {
        "loader_parameters": {
            "contrast_params": {
                "training_validation": ["T1w", "T2w", "T2star"],
                "testing": ["T1w", "T2w", "T2star"],
                "balance": {}
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "multichannel",
        "description": "Indicated if more than a contrast (e.g. ``T1w`` and ``T2w``) is
            used by the model.",
        "type": "boolean"
    }

See details in both ``train_validation`` and ``test`` for the contrasts that are input.



.. code-block:: JSON

    {
        "loader_parameters": {
            "multichannel": false
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "slice_axis",
        "description": "Sets the slice orientation for 3D NIfTI files on which the model will be used.",
        "type": "string",
        "options": {"sagittal": "plane dividing body into left/right",
                    "coronal": "plane dividing body into front/back",
                    "axial": "plane dividing body into top/bottom"
        }
    }



.. code-block:: JSON

    {
        "loader_parameters": {
            "slice_axis": "sagittal"
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "slice_filter_params",
        "description": "Discard a slice from the dataset if it meets a condition, see
            below.",
        "type": "dict",
        "options": {
            "filter_empty_input": {
                "type": "boolean",
                "description": "Discard slices where all voxel
                   intensities are zeros."
            },
            "filter_empty_mask": {
                "type": "boolean",
                "description": "Discard slices where all voxel labels are zeros."
            },
            "filter_absent_class": {
                "type": "boolean",
                "$$description": [
                    "Discard slices where all voxel labels are zero for one or more classes\n",
                    "(this is most relevant for multi-class models that need GT for all classes at train time)."
                ]
            },
            "filter_classification": {
                "type": "boolean",
                "$$description": [
                    "Discard slices where all images fail a custom classifier filter. If used,\n",
                    "``classifier_path`` must also be specified, pointing to a saved PyTorch classifier."
                ]
            }
        }
    }


.. code-block:: JSON

    {
        "loader_parameters": {
            "slice_filter_params": {
                "filter_empty_mask": false,
                "filter_empty_input": true
            }
        }
    }

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "roi_params",
        "description": "Parameters for the region of interest",
        "type": "dict",
        "options": {
            "suffix": {
                "type": "string",
                "$$description": [
                    "Suffix of the derivative file containing the ROI used to crop\n",
                    "(e.g. ``_seg-manual``) with ``ROICrop`` as transform. Please use ``null`` if",
                    "you do not want to use an ROI to crop."
                ]
            },
            "slice_filter_roi": {
                "type": "int",
                "$$description": [
                    "If the ROI mask contains less than ``slice_filter_roi`` non-zero voxels\n",
                    "the slice will be discarded from the dataset. This feature helps with\n",
                    "noisy labels, e.g., if a slice contains only 2-3 labeled voxels, we do\n",
                    "not want to use these labels to crop the image. This parameter is only\n",
                    "considered when using ``ROICrop``."
                ]
            }
        }
    }



.. code-block:: JSON

    {
        "loader_parameters": {
            "roi_params": {
                "suffix": null,
                "slice_filter_roi": null
            }
        }
    }

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "soft_gt",
        "$$description": [
            "Indicates if a soft mask will be used as ground-truth to train\n",
            "and / or evaluate a model. In particular, the masks are not binarized\n",
            "after interpolations implied by preprocessing or data-augmentation operations."
        ],
        "type": "boolean"
    }

.. code-block:: JSON

    {
        "loader_parameters": {
            "soft_gt": true
        }
    }

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "is_input_dropout",
        "$$description": [
            "Indicates if input-level dropout should be applied during training.\n",
            "This option trains a model to be robust to missing modalities by setting \n",
            "to zero input channels (from 0 to all channels - 1). Always at least one \n",
            "channel will remain. If one or more modalities are already missing, they will \n",
            "be considered as dropped."
        ],
        "type": "boolean"
    }

.. code-block:: JSON

    {
        "loader_parameters": {
            "is_input_dropout": true
        }
    }



Split Dataset
-------------

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "fname_split",
        "$$description": [
            "File name of the log (`joblib <https://joblib.readthedocs.io/en/latest/>`__)\n",
            "that contains the list of training/validation/testing filenames. This file can later\n",
            "be used to re-train a model using the same data splitting scheme. If ``null``,\n",
            "a new splitting scheme is performed. If specified, the .joblib file data splitting scheme\n",
            "bypasses all the other split dataset parameters."
        ],
        "type": "string"
    }


.. code-block:: JSON

    {
        "split_dataset": {
            "fname_split": null
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "random_seed",
        "$$description": [
            "Seed used by the random number generator to split the dataset between\n",
            "training/validation/testing sets. The use of the same seed ensures the same split between\n",
            "the sub-datasets, which is useful for reproducibility."
        ],
        "type": "int"
    }

.. code-block:: JSON

    {
        "split_dataset": {
            "random_seed": 6
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "split_method",
        "$$description": [
            "Metadata contained in a BIDS tabular file on which the files are shuffled, then split\n",
            "between train/validation/test, according to ``train_fraction`` and ``test_fraction``.\n",
            "For example, ``participant_id`` from the ``participants.tsv`` file will shuffle all participants,\n",
            "then split between train/validation/test sets."
        ],
        "type": "string"
    }

.. code-block:: JSON

    {
        "split_dataset": {
            "split_method": "participant_id"
        }
    }

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "data_testing",
        "$$description": ["(Optional) Used to specify a custom metadata to only include in the testing dataset (not validation).\n",
            "For example, to not mix participants from different institutions between the train/validation set and the test set,\n",
			"use the column ``institution_id`` from ``participants.tsv`` in ``data_type``.\n"
		],
        "type": "dict",
        "options": {
            "data_type": {
                "$$description": [
					"Metadata to include in the testing dataset.\n",
					"If specified, the ``test_fraction`` is applied to this metadata."
                ],
                "type": "string"
            },
            "data_value": {
                "$$description": [
					"(Optional) List of metadata values from the ``data_type`` column to include in\n",
                    "the testing dataset. If specified, the testing set contains only files from the\n",
                    "``data_value`` list and the ``test_fraction`` is not used."
                ],
                "type": "list"
            }
        }
    }

.. code-block:: JSON

    {
        "split_dataset": {
            "data_testing": {"data_type": "institution_id", "data_value":[]}
        }
    }

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "balance",
        "$$description": [
            "Metadata contained in ``participants.tsv`` file with categorical values. Each category\n",
            "will be evenly distributed in the training, validation and testing datasets."
        ],
        "type": "string",
        "required": "false"
    }

.. code-block:: JSON

    {
        "split_dataset": {
            "balance": null
        }
    }

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "train_fraction",
        "description": "Fraction of the dataset used as training set.",
        "type": "float",
        "range": "[0, 1]"
    }

.. code-block:: JSON

    {
        "split_dataset": {
            "train_fraction": 0.6
        }
    }

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "test_fraction",
        "$$description": [
            "Fraction of the dataset used as testing set.\n"
        ],
        "type": "float",
        "range": "[0, 1]"
    }

.. code-block:: JSON

    {
        "split_dataset": {
            "test_fraction": 0.2
        }
    }


Training Parameters
-------------------

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "batch_size",
        "type": "int",
        "range": "(0, inf)"
    }

.. code-block:: JSON

    {
        "training_parameters": {
            "batch_size": 24
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "loss",
        "$$description": [
            "Metadata for the loss function. Other parameters that could be needed in the\n",
            "Loss function definition: see attributes of the Loss function of interest\n",
            "(e.g. ``'gamma': 0.5`` for ``FocalLoss``)."
        ],
        "type": "dict",
        "options": {
            "name": {
                "type": "string",
                "description": "Name of the loss function class. See :mod:`ivadomed.losses`"
            }
        }
    }

.. code-block:: JSON

    {
        "training_parameters": {
            "loss": {
                "name": "DiceLoss"
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "training_time",
        "$$description": [
            "Metadata for the loss function. Other parameters that could be needed in the\n",
            "Loss function definition: see attributes of the Loss function of interest\n",
            "(e.g. ``'gamma': 0.5`` for ``FocalLoss``)."
        ],
        "type": "dict",
        "options": {
            "num_epochs": {
                "type": "int",
                "range": "(0, inf)"
            },
            "early_stopping_epsilon": {
                "type": "float",
                "$$description": [
                    "If the validation loss difference during one epoch\n",
                    "(i.e. ``abs(validation_loss[n] - validation_loss[n-1]`` where n is the current epoch)\n",
                    "is inferior to this epsilon for ``early_stopping_patience`` consecutive epochs,\n",
                    "then training stops."
                ]
            },
            "early_stopping_patience": {
                "type": "int",
                "range": "(0, inf)",
                "$$description": [
                    "Number of epochs after which the training is stopped if the validation loss\n",
                    "improvement is smaller than ``early_stopping_epsilon``."
                ]
            }
        }
    }

.. code-block:: JSON

    {
        "training_parameters": {
            "training_time": {
                "num_epochs": 100,
                "early_stopping_patience": 50,
                "early_stopping_epsilon": 0.001
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "scheduler",
        "type": "dict",
        "options": {
            "initial_lr": {
                "type": "float",
                "description": "Initial learning rate."
            },
            "scheduler_lr": {
                "type": "dict",
                "options": {
                    "name": {
                        "type": "string",
                        "$$description": [
                            "One of ``CosineAnnealingLR``, ``CosineAnnealingWarmRestarts``\n",
                            "and ``CyclicLR``. Please find documentation `here <https://pytorch.org/docs/stable/optim.html>`__.\n",

                        ]
                    }
                },
                "description": "Other parameters depend on the scheduler of interest"
            }
        }
    }

.. code-block:: JSON

    {
        "training_parameters": {
            "scheduler": {
                "initial_lr": 0.001,
                "scheduler_lr": {
                    "name": "CosineAnnealingLR",
                    "max_lr": 1e-2,
                    "base_lr": 1e-5
                }
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "balance_samples",
        "description": "Balance labels in both the training and the validation datasets.",
        "type": "dict",
        "options": {
          "applied": {
              "type": "boolean",
              "description": "Indicates whether to use a balanced sampler or not."
          },
          "type": {
              "type": "string",
              "$$description": [
                "Indicates which metadata to use to balance the sampler.\n",
                "Choices: ``gt`` or  the name of a column from the ``participants.tsv`` file\n",
                "(i.e. subject-based metadata)"
              ]
          }
        }
     }

.. code-block:: JSON

    {
        "training_parameters": {
            "balance_samples": {
                "applied": false,
                "type": "gt"
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "mixup_alpha",
        "description": "Alpha parameter of the Beta distribution, see `original paper on
        the Mixup technique <https://arxiv.org/abs/1710.09412>`__.",
        "type": "float"
    }

.. code-block:: JSON

    {
        "training_parameters": {
            "mixup_alpha": null
        }
    }


.. jsonschema::

   {
       "$schema": "http://json-schema.org/draft-04/schema#",
       "title": "transfer_learning",
       "type": "dict",
       "options": {
           "retrain_model": {
               "type": "string",
               "$$description": [
                   "Filename of the pretrained model (``path/to/pretrained-model``). If ``null``,\n",
                   "no transfer learning is performed and the network is trained from scratch."
               ]
           },
           "retrain_fraction": {
               "type": "float",
               "range": "[0, 1]",
               "$$description": [
                   "Controls the fraction of the pre-trained model that will be fine-tuned. For\n",
                   "instance, if set to 0.5, the second half of the model will be fine-tuned while\n",
                   "the first layers will be frozen."
               ]
           },
           "reset": {
               "type": "boolean",
               "description": "If true, the weights of the layers that are not frozen
                  are reset. If false, they are kept as loaded."
           }
       }
   }

.. code-block:: JSON

    {
        "training_parameters": {
            "transfer_learning": {
                "retrain_model": null,
                "retrain_fraction": 1.0,
                "reset": true
            }
        }
   }


Architecture
------------

Architectures for both segmentation and classification are available and
described in the :ref:`architectures` section. If the selected architecture is listed in the
`loader <https://github.com/ivadomed/ivadomed/blob/lr/fixing_documentation/ivadomed/loader/loader.py>`__ file, a
classification (not segmentation) task is run. In the case of a
classification task, the ground truth will correspond to a single label
value extracted from ``target``, instead being an array (the latter
being used for the segmentation task).


.. jsonschema::

   {
       "$schema": "http://json-schema.org/draft-04/schema#",
       "title": "default_model",
       "required": "true",
       "type": "dict",
       "$$description": [
           "Define the default model (``Unet``) and mandatory parameters that are common to all\n",
           "available :ref:`architectures`. For custom architectures (see below), the default\n",
           "parameters are merged with the parameters that are specific to the tailored architecture."
       ],
       "options": {
           "name": {
               "type": "string",
               "description": "Default: ``Unet``"
           },
           "dropout_rate": {
               "type": "float",
               "description": "Default: ``0.3``"
           },
           "bn_momentum": {
               "type": "float",
               "$$description": [
                    "Defines the importance of the running average: (1 - `bn_momentum`). A large running\n",
                    "average factor will lead to a slow and smooth learning.\n",
                    "See `PyTorch's BatchNorm classes for more details. <https://pytorch.org/docs/stable/generated/torch.nn.BatchNorm2d.html>`__ for more details. Default: ``0.1``\n"
               ]

           },
           "depth": {
               "type": "int",
               "range": "(0, inf)",
               "description": "Number of down-sampling operations. Default: ``3``"
           },
           "final_activation": {
               "type": "string",
               "required": "false",
               "$$description": [
                   "Final activation layer. Options: ``sigmoid`` (default), ``relu`` (normalized ReLU), or ``softmax``."
               ]
           },
           "length_2D": {
                "type": "[int, int]",
                "description": "(Optional) Size of the 2D patches used as model's input tensors.",
                "required": "false"
            },
            "stride_2D": {
                "type": "[int, int]",
                "$$description": [
                    "(Optional) Strictly positive integers: Pixels' shift over the input matrix to create 2D patches.\n",
                    "Ex: Stride of [1, 2] will cause a patch translation of 1 pixel in the 1st dimension and 2 pixels in\n",
                    "the 2nd dimension at every iteration until the whole input matrix is covered."
                ],
                "required": "false"
            },
           "is_2d": {
               "type": "boolean",
               "$$description": [
                   "Indicates if the model is 2D, if not the model is 3D. If ``is_2d`` is ``False``, then parameters\n",
                   "``length_3D`` and ``stride_3D`` for 3D loader need to be specified (see :ref:`Modified3DUNet <Modified3DUNet>`)."
               ]
           }
       }
   }


.. code-block:: JSON

    {
        "default_model": {
            "name": "Unet",
            "dropout_rate": 0.3,
            "bn_momentum": 0.1,
            "depth": 3,
            "final_activation": "sigmoid"
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "FiLMedUnet",
        "type": "dict",
        "required": "false",
        "options": {
            "applied": {
                "type": "boolean",
                "description": "Set to ``true`` to use this model."
            },
            "metadata": {
                "type": "string",
                "options": {
                    "mri_params": {
                        "$$description": [
                            "Vectors of ``[FlipAngle, EchoTime, RepetitionTime, Manufacturer]``\n",
                            "(defined in the json of each image) are input to the FiLM generator."
                        ]
                    },
                    "contrast": "Image contrasts (according to ``config/contrast_dct.json``) are input to the FiLM generator."
               },
               "$$description": [
                   "Choice between ``mri_params``, ``contrasts`` (i.e. image-based metadata) or the\n",
                   "name of a column from the participants.tsv file (i.e. subject-based metadata)."
               ]
           }
       }
   }

.. code-block:: JSON

    {
        "FiLMedUnet": {
            "applied": false,
            "metadata": "contrasts",
            "film_layers": [0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        }
    }


.. jsonschema::


  	{
		"$schema": "http://json-schema.org/draft-04/schema#",
      	"title": "HeMISUnet",
      	"type": "dict",
      	"required": "false",
      	"options": {
			"applied": {
				"type": "boolean",
              	"description": "Set to ``true`` to use this model."
			},
          	"missing_probability": {
				"type": "float",
                "range": "[0, 1]",
                "$$description": [
                    "Initial probability of missing image contrasts as model's input\n",
                    "(e.g. 0.25 results in a quarter of the image contrasts, i.e. channels, that\n",
                    "will not be sent to the model for training)."
                ]
            },
            "missing_probability_growth": {
                "type": "float",
                "$$description": [
                    "Controls missing probability growth at each epoch: at each epoch, the\n",
                    "``missing_probability`` is modified with the exponent ``missing_probability_growth``."
                ]
            }
         }
      }

.. code-block:: JSON

    {
        "HeMISUnet": {
            "applied": true,
            "missing_probability": 0.00001,
            "missing_probability_growth": 0.9,
            "contrasts": ["T1w", "T2w"],
            "ram": true,
            "path_hdf5": "/path/to/HeMIS.hdf5",
            "csv_path": "/path/to/HeMIS.csv",
            "target_lst": ["T2w"],
            "roi_lst": null
        }
    }

.. _Modified3DUNet:

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "Modified3DUNet",
        "type": "dict",
        "required": "false",
        "options": {
            "length_3D": {
                "type": "[int, int, int]",
                "description": "Size of the 3D patches used as model's input tensors."
            },
            "stride_3D": {
                "type": "[int, int, int]",
                "$$description": [
                    "Voxels' shift over the input matrix to create patches. Ex: Stride of [1, 2, 3]\n",
                    "will cause a patch translation of 1 voxel in the 1st dimension, 2 voxels in\n",
                    "the 2nd dimension and 3 voxels in the 3rd dimension at every iteration until\n",
                    "the whole input matrix is covered."
                ]
            },
            "attention_unet": {
                "type": "boolean",
                "description": "Use attention gates in the Unet's decoder.",
                "required": "false"
            },
            "n_filters": {
                "type": "int",
                "$$description": [
                    "Number of filters in the first convolution of the UNet.\n",
                    "This number of filters will be doubled at each convolution."
                ],
                "required": "false"
            }
       }
   }

.. code-block:: JSON

    {
        "Modified3DUNet": {
            "applied": false,
            "length_3D": [128, 128, 16],
            "stride_3D": [128, 128, 16],
            "attention": false,
            "n_filters": 8
        }
    }


Cascaded Architecture Features
------------------------------

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "object_detection_params",
        "type": "dict",
        "required": "false",
        "options": {
            "object_detection_path": {
                "type": "string",
                "$$description": [
                    "Path to object detection model and the configuration file. The folder,\n",
                    "configuration file, and model need to have the same name\n",
                    "(e.g. ``findcord_tumor/``, ``findcord_tumor/findcord_tumor.json``, and\n",
                    "``findcord_tumor/findcord_tumor.onnx``, respectively). The model's prediction\n",
                    "will be used to generate bounding boxes."
                ]
            },
            "safety_factor": {
                "type": "[int, int, int]",
                "$$description": [
                    "List of length 3 containing the factors to multiply each dimension of the\n",
                    "bounding box. Ex: If the original bounding box has a size of 10x20x30 with\n",
                    "a safety factor of [1.5, 1.5, 1.5], the final dimensions of the bounding box\n",
                    "will be 15x30x45 with an unchanged center."
                ]
            }
       }
   }

.. code-block:: JSON

    {
        "object_detection_params": {
            "object_detection_path": null,
            "safety_factor": [1.0, 1.0, 1.0]
        }
    }


Transformations
---------------

Transformations applied during data augmentation. Transformations are sorted in the order they are applied to the image samples. For each transformation, the following parameters are customizable: 

- ``applied_to``: list between ``"im", "gt", "roi"``. If not specified, then the transformation is applied to all loaded samples. Otherwise, only applied to the specified types: Example: ``["gt"]`` implies that this transformation is only applied to the ground-truth data.
- ``dataset_type``: list between ``"training", "validation", "testing"``. If not specified, then the transformation is applied to the three sub-datasets. Otherwise, only applied to the specified subdatasets. Example: ``["testing"]`` implies that this transformation is only applied to the testing sub-dataset.


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "NumpyToTensor",
        "type": "dict",
        "description": "Converts nd array to tensor object."
    }

.. code-block:: JSON

    {
        "transformation": {
            "NumpyToTensor": {
                "applied_to": ["im", "gt"]
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "CenterCrop",
        "type": "dict",
        "options": {
            "size": {
                "type": "list, int"
            },
            "applied_to": {
                "type": "list, string"
            }
        }
    }

.. code-block:: JSON

    {
        "transformation": {
            "CenterCrop": {
                "applied_to": ["im", "gt"],
                "size":  [512, 256, 16]
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "ROICrop",
        "type": "dict",
        "options": {
            "size": {
                "type": "list, int"
            },
            "applied_to": {
                "type": "list, string"
            }
        }
    }

.. code-block:: JSON

    {
        "transformation": {
            "ROICrop": {
                "size": [48, 48]
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "NormalizeInstance",
        "type": "dict",
        "$$description": [
            "Normalize a tensor or an array image with mean and standard deviation estimated from\n",
            "the sample itself."
        ]
    }

.. code-block:: JSON

    {
        "transformation": {
            "NormalizeInstance": {}
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "RandomAffine",
        "type": "dict",
        "options": {
            "degrees": {
                "type": "float or tuple of float",
                "range": "(0, inf)",
                "$$description": [
                    "Positive float or list (or tuple) of length two. Angles in degrees. If only\n",
                    "a float is provided, then rotation angle is selected within the range\n",
                    "[-degrees, degrees]. Otherwise, the tuple defines this range."
                ]
            },
            "translate": {
                "type": "list, float",
                "range": "[0, 1]",
                "$$description": [
                    "Length 2 or 3 depending on the sample shape (2D or 3D). Defines\n",
                    "the maximum range of translation along each axis."
                ]
            },
            "scale": {
                "type": "list, float",
                "range": "[0, 1]",
                "$$description": [
                    "Length 2 or 3 depending on the sample shape (2D or 3D). Defines\n",
                    "the maximum range of scaling along each axis."
                ]
            }
        }
    }

.. code-block:: JSON

    {
        "transformation": {
            "RandomAffine": {
                "translate": [0.03, 0.03],
                "applied_to": ["im"],
                "dataset_type": ["training"],
                "scale": [0.1, 0.5],
                "degrees": 180
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "RandomShiftIntensity",
        "type": "dict",
        "options": {
            "shift_range": {
                "type": "[float, float]",
                "description": "Range from which the offset applied is randomly selected."
            }
        }
    }

.. code-block:: JSON

    {
        "transformation": {
            "RandomShiftIntensity": {
                "shift_range": [28.0, 30.0]
            }
        }
     }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "ElasticTransform",
        "type": "dict",
        "$$description": [
            "Applies elastic transformation. See also:\n",
            "`Best practices for convolutional neural networks
             applied to visual document analysis <http://cognitivemedium.com/assets/rmnist/Simard.pdf>`__."
        ],
        "options": {
            "alpha_range": {
                "type": "(float, float)",
                "description": "Deformation coefficient."
            },
            "sigma_range": {
                "type": "(float, float)",
                "description": "Standard deviation."
            },
            "p": {
                "type": "float"
            }
        }
    }

.. code-block:: JSON

    {
        "transformation": {
            "ElasticTransform": {
                "alpha_range": [28.0, 30.0],
                "sigma_range":  [3.5, 4.5],
                "p": 0.1,
                "applied_to": ["im", "gt"],
                "dataset_type": ["training"]
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "Resample",
        "type": "dict",
        "options": {
            "hspace": {
                "type": "float",
                "range": "[0, 1]",
                "description": "Resolution along the first axis, in mm."
            },
            "wspace": {
                "type": "float",
                "range": "[0, 1]",
                "description": "Resolution along the second axis, in mm."
            },
            "dspace": {
                "type": "float",
                "range": "[0, 1]",
                "description": "Resolution along the third axis, in mm."
            }
        }
    }

.. code-block:: JSON

    {
        "transformation": {
            "Resample": {
                "hspace": 0.75,
                "wspace": 0.75,
                "dspace": 1
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "AdditiveGaussianNoise",
        "type": "dict",
        "options": {
            "mean": {
                "type": "float",
                "description": "Mean of Gaussian noise."
            },
            "std": {
                "type": "float",
                "description": "Standard deviation of Gaussian noise."
            }
        }
    }

.. code-block:: JSON

    {
        "transformation": {
            "AdditiveGaussianNoise": {
                "mean": 0.0,
                "std": 0.02
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "DilateGT",
        "type": "dict",
        "options": {
            "dilation_factor": {
                "type": "float",
                "$$description": [
                    "Controls the number of iterations of ground-truth dilation depending on\n",
                    "the size of each individual lesion, data augmentation of the training set.\n",
                    "Use ``0`` to disable."
                ]
            }
        }
    }

.. code-block:: JSON

    {
        "transformation": {
            "DilateGT": {
                "dilation_factor": 0
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "HistogramClipping",
        "description": "Perform intensity clipping based on percentiles.",
        "type": "dict",
        "options": {
            "min_percentile": {
                "type": "float",
                "range": "[0, 100]"
            },
            "max_percentile": {
                "type": "float",
                "range": "[0, 100]"
            }
        }
    }

.. code-block:: JSON

    {
        "transformation": {
            "HistogramClipping": {
                "min_percentile": 50,
                "max_percentile": 75
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "Clahe",
        "type": "dict",
        "options": {
            "clip_limit": {
                "type": "float"
            },
            "kernel_size": {
                "type": "list, int",
                "$$description": [
                    "Defines the shape of contextual regions used in the algorithm.\n",
                    "List length = dimension, i.e. 2D or 3D"
                ]
            }
        }
    }

.. code-block:: JSON

    {
        "transformation": {
            "Clahe": {
                "clip_limit": 0.5,
                "kernel_size": [8, 8]
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "RandomReverse",
        "type": "dict",
        "description": "Make a randomized symmetric inversion of the different values of each dimensions."
    }

.. code-block:: JSON

    {
        "transformation": {
            "RandomReverse": {
                "applied_to": ["im"]
            }
        }
    }

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "RandomGamma",
        "type": "dict",
        "$$description": [
            "Randomly changes the contrast of an image by gamma exponential."
        ],
        "options": {
            "log_gamma_range": {
                "type": "(float, float)",
                "description": "Log gamma range for changing contrast."
            },
            "p": {
                "type": "float"
            }
        }
    }

.. code-block:: JSON

    {
        "transformation": {
            "RandomGamma": {
                "log_gamma_range": [-3.0, 3.0],
                "p": 0.5,
                "applied_to": ["im"],
                "dataset_type": ["training"]
            }
        }
    }

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "RandomBiasField",
        "type": "dict",
        "$$description": [
            "Applies a random MRI bias field artifact to the image via `torchio.RandomBiasField()`"
        ],
        "options": {
            "coefficients": {
                "type": "float",
                "description": "Maximum magnitude of polynomial coefficients."
            },
            "order": {
                "type": "int",
                "description": "Order of the basis polynomial functions."
            },
            "p": {
                "type": "float"
            }
        }
    }

.. code-block:: JSON

    {
        "transformation": {
            "RandomBiasField": {
                "coefficients": 0.5,
                "order": 3,
                "p": 0.5,
                "applied_to": ["im"],
                "dataset_type": ["training"]
            }
        }
    }

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "RandomBlur",
        "type": "dict",
        "$$description": [
            "Applies a random blur to the image."
        ],
        "options": {
            "sigma_range": {
                "type": "(float, float)",
                "description": "Standard deviation range for the gaussian filter."
            },
            "p": {
                "type": "float"
            }
        }
    }

.. code-block:: JSON

    {
        "transformation": {
            "RandomBlur": {
                "sigma_range": [0.0, 2.0],
                "p": 0.5,
                "applied_to": ["im"],
                "dataset_type": ["training"]
            }
        }
    }

.. _Uncertainty:

Uncertainty
-----------

Uncertainty computation is performed if ``n_it>0`` and at least
``epistemic`` or ``aleatoric`` is ``true``. Note: both ``epistemic`` and
``aleatoric`` can be ``true``.


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "epistemic",
        "type": "boolean",
        "description": "Model-based uncertainty with `Monte Carlo Dropout <https://arxiv.org/abs/1506.02142>`__."
    }

.. code-block:: JSON

    {
        "uncertainty": {
            "epistemic": true
        }
    }

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "aleatoric",
        "type": "boolean",
        "description": "Image-based uncertainty with `test-time augmentation <https://doi.org/10.1016/j.neucom.2019.01.103>`__."
    }

.. code-block:: JSON

    {
        "uncertainty": {
            "aleatoric": true
        }
    }

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "n_it",
        "type": "int",
        "description": "Number of Monte Carlo iterations. Set to 0 for no uncertainty computation."
    }

.. code-block:: JSON

    {
        "uncertainty": {
            "n_it": 2
        }
    }


Postprocessing
--------------

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "binarize_prediction",
        "type": "dict",
        "options": {
            "thr": {
                "type": "float",
                "range": "[0, 1]",
                "$$description": [
                    "Threshold. To use soft predictions (i.e. no binarisation, float between 0 and 1)\n",
                    "for metric computation, indicate -1."
                ]
            }
        },
        "$$description": [
            "Binarizes predictions according to the given threshold ``thr``. Predictions below the\n",
            "threshold become 0, and predictions above or equal to threshold become 1."
        ]
    }



.. code-block:: JSON

    {
        "postprocessing": {
            "binarize_prediction": {
                "thr": 0.1
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "binarize_maxpooling",
        "type": "dict",
        "$$description": [
            "Binarize by setting to 1 the voxel having the maximum prediction across all classes.\n",
            "Useful for multiclass models. No parameters required (i.e., {})."
        ]
    }



.. code-block:: JSON

    {
        "postprocessing": {
            "binarize_maxpooling": {}
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "fill_holes",
        "type": "dict",
        "description": "Fill holes in the predictions. No parameters required (i.e., {})."
    }



.. code-block:: JSON

    {
        "postprocessing": {
            "fill_holes": {}
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "keep_largest",
        "type": "dict",
        "$$description": [
            "Keeps only the largest connected object in prediction. Only nearest neighbors are\n",
            "connected to the center, diagonally-connected elements are not considered neighbors.\n",
            "No parameters required (i.e., {})"
        ]
    }



.. code-block:: JSON

    {
        "postprocessing": {
            "keep_largest": {}
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "remove_noise",
        "type": "dict",
        "options": {
            "thr": {
                "type": "float",
                "range": "[0, 1]",
                "description": "Threshold. Threshold set to ``-1`` will not apply this postprocessing step."
            }
        },
        "description": "Sets to zero prediction values strictly below the given threshold ``thr``."
    }



.. code-block:: JSON

    {
        "postprocessing": {
            "remove_noise": {
                "thr": 0.1
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "remove_small",
        "type": "dict",
        "$$description": [
            "Remove small objects from the prediction. An object is defined as a group of connected\n",
            "voxels. Only nearest neighbors are connected to the center, diagonally-connected\n",
            "elements are not considered neighbors."
        ],
        "options": {
            "thr": {
                "type": "int or list",
                "$$description": [
                    "Minimal object size. If a list of thresholds is chosen, the length should\n",
                    "match the number of predicted classes."
                ]
            },
            "unit": {
                "type": "string",
                "$$description": [
                    "Either `vox` for voxels or `mm3`. Indicates the unit used to define the\n",
                    "minimal object size."
                ]
            }
        }
    }



.. code-block:: JSON

    {
        "postprocessing": {
            "remove_small": {
                "unit": "vox",
                "thr": 3
            }
        }
    }

.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "threshold_uncertainty",
        "type": "dict",
        "$$description": [
            "Removes the most uncertain predictions (set to 0) according to a threshold ``thr``\n",
            "using the uncertainty file with the suffix ``suffix``. To apply this method,\n",
            "uncertainty needs to be evaluated on the predictions with the :ref:`uncertainty <Uncertainty>` parameter."
        ],
        "options": {
            "thr": {
                "type": "float",
                "range": "[0, 1]",
                "$$description": [
                    "Threshold. Threshold set to ``-1`` will not apply this postprocessing step."
                ]
            },
            "suffix": {
                "type": "string",
                "$$description": [
                    "Indicates the suffix of an uncertainty file. Choices: ``_unc-vox.nii.gz`` for\n",
                    "voxel-wise uncertainty, ``_unc-avgUnc.nii.gz`` for structure-wise uncertainty\n",
                    "derived from mean value of ``_unc-vox.nii.gz`` within a given connected object,\n",
                    "``_unc-cv.nii.gz`` for structure-wise uncertainty derived from coefficient of\n",
                    "variation, ``_unc-iou.nii.gz`` for structure-wise measure of uncertainty\n",
                    "derived from the Intersection-over-Union of the predictions, or ``_soft.nii.gz``\n",
                    "to threshold on the average of Monte Carlo iterations."
                ]
            }
        }
    }



.. code-block:: JSON

    {
        "postprocessing": {
            "threshold_uncertainty": {
                "thr": -1,
                "suffix": "_unc-vox.nii.gz"
            }
        }
    }


Evaluation Parameters
---------------------
Dict. Parameters to get object detection metrics (true positive and false detection rates), and this, for defined
object sizes.


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "target_size",
        "type": "dict",
        "options": {
            "thr": {
                "type": "list, int",
                "$$description": [
                    "These values will create several consecutive target size bins. For instance\n",
                    "with a list of two values, we will have three target size bins: minimal size\n",
                    "to first list element, first list element to second list element, and second\n",
                    "list element to infinity."
                ]
            },
            "unit": {
                "type": "string",
                "$$description": [
                    "Either `vox` for voxels or `mm3`. Indicates the unit used to define the\n",
                    "target object sizes."
                ]
            }
        }
    }

.. code-block:: JSON

    {
        "evaluation_parameters": {
            "target_size": {
                "thr": [20, 100],
                "unit": "vox"
            }
        }
    }


.. jsonschema::

    {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "title": "overlap",
        "type": "dict",
        "options": {
            "thr": {
                "type": "int",
                "$$description": [
                    "Minimal object size overlapping to be considered a TP, FP, or FN."
                ]
            },
            "unit": {
                "type": "string",
                "$$description": [
                    "Either `vox` for voxels or `mm3`. Indicates the unit used to define the\n",
                    "overlap."
                ]
            }
        }
    }

.. code-block:: JSON

    {
        "evaluation_parameters": {
            "overlap": {
                "thr": 30,
                "unit": "vox"
            }
        }
    }
Installation
============

Supported OS
++++++++++++

    Currently, ``ivadomed`` supports GPU/CPU on ``Linux`` and ``Windows``, and CPU only on ``macOS`` and `Windows Subsystem for Linux <https://docs.microsoft.com/en-us/windows/wsl/>`_.

Step 1: Setup dedicated python environment
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    You can setup ``ivadomed`` using either Conda or Venv:

    .. tabs::

        .. tab:: Install via ``venv``

            1. Setup Python Venv Virtual Environment.

                ``ivadomed`` requires Python >= 3.6 and <3.10.

                First, make sure that a compatible version of Python 3 is installed on your system by running:

                .. tabs::

                    .. tab:: Mac/Linux

                        .. code::

                            python3 --version

                    .. tab:: Windows

                        .. code::

                            python --version

                If your system's Python is not 3.6, 3.7, 3.8, or 3.9 (or if you don't have Python 3 installed at all), please `install Python <https://wiki.python.org/moin/BeginnersGuide/Download/>`_ before continuing.

                Once you have a supported version of Python installed, run the following command:


                .. tabs::

                    .. tab:: Mac/Linux

                        .. code::

                            # Replacing ``3.X`` with the Python version number that you installed):
                            python3.X -m venv ivadomed_env

                        .. note::

                           If you use ``Debian`` or ``Ubuntu``, you may be prompted to install the ``python3-venv`` module when creating the virtual environment. This is expected, so please follow the instructions provided by Python. For other operating systems, ``venv`` will be installed by default.

                    .. tab:: Windows

                        .. code::

                            python -m venv ivadomed_env

            2. Activate the new virtual environment (default named ``ivadomed_env``)

                .. tabs::

                    .. tab:: Mac/Linux

                        .. code::

                            source ivadomed_env/bin/activate

                    .. tab:: Windows

                        .. code::

                            cd ivadomed_env/Scripts/
                            activate

        .. tab:: Install via ``conda``

            1. Create new conda environment using ``environment.yml`` file

                ::

                    conda env create --name ivadomed_env

            2. Activate the new conda environment

                ::

                    conda activate ivadomed_env


        .. tab:: Compute Canada HPC

            There are numerous constraints and limited package availabilities with ComputeCanada cluster environment.

            It is best to attempt ``venv`` based installations and follow up with ComputeCanada technical support as MANY specially compiled packages (e.g. numpy) are exclusively available for Compute Canada HPC environment.

            If you are using `Compute Canada <https://www.computecanada.ca/>`_, you can load modules as `mentioned here <https://intranet.neuro.polymtl.ca/computing-resources/compute-canada#modules>`_ and `also here <https://docs.computecanada.ca/wiki/Utiliser_des_modules/en#Loading_modules_automatically>`_.


Step 2: Install ``ivadomed``
++++++++++++++++++++++++++++


    .. tabs::

        .. tab:: PyPI Installation

            Install ``ivadomed`` and its requirements from
            `PyPI <https://pypi.org/project/ivadomed/>`__:

            ::

                pip install --upgrade pip

                pip install ivadomed

        .. tab:: Repo Installation (Advanced or Developer)

            Bleeding-edge developments are available on the project's master branch
            on Github. Install ``ivadomed`` from source:

            ::

                git clone https://github.com/ivadomed/ivadomed.git

                cd ivadomed

                pip install -e .


Step 3: Install ``torch`` and ``torchvision`` with CPU or GPU Support
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    .. tabs::

        .. tab:: CPU Support

            If you plan to run ``ivadomed`` on CPU only, install PyTorch per instructions provided below for your specific operating system:

            .. tabs::

                .. tab:: Windows/Linux

                    .. code::

                       pip install torch==1.8.0+cpu torchvision==0.9.0+cpu --find-links https://download.pytorch.org/whl/torch_stable.html

                .. tab:: Mac

                    .. code::

                       pip install torch==1.8.0 torchvision==0.9.0 --find-links https://download.pytorch.org/whl/torch_stable.html

                .. tab:: Repo Installation (Advanced or Developer)

                    Run this only if you have already downloaded/cloned the repo with access to the ``requirement.txt`` file, then run the following command while at the repository root level:

                    .. code::

                        pip install -r requirements.txt

        .. tab:: Nvidia GPU Support

            ``ivadomed`` requires CUDA11 to execute properly. If you have a nvidia GPU, try to look up its Cuda Compute Score `here <https://developer.nvidia.com/cuda-gpus>`_, which needs to be > 3.5 to support CUDA11. Then, make sure to upgrade to nvidia driver to be at least v450+ or newer.

            If you have a compatible NVIDIA GPU that supports CUDA11 and with the right driver installed try run the following command relevant to your situation:

            .. tabs::

                .. tab:: All OS

                    .. code::

                       pip install torch==1.8.1+cu111 torchvision==0.9.1+cu111 --find-links https://download.pytorch.org/whl/torch_stable.html

                .. tab:: Repo Installation (Advanced or Developer)

                    Run this only if you have already downloaded/cloned the repo with access to the ``requirement_gpu.txt`` file, then run the following command while at the repository root level:
                    .. code::

                       pip install -r requirements_gpu.txt


Developer-only Installation Steps
+++++++++++++++++++++++++++++++++

    The additional steps below are only necessary for contributors to the ``ivadomed`` project.

    The ``pre-commit`` package is used to enforce a size limit on committed files. The ``requirements_dev.txt`` also contain additional dependencies related to documentation building and testing.

    After you've installed ``ivadomed``, install the ``pre-commit`` hooks by running:

    .. code::

        pip install -r requirements_dev.txt
        pre-commit installTechnical features
==================

Physics-informed network
------------------------

CNNs can be modulated, at each layer, using the `Feature-wise Linear
Modulation (FiLM) <https://arxiv.org/abs/1709.07871>`__ technique.
FiLM permits to add priors during training/inference based on the
imaging physics (e.g. acquisition parameters), thereby improving the
performance of the output segmentations.

.. figure:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/technical_features/film_figure.png
   :alt: Figure FiLM

   Figure FiLM

.. _Uncertainty-measures:

Uncertainty measures
--------------------

At inference time, uncertainty can be estimated via two ways: -
model-based uncertainty (epistemic) based on `Monte Carlo
Dropout <https://arxiv.org/abs/1506.02142>`__. - image-based uncertainty
(aleatoric) `based on test-time
augmentation <https://doi.org/10.1016/j.neucom.2019.01.103>`__.

From the Monte Carlo samples, different measures of uncertainty can be
derived: - voxel-wise entropy - structure-wise intersection over union -
structure-wise coefficient of variation - structure-wise averaged
voxel-wise uncertainty within the structure

These measures can be used to perform some
`post-processing <https://arxiv.org/abs/1808.01200>`__ based on the
uncertainty measures.

.. figure:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/technical_features/uncertainty_measures.png
   :alt: Figure Uncertainty

   Figure Uncertainty

Two-step training scheme with class sampling
--------------------------------------------

Class sampling, coupled with a transfer learning strategy, can mitigate
class imbalance issues, while addressing the limitations of classical
under-sampling (risk of loss of information) or over-sampling (risk of
overfitting) approaches.

During a first training step, the CNN is trained on an equivalent
proportion of positive and negative samples, negative samples being
under-weighted dynamically at each epoch. During the second step, the
CNN is fine-tuned on the realistic (i.e. class-imbalanced) dataset.

Mixup
-----

`Mixup <https://arxiv.org/abs/1710.09412>`__ is a data augmentation
technique, wherein training is performed on samples that are generated
by combining two random samples from the training set and from the
associated labels. The motivation is to regularize the network while
extending the training distribution.

.. figure:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/technical_features/mixup.png
   :alt: Figure mixup

   Figure mixup

Data augmentation on lesion labels
----------------------------------

This data augmentation is motivated by the large inter-rater variability
that is common in medical image segmentation tasks. Typically, raters
disagree on the boundaries of pathologies (e.g., tumors, lesions). A
soft mask is constructed by morphological dilation of the binary
segmentation (i.e. mask provided by expert), where expert-labeled voxels
have one as value while the augmented voxels are assigned a soft value
which depends on the distance to the core of the lesion. Thus, the prior
knowledge about the subjective lesion borders is then leveraged to the
network.

.. figure:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/technical_features/dilate-gt.png
   :alt: Figure Data Augmentation on lesion ground truths

   Figure Data Augmentation on lesion ground truths

Network architectures
---------------------

-  `UNet <https://arxiv.org/abs/1505.04597>`__, with control of the
   network depth.
-  HeMIS-UNet: integrates the
   `HeMIS <https://arxiv.org/abs/1607.05194>`__ strategy to deal with
   missing modalities within a UNet training scheme.
-  FiLMed-UNet, based on `FiLM <https://arxiv.org/abs/1709.07871>`__
   strategy adapted to the `segmentation
   task <#physic-informed-network>`__.
- Countception: modified implementation of `Countception <https://arxiv.org/abs/1703.08710>`__ for keypoints detection.

Loss functions
--------------

-  `Dice Loss <https://arxiv.org/abs/1606.04797>`__. Also adapted for
   multi-label segmentation tasks, by averaging the loss for each class.
-  `Focal Loss <https://arxiv.org/abs/1708.02002>`__.
-  Focal-Dice Loss: Linear combination of the Focal and Dice losses.
-  `Generalized Dice Loss <https://arxiv.org/abs/1707.03237>`__. An
   additional feature compared to the published reference, is that the
   background volume can be weighted by the inverse of its area, which
   could be of interest in high class imbalance scenarios.
-  `Adaptive wing loss <https://arxiv.org/abs/1904.07399>`__. Loss function used to detect key points with Gaussian representation of the target.
-  Loss Combination: Linear combination of any other implemented losses. 
API Reference
=============

This document is for developers of ``ivadomed``, it contains the API functions.

Loader API
++++++++++

loader.film
^^^^^^^^^^^

.. automodule:: ivadomed.loader.film


loader.loader
^^^^^^^^^^^^^
.. automodule:: ivadomed.loader.loader


loader.utils
^^^^^^^^^^^^
.. automodule:: ivadomed.loader.utils


Object Detection API
++++++++++++++++++++

object_detection.utils
^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: ivadomed.object_detection.utils


Evaluation API
++++++++++++++

.. automodule:: ivadomed.evaluation


Losses API
++++++++++

.. automodule:: ivadomed.losses


Main API
++++++++

.. automodule:: ivadomed.main


Metrics API
+++++++++++

.. automodule:: ivadomed.metrics


Postprocessing API
++++++++++++++++++

.. automodule:: ivadomed.postprocessing


Testing API
+++++++++++

.. automodule:: ivadomed.testing


Training API
++++++++++++

.. automodule:: ivadomed.training


Transformations API
+++++++++++++++++++

.. automodule:: ivadomed.transforms


Utils API
+++++++++

.. automodule:: ivadomed.utils

Visualize API
+++++++++++++

.. automodule:: ivadomed.visualize

Inference API
+++++++++++++

.. automodule:: ivadomed.inference

Mixup API
+++++++++

.. automodule:: ivadomed.mixup

Uncertainty API
+++++++++++++++

.. automodule:: ivadomed.uncertainty

Maths API
+++++++++

.. automodule:: ivadomed.maths

.. image:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/index/overview_title.png
  :alt: Alternative text

|

``ivadomed`` is an integrated framework for medical image analysis with deep
learning, based on `PyTorch <https://pytorch.org/>`_. The name is a portmanteau between *IVADO* (The `Institute for data
valorization <https://ivado.ca/en/>`_) and *Medical*.

If you use ``ivadomed`` for your research, please cite `our paper <https://joss.theoj.org/papers/10.21105/joss.02868>`_:

.. code::

    @article{gros2021ivadomed,
      doi = {10.21105/joss.02868},
      url = {https://doi.org/10.21105/joss.02868},
      year = {2021},
      publisher = {The Open Journal},
      volume = {6},
      number = {58},
      pages = {2868},
      author = {Charley Gros and Andreanne Lemay and Olivier Vincent and Lucas Rouhier and Marie-Helene Bourget and Anthime Bucquet and Joseph Paul Cohen and Julien Cohen-Adad},
      title = {ivadomed: A Medical Imaging Deep Learning Toolbox},
      journal = {Journal of Open Source Software}
    }

Home
====

.. toctree::
   :maxdepth: 1
   :caption: Overview

   purpose.rst
   technical_features.rst
   use_cases.rst

.. toctree::
   :maxdepth: 1
   :caption: Getting started

   installation.rst
   data.rst
   configuration_file.rst
   usage.rst
   architectures.rst
   pretrained_models.rst
   scripts.rst
   help.rst

.. _tutorials:
.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   tutorials/one_class_segmentation_2d_unet.rst
   tutorials/cascaded_architecture.rst
   tutorials/uncertainty.rst
   tutorials/automate_training.rst
   tutorials/two_class_microscopy_seg_2d_unet.rst

.. toctree::
   :maxdepth: 1
   :caption: Developer section

   contributing.rst
   api_ref.rst
   contributors.rst
   license.rst
.. |br| raw:: html

   <br />

Scripts
=======

This section contains a collection of useful scripts for quality control during
the training of models.

ivadomed_visualize_transforms
"""""""""""""""""""""""""""""

.. autofunction:: ivadomed.scripts.visualize_transforms.run_visualization

ivadomed_convert_to_onnx
""""""""""""""""""""""""

.. autofunction:: ivadomed.scripts.convert_to_onnx.convert_pytorch_to_onnx

ivadomed_automate_training
""""""""""""""""""""""""""

.. autofunction:: ivadomed.scripts.automate_training.automate_training

.. autofunction:: ivadomed.scripts.automate_training.HyperparameterOption

.. autofunction:: ivadomed.scripts.automate_training.get_param_list

.. autofunction:: ivadomed.scripts.automate_training.make_config_list

.. autofunction:: ivadomed.scripts.automate_training.update_dict


ivadomed_compare_models
"""""""""""""""""""""""

.. autofunction:: ivadomed.scripts.compare_models.compute_statistics

ivadomed_prepare_dataset_vertebral_labeling
"""""""""""""""""""""""""""""""""""""""""""

.. autofunction:: ivadomed.scripts.prepare_dataset_vertebral_labeling.extract_mid_slice_and_convert_coordinates_to_heatmaps

ivadomed_extract_small_dataset
""""""""""""""""""""""""""""""

.. autofunction:: ivadomed.scripts.extract_small_dataset.extract_small_dataset

ivadomed_training_curve
"""""""""""""""""""""""

.. autofunction:: ivadomed.scripts.training_curve.run_plot_training_curves

ivadomed_download_data
""""""""""""""""""""""

.. autofunction:: ivadomed.scripts.download_data.install_data

ivadomed_visualize_and_compare_testing_models
"""""""""""""""""""""""""""""""""""""""""""""

.. autofunction:: ivadomed.scripts.visualize_and_compare_testing_models.visualize_and_compare_models
Two-class microscopy segmentation with 2D U-Net
===============================================

In this tutorial we will learn the following features:

- Training of a segmentation model (U-Net 2D) with two-class labels on a single contrast on microscopy PNG images,

- Testing of a trained model and computation of evaluation metrics,

- Visualization of the outputs of a trained model.

Download dataset
----------------

We will use a publicly available dataset consisting of 10 microscopy samples of rat spinal cord.

To download the dataset (~11MB), run the following command in your terminal:

.. code-block:: bash

   # Download data
   ivadomed_download_data -d data_axondeepseg_sem

Configuration file
------------------

In ``ivadomed``, training is orchestrated by a configuration file. Examples of configuration files are available in
the ``ivadomed/config/`` and the documentation is available in :doc:`../configuration_file`.

In this tutorial, we will use the configuration file: ``ivadomed/config/config_microscopy.json``.
First off, copy this configuration file in your local directory (to avoid modifying the source file):

.. code-block:: bash

   cp <PATH_TO_IVADOMED>/ivadomed/config/config_microscopy.json .

Then, open it with a text editor.
Below we will discuss some of the key parameters to perform a two-class 2D
microscopy segmentation training.

- ``command``: Action to perform. Here, we want to train a model, so we set the fields as follows:

  .. code-block:: xml

     "command": "train"

- ``path_output``: Folder name that will contain the output files (e.g., trained model, predictions, results).

  .. code-block:: xml

     "path_output": "log_microscopy_sem"

- ``loader_parameters:path_data``: Location of the dataset. As discussed in `Data <../data.html>`__, the dataset
  should conform to the BIDS standard. Modify the path so it points to the location of the downloaded dataset.

  .. code-block:: xml

     "path_data": ["data_axondeepseg_sem"]

- ``loader_parameters:bids_config``: Location of the custom BIDS configuration file required for microscopy
  file indexing.

  .. note::

     You will need to update the value ``"ivadomed/config/config_bids.json"`` to
     ``"<PATH_TO_IVADOMED>/ivadomed/config/config_bids.json"``

  .. code-block:: xml

     "bids_config": "<PATH_TO_IVADOMED>/ivadomed/config/config_bids.json"

- ``loader_parameters:target_suffix``: Suffix of the ground truth segmentations. The ground truths are located
  under the ``data_axondeepseg_sem/derivatives/labels`` folder. In our case, the suffix are ``_seg-axon-manual``
  and ``_seg-myelin-manual``:

  .. code-block:: xml

     "target_suffix": ["_seg-axon-manual", "_seg-myelin-manual"]

- ``loader_parameters:extensions``: File extensions of the microscopy raw data.

  .. code-block:: xml

     "extensions": [".png"]

- ``loader_parameters:contrast_params``: Contrast(s) of interest. In our case, we are training a single contrast model
  with contrast ``SEM``.

  .. code-block:: xml

     "contrast_params": {
         "training_validation": ["SEM"],
         "testing": ["SEM"],
         "balance": {}
     }

- ``loader_parameters:slice_axis``: Orientation of the 2D slice to use with the model.
  2D PNG files must use default ``axial``.

  .. code-block:: xml

     "slice_axis": "axial"

- ``split_dataset:split_method``: Describe the metadata used to split the train/validation/test sets.
  Here, ``sample_id`` from the ``samples.tsv`` file will shuffle all samples, then split them between
  train/validation/test sets.
- ``split_dataset:train_fraction``: Fraction of the dataset's ``sample_id`` in the train set. In our case ``0.6``.
- ``split_dataset:test_fraction``: Fraction of the dataset's ``sample_id`` in the test set. In our case ``0.1``.

  .. code-block:: xml

      "split_method" : "sample_id"
      "train_fraction": 0.6
      "test_fraction": 0.1

- ``training_parameters:training_time:num_epochs``: The maximum number of epochs that will be run during training. Each epoch is composed
  of a training part and a validation part. It should be a strictly positive integer. In our case, we will use
  50 epochs.

  .. code-block:: xml

     "num_epochs": 50

- ``default_model:length_2D``: Size of the 2D patches used as model’s input tensors. We recommend using patches
  between 256x256 and 512x512. In our case, we use patches of 256x256.
- ``default_model:stride_2D``: Pixels’ shift over the input matrix to create 2D patches. In our case, we use
  a stride of 244 pixels in both dimensions, resulting in an overlap of 12 pixels between patches.

  .. code-block:: xml

     "length_2D": [256, 256]
     "stride_2D": [244, 244]

- ``postprocessing:binarize_maxpooling``: Used to binarize predictions across all classes in multiclass models. For each pixel, the class, including the background class, with the highest output probability will be segmented.

  .. code-block:: xml

      "binarize_maxpooling": {}

- ``transformation:Resample``: Used to resample images to a common resolution (in mm) before splitting into patches,
  according to each image real pixel size. In our case, we resample the images to a common resolution of 0.0001 mm
  (0.1 μm) in both dimensions.

  .. code-block:: xml

     "Resample":
        {
            "hspace": 0.0001,
            "wspace": 0.0001
        },


Train model
-----------

Once the configuration file is ready, run the training:

.. code-block:: bash

   ivadomed -c config_microscopy.json

Alternatively, the "command", "path_output", and "path_data" arguments can be passed as CLI flags
in which case they supersede the configration file, see `Usage <../usage.html>`__.

.. code-block:: bash

   ivadomed --train -c config_microscopy.json --path-data path/to/bids/data --path-output path/to/output/directory

.. note::

   If a `compatible GPU <https://pytorch.org/get-started/locally/>`_ is available, it will be used by default.
   Otherwise, training will use the CPU, which will take a prohibitively long computational time (several hours).

The main parameters of the training scheme and model will be displayed on the terminal, followed by the loss value
on training and validation sets at every epoch. To know more about the meaning of each parameter, go to
:doc:`../configuration_file`. The value of the loss should decrease during the training.

.. code-block:: console

   No CLI argument given for command: ( --train | --test | --segment ). Will check config file for command...
   CLI flag --path-output not used to specify output directory. Will check config file for directory...
   CLI flag --path-data not used to specify BIDS data directory. Will check config file for directory...

   Creating output path: log_microscopy_sem
   Using GPU ID 0

   Selected architecture: Unet, with the following parameters:
   dropout_rate: 0.2
   bn_momentum: 0.1
   depth: 4
   is_2d: True
   final_activation: sigmoid
   length_2D: [256, 256]
   stride_2D: [244, 244]
   folder_name: model_seg_rat_axon-myelin_sem
   in_channel: 1
   out_channel: 3

   Dataframe has been saved in log_microscopy_sem/bids_dataframe.csv.
   After splitting: train, validation and test fractions are respectively 0.6, 0.3 and 0.1 of sample_id.

   Selected transformations for the ['training'] dataset:
   Resample: {'hspace': 0.0001, 'wspace': 0.0001}
   RandomAffine: {'degrees': 2.5, 'scale': [0.05, 0.05], 'translate': [0.015, 0.015], 'applied_to': ['im', 'gt']}
   ElasticTransform: {'alpha_range': [100.0, 150.0], 'sigma_range': [4.0, 5.0], 'p': 0.5, 'applied_to': ['im', 'gt']}
   NormalizeInstance: {'applied_to': ['im']}
   Selected transformations for the ['validation'] dataset:
   Resample: {'hspace': 0.0001, 'wspace': 0.0001}
   NormalizeInstance: {'applied_to': ['im']}

   Loading dataset: 100%|████████████████████████████████████████████████████████████████| 3/3 [00:00<00:00, 738.48it/s]
   Loaded 76 axial patches of shape [256, 256] for the validation set.
   Loading dataset: 100%|████████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 829.21it/s]
   Loaded 252 axial patches of shape [256, 256] for the training set.
   Creating model directory: log_microscopy_sem/model_seg_rat_axon-myelin_sem

   Initialising model's weights from scratch.
   Scheduler parameters: {'name': 'CosineAnnealingLR', 'base_lr': 1e-05, 'max_lr': 0.01}

   Selected Loss: DiceLoss
   with the parameters: []
   Epoch 1 training loss: -0.6894.
   Epoch 1 validation loss: -0.7908.

After 50 epochs (see ``"num_epochs"`` in the configuration file), the Dice score on the validation set should be ~85%.

.. note::

   When loading the images for training or evaluation, a temporary NIfTI file will be created for each images in the
   dataset directory (``path_data``) alongside the original PNG files.

Evaluate model
--------------

To test the trained model on the testing sub-dataset and compute evaluation metrics, run:

.. code-block:: bash

   ivadomed -c config_microscopy.json --test

If you prefer to use config files over CLI flags, set "command" to the following in you config file:

.. code-block:: xml

   "command": "test"

Then run:

.. code-block:: bash

   ivadomed -c config_microscopy.json

The model's parameters will be displayed in the terminal, followed by a preview of the results for each image.
The resulting segmentations are saved for each image in the ``<PATH_TO_OUT_DIR>/pred_masks`` while a CSV file,
saved in ``<PATH_TO_OUT_DIR>/results_eval/evaluation_3Dmetrics.csv``, contains all the evaluation metrics.
For more details on the evaluation metrics, see :mod:`ivadomed.metrics`.

.. code-block:: console

   CLI flag --path-output not used to specify output directory. Will check config file for directory...
   CLI flag --path-data not used to specify BIDS data directory. Will check config file for directory...

   Output path already exists: log_microscopy_sem
   Using GPU ID 0

   Selected architecture: Unet, with the following parameters:
   dropout_rate: 0.2
   bn_momentum: 0.1
   depth: 4
   is_2d: True
   final_activation: sigmoid
   length_2D: [256, 256]
   stride_2D: [244, 244]
   folder_name: model_seg_rat_axon-myelin_sem
   in_channel: 1
   out_channel: 3

   Dataframe has been saved in log_microscopy_sem/bids_dataframe.csv.
   After splitting: train, validation and test fractions are respectively 0.6, 0.3 and 0.1 of sample_id.

   Selected transformations for the ['testing'] dataset:
   Resample: {'hspace': 0.0001, 'wspace': 0.0001}
   NormalizeInstance: {'applied_to': ['im']}

   Loading dataset: 100%|████████████████████████████████████████████████████████████████| 1/1 [00:00<00:00, 413.48it/s]
   Loaded 16 axial patches of shape [256, 256] for the testing set.
   Loading model: log_microscopy_sem/best_model.pt

   Inference - Iteration 0: 100%|████████████████████████████████████████████████████████████████| 4/4 [00:01<00:00,  2.89it/s]
   Lossy conversion from float64 to uint8. Range [0, 1]. Convert image to uint8 prior to saving to suppress this warning.
   Lossy conversion from float64 to uint8. Range [0, 1]. Convert image to uint8 prior to saving to suppress this warning.
   {'dice_score': 0.8381376827003003, 'multi_class_dice_score': 0.8422281034034607, 'precision_score': 0.8342335786851753,
   'recall_score': 0.8420784999205466, 'specificity_score': 0.9456594910680598, 'intersection_over_union': 0.7213743575471384,
   'accuracy_score': 0.9202670087814067, 'hausdorff_score': 0.0}

   Run Evaluation on log_microscopy_sem/pred_masks

   Evaluation: 100%|████████████████████████████████████████████████████████████████| 1/1 [00:13<00:00, 13.56s/it]
   Lossy conversion from float64 to uint8. Range [0.0, 3.0]. Convert image to uint8 prior to saving to suppress this warning.
   Lossy conversion from float64 to uint8. Range [0.0, 3.0]. Convert image to uint8 prior to saving to suppress this warning.
                                avd_class0  avd_class1  dice_class0  dice_class1  ...  vol_gt_class0  vol_gt_class1  vol_pred_class0  vol_pred_class1
   image_id
   sub-rat3_sample-data9_SEM    0.082771    0.082971    0.868964     0.815492     ...  1.256960e-07   1.574890e-07   1.152920e-07     1.705560e-07

   [1 rows x 26 columns]

The test image segmentations are stored in ``<PATH_TO_OUT_DIR>/pred_masks/`` in PNG format and have the same name as
the input image with the suffix ``<class-index>_pred.png``. In our case: ``sub-rat3_sample-data9_SEM_class-0_pred.png`` and
``sub-rat3_sample-data9_SEM_class-1_pred.png`` for axons and myelin respectively (in the same order as ``target_suffix``).
A temporary NIfTI files containing the predictions for both classes with the suffix ``_pred.nii.gz`` will also be
present.

After the training for 50 epochs, the segmentations should be similar to the one presented in the following image.
The ground truth segmentations and predictions of the axons and myelin are presented in blue and red respectively for
``sub-rat3_sample-data9_SEM``):

.. image:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/tutorials/two_classes_microscopy_seg_2d_unet/axon_myelin_predictions.png
   :align: center
Estimate uncertainty
====================

    This tutorial shows how to estimate uncertainty measures on the model predictions. The uncertainty measures implemented
    in ``ivadomed`` are detailed in :ref:`Technical features <Uncertainty-measures>`.

    An interactive Colab version of this tutorial is directly accessible `here <https://colab.research.google.com/github/ivadomed/ivadomed/blob/master/testing/tutorials/tutorial_3_uncertainty_estimation.ipynb>`_.

Download dataset
----------------

    A dataset example is available for this tutorial. If not already done, download the dataset with the following line.
    For more details on this dataset see :ref:`One-class segmentation with 2D U-Net<Download dataset>`.

    .. code-block:: bash

       # Download data
       ivadomed_download_data -d data_example_spinegeneric


Configuration file
------------------
    In this tutorial we will use the configuration file: ``ivadomed/config/config.json``.
    First off, copy this configuration file in your local directory (to avoid modifying the source file):

    .. code-block:: bash

       cp <PATH_TO_IVADOMED>/ivadomed/config/config.json .

    Please open it with a text editor.
    The configuration file will be modified to be the same as the one used for
    :ref:`Technical features <Uncertainty-measures>`. As described in the tutorial
    :doc:`../tutorials/one_class_segmentation_2d_unet`, make sure ``path_data`` point to the location of the dataset.
    The parameters that are specific to this tutorial are:

    - ``path_output``: Location of the directory containing the trained model. To avoid having to train a model from
      scratch, there is a already trained model for spinal cord segmentation in the folder named `trained_model`, in the downloaded dataset.
      Modify the path so it points to the location of the trained model.

      .. code-block:: json

         "path_output": "<PATH_TO_DATASET>/data_example_spinegeneric/trained_model"

      Note that you can also pass this argument via CLI (see :ref:`Usage <usage>`)

      .. code-block:: bash

        ivadomed -c path/to/config --path-output path/to/output/directory

    - ``command``: Action to perform. Here, we want to do some inference using the previously trained model, so we set the
      field as follows:

      .. code-block:: json

         "command": "test"

      Note that you can also pass this argument via CLI (see :ref:`Usage <usage>`)

      .. code-block:: bash

        ivadomed --test -c path/to/config

    - ``uncertainty``: Type of uncertainty to estimate. Available choices are ``epistemic`` and
      ``aleatoric``. Note that both can be ``true``. More details on the implementation are available in :ref:`Technical features <Uncertainty-measures>`.
      ``n_it`` controls the number of Monte Carlo iterations that are performed to estimate the uncertainty. Set it to a
      non-zero positive integer for this tutorial (e.g. ``20``).

      .. code-block:: json

          "uncertainty": {
               "epistemic": true,
               "aleatoric": true,
               "n_it": 20
          }


    - ``transformation``: Data augmentation transformation. If you have selected the aleatoric uncertainty, the data
      augmentation that will be performed is the same as the one performed for the training. Note that only transformations
      for which a ``undo_transform`` (i.e. inverse transformation) is available will be performed since these inverse
      transformations are required to reconstruct the predicted volume.


Run uncertainty estimation
--------------------------

    Once the configuration file has been modified, run the inference with the following command:

    .. code-block:: bash

       ivadomed --test -c config.json --path-data <PATH_TO_DATASET>/data_example_spinegeneric --path-output <PATH_TO_DATASET>/data_example_spinegeneric/trained_model

    - Here, we want to do some inference using the previously trained model, so we set the
      command flag as follows:

      .. code-block:: bash

         --test

    - ``--path-data``: Location of the directory containing the dataset.

      .. code-block:: bash

         --path-data <PATH_TO_DATASET>/data_example_spinegeneric

    - ``--path-output``: Folder name that will contain the output files (e.g., trained model, predictions, results). For the purpose of this particular tutorial, since we do not train the model from scratch, we set the output path to point to a folder containing the pre-trained model for spinal cord segmentation that comes with the dataset. Hence, after running this tutorial, the corresponding output files can be found inside the `trained_model` folder.

      .. code-block:: bash

         --path-output <PATH_TO_DATASET>/data_example_spinegeneric/trained_model

    If you set the ``command``, ``path_output``, and ``path_data`` arguments in your config file, you do not need to pass the CLI flags:

    .. code-block:: bash

       ivadomed -c config.json

    If aleatoric uncertainty was selected, then data augmentation operations will be performed at testing time, as indicated
    in the terminal output (see below). Note that ``ElasticTransform`` has been deactivated because no ``undo_transform``
    function is available for it.

    .. code-block:: bash

        Selected transformations for the ['testing'] dataset:
            Resample: {'hspace': 0.75, 'wspace': 0.75, 'dspace': 1}
            CenterCrop: {'size': [128, 128]}
            RandomAffine: {'degrees': 5, 'scale': [0.1, 0.1], 'translate': [0.03, 0.03], 'applied_to': ['im', 'gt']}
            ElasticTransform: {'alpha_range': [28.0, 30.0], 'sigma_range': [3.5, 4.5], 'p': 0.1, 'applied_to': ['im', 'gt']}
            NumpyToTensor: {}
            NormalizeInstance: {'applied_to': ['im']}
        ElasticTransform transform not included since no undo_transform available for it.

    ... otherwise, only preprocessing and data normalization are performed, see below:

    .. code-block:: bash

        Selected transformations for the ['testing'] dataset:
            Resample: {'hspace': 0.75, 'wspace': 0.75, 'dspace': 1}
            CenterCrop: {'size': [128, 128]}
            NumpyToTensor: {}
            NormalizeInstance: {'applied_to': ['im']}

    For each Monte Carlo iteration, each testing image is segmented using the trained model and saved under ``pred_masks``,
    with the iteration number as suffix (e.g. ``sub-001_pred_00.nii.gz`` ... ``sub-001_pred_19.nii.gz``).

    .. code-block:: bash

        Computing model uncertainty over 20 iterations.
        Inference - Iteration 0: 100%|███████████████████████████████████████████████████████████████████████████████████| 5/5 [00:11<00:00,  2.27s/it]
        Inference - Iteration 1: 100%|███████████████████████████████████████████████████████████████████████████████████| 5/5 [00:09<00:00,  1.81s/it]
        Inference - Iteration 2: 100%|███████████████████████████████████████████████████████████████████████████████████| 5/5 [00:09<00:00,  1.96s/it]
        Inference - Iteration 3: 100%|███████████████████████████████████████████████████████████████████████████████████| 5/5 [00:08<00:00,  1.66s/it]
        Inference - Iteration 4: 100%|███████████████████████████████████████████████████████████████████████████████████| 5/5 [00:08<00:00,  1.69s/it]
        Inference - Iteration 5: 100%|███████████████████████████████████████████████████████████████████████████████████| 5/5 [00:09<00:00,  1.92s/it]
        Inference - Iteration 6: 100%|███████████████████████████████████████████████████████████████████████████████████| 5/5 [00:08<00:00,  1.74s/it]
        Inference - Iteration 7: 100%|███████████████████████████████████████████████████████████████████████████████████| 5/5 [00:08<00:00,  1.74s/it]
        Inference - Iteration 8: 100%|███████████████████████████████████████████████████████████████████████████████████| 5/5 [00:09<00:00,  1.83s/it]
        Inference - Iteration 9: 100%|███████████████████████████████████████████████████████████████████████████████████| 5/5 [00:07<00:00,  1.59s/it]
        Inference - Iteration 10: 100%|██████████████████████████████████████████████████████████████████████████████████| 5/5 [00:09<00:00,  1.85s/it]
        Inference - Iteration 11: 100%|██████████████████████████████████████████████████████████████████████████████████| 5/5 [00:09<00:00,  1.85s/it]
        Inference - Iteration 12: 100%|██████████████████████████████████████████████████████████████████████████████████| 5/5 [00:09<00:00,  1.92s/it]
        Inference - Iteration 13: 100%|██████████████████████████████████████████████████████████████████████████████████| 5/5 [00:09<00:00,  1.83s/it]
        Inference - Iteration 14: 100%|██████████████████████████████████████████████████████████████████████████████████| 5/5 [00:09<00:00,  1.84s/it]
        Inference - Iteration 15: 100%|██████████████████████████████████████████████████████████████████████████████████| 5/5 [00:09<00:00,  1.87s/it]
        Inference - Iteration 16: 100%|██████████████████████████████████████████████████████████████████████████████████| 5/5 [00:09<00:00,  1.81s/it]
        Inference - Iteration 17: 100%|██████████████████████████████████████████████████████████████████████████████████| 5/5 [00:09<00:00,  1.95s/it]
        Inference - Iteration 18: 100%|██████████████████████████████████████████████████████████████████████████████████| 5/5 [00:09<00:00,  1.82s/it]
        Inference - Iteration 19: 100%|██████████████████████████████████████████████████████████████████████████████████| 5/5 [00:08<00:00,  1.71s/it]

    The Monte Carlo samples are then used to compute uncertainty measures for each image. The results are saved under
    ``pred_masks``.

    .. code-block:: bash

        Uncertainty Computation: 100%|███████████████████████████████████████████████████████████████████████████████████| 5/5 [01:31<00:00, 18.28s/it]

    Six files are generated during this process for each testing image:

    - ``*_soft.nii.gz``: Soft segmentation (i.e. values between 0 and 1) which is generated by averaging the Monte Carlo
      samples.
    - ``*_pred.nii.gz``: Binary segmentation obtained by thresholding ``*_soft.nii.gz`` with ``1 / (Number of Monte Carlo
      iterations)``.
    - ``*_unc-vox.nii.gz``: Voxel-wise measure of uncertainty derived from the entropy of the Monte Carlo samples. The
      higher a given voxel value is, the more uncertain is the prediction for this voxel.
    - ``*_unc-avgUnc.nii.gz``: Structure-wise measure of uncertainty derived from the mean value of ``*_unc-vox.nii.gz``
      within a given connected object (e.g. a lesion, grey matter).
    - ``*_unc-cv.nii.gz``: Structure-wise measure of uncertainty derived from the coefficient of variation of the volume
      of a given connected object across the Monte Carlo samples. The higher a given voxel value is, the more uncertain is the
      prediction for this voxel.
    - ``*_unc-iou.nii.gz``: Structure-wise measure of uncertainty derived from the Intersection-over-Union of the
      predictions of a given connected object across the Monte Carlo samples. The lower a given voxel value is, the more
      uncertain is the prediction for this voxel.

    These files can further be used for post-processing to refine the segmentation. For example, the voxels
    depicted in pink are more uncertain than the ones in blue (left image): we might want to refine the model prediction by removing
    from the foreground class the voxels with low uncertainty (blue, left image) AND low prediction value (dark red, middle image).

    .. image:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/tutorials/uncertainty/uncertainty_tutorial.png
       :align: center
Optimize Hyperparameters
========================

For training, you may want to determine which hyperparameters will be the best to use. To do this,
we have the function ``ivadomed_automate_training``.


Step 1: Download Example Data
-----------------------------

To download the dataset (~490MB), run the following command in your terminal:

.. code-block:: bash

   ivadomed_download_data -d data_example_spinegeneric


Step 2: Copy Configuration File
-------------------------------

In ``ivadomed``, training is orchestrated by a configuration file. Examples of configuration files
are available in the ``ivadomed/config/`` and the documentation is available in :doc:`../configuration_file`.

In this tutorial we will use the configuration file: ``ivadomed/config/config.json``.
Copy this configuration file in your local directory (to avoid modifying the source file):

.. code-block:: bash

   cp <PATH_TO_IVADOMED>/ivadomed/config/config.json .

Then, open it with a text editor. Below we will discuss some of the key parameters to perform a one-class 2D
segmentation training.


Step 3: Create Hyperparameters Config File
------------------------------------------

The hyperparameter config file should have the same layout as the config file. To select
a hyperparameter you would like to vary, just list the different options under the
appropriate key.

In our example, we have 3 hyperparameters we would like to vary: ``batch_size``, ``loss``, and
``depth``. In your directory, create a new file called: ``config_hyper.json``. Open this
in a text editor and add the following:

.. code-block:: JSON

    {
        "training_parameters": {
            "batch_size": [2, 64],
            "loss": [
                {"name": "DiceLoss"},
                {"name": "FocalLoss", "gamma": 0.2, "alpha" : 0.5}
            ]
        },
        "default_model": {"depth": [2, 3, 4]}
    }

Step 4: (Optional) Change the Training Epochs
---------------------------------------------

The default number of training epochs in the ``config.json`` file is ``100``; however,
depending on your computer, this could be quite slow (especially if you don't have any GPUs).

To change the number of epochs, open the ``config.json`` file and change the following:

.. code-block:: JSON

    {
        "training_parameters": {
            "training_time": {
                "num_epochs": 1
            }
        }
    }


Step 5: Run the Code
--------------------

Default
^^^^^^^

If neither ``all-combin`` nor ``multi-params`` is selected, then the hyperparameters will be
combined as follows into a ``config_list``.

.. note::

    I am not showing the actual ``config_list`` here as it would take up too much space. The options
    listed below are incorporated into the base config file in ``config.json``.

.. code-block::

    batch_size = 2, loss = "DiceLoss", depth = 3
    batch_size = 64, loss = "DiceLoss", depth = 3
    batch_size = 18, loss = "FocalLoss", depth = 3
    batch_size = 18, loss = "DiceLoss", depth = 2
    batch_size = 18, loss = "DiceLoss", depth = 3
    batch_size = 18, loss = "DiceLoss", depth = 4

To run this:

.. code-block:: bash

    ivadomed_automate_training -c config.json -ch config_hyper.json \
    -n 1

All Combinations
^^^^^^^^^^^^^^^^

If the flag ``all-combin`` is selected, the hyperparameter options will be combined
combinatorically.

.. code-block::

    batch_size = 2, loss = "DiceLoss", depth = 2
    batch_size = 2, loss = "FocalLoss", depth = 2
    batch_size = 2, loss = "DiceLoss", depth = 3
    batch_size = 2, loss = "FocalLoss", depth = 3
    batch_size = 2, loss = "DiceLoss", depth = 4
    batch_size = 2, loss = "FocalLoss", depth = 4
    batch_size = 2, loss = "DiceLoss", depth = 4
    batch_size = 2, loss = "FocalLoss", depth = 4
    batch_size = 64, loss = "DiceLoss", depth = 2
    batch_size = 64, loss = "FocalLoss", depth = 2
    batch_size = 64, loss = "DiceLoss", depth = 3
    batch_size = 64, loss = "FocalLoss", depth = 3
    batch_size = 64, loss = "DiceLoss", depth = 4
    batch_size = 64, loss = "FocalLoss", depth = 4
    batch_size = 64, loss = "DiceLoss", depth = 4
    batch_size = 64, loss = "FocalLoss", depth = 4

To run:

.. code-block:: bash

    ivadomed_automate_training -c config.json -ch config_hyper.json \
    -n 1 --all-combin

Multiple Parameters
^^^^^^^^^^^^^^^^^^^

If the flag ``multi-params`` is selected, the elements from each hyperparameter list will be
selected sequentially, so all the first elements, then all the second elements, etc. If the lists
are different lengths, say ``len(list_a) = n`` and ``len(list_b) = n+m``, where ``n`` and ``m``
are strictly positive integers, then we will only use the first ``n`` elements.

.. code-block::

    batch_size = 2, loss = "DiceLoss", depth = 2
    batch_size = 64, loss = "FocalLoss", depth = 3

To run:

.. code-block:: bash

    ivadomed_automate_training -c config.json -ch config_hyper.json \
    -n 1 --multi-params

Step 6: Results
---------------

There will be an output file called ``detailed_results.csv``. This file gives an overview of the
results from all the different trials. For a more fine-grained analysis, you can also look
at each of the log directories (there is one for each config option).

An example of the ``detailed_results.csv``:

.. csv-table::
   :file: detailed_results.csv
One-class segmentation with 2D U-Net
====================================

    In this tutorial we will learn the following features:

    - Training of a segmentation model (U-Net 2D) with a single label on multiple contrasts,
    - Testing of a trained model and computation of 3D evaluation metrics.
    - Visualization of the outputs of a trained model.

    An interactive Colab version of this tutorial is directly accessible `here <https://colab.research.google.com/github/ivadomed/ivadomed/blob/master/testing/tutorials/tutorial_1_2d_segmentation_unet.ipynb>`_.

.. _Download dataset:

Download dataset
----------------

    We will use a publicly-available dataset consisting of MRI data of the spinal cord. This dataset is a subset of the
    `spine-generic multi-center dataset <https://github.com/spine-generic/data-multi-subject>`_ and has been pre-processed
    to facilitate training/testing of a new model. Namely, for each subject, all six contrasts were co-registered together.
    Semi-manual cord segmentation for all modalities and manual cerebrospinal fluid labels for T2w modality were created.
    More details `here <https://github.com/ivadomed/ivadomed/blob/master/dev/prepare_data/README.md>`_.

    In addition to the MRI data, this sample dataset also includes a trained model for spinal cord segmentation.

    To download the dataset (~490MB), run the following commands in your terminal:

    .. code-block:: bash

       # Download data
       ivadomed_download_data -d data_example_spinegeneric

Configuration file
------------------

    In ``ivadomed``, training is orchestrated by a configuration file. Examples of configuration files are available in
    the ``ivadomed/config/`` and the documentation is available in :doc:`../configuration_file`.

    In this tutorial we will use the configuration file: ``ivadomed/config/config.json``.
    First off, copy this configuration file in your local directory (to avoid modifying the source file):

    .. code-block:: bash

       cp <PATH_TO_IVADOMED>/ivadomed/config/config.json .

    Then, open it with a text editor. Below we will discuss some of the key parameters to perform a one-class 2D
    segmentation training.

    - ``command``: Action to perform. Here, we want to train a model, so we set the fields as follows:

      .. code-block:: json

         "command": "train"

    Note that you can also pass this argument via CLI (see `Usage <../usage.html>`__)

      .. code-block:: bash

        ivadomed --train -c path/to/config

    - ``path_output``: Folder name that will contain the output files (e.g., trained model, predictions, results).

      .. code-block:: json

         "path_output": "spineGeneric"

    Note that you can also pass this argument via CLI (see `Usage <../usage.html>`__)

      .. code-block:: bash

        ivadomed -c path/to/config --path-output path/to/output/directory

    - ``loader_parameters:path_data``: Location of the dataset. As discussed in `Data <../data.html>`__, the dataset
      should conform to the BIDS standard. Modify the path so it points to the location of the downloaded dataset.

      .. code-block:: json

         "path_data": "data_example_spinegeneric"

    Note that you can also pass this argument via CLI (see `Usage <../usage.html>`__)

      .. code-block:: bash

        ivadomed -c path/to/config --path-data path/to/bids/data

    - ``loader_parameters:target_suffix``: Suffix of the ground truth segmentation. The ground truth is located
      under the ``DATASET/derivatives/labels`` folder. In our case, the suffix is ``_seg-manual``:

      .. code-block:: json

         "target_suffix": ["_seg-manual"]

    - ``loader_parameters:contrast_params``: Contrast(s) of interest

      .. code-block:: json

         "contrast_params": {
             "training_validation": ["T1w", "T2w", "T2star"],
             "testing": ["T1w", "T2w", "T2star"],
             "balance": {}
         }

    - ``loader_parameters:slice_axis``: Orientation of the 2D slice to use with the model.

      .. code-block:: json

         "slice_axis": "axial"

    - ``loader_parameters:multichannel``: Turn on/off multi-channel training. If ``true``, each sample has several
      channels, where each channel is an image contrast. If ``false``, only one image contrast is used per sample.

      .. code-block:: json

         "multichannel": false

      .. note::

         The multichannel approach requires that for each subject, the image contrasts are co-registered. This implies that
         a ground truth segmentation is aligned with all contrasts, for a given subject. In this tutorial, only one channel
         will be used.

    - ``training_time:num_epochs``: the maximum number of epochs that will be run during training. Each epoch is composed
      of a training part and an evaluation part. It should be a strictly positive integer.

      .. code-block:: json

         "num_epochs": 100

Train model
-----------

    Once the configuration file is ready, run the training:

    .. code-block:: bash

       ivadomed --train -c config.json --path-data path/to/bids/data --path-output path/to/output/directory

    - We can pass other flags to execute different commands (training, testing, segmentation), see `Usage <../usage.html>`__.


    - ``--path-output``: Folder name that will contain the output files (e.g., trained model, predictions, results).

      .. code-block:: bash

         --path-output path/to/output/directory

    - ``--path-data``: Location of the dataset. As discussed in `Data <../data.html>`__, the dataset
      should conform to the BIDS standard. Modify the path so it points to the location of the downloaded dataset.

      .. code-block:: bash

         --path-data path/to/bids/data

    - If you set the ``command``, ``path_output``, and ``path_data`` arguments in your config file, you do not need to pass the CLI flags:

    .. code-block:: bash

       ivadomed -c config.json

    .. note::

       If a `compatible GPU <https://pytorch.org/get-started/locally/>`_ is available, it will be used by default.
       Otherwise, training will use the CPU, which will take a prohibitively long computational time (several hours).

    The main parameters of the training scheme and model will be displayed on the terminal, followed by the loss value
    on training and validation sets at every epoch. To know more about the meaning of each parameter, go to
    :doc:`../configuration_file`. The value of the loss should decrease during the training.

    .. code-block:: console

       Creating output path: spineGeneric
       Cuda is not available.
       Working on cpu.

       Selected architecture: Unet, with the following parameters:
       dropout_rate: 0.3
       bn_momentum: 0.1
       depth: 3
       is_2d: True
       final_activation: sigmoid
       folder_name: my_model
       in_channel: 1
       out_channel: 1
       Dataframe has been saved in spineGeneric\bids_dataframe.csv.
       After splitting: train, validation and test fractions are respectively 0.6, 0.2 and 0.2 of participant_id.

       Selected transformations for the ['training'] dataset:
       Resample: {'hspace': 0.75, 'wspace': 0.75, 'dspace': 1}
       CenterCrop: {'size': [128, 128]}
       RandomAffine: {'degrees': 5, 'scale': [0.1, 0.1], 'translate': [0.03, 0.03], 'applied_to': ['im', 'gt']}
       ElasticTransform: {'alpha_range': [28.0, 30.0], 'sigma_range': [3.5, 4.5], 'p': 0.1, 'applied_to': ['im', 'gt']}
       NumpyToTensor: {}
       NormalizeInstance: {'applied_to': ['im']}

       Selected transformations for the ['validation'] dataset:
       Resample: {'hspace': 0.75, 'wspace': 0.75, 'dspace': 1}
       CenterCrop: {'size': [128, 128]}
       NumpyToTensor: {}
       NormalizeInstance: {'applied_to': ['im']}
       Loading dataset: 100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 383.65it/s]
       Loaded 92 axial slices for the validation set.
       Loading dataset: 100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 17/17 [00:00<00:00, 282.10it/s]
       Loaded 276 axial slices for the training set.
       Creating model directory: spineGeneric\my_model

       Initialising model's weights from scratch.

       Scheduler parameters: {'name': 'CosineAnnealingLR', 'base_lr': 1e-05, 'max_lr': 0.01}

       Selected Loss: DiceLoss
       with the parameters: []
       Epoch 1 training loss: -0.0336.
       Epoch 1 validation loss: -0.0382.


    After 100 epochs (see ``num_epochs`` in the configuration file), the Dice score on the validation set should
    be ~90%.

.. _Evaluate model:

Evaluate model
--------------

    To test the trained model on the testing sub-dataset and compute evaluation metrics, run:

    .. code-block:: bash

       ivadomed --test -c config.json --path-data path/to/bids/data --path-output path/to/output/directory

    If you prefer to use config files over CLI flags, set ``command`` to the following in you config file:

    .. code-block:: json

       "command": "test"

    You can also set ``path_output``, and ``path_data`` arguments in your config file.

    Then run:

    .. code-block:: bash

       ivadomed -c config.json

    The model's parameters will be displayed in the terminal, followed by a preview of the results for each image.
    The resulting segmentation is saved for each image in the ``<PATH_TO_OUT_DIR>/pred_masks`` while a csv file,
    saved in ``<PATH_TO_OUT_DIR>/results_eval/evaluation_3Dmetrics.csv``, contains all the evaluation metrics. For more details
    on the evaluation metrics, see :mod:`ivadomed.metrics`.

    .. code-block:: console

       Output path already exists: spineGeneric
       Cuda is not available.
       Working on cpu.

       Selected architecture: Unet, with the following parameters:
       dropout_rate: 0.3
       bn_momentum: 0.1
       depth: 3
       is_2d: True
       final_activation: sigmoid
       folder_name: my_model
       in_channel: 1
       out_channel: 1
       Dataframe has been saved in spineGeneric\bids_dataframe.csv.
       After splitting: train, validation and test fractions are respectively 0.6, 0.2 and 0.2 of participant_id.

       Selected transformations for the ['testing'] dataset:
       Resample: {'hspace': 0.75, 'wspace': 0.75, 'dspace': 1}
       CenterCrop: {'size': [128, 128]}
       NumpyToTensor: {}
       NormalizeInstance: {'applied_to': ['im']}
       Loading dataset: 100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 6/6 [00:00<00:00, 373.59it/s]
       Loaded 94 axial slices for the testing set.

       Loading model: spineGeneric\best_model.pt
       Inference - Iteration 0: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 6/6 [00:29<00:00,  4.86s/it]
       {'dice_score': 0.9334570551249012, 'multi_class_dice_score': 0.9334570551249012, 'precision_score': 0.925126264682505, 'recall_score': 0.9428409070673442, 'specificity_score': 0.9999025807354961, 'intersection_over_union': 0.8756498644456311, 'accu
       racy_score': 0.9998261755671077, 'hausdorff_score': 0.05965616760384793}

       Run Evaluation on spineGeneric\pred_masks

       Evaluation: 100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 6/6 [00:05<00:00,  1.04it/s]
                         avd_class0  dice_class0  lfdr_101-INFvox_class0  lfdr_class0  ltpr_101-INFvox_class0  ltpr_class0  mse_class0  ...  n_pred_class0  precision_class0  recall_class0  rvd_class0  specificity_class0  vol_gt_class0  vol_pred_class0
       image_id                                                                                                                            ...
       sub-mpicbs06_T1w       0.086296     0.940116                     0.0          0.0                     1.0          1.0    0.002292  ...            1.0          0.902774       0.980680   -0.086296            0.999879    4852.499537      5271.249497
       sub-mpicbs06_T2star    0.038346     0.909164                     0.0          0.0                     1.0          1.0    0.003195  ...            1.0          0.892377       0.926595   -0.038346            0.999871    4563.749565      4738.749548
       sub-mpicbs06_T2w       0.032715     0.947155                     0.0          0.0                     1.0          1.0    0.001971  ...            1.0          0.932153       0.962648   -0.032715            0.999920    4852.499537      5011.249522
       sub-unf01_T1w          0.020288     0.954007                     0.0          0.0                     1.0          1.0    0.002164  ...            1.0          0.944522       0.963684   -0.020288            0.999917    6161.249412      6286.249400
       sub-unf01_T2star       0.001517     0.935124                     0.0          0.0                     1.0          1.0    0.002831  ...            1.0          0.934416       0.935834   -0.001517            0.999904    5766.249450      5774.999449

       [5 rows x 16 columns]


    The test image segmentations are stored in ``<PATH_TO_OUT_DIR>/pred_masks/`` and have the same name as the input image
    with the suffix ``_pred``. To visualize the segmentation of a given subject, you can use any Nifti image viewer.
    For `FSLeyes <https://open.win.ox.ac.uk/pages/fsl/fsleyes/fsleyes/userdoc/>`_ users, this command will open the
    input image with the overlaid prediction (segmentation) for one of the test subject:

    .. code-block:: bash

       fsleyes <PATH_TO_BIDS_DATA>/sub-mpicbs06/anat/sub-mpicbs06_T2w.nii.gz <PATH_TO_OUT_DIR>/pred_masks/sub-mpicbs06_T2w_pred.nii.gz -cm red

    After the training for 100 epochs, the segmentations should be similar to the one presented in the following image.
    The output and ground truth segmentations of the spinal cord are presented in red (subject ``sub-mpicbs06`` with
    contrast T2w):

    .. image:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/tutorials/one_class_segmentation_2d_unet/sc_prediction.png
       :align: center
Cascaded architecture
=====================

    In this tutorial we will learn the following features:

    - Design a training scheme composed of two cascaded networks.
    - Visualize the training with tensorboard.
    - Generate a GIF to visualize the learning of the model.
    - Find the optimal threshold to binarize images based on the validation sub-dataset.

    In our example, the model will first locate the spinal cord (step 1). This localisation will then be used to crop the images around this region of interest, before segmenting the cerebrospinal fluid (CSF, step 2).

Download dataset
----------------

    A dataset example is available for this tutorial. If not already done, download the dataset with the following line.
    For more details on this dataset see :ref:`One-class segmentation with 2D U-Net<Download dataset>`.

    .. code-block:: bash

       # Download data
       ivadomed_download_data -d data_example_spinegeneric

Configuration file
------------------

    In this tutorial we will use the configuration file: ``ivadomed/config/config.json``.
    First off, copy this configuration file in your local directory to avoid modifying the source file:

    .. code-block:: bash

       cp <PATH_TO_IVADOMED>/ivadomed/config/config.json .

    Then, open it with a text editor. As described in the tutorial :doc: `../tutorials/one_class_segmentation_2d_unet`, make
    sure the ``command`` is set to ``train`` and ``path_data`` point to the location of the dataset. Below, we will discuss
    some of the key parameters to use cascaded models.

    - ``debugging``: Boolean, create extended verbosity and intermediate outputs. Here we will look at the intermediate predictions
      with tensorboard, we therefore need to activate those intermediate outputs.

      .. code-block:: json

         "debugging": true

    - ``object_detection_params:object_detection_path``: Location of the object detection model. This parameter corresponds
      to the path of the first model from the cascaded architecture that segments the spinal cord. The packaged model in the
      downloaded dataset located in the folder `trained_model/seg_sc_t1-t2-t2s-mt` will be used to detect the spinal cord.
      This spinal cord segmentation model will be applied to the images and a bounding box will be created around this mask
      to crop the image.

      .. code-block:: json

         "object_detection_path": "<PATH_TO_DATASET>/trained_model/seg_sc_t1-t2-t2s-mt"

    - ``object_detection_params:safety_factor``: Multiplicative factor to apply to each dimension of the bounding box. To
      ensure all the CSF is included, a safety factor should be applied to the bounding box generated from the spinal cord.
      A safety factor of 200% on each dimension is applied on the height and width of the image. The original depth of the
      bounding box is kept since the CSF should not be present past this border.

      .. code-block:: json

         "safety_factor": [2, 2, 1]

    - ``loader_parameters:target_suffix``: Suffix of the ground truth segmentation. The ground truth are located under the
      ``DATASET/derivatives/labels`` folder. The suffix for CSF labels in this dataset is ``_csfseg-manual``:

      .. code-block:: json

         "target_suffix": ["_csfseg-manual"]

    - ``loader_parameters:contrast_params``: Contrast(s) of interest. The CSF labels are only available in T2w contrast in
      the example dataset.

      .. code-block:: json

         "contrast_params": {
             "training_validation": ["T2w"],
             "testing": ["T2w"],
             "balance": {}
         }

    - ``transformation:CenterCrop:size``: Crop size in voxel. Images will be cropped or padded to fit these dimensions. This
      allows all the images to have the same size during training. Since the images will be cropped around the spinal cord,
      the image size can be reduced to avoid large zero padding.

      .. code-block:: json

         "CenterCrop": {
             "size": [64, 64]
         }

Train model
-----------

    Once the configuration file is ready, run the training. ``ivadomed`` has an option to find a threshold value which optimized the dice score on the validation dataset. This threshold will be further used to binarize the predictions on testing data. Add the flag ``-t`` with an increment
    between 0 and 1 to perform this threshold optimization (i.e. ``-t 0.1`` will return the best threshold between 0.1,
    0.2, ..., 0.9)

    To help visualize the training, the flag ``--gif`` or ``-g`` can be used. The flag should be followed by the number of
    slices by epoch to visualize. For example, ``-g 2`` will generate 2 GIFs of 2 randomly selected slices from the
    validation set.

    Make sure to run the CLI command with the ``--train`` flag, and to point to the location of the dataset via the flag ``--path-data path/to/bids/data``.

    .. code-block:: bash

       ivadomed --train -c config.json --path-data path/to/bids/data --path-output path/to/output/directory -t 0.01 -g 1

    If you prefer to use config files over CLI flags, set ``command`` to the following in you config file:

    .. code-block:: json

       "command": "train"

    You can also set ``path_output``, and ``path_data`` arguments in your config file.

    Then run:

    .. code-block:: bash

       ivadomed -c config.json

    At the end of the training, the optimal threshold will be indicated:

    .. code-block:: console

       Running threshold analysis to find optimal threshold
        Optimal threshold: 0.01
        Saving plot: spineGeneric/roc.png


Visualize training data
-----------------------
    If the flag ``--gif`` or ``-g`` was used, the training can be visualized through gifs located in the folder specified by the ``--path-output`` flag
    ``<PATH_TO_OUT_DIR>/gifs``.

    .. figure:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/tutorials/cascaded_architecture/training.gif
       :width: 300
       :align: center

       Training visualization with GIF

    Another way to visualize the training is to use Tensorboard. Tensorboard helps to visualize the augmented input images,
    the model's prediction, the ground truth, the learning curves, and more. To access this data during or after training,
    use the following command-line:

    .. code-block:: bash

       tensorboard --logdir <PATH_TO_OUT_DIR>

    The following should be displayed in the terminal:

    .. code-block:: console

       Serving TensorBoard on localhost; to expose to the network, use a proxy or pass --bind_all
       TensorBoard 2.2.1 at http://localhost:6006/ (Press CTRL+C to quit)

    Open your browser and type the URL provided, in this case ``http://localhost:6006/``.
    In the scalars folder, the evolution of metrics, learning rate and loss through the epochs can be visualized.

    .. image:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/tutorials/cascaded_architecture/tensorboard_scalar.png
       :align: center

    In the image folder, the training and validation ground truth, input images and predictions are displayed. With this
    feature, it is possible to visualize the cropping from the first model and confirm that the spinal cord
    was correctly located.

    .. image:: https://raw.githubusercontent.com/ivadomed/doc-figures/main/tutorials/cascaded_architecture/tensorboard_images.png
       :align: center

Evaluate model
--------------
    - ``postprocessing:binarize_prediction``: Threshold at which predictions are binarized. Before testing the model,
      modify the binarization threshold to have a threshold adapted to the data:

    .. code-block:: json

        "binarize_prediction": 0.01


    To test and apply this model on the testing dataset, go to the `Evaluate model` section of the tutorial
    :ref:`One-class segmentation with 2D U-Net<evaluate model>`.
