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
