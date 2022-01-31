# Change Log

All notable changes to this project will be documented in this file. It's a team effort
to make them as straightforward as possible.

The format is based on [Keep a Changelog](http://keepachangelog.com/) and this project
adheres to [Semantic Versioning](http://semver.org/).

## [1.0.0-rc1] - In Progress

Release comment: refactoring of models means that old checkpoint files are no longer
compatible with the updates.

### Added

- Added `num_parallel_calls` option in config for data preprocessing.
- Added tests for Dice score, Jaccard Index, and cross entropy losses.
- Added statistics on inputs, DDF and TRE into tensorboard.
- Added example for using custom loss.
- Added tests on Mac OS.
- Added tests for python 3.6 and 3.7.
- Added support to custom layer channels in U-Net.
- Added support to multiple loss functions for each loss type: "image", "label" and
  "regularization".
- Added LNCC computation using separable 1-D filters for all kernels available

### Changed

- Updated pre-trained models for unpaired_ct_abdomen demo to new version
- Changed dataset config so that `format` and `labeled` are defined per split.
- Reduced TensorFlow logging level.
- Used `DEEPREG_LOG_LEVEL` to control logging in DeepReg.
- Increased all EPS to 1e-5.
- Clarify the suggestion in doc to use all-zero masks for missing labels.
- Moved contributor list to a separate page.
- Changed `no-test` flag to `full` for demo scripts.
- Renamed `neg_weight` to `background_weight`.
- Renamed `log_dir` to `exp_name` and `log_root` to `log_dir` respectively.
- Uniformed local-net, global-net, u-net under a single u-net structure.
- Simplified custom layer definitions.
- Removed multiple unnecessary custom layers and use tf.keras.layers whenever possible.
- Refactored BSplines interpolation independently of the backbone network and available
  only for DDF and DVF models.

### Fixed

- Fixed using GPU remotely
- Fixed LNCC loss regarding INF values.
- Removed loss weight checks to be more robust.
- Fixed import error under python 3.6.
- Fixed the residual module in local net architecture, compatible for previous
  checkpoints.
- Broken link in README to seminar video.

## [0.1.2] - 2021-01-31

Release comment: This is mainly a bugfix release, although some of the tasks in
1.0.0-rc1 have been included in this release, with or without public-facing
accessibility (see details below).

### Added

- Added global NCC loss
- Added the docs on registry for backbone models.
- Added backward compatible config parser.
- Added tests so that test coverage is 100%.
- Added config file docs with details on how new config works.
- Added DDF data augmentation.
- Added the registry for backbone models and losses.
- Added pylint with partial check (C0103,C0301,R1725,W0107,W9012,W9015) to CI.
- Added badges for code quality and maintainability.
- Added additional links (CoC, PyPI) and information (contributing, citing) to project
  README.md.
- Added CMIC seminar where DeepReg was introduced to the project README.md.
- Added deepreg_download entry point to access non-release folders required for Quick
  Start.

### Changed

- Refactored optimizer configuration.
- Refactored affine transform data augmentation.
- Modified the implementation of resampler to support zero boundary condition.
- Refactored loss functions into classes.
- Use CheckpointManager callback for saving and support training restore.
- Changed distribute strategy to default for <= 1 GPU.
- Migrated from Travis-CI to GitHub Actions.
- Simplified configuration for backbone models and losses.
- Simplified contributing documentation.
- Uniform kernel size for LNCC loss.
- Improved demo configurations with the updated pre-trained models for:
  grouped_mask_prostate_longitudinal, paried_mrus_prostate, unpaired_us_prostate_cv,
  grouped_mr_heart, unpaired_ct_lung, paired_ct_lung.

### Fixed

- Fixed several dead links in the documentation.
- Fixed a bug due to typo when image loss weight is zero, label loss is not applied.
- Fixed warp CLI tool by saving outputs in Nifti1 format.
- Fixed optimiser storage and loading from checkpoints.
- Fixed bias initialization for theta in GlobalNet.
- Removed invalid `first` argument in DataLoader for sample_index generator.
- Fixed build error when downloading data from the private repository.
- Fixed the typo for CLI tools in documents.

## [0.1.0] - 2020-11-02

### Added

- Added option to change the kernel size and type for LNCC image similarity loss.
- Added visualization tool for generating gifs from model outputs.
- Added the max_epochs argument for training to overwrite configuration.
- Added the log_root argument for training and prediction to customize the log file
  location.
- Added more meaningful error messages for data loading.
- Added integration tests for all demos.
- Added environment.yml file for Conda environment creation.
- Added Dockerfile.
- Added the documentation about using UCL cluster with DeepReg.

### Changed

- Updated TensorFlow version to 2.3.1.
- Updated the pre-trained models in MR brain demo.
- Updated instruction on Conda environment creation.
- Updated the documentation regarding pre-commit and unit-testing.
- Updated the issue and pull-request templates.
- Updated the instructions for all demos.
- Updated pre-commit hooks version.
- Updated JOSS paper to address reviewers' comments.
- Migrated from travis-ci.org to travis-ci.com.

### Fixed

- Fixed prediction error when number of samples cannot be divided by batch size exactly.
- Fixed division by zero handling in multiple image/label losses.
- Fixed tensor comparison in unit tests and impacted tests.
- Removed normalization of DDF/DVF when saving in Nifti formats.
- Fixed invalid link in the quick start page.

## [0.1.0b1] - 2020-09-01

Initial beta release.
<p align="center">
  <img src="https://raw.githubusercontent.com/DeepRegNet/DeepReg/main/docs/asset/deepreg_logo_purple.svg"
    alt="deepreg_logo" title="DeepReg" width="200"/>
</p>

<table align="center">
  <tr>
    <td>
      <b>Package</b>
    </td>
    <td>
      <a href="https://opensource.org/licenses/Apache-2.0">
      <img src="https://img.shields.io/badge/License-Apache%202.0-blue.svg" alt="License">
      </a>
      <a href="https://pypi.python.org/pypi/DeepReg/">
      <img src="https://img.shields.io/pypi/v/deepreg.svg" alt="PyPI Version">
      </a>
      <a href="https://pypi.python.org/pypi/DeepReg/">
      <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/deepreg">
      </a>
      <a href="https://pepy.tech/project/deepreg">
      <img src="https://static.pepy.tech/personalized-badge/deepreg?period=total&units=none&left_color=grey&right_color=orange&left_text=Downloads"
        alt="PyPI downloads">
      </a>
    </td>
  </tr>
  <tr>
    <td>
      <b>Documentation</b>
    </td>
    <td>
      <a href="https://deepreg.readthedocs.io/en/latest/?badge=latest">
      <img src="https://readthedocs.org/projects/deepreg/badge/?version=latest" alt="Documentation Status">
      </a>
    </td>
  </tr>
  <tr>
    <td>
      <b>Code</b>
    </td>
    <td>
      <a href="https://github.com/DeepRegNet/DeepReg/actions?query=workflow%3A%22Unit+Test%22">
      <img src="https://github.com/deepregnet/deepreg/workflows/Unit%20Test/badge.svg?branch=main" alt="Unit Test">
      </a>
      <a href="https://github.com/DeepRegNet/DeepReg/actions?query=workflow%3A%22Integration+Test%22">
      <img src="https://github.com/deepregnet/deepreg/workflows/Integration%20Test/badge.svg?branch=main" alt="Integration Test">
      </a>
      <a href="https://codecov.io/github/DeepRegNet/DeepReg">
      <img src="https://codecov.io/gh/DeepRegNet/DeepReg/branch/main/graph/badge.svg" alt="Coverage Status">
      </a>
      <a href="https://github.com/psf/black">
      <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Code Style">
      </a>
      <a href="https://scrutinizer-ci.com/g/DeepRegNet/DeepReg/">
      <img src="https://scrutinizer-ci.com/g/DeepRegNet/DeepReg/badges/quality-score.png" alt="Code Quality">
      </a>
      <a href="https://codeclimate.com/github/DeepRegNet/DeepReg/maintainability">
      <img src="https://api.codeclimate.com/v1/badges/65245e28aa8f2cd7c6b6/maintainability" alt="Code Maintainability">
      </a>
    </td>
  </tr>
  <tr>
    <td>
      <b>Papers</b>
    </td>
    <td>
      <a href="https://joss.theoj.org/papers/7e6de472bc82a70d7618e23f618960b3"><img
        src="https://joss.theoj.org/papers/7e6de472bc82a70d7618e23f618960b3/status.svg"
        alt="JOSS Paper"></a>
      <a href="https://zenodo.org/badge/latestdoi/269365590"><img src="https://zenodo.org/badge/269365590.svg"
        alt="DOI"></a>
    </td>
  </tr>
</table>

# DeepReg

**DeepReg is a freely available, community-supported open-source toolkit for research
and education in medical image registration using deep learning.**

- TensorFlow 2-based for efficient training and rapid deployment;
- Implementing major unsupervised and weakly-supervised algorithms, with their
  combinations and variants;
- Focusing on growing and diverse clinical applications, with all DeepReg Demos using
  open-accessible data;
- Simple built-in command line tools requiring minimal programming and scripting;
- Open, permissible and research-and-education-driven, under the Apache 2.0 license.

---

## Getting Started

- [DeepReg.Net](http://deepreg.net/)
- [Documentation, Tutorials, and Quick Start](https://deepreg.readthedocs.io/)
- [Medical Image Registration Demos using DeepReg](https://deepreg.readthedocs.io/en/latest/demo/introduction.html)
- [Issue Tracker](https://github.com/DeepRegNet/DeepReg/issues/new/choose)

## Contributing

Get involved, and help make DeepReg better! We want your help - **Really**.

**Being a contributor doesn't just mean writing code.** Equally important to the
open-source process is writing or proof-reading documentation, suggesting or
implementing tests, or giving feedback about the project. You might see the errors and
assumptions that have been glossed over. If you can write any code at all, you can
contribute code to open-source. We are constantly trying out new skills, making
mistakes, and learning from those mistakes. That's how we all improve, and we are happy
to help others learn with us.

### Code of Conduct

This project is released with a
[Code of Conduct](https://github.com/DeepRegNet/DeepReg/blob/main/docs/CODE_OF_CONDUCT.md).
By participating in this project, you agree to abide by its terms.

### Where Should I Start?

For guidance on making a contribution to DeepReg, see our
[Contribution Guidelines](https://deepreg.readthedocs.io/en/latest/contributing/guide.html).

Have a registration application with openly accessible data? Consider
[contributing a DeepReg Demo](https://deepreg.readthedocs.io/en/latest/contributing/demo.html).

## MICCAI 2020 Educational Challenge

Our [MICCAI Educational Challenge](https://miccai-sb.github.io/materials.html)
submission on DeepReg is an Award Winner!

Check it out
[here](https://github.com/DeepRegNet/DeepReg/blob/main/docs/Intro_to_Medical_Image_Registration.ipynb) -
you can also
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/DeepRegNet/DeepReg/blob/main/docs/Intro_to_Medical_Image_Registration.ipynb)

## Overview Video

Members of the DeepReg dev team presented "The Road to DeepReg" at the Centre for
Medical Imaging Computing (CMIC) seminar series at University College London on the 4th
of November 2020. You can access the talk
[here](https://www.youtube.com/watch?v=jDEyWXZM3CE&feature=youtu.be).

## Citing DeepReg

DeepReg is research software, made by a
[team of academic researchers](https://deepreg.readthedocs.io/en/latest/#contributors).
Citations and use of our software help us justify the effort which has gone into, and
will keep going into, maintaining and growing this project.

If you have used DeepReg in your research, please consider citing us:

> Fu _et al._, (2020). DeepReg: a deep learning toolkit for medical image registration.
> _Journal of Open Source Software_, **5**(55), 2705,
> https://doi.org/10.21105/joss.02705

Or with BibTex:

```
@article{Fu2020,
  doi = {10.21105/joss.02705},
  url = {https://doi.org/10.21105/joss.02705},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {55},
  pages = {2705},
  author = {Yunguan Fu and Nina Montaña Brown and Shaheer U. Saeed and Adrià Casamitjana and Zachary M. C. Baum and Rémi Delaunay and Qianye Yang and Alexander Grimwood and Zhe Min and Stefano B. Blumberg and Juan Eugenio Iglesias and Dean C. Barratt and Ester Bonmati and Daniel C. Alexander and Matthew J. Clarkson and Tom Vercauteren and Yipeng Hu},
  title = {DeepReg: a deep learning toolkit for medical image registration},
  journal = {Journal of Open Source Software}
}
```
# Description

Please include a summary of the change and which issue is fixed. Please also include
relevant motivation and context. List any dependencies that are required for this
change.

Fixes #<issue_number>

## Type of change

What types of changes does your code introduce to DeepReg?

_Please check the boxes that apply after submitting the pull request._

- [ ] Bugfix (non-breaking change which fixes an issue)
- [ ] Code style update (formatting, renaming)
- [ ] Refactoring (no functional changes, no api changes)
- [ ] Documentation Update (fix or improvement on the documentation)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Other (if none of the other choices apply)

## Checklist

_Please check the boxes that apply after submitting the pull request._

_If you're unsure about any of them, don't hesitate to ask. We're here to help! This is
simply a reminder of what we are going to look for before merging your code._

- [ ] I have
      [installed pre-commit](https://deepreg.readthedocs.io/en/latest/contributing/setup.html)
      using `pre-commit install` and formatted all changed files. If you are not
      certain, run `pre-commit run --all-files`.
- [ ] My commits' message styles matches
      [our requested structure](https://deepreg.readthedocs.io/en/latest/contributing/commit.html),
      e.g. `Issue #<issue number>: detailed message`.
- [ ] I have updated the
      [change log file](https://github.com/DeepRegNet/DeepReg/blob/main/CHANGELOG.md)
      regarding my changes.
- [ ] I have added tests that prove my fix is effective or that my feature works
- [ ] I have added necessary documentation (if appropriate)
---
name: 💡 Feature request
about: It would be nice to have some new feature.
---

## Subject of the feature

What are you trying to do and how would you want to do it differently?

Is it something you currently you cannot do?

Is this related to an issue/problem?

Has the feature been requested before? If yes, please provide a link to the issue.

If the feature request is approved, would you be willing to submit a PR? _(Help can be
provided if you need assistance submitting a PR)_

Yes / No
---
name: 😀 Other
about: Anything not covered by previous options.
---

## Subject of the issue

Please describe the issue, it can be a question or any other discussions.
---
name: 📝 Documentation request
about: Some documentation is missing or incorrect.
---

## Subject of the documentation

What are you trying to do and which part of the code is hard to understand or incorrect
or missing? Please provide a link to the readthedocs page and code if available.

If the documentation request is approved, would you be willing to submit a PR? _(Help
can be provided if you need assistance submitting a PR)_

Yes / No
---
name: 🐜 Bug report
about: Something isn't working.
---

## Subject of the issue

Describe your issue here.

If the bug is confirmed, would you be willing to submit a PR? _(Help can be provided if
you need assistance submitting a PR)_

Yes / No

## Your environment

- DeepReg version (commit hash)

  Please use `git rev-parse HEAD` to get the hash of the current commit. Using
  `pip list` will provide the fixed tag version inside `setup.py`, therefore it is not
  accurate.

  We recommend installing DeepReg using
  [Anaconda](https://docs.anaconda.com/anaconda/install/) /
  [Miniconda](https://docs.conda.io/en/latest/miniconda.html) in a separate virtual
  environment.

- OS (e.g. Ubuntu 20.04, MacOS 10.15, etc.)

  We do not officially support windows.

- Python Version (3.7, 3.8, etc.)

  We support only Python 3.7 officially.

- TensorFlow

  - TensorFlow Version (2.2, 2.3, etc.)
  - CUDA Version (10.1, etc.) if available.
  - cuDNN Version

  We support only TensorFlow 2.3 officially.

  If using GPU, please check https://www.tensorflow.org/install/source#gpu to verify the
  GPU support.

## Steps to reproduce

Tell us how to reproduce this issue. Please provide a working demo.

## Expected behaviour

Tell us what should happen.

## Actual behaviour

Tell us what happens instead.
---
name: 🧪 Test request
about: Some code is not covered by test.
---

## Subject of the code / test

What are you trying to do and which part of the code is not tested? Please provide a
link to the code if available.

If the test request is approved, would you be willing to submit a PR? _(Help can be
provided if you need assistance submitting a PR)_

Yes / No
# DeepReg Demos

This folder contains multiple self-contained registration demos for real-world clinical
scenarios. Each demo contains data preparation, training, and inference scripts.
# Classical affine registration for head-and-neck CT images

> **Note**: Please read the
> [DeepReg Demo Disclaimer](introduction.html#demo-disclaimer).

[Source Code](https://github.com/DeepRegNet/DeepReg/tree/main/demos/classical_ct_headneck_affine)

This is a special demo that uses the DeepReg package for classical affine image
registration, which iteratively solves an optimisation problem. Gradient descent is used
to minimise the image dissimilarity function of a given pair of moving anf fixed images.

## Author

DeepReg Development Team

## Application

Although in this demo the moving images are simulated using a randomly generated
transformation. The registration technique can be used in radiotherapy to compensate the
difference between CT acquired at different time points, such as pre-treatment and
intra-/post-treatment.

## Data

Data is
[an example CT volume](https://wiki.cancerimagingarchive.net/display/Public/Head-Neck-PET-CT)
with two labels.

## Instruction

### Installation

Please install DeepReg following the [instructions](../getting_started/install.html) and
change the current directory to the root directory of DeepReg project, i.e. `DeepReg/`.

### Download data

Please execute the following command to download and pre-process the data.

```bash
python demos/classical_ct_headneck_affine/demo_data.py
```

### Launch registration

Please execute the following command to register two images. The fixed image will be the
downloaded data and the moving image will be simulated by applying a random affine
transformation, such that the ground-truth is available for. The optimised
transformation will be applied to the moving images, as well as the moving labels. The
results, saved in a timestamped folder under the project directory, will compare the
warped image/labels with the ground-truth image/labels.

```bash
python demos/classical_ct_headneck_affine/demo_register.py
```

## Visualise

The following command can be executed to generate a plot of three image slices from the
the moving image, warped image and fixed image (left to right) to visualise the
registration. Please see the visualisation tool docs
[here](https://github.com/DeepRegNet/DeepReg/blob/main/docs/source/docs/visualisation_tool.md)
for more visualisation options such as animated gifs.

```bash
deepreg_vis -m 2 -i 'demos/classical_ct_headneck_affine/logs_reg/moving_image.nii.gz, demos/classical_ct_headneck_affine/logs_reg/warped_moving_image.nii.gz, demos/classical_ct_headneck_affine/logs_reg/fixed_image.nii.gz' --slice-inds '4,8,12' -s demos/classical_ct_headneck_affine/logs_reg
```

Note: The registration script must be run before running the command to generate the
visualisation.

![plot](../assets/classical_ct_headneck_affine.png)

## Contact

Please [raise an issue](https://github.com/DeepRegNet/DeepReg/issues/new/choose) for any
questions.

## Reference

[1] Vallières, M. et al. Radiomics strategies for risk assessment of tumour failure in
head-and-neck cancer. Sci Rep 7, 10117 (2017). doi: 10.1038/s41598-017-10371-5
# Unpaired lung CT registration

> **Note**: Please read the
> [DeepReg Demo Disclaimer](introduction.html#demo-disclaimer).

[Source Code](https://github.com/DeepRegNet/DeepReg/tree/main/demos/unpaired_ct_lung)

## Author

DeepReg Development Team (Shaheer Saeed)

## Application

This is a registration between CT images from different patients. The images are all
from acquired at the same timepoint in the breathing cycle. This is an inter subject
registration. This kind of registration is useful for determining how one stimulus
affects multiple patients. If a drug or invasive procedure is administered to multiple
patients, registering the images from different patients can give medical professionals
a sense of how each patient is responding in comparison to others. An example of such an
application can be seen in [2].

## Data

The dataset for this demo comes from [1] and can be downloaded from:
https://zenodo.org/record/3835682#.XsUWXsBpFhE

## Instruction

### Installation

Please install DeepReg following the [instructions](../getting_started/install.html) and
change the current directory to the root directory of DeepReg project, i.e. `DeepReg/`.

### Download data

Please execute the following command to download/pre-process the data and download the
pre-trained model. Image intensities are rescaled during pre-processing.

```bash
python demos/unpaired_ct_lung/demo_data.py
```

### Launch demo training

Please execute the following command to launch a demo training. The training logs and
model checkpoints will be saved under `demos/unpaired_ct_lung/logs_train`.

```bash
python demos/unpaired_ct_lung/demo_train.py
```

Here the training is launched using the GPU of index 0 with a limited number of steps
and reduced size. Please add flag `--full` to use the original training configuration,
such as

```bash
python demos/unpaired_ct_lung/demo_train.py --full
```

### Predict

Please execute the following command to run the prediction with pre-trained model. The
prediction logs and visualization results will be saved under
`demos/unpaired_ct_lung/logs_predict`. Check the [CLI documentation](../docs/cli.html)
for more details about prediction output.

```bash
python demos/unpaired_ct_lung/demo_predict.py
```

Optionally, the user-trained model can be used by changing the `ckpt_path` variable
inside `demo_predict.py`. Note that the path should end with `.ckpt` and checkpoints are
saved under `logs_train` as mentioned above.

## Visualise

The following command can be executed to generate a plot of three image slices from the
the moving image, warped image and fixed image (left to right) to visualise the
registration. Please see the visualisation tool docs
[here](https://github.com/DeepRegNet/DeepReg/blob/main/docs/source/docs/visualisation_tool.md)
for more visualisation options such as animated gifs.

```bash
deepreg_vis -m 2 -i 'demos/unpaired_ct_lung/logs_predict/<time-stamp>/test/<pair-number>/moving_image.nii.gz, demos/unpaired_ct_lung/logs_predict/<time-stamp>/test/<pair-number>/pred_fixed_image.nii.gz, demos/unpaired_ct_lung/logs_predict/<time-stamp>/test/<pair-number>/fixed_image.nii.gz' --slice-inds '40,48,56' -s demos/unpaired_ct_lung/logs_predict
```

Note: The prediction must be run before running the command to generate the
visualisation. The `<time-stamp>` and `<pair-number>` must be entered by the user.

![plot](../assets/unpaired_ct_lung.png)

## Contact

Please [raise an issue](https://github.com/DeepRegNet/DeepReg/issues/new/choose) for any
questions.

## Reference

[1] Hering A, Murphy K, and van Ginneken B. (2020). Lean2Reg Challenge: CT Lung
Registration - Training Data [Data set]. Zenodo. http://doi.org/10.5281/zenodo.3835682

[2] Li B, Christensen GE, Hoffman EA, McLennan G, Reinhardt JM. Establishing a normative
atlas of the human lung: intersubject warping and registration of volumetric CT images.
Acad Radiol. 2003;10(3):255-265. doi:10.1016/s1076-6332(03)80099-5
# Unpaired prostate ultrasound registration

> **Note**: Please read the
> [DeepReg Demo Disclaimer](introduction.html#demo-disclaimer).

[Source Code](https://github.com/DeepRegNet/DeepReg/tree/main/demos/unpaired_us_prostate_cv)

This DeepReg Demo is also an example of cross validation.

## Author

DeepReg Development Team

## Application

Transrectal ultrasound (TRUS) images are acquired from prostate cancer patients during
image-guided procedures. Pairwise registration between these 3D images may be useful for
intraoperative motion modelling and group-wise registration for population studies.

## Data

The 3D ultrasound images used in this demo were derived from the Prostate-MRI-US-Biopsy
dataset, hosted at the
[Cancer Imaging Archive (TCIA)](https://www.cancerimagingarchive.net/).

## Instruction

### Installation

Please install DeepReg following the [instructions](../getting_started/install.html) and
change the current directory to the root directory of DeepReg project, i.e. `DeepReg/`.

### Download data

Please execute the following command to download/pre-process the data and download the
pre-trained model. Data are split into 10 folds for cross-validation.

```bash
python demos/unpaired_us_prostate_cv/demo_data.py
```

### Launch demo training

Please execute the following command to launch a demo training (the first of the ten
runs of a 9-fold cross-validation). The training logs and model checkpoints will be
saved under `demos/unpaired_us_prostate_cv/logs_train`.

```bash
python demos/unpaired_us_prostate_cv/demo_train.py
```

Here the training is launched using the GPU of index 0 with a limited number of steps
and reduced size. Please add flag `--full` to use the original training configuration,
such as

```bash
python demos/unpaired_us_prostate_cv/demo_train.py --full
```

### Predict

Please execute the following command to run the prediction with pre-trained model. The
prediction logs and visualization results will be saved under
`demos/unpaired_us_prostate_cv/logs_predict`. Check the
[CLI documentation](../docs/cli.html) for more details about prediction output.

```bash
python demos/unpaired_us_prostate_cv/demo_predict.py
```

Optionally, the user-trained model can be used by changing the `ckpt_path` variable
inside `demo_predict.py`. Note that the path should end with `.ckpt` and checkpoints are
saved under `logs_train` as mentioned above.

## Visualise

The following command can be executed to generate a plot of three image slices from the
the moving image, warped image and fixed image (left to right) to visualise the
registration. Please see the visualisation tool docs
[here](https://github.com/DeepRegNet/DeepReg/blob/main/docs/source/docs/visualisation_tool.md)
for more visualisation options such as animated gifs.

```bash
deepreg_vis -m 2 -i 'demos/unpaired_us_prostate_cv/logs_predict/<time-stamp>/test/<pair-number>/moving_image.nii.gz, demos/unpaired_us_prostate_cv/logs_predict/<time-stamp>/test/<pair-number>/pred_fixed_image.nii.gz, demos/unpaired_us_prostate_cv/logs_predict/<time-stamp>/test/<pair-number>/fixed_image.nii.gz' --slice-inds '50,65,35' -s demos/unpaired_us_prostate_cv/logs_predict/
```

Note: The prediction must be run before running the command to generate the
visualisation. The `<time-stamp>` and `<pair-number>` must be entered by the user.

![plot](../assets/unpaired_us_prostate_cv.png)

## Contact

Please [raise an issue](https://github.com/DeepRegNet/DeepReg/issues/new/choose) for any
questions.
# Paired prostate MR-ultrasound registration

> **Note**: Please read the
> [DeepReg Demo Disclaimer](introduction.html#demo-disclaimer).

> **Warning**:
> [This demo ought to be improved in the future.](https://github.com/DeepRegNet/DeepReg/issues/621).

[Source Code](https://github.com/DeepRegNet/DeepReg/tree/main/demos/paired_mrus_brain)

This demo uses DeepReg to re-implement the algorithms described in
[Weakly-supervised convolutional neural networks for multimodal image registration](https://doi.org/10.1016/j.media.2018.07.002).
A standalone demo was hosted at https://github.com/yipenghu/label-reg.

## Author

DeepReg Development Team

## Application

Registering preoperative MR images to intraoperative transrectal ultrasound images has
been an active research area for more than a decade. The multimodal image registration
task assist a number of ultrasound-guided interventions and surgical procedures, such as
targeted biopsy and focal therapy for prostate cancer patients. One of the key
challenges in this registration task is the lack of robust and effective similarity
measures between the two image types. This demo implements a weakly-supervised learning
approach to learn voxel correspondence between intensity patterns between the multimodal
data, driven by expert-defined anatomical landmarks, such as the prostate gland
segmentaion.

## Data

This is a demo without real clinical data due to regulatory restrictions. The MR and
ultrasound images used are simulated dummy images.

## Instruction

### Installation

Please install DeepReg following the [instructions](../getting_started/install.html) and
change the current directory to the root directory of DeepReg project, i.e. `DeepReg/`.

### Download data

Please execute the following command to download/pre-process the data and download the
pre-trained model.

```bash
python demos/paired_mrus_prostate/demo_data.py
```

### Launch demo training

Please execute the following command to launch a demo training (the first of the ten
runs of a 9-fold cross-validation). The training logs and model checkpoints will be
saved under `demos/paired_mrus_prostate/logs_train`.

```bash
python demos/paired_mrus_prostate/demo_train.py
```

Here the training is launched using the GPU of index 0 with a limited number of steps
and reduced size. Please add flag `--full` to use the original training configuration,
such as

```bash
python demos/paired_mrus_prostate/demo_train.py --full
```

### Predict

Please execute the following command to run the prediction with pre-trained model. The
prediction logs and visualization results will be saved under
`demos/paired_mrus_prostate/logs_predict`. Check the
[CLI documentation](../docs/cli.html) for more details about prediction output.

```bash
python demos/paired_mrus_prostate/demo_predict.py
```

Optionally, the user-trained model can be used by changing the `ckpt_path` variable
inside `demo_predict.py`. Note that the path should end with `.ckpt` and checkpoints are
saved under `logs_train` as mentioned above.

## Visualise

The following command can be executed to generate a plot of three image slices from the
the moving image, warped image and fixed image (left to right) to visualise the
registration. Please see the visualisation tool docs
[here](https://github.com/DeepRegNet/DeepReg/blob/main/docs/source/docs/visualisation_tool.md)
for more visualisation options such as animated gifs.

```bash
deepreg_vis -m 2 -i 'demos/paired_mrus_prostate/logs_predict/<time-stamp>/test/<pair-number>/moving_image.nii.gz, demos/paired_mrus_prostate/logs_predict/<time-stamp>/test/<pair-number>/pred_fixed_image.nii.gz, demos/paired_mrus_prostate/logs_predict/<time-stamp>/test/<pair-number>/fixed_image.nii.gz' --slice-inds '12,20,36' -s demos/paired_mrus_prostate/logs_predict
```

Note: The prediction must be run before running the command to generate the
visualisation. The `<time-stamp>` and `<pair-number>` must be entered by the user.

![plot](../assets/paired_mrus_prostate.png)

## Contact

Please [raise an issue](https://github.com/DeepRegNet/DeepReg/issues/new/choose) for any
questions.
# Unpaired hippocampus MR registration

> **Note**: Please read the
> [DeepReg Demo Disclaimer](introduction.html#demo-disclaimer).

> **Warning**:
> [This demo ought to be improved in the future.](https://github.com/DeepRegNet/DeepReg/issues/620).

[Source Code](https://github.com/DeepRegNet/DeepReg/tree/main/demos/unpaired_mr_brain)

## Author

DeepReg Development Team (Adrià Casamitjana)

## Application

This is a demo targeting the alignment of hippocampal substructures (head and body)
using mono-modal MR images between different patients. The images are cropped around
those areas and manually annotated. This is a 3D intra-modal registration using a
composite loss of image and label similarity.

## Data

The dataset for this demo comes from the Learn2Reg MICCAI Challenge (Task 4) [1] and can
be downloaded from:
https://drive.google.com/uc?export=download&id=1RvJIjG2loU8uGkWzUuGjqVcGQW2RzNYA

## Instruction

### Installation

Please install DeepReg following the [instructions](../getting_started/install.html) and
change the current directory to the root directory of DeepReg project, i.e. `DeepReg/`.

### Download data

Please execute the following command to download/pre-process the data and download the
pre-trained model.

```bash
python demos/unpaired_mr_brain/demo_data.py
```

Pre-processing includes:

- Rescaling all images' intensity to 0-255.
- Creating and applying a binary mask to mask-out the padded values in images.
- Transforming label volumes using one-hot encoding (only for foreground classes)

### Launch demo training

Please execute the following command to launch a demo training. The training logs and
model checkpoints will be saved under `demos/unpaired_mr_brain/logs_train`.

```bash
python demos/unpaired_mr_brain/demo_train.py
```

Here the training is launched using the GPU of index 0 with a limited number of steps
and reduced size. Please add flag `--full` to use the original training configuration,
such as

```bash
python demos/unpaired_mr_brain/demo_train.py --full
```

### Predict

Please execute the following command to run the prediction with pre-trained model. The
prediction logs and visualization results will be saved under
`demos/unpaired_mr_brain/logs_predict`. Check the [CLI documentation](../docs/cli.html)
for more details about prediction output.

```bash
python demos/unpaired_mr_brain/demo_predict.py
```

## Visualise

The following command can be executed to generate a plot of three image slices from the
the moving image, warped image and fixed image (left to right) to visualise the
registration. Please see the visualisation tool docs
[here](https://github.com/DeepRegNet/DeepReg/blob/main/docs/source/docs/visualisation_tool.md)
for more visualisation options such as animated gifs.

```bash
deepreg_vis -m 2 -i 'demos/unpaired_mr_brain/logs_predict/<time-stamp>/test/<pair-number>/moving_image.nii.gz, demos/unpaired_mr_brain/logs_predict/<time-stamp>/test/<pair-number>/pred_fixed_image.nii.gz, demos/unpaired_mr_brain/logs_predict/<time-stamp>/test/<pair-number>/fixed_image.nii.gz' --slice-inds '20,32,44' -s demos/unpaired_mr_brain/logs_predict/
```

Note: The prediction must be run before running the command to generate the
visualisation. The `<time-stamp>` and `<pair-number>` must be entered by the user.

![plot](../assets/unpaired_mr_brain.png)

## Contact

Please [raise an issue](https://github.com/DeepRegNet/DeepReg/issues/new/choose) for any
questions.

## Reference

[1] AL Simpson et al., _A large annotated medical image dataset for the development and
evaluation of segmentation algorithms_ (2019). https://arxiv.org/abs/1902.09063
# Unpaired abdomen CT registration

> **Note**: Please read the
> [DeepReg Demo Disclaimer](introduction.html#demo-disclaimer).

> **Warning**:
> [This demo ought to be improved in the future.](https://github.com/DeepRegNet/DeepReg/issues/552).

[Source Code](https://github.com/DeepRegNet/DeepReg/tree/main/demos/unpaired_ct_abdomen)

## Author

DeepReg Development Team (Ester Bonmati)

## Application

This demo shows how to register unpaired abdominal CT data from different patients using
DeepReg. In addition, the demo demonstrates the difference between the unsupervised,
weakly-supervised and their combination, using a U-Net.

## Data

The data set is from the MICCAI Learn2Reg grand challenge
(https://learn2reg.grand-challenge.org/) task 3 [1], and can be downloaded directly from
https://learn2reg.grand-challenge.org/Datasets/.

## Instruction

### Installation

Please install DeepReg following the [instructions](../getting_started/install.html) and
change the current directory to the root directory of DeepReg project, i.e. `DeepReg/`.

### Download data

Please execute the following command to download/pre-process the data and download the
pre-trained model.

```bash
python demos/unpaired_ct_abdomen/demo_data.py
```

### Launch demo training

In this demo, three different training methods are provided: unsupervised, weakly
supervised and the combined method. Please execute one of the following commands to
launch a demo training. The training logs and model checkpoints will be saved under
`demos/unpaired_ct_abdomen/logs_train/method` with `method` be `unsup`, `weakly` or
`comb`.

```bash
python demos/unpaired_ct_abdomen/demo_train.py --method unsup
python demos/unpaired_ct_abdomen/demo_train.py --method weakly
python demos/unpaired_ct_abdomen/demo_train.py --method comb
```

Here the training is launched using the GPU of index 0 with a limited number of steps
and reduced size. Please add flag `--full` to use the original training configuration,
such as

```bash
python demos/unpaired_ct_abdomen/demo_train.py --method unsup --full
```

### Predict

Please execute one of the following commands to run the prediction with pre-trained
model. The prediction logs and visualization results will be saved under
`demos/unpaired_ct_abdomen/logs_predict/method` with `method` be `unsup`, `weakly` or
`comb`. Check the [CLI documentation](../docs/cli.html) for more details about
prediction output.

```bash
python demos/unpaired_ct_abdomen/demo_predict.py --method unsup
python demos/unpaired_ct_abdomen/demo_predict.py --method weakly
python demos/unpaired_ct_abdomen/demo_predict.py --method comb
```

Optionally, the user-trained model can be used by changing the `ckpt_path` variable
inside `demo_predict.py`. Note that the path should end with `.ckpt` and checkpoints are
saved under `logs_train` as mentioned above.

## Visualise

The following command can be executed to generate a plot of three image slices from the
the moving image, warped image and fixed image (left to right) to visualise the
registration. Please see the visualisation tool docs
[here](https://github.com/DeepRegNet/DeepReg/blob/main/docs/source/docs/visualisation_tool.md)
for more visualisation options such as animated gifs.

```bash
deepreg_vis -m 2 -i 'demos/unpaired_ct_abdomen/logs_predict/comb/<time-stamp>/test/<pair-number>/moving_image.nii.gz, demos/unpaired_ct_abdomen/logs_predict/comb/<time-stamp>/test/<pair-number>/pred_fixed_image.nii.gz, demos/unpaired_ct_abdomen/logs_predict/comb/<time-stamp>/test/<pair-number>/fixed_image.nii.gz' --slice-inds '30,50,65' -s demos/unpaired_ct_abdomen/logs_predict
```

Note: The prediction must be run before running the command to generate the
visualisation. The `<time-stamp>` and `<pair-number>` must be entered by the user.

![plot](../assets/unpaired_ct_abdomen.png)

## Contact

Please [raise an issue](https://github.com/DeepRegNet/DeepReg/issues/new/choose) for any
questions.

## Reference

[1] Adrian Dalca, Yipeng Hu, Tom Vercauteren, Mattias Heinrich, Lasse Hansen, Marc
Modat, Bob de Vos, Yiming Xiao, Hassan Rivaz, Matthieu Chabanas, Ingerid Reinertsen,
Bennett Landman, Jorge Cardoso, Bram van Ginneken, Alessa Hering, and Keelin Murphy.
(2020, March 19). Learn2Reg - The Challenge. Zenodo.
http://doi.org/10.5281/zenodo.3715652
# Pairwise registration for grouped prostate segmentation masks

> **Note**: Please read the
> [DeepReg Demo Disclaimer](introduction.html#demo-disclaimer).

[Source Code](https://github.com/DeepRegNet/DeepReg/tree/main/demos/grouped_mask_prostate_longitudinal)

This demo uses DeepReg to demonstrate a number of features:

- For grouped data in h5 files, e.g. "group-1-2" indicates the 2th visit from Subject 1;
- Use masks as the images for feature-based registration - aligning the prostate gland
  segmentation in this case - with deep learning;
- Register intra-patient longitudinal data.

This demo also implements the feature-based registration described in
[Morphological Change Forecasting for Prostate Glands using Feature-based Registration and Kernel Density Extrapolation](https://arxiv.org/abs/2101.06425).

## Author

DeepReg Development Team

## Application

Longitudinal registration detects the temporal changes and normalises the spatial
difference between images acquired at different time-points. For prostate cancer
patients under active surveillance programmes, quantifying these changes is useful for
detecting and monitoring potential cancerous regions.

## Data

This is a demo without real clinical data due to regulatory restrictions. The MR and
ultrasound images used are simulated dummy 3D Ultrasound images. Data are organized into
10 separate folds.

## Instruction

### Installation

Please install DeepReg following the [instructions](../getting_started/install.html) and
change the current directory to the root directory of DeepReg project, i.e. `DeepReg/`.

### Download data

Please execute the following command to download/pre-process the data and download the
pre-trained model.

```bash
python demos/grouped_mask_prostate_longitudinal/demo_data.py
```

### Launch demo training

Please execute the following command to launch a demo training (the first of the ten
runs of a 9-fold cross-validation). The training logs and model checkpoints will be
saved under `demos/grouped_mask_prostate_longitudinal/logs_train`.

```bash
python demos/grouped_mask_prostate_longitudinal/demo_train.py
```

Here the training is launched using the GPU of index 0 with a limited number of steps
and reduced size. Please add flag `--full` to use the original training configuration,
such as

```bash
python demos/grouped_mask_prostate_longitudinal/demo_train.py --full
```

### Predict

Please execute the following command to run the prediction with pre-trained model. The
prediction logs and visualization results will be saved under
`demos/grouped_mask_prostate_longitudinal/logs_predict`. Check the
[CLI documentation](../docs/cli.html) for more details about prediction output.

```bash
python demos/grouped_mask_prostate_longitudinal/demo_predict.py
```

Optionally, the user-trained model can be used by changing the `ckpt_path` variable
inside `demo_predict.py`. Note that the path should end with `.ckpt` and checkpoints are
saved under `logs_train` as mentioned above.

## Visualise

The following command can be executed to generate a plot of three image slices from the
the moving image, warped image and fixed image (left to right) to visualise the
registration. Please see the visualisation tool docs
[here](https://github.com/DeepRegNet/DeepReg/blob/main/docs/source/docs/visualisation_tool.md)
for more visualisation options such as animated gifs.

```bash
deepreg_vis -m 2 -i 'demos/grouped_mask_prostate_longitudinal/logs_predict/<time-stamp>/test/<pair-number>/moving_image.nii.gz, demos/grouped_mask_prostate_longitudinal/logs_predict/<time-stamp>/test/<pair-number>/pred_fixed_image.nii.gz, demos/grouped_mask_prostate_longitudinal/logs_predict/<time-stamp>/test/<pair-number>/fixed_image.nii.gz' --slice-inds '10,16,20' -s demos/grouped_mask_prostate_longitudinal/logs_predict
```

Note: The prediction must be run before running the command to generate the
visualisation. The `<time-stamp>` and `<pair-number>` must be entered by the user.

![plot](../assets/grouped_mask_prostate_longitudinal.png)

## Contact

Please [raise an issue](https://github.com/DeepRegNet/DeepReg/issues/new/choose) for any
questions.
# Pairwise registration for grouped cardiac MR images

> **Note**: Please read the
> [DeepReg Demo Disclaimer](introduction.html#demo-disclaimer).

[Source Code](https://github.com/DeepRegNet/DeepReg/tree/main/demos/grouped_mr_heart)

This demo uses the grouped dataset loader to register intra-subject multi-sequence
cardiac magnetic resonance (CMR) images.

## Author

DeepReg Development Team

## Application

Computer-assisted management for patients suffering from myocardial infraction (MI)
often requires quantifying the difference and comprising the multiple sequences, such as
the late gadolinium enhancement (LGE) CMR sequence MI, the T2-weighted CMR. They
collectively provide radiological information otherwise unavailable during clinical
practice.

## Data

This demo uses CMR images from 45 patients, acquired from the
[MyoPS2020](http://www.sdspeople.fudan.edu.cn/zhuangxiahai/0/MyoPS20/) challenge held in
conjunction with MICCAI 2020.

## Instruction

### Installation

Please install DeepReg following the [instructions](../getting_started/install.html) and
change the current directory to the root directory of DeepReg project, i.e. `DeepReg/`.

### Download data

Please execute the following command to download/pre-process the data and download the
pre-trained model. Images are re-sampled to an isotropic voxel size.

```bash
python demos/grouped_mr_heart/demo_data.py
```

### Launch demo training

Please execute the following command to launch a demo training (the first of the ten
runs of a 9-fold cross-validation). The training logs and model checkpoints will be
saved under `demos/grouped_mr_heart/logs_train`.

```bash
python demos/grouped_mr_heart/demo_train.py
```

Here the training is launched using the GPU of index 0 with a limited number of steps
and reduced size. Please add flag `--full` to use the original training configuration,
such as

```bash
python demos/grouped_mr_heart/demo_train.py --full
```

### Predict

Please execute the following command to run the prediction with pre-trained model. The
prediction logs and visualization results will be saved under
`demos/grouped_mr_heart/logs_predict`. Check the [CLI documentation](../docs/cli.html)
for more details about prediction output.

```bash
python demos/grouped_mr_heart/demo_predict.py
```

Optionally, the user-trained model can be used by changing the `ckpt_path` variable
inside `demo_predict.py`. Note that the path should end with `.ckpt` and checkpoints are
saved under `logs_train` as mentioned above.

## Visualise

The following command can be executed to generate a plot of three image slices from the
the moving image, warped image and fixed image (left to right) to visualise the
registration. Please see the visualisation tool docs
[here](https://github.com/DeepRegNet/DeepReg/blob/main/docs/source/docs/visualisation_tool.md)
for more visualisation options such as animated gifs.

```bash
deepreg_vis -m 2 -i 'demos/grouped_mr_heart/logs_predict/<time-stamp>/test/<pair-number>/moving_image.nii.gz, demos/grouped_mr_heart/logs_predict/<time-stamp>/test/<pair-number>/pred_fixed_image.nii.gz, demos/grouped_mr_heart/logs_predict/<time-stamp>/test/<pair-number>/fixed_image.nii.gz' --slice-inds '14,10,20' -s demos/grouped_mr_heart/logs_predict
```

Note: The prediction must be run before running the command to generate the
visualisation. The `<time-stamp>` and `<pair-number>` must be entered by the user.

![plot](../assets/grouped_mr_heart.png)

## Contact

Please [raise an issue](https://github.com/DeepRegNet/DeepReg/issues/new/choose) for any
questions.

## Reference

[1] Xiahai Zhuang: Multivariate mixture model for myocardial segmentation combining
multi-source images. IEEE Transactions on Pattern Analysis and Machine Intelligence (T
PAMI), vol. 41, no. 12, 2933-2946, Dec 2019. link.

[2] Xiahai Zhuang: Multivariate mixture model for cardiac segmentation from
multi-sequence MRI. International Conference on Medical Image Computing and
Computer-Assisted Intervention, pp.581-588, 2016.
# Paired brain MR-ultrasound registration

> **Note**: Please read the
> [DeepReg Demo Disclaimer](introduction.html#demo-disclaimer).

> **Warning**:
> [This demo ought to be improved in the future.](https://github.com/DeepRegNet/DeepReg/issues/620).

[Source Code](https://github.com/DeepRegNet/DeepReg/tree/main/demos/paired_mrus_brain)

## Author

DeepReg Development Team (Shaheer Saeed)

## Application

This demo aims to register pairs of brain MR and ultrasound scans. The dataset consists
of 22 subjects with low-grade brain gliomas who underwent brain tumour resection [1].
The main application for this type of registration is to better delineate brain tumour
boundaries during surgery and correct tissue shift induced by the craniotomy.

## Data

The dataset for this demo comes from Xiao et al. [1] and can be downloaded from:
https://archive.sigma2.no/pages/public/datasetDetail.jsf?id=10.11582/2020.00025.

## Instruction

Please install DeepReg following the [instructions](../getting_started/install.html) and
change the current directory to the root directory of DeepReg project, i.e. `DeepReg/`.

### Download data

Please execute the following command to download/pre-process the data and download the
pre-trained model. By default, the downloaded data is only a partial of the original
one. However the access to the original data is temporarily unavailable.

```bash
python demos/paired_mrus_brain/demo_data.py
```

### Launch demo training

Please execute the following command to launch a demo training. The training logs and
model checkpoints will be saved under `demos/paired_mrus_brain/logs_train`.

```bash
python demos/paired_mrus_brain/demo_train.py
```

Here the training is launched using the GPU of index 0 with a limited number of steps
and reduced size. Please add flag `--full` to use the original training configuration,
such as

```bash
python demos/paired_mrus_brain/demo_train.py --full
```

Note: The number of epochs and reduced dataset size for training will result in a loss
in test accuracy so please train with the full dataset and for a greater number of
epochs for improved results.

### Predict

Please execute the following command to run the prediction with pre-trained model. The
prediction logs and visualization results will be saved under
`demos/paired_mrus_brain/logs_predict`. Check the [CLI documentation](../docs/cli.html)
for more details about prediction output.

```bash
python demos/paired_mrus_brain/demo_predict.py
```

Optionally, the user-trained model can be used by changing the `ckpt_path` variable
inside `demo_predict.py`. Note that the path should end with `.ckpt` and checkpoints are
saved under `logs_train` as mentioned above.

## Visualise

The following command can be executed to generate a plot of three image slices from the
the moving image, warped image and fixed image (left to right) to visualise the
registration. Please see the visualisation tool docs
[here](https://github.com/DeepRegNet/DeepReg/blob/main/docs/source/docs/visualisation_tool.md)
for more visualisation options such as animated gifs.

```bash
deepreg_vis -m 2 -i 'demos/paired_mrus_brain/logs_predict/<time-stamp>/test/<pair-number>/moving_image.nii.gz, demos/paired_mrus_brain/logs_predict/<time-stamp>/test/<pair-number>/pred_fixed_image.nii.gz, demos/paired_mrus_brain/logs_predict/<time-stamp>/test/<pair-number>/fixed_image.nii.gz' --slice-inds '190,128,96' -s demos/paired_mrus_brain/logs_predict
```

Note: The prediction must be run before running the command to generate the
visualisation. The `<time-stamp>` and `<pair-number>` must be entered by the user.

![plot](../assets/paired_mrus_brain.png)

## Contact

Please [raise an issue](https://github.com/DeepRegNet/DeepReg/issues/new/choose) for any
questions.

## Reference

[1] Y. Xiao, M. Fortin, G. Unsgård , H. Rivaz, and I. Reinertsen, "REtroSpective
Evaluation of Cerebral Tumors (RESECT): a clinical database of pre-operative MRI and
intra-operative ultrasound in low-grade glioma surgeries". Medical Physics, Vol. 44(7),
pp. 3875-3882, 2017.
# Classical nonrigid registration for prostate MR images

> **Note**: Please read the
> [DeepReg Demo Disclaimer](introduction.html#demo-disclaimer).

[Source Code](https://github.com/DeepRegNet/DeepReg/tree/main/demos/classical_mr_prostate_nonrigid)

This is a special demo that uses the DeepReg package for classical nonrigid image
registration, which iteratively solves an optimisation problem. Gradient descent is used
to minimise the image dissimilarity function of a given pair of moving and fixed images,
often regularised by a deformation smoothness function.

## Author

DeepReg Development Team

## Application

Registering inter-subject prostate MR images may be useful to align different glands in
a common space for investigating the spatial distribution of cancer.

## Data

Data is an example MR volumes with the prostate gland segmentation from
[MICCAI Grand Challenge: Prostate MR Image Segmentation 2012](https://promise12.grand-challenge.org/).

## Instruction

### Installation

Please install DeepReg following the [instructions](../getting_started/install.html) and
change the current directory to the root directory of DeepReg project, i.e. `DeepReg/`.

### Download data

Please execute the following command to download and pre-process the data.

```bash
python demos/classical_mr_prostate_nonrigid/demo_data.py
```

### Launch registration

Please execute the following command to register two images. The optimised
transformation will be applied to the moving images, as well as the moving labels. The
results, saved in a timestamped folder under the project directory, will compare the
warped image/labels with the ground-truth image/labels.

```bash
python demos/classical_mr_prostate_nonrigid/demo_register.py
```

## Visualise

The following command can be executed to generate a plot of three image slices from the
the moving image, warped image and fixed image (left to right) to visualise the
registration. Please see the visualisation tool docs
[here](https://github.com/DeepRegNet/DeepReg/blob/main/docs/source/docs/visualisation_tool.md)
for more visualisation options such as animated gifs.

```bash
deepreg_vis -m 2 -i 'demos/classical_mr_prostate_nonrigid/logs_reg/moving_image.nii.gz, demos/classical_mr_prostate_nonrigid/logs_reg/warped_moving_image.nii.gz, demos/classical_mr_prostate_nonrigid/logs_reg/fixed_image.nii.gz' --slice-inds '4,8,12' -s demos/classical_mr_prostate_nonrigid/logs_reg
```

Note: The registration script must be run before running the command to generate the
visualisation.

![plot](../assets/classical_mr_prostate_nonrigid.png)

## Contact

Please [raise an issue](https://github.com/DeepRegNet/DeepReg/issues/new/choose) for any
questions.

## Reference

[1] Litjens, G., Toth, R., van de Ven, W., Hoeks, C., Kerkstra, S., van Ginneken, B.,
Vincent, G., Guillard, G., Birbeck, N., Zhang, J. and Strand, R., 2014. Evaluation of
prostate segmentation algorithms for MRI: the PROMISE12 challenge. Medical image
analysis, 18(2), pp.359-373.
# Paired lung CT registration

> **Note**: Please read the
> [DeepReg Demo Disclaimer](introduction.html#demo-disclaimer).

[Source Code](https://github.com/DeepRegNet/DeepReg/tree/main/demos/paired_ct_lung)

## Author

DeepReg Development Team (Shaheer Saeed)

## Application

This is a registration between CT images acquired at different time points for a single
patient. The images being registered are taken at inspiration and expiration for each
subject. This is an intra subject registration. This type of intra subject registration
is useful when there is a need to track certain features on a medical image such as
tumor location when conducting invasive procedures.

## Data

The dataset for this demo comes from
[Lean2Reg Challenge: CT Lung Registration - Training Data](https://zenodo.org/record/3835682#.XsUWXsBpFhE)
[1].

## Instruction

### Installation

Please install DeepReg following the [instructions](../getting_started/install.html) and
change the current directory to the root directory of DeepReg project, i.e. `DeepReg/`.

### Download data

Please execute the following command to download/pre-process the data and download the
pre-trained model. Image intensities are rescaled during pre-processing.

```bash
python demos/paired_ct_lung/demo_data.py
```

### Launch demo training

Please execute the following command to launch a demo training. The training logs and
model checkpoints will be saved under `demos/paired_ct_lung/logs_train`.

```bash
python demos/paired_ct_lung/demo_train.py
```

Here the training is launched using the GPU of index 0 with a limited number of steps
and reduced size. Please add flag `--full` to use the original training configuration,
such as

```bash
python demos/paired_ct_lung/demo_train.py --full
```

### Predict

Please execute the following command to run the prediction with pre-trained model. The
prediction logs and visualization results will be saved under
`demos/paired_ct_lung/logs_predict`. Check the [CLI documentation](../docs/cli.html) for
more details about prediction output.

```bash
python demos/paired_ct_lung/demo_predict.py
```

Optionally, the user-trained model can be used by changing the `ckpt_path` variable
inside `demo_predict.py`. Note that the path should end with `.ckpt` and checkpoints are
saved under `logs_train` as mentioned above.

## Visualise

The following command can be executed to generate a plot of three image slices from the
the moving image, warped image and fixed image (left to right) to visualise the
registration. Please see the visualisation tool docs
[here](https://github.com/DeepRegNet/DeepReg/blob/main/docs/source/docs/visualisation_tool.md)
for more visualisation options such as animated gifs.

```bash
deepreg_vis -m 2 -i 'demos/paired_ct_lung/logs_predict/<time-stamp>/test/<pair-number>/moving_image.nii.gz, demos/paired_ct_lung/logs_predict/<time-stamp>/test/<pair-number>/pred_fixed_image.nii.gz, demos/paired_ct_lung/logs_predict/<time-stamp>/test/<pair-number>/fixed_image.nii.gz' --slice-inds '64,50,72' -s demos/paired_ct_lung/logs_predict/
```

Note: The prediction must be run before running the command to generate the
visualisation. The `<time-stamp>` and `<pair-number>` must be entered by the user.

![plot](../assets/paired_ct_lung.png)

## Contact

Please [raise an issue](https://github.com/DeepRegNet/DeepReg/issues/new/choose) for any
questions.

## Reference

[1] Hering, Alessa, Murphy,Keelin, and van Ginneken, Bram. (2020). Lean2Reg Challenge:
CT Lung Registration : CT Lung Registration - Training Data
# DeepReg Examples

This folder contains multiple self-contained example scripts for customizing DeepReg
functionalities. The purpose is mainly to showcase some functionalities such as using
custom backbones or loss functions.
# DeepReg Contributor List

_(in unsorted order)_

DeepReg is a community-supported project, this list may not be up-to-date and only
includes people who have accepted assignments.

| Name                  | Affiliation (at time of contribution)                             |
| --------------------- | ----------------------------------------------------------------- |
| Adrià Casamitjana     | University College London                                         |
| Alexander Grimwood    | University College London                                         |
| Daniel C. Alexander   | University College London                                         |
| Dean C. Barratt       | University College London                                         |
| Ester Bonmati         | University College London                                         |
| Juan Eugenio Iglesias | University College London / Massachusetts Institute of Technology |
| Matt J. Clarkson      | University College London                                         |
| Nina Montaña Brown    | University College London                                         |
| Qianye Yang           | University College London                                         |
| Rémi Delaunay         | University College London / King’s College London                 |
| Shaheer U. Saeed      | University College London                                         |
| Stefano B. Blumberg   | University College London                                         |
| Tom Vercauteren       | King’s College London                                             |
| Yipeng Hu             | University College London                                         |
| Yiwen Li              | University of Oxford                                              |
| Yunguan Fu            | University College London / InstaDeep                             |
| Zachary M. C. Baum    | University College London                                         |
| Zhe Min               | University College London                                         |
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our community a
harassment-free experience for everyone, regardless of age, body size, visible or
invisible disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal appearance,
race, religion, or sexual identity and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming, diverse,
inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our community
include:

- Demonstrating empathy and kindness toward other people
- Being respectful of differing opinions, viewpoints, and experiences
- Giving and gracefully accepting constructive feedback
- Accepting responsibility and apologizing to those affected by our mistakes, and
  learning from the experience
- Focusing on what is best not just for us as individuals, but for the overall community

Examples of unacceptable behavior include:

- The use of sexualized language or imagery, and sexual attention or advances of any
  kind
- Trolling, insulting or derogatory comments, and personal or political attacks
- Public or private harassment
- Publishing others' private information, such as a physical or email address, without
  their explicit permission
- Other conduct which could reasonably be considered inappropriate in a professional
  setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in response to
any behavior that they deem inappropriate, threatening, offensive, or harmful.

Community leaders have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this
Code of Conduct, and will communicate reasons for moderation decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when an
individual is officially representing the community in public spaces. Examples of
representing our community include using an official e-mail address, posting via an
official social media account, or acting as an appointed representative at an online or
offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported to
the community leaders responsible for enforcement at DeepRegNet@gmail.com. All
complaints will be reviewed and investigated promptly and fairly - community leaders
will aim for a first response to complaints under 48h. The community leaders will aim to
respond within one week to the person who filed the report with either a resolution or
an explanation as to why the situation has not been resolved.

All community leaders are obligated to respect the privacy and security of the reporter
of any incident.

## Conflicts of Interest If a complaint concerns a member of the community leaders and

should the reporter not feel comfortable sharing the report with such member, reports
can be made privately to any of the members of the community leaders by contacting them
privately. Contact methods are listed at
https://github.com/DeepRegNet/DeepReg/wiki/Contact. The enforcement process will be
followed, but will exclude the member or members that the report concerns from
discussion and decision making.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining the
consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing clarity
around the nature of the violation and an explanation of why the behavior was
inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series of actions.

**Consequence**: A warning with consequences for continued behavior. No interaction with
the people involved, including unsolicited interaction with those enforcing the Code of
Conduct, for a specified period of time. This includes avoiding interactions in
community spaces as well as external channels like social media. Violating these terms
may lead to a temporary or permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including sustained
inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public communication
with the community for a specified period of time. No public or private interaction with
the people involved, including unsolicited interaction with those enforcing the Code of
Conduct, is allowed during this period. Violating these terms may lead to a permanent
ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community standards,
including sustained inappropriate behavior, harassment of an individual, or aggression
toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within the
community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 2.0,
available at https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by
[Mozilla's code of conduct enforcement ladder](https://github.com/mozilla/diversity) and
[Django's conduct reporting processes](https://www.djangoproject.com/conduct/reporting/).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
# DeepReg Demo README Example

> **Note**: Please read the
> [DeepReg Demo Disclaimer](introduction.html#demo-disclaimer).

[Source Code](https://github.com/DeepRegNet/DeepReg/tree/main/demos/)

## Author

DeepReg Development Team (optional corresponding Author Name) or Author Name for
external contributor.

## Application

Briefly describe the clinical application and the need for registration.

## Data

Describe the source of the dataset. Dataset should be open-accessible.

## Instruction

### Installation

Please install DeepReg following the [instructions](../getting_started/install.html) and
change the current directory to the root directory of DeepReg project, i.e. `DeepReg/`.

### Download data

Please execute the following command to download/pre-process the data and download the
pre-trained model.

```bash
python demos/name/demo_data.py
```

### Launch demo training

Please execute the following command to launch a demo training. The training logs and
model checkpoints will be saved under `demos/name/logs_train`.

```bash
python demos/name/demo_train.py
```

Here the training is launched using the GPU of index 0 with a limited number of steps
and reduced size. Please add flag `--full` to use the original training configuration,
such as

```bash
python demos/name/demo_train.py --full
```

### Predict

Please execute the following command to run the prediction with pre-trained model. The
prediction logs and visualization results will be saved under `demos/name/logs_predict`.
Check the [CLI documentation](../docs/cli.html) for more details about prediction
output.

```bash
python demos/name/demo_predict.py
```

Optionally, the user-trained model can be used by changing the `ckpt_path` variable
inside `demo_predict.py`. Note that the path should end with `.ckpt` and checkpoints are
saved under `logs_train` as mentioned above.

## Contact

Please [raise an issue](https://github.com/DeepRegNet/DeepReg/issues/new/choose) for any
questions.

## Reference

Related references in the form of `[1] Reference`.
# Design Experiments

DeepReg dataset loaders use a folder/directory-based file storing approach, with which
the user will be responsible for
[organising image and label files in required file formats and folders](../docs/dataset_loader.html).
This design was primarily motivated by the need to minimise the risk of data leakage (or
information leakage), both in code development and subsequent applications.

## Random-split

Every call of the `deepreg_train` or `deepreg_predict` function uses a dataset
"physically" separated by folders, including 'train', 'val' and 'test' sets used in a
random-split experiment. In this case, the user needs to randomly assign available
experiment image and label files into the three folders. Again, for more details see the
[Dataset loader](../docs/dataset_loader.html).

## Cross-validation

Experiments such as _cross-validation_ can be readily implemented by using the
"multi-folder support" in the `dataset` section of the yaml configuration files. See
details in [configuration](../docs/configuration.html).

For example, in a 3-fold cross-validation, the user may randomly partition available
experiment data files into four folders, 'fold0', 'fold1', 'fold2' and 'test'. The
'test' is a hold-out testing set. Each run of the 3-fold cross-validation then can be
specified in a different yaml file as follows.

"cv_run1.yaml":

```yaml
dataset:
  dir:
    train: # training data set
      - "data/test/h5/paired/fold0"
      - "data/test/h5/paired/fold1"
    valid: "data/test/h5/paired/fold2" # validation data set
    test: ""
```

"cv_run2.yaml":

```yaml
dataset:
  dir:
    train: # training data set
      - "data/test/h5/paired/fold0"
      - "data/test/h5/paired/fold2"
    valid: "data/test/h5/paired/fold1" # validation data set
    test: ""
```

"cv_run3.yaml":

```yaml
dataset:
  dir:
    train: # training data set
      - "data/test/h5/paired/fold1"
      - "data/test/h5/paired/fold2"
    valid: "data/test/h5/paired/fold0" # validation data set
    test: ""
```

To further facilitate flexible uses of these dataset loaders, the `deepreg_train` and
`deepreg_predict` functions also accept multiple yaml files - therefore the same `train`
section does not have to be repeated multiple times for the multiple cross-validation
folds or for the test. An example `dataset` section for configuring testing when using
`deepreg_predict` is given below.

"test.yaml":

```yaml
dataset:
  dir:
    train: ""
    valid: ""
    test: "data/test/h5/paired/test" # validation data set
```
# Running DeepReg Remotely (on the Cluster)

This tutorial gives an example of how to run DeepReg remotely (e.g. on a cluster). Our
example is specific to the UCL cluster, which has the operating system CentOS 7 (similar
to Ubuntu), with job scheduler Sun Grid Engine (SGE). More information on the specific
configuration at UCL is available [here](https://hpc.cs.ucl.ac.uk/job-submission/).

## Installing the Environment

Below is the script to install the environment for DeepReg in the cluster. If you want
to switch to the current working branch. Call `git branch origin/<branch_name>` after
downloading DeepReg.

```
git clone https://github.com/<personal_acount_id>/DeepReg.git
module load default/python/3.8.5

cd <DeepReg_Dir>

export PATH=/share/apps/anaconda3-5/bin:$PATH
conda env create -f environment.yml   #set up environment
source /share/apps/source_files/cuda/cuda-10.1.source    # set up cuda for GPU
source activate deepreg   # activate conda env

export CONDA_PIP="/home/<cs_account_id>/.conda/envs/deepreg/bin/pip"
$CONDA_PIP install -e .
```

`module` is a command to find packages and set up them in your own storge. You can get
more information by `module help`. Another way to install packages is

```
export PATH=/share/apps/<module_name>/bin:$PATH
```

You can find any available packages in cluster nodes by

```
ls /share/apps/ | grep '<package_name>*'
```

**Tip:** For now, all the packages are stored in `share/apps/`. If the path does not
exist, try `module load default/python/3.8.5`. Then, call `$PATH` to find the new
location of packages.

## Example Script

Below is the submission script for running quick start example
[here](../getting_started/quick_start.md). Change `<DeepReg_dir>` in the script to the
remote DeepReg repo location and save the below code in a `<your_name>.qsub`. Submit the
job with `qsub <your_name>.qsub` and check the status of the job with `qstat`, the saved
stdout and stderr is in `home/<cs_account_id>/logs/`.

```
#!/bin/bash
#$ -S /bin/bash   # bash for job
#$ -l gpu=true   # use gpu
#$ -l tmem=10G   # virtual mem used
#$ -l h_rt=36:0:0   # max job runtime hour:min:sec
#$ -N DeepReg_tst   # job name
#$ -wd home/<cs_account_id>/logs   # output, error log dir.
#Please call `mkdir logs` before using the script.

hostname
date

cd ../<DeepReg_dir>
export PATH=/share/apps/anaconda3-5/bin:$PATH
source activate deepreg   # activate conda env
export PATH=/share/apps/cuda-10.1/bin:/share/apps/gcc-8.3/bin:$PATH   # path for cuda, gcc
export LD_LIBRARY_PATH=/share/apps/cuda-10.1/lib64:/share/apps/gcc-8.3/lib64:$LD_LIBRARY_PATH   # path for cuda, gcc

deepreg_train \
--gpu \
--config_path config/unpaired_labeled_ddf.yaml \
--log_dir test
```

You also can directly access one of four cluster nodes reserved for development
purposes, by the command below. You can then run your code via the command line. More
information on the specific configuration at UCL is available
[here](https://hpc.cs.ucl.ac.uk/job-submission/).

```
qrsh -l tmem=14G,h_vmem=14G
```

## Contact and Version

Please contact stefano.blumberg.17@ucl.ac.uk, for information about the cluster. The
information here is likely to change as UCL updates the software on the cluster.
# Image Registration with Deep Learning

A series of scientific tutorials on deep learning for registration can be found at the
[learn2reg tutorial](https://learn2reg.github.io/), held in conjunction with
MICCAI 2019.

This document provides a practical overview for a number of algorithms supported by
DeepReg.

## Registration

Image registration is the process of mapping the coordinate system of one image into
another image. A registration method takes a pair of images as input, denoted as moving
and fixed images. In this tutorial, we register the moving image into the fixed image,
i.e. mapping the coordinates of the moving image onto the fixed image.

## Network

### Predict a dense displacement field

With deep learning, given a pair of moving and fixed images, the registration network
outputs a dense displacement field (DDF) with the same shape as the moving image. Each
value can be considered as the placement of the corresponding pixel / voxel of the
moving image. Therefore, the DDF defines a mapping from the moving image's coordinates
to the fixed image.

In this tutorial, we mainly focus on DDF-based methods.

### Predict a dense velocity field

Another option is to predict a dense (static) velocity field (DVF), such that a
diffeomorphic DDF can be numerically integrated. Read
["A fast diffeomorphic image registration algorithm"](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.474.1033&rep=rep1&type=pdf)
and
["Diffeomorphic demons: Efficient non-parametric image registration"](http://www-sop.inria.fr/asclepios/Publications/Tom.Vercauteren/DiffeoDemons-NeuroImage08-Vercauteren.pdf)
for more details.

### Predict an affine transformation

A more constrained option is to predict an affine transformation, parameterised by the
affine transformation matrix of 12 degrees of freedom. The DDF can then be computed to
resample the moving images in fixed image space.

### Predict a region of interest

Instead of outputting the transformation between coordinates, given moving image, fixed
image, and a region of interest (ROI) in the moving image, the network can predict the
ROI in the fixed image directly. Interested readers are referred to the MICCAI 2019
paper:
[Conditional segmentation in lieu of image registration](https://arxiv.org/abs/1907.00438)

## Loss

A loss function has to be defined to train a deep neural network. There are mainly three
types of losses:

### Intensity based (image based) loss

The common loss functions are normalized cross correlation (NCC), sum of squared
distance (SSD), and normalized mutual information (MI).

Intensity based losses measure the dissimilarity between a fixed image and a warped
moving image, which is an adaptation from classical image registration methods. These
losses can perform poorly on multi-modality registrations (e.g. SSD loss in CT-MRI
registration), although certain multi-modality tasks, such as registration between
different MRI sequences are handled well by MI. Generally, intensity based losses work
best when there is an inherent consistency in appearance between moving and fixed
images, which is more common in single-modality registration.

### Feature based (label based) loss

This type of loss measures the dissimilarity of the fixed image labels and warped moving
image labels. The label is often an ROI in the image, like the segmentation of an organ
in a CT image.

The common loss function is Dice loss, Jacard and average cross-entropy over all voxels.

### Deformation loss

This type of loss measures the amount of deformation in an image and penalises
non-smooth deformations . Penalising sudden or discontinuous deformations helps to
regularise the transformation between fixed and moving images.

For DDF, the common loss functions are bending energy, L1 or L2 norm of the displacement
gradient.

## Learning

Depending on the availability of the data labels, registration networks can be trained
with different approaches:

### Unsupervised

When the data label is unavailable, the training can be driven by the unsupervised loss.
The loss function often consists of the intensity based loss and deformation loss. The
following is an illustration of an unsupervised DDF-based registration network.

![Unsupervised DDF-based registration network](../_images/registration-ddf-nn-unsupervised.svg ":size=600")

### Weakly-supervised

When there is no intensity based loss that is appropriate for the image pair one would
like to register, the training can take a pair of corresponding moving and fixed labels
(in addition to the image pair), represented by binary masks, to compute a label
dissimilarity (feature based loss) to drive the registration.

Combined with the regularisation on the predicted displacement field, this forms a
weakly-supervised training. An illustration of an weakly-supervised DDF-based
registration network is provided below.

When multiple labels are available for each image, the labels can be sampled during the
training iteration, such that only one label per image is used in each iteration of the
data set (epoch).

![Weakly-supervised DDF-based registration network](../_images/registration-ddf-nn-weakly-supervised.svg ":size=600")

### Combined

When the data label is available, combining intensity based, feature based, and
deformation based losses together has shown superior registration accuracy, compared to
unsupervised and weakly supervised methods. Following is an illustration of a combined
DDF-based registration network.

![Combined DDF-based registration network](../_images/registration-ddf-nn-combined.svg ":size=600")
# DeepReg demo

The [demos](https://github.com/DeepRegNet/DeepReg/tree/main/demos) folder directly under
the DeepReg root directory contains demonstrations using DeepReg for different image
registration applications.

Contributions are welcome! Below is a set of requirements for a demo to be included as a
DeepReg Demo.

- Each demo _must_ have an independent folder directly under `demos/`;
- Name the folder as
  `[loader-type]_[image-modality]_[organ-disease]_[optional:brief-remark]`, e.g.
  `unpaired_ultrasound_prostate` or `grouped_mr_brain_longitudinal`;
- For simplicity, avoid sub-folders (other than those specified below) and separate
  files for additional functions/classes;
- Experiment using cross-validation or advanced data set sampling is NOT encouraged,
  unless the purpose of the demo is to demonstrate how to design experiments.

## Open accessible data

- Each demo _must_ have a `demo_data.py` script to automatically download and preprocess
  demo data;
- Data for training and test should be downloaded under the demo folder named `dataset`,
  such as `dataset/train` and `dataset/test`;
- Data should be hosted in a reliable and efficient (DeepReg repo will not store demo
  data or model) online storage, Kaggle, GitHub and Zendoo are all options for non-login
  access (avoid google drive for known accessibility issues);
- Relevant dataset folder structure to utilise the supported loaders can be either
  pre-arranged in data source or scripted in `demo_data.py` after downloading;
- Avoid slow and excessively large data set download. Consider downloading a subset as
  default for demonstration purpose, with options for full data set.

## Pre-trained model

- A pre-trained model _must_ be available for downloading, with
  github.com/DeepRegNet/deepreg-model-zoo being preferred for storing the models. Please
  contact the Development Team for access;
- The pre-trained model, e.g. ckpt files, should be downloaded and extracted under the
  `dataset/pretrained` folder. Avoid overwriting with user-trained models;

## Training

- Each demo _must_ have a `demo_train.py` script;
- This is accompanied by one or more config yaml files in the same folder. Please use
  the same demo folder name for the config file. Add postfix if multiple training
  methods are provided, e.g. `unpaired_ct_abdomen_comb.yaml`,
  `unpaired_ct_abdomen_unsup.yaml`.

## Predicting

- Each demo _must_ have a `demo_predict.py` script;
- By default, the pre-trained model should be used in `demo_predict.py`. However, the
  instruction should be clearly given to use the user-trained model, saved with the
  `demo_train.py`;

## A README.md file

A markdown file _must_ be provided under `demos/<demo_name>` as an entry point for each
demo, which should be based on the [template](../demo/readme_template.html). Moreover, a
`.rst` file _must_ be provided under `docs/source/demo` to link the markdown file to the
documentation page. The
[introduction.rst](https://github.com/DeepRegNet/DeepReg/blob/main/docs/source/demo/introduction.rst)
file should be updated properly as well.

Following is a checklist for modifying the README template:

- Update the link to source code;
- Update the author section;
- Update the application section;
- Update the data section, optionally, describe the used pre-processing methods;
- Update the `name` in all commands;
- Update the reference section.
- Optionally, adapt the file to custom needs.
# Documentation

We use [Sphinx](https://www.sphinx-doc.org/en/master/) to organize the documentation and
it is hosted in [ReadTheDocs](https://readthedocs.org/projects/deepreg/).

## Local Build

Please run the following command under `docs/` directory for generating the
documentation pages locally. Generated files are under `docs/build/html/`

```bash
make clean html
```

where

- `clean` removes the possible built files.
- Optionally we can add `SPHINXOPTS="-W"` to fail on any warnings, but we are currently
  having `document isn't included in any toctree` warning and no better solution has
  been found yet.

## Recommendations

There are some recommendations regarding the docs.

- **We prefer markdown files** over reStructuredText files as its linting is covered
  using [Prettier](https://prettier.io/).

  Only use reStructuredText (rst) files for some functionalities not supported by
  markdown, such as

  - `toctree`
  - warning/notes boxes in [Installation](../getting_started/install.html)

  The conversion between markdown and rst can be done automatically using free online
  tool [Pandoc](https://pandoc.org/try/).

- When linking to other pages, **please use relative paths** such as
  `../getting_started/install.html` instead of absolute paths
  `https://deepreg.readthedocs.io/en/latest/getting_started/install.html` as relative
  paths are more robust for different version of documentations.

- To refer a markdown file outside of the source folder, create an rst file and use
  `.. mdinclude:: <makrdown file path>` to include the markdown source.

  Check the source code of
  [paired lung CT image registration page](../demo/paired_ct_lung.html) as an example.
# Release

DeepReg is distributed on PyPI. To create new releases, you can follow the below
instructions to automatically submit new versions to PyPI via our GitHub Actions
workflow.

## Creating a new Release

From the main DeepReg repository page, head to
[releases](https://github.com/DeepRegNet/DeepReg/releases). From here, you can
[draft a new release](https://github.com/DeepRegNet/DeepReg/releases/new).

## Tagging & Titling

We follow [semver](https://semver.org/) naming conventions for tags, with `vX.Y.Z` where
each represents major, minor, and patch release versions.

From semver.org:

> Major version: when you make incompatible API changes,
>
> Minor version: when you add functionality in a backwards compatible manner, and
>
> Patch version: when you make backwards compatible bug fixes.

Typically, most releases will be an increment of the minor or patch versions.

Enter the new version in the format `vX.Y.Z` into the "Tag version" and "Release title"
fields of the draft a release page.

## Publish!

Click the "Publish release" button, and our GitHub Action workflow will handle the rest.

It's recommended that you check the
[output log](https://github.com/DeepRegNet/DeepReg/actions) from the GitHub Actions page
to make sure everything went as planned for the release.
# Guidelines

We welcome contributions to DeepReg. Please
[raise an issue](https://github.com/DeepRegNet/DeepReg/issues/new/choose) to report
bugs, request features, or ask questions. For code contribution, please follow the
guidelines below.

## Setup

We recommend using conda environment on Linux or Mac machines for code development. The
setup steps are:

1. [Install git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).
2. [Clone or fork](https://github.com/DeepRegNet/DeepReg) the repository.
3. [Install](../getting_started/install.html) and activate `deepreg` conda environment.
4. Run `pre-commit install` under the root of this repository `DeepReg/` to install
   pre-commit hooks.

## Resolve an issue

For resolving an issue, please

1. Create a branch.

   The branch name should start with the issue number, followed by hyphen separated
   words describing the issue, e.g. `1-update-contribution-guidelines`.

2. Implement the required features or fix the reported bugs.

   There are several guidelines for commit, coding, testing, and documentation.

   1. Please create commits with meaningful commit messages.

      The commit message should start with `Issue #<issue number>:`, for instance
      `Issue #1: add commit requirements.`

   2. Please write or update unit-tests for the added or changed functionalities.

      Pull request will not be approved if test coverage is decreased. Check
      [testing guidelines](test.html) for further details.

   3. Please write meaningful docstring and documentations for added or changed code.

      We use
      [Sphinx docstring format](https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html).

   4. Please update the
      [CHANGELOG](https://github.com/DeepRegNet/DeepReg/blob/main/CHANGELOG.md)
      regarding the changes.

3. [Create a pull request](https://github.com/DeepRegNet/DeepReg/pulls) when the branch
   is ready.

   Please resume the changes with some more details in the description and check the
   boxes after submitting the pull request. Optionally, you can create a pull request
   and add `WIP` in the name to mark is as working in progress.

## Add a DeepReg demo

Adding DeepReg demo should be done via pull request. Besides the guidelines above,
adding demo has additional requirements described in [demo guidelines](demo.html).
# Unit Test

In DeepReg, we use [pytest](https://docs.pytest.org/en/stable/) (not
[unittest](https://docs.python.org/3/library/unittest.html)) for unit tests to ensure a
certain code quality and to facilitate the code maintenance.

For testing the code locally,

- If you only need to test with python 3.7, please execute `pytest` at the repository
  root:

  ```bash
  pytest test/
  ```

- If you want to test with all supported python versions, please execute `tox` at the
  repository root:

  ```bash
  tox
  ```

Moreover, the testing is checked automatically via
[GitHub workflows](https://github.com/DeepRegNet/DeepReg/actions) and
[Codecov](https://codecov.io/gh/DeepRegNet/DeepReg) is used to monitor the test
coverage. While checking the Codecov report in file mode, generally a line highlighted
by red means it is not covered by test. Please check the
[Codecov documentation](https://docs.codecov.io/docs/viewing-source-code) for more
details.

## Test requirement

We would like to achieve 100% test coverage. In general, tests should be

- thorough, covering different scenarios.
- independent, different scenarios are not tested together.
- clean and compact, for instance,
  - Use [parameterized test](https://docs.pytest.org/en/stable/example/parametrize.html)
    to reduce code redundancy.
  - Use `is_equal_np` and `is_equal_tf` provided in
    [test/unit/util.py](https://github.com/DeepRegNet/DeepReg/blob/main/test/unit/util.py)
    to compare arrays or tensors.

The detailed requirements are as follows:

- Test all functions in python classes.
- Test the trigger of all warning and errors in functions.
- Test the correctness of output values for functions.
- Test at least the correctness of output shapes for TensorFlow functions.

## Example unit test

We provide here an example to help understanding the requirements.

```python
import pytest


def subtract(x: int) -> int:
    """
    A function subtracts one from a non-negative integer.
    :param x: a non-negative integer
    :return: x - 1
    """
    assert isinstance(x, int), f"input {x} is not int"
    assert x >= 0, f"input {x} is negative"
    return x - 1


class TestSubtract:
    @pytest.mark.parametrize("x,expected",[(0, -1), (1,0)])
    def test_value(self, x, expected):
        got = subtract(x=x)
        assert got == expected

    @pytest.mark.parametrize("x,msg", [(-1, "is negative"), (0.0, "is not int")])
    def test_err(self, x, msg):
        with pytest.raises(AssertionError) as err_info:
            subtract(x=x)
        assert msg in str(err_info.value)
```

where

- we group multiple test functions for `subtract` under the same class `TestSubtract`.
- we [parameterize test](https://docs.pytest.org/en/stable/example/parametrize.html) to
  test different inputs.
- we catch errors using `pytest.raises` and check error messages.

For further usage like [fixture](https://docs.pytest.org/en/stable/fixture.html) and
other functionalities, please check
[pytest documentation](https://docs.pytest.org/en/stable/index.html) or
[existing tests](https://github.com/DeepRegNet/DeepReg/tree/main/test/unit) in DeepReg.
You can also [raise an issue](https://github.com/DeepRegNet/DeepReg/issues/new/choose)
for any questions.
# Quick Start

This is a set of simple tests to use DeepReg command line tools. More details and other
options can be found in [Command Line Tools](../docs/cli.html).

First, [install DeepReg](install.html) and change current directory to the root
directory of DeepReg.

## Train a registration network

Train a registration network using unpaired and labeled example data with a predefined
configuration:

```bash
deepreg_train --gpu "" --config_path config/unpaired_labeled_ddf.yaml --exp_name test
```

where:

- `--gpu ""` indicates using CPU. Change to `--gpu "0"` to use the GPU at index 0.
- `--config_path <filepath>` specifies the configuration file path.
- `--log_dir test` specifies the output folder. In this case, the output is saved in
  `logs/test`.

## Evaluate a trained network

Once trained, evaluate the network using a test dataset:

```bash
deepreg_predict --gpu "" --ckpt_path logs/test/save/ckpt-2 --split test
```

where:

- `--ckpt_path <filepath>` specifies the checkpoint file path.
- `--split test` specifies prediction on the test dataset.

## Warp an image

DeepReg provides a command line interface (CLI) tool to warp an image/label with a dense
displacement field (DDF):

```bash
deepreg_warp --image data/test/nifti/unit_test/moving_image.nii.gz --ddf data/test/nifti/unit_test/ddf.nii.gz --out logs/test_warp/out.nii.gz
```

where:

- `--image <filepath>` specifies the image/label file path.
- `--ddf <filepath>` specifies the ddf file path.
- `--out <filepath>` specifies the output file path.
# Visualisation tool

DeepReg provides a visuaisation tool which allows the user to generate various
visualisations from nifti images. The tool is compatible with outputs from
deepreg_predict as well as with other nifti images.

The visualisation tool currently offers four functionalities using the command
`deepreg_vis`, with their Python functions explained in the final section of this
document. Run the examples below after changing directory to the root folder `DeepReg`.

The creation of `.gif` files requires a movie writer that is compatible with
`matplotlib` like `ffmpeg` in order to be able to write `.gif` files. Please refer to
the [matplotlib documentation](https://matplotlib.org/3.3.1/api/animation_api.html) for
more details about writers. The `ffmpeg` writer was used to test functionality of this
visualisation tool, however, other writers can also be installed and used.

## General arguments for `deepreg_vis`

- `-m` or `--mode`: This specifies which mode to use to generate the visualisation. See
  below for available modes.
- `-i` or `--image-paths`: This is the path of the image or images that need to be used
  to generate the visualisation. Multiple paths can be passed using a comma separated
  string.
- `-s` or `--save-path`: This is the path to the directory where the visualisation will
  be saved. The name of the visualisation is auto generated, however, if `--fname`
  argument is available for a specific mode, then a custom filename can be specified.
  (default results in saving to current directory)

## GIF over image slices

This creates an animation which is an iteration through the image slices of a 3D image,
can be accessed by passing `--mode 0` or `-m 0`.

In addition to the general arguments, the additional arguments applicable to this mode
are:

- `--interval`: This argument is optional and can be used to specify the time, in
  milliseconds, between successive frames of an animation. (default=50)

The output will be a file with the same name as the original with a `.gif` extension. If
multiple image paths are passed in the `-i` or `--image-paths` argument then multiple
`.gif` files will be generated, one for each unique image path passed.

A simple example, which takes an image path and saves a `.gif` animation, is shown
below:

```bash
deepreg_vis -m 0 -i ./data/test/nifti/unit_test/moving_image.nii.gz -s logs
```

## GIF that shows warping

This functionality produces an animation showing warping for a single image slice using
a ddf, can be accessed by passing `--mode 1` or `-m 1`.

In addition to the general arguments, the additional arguments applicable to this mode
are:

- `--ddf-path`: This argument is required for this mode and specifies the path of the
  ddf to use for warping the image.
- `--slice-inds`: This argument is optional and can be used to specify the indexes to be
  used to generate the visualisation. Multiple indexes can be passed by using a comma
  separated string. (default results in a random slice being used)
- `--interval`: This argument is optional and can be used to specify the time, in
  milliseconds, between successive frames of an animation. (default=50)
- `--num-interval`: The number of intervals to use for warping. (default=100)

The output will be a file with the slice number appended to the original file name, with
a `.gif` extension. For example if a file is named `moving_image.nii.gz` and slice
number 2 and 3 three are chosen, the files produced will be `moving_image_slice_2.gif`
and `moving_image_slice_3.gif`.If multiple image paths are passed in the `-i` or
`--image-paths` argument then multiple `.gif` files will be generated, if multiple slice
indexes are specified then multiple `.gif` files are generated for each unique image
path passed.

A simple example, which takes an image, slice indexes and a ddf path and saves a `.gif`
animation for each slice, is shown below:

```bash
deepreg_vis -m 1 -i ./data/test/nifti/unit_test/moving_image.nii.gz --ddf-path "./data/test/nifti/unit_test/ddf.nii.gz" --slice-inds '2,3' -s logs
```

## Plot of image slices

This functionality produces a plot of image slices from a single or multiple images.
Each column is a different image and each row is a different slice, can be accessed by
passing `--mode 2` or `-m 2`.

In addition to the general arguments, the additional arguments applicable to this mode
are:

- `--slice-inds`: This argument is optional and can be used to specify the indexes to be
  used to generate the visualisation. Multiple can be passed by using a comma separated
  string. (default results in a random slice being used)
- `--col-titles`: This is optional. The title of the column to be used in order from
  left to right column. (default results in using file name specified in image path as
  column name)
- `--fname`: Optional argument of file name to save visualisation to; should end with an
  appropriate file extension like `.png` or `.jpeg`. (default='visualisation.png')

The output is a single file which contains a static visualisation. The visualisation is
different images in the columns and different slices in the rows.

A simple example, which takes three images and three slice indexes and saves a `.png`
file, is shown below (this will create a plot with 3 columns and 3 rows):

```
deepreg_vis -m 2 -i './data/test/nifti/unit_test/moving_image.nii.gz, ./data/test/nifti/unit_test/moving_image.nii.gz, ./data/test/nifti/unit_test/moving_image.nii.gz' --slice-inds '2,3,4' -s logs
```

## Tiled GIF over image slices

This functionality produces an animation with multiple animated images that are tiled
together, can be accessed by passing `--mode 3` or `-m 3`. The images used as input must
all be of the same size.

- `--interval`: This argument is optional and can be used to specify the time, in
  milliseconds, between successive frames of an animation. (default=50)
- `--size`: This is an optional argument and can be used to specify the number of rows
  and columns for the final tiled animation. For example `'2, 3'` means 2 rows and 3
  columns, in this case 6 images must be passed for the visualisation. (default='2,2')
- `--fname`: Optional argument of file name to save visualisation to; should end with
  `.gif`. (default='visualisation.gif')

The output is a single file which contains the tiled animation. The image paths must be
passed in the order in which they are to be tiled (order is from left to right and then
next row).

A simple example, which takes four images and creates an animation over the slices of
the four images in a tiled manner, is shown below:

```bash
deepreg_vis -m 3 -i './data/test/nifti/unit_test/moving_image.nii.gz, ./data/test/nifti/unit_test/moving_image.nii.gz, ./data/test/nifti/unit_test/moving_image.nii.gz, ./data/test/nifti/unit_test/moving_image.nii.gz' --size '2,2' -s logs
```

## Using the same functionality by importing python functions

The same functionality can be accessed via python functions.

The names of the functions are:

- `gif_slices`: equivalent to `--mode 0` or `-m 0`
- `gif_warp`: equivalent to `--mode 1` or `-m 1`
- `tile_slices`: equivalent to `--mode 2` or `-m 2`
- `gif_tile_slices`: equivalent to `--mode 3` or `-m 3`

The functions can be imported in a python script as follows:

```
from deepreg.vis import function_name
```

For usage details please see the function docstrings. Function docstrings can be
accessed in a python prompt by:

```
help(func_name)
```
# Registry

DeepReg adopts the
[registry system](https://github.com/DeepRegNet/DeepReg/blob/main/deepreg/registry.py)
to facilitate the addition of custom functionalities.

## Description

The class `Registry` maintains a dictionary mapping `(category, key)` to `value`, where

- `category` is the class category, e.g. `"backbone_class""` for backbone classes.
- `key` is the name of the registered class, e.g. `"unet"` for the class `UNet`.
- `value` is the registered class, e.g. `UNet` corresponding to `"unet"`.

A global variable `REGISTRY = Registry()` is defined to provide a central control of all
classes. The supported categories and the registered classes are resumed in the
[registered classes](registered_classes.html) page.

### Register a class

To register a class into `REGISTRY`, it is recommended to use `register` as a
**decorator**. Consider the `UNet` class for backbone as an example.

```python
from deepreg.registry import REGISTRY


@REGISTRY.register(category="backbone_class", name="unet")
class UNet:
    """UNet-style backbone."""
```

The decorator automatically registers the class upon import. To ensure that this class
is registered when `import deepreg` is called, this class and related parent modules
need to be imported in the `__init__.py` files.

For the purpose of code simplicity, a specific register function is defined for each
category. For instance, we can use `register_backbone` for `UNet`:

```python
from deepreg.registry import REGISTRY


@REGISTRY.register_backbone(name="unet")
class UNet:
    """UNet-style backbone."""
```

### Instantiate a class

To instantiate a registered example in `REGISTRY`, we call `build_from_config` , which
allows creating a class instance from a config directly. The config should be a
dictionary containing the key `name` and other required keys for the class. Please check
the [configuration documentation](configuration.html) and the docstring of the related
classes for detailed configuration requirements. The path of the registered classes are
resumed in the [registered classes](registered_classes.html) page.

For instance, to instantiate a UNet class,

```python
from deepreg.registry import REGISTRY

config = dict(name="unet",
              image_size=(16,16,16),
              out_channels=3,
              num_channel_initial=2,
              out_kernel_initializer="he_normal",
              out_activation="")
unet = REGISTRY.build_from_config(category="backbone_class", config=config)
```

where `kwargs` represents other required arguments.

Similarly, for the purpose of code simplicity, a specific build function is defined for
each category. For instance, we can use `build_backbone` for `UNet`:

```python
from deepreg.registry import REGISTRY


config = dict(name="unet",
              image_size=(16,16,16),
              out_channels=3,
              num_channel_initial=2,
              out_kernel_initializer="he_normal",
              out_activation="")
unet = REGISTRY.build_backbone(config=config)
```

## Example usages

Apart from the table of [registered classes](registered_classes.html), to further
explain how to use `REGISTRY` for using customized classes, detailed examples are
provided for the following categories:

- backbone
- loss

### Custom backbone

To register a custom backbone class, the steps are as follows

1. Subclass the `Backbone` and implement a custom backbone class.
2. Import `REGISTRY` and use the decorator `@REGISTRY.register_backbone` to register the
   custom class.
3. Use the registered name in the config for using the registered custom backbone.

Please check the self-contained
[example script](https://github.com/DeepRegNet/DeepReg/blob/main/examples/custom_backbone.py)
for further details.

### Custom Loss

To register a custom loss class for images and labels, the steps are as follows

1. Subclass the `tf.keras.losses.Loss` and implement a custom backbone class.
2. Import `REGISTRY` and use the decorator `@REGISTRY.register_loss` to register the
   custom class.
3. Use the registered name in the config for using the registered custom loss.

Please check the self-contained
[example script](https://github.com/DeepRegNet/DeepReg/blob/main/examples/custom_image_label_loss.py)
for further details. There is also a more complicated
[example of parameterized custom loss](https://github.com/DeepRegNet/DeepReg/blob/main/examples/custom_parameterized_image_label_loss.py).
# Logging

When running DeepReg, there are two types of printed messages:

1. DeepReg messages, which are entirely controlled by DeepReg code.
2. TensorFlow messages, which are automatically printed by TensorFlow during certain
   workflows.

Therefore, we use two independent environment variables to control these two types of
messages. The log level controls which types of log messages would be printed. Further
details are explained below.

## DeepReg logging

DeepReg uses the
[Python module `logging`](https://docs.python.org/3/library/logging.html) to log the
messages. The log level is controlled by the environment variable `DEEPREG_LOG_LEVEL`.
The levels are given in the table below. The default level is "2".

To adjust the logging level, there are two options. Take the training as example,

- You can first define the environment variable, then run the job.

  ```bash
  export DEEPREG_LOG_LEVEL=1
  deepreg_train --gpu "" --config_path config/unpaired_labeled_ddf.yaml --exp_name test
  ```

- Alternatively, you can define the environment variable while running the job.

  ```bash
  DEEPREG_LOG_LEVEL=1 deepreg_train --gpu "" --config_path config/unpaired_labeled_ddf.yaml --exp_name test
  ```

| DEEPREG_LOG_LEVEL | Behavior                                                                                   |
| ----------------- | ------------------------------------------------------------------------------------------ |
| "0"               | Log all messages, equivalent to `logging.DEBUG`. Same as log level "1".                    |
| "1"               | Log all messages, equivalent to `logging.DEBUG`.                                           |
| "2"               | Log all messages except DEBUG, equivalent to `logging.INFO`. (default)                     |
| "3"               | Log all messages except DEBUG and INFO, equivalent to `logging.WARNING`.                   |
| "4"               | Log all messages except DEBUG, INFO, and WARNING, equivalent to `logging.ERROR`.           |
| "5"               | Log all messages except DEBUG, INFO, WARNING, and ERROR, equivalent to `logging.CRITICAL`. |

## TensorFlow logging

With TensorFlow 2.3, its log level is controlled by the environment variable
`TF_CPP_MIN_LOG_LEVEL`. The levels are given in the table below. The default level is
"2".

To adjust the logging level, there are two options. Take the training as example,

- You can first define the environment variable, then run the job.

  ```bash
  export TF_CPP_MIN_LOG_LEVEL=1
  deepreg_train --gpu "" --config_path config/unpaired_labeled_ddf.yaml --exp_name test
  ```

- Alternatively, you can define the environment variable while running the job.

  ```bash
  TF_CPP_MIN_LOG_LEVEL=1 deepreg_train --gpu "" --config_path config/unpaired_labeled_ddf.yaml --exp_name test
  ```

| TF_CPP_MIN_LOG_LEVEL | Behavior                                            |
| -------------------- | --------------------------------------------------- |
| "0"                  | Log all messages.                                   |
| "1"                  | Log all messages except INFO.                       |
| "2"                  | Log all messages except INFO and WARNING. (default) |
| "3"                  | Log all messages except INFO, WARNING, and ERROR.   |
# Dataset Loader

## Dataset type

DeepReg provides six dataset loaders to support the following three different types of
datasets:

- **Paired images**

  Images are organized into moving and fixed image pairs.

  An example case is two-modalities intra-subject registration, such as registering one
  subject's MR image to the corresponding ultrasound image.

- **Unpaired images**

  Images may be considered independent samples.

  An example case is single-modality inter-subject registration, such as registering one
  CT image to another from different subjects.

- **Grouped images**

  Images are organized into multiple groups.

  An example case is single-modality intra-subject registration, such as registering
  time-series images within individual subjects, a group is one subject in this case.

For all three above cases, the images can be either unlabeled or labeled. A label is
represented by a boolean mask on the image, such as a segmentation of an anatomical
structure or landmark.

## Dataset requirements

To use the provided dataset loaders, other detailed images and labels requirements are
described in individual dataset loader sections. General requirements are described as
follows.

- Image

  - DeepReg currently supports 3D images. But images do not have to be of the same
    shape, and it will be resized to the required shape using linear interpolation.

  - Currently, DeepReg only supports images stored in Nifti files or H5 files. Check
    [Nifti_loader](https://github.com/DeepRegNet/DeepReg/blob/main/deepreg/dataset/loader/nifti_loader.py)
    and
    [h5_loader](https://github.com/DeepRegNet/DeepReg/blob/main/deepreg/dataset/loader/h5_loader.py)
    for more details.

  - **Images are automatically normalized** at per-image level: the intensity values x
    equals to `(x-min(x)+EPS) / (max(x)-min(x)+EPS)` so that its values are between
    [0,1]. Check `GeneratorDataLoader.data_generator` in
    [loader interface](https://github.com/DeepRegNet/DeepReg/blob/main/deepreg/dataset/loader/interface.py)
    for more details.

- Label

  - If an image is labeled, the label shape is recommended to be the same as the image
    shape. Otherwise, the resize might give unexpected behaviours. But each image can
    have more than one labels.<br>

    For instance, an image of shape `(dim1, dim2, dim3)`, its label shape can be
    `(dim1, dim2, dim3)` (single label) or `(dim1, dim2, dim3, num_labels)` (multiple
    labels).

  - **All labels are assumed to have values between [0, 1].** So DeepReg accepts binary
    segmentation masks or soft labels with float values between [0,1]. This is to
    prevent accidental use of non-one-hot encoding to represent multiple class labels.
    In case of multi labels, please use one-hot encoding to transform them into multiple
    channels such that each class has its own binary label.

  - When the images are paired, the moving and fixed images must have the same number of
    labels.

  - When there are multiple labels, it is assumed that the labels are ordered, such that
    the channel of index `label_idx` is the same anatomical or pathological structure.

  - Currently, if the data are labeled, each data sample must have at least one label.
    For missing labels or partially labelled data, consider using all-zero masks as a
    workaround.

  - See further discussion in [Label sampling](exp_label_sampling.html).

.. \_paired-images:

## Paired images

For paired images, each pair contains a moving image and a fixed image. Optionally,
corresponding moving label(s) and fixed label(s).

Specifically, given a pair of images

- When the image is unlabeled,
  - moving image of shape `(m_dim1, m_dim2, m_dim3)`
  - fixed image of shape `(f_dim1, f_dim2, f_dim3)`
- When the image is labeled and there is only one label,
  - moving image of shape `(m_dim1, m_dim2, m_dim3)`
  - fixed image of shape `(f_dim1, f_dim2, f_dim3)`
  - moving label of shape `(m_dim1, m_dim2, m_dim3)`
  - fixed label of shape `(f_dim1, f_dim2, f_dim3)`
- When the image is labeled and there are multiple labels,
  - moving image of shape `(m_dim1, m_dim2, m_dim3)`
  - fixed image of shape `(f_dim1, f_dim2, f_dim3)`
  - moving label of shape `(m_dim1, m_dim2, m_dim3, num_labels)`
  - fixed label of shape `(f_dim1, f_dim2, f_dim3, num_labels)`

### Sampling

For paired images, one epoch of the dataset iterates all the image pairs sequentially
with random orders. So each image pair is sampled once in each epoch with equal chance.
For validation or testing, the random seed is fixed to ensure consistency.

When an image has multiple labels, e.g. the segmentation of different organs in a CT
image, only one label will be sampled during training. In particular, only corresponding
labels will be sampled between a pair of moving and fixed images. In case of validation
or testing, instead of sampling one label per image, all labels will be iterated.

### Configuration

An example configuration for paired dataset is provided as follows.

```yaml
dataset:
  train:
    dir: "data/test/h5/paired/train" # folder containing data
    format: "h5" # nifti or h5
    labeled: true # true or false
  valid:
    dir: "data/test/h5/unpaired/test"
    format: "h5"
    labeled: true
  test:
    dir: "data/test/h5/unpaired/test"
    format: "h5"
    labeled: true
  type: "paired" # value should be paired / unpaired / grouped
  moving_image_shape: [16, 16, 16] # value should be like [dim1, dim2, dim3]
  fixed_image_shape: [8, 8, 8] # value should be like [dim1, dim2, dim3]
```

where, the configuration can be split into common configurations that shared by all
dataset types and specific configurations for paired images:

- Common configurations
  - `dir/train` gives the directory containing training data. Same for `dir/valid` and
    `dir/test`.
  - `format` can only be Nifti or h5 currently.
  - `type` can be paired, unpaired or grouped, corresponding to the dataset type
    described above.
  - `labeled` is a boolean indicating if the data is labeled or not.
- Paired images configurations
  - `moving_image_shape` is the shape of moving images, a list of three integers.
  - `fixed_image_shape` is the shape of fixed images, a list of three integers.

Optionally, multiple dataset directories can be specified, such that the data will be
sampled from several directories, for instance:

```yaml
dataset:
  train:
    dir: # folders containing data
      - "data/test/h5/paired/train1"
      - "data/test/h5/paired/train2"
    format: "h5" # nifti or h5
    labeled: true # true or false
  valid:
    dir: "data/test/h5/unpaired/test"
    format: "h5"
    labeled: true
  test:
    dir: "data/test/h5/unpaired/test"
    format: "h5"
    labeled: true
  type: "paired" # value should be paired / unpaired / grouped
  moving_image_shape: [16, 16, 16] # value should be like [dim1, dim2, dim3]
  fixed_image_shape: [8, 8, 8] # value should be like [dim1, dim2, dim3]
```

This is particularly useful when performing an
[experiment such as cross-validation](../tutorial/cross_val.html).

### File loader

For paired data, the specific requirements for data stored in Nifti and h5 files are
described as follows.

#### Nifti

Nifti data are stored in files with suffix `.nii.gz`. Each file should contain only one
3D or 4D tensor, corresponding to an image or a label.

`obs` is short for one observation of a data sample - a 3D image volume or a 3D/4D label
volume - and the name can be any string.

All image data should be placed under `moving_images/`, `fixed_images/` with respect to
the provided directory. The label data should be placed under `moving_labels/`, and
`fixed_labels/`, if available. These are _top_ directories.

File names should be consistent between top directories, e.g.:

- moving_images/
  - obs1.nii.gz
  - obs2.nii.gz
  - ...
- fixed_images/
  - obs1.nii.gz
  - obs2.nii.gz
  - ...
- moving_labels/
  - obs1.nii.gz
  - obs2.nii.gz
  - ...
- fixed_labels/
  - obs1.nii.gz
  - obs2.nii.gz
  - ...

Check
[test paired Nifti data](https://github.com/DeepRegNet/DeepReg/tree/main/data/test/nifti/paired)
as an example.

Optionally, the data may not be all saved directly under the top directory. They can be
further grouped in subdirectories as long as the data paths are consistent.

#### H5

H5 data are stored in files with suffix `.h5`. Hierarchical multi-level indexing is not
used. Each file should contain multiple key-value pairs and values are 3D or 4D tensors.
Each file is equivalent to a top folder in Nifti cases.

All image data should be stored in `moving_images.h5`, `fixed_images.h5`. The label data
should be stored in `moving_labels.h5`, and `fixed_labels.h5`, if available.

The keys should be consistent between files, e.g.:

- moving_images.h5 has keys:
  - "obs1"
  - "obs2"
  - ...
- fixed_images.h5 has keys:
  - "obs1"
  - "obs2"
  - ...
- moving_labels.h5 has keys:
  - "obs1"
  - "obs2"
  - ...
- fixed_labels.h5 has keys:
  - "obs1"
  - "obs2"
  - ...

Check
[test paired H5 data](https://github.com/DeepRegNet/DeepReg/tree/main/data/test/h5/paired)
as an example.

## Unpaired images

For unpaired images, all images are considered as independent and they must have the
same shape. Optionally, there are corresponding labels for the images.

Specifically,

- When the image is unlabeled,
  - image of shape `(dim1, dim2, dim3)`
- When the image is labeled and there is only one label,
  - image of shape `(dim1, dim2, dim3)`
  - label of shape `(dim1, dim2, dim3)`
- When the image is labeled and there are multiple labels,
  - image of shape `(dim1, dim2, dim3)`
  - label of shape `(dim1, dim2, dim3, num_labels)`

### Sampling

During each epoch, image pairs will be sampled without replacement. Therefore, given N
images, one epoch will thereby have floor(N / 2) image pairs. For validation or testing,
the random seed is fixed to ensure consistency.

In case of multiple labels, the sampling method is the same as in
[paired data](#sampling). In particular, the only corresponding label pairs will be
sampled between the two sampled images.

### Configuration

An example configuration for unpaired dataset is provided as follows.

```yaml
dataset:
  train:
    dir: "data/test/h5/paired/train" # folder containing data
    format: "h5" # nifti or h5
    labeled: true # true or false
  valid:
    dir: "data/test/h5/unpaired/test"
    format: "h5"
    labeled: true
  test:
    dir: "data/test/h5/unpaired/test"
    format: "h5"
    labeled: true
  type: "unpaired" # value should be paired / unpaired / grouped
  image_shape: [16, 16, 16] # value should be like [dim1, dim2, dim3]
```

where

- Common configurations

  Same as [paired images](#configuration).

- Unpaired images configurations
  - `image_shape` is the shape of images, a list of three integers.

### File loader

For unpaired data, the specific requirements for data stored in nifti and h5 files are
described as follows.

#### Nifti

Nifti data are stored in files with suffix `.nii.gz` or `.nii`. Each file must contain
only one 3D or 4D tensor, corresponding to an image or a label.

`obs` is short for one observation of a data sample - a 3D image volume or a 3D/4D label
volume - and the name can be any string.

All image data should be placed under `images/`. The label data should be placed under
`labels/`, if available. These are _top_ directories.

File names should be consistent between top directories, e.g.:

- images/
  - obs1.nii.gz
  - obs2.nii.gz
  - ...
- labels/
  - obs1.nii.gz
  - obs2.nii.gz
  - ...

Check
[test unpaired Nifti data](https://github.com/DeepRegNet/DeepReg/tree/main/data/test/nifti/unpaired)
as an example.

#### H5

F5 data are stored in files with suffix `.h5`. Hierarchical multi-level indexing is not
used. Each file should contain multiple key-value pairs and values are 3D or 4D tensors.
Each file is equivalent to a top folder in Nifti cases.

All image data should be placed under `images.h5`. The label data should be placed under
`labels.h5`, if available.

The keys should be consistent between files, e.g.:

- images.h5 has keys:
  - "obs1"
  - "obs2"
  - ...
- labels.h5 has keys:
  - "obs1"
  - "obs2"
  - ...

Check
[test unpaired H5 data](https://github.com/DeepRegNet/DeepReg/tree/main/data/test/h5/unpaired)
as an example.

## Grouped images

For grouped images, images may not be paired but organized into multiple groups. Each
group must have at least two images.

The requirements are the same as unpaired images. Specifically,

- When the image is unlabeled,
  - image of shape `(dim1, dim2, dim3)`
- When the image is labeled and there is only one label,
  - image of shape `(dim1, dim2, dim3)`
  - label of shape `(dim1, dim2, dim3)`
- When the image is labeled and there are multiple labels,
  - image of shape `(dim1, dim2, dim3)`
  - label of shape `(dim1, dim2, dim3, num_labels)`

### Sampling

For sampling image pairs, DeepReg provides the following options:

- **inter-group sampling**, where the moving image and fixed image come from different
  groups.
- **intra-group sampling**, where the moving image and fixed image come from the same
  group.
- **mixed sampling**, where the image pairs are mixed from inter-group sampling and
  intra-group sampling.

For validation or testing, the random seed is fixed to ensure consistency.

In case of multiple labels, the sampling method is the same as [paired data](#sampling).
In particular, only the corresponding label pairs will be sampled between the two
sampled images.

#### Intra-group

To form image pairs, the group and image are sampled sequentially at two stages,

1. Sample a group from which the moving and fixed images will be sampled.
2. Sample two different images from the group as moving and fixed images.<br> When
   sampling images from the same group, there are multiple options, denoted by
   `intra_group_option`:
   - `forward`: the moving image always has a smaller image index than fixed image.
   - `backward`: the moving image always has a larger image index than fixed image.
   - `unconstrained`: no constraint on the image index as long as the two images are
     different.

Therefore, each epoch generates the same number of image pairs as the number of groups,
where all groups will be first shuffled and iterated. The `intra_group_option` is useful
in implementing temporal-order sensitive sampling strategy.

#### Inter-group

To form image pairs, the group and image are sampled sequentially at two stages,

1. Sample the first group, from which the moving image will be sampled.
2. Sample the second group, from which the fixed image will be sampled.
3. Sample an image from the first group as moving image.
4. Sample an image from the second group as fixed image.

Therefore, each epoch generates the same number of image pairs as the number of groups,
where all groups will be first shuffled and iterated.

#### Mixed

Optionally, it is possible to mix inter-group and intra-group sampling by specifying the
intra-group image sampling probability `intra_group_prob`=[0,1]. The value 0 means
entirely inter-group sampling and 1 means entirely intra-group sampling.

Given 0<p<1, when generating intra-group pairs, there is (1-p)\*100% chance to sample
the fixed images from a different group, after sampling the moving image from the
current intra-group images.

#### Iterated

Optionally, it is possible to generate all combinations of inter-/intra-group image
pairs, with `sample_image_in_group` set to false. This is originally designed for
evaluation. Mixing inter-/intra-group sampling is not supported with with
`sample_image_in_group` set to false.

### Configuration

An example configuration for grouped dataset is provided as follows.

```yaml
dataset:
  train:
    dir: "data/test/h5/paired/train" # folder containing data
    format: "h5" # nifti or h5
    labeled: true # true or false
  valid:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
  test:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
  type: "unpaired" # value should be paired / unpaired / grouped
  intra_group_prob: 1 # probability of intra-group sampling, value should be between 0 and 1
  intra_group_option: "forward" # option for intra-group sampling, value should be forward / backward / unconstrained
  sample_image_in_group: true # true if sampling one image per group, value should be true / false
  image_shape: [16, 16, 16] # value should be like [dim1, dim2, dim3]`
```

where

- Common configurations

  Same as [paired images](#configuration).

- Grouped images configurations
  - `intra_group_prob`, a value between 0 and 1, 0 is for inter-group only and 1 is for
    intra-group only.
  - `intra_group_option`, forward or backward or unconstrained, as described above.
  - `sample_image_in_group`, true if sampling one image at a time per group, false if
    generating all possible pairs.

### File loader

For grouped data, the specific requirements for data stored in Nifti and h5 files are
described as follows.

#### Nifti

Nifti data are stored in files with suffix `.nii.gz`. Each file should contain only one
3D or 4D tensor, corresponding to an image or a label.

`obs` is short for one observation of a data sample - a 3D image volume or a 3D/4D label
volume - and the name can be any string.

All image data should be placed under `images/`. The label data should be placed under
`labels/`, if available. These are _top_ directories.

The leaf directories will be considered as different groups, and file names should be
consistent between top directories, e.g.:

- images
  - group1
    - obs1.nii.gz
    - obs2.nii.gz
    - ...
  - ...
- labels
  - group1
    - obs1.nii.gz
    - obs2.nii.gz
    - ...
  - ...

Check
[test grouped Nifti data](https://github.com/DeepRegNet/DeepReg/tree/main/data/test/nifti/grouped)
as an example.

#### H5

H5 data are stored in files with suffix `.h5`. Hierarchical multi-level indexing is not
used. Each file should contain multiple key-value pairs and values are 3D or 4D tensors.
Each file is equivalent to a top folder in Nifti cases.

All image data should be placed under `images.h5`. The label data should be placed under
`labels.h5`, if available.

The keys must satisfy a specific format, `group-%d-%d`, where `%d` represents an integer
number. The first number corresponds to the group index, and the second number
corresponds to the observation index. For example, `group-3-2` corresponds to the second
observation from the third group.

The keys should be consistent between files, e.g.:

- images.h5 has keys:
  - "group-1-1"
  - "group-1-2"
  - ...
  - "group-2-1"
  - ...
- labels.h5 has keys:
  - "group-1-1"
  - "group-1-2"
  - ...
  - "group-2-1"
  - ...

Check
[test grouped H5 data](https://github.com/DeepRegNet/DeepReg/tree/main/data/test/h5/grouped)
as an example.
# Label sampling

Images may have multiple labels, such as with segmentation of different organs in CT
scans. In this case, for each sampled image pair, one label pair is randomly chosen by
default.

## Corresponding label pairs

When using multiple labels, ensure the labels are ordered correctly. `label_idx` in
`[width, height, depth, label_idx]` must be the same anatomical or pathological
structure; a corresponding label pair between the moving and fixed labels.

## Consistent label pairs

Consistent label pairs between a pair of moving and fixed labels requires:

1. The two images have the same number of labels, and
2. The labels have the same order

When a pair of moving and fixed images have inconsistent label pairs, label
dissimilarity cannot be defined. The following applies:

- When using the unpaired-labeled-image loader, consistent label pairs are required;
- When using the grouped-labeled-image loader, consistent label pairs are required
  between intra-group image pairs;
- When mixing intra-inter-group images in the grouped-labeled-image loader, consistent
  label pairs are required between all intra-group and inter-group image pairs.

However,

- When using the paired-labeled-image loader, consistent label pairs are not required
  between different image pairs;
- When using the grouped-labeled-image loader without mixing intra-group and inter-group
  images, consistent label pairs are not required between different image groups.

## Partially labeled image data or missing labels

When one of the label dissimilarity measures prevents accidentally missing labels. When
appropriate, enable training with missing labels with placeholder all-zero masks for the
labels.

## Option for iterating all available label pairs

This option is default for testing. All the label pairs will be sampled once for each
sampled image pair. This option is not supported when mixing intra-group and inter-group
image pairs.
# Registered Classes

> This file is generated automatically.

The following tables contain all registered classes with their categories and keys.

## Backbone

The category is `backbone_class`. Registered keys and values are as following.

| key      | value                                         |
| :------- | :-------------------------------------------- |
| "global" | `deepreg.model.backbone.global_net.GlobalNet` |
| "local"  | `deepreg.model.backbone.local_net.LocalNet`   |
| "unet"   | `deepreg.model.backbone.u_net.UNet`           |

## Model

The category is `model_class`. Registered keys and values are as following.

| key           | value                                    |
| :------------ | :--------------------------------------- |
| "conditional" | `deepreg.model.network.ConditionalModel` |
| "ddf"         | `deepreg.model.network.DDFModel`         |
| "dvf"         | `deepreg.model.network.DVFModel`         |

## Loss

The category is `loss_class`. Registered keys and values are as following.

| key             | value                                                     |
| :-------------- | :-------------------------------------------------------- |
| "bending"       | `deepreg.loss.deform.BendingEnergy`                       |
| "cross-entropy" | `deepreg.loss.label.CrossEntropyLoss`                     |
| "dice"          | `deepreg.loss.label.DiceLoss`                             |
| "gmi"           | `deepreg.loss.image.GlobalMutualInformationLoss`          |
| "gncc"          | `deepreg.loss.image.GlobalNormalizedCrossCorrelationLoss` |
| "gradient"      | `deepreg.loss.deform.GradientNorm`                        |
| "jaccard"       | `deepreg.loss.label.JaccardLoss`                          |
| "lncc"          | `deepreg.loss.image.LocalNormalizedCrossCorrelationLoss`  |
| "ssd"           | `deepreg.loss.label.SumSquaredDifferenceLoss`             |

## Data Augmentation

The category is `da_class`. Registered keys and values are as following.

| key      | value                                                |
| :------- | :--------------------------------------------------- |
| "affine" | `deepreg.dataset.preprocess.RandomAffineTransform3D` |
| "ddf"    | `deepreg.dataset.preprocess.RandomDDFTransform3D`    |

## Data Loader

The category is `data_loader_class`. Registered keys and values are as following.

| key        | value                                                       |
| :--------- | :---------------------------------------------------------- |
| "grouped"  | `deepreg.dataset.loader.grouped_loader.GroupedDataLoader`   |
| "paired"   | `deepreg.dataset.loader.paired_loader.PairedDataLoader`     |
| "unpaired" | `deepreg.dataset.loader.unpaired_loader.UnpairedDataLoader` |

## File Loader

The category is `file_loader_class`. Registered keys and values are as following.

| key     | value                                                 |
| :------ | :---------------------------------------------------- |
| "h5"    | `deepreg.dataset.loader.h5_loader.H5FileLoader`       |
| "nifti" | `deepreg.dataset.loader.nifti_loader.NiftiFileLoader` |
# Command Line Tools

With DeepReg installed, multiple command line tools are available, currently including:

- `deepreg_train`, for training a registration network.
- `deepreg_predict`, for evaluating a trained network.
- `deepreg_warp`, for warping an image with a dense displacement field.

## Train

`deepreg_train` accepts the following arguments via command line tools. More
configuration can be specified in the configuration file. Please see
[configuration file](configuration.html) for further details.

### Required arguments

- **GPU**:

  `--gpu` or `-g`, specifies the index or indices of GPUs for training.

  Example usage:

  - `--gpu ""` for CPU only
  - `--gpu "0"` for using only GPU 0
  - `--gpu "0,1"` for using GPU 0 and 1.

- **Configuration**:

  `--config_path` or `-c`, specifies the configuration file for training.

  The path must end with `.yaml`.

  Optionally, multiple paths can be specified, and the configuration will be merged. In
  case of conflicts, values are overwritten by the last config file defining them.

  Example usage:

  - `--config_path config1.yaml` for using one single configuration file.
  - `--config_path config1.yaml config2.yaml` for using multiple configuration files.

### Optional arguments

- **CPU allocation**:

  `--num_workers`, if given, TensorFlow will use limited CPUs.

  By default, it uses only 1 CPUs. Setting it to non-positive values will be using all
  CPUs.

  Example usage:

  - `--num_workers 2` for using at most 2 CPUs.

- **GPU memory allocation**:

  `--gpu_allow_growth` or `-gr`, if given, TensorFlow will only grow the memory usage as
  is needed.

  By default, it allocates all available GPU memory.

  Example usage:

  - `--gpu_allow_growth`, no extra argument is needed.

- **Load checkpoint**:

  `--ckpt_path` or `-k`, specifies the path of the saved model checkpoint, so that the
  training will be resumed from the given checkpoint.

  The path must end with `.ckpt`.

  By default, it starts training from a random initialization.

  Example usage:

  - `--ckpt_path weights-epoch2.ckpt` for reloading the given checkpoint.

- **Log directory**:

  `--log_dir`, specifies the log directory for logging output information and results.

  By default, it is `logs` under the package root.

  Example usage:

  - `--log_dir logs` for specifying the log directory `logs/` under current directory.

- **Experiment name**:

  `--exp_name` or `-n`, specifies the name of an experiment (every time a training or a
  prediction is run), which will be used together with the log directory (via `log_dir`)
  to specify the sub-folder that saves the output information and results from
  individual experiments (runs).

  If this is not provided, it creates a timestamp-named sub-folder under the `log_dir`,
  by default, e.g. `logs/20200810-194042/`.

  Example usage:

  - `--exp_name test --log_dir logs` for saving under `logs/test/`.
  - `--log_dir logs` for saving under `logs/20210101-120000/`, assuming
    `20210101-120000` is current time.
  - `--exp_name test` for saving under `DeepReg/logs/test/`, assuming `DeepReg` is the
    package root.

- **Maximum number of epochs**:

  `--max_epochs`, specifies the maximum number of epochs for training and overwrites the
  value defined in the configuration.

  By default, the value is -1, meaning the number of epochs will be defined by
  configuration.

  Example usage:

  - `--max_epochs 2` for run training only for two epochs.

### Output

During the training, multiple output files will be saved in the log directory
`logs/log_dir`, where `log_dir` is specified in the arguments, otherwise a timestamped
folder name will be used. The output files are:

- `config.yaml` is a backup of the used configuration. It can be used for prediction. In
  case of multiple configuration files, a merged configuration file will be saved.
- `train/` and `validation/` are the directories that save tensorboard logs on metrics.
- `save/` is the directory containing saved checkpoints of the trained network.

## Predict

`deepreg_predict` accepts the following arguments via command line tools. More
configuration can be specified in the configuration file. Please see
[configuration file](configuration.html) for further details.

### Required arguments

- **GPU**:

  `--gpu` or `-g`, specifies the index or indices of GPUs for training.

  Example usage:

  - `--gpu ""` for CPU only
  - `--gpu "0"` for using only GPU 0
  - `--gpu "0,1"` for using GPU 0 and 1.

- **Model checkpoint**:

  `--ckpt_path` or `-k`, specifies the path of the saved model checkpoint, so that the
  trained model will be loaded for evaluation.

  The path must end with `.ckpt`.

  Example usage:

  - `--ckpt_path weights-epoch2.ckpt` for reloading the given checkpoint.

- **Evaluation data**:

  `--split`, specifies in which data set the prediction is performed.

  It must be one of `train` / `valid` / `test`.

  Example usage:

  - `--split test` for evaluating the model on test split.

### Optional arguments

- **CPU allocation**:

  `--num_workers`, if given, TensorFlow will use limited CPUs.

  By default, it uses all available CPUs.

  Example usage:

  - `--num_workers 2` for using at most 2 CPUs.

- **GPU memory allocation**:

  `--gpu_allow_growth` or `-gr`, if given, TensorFlow will only grow the memory usage as
  is needed.

  By default, it allocates all availables in the GPU memory.

  Example usage:

  - `--gpu_allow_growth`, no extra argument is needed.

- **Log directory**:

  `--log_dir`, specifies the log directory for logging output information and results.

  By default, it is `logs` under the package root.

  Example usage:

  - `--log_dir logs` for specifying the log directory `logs/` under current directory.

- **Experiment name**:

  `--exp_name` or `-n`, specifies the name of an experiment (every time a training or a
  prediction is run), which will be used together with the log directory (via `log_dir`)
  to specify the sub-folder that saves the output information and results from
  individual experiments (runs).

  If this is not provided, it creates a timestamp-named sub-folder under the `log_dir`,
  by default, e.g. `logs/20200810-194042/`.

  Example usage:

  - `--exp_name test --log_dir logs` for saving under `logs/test/`.
  - `--log_dir logs` for saving under `logs/20210101-120000/`, assuming
    `20210101-120000` is current time.
  - `--exp_name test` for saving under `DeepReg/logs/test/`, assuming `DeepReg` is the
    package root.

- **Batch size**:

  `--batch_size` or `-b`, specifies the number of samples per step for prediction. If
  using multiple GPUs, i.e. `n` GPUs, each GPU will have mini batch size
  `batch_size / n`. Thus, `batch_size` should be divided by `n` evenly.

  The default value is 1.

  Example usage:

  - `--batch_size 2` for using a global mini-batch size of 2.

- **Save outputs in Nifti format**:

  The predicted 3D tensors can be saved in Nifti format for further calculation.

  By default, it saves outputs in Nifti format.

  Example usage:

  - `--save_nifti`, for saving the outputs in Nifti format.
  - `--no_nifti`, for not saving the outputs in Nifti format.

- **Save outputs in png format**:

  The predicted 3D tensors can be saved as a slice of 2D images for quick visualization.

  As values have to be normalized between 0~255 (or 0~1) for png files (Nifti files are
  not impacted), all images (`moving_image`, `fixed_image` and `pred_fixed_image`) and
  displacement/velocity fields (`ddf` and `dvf`) will be normalized before being saved.
  Labels (`moving_label`, `fixed_label` and `pred_fixed_label`) are not affected as they
  are already within 0~1.

  By default, it saves the outputs in png format.

  Example usage:

  - `--save_png`, for saving the outputs in png format.
  - `--no_png`, for not saving the outputs in png format.

- **Configuration**:

  `--config_path` or `-c`, specifies the configuration file for prediction.

  The path must end with `.yaml`.

  By default, it uses the configuration file saved in the directory of the given
  checkpoint.

  Example usage:

  - `--config_path config1.yaml` for using one single configuration file.

### Output

During the evaluation, multiple output files will be saved in the log directory
`logs/log_dir/mode` where

- `log_dir` is defined in arguments, or a timestamped folder name will be used;
- `mode` is `train` or `valid` or `test`, specified by the argument.

The saved files include:

- Metrics to evaluate the registration performance
  - `metrics.csv` saves the metrics on all samples. Each line corresponds to a data
    sample.
  - `metrics_stats_per_label.csv` saves the mean, median and std of each metrics on all
    samples with the same label index.
  - `metrics_stats_overall.csv` saves a set of commonly used statistics (such as mean
    and std) on the metrics over all samples.
- Inputs and predictions for each pair of image.

  Each pair has its own directory and the followings tensors are saved inside if
  available. Tensors can be saved in Nifti format (one single file) or in png format
  (one folder contains all image slices, ordered by depth) or both.

  - `ddf`, `dvf`, `affine`

    DDF stands for dense displacement field; DVF stands for dense (static) velocity
    field.

    The 12 parameters of affine transformation are saved in `affine.txt`.

  - `moving_image`, `fixed_image` and `pred_fixed_image`

    `pred_fixed_image` is the warped moving image if the network predicts a DDF or a DVF
    or an affine transformation.

  - `moving_label`, `fixed_label` and `pred_fixed_label` under directory `label_i` if
    the sample is labeled and `i` is the label index.

    `pred_fixed_label` is the predicted label in the fixed image space. In many cases,
    this is equivalent to the warped moving label, if the network predicts a DDF or a
    DVF or an affine transformation.

## Warp

`deepreg_warp` accepts the following arguments:

### Required arguments

- **Image file**:

  `--image` or `-i`, specifies the file path of the image/label.

  The image/label should be saved in a Nifti file with suffix `.nii` or `.nii.gz`. The
  image/label should be a 3D / 4D tensor, where the first three dimensions correspond to
  the moving image shape and the fourth can be a channel of features.

  Example usage:

  - `--image input_image.nii.gz`

- **DDF file**:

  `--ddf` or `-d`, specifies the file path of the DDF.

  The DDF should be saved in a Nifti file with suffix `.nii` or `.nii.gz`. The DDF
  should be a 4D tensor, where the first three dimensions correspond to the fixed image
  shape and the fourth dimension has 3 channels corresponding to x, y, z axes.

  Example usage:

  - `--image input_DDF.nii.gz`

### Optional arguments

- **Output directory**:

  `--out` or `-o`, specifies the file path for the output.

  The path should end with `.nii` or `.nii.gz`, otherwise the output path will be
  corrected automatically based on the given path.

  By default, it saves the output as `warped.nii.gz` in the current directory.

  Example usage:

  - `--out output_image.nii.gz`

### Output

The warped image is saved in the given output file path, otherwise the default file path
`warped.nii.gz` will be used.

## Visualise

In addition to the images in the output, DeepReg provides a set of tools with the
command `deepreg_vis`. See more details in
[its usage documentation](visualisation_tool.html).
# Configuration File

In addition to the arguments provided to the command line tools, detailed training and
prediction configuration is specified in a `YAML` file. The configuration file contains
two sections, `dataset` and `train`. Within `dataset` one specifies the data file
formats, sizes, as well as the data loader to use. The `train` section specifies
parameters related to the neural network.

## Dataset section

The `dataset` section specifies the path to the data to be used during training, the
data loader to use as well as the specific arguments to configure the data loader.

### Split keys - Required

The data paths, data format, and label availability of the training, validation and
testing data are specified under the corresponding sections separately:

```yaml
dataset:
  train:
    dir: "data/test/h5/paired/train"
    format: "h5"
    labeled: true
  valid:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
  test:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
```

For data paths, multiple dataset directories can be specified, such that data are
sampled across several folders:

```yaml
dataset:
  train:
    dir:
      - "data/test/h5/paired/train1"
      - "data/test/h5/paired/train2"
    format: "h5"
    labeled: true
  valid:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
  test:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
```

For data format, different formats requires different file structure and format. Thus
different file loader will be used. Check the
[data loader configuration](dataset_loader.html) for more details.

Currently, DeepReg file loaders support Nifti and H5 file types - alternate file formats
will raise errors in the data loaders. To indicate which format to use, pass a string to
this field as either "nifti" or "h5":

```yaml
dataset:
  train:
    dir:
      - "data/test/nifti/paired/train1"
      - "data/test/nifti/paired/train2"
    format: "nifti"
    labeled: true
```

The `labeled` key indicates whether segmentation labels are available for training or
evaluation. Use `true` and `false` to indicate the availability and unavailability
correspondingly. In particular, if the value passed is false, the labels will not be
used even if they are available in the associated directories.

```yaml
dataset:
  train:
    dir:
      - "data/test/nifti/paired/train1"
      - "data/test/nifti/paired/train2"
    format: "nifti"
    labeled: false # labels are not available
```

### Type key - Required

The type of data loader used will depend on how one wants to train the network.
Currently, DeepReg data loaders support the `paired`, `unpaired`, and `grouped` training
strategies. Passing a string that doesn't match any of the above would raise an error.
The data loader type would be specified using the `type` key:

```yaml
dataset:
  train:
    dir: "data/test/h5/paired/train"
    format: "h5"
    labeled: true
  valid:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
  test:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
  type: "paired" # one of "paired", "unpaired" or "grouped"
```

### Data loader dependent keys

Depending on which string is passed to the `type` key, DeepReg will initialize a
different data loader instance with different sampling strategies. These are described
in depth in the [dataset loader configuration](dataset_loader.html) documentation. Here
we outline the arguments necessary to configure the different data loaders.

#### Sample_label - Required

In the case that we have more than one label per image, we need to inform the loader
which one to use. We can use the `sample_label` argument to indicate which method to use
during training.

- `all`: for one image that has x number of labels, the loader yields x image-label
  pairs with the same image. Occurs over all images, over one epoch.
- `sample`: for one image that has x number of labels, the loader yields 1 image-label
  pair randomly sampled from all the labels. Occurs for all images in one epoch.

During validation and testing (ie for `valid` and `test` directories), data loaders will
be built to sample `all` the data-label pairs, regardless of the argument passed to
`sample_label`.

```yaml
dataset:
  train:
    dir: "data/test/h5/paired/train"
    format: "h5"
    labeled: true
  valid:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
  test:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
  type: "paired" # one of "paired", "unpaired" or "grouped"
  sample_label: "sample" # one of "sample", "all" or None
```

In the case the `labeled` argument is false, the sample_label is unused, but still must
be passed. Additionally, if the tensors in the files only have one label, regardless of
the `sample_label` argument, the data loader will only pass the one label to the
network.

For more details please refer to
[Read The Docs](https://deepreg.readthedocs.io/en/latest/docs/exp_label_sampling.html).

#### Paired

- `moving_image_shape`: Union[Tuple[int, ...], List[int]] of ints, len 3, corresponding
  to (dim1, dim2, dim3) of the 3D moving image.
- `fixed_image_shape`: Union[Tuple[int, ...], List[int]] of ints, len 3, corresponding
  to (dim1, dim2, dim3) of the 3D fixed image.

```yaml
dataset:
  train:
    dir: "data/test/h5/paired/train"
    format: "h5"
    labeled: true
  valid:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
  test:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
  type: "paired" # one of "paired", "unpaired" or "grouped"
  sample_label: "sample" # one of "sample", "all" or None
  moving_image_shape: [16, 16, 3]
  fixed_image_shape: [16, 16, 3]
```

#### Unpaired

- `image_shape`: Union[Tuple[int, ...], List[int]] of ints, len 3, corresponding to
  (dim1, dim2, dim3) of the 3D image.

```yaml
dataset:
  train:
    dir: "data/test/h5/paired/train"
    format: "h5"
    labeled: true
  valid:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
  test:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
  type: "unpaired" # one of "paired", "unpaired" or "grouped"
  sample_label: "sample" # one of "sample", "all" or None
  image_shape: [16, 16, 3]
```

#### Grouped

- `intra_group_prob`: float, between 0 and 1. Passing 0 would only generate inter-group
  samples, and passing 1 would only generate intra-group samples.
- `sample_label`: method for sampling the labels "sample", "all".
- `intra_group_option`: str, "forward", "backward, or "unconstrained"
- `sample_image_in_group`: bool, if true, only one image pair will be yielded for each
  group, so one epoch has num_groups pairs of data, if false, iterate through this
  loader will generate all possible pairs.
- `image_shape`: Union[Tuple[int, ...], List[int]] len 3, corresponding to (dim1, dim2,
  dim3) of the 3D image.

```yaml
dataset:
  train:
    dir: "data/test/h5/paired/train"
    format: "h5"
    labeled: true
  valid:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
  test:
    dir: "data/test/h5/paired/test"
    format: "h5"
    labeled: true
  type: "grouped" # one of "paired", "unpaired" or "grouped"
  sample_label: "sample" # one of "sample", "all" or None
  image_shape: [16, 16, 3]
  sample_image_in_group: true
  intra_group_prob: 0.7
  intra_group_option: "forward"
```

See the [dataset loader configuration](dataset_loader.html) for more details.

## Train section

The `train` section defines the neural network training hyper-parameters, by specifying
subsections, `method`, `backbone`, `loss`, `optimizer`, `preprocess` and other training
hyper-parameters, including `epochs` and `save_period`.

### Method - required

The `method` argument defines the registration type. It must be a string. Feasible
values are: `ddf`, `dvf`, and `conditional`, corresponding to the dense displacement
field (DDF) based model, dense velocity field (DDF) based model, and conditional model
presented in the [registration tutorial](../tutorial/registration.html).

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
```

### Backbone - required

The `backbone` subsection is used to define the network, with all the network-specific
arguments under the same indent. The first argument should be the argument `name`, which
should be string type, one of "unet", "local" or "global", to define a UNet, LocalNet or
GlobalNet backbone, respectively. With Registry functionalities, you can also define
your own networks to pass to DeepReg train via config.

The `num_channel_initial` is used to define the number of initial channels for the
network, and should be int type.

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
  backbone:
    name: "unet" # One of unet, local, global: networks currently supported by DeepReg
    num_channel_initial: 16 # Int type, number of initial channels in the network. Controls the network size.
```

#### UNet

The UNet model requires several additional arguments to define its structure:

- `depth`: int, defines the depth of the UNet from first to bottom, bottleneck layer.
- `pooling`: Boolean, pooling method used for down-sampling. True: non-parametrized
  pooling will be used, False: conv3d will be used.
- `concat_skip`: Boolean, concatenation method for skip layers in UNet. True:
  concatenation of layers, False: addition is used instead.

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
  backbone:
    name: "unet" # One of unet, local, global
    num_channel_initial: 16 # Int type, number of initial channels in the network. Controls the network size.
    depth: 3
    pooling: false
    concat_skip: true
```

#### LocalNet

The LocalNet has an encoder-decoder structure and extracts information from tensors at
one or multiple resolution levels. We can define which levels to extract info from with
the `extract_levels` argument.

- `depth`: Depth of the encoder, `depth=2` means there are in total 3 layers where 0 is
  the top layer and 2 is the bottom.
- `extract_levels`: indices of layer from which the output will be extracted, the value
  range is `[0, depth]` both side inclusive.

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
  backbone:
    name: "local" # One of unet, local, global
    num_channel_initial: 16 # Int type, number of initial channels in the network. Controls the network size.
    depth: 2
    extract_levels: [0, 1, 2]
```

#### GlobalNet

The GlobalNet has a U-net like encoder to encode the image and uses the bottleneck layer
to output an affine transformation using a CNN.

- `depth`: Depth of the encoder, `depth=2` means there are in total 3 layers where 0 is
  the top layer and 2 is the bottom.

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
  backbone:
    name: "global" # One of unet, local, global
    num_channel_initial: 16 # Int type, number of initial channels in the network. Controls the network size.
    depth: 4
```

### Loss - required

This section defines the loss in training.

There are three different categories of losses in DeepReg:

- **image loss**: loss between the fixed image and predicted fixed image (warped moving
  image).
- **label loss**: loss between the fixed label and predicted fixed label (warped moving
  label).
- **regularization loss**: loss on predicted dense displacement field (DDF).

Not all losses are applicable for all models, the details are in the following table.

|                     | DDF / DVF                      | Conditional    |
| ------------------- | ------------------------------ | -------------- |
| Image Loss          | Applicable                     | Non-applicable |
| Label Loss          | Applicable if data are labeled | Applicable     |
| Regularization Loss | Applicable                     | Non-applicable |

The configuration for non-applicable losses will be ignored without errors. The loss
will also be ignored if the weight is zero. However, each model must define at least one
loss, otherwise error will be raised by TensorFlow.

For each loss, there are multiple existing loss functions to choose. The registry
mechanism can also be used to use custom loss functions. Please read the
[registry documentation](registry.html) for more details.

#### Image

The image loss calculates dissimilarity between warped image tensors and fixed image
tensors.

- `weight`: float type, the weight of individual loss element in the total loss
  function.
- `name`: string type, one of "lncc", "ssd" or "gmi".

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
  backbone:
    name: "local" # One of unet, local, global
    num_channel_initial: 16 # Int type, number of initial channels in the network. Controls the network size.
    extract_levels: [0, 1, 2]
  loss:
    image:
      name: "lncc" # other options include "lncc", "ssd" and "gmi", for local normalised cross correlation,
      weight: 0.1
```

The following are the DeepReg image losses. Additional arguments should be added at the
same indent level:

- `lncc`: Calls a local normalized cross-correlation type loss. Requires the following
  arguments:

  - `kernel_size`: int, optional, default=9. Kernel size or kernel sigma for
    kernel_type="gaussian".
  - `kernel_type`: str, optional, default="rectangular". One of "rectangular",
    "triangular" or "gaussian"

- `ssd`: Calls a sum of squared differences loss. No additional arguments required.

- `gmi`: Calls a global mutual informatin loss. Requires the following arguments:
  - `num_bins`: int, optional, default=23. Number of bins for intensity.
  - `sigma_ratio`: float, optional, default=0.5. A hyperparameter for the Gaussian
    kernel density estimation.

#### Label

The label loss calculates dissimilarity between labels.

All default DeepReg losses can be used as multi-scale or single scale losses.
Multi-scale losses require a kernel Additionally, all losses can be weighted, so the
following two arguments are global to all provided losses:

- `weight`: float type, weight of individual loss element in total loss function.
- `scales`: list of ints, or None. Optional argument. If you do not pass this argument
  (or pass the list [0], the value `null` or an empty value pair), the loss is
  calculated at a single scale. If you pass a list of length > 1, a multi-scale loss
  will be used. WARNING: an empty list ([]) will raise an error.
- `kernel`: str, "gaussian" or "cauchy", default "gaussian". Optional argument. Defines
  the kernel to use for multi-scale losses.

EG.

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
  backbone:
    name: "local" # One of unet, local, global
    num_channel_initial: 16 # Int type, number of initial channels in the network. Controls the network size.
    extract_levels: [0, 1, 2]
  loss:
    label:
      weight: 1.0
      name: "dice" # options include "dice", "cross-entropy", "mean-squared", "generalised_dice" and "jaccard"
      scales: [1, 2]
```

The default losses require the following arguments. Additional arguments should be added
at the same indent level:

- `dice`: Calls a Dice loss on the labels, requires the following arguments:

  - `binary`: bool, default is false. If true, the tensors are thresholded at 0.5.
  - `background_weight`: float, default=0.0. `background_weight` weights the foreground
    and background classes by replacing the labels of 1s and 0s with
    `(1-background_weight)` and `background_weight`, respectively.

- `cross-entropy`: Calls a cross-entropy loss between labels, requires the following
  arguments:

  - `binary`: bool, default is false. If true, the tensors are thresholded at 0.5.
  - `background_weight`: float, default=0.0. `background_weight` weights the foreground
    and background classes by replacing the labels of 1s and 0s with
    `(1-background_weight)` and `background_weight`, respectively.

- `jaccard`: - `binary`: bool, default is false. If true, the tensors are thresholded at
0.5.
<!-- - `background_weight`: float, default=0.0. `background_weight` weights the foreground and background classes by replacing the labels of 1s and 0s with (1-background_weight) and background_weight, respectively. -->

#### Regularization

The regularization section configures the losses for the DDF. To instantiate this part
of the loss, pass "regularization" into the config file as a field.

- `weight`: float type, the weight of the regularization loss.
- `name`: string type, the type of deformation energy to compute. Options include
  "bending", "gradient"

If the `gradient` loss is used, another argument must be passed at the same indent
level: - `l1`: bool. Indicates whether to calculate the L1-norm (true) or L2-norm
(false) gradient loss of the ddf.

EG.

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
  backbone:
    name: "local" # One of unet, local, global
    num_channel_initial: 16 # Int type, number of initial channels in the network. Controls the network size.
    extract_levels: [0, 1, 2]
  loss:
    regularization:
      weight: 0.5 # weight of regularization loss
      name: "bending" # options include "bending", "gradient"
```

or

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
  backbone:
    name: "local" # One of unet, local, global
    num_channel_initial: 16 # Int type, number of initial channels in the network. Controls the network size.
    extract_levels: [0, 1, 2]
  loss:
    regularization:
      weight: 0.5 # weight of regularization loss
      name: "gradient" # options include "bending", "gradient"
      l1: false
```

#### Composite Loss

The loss function can be a composite of different loss categories by adding all fields
in the same configuration file.

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
  backbone:
    name: "local" # One of unet, local, global
    num_channel_initial: 16 # Int type, number of initial channels in the network. Controls the network size.
    extract_levels: [0, 1, 2]
  loss:
    image:
      name: "gmi"
      weight: 1.0
    label:
      weight: 1.0
      name: "dice"
```

Moreover, one may want to specify several loss functions for each category. In that
case, a dashed line (-) indicates the specification of a new loss function under each
field.

E.G.

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
  backbone:
    name: "local" # One of unet, local, global
    num_channel_initial: 16 # Int type, number of initial channels in the network. Controls the network size.
    extract_levels: [0, 1, 2]
  loss:
    image:
      - name: "lncc"
        weight: 0.5
        kernel_size: 5
      - name: "gmi"
        weight: 0.5
    label:
      name: "dice"
      weight: 0.5
```

### Optimizer - required

The optimizer can be defined by using a `name` with other required arugment. The name
must be the same to the class name under `tf.keras.optimizers`.

For instance, to use a default
[Keras Adam optimizer](https://www.tensorflow.org/api_docs/python/tf/keras/optimizers/Adam),
the configuration should be

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
  backbone:
    name: "local" # One of unet, local, global
    num_channel_initial: 16 # Int type, number of initial channels in the network. Controls the network size.
    extract_levels: [0, 1, 2]
  loss:
    regularization:
      weight: 0.5 # weight of regularization loss
      name: "bending" # options include "bending", "gradient"
  optimizer:
    name: "Adam"
```

For a
[Keras SGD optimizer](https://www.tensorflow.org/api_docs/python/tf/keras/optimizers/SGD)
with learning rate 0.001, the configuration should be

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
  backbone:
    name: "local" # One of unet, local, global
    num_channel_initial: 16 # Int type, number of initial channels in the network. Controls the network size.
    extract_levels: [0, 1, 2]
  loss:
    regularization:
      weight: 0.5 # weight of regularization loss
      name: "bending" # options include "bending", "gradient"
  optimizer:
    name: "SGD"
    learning_rate: 0.001
```

### Preprocess - required

The `preprocess` field defines how the data loader feeds data into the model.

- `batch_size`: int, specifies the number of samples per step for prediction. If using
  multiple GPUs, i.e. `n` GPUs, each GPU will have mini batch size `batch_size / n`.
  Thus, `batch_size` should be divided by `n` evenly.
- `shuffle_buffer_num_batch`: int, helps define how much data should be pre-loaded into
  memory to buffer training, such that shuffle_buffer_size = batch_size \*
  shuffle_buffer_num_batch.
- `num_parallel_calls`: int, it defines the number of cpus used during preprocessing, -1
  means unlimited and it may take all cpus and significantly more memory.

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
  backbone:
    name: "local" # One of unet, local, global
    num_channel_initial: 16 # Int type, number of initial channels in the network. Controls the network size.
    extract_levels: [0, 1, 2]
  loss:
    regularization:
      weight: 0.5 # weight of regularization loss
      name: "bending" # options include "bending", "gradient"
  optimizer:
    name: "sgd"
    sgd:
      learning_rate: 1.0e-5
      momentum: 0.9
      nesterov: false
  preprocess:
    batch_size: 32
    shuffle_buffer_num_batch: 1
    num_parallel_calls: -1 # number elements to process asynchronously in parallel during preprocessing, -1 means unlimited, heuristically it should be set to the number of CPU cores available
```

### Epochs - required

The `epochs` field defines the number of epochs to train the network for.

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
  backbone:
    name: "local" # One of unet, local, global
    num_channel_initial: 16 # Int type, number of initial channels in the network. Controls the network size.
    extract_levels: [0, 1, 2]
  loss:
    regularization:
      weight: 0.5 # weight of regularization loss
      name: "bending" # options include "bending", "gradient"
  optimizer:
    name: "sgd"
    sgd:
      learning_rate: 1.0e-5
      momentum: 0.9
      nesterov: false
  preprocess:
    batch_size: 32
    shuffle_buffer_num_batch: 1
  epochs: 1000
```

### Saving frequency - required

The `save_period` field defines the save frequency - the model will be saved every
`save_period` epochs.

```yaml
train:
  method: "ddf" # One of ddf, dvf, conditional
  backbone:
    name: "local" # One of unet, local, global
    num_channel_initial: 16 # Int type, number of initial channels in the network. Controls the network size.
    extract_levels: [0, 1, 2]
  loss:
    regularization:
      weight: 0.5 # weight of regularization loss
      name: "bending" # options include "bending", "gradient"
  optimizer:
    name: "sgd"
    sgd:
      learning_rate: 1.0e-5
      momentum: 0.9
      nesterov: false
  preprocess:
    batch_size: 32
    shuffle_buffer_num_batch: 1
  epochs: 1000
  save_period: 5
```
---
title: 'DeepReg: a deep learning toolkit for medical image registration'
tags:
  - Python
  - TensorFlow
  - medical image registration
  - image fusion
  - deep learning
  - neural networks
authors:
  - name: Yunguan Fu
    orcid: 0000-0002-1184-7421
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
  - name: Nina Montaña Brown
    orcid: 0000-0001-5685-971X
    affiliation: "1, 2"
  - name: Shaheer U. Saeed
    orcid: 0000-0002-5004-0663
    affiliation: "1, 2"
  - name: Adrià Casamitjana
    orcid: 0000-0002-0539-3638
    affiliation: 2
  - name: Zachary M. C. Baum
    affiliation: "1, 2"
    orcid: 0000-0001-6838-335X
  - name: Rémi Delaunay
    orcid: 0000-0002-0398-4995
    affiliation: "1, 4"
  - name: Qianye Yang
    orcid: 0000-0003-4401-5311
    affiliation: "1, 2"
  - name: Alexander Grimwood
    orcid: 0000-0002-2608-2580
    affiliation: "1, 2"
  - name: Zhe Min
    orcid: 0000-0002-8903-1561
    affiliation: 1
  - name: Stefano B. Blumberg
    orcid: 0000-0002-7150-9918
    affiliation: 2
  - name: Juan Eugenio Iglesias
    orcid: 0000-0001-7569-173X
    affiliation: "2, 5, 6"
  - name: Dean C. Barratt
    orcid: 0000-0003-2916-655X
    affiliation: "1, 2"
  - name: Ester Bonmati
    orcid: 0000-0001-9217-5438
    affiliation: "1, 2"
  - name: Daniel C. Alexander
    orcid: 0000-0003-2439-350X
    affiliation: 2
  - name: Matthew J. Clarkson
    orcid: 0000-0002-5565-1252
    affiliation: "1, 2"
  - name: Tom Vercauteren
    orcid: 0000-0003-1794-0456
    affiliation: 4
  - name: Yipeng Hu
    orcid: 0000-0003-4902-0486
    affiliation: "1, 2"
affiliations:
 - name: Wellcome/EPSRC Centre for Surgical and Interventional Sciences, University College London, London, UK
   index: 1
 - name: Centre for Medical Image Computing, University College London, London, UK
   index: 2
 - name: InstaDeep, London, UK
   index: 3
 - name: Department of Surgical & Interventional Engineering, King’s College London, London, UK
   index: 4
 - name: Martinos Center for Biomedical Imaging, Massachusetts General Hospital and Harvard Medical School, Boston, USA
   index: 5
 - name: Computer Science and Artificial Intelligence Laboratory, Massachusetts Institute of Technology, Boston, USA
   index: 6
date: 1 September 2020
bibliography: paper.bib
---

# Summary
Image fusion is a fundamental task in medical image analysis and computer-assisted intervention. Medical image registration, computational algorithms that align different images together [@hill2001medical], has in recent years turned the research attention towards deep learning. Indeed, the representation ability to learn from population data with deep neural networks has opened new possibilities for improving registration generalisability by mitigating difficulties in designing hand-engineered image features and similarity measures for many real-world clinical applications [@haskins2020deep; @fu2020deep]. In addition, its fast inference can substantially accelerate registration execution for time-critical tasks.

*DeepReg* is a Python package using TensorFlow [@tensorflow2015-whitepaper] that implements multiple registration algorithms and a set of predefined _dataset loaders_, supporting both labelled- and unlabelled data. DeepReg also provides command-line tool options that enable basic and advanced functionalities for model training, prediction and image warping. These implementations, together with their documentation, tutorials and demos, aim to simplify workflows for prototyping and developing novel methodology, utilising latest development and accessing quality research advances. DeepReg is unit tested and a set of customised contributor guidelines are provided to facilitate community contributions.

A submission to the MICCAI Educational Challenge has utilised the DeepReg code and demos to explore the link between classical algorithms and deep-learning-based methods [@brown2020introduction], while a recently published research work investigated temporal changes in prostate cancer imaging, by using a longitudinal registration adapted from the DeepReg code [@yang2020longitudinal].

# Statement of need
Currently, popular packages focusing on deep learning methods for medical imaging, such as NiftyNet [@gibson2018niftynet] and MONAI (https://monai.io/), do not support image registration. The existing open-sourced registration projects either implement specific published algorithms without automated testing, such as the VoxelMorph [@balakrishnan2019voxelmorph], or focus on classical methods, such as NiftiReg [@modat2010fast], SimpleElastix [@marstal2016simpleelastix] and AirLab [@sandkuhler2018airlab]. Therefore an open-sourced project focusing on image registration with deep learning is much needed for general research and education purposes.

# Implementation
DeepReg implements a framework for unsupervised learning [@de2019deep; @balakrishnan2019voxelmorph], weakly-supervised learning [@hu2018label; @hu2018weakly] and their combinations and variants, e.g. [@hu2019conditional]. Many options are included for major components of these approaches, such as different image- and label dissimilarity functions, transformation models [@ashburner2007fast; @vercauteren2009diffeomorphic; @hill2001medical], deformation regularisation [@rueckert1999nonrigid] and different neural network architectures [@hu2018weakly; @he2016deep; @simonyan2014very]. Details of the implemented methods are described in the documentation. The provided dataset loaders adopt staged random sampling strategy to ensure unbiased learning from groups, images and labels [@hu2018weakly; @yang2020longitudinal]. These algorithmic components together with the flexible dataset loaders are building blocks of many other registration tasks, such as group-wise registration and morphological template construction [@dalca2019learning; @siebert2020deep; @luo2020mvmm].

# DeepReg Demos
In addition to the tutorials and documentation, DeepReg provides a collection of demonstrations, _DeepReg Demos_, using open-accessible data with real-world clinical applications.

## Paired images
Many clinical applications for tracking organ motion and other temporal changes require _intra-subject_ _single-modality_ image registration. Registering lung CT images for the same patient, acquired at expiratory and inspiratory phases [@hering_alessa_2020_3835682], is such an example of both unsupervised (without labels) and combined supervision (trained with additional label dissimilarity based on anatomical segmentation). Furthermore, registering prostate MR, acquired before surgery, and intra-operative ultrasound images is an example of weakly-supervised learning for multimodal image registration [@hu2018weakly]. Another DeepReg Demo illustrates MR-to-ultrasound image registration is to track tissue deformation and brain tumour resection during neurosurgery [@xiao2017resect].

## Unpaired images
Unpaired images are found in applications such as _single-modality_ _inter-subject_ registration. One demo registers different brain MR images from different subjects [@simpson2019large], fundamental to population studies. Two other applications align unpaired inter-subject CT images for lung [@hering_alessa_2020_3835682] and abdominal organs [@adrian_dalca_2020_3715652]. Additionally, the support for cross-validation in DeepReg has been included in a demo, which registers 3D ultrasound images from different prostate cancer patients.

## Grouped images
Unpaired images may also be grouped in applications such as _single-modality_ _intra-subject_ registration. In this case, each subject has multiple images acquired, for instance, at two or more time points. For demonstration, multi-sequence cardiac MR images, acquired from myocardial infarction patients [@zhuang2020cardiac], are registered, where multiple images within each subject are considered as grouped images. Prostate longitudinal MR registration is proposed to track the cancer progression during active surveillance programme [@yang2020longitudinal]. Using segmentation from this application, another demo application illustrates aligning intra-patient prostate gland masks - also an example of feature-based registration based on deep learning.

# Conclusion
DeepReg provides a collection of deep learning algorithms and dataset loaders to train image registration networks, which provides a reference of basic functionalities. In its permissible open-source format, DeepReg not only provides a tool for scientific research and higher education, but also welcomes contributions from wider communities.

# Acknowledgements

This work is supported by the Wellcome/EPSRC Centre for Interventional and Surgical Sciences (203145Z/16/Z). Support was also from the Engineering and Physical Sciences Research Council (EPSRC) (EP/M020533/1, NS/A000049/1), National Institute for Health Research University College London Hospitals Biomedical Research Centre and Wellcome Trust (203148/Z/16/Z). TV is supported by a Medtronic / Royal Academy of Engineering Research Chair (RCSRF1819\7\34). NMB, ZB, RD are also supported by the EPSRC CDT i4health (EP/S021930/1). ZB is supported by the Natural Sciences and Engineering Research Council of Canada Postgraduate Scholarships-Doctoral Program and the University College London Overseas and Graduate Research Scholarships.

# References
<!-- This will be filled in by references in paper.bib -->
.. include:: introduction.rst

.. toctree::
    :hidden:
    :maxdepth: 2
    :caption: Getting Started

    getting_started/install
    getting_started/quick_start

.. toctree::
    :hidden:
    :maxdepth: 2
    :caption: Tutorial

    tutorial/registration
    tutorial/cross_val
    tutorial/run_cluster

.. toctree::
    :hidden:
    :maxdepth: 2
    :caption: DeepReg Demo

    demo/introduction

.. toctree::
    :hidden:
    :maxdepth: 2
    :caption: Documentation

    docs/cli
    docs/logging
    docs/configuration
    docs/dataset_loader
    docs/registry
    docs/experimental

.. toctree::
    :hidden:
    :maxdepth: 2
    :caption: API Reference

    api/entry_point
    api/loader
    api/registry
    api/network
    api/backbone
    api/layer
    api/loss
    api/optimizer

.. toctree::
    :hidden:
    :maxdepth: 2
    :caption: Contributing to DeepReg

    contributing/guide
    contributing/test
    contributing/demo
    contributing/docs
    contributing/release
DeepReg
=======

DeepReg is a freely available, community-supported open-source toolkit
for research and education in medical image registration using deep
learning.

The current version is implemented as a `TensorFlow2`_-based framework,
and contains implementations for unsupervised- and weakly-supervised
algorithms with their combinations and variants. DeepReg has a practical
focus on growing and diverse clinical applications, as seen in the
provided examples - `DeepReg Demos`_.

`Get involved`_ and help make DeepReg better!


Features
--------

DeepReg extends and simplifies workflows for medical imaging researchers
working in TensorFlow 2, and can be easily installed and used for
efficient training and rapid deployment of deep-learning registration
algorithms.

DeepReg is designed to be used with minimal programming or scripting,
owing to its built-in command line tools.

Our development and all related work involved in the project is public,
and released under the Apache 2.0 license.


Contact
-------

For development matters, please `raise an issue`_.

For matters regarding the `Code of Conduct`_, such as a complaint,
please email the DeepReg Development Team: DeepRegNet@gmail.com.

Alternatively, please contact one or more members of the CoC Committee as appropriate: Nina Montana Brown (nina.brown.15@ucl.ac.uk), Ester Bonmati (e.bonmati@ucl.ac.uk), Matt Clarkson (m.clarkson@ucl.ac.uk).


.. image:: ../asset/weiss.jpg
    :width: 300
    :alt: WEISS Logo


.. image:: ../asset/medicalengineering.svg
    :width: 250
    :alt: CME Logo

.. _TensorFlow2: https://www.tensorflow.org/
.. _DeepReg Demos: https://deepreg.readthedocs.io/en/latest/demo/introduction.html
.. _Get involved: https://deepreg.readthedocs.io/en/latest/contributing/guide.html
.. _WEISS: https://www.ucl.ac.uk/interventional-surgical-sciences/
.. _CME: https://medicalengineering.org.uk/
.. _Code of Conduct: https://github.com/DeepRegNet/DeepReg/blob/main/docs/CODE_OF_CONDUCT.md
.. _raise an issue: https://github.com/DeepRegNet/DeepReg/issues/new


Contributors
------------

DeepReg is maintained by a team of developers and researchers.
People with significant contributions to DeepReg are acknowledged in the `Contributor List`_.

This open-source initiative started within University College London,
with support from the Wellcome/EPSRC Centre for Interventional and
Surgical Sciences (`WEISS`_), and partial support from the
Wellcome/EPSRC Centre for Medical Engineering (`CME`_).

.. _Contributor List: https://github.com/DeepRegNet/DeepReg/blob/main/docs/Contributors.md
.. mdinclude:: ../../../demos/unpaired_us_prostate_cv/README.md
.. mdinclude:: ../../../demos/grouped_mr_heart/README.md
.. mdinclude:: ../../../demos/unpaired_ct_lung/README.md
.. mdinclude:: ../../../demos/paired_mrus_prostate/README.md
.. mdinclude:: ../../../demos/paired_ct_lung/README.md
.. mdinclude:: ../../../demos/classical_mr_prostate_nonrigid/README.md
.. mdinclude:: ../../../demos/classical_ct_headneck_affine/README.md
.. mdinclude:: ../../../demos/unpaired_ct_abdomen/README.md
.. mdinclude:: ../../../demos/grouped_mask_prostate_longitudinal/README.md
.. mdinclude:: ../../../demos/unpaired_mr_brain/README.md
Introduction to DeepReg Demos
=============================

DeepReg offers multiple built-in dataset loaders to support real-world clinical scenarios, in which images may be paired, unpaired or grouped. Images may also be labeled with segmented regions of interest to assist registration.

A typical workflow to develop a `registration network`_ using DeepReg
includes:

- Select a dataset loader, among the `unpaired, paired and grouped`_,
  and prepare data into folders as required;
- Configure the network training in the configuration yaml file(s), as
  specified in `supported configuration details`_;
- Train and tune the registration network with the `command line tool`_
  ``deepreg_train``;
- Test or use the final trained registration network with the `command line tool`_
  ``deepreg_predict``.

Besides the tutorials, a series of DeepReg Demos are provided to
showcase a wide range of applications with real clinical image and label data.
These applications range from ultrasound, CT and MR images,
covering many clinical specialties such as neurology, urology,
gastroenterology, oncology, respiratory and cardiovascular diseases.

Each DeepReg Demo provides a step-by-step instruction
to explain how different scenarios can be implemented with DeepReg.
All data sets used are open-accessible.
Pre-trained models with numerical and graphical inference results are also available.

.. _demo-disclaimer:

.. note::

   DeepReg Demos are provided to demonstrate functionalities in DeepReg.
   Although effort has been made to ensure these demos are representative
   of real-world applications, the implementations and the results are not
   peer-reviewed or tested for clinical efficacy. Substantial further
   adaptation and development may be required for any potential clinical
   adoption.

.. _registration network: ../tutorial/registration.html
.. _unpaired, paired and grouped: ../docs/dataset_loader.html
.. _supported configuration details: ../docs/configuration.html
.. _command line tool: ../docs/cli.html

Paired Images
=============

The following DeepReg Demos provide examples of
using paired images.

- `Paired lung CT registration <paired_ct_lung.html>`__

  This demo registers paired CT lung images, with optional weak supervision.

- `Paired brain MR-ultrasound registration <paired_mrus_brain.html>`__

  This demo registers paired preoperative MR images and 3D tracked ultrasound images for
  locating brain tumours during neurosurgery, with optional weak supervision.

- `Paired prostate MR-ultrasound registration <paired_mrus_prostate.html>`__

  This demo registers paired MR-to-ultrasound prostate images, an example of
  weakly-supervised multimodal image registration.

.. toctree::
    :glob:
    :hidden:
    :maxdepth: 2

    paired_*

Unpaired Images
===============

The following DeepReg Demos provide examples of
using unpaired images.

- `Unpaired abdominal CT registration <unpaired_ct_abdomen.html>`__

  This demo compares three training strategies, using unsupervised, weakly-supervised and
  combined losses, to register inter-subject abdominal CT images.

- `Unpaired lung CT registration <unpaired_ct_lung.html>`__

  This demo registers unpaired CT lung images, with optional weak supervision.

- `Unpaired hippocampus MR registration <unpaired_mr_brain.html>`__

  This demo aligns hippocampus on MR images between different patients, with optional weak
  supervision.

- `Unpaired prostate ultrasound registration <unpaired_us_prostate_cv.html>`__

  This demo registers 3D ultrasound images with a 9-fold cross-validation. This strategy
  is applicable for any of the available dataset loaders.

.. toctree::
    :glob:
    :hidden:
    :maxdepth: 2

    unpaired_*

Grouped Images
==============

The following DeepReg Demos provide examples of using grouped images.

- `Pairwise registration for grouped prostate segmentation masks <grouped_mask_prostate_longitudinal.html>`__

  This demo registers grouped masks (as input images) of prostate glands from MR images,
  an example of feature-based registration.

- `Pairwise registration for grouped cardiac MR images <grouped_mr_heart.html>`__

  This demo registers grouped CMR images, where each group has multi-sequence CMR images from a single patient.

.. toctree::
    :glob:
    :hidden:
    :maxdepth: 2

    grouped_*

Classical Registration
======================

The following DeepReg Demos provide examples of
using classical registration methods.

- `Classical affine registration for head-and-neck CT images <classical_ct_headneck_affine.html>`__

  This demo registers head-and-neck CT images using iterative affine registration.

- `Classical nonrigid registration for prostate MR images <classical_mr_prostate_nonrigid.html>`__

  This demo registers prostate MR images using iterative nonrigid registration.

.. toctree::
    :hidden:
    :glob:
    :maxdepth: 2

    classical_*
.. mdinclude:: ../../../demos/paired_mrus_brain/README.md
Layer
=====

Layer
-----

Activation
^^^^^^^^^^

.. autoclass:: deepreg.model.layer.Activation
    :members:

AdditiveUpSampling
^^^^^^^^^^^^^^^^^^

.. autoclass:: deepreg.model.layer.AdditiveUpSampling
    :members:

Conv3d
^^^^^^

.. autoclass:: deepreg.model.layer.Conv3d
    :members:

Conv3dBlock
^^^^^^^^^^^

.. autoclass:: deepreg.model.layer.Conv3dBlock
    :members:

Conv3dWithResize
^^^^^^^^^^^^^^^^

.. autoclass:: deepreg.model.layer.Conv3dWithResize
    :members:

Deconv3d
^^^^^^^^

.. autoclass:: deepreg.model.layer.Deconv3d
    :members:

Deconv3dBlock
^^^^^^^^^^^^^

.. autoclass:: deepreg.model.layer.Deconv3dBlock
    :members:

Dense
^^^^^

.. autoclass:: deepreg.model.layer.Dense
    :members:

DownSampleResnetBlock
^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: deepreg.model.layer.DownSampleResnetBlock
    :members:

IntDVF
^^^^^^

.. autoclass:: deepreg.model.layer.IntDVF
    :members:

LocalNetResidual3dBlock
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: deepreg.model.layer.LocalNetResidual3dBlock
    :members:

LocalNetUpSampleResnetBlock
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: deepreg.model.layer.LocalNetUpSampleResnetBlock
    :members:

MaxPool3d
^^^^^^^^^

.. autoclass:: deepreg.model.layer.MaxPool3d
    :members:

Norm
^^^^

.. autoclass:: deepreg.model.layer.Norm
    :members:

Residual3dBlock
^^^^^^^^^^^^^^^

.. autoclass:: deepreg.model.layer.Residual3dBlock
    :members:

UpSampleResnetBlock
^^^^^^^^^^^^^^^^^^^

.. autoclass:: deepreg.model.layer.UpSampleResnetBlock
    :members:

Warping
^^^^^^^

.. autoclass:: deepreg.model.layer.Warping
    :members:

Util
----

.. automodule:: deepreg.model.layer_util
    :members:
Registry
========

.. autoclass:: deepreg.registry.Registry
   :members:
   :private-members:
Loss
====

Image Loss
----------

.. automodule:: deepreg.loss.image
    :members:

Label Loss
----------

.. automodule:: deepreg.loss.label
    :members:

Deformation Loss
----------------

.. automodule:: deepreg.loss.deform
    :members:

Loss Util
---------

.. automodule:: deepreg.loss.util
    :members:
Network
=======

.. automodule:: deepreg.model.network
    :members:
Entry Point
===========

Train
-----

.. automodule:: deepreg.train
    :members:

Predict
-------

.. automodule:: deepreg.predict
    :members:

Warp
----

.. automodule:: deepreg.warp
    :members:
Optimizer
=========

.. automodule:: deepreg.model.optimizer
    :members:
Backbone
========

Interface
---------

.. autoclass:: deepreg.model.backbone.interface.BackboneInterface
    :members:

Local Net
---------

.. autoclass:: deepreg.model.backbone.local_net.LocalNet
    :members:

Global Net
----------

.. autoclass:: deepreg.model.backbone.global_net.GlobalNet
    :members:

U-Net
-----

.. autoclass:: deepreg.model.backbone.u_net.UNet
    :members:
Dataset Loader
==============

Paired Loader
-------------

.. automodule:: deepreg.dataset.loader.paired_loader
    :members:

Unpaired Loader
---------------

.. automodule:: deepreg.dataset.loader.unpaired_loader
    :members:

Grouped Loader
--------------

.. automodule:: deepreg.dataset.loader.grouped_loader
    :members:

File Loader
===========

Interface
---------

.. autoclass:: deepreg.dataset.loader.interface.FileLoader
    :members:

Nifti Loader
------------

.. automodule:: deepreg.dataset.loader.nifti_loader
    :members:

H5 Loader
---------

.. automodule:: deepreg.dataset.loader.h5_loader
    :members:
Installation
============

DeepReg can be installed in Python 3.7 and external python dependencies are mainly defined in `requirements`_.
DeepReg primarily supports and is regularly tested with Ubuntu and Mac OS.

There are multiple different methods to install DeepReg:

1. Clone `DeepReg`_ and create a virtual environment using `Anaconda`_ / `Miniconda`_ (**recommended**).
2. Clone `DeepReg`_ and build a docker image using the provided docker file.
3. Install directly from PyPI release without cloning `DeepReg`_.

Install via Conda
-----------------

The recommended method is to install DeepReg in a dedicated virtual
environment using `Anaconda`_ / `Miniconda`_.

Please clone `DeepReg`_ first and change current directory to the DeepReg root directory:

.. code:: bash

    git clone https://github.com/DeepRegNet/DeepReg.git
    cd DeepReg

Then, install or update the conda environment following the instructions below.
Please see the `official conda documentation <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`__
for more details.

.. tabs::

    .. tab:: Linux

        Install prerequisites (Optional).

        .. code:: bash

            sudo apt-get update
            sudo apt-get install graphviz

        Install DeepReg without GPU support.

        .. code:: bash

            conda env create -f environment_cpu.yml
            conda activate deepreg

        Install DeepReg with GPU support.

        .. code:: bash

            conda env create -f environment.yml
            conda activate deepreg

    .. tab:: Mac OS

        Install prerequisites (Optional).

        .. code:: bash

            brew install graphviz

        Install DeepReg without GPU support.

        .. code:: bash

            conda env create -f environment_cpu.yml
            conda activate deepreg

        Install DeepReg with GPU support.

        .. warning::

            Not supported or tested.

    .. tab:: Windows

        Install DeepReg without GPU support.

        .. warning::

            DeepReg on Windows is not fully supported.
            However, you can use the `Windows Subsystem for Linux <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`__.
            Set up WSL and follow the DeepReg setup instructions for Linux.

        Install DeepReg with GPU support.

        .. warning::

            Not supported or tested.


After activating the conda environment, please install DeepReg locally:

.. code:: bash

    pip install -e .

Install via docker
------------------

We also provide the docker file for building the docker image.
Please clone `DeepReg`_ repository first:

.. code:: bash

    git clone https://github.com/DeepRegNet/DeepReg.git

Then, install DeepReg following the instructions below.

Install docker
^^^^^^^^^^^^^^

Docker can be installed following the `official documentation <https://docs.docker.com/get-docker/>`__.

For Linux based OS, there are some `additional setup <https://docs.docker.com/engine/install/linux-postinstall/>`__ after the installation.
Otherwise you might have permission errors.

Build docker image
^^^^^^^^^^^^^^^^^^

.. code:: bash

    docker build . -t deepreg -f Dockerfile

where

- :code:`-t` names the built image as :code:`deepreg`.
- :code:`-f` provides the docker file for configuration.

Create a container
^^^^^^^^^^^^^^^^^^

.. code:: bash

    docker run --name <container_name> --privileged=true -ti deepreg bash

where
- :code:`--name` names the created container.
- :code:`--privileged=true` is required to solve the permission issue linked to TensorFlow profiler.
- :code:`-it` allows interaction with container and enters the container directly,
check more info on `stackoverflow <https://stackoverflow.com/questions/48368411/what-is-docker-run-it-flag>`__.

Remove a container
^^^^^^^^^^^^^^^^^^

.. code:: bash

    docker rm -v <container_name>

which removes a created container and its volumes, check more info on `docker documentation <https://docs.docker.com/engine/reference/commandline/rm/)>`__.

Install via PyPI
----------------

Please use the following command to install DeepReg directly from the PyPI release:

.. code:: bash

    pip install deepreg

The PyPI release currently does not ship with test data and demos.
Running examples, such as those in `Quick Start`_ and `DeepReg Demo`_,
in this documentation may require downloading additional test data.

Once you have installed DeepReg via :code:`pip`, you can run the following
command to download the necessary files to run all examples by:

.. code:: bash

    deepreg_download

The above will download the files to the current working directory.
If you need to download to a specific directory, use the
:code:`--output_dir` or :code:`-d` flag to specify this.

**Note**

1. All dependencies, APIs and command-line tools will be installed automatically via each installation method.
2. Only released versions of DeepReg are available via PyPI release.
   Therefore it is different from the `latest (unstable) version <https://github.com/DeepRegNet/DeepReg>`__ on GitHub.

.. _Quick Start: quick_start.html
.. _DeepReg Demo: ../demo/introduction.html
.. _Anaconda: https://docs.anaconda.com/anaconda/install
.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _DeepReg: https://github.com/DeepRegNet/DeepReg
.. _requirements: https://github.com/DeepRegNet/DeepReg/blob/main/requirements.txt
Experimental Features
=====================

DeepReg provides some experimental features.
These are still in development with variable levels of readiness.

The following tutorials provide an overview of these features. To submit feedback, open a `new issue <https://github.com/DeepRegNet/DeepReg/issues/new>`__.

-  `Label sampling`_

.. _Label sampling: exp_label_sampling.html

.. toctree::
    :hidden:
    :maxdepth: 2

    exp_label_sampling
