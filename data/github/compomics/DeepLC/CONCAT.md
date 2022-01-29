# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2022-01-21
- Bug fix duplicate peptide+mod for DeepCALLC
- New feature DeepCALLC

## [0.1.39] - 2022-01-12
- Version bump

## [0.1.38] - 2022-01-10
- Deep(CAL)LC functionality

## [0.1.37] - 2021-09-09
- Pygam as calibration function

## [0.1.36] - 2021-09-09
- Update to Streamlit webserver: Use `st.form` and new official download button

## [0.1.35] - 2021-08-09
- More elegant solution for library call as global

## [0.1.34] - 2021-07-29
- Fix var call dict object

## [0.1.33] - 2021-07-29
- Fix library delete call

## [0.1.32] - 2021-07-29
- Temporary fix suggested by markmipt for the library

## [0.1.31] - 2021-07-01
- Change logging to specified logger object instead of standard logger

## [0.1.30] - 2021-06-09
- GUI: Fix small font through starting jar with cmd to increase font size
- Added testing for Python 3.9
- Relax h5py requirement to allow v3
- Fixed GitHub Action workflow for Streamlit docker image build

## [0.1.29] - 2021-03-24
- Bug in writing library where a list was assumed so library only partially filled

## [0.1.28] - 2021-03-17
- Make it optional to reload library

## [0.1.27] - 2021-03-16
- Ignore library messages

## [0.1.26] - 2021-03-15
- Force older library h5py for compatability

## [0.1.25] - 2021-03-12
- Library gets appended for non-calibration peptides too

## [0.1.24] - 2021-03-12
- Change the default windows install to pip instead of conda
- Change reporting of identifiers used from library

## [0.1.23] - 2021-03-04
- Log the amount of identifiers in library used

## [0.1.22] - 2021-03-04
- Publish PyPI and GitHub release

## [0.1.21] - 2021-03-04
- Add library functionality that allows for storing and retrieving predictions (without running the model)

## [0.1.20] - 2021-02-19
- Describe hyperparameters and limit CPU threads
- Additional modfications, including those exclusive to pFind
- Change calibration error to warning (since it is a warning if it is out of range...)

## [0.1.18] - 2021-01-11
- Limit CPU usage by tensorflow by connecting to n_jobs

## [0.1.17] - 2020-08-07
- Support for Python 3.8

## [0.1.16] - 2020-05-18
- Bug fix in the calibration function
- Had to order the predicted instead of observed retention times of the calibration analytes
- Thanks to @courcelm for both finding and fixing the issue

## [0.1.15] - 2020-05-15
- Different calibration function, should not contain gaps anymore
- Changed to more accurate rounding
- Changed to splitting in groups of retention time instead of groups of peptides

## [0.1.14] - 2020-02-21
- Changed default model to a different data set

## [0.1.13] - 2020-02-21
- Duplicate peptides charge fix

## [0.1.12] - 2020-02-15
- Support for charges and spaces in peprec

## [0.1.11] - 2020-02-13
- Fixes in GUI

## [0.1.10] - 2020-02-10
- Include less models in package to meet PyPI 60MB size limitation

## [0.1.9] - 2020-02-09
- Bugfix: Pass custom activation function

## [0.1.8] - 2020-02-07
- Fixed support for averaging predictions of groups of models (ensemble) when no models were passed
- New models for ensemble

## [0.1.7] - 2020-02-07
- Support for averaging predictions of groups of models (ensemble)

## [0.1.6] - 2020-01-21
- Fix the latest release

## [0.1.5] - 2020-01-21
- Spaces in paths to files and installation allowed
- References to other CompOmics tools removed in GUI

## [0.1.5] - 2020-02-13
- Fixes in GUI

## [0.1.4] - 2020-01-17
- Fix the latest release

## [0.1.3] - 2020-01-17
- Fixed the .bat installer (now uses bioconda)

## [0.1.2] - 2019-12-19
- Example files in GUI folder
- Unnecesary bat and sh for running GUI removed

## [0.1.1] - 2019-12-18
- Switch to setuptools
- Reorder publish workflow; build wheels

## [0.1.1.dev8] - 2019-12-16
- Remove xgboost dependancy

## [0.1.1.dev7] - 2019-12-16
- Use dot instead of dash in versioning for bioconda

## [0.1.1-dev6] - 2019-12-15
- Fix publish action branch specification (2)

## [0.1.1-dev5] - 2019-12-15
- Fix publish action branch specification

## [0.1.1-dev4] - 2019-12-15
- Test other trigger for publish action

## [0.1.1-dev3] - 2019-12-15
- Update documentation, specify branch in publish action

## [0.1.1-dev2] - 2019-12-15
- Add long description to setup.py

## [0.1.1-dev1] - 2019-12-15
- Initial pre-release
<img src="https://github.com/compomics/DeepLC/raw/master/img/deeplc_logo.png"
width="150" height="150" /> <br/><br/>

[![GitHub release](https://flat.badgen.net/github/release/compomics/deeplc)](https://github.com/compomics/DeepLC/releases/latest/)
[![PyPI](https://flat.badgen.net/pypi/v/deeplc)](https://pypi.org/project/deeplc/)
[![Conda](https://img.shields.io/conda/vn/bioconda/deeplc?style=flat-square)](https://bioconda.github.io/recipes/deeplc/README.html)
[![GitHub Workflow Status](https://flat.badgen.net/github/checks/compomics/deeplc/)](https://github.com/compomics/deeplc/actions/)
[![License](https://flat.badgen.net/github/license/compomics/deeplc)](https://www.apache.org/licenses/LICENSE-2.0)
[![Twitter](https://flat.badgen.net/twitter/follow/compomics?icon=twitter)](https://twitter.com/compomics)

DeepLC: Retention time prediction for (modified) peptides using Deep Learning.

---

- [Introduction](#introduction)
- [Usage](#usage)
  - [Web application](#web-application)
  - [Graphical user interface](#graphical-user-interface)
  - [Python package](#python-package)
    - [Installation](#installation)
    - [Command line interface](#command-line-interface)
    - [Python module](#python-module)
  - [Input files](#input-files)
  - [Prediction models](#prediction-models)
- [Citation](#citation)
- [Q&A](#qa)

---

## Introduction

DeepLC is a retention time predictor for (modified) peptides that employs Deep
Learning. Its strength lies in the fact that it can accurately predict
retention times for modified peptides, even if hasn't seen said modification
during training.

DeepLC can be used through the
[web application](https://iomics.ugent.be/deeplc/),
locally with a graphical user interface (GUI), or as a Python package. In the
latter case, DeepLC can be used from the command line, or as a Python module.

## Usage

### Web application
[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://iomics.ugent.be/deeplc/)

Just go to [iomics.ugent.be/deeplc](https://iomics.ugent.be/deeplc/) and get started!


### Graphical user interface
#### Installation

[![Download GUI](https://flat.badgen.net/badge/download/GUI/green)](https://github.com/compomics/DeepLC/releases/latest/)

1. Download `deeplc_gui.zip` from the
[latest release](https://github.com/compomics/DeepLC/releases/latest/) and
unzip.
2. Install DeepLC GUI with `install_gui_windows.bat` or `install_gui_linux.sh`,
depending on your operating system.
3. Run DeepLC GUI by running the `deeplc_gui.jar`.

#### Run parameters

Dictionary divider - this parameter defines the precision to use for fast-lookup of retention times for calibration. A value of 10 means a precision of 0.1 (and 100 a precision of 0.01) between the calibration anchor points. This parameter does not influence the precision of the calibration, but setting it too high might mean that there is bad selection of the models between anchor points. A safe value is usually higher than 10.

Split calibration - the number of divisions for the chromatogram. If the value is set to 10 the chromatogram is split up into 10 equidistant parts. For each part the median value of the calibration peptides is selected. These are the anchor points. Between each anchor point a linear fit is made.

Batch number - define the number of peptides to make predictions for in a single go (reduce to fit into memory, increase for faster prediction speeds).

### Python package

#### Installation

[![install with bioconda](https://flat.badgen.net/badge/install%20with/bioconda/green)](http://bioconda.github.io/recipes/deeplc/README.html)
[![install with pip](https://flat.badgen.net/badge/install%20with/pip/green)](http://bioconda.github.io/recipes/deeplc/README.html)
[![container](https://flat.badgen.net/badge/pull/biocontainer/green)](https://quay.io/repository/biocontainers/deeplc)

Install with conda, using the bioconda and conda-forge channels:
`conda install -c bioconda -c conda-forge deeplc`

Or install with pip:
`pip install deeplc`

#### Command line interface

To use the DeepLC CLI, run:

```sh
deeplc --file_pred <path/to/peptide_file.csv>
```

We highly recommend to add a peptide file with known retention times for
calibration:

```sh
deeplc --file_pred  <path/to/peptide_file.csv> --file_cal <path/to/peptide_file_with_tr.csv>
```

For an overview of all CLI arguments, run `deeplc --help`.

#### Python module

Minimal example:

```python
import pandas as pd
from deeplc import DeepLC

peptide_file = "datasets/test_pred.csv"
calibration_file = "datasets/test_train.csv"

pep_df = pd.read_csv(peptide_file, sep=",")
pep_df['modifications'] = pep_df['modifications'].fillna("")

cal_df = pd.read_csv(calibration_file, sep=",")
cal_df['modifications'] = cal_df['modifications'].fillna("")

dlc = DeepLC()
dlc.calibrate_preds(seq_df=cal_df)
preds = dlc.make_preds(seq_df=pep_df)
```

For a more elaborate example, see
[examples/deeplc_example.py](https://github.com/compomics/DeepLC/blob/master/examples/deeplc_example.py)
.

### Input files

DeepLC expects comma-separated values (CSV) with the following columns:

- `seq`: unmodified peptide sequences
- `modifications`: MS2PIP-style formatted modifications: Every modification is
  listed as `location|name`, separated by a pipe (`|`) between the location, the
  name, and other modifications. `location` is an integer counted starting at 1
  for the first AA. 0 is reserved for N-terminal modifications, -1 for
  C-terminal modifications. `name` has to correspond to a Unimod (PSI-MS) name.
- `tr`: retention time (only required for calibration)

For example:

```csv
seq,modifications,tr
AAGPSLSHTSGGTQSK,,12.1645
AAINQKLIETGER,6|Acetyl,34.095
AANDAGYFNDEMAPIEVKTK,12|Oxidation|18|Acetyl,37.3765
```

See
[examples/datasets](https://github.com/compomics/DeepLC/tree/master/examples/datasets)
for more examples.

### Prediction models

DeepLC comes with multiple CNN models trained on data from various experimental
settings:

| Model filename | Experimental settings | Publication |
| - | - | - |
| full_hc_dia_fixed_mods.hdf5 | Reverse phase | [Rosenberger et al. 2014](https://doi.org/10.1038/sdata.2014.31) |
| full_hc_LUNA_HILIC_fixed_mods.hdf5 | HILIC | [Spicer et al. 2018](https://doi.org/10.1016/j.chroma.2017.12.046) |
| full_hc_LUNA_SILICA_fixed_mods.hdf5 | HILIC | [Spicer et al. 2018](https://doi.org/10.1016/j.chroma.2017.12.046) |
| full_hc_PXD000954_fixed_mods.hdf5 | Reverse phase | [Rosenberger et al. 2014](https://doi.org/10.1038/sdata.2014.31) |

By default, DeepLC selects the best model based on the calibration dataset. If
no calibration is performed, the first default model is selected. Always keep
note of the used models and the DeepLC version.

The table above is for an old version of DeepLC, the current version comes with:

| Model filename | Experimental settings | Publication |
| - | - | - |
| full_hc_hela_hf_psms_aligned_1fd8363d9af9dcad3be7553c39396960.hdf5 | Reverse phase | [Kelstrup et al. 2018](https://doi.org/10.1021/acs.jproteome.7b006021) |
| full_hc_hela_hf_psms_aligned_8c22d89667368f2f02ad996469ba157e.hdf5 | Reverse phase | [Kelstrup et al. 2018](https://doi.org/10.1021/acs.jproteome.7b00602) |
| full_hc_hela_hf_psms_aligned_cb975cfdd4105f97efa0b3afffe075cc.hdf5 | Reverse phase | [Kelstrup et al. 2018](https://doi.org/10.1021/acs.jproteome.7b00602) |
| full_hc_PXD005573_mcp_cb975cfdd4105f97efa0b3afffe075cc.hdf5 | Reverse phase | [Bruderer et al. 2017](https://pubmed.ncbi.nlm.nih.gov/29070702/) |

For all the full models that can be used in DeepLC (including some TMT models!) please see:

[https://github.com/RobbinBouwmeester/DeepLCModels](https://github.com/RobbinBouwmeester/DeepLCModels)

## Citation

If you use DeepLC for your research, please use the following citation:
>**DeepLC can predict retention times for peptides that carry as-yet unseen modifications**  
>Robbin Bouwmeester, Ralf Gabriels, Niels Hulstaert, Lennart Martens, Sven Degroeve  
>bioRxiv 2020.03.28.013003; [doi: 10.1101/2020.03.28.013003](https://doi.org/10.1101/2020.03.28.013003)


## Q&A

**__Q: Is it required to indicate fixed modifications in the input file?__**

Yes, even modifications like carbamidomethyl should be in the input file.

**__Q: So DeepLC is able to predict the retention time for any modification?__**

Yes, DeepLC can predict the retention time of any modification. However, if the
modification is **very** different from the peptides the model has seen during
training the accuracy might not be satisfactory for you. For example, if the model
has never seen a phosphor atom before, the accuracy of the prediction is going to
be low.

**__Q: Installation fails. Why?__**

Please make sure to install DeepLC in a path that does not contain spaces. Run
the latest LTS version of Ubuntu or Windows 10. Make sure you have enough disk
space available, surprisingly TensorFlow needs quite a bit of disk space. If
you are still not able to install DeepLC, please feel free to contact us:

Robbin.Bouwmeester@ugent.be and Ralf.Gabriels@ugent.be

**__Q: I have a special usecase that is not supported. Can you help?__**

Ofcourse, please feel free to contact us:

Robbin.Bouwmeester@ugent.be and Ralf.Gabriels@ugent.be

**__Q: DeepLC runs out of memory. What can I do?__**

You can try to reduce the batch size. DeepLC should be able to run if the batch size is low
enough, even on machines with only 4 GB of RAM.

**__Q: I have a graphics card, but DeepLC is not using the GPU. Why?__**

For now DeepLC defaults to the CPU instead of the GPU. Clearly, because you want
to use the GPU, you are a power user :-). If you want to make the most of that expensive
GPU, you need to change or remove the following line (at the top) in __deeplc.py__:

```
# Set to force CPU calculations
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
```

Also change the same line in the function __reset_keras()__:

```
# Set to force CPU calculations
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
```

Either remove the line or change to (where the number indicates the number of GPUs):

```
# Set to force CPU calculations
os.environ['CUDA_VISIBLE_DEVICES'] = '1'
```

**__Q: What modification name should I use?__**

The names from unimod are used. The PSI-MS name is used by default, but the Interim name
is used as a fall-back if the PSI-MS name is not available. Please also see __unimod_to_formula.csv__
in the folder __unimod/__ for the naming of specific modifications.

**__Q: I have a modification that is not in unimod. How can I add the modification?__**

In the folder __unimod/__ there is the file __unimod_to_formula.csv__ that can be used to
add modifications. In the CSV file add a name (**that is unique and not present yet**) and
the change in atomic composition. For example:

```
Met->Hse,O,H(-2) C(-1) S(-1)
```

Make sure to use negative signs for the atoms subtracted.

**__Q: Help, all my predictions are between [0,10]. Why?__**

It is likely you did not use calibration. No problem, but the retention times for training
purposes were normalized between [0,10]. This means that you probably need to adjust the
retention time yourselve after analysis or use a calibration set as the input.

**__Q: How does the ensemble part of DeepLC work?__**

Models within the same directory are grouped if they overlap in their name. The overlap
has to be in their full name, except for the last part of the name after a "_"-character.

The following models will be grouped:

```
full_hc_dia_fixed_mods_a.hdf5
full_hc_dia_fixed_mods_b.hdf5
```

None of the following models will not be grouped:

```
full_hc_dia_fixed_mods2_a.hdf5
full_hc_dia_fixed_mods_b.hdf5
full_hc_dia_fixed_mods_2_b.hdf5
```

**__Q: I would like to take the ensemble average of multiple models, even if they are trained on different datasets. How can I do this?__**

Feel free to experiment! Models within the same directory are grouped if they overlap in
their name. The overlap has to be in their full name, except for the last part of the
name after a "_"-character.

The following models will be grouped:

```
model_dataset1.hdf5
model_dataset2.hdf5
```

So you just need to rename your models.
# Contributing

This document briefly describes how to contribute to
[DeepLC](https://github.com/compomics/DeepLC).

## Before you begin

If you have an idea for a feature, use case to add or an approach for a bugfix,
it is best to communicate with the community by creating an issue in
[GitHub issues](https://github.com/compomics/DeepLC/issues).

## How to contribute

- Fork [DeepLC](https://github.com/compomics/DeepLC) on GitHub to
make your changes.
- Commit and push your changes to your
[fork](https://help.github.com/articles/pushing-to-a-remote/).
- Open a
[pull request](https://help.github.com/articles/creating-a-pull-request/)
with these changes. You pull request message ideally should include:
   - A description of why the changes should be made.
   - A description of the implementation of the changes.
   - A description of how to test the changes.
- The pull request should pass all the continuous integration tests which are
  automatically run by
  [GitHub Actions](https://github.com/compomics/DeepLC/actions).


## Development workflow

- When a new version is ready to be published:

    1. Change the version number in `setup.py` using
    [semantic versioning](https://semver.org/).
    2. Update the changelog (if not already done) in `CHANGELOG.md` according to
    [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
    3. Set a new tag with the version number, e.g. `git tag v0.1.5`.
    4. Push to GitHub, with the tag: `git push; git push --tags`.

- When a new tag is pushed to (or made on) GitHub that matches `v*`, the
following GitHub Actions are triggered:

    1. The Python package is build and published to PyPI.
    2. A zip archive is made of the `./deeplc_gui/` directory, excluding
    `./deeplc_gui/src` with
    [Zip Release](https://github.com/marketplace/actions/zip-release).
    3. A GitHub release is made with the zipped GUI files as assets and the new
    changes listed in `CHANGELOG.md` with
    [Git Release](https://github.com/marketplace/actions/git-release).
    4. After some time, the bioconda package should get updated automatically.
# DeepLC Streamlit-based web server


## Usage

Pull and run with

```sh
docker run -p 8501 ghcr.io/compomics/deeplc-streamlit
```

Streamlit can be further configured using environment variables:

```sh
docker run \
    -p 8501 \
    -e STREAMLIT_SERVER_MAX_UPLOAD_SIZE=200 \
    ghcr.io/compomics/deeplc-streamlit
```
See
[Streamlit configuration](https://docs.streamlit.io/en/stable/streamlit_configuration.html)
for more info.
