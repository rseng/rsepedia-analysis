# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## Changed

- Minor linting [#93](https://github.com/matchms/ms2deepscore/pull/93)

## Fixed

- Handled numby dependency issues [#94](https://github.com/matchms/ms2deepscore/issues/94) and [#95](https://github.com/matchms/ms2deepscore/issues/95)

## [0.2.2] - 2021-08-19

## Fixed

- now compatible with new Tensorflow 2.6, also checked by additional CI runs for Tensorflow 2.4, 2.5 and 2.6 [#92](https://github.com/matchms/ms2deepscore/pull/92)

## [0.2.1] - 2021-07-20

## Changed

- Speed improvement of spectrum binning step [#90](https://github.com/matchms/ms2deepscore/pull/90)

## [0.2.0] - 2021-04-01

## Added

- `MS2DeepScoreMonteCarlo` Monte-Carlo dropout based ensembling do obtain mean/median score and STD [#65](https://github.com/matchms/ms2deepscore/pull/65)
- choice between `median` (default) and `mean` ensemble score which come with `IQR` and `STD` as uncertainty measures [#86](https://github.com/matchms/ms2deepscore/pull/86)
- `dropout_in_first_layer` option for SiameseModel (default is False) [#86](https://github.com/matchms/ms2deepscore/pull/86)
- `use_fixed_set` option for data generators to create deterministic training/testing data with fixed random seed [#73](https://github.com/matchms/ms2deepscore/issues/73)

## Changed

- small update of `create_histograms_plot` to make the plot prettier/better to read [#85](https://github.com/matchms/ms2deepscore/pull/85)

## Fixed

- solved minor unclarity with the pair selection for non-available reference scores [#79](https://github.com/matchms/ms2deepscore/pull/79)
- solved minor unclarity with the addition of noise peaks during data augmentation [#78](https://github.com/matchms/ms2deepscore/pull/78)

## [0.1.3] - 2021-03-09

## Changed

- Allow users to define L1 and L2 regularization of `SiameseModel` [#67](https://github.com/matchms/ms2deepscore/issues/67)
- Allow users to define number and size of `SiameseModel` [#64](https://github.com/matchms/ms2deepscore/pull/64)

## [0.1.2] - 2021-03-05

## Added

- `create_confusion_matrix_plot` in `plotting` [#58](https://github.com/matchms/ms2deepscore/pull/58)

## [0.1.1] - 2021-02-09

## Added

- noise peak addition during training via data generators [#55](https://github.com/matchms/ms2deepscore/pull/55)
- L1 and L2 regularization for first dense layer [#55](https://github.com/matchms/ms2deepscore/pull/55)

## Changed

- move vector calculation to separate calculate_vectors method [#52](https://github.com/matchms/ms2deepscore/pull/52)

## [0.1.0] - 2021-02-08

### Added

- This is the initial version of MS2DeepScore

[Unreleased]: https://github.com/matchms/ms2deepscore/compare/0.2.2...HEAD
[0.2.2]: https://github.com/matchms/ms2deepscore/compare/0.2.1...0.2.2
[0.2.1]: https://github.com/matchms/ms2deepscore/compare/0.2.0...0.2.1
[0.2.0]: https://github.com/matchms/ms2deepscore/compare/0.1.3...0.2.0
[0.1.3]: https://github.com/matchms/ms2deepscore/compare/0.1.2...0.1.3
[0.1.2]: https://github.com/matchms/ms2deepscore/compare/0.1.1...0.1.2
[0.1.1]: https://github.com/matchms/ms2deepscore/compare/0.1.0...0.1.1
[0.1.0]: https://github.com/matchms/ms2deepscore/releases/tag/0.1.0
![GitHub](https://img.shields.io/github/license/matchms/ms2deepscore)
[![PyPI](https://img.shields.io/pypi/v/ms2deepscore)](https://pypi.org/project/ms2deepscore/)
![GitHub Workflow Status](https://img.shields.io/github/workflow/status/matchms/ms2deepscore/CI%20Build)
[![SonarCloud Quality Gate](https://sonarcloud.io/api/project_badges/measure?project=matchms_ms2deepscore&metric=alert_status)](https://sonarcloud.io/dashboard?id=matchms_ms2deepscore)
[![SonarCloud Coverage](https://sonarcloud.io/api/project_badges/measure?project=matchms_ms2deepscore&metric=coverage)](https://sonarcloud.io/component_measures?id=matchms_ms2deepscore&metric=Coverage&view=list)  
[![DOI](https://zenodo.org/badge/310047938.svg)](https://zenodo.org/badge/latestdoi/310047938)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu)

# ms2deepscore
ms2deepscore provides a Siamese neural network that is trained to predict molecular structural similarities (Tanimoto scores) 
from pairs of mass spectrometry spectra. 

The library provides an intuitive classes to prepare data, train a siamese model,
and compute similarities between pairs of spectra.

In addition to the prediction of a structural similarity, 
MS2DeepScore can also make use of Monte-Carlo dropout to assess the model uncertainty.

## Reference
If you use MS2DeepScore for your research, please cite the following:

**"MS2DeepScore - a novel deep learning similarity measure to compare tandem mass spectra"**
Florian Huber, Sven van der Burg, Justin J.J. van der Hooft, Lars Ridder, 13, Article number: 84 (2021), Journal of Cheminformatics, doi: https://doi.org/10.1186/s13321-021-00558-4


## Setup
### Requirements

Python 3.7 or higher

### Installation
Simply install using pip: `pip install ms2deepscore`

### Prepare environment
We recommend to create an Anaconda environment with

```
conda create --name ms2deepscore python=3.8
conda activate ms2deepscore
pip install ms2deepscore
```
Alternatively, simply install in the environment of your choice by .


Or, to also include the full [matchms](https://github.com/matchms/matchms) functionality:
```
conda create --name ms2deepscore python=3.8
conda activate ms2deepscore
conda install --channel bioconda --channel conda-forge matchms
pip install ms2deepscore
```

## Getting started: How to prepare data, train a model, and compute similarities.
See [notebooks/MS2DeepScore_tutorial.ipynb](https://github.com/matchms/ms2deepscore/blob/main/notebooks/MS2DeepScore_tutorial.ipynb) 
for a more extensive fully-working example on test data.
If you are not familiar with `matchms` yet, then we also recommand our [tutorial on how to get started using matchms](https://blog.esciencecenter.nl/build-your-own-mass-spectrometry-analysis-pipeline-in-python-using-matchms-part-i-d96c718c68ee).

There are two different ways to use MS2DeepScore to compute spectral similarities. You can train a new model on a dataset of your choice. That, however, should preferentially contain a substantial amount of spectra to learn relevant features, say > 10,000 spectra of sufficiently diverse types.
The second way is much simpler: Use a model that was pretrained on a large dataset. 

## 1) Use a pretrained model to compute spectral similarities
We provide a model which was trained on > 100,000 MS/MS spectra from [GNPS](https://gnps.ucsd.edu/), which can simply be downloaded [from zenodo here](https://zenodo.org/record/4699356).
To then compute the similarities between spectra of your choice you can run something like:
```python
from matchms import calculate_scores()
from matchms.importing import load_from_msp
from ms2deepscore import MS2DeepScore
from ms2deepscore.models import load_model

# Import data
references = load_from_msp("my_reference_spectra.msp")
queries = load_from_msp("my_query_spectra.msp")

# Load pretrained model
model = load_model("MS2DeepScore_allGNPSpositive_10k_500_500_200.hdf5")

similarity_measure = MS2DeepScore(model)
# Calculate scores and get matchms.Scores object
scores = calculate_scores(references, queries, similarity_measure)
```

If you want to calculate all-vs-all spectral similarities, e.g. to build a network, than you can run:
```python
scores = calculate_scores(references, references, similarity_measure, is_symmetric=True)
```

To use Monte-Carlo Dropout to also get a uncertainty measure with each score, run the following:
```python
from matchms import calculate_scores()
from matchms.importing import load_from_msp
from ms2deepscore import MS2DeepScoreMonteCarlo
from ms2deepscore.models import load_model

# Import data
references = load_from_msp("my_reference_spectra.msp")
queries = load_from_msp("my_query_spectra.msp")

# Load pretrained model
model = load_model("MS2DeepScore_allGNPSpositive_10k_500_500_200.hdf5")

similarity_measure = MS2DeepScoreMonteCarlo(model, n_ensembles=10)
# Calculate scores and get matchms.Scores object
scores = calculate_scores(references, queries, similarity_measure)
```
In that scenario, `scores["score"]` contains the similarity scores (median of the ensemble of 10x10 scores) and `scores["uncertainty"]` give an uncertainty estimate (interquartile range of ensemble of 10x10 scores.

## 2) Train an own MS2DeepScore model
### Data preperation
Bin spectrums using `ms2deepscore.SpectrumBinner`. 
In this binned form we can feed spectra to the model.
```python
from ms2deepscore import SpectrumBinner
spectrum_binner = SpectrumBinner(1000, mz_min=10.0, mz_max=1000.0, peak_scaling=0.5)
binned_spectrums = spectrum_binner.fit_transform(spectrums)
```
Create a data generator that will generate batches of training examples.
Each training example consists of a pair of binned spectra and the corresponding reference similarity score.
```python
from ms2deepscore.data_generators import DataGeneratorAllSpectrums
dimension = len(spectrum_binner.known_bins)
data_generator = DataGeneratorAllSpectrums(binned_spectrums, tanimoto_scores_df,
                                           dim=dimension)
```
### Train a model
Initialize and train a SiameseModel. 
It consists of a dense 'base' network that produces an embedding for each of the 2 inputs.
The 'head' model computes the cosine similarity between the embeddings.
```python
from tensorflow import keras
from ms2deepscore.models import SiameseModel
model = SiameseModel(spectrum_binner, base_dims=(200, 200, 200), embedding_dim=200,
                     dropout_rate=0.2)
model.compile(loss='mse', optimizer=keras.optimizers.Adam(lr=0.001))
model.fit(data_generator,
          validation_data=data_generator,
          epochs=2)
```
### Predict similarity scores
Calculate similariteis for a pair of spectra
```python
from ms2deepscore import MS2DeepScore
similarity_measure = MS2DeepScore(model)
score = similarity_measure.pair(spectrums[0], spectrums[1])
```

## Contributing
We welcome contributions to the development of ms2deepscore! Have a look at the [contribution guidelines](https://github.com/matchms/ms2deepscore/blob/main/CONTRIBUTING.md).
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/). Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.rst).

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation);
1. you want to make a new release of the code base.

The sections below outline the steps in each case.

## You have a question

1. use the search functionality [here](https://github.com/matchms/ms2deepscore/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/matchms/ms2deepscore/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. make sure the existing tests still work by running ``python setup.py test``;
1. add your own tests (if necessary);
1. update or expand the documentation;
1. update the `CHANGELOG.md` file with change;
1. [push](http://rogerdudler.github.io/git-guide/>) your feature branch to (your fork of) the ms2deepscore repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.

## You want to make a new release of the code base

To create release you need write permission on the repository.

1. Check author list in `citation.cff` and `.zenodo.json` files
1. Bump the version using `bump2version <major|minor|patch>`. For example, `bump2version major` will increase major version numbers everywhere its needed (code, meta, etc.) in the repo.
1. Update the `CHANGELOG.md` to include changes made
1. Goto [GitHub release page](https://github.com/matchms/ms2deepscore/releases)
1. Press draft a new release button
1. Fill version, title and description field
1. Press the Publish Release button

A GitHub action will run which will publish the new version to [pypi](https://pypi.org/project/ms2deepscore).
Also a Zenodo entry will be made for the release with its own DOI.
