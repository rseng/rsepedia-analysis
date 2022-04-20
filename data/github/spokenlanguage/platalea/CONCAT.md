# Changelog
<!--
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
-->

## [Unreleased]
<!-- track upcoming changes here; move to new versioned section at release time -->
- added support for HowTo100M dataset

## [2.0] - 16 July 2021

Changes related to the ZeroSpeech challenge:
 - added support for SpokenCOCO dataset
 - added code to support the use of pretrained features + utility script to extract CPC features
 - refactored tokenization helpers making the tokenizer a global variable of dataset.py
 - changed platalea default config path ~/.platalea -> ~/.config/platalea
 - disabled use of wandb by default in basic.py and transformer.py experiments
 - pinning down pytorch version

Resolves issues #53, #103, #104 and (temporarily) solves #116.

## [1.0] - 9 December 2020

### Added
- Introducing an attention-based encoder-decoder architecture for speech recognition.
- Multitask training with multiple objectives (e.g. cross-modality retrieval and speech transcription) is also possible now.

<!--
### Removed

### Changed
-->

## [0.9] - 20 January 2020

State of the repo before @bhigy's merge leading to version 1.0.

[Unreleased]: https://github.com/spokenlanguage/platalea/compare/v1.0...HEAD
[2.0]: https://github.com/spokenlanguage/platalea/releases/tag/v2.0
[1.0]: https://github.com/spokenlanguage/platalea/releases/tag/v1.0
[0.9]: https://github.com/gchrupala/platalea/releases/tag/v0.9
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/).

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation);
1. you want to make a new release of the code base.

The sections below outline the steps in each case.

## You have a question

1. use the search functionality [here](https://github.com/spokenlanguage/platalea/issues) to see if someone already filed the same issue;
2. if your issue search did not yield any relevant results, make a new issue;
3. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/spokenlanguage/platalea/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. make sure the existing tests still work by running ``tox``;
1. also make sure you do not diverge from PEP8 standards by fixing any issues flake8 reports in the ``tox`` run; if you do diverge, make sure you clearly explain why;
1. add your own tests (if necessary);
1. update or expand the documentation;
1. push your feature branch to (your fork of) the platalea repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
# Platalea
Understanding visually grounded spoken language via multi-tasking

[![DOI](https://zenodo.org/badge/239750248.svg)](https://zenodo.org/badge/latestdoi/239750248)
[![install and run tests](https://github.com/egpbos/platalea/workflows/install%20and%20run%20tests/badge.svg?branch=master)](https://github.com/spokenlanguage/platalea/actions/workflows/pythonapp.yml)
[![codecov](https://codecov.io/gh/spokenlanguage/platalea/branch/master/graph/badge.svg)](https://codecov.io/gh/spokenlanguage/platalea)

## Installation

Clone this repo and cd into it:

```sh
git clone https://github.com/spokenlanguage/platalea.git
cd platalea
```

To install in a conda environment, assuming conda has already been installed, run the following to download and install dependencies:

```sh
conda create -n platalea python==3.8 pytorch -c conda-forge -c pytorch
conda activate platalea
pip install torchvision
```

Then install platalea with:

```sh
pip install .
```

### Experiment dependencies
Different experiments may have different additional dependencies.
The `basic` experiment needs the following:

```sh
pip install sklearn python-Levenshtein
```

## Datasets

### Flickr8K
The repository has been developed to work with Flickr8K dataset. The code can
be made to work with other datasets but this will require some adaptations.

To use Flickr8K, you need to download:
* [Flickr8K](http://hockenmaier.cs.illinois.edu/Framing_Image_Description/KCCA.html) [1].
  Note that downloading from the official website seems broken at the moment.
  Alternatively, the dataset can be obtained from
  [here](https://github.com/jbrownlee/Datasets/blob/master/Flickr8k_Dataset.names).
* The [Flickr Audio Caption Corpus](https://groups.csail.mit.edu/sls/downloads/flickraudio/) [2].
* Some additional [metadata files](https://surfdrive.surf.nl/files/index.php/s/EF1bA9YYfhiBxoN).

Create a folder to store the dataset (we will assume here that the folder is
`~/corpora/flickr8k`)  and move all the files you downloaded there, then
extract the content of the archives.  You can now setup the environment and
start preprocessing the data.

#### Configuration

We use ConfigArgParse for setting necessary input variables, including the
location of the dataset.  This means you can use either a configuration file
(config.ini or config.yml), environment variables or command line arguments to
specify the necessary configuration parameters.

To specify the location of the dataset, one option is to create a configuration
file under your home directory (`~/.config/platalea/config.yml`)., with
follwing content:

```
flickr8k_root   /home/<user>/corpora/flickr8k
```

The same result can be achieved with an environment variable:

```sh
export FLICKR8K_ROOT=/home/<user>/corpora/flickr8k
```

You could also specify this option directly on the command line when running
an experiment (the respective options would be `--flickr8k_root=...`).

#### Preprocessing

Run the preprocessing script to extract input features:

```bash
python platalea/utils/preprocessing.py flickr8k
```

### Howto100Men-cc

This repository has support for a subset of the [howto100M dataset](https://github.com/antoine77340/howto100m). The subset contains all videos with creative commons license that claim to be in english language according to their metadata.

- Sampling rate of the audio features is 100Hz
- Sampling rate of the video features is 1Hz.
- A dataset item is defined as the combination of video and audio for a fragment of 3 seconds.

#### Preprocessing

This code contains functionality for extracting audio features from the videos. These files need to be in the datafolder before preprocessing the dataset. Preprocessing will create an index file with references to the video feature files. The video features (S3D) need to be acquired elsewhere. Videos from howto100m need to be downloaded from youtube using the metadata from the howto100m dataset provided above.

To start preprocessing, run the following:

```sh
python -m platalea.utils.preprocessing howto100m-encc --howto100m_root /corpora/howto100m/
```

#### Running/Training

Running any experiment using the howto100M dataset has not yet been implemented. 


## Training

You can now train a model using one of the examples provided under
`platalea/experiments`, e.g.:

```sh
cd platalea/experiments/flickr8k
mkdir -p runs/test
cd runs/test
python -m platalea.experiments.flickr8k.basic
```

After the model is trained, results are available in `results.json`.

### Weights and Biases (wandb)

Some experiments support the use of wandb for cloud logging of results.
In the examples we provide under `platalea/experiments`, this option is disabled by default.
To force-enable it, the call to `experiment()` should be changed from `experiment(..., wandb_mode='disabled')` to `experiment(..., wandb_mode='online')`. To default back to wandb normal behavior (where the mode can be set through command line or environment variable), use `wandb_mode=None`.

## Contributing

If you want to contribute to the development of platalea, have a look at the [contribution guidelines](CONTRIBUTING.md).

## Changelog

We keep track of what is added, changed and removed in releases in the [changelog](CHANGELOG.md).

## References

[1] Hodosh, M., Young, P., & Hockenmaier, J. (2013). Framing Image Description
as a Ranking Task: Data, Models and Evaluation Metrics. Journal of Artificial
Intelligence Research, 47, 853–899. https://doi.org/10.1613/jair.3994.

[2] Harwath, D., & Glass, J. (2015). Deep multimodal semantic embeddings for
speech and images. 2015 IEEE Workshop on Automatic Speech Recognition and
Understanding (ASRU), 237–244. https://doi.org/10.1109/ASRU.2015.7404800.
<!-- Describe your PR here -->

<!-- Please, make sure the following items are checked -->
Checklist before merging:

- [ ] Does this PR warrant a version bump? If so, make sure you add a git tag and GitHub release.
- [ ] Did you update the changelog with a user-readable summary? Please, include references to relevant issues or PR discussions. See [CII [release_notes] criterion](https://bestpractices.coreinfrastructure.org/en/criteria/0#0.release_notes_vulns).
- [ ] When adding new functionality: did you also add a test for it?
