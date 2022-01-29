# AUDIOMATE

[![PyPI](https://img.shields.io/pypi/v/audiomate.svg)](https://pypi.python.org/pypi/audiomate)
[![Build Status](https://travis-ci.com/ynop/audiomate.svg?branch=master)](https://travis-ci.com/ynop/audiomate)
[![Documentation Status](https://readthedocs.org/projects/audiomate/badge/?version=latest)](https://audiomate.readthedocs.io/en/latest/?badge=latest)
[![DeepSource](https://static.deepsource.io/deepsource-badge-light-mini.svg)](https://deepsource.io/gh/ynop/audiomate/?ref=repository-badge)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02135/status.svg)](https://doi.org/10.21105/joss.02135)

Audiomate is a library for easy access to audio datasets.
It provides the datastructures for accessing/loading different datasets in a generic way.
This should ease the use of audio datasets for example for machine learning tasks.

```python
import audiomate
from audiomate.corpus import io

# Download a dataset
esc_downloader = io.ESC50Downloader()
esc_downloader.download('/local/path')

# Load and work with the dataset
esc50 = audiomate.Corpus.load('/local/path', reader='esc-50')

# e.g. Read the audio signal and the label of specific sample/utterance
utterance = esc50.utterances['1-100032-A-0']
samples = utterance.read_samples()
label_list = utterance.label_lists[audiomate.corpus.LL_SOUND_CLASS]

for label in label_list:
  print(label.start, label.value)
```

Furthermore it provides tools for interacting with datasets
(validation, splitting, subsets, merge, filter), extracting features,
feeding samples for training ML models and more.

* [Documentation](https://audiomate.readthedocs.io)
* [Examples](https://github.com/ynop/audiomate/tree/master/examples)
* [Changelog](https://audiomate.readthedocs.io/en/latest/notes/changelog.html)

Currently supported datasets:
* [Acoustic Event Dataset](https://arxiv.org/pdf/1604.07160.pdf)
* [AudioMNIST](https://github.com/soerenab/AudioMNIST)
* [Mozilla Common Voice](https://voice.mozilla.org/)
* [ESC-50](https://github.com/karoldvl/ESC-50)
* [Fluent Speech Commands](http://www.fluent.ai/research/fluent-speech-commands/)
* [Free Spoken Digit Dataset](https://github.com/Jakobovski/free-spoken-digit-dataset)
* [German Distant Speech Corpus](https://www.inf.uni-hamburg.de/en/inst/ab/lt/resources/data/acoustic-models.html)
* [Google Speech Commands](https://research.googleblog.com/2017/08/launching-speech-commands-dataset.html)
* [GTZAN](http://marsyas.info/downloads/datasets.html)
* [LibriSpeech](https://www.openslr.org/12/)
* [M-AILABS Speech Dataset](https://www.caito.de/2019/01/the-m-ailabs-speech-dataset/)
* [MUSAN](http://www.openslr.org/17/)
* [LITIS Rouen Audio scene dataset](https://sites.google.com/site/alainrakotomamonjy/home/audio-scene)
* [Spoken Wikipedia Corpora](https://nats.gitlab.io/swc/)
* [Tatoeba](https://tatoeba.org/)
* [TIMIT](https://github.com/philipperemy/timit)
* [Urbansound8k](http://urbansounddataset.weebly.com/urbansound8k.html)
* [Voxforge](http://www.voxforge.org/de)

Currently supported formats:
* [Kaldi](http://kaldi-asr.org/)
* [Mozilla DeepSpeech](https://github.com/mozilla/DeepSpeech)
* [Wav2Letter](https://github.com/facebookresearch/wav2letter)
* [NVIDIA Jasper](https://github.com/NVIDIA/DeepLearningExamples/tree/master/PyTorch/SpeechRecognition/Jasper)
* [Custom Formats](https://audiomate.readthedocs.io/en/latest/documentation/formats.html)

## Installation

```sh
pip install audiomate
```

Install the latest development version:

```sh
pip install git+https://github.com/ynop/audiomate.git
```

### Dependencies

#### sox
For parts of the functionality (e.g. audio format conversion) [sox](http://sox.sourceforge.net) is used. In order to use it, you have to install sox.

```sh
# macos
brew install sox

# with support for specific formats
brew install sox --with-lame --with-flac --with-libvorbis

# linux
apt-get install sox

# anaconda for macOS/windows/linux:
conda install -c conda-forge sox
```

## Development

### Prerequisites

* [A supported version of Python > 3.5](https://docs.python.org/devguide/index.html#status-of-python-branches)

It's recommended to use a virtual environment when developing audiomate.
To create one, execute the following command in the project's root directory:

```
python -m venv .
```

To install audiomate and all it's dependencies, execute:

```
pip install -e .
```

### Running the test suite

```
pip install -e .[dev]
pytest
```

With PyCharm you might have to change the default test runner. Otherwise, it might only suggest to use nose. To do so, go to File > Settings > Tools > Python Integrated Tools (on the Mac it's PyCharm > Preferences > Settings > Tools > Python Integrated Tools) and change the test runner to py.test.

### Benchmarks

In order to check the runtime of specific parts, ``pytest-benchmark`` is used. Benchmarks are normal test functions, but call the benchmark fixture for the code under test.

To run benchmarks:

```
# Run all
pytest bench

# Specific benchmark
pytest bench/corpus/test_merge_corpus.py
```

To compare between different runs:

```
pytest-benchmark compare
```

### Editing the Documentation

The documentation is written in [reStructuredText](http://docutils.sourceforge.net/rst.html) and transformed into various output formats with the help of [Sphinx](http://www.sphinx-doc.org/).

* [Syntax reference reStructuredText](http://docutils.sourceforge.net/docs/user/rst/quickref.html)
* [Sphinx-specific additions to reStructuredText](http://www.sphinx-doc.org/en/stable/markup/index.html)

To generate the documentation, execute:

```
pip install -e .[dev]
cd docs
make html
```

The generated files are written to `docs/_build/html`.

### Versions

Versions is handled using [bump2version](https://github.com/c4urself/bump2version). To bump the version:

```
bump2version [major,minor,patch,release,num]
```

In order to directly go to a final relase version (skip .dev/.rc/...):

```
bump2version [major,minor,patch] --new-version x.x.x
```

### Release

Commands to create a new release on pypi.

```
rm -rf build
rm -rf dist

python setup.py sdist
python setup.py bdist_wheel
twine upload dist/*
```
# Contributing to audiomate

If you would like to contribute code or documentation to audiomate you can do so through GitHub by forking the repository and starting a pull request. Please follow the guidelines in this document when preparing the pull request. If you need help compiling the code, head to the [README](README.md) for quick instructions. All the details are outlined in the [separate documentation](docs).

## Code Conventions and Housekeeping

Please make every effort to follow the existing conventions and style in order to keep the code as readable as possible. This makes it easier to review your contribution and reduces the amount of work necessary before a merge.

* When writing a commit message please follow [these conventions](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html). If you are addressing an existing issue, add `Fixes GH-XXXX` at the end of the commit message (where `XXXX` denotes the issue number).
* Run `flake8` to ensure that the formatting of your code matches the project's code style. Audiomate's code is formatted according to [PEP 8 -- Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/), but the maximum line length is 100 characters and strings have to be delimited by single quotes.
* Add docstrings to the public API in accordance with [PEP 257 -- Docstring Conventions](https://www.python.org/dev/peps/pep-0257/).
* Cover your changes with unit tests.
* Before opening the pull request, rebase it onto the HEAD of the current master (or the relevant target branch).
---
title: 'audiomate: A Python package for working with audio datasets'
tags:
    - Python
    - audio
    - speech
    - music
    - corpus
    - dataset
authors:
    - name: Matthias Büchi
      orcid: 0000-0003-0207-5711
      affiliation: 1
    - name: Andreas Ahlenstorf
      affiliation: 1
affiliations:
    - name: ZHAW Zurich University of Applied Sciences, Winterthur, Switzerland
      index: 1
date: 20 December 2019
bibliography: paper.bib
---

# Summary

Machine learning tasks in the audio domain frequently require large datasets with training data.
Over the last years, numerous datasets have been made available for various purposes, for example, [@musan2015] and [@ardila2019common].
Unfortunately, most of the datasets are stored in widely differing formats.
As a consequence, machine learning practitioners have to convert datasets into other formats before they can be used or combined.
Furthermore, common tasks like reading, partitioning, or shuffling of datasets have to be developed over and over again for each format and require intimate knowledge of the formats.
We purpose Audiomate, a Python toolkit, to solve this problem.

Audiomate provides a uniform programming interface to work with numerous datasets.
Knowledge about the structure or on-disk format of the datasets is not necessary.
Audiomate facilitates and simplifies a wide range of tasks:

* Reading and writing of numerous dataset formats using a uniform programming interface, for example [@musan2015], [@Panayotov2015LibrispeechAA] and [@ardila2019common]
* Accessing metadata, like speaker information and labels
* Reading audio data (single files, batches of files)
* Retrieval of information about the data (e.g., number of speakers, total duration).
* Merging of multiple datasets (e.g., combine two speech datasets).
* Splitting data into smaller subsets (e.g., create training, validation, and test sets with a reasonable distribution of classes).
* Validation of data for specific requirements (e.g., check whether all samples were assigned a label)

# Use Cases

To illustrate Audiomate's capabilities, we present two typical applications where Audiomate significantly simplifies the task of a developer: Training a speech recognition model with Mozilla's implementation of DeepSpeech and training a deep neural network to recognize music.

## Converting Datasets

In this example, we illustrate how to employ Audiomate to convert the LibriSpeech dataset [@Panayotov2015LibrispeechAA] into the CSV-format expected by Mozilla's implementation (https://github.com/mozilla/DeepSpeech) of DeepSpeech [@deepspeech] which can, in turn, be used to train an automatic speech recognition model.

```python
import audiomate
from audiomate.corpus import io

# Download LibriSpeech corpus
downloader = io.LibriSpeechDownloader()
downloader.download('/local/data/librispeech')

# Read LibriSpeech
reader = io.LibriSpeechReader()
librispeech = reader.load('/local/data/librispeech')

# Save in DeepSpeech format
writer = io.MozillaDeepSpeechWriter()
writer.save(librispeech, '/local/data/librispeech_ds')
```

Knowledge of the on-disk formats of the datasets is not required.

Some datasets contain invalid or corrupted files.
If those are known, Audiomate tries to rectify the problems or automatically excludes those files before processing any data.

## Merging and Partitioning Datasets

Another area where Audiomate excels is mixing datasets and partitioning them into training, test, and validation sets.
Assume that the task is to train a neural network to detect segments in audio streams that are music.
MUSAN [@musan2015] and GTZAN [@GTZAN] are two suitable datasets for this task because they provide a wide selection of music, speech, and noise samples.
In the example below, we first download MUSAN and GTZAN to the local disk before creating `Loader` instances for each format that allow Audiomate to access both datasets using a unified interface. Then, we instruct Audiomate to merge both datasets.
Afterwards, we use a `Splitter` to partition the merged dataset into a train and test set.
By merely creating views, Audiomate avoids creating unnecessary disk I/O and is therefore ideally suited to work with large datasets in the range of tens or hundreds of gigabytes.
Ultimately, we load the samples and labels by iterating over all utterances.
Audio samples are numpy arrays. They allow for fast access, high processing speed and ensure interoperability with
third-party programs that can operate on numpy arrays, for example TensorFlow or PyTorch.
Alternatively, it is possible to load the samples in batches, which is ideal for feeding them to a deep learning toolkit like PyTorch.

```python
import audiomate
from audiomate.corpus import io
from audiomate.corpus import subset

musan_dl = io.MusanDownloader()
musan_dl.download('/local/data/musan')

gtzan_dl = io.GtzanDownloader()
gtzan_dl.download('/local/data/gtzan')

musan = audiomate.Corpus.load('/local/data/musan', reader='musan')
gtzan = audiomate.Corpus.load('/local/data/gtzan', reader='gtzan')

full = audiomate.Corpus.merge_corpora([musan, gtzan])

splitter = subset.Splitter(full, random_seed=222)
subviews = splitter.split(proportions={
    'train': 0.8,
    'test': 0.2,
})

for utterance in subviews['train'].utterances.values():
    samples = utterance.read_samples()
    labels = utterance.label_lists[audiomate.corpus.LL_DOMAIN]
```

# Implementation

Audiomate was designed with extensibility in mind.
Therefore, it is straightforward to add support for additional data formats.
Support for another format can be added by implementing at least one of three available abstract interfaces.

* Reader: A `Reader` defines the procedure to load data that is structured in a specific format.
          It converts it into a Audiomate-specific data structure.

* Writer: A `Writer` defines the procedure to store data in a specific format.
          It does that by converting the data from the Audiomate-specific data structure into the target format.

* Downloader: A `Downloader` can be used to download a dataset.
              It downloads all required files automatically.

Rarely all interfaces are implemented for a particular format.
Usually, `Reader` and `Downloader` are implemented for datasets, while `Writer` is implemented for machine learning toolkits.

Audiomate supports more than a dozen datasets and half as many toolkits.

# Related Work

A variety of frameworks and tools offer functionality similar to Audiomate.

**Data loaders** Data loaders are libraries that focus on downloading and preprocessing data sets to make them easily accessible without requiring a specific tool or framework.
In contrast to Audiomate, they cannot convert between formats, split or merge data sets.
Examples of libraries in that category are [@mirdata], [@speechcorpusdownloader], and [@audiodatasets].
Furthermore, some of these libraries focus on a particular kind of data, such as music, and do not assist with speech data sets.

**Tools for specific frameworks** Various machine learning tools and deep learning frameworks include the necessary infrastructure to make various datasets readily available to their users.
One notable example is TensorFlow [@tensorflow], which includes data loaders for different kinds of data, including image, speech, and music data sets, such as [@ardila2019common].
Another one is torchaudio [@torchaudio] for PyTorch, which not only offers data loaders but is also capable of converting between various formats.
In contrast to Audiomate, those tools or libraries support a specific machine learning or deep learning framework (TensorFlow or PyTorch, respectively), whereas Audiomate is framework agnostic.

# References
