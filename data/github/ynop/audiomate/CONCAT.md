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
Welcome to audiomate's documentation!
=====================================

Audiomate is a library for easy access to audio datasets. It provides the datastructures for accessing/loading different datasets in a generic way.
This should ease the use of audio datasets for example for machine learning tasks.

.. image:: documentation/idea.*


.. toctree::
    :maxdepth: 1
    :caption: Notes

    notes/installation
    notes/changelog

.. toctree::
    :maxdepth: 2
    :caption: Documentation

    documentation/structure
    documentation/formats
    documentation/new_dataset_format
    documentation/data_mapping
    documentation/indirect_support
    documentation/logging.rst

.. toctree::
    :maxdepth: 2
    :caption: Package Reference

    reference/tracks
    reference/annotations
    reference/issuers
    reference/containers
    reference/corpus
    reference/io
    reference/subset
    reference/validation
    reference/conversion
    reference/processing
    reference/encoding
    reference/feeding
    reference/formats
    reference/utils

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
.. _section_broadcast_format:

Broadcast Format
================

The broadcast format is basically the same as :ref:`section_default_format`, except it uses another format to store labels.
This format is meant for data where not many utterances are given, but with a lot of labels. So instead to have all labels per label-list in one file,
a label-file per utterance is used.

**labels.txt**

This files defines where to find the effective label files. It stores the label-file path per utterance. Additionaly a label-list-id can be given, if there are multiple label-lists per utterance.

.. code-block:: bash

    <utt-id> <label-file-path> <label-list-idx>

Example:

.. code-block:: bash

    utt-1 files/a/labels.txt
    utt-2 files/b/music.txt music
    utt-2 files/b/jingles.txt jingles
    utt-3 files/c/trailers.txt

**[label-file].txt**

The label files reference by the *labels.txt* are in the following format. It contains the start and end in seconds.
The values are **Tab-separated**.
Optionally additional meta-information can be stored per label.
This has to be a json string in square brackets with a space separated after the label-value.


.. code-block:: bash

    <start> <end>   <value> [<label-meta>]

Example:

.. code-block:: bash

    0	40  hallo
    40.5    100 velo
    102.4   109.2   auto [{"lang": "de", "type": 2}]
Add Dataset/Format
==================

In this section it is described how to a downloader, reader or writer for a new dataset or another corpus format.
The implementation is pretty straight-forward. For examples checkout some of the existing implementations at
:mod:`audiomate.corpus.io`.

.. IMPORTANT::

    * Use the same name (``type()`` method) for downloader/reader/writer.
    * Import your components in :mod:`audiomate.corpus.io.__init___`. So all components are available from the io module.
    * Checkout :ref:`data-mapping` on what and how to add info/data when reading a corpus.

Corpus Downloader
-----------------
If we aim to load some specific dataset/corpus, a downloader can be implemented,
if it is possible to automate the whole download process. First we create a new class that inherits from
:class:`audiomate.corpus.io.CorpusDownloader`. There we have to implement two methods.
The ``type`` method just has to return a string with the name of the new dataset/format.
The ``_download`` method will do the heavy work of download all the files to the path ``target_path``.

.. code-block:: python

    from audiomate.corpus.io import base


    class MyDownloader(base.CorpusDownloader):

        @classmethod
        def type(cls):
            return 'MyDataset'

        def _download(self, target_path):
            # Download the data to target_path

In the module :mod:`audiomate.corpus.io.downloader`, common base classes for downloaders are implemented. This is useful since for a lot of corpora the way of downloading is similar.

* :class:`audiomate.corpus.io.ArchiveDownloader`: For corpora based on a single archive.

Corpus Reader
-------------

The reader is the one component that is mostly used. Either for a specific dataset/corpus or a custom format,
a reader is most likely to be required. First we create a new class that inherits from
:class:`audiomate.corpus.io.CorpusReader`. There we have to implement three methods.
The ``type`` method just has to return a string with the name of the new dataset/format.
The ``_check_for_missing_files`` method can be used to check if the given path is a valid input.
For example if the format/dataset requires some specific meta-files it can be check here if they are available.
Finally in the ``_load`` method the actual loading is done and the loaded corpus is returned.

.. code-block:: python

    from audiomate.corpus.io import base


    class MyReader(base.CorpusReader):

        @classmethod
        def type(cls):
            return 'MyDataset'

        def _check_for_missing_files(self, path):
            # Check the path for missing files that are required to read with this reader.
            # Return a list of missing files
            return []

        def _load(self, path):
            # Create a new corpus
            corpus = audiomate.Corpus(path=path)

            # Create files ...
            corpus.new_file(file_path, file_idx)

            # Issuers ...
            issuer = assets.Speaker(issuer_idx)
            corpus.import_issuers(issuer)

            # Utterances with labels ...
            utterance = corpus.new_utterance(file_idx, file_idx, issuer_idx)
            utterance.set_label_list(annotations.LabelList(idx='transcription', labels=[
                annotations.Label(str(digit))
            ]))

            return corpus

For some datasets there are files/utterances that are not valid.
(This can be due to a corrupt file, invalid transcription, ...)
For this case a json-file ``audiomate/corpus/io/data/[reader-type]/invalid_utterances.json`` can be created
that contains a list with ids of invalid utterances.
The ids correspond to id of the utterance, if it would be loaded anyway.

Testing
^^^^^^^
For testing a reader the :class:`tests.corpus.io.reader_test.CorpusReaderTest` can be used.
It provides base test methods for checking the correctness/existence of the basic components (tracks, utterances, labels, ...).

.. code-block:: python

   from tests.corpus.io import reader_test as rt

   class TestMyReader(rt.CorpusReaderTest):

      #
      # Define via EXPECTED_* variables, what components are expected to be loaded
      #
      EXPECTED_NUMBER_OF_TRACKS = 3
      EXPECTED_TRACKS = [
         rt.ExpFileTrack('file-id', '/path/to/file'),
      ]

      #
      # Override the load method, that loads the sample-corpus.
      #
      def load(self):
         return MyReader().load('/path/to/sample/corpus')


For testing any custom functionality specific test-methods can be added as well.

Corpus Writer
-------------

A writer is only useful for custom formats. For a specific dataset a writer is most likely not needed.
First we create a new class that inherits from :class:`audiomate.corpus.io.CorpusWriter`.
There we have to implement two methods.
The ``type`` method just has to return a string with the name of the new dataset/format.
The ``_save`` method does the serialization of the given corpus to the given path.


.. code-block:: python

    from audiomate.corpus.io import base


    class DefaultWriter(base.CorpusWriter):

        @classmethod
        def type(cls):
            return 'MyDataset'

        def _save(self, corpus, path):
            # Do the serialization
Indirectly Supported Corpora
============================

Some corpora are hard to integrate directly, e.g. due to necessary preprocessing steps.
Corpus Formats
==============

.. toctree::
    :hidden:

    default_format
    broadcast_format


A corpus format defines how a corpus is saved on disk. For the use with this library some formats were specifically developed:

* Default :class:`audiomate.corpus.io.DefaultReader` / :class:`audiomate.corpus.io.DefaultWriter`
* Broadcast :class:`audiomate.corpus.io.BroadcastReader`

Furthermore there exist downloaders, readers and writers for other formats or specific datasets.
For a list of available downloaders, readers and writers check :ref:`io_implementations`.


.. _section_default_format:

Default Format
==============

This describes, how a corpus with the default format is saved on disk. Every corpus is a folder with a bunch of files.

**files.txt**

This file contains a list of every audio file in the corpus. Every file is identified by a unique id.
Every line in the file contains the mapping from file-id to the file-path for a single file. The filepath is the path to the audio file relative to the corpus folder.

.. code-block:: bash

    <recording-id> <wav-file-path>

Example:

.. code-block:: bash

    2014-03-17-09-45-16_Kinect-Beam train/2014-03-17-09-45-16_Kinect-Beam.wav
    2014-03-17-09-45-16_Realtek train/2014-03-17-09-45-16_Realtek.wav
    2014-03-17-09-45-16_Yamaha train/2014-03-17-09-45-16_Yamaha.wav
    2014-03-17-10-26-07_Realtek train/2014-03-17-10-26-07_Realtek.wav


**utterances.txt**

This file contains all utterances in the corpus. An utterance is a part of a file (A file can contain one or more utterances).
Every line in this file defines a single utterance, which consists of utterance-id, file-id, start and end. Start and end are measured in seconds within the file.
If end is -1 it is considered to be the end of the file (If the utterance is the full length of the file, start and end are 0/-1).

.. code-block:: bash

    <utterance-id> <recording-id> <start> <end>

Example:

.. code-block:: bash

    1_hello 2014-03-17-09-45-16_Kinect-Beam
    1_hello_sam 2014-03-17-09-45-16_Realtek 0 -1
    2_this_is 2014-03-17-09-45-16_Yamaha 0 5
    3_goto 2014-03-17-09-45-16_Yamaha 5 -1

**utt_issuers.txt**

This file contains the mapping from utterance to issuers, which gives the information who/what is the origin of a given utterance (e.g. the speaker).
Every line contains one mapping from utterance-id to issuer-id.

.. code-block:: bash

    <utterance-id> <issuer-id>

Example:

.. code-block:: bash

    1_hello marc
    1_hello_sam marc
    2_this_is sam
    3_goto jenny


**issuers.json**

This file contains additional information about the issuers.
Depending on the type of the issuer (Issuer, Speaker, Artist which are defined via `type` parameter) different parameters can be set.

Issuer:
    * info: An arbitrary info dictionary.

Speaker:
    * gender (MALE/FEMALE): The gender of a speaker
    * age_group (child, youth, adult, senior): The age group of a speaker
    * native_language (language code ISO 639-3): The language natively spoken

Artist:
    * name: The name of the artist/band/group

.. code-block:: json

    {
      "speaker-1": {
        "info": {},
        "type": "speaker",
        "gender": "MALE",
        "age_group": "ADULT",
        "native_language": "deu"
      },
      "speaker-2": {
        "info": {},
        "type": "artist",
        "name": "Ohooo"
      },
      "speaker-3": {
        "info": {
          "region": "zh"
        }
      }
    }


**labels_[x].txt**

There can be multiple label-lists in a corpus (e.g. text-transcription, raw-text-transcription - with punctuation, audio classification type, ...).
Every label-list is saved in a separate file with the prefix *labels_*.
A single file contains labels of a specific type for all utterances. A label-list of an utterance can contain one or more labels (e.g. in a text segmentation every word could be a label).
A label optionally can have a start and end time (in seconds within the utterance). For labels without start/end defined 0/-1 is set.
Every line in the file defines one label. The labels are stored in order per utterance (e.g. 1. word, 2. word, 3. word, ...).
Optionally addtional meta-information can be stored per label. This has to be a json string in square brackets.

.. code-block:: bash

    <utterance-id> <start> <end> <label-value> [<label-meta>]

Example:

.. code-block:: bash

    1_hello 0 -1 hi
    1_hello 0 -1 this
    1_hello 0 -1 is
    1_hello_sam 0 -1 hello
    1_hello_sam 0 -1 sam
    2_this_is 0 -1 this
    2_this_is 0 -1 is [{"prio": 3}]
    2_this_is 0 -1 me [{"stress": true}]
    3_goto 0 -1 go
    3_goto 0 -1 to
    3_goto 0 -1 the
    3_goto 0 -1 mall

**features.txt**

Contains a list of stored features. A corpus can have different feature containers. Every container contains the features of all utterances of a given type (e.g. MFCC features).
A feature container is a h5py file which contains a dataset per utterance. Every line contains one container of features.

.. code-block:: bash

    <feature-name> <relative-path>

Example:

.. code-block:: bash

    mfcc mfcc_features
    fbank fbank_features

**audio.txt**

Contains a list of tracks that are stored in audio-containers.
Every entry consists of a track-id, the relative path to the container and
a key that identifies the track in the audio-container.

.. code-block:: bash

    <track-id> <audio-container-path> <audio-container-key>

Example:

.. code-block:: bash

    track-1 ../audio.hdf5 track-1
    track-2 ../audio.hdf5 track-2x
Corpus Structure
================

To represent any corpus/dataset in a generic way, a structure
is needed that can represent the data of any audio dataset as good as possible.
The basic structure consists of the following components.

.. image:: basic_structure.*

Corpus
------
The Corpus is the main object that represents a dataset/corpus.

Track
-----
A track is an abstract representation of an audio signal.
There are currently two implementations.
One that reads the audio signal from a file
and one that read the audio signal from a HDF5 container.

Utterance
---------
An utterance represents a segment of a track.
It is used to divide a track into independent segments.
A track can have one or more utterances.
The utterances are basically the samples in terms of machine learning.

Issuer
------
The issuer is defined as the person/thing/... who generate/produced the utterance (e.g. The speaker who read a given utterance).

An issuer can be further distinguished into different types.
The current implementation provides classes for speaker (for spoken audio content)
and for artists (for musical content).

LabelList
---------
The label-list is a container for holding all labels of a given type for one utterance.
For example there is a label-list containing the textual transcription of recorded speech.
Another possible type of label-list could hold all labels classifying the audio type (music, speech, noise) of every part of a radio broadcast recording.

Label
-----
The label is defining any kind of annotation for a part of or the whole utterance.

FeatureContainer
----------------
A feature-container is a container holding the feature matrices of a given type (e.g. mfcc) for all utterances.
A corpus can contain multiple feature-containers.
.. _data-mapping:

Data Mapping
============

Since we want to have a consistent abstraction of different formats and datasets,
it is important that all data and information is mapped correctly into the python classes.

Issuer
------

The issuer holds information about the source of the audio content.
Depending on the audio content different attributes are important.
Therefore different types of issuers can be used.

Speech
    For audio content that mainly contains spoken content the :class:`audiomate.issuers.Speaker` has to be used.
    This is most common for datasets regarding speech recognition/synthesis etc.

Music
    For audio content that contains music, the :class:`audiomate.issuers.Artist` has to be used.

Labels
------

In the corpus data structures an utterance can have multiple label-lists. In order to access a label-list a key is used.

.. code-block:: python

    utterance = ...
    label_list = utterance.label_lists['word-transcription']

The used key should be consistent for all datasets. Therefore the identifiers/keys should be selected from below
if possible. For these predefined keys, constants are defined in :mod:`audiomate.corpus`.

general
^^^^^^^

domain
    A high-level category for a given audio excerpt. Should be one of the following values:

        * speech
        * music
        * noise


speech
^^^^^^

word-transcript
    Non-aligned transcription of speech.

word-transcript-raw
    Non-aligned transcription of speech. Used for unprocessed transcriptions (e.g. containing punctuation, ...).

word-transcript-aligned
    Aligned transcription of speech. The begin and end of the words is defined.
    Every word is a single label in the label-list.

phone-transcript
    Non-aligned transcription of phones.

phone-transcript-aligned
    Aligned transcription of phones. Begin and end of phones is defined.

music
^^^^^

genre
    The genre of the music.

noise
^^^^^

sound-class
    Labels defining any sound-event, acoustic-scene, environmental noise, ...
    e.g. siren, dog_bark, train, car, snoring ...


This list isn't complete. Please open an issue for any additional domains/classes that maybe needed.
.. _logging:

Logging
=======

Logging in audiomate is done using the standard Python logging facilities.

Enable Logging
--------------

By default, only messages of severity ``Warning`` or higher are printed to ``sys.stderr``.
Audiomate provides detailed information about progress of long-running tasks with messages of severity ``Info``.
To enable logging of messages of lower severity, configure Python's logging system as follows:

.. code-block:: python

    import logging

    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)-15s  %(name)s  %(message)s'
    )

For further information check the python `logging documentation <https://docs.python.org/3/howto/logging.html>`_.

Create log messages in audiomate
--------------------------------

Logging in audiomate is done with a single logger.
The logger is available in :mod:`audiomate.logutil`.

.. code-block:: python

    from audiomate import logutil

    logger = logutil.getLogger()

    def some_functionality():
        logger.debug('message')

Since audiomate has a lot of long-running tasks,
a special function for logging the progress of a loop can be used.
It basically is a wrapper around an iterable to check and log the progress.
In order to keep the logs as small as possible,
progress is logged in steps of 5 minutes.


.. code-block:: python

    from audiomate import logutil

    logger = logutil.getLogger()

    for utterance in logger.progress(
            corpus.utterances.values(),
            total=corpus.num_utterances,
            description='Process utterances'):

        # Do something with the utterance,
        # that takes up some time.
Installation
============

Install the latest stable version::

    pip install audiomate

Install the latest development version::

    pip install git+https://github.com/ynop/audiomate.git

Dependencies
------------

**sox**

For parts of the functionality (e.g. audio format conversion) `sox <http://sox.sourceforge.net>`_ is used. In order to use it, you have to install sox.

.. code-block:: bash

   # macos
   brew install sox

   # with support for specific formats
   brew install sox --with-lame --with-flac --with-libvorbis

   # linux
   apt-get install sox

   # anaconda for macOS/windows/linux:
   conda install -c conda-forge soxChangelog
=========

Next Version
------------

v6.0.0
------

**Breaking Changes**

* Drop support of Python 3.5 because a required dependency (llvmlite) does not support it anymore.

**New Features**

* Setup consistent way for logging. (:ref:`logging`)

* Added downloader (:class:`audiomate.corpus.io.CommonVoiceDownloader`) for the `Common Voice Corpora <https://voice.mozilla.org/de/datasets>`_.

* Add existence checks for reader (:class:`audiomate.corpus.io.CorpusReader`) to see if folder exists.

* Add existence checks and a option for forcing redownload for downloader (:class:`audiomate.corpus.io.CorpusDownloader`).

v5.2.0
------

**New Features**

* Added reader (:class:`audiomate.corpus.io.LibriSpeechReader`) and
  downloader (:class:`audiomate.corpus.io.LibriSpeechDownloader`) for the
  `LibriSpeech Dataset <https://www.openslr.org/12/>`_.

v5.1.0
------

**New Features**

* Added Downloader for SWC Corpus ((:class:`audiomate.corpus.io.SWCDownloader`).

* Updated SWC-Reader (:class:`audiomate.corpus.io.SWCReader`) with an own implementation,
  so no manual preprocessing is needed anymore.

* Added conversion class (:class:`audiomate.corpus.conversion.WavAudioFileConverter`) to convert
  all files (or files that do not meet the requirements) of a corpus.

* Added writer (:class:`audiomate.corpus.io.NvidiaJasperWriter`) for
  `NVIDIA Jasper <https://github.com/NVIDIA/DeepLearningExamples/tree/master/PyTorch/SpeechRecognition/Jasper>`_.

* Create a consistent way to define invalid utterances of a dataset.
  Invalid utterance ids are defined in a json-file (e.g. ``audiomate/corpus/io/data/tuda/invalid_utterances.json``).
  Those are loaded automatically in the base-reader and can be accessed in the concrete implementation.

v5.0.0
------

**Breaking Changes**

* Changed :class:`audiomate.corpus.validation.InvalidItemsResult` to use it not only for Utterances, but also for Tracks for example.

* Refactoring and addition of splitting functions in the :class:`audiomate.corpus.subset.Splitter`.

**New Features**

* Added :class:`audiomate.corpus.validation.TrackReadValidator` to check for corrupt audio tracks/files.

* Added reader (:class:`audiomate.corpus.io.FluentSpeechReader`) for the
  `Fluent Speech Commands Dataset <http://www.fluent.ai/research/fluent-speech-commands/>`_.

* Added functions to check for contained tracks and issuers (:meth:`audiomate.corpus.CorpusView.contains_track`, :meth:`audiomate.corpus.CorpusView.contains_issuer`).

* Multiple options for controlling the behavior of the :class:`audiomate.corpus.io.KaldiWriter`.

* Added writer (:class:`audiomate.corpus.io.Wav2LetterWriter`) for the
  `wav2letter engine <https://github.com/facebookresearch/wav2letter/>`_.

* Added module with functions to read/write sclite trn files (:mod:`audiomate.formats.trn`).

**Fixes**

* Improved performance of Tuda-Reader (:class:`audiomate.corpus.io.TudaReader`).

* Added wrapper for the ```audioread.audio_open``` function (:mod:`audiomate.utils.audioread`) to cache available
  backends. This speeds up audioopen operations a lot.

* Performance improvements, especially for importing utterances, merging, subviews.

v4.0.1
------

**Fixes**

* Fix :class:`audiomate.corpus.io.CommonVoiceReader` to use correct file-extension of the audio files.

v4.0.0
------

**Breaking Changes**

* For utterances and labels ``-1`` was used for representing that the end is the same as the end of the parent utterance/track.
  In order to prevent ``-1`` checks in different methods/places ``float('inf')`` is now used.
  This makes it easier to implement stuff like label overlapping.

* :class:`audiomate.annotations.LabelList` is now backed by an interval-tree instead of a simple list. Therefore the labels have no fixed order anymore. The interval-tree provides functionality for operations like merging, splitting, finding overlaps with much lower code complexity.

* Removed module :mod:`audiomate.annotations.label_cleaning`, since those methods are available on :class:`audiomate.annotations.LabelList` directly.

**New Features**

* Added reader (:class:`audiomate.corpus.io.RouenReader`) and
  downloader (:class:`audiomate.corpus.io.RouenDownloader`) for the
  `LITIS Rouen Audio scene dataset <https://sites.google.com/site/alainrakotomamonjy/home/audio-scene>`_.

* Added downloader (:class:`audiomate.corpus.io.AEDDownloader`) for the
  `Acoustic Event Dataset <https://data.vision.ee.ethz.ch/cvl/ae_dataset/>`_.

* [`#69 <https://github.com/ynop/audiomate/issues/69>`_] Method to get labels within range: :meth:`audiomate.annotations.LabelList.labels_in_range`.

* [`#68 <https://github.com/ynop/audiomate/issues/68>`_] Add convenience method to create Label-List with list of label values: :meth:`audiomate.annotations.LabelList.with_label_values`.

* [`#61 <https://github.com/ynop/audiomate/issues/61>`_] Added function to split utterances of a corpus into multiple utterances with a maximal duration:
  :meth:`audiomate.corpus.CorpusView.split_utterances_to_max_time`.

* Add functions to check for overlap between labels: :meth:`audiomate.annotations.Label.do_overlap` and
  :meth:`audiomate.annotations.Label.overlap_duration`.

* Add function to merge equal labels that overlap within a label-list:
  :meth:`audiomate.annotations.LabelList.merge_overlapping_labels`.

* Added reader (:class:`audiomate.corpus.io.AudioMNISTReader`) and
  downloader (:class:`audiomate.corpus.io.AudioMNISTDownloader`) for the
  `AudioMNIST dataset <https://github.com/soerenab/AudioMNIST>`_.


**Fixes**

* [`#76 <https://github.com/ynop/audiomate/issues/76>`_][`#77 <https://github.com/ynop/audiomate/issues/77>`_][`#78 <https://github.com/ynop/audiomate/issues/78>`_] Multiple fixes on KaldiWriter


v3.0.0
------

**Breaking Changes**

* Moved label-encoding to its own module (:mod:`audiomate.encoding`).
  It now provides the processing of full corpora and store it in containers.

* Moved :class:`audiomate.feeding.PartitioningFeatureIterator` to the :mod:`audiomate.feeding` module.

* Added :class:`audiomate.containers.AudioContainer` to store audio tracks
  in a single file. All container classes are now in a separate module
  :mod:`audiomate.containers`.

* Corpus now contains Tracks not Files anymore. This makes it possible to
  different kinds of audio sources. Audio from a file is now included using
  :class:`audiomate.tracks.FileTrack`. New is the
  :class:`audiomate.tracks.ContainerTrack`, which reads data stored in
  a container.

* The :class:`audiomate.corpus.io.DefaultReader` and the
  :class:`audiomate.corpus.io.DefaultWriter` now load and store tracks,
  that are stored in a container.

* All functionality regarding labels was moved to its own module
  :mod:`audiomate.annotations`.

* The class :class:`audiomate.tracks.Utterance` was moved to the tracks module.

**New Features**

* Introducing the :mod:`audiomate.feeding` module. It provides different tools for accessing container data.
  Via a :class:`audiomate.feeding.Dataset` data can be accessed by indices.
  With a :class:`audiomate.feeding.DataIterator` one can easily iterate over data, such as frames.

* Added processing steps for computing Onset-Strength (:class:`audiomate.processing.pipeline.OnsetStrength`))
  and Tempogram (:class:`audiomate.processing.pipeline.Tempogram`)).

* Introduced :class:`audiomate.corpus.validation` module, that is used to validate a corpus.

* Added reader (:class:`audiomate.corpus.io.SWCReader`) for the
  `SWC corpus <https://audiomate.readthedocs.io/en/latest/documentation/indirect_support.html>`_.
  But it only works for the prepared corpus.

* Added function (:func:`audiomate.corpus.utils.label_cleaning.merge_consecutive_labels_with_same_values`)
  for merging consecutive labels with the same value

* Added downloader (:class:`audiomate.corpus.io.GtzanDownloader`) for the
  `GTZAN Music/Speech <https://marsyasweb.appspot.com/download/data_sets/>`_.

* Added :meth:`audiomate.corpus.assets.Label.tokenized` to get a list of tokens from a label.
  It basically splits the value and trims whitespace.

* Added methods on :class:`audiomate.corpus.CorpusView`, :class:`audiomate.corpus.assets.Utterance`
  and :class:`audiomate.corpus.assets.LabelList` to get a set of occurring tokens.

* Added :class:`audiomate.encoding.TokenOrdinalEncoder` to encode labels of an utterance
  by mapping every token of the label to a number.

* Create container base class (:class:`audiomate.corpus.assets.Container`), that can be used to store arbitrary data
  per utterance. The :class:`audiomate.corpus.assets.FeatureContainer` is now an extension of the container,
  that provides functionality especially for features.

* Added functions to split utterances and label-lists into multiple parts.
  (:meth:`audiomate.corpus.assets.Utterance.split`, :meth:`audiomate.corpus.assets.LabelList.split`)

* Added :class:`audiomate.processing.pipeline.AddContext` to add context to frames,
  using previous and subsequent frames.

* Added reader (:class:`audiomate.corpus.io.MailabsReader`) and
  downloader (:class:`audiomate.corpus.io.MailabsDownloader`) for the
  `M-AILABS Speech Dataset <http://www.m-ailabs.bayern/en/the-mailabs-speech-dataset/>`_.

**Fixes**

* [`#58 <https://github.com/ynop/audiomate/issues/58>`_] Keep track of number of samples per frame and between frames.
  Now the correct values will be stored in a Feature-Container, if the processor implements it correctly.

* [`#72 <https://github.com/ynop/audiomate/issues/72>`_] Fix bug, when reading samples from utterance,
  using a specific duration, while the utterance end is not defined.

v2.0.0
------

**Breaking Changes**

* Update various readers to use the correct label-list identifiers as defined
  in :ref:`data-mapping`.

**New Features**

* Added downloader (:class:`audiomate.corpus.io.TatoebaDownloader`) and
  reader (:class:`audiomate.corpus.io.TatoebaReader`) for the
  `Tatoeba platform <https://tatoeba.org/>`_.

* Added downloader (:class:`audiomate.corpus.io.CommonVoiceDownloader`) and
  reader (:class:`audiomate.corpus.io.CommonVoiceReader`) for the
  `Common Voice Corpus <https://voice.mozilla.org/>`_.

* Added processing steps :class:`audiomate.processing.pipeline.AvgPool` and
  :class:`audiomate.processing.pipeline.VarPool` for computing average and variance over
  a given number of sequential frames.

* Added downloader (:class:`audiomate.corpus.io.MusanDownloader`) for the
  `Musan Corpus <http://www.openslr.org/17/>`_.

* Added constants for common label-list identifiers/keys in :mod:`audiomate.corpus`.

v1.0.0
------

**Breaking Changes**

* The (pre)processing module has moved to :mod:`audiomate.processing`. It now supports online processing in chunks.
  For this purpose a pipeline step can require context.
  The pipeline automatically buffers data, until enough frames are ready.

**New Features**

* Added downloader (:class:`audiomate.corpus.io.FreeSpokenDigitDownloader`) and
  reader (:class:`audiomate.corpus.io.FreeSpokenDigitReader`) for the
  `Free-Spoken-Digit-Dataset <https://github.com/Jakobovski/free-spoken-digit-dataset>`_.


v0.1.0
------

Initial release
audiomate.encoding
==================

.. automodule:: audiomate.encoding
.. currentmodule:: audiomate.encoding

Encoder
-------

.. autoclass:: Encoder
   :members:

Frame-Based
-----------

.. autoclass:: FrameHotEncoder
   :members:

.. autoclass:: FrameOrdinalEncoder
   :members:

Utterance-Based
---------------

.. autoclass:: TokenOrdinalEncoder
   :members:
audiomate.processing
====================

.. automodule:: audiomate.processing
.. currentmodule:: audiomate.processing

Processor
---------

.. autoclass:: Processor
   :members:

Pipeline
--------

.. automodule:: audiomate.processing.pipeline
.. currentmodule:: audiomate.processing.pipeline

.. autoclass:: audiomate.processing.pipeline.Chunk
   :members:

.. autoclass:: audiomate.processing.pipeline.Step
   :members:

.. autoclass:: audiomate.processing.pipeline.Computation
   :members:

.. autoclass:: audiomate.processing.pipeline.Reduction
   :members:

Implementations
---------------

Some processing pipeline steps are already implemented.

.. _table-processing-step-implementations:

.. table:: Implementations of processing pipeline steps.

  ==============================  ===========
  Name                            Description
  ==============================  ===========
  MeanVarianceNorm                Normalizes features with given mean and variance.
  MelSpectrogram                  Exctracts MelSpectrogram features.
  MFCC                            Extracts MFCC features.
  PowerToDb                       Convert power spectrum to Db.
  Delta                           Compute delta features.
  AddContext                      Add previous and subsequent frames to the current frame.
  Stack                           Reduce multiple features into one by stacking them on top of each other.
  AvgPool                         Compute the average (per dimension) over a given number of sequential frames.
  VarPool                         Compute the variance (per dimension) over a given number of sequential frames.
  OnsetStrength                   Compute onset strengths.
  Tempogram                       Compute tempogram features.
  ==============================  ===========

.. autoclass:: audiomate.processing.pipeline.MeanVarianceNorm
   :members:

.. autoclass:: audiomate.processing.pipeline.MelSpectrogram
   :members:

.. autoclass:: audiomate.processing.pipeline.MFCC
   :members:

.. autoclass:: audiomate.processing.pipeline.PowerToDb
   :members:

.. autoclass:: audiomate.processing.pipeline.Delta
   :members:

.. autoclass:: audiomate.processing.pipeline.AddContext
   :members:

.. autoclass:: audiomate.processing.pipeline.Stack
   :members:

.. autoclass:: audiomate.processing.pipeline.AvgPool
   :members:

.. autoclass:: audiomate.processing.pipeline.VarPool
   :members:

.. autoclass:: audiomate.processing.pipeline.OnsetStrength
   :members:

.. autoclass:: audiomate.processing.pipeline.Tempogram
   :members:
audiomate.corpus.validation
===========================

.. automodule:: audiomate.corpus.validation
    :members:

Base
----

.. autoclass:: Validator
    :members:

.. autoclass:: ValidationResult
    :members:

.. autoclass:: InvalidItemsResult
    :members:

Combination
-----------

.. autoclass:: CombinedValidator
    :members:

.. autoclass:: CombinedValidationResult
    :members:

Label-List
----------

.. autoclass:: UtteranceTranscriptionRatioValidator
    :members:

.. autoclass:: LabelCountValidator
    :members:

.. autoclass:: LabelCoverageValidator
    :members:

.. autoclass:: LabelCoverageValidationResult
    :members:

.. autoclass:: LabelOverflowValidator
    :members:

.. autoclass:: LabelOverflowValidationResult
    :members:

Track
-----

.. autoclass:: TrackReadValidator
    :members:
audiomate.utils
===============

Audio
-----

.. automodule:: audiomate.utils.audio
    :members:

Audioread
---------

.. automodule:: audiomate.utils.audioread
    :members:

JSON File
---------

.. automodule:: audiomate.utils.jsonfile
    :members:

Naming
------

.. automodule:: audiomate.utils.naming
    :members:

Text
----

.. automodule:: audiomate.utils.text
    :members:

Text File
---------

.. automodule:: audiomate.utils.textfile
    :members:

Units
-----

.. automodule:: audiomate.utils.units
    :members:

Misc
----

.. automodule:: audiomate.utils.misc
    :members:
audiomate.corpus.conversion
===========================

.. automodule:: audiomate.corpus.conversion
.. currentmodule:: audiomate.corpus.conversion
audiomate.formats
=================

.. automodule:: audiomate.formats
.. currentmodule:: audiomate.formats

Audacity Labels
---------------

.. automodule:: audiomate.formats.audacity
   :members:

CTM Files
---------

.. automodule:: audiomate.formats.ctm
   :members:

TRN Files
---------

.. automodule:: audiomate.formats.trn
   :members:
audiomate.feeding
=================

.. automodule:: audiomate.feeding
.. currentmodule:: audiomate.feeding


Datasets
--------

.. autoclass:: Dataset
    :members:
    :inherited-members:

.. autoclass:: FrameDataset
    :members:
    :inherited-members:

.. autoclass:: MultiFrameDataset
    :members:
    :inherited-members:

Iterator
--------

.. autoclass:: DataIterator
    :members:
    :inherited-members:

.. autoclass:: FrameIterator
    :members:
    :inherited-members:

.. autoclass:: MultiFrameIterator
    :members:
    :inherited-members:

Partitioning
------------

.. autoclass:: PartitioningContainerLoader
    :members:
    :inherited-members:

.. autoclass:: PartitionInfo
    :members:
    :inherited-members:

.. autoclass:: PartitionData
    :members:
    :inherited-members:

.. autoclass:: PartitioningFeatureIterator
    :members:
    :inherited-members:
audiomate.tracks
================

.. automodule:: audiomate.tracks
.. currentmodule:: audiomate.tracks

Track
-----

.. autoclass:: Track
   :members:
   :inherited-members:

FileTrack
---------

.. autoclass:: FileTrack
   :members:
   :inherited-members:

ContainerTrack
--------------

.. autoclass:: ContainerTrack
   :members:
   :inherited-members:

Utterance
---------

.. autoclass:: Utterance
   :members:
   :inherited-members:
audiomate.corpus
================

.. automodule:: audiomate.corpus
.. currentmodule:: audiomate.corpus

CorpusView
----------

.. autoclass:: CorpusView
   :members:
   :inherited-members:

Corpus
------

.. autoclass:: Corpus
   :members:
audiomate.corpus.subset
=======================

.. automodule:: audiomate.corpus.subset
.. currentmodule:: audiomate.corpus.subset
audiomate.corpus.io
===================

.. automodule:: audiomate.corpus.io
    :members:
.. currentmodule:: audiomate.corpus.io

Base Classes
------------

.. autoclass:: CorpusDownloader
   :members:
   :inherited-members:
   :private-members:

.. autoclass:: ArchiveDownloader
   :members:
   :inherited-members:
   :private-members:

.. autoclass:: CorpusReader
   :members:
   :inherited-members:
   :private-members:

.. autoclass:: CorpusWriter
   :members:
   :inherited-members:
   :private-members:

.. _io_implementations:

Implementations
---------------

.. _table-format-support-of-readers-writers-by-format:

.. table:: Support for Reading and Writing by Format


  ================================  ========  =====  =======
  Format                            Download  Read   Write
  ================================  ========  =====  =======
  Acoustic Event Dataset            x         x
  AudioMNIST                        x         x
  Broadcast                                   x
  Common Voice                      x         x
  Default                                     x      x
  ESC-50                            x         x
  Free-Spoken-Digit-Dataset         x         x
  Folder                                      x
  Fluent Speech Commands Dataset              x
  Google Speech Commands                      x
  GTZAN                             x         x
  Kaldi                                       x      x
  LibriSpeech                       x         x
  Mozilla DeepSpeech                                 x
  MUSAN                             x         x
  M-AILABS Speech Dataset           x         x
  LITIS Rouen Audio scene dataset   x         x
  Spoken Wikipedia Corpora          x         x
  Tatoeba                           x         x
  TIMIT                                       x
  TUDA German Distant Speech        x         x
  Urbansound8k                                x
  VoxForge                          x         x
  Wav2Letter                                         x
  ================================  ========  =====  =======

Acoustic Event Dataset
^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: AEDDownloader
   :members:

.. autoclass:: AEDReader
   :members:

AudioMNIST
^^^^^^^^^^
.. autoclass:: AudioMNISTDownloader
   :members:

.. autoclass:: AudioMNISTReader
   :members:

Broadcast
^^^^^^^^^
.. autoclass:: BroadcastReader
   :members:

Common-Voice
^^^^^^^^^^^^
.. autoclass:: CommonVoiceDownloader
   :members:

.. autoclass:: CommonVoiceReader
   :members:

Default
^^^^^^^
.. autoclass:: DefaultReader
   :members:

.. autoclass:: DefaultWriter
   :members:

ESC-50
^^^^^^
.. autoclass:: ESC50Downloader
   :members:

.. autoclass:: ESC50Reader
   :members:

Folder
^^^^^^
.. autoclass:: FolderReader
   :members:

Free-Spoken-Digit-Dataset
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: FreeSpokenDigitDownloader
   :members:

.. autoclass:: FreeSpokenDigitReader
   :members:

Fluent Speech Commands Dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: FluentSpeechReader
   :members:

Google Speech Commands
^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: SpeechCommandsReader
   :members:

GTZAN
^^^^^
.. autoclass:: GtzanDownloader
   :members:

.. autoclass:: GtzanReader
   :members:

Kaldi
^^^^^
.. autoclass:: KaldiReader
   :members:

.. autoclass:: KaldiWriter
   :members:

LibriSpeech
^^^^^^^^^^^
.. autoclass:: LibriSpeechDownloader
   :members:

.. autoclass:: LibriSpeechReader
   :members:

Mozilla DeepSpeech
^^^^^^^^^^^^^^^^^^
.. autoclass:: MozillaDeepSpeechWriter
   :members:

MUSAN
^^^^^
.. autoclass:: MusanDownloader
   :members:

.. autoclass:: MusanReader
   :members:

M-AILABS Speech Dataset
^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: MailabsDownloader
   :members:

.. autoclass:: MailabsReader
   :members:

NVIDIA Jasper
^^^^^^^^^^^^^
.. autoclass:: NvidiaJasperWriter
   :members:

LITIS Rouen Audio scene dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: RouenDownloader
   :members:

.. autoclass:: RouenReader
   :members:

SWC - Spoken Wikipedia Corpora
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: SWCDownloader
   :members:

.. autoclass:: SWCReader
   :members:

Tatoeba
^^^^^^^
.. autoclass:: TatoebaDownloader
   :members:

.. autoclass:: TatoebaReader
   :members:

TIMIT DARPA Acoustic-Phonetic Continuous Speech Corpus
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: TimitReader
   :members:

TUDA German Distant Speech
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: TudaDownloader
   :members:

.. autoclass:: TudaReader
   :members:

Urbansound8k
^^^^^^^^^^^^
.. autoclass:: Urbansound8kReader
   :members:

VoxForge
^^^^^^^^

.. autoclass:: VoxforgeDownloader
   :members:

.. autoclass:: VoxforgeReader
   :members:

Wav2Letter
^^^^^^^^^^
.. autoclass:: Wav2LetterWriter
   :members:
audiomate.issuers
=================

.. automodule:: audiomate.issuers
.. currentmodule:: audiomate.issuers

Issuer
------

.. autoclass:: Issuer
   :members:
   :inherited-members:

Speaker
-------

.. autoclass:: Speaker
   :members:

Artist
------

.. autoclass:: Artist
   :members:

audiomate.annotations
=====================

.. automodule:: audiomate.annotations
.. currentmodule:: audiomate.annotations

Label
-----

.. autoclass:: Label
   :members:
   :inherited-members:

LabelList
---------

.. autoclass:: LabelList
   :members:
   :inherited-members:

Relabeling
----------

.. automodule:: audiomate.annotations.relabeling
    :members:

Exceptions
----------
.. autoexception:: audiomate.annotations.relabeling.UnmappedLabelsException
audiomate.containers
====================

.. automodule:: audiomate.containers
.. currentmodule:: audiomate.containers

Container
---------

.. autoclass:: Container
   :members:
   :inherited-members:

FeatureContainer
----------------

.. autoclass:: FeatureContainer
   :members:
   :inherited-members:

AudioContainer
----------------

.. autoclass:: AudioContainer
   :members:
   :inherited-members:
