---
title: 'Spleeter: a fast and efficient music source separation tool with pre-trained models'
tags:
  - Python
  - musical signal processing
  - source separation
  - vocal isolation
authors:
  - name:  Romain Hennequin
    orcid: 0000-0001-8158-5562
    affiliation: 1
  - name: Anis Khlif
    affiliation: 1
  - name: Felix Voituret
    affiliation: 1
  - name: Manuel Moussallam
    orcid: 0000-0003-0886-5423
    affiliation: 1
affiliations:
 - name: Deezer Research, Paris
   index: 1
date: 04 March 2020
bibliography: paper.bib

---

## Summary

We present and release a new tool for music source separation with pre-trained models called Spleeter. Spleeter was designed with ease of use, separation performance, and speed in mind. Spleeter is based on Tensorflow [@tensorflow2015-whitepaper] and makes it possible to:

- split music audio files into several stems with a single command line using pre-trained models. A music audio file can be separated into $2$ stems (vocals and accompaniments), $4$ stems (vocals, drums, bass, and other) or $5$ stems (vocals, drums, bass, piano and other).
- train source separation models or fine-tune pre-trained ones with Tensorflow (provided you have a dataset of isolated sources).

The performance of the pre-trained models are very close to the published state-of-the-art and is one of the best performing $4$ stems separation model on the common musdb18 benchmark [@musdb18] to be publicly released. Spleeter is also very fast as it can separate a mix audio file into $4$ stems $100$ times faster than real-time (we note, though, that the model cannot be applied in real-time as it needs buffering) on a single Graphics Processing Unit (GPU) using the pre-trained $4$-stems model.

## Purpose

We release Spleeter with pre-trained state-of-the-art models in order to help the Music Information Retrieval (MIR) research community leverage the power of source separation in various MIR tasks, such as vocal lyrics analysis from audio (audio/lyrics alignment, lyrics transcription...), music transcription (chord transcription, drums transcription, bass transcription, chord estimation, beat tracking), singer identification, any type of multilabel classification (mood/genre...), vocal melody extraction or cover detection.
We believe that source separation has reached a level of maturity that makes it worth considering for these tasks and that specific features computed from isolated vocals, drums or bass may help increase performances, especially in low data availability scenarios (small datasets, limited annotation availability) for which supervised learning might be difficult.
Spleeter also makes it possible to fine-tune the provided state-of-the-art models in order to adapt the system to a specific use-case.
Finally, having an available source separation tool such as Spleeter will allow researchers to compare performances of their new models to a state-of-the-art one on their private datasets instead of musdb18, which is usually the only used dataset for reporting separation performances for unreleased models.
Note that we cannot release the training data for copyright reasons, and thus, sharing pre-trained models were the only way to make these results available to the community.

## Implementation details

Spleeter contains pre-trained models for:

- vocals/accompaniment separation.
- $4$ stems separation as in SiSec [@SISEC18]  (vocals, bass, drums and other).
- $5$ stems separation with an extra piano stem (vocals, bass, drums, piano, and other). It is, to the authors' knowledge, the first released model to perform such a separation.

The pre-trained models are U-nets [@unet2017] and follow similar specifications as in [@deezerICASSP2019]. The U-net is an encoder/decoder Convolutional Neural Network (CNN) architecture with skip connections. We used $12$-layer U-nets ($6$ layers for the encoder and $6$ for the decoder). A U-net is used for estimating a soft mask for each source (stem). Training loss is a $L_1$-norm between masked input mix spectrograms and source-target spectrograms. The models were trained on Deezer's internal datasets (noteworthily the Bean dataset that was used in [@deezerICASSP2019]) using Adam [@Adam]. Training time took approximately a full week on a single GPU. Separation is then done from estimated source spectrograms using soft masking or multi-channel Wiener filtering.

Training and inference are implemented in Tensorflow which makes it possible to run the code on Central Processing Unit (CPU) or GPU.

## Speed

As the whole separation pipeline can be run on a GPU and the model is based on a CNN, computations are efficiently parallelized and model inference is very fast. For instance, Spleeter is able to separate the whole musdb18 test dataset (about $3$ hours and $27$ minutes of audio) into $4$ stems in less than $2$ minutes, including model loading time (about $15$ seconds), and audio wav files export, using a single GeForce RTX 2080 GPU, and a double Intel Xeon Gold 6134 CPU @ 3.20GHz (CPU is used for mix files loading and stem files export only). In this setup, Spleeter is able to process $100$ seconds of stereo audio in less than $1$ second, which makes it very useful for efficiently processing large datasets.

## Separation performances

The models compete with the state-of-the-art on the standard musdb18 dataset [@musdb18] while it was not trained, validated or optimized in any way with musdb18 data. We report results in terms of standard source separation metrics [@separation_metrics], namely Signal to Distortion Ratio (SDR), Signal to Artifacts Ratio (SAR), Signal to Interference Ratio (SIR) and source Image to Spatial distortion Ratio (ISR), are presented in the following table compared to Open-Unmix [@Open-Unmix] and Demucs [@demucs] (only SDR are reported for Demucs since other metrics are not available in the paper) which are, to the authors' knowledge, the only released system that performs near state-of-the-art performances.
We present results for soft masking and for multi-channel Wiener filtering (applied using Norbert [@Norbert]). As can be seen, for most metrics Spleeter is competitive with Open-Unmix and especially on SDR for all instruments, and is almost on par with Demucs.


|           |Spleeter Mask  |Spleeter MWF   |Open-Unmix |Demucs|
|-----------|---------------|---------------|-----------|------|
| Vocals SDR|6.55           |6.86           |6.32       |7.05  |
| Vocals SIR|15.19          |15.86          |13.33      |13.94 |
| Vocals SAR|6.44           |6.99           |6.52       |7.00  |
| Vocals ISR|12.01          |11.95          |11.93      |12.04 |
| Bass SDR  |5.10           |5.51           |5.23       |6.70  |
| Bass SIR  |10.01          |10.30          |10.93      |13.03 |
| Bass SAR  |5.15           |5.96           |6.34       |6.68  |
| Bass ISR  |9.18           |9.61           |9.23       |9.99  |
| Drums SDR |5.93           |6.71           |5.73       |7.08  |
| Drums SIR |12.24          |13.67          |11.12      |13.74 |
| Drums SAR |5.78           |6.54           |6.02       |7.04  |
| Drums ISR |10.50          |10.69          |10.51      |11.96 |
| Other SDR |4.24           |4.55           |4.02       |4.47  |
| Other SIR |7.86           |8.16           |6.59       |7.11  |
| Other SAR |4.63           |4.88           |4.74       |5.26  |
| Other ISR |9.83           |9.87           |9.31       |10.86 |

Spleeter [@spleeter] source code and pre-trained models are available on [github](https://www.github.com/deezer/spleeter) and distributed under a MIT license. This repository will eventually be used for releasing other models with improved performances or models separating into more than $5$ stems in the future.

## Distribution

Spleeter is available as a standalone Python package, and also provided as a [conda](https://github.com/conda-forge/spleeter-feedstock) recipe and self-contained [Dockers](https://hub.docker.com/r/researchdeezer/spleeter) which makes it usable as-is on various platforms.

## Acknowledgements

We acknowledge contributions from Laure Pretet who trained first models and wrote the first piece of code that lead to Spleeter.

## References
# Changelog History

## 2.3.0

Updating dependencies to enable TensorFlow 2.5 support (and Python 3.9 overall)
Removing the destructor from the `Separator` class

## 2.2.0

Minor changes mainly fixing some issues:
* mono training was not working due to hardcoded filters in the dataset
* default argument of `separate` was of wrong type
* added a way to request spleeter version with the `--version` argument in the CLI

## 2.1.0

This version introduce design related changes, especially transition to Typer for CLI managment and Poetry as
library build backend.

* `-i` option is now deprecated and replaced by traditional CLI input argument listing
* Project is now built using Poetry
* Project requires code formatting using Black and iSort
* Dedicated GPU package `spleeter-gpu` is not supported anymore, `spleeter` package will support both CPU and GPU hardware

### API changes:

* function `get_default_audio_adapter` is now available as `default()` class method within `AudioAdapter` class
* function `get_default_model_provider` is now available as `default()` class method within `ModelProvider` class
* `STFTBackend` and `Codec` are now string enum
* `GithubModelProvider` now use `httpx` with HTTP/2 support
* Commands are now located in `__main__` module, wrapped as simple function using Typer options module provide specification for each available option and argument
* `types` module provide custom type specification and must be enhanced in future release to provide more robust typing support with MyPy
* `utils.logging` module has been cleaned, logger instance is now a module singleton, and a single function is used to configure it with verbose parameter
* Added a custom logger handler (see tiangolo/typer#203 discussion)


## 2.0

First release, October 9th 2020

Tensorflow-2 compatible version, allowing uses in python 3.8.

## 1.5.4

First release, July 24th 2020

Add some padding of the input waveform to avoid separation artefacts on the edges due to unstabilities in the inverse fourier transforms.
Also add tests to ensure both librosa and tensorflow backends have same outputs.

## 1.5.2

First released, May 15th 2020

### Major changes

* PR #375 merged to avoid mutliple tf.graph instantiation failures

### Minor changes

* PR #362 use tf.abs instead of numpy
* PR #352 tempdir cleaning


## 1.5.1

First released, April 15th 2020

### Major changes

* Bugfixes on the LibRosa STFT backend

### Minor changes

* Typos, and small bugfixes

## 1.5.0

First released, March 20th 2020

### Major changes

* Implement a new STFT backend using LibRosa, faster on CPU than TF implementation
* Switch tensorflow version to 1.15.2

### Minor changes

* Typos, and small bugfixes

## 1.4.9

First released, Dec 27th 2019

### Major changes

* Add new configuration for processing until 16Khz

### Minor changes

* Typos, and small bugfixes
<img src="https://github.com/deezer/spleeter/raw/master/images/spleeter_logo.png" height="80" />

[![Github actions](https://github.com/deezer/spleeter/workflows/pytest/badge.svg)](https://github.com/deezer/spleeter/actions) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/spleeter) [![PyPI version](https://badge.fury.io/py/spleeter.svg)](https://badge.fury.io/py/spleeter) [![Conda](https://img.shields.io/conda/vn/deezer-research/spleeter)](https://anaconda.org/deezer-research/spleeter) [![Docker Pulls](https://img.shields.io/docker/pulls/deezer/spleeter)](https://hub.docker.com/r/researchdeezer/spleeter) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/deezer/spleeter/blob/master/spleeter.ipynb) [![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/spleeter/community) [![status](https://joss.theoj.org/papers/259e5efe669945a343bad6eccb89018b/status.svg)](https://joss.theoj.org/papers/259e5efe669945a343bad6eccb89018b)

> :warning: [Spleeter 2.1.0](https://pypi.org/project/spleeter/) release introduces some breaking changes, including new CLI option naming for input, and the drop
> of dedicated GPU package. Please read [CHANGELOG](CHANGELOG.md) for more details.

## About

**Spleeter** is [Deezer](https://www.deezer.com/) source separation library with pretrained models
written in [Python](https://www.python.org/) and uses [Tensorflow](https://tensorflow.org/). It makes it easy
to train source separation model (assuming you have a dataset of isolated sources), and provides
already trained state of the art model for performing various flavour of separation :

* Vocals (singing voice) / accompaniment separation ([2 stems](https://github.com/deezer/spleeter/wiki/2.-Getting-started#using-2stems-model))
* Vocals / drums / bass / other separation ([4 stems](https://github.com/deezer/spleeter/wiki/2.-Getting-started#using-4stems-model))
* Vocals / drums / bass / piano / other separation ([5 stems](https://github.com/deezer/spleeter/wiki/2.-Getting-started#using-5stems-model))

2 stems and 4 stems models have [high performances](https://github.com/deezer/spleeter/wiki/Separation-Performances) on the [musdb](https://sigsep.github.io/datasets/musdb.html) dataset. **Spleeter** is also very fast as it can perform separation of audio files to 4 stems 100x faster than real-time when run on a GPU.

We designed **Spleeter** so you can use it straight from [command line](https://github.com/deezer/spleeter/wiki/2.-Getting-started#usage)
as well as directly in your own development pipeline as a [Python library](https://github.com/deezer/spleeter/wiki/4.-API-Reference#separator). It can be installed with [pip](https://github.com/deezer/spleeter/wiki/1.-Installation#using-pip) or be used with
[Docker](https://github.com/deezer/spleeter/wiki/2.-Getting-started#using-docker-image).

### Projects and Softwares using **Spleeter**

Since it's been released, there are multiple forks exposing **Spleeter** through either a Guided User Interface (GUI) or a standalone free or paying website. Please note that we do not host, maintain or directly support any of these initiatives.

That being said, many cool projects have been built on top of ours. Notably the porting to the *Ableton Live* ecosystem through the [Spleeter 4 Max](https://github.com/diracdeltas/spleeter4max#spleeter-for-max) project.

**Spleeter** pre-trained models have also been used by professionnal audio softwares. Here's a non-exhaustive list:

* [iZotope](https://www.izotope.com/en/shop/rx-8-standard.html) in its *Music Rebalance* feature within **RX 8**
* [SpectralLayers](https://new.steinberg.net/spectralayers/) in its *Unmix* feature in **SpectralLayers 7**
* [Acon Digital](https://acondigital.com/products/acoustica-audio-editor/) within **Acoustica 7**
* [VirtualDJ](https://www.virtualdj.com/stems/) in their stem isolation feature
* [Algoriddim](https://www.algoriddim.com/apps) in their **NeuralMix** and **djayPRO** app suite

🆕 **Spleeter** is a baseline in the ongoing [Music Demixing Challenge](https://www.aicrowd.com/challenges/music-demixing-challenge-ismir-2021)!

## Quick start

Want to try it out but don't want to install anything ? We have set up a [Google Colab](https://colab.research.google.com/github/deezer/spleeter/blob/master/spleeter.ipynb).

Ready to dig into it ? In a few lines you can install **Spleeter**  and separate the vocal and accompaniment parts from an example audio file.
You need first to install `ffmpeg` and `libsndfile`. It can be done on most platform using [Conda](https://github.com/deezer/spleeter/wiki/1.-Installation#using-conda):

```bash
# install dependencies using conda
conda install -c conda-forge ffmpeg libsndfile
# install spleeter with pip
pip install spleeter
# download an example audio file (if you don't have wget, use another tool for downloading)
wget https://github.com/deezer/spleeter/raw/master/audio_example.mp3
# separate the example audio into two components
spleeter separate -p spleeter:2stems -o output audio_example.mp3
```

> :warning: Note that we no longer recommend using `conda` for installing spleeter.

> ⚠️ There are known issues with Apple M1 chips, mostly due to TensorFlow compatibility. Until these are fixed, you can use [this workaround](https://github.com/deezer/spleeter/issues/607#issuecomment-828352392)

You should get two separated audio files (`vocals.wav` and `accompaniment.wav`) in the `output/audio_example` folder.

For a detailed documentation, please check the [repository wiki](https://github.com/deezer/spleeter/wiki/1.-Installation)

## Development and Testing

This project is managed using [Poetry](https://python-poetry.org/docs/basic-usage/), to run test suite you
can execute the following set of commands:

```bash
# Clone spleeter repository
git clone https://github.com/Deezer/spleeter && cd spleeter
# Install poetry
pip install poetry
# Install spleeter dependencies
poetry install
# Run unit test suite
poetry run pytest tests/
```

## Reference

* Deezer Research - Source Separation Engine Story - deezer.io blog post:
  * [English version](https://deezer.io/releasing-spleeter-deezer-r-d-source-separation-engine-2b88985e797e)
  * [Japanese version](http://dzr.fm/splitterjp)
* [Music Source Separation tool with pre-trained models / ISMIR2019 extended abstract](http://archives.ismir.net/ismir2019/latebreaking/000036.pdf)

If you use **Spleeter** in your work, please cite:

```BibTeX
@article{spleeter2020,
  doi = {10.21105/joss.02154},
  url = {https://doi.org/10.21105/joss.02154},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {50},
  pages = {2154},
  author = {Romain Hennequin and Anis Khlif and Felix Voituret and Manuel Moussallam},
  title = {Spleeter: a fast and efficient music source separation tool with pre-trained models},
  journal = {Journal of Open Source Software},
  note = {Deezer Research}
}
```

## License

The code of **Spleeter** is [MIT-licensed](LICENSE).

## Disclaimer

If you plan to use **Spleeter** on copyrighted material, make sure you get proper authorization from right owners beforehand.

## Troubleshooting

**Spleeter** is a complex piece of software and although we continously try to improve and test it you may encounter unexpected issues running it. If that's the case please check the [FAQ page](https://github.com/deezer/spleeter/wiki/5.-FAQ) first as well as the list of [currently open issues](https://github.com/deezer/spleeter/issues)

### Windows users

   It appears that sometimes the shortcut command `spleeter` does not work properly on windows. This is a known issue that we will hopefully fix soon. In the meantime replace `spleeter separate` by `python -m spleeter separate` in command line and it should work.

## Contributing

If you would like to participate in the development of **Spleeter** you are more than welcome to do so. Don't hesitate to throw us a pull request and we'll do our best to examine it quickly. Please check out our [guidelines](.github/CONTRIBUTING.md) first.

## Note

This repository include a demo audio file `audio_example.mp3` which is an excerpt
from Slow Motion Dream by Steven M Bryant (c) copyright 2011 Licensed under a Creative
Commons Attribution (3.0) [license](http://dig.ccmixter.org/files/stevieb357/34740)
Ft: CSoul,Alex Beroza & Robert Siekawitch
# How-to contribute

Those are the main contributing guidelines for contributing to this project:

- Verify that your contribution does not embark proprietary code or infringe any copyright of any sort.
- Avoid adding any unnecessary dependencies to the project, espcially of those are not easily packaged and installed through `conda` or `pip`.
- Python contributions must follow the [PEP 8 style guide](https://www.python.org/dev/peps/pep-0008/).
- Use [Pull Request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests) mechanism and please be patient while waiting for reviews.
- Remain polite and civil in all exchanges with the maintainers and other contributors.
- Any issue submitted which does not respect provided template, or lack of information, will be considered as invalid and automatically closed.

## Get started

This project is managed using [Poetry](https://python-poetry.org/docs/basic-usage/),
in order to contribute, the safest is to create your
[own fork of spleeter](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) first and then setup your development environment:

```bash
# Clone spleeter repository fork
git clone https://github.com/<your_name>/spleeter && cd spleeter
# Install poetry
pip install poetry
# Install spleeter dependencies
poetry install
# Run unit test suite
poetry run pytest tests/
```

You can then make your changes and experiment freely. Once you're done, remember to check that the tests still run. If you've added a new feature, add tests!

Then finally, you're more than welcome to create a [Pull Request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork) in **Spleeter** main repo. We will look at it as soon as possible and eventually integrate your changes in the project.

## PR requirements

Following command should be ran successfully before to consider a PR for merging:

```bash
poetry run pytest tests/
poetry run black spleeter
poetry run isort spleeter
```
# Pull request title

- [ ] I read [contributing guideline](https://github.com/deezer/spleeter/blob/master/.github/CONTRIBUTING.md)
- [ ] I didn't find a similar pull request already open.
- [ ] My PR is related to Spleeter only, not a derivative product (such as Webapplication, or GUI provided by others)

## Description

A few sentences describing the overall goals of the pull request's commits.

## How this patch was tested

You tested it, right?

- [ ] I implemented unit test whicn ran successfully using `poetry run pytest tests/`
- [ ] Code has been formatted using `poetry run black spleeter`
- [ ] Imports has been formatted using `poetry run isort spleeter``

## Documentation link and external references

Please provide any info that may help us better understand your code.
---
name: Bug
about: Report a bug
title: "[Bug] name your bug"
labels: bug, invalid
---

- [ ] I didn't find a similar issue already open.
- [ ] I read the documentation (README AND Wiki)
- [ ] I have installed FFMpeg
- [ ] My problem is related to Spleeter only, not a derivative product (such as Webapplication, or GUI provided by others)

## Description

<!-- Give us a clear and concise description of the bug you are reporting. -->

## Step to reproduce

<!-- Indicates clearly steps to reproduce the behavior: -->

1. Installed using `...`
2. Run as `...`
3. Got `...` error

## Output

```bash
Share what your terminal says when you run the script (as well as what you would expect).
```

## Environment

<!-- Fill the following table -->

|                   |                                 |
| ----------------- | ------------------------------- |
| OS                | Windows / Linux / MacOS / other |
| Installation type | Conda / pip / other             |
| RAM available     | XGo                             |
| Hardware spec     | GPU / CPU / etc ...             |

## Additional context

<!-- Add any other context about the problem here, references, cites, etc.. -->
---
name: Feature request
about: Submit idea for new feature
labels: feature, enhancement
title: "[Feature] your feature name"
---

## Description

<!-- Describe your feature request here. -->

## Additional information

<!-- Add any additional description -->
---
name: Discussion
about: Ideas sharing or theorical question solving 
labels: question
title: "[Discussion] your question"
---

<!-- Please respect the title [Discussion] tag. -->
