---
title: "Omnizart: A General Toolbox for Automatic Music Transcription"
tags:
  - Python
  - automatic music transcription
  - music information retrieval
  - audio signal processing
  - artificial intelligence
authors:
  - name: Yu-Te Wu
    affiliation: 1
  - name: Yin-Jyun Luo
    affiliation: 1
  - name: Tsung-Ping Chen
    affiliation: 1
  - name: I-Chieh Wei
    affiliation: 1
  - name: Jui-Yang Hsu
    affiliation: 1
  - name: Yi-Chin Chuang
    affiliation: 1
  - name: Li Su
    affiliation: 1
affiliations:
  - name: Music and Culture Technology Lab, Institute of Information Science, Academia Sinica, Taipei, Taiwan
    index: 1
date: 18 April 2021
bibliography: paper.bib
# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

We present and release Omnizart, a new Python library that provides a streamlined solution to automatic music transcription (AMT).
Omnizart encompasses modules that construct the life-cycle of deep learning-based AMT, and is designed for ease of use with a compact command-line interface.
To the best of our knowledge, Omnizart is the first toolkit that offers transcription models for various music content including piano solo, instrument ensembles, percussion and vocal. Omnizart also supports models for chord recognition and beat/downbeat tracking, which are highly related to AMT.

In summary, Omnizart incorporates:

- Pre-trained models for frame-level and note-level transcription of multiple pitched instruments, vocal melody, and drum events;
- Pre-trained models of chord recognition and beat/downbeat tracking;
- The main functionalities in the life-cycle of AMT research, covering dataset downloading, feature pre-processing, model training, to the sonification of the transcription result.

Omnizart is based on Tensorflow [@abadi2016tensorflow].
The complete code base, command-line interface, documentation, as well as demo examples can all be accessed from the [project website](https://github.com/Music-and-Culture-Technology-Lab/omnizart).

# Statement of need

AMT of polyphonic music is a complicated MIR task because the note-, melody-, timbre-, and rhythm-level attributes of music are overlapped with each other in music signals. A unified solution of AMT is therefore in eager demand. AMT is also strongly related to other MIR tasks such as source separation and music generation with transcribed data needed as supervisory resources.
Omnizart considers multi-instrument transcription and collects several state-of-the-art models for transcribing pitched and percussive instruments, as well as singing voice, within polyphonic music signals. Omnizart is an AMT tool that unifies multiple transcription utilities and enables further productivity. Omnizart can save one's time and labor in generating a massive number of multi-track MIDI files, which could have a large impact on music production, music generation, education, and musicology research.

# Implementation Details

## Piano solo transcription

The piano solo transcription model in Omnizart reproduces the implementation of @wu2020multi.
The model features a U-net that takes as inputs the audio spectrogram, generalized cepstrum (GC) [@su2015combining], and GC of spectrogram (GCoS) [@wu2018automatic], and outputs a multi-channel time-pitch representation with time- and pitch-resolution of 20 ms and 25 cents, respectively.
For the U-net, implementation of the encoder and the decoder follows DeepLabV3+ [@Chen2018DeepLabV3], and the bottleneck layer is adapted from the Image Transformer [@parmar2018image].

The model is trained on the MAESTRO dataset [@hawthorne2018enabling], an external dataset containing 1,184 real piano performance recordings with a total length of 172.3 hours.
The model achieves 72.50\% and 79.57\% for frame- and note-level F1-scores, respectively, on the Configuration-II test set of the MAPS dataset [@kelz2016potential].

## Multi-instrument polyphonic transcription

The multi-instrument transcription model extends the piano solo model to support 11 output classes, namely piano, violin, viola, cello, flute, horn, bassoon, clarinet, harpsichord, contrabass, and oboe, accessed from MusicNet [@thickstun2017learning].
Detailed characteristics of the model can be seen in @wu2020multi.
The evaluation on the test set from MusicNet [@thickstun2018invariances] yields 66.59\% for the note streaming task.

## Drum transcription

The model for drum transcription is a re-implementation of @wei2021improving.
Building blocks of the network include convolutional layers and the attention mechanism.

The model is trained on a dataset with 1,454 audio clips of polyphonic music with synchronized drum events [@wei2021improving].
The model demonstrates SoTA performance on two commonly used benchmark datasets, i.e., 74\% for ENST [@gillet2006enst] and 71\% for MDB-Drums [@southall2017mdb] in terms of the note-level F1-score.

## Vocal transcription in polyphonic music

The system for vocal transcription features a pitch extractor and a module for note segmentation.
The inputs to the model are composed of spectrogram, GS, and GCoS derived from polyphonic music recordings [@wu2018automatic].

A pre-trained Patch-CNN [@su2018vocal] is leveraged as the pitch extractor.
The module for note segmentation is implemented with PyramidNet-110 and ShakeDrop regularization [@yamada2019shakedrop], which is trained using Virtual Adversarial Training [@miyato2018virtual] enabling semi-supervised learning.

The training data includes labeled data from TONAS [@mora2010characterization] and unlabeled data from MIR-1K [@hsu2009improvement].
The model yields the SoTA F1-score of 68.4\% evaluated with the ISMIR2014 dataset [@molina2014evaluation].

## Chord recognition

The harmony recognition model of Omnizart is implemented using the Harmony Transformer (HT) [@chen2019harmony].
The HT model is based on an encoder-decoder architecture,
where the encoder performs chord segmentation on the input, and the decoder recognizes the chord progression based on the segmentation result.

The original HT supports both audio and symbolic inputs.
Currently, Omnizart supports only audio inputs.
A given audio input is pre-processed using Chordino VAMP plugin [@mauch2010approximate] as the non-negative-least-squares chromagram.
The outputs of the model include 25 chord types, covering 12 major and minor chords together with a class referred to the absence of chord, with a time resolution of 230 ms.

In an experiment with evaluations on the McGill Billboard dataset [@burgoyne2011anexpert], the HT outperforms the previous state of the art [@chen2019harmony].

## Beat/downbeat tracking

The model for beat and downbeat tracking provided in Omnizart is a reproduction of @chuang2020beat.
Unlike most of the available open-source projects such as \texttt{madmom} [@bock2016madmom] and \texttt{librosa} [@mcfee2015librosa] which focus on audio, the provided model targets symbolic data.

The input and output of the model are respectively MIDI and beat/downbeat positions with the time resolution of 10 ms.
The input representation combines piano-roll, spectral flux, and inter-onset interval extracted from MIDI.
The model composes a two-layer BLSTM network with the attention mechanism, and predicts probabilities of the presence of beat and downbeat per time step.

Experiments on the MusicNet dataset [@thickstun2018invariances] with the synchronized beat annotation show that the proposed model outperforms the state-of-the-art beat trackers which operate on synthesized audio [@chuang2020beat].

# Conclusion

Omnizart represents the first systematic solution for the polyphonic AMT of general music contents ranging from pitched instruments, percussion instruments, to voices.
In addition to note transcription, Omnizart also includes high-level MIR tasks such as chord recognition and beat/downbeat tracking.
As an ongoing project, the research group will keep refining the package and extending the scope of transcription in the future.

# References
# Changelog

## 0.5.0 - 2021-12-09

Official Open JOSS reviewed version.

## Bugs
- Fix bug of name conflict while loading chord model.


## 0.4.2 - 2021-11-16

Accumulated release. Various improvements and bug fix. See details below.

## Feature
- Migrate checkpoints from private Google Drive to Github release. 
See [here](https://github.com/Music-and-Culture-Technology-Lab/omnizart/releases/tag/checkpoints-20211001)
- Replace opencv 

## Dependency
- Upgrade Tensorflow version to 2.5.0 for Nvidia 30 series GPU compatibility.
- Upgrade Spleeter version to 2.3.0 for new TF version compatibility.
- Replace Opencv with PIL for drum feature resizing and remove opencv from the dependency.

## Enhancement
- Simplify model loading mechanism by unifying the all checkpoint format to use TF format.
- Lazy import extraction functions to boost loading time.
- Change the order of Dockerfile commands for better utilizing cache.

## Documentation
- Add notice about compatibility issue of running on certain CPU architecture.
- Add explaination about enabling auto completion.
- Rephrase sentences in paper according to JOSS review feedback.
- Add explaination about installing development dependencies.
- Use pepy as the alternative source for 'download' state badge.


## Bugs
- Fix bug of unable to find vocal contour checkpoint.
- Fix bug of fail to custom layers of chord module.
- Fix various unit tests bugs.
- Fix minor linter errors.



## 0.4.1 - 2021-06-04
Hotfix version according to issue [#23](https://github.com/Music-and-Culture-Technology-Lab/omnizart/issues/23)

## Feature
- Add a new piano transcription model and set it as the default model while using `music` module.

## Bugs
- Fix bug while parsing weight files in the checkpoint folder.

---

## 0.4.0 - 2021-05-31
Various improvements on music module and some critical bug fixes.

## Enhancement
- Improve the peak finding and thresholding strategy for more stable and better performance.
- Modify the feeding strategy of feature slices with adjustable overlapping rate while making predictions.
- Apply learning rate scheduler for music module.
- Replace the usage of custom training loop of music module with the built-in TF `.fit()` function.

## Bugs
- Fix a critical bug of inference of music module that would lead to missing onsets.
- Fix generation of pertubation of vocal module while training.

## Documentation
- Merge the demo page into master from `build_doc` branch.

---

## 0.3.4 - 2021-05-10
Hotifx version according to issue #19.

## Bugs
- Fix bug of treating numpy array as list while appending elements.

---

## 0.3.3 - 2021-05-07
Hotfix version according to issue #19.

## Bugs
- Fix column inconsistency of `aggregate_f0_info` and `write_agg_f0_results`.
- Update version of dependencies according to the security alert.

---

## 0.3.2 - 2021-02-13

### Enhancement
- Move `load_label` functions of different datasets into dataset structure classes.
- Add custom exception on fail downloading GD file due to access limit.
- Add unit tests on parsing label files into shared intermediate format.

### Bugs
- Fix wrong access name of the dict in vocal midi inference function.
- Fix bug of generating beat module training labels.

---

## 0.3.1 - 2021-01-18

Hotfix release of spleeter error.

### Bugs
- Call Spleeter in CLI mode instead of using python class.

---

## 0.3.0 - 2021-01-17

Release the `beat` module for symbolic domain beat transcription.

### Features
- Release `beat` module.
- Add an example `patch-cnn` module for demonstrating the implementation progress.

### Enhancement
- Refactor the flow of chord module for parsing the feature and label files.
- Modularize F0 information aggragation functions to *utils.py* and *io.py*.
- Improve verbosity on fail to open hdf files.

### Documentation
- Re-arrange the side bar with an additional group of CLI.
- Add custom CSS style for adjusting the width of audio and video elements.

### Bugs
- Fix Spleeter import errors after upgrading to v2.1.2.

---

## 0.2.0 - 2020-12-13

### Vocal melody transcription in both frame- and note-level are live!
We release the modules for vocal melody transcription after a decent amount of effort. 
Now you can transcribe your favorite singing voice.

### Features
- Release `vocal` and `vocal-contour` submodules.

### Enhancement
- Improve chord transcription results by filtering out chord predictions with short duration.
- Resolve the path for transcription output in a consistent way.

### Documentation
- Re-organize Quick Start and Tutorial pages to improve accessibility.
- Move the section for development from README.md to CONTRIBUTING.md.

### Bug Fix
- Fix bugs of passing the wrong parameter to vamp for chroma feature extraction.

---

## 0.1.1 - 2020-12-01
### Features
- Add more supported datasets for download and process.
- Supports to save checkpoints in .pb format with customized model checkpoint callback.

### Enhancement
- Huge refactor of constants.dataset. Improves reusability and add more useful common utilities.
- Modularize common parts of app classes.
- Construct base class of loading dataset samples. Reduce duplicate code and reuse the same functionalities.
- Filter out messy Tensorflow warnings when using CLI.

### Bug Fix
- Resolved bugs of some function parameters not actually being used inside functions.
- Fix CFP extraction down_fs don't actually work.

---

## 0.1.0 - 2020-11-16
### Features
- Add command for synthesizing MIDI file.
- Provides colab for quick start now!

### Enhancement
- Lazy import application instance for avoiding pulling large amount of dependencies.
- Group sub-commands into different sections when showing help message.

---

## 0.1.0-beta.2 - 2020-11-10

### Enhancement
- Better dealing with the input model path.
- Better approach for resolving dataset path when given with "./".
- Add documentation for Conda user for manually install omnizart.

### Bug Fix
- Fix wrong save path of checkpoints.
- Fix installation script for not upgrading pip after activating virtual environment.

---

## 0.1.0-beta.1 - 2020-11-08

First release of `omnizart` CLI tool, as well as a python package.

### Features
- Multi-instrument transcription
- Drum transcription
- Chord transcription
- Download datasets
- Extract feature of datasets for each module
- Train models for each module

---
# OMNIZART

[![build](https://github.com/Music-and-Culture-Technology-Lab/omnizart/workflows/general-check/badge.svg)](https://github.com/Music-and-Culture-Technology-Lab/omnizart/actions?query=workflow%3Ageneral-check)
[![docs](https://github.com/Music-and-Culture-Technology-Lab/omnizart/workflows/docs/badge.svg?branch=build_doc)](https://music-and-culture-technology-lab.github.io/omnizart-doc/)
[![PyPI version](https://badge.fury.io/py/omnizart.svg)](https://pypi.org/project/omnizart/)
![PyPI - License](https://img.shields.io/pypi/l/omnizart)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/omnizart)](https://pypistats.org/packages/omnizart)
[![Docker Pulls](https://img.shields.io/docker/pulls/mctlab/omnizart)](https://hub.docker.com/r/mctlab/omnizart)

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03391/status.svg)](https://doi.org/10.21105/joss.03391)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5769022.svg)](https://doi.org/10.5281/zenodo.5769022)


Omnizart is a Python library that aims for democratizing automatic music transcription.
Given polyphonic music, it is able to transcribe pitched instruments, vocal melody, chords, drum events, and beat.
This is powered by the research outcomes from [Music and Culture Technology (MCT) Lab](https://sites.google.com/view/mctl/home). The paper has been published to [Journal of Open Source Software (JOSS)](https://doi.org/10.21105/joss.03391).

### Transcribe your favorite songs now in Colab [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://bit.ly/OmnizartColab) or [![Replicate](https://replicate.com/breezewhite/omnizart/badge)](https://replicate.ai/breezewhite/omnizart)

# Quick start

Visit the [complete document](https://music-and-culture-technology-lab.github.io/omnizart-doc/) for detailed guidance.

## Pip
``` bash
# Install omnizart
pip install omnizart

# Download the checkpoints
omnizart download-checkpoints

# Transcribe your songs
omnizart drum transcribe <path/to/audio.wav>
omnizart chord transcribe <path/to/audio.wav>
omnizart music transcribe <path/to/audio.wav>
```

## Docker
```
docker pull mctlab/omnizart:latest
docker run -it mctlab/omnizart:latest bash
```

# Supported applications
| Application      | Transcription      | Training           | Evaluation | Description                                      |
|------------------|--------------------|--------------------|------------|--------------------------------------------------|
| music            | :heavy_check_mark: | :heavy_check_mark: |            | Transcribe musical notes of pitched instruments. |
| drum             | :heavy_check_mark: | :interrobang:      |            | Transcribe events of percussive instruments.     |
| vocal            | :heavy_check_mark: | :heavy_check_mark: |            | Transcribe note-level vocal melody.              |
| vocal-contour    | :heavy_check_mark: | :heavy_check_mark: |            | Transcribe frame-level vocal melody (F0).        |
| chord            | :heavy_check_mark: | :heavy_check_mark: |            | Transcribe chord progressions.                   |
| beat             | :heavy_check_mark: | :heavy_check_mark: |            | Transcribe beat position.                        |

**NOTES**
The current implementation for the drum model has unknown bugs, preventing loss convergence when training from scratch.
Fortunately, you can still enjoy drum transcription with the provided checkpoints.

## Compatibility Issue
Currently, Omnizart is **incompatible for ARM-based MacOS** system due to the underlying dependencies.
More details can be found in the [issue #38](https://github.com/Music-and-Culture-Technology-Lab/omnizart/issues/38).

## Citation
If you use this software in your work, please cite:

```
@article{Wu2021,
  doi = {10.21105/joss.03391},
  url = {https://doi.org/10.21105/joss.03391},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {68},
  pages = {3391},
  author = {Yu-Te Wu and Yin-Jyun Luo and Tsung-Ping Chen and I-Chieh Wei and Jui-Yang Hsu and Yi-Chin Chuang and Li Su},
  title = {Omnizart: A General Toolbox for Automatic Music Transcription},
  journal = {Journal of Open Source Software}
}
```

# Development
Describes the neccessary background of how to develop this project.

## Download and install
``` bash
git clone https://github.com/Music-and-Culture-Technology-Lab/omnizart.git

# Install dependenies. For more different installation approaches, please refer to the official documentation page.
# The following command will download the checkpoints automatically.
cd omnizart
make install

# For developers, you have to install Dev dependencies as well, since they will not be installed by default.
make install-dev
```

## Package management
Uses [poetry](https://python-poetry.org/) for package management, instead of writing `requirements.txt` and `setup.py` manually.
We still provide the above two files for convenience. You can also generate them by executing ``make export``.

### ATTENTION! MUST SEE!
There is a major difference between install with `poetry install` and `python setup.py install`. When using poetry for installation, which
is the default approach when running `make insatll`, the site-pacakges and resource files are placed in the **current** folder.
This is different from executing `python setup.py install`, which resource files are installed in where a normal package you download through `pip install` will be placed (e.g. ~/.local/lib/python3.6/site-packages) .

And why things aren't placed in the normal path, but the command still can be executed? The answer is that poetry add an additional package path to your *PATH*  environment variable, and guess what is that path? Bingo! Your current path where you can execute `poetry install`! The difference has a major impact on running `omnizart download-checkpoints`. The default save path of checkpoints is to where omnizart being installed. That would be fine for end users, but not good news for developers though. That means after you git clone this project, and installed with `setup.py` approach, the **checkpoints are stored under ~/.local/.../site-packages/**, not your current development path. Therefore, it is strongly suggested that developers should use the default installation approach for a more comfortable developing experience^^.

Feedback: what a big trap there is...


## Documentation
Automatically generate documents from inline docstrings of module, class, and function. 
[Hosted document page](https://music-and-culture-technology-lab.github.io/omnizart-doc/index.html)

Documentation style: Follows `numpy` document flavor. Learn more from [numpydoc](https://numpydoc.readthedocs.io/en/latest/format.html).

Document builder: [sphinx](https://www.sphinx-doc.org/en/master/)

To generate documents, `cd docs/` and execute `make html`. To see the rendered results, run `make serve` and view from the browser.
All documents and docstrings use **reStructured Text** format. More informations about this format can be found from 
[Sphinx's Document](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html).

## Linters
Uses flake8 and pylint for coding style check.

To check with linters, execute `make check`.

You don't have to achieve a perfect score on pylint check, just pass 9.5 points still counted as a successful check.

### Caution!
There is convenient make command for formating the code, but it should be used very carefully.
Not only it could format the code for you, but also could mess up the code, and eventually you should still need
to check the format manually after refacorting the code with tools. 

To format the code with black and yapf, enter `make format`.

## Unittest
Uses `pytest` for unittest. The overall coverage rate should pass 25%, or CI would fail.

## CI/CD
Uses github actions for automatic linting, unittesting, document building, and package release.
Currently supports two workflows:
* General check
* Documentation page publishing
* Publish PyPI package and docker image

### General Check
Everytime you push to the master branch, file a pull request, and merge into master branch, will trigger
this action. This will do checks like code format, and unittests by leveraging the above mentioned
tools. If the check fails, you will not be able to merge the feature branch into master branch.

### Documentation Page Publishing
We use [github page](https://pages.github.com/) to host our documentation, and is separated as an [independent
repository](https://github.com/Music-and-Culture-Technology-Lab/omnizart-doc). 

**Please do not directly modify the content of the omnizart-doc repository!!**

The only permitted way to update the documentation page is by updating the `build_doc` branch, and
let the workflow do the rest of things.

Steps to update the documentation page:
* Clone **this** repo
* Create a new branch. **DO NOT UPDATE THE `build_doc` BRANCH DIRECTLY!!**
* File a pull request
* Merge into master (by admin)
* Merge into `build_doc` branch (by admin)
* Push to this repo (by admin)

### Publish PyPI Package and Docker Image
Publish the python package to PyPI and also the docker image to dockerhub when push tags to the repository.
The publish process will be automatically done by the github actions. There are several steps in the process:

1. Pack and publish the python package.
2. Build the docker image and publish to Docker Hub.
3. Create release -> this will also trigger the automation of documentation publishment.


## Docker
We provide both the Dockerfile for local image build and also the pre-build image on Docker Hub.

To build the image, run the following:
```
docker build -t omnizart:my-image .
```

To use the pre-build image, follow below steps:
```
# Download from Docker Hub
docker pull mctlab/omnizart

# Execute the image
docker run -it mctlab/omnizart:latest

### For those who want to leverage the power of GPU for acceleration, make sure
### you have installed docker>=19.03 and the 'nvidia-container-toolkit' package.
# Execute the docker with GPU support
docker run --gpus all -it mctlab/omnizart:latest
```


## Command Test
To actually install and test the `omnizart` command, execute `make install`. This will automatically create a virtual environment and install everything needed inside it. After installation, just follow the instruction showing on the screen to activate the environment, then type `omnizart --help` to check if it works. After testing the command, type `deactivate` to leave the virtual environment.

## Others
### Log Level
The default log level is set to `warn`. You can change it by exporting environment variable *LOG_LEVEL* to one of `debug`, `info`, `warning`, `error`, or `critical`. The verbosity is sorted from high to low (debug -> critical). For the consideration behind the log level design, please refer to the [soruce code](https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/master/omnizart/utils.py#L20) or the [documentation page](https://music-and-culture-technology-lab.github.io/omnizart-doc/utils.html#omnizart.utils.get_logger)
### NNLS Chroma ###

System identifier – vamp:nnls-chroma:nnls-chroma
RDF URI – http://vamp-plugins.org/rdf/plugins/nnls-chroma#nnls-chroma

#### General Description ####

NNLS Chroma analyses a single channel of audio using frame-wise spectral input from the Vamp host. The plugin was originally developed to extract treble and bass chromagrams for subsequent use in chord extraction methods. The spectrum is transformed to a log-frequency spectrum (constant-Q) with three bins per semitone. On this representation, two processing steps are performed:
* tuning, after which each centre bin (i.e. bin 2, 5, 8, ...) corresponds to a semitone, even if the tuning of the piece deviates from 440 Hz standard pitch.
* running standardisation: subtraction of the running mean, division by the running standard deviation. This has a spectral whitening effect.

The processed log-frequency spectrum is then used as an input for NNLS approximate transcription (using a dictionary of harmonic notes with geometrically decaying harmonics magnitudes). The output of the NNLS approximate transcription is semitone-spaced. To get the chroma, this semitone spectrum is multiplied (element-wise) with the desired profile (chroma or bass chroma) and then mapped to 12 bins. The resulting chroma frames can be normalised by (dividing by) their norm (L1, L2 and maximum norm available).

#### Parameters ####

The default settings (in brackets, below) are those used for Matthias Mauch's 2010 MIREX submissions.

* use approximate transcription (NNLS) (on or off; default: on): toggle between NNLS approximate transcription and linear spectral mapping.
* spectral roll on spectral roll on (0 % -- 5 %; default: 0 %): consider the cumulative energy spectrum (from low to high frequencies). All bins below the first bin whose cumulative energy exceeds the quantile [spectral roll on] x [total energy] will be set to 0. A value of 0 means that no bins will be changed.
* tuning mode (global or local; default: global): local uses a local average for tuning, global uses all audio frames. Local tuning is only advisable when the tuning is likely to change over the audio, for example in podcasts, or in a cappella singing.
* spectral whitening (0.0 -- 1.0; default: 1.0): determines how much the log-frequency spectrum is whitened. A value of 0.0 means no whitening. For values other than 0.0 the log-freq spectral bins are divided by  [standard deviation of their neighbours]^[spectral whitening], where "^" means "to the power of".
* spectral shape (0.5 -- 0.9; default: 0.7): the shape of the notes in the NNLS dictionary. Their harmonic amplitude follows a geometrically decreasing pattern, in which the i-th harmonic has an amplitude of [spectral shape]^[i-1], where "^" means "to the power of".
* chroma normalisation (none, maximum norm, L1 norm, L2 norm; default: none): determines whether or how the chromagrams are normalised. If the setting is not 'none', then each chroma frame separately is divided by the chosen vector norm. Note that normalisation implies that the joint 24-dim. "Chroma and Bass Chromagram" output will be different from the individual 12-dim. "Chromagram" and "Bass Chromagram" outputs.

#### Outputs ####

* Log-frequency Spectrum: a spectrum similar to the well-known constant Q spectrum, in which bins are linear in log-frequency. Three bins per semitone.
* Tuned Log-frequency Spectrum: has the same format as Log-frequency Spectrum, but has been processed by the following processes: tuning, subtraction of background spectrum, spectral whitening.
* Semitone Spectrum: a spectral representation with one bin per semitone. If NNLS is selected in the parameters, this is the note activation, otherwise just a linear mapping to semitones.
* Bass Chromagram: a 12-dimensional chromagram, restricted to the bass range. At each frame the Semitone Spectrum is multiplied by a bass pattern and then mapped to the 12 chroma bins. 
* Chromagram: a 12-dimensional chromagram, restricted with mid-range emphasis. At each frame the Semitone Spectrum is multiplied by a mid-range pattern and then mapped to the 12 chroma bins.
* Chromagram and Bass Chromagram: a 24-dimensional chromagram, consisting of the both Bass Chromgram and Chromagram, see above. When normalisation is used, this representation will however be scaled differently, and hence be different from the individual chromagrams.
* Consonance estimate: A simple consonance value based on the convolution of a consonance profile with the Semitone Spectrum. Experimental status. Compare two pieces of audio in terms of consonance if the instrumentation is similar. Instruments with fluctuating pitches (also: voice) will decrease the consonance value.

### Chordino ###

System identifier – vamp:nnls-chroma:chordino
RDF URI – http://vamp-plugins.org/rdf/plugins/nnls-chroma#chordino

#### General Description ####

Chordino provides a simple chord transcription based on NNLS Chroma (described above). Chord profiles given by the user in the file "chord.dict" are used to calculate frame-wise chord similarities. Two simple (non-state-of-the-art!) algorithms are available that smooth these to provide a chord transcription: a simple chord change method, and a standard HMM/Viterbi approach.

#### Parameters ####

* use approximate transcription (NNLS) (on or off; default: on): toggle between NNLS approximate transcription and linear spectral mapping.
* HMM (Viterbi decoding) (on or off; default: on): uses HMM/Viterbi smoothing. Otherwise: heuristic chord change smoothing.
* spectral roll on (0 % -- 5 %; default: 0 %): consider the cumulative energy spectrum (from low to high frequencies). All bins below the first bin whose cumulative energy exceeds the quantile [spectral roll on] x [total energy] will be set to 0. A value of 0 means that no bins will be changed.
* tuning mode (global or local; default: global): local uses a local average for tuning. Local tuning is only advisable when the tuning is likely to change over the audio, for example in podcasts, or in a cappella singing.
* spectral whitening (0.0 -- 1.0; default: 1.0): determines how much the log-frequency spectrum is whitened. A value of 0.0 means no whitening. For values other than 0.0 the log-freq spectral bins are divided by  [standard deviation of their neighbours]^[spectral whitening], where "^" means "to the power of".
* spectral shape (0.5 -- 0.9; default: 0.7): the shape of the notes in the NNLS dictionary. Their harmonic amplitude follows a geometrically decreasing pattern, in which the i-th harmonic has an amplitude of [spectral shape]^[i-1], where "^" means "to the power of".
* chroma normalisation (none, maximum norm, L1 norm, L2 norm; default: none): determines whether or how the chromagrams are normalised. If the setting is not 'none', then each chroma frame separately is divided by the chosen vector norm. Note that normalisation implies that the joint 24-dim. "Chroma and Bass Chromagram" output will be different from the individual 12-dim. "Chromagram" and "Bass Chromagram" outputs.
* boost likelihood of the N (no chord) label (0.0 -- 1.0; default: 0.1): leads to greater values in the profile of the "no chord" chord, hence non-harmonic parts of audio files are more likely to be recogised as such. Warning: for values above the default, it quickly leads to many chords being misclassified as N.

#### Outputs ####

* Chord Estimate: estimated chord times and labels.
* Harmonic Change Value: an indication of the likelihood of harmonic change. Depends on the chord dictionary. Calculation is different depending on whether the Viterbi algorithm is used for chord estimation, or the simple chord estimate.
* Note Representation of Chord Estimate: a simple MIDI-like represenation of the estimated chord with bass note (if applicable) and chord notes. Can be used, for example, to export MIDI chords from Sonic Visuliser.

### Tuning ###

System identifier – vamp:nnls-chroma:tuning
RDF URI – http://vamp-plugins.org/rdf/plugins/nnls-chroma#tuning

#### General Description ####

The tuning plugin can estimate the local and global tuning of piece. The same tuning method is used for the NNLS Chroma and Chordino plugins.

#### Parameter ####

* spectral roll on spectral roll on (0 % -- 5 %; default: 0 %): consider the cumulative energy spectrum (from low to high frequencies). All bins below the first bin whose cumulative energy exceeds the quantile [spectral roll on] x [total energy] will be set to 0. A value of 0 means that no bins will be changed.

#### Outputs ####

* Tuning: returns a single label (at time 0 seconds) containing an estimate of the concert pitch in Hz.
* Local Tuning: returns a tuning estimate at every analysis frame, an average of the (recent) previous frame-wise estimates of the concert pitch in Hz.

### References and Credits ###

If you make use of this software for any public or commercial purpose,
we ask you to kindly mention the authors and Queen Mary, University of
London in your user-visible documentation. We're very happy to see
this sort of use, but would much appreciate being credited, separately
from the requirements of the software license itself (see below).

If you make use of this software for academic purposes, please cite:

Mauch, Matthias and Dixon, Simon: [*Approximate Note Transcription for the Improved Identification of Difficult Chords*](http://schall-und-mauch.de/artificialmusicality/?p=89), Proceedings of the 11th International Society for Music Information Retrieval Conference (ISMIR 2010), 2010.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
Add New Modules
===============

This page describes how to add new modules into Omnizart project, adapt the original implementations
to omnizart's architecture.

Before starting walking through the integration process, be sure you have already read the
`CONTRIBUTING.md <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/master/CONTRIBUTING.md>`_,
and the `slides of omnziart <https://drive.google.com/file/d/1IO1lh07nMvSi0X0nzRDT7kuE1f468Rl1/view?usp=sharing>`_
could also be helpful for your understanding of this project.
Additionally, there are few more things to be always kept in mind while developing omnizart.

Principles
##########

* **Find an existing module and start development** - There are already several implemented modules
  that are fully functional, and being great examples that give you hints on your way developing
  new modules. Most of them are very similar of their overall architecture, but vary in detail.
  Most the time, you could just copy and paste the small pieces to your module, and modify just a
  small part of them to adapt to your task.
* **Try not to make your own wheels** - There have been many useful and validated functions that are
  developed to deal with the daily works. They are already there to cover 90% of every details of a
  module, thus new logics are in very small chances being needed.
  Most of the time you need to implement the most would be the part of feature and label extraction,
  which will be explained in the upcoming sections.
* **Check with linters frequently** - You should always do ``make lint`` before you push to github,
  checking that there aren't any errors with the code format, or the build process would fail.
* **Don't permit linter errors easily** - You may find some comments that permits the linter errors
  while surfing the code. Those are quick solutions while in the early development of omnizart, which
  saves lots of time fixing those lint errors. But it should not be the main concern now, as the
  architecture is more stable and less error prone. You should follow every hints by the linters
  and fix them before you file a pull request.


----

So now we are all set and ready to add a new module to omnizart. Here we will take the
`PR #11 <https://github.com/Music-and-Culture-Technology-Lab/omnizart/pull/11>`_ as the example.

Setup
#####

1. **IMPORTANT** - Give your module a short, yet descriptive name. In the example, the name is
   ``PatchCNN`` (camel-case), ``patch_cnn`` (snake-case), and ``patch-cnn`` (for CLI use).

2. Create a folder named after your module under ``omnizart``. There should be at least two files:
   ``app.py`` and ``__init__.py``.

Implement Feature Generation
############################

The process unit should be **a dataset**, means the function accepts the path to the dataset itself, and will handle the rest
of the things like dataset type inferring, folder structure handling, file parsing, feature extraction, and output storage.

Commits
*******

* `3ff6c4a <https://github.com/Music-and-Culture-Technology-Lab/omnizart/pull/11/commits/3ff6c4abe5ab98242d33c146353b5282ce5f6b66>`_
  - Builds the main structure of feature extraction.
* `f3138eb <https://github.com/Music-and-Culture-Technology-Lab/omnizart/pull/11/commits/f3138eb4a0650c91692f70e09bab1578be11c132>`_
  - Contains the patch of label extraction function.
* `0190f18 <https://github.com/Music-and-Culture-Technology-Lab/omnizart/pull/11/commits/0190f1895027cf859647c2099d3c03a24f73246a>`_
  - Contains the patch of label extraction function.

Critical Files/Functions
************************

* `omnizart.patch_cnn.app.PatchCNNTranscription.generate_feature <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L106-L193>`_
  - The main function for managing the process of feature generation.

* `omnizart.feature.cfp.extract_patch_cfp <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/feature/cfp.py#L355-L451>`_
  - The function for feature extraction, which takes audio path as the input and outputs the required feature representations.

* `omnizart.patch_cnn.app.extract_label <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L278-L327>`_
  - The function for label extraction, generating the representation of ground-truth. Accepts the path to the ground-truth file, parses the contents
  into intermediate format (see :class:`omnizart.base.Label`), and extracts necessary informations.

  Normally, it should be defined in a separate file called ``labels.py`` under *omnizart/<module>/* when the extraction process contains lots of logics.
  Since the label extraction in this example is relatively simple, it is okay to put it under ``omnizart.patch_cnn.app``.
  See the conventional case :class:`omnizart.drum.labels`.

* `omnizart.patch_cnn._parallel_feature_extraction <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L336-L373>`_

  - To boost the process of feature extraction, files are processed in parallel. You can use the function
  :class:`omnizart.utils.parallel_generator` to accelerate the process.

* `omnizart.setting_loaders.PatchCNNSettings <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/setting_loaders.py#L330-L357>`_
  - The data class that holds the necessary hyper parameters that will be used by different functions of this module. For feature extraction, the
  parameters are registered under ``PatchCNNSettings.feature`` attribute.

* `omnizart/defatuls/patch_cnn.yaml <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/checkpoints/patch_cnn/patch_cnn_melody/configurations.yaml#L13-L53>`_
  - The configuration file of the module, records the values of hyper parameters and will be consumed by the data class (i.e. PatchCNNSettings).

Overall Process Flow
********************

1. Determine the dataset type from the given dataset path.
    * `music module <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/master/omnizart/music/app.py#L169-L179>`_
2. Choose the corresponding dataset structure class.
    * `patch-cnn <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L135>`_
    * `music module <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/master/omnizart/music/app.py#L182-L186>`_
3. Parse audio/ground-truth file pairs.
    * `patch-cnn <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L163-L167>`_
4. Make sure feature output path exists.
    * `patch-cnn <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L169-L172>`_
5. Parallel generate feature and label representation.
    * `patch-cnn <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L174-L188>`_
6. Write the settings to the output path, named as *.success.yaml*.
    * `patch-cnn <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L190-L193>`_


Implement Model Training
########################

All the training should happen in the ``.fit()`` function to fine-tune the model. There is supposed no need to manually
write the training loop.

Commits
*******

* `2d6f74d <https://github.com/Music-and-Culture-Technology-Lab/omnizart/pull/11/commits/2d6f74da88e52cef7ef6e96f3b93be97771bdf31>`_

Critical Files/Functions
************************

* `omnizart.patch_cnn.app.PatchCNNTranscription.train <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L195-L275>`_
  - The main function for managing the training flow.

* `omnizart.models.patch_cnn <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/master/omnizart/models/patch_cnn.py>`_
  - Definition of the model. You can also customize the ``train_step`` function to do more sophisticated loss computation. See examples
  in `vocal <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/master/omnizart/models/pyramid_net.py#L233-L284>`_
  and `chord <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/master/omnizart/models/chord_model.py#L547-L600>`_
  modules.

* `omnizart.patch_cnn.app.PatchCNNDatasetLoader <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L376-L380>`_
  - The dataset loader for feeding data to models. Dealing with listing files, iterating through all feature/label pairs,
  indexing, and additionally augmenting, clipping, or transforming the feature/label on the fly.

* `omnizart.setting_loaders.PatchCNNSettings <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/setting_loaders.py#L366-L386>`_
  - The data class that holds the necessary hyper parameters that will be used by different functions of this module. For model training,
  related hyper parameters are registered under ``PatchCNNSettings.dataset``, ``PatchCNNSettings.model``, and
  ``PatchCNNSettings.training`` attributes.

* `omnizart/defatuls/patch_cnn.yaml (1) <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/checkpoints/patch_cnn/patch_cnn_melody/configurations.yaml#L54-L75>`_ /
  `omnizart/defatuls/patch_cnn.yaml (2) <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/checkpoints/patch_cnn/patch_cnn_melody/configurations.yaml#L88-L118>`_
  - The configuration file of the module, records the values of hyper parameters and will be consumed by the data class (i.e. PatchCNNSettings).

Overall Process Flow
********************

1. Check whether there is an input model or not. If given input model path, this indicating the user wants to fine-tune on a previously trained model. The coressponding settings should also be updated.
    * `patch-cnn <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L216-L219>`_
    * `drum <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/master/omnizart/drum/app.py#L167-L172>`_
2. Decide the portion of training and validation set.
    * `patch-cnn <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L221-L223>`_
    * `drum <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/master/omnizart/drum/app.py#L174-L176>`_
3. Construct dataset loader instances for training and validation.
    * `patch-cnn <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L225-L236>`_
4. Construct a fresh model if there is no input model.
    * `patch-cnn <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L238-L240>`_
5. Compile the model with loss function
    * `patch-cnn <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L242-L244>`_
    * `music <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/master/omnizart/drum/app.py#L167-L172>`_
6. Resolve the output path of the model
    * `patch-cnn <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L246-L255>`_
7. Construct the callbacks for storing the checkpoints, early stopping the training, and others.
    * `patch-cnn <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L257-L262>`_
8. Start training
    * `patch-cnn <https://github.com/Music-and-Culture-Technology-Lab/omnizart/blob/273fc60fbc6e3728c07abf71e06cf8f092bfabeb/omnizart/patch_cnn/app.py#L264-L274>`_


Implement Transcription
#######################

Add Unit Tests
##############

Commit Checkpoints
##################

Implement CLI
#############

Add Documentation
#################

----

Optional
########

This section holds the optional actions you can do, while it is not necessary to be done
during implementing a new module.

Add new supported datasets
**************************

If you want to add a new dataset that is currently not supported by ``omnizart`` (which is defined in
:class:`omnizart.constants.datasets`), things should be noticed are explained in this section.

(To be continued...)
Feature
=======


.. Introduction

.. automodule:: omnizart.feature
    :members:


CFP
###

.. automodule:: omnizart.feature.cfp
    :members:
    :undoc-members:


HCFP
####

.. automodule:: omnizart.feature.hcfp
    :members:
    :undoc-members:


CQT
###

.. automodule:: omnizart.feature.cqt
    :members:
    :undoc-members:


Beat Tracking
#############

.. automodule:: omnizart.feature.beat_for_drum
    :members:
    :undoc-members:
Utilities
=========

Some common utility functions.


Remote
######

.. automodule:: omnizart.remote
    :members:
    :undoc-members:


Utility Functions
#################

.. automodule:: omnizart.utils
    :members:
    :undoc-members:
Models
======

.. Introduction

.. automodule:: omnizart.models
    :members:


U-Net
#####

.. automodule:: omnizart.models.u_net
    :members:
    :undoc-members:
    :show-inheritance:


Tensor2Tensor
#############

.. automodule:: omnizart.models.t2t
    :members:
    :undoc-members:


Spectral Normalization Model
############################

.. automodule:: omnizart.models.spectral_norm_net
    :members:
    :undoc-members:
    :show-inheritance:


Chord Transformer
#################

.. automodule:: omnizart.models.chord_model
    :members:
    :undoc-members:
    :show-inheritance:


Pyramid Net
###########

.. automodule:: omnizart.models.pyramid_net
    :members:
    :undoc-members:
    :show-inheritance:


Utils
#####

.. automodule:: omnizart.models.utils
    :members:
    :undoc-members:
Base Classes
============

.. automodule:: omnizart.base


Transcription
-------------
.. autoclass:: omnizart.base.BaseTranscription
    :members:


Label
-----
.. autoclass:: omnizart.base.Label
    :members:


Dataset Loader
--------------
.. autoclass:: omnizart.base.BaseDatasetLoader
    :members:

Training
========

Including logics of main training loop, progress visualization, and callback
functions.


Training
########

.. autodata:: omnizart.train.PROGRESS_BAR_FORMAT
    :annotation: = Format of the training progress bar

.. automodule:: omnizart.train
    :members:
    :undoc-members:


Callbacks
#########
.. automodule:: omnizart.callbacks
    :members:
.. Documents are written in reStructured Text (.rst) format.
   Learn the syntax from: https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html

   Heading Level (most significant to least):
     Underline with '='
     Underline with '#'
     Underline with '*'


Tutorial
========

This page describes the workflow and usage of ``omnizart`` command-line interface,
covering core and utility.

The root entry is ``omnizart`` followed by sub-commands.
The available sub-commands can be found by typing ``omnizart --help``.

Core
####

In general, the core sub-commands follow a pipeline of ``application``-``action``-``arguments``:

.. code-block:: bash

   omnizart application action --arguments

where we apply an ``action`` to the ``application`` of interest, with corresponding ``arguments``.
Detailed descriptions for the usage of each sub-command can be found in the dedicated pages for each ``application``:

* :doc:`music/cli`
* :doc:`drum/cli`
* :doc:`chord/cli`
* :doc:`vocal-contour/cli`
* :doc:`vocal/cli`
* :doc:`beat/cli`

All the applications share a same set of actions: **transcribe**, **generate-feature**, and **train-model**.
Let's have a walkthrough of each ``action``.

Transcribe
**********

As the name suggests, this action transcribes a given input.
The supported applications are as follows:

* ``music`` - Transcribe musical notes of pitched instruments in MIDI.
* ``drum`` - Transcribe events of percussive instruments in MIDI.
* ``chord`` - Transcribe chord progressions in MIDI and CSV.
* ``vocal`` - Transcribe note-level vocal melody in MIDI.
* ``vocal-contour`` - Transcribe frame-level vocal melody (F0) in text.
* ``beat`` - Transcribe beat position.

Note that all the applications receive polyphonic music in WAV, except ``beat`` receives inputs in MIDI.

Example usage:

.. code-block:: bash

   # Transcribe percussive events given pop.wav, with specified model path and output directory
   omnizart drum transcribe pop.wav --model-path ./my-model --output ./trans_pop.mid

Note: ``--model-path`` can be left unspecified, and the default will be the downloaded checkpoints.
Execute ``omnizart download-checkpoints`` if you have not done in the installation from :doc:`quick-start`.


Generate Feature
****************

This action generates the features that are necessary for training and testing.
You can definitely skip this if you are only into transcribing with the given checkpoints.
The processed features will be stored in *<path/to/dataset>/train_feature* and *<path/to/dataset>/test_feature*.

The supported datasets for feature processing are application-dependent, summarized as follows:

+-------------+-------+------+-------+------+-------+---------------+------+
| Module      | music | drum | chord | beat | vocal | vocal-contour | beat |
+=============+=======+======+=======+======+=======+===============+======+
| Maestro     |   O   |      |       |      |       |               |      |
+-------------+-------+------+-------+------+-------+---------------+------+
| Maps        |   O   |      |       |      |       |               |      |
+-------------+-------+------+-------+------+-------+---------------+------+
| MusicNet    |   O   |      |       |      |       |               |  O   |
+-------------+-------+------+-------+------+-------+---------------+------+
| Pop         |   O   |  O   |       |      |       |               |      |
+-------------+-------+------+-------+------+-------+---------------+------+
| Ext-Su      |   O   |      |       |      |       |               |      |
+-------------+-------+------+-------+------+-------+---------------+------+
| BillBoard   |       |      |   O   |      |       |               |      |
+-------------+-------+------+-------+------+-------+---------------+------+
| BPS-FH      |       |      |       |      |       |               |      |
+-------------+-------+------+-------+------+-------+---------------+------+
| MIR-1K      |       |      |       |      |   O   |       O       |      |
+-------------+-------+------+-------+------+-------+---------------+------+
| MedleyDB    |       |      |       |      |       |       O       |      |
+-------------+-------+------+-------+------+-------+---------------+------+
| Tonas       |       |      |       |      |   O   |               |      |
+-------------+-------+------+-------+------+-------+---------------+------+

Before running the commands below, make sure to download the corresponding datasets first.
This can be easily done in :ref:`Download Datasets`.

.. code-block:: bash

   # Generate features for the music application
   omnizart music generate-feature --dataset-path <path/to/dataset>

   # Generate features for the drum application
   omnizart drum generate-feature --dataset-path <path/to/dataset>


Train Model
***********

This action trains a model from scratch given the generated features from :ref:`Generate Feature`.
Once again, you can skip this if you are only up to transcribing music, and use the provided checkpoints.

.. code-block:: bash

   omnizart music train-model -d <path/to/feature/folder> --model-name My-Music
   omnizart drum train-model -d <path/to/feature/folder> --model-name My-Drum
   omnizart chord train-model -d <path/to/feature/folder> --model-name My-Chord


Utility
#######


Download Datasets
*****************

This sub-command belongs to the utility, used to download the datasets for training and testing the models.
Current supported datasets are:

* `Maestro <https://magenta.tensorflow.org/datasets/maestro>`_ - MIDI and Audio Edited for Synchronous TRacks and Organization dataset.
* `MusicNet <https://homes.cs.washington.edu/~thickstn/musicnet.html>`_ - MusicNet dataset with a collection of 330 freely-licensed classical music recordings.
* `McGill <https://ddmal.music.mcgill.ca/research/The_McGill_Billboard_Project_(Chord_Analysis_Dataset)/>`_ - McGill BillBoard dataset.
* `BPS-FH <https://github.com/Tsung-Ping/functional-harmony>`_ - Beethoven Piano Sonata with Function Harmony dataset.
* Ext-Su - Extended Su dataset.
* `MIR-1K <https://sites.google.com/site/unvoicedsoundseparation/mir-1k>`_ - 1000 short clips of Mandarin pop songs.
* `MedleyDB <http://medleydb.weebly.com/>`_ - 122 multitracks.

Example usage:

.. code-block:: bash

   # Download the MAESTRO dataset and output to the */data* folder.
   omnizart download-dataset Maestro --output /data

   # Download the MusicNet dataset and unzip the dataset after download.
   omnizart download-dataset MusicNet --unzip

   # To see a complete list of available datasets, execute the following command
   omnizart download-dataset --help


Download Checkpoints
********************

This is the other sub-command for the utility, used to download the archived checkpoints of pre-trained models.

.. code-block:: bash

   # Simply run the following command, and no other options are needed to be specified.
   omnizart download-checkpoints
.. omnizart documentation master file, created by
   sphinx-quickstart on Tue Aug 25 10:43:56 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


OMNIZART: MUSIC TRANSCRIPTION MADE EASY
=======================================

.. figure:: ../../figures/features2.png
   :align: center


Omnizart is a Python library and a streamlined solution for automatic music transcription.
This library gathers the research outcomes from `Music and Cultural Technology Lab <https://sites.google.com/view/mctl/home>`_, 
analyzing polyphonic music and transcribes 
**musical notes of instruments** :cite:`music`,
**chord progression** :cite:`chord`,
**drum events** :cite:`drum`,
**frame-level vocal melody** :cite:`vocalcontour`,
**note-level vocal melody**  :cite:`vocal`, and
**beat** :cite:`beat`.

Omnizart provides the main functionalities that construct the life-cycle of deep learning-based music transcription,
covering from *dataset downloading*, *feature pre-processing*, *model training*, to *transcription* and *sonification*.
Pre-trained checkpoints are also provided for the immediate usage of transcription. The paper can be found from
`Journal of Open Source Software (JOSS) <https://doi.org/10.21105/joss.03391>`_.

Demonstration
#############

Colab
*****

Play with the `Colab notebook <https://bit.ly/OmnizartColab>`_ to transcribe your favorite song almost immediately!

Replicate web demo
******************

Transcribe music with `Replicate web UI <https://replicate.com/breezewhite/omnizart>`_.

Sound samples
*************

Original song

.. raw:: html

   <iframe src="https://www.youtube-nocookie.com/embed/hjJhweRlE-A" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

Chord transcription

.. raw:: html

   <audio controls="controls">
      <source src="_audio/high_chord_synth.mp3" type="audio/mpeg">
      Your browser does not support the <code>audio</code> element.
   </audio>


Drum transcription

.. raw:: html

   <audio controls="controls">
      <source src="_audio/high_drum_synth.mp3" type="audio/mpeg">
      Your browser does not support the <code>audio</code> element.
   </audio>


Note-level vocal transcription

.. raw:: html

   <audio controls="controls">
      <source src="_audio/high_vocal_synth.mp3" type="audio/mpeg">
      Your browser does not support the <code>audio</code> element.
   </audio>


Frame-level vocal transcription

.. raw:: html

   <audio controls="controls">
      <source src="_audio/high_vocal_contour.mp3" type="audio/mpeg">
      Your browser does not support the <code>audio</code> element.
   </audio>


Source files can be downloaded `here <https://drive.google.com/file/d/15VqHearznV9L83cyl61ccACsXXJ4vBHo/view?usp=sharing>`_.
You can use *Audacity* to open the files.


.. toctree::
   :maxdepth: 2
   :caption: Contents

   quick-start.rst
   tutorial.rst
   demo.rst

.. toctree::
   :maxdepth: 2
   :caption: Command Line Interface

   music/cli.rst
   drum/cli.rst
   chord/cli.rst
   vocal/cli.rst
   vocal-contour/cli.rst
   beat/cli.rst
   patch-cnn/cli.rst


.. toctree::
   :maxdepth: 2
   :caption: API Reference

   music/api.rst
   drum/api.rst
   chord/api.rst
   vocal/api.rst
   vocal-contour/api.rst
   patch-cnn/api.rst
   beat/api.rst
   feature.rst
   models.rst
   training.rst
   base.rst
   constants.rst
   utils.rst

.. Indices and tables
  ==================
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`


References
##########

.. bibliography::
   refs.bib
Constants
=========

.. Introduction

.. automodule:: omnizart.constants
    :members:
    :undoc-members:


Feature
#######

Records the settings of the feature extraction process.
There are also default settings that can be adjusted, and are records in `defaults/*.yaml` files.

.. automodule:: omnizart.constants.feature
    :members:
    :undoc-members:


Datasets
########

Records the directory structure of each dataset, and will be used for extracting the feature
of the whole dataset (not yet supported...).

.. automodule:: omnizart.constants.datasets
    :members:
    :show-inheritance:


Midi
####

Records MIDI related settings, including mapping of program number and the corresponding 
instrument name.

.. automodule:: omnizart.constants.midi
    :members:
    :undoc-members:
Demonstration
=============


葬予規路火烌猶在 - 柯拉琪 Collage
-------------------------------

.. raw:: html

    <iframe src="_static/demo-collage.html" style="border: none; height: 500px;"></iframe>


Last Stardust Piano Ver. - Animenz
---------------------------------

.. raw:: html

    <iframe src="_static/demo-stardust.html" style="border: none; height: 500px;"></iframe>


Self Spiral - The Surrealist
----------------------------
A pure instrumental piece.

.. raw:: html

    <iframe src="_static/demo-surrealist.html" style="border: none; height: 500px;"></iframe>


Question Everything - Dreamshade
--------------------------------

.. raw:: html

    <iframe src="_static/demo-dreamshade.html" style="border: none; height: 700px;"></iframe>


Cosmo Funk - Snail's House
--------------------------

.. raw:: html

    <iframe src="_static/demo-cosmo.html" style="border: none; height: 800px;"></iframe>
Quick Start
===========

Colab
#####

Play with the `Colab notebook <https://bit.ly/omnizart-colab>`_  to transcribe your favorite song without hassles.
You can also follow the installation below to enjoy Omnizart locally.

Installation
############

Using pip
*********

Omnizart is under development and will be updated regularly on PyPI.
Use ``pip`` to install the latest stable version.

.. code-block:: bash

    # Install the prerequisites manually since there are some dependencies can't be
    # resolved automatically.
    pip install numpy Cython

    # Additional system packages are required to fully use Omnizart.
    sudo apt-get install libsndfile-dev fluidsynth ffmpeg

    # Install Omnizart
    pip install omnizart

    # Then download the checkpoints
    omnizart download-checkpoints


Development installation
************************

For the development installation, clone the git repo and the installation
creates a virtual environment under the directory *omnizart/* by default.

.. code-block:: bash

    # Clone the omnizart repository from GitHub
    git clone https://github.com/Music-and-Culture-Technology-Lab/omnizart.git

    # Install dependencies, with checkpoints automatically downloaded
    cd omnizart
    make install

    # Install Dev dependencies, since they will not be installed by default
    poetry install


CLI
###

Below is an example usage of pitched instrument transcription with the command-line interface,
first transcribing a piece of music and then synthesizing the results.
For more details and other types of transcription, refer to :doc:`tutorial`.

Transcription
*************

The example transcribes a piece of music, being monophonic or polyphonic,
to a MIDI file of the transcribed pitched notes and a CSV file with more information.

.. code-block:: bash

    omnizart music transcribe <path/to/example.wav>


Sonification
************

Omnizart renders the transcribed MIDI file with default soundfonts,
synthesizing an audio in WAV by the command below.
For the first-time execution, it is expected to take a bit for downloading the free-licensed soundfonts.

.. code-block:: bash

    omnizart synth example.mid


Auto Completion
***************

To enable auto completion, type the following according to your environment type.

.. code-block:: bash

    # For bash
    _OMNIZART_COMPLETE=source_bash omnizart > omnizart-complete.sh

    # For zsh
    _OMNIZART_COMPLETE=source_zsh omnizart > omnizart-complete.sh

    # Source the generated script to enable
    source omnizart-complete.sh
omnizart vocal-contour
======================

Lists the available options of each sub-command.


transcribe
##########

.. click:: omnizart.cli.vocal_contour.transcribe:transcribe
    :prog: omnizart vocal-contour transcribe


generate-feature
################

.. click:: omnizart.cli.vocal_contour.generate_feature:generate_feature
    :prog: omnizart vocal_contour generate-feature


train-model
###########

.. click:: omnizart.cli.vocal_contour.train_model:train_model
    :prog: omnizart vocal_contour train-model
Vocal-Contour Transcription
===========================


.. automodule:: omnizart.vocal_contour


App
###
.. automodule:: omnizart.vocal_contour.app
    :members:
    :show-inheritance:


Inference
#########
.. automodule:: omnizart.vocal_contour.inference
    :members:


Loss Functions
##############
.. automodule:: omnizart.music.losses
    :members:


Settings
########
Below are the default settings for frame-level vocal transcription. 
It will be loaded by the class :class:`omnizart.setting_loaders.VocalContourSettings`. 
The name of the attributes will be converted to snake-case (e.g. HopSize -> hop_size). 
There is also a path transformation when applying the settings into the ``VocalContourSettings`` instance. 
For example, the attribute ``BatchSize`` defined in the yaml path *General/Training/Settings/BatchSize* is transformed 
to *VocalContourSettings.training.batch_size*. 
The level of */Settings* is removed among all fields.

.. literalinclude:: ../../../omnizart/defaults/vocal_contour.yaml
    :language: yaml
omnizart beat
=============

Lists the detailed available options of each sub-commands.


transcribe
##########

.. click:: omnizart.cli.beat.transcribe:transcribe
    :prog: omnizart beat transcribe


generate-feature
################

.. click:: omnizart.cli.beat.generate_feature:generate_feature
    :prog: omnizart beat generate-feature


train-model
###########

.. click:: omnizart.cli.beat.train_model:train_model
    :prog: omnizart beat train-model
Beat Transcription
===================


.. automodule:: omnizart.beat


App
###
.. autoclass:: omnizart.beat.app.BeatTranscription
    :members:
    :show-inheritance:


Dataset
#######
.. autoclass:: omnizart.beat.app.BeatDatasetLoader
    :members:
    :show-inheritance:


Inference
#########
.. automodule:: omnizart.beat.inference
    :members:


Loss Functions
##############
.. autofunction:: omnizart.beat.app.weighted_binary_crossentropy


Features
########
.. automodule:: omnizart.beat.features
    :members:


Prediction
##########
.. automodule:: omnizart.beat.prediction
    :members:


Settings
########
Below are the default settings for building the beat model. It will be loaded
by the class :class:`omnizart.setting_loaders.BeatSettings`. The name of the
attributes will be converted to snake-case (e.g., HopSize -> hop_size). There
is also a path transformation process when applying the settings into the
``BeatSettings`` instance. For example, if you want to access the attribute
``BatchSize`` defined in the yaml path *General/Training/Settings/BatchSize*,
the corresponding attribute will be *BeatSettings.training.batch_size*.
The level of */Settings* is removed among all fields.

.. literalinclude:: ../../../omnizart/defaults/beat.yaml
    :language: yaml
omnizart chord
==============

Lists the detailed available options of each sub-commands.


transcribe
##########

.. click:: omnizart.cli.chord.transcribe:transcribe
    :prog: omnizart chord transcribe


generate-feature
################

.. click:: omnizart.cli.chord.generate_feature:generate_feature
    :prog: omnizart chord generate-feature


train-model
###########

.. click:: omnizart.cli.chord.train_model:train_model
    :prog: omnizart chord train-model
Chord Transcription
===================


.. automodule:: omnizart.chord


App
###
.. autoclass:: omnizart.chord.app.ChordTranscription
    :members:
    :show-inheritance:


Feature
#######
.. automodule:: omnizart.chord.features
    :members:
    :undoc-members:


Dataset
#######
.. autoclass:: omnizart.chord.app.McGillDatasetLoader
    :members:
    :show-inheritance:


Inference
#########
.. automodule:: omnizart.chord.inference
    :members:
    :undoc-members:


Settings
########
Below are the default settings for building the chord model. It will be loaded
by the class :class:`omnizart.setting_loaders.ChordSettings`. The name of the
attributes will be converted to snake-case (e.g., HopSize -> hop_size). There
is also a path transformation process when applying the settings into the
``ChordSettings`` instance. For example, if you want to access the attribute
``BatchSize`` defined in the yaml path *General/Training/Settings/BatchSize*,
the corresponding attribute will be *ChordSettings.training.batch_size*.
The level of */Settings* is removed among all fields.

.. literalinclude:: ../../../omnizart/defaults/chord.yaml
    :language: yaml
omnizart drum
=============

Lists the detailed available options of each sub-commands.


transcribe
##########

.. click:: omnizart.cli.drum.transcribe:transcribe
    :prog: omnizart drum transcribe


generate-feature
################

.. click:: omnizart.cli.drum.generate_feature:generate_feature
    :prog: omnizart drum generate-feature


train-model
###########

.. click:: omnizart.cli.drum.train_model:train_model
    :prog: omnizart drum train-model
Drum Transcription
==================


.. automodule:: omnizart.drum


App
###
.. autoclass:: omnizart.drum.app.DrumTranscription
    :members:
    :show-inheritance:


Dataset
#######
.. autoclass:: omnizart.drum.app.PopDatasetLoader
    :members:
    :show-inheritance:


Inference
#########
.. automodule:: omnizart.drum.inference
    :members:
    :undoc-members:


Labels
######
.. automodule:: omnizart.drum.labels
    :members:
    :undoc-members:


Prediction
##########
.. automodule:: omnizart.drum.prediction
    :members:
    :undoc-members:


Settings
########
Below are the default settings for building the drum model. It will be loaded
by the class :class:`omnizart.setting_loaders.DrumSettings`. The name of the
attributes will be converted to snake-case (e.g., HopSize -> hop_size). There
is also a path transformation process when applying the settings into the
``DrumSettings`` instance. For example, if you want to access the attribute
``BatchSize`` defined in the yaml path *General/Training/Settings/BatchSize*,
the corresponding attribute will be *DrumSettings.training.batch_size*.
The level of */Settings* is removed among all fields.

.. literalinclude:: ../../../omnizart/defaults/drum.yaml
    :language: yaml
omnizart music
==============

Lists the detailed available options of each sub-commands.


transcribe
##########

.. click:: omnizart.cli.music.transcribe:transcribe
    :prog: omnizart music transcribe


generate-feature
################

.. click:: omnizart.cli.music.generate_feature:generate_feature
    :prog: omnizart music generate-feature


train-model
###########

.. click:: omnizart.cli.music.train_model:train_model
    :prog: omnizart music train-model
Music Transcription
===================


.. automodule:: omnizart.music


App
###
.. autoclass:: omnizart.music.app.MusicTranscription
    :members:
    :show-inheritance:


Dataset
#######
.. autoclass:: omnizart.music.app.MusicDatasetLoader
    :members:
    :show-inheritance:


Inference
#########
.. automodule:: omnizart.music.inference
    :members:


Loss Functions
##############
.. automodule:: omnizart.music.losses
    :members:


Labels
######
.. automodule:: omnizart.music.labels
    :members:
    :undoc-members:


Prediction
##########
.. automodule:: omnizart.music.prediction
    :members:


Settings
########
Below are the default settings for building the music model. It will be loaded
by the class :class:`omnizart.setting_loaders.MusicSettings`. The name of the
attributes will be converted to snake-case (e.g., HopSize -> hop_size). There
is also a path transformation process when applying the settings into the
``MusicSettings`` instance. For example, if you want to access the attribute
``BatchSize`` defined in the yaml path *General/Training/Settings/BatchSize*,
the corresponding attribute will be *MusicSettings.training.batch_size*.
The level of */Settings* is removed among all fields.

.. literalinclude:: ../../../omnizart/defaults/music.yaml
    :language: yaml
omnizart patch-cnn
==================

Lists the detailed available options of each sub-commands.


transcribe
##########
.. click:: omnizart.cli.patch_cnn.transcribe:transcribe
    :prog: omnizart patch-cnn transcribe



generate-feature
################
.. click:: omnizart.cli.patch_cnn.generate_feature:generate_feature
    :prog: omnizart patch-cnn generate-feature


train-model
###########
.. click:: omnizart.cli.patch_cnn.train_model:train_model
    :prog: omnizart patch-cnn train-model
Patch-CNN Transcription
=======================


.. automodule:: omnizart.patch_cnn


App
###
.. autoclass:: omnizart.patch_cnn.app.PatchCNNTranscription
    :members:
    :show-inheritance:


Dataset
#######
.. autoclass:: omnizart.patch_cnn.app.PatchCNNDatasetLoader
    :members:
    :show-inheritance:


Inference
#########
.. automodule:: omnizart.patch_cnn.inference
    :members:


Labels
######
.. autofunction:: omnizart.patch_cnn.app.extract_label



Settings
########
Below are the default settings for building the PatchCNN model. It will be loaded
by the class :class:`omnizart.setting_loaders.PatchCNNSettings`. The name of the
attributes will be converted to snake-case (e.g. HopSize -> hop_size). There
is also a path transformation process when applying the settings into the
``PatchCNNSettings`` instance. For example, if you want to access the attribute
``BatchSize`` defined in the yaml path *General/Training/Settings/BatchSize*,
the coressponding attribute will be *MusicSettings.training.batch_size*.
The level of */Settings* is removed among all fields.

.. literalinclude:: ../../../omnizart/defaults/patch_cnn.yaml
    :language: yaml
omnizart vocal
==============

Lists the detailed available options of each sub-commands.


transcribe
##########

.. click:: omnizart.cli.vocal.transcribe:transcribe
    :prog: omnizart vocal transcribe


generate-feature
################

.. click:: omnizart.cli.vocal.generate_feature:generate_feature
    :prog: omnizart vocal generate-feature


train-model
###########

.. click:: omnizart.cli.vocal.train_model:train_model
    :prog: omnizart vocal train-model
Vocal Transcription
===================


.. automodule:: omnizart.vocal


App
###
.. autoclass:: omnizart.vocal.app.VocalTranscription
    :members:
    :show-inheritance:


Dataset
#######
.. autoclass:: omnizart.vocal.app.VocalDatasetLoader
    :members:
    :show-inheritance:


Inference
#########
.. automodule:: omnizart.vocal.inference
    :members:


Labels
######
.. automodule:: omnizart.vocal.labels
    :members:
    :undoc-members:


Prediction
##########
.. automodule:: omnizart.vocal.prediction
    :members:
    :undoc-members:


Settings
########
Below are the default settings for building the vocal model. It will be loaded
by the class :class:`omnizart.setting_loaders.VocalSettings`. The name of the
attributes will be converted to snake-case (e.g. HopSize -> hop_size). There
is also a path transformation process when applying the settings into the
``VocalSettings`` instance. For example, if you want to access the attribute
``BatchSize`` defined in the yaml path *General/Training/Settings/BatchSize*,
the coressponding attribute will be *VocalSettings.training.batch_size*.
The level of */Settings* is removed among all fields.

.. literalinclude:: ../../../omnizart/defaults/vocal.yaml
    :language: yaml
