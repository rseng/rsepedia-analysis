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
