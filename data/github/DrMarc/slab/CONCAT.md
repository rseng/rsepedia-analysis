![Package](https://github.com/DrMarc/slab/workflows/Python%20package/badge.svg)
[![PyPI](https://github.com/DrMarc/slab/workflows/PyPi/badge.svg)](https://pypi.org/project/slab/)
[![Documentation Status](https://readthedocs.org/projects/slab/badge/?version=latest)](https://slab.readthedocs.io/en/latest/?badge=latest)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-brightgreen.svg)](https://github.com/DrMarc/slab/graphs/commit-activity)
![PyPI pyversions](https://img.shields.io/badge/python-%3E%3D3.6-blue)
![PyPI license](https://img.shields.io/badge/license-MIT-brightgreen)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03284/status.svg)](https://doi.org/10.21105/joss.03284)

**slab**: easy manipulation of sounds and psychoacoustic experiments in Python
======================

**Slab** ('es-lab', or sound laboratory) is an open source project and Python package that makes working with sounds and running psychoacoustic experiments simple, efficient, and fun! For instance, it takes just eight lines of code to run a pure tone audiogram using an adaptive staircase:
```python
import slab
stimulus = slab.Sound.tone(frequency=500, duration=0.5) # make a 0.5 sec pure tone of 500 Hz
stairs = slab.Staircase(start_val=50, n_reversals=10) # set up the adaptive staircase
for level in stairs: # the staircase object returns a value between 0 and 50 dB for each trial
    stimulus.level = level
    stairs.present_tone_trial(stimulus) # plays the tone and records a keypress (1 for 'heard', 2 for 'not heard')
print(stairs.threshold()) # print threshold when done
```

Why slab?
---------
The package aims to lower the entrance barrier for working with sounds in Python and provide easy access to typical operations in psychoacoustics, specifically for students and researchers in the life sciences. The typical BSc or MSc student entering our lab has limited programming and signal processing training and is unable to implement a psychoacoustic experiment from scratch within the time limit of a BSc or MSc thesis. Slab solves this issue by providing easy-to-use building blocks for such experiments. The implementation is well documented and sufficiently simple for curious students to understand. All functions provide sensible defaults and will many cases 'just work' without arguments (vowel = slab.Sound.vowel() gives you a 1-second synthetic vowel 'a', vowel.spectrogram() plots the spectrogram). This turned out to be useful for teaching and demonstrations. Many students in our lab have now used the package to implement their final projects and exit the lab as proficient Python programmers.

Features
--------
Slab represents sounds as [Numpy](https://www.numpy.org) arrays and provides classes and methods to perform typical sound manipulation tasks and psychoacoustic procedures. The main classes are:

**Signal**: Provides a generic signal object with properties duration, number of samples, sample times, number of channels. Keeps the data in a 'data' property and implements slicing, arithmetic operations, and conversion between sample points and time points.
```python
sig = slab.Sound.pinknoise(n_channels=2) # make a 2-channel pink noise
sig.duration
# 1.0
sig.n_samples
# 8000
sig2 = sig.resample(samplerate=4000) # resample to 4 kHz
env = sig2.envelope() # returns a new signal containing the lowpass Hilbert envelopes of both channels
sig.delay(duration=0.0006, channel=0) # delay the first channel by 0.6 ms
```

**Sound**: Inherits from Signal and provides methods for generating, manipulating, displaying, and analysing sound stimuli. Can compute descriptive sound features and apply manipulations to all sounds in a folder.<sup id="a1">[1](#f1)</sup>
```python
vowel = slab.Sound.vowel(vowel='a', duration=.5) # make a 0.5-second synthetic vowel sound
vowel.play() # play the sound
vowel = vowel.ramp() # apply default raised-cosine onset and offset ramps
vowel = vowel.filter(kind='bp', frequency=[50, 3000]) # apply bandpass filter between 50 and 3000 Hz
vowel.spectrogram() # plot the spectrogram
vowel.spectrum(low_cutoff=100, high_cutoff=4000, log_power=True) # plot a band-limited spectrum
vowel.waveform(start=0, end=.1) # plot the waveform
vowel.write('vowel.wav') # save the sound to a WAV file
vocoded_vowel = vowel.vocode() # run a vocoding algorithm
vocoded_vowel.play() # play the vocoded sound
vowel.spectral_feature(feature='centroid') # compute the spectral centroid of the sound in Hz
# [1016.811]
```

**Binaural**: Inherits from Sound and provides methods for generating and manipulating binaural sounds, including advanced interaural time and intensity manipulation. Binaural sounds have left and a right channel properties.
```python
sig = slab.Binaural.pinknoise()
sig = sig.pulse() # make a 2-channel pulsed pink noise
sig.n_channels
# 2
right_lateralized = sig.itd(duration=600e-6) # add an interaural time difference of 600 µsec, right channel leading
# apply a linearly increasing or decreasing interaural time difference.
# This is achieved by sinc interpolation of one channel with a dynamic delay:
moving = sig.itd_ramp(from_itd=-0.001, to_itd=0.01)
level_spectrum = slab.Binaural.make_interaural_level_spectrum(hrtf) # compute frequency-band-specific ILDs from KEMAR
lateralized = sig.at_azimuth(azimuth=-45, ils=level_spectrum) # add frequency-dependent ITD and ILD corresponding to a sound at 45 deg
external = lateralized.externalize() # add an under-sampled HRTF filter that results in the percept of an external source
# (i.e. outside of the head), defaults to the KEMAR HRTF recordings, but any HRTF can be supplied
```

**Filter**: Inherits from Signal and provides methods for generating, measuring, and manipulating FIR and FFT filters, filter banks, and transfer functions.
```python
sig = slab.Sound.whitenoise()
filt = slab.Filter.band(frequency=2000, kind='hp') # make a highpass filter
filt.tf() # plot the transfer function
sig_filt = filt.apply(sig) # apply it to a sound
# applying a whole filterbank is equally easy:
fbank = slab.Filter.cos_filterbank(length=sig.n_samples, bandwidth=1/10, low_cutoff=100) # make a cosine filter bank
fbank.tf() # plot the transfer function of all filters in the bank
subbands = fbank.apply(sig) # make a multi-channel sound containing the passbands of the filters in the filter bank
subbands.spectrum(low_cutoff=90) # each band is limited by the corresponding fbank filter
# the subbands could now be manipulated and then combined with the collapse_subbands method
fbank.filter_bank_center_freqs() # return the centre frequencies of the filters in the filter bank
fbank = slab.Filter.equalizing_filterbank(reference, measured) # generate inverse filters to minimize the difference
# between measured signals and a reference sound. Used to equalize loudspeakers, microphones, or speaker arrays.
# measured is typically a recorded signal (potentially multi-channel), and reference for instance a flat white noise.
fbank.save('equalizing_filters.npy') # saves the filter bank as .npy file.
```

**HRTF**: Inherits from Filter, reads .sofa format HRTFs and provides methods for manipulating, plotting, and applying head-related transfer functions.
```python
hrtf = slab.HRTF.kemar() # load in-built KEMAR HRTF
print(hrtf) # print information
# <class 'hrtf.HRTF'> sources 710, elevations 14, samples 710, samplerate 44100.0
sourceidx = hrtf.cone_sources(20) # select sources on a cone of confusion at 20 deg from midline
hrtf.plot_sources(sourceidx) # plot the sources in 3D, highlighting the selected sources
hrtf.plot_tf(sourceidx,ear='left') # plot transfer functions of selected sources in a waterfall plot
dtf = hrtf.diffuse_field_equalization() # apply diffuse field equalization to remove non-spatial components of the HRTF
```

**Psychoacoustics**: A collection of classes for working trial sequences, adaptive staircases, forced-choice procedures, stimulus presentation and response recording from the keyboard and USB button boxes, handling of precomputed stimulus lists, results files, and experiment configuration files.
```python
# set up an 1up-2down adaptive weighted staircase with dynamic step sizes:
stairs = slab.Staircase(start_val=30, max_val=40, n_up=1, n_down=2,
                            step_sizes=[3, 1], step_up_factor=1.5)
for trial in stairs: # draw a value from the staircase; the loop terminates with the staircase
    response = stairs.simulate_response(25) # simulate a response from a participant using a psychometric function
    print(f'trial # {stairs.this_trial_n}: intensity {trial}, response {response}')
    stairs.add_response(response) # logs the response and advances the staircase
    stairs.plot() # updates a plot of the staircase in each trial to keep an eye on the performance of the listener
stairs.reversal_intensities # returns a list of stimulus values at the reversal points of the staircase
stairs.threshold() # computes and returns the final threshold
stairs.save_json('stairs.json') # the staircase object can be saved as a human readable json file

# for non-adaptive experiments and all other cases where you need a controlled sequence of stimulus values:
trials = slab.Trialsequence(conditions=5, n_reps=2) # sequence of 5 conditions, repeated twice, without direct repetitions
trials = slab.Trialsequence(conditions=['red', 'green', 'blue'], kind='infinite') # infinite sequence of color names
trials = slab.Trialsequence(conditions=3, n_reps=20, deviant_freq=0.12) # stimulus sequence for an oddball design
trials.transitions() # return the array of transition probabilities between all combinations of conditions.
trials.condition_probabilities() # return a list of frequencies of conditions
for trial in trials: # use the trials object in a loop to go through the trials
    print(trial) # here you would generate or select a stimulus according to the condition
    trials.present_afc_trial(target, distractor, isi=0.2) # present a 2-alternative forced-choice trial and record the response

stims = slab.Precomputed(lambda: slab.Sound.pinknoise(), n=10) # make 10 instances of noise as one Sound-like object
stims = slab.Precomputed([stim1, stim2, stim3, stim4, stim5]) # or use a list of sound objects, or a list comprehension
stims.play() # play a random instance
stims.play() # play another one, guaranteed to be different from the previous one
stims.sequence # the sequence of instances played so far
stims.write('stims.zip') # save the sounds as zip file of wavs
stims = slab.Precomputed.read('stims.zip') # reloads the file into a Precomputed object
```

<b id="f1">1)</b> The basic functionality of the Signal class and some of the sound generation methods in the Sound class were based on the brian.hears Sound class (now [brain2hears](https://brian2hears.readthedocs.io/en/stable/), an auditory modelling package). [↩](#a1)

Installation
------------

Install the current stable release from the python package index with pip:
```pip install slab```

### Other requirements ###

On *Linux*, there is only one requirement outside of Python: you may need to install *libsndfile* using your distribution’s package manager, for instance:

```sudo apt-get install libsndfile1```

On Macs with M1 processors, the SoundCard module that slab uses to play and record sounds is currently not working. You can workaround this issue by uninstalling SoundCard:

```pip uninstall soundcard```

Slab will fall back to `afplay` to play sounds. Recording sounds directly from slab is not possible in this case.

Other optional requirements can be installed by telling pip which extras you want:

```pip install slab[name_of_extra]```

The options for `name_of_extra` are:
- `windows`: if you are running Windows - this will install windows-curses for you, which is required for getting button presses in the psychoacoustics classes,
- `hrtf`: if you want to use spatial stimuli with the `Binaural` and `HRTF` classes,
- `testing`: (for developers) if you want to run the unit tests for slab, and
- `docs`: (for developers) if you want to build the documentation locally.

You can combine these options: `pip install slab[windows, hrtf]` if you are on Windows and use spatial sounds.

Detailed installation instructions can be found [here](https://slab.readthedocs.io/en/latest/index.html#installation).

You can also get the latest development version directly from GitHub (if you have [git](https://git-scm.com)) by running:
```pip install git+https://github.com/DrMarc/slab.git```

The releases use [semantic versioning](https://semver.org): major.minor.patch, where `major` increments for changes that break backwards compatibility, `minor` increments for added functionality, and `patch` increments for internal bug fixes.
```slab.__version__``` prints the installed version.

Documentation
-------------

Read the tutorial-style documentation on [ReadTheDocs](https://slab.readthedocs.io/).

Citing slab
-----------

Schönwiesner et al., (2021). s(ound)lab: An easy to learn Python package for designing and running psychoacoustic experiments. Journal of Open Source Software, 6(62), 3284, https://doi.org/10.21105/joss.03284

```
@article{Schönwiesner2021,
  doi = {10.21105/joss.03284},
  url = {https://doi.org/10.21105/joss.03284},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {62},
  pages = {3284},
  author = {Marc Schönwiesner and Ole Bialas},
  title = {s(ound)lab: An easy to learn Python package for designing and running psychoacoustic experiments.},
  journal = {Journal of Open Source Software}
}
```

Contributing to this project
----------------------------

Anyone and everyone is welcome to contribute. Please take a moment to
review the [guidelines for contributing](CONTRIBUTING.md).

* [Bug reports](CONTRIBUTING.md#bugs)
* [Feature requests](CONTRIBUTING.md#features)
* [Pull requests](CONTRIBUTING.md#pull-requests)

License
-------

The project is licensed under the MIT license.

[![forthebadge made-with-python](http://ForTheBadge.com/images/badges/made-with-python.svg)](https://www.python.org/)
# Contributing to this project

Please take a moment to review this document in order to make the contribution
process easy and effective for everyone involved.

Following these guidelines helps to communicate that you respect the time of
the developers managing and developing this open source project. In return,
they should reciprocate that respect in addressing your issue or assessing
patches and features.


## Using the issue tracker

The issue tracker is the preferred channel for [bug reports](#bugs),
[features requests](#features) and [submitting pull
requests](#pull-requests), but please respect the following restrictions:

* Please **do not** use the issue tracker for personal support requests (use
  [Stack Overflow](http://stackoverflow.com)).

* Please **do not** derail or troll issues. Keep the discussion on topic and
  respect the opinions of others.


<a name="bugs"></a>
## Bug reports

A bug is a _demonstrable problem_ that is caused by the code in the repository.
Good bug reports are extremely helpful - thank you!

Guidelines for bug reports:

1. **Use the GitHub issue search** &mdash; check if the issue has already been
   reported.

2. **Check if the issue has been fixed** &mdash; try to reproduce it using the
   latest `master` or development branch in the repository.

3. **Isolate the problem** &mdash; create a [reduced test
   case](http://css-tricks.com/reduced-test-cases/) and a live example.

A good bug report shouldn't leave others needing to chase you up for more
information. Please try to be as detailed as possible in your report. What is
your environment? What steps will reproduce the issue? What browser(s) and OS
experience the problem? What would you expect to be the outcome? All these
details will help people to fix any potential bugs.

Example:

> Short and descriptive example bug report title
>
> A summary of the issue and the browser/OS environment in which it occurs. If
> suitable, include the steps required to reproduce the bug.
>
> 1. This is the first step
> 2. This is the second step
> 3. Further steps, etc.
>
> `<url>` - a link to the reduced test case
>
> Any other information you want to share that is relevant to the issue being
> reported. This might include the lines of code that you have identified as
> causing the bug, and potential solutions (and your opinions on their
> merits).


<a name="features"></a>
## Feature requests

Feature requests are welcome. But take a moment to find out whether your idea
fits with the scope and aims of the project. It's up to *you* to make a strong
case to convince the project's developers of the merits of this feature. Please
provide as much detail and context as possible.


<a name="pull-requests"></a>
## Pull requests

Good pull requests - patches, improvements, new features - are a fantastic
help. They should remain focused in scope and avoid containing unrelated
commits.

**Please ask first** before embarking on any significant pull request (e.g.
implementing features, refactoring code, porting to a different language),
otherwise you risk spending a lot of time working on something that the
project's developers might not want to merge into the project.

Please adhere to the coding conventions used throughout a project (indentation,
accurate comments, etc.) and any other requirements (such as test coverage).

Follow this process if you'd like your work considered for inclusion in the
project:

1. [Fork](http://help.github.com/fork-a-repo/) the project, clone your fork,
   and configure the remotes:

   ```bash
   # Clone your fork of the repo into the current directory
   git clone https://github.com/<your-username>/<repo-name>
   # Navigate to the newly cloned directory
   cd <repo-name>
   # Assign the original repo to a remote called "upstream"
   git remote add upstream https://github.com/<upstream-owner>/<repo-name>
   ```

2. If you cloned a while ago, get the latest changes from upstream:

   ```bash
   git checkout <dev-branch>
   git pull upstream <dev-branch>
   ```

3. Create a new topic branch (off the main project development branch) to
   contain your feature, change, or fix:

   ```bash
   git checkout -b <topic-branch-name>
   ```

4. Commit your changes in logical chunks. Please adhere to these [git commit
   message guidelines](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html)
   or your code is unlikely be merged into the main project. Use Git's
   [interactive rebase](https://help.github.com/articles/interactive-rebase)
   feature to tidy up your commits before making them public.

5. Locally merge (or rebase) the upstream development branch into your topic branch:

   ```bash
   git pull [--rebase] upstream <dev-branch>
   ```

6. Push your topic branch up to your fork:

   ```bash
   git push origin <topic-branch-name>
   ```

7. [Open a Pull Request](https://help.github.com/articles/using-pull-requests/)
    with a clear title and description.

**IMPORTANT**: By submitting a patch, you agree to allow the project owner to
license your work under the same license as that used by the project.
---
title: 's(ound)lab: An easy to learn Python package for designing and running psychoacoustic experiments.'
tags:
  - Python
  - psychoacoustics
  - audio
  - signal processing
  - teaching
authors:
  - name: Marc Schönwiesner^[co-first author, corresponding author]
    orcid: 0000-0002-2023-1207
    affiliation: "1, 2"
  - name: Ole Bialas^[co-first author]
    affiliation: 1
affiliations:
- index: 1
  name: Institute of Biology, Faculty of Life sciences, Leipzig University, Germany
- index: 2
  name: Institute of Psychology, Faculty of Arts and Sciences, University of Montreal, Canada

date: 30 September 2020
bibliography: paper.bib

---
# Summary
Slab enables researchers and students to prototype and implement psychoacoustic experiments quickly. Slab implements many of the procedures for psychoacoustic research and experiment control and is easily combined with other Python software. A secondary aim of slab is to enable researchers and students without prior training in computer programming and digital signal processing to implement and conduct these experiments. Slab provides building blocks rather than ready-made solutions, so that experimenters still need to carefully consider stimulation, sequencing and data management. This also makes slab very flexible and easy to customise. In the documentation (see slab.readthedocs.io), we provide tutorials suitable for beginners. We also provide experiments conducted in our lab as worked examples.

Slab can:

- generate and manipulate single- and multi-channel sounds
- analyse sound by extracting basic sound features
- aid experimental design through stimulus sequence management and response simulation
- calibrate the experimental setup (loudness calibration and frequency equalisation)
- display and manipulate head-related transfer functions

Below is an example script that estimates the detection threshold for a small change in the location of a sound source (minimum audible angle) with an amplitude-modulated pink noise and a staircase procedure. It illustrates the use of the `Binaural` and `Staircase` classes. The method `present_afc_trial` is a higher-level convenience function to present several sounds and acquire a response from the participant, but each of these steps can be performed separately in one or two lines when implementing non-standard paradigms.
```
# 1up-2down staircase, starting at 20˚ separation:
stairs = slab.Staircase(start_val=20, min_val=0, n_reversals=18)
for angle in stairs:
  # generate fresh noise in each trial:
  midline = slab.Binaural.pinknoise(duration=0.25)
  midline =  midline.am()  # apply amplitude modulation
  midline = midline.ramp()
  # get ITD equivalent to current angle:
  itd = slab.Binaural.azimuth_to_itd(angle, head_radius=10)
  stimulus_left = midline.itd(itd)  # apply the itd
  # apply smoothed head-related transfer function to
  # evoke externalized percept:
  midline = midline.externalize()
  stimulus_left = stimulus_left.externalize()
  # 3 alternatives (target and 2 distractors):
  stairs.present_afc_trial(target=stimulus_left, distractors=[midline]*2)
print(stairs.threshold(n=14))
```

# Statement of need
Slab was written to address our own need for a Python package that allows incoming students and new researchers to implement their own experiments with clean and maintainable code. Experimenters should be able to write and understand the code that they are using. Many students in our lab have very quickly progressed from programming novices to completing psychoacoustic research theses using slab, and we think the package may be useful to others in the same situation. Our approach differs from existing software packages for running behavioural experiments, which provide a high level graphical user interface to customise the parameters of pre-made experiments [@psychopy2; @pychoacoustics]. In our experience, this leads to very little generalisable learning of Python and experimental control. Slab facilitates this learning by providing basic building blocks, implemented concisely in pure Python, that can be used to construct experiments of various levels of complexity.
There is some overlap with librosa [@mcfee2015librosa], a Python package for music analysis, but that package focuses on feature extraction and does not support psychoacoustic experimentation.
Slab is also one of very few packages that features manipulation of head-related transfer functions and a simple API for reading a standard file format (SOFA) for such data. There is overlap with more recent implementations of the complete SOFA API [@pysofaconventions; @python-sofa], but these packages provide no methods for typical experimental manipulations of head-related transfer functions. We will likely use `pysofaconventions` internally for handling SOFA files within `slab` in the near future.
The architecture of the `Signal` class and some of the sound generation methods in the `Sound` class are inspired on the 1.4 version of Brian.hears [@brian2hears], but we made several simplifications based on learning reports from students. For instance, signal objects do not implement buffering and do not inherit from Numpy arrays [@harris2020array] directly, because these features significantly hindered students' understanding of the code.

# Audience
Slab is directed towards students and researchers of all levels studying the perception of sound.
Researchers and incoming students at our lab use it routinely in behavioural and neuroimaging experiments, and the package and has been used in several graduate courses psychophysics and auditory neuroscience.

# References
---
name: Bug report
about: Report a bug and help us improve slab
title: "[BUG in class ...]"
labels: ''
assignees: DrMarc

---

**Describe the bug**
A clear and concise description of what the bug is.
For instance: "Binaural.pinknoise() should return a stereo sound (two channels), but instead returns a mono sound (one channel)."

**To Reproduce**
Minimal code example. This should ideally run as is, without additional data. Include all imports.
For instance:
> import slab
> sig = slab.Binaural.pinknoise()
> sig.n_channels

**Expected behaviour**
The expected output of the minimal code example. For instance, in the example above:
> 2

**Actual behaviour**
The output of the minimal code example. For instance, in the example above:
> 1

**OS and Version:**
 - OS: [e.g. MacOS, Win, Linux]
 - Version [e.g. v0.9]
