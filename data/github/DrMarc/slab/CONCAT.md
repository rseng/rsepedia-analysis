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
.. _Hrtfs:

HRTFs
=====

A head-related transfer function (HRTF) describes the impact of the listeners ears, head and torso on incoming sound
for every position in space. Knowing the listeners HRTF, you can simulate a sound source at any position by filtering
it with the transfer function corresponding to that position. The :class:`HRTF` class provides methods for
manipulating, plotting, and applying head-related transfer functions.

Reading HRTF data
-----------------
Typically the :class:`HRTF` class is instantiated by loading a file. The canonical format for HRTF-data is called
sofa (Spatially Oriented Format for Acoustics). To read sofa files, you need to install the h5netcdf module:
`pip install h5netcdf`. The module includes a set of standard HRTF recordings from the KEMAR (a mannequin for acoustic
recordings). You can get the path to the folder containing the recordings with the :func:`data_path` function. The
first time you call this function, the recordings will be downloaded from the sofa website. You can read them by
calling the :class:`HRTF` class with the name of the file as an argument. Print the resulting object to obtain
information about the structure of the HRTF data ::

    hrtf = slab.HRTF.kemar()
    print(hrtf)
    # <class 'hrtf.HRTF'> sources 710, elevations 14, samples 710, samplerate 44100.0

Libraries of many other recordings can be found on the `website of the sofa file format <https://www.sofaconventions.org/>`_.

.. note:: The class is at the moment geared towards plotting and analysis of HRTF files in the `sofa format <https://www.sofaconventions.org/>`_, because we needed that functionality for grant applications. The functionality will grow as we start to record and manipulate HRTFs more often.

.. note:: When we started to writing this code, there was no python module for reading and writing sofa files. Now that `pysofaconventions <https://github.com/andresperezlopez/pysofaconventions>`_ is available, we will at some point switch internally to using that module as backend for reading sofa files, instead of our own limited implementation.

Plotting sources
--------------------
The HRTF is a set of many transfer functions, each belonging to a certain sound source position (for example,
there are 710 sources in the KEMAR recordings). You can plot the source positions in 3D with the :meth:`.plot_sources`
to get an impression of the density of the recordings. The red dot indicates the position of the listener and the red
arrow indicates the lister's gaze direction. Optionally, you can supply a list of source indices which will be
highlighted in red. This can be useful when you are selecting source locations for an experiment and want to confirm
that you chose correctly. In the example below, we select sources using the methods :meth:`.elevation_sources`, which
selects sources along a horizontal sphere slice at a given elevation and :meth:`.cone_sources`, which selects sources
along a vertical slice through the source sphere in front of the listener a given angular distance away from the
midline:

.. plot::
    :include-source:
    :context:

    # cone_sources and elevation_sources return lists of indices which are concatenated by adding:
    hrtf = slab.HRTF.kemar()
    sourceidx = hrtf.cone_sources(0) + hrtf.elevation_sources(0)
    hrtf.plot_sources(sourceidx) # plot the sources in 3D, highlighting the selected sources

Try a few angles for the :meth:`.elevation_sources` and :meth:`.cone_sources` methods to understand how selecting
the sources works!

Plotting transfer functions
---------------------------
As mentioned before, a HRTF is collection of transfer functions. Each single transfer function is an instance of the
:class:`slab.Filter` with two channels - one for each ear. The transfer functions are located in the :attr:`data`
list and the coordinates of the corresponding sources in the :attr:`sources` list. In the example below, we select a
source, print it's coordinates and plot the corresponding transfer function.

.. plot::
    :include-source:
    :context: close-figs

    from matplotlib import pyplot as plt
    hrtf = slab.HRTF.kemar()
    fig, ax = plt.subplots(1)
    idx = 10
    source = hrtf.sources[idx]  # the source's azimuth, elevation and distance
    filt = hrtf.data[idx] # the corresponding filter
    fig.suptitle(f"source at azimuth {source[0].round(2)} and elevation {source[1]}")
    filt.channel(0).tf(axis=ax, show=False)
    filt.channel(1).tf(axis=ax, show=Fals)
    plt.legend()
    plt.show()

The :class:`HRTF` class also has a :meth:`.plot_tf` method to plot transfer functions as either `waterfall`
(as is Wightman and Kistler, 1989), `image` plot (as in Hofman 1998). The function takes a list of source indices as an
argument which will be included in the plot. The function below shows how to generate a `waterfall` and `image` plot
for the sources along the central cone. Before plotting, we apply a diffuse field equalization to remove non-spatial
components of the HRTF, which makes the features of the HRTF that change with direction easier to see:

.. plot::
    :include-source:
    :context: close-figs

    from slab import data_path
    from matplotlib import pyplot as plt
    hrtf = slab.HRTF.kemar()
    fig, ax = plt.subplots(2)
    dtf = hrtf.diffuse_field_equalization()
    sourceidx = hrtf.cone_sources(0)
    ax[0].set_title("waterfall plot")
    ax[1].set_title("image plot")
    hrtf.plot_tf(sourceidx, ear='left', axis=ax[0], show=False, kind="waterfall")
    hrtf.plot_tf(sourceidx, ear='left', axis=ax[1], show=False, kind="image")
    plt.tight_layout()
    plt.show()


As you can see the HRTF changes systematically with the elevation of the sound source, especially for frequencies above
6 kHz. Individual HRTFs vary in the amount of spectral change across elevations, mostly due to differences in the
shape of the ears. You can compute a measure of the HRTFs spectral dissimilarity the vertical axis, called vertical
spatial information (VSI, `Trapeau and Schönwiesner, 2016 <https://pubmed.ncbi.nlm.nih.gov/27586720/>`_).
The VSI relates to behavioral localization accuracy in the vertical dimension: listeners with acoustically more
informative spectral cues tend to localize sounds more accurately in the vertical axis. Identical filters give a VSI
of zero, highly dissimilar filters give a VSI closer to one. The hrtf has to be diffuse-field equalized for this
measure to be sensible, and the :meth:`.vsi` method will apply the equalization. The KEMAR mannequin have a VSI
of about 0.73::

    hrtf.vsi()
    # .73328

The :meth:`.vsi` method accepts arbitrary lists of source indices for the dissimilarity computation.
We can for instance check how the VSI changes when sources further off the midline are used. There are some reports
in the literature that listeners can perceive the elevation of a sound source better if it is a few degrees to the
side. We can check whether this is due to more dissimilar filters at different angles (we'll reuse the `dtf` from above
to avoid recalculation of the diffuse-field equalization in each iteration)::

    for cone in range(0,51,10):
        sources = dtf.cone_sources(cone)
        vsi = dtf.vsi(sources=sources, equalize=False)
        print(f'{cone}˚: {vsi:.2f}')
        # 0˚: 0.73
        # 10˚: 0.63
        # 20˚: 0.69
        # 30˚: 0.74
        # 40˚: 0.76
        # 50˚: 0.73

The effect seems to be weak for KEMAR, (VSI falls off for directions slightly off the midline and then increases again at around 30-40˚).


Virtually displaying 3D sound
-----------------------------
The HRTF describes the directional filtering of incoming sounds by the listeners ears, head and torso. Since this is the
basis for localizing sounds in three dimensions, we can apply the HRTF to a sound to evoke the impression of it coming
from a certain direction in space when played through headphones. The :meth:`HRTF.apply` method returns an instance of
the :class:`slab.Binaural`. It is important to use the :meth:`~HRTF.apply` method of the :class:`HRTF` class instead of
the :meth:`~Filter.apply` method of the individual :class:`Filter` class objects in the HRTF, because only the
meth:`HRTF.apply` method conserves ITDs. The :meth:`Filter.apply` method does not do that, because when applying a
generic filter, you normally do not want to introduce delays.

In the example below we apply the transfer functions corresponding to three sound sources at
different elevations along the vertical midline to white noise.

.. plot::
    :include-source:
    :context: close-figs

    from slab import data_path, Sound
    from matplotlib import pyplot as plt
    hrtf = slab.HRTF.kemar()
    sound = slab.Sound.pinknoise(samplerate=hrtf.samplerate)  # the sound to be spatialized
    fig, ax = plt.subplots(3)
    sourceidx = [0, 260, 536]  # sources at elevations -40, 0 and 40
    spatial_sounds = []
    for i, index in enumerate(sourceidx):
        spatial_sounds.append(hrtf.apply(index, sound))
        # only plot frequencies above 5kHz because low frequencies are unaffected by the HRTF
        spatial_sounds[i].spectrum(axis=ax[i], low_cutoff=5000, show=False)
    plt.show()

You can use the :meth:`~Sound.play` method of the sounds to listen to them - see if you can identify the virtual sound
source position. Your ability to do so depends on how similar your own HRTF is to that of the the KEMAR artificial head.
Your auditory system can get used to new HRTFs, so if you listen to the KEMAR recordings long enough they will eventually
produce virtual sound sources at the correct locations.

Binaural filters from the KEMAR HRTF will impose the correct spectral profile, but no ITD. After applying an HRTF filter
corresponding to an off-center direction, you should also apply an ITD corresponding to the direction using the
:meth:`Binaural.azimuth_to_itd` and :meth:`Binaural.itd` methods.

Finally, the HRTF filters are recorded only at certain locations (710, in case of KEMAR - plot the source locations to
inspect them). You can interpolate a filter for any location covered by these sources with the :meth:`HRTF.interpolate`
method. It triangulates the source locations and finds three sources that form a triangle around the requested location
and interpolate a filter with a (barycentric) weighted average in the spectral domain. The resulting filter may not have
the same overall gain, so remember to set the level of your stimulus after having applied the interpolated HRTF.

.. currentmodule:: slab.experiments

Worked examples
===============

The folder slab.experiments contains the full code from actual psychoacoustic experiments in our lab. We use this folder mainly to make the code available and enable easy replication. The examples are well documented and may give you and idea of the typical structure of such experiments. To run these experiments, import them from slab.experiments::

    from slab.experiments import motion_speed
    motion_speed.main_experiment(subject='test')

Currently available are:

.. autofunction:: slab.experiments.room_voice_interference.main_experiment

.. autofunction:: slab.experiments.motion_speed.main_experiment


Quick standard experiments:
---------------------------

.. _audiogram:

Audiogram
^^^^^^^^^
Run a pure tone audiogram at the standard frequencies 125, 250, 500, 1000, 2000, 4000 Hz using an adaptive staircase: ::

    from matplotlib import pyplot as plt
    freqs = [125, 250, 500, 1000, 2000, 4000]
    threshs = []
    for frequency in freqs:
        stimulus = slab.Sound.tone(frequency=frequency, duration=0.5)
        stairs = slab.Staircase(start_val=50, n_reversals=18)
        print(f'Starting staircase with {frequency} Hz:')
        for level in stairs:
            stimulus.level = level
            stairs.present_tone_trial(stimulus)
            stairs.print_trial_info()
        threshs.append({stairs.threshold())
        print(f'Threshold at {frequency} Hz: {stairs.threshold()} dB')
    plt.plot(freqs, threshs) # plot the audiogram


Temporal modulation transfer function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Measure temporal modulation transfer functions via detection thresholds for amplitude modulations. The parameters of the test replicate Fig. 2 in Viemeister [1979]_ and present sinusoidal 2 to 4000 Hz modulations in a 77-dB wideband noise carrier using an adaptive staircase. ::

    from matplotlib import pyplot as plt
    mod_freqs = [2, 4, 8, 16, 32, 64, 125, 250, 500, 1000, 2000, 4000]
    threshs = []
    base_stimulus = slab.Sound.pinknoise(duration=1.)
    base_stimulus.level = 77
    for frequency in mod_freqs:
    stairs = slab.Staircase(start_val=0.8, n_reversals=16, step_type='db',
                step_sizes=[4,2], min_val=0, max_val=1, nup=1, ndown=2)
        print(f'Starting staircase with {frequency} Hz:')
        for depth in stairs:
            stimulus = base_stimulus.am(frequency=frequency, depth=depth)
            stairs.present_afc_trial(stimulus, base_stimulus)
        threshs.append(stairs.threshold(n=14))
        print(f'Threshold at {frequency} Hz: {stairs.threshold(n=14)} modulation depth')
    plt.plot(freqs, threshs) # plot the transfer function


.. [1979] Viemeister (1979) Temporal modulation transfer functions based upon modulation thresholds. JASA 66(5), 1364–1380
.. currentmodule:: slab

.. _Reference:

Reference documentation
=======================

.. note:: This reference documentation is auto-generated from the doc strings in the module. For a tutorial-like overview of the functionality of slab, please see the previous sections.

Sounds
^^^^^^
Inherits from :class:`slab.Signal`.

.. autoclass:: Sound
   :members:
   :member-order: bysource

.. automethod:: slab.sound.apply_to_path

Signal
------
:class:`slab.Sound` inherits from Signal, which provides basic methods to handle signals:

.. autoclass:: Signal
   :members:
   :member-order: bysource

Binaural sounds
---------------
Binaural sounds inherit from Sound and provide methods for manipulating interaural parameters of two-channel sounds.

.. autoclass:: Binaural
  :members:
  :member-order: bysource

Psychoacoustic procedures
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoclass:: Trialsequence
   :members:
   :inherited-members:
   :member-order: bysource

.. autoclass:: Staircase
   :members:
   :inherited-members:
   :member-order: bysource

.. autoclass:: Precomputed
   :members:
   :member-order: bysource

.. autoclass:: ResultsFile
   :members:
   :member-order: bysource

.. automethod:: slab.psychoacoustics.key

.. automethod:: slab.psychoacoustics.load_config


Filters
-------
.. autoclass:: Filter
   :members:
   :member-order: bysource

HRTFs
-----
.. autoclass:: HRTF
   :members:
   :member-order: bysource
.. currentmodule:: slab

.. _Filters:

Filters
=======
The :class:`Filter` class can be used to generate, manipulate and save filter banks and transfer functions. Filters are represented internally as :class:`Signal` and come in two flavours: finite impulse responses (FIR) and frequency bin amplitudes (FFT). The :attr:`fir` (True or False).

Simple Filters
--------------
Simple low-, high-, bandpass, and bandstop filters can be used to suppress selected frequency bands in a sound. For example, if you don't want the sound to contain power above 1 kHz, apply a 1 kHz lowpass filter:

.. plot::
    :include-source:
    :context:

    from matplotlib import pyplot as plt
    sound = slab.Sound.whitenoise()
    filt = slab.Filter.band(frequency=1000, kind='lp')
    sound_filt = filt.apply(sound)
    _, [ax1, ax2] = plt.subplots(2, sharex=True)
    sound.spectrum(axis=ax1, color="blue")
    sound_filt.spectrum(axis=ax2, color="red")

The :meth:`~Sound.filter` of the :class:`Sound` class wraps around :meth:`Filter.cutoff_filter` and :meth:`Filter.apply` so that you can use these filters conveniently from within the :class:`Sound` class.

Filter design is tricky and it is good practice to plot and inspect the transfer function of the filter:

.. plot::
    :include-source:
    :context: close-figs

    filt.tf()


Inspecting the waveform of the sound (using the :meth:`~slab.Sound.waveform` method) or the rms level (using the :attr:`slab.Sound/level` attribute) shows that the amplitude of the filtered signal is smaller than that of the original, because the filter has removed power. You might be tempted to correct this difference by increasing the level of the filtered sound, but this is not recommended because the perception of intensity (loudness) depends non-linearly on the frequency content of the sound.

Filter banks
------------
A :class:`Filter` objects can hold multiple channels, just like a :class:`Sound` object. In the following, we will refer
to filters with multiple channels as filter banks. You can create filter banks the same way you create multi-channel
sound (since the :class:`Filter` and :class:`Sound` class both inherit from the parent :class:`Signal` class).
When you apply a filter bank to a single sound, each filter will be applied to a separate copy of the sound and the
:func:`apply` function will return a sound with a number of channels equal to the number of filters in the bank.
This way you can create, for example, a series of sounds with different frequency bands:

.. plot::
    :include-source:
    :context: close-figs

    filters = []
    low_cutoff_freqs = [500, 1000, 1500]
    high_cutoff_freqs = [1000, 1500, 2000]
    for low, high in zip(low_cutoff_freqs, high_cutoff_freqs):
        filters.append(slab.Filter.band(frequency=(low, high), kind='bp'))
    fbank = slab.Filter(filters)  # put the list into a single filter object
    sound_filt = fbank.apply(sound)  # apply each filter to a copy of sound
    # plot the spectra, each color represents one channel of the filtered sound
    _, ax = plt.subplots(1)
    sound_filt.spectrum(axis=ax, show=False)
    ax.set_xlim(100, 5000)
    plt.show()

.. plot::
    :context: close-figs

The channels, or subbands, of the filtered sound can be modified and re-combine with the :meth:`.combine_subbands`
method. An example of this process is the vocoder implementation in the :class:`Sound` class, which uses these features
of the :class:`Filter` class. The multi-channel filter is generated with :meth:`.cos_filterbank`, which produces
cosine-shaped filters that divide the sound into small frequency bands which are spaced in a way that mimics the
filters of the human auditory periphery (`equivalent rectangular bandwidth, ERB <https://en.wikipedia.org/wiki/
Equivalent_rectangular_bandwidth>`_). Here is an example of the transfer functions of this filter bank:

.. plot::
    :include-source:
    :context: close-figs

    fbank = slab.Filter.cos_filterbank()
    fbank.tf()

A speech signal is filtered with this bank, and the envelopes of the subbands are computed using the
:meth:`envelope` method of the :class:`Signal` class. The envelopes are filled with noise, and the
subbands are collapsed back into one sound. This removes most spectral information but retains temporal information
in a speech signal and sound a bit like whispered speech. Here are the essential bits of code from the
:meth:`~slab.Sound.vocode` method to illustrate the use of a filter bank. The first line records a speech sample
from the microphone (say something!)::

    signal = slab.Sound.record() # record a 1 s speech sample from the microphone
    fbank = slab.Filter.cos_filterbank(length=signal.n_samples) # make the filter bank
    subbands = fbank.apply(signal) # get a sound channel for each filter channel
    envs = subbands.envelope() # now get the envelope of each frequency band...
    noise = slab.Sound.whitenoise()
    subbands_noise = fbank.apply(noise)
    subbands_noise *= envs  # ... and fill them with noise
    subbands_noise.level = subbands.level # keep subband level of original
    vocoded = slab.Filter.collapse_subbands(subbands_noise, filter_bank=fbank)
    vocoded.play()


If you want to apply multiple filters to the same sound in sequence without creating subbands, you can simply use a for
loop. For example, you could remove different parts of the spectrum using bandstop filters::

    sound = slab.Sound.whitenoise()
    # create a filter bank which consists of three separate bandstop filters
    filter_bank = slab.Filter([slab.Filter.band(kind="bs", frequency=f) for f in [(200, 300), (500, 600), (800, 900)]])
    for i in range(filter_bank.n_channels):
        sound = filter_bank.channel(i).apply(sound)


If the a one-channel filter is applied to a multi-channel sound, the filter will be applied to each
channel individually. This can be used, for example, to easily pre-process a set of recordings (where
every recordings is represented by a channel in the :class:`slab.Sound` object). If a multi-channel filter
is applied to a multi-channel signal with the same number of channels each filter channel is applied to
the corresponding signal channel. This mechanism is used, for example, during the equalization of a set of loudspeakers.

Equalization
------------
In Psychoacoustic experiments, we are often interested in the effect of a specific feature. One could,
for example, take the bandpass filtered sounds from the example above and investigate how well listeners
can discriminate them from a noisy background - a typical cocktail-party task. However, if the transfer
function of the loudspeakers or headphones used in the experiment is not flat, the findings will be biased.
Imagine that the headphones used were bad at transmitting frequencies below 1000 Hz. This would make a sound
with center frequency of 550 Hz harder to detect than one with a center frequency of 1550 Hz. To prevent this from
happening, we have to equalize the headphones' transfer function. You can measure the
transfer function of your system by playing a wide-band sound, like a chirp, and recording it with a probe microphone
(which itself must have a flat transfer function). From this recording, you can calculate the transfer function, which
is basically the difference in the power spectrum of the played sound and the recording. We can take the opposite of
that difference to create an inverse filter. Apply the inverse filter to a sound before playing it through that system
to compensate for the uneven transfer, because the inverse filter and the actual transfer function cancel each other.
The :meth:`~slab.Filter.equalizing_filterbank` method does most of this work for you. For a demonstration,
we simulate a (pretty bad) loudspeaker transfer function by applying a random filter:

.. plot::
    :include-source:
    :context: close-figs

    import random
    freqs = [f * 400 for f in range(10)]
    gain = [random.random()+.4 for _ in range(10)]
    tf = slab.Filter.band(frequency=freqs, gain=gain)
    sound = slab.Sound.whitenoise()
    recording = tf.apply(sound)
    recording.spectrum()

With the original sound and the simulated recording we can compute an inverse filter und pre-filter the sound
(or in this case, just filter the recording) to achieve a nearly flat playback through our simulated bad loudspeaker:

.. plot::
    :include-source:
    :context: close-figs

    inverse = slab.Filter.equalizing_filterbank(reference=sound, sound=recording)
    equalized = inverse.apply(recording)
    equalized.spectrum()

If there are multiple channels in your recording (assembled from recordings of the same white noise through several
loudspeakers, for instance) then the :meth:`~slab.Filter.equalizing_filterbank` method returns a filter bank with one
inverse filter for each signal channel, which you can :meth:`~slab.Filter.apply` just as in the example above.

**slab**: easy manipulation of sounds and psychoacoustic experiments in Python
==============================================================================


**Slab** ('es-lab', or sound laboratory) is an open source project and Python package that makes working with sounds and running psychoacoustic experiments simple, efficient, and fun! For instance, it takes just eight lines of code to run a pure tone audiogram using an adaptive staircase: :ref:`audiogram`



Why slab?
---------
The package aims to lower the entrance barrier for working with sounds in Python and provide easy access to typical operations in psychoacoustics, specifically for students and researchers in the life sciences. The typical BSc or MSc student entering our lab has limited programming and signal processing training and is unable to implement a psychoacoustic experiment from scratch within the time limit of a BSc or MSc thesis. Slab solves this issue by providing easy-to-use building blocks for such experiments. The implementation is well documented and sufficiently simple for curious students to understand. All functions provide sensible defaults and will in many cases 'just work' without arguments (`vowel = slab.Sound.vowel()` gives you a 1-second synthetic vowel 'a', `vowel.spectrogram()` plots the spectrogram). This turned out to be useful for teaching and demonstrations. Many students in our lab have now used the package to implement their final projects and exit the lab as proficient Python programmers.

.. _installation:

Installation
------------

Install the current stable release from the python package index with pip::

    pip install slab

or get the latest development version directly from GitHub (if you have `git <https://git-scm.com>`_) by running::

    pip git+https://github.com/DrMarc/slab.git

**The current version of slab is** |version|.

The releases use `semantic versioning <https://semver.org>`_: ``major.minor.patch``, where ``major`` increments for changes that break backwards compatibility, ``minor`` increments of added functionality, and ``patch`` increases for internal bug fixes.
```slab.__version__``` prints the installed version.

To run the tests::

    pip install slab[testing]

Then go to the installation directory and run::

    pytest

On Linux, you may need to install libsndfile (required by SoundFile) using your distribution's package manager, for instance::

    sudo apt-get install libsndfile1

On Windows, you may need to install `windows-curses <https://pypi.org/project/windows-curses/>`_ (required for getting button presses in the psychoacoustics classes)::

    pip install windows-curses

Working with head related transfer functions requires the h5netcdf module (trying to load a hrtf file will raise an error and tell you to install::

    pip install h5netcdf

All other dependencies should have been automatically installed, and you should see meaningful errors if that did not happen for some reason. The dependencies are: numpy, scipy.signal (for filtering and several other DSP functions), matplotlib (for all plotting), SoundFile (for reading and writing wav files), curses or windows-curses (for getting key presses), and SoundCard (for playing and recording sounds). We have seen a hard-to-replicate problem on some Macs with the SoundCard module: a pause of several seconds after a sound is played. *Macs with M1 processors* also have an issue with SoundCard. If you experience these issues, just uninstall SoundCard::

    pip uninstall SoundCard

Slab will then use another method to play sounds (winsound on Windows, afplay on Macs, and `SoX <http://sox.sourceforge.net>`_ on Linux), and will record sounds from the microphone using SoX. There are many other packages to play sounds, depending on our operating system. If you prefer a different one, you can easily modify or replace the :meth:`~slab.Sound.play` method.


Citing slab
-----------

Schönwiesner et al., (2021). s(ound)lab: An easy to learn Python package for designing and running psychoacoustic experiments. Journal of Open Source Software, 6(62), 3284, https://doi.org/10.21105/joss.03284

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

To return the citations as strings, you can also call::

    slab.cite()


.. toctree::
  :caption: Contents
  :maxdepth: 2

  introduction
  sounds
  psychoacoustics
  filter
  hrtf
  Worked examples <examples>
  Reference documentation <reference>

**Index of functions and classes:** :ref:`genindex`

**Search the documentation:** :ref:`search`
.. _Psychoacoustics:

Psychoacoustics
===============
The :class:`Psychoacoustics` class simplifies psychoacoustic experiments by providing classes and methods
for trial sequences and adaptive staircases, results and configuration files, response collection via keyboard and
button boxes, and handling of collections of precomputed stimuli. This all-in-one approach makes for clean code and
easy data management.

Trial sequences
---------------
Experiments are often defined by a sequence of trials of different conditions. This sequence is generated before the
experiment according to certain rules. In the most basic case, a set of experimental conditions are repeated a number of
times pseudorandom order. Such experiments can be handled by the :class:`Trialsequence` class. To generate an
instance of :class:`Trialsequence` you define a list of ``conditions`` and specify how often each of them is
repeated (``n_reps``). You can also specify the ``kind`` of list you want to generate: "non_repeating" means that
the same condition will not appear twice in a row, "random_permutation" means that the order is completely randomised.
For example, generate pure tones with different frequencies and play them in non-repeating, randomised order.::

  freqs = [495, 498, 501, 504]  # frequencies of the tones
  seq = slab.Trialsequence(conditions=freqs, n_reps=10)  # 10 repetitions per condition
  # now we draw elements from the list, generate a tone and play it until we reach the end:
  for freq in seq:
    stimulus = slab.Sound.tone(frequency=freq)
    stimulus.play()

Usually, we do not only want to play sounds to the participants in our experiment. Instead, we want them to perform some
kind of task and give a response. In the example above we could, for instance, ask after every tone if that tone was
higher or lower in frequency than the previous one. The response is captured with the :meth:`~slab.psychoacoustics.key`
context manager which can record single button presses (using either the :mod:`curses` module or the :meth:`key_press_event`
of the `stairs` plot, see :ref:`_responses`). In our example, we instruct the subject to press "y" (yes) if the played
tone was higher then the previous and "n" (no) if it was lower (a 1-back task). After each trial we check if the
response was correct and store that information as 1 (correct) or 0 (wrong) in the trial sequence.::

  for freq in seq:
    stimulus = slab.Sound.tone(frequency=freq)
    stimulus.play()
    if seq.this_n > 0:  # don't get response for first trial
      previous = seq.get_future_trial(-1)
      with slab.key() as key:  # wait for a key press
        response = key.getch()
      # check if the response was correct, if so store a 1, else store 0
      if (freq > previous and response == ord('y')) or (freq<previous and response == ord('n')):
        seq.add_response(1)
      else:
        seq.add_response(0)
  seq.save_json("sequence.json")  # save the trial sequence and response

There are two ways for a response to be correct in this experiment. Either the frequency of the stimulus was higher
than the last one and the 'y' key was pressed, or it was lower and the 'n' key was pressed. (:func:`ord()` is used to
get the key codes of the 'y' and 'n' keys (112 and 110, respectively). All other options, including missed responses,
are counted as wrong answers. Since we encoded correct responses as 1 and wrong responses as 0, we could just sum over
the list of responses and divide by the length of the list to get the fraction of trials that was answered correctly.

Kinds of trial sequences
^^^^^^^^^^^^^^^^^^^^^^^^
Trial sequences are useful for non-adaptive testing (the current stimulus does not depend on the listeners previous
responses) and other situations where you need a controlled sequence of stimulus values. The :class:`Trialsequence`
class constructs several controlled sequences (random permutation, non-repeating, infinite, oddball), computes
transition probabilities and condition frequencies, and can keep track of responses::

    # sequence of 5 conditions, repeated twice, without direct repetitions:
    seq = slab.Trialsequence(conditions=5, n_reps=2)

    # infinite sequence of color names:
    seq = slab.Trialsequence(conditions=['red', 'green', 'blue'], kind='infinite')

    # stimulus sequence for an oddball design:
    seq = slab.Trialsequence(conditions=1, deviant_freq=0.12, n_reps=60)

The list of trials is contained in the :attr:`trials` of the :class:`Trialsequence` object, but you don't normally need
to access this list directly. A :class:`Trialsequence` object can be used like a :class:`Staircase` object in a
listening experiment and will return the current stimulus value when used in a loop. Below is
:ref:`the detection threshold task <detection_example>` from the :class:`Staircase`, rewritten using Fechner's method of
constant stimuli with a :class:`Trialsequence`::

    stimulus = slab.Sound.tone(duration=0.5)
    levels = list(range(0, 50, 10)) # the sound levels to test
    trials = slab.Trialsequence(conditions=levels, n_reps=10) # each repeated 10 times
    for level in trials:
        stimulus.level = level
        stimulus.play()
        with slab.key() as key:
            response = key.getch()
        trials.add_response(response)
    trials.response_summary()

Because there is no simple threshold, the :class:`Trialsequence` class provides a :meth:`.response_summary`, which
tabulates responses by condition index in a nested list.

The infinite kind of :class:`Trialsequence` is perhaps less suitable for controlling the stimulus parameter of interest,
but it is very useful for varying other stimulus attributes in a controlled fashion from trial to trial (think of
'roving' paradigms). Unlike when selecting a random value in each trial, the infinite :class:`Trialsequence` guarantees
locally equal value frequencies, avoids direct repetition, and keeps a record in case you want to include the sequence as
nuisance covariate in the analysis later on. Here is a real-world example from an experiment with pseudo-words, in which
several words without direct repetition were needed in each trial. word_list contained the words as strings, later used
to load the correct stimulus file::

    word_seq = slab.Trialsequence(conditions=word_list, kind='infinite')
    word = next(word_seq) # draw a word from the list

This is one of the very few cases where it makes sense to get the next trial by calling Python's :func:`next` function,
because this is not the main trial sequence. The main trial sequence (the one determining the values of your main
experimental parameter) should normally be used in a `for` loop as in the previous example.

Controlling transitions
^^^^^^^^^^^^^^^^^^^^^^^^
While randomized sequences do the job most of the time, in some cases it is necessary to control the transitions
between the individual conditions more tightly. For instance, you may want to ensure nearly equal transitions,
or avoid certain combinations of subsequent conditions entirely. The :meth:`.transitions` method counts, for each
condition, how often every other condition follows this one. You can divide the count by the number of repetitions in
the sequence to get the transitional probabilities::

    trials = slab.Trialsequence(conditions=4, n_reps=10)
    trials.transitions()
    out:
    array([[0., 2., 6., 2.],
           [3., 0., 0., 7.],
           [2., 6., 0., 1.],
           [4., 2., 4., 0.]])
    trials.transitions() / 10  # divide by n_reps to get the probability
    out:
    array([[0. , 0.2, 0.6, 0.2],
           [0.3, 0. , 0. , 0.7],
           [0.2, 0.6, 0. , 0.1],
           [0.4, 0.2, 0.4, 0. ]])

The diagonal of this array contains only zeroes, because a condition cannot follow itself in the default
``non_repeating`` trial sequence. The other entries are uneven; for instance, condition 1 is followed by condition
3 seven times, but never by condition 2. If you want near-equal transitions, then you could generate sequences in a
loop until a set condition is fulfilled, for instance, no transition > 4::

    import numpy
    trans = 5
    while numpy.any(trans>4):
        trials = slab.Trialsequence(conditions=4, n_reps=10)
        trans = trials.transitions()
    print(trans)
    out:
    array([[0., 3., 3., 3.],
           [4., 0., 3., 3.],
           [3., 4., 0., 3.],
           [3., 3., 4., 0.]])

If your condition is more complicated, you can perform several tests in the loop body and set a flag that determines
when all have been satisfied and the loop should be end.
But be careful, setting these constraints too tightly may result in an infinite loop.

Alternative Choices
^^^^^^^^^^^^^^^^^^^
Often, an experimental paradigm requires more complex responses than yes or no. A common option is the classical
"forced choice" paradigm, in which the subject has to pick a response from a defined set of responses. Since this is a
common paradigm, the :class:`Trialsequence` and :class:`Staircase` class have a method for it called
:meth:`present_afc_trial` (afc stands for alternative forced choice). With this function we can make our frequency
discrimination task from the example above a bit more elaborate. We define the frequencies of our target tones and add
two distractor tones with a frequency of 500 Hz. In each trial, all three tones (target + 2 x distractor) are played in
random order. The participant answers the question: "which tone was different from the others?" and responds by pressing
the key "1", "2" or "3". All of this can be done in only 6 lines of code: ::

    distractor = slab.Sound.tone(duration=0.5, frequency=500)
    freqs = list(range(495, 505))
    trials = slab.Trialsequence(conditions=freqs, n_reps=2)
    for freq in trials:
        target = slab.Sound.tone(frequency=freq, duration=0.5)
        trials.present_afc_trial(target, [distractor, distractor], isi=0.2)

Adaptive staircases
-------------------
In many cases, you do not want to test every condition with the same frequency, but adapt the stimulus presentation to
the responses of the participant. For example, when measuring an audiogram, you want to spend most of the testing time
around the threshold to make the testing efficient. The :class:`Staircase` class lets you do that. You pick an initial
value for the stimulus parameter (``start_val``) and a step size (``step_sizes``). With each trial, the starting value
is decreased by one step size until the subject is not able to respond correctly anymore. Then it is increased step wise
until the response is correct again, then decreased again and so on. This procedure is repeated until the given number
of reversals (``n_reversals``) is reached. The step size can be a list in which case the current step size moves one
index in the list by each reversal until the end of the list is reached.
For example, we could use a step size of 4 until we crossed the threshold for the first time, then use a step size of
1 for the rest of the experiment. This ensures that we get to the threshold quickly and, once we are there, measure
it precisely. (The :meth:`simulate_response` method used here is explained under :ref:`_simulating`.)

.. plot::
    :include-source:

    stairs = slab.Staircase(start_val=10, n_reversals=18, step_sizes=[4,1])
    for stimulus_value in stairs:
        response = stairs.simulate_response(threshold=3) # simulate subject's response
        stairs.add_response(response) # initiates calculation of next stimulus value
        stairs.plot()

Calling the plot function in the for loop (*after* :meth:`Staircase.add_response`) will update the plot each
trial and let you monitor the performance of the participant, including the current stimulus value (grey dot), and
correct/incorrect responses (green and red dots). (On some Windows systems, the plot captures the focus and may prevent
you from entering responses in the terminal window. In that case, switch the :data:`slab.psychoacoustics.input_method`
to 'figure'. This will get a button press through the stairs figure's :meth:`key_press_event`.)

An audiogram is a typical example for a staircase procedure. We can define a list of frequencies and run a
staircase for each one. Afterwards we can print out the result using the :meth:`thresh()` method.::

    from matplotlib import pyplot as plt
    freqs = [125, 250, 500, 1000, 2000, 4000]
    threshs = []
    for frequency in freqs:
        stimulus = slab.Sound.tone(frequency=frequency, duration=0.5)
        stairs = slab.Staircase(start_val=50, n_reversals=18)
        print(f'Starting staircase with {frequency} Hz:')
        for level in stairs:
            stimulus.level = level
            stairs.present_tone_trial(stimulus)
        threshs.append(stairs.threshold())
        print(f'Threshold at {frequency} Hz: {stairs.threshold()} dB')
    plt.plot(freqs, threshs) # would plot the audiogram

:meth:`present_tone_trial()` is a convenience method that presents the trial, acquires a response, and optionally prints
trial information. All of this can be done explicitly, as shown in the :class:`Trialsequence` example.

Staircase Parameters
^^^^^^^^^^^^^^^^^^^^
Setting up a near optimal staircase requires some expertise and pilot data. Practical recommendations can be found in
`García-Pérez (1998) <https://pubmed.ncbi.nlm.nih.gov/9797963/>`_. ``start_val`` sets the stimulus value presented in
the first trial and the starting point of the staircase. This stimulus should in general be easy to detect/discriminate
for all participants. You can limit the range of stimulus values between ``min_val`` and ``max_val`` (the default is
infinity in both directions). ``step_sizes`` determines how far to go up or down when changing the stimulus value
adaptively. If it is a list of values, then the first element is used until the first reversal, the second until the
second reversal, etc. ``step_type`` determines what kind of steps are taken: 'lin' adds/subtracts the step size from
the current stimulus value, 'db' and 'log' will step by a certain number of decibels or log units.
Typically you would start with a large step size to quickly get close to the threshold, and then switch to a smaller
step size. Steps going up are multiplied with ``step_up_factor`` to allow unequal step sizes and weighted up-down
procedures (`Kaernbach (1991) <https://pubmed.ncbi.nlm.nih.gov/2011460/>`_).
Optimal step sizes are a bit smaller than the spread of the psychometric function for the parameter you are testing.
You can set the number of correct responses required to reduce the stimulus value with ``ndown`` and the number of
incorrect responses required to increase the value with ``nup``. The default is a 1up-2down procedure.
You can also add a number of training trials, in which the stimulus value does not change, with ``n_pretrials``.

.. _simulating:

Simulating responses
^^^^^^^^^^^^^^^^^^^^
For testing and comparing different staircase settings it can be useful to simulate responses. The first staircase
example uses :meth:`.simulate_responses` to draw responses from a logistic psychometric function with a given threshold
and width (expressed as the stimulus range in which the function increases from 20% to 80% hitrate).
For instance, if the current stimulus value is at the threshold, then the function returns a hit with 50% probability.
This is useful to simulate and compare different staircase settings and determine to which hit rate they converge.
For instance, let's get a feeling for the effect of the length of the measurement (number of reversals required to
end the staircase) and the accuracy of the threshold (standard deviation of thresholds across 100 simulated runs).
We test from 10 to 40 reversals and run 100 staircases in the inner loop, each time saving the threshold,
then computing the interquartile range and plotting it against the number of reversals. Longer measurements
should reduce the variability:

.. plot::
    :include-source:

    from matplotlib import pyplot as plt
    stairs_iqr =[]
    for reversals in range(10,41,5):
        threshs = []
        for _ in range(100):
            stairs = slab.Staircase(start_val=10, n_reversals=reversals)
            for trial in stairs:
                resp = stairs.simulate_response(3)
                stairs.add_response(resp)
            threshs.append(stairs.threshold())
        threshs.sort()
        stairs_iqr.append(threshs[74] - threshs[24]) # 75th-25th percentile
    plt.plot(range(10,41,5), stairs_iqr)
    plt.gca().set(xlabel='reversals', ylabel='threshold IQR')

Many other useful simulations are possible. You could check whether a 1up-3down procedure procedure would arrive at a
similar accuracy in fewer trials, what the best step size for a given psychometric function is, or how much a wider than
expected psychometric function increases experimental time. Simulations are a good starting point, but the psychometric
function is a very simplistic model for human behaviour. Check the results with pilot data.

Simulation is also useful for finding the hitrate (or point on the psychometric function) that a staircase converges on
in cases that are difficult for calculate. For instance, it is not immediately obvious on what threshold a 1up-4down
staircase with step_up_factor 1.5 and a 3-alternative forced choice presentation converges on::

    import numpy
    threshs = []
    width = 2
    thresh = 3
    for _ in range(100):
        stairs = slab.Staircase(start_val=10, n_reversals=30, n_down=4, step_up_factor=1.5)
        for trial in stairs:
            resp = stairs.simulate_response(threshold=thresh, transition_width=width, intervals=3)
            stairs.add_response(resp)
        threshs.append(stairs.threshold())
    # now we have 100 thresholds, take mean and convert to equivalent hitrate:
    hitrate = 1 / (1 + numpy.exp(4 * (0.5/width)  * (thresh - numpy.mean(threshs))))
    print(hitrate)
    # 0.83

As you can see, even through the threshold in the response simulation is 3 (that is, the rate of correct responses is
> 0.5 above this value; how fast it increases from there depends on the transition_width), the mean threshold returned
from the procedure is over 4.5. The last line translates this value in relation to the width of the simulated
psychometric function into a hitrate of about 0.83.

.. _responses:

Acquiring key presses
---------------------
When you use a staircase in a listening experiment, you need to record responses from the participant, usually in the
form of button presses. The :meth:`~slab.psychoacoustics.key` context manager can record single button presses
from the computer keyboard (or an attached USB number pad), or via the key press event handler of a matplotlib figure,
or from a custom USB buttonbox. The input is selected by setting :data:`slab.psychoacoustics.input_method` to 'keyboard',
'buttonbox', or 'figure'. This allow you to test your code on your laptop and switch to button box input at the lab
computer by changing a single line of code. Getting a button press from the keyboard will clear your terminal while
waiting for the response, and restore it afterwards. The the lab, you may not want to use a keyboard, which can be
distracting. A simple response box with the required number of buttons can be constructed easily with an
Arduino-compatible micro-controller that can send key codes to the computer via USB. Check for a press of a button
attached to a digital input and send a string corresponding to the key code of the desired key followed by the Enter key.
If you use the :meth:`~Staicase.plot` method of the :class:`Staircase` class to show the progress of the test, you can
set the :data:`~slab.psychoacoustics.input_method` to 'figure' to get a keypress via the figure's key press event
handler.

The :meth:`~slab.psychoacoustics.key` method uses the key code of a button, rather than the string character it produces
when pressed. You can find the code of a key by calling Python's :func:`ord` function. For instance, `ord('y')` returns
121, the code of the 'y' key.

The :class:`Trialsequence` and :class:`Staircase` classes have two convenience methods to present tones and acquire a
response from the listener in one step: :meth:`present_tone_trial` and :meth:`present_afc_trial`. Both take a list of
key codes that are considered valid responses (:param:`key_codes`). The list defaults to the number keys from 1 to 9.
If you use any of these keys in :meth:`present_tone_trial`, then you just need to specify which of them is counted as a
correct response by setting the argument `correct_key_idx` to the list index that contains the correct key (instead of a
single index you can specify a list of indices if you want to count several keys as correct). In
:meth:`present_afc_trial`, the order of the keys in :param:`key_codes` should correspond to the keys that should be
pressed to indicate interval 1, 2, etc. In this case, the correct key is different in each trial, depending on the
interval that contains the target stimulus.

Here is an example of how to use the :meth:`~slab.psychoacoustics.key` in a staircase that finds the detection threshold
for a 500 Hz tone, after every trial you have to indicate whether you could or could not hear the sound by pressing "y"
for yes or any other button for no:

.. _detection_example:

::

    stimulus = slab.Sound.tone(duration=0.5)
    stairs = slab.Staircase(start_val=60, step_sizes=[10, 3])
    for level in stairs:
        stimulus.level = level
        stimulus.play()
        with slab.key('Press y for yes or n for no.') as key:
            response = key.getch()
        if response == 121:  # 121 is the unicode for the "y" key
            stairs.add_response(True) # initiates calculation of next stimulus value
        else:
            stairs.add_response(False)
    stairs.plot()
    stairs.threshold()

Note that slab is not optimal for measuring reaction times due to the timing uncertainties in the millisecond range
introduced by modern multi-tasking operating systems. If you are serious about reaction times, you should use an
external DSP device to ensure accurate timing. Ubiquitous in auditory research are the realtime processors from
Tucker-Davies Technologies (our module `freefield` module works with these devices).

Precomputed sounds
------------------
If you present white noise in an experiment, you probably do not want to play the exact same noise in each trial
('frozen' noise), but different random instances of noise. The :class:`Precomputed` class manages a list of
pre-generated stimuli, but behave like a single sound. You can pass a list of sounds, a function to generate sounds
together with an indication of how many you want, or a generator expression to initialize the :class:`Precomputed`
object. The object has a :meth:`~Precomputed.play` method that plays a random stimulus from the list (but never the
stimulus played just before), and remembers all previously played stimuli in the :attr:`sequence`. The
:class:`Precomputed` object can be saved to a zip file and loaded back later on::

    # generate 10 instances of pink noise::
    stims = slab.Precomputed(lambda: slab.Sound.pinknoise(), n=10)
    stims.play() # play a random instance
    stims.play() # play another one, guaranteed to be different from the previous one
    stims.sequence # the sequence of instances played so far
    stims.write('stims.zip') # save the sounds as zip file
    stims = slab.Precomputed.read('stims.zip') # reloads the file into a Precomputed object


Results files
-------------
In most experiments, the performance of the listener, experimental settings, the presented stimuli, and other
information need to be saved to disk during the experiment. The :class:`ResultsFile` class helps with several typical
functions of these files, like generating timestamps, creating the necessary folders, and ensuring that the file is
readable if the experiment is interrupted writing to the file after each trial. Information is written incrementally to
the file in single lines of JSON (a `JSON Lines <http://jsonlines.org>`_ file).

Set the folder that will hold results files from all participants for the experiment somewhere at the top of your script
with the :data:`.results_folder`. Then you can create a file by initializing a class instance with a subject name::

    subject_ID = 'MS01'
    slab.ResultsFile.results_folder = 'MyResults'
    file = slab.ResultsFile(subject='MS')
    print(file.name)
    file.write(subject_ID)

You can now use the :meth:`~ResultsFile.write` method to write any information to the file, to be precise, you can write
any object that can be converted to JSON, like strings, lists, or dictionaries. Numpy data types need to be converted to
python types. A numpy array can be converted to a list before saving by calling its :meth:`numpy.ndarray.tolist` method,
and numpy ints or floats need to be converted by calling their :meth:`~numpy.int64.item` method. You can try out what
the JSON representation of an item is by calling::

    import json
    import numpy
    a = 'a string'
    b = [1, 2, 3, 4]
    c = {'frequency': 500, 'duration': 1.5}
    d = numpy.array(b)
    for item in [a, b, c]:
        json.dumps(item)
    json.dumps(d.tolist())

:class:`Trialsequence` and :class:`Staircase` objects can pass their entire current state to the write method, which
makes it easy to save all settings and responses from these objects::

    trials = slab.Trialsequence(conditions=4, n_reps=10)
    file.write(trials, tag='trials')

The :meth:`~ResultsFile.write` method writes a dictionary with a single key-value pair, where the key is supplied as
``tag`` argument argument (default is a time stamp in the format '%Y-%m-%d-%H-%M-%S'), and the value is the
json-serialized data you want to save. The information can be read back from the file, either while the experiment is
running and you need to access a previously saved result (:meth:`~ResultsFile.read`), or for later data analysis (:meth:`ResultsFile.read_file`). Both methods can take a ``tag`` argument to extract all instances saved under that tag
in a list.

Configuration files
-------------------
Another recurring issue when implementing experiments is loading configuration settings from a text file. Experiments
sometimes use configuration files when experimenters (who might not by Python programmers) need to set parameters
without changing the code. The format is a plain text file with a variable assignment on each line, because it is meant
to be written and changed by humans. The function :func:`~slab.psychoacoustics.load_config` reads the text file and
return a :func:`~collections.namedtuple` with the variable names and values. If you have a text file with the following
content::

    samplerate = 32000
    pause_duration = 30
    speeds = [60,120,180]

you can make all variables available to your script as attributes of the named tuple object::

    conf = slab.load_config('example.txt')
    conf.speeds
    % [60, 120, 180]
Introduction
============

Overview
--------

In this documentation we do not aim at providing a comprehensive explanation of
every single slab function (a complete description can be found in the :ref:`reference` section).
Rather, we want to provide some guidance for you to start generating sounds and running experiments.

For starters, you should have a look at the :ref:`Sounds` section. There, you will learn how to
generate, manipulate and write/read Sounds in slab. Next, you should see the :ref:`Psychoacoustics`
section which is about generating trial sequences and running experiments. With these tools you can
already do plenty of things! For example...

The :ref:`Filters` section contains some more advanced, but powerful, methods for processing
digital signals. The :ref:`hrtfs` section describes the handling of head related transfer functions and
will only be relevant if you are interested in spatial audio.


Frequently Asked Questions
--------------------------

* **Where can I learn enough Python to use this module?**

You can find many free courses online. We usually point our students to `Google's Python class <https://developers.google.com/edu/python>`_. For those of you who prefer video, Coursera has two suitable courses: `Python for Everybody <https://www.coursera.org/learn/python>`_ and `An Introduction to Interactive Programming with Python <https://www.coursera.org/learn/interactive-python-1?trk=profile_certification_title>`_.
There are also courses specifically for sound and signal processing, for instance `this one <https://www.coursera.org/learn/audio-signal-processing>`_.


* **Which Python environment do you use in the lab?**

We recommend `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_, which bundles Python and the conda package manager and installs quickly. You can then install only the packages that you need for your work, like IPython, numpy, scipy, and matplotlib, with a single command::

    conda install ipython numpy scipy matplotlib

When programming we use the command line with IPython and the Atom text editor with a package for syntax highlighting. Some lab members use `PyCharm <https://www.jetbrains.com/pycharm/>`_ or `Spyder <https://www.spyder-ide.org>`_ as integrated development environments. We don't recommend IDEs for beginners, because in our experience, students tend to conflate the IDE with Python itself and develop programming habits that they need to unlearn when they want to get productive.


* **I get import errors when using certain functions!**

Slab requires additional modules for some functionality. These modules are not installed automatically because not everyone may need them (such as HRTF file reading) or the installation is OS-dependent (such as SoundFile and curses). Please see :ref:`installation` for how and what to install should you need it. The import error messages will in most cases give you the necessary installation command for Mac/Linux systems.


* **I have set the level of a sound to 70 dB but it is way louder, why?**

This is because slab does not know the hardware you are using to play sound. For example, white noise is generated so that the maximum value in the time series is +1 and the minimum minus one ("full scale"). The RMS of this signal, expressed in deciBels happens to be about 82 dB, but you need to calibrate your system (see :ref:`calibration`) so that the calculated intensity is meaningful. Relative intensities are correct without calibration---so decreasing the intensity by 10 dB (`sound.level -= 10`) will work as expected.


* **What is the difference between white noise and pink noise?**

White noise is a signal that consists of random numbers. This signal has equal power at all frequencies. However, our auditory system does not perceive it that way, which is why white noise appears high-pitched. In the pink noise signal, the power decreases with frequency to correct for this effect. Pink noise is thus a more appropriate choice for a masking or background noise, because it has the same power in each octave. However, there are even better options. The :meth:`~slab.Sound.erb_noise` method constructs a noise with equal energy not in octaves, but in fractions of approximated auditory filters widths (equivalent rectangular bandwidths, ERB). Or the :meth:`~slab.Sound.multitone_masker`, which is a noise-like combination of many pure tones at ERB intervals. This noise does not have random amplitude variations and masks evenly across frequency and time.


* **I think I found a bug!**

Please see the `bug reports <https://github.com/user/DrMarc/soundlab/CONTRIBUTING.md#bugs>`_ section in the contribution guidelines.


* **How can I contribute to the project?**

Please see the `pull request <https://github.com/user/DrMarc/soundlab/CONTRIBUTING.md#pull-requests>`_ section in the contribution guidelines if you want to contribute code or useful examples for the documentation.
.. _Sounds:

Sound
=====

Generating sounds
-----------------
The :class:`Sound` class provides methods for generating, manipulating, displaying, and analysing sound stimuli.
You can generate typical experimental stimuli with this class, including tones, noises, and click trains, and also
more specialized stimuli, like equally-masking noises, Schroeder-phase harmonics, iterated ripple noise and synthetic
vowels.
Slab methods assume sensible defaults where possible. You can call most methods without arguments to get an impression
of what they do (f.i. :meth:`slab.Sound.tone()` returns a 1s-long 1kHz tone at 70 dB sampled at 8 kHz) and then
customise from there.
For instance, let's make a 500 ms long 500 Hz pure tone signal with a band-limited (one octave below and above
the tone) pink noise background with a 10 dB signal-to-noise ratio: ::

  tone = slab.Sound.tone(frequency=500, duration=0.5)
  tone.level = 80 # setting the intensity to 80 dB
  noise = slab.Sound.pinknoise(duration=0.5)
  noise.filter(frequency=(250, 1000), kind='bp') # bandpass .25 to 1 kHz
  noise.level = 70 # 10 dB lower than the tone
  stimulus = tone + noise # combine the two signals
  stimulus = stimulus.ramp() # apply on- and offset ramps to avoid clicks
  stimulus.play()

:class:`Sound` objects have many useful methods for manipulating (like :meth:`.ramp`, :meth:`.filter`,
and :meth:`.pulse`) or inspecting them (like :meth:`.waveform`, :meth:`.spectrum`, and :meth:`.spectral_feature`).
A complete list is in the :ref:`Reference` section, and the majority is also discussed here. If you use IPython,
you can tap the `tab` key after typing ``slab.Sound.``, or the name of any Sound object followed by a full stop,
to get an interactive list of the possibilities.

Sounds can also be created by recording them with :meth:`slab.Sound.record`. For instance
``recording = slab.Sound.record(duration=1.0, samplerate=44100)`` will record a 1-second sound at 44100 Hz from the
default audio input (usually the microphone). The ``record`` method uses
`SoundCard <https://github.com/bastibe/SoundCard>`_ if installed, or `SoX <http://sox.sourceforge.net>`_
(via a temporary file) otherwise. Both are cross-platform and easy to install. If neither tool is installed,
you won't be able to record sounds.

Specifying durations
--------------------
Sometimes it is useful to specify the duration of a stimulus in samples rather than seconds. All methods that generate
sounds have a :attr:`duration` argument that accepts floating point numbers or integers. Floating point numbers are
interpreted as durations in seconds (``slab.Sound.tone(duration=1.0)`` results in a 1 second tone). Integers are
interpreted as number of samples (``slab.Sound.tone(duration=1000)`` gives you 1000 samples of a tone).

Setting the sample rate
-----------------------
We did not specify a sample rate for any of the stimuli in the examples above. When the :attr:`samplerate` argument of
a sound-generating method is not specified, the default sample rate (8 kHz if not set otherwise) is used. It is possible
to set a sample rate separately for each Sound object, but it is usually better to set a suitable default sample rate
at the start of your script or Python session using :func:`slab.set_default_samplerate`. This rate is kept in the class
variable :data:`_default_samplerate` and is used whenever you call a sound generating method without specifying a rate.
This rate depends on the frequency content of your stimuli and should be at least double the highest frequency of
interest. For some speech sounds or narrow bad noises you might get away with 8 kHz; for spatial sounds you may need 48
kHz or more.

Specifying levels
--------------------------
Same as for the sample rate, sounds are generated at a default level (70 dB if not set otherwise). The default is kept
in the class variable :data:`_default_level` and you can set set it to a different value using
:func:`slab.set_default_level`. Level are not specified directly when generating sounds, but rather afterwards by
setting the :attr:`level` property::

    sig = slab.Sound.pinknoise()
    sig.level # return the current level
    sig.level = 85 # set a new level

Note that the returned level will *not* be the actual physical playback level, because that depends on the playback
hardware (soundcard, amplifiers, headphones, speakers). Calibrate your system if you need to play stimuli at a known
level (see :ref:`calibration`).

.. _calibration:

Calibrating the output
----------------------
Analogous to setting the default level at which sounds are generated with ``slab.set_default_level()``. Each sound's
level can be set individually by changing its :attr:`level` property. Setting the :attr:`level` property of a
stimulus changes the root-mean-square of the waveform and relative changes are correct (reducing the level attribute by
10 dB will reduce the sound output by the same amount), but the *absolute* intensity is only correct if you calibrate
your output. The recommended procedure it to set your system volume to maximum, connect the listening hardware
(headphone or loudspeaker) and set up a sound level meter. Then call :func:`slab.calibrate`. The :func:`.calibrate`
function will play a 1 kHz tone for 5 seconds. Note the recorded intensity on the meter and enter it when requested. The
function returns a calibration intensity, i.e. difference between the tone's level attribute and the recorded level.
Pass this value to :func:`slab.set_calibration_intensity` to to correct the intensities returned by the :attr:`level`
property all sounds. The calibration intensity is saved in the class variable :data:`_calibration_intensity`.
It is applied to all level calculations so that a sound's level attribute now roughly corresponds to the actual output
intensity in dB SPL---'roughly' because your output hardware may not have a flat frequency transfer function
(some frequencies play louder than others). See :ref:`Filters` for methods to equalize transfer functions.

Experiments sometimes require you to play different stimuli at comparable loudness. Loudness is the perception of sound
intensity and it is difficult to calculate. You can use the :meth:`Sound.aweight` method of a sound to filter it so that
frequencies are weighted according to the typical human hearing thresholds. This will increase the correspondence
between the rms intensity measure returned by the :attr:`level` attribute and the perceived loudness. However, in most
cases, controlling relative intensities is sufficient.

To increase the accuracy of the calibration for your experimental stimuli, pass a sound with a similar spectrum to
:func:`slab.calibrate`. For instance, if your stimuli are wide band pink noises, then you may want to use a pink noise
for calibration. The `level` of the noise should be high, but not cause clipping.

If you do not have a sound level meter, then you can present sounds in dB HL (hearing level). For that, measure the
hearing threshold of the listener at the frequency or frequencies that are presented in your experiment and play your
stimuli at a set level above that threshold. You can measure the hearing threshold at one frequency (or for any
broadband sound) with the few lines of code (see :ref:`audiogram`).

Saving and loading sounds
-------------------------
You can save sounds to wav files by calling the object's :meth:`.Sound.write` method (``signal.write('signal.wav')``).
By default, sounds are normalized to have a maximal amplitude of 1 to avoid clipping when writing the file.
You should set :attr:`signal.level` to the intended level when loading a sound from file or disable normalization
if you know what you are doing. You can load a wav file by initializing a Sound object with the filename:
``signal = slab.Sound('signal.wav')``.

Combining sounds
----------------
Several functions allow you to string stimuli together. For instance, in a forward masking experiment [#f1]_ we need a
masking noise followed by a target sound after a brief silent interval. An example implementation of a complete
experiment is discussed in the :ref:`Psychoacoustics` section, but here, we will construct the stimulus: ::

    masker = slab.Sound.tone(frequency=550, duration=0.5) # a 0.5s 550 Hz tone
    masker.level = 80 # at 80 dB
    masker.ramp() # default 10 ms raised cosine ramps
    silence = slab.Sound.silence(duration=0.01) # 10 ms silence
    signal = slab.Sound.tone(duration=0.05) # using the default 500 Hz
    signal.level = 80 # let's start at the same intensity as the masker
    signal.ramp(duration=0.005) # short signal, we'll use 5 ms ramps
    stimulus = slab.Sound.sequence(masker, silence, signal)
    stimulus.play()

We can make a classic non-interactive demonstration of forward masking by playing these stimuli with decreasing signal
level in a loop, once without the masker, and once with the masker.
Count for how many steps you can hear the signal tone: ::

    import time # we need the sleep function
    for level in range(80, 10, -5): # down from 80 in steps of 5 dB
        signal.level = level
        signal.play()
        time.sleep(0.5)
    # now with the masker
    for level in range(80, 10, -5): # down from 80 in steps of 5 dB
        signal.level = level
        stimulus = slab.Sound.sequence(masker, silence, signal)
        stimulus.play()
        time.sleep(0.5)

Many listeners can hear all of the steps without the masker, but only the first 6 or 7 steps with the masker. This
depends on the intensity at which you play the demo (see :ref:`Calibrating the output<calibration>` below).
The :meth:`.sequence` method is an example of list unpacking---you can provide any number of sounds to be concatenated.
If you have a list of sounds, call the method like so: ``slab.Sound.sequence(*[list_of_sound_objects])``
to unpack the list into function arguments.

Another method to put sounds together is :meth:`.crossfade`, which applies a crossfading between two sounds with a
specified :attr:`overlap` in seconds. An interesting experimental use is in adaptation designs, in which one longer
stimulus is played to adapt neuronal responses to its sound features, and then a new stimulus feature is introduced
(but nothing else changes). Responses (measured for instance with EEG) at that point will be mostly due to that feature.
A classical example is the pitch onset response, which is evoked when the temporal fine structure of a continuous noise
is regularized to produce a pitch percept without altering the sound spectrum
(see `Krumbholz et al. (2003) <https://pubmed.ncbi.nlm.nih.gov/12816892/>`_).
It is easy to generate the main stimulus of that study, a noise transitioning to an iterates ripple noise after two
seconds, with 5 ms crossfade overlap, then filtered between 0.8 and 3.2 kHz: ::

    slab.set_default_samplerate(16000) # we need a higher sample rate
    slab.set_default_level(80)  # set the level for all sounds to 80 dB
    adapter = slab.Sound.whitenoise(duration=2.0)
    irn = slab.Sound.irn(frequency=125, n_iter=2, duration=1.0) # pitched sound
    stimulus = slab.Sound.crossfade(adapter, irn, overlap=0.005) # crossfade
    stimulus.filter(frequency=[800, 3200], kind='bp') # filter
    stimulus.ramp(duration=0.005) # 5 ms on- and offset ramps
    stimulus.spectrogram() # note that there is no change at the transition
    stimulus.play() # but you can hear the onset of the regularity (pitch)

Plotting and analysis
---------------------
You can inspect sounds by plotting the :meth:`.waveform`, :meth:`.spectrum`, or :meth:`.spectrogram`:

.. plot::
    :include-source:

    from matplotlib import pyplot as plt
    a = slab.Sound.vowel(vowel='a')
    e = slab.Sound.vowel(vowel='e')
    i = slab.Sound.vowel(vowel='i')
    signal = slab.Sound.sequence(a,e,i)
    import matplotlib.pyplot as plt # preparing a 2-by-2 figure
    _, [[ax1, ax2], [ax3, ax4]] = plt.subplots(
                    nrows=2, ncols=2, constrained_layout=True)
    signal.waveform(axis=ax1, show=False)
    signal.waveform(end=0.05, axis=ax2, show=False) # first 50ms
    signal.spectrogram(upper_frequency=5000, axis=ax3, show=False)
    signal.spectrum(axis=ax4)

Instead of plotting, :meth:`.spectrum` and :meth:`.spectrogram` will return the time frequency bins and spectral power
values for further analysis if you set the :attr:`show` argument to False. All plotting functions can draw into an
existing matplotlib.pyplot axis supplied with the :attr:`axis` argument.

.. _spectral_features:

You can also extract common features from sounds, such as the :meth:`.crest_factor` (a measure of how 'peaky'
the waveform is), or the average :meth:`.onset_slope` (a measure of how fast the on-ramps in the sound are---important
for sound localization). Features of the spectral content are bundled in the :meth:`.spectral_feature` method.
It can compute spectral centroid, flux, flatness, and rolloff, either for an entire sound (suitable for stationary
sounds), or for successive time windows (frames, suitable for time-varying sounds).
* The centroid is a measure of the center of mass of a spectrum (i.e. the 'center' frequency).
* The flux measures how quickly the power spectrum is changing by comparing the power spectrum for one frame against the
power spectrum from the previous frame; flatness measures how tone-like a sound is, as opposed to being noise-like, and
is calculated by dividing the geometric mean of the power spectrum by the arithmetic mean (see `Dubnov (2004) <https://ieeexplore.ieee.org/document/1316889>`_).
* The rolloff measures the frequency at which the spectrum rolls off, typically used to find a suitable low-cutoff
frequency that retains most of the sound power.
These particular features are integrated in slab because we find them useful in our daily work. Many more features are
available in packages specialised on audio processing, for instance `librosa <https://librosa.org>`_. librosa interfaces
easily with slab, you can just hand the sample data and the sample rate of an slab object separately to most of its
methods::

    import librosa
    sig = slab.Sound('music.wav') # load wav file into slab.Sound object
    librosa.beat.beat_track(y=sig.data, sr=sig.samplerate)

When working with environmental sounds or other recorded stimuli, one often needs to compute relevant features for
collections of recordings in different experimental conditions. The slab module contains a function
:func:`slab.apply_to_path`, which applies a function to all sound files in a given folder and returns a dictionary of file
names and computed features. In fact, you can also use that function to modify (for instance ramp and filter) all files
in a folder.

For other time-frequency processing, the :meth:`.frames` provides an easy way to step through the signal in short
windowed frames and compute some values from it. For instance, you could detect on- and offsets in the signal
by computing the crest factor in each frame: ::

    from matplotlib import pyplot as plt
    signal.pulse() # apply a 4 Hz pulse to the 3 vowels from above
    signal.waveform() # note the pulses
    crest = [] # the short-term crest factor will show on- and offsets
    frames = signal.frames(duration=64)
    for f in frames:
        crest.append(f.crest_factor())
    times = signal.frametimes(duration=64) # frame center times
    import matplotlib.pyplot as plt
    plt.plot(times, crest) # peaks in the crest factor mark intensity ramps

Binaural sounds
---------------
For experiments in spatial hearing, or any other situation that requires differential manipulation of the left and
right channel of a sound, you can use the :class:`Binaural` class. It inherits all methods from :class:`Sound` and
provides additional methods for generating and manipulating binaural sounds, including advanced interaural time
and intensity manipulation.

Generating binaural sounds
^^^^^^^^^^^^^^^^^^^^^^^^^^
Binaural sounds support all sound generating functions with a :attr:`n_hannels` attribute of the :class:`Sound` class,
but automatically set :attr:`n_channels` to 2. Noises support an additional :attr:`kind` argument,
which can be set to 'diotic' (identical noise in both channels) or 'dichotic' (uncorrelated noise). Other methods just
return 2-channel versions of the stimuli. You can recast any Sound object as Binaural sound, which duplicates the first
channel if :attr:`n_channels` is 1 or greater than 2: ::

    monaural = slab.Sound.tone()
    monaural.n_channels
    out: 1
    binaural = slab.Binaural(monaural)
    binaural.n_channels
    out: 2
    binaural.left # access to the left channel
    binaural.right # access to the right channel

Loading a wav file with ``slab.Binaural('file.wav')`` returns a Binaural sound object with two channels (even if the
wav file contains only one channel).

Manipulating ITD and ILD
^^^^^^^^^^^^^^^^^^^^^^^^
The easiest manipulation of a binaural parameter may be to change the interaural level difference (ILD).
This can be achieved by setting the :attr:`level` attributes of both channels: ::

    noise = slab.Binaural.pinknoise()
    noise.left.level = 75
    noise.right.level = 85
    noise.level
    out: array([75., 85.])

The :meth:`.ild` makes this easier and keeps the overall level constant: ``noise.ild(10)`` amplifies the right channel
by 5 dB and attenuates the left channel by the same amount to achieve a 10dB level difference. Positive dB values
move the virtual sound source to the right and negative values move the source to the left. The pink noise in the
example is a broadband signal, and the ILD is frequency dependent and should not be the same for all frequencies. A
frequency-dependent level difference can be computed and applied with :meth:`.interaural_level_spectrum`. The level
spectrum is computed from a head-related transfer function (HRTF) and can be customised for individual listeners.
See :ref:`hrtfs` for how to handle these functions. The default level spectrum is computed form the HRTF of the KEMAR
binaural recording mannequin (as measured by
`Gardener and Martin (1994) <https://sound.media.mit.edu/resources/KEMAR.html>`_ at the MIT Media Lab).
The level spectrum takes a while to compute and it may be useful to save it. It is a Python dict containing the level
differences in a numpy array along with a frequency vector, an azimuth vector, and the sample rate. You can save it for
instance with pickle: ::

    import pickle
    ils = slab.Binaural.make_interaural_level_spectrum()
    pickle.dump(ils, open('ils.pickle', 'wb')) # save using pickle
    ils = pickle.load(open('ils.pickle', 'rb')) # load pickle

If the limitations of pickle worry you, you can use numpy.save with a small caveat when loading: numpy.save wraps the
dict in an object and we need to remove that after loading with the somewhat strange index `[()]`: ::

    import numpy
    numpy.save('ils.npy', ils) # save using numpy
    ils = numpy.load('ils.npy, allow_pickle=True)[()] # load and get the original dict from the wrapping object

If you are unsure which ILD value is appropriate, :meth:`.azimuth_to_ild` can compute ILDs corresponding to an azimuth
angle, for instance 45 degrees, and a frequency: ::

    slab.Binaural.azimuth_to_ild(45)
    # -9.12  # correct ILD in dB
    noise.ild(-9.12)  # apply the ILD

A dynamic ILD, which evokes the perception of a moving sound source, can be applied with
:meth:`.ild_ramp`. The ramp is linear from and to a given ILD.

Similar functions exist to manipulate interaural time differences (ITD): :meth:`.itd`, :meth:`.azimuth_to_ild`
(using a given head radius), and :meth:`.itd_ramp`. To present a signal from a given azimuth using both cues,
use the :meth:`.at_azimuth`, which calculates the correct ILD and ITD for you and applies it.

ITD and ILD manipulation leads to the percept of *lateralization*, that is, a source somewhere between the
ears inside the head. Additional spectral shaping is necessary to generate an externalized percept (outside the head).
This shaping can be achieved with the :meth:`.externalize`, which applies a low-resolution HRTF filter
(KEMAR by default). Using both ramp functions and externalization, it is easy to generate a convincing sound source
movement with pulsed pink noise: ::

    noise = slab.Binaural.pinknoise(samplerate=44100)
    from_ild = slab.Binaural.azimuth_to_ild(-90)
    from_itd = slab.Binaural.azimuth_to_itd(-90)
    to_ild = slab.Binaural.azimuth_to_ild(90)
    to_itd = slab.Binaural.azimuth_to_itd(90)
    noise_moving = noise.ild_ramp(from_ild, to_ild)
    noise_moving = noise_moving.itd_ramp(from_itd, to_itd)
    noise_moving.externalize() # apply filter in place
    noise_moving.play() # best through headphones


Signals
-------
Sounds inherit from the :class:`Signal` class, which provides a generic signal object with properties duration,
number of samples, sample times, number of channels. The actual samples are kept as numpy array in the :attr:`data`
property and can be accessed, if necessary as for instance :attr:`signal.data`. Signals support slicing, arithmetic
operations, and conversion between sample points and time points directly, without having to access the :attr:`data`
property. The methods :meth:`.resample`, :meth:`.envelope`, and :meth:`.delay` are also implemented in Signal and
passed to the child classes :class:`Sound`, :class:`Binaural`, and :class:`Filter`. You do not normally need to use
the Signal class directly. ::

    sig = slab.Sound.pinknoise(n_channels=3)
    sig.duration
    out: 1.0
    sig.n_samples
    out: 8000
    sig.data.shape # accessing the sample array
    out: (8000, 3) # which has shape (n_samples x n_channels)
    sig2 = sig.resample(samplerate=4000) # resample to 4 kHz
    env = sig2.envelope() # returns a new signal containing the lowpass Hilbert envelopes of both channels
    sig.delay(duration=0.0006, channel=0) # delay the first channel by 0.6 ms

.. rubric:: Footnotes

.. [#f1] Forward masking occurs when a signal cannot be heard due to a preceding masking sound. Typically, three intervals are presented to the listener, two contain only the masker and one contains the masker followed by the signal. The listener has to identify the interval with the signal. The level of the masker is fixed and the signal level is varied adaptively to obtain the masked threshold.
