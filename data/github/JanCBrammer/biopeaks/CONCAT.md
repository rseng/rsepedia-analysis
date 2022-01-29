<img src="https://github.com/JanCBrammer/biopeaks/raw/master/docs/images/logo.png" alt="logo" style="width:600px;"/>

![GH Actions](https://github.com/JanCBrammer/biopeaks/workflows/CI/badge.svg?branch=master)
[![codecov](https://codecov.io/gh/JanCBrammer/biopeaks/branch/master/graph/badge.svg)](https://codecov.io/gh/JanCBrammer/biopeaks)
[![DOI](https://www.zenodo.org/badge/172897525.svg)](https://www.zenodo.org/badge/latestdoi/172897525)
[![PyPI version](https://img.shields.io/pypi/v/biopeaks.svg)](https://pypi.org/project/biopeaks/)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.02621/status.svg)](https://doi.org/10.21105/joss.02621)


# General Information

`biopeaks` is a straightforward graphical user interface for feature extraction from electrocardiogram (ECG), photoplethysmogram (PPG) and breathing biosignals.
It processes these biosignals semi-automatically with sensible defaults and offers the following functionality:

+ processes files in the open biosignal formats [EDF](https://en.wikipedia.org/wiki/European_Data_Format), [OpenSignals (Bitalino)](https://bitalino.com/en/software)
as well as plain text files (.txt, .csv, .tsv)
+ interactive biosignal visualization
+ biosignal segmentation
+ benchmarked, automatic extrema detection (R-peaks in ECG, systolic peaks in PPG, exhalation troughs and inhalation
peaks in breathing signals) with signal-specific, sensible defaults
+ automatic state-of-the-art [artifact correction](https://www.tandfonline.com/doi/full/10.1080/03091902.2019.1640306)
 for ECG and PPG extrema
+ manual editing of extrema
+ extraction of instantaneous features: (heart- or breathing-) rate and period, as well as breathing amplitude
+ .csv export of extrema and instantaneous features for further analysis (e.g., heart rate variability)
+ automatic analysis of multiple files (batch processing)


![GUI](https://github.com/JanCBrammer/biopeaks/raw/master/docs/images/screenshot_statistics.png)


# Installation

`biopeaks` can be installed from PyPI:

```
pip install biopeaks
```

Alternatively, on Windows, download [biopeaks.exe](https://github.com/JanCBrammer/biopeaks/releases/latest)
and run it. Running the executable does not require a Python installation.

You can find more details on the installation [here](https://jancbrammer.github.io/biopeaks/installation.html).


# Documentation

Have a look at the [user guide](https://jancbrammer.github.io/biopeaks/user_guide.html) to get started with `biopeaks`.


# Contributors welcome!

Improvements or additions to the repository (documentation, tests, code) are welcome and encouraged.
Spotted a typo in the documentation? Caught a bug in the code? Ideas for improving the documentation,
increase test coverage, or adding features to the GUI? Get started with the [contributor guide](https://jancbrammer.github.io/biopeaks/contributor_guide.html).


# Citation

Please refer to the [biopeaks paper](https://joss.theoj.org/papers/10.21105/joss.02621) in The Journal of Open Source Software.


# Changelog

Have a look at the [changelog](https://jancbrammer.github.io/biopeaks/changelog.html) to get an overview of what has changed throughout the versions of `biopeaks`.





# Contributor Covenant Code of Conduct

## Our Pledge

Contributors pledge to make participation in the `biopeaks` community a
harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse and inclusive community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
jan.c.brammer@gmail.com.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
---
title: 'biopeaks: a graphical user interface for feature extraction from heart- and breathing biosignals'
tags:
  - Python
  - GUI
  - biosignals
  - heart
  - breathing
  - PPG
  - ECG
  - feature extraction
authors:
  - name: Jan C. Brammer
    orcid: 0000-0002-7664-3753
    affiliation: 1
affiliations:
  - name: Behavioral Science Institute, Radboud University Nijmegen, Nijmegen, The Netherlands
    index: 1
date: 23 August 2020
bibliography: paper.bib

---


# Statement of need

Heart- and breathing biosignals increasingly gain popularity in academia and industry, sparked by the availability of
easy-to-use, low-cost biosignal sensors along with a growing ecosystem of free,
open-source software for the analysis of biosignals [@neurokit; @biosppy]. However, an open-source, freely
available graphical user interface (GUI) for biosignal analysis is currenly lacking. `biopeaks` addresses this need.
Compared to application programming interfaces [@neurokit; @biosppy], its GUI allows for
more intuitive and immediate visual interaction with the biosignal, which is especially valuable
during data preprocessing and exploration. At the time of writing, `biopeaks` is used in multiple projects at the [Gemhlab](https://gemhlab.com/).


# Functionality

+ processing of open biosignal formats EDF [@edf], OpenSignals [@opensignals], as well as plain text files (.txt, .csv, .tsv)
+ interactive biosignal visualization
+ biosignal segmentation
+ automatic extrema detection (R-peaks in electrocardiogram (ECG), systolic peaks in photoplethysmogram (PPG), as well as exhalation troughs and inhalation peaks in breathing biosignals)
with biosignal-specific, sensible defaults
+ automatic state-of-the-art artifact correction for ECG and PPG extrema
+ point-and-click editing of extrema
+ extraction of instantaneous features: rate and period for heart and breathing biosignals, as well as breathing amplitude
+ automatic analysis of multiple files (batch processing)
+ .csv export of extrema and instantaneous features for further analysis (e.g., heart rate variability)

An analyst who wants to extract information from heart or breathing biosignals performs multiple analysis steps.
First, they verify if the biosignal's quality is sufficient for analysis, since biosignals can be corrupted
for a number of reasons (e.g., movement artifacts, poor sensor placement). `biopeaks` allows
the analyst to quickly visualize a biosignal and interact with it (panning, zooming) to evaluate its quality.
If the analyst deems the biosignal's quality sufficient, they proceed to identify local extrema in the physiological time series.
Local extrema include R-peaks in ECG and systolic peaks in PPG, representing
the contraction of the ventricular heart muscle. In breathing biosignals,
the relevant local extrema are inhalation peaks and exhalation troughs. `biopeaks` detects these extrema automatically
using three biosignal-specific algorithms. Breathing extrema are detected using a variant of the "zero-crossing algorithm
with amplitude threshold" [@khodadad]. Systolic peaks in PPG signals are identified using an implementation of "Method IV;
Event-Related Moving Averages with Dynamic Threshold" introduced by Elgendi et al. [@elgendi]. Lastly, the ECG R-peak detector is a
custom algorithm that has been evaluated on the Glasgow University Database (GUDB) [@gudb] which contains ECG signals along with R-peak annotations. The performance of the R-peak detector has been evaluated in terms of sensitivity (aka recall; i.e., how many of the correct extrema were detected?) and precision (i.e., how many of the detected extrema are correct extrema?). Peak detection has been evaluated on the records of all 25 participants included in the GUDB using the ECG channel corresponding to Einthoven lead II. The tolerance for true positive peak detection was set to one sample. The GUDB has not been used to optimize the R-peak detector prior to the performance evaluation. The performance at rest (sitting, 25 records) and in dynamic conditions (handbike, 24 records due to the missing R-peak annotations of participant 04) is as follows:

|           |    |sitting|handbike|
|:---------:|:--:|:-----:|:------:|
|precision  |mean|.9995  |.9855   |
|           |std |.0017  |.0234   |
|sensitivity|mean|.9974  |.9853   |
|           |std |.0037  |.0250   |

The code for performance evaluation is included in the `biopeaks` installation and can be run without downloading the GUDB (the database is streamed).
Despite the robust performance of the extrema detectors, algorithmically identified extrema can be misplaced (false positives) or extrema might be missed (false negatives),
if there are noisy segments in the biosignal. If left uncorrected, these errors can significantly distort subsequent analysis steps [@berntson]. To address this problem and ensure the correct placement of extrema, `biopeaks` offers intuitive
point-and-click extrema editing (i.e., removing and adding extrema). Additionally, for cardiac biosignals,
`biopeaks` offers state-of-the-art automatic extrema correction [@lipponen]. Finally, based on the local extrema, the analyst can extract features
from the biosignal. The features are based on temporal or amplitude differences between the local extrema.
\autoref{Figure 1}, \autoref{Figure 2}, and \autoref{Figure 3} illustrate the extraction of instantaneous heart period, breathing period, and breathing (inhalation) amplitude respectively.

![Extraction of heart period (panel B) based on R-peaks in an ECG (panel A). Note that this is conceptually identical to the extraction of heart period based on systolic peaks in PPG.\label{Figure 1}](fig_heartperiod.png)

![Extraction of breathing period (panel B) based on inhalation peaks in a breathing biosignal (panel A).\label{Figure 2}](fig_breathingperiod.png)

![Extraction of inhalation amplitude (panel B) based on breathing extrema in a breathing biosignal (panel A).\label{Figure 3}](fig_breathingamplitude.png)

In summary, `biopeaks` is designed to make heart- and breathing biosignal inspection, extrema detection and editing, as well as feature
extraction fast and intuitive. It does not aim to offer a suite of low-level signal processing tools (e.g., digital filter design),
and abstracts these details away with sensible biosignal-specific defaults. `biopeaks` is implemented in Python using the cross-platform
PySide2 framework (official Python bindings for Qt) [@pyside2], and leverages Matplotlib [@matplotlib], NumPy [@numpy], SciPy [@scipy] and pandas [@pandas] for visualization and signal processing.
There are freely available alternatives to `biopeaks` that are implemented in MATLAB [@artiifact; @physiodatatoolbox] or C# [@signalplant].
However, the source code of these tools is not available [@artiifact; @signalplant] or they are not released under an open
source license [@artiifact; @physiodatatoolbox; @signalplant].


# Acknowledgements
Nastasia Griffioen, Babet Halberstadt, Joanneke Weerdmeester, and Abele Michela provided invaluable feedback throughout the development of `biopeaks`.


# References
# Changelog

### Version 1.4.4 (January 05, 2022)
+ enhancement: ported from PySide2 to PySide6.

### Version 1.4.3 (January 03, 2022)
Removed due to build error in Windows executable (also yanked from PyPI).

### Version 1.4.2 (June 07, 2021)
+ enhancement: [improved stopping criterion for iterative artifact correction](https://github.com/JanCBrammer/biopeaks/commit/3a25e7c4f9cef3cab28afe067449f280340e71ee). 
+ enhancement: [using sos format instead of ba format for butterworth filters.](https://github.com/JanCBrammer/biopeaks/commit/8f52909cebafd3b162c943dddf9e4ca1d8838cab). 
  
### Version 1.4.1 (October 26, 2020)
+ enhancement: updated documentation and docstrings.
+ enhancement: more convenient PPG benchmarking.
+ bugfix: `resp.ensure_peak_trough_alternation()` now considers neighboring extrema with equal amplitude.

### Version 1.4.0 (August 04, 2020)
+ enhancement: added support for plain text files (.txt, .csv, .tsv).
+ enhancement: stream [Glasgow University Database (GUDB)](http://researchdata.gla.ac.uk/716/) for ECG benchmarking (download is no longer required).

### Version 1.3.2 (June 07, 2020)
+ enhancement: visibility of configuration panel can now be toggled (more screen space for signals).
+ bugfix: fixed index-out-of-range error in `heart._correct_misaligned()`.

### Version 1.3.1 (May 06, 2020)
+ bugfix: corrected control-flow of artifact classification during auto-correction of ECG and PPG peaks.
+ bugfix: using PySide2 resources instead of PyQt5 resources

### Version 1.3.0 (April 27, 2020)
+ enhancement: using the official Qt Python bindings (PySide2) instead of PyQt5.
+ bugfix: enforcing minimal figure height to comply with requirements of Matplotlib > 3.2.0.

### Version 1.2.2 (April 15, 2020)
+ enhancement: faster auto-correction of ECG and PPG peaks with pandas rolling window.
+ bugfix: corrected error in calculation of subspace 22 during auto-correction of ECG and PPG peaks.
+ bugfix: during segmentation, extrema at segment boundaries are now excluded from segment.

### Version 1.2.1 (April 01, 2020)
+ enhancement: auto-correction of ECG and PPG peaks is now optional (instead
of being applied by default during the calculation of the statistics).
+ enhancement: improved baseline removal for respiration signals.

### Version 1.2.0 (March 20, 2020)
+ enhancement: added peak detection for photoplethysmogram (PPG) (heart.ppg_peaks()), based on
[Elgendi et al., (2013)](https://journals.plos.org/plosone/article/comments?id=10.1371/journal.pone.0076585).
The performance of heart.ppg_peaks() has been evaluated on the PPG signals of all
42 subjects in the [Capnobase IEEE TBME benchmark dataset](http://www.capnobase.org/index.php?id=857).
The dataset has not been used to optimize `heart.ppg_peaks()` in any way prior to
the performance evaluation. The tolerance for peak detection was set to 50 milliseconds in
accordance with [Elgendi et al., (2013)](https://journals.plos.org/plosone/article/comments?id=10.1371/journal.pone.0076585).

|metric     |summary|version 1.2.0
|:---------:|:-----:|:-----------:
|precision  |mean   |.996
|           |std    |.004
|sensitivity|mean   |.999
|           |std    |.001

+ bugfix: the PATCH version has been reset to 0 after incrementing MINOR version (https://semver.org/)

### Version 1.1.6 (March 06, 2020)
+ enhancement: some small improvements of the statistics panel in **datadisplay**.

### Version 1.1.5 (March 01, 2020)
+ enhancement: added support for [EDF files](https://en.wikipedia.org/wiki/European_Data_Format).
+ enhancement: `ecg.ecg_peaks()` now filters out power-line noise at 50Hz. This
further increases the performance on the [Glasgow University Database (GUDB)](http://researchdata.gla.ac.uk/716/).
Again, the GUBD has not been used to optimize `ecg.ecg_peaks()` in any way prior to
the performance evaluation. The tolerance for peak detection was set to one
sample.

|condition|metric     |summary|version 1.0.2|version 1.0.3|version 1.1.5
|:-------:|:---------:|:-----:|:-----------:|:-----------:|:-----------:
|sitting  |precision  |mean   |.999         |.998         |.998
|         |           |std    |.002         |.005         |.002
|         |sensitivity|mean   |.996         |.996         |.998
|         |           |std    |.008         |.004         |.004
|handbike |precision  |mean   |.904         |.930         |.984
|         |           |std    |.135         |.127         |.022
|         |sensitivity|mean   |.789         |.857         |.984
|         |           |std    |.281         |.247         |.025

### Version 1.0.5 (February 09, 2020)
+ enhancement: improved ECG artifact detection and correction.

### Version 1.0.4 (January 08, 2020)
+ bugfix: `controller.edit_peaks()` works properly again.
+ enhancement: moved the modality menu to processing options.

### Version 1.0.3 (December 26, 2019)
+ enhancement: improved sensitivity of `ecg.ecg_peaks()` without decreasing
precision in moderately dynamic conditions (handbike) while maintaining
high performance in resting conditions (sitting). Performance has been
evaluated on lead 2 of all 25 subjects in the [Glasgow University Database (GUDB)](http://researchdata.gla.ac.uk/716/).
The GUBD has not been used to optimize `ecg.ecg_peaks()` in any way prior to
the performance evaluation. The tolerance for peak detection was set to one
sample.

|condition|metric     |summary|version 1.0.2|version 1.0.3
|:-------:|:---------:|:-----:|:-----------:|:-----------:
|sitting  |precision  |mean   |.999         |.998
|         |           |std    |.002         |.005
|         |sensitivity|mean   |.996         |.996
|         |           |std    |.008         |.004
|handbike |precision  |mean   |.904         |.930
|         |           |std    |.135         |.127
|         |sensitivity|mean   |.789         |.857
|         |           |std    |.281         |.247

### Version 1.0.2 (December 1, 2019)
+ enhancement: `resp.resp_extrema()` is now based on zerocrossings and makes
fewer assumptions about breathing rate.

### Version 1.0.1 (November 28, 2019)
+ bugfix: `controller.save_signal()` now preserves the header if the data are
saved to the same location (i.e., if `model.rpathsignal` and `model.wpathsignal` are
identical).
# User Guide

## Layout of the interface

![blank](images/screenshot_blank.png)

In the **menubar**, you can find the sections **_biosignal_**, **_peaks_**,
and **_statistics_**. These contain methods for the interaction with your biosignals. On the left
side, there's a panel containing **configurations** that allows you to customize your workflow.
You can toggle the visibility of the **configurations** with **_show/hide configurations_** in the **menubar**.
To the right of the **configurations** is the **datadisplay** which consists of three panels. The upper
panel contains the biosignal as well as peaks identified in the biosignal, while the middle panel can be used
to optionally display a marker channel. The lower panel contains any statistics
derived from the peaks. Adjust the height of the panels by dragging the segmenters between them up or down.
Beneath the **datadisplay**, in the lower left corner, you find the **displaytools**. These allow you to interact with the
biosignal. Have a look in the [functionality section](#functionality) for details on these elements.


## Getting started
The following work-flows are meant as an introduction to the interface. Many other
work-flows are possible. Note that `biopeaks` works with the OpenSignals file format
as well as EDF files, and plain text files (.txt, .csv, .tsv) containing biosignal channels as columns. The functions
used in the exemplary work-flow are described in detail in the [functionality section](#functionality).

### extracting ECG features from an OpenSignals file

1. Download the [example data](https://github.com/JanCBrammer/biopeaks/blob/master/docs/example_data.txt) (click the `Download` button on the upper right).

2. Configure the processing options in the **configurations**:
    * Since we want to analyze ECG data, make sure that **_processing options_** -> _modality_ is set to "ECG".
    * Set **_processing options_** -> _mode_  to "single file".
    * Set **_channels_** -> _biosignal_ to "A3" (ECG has been recorded on analog channel 3) and _marker_ to "I1" (events have been marked in input channel 1).
3. Load the example data using **menubar** -> **_biosignal_** -> _load_ -> _Opensignals_. More details on loading data can be found [here](#load-biosignal).

4. Once the data has been loaded, lets select a segment. **menubar** -> **_biosignal_** -> _select segment_ opens the **segmentdialog** on the right side of the **datadisplay**. Select the segment from the second marker at 150 seconds to the end of the signal at 340 seconds by entering "150" and "340" in the `start` and `end` fields respectively. Click `preview segment` to verify that the correct segment has been selected. Now, click `confirm segment` to cut out the selection. Note that the original file is not affected by this. More details on segmenting can be found [here](#segment-biosignal).

5. Now you can identify the R-peaks using **menubar** -> **_peaks_** -> _find_. Around second 10, you'll notice irregularities in the placement of the R-peaks. Zoom in on the biosignal by clicking on the magnifying glass in the **displaytools**. You can now use the mouse cursor to draw a rectangle around the biosignal between seconds 5 and 20. This will result in the magnification of that segment. Note that zooming is easier when you enlarge the panel containing the biosignal by dragging down the gray segmenter underneath the panel. You'll see a slightly misplaced R-peak around second 11 as well as a missing R-peak at second 13. In the next two steps, we'll correct these R-peaks using manual and automatic correction. The functionality section contains more details on how to [find extrema](#find-peaks) and [use the displaytools](#displaytools).

6. First, lets use autocorrection: **menubar** -> **_peaks_** -> _autocorrect_. This corrects pronounced irregularities in the R-peaks, such as missing R-peaks or relatively large misplacements. You'll see that the autocorrection adds the missing R-peak at 13 seconds. However, it does not correct the subtle misplacement of the R-peak at 11 seconds. More details on the autocorrection can be found [here](#auto-correct-peaks).

7. To correct the misplacement at 11 seconds you can use manual peak editing. Check the box next to **configurations** -> **peak** -> _editable_. Now click on the biosignal panel once to enable peak editing. Delete the misplaced R-peak by placing the mouse cursor in it's vicinity and pressing "d" on your keyboard. Next, you can insert the R-peak at the correct position by placing the mouse cursor on the R-peak at 11 seconds and pressing "a". You can play around with adding and deleting peaks to get a feeling for how it works. More details on editing peaks can be found [here](#edit-peaks).

8. Now you're ready to extract heart period and heart rate by clicking **menubar** -> **_statistics_** -> _calculate_. Both statistics will be displayed in the lower panel of the **datadisplay**. More details on calculating statistics can be found [here](#calculate-statistics).

9.  Finally, to be able to reproduce your results, save the segment (**menubar** -> **_biosignal_** -> _save_), peaks (**menubar** -> **_peaks_** -> _save_), and statistics (**menubar** -> **_statistics_** -> _save_). Prior to saving the statistics, check the boxes next to the statistics that you'd like to save: **configurations** -> **_select statictics for saving_**. You can now use the statistics and/or peaks for further analysis such as averages or heart rate variability.


### extracting respiration features from a custom file

You can analyze any plain text file (.txt, .csv, .tsv) as long as it contains biosignal(s) as column(s).

1. Download the [example data](https://github.com/JanCBrammer/biopeaks/blob/master/docs/example_data.txt) (click the `Download` button on the upper right).

2. Configure the processing options in the **configurations**:
    * Since we want to analyze respiration data, make sure that **_processing options_** -> _modality_ is set to "RESP".
    * Set **_processing options_** -> _mode_  to "single file".

3. Load the example data using **menubar** -> **_biosignal_** -> _load_ -> _Custom_. This will pop up a dialog that prompts you for information about the custom file. First, you need to specify which column contains the biosignal. Since the respiration channel is recorded in column 6 of the example data, fill in "6" in the `biosignal column` field. Similarly, the marker we'd like to display is recorded in column 2. Since the columns of the example data file are separated by tabs, you can select "tab" from the `column separator` dropdown menu. If the data file contains a header, specify the `number of header rows`. The example data has 3 header rows. Finally, you need to specify the `sampling rate`. The example data has been recorded at 1000 Hz. Click `continue loading file` to select the example data file from your file system. More details on how to load a custom file can be found [here](#plain-text-files).

4. Once the biosignal and marker are displayed in the upper and lower **datadisplay** panels respectively, you can select a segment. **menubar** -> **_biosignal_** -> _select segment_ opens the **segmentdialog** on the right side of the **datadisplay**. Select the segment from the start of the recording (0 seconds) to the first marker at 95 seconds by entering "0" and "95" in the `start` and `end` fields respectively. Click `preview segment` to verify that the correct segment has been selected. Now, click `confirm segment` to cut out the selection. Note that the original file is not affected by this. More details on segmenting can be found [here](#segment-biosignal).

5. Now you can use **menubar** -> **_peaks_** -> _find_ to identify the breathing extrema.

6. Once the extrema are displayed on top of the biosignal, you're ready to extract breathing period, -rate, and -amplitude by clicking **menubar** -> **_statistics_** -> _calculate_. All three statistics will be displayed in the lower panel of the **datadisplay**. More details on calculating statistics can be found [here](#calculate-statistics).

7. Finally, to be able to reproduce your results, save the segment (**menubar** -> **_biosignal_** -> _save_), breathing extrema (**menubar** -> **_peaks_** -> _save_), and statistics (**menubar** -> **_statistics_** -> _save_). Prior to saving the statistics, check the boxes next to the statistics that you'd like to save: **configurations** -> **_select statictics for saving_**. You can now use the statistics and/or peaks for further analysis.

## Functionality

- [load biosignal](#load-biosignal)
  - [OpenSignals and EDF](#opensignals-and-edf)
  - [Plain text files](#plain-text-files)
- [segment biosignal](#segment-biosignal)
- [save biosignal](#save-biosignal)
- [find peaks](#find-peaks)
- [save peaks](#save-peaks)
- [load peaks](#load-peaks)
- [calculate statistics](#calculate-statistics)
- [save statistics](#save-statistics)
- [edit peaks](#edit-peaks)
- [auto-correct peaks](#auto-correct-peaks)
- [batch processing](#batch-processing)
- [displaytools](#displaytools)


### load biosignal
#### OpenSignals and EDF
Under **configurations** -> **_channels_**, specify which channel contains the
_biosignal_ corresponding to your modality. Optionally, in addition to the
_biosignal_, you can select a _marker_. This is useful if you recorded a channel
that marks interesting events such as the onset of an experimental condition,
button presses etc.. You can use the _marker_ to display any other
channel alongside your _biosignal_. Once these options are selected,
you can load the biosignal: **menubar** -> **_biosignal_** -> _load_ -> _Opensignals_ or _EDF_. A
dialog will let you select a file.
#### Plain text files
If you have a .txt, .csv, or .tsv file that contains biosignal channels as columns,
you can load a biosignal channel and optionally a marker channel using **menubar** -> **_biosignal_** -> _load_ -> _Custom_.
This will open a dialog that prompts you for some information.

![customfile](images/screenshot_customfile.png)

Note that the fields _biosignal_column_, _number of header rows_, and _sampling rate_
are required, whereas _marker column_ is optional. In case your file doesn't have a header,
simply put "0" into the _number of header rows_ field. Also make sure to select the
correct _column separator_ ("comma" for .csv, "tab" for "tsv", "colon",
or "space" in case your columns are separated with these characters). Once you're done,
pressing _continue loading file_ opens another dialog that will let you select the file.


Once the biosignal has been loaded successfully it is displayed in the upper **datadisplay**. If
you selected a marker, it will be displayed in the middle **datadisplay**. The current file name is displayed in the lower right
corner of the interface. You can load a new biosignal from either the same file (i.e., another channel)
or a different file at any time. Note however that this will remove all data that is currently in the interface.

![biosignal](images/screenshot_biosignal.png)

### segment biosignal
**menubar** -> **_biosignal_** -> _select segment_ opens the **segmentdialog**
on the right side of the **datadisplay**.

![segmentdialog](images/screenshot_segmentdialog.png)

Specify the start and end of the segment in seconds either by entering values,
or with the mouse. For the latter option, first click on
the mouse icon in the respective field and then left-click anywhere on the
upper **datadisplay** to select a time point. To see which time point is
currently under the mouse cursor have a look at the x-coordinate
displayed in the lower right corner of the **datadisplay** (displayed when you hover
the mouse over the upper **datadisplay**). If you click **_preview segment_**
the segment will be displayed as a shaded region in the upper **datadisplay**
but the segment won't be cut out yet.

![segmenthighlight](images/screenshot_segmenthighlight.png)

You can change the start and end values and preview the segment until you are certain that the desired segment is
selected. Then you can cut out the segment with **_confirm segment_**. This also closes the **segmentdialog**. Alternatively, the
**segmentdialog** can be closed any time by clicking **_abort segmentation_**.
Clicking **_abort segmentation_** discards any values that might have been
selected. You can segment the biosignal multiple times. Other data (peaks,
statistics) will be also be segmented if they are already computed. Note that
the selected segment must have a minimum duration of five seconds. Also, after
the segmentation, the signal starts at second 0 again. That is, relative timing
is not preserved during segmentation. The original file is not affected by the segmentation.

### save biosignal
**menubar** -> **_biosignal_** -> _save_ opens a dialog that lets you
select a directory and file name for saving the biosignal.
Note that saving the biosignal is only possible after segmentation. The file is
saved in its original format containing all channels.

### find peaks
Before identifying peaks, you need to select the modality of your biosignal
in **configurations** -> **_processing options_** -> _modality_ ("ECG" for
electrocardiogram, "PPG" for photoplethysmogram, and "RESP" for breathing).
This is important since `biopeaks` uses modality-specific peak detectors.
Then, **menubar** -> **_peaks_** -> _find_ automatically identifies the peaks in the
biosignal. The peaks appear as dots displayed on top of the biosignal.

![peaks](images/screenshot_peaks.png)

### save peaks
**menubar** -> **_peaks_** -> _save_ opens a file dialog that lets you select a
directory and file name for saving the peaks to a CSV file. The format of the
file depends on the _modality_. For ECG and PPG, `biopeaks` saves a column containing
the occurrences of R-peaks or systolic peaks respectively in seconds. The first element contains the header
"peaks". For breathing, `biopeaks` saves two columns containing the occurrences
of inhalation peaks and exhalation troughs respectively in seconds. The first
row contains the header "peaks, troughs". Note that if there are less peaks
than troughs or vice versa, the column with less elements will be padded with
a NaN.

### load peaks
**menubar** -> **_peaks_** -> _load_ opens a file dialog that lets you select
the file containing the peaks. Note that prior to loading the peaks, you have
to load the associated biosignal. Also, loading peaks won't work if there are
already peaks in memory (i.e., if there are already peaks displayed in the
upper **datadisplay**). Note that it's only possible to load peaks that have
been saved with `biopeaks` or adhere to the same format. The peaks appear as
dots displayed on top of the biosignal.

### calculate statistics
**menubar** -> **_statistics_** -> _calculate_ automatically calculates all
possible statistics for the selected _modality_. The statistics will be
displayed in the lowest **datadisplay**.

![statistics](images/screenshot_statistics.png)

### save statistics
First select the statistics that you'd like to save: **configurations** ->
**_select statictics for saving_**. Then,
**menubar** -> **_statistics_** -> _save_, opens a file dialog
that lets you choose a directory and file name for saving a CSV file. The
format of the file depends on the _modality_. Irrespective
of the modality the first two columns contain period and rate (if both have
been chosen for saving).
For breathing, there will be an additional third column containing the tidal
amplitude (if it has been chosen for saving). The first row contains the
header. Note that the statistics are linearly interpolated to match the biosignal's
timescale (i.e., they represent instantaneous statistics sampled at the biosignal's sampling rate).

### edit peaks
It happens that the automatic peak detection places peaks wrongly or fails to
detect some peaks. You can
catch these errors by visually inspecting the peak placement. If you spot
errors in peak placement you can correct those manually. To do so make sure to
select **configurations** -> **peak** -> _editable_. Now click on the
upper **datadisplay** once to enable peak editing. To delete a peak place the
mouse cursor in it's vicinity and press "d". To add a peak,
press "a". Editing peaks is most convenient if you zoom in on the biosignal
region that you want to edit using the [**displaytools**](#displaytools).
The statistics in the lowest **datadisplay**
can be a useful guide when editing peaks. Isolated, unusually large or small
values in period or rate can indicate misplaced peaks. Note, that when editing breathing
extrema, any edits that break the alternation of peaks and troughs
(e.g., two consecutive peaks) will automatically be discarded when you save
the extrema. If you already calculated statistics, don't forget to calculate
them again after peak editing.

### auto-correct peaks
If the _modality_ is ECG or PPG, you can automatically correct the peaks with
**menubar** -> **_peaks_** -> _autocorrect_. Note that the auto-correction tries
to spread the peaks evenly across the signal which can lead to peaks that are
slightly misplaced. Also, the auto-correction does not guarantee that
all errors in peak placement will be caught. It is always good to check for errors manually!

### batch processing
> There is no substitute for manually checking the biosignal's
> quality as well as the placement of the peaks. Manually checking and editing
> peak placement is the only way to guarantee sensible statistics. Only use
> batch processing if you are sure that the biosignal's quality is sufficient!

You can configure the batch processing in the **configurations**.
To enable batch processing, select
**_processing options_** -> _mode_ -> "multiple files". Make sure to
select the correct _modality_ in the **_processing options_** as well. Also select
the desired _biosignal channel_ in **_channels_** (for custom files you can
specify the _biosignal_column_ in the file dialog while loading the files). Further, indicate if you'd
like to save the peaks during batch processing: **_peak options_** ->
_save during batch processing_. You can also choose to apply the auto-correction
to the peaks by selecting **_peak options_** ->
_correct during batch processing_. Also, select the statistics you'd like
to save: **_select statictics for saving_**. Now, select
all files that should be included in the batch: **menubar** -> **_biosignal_**
-> _load_. A dialog will let you select the files (select multiple files with
the appropriate keyboard commands of your operating system). Next, a dialog
will ask you to choose a directory for saving the peaks (if you enabled that
option). The peaks will be saved to a file with the same name as the biosignal
file, with a "_peaks.csv" extension.
Finally, a dialog will ask you to select a directory for saving the statistics
(if you chose any statistics for saving). The statistics will be saved to a
file with the same name as the biosignal file, with a "_stats.csv" extension. Once all
dialogs are closed, `biopeaks` carries out the following actions for each file
in the batch: loading the biosignal, identifying the
peaks, calculating the statistics and finally saving the desired data (peaks
and/or statistics). Note that nothing will be shown in the **datadisplay**
while the batch is processed. You can keep track of the progress by looking
at the file name displayed in the lower right corner of the interface.
Note that segmentation or peak editing are not possible during batch
processing.

### displaytools
The **displaytools** allow you to interact with the biosignal. Have a look
[here](https://matplotlib.org/3.1.1/users/navigation_toolbar.html) for a
detailed description of how to use them.
Welcome to `biopeaks`, a straightforward graphical user interface for feature extraction from electrocardiogram (ECG), photoplethysmogram (PPG) and breathing biosignals.
It processes these biosignals semi-automatically with sensible defaults and offers the following functionality:


+ processes files in the open biosignal formats [EDF](https://en.wikipedia.org/wiki/European_Data_Format), [OpenSignals (Bitalino)](https://bitalino.com/en/software)
as well as plain text files (.txt, .csv, .tsv)
+ interactive biosignal visualization
+ biosignal segmentation
+ benchmarked, automatic extrema detection (R-peaks in ECG, systolic peaks in PPG, exhalation troughs and inhalation
peaks in breathing signals) with signal-specific, sensible defaults
+ automatic state-of-the-art [artifact correction](https://www.tandfonline.com/doi/full/10.1080/03091902.2019.1640306)
 for ECG and PPG extrema
+ manual editing of extrema
+ extraction of instantaneous features: (heart- or breathing-) rate and period, as well as breathing amplitude
+ .csv export of extrema and instantaneous features for further analysis (e.g., heart rate variability)
+ automatic analysis of multiple files (batch processing)


Have a look at the documentation to find out more:

+ [Installation](installation.md)
+ [User Guide](user_guide.md)
+ [Contributor Guide](contributor_guide.md)
+ [Changelog](changelog.md)
+ [Further Resources](additional_resources.md)
+ [Citation](citation.md)
+ [GitHub repository](https://github.com/JanCBrammer/biopeaks)
# Installation

## As executable

For Windows you can download [biopeaks.exe](https://github.com/JanCBrammer/biopeaks/releases/latest). Double-click the 
executable to run it. You don't need a Python installation on your computer to run the executable.
Currently, there are no executables available for macOS or Linux
(please [open an issue](https://help.github.com/en/github/managing-your-work-on-github/creating-an-issue) if you're interested).

## As Python package

### Instructions for users without a Python installation
If you don't have experience with installing Python packages and/or if you
aren't sure if you have Python on your computer start by setting up Python.
Go to <https://docs.conda.io/en/latest/miniconda.html> and install the latest
miniconda distribution for your operating system.
Follow these [instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
in case you're unsure about the installation. Once you've installed miniconda, open the
[Anaconda Prompt](https://docs.anaconda.com/anaconda/user-guide/getting-started/)
and run the following commands (hit enter once you've typed each of the lines below and wait for
the commands to be executed):

```
conda create -y -n biopeaks python=3.9
conda activate biopeaks
pip install scipy numpy matplotlib pandas PySide6 biopeaks 
```

After the successful installation, open the application by typing
```
biopeaks
```
Note that every time you open the Anaconda Prompt, you need to activate the
biopeaks environment before starting the application:
```
conda activate biopeaks
biopeaks
```

### Instructions for users who already have a Python installation
Have a look at the project's [pyproject.toml file](https://github.com/JanCBrammer/biopeaks/blob/master/pyproject.toml)
for an up-to-date list of the dependencies. In order to manage the dependencies, it is highly recommended to install
`biopeaks` into an isolated environment using [miniconda](https://docs.conda.io/en/latest/miniconda.html),
[Poetry](https://python-poetry.org/), or other tools for creating and managing virtual environments.

Once you've set up an environment containing all the dependencies, install `biopeaks` with

```
pip install biopeaks
```

You can then open the application by typing

```
biopeaks
```
# Further Resources
Check out these free alternatives:
+ <https://www.medisig.com/signalplant/>
+ <https://physiodatatoolbox.leidenuniv.nl/>
+ <http://www.artiifact.de/># Contributor Guide

Thanks for your interest in contributing to `biopeaks`! Please have a look at the [code of conduct](https://github.com/JanCBrammer/biopeaks/blob/master/code_of_conduct.md).


## Ways to contribute

### Reporting bugs or asking questions

Please report bugs or ask questions by [opening an issue](https://help.github.com/en/github/managing-your-work-on-github/creating-an-issue)
in the [`biopeaks` repository](https://github.com/JanCBrammer/biopeaks).


### Improve documentation, tests, or code

If you plan to contribute relatively large changes, please [open an issue](https://help.github.com/en/github/managing-your-work-on-github/creating-an-issue)
in the [`biopeaks` repository](https://github.com/JanCBrammer/biopeaks) before
you start working on your contribution. This way we can discuss your plans before you start writing/coding.

You can follow these steps to contribute documentation, tests, or code:

1. [Fork](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo) the [`biopeaks` repository](https://github.com/JanCBrammer/biopeaks).
2. Add a `topic` branch with a descriptive name to your fork. For example, if you want to contribute an improvement to the documentation you could call the `topic` branch `improve_docs`.
3. Install `biopeaks` in development mode:
   1. Make a [local clone of your fork](https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github/fork-a-repo#step-2-create-a-local-clone-of-your-fork).
   2. Navigate to the directory containing the cloned fork.
   3. Install `biopeaks` with `pip install -e .` The `e` stand for editable, meaning all the changes you make to the cloned fork take immediate effect. The `.` simply tells pip to install the content of the current directory.
4. Implement your contribution in the `topic` branch, following the [conventions](#conventions).
5. [Make a pull request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork) from the `topic` branch on your fork to the [`dev` branch of the `biopeaks` repository](https://github.com/JanCBrammer/biopeaks/tree/dev).
6. Once all CI tests pass and your changes have been reviewed, your PR will be merged and you're a contributor!


## Conventions

### General

* avoid introducing new dependencies
* write [numpydoc](https://numpydoc.readthedocs.io/en/latest/format.html) docstrings
  for every (non-private) new function
* add tests if you contribute code that is not covered by the existing [tests](#tests)

### Code style

* aim for simplicity and readability
* follow [PEP8 guidelines](https://www.python.org/dev/peps/pep-0008/)
* write code in Python 3.6+

### Architecture

The GUI is structured according to a variant of the
[model-view-controller architecture](https://martinfowler.com/eaaDev/uiArchs.html).
To understand the relationship of the `model`, `view`, and `controller` have a look
at how each of them is instantiated in [`__main__.py`](https://github.com/JanCBrammer/biopeaks/blob/master/biopeaks/__main__.py).
For example, the `view` has references to the `model` as well as the
`controller`, whereas the `model` has no reference to any of the other
components of the architecture (i.e., the `model` is agnostic to the `view` and
`controller`).


## Documentation

The documentation is hosted on GitHub pages, a static website associated
with the `biopeaks` repository: <https://jancbrammer.github.io/biopeaks/>. It
is automatically build from the `/docs` folder in the root of the `biopeaks` repository.
The website is re-build every very time the content of `/docs` changes on the
master branch (pushes, merged pull requests). `/docs` includes an `index.md` file
that constitutes the "landing page". It contains links to all other parts of the
documentation. The layout of the website is defined in `/docs/layouts`. For
additional information, head over to the [GitHub pages documentation](https://docs.github.com/en/free-pro-team@latest/github/working-with-github-pages).


## Tests

The OpenSignals test data have been recorded with<br/>
software: opensignals v2.0.0, 20190805<br/>
hardware: BITalino (r)evolution (firmware 1281)

The EDF test data have been downloaded from https://www.teuniz.net/edf_bdf_testfiles/

All test data are part of the biopeaks installation and do not have to be downloaded.

Please make sure to have [pytest](https://docs.pytest.org/en/latest/) as well as
[pytest-qt](https://pypi.org/project/pytest-qt/) installed before running the tests.

The tests can then be run in the test directory with [pytest](https://docs.pytest.org/en/latest/):
```
pytest -v
```


## Algorithm benchmarks

### ECG
To validate the performance of the ECG peak detector `heart.ecg_peaks()`, please install the [wfdb](https://github.com/MIT-LCP/wfdb-python) and [aiohttp](https://github.com/aio-libs/aiohttp).

You can then run the `benchmark_ECG_stream` script in the `benchmarks` folder. The script streams ECG and annotation files from the [Glasgow University Database (GUDB)](http://researchdata.gla.ac.uk/716/).
You can select an experiment, ECG channel, and annotation file.

Alternatively, you can download the GUDB and run the `benchmark_ECG_local` script in the `benchmarks` folder. In the script, replace the `data_dir` with your local directory (see comments in the script).

### PPG

To validate the performance of the PPG peak detector `heart.ppg_peaks()`
please download the [Capnobase IEEE TBME benchmark dataset](http://www.capnobase.org/index.php?id=857) and install [wfdb](https://github.com/MIT-LCP/wfdb-python) and [h5py](https://www.h5py.org/).

You can then run the `benchmark_PPG_local` script in the `benchmarks` folder. In the script, replace the `data_dir` with your local directory (see comments in the script).


## Resources

### [Using git](https://github.com/dictcp/awesome-git)

### [Using GitHub](https://docs.github.com/en)


## Local development on Windows

### Editable development
To develop and build `biopeaks` locally on Windows I found the following to be an ok solution (albeit somewhat hacky):
Set up a minimal Python environment using [miniconda](https://docs.conda.io/en/latest/miniconda.html). This environment merely contains Python and pip.
```
conda create --name biopeaks_dev python=3.9
```
Within the conda environment, use [Poetry](https://python-poetry.org/) to manage dependencies.
```
conda activate biopeaks_dev
pip install poetry
```
Note that Poetry usually creates its own virtual environments. However, we'll use it
inside an existing environment. To make sure Poetry doesn't redundantly nest a virual environment within the conda environment we run
```
poetry config virtualenvs.create false --local
```
Now we can install an editable version of `biopeaks` alongside its depencencies with
```
poetry install --extras pyinstaller
```

### Building executable with PyInstaller
Create an additional environment that contains only the build dependencies in
order to reduce the build size. Configure the build environment just like the [editable development environment](#editable-development).
```
conda create --name biopeaks_build python=3.9
conda activate biopeaks_build
pip install poetry
poetry config virtualenvs.create false --local
```
Now we use Poetry to only install the build dependencies, leaving out the development dependencies.
The latter would unnecessarily increase the build size.
```
poetry install --no-root --no-dev --extras "pyinstaller"
```
Now we can build the application from the root of the repository using PyInstaller.
Note that PyInstaller needs access to the `__main__.py` entry-point as if the file
would be located outside the `biopeaks` sub-directory (since `biopeaks` is imported
using absolute imports inside `__main__.py`). This is why we need to pass the root (`.`)
to the PyInstaller paths. For more details see https://pyinstaller.readthedocs.io/en/stable/runtime-information.html.
```
pyinstaller --onefile --windowed --name=biopeaks --paths=. \
--icon=biopeaks\images\python_icon.ico biopeaks\__main__.py 
```
# Citation

Please refer to the [biopeaks paper](https://joss.theoj.org/papers/10.21105/joss.02621) in The Journal of Open Source Software.
