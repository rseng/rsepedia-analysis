# *k*-Nearest Neighbors Classifier-Based Sleep Staging (*k*NNSS)

## Welcome to *k*NNSS software

*k*NNSS is a Matlab package for automated sleep stage scoring using the *k*-nearest neighbors algorithm. Compared to other automated sleep scoring packages its main advantage is simplicity and the use of physiologically relevant, human-interpretable features. The code is documented in the following JOSS publication: [![DOI](https://joss.theoj.org/papers/10.21105/joss.02377/status.svg)](https://doi.org/10.21105/joss.02377)

## LICENSE

    /* This Source Code Form is subject to the terms of the Mozilla Public
     * License, v. 2.0. If a copy of the MPL was not distributed with this
     * file, You can obtain one at http://mozilla.org/MPL/2.0/.

## Using *k*NNSS
### Compatibility and dependencies
*k*NNSS was tested in Matlab R2014b and R2018a. The software uses the Signal Processing Toolbox and the Statistics and Machine Learning Toolbox from Matlab. Also, some functions from the [Chronux Toolbox](http://chronux.org)  are used to calculate spectral properties of data (see below for a list of files). European Data Format (EDF; https://www.edfplus.info/) files are read using Brett Shoelson's (brett.shoelson@mathworks.com) edfread() function. Hjorth parameters are calculated by hjorth() a function written by Alois Schloegl (a.schloegl@ieee.org) as part of the [BIOSIG-toolbox](http://biosig.sf.net). This GitHub project is self-containing, no further downloads are required.

#### Chronux Toolbox file list
change_row_to_column.m, dpsschk.m, getfgrid.m, getparams.m, mtfftc.m, mtspecgramc.m, mtspectrumc.m

### Folder content

- Example_Data: this directory contains three sub-directories. IntermRes is initially empty, intermediate files generated during the execution of the test scripts will be stored here. RawData contains polysomnographycal recordings from Thomas Kilduff's laboratory at SRI International in EDF format. Sub-folder EDF contain raw recordings, sub-folder FFT stores associated manual sleep scores used to train classifiers and test prediction efficacy.

- Function_Library: this directory stores all the necessary Matlab functions required by *k*NNSS. Functions inputs and outputs are described in detail in the function header. Use Matlab's help to read about details of each function. A list of high-level purpose of functions is found in the Readme in this folder.

- Software_Verification: a set of scripts to test the code on your system are in this folder. For a detailed description see the Readme in this folder.

- JOSS_Paper: contains the manuscript submitted to the Journal of Open Source Software. The journal paper will provide background to the software, a summary of our work, applications and references.

### Installation
All functions are stored in the Function_Library folder. Simply adding this folder to Matlab's function path will take care of installation. Alternatively, check out contents of the Software_Verification folder.

## Top-level inputs and outputs
*k*NNSS reads raw electrophysiological data from file. Currently the European Data Format and Matlab files are supported but see generate_statespace.m in Function_Library to see how additional formats can be added. Manual scores are read together with the data. These are typically text files that come in many format.

Once the algorithm finishes predicting labels the output is a Matlab structure type variable, containing a cell array of strings in its fields. Strings correspond to predicted labels.

## Getting help
For more information please do not hesitate to contact Tamás Kiss (kiss.t (at) wigner.hu).

Thank you!
# Contributing

Thank you for helping to make this package better!

If you use our code and run into a bug or any problem, please open an issue here. Any questions, comments or suggestions can be sent via GitHub or send an e-mail to Tamás Kiss at kiss.t (at) wigner.hu.

Contributions to the code are very welcome. Since this repository follows the standard open source protocol, please use the basic [GitHub workflow](https://guides.github.com/introduction/flow/) to suggest a modification or addition. TL;DR: fork this repository and clone it to your own machine. Then make your changes, push your work back up to the fork on GitHub and open a [pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests) so that we can review your changes and discuss them with you before accepting and merging them to the main project.

Please restrict your commit to one modification per commit and ideally form a legible commit history. In your pull request, please describe clearly and concisely your motivation and the feature you developed.

Also, please try to keep your code easily readable and tidy using [smart indenting](https://www.mathworks.com/help/matlab/matlab_prog/improve-code-readability.html). We feel that writing help to functions greatly improves understanding of the code and its function so please follow the MATLAB basic structure for `help` text: [Add Help for Your Program](https://www.mathworks.com/help/matlab/matlab_prog/add-help-for-your-program.html).

## Code of conduct

Please be polite :-) For a more detailed version, please see [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/version/2/0/code_of_conduct).
# Repository for Example Dataset

This folder has two subfolders:

- RawData: This folder stores the abridged version of two datasets allowing the user to assess software capabilities: One from male Trace Amine-Associated Receptor 1 (TAAR1) knockout mice (file names starting with 'tko') described in detail in [Schwartz *et al*. (2018)](https://doi.org/10.3389/fphar.2018.00035), and the other from male Sprague-Dawley rats (file names starting with 'A') collected in the Sleep Neurobiology Laboratory at [SRI International](https://www.sri.com/). The EDF subfolder contains electrophysiological recordings in the [European Data Format](https://www.edfplus.info), while the associated FFT folder presents the corresponding state designations assigned to each 10-sec recording epoch.

- IntermRes, which is initially empty, is the folder that will store results of intermediate calculations, like the feature table file (features & manual scores), the partitioned training and test sets, as well as the trained classifiers.
This folder is needed to store intermediate results temporarily to
avoid multiple re-calculations.
# High-level function description for *k*NNSS

This readme describes what each function of the *k*NNSS package is used for. A detailed description of mandatory and optional inputs for functions, as well as their outputs can be found in the function headers.

## Top-level functions

generate_statespace.m
:  This is the first function to use in the processing chain. It loads raw data, cuts it into epochs and calculates feature describing each epoch. It also loads manual scores associated with the raw data and labels each epoch using manual scores (optional). Feature data (since its calculation can be time consuming) can be saved in Matlab files. This function is written with modularity in mind: if further features are to be included their calculation can be implemented in a separate function and values passed to generate_statespace that merges them and saves a unified feature file. An example use is shown in the FirstStep.m script in SoftwareVerification folder.

train_one_model.m
: This function uses features calculated by generate_statespace.m and trains a *k*-Nearest Neighbors (*k*NN) classifier. This function gathers feature information from multiple experiments (could be from the same or from different subjects) and uses the merged data to train a single classifier for all experiments. This way a general classifier is created that can be used across multiple experiments. Its use is presented in the SecondStep.m script in SoftwareVerification folder.

train_many_models.m
: Similarly to train_one_model.m, this function uses pre-calculated feature values to train a *k*NN classifier. It differs from train_one_model.m in that it trains a classifier for each experiment it gets. Individual classifiers are more accurate for any given subject, however, they do not generalize as well as models trained on multiple subjects or conditions. Its use is presented in the SecondStep.m script in SoftwareVerification folder.

evaluate_model_goodness.m
: This function calculates confusion matrices from manually and automatically determined label sets and displays a set of measures to describe how well the classifiers perform.

## Helper functions
### Functions from the [Chronux toolbox](http://chronux.org)  used for spectrum  calculations
change_row_to_column.m
: Helper routine to transform 1d arrays into column vectors that are needed by other routines in Chronux

dpsschk.m
: Helper function to calculate tapers

getfgrid.m
:  Helper function that gets the frequency grid associated with a given FFT based computation

getparams.m
: Helper function to convert structure params to variables used by the various routines in Chronux

mtfftc.m
: Multi-taper fourier transform - continuous data

mtspecgramc.m
: Multi-taper time-frequency spectrum - continuous process

mtspectrumc.m
: Multi-taper spectrum - continuous process

### Other helper functions

edfread.m
: Read European Data Format file into MATLAB. This function was written by Brett Shoelson, PhD (brett.shoelson@mathworks.com).

canonize_fieldname.m
: This function attempts to create properly formatted field names.

confusion2PerformanceMetrics.m
: This function converts the confusion matrix into different preformance metrics including SEN, FPR, SPC, ACC, PPV, NPR. For details see the [Wikipedia page for confusion matrix](https://en.wikipedia.org/wiki/Confusion_matrix).

estimate_feature_goodness.m
: This function uses [filter methods](https://en.wikipedia.org/wiki/Feature_selection#Filter_method) to select features used for classifier training that are most dissimilar.

hjorth.m
: Calculates the [Hjorth parameters](https://en.wikipedia.org/wiki/Hjorth_parameters) activity, mobility, and complexity. This function is part of the [BIOSIG-toolbox](http://biosig.sf.net/) and was written by Alois Schloegl (a.schloegl@ieee.org).

select_features.m
: The pair of estimate_feature_goodness.m, this function uses the [wrapper method](https://en.wikipedia.org/wiki/Feature_selection#Wrapper_method) (sequential feature selection here) to find best features to be used for training the classifier.

flag_outliers.m
: This function finds outliers in a feature table using median and mad.

metric_band_power.m
:  This function calculates power in given frequency bands using multitaper power estimation.

preprocess_features.m
: Feature pre-processing is performed by this function. It is called before fitting a classifier or before predicting new data. Pre-processing steps include removal of data with artifacts, missing data segments, removal of manually selected epochs, normalization, and log transformation.

train_classifier.m
: This function trains the KNN classifier.# Code to Test *k*NNSS on Your System

## Outline
The short example stored in this directory guides you through typical steps of sleep stage scoring. There are two independent datasets you can work on, one was recorded in transgenic mice, the other in control rats (see Readme in the Example_Data folder for further information on the data).

- To begin, please edit the `par.ROOTDIR` variable in setup.m file to specify the folder in which *k*NNSS is placed. setup.m will set all other path to data and functions.

- Using our supervised classification method the first step to take is the calculation of features that will describe each epoch (time window) of the recording. The FirstStep.m script will do this. Also, features are combined with manual scores to form a training set. Manual scores in this example are simplified and some labels are merged. Once features are calculated they are saved in files since their calculation might be quite time consuming in case of many and long recordings. The calculated features (plus manual scores) are partitioned into a training and a test set. The former will be used to label feature space of a *k*NN classifier, while the latter will serve to test reliability of prediction.

- Model fitting takes place in SecondStep.m. This script fits (and stores in Matlab files) two types of models. On one hand, when using the  `train_one_model` function, data from all mice are combined into a single training set and used to train the classifier. On the other hand, `train_many_models` labels one separate classifier for each animal. These models use data in the training set created during the previous step.

- Prediction of the unseen part of the data (the test set) happens in ThirdStep.m. Both types of models are tested: the single model is used to consecutively predict data from all animals, and subsequently data from each animal is predicted by its dedicated classifier as well. Since rapid eye movement (REM) sleep epochs are very scarce in the last step the single model is "deflated" by removing a number or non-REM sleep epochs. Successful testing of the software will output three named figures similar to the ones below (since feature selection has a random part to it, results might slightly vary). For the transgenic mouse data the prediction accuracy turns out as shown in Figs 1-3.

![Output figures for mouse](ExampleResultFigMouse.png)

For the rat data please see Figs 4-6 below for performance.

![Output figures for rat](ExampleResultFigRat.png)
---
title: 'Automated Sleep Stage Scoring Using *k*-Nearest Neighbors Classifier'

tags:
  - Supervised clustering
  - Polysomnography
  - Power spectrum
  - MATLAB

authors:
  - name: Tamás Kiss
    orcid: 0000-0001-6360-0714
    affiliation: "1, 2"

  - name: Stephen Morairty
    orcid: 0000-0002-0781-1645
    affiliation: 3

  - name: Michael Schwartz
    orcid: 0000-0002-5464-638X
    affiliation: 3

  - name: Thomas S.\ Kilduff
    orcid: 0000-0002-6823-0094
    affiliation: 3

  - name: Derek L.\ Buhl
    orcid: 0000-0003-4433-7150
    affiliation: "1, 4"

  - name: Dmitri Volfson
    orcid: 0000-0002-5167-7834
    affiliation: "1, 4"

affiliations:
 - name: Global Research and Development, Pfizer Inc, Groton, CT, USA
   index: 1

 - name: Department of Computational Sciences, Wigner Research Centre for Physics, Budapest, Hungary
   index: 2

 - name: Center for Neuroscience, SRI International, Menlo Park, CA, USA
   index: 3

 - name: Current affiliation -- Takeda Pharmaceuticals, Inc., Cambridge, MA, USA
   index: 4

date: 1 May 2020
bibliography: paper.bib
---

# Polysomnographic Sleep Stage Scoring

Many features of sleep, such as the existence of rapid eye movement
(REM) sleep or non-REM sleep stages, as well as some of the underlying
physiological mechanisms controlling sleep, are conserved across
different mammalian species. Sleep research is important to
understanding the impact of disease on circadian biology and optimal
waking performance, and to advance treatments for sleep disorders,
such as narcolepsy, shift work disorder, non-24 sleep-wake disorder,
and neurodegenerative disease.  Given the evolutionary relatedness of
mammalian species, sleep architecture and changes therein may provide
reliable translational biomarkers for pharmacological engagement in
proof-of-mechanism clinical studies.

Key physiological indicators in sleep include electroencephalography
(EEG) or electrocorticography, electrooculography (EOG), and
electromyography (EMG). Polysomnography (PSG) is the simultaneous
collection of some or all of these measurements and is typically
performed in a specialized sleep laboratory. Determination of the wake
or sleep stage someone is in (i.e., wake, REM sleep, or non-REM sleep,
which is broken down into stages 1, 2, or 3), relies on the judgment
of a trained professional who scores the data based on the
standardized criteria for the recording and staging of human PSG set
forth by @AASMM. Disagreement between individual recordings might
arise due to differences in instrumentation or to the subjective
opinion of the individual scoring the stages. Animal sleep studies
show even greater variability [@ROBERT1999111], as each laboratory
uses methods that best suit their individual needs (e.g.,
electrode/reference positions, muscle choice for EMG implantation, use
of EOG, etc.). While these technical differences make it difficult to
compare studies, the variability in scoring of sleep stages makes it
even more challenging. Although numerous scoring algorithms exist,
most are unreliable, especially following drug treatment. After nearly
half a century of PSG studies, the gold standard of scoring sleep
architecture remains a complete and thorough examination of the PSG
signals, which are scored in 4-, 10-, or 12-second epochs in animal
studies and 30-second epochs in human studies, making it very
difficult to screen through drugs in animal studies and cumbersome to
implement large clinical trials.

# Applications and Advantage

To expedite the tedious process of visually analyzing PSG signals and
to further objectivity in the scoring procedure, a number of sleep
staging algorithms have been developed both for animals
[@STEPHENSON2009263;@BASTIANINI2014277;@barger2019;@vladimir2020] and human
subjects [@PENZEL2000131;@gunn2020;@zhang2020] as reviewed most
recently by @fiorillo2019 and @faust2019. However, computer-based
methods are typically tested on data obtained from healthy subjects or
control animals, and performance is assessed only in a few cases in
subjects with sleep disorders or following drug treatment
[@BOOSTANI201777;@allocca2019]. Furthermore, scoring sleep for
hundreds of animals in a typical preclinical drug discovery effort
often becomes a bottleneck and a potential source of
subjectivity affecting research outcomes.

In this paper, we present an automated approach intended to eliminate
these potential issues. The initial application of our approach is for
basic and discovery research in which experiments are conducted in
large cohorts of rodents, with the expectation that results can be
translated to higher-order mammals or even humans. Building on
features classically extracted from EEG and EMG data and machine
learning-based classification of PSG, this approach is capable of
staging sleep in multiple species under control and drug-treated
conditions, facilitating the detection of treatment-induced changes or
other manipulations (e.g., genetic). Using human interpretable
features calculated from EEG and EMG will be important to understand
drug mechanisms, for prediction of treatment outcomes, and as
biomarkers or even translational biomarkers. For example, one of the
features used by the algorithm is the power in the theta frequency
band (called `eeg_theta` in the code), which is the 4 Hz to 12 Hz
range and it is known that an increase of theta activity together with
low EMG activity (our relevant features are called `emg_high` and
`emg_RMS`) are the hallmark of REM sleep (see the
[figure](https://en.wikipedia.org/wiki/Rapid_eye_movement_sleep#/media/File:Normal_EEG_of_mouse.png)
in @wiki:REM). However, theta power is also associated with other
phenomena, like anxiety [@John2014], thus our `eeg_theta` feature,
besides being used for sleep scoring can also be used as a biomarker
of drug effect.

Multiple software applications have been developed to address the
problem of automated sleep stage scoring. In their comparative review,
@BOOSTANI201777 found that the best results could be achieved when
entropy of wavelet coefficients along with a random forest classifier
were chosen as feature and classifier, respectively. Another recent
method [@miladinovic2019] used cutting-edge machine learning methods
combining a convolutional neural network-based architecture to produce
domain invariant predictions integrated with a hidden Markov model to
constrain state dynamics based upon known sleep physiology. While our
method also builds on machine learning techniques, it is based on
interpretable features and uses a simpler algorithm for classification
-- which should make it an ideal choice for the broader community as
well as for sleep experts who might not be too familiar with complex
machine learning approaches. Furthermore, we chose not to constrain
the number of identifiable sleep/wake states or the probability of
transition from one state to another, as we and others have found that
drug interventions [@Harvey2013] and disease processes [@Mooij2020]
tend to change not only the amount of time spent in different sleep
stages but their transition probabilities as well. Finally, our method
is a supervised method that requires a training set. While this might
seem to be a disadvantage over non-supervised methods, we have found
that drug treatment or pathological conditions can result in sleep
stages not observed in healthy controls. Thus, the algorithm must be
trained to these new stages.

# Brief Software Description

Our software package, implemented in Matlab, is available for download
on GitHub [@kNNSS]. Automatic sleep staging consists of the classical
consecutive steps of machine learning-based sleep scoring algorithms
\autoref{fig:method}. First, offline stored EEG and EMG data are
loaded into memory to allow for the uniform processing of time-series
data and segmented into consecutive 10-second, non-overlapping epochs
that correspond to manually scored epochs. Second, features are
extracted from the raw signal for all epochs. Features consist of the
power contained in physiologically-relevant frequency bands, as well
as Hjorth parameters for both EEG and EMG data. Third, features
undergo a pre-processing step including the following operations:
unusable epochs that contain too much noise or contain no signal are
removed.  Features are then transformed using the logarithm function
making feature distributions more Gaussian-like, thereby
facilitating subsequent machine classification.  Finally, each feature
is normalized to its median wake value within an animal to enable
usability of the algorithm across laboratories.  Wake periods can be
identified before running the algorithm using the manually-scored
training set or an experiment can be performed such that a given
period is expected to be comprised of an extended period of
wakefulness.  Following feature extraction, a combined filter and
wrapper method-based feature selection step is applied.  This step
ensures that features with the most predictive value are chosen and
also helps to prevent over-fitting. For classification, the
*k*-nearest neighbors classifier is used on data pre-processed
following the procedure described above.

![Summary of training and using the *k*-nearest neighbors algorithm
 for predicting sleep stage labels.\label{fig:method}](fig1.png)

The algorithm was used to predict sleep stages in mice
(\autoref{fig:efficM}), rats (\autoref{fig:efficR}) and non-human
primates (data not shown). Prediction accuracy was found to depend on
a number of parameters of the input data, including consistency of
manual scores and physiological signals, as well as the amount of
artifacts. Furthermore, relative frequency of predicted labels can
influence efficacy, with rare labels being harder to predict. The code
on GitHub [@kNNSS] accompanying this paper contains the abridged
version of two datasets, one from male Trace Amine-Associated Receptor
1 (TAAR1) knockout mice described in detail in @schwartz2018
(\autoref{fig:efficM}) and the other from male Sprague-Dawley rats
collected in the Sleep Neurobiology Laboratory at SRI International
(\autoref{fig:efficR}). The rodents in both datasets received an oral
dosing of a water-based vehicle solution.

Three labels were predicted: wake (W), non-REM
sleep (NR), and REM sleep (R), and prediction efficacy was
calculated. (However, note that any number of stages can be trained
depending on how elaborate the manual scoring is.) The model was first
used to train a single classifier merging training data from all
animals (\autoref{fig:efficM} A, \autoref{fig:efficR} A), then
individual models were trained, one for each animal
(\autoref{fig:efficM} B, \autoref{fig:efficR} B). The GitHub repository
includes additional information on prediction accuracy, including
detailed values of true and false positive rates, as well as a method
to deal with imbalanced data.

![Estimation of prediction accuracy for the transgenic mouse data. For
 each state (wake -- W, non-REM -- NR, REM -- R) and animal (points on
 plots) true and false positive rates are calculated. Red crosses
 denote mean and SEM.  In A, training data was merged and one single
 classifier was trained to predict sleep stages of all animals. In B,
 an individual classifier was trained for each animal
 separately.\label{fig:efficM}](fig2.png)

State labels were predicted the same way for the rat data (the same
set of GitHub scripts were run) and prediction accuracy represented on
\autoref{fig:efficR} shows very similar results.

![Estimation of prediction accuracy for the rat data. Prediction and
 figure set up as in
 \autoref{fig:efficM}.\label{fig:efficR}](fig3.png)

# Acknowledgments

TK, DV, and DLB were full time employees and shareholders of Pfizer
Inc.\ during development of this software package. This work was supported
by Pfizer Inc.\ and SRI International.

# Author contributions

DV and TK developed the software, SM, MS, TSK, and DLB contributed data for
development and testing, all authors took part in debugging and testing the
software, and all authors wrote or contributed to writing the manuscript.

# References
