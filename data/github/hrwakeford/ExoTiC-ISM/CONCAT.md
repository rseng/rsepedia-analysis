<!-- PROJECT SHIELDS -->
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02281/status.svg)](https://doi.org/10.21105/joss.02281)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3923986.svg)](https://doi.org/10.5281/zenodo.3923986)
[![MIT License][license-shield]][license-url]
![Python version][python-version-url]

<img src="logo.png" align="left" />

# ExoTiC-ISM
**Exoplanet Timeseries Characterisation - Instrument Systematic Marginalisation**

This code performs Levenberg-Marquardt least-squares minimisation across a grid of pseudo-stochastic instrument systematic models to produce marginalised transit parameters given a lightcurve for a specified wavelength range.

This was developed and tested for data from Wide Field Camera 3 (WFC3) on the Hubble Space Telescope (HST), specifically with the G141 spectroscopic grism, as published in [Wakeford et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016ApJ...819...10W/abstract). This method can also be applied to the WFC3 IR G102 grism, and UVIS G280 grism by selecting the correct parameters.
Future work includes plans to extend this to Space Telescope Imaging Spectrograph (STIS) instrument data, and eventually data from the James Webb Space Telescope (JWST).

This code follows the method outlined in [Wakeford et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016ApJ...819...10W/abstract), using marginalisation across a stochastic grid of 50 polynomial models (see bottom for full grid). 
These 50 instrument systematic models contain a combination of corrective factors for likely HST systematics. These include a linear trend in time across the whole lightcurve, accounting for HST breathing effects caused by thermal changes in the telescope with up to a 4th order polynomial, and correcting for positional shifts of the target spectrum on the detector fitting up to a 4th order polynomial. See [Wakeford et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016ApJ...819...10W/abstract) section 2.2 for details and Table 2 therein for the full grid of systematic models included. 

The evidence (marginal liklihood) is calculated from the AIC for each model when fit with the data and converted to a normalised weighting that is used to marginalise each of the global fit parameters. See equations 15 and 16 in [Wakeford et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016ApJ...819...10W/abstract) to marginalise over the parameters and their uncertainties.

Additional functions in this package that can be utilized independently include the following:
- The program makes use of the analytic transit model in [Mandel & Agol (2002)](https://ui.adsabs.harvard.edu/abs/2002ApJ...580L.171M/abstract) that has been translated into python and can be used independently to fit any transit lightcurve once exotic-ism has been installed. 
- It utilizes Levenberg-Marquardt least squares minimisation using [Sherpa](https://sherpa.readthedocs.io/en/latest/), a Python package for modeling and fitting data. 
- The transit model uses a 4-parameter limb darkening law, as outlined in [Claret (2010)](https://ui.adsabs.harvard.edu/abs/2000A%26A...363.1081C/abstract) and [Sing (2010)](https://ui.adsabs.harvard.edu/abs/2010A%26A...510A..21S/abstract) using 1D Kurucz stellar models (provided on install of this package) or 3D stellar models for a smaller subset of parameters from [Magic et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015A&A...573A..90M/abstract).

This package was built from the original IDL code used for the analysis in [Wakeford et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016ApJ...819...10W/abstract), initially translated by Matthew Hill and then further adapted and transformed into a full astronomy Python package with the help of Iva Laginja.

Note how this package is not distributed through pip or conda, so you will  always need to clone it if you want to work with it.


## Table of contents

* [Supported instruments and gratings](#supported-instruments-and-gratings)
* [Quickstart](#quickstart)
  * [Clone the repo and create conda environment](#clone-the-repo-and-create-conda-environment)
  * [Set up local configfile](#set-up-local-configfile)
  * [Run the main script](#run-the-main-script)
* [Full setup](#full-setup)
  * [Prerequisites](#prerequisites)
  * [Configuration file](#configuration-file)
  * [Output data](#output-data)
  * [Changing input data and/or input parameters](#changing-input-data-andor-input-parameters)
  * [Input filenames](#input-filenames)
  * [The Systematic Model Grid](#the-systematic-model-grid)
  * [Running the testing suite](#running-the-testing-suite)
* [About this repository](#about-this-repository)
  * [Contributing and code of conduct](#contributing-and-code-of-conduct)
  * [Authors](#authors)
  * [Licence](#license)
  * [Acknowledgments](#acknowledgments)
  

## Supported instruments and gratings
Current supported instruments and gratings are:  
**HST** *WFC3* IR/G102, IR/G141 grisms


##  Quickstart

*This section will you give all the necessary terminal commands to go from opening our GitHub page in the browser to having 
reduced results of the template data on your local machine. For a more thorough description of the individual steps, please continue to the section 
**Prerequisites** and beyond.*

For a tutorial on ExoTiC-ISM, please see `notebooks` --> `tutorials` --> `1_Intro-tutorial.ipynb`.

We assume that you have `conda` and `git` installed and that you're using `bash`.

### Clone the repo and create conda environment

- Navigate to the directory you want to clone the repository into:  
```bash
$ cd /User/<YourUser>/repos/
```

- Clone the repository:  
```bash
$ git clone https://github.com/hrwakeford/ExoTiC-ISM.git
```
or use SSH if that is your preferred way of cloning repositories:
```bash
$ git clone git@github.com:hrwakeford/ExoTiC-ISM.git
```

- Navigate into the cloned repository:  
```bash
$ cd ExoTiC-ISM
```

- Create the `exoticism` conda environment:  
```bash
$ conda env create --file environment.yml
```

- Activate the environment:
```bash
$ conda activate exoticism
```

- Install the package into this environment, in editable mode:
```bash
$ python setup.py develop
```

### Set up local configfile

- Go into the code directory:  
```bash
$ cd exoticism
```

- Copy the file `config.ini` and name the copy `config_local.ini`.

- Open your local configfile `config_local.ini` and edit the entry `[data_paths][local_path]` to point to your local repo clone that you just created, e.g.:  
```ini
[data_paths]
local_path = /Users/<YourUser>/repos/ExoTiC-ISM
```

- In the same file, define with `[data_paths][output_path]` where your output data should be saved to, e.g.:  
```ini
[data_paths]
...
output_path = /Users/<YourUser>/<path-to-data>
```

### Run the main script

- Activate the conda environment you just created:
```bash
$ conda activate exoticism
```

- Run the marginalisation on the demo data from the template:  
```bash
$ python marginalisation.py
```

The script takes a short while to run and will output messages to the terminal and save the final data to the path you 
specified under `[data_paths][output_path]` in your `config_local.ini`!

## Full setup

### Prerequisites

This is not an installable package (yet), so you will need to clone it if you want to work with it.
Sherpa is distributed for Mac and Linux, this means Windows users will have to use a Linux virtual machine or find an alternative solution. 

We highly recommend the usage of the package and environment manager [Conda](https://docs.conda.io/projects/conda/en/latest/index.html), 
which is free and runs on Windows, macOS and Linux. We have included an [environment](environment.yml) file in our repository 
from which you can directly build a new conda environment in which we have tested our package. We developed and tested our 
 package with **Python 3.7.3** in **conda 4.6.7**.
 
After cloning the repository, run

```bash
$ conda env create --file environment.yml
```

to build the environment, or optionally

```bash
$ conda env create --name <myEnvName> --file environment.yml
```

to give the environment your own name.

The last step is to install the `exoticism` package into your newly created environment. We do currently not support a
plain install ($ python setup.py install), instead please install it in editable mode:
```bash
$ python setup.py develop
```


### Configuration file

The main configuration file is `config.ini`, which holds all of your simulation parameters. This specific file,
however, is version controlled, and the paths to local directories will get messed up if you push or pull this
file; you might also lose the changes you made to the parameters. This is why `config.ini` is initially supposed to be used as a **template**.

In order to make it work for you, copy `config.ini` and rename the copy to `config_local.ini`. In this **local configfile**, 
you can set all your parameters, and it will override the `config.ini` at runtime. Whichever configfile is used in the end, 
the version controlled one or the local one, a copy of it is always saved together with the output data. In the case you 
want to version control the configfile you use, we recommend that you **fork** the repository and simply use the `config.ini` file directly.

**The configfile** has the following structure, except here we added some extra comments for clarity:
```ini
[data_paths]
local_path = /Users/MyUser/repos/ExoTiC-ISM           ; your global path to the repo clone
input_path = ${local_path}/data                       ; global path to the input data, defaults to template data in repo
output_path = /Users/MyUser/outputs                   ; global path ot the output directory 
run_name = testing                                    ; suffix for output data directory

[setup]
data_set = W17                                   ; data selection; refers to section in configfile
instrument = WFC3
grating = G141
grid_selection = fit_time
ld_model = 3D                     ; 3D or 2D limb darkening model
plotting = True
report = True

[smooth_model]
resolution = 0.0001
half_range = 0.2


; Stellar and planet system parameters - make a new section for each new data set

[W17]
lightcurve_file = W17_${setup:grating}_lightcurve_test_data.txt         ; lightcurve data file
wvln_file = W17_${setup:grating}_wavelength_test_data.txt               ; wavelength data file
rl = 0.12169232                                             ; Rp/R* estimate - the transit depth
epoch = 57957.970153390                                     ; in MJD
inclin = 87.34635                                           ; inclination in deg
ecc = 0.0                                                   ; eccentricity in deg
omega = 0.0                                                 ; deg
Per = 3.73548535                                            ; planet period in days
aor = 7.0780354                                             ;a/r* (unitless) --> "distance of the planet from the star (meters)/stellar radius (meters)"

; limb darkening parameters
metallicity = -1.0                ; stellar metallicity
Teff = 6550                       ; stellar effective temperature
logg = 4.5                        ; log gravity of star

[constants]
dtosec = 86400                    ; conversion factor from days to seconds
HST_period = 0.06691666           ; Hubbe Space Telescope period in days
```

### Output data

The relevant data files and plots from your run, together with the used configfile, will be saved to the directory you specify under **`output_path`** in your 
local configfile. The results of each new run will be saved in a subdirectory under `[data_paths] -> output_path` that is labelled with a time stamp, the
name of the stellar system data and a custom suffix, which you set in the configfile.

### Changing input data and/or input parameters

We provide demo data for the exoplanet WASP-17b, which is one of the datasets analyzed in [Wakeford et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016ApJ...819...10W/abstract).
Please refer to section "Supported instruments and gratings" for a list of currently supported instruments and gratings. If you want to perform the marginalisation on a different 
transit dataset, you have to point the configfile to your input data folder and also update the planetary parameters by adding a new section to the configfile.

#### Input filenames

Due to the structure of the configfile (see above), we follow this naming convention for input files:  
star + grating + arbitrary string + .txt  
*E.g.: "W17_G141_lightcurve_test_data.txt"*

The star and grating can then be set once in the config section '[setup]', while the full filename needs to be added in 
the respective stellar and planetary parameters section, with a placeholder for the grating name.  
*E.g.: "W17_${setup:grating}_lightcurve_test_data.txt"*

All header lines should start with '#'.

### The Systematic Model Grid
<img src="Systematic_model_tabel.png" align="left" />
This table shows the functional form of the systematic models as presented in Wakeford et al. (2016). Each check mark shows which of the parameters are thawed in the model and fit to the data in combination. The grid contains corrections for a linear trend in time across the whole observation, corrections for thermal variations on the time scale of a HST orbit around the Earth, and positional shifts of the observed spectrum on the detector.


### Running the testing suite
Please see [here](CONTRIBUTING.md/#running-the-test-suite) how to run our units tests locally.


## About this repository

### Contributing and code of conduct

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for contribution guidelines, and the process for submitting issues and pull requests to us.
Please also see our [CODE OF CONDUCT](CODE_OF_CONDUCT.md).

If you use this code in your work, please find citation snippets to give us credits with in [CITATION.txt](CITATION.txt).

### Authors

* **Hannah R. Wakeford** - *Method author* - [@hrwakeford](https://github.com/hrwakeford)
* **Iva Laginja** - *Turning the code into a functional Python repository* - [@ivalaginja](https://github.com/ivalaginja)
* **Matthew Hill** - *Translation from IDL to Python* - [@mattjhill](https://github.com/mattjhill)

### License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.txt) file for details.

### Acknowledgments

* Tom J. Wilson for statistical testing of the code
* Matthew Hill for a functional translation from IDL to Python
* Iva Laginja for finding `Sherpa`, making the clunky `mpfit` dispensable
* The [`Sherpa` team](https://github.com/sherpa/sherpa), providing a fantastic package and answering fast to GitHub issues
* This work is  based  on  observations  made  with  the  NASA/ESA Hubble Space Telescope, HST-GO-14918, that were obtained at the Space Telescope Science Institute, which isoperated by the Association of Universities for Research in Astronomy, Inc.  


<!-- MARKDOWN LINKS & IMAGES -->
[license-shield]: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square
[license-url]: https://choosealicense.com/licenses/mit
[python-version-url]: https://img.shields.io/badge/Python-3.6-green.svg?style=flat
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
 advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
 address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
 professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at stellarplanet+exoticism@gmail.com. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
# Contribution guidelines

## General guidelines
Hello and thank you for contributing to our project! Go ahead open a new issue or new pull request (PR) for bugs,
feedback, or new features you would like to see by following these guidelines.
Leave a comment anywhere and we will be happy to assist. New contributions and contributors are very welcome!
We recommend that you **fork** the repository and submit your requests according to these guidelines:

1) Search previous issues and pull requests as your question or suggestion might be a duplicate.
2) Open a separate issue or pull request for each suggestion or bug report.
3) Pull requests, issues and commits should have a useful title.
4) Give every PR a description and open it against **develop**, the default branch of the repository.
5) Every PR that adds or changes a feature should also add or change an/the according **test**, if applicable, as well
as **documentation**.

## Running the test suite
We use the [pytest](https://docs.pytest.org/en/stable/) framework for the unit and regression tests on our repository.
To run the tests locally to make sure that they all pass, you just have to navigate into the top-level repository
directory, activate your `exoticism` conda environment (`$ conda activate exoticism`) and run:
```bash
$ pytest
```
This will run all tests and give you a test report with passing and failing tests. The package `pytest` is automatically
installed in the `exoticism` conda environment when creating it from our `environment.yml` file as described in the
[Quickstart section of our README](README.md/#quickstart).

## Code of conduct
This package follows the Contributor Covenant [Code of Conduct](CODE_OF_CONDUCT.md) to provide a welcoming community to everybody.
---
title: 'ExoTiC-ISM: A Python package for marginalised exoplanet transit parameters across a grid of systematic instrument models'
tags:
  - Python
  - astronomy
  - exoplanets
  - transit spectroscopy
  - marginalisation
  - instrument systematics
authors:
  - name: Iva Laginja
    orcid: 0000-0003-1783-5023
    affiliation: "1, 2, 3"
  - name: Hannah R. Wakeford
    orcid: 0000-0003-4328-3867
    affiliation: 4
affiliations:
 - name: Space Telescope Science Institute, Baltimore, USA
   index: 1
 - name: DOTA, ONERA, Université Paris Saclay, F-92322 Châtillon, France
   index: 2
 - name: Aix Marseille Université, CNRS, CNES, LAM, Marseille, France
   index: 3
 - name: School of Physics, University of Bristol, HH Wills Physics Laboratory, Tyndall Avenue, Bristol BS8 1TL, UK
   index: 4
date: 27 May 2020
bibliography: paper.bib

---

# Science background

## Transit spectroscopy of exoplanets

There have been a slew of planet detections outside our own solar system over the past two decades and several 
characterisation methods can be used to determine their chemical compositions. One of them is transit spectroscopy. 
With this technique, astronomers measure the star light passing through an exoplanet's atmosphere while it is 
transiting in front of its host star. Imprinted on this light are the absorption signatures of different 
materials - atoms and molecules in the gas phase, or solid or liquid aerosols - in the transiting planet's atmosphere. 
Using a spectrograph the flux is recorded as a function of wavelength, allowing scientists to construct 
absorption/transmission spectra, with the goal of identifying the chemical composition of the atmosphere.

There are many different chemical components that transit spectroscopy can reveal in an exoplanet, but a majority of the 
exoplanets studied via transmission spectroscopy are on close-in orbits around their 
stars lasting only several days, most of them giant Jupiter- or Neptune-sized worlds. For these giant, 
close-in exoplanets, the most dominant source of 
absorption will be from water vapour, which is expected to be well-mixed throughout their atmosphere. H$_2$O has 
strong absorption in the near-infrared (IR) with broad peaks at 0.9, 1.4, 1.9, and 2.7 $\mu$m. However, these 
absorption features cannot be measured from the ground as the Earth's atmosphere, filled with water vapour, gets in 
the way. To measure H$_2$O in the atmospheres of exoplanets, astronomers use the Hubble Space Telescope's Wide Field 
Camera 3 (HST WFC3) infrared capabilities to detect the absorption signatures of H$_2$O at 0.9 $\mu$m with the G102 
grism, and at 1.4 $\mu$m with the G141 grism [e.g. @kreidberg2015; @sing2016; @wakeford2017; @wakeford2018; @spake2018].

## Calibration of instrument systematics with marginalisation

While all telescopes aim at recording the light of distant worlds as accurately as possible, every 
instrument will have its own signature that it superimposes on the collected signal. These effects that contaminate our 
data are called "instrument systematics" and need to be calibrated out before an observation can be interpreted 
truthfully. All the different systematics that influence the result recorded by an instrument are combined into 
a systematic instrument model to be used during the calibration. These models are not always obvious and different 
authors often suggest differing systematic instrument models being applied to data from the same instrument, often 
favouring models that work particularly well for any given data set. This makes it hard to compare data sets to each 
other, yielding moderately different results for the parameters of a transiting planet when a different systematic 
model was applied to the data.

The solution that was applied to WFC3 data in @wakeford2016 performs a marginalisation across a grid of systematic 
models that take different corrections across an exoplanet transit data set into account. Following the method proposed 
by @gibson2014, a Levenberg-Marquardt least-squares minimisation is performed across all systematic models, which yields 
a set of fitted transit parameters for each systematic model. We then use the resulting Akaike Information 
Criterion (AIC) to calculate each model’s evidence (marginal likelihood) and normalised weight. These weights are then 
used to calculate the marginalised fit parameters, leading to results that will not depend as heavily on the individual 
choice of systematic model. Finally, performing this for each lightcurve constructed at each wavelength from 
the measured spectrum results in the measured transmission spectrum of the exoplanet.

# The ``ExoTiC-ISM`` package

## Functionality

``ExoTiC-ISM`` (Exoplanet Timeseries Characterisation - Instrument Systematic Marginalisation) is an open-source Python 
package that computes the transit depth from a timeseries lightcurve, while sampling a grid of pseudo-stochastic models 
to account for instrument based systematics that may impact the measurement, following the method proposed by 
@gibson2014 and implemented by @wakeford2016. While there are a number of Python solutions to create and fit transiting 
planet light curves, ``ExoTiC-ISM`` is a lightcurve fitting tool that focuses particularly 
on the statistical method of marginalisation. It allows for a transparent method of data analysis of systematics
impacting a measurement. While other methods, such as Gaussian processes (GP), can account for the likelyhood of 
systematics impacting your measurement, these methods can typically not easily determine which systematics are the most 
important, and which combination of systematics is specifically affecting your data set. ``ExoTiC-ISM`` allows you to 
evaluate a grid of instrument systematic models to obtain the needed information on the dominant systematics enabling 
you to design the next observation to be more efficient and precise. As the 
authors of the original method paper state [@wakeford2016]: “The use of marginalisation 
allows for transparent interpretation and understanding of the instrument and the impact of each systematic [model] 
evaluated statistically for each data set, expanding the ability to make true and comprehensive comparisons between 
exoplanet atmospheres.”

The currently implemented instrument systematic grid is composed of a series of 49 polynomial functions 
that are specifically designed to account for systematics associated with the detectors on HST WFC3 [@wakeford2016]. 
However, this can be adapted to other instruments.
The package performs the Levenberg-Marquardt least-squares minimisation across all models with the 
``sherpa`` package [@sherpa.v4.11.0; @sherpa_paper_1; @sherpa_paper_2] for modeling and fitting data, and then 
calculates the AIC and normalised weight to 
marginalise over the fit parameters (e.g. transit depth $rl$, inclination $i$, a/R$_*$, center of transit time) using 
each systematic model. This method is different from evaluating each systematic model independently 
and selecting the “best” one purely by minimising the scatter of its residuals as that would not include a 
penalisation for increased model complexity nor information from similarly likely systematic corrections.

The original code was written in IDL, which was used to publish marginalised transit parameters for five different 
exoplanets [@wakeford2016] observed in the IR with the G141 grism on HST's WFC3. The ``ExoTiC-ISM`` package described 
in this paper implements a marginalisation for that same grism and extends its functionality to the G102 grism, which 
uses the same grid of systematic models [see results by @wakeford2017; @wakeford2018].

## Dependencies and usage

``ExoTiC-ISM`` is written in Python with support for Python 3.6 and 3.7 on MacOS and Linux. It makes use of the packages 
``numpy`` [@numpy1; @numpy2], ``astropy`` [@astropy2013; @astropy2018], ``pandas`` [@pandas-paper; @pandas-zenodo], 
``matplotlib`` [@matplotlib; @matplotlib-zenodo], ``sherpa`` [@sherpa.v4.11.0; @sherpa_paper_1; @sherpa_paper_2], as well as some custom functions, 
like an implementation of the transit function by @mandel2002. It applies a 4-parameter limb darkening law as outlined 
in @claret2000 and @sing2010, using either the 1D Kurucz stellar models or the 3D stellar atmosphere models by @magic2015.

The required inputs for the analysis are two text files that contain the uncorrected lightcurve of the observed object, and a 
wavelength array. Input parameters are rendered from a ``config.ini`` file and we provide an ``environment.yml`` file 
to build a ``conda`` environment to run the package in. Development in Python and hosting the repository on GitHub 
will facilitate the usage of the package by researchers, as well as further functional development; an introductory 
tutorial is provided in the form of a Jupyter Notebook.

## Outlook

While its current capabilities are limited to WFC3 data taken with the G141 and G102 grism, the package’s 
functionality will be extended to the UVIS G280 grism [@wakeford2020] and the G430L and G750L gratings of the Space Telescope 
Imaging Spectrograph (STIS) on HST. This will lay the groundwork for the envisioned future extension to implement 
systematic grids for select instruments on the James Webb Space Telescope (JWST) and obtain robust transit spectra 
for JWST data.

# Acknowledgements

The authors would like to acknowledge Matthew Hill who translated the first part of the IDL version of this code to 
Python. We also acknowledge contributions by Tom J. Wilson to the statistical testing of the code. 
We also thank the Sherpa team for their fast and detailed responses to questions we had during the 
implementation of their package. This work is based on observations made with the NASA/ESA Hubble Space Telescope, 
HST-GO-14918, that were obtained at the Space Telescope Science Institute, which is operated by the Association of 
Universities for Research in Astronomy, Inc.

# References
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
