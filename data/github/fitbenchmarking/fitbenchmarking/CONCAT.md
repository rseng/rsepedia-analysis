[![Build Status](https://img.shields.io/github/workflow/status/fitbenchmarking/fitbenchmarking/Build%20and%20Publish?style=flat-square)](https://github.com/fitbenchmarking/fitbenchmarking/actions/workflows/release.yml)
[![Tests Status](https://img.shields.io/github/workflow/status/fitbenchmarking/fitbenchmarking/Tests?label=tests&style=flat-square)](https://github.com/fitbenchmarking/fitbenchmarking/actions/workflows/main.yml)
[![Documentation Status](https://img.shields.io/readthedocs/fitbenchmarking?style=flat-square)](https://fitbenchmarking.readthedocs.io/en/latest)
[![Coverage Status](https://img.shields.io/coveralls/github/fitbenchmarking/fitbenchmarking.svg?style=flat-square)](https://coveralls.io/github/fitbenchmarking/fitbenchmarking)
![Windows Supported](https://img.shields.io/badge/win10-support-blue.svg?style=flat-square&logo=windows)
![Ubuntu Supported](https://img.shields.io/badge/18.04-support-orange.svg?style=flat-square&logo=ubuntu)
[![Chat](https://img.shields.io/badge/chat-CompareFitMinimizers-lightgrey.svg?style=flat-square&logo=slack)](https://slack.com/)
# FitBenchmarking

FitBenchmarking is an open source tool for comparing different minimizers/fitting frameworks. FitBenchmarking is cross platform and we support Windows, Linux and Mac OS. For questions, feature requests or any other inquiries, please open an issue on GitHub, or send us an e-mail at support@fitbenchmarking.com.

- **Installation Instructions:** https://fitbenchmarking.readthedocs.io/en/latest/users/install_instructions/index.html
- **User Documentation & Example Usage:** https://fitbenchmarking.readthedocs.io/en/latest/users/index.html
- **Community Guidelines:** https://fitbenchmarking.readthedocs.io/en/latest/contributors/guidelines.html
- **Automated Tests:** Run via GitHub Actions, https://github.com/fitbenchmarking/fitbenchmarking/actions, and tests are documented at https://fitbenchmarking.readthedocs.io/en/latest/users/tests.html

The package is the result of a collaboration between STFC’s Scientific Computing Department and ISIS Neutron and Muon Facility and the Diamond Light Source. We also would like to acknowledge support from:

* EU SINE2020 WP-10, which received funding from the European Union’s Horizon2020 research and innovation programme under grant agreement No 654000.
* EPSRC Grant EP/M025179/1  Least Squares: Fit for the Future.
* The Ada Lovelace Centre (ALC). ALC is an integrated, cross-disciplinary data intensive science centre, for better exploitation of research carried out at our large scale National Facilities including the Diamond Light Source (DLS), the ISIS Neutron and Muon Facility, the Central Laser Facility (CLF) and the Culham Centre for Fusion Energy (CCFE).
---
title: '`FitBenchmarking`: an open source `Python` package comparing data fitting software'
tags:
  - Python
  - fitting
  - non-linear least squares
authors:
  - name: Anders Markvardsen
    affiliation: 1
  - name: Tyrone Rees
    affiliation: 1
  - name: Michael Wathen
    affiliation: 1
  - name: Andrew Lister
    affiliation: 1
  - name: Patrick Odagiu
    affiliation: 1
  - name: Atijit Anuchitanukul
    affiliation: 1
  - name: Tom Farmer
    affiliation: 1
  - name: Anthony Lim
    affiliation: 1
  - name: Federico Montesino
    affiliation: 1
  - name: Tim Snow
    affiliation: 2
  - name:  Andrew McCluskey
    affiliation: 2
affiliations:
 - name: Science and Technology Facilities Council, Rutherford Appleton Laboratory, Harwell Campus, Didcot, Oxfordshire, OX11 0QX
   index: 1
 - name: Diamond Light Source Ltd, Diamond House, Harwell Campus, Didcot, Oxfordshire, OX11 0DE
   index: 2
date: October 2020
bibliography: paper.bib
---
# Summary

Fitting a mathematical model to data is a fundamental task across all scientific disciplines. [`FitBenchmarking`](https://fitbenchmarking.com/) has been designed to help:

* Scientists, who want to know the best algorithm for fitting their data to a given model using specific hardware.
* Scientific software developers, who want to identify the best fitting algorithms and implementations. This allows them to recommend a default solver, to see if it is worth adding a new minimizer, and to test their implementation.
* Mathematicians and numerical software developers, who want to understand the types of problems on which current algorithms do not perform well, and to have a route to expose newly developed methods to users.

Representatives of each of these communities have got together to build `FitBenchmarking`. We hope this tool will help foster fruitful interactions and collaborations across the disciplines.

![Benchmarking paradigm: associating fitting problems represented in individual scientific software packages (top cycle) to optimization software packages (bottom cycle), and bringing these closer together. \label{fig:concept}](figures/FitBenchmarkingConcept.png){width=60%}

`FitBenchmarking` is easy to install via `pip` and our [documentation](https://fitbenchmarking.com/) guides users through the installation of some external packages we support. We provide several data sets from a range of applications and adding new data in these formats is as easy as dropping the data into a new folder. The data and fitting packages currently supported are shown in Figure \ref{fig:concept}. A key part of `FitBenchmarking` is the ease of which a user, with a basic knowledge of `Python`, can add new fitting software, data formats and different fitting comparison output metrics.

# State of the field

Fitting data to models is a form of optimization, and 
`CUTEst` [@cutest], and its predecessors, has been the standard tool to
benchmark optimization packages for some time. `CUTEst` can benchmark any problem
written in a custom `SIF` format.  However, only the hooks to run the same problem are
provided, the user must provide their own data analysis.  Tools such as
`Paver` [@paver], part of the COIN-OR initiative, can be used
alongside `CUTEst` (or other tools) for this purpose.
The packages `Olympus` [@olympus] and `Benchopt` [@benchopt] have been recently
developed as benchmarking and analysis frameworks for optimization problems.
`Olympus` is designed for experiment planning and provides analytic benchmark problems,
experimental datasets, and emulated datasets, but could be adapted to be applied to
any optimization (or data-fitting) problem.
`Benchopt`, on the other hand, is currently primarily used to benchmark data fitting
using a range of cost functions.  `Benchopt` ships with a limited number of example data
sets, but it is well documented how to write new benchmarks with custom data and
objective functions.


# Statement of need

While there is some overlap between `FitBenchmarking` and the rest of the field,
what makes our software unique is:

* It is designed to interface directly to the source of data, be that a scientific
  software package or an academic data set.  Our `parser` class can be extended
  to make it clear what a developer needs to do to get data into `FitBenchmarking`.
* While being easy to extend using new software, or new data from currently supported
  packages, `FitBechmarking` ships with open datasets that all can use for testing.
* FitBenchmarking tests implementations of algorithms, not just algorithms.
  A growing number of optimization packages that can be used for data fitting are
  supported, and it is straightforward to extend our `controller` class to add new
  software.
* `FitBenchmarking` performs its own data processing and analysis and, if needed,
  the output generated can be customized for new data sets and/or minimizers.  

As far as we are aware, `FitBenchmarking` is the only package that is designed
specifically to interface directly with optimization  packages and individual
scientific software packages to test different implementations of fitting algorithms.
`FitBenchmarking` originally started as a tool to benchmark fitting algorithms in the data reduction package `Mantid` [@mantid], which is used to process neutron scattering and muon spectroscopy data. `FitBenchmarking` has since been significantly extended to take data and models from other real world applications and data analysis / modelling / treatment packages, such as `SasView` [@sasview] and `CUTEst` [@cutest]. It fits models to the data by using a range of data fitting and nonlinear optimization software packages, and present comparisons through a variety of different metrics. These include comparison tables and performance profile plots.

`FitBenchmarking` compares how different fitting algorithms perform for the same data, model and initial guess. The best parameters for the model are found by solving a nonlinear least-squares problem, which can either be solved using a dedicated optimisation software package or using a fitting algorithm implementation within a scientific software package. Figure \ref{fig:sample} displays a data set from `FitBenchmarking` where the crosses are the data points and the two curves are the fits found by two optimization algorithms implemented in `GSL` [@gsl]. From Figure \ref{fig:sample}, it is clear that the solution given by lmsder is better. As the volume of data increases, and we do more and more scientific analysis algorithmically, it is increasingly important that we apply the best available algorithm for a given category of fitting problems. `FitBenchmarking` generates HTML output that makes it easy to compare minimizers on a given problem set.

![A sample fit: this problem is shipped with `FitBenchmarking`. The data was collected from an instrument named VESUVIO at the ISIS Neutron and Muon Source and has a difficult initial guess. \label{fig:sample}](figures/nmsimplex2_fit_for_EVS14188-90_processed_Gaussian_peaks_1_1.png){width=70%}

`FitBenchmarking` will help the scientist make an informed choice by comparing runtime and accuracy of all available minimizers, on their specific hardware, on problems from their science area.

`FitBenchmarking` will help the scientific software developer ensure that the most robust and quickest algorithms for the type of data analysis they support are available in their software.

`FitBenchmarking` will help mathematicians see what the state of the art is, and what kinds of data are problematic. It will give them access to real data, and will give a route for novel methods to quickly make it into production.

# Acknowledgements

We would like to acknowledge funding support from:

* European Union’s Horizon2020 research and innovation programme, EU SINE2020 WP-10,
* EPSRC Grant EP/M025179/1 -- Least Squares: Fit for the Future.
* The Ada Lovelace Centre (ALC).

We would also like to thank Nick Draper, Roman Tolchenov, Nick Gould and Jaroslav Fowkes for their helpful comments and advice.

# References
#### Description of Work

Fixes


#### Testing Instructions

1.
2.
3.

Function: Does the change do what it's supposed to?

Tests: Does it pass? Is there adequate coverage for new code?

Style: Is the coding style consistent? Is anything overly confusing?

Documentation: Is there a suitable change to documentation for this change?
---
name: Test
about: Specify a test that needs to be implemented
title: ''
labels: Testing

---

**Which module and class/method/function does this relate to?**
Please provide the module and class, method or function name(s) that the tests will apply to.

**What aspect requires additional tests?**
Please provide a description of the specifics of the test - what is the behaviour that the tests will inspect.  Ideally this will include a description of the setup (e.g. Pytest fixtures) and the expected results.

**Is this a unit, system or functional test?**
Simply state what type of test you are expecting is required.

**Additional context**
Add any other context about the tests here.
---
name: Documentation
about: Suggest a improvement to user, developer or other documentation
title: ''
labels: Documentation

---

**Description of the documentation**
A description of what area can benefit from better documentation.

**Additional details**
Add any other context or screenshots about this request.
---
name: Bug report
about: Report an error which requires fixing
title: ''
labels: Bug

---

**Description of the error**
A clear and concise description of what the problem is, including steps to reproduce it and the environment you are running FitBenchmarking in.

**Describe the expected result**
What is the result you expect from running the steps described above?

**Describe the actual result**
What was the actual result.

**Suggested fix**

**Additional details**
Add any other context or screenshots about the feature request here.
---
name: Feature request
about: Suggest a new feature
title: ''
labels: Enhancement

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
This folder contains problem definition files each defining a
fit benchmarking problem that minimizers can be benchmarked against.

More specifically a folder or sub-folder of this directory contain a problem set.

See https://fitbenchmarking.readthedocs.io and section on Problem Definition Files
for the formats supported.

A data file name specified within a definition file is recommended to be stored in a sub-folder named `data_files` relative to the location of the definition file.

Examples of Problem Definition folders include (please note this list is stadily changes and hence below may get out of sync from time to time):

* SAS_modelling : fitting problems specific relevant to fitting SAS (Small Angle Scattering) data
* Muon : generic fitting problems relevant to fitting Muon data collected at a Muon facility such as ISIS Neutron and Muon Facility
* Neutron : generic fitting problems relevant to fitting Neutron data collected at a Neutron facility such as ISIS Neutron and Muon Facility
* NIST : set of made up fitting problem (not against measured data with error bars) as described [here](https://www.itl.nist.gov/div898/strd/nls/nls_main.shtml)
* CUTEes : fitting problems relevant to this tool included in [CUTEet](http://epubs.stfc.ac.uk/bitstream/9327/RAL-TR-2013-005.pdf)

This folder contains scripts that have been used at various times to create more generic type data formatted files from more software specific type ones,
thereby allowing problem definition files referencing such data to be made easier available for benchmarking across minimizer libraries.
Of course, this is not aways feasible, but where it is it is recommended.
A particular simple format is column ascii format, for 2D data this is X and Y columns with an optional E column,
where the E column contains the errors of the Y values.This folder includes images etc. used in the main FitBenchmarking [README](../README.md). Other documentation can also be found on the FitBenchmarking [Wiki](https://github.com/fitbenchmarking/fitbenchmarking/wiki).

Furthermore, sphinx documentation is being setup within the source folder and which can be build with Makefile and make.bat. Also this docs is linked up with [Read the Docs](https://readthedocs.org/) and the automatically build sphinx is viewable from https://fitbenchmarking.readthedocs.io
