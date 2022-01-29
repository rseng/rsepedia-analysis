# bfit

<a href="https://pypi.org/project/bfit/" alt="PyPI Version"><img src="https://img.shields.io/pypi/v/bfit?label=PyPI%20Version"/></a>
<img src="https://img.shields.io/pypi/format/bfit?label=PyPI%20Format"/>
<img src="https://img.shields.io/github/languages/code-size/dfujim/bfit"/>
<img src="https://img.shields.io/tokei/lines/github/dfujim/bfit"/>
<img src="https://img.shields.io/pypi/l/bfit"/>
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03598/status.svg)](https://doi.org/10.21105/joss.03598)

<a href="https://github.com/dfujim/bfit/commits/master" alt="Commits"><img src="https://img.shields.io/github/commits-since/dfujim/bfit/latest/master"/></a>
<a href="https://github.com/dfujim/bfit/commits/master" alt="Commits"><img src="https://img.shields.io/github/last-commit/dfujim/bfit"/></a>

[bfit] is a [Python] application aimed to aid in the analysis of β-detected
nuclear magnetic/quadrupole resonance (β-NMR and β-NQR) data taken at [TRIUMF].
These techniques are similar to muon spin rotation ([μSR]) and "conventional"
nuclear magnetic resonance ([NMR]), but use radioactive nuclei as their [NMR]
probe in place of the [muon] or a stable isotope.
The instruments and research program are governed through [TRIUMF]'s [CMMS],
with more information given at <https://bnmr.triumf.ca>.
An overview of instrumentation details and scientific applications of the
β-NMR/β-NQR techniques can be found in several recent journal articles:

- W. A. MacFarlane.
  <i>Implanted-ion βNMR: a new probe for nanoscience</i>.
  <a href="https://doi.org/10.1016/j.ssnmr.2015.02.004">
  Solid State Nucl. Magn. Reson. <b>68-69</b>, 1-12 (2015)</a>.
- G. D. Morris.
  <i>β-NMR</i>.
  <a href="https://doi.org/10.1007/s10751-013-0894-6">
  Hyperfine Interact. <b>225</b>, 173-182 (2014)</a>.

The intended user of [bfit] is anyone performing experiments with or analyzing
data taken from [TRIUMF]'s β-NMR or β-NQR spectrometers - independent of whether
they are a new student, visiting scientist, or someone with decades of experience.
(e.g., someone from the "local" [TRIUMF]/[CMMS]/[UBC] group).
A key goal of the project is to alleviate much of the technical tedium that is
often encountered during any analysis.
More generally, [bfit] has been written to fulfill the following needs:

* Provide the means for quick on-line analyses during beam time.
* Provide a useful and flexible API for refined analyses in [Python],
  to be used in conjunction with [bdata] and the [SciPy] ecosystem.
* Provide an intuitive, user-friendly interface for non-programmers.
* Be easily maintainable and distributable.

## Citing

If you use [bfit] in your work, please cite:

D. Fujimoto, "bfit: A Python Application For Beta-Detected NMR," [J. Open Source Softw. <b>6</b>, 3598 (2021).](https://doi.org/10.21105/joss.03598)

## Useful Links

* [bfit]
  * [Wiki]
    * [API Reference]
    * [API Tutorial]
    * [GUI Tutorial]
* [mudpy]
* [bdata]

## Community Guidelines

* Contributing:
  * Please submit your contribution to [bfit] through the list of
    [Pull Requests]!
* Reporting issues and/or seeking support:
  * Please file a new ticket in [bfit]'s list of [Issues] - I will get an email
    notification of your problem and try to fix it ASAP!

## Installation and Use

### Dependencies

The following packages/applications are needed _prior_ to [bfit] installation:
- [Python] 3.6 or higher: a dynamically typed programming language. [[install](https://wiki.python.org/moin/BeginnersGuide/Download)]
- [Tkinter] : [Python]'s de facto standard GUI package. [[install](https://tkdocs.com/tutorial/install.html)]

### Install Instructions

|  | Command |
|:-- | :--|
From the [PyPI] as user (recommended) | `pip install --user bfit` |
From the [PyPI] as root | `pip install bfit` |
From source | `python3 setup.py install` |

Note that `pip` should point to a (version 3) [Python] executable
(e.g., `python3`, `python3.8`, etc.).
If the above does not work, try using `pip3` or `python3 -m pip` instead.

### Optional Setup

For convenience,
you may want to tell [bfit] where the data is stored on your machine.
This is done by defining two environment variables:
`BNMR_ARCHIVE` and `BNQR_ARCHIVE`.
This can be done, for example, in your `.bashrc` script.
Both variables expect the data to be stored in directories with a particular
heirarchy:

```
/path/
|---bnmr/
|---bnqr/
|-------2017/
|-------2018/
|-----------045123.msr
```

Here, the folders `/path/bnmr/` and `/path/bnqr/` both contain runs
(i.e., `.msr` files) organized into subdirectories by year of aquasition.
In this case, you would set (in your `.bashrc`):

```bash
export BNMR_ARCHIVE=/path/bnmr/
export BNQR_ARCHIVE=/path/bnqr/
```

If [bfit] cannot find the data, it will attempt to download the relavent [MUD]
(i.e., `.msr`) files from <https://cmms.triumf.ca/mud/runSel.html>.
This is the default behaviour for [bfit] installed from [PyPI].

### First Startup 

To launch the GUI from a terminal simply call `bfit`, if this fails, one can also use the alternative syntax `python3 -m bfit`, where `python3` may be replaced with any (version 3) [Python] executable.

### Testing

Testing your installation of [bfit] is accomplished by running `pytest` within the installation folder. Note that some tests, particularly those involving drawing, fail when run as group in this environment, but they should pass on a subsequent attempts: `pytest --lf`. Further testing information can be found [here](https://github.com/dfujim/bfit/wiki/Installation-and-first-startup).

[Python]: https://www.python.org/
[SciPy]: https://www.scipy.org/
[Cython]: https://cython.org/
[NumPy]: https://numpy.org/
[pandas]: https://pandas.pydata.org/
[Matplotlib]: https://matplotlib.org/
[Tkinter]: https://wiki.python.org/moin/TkInter
[PyYAML]: https://pyyaml.org/
[pytest]: https://docs.pytest.org/en/6.2.x/
[tqdm]: https://github.com/tqdm/tqdm
[requests]: https://requests.readthedocs.io/en/master/
[Jupyter]: https://jupyter.org/
[argparse]: https://docs.python.org/3/library/argparse.html

[YAML]: https://yaml.org/
[C]: https://en.wikipedia.org/wiki/C_(programming_language)
[HTTP]: https://en.wikipedia.org/wiki/Hypertext_Transfer_Protocol

[TRIUMF]: https://www.triumf.ca/
[CMMS]: https://cmms.triumf.ca
[MUD]: https://cmms.triumf.ca/mud/
[archive]: https://cmms.triumf.ca/mud/runSel.html
[`data/BNMR/2020/040123.msr`]: https://cmms.triumf.ca/mud/mud_hdrs.php?ray=Run%2040123%20from%20BNMR%20in%202020&cmd=heads&fn=data/BNMR/2020/040123.msr

[PHYSICA]: https://computing.triumf.ca/legacy/physica/
[UBC]: https://www.ubc.ca/
[μSR]: https://en.wikipedia.org/wiki/Muon_spin_spectroscopy
[NMR]: https://en.wikipedia.org/wiki/Nuclear_magnetic_resonance
[muon]: https://en.wikipedia.org/wiki/Muon

[bnmr_1f]: https://gitlab.com/rmlm/bnmr_1f
[bnmr_2e]: https://gitlab.com/rmlm/bnmr_2e
[bnmrfit]: https://gitlab.com/rmlm/bnmrfit
[bnmroffice]: https://github.com/hsaadaoui/bnmroffice
[musrfit]: https://bitbucket.org/muonspin/musrfit
[musrfit documentation]: https://lmu.web.psi.ch/musrfit/user/html/index.html

[mudpy]: https://github.com/dfujim/mudpy
[bdata]: https://github.com/dfujim/bdata

[bfit]: https://github.com/dfujim/bfit
[Pull Requests]: https://github.com/dfujim/bfit/pulls
[Issues]: https://github.com/dfujim/bfit/issues
[PyPI]: https://pypi.org/project/bfit/
[API Reference]: https://github.com/dfujim/bfit/wiki/API-Reference
[API Tutorial]: https://github.com/dfujim/bfit/wiki/API-Tutorial
[GUI Tutorial]: https://github.com/dfujim/bfit/wiki/GUI-Tutorial
[Wiki]: https://github.com/dfujim/bfit/wiki

[ROOT]: https://github.com/root-project/root
[MINUIT]: https://doi.org/10.1016/0010-4655(75)90039-9
[MINUIT2]: https://root.cern/doc/master/Minuit2Page.html
[iminuit]: https://github.com/scikit-hep/iminuit
---
title: 'bfit: A Python Application For Beta-Detected NMR'
tags:
  - Python
  - beta-detected NMR
authors:
  - name: Derek Fujimoto
    orcid: 0000-0003-2847-2053
    affiliation: "1,2"
affiliations:
 - name: Stewart Blusson Quantum Matter Institute, University of British Columbia, Vancouver, BC V6T 1Z4, Canada
   index: 1
 - name: Department of Physics and Astronomy, University of British Columbia, Vancouver, BC V6T 1Z1, Canada
   index: 2
date: 17 May 2021
bibliography: paper.bib
---

# Summary

Beta-detected nuclear magnetic resonance ($\beta$-NMR) measures the beta-decay of probe radioactive nuclei to infer the electromagnetic character of the probe's local environment. Similar to muon spin rotation ($\mu$SR), this technique allows for unique insight of material properties not easily measured by conventional NMR. The [`bfit`] package provides a graphical user interface (GUI) and application programming interface (API) to facilitate the analysis of implanted-ion $\beta$-NMR measurements taken at TRIUMF.

# Background

$\beta$-NMR leverages the parity-violating nuclear weak interaction to measure the spin precession of a ensemble of radioactive probe nuclei [@MacFarlane2015]. These nuclei can either be activated by neutrons or implanted as a foreign species in the form of a low-energy particle beam. Upon decay, the direction of the emitted electron is correlated with the nuclear spin orientation. As with many nuclear and particle physics experiments, the data collected is the counted number of electrons emitted in a given direction. These counts are then histogrammed and processed to yield a signal of interest.

The activation or implantation of the probe nuclei require high-intensity particle beams, restricting the technique to large nationally-supported facilities. Even today, there are only a handful of locations capable of conducting $\beta$-NMR measurements, such as TRIUMF, which is situated in Vancouver, Canada. This facility has been running $\beta$-NMR experiments for the past 20 years, and has developed the Muon Data (MUD) file format [@Whidden1994] as a means of storing $\mu$SR and $\beta$-NMR data.

# Statement of need

At TRIUMF, $\beta$-NMR receives approximately 5 weeks of radioactive beam time per year. As with other large-facility experiments employing particle beams, this data is extremely limited and expensive to generate. Having the tools for rapid on-line analysis is therefore crucial for efficient and informed measurement. Additionally, many of the experimenters using the $\beta$-NMR spectrometer are visiting scientists or students who have little experience with the technical aspects of the measurement.

As with many older science applications, the MUD API is written in C and FORTRAN. These statically-typed and compiled languages are known for their computational efficiency, but are accompanied by long development times, relative to modern languages. In many communities, scientific computing has shifted to languages such as Python: a dynamically-typed and interpreted language. As a result, Python has amassed a massive library of data analysis tools [@Virtanen2020]. The short development time of Python programs is particularly important in the context of scientific analyses, which are typically run only a few times by select individuals. As a result, the development time of the analysis code comprises a large part of the program's effective run time. The aim of this work is to bring this rapid prototyping style of analysis to $\beta$-NMR. To further streamline on-line analyses, [`bfit`] provides an intuitive GUI capable of a moderately high degree of sophistication.

It should be acknowledged that, while a large body of analysis software exists to support $\mu$SR workers (such as WIMDA [@Pratt2000], MANTID [@Arnold2014], and Musrfit [@Suter2012]), $\beta$-NMR does not have a comparably extensive suite of maintained analysis programs. While there have been some recent improvements to this situation [@Saadaoui2018], the analysis required for any non-trivial $\beta$-NMR experiment necessitates the development of new code to meet the individual requirements of each experiment. While such code may employ Musrfit, which is compatible with the MUD file format, this approach may be cumbersome for complex or rapid analyses, and presents a entry high entry barrier for new users. The Python API of [`bfit`] is well suited for addressing these issues.

# Usage and features

The [`bfit`] GUI has three primary functions which are contained in the _Inspect_, _Fetch_, and _Fit_ tabs. The purpose of the _Inspect_ tab (shown below) is to quickly view the file headers and plot the data in order to detect and solve problems as they may arise during measurement. The _Fetch_ tab has been designed to prepare the data for analysis, loading runs in batch and allowing the user to draw and compare each run. The _Fit_ tab provides the tools needed to fit a model to the data, and to view and analyze the result. These tools include global fitting (i.e., sharing fit parameters between data sets), constrained fitting (i.e., constraining a parameter to follow a specific model dependent on the experimental conditions, such as temperature), non-trivial fitting functions specific to pulsed-beam operation (leveraging double exponential integration [@Cook2014]), multiple minimization routines, and more.

![The inspection tab of the `bfit` GUI.](inspect_tab.png){ width=80% }

While the GUI greatly facilitates rapid on-line analysis, the [`bfit`] API provides the flexibility needed for publishable analyses. The analysis tools and functions utilized in the GUI are readily accessible via the API, and documented in the [wiki]. Many of these tools are very general, easily interfacing with other Python packages, and can accommodate a great deal of complexity and sophistication.

# Acknowledgements

The author would like to thank the members of the $\beta$-NMR group at TRIUMF for their useful input and feedback. In particular, discussions with R. M. L. McFadden have been particularly useful. The author additionally acknowledges the support of a SBQMI QuEST fellowship.

# References

[`bfit`]: https://github.com/dfujim/bfit
[wiki]: https://github.com/dfujim/bfit/wiki
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is. 

```
Please paste terminal output within these quotes 

```

**To Reproduce**
Steps to reproduce the behaviour:
1. Which runs were loaded
2. Click on '....'
3. See error

**Desktop (please complete the following information):**
 - Machine or OS: [e.g. isdaq01 or OpenSUSE]
 - Version [e.g. 4.2.5]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
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
---
name: Documentation needed
about: Something is unclear in the wiki
title: ''
labels: documentation
assignees: ''

---

Please indicate if

**More content is needed**
- Please describe your issue and provide some context
- What information do you need?

**Content is unclear or inaccurate**
- Provide a link to the wiki page
- What is the issue with the page? 
- How could we make the page clearer or more accurate?
