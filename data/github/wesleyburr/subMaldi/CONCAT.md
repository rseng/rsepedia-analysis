---
title: "subMALDI: an open framework R package for processing irregularly-spaced mass spectrometry data"
tags:
  - mass spectrometry
  - spectral processing
  - MALDI-MS
authors: 
  - name: Kristen Yeh
    orcid: 0000-0002-3411-6816
    affiliation: 1
  - name: Sophie Castel
    orcid: 0000-0001-9086-0917
    affiliation: 4
  - name: Naomi L. Stock
    orcid: 0000-0002-3472-9284
    affiliation: 2
  - name: Theresa Stotesbury
    orcid: 0000-0001-6452-4389
    affiliation: 3
  - name: Wesley Burr
    orcid: 0000-0002-2058-1899
    affiliation: 4
affiliations:
  - name: Forensic Science Program, Trent University
    index: 1
  - name: Water Quality Center, Trent University
    index: 2
  - name: Faculty of Science, Forensic Science & Applied Bioscience, Ontario Tech University
    index: 3
  - name: Faculty of Science, Mathematics, Trent University
    index: 4
date: August 4, 2020
bibliography: paper.bib
output: pdf_document
---


# Summary

Mass spectrometry (MS) is an essential analytical technique used in many
fields of science, including chemistry, biology, medicine, and more
[@Gross:2011]. Its uses are varied, from biotechnology studies of
biomolecular sequencing [@Maux:2001], genetic analysis of human DNA
[@Null:2001], exploration of the structure of single cells [@Jones:2003]
and even examination of extraterrestrial objects [@Fenselau:2003]. 
This incredible breadth of applications using MS results in highly complex
data, which often requires significant processing in order to obtain
actionable insights.

Modern instrumentation often includes proprietary software for spectral
processing and analysis (e.g. Bruker Daltonics’ Data Analysis). These
tools, though convenient, often fail to provide sufficient documentation
of the algorithms employed in the software and have limited analytical
capabilities. Other commercial tools are available to supplement these
programs (e.g. Agilent Technologies’ MassHunter Profinder and Thermo
Scientific’s SIEVE^TM^), however, they come at a cost.  Open source
software for analysis of MS data is also available online. These
applications are often implemented in a variety of statistical computing
languages, including Python (e.g. pyOpenMS) [@Rost:2014], Matlab
(e.g. LIMPIC) [@Mantini:2007], C++ (e.g. ProteoWizard) [@Chambers:2012]
and R (e.g. MSnbase, MALDIquant) [@Gatto:2012; @Gibb:2012]. While more
accessible and well-documented than proprietary software, these available
open source applications [@Gibb:2016] often utilize complex data structures (e.g. S3 and
S4 class objects in R), which can make it difficult for researchers without
strong coding backgrounds to access their raw spectral data. In order to 
simplify the organization and processing of mass spectrometry data, we propose 
the R package `subMALDI`.

`subMALDI` is an open framework tool that permits organization,
pre-processing (smoothing, baseline correction, peak detection), and
normalization of spectral data sets without masking into S3 or S4 class
objects. After every step of processing, the *m/z* and intensity data
of each spectrum is readily accessible, providing researchers with a
more thorough understanding of the data manipulation that occurs during
analysis. As a result of the package's open framework, subMALDI data sets
are compatible with functions from a wide variety of other R packages,
and user-defined functions are easier to implement and test.

# Statement of Need

`subMALDI` permits the direct comparison of irregularly spaced
spectral replicates in an open framework, an important feature that
other open source tools do not contain. While matrix-assisted laser
desoprtion/ionization (MALDI) mass spectra (and, also, any single spectra
aquired data) are often visualized on a
continuous scale, the data observed are positive intensity values,
corresponding to discretely measured mass-to-charge (*m/z*) values
[@Stanford:2016]. When spectral replicates are acquired of a sample, there
is variation in the number and value of *m/z* responses with accompanying
peaks due to spectra centroiding in the mass analyzer. This results in 
irregularly spaced data. This has implications for the
statistical interpretations of inter-and intra-sample comparisons. In
order to generate meaningful results from unevenly spaced data, it
is essential that the data set be standardized by some means. In 
statistical computing languages, replicates often must be aligned against
the same data structure: for our purposes, this will be the default
data structure in R, the data.frame [@Wickham:2014].

`subMALDI` processes each raw spectrum with one of several smoothing
filters, baseline correction methods, and peak detection algorithms
included in the package. The processed spectral intensity values are then
aligned to an array of all the theoretically possible *m/z* values in
the observed mass range, at a specified resolution. The resulting data
frame contains all *m/z* data in the first column, with the intensity
data of each spectral replicate in adjacent columns.

`subMALDI` was designed for use by researchers who wish to organize,
process, and analyze single spectra data, particularly MS data, while still 
being able to access their raw data at various points throughout the process. 
It has been utilized in a scientific article in the *Journal of Forensic Chemistry*
[@Yeh:2020] and in our laboratory for analysis of MALDI-MS and electrospray-ionization
(ESI) MS data. The open framework format and data structures of `subMALDI` create a more
transparent pipeline for processing of MS data, where users can easily
access their raw data and better understand the processing algorithms
that are being executed on their data sets. The `subMALDI` framework is
intended to reduce the "black-box" characteristics of MS data analysis and
assist students and researchers in obtaining a more thorough understanding
of MS and the complex, diverse data sets that it is used to produce.

# Acknowledgement

We are grateful to the Canadian Foundation for Innovation, and the Ontario
Research Fund for funding the Bruker SolariX XR MALDI FT-ICR-MS in the
Water Quality Centre at Trent University.

# Funding

This work was supported by a Natural Sciences and Engineering Research
Council (NSERC) Discovery Grant to W. Burr (2017-04741) and a Vice
President Research Fund to T. Stotesbury (2019). Author K. Yeh was
supported by two NSERC Undergraduate Student Research Awards (USRA),
2019 and 2020.

# References
# subMaldi

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02694/status.svg)](https://doi.org/10.21105/joss.02694)

subMALDI is an open framework package for the R programming environment that permits organization,  pre-processing (smoothing, baseline correction, peak detection), and normalization of spectral data sets. As a result of the package's design, subMALDI data objects are compatible with functions from a wide variety of other R packages, and user-defined functions are easier to implement and test.

<center><img src="subMALDIprocessing.png" width="650"></center>

# Installing in R

To install this package in R, begin by downloading the source copy of the package from the [releases folder](https://github.com/wesleyburr/subMaldi/releases) and 
install it using type = "source". But, before you do this, a dependency from the BioConductor project is required. 
Begin by installing the BioConductor Manager:

    install.packages("BiocManager")
    
and then install the package MassSpecWavelet:

    BiocManager::install("MassSpecWavelet")
    
Then install the CRAN-supported packages required:

    install.packages(c("tidyr", "ggplot2", "ggpmisc", "reshape2", "RColorBrewer", "signal"))
    
Once the packages are all installed, you can then install the package:

    install.packages("subMALDI-1.0-2.tar.gz", type = "source")
   
Update the file name as necessary to whatever the copy you downloaded was set to for versioning.

Note: if you install from source, and set **dependencies = TRUE**, this should all be moot. This
assumes you run the commands exactly as given, in which case, the BioConductor sequence doesn't always
seem to work. For example,

    install.packages("subMALDI-1.0-2.tar.gz", type = "source", dependencies = TRUE)

# Getting Started

Refer to the [vignette on workflow](https://wesleyburr.github.io/subMaldi/articles/subMALDI_workflow.html)
for a brief summary of the package functionality, and the 
[vignette on processing](https://wesleyburr.github.io/subMaldi/articles/subMALDIprocessing.html) for a simple example.

# Contributing?

Examine our [contribution guidelines](Contributing.md) for more.

# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, caste, color, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

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
reported to the community leaders responsible for enforcement via [email](mailto:wesleyburr@trentu.ca).
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
[https://www.contributor-covenant.org/version/2/0/code_of_conduct.html][v2.0].

Community Impact Guidelines were inspired by 
[Mozilla's code of conduct enforcement ladder][Mozilla CoC].

For answers to common questions about this code of conduct, see the FAQ at
[https://www.contributor-covenant.org/faq][FAQ]. Translations are available 
at [https://www.contributor-covenant.org/translations][translations].

[homepage]: https://www.contributor-covenant.org
[v2.0]: https://www.contributor-covenant.org/version/2/0/code_of_conduct.html
[Mozilla CoC]: https://github.com/mozilla/diversity
[FAQ]: https://www.contributor-covenant.org/faq
[translations]: https://www.contributor-covenant.org/translations

# Contributing to subMALDI

<!-- This CONTRIBUTING.md is adapted from https://gist.github.com/peterdesmet/e90a1b0dc17af6c12daf6e8b2f044e7c -->

First of all, thanks for considering contributing to subMALDI! 👍 It's people like you that make it rewarding for us - the project maintainers - to work on subMALDI. 😊

subMALDI is an open source project, maintained by people who care. We are not directly funded to do so.

[repo]: https://github.com/wesleyburr/subMaldi
[issues]: https://github.com/wesleyburr/subMaldi/issues
[new_issue]: https://github.com/wesleyburr/subMaldi/issues/new
[email]: mailto:sophie.castel@ontariotechu.net
[website]: https://wesleyburr.github.io/subMaldi
[citation]: https://joss.theoj.org/papers/10.21105/joss.02694

## Code of conduct

Please note that this project is released using a [Contributor Code of Conduct](code_of_conduct.md), via [Contributor Convenant](https://www.contributor-covenant.org/). By participating in this project you agree to abide by its terms.

## How you can contribute

There are several ways you can contribute to this project. If you want to know more about why and how to contribute to open source projects like this one, see this [Open Source Guide](https://opensource.guide/how-to-contribute/).

### Share the love ❤️

Think subMALDI is useful? Let others discover it, by telling them in person, via Twitter or a blog post.

Using subMALDI for a paper you are writing? Consider [citing it][citation].

### Ask a question ⁉️

Using subMALDI and got stuck? Browse the [documentation][website] to see if you can find a solution. Still stuck? Post your question as an [issue on GitHub][new_issue]. While we cannot offer user support, we'll try to do our best to address it, as questions often lead to better documentation or the discovery of bugs.

Want to ask a question in private? Contact the package maintainer by [email][email].

### Propose an idea 💡

Have an idea for a new subMALDI feature? Take a look at the [documentation][website] and [issue list][issues] to see if it isn't included or suggested yet. If not, suggest your idea as an [issue on GitHub][new_issue]. While we can't promise to implement your idea, it helps to:

* Explain in detail how it would work.
* Keep the scope as narrow as possible.

See below if you want to contribute code for your idea as well.

### Report a bug 🐛

Using subMALDI and discovered a bug? That's annoying! Don't let others have the same experience and report it as an [issue on GitHub][new_issue] so we can fix it. A good bug report makes it easier for us to do so, so please include:

* Your operating system name and version (e.g., Mac OS 10.13.6).
* Your version of R, RStudio, and subMALDI
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug. A [reprex](https://rstudio.com/resources/webinars/help-me-help-you-creating-reproducible-examples/) would be lovely!

### Improve the documentation 📖

Noticed a typo on the website? Think a function could use a better example? Good documentation makes all the difference, so your help to improve it is very welcome!

#### The website

[This website][website] is generated with [`pkgdown`](http://pkgdown.r-lib.org/). That means we don't have to write any html: content is pulled together from documentation in the code, vignettes, [Markdown](https://guides.github.com/features/mastering-markdown/) files, the package `DESCRIPTION` and `_pkgdown.yml` settings. If you know your way around `pkgdown`, you can [propose a file change](https://help.github.com/articles/editing-files-in-another-user-s-repository/) to improve documentation. If not, [report an issue][new_issue] and we can point you in the right direction.

#### Function documentation

Functions are described as comments near their code and translated to documentation using [`roxygen2`](https://klutometis.github.io/roxygen/). If you want to improve a function description:

1. Go to `R/` directory in the [code repository][repo].
2. Look for the file with the name of the function.
3. [Propose a file change](https://help.github.com/articles/editing-files-in-another-user-s-repository/) to update the function documentation in the roxygen comments (starting with `#'`).

### Contribute code 📝

Care to fix bugs or implement new functionality for subMALDI? Awesome! 👏 Have a look at the [issue list][issues] and leave a comment on the things you want to work on. See also the development guidelines below.

## Development guidelines

We try to follow the [GitHub flow](https://guides.github.com/introduction/flow/) for development.

1. Fork [this repo][repo] and clone it to your computer. To learn more about this process, see [this guide](https://guides.github.com/activities/forking/).
2. If you have forked and cloned the project before and it has been a while since you worked on it, [pull changes from the original repo](https://help.github.com/articles/merging-an-upstream-repository-into-your-fork/) to your clone by using `git pull upstream master`.
3. Open the RStudio project file (`.Rproj`).
4. Make your changes:
    * Write your code.
    * Test your code (bonus points for adding unit tests).
    * Document your code (see function documentation above).
    * Check your code with `devtools::check()` and aim for 0 errors and warnings.
5. Commit and push your changes.
6. Submit a [pull request](https://guides.github.com/activities/forking/#making-a-pull-request).
---
name: Bug report
about: Create a report to help us improve
title: "[BUG] "
labels: ''
assignees: castels

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
 - R Version
 - Package Version [e.g. 22]

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
