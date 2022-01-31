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
---
title: "dev"
author: "Sophie Castel"
date: "1/27/2021"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

#########################################################################


## JOSS 11 (Closed)

```{r}
# loading uncompressed data files

load(file = "/home/sophie/subMaldi/data/After1.rda")
load(file = "/home/sophie/subMaldi/data/After2.rda")
load(file = "/home/sophie/subMaldi/data/Before1.rda")
load(file = "/home/sophie/subMaldi/data/Before2.rda")
load(file = "/home/sophie/subMaldi/data/Blank1.rda")
load(file = "/home/sophie/subMaldi/data/Blank2.rda")
load(file = "/home/sophie/subMaldi/data/bsline.rda")
load(file = "/home/sophie/subMaldi/data/Master.rda")
load(file = "/home/sophie/subMaldi/data/Master2.rda")

# test compress, comparing methods
save(Master, file = "/home/sophie/subMaldi/data/compressed/Master_T.rda", compress = TRUE)
save(Master, file = "/home/sophie/subMaldi/data/compressed/Master_xz.rda", compress = "xz")

file.info("/home/sophie/subMaldi/data/compressed/Master_T.rda")
file.info("/home/sophie/subMaldi/data/compressed/Master_xz.rda")

# opening test compress
load("/home/sophie/subMaldi/data/compressed/Master_xz.rda")

# saving else compressed files
save(After1, file = "/home/sophie/subMaldi/data/compressed/After1.rda", compress = "xz")
save(After2, file = "/home/sophie/subMaldi/data/compressed/After2.rda", compress = "xz")
save(Before1, file = "/home/sophie/subMaldi/data/compressed/Before1.rda", compress = "xz")
save(Before2, file = "/home/sophie/subMaldi/data/compressed/Before2.rda", compress = "xz")
save(Blank1, file = "/home/sophie/subMaldi/data/compressed/Blank1.rda", compress = "xz")
save(Blank2, file = "/home/sophie/subMaldi/data/compressed/Blank2.rda", compress = "xz")
save(bsline, file = "/home/sophie/subMaldi/data/compressed/bsline.rda", compress = "xz")
save(Master, file = "/home/sophie/subMaldi/data/compressed/Master.rda", compress = "xz")
save(Master2, file = "/home/sophie/subMaldi/data/compressed/Master2.rda", compress = "xz")

```


#########################################################################





## JOSS 5 (Closed)
Wrote biocViews: into DESCRIPTION install_github() can search BioConductor for MassSpecWavelet




#########################################################################




## JOSS 12 - generalize and fix redundancies in code (Closed)

## avgSpectra.R (done)

```{r}
avgSpectra <- function(dat, method = "mean", spectra_cols){
  
  # --------------
  # LOGICAL CHECKS
  # --------------
  
  if(length(spectra_cols) < 2){
    stop("Only one spectrum input. Please enter two spectra for averaging.")
  }
  
  if(!all(spectra_cols %in% colnames(dat))){
    logic <- which(!(spectra_cols %in% colnames(dat)))
    stop(c("Columns '",paste0(as.character(spectra_cols[logic]), sep = "', "), " not found in specified dataframe."))
  }
  
  if(method == "sum"){
    .avg_sum(dat = dat, spectra_cols) } 
  else if(method == "mean"){
    .avg_mean(dat = dat, spectra_cols) }
}

# --------------
# METHOD = SUM
# --------------

.avg_sum <- function(dat, spectra_cols){
  
  mz <- dat$full_mz
  
  spectra <- lapply(spectra_cols, function(x){dat[x]})
  i <- do.call(what = data.frame, args = c(spectra))
  
  dat <- cbind(mz, dat, Sum = apply(i, 1, sum, na.rm = TRUE))
  
  return(dat)
}

# --------------
# METHOD = MEAN
# --------------

.avg_mean <- function(dat, spectra_cols){
  
  mz <- dat$full_mz
  
  spectra <- lapply(spectra_cols, function(x){dat[x]})
  i <- do.call(what = data.frame, args = c(spectra))
  
  dat <- cbind(mz, dat, Average = apply(i, 1, mean, na.rm = TRUE))
  
  return(dat)
}

```

## find_max.R (done)

```{r}
# -----------------
# FIND PEAK MAXIMA
# -----------------

# Find peak maximum and associated m/z of spectra
find_max <- function (dat, mass_dat, spectra_cols){
  
  mz <- dat[[mass_dat]]
  
  spectra <- lapply(spectra_cols, function(x){dat[x]})
  i <- do.call(what = data.frame, args = c(spectra))
  rownames(i) <- mz
  
  max_i <- apply(i, 2, max)
  which_max <- apply(i, 2, which.max)
  max_mz <- mz[which_max]
  
  max_spec <- cbind(max_i, max_mz)
  colnames(max_spec) <- c("Intensity (Max)", "Mass")
  
  return(max_spec)
}

# -----------------------------------------------------------------------

# ------------------------
# FIND PEAK MAXIMA OF SET
# ------------------------

find_max_set <- function(dat, mass_dat, spectra_cols){

  # --------------
  # LOGICAL CHECKS
  # --------------
  
  stopifnot(
    is.character(mass_dat),
    is.character(spectra_cols),
    mass_dat %in% colnames(dat),
    all(spectra_cols %in% colnames(dat))
  )
  
  
  mz <- dat[[mass_dat]]
  
  spectra <- lapply(spectra_cols, function(x){dat[x]})
  i <- do.call(what = data.frame, args = c(spectra))
  rownames(i) <- mz
  
  which_max <- which(i == max(i), arr.ind = TRUE)
  
  max_spec <- data.frame(Intensity = as.numeric(i[which_max[1], which_max[2]]), Mass = as.numeric(rownames(which_max)))
  
  colnames(max_spec) <- c("Intensity (Max)", "Mass")
  rownames(max_spec) <- colnames(i)[which_max[2]]
  
  return(max_spec)
  
}
```

## norm_custimp.R (done)

```{r}
# --------------------------------------
# METHOD: CUSTOM M/Z VALUE, NOT PRECISE
# --------------------------------------

norm_custimp <- function(dat, mass_dat, norm_mz, spectra_cols, showHI = FALSE){
  
  # ---------------------
  # LOGICAL CHECKS
  # ---------------------
  
  stopifnot(
    is.character(mass_dat),
    is.character(spectra_cols),
    mass_dat %in% colnames(dat),
    all(spectra_cols %in% colnames(dat))
  )
  
  if(is.null(norm_mz)){
    stop('Please select a m/z value. norm_mz is NULL.')
  }   
  
  
  mz <- dat[[mass_dat]]
  
  logic_i <- which(startsWith(as.character(mz), as.character(norm_mz)))
  
  if(length(logic_i) == 0){
    stop(c("The selected value norm_mz = ", norm_mz, " not found in '", mass_dat,"'. Please try another value."))
  }
  
  spectra <- lapply(spectra_cols, function(x){dat[x]})
  i <- do.call(what = data.frame, args = c(spectra))
  
  i_cust <- i[logic_i,]
  
  if(nrow(i_cust) > 1){
    stop("More than one peak per spectra for the given norm_mz value. Please be more precise.")    
  }
  
  i_cust <- as.numeric(i_cust)
  
  logic <- i_cust == 0 
  
  if(any(logic)){
    stop(c("The selected maximum intensity is 0 in spectra labeled: ", paste0(colnames(i)[logic], sep = " ")))
  }
  
  i_norm <- t(t( i - min(i) )* ( i_cust - min(i) )^-1) # (i-min(i))/(i_cust-min(i)) Matrix multiplication
  
  if(!showHI){
    i_norm[i_norm > 1] <- 1
  }
  
  dat <- data.frame(mz, i_norm)
  colnames(dat) <- c("full_mz", spectra_cols)
  
  return(dat)
}

```

## norm_custom.R (done)

```{r}

# -------------------------
# METHOD: CUSTOM M/Z VALUE
# -------------------------

.norm_custom <- function(y, custom_y){
  return((y - min(y)) / (custom_y - min(y))) # single value
}

.normMethod_custom <- function(dat, mass_dat, norm_mz, spectra_cols, showHI = FALSE){
  
  # ---------------------
  # LOGICAL CHECKS
  # ---------------------
    
  stopifnot(
    is.character(mass_dat),
    is.character(spectra_cols),
    mass_dat %in% colnames(dat),
    all(spectra_cols %in% colnames(dat))
  )
    

    
  if(is.null(norm_mz)){
    stop('Please select a m/z value. norm_mz is NULL.')
  } 
  
    
  mz <- dat[[mass_dat]]
  
  spectra <- lapply(spectra_cols, function(x){dat[x]})
  i <- do.call(what = data.frame, args = c(spectra))
    
  i_cust <- as.numeric(i[which(mz == norm_mz),])
  
  if(all(is.na(i_cust))){
    stop(c("The selected value norm_mz = ", norm_mz, " not found in '", mass_dat,"'. Please try another value."))
  }
  
  logic <- i_cust == 0 
    
  if(any(logic)){
    stop(c("The selected maximum intensity is 0 in spectra labeled: ", paste0(colnames(i[logic]), sep = " ")))
  }
    
  i_norm <- t(t( i - min(i) )* ( i_cust - min(i) )^-1) # (i-min(i))/(i_cust-min(i)) Matrix multiplication
    
  if(!showHI){
    i_norm[i_norm > 1] <- 1
  }
    
  dat <- data.frame(mz, i_norm)
  colnames(dat) <- c("full_mz", spectra_cols)
    
  return(dat)
    
}

```

## norm_max.R (done)

```{r}
# -----------------------
# NORMALIZATION FUNCTION
# -----------------------

# Function rescales intensity to 0,1
.normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

# ------------------------------
# METHOD: MAX. OF EACH SPECTRUM
# ------------------------------

norm_max <- function(dat, mass_dat, spectra_cols){

  # ---------------------
  # LOGICAL CHECKS
  # ---------------------
  
  stopifnot(
    is.character(mass_dat),
    is.character(spectra_cols),
    mass_dat %in% colnames(dat),
    all(spectra_cols %in% colnames(dat))
  )
  
  
  
  mz <- dat[[mass_dat]]
  
  spectra <- lapply(spectra_cols, function(x){dat[x]})
  i <- do.call(what = data.frame, args = c(spectra))
  
  i_n <- apply(i, 2, FUN = .normalize)
  dat <- cbind(mz, i_n)
  colnames(dat) <- c("full_mz", spectra_cols)
  
  return(dat)
}
  
# ---------------------
# METHOD: MAX. OF SET
# ---------------------

# Normalization function specific to custom max_y
.normalize_set <- function(y, max_y){
  return((y - min(y)) / (max_y - min(y)))
}

.normMethod_max_set <- function(dat, mass_dat, spectra_cols){
  
  # ---------------------
  # LOGICAL CHECKS
  # ---------------------
  
  stopifnot(
    is.character(mass_dat),
    is.character(spectra_cols),
    mass_dat %in% colnames(dat),
    all(spectra_cols %in% colnames(dat))
  )

  if(length(spectra_cols) < 2){
    stop("Only one spectrum input. Please enter two spectra.")
  }
  
  
  mz <- dat[[mass_dat]]
  
  spectra <- lapply(spectra_cols, function(x){dat[x]})
  i <- do.call(what = data.frame, args = c(spectra))
  
  which_max <- which(i == max(i), arr.ind = TRUE)
  max_i <- i[which_max[1], which_max[2]]
  
  i_n <- apply(i, 2, FUN = function(x) {.normalize_set(y = x, max_y = max_i)} )
  dat <- cbind(mz, i_n)
  colnames(dat) <- c("full_mz", spectra_cols)
  
  return(dat)
}

```

## norm_median.R (done)

```{r}

# ---------------
# METHOD: MEDIAN
# ---------------

norm_med<- function(dat, mass_dat, spectra_cols){

  # ---------------------
  # LOGICAL CHECKS
  # ---------------------
  
  stopifnot(
    is.character(mass_dat),
    is.character(spectra_cols),
    mass_dat %in% colnames(dat),
    all(spectra_cols %in% colnames(dat))
  )

  if(length(spectra_cols) < 2){
    stop("Only one spectrum input. Please enter two spectra.")
  }
  
  
  mz <- dat[[mass_dat]]
  
  spectra <- lapply(spectra_cols, function(x){dat[x]})
  i <- do.call(what = data.frame, args = c(spectra))
  i[i == 0 ] <- NA
  
  i_med <- apply(i,2, FUN = median, na.rm = TRUE)
  max_med <- max(i_med)
  
  for(j in 1:length(spectra_cols)){
    
    # multiplicative factors
    factors <- c()
    
      if(i_med[j] == max_med){
        
        other <- i_med[-j]
        other_index <- which(i_med %in% other)
          
        factors <- rep(1,length(spectra_cols))
        
        factors[other_index] = i_med[j]/i_med[other_index]
        
      }
    }
    
  #i_mult <- i*factors
  
  i_mult <- t(t(i)*factors)
  
  dat <- cbind(mz, i_mult)
  colnames(dat) <- c("full_mz", spectra_cols)
  
  return(dat)
  
}
```

## norm_quantile.R (done)

```{r}
# --------------------------------------------------------------------------------------------
# Date: January 28, 2021
# Author: Kristen Yeh, Sophie Castel
# Title: Normalization Method - Quantile
# --------------------------------------------------------------------------------------------


# Quantile normalization consists of two steps:
    # 1. Create a mapping between ranks and values. For rank 1, find the n values,
       # one per spectrum, that are the smallest value on the spectrum, and save their
       # averages. Similarly for rank 2 and the second smallest values, and on up to
       # the n largest values, one per spectrum.
    # 2. For each spectrum, replace the actual values with these averages.

# Basically, index each intensity value in each spectrum then sort
# by intensity (value). Probably need 2 truncate. Average spectrum for the sorted spectra
# is calculated. The averages are then inserted into the sorted spectra.
# The spectra are then reverted to their unsorted order using the index.


norm_quantile <- function(dat, mass_dat, spectra_cols){
  
  # ---------------------
  # LOGICAL CHECKS
  # ---------------------
  
  stopifnot(
    is.character(mass_dat),
    is.character(spectra_cols),
    mass_dat %in% colnames(dat),
    all(spectra_cols %in% colnames(dat)),
    is.numeric(lower),
    is.numeric(upper)
  )    
  
  if(length(spectra_cols) < 2){
    stop("Only one spectrum input. Please enter two spectra for comparison.")
    }
 
  mz <- dat[[mass_dat]]
  
  spectra <- lapply(spectra_cols, function(x){dat[x]})
  i <- do.call(what = data.frame, args = c(spectra))
  
  i <- data.frame(mz, i)
  colnames(i) <- c("full_mz", spectra_cols)
  
  i_melt <- reshape2::melt(i, id.vars = "full_mz")
  colnames(i_melt) <- c("full_mz","Spectrum","Intensity")
  
  i_melt <- dplyr::filter(i_melt, Intensity != 0)
  
  i_melt <- i_melt %>% group_by(Spectrum)  %>% dplyr::mutate(id = row_number())
  ordr <- i_melt[order(i_melt$Intensity, decreasing = TRUE),]
  
  r <- rep(NA, length(spectra_cols))
  
  for(j in 1:length(spectra_cols)){
    r[j] <- max(as.numeric(i_melt$id[which(i_melt$Spectrum == spectra_cols[j])]))
  }
  
  r_min <- min(r)
  
  logic <-which(r == r_min)
  
  trunc <- ordr[(which(ordr$id <= r[logic])),]
  
  i_trunc <- data.frame(matrix(ncol = length(unique(trunc$Spectrum)), nrow = max(trunc$id)))
  
  for(j in 1:length(unique(trunc$Spectrum))){
    
    i_trunc[j] <- trunc$Intensity[which(trunc$Spectrum == spectra_cols[j])]
    colnames(i_trunc)[j] <- spectra_cols[j]
    
  }
  
  i_trunc <- data.frame(i_trunc, avg = apply(i_trunc, 1, FUN = mean, na.rm = TRUE))
  
  for(j in 1:length(unique(trunc$Spectrum))){
    trunc$Intensity[which(trunc$Spectrum == colnames(i_trunc)[j])] <- i_trunc$avg
  }
  
  dat <- trunc[order(trunc[mass_dat], decreasing = FALSE), ] %>% spread(Spectrum, Intensity)
  
  dat[is.na(dat)] <- 0
  dat <- dat[-2]
  return(as.data.frame(dat))
}

```

## norm_RMS.R (done)

```{r}

# -------------------------
# METHOD: ROOT MEAN SQUARE
# -------------------------

.RMS <- function(y, n){
  right <- 1 / (n - 1)
  left <- sum(y^2)
  prod <- prod(left, right)
  sqrt(prod)
}
  
.norm_RMS <- function(y, RMS){
  return(y/RMS) # single value
}

norm_RMS <- function(dat, mass_dat, spectra_cols){
  
  mz <- dat[[mass_dat]]
  
  spectra <- lapply(spectra_cols, function(x){dat[x]})
  i <- do.call(what = data.frame, args = c(spectra))
  
  n <- nrow(dat)
  
  rms <- apply(i, 2, FUN = function(x) { .RMS(y = x, n = n)})
  
  i_rms <- t(t(i)*rms^-1)  # i / rms matrix multiplication
  
  dat <- data.frame(mz, i_rms)
  
  colnames(dat) <- c("full_mz", spectra_cols)

  return(dat)
  
}

```

## norm_stdev.R (done)

```{r}
# ------------------------------------
# METHOD: STANDARD DEVIATION OF NOISE
# ------------------------------------

# Find a region of the spectrum with only noise and minimal peaks
# (Should be same region for all spectra)
# Evaluate standard deviation of intensity in that region for each spec
# Divide intensity of EACH PEAK IN SPEC by its st. dev.


norm_stdev <- function(dat, mass_dat, lower = 900, upper= 1100 , spectra_cols){
  
  # ---------------------
  # LOGICAL CHECKS
  # ---------------------
  
  stopifnot(
    is.character(mass_dat),
    is.character(spectra_cols),
    mass_dat %in% colnames(dat),
    all(spectra_cols %in% colnames(dat)),
    is.numeric(lower),
    is.numeric(upper)
  )
  
  mz <- dat[[mass_dat]]
  
  spectra <- lapply(spectra_cols, function(x){dat[x]})
  i <- do.call(what = data.frame, args = c(spectra))
  
  noise <- i[which(mz > lower & upper > mz),]
  noise[noise == 0] <- NA
  
  std_dev <- apply(noise, 2, FUN = sd, na.rm = TRUE)
  
  i_sd <- t(t(i)*std_dev^-1) # i / st_dev matrix multiplication
      
  dat <- data.frame(mz, i_sd)
  
  colnames(dat) <- c("full_mz", spectra_cols)
  
  return(dat)
}

```

## norm_TIC.R (done)

```{r}
# ----------------------
# METHOD: NORMALIZE TIC
# ----------------------

# This method evaluates the sum of all intensities of each spectrum in a set.
# If the sum is not equivalent between spectra, multiplies the spectra by
# normalization factor.


norm_TIC <- function(dat, mass_dat, spectra_cols){
  
  # ---------------------
  # LOGICAL CHECKS
  # ---------------------
  
  stopifnot(
    is.character(mass_dat),
    is.character(spectra_cols),
    mass_dat %in% colnames(dat),
    all(spectra_cols %in% colnames(dat)),
    is.numeric(lower),
    is.numeric(upper)
  )    
  
  if(length(spectra_cols) < 2){
    stop("Only one spectrum input. Please enter two spectra for comparison.")
    }
  
  
    mz <- dat[[mass_dat]]
    
    spectra <- lapply(spectra_cols, function(x){dat[x]})
    i <- do.call(what = data.frame, args = c(spectra))

    i_sum <- apply(i, 2, FUN = sum)
  
    max_sum <- max(i_sum)
    
    for(j in 1:length(spectra_cols)){
      
      if(i_sum[j] == max_sum){
      
        other <- i_sum[-j]
        other_index <- which(i_sum %in% other)
        
        factors <- rep(1,length(spectra_cols))
        
        factors[other_index] = i_sum[j]/i_sum[other_index]
        
      }
    }
    
    i_mult <- t(t(i)*factors)
    
    dat <- data.frame(mz, i_mult)
    colnames(dat) <- c("full_mz", spectra_cols)
    
    return(dat)
    
  }

```

## norm_rel_TIC (done)
```{r}
# ---------------------
# METHOD: RELATIVE TIC
# ---------------------

# This method evaluates the sum of all intensities of each spectrum in a set.
# If the sum is not equivalent between spectra, multiplies the spectra by
# normalization factor.
# Each peak in each spectrum is then divided by the normalized TIC.

# Divide intensity of each peak in spectrum by sum of all intensities
.norm_TIC <- function(y, all_y){
  return(y/all_y) # single value
}

norm_rel_TIC <- function(dat, mass_dat, spectra_cols){
  
  # ---------------------
  # LOGICAL CHECKS
  # ---------------------
  
  stopifnot(
    is.character(mass_dat),
    is.character(spectra_cols),
    mass_dat %in% colnames(dat),
    all(spectra_cols %in% colnames(dat)),
    is.numeric(lower),
    is.numeric(upper)
  )    
  
  
  mz <- dat[[mass_dat]]
    
  spectra <- lapply(spectra_cols, function(x){dat[x]})
  i <- do.call(what = data.frame, args = c(spectra))

  i_sum <- apply(i, 2, FUN = sum)
  
  max_sum <- max(i_sum)
  
  for(j in 1:length(spectra_cols)){
      
    if(i_sum[j] == max_sum){
      
      other <- i_sum[-j]
      other_index <- which(i_sum %in% other)
        
      factors <- rep(1,length(spectra_cols))
        
      factors[other_index] = i_sum[j]/i_sum[other_index]
        
      }
    }
  
    i_mult <- t(t(i)*factors)
    
    i_rel_TIC<- t(t(i_mult)*max_sum^-1) # i_mult / max_sum matrix multiplication
    
    dat <- data.frame(mz, i_rel_TIC)
    colnames(dat) <- c("full_mz", spectra_cols)
    
    return(dat)
}

```










## JOSS 6: roxygen2 (Closed)

```{r}
library(Rd2roxygen)

oxygnze <- function(file){
  
  options(roxygen.comment = "##' ")
  str(info <- parse_file(file))
  
  # parse_and_save() combines these two steps
  cat(create_roxygen(info), sep = "\n")
  
}

```



## avgSpectra.Rd (done)
```{r}
oxygnze(file = "/home/sophie/subMaldi/man/avgSpectra.Rd")
```

## baselineCorr.Rd (done)

```{r}
oxygnze(file = "/home/sophie/subMaldi/man/baselineCorr.Rd")
```


## createSpecDF.Rd (done)

```{r}
oxygnze(file = "/home/sophie/subMaldi/man/createSpecDF.Rd")
```


## find_max.Rd (done)

```{r}
oxygnze(file = "/home/sophie/subMaldi/man/find_max.Rd")
```


## find_max_set.Rd (done)

```{r}
oxygnze(file = "/home/sophie/subMaldi/man/find_max_set.Rd")
```


## mapSpectrum.Rd (done)

```{r}
oxygnze(file = "/home/sophie/subMaldi/man/mapSpectrum.Rd")
```


## normSpectra.Rd (done)

```{r}
oxygnze(file = "/home/sophie/subMaldi/man/normSpectra.Rd")
```


## peakDet.Rd (done)

```{r}
oxygnze(file = "/home/sophie/subMaldi/man/peakDet.Rd")
```


## readcsvDir.Rd (done)

```{r}
oxygnze(file = "/home/sophie/subMaldi/man/readcsvDir.Rd")
```


## readcsvSpec.Rd (done)

```{r}
oxygnze(file = "/home/sophie/subMaldi/man/readcsvSpec.Rd")
```


## rmveEmpty.Rd (done)

```{r}
oxygnze(file = "/home/sophie/subMaldi/man/rmveEmpty.Rd")
```


## smoothSpectrum.Rd (done)

```{r}
oxygnze(file = "/home/sophie/subMaldi/man/smoothSpectrum.Rd")
```


## subSpectra.Rd (done)
```{r}
oxygnze(file = "/home/sophie/subMaldi/man/subSpectra.Rd")
```




## JOSS 10: test code

## baselineCorr.R

```{r}

```

## baseMALDI.R

```{r}

```

## avgSpectra.R

```{r}

```

## baseline_loess.R

```{r}

```

## baseline_sub.R

```{r}

```

## createSpecDF.R

```{r}

```

## find_max.R

```{r}

```

## find_max_set.R

```{r}

```

## mapSpectrum.R

```{r}

```

## norm_RMS.R

```{r}

```

## norm_TIC.R

```{r}

```

## norm_rel_TIC.R

```{r}

```

## norm_custimp.R

```{r}

```

## norm_custom.R

```{r}

```

## norm_max.R

```{r}

```

## norm_median.R

```{r}

```

## norm_quantile.R

```{r}

```

## norm_stdev.R

```{r}

```

## normSpectra.R

```{r}

library(testthat)

test_that("Method = NULL yields error", {
  expect_error(normSpectra(dat = Master, 
                           mass_dat = "full_mz",
                           method = "apple",
                           norm_mz = NULL,
                           upper = NULL,
                           lower = NULL,
                           spectra_cols = NULL,
                           showHI = FALSE))
})



```

## peakDet.R

```{r}

```

## peak_sub.R

```{r}

```

## plotSpectra.R

```{r}
library(testthat)





```

## readcsvDir.R

```{r}

```

## readcsvSpec.R

```{r}

```

## rmveEmpty.R

```{r}

```

## smoothSpectrum.R

```{r}

```

## smooth_sub.R

```{r}

```

## subSpectra.R

```{r}

```


---
title: "subMALDI: Sample Workflow"
author: Kristen Yeh^1^, Sophie Castel^2^ and Wesley Burr ^2^
output: 
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 1
    df_print: kable
    number_sections: true
vignette: >
  %\VignetteIndexEntry{subMALDI: Sample Workflow}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

^1^ Forensic Science Program, Trent University, Peterborough, ON, Canada

^2^ Faculty of Science, Mathematics, Trent University, Peterborough, ON, Canada

# Introduction

In this vignette we demonstrate a processing workflow for irregularly-spaced mass spectrometry data using __*subMALDI*__. This package is freely available from [GitHub](https://github.com/wesleyburr/subMaldi) and was created using MALDI FT-ICR mass spectrometry data. 

---


# Installation and Import

## Installation

To download the package, select the '.tar.gz' file from the git repository. This file can then used to install the package in R. Below, we load the package and its dependencies.

```{r setup}
library(subMALDI)
```

## Import
Two import functions are included in the package. Both read .CSV files which contain spectral data from a single spectrum per file. While `readcsvSpec()` imports one spectrum at a time and turns them into a pairwise data frame of *m/z* and intensity data, `readcsvDir()` is capable of reading all the .CSV files in a folder. The pairwise data frames generated for all of these spectra are then output as .rda files into a designated output folder.

For example, let's say we have a set of 5 spectral replicates in our working directory. These are stored as .CSV files in a folder called "my_spec". In the original file, mass data is stored in a column called "m.z", while intensity data is stored in a column called "I". If we wanted to import only one of the spectra, specifically "Sample1.csv", then we would run the following code. 

```{r eval = FALSE}
my_spec1 <- readcsvSpec(spec_file = "~/my_spec/Sample1.csv", massCol = "m.z", intenseCol = "I")
```

Conversely, if we wanted to import all spectra in the folder "my_spec" into a new folder called "out", we would run the following.

```{r eval = FALSE}
readcsvDir(direct = "~/my_spec/", massCol = "m.z", intenseCol = "I", output = "~/my_spec/out/")

load("~/my_spec/out/Sample1.rda")
```

This will create a pairwise data frame of *m/z* and intensity data, much like the one shown below.

```{r }
data("Blank1")
head(Blank1)
```


---

## Reading in an ASCII File
The mass spectrometer may produce raw output data in \code{.ascii} format. In this case, the \code{read_ascii()} function can be used to read in the data file and convert it to a \code{data.frame} object.  Simply specify the filepath, as shown below:

```{r}
# spec_data <- read_ascii(filename = "~/filename.ascii")
```

Once loaded into the global environment, the object can be saved to the local directory as a \code{.csv} file using the \code{save()} function:

```{r}
# write.csv(spec_data, file = "~/filename.csv")
```

# Creating a Mapped Data Frame

In order to compare unevenly spaced spectral data, it is essential that the spectra be standardized by some means. When performing R analysis, replicates must be aligned against the same data structure, requiring that each spectrum has a consistent number of data points. 

## Create an Empty Data Frame

__*subMALDI*__ permits the comparison of unevenly spaced mass spectrometry data by mapping each peak in a spectrum to all of the theoretically possible *m/z* values in the analyzed mass range.

First, an empty data frame containing all of the theoretically possible *m/z* values is created. As the sample data sets included with this package were acquired between the mass range of *m/z*  53.76 to 1100 with a resolution of four decimal places, the empty data frame we create will be customized for this data. This is performed using the following code.

```{r }
spec_df <- createSpecDF(min_mz = 53.76, max_mz = 1100, res = 1e-04, dig = 4)

head(spec_df)
```

We have created a long data frame containing *m/z* data in the first column. This data frame can now be manipulated with `dplyr`'s `select()` function to create columns specific to our spectral data. 

We will add 4 columns of sample data, corresponding to 2 matrix blanks and 2 sample spectra.

```{r }
spec_df <- dplyr::select(spec_df, "full_mz")
spec_df <- transform(spec_df, "Blank1" = 0, "Blank2" = 0, "Before1" = 0, "Before2" = 0)

head(spec_df)
```

Once the data frame is customized for our sample data, each spectrum can now be mapped to its respective column.

## Map Spectral Data to Empty Data Frame

Before any spectra are mapped to our data frame, the spectral data must be loaded. The pairwise spectra included with the package are loaded below.

```{r }
data("Blank1")
data("Blank2")
data("Before1")
data("Before2")
```

Now that the pairwise spectra have been loaded, we can map them to our empty data frame. It is important to note that the `thresh` and `dig` values must match the `res` and `dig` values used in `createSpecDF()`, otherwise data will fail to map.

```{r }
spec_df <- mapSpectrum(dat = Blank1, massCol = "mass", intenseCol = "Intensity", 
            spec_df = spec_df, colName = "Blank1", thresh = 1e-04, dig = 4)
spec_df <- mapSpectrum(dat = Blank2, massCol = "mass", intenseCol = "Intensity", 
            spec_df = spec_df, colName = "Blank2", thresh = 1e-04, dig = 4)
spec_df <- mapSpectrum(dat = Before1, massCol = "mass", intenseCol = "Intensity", 
            spec_df = spec_df, colName = "Before1", thresh = 1e-04, dig = 4)
spec_df <- mapSpectrum(dat = Before2, massCol = "mass", intenseCol = "Intensity", 
            spec_df = spec_df, colName = "Before2", thresh = 1e-04, dig = 4)

head(spec_df)
```

## Remove Empty Rows

From the output above, we can see there are several rows which contain only zero intensity values. The data frame is `r nrow(spec_df)` rows long, and should be reduced in size to minimize the computational load of further processing functions.

It is important to note that *ALL* samples or spectral replicates to be compared must be mapped to the data frame prior to removal of empty rows, as new samples cannot be mapped after the original mapping *m/z* vector has been truncated.

Empty rows are removed from our mapped spectral data frame below. 

```{r }
spec <- rmveEmpty(spec_df)

head(spec)
nrow(spec)
```

The number of rows has been reduced by `r nrow(spec_df)-nrow(spec)`. We can now proceed with further processing.


---


# Manipulate Raw Data

In our spectral data frame we have three sets of spectral replicates: 2 matrix blanks, 2 sample spectra acquired prior to chemical intervention, and 2 sample spectra acquired after chemical intervention.

To distinguish between peaks originating from the matrix and peaks originating from the sample, we will average our matrix blank spectra, and then subtract the average blank spectrum from each sample.

## Average Intensity Data

There are two methods available for averaging our intensity data: sum or mean. In this example we will use the method "mean" on the matrix blanks "Blank1" and "Blank2". This is performed using the code below.

```{r }
avg_spec <- avgSpectra(spec, method = "mean", spectra_cols = c("Blank1", "Blank2"))

head(avg_spec)
```

A new column has now been added to our original data frame, containing the averaged blank intensity data. Now that the data has been averaged, we can subtract our blank spectrum from each sample.

## Subtract Intensity Data

We will subtract our averaged matrix blank from samples "Before1" and "Before2" below. First, we will restructure our data frame to remove the matrix blank columns, whose intensity data is now contained in the "Average" column. Next, we will add new columns to the data frame for the new subtracted intensity data.

```{r}
spec <- dplyr::select(avg_spec, "full_mz", "Before1", "Before2", "Average")
sub_spec <- transform(spec, "Sub_Before1" = 0, "Sub_Before2" = 0)
```

Finally, we will subtract our averaged matrix blank from each sample spectrum.

```{r }
sub_spec <- subSpectra(dat = sub_spec, Blank_Var = "Average", 
                       Sample = "Before1", Sub_Sample = "Sub_Before1")
sub_spec <- subSpectra(dat = sub_spec, Blank_Var = "Average", 
                       Sample = "Before2", Sub_Sample = "Sub_Before2")

head(sub_spec)
```

---


# Normalization Methods

In addition to simple subtraction and averaging functions, __*subMALDI*__ also offers several normalization methods. These are listed below:

1. Custom *m/z*
1. Maximum Intensity
1. Root Mean Square (RMS)
1. Total Ion Current (TIC)
1. Median
1. Standard Deviation of Noise
1. Quantile

A sample normalization for each method will be shown in order, below. We will begin by loading a smaller sample data set to work with, and subsetting it even further.

```{r }
data("Master2")
norm_spec <- dplyr::select(Master2, "full_mz", "Before1", "Before2")
head(norm_spec)
```

```{r}
norm_spec <- subset(norm_spec, full_mz >= 53.76 & full_mz <= 800)
```


The raw data is plotted below.

```{r warning = FALSE, fig.align = "center", fig.width = 6, fig.height = 3}
plotSpectra(norm_spec, mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), x_ticks = 5)
```

Please note that while only two spectra are normalized in this vignette, __*subMALDI*__ is capable of normalizing several spectral replicates at once, the resulting plots of which can be arranged using the __*nrows*__ argument.

## Custom *m/z* Normalization

Custom *m/z* normalization is one of the simpler methods included with __*subMALDI*__. This method allows the user to select a specific *m/z* peak and normalizes the intensity of the spectrum to the intensity of the selected peak. In other words, the maximum intensity in the spectrum is set to the intensity of the selected peak.

```{r warning = FALSE, fig.show = "hold"}
custom_ex <- normSpectra(norm_spec, mass_dat = "full_mz", method = "custom", 
                         norm_mz = 255.23, spectra_cols = c("Before1", "Before2"))

plotSpectra(norm_spec, mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), 
            x_ticks = 5)
plotSpectra(custom_ex, mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), 
            x_ticks = 5, palette = "Accent")
```

The original spectra are shown to the right, while the normalized spectra are displayed on the left.

## Maximum Intensity Normalization

Maximum intensity normalization is another simple method included with __*subMALDI*__. This method takes the most intense peak in a spectrum and re-scales the intensity data from an exponential scale (ex. 0 to 1e11) to a scale of 0,1. 

A sample normalization is performed below.

```{r warning = FALSE, fig.show = "hold"}
max_ex <- normSpectra(norm_spec, mass_dat = "full_mz", method = "max", 
                      spectra_cols = c("Before1", "Before2"))

plotSpectra(norm_spec,mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), x_ticks = 5)
plotSpectra(max_ex, mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), x_ticks = 5, palette = "Accent")
```

The original spectra are shown to the right, while the normalized spectra are displayed on the left.

### Maximum Intensity of Set Normalization

__*subMALDI*__ also offers a maximum intensity normalization for a set of spectral replicates. In this method, the intensity of the most intense peak in the most intense spectrum in the data set is set as the maximum intensity. Again, the data is re-scaled to 0,1. 

In this case, 1 represents the same value for all spectra in a data set, whereas in the previous method, 1 represents different intensity values in different spectra.

```{r warning = FALSE, fig.show = "hold"}
maxset_ex <- normSpectra(norm_spec, mass_dat = "full_mz", method = "max_set", 
                         spectra_cols = c("Before1", "Before2"))

plotSpectra(norm_spec,mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), x_ticks = 5)
plotSpectra(maxset_ex, mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), 
            x_ticks = 5, palette = "Accent")
```

The original spectra are shown to the right, while the normalized spectra are displayed on the left.

## Root Mean Squared (RMS) Normalization

In this normalization method, the intensity of each peak in a spectrum is divided by the spectrum's root mean squared error. The equation for RMS error is given below, where *n* represents the number of peaks in each spectrum and *I* represents the intensity of each individual peak in a spectrum.

$RMS = \sqrt{\frac{1}{n-1}*\sum{I^2}}$

A sample RMS normalization is performed below.

```{r warning = FALSE, fig.show = "hold"}
RMS_ex <- normSpectra(norm_spec, mass_dat = "full_mz", method = "RMS", spectra_cols = c("Before1", "Before2"))

plotSpectra(norm_spec,mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), x_ticks = 5)
plotSpectra(RMS_ex, mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), 
            x_ticks = 5, palette = "Accent")
```

The original spectra are shown to the right, while the normalized spectra are displayed on the left.

## Total Ion Current (TIC) Normalization

Total ion current (TIC) is one of the most common normalization techniques for MALDI-MS data^1^. 

In this method, the total intensity of each spectrum in a data set is evaluated. If the total intensity of each spectrum in the data set is not already equal, the method multiplies each spectrum's intensities by a normalization factor to equalize the TIC.

A sample normalization is performed below.

```{r warning = FALSE, fig.show = "hold"}
TIC_ex <- normSpectra(norm_spec, mass_dat= "full_mz", method = "TIC", spectra_cols = c("Before1", "Before2"))

plotSpectra(norm_spec,mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), x_ticks = 5)
plotSpectra(TIC_ex, mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"),
            x_ticks = 5, palette = "Accent")
```

The original spectra are shown to the right, while the normalized spectra are displayed on the left.

### Relative TIC Normalization

__*subMALDI*__ also includes a TIC normalization method which re-scales the intensity data to a relative intensity scale. This technique is identical to the method mentioned above, however, it divides each peak in a data set by the normalized TIC. 

A sample normalization is performed below.

```{r warning = FALSE, fig.show = "hold"}
rel_TIC_ex <- normSpectra(norm_spec, mass_dat = "full_mz", method = "rel_TIC", spectra_cols = c("Before1", "Before2"))

plotSpectra(norm_spec, mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), x_ticks = 5)
plotSpectra(rel_TIC_ex, mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"),
            x_ticks = 5, palette = "Accent")
```

The original spectra are shown to the right, while the normalized spectra are displayed on the left. 

The TIC of each spectrum in the output of the relative TIC normalization are equal to one. The TIC of the data set should be equal to the number of spectra in the set. An example is shown below.

In our output, we have two spectra. 

```{r}
head(rel_TIC_ex)
```

The sum of each spectrum's respective intensities is evaluated below, in addition to the total sum of intensities in the output data set. We can see that both individual spectra have normalized TICs of 1, and that the total intensity of the spectrum is 2, the number of spectra in the set.

```{r}
sum(rel_TIC_ex$Before1)
sum(rel_TIC_ex$Before2)
sum(rel_TIC_ex[,-1])
```

## Median Normalization

Median normalization is similar to the TIC normalization method^1^. In both techniques, a simple operation is performed on the intensity data of each spectrum in a data set. If the result of the operation is not equal between all spectra, then the intensities are multiplied by a correction factor to make them so^1^.

Unlike TIC normalization, which evaluates the sum of all intensities in a spectrum, median normalization evaluates the median of each spectrum. 

First, the intensity data from each spectrum is sorted from most to least intense. Next, the data set is truncated so that all spectra contain the same number of data points^1^. The number of peaks removed is dependent on the variation in size between each spectrum in the data set. The spectrum containing the least amount of non-zero peaks acts as the guide for truncation.

Once the spectra have been corrected to contain the same number of data points, the median intensity value in the sorted data set is evaluated. If the median is not consistent between spectral replicates, then the intensity data of each spectrum is multiplied by a normalization factor.

A sample median normalization is performed below.

```{r warning = FALSE, fig.show = "hold"}
median_ex <- normSpectra(norm_spec, mass_dat = "full_mz", method = "median", spectra_cols = c("Before1", "Before2"))

plotSpectra(norm_spec,mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), x_ticks = 5)
plotSpectra(median_ex, mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), 
            x_ticks = 5, palette = "Accent")
```

The original spectra are shown to the right, while the normalized spectra are displayed on the left.

## Standard Deviation of Noise Normalization

Standard deviation of noise normalization is focused on normalizing noisy regions of each spectrum. The idea is that the standard deviation of noisy regions should be equal in all spectra^1^. 

In order to use this method, a large region of the spectra must be identified wherein there is *only* noise and no peaks. This region must be the same between all spectral replicates to be normalized.

A plot of our raw data is shown below. Here we can see that there is a noisy region from *m/z* 900 to 1000. We will use this region for our sample normalization.

```{r warning = FALSE, fig.width = 6, fig.height = 3, fig.align = "center"}
plotSpectra(norm_spec, mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), 
            x_ticks = 5, min_mz = 750, max_mz = 1100)
```

Once a spectral region of noise has been identified, the standard deviation of the intensities in this region is evaluated. The intensity of each peak in a spectrum is then divided by that spectrum's standard deviation of noise. The output spectra should have a standard deviation of 1 in the selected noise region^1^.

A sample standard deviation of noise normalization is performed below, using a lower *m/z* bound of 900 and an upper bound of *m/z* 1000.

```{r warning = FALSE, fig.show = "hold"}
stdev_ex <- normSpectra(norm_spec, mass_dat = "full_mz", method = "stdev", 
                        lower = 900, upper = 1000, spectra_cols = c("Before1", "Before2"))

plotSpectra(norm_spec,mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), x_ticks = 5)
plotSpectra(stdev_ex, mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), x_ticks = 5, palette = "Accent")
```

The original spectra are shown to the right, while the normalized spectra are displayed on the left.

## Quantile Normalization

Quantile normalization is the final method included in the package. The method was originally developed for multiple high-density arrays^1^.

Quantile normalization consists of two steps^2^:

1. An index is created for all intensity values in a spectrum to document their original order. The data is then sorted by decreasing intensity. For each rank, the average of the intensities (one per spectrum) is saved.
2. The actual intensity values in each spectrum are replaced with these averages. The data is then resorted by the index.


This makes the distributions of the values equal in all spectra^1^. The quantile normalization technique is demonstrated below.

```{r warning = FALSE, fig.show = "hold"}
quantile_ex <- normSpectra(norm_spec, mass_dat = "full_mz", method = "quantile", spectra_cols = c("Before1", "Before2"))

plotSpectra(norm_spec,mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), x_ticks = 5)
plotSpectra(quantile_ex, mass_dat = "full_mz", spectra_cols = c("Before1", "Before2"), x_ticks = 5, palette = "Accent")
```

The original spectra are shown to the right, while the normalized spectra are displayed on the left.

---

# Plotting Spectra

Several sample plots were shown above to illustrate the different normalization methods. In this section we will show the features of the \code{plotSpectra()} function, and demonstrate its capabilities.

We can plot data from the original pairwise data frame, or from a mapped spectral data frame. An example of the former is shown below. 

In this example, we use the arguments `lbls = TRUE` and `span = 9`. This labels the peak maxima in the spectrum. The span argument determines how many maxima will be labeled. In this case, we've selected a value of 9. This means that four peaks on either side of each labeled peak are not labelled. 

```{r warning = FALSE, fig.width = 6, fig.height = 3, fig.align = "center"}
data("Blank1")
plotSpectra(dat = Blank1, mass_dat = "mass", spectra_cols = "Intensity", 
             lbls = TRUE, span = 9, x_ticks = 6)
```

We can also customize the color and the *m/z* range of the plot. In this example we will plot a single spectrum from a mapped spectral data frame without labels. We will change the color of the spectrum to `"black"` and zoom in on the *m/z* region from *m/z* 100 to 500.

```{r warning = FALSE, fig.width = 6, fig.height = 3, fig.align = "center"}
data("Master")
plotSpectra(dat = Master, mass_dat = "full_mz", spectra_cols = "After2", 
             colours = "black", min_mz = 100, max_mz = 500, x_ticks = 7)
```

Not all spectral data has the same resolution. As such, the labels on each mass spectrum should accurately reflect the resolution of the raw data.

By default, __*subMALDI*__ displays all labels in the format XXX.XXXX, as this is the resolution of most of the data included with the package. In this example, we will customize the resolution of our *m/z* labels to match a lower resolution data set.

Here we load the low resolution spectral data set.

```{r}
data("Master2")
```

We will customize our label output using the argument \code{lbl.fmt}. The default value for this argument is "%3.4f". This means our output will have three digits on the left of the decimal point and four to the right. 

Our low resolution data only has two decimal places, therefore our \code{lbl.fmt} argument will be "%3.2f".

```{r  warning = FALSE, fig.width = 6, fig.height = 3, fig.align = "center"}
plotSpectra(dat = Master2, mass_dat = "full_mz", spectra_cols = "Before1", lbls = TRUE, lbl.fmt = "%3.2f", x_ticks = 7)
```

## Plotting Multiple Spectra

When plotting multiple spectra on the same *m/z* axis, the `nrows` argument should be used. Here, we set `nrows = 1` such that the four plots are stacked on top of one another, and sharing the same *m/z* axis. 

```{r warning = FALSE, fig.width = 6, fig.height = 6, fig.align = "center"}
plotSpectra(dat = Master, mass_dat = "full_mz", spectra_cols = c("Before1", "Before2", "After1", "After2"), 
            lbls = TRUE, span = 15, nrows = 4, x_ticks = 7)
```

If a grid format is preferred, the `nrows` argument can be adjusted and will produce the best-fitting layout.

---

# References
1. O. Haglund, Qualitative comparison of normalization approaches in MALDI-MS, M.Sc Thesis, Royal Institute of Technology, Stockholm, Sweden, 2008.
2. K. V. Ballman, D. E. Grill, A. L. Oberg, T. M. Therneau, Faster Cyclic Loess: normalizing RNA analysis via linear models, *Bioinformatics*, 2004, **20**, 2778-2786. 







---
title: "Processing with subMALDI"
author: Kristen Yeh^1^, Sophie Castel^2^ and Wesley Burr ^2^
output: 
  rmarkdown::html_vignette:
    toc: true
    toc_depth: 1
    number_sections: true
vignette: >
  %\VignetteIndexEntry{Processing with subMALDI}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```


^1^ Forensic Science Program, Trent University, Peterborough, ON, Canada

^2^ Faculty of Science, Mathematics, Trent University, Peterborough, ON, Canada

# Introduction

In this vignette we demonstrate a pre-processing workflow for smoothing, baseline correction, and peak detection of raw, irregularly-spaced mass spectrometry data using `subMALDI`. This package is freely available from [GitHub](https://github.com/wesleyburr/subMaldi) and was created using MALDI FT-ICR mass spectrometry data. 

It is recommended that users first read through the *"`subMALDI` Workflow"* vignette before proceeding with pre-processing functions. This will assist users in understanding the basic functions included in `subMALDI`, including import, plotting, and mapping functions for comparison of multiple spectra.

## Loading Sample Data

We'll begin this vignette by installing the package and loading a sample data set included with the package, called `"bsline.rda"`.


```{r, warning = FALSE, error = FALSE, message = FALSE}
library("subMALDI")
data("bsline")
```

Although most modern-day high resolution mass spectrometry instruments come with their own spectral processing software which are capable of automatic baseline correction, this software can fail, producing a spectrum with an uneven baseline. 

Our data set `bsline` is an example of this scenario. A plot of the raw spectrum is shown below.

```{r warning = FALSE, fig.width = 6, fig.height = 3, fig.align = "center"}
plotSpectra(dat = bsline, mass_dat = "mass", spectra_cols = "raw",
             colours = "black", min_mz = 200, max_mz = 5000, x_ticks = 7)
```


Observing this plot, we can clearly see the baseline of the spectrum is increasing as the x-axis increases, with the exception of a nearly Gaussian-shaped curve from the *m/z* range 700 to 1700. Before we can extract meaningful data from this raw spectrum, this irregular baseline must be corrected.

We will work through the pre-processing and baseline correction below, demonstrating several different methods for each processing step.

---

# Smoothing Filters

As the sensitivity of high-resolution mass spectrometry instruments increases, so does the amount of detectable noise. This produces spectra which appear thicker around the baseline than at high intensity peaks, without proper spectral processing. We can see this phenomenon in our sample data set `bsline` from the *m/z* range 200 to 1700. Here, the baseline of the spectrum is very thick relative to that at higher *m/z* values and high intensity peaks. 

In order to minimize the noise in our spectrum and create a more uniform signal throughout the analyzed *m/z* range, we will use the `subMALDI` function `smoothSpectrum`. This function offers two different filters for spectral smoothing: the Savitzky-Golay and moving average filters.

## Savitzky-Golay

The Savitzy-Golay filter is a least-squares smoothing method [1]. It performs a least squares fit of a small set of consecutive data points to a polynomial and takes the central point of the fitted polynomial curve as the output [1,2]. 

The output of a Savitzky-Golay smoothing filter is given by the following equation  [1,2]: 
$y[n] = \frac{\sum{A*x[n-i]}}{\sum{A}}$

Where `n` denotes signal index and `k` indicates filter width [1]. `A` controls the polynomial order [1,2].

Use of `subMALDI`'s Savitzy-Golay filter is demonstrated below. When using `smoothSpectrum(method = "sgolay"`, the following arguments are required:

* `p` = Filter order
* `n` = Filter length, must be odd
* `m` = Return the m-th derivative of the filter coefficients (Default = 0)
* `ts` = Time scaling factor (Default = 1)


```{r}
sgolay <- smoothSpectrum(dat = bsline, mass_dat = "mass", intensity_dat = "raw",
                       method = "sgolay", p = 4, n = 7, m = 0, ts = 1)
```


We'll now compare our smoothed version of `bsline` to the raw spectrum by merging the intensity data into a single data frame.


```{r warning = FALSE, fig.width = 6, fig.height = 4, fig.align = "center"}
test <- data.frame(bsline$mass, bsline$raw, sgolay$sg)
names(test) <- c("mass", "Raw Spectrum", "Savitzky-Golay Smooth")

plotSpectra(dat = test, mass_dat = "mass", spectra_cols = c("Raw Spectrum", "Savitzky-Golay Smooth"), 
            min_mz = 200, max_mz = 5000, x_ticks = 7)
```


The noise of the raw spectrum has now been reduced.

## Moving Average

The output of a moving average smoothing filter is denoted by the following equation [1]: $y[n] = \frac{1}{2k+1}*\sum{x[n-i]}$

Again, `k` denotes filter width and `n` indicates signal index. The larger the filter width, the more intense the smoothing filter [1].

We will use the `"mov_avg"` smoothing filter on our sample data set below. When using `smoothSpectrum(method = "sgolay"`, the `n` argument, window size/filter width, is required. Since `bsline` is only extraneously noisy in the first portion of the spectrum, the value of `n` should not be too large (< 10), as this will risk distorting the shape of peaks at higher *m/z* values. We will use `n = 7` in this example.


```{r}
mov_avg <- smoothSpectrum(dat = bsline, mass_dat = "mass", intensity_dat = "raw",
                          method = "mov_avg", n = 7)
```


We'll now compare our smoothed version of `bsline` to the raw spectrum by mapping the raw and smoothed intensities to a spectral data frame called `test`. The code that follows is explained in depth in the vignette *"`subMALDI` Workflow".*

Once our raw and smoothed spectra have been mapped to the data frame, we can remove *m/z* rows that do not correspond to any intensities in either spectra with `rmveEmpty()`.


```{r}
test <- createSpecDF(min_mz = 200, max_mz = 5000)
test <- dplyr::select(test, full_mz)
test <- transform(test, "RawSpectrum" = 0, "MovingAverage" = 0)
test <- mapSpectrum(bsline, massCol = "mass", intenseCol = "raw", 
                 colName = "RawSpectrum", spec_df = test)
test <- mapSpectrum(mov_avg, massCol = "mass", intenseCol = "mov_avg", 
                 colName = "MovingAverage", spec_df = test)
test <- rmveEmpty(test)
```


The raw and smoothed spectra are plotted below.


```{r warning = FALSE, fig.width = 6, fig.height = 4, fig.align = "center"}
plotSpectra(dat = test, mass_dat = "full_mz", spectra_cols = c("RawSpectrum", "MovingAverage"), min_mz = 200, max_mz = 5000,
            x_ticks = 5, intensity_scale = "fixed")
```


The noise in the smoothed spectrum is now noticeably reduced relative to the raw spectrum.

---

# Baseline Correction Tools

Baseline correction often consists of two steps. First, estimating the baseline of the spectrum, and second, subtracting the baseline from the signal [1]. `subMALDI` offers three methods for baseline estimation: local monotone minimum, linear interpolation, and LOESS curve fitting.

Depending on the level of noise in your data, you may choose to skip spectral smoothing and proceed to baseline correction and peak detection. `subMALDI`'s baseline correction functions are compatible with both raw and smoothed spectra. 

In the examples that follow, we will correct the baseline of our smoothed spectra from section 2, `sgolay` and `mov_avg`, and the baseline of our raw spectrum, `bsline`, using `baselineCorr()`.

## Monotone Minimum

The monotone minimum method computes the difference in intensity between adjacent peaks to determine the slope of each point [1]. The baseline is then estimated as follows, moving from the leftmost to rightmost point:

* If slope of point *A* is greater than 0, the nearest point *B* with slope of less than 0 is located. Points in between *A* and *B* are recorded as potential peaks. 
* If slope of *A* < 0, nearest point *B* with slope > 0 is located. Points between *A* and *B* are recorded as valleys. 

Each intensity value in the spectrum is then subtracted by the intensity of the nearest valley, correcting any irregularities in the baseline of the data.

Use of monotone minimum method of baseline correction is demonstrated below, on the Savitzky-Golay smoothed spectrum `sgolay` from section 2.1.


```{r}
mono_min <- baselineCorr(dat = sgolay, mass_dat = "mass", 
                        intensity_dat = "sg", method = "monotone_min")
```


The smoothed and baseline corrected intensities are merged into a single data frame below, so that they can be plotted for comparison.


```{r warning = FALSE, fig.width = 6, fig.height = 4, fig.align = "center"}
test <- data.frame(sgolay$mass, sgolay$sg, mono_min$baseline)
names(test) <- c("mz", "Smoothed", "BaselineCorrected")

plotSpectra(dat = test, mass_dat = "mz", spectra_cols = c("Smoothed", "BaselineCorrected"), min_mz = 200, max_mz = 5000,
            x_ticks = 5, intensity_scale = "fixed")
```


Observing the plot above, we can see that the baseline of the spectrum now runs along the *x*-axis, and is no longer irregular. 

__NOTE:__ Depending on the size of your data set, this function may take quite a while to run. Spectral data frames longer than 100,000 rows may take an hour to run. Methods `"linear"` and `"loess"` are capable of processing large data sets more efficiently.

## Linear Interpolation

Linear interpolation estimates the baseline of a spectrum by dividing it into small segments [1]. In each segment, the mean of the points is evaluated and recorded as a baseline predictor for that segment [1].  

A baseline is then generated by linearly interpolating baseline predictors across all small segments [1]. In each segment, the intensity of each peak is subtracted by the intensity of the baseline predictor in that segment.

We will use this method of baseline correction on our `mov_avg` spectrum from section 2.2 below. When using `baselineCorr(method = "linear")`, the `n` argument, segment size, is required. To be consistent with previous examples, we will use `n = 7`.


```{r}
linear_int <- baselineCorr(dat = mov_avg, mass_dat = "mass", 
                           intensity_dat = "mov_avg", method = "linear", n = 7)
```


In order to compare the smoothed and baseline-corrected spectra, we will once again merge their intensities into a single data frame. Once the data has been organized, we can plot it for comparison.


```{r warning = FALSE, fig.width = 6, fig.height = 4, fig.align = "center"}
test <- data.frame(mov_avg$mass, mov_avg$mov_avg, linear_int$baseline)
names(test) <- c("mz", "Smoothed", "BaselineCorrected")

plotSpectra(dat = test, mass_dat = "mz", spectra_cols = c("Smoothed", "BaselineCorrected"), min_mz = 200, max_mz = 5000,
            x_ticks = 7)
```


Again, we can see the baseline now falls along the *x*-axis and is no longer irregular.

## LOESS

Much like the linear interpolation method, LOESS curve fitting begins by dividing a spectrum into small segments [1]. The quantile is then evaluated in each segment [1]. The baseline is estimated as follows:

* If the intensity of point *A* in a segment is less than the quantile, the intensity of point *A* is recorded as a baseline predictor for that point [1]. 
* If the intensity of point *A* is greater than or equal to the quantile in the segment, the quantile for the segment is recorded the baseline predictor for that point [1].

Local polynomial regression fitting is applied to the predictor points [1]. The baseline is then subtracted from the raw spectral signal.

The LOESS method of baseline correction is demonstrated below on our sample raw spectrum, `bsline`.


```{r}
loess <- baselineCorr(dat = bsline, mass_dat = "mass", 
                      intensity_dat = "raw", method = "loess")
```


The raw and baseline corrected spectra are merged into a single data frame and plotted below.


```{r warning = FALSE, fig.width = 6, fig.height = 4, fig.align = "center"}
test <- data.frame(bsline$mass, bsline$raw, loess$baseline)
names(test) <- c("mz", "RawSpectrum", "BaselineCorrected")

plotSpectra(dat = test, mass_dat = "mz", spectra_cols = c("RawSpectrum", "BaselineCorrected"), min_mz = 200, max_mz = 5000,
            x_ticks = 7)
```


Once more, we see the baseline of the corrected spectrum running along the *x*-axis, with the irregularities of the raw spectral baseline removed.

---

# Peak Detection Methods

Spectral data sets can be further processed using peak detection methods. In this step, false peak candidates are removed from a spectrum if they do not meet the criteria in the method.

`subMALDI` offers two methods for peak detection: signal-to-noise ratio and slopes of peaks. Examples can be found below, in sections 4.1 and 4.2.

`subMALDI` is also compatible with other R packages designed for processing mass spectra. To demonstrate this, we will use the continuous wavelet transform (CWT) algorithm from `MassSpecWavelet` [3] to baseline correct and detect peaks in `subMALDI`'s sample data. This is shown in section 4.3.


## Signal-to-Noise (SNR)

Each spectrum is divided into segments of size `n`. Noise is calculated as the median absolute deviation of points within each segment [1]. 

The intensity of each peak is then divided by the noise in that segment, yielding an SNR value for each peak. If the SNR value of a peak is lower than the defined `SNR_thresh`, the peak candidate is discarded [1]. 

We will perform SNR peak detection on our monotone minimum baseline corrected data set from section 3.1, `monotone_min`, below.


```{r}
snr <- peakDet(mono_min, "mz", "baseline", method = "snr", 
                  n = 7, SNR_thresh = 3)
```


As performed previously, we will merge the intensity data from the baseline corrected and peak-filtered spectra. Once the data is in a single frame, we can plot the spectra for comparison.


```{r warning = FALSE, fig.width = 6, fig.height = 4, fig.align = "center"}
test <- data.frame(mono_min$mz, mono_min$baseline, snr$peaks)
names(test) <- c("mz", "Baseline", "Peaks")

plotSpectra(dat = test, mass_dat = "mz", spectra_cols = c("Baseline", "Peaks"), min_mz = 200, max_mz = 5000,
            x_ticks = 7)
```


We can see that some peaks present in the baseline corrected spectrum are not present in the peak-filtered spectrum. These peaks did not meet the SNR criteria.

## Slopes of Peaks

This method uses the shapes of peaks to remove false peak candidates [1]. First, the left and right endpoints of each peak are identified on the baseline. Next, the slopes of each endpoint are evaluated. 

If the either the left or right slope are less than a defined threshold, the peak candidate is discarded [1]. The threshold is defined as half of the local noise level, or half of the median absolute deviation in a window of size `n` [1].

We will demonstrate this method of peak detection below, using the linear interpolated data set generated in section 3.2, `linear_int`.


```{r}
slopes <- peakDet(linear_int, "mz", "baseline", method = "slopes",
                    n = 7)
```



The filtered peaks and baseline-corrected intensities are merged and plotted below.


```{r warning = FALSE, fig.width = 6, fig.height = 4, fig.align = "center"}
test <- data.frame(linear_int$mz, linear_int$baseline, slopes$peaks)
names(test) <- c("mz", "Baseline", "Peaks")

plotSpectra(dat = test, mass_dat = "mz", spectra_cols = c("Baseline", "Peaks"), min_mz = 200, max_mz = 5000,
            x_ticks = 7)
```


Once more, we can see the difference in peak distribution between the baseline corrected and peak-filtered spectra. Peaks that are not shown in the peak-filtered spectrum were discarded as the slopes of their endpoints failed to meet a threshold.

## Continuous Wavelet Transform (CWT)

Several other R packages are available for processing mass spectrometry data [3-5]. While these differ from `subMALDI` in their capabilities and organization, several of the functions from these packages can be used with `subMALDI`'s data structure.

Here we will demonstrate the compatibility of `subMALDI` with the continuous wavelet transform (CWT) algorithm included with `MassSpecWavelet`. The CWT algorithm is a more sophisticated method of pre-processing which is capable of both baseline correction and peak detection in a single step [1].

The package `MassSpecWavelet` is loaded below.

```{r}
library(MassSpecWavelet)
```


We will begin the CWT method by selecting our scales, or scaling factors. The scale factor compresses or dilates a signal. 

When the scale factor is low, the signal is compressed, resulting in a more detailed spectrum. Conversely, when the scale factor is high, the signal is elongated and the resulting spectrum is less detailed.

In the continuous wavelet transform, the scaling factor is usually a vector. This allows the algorithm to compare each spectral signal to different sized wavelets for peak identification. While a wider variety of scaling factors generally results in more identified peaks, a long scaling factor vector may significantly impact the run time of the analysis. 

In the example below, our scales have been selected to ensure maximum peak identification with minimum run time.


```{r}
scales <- seq(1,40, by = 1)
```


`MassSpecWavelet` recognizes spectral data sets as a vector of intensity values. To match this data structure, we will grab only the intensity column of `bsline` and assign it to the variable `y`. For later use, we will also grab the *m/z* column from `bsline` and assign it to the variable `x`.


```{r}
y <- bsline$raw
x <- bsline$mass
```


Once our spectral data has been vectorized for use with `MassSpecWavelet`, we can generate our CWT coefficients using `scales`. This is performed using the code below.


```{r}
coeff <- cwt(y, scales = scales, wavelet = "mexh")
coeff <- cbind(y,coeff)
colnames(coeff) <- c(0,scales)
```


The next step is to identify local maxima in our coefficient matrix. This is done as follows.


```{r}
localMax <- getLocalMaximumCWT(coeff) 
```
  

Once identified, the local maximal coefficients are used to identify ridge lines [1]. To generate each line, the local maximum coefficient is connected to its adjacent scale [1].
  
  
```{r}
ridgeList <- getRidge(localMax)
```


Ridge lines are then used to discard false peak candidates [1]. Criteria for exclusion are as follows:
* If the length of the ridge line is smaller than a user-defined threshold, the peak is discarded.
* If the width of a peak is not within a given range, it is discarded. [3]

We will use these lines to generate a peak index below.


```{r}
majorPeakInfo <- identifyMajorPeaks(y, ridgeList, coeff, 
                                    scales = as.numeric(colnames(coeff)),
                                    SNR.Th = 3)
peakIndex <- majorPeakInfo$peakIndex
```


Now that our peak index has been generated, we can flag unwanted peaks in our original data set for removal. We will create a new data frame containing the original *m/z* data in the first column, intensity data in the second column, and a peak index in the third column, which has been assigned a value of `FALSE`. We will call this data frame `peaks`.


```{r}
peaks <- data.frame(x,y,"Peaks" = FALSE)
```


Using the `for` loop below, data points present in `peakIndex` are flagged as true peaks in our `peaks` data frame. 

This is done by changing the value in the "Peaks" column from `FALSE `to `TRUE`, at rows indicated by `peakIndex`.


```{r}
for(i in peakIndex){
  peaks[i, "Peaks"] <- TRUE
}
```


Once all peaks are identified and false candidates are flagged for removal, we can begin to discard unwanted signals. This is performed below, by making the intensity value ("y") 0 wherever the "Peaks" column is `FALSE`. 


```{r}
peaks[which(peaks$Peaks == FALSE), "y"] <- 0
```


We can now remove the "Peaks" column, as our data has already been manipulated. We will do this by selecting only our `x` and `y` columns, and assigning them to a new data frame called `cwt`.


```{r}
cwt <- dplyr::select(peaks, x, y)
```


Finally, we will merge the raw and CWT-processed data for comparison by plotting.

```{r warning = FALSE, fig.width = 6, fig.height = 4, fig.align = "center"}
test <- data.frame(bsline$mass, bsline$raw, cwt$y)
names(test) <- c("mz", "RawSpectrum", "CWT")

plotSpectra(dat = test, mass_dat = "mz", spectra_cols = c("RawSpectrum", "CWT"), min_mz = 200, max_mz = 5000,
            x_ticks = 7)
```


We can see that the algorithm has successfully corrected the irregular baseline in the raw spectrum, and has reduced the the amount of peaks relative to previous plots. 

---

# Conclusion

After users have completed the pre-processing pipeline on all individual raw spectra, they can now proceed with mapping spectral replicates to a single spectral data frame, as is demonstrated in the *"`subMALDI` Workflow"* vignette. These replicates can then be normalized, averaged, subtracted, and plotted using `subMALDI`'s suite of functions, or further processed or analyzed using other compatible R pacakges!

---

# References
1. C. Yang, Z. He, W. Yu, Comparison of public peak detection algorithms for MALDI mass spectrometry data analysis, *BMC Bioinformatics*, 2009, **10**, 4. 
2. A. Savitzky, M.J.E. Golay, Smoothing and differentiation of data by simplified least squares procedures, *Analytical Chemsitry*, 1964, **36(8)**, 1627-1639.
3. P. Du, W.A. Kibbe, S.M. Lin, Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching, *Bioinformatics*, 2006, **22**, 2059-2065.
4. S. Gibb, K. Strimmer, MALDIquant: a versatile R pacakge for the analysis of mass spectrometry data, *Bioinformatics*, 2012, **28(17)**, 2270-2271.
5. L. Gatto, K.S. Lilley, MSnbase - an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing, and quantitation, *Bioinformatics*, 2012, **28(2)**, 288-289.


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapSpectrum.R
\name{mapSpectrum}
\alias{mapSpectrum}
\title{Map Spectrum to \emph{m/z} Vector}
\usage{
mapSpectrum(
  dat,
  massCol,
  intenseCol,
  dig = 4,
  thresh = 1e-04,
  spec_df,
  colName
)
}
\arguments{
\item{dat}{A pairwise data frame containing your spectral data. Should
contain two columns: one for \emph{m/z} and one for intensity.}

\item{massCol}{A character string; the name of the \emph{m/z} column in
\code{dat}.}

\item{intenseCol}{A character string; the name of the intensity column in
\code{dat}.}

\item{dig}{Number of decimal places to round the \emph{m/z} data to; must
match the same value used in \code{createSpecDF} in order for the columns
to fill. Default = 4.}

\item{thresh}{Single numeric value; all \emph{m/z} values within
\code{thresh} of each other are binned under that with the maximum
intensity. Default = 5e-5.}

\item{spec_df}{An empty data frame created using \code{createSpecDF}.}

\item{colName}{A character string; the name of the column that should be
filled with the spectral data.}
}
\value{
Returns a vector that is used to fill in the \code{colName} column
of the mapped spectral data frame.
}
\description{
Fills in the columns of the empty data frame created using
\code{createSpecDF}. Mass to charge values from pairwise spectral data are
compared to the \code{full_mz} vector. All peaks within the \code{thresh}
of one another are binned, and only the maximum intensity of that bin is
filled into the mapped spectrum.
}
\section{Warning}{
 It is important that the values for \code{thresh} and
\code{dig} are equal to that of \code{res} and \code{dig} used in
\code{createSpecDF()}. Otherwise the data will fail to map.
}

\examples{

## Load sample dataset "Blank1.rda"
data("lank1")

## Create empty spectral data frame to map to
spec_df <- createSpecDF(min_mz = 53.76, max_mz = 1100, res = 0.0001, dig = 4)

## Map binary spectral data to empty spectral data frame
spec_df <- mapSpectrum(dat = Blank1, massCol = "mass", intenseCol = "Intensity", 
            spec_df = spec_df, colName = "Sample", thresh = 1e-04, dig = 4)


}
\references{
https://github.com/wesleyburr/subMaldi
}
\seealso{
\code{\link{createSpecDF}}
}
\author{
Kristen Yeh <kristenyeh@trentu.ca> Wesley Burr <wesleyburr@trentu.ca>
}
\keyword{array}
\keyword{methods}
\name{norm_custom}
\alias{norm_custom}
\title{
Normalize Peak Intensity to Custom \emph{m/z}
}
\description{
Called internally by \code{normSpectra}. Normalizes the intensity data of each input spectrum to the intensity of the selected \emph{m/z} peak. Capable of normalizing 6 spectra to the same \emph{m/z} value at once.
}
\usage{
norm_custom(dat, mass_dat, norm_mz, spectra_cols, showHI = FALSE)
}
\arguments{
  \item{dat}{The name of the spectral data frame, containing \code{m/z} data in the first column and spectral intensity data in subsequent columns.
}
  \item{mass_dat}{A character string; the name of the column in \code{dat} containing the \emph{m/z} data for the spectrum.
}
  \item{norm_mz}{Numeric; the \emph{m/z} peak to which the spectral intensity should be normalized to. Value should have the same number of decimal places as the \emph{m/z} data in \code{dat}.
}
  \item{spectra_cols}{A character string; the names of the column in \code{dat} containing the intensity data for the spectrs to be analyzed.}
  \item{showHI}{A Boolean variable; clipping argument for results greater than 1 after normalization.}
}
\value{Returns a new data frame including the original \emph{m/z} data and normalized intensity data.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca>
Wesley Burr <wesleyburr@trentu.ca>
}
\seealso{
\code{\link{normSpectra}}
}
\examples{
## Load sample dataset "Master2.rda"
data("Master2")

## Normalize the intensity of peaks in "Before1" to the intensity of the peak at 255.23
ex <- normSpectra(dat = Master2, mass_dat = "full_mz", method = "custom", norm_mz = 255.23, 
            spectra_cols = "Before1")

## Normalize the intensity of peaks in 6 spectra to the intensity of the peak at 255.23
ex <- normSpectra(dat = Master2, mass_dat = "full_mz", method = "custom", 
            norm_mz = 255.23, spectra_cols = c("Blank1", "Blank2", "Before1", 
            "Before2", "After1", "After2"))

}
\keyword{ methods }
\keyword{ manip }
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peakDet.R
\name{peakDet}
\alias{peakDet}
\title{Peak Detection}
\usage{
peakDet(
  dat,
  mass_dat,
  intensity_dat,
  method = NULL,
  n = NULL,
  SNR_thresh = NULL
)
}
\arguments{
\item{dat}{The name of the spectral data frame, containing \code{m/z} data
in the first column and spectral intensity data in subsequent columns.}

\item{mass_dat}{Character string. The name of the column in \code{dat}
containing the \emph{m/z} data for the spectrum.}

\item{intensity_dat}{Character string. The name of the column in \code{dat}
containing the intensity data for the spectrum.}

\item{method}{Character string. The method of peak detection. Either
\code{"snr"} for signal-to-noise ratio, or \code{"slopes"} for slopes of
peaks.}

\item{n}{Single numeric value. For both \code{method = "snr"} and
\code{method = "slopes"}, the window size used to calculate noise. Noise is
defined as the median of the absolute deviation (MAD) of points within a
window [1].}

\item{SNR_thresh}{Single numeric value. When \code{method = "snr"}, the
signal-to-noise ratio (SNR) threshold for discarding peaks. If the SNR of a
peak falls below this threshold, it will be discarded.}
}
\value{
Returns a new data frame containing only the peaks which have
passed the detection criteria.
}
\description{
Two methods for peak detection in baseline corrected spectral data. Methods
include signal-to-noise ratio and slopes of peaks.
}
\section{Methods}{
 \describe{ \item{snr}{ Each spectrum is divided into
segments of size \code{n}. Noise is calculated as the median absolute
deviation of points within each segment [1]. If the intensity of a peak
divided by the noise in that segment is less than the indicated
\code{SNR_thresh}, the peak is discarded. } \item{slopes}{ Uses the shapes
of peaks to remove false peak candidates [1]. First, the left and right
endpoints of each peak are identified on the baseline. Next, the slopes of
each endpoint are evaluated. If the either the left or right slope are less
than a defined threshold, the peak candidate is discarded [1]. The
threshold is defined as half of the local noise level, or half of the
median absolute deviation in a window of size \code{n}.} }
}

\examples{

## Load sample dataset "bsline"
data("bsline")

## Baseline correct using method "linear"
linear <- baselineCorr(bsline, "mass", "raw", 
                        method = "linear", n = 7)

## Detect peaks using method "snr"
snr <- peakDet(linear, "mz", "baseline", method = "snr", 
                  n = 7, SNR_thresh = 3)

## Detect peaks using method "slopes"
slopes <- peakDet(linear, "mz", "baseline", method = "slopes",
                    n = 7)

}
\references{
https://github.com/wesleyburr/subMaldi (1) Yang, C., He, Z. &
Yu, W. Comparison of public peak detection algorithms for MALDI mass
spectrometry data analysis. BMC Bioinformatics 10, 4 (2009).
https://doi.org/10.1186/1471-2105-10-4
}
\seealso{
\code{\link{smoothSpectrum}}, \code{\link{baselineCorr}}
}
\author{
Kristen Yeh <kristenyeh@trentu.ca> Wesley Burr <wesleyburr@trentu.ca>
}
\keyword{manip}
\keyword{methods}
\docType{package}
\name{subMALDI-package}
\alias{subMALDI}
\alias{subMALDI-package}
\title{subMALDI: Organization and Processing of MALDI-MS Datasets}
\description{
To learn more about subMALDI, start with the vignettes:
\code{browseVignettes(package = "subMALDI")}
}
\seealso{
Useful links:
\itemize{
  \item Report bugs at \url{https://github.com/wesleyburr/subMALDI/issues}
}

}
\author{
\strong{Maintainer}: Wesley Burr \email{wesleyburr@trentu.ca} (\href{https://orcid.org/0000-0002-2058-1899}{ORCID})

Authors:
\itemize{
  \item Kristen Yeh (\href{https://orcid.org/0000-0002-3411-6816}{ORCID})
  \item Wesley Burr (\href{https://orcid.org/0000-0002-2058-1899}{ORCID})
  \item Sophie Castel (\href{https://orcid.org/0000-0001-9086-0917}{ORCID})
}
}

\name{norm_median}
\alias{norm_median}
\title{
Normalize a Spectral Dataset to its Median Intensity
}
\description{
Called internally by \code{normSpectra()}. Evaluates the median of each spectrum in a dataset after removing 0 values introduced by mapping. If the medians are not equal between spectra, a normalization factor is applied to the spectral intensities until all median intensities in the dataset are equal. Capable of normalizing 2-6 spectra at once.
}
\usage{
norm_median(dat, mass_dat, spectra_cols)
}
\arguments{
  \item{dat}{The name of the spectral data frame, containing \code{m/z} data in the first column and spectral intensity data in subsequent columns.}
  \item{mass_dat}{A character string; the name of the column in \code{dat} containing the \emph{m/z} data for the spectrum.}
  \item{spectra_cols}{A character string; the names of the column in \code{dat} containing the intensity data for the spectrs to be analyzed.}
}
\value{
Returns a new data frame including the original \emph{m/z} data and normalized intensity data.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca>
Wesley Burr <wesleyburr@trentu.ca>
}
\seealso{
\code{\link{normSpectra}}
}
\examples{
## Load sample dataset "Master2.rda"
data("Master2")

## Normalize spectra "Before1" and "Before2" to their median intensity
ex <- normSpectra(dat = Master, mass_dat = "full_mz", method = "median", 
            spectra_cols = c("Before1", "Before2"))

}

\keyword{ methods }
\keyword{ manip }
\name{Before2}
\alias{Before2}
\docType{data}
\title{
Pairwise Spectral Data Frame from Sample
}
\description{
Raw spectral data in binary dataframe. The first column contains the \emph{m/z} data for the spectrum, while the second contains the intensity data. Spectrum is obtained from a biological sample, prior to chemical intervention.
}
\usage{data("Before2")}
\format{
  A data frame with 213593 observations on the following 2 variables.
  \describe{
    \item{\code{mass}}{A numeric vector; the \emph{m/z} data for the spectrum.}
    \item{\code{Intensity}}{A numeric vector; the intensity data for the spectrum.}
  }
}
\source{
Yeh, K., Stock N. L., Burr, W. & Stotesbury, T. Preliminary analysis of latent fingerprints recovered from underneath bloodstains using Matrix-Assisted Laser Desportion/Ionization Fourier-Transform Ion Cyclotron Resonance Mass Spectrometry Imaging (MALDI FT-ICR MSI). Forensic Chemistry (2020). In press.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\examples{
data(Before2)
}
\keyword{datasets}
\name{norm_max}
\alias{norm_max}
\title{
Normalization Method: Maximum Intensity
}
\description{
Called internally by \code{normSpectra}. Normalizes the intensity data of spectra to a scale of 0,1. Capable of normalizing 1-6 spectra, to their individual maxima, at once.
}
\usage{
norm_max(dat, mass_dat, spectra_cols)
}
\arguments{
  \item{dat}{The name of the spectral data frame, containing \code{m/z} data in the first column and spectral intensity data in subsequent columns.}
  \item{mass_dat}{A character string; the name of the column in \code{dat} containing the \emph{m/z} data for the spectrum.}
  \item{spectra_cols}{A character string; the names of the column in \code{dat} containing the intensity data for the spectrs to be analyzed.}
}
\value{Returns a new data frame including the original \emph{m/z} data and normalized intensity data.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca>
Wesley Burr <wesleyburr@trentu.ca>
}
\seealso{
\code{\link{normSpectra}}
}
\examples{
## Load sample dataset "Master2.rda"
data("Master2")

## Normalize intensity data for three spectra
normSpectra(dat = Master2, mass_dat = "full_mz", method = "max", 
            spectra_cols = c("Blank1", "Before1", "After1"))
}
\keyword{ methods }
\keyword{ manip }
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotSpectra.R
\name{plotSpectra}
\alias{plotSpectra}
\title{Plot and Compare Spectra}
\usage{
plotSpectra(
  dat,
  mass_dat,
  spectra_cols,
  palette = NULL,
  colours = (grDevices::colorRampPalette(c("#cc6600",
    "#33ccff")))(length(spectra_cols)),
  span = 5,
  thresh = 0.1,
  lbls = FALSE,
  lbl.fmt = "\%3.4f",
  min_mz = min(dat[[mass_dat]]),
  max_mz = max(dat[[mass_dat]]),
  min_I = 0,
  max_I = max(dat[spectra_cols]),
  x_ticks = 100,
  nrows = ceiling(length(spectra_cols)/2),
  intensity_scale = "free_y"
)
}
\arguments{
\item{dat}{The name of the spectral data frame, containing the \emph{m/z}
data in the first column.}

\item{mass_dat}{A character string; the name of the column in \code{dat}
containing the \emph{m/z} data.}

\item{spectra_cols}{A character string; the name(s) of the column(s) in \code{dat}
containing the intensity data for the spectra-of-interest.}

\item{palette}{A character element; the RColorBrewer palette to use. See below for available palettes.}

\item{colours}{A character string indicating the desired colour(s)}

\item{span}{Single numeric value; the span of peak maxima in between each
label. Default = 5 (ignores two peak maxima on either side of each label).}

\item{thresh}{Single numeric value (0-100); the threshold of peak
intensities which should be labeled.}

\item{lbls}{Logical. If \code{lbls = TRUE}, labels indicating the
\emph{m/z} value of each peak maxima within the indicated \code{span} will
be included in the output plot. If \code{lbls = FALSE}, labels are not
shown.}

\item{lbl.fmt}{Character string in the format "\%a.bf", where \emph{a} is
the number of figures displayed to the left of decimal places in labels,
and \emph{b} is the number of figures displayed to the right of decimal
places in labels.  Default = "\%3.4f".}

\item{min_mz}{Single numeric value; minimum \emph{m/z} value of the
observed range.}

\item{max_mz}{Single numeric value; upper end of \emph{m/z} range observed
in spectra.}

\item{min_I}{Single numeric value; minimum intensity value of the observed
range.}

\item{max_I}{Single numeric value; upper end of the intensity range
observed in spectra.}

\item{x_ticks}{Single numeric value; the number of ticks on the x-axis.}

\item{nrows}{Single integer value; the number of rows in plot arrangement.}

\item{intensity_scale}{A character string; the method that should be used
for determining the y-axis scales for each spectrum. If \code{method =
"free_y"}, each spectrum will be plotted with its own intensity scale. If
\code{method = "fixed"}, each spectrum will be plotted with the y-axis of
the most intense spectrum in the set.}
}
\value{
Returns a line plot of the input spectra.
}
\description{
Plot single or multiple spectra. For multiple spectra, a grid layout can be called using 'nrows'.
}
\section{RColorBrewer Palettes}{
 
\describe{\itemize{\item{\code{"Accent"}}
                   \item{\code{"Dark2"}} 
                   \item{\code{"Paired"}}  
                   \item{\code{"Pastel1"}}  
                   \item{\code{"Pastel2"}}  
                   \item{\code{"Set1"}}  
                   \item{\code{"Set2"} (default)} 
                   \item{\code{"Set3"}}
                   } 
          }
}

\examples{

## Plotting using the sample dataset "Master.rda"
data("Master")
plotSpectra(dat = Master, mass_dat = "full_mz",
            spectra_cols = c("Blank1", "Blank2"),
            intensity_scale = "free_y", lbls = TRUE, nrows = 2, x_ticks = 10)
}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca> Wesley Burr <wesleyburr@trentu.ca> Sophie Castel <sophie.castel@ontariotechu.net>
}
\keyword{aplot}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/avgSpectra.R
\name{avgSpectra}
\alias{avgSpectra}
\title{Average Spectral Replicates}
\usage{
avgSpectra(dat, method = "mean", spectra_cols)
}
\arguments{
\item{dat}{The mapped spectral data frame, containing \code{full_mz} in the
first column.}

\item{method}{A character string; the method used to combine the spectra.
Methods include "sum" and "mean". Default = "mean."}

\item{spectra_cols}{A character vector; the names of the column in
\code{dat} containing the intensity data for the spectra-of-interest.}
}
\value{
Returns a new column in the input data frame containing the
averaged intensity data.
}
\description{
Combines spectral replicates either by averaging (method = "mean") or
summing (method = "sum") the intensity values across each row representing a
mass-to-charge value in \code{full_mz}.
}
\examples{

## Load sample dataset "Master.rda"
data("Master")

## Average blank spectrum 1 and 2 using the method "mean"
ex <- avgSpectra(Master, method = "mean", spectra_cols = c("Blank1", "Blank2"))

## Average blank spectrum 1 and 2 using the method "sum"
ex <- avgSpectra(Master, method = "sum", spectra_cols = c("Blank1", "Blank2"))

}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca> Wesley Burr <wesleyburr@trentu.ca>
Sophie Castel <sophie.castel@ontariotechu.net>
}
\keyword{manip}
\keyword{methods}
\name{norm_stdev}
\alias{norm_stdev}
\title{
Normalize Spectrum to Standard Deviation of Noise
}
\description{
Called internally by \code{normSpectra()}.Evaluates the standard deviation of intensity values within the same noisy region (a region lacking peaks) of each spectrum. All intensities in each spectrum are then divided by the spectrum's standard deviation in the noise region. All output spectra should have a standard deviation of 1 in the selected region. Capable of normalizing 1-6 spectra at once.
}
\usage{
norm_stdev(dat, mass_dat, lower = 900, upper= 1100 , 
                  spectra_cols)
}
\arguments{
  \item{dat}{The name of the spectral data frame, containing \code{m/z} data in the first column and spectral intensity data in subsequent columns.}
  \item{mass_dat}{A character string; the name of the column in \code{dat} containing the \emph{m/z} data for the spectrum.}
  \item{lower}{Single numeric value; the lower \emph{m/z} bound of the noisy region in all spectra.}
  \item{upper}{Single numeric value; the lower \emph{m/z} bound of the noisy region in all spectra.}
  \item{spectra_cols}{A character string; the names of the column in \code{dat} containing the intensity data for the spectrs to be analyzed.}
}
\value{Returns a new data frame including the original \emph{m/z} data and normalized intensity data.
}
\references{
https://github.com/wesleybur/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca>
Wesley Burr <wesleyburr@trentu.ca>
}
\seealso{
code{\link{normSpectra}}
}
\examples{
## Load sample dataset "Master.rda"
data("Master")

## Normalize spectra "After1" and "After2" 
# using the noisy region from m/z 900 to 1100
ex <- normSpectra(Master, mass_dat = "full_mz", method = "stdev",
            lower = 900, upper = 1100, spectra_cols = c("After1", "After2"))
}

\keyword{ methods }
\keyword{ manip }
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normSpectra.R
\name{normSpectra}
\alias{normSpectra}
\title{Normalize Spectral Data}
\usage{
normSpectra(
  dat,
  mass_dat,
  method = NULL,
  norm_mz = NULL,
  upper = NULL,
  lower = NULL,
  spectra_cols = NULL,
  showHI = FALSE
)
}
\arguments{
\item{dat}{The name of the spectral data frame, containing \code{m/z} data
in the first column and spectral intensity data in subsequent columns.}

\item{mass_dat}{A character string; the name of the column in \code{dat}
containing the \emph{m/z} data for the spectrum.}

\item{method}{A character string; the normalization method that should be
used to process the data. See 'Methods' below for list of methods. Default
= NULL.}

\item{norm_mz}{Numeric. If \code{method = "custom"}, the \emph{m/z} peak to
which the spectral intensity should be normalized to. Value should have the
same number of decimal places as the \emph{m/z} data in \code{dat}. If
\code{method = "custom_imprecise", this value must be given as a character
string of numbers.}}

\item{upper}{Numeric. If \code{method = "stdev"}, the upper \emph{m/z}
bound of the noise region of the spectrum.}

\item{lower}{Numeric. If \code{method = "stdev"}, the lower \emph{m/z}
bound of the noise region of the spectrum.}

\item{spectra_cols}{A character vector; the names of the column in
\code{dat} containing the intensity data for the spectra-of-interest.}

\item{showHI}{Logical. To be used with \code{method = "custom"}. If
\code{TRUE}, intensity values greater than \code{norm_mz} will be kept. If
\code{FALSE}, any peaks with intensity greater than \code{norm_mz} will be
truncated to the intensity of \code{norm_mz}. Default = FALSE.}
}
\value{
Returns a new data frame including the original \emph{m/z} data and
normalized intensity data.
}
\description{
Normalize spectral data to a common scale using several different methods.
}
\section{Methods}{
 \describe{ \item{max}{ Normalizes the intensity data of
spectra to a scale of 0,1. } \item{max_set}{ Normalizes the intensity data
of spectra to a scale of 0,1, where 1 is the single most intense peak of
the spectral set.} \item{custom}{ Normalizes the intensity data of each
input spectrum to the intensity of the selected \emph{m/z} peak. }
\item{custom_imprecise}{Normalizes the intensity data of each input
spectrum to the intensity of the selected \emph{m/z} peak. Allows for less
precise normalization; if data contains four decimal places, input
\emph{m/z} values can be input to 2 or 3 decimal places.} \item{TIC}{
Evaluates the sum of all intensities (TIC) of each spectrum in a dataset.
If the TIC of all spectra are not equal, their intensities are multiplied
by a normalization factor. } \item{rel_TIC}{ Evaluates the sum of all
intensities (TIC) of each spectrum in a dataset. If the TIC of all spectra
are not equal, their intensities are multiplied by a normalization factor.
Each peak intensity is then divided by the normalized peak intensity so
that each spectrum in the dataset has a TIC of 1.} \item{RMS}{ Normalizes
each spectrum by dividing each intensity by the spectrum's RMS. }
\item{median}{ Evaluates the median of each spectrum in a dataset after
removing 0 values introduced by mapping. If the medians are not equal
between spectra, a normalization factor is applied to the spectral
intensities until all median intensities in the dataset are equal.}
\item{stdev}{ Evaluates the standard deviation of intensity values within
the same noisy region (a region lacking peaks) of each spectrum. All
intensities in each spectrum are then divided by the spectrum's standard
deviation in the noise region. All output spectra should have a standard
deviation of 1 in the selected region. } \item{quantile}{ Normalizes the
distributions of the values in each spectrum in a set.  Sorts the intensity
data of each spectrum and evaluates the average intensity for each rank.
The intensity values are then replaces with the averaged intensities,
rearranged in their original order. } }
}

\examples{

## Load sample dataset "Master2.rda"
data("Master2")

## Normalize spectrum "Before1" to its maximum intensity
ex <- normSpectra(dat = Master2, mass_dat = "full_mz", 
            method = "max", spectra_cols = "Before1")

## Normalize the spectra "Before1" and "Before2" to the TIC
ex <- normSpectra(dat = Master2, mass_dat = "full_mz",
            method = "TIC", spectra_cols = c("Before1", "Before2"))

## Normalize spectrum "After1" to the intensity of the peak at m/z 253.22
ex <- normSpectra(dat = Master2, mass_dat = "full_mz", 
            method = "custom", norm_mz = 253.22, spectra_cols = "After1")


}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca> Wesley Burr <wesleyburr@trentu.ca>
Sophie Castel <sophie.castel@ontariotechu.net>
}
\keyword{manip}
\keyword{methods}
\name{norm_max_set}
\alias{norm_max_set}
\title{
Normalization Method: Maximum Intensity of Set
}
\description{
Called internally by \code{normSpectra}. Normalizes the intensity data of spectra to a scale of 0,1, where 1 is the single most intense peak of the spectral set. Capable of normalizing 1-6 spectra at once.
}
\usage{
norm_max_set(dat, mass_dat, spectra_cols)
}
\arguments{
  \item{dat}{The name of the spectral data frame, containing \code{m/z} data in the first column and spectral intensity data in subsequent columns.}
  \item{mass_dat}{A character string; the name of the column in \code{dat} containing the \emph{m/z} data for the spectrum.}
  \item{spectra_cols}{A character string; the names of the column in \code{dat} containing the intensity data for the spectrs to be analyzed.}
}
\value{
Returns a new data frame including the original \emph{m/z} data and normalized intensity data of each input spectrum.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca>
Wesley Burr <wesleyburr@trentu.ca>
}
\seealso{
\code{\link{normSpectra}}
}
\examples{
## Load sample datasat "Master2.rda"
data("Master2")

## Normalize the intensity data of two samples to the maximum of the blank
ex <- normSpectra(dat = Master2, mass_dat = "full_mz", method = "max_set", 
            spectra_cols = c("Blank1", "Before1", "Before2"))
}

\keyword{ methods }
\keyword{ manip }
\name{After2}
\alias{After2}
\docType{data}
\title{
Pairwise Spectral Data Frame from Sample
}
\description{
Raw spectral data in binary dataframe. The first column contains the \emph{m/z} data for the spectrum, while the second contains the intensity data. Spectrum is obtained from a biological sample, after chemical intervention.
}
\usage{data("After2")}
\format{
  A data frame with 176209 observations on the following 2 variables.
  \describe{
    \item{\code{mass}}{A numeric vector; the \emph{m/z} data for the spectrum.}
    \item{\code{Intensity}}{A numeric vector; the intensity data for the spectrum.}
  }
}
\source{
Yeh, K., Stock N. L., Burr, W. & Stotesbury, T. Preliminary analysis of latent fingerprints recovered from underneath bloodstains using Matrix-Assisted Laser Desportion/Ionization Fourier-Transform Ion Cyclotron Resonance Mass Spectrometry Imaging (MALDI FT-ICR MSI). Forensic Chemistry (2020). In press.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\examples{
data(After2)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readcsvDir.R
\name{readcsvDir}
\alias{readcsvDir}
\title{Import .csv Spectra from a Directory and Export}
\usage{
readcsvDir(direct, massCol, intenseCol, output)
}
\arguments{
\item{direct}{The path to the directory where the \code{.csv} files are
held.}

\item{massCol}{A character string; the name of the mass column in the
\code{.csv} files.}

\item{intenseCol}{A character string; the name of the intensity column in
the \code{.csv} files.}

\item{output}{The path to the directory where the exported \code{.rda}
files containing the binary spectral data frames should go.}
}
\value{
Returns a directory of \code{.rda} files containing binary spectral
data frames into the \code{output} file path.
}
\description{
Imports all \code{.csv} files in a directory, turns them into a binary data
frame containing a "mass" and "Intensity" column, and outputs them as
\code{.rda} files into an output directory. Allows for rapid import of many
spectral datasets at once.
}
\details{
Assumes that the \code{.csv} files are organized in tidy format, i.e.,
columns are variables and rows are individual m/z observations.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca> Wesley Burr <wesleyburr@trentu.ca>
}
\keyword{methods}
\name{bsline}
\alias{bsline}
\docType{data}
\title{
Raw Spectrum with Irregular Baseline
}
\description{
A spectrum acquired using ESI FT-ICR-MS with an irregular baseline.
}
\usage{data("bsline")}
\format{
  A data frame with 23807 observations on the following 2 variables.
  \describe{
    \item{\code{mass}}{a numeric vector. The \emph{m/z} data of the spectrum.}
    \item{\code{raw}}{A numeric vector. The raw, unfiltered or corrected intensity of the spectrum.}
  }
}

\references{
https://github.com/wesleyburr/subMaldi/
}
\examples{
data(bsline)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subSpectra.R
\name{subSpectra}
\alias{subSpectra}
\title{Subtract Blank Peaks from Sample Spectra}
\usage{
subSpectra(dat, Blank_Var, Sample, Sub_Sample, showNeg = FALSE)
}
\arguments{
\item{dat}{The spectral data frame, containing \code{full_mz} in the first
column and intensity data in the subsequent columns.}

\item{Blank_Var}{A character string; the name of the blank column that will
be subtracted from the sample column.}

\item{Sample}{A character string; the name of the sample column that the
blank data will be subtracted from.}

\item{Sub_Sample}{A character string; the name of the column to be filled
with the subtracted spectrum.}

\item{showNeg}{Logical; if \code{showNeg = TRUE}, then negative values
produced by the subtraction will be kept in the dataset. If \code{showNeg =
FALSE}, then negative values created by the subtraction are set to 0 in the
output dataset.}
}
\value{
Returns a vector of intensity values that are filled into the
\code{Sub_Sample} column of the mapped data frame at the corresponding rows
of \code{full_mz}.
}
\description{
This function takes the intensity values in the blank column of the mapped
spectral data frame and subtracts them from the intensity values in the
sample column. If the intensity value in the blank spectrum is greater than
that in the sample spectrum, the intensity of the sample peak is given an
intensity of 0.
}
\examples{

## Load sample datasets "Master.rda"
data("Master")

## Make subtraction column for new data
Master <- transform(Master, "Subtracted" = 0)

## Subtract Blank1 from Before1
Master <- subSpectra(dat = Master, Blank_Var = "Blank1", 
                      Sample = "Before1", Sub_Sample = "Subtracted")


}
\references{
https://github.com/wesleyburr/subMaldi
}
\seealso{
\code{\link{createSpecDF}}, \code{\link{mapSpectrum}}
}
\author{
Kristen Yeh <kristenyeh@trentu.ca> Wesley Burr <wesleyburr@trentu.ca>
}
\keyword{manip}
\keyword{methods}
\name{Blank1}
\alias{Blank1}
\docType{data}
\title{
Pairwise Spectral Data from Matrix Blank
}
\description{
Raw spectral data in binary dataframe. The first column contains the \emph{m/z} data for the spectrum, while the second contains the intensity data. Spectrum is a MALDI matrix blank.
}
\usage{data("Blank1")}
\format{
  A data frame with 105717 observations on the following 2 variables.
  \describe{
    \item{\code{mass}}{A numeric vector; the \emph{m/z} data for the spectrum.}
    \item{\code{Intensity}}{A numeric vector; the intensity data for the spectrum.}
  }
}
\source{
Yeh, K., Stock N. L., Burr, W. & Stotesbury, T. Preliminary analysis of latent fingerprints recovered from underneath bloodstains using Matrix-Assisted Laser Desportion/Ionization Fourier-Transform Ion Cyclotron Resonance Mass Spectrometry Imaging (MALDI FT-ICR MSI). Forensic Chemistry (2020). In press.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\examples{
data(Blank1)
}
\keyword{datasets}
\name{After1}
\alias{After1}
\docType{data}
\title{
Pairwise Spectral Data Frame from Sample
}
\description{
Raw spectral data in binary dataframe. The first column contains the \emph{m/z} data for the spectrum, while the second contains the intensity data. Spectrum is obtained from a biological sample, after chemical intervention.
}
\usage{data("After1")}
\format{
  A data frame with 134674 observations on the following 2 variables.
  \describe{
    \item{\code{mass}}{A numeric vector; the \emph{m/z} data for the spectrum.}
    \item{\code{Intensity}}{A numeric vector; the intensity data for the spectrum.}
  }
}
\source{
Yeh, K., Stock N. L., Burr, W. & Stotesbury, T. Preliminary analysis of latent fingerprints recovered from underneath bloodstains using Matrix-Assisted Laser Desportion/Ionization Fourier-Transform Ion Cyclotron Resonance Mass Spectrometry Imaging (MALDI FT-ICR MSI). Forensic Chemistry (2020). In press.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\examples{
data(After1)
}
\keyword{datasets}
\name{norm_custimp}
\alias{norm_custimp}
\title{
Normalization Method: Maximum Intensity of Set
}
\description{
Called internally by \code{normSpectra}. Normalizes the intensity data of spectra to a scale of 0,1, where 1 is the single most intense peak of the spectral set; does so against a specified m/z. Capable of normalizing 1-6 spectra at once.
}
\usage{
norm_custimp(dat, mass_dat, norm_mz, spectra_cols, showHI = FALSE)
}
\arguments{
  \item{dat}{The name of the spectral data frame, containing \code{m/z} data in the first column and spectral intensity data in subsequent columns.}
  \item{mass_dat}{A character string; the name of the column in \code{dat} containing the \emph{m/z} data for the spectrum.}
  \item{norm_mz}{Numeric m/z reference to be normalized against.}
  \item{spectra_cols}{A character string; the names of the column in \code{dat} containing the intensity data for the spectrs to be analyzed.}
  \item{showHI}{A Boolean variable; clipping argument for results greater than 1 after normalization.}
}
\value{
Returns a new data frame including the original \emph{m/z} data and normalized intensity data of each input spectrum.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca>
Wesley Burr <wesleyburr@trentu.ca>
}
\seealso{
\code{\link{normSpectra}}
}
\examples{
## Load sample datasat "Master2.rda"
data("Master2")

## Normalize the intensity data of two samples to the maximum of the blank
ex <- normSpectra(dat = Master2, mass_dat = "full_mz", norm_mz = 255.23, 
                  method = "custom_imprecise", spectra_cols = c("Blank1", "Before1", "Before2"))
}

\keyword{ methods }
\keyword{ manip }
\name{norm_TIC}
\alias{norm_TIC}
\title{
Normalize Peak Intensity by Total Ion Current
}
\description{
Called internally by \code{normSpectra}. Evaluates the sum of all intensities (TIC) of each spectrum in a dataset. If the TIC of all spectra are not equal, their intensities are multiplied by a normalization factor. Capable of normalizing 2-6 spectra at once.
}
\usage{
norm_TIC(dat, mass_dat, spectra_cols)
}
\arguments{
  \item{dat}{The name of the spectral data frame, containing \code{m/z} data in the first column and spectral intensity data in subsequent columns.}
  \item{mass_dat}{A character string; the name of the column in \code{dat} containing the \emph{m/z} data for the spectrum.}
  \item{spectra_cols}{A character string; the names of the column in \code{dat} containing the intensity data for the spectrs to be analyzed.}
  }
\value{Returns a new data frame including the original \emph{m/z} data and normalized intensity data.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca>
Wesley Burr <wesleyburr@trentu.ca>
}

\seealso{
\code{\link{normSpectra}}
}
\examples{
## Load sample dataset "Master2.rda"
data("Master2")

# Normalize the TIC of "Before1" and "Before2"
ex <- normSpectra(dat = Master2, mass_dat = "full_mz", method = "TIC", 
            spectra_cols = c("Before1", "Before2"))

}

\keyword{ methods }
\keyword{ manip }
\name{Before1}
\alias{Before1}
\docType{data}
\title{
Pairwise Spectral Data Frame from Sample
}
\description{
Raw spectral data in binary dataframe. The first column contains the \emph{m/z} data for the spectrum, while the second contains the intensity data. Spectrum is obtained from a biological sample, prior to chemical intervention.
}
\usage{data("Before1")}
\format{
  A data frame with 437372 observations on the following 2 variables.
  \describe{
    \item{\code{mass}}{A numeric vector; the \emph{m/z} data for the spectrum.}
    \item{\code{Intensity}}{A numeric vector; the intensity data for the spectrum.}
  }
}
\source{
Yeh, K., Stock N. L., Burr, W. & Stotesbury, T. Preliminary analysis of latent fingerprints recovered from underneath bloodstains using Matrix-Assisted Laser Desportion/Ionization Fourier-Transform Ion Cyclotron Resonance Mass Spectrometry Imaging (MALDI FT-ICR MSI). Forensic Chemistry (2020). In press.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\examples{
data(Before1)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_max_set.R
\name{find_max_set}
\alias{find_max_set}
\title{Find Most Intense Peak of a Set}
\usage{
find_max_set(dat, mass_dat, spectra_cols)
}
\arguments{
\item{dat}{The name of the spectral data frame, containing \code{m/z} data
in the first column and spectral intensity data in subsequent columns.}

\item{mass_dat}{A character string; the name of the column in \code{dat}
containing the \emph{m/z} data for the spectrum.}

\item{spectra_cols}{A character vector; the names of the column in
\code{dat} containing the intensity data for the spectra-of-interest.}
}
\value{
Returns a data frame indidcating the peak maximum of the spectral
set. The sample number, \emph{m/z} value, and intensity data are returned.
}
\description{
Analyzes spectral data and returns information about the most intense peak
in all of the spectral set. Indicates the most intense spectrum, the
\emph{m/z} value, and intesnity of the maximum. User must specify at least
two spectra.
}
\examples{

## Load sample dataset "Master.rda"
data("Master")

## Find the most intense peak of six spectra
find_max_set(dat = Master, mass_dat = "full_mz", spectra_cols = c("After1", "Blank2", 
            "Before1", "Blank1", "After2", "Before2"))

}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca> Wesley Burr <wesleyburr@trentu.ca>
Sophie Castel <sophie.castel@ontariotechu.net>
}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_max.R
\name{find_max}
\alias{find_max}
\title{Find Peak Maxima}
\usage{
find_max(dat, mass_dat, spectra_cols)
}
\arguments{
\item{dat}{The name of the spectral data frame, containing \code{m/z} data
in the first column and spectral intensity data in subsequent columns.}

\item{mass_dat}{A character string; the name of the column in \code{dat}
containing the \emph{m/z} data for the spectrum.}

\item{spectra_cols}{A character vector; the names of the column in
\code{dat} containing the intensity data for the spectra-of-interest.}
}
\value{
Returns a data frame indidcating the most intense peaks of each
input spectrum. Indicates the spectrum the data is from, the \emph{m/z}
value associated with the peak, and the intensity of the maxima.
}
\description{
Analyzes spectral data and returns a list of the most intense peak in each
spectrum, including the \emph{m/z} value associated with the peak.
}
\examples{

## Load sample dataset "Master.rda"
data("Master")


## Find maxima of four spectra
find_max(dat = Master, mass_dat = "full_mz", spectra_cols = c("Blank1", "Before1", "After1", "After2"))

}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca> Wesley Burr <wesleyburr@trentu.ca>
Sophie Castel <sophie.castel@ontariotechu.net>
}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baselineCorr.R
\name{baselineCorr}
\alias{baselineCorr}
\title{Baseline Correction}
\usage{
baselineCorr(dat, mass_dat, intensity_dat, method = NULL, n = NULL)
}
\arguments{
\item{dat}{The name of the spectral data frame, containing \code{m/z} data
in the first column and spectral intensity data in subsequent columns.}

\item{mass_dat}{Character string. The name of the column in \code{dat}
containing the \emph{m/z} data for the spectrum that will be corrected.}

\item{intensity_dat}{Character string. The name of the column in \code{dat}
containing the intensity data for the spectrum to be baseline corrected.}

\item{method}{Character string. The method that is to be used to perform
the baseline correction. Either \code{"monotone"} for monotone minimum
correction, \code{"linear"} for linear interpolation, or \code{"loess"} for
LOESS curve fitting.}

\item{n}{Single odd numeric value. If \code{method = "linear"}, the size of
the window that should be used for linear interpolation of the spectral
baseline.}
}
\value{
Returns a new data frame containing the baseline corrected spectral
data.
}
\description{
Offers three different methods for baseline correction of raw spectral
data. Methods include monotone minimum, linear interpolation, and LOESS
curve fitting.
}
\section{Methods}{
 \describe{ \item{monotone_min}{ Identifies valleys in the
spectrum by evaluating the slopes between adjacent intensity values. Moving
from left to right, if slope of point A is greater than 0, the nearest
point B with slope of less than 0 is located. Points in between A and B are
recorded as potential peaks. If slope of A < 0, nearest point B with slope
> 0 is located. Points between A and B are recorded as valleys. Each
intensity value in the spectrum is then subtracted by the intensity of the
nearest valley, correcting any irregularities in the baseline of the data.
} \item{linear}{ Divides a spectrum into small segments and evaluates the
mean of points in each segment as a valley (1). A baseline is then
generated by linearly interpolating valleys across all small segments. In
each segment, the intensity of each peak is subtracted by the intensity of
the valley in that segment. } \item{loess}{ Divides a spectrum into small
segments and evaluates the quantile in each segment. The baseline is
estimated as follows. If the intensity of point A in a segment is less than
the quantile, the intensity of point A is recorded as a baseline predictor
for that point. If the intensity of point A is greater than or equal to the
quantile in the segment, the quantile for the segment is recorded the
baseline predictor for that point. Local polynomial regression fitting is
applied to the predictor points. The baseline is then subtracted from the
raw spectral signal.  } }
}

\examples{

## Load sample dataset "bsline"
data("bsline")

## Correct the baseline using the monotone minimum method
bsline_mono <- baselineCorr(bsline, "mass", "raw", 
                      method = "monotone_min")

## Correct the baseline using the linear interpolation method
## Window size of 7
bsline_linear <- baselineCorr(bsline, "mass", "raw", 
                      method = "linear", n = 7)

}
\references{
https://github.com/wesleyburr/subMaldi (1) Yang, C., He, Z. &
Yu, W. Comparison of public peak detection algorithms for MALDI mass
spectrometry data analysis. BMC Bioinformatics 10, 4 (2009).
https://doi.org/10.1186/1471-2105-10-4
}
\seealso{
\code{\link{smoothSpectrum}}, \code{\link{peakDet}}
}
\author{
Kristen Yeh <kristenyeh@trentu.ca> Wesley Burr <wesleyburr@trentu.ca>
}
\keyword{manip}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createSpecDF.R
\name{createSpecDF}
\alias{createSpecDF}
\title{Create Empty Spectral Data Frame}
\usage{
createSpecDF(min_mz = 53.76, max_mz = 1100, res = 1e-04, dig = 4)
}
\arguments{
\item{min_mz}{Single numeric value; minimum \emph{m/z} value of the
observed range. Default = 53.76.}

\item{max_mz}{Single numeric value; upper end of \emph{m/z} range observed
in spectra. Default = 1100.}

\item{res}{The resolution of peaks; a single numeric value indicating the
step size from the minimum to maximum \emph{m/z} bin. Also known as bin
width. Default = 0.0001.}

\item{dig}{Number of decimal places to round the \emph{m/z} vector to;
required for \emph{m/z} values to match full \emph{m/z} vector. Default =
4.}
}
\value{
Returns a data frame which can be used to map irregularly spaced
spectral data to a set range of \emph{m/z} values, contained in the first
column of the frame.
}
\description{
Creates an empty data frame to be used for mapping spectral data to a
vector of \emph{m/z} values. The vector of \emph{m/z} values is stored in
the first column and can be called by "\code{full_mz}".
}
\examples{

## Creating an empty spectrum with an m/z range of 500 to 2000 m/z, with a step size of 0.001
spec_df <- createSpecDF(min_mz = 500, max_mz = 2000, res = 0.001, dig = 3)

## Creating an empty spectrum with an m/z range of 100 to 1000, 
# with a step size of 0.0001 and sample names
spec_df <- createSpecDF(min_mz = 100, max_mz = 1000, res = 0.0001, dig = 4)

}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca>
}
\keyword{array}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmveEmpty.R
\name{rmveEmpty}
\alias{rmveEmpty}
\title{Remove Empty Rows}
\usage{
rmveEmpty(dat)
}
\arguments{
\item{dat}{The spectral data frame, containing \code{full_mz} in the first
column, to be trimmed of all-zero rows.}
}
\value{
Returns the input data frame without its all-zero rows.
}
\description{
Data frames created using createSpecDF can have thousands to millions of
rows, especially when dealing with high resolution mass spectrometry data.
This can be quite taxing to the speed of executing certain functions,
especially when it comes to visualizing data. This function removes empty
rows to reduce the computational load and increase the ease of use of
functions on the data frame. Any row where all elements (but the first,
which is \code{full_mz}) equal 0 are removed from the data frame.
}
\note{
Avoid truncating a dataframe until all samples to be compared have
been mapped.
}
\examples{

## Load sample dataset "Mastre.rda"
data("Master")

## Select only the spectra "Before1" and "Before2"
ex <- dplyr::select(Master, "full_mz", "Before1", "Before2")

## Use rmveEmpty(x) on those data frames to reduce computational load for use with
# other functions and packages
ex <- rmveEmpty(dat = ex)

}
\references{
https://github.com/wesleyburr/subMaldi
}
\seealso{
\code{\link{createSpecDF}},\code{\link{mapSpectrum}}
}
\author{
Kristen Yeh <kristenyeh@trentu.ca>
}
\keyword{manip}
\keyword{methods}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readcsvSpec.R
\name{readcsvSpec}
\alias{readcsvSpec}
\title{Import a Single .csv Spectrum}
\usage{
readcsvSpec(spec_file, massCol, intenseCol)
}
\arguments{
\item{spec_file}{A character string; the file path to the \code{.csv}
spectrum that is to be imported.}

\item{massCol}{A character string; the name of the mass column in the
spectrum's \code{.csv} file.}

\item{intenseCol}{A character string; the name of the intensity column in
the spectrum's \code{.csv} file.}
}
\value{
Returns a binary data frame containing the imported data in a
\emph{m/z} column denoted "mass", and an intensity column denoted
"Intensity".
}
\description{
Imports a single spectrum in .csv format and turns it into a binary data
frame containing a "mass" and "Intensity" column.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca> Wesley Burr <wesleyburr@trentu.ca>
}
\keyword{methods}
\name{norm_RMS}
\alias{norm_RMS}
\title{
Normalize Spectral Intensity by Root Mean Square
}
\description{
Called internally by \code{normSpectra}. Normalizes the intensity data of each input spectrum by dividing each peak's intensity by RMS error. Capable of normalizing 6 spectra to the same \emph{m/z} value at once.
}
\usage{
norm_RMS(dat, mass_dat, spectra_cols)
}
\arguments{
  \item{dat}{The name of the spectral data frame, containing \code{m/z} data in the first column and spectral intensity data in subsequent columns.}
  \item{mass_dat}{A character string; the name of the column in \code{dat} containing the \emph{m/z} data for the spectrum.}
  \item{spectra_cols}{A character string; the names of the column in \code{dat} containing the intensity data for the spectrs to be analyzed.}
}
\value{Returns a new data frame including the original \emph{m/z} data and normalized intensity data.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca>
Wesley Burr <wesleyburr@trentu.ca>
}
\seealso{
\code{\link{normSpectra}}
}
\examples{
## Load sample dataset "Master2.rda"
data("Master2")

## Normalize spectrum "Before1" using the RMS method
ex <- normSpectra(dat = Master2, mass_dat = "full_mz", method = "RMS", spectra_cols = "Before1")


}

\keyword{ methods }
\keyword{ manip }
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smoothSpectrum.R
\name{smoothSpectrum}
\alias{smoothSpectrum}
\title{Smooth Noise in Spectrum}
\usage{
smoothSpectrum(
  dat,
  mass_dat,
  intensity_dat,
  method = NULL,
  p = NULL,
  n = NULL,
  m = 0,
  ts = 1
)
}
\arguments{
\item{dat}{The name of the spectral data frame, containing \code{m/z} data
in the first column and spectral intensity data in subsequent columns.}

\item{mass_dat}{A character string; the name of the column in \code{dat}
containing the \emph{m/z} data for the spectrum.}

\item{intensity_dat}{A character string; the name of the column in
\code{dat} containing the intensity data for the spectrum to be smoothed.}

\item{method}{A character string; the method to be used for smoothing.
Available methods include a Savitzky-Golay filter (\code{"sgolay"}) and a
moving average filter (\code{"mov_avg."})}

\item{p}{Single numeric value. If \code{method = "sgolay"}, the filter
order of smoothing. Default = NULL.}

\item{n}{Single odd numeric value. If \code{method = "sgolay"}, the length
of the smoothing filter. If \code{method = "mov_avg"}, the window span
size. Default = NULL.}

\item{m}{Single numeric value. If \code{method = "sgolay"}, returns the
m-th derivative of the filter coefficients. Default = 0.}

\item{ts}{Single numeric value. If \code{method = "sgolay"}, the time
scaling factor. Default = 1.}
}
\value{
Returns a new data frame containing the smoothed spectral data.
}
\description{
Offers two different methods for smoothing noise in raw spectral data: a
moving average filter and the Savitzky-Golay filter (1).
}
\examples{

## Load sample dataset "Before1.rda"
data("Before1")

## Testing method "sgolay"
# test <- smoothSpectrum(dat = Before1, mass_dat = "mass",
#                       intensity_dat = "Intensity", 
#                       method = "sgolay", p = 4, 
#                       n = 25, m = 0)

## Testing method "mov_avg"                        
# test <- smoothSpectrum(dat = Before1, mass_dat = Before1$mass,
#                       intensity_dat = Before1$Intensity,
#                       method = "mov_avg", n = 25)
#                  

}
\references{
https://github.com/wesleyburr/subMaldi (1) A. Savitzky, M.J.E.
Golay, Smoothing and differentiation of data by simplified least-squares
procedures, Anal. Chem. 36 (8) (1964) 1627-1639.
}
\author{
Kristen Yeh <kristenyeh@trentu.ca> Wesley Burr <wesleyburr@trentu.ca>
}
\keyword{manip}
\keyword{methods}
\name{norm_quantile}
\alias{norm_quantile}
\title{
Normalize Distribution of Intensities
}
\description{
Called internally by \code{normSpectra()}. Normalizes the distributions of the values in each spectrum in a set.  Sorts the intensity data of each spectrum and evaluates the average intensity for each rank. The intensity values are then replaces with the averaged intensities, rearranged in their original order. Capable of analyzing 2-6 spectra at once.
}
\usage{
norm_quantile(dat, mass_dat, spectra_cols)
}
\arguments{
  \item{dat}{The name of the spectral data frame, containing \code{m/z} data in the first column and spectral intensity data in subsequent columns.}
  \item{mass_dat}{A character string; the name of the column in \code{dat} containing the \emph{m/z} data for the spectrum.}
  \item{spectra_cols}{A character string; the names of the column in \code{dat} containing the intensity data for the spectrs to be analyzed.}
}
\value{
Returns a new data frame including the original \emph{m/z} data and normalized intensity data.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca>
Wesley Burr <wesleyburr@trentu.ca>
}
\examples{
## Load sample dataset "Master.rda"
data("Master")

## Normalize the spectra "Before1" and "Before2"
# using the quantile method
ex <- normSpectra(dat = Master, mass_dat = "full_mz", method = "quantile",
                  spectra_cols = c("Before1", "Before2"))

}

\keyword{ methods }
\keyword{ manip }
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_ascii.R
\name{read_ascii}
\alias{read_ascii}
\title{Read Raw ASCII Files}
\usage{
read_ascii(filename)
}
\arguments{
\item{filename}{character; The filepath to the local \code{.ascii} data that is to be converted.}
}
\value{
Returns an object of class \code{data.frame} with two numeric columns: \code{"mass"} and \code{"intensity"}
}
\description{
A function to convert raw ASCII files obtained from the Bruker SolariX XR (eXtreme Resolution) FT-ICR Mass Spectrometer.
}
\examples{

## Converting sample ASCII file 'raw_ascii.ascii' to 'data.frame'
file_loc <- system.file("extdata", "raw_test_file.ascii", package = "subMALDI")
converted <- read_ascii(filename = file_loc)

}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Sophie Castel <sophie.castel@ontariotechu.net>
}
\keyword{ascii}
\keyword{convert}
\keyword{data}
\keyword{read}
\name{norm_rel_TIC}
\alias{norm_rel_TIC}
\title{
Normalize Relative Peak Intensity by Total Ion Current
}
\description{
Called internally by \code{normSpectra}. Evaluates the sum of all intensities (TIC) of each spectrum in a dataset. If the TIC of all spectra are not equal, their intensities are multiplied by a normalization factor. Once the TIC of each spectrum has been standardized, the intensity of each peak in all spectra are divided by the normalized TIC. The resulting spectra will have individual TICs of 1. Capable of normalizing 1-6 spectra at once.
}
\usage{
norm_rel_TIC(dat, mass_dat, spectra_cols)
}
\arguments{
  \item{dat}{The name of the spectral data frame, containing \code{m/z} data in the first column and spectral intensity data in subsequent columns.}
  \item{mass_dat}{A character string; the name of the column in \code{dat} containing the \emph{m/z} data for the spectrum.}
  \item{spectra_cols}{A character string; the names of the column in \code{dat} containing the intensity data for the spectrs to be analyzed.}
}
\value{ Returns a new data frame including the original \emph{m/z} data and normalized intensity data.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\author{
Kristen Yeh <kristenyeh@trentu.ca>
Wesley Burr <wesleyburr@trentu.ca>
}
\note{
\code{\link{normSpectra}}
}
\examples{
## Load sample dataset "Master2.rda"
data("Master2")

# Normalize the TIC of "Before1" and "Before2"
ex <- normSpectra(dat = Master2, mass_dat = "full_mz", 
            method = "rel_TIC", spectra_cols = c("Before1", "Before2"))

}

\keyword{ methods }
\keyword{ manip }
\name{Master}
\alias{Master}
\docType{data}
\title{
Master Data Frame (High Resolution)
}
\description{
A mapped spectral data frame created using \code{createSpecDF()} and \code{mapSpectrum()}. Contains the baseline corrected (by linear interpolation) intensity data for 6 spectra, with \emph{m/z} data recorded up to 4 decimal places. 
}
\usage{data("Master")}
\format{
  A data frame with 980506 observations on the following 7 variables.
  \describe{
    \item{\code{full_mz}}{A numeric vector; the \emph{m/z} data for all spectra in the data frame.}
    \item{\code{Blank1}}{A numeric vector; the intensity data for the first spectrum of the dataset. A MALDI matrix blank spectrum.}
    \item{\code{Blank2}}{A numeric vector; the intensity data for the second spectrum in the dataset. A MALDI matrix blank spectrum.}
    \item{\code{Before1}}{A numeric vector; the intensity data for the third spectrum in the dataset. Spectrum obtained prior to chemical intervention.}
    \item{\code{Before2}}{A numeric vector; the intensity data for the fourth spectrum in the dataset. Spectrum obtained prior to chemical intervention.}
    \item{\code{After1}}{A numeric vector; the intensity data for the fifth spectrum in the dataset. Spectrum obtained after chemical intervention.}
    \item{\code{After2}}{A numeric vector; the intensity data for the sixth spectrum in the dataset. Spectrum obtained after chemical intervention.}
  }
}

\source{
Yeh, K., Stock N. L., Burr, W. & Stotesbury, T. Preliminary analysis of latent fingerprints recovered from underneath bloodstains using Matrix-Assisted Laser Desportion/Ionization Fourier-Transform Ion Cyclotron Resonance Mass Spectrometry Imaging (MALDI FT-ICR MSI). Forensic Chemistry (2020). In press.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\examples{
data(Master)

## Separate samples
blanks <- dplyr::select(Master, "full_mz", "Blank1", "Blank2")
precon <- dplyr::select(Master, "full_mz", "Before1", "Before2")
postcon <- dplyr::select(Master, "full_mz", "After1", "After2")
}
\keyword{datasets}
\name{Blank2}
\alias{Blank2}
\docType{data}
\title{
Pairwise Spectral Data from Matrix Blank
}
\description{
Raw spectral data in binary dataframe. The first column contains the \emph{m/z} data for the spectrum, while the second contains the intensity data. Spectrum is a MALDI matrix blank.
}
\usage{data("Blank2")}
\format{
  A data frame with 276704 observations on the following 2 variables.
  \describe{
    \item{\code{mass}}{A numeric vector; the \emph{m/z} data for the spectrum.}
    \item{\code{Intensity}}{A numeric vector; the intensity data for the spectrum.}
  }
}
\source{
Yeh, K., Stock N. L., Burr, W. & Stotesbury, T. Preliminary analysis of latent fingerprints recovered from underneath bloodstains using Matrix-Assisted Laser Desportion/Ionization Fourier-Transform Ion Cyclotron Resonance Mass Spectrometry Imaging (MALDI FT-ICR MSI). Forensic Chemistry (2020). In press.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\examples{
data(Blank2)
}
\keyword{datasets}
\name{Master2}
\alias{Master2}
\docType{data}
\title{
Master Data Frame (Low Resolution)
}
\description{
A mapped spectral data frame created using \code{createSpecDF()} and \code{mapSpectrum()}. Contains the baseline corrected (by linear interpolation) intensity data for 6 spectra, with \emph{m/z} data recorded up to 2 decimal places. 
}
\usage{data("Master2")}
\format{
  A data frame with 91346 observations on the following 7 variables.
  \describe{
    \item{\code{full_mz}}{A numeric vector; the \emph{m/z} data for all spectra in the data frame.}
    \item{\code{Blank1}}{A numeric vector; the intensity data for the first spectrum of the dataset. A MALDI matrix blank spectrum.}
    \item{\code{Blank2}}{A numeric vector; the intensity data for the second spectrum in the dataset. A MALDI matrix blank spectrum.}
    \item{\code{Before1}}{A numeric vector; the intensity data for the third spectrum in the dataset. Spectrum obtained prior to chemical intervention.}
    \item{\code{Before2}}{A numeric vector; the intensity data for the fourth spectrum in the dataset. Spectrum obtained prior to chemical intervention.}
    \item{\code{After1}}{A numeric vector; the intensity data for the fifth spectrum in the dataset. Spectrum obtained after chemical intervention.}
    \item{\code{After2}}{A numeric vector; the intensity data for the sixth spectrum in the dataset. Spectrum obtained after chemical intervention.}
  }
}
\source{
Yeh, K., Stock N. L., Burr, W. & Stotesbury, T. Preliminary analysis of latent fingerprints recovered from underneath bloodstains using Matrix-Assisted Laser Desportion/Ionization Fourier-Transform Ion Cyclotron Resonance Mass Spectrometry Imaging (MALDI FT-ICR MSI). Forensic Chemistry (2020). In press.
}
\references{
https://github.com/wesleyburr/subMaldi
}
\examples{
data(Master2)
data(Master)

## Separate samples
blanks <- dplyr::select(Master2, "full_mz", "Blank1", "Blank2")
precon <- dplyr::select(Master2, "full_mz", "Before1", "Before2")
postcon <- dplyr::select(Master2, "full_mz", "After1", "After2")
}
\keyword{datasets}
