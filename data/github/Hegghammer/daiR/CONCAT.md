---
title: 'daiR: an R package for OCR with Google Document AI'
tags:
  - R
  - optical character recognition
  - cloud computing
  - text mining
  - natural language processing
authors:
  - name: Thomas Hegghammer^[corresponding author]
    orcid: 0000-0001-6253-1518
    affiliation: 1
affiliations:
  - name: Senior Research Fellow, Norwegian Defence Research Establishment (FFI)
    index: 1
date: 31 August 2021
bibliography: paper.bib
---

# Statement of need

Optical character recognition (OCR) promises to open up centuries worth of text to computational analysis. But OCR software has long been sensitive to visual noise and weak on non-Western languages. In April 2021, Google launched Document AI (DAI), a server-based processor offering high-accuracy OCR for over sixty languages [@vanguri:2021]. The `daiR` [@hegghammer:2021] package provides an R interface to the Document AI API along with additional tools for output parsing and visualization.

# Summary

Text as data is a growing field in the social sciences and digital humanities, but computational access to text produced before the late 20th century has been limited by the difficulty of extracting text from document scans. Established OCR libraries such as Tesseract [@tesseract:2021] are highly sensitive to noise and often require extensive corpus-specific adaptations to render text accurately.

The past two years have seen the introduction of server-based OCR processors, such as Amazon Textract [@amazon:2021] and Google Document AI (DAI), which offer very high accuracy out of the box [@hegghammer:2021b]. Of the two, DAI performs better in benchmarking tests and offers broader language support.

In R, where many scholars do their text analysis work, there are packages for Tesseract [@ooms:2021] and Amazon Textract [@kretch:2021], but not for Document AI. The primary objective of `daiR` is therefore to provide access, from within R, to all the main functionalities of the Document AI API. The secondary aim is to offer tools to help parse the output returned by
the DAI processor.

DAI is part of Google Cloud Services (GCS), a suite of cloud computing services for storage, analytics, and machine learning. `daiR` joins a family of existing R packages that interface with GCS, such as `googleLanguageR` [@edmondson:2020] and `googleCloudStorageR` [@edmondson:2021], that together allow for the implementation of multiple GCS tools into an R-based text mining workflow.

`daiR` also includes a range of tools to process DAI's output, which comes in complex JSON files. One set of functions extracts text and table data from the JSON files and brings them into R as character vectors or data frames. Another set draws block, paragraph, line, and token boundary boxes on images of the submitted documents, to help with visual inspection. A third group of functions helps rearrange text blocks in the cases where Document AI has misread their order. Document AI has near-perfect character recognition, but
its parsing of complex page layouts is fallible. This problem is likely to diminish over time as Document AI's algorithm trains on ever larger document data sets. In the meantime, `daiR` makes it relatively easy to correct DAI's errors and obtain an accurately rendered text.

`daiR` is the first R tool to offer high-accuracy text extraction from noisy historical documents out of the box. Until now, scholars have often dealt with Tesseract's high error rates by treating error as noise and using bag-of-words techniques such as topic modeling. Low-error OCR opens up for wider use of natural language processing and other methods that require correctly parsed and ordered text. DAI's improved language coverage may also help reduce the prevalence of English-language data in computational text analysis.

# Acknowledgements
I am grateful to Mark Edmondson, Trond Arne Sørby, Neil Ketchley, and Hallvar Gisnås for contributions to this project and to Christopher Barrie for valuable reviewer comments.

# References
# daiR: OCR with Google Document AI in R

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/daiR)](https://CRAN.R-project.org/package=daiR)
[![R-CMD-check](https://github.com/Hegghammer/daiR/actions/workflows/package-check.yml/badge.svg)](https://github.com/Hegghammer/daiR/actions/workflows/package-check.yml)
[![Codecov test coverage](https://codecov.io/gh/Hegghammer/daiR/branch/master/graph/badge.svg)](https://codecov.io/gh/Hegghammer/daiR?branch=master)
<!-- badges: end -->

**daiR** is an R package for [Google Document AI](https://cloud.google.com/document-ai), a powerful server-based OCR processor with support for over 60 languages. The package provides an inferface for the Document AI API and comes with additional tools for output file parsing and text reconstruction.

<img src="man/figures/frontpage_image.png" width="400" class="center">

## Use

Quick OCR short documents:

```R
## NOT RUN
library(daiR)
response <- dai_sync("file.pdf")
text <- text_from_dai_response(response)
cat(text)
```

Batch process asynchronously via Google Storage:

```R
## NOT RUN
library(googleCloudStorageR)
library(purrr)
my_files <- c("file1.pdf", "file2.pdf", "file3.pdf")
map(my_files, gcs_upload)
dai_async(my_files)
contents <- gcs_list_objects()
output_files <- grep("json$", contents$name, value = TRUE)
map(output_files, ~ gcs_get_object(.x, saveToDisk = file.path(tempdir(), .x)))
sample_text <- text_from_dai_file(file.path(tempdir(), output_files[1]))
cat(sample_text)
```

Turn images of tables into R dataframes:

```R
## NOT RUN:
response <- dai_sync_tab("tables.pdf")
dfs <- tables_from_dai_response(response) 
```

## Requirements

Google Document AI is a [paid service](https://cloud.google.com/document-ai/pricing) that requires a [Google Cloud](https://console.cloud.google.com/) account and a [Google Storage](https://cloud.google.com/storage) bucket. I recommend using Mark Edmondson's `googleCloudStorageR` [package](https://github.com/cloudyr/googleCloudStorageR) in combination with `daiR`. See [vignettes](http://dair.info/) for more on authentication and setup.

## Installation

Download from CRAN:

```R
utils::install.packages("daiR")
```

Or install the latest development version from Github:

```R
devtools::install_github("hegghammer/daiR")
```

# daiR: OCR with Google Document AI in R

<img align="right" src="man/figures/logo.png" width="120">

**daiR** is an R package for [Google Document AI](https://cloud.google.com/document-ai), a powerful server-based OCR processor with support for over 60 languages. The package provides an interface for the Document AI API and comes with additional tools for output file parsing and text reconstruction. See the `daiR` [website](https://dair.info/) and this [journal article](https://joss.theoj.org/papers/10.21105/joss.03538#) for more details.

## Use

Quick OCR short documents:

```R
## NOT RUN
library(daiR)
response <- dai_sync("file.pdf")
text <- text_from_dai_response(response)
cat(text)
```

Batch process asynchronously via Google Storage:

```R
## NOT RUN
library(googleCloudStorageR)
library(purrr)
my_files <- c("file1.pdf", "file2.pdf", "file3.pdf")
map(my_files, gcs_upload)
dai_async(my_files)
contents <- gcs_list_objects()
output_files <- grep("json$", contents$name, value = TRUE)
map(output_files, ~ gcs_get_object(.x, saveToDisk = file.path(tempdir(), .x)))
sample_text <- text_from_dai_file(file.path(tempdir(), output_files[1]))
cat(sample_text)
```

Turn images of tables into R dataframes:

```R
## NOT RUN:
response <- dai_sync_tab("tables.pdf")
dfs <- tables_from_dai_response(response)
```

## Requirements

Google Document AI is a [paid service](https://cloud.google.com/document-ai/pricing) that requires a [Google Cloud](https://console.cloud.google.com/) account and a [Google Storage](https://cloud.google.com/storage) bucket. I recommend using Mark Edmondson's `googleCloudStorageR` [package](https://github.com/cloudyr/googleCloudStorageR) in combination with `daiR`.

## Installation

Download from CRAN:

```R
utils::install.packages("daiR")
```

Or install the latest development version from Github:

```R
devtools::install_github("hegghammer/daiR")
```

## Citation

To cite `daiR` in publications, please use

>Hegghammer, T., (2021). daiR: an R package for OCR with Google Document AI. *Journal of Open Source Software*, 6(68), 3538, https://doi.org/10.21105/joss.03538

Bibtex:
```
@article{Hegghammer2021,
  doi = {10.21105/joss.03538},
  url = {https://doi.org/10.21105/joss.03538},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {68},
  pages = {3538},
  author = {Thomas Hegghammer},
  title = {daiR: an R package for OCR with Google Document AI},
  journal = {Journal of Open Source Software}
}
```

## Acknowledgments

Thanks to Mark Edmondson, Hallvar Gisnås, Will Hanley, Neil Ketchley, Trond Arne Sørby, Chris Barrie, and Geraint Palmer for contributions to the project.

## Code of conduct

Please note that the daiR project is released with a [Contributor Code of Conduct](https://www.contributor-covenant.org/version/2/0/code_of_conduct/). By contributing to this project, you agree to abide by its terms.

<!-- badges: start -->
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03538/status.svg)](https://doi.org/10.21105/joss.03538)
[![CRAN status](https://www.r-pkg.org/badges/version/daiR)](https://CRAN.R-project.org/package=daiR)
[![R-CMD-check](https://github.com/Hegghammer/daiR/actions/workflows/package-check.yml/badge.svg)](https://github.com/Hegghammer/daiR/actions/workflows/package-check.yml)
[![Codecov test coverage](https://codecov.io/gh/Hegghammer/daiR/branch/master/graph/badge.svg)](https://codecov.io/gh/Hegghammer/daiR?branch=master)
<!-- badges: end -->
# daiR 0.9.5

- Added two new functions: `dai_notify()` and `merge_shards()`
- Modified `text_from_dai_response()` and `text_from_dai_file()` to allow saving the output straight to a text file. 
- Fixed a bug in `dai_status()` that caused an error when processing responses from the v1beta2 endpoint (`dai_tab_async()`).  

# daiR 0.9.3

- Changed the "draw_*" functions to allow custom output filenames. 

# daiR 0.9.2

- Added more support documentation in connection with JOSS release.

# daiR 0.9.1

- Added new function (draw_blocks_new()) to inspect block bounding boxes after reprocessing.

# daiR 0.9.0 

- Initial CRAN version.

# daiR 0.8.0 

- Simplified the auth functions to avoid modifying the global environment. Tokens are no longer held in an internal .auth object. Removed dai_has_token(), dai_deauth(), and create_folder().  

# daiR 0.7.0

- Revised draw_blocks(), draw_paragraphs(), draw_lines(), and draw_tokens() functions. These functions no longer require supplying a pdf file, as they get the images from a base64-encoded string in the json file. 

# daiR 0.6.0

- New processing functions adapted to the new stable release of Document AI API (v1). The `dai_sync()`/`dai_async()` functions now access the new v1 endpoint, which has foreign language support. However, the v1 endpoint currently does not support table extraction, so the old processing functions (which access the v1beta2 endpoint) are kept under the new names `dai_sync_tab()` and `dai_async_tab()`. I expect this to be a temporary solution until DAI's capabilities are consolidated in a single endpoint, at which stage the `*_tab()` functions will be phased out.

# daiR 0.4.0

- New table extraction functions. 

# daiR 0.2.0

- New helper functions and more robust code.

# daiR 0.1.0

- Revised and expanded auth functions.

# daiR 0.0.1

- Initial Github release.
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
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
reported to the community leaders responsible for enforcement at
hegghammer@gmail.com (Thomas Hegghammer). All complaints will be reviewed and
investigated promptly and fairly.

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
## How to contribute

### Reporting bugs
If you find a bug in daiR, please report it on the [issue tracker](https://github.com/Hegghammer/daiR/issues/new?assignees=Hegghammer&labels=bug&template=bug_report.md&title=%5BBUG%5D%3A).

### Suggesting enhancements
If you want to suggest a new feature or an improvement of a current feature, you can submit this
on the [issue tracker](https://github.com/Hegghammer/daiR/issues/new?assignees=Hegghammer&labels=enhancement&template=feature_request.md&title=%5BFEATURE%5D%3A).

### Reporting other problems
If you encounter a problem using daiR and are not sure what the source of the problem is, you can report this too on the [issue tracker](https://github.com/Hegghammer/daiR/issues/new?assignees=Hegghammer&labels=help+wanted&template=technical-problem.md&title=%5BPROBLEM%5D%3A)

### Submitting a pull request
If you want to directly submit code to daiR, you can do this by forking the daiR repo and then submitting a pull request.

### Code of conduct
We expect all our contributors to follow [the Contributor Covenant](CODE_OF_CONDUCT.md). Any unacceptable
behaviour can be reported to Thomas (hegghammer@gmail.com).
## Resubmission

This is a resubmission. In this version I have addressed Gregor Seyers comments. I have:

* put 'daiR' in single quotes throughout, except in the top line of the DESCRIPTION
file, as `devtools::check(cran=TRUE)` throws an error if I do. 
* put 'Document AI' in single quotes throughout in the DESCRIPTION file.
* added a web reference for the API in the DESCRIPTION file.
* added `\value` to all .Rd files that didn't have it and reviewed all `\value`
entries to make sure they communicate the structure/class and meaning of the output,
including in the places where no value is returned. 
* removed all instances I could find of functions writing to the user's homespace. 
I checked all the examples, tests, vignettes, as well as readme.md and changed to tempdir()
throughout.
* removed the function that wrote to the global environment. I should mention that 
the function --- which creates an `.auth` object on load to store access tokens --- 
was borrowed from a set of large R packages currently on CRAN, notably 
['bigRQuery'](https://github.com/r-dbi/bigrquery/blob/main/R/zzz.R) and 
['googledrive'](https://github.com/tidyverse/googledrive/blob/master/R/zzz.R).
This led me to believe that CRAN makes exceptions for credential-storing functions.
My new authentication solution works, but in case it breaks, it would be useful to know 
whether CRAN does indeed allow this particular operation. (I'm assuming the
maintainers of the other packages use it for good reason.)

I also made some additional changes. I have:

* removed two functions (`dai_has_token()` and `dai_deauth`) that are redundant under 
the new authentication solution.
* removed one function (`create_folder()`) that I found on closer inspection to be  
unnecessary.
* rewritten several function descriptions (in the .Rd files) for improved clarity 
and consistency.
* revised news.md and the vignettes to reflect the above changes. 
* changed the new version number to 0.9.0 in view of the scale of the combined changes.  

## Test environments
* local Win 10 Enterprise install, R 4.1.0
* windows 10.0.17763 (on Github actions), R 4.1.0
* ubuntu 20.04 (on Github actions), R 4.1.0
* mac OS 10.15 (on Github actions), R 4.1.0
* windows (on WinBuilder), R Devel
* fedora 24 (on rhub), R Devel

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE on rhub and WinBuilder:

* New submission 

## Downstream dependencies
I am not aware of any downstream dependencies.

################################################

## Package history 
This is a first submission.

## Test environments
* local Win 10 Enterprise install, R 4.1.0
* windows 10.0.17763 (on Github actions), R 4.1.0
* ubuntu 20.04 (on Github actions), R 4.1.0
* mac OS 10.15 (on Github actions), R 4.1.0
* windows (on WinBuilder), R Devel
* fedora 24 (on rhub), R Devel

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE on rhub and WinBuilder:

* Possibly mis-spelled words in DESCRIPTION:
  JSON (14:39)
  daiR (13:15, 14:77)
  
  These are proper names. 

## Downstream dependencies
I am not aware of any downstream dependencies.
---
name: Bug report
about: Report what you think is a problem in the code
title: "[BUG]:"
labels: bug
assignees: Hegghammer

---

**Describe the problem**
A concise description of the bug

**To Reproduce**
Steps to reproduce the behavior

** Your configuration**
 - OS: [e.g. Windows 10]
 - R version [e.g. R 4.1.2]
 - IDE type and version: [e.g. RStudio 2021.09.1]

**Additional context and screenshots**
Add any other context about the problem here. If applicable, add screenshots to help explain your problem.
---
name: Feature request
about: Suggest an idea for this project
title: "[FEATURE]:"
labels: enhancement
assignees: Hegghammer

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
name: Technical problem
about: Something's not working and you're not sure why
title: "[PROBLEM]:"
labels: help wanted
assignees: Hegghammer

---

**Describe the problem**
A concise description of the problem

**Summarize the procedure**
A description of the steps you followed until the problem occurred (the more detailed the better)

** Your configuration**
 - OS: [e.g. Windows 10]
 - R version [e.g. R 4.1.2]
 - IDE type and version: [e.g. RStudio 2021.09.1]

**Additional context and screenshots**
Add any other context about the problem here. If applicable, add screenshots to help explain your problem.
