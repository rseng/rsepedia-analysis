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
---
title: "Using Google Document AI with R"
author: "Thomas Hegghammer"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using Google Document AI with R}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
options(rmarkdown.html_vignette.check_title = FALSE)
library(knitr)
opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

**Last updated 14 April 2021**
\
\
\

## About Document AI

[Google Document AI](https://cloud.google.com/document-ai) (DAI) is a server-based OCR engine that extracts text from pdf files. Released in November 2020, it is much more powerful than static libraries such as [`tesseract`](https://github.com/tesseract-ocr/tesseract). Short of corpus-specific, self-trained processors, DAI offers some of the best OCR capabilities currently available to the general public. At the time of writing, DAI is more expensive than Amazon's [Textract](https://aws.amazon.com/textract/), but promises to support many more languages.

DAI is accessed through an API, but this API currently has no official R [client library](https://cloud.google.com/document-ai/docs/libraries). This is where the `daiR` package comes in; it provides a light wrapper for DAI's [REST API](https://cloud.google.com/document-ai/docs/reference/rest), making it possible to submit documents to DAI from within R. In addition, `daiR` comes with pre- and postprocessing tools intended to make the whole text extraction process easier. 

Google Document AI is closely connected with [Google Storage](https://cloud.google.com/storage), as the latter serves as a drop-off and pick-up point for files you want processed in DAI. An R workflow for DAI processing consists of three core steps: 

1. Upload your files to a Google Storage bucket. This can be done manually in the [Google Cloud Console](https://console.cloud.google.com/storage/) or programmatically with the package [`googleCloudStorager`](https://code.markedmondson.me/googleCloudStorageR/index.html). 
2. Using `daiR`, tell DAI to process the files in your bucket. DAI will return its output to your Storage bucket in the form of json files.
3. Download the json files from your Storage bucket to your hard drive. Again you can use either the Cloud Console or `googleCloudStorager`.

## Setup

A [previous vignette](https://dair.info/articles/setting_up_google_storage.html) covered the setting up of a Google Cloud service account and interacting with Google Storage. Here we pick up from where that vignette left off, and assume that the following things are in place:

1. A Google Cloud Services (GCS) **project** linked to your billing account and with the Document AI API enabled.
2. A **service account** with the role "Owner".
3. A **json file** with the service account key, the path to which is stored in an environment variable called `GCS_AUTH_FILE`.
4. The name of your default bucket stored in an environment variable called `GCS_DEFAULT_BUCKET`.

To use Document AI, we need to complete a few more steps. 

### Step 1: Activate Document AI
First, we must activate the API. Go to the [Google Cloud Console](https://console.cloud.google.com/) and open the navigation menu on the left hand side. Click on "APIs and services". Then click on "Enable APIs and Services", type "document ai" in the search field, click on "cloud document ai API", and then "Enable". 

### Step 2: Create a processor
Open the navigation menu on the left again. Scroll down, almost to the bottom, till you see "Document AI" (under the group heading "Artificial intelligence"). Click on "Document AI".

Now click the blue button labelled "Create processor". On the next page, choose the "Document OCR" processor type. A pane should open on your right where you can choose a name for the processor. Call it what you like; the name is mainly for your own reference. Select a location (where you want your files to be processed), then click create.

You should now see a page listing the processor's Name, ID, Status and other attributes. The main thing you want here is the **ID**. Select it and copy it to the clipboard.  

### Step 3: Store the processor id as an environment variable

Open your `.Renviron` file by calling `usethis::edit_r_environ()`. Add `DAI_PROCESSOR_ID="<your processor id>"` on a separate line. Save `.Renviron` and restart RStudio.   

That's it. If these things are in place, you can start processing right after loading the package.

**A note on access tokens**: Unlike some other GCS wrappers, daiR does not authenticate on startup and store access tokens in the environment. Instead it generates tokens on a per-call basis. If you prefer to generate one token per session, you can use `dai_token()` to store your token in an object and pass that object directly into the API call functions using the latter's `token =` parameter. This also means you can use auth functions from pretty much any other GCS wrapper to generate your token.

```{r, eval=FALSE}
library(daiR)
```

Now let's try this thing out.

## Synchronous processing

The quickest and easiest way to OCR with DAI is through synchronous processing. You simply pass an image file or a pdf (of up to 5 pages) to the processor and get the result into your R environment within seconds.

We can try with a sample pdf from the CIA's Freedom of Information Act Electronic Reading Room: 

```{r, eval=FALSE}
setwd(tempdir())
download.file("https://www.cia.gov/readingroom/docs/AGH%2C%20LASLO_0011.pdf", 
              destfile = "CIA.pdf", 
              mode = "wb")
```

We send it to Document AI with `dai_sync()` and store the HTTP response in an object, for example `response`.

```{r, eval=FALSE}
response1 <- dai_sync("CIA.pdf")
```

Then we extract the text with `text_from_dai_response()`:

```{r, eval=FALSE}
text <- text_from_dai_response(response1)
cat(text)
```

Synchronous processing is very convenient, but has two limitations. One is that OCR accuracy may be slightly reduced compared with asynchronous processing, because `dai_sync()` converts the source file to a lightweight, grayscale image before passing it to DAI. The other is scaling; If you have a large pdf or many files, it is usually easier to process them asynchronously.

## Asynchronous processing

In asynchronous (offline) processing, you don't send DAI the actual document, but rather its location on Google Storage so that DAI can process it "in its own time". While slower than synchronous OCR, it allows for batch processing. The `daiR` function `dai_async()` is vectorized, so you can send multiple files with a single call. For this vignette, however, we'll just use a single document; the same as in the previous example.

The first step is to upload the source file(s) to a Google Storage bucket where DAI can find it.^[Note that if you do not have a `GCS_DEFAULT_BUCKET` variable in your .Renviron file, you will need to either set a default bucket for the current session with `gcs_global_bucket("<a bucket name>")` or supply a `bucket = "<a bucket name>"` parameter explicitly inside `gcs_upload()`.]

```{r, eval=FALSE}
library(googleCloudStorageR)
gcs_upload("CIA.pdf")
```

Let's check that our file made it safely: 

```{r, eval=FALSE}
gcs_list_objects()
```

We're now ready to send it off to Document AI with `daiR`'s workhorse function, `dai_async()`, as follows: 

```{r, eval=FALSE}
response2 <- dai_async("CIA.pdf")
```

A few words about this function. Its core parameter, `files`, tells DAI what to process. You can submit either .pdf, .gif, or .tiff files, and your `files` vector can contain a mixture of these three file formats.

You can also specify a `dest_folder`: the name of the bucket folder where you want the output. It defaults to the root of the bucket, but you can specify another subfolder. If the folder does not exist already, it will be created.

The function also takes a location parameter (`loc`), which defaults to "eu" but can be set to "us". It has nothing to do with where you are based, but with which of Google's servers will process your files. The parameter `skip_rev` can be ignored by most; it is for passing selected documents to human review in business workflows. The remaining parameters default to things that are defined by your environment variables (provided you followed the recommendations above).   

Back to our processing. If your call returned "status: 200", it was accepted by the API. This does not necessarily mean that the processing was successful, because the API has no way of knowing right away if the filepaths you provided exist in your bucket. If there were errors in your filepaths, your HTTP request would get a 200, but your files would not actually process. They would turn up as empty files in the folder you provided. So if you see json files of around 70 bytes each in the destination folder, you know there was something wrong with your filenames.

You can check the status of a job with `dai_status()`. Just pass the response object from your `dai_async()` into the parentheses, and it will tell you whether the job is finished. It won't tell you how much time remains, but in my experience, processing takes about 5-20 seconds per page. 

```{r, eval=FALSE}
dai_status(response2)
```

When `dai_status()` says "SUCCEEDED", the json output files are waiting for you in the bucket. 

```{r, eval=FALSE}
gcs_list_objects()
```

Output file names look cryptic, but there's a logic to them, namely: `"<job_number>/<document_number>/<filename>-<shard_number>.json"`
Our file will thus take the form `"<job_number>/0/CIA-0.json"`, with `<job_number>` changing from one processing call to the next. Let us store the name in a vector for simplicity:

```{r}
## NOT RUN
our_file <- "<job_number>/0/CIA-0.json"
```

Now let's download it and save it under a simpler name:

```{r, eval=FALSE}
gcs_get_object(our_file, saveToDisk = "CIA.json", overwrite = TRUE)
```

Finally we extract the text using `text_from_dai_file`:

```{r, eval=FALSE}
text <- text_from_dai_file("CIA.json")
cat(text)
```

# Large batches

Although `dai_async()` takes batches of files, it is constrained by Google's [rate limits](https://cloud.google.com/document-ai/quotas). Currently, a `dai_async()` call can contain maximum 50 files (a multi-page pdf counts as one file), and you can not have more than 5 batch requests and 10 000 pages undergoing processing at any one time.

Therefore, if you're looking to process a large batch, you need to spread the `dai_async()` calls out over time. The simplest solution is to make a function that sends files off individually with a small wait in between. Say we have a vector called `big_batch` containing thousands of filenames. First we would make a function like this:

```{r, eval=FALSE}
process_slowly <- function(file) {
  dai_async(file)
  Sys.sleep(15)
}
```

Then we would iterate it over our file vector:

```{r, eval=FALSE}
## NOT RUN
map(big_batch, process_slowly)
```

This will hold up your console for a while, so it may be worth doing in the background as an RStudio [job](https://blog.rstudio.com/2019/03/14/rstudio-1-2-jobs/). 

Finding the optimal wait time for the `Sys-sleep()` may require some trial and error. As a rule of thumb, it should approximate the time it takes for DAI to process *one* of your files. This, in turn, depends on the size of the files, for a 100-page pdf will take a lot longer to process than a single-page one. In my experience, a 10-second interval works fine for a batch of single-page pdfs. Multi-page pdfs require proportionally more time. If your files vary in size, calibrate the wait time to the largest file, or you may get 429s (HTTP code for "rate limit exceeded") half way through the iteration.

Although this procedure is relatively slow, it need not add much to the overall processing time. DAI starts processing the first files it receives right away, so when your loop ends, DAI will be mostly done with the OCR as well. 

```{r, echo=FALSE, message=FALSE, warning=FALSE, eval=FALSE}
#cleanup
contents <- gcs_list_objects()
map(contents$name, gcs_delete_object)
```
---
title: "Correcting text output from Google Document AI"
author: "Thomas Hegghammer"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Correcting text output from Google Document AI}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
options(rmarkdown.html_vignette.check_title = FALSE)
library(knitr)
opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

**Last updated 4 April 2021**
\
\
\

Google Document AI (DAI) has excellent character recognition, but often reads columns wrong. This vignette will show you how to identify and reorder jumbled text with the tools in the `daiR` package.  

## The problem

Server-based OCR engines such as Google Document AI and Amazon Textract represent a major advance in OCR technology. They handle visual noise extremely well and effectively eliminate the need for image preprocessing, the most agonizing part of OCR in `tesseract` and other standalone libraries. DAI also reads non-Western languages such as Arabic better than any other general engine I have seen. 

But DAI and Textract still struggle with text columns and irregular page layouts. In my experience, DAI will misread a multi-column page about half the time, and the error rate increases with the complexity of the layout. This is not a problem if you plan to apply "bag-of-words" text mining techniques, but if you're looking at Natural Language Processing or actually reading the text, you cannot trust Document AI or Textract to always return accurate text. 

DAI column-reading errors are of two main types. The first is to put text blocks in the wrong order, and the second is to merge blocks that shouldn't be merged. Both errors can be corrected programmatically with the tools in the `daiR` package.    

## Reordering blocks

To illustrate, let's feed DAI a simple two-column text. This one is from the CIA's archive of declassified intelligence documents: 

```{r, message=FALSE, eval=FALSE}
setwd(tempdir())
download.file("https://www.cia.gov/readingroom/docs/1968-03-08.pdf", 
              destfile = "CIA_columns.pdf", 
              mode = "wb")
```

```{r, echo=FALSE, out.width = "50%"}
include_graphics("CIA_columns.jpg")
```

We start by uploading it to Google Storage and passing it to Document AI.  

```{r, echo=FALSE}
library(daiR)
suppressMessages(library(googleCloudStorageR))
```

```{r, eval=FALSE}
library(daiR)
library(googleCloudStorageR)
```

```{r, eval=FALSE}
resp <- gcs_upload("CIA_columns.pdf")
resp <- dai_async("CIA_columns.pdf")
```

Note that this code assumes you have set `GCS_AUTH_FILE` and `GCS_DEFAULT_BUCKET` variables in your .Renviron file (see the [two](https://dair.info/articles/setting_up_google_storage.html) [previous](https://dair.info/articles/using_document_ai.html) vignettes for details).

We check our bucket for the json output and download it when it's ready:

```{r, echo=FALSE, eval=FALSE}
content <- gcs_list_objects()
count <- 0
while (count < 150 && nrow(content) < 2){
  Sys.sleep(2)
  content <- gcs_list_objects()
  count <- count + 1
}
```

```{r, eval=FALSE}
gcs_list_objects()

## NOT RUN
our_file <- "<job_number>/0/CIA_columns-0.json"
```

```{r, eval=FALSE}
gcs_get_object(our_file, saveToDisk = "CIA_columns.json")
```


Finally we extract the text:

```{r, eval=FALSE}
text <- text_from_dai_file("CIA_columns.json")
cat(text)
```

On first inspection, this does not look so bad. But notice the transition from the first to the second paragraph: 

> ... they might reduce public support for the new Dubcek administration. Czech consumer in connection with these price increases.

Something's not right. Could it be a column-reading error?

We can find out with the function `draw_blocks()`, which extracts boundary box data from the `.json` file and draws numbered rectangles on an image of each page of the source document. (The images come from the json file, where they are stored as base64-encoded strings.) 

```{r, eval=FALSE}
draw_blocks(dest_path2, dir = tempdir())
```

Check your temporary directory (type `tempdir()` for the path) for a file ending in `_blocks.png` and pull it up:

```{r, echo=FALSE, out.width = "50%"}
include_graphics("CIA_columns1_blocks.png")
```

We can immediately see that the blocks are in the wrong order. How to fix this?

Fortunately, the `.json` file from DAI comes with a ton of data that allow us to programmatically reorder the text. The key is to generate a token dataframe with page location data and then filter and reorder as necessary. We create the dataframe with `build_token_df()`:

```{r, eval=FALSE}
token_df <- build_token_df(dest_path2)
str(token_df)
```

The dataframe has the words in the order in which DAI proposes to read them, and the `block` column has the number of the block to which each word belongs. This allows us to reorder the blocks while keeping the within-block word order intact. 

We see from the annotated image that the real order of the blocks should be 1 - 2 - 3 - 5 - 7 - 4 - 6. We can store this in a vector that we use to reorder the dataframe. 

```{r, eval=FALSE}
order <- c(1, 2, 3, 5, 7, 4, 6)
token_df$block <- factor(token_df$block, levels = order)
token_df_correct <- token_df[order(token_df$block),]
```

We get the correct text from the `token_df_correct$token` column:

```{r, message = FALSE, warning = FALSE, eval=FALSE}
library(dplyr)
text <- token_df_correct$token %>% 
  paste(collapse="")
```

Now the transition from the first to the second paragraph makes more sense:

```{r, eval=FALSE}
snippet <- substr(text, start = 1, stop = 700)
cat(snippet)
```

## Splitting blocks

A more complex --- and, unfortunately, more common --- situation is when DAI fails to distinguish between columns. This means that lines do not end where they should, resulting in long stretches of incomprehensible text. We can illustrate this with an article about the great [Peshtigo forest fire](https://en.wikipedia.org/wiki/Peshtigo_fire) in Wisconsin in 1871, available on the Internet Archive.

```{r, echo=FALSE, out.width = "50%"}
include_graphics("peshtigo.jpg")
```

We do our processing routine again:

```{r, eval=FALSE}
download.file("https://archive.org/download/themarinetteandpeshtigoeagleoct141871/The%20Marinette%20and%20Peshtigo%20Eagle%20-%20Oct%2014%201871.pdf", 
              destfile = "peshtigo.pdf", 
              mode = "wb")
resp <- gcs_upload("peshtigo.pdf")
resp <- dai_async("peshtigo.pdf")
```

```{r, eval=FALSE}
# wait till ready
gcs_list_objects()

## NOT RUN
our_file <- "<job_number>/0/peshtigo-0.json"
```

```{r, eval=FALSE}
gcs_get_object(our_file, saveToDisk = "peshtigo.json")
```


This time we'll skip the text printout and go straight to inspecting the boundary boxes:

```{r, eval=FALSE}
draw_blocks("peshtigo.json", dir = tempdir())
```

```{r, echo=FALSE, out.width = "50%"}
include_graphics("peshtigo1_blocks.png")
```

As we can see, this time DAI has failed to distinguish between the two main columns. We can verify this by checking the beginning of the text:

```{r, eval=FALSE}
text <- text_from_dai_file("peshtigo.json")
snippet <- substr(text, start = 1, stop = 1000)
cat(snippet)
```
This means that we must find a way of splitting block 12 vertically. 

What we will do is create a new boundary box that captures only the right-hand column. Then we will feed the location coordinates of the new box back into the token dataframe so that the tokens that fall within it are assigned a new block number. We can then reorder the blocks as we did in the previous example.

There are two main ways to obtain the coordinates of a new block: mathematically or through image annotation.

### Mathematical splitting

We can split blocks mathematically by using the location data for existing blocks in the json file. We start by building a block dataframe to keep track of the blocks. 

```{r, eval=FALSE}
block_df <- build_block_df(dest_path4)
```

Then we use the function `split_block()` to cut block 12 vertically in half. This function takes as input a block dataframe, the page and number of the block to split, and a parameter `cut_point`, which is a number from 1 to 99 for the relative location of the cut point. `split_block()` returns a new block dataframe that includes the new block and revised coordinates for the old one. 

```{r, eval=FALSE}
new_block_df <- split_block(block_df, block = 12, cut_point = 50)
```

If we had more blocks to split, we could repeat the procedure as many times as necessary. We just have to make sure to feed the latest version of the block dataframe into the `split_block()` function. 

When we have a block dataframe that captures the layout fairly accurately, we can use the `reassign_tokens()` function to assign new block values to the words in the *token* dataframe. `reassign_tokens()` takes as input the token dataframe and the new block dataframe and returns a revised token dataframe.

```{r, eval=FALSE}
token_df <- build_token_df(dest_path4)
token_df_correct <- reassign_tokens(token_df, new_block_df)  
```

In this particular case, the blocks are in the right order after splitting, so we can extract a correct text right away. In other cases the blocks may need reordering, in which case we use the procedure from the previous section.  

```{r, eval=FALSE}
text <- token_df_correct$token %>% 
  paste(collapse="")
snippet <- substr(text, start = 1, stop = 1000)
cat(snippet)
```

Mathematical splitting will often be the easiest method, and it can be particularly efficient when you have a lot of documents with the exact same column structure. However, it may sometimes be difficult to tell with the naked eye where the cut point should be. At other times the space between columns may be so narrow as to make precision important. For these situations we can use manual image annotation. 

### Manual splitting

In principle you can use any image annotation tool, so long as you format the resulting coordinates in a way that `daiR`'s processing functions understand. In the following, I will use [labelme](https://github.com/wkentaro/labelme) because it's easy to use and `daiR` has a helper function for it.

Labelme opens from the command line, but has a fairly intuitive graphical user interface. We load the annotated image generated by `draw_blocks()`, click "create polygons" in the left pane, right-click while the cursor is in the page pane, and choose "create rectangle". 

```{r, echo=FALSE, out.width = "50%"}
include_graphics("labelme1.png")
```

Then we mark the right-hand column and label it 13 (for the number of the new block).

```{r, echo=FALSE, out.width = "50%"}
include_graphics("labelme2.png")
```

Click "save" and store the json file, for example as `peshtigo1_blocks.json`. Now we can load it in R and get the coordinates of the new block 13 with the function `from_labelme()`. This function returns a one-row dataframe formatted like block dataframes generated with `build_block_df`.

```{r, eval=FALSE}
block13 <- from_labelme("peshtigo1_blocks.json")
```

We can then assign a new block number to the tokens that fall within block13. For this we use `reassign_tokens2()`, which reassigns tokens on a specified page according to the coordinates of a single new block. 

```{r, eval=FALSE}
token_df_new <- reassign_tokens2(token_df, block13) 
```

Now we just need to reorder the token data frame by blocks, and the words will be in the right order. In this particular case, we do not need to supply a custom block order, since the block numbering reflects the right order of the text.

```{r, eval=FALSE}
token_df_correct <- token_df_new[order(token_df_new$block), ]
```

And again we have a text in the right order. 

```{r, eval=FALSE}
text <- token_df_correct$token %>% 
  paste(collapse="")
snippet <- substr(text, start = 1, stop = 1000)
cat(snippet)
```

```{r, echo=FALSE, message=FALSE, warning=FALSE, eval=FALSE}
#cleanup
contents <- gcs_list_objects()
map(contents$name, gcs_delete_object)
```
---
title: "Complex file and folder management"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Complex file and folder management}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
options(rmarkdown.html_vignette.check_title = FALSE)

library(knitr)

opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

**Last updated 4 April 2021**
\
\
\

The main `daiR` vignettes use deliberately simple examples involving uploads of pdf files straight into the root of the bucket and down again. In real life you may be dealing with slightly more complex scenarios. 

## Image files

Document AI accepts only PDFs, GIFs and TIFFs, but sometimes your source documents are in other formats. `daiR`'s helper function `image_to_pdf()` is designed to help with this. Based as it is on `imagemagick`, it converts almost any image file format to pdf. You can also pass a vector of image files and ask for a single `.pdf` output, which is useful for collating several pagewise images to a single, multipage `.pdf`.

To illustrate, we can take this image of an old text from the National Park Service Website:

```{r, eval=FALSE}
setwd(tempdir())
download.file("https://www.nps.gov/articles/images/dec-of-sentiments-loc-copy.jpg", 
              destfile = "nps.jpg", 
              mode = "wb")
```

And convert it to a pdf like so:

```{r, eval=FALSE}
library(daiR)
image_to_pdf("nps.jpg", "nps.pdf")
```

And the file is ready for processing with Document AI.

## Processing a folder tree

At other times you may want to have folders inside your bucket. A typical scenario is when your source documents are stored in a folder tree and you want to batch process everything without losing the original folder structure. 

Problem is, it's technically not possible to have folders in Google Storage; files in a bucket are kept side by side in a flat structure. We can, however, *imitate* a folder structure by adding prefixes with forward slashes to filenames. This is not complicated, but requires paying attention to filenames at the upload and download stage. 

To illustrate, let's create two folders in our working directory: `folder1` and `folder2`: 

```{r, eval=FALSE}
library(fs)
dir1 <- file.path(tempdir(), "folder1")
dir2 <- file.path(tempdir(), "folder2")
dir_create(c(dir1, dir2))
```

Then we create three duplicates of the file `nps.pdf` and put two pdfs in each folder.  

```{r, eval=FALSE}
dest_path3 <- file.path(dir1, "nps.pdf")
dest_path4 <- file.path(dir1, "nps2.pdf")
dest_path5 <- file.path(dir2, "nps3.pdf")
dest_path6 <- file.path(dir2, "nps4.pdf")
file_copy(dest_path2, dest_path3)
file_copy(dest_path2, dest_path4)
file_copy(dest_path2, dest_path5)
file_copy(dest_path2, dest_path6)
```

To upload this entire structure to Google Storage, we create a vector of files in all subfolders with the parameter `recurse = TRUE` in the `dir_ls()` function. I'm assuming here that the working directory is otherwise empty of pdf files. 

```{r, eval=FALSE}
pdfs <- dir_ls(tempdir(), glob = "*.pdf", recurse = TRUE)
```

We then iterate the `gcs_upload()` function over our vector:

```{r, eval=FALSE}
library(googleCloudStorageR)
library(purrr)
resp <- map(pdfs, ~ gcs_upload(.x, name = .x))
```

If we now check the bucket contents, we see that the files are in their respective "folders".

```{r, eval=FALSE}
gcs_list_objects()
```

Bear in mind, though, that this is an optical illusion; the files are technically still on the same level. In reality, the `folder1/` and `folder2/` elements are an integral part of the filenames. 

We can process these files as they are with the following command:

```{r, eval=FALSE}
resp <- dai_async(pdfs) 
```

In which case DAI returns `.json` files titled `folder1/<job_number>/0/nps-0.json` and so forth. We can download these the usual way:

```{r, eval=FALSE}
content <- gcs_list_objects()
jsons <- grep("*.json$", content$name, value = TRUE)
resp <- map(jsons, ~ gcs_get_object(.x, saveToDisk = file.path(tempdir(), .x)))
```

And the json files will be stored in their respective subfolders alongside the source pdfs. 

Note, however, that this last script only worked because there already were folders titled `folder1` and `folder2` in our temporary directory. If there hadn't been, R would have returned an error, because the `gcs_get_object()` function cannot create new folders on your hard drive.

If you wanted to download the files to another folder where there wasn't a corresponding folder tree to "receive" them, you would have to use a workaround such as changing the forward slash in the bucket filepaths for an underscore (or something else) as follows:

```{r, eval=FALSE}
resp <- map(jsons, ~ gcs_get_object(.x, 
                                    saveToDisk = file.path(tempdir(), "folder3", gsub("/", "_", .x))
```

```{r, echo=FALSE, message=FALSE, warning=FALSE, eval=FALSE}
#cleanup
contents <- gcs_list_objects()
resp <- map(contents$name, gcs_delete_object)
```
---
title: "Setting up a Google Storage bucket"
author: "Thomas Hegghammer"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Setting up a Google Storage bucket}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
options(rmarkdown.html_vignette.check_title = FALSE)

library(knitr)

opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

**Last updated 4 April 2021**
\
\
\

This is a rundown for beginners of how to set up and interact with Google Storage in R. You need Google Storage in order to use Document AI and other Google APIs on scale, because these services do not accept bulk file submissions directly. Instead they use Google Storage as an intermediary, so you need to know how to get files in and out of it.    

It is possible to bulk upload and download files to Google Storage in the [Google Cloud Console](https://console.cloud.google.com/storage). In fact, for uploads it can sometimes be easier than doing it programmatically. But downloads and deletions are cumbersome if you have a lot of files. And since bulk processing in DAI can only be done with code, you might as well keep the whole workflow in R. 

## Authentication

The biggest hurdle to using any Google API is authentication, which is daunting for several reasons. It involves abstract new concepts like "service accounts", "Oauth2.0", and "scopes". Moreover, the Google Cloud Console is so crowded it's a nightmare to navigate as a beginner. In addition, different R packages have different procedures for authenticating with Google Cloud Services (GCS).

A full explanation of Google API authentication could fill a book, but suffice to say here that there are several different ways to authenticate to GCS from R. In the following I will walk you through one such way, the one I think is the simplest and most robust if you are primarily planning to use Google Storage and Google Document AI and to do so in an interactive setting. 

The procedure outlined below relies on service accounts and json key files. It involves a little bit more hands-on configuration than some of the setup wizards in other Google API wrappers, but the solution is robust and secure. My thinking is: if you are computer-literate enough to consider working programmatically with an API, you are also able to find your way around the Google Cloud Console. But know that there are [several](https://cran.r-project.org/package=gargle/vignettes/get-api-credentials.html) [other](https://cran.r-project.org/package=googleAuthR/vignettes/google-authentication-types.html) authentication strategies; which is better depends on where and how you intend to use the R wrapper.

### Step 1: Get a Gmail account

If you have one already, you can use that. Or you can [create](https://accounts.google.com/signup/v2/webcreateaccount) a burner account for your GCS work.

### Step 2: Activate the Google Cloud Console

While logged in to your gmail account, go to the [Google Cloud Console](https://console.cloud.google.com/). Agree to the terms of service and click "Try for free". 

```{r, echo=FALSE, out.width = "50%"}
include_graphics("storage1.PNG")
```

Accept the terms again, and add an address and a **credit card**. This last part is a prerequisite for using GCS.

### Step 3: Link your project to your billing account

The largest "unit" of your GCS activities is your *project*. You can think of it as your root folder, since you will most likely only ever need one unless you are a business or a developer (in principle, though, you can have as many projects as you like). 

When you activate GCS, you are assigned a project named "My first project". Click on "My first project" in the top blue bar, just to the right of "Google cloud services". You'll see a screen like this: 

```{r, echo=FALSE, out.width = "50%"}
include_graphics("storage2.PNG")
```

Note that your project has an *ID*, usually consisting of an adjective, a noun, and a number. You'll need this soon, so I recommend opening RStudio and storing it as a vector with `project_id <- "<your project id>"`. 

Return to the Google Cloud Console and look at the left column. Toward the top you see an entry called "Billing". Click it. You'll get to a screen saying "This project has no billing account". Click "link a billing account" and set the billing account to "My billing account". 

```{r, echo=FALSE, out.width = "50%"}
include_graphics("storage3.PNG")
```

All this is necessary for you to be able to access Google Storage and other Google tools programmatically. Both Google Storage and DAI are paid services, although for Google Storage the cost is negligible unless you plan to keep very large amounts of data there for a long time. For DAI, you're looking at around EUR 0.06 per processed page, though at the current time of writing, you get 300$ worth of free credits. 

### Step 4: Set up a service account

Now we need to create a service account. Bring out the navigation menu on the left hand side by clicking the little circle with the three horizontal lines in the top left of the screen. Click on "APIs and services." Then click on "credentials" in the left pane. You should see this: 

```{r, echo=FALSE, out.width = "50%"}
include_graphics("storage4.PNG")
```

Click on "create credentials" in the top middle, then choose service account. Give it any name you like (e.g. "my_rstudio_service_account") and a description (e.g. "Interacting with GCS through R") and click "create". 

In section 2 titled "Grant this service account access to project", add "Basic > **Owner**" to the service account's roles.  

```{r, echo=FALSE, out.width = "50%"}
include_graphics("storage5.PNG")
```

Click "continue", then "done" at the bottom. You should now see your service account listed at the bottom. 

```{r, echo=FALSE, out.width = "50%"}
include_graphics("storage6.PNG")
```

### Step 5: Download a json file with the service account key

Now we need to generate a json file containing the login details for this service account. Click the small edit icon on the bottom right. On the next page, click "add key", choose "create new key", select JSON format, and click "create". This should prompt a save file window. Save the file to your hard drive. You can change the name to something more memorable if you like (but keep the ".json" extension). Also, take note of where you stored it. Now we are done in the Google Cloud Console and can finally start working in RStudio.

### Step 6: Tell RStudio where to find the json file

The last step is to store the path to the json file in your .Renviron file so that RStudio can authenticate you whenever you are working with GCS from R. Start by writing the following in the console:

```{r, message=FALSE, eval=FALSE}
usethis::edit_r_environ()
```

This will open a pane with your .Renviron file. If you haven't modified it before, it is probably empty. 

All you need to do is add a line with the following: **`GCS_AUTH_FILE='<full path to the json file you stored earlier>'`**. Make sure all the slashes in the filepath are forward slashes. Save the file, close it, and restart RStudio.

Now when you load the library `googleCloudStorageR`, you will be auto-authenticated and ready to communicate with your Google Storage account from within R. 

```{r, message=FALSE}
library(googleCloudStorageR)
```

## Working with googleCloudStorageR

`googleCloudStorageR` is a so-called wrapper for the Google Storage API, which means it translates your R input into URLs that the API can understand. When you execute `googleCloudStorageR` functions, you are really sending GET and POST requests to Google and receiving responses in return.

### Creating and inspecting buckets

Google Storage is a file repository, and it keeps your files in so-called "buckets". You need at least one bucket to store files. To inspect your Storage account, first bring out your project id. If you did not store it in step 3 above, you can look it up in the Google Cloud Console or use the `daiR` function `get_project_id()`. 

```{r, eval=FALSE}
project_id <- daiR::get_project_id()
```

Now let's see how many buckets we have:
```{r, eval = FALSE}
gcs_list_buckets(project_id)
```

Answer: zero, because we haven't created one yet. We can do this with `gcs_create_bucket()`. Note that it has to be globally unique ("my_bucket" won't work because someone's already taken it). For this example, let's use "dair-example-bucket". Also add a location ("EU" or "US").

```{r, eval = FALSE}
gcs_create_bucket("dair-example-bucket", project_id, location = "EU")
```

Now we can see the bucket listed:
```{r, eval=FALSE}
gcs_list_buckets(project_id)
```

You will need to supply a bucket name with every call to Google Storage (and Document AI), so you may want to store the name of a default bucket in the environment. You have two options:

1) Store it permanently in you .Renviron file by calling `usethis::edit_r_environ()` and adding `GCS_DEFAULT_BUCKET=<your bucket name>` to the list of variables, just as you did with the json key file path (`GCS_AUTH_FILE`) earlier. Note that adding a default bucket to .Renviron will not prevent you from supplying other bucket names in individual calls to Google Storage when necessary. 

2) Set it for the current session with `gcs_global_bucket("<your bucket name>")`


To get a bucket's file inventory, we use `gcs_list_objects()`. Leaving the parentheses empty will get information about the default bucket if you have set it.
```{r, eval=FALSE}
gcs_list_objects()
```

At this point it's obviously empty, so let's upload something. 

### Uploading files

This we do with `gcs_upload()`. If the file is in your working directory, just write the filename; otherwise provide the full file path. If you want, you can store the file under another name in Google Storage with the `name` parameter. Otherwise, just leave the parameter out.  

```{r, eval=FALSE}
setwd(tempdir())
write.csv(mtcars, "mtcars.csv")
resp <- gcs_upload("mtcars.csv")
```

`resp` here is just a container for the response that the API returns, to prevent it from printing in the console. 

Now let's check the contents:
```{r, eval=FALSE}
contents <- gcs_list_objects()
contents
```

The Google Storage API handles only one file at a time, so for bulk uploads you need to use iteration.  

```{r, eval=FALSE}
library(purrr)
library(fs)
write.csv(iris, "iris.csv")
my_files <- dir_ls(glob = "*.csv")
resp <- map(my_files, ~ gcs_upload(.x, name = .x))
```

Let's check the contents again:
```{r, eval=FALSE}
contents <- gcs_list_objects()
contents
```

Note that there's a file size limit of 5Mb, but you can change it with `gcs_upload_set_limit()`. 
```{r, eval=FALSE}
gcs_upload_set_limit(upload_limit = 20000000L)
```

### Downloading files

Downloads are performed with gcs_get_object(). Here, too, you can save it under a different name, but this time the parameter is `saveToDisk`.

```{r, eval=FALSE}
resp <- gcs_get_object("mtcars.csv", saveToDisk = "mtcars_duplicate.csv")
```

To download multiple files we need to loop or map. Let's say we wanted to download all the files in the bucket:

```{r, eval=FALSE}
contents <- gcs_list_objects()
resp <- map(contents$name, ~ gcs_get_object(.x, saveToDisk = .x, overwrite = TRUE))
```

### Deleting

We can delete files in the bucket with `gcs_delete_object()`: 
```{r, eval=FALSE}
gcs_delete_object("mtcars.csv")
```

To delete several, we again need to loop or map. Let's try to delete everything in the bucket:

```{r, eval=FALSE}
contents <- gcs_list_objects()
map(contents$name, gcs_delete_object)
```

And the bucket is empty. With this you should be ready to process files with Document AI and other APIs that involve Google Storage as an intermediary.

```{r, echo=FALSE, message=FALSE, warning=FALSE, eval=FALSE}
# cleanup
file_delete(dir_ls(tempdir(), glob = "*.csv"))
```
---
title: "Basic processing"
output: html_document
vignette: >
  %\VignetteIndexEntry{Basic processing}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

### Process synchronously

Pass a single-page pdf or image file to Document AI and get the output immediately:

```{r, eval=FALSE}
library(daiR) 
## Not run:
myfile <- "<sample.pdf>"
response <- dai_sync(myfile)
text <- text_from_dai_response(response)
cat(text)
```

### Process asynchronously

Send larger batches for offline processing in three steps:

#### 1. Upload files to your Google Storage bucket

```{r, eval=FALSE}
## Not run:
library(googleCloudStorageR)

my_pdfs <- c("<sample1.pdf>", "<sample2.pdf>")
purrr::map(my_pdfs, ~ gcs_upload(.x, name = .x))
```

#### 2. Tell Document AI to process them:

```{r, eval=FALSE}
## Not run:
dai_async(my_pdfs)
```

#### 3. Download the json output and extract the text:

```{r, eval=FALSE}
## Not run:
bucket_contents <- gcs_list_objects()
only_jsons <- grep("*.json", bucket_contents$name, value = TRUE)
map(only_jsons, ~ gcs_get_object(.x, saveToDisk = .x))
text <- text_from_dai_file(only_jsons[1])
cat(text)
```
---
title: "Extracting tables"
author: "Thomas Hegghammer"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Extracting tables}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
options(rmarkdown.html_vignette.check_title = FALSE)
library(knitr)
opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

**Last updated 4 April 2021**
\
\
\

In optical character recognition (OCR) of historical documents, we're often interested in extracting tables and getting the data into R in tidy format. Having been designed primarily for processing business documents, Google Document AI has very powerful table extraction capabilities. With the help of the `daiR` package, you can get these tables into R dataframes with a single line of code.

To extract table data, we use the `dai_sync_tab()` and `dai_async_tab()` commands. They work much like their siblings `dai_sync()`/`dai_async()` except that their responses also include table data. In return, the `*sync_tab()` functions have weaker language support (they currently only work for English-language documents).^[Under the hood, the `*sync()` and `*sync_tab()` functions access different endpoints in the Document AI API, each with its own capabilities. I expect these capabilities to be consolidated in a single endpoint in the future, in which case the `*sync_tab()` functions will be deprecated.] 

To get the table data into R, we use the functions `tables_from_dai_response()` and `tables_from_dai_file()`. They do the same thing, but take different inputs. You use the former for response objects obtained with `dai_sync_tab()` and the latter for .json files obtained with `dai_async_tab()`.

## Tables from response objects

Let's try sync-processing a small pdf from the so-called [Truth Tobacco Industry Documents](https://en.wikipedia.org/wiki/Truth_Tobacco_Industry_Documents). It contains four tables like this: 

```{r, echo=FALSE, out.width = "100%"}
include_graphics("tobacco.PNG")
```

First we download it and send it to Document AI for processing:
```{r, eval=FALSE}
library(daiR)
setwd(tempdir())
download.file("https://archive.org/download/tobacco_lpnn0000/lpnn0000.pdf", 
              destfile = "tobacco.pdf",
              mode = "wb")

resp <- dai_sync_tab("tobacco.pdf")
```

To get the tables into R, we pass the response object (`resp`) to our extraction function:

```{r, eval=FALSE}
tables <- tables_from_dai_response(resp)
```

Note that the output is not a dataframe but a *list* of dataframes. We access the individual tables either by their index, like so:

```{r, eval=FALSE}
table1 <- tables[[1]]
```

... or by importing them all into our global environment with `assign()`, as follows:

```{r, eval=FALSE}
for(i in 1:length(tables)) {
  assign(paste0("table", i), tables[[i]])
}
```

Either way, the result is pretty good. Here's what the first table looks like in my RStudio before any cleaning:

```{r, echo=FALSE, out.width = "100%"}
include_graphics("table.PNG")
```

## Tables from json files

With json files obtained from `dai_async_tab()` the process is very similar, except that we pass a filepath instead of a response object to the extraction function. If we processed the same document asynchronously and obtained a file named `"tobacco.json"`, we would extract the tables like this:

```{r, eval=FALSE}
## NOT RUN
tables <- tables_from_dai_file("tobacco.json")
```

That's all there is to it.

Note that the quality of the output depends on the structure of the original table. Like all other table extraction engines, DAI struggles with complex tables that contain row and column spans. In such cases you can usually turn DAI's initial output into a decent dataframe with a bit of data cleaning inside R. In most cases, though, Document AI delivers very decent output that `daiR` can turn into workable R dataframes in a matter of seconds. 
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_output.R
\name{build_block_df}
\alias{build_block_df}
\title{Build block dataframe}
\usage{
build_block_df(json)
}
\arguments{
\item{json}{filepath of a JSON file obtained using \code{dai_async()}}
}
\value{
a block data frame
}
\description{
Creates a dataframe with the block bounding boxes
identified by Document AI (DAI) in an asynchronous request.
Rows are blocks, in the order DAI proposes to read them. Columns
are location variables such as page coordinates and page numbers.
}
\details{
The location variables are: page number, left boundary,
right boundary, top boundary, and bottom boundary.
}
\examples{
\dontrun{
block_df <- build_block_df("pdf_output.json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_output.R
\name{reassign_tokens2}
\alias{reassign_tokens2}
\title{Assign tokens to a single new block}
\usage{
reassign_tokens2(token_df, block, page = 1)
}
\arguments{
\item{token_df}{a data frame generated by \code{dair::build_token_df}}

\item{block}{a one-row data frame of the same format as \code{token_df}}

\item{page}{the number of the page on which the block belongs}
}
\value{
a token data frame
}
\description{
This is a specialized function for use in connection
with text reordering. It is designed to facilitate manual splitting
of block boundary boxes and typically takes a one-row block dataframe
generated by \code{from_labelme()}.
}
\examples{
\dontrun{
new_token_df <- reassign_tokens2(token_df, new_block_df)
new_token_df <- reassign_tokens2(token_df, new_block_df, 5)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{pdf_to_binbase}
\alias{pdf_to_binbase}
\title{PDF to base64 tiff}
\usage{
pdf_to_binbase(file)
}
\arguments{
\item{file}{path to a single-page pdf file}
}
\value{
a base64-encoded string
}
\description{
Converts a PDF file to a base64-encoded binary .tiff file.
}
\examples{
\dontrun{
doc_encoded <- pdf_to_binbase("document.pdf")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inspect_output.R
\name{draw_blocks}
\alias{draw_blocks}
\title{Inspect block bounding boxes}
\usage{
draw_blocks(json, prefix = NULL, dir = getwd())
}
\arguments{
\item{json}{filepath of a JSON file obtained using \code{dai_async()}}

\item{prefix}{string to be prepended to output filename}

\item{dir}{path to the desired output directory}
}
\value{
no return value, called for side effects
}
\description{
Plots the block bounding boxes identified by
Document AI (DAI) onto images of the submitted document.
Generates an annotated .png file for each page in the original
document.
}
\details{
Not vectorized, but documents can be multi-page.
}
\examples{
\dontrun{
draw_blocks("pdf_output.json", dir = tempdir())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/send_to_dai.R
\name{dai_async}
\alias{dai_async}
\title{OCR documents asynchronously}
\usage{
dai_async(
  files,
  dest_folder = NULL,
  bucket = Sys.getenv("GCS_DEFAULT_BUCKET"),
  proj_id = get_project_id(),
  proc_id = Sys.getenv("DAI_PROCESSOR_ID"),
  skip_rev = "true",
  loc = "eu",
  token = dai_token()
)
}
\arguments{
\item{files}{a vector or list of pdf filepaths in a GCS Storage bucket
Filepaths must include all parent bucket folder(s) except the bucket name}

\item{dest_folder}{the name of the GCS Storage bucket subfolder where
you want the json output}

\item{bucket}{the name of the GCS Storage bucket where the files
to be processed are located}

\item{proj_id}{a GCS project id}

\item{proc_id}{a Document AI processor id}

\item{skip_rev}{whether to skip human review; "true" or "false"}

\item{loc}{a two-letter region code; "eu" or "us"}

\item{token}{an access token generated by \code{dai_auth()} or another
auth function}
}
\value{
A list of HTTP responses
}
\description{
Sends files from a Google Cloud Services (GCS) Storage
bucket to the GCS Document AI v1 API for asynchronous (offline) processing.
The output is delivered to the same bucket as JSON files containing
the OCRed text and additional data.
}
\details{
Requires a GCS access token and some configuration of the
.Renviron file; see package vignettes for details. Currently, a
\code{dai_async()} call can contain a maximum of 50 files (but a
multi-page pdf counts as one file). You can not have more than
5 batch requests and 10,000 pages undergoing processing at any one time.
Maximum pdf document length is 2,000 pages. With long pdf documents,
Document AI divides the JSON output into separate files ('shards') of
20 pages each. If you want longer shards, use \code{dai_tab_async()},
which accesses another API endpoint that allows for shards of up to
100 pages.
}
\examples{
\dontrun{
# with daiR configured on your system, several parameters are automatically provided,
# and you can pass simple calls, such as:
dai_async("my_document.pdf")

# NB: Include all parent bucket folders (but not the bucket name) in the filepath:
dai_async("for_processing/pdfs/my_document.pdf")

# Bulk process by passing a vector of filepaths in the files argument:
dai_async(my_files)

# Specify a bucket subfolder for the json output:
dai_async(my_files, dest_folder = "processed")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_output.R
\name{split_block}
\alias{split_block}
\title{Split a block bounding box}
\usage{
split_block(block_df, page = 1, block, cut_point, direction = "v")
}
\arguments{
\item{block_df}{A dataframe generated by \code{build_block_df()}.}

\item{page}{The number of the page where the split will be made.
Defaults to 1.}

\item{block}{The number of the block to be split.}

\item{cut_point}{A number between 0 and 100, where 0 is the
existing left/top limit and 100 is the existing right/bottom limit.}

\item{direction}{"V" for vertical split or "H" for horizontal split.
Defaults to "V".}
}
\value{
a block data frame
}
\description{
This function 'splits' (in the sense of changing the
coordinates) of an existing block bounding box vertically or
horizontally at a specified point. It takes a block data frame as
input and modifies it. The splitting produces a new block, which
is added to the data frame while the old block's coordinates are
updated. The function returns a revised block data frame.
}
\examples{
\dontrun{
new_block_df <- split_block(df = old_block_df, block = 7, cut_point = 33)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_tables.R
\name{tables_from_dai_response}
\alias{tables_from_dai_response}
\title{Get tables from response object}
\usage{
tables_from_dai_response(object)
}
\arguments{
\item{object}{an HTTP response object returned by \code{dai_sync_tab()}}
}
\value{
a list of data frames
}
\description{
Extracts all the tables that Document AI (DAI)
identified in a synchronous processing request.
}
\examples{
\dontrun{
tables <- tables_from_dai_response(response)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auth.R
\name{dai_token}
\alias{dai_token}
\title{Produce access token}
\usage{
dai_token(
  path = Sys.getenv("GCS_AUTH_FILE"),
  scopes = "https://www.googleapis.com/auth/cloud-platform"
)
}
\arguments{
\item{path}{path to a JSON file with a service account key}

\item{scopes}{GCS auth scopes for the token}
}
\value{
a GCS access token object (if credentials are valid) or a message (if not).
}
\description{
Produces an access token for Google Cloud Services (GCS)
}
\examples{
\dontrun{
token <- dai_token()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inspect_output.R
\name{draw_paragraphs}
\alias{draw_paragraphs}
\title{Inspect paragraph bounding boxes}
\usage{
draw_paragraphs(json, prefix = NULL, dir = getwd())
}
\arguments{
\item{json}{filepath of a JSON file obtained using \code{dai_async()}}

\item{prefix}{string to be prepended to output filename}

\item{dir}{path to the desired output directory.}
}
\value{
no return value, called for side effects
}
\description{
Plots the paragraph bounding boxes identified by
Document AI (DAI) onto images of the submitted document.
Generates an annotated .png file for each page in the original
document.
}
\details{
Not vectorized, but documents can be multi-page.
}
\examples{
\dontrun{
draw_paragraphs("pdf_output.json", dir = tempdir())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/send_to_dai.R
\name{dai_sync}
\alias{dai_sync}
\title{OCR document synchronously}
\usage{
dai_sync(
  file,
  proj_id = get_project_id(),
  proc_id = Sys.getenv("DAI_PROCESSOR_ID"),
  skip_rev = "true",
  loc = "eu",
  token = dai_token()
)
}
\arguments{
\item{file}{path to a single-page pdf or image file}

\item{proj_id}{a GCS project id.}

\item{proc_id}{a Document AI processor id}

\item{skip_rev}{whether to skip human review; "true" or "false".}

\item{loc}{a two-letter region code; "eu" or "us".}

\item{token}{an authentication token generated by \code{dai_auth()} or
another auth function.}
}
\value{
a HTTP response object
}
\description{
Sends a single document to the Google Cloud Services (GCS)
Document AI v1 API for synchronous (immediate) processing. Returns a
HTTP response object containing the OCRed text and additional data.
}
\details{
Requires a GCS access token and some configuration of the
.Renviron file; see package vignettes for details.Input files can be in
either .pdf, .bmp, .gif, .jpeg, .jpg, .png, or .tiff format. PDF files
can be up to five pages long. Extract the text from the response object with
\code{text_from_dai_response()}. Inspect the entire response object with
\code{httr::content()}.
}
\examples{
\dontrun{
response <- dai_sync("doc_page.pdf")

my_page_scan <- "001.png"
response <- dai_sync(my_page_scan)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_output.R
\name{reassign_tokens}
\alias{reassign_tokens}
\title{Assign tokens to new blocks}
\usage{
reassign_tokens(token_df, block_df)
}
\arguments{
\item{token_df}{a dataframe generated by \code{build_token_df()}}

\item{block_df}{a dataframe generated by \code{dair::split_block()}
or \code{dair::build_block_df()}}
}
\value{
a token data frame
}
\description{
This is a specialized function for use in connection
with text reordering. It modifies a token dataframe by assigning
new block bounding box values to a subset of tokens based on
prior modifications made to a block dataframe.
}
\details{
The token and block data frames provided as input must be
from the same JSON output file.
}
\examples{
\dontrun{
new_token_df <- reassign_tokens(token_df, new_block_df)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_tables.R
\name{tables_from_dai_file}
\alias{tables_from_dai_file}
\title{Get tables from output file}
\usage{
tables_from_dai_file(file)
}
\arguments{
\item{file}{filepath of a JSON file obtained using \code{dai_async_tab()}}
}
\value{
a list of data frames
}
\description{
Extracts all the tables that Document AI (DAI)
identified in an asynchronous processing request.
}
\examples{
\dontrun{
tables <- tables_from_dai_file("document.json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inspect_output.R
\name{text_from_dai_response}
\alias{text_from_dai_response}
\title{Get text from HTTP response object}
\usage{
text_from_dai_response(
  response,
  save_to_file = FALSE,
  dest_dir = getwd(),
  filename = "output"
)
}
\arguments{
\item{response}{an HTTP response object returned by \code{dai_sync()}}

\item{save_to_file}{boolean; whether to save the text as a .txt file}

\item{dest_dir}{folder path for the .txt output file if \code{save_to_file = TRUE}}

\item{filename}{string to form the stem of the .txt output file}
}
\value{
a string (if \code{save_to_file = FALSE})
}
\description{
Extracts the text OCRed by Document AI (DAI) in a
synchronous processing request.
}
\examples{
\dontrun{
text <- text_from_dai_response(response)

text_from_dai_response(response, save_to_file = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inspect_output.R
\name{text_from_dai_file}
\alias{text_from_dai_file}
\title{Get text from output file}
\usage{
text_from_dai_file(file, save_to_file = FALSE, dest_dir = getwd())
}
\arguments{
\item{file}{filepath of a JSON file obtained using \code{dai_async()}}

\item{save_to_file}{boolean; whether to save the text as a .txt file}

\item{dest_dir}{folder path for the .txt output file if save_to_file=TRUE}
}
\value{
a string (if \code{save_to_file = FALSE})
}
\description{
Extracts the text OCRed by Document AI (DAI) in an
asynchronous processing request.
}
\examples{
\dontrun{
text <- text_from_dai_file("mydoc-0.json")
text_from_dai_file("mydoc-0.json", save_to_file = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/send_to_dai.R
\name{dai_async_tab}
\alias{dai_async_tab}
\title{OCR asynchronously and get table data}
\usage{
dai_async_tab(
  files,
  filetype = "pdf",
  dest_folder = NULL,
  bucket = Sys.getenv("GCS_DEFAULT_BUCKET"),
  proj_id = get_project_id(),
  loc = "eu",
  token = dai_token(),
  pps = 100
)
}
\arguments{
\item{files}{A vector or list of pdf filepaths in a GCS Storage bucket.
Filepaths must include all parent bucket folder(s) except the bucket name.}

\item{filetype}{Either "pdf", "gif", or "tiff". If \code{files} is a
vector, all elements must be of the same type.}

\item{dest_folder}{The name of the bucket subfolder where you want the
JSON output.}

\item{bucket}{The name of the GCS Storage bucket. Not necessary if
you have set a default bucket as a .Renviron variable named
\code{GCS_DEFAULT_BUCKET} as described in the package vignette}

\item{proj_id}{a GCS project id}

\item{loc}{a two-letter region code ("eu" or "us")}

\item{token}{an access token generated by \code{dai_auth()} or another
auth function.}

\item{pps}{an integer from 1 to 100 for the desired number of pages per
shard in the JSON output}
}
\value{
A list of HTTP responses
}
\description{
Sends files from a Google Cloud Services (GCS) Storage
bucket to the GCS Document AI v1beta2 API for asynchronous (offline)
processing. The output is delivered to the same bucket as JSON files
containing the OCRed text and additional information, including
table-related data.
}
\details{
This function accesses a different API endpoint than the main
\code{dai_async()} function, one that has less language support, but
returns table data in addition to parsed text (which \code{dai_async()}
currently does not). This function may be deprecated if/when the v1
API endpoint incorporates table extraction. Use of this service
requires a GCS access token and some configuration of the .Renviron file;
see vignettes for details. Note that this API endpoint does not require
a Document AI processor id. Maximum pdf document length is 2,000 pages,
and the maximum number of pages in active processing is 10,000. Also note
that this function does not provide 'true' batch processing; instead it
successively submits single requests at 10-second intervals.
}
\examples{
\dontrun{
# with daiR configured on your system, several parameters are automatically provided,
# and you can pass simple calls, such as:
dai_async_tab("my_document.pdf")

# NB: Include all parent bucket folders (but not the bucket name) in the filepath:
dai_async_tab("for_processing/pdfs/my_document.pdf")

# Bulk process by passing a vector of filepaths in the files argument:
dai_async_tab(my_files)

# Specify a bucket subfolder for the json output:
dai_async_tab(my_files, dest_folder = "processed")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/send_to_dai.R
\name{dai_notify}
\alias{dai_notify}
\title{Notify on job completion}
\usage{
dai_notify(response, loc = "eu", token = dai_token(), sound = 2)
}
\arguments{
\item{response}{a HTTP response object generated by \code{dai_async()}}

\item{loc}{A two-letter region code; "eu" or "us"}

\item{token}{An authentication token generated by \code{dai_auth()} or
another auth function}

\item{sound}{A number from 1 to 10 for the Beepr sound selection
(https://www.r-project.org/nosvn/pandoc/beepr.html).}
}
\value{
no return value, called for side effects
}
\description{
Queries to the Google Cloud Services (GCS) Document AI API
about the status of a previously submitted asynchronous job
and emits a sound notification when the job is complete.
}
\examples{
\dontrun{
response <- dai_async(myfiles)
dai_notify(response)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{is_pdf}
\alias{is_pdf}
\title{Check that a file is PDF}
\usage{
is_pdf(file)
}
\arguments{
\item{file}{a filepath}
}
\value{
a boolean
}
\description{
Checks whether a file is a PDF file.
}
\examples{
\dontrun{
is_pdf("document.pdf")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/send_to_dai.R
\name{dai_sync_tab}
\alias{dai_sync_tab}
\title{OCR synchronously and get table data}
\usage{
dai_sync_tab(file, proj_id = get_project_id(), loc = "eu", token = dai_token())
}
\arguments{
\item{file}{path to a single pdf or image file}

\item{proj_id}{a GCS project id}

\item{loc}{a two-letter region code ("eu" or "us")}

\item{token}{An access token generated by \code{dai_auth()} or another
auth function.}
}
\value{
a HTTP response object
}
\description{
Sends a single document to the Google Cloud Services (GCS)
Document AI v1beta2 API for synchronous (immediate) processing. Returns
a response object containing the OCRed text and additional information,
including table-related data.
}
\details{
This function accesses a different API endpoint than the main
\code{dai_sync()} function, one that has less language support, but
returns table data in addition to parsed text (which \code{dai_sync()}
currently does not). This function may be deprecated if/when the v1
endpoint incorporates table extraction. Use of this service requires
a GCS access token and some configuration of the .Renviron file; see
vignettes for details. Input files can be in either .pdf, .bmp, .gif,
.jpeg, .jpg, .png, or .tiff format. PDFs can be up to five pages long.
Extract the text from the response object with
\code{text_from_dai_response()}. Inspect the entire response object
with \code{httr::content()}.
}
\examples{
\dontrun{
response <- dai_sync("doc_page.pdf")

my_page_scan <- "001.png"
response <- dai_sync(my_page_scan)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inspect_output.R
\name{merge_shards}
\alias{merge_shards}
\title{Merge shards}
\usage{
merge_shards(source_dir, dest_dir)
}
\arguments{
\item{source_dir}{folder path for input files}

\item{dest_dir}{folder path for output files}
}
\value{
no return value, called for side effects
}
\description{
Merges text files from Document AI output shards into a
single text file corresponding to the parent document.
}
\details{
The function works on .txt files generated from .json output files,
not on .json files directly. It also presupposes that the .txt filenames
have the same name stems as the .json files from which they were extracted.
For the v1 API, this means files ending with "-0.txt", "-1.txt", "-2.txt",
and so forth. For the v1beta2 API, it means files ending with
"-page-1-to-100.txt", "-page-101-to-200.txt", etc. The safest approach is
to generate .txt files using \code{text_from_dai_file()} with the \code{save_to_file}
parameter set to TRUE.
}
\examples{
\dontrun{
merge_shards(source_dir = getwd(), dest_dir = tempdir())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auth.R
\name{dai_auth}
\alias{dai_auth}
\title{Check authentication}
\usage{
dai_auth(
  path = Sys.getenv("GCS_AUTH_FILE"),
  scopes = "https://www.googleapis.com/auth/cloud-platform"
)
}
\arguments{
\item{path}{path to a JSON file with a service account key}

\item{scopes}{GCS auth scopes for the token}
}
\value{
no return value, called for side effects
}
\description{
Checks whether the user can obtain an access token for
Google Cloud Services (GCS) using a service account key stored on file.
}
\details{
daiR takes a very parsimonious approach to authentication,
with the native auth functions only supporting service account files.
Those who prefer other authentication methods can pass those directly
to the \code{token} parameter in the various functions that call the
Document AI API.
}
\examples{
\dontrun{
dai_auth()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/send_to_dai.R
\name{dai_status}
\alias{dai_status}
\title{Check job status}
\usage{
dai_status(response, loc = "eu", token = dai_token(), verbose = FALSE)
}
\arguments{
\item{response}{A HTTP response object generated by \code{dai_async()} or
\code{dai_tab_async()}}

\item{loc}{A two-letter region code; "eu" or "us"}

\item{token}{An authentication token generated by \code{dai_auth()} or
another auth function}

\item{verbose}{boolean; Whether to output the full response}
}
\value{
If verbose was set to \code{TRUE}, a HTTP response object.
If verbose was set to \code{FALSE}, a string summarizing the status.
}
\description{
Queries the Google Cloud Services (GCS) Document AI API
about the status of a previously submitted asynchronous job.
}
\examples{
\dontrun{
# Short status message:
response <- dai_async(myfiles)
dai_status(response)

# Full status details:
response <- dai_async(myfiles)
status <- dai_status(response, verbose = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_output.R
\name{build_token_df}
\alias{build_token_df}
\title{Build token dataframe}
\usage{
build_token_df(json)
}
\arguments{
\item{json}{filepath of a JSON file obtained using \code{dai_async()}}
}
\value{
a token data frame
}
\description{
Builds a token dataframe from the text OCRed by
Document AI (DAI) in an asynchronous request. Rows are tokens, in the
order DAI proposes to read them. Columns are location variables
such as page coordinates and block bounding box numbers.
}
\details{
The location variables are: start index, end index,
left boundary, right boundary, top boundary, bottom boundary,
page number, and block number. Start and end indices refer to
character position in the string containing the full text.
}
\examples{
\dontrun{
token_df <- build_token_df("pdf_output.json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inspect_output.R
\name{draw_tokens}
\alias{draw_tokens}
\title{Inspect token bounding boxes}
\usage{
draw_tokens(json, prefix = NULL, dir = getwd())
}
\arguments{
\item{json}{filepath of a JSON file obtained using \code{dai_async()}}

\item{prefix}{string to be prepended to output filename}

\item{dir}{path to the desired output directory.}
}
\value{
no return value, called for side effects
}
\description{
Plots the token (i.e., word) bounding boxes identified
by Document AI (DAI) onto images of the submitted document.
Generates an annotated .png file for each page in the original
document.
}
\details{
Not vectorized, but documents can be multi-page.
}
\examples{
\dontrun{
draw_tokens("pdf_output.json", dir = tempdir())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_output.R
\name{redraw_blocks}
\alias{redraw_blocks}
\title{Inspect revised block bounding boxes}
\usage{
redraw_blocks(json, token_df, dir = getwd())
}
\arguments{
\item{json}{filepath of a JSON file obtained using \code{dai_async()}}

\item{token_df}{a token data frame generated with \code{build_token_df()},
\code{reassign_tokens()}, or \code{reassign_tokens2()}.}

\item{dir}{path to the desired output directory.}
}
\value{
no return value, called for side effects
}
\description{
Tool to visually check the order of block bounding boxes after
manual processing (e.g. block reordering or splitting). Takes as its main
input a token dataframe generated with \code{build_token_df()},
\code{reassign_tokens()}, or \code{reassign_tokens2()}.
The function plots the block bounding boxes onto images of the submitted
document. Generates an annotated .png file for each page in the
original document.
}
\details{
Not vectorized, but documents can be multi-page.
}
\examples{
\dontrun{
redraw_blocks("pdf_output.json", revised_token_df, dir = tempdir())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{image_to_pdf}
\alias{image_to_pdf}
\title{Convert images to PDF}
\usage{
image_to_pdf(files, pdf_name)
}
\arguments{
\item{files}{a vector of image files}

\item{pdf_name}{a string with the name of the new PDF}
}
\value{
no return value, called for side effects
}
\description{
This helper function converts a vector of images to a
single PDF.
}
\details{
Combines any number of image files of almost any type
to a single PDF. The vector can consist of different image file types.
See the 'Magick' package documentation \url{https://cran.r-project.org/package=magick}
for details on supported file types. Note that on Linux, ImageMagick may
not allow conversion to pdf for security reasons.
}
\examples{
\dontrun{
# Single file
new_pdf <- file.path(tempdir(), "document.pdf")
image_to_pdf("document.jpg", new_pdf)

# A vector of image files:
image_to_pdf(images)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{.onAttach}
\alias{.onAttach}
\title{Run when daiR is attached}
\usage{
.onAttach(libname, pkgname)
}
\arguments{
\item{libname}{name of library}

\item{pkgname}{name of package}
}
\description{
Run when daiR is attached
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{is_json}
\alias{is_json}
\title{Check that a file is JSON}
\usage{
is_json(file)
}
\arguments{
\item{file}{a filepath}
}
\value{
a boolean
}
\description{
Checks whether a file is a JSON file.
}
\examples{
\dontrun{
is_json("file.json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_output.R
\name{from_labelme}
\alias{from_labelme}
\title{Extract block coordinates from labelme files}
\usage{
from_labelme(json, page = 1)
}
\arguments{
\item{json}{a json file generated by 'Labelme'}

\item{page}{the number of the annotated page}
}
\value{
a data frame with location coordinates for the rectangle
marked in 'Labelme'.
}
\description{
This is a specialized function for use in connection
with text reordering. It takes the output from the image
annotation tool 'Labelme' \url{https://github.com/wkentaro/labelme}
and turns it into a one-row data frame compatible with other
'daiR' functions for text reordering such as
\code{reassign_tokens2()}. See package vignette on text reconstruction
for details.
}
\examples{
\dontrun{
new_block <- from_labelme("document1_blocks.json")
new_block <- from_labelme("document5_blocks.json", 5)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auth.R
\name{dai_user}
\alias{dai_user}
\title{Get user information}
\usage{
dai_user()
}
\value{
a list of user information elements
}
\description{
Fetches the Google Cloud Services (GCS) user information
associated with a service account key.
}
\examples{
\dontrun{
dai_user()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inspect_output.R
\name{draw_lines}
\alias{draw_lines}
\title{Inspect line bounding boxes}
\usage{
draw_lines(json, prefix = NULL, dir = getwd())
}
\arguments{
\item{json}{filepath of a JSON file obtained using \code{dai_async()}}

\item{prefix}{string to be prepended to output filename}

\item{dir}{path to the desired output directory.}
}
\value{
no return value, called for side effects
}
\description{
Plots the line bounding boxes identified by
Document AI (DAI) onto images of the submitted document.
Generates an annotated .png file for each page in the original
document.
}
\details{
Not vectorized, but documents can be multi-page.
}
\examples{
\dontrun{
draw_lines("pdf_output.json", dir = tempdir())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auth.R
\name{get_project_id}
\alias{get_project_id}
\title{Get project id}
\usage{
get_project_id(path = Sys.getenv("GCS_AUTH_FILE"))
}
\arguments{
\item{path}{path to the JSON file with your service account key}
}
\value{
a string with a GCS project id
}
\description{
Fetches the Google Cloud Services (GCS) project id
associated with a service account key.
}
\examples{
\dontrun{
project_id <- get_project_id()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{img_to_binbase}
\alias{img_to_binbase}
\title{Image to base64 tiff}
\usage{
img_to_binbase(file)
}
\arguments{
\item{file}{path to an image file}
}
\value{
a base64-encoded string
}
\description{
Converts an image file to a base64-encoded binary .tiff file.
}
\examples{
\dontrun{
img_encoded <- pdf_to_binbase("image.png")
}
}
