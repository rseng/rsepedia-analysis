# phonfieldwork

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![](https://badges.ropensci.org/385_status.svg)](https://github.com/ropensci/software-review/issues/385)
[![CRAN version](http://www.r-pkg.org/badges/version/phonfieldwork)](https://cran.r-project.org/package=phonfieldwork)
[![](http://cranlogs.r-pkg.org/badges/grand-total/phonfieldwork)](https://CRAN.R-project.org/package=phonfieldwork)
[![R build status](https://github.com/ropensci/phonfieldwork/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/phonfieldwork/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/ropensci/phonfieldwork/master.svg)](https://codecov.io/github/ropensci/phonfieldwork?branch=master)
[![DOI](https://zenodo.org/badge/194053227.svg)](https://zenodo.org/badge/latestdoi/194053227)

`phonfieldwork` is a package for phonetic fieldwork research and experiments. This package make it easier:

- creating a html/pptx presentation from stimuli-translation list, 
- renaming soundfiles according to the list of stimuli, 
- concatenating multiple soundfiles and create a Praat TextGrid whose interval labels are the original names of the sound
- extracting sounds according to annotation
- extracting annotation from multiple linguistic formats (Praat `.TextGrid`, ELAN `.eaf`, EXMARaLDA `.exb`, Audacity `.txt` and subtitles `.srt`)
- visualising an oscilogram, a spectrogram and an annotation
- creating an html viewer [like this](https://ropensci.github.io/phonfieldwork/s1/stimuli_viewer.html), ethical problems of this kind of viewer in linguistic research are covered in the vignette `vignette("ethical_research_with_phonfieldwork")`.

For more ditails see [tutorial](https://docs.ropensci.org/phonfieldwork/).

The main goal of the `phonfieldwork` package is to make a full research workflow from data collection to data extraction and data representation easier for people that are not familiar with programming. Hovewer most of the `phonfieldwork` funnctionality can be found in other software and packages:

* stimuli presentation creation could be done with any programming language and probably without them
* automatic file renaming and automatic merge could be done with any programming language
* Praat `.TextGrid` manipulation could be done with Praat, R packages [`rPraat`](https://cran.r-project.org/package=rPraat) and [`textgRid`](https://cran.r-project.org/package=textgRid), the Python package ['pympi'](https://dopefishh.github.io/pympi/index.html))
* ELAN `.eaf` manipulationcould be done with ELAN, the R package [`FRelan`](https://github.com/langdoc/FRelan) and the Python package [`pympi`](https://dopefishh.github.io/pympi/index.html)
* Praat `.TextGrid`, ELAN `.eaf`, and 'EXMARaLDA .exb import and export could be done with the R package [`act`](https://cran.r-project.org/package=act)
* cut sounds according to annotation could be done with Praat and the R package`tuneR`
* spectrogram visualisation could be done with multiple R packages [`signal`](https://cran.r-project.org/package=signal), [`tuneR`](https://cran.r-project.org/package=tuneR), [`seewave`](https://cran.r-project.org/package=seewave), [`phonTools`](https://cran.r-project.org/package=phonTools), [`monitor`](https://cran.r-project.org/package=monitor), [`warbleR`](https://cran.r-project.org/package=warbleR), [`soundgen`](https://cran.r-project.org/package=soundgen) and many others

## Installation

Install from CRAN:

```
install.packages("phonfieldwork")
```

Get the development version from GitHub:

```
install.packages("remotes")
remotes::install_github("ropensci/phonfieldwork")
```
Load a library:
```
library(phonfieldwork)
```

In order to work with some `rmarkdown` functions you will need to install `pandoc`, see `vignette("pandoc")` for the details.

## To do:

* export to ELAN and EXMARALDA files
* use ELAN and EXMARALDA files in the whole pipline discribed in docs
* use the same pipline with video (for Sign Languages)
* make [TECkit](https://scripts.sil.org/cms/scripts/render_download.php?format=file&media_id=BeyondUTR22_pdf&filename=BeyondUTR22_pdf.pdf) to df and back

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# phonfieldwork 0.0.12

- add a `separate_duration` argument to the `concatenate_soundfiles()` function that makes it possible to use some silent separator during the file concatenation.
- make `rename_soundfiles()` function to work with mp3 files.
- fix the bug with `"` sign in textgrids.
- make pictures optional in the `create_viewer()` function

# phonfieldwork 0.0.11

- correct empty tiers behavior #34 (thanks to Shungo Suzuki)
- add possibility to have different values in the `n_of_annotations` argument of `create_subannotation()` (thanks to Jenya Korovina for the idea)
- rename `tier` argument of the `create_empty_textgrid()` to `tier_name`.
- create the `remove_textgrid_tier()` function.

# phonfieldwork 0.0.10

- add tryCatch to the `read_from_folder()` function

# phonfieldwork 0.0.8

- add possibility to read short format of `.TextGrid`s and fix `textgrid_to_df()` and `tier_to_df()` functions
- add `encoding`, `formant_df`, `intensity`, `picth` and `pitch_range` arguments to the `draw_sound()` function
- add the `formant_to_df()` function
- add the `picth_to_df()` function
- add the `intensity_to_df()` function
- add an argument `external` to the `create_presentation()` function in order to mark external images or gifs
- remove all `encoding` arguments and replace it with encoding autodetection from `uchardet` (thanks to Artem Klevtsov for help)
- add `autonumber`, `loging` and `missing` arguments to `rename_soundfiles()` function (thanks to Niko Partanen)
- a lot of minor style changes (thanks to Jonathan Keane)
- add the `create_empty_textgrid()` function (thanks to Niko Partanen)
- add the `data_manipulation_with_tidyverse` vignette (thanks to Niko Partanen)
- add the `concatenate_textgrids()` function
- add the `read_from_folder()` function and remove `..._from_folder` arguments
- pass rOpenSci review! Move tutorial to <https://ropensci.github.io/phonfieldwork/>

# phonfieldwork 0.0.7

- add a vigniettes about ethical research and introduction to work with phonfieldwork
- add an argument `textgrids_from_folder` to the `textgrid_to_df()` function
- add an argument `exbs_from_folder` to the `exb_to_df()` function
- add an argument `eafs_from_folder` to the `eaf_to_df()` function
- add Raven style annotations8
- replace freqmax with frequency_range argument
- replace example_textgrid with systemfile() call
- add window annotation to spectrograms
- add bridge to lingtypology package: `map` argument in the `create_viewer()` function
- make it possible to visualise all types of annotations with the `draw_sound()` function
- add the `source` column to all `..._to_df()` functions
- add the `audacity_to_df` function
- add the `srt_to_df` function
- chage textgrid related functions' output from `start`, `end`, `annotation` to `time_start`, `time_end`, `content`
- correct point tier visualization

# phonfieldwork 0.0.6

- add `encoding` arguments to functions for working with TextGrids
- add `textgrid` argument to the `draw_sound()` function
- change subgraphs alignment in the `draw_sound()` function including textgrid annotation
- add `from` and `to` arguments to the `draw_sound()` function
- add .mp3 format reading options to all functions that work with sounds
- add `text_size` argument to the `draw_sound()` function
- add `zoom`` argument to the `draw_sound()` function
- add the `get_sound_duration()` function
- fix ploting of multiple sounds with multiple .TextGrids
- add an argument `title_as_filename` to the `draw_sound()` function
- change .TextGrid associated arguments of the `create_viewer()` function to `table` argument; as a result users now need to provide a table for the annotation viewer and not a .TextGrid

# phonfieldwork 0.0.5

- add `textgrid_to_df()` function for reading Praat files
- add `create_glossed_document()` function for converting .flextext files into a glossed document
- add `flextext_to_df()` function for reading FLEx files
- add `eaf_to_df()` function for reading ELAN files
- add `exb_to_df()` function for reading EXMARaLDA files
- rename `textgrid` argument into `annotation` argument in `concatenate_soundfiles()` function adding new possible values

# phonfieldwork 0.0.4

- add `create_subannotation()` function

# phonfieldwork 0.0.3

- vertically and horisontally center text in presentations created by `create_presentation()`; thx @Pandaklez #1
- add the `font_size` argument to the `create_presentation()` function
- add `rename_videofiles()` function
- rebuild html viewer for sounds with JavaScript with the help of new functions `create_image_look_up()` and `create_sound_play()`.

# phonfieldwork 0.0.2

- make the `create_presentation()` function render silently
- add a new function `draw_sound()` for creating spectrogram and oscilogram
- add a new function `create_viewer()` for creating an html viewer with sound and spectrograms
- correct work of `autonumbering` function in `extract_intervals()` function
- add new functions  `get_textgrid_names()` and `set_textgrid_names()`
- finish tutorial <https://ropensci.github.io/phonfieldwork/>

# phonfieldwork 0.0.1

- initial release
# Contributing to `phonfieldwork`

## Issues

When filing an issue, the most important thing is to include a minimal reproducible example so that I can quickly verify the problem. So please include:

* required packages
* package versions
* data
* code

There are some additional information that I put in the issue template.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''
---

Please, put the results of these functions below, if your issue is related to the technical bugs. Thank you!

<details> <summary> info about OS and package versions </summary>
```
sessionInfo()$R.version$platform
sessionInfo()$R.version$version.string
packageVersion("rmarkdown")
packageVersion("phonfieldork")
```
</details>
---
name: Feature request
about: 'Suggest an idea for this project '
title: ''
labels: enhancement
assignees: ''

---


