---
title: "Bioassays: a new package in R for analyzing cellular assays"

tags:
  - R
  - Multi-well plate
  - bioassays

authors:
  - name: Anwar Azad Palakkan
    orcid: 0000-0002-1693-1219
    affiliation: 1

affiliations:
  - name: Deanery of Biomedical Sciences, University of Edinburgh, UK
    index: 1

date: 11 June 2020

bibliography: paper.bib
---

# Summary

  Experimental techniques that acquire data from multi-well plates are very common in life science research; today they account for half of all assays used in drug development and toxicity testing. These experiments include colorimetric, fluorometric and luminometric assays for cytotoxicity / cell viability [@Markossian:2004], bacterial assays [@Burton:2007; @Pant:2016] and immunoassays [@Lin:2015; @Song:2015]. Such assays feature in small-scale experiments and also in high-throughput screens. High-throughput imaging and spectrofluorometric quantification-based methods are increasingly used in the pharmaceutical industry and in academic research aiming to screen large libraries (e.g. potential drug compounds, RNAi, CRISPR/Cas9) or to quantify rare biological events.

Plate-reading equipment is usually supplied with accompanying analytic software, such as  Gen5 (BioTek), Multi-Mode Analysis Software (Molecular Devices), and Multiwell-Analyzer (Multi Channel Systems). There are, however, drawbacks to its use. Usually, the software is bound to the data coming from a specific piece of equipment, it is seldom possible to use the software on multiple computers without making additional license purchases. This is a problem in academic settings in which a single machine is shared and used intensively by several research groups. Also, proprietary software has closed code that cannot be checked or understood in detail by the  user and cannot be expanded or adapted for analyses not envisaged by the equipment manufacturer.

In response to these limitations, some open-source software packages have been developed to serve particular unusual applications such as single-cell migration [@Masuzzo:2017] or differential scanning fluorimetry [@Wang:2012]. Open software with more general capabilities like DRfit [@Hofmann:2019] that run on JAVA platforms has been developed, but they are not ideal for high throughput screening. Many biologists are familiar with the statistical language, R, and R libraries have been created to assist with analysis of data from plate readers. Examples include ‘platetools’ (https://cran.r-project.org/package=platetools), ‘plater’ (https://cran.r-project.org/package=plater) and ‘phenoScreen’ (https://rpubs.com/Swarchal/phenoScreen), each of which addresses particular functions. ‘Plater’, for example, is focused on data formatting, ‘phenoScreen’ on data visualization, and ‘platetools’ on both data visualizing and formatting. The packages are useful, but can be difficult for a beginner to grasp, partly because their documentation seems not to have been written with beginners in mind.

Here, we introduce a new open-source package, 'bioassays', designed to provide a wide range of functions relevant to multi-well plate assays. ‘Bioassays’ is an R package freely available on CRAN [@R:2014], and has supports for formatting, visualizing and analyzing multi-well plate assays. It has functions for handling outliers, for handling multiple data sets with separate blanks, for estimating values from standard curves, for summarizing data and for doing statistical analysis. Moreover, we have provided examples (bioassays-examples) in the vignettes, which will be very helpful for beginners to grasp how this package can be used. Because this is a R package, more experienced users will be able add to its functionality using other R packages.


# Bioassays
Bioassays can handle data from any of the standard multi-well plate format: 6, 12, 24, 96 or 384 well plates. A prerequisite for the package is the need for both input data (\autoref{figure 1}) and metadata (\autoref{figure 2}) to be in comma-separated variable (csv) format; most plate readers can export data as a .csv file. Metadata need  “row” and “col” columns to indicate the location of well.  The package has functions for data extraction (extract_filename), formatting (data2plateformat, plate2df, matrix96 and plate_metadata), visualization (heatplate) and data analysis (reduceblank, estimate, dfsummary, pvalue).

The function extract_filename is useful for extracting specific information from file names, such as compound name, plate number etc., which can be used for automated data analysis. This function provides a very easy way to pass data of this type into the analysis by simply editing file names.

The function data2plateformat converts the data (eg: readings from a 96 well plate) to an appropriate matrix format with suitable row and column names. The plate2df function can convert such data into a data frame with ‘row’ and ‘col’ columns to indicate well position. The function matrix96 can convert both character and numeric columns to a matrix. This function can automatically determine the type of multi-well from which the data is coming, and it have provision to convert negative and NA values to 0, if needed. The plate_metadata function combines plate-specific information (such as compound used, standard concentration, dilution of samples, etc.), to produce unique plate metadata.

For data visualization, the ‘heatplate’ function can be used. This can plot both heat plots (\autoref{figure 3}) and categorical plots (\autoref{figure 4}) automatically, depending on the data type.  The ‘heatplate’ function also has provisions for adjusting the displayed well size, for visual esthetics. Any multi-well plate data can be structured for the plots (function 'heatplate') by using functions 'data2plateformat', 'plate2df', 'matrix96' alone or in combinations (depending on the input data).

The function ‘reduceblank’ can help to reduce blank values from the readings. This function can handle separate blanks for different datasets. The function ‘estimate’ can be used to estimate unknown variable (eg: concentration) based on standard curve. The function ‘dfsummary’ is really versatile for handling multiple data sets. It can group samples and summarize data sets separately. It has additional controls for handling outliers and omitting unwanted data sets. Function ‘pvalue’ can be used to test significance by t-test. It has provisions for asserting control group and level of significance.

# Availability and installation

Bioassays is supported on Windows and macOS. The package can be installed using install.packages command in R. The source code, vignette, datasets and detailed examples on how to use the package are available on CRAN ( https://CRAN.R-project.org/package=bioassays) and GitHub (https://github.com/anwarbio/bioassays).

# Acknowledgment
The author acknowledges Prof Jamie Davies for suggestions and correcting the article. Support from lab members (Jamie Davies Group, Centre for Discovery Brain Sciences, University of Ediburgh) are acknowledged.

# Figures
![Input data format from a 96 well plate reading.\label{figure 1}](figure1.png)

![Metafile data format.\label{figure 2}](figure2.png)

![Heat map of 384 well (normalized values).\label{figure 3}](figure3.png)

![Categorical plot of 384 well plate.\label{figure 4}](figure4.png)

# References
# Contributing

Thanks for your interest in contributing to bioassays.
Here are some guidelines to help make it easier to merge your Pull Request:

* For potentially large changes, please open an Issue first to discuss
* Please follow the [Hadley style guide][style]
* Run `devtools::test()` to run the tests
* (Optional) Add new test(s) in `tests/testthat/`

If you're new to submitting Pull Requests, please read the section [Contribute
to other projects][contribute] in the tutorial [A quick introduction to version
control with Git and GitHub][git-tutorial].

## More about this repository

For the most part, I tried to follow the guidelines from [R packages][r-pkg] by
[Hadley Wickham][hadley]. The unit tests are performed with [testthat][], the
documentation is built with [roxygen2][], and the online package documentation
is created with [rmarkdown](https://github.com/rstudio/rmarkdown). Continuous integration testing is performed for
for macOS by [Travis CI](https://travis-ci.com/github/anwarbio/bioassays/builds), and for Windows
by [AppVeyor](https://ci.appveyor.com/project/anwarbio/bioassays).

The repository contains the files `LICENSE` and `LICENSE.md` to both adhere to
[R package conventions for defining the license][r-exts-licensing] and also to
make the license clear in a more conventional manner (suggestions for
improvement welcome). Directories are standard for R packages
as described in the manual [Writing R Extensions][r-exts].

## Release checklist (for maintainers)

* Bump version with [scripts/bump-version.R](scripts/bump-version.R)
* Update [NEWS.md](NEWS.md): Check `git log` and make sure to reference GitHub
Issues/PRs
* Run [scripts/build.sh](scripts/build.sh) to confirm tests pass locally
* Test on [rhub][]:
    * Have to validate email first with `rhub::validate_email()`. Copy-paste
    token from email into R console.
    * Check on Ubuntu with `rhub::check_on_ubuntu()`
    * Check on macOS with `rhub::check_on_macos()`
    * Check on Windows with `rhub::check_on_windows()`
* Test on [winbuilder][]:
    * Check with R release with `devtools::check_win_release()`
    * Check with R devel with `devtools::check_win_devel()`
* Update [cran-comments.md](cran-comments.md)
* Commit with `git commit -am "Bump version: x.x.x.9xxx -> x.x.x and re-build
docs."`
* Push with `git push origin master` and wait for CI builds to pass
* Tag with `git tag -a vx.x.x`. Summarize [NEWS.md](NEWS.md) entry into bullet
points. Run ` git tag -l -n9` for past examples. Push with `git push origin
--tags`.
* Make a release. On GitHub, go to Releases -> Tags -> Edit release notes. Name
the release "workflowr x.x.x" and copy-paste the Markdown entry from
[NEWS.md](NEWS.md).
* Build tarball with `R CMD build .` and upload to [CRAN submission
site][cran-submit]. You will receive an email to request confirmation, then an
email confirming the package was submitted, and then an email with the test
results. Once it is accepted to CRAN, monitor the [check results][check-results]
for any surprise errors. Also, these builds are when the binaries are built for
Windows and macOS, so they aren't available until they are finished. You will
receive an email once all the Windows binaries are available for download
(devel, release, oldrel).

[appveyor]: https://ci.appveyor.com
[check-results]: https://cran.r-project.org/web/checks/check_results_workflowr.html
[circleci]: https://circleci.com
[Codecov]: https://codecov.io/
[contribute]: http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004668#sec011
[covr]: https://github.com/jimhester/covr
[cran-submit]: https://cran.r-project.org/submit.html
[foghorn]: https://cran.r-project.org/package=foghorn
[git-tutorial]: http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004668
[hadley]: http://hadley.nz/
[pkgdown]: https://github.com/r-lib/pkgdown
[pt]: https://rstudio.github.io/rstudio-extensions/rstudio_project_templates.html
[r-exts]: https://cran.r-project.org/doc/manuals/R-exts.html
[r-exts-licensing]: https://cran.r-project.org/doc/manuals/R-exts.html#Licensing
[r-pkg]: http://r-pkgs.had.co.nz/
[rhub]: https://r-hub.github.io/rhub/
[roxygen2]: https://github.com/klutometis/roxygen
[style]: http://adv-r.had.co.nz/Style.html
[testthat]: https://github.com/hadley/testthat
[travis]: https://travis-ci.org/
[winbuilder]: https://win-builder.r-project.org/

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Bioassays

<!-- badges: start -->

<!-- badges: end -->

'Bioassays' is an R package, designed to provide a wide range of functions relevant to multi-well plate assays. It can handle data from any of the standard multi-well plate format: 6, 12, 24, 96 or 384 well plates. ‘Bioassays’ can help in formatting, visualizing and analyzing multi-well plate data. It has functions for handling outliers, for handling multiple data sets with separate blanks, for estimating values from standard curves, for summarizing data and for doing statistical analysis. Moreover, it is strongly documented in a manner designed to be easy for even beginners to grasp.
## Installation

You can install the released version of bioassays from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("bioassays")
```
The latest development version can be installed from github with the remotes package
``` r
library(remotes)
install_github("anwarbio/bioassays")
```
## Example

Detail’s of various funtions in this package is provided in
‘bioassays-vignette’. Examples on how to use this package is provided
in ‘bioassays-example’ in Vignette.

## License: GPL-3.0
https://github.com/anwarbio/bioassays/blob/master/LICENSE

## Contribute

We love your input! Users may request new features by opening a [GitHub Issue](https://github.com/anwarbio/bioassays/issues), or may contribute their own additions and improvements via a pull request. Similarly, if you run into problems while using this package, or require technical support, do not hesitate to request support through a [GitHub Issue](https://github.com/anwarbio/bioassays/issues). If you use ‘Bioassays’ in your work and would like to further collaborate, I would be more than willing to discuss it over email or [GitHub Issue](https://github.com/anwarbio/bioassays/issues). When you submit code changes, your submissions are understood to be under the same [GPL-3.0](https://github.com/anwarbio/bioassays/blob/master/LICENSE) that covers this project.

An incomplete list of possible improvements:
* Include support for qPCR data analysis
* Include provisions for user template
* Include capability for plotting data

## version 1.0.0
- OCT 2020
---

### NEWS.md setup

- added NEWS.md creation with newsmd

