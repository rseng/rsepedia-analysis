## FielDHub: A Shiny App for Design of Experiments in Life Sciences


[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)

## Installation

### Development version from GitHub

``` r
devtools::install_github("DidierMurilloF/FielDHub")
```

### Stable version from R CRAN

``` r
install.packages("FielDHub")
```

## FielDHub Paper

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03122/status.svg)](https://doi.org/10.21105/joss.03122)

## Overview

A shiny design of experiments (DOE) app that aids in the creation of traditional, un-replicated, augmented and partially-replicated designs applied to agriculture, plant breeding, forestry, animal and biological sciences. 

For more details and examples of all functions present in the FielDHub package. Please, go to <https://didiermurillof.github.io/FielDHub/reference/index.html>.


![](FielDHub_Overview.png)

## Usage

This is a basic example which shows you how to launch the app:

``` r
library(FielDHub)
run_app()
```
### Diagonal Arrangement Example

A project needs to test 270 genotypes in a field containing 20 rows and 15 columns of plots. In this example, these 270 genotypes are divided among three different experiments. In addition, four checks are included in a systematic diagonal arrangement across experiments to fill 30 plots representing 10% of the total number of experimental plots. An option to include filler plots is also available for fields where the number of experimental plots does not equal the number of available field plots.

![](DExample.PNG)

The figure above shows a map of an experiment randomized as a Decision Block Unreplicated Design with Checks on Diagonals. Yellow, gray, and green shade the blocks of unreplicated experiments, while distinctively colored check plots are replicated throughout the field in a systematic diagonal arrangement.

To illustrate using FielDHub to build experimental designs through R code, the design produced in the R Shiny interface described above can also be created using the function `diagonal_arrangement()` in the R script below. Note, that to obtain identical results, users must include the same seed number in the script as was used in the Shiny app. In this case, the seed number is 1249. 

``` r
diagonalExample <- diagonal_arrangement(nrows = 20, 
                                        ncols = 15, 
                                        lines = 270, 
                                        checks = 4, 
                                        plotNumber = 101, 
                                        splitBy = "row", 
                                        seed = 1249, 
                                        kindExpt = "DBUDC", 
                                        blocks = c(100, 100, 70))
```

Users can access the returned values from `diagonal_arrangement()` as follow,

``` r
> diagonalExample$infoDesign
$Lines
[1] 100 100  70

$checks
[1] 4

$RepChecks
1 2 3 4 
7 7 8 8 

$percentChecks
[1] "10%"

$Fillers
[1] 0

$seed
[1] 1249

> diagonalExample$layoutRandom
      Col1 Col2 Col3 Col4 Col5 Col6 Col7 Col8 Col9 Col10 Col11 Col12 Col13 Col14 Col15
Row20  248    2  268  225  258  223  255  226  264   239   274     1   232   235   256
Row19  214  273  244  220    3  269  247  207  243   215   245   221   208   229     4
Row18  240  252  251  272  238  212  250    3  206   242   259   253   222   213   227
Row17    2  261  210  230  254  262  266  237  236   231     4   218   219   257   265
Row16  224  263  211    1  216  209  217  270  271   205   233   260   234     2   246
Row15  166  135  201  127  194  203    1  171  190   116   156   228   241   249   267
Row14  126  191  193  162  204  137  114  178  173     1   115   130   167   155   199
Row13  111  143    4  182  150  154  117  160  136   192   122   119     3   177   131
Row12  105  197  159  157  181    1  147  107  153   158   170   152   141   148   187
Row11  161  186  112  175  183  124  200  195    2   129   120   144   146   188   118
Row10  128    3  138  165  142  184  196  121  169   132   168     4   113   172   133
Row9   145  198  139  179    2  109  180  140  125   134   149   164   176   106     3
Row8   202  189  174  163  108  185  110    4  151   123    84    70    43    57   101
Row7     4   24   15   93   22   78   87   96   18    47     3    44   100    99    53
Row6    92   19    8    1   62   49   36   85   38     5    59    33     7     2    39
Row5    66   82  104   26   34   14    4   97   73    16    65    41    77   103    42
Row4    54   86   46   75   13   69   98   10   63     2    37    71    29    45    40
Row3    89   28    3   88   17   20   72   31   35    25    56    74     3    68    55
Row2    91   23   11   51   30    1   83   64   50    90    79    21    12     6    95
Row1    48  102    9   94   52   80   76   60    4    61    32    58    67    27    81

> head(diagonalExample$fieldBook, 12)
   ID  EXPT LOCATION YEAR PLOT ROW COLUMN CHECKS ENTRY TREATMENT
1   1 Expt1        1 2021  101   1      1      0    48    gen 48
2   2 Expt1        1 2021  102   1      2      0   102   gen 102
3   3 Expt1        1 2021  103   1      3      0     9     gen 9
4   4 Expt1        1 2021  104   1      4      0    94    gen 94
5   5 Expt1        1 2021  105   1      5      0    52    gen 52
6   6 Expt1        1 2021  106   1      6      0    80    gen 80
7   7 Expt1        1 2021  107   1      7      0    76    gen 76
8   8 Expt1        1 2021  108   1      8      0    60    gen 60
9   9 Expt1        1 2021  109   1      9      4     4   Check 4
10 10 Expt1        1 2021  110   1     10      0    61    gen 61
11 11 Expt1        1 2021  111   1     11      0    32    gen 32
12 12 Expt1        1 2021  112   1     12      0    58    gen 58
```

The main difference between using the FielDHub Shiny app and using the standalone function `diagonal_arrangement()` is that the standalone function will allocate filler only if it is necessary, while in R Shiny, filler plots are generated automatically. In cases where users include fillers, either between or after experiments, the Shiny app is preferable for filling and visualizing all field plots.

### Partially Replicated Design Example

Partially replicated designs are commonly employed in early generation field trials. This type of design is characterized by replication of a portion of the entries, with the remaining entries only appearing once in the experiment. As an example, considered a field trial with 288 plots containing 75 entries appearing two times each, and 138 entries only appearing once. This field trials is arranged in a field of 16 rows by 18 columns.

![](pREPExample.PNG)

In the figure above, green plots contain replicated entries, and yellow plots contain entries that only appear once.

Instead of using the Shiny FielDHub app, users can use the standalone FielDHub function `partially_replicated()`. The partially replicated layout described above can be produced through scripting as follows. As noted in the previous example, to obtain identical results between the script and the Shiny app, users need to use the same seed number, which, in this case, is 77.

``` r
pREPexample <- partially_replicated(nrows = 16, 
                                    ncols = 18,  
                                    repGens = c(138,75),
                                    repUnits = c(1,2),
                                    planter = "serpentine", 
                                    plotNumber = 1,
                                    exptName = "ExptA",
                                    locationNames = "FARGO",
                                    seed = 77)

```
Users can access returned values from `partially_replicated()` as follows,

``` r
> pREPexample$layoutRandom
      Col1 Col2 Col3 Col4 Col5 Col6 Col7 Col8 Col9 Col10 Col11 Col12 Col13 Col14 Col15 Col16 Col17 Col18
Row16  176   45  192  127   71   62  102  120   95    79    61   158   141    32    21    25    14    71
Row15   18  121    9  123  171   60   11  204    2   182   119    65   198    63   188    34   124    55
Row14   68  108  200   16    3    5  136  183   78    43   149    89    23    47   160   168   152    10
Row13  187   66  153   24   21   23   65  199   34    64    11   191    94    96   131   181   186   103
Row12  211   15    8  174   56  169   28   84   90    46   156    57    51    42   128   170    27   114
Row11   19    7   87    7   93   86   26   22  151    67    30     5    92    48    51   150    59    41
Row10  144  101  104  134  213  167   73   33    6    74    25   122    32   107    18    40   139    59
Row9    70   38   70   12   88   99  179  210   37    22   155   142   116   161   137    10    19    30
Row8     1  173   57   20  157   69   13  110  195   147    38   163     1    83    40   162     4   209
Row7    97   82   36   29   56   35   50    6   75    85    58    17    66   172   145    63    44     3
Row6    24    8   43  177   31  143   91   64   98   194    74    55   126   115   117   159   100    45
Row5    77  205   37   62  106   29  138   53  118    35   132    54   166    81    47    41   109    76
Row4    73  129  125  135   60  164   14   39   13   165   112   105   140    12    72   212    46   203
Row3    69   28  196  154  146   75  189   36   54    72    31    33    53    67   175    27     2   130
Row2    16   68  206   49   52  208   26  184  180   193   178    61    48   197    52    50   201   111
Row1    80   49   15   44  185    9  133   58    4   113   190   148    39    17    20   207   202    42

> head(pREPexample$fieldBook, 12)
   ID  EXPT LOCATION YEAR PLOT ROW COLUMN CHECKS ENTRY TREATMENT
1   1 ExptA    FARGO 2021    1   1      1      0    80       G80
2   2 ExptA    FARGO 2021    2   1      2      1    49       G49
4   3 ExptA    FARGO 2021    3   1      3      1    15       G15
6   4 ExptA    FARGO 2021    4   1      4      1    44       G44
8   5 ExptA    FARGO 2021    5   1      5      0   185      G185
10  6 ExptA    FARGO 2021    6   1      6      1     9        G9
11  7 ExptA    FARGO 2021    7   1      7      0   133      G133
13  8 ExptA    FARGO 2021    8   1      8      1    58       G58
15  9 ExptA    FARGO 2021    9   1      9      1     4        G4
16 10 ExptA    FARGO 2021   10   1     10      0   113      G113
17 11 ExptA    FARGO 2021   11   1     11      0   190      G190
18 12 ExptA    FARGO 2021   12   1     12      0   148      G148
```

To see more examples, please go to <https://didiermurillof.github.io/FielDHub/reference/index.html>.


# FielDHub 0.1.0

* Added a `NEWS.md` file to track changes to the package.
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

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
* Focusing on what is best not just for us as individuals, but for the overall
community

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

Community leaders are responsible for clarifying and enforcing our standards
of acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies
when an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at [email][didier.murilloflorez@gmail.com]. All complaints will be reviewed and investigated promptly and fairly.

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

**Community Impact**: A violation through a single incident or series of
actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or permanent
ban.

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
standards, including sustained inappropriate behavior, harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within the
community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at <https://www.contributor-covenant.org/version/2/0/code_of_conduct.html>.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
<https://www.contributor-covenant.org/faq>. Translations are available at <https://www.contributor-covenant.org/translations>.
# Contributing to FielDHub

<!-- This CONTRIBUTING.md is adapted from https://gist.github.com/peterdesmet/e90a1b0dc17af6c12daf6e8b2f044e7c -->

First of all, thanks for considering contributing to FielDHub! 👍 It's people like you that make it rewarding for us - the project maintainers - to work on FielDHub. 😊

FielDHub is an open-source project, maintained by people who care.

[repo]: https://github.com/DidierMurilloF/FielDHub
[issues]: https://github.com/DidierMurilloF/FielDHub/issues
[new_issue]: https://github.com/DidierMurilloF/FielDHub/issues/new
[website]: https://DidierMurilloF.github.io/FielDHub
[citation]: https://DidierMurilloF.github.io/FielDHub/authors.html
[email]: didier.murilloflorez@ndsu.edu

## Code of conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## How you can contribute

There are several ways you can contribute to this project. If you want to know more about why and how to contribute to open source projects like this one, see this [Open Source Guide](https://opensource.guide/how-to-contribute/).

### Share the love ❤️

Think FielDHub is useful? Let others discover it, by telling them in person, via Twitter or a blog post.

Using FielDHub for a paper you are writing? Consider [citing it][citation].

### Ask a question ⁉️

Using FielDHub and got stuck? Browse the [documentation][website] to see if you can find a solution. Still stuck? Post your question as an [issue on GitHub][new_issue]. While we cannot offer user support, we'll try to do our best to address it, as questions often lead to better documentation or the discovery of bugs.

Want to ask a question in private? Contact the package maintainer by [email][didier.murilloflorez@ndsu.edu].

### Propose an idea 💡

Have an idea for a new FielDHub feature? Take a look at the [documentation][website] and [issue list][issues] to see if it isn't included or suggested yet. If not, suggest your idea as an [issue on GitHub][new_issue]. While we can't promise to implement your idea, it helps to:

* Explain in detail how it would work.
* Keep the scope as narrow as possible.

See below if you want to contribute code for your idea as well.

### Report a bug 🐛

Using FielDHub and discovered a bug? That's annoying! Don't let others have the same experience and report it as an [issue on GitHub][new_issue] so we can fix it. A good bug report makes it easier for us to do so, so please include:

* Your operating system name and version (e.g. Mac OS 10.13.6).
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

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

Care to fix bugs or implement new functionality for FielDHub? Awesome! 👏 Have a look at the [issue list][issues] and leave a comment on the things you want to work on. See also the development guidelines below.

## Development guidelines

We try to follow the [GitHub flow](https://guides.github.com/introduction/flow/) for development.

1. Fork [this repo][repo] and clone it to your computer. To learn more about this process, see [this guide](https://guides.github.com/activities/forking/).
2. If you have forked and cloned the project before and it has been a while since you worked on it, [pull changes from the original repo](https://help.github.com/articles/merging-an-upstream-repository-into-your-fork/) to your clone by using `git pull upstream master`.
3. Open the RStudio project file (`.Rproj`).
4. Make your changes:
    * Write your code.
    * Test your code (bonus points for adding unit tests).
    * Document your code (see function documentation above).
    * Check your code with `devtools::check()` and aim for 0 errors, warnings and notes.
5. Commit and push your changes.
6. Submit a [pull request](https://guides.github.com/activities/forking/#making-a-pull-request).

## Attribution

This Contributing is adapted from [CONTRIBUTING](https://gist.github.com/peterdesmet/e90a1b0dc17af6c12daf6e8b2f044e7c).---
title: '`FielDHub`: A Shiny App for Design of Experiments in Life Sciences'
tags:
  - R
  - Shiny
  - DOE
  - Unreplicated Designs
  - Partially-Replicated Designs
  - Plant Breeding
authors:
  - name: Didier A. Murillo
    affiliation: 1
  - name: Salvador A. Gezan
    affiliation: 2
  - name: Ana M. Heilman
    affiliation: 1
  - name: Thomas C. Walk
    affiliation: 1
  - name: Johan S. Aparicio
    affiliation: 3
  - name: Richard D. Horsley
    affiliation: 1
affiliations:
 - name: North Dakota State University
   index: 1
 - name: VSN International
   index: 2
 - name: CIAT (International Center for Tropical Agriculture)
   index: 3
date: 4 March 2021
bibliography: paper.bib
---

# Summary

`FielDHub` is an R Shiny design of experiments (DOE) app that aids in the creation of traditional, unreplicated, augmented and partially-replicated [@Cullis2006_cit] designs  applied to agriculture, plant breeding, forestry, animal and biological sciences. One of the problems that life scientists often face is the lack of freely available and user-friendly interactive tools to create designs that fit their needs. A few open-source DOE R packages options exist including agricolae [@agricolae_cit] and blocksdesign [@blocksdesign_cit], but they require users to be familiar with the R programming language and do not have a graphical user interface (GUI).

# Statement of need

`FielDHub` allows users to perform randomizations of field, laboratory, and greenhouse experiments, while providing output via interactive field layouts and tables that can be extracted and saved. This app has a novel design that offers DOE options and features that are not currently available in most software applications.  Users are guided in each step of the DOE platform in an interactive interface, which includes a feature that helps to generate randomizations with an option to simulate data for a response variable. This last feature makes it suitable for teaching and evaluation purposes, where instructors can use the graphical dynamic user interface and/or use the functions included in the R package for teaching R scripting courses. This app also provides a graphical workflow to import treatment lists and export field books. For field experiments with a strict spatial arrangement, it allows users to specify the dimensions of the field (the number of rows and columns), while controlling the percentage of check plots, and obtaining field maps and field books that can be used directly as templates and input files for centralized databases.

`FielDHub` is currently being used by different breeding programs at NDSU and in graduate courses to teach the concepts of randomization, blocking, replication and simulations. The combination, in a single application, of novel and traditional designs, an interactive user interface, visualizations, and generation of templates and field layouts will enable the discovery of outstanding genotypes, while using efficient experimental designs that meet the requirements of the research being conducted.

Some of the features and designs implemented in `FielDHub` are summarized below:

1. Novel Designs: `FielDHub` has implemented a class of experimental designs known as augmented designs, partially replicated, and unreplicated designs. Examples are provided for each of the options with a default input data to demonstrate the functionalities of the app.

2. Reactive Interface: `FielDHub` provides output via an interactive interface, where users enter values that automatically generate tables, layouts, and output files within seconds.

3. Modularization: `FielDHub` was built in Shiny modules using the golem framework [@golem_cit]. Modularity makes the app easy to test, maintain, and deploy. 

4. Local and Remote Deployment: `FielDHub` can be deployed either to a local computer or to a server for online use. Currently it has been used within a server instance that has been utilized by graduate students and researchers alike in NDSU.

5. Simulations: `FielDHub` allows users to simulate a response variable along with the randomization. This feature can be used to define the corresponding linear model and to assess the efficiency of the experimental design, particularly in relation to its spatial components, or it can also be used to teach statistical concepts.

![\label{fig:Fig}](FielDHub_Overview_Map.png)
<div align="center"> Fig 1. Overview of `FielDHub` main features.</div>

# Usage

Plant breeding field research projects at NDSU are using `FielDHub` to refine experimental techniques in order to obtain unbiased and more precise estimates of the true treatment effects and their differences using unreplicated designs [@estefanova_cit; @Federer_cit]. Often these projects face limitations of seed quantity and available field space in conducting trials with large numbers of genotypes and opt for the use of partially replicated or unreplicated designs [@estefanova_cit]. As an example (Fig. 2), we consider here 270 genotypes arranged in a field of 15 rows by 20 columns. These genotypes are grouped in three different experiments/sites. In addition, we used four checks that are replicated in a systematic diagonal arrangement to fill 27 plots that represent 9% of the total experimental plots. An option to include filler plots is also available.

![\label{fig:Fig2}](Example_FielDHub4.png)
<div align="center"> Fig 2. Unreplicated design with checks in a systematic diagonal arrangement. </div>

# Acknowledgments

This application was developed as part of a collaboration between North Dakota State University (ND, USA), VSN International (Hemel Hempstead, UK) and CIAT (Cali, Colombia). The authors want to express the support from their respective institutions for allowing us to dedicate time to this exciting product. Special thanks to Dr. Blaine Johnson and Dr. Andrew Green for their useful contributions. 

# Availability and Community Guidelines

The software is available at the [GitHub](https://github.com/DidierMurilloF/FielDHub/) repository. The GitHub repository also contains the source code for this paper. Users and contributors are welcome to contribute, request features, and report bugs in this GitHub repository.

# References
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# FielDHub

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Travis build status](https://travis-ci.com/DidierMurilloF/FielDHub.svg?branch=master)](https://travis-ci.com/DidierMurilloF/FielDHub)
[![R-CMD-check](https://github.com/DidierMurilloF/FielDHub/workflows/R-CMD-check/badge.svg)](https://github.com/DidierMurilloF/FielDHub/actions)
<!-- badges: end -->

FielDHub is a customer centric R Shiny design of experiment (DOE) app that aids in the 
creation of non-traditional and traditional DOE. Our app provides a graphical workflow to import 
treatment lists and export field books. For field experiments, it allows users to specify the 
dimensions of the field as row and columns, while controlling the percentage of check plots, and 
obtaining randomized field maps and field books that can be used directly as templates and input 
files for central databases.

## Installation

You can install the dev version of FielDHub from:

``` r
devtools::install_bitbucket("DidierMurillo/fieldhub-package", 
                            auth_user = "DidierMurillo", 
                            password = rstudioapi::askForPassword())
```

## Example

This is a basic example which shows you how to launch the app:

```{r example}
#library(FielDHub)
#run_app()
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_app.R
\name{run_app}
\alias{run_app}
\title{Run the Shiny Application}
\usage{
run_app(...)
}
\arguments{
\item{...}{Unused, for extensibility}
}
\value{
A shiny app object
}
\description{
Run the Shiny Application
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_split_split_plot.R
\name{split_split_plot}
\alias{split_split_plot}
\title{Generates a Split Split Plot Design}
\usage{
split_split_plot(
  wp = NULL,
  sp = NULL,
  ssp = NULL,
  reps = NULL,
  type = 2,
  l = 1,
  plotNumber = 101,
  seed = NULL,
  locationNames = NULL,
  factorLabels = TRUE,
  data = NULL
)
}
\arguments{
\item{wp}{Number of whole plots, as an integer or a vector.}

\item{sp}{Number of sub plots per whole plot, as an integer or a vector.}

\item{ssp}{Number of sub-sub plots, as an integer or a vector.}

\item{reps}{Number of blocks (full replicates).}

\item{type}{Option for CRD or RCBD designs. Values are \code{type = 1} (CRD) or \code{type = 2} (RCBD). By default \code{type = 2}.}

\item{l}{Number of locations. By default \code{l = 1}.}

\item{plotNumber}{Numeric vector with the starting plot number for each location. By default \code{plotNumber = 101}.}

\item{seed}{(optional) Real number that specifies the starting seed to obtain reproducible designs.}

\item{locationNames}{(optional) Names for each location.}

\item{factorLabels}{(optional) If \code{TRUE} retain the levels
labels from the original data set otherwise, numeric labels will be
assigned. Default is \code{factorLabels =TRUE}.}

\item{data}{(optional) Data frame with label list of treatments.}
}
\value{
A list with two elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{fieldBook} is a data frame with the split split plot field book.
}
}
\description{
It randomly generates a split split plot design (SSPD) across locations.
}
\examples{
# Example 1: Generates a split split plot design SSPD with 5 whole plots, 2 sub-plots,
# 3 sub-sub plots, and 3 reps in an RCBD arrangement. This is for one location.
SSPD1 <- split_split_plot(wp = 4, sp = 2, ssp = 3, reps = 5, l = 1, 
                          plotNumber = 101, 
                          seed = 23, 
                          type = 2, 
                          locationNames = "FARGO")
SSPD1$infoDesign
head(SSPD1$fieldBook,12)

# Example 2: Generates a split split plot design SSPD with 2 whole plost 
# (Irrigation, No irrigation), 5 sub plots (4 types of fungicide + one control), and 
# 10 sub-sub plots (Ten varieties of beans), and 4 reps in an RCBD arrangement.
# This is for 3 locations. In this case, we show how to use the option data.
wp <- paste("IRR_", c("NO", "Yes"), sep = "") #Irrigation (2 Whole plots)
sp <- c("NFung", paste("Fung", 1:4, sep = "")) #Fungicides (5 Sub plots)
ssp <- paste("Beans", 1:10, sep = "") #Beans varieties (10 Sub-sub plots)
split_split_plot_Data <- data.frame(list(WHOLPLOT = c(wp, rep(NA, 8)), 
                                         SUBPLOT = c(sp, rep(NA, 5)),
                                         SUB_SUBPLOTS = ssp))
head(split_split_plot_Data, 10)
SSPD2 <- split_split_plot(reps = 4, l = 3, 
                          plotNumber = c(101, 1001, 2001),
                          seed = 23, 
                          type = 2, 
                          locationNames = c("A", "B", "C"),
                          data = split_split_plot_Data)
SSPD2$infoDesign
head(SSPD2$fieldBook,12)
             
}
\references{
Federer, W. T. (1955). Experimental Design. Theory and Application. New York, USA. The
Macmillan Company.
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_diagonal_arrangement.R
\name{diagonal_arrangement}
\alias{diagonal_arrangement}
\title{Spatial Un-replicated Diagonal Arrangement Design}
\usage{
diagonal_arrangement(
  nrows = NULL,
  ncols = NULL,
  lines = NULL,
  checks = NULL,
  planter = "serpentine",
  l = 1,
  plotNumber = 101,
  kindExpt = "SUDC",
  splitBy = "row",
  seed = NULL,
  blocks = NULL,
  exptName = NULL,
  locationNames = NULL,
  data = NULL
)
}
\arguments{
\item{nrows}{Number of rows in the field.}

\item{ncols}{Number of columns in the field.}

\item{lines}{Number of genotypes, experimental lines or treatments.}

\item{checks}{Number of genotypes checks.}

\item{planter}{Option for \code{serpentine} or \code{cartesian} plot arrangement. 
By default  \code{planter = 'serpentine'}.}

\item{l}{Number of locations or sites. By default  \code{l = 1}.}

\item{plotNumber}{Numeric vector with the starting plot number for each location. 
By default \code{plotNumber = 101}.}

\item{kindExpt}{Type of diagonal design, with single options: Single Un-replicated Diagonal Checks
\code{'SUDC'} and Decision Blocks Un-replicated Design with Diagonal Checks \code{'DBUDC'} 
for multiple experiments. By default \code{kindExpt = 'SUDC'}.}

\item{splitBy}{Option to split the field when \code{kindExpt = 'DBUDC'} is selected. 
By default \code{splitBy = 'row'}.}

\item{seed}{(optional) Real number that specifies the starting seed to obtain reproducible designs.}

\item{blocks}{Number of experiments or blocks to generate an \code{DBUDC} design. 
If \code{kindExpt = 'DBUDC'} and data is null, \code{blocks} are mandatory.}

\item{exptName}{(optional) Name of the experiment.}

\item{locationNames}{(optional) Names each location.}

\item{data}{(optional) Data frame with 3 columns: \code{ENTRY | NAME | BLOCK} or only 2 
columns \code{ENTRY | NAME} if \code{kindExpt = 'SUDC'}.}
}
\value{
A list with five elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{layoutRandom} is a matrix with the randomization layout.
  \item \code{plotsNumber} is a matrix with the layout plot number.
  \item \code{data_entry} is a data frame with the data input.
  \item \code{fieldBook} is a data frame with field book design. This includes the index (Row, Column).
}
}
\description{
Randomly generates an spatial un-replicated diagonal arrangement design.
}
\examples{

# Example 1: Generates a spatial single diagonal arrangement design in one location
# with 270 treatments and 30 check plots for a field with dimensions 15 rows x 20 cols
# in a serpentine arrangement.
spatd <- diagonal_arrangement(nrows = 15, ncols = 20, lines = 270, 
                              checks = 4, 
                              plotNumber = 101, 
                              kindExpt = "SUDC", 
                              planter = "serpentine", 
                              seed = 1987,
                              exptName = "20WRY1", 
                              locationNames = "MINOT")
spatd$infoDesign
spatd$layoutRandom
spatd$plotsNumber
head(spatd$fieldBook, 12)

# Example 2: Generates a spatial decision block diagonal arrangement design in one location
# with 720 treatments allocated in 5 experiments or blocks for a field with dimensions
# 30 rows x 26 cols in a serpentine arrangement. In this case, we show how to set up the data 
# option with the entries list.
checks <- 5;expts <- 5
list_checks <- paste("CH", 1:checks, sep = "")
treatments <- paste("G", 6:725, sep = "")
BLOCK <- c(rep("ALL", checks), rep(1:expts, c(150,155,95,200,120)))
treatment_list <- data.frame(list(ENTRY = 1:725, NAME = c(list_checks, treatments), BLOCK = BLOCK))
head(treatment_list, 12) 
tail(treatment_list, 12)
spatDB <- diagonal_arrangement(nrows = 30, ncols = 26,
                               checks = 5, 
                               plotNumber = 1, 
                               kindExpt = "DBUDC", 
                               planter = "serpentine", 
                               splitBy = "row", 
                               data = treatment_list)
spatDB$infoDesign
spatDB$layoutRandom
spatDB$plotsNumber
head(spatDB$fieldBook,12)

# Example 3: Generates a spatial decision block diagonal arrangement design in one location
# with 270 treatments allocated in 3 experiments or blocks for a field with dimensions
# 20 rows x 15 cols in a serpentine arrangement. Which in turn is an augmented block (3 blocks).
spatAB <- diagonal_arrangement(nrows = 20, ncols = 15, lines = 270, 
                               checks = 4, 
                               plotNumber = c(1,1001,2001), 
                               kindExpt = "DBUDC", 
                               planter = "serpentine",
                               exptName = c("20WRA", "20WRB", "20WRC"), 
                               blocks = c(90, 90, 90),
                               splitBy = "column")
spatAB$infoDesign
spatAB$layoutRandom
spatAB$plotsNumber
head(spatAB$fieldBook,12)

}
\references{
Clarke, G. P. Y., & Stefanova, K. T. (2011). Optimal design for early-generation plant
breeding trials with unreplicated or partially replicated test lines. Australian & New
Zealand Journal of Statistics, 53(4), 461–480.
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_partially_replicated.R
\name{partially_replicated}
\alias{partially_replicated}
\title{Generates a Spatial Partially Replicated Arrangement Design}
\usage{
partially_replicated(
  nrows = NULL,
  ncols = NULL,
  repGens = NULL,
  repUnits = NULL,
  planter = "serpentine",
  l = 1,
  plotNumber = 101,
  seed = NULL,
  exptName = NULL,
  locationNames = NULL,
  data = NULL
)
}
\arguments{
\item{nrows}{Number of rows field.}

\item{ncols}{Number of columns field.}

\item{repGens}{Numeric vector with the amount genotypes to replicate.}

\item{repUnits}{Numeric vector with the number of reps of each genotype.}

\item{planter}{Option for \code{serpentine} or \code{cartesian} movement. By default  \code{planter = 'serpentine'}.}

\item{l}{Number of locations. By default \code{l = 1}.}

\item{plotNumber}{Numeric vector with the starting plot number for each location. By default \code{plotNumber = 101}.}

\item{seed}{(optional) Real number that specifies the starting seed to obtain reproducible designs.}

\item{exptName}{(optional) Name of the experiment.}

\item{locationNames}{(optional) Name for each location.}

\item{data}{(optional) Dataframe with 3 columns: \code{ENTRY | NAME | REPS}.}
}
\value{
A list with five elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{layoutRandom} is a matrix with the randomization layout.
  \item \code{plotNumber} is a matrix with the layout plot number.
  \item \code{data_entry} is a data frame with the data input.
  \item \code{fieldBook} is a data frame with field book design. This includes the index (Row, Column).
}
}
\description{
Randomly generates a spatial partially replicated design, where the distance 
between checks is maximized in such a way that each row and column have control plots. 
Note that design generation needs the dimension of the field (number of rows and columns).
}
\examples{
# Example 1: Generates a spatial optimized partially replicated arrangement design in one 
# location with 342 genotypes for a field with dimensions 25 rows x 18 cols. 
# Note that there are 280 genotypes unreplicated (only one time), 50 genotypes replicated 
# two times, and 10 genotypes replicated three times, and two checks 20 times each one.
SpatpREP1 <- partially_replicated(nrows = 25, 
                                  ncols = 18,  
                                  repGens = c(280,50,10,1,1),
                                  repUnits = c(1,2,3,20,20),
                                  planter = "cartesian", 
                                  plotNumber = 101,
                                  seed = 77)
SpatpREP1$infoDesign
SpatpREP1$layoutRandom
SpatpREP1$plotNumber
head(SpatpREP1$fieldBook,12)

# Example 2: Generates a spatial optimized partially replicated arrangement design with 492 
# genotypes in a field with dimensions 30 rows x 20 cols. Note that there 384 genotypes 
# unreplicated (only one time), 108 genotypes replicated two times. 
# In this case we don't have check plots.
# As example, we set up the data option with the entries list.
NAME <- paste("G", 1:492, sep = "")
repGens = c(108, 384);repUnits = c(2,1)
REPS <- rep(repUnits, repGens)
treatment_list <- data.frame(list(ENTRY = 1:492, NAME = NAME, REPS = REPS))
head(treatment_list, 12) 
tail(treatment_list, 12)
SpatpREP2 <- partially_replicated(nrows = 30, 
                                  ncols = 20, 
                                  planter = "serpentine", 
                                  plotNumber = 101,
                                  seed = 41,
                                  data = treatment_list)
SpatpREP2$infoDesign
SpatpREP2$layoutRandom
SpatpREP2$plotNumber
head(SpatpREP2$fieldBook,10)

}
\references{
Cullis, S., B. R., & Coombes, N. E. (2006). On the design of early generation variety trials
with correlated data. Journal of Agricultural, Biological, and Environmental Statistics, 11,
381–393. https://doi.org/10.1198/108571106X154443
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_S3_methods.R
\name{print.FielDHub}
\alias{print.FielDHub}
\title{Print a \code{FielDHub} object}
\usage{
\method{print}{FielDHub}(x, n, ...)
}
\arguments{
\item{x}{an object inheriting from class}

\item{n}{a single integer. If positive or zero, size for the
resulting object: number of elements for a vector (including
lists), rows for a matrix or data frame or lines for a function. If
negative, all but the n last/first number of elements of x.}

\item{...}{further arguments passed to \code{\link{head}}.}
}
\value{
an object inheriting from class \code{FielDHub}
}
\description{
Prints information about any \code{FielDHub} function.
}
\examples{
# Example 1: Generates a CRD design with 5 treatments and 5 reps each.
crd1 <- CRD(t = 5, reps = 5, plotNumber = 101,
seed = 1985, locationName = "Fargo")
crd1$infoDesign
print(crd1)

}
\author{
Thiago de Paula Oliveira,
  \email{thiago.paula.oliveira@alumni.usp.br} [aut],
  Didier Murillo [aut]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_RCDB_augmented.R
\name{RCBD_augmented}
\alias{RCBD_augmented}
\title{Generates an Augmented Randomized Complete Block Design (ARCBD)}
\usage{
RCBD_augmented(
  lines = NULL,
  checks = NULL,
  b = NULL,
  l = 1,
  planter = "serpentine",
  plotNumber = 101,
  exptName = NULL,
  seed = NULL,
  locationNames = NULL,
  repsExpt = 1,
  random = TRUE,
  data = NULL
)
}
\arguments{
\item{lines}{Treatments, number of lines for test.}

\item{checks}{Number of checks per augmented block.}

\item{b}{Number of augmented blocks.}

\item{l}{Number of locations. By default \code{l = 1}.}

\item{planter}{Option for \code{serpentine} or \code{cartesian} arrangement. By default \code{planter = 'serpentine'}.}

\item{plotNumber}{Numeric vector with the starting plot number for each location. By default \code{plotNumber = 101}.}

\item{exptName}{(optional) Name of experiment.}

\item{seed}{(optional) Real number that specifies the starting seed to obtain reproducible designs.}

\item{locationNames}{(optional) Name for each location.}

\item{repsExpt}{(optional) Number of reps of experiment. By default \code{repsExpt = 1}.}

\item{random}{Logical value to randomize treatments or not. By default \code{random = TRUE}.}

\item{data}{(optional) Data frame with the labels of treatments.}
}
\value{
A list with five elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{layoutRandom} is the ARCBD layout randomization for the first location.
  \item \code{plotNumber} is the plot number layout for the first location.
  \item \code{exptNames} is the experiment names layout.
  \item \code{data_entry} is a data frame with the data input.
  \item \code{fieldBook} is a data frame with the ARCBD field book.
}
}
\description{
It randomly generates an augmented randomized complete block design across locations (ARCBD).
}
\examples{
# Example 1: Generates an ARCBD with 6 blocks, 3 checks for each, and 50 treatments 
# in two locations.
ARCBD1 <- RCBD_augmented(lines = 50, checks = 3, b = 6, l = 2, 
                         planter = "cartesian", 
                         plotNumber = c(1,1001),
                         seed = 23, 
                         locationNames = c("FARGO", "MINOT"))
ARCBD1$infoDesign
ARCBD1$layoutRandom
ARCBD1$exptNames
ARCBD1$plotNumber
head(ARCBD1$fieldBook, 12)
                   
# Example 2: Generates an ARCBD with 17 blocks, 4 checks for each, and 350 treatments 
# in 3 locations.
# In this case, we show how to use the option data.
checks <- 4;
list_checks <- paste("CH", 1:checks, sep = "")
treatments <- paste("G", 5:354, sep = "")
treatment_list <- data.frame(list(ENTRY = 1:354, NAME = c(list_checks, treatments)))
head(treatment_list, 12)
ARCBD2 <- RCBD_augmented(lines = 350, checks = 4, b = 17, l = 3, 
                         planter = "serpentine", 
                         plotNumber = c(101,1001,2001), 
                         seed = 24, 
                         locationNames = LETTERS[1:3],
                         data = treatment_list)
ARCBD2$infoDesign
ARCBD2$layoutRandom
ARCBD2$exptNames
ARCBD2$plotNumber
head(ARCBD2$fieldBook, 12)
                                       
}
\references{
Federer, W. T. (1955). Experimental Design. Theory and Application. New York, USA. The
Macmillan Company.
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_latin_square.R
\name{latin_square}
\alias{latin_square}
\title{Generates a Latin Square Design}
\usage{
latin_square(
  t = NULL,
  reps = 1,
  plotNumber = 101,
  planter = "serpentine",
  seed = NULL,
  locationNames = NULL,
  data = NULL
)
}
\arguments{
\item{t}{Number of treatments.}

\item{reps}{Number of full resolvable squares. By default \code{reps = 1}.}

\item{plotNumber}{Starting plot number. By default \code{plotNumber = 101}.}

\item{planter}{Option for \code{serpentine} or \code{cartesian} arrangement. By default \code{planter = 'serpentine'}.}

\item{seed}{(optional) Real number that specifies the starting seed to obtain reproducible designs.}

\item{locationNames}{(optional) Name for the location.}

\item{data}{(optional) Data frame with label list of treatments.}
}
\value{
A list with information on the design parameters.

Data frame with the latin square field book.

A list with two elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{fieldBook} is a data frame with the latin square field book.
}
}
\description{
Randomly generates a latin square design of up 10 treatments.
}
\examples{
# Example 1: Generates a latin square design with 4 treatments and 2 reps.
latinSq1 <- latin_square(t = 4,
                         reps = 2,
                         plotNumber = 101,
                         planter = "cartesian",
                         seed = 1980)
print(latinSq1)
summary(latinSq1)
head(latinSq1$fieldBook)

# Example 2: Generates a latin square design with 5 treatments and 3 reps.
latin_data <- data.frame(list(ROW = paste("Period", 1:5, sep = ""),
                              COLUMN = paste("Cow", 1:5, sep = ""),
                              TREATMENT = paste("Diet", 1:5, sep = "")))
print(latin_data)
latinSq2 <- latin_square(t = NULL,
                         reps = 3,
                         plotNumber = 101,
                         planter = "cartesian",
                         seed = 1981,
                         data = latin_data)
latinSq2$squares
latinSq2$plotSquares
head(latinSq2$fieldBook)

}
\references{
Federer, W. T. (1955). Experimental Design. Theory and Application. New York, USA. The
Macmillan Company.
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb],
        Thiago de Paula Oliveira[ctb] 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_square_lattice.R
\name{square_lattice}
\alias{square_lattice}
\title{Generates a Square Lattice Design.}
\usage{
square_lattice(
  t = NULL,
  k = NULL,
  r = NULL,
  l = 1,
  plotNumber = 101,
  locationNames = NULL,
  seed = NULL,
  data = NULL
)
}
\arguments{
\item{t}{Number of  treatments.}

\item{k}{Size of incomplete blocks (number of units per incomplete block).}

\item{r}{Number of blocks (full resolvable replicates).}

\item{l}{Number of locations. By default \code{l = 1}.}

\item{plotNumber}{Numeric vector with the starting plot number for each location. By default \code{plotNumber = 101}.}

\item{locationNames}{(optional) Names for each location.}

\item{seed}{(optional) Real number that specifies the starting seed to obtain reproducible designs.}

\item{data}{(optional) Data frame with label list of treatments.}
}
\value{
A list with two elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{fieldBook} is a data frame with the square lattice design field book.
}
}
\description{
It randomly generates a square lattice design across locations.
}
\examples{
# Example 1: Generates a square lattice design with 5 full blocks, 8 units per IBlock,
# 8 IBlocks for a square number of treatmens of 64 in two locations.
squareLattice1 <- square_lattice(t = 64, k = 8, r = 5, l = 2, 
                                 plotNumber = c(1001, 2001),
                                 locationNames = c("FARGO", "MINOT"), 
                                 seed = 1986)
squareLattice1$infoDesign
head(squareLattice1$fieldBook,12)

# Example 2: Generates a square lattice design with 3 full blocks, 7 units per IBlock,
# 7 IBlocks for a square number of treatmens of 49 in one location.
# In this case, we show how to use the option data.
treatments <- paste("G", 1:49, sep = "")
ENTRY <- 1:49
treatment_list <- data.frame(list(ENTRY = ENTRY, TREATMENT = treatments))
head(treatment_list) 
squareLattice2 <- square_lattice(t = 49, k = 7, r = 3, l = 1, 
                                 plotNumber = 1001,
                                 locationNames = "CASSELTON", 
                                 seed = 1986,
                                 data = treatment_list)
squareLattice2$infoDesign
head(squareLattice2$fieldBook,12)

}
\references{
Edmondson., R. N. (2021). blocksdesign: Nested and crossed block designs for factorial and
unstructured treatment sets. https://CRAN.R-project.org/package=blocksdesign
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_alpha_lattice.R
\name{alpha_lattice}
\alias{alpha_lattice}
\title{Generates an Alpha Design}
\usage{
alpha_lattice(
  t = NULL,
  k = NULL,
  r = NULL,
  l = 1,
  plotNumber = 101,
  locationNames = NULL,
  seed = NULL,
  data = NULL
)
}
\arguments{
\item{t}{Number of  treatments.}

\item{k}{Size of incomplete blocks (number of units per incomplete block).}

\item{r}{Number of full blocks (or resolvable replicates) (also number of replicates per treatment).}

\item{l}{Number of locations. By default \code{l = 1}.}

\item{plotNumber}{Numeric vector with the starting plot number for each location. By default \code{plotNumber = 101}.}

\item{locationNames}{(optional) String with names for each of the \code{l} locations.}

\item{seed}{(optional) Real number that specifies the starting seed to obtain reproducible designs.}

\item{data}{(optional) Data frame with label list of treatments.}
}
\value{
A list with two elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{fieldBook} is a data frame with the alpha design field book.
}
}
\description{
Randomly generates an alpha design like \code{alpha(0,1)} across multiple locations.
}
\examples{
# Example 1: Generates an alpha design with 7 full blocks and 15 treatments.
# Size of IBlocks k = 3.
alphalattice1 <- alpha_lattice(t = 15, k = 3, r = 7, 
                               l = 1, 
                               plotNumber = 101, 
                               locationNames = "GreenHouse", 
                               seed = 1247)
alphalattice1$infoDesign
head(alphalattice1$fieldBook, 10)

# Example 2: Generates an alpha design with 5 full blocks and 50 treatment.
# Size of IBlocks k = 10. 
# In this case, we show how to use the option data.
treatments <- paste("G-", 1:50, sep = "")
ENTRY <- 1:50
treatment_list <- data.frame(list(ENTRY = ENTRY, TREATMENT = treatments))
head(treatment_list) 
alphalattice2 <- alpha_lattice(t = 50, k = 10, r = 5, 
                               l = 1, 
                               plotNumber = 1001, 
                               locationNames = "A", 
                               seed = 1945,
                               data = treatment_list)
alphalattice2$infoDesign
head(alphalattice2$fieldBook, 10)

}
\references{
Edmondson., R. N. (2021). blocksdesign: Nested and crossed block designs for factorial and
unstructured treatment sets. https://CRAN.R-project.org/package=blocksdesign
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_optimized_arrangement.R
\name{optimized_arrangement}
\alias{optimized_arrangement}
\title{Generates an Spatial Un-replicated Optimized Arrangement Design}
\usage{
optimized_arrangement(
  nrows = NULL,
  ncols = NULL,
  lines = NULL,
  amountChecks = NULL,
  checks = NULL,
  planter = "serpentine",
  l = 1,
  plotNumber = 101,
  seed = NULL,
  exptName = NULL,
  locationNames = NULL,
  data = NULL
)
}
\arguments{
\item{nrows}{Number of rows in the field.}

\item{ncols}{Number of columns in the field.}

\item{lines}{Number of genotypes, experimental lines or treatments.}

\item{amountChecks}{Integer with the amount total of checks or a numeric vector with the replicates of each check label.}

\item{checks}{Number of genotypes as checks.}

\item{planter}{Option for \code{serpentine} or \code{cartesian} arrangement. By default  \code{planter = 'serpentine'}.}

\item{l}{Number of locations. By default \code{l = 1}.}

\item{plotNumber}{Numeric vector with the starting plot number for each location. By default \code{plotNumber = 101}.}

\item{seed}{(optional) Real number that specifies the starting seed to obtain reproducible designs.}

\item{exptName}{(optional) Name of the experiment.}

\item{locationNames}{(optional) Name for each location.}

\item{data}{(optional) Data frame with 3 columns: \code{ENTRY | NAME | REPS}.}
}
\value{
A list with five elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{layoutRandom} is a matrix with the randomization layout.
  \item \code{plotNumber} is a matrix with the layout plot number.
  \item \code{data_entry} is a data frame with the data input.
  \item \code{fieldBook} is a data frame with field book design. This includes the index (Row, Column).
}
}
\description{
Randomly generates a spatial un-replicated optimized arrangement design, where the distance
between checks is maximized in such a way that each row and column have control plots.
Note that design generation needs the dimension of the field (number of rows and columns).
}
\examples{
# Example 1: Generates a spatial unreplicated optimized arrangement design in one location
# with 362 genotypes + 38 check plots (5 checks) for a field with dimension 20 rows x 20 cols.
OptimAd1 <- optimized_arrangement(nrows = 20, ncols = 20, lines = 362, 
                                  amountChecks = 38, 
                                  checks = 1:5,
                                  planter = "cartesian", 
                                  plotNumber = 101,
                                  seed = 14,
                                  exptName = "20RW1",
                                  locationNames = "CASSELTON")
OptimAd1$infoDesign
OptimAd1$layoutRandom
OptimAd1$plotNumber
head(OptimAd1$fieldBook,12)
                  
# Example 2: Generates a spatial unreplicated optimized arrangement design in one location
# with 635 genotypes + 65 check plots (4 checks) for a field with dimension 20 rows x 35 cols.
# As example, we set up the data option with the entries list.
checks <- 4
list_checks <- paste("CH", 1:checks, sep = "")
treatments <- paste("G", 5:639, sep = "")
REPS <- c(17, 16, 16, 16, rep(1, 635))
treatment_list <- data.frame(list(ENTRY = 1:639, NAME = c(list_checks, treatments), REPS = REPS))
head(treatment_list, 12) 
tail(treatment_list, 12)
OptimAd2 <- optimized_arrangement(nrows = 20, ncols = 35, 
                                  planter = "serpentine", 
                                  plotNumber = 101,
                                  seed = 12,
                                  exptName = "20YWA2",
                                  locationNames = "MINOT",
                                  data = treatment_list)
OptimAd2$infoDesign
OptimAd2$layoutRandom
OptimAd2$plotNumber
head(OptimAd2$fieldBook,12)
                  
}
\references{
Clarke, G. P. Y., & Stefanova, K. T. (2011). Optimal design for early-generation plant
breeding trials with unreplicated or partially replicated test lines. Australian & New
Zealand Journal of Statistics, 53(4), 461–480.
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_rectangular_lattice.R
\name{rectangular_lattice}
\alias{rectangular_lattice}
\title{Generates a Rectangular Lattice Design.}
\usage{
rectangular_lattice(
  t = NULL,
  k = NULL,
  r = NULL,
  l = 1,
  plotNumber = 101,
  locationNames = NULL,
  seed = NULL,
  data = NULL
)
}
\arguments{
\item{t}{Number of  treatments.}

\item{k}{Size of incomplete blocks (number of units per incomplete block).}

\item{r}{Number of blocks (full resolvable replicates).}

\item{l}{Number of locations. By default \code{l = 1}.}

\item{plotNumber}{Numeric vector with the starting plot number for each location. By default \code{plotNumber = 101}.}

\item{locationNames}{(optional) Names for each location.}

\item{seed}{(optional) Real number that specifies the starting seed to obtain reproducible designs.}

\item{data}{(optional) Data frame with label list of treatments.}
}
\value{
A list with two elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{fieldBook} is a data frame with the rectangular lattice design field book.
}
}
\description{
It randomly generates a rectangular lattice design across locations.
}
\examples{
# Example 1: Generates a rectangular lattice design with 6 full blocks, 4 units per IBlock (k)
# and 20 treatments in one location.
rectangularLattice1 <- rectangular_lattice(t = 20, k = 4, r = 6, l = 1, 
                                           plotNumber = 101,
                                           locationNames = "FARGO", 
                                           seed = 126)
rectangularLattice1$infoDesign
head(rectangularLattice1$fieldBook,12)

# Example 2: Generates a rectangular lattice design with 5 full blocks, 7 units per IBlock (k)
# and 56 treatments across 2 locations.
# In this case, we show how to use the option data.
treatments <- paste("ND-", 1:56, sep = "")
ENTRY <- 1:56
treatment_list <- data.frame(list(ENTRY = ENTRY, TREATMENT = treatments))
head(treatment_list) 
rectangularLattice2 <- rectangular_lattice(t = 56, k = 7, r = 5, l = 2, 
                                           plotNumber = c(1001,2001),
                                           locationNames = c("Loc1", "Loc2"), 
                                           seed = 127,
                                           data = treatment_list)
rectangularLattice2$infoDesign
head(rectangularLattice2$fieldBook,12)

}
\references{
Edmondson., R. N. (2021). blocksdesign: Nested and crossed block designs for factorial and
unstructured treatment sets. https://CRAN.R-project.org/package=blocksdesign
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_row_column.R
\name{row_column}
\alias{row_column}
\title{Generates a Resolvable Row-Column Design (RowColD)}
\usage{
row_column(
  t = NULL,
  nrows = NULL,
  r = NULL,
  l = 1,
  plotNumber = 101,
  locationNames = NULL,
  seed = NULL,
  data = NULL
)
}
\arguments{
\item{t}{Number of  treatments.}

\item{nrows}{Number of rows of a full resolvable replicate.}

\item{r}{Number of blocks (full resolvable replicates).}

\item{l}{Number of locations. By default \code{l = 1}.}

\item{plotNumber}{Numeric vector with the starting plot number for each location. By default \code{plotNumber = 101}.}

\item{locationNames}{(optional) Names for each location.}

\item{seed}{(optional) Real number that specifies the starting seed to obtain reproducible designs.}

\item{data}{(optional) Data frame with label list of treatments}
}
\value{
A list with four elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{resolvableBlocks} a list with the resolvable row columns blocks. 
  \item \code{concurrence} is the concurrence matrix.
  \item \code{fieldBook} is a data frame with the row-column field book.
}
}
\description{
It randomly generates a resolvable row-column designs (RowColD). 
Note that design optimization is only done at the level of rows and not columns; 
hence, design is suboptimal. The randomization can be done across locations.
}
\examples{

# Example 1: Generates a row-column design with 3 full blocks and 36 treatments
# and 6 rows. This for one location.
rowcold1 <- row_column(t = 36, nrows = 6, r = 3, l = 1, 
                       plotNumber= 101, 
                       locationNames = "Loc1",
                       seed = 21)
rowcold1$infoDesign
rowcold1$resolvableBlocks
head(rowcold1$fieldBook,12)

# Example 2: Generates a row-column design with 3 full blocks and 30 treatments
# and 5 rows, for two locations.
# In this case, we show how to use the option data.
treatments <- paste("ND-", 1:30, sep = "")
ENTRY <- 1:30
treatment_list <- data.frame(list(ENTRY = ENTRY, TREATMENT = treatments))
head(treatment_list)
rowcold2 <- row_column(t = 30, nrows = 5, r = 3, l = 2, 
                       plotNumber= c(101,1001), 
                       locationNames = c("A", "B"),
                       seed = 15,
                       data = treatment_list)
rowcold2$infoDesign
rowcold2$resolvableBlocks
head(rowcold2$fieldBook,12)
  

}
\references{
Edmondson., R. N. (2021). blocksdesign: Nested and crossed block designs for factorial and
unstructured treatment sets. https://CRAN.R-project.org/package=blocksdesign
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_split_families.R
\name{split_families}
\alias{split_families}
\title{Split a population of genotypes randomly into several locations.}
\usage{
split_families(l = NULL, data = NULL)
}
\arguments{
\item{l}{Number of locations.}

\item{data}{Data frame with the entry (ENTRY) and the labels of each treatment (NAME)
and number of individuals per family group (FAMILY).}
}
\value{
A list with two elements.
\itemize{
  \item \code{rowsEachlist} is a table with a summary of cases.
  \item \code{data_locations} is a data frame with the entries for each location
}
}
\description{
Split a population of genotypes randomly into several locations, with the
aim of having approximatelly the same number of replicates of each genotype, line or
treatment per location.
}
\examples{
# Example 1: Split a population of 3000 and 200 families into 8 locations. 
# Original dataset is been simulated.
set.seed(77)
N <- 3000; families <- 200
ENTRY <- 1:N
NAME <- paste0("SB-", 1:N)
FAMILY <- vector(mode = "numeric", length = N)
x <- 1:N
for (i in x) { FAMILY[i] <- sample(1:families, size = 1, replace = TRUE) }
gen.list <- data.frame(list(ENTRY = ENTRY, NAME = NAME, FAMILY = FAMILY))
head(gen.list)
# Now we are going to use the split_families() function.
split_population <- split_families(l = 8, data = gen.list)
print(split_population)
summary(split_population)
head(split_population$data_locations,12)

}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_strip_plot.R
\name{strip_plot}
\alias{strip_plot}
\title{Strip Plot Design}
\usage{
strip_plot(
  Hplots = NULL,
  Vplots = NULL,
  b = 1,
  l = 1,
  plotNumber = NULL,
  planter = "serpentine",
  locationNames = NULL,
  seed = NULL,
  factorLabels = TRUE,
  data = NULL
)
}
\arguments{
\item{Hplots}{Number of horizontal factors, as an integer or a vector.}

\item{Vplots}{Number of vertical factors, as an integer or a vector.}

\item{b}{Number of blocks (full replicates).}

\item{l}{Number of locations. By default \code{l = 1}.}

\item{plotNumber}{Numeric vector with the starting plot number for each location. By default \code{plotNumber = 101}.}

\item{planter}{Option for \code{serpentine} or \code{cartesian} arrangement. By default \code{planter = 'serpentine'}.}

\item{locationNames}{(optional) Names for each location.}

\item{seed}{(optional) Real number that specifies the starting seed to obtain reproducible designs.}

\item{factorLabels}{(optional) If \code{TRUE} retain the levels
labels from the original data set otherwise, numeric labels will be
assigned. Default is \code{factorLabels =TRUE}.}

\item{data}{(optional) data frame with the labels of vertical and hirizontal plots.}
}
\value{
A list with four elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{stripsBlockLoc} is a list with the strip blocks for each location.
  \item \code{plotLayouts} is a list with the layout plot numbers for each location.
  \item \code{fieldBook} is a data frame  with the strip plot field book.
}
}
\description{
It randomly generates a strip plot design across locations.
}
\examples{
# Example 1: Generates a strip plot design with 5 vertical strips and 4 horizontal strips,
# with 3 reps in one location.
H <- paste("H", 1:4, sep = "")
V <- paste("V", 1:5, sep = "")
strip1 <- strip_plot(Hplots = H, 
                     Vplots = V, 
                     b = 3, 
                     l = 1, 
                     plotNumber = 101,
                     planter = "serpentine",
                     locationNames = "A", 
                     seed = 333)
strip1$infoDesign                  
strip1$stripsBlockLoc
strip1$plotLayouts
head(strip1$fieldBook,12)                     

# Example 2: Generates a strip plot design with 5 vertical strips and 5 horizontal strips,
# with 6 reps across to 3 locations. In this case, we show how to use the option data.
Hplots <- LETTERS[1:5]
Vplots <- LETTERS[1:4]
strip_data <- data.frame(list(HPLOTS = Hplots, VPLOTS = c(Vplots, NA)))
head(strip_data)
strip2 <- strip_plot(Hplots = 5, 
                     Vplots = 5, 
                     b = 6, 
                     l = 3, 
                     plotNumber = c(101,1001,2001),
                     planter = "cartesian",
                     locationNames = c("A", "B", "C"), 
                     seed = 222,
                     data = strip_data)
strip2$infoDesign                  
strip2$stripsBlockLoc
strip2$plotLayouts
head(strip2$fieldBook,12)

}
\references{
Federer, W. T. (1955). Experimental Design. Theory and Application. New York, USA. The
Macmillan Company.
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\arguments{
\item{lhs}{A value or the magrittr placeholder.}

\item{rhs}{A function call using the magrittr semantics.}
}
\value{
No return value, called for side effects.
}
\description{
See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_full_factorial.R
\name{full_factorial}
\alias{full_factorial}
\title{Generates a Full Factorial Design}
\usage{
full_factorial(
  setfactors = NULL,
  reps = NULL,
  l = 1,
  type = 2,
  plotNumber = 101,
  continuous = FALSE,
  planter = "serpentine",
  seed = NULL,
  locationNames = NULL,
  factorLabels = TRUE,
  data = NULL
)
}
\arguments{
\item{setfactors}{Numeric vector with levels of each factor.}

\item{reps}{Number of replicates (full blocks).}

\item{l}{Number of locations. By default \code{l = 1}.}

\item{type}{Option for CRD or RCBD designs. Values are \code{type =
1} (CRD) or \code{type = 2} (RCBD). By default \code{type = 2}.}

\item{plotNumber}{Numeric vector with the starting plot number for
each location. By default \code{plotNumber = 101}.}

\item{continuous}{Logical for plot number continuous or not. By
default \code{continuous = FALSE}.}

\item{planter}{Option for \code{serpentine} or \code{cartesian} plot
arrangement. By default \code{planter = 'serpentine'}.}

\item{seed}{(optional) Real number that specifies the starting seed
to obtain reproducible designs.}

\item{locationNames}{(optional) Names for each location.}

\item{factorLabels}{(optional) If \code{TRUE} retain the levels
labels from the original data set otherwise, numeric labels will be
assigned. Default is \code{factorLabels =TRUE}.}

\item{data}{(optional) Data frame with the labels of factors.}
}
\value{
A list with two elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{fieldBook} is a data frame with the full factorial field book.
}
}
\description{
It randomly generates a full factorial design across locations.
}
\examples{
# Example 1: Generates a full factorial with 3 factors each with 2 levels.
# This in an RCBD arrangement with 3 reps.
fullFact1 <- full_factorial(setfactors = c(2,2,2), reps = 3, l = 1, type = 2,
                            plotNumber = 101,
                            continuous = TRUE,
                            planter = "serpentine",
                            seed = 325,
                            locationNames = "FARGO")
fullFact1$infoDesign
head(fullFact1$fieldBook,10)

# Example 2: Generates a full factorial with 3 factors and each with levels: 2,3,
# and 2, respectively. In this case, we show how to use the option data
FACTORS <- rep(c("A", "B", "C"), c(2,3,2))
LEVELS <- c("a0", "a1", "b0", "b1", "b2", "c0", "c1")
data_factorial <- data.frame(list(FACTOR = FACTORS, LEVEL = LEVELS))
print(data_factorial)
# This in an RCBD arrangement with 5 reps in 3 locations.
fullFact2 <- full_factorial(setfactors = NULL, reps = 5, l = 3, type = 2,
                            plotNumber = c(101,1001,2001),
                            continuous = FALSE,
                            planter = "serpentine",
                            seed = 326,
                            locationNames = c("Loc1","Loc2","Loc3"),
                            data = data_factorial)
fullFact2$infoDesign
head(fullFact2$fieldBook,10)

}
\references{
Federer, W. T. (1955). Experimental Design. Theory and Application. New York, USA. The
Macmillan Company.
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_incomplete_blocks.R
\name{incomplete_blocks}
\alias{incomplete_blocks}
\title{Generates a Resolvable Incomplete Block Design}
\usage{
incomplete_blocks(
  t = NULL,
  k = NULL,
  r = NULL,
  l = 1,
  plotNumber = 101,
  locationNames = NULL,
  seed = NULL,
  data = NULL
)
}
\arguments{
\item{t}{Number of  treatments.}

\item{k}{Size of incomplete blocks (number of units per incomplete block).}

\item{r}{Number of full blocks (or resolvable replicates) (also number of replicates per treatment).}

\item{l}{Number of locations. By default \code{l = 1}.}

\item{plotNumber}{Numeric vector with the starting plot number for each location. By default \code{plotNumber = 101}.}

\item{locationNames}{(optional) Names for each location.}

\item{seed}{(optional) Real number that specifies the starting seed to obtain reproducible designs.}

\item{data}{(optional) Data frame with label list of treatments.}
}
\value{
A list with two elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{fieldBook} is a data frame with the incomplete block design field book.
}
}
\description{
Randomly generates a resolvable incomplete block design (IBD) of characteristics (t, k, r).
The randomization can be done across locations.
}
\examples{
# Example 1: Generates a resolvable IBD of characteristics (t,k,r) = (12,4,2).
# 1-resolvable IBDs
ibd1 <- incomplete_blocks(t = 12,
                          k = 4,
                          r = 2,
                          seed = 1984)
ibd1$infoDesign
head(ibd1$fieldBook)

# Example 2: Generates a balanced resolvable IBD of characteristics (t,k,r) = (15,3,7).
# In this case, we show how to use the option data.
treatments <- paste("TX-", 1:15, sep = "")
ENTRY <- 1:15
treatment_list <- data.frame(list(ENTRY = ENTRY, TREATMENT = treatments))
head(treatment_list)
ibd2 <- incomplete_blocks(t = 15,
                          k = 3,
                          r = 7,
                          seed = 1985,
                          data = treatment_list)
ibd2$infoDesign
head(ibd2$fieldBook)

}
\references{
Edmondson., R. N. (2021). blocksdesign: Nested and crossed block designs for factorial and
unstructured treatment sets. https://CRAN.R-project.org/package=blocksdesign
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_S3_methods.R
\name{summary.FielDHub}
\alias{summary.FielDHub}
\title{Summary a \code{FielDHub} object}
\usage{
\method{summary}{FielDHub}(object, ...)
}
\arguments{
\item{object}{an object inheriting from class
\code{FielDHub}}

\item{...}{Unused, for extensibility}
}
\value{
an object inheriting from class \code{summary.FielDHub}
}
\description{
Summarise information on the design parameters, and data
  frame structure
}
\examples{
# Example 1: Generates a CRD design with 5 treatments and 5 reps each.
crd1 <- CRD(t = 5, reps = 5, plotNumber = 101,
seed = 1985, locationName = "Fargo")
crd1$infoDesign
summary(crd1)

}
\author{
Thiago de Paula Oliveira,
  \email{thiago.paula.oliveira@alumni.usp.br}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_split_plot.R
\name{split_plot}
\alias{split_plot}
\title{Generates a Split Plot Design}
\usage{
split_plot(
  wp = NULL,
  sp = NULL,
  reps = NULL,
  type = 2,
  l = 1,
  plotNumber = 101,
  seed = NULL,
  locationNames = NULL,
  factorLabels = TRUE,
  data = NULL
)
}
\arguments{
\item{wp}{Number of whole plots, as an integer or a vector.}

\item{sp}{Number of sub plots per whole plot, as an integer or a vector.}

\item{reps}{Number of blocks (full replicates).}

\item{type}{Option for CRD or RCBD designs. Values are \code{type = 1} (CRD) or \code{type = 2} (RCBD). By default \code{type = 2}.}

\item{l}{Number of locations. By default \code{l = 1}.}

\item{plotNumber}{Numeric vector with the starting plot number for each location. By default \code{plotNumber = 101}.}

\item{seed}{(optional) Real number that specifies the starting seed to obtain reproducible designs.}

\item{locationNames}{(optional) Names for each location.}

\item{factorLabels}{(optional) If \code{TRUE} retain the levels
labels from the original data set otherwise, numeric labels will be
assigned. Default is \code{factorLabels =TRUE}.}

\item{data}{(optional) Data frame with label list of treatments.}
}
\value{
A list with two elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{fieldBook} is a data frame with the split plot field book.
}
}
\description{
It randomly generates a split plot design (SPD) across locations.
}
\examples{
# Example 1: Generates a split plot design SPD with 4 whole plots, 2 sub plots per whole plot,
# and 4 reps in an RCBD arrangement. This in for a single location.
SPDExample1 <- split_plot(wp = 4, sp = 2, reps = 5, l = 1, 
                          plotNumber = 101, 
                          seed = 14,
                          type = 2, 
                          locationNames = "FARGO")
SPDExample1$infoDesign
SPDExample1$layoutlocations
head(SPDExample1$fieldBook,12)

# Example 2: Generates a split plot design SPD with 5 whole plots 
# (4 types of fungicide + one control), 10 sub plots per whole plot (10 bean varieties), 
# and 6 reps in an RCBD arrangement. This in 3 locations or sites.
# In this case, we show how to use the option data.
wp <- c("NFung", paste("Fung", 1:4, sep = ""))  # Fungicides (5 Whole plots)
sp <- paste("Beans", 1:10, sep = "")            # Beans varieties (10 sub plots)
split_plot_Data <- data.frame(list(WHOLPLOT = c(wp, rep(NA, 5)), SUBPLOT = sp))
head(split_plot_Data, 12)
SPDExample2 <- split_plot(reps = 6, l = 3, 
                          plotNumber = c(101, 1001, 2001),
                          seed = 23, 
                          type = 2, 
                          locationNames = c("A", "B", "C"),
                          data = split_plot_Data)
SPDExample2$infoDesign
SPDExample2$layoutlocations
head(SPDExample2$fieldBook,12)
             
                  
}
\references{
Federer, W. T. (1955). Experimental Design. Theory and Application. New York, USA. The
Macmillan Company.
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_S3_methods.R
\name{print.summary.FielDHub}
\alias{print.summary.FielDHub}
\title{Print the summary of a \code{FielDHub} object}
\usage{
\method{print}{summary.FielDHub}(x, ...)
}
\arguments{
\item{x}{an object inheriting from class \code{FielDHub}}

\item{...}{Unused, for extensibility}
}
\value{
an object inheriting from class \code{FielDHub}
}
\description{
Print summary information on the design parameters, and
  data frame structure
}
\author{
Thiago de Paula Oliveira,
  \email{thiago.paula.oliveira@alumni.usp.br} [aut],
  Didier Murillo [aut]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_RCBD.R
\name{RCBD}
\alias{RCBD}
\title{Generates a Randomized Complete Block Design (RCBD)}
\usage{
RCBD(
  t = NULL,
  reps = NULL,
  l = 1,
  plotNumber = 101,
  continuous = FALSE,
  planter = "serpentine",
  seed = NULL,
  locationNames = NULL,
  data = NULL
)
}
\arguments{
\item{t}{An integer number with total number of treatments or a vector of dimension t with labels.}

\item{reps}{Number of replicates (full blocks) of each treatment.}

\item{l}{Number of locations. By default \code{l = 1}.}

\item{plotNumber}{Numeric vector with the starting plot number for each location. By default \code{plotNumber = 101}.}

\item{continuous}{Logical value for plot number continuous or not. By default \code{continuous = FALSE}.}

\item{planter}{Option for \code{serpentine} or \code{cartesian} arrangement. By default \code{planter = 'serpentine'}.}

\item{seed}{(optional) Real number that specifies the starting seed to obtain reproducible designs.}

\item{locationNames}{(optional) Names for each location.}

\item{data}{(optional) Data frame with the labels of treatments.}
}
\value{
A list with five elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{layoutRandom} is the RCBD layout randomization for each location.
  \item \code{plotNumber} is the plot number layout for each location.
  \item \code{fieldBook} is a data frame with the RCBD field book design.
}
}
\description{
It randomly generates a randomized complete block design (RCBD) across locations.
}
\examples{
# Example 1: Generates a RCBD design with 3 blocks and 20 treatments across 3 locations.
rcbd1 <- RCBD(t = LETTERS[1:20], reps = 5, l = 3, 
              plotNumber = c(101,1001, 2001), 
              continuous = TRUE,
              planter = "serpentine", 
              seed = 1020, 
              locationNames = c("FARGO", "MINOT", "CASSELTON"))
rcbd1$infoDesign                  
rcbd1$layoutRandom
rcbd1$plotNumber
head(rcbd1$fieldBook)

# Example 2: Generates a RCBD design with 6 blocks and 18 treatments in one location.
# In this case, we show how to use the option data.
treatments <- paste("ND-", 1:18, sep = "")
treatment_list <- data.frame(list(TREATMENT = treatments))
head(treatment_list)
rcbd2 <- RCBD(reps = 6, l = 1, 
              plotNumber = 101, 
              continuous = FALSE, 
              planter = "serpentine", 
              seed = 13, 
              locationNames = "IBAGUE",
              data = treatment_list)
rcbd2$infoDesign                  
rcbd2$layoutRandom
rcbd2$plotNumber
head(rcbd2$fieldBook)


}
\references{
Federer, W. T. (1955). Experimental Design. Theory and Application. New York, USA. The
Macmillan Company.
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_CRD.R
\name{CRD}
\alias{CRD}
\title{Generates a Completely Randomized Design (CRD)}
\usage{
CRD(
  t = NULL,
  reps = NULL,
  plotNumber = 101,
  locationName = NULL,
  seed = NULL,
  data = NULL
)
}
\arguments{
\item{t}{An integer number with total number of treatments or a vector of dimension t with labels.}

\item{reps}{Number of replicates of each treatment.}

\item{plotNumber}{Starting plot number. By default \code{plotNumber = 101}.}

\item{locationName}{(optional) Name of the location.}

\item{seed}{(optional) Real number that specifies the starting seed to obtain reproducible designs.}

\item{data}{(optional) Data frame with the 2 columns with labels of each treatments and its number of replicates.}
}
\value{
A list with two elements.
\itemize{
  \item \code{infoDesign} is a list with information on the design parameters.
  \item \code{fieldBook} is a data frame with the CRD field book.
}
}
\description{
It randomly generates a completely randomized design.
}
\examples{
# Example 1: Generates a CRD design with 10 treatments and 5 reps each.
crd1 <- CRD(t = 10,
            reps = 5,
            plotNumber = 101,
            seed = 1987,
            locationName = "Fargo")
crd1$infoDesign
head(crd1$fieldBook,10)

# Example 2: Generates a CRD design with 15 treatments and 6 reps each.
Gens <- paste("Wheat", 1:15, sep = "")
crd2 <- CRD(t = Gens,
            reps = 6,
            plotNumber = 1001,
            seed = 1654,
            locationName = "Fargo")
crd2$infoDesign
head(crd2$fieldBook,10)

# Example 3: Generates a CRD design with 12 treatments and 4 reps each.
# In this case, we show how to use the option data.
treatments <- paste("ND-", 1:12, sep = "")
treatment_list <- data.frame(list(TREATMENT = treatments, REP = 4))
head(treatment_list)
crd3 <- CRD(t = NULL,
            reps = NULL,
            plotNumber = 2001,
            seed = 1655,
            locationName = "Cali",
            data = treatment_list)
crd3$infoDesign
head(crd3$fieldBook,10)

}
\references{
Federer, W. T. (1955). Experimental Design. Theory and Application. New York, USA. The
Macmillan Company.
}
\author{
Didier Murillo [aut],
        Salvador Gezan [aut],
        Ana Heilman [ctb],
        Thomas Walk [ctb], 
        Johan Aparicio [ctb], 
        Richard Horsley [ctb]
}
