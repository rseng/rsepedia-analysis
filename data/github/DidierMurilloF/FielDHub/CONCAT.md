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
