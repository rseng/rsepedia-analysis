
<!-- README.md is generated from README.Rmd. Please edit that file -->

# oolong <img src="man/figures/oolong_logo.svg" align="right" height="200" />

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/chainsawriot/oolong.svg?branch=master)](https://travis-ci.org/chainsawriot/oolong)
[![Codecov test
coverage](https://codecov.io/gh/chainsawriot/oolong/branch/master/graph/badge.svg)](https://codecov.io/gh/chainsawriot/oolong?branch=master)
[![joss
stataus](https://joss.theoj.org/papers/6e535564e7142d705f4f3d68b18dac62/status.svg)](https://joss.theoj.org/papers/6e535564e7142d705f4f3d68b18dac62)
[![CRAN
status](https://www.r-pkg.org/badges/version/oolong)](https://CRAN.R-project.org/package=oolong)
[![R-CMD-check](https://github.com/chainsawriot/oolong/workflows/R-CMD-check/badge.svg)](https://github.com/chainsawriot/oolong/actions)
<!-- badges: end -->

<img src="man/figures/oolong_demo.gif" align="center" height="400" />

The goal of oolong \[1\] is to generate and administrate validation
tests easily for typical automated content analysis tools such as topic
models and dictionary-based tools.

Please refer to the [overview](overview_gh.md) for an introduction to
this package. If you need to deploy the test online, please refer to the
[Deployment Vignette](deploy_gh.md). If you use BTM, please refer to the
[BTM Vignette](btm_gh.md).

## Citation

Please cite this package as:

Chan C-h. & Sältzer M., (2020). oolong: An R package for validating
automated content analysis tools. Journal of Open Source Software,
5(55), 2461, <https://doi.org/10.21105/joss.02461>

For a BibTeX entry, use the output from `citation(package = "oolong")`.

## Contributing

Contributions in the form of feedback, comments, code, and bug report
are welcome.

  - Fork the source code, modify, and issue a [pull
    request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork).
  - Issues, bug reports: [File a Github
    issue](https://github.com/chainsawriot/oolong).
  - Github is not your thing? Contact Chung-hong Chan by e-mail, post,
    or other methods listed on this
    [page](https://www.mzes.uni-mannheim.de/d7/en/profiles/chung-hong-chan).

## Code of Conduct

Please note that the oolong project is released with a [Contributor Code
of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

-----

1.  /ˈuːlʊŋ/ 烏龍, literally means “Dark Dragon”, is a semi-oxidized tea
    from Asia. It is very popular in Taiwan, Japan and Hong Kong. In
    Cantonese and Taiwanese Mandarin, the same word can also mean
    “confused”. It perfectly captures the spirit of human-in-the-loop
    validation.
# oolong 0.4.1

* Eliminate `miniUI` as a dependency.
* Update the documentation to reflect newly published papers, e.g. Ying et al.

# oolong 0.4.0

* Add `export_oolong` and `deploy_oolong` for online deployment [thanks Marius Sältzer, Daniel Braby (and his friend Louis), Johannes Gruber and Felicia Loecherbach for testing this feature; thanks SAGE Ocean for the concept grant to support the development of this feature]
* Support models from `seededlda` [thanks Marius Sältzer]
* Support Naive Bayes models from `quanteda.textmodels` [thanks Marius Sältzer]
* Support generation of word set intrusion test (Ying et al. forthcoming)
* Support generation of oolong object with only topic intrusion test
* Add new wrappers: `wi`, `ti`, `witi`, `wsi`, and `gs`
* Add `userid` as an suggested parameter
* Total revamp of the object of all oolong tests; add more meta data. Add `update_oolong` for updating object created by older versions of oolong
* Update the print method of all oolong tests; it is now based on `cli`
* Various bug fixes; all Shiny components are now automatically tested

# oolong 0.3.11

* Support BTM [thanks Marius Sältzer]
* Update Shiny UI (with jump button)
* Various bug fixes

# oolong 0.3.4

* Initial CRAN version.
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
reported to the community leaders responsible for enforcement at [INSERT CONTACT
METHOD]. All complaints will be reviewed and investigated promptly and fairly.

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
version 2.0,
available at https://www.contributor-covenant.org/version/2/0/
code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at https://
www.contributor-covenant.org/translations.
## Submission v0.4.1 (Nov 8 2021)

This is a maintenance release

Deploy
================
Chung-hong Chan

In oolong 0.3.22, functions for deploying oolong tests were added
(`export_oolong`, `revert_oolong` etc.). These functions make it
possible for the coders to conduct validation tests online using their
browser, rather than having to install R on their computer.

The basic workflow is simple: 1) create the oolong test object as usual;
2) deploy the test online and obtain the URL to the test; 3) ask your
coders to conduct the test online and send back the data file; 4) revert
back from the data file to an oolong object.

# Create an oolong test

Please note that one cannot deploy oolong test objects with *both* word
and topic intrusion tests, i.e. those created using `witi()` online. If
you need to do both tests, you need to deploy them as two separate
instances: one created using `wi()` and another created using `ti()`.

In this guide, we assume you want to deploy a word set intrusion test
online.

``` r
library(oolong)
wsi_test <- wsi(abstracts_keyatm)
wsi_test
#> 
#> ── oolong (topic model) ────────────────────────────────────────────────────────
#> ✖ WI ✖ TI ✔ WSI
#> ℹ WSI: n = 10, 0 coded.
#> 
#> ── Methods ──
#> 
#> • <$do_word_set_intrusion_test()>: do word set intrusion test
#> • <$lock()>: finalize and see the results
```

# Deploy the test online

First, you need to export the oolong test object as a stand alone Shiny
app. This stand alone Shiny app will be in a directory.

``` r
export_oolong(wsi_test, dir = "./wsi_test", use_full_path = FALSE)
#> ℹ The Shiny has been written to the directory: ./wsi_test
#> ℹ You can test the app with: shiny::runApp("./wsi_test")
```

The directory has only two files

``` r
fs::dir_tree("./wsi_test")
#> ./wsi_test
#> ├── app.R
#> └── oolong.RDS
```

This structure is called [“Single-file Shiny
app.”](https://shiny.rstudio.com/articles/app-formats.html)
Experienced Shiny users might have their preferred method of deploying
this app to whatever Shiny server they can master.

For less experienced users, the simplest way to deploy this app online
is to use [shinyapps.io](https://www.shinyapps.io/) (free tier available
with 25 hours of computational time per month). Please register for an
account at shinyapps.io and configure rsconnect. Please refer to [this
guide](https://shiny.rstudio.com/articles/shinyapps.html) for more
information. Please remember to configure the tokens.

``` r
## replace <ACCOUNT>, <TOKEN>, <SECRET> with the information from your profile on Shinyapps.io: click Your name -> Tokens
rsconnect::setAccountInfo(name="<ACCOUNT>", token="<TOKEN>", secret="<SECRET>")
```

For RStudio users, the simplest way to deploy the app to shinyapps.io is
to first launch the app.

``` r
library(shiny)
runApp("./wsi_test")
```

And then click the **Publish** button at the right corner of the
launched window.

You will be asked for the title of the app, just give it a name,
e.g. *wsi\_test*. You probably can keep other default settings and push
the **Publish** button to initialize the deployment process.

<img src="vignettes/figures/deploying_shinyappsio.png" align="center" height="400" />

If there is no hiccup, you will get a URL to your deployed oolong test.
Something like: *<https://yourname.shinyapps.io/wsi_test/>*

# Conduct the test

You can give the URL to your coders and they conduct the test with their
browser online. The only difference of the deployed version is that,
there will be a userid prompt and download button after the coding.

<img src="vignettes/figures/oolong_download.png" align="center"/>

You should instruct your coders to download the data file after coding
and return it to you. \[1\]

# Revert

You can then obtain a locked oolong object from the original oolong and
the downloaded data file. `revert_oolong` will do verifications with the
original oolong object to make sure no error and no cheating.

``` r
revert_oolong(wsi_test, "oolong_2021-05-22 20 51 26 Hadley Wickham.RDS")
```

    #> 
    #> ── oolong (topic model) ────────────────────────────────────────────────────────
    #> ✖ WI ✖ TI ✔ WSI
    #> ☺ Hadley Wickham
    #> ℹ WSI: n = 10, 10 coded.
    #> 
    #> ── Results: ──
    #> 
    #> ℹ 80%  precision (WSI)

1.  Future versions might provide permanent storage
Overview
================
Chung-hong Chan

The validation test is called “oolong test” (for reading tea leaves).
This package provides several functions for generating different types
of oolong test.

| function | purpose                                                                                                                           |
| -------: | :-------------------------------------------------------------------------------------------------------------------------------- |
|   `wi()` | validating a topic model with [word intrusion test](#word-intrusion-test) (Chang et al., 2008)                                    |
|   `ti()` | validating a topic model with [topic intrusion test](#topic-intrusion-test) (Chang et al., 2008; aka “T8WSI” in Ying et al. 2021) |
| `witi()` | validating a topic model with [word intrusion test](#word-intrusion-test) and [topic intrusion test](#topic-intrusion-test)       |
|  `wsi()` | validating a topic model with [word set intrusion test](#word-set-intrusion-test) (Ying et al. 2021)                              |
|   `gs()` | oolong test for [creating gold standard](#creating-gold-standard) (see Song et al., 2020)                                         |

All of these tests can also be generated with the function
[`create_oolong`](#backward-compatibility). As of version 0.3.20, it is
no longer recommended.

## Installation

Because the package is constantly changing, we suggest using the
development version from [GitHub](https://github.com/):

``` r
# install.packages("devtools")
devtools::install_github("chainsawriot/oolong")
```

You can also install the “stable” (but slightly older) version from
CRAN:

``` r
install.packages("oolong")
```

## Validating Topic Models

#### Word intrusion test

`abstracts_keyatm` is an example topic model trained with the data
`abstracts` using the `keyATM` package. Currently, this package supports
structural topic models / correlated topic models from `stm`, Warp LDA
models from `text2vec` , LDA/CTM models from `topicmodels`, Biterm Topic
Models from `BTM`, Keyword Assisted Topic Models from `keyATM`, and
seeded LDA models from `seededlda`. Although not strictly a topic model,
Naive Bayes models from `quanteda.textmodels` are also supported. See
the section on [Naive Bayes](#about-naive-bayes) for more information.

``` r
library(oolong)
library(keyATM)
#> keyATM 0.4.0 successfully loaded.
#>  Papers, examples, resources, and other materials are at
#>  https://keyatm.github.io/keyATM/
library(quanteda)
#> Package version: 3.2.0
#> Unicode version: 13.0
#> ICU version: 66.1
#> Parallel computing: 8 of 8 threads used.
#> See https://quanteda.io for tutorials and examples.
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

``` r
abstracts_keyatm
#> keyATM_output object for the base model.
```

To create an oolong test with word intrusion test, use the function
`wi`. It is recommended to provide a user id of coder who are going to
be doing the test.

``` r
oolong_test <- wi(abstracts_keyatm, userid = "Hadley")
oolong_test
#> 
#> ── oolong (topic model) ────────────────────────────────────────────────────────
#> ✔ WI ✖ TI ✖ WSI
#> ☺ Hadley
#> ℹ WI: k = 10, 0 coded.
#> 
#> ── Methods ──
#> 
#> • <$do_word_intrusion_test()>: do word intrusion test
#> • <$lock()>: finalize and see the results
```

As instructed, use the method `$do_word_intrusion_test()` to start
coding.

``` r
oolong_test$do_word_intrusion_test()
```

After the coding, you need to first lock the test. Then, you can look at
the model precision by printing the oolong test.

``` r
oolong_test$lock()
oolong_test
#> 
#> ── oolong (topic model) ────────────────────────────────────────────────────────
#> ✔ WI ✖ TI ✖ WSI
#> ☺ Hadley
#> ℹ WI: k = 10, 10 coded.
#> 
#> ── Results: ──
#> 
#> ℹ 90%  precision
```

#### Word set intrusion test

Word set intrusion test is a variant of word intrusion test (Ying et
al., 2021), in which multiple word sets generated from top terms of one
topic are juxtaposed with one intruder word set generated similarly from
another topic. In Ying et al., this test is called “R4WSI” because 4
word sets are displayed. By default, oolong generates also R4WSI.
However, it is also possible to generate R(N)WSI by setting the
parameter `n_correct_ws` to N - 1.

``` r
oolong_test <- wsi(abstracts_keyatm, userid = "Garrett")
oolong_test
#> 
#> ── oolong (topic model) ────────────────────────────────────────────────────────
#> ✖ WI ✖ TI ✔ WSI
#> ☺ Garrett
#> ℹ WSI: n = 10, 0 coded.
#> 
#> ── Methods ──
#> 
#> • <$do_word_set_intrusion_test()>: do word set intrusion test
#> • <$lock()>: finalize and see the results
```

Use the method `$do_word_set_intrusion_test()` to start coding.

``` r
oolong_test$do_word_set_intrusion_test()
```

``` r
oolong_test$lock()
oolong_test
#> 
#> ── oolong (topic model) ────────────────────────────────────────────────────────
#> ✖ WI ✖ TI ✔ WSI
#> ☺ Garrett
#> ℹ WSI: n = 10, 10 coded.
#> 
#> ── Results: ──
#> 
#> ℹ 90%  precision (WSI)
```

#### Topic intrusion test

For example, `abstracts_keyatm` was generated with the corpus
`abstracts$text`

``` r
library(tibble)
abstracts
#> # A tibble: 2,500 × 1
#>    text                                                                         
#>    <chr>                                                                        
#>  1 This study explores the benefits and risks featured in medical tourism broke…
#>  2 This article puts forth the argument that with the transfer of stock trading…
#>  3 The purpose of this study was to evaluate the effect the visual fidelity of …
#>  4 Among the many health issues relevant to college students, overconsumption o…
#>  5 This address, delivered at ICA's 50th anniversary conference, calls on the a…
#>  6 The Internet has often been used to reach men who have sex with men (MSMs) i…
#>  7 This article argues that the literature describing the internet revolution i…
#>  8 This research study examined Bud Goodall's online health narrative as a case…
#>  9 Information technology and new media allow for collecting and sharing person…
#> 10 Using a national, telephone survey of 1,762 adolescents aged 12-17 years, th…
#> # … with 2,490 more rows
```

Creating the oolong test object with the corpus used for training the
topic model will generate topic intrusion test cases.

``` r
oolong_test <- ti(abstracts_keyatm, abstracts$text, userid = "Julia")
oolong_test
#> 
#> ── oolong (topic model) ────────────────────────────────────────────────────────
#> ✖ WI ✔ TI ✖ WSI
#> ☺ Julia
#> ℹ TI: n = 25, 0 coded.
#> 
#> ── Methods ──
#> 
#> • <$do_topic_intrusion_test()>: do topic intrusion test
#> • <$lock()>: finalize and see the results
```

Similarly, use the `$do_topic_intrusion_test` to code the test cases,
lock the test with `$lock()` and then you can look at the TLO (topic log
odds) value by printing the oolong test.

``` r
oolong_test$do_topic_intrusion_test()
oolong_test$lock()
```

``` r
oolong_test
#> 
#> ── oolong (topic model) ────────────────────────────────────────────────────────
#> ✖ WI ✔ TI ✖ WSI
#> ☺ Julia
#> ℹ TI: n = 25, 25 coded.
#> 
#> ── Results: ──
#> 
#> ℹ TLO: -0.009
```

### Suggested workflow

The test makes more sense if more than one coder is involved. A
suggested workflow is to create the test, then clone the oolong object.
Ask multiple coders to do the test(s) and then summarize the results.

Preprocess and create a document-feature matrix

``` r
dfm(abstracts$text, tolower = TRUE, stem = TRUE, remove = stopwords('english'), remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE, remove_hyphens = TRUE) %>% dfm_trim(min_docfreq = 3, max_docfreq = 500) %>% dfm_select(min_nchar = 3, pattern = "^[a-zA-Z]+$", valuetype = "regex") -> abstracts_dfm
```

Train a topic model.

``` r
require(keyATM)
abstracts_keyatm <- keyATM(keyATM_read(abstracts_dfm), no_keyword_topics = 0, keywords = abstracts_dictionary, model = "base", options = list(seed = 46709394))
```

Create a new oolong object.

``` r
oolong_test_rater1 <- witi(abstracts_keyatm, abstracts$text, userid = "Yihui")
```

Clone the oolong object to be used by other raters.

``` r
oolong_test_rater2 <- clone_oolong(oolong_test_rater1, userid = "Jenny")
```

Ask different coders to code each object and then lock the object.

``` r
## Let Yihui do the test.
oolong_test_rater1$do_word_intrusion_test()
oolong_test_rater1$do_topic_intrusion_test()
oolong_test_rater1$lock()

## Let Jenny do the test.
oolong_test_rater2$do_word_intrusion_test()
oolong_test_rater2$do_topic_intrusion_test()
oolong_test_rater2$lock()
```

Get a summary of the two objects.

``` r
summarize_oolong(oolong_test_rater1, oolong_test_rater2)
#> 
#> ── Summary (topic model): ──────────────────────────────────────────────────────
#> 
#> ── Word intrusion test ──
#> 
#> ℹ Mean model precision: 0.3
#> ℹ Quantiles of model precision: 0.1, 0.2, 0.3, 0.4, 0.5
#> ℹ P-value of the model precision
#> (H0: Model precision is not better than random guess): 0.0693
#> ℹ Krippendorff's alpha: 0.095
#> ℹ K Precision:
#> 0, 0.5, 1, 0, 0, 0.5, 0, 0.5, 0, 0.5
#> 
#> ── Topic intrusion test ──
#> 
#> ℹ Mean TLO: -3.2
#> ℹ Median TLO: -4.07
#> ℹ Quantiles of TLO: -8.26, -6.09, -4.07, 0, 0
#> ℹ P-Value of the median TLO 
#> (H0: Median TLO is not better than random guess): 0.096
```

### About the p-values

The test for model precision (MP) is based on an one-tailed, one-sample
binomial test for each rater. In a multiple-rater situation, the
p-values from all raters are combined using the Fisher’s method (a.k.a.
Fisher’s omnibus test).

H0: MP is not better than 1/ n\_top\_terms

H1: MP is better than 1/ n\_top\_terms

The test for the median of TLO is based on a permutation test.

H0: Median TLO is not better than random guess.

H1: Median TLO is better than random guess.

One must notice that the two statistical tests are testing the bear
minimum. A significant test only indicates the topic model can make the
rater(s) perform better than random guess. It is not an indication of
good topic interpretability. Also, one should use a very conservative
significant level, e.g. \(\alpha < 0.001\).

### About Warp LDA

There is a subtle difference between the support for `stm` and for
`text2vec`.

`abstracts_warplda` is a Warp LDA object trained with the same dataset
as the `abstracts_stm`

``` r
abstracts_warplda
#> <WarpLDA>
#>   Inherits from: <LDA>
#>   Public:
#>     clone: function (deep = FALSE) 
#>     components: 0 1 0 46 0 95 0 20 42 8 31 36 50 23 0 0 0 58 0 43 0 0 0  ...
#>     fit_transform: function (x, n_iter = 1000, convergence_tol = 0.001, n_check_convergence = 10, 
#>     get_top_words: function (n = 10, topic_number = 1L:private$n_topics, lambda = 1) 
#>     initialize: function (n_topics = 10L, doc_topic_prior = 50/n_topics, topic_word_prior = 1/n_topics, 
#>     plot: function (lambda.step = 0.1, reorder.topics = FALSE, doc_len = private$doc_len, 
#>     topic_word_distribution: 0 9.41796948577887e-05 0 0.00446992517733942 0 0.0086837 ...
#>     transform: function (x, n_iter = 1000, convergence_tol = 0.001, n_check_convergence = 10, 
#>   Private:
#>     calc_pseudo_loglikelihood: function (ptr = private$ptr) 
#>     check_convert_input: function (x) 
#>     components_: 0 1 0 46 0 95 0 20 42 8 31 36 50 23 0 0 0 58 0 43 0 0 0  ...
#>     doc_len: 80 68 85 88 69 118 99 50 57 88 70 67 53 62 66 92 89 79 1 ...
#>     doc_topic_distribution: function () 
#>     doc_topic_distribution_with_prior: function () 
#>     doc_topic_matrix: 0 0 0 0 0 3 111 0 0 0 0 0 90 134 0 174 0 321 0 0 109 38  ...
#>     doc_topic_prior: 0.1
#>     fit_transform_internal: function (model_ptr, n_iter, convergence_tol, n_check_convergence, 
#>     get_c_all: function () 
#>     get_c_all_local: function () 
#>     get_doc_topic_matrix: function (prt, nr) 
#>     get_topic_word_count: function () 
#>     init_model_dtm: function (x, ptr = private$ptr) 
#>     internal_matrix_formats: list
#>     is_initialized: FALSE
#>     n_iter_inference: 10
#>     n_topics: 20
#>     ptr: externalptr
#>     reset_c_local: function () 
#>     run_iter_doc: function (update_topics = TRUE, ptr = private$ptr) 
#>     run_iter_word: function (update_topics = TRUE, ptr = private$ptr) 
#>     seeds: 135203513.874082 471172603.061186
#>     set_c_all: function (x) 
#>     set_internal_matrix_formats: function (sparse = NULL, dense = NULL) 
#>     topic_word_distribution_with_prior: function () 
#>     topic_word_prior: 0.01
#>     transform_internal: function (x, n_iter = 1000, convergence_tol = 0.001, n_check_convergence = 10, 
#>     vocabulary: explor benefit risk featur medic broker websit well type ...
```

All the API endpoints are the same, except the one for the creation of
topic intrusion test cases. You must supply also the `input_dfm`.

``` r
### Just word intrusion test.
oolong_test <- wi(abstracts_warplda, userid = "Lionel")
oolong_test
#> 
#> ── oolong (topic model) ────────────────────────────────────────────────────────
#> ✔ WI ✖ TI ✖ WSI
#> ☺ Lionel
#> ℹ WI: k = 20, 0 coded.
#> 
#> ── Methods ──
#> 
#> • <$do_word_intrusion_test()>: do word intrusion test
#> • <$lock()>: finalize and see the results
```

``` r
abstracts_dfm
#> Document-feature matrix of: 2,500 documents, 3,998 features (98.61% sparse) and 0 docvars.
#>        features
#> docs    explor benefit risk featur medic broker websit well type persuas
#>   text1      1       2    2      2     6      3      6    1    3       1
#>   text2      0       0    1      0     0      0      0    0    1       0
#>   text3      0       1    0      0     0      0      0    0    0       0
#>   text4      1       0    0      0     0      0      0    0    0       0
#>   text5      1       0    0      0     0      0      0    0    0       0
#>   text6      0       1    1      0     0      0      0    0    0       0
#> [ reached max_ndoc ... 2,494 more documents, reached max_nfeat ... 3,988 more features ]
```

``` r
oolong_test <- witi(abstracts_warplda, abstracts$text, input_dfm = abstracts_dfm, userid = "Mara")
```

``` r
oolong_test
#> 
#> ── oolong (topic model) ────────────────────────────────────────────────────────
#> ✔ WI ✔ TI ✖ WSI
#> ☺ Mara
#> ℹ WI: k = 20, 0 coded.
#> ℹ TI: n = 25, 0 coded.
#> 
#> ── Methods ──
#> 
#> • <$do_word_intrusion_test()>: do word intrusion test
#> • <$do_topic_intrusion_test()>: do topic intrusion test
#> • <$lock()>: finalize and see the results
```

## About Biterm Topic Model

Please refer to the vignette about BTM.

## About Naive Bayes

Naive Bayes model is a supervised machine learning model. This package
supports Naive Bayes models trained using `quanteda.textmodels`.

Suppose `newsgroup_nb` is a Naive Bayes model trained on a subset of the
classic \[20 newsgroups\] dataset.

``` r
tokens(newsgroup5$text, remove_punct = TRUE, remove_symbols = TRUE, remove_numbers = TRUE, remove_url = TRUE, spilit_hyphens = TRUE) %>% tokens_wordstem %>% tokens_remove(stopwords("en")) %>% dfm(tolower = TRUE) %>% dfm_trim(min_termfreq = 3, max_docfreq = 0.06, docfreq_type = "prop") -> newsgroup_dfm
docvars(newsgroup_dfm, "group") <- newsgroup5$title
newsgroup_nb <- textmodel_nb(newsgroup_dfm, docvars(newsgroup_dfm, "group"), distribution = "Bernoulli")
```

You can still generate word intrusion and word set intrusion tests.

``` r
wi(newsgroup_nb)
#> 
#> ── oolong (topic model) ────────────────────────────────────────────────────────
#> ✔ WI ✖ TI ✖ WSI
#> ℹ WI: k = 20, 0 coded.
#> 
#> ── Methods ──
#> 
#> • <$do_word_intrusion_test()>: do word intrusion test
#> • <$lock()>: finalize and see the results
```

``` r
wsi(newsgroup_nb)
#> 
#> ── oolong (topic model) ────────────────────────────────────────────────────────
#> ✖ WI ✖ TI ✔ WSI
#> ℹ WSI: n = 20, 0 coded.
#> 
#> ── Methods ──
#> 
#> • <$do_word_set_intrusion_test()>: do word set intrusion test
#> • <$lock()>: finalize and see the results
```

## Validating Dictionary-based Methods

### Creating gold standard

`trump2k` is a dataset of 2,000 tweets from @realdonaldtrump.

``` r
tibble(text = trump2k)
#> # A tibble: 2,000 × 1
#>    text                                                                         
#>    <chr>                                                                        
#>  1 "In just out book, Secret Service Agent Gary Byrne doesn't believe that Croo…
#>  2 "Hillary Clinton has announced that she is letting her husband out to campai…
#>  3 "\"@TheBrodyFile: Always great to visit with @TheBrodyFile one-on-one with \…
#>  4 "Explain to @brithume and @megynkelly, who know nothing, that I will beat Hi…
#>  5 "Nobody beats me on National Security. https://t.co/sCrj4Ha1I5"              
#>  6 "\"@realbill2016: @realDonaldTrump @Brainykid2010 @shl Trump leading LA Time…
#>  7 "\"@teapartynews: Trump Wins Tea Party Group's 'Nashville Straw Poll' - News…
#>  8 "Big Republican Dinner tonight at Mar-a-Lago in Palm Beach. I will be there!"
#>  9 ".@HillaryClinton loves to lie. America has had enough of the CLINTON'S! It …
#> 10 "\"@brianstoya: @realDonaldTrump For POTUS #2016\""                          
#> # … with 1,990 more rows
```

For example, you are interested in studying the sentiment of these
tweets. One can use tools such as AFINN to automatically extract
sentiment in these tweets. However, oolong recommends to generate gold
standard by human coding first using a subset. By default, oolong
selects 1% of the origin corpus as test cases. The parameter `construct`
should be an adjective, e.g. positive, liberal, populistic, etc.

``` r
oolong_test <- gs(input_corpus = trump2k, construct = "positive", userid = "Joe")
oolong_test
#> 
#> ── oolong (gold standard generation) ───────────────────────────────────────────
#> ☺ Joe
#> ℹ GS: n = 20, 0 coded.
#> ℹ Construct:  positive.
#> 
#> ── Methods ──
#> 
#> • <$do_gold_standard_test()>: generate gold standard
#> • <$lock()>: finalize this object and see the results
```

As instructed, use the method `$do_gold_standard_test()` to start
coding.

``` r
oolong_test$do_gold_standard_test()
```

After the coding, you need to first lock the test and then the
`$turn_gold()` method is available.

``` r
oolong_test$lock()
oolong_test
#> 
#> ── oolong (gold standard generation) ───────────────────────────────────────────
#> ☺ Joe
#> ℹ GS: n = 20, 20 coded.
#> ℹ Construct:  positive.
#> 
#> ── Methods ──
#> 
#> • <$turn_gold()>: convert the test results into a quanteda corpus
```

### Example: Validating AFINN using the gold standard

A locked oolong test can be converted into a quanteda-compatible corpus
for further analysis. The corpus contains two `docvars`, ‘answer’.

``` r
oolong_test$turn_gold()
#> Corpus consisting of 20 documents and 1 docvar.
#> text1 :
#> "Thank you Eau Claire, Wisconsin.  #VoteTrump on Tuesday, Apr..."
#> 
#> text2 :
#> ""@bobby990r_1: @realDonaldTrump would lead polls the second ..."
#> 
#> text3 :
#> ""@KdanielsK: @misstcassidy @AllAboutTheTea_ @realDonaldTrump..."
#> 
#> text4 :
#> "Thank you for a great afternoon Birmingham, Alabama! #Trump2..."
#> 
#> text5 :
#> ""@THETAINTEDT: @foxandfriends @realDonaldTrump Trump 2016 ht..."
#> 
#> text6 :
#> "People believe CNN these days almost as little as they belie..."
#> 
#> [ reached max_ndoc ... 14 more documents ]
#> ℹ Access the answer from the coding with quanteda::docvars(obj, 'answer')
```

In this example, we calculate the AFINN score for each tweet using
quanteda. The dictionary `afinn` is bundle with this package.

``` r
gold_standard <- oolong_test$turn_gold()
dfm(gold_standard, remove_punct = TRUE) %>% dfm_lookup(afinn) %>% quanteda::convert(to = "data.frame") %>%
    mutate(matching_word_valence = (neg5 * -5) + (neg4 * -4) + (neg3 * -3) + (neg2 * -2) + (neg1 * -1)
           + (zero * 0) + (pos1 * 1) + (pos2 * 2) + (pos3 * 3) + (pos4 * 4) + (pos5 * 5),
           base = ntoken(gold_standard, remove_punct = TRUE), afinn_score = matching_word_valence / base) %>%
    pull(afinn_score) -> all_afinn_score
#> Warning: 'dfm.corpus()' is deprecated. Use 'tokens()' first.
#> Warning: '...' should not be used for tokens() arguments; use 'tokens()' first.
all_afinn_score
#>       text1       text2       text3       text4       text5       text6 
#>  0.33333333 -0.09090909 -0.16666667  0.45454545  0.00000000  0.00000000 
#>       text7       text8       text9      text10      text11      text12 
#>  0.16666667  0.38461538  0.00000000  0.38461538 -0.29166667  0.00000000 
#>      text13      text14      text15      text16      text17      text18 
#>  0.50000000  0.07142857  0.00000000 -0.12000000  0.28571429  0.16000000 
#>      text19      text20 
#>  0.36842105  0.38888889
```

Put back the vector of AFINN score into the respective `docvars` and
study the correlation between the gold standard and AFINN.

``` r
summarize_oolong(oolong_test, target_value = all_afinn_score)
#> New names:
#> * NA -> ...1
#> `geom_smooth()` using formula 'y ~ x'
#> `geom_smooth()` using formula 'y ~ x'
#> 
#> ── Summary (gold standard generation): ─────────────────────────────────────────
#> ℹ Correlation: 0.718 (p = 4e-04)
#> ℹ Effect of content length: -0.323 (p = 0.1643)
```

### Suggested workflow

Create an oolong object, clone it for another coder. According to Song
et al. (Forthcoming), you should at least draw 1% of your data.

``` r
trump <- gs(input_corpus = trump2k, exact_n = 40, userid = "JJ")
trump2 <- clone_oolong(trump, userid = "Winston")
```

Instruct two coders to code the tweets and lock the objects.

``` r
trump$do_gold_standard_test()
trump2$do_gold_standard_test()
trump$lock()
trump2$lock()
```

Calculate the target value (in this case, the AFINN score) by turning
one object into a corpus.

``` r
gold_standard <- trump$turn_gold()
dfm(gold_standard, remove_punct = TRUE) %>% dfm_lookup(afinn) %>% quanteda::convert(to = "data.frame") %>%
    mutate(matching_word_valence = (neg5 * -5) + (neg4 * -4) + (neg3 * -3) + (neg2 * -2) + (neg1 * -1)
           + (zero * 0) + (pos1 * 1) + (pos2 * 2) + (pos3 * 3) + (pos4 * 4) + (pos5 * 5),
           base = ntoken(gold_standard, remove_punct = TRUE), afinn_score = matching_word_valence / base) %>%
    pull(afinn_score) -> target_value
#> Warning: 'dfm.corpus()' is deprecated. Use 'tokens()' first.
#> Warning: '...' should not be used for tokens() arguments; use 'tokens()' first.
```

Summarize all oolong objects with the target value.

``` r
res <- summarize_oolong(trump, trump2, target_value = target_value)
#> New names:
#> * NA -> ...1
#> * NA -> ...2
#> `geom_smooth()` using formula 'y ~ x'
#> `geom_smooth()` using formula 'y ~ x'
```

Read the results. The diagnostic plot consists of 4 subplots. It is a
good idea to read Bland & Altman (1986) on the difference between
correlation and agreement.

  - Subplot (top left): Raw correlation between human judgement and
    target value. One should want to have a good correlation between the
    two.
  - Subplot (top right): Bland-Altman plot. One should want to have no
    correlation. Also, the dots should be randomly scattering around the
    mean value. If it is so, the two measurements (human judgement and
    target value) are in good agreement.
  - Subplot (bottom left): Raw correlation between target value and
    content length. One should want to have no correlation, as an
    indication of good reliability against the influence of content
    length. (See Chan et al. 2020)
  - Subplot (bottom right): Cook’s distance of all data point. One
    should want to have no dot (or at least very few dots) above the
    threshold. It is an indication of how the raw correlation between
    human judgement and target value can or cannot be influenced by
    extreme values in your data.

The textual output contains the Krippendorff’s alpha of the codings by
your raters. In order to claim validity of your target value, you must
first establish the reliability of your gold standard. Song et
al. \[Forthcoming\] suggest Krippendorff’s Alpha \> 0.7 as an
acceptable cut-off.

``` r
res
#> 
#> ── Summary (gold standard generation): ─────────────────────────────────────────
#> ℹ Krippendorff's Alpha: 0.931
#> ℹ Correlation: 0.744 (p = 2e-04)
#> ℹ Effect of content length: -0.323 (p = 0.1643)
```

``` r
plot(res)
```

<img src="man/figures/README-diagnosis-1.png" width="100%" />

## Backward compatibility

Historically, oolong test objects could only be generated with only one
function: `create_oolong`. It is no longer the case and no longer
recommended anymore. It is still retained for backward compatibility
purposes. If you still need to use `create_oolong()`, the most important
parameters are `input_model` and `input_corpus`. Setting each of them to
`NULL` generates different tests.

| input\_model | input\_corpus | output                                                                                                                                      |
| ------------ | :-----------: | ------------------------------------------------------------------------------------------------------------------------------------------- |
| Not NULL     |     NULL      | oolong test for validating a topic model with [word intrusion test](#word-intrusion-test)                                                   |
| Not NULL     |   Not NULL    | oolong test for validating a topic model with [word intrusion test](#word-intrusion-test) and [topic intrusion test](#topic-intrusion-test) |
| NULL         |   Not NULL    | oolong test for [creating gold standard](#creating-gold-standard)                                                                           |
| NULL         |     NULL      | error                                                                                                                                       |

## References

1.  Chang, J., Gerrish, S., Wang, C., Boyd-Graber, J. L., & Blei, D. M.
    (2009). Reading tea leaves: How humans interpret topic models. In
    Advances in neural information processing systems (pp. 288-296).
    [link](https://papers.nips.cc/paper/3700-reading-tea-leaves-how-humans-interpret-topic-models)
2.  Ying, L., Montgomery, J. M., & Stewart, B. M. (2021). Inferring
    concepts from topics: Towards procedures for validating topics as
    measures. Political Analysis.
    [link](https://doi.org/10.1017/pan.2021.33)
3.  Song et al. (2020) In validations we trust? The impact of imperfect
    human annotations as a gold standard on the quality of validation of
    automated content analysis. Political Communication.
    [link](https://doi.org/10.1080/10584609.2020.1723752)
4.  Bland, J. M., & Altman, D. (1986). Statistical methods for assessing
    agreement between two methods of clinical measurement. The lancet,
    327(8476), 307-310.
5.  Chan et al. (2020) Four best practices for measuring news sentiment
    using ‘off-the-shelf’ dictionaries: a large-scale p-hacking
    experiment. Computational Communication Research.
    [link](https://osf.io/preprints/socarxiv/np5wa/)
6.  Nielsen, F. Å. (2011). A new ANEW: Evaluation of a word list for
    sentiment analysis in microblogs. arXiv preprint arXiv:1103.2903.
    [link](https://arxiv.org/abs/1103.2903)

-----
BTM
================
Chung-hong Chan

The package BTM by Jan Wijffels et al. finds “topics in collections of
short text”. Compared to other topic model packages, BTM requires a
special data format for training. Oolong has no problem generating word
intrusion tests with BTM. However, that special data format can make
creation of topic intrusion tests very tricky.

This guide provides our recommendations on how to use BTM, so that the
model can be used for generating topic intrusion tests.

# Requirement \#1: Keep your quanteda corpus

It is because every document has a unique document id.

``` r
require(BTM)
#> Loading required package: BTM
require(quanteda)
#> Loading required package: quanteda
#> Package version: 3.2.0
#> Unicode version: 13.0
#> ICU version: 66.1
#> Parallel computing: 8 of 8 threads used.
#> See https://quanteda.io for tutorials and examples.
require(oolong)
#> Loading required package: oolong
trump_corpus <- corpus(trump2k)
```

And then you can do regular text cleaning, stemming procedure with
`quanteda`. Instead of making the product a `DFM` object, make it a
`token` object. You may read [this
issue](https://github.com/quanteda/quanteda/issues/1404) by Benoit et
al.

``` r
tokens(trump_corpus, remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE, split_hyphens = TRUE, remove_url = TRUE) %>% tokens_tolower() %>% tokens_remove(stopwords("en")) %>% tokens_remove("@*")  -> trump_toks
```

# Requirement \#2: Keep your data frame

Use this function to convert the `token` object to a data frame.

``` r
as.data.frame.tokens <- function(x) {
  data.frame(
    doc_id = rep(names(x), lengths(x)),
    tokens = unlist(x, use.names = FALSE)
  )
}

trump_dat <- as.data.frame.tokens(trump_toks)
```

Train a BTM model

``` r
trump_btm <- BTM(trump_dat, k = 8, iter = 500, trace = 10)
```

## Pecularities of BTM

This is how you should generate \(\theta_{t}\) . However, there are many
NaN and there are only 1994 rows (`trump2k` has 2000 tweets) due to
empty documents.

``` r
theta <- predict(trump_btm, newdata = trump_dat)
dim(theta)
#> [1] 1994    8
```

``` r
setdiff(docid(trump_corpus), row.names(theta))
#> [1] "text604"  "text633"  "text659"  "text1586" "text1587" "text1761"
```

``` r
trump_corpus[604]
#> Corpus consisting of 1 document.
#> text604 :
#> "http://t.co/PtViAyrO4A"
```

Also, the row order is messed up.

``` r
head(row.names(theta), 100)
#>   [1] "text1"    "text10"   "text100"  "text1000" "text1001" "text1002"
#>   [7] "text1003" "text1004" "text1005" "text1006" "text1007" "text1008"
#>  [13] "text1009" "text101"  "text1010" "text1011" "text1012" "text1013"
#>  [19] "text1014" "text1015" "text1016" "text1017" "text1018" "text1019"
#>  [25] "text102"  "text1020" "text1021" "text1022" "text1023" "text1024"
#>  [31] "text1025" "text1026" "text1027" "text1028" "text1029" "text103" 
#>  [37] "text1030" "text1031" "text1032" "text1033" "text1034" "text1035"
#>  [43] "text1036" "text1037" "text1038" "text1039" "text104"  "text1040"
#>  [49] "text1041" "text1042" "text1043" "text1044" "text1045" "text1046"
#>  [55] "text1047" "text1048" "text1049" "text105"  "text1050" "text1051"
#>  [61] "text1052" "text1053" "text1054" "text1055" "text1056" "text1057"
#>  [67] "text1058" "text1059" "text106"  "text1060" "text1061" "text1062"
#>  [73] "text1063" "text1064" "text1065" "text1066" "text1067" "text1068"
#>  [79] "text1069" "text107"  "text1070" "text1071" "text1072" "text1073"
#>  [85] "text1074" "text1075" "text1076" "text1077" "text1078" "text1079"
#>  [91] "text108"  "text1080" "text1081" "text1082" "text1083" "text1084"
#>  [97] "text1085" "text1086" "text1087" "text1088"
```

# Oolong’s support for BTM

Oolong has no problem generating word intrusion test for BTM like you do
with other topic models.

``` r
oolong <- create_oolong(trump_btm)
oolong
#> 
#> ── oolong (topic model) ────────────────────────────────────────────────────────
#> ✔ WI ✖ TI ✖ WSI
#> ℹ WI: k = 8, 0 coded.
#> 
#> ── Methods ──
#> 
#> • <$do_word_intrusion_test()>: do word intrusion test
#> • <$lock()>: finalize and see the results
```

For generating topic intrusion tests, however, you must provide the data
frame you used for training (in this case `trump_dat`). Your
`input_corpus` must be a quanteda corpus too.

``` r
oolong <- create_oolong(trump_btm, trump_corpus, btm_dataframe = trump_dat)
oolong
#> 
#> ── oolong (topic model) ────────────────────────────────────────────────────────
#> ✔ WI ✔ TI ✖ WSI
#> ℹ WI: k = 8, 0 coded.
#> ℹ TI: n = 20, 0 coded.
#> 
#> ── Methods ──
#> 
#> • <$do_word_intrusion_test()>: do word intrusion test
#> • <$do_topic_intrusion_test()>: do topic intrusion test
#> • <$lock()>: finalize and see the results
```

`btm_dataframe` must not be NULL.

``` r
oolong <- create_oolong(trump_btm, trump_corpus)
#> Error: You need to provide input_corpus (in quanteda format) and btm_dataframe for generating topic intrusion tests.
```

`input_corpus` must be a quanteda corpus.

``` r
oolong <- create_oolong(trump_btm, trump2k, btm_dataframe = trump_dat)
#> Error: You need to provide input_corpus (in quanteda format) and btm_dataframe for generating topic intrusion tests.
```
---
title: 'oolong: An R package for validating  automated content analysis tools'
tags:
  - R
  - text analysis
  - topic model
  - sentiment analysis
  - validation
authors:
  - name: Chung-hong Chan
    orcid: 0000-0002-6232-7530
    affiliation: 1
  - name: Marius Sältzer
    orcid: 0000-0002-8604-4666
    affiliation: 1
affiliations:
 - name: Mannheimer Zentrum für Europäische Sozialforschung, Universität Mannheim
   index: 1
citation_author: Chan & Sältzer.
date: 15 October 2020
year: 2020
bibliography: paper.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS
---

# Statement of need

Oolong is an R package providing functions for semantic validation of topic modeling and dictionary-based methods, two main tools for doing automated content analysis [@boumans2016taking;@gunther2016word]. 

While the validation of statistical properties of topics models is well established, the substantive meaning of categories uncovered is often less clear and their interpretation reliant on "intuition" or "eyeballing". As @chang2009reading [p. 1] put it: "qualitative evaluation of the latent space" or figuratively, reading tea leaves.

The story for dictionary-based methods is not better. Researchers usually assume these dictionaries have built-in validity and use them directly in their research. However, multiple validation studies [@boukes2020whatsthetone;@gonzalez2015signals;@ribeiro2016sentibench] demonstrate these dictionaries have very limited criterion validity.

Oolong provides a set of tools to objectively judge substantive interpretability to applied users in disciplines such as political science and communication science. It allows standardized content based testing of topic models as well as dictionary-based methods with clear numeric indicators of semantic validity. Oolong makes it easy to generate standard validation tests suggested by @chang2009reading and @song2020validations.

# Validation of automated content analysis

Validity is a requirement of content analysis [@krippendorff2018content; @neuendorf2016content]. Validation of automated methods has been called for by many scholars, e.g. @grimmer2013text; @ribeiro2016sentibench; @van2018communication. But how to validate these methods? The paper by @dimaggio2013exploiting conceptualizes validation of automated methods as three different operations and the three operations supplement each other. These three operations are: 1) *statistical* validation --to see if the model results agree with the assumptions of the model. Examples of statistical validation are calculation of pointwise mutual information, perplexity or semantic coherence of a topic model; 2) *semantic* validation --to see if the model results are semantically making sense. This procedure involves comparing model results with human judgment [@grimmer2011general]; 3) *predictive* validation --to see if the model results can predict external events [@quinn2010analyze]. For example, one can study whether external events can explain surges in attention to a topic extracted by a topic model.

This package focuses on semantic validation for three reasons: 
First, there is existing architecture for conducting statistical validation and predictive validation. Topic modeling packages such as `text2vec` [@selivanov2020tex2vec], `topicmodels` [@bettina2011topicmodels], and `textmineR` [@jones2019textminer] provide functions to calculate metrics such as perplexity and semantic coherence. Packages such as `stminsights` [@schwemmer2018stminsights] and `LDAvis` [@sievert2015ldavis] offer additional qualitative methods for predictive validation. As of writing, `tosca` [@koppers2020tosca] is the only package dealing with semantic validation. But the text-based interface might pose challenges to human annotators and it can only support topic models from the `lda` package [@change2015lda].

Second, results from statistical validation do not always agree with those from semantic validation. For example, a topic model with a lower perplexity does not have a better interpretability [@chang2009reading]. Of course, there are also metrics from statistical validation that are shown to be correlated with semantic validity, e.g. semantic coherence [@mimno2011optimizing]. But this correlation is also dependent on the text material. For example, @fan2019assessing show that semantic coherence is weakly correlated at best with human assessment, when the text material used for training a topic model has some frequent terms. But still, calculation of semantic coherence is recommended in the best practice paper by @maier2018applying. Nonetheless, conducting only statistical validation is not adequate because these three validation operations supplement each other.

Finally, predictive validation is dependent on research questions and thus it is difficult to be generalized as a reusable software framework. Additionally, the relationship between external (sociopolitical) events and the results from automated content analysis tools is usually what social scientists are eager to study, cf. using topic models for information retrieval [@yi2008evaluating]. We do not believe social scientists would ignore conducting any form of predictive validation.

Oolong focuses on semantic validation. The package provides the "human-in-the-loop" semantic validation procedures suggested by @chang2009reading and @song2020validations. The procedure proposed by @chang2009reading has been adopted in subsequent social science studies as the gold standard to validate topic models, e.g. @bohr2020reporting, @chuang2015topiccheck, and @miller2017australia. The procedure proposed by @song2020validations emphasizes both criterion validity and interrater reliability.

# Semantic validation of topic models

Topic models can be validated by word intrusion test and topic intrusion test [@chang2009reading]. In these tests, a human rater is asked to pick an odd word from a bunch of words (word intrusion test) or pick an odd topic from a bunch of topics for a document (topic intrusion test). Oolong provides an easy-to-use Shiny interface for these tests (Figure 1).

Currently, oolong supports a variety of topic models, e.g. structural topic models / correlated topic models from `stm` [@roberts2019stm], warp-LDA models from `text2vec` [@selivanov2020tex2vec], latent dirichlet allocation / correlated-topic models from `topicmodels` [@bettina2011topicmodels], biterm topic models from `BTM` [@wijffels2020btm] and keyword-assisted topic models from `keyATM` [@eshima2020keyatm].

For instance, `abstracts_stm` is a structural topic model trained with the text data from `abstracts$text` [@chan2020high].

\begin{figure}
\includegraphics[width=0.5\linewidth]{paper_files/fig1} \caption{A screenshot of word intrusion test}\label{fig:unnamed-chunk-1}
\end{figure}


```r
library(stm)
library(tibble)
library(dplyr)
library(quanteda)
library(oolong)
```


```r
abstracts_stm
```

```
## A topic model with 20 topics, 2500 documents and a 3998 word dictionary.
```

The function `create_oolong` creates a test object with both word intrusion test and topic intrusion test.


```r
oolong_test <- create_oolong(input_model = abstracts_stm,
                             input_corpus = abstracts$text)
oolong_test
```

```
## An oolong test object with k = 20, 0 coded.
## Use the method $do_word_intrusion_test() to do word intrusion test.
## With 25 cases of topic intrusion test. 0 coded.
## Use the method $do_topic_intrusion_test() to do topic intrusion test.
## Use the method $lock() to finalize this object and see the results.
```

The tests can be administered with methods `do_word_intrusion_test` and `do_topic_intrusion_test`.

```r
oolong_test$do_word_intrusion_test()
oolong_test$do_topic_intrusion_test()
```

After both tests has been done by a human rater, the test object must be locked and then accuracy metrics such as model precision (MP) and TLO (topic log odd) are displayed. 





```r
oolong_test$lock()
oolong_test
```

```
## An oolong test object with k = 20, 20 coded.
## 95%  precision
## With 25 cases of topic intrusion test. 25 coded.
## TLO: -0.235
```

The suggested workflow is to have at least two human raters to do the same set of tests. Test object can be cloned to allow multiple raters to do the test. More than one test object can be studied together using the function `summarize_oolong()`.


```r
oolong_test_rater1 <- create_oolong(abstracts_stm, abstracts$text)
oolong_test_rater2 <- clone_oolong(oolong_test_rater1)
```

```r
## Let rater 1 do the test.
oolong_test_rater1$do_word_intrusion_test()
oolong_test_rater1$do_topic_intrusion_test()
oolong_test_rater1$lock()

## Let rater 2 do the test.
oolong_test_rater2$do_word_intrusion_test()
oolong_test_rater2$do_topic_intrusion_test()
oolong_test_rater2$lock()
```



Get a summary of the two objects.




```r
summarize_oolong(oolong_test_rater1, oolong_test_rater2)
```

```
## Mean model precision: 0.35
## Quantiles of model precision: 0.3, 0.325, 0.35, 0.375, 0.4
## P-value of the model precision
##  (H0: Model precision is not better than random guess): 0.0089
## Krippendorff's alpha: 0.143
## K Precision:
## 0.5, 0, 0, 0.5, 0, 0.5, 0, 0, 1, 0, 0, 1, 0, 1, 0.5, 0, 0.5, 0.5, 0.5, 0.5
## Mean TLO: -1.58
## Median TLO: -0.08
## Quantiles of TLO: -6.64, -2.9, -0.08, 0, 0
## P-Value of the median TLO 
## (H0: Median TLO is not better than random guess): 0
```

Two key indicators of semantic validity are mean model precision and median TLO. Please interpret the magnitude of the two values [see @chang2009reading] rather than the two statisical tests. The two statistical tests are testing whether the raters did better than random guess. Therefore, rejection of the null hypothesis is just the bare minimum of topic interpretability, *not* an indicator of adquate semantic validity of the topic model. Besides, please a very conservative significant level, e.g. alpha < 0.001.

# Semantic validation of dictionary-based methods

Dictionary-based methods such as AFINN [@nielsen2011new] can be validated by creating a gold standard dataset [@song2020validations]. Oolong provides a workflow for generating such gold standard dataset.

For example, you are interested in studying the sentiment of tweets from Donald Trump. `trump2k` is a random subset of 2,000 tweets from Donald Trump. And you would like to use AFINN to extract sentiment from these tweets. In this analysis, AFINN sentiment score is the *target value*.

A test object can be generated also with `create_oolong`. The argument `construct` should be an adjective, e.g. "positive" or "liberal".


```r
trump <- create_oolong(input_corpus = trump2k,
                       construct = "positive",
                       exact_n = 20)
trump
```

```
## An oolong test object (gold standard generation) with 20 cases, 0 coded.
## Use the method $do_gold_standard_test() to generate gold standard.
## Use the method $lock() to finalize this object and see the results.
```

Similarly, we suggest to have at least two human coders to do the same set of tests.


```r
trump2 <- clone_oolong(trump)
```

Instruct two coders to code the tweets and lock the objects.

```r
trump$do_gold_standard_test()
trump2$do_gold_standard_test()
trump$lock()
trump2$lock()
```



The method `turn_gold` converts a test object into a quanteda corpus [@benoit2018quanteda]. 


```r
gold_standard <- trump$turn_gold()
gold_standard
```

```
## Corpus consisting of 20 documents and 1 docvar.
## text1 :
## "Thank you Eau Claire, Wisconsin.  #VoteTrump on Tuesday, Apr..."
## 
## text2 :
## ""@bobby990r_1: @realDonaldTrump would lead polls the second ..."
## 
## text3 :
## ""@KdanielsK: @misstcassidy @AllAboutTheTea_ @realDonaldTrump..."
## 
## text4 :
## "Thank you for a great afternoon Birmingham, Alabama! #Trump2..."
## 
## text5 :
## ""@THETAINTEDT: @foxandfriends @realDonaldTrump Trump 2016 ht..."
## 
## text6 :
## "People believe CNN these days almost as little as they belie..."
## 
## [ reached max_ndoc ... 14 more documents ]
## Access the answer from the coding with quanteda::docvars(obj, 'answer')
```

This corpus can be used to calculate the target value, e.g. AFINN.


```r
dfm(gold_standard, remove_punct = TRUE) %>% dfm_lookup(afinn) %>%
    quanteda::convert(to = "data.frame") %>%
    mutate(matching_word_valence = (neg5 * -5) + (neg4 * -4) +
               (neg3 * -3) + (neg2 * -2) + (neg1 * -1) +
               (zero * 0) + (pos1 * 1) + (pos2 * 2) + (pos3 * 3) +
               (pos4 * 4) + (pos5 * 5),
           base = ntoken(gold_standard, remove_punct = TRUE),
           afinn_score = matching_word_valence / base) %>% 
		   pull(afinn_score) -> afinn_score
```

Summarize all oolong objects with the target value.


```r
res <- summarize_oolong(trump, trump2, target_value = afinn_score)
```

Printing the summary shows Krippendorff's Alpha, an indicator of interrater reliability. The validity metrics of a text analytic method can be tinted by poor interrater reliability of manual annotations [@song2020validations]. It is important to ensure high interrater reliability first.


```r
res
```

```
## Krippendorff's Alpha: 0.931
## Correlation: 0.744 (p = 2e-04)
## Effect of content length: -0.323 (p = 0.1643)
```

Additional diagnostic plots can also be displayed (Figure 2).


```r
plot(res)
```

![Diagnostic plots generated by oolong](paper_files/figure-latex/diagplot-1.pdf) 

The 4 subplots from left to right, top to bottom are: 

1. Correlation between human judgement and target value - A strong correlation between the two is an indicator of criterion validity of the target value.
2. Bland-Altman plot - If the dots are randomly scattering around the mean value (solid line), it is an indicator of good agreement between human judgement and the target value.
3. Correlation between target value and content length - If there is no strong correlation between the target value and content length, it is an indicator of robustness against the influence of content length [see @chan_4b].
4. Cook's distance of all data points - if there are only a few dots above the threshold (dotted line), it is an indicator of robustness against the influence of outliers.

# Acknowledgements

The development of oolong is partially supported by SAGE Concept Grant.

# References
# gs_creation

    Code
      create_oolong(input_corpus = abstracts$text)
    Message <cliMessage>
      -- oolong (gold standard generation) -------------------------------------------
      i GS: n = 25, 0 coded.
      i Construct:  positive.
      -- Methods --
      * <$do_gold_standard_test()>: generate gold standard
      * <$lock()>: finalize this object and see the results

# gs_turngold

    Code
      x <- create_oolong(input_corpus = abstracts$text)
      x$lock(force = TRUE)
      x$turn_gold()
    Output
      Corpus consisting of 25 documents and 1 docvar.
      text1 :
      "This study examined the relationship among personal network ..."
      
      text2 :
      "Privacy seals were developed to address concerns about onlin..."
      
      text3 :
      "Using a questionnaire validated by the project Biohead-Citiz..."
      
      text4 :
      "Communication and language barriers isolate Deaf American Si..."
      
      text5 :
      "Media use and aging is an important interdisciplinary topic ..."
      
      text6 :
      "In the last decade, research has provided a series of insigh..."
      
      [ reached max_ndoc ... 19 more documents ]
    Message <cliMessage>
      i Access the answer from the coding with quanteda::docvars(obj, 'answer')

# check_calculation_topic_intrusion_multiobject (Printing)

    Code
      res
    Message <cliMessage>
      -- Summary (topic model): ------------------------------------------------------
      -- Word intrusion test --
      i Mean model precision: 1
      i Quantiles of model precision: 1, 1, 1, 1, 1
      i P-value of the model precision
      (H0: Model precision is not better than random guess): 0
      i Krippendorff's alpha: 1
      i K Precision:
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1
      -- Topic intrusion test --
      i Mean TLO: 0
      i Median TLO: 0
      i Quantiles of TLO: 0, 0, 0, 0, 0
      i P-Value of the median TLO 
      (H0: Median TLO is not better than random guess): 0

# ti only

    Code
      create_oolong(input_model = abstracts_keyatm, input_corpus = abstracts$text,
      type = "ti")
    Message <cliMessage>
      -- oolong (topic model) --------------------------------------------------------
      x WI v TI x WSI
      i TI: n = 25, 0 coded.
      -- Methods --
      * <$do_topic_intrusion_test()>: do topic intrusion test
      * <$lock()>: finalize and see the results

# wsi only

    Code
      create_oolong(input_model = abstracts_keyatm, input_corpus = abstracts$text,
      type = "wsi", wsi_n_top_terms = 100)
    Message <cliMessage>
      -- oolong (topic model) --------------------------------------------------------
      x WI x TI v WSI
      i WSI: n = 10, 0 coded.
      -- Methods --
      * <$do_word_set_intrusion_test()>: do word set intrusion test
      * <$lock()>: finalize and see the results

# check_calculation_wsi_multiobject (printing)

    Code
      res
    Message <cliMessage>
      -- Summary (topic model): ------------------------------------------------------
      -- Word set intrusion test --
      i Mean model precision: 0.833333333333333
      i K Precision:
      0.3, 0.7, 0.7, 0.7, 1, 1, 1, 1, 1, 1
      i Krippendorff's alpha: 0.056

---

    Code
      res
    Message <cliMessage>
      -- Summary (topic model): ------------------------------------------------------
      -- Word set intrusion test --
      i Mean model precision: 1
      i K Precision:
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1

# export printing

    Code
      export_oolong(obj1, dir = newdir, verbose = TRUE, use_full_path = FALSE)
    Message <cliMessage>
      i The Shiny has been written to the directory: ~/oolong_testing
      i You can test the app with: shiny::runApp("~/oolong_testing")

---

    Code
      export_oolong(obj1, dir = newdir, verbose = FALSE, use_full_path = FALSE)

# update

    Code
      update_oolong(y)
    Warning <simpleWarning>
      Please consider setting the userid by assigning the userid to the slot $userid, e.g. oolong$userid <- "myname"
      The oolong object is too old. Some security features might not be available in the updated oolong object.
    Message <cliMessage>
      -- oolong (topic model) --------------------------------------------------------
      v WI v TI x WSI
      i WI: k = 10, 0 coded.
      i TI: n = 10, 0 coded.
      -- Methods --
      * <$do_word_intrusion_test()>: do word intrusion test
      * <$do_topic_intrusion_test()>: do topic intrusion test
      * <$lock()>: finalize and see the results

---

    Code
      update_oolong(y)
    Warning <simpleWarning>
      Please consider setting the userid by assigning the userid to the slot $userid, e.g. oolong$userid <- "myname"
      The oolong object is too old. Some security features might not be available in the updated oolong object.
    Message <cliMessage>
      -- oolong (gold standard generation) -------------------------------------------
      i GS: n = 25, 0 coded.
      i Construct:  positive.
      -- Methods --
      * <$do_gold_standard_test()>: generate gold standard
      * <$lock()>: finalize this object and see the results

---

    Code
      update_oolong(y)
    Warning <simpleWarning>
      Please consider setting the userid by assigning the userid to the slot $userid, e.g. oolong$userid <- "myname"
      The oolong object is too old. Some security features might not be available in the updated oolong object.
    Message <cliMessage>
      -- oolong (gold standard generation) -------------------------------------------
      i GS: n = 25, 0 coded.
      i Construct:  positive.
      -- Methods --
      * <$turn_gold()>: convert the test results into a quanteda corpus

---

    

---

    Code
      update_oolong(y)
    Message <cliMessage>
      -- oolong (topic model) --------------------------------------------------------
      v WI x TI x WSI
      i WI: k = 20, 0 coded.
      -- Methods --
      * <$do_word_intrusion_test()>: do word intrusion test
      * <$lock()>: finalize and see the results

---
title: "{title}"
output: rmarkdown::html_vignette
author:
  - Chung-hong Chan ^[University of Mannheim]
vignette: >
  %\VignetteIndexEntry{{title}}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  )
set.seed(46709394)
```

---
title: "{title}"
output: github_document
author:
  - Chung-hong Chan ^[University of Mannheim]
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
  )
set.seed(42)
```
---
title: "BTM"
output: github_document
author:
  - Chung-hong Chan ^[University of Mannheim]
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
  )
set.seed(42)
```

The package BTM by Jan Wijffels et al. finds "topics in collections of short text". Compared to other topic model packages, BTM requires a special data format for training. Oolong has no problem generating word intrusion tests with BTM. However, that special data format can make creation of topic intrusion tests very tricky.

This guide provides our recommendations on how to use BTM, so that the model can be used for generating topic intrusion tests.

# Requirement #1: Keep your quanteda corpus

It is because every document has a unique document id.

```{r}
require(BTM)
require(quanteda)
require(oolong)
trump_corpus <- corpus(trump2k)
```

And then you can do regular text cleaning, stemming procedure with `quanteda`. Instead of making the product a `DFM` object, make it a `token` object. You may read [this issue](https://github.com/quanteda/quanteda/issues/1404) by Benoit et al.

```{r}
tokens(trump_corpus, remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE, split_hyphens = TRUE, remove_url = TRUE) %>% tokens_tolower() %>% tokens_remove(stopwords("en")) %>% tokens_remove("@*")  -> trump_toks
```

# Requirement #2: Keep your data frame

Use this function to convert the `token` object to a data frame.

```{r}
as.data.frame.tokens <- function(x) {
  data.frame(
    doc_id = rep(names(x), lengths(x)),
    tokens = unlist(x, use.names = FALSE)
  )
}

trump_dat <- as.data.frame.tokens(trump_toks)
```

Train a BTM model

```{r, message = FALSE, results = 'hide', warning = FALSE}
trump_btm <- BTM(trump_dat, k = 8, iter = 500, trace = 10)
```

## Pecularities of BTM

This is how you should generate $\theta_{t}$ . However, there are many NaN and there are only 1994 rows (`trump2k` has 2000 tweets) due to empty documents.

```{r}
theta <- predict(trump_btm, newdata = trump_dat)
dim(theta)
```

```{r}
setdiff(docid(trump_corpus), row.names(theta))
```

```{r}
trump_corpus[604]
```

Also, the row order is messed up.

```{r}
head(row.names(theta), 100)
```


# Oolong's support for BTM

Oolong has no problem generating word intrusion test for BTM like you do with other topic models.

```{r}
oolong <- create_oolong(trump_btm)
oolong
```

For generating topic intrusion tests, however, you must provide the data frame you used for training (in this case `trump_dat`). Your `input_corpus` must be a quanteda corpus too.

```{r}
oolong <- create_oolong(trump_btm, trump_corpus, btm_dataframe = trump_dat)
oolong
```

`btm_dataframe` must not be NULL.

```{r, error = TRUE}
oolong <- create_oolong(trump_btm, trump_corpus)
```

`input_corpus` must be a quanteda corpus.

```{r, error = TRUE}
oolong <- create_oolong(trump_btm, trump2k, btm_dataframe = trump_dat)
```

The package BTM by Jan Wijffels et al. finds "topics in collections of short text". Compared to other topic model packages, BTM requires a special data format for training. Oolong has no problem generating word intrusion tests with BTM. However, that special data format can make creation of topic intrusion tests very tricky.

This guide provides our recommendations on how to use BTM, so that the model can be used for generating topic intrusion tests.

# Requirement #1: Keep your quanteda corpus

It is because every document has a unique document id.

```{r}
require(BTM)
require(quanteda)
require(oolong)
trump_corpus <- corpus(trump2k)
```

And then you can do regular text cleaning, stemming procedure with `quanteda`. Instead of making the product a `DFM` object, make it a `token` object. You may read [this issue](https://github.com/quanteda/quanteda/issues/1404) by Benoit et al.

```{r}
tokens(trump_corpus, remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE, split_hyphens = TRUE, remove_url = TRUE) %>% tokens_tolower() %>% tokens_remove(stopwords("en")) %>% tokens_remove("@*")  -> trump_toks
```

# Requirement #2: Keep your data frame

Use this function to convert the `token` object to a data frame.

```{r}
as.data.frame.tokens <- function(x) {
  data.frame(
    doc_id = rep(names(x), lengths(x)),
    tokens = unlist(x, use.names = FALSE)
  )
}

trump_dat <- as.data.frame.tokens(trump_toks)
```

Train a BTM model

```{r, message = FALSE, results = 'hide', warning = FALSE}
trump_btm <- BTM(trump_dat, k = 8, iter = 500, trace = 10)
```

## Pecularities of BTM

This is how you should generate $\theta_{t}$ . However, there are many NaN and there are only 1994 rows (`trump2k` has 2000 tweets) due to empty documents.

```{r}
theta <- predict(trump_btm, newdata = trump_dat)
dim(theta)
```

```{r}
setdiff(docid(trump_corpus), row.names(theta))
```

```{r}
trump_corpus[604]
```

Also, the row order is messed up.

```{r}
head(row.names(theta), 100)
```


# Oolong's support for BTM

Oolong has no problem generating word intrusion test for BTM like you do with other topic models.

```{r}
oolong <- create_oolong(trump_btm)
oolong
```

For generating topic intrusion tests, however, you must provide the data frame you used for training (in this case `trump_dat`). Your `input_corpus` must be a quanteda corpus too.

```{r}
oolong <- create_oolong(trump_btm, trump_corpus, btm_dataframe = trump_dat)
oolong
```

`btm_dataframe` must not be NULL.

```{r, error = TRUE}
oolong <- create_oolong(trump_btm, trump_corpus)
```

`input_corpus` must be a quanteda corpus.

```{r, error = TRUE}
oolong <- create_oolong(trump_btm, trump2k, btm_dataframe = trump_dat)
```
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
set.seed(42)
```

# oolong <img src="man/figures/oolong_logo.svg" align="right" height="200" />

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/chainsawriot/oolong.svg?branch=master)](https://travis-ci.org/chainsawriot/oolong)
[![Codecov test coverage](https://codecov.io/gh/chainsawriot/oolong/branch/master/graph/badge.svg)](https://codecov.io/gh/chainsawriot/oolong?branch=master)
[![joss stataus](https://joss.theoj.org/papers/6e535564e7142d705f4f3d68b18dac62/status.svg)](https://joss.theoj.org/papers/6e535564e7142d705f4f3d68b18dac62)
[![CRAN status](https://www.r-pkg.org/badges/version/oolong)](https://CRAN.R-project.org/package=oolong)
[![R-CMD-check](https://github.com/chainsawriot/oolong/workflows/R-CMD-check/badge.svg)](https://github.com/chainsawriot/oolong/actions)
<!-- badges: end -->

<img src="man/figures/oolong_demo.gif" align="center" height="400" />

The goal of oolong ^[/ˈuːlʊŋ/ 烏龍, literally means "Dark Dragon", is a semi-oxidized tea from Asia. It is very popular in Taiwan, Japan and Hong Kong. In Cantonese and Taiwanese Mandarin, the same word can also mean "confused". It perfectly captures the spirit of human-in-the-loop validation. ] is to generate and administrate validation tests easily for typical automated content analysis tools such as topic models and dictionary-based tools.

Please refer to the [overview](overview_gh.md) for an introduction to this package. If you need to deploy the test online, please refer to the [Deployment Vignette](deploy_gh.md). If you use BTM, please refer to the [BTM Vignette](btm_gh.md).

## Citation

Please cite this package as:

Chan C-h. & Sältzer M., (2020). oolong: An R package for validating  automated content analysis tools. Journal of Open Source Software, 5(55), 2461, [https://doi.org/10.21105/joss.02461](https://doi.org/10.21105/joss.02461)

For a BibTeX entry, use the output from `citation(package = "oolong")`.

## Contributing

Contributions in the form of feedback, comments, code, and bug report are welcome.

* Fork the source code, modify, and issue a [pull request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork).
* Issues, bug reports: [File a Github issue](https://github.com/chainsawriot/oolong).
* Github is not your thing? Contact Chung-hong Chan by e-mail, post, or other methods listed on this [page](https://www.mzes.uni-mannheim.de/d7/en/profiles/chung-hong-chan).

## Code of Conduct

Please note that the oolong project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

---
---
title: "Deploy"
output: github_document
author:
  - Chung-hong Chan ^[University of Mannheim]
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/vignettes/figures/README-",
  out.width = "100%"
  )
set.seed(42)
```
In oolong 0.3.22, functions for deploying oolong tests were added (`export_oolong`, `revert_oolong` etc.). These functions make it possible for the coders to conduct validation tests online using their browser, rather than having to install R on their computer.

The basic workflow is simple: 1) create the oolong test object as usual; 2) deploy the test online and obtain the URL to the test; 3) ask your coders to conduct the test online and send back the data file; 4) revert back from the data file to an oolong object.

# Create an oolong test

Please note that one cannot deploy oolong test objects with *both* word and topic intrusion tests, i.e. those created using `witi()` online. If you need to do both tests, you need to deploy them as two separate instances: one created using `wi()` and another created using `ti()`.

In this guide, we assume you want to deploy a word set intrusion test online.

```{r}
library(oolong)
wsi_test <- wsi(abstracts_keyatm)
wsi_test
```

# Deploy the test online

First, you need to export the oolong test object as a stand alone Shiny app. This stand alone Shiny app will be in a directory.

```{r}
export_oolong(wsi_test, dir = "./wsi_test", use_full_path = FALSE)
```

The directory has only two files

```{r}
fs::dir_tree("./wsi_test")
```

This structure is called ["Single-file Shiny app."](https://shiny.rstudio.com/articles/app-formats.html) Experienced Shiny users might have their preferred method of deploying this app to whatever Shiny server they can master.

For less experienced users, the simplest way to deploy this app online is to use [shinyapps.io](https://www.shinyapps.io/) (free tier available with 25 hours of computational time per month). Please register for an account at shinyapps.io and configure rsconnect. Please refer to [this guide](https://shiny.rstudio.com/articles/shinyapps.html) for more information. Please remember to configure the tokens.

```r
## replace <ACCOUNT>, <TOKEN>, <SECRET> with the information from your profile on Shinyapps.io: click Your name -> Tokens
rsconnect::setAccountInfo(name="<ACCOUNT>", token="<TOKEN>", secret="<SECRET>")
```

For RStudio users, the simplest way to deploy the app to shinyapps.io is to first launch the app.

```{r, eval = FALSE}
library(shiny)
runApp("./wsi_test")
```

And then click the **Publish** button at the right corner of the launched window.

You will be asked for the title of the app, just give it a name, e.g. *wsi_test*. You probably can keep other default settings and push the **Publish** button to initialize the deployment process.

<img src="vignettes/figures/deploying_shinyappsio.png" align="center" height="400" />

If there is no hiccup, you will get a URL to your deployed oolong test. Something like: *https://yourname.shinyapps.io/wsi_test/*

# Conduct the test

You can give the URL to your coders and they conduct the test with their browser online. The only difference of the deployed version is that, there will be a userid prompt and download button after the coding.

<img src="vignettes/figures/oolong_download.png" align="center"/>

You should instruct your coders to download the data file after coding and return it to you. ^[Future versions might provide permanent storage]

# Revert

You can then obtain a locked oolong object from the original oolong and the downloaded data file. `revert_oolong` will do verifications with the original oolong object to make sure no error and no cheating.

```{r, include = FALSE}
wsi_test <- readRDS(system.file("extdata", "wsi_test.RDS", package = "oolong"))
```

```r
revert_oolong(wsi_test, "oolong_2021-05-22 20 51 26 Hadley Wickham.RDS")
```

```{r, echo = FALSE}
revert_oolong(wsi_test, system.file("extdata", "hadley.RDS", package = "oolong"))
```

```{r, include = FALSE}
unlink("./wsi_test", recursive = TRUE)
```
---
title: "Overview"
output: github_document
author:
  - Chung-hong Chan ^[University of Mannheim]
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
  )
set.seed(42)
```
The validation test is called "oolong test" (for reading tea leaves). This package provides several functions for generating different types of oolong test.

| function | purpose                                                                                                                           |
|---------:|:----------------------------------------------------------------------------------------------------------------------------------|
|   `wi()` | validating a topic model with [word intrusion test](#word-intrusion-test) (Chang et al., 2008)                                    |
|   `ti()` | validating a topic model with [topic intrusion test](#topic-intrusion-test) (Chang et al., 2008; aka "T8WSI" in Ying et al. 2021) |
| `witi()` | validating a topic model with [word intrusion test](#word-intrusion-test) and [topic intrusion test](#topic-intrusion-test)       |
|  `wsi()` | validating a topic model with [word set intrusion test](#word-set-intrusion-test) (Ying et al. 2021)                              |
|   `gs()` | oolong test for [creating gold standard](#creating-gold-standard) (see Song et al., 2020)                                         |

All of these tests can also be generated with the function [`create_oolong`](#backward-compatibility). As of version 0.3.20, it is no longer recommended.

## Installation

Because the package is constantly changing, we suggest using the development version from [GitHub](https://github.com/):

``` r
# install.packages("devtools")
devtools::install_github("chainsawriot/oolong")
```

You can also install the "stable" (but slightly older) version from CRAN:

```r
install.packages("oolong")
```

## Validating Topic Models

#### Word intrusion test

`abstracts_keyatm` is an example topic model trained with the data `abstracts` using the `keyATM` package. Currently, this package supports structural topic models / correlated topic models from `stm`, Warp LDA models from `text2vec` , LDA/CTM models from `topicmodels`,  Biterm Topic Models from `BTM`, Keyword Assisted Topic Models from `keyATM`, and seeded LDA models from `seededlda`. Although not strictly a topic model, Naive Bayes models from `quanteda.textmodels` are also supported. See the section on [Naive Bayes](#about-naive-bayes) for more information.

```{r}
library(oolong)
library(keyATM)
library(quanteda)
library(dplyr)
```

```{r example}
abstracts_keyatm
```

To create an oolong test with word intrusion test, use the function `wi`. It is recommended to provide a user id of coder who are going to be doing the test.

```{r createtest}
oolong_test <- wi(abstracts_keyatm, userid = "Hadley")
oolong_test
```

As instructed, use the method `$do_word_intrusion_test()` to start coding. 

```{r, eval = FALSE}
oolong_test$do_word_intrusion_test()
```

After the coding, you need to first lock the test. Then, you can look at the model precision by printing the oolong test.

```{r, include = FALSE}
### Mock this process
oolong_test$.__enclos_env__$private$test_content$wi$answer <- oolong_test$.__enclos_env__$private$test_content$wi$intruder
oolong_test$.__enclos_env__$private$test_content$wi$answer[1] <- "wronganswer"
```

```{r lock}
oolong_test$lock()
oolong_test
```

#### Word set intrusion test

Word set intrusion test is a variant of word intrusion test (Ying et al., 2021), in which multiple word sets generated from top terms of one topic are juxtaposed with one intruder word set generated similarly from another topic. In Ying et al., this test is called "R4WSI" because 4 word sets are displayed. By default, oolong generates also R4WSI. However, it is also possible to generate R(N)WSI by setting the parameter `n_correct_ws` to N - 1.

```{r wsi1}
oolong_test <- wsi(abstracts_keyatm, userid = "Garrett")
oolong_test
```

Use the method `$do_word_set_intrusion_test()` to start coding. 

```{r wsi2, eval = FALSE}
oolong_test$do_word_set_intrusion_test()
```

```{r, include = FALSE}
### Mock this process
oolong_test$.__enclos_env__$private$test_content$wsi$answer <- oolong_test$.__enclos_env__$private$test_content$wsi$intruder
oolong_test$.__enclos_env__$private$test_content$wsi$answer[1] <- "wronganswer"
```

```{r wsi3}
oolong_test$lock()
oolong_test
```

#### Topic intrusion test

For example, `abstracts_keyatm` was generated with the corpus `abstracts$text`

```{r newgroup5}
library(tibble)
abstracts
```

Creating the oolong test object with the corpus used for training the topic model will generate topic intrusion test cases. 

```{r createtest2}
oolong_test <- ti(abstracts_keyatm, abstracts$text, userid = "Julia")
oolong_test
```

Similarly, use the `$do_topic_intrusion_test` to code the test cases, lock the test with `$lock()` and then you can look at the TLO (topic log odds) value by printing the oolong test.

```{r, eval = FALSE}
oolong_test$do_topic_intrusion_test()
oolong_test$lock()
```

```{r, include = FALSE}
genius_topic <- function(obj1) {
    obj1$.__enclos_env__$private$test_content$ti$answer <- obj1$.__enclos_env__$private$test_content$ti$intruder
    return(obj1)
}
genius_word <- function(obj1) {
    obj1$.__enclos_env__$private$test_content$wi$answer <- obj1$.__enclos_env__$private$test_content$wi$intruder
    return(obj1)
}
oolong_test <- genius_word(genius_topic(oolong_test))
oolong_test$.__enclos_env__$private$test_content$ti$answer[2] <- sample(oolong_test$.__enclos_env__$private$test_content$ti$candidates[[2]], 1)
oolong_test$lock()
```

```{r topic_res}
oolong_test
```

### Suggested workflow

The test makes more sense if more than one coder is involved. A suggested workflow is to create the test, then clone the oolong object. Ask multiple coders to do the test(s) and then summarize the results.

Preprocess and create a document-feature matrix

```{r, eval = FALSE}
dfm(abstracts$text, tolower = TRUE, stem = TRUE, remove = stopwords('english'), remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE, remove_hyphens = TRUE) %>% dfm_trim(min_docfreq = 3, max_docfreq = 500) %>% dfm_select(min_nchar = 3, pattern = "^[a-zA-Z]+$", valuetype = "regex") -> abstracts_dfm
```

Train a topic model.

```{r step0, eval = FALSE}
require(keyATM)
abstracts_keyatm <- keyATM(keyATM_read(abstracts_dfm), no_keyword_topics = 0, keywords = abstracts_dictionary, model = "base", options = list(seed = 46709394))
```

Create a new oolong object.

```{r step1}
oolong_test_rater1 <- witi(abstracts_keyatm, abstracts$text, userid = "Yihui")
```

Clone the oolong object to be used by other raters.

```{r step2}
oolong_test_rater2 <- clone_oolong(oolong_test_rater1, userid = "Jenny")
```

Ask different coders to code each object and then lock the object.

```{r, eval = FALSE}
## Let Yihui do the test.
oolong_test_rater1$do_word_intrusion_test()
oolong_test_rater1$do_topic_intrusion_test()
oolong_test_rater1$lock()

## Let Jenny do the test.
oolong_test_rater2$do_word_intrusion_test()
oolong_test_rater2$do_topic_intrusion_test()
oolong_test_rater2$lock()
```

```{r, include = FALSE}
### Mock this process
set.seed(46709394)
oolong_test_rater1 <- oolong:::.monkey_test(oolong_test_rater1, intelligent = 0.3)
oolong_test_rater2 <- oolong:::.monkey_test(oolong_test_rater2, intelligent = 0)
oolong_test_rater1$lock()
oolong_test_rater2$lock()
```

Get a summary of the two objects.

```{r, step3}
summarize_oolong(oolong_test_rater1, oolong_test_rater2)
```

### About the p-values

The test for model precision (MP) is based on an one-tailed, one-sample binomial test for each rater. In a multiple-rater situation, the p-values from all raters are combined using the Fisher's method (a.k.a. Fisher's omnibus test).

H0: MP is not better than 1/ n\_top\_terms

H1: MP is better than 1/ n\_top\_terms


The test for the median of TLO is based on a permutation test.

H0: Median TLO is not better than random guess.

H1: Median TLO is better than random guess.

One must notice that the two statistical tests are testing the bear minimum. A significant test only indicates the topic model can make the rater(s) perform better than random guess. It is not an indication of good topic interpretability. Also, one should use a very conservative significant level, e.g. $\alpha < 0.001$.

### About Warp LDA

There is a subtle difference between the support for `stm` and for `text2vec`.

`abstracts_warplda` is a Warp LDA object trained with the same dataset as the `abstracts_stm`

```{r warplda}
abstracts_warplda
```

All the API endpoints are the same, except the one for the creation of topic intrusion test cases. You must supply also the `input_dfm`.

```{r warplda2}
### Just word intrusion test.
oolong_test <- wi(abstracts_warplda, userid = "Lionel")
oolong_test
```


```{r warplda3}
abstracts_dfm
```

```{r warplda4, , message = FALSE, results = 'hide', warning = FALSE}
oolong_test <- witi(abstracts_warplda, abstracts$text, input_dfm = abstracts_dfm, userid = "Mara")
```

```{r warplda5}
oolong_test
```

## About Biterm Topic Model

Please refer to the vignette about BTM.

## About Naive Bayes

Naive Bayes model is a supervised machine learning model. This package supports Naive Bayes models trained using `quanteda.textmodels`.

Suppose `newsgroup_nb` is a Naive Bayes model trained on a subset of the classic [20 newsgroups] dataset.

```r
tokens(newsgroup5$text, remove_punct = TRUE, remove_symbols = TRUE, remove_numbers = TRUE, remove_url = TRUE, spilit_hyphens = TRUE) %>% tokens_wordstem %>% tokens_remove(stopwords("en")) %>% dfm(tolower = TRUE) %>% dfm_trim(min_termfreq = 3, max_docfreq = 0.06, docfreq_type = "prop") -> newsgroup_dfm
docvars(newsgroup_dfm, "group") <- newsgroup5$title
newsgroup_nb <- textmodel_nb(newsgroup_dfm, docvars(newsgroup_dfm, "group"), distribution = "Bernoulli")
```

You can still generate word intrusion and word set intrusion tests.

```{r}
wi(newsgroup_nb)
```

```{r}
wsi(newsgroup_nb)
```

## Validating Dictionary-based Methods

### Creating gold standard

`trump2k` is a dataset of 2,000 tweets from \@realdonaldtrump.

```{r trump2k}
tibble(text = trump2k)
```

For example, you are interested in studying the sentiment of these tweets. One can use tools such as AFINN to automatically extract sentiment in these tweets. However, oolong recommends to generate gold standard by human coding first using a subset. By default, oolong selects 1% of the origin corpus as test cases. The parameter `construct` should be an adjective, e.g. positive, liberal, populistic, etc.

```{r goldstandard}
oolong_test <- gs(input_corpus = trump2k, construct = "positive", userid = "Joe")
oolong_test
```

As instructed, use the method `$do_gold_standard_test()` to start coding.

```{r, eval = FALSE}
oolong_test$do_gold_standard_test()
```

After the coding, you need to first lock the test and then the `$turn_gold()` method is available.

```{r, include = FALSE}
oolong_test$.__enclos_env__$private$test_content$gs <- 
structure(list(case = 1:20, text = c("Thank you Eau Claire, Wisconsin. \n#VoteTrump on Tuesday, April 5th!\nMAKE AMERICA GREAT AGAIN! https://t.co/JI5JqwHnMC", 
"\"@bobby990r_1: @realDonaldTrump would lead polls the second he announces candidacy! America is waiting for him to LEAD us out of this mess!", 
"\"@KdanielsK: @misstcassidy @AllAboutTheTea_ @realDonaldTrump My money is on Kenya getting fired first.\"", 
"Thank you for a great afternoon Birmingham, Alabama! #Trump2016 #MakeAmericaGreatAgain https://t.co/FrOkqCzBoD", 
"\"@THETAINTEDT: @foxandfriends @realDonaldTrump Trump 2016 http://t.co/UlQWGKUrCJ\"", 
"People believe CNN these days almost as little as they believe Hillary....that's really saying something!", 
"It was great being in Michigan. Remember, I am the only presidential candidate who will bring jobs back to the U.S.and protect car industry!", 
"\"@DomineekSmith: @realDonaldTrump is the best Republican presidential candidate of all time.\"  Thank you.", 
"Word is that little Morty Zuckerman’s @NYDailyNews loses more than $50 million per year---can that be possible?", 
"\"@Chevy_Mama: @realDonaldTrump I'm obsessed with @celebrityapprenticeNBC. Honestly,  Mr Trump, you are very inspiring\"", 
"President Obama said \"ISIL continues to shrink\" in an interview just hours before the horrible attack in Paris. He is just so bad! CHANGE.", 
".@HillaryClinton loves to lie. America has had enough of the CLINTON'S! It is time to #DrainTheSwamp! Debates https://t.co/3Mz4T7qTTR", 
"\"@jerrimoore: @realDonaldTrump awesome. A treat to get to see the brilliant Joan Rivers once more #icon\"", 
"\"@shoegoddesss: @realDonaldTrump Will definitely vote for you. Breath of fresh air. America needs you!\"", 
"Ted is the ultimate hypocrite. Says one thing for money, does another for votes. \nhttps://t.co/hxdfy0mjVw", 
"\"@Lisa_Milicaj: Truth be told, I  never heard of The National Review until they \"tried\" to declare war on you. No worries, you got my vote!\"", 
"THANK YOU Daytona Beach, Florida!\n#MakeAmericaGreatAgain https://t.co/IAcLfXe463", 
"People rarely say that many conservatives didn’t vote for Mitt Romney. If I can get them to vote for me, we win in a landslide.", 
"Trump National Golf Club, Washington, D.C. is on 600 beautiful acres fronting the Potomac River. A fantastic setting! http://t.co/pYtkbyKwt5", 
"\"@DRUDGE_REPORT: REUTERS 5-DAY ROLLING POLL: TRUMP 34%, CARSON 19.6%, RUBIO 9.7%, CRUZ 7.7%...\" Thank you - a great honor!"
), answer = c(4L, 4L, 2L, 5L, 3L, 2L, 4L, 5L, 2L, 4L, 1L, 1L, 
4L, 4L, 2L, 4L, 4L, 4L, 4L, 4L), target_value = c(NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
NA)), row.names = c(NA, -20L), class = c("tbl_df", "tbl", "data.frame"
))
```

```{r gs_locking}
oolong_test$lock()
oolong_test
```

### Example: Validating AFINN using the gold standard

A locked oolong test can be converted into a quanteda-compatible corpus for further analysis. The corpus contains two `docvars`, 'answer'.

```{r}
oolong_test$turn_gold()
```

In this example, we calculate the AFINN score for each tweet using quanteda. The dictionary `afinn` is bundle with this package.

```{r}
gold_standard <- oolong_test$turn_gold()
dfm(gold_standard, remove_punct = TRUE) %>% dfm_lookup(afinn) %>% quanteda::convert(to = "data.frame") %>%
    mutate(matching_word_valence = (neg5 * -5) + (neg4 * -4) + (neg3 * -3) + (neg2 * -2) + (neg1 * -1)
           + (zero * 0) + (pos1 * 1) + (pos2 * 2) + (pos3 * 3) + (pos4 * 4) + (pos5 * 5),
           base = ntoken(gold_standard, remove_punct = TRUE), afinn_score = matching_word_valence / base) %>%
    pull(afinn_score) -> all_afinn_score
all_afinn_score
```

Put back the vector of AFINN score into the respective `docvars` and study the correlation between the gold standard and AFINN.

```{r}
summarize_oolong(oolong_test, target_value = all_afinn_score)
```

### Suggested workflow

Create an oolong object, clone it for another coder. According to Song et al. (Forthcoming), you should at least draw 1% of your data.

```{r}
trump <- gs(input_corpus = trump2k, exact_n = 40, userid = "JJ")
trump2 <- clone_oolong(trump, userid = "Winston")
```

Instruct two coders to code the tweets and lock the objects.

```{r, eval = FALSE}
trump$do_gold_standard_test()
trump2$do_gold_standard_test()
trump$lock()
trump2$lock()
```

```{r, include = FALSE}
trump$.__enclos_env__$private$test_content$gs <- 
structure(list(case = 1:20, text = c("Thank you Eau Claire, Wisconsin. \n#VoteTrump on Tuesday, April 5th!\nMAKE AMERICA GREAT AGAIN! https://t.co/JI5JqwHnMC", 
"\"@bobby990r_1: @realDonaldTrump would lead polls the second he announces candidacy! America is waiting for him to LEAD us out of this mess!", 
"\"@KdanielsK: @misstcassidy @AllAboutTheTea_ @realDonaldTrump My money is on Kenya getting fired first.\"", 
"Thank you for a great afternoon Birmingham, Alabama! #Trump2016 #MakeAmericaGreatAgain https://t.co/FrOkqCzBoD", 
"\"@THETAINTEDT: @foxandfriends @realDonaldTrump Trump 2016 http://t.co/UlQWGKUrCJ\"", 
"People believe CNN these days almost as little as they believe Hillary....that's really saying something!", 
"It was great being in Michigan. Remember, I am the only presidential candidate who will bring jobs back to the U.S.and protect car industry!", 
"\"@DomineekSmith: @realDonaldTrump is the best Republican presidential candidate of all time.\"  Thank you.", 
"Word is that little Morty Zuckerman’s @NYDailyNews loses more than $50 million per year---can that be possible?", 
"\"@Chevy_Mama: @realDonaldTrump I'm obsessed with @celebrityapprenticeNBC. Honestly,  Mr Trump, you are very inspiring\"", 
"President Obama said \"ISIL continues to shrink\" in an interview just hours before the horrible attack in Paris. He is just so bad! CHANGE.", 
".@HillaryClinton loves to lie. America has had enough of the CLINTON'S! It is time to #DrainTheSwamp! Debates https://t.co/3Mz4T7qTTR", 
"\"@jerrimoore: @realDonaldTrump awesome. A treat to get to see the brilliant Joan Rivers once more #icon\"", 
"\"@shoegoddesss: @realDonaldTrump Will definitely vote for you. Breath of fresh air. America needs you!\"", 
"Ted is the ultimate hypocrite. Says one thing for money, does another for votes. \nhttps://t.co/hxdfy0mjVw", 
"\"@Lisa_Milicaj: Truth be told, I  never heard of The National Review until they \"tried\" to declare war on you. No worries, you got my vote!\"", 
"THANK YOU Daytona Beach, Florida!\n#MakeAmericaGreatAgain https://t.co/IAcLfXe463", 
"People rarely say that many conservatives didn’t vote for Mitt Romney. If I can get them to vote for me, we win in a landslide.", 
"Trump National Golf Club, Washington, D.C. is on 600 beautiful acres fronting the Potomac River. A fantastic setting! http://t.co/pYtkbyKwt5", 
"\"@DRUDGE_REPORT: REUTERS 5-DAY ROLLING POLL: TRUMP 34%, CARSON 19.6%, RUBIO 9.7%, CRUZ 7.7%...\" Thank you - a great honor!"
), answer = c(4L, 4L, 2L, 5L, 3L, 2L, 4L, 5L, 2L, 4L, 1L, 1L, 
4L, 4L, 2L, 4L, 4L, 4L, 4L, 4L), target_value = c(NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
NA)), row.names = c(NA, -20L), class = c("tbl_df", "tbl", "data.frame"
                                         ))

trump2$.__enclos_env__$private$test_content$gs <- 
structure(list(case = 1:20, text = c("Thank you Eau Claire, Wisconsin. \n#VoteTrump on Tuesday, April 5th!\nMAKE AMERICA GREAT AGAIN! https://t.co/JI5JqwHnMC", 
"\"@bobby990r_1: @realDonaldTrump would lead polls the second he announces candidacy! America is waiting for him to LEAD us out of this mess!", 
"\"@KdanielsK: @misstcassidy @AllAboutTheTea_ @realDonaldTrump My money is on Kenya getting fired first.\"", 
"Thank you for a great afternoon Birmingham, Alabama! #Trump2016 #MakeAmericaGreatAgain https://t.co/FrOkqCzBoD", 
"\"@THETAINTEDT: @foxandfriends @realDonaldTrump Trump 2016 http://t.co/UlQWGKUrCJ\"", 
"People believe CNN these days almost as little as they believe Hillary....that's really saying something!", 
"It was great being in Michigan. Remember, I am the only presidential candidate who will bring jobs back to the U.S.and protect car industry!", 
"\"@DomineekSmith: @realDonaldTrump is the best Republican presidential candidate of all time.\"  Thank you.", 
"Word is that little Morty Zuckerman’s @NYDailyNews loses more than $50 million per year---can that be possible?", 
"\"@Chevy_Mama: @realDonaldTrump I'm obsessed with @celebrityapprenticeNBC. Honestly,  Mr Trump, you are very inspiring\"", 
"President Obama said \"ISIL continues to shrink\" in an interview just hours before the horrible attack in Paris. He is just so bad! CHANGE.", 
".@HillaryClinton loves to lie. America has had enough of the CLINTON'S! It is time to #DrainTheSwamp! Debates https://t.co/3Mz4T7qTTR", 
"\"@jerrimoore: @realDonaldTrump awesome. A treat to get to see the brilliant Joan Rivers once more #icon\"", 
"\"@shoegoddesss: @realDonaldTrump Will definitely vote for you. Breath of fresh air. America needs you!\"", 
"Ted is the ultimate hypocrite. Says one thing for money, does another for votes. \nhttps://t.co/hxdfy0mjVw", 
"\"@Lisa_Milicaj: Truth be told, I  never heard of The National Review until they \"tried\" to declare war on you. No worries, you got my vote!\"", 
"THANK YOU Daytona Beach, Florida!\n#MakeAmericaGreatAgain https://t.co/IAcLfXe463", 
"People rarely say that many conservatives didn’t vote for Mitt Romney. If I can get them to vote for me, we win in a landslide.", 
"Trump National Golf Club, Washington, D.C. is on 600 beautiful acres fronting the Potomac River. A fantastic setting! http://t.co/pYtkbyKwt5", 
"\"@DRUDGE_REPORT: REUTERS 5-DAY ROLLING POLL: TRUMP 34%, CARSON 19.6%, RUBIO 9.7%, CRUZ 7.7%...\" Thank you - a great honor!"
), answer = c(5L, 3L, 2L, 5L, 3L, 1L, 4L, 5L, 2L, 4L, 1L, 1L, 
4L, 4L, 2L, 4L, 4L, 4L, 4L, 4L), target_value = c(NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
NA)), row.names = c(NA, -20L), class = c("tbl_df", "tbl", "data.frame"
                                         ))
trump$lock()
trump2$lock()
```

Calculate the target value (in this case, the AFINN score) by turning one object into a corpus.

```{r}
gold_standard <- trump$turn_gold()
dfm(gold_standard, remove_punct = TRUE) %>% dfm_lookup(afinn) %>% quanteda::convert(to = "data.frame") %>%
    mutate(matching_word_valence = (neg5 * -5) + (neg4 * -4) + (neg3 * -3) + (neg2 * -2) + (neg1 * -1)
           + (zero * 0) + (pos1 * 1) + (pos2 * 2) + (pos3 * 3) + (pos4 * 4) + (pos5 * 5),
           base = ntoken(gold_standard, remove_punct = TRUE), afinn_score = matching_word_valence / base) %>%
    pull(afinn_score) -> target_value
```

Summarize all oolong objects with the target value.

```{r}
res <- summarize_oolong(trump, trump2, target_value = target_value)
```

Read the results. The diagnostic plot consists of 4 subplots. It is a good idea to read Bland & Altman (1986) on the difference between correlation and agreement.

* Subplot (top left): Raw correlation between human judgement and target value. One should want to have a good correlation between the two.
* Subplot (top right): Bland-Altman plot. One should want to have no correlation. Also, the dots should be randomly scattering around the mean value. If it is so, the two measurements (human judgement and target value) are in good agreement.
* Subplot (bottom left): Raw correlation between target value and content length. One should want to have no correlation, as an indication of good reliability against the influence of content length. (See Chan et al. 2020)
* Subplot (bottom right): Cook's distance of all data point. One should want to have no dot (or at least very few dots) above the threshold. It is an indication of how the raw correlation between human judgement and target value can or cannot be influenced by extreme values in your data. 

The textual output contains the Krippendorff's alpha of the codings by your raters. In order to claim validity of your target value, you must first establish the reliability of your gold standard. Song et al. [Forthcoming] suggest Krippendorff's Alpha > 0.7 as an acceptable cut-off.

```{r}
res
```

```{r diagnosis}
plot(res)
```

## Backward compatibility

Historically, oolong test objects could only be generated with only one function: `create_oolong`. It is no longer the case and no longer recommended anymore. It is still retained for backward compatibility purposes. If you still need to use  `create_oolong()`, the most important parameters are `input_model` and `input_corpus`. Setting each of them to `NULL` generates different tests.

| input\_model | input\_corpus | output                                                                                                                                      |
|--------------|:-------------:|---------------------------------------------------------------------------------------------------------------------------------------------|
| Not NULL     | NULL          | oolong test for validating a topic model with [word intrusion test](#word-intrusion-test)                                                   |
| Not NULL     | Not NULL      | oolong test for validating a topic model with [word intrusion test](#word-intrusion-test) and [topic intrusion test](#topic-intrusion-test) |
| NULL         | Not NULL      | oolong test for [creating gold standard](#creating-gold-standard)                                                                           |
| NULL         | NULL          | error                                                                                                                                       |


## References

1. Chang, J., Gerrish, S., Wang, C., Boyd-Graber, J. L., & Blei, D. M. (2009). Reading tea leaves: How humans interpret topic models. In Advances in neural information processing systems (pp. 288-296). [link](https://papers.nips.cc/paper/3700-reading-tea-leaves-how-humans-interpret-topic-models)
2. Ying, L., Montgomery, J. M., & Stewart, B. M. (2021). Inferring concepts from topics: Towards procedures for validating topics as measures. Political Analysis. [link](https://doi.org/10.1017/pan.2021.33)
3. Song et al. (2020) In validations we trust? The impact of imperfect human annotations as a gold standard on the quality of validation of automated content analysis. Political Communication. [link](https://doi.org/10.1080/10584609.2020.1723752)
4. Bland, J. M., & Altman, D. (1986). Statistical methods for assessing agreement between two methods of clinical measurement. The lancet, 327(8476), 307-310.
5. Chan et al. (2020) Four best practices for measuring news sentiment using ‘off-the-shelf’ dictionaries: a large-scale p-hacking experiment. Computational Communication Research. [link](https://osf.io/preprints/socarxiv/np5wa/)
6. Nielsen, F. Å. (2011). A new ANEW: Evaluation of a word list for sentiment analysis in microblogs. arXiv preprint arXiv:1103.2903. [link](https://arxiv.org/abs/1103.2903)

---
In oolong 0.3.22, functions for deploying oolong tests were added (`export_oolong`, `revert_oolong` etc.). These functions make it possible for the coders to conduct validation tests online using their browser, rather than having to install R on their computer.

The basic workflow is simple: 1) create the oolong test object as usual; 2) deploy the test online and obtain the URL to the test; 3) ask your coders to conduct the test online and send back the data file; 4) revert back from the data file to an oolong object.

# Create an oolong test

Please note that one cannot deploy oolong test objects with *both* word and topic intrusion tests, i.e. those created using `witi()` online. If you need to do both tests, you need to deploy them as two separate instances: one created using `wi()` and another created using `ti()`.

In this guide, we assume you want to deploy a word set intrusion test online.

```{r}
library(oolong)
wsi_test <- wsi(abstracts_keyatm)
wsi_test
```

# Deploy the test online

First, you need to export the oolong test object as a stand alone Shiny app. This stand alone Shiny app will be in a directory.

```{r}
export_oolong(wsi_test, dir = "./wsi_test", use_full_path = FALSE)
```

The directory has only two files

```{r}
fs::dir_tree("./wsi_test")
```

This structure is called ["Single-file Shiny app."](https://shiny.rstudio.com/articles/app-formats.html) Experienced Shiny users might have their preferred method of deploying this app to whatever Shiny server they can master.

For less experienced users, the simplest way to deploy this app online is to use [shinyapps.io](https://www.shinyapps.io/) (free tier available with 25 hours of computational time per month). Please register for an account at shinyapps.io and configure rsconnect. Please refer to [this guide](https://shiny.rstudio.com/articles/shinyapps.html) for more information. Please remember to configure the tokens.

```r
## replace <ACCOUNT>, <TOKEN>, <SECRET> with the information from your profile on Shinyapps.io: click Your name -> Tokens
rsconnect::setAccountInfo(name="<ACCOUNT>", token="<TOKEN>", secret="<SECRET>")
```

For RStudio users, the simplest way to deploy the app to shinyapps.io is to first launch the app.

```{r, eval = FALSE}
library(shiny)
runApp("./wsi_test")
```

And then click the **Publish** button at the right corner of the launched window.

You will be asked for the title of the app, just give it a name, e.g. *wsi_test*. You probably can keep other default settings and push the **Publish** button to initialize the deployment process.

<img src="figures/deploying_shinyappsio.png" align="center" height="400" />

If there is no hiccup, you will get a URL to your deployed oolong test. Something like: *https://yourname.shinyapps.io/wsi_test/*

# Conduct the test

You can give the URL to your coders and they conduct the test with their browser online. The only difference of the deployed version is that, there will be a userid prompt and download button after the coding.

<img src="figures/oolong_download.png" align="center"/>

You should instruct your coders to download the data file after coding and return it to you. ^[Future versions might provide permanent storage]

# Revert

You can then obtain a locked oolong object from the original oolong and the downloaded data file. `revert_oolong` will do verifications with the original oolong object to make sure no error and no cheating.

```{r, include = FALSE}
wsi_test <- readRDS(system.file("extdata", "wsi_test.RDS", package = "oolong"))
```

```r
revert_oolong(wsi_test, "oolong_2021-05-22 20 51 26 Hadley Wickham.RDS")
```

```{r, echo = FALSE}
revert_oolong(wsi_test, system.file("extdata", "hadley.RDS", package = "oolong"))
```

```{r, include = FALSE}
unlink("./wsi_test", recursive = TRUE)
```
The validation test is called "oolong test" (for reading tea leaves). This package provides several functions for generating different types of oolong test.

| function | purpose                                                                                                                           |
|---------:|:----------------------------------------------------------------------------------------------------------------------------------|
|   `wi()` | validating a topic model with [word intrusion test](#word-intrusion-test) (Chang et al., 2008)                                    |
|   `ti()` | validating a topic model with [topic intrusion test](#topic-intrusion-test) (Chang et al., 2008; aka "T8WSI" in Ying et al. 2021) |
| `witi()` | validating a topic model with [word intrusion test](#word-intrusion-test) and [topic intrusion test](#topic-intrusion-test)       |
|  `wsi()` | validating a topic model with [word set intrusion test](#word-set-intrusion-test) (Ying et al. 2021)                              |
|   `gs()` | oolong test for [creating gold standard](#creating-gold-standard) (see Song et al., 2020)                                         |

All of these tests can also be generated with the function [`create_oolong`](#backward-compatibility). As of version 0.3.20, it is no longer recommended.

## Installation

Because the package is constantly changing, we suggest using the development version from [GitHub](https://github.com/):

``` r
# install.packages("devtools")
devtools::install_github("chainsawriot/oolong")
```

You can also install the "stable" (but slightly older) version from CRAN:

```r
install.packages("oolong")
```

## Validating Topic Models

#### Word intrusion test

`abstracts_keyatm` is an example topic model trained with the data `abstracts` using the `keyATM` package. Currently, this package supports structural topic models / correlated topic models from `stm`, Warp LDA models from `text2vec` , LDA/CTM models from `topicmodels`,  Biterm Topic Models from `BTM`, Keyword Assisted Topic Models from `keyATM`, and seeded LDA models from `seededlda`. Although not strictly a topic model, Naive Bayes models from `quanteda.textmodels` are also supported. See the section on [Naive Bayes](#about-naive-bayes) for more information.

```{r}
library(oolong)
library(keyATM)
library(quanteda)
library(dplyr)
```

```{r example}
abstracts_keyatm
```

To create an oolong test with word intrusion test, use the function `wi`. It is recommended to provide a user id of coder who are going to be doing the test.

```{r createtest}
oolong_test <- wi(abstracts_keyatm, userid = "Hadley")
oolong_test
```

As instructed, use the method `$do_word_intrusion_test()` to start coding. 

```{r, eval = FALSE}
oolong_test$do_word_intrusion_test()
```

After the coding, you need to first lock the test. Then, you can look at the model precision by printing the oolong test.

```{r, include = FALSE}
### Mock this process
oolong_test$.__enclos_env__$private$test_content$wi$answer <- oolong_test$.__enclos_env__$private$test_content$wi$intruder
oolong_test$.__enclos_env__$private$test_content$wi$answer[1] <- "wronganswer"
```

```{r lock}
oolong_test$lock()
oolong_test
```

#### Word set intrusion test

Word set intrusion test is a variant of word intrusion test (Ying et al., 2021), in which multiple word sets generated from top terms of one topic are juxtaposed with one intruder word set generated similarly from another topic. In Ying et al., this test is called "R4WSI" because 4 word sets are displayed. By default, oolong generates also R4WSI. However, it is also possible to generate R(N)WSI by setting the parameter `n_correct_ws` to N - 1.

```{r wsi1}
oolong_test <- wsi(abstracts_keyatm, userid = "Garrett")
oolong_test
```

Use the method `$do_word_set_intrusion_test()` to start coding. 

```{r wsi2, eval = FALSE}
oolong_test$do_word_set_intrusion_test()
```

```{r, include = FALSE}
### Mock this process
oolong_test$.__enclos_env__$private$test_content$wsi$answer <- oolong_test$.__enclos_env__$private$test_content$wsi$intruder
oolong_test$.__enclos_env__$private$test_content$wsi$answer[1] <- "wronganswer"
```

```{r wsi3}
oolong_test$lock()
oolong_test
```

#### Topic intrusion test

For example, `abstracts_keyatm` was generated with the corpus `abstracts$text`

```{r newgroup5}
library(tibble)
abstracts
```

Creating the oolong test object with the corpus used for training the topic model will generate topic intrusion test cases. 

```{r createtest2}
oolong_test <- ti(abstracts_keyatm, abstracts$text, userid = "Julia")
oolong_test
```

Similarly, use the `$do_topic_intrusion_test` to code the test cases, lock the test with `$lock()` and then you can look at the TLO (topic log odds) value by printing the oolong test.

```{r, eval = FALSE}
oolong_test$do_topic_intrusion_test()
oolong_test$lock()
```

```{r, include = FALSE}
genius_topic <- function(obj1) {
    obj1$.__enclos_env__$private$test_content$ti$answer <- obj1$.__enclos_env__$private$test_content$ti$intruder
    return(obj1)
}
genius_word <- function(obj1) {
    obj1$.__enclos_env__$private$test_content$wi$answer <- obj1$.__enclos_env__$private$test_content$wi$intruder
    return(obj1)
}
oolong_test <- genius_word(genius_topic(oolong_test))
oolong_test$.__enclos_env__$private$test_content$ti$answer[2] <- sample(oolong_test$.__enclos_env__$private$test_content$ti$candidates[[2]], 1)
oolong_test$lock()
```

```{r topic_res}
oolong_test
```

### Suggested workflow

The test makes more sense if more than one coder is involved. A suggested workflow is to create the test, then clone the oolong object. Ask multiple coders to do the test(s) and then summarize the results.

Preprocess and create a document-feature matrix

```{r, eval = FALSE}
dfm(abstracts$text, tolower = TRUE, stem = TRUE, remove = stopwords('english'), remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE, remove_hyphens = TRUE) %>% dfm_trim(min_docfreq = 3, max_docfreq = 500) %>% dfm_select(min_nchar = 3, pattern = "^[a-zA-Z]+$", valuetype = "regex") -> abstracts_dfm
```

Train a topic model.

```{r step0, eval = FALSE}
require(keyATM)
abstracts_keyatm <- keyATM(keyATM_read(abstracts_dfm), no_keyword_topics = 0, keywords = abstracts_dictionary, model = "base", options = list(seed = 46709394))
```

Create a new oolong object.

```{r step1}
oolong_test_rater1 <- witi(abstracts_keyatm, abstracts$text, userid = "Yihui")
```

Clone the oolong object to be used by other raters.

```{r step2}
oolong_test_rater2 <- clone_oolong(oolong_test_rater1, userid = "Jenny")
```

Ask different coders to code each object and then lock the object.

```{r, eval = FALSE}
## Let Yihui do the test.
oolong_test_rater1$do_word_intrusion_test()
oolong_test_rater1$do_topic_intrusion_test()
oolong_test_rater1$lock()

## Let Jenny do the test.
oolong_test_rater2$do_word_intrusion_test()
oolong_test_rater2$do_topic_intrusion_test()
oolong_test_rater2$lock()
```

```{r, include = FALSE}
### Mock this process
set.seed(46709394)
oolong_test_rater1 <- oolong:::.monkey_test(oolong_test_rater1, intelligent = 0.3)
oolong_test_rater2 <- oolong:::.monkey_test(oolong_test_rater2, intelligent = 0)
oolong_test_rater1$lock()
oolong_test_rater2$lock()
```

Get a summary of the two objects.

```{r, step3}
summarize_oolong(oolong_test_rater1, oolong_test_rater2)
```

### About the p-values

The test for model precision (MP) is based on an one-tailed, one-sample binomial test for each rater. In a multiple-rater situation, the p-values from all raters are combined using the Fisher's method (a.k.a. Fisher's omnibus test).

H0: MP is not better than 1/ n\_top\_terms

H1: MP is better than 1/ n\_top\_terms


The test for the median of TLO is based on a permutation test.

H0: Median TLO is not better than random guess.

H1: Median TLO is better than random guess.

One must notice that the two statistical tests are testing the bear minimum. A significant test only indicates the topic model can make the rater(s) perform better than random guess. It is not an indication of good topic interpretability. Also, one should use a very conservative significant level, e.g. $\alpha < 0.001$.

### About Warp LDA

There is a subtle difference between the support for `stm` and for `text2vec`.

`abstracts_warplda` is a Warp LDA object trained with the same dataset as the `abstracts_stm`

```{r warplda}
abstracts_warplda
```

All the API endpoints are the same, except the one for the creation of topic intrusion test cases. You must supply also the `input_dfm`.

```{r warplda2}
### Just word intrusion test.
oolong_test <- wi(abstracts_warplda, userid = "Lionel")
oolong_test
```


```{r warplda3}
abstracts_dfm
```

```{r warplda4, , message = FALSE, results = 'hide', warning = FALSE}
oolong_test <- witi(abstracts_warplda, abstracts$text, input_dfm = abstracts_dfm, userid = "Mara")
```

```{r warplda5}
oolong_test
```

## About Biterm Topic Model

Please refer to the vignette about BTM.

## About Naive Bayes

Naive Bayes model is a supervised machine learning model. This package supports Naive Bayes models trained using `quanteda.textmodels`.

Suppose `newsgroup_nb` is a Naive Bayes model trained on a subset of the classic [20 newsgroups] dataset.

```r
tokens(newsgroup5$text, remove_punct = TRUE, remove_symbols = TRUE, remove_numbers = TRUE, remove_url = TRUE, spilit_hyphens = TRUE) %>% tokens_wordstem %>% tokens_remove(stopwords("en")) %>% dfm(tolower = TRUE) %>% dfm_trim(min_termfreq = 3, max_docfreq = 0.06, docfreq_type = "prop") -> newsgroup_dfm
docvars(newsgroup_dfm, "group") <- newsgroup5$title
newsgroup_nb <- textmodel_nb(newsgroup_dfm, docvars(newsgroup_dfm, "group"), distribution = "Bernoulli")
```

You can still generate word intrusion and word set intrusion tests.

```{r}
wi(newsgroup_nb)
```

```{r}
wsi(newsgroup_nb)
```

## Validating Dictionary-based Methods

### Creating gold standard

`trump2k` is a dataset of 2,000 tweets from \@realdonaldtrump.

```{r trump2k}
tibble(text = trump2k)
```

For example, you are interested in studying the sentiment of these tweets. One can use tools such as AFINN to automatically extract sentiment in these tweets. However, oolong recommends to generate gold standard by human coding first using a subset. By default, oolong selects 1% of the origin corpus as test cases. The parameter `construct` should be an adjective, e.g. positive, liberal, populistic, etc.

```{r goldstandard}
oolong_test <- gs(input_corpus = trump2k, construct = "positive", userid = "Joe")
oolong_test
```

As instructed, use the method `$do_gold_standard_test()` to start coding.

```{r, eval = FALSE}
oolong_test$do_gold_standard_test()
```

After the coding, you need to first lock the test and then the `$turn_gold()` method is available.

```{r, include = FALSE}
oolong_test$.__enclos_env__$private$test_content$gs <- 
structure(list(case = 1:20, text = c("Thank you Eau Claire, Wisconsin. \n#VoteTrump on Tuesday, April 5th!\nMAKE AMERICA GREAT AGAIN! https://t.co/JI5JqwHnMC", 
"\"@bobby990r_1: @realDonaldTrump would lead polls the second he announces candidacy! America is waiting for him to LEAD us out of this mess!", 
"\"@KdanielsK: @misstcassidy @AllAboutTheTea_ @realDonaldTrump My money is on Kenya getting fired first.\"", 
"Thank you for a great afternoon Birmingham, Alabama! #Trump2016 #MakeAmericaGreatAgain https://t.co/FrOkqCzBoD", 
"\"@THETAINTEDT: @foxandfriends @realDonaldTrump Trump 2016 http://t.co/UlQWGKUrCJ\"", 
"People believe CNN these days almost as little as they believe Hillary....that's really saying something!", 
"It was great being in Michigan. Remember, I am the only presidential candidate who will bring jobs back to the U.S.and protect car industry!", 
"\"@DomineekSmith: @realDonaldTrump is the best Republican presidential candidate of all time.\"  Thank you.", 
"Word is that little Morty Zuckerman’s @NYDailyNews loses more than $50 million per year---can that be possible?", 
"\"@Chevy_Mama: @realDonaldTrump I'm obsessed with @celebrityapprenticeNBC. Honestly,  Mr Trump, you are very inspiring\"", 
"President Obama said \"ISIL continues to shrink\" in an interview just hours before the horrible attack in Paris. He is just so bad! CHANGE.", 
".@HillaryClinton loves to lie. America has had enough of the CLINTON'S! It is time to #DrainTheSwamp! Debates https://t.co/3Mz4T7qTTR", 
"\"@jerrimoore: @realDonaldTrump awesome. A treat to get to see the brilliant Joan Rivers once more #icon\"", 
"\"@shoegoddesss: @realDonaldTrump Will definitely vote for you. Breath of fresh air. America needs you!\"", 
"Ted is the ultimate hypocrite. Says one thing for money, does another for votes. \nhttps://t.co/hxdfy0mjVw", 
"\"@Lisa_Milicaj: Truth be told, I  never heard of The National Review until they \"tried\" to declare war on you. No worries, you got my vote!\"", 
"THANK YOU Daytona Beach, Florida!\n#MakeAmericaGreatAgain https://t.co/IAcLfXe463", 
"People rarely say that many conservatives didn’t vote for Mitt Romney. If I can get them to vote for me, we win in a landslide.", 
"Trump National Golf Club, Washington, D.C. is on 600 beautiful acres fronting the Potomac River. A fantastic setting! http://t.co/pYtkbyKwt5", 
"\"@DRUDGE_REPORT: REUTERS 5-DAY ROLLING POLL: TRUMP 34%, CARSON 19.6%, RUBIO 9.7%, CRUZ 7.7%...\" Thank you - a great honor!"
), answer = c(4L, 4L, 2L, 5L, 3L, 2L, 4L, 5L, 2L, 4L, 1L, 1L, 
4L, 4L, 2L, 4L, 4L, 4L, 4L, 4L), target_value = c(NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
NA)), row.names = c(NA, -20L), class = c("tbl_df", "tbl", "data.frame"
))
```

```{r gs_locking}
oolong_test$lock()
oolong_test
```

### Example: Validating AFINN using the gold standard

A locked oolong test can be converted into a quanteda-compatible corpus for further analysis. The corpus contains two `docvars`, 'answer'.

```{r}
oolong_test$turn_gold()
```

In this example, we calculate the AFINN score for each tweet using quanteda. The dictionary `afinn` is bundle with this package.

```{r}
gold_standard <- oolong_test$turn_gold()
dfm(gold_standard, remove_punct = TRUE) %>% dfm_lookup(afinn) %>% quanteda::convert(to = "data.frame") %>%
    mutate(matching_word_valence = (neg5 * -5) + (neg4 * -4) + (neg3 * -3) + (neg2 * -2) + (neg1 * -1)
           + (zero * 0) + (pos1 * 1) + (pos2 * 2) + (pos3 * 3) + (pos4 * 4) + (pos5 * 5),
           base = ntoken(gold_standard, remove_punct = TRUE), afinn_score = matching_word_valence / base) %>%
    pull(afinn_score) -> all_afinn_score
all_afinn_score
```

Put back the vector of AFINN score into the respective `docvars` and study the correlation between the gold standard and AFINN.

```{r}
summarize_oolong(oolong_test, target_value = all_afinn_score)
```

### Suggested workflow

Create an oolong object, clone it for another coder. According to Song et al. (Forthcoming), you should at least draw 1% of your data.

```{r}
trump <- gs(input_corpus = trump2k, exact_n = 40, userid = "JJ")
trump2 <- clone_oolong(trump, userid = "Winston")
```

Instruct two coders to code the tweets and lock the objects.

```{r, eval = FALSE}
trump$do_gold_standard_test()
trump2$do_gold_standard_test()
trump$lock()
trump2$lock()
```

```{r, include = FALSE}
trump$.__enclos_env__$private$test_content$gs <- 
structure(list(case = 1:20, text = c("Thank you Eau Claire, Wisconsin. \n#VoteTrump on Tuesday, April 5th!\nMAKE AMERICA GREAT AGAIN! https://t.co/JI5JqwHnMC", 
"\"@bobby990r_1: @realDonaldTrump would lead polls the second he announces candidacy! America is waiting for him to LEAD us out of this mess!", 
"\"@KdanielsK: @misstcassidy @AllAboutTheTea_ @realDonaldTrump My money is on Kenya getting fired first.\"", 
"Thank you for a great afternoon Birmingham, Alabama! #Trump2016 #MakeAmericaGreatAgain https://t.co/FrOkqCzBoD", 
"\"@THETAINTEDT: @foxandfriends @realDonaldTrump Trump 2016 http://t.co/UlQWGKUrCJ\"", 
"People believe CNN these days almost as little as they believe Hillary....that's really saying something!", 
"It was great being in Michigan. Remember, I am the only presidential candidate who will bring jobs back to the U.S.and protect car industry!", 
"\"@DomineekSmith: @realDonaldTrump is the best Republican presidential candidate of all time.\"  Thank you.", 
"Word is that little Morty Zuckerman’s @NYDailyNews loses more than $50 million per year---can that be possible?", 
"\"@Chevy_Mama: @realDonaldTrump I'm obsessed with @celebrityapprenticeNBC. Honestly,  Mr Trump, you are very inspiring\"", 
"President Obama said \"ISIL continues to shrink\" in an interview just hours before the horrible attack in Paris. He is just so bad! CHANGE.", 
".@HillaryClinton loves to lie. America has had enough of the CLINTON'S! It is time to #DrainTheSwamp! Debates https://t.co/3Mz4T7qTTR", 
"\"@jerrimoore: @realDonaldTrump awesome. A treat to get to see the brilliant Joan Rivers once more #icon\"", 
"\"@shoegoddesss: @realDonaldTrump Will definitely vote for you. Breath of fresh air. America needs you!\"", 
"Ted is the ultimate hypocrite. Says one thing for money, does another for votes. \nhttps://t.co/hxdfy0mjVw", 
"\"@Lisa_Milicaj: Truth be told, I  never heard of The National Review until they \"tried\" to declare war on you. No worries, you got my vote!\"", 
"THANK YOU Daytona Beach, Florida!\n#MakeAmericaGreatAgain https://t.co/IAcLfXe463", 
"People rarely say that many conservatives didn’t vote for Mitt Romney. If I can get them to vote for me, we win in a landslide.", 
"Trump National Golf Club, Washington, D.C. is on 600 beautiful acres fronting the Potomac River. A fantastic setting! http://t.co/pYtkbyKwt5", 
"\"@DRUDGE_REPORT: REUTERS 5-DAY ROLLING POLL: TRUMP 34%, CARSON 19.6%, RUBIO 9.7%, CRUZ 7.7%...\" Thank you - a great honor!"
), answer = c(4L, 4L, 2L, 5L, 3L, 2L, 4L, 5L, 2L, 4L, 1L, 1L, 
4L, 4L, 2L, 4L, 4L, 4L, 4L, 4L), target_value = c(NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
NA)), row.names = c(NA, -20L), class = c("tbl_df", "tbl", "data.frame"
                                         ))

trump2$.__enclos_env__$private$test_content$gs <- 
structure(list(case = 1:20, text = c("Thank you Eau Claire, Wisconsin. \n#VoteTrump on Tuesday, April 5th!\nMAKE AMERICA GREAT AGAIN! https://t.co/JI5JqwHnMC", 
"\"@bobby990r_1: @realDonaldTrump would lead polls the second he announces candidacy! America is waiting for him to LEAD us out of this mess!", 
"\"@KdanielsK: @misstcassidy @AllAboutTheTea_ @realDonaldTrump My money is on Kenya getting fired first.\"", 
"Thank you for a great afternoon Birmingham, Alabama! #Trump2016 #MakeAmericaGreatAgain https://t.co/FrOkqCzBoD", 
"\"@THETAINTEDT: @foxandfriends @realDonaldTrump Trump 2016 http://t.co/UlQWGKUrCJ\"", 
"People believe CNN these days almost as little as they believe Hillary....that's really saying something!", 
"It was great being in Michigan. Remember, I am the only presidential candidate who will bring jobs back to the U.S.and protect car industry!", 
"\"@DomineekSmith: @realDonaldTrump is the best Republican presidential candidate of all time.\"  Thank you.", 
"Word is that little Morty Zuckerman’s @NYDailyNews loses more than $50 million per year---can that be possible?", 
"\"@Chevy_Mama: @realDonaldTrump I'm obsessed with @celebrityapprenticeNBC. Honestly,  Mr Trump, you are very inspiring\"", 
"President Obama said \"ISIL continues to shrink\" in an interview just hours before the horrible attack in Paris. He is just so bad! CHANGE.", 
".@HillaryClinton loves to lie. America has had enough of the CLINTON'S! It is time to #DrainTheSwamp! Debates https://t.co/3Mz4T7qTTR", 
"\"@jerrimoore: @realDonaldTrump awesome. A treat to get to see the brilliant Joan Rivers once more #icon\"", 
"\"@shoegoddesss: @realDonaldTrump Will definitely vote for you. Breath of fresh air. America needs you!\"", 
"Ted is the ultimate hypocrite. Says one thing for money, does another for votes. \nhttps://t.co/hxdfy0mjVw", 
"\"@Lisa_Milicaj: Truth be told, I  never heard of The National Review until they \"tried\" to declare war on you. No worries, you got my vote!\"", 
"THANK YOU Daytona Beach, Florida!\n#MakeAmericaGreatAgain https://t.co/IAcLfXe463", 
"People rarely say that many conservatives didn’t vote for Mitt Romney. If I can get them to vote for me, we win in a landslide.", 
"Trump National Golf Club, Washington, D.C. is on 600 beautiful acres fronting the Potomac River. A fantastic setting! http://t.co/pYtkbyKwt5", 
"\"@DRUDGE_REPORT: REUTERS 5-DAY ROLLING POLL: TRUMP 34%, CARSON 19.6%, RUBIO 9.7%, CRUZ 7.7%...\" Thank you - a great honor!"
), answer = c(5L, 3L, 2L, 5L, 3L, 1L, 4L, 5L, 2L, 4L, 1L, 1L, 
4L, 4L, 2L, 4L, 4L, 4L, 4L, 4L), target_value = c(NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
NA)), row.names = c(NA, -20L), class = c("tbl_df", "tbl", "data.frame"
                                         ))
trump$lock()
trump2$lock()
```

Calculate the target value (in this case, the AFINN score) by turning one object into a corpus.

```{r}
gold_standard <- trump$turn_gold()
dfm(gold_standard, remove_punct = TRUE) %>% dfm_lookup(afinn) %>% quanteda::convert(to = "data.frame") %>%
    mutate(matching_word_valence = (neg5 * -5) + (neg4 * -4) + (neg3 * -3) + (neg2 * -2) + (neg1 * -1)
           + (zero * 0) + (pos1 * 1) + (pos2 * 2) + (pos3 * 3) + (pos4 * 4) + (pos5 * 5),
           base = ntoken(gold_standard, remove_punct = TRUE), afinn_score = matching_word_valence / base) %>%
    pull(afinn_score) -> target_value
```

Summarize all oolong objects with the target value.

```{r}
res <- summarize_oolong(trump, trump2, target_value = target_value)
```

Read the results. The diagnostic plot consists of 4 subplots. It is a good idea to read Bland & Altman (1986) on the difference between correlation and agreement.

* Subplot (top left): Raw correlation between human judgement and target value. One should want to have a good correlation between the two.
* Subplot (top right): Bland-Altman plot. One should want to have no correlation. Also, the dots should be randomly scattering around the mean value. If it is so, the two measurements (human judgement and target value) are in good agreement.
* Subplot (bottom left): Raw correlation between target value and content length. One should want to have no correlation, as an indication of good reliability against the influence of content length. (See Chan et al. 2020)
* Subplot (bottom right): Cook's distance of all data point. One should want to have no dot (or at least very few dots) above the threshold. It is an indication of how the raw correlation between human judgement and target value can or cannot be influenced by extreme values in your data. 

The textual output contains the Krippendorff's alpha of the codings by your raters. In order to claim validity of your target value, you must first establish the reliability of your gold standard. Song et al. [Forthcoming] suggest Krippendorff's Alpha > 0.7 as an acceptable cut-off.

```{r}
res
```

```{r diagnosis}
plot(res)
```

## Backward compatibility

Historically, oolong test objects could only be generated with only one function: `create_oolong`. It is no longer the case and no longer recommended anymore. It is still retained for backward compatibility purposes. If you still need to use  `create_oolong()`, the most important parameters are `input_model` and `input_corpus`. Setting each of them to `NULL` generates different tests.

| input\_model | input\_corpus | output                                                                                                                                      |
|--------------|:-------------:|---------------------------------------------------------------------------------------------------------------------------------------------|
| Not NULL     | NULL          | oolong test for validating a topic model with [word intrusion test](#word-intrusion-test)                                                   |
| Not NULL     | Not NULL      | oolong test for validating a topic model with [word intrusion test](#word-intrusion-test) and [topic intrusion test](#topic-intrusion-test) |
| NULL         | Not NULL      | oolong test for [creating gold standard](#creating-gold-standard)                                                                           |
| NULL         | NULL          | error                                                                                                                                       |


## References

1. Chang, J., Gerrish, S., Wang, C., Boyd-Graber, J. L., & Blei, D. M. (2009). Reading tea leaves: How humans interpret topic models. In Advances in neural information processing systems (pp. 288-296). [link](https://papers.nips.cc/paper/3700-reading-tea-leaves-how-humans-interpret-topic-models)
2. Ying, L., Montgomery, J. M., & Stewart, B. M. (2021). Inferring concepts from topics: Towards procedures for validating topics as measures. Political Analysis. [link](https://doi.org/10.1017/pan.2021.33)
3. Song et al. (2020) In validations we trust? The impact of imperfect human annotations as a gold standard on the quality of validation of automated content analysis. Political Communication. [link](https://doi.org/10.1080/10584609.2020.1723752)
4. Bland, J. M., & Altman, D. (1986). Statistical methods for assessing agreement between two methods of clinical measurement. The lancet, 327(8476), 307-310.
5. Chan et al. (2020) Four best practices for measuring news sentiment using ‘off-the-shelf’ dictionaries: a large-scale p-hacking experiment. Computational Communication Research. [link](https://osf.io/preprints/socarxiv/np5wa/)
6. Nielsen, F. Å. (2011). A new ANEW: Evaluation of a word list for sentiment analysis in microblogs. arXiv preprint arXiv:1103.2903. [link](https://arxiv.org/abs/1103.2903)

---
---
title: "BTM"
output: rmarkdown::html_vignette
author:
  - Chung-hong Chan ^[University of Mannheim]
vignette: >
  %\VignetteIndexEntry{BTM}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  )
set.seed(46709394)
```


The package BTM by Jan Wijffels et al. finds "topics in collections of short text". Compared to other topic model packages, BTM requires a special data format for training. Oolong has no problem generating word intrusion tests with BTM. However, that special data format can make creation of topic intrusion tests very tricky.

This guide provides our recommendations on how to use BTM, so that the model can be used for generating topic intrusion tests.

# Requirement #1: Keep your quanteda corpus

It is because every document has a unique document id.

```{r}
require(BTM)
require(quanteda)
require(oolong)
trump_corpus <- corpus(trump2k)
```

And then you can do regular text cleaning, stemming procedure with `quanteda`. Instead of making the product a `DFM` object, make it a `token` object. You may read [this issue](https://github.com/quanteda/quanteda/issues/1404) by Benoit et al.

```{r}
tokens(trump_corpus, remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE, split_hyphens = TRUE, remove_url = TRUE) %>% tokens_tolower() %>% tokens_remove(stopwords("en")) %>% tokens_remove("@*")  -> trump_toks
```

# Requirement #2: Keep your data frame

Use this function to convert the `token` object to a data frame.

```{r}
as.data.frame.tokens <- function(x) {
  data.frame(
    doc_id = rep(names(x), lengths(x)),
    tokens = unlist(x, use.names = FALSE)
  )
}

trump_dat <- as.data.frame.tokens(trump_toks)
```

Train a BTM model

```{r, message = FALSE, results = 'hide', warning = FALSE}
trump_btm <- BTM(trump_dat, k = 8, iter = 500, trace = 10)
```

## Pecularities of BTM

This is how you should generate $\theta_{t}$ . However, there are many NaN and there are only 1994 rows (`trump2k` has 2000 tweets) due to empty documents.

```{r}
theta <- predict(trump_btm, newdata = trump_dat)
dim(theta)
```

```{r}
setdiff(docid(trump_corpus), row.names(theta))
```

```{r}
trump_corpus[604]
```

Also, the row order is messed up.

```{r}
head(row.names(theta), 100)
```


# Oolong's support for BTM

Oolong has no problem generating word intrusion test for BTM like you do with other topic models.

```{r}
oolong <- create_oolong(trump_btm)
oolong
```

For generating topic intrusion tests, however, you must provide the data frame you used for training (in this case `trump_dat`). Your `input_corpus` must be a quanteda corpus too.

```{r}
oolong <- create_oolong(trump_btm, trump_corpus, btm_dataframe = trump_dat)
oolong
```

`btm_dataframe` must not be NULL.

```{r, error = TRUE}
oolong <- create_oolong(trump_btm, trump_corpus)
```

`input_corpus` must be a quanteda corpus.

```{r, error = TRUE}
oolong <- create_oolong(trump_btm, trump2k, btm_dataframe = trump_dat)
```
---
title: "Overview"
output: rmarkdown::html_vignette
author:
  - Chung-hong Chan ^[University of Mannheim]
vignette: >
  %\VignetteIndexEntry{Overview}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  )
set.seed(46709394)
```

The validation test is called "oolong test" (for reading tea leaves). This package provides several functions for generating different types of oolong test.

| function | purpose                                                                                                                           |
|---------:|:----------------------------------------------------------------------------------------------------------------------------------|
|   `wi()` | validating a topic model with [word intrusion test](#word-intrusion-test) (Chang et al., 2008)                                    |
|   `ti()` | validating a topic model with [topic intrusion test](#topic-intrusion-test) (Chang et al., 2008; aka "T8WSI" in Ying et al. 2021) |
| `witi()` | validating a topic model with [word intrusion test](#word-intrusion-test) and [topic intrusion test](#topic-intrusion-test)       |
|  `wsi()` | validating a topic model with [word set intrusion test](#word-set-intrusion-test) (Ying et al. 2021)                              |
|   `gs()` | oolong test for [creating gold standard](#creating-gold-standard) (see Song et al., 2020)                                         |

All of these tests can also be generated with the function [`create_oolong`](#backward-compatibility). As of version 0.3.20, it is no longer recommended.

## Installation

Because the package is constantly changing, we suggest using the development version from [GitHub](https://github.com/):

``` r
# install.packages("devtools")
devtools::install_github("chainsawriot/oolong")
```

You can also install the "stable" (but slightly older) version from CRAN:

```r
install.packages("oolong")
```

## Validating Topic Models

#### Word intrusion test

`abstracts_keyatm` is an example topic model trained with the data `abstracts` using the `keyATM` package. Currently, this package supports structural topic models / correlated topic models from `stm`, Warp LDA models from `text2vec` , LDA/CTM models from `topicmodels`,  Biterm Topic Models from `BTM`, Keyword Assisted Topic Models from `keyATM`, and seeded LDA models from `seededlda`. Although not strictly a topic model, Naive Bayes models from `quanteda.textmodels` are also supported. See the section on [Naive Bayes](#about-naive-bayes) for more information.

```{r}
library(oolong)
library(keyATM)
library(quanteda)
library(dplyr)
```

```{r example}
abstracts_keyatm
```

To create an oolong test with word intrusion test, use the function `wi`. It is recommended to provide a user id of coder who are going to be doing the test.

```{r createtest}
oolong_test <- wi(abstracts_keyatm, userid = "Hadley")
oolong_test
```

As instructed, use the method `$do_word_intrusion_test()` to start coding. 

```{r, eval = FALSE}
oolong_test$do_word_intrusion_test()
```

After the coding, you need to first lock the test. Then, you can look at the model precision by printing the oolong test.

```{r, include = FALSE}
### Mock this process
oolong_test$.__enclos_env__$private$test_content$wi$answer <- oolong_test$.__enclos_env__$private$test_content$wi$intruder
oolong_test$.__enclos_env__$private$test_content$wi$answer[1] <- "wronganswer"
```

```{r lock}
oolong_test$lock()
oolong_test
```

#### Word set intrusion test

Word set intrusion test is a variant of word intrusion test (Ying et al., 2021), in which multiple word sets generated from top terms of one topic are juxtaposed with one intruder word set generated similarly from another topic. In Ying et al., this test is called "R4WSI" because 4 word sets are displayed. By default, oolong generates also R4WSI. However, it is also possible to generate R(N)WSI by setting the parameter `n_correct_ws` to N - 1.

```{r wsi1}
oolong_test <- wsi(abstracts_keyatm, userid = "Garrett")
oolong_test
```

Use the method `$do_word_set_intrusion_test()` to start coding. 

```{r wsi2, eval = FALSE}
oolong_test$do_word_set_intrusion_test()
```

```{r, include = FALSE}
### Mock this process
oolong_test$.__enclos_env__$private$test_content$wsi$answer <- oolong_test$.__enclos_env__$private$test_content$wsi$intruder
oolong_test$.__enclos_env__$private$test_content$wsi$answer[1] <- "wronganswer"
```

```{r wsi3}
oolong_test$lock()
oolong_test
```

#### Topic intrusion test

For example, `abstracts_keyatm` was generated with the corpus `abstracts$text`

```{r newgroup5}
library(tibble)
abstracts
```

Creating the oolong test object with the corpus used for training the topic model will generate topic intrusion test cases. 

```{r createtest2}
oolong_test <- ti(abstracts_keyatm, abstracts$text, userid = "Julia")
oolong_test
```

Similarly, use the `$do_topic_intrusion_test` to code the test cases, lock the test with `$lock()` and then you can look at the TLO (topic log odds) value by printing the oolong test.

```{r, eval = FALSE}
oolong_test$do_topic_intrusion_test()
oolong_test$lock()
```

```{r, include = FALSE}
genius_topic <- function(obj1) {
    obj1$.__enclos_env__$private$test_content$ti$answer <- obj1$.__enclos_env__$private$test_content$ti$intruder
    return(obj1)
}
genius_word <- function(obj1) {
    obj1$.__enclos_env__$private$test_content$wi$answer <- obj1$.__enclos_env__$private$test_content$wi$intruder
    return(obj1)
}
oolong_test <- genius_word(genius_topic(oolong_test))
oolong_test$.__enclos_env__$private$test_content$ti$answer[2] <- sample(oolong_test$.__enclos_env__$private$test_content$ti$candidates[[2]], 1)
oolong_test$lock()
```

```{r topic_res}
oolong_test
```

### Suggested workflow

The test makes more sense if more than one coder is involved. A suggested workflow is to create the test, then clone the oolong object. Ask multiple coders to do the test(s) and then summarize the results.

Preprocess and create a document-feature matrix

```{r, eval = FALSE}
dfm(abstracts$text, tolower = TRUE, stem = TRUE, remove = stopwords('english'), remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE, remove_hyphens = TRUE) %>% dfm_trim(min_docfreq = 3, max_docfreq = 500) %>% dfm_select(min_nchar = 3, pattern = "^[a-zA-Z]+$", valuetype = "regex") -> abstracts_dfm
```

Train a topic model.

```{r step0, eval = FALSE}
require(keyATM)
abstracts_keyatm <- keyATM(keyATM_read(abstracts_dfm), no_keyword_topics = 0, keywords = abstracts_dictionary, model = "base", options = list(seed = 46709394))
```

Create a new oolong object.

```{r step1}
oolong_test_rater1 <- witi(abstracts_keyatm, abstracts$text, userid = "Yihui")
```

Clone the oolong object to be used by other raters.

```{r step2}
oolong_test_rater2 <- clone_oolong(oolong_test_rater1, userid = "Jenny")
```

Ask different coders to code each object and then lock the object.

```{r, eval = FALSE}
## Let Yihui do the test.
oolong_test_rater1$do_word_intrusion_test()
oolong_test_rater1$do_topic_intrusion_test()
oolong_test_rater1$lock()

## Let Jenny do the test.
oolong_test_rater2$do_word_intrusion_test()
oolong_test_rater2$do_topic_intrusion_test()
oolong_test_rater2$lock()
```

```{r, include = FALSE}
### Mock this process
set.seed(46709394)
oolong_test_rater1 <- oolong:::.monkey_test(oolong_test_rater1, intelligent = 0.3)
oolong_test_rater2 <- oolong:::.monkey_test(oolong_test_rater2, intelligent = 0)
oolong_test_rater1$lock()
oolong_test_rater2$lock()
```

Get a summary of the two objects.

```{r, step3}
summarize_oolong(oolong_test_rater1, oolong_test_rater2)
```

### About the p-values

The test for model precision (MP) is based on an one-tailed, one-sample binomial test for each rater. In a multiple-rater situation, the p-values from all raters are combined using the Fisher's method (a.k.a. Fisher's omnibus test).

H0: MP is not better than 1/ n\_top\_terms

H1: MP is better than 1/ n\_top\_terms


The test for the median of TLO is based on a permutation test.

H0: Median TLO is not better than random guess.

H1: Median TLO is better than random guess.

One must notice that the two statistical tests are testing the bear minimum. A significant test only indicates the topic model can make the rater(s) perform better than random guess. It is not an indication of good topic interpretability. Also, one should use a very conservative significant level, e.g. $\alpha < 0.001$.

### About Warp LDA

There is a subtle difference between the support for `stm` and for `text2vec`.

`abstracts_warplda` is a Warp LDA object trained with the same dataset as the `abstracts_stm`

```{r warplda}
abstracts_warplda
```

All the API endpoints are the same, except the one for the creation of topic intrusion test cases. You must supply also the `input_dfm`.

```{r warplda2}
### Just word intrusion test.
oolong_test <- wi(abstracts_warplda, userid = "Lionel")
oolong_test
```


```{r warplda3}
abstracts_dfm
```

```{r warplda4, , message = FALSE, results = 'hide', warning = FALSE}
oolong_test <- witi(abstracts_warplda, abstracts$text, input_dfm = abstracts_dfm, userid = "Mara")
```

```{r warplda5}
oolong_test
```

## About Biterm Topic Model

Please refer to the vignette about BTM.

## About Naive Bayes

Naive Bayes model is a supervised machine learning model. This package supports Naive Bayes models trained using `quanteda.textmodels`.

Suppose `newsgroup_nb` is a Naive Bayes model trained on a subset of the classic [20 newsgroups] dataset.

```r
tokens(newsgroup5$text, remove_punct = TRUE, remove_symbols = TRUE, remove_numbers = TRUE, remove_url = TRUE, spilit_hyphens = TRUE) %>% tokens_wordstem %>% tokens_remove(stopwords("en")) %>% dfm(tolower = TRUE) %>% dfm_trim(min_termfreq = 3, max_docfreq = 0.06, docfreq_type = "prop") -> newsgroup_dfm
docvars(newsgroup_dfm, "group") <- newsgroup5$title
newsgroup_nb <- textmodel_nb(newsgroup_dfm, docvars(newsgroup_dfm, "group"), distribution = "Bernoulli")
```

You can still generate word intrusion and word set intrusion tests.

```{r}
wi(newsgroup_nb)
```

```{r}
wsi(newsgroup_nb)
```

## Validating Dictionary-based Methods

### Creating gold standard

`trump2k` is a dataset of 2,000 tweets from \@realdonaldtrump.

```{r trump2k}
tibble(text = trump2k)
```

For example, you are interested in studying the sentiment of these tweets. One can use tools such as AFINN to automatically extract sentiment in these tweets. However, oolong recommends to generate gold standard by human coding first using a subset. By default, oolong selects 1% of the origin corpus as test cases. The parameter `construct` should be an adjective, e.g. positive, liberal, populistic, etc.

```{r goldstandard}
oolong_test <- gs(input_corpus = trump2k, construct = "positive", userid = "Joe")
oolong_test
```

As instructed, use the method `$do_gold_standard_test()` to start coding.

```{r, eval = FALSE}
oolong_test$do_gold_standard_test()
```

After the coding, you need to first lock the test and then the `$turn_gold()` method is available.

```{r, include = FALSE}
oolong_test$.__enclos_env__$private$test_content$gs <- 
structure(list(case = 1:20, text = c("Thank you Eau Claire, Wisconsin. \n#VoteTrump on Tuesday, April 5th!\nMAKE AMERICA GREAT AGAIN! https://t.co/JI5JqwHnMC", 
"\"@bobby990r_1: @realDonaldTrump would lead polls the second he announces candidacy! America is waiting for him to LEAD us out of this mess!", 
"\"@KdanielsK: @misstcassidy @AllAboutTheTea_ @realDonaldTrump My money is on Kenya getting fired first.\"", 
"Thank you for a great afternoon Birmingham, Alabama! #Trump2016 #MakeAmericaGreatAgain https://t.co/FrOkqCzBoD", 
"\"@THETAINTEDT: @foxandfriends @realDonaldTrump Trump 2016 http://t.co/UlQWGKUrCJ\"", 
"People believe CNN these days almost as little as they believe Hillary....that's really saying something!", 
"It was great being in Michigan. Remember, I am the only presidential candidate who will bring jobs back to the U.S.and protect car industry!", 
"\"@DomineekSmith: @realDonaldTrump is the best Republican presidential candidate of all time.\"  Thank you.", 
"Word is that little Morty Zuckerman’s @NYDailyNews loses more than $50 million per year---can that be possible?", 
"\"@Chevy_Mama: @realDonaldTrump I'm obsessed with @celebrityapprenticeNBC. Honestly,  Mr Trump, you are very inspiring\"", 
"President Obama said \"ISIL continues to shrink\" in an interview just hours before the horrible attack in Paris. He is just so bad! CHANGE.", 
".@HillaryClinton loves to lie. America has had enough of the CLINTON'S! It is time to #DrainTheSwamp! Debates https://t.co/3Mz4T7qTTR", 
"\"@jerrimoore: @realDonaldTrump awesome. A treat to get to see the brilliant Joan Rivers once more #icon\"", 
"\"@shoegoddesss: @realDonaldTrump Will definitely vote for you. Breath of fresh air. America needs you!\"", 
"Ted is the ultimate hypocrite. Says one thing for money, does another for votes. \nhttps://t.co/hxdfy0mjVw", 
"\"@Lisa_Milicaj: Truth be told, I  never heard of The National Review until they \"tried\" to declare war on you. No worries, you got my vote!\"", 
"THANK YOU Daytona Beach, Florida!\n#MakeAmericaGreatAgain https://t.co/IAcLfXe463", 
"People rarely say that many conservatives didn’t vote for Mitt Romney. If I can get them to vote for me, we win in a landslide.", 
"Trump National Golf Club, Washington, D.C. is on 600 beautiful acres fronting the Potomac River. A fantastic setting! http://t.co/pYtkbyKwt5", 
"\"@DRUDGE_REPORT: REUTERS 5-DAY ROLLING POLL: TRUMP 34%, CARSON 19.6%, RUBIO 9.7%, CRUZ 7.7%...\" Thank you - a great honor!"
), answer = c(4L, 4L, 2L, 5L, 3L, 2L, 4L, 5L, 2L, 4L, 1L, 1L, 
4L, 4L, 2L, 4L, 4L, 4L, 4L, 4L), target_value = c(NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
NA)), row.names = c(NA, -20L), class = c("tbl_df", "tbl", "data.frame"
))
```

```{r gs_locking}
oolong_test$lock()
oolong_test
```

### Example: Validating AFINN using the gold standard

A locked oolong test can be converted into a quanteda-compatible corpus for further analysis. The corpus contains two `docvars`, 'answer'.

```{r}
oolong_test$turn_gold()
```

In this example, we calculate the AFINN score for each tweet using quanteda. The dictionary `afinn` is bundle with this package.

```{r}
gold_standard <- oolong_test$turn_gold()
dfm(gold_standard, remove_punct = TRUE) %>% dfm_lookup(afinn) %>% quanteda::convert(to = "data.frame") %>%
    mutate(matching_word_valence = (neg5 * -5) + (neg4 * -4) + (neg3 * -3) + (neg2 * -2) + (neg1 * -1)
           + (zero * 0) + (pos1 * 1) + (pos2 * 2) + (pos3 * 3) + (pos4 * 4) + (pos5 * 5),
           base = ntoken(gold_standard, remove_punct = TRUE), afinn_score = matching_word_valence / base) %>%
    pull(afinn_score) -> all_afinn_score
all_afinn_score
```

Put back the vector of AFINN score into the respective `docvars` and study the correlation between the gold standard and AFINN.

```{r}
summarize_oolong(oolong_test, target_value = all_afinn_score)
```

### Suggested workflow

Create an oolong object, clone it for another coder. According to Song et al. (Forthcoming), you should at least draw 1% of your data.

```{r}
trump <- gs(input_corpus = trump2k, exact_n = 40, userid = "JJ")
trump2 <- clone_oolong(trump, userid = "Winston")
```

Instruct two coders to code the tweets and lock the objects.

```{r, eval = FALSE}
trump$do_gold_standard_test()
trump2$do_gold_standard_test()
trump$lock()
trump2$lock()
```

```{r, include = FALSE}
trump$.__enclos_env__$private$test_content$gs <- 
structure(list(case = 1:20, text = c("Thank you Eau Claire, Wisconsin. \n#VoteTrump on Tuesday, April 5th!\nMAKE AMERICA GREAT AGAIN! https://t.co/JI5JqwHnMC", 
"\"@bobby990r_1: @realDonaldTrump would lead polls the second he announces candidacy! America is waiting for him to LEAD us out of this mess!", 
"\"@KdanielsK: @misstcassidy @AllAboutTheTea_ @realDonaldTrump My money is on Kenya getting fired first.\"", 
"Thank you for a great afternoon Birmingham, Alabama! #Trump2016 #MakeAmericaGreatAgain https://t.co/FrOkqCzBoD", 
"\"@THETAINTEDT: @foxandfriends @realDonaldTrump Trump 2016 http://t.co/UlQWGKUrCJ\"", 
"People believe CNN these days almost as little as they believe Hillary....that's really saying something!", 
"It was great being in Michigan. Remember, I am the only presidential candidate who will bring jobs back to the U.S.and protect car industry!", 
"\"@DomineekSmith: @realDonaldTrump is the best Republican presidential candidate of all time.\"  Thank you.", 
"Word is that little Morty Zuckerman’s @NYDailyNews loses more than $50 million per year---can that be possible?", 
"\"@Chevy_Mama: @realDonaldTrump I'm obsessed with @celebrityapprenticeNBC. Honestly,  Mr Trump, you are very inspiring\"", 
"President Obama said \"ISIL continues to shrink\" in an interview just hours before the horrible attack in Paris. He is just so bad! CHANGE.", 
".@HillaryClinton loves to lie. America has had enough of the CLINTON'S! It is time to #DrainTheSwamp! Debates https://t.co/3Mz4T7qTTR", 
"\"@jerrimoore: @realDonaldTrump awesome. A treat to get to see the brilliant Joan Rivers once more #icon\"", 
"\"@shoegoddesss: @realDonaldTrump Will definitely vote for you. Breath of fresh air. America needs you!\"", 
"Ted is the ultimate hypocrite. Says one thing for money, does another for votes. \nhttps://t.co/hxdfy0mjVw", 
"\"@Lisa_Milicaj: Truth be told, I  never heard of The National Review until they \"tried\" to declare war on you. No worries, you got my vote!\"", 
"THANK YOU Daytona Beach, Florida!\n#MakeAmericaGreatAgain https://t.co/IAcLfXe463", 
"People rarely say that many conservatives didn’t vote for Mitt Romney. If I can get them to vote for me, we win in a landslide.", 
"Trump National Golf Club, Washington, D.C. is on 600 beautiful acres fronting the Potomac River. A fantastic setting! http://t.co/pYtkbyKwt5", 
"\"@DRUDGE_REPORT: REUTERS 5-DAY ROLLING POLL: TRUMP 34%, CARSON 19.6%, RUBIO 9.7%, CRUZ 7.7%...\" Thank you - a great honor!"
), answer = c(4L, 4L, 2L, 5L, 3L, 2L, 4L, 5L, 2L, 4L, 1L, 1L, 
4L, 4L, 2L, 4L, 4L, 4L, 4L, 4L), target_value = c(NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
NA)), row.names = c(NA, -20L), class = c("tbl_df", "tbl", "data.frame"
                                         ))

trump2$.__enclos_env__$private$test_content$gs <- 
structure(list(case = 1:20, text = c("Thank you Eau Claire, Wisconsin. \n#VoteTrump on Tuesday, April 5th!\nMAKE AMERICA GREAT AGAIN! https://t.co/JI5JqwHnMC", 
"\"@bobby990r_1: @realDonaldTrump would lead polls the second he announces candidacy! America is waiting for him to LEAD us out of this mess!", 
"\"@KdanielsK: @misstcassidy @AllAboutTheTea_ @realDonaldTrump My money is on Kenya getting fired first.\"", 
"Thank you for a great afternoon Birmingham, Alabama! #Trump2016 #MakeAmericaGreatAgain https://t.co/FrOkqCzBoD", 
"\"@THETAINTEDT: @foxandfriends @realDonaldTrump Trump 2016 http://t.co/UlQWGKUrCJ\"", 
"People believe CNN these days almost as little as they believe Hillary....that's really saying something!", 
"It was great being in Michigan. Remember, I am the only presidential candidate who will bring jobs back to the U.S.and protect car industry!", 
"\"@DomineekSmith: @realDonaldTrump is the best Republican presidential candidate of all time.\"  Thank you.", 
"Word is that little Morty Zuckerman’s @NYDailyNews loses more than $50 million per year---can that be possible?", 
"\"@Chevy_Mama: @realDonaldTrump I'm obsessed with @celebrityapprenticeNBC. Honestly,  Mr Trump, you are very inspiring\"", 
"President Obama said \"ISIL continues to shrink\" in an interview just hours before the horrible attack in Paris. He is just so bad! CHANGE.", 
".@HillaryClinton loves to lie. America has had enough of the CLINTON'S! It is time to #DrainTheSwamp! Debates https://t.co/3Mz4T7qTTR", 
"\"@jerrimoore: @realDonaldTrump awesome. A treat to get to see the brilliant Joan Rivers once more #icon\"", 
"\"@shoegoddesss: @realDonaldTrump Will definitely vote for you. Breath of fresh air. America needs you!\"", 
"Ted is the ultimate hypocrite. Says one thing for money, does another for votes. \nhttps://t.co/hxdfy0mjVw", 
"\"@Lisa_Milicaj: Truth be told, I  never heard of The National Review until they \"tried\" to declare war on you. No worries, you got my vote!\"", 
"THANK YOU Daytona Beach, Florida!\n#MakeAmericaGreatAgain https://t.co/IAcLfXe463", 
"People rarely say that many conservatives didn’t vote for Mitt Romney. If I can get them to vote for me, we win in a landslide.", 
"Trump National Golf Club, Washington, D.C. is on 600 beautiful acres fronting the Potomac River. A fantastic setting! http://t.co/pYtkbyKwt5", 
"\"@DRUDGE_REPORT: REUTERS 5-DAY ROLLING POLL: TRUMP 34%, CARSON 19.6%, RUBIO 9.7%, CRUZ 7.7%...\" Thank you - a great honor!"
), answer = c(5L, 3L, 2L, 5L, 3L, 1L, 4L, 5L, 2L, 4L, 1L, 1L, 
4L, 4L, 2L, 4L, 4L, 4L, 4L, 4L), target_value = c(NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
NA)), row.names = c(NA, -20L), class = c("tbl_df", "tbl", "data.frame"
                                         ))
trump$lock()
trump2$lock()
```

Calculate the target value (in this case, the AFINN score) by turning one object into a corpus.

```{r}
gold_standard <- trump$turn_gold()
dfm(gold_standard, remove_punct = TRUE) %>% dfm_lookup(afinn) %>% quanteda::convert(to = "data.frame") %>%
    mutate(matching_word_valence = (neg5 * -5) + (neg4 * -4) + (neg3 * -3) + (neg2 * -2) + (neg1 * -1)
           + (zero * 0) + (pos1 * 1) + (pos2 * 2) + (pos3 * 3) + (pos4 * 4) + (pos5 * 5),
           base = ntoken(gold_standard, remove_punct = TRUE), afinn_score = matching_word_valence / base) %>%
    pull(afinn_score) -> target_value
```

Summarize all oolong objects with the target value.

```{r}
res <- summarize_oolong(trump, trump2, target_value = target_value)
```

Read the results. The diagnostic plot consists of 4 subplots. It is a good idea to read Bland & Altman (1986) on the difference between correlation and agreement.

* Subplot (top left): Raw correlation between human judgement and target value. One should want to have a good correlation between the two.
* Subplot (top right): Bland-Altman plot. One should want to have no correlation. Also, the dots should be randomly scattering around the mean value. If it is so, the two measurements (human judgement and target value) are in good agreement.
* Subplot (bottom left): Raw correlation between target value and content length. One should want to have no correlation, as an indication of good reliability against the influence of content length. (See Chan et al. 2020)
* Subplot (bottom right): Cook's distance of all data point. One should want to have no dot (or at least very few dots) above the threshold. It is an indication of how the raw correlation between human judgement and target value can or cannot be influenced by extreme values in your data. 

The textual output contains the Krippendorff's alpha of the codings by your raters. In order to claim validity of your target value, you must first establish the reliability of your gold standard. Song et al. [Forthcoming] suggest Krippendorff's Alpha > 0.7 as an acceptable cut-off.

```{r}
res
```

```{r diagnosis}
plot(res)
```

## Backward compatibility

Historically, oolong test objects could only be generated with only one function: `create_oolong`. It is no longer the case and no longer recommended anymore. It is still retained for backward compatibility purposes. If you still need to use  `create_oolong()`, the most important parameters are `input_model` and `input_corpus`. Setting each of them to `NULL` generates different tests.

| input\_model | input\_corpus | output                                                                                                                                      |
|--------------|:-------------:|---------------------------------------------------------------------------------------------------------------------------------------------|
| Not NULL     | NULL          | oolong test for validating a topic model with [word intrusion test](#word-intrusion-test)                                                   |
| Not NULL     | Not NULL      | oolong test for validating a topic model with [word intrusion test](#word-intrusion-test) and [topic intrusion test](#topic-intrusion-test) |
| NULL         | Not NULL      | oolong test for [creating gold standard](#creating-gold-standard)                                                                           |
| NULL         | NULL          | error                                                                                                                                       |


## References

1. Chang, J., Gerrish, S., Wang, C., Boyd-Graber, J. L., & Blei, D. M. (2009). Reading tea leaves: How humans interpret topic models. In Advances in neural information processing systems (pp. 288-296). [link](https://papers.nips.cc/paper/3700-reading-tea-leaves-how-humans-interpret-topic-models)
2. Ying, L., Montgomery, J. M., & Stewart, B. M. (2021). Inferring concepts from topics: Towards procedures for validating topics as measures. Political Analysis. [link](https://doi.org/10.1017/pan.2021.33)
3. Song et al. (2020) In validations we trust? The impact of imperfect human annotations as a gold standard on the quality of validation of automated content analysis. Political Communication. [link](https://doi.org/10.1080/10584609.2020.1723752)
4. Bland, J. M., & Altman, D. (1986). Statistical methods for assessing agreement between two methods of clinical measurement. The lancet, 327(8476), 307-310.
5. Chan et al. (2020) Four best practices for measuring news sentiment using ‘off-the-shelf’ dictionaries: a large-scale p-hacking experiment. Computational Communication Research. [link](https://osf.io/preprints/socarxiv/np5wa/)
6. Nielsen, F. Å. (2011). A new ANEW: Evaluation of a word list for sentiment analysis in microblogs. arXiv preprint arXiv:1103.2903. [link](https://arxiv.org/abs/1103.2903)

---
---
title: "Deploy"
output: rmarkdown::html_vignette
author:
  - Chung-hong Chan ^[University of Mannheim]
vignette: >
  %\VignetteIndexEntry{Deploy}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  )
set.seed(46709394)
```

In oolong 0.3.22, functions for deploying oolong tests were added (`export_oolong`, `revert_oolong` etc.). These functions make it possible for the coders to conduct validation tests online using their browser, rather than having to install R on their computer.

The basic workflow is simple: 1) create the oolong test object as usual; 2) deploy the test online and obtain the URL to the test; 3) ask your coders to conduct the test online and send back the data file; 4) revert back from the data file to an oolong object.

# Create an oolong test

Please note that one cannot deploy oolong test objects with *both* word and topic intrusion tests, i.e. those created using `witi()` online. If you need to do both tests, you need to deploy them as two separate instances: one created using `wi()` and another created using `ti()`.

In this guide, we assume you want to deploy a word set intrusion test online.

```{r}
library(oolong)
wsi_test <- wsi(abstracts_keyatm)
wsi_test
```

# Deploy the test online

First, you need to export the oolong test object as a stand alone Shiny app. This stand alone Shiny app will be in a directory.

```{r}
export_oolong(wsi_test, dir = "./wsi_test", use_full_path = FALSE)
```

The directory has only two files

```{r}
fs::dir_tree("./wsi_test")
```

This structure is called ["Single-file Shiny app."](https://shiny.rstudio.com/articles/app-formats.html) Experienced Shiny users might have their preferred method of deploying this app to whatever Shiny server they can master.

For less experienced users, the simplest way to deploy this app online is to use [shinyapps.io](https://www.shinyapps.io/) (free tier available with 25 hours of computational time per month). Please register for an account at shinyapps.io and configure rsconnect. Please refer to [this guide](https://shiny.rstudio.com/articles/shinyapps.html) for more information. Please remember to configure the tokens.

```r
## replace <ACCOUNT>, <TOKEN>, <SECRET> with the information from your profile on Shinyapps.io: click Your name -> Tokens
rsconnect::setAccountInfo(name="<ACCOUNT>", token="<TOKEN>", secret="<SECRET>")
```

For RStudio users, the simplest way to deploy the app to shinyapps.io is to first launch the app.

```{r, eval = FALSE}
library(shiny)
runApp("./wsi_test")
```

And then click the **Publish** button at the right corner of the launched window.

You will be asked for the title of the app, just give it a name, e.g. *wsi_test*. You probably can keep other default settings and push the **Publish** button to initialize the deployment process.

<img src="figures/deploying_shinyappsio.png" align="center" height="400" />

If there is no hiccup, you will get a URL to your deployed oolong test. Something like: *https://yourname.shinyapps.io/wsi_test/*

# Conduct the test

You can give the URL to your coders and they conduct the test with their browser online. The only difference of the deployed version is that, there will be a userid prompt and download button after the coding.

<img src="figures/oolong_download.png" align="center"/>

You should instruct your coders to download the data file after coding and return it to you. ^[Future versions might provide permanent storage]

# Revert

You can then obtain a locked oolong object from the original oolong and the downloaded data file. `revert_oolong` will do verifications with the original oolong object to make sure no error and no cheating.

```{r, include = FALSE}
wsi_test <- readRDS(system.file("extdata", "wsi_test.RDS", package = "oolong"))
```

```r
revert_oolong(wsi_test, "oolong_2021-05-22 20 51 26 Hadley Wickham.RDS")
```

```{r, echo = FALSE}
revert_oolong(wsi_test, system.file("extdata", "hadley.RDS", package = "oolong"))
```

```{r, include = FALSE}
unlink("./wsi_test", recursive = TRUE)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong_shiny.R
\name{export_oolong}
\alias{export_oolong}
\title{Export a deployable Shiny app from an oolong object into a directory}
\usage{
export_oolong(
  oolong,
  dir = base::tempdir(),
  verbose = TRUE,
  use_full_path = TRUE
)
}
\arguments{
\item{oolong}{an oolong object to be exported. Please note that the "witi" type, i.e. oolong object with both word and topic intrusion tests, cannot be exported. Also the object must not be locked and ever coded.}

\item{dir}{character string, the directory to be exported. Default to a temporary directory}

\item{verbose}{logical, whether to display information after exporting}

\item{use_full_path}{logical, whether to expand dir into full path}
}
\value{
directory exported, invisible
}
\description{
This function exports your oolong test into a launched Shiny app that is ideal for online deployment. Deploying the Shiny app online allows coders to conduct the test online with their browser, rather than having to install R on their own computer. In contrast to the testing interfaces launched with methods such as \code{$do_word_intrusion_test()}, the deployable version provides data download after the coder finished coding. Downloaded data can then revert back to a locked oolong object using \code{revert_oolong}. Further version might provide solutions to permanent storage. The deployable Shiny app will be in a directory. The Shiny app is both launchable with shiny::runApp() and deployable with rsconnect::deployApp(). Please refer to \code{vignette("deploy", package = "oolong")} for more details.
}
\examples{
# Please try this example in interactive R sessions only.
if (interactive()) {
   data(abstracts_stm)
   x <- wi(abstracts_stm)
   export_oolong(x)
}
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong.R, R/oolong_ui.R
\name{create_oolong}
\alias{create_oolong}
\alias{wi}
\alias{witi}
\alias{ti}
\alias{wsi}
\alias{gs}
\title{Generate an oolong test}
\usage{
create_oolong(
  input_model = NULL,
  input_corpus = NULL,
  n_top_terms = 5,
  bottom_terms_percentile = 0.6,
  exact_n = NULL,
  frac = 0.01,
  n_top_topics = 3,
  n_topiclabel_words = 8,
  use_frex_words = FALSE,
  difficulty = 1,
  input_dfm = NULL,
  construct = "positive",
  btm_dataframe = NULL,
  n_correct_ws = 3,
  wsi_n_top_terms = 20,
  userid = NA,
  type = "witi"
)

wi(
  input_model = NULL,
  userid = NA,
  n_top_terms = 5,
  bottom_terms_percentile = 0.6,
  difficulty = 1,
  use_frex_words = FALSE
)

witi(
  input_model = NULL,
  input_corpus = NULL,
  userid = NA,
  n_top_terms = 5,
  bottom_terms_percentile = 0.6,
  exact_n = NULL,
  frac = 0.01,
  n_top_topics = 3,
  n_topiclabel_words = 8,
  use_frex_words = FALSE,
  difficulty = 1,
  input_dfm = NULL,
  btm_dataframe = NULL
)

ti(
  input_model = NULL,
  input_corpus = NULL,
  userid = NA,
  exact_n = NULL,
  frac = 0.01,
  n_top_topics = 3,
  n_topiclabel_words = 8,
  use_frex_words = FALSE,
  difficulty = 1,
  input_dfm = NULL,
  btm_dataframe = NULL
)

wsi(
  input_model = NULL,
  userid = NA,
  n_topiclabel_words = 4,
  n_correct_ws = 3,
  wsi_n_top_terms = 20,
  difficulty = 1,
  use_frex_words = FALSE
)

gs(
  input_corpus = NULL,
  userid = NA,
  construct = "positive",
  exact_n = NULL,
  frac = 0.01
)
}
\arguments{
\item{input_model}{(wi, ti, witi, wsi) a STM, WarpLDA, topicmodels, KeyATM, seededlda, textmodel_nb, or BTM object; if it is NULL, create_oolong assumes that you want to create gold standard.}

\item{input_corpus}{(wi, ti, witi, wsi, gs) if input_model is not null, it should be the corpus (character vector or quanteda::corpus object) to generate the model object. If input_model and input_corpus are not NULL, topic intrusion test cases are generated. If input_model is a BTM object, this argument is ignored. If input_model is null, it generates gold standard test cases.}

\item{n_top_terms}{(wi, witi) integer, number of top topic words to be included in the candidates of word intrusion test.}

\item{bottom_terms_percentile}{(wi, witi) double, a term is considered to be an word intruder when its theta less than the percentile of this theta, must be within the range of 0 to 1}

\item{exact_n}{(ti, witi, gs) integer, number of topic intrusion test cases to generate, ignore if frac is not NULL}

\item{frac}{(ti, witi, gs) double, fraction of test cases to be generated from the corpus}

\item{n_top_topics}{(wi, witi) integer, number of most relevant topics to be shown alongside the intruder topic}

\item{n_topiclabel_words}{(witi, ti, wsi) integer, number of topic words to be shown as the topic ("ti" and "witi") / word set ("wsi") label}

\item{use_frex_words}{(wi, witi, ti, wsi) logical, for a STM object, use FREX words if TRUE, use PROB words if FALSE}

\item{difficulty}{(wi, witi, ti, wsi) double, adjust the difficulty of the test. Higher value indicates higher difficulty and must be within the range of 0 to 1, no effect for STM if use_frex_words is FALSE. Ignore for topicmodels objects}

\item{input_dfm}{(wi, witi, ti, wsi) a dfm object used for training the input_model, if input_model is a WarpLDA object}

\item{construct}{(gs) string, an adjective to describe the construct you want your coders to code the the gold standard test cases}

\item{btm_dataframe}{(witi, ti) dataframe used for training the input_model, if input_model is a BTM object}

\item{n_correct_ws}{(wsi) number of word sets to be shown alongside the intruder word set}

\item{wsi_n_top_terms}{(wsi) number of top topic words from each topic to be randomized selected as the word set label}

\item{userid}{a character string to denote the name of the coder. Default to NA (no userid); not recommended}

\item{type}{(create_oolong) a character string to denote what you want to create. "wi": word intrusion test; "ti": topic intrusion test; "witi": both word intrusion test and topic intrusion test; "gs": gold standard generation}
}
\value{
an oolong test object.
}
\description{
\code{create_oolong} generates an oolong test object that can either be used for validating a topic model or for creating ground truth (gold standard) of a text corpus. \code{wi} (word intrusion test), \code{ti} (topic intrusion test), \code{witi} (word and topic intrusion tests), \code{wsi} (word set intrusion test) and \code{gs} are handy wrappers to \code{create_oolong}. It is recommended to use these wrappers instead of \code{create_oolong}.
}
\section{Usage}{


Use \code{wi}, \code{ti}, \code{witi}, \code{wsi} or \code{gs} to generate an oolong test of your choice. It is recommended to supply also \code{userid} (current coder).
The names of the tests (word intrusion test and topic intrusion test) follow Chang et al (2009). In Ying et al. (forthcoming), topic intrusion test is named "T8WSI" (Top 8 Word Set Intrusion). Word set intrusion test in this package is actually the "R4WSI" (Random 4 Word Set Intrusion) in Lu et al (forthcoming). The default settings of \code{wi}, \code{witi}, and \code{ti} follow Chang et al (2009), e.g. \code{n_top_terms} = 5; instead of \code{n_top_terms} = 4 as in Lu et al (forthcoming). The default setting of \code{wsi} follows Ying et al. (forthcoming), e.g. \code{n_topiclabel_words} = 4.
As suggested by Song et al. (2020), 1% of the articles from \code{input_corpus} are randomly selected as the test cases of both \code{ti} and \code{gs}, i.e. \code{frac} = 0.01. However, it is generally believed that this proportion is dependent of the size of \code{input_corpus}, e.g. it does not make sense to draw 1% of the articles from only 100 articles. Use \code{exact_n} in these cases.
}

\section{About create_oolong}{


Because \code{create_oolong} is not intuitive to use, it is no longer recommended to use \code{create_oolong} to generate oolong test. \code{create_oolong} is retained only for backward compatibility purposes. This function generates an oolong test object based on \code{input_model} and \code{input_corpus}. If \code{input_model} is not NULL, it generates oolong test for a topic model (tm). If \code{input_model} is NULL but input_corpus is not NULL, it generates oolong test for generating gold standard (gs).
}

\section{Methods}{

An oolong object, depends on its purpose, has the following methods:
\describe{
  \item{\code{$do_word_intrusion_test()}}{(tm) launch the shiny-based word intrusion test. The coder should find out the intruder word that is not related to other words.}
  \item{\code{$do_topic_intrusion_test()}}{(tm) launch the shiny-based topic intrusion test. The coder should find out the intruder topic that is least likely to be the topic of the document.}
  \item{\code{$do_word_set_intrusion_test()}}{(tm) launch the shiny-based word set intrusion test. The coder should find out the intruder word set that is not related to other word sets.}
  \item{\code{$do_gold_standard_test()}}{(gs) launch the shiny-based test for generating gold standard. The coder should determine the level of the predetermined constructs with a 5-point Likert scale.}
  \item{\code{$lock(force = FALSE)}}{(gs/tm) lock the object so that it cannot be changed anymore. It enables \code{\link{summarize_oolong}} and the following method.}
  \item{\code{$turn_gold()}}{(gs) convert the oolong object into a quanteda compatible corpus.}
}
For more details, please see the overview vignette: \code{vignette("overview", package = "oolong")}
}

\examples{
## Creation of oolong test with only word intrusion test
data(abstracts_keyatm)
data(abstracts)
oolong_test <- wi(input_model = abstracts_keyatm, userid = "Hadley")
## Creation of oolong test with both word intrusion test and topic intrusion test
oolong_test <- witi(input_model = abstracts_keyatm, input_corpus = abstracts$text, userid = "Julia")
## Creation of oolong test with topic intrusion test
oolong_test <- ti(input_model = abstracts_keyatm, input_corpus = abstracts$text, userid = "Jenny")
## Creation of oolong test with word set intrusion test
oolong_test <- wsi(input_model = abstracts_keyatm, userid = "Garrett")
## Creation of gold standard
oolong_test <- gs(input_corpus = trump2k, userid = "Yihui")
## Using create_oolong(); not recommended
oolong_test <- create_oolong(input_model = abstracts_keyatm,
input_corpus = abstracts$text, userid = "JJ")
oolong_test <- create_oolong(input_model = abstracts_keyatm,
input_corpus = abstracts$text, userid = "Mara", type = "ti")
oolong_test <- create_oolong(input_corpus = abstracts$text, userid = "Winston", type = "gs")
}
\references{
Chang, J., Gerrish, S., Wang, C., Boyd-Graber, J. L., & Blei, D. M. (2009). Reading tea leaves: How humans interpret topic models. In Advances in neural information processing systems (pp. 288-296).

  Song et al. (2020) In validations we trust? The impact of imperfect human annotations as a gold standard on the quality of validation of automated content analysis. Political Communication.

  Ying, L., Montgomery, J. M., & Stewart, B. M. (Forthcoming). Inferring concepts from topics: Towards procedures for validating topics as measures. Political Analysis
}
\author{
Chung-hong Chan, Marius Sältzer
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong_gs.R
\name{print.oolong_gold_standard}
\alias{print.oolong_gold_standard}
\title{Print oolong gold standard object}
\usage{
\method{print}{oolong_gold_standard}(x, ...)
}
\arguments{
\item{x}{an oolong gold standard object}

\item{...}{other parameters}
}
\value{
None, a summary of the quanteda::corpus and what you should do are displayed
}
\description{
This function prints a summary of the oolong gold standard object. An oolong gold standard object is a result of $turn_gold() method. It is a quanteda::corpus compatible object.
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong_summary.R
\name{summarize_oolong}
\alias{summarize_oolong}
\alias{summarise_oolong}
\title{Summarize oolong objects}
\usage{
summarize_oolong(..., target_value = NULL, n_iter = 1500)

summarise_oolong(..., target_value = NULL, n_iter = 1500)
}
\arguments{
\item{...}{(tm/gs) one or more oolong objects to be summarized}

\item{target_value}{(gs) a vector of numeric values, the value you want to validate against the human-coded gold standard. One example of this target value is sentiment score extracted automatically from text}

\item{n_iter}{(ti) number of iterations to calculate the median test}
}
\value{
An oolong summary.
Depends on purpose, an oolong summary object has the following values:
\describe{
  \item{\code{$type}}{(gs/tm) type of analysis, either 'gs' or 'tm'}
  \item{\code{$kripp_alpha}; \code{$kripp_alpha_wsi}}{(wi, wsi) Krippendorff's Alpha, if more than one oolong object is analyzed.}
  \item{\code{$rater_precision}; \code{$rater_precision_wsi}}{(wi, wsi) Model precision}
  \item{\code{$res$rater_precision_p_value}}{(wi) Model precision's p-value calculated by one-sample binomial test and Fisher's Omnibus method.}
  \item{\code{$k_precision}; \code{$k_precision_wsi}}{(wi, wsi) precision for each topic}
  \item{\code{$tlo}}{(ti) vector of topic log odds}
  \item{\code{$tlo_pvalue}}{(ti) Median topic log odds's p-value calculated by permutation test.}
  \item{\code{$cor}}{(gs) Pearson's correlation between average answer and target value}
  \item{\code{$cor_length}}{(gs) Pearson's correlation between content length and target value}
  \item{\code{$diag_plot}}{(gs) diagnostic plot.}
}
A useful summary of an object can be obtained either by \code{\link{print.oolong_summary}} or \code{\link{plot.oolong_summary}}. For details, please see the overview vignette: \code{vignette("overview", package = "oolong")}
}
\description{
This function summarizes one or more oolong objects. All oolong objects must be locked.
}
\examples{
# Please try this example in interactive R sessions only.
if (interactive()) {
   data(abstracts_stm)
   oolong_test1 <- create_oolong(abstracts_stm)
   oolong_test2 <- clone_oolong(oolong_test1)
   oolong_test1$do_word_intrusion_test()
   oolong_test2$do_word_intrusion_test()
   oolong_test1$lock()
   oolong_test2$lock()
   summarize_oolong(oolong_test1, oolong_test2)
}
}
\references{
Chang, J., Gerrish, S., Wang, C., Boyd-Graber, J. L., & Blei, D. M. (2009). Reading tea leaves: How humans interpret topic models. In Advances in neural information processing systems (pp. 288-296).

  Song et al. (2020) In validations we trust? The impact of imperfect human annotations as a gold standard on the quality of validation of automated content analysis. Political Communication.

  Ying, L., Montgomery, J. M., & Stewart, B. M. (Forthcoming). Inferring concepts from topics: Towards procedures for validating topics as measures. Political Analysis.
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong_update.R
\name{update_oolong}
\alias{update_oolong}
\title{Update an oolong object to the latest version}
\usage{
update_oolong(oolong, verbose = TRUE)
}
\arguments{
\item{oolong}{an oolong object to be updated}

\item{verbose, }{logical, display messages}
}
\value{
an updated oolong object
}
\description{
This function update an old oolong object to the latest version.
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong.R
\name{clone_oolong}
\alias{clone_oolong}
\title{Clone an oolong object}
\usage{
clone_oolong(oolong, userid = NA)
}
\arguments{
\item{oolong}{an oolong object.}

\item{userid}{a character string to denote the name of the coder}
}
\value{
an oolong object
}
\description{
Clone a new oolong object. The oolong must not be locked and ever coded.
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong_update.R
\name{check_oolong}
\alias{check_oolong}
\title{Check whether the oolong needs to be updated}
\usage{
check_oolong(oolong, verbose = FALSE)
}
\arguments{
\item{oolong}{an oolong object to be checked}

\item{verbose, }{logical, display messages}
}
\value{
Nothing
}
\description{
This function raises an error when the input oolong object needs to be updated. Oolong objects generated with an old version of oolong need to be updated to use the functionalities from the recent versions of oolong.
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong.R
\name{revert_oolong}
\alias{revert_oolong}
\title{Obtain a locked oolong from a downloaded data file}
\usage{
revert_oolong(oolong, rds_file)
}
\arguments{
\item{oolong}{an oolong object used for deployment}

\item{rds_file}{path to the downloaded RDS file}
}
\value{
a locked oolong object based on the data in the downloaded RDS file
}
\description{
To generate a locked oolong object with the original oolong object and the RDS file. The RDS file should have been downloaded from a deployed Shiny app.
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong_data_misc.R
\docType{data}
\name{abstracts_keyatm}
\alias{abstracts_keyatm}
\alias{abstracts_warplda}
\alias{abstracts_btm}
\title{Topic models trained with the abstracts dataset.}
\format{
An object of class \code{keyATM_output} (inherits from \code{base}, \code{list}) of length 18.

An object of class \code{WarpLDA} (inherits from \code{LDA}, \code{TopicModel}, \code{mlapiDecomposition}, \code{mlapiTransformation}, \code{mlapiBase}, \code{R6}) of length 9.

An object of class \code{BTM} of length 9.
}
\usage{
abstracts_keyatm

abstracts_warplda

abstracts_btm
}
\description{
These are topic models trained with different topic model packages.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong_data_misc.R
\docType{data}
\name{abstracts}
\alias{abstracts}
\alias{abstracts_dfm}
\alias{abstracts_dictionary}
\title{Abstracts of communication journals dataset}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 2500 rows and 1 columns.

An object of class \code{dfm} with 2500 rows and 3998 columns.

An object of class \code{list} of length 10.
}
\usage{
abstracts

abstracts_dfm

abstracts_dictionary
}
\description{
This is a random sample of all abstracts of papers published in high-impact communication journals from 2000 to 2017. abstracts_dictionary is a list of terms that can be used for semisupervised techniques such as keyATM.
}
\references{
Chan, C-h, & Grill, C. (2020). [The Highs in Communication Research: Research Topics With High Supply, High Popularity, and High Prestige in High-Impact Journals.](https://doi.org/10.1177/0093650220944790) Communication Research.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong_data_misc.R
\docType{data}
\name{trump2k}
\alias{trump2k}
\title{Trump's tweets dataset}
\format{
An object of class \code{character} of length 2000.
}
\usage{
trump2k
}
\description{
This is a random sample of 2000 tweets from @realdonaldtrump account before his assumption of duty as the president of the United States.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong_shiny.R
\name{deploy_oolong}
\alias{deploy_oolong}
\title{Deploy an oolong test}
\usage{
deploy_oolong(oolong)
}
\arguments{
\item{oolong}{an oolong object to be deployed. Please note that the "witi" type, i.e. oolong object with both word and topic intrusion tests, cannot be deployed. Also the object must not be locked and ever coded.}
}
\value{
Nothing, it launches a deployable version of the coding interface
}
\description{
In most of the time, you should not use this function. You should write the deployable version of your app into a directory using \code{export_oolong} instead. Please refer to \code{vignette("deploy", package = "oolong")} for more details.
}
\examples{
# Please try this example in interactive R sessions only.
if (interactive()) {
   data(abstracts_stm)
   x <- wi(abstracts_stm)
   deploy_oolong(x)
}
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong.R
\docType{package}
\name{oolong}
\alias{oolong}
\title{Oolong: create and administrate validation tests for typical automated content analysis tools}
\description{
The oolong package makes it very easy to create, administrate and evaluate typical automated content analysis tools by providing a framework for creating human-in-the-loop validation tests. For topic models, oolong can generate tests such as word intrusion test, topic intrusion test (Chang et al.) and word set intrusion test (Ying et al.) For dictionary-based methods, oolong can generate standardized interface for making gold standard ('Ground truth') data.
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong_summary.R
\name{print.oolong_summary}
\alias{print.oolong_summary}
\alias{plot.oolong_summary}
\title{Print and plot oolong summary}
\usage{
\method{print}{oolong_summary}(x, ...)

\method{plot}{oolong_summary}(x, ...)
}
\arguments{
\item{x}{an oolong_summary}

\item{...}{other parameters}
}
\value{
None
}
\description{
These functions print or plot a useful summary of the results from \code{\link{summarize_oolong}}. For details, please see the overview vignette: \code{vignette("overview", package = "oolong")}
}
\section{Summary}{


Print function displays the following information:
\describe{
  \item{Mean model precision}{(wi, wsi) Higher value indicates better topic interpretability}
  \item{Quantiles of model precision}{(wi) Higher value indicates better topic interpretability}
  \item{P-value of the model precision}{(wi) Model precision's p-value calculated by one-sample binomial test and Fisher's Omnibus method.}
  \item{Krippendorff's alpha}{(wi, wsi, gs) Krippendorff's Alpha, if more than one oolong object is analyzed.}
  \item{K Precision}{(wi, wsi) Model precision for each topic.}
  \item{Mean TLO}{(ti) Mean topic log odds, higher value indicates better interpretability}
  \item{Median TLO}{(ti) Median topic log odds, higher value indicates better interpretability}
  \item{Quantiles of TLO}{(ti) Quantiles of topic log odds}
  \item{P-Value of the median TLO}{(ti) Median topic log odds's p-value calculated by permutation test.}
  \item{Correlation (average answer)}{(gs) Pearson's correlation between average answer and target value}
  \item{Corrlation (content length)}{(gs) Pearson's correlation between content length and target value}
}
}

\section{Diagnostic plot}{


Plot function displays a diagnostic plot with the following subplots (gs only).
\describe{
  \item{Top left}{Correlation between answer from coders and target value to check for correlation between two values. Both axes are minmax transformed.}
  \item{Top right}{Bland-altman plot of answer from coders and target value to check for agreement between two values.}
  \item{Bottom left}{Correlation between target value and content length to check for the influence of content length.}
  \item{Bottom right}{Cook's distance to check for influential observations.}
}
}

\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong_data_misc.R
\docType{data}
\name{newsgroup_nb}
\alias{newsgroup_nb}
\title{Naive Bayes model trained on 20 newsgroups data}
\format{
An object of class \code{textmodel_nb} (inherits from \code{textmodel}, \code{list}) of length 7.
}
\usage{
newsgroup_nb
}
\description{
This is a Naive Bayes model (of the class 'textmodel_nb') trained on 20 newsgroups data.
}
\references{
Lang, K. (1995). Newsweeder: Learning to filter netnews. In Machine Learning Proceedings 1995 (pp. 331-339). Morgan Kaufmann.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oolong_data_misc.R
\docType{data}
\name{afinn}
\alias{afinn}
\title{AFINN dictionary}
\format{
An object of class \code{dictionary2} of length 11.
}
\usage{
afinn
}
\description{
This is the AFINN sentiment dictionary in quanteda::dictionary format.
}
\references{
Nielsen, F. Å. (2011). A new ANEW: Evaluation of a word list for sentiment analysis in microblogs. arXiv preprint arXiv:1103.2903.
}
\keyword{datasets}
