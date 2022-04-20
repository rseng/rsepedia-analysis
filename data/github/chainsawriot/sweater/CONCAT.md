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

<!-- README.md is generated from README.Rmd. Please edit that file -->

# sweater <img src="man/figures/sweater_logo.svg" align="right" height="200" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/chainsawriot/sweater/workflows/R-CMD-check/badge.svg)](https://github.com/chainsawriot/sweater/actions)
[![Codecov test
coverage](https://codecov.io/gh/chainsawriot/sweater/branch/master/graph/badge.svg)](https://app.codecov.io/gh/chainsawriot/sweater?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/sweater)](https://CRAN.R-project.org/package=sweater)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04036/status.svg)](https://doi.org/10.21105/joss.04036)
<!-- badges: end -->

The goal of sweater (**S**peedy **W**ord **E**mbedding **A**ssociation
**T**est & **E**xtras using **R**) is to test for associations among
words in word embedding spaces. The methods provided by this package can
also be used to test for unwanted associations, or biases.

The package provides functions that are speedy. They are either
implemented in C++, or are speedy but accurate approximation of the
original implementation proposed by Caliskan et al (2017). See the
benchmark
[here](https://github.com/chainsawriot/sweater/blob/master/paper/benchmark.md).

This package provides extra methods such as Relative Norm Distance,
Embedding Coherence Test, SemAxis and Relative Negative Sentiment Bias.

If your goal is to reproduce the analysis in Caliskan et al (2017),
please consider using the [original Java
program](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/DX4VWP&version=2.0)
or the R package [cbn](https://github.com/conjugateprior/cbn) by Lowe.
To reproduce the analysis in Garg et al (2018), please consider using
the [original Python
program](https://github.com/nikhgarg/EmbeddingDynamicStereotypes). To
reproduce the analysis in Mazini et al (2019), please consider using the
[original Python
program](https://github.com/TManzini/DebiasMulticlassWordEmbedding/).

Please cite this software as:

Chan, C., (2022). sweater: Speedy Word Embedding Association Test and
Extras Using R. Journal of Open Source Software, 7(72), 4036,
<https://doi.org/10.21105/joss.04036>

For a BibTeX entry, use the output from `citation(package = "sweater")`.

## Installation

Recommended: install the latest development version

``` r
remotes::install_github("chainsawriot/sweater")
```

or the “stable” release

``` r
install.packages("sweater")
```

## Notation of a query

All tests in this package use the concept of queries (see Badilla et
al., 2020) to study associations in the input word embeddings `w`. This
package uses the “STAB” notation from Brunet et al (2019). \[1\]

All tests depend on two types of words. The first type, namely,
`S_words` and `T_words`, is *target words* (or *neutral words* in Garg
et al). In the case of studying biases, these are words that **should**
have no bias. For instance, the words such as “nurse” and “professor”
can be used as target words to study the gender bias in word embeddings.
One can also separate these words into two sets, `S_words` and
`T_words`, to group words by their perceived bias. For example, Caliskan
et al. (2017) grouped target words into two groups: mathematics (“math”,
“algebra”, “geometry”, “calculus”, “equations”, “computation”,
“numbers”, “addition”) and arts (“poetry”, “art”, “dance”,
“literature”, “novel”, “symphony”, “drama”, “sculpture”). Please note
that also `T_words` is not always required.

The second type, namely `A_words` and `B_words`, is *attribute words*
(or *group words* in Garg et al). These are words with known properties
in relation to the bias that one is studying. For example, Caliskan et
al. (2017) used gender-related words such as “male”, “man”, “boy”,
“brother”, “he”, “him”, “his”, “son” to study gender bias. These words
qualify as attribute words because we know they are related to a certain
gender.

It is recommended to use the function `query()` to make a query and
`calculate_es()` to calculate the effect size.

## Available methods

| Target words       | Attribute words    | Method                                                      | `method` argument | Suggested by `query`? | legacy functions \[2\]                                |
| ------------------ | ------------------ | ----------------------------------------------------------- | ----------------- | --------------------- | ----------------------------------------------------- |
| S\_words           | A\_words           | Mean Average Cosine Similarity (Mazini et al. 2019)         | “mac”             | yes                   | mac(), mac\_es()                                      |
| S\_words           | A\_words, B\_words | Relative Norm Distance (Garg et al. 2018)                   | “rnd”             | yes                   | rnd(), rnd\_es()                                      |
| S\_words           | A\_words, B\_words | Relative Negative Sentiment Bias (Sweeney & Najafian. 2019) | “rnsb”            | no                    | rnsb(), rnsb\_es()                                    |
| S\_words           | A\_words, B\_words | Embedding Coherence Test (Dev & Phillips. 2019)             | “ect”             | no                    | ect(), ect\_es(), plot\_ect()                         |
| S\_words           | A\_words, B\_words | SemAxis (An et al. 2018)                                    | “semaxis”         | no                    | semaxis()                                             |
| S\_words           | A\_words, B\_words | Normalized Association Score (Caliskan et al. 2017)         | “nas”             | no                    | nas()                                                 |
| S\_words, T\_words | A\_words, B\_words | Word Embedding Association Test (Caliskan et al. 2017)      | “weat”            | yes                   | weat(), weat\_es(), weat\_resampling(), weat\_exact() |
| S\_words, T\_words | A\_words, B\_words | Word Embeddings Fairness Evaluation (Badilla et al. 2020)   | To be implemented |                       |                                                       |

## Example: Mean Average Cosine Similarity

The simplest form of bias detection is Mean Average Cosine Similarity
(Mazini et al. 2019). The same method is used also in Kroon et
al. (2020). `googlenews` is a subset of [the pretrained word2vec word
embeddings provided by
Google](https://code.google.com/archive/p/word2vec/).

By default, the `query()` function guesses the method you want to use
based on the combination of target words and attribute words provided
(see the “Suggested?” column in the above table). You can also make this
explicit by specifying the `method` argument. Printing the returned
object shows the effect size (if available) as well as the functions
that can further process the object: `calculate_es` and `plot`. Please
read the help file of `calculate_es` (`?calculate_es`) on what is the
meaning of the effect size for a specific test.

``` r
require(sweater)
```

``` r
S1 <- c("janitor", "statistician", "midwife", "bailiff", "auctioneer", 
"photographer", "geologist", "shoemaker", "athlete", "cashier", 
"dancer", "housekeeper", "accountant", "physicist", "gardener", 
"dentist", "weaver", "blacksmith", "psychologist", "supervisor", 
"mathematician", "surveyor", "tailor", "designer", "economist", 
"mechanic", "laborer", "postmaster", "broker", "chemist", "librarian", 
"attendant", "clerical", "musician", "porter", "scientist", "carpenter", 
"sailor", "instructor", "sheriff", "pilot", "inspector", "mason", 
"baker", "administrator", "architect", "collector", "operator", 
"surgeon", "driver", "painter", "conductor", "nurse", "cook", 
"engineer", "retired", "sales", "lawyer", "clergy", "physician", 
"farmer", "clerk", "manager", "guard", "artist", "smith", "official", 
"police", "doctor", "professor", "student", "judge", "teacher", 
"author", "secretary", "soldier")

A1 <- c("he", "son", "his", "him", "father", "man", "boy", "himself", 
"male", "brother", "sons", "fathers", "men", "boys", "males", 
"brothers", "uncle", "uncles", "nephew", "nephews")

## The same as:
## mac_neg <- query(googlenews, S_words = S1, A_words = A1, method = "mac")
mac_neg <- query(googlenews, S_words = S1, A_words = A1)
mac_neg
#> 
#> ── sweater object ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> Test type:  mac 
#> Effect size:  0.1375856
#> 
#> ── Functions ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> • `calculate_es()`: Calculate effect size
#> • `plot()`: Plot the bias of each individual word
```

The returned object is an S3 object. Please refer to the help file of
the method for the definition of all slots (in this case: `?mac`). For
example, the magnitude of bias for each word in `S1` is available in the
`P` slot.

``` r
sort(mac_neg$P)
#>         sales      designer     economist       manager      clerical 
#>  -0.002892495   0.039197285   0.046155954   0.047322071   0.048912403 
#>      operator administrator        author    auctioneer        tailor 
#>   0.050275206   0.050319552   0.051470909   0.065440629   0.074771460 
#>     secretary     librarian     scientist  statistician         pilot 
#>   0.077506781   0.079040760   0.082535536   0.088000351   0.088337791 
#>     geologist      official     architect        broker     professor 
#>   0.088567238   0.090706054   0.091598653   0.098761198   0.101847166 
#>      engineer     collector         smith       chemist      surveyor 
#>   0.103448025   0.104596505   0.104956871   0.110798023   0.112098241 
#>     inspector        weaver     physicist       midwife    supervisor 
#>   0.112383017   0.113221694   0.114302092   0.115791724   0.118784135 
#>     physician        artist     conductor        clergy         guard 
#>   0.118990813   0.119571390   0.120602413   0.123313906   0.128804364 
#>    accountant    instructor         judge    postmaster         nurse 
#>   0.131700192   0.133135210   0.135238197   0.138497652   0.143781092 
#>          cook     attendant       sheriff        dancer  photographer 
#>   0.145019382   0.149134946   0.149992633   0.150637430   0.151388282 
#>  psychologist       cashier       surgeon mathematician       retired 
#>   0.151908676   0.153591372   0.158348402   0.158969004   0.165010593 
#>         clerk       student        porter      gardener       dentist 
#>   0.165903226   0.167006052   0.172551327   0.173346664   0.174776368 
#>       teacher       athlete       bailiff       painter        driver 
#>   0.175027901   0.176353551   0.176440157   0.176625091   0.181269327 
#>         baker     shoemaker        lawyer    blacksmith        farmer 
#>   0.183320490   0.183548112   0.189963886   0.198764788   0.199243319 
#>         mason        police   housekeeper        sailor      musician 
#>   0.203577329   0.206264491   0.208280255   0.208689761   0.219184802 
#>       janitor      mechanic        doctor       soldier       laborer 
#>   0.220953800   0.224008333   0.226657160   0.238053858   0.251032714 
#>     carpenter 
#>   0.259775292
```

## Example: Relative Norm Distance

This analysis reproduces the analysis in Garg et al (2018), namely
Figure 1.

``` r
B1 <- c("she", "daughter", "hers", "her", "mother", "woman", "girl", 
"herself", "female", "sister", "daughters", "mothers", "women", 
"girls", "females", "sisters", "aunt", "aunts", "niece", "nieces"
)

garg_f1 <- query(googlenews, S_words = S1, A_words = A1, B_words = B1)
garg_f1
#> 
#> ── sweater object ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> Test type:  rnd 
#> Effect size:  -6.341598
#> 
#> ── Functions ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> • `calculate_es()`: Calculate effect size
#> • `plot()`: Plot the bias of each individual word
```

The object can be plotted by the function `plot` to show the bias of
each word in S. Words such as “nurse”, “midwife” and “librarian” are
more associated with female, as indicated by the positive relative norm
distance.

``` r
plot(garg_f1)
```

<img src="man/figures/README-rndplot-1.png" width="100%" />

The effect size is simply the sum of all relative norm distance values
(Equation 3 in Garg et al. 2018). It is displayed simply by printing the
object. You can also use the function `calculate_es` to obtain the
numeric result.

The more positive effect size indicates that words in `S_words` are more
associated with `B_words`. As the effect size is negative, it indicates
that the concept of occupation is more associated with `A_words`,
i.e. male.

``` r
calculate_es(garg_f1)
#> [1] -6.341598
```

## Example: SemAxis

This analysis attempts to reproduce the analysis in An et al. (2018).

You may obtain the word2vec word vectors trained with Trump supporters
Reddit from [here](https://github.com/ghdi6758/SemAxis). This package
provides a tiny version of the data `small_reddit` for reproducing the
analysis.

``` r
S2 <- c("mexicans", "asians", "whites", "blacks", "latinos")
A2 <- c("respect")
B2 <- c("disrespect")
res <- query(small_reddit, S_words = S2, A_words = A2, B_words = B2, method = "semaxis", l = 1)
plot(res)
```

<img src="man/figures/README-semxaxisplot-1.png" width="100%" />

## Example: Embedding Coherence Test

Embedding Coherence Test (Dev & Phillips, 2019) is similar to SemAxis.
The only significant different is that no “SemAxis” is calculated (the
difference between the average word vectors of `A_words` and `B_words`).
Instead, it calculates two separate axes for `A_words` and `B_words`.
Then it calculates the proximity of each word in `S_words` with the two
axes. It is like doing two separate `mac`, but `ect` averages the word
vectors of `A_words` / `B_words` first.

It is important to note that `P` is a 2-D matrix. Hence, the plot is
2-dimensional. Words above the equality line are more associated with
`B_words` and vice versa.

``` r
res <- query(googlenews, S_words = S1, A_words = A1, B_words = B1, method = "ect")
res$P
#>           janitor statistician   midwife   bailiff auctioneer photographer
#> A_words 0.3352883   0.13495237 0.1791162 0.2698131 0.10123085    0.2305419
#> B_words 0.2598501   0.08300127 0.3851766 0.2331852 0.06957685    0.2077952
#>          geologist shoemaker   athlete   cashier    dancer housekeeper
#> A_words 0.13817054 0.2842002 0.2607956 0.2340296 0.2282981   0.3205498
#> B_words 0.05101061 0.1850456 0.2570477 0.3171645 0.3508183   0.4610773
#>         accountant  physicist  gardener   dentist    weaver blacksmith
#> A_words  0.2029543 0.17446868 0.2657907 0.2672548 0.1767915  0.3080301
#> B_words  0.1789482 0.08362829 0.2873140 0.2623802 0.2475565  0.1603038
#>         psychologist supervisor mathematician   surveyor     tailor   designer
#> A_words    0.2322444  0.1852041     0.2423898 0.17124643 0.11379186 0.06231389
#> B_words    0.2418605  0.1920407     0.1332954 0.09133125 0.07585015 0.14343468
#>          economist  mechanic   laborer postmaster    broker   chemist librarian
#> A_words 0.07450962 0.3435494 0.3904412  0.2128712 0.1525395 0.1696522 0.1237070
#> B_words 0.04008006 0.1882135 0.3011930  0.2223472 0.1112061 0.1440956 0.3147546
#>         attendant   clerical  musician    porter scientist carpenter    sailor
#> A_words 0.2278508 0.07601974 0.3349666 0.2642203 0.1263250 0.4006367 0.3169384
#> B_words 0.2495253 0.15137979 0.2735083 0.1957056 0.1023058 0.2425019 0.3083380
#>         instructor   sheriff     pilot inspector     mason     baker
#> A_words  0.2034101 0.2256034 0.1339011 0.1741268 0.3154815 0.2847909
#> B_words  0.1903228 0.2029597 0.1112940 0.1272682 0.1585883 0.2981460
#>         administrator    architect collector   operator   surgeon    driver
#> A_words    0.08028339 0.1397101748 0.1572854 0.07317863 0.2337787 0.2733306
#> B_words    0.10544115 0.0008324421 0.1341877 0.08706450 0.1926543 0.2363398
#>           painter conductor     nurse      cook   engineer   retired
#> A_words 0.2703030 0.1832604 0.2187359 0.2278016 0.16052771 0.2494770
#> B_words 0.2413599 0.1034218 0.4470728 0.2849471 0.03511008 0.1146753
#>                sales    lawyer    clergy physician    farmer     clerk
#> A_words -0.006505338 0.2937436 0.1920894 0.1777700 0.3090903 0.2519372
#> B_words  0.032652565 0.2345743 0.2081210 0.1555298 0.2220792 0.3146901
#>            manager     guard    artist      smith  official    police    doctor
#> A_words 0.07080773 0.1948853 0.1819504 0.15938222 0.1300515 0.3116599 0.3413265
#> B_words 0.03393879 0.1344678 0.2274278 0.09691327 0.0743546 0.2590763 0.3390124
#>         professor   student     judge   teacher    author secretary   soldier
#> A_words 0.1604224 0.2540493 0.2008630 0.2675705 0.0828586 0.1211243 0.3599860
#> B_words 0.1368013 0.3299938 0.2493299 0.3567416 0.1224295 0.1220939 0.3076572
plot(res)
```

<img src="man/figures/README-ectplot-1.png" width="100%" />

Effect size can also be calculated. It is the Spearman Correlation
Coefficient of the two rows in `P`. Higher value indicates more
“coherent”, i.e. less bias.

``` r
res
#> 
#> ── sweater object ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> Test type:  ect 
#> Effect size:  0.7001504
#> 
#> ── Functions ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> • `calculate_es()`: Calculate effect size
#> • `plot()`: Plot the bias of each individual word
```

## Example: Relative Negative Sentiment Bias

This analysis attempts to reproduce the analysis in Sweeney & Najafian
(2019).

Please note that the datasets `glove_sweeney`, `bing_pos` and `bing_neg`
are not included in the package. If you are interested in reproducing
the analysis, the 3 datasets are available from
[here](https://github.com/chainsawriot/sweater/tree/master/tests/testdata).

``` r
load("tests/testdata/bing_neg.rda")
load("tests/testdata/bing_pos.rda")
load("tests/testdata/glove_sweeney.rda")

S3 <- c("swedish", "irish", "mexican", "chinese", "filipino",
        "german", "english", "french", "norwegian", "american",
        "indian", "dutch", "russian", "scottish", "italian")
sn <- query(glove_sweeney, S_words = S3, A_words = bing_pos, B_words = bing_neg, method = "rnsb")
```

The analysis shows that `indian`, `mexican`, and `russian` are more
likely to be associated with negative sentiment.

``` r
plot(sn)
```

<img src="man/figures/README-rnsbplot-1.png" width="100%" />

The effect size from the analysis is the Kullback–Leibler divergence of
P from the uniform distribution. It is extremely close to the value
reported in the original paper (0.6225).

``` r
sn
#> 
#> ── sweater object ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> Test type:  rnsb 
#> Effect size:  0.6228853
#> 
#> ── Functions ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> • `calculate_es()`: Calculate effect size
#> • `plot()`: Plot the bias of each individual word
```

## Support for Quanteda Dictionaries

`rnsb` supports [quanteda](https://github.com/quanteda/quanteda)
dictionaries as `S_words`. This support will be expanded to other
methods later.

This analysis uses the data from
[here](https://github.com/chainsawriot/sweater/tree/master/tests/testdata).

For example, `newsmap_europe` is an abridged dictionary from the package
newsmap (Watanabe, 2018). The dictionary contains keywords of European
countries and has two levels: regional level (e.g. Eastern Europe) and
country level (e.g. Germany).

``` r
load("tests/testdata/newsmap_europe.rda")
load("tests/testdata/dictionary_demo.rda")

require(quanteda)
#> Loading required package: quanteda
#> Package version: 3.2.0
#> Unicode version: 13.0
#> ICU version: 66.1
#> Parallel computing: 8 of 8 threads used.
#> See https://quanteda.io for tutorials and examples.
newsmap_europe
#> Dictionary object with 4 primary key entries and 2 nested levels.
#> - [EAST]:
#>   - [BG]:
#>     - bulgaria, bulgarian*, sofia
#>   - [BY]:
#>     - belarus, belarusian*, minsk
#>   - [CZ]:
#>     - czech republic, czech*, prague
#>   - [HU]:
#>     - hungary, hungarian*, budapest
#>   - [MD]:
#>     - moldova, moldovan*, chisinau
#>   - [PL]:
#>     - poland, polish, pole*, warsaw
#>   [ reached max_nkey ... 4 more keys ]
#> - [NORTH]:
#>   - [AX]:
#>     - aland islands, aland island*, alandish, mariehamn
#>   - [DK]:
#>     - denmark, danish, dane*, copenhagen
#>   - [EE]:
#>     - estonia, estonian*, tallinn
#>   - [FI]:
#>     - finland, finnish, finn*, helsinki
#>   - [FO]:
#>     - faeroe islands, faeroe island*, faroese*, torshavn
#>   - [GB]:
#>     - uk, united kingdom, britain, british, briton*, brit*, london
#>   [ reached max_nkey ... 10 more keys ]
#> - [SOUTH]:
#>   - [AD]:
#>     - andorra, andorran*
#>   - [AL]:
#>     - albania, albanian*, tirana
#>   - [BA]:
#>     - bosnia, bosnian*, bosnia and herzegovina, herzegovina, sarajevo
#>   - [ES]:
#>     - spain, spanish, spaniard*, madrid, barcelona
#>   - [GI]:
#>     - gibraltar, gibraltarian*, llanitos
#>   - [GR]:
#>     - greece, greek*, athens
#>   [ reached max_nkey ... 11 more keys ]
#> - [WEST]:
#>   - [AT]:
#>     - austria, austrian*, vienna
#>   - [BE]:
#>     - belgium, belgian*, brussels
#>   - [CH]:
#>     - switzerland, swiss*, zurich, bern
#>   - [DE]:
#>     - germany, german*, berlin, frankfurt
#>   - [FR]:
#>     - france, french*, paris
#>   - [LI]:
#>     - liechtenstein, liechtenstein*, vaduz
#>   [ reached max_nkey ... 3 more keys ]
```

Country-level analysis

``` r
country_level <- rnsb(w = dictionary_demo, S_words = newsmap_europe, A_words = bing_pos, B_words = bing_neg, levels = 2)
plot(country_level)
```

<img src="man/figures/README-rnsb2-1.png" width="100%" />

Region-level analysis

``` r
region_level <- rnsb(w = dictionary_demo, S_words = newsmap_europe, A_words = bing_pos, B_words = bing_neg, levels = 1)
plot(region_level)
```

<img src="man/figures/README-rnsb3-1.png" width="100%" />

Comparison of the two effect sizes. Please note the much smaller effect
size from region-level analysis. It reflects the evener distribution of
P across regions than across countries.

``` r
calculate_es(country_level)
#> [1] 0.0796689
calculate_es(region_level)
#> [1] 0.00329434
```

## Example: Normalized Association Score

Normalized Association Score (Caliskan et al., 2017) is similar to
Relative Norm Distance above.

``` r
S3 <- c("janitor", "statistician", "midwife", "bailiff", "auctioneer", 
"photographer", "geologist", "shoemaker", "athlete", "cashier", 
"dancer", "housekeeper", "accountant", "physicist", "gardener", 
"dentist", "weaver", "blacksmith", "psychologist", "supervisor", 
"mathematician", "surveyor", "tailor", "designer", "economist", 
"mechanic", "laborer", "postmaster", "broker", "chemist", "librarian", 
"attendant", "clerical", "musician", "porter", "scientist", "carpenter", 
"sailor", "instructor", "sheriff", "pilot", "inspector", "mason", 
"baker", "administrator", "architect", "collector", "operator", 
"surgeon", "driver", "painter", "conductor", "nurse", "cook", 
"engineer", "retired", "sales", "lawyer", "clergy", "physician", 
"farmer", "clerk", "manager", "guard", "artist", "smith", "official", 
"police", "doctor", "professor", "student", "judge", "teacher", 
"author", "secretary", "soldier")
A3 <- c("he", "son", "his", "him", "father", "man", "boy", "himself", 
"male", "brother", "sons", "fathers", "men", "boys", "males", 
"brothers", "uncle", "uncles", "nephew", "nephews")
B3 <- c("she", "daughter", "hers", "her", "mother", "woman", "girl", 
"herself", "female", "sister", "daughters", "mothers", "women", 
"girls", "females", "sisters", "aunt", "aunts", "niece", "nieces"
)

nas_f1 <- query(googlenews, S_words= S3, A_words = A3, B_words = B3, method = "nas")
plot(nas_f1)
```

<img src="man/figures/README-nasplot-1.png" width="100%" />

There is a very strong correlation between NAS and RND.

``` r
cor.test(nas_f1$P, garg_f1$P)
#> 
#>  Pearson's product-moment correlation
#> 
#> data:  nas_f1$P and garg_f1$P
#> t = -24.93, df = 74, p-value < 2.2e-16
#> alternative hypothesis: true correlation is not equal to 0
#> 95 percent confidence interval:
#>  -0.9650781 -0.9148179
#> sample estimates:
#>        cor 
#> -0.9453038
```

## Example: Word Embedding Association Test

This example reproduces the detection of “Math. vs Arts” gender bias in
Caliskan et al (2017).

``` r
data(glove_math) # a subset of the original GLoVE word vectors

S4 <- c("math", "algebra", "geometry", "calculus", "equations", "computation", "numbers", "addition")
T4 <- c("poetry", "art", "dance", "literature", "novel", "symphony", "drama", "sculpture")
A4 <- c("male", "man", "boy", "brother", "he", "him", "his", "son")
B4 <- c("female", "woman", "girl", "sister", "she", "her", "hers", "daughter")
sw <- query(glove_math, S4, T4, A4, B4)

# extraction of effect size
sw
#> 
#> ── sweater object ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> Test type:  weat 
#> Effect size:  1.055015
#> 
#> ── Functions ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> • `calculate_es()`: Calculate effect size
#> • `weat_resampling()`: Conduct statistical test
```

## A note about the effect size

By default, the effect size from the function `weat_es` is adjusted by
the pooled standard deviaion (see Page 2 of Caliskan et al. 2007). The
standardized effect size can be interpreted the way as Cohen’s d (Cohen,
1988).

One can also get the unstandardized version (aka. test statistic in the
original paper):

``` r
## weat_es
calculate_es(sw, standardize = FALSE)
#> [1] 0.02486533
```

The original implementation assumes equal size of `S` and `T`. This
assumption can be relaxed by pooling the standard deviaion with sample
size adjustment. The function `weat_es` does it when `S` and `T` are of
different length.

Also, the effect size can be converted to point-biserial correlation
(mathematically equivalent to the Pearson’s product moment correlation).

``` r
weat_es(sw, r = TRUE)
#> [1] 0.4912066
```

## Exact test

The exact test described in Caliskan et al. (2017) is also available.
But it takes a long time to calculate.

``` r
## Don't do it. It takes a long time and is almost always significant.
weat_exact(sw)
```

Instead, please use the resampling approximation of the exact test. The
p-value is very close to the reported 0.018.

``` r
weat_resampling(sw)
#> 
#>  Resampling approximation of the exact test in Caliskan et al. (2017)
#> 
#> data:  sw
#> bias = 0.024865, p-value = 0.0171
#> alternative hypothesis: true bias is greater than 7.245425e-05
#> sample estimates:
#>       bias 
#> 0.02486533
```

## How to get help

  - Read the
    [documentation](https://rdrr.io/github/chainsawriot/sweater/man/)
  - Search for [issues](https://github.com/chainsawriot/sweater/issues)
  - If you have further questions about the package, please contact
    Chung-hong Chan by e-mail, post, or other methods listed on this
    [page](https://www.mzes.uni-mannheim.de/d7/en/profiles/chung-hong-chan).

## Contributing

Contributions in the form of feedback, comments, code, and bug report
are welcome.

  - Fork the source code, modify, and issue a [pull
    request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork).
  - Issues, bug reports: [File a Github
    issue](https://github.com/chainsawriot/sweater/issues).

## Code of Conduct

Please note that the sweater project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## References

1.  An, J., Kwak, H., & Ahn, Y. Y. (2018). SemAxis: A lightweight
    framework to characterize domain-specific word semantics beyond
    sentiment. arXiv preprint arXiv:1806.05521.
2.  Badilla, P., Bravo-Marquez, F., & Pérez, J. (2020). WEFE: The word
    embeddings fairness evaluation framework. In Proceedings of the 29
    th Intern. Joint Conf. Artificial Intelligence.
3.  Brunet, M. E., Alkalay-Houlihan, C., Anderson, A., & Zemel, R.
    (2019, May). Understanding the origins of bias in word embeddings.
    In International Conference on Machine Learning (pp. 803-811). PMLR.
4.  Caliskan, Aylin, Joanna J. Bryson, and Arvind Narayanan. “Semantics
    derived automatically from language corpora contain human-like
    biases.” Science 356.6334 (2017): 183-186.
5.  Cohen, J. (1988), Statistical Power Analysis for the Behavioral
    Sciences, 2nd Edition. Hillsdale: Lawrence Erlbaum.
6.  Dev, S., & Phillips, J. (2019, April). Attenuating bias in word
    vectors. In The 22nd International Conference on Artificial
    Intelligence and Statistics (pp. 879-887). PMLR.
7.  Garg, N., Schiebinger, L., Jurafsky, D., & Zou, J. (2018). Word
    embeddings quantify 100 years of gender and ethnic stereotypes.
    Proceedings of the National Academy of Sciences, 115(16),
    E3635-E3644.
8.  Manzini, T., Lim, Y. C., Tsvetkov, Y., & Black, A. W. (2019). Black
    is to criminal as caucasian is to police: Detecting and removing
    multiclass bias in word embeddings. arXiv preprint arXiv:1904.04047.
9.  McGrath, R. E., & Meyer, G. J. (2006). When effect sizes disagree:
    the case of r and d. Psychological methods, 11(4), 386.
10. Rosenthal, R. (1991), Meta-Analytic Procedures for Social Research.
    Newbury Park: Sage
11. Sweeney, C., & Najafian, M. (2019, July). A transparent framework
    for evaluating unintended demographic bias in word embeddings. In
    Proceedings of the 57th Annual Meeting of the Association for
    Computational Linguistics (pp. 1662-1667).
12. Watanabe, K. (2018). Newsmap: A semi-supervised approach to
    geographical news classification. Digital Journalism, 6(3), 294-309.

-----

1.  In the pre 0.1.0 version of this package, the package used `S`, `T`,
    `A`, and `B` as the main parameters. It was later rejected because
    the symbol `T` is hardlinked to the logical value `TRUE` [as a
    global
    variable](https://stat.ethz.ch/R-manual/R-devel/library/base/html/logical.html);
    and it is considered to be a [bad
    style](https://style.tidyverse.org/syntax.html) to use the symbol
    `T`. Accordingly, they were renamed to `S_words`, `T_words`,
    `A_words`, and `B_words` respectively. But in general, please stop
    using the symbol `T` to represent `TRUE`\!

2.  Please use the `query` function. These functions are kept for
    backward compatibility.

## Benchmark

This is a version of WEAT written entirely in R.

``` r
require(purrr)
```

    ## Loading required package: purrr

``` r
require(lsa)
```

    ## Loading required package: lsa

    ## Loading required package: SnowballC

``` r
take <- function(word, w) {
    return(as.vector(w[word, , drop = FALSE]))
}

get_x <- function(w, words) {
    purrr::map(words, take, w = w)
}

g <- function(c, A, B, w) {
    A_emb <- get_x(w, A)
    B_emb <- get_x(w, B)
    c_emb <- get_x(w, c)[[1]]
    a_cos_diff <- mean(purrr::map_dbl(A_emb, ~cosine(., c_emb)))
    b_cos_diff <- mean(purrr::map_dbl(B_emb, ~cosine(., c_emb)))
    return(a_cos_diff - b_cos_diff)
}

.clean <- function(x, w_lab, verbose = FALSE) {
    new_x <- intersect(x, w_lab)
    if (length(new_x) < length(x) & verbose) {
        print("Some word(s) are not available in w.")
    }
    return(new_x)
}


r_weat <- function(w, S, T, A, B, verbose = FALSE) {
    w_lab <- rownames(w)
    A <- .clean(A, w_lab, verbose = verbose)
    B <- .clean(B, w_lab, verbose = verbose)
    S <- .clean(S, w_lab, verbose = verbose)
    T <- .clean(T, w_lab, verbose = verbose)
    S_diff <- purrr::map_dbl(S, g, A, B, w)
    T_diff <- purrr::map_dbl(T, g, A, B, w)
    ## union_diff <- purrr::map_dbl(union(S, T), g, A, B, w)
    return((mean(S_diff) - mean(T_diff)) / sd(c(S_diff, T_diff)))
}
require(compiler)
```

    ## Loading required package: compiler

``` r
r_weat_c <- cmpfun(r_weat)
```

## The Calikskan et al. example.

``` r
require(sweater)
```

    ## Loading required package: sweater

    ## 
    ## Attaching package: 'sweater'

    ## The following object is masked from 'package:lsa':
    ## 
    ##     query

``` r
S2 <- c("math", "algebra", "geometry", "calculus", "equations",
        "computation", "numbers", "addition")
T2 <- c("poetry", "art", "dance", "literature", "novel", "symphony",
        "drama", "sculpture")
A2 <- c("male", "man", "boy", "brother", "he", "him", "his", "son")
B2 <- c("female", "woman", "girl", "sister", "she", "her", "hers",
        "daughter")
r_weat(glove_math, S2, T2, A2, B2)
```

    ## [1] 1.055015

``` r
r_weat_c(glove_math, S2, T2, A2, B2)
```

    ## [1] 1.055015

The same implementation in C++ from `sweater`

``` r
calculate_es(query(glove_math, S2, T2, A2, B2))
```

    ## [1] 1.055015

``` r
cpp_weat <- function(w, S, T, A, B) {
     calculate_es(query(w, S, T, A, B))
}
```

The C++ implementation in `sweater` is \>10x faster. Byte-code
compilation (`r_weat_c`) can bring about almost no little improvement.

``` r
require(bench)
```

    ## Loading required package: bench

``` r
benchmark_res <- bench::mark(
                            r_weat(glove_math, S2, T2, A2, B2),
                            r_weat_c(glove_math, S2, T2, A2, B2),
                            cpp_weat(glove_math, S2, T2, A2, B2),
                            relative = TRUE)
benchmark_res
```

    ## # A tibble: 3 × 6
    ##   expression                             min median `itr/sec` mem_alloc `gc/sec`
    ##   <bch:expr>                           <dbl>  <dbl>     <dbl>     <dbl>    <dbl>
    ## 1 r_weat(glove_math, S2, T2, A2, B2)    13.5   12.2      1.12      2.35     1   
    ## 2 r_weat_c(glove_math, S2, T2, A2, B2)  13.5   12.6      1         2.35     1.02
    ## 3 cpp_weat(glove_math, S2, T2, A2, B2)   1      1       13.9       1        2.54

### Random benchmark

In this benchmark, we test how the lengths of S/T/A/B affect the
performance. `sweater` is at least 7x faster.

``` r
set.seed(12121)
stab_length <- seq(10, 100, 10)
r_bench <- function(stab_n) {
    w_lab <- rownames(googlenews)
    S <- sample(w_lab, stab_n)
    T <- sample(w_lab, stab_n)
    A <- sample(w_lab, stab_n)
    B <- sample(w_lab, stab_n)
    bench::mark(r_weat(googlenews, S, T, A, B),
                r_weat_c(googlenews, S, T, A, B),
                cpp_weat(googlenews, S, T, A, B),
                relative = TRUE)
}

res <- map(stab_length, r_bench)
```

    ## Warning: Some expressions had a GC in every iteration; so filtering is disabled.
    
    ## Warning: Some expressions had a GC in every iteration; so filtering is disabled.
    
    ## Warning: Some expressions had a GC in every iteration; so filtering is disabled.
    
    ## Warning: Some expressions had a GC in every iteration; so filtering is disabled.
    
    ## Warning: Some expressions had a GC in every iteration; so filtering is disabled.
    
    ## Warning: Some expressions had a GC in every iteration; so filtering is disabled.
    
    ## Warning: Some expressions had a GC in every iteration; so filtering is disabled.

``` r
res %>% map_dfr(~.[1,3]) %>% dplyr::mutate(stab_length = stab_length)
```

    ## # A tibble: 10 × 2
    ##    median stab_length
    ##     <dbl>       <dbl>
    ##  1   9.93          10
    ##  2  11.4           20
    ##  3  11.6           30
    ##  4   7.67          40
    ##  5  12.0           50
    ##  6   9.10          60
    ##  7   9.60          70
    ##  8   8.54          80
    ##  9   9.94          90
    ## 10   9.16         100

### versus WEFE

The hidden gem of this package is the function `read_word2vec`. It is so
flexible and yet speedy. In the following example, we are going to
compare the typical task of the bias detection pipeline:

1.  read a pretrained word embedding file (In this case GloVE, 5.3 GB)
2.  do WEAT

`read_word2vec` is based on the C++ based function `data.table::fread`
and it can read almost all formats (word2vec, glove, fastText, etc). The
entire workflow can be finished in less than a minute.

``` bash
time Rscript bench.R
```

    ## Loading required package: sweater
    ## 
    ## ── sweater object ──────────────────────────────────────────────────────────────
    ## Test type:  weat 
    ## Effect size:  1.055015 
    ## 
    ## ── Functions ───────────────────────────────────────────────────────────────────
    ## • `calculate_es()`: Calculate effect size
    ## • `weat_resampling()`: Conduct statistical test
    ## 
    ## real 0m23.101s
    ## user 0m55.114s
    ## sys  0m6.360s

The Python workflow, however, needs to use `gensim` to read the
pretained word embedding file and it can’t read GLoVE format directly
and the file needs to first convert to word2vec \[1\]. The
`KeyedVectors.load_word2vec_format` is not written in a low-level
language and thus is many times slower than `read_word2vec`. And the
result reported is not the same as the numbers reported in Caliskan et
al.

``` bash
time python3 bench.py
```

    ## {'query_name': 'S and T wrt A and B', 'result': 0.19892264262307435, 'weat': 0.19892264262307435, 'effect_size': 1.0896144272748516, 'p_value': nan}
    ## 
    ## real 12m6.472s
    ## user 11m46.732s
    ## sys  0m19.764s

### versus the original Java code

The original Java code by Caliskan et al. is extremely fast because the
code is highly optimized. If all you need to do is WEAT and you know how
to write Java, it is recommended using the Java code. See this footnote
\[2\] on how to run the code.

``` bash
## javac -cp ./lib/commons-lang3-3.3.2.jar:./lib/commons-math3-3.6.1.jar WeatBenchmark.java Utils.java
time java -classpath .:./lib/commons-lang3-3.3.2.jar:./lib/commons-math3-3.6.1.jar WeatBenchmark
```

    ## Array before check is:[math, algebra, geometry, calculus, equations, computation, numbers, addition]
    ## Array length before check is:8
    ## Array after check is:[math, algebra, geometry, calculus, equations, computation, numbers, addition]
    ## Array length after check is:8
    ## Array before check is:[poetry, art, sculpture, dance, literature, novel, symphony, drama]
    ## Array length before check is:8
    ## Array after check is:[poetry, art, sculpture, dance, literature, novel, symphony, drama]
    ## Array length after check is:8
    ## Array before check is:[brother, male, man, boy, son, he, his, him]
    ## Array length before check is:8
    ## Array after check is:[brother, male, man, boy, son, he, his, him]
    ## Array length after check is:8
    ## Array before check is:[sister, female, woman, girl, daughter, she, hers, her]
    ## Array length before check is:8
    ## Array after check is:[sister, female, woman, girl, daughter, she, hers, her]
    ## Array length after check is:8
    ## The differenceOfMeans is: 0.024865325959943483
    ## Getting the entire distribution...
    ## effectSize: 1.0550147873162645
    ## 
    ## real 0m10.974s
    ## user 0m11.003s
    ## sys  0m1.186s

## Testing environment

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3.10.3
    ## LAPACK: /usr/lib/x86_64-linux-gnu/atlas/liblapack.so.3.10.3
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] compiler  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] bench_1.1.2     sweater_0.1.5   lsa_0.73.2      SnowballC_0.7.0
    ## [5] purrr_0.3.4    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.8       knitr_1.37       magrittr_2.0.2   tidyselect_1.1.1
    ##  [5] R6_2.5.1         rlang_1.0.0      fastmap_1.1.0    fansi_1.0.2     
    ##  [9] dplyr_1.0.7      stringr_1.4.0    tools_4.1.2      xfun_0.29       
    ## [13] utf8_1.2.2       DBI_1.1.2        cli_3.1.1        htmltools_0.5.2 
    ## [17] ellipsis_0.3.2   assertthat_0.2.1 yaml_2.2.2       digest_0.6.29   
    ## [21] tibble_3.1.6     lifecycle_1.0.1  crayon_1.4.2     vctrs_0.3.8     
    ## [25] glue_1.6.1       evaluate_0.14    rmarkdown_2.11   stringi_1.7.6   
    ## [29] pillar_1.6.5     generics_0.1.1   profmem_0.6.0    pkgconfig_2.0.3

``` bash
neofetch --stdout
```

    ## chainsawriot@mzes153 
    ## -------------------- 
    ## OS: Ubuntu 20.04.3 LTS x86_64 
    ## Host: LIFEBOOK U749 10601736746 
    ## Kernel: 5.13.0-27-generic 
    ## Uptime: 9 days, 7 hours, 50 mins 
    ## Packages: 2985 (dpkg), 22 (snap) 
    ## Shell: zsh 5.8 
    ## Resolution: , 1920x1080 
    ## DE: GNOME 
    ## WM: stumpwm 
    ## Theme: Yaru-dark [GTK2/3] 
    ## Icons: Yaru [GTK2/3] 
    ## Terminal: R 
    ## CPU: Intel i5-8365U (8) @ 4.100GHz 
    ## GPU: Intel UHD Graphics 620 
    ## Memory: 4943MiB / 31780MiB

-----

1.  Actually the only difference between the two format is the GLoVE
    format doesn’t record the dimensionality of the matrix in the first
    line.

2.  Please download the [original Java
    code](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/DX4VWP&version=2.0)
    and keep the `Utils.java` in the same directory. Also, put the 2 jar
    files provided inside the `lib` directory. The commented-out line of
    command compiles the `WeatBenchmark.java` to JVM Byecode. It should
    be finished in 1s.
---
title: 'sweater: Speedy Word Embedding Association Test and Extras Using R'
tags:
  - R
  - word embedding
  - implicit bias
  - media bias
  - algorithmic accountability
authors:
  - name: Chung-hong Chan
    orcid: 0000-0002-6232-7530
    affiliation: 1
affiliations:
 - name: Mannheimer Zentrum für Europäische Sozialforschung, Universität Mannheim
   index: 1
citation_author: Chan.
date: 28 January 2022
year: 2022
bibliography: paper.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS
---



# Statement of need

The goal of this R package is to detect associations among words in word embedding spaces. Word embeddings can capture how similar or different two words are in terms of implicit and explicit meanings. Using the example in @collobert2011natural, the word vector for "XBox" is close to that for "PlayStation", as measured by a distance measure such as cosine distance. Word embeddings can also be used to study associations among words that are otherwise difficult to detect. For instance, @jing2021characterizing used word embeddings to study how Democrats and Republicans are divided along party lines about COVID-19.

The same technique can also be used to detect unwanted implicit associations, or biases. For example, @kroon2020guilty measure how close the word vectors for various ethnic group names (e.g. "Dutch", "Belgian" , and "Syrian") are to those for various nouns related to threats (e.g. "terrorist", "murderer", and "gangster"). These biases in word embedding can be understood through the implicit social cognition model of media priming [@arendt:2013:DDM]. In this model, implicit stereotypes are defined as the "strength of the automatic association between a group concept (e.g., minority group) and an attribute (e.g., criminal)." [@arendt:2013:DDM, p. 832] All of these bias detection methods are based on the strength of association between a concept (or a target) and an attribute in embedding spaces.

The importance of detecting biases in word embeddings is twofold. First, pretrained, biased word embeddings deployed in real-life machine learning systems can pose fairness concerns [@packer2018text;@boyarskaya2020overcoming]. Second, biases in word embeddings reflect the biases in the original training material. Social scientists, communication researchers included, have exploited these methods to quantify (implicit) media biases by extracting biases from word embeddings locally trained on large text corpora [e.g. @kroon2020guilty;@knoche2019identifying;@sales2019media].

Previously, the software of these methods is only available piecemeal as the addendum of the original papers and was implemented in different languages (Java, Python, etc.). `sweater` provides several of these bias detection methods in one unified package with a consistent R interface [@rcore]. Also, `sweater` provides several methods that are implemented in C++ for speed and interfaced to R using the `Rcpp` package [@eddelbuettel:2013:SRC] [^BENCH].

[^BENCH]: Compared with a pure R implementation, the C++ implementation of Word Embedding Association Test in `sweater` is at least 7 times faster. See the benchmark [here](https://github.com/chainsawriot/sweater/blob/master/paper/benchmark.md).

# Related work

The R package `cbn` (https://github.com/conjugateprior/cbn) by Will Lowe provides tools for replicating the study by @caliskan:2017:S. The Python package `wefe` [@badilla2020wefe] provides several methods for bias evaluation in a unified (Python) interface.

# Usage

In this section, I demonstrate how the package can be used to detect biases and reproduce some published findings.

## Word Embeddings

The input word embedding $w$ is a dense $m\times n$ matrix, where $m$ is the total size of the vocabulary in the training corpus and $n$ is the vector dimension size.

`sweater` supports input word embeddings, $w$, in several formats. For locally trained word embeddings, output from the R packages `word2vec` [@wijffelsword2vec], `rsparse` [@rsparse] and `text2vec` [@selivanov2020tex2vec] can be used directly with the packages primary functions, such as `query` [^TRAIN]. Pretrained word embeddings in the so-called "word2vec" file format, such as those obtained online [^SOURCE], can be converted to the dense numeric matrix format required with the `read_word2vec` function.

The package also provides three trimmed word embeddings for experimentation: `googlenews` [@mikolov2013distributed], `glove_math` [@pennington:2014:G] , and `small_reddit` [@an2018semaxis].

[^TRAIN]: The vignette of `text2vec` provides a guide on how to locally train word embeddings using the GLoVE algorithm [@pennington:2014:G]. https://cran.r-project.org/web/packages/text2vec/vignettes/glove.html

[^SOURCE]: For example, the [pretrained GLoVE word embeddings](https://nlp.stanford.edu/projects/glove/), [pretrained word2vec word embeddings](https://wikipedia2vec.github.io/wikipedia2vec/pretrained/)  and pretrained [fastText word embeddings](https://fasttext.cc/docs/en/english-vectors.html).

## Query

`sweater` uses the concept of a *query* [@badilla2020wefe] to study associations in $w$ and the $\mathcal{S}\mathcal{T}\mathcal{A}\mathcal{B}$ notation from @brunet2019understanding to form a query. A query contains two or more sets of seed words (wordsets selected by the individual administering the test, sometimes called "seed lexicons" or "dictionaries"). Among these seed wordsets, there should be at least one set of *target words* and one set of *attribute words*.

In the situation of bias detection, target words are words that **should** have no bias and usually represent the concept one would like to probe for biases. For instance, @garg:2018:W investigated the "women bias" of occupation-related words and their target words contain "nurse", "mathematician", and "blacksmith". These words can be used as target words because in an ideal world with no "women bias" associated with occupations, these occupation-related words should have no gender association.

Target words are denoted as wordsets $\mathcal{S}$ and $\mathcal{T}$. All methods require $\mathcal{S}$ while $\mathcal{T}$ is only required for Word Embedding Association Test (WEAT). For instance, the study of gender stereotypes in academic pursuits by @caliskan:2017:S used $\mathcal{S} = \{math, algebra, geometry, calculus, equations, ...\}$ and $\mathcal{T}= \{poetry, art, dance, literature, novel, ...\}$.

In the situation of bias detection, attribute words are words that have known properties in relation to the bias. They are denoted as wordsets $\mathcal{A}$ and $\mathcal{B}$. All methods require both wordsets except Mean Average Cosine Similarity [@manzini2019black]. For instance, the study of gender stereotypes by @caliskan:2017:S used $\mathcal{A} = \{he, son, his, him, ...\}$ and $\mathcal{B} = \{she, daughter, hers, her, ...\}$. In some applications, popular off-the-shelf sentiment dictionaries can also be used as $\mathcal{A}$ and $\mathcal{B}$ [e.g. @sweeney2020reducing]. That being said, it is up to the researchers to select and derive these seed words in a query. However, the selection of seed words has been shown to be the most consequential part of the entire analysis [@antoniak2021bad;@du2021assessing]. Please read @antoniak2021bad for recommendations.

## Supported methods

Table 1 lists all methods supported by sweater. The function `query` is used to conduct a query. The function `calculate_es` can be used for some methods to calculate the effect size representing the overall bias of $w$ from the query.

Table: All methods supported by sweater

| Method                                                          | Target words                 | Attribute words              |
|-----------------------------------------------------------------|------------------------------|------------------------------|
| Mean Average Cosine Similarity [@manzini2019black]              | $\mathcal{S}$                | $\mathcal{A}$                |
| Relative Norm Distance [@garg:2018:W]                           | $\mathcal{S}$                | $\mathcal{A}$, $\mathcal{B}$ |
| Relative Negative Sentiment Bias [@sweeney2020reducing] [^DICT] | $\mathcal{S}$                | $\mathcal{A}$, $\mathcal{B}$ |
| SemAxis [@an2018semaxis]                                        | $\mathcal{S}$                | $\mathcal{A}$, $\mathcal{B}$ |
| Normalized Association Score [@caliskan:2017:S]                 | $\mathcal{S}$                | $\mathcal{A}$, $\mathcal{B}$ |
| Embedding Coherence Test [@dev2019attenuating]                  | $\mathcal{S}$                | $\mathcal{A}$, $\mathcal{B}$ |
| Word Embedding Association Test [@caliskan:2017:S]              | $\mathcal{S}$, $\mathcal{T}$ | $\mathcal{A}$, $\mathcal{B}$ |

[^DICT]: Experimental support for quanteda dictionaries [@benoit2018quanteda] is current available for this method. The support will be expanded to all methods later.

## Example 1

Relative Norm Distance (RND) [@garg:2018:W] is calculated with two sets of attribute words. The following analysis reproduces the calculation of "women bias" values in @garg:2018:W. The publicly available word2vec word embeddings trained on the Google News corpus is used [@mikolov2013distributed]. Words such as "nurse", "midwife" and "librarian" are more associated with female, as indicated by the positive relative norm distance (Figure 1).




```r
library(sweater)
```


```r
S1 <- c("janitor", "statistician", "midwife", "bailiff", "auctioneer",
       "photographer", "geologist", "shoemaker", "athlete", "cashier",
       "dancer", "housekeeper", "accountant", "physicist", "gardener",
       "dentist", "weaver", "blacksmith", "psychologist", "supervisor",
       "mathematician", "surveyor", "tailor", "designer", "economist",
       "mechanic", "laborer", "postmaster", "broker", "chemist",
       "librarian", "attendant", "clerical", "musician", "porter",
       "scientist", "carpenter", "sailor", "instructor", "sheriff",
       "pilot", "inspector", "mason", "baker", "administrator",
       "architect", "collector", "operator", "surgeon", "driver",
       "painter", "conductor", "nurse", "cook", "engineer", "retired",
       "sales", "lawyer", "clergy", "physician", "farmer", "clerk",
       "manager", "guard", "artist", "smith", "official", "police",
       "doctor", "professor", "student", "judge", "teacher", "author",
       "secretary", "soldier")
A1 <- c("he", "son", "his", "him", "father", "man", "boy", "himself",
        "male", "brother", "sons", "fathers", "men", "boys", "males",
        "brothers", "uncle", "uncles", "nephew", "nephews")
B1 <- c("she", "daughter", "hers", "her", "mother", "woman", "girl",
       "herself", "female", "sister", "daughters", "mothers", "women",
       "girls", "females", "sisters", "aunt", "aunts", "niece", "nieces")
res_rnd_male <- query(w = googlenews, S_words = S1,
                      A_words = A1, B_words= B1,
                      method = "rnd")
plot(res_rnd_male)
```

![Bias of words in the target wordset according to relative norm distance](paper_files/figure-latex/rnd-1.pdf) 

## Example 2

Word Embedding Association Test (WEAT) [@caliskan:2017:S] requires all four wordsets of $\mathcal{S}$, $\mathcal{T}$, $\mathcal{A}$, and $\mathcal{B}$. The method is modeled after the Implicit Association Test (IAT) [@nosek:2005:UUI] and it measures the relative strength of $\mathcal{S}$'s association with $\mathcal{A}$ to $\mathcal{B}$ against the same of $\mathcal{T}$. The effect sizes calculated from a large corpus, as shown by @caliskan:2017:S, are comparable to the published IAT effect sizes obtained from volunteers.

In this example, the publicly available GLoVE embeddings made available by the original Stanford Team [@pennington:2014:G] were used. In the following example, the calculation of "Math vs Arts" gender bias in @caliskan:2017:S is reproduced. In this example, the positive effect size indicates the words in the wordset $\mathcal{S}$ are more associated with males than are the words in wordset $\mathcal{T}$.


```r
S2 <- c("math", "algebra", "geometry", "calculus", "equations",
        "computation", "numbers", "addition")
T2 <- c("poetry", "art", "dance", "literature", "novel", "symphony",
        "drama", "sculpture")
A2 <- c("male", "man", "boy", "brother", "he", "him", "his", "son")
B2 <- c("female", "woman", "girl", "sister", "she", "her", "hers",
        "daughter")
sw <- query(w = glove_math,
            S_words = S2, T_words = T2,
            A_words = A2, B_words = B2)
sw
```

```
## 
```

```
## -- sweater object ------------------------------------------------------------------------------------------------------------------------------------------
```

```
## Test type:  weat 
## Effect size:  1.055015
```

```
## 
```

```
## -- Functions -----------------------------------------------------------------------------------------------------------------------------------------------
```

```
## * `calculate_es()`: Calculate effect size
```

```
## * `weat_resampling()`: Conduct statistical test
```

The statistical significance of the effect size can be evaluated using the function `weat_resampling`.


```r
weat_resampling(sw)
```

```
## 
## 	Resampling approximation of the exact test in Caliskan et al. (2017)
## 
## data:  sw
## bias = 0.024865, p-value = 0.0171
## alternative hypothesis: true bias is greater than 7.245425e-05
## sample estimates:
##       bias 
## 0.02486533
```

# Acknowledgements

The development of this package was supported by the Federal Ministry for Family Affairs, Senior Citizens, Women and Youth (*Bundesministerium für Familie, Senioren, Frauen und Jugend*), the Federal Republic of Germany -- Research project: "*Erfahrungen von Alltagsrassismus und medienvermittelter Rassismus in der (politischen) Öffentlichkeit*".

# References
# print method

    Code
      sw
    Message <cliMessage>
      
      -- sweater object --------------------------------------------------------------
    Output
      Test type:  weat 
      Effect size:  1.055015 
    Message <cliMessage>
      
      -- Functions -------------------------------------------------------------------
      * `calculate_es()`: Calculate effect size
      * `weat_resampling()`: Conduct statistical test

---

    Code
      query(glove_math, S_words = S1, A_words = A1, method = "mac")
    Message <cliMessage>
      
      -- sweater object --------------------------------------------------------------
    Output
      Test type:  mac 
      Effect size:  0.09579005 
    Message <cliMessage>
      
      -- Functions -------------------------------------------------------------------
      * `calculate_es()`: Calculate effect size
      * `plot()`: Plot the bias of each individual word

---

    Code
      query(glove_math, S_words = S1, A_words = A1, B_words = B1, method = test_types)
    Message <cliMessage>
      
      -- sweater object --------------------------------------------------------------
    Output
      Test type:  rnd 
      Effect size:  -1.147534 
    Message <cliMessage>
      
      -- Functions -------------------------------------------------------------------
      * `calculate_es()`: Calculate effect size
      * `plot()`: Plot the bias of each individual word

---

    Code
      query(glove_math, S_words = S1, A_words = A1, B_words = B1, method = test_types)
    Message <cliMessage>
      
      -- sweater object --------------------------------------------------------------
    Output
      Test type:  semaxis 
    Message <cliMessage>
      
      -- Functions -------------------------------------------------------------------
      * `plot()`: Plot the bias of each individual word

---

    Code
      query(glove_math, S_words = S1, A_words = A1, B_words = B1, method = test_types)
    Message <cliMessage>
      
      -- sweater object --------------------------------------------------------------
    Output
      Test type:  nas 
    Message <cliMessage>
      
      -- Functions -------------------------------------------------------------------
      * `plot()`: Plot the bias of each individual word

---

    Code
      query(glove_math, S_words = S1, A_words = A1, B_words = B1, method = test_types)
    Message <cliMessage>
      
      -- sweater object --------------------------------------------------------------
    Output
      Test type:  rnsb 
      Effect size:  0.01772161 
    Message <cliMessage>
      
      -- Functions -------------------------------------------------------------------
      * `calculate_es()`: Calculate effect size
      * `plot()`: Plot the bias of each individual word

---

    Code
      query(glove_math, S_words = S1, A_words = A1, B_words = B1, method = test_types)
    Message <cliMessage>
      
      -- sweater object --------------------------------------------------------------
    Output
      Test type:  ect 
      Effect size:  0.9761905 
    Message <cliMessage>
      
      -- Functions -------------------------------------------------------------------
      * `calculate_es()`: Calculate effect size
      * `plot()`: Plot the bias of each individual word

# Display mistake

    Code
      weat_resampling(sw)
    Output
      
      	Resampling approximation of the exact test in Caliskan et al. (2017)
      
      data:  sw
      bias = 0.024865, p-value = 0.0169
      alternative hypothesis: true bias is greater than -3.669546e-05
      sample estimates:
            bias 
      0.02486533 
      

---

    Code
      weat_resampling(whatever)
    Output
      
      	Resampling approximation of the exact test in Caliskan et al. (2017)
      
      data:  whatever
      bias = 0.024865, p-value = 0.0157
      alternative hypothesis: true bias is greater than 8.26938e-05
      sample estimates:
            bias 
      0.02486533 
      

# exact

    Code
      weat_exact(sw2)
    Output
      
      	The exact test in Caliskan et al. (2017)
      
      data:  sw2
      p-value = 0.15
      alternative hypothesis: greater
      

# exact., reverse

    Code
      weat_exact(sw3)
    Output
      
      	The exact test in Caliskan et al. (2017)
      
      data:  sw3
      p-value = 0.8
      alternative hypothesis: greater
      

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
set.seed(46709394)
devtools::load_all()
```

# sweater <img src="man/figures/sweater_logo.svg" align="right" height="200" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/chainsawriot/sweater/workflows/R-CMD-check/badge.svg)](https://github.com/chainsawriot/sweater/actions)
[![Codecov test coverage](https://codecov.io/gh/chainsawriot/sweater/branch/master/graph/badge.svg)](https://app.codecov.io/gh/chainsawriot/sweater?branch=master)
[![CRAN status](https://www.r-pkg.org/badges/version/sweater)](https://CRAN.R-project.org/package=sweater)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04036/status.svg)](https://doi.org/10.21105/joss.04036)
<!-- badges: end -->

The goal of sweater (**S**peedy **W**ord **E**mbedding **A**ssociation **T**est & **E**xtras using **R**) is to test for associations among words in word embedding spaces. The methods provided by this package can also be used to test for unwanted associations, or biases.

The package provides functions that are speedy. They are either implemented in C++, or are speedy but accurate approximation of the original implementation proposed by Caliskan et al (2017). See the benchmark [here](https://github.com/chainsawriot/sweater/blob/master/paper/benchmark.md).

This package provides extra methods such as Relative Norm Distance, Embedding Coherence Test, SemAxis and Relative Negative Sentiment Bias.

If your goal is to reproduce the analysis in Caliskan et al (2017), please consider using the [original Java program](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/DX4VWP&version=2.0) or the R package [cbn](https://github.com/conjugateprior/cbn) by Lowe. To reproduce the analysis in Garg et al (2018), please consider using the [original Python program](https://github.com/nikhgarg/EmbeddingDynamicStereotypes). To reproduce the analysis in Mazini et al (2019), please consider using the [original Python program](https://github.com/TManzini/DebiasMulticlassWordEmbedding/).

Please cite this software as:

Chan, C., (2022). sweater: Speedy Word Embedding Association Test and Extras Using R. Journal of Open Source Software, 7(72), 4036, https://doi.org/10.21105/joss.04036

For a BibTeX entry, use the output from `citation(package = "sweater")`.

## Installation

Recommended: install the latest development version

``` r
remotes::install_github("chainsawriot/sweater")
```

or the "stable" release

```r
install.packages("sweater")
```

## Notation of a query

All tests in this package use the concept of queries (see Badilla et al., 2020) to study associations in the input word embeddings `w`. This package uses the "STAB" notation from Brunet et al (2019). [^1]

[^1]: In the pre 0.1.0 version of this package, the package used `S`, `T`, `A`, and `B` as the main parameters. It was later rejected because the symbol `T` is hardlinked to the logical value `TRUE` [as a global variable](https://stat.ethz.ch/R-manual/R-devel/library/base/html/logical.html); and it is considered to be a [bad style](https://style.tidyverse.org/syntax.html) to use the symbol `T`. Accordingly, they were renamed to `S_words`, `T_words`, `A_words`, and `B_words` respectively. But in general, please stop using the symbol `T` to represent `TRUE`!

All tests depend on two types of words. The first type, namely, `S_words` and `T_words`, is *target words* (or *neutral words* in Garg et al). In the case of studying biases, these are words that **should** have no bias. For instance, the words such as "nurse" and "professor" can be used as target words to study the gender bias in word embeddings. One can also separate these words into two sets, `S_words` and `T_words`, to group words by their perceived bias. For example, Caliskan et al. (2017) grouped target words into two groups: mathematics ("math", "algebra", "geometry", "calculus", "equations", "computation", "numbers", "addition") and arts ("poetry", "art", "dance", "literature", "novel", "symphony", "drama", "sculpture"). Please note that also `T_words` is not always required.

The second type, namely `A_words` and `B_words`, is *attribute words* (or *group words* in Garg et al). These are words with known properties in relation to the bias that one is studying. For example, Caliskan et al. (2017) used gender-related words such as "male", "man", "boy", "brother", "he", "him", "his", "son" to study gender bias. These words qualify as attribute words because we know they are related to a certain gender.

It is recommended to use the function `query()` to make a query and `calculate_es()` to calculate the effect size.

## Available methods

| Target words     | Attribute words  | Method                                                      | `method` argument | Suggested by `query`? | legacy functions [^legacy]                         |
|------------------|------------------|-------------------------------------------------------------|-------------------|-----------------------|----------------------------------------------------|
| S_words          | A_words          | Mean Average Cosine Similarity (Mazini et al. 2019)         | "mac"             | yes                   | mac(), mac_es()                                    |
| S_words          | A_words, B_words | Relative Norm Distance (Garg et al. 2018)                   | "rnd"             | yes                   | rnd(), rnd_es()                                    |
| S_words          | A_words, B_words | Relative Negative Sentiment Bias (Sweeney & Najafian. 2019) | "rnsb"            | no                    | rnsb(), rnsb_es()                                  |
| S_words          | A_words, B_words | Embedding Coherence Test (Dev & Phillips. 2019)             | "ect"             | no                    | ect(), ect_es(), plot_ect()                        |
| S_words          | A_words, B_words | SemAxis (An et al. 2018)                                    | "semaxis"         | no                    | semaxis()                                          |
| S_words          | A_words, B_words | Normalized Association Score (Caliskan et al. 2017)         | "nas"             | no                    | nas()                                              |
| S_words, T_words | A_words, B_words | Word Embedding Association Test (Caliskan et al. 2017)      | "weat"            | yes                   | weat(), weat_es(), weat_resampling(), weat_exact() |
| S_words, T_words | A_words, B_words | Word Embeddings Fairness Evaluation (Badilla et al. 2020)   | To be implemented |                       |                                                    |

## Example: Mean Average Cosine Similarity

The simplest form of bias detection is Mean Average Cosine Similarity (Mazini et al. 2019). The same method is used also in Kroon et al. (2020). `googlenews` is a subset of [the pretrained word2vec word embeddings provided by Google](https://code.google.com/archive/p/word2vec/).

By default, the `query()` function guesses the method you want to use based on the combination of target words and attribute words provided (see the "Suggested?" column in the above table). You can also make this explicit by specifying the `method` argument. Printing the returned object shows the effect size (if available) as well as the functions that can further process the object: `calculate_es` and `plot`. Please read the help file of `calculate_es` (`?calculate_es`) on what is the meaning of the effect size for a specific test.

```{r, eval = FALSE}
require(sweater)
```

```{r mac_neg}
S1 <- c("janitor", "statistician", "midwife", "bailiff", "auctioneer", 
"photographer", "geologist", "shoemaker", "athlete", "cashier", 
"dancer", "housekeeper", "accountant", "physicist", "gardener", 
"dentist", "weaver", "blacksmith", "psychologist", "supervisor", 
"mathematician", "surveyor", "tailor", "designer", "economist", 
"mechanic", "laborer", "postmaster", "broker", "chemist", "librarian", 
"attendant", "clerical", "musician", "porter", "scientist", "carpenter", 
"sailor", "instructor", "sheriff", "pilot", "inspector", "mason", 
"baker", "administrator", "architect", "collector", "operator", 
"surgeon", "driver", "painter", "conductor", "nurse", "cook", 
"engineer", "retired", "sales", "lawyer", "clergy", "physician", 
"farmer", "clerk", "manager", "guard", "artist", "smith", "official", 
"police", "doctor", "professor", "student", "judge", "teacher", 
"author", "secretary", "soldier")

A1 <- c("he", "son", "his", "him", "father", "man", "boy", "himself", 
"male", "brother", "sons", "fathers", "men", "boys", "males", 
"brothers", "uncle", "uncles", "nephew", "nephews")

## The same as:
## mac_neg <- query(googlenews, S_words = S1, A_words = A1, method = "mac")
mac_neg <- query(googlenews, S_words = S1, A_words = A1)
mac_neg
```

The returned object is an S3 object. Please refer to the help file of the method for the definition of all slots (in this case: `?mac`). For example, the magnitude of bias for each word in `S1` is available in the `P` slot. 

```{r mac_neg2}
sort(mac_neg$P)
```

## Example: Relative Norm Distance

This analysis reproduces the analysis in Garg et al (2018), namely Figure 1.

```{r}
B1 <- c("she", "daughter", "hers", "her", "mother", "woman", "girl", 
"herself", "female", "sister", "daughters", "mothers", "women", 
"girls", "females", "sisters", "aunt", "aunts", "niece", "nieces"
)

garg_f1 <- query(googlenews, S_words = S1, A_words = A1, B_words = B1)
garg_f1
```

The object can be plotted by the function `plot` to show the bias of each word in S. Words such as "nurse", "midwife" and "librarian" are more associated with female, as indicated by the positive relative norm distance.

```{r rndplot, fig.height = 12}
plot(garg_f1)
```

The effect size is simply the sum of all relative norm distance values (Equation 3 in Garg et al. 2018). It is displayed simply by printing the object. You can also use the function `calculate_es` to obtain the numeric result.

The more positive effect size indicates that words in `S_words` are more associated with `B_words`. As the effect size is negative, it indicates that the concept of occupation is more associated with `A_words`, i.e. male.

```{r}
calculate_es(garg_f1)
```

## Example: SemAxis

This analysis attempts to reproduce the analysis in An et al. (2018).

You may obtain the word2vec word vectors trained with Trump supporters Reddit from [here](https://github.com/ghdi6758/SemAxis). This package provides a tiny version of the data `small_reddit` for reproducing the analysis.

```{r semxaxisplot}
S2 <- c("mexicans", "asians", "whites", "blacks", "latinos")
A2 <- c("respect")
B2 <- c("disrespect")
res <- query(small_reddit, S_words = S2, A_words = A2, B_words = B2, method = "semaxis", l = 1)
plot(res)
```

## Example: Embedding Coherence Test

Embedding Coherence Test (Dev & Phillips, 2019) is similar to SemAxis. The only significant different is that no "SemAxis" is calculated (the difference between the average word vectors of `A_words` and `B_words`). Instead, it calculates two separate axes for `A_words` and `B_words`. Then it calculates the proximity of each word in `S_words` with the two axes. It is like doing two separate `mac`, but `ect` averages the word vectors of `A_words` / `B_words` first.

It is important to note that `P` is a 2-D matrix. Hence, the plot is 2-dimensional. Words above the equality line are more associated with `B_words` and vice versa.

```{r ectplot}
res <- query(googlenews, S_words = S1, A_words = A1, B_words = B1, method = "ect")
res$P
plot(res)
```

Effect size can also be calculated. It is the Spearman Correlation Coefficient of the two rows in `P`. Higher value indicates more "coherent", i.e. less bias.

```{r}
res
```


## Example: Relative Negative Sentiment Bias

This analysis attempts to reproduce the analysis in Sweeney & Najafian (2019).

Please note that the datasets `glove_sweeney`, `bing_pos` and `bing_neg` are not included in the package. If you are interested in reproducing the analysis, the 3 datasets are available from [here](https://github.com/chainsawriot/sweater/tree/master/tests/testdata).

```{r}
load("tests/testdata/bing_neg.rda")
load("tests/testdata/bing_pos.rda")
load("tests/testdata/glove_sweeney.rda")

S3 <- c("swedish", "irish", "mexican", "chinese", "filipino",
        "german", "english", "french", "norwegian", "american",
        "indian", "dutch", "russian", "scottish", "italian")
sn <- query(glove_sweeney, S_words = S3, A_words = bing_pos, B_words = bing_neg, method = "rnsb")
```

The analysis shows that `indian`, `mexican`, and `russian` are more likely to be associated with negative sentiment.

```{r rnsbplot}
plot(sn)
```

The effect size from the analysis is the Kullback–Leibler divergence of P from the uniform distribution. It is extremely close to the value reported in the original paper (0.6225).

```{r}
sn
```

## Support for Quanteda Dictionaries

`rnsb` supports [quanteda](https://github.com/quanteda/quanteda) dictionaries as `S_words`. This support will be expanded to other methods later.

This analysis uses the data from [here](https://github.com/chainsawriot/sweater/tree/master/tests/testdata).

For example, `newsmap_europe` is an abridged dictionary from the package newsmap (Watanabe, 2018). The dictionary contains keywords of European countries and has two levels: regional level (e.g. Eastern Europe) and country level (e.g. Germany).

```{r}
load("tests/testdata/newsmap_europe.rda")
load("tests/testdata/dictionary_demo.rda")

require(quanteda)
newsmap_europe
```

Country-level analysis

```{r rnsb2, fig.height = 10}
country_level <- rnsb(w = dictionary_demo, S_words = newsmap_europe, A_words = bing_pos, B_words = bing_neg, levels = 2)
plot(country_level)
```

Region-level analysis

```{r rnsb3}
region_level <- rnsb(w = dictionary_demo, S_words = newsmap_europe, A_words = bing_pos, B_words = bing_neg, levels = 1)
plot(region_level)
```

Comparison of the two effect sizes. Please note the much smaller effect size from region-level analysis. It reflects the evener distribution of P across regions than across countries.

```{r}
calculate_es(country_level)
calculate_es(region_level)
```

## Example: Normalized Association Score

Normalized Association Score (Caliskan et al., 2017) is similar to Relative Norm Distance above.

```{r nasplot, fig.height = 12}
S3 <- c("janitor", "statistician", "midwife", "bailiff", "auctioneer", 
"photographer", "geologist", "shoemaker", "athlete", "cashier", 
"dancer", "housekeeper", "accountant", "physicist", "gardener", 
"dentist", "weaver", "blacksmith", "psychologist", "supervisor", 
"mathematician", "surveyor", "tailor", "designer", "economist", 
"mechanic", "laborer", "postmaster", "broker", "chemist", "librarian", 
"attendant", "clerical", "musician", "porter", "scientist", "carpenter", 
"sailor", "instructor", "sheriff", "pilot", "inspector", "mason", 
"baker", "administrator", "architect", "collector", "operator", 
"surgeon", "driver", "painter", "conductor", "nurse", "cook", 
"engineer", "retired", "sales", "lawyer", "clergy", "physician", 
"farmer", "clerk", "manager", "guard", "artist", "smith", "official", 
"police", "doctor", "professor", "student", "judge", "teacher", 
"author", "secretary", "soldier")
A3 <- c("he", "son", "his", "him", "father", "man", "boy", "himself", 
"male", "brother", "sons", "fathers", "men", "boys", "males", 
"brothers", "uncle", "uncles", "nephew", "nephews")
B3 <- c("she", "daughter", "hers", "her", "mother", "woman", "girl", 
"herself", "female", "sister", "daughters", "mothers", "women", 
"girls", "females", "sisters", "aunt", "aunts", "niece", "nieces"
)

nas_f1 <- query(googlenews, S_words= S3, A_words = A3, B_words = B3, method = "nas")
plot(nas_f1)
```

There is a very strong correlation between NAS and RND.

```{r}
cor.test(nas_f1$P, garg_f1$P)
```

## Example: Word Embedding Association Test

This example reproduces the detection of "Math. vs Arts" gender bias in Caliskan et al (2017).

```{r maths}
data(glove_math) # a subset of the original GLoVE word vectors

S4 <- c("math", "algebra", "geometry", "calculus", "equations", "computation", "numbers", "addition")
T4 <- c("poetry", "art", "dance", "literature", "novel", "symphony", "drama", "sculpture")
A4 <- c("male", "man", "boy", "brother", "he", "him", "his", "son")
B4 <- c("female", "woman", "girl", "sister", "she", "her", "hers", "daughter")
sw <- query(glove_math, S4, T4, A4, B4)

# extraction of effect size
sw
```

## A note about the effect size

By default, the effect size from the function `weat_es` is adjusted by the pooled standard deviaion (see Page 2 of Caliskan et al. 2007). The standardized effect size can be interpreted the way as Cohen's d (Cohen, 1988).

One can also get the unstandardized version (aka. test statistic in the original paper):

```{r}
## weat_es
calculate_es(sw, standardize = FALSE)
```

The original implementation assumes equal size of `S` and `T`. This assumption can be relaxed by pooling the standard deviaion with sample size adjustment. The function `weat_es` does it when `S` and `T` are of different length.

Also, the effect size can be converted to point-biserial correlation (mathematically equivalent to the Pearson's product moment correlation).

```{r}
weat_es(sw, r = TRUE)
```

## Exact test

The exact test described in Caliskan et al. (2017) is also available. But it takes a long time to calculate.

```r
## Don't do it. It takes a long time and is almost always significant.
weat_exact(sw)
```

Instead, please use the resampling approximation of the exact test. The p-value is very close to the reported 0.018.

```{r}
weat_resampling(sw)
```

## How to get help

* Read the [documentation](https://rdrr.io/github/chainsawriot/sweater/man/)
* Search for [issues](https://github.com/chainsawriot/sweater/issues)
* If you have further questions about the package, please contact Chung-hong Chan by e-mail, post, or other methods listed on this [page](https://www.mzes.uni-mannheim.de/d7/en/profiles/chung-hong-chan).

## Contributing

Contributions in the form of feedback, comments, code, and bug report are welcome.

* Fork the source code, modify, and issue a [pull request](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork).
* Issues, bug reports: [File a Github issue](https://github.com/chainsawriot/sweater/issues).

## Code of Conduct

Please note that the sweater project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

## References

1. An, J., Kwak, H., & Ahn, Y. Y. (2018). SemAxis: A lightweight framework to characterize domain-specific word semantics beyond sentiment. arXiv preprint arXiv:1806.05521.
2. Badilla, P., Bravo-Marquez, F., & Pérez, J. (2020). WEFE: The word embeddings fairness evaluation framework. In Proceedings of the 29 th Intern. Joint Conf. Artificial Intelligence.
3. Brunet, M. E., Alkalay-Houlihan, C., Anderson, A., & Zemel, R. (2019, May). Understanding the origins of bias in word embeddings. In International Conference on Machine Learning (pp. 803-811). PMLR.
4. Caliskan, Aylin, Joanna J. Bryson, and Arvind Narayanan. "Semantics derived automatically from language corpora contain human-like biases." Science 356.6334 (2017): 183-186.
5. Cohen, J. (1988), Statistical Power Analysis for the Behavioral Sciences, 2nd Edition. Hillsdale: Lawrence Erlbaum.
6. Dev, S., & Phillips, J. (2019, April). Attenuating bias in word vectors. In The 22nd International Conference on Artificial Intelligence and Statistics (pp. 879-887). PMLR.
7. Garg, N., Schiebinger, L., Jurafsky, D., & Zou, J. (2018). Word embeddings quantify 100 years of gender and ethnic stereotypes. Proceedings of the National Academy of Sciences, 115(16), E3635-E3644.
8. Manzini, T., Lim, Y. C., Tsvetkov, Y., & Black, A. W. (2019). Black is to criminal as caucasian is to police: Detecting and removing multiclass bias in word embeddings. arXiv preprint arXiv:1904.04047.
9. McGrath, R. E., & Meyer, G. J. (2006). When effect sizes disagree: the case of r and d. Psychological methods, 11(4), 386.
10. Rosenthal, R. (1991), Meta-Analytic Procedures for Social Research. Newbury Park: Sage
11. Sweeney, C., & Najafian, M. (2019, July). A transparent framework for evaluating unintended demographic bias in word embeddings. In Proceedings of the 57th Annual Meeting of the Association for Computational Linguistics (pp. 1662-1667).
12. Watanabe, K. (2018). Newsmap: A semi-supervised approach to geographical news classification. Digital Journalism, 6(3), 294-309.

---

[^legacy]: Please use the `query` function. These functions are kept for backward compatibility.
---
title: "sweater demo"
output: learnr::tutorial
runtime: shiny_prerendered
---

```{r setup, include=FALSE}
library(learnr)
knitr::opts_chunk$set(echo = FALSE)
```


## read_word2vec

<img src="https://github.com/chainsawriot/sweater/raw/master/man/figures/sweater_logo.svg" width="200">

The main purpose of sweater is to evaluate associations among words in word embedding spaces.

Unwanted associations = biases

You can train your own word embeddings. `sweater` also provides the function `read_word2vec` to read pretrained word embeddings efficiently.

### a demo of read_word2vec

The word embedding file is 5.3GB! gensim needs 10 minutes to read it.

Here < 1min.

```{r read, exercise=TRUE}
require(sweater)
big_glove <- read_word2vec("~/dev/sweater/raw_data/glove.840B.300d.txt")
dim(big_glove)
```

## query

We use the concept of query to look for biases.

A query requires two sets of words: Target words ($\mathcal{S}$, $\mathcal{T}$) and Attribute words ($\mathcal{A}$, $\mathcal{B}$).

**Target words** ($\mathcal{S}$, $\mathcal{T}$) should have no bias.

**Attribute words** ($\mathcal{A}$, $\mathcal{B}$) are attributes in relation to the bias.

| Method                           | Target words                 | Attribute words              |
|----------------------------------|------------------------------|------------------------------|
| Mean Average Cosine Similarity   | $\mathcal{S}$                | $\mathcal{A}$                |
| Relative Norm Distance           | $\mathcal{S}$                | $\mathcal{A}$, $\mathcal{B}$ |
| Relative Negative Sentiment Bias | $\mathcal{S}$                | $\mathcal{A}$, $\mathcal{B}$ |
| SemAxis                          | $\mathcal{S}$                | $\mathcal{A}$, $\mathcal{B}$ |
| Normalized Association Score     | $\mathcal{S}$                | $\mathcal{A}$, $\mathcal{B}$ |
| Embedding Coherence Test         | $\mathcal{S}$                | $\mathcal{A}$, $\mathcal{B}$ |
| Word Embedding Association Test  | $\mathcal{S}$, $\mathcal{T}$ | $\mathcal{A}$, $\mathcal{B}$ |

```r
query(w, S_words, T_words, A_words, B_words, method = "guess", verbose = FALSE)
```

## Example: Garg et al.

Gender biases of occupation words in the pretrained Google News Word Embeddings 

```{r garg, exercise=TRUE, exercise.lines = 30}
require(sweater)
S1 <- c("janitor", "statistician", "midwife", "bailiff", "auctioneer", 
"photographer", "geologist", "shoemaker", "athlete", "cashier", 
"dancer", "housekeeper", "accountant", "physicist", "gardener", 
"dentist", "weaver", "blacksmith", "psychologist", "supervisor", 
"mathematician", "surveyor", "tailor", "designer", "economist", 
"mechanic", "laborer", "postmaster", "broker", "chemist", "librarian", 
"attendant", "clerical", "musician", "porter", "scientist", "carpenter", 
"sailor", "instructor", "sheriff", "pilot", "inspector", "mason", 
"baker", "administrator", "architect", "collector", "operator", 
"surgeon", "driver", "painter", "conductor", "nurse", "cook", 
"engineer", "retired", "sales", "lawyer", "clergy", "physician", 
"farmer", "clerk", "manager", "guard", "artist", "smith", "official", 
"police", "doctor", "professor", "student", "judge", "teacher", 
"author", "secretary", "soldier")

A1 <- c("he", "son", "his", "him", "father", "man", "boy", "himself", 
"male", "brother", "sons", "fathers", "men", "boys", "males", 
"brothers", "uncle", "uncles", "nephew", "nephews")

B1 <- c("she", "daughter", "hers", "her", "mother", "woman", "girl", 
"herself", "female", "sister", "daughters", "mothers", "women", 
"girls", "females", "sisters", "aunt", "aunts", "niece", "nieces"
)
garg_f1 <- query(googlenews, S_words = S1, A_words = A1, B_words = B1)
garg_f1
```

## Plot

```{r garg2, exercise=TRUE, exercise.lines = 30}
require(sweater)
S1 <- c("janitor", "statistician", "midwife", "nurse",
"engineer", "teacher", "author", "secretary", "soldier")

A1 <- c("he", "son", "his", "him", "father", "man", "boy", "himself", 
"male", "brother", "sons", "fathers", "men", "boys", "males", 
"brothers", "uncle", "uncles", "nephew", "nephews")

B1 <- c("she", "daughter", "hers", "her", "mother", "woman", "girl", 
"herself", "female", "sister", "daughters", "mothers", "women", 
"girls", "females", "sisters", "aunt", "aunts", "niece", "nieces"
)
garg_f1 <- query(googlenews, S_words = S1, A_words = A1, B_words = B1)

plot(garg_f1)
```

## Example: Caliskan et al.

Gender biases of academic pursuits in the pretrained GLoVE Word Embeddings 

```{r caliskan, exercise = TRUE, exercise.lines = 30}
require(sweater)
S4 <- c("math", "algebra", "geometry", "calculus", "equations", "computation", "numbers", "addition")
T4 <- c("poetry", "art", "dance", "literature", "novel", "symphony", "drama", "sculpture")
A4 <- c("male", "man", "boy", "brother", "he", "him", "his", "son")
B4 <- c("female", "woman", "girl", "sister", "she", "her", "hers", "daughter")
sw <- query(glove_math, S4, T4, A4, B4)
sw
```

## Effect size

```{r effectsize, exercise = TRUE, exercise.lines = 30}
require(sweater)
S4 <- c("math", "algebra", "geometry", "calculus", "equations", "computation", "numbers", "addition")
T4 <- c("poetry", "art", "dance", "literature", "novel", "symphony", "drama", "sculpture")
A4 <- c("male", "man", "boy", "brother", "he", "him", "his", "son")
B4 <- c("female", "woman", "girl", "sister", "she", "her", "hers", "daughter")
sw <- query(glove_math, S4, T4, A4, B4)
calculate_es(sw)
```

## Significance testing

```{r sig, exercise = TRUE, exercise.lines = 30}
require(sweater)
S4 <- c("math", "algebra", "geometry", "calculus", "equations", "computation", "numbers", "addition")
T4 <- c("poetry", "art", "dance", "literature", "novel", "symphony", "drama", "sculpture")
A4 <- c("male", "man", "boy", "brother", "he", "him", "his", "son")
B4 <- c("female", "woman", "girl", "sister", "she", "her", "hers", "daughter")
sw <- query(glove_math, S4, T4, A4, B4)
weat_resampling(sw)
```

## Supers

* Super efficient (C++, RCpp)
* Super easy to use
* Super nice documentation
* Super nice interface (Guide you what to do next)
* Super nice integration with existing R packages (rsparse, word2vec, text2vec,
quanteda)

### R

```r
require(sweater)

big_glove <- read_word2vec("glove.840B.300d.txt")

S2 <- c("math", "algebra", "geometry", "calculus", "equations",
        "computation", "numbers", "addition")
T2 <- c("poetry", "art", "dance", "literature", "novel", "symphony",
        "drama", "sculpture")
A2 <- c("male", "man", "boy", "brother", "he", "him", "his", "son")
B2 <- c("female", "woman", "girl", "sister", "she", "her", "hers",
        "daughter")
query(big_glove, S2, T2, A2, B2, method = "weat")
```

### Python

```python
from wefe.query import Query
from wefe.word_embedding_model import WordEmbeddingModel
from wefe.metrics.WEAT import WEAT
from gensim.models import KeyedVectors
from gensim.scripts.glove2word2vec import glove2word2vec
from gensim.test.utils import datapath, get_tmpfile

# load_word2vec_format can't read glove file directly.
# if you want to do it "directly"

glove_file = datapath('/home/chainsawriot/dev/sweater/paper/glove.840B.300d.txt')

tmp_file = get_tmpfile("test_word2vec.txt")

_ = glove2word2vec(glove_file, tmp_file)

glove_embeddings = WordEmbeddingModel(KeyedVectors.load_word2vec_format(tmp_file), "word2vec")

target_sets = [["math", "algebra", "geometry", "calculus", "equations",
        "computation", "numbers", "addition"], ["poetry", "art", "dance", "literature", "novel", "symphony",
        "drama", "sculpture"]]
attribute_sets = [["male", "man", "boy", "brother", "he", "him", "his", "son"], ["female", "woman", "girl", "sister", "she", "her", "hers",
        "daughter"]]

attribute_sets_names = ['A', 'B']
target_sets_names = ['S', 'T']

query = Query(target_sets, attribute_sets, target_sets_names, attribute_sets_names)

weat = WEAT()
result = weat.run_query(query, glove_embeddings, calculate_p_value = False)
print(result)
```

### Java

```java
import java.util.Arrays;
public class WeatBenchmark {
	public static void main(String[] args) throws Exception{
    	String semanticModel="glove.840B.300d.txt";
    	int wordDimension =300;
    	String delimiter =" ";	//dimension delimiter in the word embeddings
    	boolean caseSensitive=true; //prefer case sensitivity
    	boolean checkWordPresence=true;
    	int iterations = 1000000;
    	String distribution = "normal";
    	    	
    	String[] target1 = {"math" , "algebra" , "geometry" , "calculus" , "equations" , "computation" , "numbers" , "addition"};
    	String[] target2 = {"poetry" , "art" , "sculpture" , "dance" , "literature" , "novel" , "symphony" , "drama"};
    	String[] attribute1 = {"brother" , "male" , "man" , "boy" , "son" , "he" , "his" , "him"};
    	String[] attribute2 = {"sister" , "female" , "woman" , "girl" , "daughter" , "she" , "hers" , "her"};
	if(checkWordPresence == true){
	    //remove words from categories if they do not exist
	    target1 = Utils.removeCategoryWordsIfNotInDictionary(target1, semanticModel, wordDimension, delimiter, caseSensitive);
	    target2 = Utils.removeCategoryWordsIfNotInDictionary(target2, semanticModel, wordDimension, delimiter, caseSensitive);
	    attribute1 = Utils.removeCategoryWordsIfNotInDictionary(attribute1, semanticModel, wordDimension, delimiter, caseSensitive);
	    attribute2 = Utils.removeCategoryWordsIfNotInDictionary(attribute2, semanticModel, wordDimension, delimiter, caseSensitive);
	}
	double ts = Utils.getTestStatistic(target1, target2, attribute1, attribute2,  caseSensitive,  semanticModel,  wordDimension,  delimiter);
	// Actually this function doesn't need the iterations param,
	// but it is included in the params; but I don't want to modify Utils.java
	double [] entireDistribution = Utils.getEntireDistribution(target1, target2, attribute1, attribute2,  caseSensitive,  semanticModel,  wordDimension,  delimiter, iterations); 
	double es = Utils.effectSize(entireDistribution, ts);
	System.out.println("effectSize: "+ es );
    	}
}
```

## Use it now

```{r install, exercise = TRUE}
install.packages("sweater")
remotes::install_github("chainsawriot/sweater")
```

### Citation:

Chan, C.H. (2022) sweater: Speedy Word Embedding Association Test and Extras Using R. Journal Of Open Source Software (Accepted)

### Funding

The development of this package was supported by the Federal Ministry for Family Affairs, Senior Citizens, Women and Youth (*Bundesministerium für Familie, Senioren, Frauen und Jugend*), the Federal Republic of Germany -- Research project: "*Erfahrungen von Alltagsrassismus und medienvermittelter Rassismus in der (politischen) Öffentlichkeit*".---
title: "resampling"
output: rmarkdown::html_vignette

vignette: >
  %\VignetteIndexEntry{resampling}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The original significance test proposed by Caliskan et al. in their paper is based on exact inference using a permutaion approach. This package provides such an exact test. But like any permutation test (such as Fisher's exact test), the number of permuations, as well as the running time, grows *factorially* --- O(n!) ---  with the number of samples (in the case, the total length of S and T). In the "math vs poetry" example used extensively in this package, the total length of S and T is 16. The total number of iterations needed is 16! = 20922789888000. The number grows 306x to 18! = 6402373705728000 by just adding one item each to S and T. It is unrealistic to expect an ordinary computer can do such computation within a reasonable time frame.

In the following benchmark, we studied how much time does it take to finish an exact test with the word set with a length of 2 to 4. I actually don't recommend trying anything beyond 6; it would take probably a few hours, days even. And actually, the function `weat_exact` will warn you about that.

```{r setup}
library(sweater)
S1 <- c("math", "algebra", "geometry", "calculus", "equations",
 "computation", "numbers", "addition")
T1 <- c("poetry", "art", "dance", "literature", "novel", "symphony", "drama", "sculpture")
A1 <- c("male", "man", "boy", "brother", "he", "him", "his", "son")
B1 <- c("female", "woman", "girl", "sister", "she", "her", "hers", "daughter")

do_exact_weat <- function(n = 4) {
    sw <- weat(glove_math, S1[1:n], T1[1:n], A1[1:n], B1[1:n])
    weat_exact(sw)
}

require(bench)
benchmark_res <- bench::mark(do_exact_weat(2),
                             do_exact_weat(3),
                             do_exact_weat(4),
                             relative = TRUE)
benchmark_res
```

Instead, I recommend using the resampling approximation of the exact inference. Actually, in the [Java code](https://dataverse.harvard.edu/file.xhtml?persistentId=doi:10.7910/DVN/DX4VWP/RQ1CIW&version=2.0) provided by Caliskan et al., the same resampling approximation is used.

```{r}
set.seed(11111)
S1 <- c("math", "algebra", "geometry", "calculus", "equations", "computation", "numbers", "addition")
T1 <- c("poetry", "art", "dance", "literature", "novel", "symphony", "drama", "sculpture")
A1 <- c("male", "man", "boy", "brother", "he", "him", "his", "son")
B1 <- c("female", "woman", "girl", "sister", "she", "her", "hers", "daughter")

sw <- weat(glove_math, S1, T1, A1, B1)
weat_resampling(sw)
```

The method is extremely accurate. In the following simulation, we repeat the same resampling test for 50 times (which is still faster than conducting one exact test). The p-values are extremely close to the 0.018 repeated by Caliskan et al.</p>
<div class="sourceCode" id="cb3"><pre class="sourceCode r"><code class="sourceCode r"><a class="sourceLine" id="cb3-1" title="1">res &lt;-<span class="st"> </span><span class="kw">replicate</span>(<span class="dv">50</span>, <span class="kw">weat_resampling</span>(sw)<span class="op">$</span>p.value)</a>
<a class="sourceLine" id="cb3-2" title="2"><span class="kw">hist</span>(res)</a></code></pre></div>
<p><img src="data:image/png;base64,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" /><!-- --></p>
<p>The approximate can be increased by increasing the <code>n_resampling</code> parameter (default to 9999). In the Java code provided by Caliskan et al., this sets to 100000.</p>
<div class="sourceCode" id="cb4"><pre class="sourceCode r"><code class="sourceCode r"><a class="sourceLine" id="cb4-1" title="1">res &lt;-<span class="st"> </span><span class="kw">replicate</span>(<span class="dv">50</span>, <span class="kw">weat_resampling</span>(sw, <span class="dt">n_resampling =</span> <span class="dv">15999</span>)<span class="op">$</span>p.value)</a>
<a class="sourceLine" id="cb4-2" title="2"><span class="kw">hist</span>(res)</a></code></pre></div>
<p><img src="data:image/png;base64,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" /><!-- --></p>
<p>In the following example, we compare the exact inference with the resampling inference, when n = 5. On a regular computer, the exact test would take 1 minute.</p>
<div class="sourceCode" id="cb5"><pre class="sourceCode r"><code class="sourceCode r"><a class="sourceLine" id="cb5-1" title="1">S2 &lt;-<span class="st"> </span><span class="kw">c</span>(<span class="st">&quot;math&quot;</span>, <span class="st">&quot;algebra&quot;</span>, <span class="st">&quot;geometry&quot;</span>, <span class="st">&quot;calculus&quot;</span>, <span class="st">&quot;equations&quot;</span>)</a>
<a class="sourceLine" id="cb5-2" title="2">T2 &lt;-<span class="st"> </span><span class="kw">c</span>(<span class="st">&quot;poetry&quot;</span>, <span class="st">&quot;art&quot;</span>, <span class="st">&quot;dance&quot;</span>, <span class="st">&quot;literature&quot;</span>, <span class="st">&quot;novel&quot;</span>)</a>
<a class="sourceLine" id="cb5-3" title="3">A2 &lt;-<span class="st"> </span><span class="kw">c</span>(<span class="st">&quot;male&quot;</span>, <span class="st">&quot;man&quot;</span>, <span class="st">&quot;boy&quot;</span>, <span class="st">&quot;brother&quot;</span>, <span class="st">&quot;he&quot;</span>, <span class="st">&quot;him&quot;</span>, <span class="st">&quot;his&quot;</span>, <span class="st">&quot;son&quot;</span>)</a>
<a class="sourceLine" id="cb5-4" title="4">B2 &lt;-<span class="st"> </span><span class="kw">c</span>(<span class="st">&quot;female&quot;</span>, <span class="st">&quot;woman&quot;</span>, <span class="st">&quot;girl&quot;</span>, <span class="st">&quot;sister&quot;</span>, <span class="st">&quot;she&quot;</span>, <span class="st">&quot;her&quot;</span>, <span class="st">&quot;hers&quot;</span>, <span class="st">&quot;daughter&quot;</span>)</a>
<a class="sourceLine" id="cb5-5" title="5"></a>
<a class="sourceLine" id="cb5-6" title="6">sw3 &lt;-<span class="st"> </span><span class="kw">weat</span>(glove_math, T2, S2, A2, B2)</a>
<a class="sourceLine" id="cb5-7" title="7"></a>
<a class="sourceLine" id="cb5-8" title="8"><span class="kw">system.time</span>(exact_res &lt;-<span class="st"> </span><span class="kw">weat_exact</span>(sw3))</a>
<a class="sourceLine" id="cb5-9" title="9"><span class="co">#&gt;    user  system elapsed </span></a>
<a class="sourceLine" id="cb5-10" title="10"><span class="co">#&gt;  74.011   0.748  74.777</span></a></code></pre></div>
<div class="sourceCode" id="cb6"><pre class="sourceCode r"><code class="sourceCode r"><a class="sourceLine" id="cb6-1" title="1">exact_res</a>
<a class="sourceLine" id="cb6-2" title="2"><span class="co">#&gt; </span></a>
<a class="sourceLine" id="cb6-3" title="3"><span class="co">#&gt;  The exact test in Caliskan et al. (2017)</span></a>
<a class="sourceLine" id="cb6-4" title="4"><span class="co">#&gt; </span></a>
<a class="sourceLine" id="cb6-5" title="5"><span class="co">#&gt; data:  sw3</span></a>
<a class="sourceLine" id="cb6-6" title="6"><span class="co">#&gt; p-value = 0.9802</span></a>
<a class="sourceLine" id="cb6-7" title="7"><span class="co">#&gt; alternative hypothesis: greater</span></a></code></pre></div>
<div class="sourceCode" id="cb7"><pre class="sourceCode r"><code class="sourceCode r"><a class="sourceLine" id="cb7-1" title="1">res &lt;-<span class="st"> </span><span class="kw">replicate</span>(<span class="dv">100</span>, <span class="kw">weat_resampling</span>(sw3)<span class="op">$</span>p.value)</a>
<a class="sourceLine" id="cb7-2" title="2"><span class="kw">hist</span>(res, <span class="dt">main =</span> <span class="st">&quot;Distribution of p-values&quot;</span>)</a>
<a class="sourceLine" id="cb7-3" title="3"><span class="kw">abline</span>(<span class="dt">v =</span> exact_res<span class="op">$</span>p.value, <span class="dt">col =</span> <span class="st">&quot;red&quot;</span>, <span class="dt">lwd =</span> <span class="dv">3</span>)</a></code></pre></div>
<p><img src="data:image/png;base64,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" /><!-- --></p>


<!-- The method is extremely accurate. In the following simulation, we repeat the same resampling test for 50 times (which is still faster than conducting one exact test). The p-values are extremely close to the 0.018 repeated by Caliskan et al. -->

<!-- ```{r} -->
<!-- res <- replicate(50, weat_resampling(sw)$p.value) -->
<!-- hist(res) -->
<!-- ``` -->

<!-- The approximate can be increased by increasing the `n_resampling` parameter (default to 9999). In the Java code provided by Caliskan et al., this sets to 100000. -->

<!-- ```{r} -->
<!-- res <- replicate(50, weat_resampling(sw, n_resampling = 15999)$p.value) -->
<!-- hist(res) -->
<!-- ``` -->

<!-- In the following example, we compare the exact inference with the resampling inference, when n = 5. On a regular computer, the exact test would take 1 minute. -->

<!-- ```{r, eval = FALSE} -->
<!-- S2 <- c("math", "algebra", "geometry", "calculus", "equations") -->
<!-- T2 <- c("poetry", "art", "dance", "literature", "novel") -->
<!-- A2 <- c("male", "man", "boy", "brother", "he", "him", "his", "son") -->
<!-- B2 <- c("female", "woman", "girl", "sister", "she", "her", "hers", "daughter") -->

<!-- sw3 <- weat(glove_math, T2, S2, A2, B2) -->

<!-- system.time({exact_res <- weat_exact(sw3)}) -->
<!-- ``` -->

<!-- ```{r} -->
<!-- exact_res -->
<!-- ``` -->

<!-- ```{r} -->
<!-- res <- replicate(100, weat_resampling(sw3)$p.value) -->
<!-- hist(res, main = "Distribution of p-values") -->
<!-- abline(v = exact_res$p.value, col = "red", lwd = 3) -->
<!-- ``` -->
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/semaxis.R
\docType{data}
\name{small_reddit}
\alias{small_reddit}
\title{A subset of the pretrained word2vec word vectors on Reddit}
\format{
An object of class \code{matrix} (inherits from \code{array}) with 106 rows and 300 columns.
}
\usage{
small_reddit
}
\description{
This is a subset of the pretrained word2vec word vectors on Reddit provided by An et al. (2018). With this dataset, you can try with the "l" parameter of \code{\link[=semaxis]{semaxis()}} up to 10.
}
\references{
An, J., Kwak, H., & Ahn, Y. Y. (2018). \href{https://arxiv.org/abs/1806.05521}{SemAxis: A lightweight framework to characterize domain-specific word semantics beyond sentiment.} arXiv preprint arXiv:1806.05521.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\docType{data}
\name{googlenews}
\alias{googlenews}
\title{A subset of the pretrained word2vec word vectors}
\format{
An object of class \code{matrix} (inherits from \code{array}) with 116 rows and 300 columns.
}
\usage{
googlenews
}
\description{
This is a subset of the original pretrained word2vec word vectors trained on Google News. The same word vectors were used in Garg et al. (2018) to study biases.
}
\references{
Garg, N., Schiebinger, L., Jurafsky, D., & Zou, J. (2018). Word embeddings quantify 100 years of gender and ethnic stereotypes. Proceedings of the National Academy of Sciences, 115(16), E3635-E3644. \doi{10.1073/pnas.1720347115}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{calculate_es}
\alias{calculate_es}
\title{Calculate the effect size of a query}
\usage{
calculate_es(x, ...)
}
\arguments{
\item{x}{an S3 object returned from a query, either by the function \code{\link[=query]{query()}} or underlying functions such as \code{\link[=mac]{mac()}}}

\item{...}{additional parameters for the effect size functions
\describe{
\item{\code{r}}{for \code{weat}: a boolean to denote whether convert the effect size to biserial correlation coefficient.}
\item{\code{standardize}}{for \code{weat}: a boolean to denote whether to correct the difference by the standard division. The standardized version can be interpreted the same way as Cohen's d. }
}}
}
\value{
effect size
}
\description{
This function calculates the effect of a query.
}
\details{
The following methods are supported.
\describe{
\item{\code{mac}}{mean cosine distance value. The value makes sense only for comparison (e.g. before and after debiasing). But a lower value indicates greater association between the target words and the attribute words.}
\item{\code{rnd}}{sum of all relative norm distances. It equals to zero when there is no bias.}
\item{\code{rnsb}}{Kullback-Leibler divergence of the predicted negative probabilities, P, from the uniform distribution. A lower value indicates less bias.}
\item{\code{ect}}{Spearman Coefficient of an Embedding Coherence Test. The value ranges from -1 to +1 and a larger value indicates less bias.}
\item{\code{weat}}{The standardized effect size (default) can be interpreted the same way as Cohen's D.}
}
}
\references{
Caliskan, A., Bryson, J. J., & Narayanan, A. (2017). Semantics derived automatically from language corpora contain human-like biases. Science, 356(6334), 183-186. \doi{10.1126/science.aal4230}
Dev, S., & Phillips, J. (2019, April). \href{https://proceedings.mlr.press/v89/dev19a.html}{Attenuating bias in word vectors.} In The 22nd International Conference on Artificial Intelligence and Statistics (pp. 879-887). PMLR.
Garg, N., Schiebinger, L., Jurafsky, D., & Zou, J. (2018). Word embeddings quantify 100 years of gender and ethnic stereotypes. Proceedings of the National Academy of Sciences, 115(16), E3635-E3644. \doi{10.1073/pnas.1720347115}
Manzini, T., Lim, Y. C., Tsvetkov, Y., & Black, A. W. (2019). \href{https://arxiv.org/abs/1904.04047}{Black is to criminal as caucasian is to police: Detecting and removing multiclass bias in word embeddings.} arXiv preprint arXiv:1904.04047.
Sweeney, C., & Najafian, M. (2019, July). \href{https://aclanthology.org/P19-1162/}{A transparent framework for evaluating unintended demographic bias in word embeddings.} In Proceedings of the 57th Annual Meeting of the Association for Computational Linguistics (pp. 1662-1667).
}
\seealso{
\code{\link[=weat_es]{weat_es()}}, \code{\link[=mac_es]{mac_es()}}, \code{\link[=rnd_es]{rnd_es()}}, \code{\link[=rnsb_es]{rnsb_es()}}, \code{\link[=ect_es]{ect_es()}}
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnd.R
\name{rnd}
\alias{rnd}
\title{Relative Norm Distance}
\usage{
rnd(w, S_words, A_words, B_words, verbose = FALSE)
}
\arguments{
\item{w}{a numeric matrix of word embeddings (e.g. from \code{\link[rsparse:GloVe]{rsparse::GloVe()}})}

\item{S_words}{a character vector of the first set of target words. In an example of studying gender stereotype, it can include occupations such as programmer, engineer, scientists...}

\item{A_words}{a character vector of the first set of attribute words. In an example of studying gender stereotype, it can include words such as man, male, he, his.}

\item{B_words}{a character vector of the second set of attribute words. In an example of studying gender stereotype, it can include words such as woman, female, she, her.}

\item{verbose}{logical, whether to display information}
}
\value{
A list with class \code{"rnd"} containing the following components:
\describe{
\item{\code{$norm_diff}}{a vector of relative norm distances for every word in S_words}
\item{\code{$S_words}}{the input S_words}
\item{\code{$A_words}}{the input A_words}
\item{\code{$B_words}}{the input B_words}
}
\code{\link{rnd_es}} can be used to obtain the effect size of the test.
}
\description{
This function calculate the relative norm distance (RND) of word embeddings. If possible, please use \code{\link[=query]{query()}} instead.
}
\examples{
data(googlenews)
S1 <- c("janitor", "statistician", "midwife", "bailiff", "auctioneer",
"photographer", "geologist", "shoemaker", "athlete", "cashier", "dancer",
"housekeeper", "accountant", "physicist", "gardener", "dentist", "weaver",
"blacksmith", "psychologist", "supervisor", "mathematician", "surveyor",
"tailor", "designer", "economist", "mechanic", "laborer", "postmaster",
"broker", "chemist", "librarian", "attendant", "clerical", "musician",
"porter", "scientist", "carpenter", "sailor", "instructor", "sheriff",
"pilot", "inspector", "mason", "baker", "administrator", "architect",
"collector", "operator", "surgeon", "driver", "painter", "conductor",
"nurse", "cook", "engineer", "retired", "sales", "lawyer", "clergy",
"physician", "farmer", "clerk", "manager", "guard", "artist", "smith",
"official", "police", "doctor", "professor", "student", "judge",
"teacher", "author", "secretary", "soldier")
A1 <- c("he", "son", "his", "him", "father", "man", "boy", "himself",
"male", "brother", "sons", "fathers", "men", "boys", "males", "brothers",
"uncle", "uncles", "nephew", "nephews")
B1 <- c("she", "daughter", "hers", "her", "mother", "woman", "girl",
"herself", "female", "sister", "daughters", "mothers", "women", "girls",
"females", "sisters", "aunt", "aunts", "niece", "nieces")
garg_f1 <- rnd(googlenews, S1, A1, B1)
plot_bias(garg_f1)
}
\references{
Garg, N., Schiebinger, L., Jurafsky, D., & Zou, J. (2018). Word embeddings quantify 100 years of gender and ethnic stereotypes. Proceedings of the National Academy of Sciences, 115(16), E3635-E3644. \doi{10.1073/pnas.1720347115}
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ect.R
\name{plot_ect}
\alias{plot_ect}
\title{Plot an ECT result on a two-dimensional plane}
\usage{
plot_ect(x, ...)
}
\arguments{
\item{x}{an ect object from the \link{ect} function.}

\item{...}{additional parameters to the underlying \code{\link[=plot]{plot()}} function}
}
\value{
a plot
}
\description{
This functions plot the words in \code{S_words} on a 2D plane according to their association with the average vectors of \code{A_words} and \code{B_words}. A equality line is also added. Words along the equality line have less bias. Words located on the upper side of the equality line have a stronger association with \code{A_words} and vice versa.
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ect.R
\name{ect}
\alias{ect}
\title{Embedding Coherence Test}
\usage{
ect(w, S_words, A_words, B_words, verbose = FALSE)
}
\arguments{
\item{w}{a numeric matrix of word embeddings (e.g. from \code{\link[rsparse:GloVe]{rsparse::GloVe()}})}

\item{S_words}{a character vector of the first set of target words. In an example of studying gender stereotype, it can include occupations such as programmer, engineer, scientists...}

\item{A_words}{a character vector of the first set of attribute words. In an example of studying gender stereotype, it can include words such as man, male, he, his.}

\item{B_words}{a character vector of the second set of attribute words. In an example of studying gender stereotype, it can include words such as woman, female, she, her.}

\item{verbose}{logical, whether to display information}
}
\value{
A list with class \code{"ect"} containing the following components:
\describe{
\item{\code{$A_words}}{the input A_words}
\item{\code{$B_words}}{the input B_words}
\item{\code{$S_words}}{the input S_words}
\item{\code{$u_a}}{Cosine similarity between each word vector of S_words and average vector of A_words}
\item{\code{$u_b}}{Cosine similarity between each word vector of S_words and average vector of B_words}
}
}
\description{
This function estimate the Embedding Coherence Test (ECT) of word embeddings (Dev & Philips, 2019). If possible, please use \code{\link[=query]{query()}} instead.
}
\examples{
data(googlenews)
S1 <- c("janitor", "statistician", "midwife", "bailiff", "auctioneer",
"photographer", "geologist", "shoemaker", "athlete", "cashier", "dancer",
"housekeeper", "accountant", "physicist", "gardener", "dentist", "weaver",
"blacksmith", "psychologist", "supervisor", "mathematician", "surveyor",
"tailor", "designer", "economist", "mechanic", "laborer", "postmaster",
"broker", "chemist", "librarian", "attendant", "clerical", "musician",
"porter", "scientist", "carpenter", "sailor", "instructor", "sheriff",
"pilot", "inspector", "mason", "baker", "administrator", "architect",
"collector", "operator", "surgeon", "driver", "painter", "conductor",
"nurse", "cook", "engineer", "retired", "sales", "lawyer", "clergy",
"physician", "farmer", "clerk", "manager", "guard", "artist", "smith",
"official", "police", "doctor", "professor", "student", "judge",
"teacher", "author", "secretary", "soldier")
A1 <- c("he", "son", "his", "him", "father", "man", "boy", "himself",
"male", "brother", "sons", "fathers", "men", "boys", "males", "brothers",
"uncle", "uncles", "nephew", "nephews")
B1 <- c("she", "daughter", "hers", "her", "mother", "woman", "girl",
"herself", "female", "sister", "daughters", "mothers", "women", "girls",
"females", "sisters", "aunt", "aunts", "niece", "nieces")
garg_f1 <- ect(googlenews, S1, A1, B1)
plot_ect(garg_f1)
}
\references{
Dev, S., & Phillips, J. (2019, April). \href{https://proceedings.mlr.press/v89/dev19a.html}{Attenuating bias in word vectors.} In The 22nd International Conference on Artificial Intelligence and Statistics (pp. 879-887). PMLR.
}
\seealso{
\code{\link[=ect_es]{ect_es()}} can be used to obtain the effect size of the test.
\code{\link[=plot_ect]{plot_ect()}} can be used to visualize the result.
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sweater.R
\docType{data}
\name{glove_math}
\alias{glove_math}
\title{A subset of the pretrained GLoVE word vectors}
\format{
An object of class \code{matrix} (inherits from \code{array}) with 32 rows and 300 columns.
}
\usage{
glove_math
}
\description{
This is a subset of the original pretrained GLoVE word vectors provided by Pennington et al (2017). The same word vectors were used in Caliskan et al. (2017) to study biases.
}
\references{
Pennington, J., Socher, R., & Manning, C. D. (2014, October). \href{https://aclanthology.org/D14-1162/}{Glove: Global vectors for word representation.} In Proceedings of the 2014 conference on empirical methods in natural language processing (EMNLP) (pp. 1532-1543).
Caliskan, A., Bryson, J. J., & Narayanan, A. (2017). Semantics derived automatically from language corpora contain human-like biases. Science, 356(6334), 183-186. \doi{10.1126/science.aal4230}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/semaxis.R
\name{semaxis}
\alias{semaxis}
\title{Characterise word semantics using the SemAxis framework}
\usage{
semaxis(w, S_words, A_words, B_words, l = 0, verbose = FALSE)
}
\arguments{
\item{w}{a numeric matrix of word embeddings (e.g. from \code{\link[rsparse:GloVe]{rsparse::GloVe()}})}

\item{S_words}{a character vector of the first set of target words. In an example of studying gender stereotype, it can include occupations such as programmer, engineer, scientists...}

\item{A_words}{a character vector of the first set of attribute words. In an example of studying gender stereotype, it can include words such as man, male, he, his.}

\item{B_words}{a character vector of the second set of attribute words. In an example of studying gender stereotype, it can include words such as woman, female, she, her.}

\item{l}{an integer indicates the number of words to augment each word in A and B based on cosine , see An et al (2018). Default to 0 (no augmentation).}

\item{verbose}{logical, whether to display information}
}
\value{
A list with class \code{"semaxis"} containing the following components:
\describe{
\item{\code{$P}}{for each of words in S, the score according to SemAxis}
\item{\code{$V}}{the semantic axis vector}
\item{\code{$S_words}}{the input S_words}
\item{\code{$A_words}}{the input A_words}
\item{\code{$B_words}}{the input B_words}
}
}
\description{
This function calculates the axis and the score using the SemAxis framework proposed in An et al (2018). If possible, please use \code{\link[=query]{query()}} instead.
}
\examples{
data(glove_math)
S1 <- c("math", "algebra", "geometry", "calculus", "equations",
"computation", "numbers", "addition")
A1 <- c("male", "man", "boy", "brother", "he", "him", "his", "son")
B1 <- c("female", "woman", "girl", "sister", "she", "her", "hers", "daughter")
semaxis(glove_math, S1, A1, B1, l = 0)$P
}
\references{
An, J., Kwak, H., & Ahn, Y. Y. (2018). \href{https://arxiv.org/abs/1806.05521}{SemAxis: A lightweight framework to characterize domain-specific word semantics beyond sentiment.} arXiv preprint arXiv:1806.05521.
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sweater.R
\name{weat_es}
\alias{weat_es}
\title{Calculation of WEAT effect size}
\usage{
weat_es(x, standardize = TRUE, r = FALSE)
}
\arguments{
\item{x}{an object from the \link{weat} function.}

\item{standardize}{a boolean to denote whether to correct the difference by the standard division. The standardized version can be interpreted the same way as Cohen's d.}

\item{r}{a boolean to denote whether convert the effect size to biserial correlation coefficient.}
}
\value{
the effect size of the query
}
\description{
This function calculates the effect size from a sweater object. The original implementation in Caliskan et al. (2017) assumes the numbers of words in S and in T must be equal. The current implementation eases this assumption by adjusting the variance with the difference in sample sizes. It is also possible to convert the Cohen's d to Pearson's correlation coefficient (r). If possible, please use \code{\link[=calculate_es]{calculate_es()}} instead.
}
\examples{
# Reproduce the number in Caliskan et al. (2017) - Table 1, "Math vs. Arts"
data(glove_math)
S1 <- c("math", "algebra", "geometry", "calculus", "equations",
"computation", "numbers", "addition")
T1 <- c("poetry", "art", "dance", "literature", "novel", "symphony", "drama", "sculpture")
A1 <- c("male", "man", "boy", "brother", "he", "him", "his", "son")
B1 <- c("female", "woman", "girl", "sister", "she", "her", "hers", "daughter")
sw <- weat(glove_math, S1, T1, A1, B1)
weat_es(sw)
}
\references{
Caliskan, A., Bryson, J. J., & Narayanan, A. (2017). Semantics derived automatically from language corpora contain human-like biases. Science, 356(6334), 183-186. \doi{10.1126/science.aal4230}
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnsb.R
\name{rnsb_es}
\alias{rnsb_es}
\title{Calculation the Kullback-Leibler divergence}
\usage{
rnsb_es(x)
}
\arguments{
\item{x}{an rnsb object from the \link{rnsb} function.}
}
\value{
the Kullback-Leibler divergence.
}
\description{
This function calculates the Kullback-Leibler divergence of the predicted negative probabilities, P, from the uniform distribution. If possible, please use \code{\link[=calculate_es]{calculate_es()}} instead.
}
\references{
Sweeney, C., & Najafian, M. (2019, July). \href{https://aclanthology.org/P19-1162/}{A transparent framework for evaluating unintended demographic bias in word embeddings.} In Proceedings of the 57th Annual Meeting of the Association for Computational Linguistics (pp. 1662-1667).
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nas.R
\name{nas}
\alias{nas}
\title{Calculate Normalized Association Score}
\usage{
nas(w, S_words, A_words, B_words, verbose = FALSE)
}
\arguments{
\item{w}{a numeric matrix of word embeddings (e.g. from \code{\link[rsparse:GloVe]{rsparse::GloVe()}})}

\item{S_words}{a character vector of the first set of target words. In an example of studying gender stereotype, it can include occupations such as programmer, engineer, scientists...}

\item{A_words}{a character vector of the first set of attribute words. In an example of studying gender stereotype, it can include words such as man, male, he, his.}

\item{B_words}{a character vector of the second set of attribute words. In an example of studying gender stereotype, it can include words such as woman, female, she, her.}

\item{verbose}{logical, whether to display information}
}
\value{
A list with class \code{"nas"} containing the following components:
\describe{
\item{\code{$P}}{a vector of normalized association score for every word in S}
\item{\code{$raw}}{a list of raw results used for calculating normalized association scores}
\item{\code{$S_words}}{the input S_words}
\item{\code{$A_words}}{the input A_words}
\item{\code{$B_words}}{the input B_words}
}
}
\description{
This functions quantifies the bias in a set of word embeddings by Caliskan et al (2017). In comparison to WEAT introduced in the same paper, this method is more suitable for continuous ground truth data. See Figure 1 and Figure 2 of the original paper. If possible, please use \code{\link[=query]{query()}} instead.
}
\references{
Caliskan, A., Bryson, J. J., & Narayanan, A. (2017). Semantics derived automatically from language corpora contain human-like biases. Science, 356(6334), 183-186. \doi{10.1126/science.aal4230}
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sweater.R
\name{weat}
\alias{weat}
\title{Speedy Word Embedding Association Test}
\usage{
weat(w, S_words, T_words, A_words, B_words, verbose = FALSE)
}
\arguments{
\item{w}{a numeric matrix of word embeddings (e.g. from \code{\link[rsparse:GloVe]{rsparse::GloVe()}})}

\item{S_words}{a character vector of the first set of target words. In an example of studying gender stereotype, it can include occupations such as programmer, engineer, scientists...}

\item{T_words}{a character vector of the second set of target words. In an example of studying gender stereotype, it can include occupations such as nurse, teacher, librarian...}

\item{A_words}{a character vector of the first set of attribute words. In an example of studying gender stereotype, it can include words such as man, male, he, his.}

\item{B_words}{a character vector of the second set of attribute words. In an example of studying gender stereotype, it can include words such as woman, female, she, her.}

\item{verbose}{logical, whether to display information}
}
\value{
A list with class \code{"weat"} containing the following components:
\describe{
\item{\code{$S_diff}}{for each of words in S_words, mean of the mean differences in cosine similarity between words in A_words and words in B_words}
\item{\code{$T_diff}}{for each of words in T_words, mean of the mean differences in cosine similarity between words in A_words and words in B_words}
\item{\code{$S_words}}{the input S_words}
\item{\code{$T_words}}{the input T_words}
\item{\code{$A_words}}{the input A_words}
\item{\code{$B_words}}{the input B_words}
}
\code{\link{weat_es}} can be used to obtain the effect size of the test; \code{\link{weat_resampling}} for a test of significance.
}
\description{
This functions test the bias in a set of word embeddings using the method by Caliskan et al (2017). If possible, please use \code{\link[=query]{query()}} instead.
}
\examples{
# Reproduce the number in Caliskan et al. (2017) - Table 1, "Math vs. Arts"
data(glove_math)
S1 <- c("math", "algebra", "geometry", "calculus", "equations",
"computation", "numbers", "addition")
T1 <- c("poetry", "art", "dance", "literature", "novel", "symphony", "drama", "sculpture")
A1 <- c("male", "man", "boy", "brother", "he", "him", "his", "son")
B1 <- c("female", "woman", "girl", "sister", "she", "her", "hers", "daughter")
sw <- weat(glove_math, S1, T1, A1, B1)
weat_es(sw)
}
\references{
Caliskan, A., Bryson, J. J., & Narayanan, A. (2017). Semantics derived automatically from language corpora contain human-like biases. Science, 356(6334), 183-186. \doi{10.1126/science.aal4230}
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R, R/s3.R
\name{plot_bias}
\alias{plot_bias}
\alias{plot.sweater}
\title{Visualize the bias of words in S}
\usage{
plot_bias(x)

\method{plot}{sweater}(x, ...)
}
\arguments{
\item{x}{an S3 object returned from mac, rnd, semaxis, nas or rnsb}

\item{...}{other parameters}
}
\value{
a plot
}
\description{
For \code{ect}, this function calls \code{\link[=plot_ect]{plot_ect()}}. For other tests (except \code{weat}), this function plots the bias of words in \code{S} as a Cleveland Dot Plot. Plotting the result of \code{weat} is not supported.
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnsb.R
\name{rnsb}
\alias{rnsb}
\title{Relative Negative Sentiment Bias}
\usage{
rnsb(w, S_words, A_words, B_words, levels = 1, verbose = FALSE)
}
\arguments{
\item{w}{a numeric matrix of word embeddings (e.g. from \code{\link[rsparse:GloVe]{rsparse::GloVe()}})}

\item{S_words}{a character vector of the first set of target words. In an example of studying gender stereotype, it can include occupations such as programmer, engineer, scientists...}

\item{A_words}{a character vector of the first set of attribute words. In an example of studying gender stereotype, it can include words such as man, male, he, his.}

\item{B_words}{a character vector of the second set of attribute words. In an example of studying gender stereotype, it can include words such as woman, female, she, her.}

\item{levels}{levels of entries in a hierarchical dictionary that will be applied (see \code{\link[quanteda:dfm_lookup]{quanteda::dfm_lookup()}})}

\item{verbose}{logical, whether to display information}
}
\value{
A list with class \code{"rnsb"} containing the following components:
\describe{
\item{\code{$classifer}}{ a logistic regression model with L2 regularization trained with LiblineaR}
\item{\code{$A_words}}{the input A_words}
\item{\code{$B_words}}{the input B_words}
\item{\code{$S_words}}{the input S_words}
\item{\code{$P}}{the predicted negative sentiment probabilities}
}
\code{\link{rnsb_es}} can be used to obtain the effect size of the test.
}
\description{
This function estimate the Relative Negative Sentiment Bias (RNSB) of word embeddings (Sweeney & Najafian, 2 019). If possible, please use \code{\link[=query]{query()}} instead.
}
\examples{
data(googlenews)
S1 <- c("janitor", "statistician", "midwife", "bailiff", "auctioneer",
"photographer", "geologist", "shoemaker", "athlete", "cashier", "dancer",
"housekeeper", "accountant", "physicist", "gardener", "dentist", "weaver",
"blacksmith", "psychologist", "supervisor", "mathematician", "surveyor",
"tailor", "designer", "economist", "mechanic", "laborer", "postmaster",
"broker", "chemist", "librarian", "attendant", "clerical", "musician",
"porter", "scientist", "carpenter", "sailor", "instructor", "sheriff",
"pilot", "inspector", "mason", "baker", "administrator", "architect",
"collector", "operator", "surgeon", "driver", "painter", "conductor",
"nurse", "cook", "engineer", "retired", "sales", "lawyer", "clergy",
"physician", "farmer", "clerk", "manager", "guard", "artist", "smith",
"official", "police", "doctor", "professor", "student", "judge",
"teacher", "author", "secretary", "soldier")
A1 <- c("he", "son", "his", "him", "father", "man", "boy", "himself",
"male", "brother", "sons", "fathers", "men", "boys", "males", "brothers",
"uncle", "uncles", "nephew", "nephews")
B1 <- c("she", "daughter", "hers", "her", "mother", "woman", "girl",
"herself", "female", "sister", "daughters", "mothers", "women", "girls",
"females", "sisters", "aunt", "aunts", "niece", "nieces")
garg_f1 <- rnsb(googlenews, S1, A1, B1)
plot_bias(garg_f1)
}
\references{
Sweeney, C., & Najafian, M. (2019, July). \href{https://aclanthology.org/P19-1162/}{A transparent framework for evaluating unintended demographic bias in word embeddings.} In Proceedings of the 57th Annual Meeting of the Association for Computational Linguistics (pp. 1662-1667).
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sweater.R
\name{weat_exact}
\alias{weat_exact}
\alias{weat_resampling}
\title{Test of significance for WEAT}
\usage{
weat_exact(x)

weat_resampling(x, n_resampling = 9999)
}
\arguments{
\item{x}{an object from the \link{weat} function.}

\item{n_resampling}{an integer specifying the number of replicates used to estimate the exact test}
}
\value{
A list with class \code{"htest"}
}
\description{
This function conducts the test of significance for WEAT as described in Caliskan et al. (2017). The exact test (proposed in Caliskan et al.) takes an unreasonably long time, if the total number of words in S and T is larger than 10. The resampling test is an approximation of the exact test.
}
\examples{
# Reproduce the number in Caliskan et al. (2017) - Table 1, "Math vs. Arts"
data(glove_math)
S1 <- c("math", "algebra", "geometry", "calculus", "equations",
"computation", "numbers", "addition")
T1 <- c("poetry", "art", "dance", "literature", "novel", "symphony", "drama", "sculpture")
A1 <- c("male", "man", "boy", "brother", "he", "him", "his", "son")
B1 <- c("female", "woman", "girl", "sister", "she", "her", "hers", "daughter")
sw <- weat(glove_math, S1, T1, A1, B1)
weat_resampling(sw)
}
\references{
Caliskan, A., Bryson, J. J., & Narayanan, A. (2017). Semantics derived automatically from language corpora contain human-like biases. Science, 356(6334), 183-186. \doi{10.1126/science.aal4230}
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mac.R
\name{mac_es}
\alias{mac_es}
\title{Calculation of MAC Effect Size}
\usage{
mac_es(x)
}
\arguments{
\item{x}{an object from the function \link{mac}}
}
\value{
Mean of all cosine similarity values
}
\description{
This function calculates the mean of cosine distance values. If possible, please use \code{\link[=calculate_es]{calculate_es()}} instead.
}
\references{
Manzini, T., Lim, Y. C., Tsvetkov, Y., & Black, A. W. (2019). \href{https://arxiv.org/abs/1904.04047}{Black is to criminal as caucasian is to police: Detecting and removing multiclass bias in word embeddings.} arXiv preprint arXiv:1904.04047.
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R, R/s3.R
\name{query}
\alias{query}
\alias{print.sweater}
\title{A common interface for making query}
\usage{
query(
  w,
  S_words,
  T_words,
  A_words,
  B_words,
  method = "guess",
  verbose = FALSE,
  ...
)

\method{print}{sweater}(x, ...)
}
\arguments{
\item{w}{a numeric matrix of word embeddings (e.g. from \code{\link[rsparse:GloVe]{rsparse::GloVe()}})}

\item{S_words}{a character vector of the first set of target words. In an example of studying gender stereotype, it can include occupations such as programmer, engineer, scientists...}

\item{T_words}{a character vector of the second set of target words. In an example of studying gender stereotype, it can include occupations such as nurse, teacher, librarian...}

\item{A_words}{a character vector of the first set of attribute words. In an example of studying gender stereotype, it can include words such as man, male, he, his.}

\item{B_words}{a character vector of the second set of attribute words. In an example of studying gender stereotype, it can include words such as woman, female, she, her.}

\item{method}{string, the method to be used to make the query. Available options are: \code{weat}, \code{mac}, \code{nas}, \code{semaxis}, \code{rnsb}, \code{rnd}, \code{nas}, \code{ect} and \code{guess}. If "guess", the function selects one of the following methods based on your provided wordsets.
\itemize{
\item{S_words & A_words - }{"mac"}
\item{S_words, A_words & B_words - }{"rnd"}
\item{S_words, T_words, A_words & B_words - }{"weat"}
}}

\item{verbose}{logical, whether to display information}

\item{...}{additional parameters for the underlying function
\describe{
\item{\code{l}}{for "semaxis": an integer indicates the number of words to augment each word in A and B based on cosine , see An et al (2018). Default to 0 (no augmentation).}
\item{\code{levels}}{for "rnsb": levels of entries in a hierarchical dictionary that will be applied (see \code{\link[quanteda:dfm_lookup]{quanteda::dfm_lookup()}})}
}}

\item{x}{a sweater S3 object}
}
\value{
a sweater S3 object
}
\description{
This function makes a query based on the supplied parameters. The object can then be displayed by the S3 method \code{\link[=print.sweater]{print.sweater()}} and plotted by \code{\link[=plot.sweater]{plot.sweater()}}.
}
\examples{
data(googlenews)
S1 <- c("janitor", "statistician", "midwife", "bailiff", "auctioneer",
"photographer", "geologist", "shoemaker", "athlete", "cashier", "dancer",
"housekeeper", "accountant", "physicist", "gardener", "dentist", "weaver",
"blacksmith", "psychologist", "supervisor", "mathematician", "surveyor",
"tailor", "designer", "economist", "mechanic", "laborer", "postmaster",
"broker", "chemist", "librarian", "attendant", "clerical", "musician",
"porter", "scientist", "carpenter", "sailor", "instructor", "sheriff",
"pilot", "inspector", "mason", "baker", "administrator", "architect",
"collector", "operator", "surgeon", "driver", "painter", "conductor",
"nurse", "cook", "engineer", "retired", "sales", "lawyer", "clergy",
"physician", "farmer", "clerk", "manager", "guard", "artist", "smith",
"official", "police", "doctor", "professor", "student", "judge",
"teacher", "author", "secretary", "soldier")
A1 <- c("he", "son", "his", "him", "father", "man", "boy", "himself",
"male", "brother", "sons", "fathers", "men", "boys", "males", "brothers",
"uncle", "uncles", "nephew", "nephews")
B1 <- c("she", "daughter", "hers", "her", "mother", "woman", "girl",
"herself", "female", "sister", "daughters", "mothers", "women", "girls",
"females", "sisters", "aunt", "aunts", "niece", "nieces")
garg_f1 <- query(googlenews, S_words = S1, A_words = A1, B_words = B1)
garg_f1
plot(garg_f1)
}
\seealso{
\code{\link[=weat]{weat()}}, \code{\link[=mac]{mac()}}, \code{\link[=nas]{nas()}}, \code{\link[=semaxis]{semaxis()}}, \code{\link[=rnsb]{rnsb()}}, \code{\link[=rnd]{rnd()}}, \code{\link[=nas]{nas()}}, \code{\link[=ect]{ect()}}
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mac.R
\name{mac}
\alias{mac}
\title{Mean average cosine similarity}
\usage{
mac(w, S_words, A_words, verbose = FALSE)
}
\arguments{
\item{w}{a numeric matrix of word embeddings (e.g. from \code{\link[rsparse:GloVe]{rsparse::GloVe()}})}

\item{S_words}{a character vector of the first set of target words. In an example of studying gender stereotype, it can include occupations such as programmer, engineer, scientists...}

\item{A_words}{a character vector of the first set of attribute words. In an example of studying gender stereotype, it can include words such as man, male, he, his.}

\item{verbose}{logical, whether to display information}
}
\value{
A list with class \code{"mac"} containing the following components:
\describe{
\item{\code{$P}}{a vector of cosine similarity values for every word in S_words}
\item{\code{$S_words}}{the input S_words}
\item{\code{$A_words}}{the input A_words}
}
\code{\link{mac_es}} can be used to obtain the effect size of the test.
}
\description{
This function calculates the mean average cosine similarity (MAC) score proposed in Manzini et al (2019). If possible, please use \code{\link[=query]{query()}} instead.
}
\examples{
data(googlenews)
S1 <- c("janitor", "statistician", "midwife", "bailiff", "auctioneer",
"photographer", "geologist", "shoemaker", "athlete", "cashier", "dancer",
"housekeeper", "accountant", "physicist", "gardener", "dentist", "weaver",
"blacksmith", "psychologist", "supervisor", "mathematician", "surveyor",
"tailor", "designer", "economist", "mechanic", "laborer", "postmaster",
"broker", "chemist", "librarian", "attendant", "clerical", "musician",
"porter", "scientist", "carpenter", "sailor", "instructor", "sheriff",
"pilot", "inspector", "mason", "baker", "administrator", "architect",
"collector", "operator", "surgeon", "driver", "painter", "conductor",
"nurse", "cook", "engineer", "retired", "sales", "lawyer", "clergy",
"physician", "farmer", "clerk", "manager", "guard", "artist", "smith",
"official", "police", "doctor", "professor", "student", "judge", "teacher",
"author", "secretary", "soldier")
A1 <- c("he", "son", "his", "him", "father", "man", "boy", "himself",
"male", "brother", "sons", "fathers", "men", "boys", "males", "brothers",
"uncle", "uncles", "nephew", "nephews")
x <- mac(googlenews, S1, A1)
x$P
}
\references{
Manzini, T., Lim, Y. C., Tsvetkov, Y., & Black, A. W. (2019). \href{https://arxiv.org/abs/1904.04047}{Black is to criminal as caucasian is to police: Detecting and removing multiclass bias in word embeddings.} arXiv preprint arXiv:1904.04047.
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ect.R
\name{ect_es}
\alias{ect_es}
\title{Calculate the Spearman Coefficient of an ECT result}
\usage{
ect_es(x)
}
\arguments{
\item{x}{an ect object from the \code{\link[=ect]{ect()}} function.}
}
\value{
Spearman Coefficient
}
\description{
This functions calculates the Spearman Coefficient of an Embedding Coherence Test. The value ranges from -1 to +1 and a larger value indicates less bias. If possible, please use \code{\link[=calculate_es]{calculate_es()}} instead.
}
\references{
Dev, S., & Phillips, J. (2019, April). \href{https://proceedings.mlr.press/v89/dev19a.html}{Attenuating bias in word vectors.} In The 22nd International Conference on Artificial Intelligence and Statistics (pp. 879-887). PMLR.
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnd.R
\name{rnd_es}
\alias{rnd_es}
\title{Calculation of sum of all relative norm distances}
\usage{
rnd_es(x)
}
\arguments{
\item{x}{an object from the function \link{rnd}}
}
\value{
Sum of all relative norm distances
}
\description{
This function calculates the sum of all relative norm distances from the relative norm distance test. If possible, please use \code{\link[=calculate_es]{calculate_es()}} instead.
}
\references{
Garg, N., Schiebinger, L., Jurafsky, D., & Zou, J. (2018). Word embeddings quantify 100 years of gender and ethnic stereotypes. Proceedings of the National Academy of Sciences, 115(16), E3635-E3644. \doi{10.1073/pnas.1720347115}
}
\author{
Chung-hong Chan
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{read_word2vec}
\alias{read_word2vec}
\title{A helper function for reading word2vec format}
\usage{
read_word2vec(x)
}
\arguments{
\item{x}{path to your text file}
}
\value{
a dense matrix
}
\description{
This function reads word2vec text format and return a dense matrix that can be used by this package.
The file can have or have not the "verification line", i.e. the first line contains the dimensionality of the matrix. If the verification line exists, the function will check the returned matrix for correctness.
}
\author{
Chung-hong Chan
}
