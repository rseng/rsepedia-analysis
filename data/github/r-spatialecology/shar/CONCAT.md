---
title: '`shar`: An `R` package to analyze species-habitat associations using point pattern analysis'
tags:
  - R
  - open source software
  - spatial point pattern analysis
  - gamma-test
  - pattern reconstruction
  - torus-translation test
  - randomized-habitats procedure;
authors:
  - name: Maximilian H. K. Hesselbarth^[corresponding author]
    orcid: 0000-0003-1125-9918
    affiliation: "1, 2"
affiliations:
 - name: Department of Ecosystem Modelling, University of Goettingen, Buesgenweg 4, 37077, Goettingen
   index: 1
 - name: Department of Ecology and Evolutionary Biology, University of Michigan, 1105 N University Ave, Ann Arbor, Michigan 48109, USA
   index: 2
date: 2021-09-28
bibliography: paper.bib
---

# Summary

Studying species-habitat associations is one approach to reveal the importance of abiotic processes in shaping the spatial distribution of species.
Even though the `R` programming language offers many packages for spatial point pattern analysis, currently there is no comprehensive package specifically designed to analyze species-habitat associations.
The `shar` package builds on widely used `R` packages for spatial analyses and provides an easy and straightforward way to uncover species-habitat associations for discrete environmental data.

# Statement of need

Species-habitat associations are a result of certain species being specialized to certain environmental conditions [@Tilman1993] and typically result in a clustering of individuals at suitable habitats [@Harms2001; @Comita2007].
Thus, analyzing species-habitat associations can help to understand the importance of abiotic processes shaping the spatial distribution of species [@Garzon-Lopez2014].
However, since biotic processes (e.g., competition, limited dispersal) can also lead to a clustering of individuals, analyses of species-habitat associations need to control for potential biotic processes because they result in a violation of the independence assumption of similar statistical tests, such as the $\chi^2$ test, ordination methods, or canonical correspondence analysis [@Plotkin2000; @Harms2001].
Spatial autocorrelation of the environmental conditions could also violate the independence assumption of the previously mentioned statistical tests [@Harms2001].
For example, previous research has shown that violating the independence assumptions of the $\chi^2$ test resulted in more (i.e., possible false-positive) species-habitat associations than the more conservative methods that are provided by the `shar` package [@Plotkin2000; @Harms2001].

Ecologists use the spatial distribution of ecological objects to infer the underlying processes that shaped their distribution because spatial patterns can act as a "memory" of the processes that shaped them [@Velazquez2016].
A spatial point pattern contains the fully mapped locations (in terms of *x*,*y* coordinates) of all individual objects in a normally two-dimensional study area and assumes that the object locations can be described by discrete points [@Wiegand2014; @Velazquez2016].
For example, many studies use individual tree locations to infer the processes that determined their distribution [@Velazquez2016], and further examples include gopher mounds [@Klaas2000], gorilla nests [@Funwi-Gabga2012], "fairy circles" [@Getzin2015], or bacteria on leaves [@Esser2015a].

The `spatstat` package [@Baddeley2015] allows ecologists to access many methods of spatial point pattern analysis, such as summary functions and null model simulations, and to simulate heterogeneous point process models that can be used to show the importance of abiotic processes for continuous environmental data [@Getzin2008].
However, even though many ecological studies on species-habitat associations using discrete environmental data can be found in the literature [@John2007; @Garzon-Lopez2014; @Guo2016; @Yang2016; @Du2017; @Furniss2017], `spatstat` cannot be used to reveal such associations without larger programming efforts by the users.
The `inlabru` package provides an approach to  analyze the importance of abiotic processes, mostly for continuous environmental data, using Bayesian spatial modelling [@Bachl2019].
The `fgeo` package [@Lepore2019] allows to visualize and analyze forest diversity, including species-habitat associations.
But, the `fgeo` was designed to specifally handle ForestGEO data (<https://forestgeo.si.edu>) and furthermore includes only a subset of methods available to analyze species-habitat associations.

Thus, the `shar` package was developed to provide a comprehensive tool set to analyze species-habitat associations of spatial point patterns.
All methods in the `shar` package are designed for discrete environmental data and have the advantage of very few assumptions about the data characteristics.
In order to make the `shar` package as accessible for as many people as possible, it builds on two of the most commonly used `R` packages for spatial data, namely the `spatstat` and `raster` packages [@Hijmans2019].

# Methodological background

To analyze species-habitat associations, potential interdependence between the object locations and the environmental conditions needs to be broken by randomizing the data as a null model.
Within the field of spatial point pattern analysis, there are two related approaches to break potential dependencies [@Plotkin2000; @Harms2001].
Both require the spatial location of all objects, as well as discrete raster data for environmental conditions.

The first approach to simulate null model data is to randomize the environmental data, while keeping the object locations stable.
This can be achieved by shifting the raster data around a torus ("torus translation test") or using a random walk algorithm ("randomized-habitats procedure") [@Harms2001].
The second approach is to randomize the object locations, while keeping the environmental data stable.
This can be achieved by fitting point process models to the object locations ("gamma test") [@Plotkin2000] or using a simulated annealing approach ("pattern reconstruction") [@Kirkpatrick1983].

The two approaches differ in how they randomize the null model data, but both control for potential biotic processes by preserving the spatial structure of the data [@Plotkin2000; @Wiegand2014] and result in similar results.
Finally, species-habitat associations are present if species are found in certain habitats in the data more often than expected compared to the randomized null model data [@Plotkin2000; @Harms2001].
Given the characteristics of the method, a positive association to one habitat inevitably leads to a negative association to at least one of the other habitats (and vice versa) [@Yamada2006].

# How to use the package

Analyzing species-habitat associations is straightforward with the `shar` package.
Only two objects are needed to quantify species-habitat associations, namely a `spatstat` object that includes all object locations within the study area and a `raster` object with discrete habitat classes.

However, all methods require "fully mapped data" in the sense that NA cells of the environmental data are allowed only if simultaneously these areas cannot accommodate any locations of the point pattern (e.g., a water body within a forest area).
This needs to be reflected in the observation window of the point pattern.
For the torus translation method, no NA values are allowed at all.

To randomize the environmental data, either `translate_raster()` or `randomize_raster()` can be used.
For the former, the number of randomizations of the null model data is automatically determined by the number of rows and columns of the `raster` object.
For the later, the number of randomizations must be specified using the `n_random` argument.

```
torus_trans <- translate_raster(raster = landscape_discrete)

random_walk <- randomize_raster(raster = landscape_discrete, n_random = 39)
```

Alternatively, to randomize the object locations, either `fit_point_process()`, `reconstruct_pattern()`, or `reconstruct_pattern_marks()` can be used.
In all cases, the number of randomization must be specified using the `n_random` argument.
In order to preserve the spatial structure of the input as detailed as possible, several options are present to acknowledge for example if the input object locations are clustered or heterogeneously distributed in the study area.

```
gamma_test <- fit_point_process(pattern = species_a, n_random = 39,
                                process = "cluster")

reconstruction <- reconstruct_pattern(pattern = species_a, n_random = 39,
                                      method = "cluster")
```

Lastly, the input data and the randomized null model data are used to test if significant species-habitat associations are present.
The `results_habitat_association()` function automatically detects which part of the data was randomized and can be used identically with either of the used randomization approach.

```
results_habitat_association(pattern = species_a, raster = random_walk,
                            significance_level = 0.01)

> Input: randomized raster
> Quantile thresholds: negative < 0.005 || positive > 0.995
  habitat count lo hi significance
1       1    35 10 35         n.s.
2       2    44 19 53         n.s.
3       3    36 15 49         n.s.
4       4     4 15 58     negative
5       5    73 48 90         n.s.
```

The `shar` package also provides several utility and plotting functions such as a generic `plot()` function to plot the null model data, `calculate_energy()` to calculate the difference between the input object locations and the randomized null model data object locations, or `classify_habitats()` to classify continuous environmental data into discrete habitats. For all functions, please see the "Functions" article on the `shar` homepage (<https://r-spatialecology.github.io/shar>).

# Parallelization

One major drawback of the `shar` package is the computation time related to some of the randomization methods for null model data.
This is the case especially for pattern reconstruction, even though most point pattern analysis studies use less than 1000 null model randomizations [@Velazquez2016].
However, since the randomizations of the null model data are independent of each other, this could be parallized using available frameworks, such as the `future` [@Bengtsson2021] or
`parallel` [@RCoreTeam2021] package.
The `shar` package does not allow to run code in parallel internally to not limit users to a specific parallelization framework.
For a short example how to simulate null model data in parallel using the `future` package, please see the "Parallelization" article on the `shar` homepage (<https://r-spatialecology.github.io/shar>).
However, the presented approach could be used with any other parallelization framework as well.

# Acknowledgements

Support was provided by the German Research Association (DFG) Research Training Group 1644 "Scaling Problems in Statistics", grant number 152112243.
M.H.K.H. is thankful to Sebastian Hanss und Marco Sciaini for their help during the development of the `shar` package and Katrina Munsterman, Bridget Shayka and Samantha Iliff for comments on earlier drafts of the manuscript.
Thomas Etherington and Lionel Hertzog provided valuable feedback during the review process that improved the manuscript and the `shar` package.

# References

<!-- README.md is generated from README.Rmd. Please edit that file -->

# **shar** | **S**pecies **h**abitat **a**ssociations in **R** <img src="man/figures/logo.png" align="right" alt="" width="150" />

<!-- badges: start -->

| CI                                                                                                                                                                                   | Development                                                                                                                        | CRAN                                                                                                                    | License                                                                                                                                              |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| [![R-CMD-check](https://github.com/r-spatialecology/shar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/r-spatialecology/shar/actions/workflows/R-CMD-check.yaml) | [![Project Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)                       | [![CRAN status](https://www.r-pkg.org/badges/version/shar)](https://cran.r-project.org/package=shar)                    | [![JOSS](https://joss.theoj.org/papers/1b786c028a5425858cb0e5428bd9173b/status.svg)](https://joss.theoj.org/papers/1b786c028a5425858cb0e5428bd9173b) |
| [![codecov](https://codecov.io/gh/r-spatialecology/shar/branch/main/graph/badge.svg?token=XMo844ABs4)](https://codecov.io/gh/r-spatialecology/shar)                                  | [![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable) | [![CRAN logs](http://cranlogs.r-pkg.org/badges/grand-total/shar)](http://cran.rstudio.com/web/packages/shar/index.html) | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)                                      |

<!-- badges: end -->

**S**pecies-**h**abitat **a**ssociations in **R** is a `R` package to
analyze species-habitat associations. Therefore, information about the
location of the species is needed (as a point pattern) and about the
environmental conditions (as a raster map). In order to analyse the data
for significant habitat associations either the location data or the
environmental data is randomized n-times. Then, counts within the
habitats are compared between the randomized data and the observed data.
Positive or negative associations are present if the observed counts is
higher or lower than the randomized counts (using quantile thresholds).
Methods are mainly described in Plotkin et al. (2000), Harms et
al. (2001) and Wiegand & Moloney (2014). **shar** is mainly based on
the [`spatstat`](http://spatstat.org) (Baddeley et al. 2015) and
[`raster`](https://rspatial.org/raster/) (Hijmans 2017) package.

#### Citation

The **shar** package is part of our academic work. To cite the package
or acknowledge its use in publications, please cite the following paper.

> Hesselbarth, M.H.K., (2021). shar: A R package to analyze
> species-habitat associations using point pattern analysis. Journal of
> Open Source Software, 6(67), 3811.
> <https://doi.org/10.21105/joss.03811>

The get a BibTex entry, please use `citation("shar")`.

## Installation

You can install the released version of **shar** from
[CRAN](https://cran.r-project.org/web/packages/shar/index.html) with:

``` r
install.packages("shar")
```

And the development version from
[GitHub](https://github.com/r-spatialecology/shar) with:

``` r
install.packages("remotes")

remotes::install_github("r-spatialecology/shar")
```

This also automatically installs all non-base `R` package dependencies,
namely the following packages: `classInt`, `raster`, `spatstat.core`,
`spatstat.geom`.

## How to use shar

``` r
library(shar)
library(spatstat)
library(raster)

set.seed(42)
```

**shar** comes with build-in example data sets. `species_a` and
`species_b` are exemplary location of species, e.g. trees, as
`ppp`-objects from the `spatstat` package. `landscape` contains
exemplary continuous environmental data. However, all methods depend on
discrete data. Therefore we need to classify the data first. However,
all methods require “fully mapped data” in a sense that NA cells of the
environmental data are allowed only if simultaneously these areas cannot
accommodate any locations of the point pattern (e.g., a water body
within a forest area). This needs to be reflected in the observation
window of the point pattern. For the torus translation method, no NA
values are allowed at all.

``` r
landscape_classified <- classify_habitats(raster = landscape, n = 5, style = "fisher")
```

There are two possibilities to randomize the environmental data, both
described in Harms et al. (2001). The first shifts the habitat map in
all 4 cardinal directions around a torus. The second one assigns the
habitat values to an empty map using a random walk algorithm. Both
functions return a list with randomized rasters and the observed one.
For more information on the methods, please click
[here](https://r-spatialecology.github.io/shar/articles/articles/background.html).

``` r
torus_trans <- translate_raster(raster = landscape_classified)

random_walk <- randomize_raster(raster = landscape_classified, n_random = 99)
```

To plot the randomized raster, you can use the plot function and specify
the number of raster as as well as the color palette used for the
discrete environmental data.

``` r
col_palette <- c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")

plot(torus_trans, n = 3, col = col_palette)
```

To randomize the point pattern, either use the Gamma test described by
Plotkin et al. (2000) or pattern reconstruction (Kirkpatrick et
al. 1983; Tscheschel & Stoyan 2006).

``` r
gamma_test <- fit_point_process(pattern = species_b, process = "cluster", n_random = 99)

# (this can takes some time)
reconstruction <- reconstruct_pattern(pattern = species_b, n_random = 99, e_threshold = 0.05)
```

Of course, there are several utility functions. For example, you can
plot the summary function of the observed pattern and the simulation
envelopes of randomized patterns (`what = "sf"`) or some randomized and
the observed pattern (`what = "pp"`) using the plot function.

``` r
plot(reconstruction, what = "pp")
```

<img src="man/figures/README-plot-random_pattern-1.png" width="100%" height="100%" style="display: block; margin: auto;" />

Another utility functions allows to calculate the differences between
the observed pattern and the randomized patterns (also called energy
using summary functions).

``` r
calculate_energy(reconstruction, return_mean = TRUE)
## [1] 0.04908566
```

The data was created that `species_a` has a negative association to
habitat 4 and `species_b` has a positive association to habitat 5, which
is reflected in the results.

Given the characteristics of the method, a positive association to one
habitat inevitably leads to a negative association to at least one of
the other habitats (and vice versa; Yamada et al. 2006). For example, a
high amount of individual points in the positively associated habitat
simultaneously mean that less individual points can be present in the
other habitats.

Furthermore, please be aware that due to the randomization of the null
model data, results might slightly differ between different
randomization approaches (e.g., `fit_point_process()`
vs. `translate_raster()`) and even for repetitions of the same
approach. Thus, the exact `lo` and `hi` thresholds might be slightly
different when re-running the examples. However, the counts of the
observed data should be identical, and general results and trends should
be similar.

``` r
significance_level <- 0.01

results_habitat_association(pattern = species_a, raster = torus_trans, significance_level = significance_level)
## > Input: randomized raster
## > Quantile thresholds: negative < 0.005 || positive > 0.995
##   habitat breaks count lo hi significance
## 1       1     NA    35 10 35         n.s.
## 2       2     NA    44 19 53         n.s.
## 3       3     NA    36 15 49         n.s.
## 4       4     NA     4 15 58     negative
## 5       5     NA    73 48 90         n.s.

results_habitat_association(pattern = reconstruction, raster = landscape_classified, significance_level = significance_level)
## > Input: randomized pattern
## > Quantile thresholds: negative < 0.005 || positive > 0.995
##   habitat breaks count    lo    hi significance
## 1       1     NA     6 21.96 49.02     negative
## 2       2     NA    18 32.47 64.51     negative
## 3       3     NA    18 26.98 56.10     negative
## 4       4     NA    21 17.98 40.00         n.s.
## 5       5     NA   129 24.96 52.02     positive
```

## Contributing and Code of Conduct

Contributions to **shar** are highly welcomed and appreciated. This
includes any form of feedback, bug reports, feature
requests/suggestions, or general questions about the usage.

Please note that the **shar** package is released with a [Contributor
Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project,
you agree to abide by its terms.

To see how to contribute to this project, please see the [Contributing
guidelines](CONTRIBUTING.md).

### References

Baddeley, A., Rubak, E., Turner, R., 2015. Spatial point patterns:
Methodology and applications with R. Chapman and Hall/CRC Press, London.
ISBN 978-1-4822-1020-0

Harms, K.E., Condit, R., Hubbell, S.P., Foster, R.B., 2001. Habitat
associations of trees and shrubs in a 50-ha neotropical forest plot.
Journal of Ecology 89, 947–959.
<https://doi.org/10.1111/j.1365-2745.2001.00615.x>

Hijmans, R.J., 2019. raster: Geographic data analysis and modeling. R
package version 2.9-5. <https://cran.r-project.org/package=raster>.

Kirkpatrick, S., Gelatt, C.D.Jr., Vecchi, M.P., 1983. Optimization by
simulated annealing. Science 220, 671–680.
<https://doi.org/10.1126/science.220.4598.671>

Plotkin, J.B., Potts, M.D., Leslie, N., Manokaran, N., LaFrankie, J.V.,
Ashton, P.S., 2000. Species-area curves, spatial aggregation, and
habitat specialization in tropical forests. Journal of Theoretical
Biology 207, 81–99. <https://doi.org/10.1006/jtbi.2000.2158>

Tscheschel, A., Stoyan, D., 2006. Statistical reconstruction of random
point patterns. Computational Statistics and Data Analysis 51, 859–871.
<https://doi.org/10.1016/j.csda.2005.09.007>

Wiegand, T., Moloney, K.A., 2014. Handbook of spatial point-pattern
analysis in ecology. Chapman and Hall/CRC Press, Boca Raton. ISBN
978-1-4200-8254-8

Yamada, T., Tomita, A., Itoh, A., Yamakura, T., Ohkubo, T., Kanzaki, M.,
Tan, S., Ashton, P.S., 2006. Habitat associations of Sterculiaceae trees
in a Bornean rain forest plot. Journal of Vegetation Science 17,
559–566. <https://doi.org/10.1111/j.1654-1103.2006.tb02479.x>
# shar 1.3.2
* Improvements
  * Improvement of `classify_habitats()` to be more variable
  * Adding breaks argument to `results_habitat_association()`
  * Adding `classint_to_vector()` helper function

# shar 1.3.1
* Improvements
  * Bugfix in `plot.rd_pat()` and `plot.rd_mar()`
  * Adding internal `sample_randomized()` function used during plotting

# shar 1.3
* Improvements
  * Better documentation
  * Combine `reconstruct_pattern_homo()`, `reconstruct_pattern_cluster()`, and `reconstruct_pattern_hetero()` to `reconstruct_pattern()`
  * Replaced `plot_randomized_*()` function with generic `plot()` methods
  * Adding warnings and errors if `NA` values are present  
* New functionality
  * `list_to_randomized()` function
  * Parallelization article
* Adding JOSS paper as reference

# shar 1.2.1
* Improvements
  * `reconstruct_pattern_homo()` has arguments to specify number of points and window
  * `reconstruct_pattern_marks()` allows to have different number of points for `pattern` and `marked_pattern` argument

# shar 1.2
* Improvements
  * Include new `spatstat` package structure
  * Use GPL3 License

# shar 1.1.1
* Improvements
   * Add logo
   * Update to MIT License
   * renamed `master` to `main` branch

# shar 1.1
* Improvements
  * Use `energy_df` to get energy for printing if available
  * Updated tests
  * More stable progress printing

# shar 1.0.1
* Improvements
  * No calculation of energy for printing (too slow)

# shar 1.0
* Improvements
  * Printing methods for most objects
  * Possibility to specify intervals of r for all reconstruction functions
* New functionality
  * `plot_energy()` to plot energy over iterations for reconstructed patterns
  * `reconstruct_pattern_hetero()` allows to reconstruct heterogeneous patterns
* Renameing/Structure
  * `reconstruct_pattern()` was split to three functions: `reconstruct_pattern_homo()`, `reconstruct_pattern_hetero()`, `reconstruct_pattern_cluster()`,
  * `reconstruct_marks()` is now called `reconstruct_pattern_marks()`

# shar 0.5
* Improvements
  * Annealing probability can be specified for reconstruction

# shar 0.4
* Improvements
  * Easier user experience because classes are used to specify provided input
  * `results_habitat_associations()` checks if extent of inputs is identical
  * `reconstruct_marks()` and `calculate_energy()` use now weights for the summary functions
* Bugfixes
  * Bug in `calculate_energy()` for reconstructed marks

# shar 0.3.1
* Improvements
  * Better structure of vignettes
  * Adding CONTRIBUTING.md
  * Trying to fix some failing tests for older R versions
 * New functionality
  * Allowing to translate raster only in n steps

# shar 0.3
* Improvements
  * `plot_randomized_pattern()` now uses envelopes to plot randomized summary functions
  * `plot_randomized_pattern()` includes a quantum bar
  * `plot_randomized_pattern()` now can return plots after each other (Press <Enter>)
  * `calculate_engery()` can also calculate the energy for marked reconstructions
  * Improved warning messages
* Bugfixes
  * Explicitly C++11 compiler
* New functionality
  * `plot_randomized_pattern()` to plot randomized rasters

# shar 0.2.1
* Improvements
  * minor speed improvment for `reconstruct_pattern()`, `reconstruct_marks()` and `calculate_energy()`
    * The starting pattern is now identical for all n_random and only the relocation process differs between randomizations
    * All summary functions are only calculated for 250 steps from 0 to `rmax`
* Bugfixes
* New functionality
  * `rcpp_sample()` as a faster Rcpp implementation of `sample()`

# shar 0.2
* Improvements
  * Replaced `cat()` with `message()` for all printing to console
  * All defaults set  to `n_random = 1`
  * `comp_fast` argument equals TRUE if number of points exceed threshold
  * `reconstruct_pattern()` stops if energy did not decrease for n iterations
  * `reconstruct_marks()` stops if energy did not decrease for n iterations
  * `plot_randomized_pattern()` can also plot point patterns
* Bugfixes
  * Bug in `fit_point_process()` that more points as present could be removed from simulated pattern
  * Bug in `reconstruct_pattern()` that more points as present could be removed from simulated pattern

# shar 0.1
* First submission to CRAN
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
reported to the community leaders responsible for enforcement at mhk.hesselbarth@gmail.com. 
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
available at <https://www.contributor-covenant.org/version/2/0/code_of_conduct.html>.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
<https://www.contributor-covenant.org/faq>. Translations are available at <https://www.contributor-covenant.org/translations>.
# Contributing to **shar**

This outlines how to propose a change to **shar**. There are several ways you can contribute to this project. If you want to know more about why and how to contribute to open source projects like this one, see this [Open Source Guide](https://opensource.guide/how-to-contribute/).

#### Code of Conduct

Please note that the **shar** project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project you agree to abide by its terms.

## How to contribute

### Ask a question :interrobang:

Browse the [documentation](https://r-spatialecology.github.io/shar/) to see if you can find a solution. Still stuck? Open an [issue on GitHub](https://github.com/r-spatialecology/shar/issues) on GitHub. We'll try to do our best to address it, as questions often lead to better documentation or the discovery of bugs.

Want to ask a question in private? Contact the package maintainer by [mhk.hesselbarth\<at\>gmail.com](mailto:mhk.hesselbarth@gmail.com).

Please try to include a reproducible example using for example the [`reprex`](https://reprex.tidyverse.org) package.

### Propose an idea :bulb:

Take a look at the [documentation](https://r-spatialecology.github.io/shar/) and [issue on GitHub](https://github.com/r-spatialecology/shar/issues) list to see if it isn't included or suggested yet. If not, please open a new issue!

While we can't promise to implement your idea, it helps to:

* Explain in detail how it would work.
* Keep the scope as narrow as possible.

### Report a bug :bug:

Report it as an [issue on GitHub](https://github.com/r-spatialecology/shar/issues) so we can fix it. A good bug report makes it easier for us to do so, so please include:

* The content of `utils::sessionInfo()`.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug

Again, please try to include a reproducible example using for example the [`reprex`](https://reprex.tidyverse.org) package.

### Improve the documentation :book:

Good documentation makes all the difference, so your help to improve it is very welcome!

We use [roxygen2](https://cran.r-project.org/package=roxygen2), with [Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), for documentation.

### Pull request process :arrow_up_down:

We try to follow the [GitHub flow](https://guides.github.com/introduction/flow/) for development.

1. Fork [the repo](https://github.com/r-spatialecology/shar) and clone it to your computer. To learn more about this process, see [this guide](https://guides.github.com/activities/forking/). Don't forget to pull all new changes before starting to work!

2. Open the RStudio project file (`.Rproj`) and install all development dependencies with (e.g., using `devtools::install_dev_deps()`). Make sure the package passes R CMD check by running `devtools::check()`.

3. Create a new Git branch and use a name that briefly describes the proposed changes.

4. Make your changes:
    * Write your code.
    * Test your code (bonus points for adding unit tests using the [`testthat`](https://testthat.r-lib.org) package).
    * Document your code (see function documentation above).
    * Check your code with `devtools::check()` and aim for 0 errors, warnings and notes.
5. Commit and push your changes.
6. Submit a [pull request](https://guides.github.com/activities/forking/#making-a-pull-request).

New code should follow the tidyverse [style guide](https://style.tidyverse.org).

### References

This CONTRIBUTING.md is adapted from [here](https://gist.github.com/peterdesmet/e90a1b0dc17af6c12daf6e8b2f044e7c) and [here](https://github.com/r-lib/usethis/blob/main/inst/templates/tidy-contributing.md).
For details changes, please see NEWS.md.

## shar 1.3.2
Improvements of existing functions

## shar 1.3.1
Bug fixes

## shar 1.3
Improvements of general package structure

## shar 1.2.1
Improvement of existing functions

## shar 1.2
Update spatstat dependencies

## shar 1.1.1
Minor improvments and new license

## shar 1.1
Improvments of existing functions

## shar 1.0.1
Some minor improvments to existing functions

## shar 1.0
New functionality and renameing/splitting of some functions

## shar 0.5
Improvments of existing functions

## shar 0.4
Improvments of general package structure

## shar 0.3.1
Some minor bugfixes

## shar 0.3
Some bugfixes and improvements of existing functions as well as new functions

## shar 0.2 
Some bugfixes and improvements of existing functions

## Review CRAN submission
1. Thanks. Please omit the redundant "in R". 

* Done as suggested
  
2. Is there some reference about the method you can add in the Description field in the form Authors (year) <doi:.....>? 

* Added "Methods are mainly based on Plotkin et al. (2000) <doi:10.1006/jtbi.2000.2158> and Harms et al. (2001) <doi:10.1111/j.1365-2745.2001.00615.x>." to the description field

3. Thanks, we see: 
  running tests for arch 'i386' ... [320s] OK 
  running tests for arch 'x64' ... [356s] OK
which is together alreadyv more than 10 min which is the CRAN threshold for a package check. Can you pls simplify the test cases or run less important tests only conditionally if some env var is set that you only define on your machine?

* The test run faster now (checked on win-builder.r-project)
  running tests for arch 'i386' ... [152s] OK
  running tests for arch 'x64' ... [164s] OK

* Renamed package from `SHAR` to `shar`

## Test environments
* Windows 10, R 3.5.1
* macOS Mojave, R 3.5.1
* https://win-builder.r-project.org (devel and release)

## R CMD check results
0 errors | 0 warnings | 0 note

## Reverse dependencies
There are currently no reverse dependencies.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE, 
  fig.path = "man/figures/README-"
)
```

# **shar** | **S**pecies **h**abitat **a**ssociations in **R** <img src="man/figures/logo.png" align="right" alt="" width="150" />

<!-- badges: start -->

| CI | Development | CRAN | License |
|----|-------------|------|---------|
| [![R-CMD-check](https://github.com/r-spatialecology/shar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/r-spatialecology/shar/actions/workflows/R-CMD-check.yaml) | [![Project Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) | [![CRAN status](https://www.r-pkg.org/badges/version/shar)](https://cran.r-project.org/package=shar) |[![JOSS](https://joss.theoj.org/papers/1b786c028a5425858cb0e5428bd9173b/status.svg)](https://joss.theoj.org/papers/1b786c028a5425858cb0e5428bd9173b) |
| [![codecov](https://codecov.io/gh/r-spatialecology/shar/branch/main/graph/badge.svg?token=XMo844ABs4)](https://codecov.io/gh/r-spatialecology/shar) | [![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable) | [![CRAN logs](http://cranlogs.r-pkg.org/badges/grand-total/shar)](http://cran.rstudio.com/web/packages/shar/index.html) | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |

<!-- badges: end -->

**S**pecies-**h**abitat **a**ssociations in **R** is a `R` package to analyze species-habitat associations. Therefore, information about the location of the species is needed (as a point pattern) and about the environmental conditions (as a raster map). In order to analyse the data for significant habitat associations either the location data or the environmental data is randomized n-times. Then, counts within the habitats are compared between the randomized data and the observed data. Positive or negative associations are present if the observed counts is higher or lower than the randomized counts (using quantile thresholds). Methods are mainly described in Plotkin et al. (2000), Harms et al. (2001) and Wiegand & Moloney (2014). **shar** is mainly based on the [`spatstat`](http://spatstat.org) (Baddeley et al. 2015) and [`raster`](https://rspatial.org/raster/) (Hijmans 2017) package.

#### Citation

The **shar** package is part of our academic work. To cite the package or acknowledge its use in publications, please cite the following paper.

> Hesselbarth, M.H.K., (2021). shar: A R package to analyze species-habitat associations using point pattern analysis. Journal of Open Source Software, 6(67), 3811. https://doi.org/10.21105/joss.03811

The get a BibTex entry, please use `citation("shar")`.

## Installation

You can install the released version of **shar** from [CRAN](https://cran.r-project.org/web/packages/shar/index.html) with:

```{r install-CRAN, eval = FALSE}
install.packages("shar")
```

And the development version from [GitHub](https://github.com/r-spatialecology/shar) with:

```{r install-github, eval = FALSE}
install.packages("remotes")

remotes::install_github("r-spatialecology/shar")
```

This also automatically installs all non-base `R` package dependencies, namely the following packages: `classInt`, `raster`, `spatstat.core`, `spatstat.geom`.

## How to use shar

```{r import-libs, message = FALSE, warning = FALSE}
library(shar)
library(spatstat)
library(raster)

set.seed(42)
```

**shar** comes with build-in example data sets. `species_a` and `species_b` are exemplary location of species, e.g. trees, as `ppp`-objects from the `spatstat` package. `landscape` contains exemplary continuous environmental data. However, all methods depend on discrete data. Therefore we need to classify the data first. 
However, all methods require "fully mapped data" in a sense that NA cells of the environmental data are allowed only if simultaneously these areas cannot accommodate any locations of the point pattern (e.g., a water body within a forest area). This needs to be reflected in the observation window of the point pattern. For the torus translation method, no NA values are allowed at all.

```{r environmental-data}
landscape_classified <- classify_habitats(raster = landscape, n = 5, style = "fisher")
```

There are two possibilities to randomize the environmental data, both described in Harms et al. (2001). The first shifts the habitat map in all 4 cardinal directions around a torus. The second one assigns the habitat values to an empty map using a random walk algorithm. Both functions return a list with randomized rasters and the observed one. For more information on the methods, please click [here](https://r-spatialecology.github.io/shar/articles/articles/background.html).

```{r habitat_random, eval = FALSE}
torus_trans <- translate_raster(raster = landscape_classified)

random_walk <- randomize_raster(raster = landscape_classified, n_random = 99)
```

To plot the randomized raster, you can use the plot function and specify the number of raster as as well as the color palette used for the discrete environmental data.

```{r plot_habitat-random, eval = FALSE}
col_palette <- c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")

plot(torus_trans, n = 3, col = col_palette)
```

To randomize the point pattern, either use the Gamma test described by Plotkin et al. (2000) or pattern reconstruction (Kirkpatrick et al. 1983; Tscheschel & Stoyan 2006). 

```{r pattern-random, eval = FALSE}
gamma_test <- fit_point_process(pattern = species_b, process = "cluster", n_random = 99)

# (this can takes some time)
reconstruction <- reconstruct_pattern(pattern = species_b, n_random = 99, e_threshold = 0.05)
```

Of course, there are several utility functions. For example, you can plot the summary function of the observed pattern and the simulation envelopes of randomized patterns (`what = "sf"`) or some randomized and the observed pattern (`what = "pp"`) using the plot function. 

```{r plot-random_pattern, fig.align = "center", out.height = "100%", out.width = "100%", message = FALSE}
plot(reconstruction, what = "pp")
```

Another utility functions allows to calculate the differences between the observed pattern and the randomized patterns (also called energy using summary functions). 

```{r calculate-energy, message = FALSE}
calculate_energy(reconstruction, return_mean = TRUE)
```

The data was created that `species_a` has a negative association to habitat 4 and `species_b` has a positive association to habitat 5, which is reflected in the results. 

Given the characteristics of the method, a positive association to one habitat inevitably leads to a negative association to at least one of the other habitats (and vice versa; Yamada et al. 2006). For example, a high amount of individual points in the positively associated habitat simultaneously mean that less individual points can be present in the other habitats.

Furthermore, please be aware that due to the randomization of the null model data, results might slightly differ between different randomization approaches (e.g., `fit_point_process()` vs. `translate_raster()`) and even for repetitions of the same approach. Thus, the exact `lo` and `hi` thresholds might be slightly different when re-running the examples. However, the counts of the observed data should be identical, and general results and trends should be similar.

```{r results}
significance_level <- 0.01

results_habitat_association(pattern = species_a, raster = torus_trans, significance_level = significance_level)

results_habitat_association(pattern = reconstruction, raster = landscape_classified, significance_level = significance_level)
```

## Contributing and Code of Conduct

Contributions to **shar** are highly welcomed and appreciated. This includes any form of feedback, bug reports, feature requests/suggestions, or general questions about the usage.

Please note that the **shar** package is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.

To see how to contribute to this project, please see the [Contributing guidelines](CONTRIBUTING.md).

### References 

Baddeley, A., Rubak, E., Turner, R., 2015. Spatial point patterns: Methodology and applications with R. Chapman and Hall/CRC Press, London. ISBN 978-1-4822-1020-0

Harms, K.E., Condit, R., Hubbell, S.P., Foster, R.B., 2001. Habitat associations of trees and shrubs in a 50-ha neotropical forest plot. Journal of Ecology 89, 947–959. <https://doi.org/10.1111/j.1365-2745.2001.00615.x>

Hijmans, R.J., 2019. raster: Geographic data analysis and modeling. R package version 2.9-5. <https://cran.r-project.org/package=raster>.

Kirkpatrick, S., Gelatt, C.D.Jr., Vecchi, M.P., 1983. Optimization by simulated annealing. Science 220, 671–680. <https://doi.org/10.1126/science.220.4598.671>

Plotkin, J.B., Potts, M.D., Leslie, N., Manokaran, N., LaFrankie, J.V., Ashton, P.S., 2000. Species-area curves, spatial aggregation, and habitat specialization in tropical forests. Journal of Theoretical Biology 207, 81–99. <https://doi.org/10.1006/jtbi.2000.2158>

Tscheschel, A., Stoyan, D., 2006. Statistical reconstruction of random point patterns. Computational Statistics and Data Analysis 51, 859–871. <https://doi.org/10.1016/j.csda.2005.09.007>

Wiegand, T., Moloney, K.A., 2014. Handbook of spatial point-pattern analysis in ecology. Chapman and Hall/CRC Press, Boca Raton. ISBN 978-1-4200-8254-8

Yamada, T., Tomita, A., Itoh, A., Yamakura, T., Ohkubo, T., Kanzaki, M., Tan, S., Ashton, P.S., 2006. Habitat associations of Sterculiaceae trees in a Bornean rain forest plot. Journal of Vegetation Science 17, 559–566. <https://doi.org/10.1111/j.1654-1103.2006.tb02479.x>
---
title: "Get started"
author: "Maximilian H.K. Hesselbarth"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Get started}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r import-libs, message = FALSE, warning = FALSE}
library(shar)
library(spatstat)
library(raster)

set.seed(42)
```

## Design

The core of **shar** are functions to to simulate null model data by randomizing either the environmental data (i.e. raster data) or the locations of species (i.e. point pattern data). The null data is then used to analyse if significant species-habitat associations are present. Additionally, functions to visualize and analyse the results are available as well as some utility functions. The methods are mainly described in Harms et al. (2001), Plotkin et al. (2000) and Wiegand & Moloney (2014). The methods are not necessary complementary, but are rather different approaches for the same result.

## Preprocessing of input data

All functions are designed for discrete habitat classes. Following, in case only continuous data is available, this has to be classified to discrete classes. `classify_habitats` provides several ways to classify the data such as the Fisher-Jenks algorithm (Fisher 1958, Jenks & Caspall 1971)

```{r classify_habitats}
landscape_discrete <- classify_habitats(raster = landscape, n = 5, style = "fisher")
```

## Randomize environmental data

There are two functions to randomize the environmental data: `translate_raster()` and `randomize_raster()`. The first function is a torus translation of the raster, shifting the habitat map in all four cardinal directions. This is only possible for rectangular observation areas and results in `n_random <- (raster::nrow(landscape) + 1) * (raster::ncol(landscape) + 1)  - 4` randomized rasters. The other function randomizes the environmental data using a random-walk algorithm. Here, the number of randomized rasters can be specified using the `n_random` argument.

However, all methods require "fully mapped data" in a sense that NA cells of the environmental data are allowed only if simultaneously these areas cannot accommodate any locations of the point pattern (e.g., a water body within a forest area). This needs to be reflected in the observation window of the point pattern. For the torus translation method, no NA values are allowed at all.

```{r randomize_raster, eval = FALSE}
torus_trans <- translate_raster(raster = landscape_discrete)

random_walk <- randomize_raster(raster = landscape_discrete, n_random = 99)
```

## Randomize location data

To randomize the location data (i.e. the point pattern) either `fit_point_process()` or `reconstruct_pattern()` are available. The first fits either a Poisson process or a cluster process to the data. The difference to solutions from the `spatstat` package is that the number of points is always identical. The second functions reconstructs the spatial characteristics of the data using pattern reconstruction (Tscheschel & Stoyan 2006). This is advantageous for point patterns not describable by simple point process models. For both function, the number of patterns can be specified by the `n_random` argument. 

```{r randomize_pp, eval = FALSE}
gamma_test <- fit_point_process(pattern = species_b, process = "cluster", n_random = 99)

# (this can takes some time)
reconstruction <- reconstruct_pattern(pattern = species_b, n_random = 99, 
                                      e_threshold = 0.05, method = "cluster")
``` 

## Analyse results

The most important function to analyse results is `results_habitat_association()`. This function compares the observed data to the null model data and by that is able to show significant species-habitat associations. The functions work for both, randomized environmental data or randomized location data. 

Please be aware that due to the randomization of the null model data, results might slightly differ between different randomization approaches (e.g., `fit_point_process()` vs. `translate_raster()`) and even for repetitions of the same approach. However, the counts of the observed data should be identical, and general results and trends should be similar.

```{r results}
results_habitat_association(pattern = species_a, raster = random_walk)

results_habitat_association(pattern = reconstruction, raster = landscape_discrete)
```

There is also the possibility to visualize the randomized data using the `plot()` function.

```{r plotting, fig.align = "center", out.height = "100%", out.width = "100%", message = FALSE}
plot(random_walk)

plot(reconstruction, ask = FALSE)
```

For the randomized point pattern data, it is also possible to show the "difference" in terms of the energy (Tscheschel & Stoyan 2006) between the patterns.

```{r energy, message = FALSE}
calculate_energy(pattern = gamma_test, return_mean = TRUE)

calculate_energy(pattern = reconstruction, return_mean = TRUE)
```

### References

Fisher, W.D., 1958. On grouping for maximum homogeneity. Journal of the American Statistical Association 53, 789–798. <https://doi.org/10.1080/01621459.1958.10501479>

Harms, K.E., Condit, R., Hubbell, S.P., Foster, R.B., 2001. Habitat associations of trees and shrubs in a 50-ha neotropical forest plot. Journal of Ecology 89, 947–959. <https://doi.org/10.1111/j.1365-2745.2001.00615.x>

Jenks, G.F., Caspall, F.C., 1971. Error in choroplethic maps: Definition, measurement, reduction. Annals of the Association of American Geographers 61, 217–244. <https://doi.org/10.1111/j.1467-8306.1971.tb00779.x>

Plotkin, J.B., Potts, M.D., Leslie, N., Manokaran, N., LaFrankie, J.V., Ashton, P.S., 2000. Species-area curves, spatial aggregation, and habitat specialization in tropical forests. Journal of Theoretical Biology 207, 81–99. <https://doi.org/10.1006/jtbi.2000.2158>

Tscheschel, A., Stoyan, D., 2006. Statistical reconstruction of random point patterns. Computational Statistics and Data Analysis 51, 859–871. <https://doi.org/10.1016/j.csda.2005.09.007>

Wiegand, T., Moloney, K.A., 2014. Handbook of spatial point-pattern analysis in ecology. Chapman and Hall/CRC Press, Boca Raton. ISBN 978-1-4200-8254-8
---
title: "Parallelization"
author: "Maximilian H.K. Hesselbarth"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Parallelization}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The `shar` packages has no build-in parallelization, however, there are many `R` frameworks that allow to run code in parallel or on high-performance clusters (see e.g., [`future`]( https://CRAN.R-project.org/package=future), [`clustermq`]( https://CRAN.R-project.org/package=clustermq) or [`rslurm`]( https://CRAN.R-project.org/package=rslurm)). Thus, `shar` provides utility functions that facilitate its usage together with any parallelization framework. 

The following examples illustrates how to use `future` package to randomize patterns using `fit_point_process` in parallel using all available cores on a local machine. Similarly, the core idea of the following code could be used to run `shar` on a high performance cluster.

First, we need to load all required packages. This includes `future` and `future.apply`.

```{r load-packages, message = FALSE, warning = FALSE}
library(shar)
library(spatstat)
library(raster)

library(future)
library(future.apply)
```

The `future` packages allows to run code in parallel using only a few lines of code. By setting the `future` plan to `multisession`, the package automatically resolves all following `futures` in parallel.

Importantly with this approach, you need one randomization per core (`n_random = 1`) and set `simplify = TRUE` to return the point pattern only. This results in a list of randomized point patterns.

```{r parallel}
future::plan(multisession)

fitted_list <- future.apply::future_lapply(X = 1:39, FUN = function(i) {
   shar::fit_point_process(pattern = shar::species_b, n_random = 1,
                          return_input = FALSE, simplify = TRUE, verbose = FALSE)
}, future.seed = 42)
```

Next, you can use the `list_to_randomized()` function to convert this list of randomized pattern to a `rd_pat` object that will work will all other functions of the `shar` package.

```{r convert-list}
fitted_rd <- shar::list_to_randomized(list = fitted_list, observed = shar::species_b)
```

Lastly, the created objects can be used to analyse if species-habitat associations are present as usual.

```{r results}
landscape_classified <- shar::classify_habitats(raster = landscape, n = 5, style = "fisher")

results_habitat_association(pattern = fitted_rd, raster = landscape_classified)
```

Of course, this idea can be used to randomize the raster data as well. Furthermore, any other parallelization framework could be used.
---
title: "Analysing the climatic niche of Cormus domestica"
author: "Zeke Marshall"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Analysing the climatic niche of Cormus domestica}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
editor_options: 
  chunk_output_type: console
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Outline

This vignette demonstrates how to use `shar` to analyse species occurrence data obtained from the Global Biodiversity Information Facility [(GBIF)](https://www.gbif.org) and environmental raster data obtained from the Climate Research Unit [(CRU)](https://www.uea.ac.uk/groups-and-centres/climatic-research-unit) entirely in `R`. The "Gamma test" approach as detailed in the `vignette("background")` is used. The distribution of the tree species *Cormus domestica* in Europe is selected, a tree which tolerates a wide range of conditions but favors warm to mild climates, occurring in the "Subtropical dry forest" and "Temperate Continental" FAO ecological zones. *Cormus domestica* is most commonly found in Southern Europe, though there it's natural range is uncertain owing to it's cultivation and distribution by the Roman Empire (De Rigo et al., 2016; Rotach, 2003).

## Load required packages

```{r packages, message=FALSE, warning=FALSE}
library(dplyr) # For data wrangling
library(magrittr) # For the non-base R pipe function (%>%)
library(rgbif) # For retrieving species occurrence data 
library(rnaturalearth) # For retrieving geographical data
library(getCRUCLdata) # For retrieving climate raster data
library(sf) # For spatial data operations
library(raster) # For spatial data operations
library(terra) # For spatial data operations
library(shar) # For species-habitat association analysis
library(spatstat) # For spatial point pattern analysis
library(patchwork) # For composing multiple plots
```

## Download occurrence data

To retrieve species occurrence data the `R` package `rgbif` (Chamberlain & Boettiger, 2017) is used, which provides an interface to access the GBIF database.

```{r gbif, echo=TRUE, message=FALSE, warning=FALSE}
# Retrieve key for Cormus domestica
key <- rgbif::name_backbone(name = 'Cormus domestica', kingdom = 'plants')

# Retrieve occurrences
res <- rgbif::occ_search(taxonKey = as.numeric(key$usageKey), limit = 99999)

# Create a simple data frame containing only the unique identifier (id),
# latitude (lat), and longtitude (lon).
data_simp <- data.frame(id = res$data$key, 
                        lat = res$data$decimalLatitude, lon = res$data$decimalLongitude) %>% 
  dplyr::filter(!is.na(lat) | !is.na(lon))
```

## Download map data

Spatial polygon data for the world is obtained from the `rnaturalearth` package (South, 2022), the map is then restricted to the European region. The `spatstat` package requires geospatial data in the format of a projected coordinate system; the data is therefore converted from the geographic coordinate system [4336](https://epsg.io/4326) to the projected coordinate system [3395](https://epsg.io/3395). The `shar` function `fit_point_process` requires a spatial point pattern (`ppp`) object bounded within an observation window of the class `owin`, which is then created.

```{r maps, echo=TRUE, message=FALSE, warning=FALSE}
# Retrieve data from rnaturalearth
worldmap <- rnaturalearth::ne_countries(returnclass = "sf", scale = 50) %>%
  sf::st_transform(crs = sf::st_crs(3395))

# Manually establish bounding box
eur_box <- sf::st_sfc(sf::st_point(c(-20, 30)), sf::st_point(c(45, 73)), crs = 4326) %>%
  sf::st_transform(crs = sf::st_crs(3395))

# Crop world map to include polygons within Europe extent
eur <- sf::st_crop(x = worldmap, y = eur_box)

# Define observation window
eur_owin <- spatstat.geom::as.owin(eur$geometry)
```

## Download climate data

The environmental variable selected for demonstrative purposes is the mean temperature in June over the 1961-1990 period. Data is obtained through the `getCRUCLdata` package (Sparks, 2017) which provides access to the datasets described in New et al. (2002).

```{r cru_data, echo=TRUE, message=FALSE, warning=FALSE}
# Download data as a raster brick through the getCRUCLdata package
# Mean temperature (tmn) data should be 180.4MB
cru_data <- getCRUCLdata::get_CRU_stack(pre = FALSE, pre_cv = FALSE, rd0 = FALSE,
                                        tmp = TRUE, dtr = FALSE, reh = FALSE,
                                        tmn = FALSE, tmx = FALSE, sunp = FALSE,
                                        frs = FALSE, wnd = FALSE, elv = FALSE,
                                        cache = FALSE)
```

## Prepare landscape raster

The climate data obtained above is restricted to the European region. It is then classified into 10 habitats based on temperature ranges, achieved by setting the lower and upper bounds of these ranges in the `fixedBreaks` argument of the `classify_habitats` function.

```{r landscape_ras, echo=TRUE, message=FALSE, warning=FALSE}
# Select temperature variable and the month of June
tmp_raster_jun <- cru_data$tmp$jun

# Crop tmp raster
tmp_raster_jun_eur <- terra::crop(x = tmp_raster_jun, 
                                  y = c(xmin = -20, xmax = 45, ymin = 30, ymax = 73))

# Reproject raster
tmp_raster_jun_eur_3395 <- raster::projectRaster(tmp_raster_jun_eur, crs = 3395)

# Classify landscape
landscape_classified <- shar::classify_habitats(raster = tmp_raster_jun_eur_3395,
                                                return_breaks = TRUE, style = "fixed",
                                                fixedBreaks = c(0, 5, 7.5,
                                                                10, 12.5, 15, 
                                                                17.5, 20, 25,
                                                                30, 35))
```

```{r land_plots, echo=FALSE, message=FALSE, warning=FALSE, fig.retina=FALSE, out.width="100%", dpi=400}
raster_unclassed_df <- terra::as.data.frame(tmp_raster_jun_eur_3395, xy = TRUE) %>%
  dplyr::rename("value" = "jun") %>% 
  dplyr::mutate("type" = "Unclassified")

raster_classed_df <- terra::as.data.frame(landscape_classified$raster, xy = TRUE) %>%
  dplyr::rename("value" = "layer") %>%
  dplyr::mutate("type" = "Classified")

plot_unclassed <- ggplot2::ggplot() +
  ggplot2::geom_raster(data = raster_unclassed_df,
                       mapping = ggplot2::aes(x = x, y = y, fill = value)) +
  ggplot2::geom_sf(data = eur$geometry,
                   mapping = ggplot2::aes(),
                   colour = "black", fill = NA, size = 0.1) +
  ggplot2::theme_minimal() +
  ggplot2::xlab(label = NULL) +
  ggplot2::ylab(label = NULL) +
  ggplot2::labs(fill = NULL) +
  ggplot2::scale_fill_distiller(palette = "RdBu", 
                                na.value = "transparent") +
  ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "#c9c9c9", 
                                                          linetype = "dashed", 
                                                          size = 0.075), 
                 panel.background = ggplot2::element_rect(fill = "#f0f8ff"), 
                 panel.border = ggplot2::element_rect(fill = NA),
                 text = ggplot2::element_text(size = 12),
                 axis.text.x = ggplot2::element_text(size = 9),
                 axis.text.y = ggplot2::element_text(size = 9),
                 plot.margin = ggplot2::margin(t = 0,  # Top margin
                                               r = 0,  # Right margin
                                               b = 0,  # Bottom margin
                                               l = 0),
                 legend.position = "bottom",
                 # legend.position = c(0.5, -0.2), 
                 legend.direction = "horizontal",
                 legend.justification = "center",
                 legend.text = ggplot2::element_text(size = 8),
                 legend.key.height = ggplot2::unit(0.25, 'cm'),
                 legend.key.width = ggplot2::unit(0.75, "cm")
                 )

plot_classed <- ggplot2::ggplot() +
  ggplot2::geom_raster(data = raster_classed_df,
                       mapping = ggplot2::aes(x = x, y = y,
                                              fill = factor(value))) +
  ggplot2::geom_sf(data = eur$geometry,
                   mapping = ggplot2::aes(),
                   colour = "black", fill = NA, size = 0.1) +
  ggplot2::theme_minimal() +
  ggplot2::xlab(label = NULL) +
  ggplot2::ylab(label = NULL) +
  ggplot2::labs(fill = NULL) +
  ggplot2::scale_fill_brewer(palette = "RdBu", 
                             direction = -1,
                             na.value = "transparent", 
                             guide = ggplot2::guide_legend()) +
  ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "#c9c9c9", 
                                                          linetype = "dashed", 
                                                          size = 0.075), 
                 panel.background = ggplot2::element_rect(fill = "#f0f8ff"), 
                 panel.border = ggplot2::element_rect(fill = NA),
                 text = ggplot2::element_text(size = 12),
                 axis.text.x = ggplot2::element_text(size = 9),
                 axis.text.y = ggplot2::element_blank(), # ggplot2::element_text(size = 4),
                 plot.margin = ggplot2::margin(t = 0,  # Top margin
                                               r = 0,  # Right margin
                                               b = 0,  # Bottom margin
                                               l = 0),
                 legend.position = "bottom",
                 # legend.position = c(0.5, -0.2),
                 legend.direction = "horizontal",
                 legend.justification = "center",
                 legend.text = ggplot2::element_text(size = 8),
                 legend.key.height = ggplot2::unit(0.25, 'cm'),
                 legend.key.width = ggplot2::unit(0.75, "cm")
                 )

plot_unclassed + plot_classed
```

## Prepare occurrence data

The occurrence data is prepared, then the `shar` function `fit_point_process` is called, yielding the randomized occurrence data within the observation window as required by the `results_habitat_association` function.

```{r occ_prep, echo=TRUE, message=FALSE, warning=FALSE}
# Convert occurrence data to a simple features object
data_sf <- sf::st_as_sf(data_simp, coords = c("lon", "lat"), crs = 4326)

# Restrict occurrences to those within the European region, then re-project the data
data_sf_eur_3395 <- sf::st_crop(x = data_sf,
                                y = c(xmin = -20, xmax = 45, ymin = 30, ymax = 73)) %>%
  sf::st_transform(crs = sf::st_crs(3395))

# Extract the coordinates as a matrix from the sf occurrences object
data_sf_eur_coords <- sf::st_coordinates(data_sf_eur_3395)

# Create a spatial points pattern object containing the occurrence data
data_sf_eur_ppp <- spatstat.geom::as.ppp(X = data_sf_eur_coords, W = eur_owin)

# Fit point pattern process to data
rand_pattern <- shar::fit_point_process(pattern = data_sf_eur_ppp, n_random = 19)
```

```{r occ_plots, echo=FALSE, message=FALSE, warning=FALSE, fig.retina=FALSE, out.width="100%", dpi=400}
recon_occ_df <- as.data.frame(rand_pattern$randomized$randomized_1) %>% 
  dplyr::mutate("type" = "Randomised Occurrences")

real_occ_df <- data_sf_eur_coords %>%
  as.data.frame() %>% 
  dplyr::rename(x = "X", y = "Y") %>% 
  dplyr::mutate("type" = "Real Occurrences")

real_occ_plot <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = eur$geometry,
                   mapping = ggplot2::aes(),
                   colour = "black", fill = "white", size = 0.1) +
  ggplot2::geom_point(real_occ_df, 
                      mapping = ggplot2::aes(x = x, y = y),
                      size = 0.2, stroke = 0, shape = 16, color = "Red") +
  ggplot2::theme_minimal() +
  ggplot2::xlab(label = NULL) +
  ggplot2::ylab(label = NULL) +
  ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "#c9c9c9", 
                                                          linetype = "dashed", 
                                                          size = 0.075), 
                 panel.background = ggplot2::element_rect(fill = "#f0f8ff"), 
                 panel.border = ggplot2::element_rect(fill = NA),
                 text = ggplot2::element_text(size = 12),
                 axis.text.x = ggplot2::element_text(size = 9),
                 axis.text.y = ggplot2::element_text(size = 9),
                 plot.margin = ggplot2::margin(t = 0,  
                                               r = 0,  
                                               b = 0,  
                                               l = 0))

recon_occ_plot <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = eur$geometry,
                   mapping = ggplot2::aes(),
                   colour = "black", fill = "white", size = 0.1) +
  ggplot2::geom_point(recon_occ_df, 
                      mapping = ggplot2::aes(x = x, y = y),
                      size = 0.2, stroke = 0, shape = 16, color = "Red") +
  ggplot2::theme_minimal() +
  ggplot2::xlab(label = NULL) +
  ggplot2::ylab(label = NULL) +
  ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "#c9c9c9", 
                                                          linetype = "dashed", 
                                                          size = 0.075), 
                 panel.background = ggplot2::element_rect(fill = "#f0f8ff"), 
                 panel.border = ggplot2::element_rect(fill = NA),
                 text = ggplot2::element_text(size = 12),
                 axis.text.x = ggplot2::element_text(size = 9),
                 axis.text.y = ggplot2::element_blank(), # ggplot2::element_text(size = 4),
                 plot.margin = ggplot2::margin(t = 0,  
                                               r = 0,  
                                               b = 0,  
                                               l = 0))

real_occ_plot + recon_occ_plot
```

## Results

The analysis function `results_habitat_association` is then called. The results of the analysis show that *Cormus domestica* is positively associated with locations which experience a mean June temperature of 15C - 17.5C (habitat 6) & 17.5C - 20C (habitat 7). Furthermore, *Cormus domestica* is negatively associated with all other locations classified by temperature.

```{r model_run, echo=TRUE, message=FALSE, warning=FALSE}
# Establish significance level
sig_level <- 0.01

# Run analysis
results <- shar::results_habitat_association(pattern = rand_pattern, 
                                             raster = landscape_classified$raster,
                                             breaks = landscape_classified$breaks,
                                             significance_level = sig_level) %>% 
  dplyr::arrange(habitat)

results
```

## References

Chamberlain SA, Boettiger C. 2017. R Python, and Ruby clients for GBIF species occurrence data. PeerJ Preprints 5:e3304v1 <doi:10.7287/peerj.preprints.3304v1>

De Rigo, D., Caudullo, G., Houston Durrant, T. and San-Miguel-Ayanz, J., 2016. The European Atlas of Forest Tree Species: modelling, data and information on forest tree species. *European Atlas of Forest Tree Species*, p.e01aa69. <doi:10.2788/4251>

New, M., Lister, D., Hulme, M. and Makin, I., 2002. A high-resolution data set of surface climate over global land areas. *Climate research*, *21*(1), pp.1-25. <doi:10.3354/cr021001>

Rotach, P., 2003. EUFORGEN Technical Guidelines for genetic conservation and use for service tree (Sorbus domestica). Bioversity International.

South A (2022). *rnaturalearth: World Map Data from Natural Earth*. <https://docs.ropensci.org/rnaturalearth> (website) <https://github.com/ropensci/rnaturalearth.>

Sparks, (2017). getCRUCLdata: Use and Explore CRU CL v. 2.0 Climatology Elements in R. Journal of Open Source Software, 2(12), 230, <doi:10.21105/joss.00230>
---
title: "Publication record"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Publication record}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

We try to maintain a list of publications that have used the **shar** package. If you have used **shar** in your own work, we would love to hear from you. Please [open an issue on GitHub](https://github.com/r-spatialecology/shar/issues) or [drop us an e-mail](mailto:mhk.hesselbarth@gmail.com) so we can add your work to this list.

The **shar** package is part of our academic work. To cite the package or acknowledge its use in publications, please cite the following paper.

> Hesselbarth, M.H.K., (2021). shar: A R package to analyze species-habitat associations using point pattern analysis. Journal of Open Source Software, 6(67), 3811. https://doi.org/10.21105/joss.03811

The get a BibTex entry, please use `citation("shar")`.

## A list of publications

...Work-in-Progress...
---
title: "Background"
author: "Maximilian H.K. Hesselbarth"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Background}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

Species-habitat associations are present if species are specialized to small-scale environmental conditions (Tilman & Pacala, 1993) and are more common within suitable habitats (Comita et al., 2007; Harms et al., 2001). Following, species-habitat associations show the importance of abiotic processes on the spatial patterning of plant populations (Garzon-Lopez et al., 2014).

There are mainly two methods that can be found in the literature to analyse such small-scale species-habitat associations, namely the gamma-test (Plotkin et al., 2000) and the torus-translation test (Harms et al., 2001). Both methods require data on the locations of plants (as `spatstat::ppp` point pattern) and on the environmental conditions (classified into discrete habitats as `raster::raster`). To show significance species-habitat associations, both methods randomize one of the components. Because the locations of plants as well as the habitats are most likely auto-correlated, the spatial structure must be kept while randomizing the data (Wiegand & Moloney, 2014).

All methods compare the abundance within a habitat between the observed data and the randomized null model data. If the count is below or above a pre-set threshold (e.g. the 2.5th and 97.5th quantiles), negative or positive associations, respectively, are present.

```{r overview-methods, echo = FALSE, out.height = "100%", out.width = "100%", fig.align = "center"}
knitr::include_graphics("methods_overview.png")
```

## Randomize environmental data

The following two methods randomize the environmental data, i.e. the RasterLayer data, while keeping the point pattern data fixed.

### Torus-translation-test
The torus-translation test (Harms et al. 2001) shifts the habitat map in all four cardinal directions around a torus. This is only possible for square rasters. To use this method in **shar** use `translate_raster()`.

### Randomized-habitats procedure
The randomze-habitats procedure (Harms et al. 2001) is also possible for non-square raster and randomizes the habitats using a random-walk algorithm. To use this method in **shar** use `randomize_raster()`.

```{r random-walk, echo = FALSE, out.height = "35%", out.width = "35%", fig.align = "center"}
knitr::include_graphics("random_walk.gif")
```

## Randomize point pattern data

Contrastingly to the two methods described above, the following two methods randomize the point pattern data, while keeping the environmental data fixed.

```{r gamma-vs-reconstrution, echo = FALSE, out.height = "75%", out.width = "75%", fig.align = "center"}
knitr::include_graphics("gamma_vs_reconstruction.png")
```

### Gamma-test
The gamma-test (Plotkin et al. 2000) randomizes the data by fitting a point process model to the observed data and simulation n random point patterns using the fitted point process model. However, the method only works for point patterns that can be described by a theoretical point process model. To use this method in **shar** use `fit_point_process()`.

### Pattern reconstruction
Pattern reconstruction (Tscheschel & Stoyan 2006) randomizes the point pattern using simulated annealing (Kirkpatrick et al. 1983). This allows to randomize also complex point patterns without a theoretical point process model. To use this method in **shar** use `reconstruct_pattern()`.

```{r pattern-reconstruction, echo = FALSE, out.height = "75%", out.width = "75%", fig.align = "center"}
knitr::include_graphics("pattern_reconstruction.gif")
```

### References

Comita, L.S., Condit, R., Hubbell, S.P., 2007. Developmental changes in habitat associations of tropical trees. Journal of Ecology 95, 482–492. <https://doi.org/10.1111/j.1365-2745.2007.01229.x>

Garzon-Lopez, C.X., Jansen, P.A., Bohlman, S.A., Ordonez, A., Olff, H., 2014. Effects of sampling scale on patterns of habitat association in tropical trees. Journal of Vegetation Science 25, 349–362. <https://doi.org/10.1111/jvs.12090>

Harms, K.E., Condit, R., Hubbell, S.P., Foster, R.B., 2001. Habitat associations of trees and shrubs in a 50-ha neotropical forest plot. Journal of Ecology 89, 947–959. <https://doi.org/10.1111/j.1365-2745.2001.00615.x>

Kirkpatrick, S., Gelatt, C.D.Jr., Vecchi, M.P., 1983. Optimization by simulated annealing. Science 220, 671–680. <https://doi.org/10.1126/science.220.4598.671>

Plotkin, J.B., Potts, M.D., Leslie, N., Manokaran, N., LaFrankie, J.V., Ashton, P.S., 2000. Species-area curves, spatial aggregation, and habitat specialization in tropical forests. Journal of Theoretical Biology 207, 81–99. <https://doi.org/10.1006/jtbi.2000.2158>

Tilman, D., Pacala, S.W., 1993. The maintenance of species richness in plant communities, in: Ricklefs, R.E., Schluter, D. (Eds.), Species Diversity in Ecological Communities. University of Chicago Press, Chicago, pp. 13–25. ISBN 978-0-226-71823-1

Tscheschel, A., Stoyan, D., 2006. Statistical reconstruction of random point patterns. Computational Statistics and Data Analysis 51, 859–871. <https://doi.org/10.1016/j.csda.2005.09.007>

Wiegand, T., Moloney, K.A., 2014. Handbook of spatial point-pattern analysis in ecology. Chapman and Hall/CRC Press, Boca Raton. ISBN 978-1-4200-8254-8
---
title: "How to reconstruct multiple patterns"
author: "Maximilian H.K. Hesselbarth"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{How to reconstruct multiple patterns}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  chunk_output_type: console
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

In case you want to reconstruct several patterns at once (e.g. for different points in time if repeated censuses are available), you can use the following code. Please be aware that the maximum number of iterations was set to `max_runs = 10` to keep the computational time low for this example. For an applied case, this value should be increased.

```{r load-packages, message = FALSE, warning = FALSE}
library(shar)
library(spatstat)
library(raster)
```

In case you want to only create the spatial characteristics, this is straightforward using `lapply()`. 

```{r several-patterns}
# create list with patterns
list_pattern <- list(species_a, species_b)

# reconstruct all patterns in list
result <- lapply(list_pattern, function(x) reconstruct_pattern(pattern = x, n_random = 3, 
                                                               max_runs = 10, verbose = FALSE))

```

The result will be a nested list including all *m* randomization (including the observed pattern) of the *n* provided input patterns. 

```{r result-spatial}
# get mean energy
lapply(result, function(x) calculate_energy(pattern = x,
                                            verbose = FALSE))
```

Another possible would be to first reconstruct *n* times the spatial characteristics and afterwards reconstruct the marks *m* times for each of the *n* spatial reconstructions.

Firstly, reconstruct only the spatial characteristics *n* times. The observed pattern is not needed in this case, so you can put `return_input = FALSE`.

```{r reconstruct-pattern}
# reconstruct spatial strucutre
reconstructed_pattern <- reconstruct_pattern(species_a, n_random = 3, 
                                             max_runs = 10, return_input = FALSE,
                                             verbose = FALSE)
```

Secondly, to reconstruct the (numeric) marks of the observed pattern for each of the spatially reconstructed patterns, just use `lapply()` in combination with `reconstruct_pattern_marks()`.

```{r reconstruct-marks}
# get only selected marks of input (numeric marks)
species_a_marks <- subset(species_a, select = dbh)

# reconstruct marks 3 times for each input pattern
result_marks <- lapply(reconstructed_pattern$randomized, 
                       function(x) reconstruct_pattern_marks(pattern = x, 
                                                             marked_pattern = species_a_marks, 
                                                             max_runs = 10,
                                                             n_random = 3, verbose = FALSE))
```

Again, the result is a nested list with the same dimensions as provided input patterns and reconstructions.

```{r result-marks}
# get energy
lapply(result_marks, function(x) calculate_energy(pattern = x, 
                                                  verbose = FALSE))
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classint_to_vector.R
\name{classint_to_vector}
\alias{classint_to_vector}
\title{classint_to_vector}
\usage{
classint_to_vector(x, digits = NULL)
}
\arguments{
\item{x}{classIntervals object}

\item{digits}{Integer with digits used for rounding.}
}
\value{
vector
}
\description{
Convert classIntervals to vector
}
\details{
Returns a character vector with breaks of a \code{classIntervals} object. If
\code{digits = NULL}, results will not be rounded
}
\examples{
\dontrun{
classint_to_vector(x = landscape_classified$breaks, digits = 4)
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/translate_raster.R
\name{translate_raster}
\alias{translate_raster}
\title{translate_raster}
\usage{
translate_raster(
  raster,
  steps_x = NULL,
  steps_y = NULL,
  return_input = TRUE,
  simplify = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{raster}{RasterLayer with discrete habitat classes.}

\item{steps_x, steps_y}{Integer with number of steps (cells) the raster is translated
into the corresponding direction. If both are null, all possible combinations are used
resulting in n = ((50 + 1) * (50 + 1)) - 4 rasters.}

\item{return_input}{Logical if the original input data is returned.}

\item{simplify}{Logical if only the raster will be returned if \code{n_random = 1}
and \code{return_input = FALSE}.}

\item{verbose}{Logical if progress report is printed.}
}
\value{
rd_ras
}
\description{
Torus translation
}
\details{
Torus translation test as described in Harms et al. (2001). The raster is shifted
in all four cardinal directions by steps equal to the raster resolution. If a cell
exits the extent on one side, it enters the extent on the opposite side.

The method does not allow any NA values to be present in the RasterLayer.
}
\examples{
\dontrun{
landscape_classified <- classify_habitats(landscape, n = 5, style = "fisher")

landscape_random <- translate_raster(landscape_classified)
landscape_random_sub <- translate_raster(landscape_classified,
steps_x = 1:10, steps_y = 1:5)
}

}
\references{
Harms, K.E., Condit, R., Hubbell, S.P., Foster, R.B., 2001. Habitat associations
of trees and shrubs in a 50-ha neotropical forest plot. Journal of Ecology 89, 947–959.
<https://doi.org/10.1111/j.1365-2745.2001.00615.x>
}
\seealso{
\code{\link{randomize_raster}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{landscape}
\alias{landscape}
\title{Example landscape (random cluster neutral landscape model).}
\format{
A RasterLayer object.
}
\source{
Simulated neutral landscape model with R. https://github.com/ropensci/NLMR/
}
\usage{
landscape
}
\description{
An example map to show landscapetools functionality
generated with the \code{NLMR::nlm_fbm()} algorithm.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.rd_pat.R
\name{plot.rd_pat}
\alias{plot.rd_pat}
\title{plot.rd_pat}
\usage{
\method{plot}{rd_pat}(
  x,
  what = "sf",
  n = NULL,
  probs = c(0.025, 0.975),
  comp_fast = 1000,
  ask = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{rd_pat bject with randomized patterns.}

\item{what}{Character specifying to plot summary functions of point patterns
(\code{what = "sf"}) or actual patterns (\code{what = "pp"}).}

\item{n}{Integer with number or vector of ids of randomized pattern to plot.
See Details section for more information.}

\item{probs}{Vector with quantiles of randomized data used for envelope construction.}

\item{comp_fast}{Integer with threshold at which summary functions are estimated
in a computational fast way.}

\item{ask}{Logical if the user is asked to press <RETURN> before second summary function
is plotted (only used if \code{what = "sf"}).}

\item{verbose}{Logical if progress report is printed.}

\item{...}{Not used.}
}
\value{
void
}
\description{
Plot method for rd_pat object
}
\details{
The function plots the pair correlation function and the nearest neighbour function of
the observed pattern and the reconstructed patterns (as "simulation envelopes").
For large patterns \code{comp_fast = TRUE} decreases the computational demand because no edge
correction is used and the pair correlation function is estimated based on Ripley's
K-function. For more information see \code{\link{estimate_pcf_fast}}.

It is also possible to plot n randomized patterns and the observed pattern
using \code{what = "pp"}. If \code{n} is a single number, \code{n} randomized
patterns will be sampled to plot. If \code{n} is a vector, the corresponding patterns
will be plotted.
}
\examples{
\dontrun{
pattern_random <- fit_point_process(species_a, n_random = 39)
plot(pattern_random)

pattern_recon <- reconstruct_pattern(species_b, n_random = 19,
max_runs = 1000, method = "hetero")
plot(pattern_recon)
}

}
\seealso{
\code{\link{reconstruct_pattern}} \cr
\code{\link{fit_point_process}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gamma_test}
\alias{gamma_test}
\title{Gamma test}
\format{
rd_pat object.
}
\usage{
gamma_test
}
\description{
Randomized data for species b using the gamma test.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/list_to_randomized.R
\name{list_to_randomized}
\alias{list_to_randomized}
\title{list_to_randomized}
\usage{
list_to_randomized(list, observed = NULL)
}
\arguments{
\item{list}{List}

\item{observed}{Observed}
}
\value{
rd_pat, rd_ras
}
\description{
Convert list to rd_* object.
}
\details{
Convert list of randomized point pattern or raster layer to a rd_* object that
can be used with all functions of the package. The main purpose of this utility function
is to allow an easy parallelization of the randomization approach.

For more information, please see the "Parallelization" article.
}
\examples{
\dontrun{
fit_list <- lapply(X = 1:39, FUN = function(i) {fit_point_process(pattern = species_a,
n_random = 1, simplify = TRUE, return_input = FALSE, verbose = FALSE)})

list_to_randomized(list = fit_list, observed = species_a)
}

}
\seealso{
\code{\link{randomize_raster}} \cr
\code{\link{translate_raster}} \cr
\code{\link{reconstruct_pattern}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{reconstruction}
\alias{reconstruction}
\title{Reconstruction}
\format{
rd_pat object.
}
\usage{
reconstruction
}
\description{
Randomized data for species b using pattern reconstruction.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{species_b}
\alias{species_b}
\title{Species b}
\format{
A spatstat ppp object.
}
\usage{
species_b
}
\description{
A species with positive associations to habitat 5 of \code{landscape}. Please be
aware that a positive association to one habitat will inevitable lead to negative
associations to other habitats (Yamada et al. 2006)
}
\references{
Yamada, T., Tomita, A., Itoh, A., Yamakura, T., Ohkubo, T., Kanzaki, M., Tan, S.,
Ashton, P.S., 2006. Habitat associations of Sterculiaceae trees in a Bornean rain
forest plot. Journal of Vegetation Science 17, 559–566.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_pcf_fast.R
\name{estimate_pcf_fast}
\alias{estimate_pcf_fast}
\title{estimate_pcf_fast}
\usage{
estimate_pcf_fast(pattern, ...)
}
\arguments{
\item{pattern}{ppp object with point pattern.}

\item{...}{Arguments passed down to \code{link{Kest}} or \code{\link{pcf.fv}}.}
}
\value{
fv.object
}
\description{
Fast estimation of the pair correlation function
}
\details{
The functions estimates the pair correlation functions based on an estimation
of Ripley's K-function. This makes it computationally faster than estimating the
pair correlation function directly. It is a wrapper around \code{\link{Kest}} and
\code{\link{pcf.fv}}.
}
\examples{
pcf_species_b <- estimate_pcf_fast(species_a)

}
\references{
Chiu, S.N., Stoyan, D., Kendall, W.S., Mecke, J., 2013. Stochastic geometry and
its applications, 3rd ed, Wiley Series in Probability and Statistics.
John Wiley & Sons Inc, Chichester, UK. ISBN 978-0-470-66481-0

Ripley, B.D., 1977. Modelling spatial patterns. Journal of the Royal Statistical
Society. Series B (Methodological) 39, 172–192.
<https://doi.org/10.1111/j.2517-6161.1977.tb01615.x>

Stoyan, D., Stoyan, H., 1994. Fractals, random shapes and point fields.
John Wiley & Sons, Chichester. ISBN 978-0-471-93757-9
}
\seealso{
\code{\link{Kest}} \cr
\code{\link{pcf.fv}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shar-package.R
\docType{package}
\name{shar}
\alias{shar}
\alias{shar-package}
\title{Package description}
\description{
Analyse species-habitat associations in R. Therefore, information about the
location of the species is needed and about the environmental conditions. To test
for significance habitat associations, one of the two components is randomized.
Methods are mainly based on Plotkin et al. (2000) <doi:10.1006/jtbi.2000.2158> and
Harms et al. (2001) <doi:10.1111/j.1365-2745.2001.00615.x>.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://r-spatialecology.github.io/shar/}
  \item Report bugs at \url{https://github.com/r-spatialecology/shar/issues/}
}

}
\author{
\strong{Maintainer}: Maximillian H.K. Hesselbarth \email{mhk.hesselbarth@gmail.com} (\href{https://orcid.org/0000-0003-1125-9918}{ORCID})

Authors:
\itemize{
  \item Marco Sciaini \email{marco.sciaini@posteo.net} (\href{https://orcid.org/0000-0002-3042-5435}{ORCID})
}

Other contributors:
\itemize{
  \item Zeke Marshall \email{ee18zm@leeds.ac.uk} (\href{https://orcid.org/0000-0001-9260-7827}{ORCID}) [contributor]
  \item Thomas Etherington \email{teth001@aucklanduni.ac.nz} (\href{https://orcid.org/0000-0002-3187-075X}{ORCID}) [contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reconstruct_pattern_hetero.R
\name{reconstruct_pattern_hetero}
\alias{reconstruct_pattern_hetero}
\title{reconstruct_pattern_hetero}
\usage{
reconstruct_pattern_hetero(
  pattern,
  n_random = 1,
  e_threshold = 0.01,
  max_runs = 1000,
  no_change = Inf,
  annealing = 0.01,
  comp_fast = 1000,
  weights = c(0.5, 0.5),
  r_length = 250,
  return_input = TRUE,
  simplify = FALSE,
  verbose = TRUE,
  plot = FALSE
)
}
\arguments{
\item{pattern}{ppp object with pattern.}

\item{n_random}{Integer with number of randomizations.}

\item{e_threshold}{Double with minimum energy to stop reconstruction.}

\item{max_runs}{Integer with maximum number of iterations if \code{e_threshold}
is not reached.}

\item{no_change}{Integer with number of iterations at which the reconstruction will
stop if the energy does not decrease.}

\item{annealing}{Double with probability to keep relocated point even if energy
did not decrease.}

\item{comp_fast}{Integer with threshold at which summary functions are estimated
in a computational fast way.}

\item{weights}{Vector with weights used to calculate energy.
The first number refers to Gest(r), the second number to pcf(r).}

\item{r_length}{Integer with number of intervals from \code{r = 0} to \code{r = rmax} for which
the summary functions are evaluated.}

\item{return_input}{Logical if the original input data is returned.}

\item{simplify}{Logical if only pattern will be returned if \code{n_random = 1}
and \code{return_input = FALSE}.}

\item{verbose}{Logical if progress report is printed.}

\item{plot}{Logical if pcf(r) function is plotted and updated during optimization.}
}
\value{
rd_pat
}
\description{
Pattern reconstruction for heterogeneous patterns
}
\examples{
\dontrun{
input_pattern <- spatstat.core::rpoispp(lambda = function(x, y) {100 * exp(-3 * x)},
nsim = 1)
pattern_recon <- reconstruct_pattern_hetero(input_pattern, n_random = 19, max_runs = 1000)
}

}
\references{
Kirkpatrick, S., Gelatt, C.D.Jr., Vecchi, M.P., 1983. Optimization by simulated
annealing. Science 220, 671–680. <https://doi.org/10.1126/science.220.4598.671>

Tscheschel, A., Stoyan, D., 2006. Statistical reconstruction of random point
patterns. Computational Statistics and Data Analysis 51, 859–871.
<https://doi.org/10.1016/j.csda.2005.09.007>

Wiegand, T., Moloney, K.A., 2014. Handbook of spatial point-pattern analysis in
ecology. Chapman and Hall/CRC Press, Boca Raton. ISBN 978-1-4200-8254-8
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.rd_mar.R
\name{plot.rd_mar}
\alias{plot.rd_mar}
\title{plot.rd_mar}
\usage{
\method{plot}{rd_mar}(
  x,
  what = "sf",
  n = NULL,
  probs = c(0.025, 0.975),
  comp_fast = 1000,
  ask = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{rd_mar object with randomized patterns.}

\item{what}{Character specifying to plot summary functions of point patterns
(\code{what = "sf"}) or actual patterns (\code{what = "pp"}).}

\item{n}{Integer with number or vector of ids of randomized pattern to plot.
See Details section for more information.}

\item{probs}{Vector with quantiles of randomized data used for envelope construction.}

\item{comp_fast}{Integer with threshold at which summary functions are estimated
in a computational fast way.}

\item{ask}{Logical if the user is asked to press <RETURN> before second summary function
is plotted (only used if \code{what = "sf"}).}

\item{verbose}{Logical if progress report is printed.}

\item{...}{Not used.}
}
\value{
void
}
\description{
Plot method for rd_pat object
}
\details{
The function plots the pair correlation function and the nearest neighbour function of
the observed pattern and the reconstructed patterns (as "simulation envelopes").
For large patterns \code{comp_fast = TRUE} decreases the computational demand because no edge
correction is used and the pair correlation function is estimated based on Ripley's
K-function. For more information see \code{\link{estimate_pcf_fast}}.

It is also possible to plot n randomized patterns and the observed pattern
using \code{what = "pp"}. If \code{n} is a single number, \code{n} randomized
patterns will be sampled to plot. If \code{n} is a vector, the corresponding patterns
will be plotted.
}
\examples{
\dontrun{
pattern_recon <- reconstruct_pattern(species_a, n_random = 1, max_runs = 1000,
simplify = TRUE, return_input = FALSE)
marks_sub <- spatstat.geom::subset.ppp(species_a, select = dbh)
marks_recon <- reconstruct_pattern_marks(pattern_recon, marks_sub,
n_random = 19, max_runs = 1000)
plot(marks_recon)
}


}
\seealso{
\code{\link{reconstruct_pattern}} \cr
\code{\link{fit_point_process}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_neighbourhood.R
\name{create_neighbourhood}
\alias{create_neighbourhood}
\title{create_neighbourhood}
\usage{
create_neighbourhood(cells, matrix, directions = 4)
}
\arguments{
\item{cells}{Matrix with cell ids of focal cells.}

\item{matrix}{Matrix in which cells are located.}

\item{directions}{Integer with cells neighbourhood rule: 4 (rook's case), 8 (queen's case).}
}
\value{
matrix
}
\description{
Create neighbourhood
}
\details{
Get cell ids of all neighbouring cells. The neighbourhoood rule can be specified
and is either rook's case (4 neighbours) or queen's case (8 neighbours).
}
\examples{
\dontrun{
mat <- matrix(1, nrow= 10, ncol = 10)
cell_id <- rbind(cbind(3,5), cbind(7,1))
create_neighbourhood(cell_id, mat)
}

}
\seealso{
\code{\link{randomize_raster}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.rd_pat.R
\name{print.rd_pat}
\alias{print.rd_pat}
\title{print.rd_pat}
\usage{
\method{print}{rd_pat}(x, digits = 4, ...)
}
\arguments{
\item{x}{rd_pat object with randomized patterns.}

\item{digits}{Integer with number of decimal places (round).}

\item{...}{Arguments passed to \code{cat}.}
}
\value{
void
}
\description{
Print method for rd_pat object
}
\details{
Printing method for random patterns created with \code{reconstruct_pattern_*}.
}
\examples{
pattern_random <- fit_point_process(species_a, n_random = 199)
print(pattern_random)

\dontrun{
pattern_recon <- reconstruct_pattern(species_b, n_random = 19, max_runs = 1000,
method = "hetero")
print(pattern_recon)
}

}
\seealso{
\code{\link{reconstruct_pattern}} \cr
\code{\link{fit_point_process}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_energy.R
\name{calculate_energy}
\alias{calculate_energy}
\title{calculate_energy}
\usage{
calculate_energy(
  pattern,
  weights = c(0.5, 0.5),
  return_mean = FALSE,
  comp_fast = 1000,
  verbose = TRUE
)
}
\arguments{
\item{pattern}{List with reconstructed patterns.}

\item{weights}{Vector with weights used to calculate energy.
The first number refers to Gest(r), the second number to pcf(r).}

\item{return_mean}{Logical if the mean energy is returned.}

\item{comp_fast}{Integer with threshold at which summary functions are estimated
in a computational fast way.}

\item{verbose}{Logical if progress report is printed.}
}
\value{
vector
}
\description{
Calculate mean energy
}
\details{
The function calculates the mean energy (or deviation) between the observed
pattern and all reconstructed patterns (for more information see Tscheschel &
Stoyan (2006) or Wiegand & Moloney (2014)). The pair correlation function and the
nearest neighbour distance function are used to describe the patterns. For large
patterns \code{comp_fast = TRUE} decreases the computational demand, because no edge
correction is used and the pair correlation function is estimated based on Ripley's
K-function. For more information see \code{\link{estimate_pcf_fast}}.
}
\examples{
pattern_random <- fit_point_process(species_a, n_random = 19)
calculate_energy(pattern_random)
calculate_energy(pattern_random, return_mean = TRUE)

\dontrun{
marks_sub <- spatstat.geom::subset.ppp(species_a, select = dbh)
marks_recon <- reconstruct_pattern_marks(pattern_random$randomized[[1]], marks_sub,
n_random = 19, max_runs = 1000)
calculate_energy(marks_recon, return_mean = FALSE)
}

}
\references{
Kirkpatrick, S., Gelatt, C.D.Jr., Vecchi, M.P., 1983. Optimization by simulated
annealing. Science 220, 671–680. <https://doi.org/10.1126/science.220.4598.671>

Tscheschel, A., Stoyan, D., 2006. Statistical reconstruction of random point
patterns. Computational Statistics and Data Analysis 51, 859–871.
<https://doi.org/10.1016/j.csda.2005.09.007>

Wiegand, T., Moloney, K.A., 2014. Handbook of spatial point-pattern analysis in
ecology. Chapman and Hall/CRC Press, Boca Raton. ISBN 978-1-4200-8254-8
}
\seealso{
\code{\link{plot_energy}} \cr
\code{\link{reconstruct_pattern}} \cr
\code{\link{fit_point_process}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.rd_ras.R
\name{print.rd_ras}
\alias{print.rd_ras}
\title{print.rd_ras}
\usage{
\method{print}{rd_ras}(x, ...)
}
\arguments{
\item{x}{rd_ras object with randomized raster.}

\item{...}{Arguments passed to \code{cat}.}
}
\value{
void
}
\description{
Print method for rd_ras object
}
\details{
Printing method for random patterns created with \code{\link{randomize_raster}} or
\code{\link{translate_raster}}.
}
\examples{
\dontrun{
landscape_classified <- classify_habitats(landscape, n = 5, style = "fisher")
landscape_random <- randomize_raster(landscape_classified, n_random = 19)

print(landscape_random)
}

}
\seealso{
\code{\link{randomize_raster}} \cr
\code{\link{translate_raster}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{species_a}
\alias{species_a}
\title{Species a}
\format{
A spatstat ppp object.
}
\usage{
species_a
}
\description{
A species with negative associations to habitat 4 of \code{landscape}. Please be
aware that a negative association to one habitat will inevitable lead to positive
associations to other habitats (Yamada et al. 2006).
}
\references{
Yamada, T., Tomita, A., Itoh, A., Yamakura, T., Ohkubo, T., Kanzaki, M., Tan, S.,
Ashton, P.S., 2006. Habitat associations of Sterculiaceae trees in a Bornean rain
forest plot. Journal of Vegetation Science 17, 559–566.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reconstruct_pattern_marks.R
\name{reconstruct_pattern_marks}
\alias{reconstruct_pattern_marks}
\title{reconstruct_pattern_marks}
\usage{
reconstruct_pattern_marks(
  pattern,
  marked_pattern,
  n_random = 1,
  e_threshold = 0.01,
  max_runs = 10000,
  no_change = Inf,
  annealing = 0.01,
  r_length = 250,
  return_input = TRUE,
  simplify = FALSE,
  verbose = TRUE,
  plot = FALSE
)
}
\arguments{
\item{pattern}{ppp object with pattern.}

\item{marked_pattern}{ppp  object with marked pattern. See Details section for more information.}

\item{n_random}{Integer with number of randomizations.}

\item{e_threshold}{Double with minimum energy to stop reconstruction.}

\item{max_runs}{Integer with maximum number of iterations if \code{e_threshold}
is not reached.}

\item{no_change}{Integer with number of iterations at which the reconstruction will
stop if the energy does not decrease.}

\item{annealing}{Double with probability to keep relocated point even if energy
did not decrease.}

\item{r_length}{Integer with number of intervals from \code{r = 0} to \code{r = rmax} for which
the summary functions are evaluated.}

\item{return_input}{Logical if the original input data is returned.}

\item{simplify}{Logical if only pattern will be returned if \code{n_random = 1}
and \code{return_input = FALSE}.}

\item{verbose}{Logical if progress report is printed.}

\item{plot}{Logical if pcf(r) function is plotted and updated during optimization.}
}
\value{
rd_mar
}
\description{
Pattern reconstruction of marked pattern
}
\details{
The function randomizes the numeric marks of a point pattern using pattern reconstruction
as described in Tscheschel & Stoyan (2006) and Wiegand & Moloney (2014). Therefore,
an unmarked as well as a marked pattern must be provided. The unmarked pattern must have
the spatial characteristics and the same observation window and number of points
as the marked one (see \code{reconstruct_pattern_*} or \code{\link{fit_point_process}}).
Marks must be numeric because the mark-correlation function is used as summary function.
Two randomly chosen marks are switch each iterations and changes only kept if the
deviation between the observed and the reconstructed pattern decreases.

\code{spatstat} sets \code{r_length} to 513 by default. However, a lower value decreases
the computational time while increasing the "bumpiness" of the summary function.
}
\examples{
\dontrun{
pattern_recon <- reconstruct_pattern(species_a, n_random = 1, max_runs = 1000,
simplify = TRUE, return_input = FALSE)
marks_sub <- spatstat.geom::subset.ppp(species_a, select = dbh)
marks_recon <- reconstruct_pattern_marks(pattern_recon, marks_sub,
n_random = 19, max_runs = 1000)
}

}
\references{
Kirkpatrick, S., Gelatt, C.D.Jr., Vecchi, M.P., 1983. Optimization by simulated
annealing. Science 220, 671–680. <https://doi.org/10.1126/science.220.4598.671>

Tscheschel, A., Stoyan, D., 2006. Statistical reconstruction of random point
patterns. Computational Statistics and Data Analysis 51, 859–871.
<https://doi.org/10.1016/j.csda.2005.09.007>

Wiegand, T., Moloney, K.A., 2014. Handbook of spatial point-pattern analysis in
ecology. Chapman and Hall/CRC Press, Boca Raton. ISBN 978-1-4200-8254-8
}
\seealso{
\code{\link{fit_point_process}} \cr
\code{\link{reconstruct_pattern}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reconstruct_pattern.R
\name{reconstruct_pattern}
\alias{reconstruct_pattern}
\title{reconstruct_pattern}
\usage{
reconstruct_pattern(
  pattern,
  method = "homo",
  n_random = 1,
  e_threshold = 0.01,
  max_runs = 1000,
  no_change = Inf,
  annealing = 0.01,
  comp_fast = 1000,
  n_points = NULL,
  window = NULL,
  weights = c(0.5, 0.5),
  r_length = 250,
  return_input = TRUE,
  simplify = FALSE,
  verbose = TRUE,
  plot = FALSE
)
}
\arguments{
\item{pattern}{ppp object with pattern.}

\item{method}{Character with specifying the method. Either \code{"homo"},
\code{"cluster"} or \code{"hetero"}.}

\item{n_random}{Integer with number of randomizations.}

\item{e_threshold}{Double with minimum energy to stop reconstruction.}

\item{max_runs}{Integer with maximum number of iterations if \code{e_threshold}
is not reached.}

\item{no_change}{Integer with number of iterations at which the reconstruction will
stop if the energy does not decrease.}

\item{annealing}{Double with probability to keep relocated point even if energy
did not decrease.}

\item{comp_fast}{Integer with threshold at which summary functions are estimated
in a computational fast way.}

\item{n_points}{Integer with number of points to be simulated.}

\item{window}{owin object with window of simulated pattern.}

\item{weights}{Vector with weights used to calculate energy.
The first number refers to Gest(r), the second number to pcf(r).}

\item{r_length}{Integer with number of intervals from \code{r=0} to \code{r=rmax} for which
the summary functions are evaluated.}

\item{return_input}{Logical if the original input data is returned.}

\item{simplify}{Logical if only pattern will be returned if \code{n_random=1}
and \code{return_input=FALSE}.}

\item{verbose}{Logical if progress report is printed.}

\item{plot}{Logical if pcf(r) function is plotted and updated during optimization.}
}
\value{
rd_pat
}
\description{
Pattern reconstruction
}
\details{
The functions randomizes the observed pattern by using pattern reconstruction
as described in Tscheschel & Stoyan (2006) and Wiegand & Moloney (2014). The
algorithm shifts a point to a new location and keeps the change only, if the
deviation between the observed and the reconstructed pattern decreases.
The pair correlation function and the nearest neighbour distance function are
used to describe the patterns.

For large patterns (\code{n > comp_fast}) the pair correlation function can be estimated
from Ripley's K-function without edge correction. This decreases the computational
time. For more information see \code{\link{estimate_pcf_fast}}.

The reconstruction can be stopped automatically if for n steps the energy does not
decrease. The number of steps can be controlled by \code{no_change} and is set to
\code{no_change = Inf} as default to never stop automatically.

The weights must be 0 < sum(weights) <= 1. To weight both summary functions identical,
use \code{weights = c(0.5, 0.5)}.

\code{spatstat} sets \code{r_length} to 513 by default. However, a lower value decreases
the computational time, while increasing the "bumpiness" of the summary function.

The arguments \code{n_points} and \code{window} are used for \code{method="homo"} only.

\subsection{method="homo":}{
The algorithm starts with a random pattern.
}

\subsection{method="cluster":}{
The algorithm starts with a random but clustered pattern.
}

\subsection{method="hetero":}{
The algorithm starts with a random but heterogeneous pattern.
}
}
\examples{
\dontrun{
pattern_recon <- reconstruct_pattern(species_b, n_random = 19, max_runs = 1000)
}

}
\references{
Kirkpatrick, S., Gelatt, C.D.Jr., Vecchi, M.P., 1983. Optimization by simulated
annealing. Science 220, 671–680. <https://doi.org/10.1126/science.220.4598.671>

Tscheschel, A., Stoyan, D., 2006. Statistical reconstruction of random point
patterns. Computational Statistics and Data Analysis 51, 859–871.
<https://doi.org/10.1016/j.csda.2005.09.007>

Wiegand, T., Moloney, K.A., 2014. Handbook of spatial point-pattern analysis in
ecology. Chapman and Hall/CRC Press, Boca Raton. ISBN 978-1-4200-8254-8
}
\seealso{
\code{\link{calculate_energy}} \cr
\code{\link{reconstruct_pattern_marks}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randomize_raster.R
\name{randomize_raster}
\alias{randomize_raster}
\title{randomize_raster}
\usage{
randomize_raster(
  raster,
  n_random = 1,
  directions = 4,
  return_input = TRUE,
  simplify = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{raster}{RasterLayer with discrete habitat classes.}

\item{n_random}{Integer with number of randomizations.}

\item{directions}{Interger with cells neighbourhood rule: 4 (rook's case), 8 (queen's case).}

\item{return_input}{Logical if the original input data is returned.}

\item{simplify}{Logical if only the raster will be returned if \code{n_random = 1}
and \code{return_input = FALSE}.}

\item{verbose}{Logical if progress report is printed.}
}
\value{
rd_ras
}
\description{
Randomized-habitats procedure
}
\details{
The function randomizes a habitat map with discrete classes (as RasterLayer) as proposed
by Harms et al. (2001) as “randomized-habitats procedure”. The algorithm starts with an
empty habitat map and starts to assign random neighbouring cells to each habitat
(in increasing order of abundance in observed map). We modified the procedure
slightly by increasing a probability to jump to a non-neighbouring cell as the
current patch becomes larger.

In case the RasterLayer contains NA cells, this needs to be reflected in the observation
window of the point pattern as well (i.e., no point locations possible in these areas).
}
\examples{
\dontrun{
landscape_classified <- classify_habitats(landscape, n = 5, style = "fisher")
landscape_random <- randomize_raster(landscape_classified, n_random = 19)
}

}
\references{
Harms, K.E., Condit, R., Hubbell, S.P., Foster, R.B., 2001. Habitat associations of
trees and shrubs in a 50-ha neotropical forest plot. Journal of Ecology 89, 947–959.
<https://doi.org/10.1111/j.1365-2745.2001.00615.x>
}
\seealso{
\code{\link{translate_raster}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.rd_mar.R
\name{print.rd_mar}
\alias{print.rd_mar}
\title{print.rd_mar}
\usage{
\method{print}{rd_mar}(x, digits = 4, ...)
}
\arguments{
\item{x}{rd_mar object with randomized patterns.}

\item{digits}{Integer with number of decimal places (round) to be printed.}

\item{...}{Arguments passed to \code{cat}.}
}
\value{
void
}
\description{
Print method for rd_mar object
}
\details{
Printing method for random patterns created with \code{\link{reconstruct_pattern_marks}}.
}
\examples{
\dontrun{
pattern_recon <- reconstruct_pattern(species_a, n_random = 1, max_runs = 1000,
simplify = TRUE, return_input = FALSE)
marks_sub <- spatstat.geom::subset.ppp(species_a, select = dbh)
marks_recon <- reconstruct_pattern_marks(pattern_recon, marks_sub,
n_random = 19, max_runs = 1000)
print(marks_recon)
}

}
\seealso{
\code{\link{reconstruct_pattern_marks}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_randomized.R
\name{sample_randomized}
\alias{sample_randomized}
\title{sample_randomized}
\usage{
sample_randomized(randomized, n = NULL, verbose = TRUE)
}
\arguments{
\item{randomized}{List with randomized raster or patterns.}

\item{n}{Integer with number or vector of ids of randomized pattern to plot.}

\item{verbose}{Logical if progress report is printed.}
}
\value{
list
}
\description{
Sample randomized list
}
\details{
Get list with \code{n} randomized raster or patterns. If \code{n} is a single number,
\code{n} randomized elements will be sampledt. If \code{n} is a vector, the
corresponding elements will be returned.
}
\examples{
\dontrun{
sample_randomized(randomized = reconstruction$randomized, n = c(5, 10, 15))
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{random_walk}
\alias{random_walk}
\title{Random walk}
\format{
rd_ras object.
}
\usage{
random_walk
}
\description{
Randomization of the \code{landscape} using the habitat randomization algorithm.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reconstruct_pattern_homo.R
\name{reconstruct_pattern_homo}
\alias{reconstruct_pattern_homo}
\title{reconstruct_pattern_homo}
\usage{
reconstruct_pattern_homo(
  pattern,
  n_random = 1,
  e_threshold = 0.01,
  max_runs = 1000,
  no_change = Inf,
  annealing = 0.01,
  n_points = NULL,
  window = NULL,
  comp_fast = 1000,
  weights = c(0.5, 0.5),
  r_length = 250,
  return_input = TRUE,
  simplify = FALSE,
  verbose = TRUE,
  plot = FALSE
)
}
\arguments{
\item{pattern}{ppp object with pattern.}

\item{n_random}{Integer with number of randomizations.}

\item{e_threshold}{Double with minimum energy to stop reconstruction.}

\item{max_runs}{Integer with maximum number of iterations if \code{e_threshold}
is not reached.}

\item{no_change}{Integer with number of iterations at which the reconstruction will
stop if the energy does not decrease.}

\item{annealing}{Double with probability to keep relocated point even if energy
did not decrease.}

\item{n_points}{Integer with number of points to be simulated.}

\item{window}{owin object with window of simulated pattern.}

\item{comp_fast}{Integer with threshold at which summary functions are estimated
in a computational fast way.}

\item{weights}{Vector with weights used to calculate energy.
The first number refers to Gest(r), the second number to pcf(r).}

\item{r_length}{Integer with number of intervals from \code{r = 0} to \code{r = rmax} for which
the summary functions are evaluated.}

\item{return_input}{Logical if the original input data is returned.}

\item{simplify}{Logical if only pattern will be returned if \code{n_random = 1}
and \code{return_input = FALSE}.}

\item{verbose}{Logical if progress report is printed.}

\item{plot}{Logical if pcf(r) function is plotted and updated during optimization.}
}
\value{
rd_pat
}
\description{
Pattern reconstruction for homogeneous pattern
}
\examples{
\dontrun{
pattern_recon_a <- reconstruct_pattern_homo(species_a, n_random = 19,
max_runs = 1000)

pattern_recon_b <- reconstruct_pattern_homo(species_a, n_points = 70,
n_random = 19, max_runs = 1000)
}

}
\references{
Kirkpatrick, S., Gelatt, C.D.Jr., Vecchi, M.P., 1983. Optimization by simulated
annealing. Science 220, 671–680. <https://doi.org/10.1126/science.220.4598.671>

Tscheschel, A., Stoyan, D., 2006. Statistical reconstruction of random point
patterns. Computational Statistics and Data Analysis 51, 859–871.
<https://doi.org/10.1016/j.csda.2005.09.007>

Wiegand, T., Moloney, K.A., 2014. Handbook of spatial point-pattern analysis in
ecology. Chapman and Hall/CRC Press, Boca Raton. ISBN 978-1-4200-8254-8
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify_habitats.R
\name{classify_habitats}
\alias{classify_habitats}
\title{classify_habitats}
\usage{
classify_habitats(raster, return_breaks = FALSE, ...)
}
\arguments{
\item{raster}{RasterLayer with continuous environmental values.}

\item{return_breaks}{Logical if breaks should be returned as well.}

\item{...}{Arguments passed on to \code{classIntervals}.}
}
\value{
RasterLayer
}
\description{
Classify habitats
}
\details{
Classifies a RasterLayer from the \code{raster} packages with continuous
values into n discrete classes. The \code{cut} function used to classify the raster,
uses \code{include.lowest = TRUE}.

For more information about the classification methods, see \code{classIntervals} from
the \code{classInt} package and/or the provided References. The help page of \code{classIntervals}
also includes further possible arguments to find  breaks (e.g., different styles, number
of classes, fixed breaks, etc.).
}
\examples{
landscape_classified <- classify_habitats(landscape, n = 5, style = "fisher")

landscape_classified <- classify_habitats(landscape, style = "fixed",
fixedBreaks = c(0, 0.25, 0.75, 1.0), return_breaks = TRUE)

}
\references{
Armstrong, M.P., Xiao, N., Bennett, D.A., 2003. Using genetic algorithms to create
multicriteria class intervals for choropleth maps. Annals of the Association of
American Geographers 93, 595–623. <https://doi.org/10.1111/1467-8306.9303005>

Dent, B.D., 1999. Cartography: Thematic map design, 5th ed. WCB/McGraw-Hill, Boston, USA.
ISBN 978-0-697-38495-9

Fisher, W.D., 1958. On grouping for maximum homogeneity. Journal of the American
Statistical Association 53, 789–798. <https://doi.org/10.1080/01621459.1958.10501479>

Jenks, G.F., Caspall, F.C., 1971. Error in choroplethic maps: Definition, measurement,
reduction. Annals of the Association of American Geographers 61, 217–244.
<https://doi.org/10.1111/j.1467-8306.1971.tb00779.x>

Jiang, B., 2013. Head/tail breaks: A new classification scheme for data with a
heavy-tailed distribution. The Professional Geographer 65, 482-494.
<https://doi.org/10.1080/00330124.2012.700499>

Slocum, T.A., McMaster, R.B., Kessler, F.C., Howard, H.H., 2009. Thematic cartography
and geovisualization, 3rd ed. ed, Prentice Hall Series in Geographic Information Science.
Pearson Prentice Hall, Upper Saddle River, USA. ISBN 978-0-13-229834-6

Wand, M. P., 1995. Data-based choice of histogram binwidth. The American
Statistician 51, 59-64. <https://doi.org/10.1080/00031305.1997.10473591>
}
\seealso{
\code{\link{classIntervals}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/results_habitat_association.R
\name{results_habitat_association}
\alias{results_habitat_association}
\title{results_habitat_association}
\usage{
results_habitat_association(
  pattern,
  raster,
  significance_level = 0.05,
  breaks = NULL,
  digits = NULL,
  verbose = TRUE
)
}
\arguments{
\item{pattern}{ppp object with original point pattern data or rd_pat or rd_mar
object with randomized point pattern.}

\item{raster}{RasterLayer with original discrete habitat data or rd_ras object with
randomized environmental data.}

\item{significance_level}{Double with significance level.}

\item{breaks}{Vector with breaks of habitat classes.}

\item{digits}{Integer with digits used during rounding.}

\item{verbose}{Logical if messages should be printed.}
}
\value{
data.frame
}
\description{
Results habitat association
}
\details{
The functions shows significant habitat associations by comparing the number of
points within a habitat between the observed data and randomized data as described in
Plotkin et al. (2000) and Harms et al. (2001). Significant positive or associations are present
if the observed count in a habitat is above or below a certain threshold of the
randomized count, respectively.

In case the RasterLayer contains NA cells, this needs to be reflected in the observation
window of the point pattern as well (i.e., no point locations possible in these areas).

If \code{breaks = NULL} (default), only habitat labels (but not breaks) will be
returned. If a vector with \code{breaks} is provided (same order as increasing habitat values),
the breaks will be included as well.
}
\examples{
landscape_classified <- classify_habitats(landscape, n = 5, style = "fisher")
species_a_random <- fit_point_process(species_a, n_random = 199)
results_habitat_association(pattern = species_a_random, raster = landscape_classified)

}
\references{
Harms, K.E., Condit, R., Hubbell, S.P., Foster, R.B., 2001. Habitat associations of
trees and shrubs in a 50-ha neotropical forest plot. Journal of Ecology 89, 947–959.
<https://doi.org/10.1111/j.1365-2745.2001.00615.x>

Plotkin, J.B., Potts, M.D., Leslie, N., Manokaran, N., LaFrankie, J.V.,
Ashton, P.S., 2000. Species-area curves, spatial aggregation, and habitat specialization
in tropical forests. Journal of Theoretical Biology 207, 81–99.
<https://doi.org/10.1006/jtbi.2000.2158>
}
\seealso{
\code{\link{reconstruct_pattern}} \cr
\code{\link{fit_point_process}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_points.R
\name{extract_points}
\alias{extract_points}
\title{extract_points}
\usage{
extract_points(raster, pattern)
}
\arguments{
\item{raster}{RasterLayer with environmental data}

\item{pattern}{ppp object with point pattern.}
}
\value{
data.frame
}
\description{
Extract points
}
\details{
The function extracts the number of points within each discrete habitat.
}
\examples{
\dontrun{
landscape_classified <- classify_habitats(landscape, n = 5, style = "fisher")
extract_points(raster = landscape_classified, pattern = species_b)
}

}
\seealso{
\code{\link{results_habitat_association}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{torus_trans}
\alias{torus_trans}
\title{Torus trans}
\format{
rd_ras object.
}
\usage{
torus_trans
}
\description{
Torus translation of the classified \code{landscape}.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reconstruct_pattern_cluster.R
\name{reconstruct_pattern_cluster}
\alias{reconstruct_pattern_cluster}
\title{reconstruct_pattern_cluster}
\usage{
reconstruct_pattern_cluster(
  pattern,
  n_random = 1,
  e_threshold = 0.01,
  max_runs = 1000,
  no_change = Inf,
  annealing = 0.01,
  comp_fast = 1000,
  weights = c(0.5, 0.5),
  r_length = 250,
  return_input = TRUE,
  simplify = FALSE,
  verbose = TRUE,
  plot = FALSE
)
}
\arguments{
\item{pattern}{ppp object with pattern.}

\item{n_random}{Integer with number of randomizations.}

\item{e_threshold}{Double with minimum energy to stop reconstruction.}

\item{max_runs}{Integer with maximum number of iterations if \code{e_threshold}
is not reached.}

\item{no_change}{Integer with number of iterations at which the reconstruction will
stop if the energy does not decrease.}

\item{annealing}{Double with probability to keep relocated point even if energy
did not decrease.}

\item{comp_fast}{Integer with threshold at which summary functions are estimated
in a computational fast way.}

\item{weights}{Vector with weights used to calculate energy.
The first number refers to Gest(r), the second number to pcf(r).}

\item{r_length}{Integer with number of intervals from \code{r = 0} to \code{r = rmax} for which
the summary functions are evaluated.}

\item{return_input}{Logical if the original input data is returned.}

\item{simplify}{Logical if only pattern will be returned if \code{n_random = 1}
and \code{return_input = FALSE}.}

\item{verbose}{Logical if progress report is printed.}

\item{plot}{Logical if pcf(r) function is plotted and updated during optimization.}
}
\value{
rd_pat
}
\description{
Pattern reconstruction for clustered patterns
}
\examples{
\dontrun{
pattern_recon <- reconstruct_pattern_cluster(species_b, n_random = 19, max_runs = 1000)
}

}
\references{
Kirkpatrick, S., Gelatt, C.D.Jr., Vecchi, M.P., 1983. Optimization by simulated
annealing. Science 220, 671–680. <https://doi.org/10.1126/science.220.4598.671>

Tscheschel, A., Stoyan, D., 2006. Statistical reconstruction of random point
patterns. Computational Statistics and Data Analysis 51, 859–871.
<https://doi.org/10.1016/j.csda.2005.09.007>

Wiegand, T., Moloney, K.A., 2014. Handbook of spatial point-pattern analysis in
ecology. Chapman and Hall/CRC Press, Boca Raton. ISBN 978-1-4200-8254-8
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_energy.R
\name{plot_energy}
\alias{plot_energy}
\title{plot_energy}
\usage{
plot_energy(pattern, col = NULL)
}
\arguments{
\item{pattern}{rd_pat or rd_mar object with randomized patterns.}

\item{col}{Vector with colors. Must be the same length as \code{n_random}.}
}
\value{
void
}
\description{
Plot energy of pattern reconstruction
}
\details{
The function plots the decrease of the energy over time, i.e. the iterations.
This can help to identify if the chosen \code{max_runs} for the reconstruction
were sufficient. The \code{pattern} object must have been created using
\code{reconstruct_pattern_*} .
}
\examples{
\dontrun{
pattern_recon <- reconstruct_pattern(species_a, n_random = 3, max_runs = 1000)
plot_energy(pattern_recon)

marks_sub <- spatstat.geom::subset.ppp(species_a, select = dbh)
marks_recon <- reconstruct_pattern_marks(pattern_recon$randomized[[1]], marks_sub,
n_random = 1, max_runs = 1000)
plot_energy(marks_recon)
}

}
\seealso{
\code{\link{reconstruct_pattern}} \cr
\code{\link{fit_point_process}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.rd_ras.R
\name{plot.rd_ras}
\alias{plot.rd_ras}
\title{plot.rd_ras}
\usage{
\method{plot}{rd_ras}(x, n = NULL, col, verbose = TRUE, nrow, ncol, ...)
}
\arguments{
\item{x}{rd_ras object with randomized raster.}

\item{n}{Integer with number or vector of ids of randomized raster to plot.
See Details section for more information.}

\item{col}{Vector with color palette used for plotting.}

\item{verbose}{Logical if messages are printed.}

\item{nrow, ncol}{Integer with number of rows and columns of plot grid.}

\item{...}{Not used.}
}
\value{
void
}
\description{
Plot method for rd_ras object
}
\details{
Function to plot randomized raster. If \code{n} is a single number, \code{n} randomized
raster will be sampled to plot. If \code{n} is a vector, the corresponding raster
will be plotted. \code{col, nrow, ncol} are passed to \code{plot}.
}
\examples{
\dontrun{
landscape_classified <- classify_habitats(landscape, n = 5, style = "fisher")
landscape_random <- randomize_raster(landscape_classified, n_random = 19)
plot(landscape_random)
}

}
\seealso{
\code{\link{randomize_raster}} \cr
\code{\link{translate_raster}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_point_process.R
\name{fit_point_process}
\alias{fit_point_process}
\title{fit_point_process}
\usage{
fit_point_process(
  pattern,
  n_random = 1,
  process = "poisson",
  return_input = TRUE,
  simplify = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{pattern}{ppp object with point pattern}

\item{n_random}{Integer with number of randomizations.}

\item{process}{Character specifying which point process model to use.
Either \code{"poisson"} or \code{"cluster"}.}

\item{return_input}{Logical if the original input data is returned.}

\item{simplify}{Logical if only pattern will be returned if \code{n_random = 1}
and \code{return_input = FALSE}.}

\item{verbose}{Logical if progress report is printed.}
}
\value{
rd_pat
}
\description{
Fit point process to randomize data
}
\details{
The functions randomizes the observed point pattern by fitting a point process to
the data and simulating \code{n_random} patterns using the fitted point process.
It is possible to choose between a Poisson process or a Thomas cluster process model.
For more information about the point process models, see e.g. Wiegand & Moloney (2014).
}
\examples{
pattern_fitted <- fit_point_process(pattern = species_a, n_random = 39)

}
\references{
Plotkin, J.B., Potts, M.D., Leslie, N., Manokaran, N., LaFrankie, J.V.,
Ashton, P.S., 2000. Species-area curves, spatial aggregation, and habitat specialization
in tropical forests. Journal of Theoretical Biology 207, 81–99.
<https://doi.org/10.1006/jtbi.2000.2158>

Wiegand, T., Moloney, K.A., 2014. Handbook of spatial point-pattern analysis in
ecology. Chapman and Hall/CRC Press, Boca Raton. ISBN 978-1-4200-8254-8
}
