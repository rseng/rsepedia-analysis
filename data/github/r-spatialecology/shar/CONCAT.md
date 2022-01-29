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
