
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stabm

[![R-CMD-check](https://github.com/bommert/stabm/workflows/R-CMD-check/badge.svg)](https://github.com/bommert/stabm/actions)
[![CRAN
Status](https://www.r-pkg.org/badges/version-ago/stabm)](https://cran.r-project.org/package=stabm)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03010/status.svg)](https://doi.org/10.21105/joss.03010)

`stabm` provides an implementation of many measures which assess the
stability of feature selection. The following stability measures are
currently included:

``` r
stabm::listStabilityMeasures()
#>                           Name Corrected Adjusted Minimum Maximum
#> 1               stabilityDavis     FALSE    FALSE       0       1
#> 2                stabilityDice     FALSE    FALSE       0       1
#> 3             stabilityHamming     FALSE    FALSE       0       1
#> 4   stabilityIntersectionCount      TRUE     TRUE    <NA>       1
#> 5  stabilityIntersectionGreedy      TRUE     TRUE    <NA>       1
#> 6     stabilityIntersectionMBM      TRUE     TRUE    <NA>       1
#> 7    stabilityIntersectionMean      TRUE     TRUE    <NA>       1
#> 8             stabilityJaccard     FALSE    FALSE       0       1
#> 9               stabilityKappa      TRUE    FALSE      -1       1
#> 10         stabilityLustgarten      TRUE    FALSE      -1       1
#> 11           stabilityNogueira      TRUE    FALSE      -1       1
#> 12         stabilityNovovicova     FALSE    FALSE       0       1
#> 13             stabilityOchiai     FALSE    FALSE       0       1
#> 14                stabilityPhi      TRUE    FALSE      -1       1
#> 15           stabilitySechidis     FALSE     TRUE    <NA>      NA
#> 16              stabilitySomol      TRUE    FALSE       0       1
#> 17         stabilityUnadjusted      TRUE    FALSE      -1       1
#> 18               stabilityWald      TRUE    FALSE     1-p       1
#> 19                 stabilityYu      TRUE     TRUE    <NA>       1
#> 20           stabilityZucknick     FALSE     TRUE       0       1
```

## Installation

You can install the released version of stabm from
[CRAN](https://cran.r-project.org/package=stabm) with:

``` r
install.packages("stabm")
```

For the development version, use
[devtools](https://cran.r-project.org/package=devtools):

``` r
devtools::install_github("bommert/stabm")
```

## Contributions

This R package is licensed under the
[LGPL-3](https://www.gnu.org/licenses/lgpl-3.0.en.html). If you
encounter problems using this software (lack of documentation,
misleading or wrong documentation, unexpected behaviour, bugs, …) or
just want to suggest features, please open an issue in the [issue
tracker](https://github.com/bommert/stabm/issues). Pull requests are
welcome and will be included at the discretion of the author.

## Code of Conduct

Please note that the `stabm` project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Citation

If you use stabm, please cite our [JOSS
article](https://doi.org/10.21105/joss.03010):

    @Article{stabm,
      title = {{stabm}: Stability Measures for Feature Selection},
      author = {Andrea Bommert and Michel Lang},
      journal = {Journal of Open Source Software},
      year = {2021},
      doi = {10.21105/joss.03010},
      publisher = {The Open Journal},
      volume = {6},
      number = {59},
      pages = {3010},
    }
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
---
title: 'stabm: Stability Measures for Feature Selection'
tags:
  - R
  - feature selection stability
  - stability measures
  - similarity measures
authors:
  - name: Andrea Bommert
    orcid: 0000-0002-1005-9351
    affiliation: 1
  - name: Michel Lang
    orcid: 0000-0001-9754-0393
    affiliation: 1
affiliations:
 - name: Faculty of Statistics, TU Dortmund University, 44221 Dortmund, Germany
   index: 1
date: 05 February 2021
bibliography: paper.bib
---

# Summary
The R [@R] package *stabm* provides functionality for quantifying the similarity of two or more sets.
For example, consider the two sets $\{A, B, C, D\}$ and $\{A, B, C, E\}$.
Intuitively, these sets are quite similar because their overlap is large compared to the cardinality of the two sets.
The R package *stabm* implements functions to express the similarity of sets by a real valued score.
Quantifying the similarity of sets is useful for comparing sets of selected features.
But also for many other tasks like similarity analyses of gene sets or text corpora, the R package *stabm* can be employed.

In the context of feature selection, the similarity of sets of selected features is assessed in order to determine the stability of a feature selection algorithm.
The stability of a feature selection algorithm is defined as the robustness of the set of selected features towards different data sets from the same data generating distribution [@kalousis2007stability].
For stability assessment, either *m* data sets from the same data generating process are available or *m* data sets are created from one data set.
The latter is often achieved with subsampling or random perturbations [@awada2012review].
Then, the feature selection algorithm of interest is applied to each of the *m* data sets, resulting in *m* feature sets.
To quantify the stability of the feature selection algorithm, the similarity of the *m* sets is calculated.
In the context of feature selection stability, set similarity measures are called stability measures.

The R package *stabm* provides an open-source implementation of the 20 stability measures displayed in the table below.
Argument checks are performed with checkmate [@lang2017checkmate] to provide helpful error messages.
It is publicly available on CRAN and on Github and it has only a few dependencies.

|Name | Reference|
|-----|----------|
|stabilityDavis | @davis2006reliable|
|stabilityDice | @dice1945measures|
|stabilityHamming | @dunne2002solutions|
|stabilityIntersectionCount | @bommert2020adjusted|
|stabilityIntersectionGreedy | @bommert2020adjusted|
|stabilityIntersectionMBM | @bommert2020adjusted|
|stabilityIntersectionMean | @bommert2020adjusted|
|stabilityJaccard | @jaccard1901etude|
|stabilityKappa | @carletta1996assessing|
|stabilityLustgarten | @lustgarten2009measuring|
|stabilityNogueira | @nogueira2018stability|
|stabilityNovovicova | @novovicova2009new|
|stabilityOchiai | @ochiai1957zoogeographic|
|stabilityPhi | @nogueira2016measuring|
|stabilitySechidis | @sechidis2020stability|
|stabilitySomol | @somol2008evaluating|
|stabilityUnadjusted | @bommert2020adjusted|
|stabilityWald | @wald2013stability|
|stabilityYu | @yu2012stable|
|stabilityZucknick | @zucknick2008comparing|

# Statement of Need
The R package *stabm* provides an implementation of many stability measures.
For theoretical and empirical comparative studies of the stability measures implemented in *stabm*, we refer to @bommert2017multicriteria, @bommert2020adjusted, @bommert2020integration, and @nogueira2018stability.
It has been demonstrated that considering the feature selection stability when fitting a predictive model often is beneficial for obtaining models with high predictive accuracy [@bommert2017multicriteria; @bommert2020integration; @schirra2016selection].
The stability measures implemented in the R package *stabm* have been employed in @bommert2017multicriteria, @bommert2020benchmark, @bommert2020adjusted, and @bommert2020integration.

# Related Software
A subset of the implemented stability measures is also available in other R or Python packages.
The R package *sets* [@meyer2009sets] and the Python package *scikit-learn* [@pedregosa2011scikit] provide an implementation of the Jaccard index [@jaccard1901etude] to assess the similarity of two sets.
The Python package *GSimPy* [@zhang2020gsimpy] implements the Jaccard index, the Dice index [@dice1945measures], and the Ochiai index [@ochiai1957zoogeographic].
The source code for the publication @nogueira2018stability provides an implementation of their stability measure in R, Python, and Matlab.

# Acknowledgements

This work was supported by the German Research Foundation (DFG), Project RA 870/7-1, and Collaborative Research Center SFB 876, A3.

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
# stabm

[![R-CMD-check](https://github.com/bommert/stabm/workflows/R-CMD-check/badge.svg)](https://github.com/bommert/stabm/actions)
[![CRAN Status](https://www.r-pkg.org/badges/version-ago/stabm)](https://cran.r-project.org/package=stabm)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03010/status.svg)](https://doi.org/10.21105/joss.03010)


`stabm` provides an implementation of many measures which assess the stability of feature selection.
The following stability measures are currently included:

```{r}
stabm::listStabilityMeasures()
```


## Installation

You can install the released version of stabm from [CRAN](https://cran.r-project.org/package=stabm) with:

```{r, eval=FALSE}
install.packages("stabm")
```
For the development version, use [devtools](https://cran.r-project.org/package=devtools):

```{r, eval = FALSE}
devtools::install_github("bommert/stabm")
```


## Contributions

This R package is licensed under the [LGPL-3](https://www.gnu.org/licenses/lgpl-3.0.en.html).
If you encounter problems using this software (lack of documentation, misleading or wrong documentation, unexpected behaviour, bugs, ...) or just want to suggest features,
please open an issue in the [issue tracker](https://github.com/bommert/stabm/issues).
Pull requests are welcome and will be included at the discretion of the author.


## Code of Conduct

Please note that the `stabm` project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.


## Citation

If you use stabm, please cite our [JOSS article](https://doi.org/10.21105/joss.03010):
```{r echo = FALSE, comment = ""}
toBibtex(citation("stabm"))
```---
title: "stabm"
author: "Andrea Bommert"
date: "`r Sys.Date()`"
output: rmarkdown::pdf_document
vignette: >
  %\VignetteIndexEntry{stabm}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: references.bib
---
```{r,include=FALSE}
library(stabm)
```

# Introduction
The `R` package `stabm` provides functionality for quantifying the similarity of two or more sets.
The anticipated usecase is comparing sets of selected features, but other sets, e.g. gene list, can be analyzed as well.
Quantifying the similarity of feature sets is necessary when assessing the feature selection stability.
The stability of a feature selection algorithm is defined as the robustness of the set of selected features towards different data sets from the same data generating distribution [@kalousis2007stability].
Stability measures quantify the similarity of the sets of selected features for different training data sets.
Many stability measures have been proposed in the literature, see for example @bommert2017multicriteria, @bommert2020adjusted, @bommert2020integration and @nogueira2018stability for comparative studies.
The `R` package `stabm` provides an implementation of many stability measures.
Detailed definitions and analyses of all stability measures implemented in `stabm` are given in @bommert2020integration.


# Usage
A list of all stability measures implemented in `stabm` is available with:

```{r}
listStabilityMeasures()
```

This list states the names of the stability measures and some information about them.

- Corrected: Does a measure fulfill the property *correction for chance* as defined in @nogueira2018stability? This property indicates whether the expected value of the stability measure is independent of the number of selected features. Stability measures not fulfilling this property, usually attain the higher values, the more features are selected. For the measures that are not corrected for chance in their original definition, `stabm` provides the possibility to transform these measures, such that they are corrected for chance.
- Adjusted: Does a measure consider similarities between features when evaluating the feature selection stability? Adjusted measures have been created based on traditional stability measures by including an adjustment term that takes into account feature similarities, see @bommert2020integration for details.
- Minimum and Maximum: Bounds for the stability measures, useful for interpreting obtained stability values.

Now, let us consider an example with 3 sets of selected features

- $V_1 = \{X_1, X_2, X_3\}$
- $V_2 = \{X_1, X_2, X_3, X_4\}$
- $V_3 = \{X_1, X_2, X_3, X_5, X_6, X_7\}$

and a total number of 10 features. We can evaluate the feature selection stability with stability measures of our choice.

```{r}
feats = list(1:3, 1:4, c(1:3, 5:7))
stabilityJaccard(features = feats)
stabilityNogueira(features = feats, p = 10)
```

For adjusted stability measures, a matrix indicating the similarities between the features has to be specified.

```{r}
mat = 0.92 ^ abs(outer(1:10, 1:10, "-"))
set.seed(1)
stabilityIntersectionCount(features = feats, sim.mat = mat, N = 1000)
```

Finally, `stabm` also provides a visualization of the feature sets.

```{r, fig.width=4.5, fig.height=3, fig.align="center", message=FALSE}
plotFeatures(feats)
```


# Example: Feature Selection

In this example, we will analyze the stability of the feature selection of regression trees
on the `BostonHousing2` data set from the [mlbench package](https://cran.r-project.org/package=mlbench).

```{r}
library(rpart) # for classification trees
data("BostonHousing2", package = "mlbench")

# remove feature that is a version of the target variable
dataset = subset(BostonHousing2, select = -cmedv)
```

We write a small function which subsamples the `BostonHousing2` data frame to `ratio` percent of the observations, fits a regression tree and then returns the used features as character vector:

```{r}
fit_tree = function(target = "medv", data = dataset, ratio = 0.67, cp = 0.01) {
    n = nrow(data)
    i = sample(n, n * ratio)
    formula = as.formula(paste(target,  "~ ."))
    model = rpart::rpart(formula = formula, data = data, subset = i, 
      control = rpart.control(maxsurrogate = 0, cp = cp))
    names(model$variable.importance)
}

set.seed(1)
fit_tree()
```


We repeat this step 30 times, resulting in a list of character vectors of selected features:
```{r}
set.seed(1)
selected_features = replicate(30, fit_tree(), simplify = FALSE)
```

A quick analysis of the list reveals that three features are selected in all repetitions while six other features are only selected in some of the repetitions:
```{r}
# Selected in each repetition:
Reduce(intersect, selected_features)

# Sorted selection frequency across all 30 repetitions:
sort(table(unlist(selected_features)), decreasing = TRUE)
```

The selection frequency can be visualized with the `plotFeatures()` function:
```{r}
plotFeatures(selected_features)
```

To finally express the selection frequencies with one number, e.g. to compare the stability of regression trees to the stability of a different modeling approach, any of the implemented stability measures can be calculated:
```{r}
stabilityJaccard(selected_features)
```

We consider a second parametrization of regression trees and observe that this parametrization provides a more stable feature selection than the default parametrization: the value of the Jaccard stability measure is higher here:
```{r}
set.seed(1)
selected_features2 = replicate(30, fit_tree(cp = 0.02), simplify = FALSE)
stabilityJaccard(selected_features2)
plotFeatures(selected_features2)
```

Now, we consider a different regression problem, for which there are highly correlated features.
Again, we repeatedly select features using regression trees:
```{r}
dataset2 = subset(BostonHousing2, select = -town)
dataset2$chas = as.numeric(dataset2$chas)

set.seed(1)
selected_features3 = replicate(30, fit_tree(target = "rm", data = dataset2, cp = 0.075), 
  simplify = FALSE)
```

We choose to assess the similarities between the features with absolute Pearson correlations, but other similarity measures could be used as well.
The similarity values of the selected features show that the two features *medv* and *cmedv* are almost perfectly correlated:
```{r}
# similarity matrix
sim.mat = abs(cor(subset(dataset2, select = -rm)))

sel.feats = unique(unlist(selected_features3))
sim.mat[sel.feats, sel.feats]
```
Also, each of the 30 feature sets includes either *medv* or *cmedv*:
```{r}
plotFeatures(selected_features3, sim.mat = sim.mat)
```

When evaluating the feature selection stability, we want that the choice of *medv* instead of *cmedv* or vice versa is not seen as a lack of stability, because they contain almost the same information.
Therefore, we use one of the *adjusted* stability measures, see `listStabilityMeasures()` in Section *Usage*.
```{r}
stabilityIntersectionCount(selected_features3, sim.mat = sim.mat, N = 1000)
```

The effect of the feature similarities for stability assessment can be quantified by considering the identity matrix as similarity matrix and thereby neglecting all similarities.
Without taking into account the feature similarities, the stability value is much lower:
```{r}
no.sim.mat = diag(nrow(sim.mat))
colnames(no.sim.mat) = row.names(no.sim.mat) = colnames(sim.mat)
stabilityIntersectionCount(selected_features3, sim.mat = no.sim.mat)
```


# Example: Clustering

As a second example, we analyze the stability of the clusters resulting from k-means clustering.
```{r}
set.seed(1)

# select a subset of instances for visualization purposes
inds = sample(nrow(dataset2), 50)
dataset.cluster = dataset2[inds, ]

# run k-means clustering with k = 3 30 times
km = replicate(30, kmeans(dataset.cluster, centers = 3), simplify = FALSE)

# change cluster names for comparability
best = which.min(sapply(km, function(x) x$tot.withinss))
best.centers = km[[best]]$centers
km.clusters = lapply(km, function(kmi) {
  dst = as.matrix(dist(rbind(best.centers, kmi$centers)))[4:6, 1:3]
  rownames(dst) = colnames(dst) = 1:3
  # greedy choice of best matches of clusters
  new.cluster.names = numeric(3)
  while(nrow(dst) > 0) {
    min.dst = which.min(dst)
    row = (min.dst - 1) %% nrow(dst) + 1
    row.o = as.numeric(rownames(dst)[row])
    col = ceiling(min.dst / nrow(dst))
    col.o = as.numeric(colnames(dst)[col])
    new.cluster.names[row.o] = col.o
    dst = dst[-row, -col, drop = FALSE]
  }
  new.cluster.names[kmi$cluster]
})

# for each cluster, create a list containing the instances 
# belonging to this cluster over the 30 repetitions
clusters = lapply(1:3, function(i) {
  lapply(km.clusters, function(kmc) {
    which(kmc == i)
  })
})
```

For each cluster, we evaluate the stability of the instances assigned to this cluster:
```{r}
stab.cl = sapply(clusters, stabilityJaccard)
stab.cl
```

We average these stability values with a weighted mean based on the average cluster sizes:
```{r}
w = sapply(clusters, function(cl) {
  mean(lengths(cl))
})

sum(stab.cl * w) / sum(w)
```

# References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/documentation_stability.R
\name{uncorrectedDocumentation}
\alias{uncorrectedDocumentation}
\title{Uncorrected Stability Measures}
\arguments{
\item{p}{\code{numeric(1)}\cr
Total number of features in the datasets.
Required, if \code{correction.for.chance} is set to "estimate" or "exact".}

\item{correction.for.chance}{\code{character(1)}\cr
Should a correction for chance be applied? Correction for chance means that if
features are chosen at random, the expected value must be independent of the number
of chosen features. To correct for chance, the original score is transformed by
\eqn{(score - expected) / (maximum - expected)}. For stability measures whose
score is the average value of pairwise scores, this transformation
is done for all components individually.
Options are "none", "estimate" and "exact".
For "none", no correction is performed, i.e. the original score is used.
For "estimate", \code{N} random feature sets of the same sizes as the input
feature sets (\code{features}) are generated.
For "exact", all possible combinations of feature sets of the same
sizes as the input feature sets are used. Computation is only feasible for very
small numbers of features (\code{p}) and numbers of considered datasets
(\code{length(features)}).}
}
\description{
Uncorrected Stability Measures
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_uncorrected.R
\encoding{UTF-8}
\name{stabilityDice}
\alias{stabilityDice}
\title{Stability Measure Dice}
\usage{
stabilityDice(
  features,
  p = NULL,
  correction.for.chance = "none",
  N = 10000,
  impute.na = NULL
)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{p}{\code{numeric(1)}\cr
Total number of features in the datasets.
Required, if \code{correction.for.chance} is set to "estimate" or "exact".}

\item{correction.for.chance}{\code{character(1)}\cr
Should a correction for chance be applied? Correction for chance means that if
features are chosen at random, the expected value must be independent of the number
of chosen features. To correct for chance, the original score is transformed by
\eqn{(score - expected) / (maximum - expected)}. For stability measures whose
score is the average value of pairwise scores, this transformation
is done for all components individually.
Options are "none", "estimate" and "exact".
For "none", no correction is performed, i.e. the original score is used.
For "estimate", \code{N} random feature sets of the same sizes as the input
feature sets (\code{features}) are generated.
For "exact", all possible combinations of feature sets of the same
sizes as the input feature sets are used. Computation is only feasible for very
small numbers of features (\code{p}) and numbers of considered datasets
(\code{length(features)}).}

\item{N}{\code{numeric(1)}\cr
Number of random feature sets to consider. Only relevant if \code{correction.for.chance}
is set to "estimate".}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as (see Notation)
\deqn{\frac{2}{m (m - 1)} \sum_{i=1}^{m-1} \sum_{j = i+1}^m
\frac{2 |V_i \cap V_j|}{|V_i| + |V_j|}.}
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
stabilityDice(features = feats)
}
\references{
Dice LR (1945).
\dQuote{Measures of the Amount of Ecologic Association Between Species.}
\emph{Ecology}, \bold{26}(3), 297--302.
\doi{10.2307/1932409}.

Bommert A, Rahnenführer J, Lang M (2017).
\dQuote{A Multicriteria Approach to Find Predictive and Sparse Models with Stable Feature Selection for High-Dimensional Data.}
\emph{Computational and Mathematical Methods in Medicine}, \bold{2017}, 1--18.
\doi{10.1155/2017/7907163}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/documentation_stability.R
\name{adjustedDocumentation}
\alias{adjustedDocumentation}
\title{Adjusted Stability Measures}
\arguments{
\item{correction.for.chance}{\code{character(1)}\cr
How should the expected value of the stability score (see Details)
be assessed? Options are "estimate", "exact" and "none".
For "estimate", \code{N} random feature sets of the same sizes as the input feature
sets (\code{features}) are generated.
For "exact", all possible combinations of feature sets of the same
sizes as the input feature sets are used. Computation is only feasible for very
small numbers of features and numbers of considered datasets (\code{length(features)}).
For "none", the transformation \eqn{(score - expected) / (maximum - expected)}
is not conducted, i.e. only \eqn{score} is used.
This is not recommended.}
}
\description{
Adjusted Stability Measures
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_uncorrected.R
\encoding{UTF-8}
\name{stabilityNovovicova}
\alias{stabilityNovovicova}
\title{Stability Measure Novovičová}
\usage{
stabilityNovovicova(
  features,
  p = NULL,
  correction.for.chance = "none",
  N = 10000,
  impute.na = NULL
)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{p}{\code{numeric(1)}\cr
Total number of features in the datasets.
Required, if \code{correction.for.chance} is set to "estimate" or "exact".}

\item{correction.for.chance}{\code{character(1)}\cr
Should a correction for chance be applied? Correction for chance means that if
features are chosen at random, the expected value must be independent of the number
of chosen features. To correct for chance, the original score is transformed by
\eqn{(score - expected) / (maximum - expected)}. For stability measures whose
score is the average value of pairwise scores, this transformation
is done for all components individually.
Options are "none", "estimate" and "exact".
For "none", no correction is performed, i.e. the original score is used.
For "estimate", \code{N} random feature sets of the same sizes as the input
feature sets (\code{features}) are generated.
For "exact", all possible combinations of feature sets of the same
sizes as the input feature sets are used. Computation is only feasible for very
small numbers of features (\code{p}) and numbers of considered datasets
(\code{length(features)}).}

\item{N}{\code{numeric(1)}\cr
Number of random feature sets to consider. Only relevant if \code{correction.for.chance}
is set to "estimate".}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as (see Notation)
\deqn{\frac{1}{q \log_2(m)} \sum_{j: X_j \in V} h_j \log_2(h_j).}
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
stabilityNovovicova(features = feats)
}
\references{
Novovičová J, Somol P, Pudil P (2009).
\dQuote{A New Measure of Feature Selection Algorithms' Stability.}
In \emph{2009 IEEE International Conference on Data Mining Workshops}.
\doi{10.1109/icdmw.2009.32}.

Bommert A, Rahnenführer J, Lang M (2017).
\dQuote{A Multicriteria Approach to Find Predictive and Sparse Models with Stable Feature Selection for High-Dimensional Data.}
\emph{Computational and Mathematical Methods in Medicine}, \bold{2017}, 1--18.
\doi{10.1155/2017/7907163}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_adjusted.R
\encoding{UTF-8}
\name{stabilityIntersectionMBM}
\alias{stabilityIntersectionMBM}
\title{Stability Measure Adjusted Intersection MBM}
\usage{
stabilityIntersectionMBM(
  features,
  sim.mat,
  threshold = 0.9,
  correction.for.chance = "estimate",
  N = 10000,
  impute.na = NULL
)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{sim.mat}{\code{numeric matrix}\cr
Similarity matrix which contains the similarity structure of all features based on
all datasets. The similarity values must be in the range of [0, 1] where 0 indicates
very low similarity and 1 indicates very high similarity. If the list elements of
\code{features} are integerish vectors, then the feature numbering must correspond to the
ordering of \code{sim.mat}. If the list elements of \code{features} are character
vectors, then \code{sim.mat} must be named and the names of \code{sim.mat} must correspond
to the entries in \code{features}.}

\item{threshold}{\code{numeric(1)}\cr
Threshold for indicating which features are similar and which are not. Two features
are considered as similar, if and only if the corresponding entry of \code{sim.mat} is greater
than or equal to \code{threshold}.}

\item{correction.for.chance}{\code{character(1)}\cr
How should the expected value of the stability score (see Details)
be assessed? Options are "estimate", "exact" and "none".
For "estimate", \code{N} random feature sets of the same sizes as the input feature
sets (\code{features}) are generated.
For "exact", all possible combinations of feature sets of the same
sizes as the input feature sets are used. Computation is only feasible for very
small numbers of features and numbers of considered datasets (\code{length(features)}).
For "none", the transformation \eqn{(score - expected) / (maximum - expected)}
is not conducted, i.e. only \eqn{score} is used.
This is not recommended.}

\item{N}{\code{numeric(1)}\cr
Number of random feature sets to consider. Only relevant if \code{correction.for.chance}
is set to "estimate".}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as (see Notation)
\deqn{\frac{2}{m(m-1)}\sum_{i=1}^{m-1} \sum_{j=i+1}^{m}
\frac{I(V_i, V_j) - E(I(V_i, V_j))}{\sqrt{|V_i| \cdot |V_j|} - E(I(V_i, V_j))}}
with \deqn{I(V_i, V_j) = |V_i \cap V_j| + \mathop{\mathrm{MBM}}(V_i \setminus V_j, V_j \backslash V_i).}
\eqn{\mathop{\mathrm{MBM}}(V_i \setminus V_j, V_j \backslash V_i)} denotes the size of the
maximum bipartite matching based on the graph whose vertices are the features
of \eqn{V_i \setminus V_j} on the one side and the features of \eqn{V_j \backslash V_i}
on the other side. Vertices x and y are connected if and only if \eqn{\mathrm{Similarity}(x, y)
\geq \mathrm{threshold}.}
Requires the package \CRANpkg{igraph}.
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
mat = 0.92 ^ abs(outer(1:10, 1:10, "-"))
stabilityIntersectionMBM(features = feats, sim.mat = mat, N = 1000)
}
\references{
Bommert A, Rahnenführer J (2020).
\dQuote{Adjusted Measures for Feature Selection Stability for Data Sets with Similar Features.}
In \emph{Machine Learning, Optimization, and Data Science}, 203--214.
\doi{10.1007/978-3-030-64583-0_19}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_adjusted.R
\name{stabilityYu}
\alias{stabilityYu}
\title{Stability Measure Yu}
\usage{
stabilityYu(
  features,
  sim.mat,
  threshold = 0.9,
  correction.for.chance = "estimate",
  N = 10000,
  impute.na = NULL
)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{sim.mat}{\code{numeric matrix}\cr
Similarity matrix which contains the similarity structure of all features based on
all datasets. The similarity values must be in the range of [0, 1] where 0 indicates
very low similarity and 1 indicates very high similarity. If the list elements of
\code{features} are integerish vectors, then the feature numbering must correspond to the
ordering of \code{sim.mat}. If the list elements of \code{features} are character
vectors, then \code{sim.mat} must be named and the names of \code{sim.mat} must correspond
to the entries in \code{features}.}

\item{threshold}{\code{numeric(1)}\cr
Threshold for indicating which features are similar and which are not. Two features
are considered as similar, if and only if the corresponding entry of \code{sim.mat} is greater
than or equal to \code{threshold}.}

\item{correction.for.chance}{\code{character(1)}\cr
How should the expected value of the stability score (see Details)
be assessed? Options are "estimate", "exact" and "none".
For "estimate", \code{N} random feature sets of the same sizes as the input feature
sets (\code{features}) are generated.
For "exact", all possible combinations of feature sets of the same
sizes as the input feature sets are used. Computation is only feasible for very
small numbers of features and numbers of considered datasets (\code{length(features)}).
For "none", the transformation \eqn{(score - expected) / (maximum - expected)}
is not conducted, i.e. only \eqn{score} is used.
This is not recommended.}

\item{N}{\code{numeric(1)}\cr
Number of random feature sets to consider. Only relevant if \code{correction.for.chance}
is set to "estimate".}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
Let \eqn{O_{ij}} denote the number of features in \eqn{V_i} that are not
shared with \eqn{V_j} but that have a highly simlar feature in \eqn{V_j}:
\deqn{O_{ij} = |\{ x \in (V_i \setminus V_j) : \exists y \in (V_j \backslash V_i) \ with \
Similarity(x,y) \geq threshold \}|.}
Then the stability measure is defined as (see Notation)
\deqn{\frac{2}{m(m-1)}\sum_{i=1}^{m-1} \sum_{j=i+1}^{m}
\frac{I(V_i, V_j) - E(I(V_i, V_j))}{\frac{|V_i| + |V_j|}{2} - E(I(V_i, V_j))}} with
\deqn{I(V_i, V_j) = |V_i \cap V_j| + \frac{O_{ij} + O_{ji}}{2}.}
Note that this definition slightly differs from its original in order to make it suitable
for arbitrary datasets and similarity measures and applicable in situations with \eqn{|V_i| \neq |V_j|}.
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
mat = 0.92 ^ abs(outer(1:10, 1:10, "-"))
stabilityYu(features = feats, sim.mat = mat, N = 1000)
}
\references{
Yu L, Han Y, Berens ME (2012).
\dQuote{Stable Gene Selection from Microarray Data via Sample Weighting.}
\emph{IEEE/ACM Transactions on Computational Biology and Bioinformatics}, \bold{9}(1), 262--272.
\doi{10.1109/tcbb.2011.47}.

Zhang M, Zhang L, Zou J, Yao C, Xiao H, Liu Q, Wang J, Wang D, Wang C, Guo Z (2009).
\dQuote{Evaluating reproducibility of differential expression discoveries in microarray studies by considering correlated molecular changes.}
\emph{Bioinformatics}, \bold{25}(13), 1662--1668.
\doi{10.1093/bioinformatics/btp295}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\docType{package}
\name{stabm-package}
\alias{stabm}
\alias{stabm-package}
\title{stabm: Stability Measures for Feature Selection}
\description{
An implementation of many measures for the
    assessment of the stability of feature selection. Both simple measures
    and measures which take into account the similarities between features
    are available, see Bommert (2020) <doi:10.17877/DE290R-21906>.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://bommert.github.io/stabm/}
  \item \url{https://github.com/bommert/stabm}
  \item Report bugs at \url{https://github.com/bommert/stabm/issues}
}

}
\author{
\strong{Maintainer}: Andrea Bommert \email{bommert@statistik.tu-dortmund.de} (\href{https://orcid.org/0000-0002-1005-9351}{ORCID})

Authors:
\itemize{
  \item Michel Lang \email{michellang@gmail.com} (\href{https://orcid.org/0000-0001-9754-0393}{ORCID})
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_adjusted.R
\encoding{UTF-8}
\name{stabilitySechidis}
\alias{stabilitySechidis}
\title{Stability Measure Sechidis}
\usage{
stabilitySechidis(features, sim.mat, threshold = 0.9, impute.na = NULL)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{sim.mat}{\code{numeric matrix}\cr
Similarity matrix which contains the similarity structure of all features based on
all datasets. The similarity values must be in the range of [0, 1] where 0 indicates
very low similarity and 1 indicates very high similarity. If the list elements of
\code{features} are integerish vectors, then the feature numbering must correspond to the
ordering of \code{sim.mat}. If the list elements of \code{features} are character
vectors, then \code{sim.mat} must be named and the names of \code{sim.mat} must correspond
to the entries in \code{features}.}

\item{threshold}{\code{numeric(1)}\cr
Threshold for indicating which features are similar and which are not. Two features
are considered as similar, if and only if the corresponding entry of \code{sim.mat} is greater
than or equal to \code{threshold}.}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as
\deqn{1 - \frac{\mathop{\mathrm{trace}}(CS)}{\mathop{\mathrm{trace}}(C \Sigma)}} with (\eqn{p \times p})-matrices
\deqn{(S)_{ij} = \frac{m}{m-1}\left(\frac{h_{ij}}{m} - \frac{h_i}{m} \frac{h_j}{m}\right)} and
\deqn{(\Sigma)_{ii} = \frac{q}{mp} \left(1 - \frac{q}{mp}\right),}
\deqn{(\Sigma)_{ij} = \frac{\frac{1}{m} \sum_{i=1}^{m} |V_i|^2 - \frac{q}{m}}{p^2 - p} - \frac{q^2}{m^2 p^2}, i \neq j.}
The matrix \eqn{C} is created from matrix \code{sim.mat} by setting all values of \code{sim.mat} that are smaller
than \code{threshold} to 0. If you want to \eqn{C} to be equal to \code{sim.mat}, use \code{threshold = 0}.
}
\note{
This stability measure is not corrected for chance.
Unlike for the other stability measures in this R package, that are not corrected for chance,
for \code{stabilitySechidis}, no \code{correction.for.chance} can be applied.
This is because for \code{stabilitySechidis}, no finite upper bound is known at the moment,
see \link{listStabilityMeasures}.
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
mat = 0.92 ^ abs(outer(1:10, 1:10, "-"))
stabilitySechidis(features = feats, sim.mat = mat)
}
\references{
Sechidis K, Papangelou K, Nogueira S, Weatherall J, Brown G (2020).
\dQuote{On the Stability of Feature Selection in the Presence of Feature Correlations.}
In \emph{Machine Learning and Knowledge Discovery in Databases}, 327--342.
Springer International Publishing.
\doi{10.1007/978-3-030-46150-8_20}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_corrected.R
\encoding{UTF-8}
\name{stabilityWald}
\alias{stabilityWald}
\title{Stability Measure Wald}
\usage{
stabilityWald(features, p, impute.na = NULL)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{p}{\code{numeric(1)}\cr
Total number of features in the datasets.}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as (see Notation)
\deqn{\frac{2}{m (m - 1)} \sum_{i=1}^{m-1} \sum_{j = i+1}^m
\frac{|V_i \cap V_j| - \frac{|V_i| \cdot |V_j|}{p}}
{\min \{|V_i|, |V_j|\} - \frac{|V_i| \cdot |V_j|}{p}}.}
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
stabilityWald(features = feats, p = 10)
}
\references{
Wald R, Khoshgoftaar TM, Napolitano A (2013).
\dQuote{Stability of Filter- and Wrapper-Based Feature Subset Selection.}
In \emph{2013 IEEE 25th International Conference on Tools with Artificial Intelligence}.
\doi{10.1109/ictai.2013.63}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_corrected.R
\encoding{UTF-8}
\name{stabilityPhi}
\alias{stabilityPhi}
\title{Stability Measure Phi}
\usage{
stabilityPhi(features, p, impute.na = NULL)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{p}{\code{numeric(1)}\cr
Total number of features in the datasets.}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as the average
phi coefficient between all pairs of feature sets.
It can be rewritten as (see Notation)
\deqn{\frac{2}{m (m - 1)} \sum_{i=1}^{m-1} \sum_{j = i+1}^m
\frac{|V_i \cap V_j| - \frac{|V_i| \cdot |V_j|}{p}}
{\sqrt{|V_i| (1 - \frac{|V_i|}{p}) \cdot |V_j| (1 - \frac{|V_j|}{p})}}.}
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
stabilityPhi(features = feats, p = 10)
}
\references{
Nogueira S, Brown G (2016).
\dQuote{Measuring the Stability of Feature Selection.}
In \emph{Machine Learning and Knowledge Discovery in Databases}, 442--457.
Springer International Publishing.
\doi{10.1007/978-3-319-46227-1_28}.

Bommert A, Rahnenführer J, Lang M (2017).
\dQuote{A Multicriteria Approach to Find Predictive and Sparse Models with Stable Feature Selection for High-Dimensional Data.}
\emph{Computational and Mathematical Methods in Medicine}, \bold{2017}, 1--18.
\doi{10.1155/2017/7907163}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_corrected.R
\encoding{UTF-8}
\name{stabilityKappa}
\alias{stabilityKappa}
\title{Stability Measure Kappa}
\usage{
stabilityKappa(features, p, impute.na = NULL)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{p}{\code{numeric(1)}\cr
Total number of features in the datasets.}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as the average
kappa coefficient between all pairs of feature sets.
It can be rewritten as (see Notation)
\deqn{\frac{2}{m (m - 1)} \sum_{i=1}^{m-1} \sum_{j = i+1}^m
\frac{|V_i \cap V_j| - \frac{|V_i| \cdot |V_j|}{p}}
{\frac{|V_i| + |V_j|}{2} - \frac{|V_i| \cdot |V_j|}{p}}.}
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
stabilityKappa(features = feats, p = 10)
}
\references{
Carletta, Jean (1996).
\dQuote{Assessing Agreement on Classification Tasks: The Kappa Statistic.}
\emph{Computational Linguistics}, \bold{22}(2), 249--254.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/documentation_stability.R
\name{correctedDocumentation}
\alias{correctedDocumentation}
\title{Corrected Stability Measures}
\arguments{
\item{p}{\code{numeric(1)}\cr
Total number of features in the datasets.}
}
\description{
Corrected Stability Measures
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_adjusted.R
\encoding{UTF-8}
\name{stabilityIntersectionCount}
\alias{stabilityIntersectionCount}
\title{Stability Measure Adjusted Intersection Count}
\usage{
stabilityIntersectionCount(
  features,
  sim.mat,
  threshold = 0.9,
  correction.for.chance = "estimate",
  N = 10000,
  impute.na = NULL
)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{sim.mat}{\code{numeric matrix}\cr
Similarity matrix which contains the similarity structure of all features based on
all datasets. The similarity values must be in the range of [0, 1] where 0 indicates
very low similarity and 1 indicates very high similarity. If the list elements of
\code{features} are integerish vectors, then the feature numbering must correspond to the
ordering of \code{sim.mat}. If the list elements of \code{features} are character
vectors, then \code{sim.mat} must be named and the names of \code{sim.mat} must correspond
to the entries in \code{features}.}

\item{threshold}{\code{numeric(1)}\cr
Threshold for indicating which features are similar and which are not. Two features
are considered as similar, if and only if the corresponding entry of \code{sim.mat} is greater
than or equal to \code{threshold}.}

\item{correction.for.chance}{\code{character(1)}\cr
How should the expected value of the stability score (see Details)
be assessed? Options are "estimate", "exact" and "none".
For "estimate", \code{N} random feature sets of the same sizes as the input feature
sets (\code{features}) are generated.
For "exact", all possible combinations of feature sets of the same
sizes as the input feature sets are used. Computation is only feasible for very
small numbers of features and numbers of considered datasets (\code{length(features)}).
For "none", the transformation \eqn{(score - expected) / (maximum - expected)}
is not conducted, i.e. only \eqn{score} is used.
This is not recommended.}

\item{N}{\code{numeric(1)}\cr
Number of random feature sets to consider. Only relevant if \code{correction.for.chance}
is set to "estimate".}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as (see Notation)
\deqn{\frac{2}{m(m-1)}\sum_{i=1}^{m-1} \sum_{j=i+1}^{m}
\frac{I(V_i, V_j) - E(I(V_i, V_j))}{\sqrt{|V_i| \cdot |V_j|} - E(I(V_i, V_j))}}
with \deqn{I(V_i, V_j) = |V_i \cap V_j| + \min (C(V_i, V_j), C(V_j, V_i))} and
\deqn{C(V_k, V_l) = |\{x \in  V_k \setminus V_l : \exists y \in
V_l \setminus V_k \ \mathrm{with Similarity} (x,y) \geq \mathrm{threshold} \}|.}
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
mat = 0.92 ^ abs(outer(1:10, 1:10, "-"))
stabilityIntersectionCount(features = feats, sim.mat = mat, N = 1000)
}
\references{
Bommert A, Rahnenführer J (2020).
\dQuote{Adjusted Measures for Feature Selection Stability for Data Sets with Similar Features.}
In \emph{Machine Learning, Optimization, and Data Science}, 203--214.
\doi{10.1007/978-3-030-64583-0_19}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_uncorrected.R
\encoding{UTF-8}
\name{stabilityOchiai}
\alias{stabilityOchiai}
\title{Stability Measure Ochiai}
\usage{
stabilityOchiai(
  features,
  p = NULL,
  correction.for.chance = "none",
  N = 10000,
  impute.na = NULL
)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{p}{\code{numeric(1)}\cr
Total number of features in the datasets.
Required, if \code{correction.for.chance} is set to "estimate" or "exact".}

\item{correction.for.chance}{\code{character(1)}\cr
Should a correction for chance be applied? Correction for chance means that if
features are chosen at random, the expected value must be independent of the number
of chosen features. To correct for chance, the original score is transformed by
\eqn{(score - expected) / (maximum - expected)}. For stability measures whose
score is the average value of pairwise scores, this transformation
is done for all components individually.
Options are "none", "estimate" and "exact".
For "none", no correction is performed, i.e. the original score is used.
For "estimate", \code{N} random feature sets of the same sizes as the input
feature sets (\code{features}) are generated.
For "exact", all possible combinations of feature sets of the same
sizes as the input feature sets are used. Computation is only feasible for very
small numbers of features (\code{p}) and numbers of considered datasets
(\code{length(features)}).}

\item{N}{\code{numeric(1)}\cr
Number of random feature sets to consider. Only relevant if \code{correction.for.chance}
is set to "estimate".}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as (see Notation)
\deqn{\frac{2}{m (m - 1)} \sum_{i=1}^{m-1} \sum_{j = i+1}^m
\frac{|V_i \cap V_j|}{\sqrt{|V_i| \cdot |V_j|}}.}
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
stabilityOchiai(features = feats)
}
\references{
Ochiai A (1957).
\dQuote{Zoogeographical Studies on the Soleoid Fishes Found in Japan and its Neighbouring Regions-III.}
\emph{Nippon Suisan Gakkaishi}, \bold{22}(9), 531-535.
\doi{10.2331/suisan.22.531}.

Bommert A, Rahnenführer J, Lang M (2017).
\dQuote{A Multicriteria Approach to Find Predictive and Sparse Models with Stable Feature Selection for High-Dimensional Data.}
\emph{Computational and Mathematical Methods in Medicine}, \bold{2017}, 1--18.
\doi{10.1155/2017/7907163}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_uncorrected.R
\encoding{UTF-8}
\name{stabilityHamming}
\alias{stabilityHamming}
\title{Stability Measure Hamming}
\usage{
stabilityHamming(
  features,
  p,
  correction.for.chance = "none",
  N = 10000,
  impute.na = NULL
)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{p}{\code{numeric(1)}\cr
Total number of features in the datasets.
Required, if \code{correction.for.chance} is set to "estimate" or "exact".}

\item{correction.for.chance}{\code{character(1)}\cr
Should a correction for chance be applied? Correction for chance means that if
features are chosen at random, the expected value must be independent of the number
of chosen features. To correct for chance, the original score is transformed by
\eqn{(score - expected) / (maximum - expected)}. For stability measures whose
score is the average value of pairwise scores, this transformation
is done for all components individually.
Options are "none", "estimate" and "exact".
For "none", no correction is performed, i.e. the original score is used.
For "estimate", \code{N} random feature sets of the same sizes as the input
feature sets (\code{features}) are generated.
For "exact", all possible combinations of feature sets of the same
sizes as the input feature sets are used. Computation is only feasible for very
small numbers of features (\code{p}) and numbers of considered datasets
(\code{length(features)}).}

\item{N}{\code{numeric(1)}\cr
Number of random feature sets to consider. Only relevant if \code{correction.for.chance}
is set to "estimate".}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as (see Notation)
\deqn{\frac{2}{m (m - 1)} \sum_{i=1}^{m-1} \sum_{j = i+1}^m
\frac{|V_i \cap V_j| + |V_i^c \cap V_j^c|}{p}.}
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
stabilityHamming(features = feats, p = 10)
}
\references{
Dunne, Kevin, Cunningham, Padraig, Azuaje, Francisco (2002).
\dQuote{Solutions to instability problems with sequential wrapper-based approaches to feature selection.}
Machine Learning Group, Department of Computer Science, Trinity College, Dublin.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_corrected.R
\encoding{UTF-8}
\name{stabilityNogueira}
\alias{stabilityNogueira}
\title{Stability Measure Nogueira}
\usage{
stabilityNogueira(features, p, impute.na = NULL)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{p}{\code{numeric(1)}\cr
Total number of features in the datasets.}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as (see Notation)
\deqn{1 - \frac{\frac{1}{p} \sum_{j=1}^p \frac{m}{m-1} \frac{h_j}{m} \left(1 - \frac{h_j}{m}\right)}
{\frac{q}{mp} (1 - \frac{q}{mp})}.}
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
stabilityNogueira(features = feats, p = 10)
}
\references{
Nogueira S, Sechidis K, Brown G (2018).
\dQuote{On the Stability of Feature Selection Algorithms.}
\emph{Journal of Machine Learning Research}, \bold{18}(174), 1--54.
\url{https://jmlr.org/papers/v18/17-514.html}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_uncorrected.R
\encoding{UTF-8}
\name{stabilityDavis}
\alias{stabilityDavis}
\title{Stability Measure Davis}
\usage{
stabilityDavis(
  features,
  p,
  correction.for.chance = "none",
  N = 10000,
  impute.na = NULL,
  penalty = 0
)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{p}{\code{numeric(1)}\cr
Total number of features in the datasets.
Required, if \code{correction.for.chance} is set to "estimate" or "exact".}

\item{correction.for.chance}{\code{character(1)}\cr
Should a correction for chance be applied? Correction for chance means that if
features are chosen at random, the expected value must be independent of the number
of chosen features. To correct for chance, the original score is transformed by
\eqn{(score - expected) / (maximum - expected)}. For stability measures whose
score is the average value of pairwise scores, this transformation
is done for all components individually.
Options are "none", "estimate" and "exact".
For "none", no correction is performed, i.e. the original score is used.
For "estimate", \code{N} random feature sets of the same sizes as the input
feature sets (\code{features}) are generated.
For "exact", all possible combinations of feature sets of the same
sizes as the input feature sets are used. Computation is only feasible for very
small numbers of features (\code{p}) and numbers of considered datasets
(\code{length(features)}).}

\item{N}{\code{numeric(1)}\cr
Number of random feature sets to consider. Only relevant if \code{correction.for.chance}
is set to "estimate".}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}

\item{penalty}{\code{numeric(1)}\cr
Penalty parameter, see Details.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as (see Notation)
\deqn{\max \left\{ 0, \frac{1}{|V|} \sum_{j=1}^p \frac{h_j}{m} - \frac{penalty}{p}
\cdot \mathop{\mathrm{median}} \{ |V_1|, \ldots, |V_m| \}  \right\}.}
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
stabilityDavis(features = feats, p = 10)
}
\references{
Davis CA, Gerick F, Hintermair V, Friedel CC, Fundel K, Kuffner R, Zimmer R (2006).
\dQuote{Reliable gene signatures for microarray classification: assessment of stability and performance.}
\emph{Bioinformatics}, \bold{22}(19), 2356--2363.
\doi{10.1093/bioinformatics/btl400}.

Bommert A, Rahnenführer J, Lang M (2017).
\dQuote{A Multicriteria Approach to Find Predictive and Sparse Models with Stable Feature Selection for High-Dimensional Data.}
\emph{Computational and Mathematical Methods in Medicine}, \bold{2017}, 1--18.
\doi{10.1155/2017/7907163}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_corrected.R
\encoding{UTF-8}
\name{stabilityLustgarten}
\alias{stabilityLustgarten}
\title{Stability Measure Lustgarten}
\usage{
stabilityLustgarten(features, p, impute.na = NULL)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{p}{\code{numeric(1)}\cr
Total number of features in the datasets.}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as (see Notation)
\deqn{\frac{2}{m (m - 1)} \sum_{i=1}^{m-1} \sum_{j = i+1}^m
\frac{|V_i \cap V_j| - \frac{|V_i| \cdot |V_j|}{p}}
{\min \{|V_i|, |V_j|\} - \max \{ 0, |V_i| + |V_j| - p \}}.}
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
stabilityLustgarten(features = feats, p = 10)
}
\references{
Lustgarten, L J, Gopalakrishnan, Vanathi, Visweswaran, Shyam (2009).
\dQuote{Measuring stability of feature selection in biomedical datasets.}
In \emph{AMIA annual symposium proceedings}, volume 2009, 406.
American Medical Informatics Association.

Bommert A, Rahnenführer J, Lang M (2017).
\dQuote{A Multicriteria Approach to Find Predictive and Sparse Models with Stable Feature Selection for High-Dimensional Data.}
\emph{Computational and Mathematical Methods in Medicine}, \bold{2017}, 1--18.
\doi{10.1155/2017/7907163}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/documentation_stability.R
\name{stabilityDocumentation}
\alias{stabilityDocumentation}
\title{Stability of Feature Selection}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{penalty}{\code{numeric(1)}\cr
Penalty parameter, see Details.}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}

\item{N}{\code{numeric(1)}\cr
Number of random feature sets to consider. Only relevant if \code{correction.for.chance}
is set to "estimate".}

\item{sim.mat}{\code{numeric matrix}\cr
Similarity matrix which contains the similarity structure of all features based on
all datasets. The similarity values must be in the range of [0, 1] where 0 indicates
very low similarity and 1 indicates very high similarity. If the list elements of
\code{features} are integerish vectors, then the feature numbering must correspond to the
ordering of \code{sim.mat}. If the list elements of \code{features} are character
vectors, then \code{sim.mat} must be named and the names of \code{sim.mat} must correspond
to the entries in \code{features}.}

\item{threshold}{\code{numeric(1)}\cr
Threshold for indicating which features are similar and which are not. Two features
are considered as similar, if and only if the corresponding entry of \code{sim.mat} is greater
than or equal to \code{threshold}.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\seealso{
\link{listStabilityMeasures}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_uncorrected.R
\encoding{UTF-8}
\name{stabilityJaccard}
\alias{stabilityJaccard}
\title{Stability Measure Jaccard}
\usage{
stabilityJaccard(
  features,
  p = NULL,
  correction.for.chance = "none",
  N = 10000,
  impute.na = NULL
)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{p}{\code{numeric(1)}\cr
Total number of features in the datasets.
Required, if \code{correction.for.chance} is set to "estimate" or "exact".}

\item{correction.for.chance}{\code{character(1)}\cr
Should a correction for chance be applied? Correction for chance means that if
features are chosen at random, the expected value must be independent of the number
of chosen features. To correct for chance, the original score is transformed by
\eqn{(score - expected) / (maximum - expected)}. For stability measures whose
score is the average value of pairwise scores, this transformation
is done for all components individually.
Options are "none", "estimate" and "exact".
For "none", no correction is performed, i.e. the original score is used.
For "estimate", \code{N} random feature sets of the same sizes as the input
feature sets (\code{features}) are generated.
For "exact", all possible combinations of feature sets of the same
sizes as the input feature sets are used. Computation is only feasible for very
small numbers of features (\code{p}) and numbers of considered datasets
(\code{length(features)}).}

\item{N}{\code{numeric(1)}\cr
Number of random feature sets to consider. Only relevant if \code{correction.for.chance}
is set to "estimate".}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as (see Notation)
\deqn{\frac{2}{m (m - 1)} \sum_{i=1}^{m-1} \sum_{j = i+1}^m
\frac{|V_i \cap V_j|}{|V_i \cup V_j|}.}
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
stabilityJaccard(features = feats)
}
\references{
Jaccard, Paul (1901).
\dQuote{Étude comparative de la distribution florale dans une portion des Alpes et du Jura.}
\emph{Bulletin de la Société Vaudoise des Sciences Naturelles}, \bold{37}, 547-579.
\doi{10.5169/SEALS-266450}.

Bommert A, Rahnenführer J, Lang M (2017).
\dQuote{A Multicriteria Approach to Find Predictive and Sparse Models with Stable Feature Selection for High-Dimensional Data.}
\emph{Computational and Mathematical Methods in Medicine}, \bold{2017}, 1--18.
\doi{10.1155/2017/7907163}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/list_stability_measures.R
\name{listStabilityMeasures}
\alias{listStabilityMeasures}
\title{List All Available Stability Measures}
\usage{
listStabilityMeasures()
}
\value{
\code{data.frame} \cr
For each stability measure, its name,
the information, whether it is corrected for chance by definition,
the information, whether it is adjusted for similar features,
its minimal value and its maximal value are displayed.
}
\description{
Lists all stability measures of package \emph{stabm} and
provides information about them.
}
\section{Note}{
 The given minimal values might only be reachable
in some scenarios, e.g. if the feature sets have a certain size.
The measures which are not corrected for chance by definition can
be corrected for chance with \code{correction.for.chance}.
This however changes the minimal value.
For the adjusted stability measures, the minimal value depends
on the similarity structure.
}

\examples{
listStabilityMeasures()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_features.R
\name{plotFeatures}
\alias{plotFeatures}
\title{Plot Selected Features}
\usage{
plotFeatures(features, sim.mat = NULL)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{sim.mat}{\code{numeric matrix}\cr
Similarity matrix which contains the similarity structure of all features based on
all datasets. The similarity values must be in the range of [0, 1] where 0 indicates
very low similarity and 1 indicates very high similarity. If the list elements of
\code{features} are integerish vectors, then the feature numbering must correspond to the
ordering of \code{sim.mat}. If the list elements of \code{features} are character
vectors, then \code{sim.mat} must be named and the names of \code{sim.mat} must correspond
to the entries in \code{features}.}
}
\value{
Object of class \code{ggplot}.
}
\description{
Creates a heatmap of the features which are selected in at least one feature set.
The sets are ordered according to average linkage hierarchical clustering based on the Manhattan
distance. If \code{sim.mat} is given, the features are ordered according to average linkage
hierarchical clustering based on \code{1 - sim.mat}. Otherwise, the features are ordered in
the same way as the feature sets.

Note that this function needs the packages \CRANpkg{ggplot2}, \CRANpkg{cowplot} and
\CRANpkg{ggdendro} installed.
}
\examples{
feats = list(1:3, 1:4, 1:5)
mat = 0.92 ^ abs(outer(1:10, 1:10, "-"))
plotFeatures(features = feats)
plotFeatures(features = feats, sim.mat = mat)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_corrected.R
\encoding{UTF-8}
\name{stabilityUnadjusted}
\alias{stabilityUnadjusted}
\title{Stability Measure Unadjusted}
\usage{
stabilityUnadjusted(features, p, impute.na = NULL)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{p}{\code{numeric(1)}\cr
Total number of features in the datasets.}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as (see Notation)
\deqn{\frac{2}{m (m - 1)} \sum_{i=1}^{m-1} \sum_{j = i+1}^m
\frac{|V_i \cap V_j| - \frac{|V_i| \cdot |V_j|}{p}}
{\sqrt{|V_i| \cdot |V_j|} - \frac{|V_i| \cdot |V_j|}{p}}.}
This is what \link{stabilityIntersectionMBM}, \link{stabilityIntersectionGreedy},
\link{stabilityIntersectionCount} and \link{stabilityIntersectionMean}
become, when there are no similar features.
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
stabilityUnadjusted(features = feats, p = 10)
}
\references{
Bommert A, Rahnenführer J (2020).
\dQuote{Adjusted Measures for Feature Selection Stability for Data Sets with Similar Features.}
In \emph{Machine Learning, Optimization, and Data Science}, 203--214.
\doi{10.1007/978-3-030-64583-0_19}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_corrected.R
\encoding{UTF-8}
\name{stabilitySomol}
\alias{stabilitySomol}
\title{Stability Measure Somol}
\usage{
stabilitySomol(features, p, impute.na = NULL)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{p}{\code{numeric(1)}\cr
Total number of features in the datasets.}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as (see Notation)
\deqn{\frac{\left(\sum\limits_{j=1}^p \frac{h_j}{q} \frac{h_j - 1}{m-1}\right) -
c_{\min}}{c_{\max} - c_{\min}}} with
\deqn{c_{\min} = \frac{q^2 - p(q - q \ \mathop{mod} \ p) - \left(q \ \mathop{mod} \ p\right)^2}{p q (m-1)},}
\deqn{c_{\max} = \frac{\left(q \ \mathop{mod} \ m\right)^2 + q(m-1) - \left(q \ \mathop{mod} \ m\right)m}{q(m-1)}.}
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
stabilitySomol(features = feats, p = 10)
}
\references{
Somol P, Novovičová J (2010).
\dQuote{Evaluating Stability and Comparing Output of Feature Selectors that Optimize Feature Subset Cardinality.}
\emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, \bold{32}(11), 1921--1939.
\doi{10.1109/tpami.2010.34}.

Bommert A, Rahnenführer J, Lang M (2017).
\dQuote{A Multicriteria Approach to Find Predictive and Sparse Models with Stable Feature Selection for High-Dimensional Data.}
\emph{Computational and Mathematical Methods in Medicine}, \bold{2017}, 1--18.
\doi{10.1155/2017/7907163}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_adjusted.R
\encoding{UTF-8}
\name{stabilityIntersectionGreedy}
\alias{stabilityIntersectionGreedy}
\title{Stability Measure Adjusted Intersection Greedy}
\usage{
stabilityIntersectionGreedy(
  features,
  sim.mat,
  threshold = 0.9,
  correction.for.chance = "estimate",
  N = 10000,
  impute.na = NULL
)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{sim.mat}{\code{numeric matrix}\cr
Similarity matrix which contains the similarity structure of all features based on
all datasets. The similarity values must be in the range of [0, 1] where 0 indicates
very low similarity and 1 indicates very high similarity. If the list elements of
\code{features} are integerish vectors, then the feature numbering must correspond to the
ordering of \code{sim.mat}. If the list elements of \code{features} are character
vectors, then \code{sim.mat} must be named and the names of \code{sim.mat} must correspond
to the entries in \code{features}.}

\item{threshold}{\code{numeric(1)}\cr
Threshold for indicating which features are similar and which are not. Two features
are considered as similar, if and only if the corresponding entry of \code{sim.mat} is greater
than or equal to \code{threshold}.}

\item{correction.for.chance}{\code{character(1)}\cr
How should the expected value of the stability score (see Details)
be assessed? Options are "estimate", "exact" and "none".
For "estimate", \code{N} random feature sets of the same sizes as the input feature
sets (\code{features}) are generated.
For "exact", all possible combinations of feature sets of the same
sizes as the input feature sets are used. Computation is only feasible for very
small numbers of features and numbers of considered datasets (\code{length(features)}).
For "none", the transformation \eqn{(score - expected) / (maximum - expected)}
is not conducted, i.e. only \eqn{score} is used.
This is not recommended.}

\item{N}{\code{numeric(1)}\cr
Number of random feature sets to consider. Only relevant if \code{correction.for.chance}
is set to "estimate".}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as (see Notation)
\deqn{\frac{2}{m(m-1)}\sum_{i=1}^{m-1} \sum_{j=i+1}^{m}
\frac{I(V_i, V_j) - E(I(V_i, V_j))}{\sqrt{|V_i| \cdot |V_j|} - E(I(V_i, V_j))}} with
\deqn{I(V_i, V_j) = |V_i \cap V_j| + \mathop{\mathrm{GMBM}}(V_i \setminus V_j, V_j \backslash V_i).}
\eqn{\mathop{\mathrm{GMBM}}(V_i \setminus V_j, V_j \backslash V_i)} denotes a greedy approximation
of \eqn{\mathop{\mathrm{MBM}}(V_i \setminus V_j, V_j \backslash V_i)}, see \link{stabilityIntersectionMBM}.
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
mat = 0.92 ^ abs(outer(1:10, 1:10, "-"))
stabilityIntersectionGreedy(features = feats, sim.mat = mat, N = 1000)
}
\references{
Bommert A, Rahnenführer J (2020).
\dQuote{Adjusted Measures for Feature Selection Stability for Data Sets with Similar Features.}
In \emph{Machine Learning, Optimization, and Data Science}, 203--214.
\doi{10.1007/978-3-030-64583-0_19}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_adjusted.R
\encoding{UTF-8}
\name{stabilityZucknick}
\alias{stabilityZucknick}
\title{Stability Measure Zucknick}
\usage{
stabilityZucknick(
  features,
  sim.mat,
  threshold = 0.9,
  correction.for.chance = "none",
  N = 10000,
  impute.na = NULL
)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{sim.mat}{\code{numeric matrix}\cr
Similarity matrix which contains the similarity structure of all features based on
all datasets. The similarity values must be in the range of [0, 1] where 0 indicates
very low similarity and 1 indicates very high similarity. If the list elements of
\code{features} are integerish vectors, then the feature numbering must correspond to the
ordering of \code{sim.mat}. If the list elements of \code{features} are character
vectors, then \code{sim.mat} must be named and the names of \code{sim.mat} must correspond
to the entries in \code{features}.}

\item{threshold}{\code{numeric(1)}\cr
Threshold for indicating which features are similar and which are not. Two features
are considered as similar, if and only if the corresponding entry of \code{sim.mat} is greater
than or equal to \code{threshold}.}

\item{correction.for.chance}{\code{character(1)}\cr
Should a correction for chance be applied? Correction for chance means that if
features are chosen at random, the expected value must be independent of the number
of chosen features. To correct for chance, the original score is transformed by
\eqn{(score - expected) / (maximum - expected)}. For stability measures whose
score is the average value of pairwise scores, this transformation
is done for all components individually.
Options are "none", "estimate" and "exact".
For "none", no correction is performed, i.e. the original score is used.
For "estimate", \code{N} random feature sets of the same sizes as the input
feature sets (\code{features}) are generated.
For "exact", all possible combinations of feature sets of the same
sizes as the input feature sets are used. Computation is only feasible for very
small numbers of features (\code{p}) and numbers of considered datasets
(\code{length(features)}).}

\item{N}{\code{numeric(1)}\cr
Number of random feature sets to consider. Only relevant if \code{correction.for.chance}
is set to "estimate".}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as
\deqn{\frac{2}{m(m-1)}\sum_{i=1}^{m-1} \sum_{j=i+1}^{m}
\frac{|V_i \cap V_j| + C(V_i, V_j) + C(V_j, V_i)}{|V_i \cup V_j|}} with
\deqn{C(V_k, V_l) = \frac{1}{|V_l|} \sum_{(x, y) \in V_k \times (V_l \setminus V_k) \ \mathrm{with Similarity}(x,y) \geq \mathrm{threshold}} \mathop{\mathrm{Similarity}}(x,y).}
Note that this definition slightly differs from its original in order to make it suitable
for arbitrary similarity measures.
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
mat = 0.92 ^ abs(outer(1:10, 1:10, "-"))
stabilityZucknick(features = feats, sim.mat = mat)
}
\references{
Zucknick M, Richardson S, Stronach EA (2008).
\dQuote{Comparing the Characteristics of Gene Expression Profiles Derived by Univariate and Multivariate Classification Methods.}
\emph{Statistical Applications in Genetics and Molecular Biology}, \bold{7}(1).
\doi{10.2202/1544-6115.1307}.

Bommert A, Rahnenführer J, Lang M (2017).
\dQuote{A Multicriteria Approach to Find Predictive and Sparse Models with Stable Feature Selection for High-Dimensional Data.}
\emph{Computational and Mathematical Methods in Medicine}, \bold{2017}, 1--18.
\doi{10.1155/2017/7907163}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_functions_adjusted.R
\encoding{UTF-8}
\name{stabilityIntersectionMean}
\alias{stabilityIntersectionMean}
\title{Stability Measure Adjusted Intersection Mean}
\usage{
stabilityIntersectionMean(
  features,
  sim.mat,
  threshold = 0.9,
  correction.for.chance = "estimate",
  N = 10000,
  impute.na = NULL
)
}
\arguments{
\item{features}{\code{list (length >= 2)}\cr
Chosen features per dataset. Each element of the list contains the features for one dataset.
The features must be given by their names (\code{character}) or indices (\code{integerish}).}

\item{sim.mat}{\code{numeric matrix}\cr
Similarity matrix which contains the similarity structure of all features based on
all datasets. The similarity values must be in the range of [0, 1] where 0 indicates
very low similarity and 1 indicates very high similarity. If the list elements of
\code{features} are integerish vectors, then the feature numbering must correspond to the
ordering of \code{sim.mat}. If the list elements of \code{features} are character
vectors, then \code{sim.mat} must be named and the names of \code{sim.mat} must correspond
to the entries in \code{features}.}

\item{threshold}{\code{numeric(1)}\cr
Threshold for indicating which features are similar and which are not. Two features
are considered as similar, if and only if the corresponding entry of \code{sim.mat} is greater
than or equal to \code{threshold}.}

\item{correction.for.chance}{\code{character(1)}\cr
How should the expected value of the stability score (see Details)
be assessed? Options are "estimate", "exact" and "none".
For "estimate", \code{N} random feature sets of the same sizes as the input feature
sets (\code{features}) are generated.
For "exact", all possible combinations of feature sets of the same
sizes as the input feature sets are used. Computation is only feasible for very
small numbers of features and numbers of considered datasets (\code{length(features)}).
For "none", the transformation \eqn{(score - expected) / (maximum - expected)}
is not conducted, i.e. only \eqn{score} is used.
This is not recommended.}

\item{N}{\code{numeric(1)}\cr
Number of random feature sets to consider. Only relevant if \code{correction.for.chance}
is set to "estimate".}

\item{impute.na}{\code{numeric(1)}\cr
In some scenarios, the stability cannot be assessed based on all feature sets.
E.g. if some of the feature sets are empty, the respective pairwise comparisons yield NA as result.
With which value should these missing values be imputed? \code{NULL} means no imputation.}
}
\value{
\code{numeric(1)} Stability value.
}
\description{
The stability of feature selection is defined as the robustness of
the sets of selected features with respect to small variations in the data on which the
feature selection is conducted. To quantify stability, several datasets from the
same data generating process can be used. Alternatively, a single dataset can be
split into parts by resampling. Either way, all datasets used for feature selection must
contain exactly the same features. The feature selection method of interest is
applied on all of the datasets and the sets of chosen features are recorded.
The stability of the feature selection is assessed based on the sets of chosen features
using stability measures.
}
\details{
The stability measure is defined as (see Notation)
\deqn{\frac{2}{m(m-1)}\sum_{i=1}^{m-1} \sum_{j=i+1}^{m}
\frac{I(V_i, V_j) - E(I(V_i, V_j))}{\sqrt{|V_i| \cdot |V_j|} - E(I(V_i, V_j))}}
with \deqn{I(V_i, V_j) = |V_i \cap V_j| + \min (C(V_i, V_j), C(V_j, V_i)),}
\deqn{C(V_k, V_l) = \sum_{x \in V_k \setminus V_l : |G^{kl}_x| > 0}
\frac{1}{|G^{kl}_x|} \sum_{y \in G^{kl}_x} \ \mathrm{Similarity} (x,y)} and
\deqn{G^{kl}_x = \{y \in V_l \setminus V_k: \ \mathrm{Similarity} (x, y) \geq \mathrm{threshold} \}.}
}
\section{Notation}{
 For the definition of all stability measures in this package,
the following notation is used:
Let \eqn{V_1, \ldots, V_m} denote the sets of chosen features
for the \eqn{m} datasets, i.e. \code{features} has length \eqn{m} and
\eqn{V_i} is a set which contains the \eqn{i}-th entry of \code{features}.
Furthermore, let \eqn{h_j} denote the number of sets that contain feature
\eqn{X_j} so that \eqn{h_j} is the absolute frequency with which feature \eqn{X_j}
is chosen.
Analogously, let \eqn{h_{ij}} denote the number of sets that include both \eqn{X_i} and \eqn{X_j}.
Also, let \eqn{q = \sum_{j=1}^p h_j = \sum_{i=1}^m |V_i|} and \eqn{V = \bigcup_{i=1}^m V_i}.
}

\examples{
feats = list(1:3, 1:4, 1:5)
mat = 0.92 ^ abs(outer(1:10, 1:10, "-"))
stabilityIntersectionMean(features = feats, sim.mat = mat, N = 1000)
}
\references{
Bommert A, Rahnenführer J (2020).
\dQuote{Adjusted Measures for Feature Selection Stability for Data Sets with Similar Features.}
In \emph{Machine Learning, Optimization, and Data Science}, 203--214.
\doi{10.1007/978-3-030-64583-0_19}.

Bommert A (2020).
\emph{Integration of Feature Selection Stability in Model Fitting}.
Ph.D. thesis, TU Dortmund University, Germany.
\doi{10.17877/DE290R-21906}.
}
\seealso{
\link{listStabilityMeasures}
}
