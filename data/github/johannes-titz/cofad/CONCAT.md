Cofad User Guide
================

[![Build
Status](https://travis-ci.com/johannes-titz/cofad.svg?branch=main)](https://travis-ci.com/johannes-titz/cofad)
[![CRAN
status](https://www.r-pkg.org/badges/version/cofad)](https://CRAN.R-project.org/package=cofad)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03822/status.svg)](https://doi.org/10.21105/joss.03822)

<!-- [![DOI](https://joss.theoj.org/papers/10.21105/joss.02116/status.svg)](https://doi.org/10.21105/joss.02116) -->
<!-- To cite cofad in publications use: -->
<!-- Titz, J. (2020). cofad: A modern graphical user interface for 2-level mixed models. *Journal of Open Source Software, 5*(49), 2116. https://doi.org/10.21105/joss.02116 -->
<!-- A BibTeX entry for LaTeX users is -->
<!-- ``` -->
<!-- @article{titz2020, -->
<!--   title = {cofad: {A} }, -->
<!--   author = {Titz, Johannes and Burkhardt, Markus}, -->
<!--   year = {2020}, -->
<!--   journal = {Journal of Open Source Software}, -->
<!--   volume = {5}, -->
<!--   pages = {2116}, -->
<!--   doi = {10.21105/joss.02116}, -->
<!--   number = {49} -->
<!-- } -->
<!-- ``` -->

## Introduction

Cofad is an R package for conducting COntrast analysis in FActorial
Designs like ANOVAs. If contrast analysis was to win a price it would be
the one for the most underestimated, underused statistical technique.
This is unfortunate because in every case a contrast analysis is at
least as good as an ANOVA, but in most cases it is better. Contrast
analysis gets rid off the unspecific omnibus-hypothesis *there are
differences somewhere* and replaces it with a very specific numerical
hypothesis. Furthermore, contrast analysis focuses on effects instead of
significance. This is expressed doubly: First, there are three different
effect sizes for contrast analysis:
![r\_\\mathrm{effectsize}](https://latex.codecogs.com/png.latex?r_%5Cmathrm%7Beffectsize%7D "r_\mathrm{effectsize}"),
![r\_\\mathrm{contrast}](https://latex.codecogs.com/png.latex?r_%5Cmathrm%7Bcontrast%7D "r_\mathrm{contrast}")
and
![r\_\\mathrm{alerting}](https://latex.codecogs.com/png.latex?r_%5Cmathrm%7Balerting%7D "r_\mathrm{alerting}").
Second, the effect size refers not to the data but to the tested
hypothesis. The larger the effect, the more this speaks for the
hypothesis. One can even compare different hypotheses against each other
(experimentum crucis!) by looking at the effect size for each
hypothesis.

Sounds interesting? Then check out some introductory literature such as
Furr (2004), Rosenthal & Rosnow (1985), Rosenthal, Rosnow, & Rubin
(2000), or, for the German-speaking audience, Sedlmeier & Renkewitz
(2018). Contrast analysis is fairly easy to understand if you know what
an ANOVA and a correlation is. In this vignette we assume you are
familiar with the basics of contrast analysis and want to apply it to a
specific data set. First we show how to install cofad and use the
graphical user interface. Then we demonstrate some exemplary analyses
for between, within and mixed designs in R.

## Installation

Cofad has two components, the plain R package and a shiny-app that
offers an intuitive graphical user interface.

If you just want to use the cofad-app, you do not need to install it.
Just go to <https://cofad.titz.science> and use it there. An example
data file is loaded when you go to <https://cofad.titz.science/example>.

If you prefer the command line interface or want to use the cofad-app
locally, install it from github (you need the package devtools for
this):

``` r
# install.packages("devtools") # uncomment if you do not have devtools installed
devtools::install_github("johannes-titz/cofad")
```

Now you can load cofad and use it in your R scripts.

You can also run the app:

``` r
cofad::run_app()
```

If you have any problems installing cofad, check that your R version is
up to date (currently R version 4.1.2 (2021-11-01)). If you are using
Windows, enable TLS 1.2 in the Internet Options Advanced tab (see
<https://github.com/r-lib/remotes/issues/130#issuecomment-423830669>).
Under Windows, you will also need Rtools to build the package:
<https://cran.r-project.org/bin/windows/Rtools/>.

If it still does not work drop an e-mail at johannes at titz.science or
at johannes.titz at gmail.com.

## Using cofad

Before we start: Your data has to be in the long-format (also referred
to as narrow or tidy)! If you do not know what this means, please check
the short description of the Wikipedia-article:
<https://en.wikipedia.org/wiki/Wide_and_narrow_data>

### Graphical-User-Interface

The graphical-user-interface is self-explanatory. Just load your data
and drag the variables to the correct position. At the moment you can
only read .sav (SPSS) and .csv files.

As an example go to <https://cofad.titz.science/example> which will load
a data set from Rosenthal et al. (2000) (Table 5.3). The cognitive
ability of nine children belonging to different age groups (between) was
measured four times (within).

There are two hypotheses:

1.  cognitive ability linearly increases over time (within)
    (![\\lambda\_\\mathrm{1} = -3, \\lambda\_\\mathrm{2} = -1, \\lambda\_\\mathrm{3} = 1, \\lambda\_\\mathrm{4} = 3](https://latex.codecogs.com/png.latex?%5Clambda_%5Cmathrm%7B1%7D%20%3D%20-3%2C%20%5Clambda_%5Cmathrm%7B2%7D%20%3D%20-1%2C%20%5Clambda_%5Cmathrm%7B3%7D%20%3D%201%2C%20%5Clambda_%5Cmathrm%7B4%7D%20%3D%203 "\lambda_\mathrm{1} = -3, \lambda_\mathrm{2} = -1, \lambda_\mathrm{3} = 1, \lambda_\mathrm{4} = 3"))
2.  cognitive ability linearly increase over age groups (between)
    (![\\lambda\_\\mathrm{Age 8} = -1, \\lambda\_\\mathrm{Age 10} = 0, \\lambda\_\\mathrm{Age12} = 1](https://latex.codecogs.com/png.latex?%5Clambda_%5Cmathrm%7BAge%208%7D%20%3D%20-1%2C%20%5Clambda_%5Cmathrm%7BAge%2010%7D%20%3D%200%2C%20%5Clambda_%5Cmathrm%7BAge12%7D%20%3D%201 "\lambda_\mathrm{Age 8} = -1, \lambda_\mathrm{Age 10} = 0, \lambda_\mathrm{Age12} = 1"))

Now drag the variables to the correct position and set the lambdas
accordingly:

![cofad GUI](gui1b.png)

The result should look like this:

![cofad GUI](gui2b.png)

A mixed design is ideal for testing out the cofad-app. You can now
construct a separate within-model by removing the between variable
“age”. Then you can construct a separate between-model by removing
“time” from within and dragging “age” back into the between panel.

The graphical user interface will suffice for most users, but some will
prefer to use the scripting capabilities of R. In the next sections we
will look at several script examples for different designs.

### Between-Subjects Designs

Let us first load the package:

``` r
library(cofad)
```

Now we need some data and hypotheses. We can simply take the data from
Furr (2004), where we have different empathy ratings of students from
different majors. This data set is available in the cofad package:

``` r
data("furr_p4")
furr_p4
#>    empathy      major
#> 1       51 psychology
#> 2       56 psychology
#> 3       61 psychology
#> 4       58 psychology
#> 5       54 psychology
#> 6       62  education
#> 7       67  education
#> 8       57  education
#> 9       65  education
#> 10      59  education
#> 11      50   business
#> 12      49   business
#> 13      47   business
#> 14      45   business
#> 15      44   business
#> 16      50  chemistry
#> 17      45  chemistry
#> 18      40  chemistry
#> 19      49  chemistry
#> 20      41  chemistry
```

Furr states three hypotheses:

-   Contrast A: Psychology majors have higher empathy scores than
    Education majors
    (![\\lambda\_\\mathrm{psych} = 1, \\lambda\_\\mathrm{edu} = -1](https://latex.codecogs.com/png.latex?%5Clambda_%5Cmathrm%7Bpsych%7D%20%3D%201%2C%20%5Clambda_%5Cmathrm%7Bedu%7D%20%3D%20-1 "\lambda_\mathrm{psych} = 1, \lambda_\mathrm{edu} = -1")).
-   Contrast B: Business majors have higher empathy scores than
    Chemistry majors
    (![\\lambda\_\\mathrm{bus} = 1, \\lambda\_\\mathrm{chem} = -1](https://latex.codecogs.com/png.latex?%5Clambda_%5Cmathrm%7Bbus%7D%20%3D%201%2C%20%5Clambda_%5Cmathrm%7Bchem%7D%20%3D%20-1 "\lambda_\mathrm{bus} = 1, \lambda_\mathrm{chem} = -1")).
-   Contrast C: On average, Psychology and Education majors have higher
    empathy scores than Business and Chemistry majors
    (![\\lambda\_\\mathrm{psych} = 1, \\lambda\_\\mathrm{edu} = 1, \\lambda\_\\mathrm{bus} = -1, \\lambda\_\\mathrm{chem} = -1](https://latex.codecogs.com/png.latex?%5Clambda_%5Cmathrm%7Bpsych%7D%20%3D%201%2C%20%5Clambda_%5Cmathrm%7Bedu%7D%20%3D%201%2C%20%5Clambda_%5Cmathrm%7Bbus%7D%20%3D%20-1%2C%20%5Clambda_%5Cmathrm%7Bchem%7D%20%3D%20-1 "\lambda_\mathrm{psych} = 1, \lambda_\mathrm{edu} = 1, \lambda_\mathrm{bus} = -1, \lambda_\mathrm{chem} = -1")).

These hypotheses are only mean comparisons, but this is a good way to
start. Let’s use cofad to conduct the contrast analysis:

``` r
ca <- calc_contrast(dv = empathy, between = major,
                    lambda_between = c("psychology" = 1, "education" = -1,
                                       "business" = 0, "chemistry" = 0),
                    data = furr_p4)
ca
#> 
#> We ran a contrast analysis for the following between contrasts: business = 0; chemistry = 0; education = -1; psychology = 1. This resulted in statistics of F(1,16) = 6.154; p = 0.02461 and an effect magnitude of r_effectsize = -0.276. Attention: Contrast fits in the opposite direction!
```

The print method shows some basic information that can be directly used
in a publication. With the summary method some more details are shown:

``` r
summary(ca)
#> $`F-Table`
#>            SS df     MS     F     p
#> contrast   90  1 90.000 6.154 0.025
#> within    234 16 14.625    NA    NA
#> total    1179 19     NA    NA    NA
#> 
#> $Effects
#>              effects
#> r_effectsize  -0.276
#> r_contrast    -0.527
#> r_alerting    -0.309
```

From this table,
![r\_\\mathrm{effectsize}](https://latex.codecogs.com/png.latex?r_%5Cmathrm%7Beffectsize%7D "r_\mathrm{effectsize}")
is probably the most useful statistic. It is just the correlation
between the lambdas and the dependent variable, which can also be
calculated by hand:

``` r
lambdas <- rep(c(1, -1, 0, 0), each = 5)
cor(furr_p4$empathy, lambdas)
#> [1] -0.2762895
```

The other two hypotheses can be tested accordingly:

``` r
ca <- calc_contrast(dv = empathy, between = major,
                    lambda_between = c("psychology" = 0, "education" = 0,
                                       "business" = 1, "chemistry" = -1),
                    data = furr_p4)
ca
#> 
#> We ran a contrast analysis for the following between contrasts: business = 1; chemistry = -1; education = 0; psychology = 0. This resulted in statistics of F(1,16) = 0.684; p = 0.4205 and an effect magnitude of r_effectsize = 0.092.
ca <- calc_contrast(dv = empathy, between = major,
                    lambda_between = c("psychology" = 1, "education" = 1,
                                       "business" = -1, "chemistry" = -1),
                    data = furr_p4)
ca
#> 
#> We ran a contrast analysis for the following between contrasts: business = -1; chemistry = -1; education = 1; psychology = 1. This resulted in statistics of F(1,16) = 57.778; p = 1.07e-06 and an effect magnitude of r_effectsize = 0.847.
```

You will find that the numbers are identical to the ones presented in
Furr (2004). Now, imagine we have a more fun hypothesis and not just
mean differences. From an elaborate theory we could derive that the
means should be 73, 61, 51 and 38. We can test this with cofad directly
because cofad will center the lambdas (the mean of the lambdas has to be
0):

``` r
ca <- calc_contrast(dv = empathy, between = major,
                    lambda_between = c("psychology" = 73, "education" = 61,
                                       "business" = 51, "chemistry" = 38),
                    data = furr_p4)
#> Warning in check_lambda(lambda_between): lambdas are centered and rounded to 3
#> digits
ca
#> 
#> We ran a contrast analysis for the following between contrasts: business = -4.75; chemistry = -17.75; education = 5.25; psychology = 17.25. This resulted in statistics of F(1,16) = 37.466; p = 1.475e-05 and an effect magnitude of r_effectsize = 0.682.
```

The manual test gives the same effect size:

``` r
lambdas <- rep(c(73, 61, 51, 38), each = 5)
cor(furr_p4$empathy, lambdas)
#> [1] 0.6817294
```

Let us now run an analysis for within-subjects designs.

## Within-Subjects Designs

For within designs the calculations are quite different, but cofad takes
care of the details. We just have to use the within parameters *within*
and *lambda_within* instead of the between equivalents. As an example we
use Table 16.5 from Sedlmeier & Renkewitz (2018). Reading ability was
assessed for eight participants under four different conditions. The
hypothesis is that you can read best without music, white noise reduces
your reading ability and music (independent of type) reduces it even
further.

``` r
data("sedlmeier_p537")
head(sedlmeier_p537)
#>   reading_test participant         music
#> 1           27           1 without music
#> 2           25           2 without music
#> 3           30           3 without music
#> 4           29           4 without music
#> 5           30           5 without music
#> 6           33           6 without music
calc_contrast(dv = reading_test, within = music,
              lambda_within = c("without music" = 1.25, 
                                "white noise" = 0.25,
                                "classic" = -0.75,
                                "jazz" = -0.75),
              id = participant, data = sedlmeier_p537)
#> 
#> We ran a contrast analysis for the following within contrasts: classic = -0.75; jazz = -0.75; white noise = 0.25; without music = 1.25. This resulted in statistics of t(7) = 5.269; p = 0.000581 and an effect magnitude of g_effectsize = 1.863.
```

You can see that the significance test is just a
![t](https://latex.codecogs.com/png.latex?t "t")-test and the reported
effect size is referring to a mean comparison
(![g](https://latex.codecogs.com/png.latex?g "g")). (The
![t](https://latex.codecogs.com/png.latex?t "t")-test is one-tailed,
because contrast analysis has always a specific hypotheses.) When
conducting the analysis by hand, we can see why:

``` r
mtr <- matrix(sedlmeier_p537$reading_test, ncol = 4)
lambdas <- c(1.25, 0.25, -0.75, -0.75)
lc1 <- mtr %*% lambdas
t.test(lc1)
#> 
#>  One Sample t-test
#> 
#> data:  lc1
#> t = 5.2689, df = 7, p-value = 0.001162
#> alternative hypothesis: true mean is not equal to 0
#> 95 percent confidence interval:
#>  3.238361 8.511639
#> sample estimates:
#> mean of x 
#>     5.875
```

Only the linear combination of the dependent variable and the contrast
weights for each participant is needed. With these values a normal
![t](https://latex.codecogs.com/png.latex?t "t")-test against 0 is
conducted. While you can do this manually, using cofad is quicker and it
also gives you more information such as the different effect sizes.

## Mixed Designs

A mixed design combines between and within factors. In this case cofad
first calculates the linear combination (*L*-Values) for the within
factor. This new variable serves as the dependent variable for a between
contrast analysis. We will again look at the example presented in
Rosenthal et al. (2000) (see the section graphical user interface). The
cognitive ability of nine children belonging to different age groups
(between) was measured four times (within).

There are two hypotheses:

1.  cognitive ability linearly increases over time (within)
    (![\\lambda\_\\mathrm{1} = -3, \\lambda\_\\mathrm{2} = -1, \\lambda\_\\mathrm{3} = 1, \\lambda\_\\mathrm{4} = 3](https://latex.codecogs.com/png.latex?%5Clambda_%5Cmathrm%7B1%7D%20%3D%20-3%2C%20%5Clambda_%5Cmathrm%7B2%7D%20%3D%20-1%2C%20%5Clambda_%5Cmathrm%7B3%7D%20%3D%201%2C%20%5Clambda_%5Cmathrm%7B4%7D%20%3D%203 "\lambda_\mathrm{1} = -3, \lambda_\mathrm{2} = -1, \lambda_\mathrm{3} = 1, \lambda_\mathrm{4} = 3"))
2.  cognitive ability linearly increase over age groups (between)
    (![\\lambda\_\\mathrm{Age 8} = -1, \\lambda\_\\mathrm{Age 10} = 0, \\lambda\_\\mathrm{Age12} = 1](https://latex.codecogs.com/png.latex?%5Clambda_%5Cmathrm%7BAge%208%7D%20%3D%20-1%2C%20%5Clambda_%5Cmathrm%7BAge%2010%7D%20%3D%200%2C%20%5Clambda_%5Cmathrm%7BAge12%7D%20%3D%201 "\lambda_\mathrm{Age 8} = -1, \lambda_\mathrm{Age 10} = 0, \lambda_\mathrm{Age12} = 1"))

Let’s have a look at the data and calculation:

``` r
data("rosenthal_tbl53")
head(rosenthal_tbl53)
#>   dv between id within
#> 1  3    age8  1      1
#> 2  1    age8  2      1
#> 3  4    age8  3      1
#> 4  4   age10  4      1
#> 5  5   age10  5      1
#> 6  5   age10  6      1
lambda_within <- c("1" = -3, "2" = -1, "3" = 1, "4" = 3)
lambda_between <-c("age8" = -1, "age10" = 0, "age12" = 1)

contr_mx <- calc_contrast(dv = dv, 
                          between = between,
                          lambda_between = lambda_between,
                          within = within,
                          lambda_within = lambda_within,
                          id = id, 
                          data = rosenthal_tbl53)
contr_mx
#> 
#> We ran a contrast analysis for the following between contrasts: age10 = 0; age12 = 1; age8 = -1 and within contrasts: 1 = -3; 2 = -1; 3 = 1; 4 = 3. This resulted in statistics of F(1,6) = 20.211; p = 0.004123 and an effect magnitude of r_effectsize = 0.871.
```

The results look like a contrast analysis for between-subject designs.
The summary gives some more details: The effect sizes, within group
means and standard errors of the *L*-values.

``` r
summary(contr_mx)
#> $F_Table
#>              SS df     MS      F     p
#> contrast 42.667  1 42.667 20.211 0.004
#> within   12.667  6  2.111     NA    NA
#> total    56.222  8     NA     NA    NA
#> 
#> $Effects
#>              effect
#> r_effectsize  0.871
#> r_contrast    0.878
#> r_alerting    0.990
#> 
#> $Within_Groups
#>              M        SE
#> age10 4.000000 1.0000000
#> age12 7.333333 0.8819171
#> age8  2.000000 0.5773503
```

## Aggregated Data

Sometimes you would like to run a contrast analysis on aggregated data
(e.g. when no raw data is available). If you have the means, standard
deviations and sample sizes for every condition, you can do this with
cofad. For instance, if we take our first example and aggregate it, we
can still calculate the contrast analysis:

``` r
library(dplyr)
furr_agg <- furr_p4 %>% 
  group_by(major) %>% 
  summarize(mean = mean(empathy), sd = sd(empathy), n = n())
lambdas = c("psychology" = 1, "education" = -1, "business" = 0, "chemistry" = 0)
calc_contrast_aggregated(mean, sd, n, major, lambdas, furr_agg)
#> 
#> We ran a contrast analysis for the following between contrasts: business = 0; chemistry = 0; education = -1; psychology = 1. This resulted in statistics of F(1,16) = 6.154; p = 0.02461 and an effect magnitude of r_effectsize = -0.276. Attention: Contrast fits in the opposite direction!
```

And the result is indeed the same when compared to the analysis with the
raw data:

``` r
ca <- calc_contrast(dv = empathy, between = major,
                    lambda_between = c("psychology" = 1, "education" = -1,
                                       "business" = 0, "chemistry" = 0),
                    data = furr_p4)
ca
#> 
#> We ran a contrast analysis for the following between contrasts: business = 0; chemistry = 0; education = -1; psychology = 1. This resulted in statistics of F(1,16) = 6.154; p = 0.02461 and an effect magnitude of r_effectsize = -0.276. Attention: Contrast fits in the opposite direction!
```

Note that this will only work for between-subjects designs.

## Issues and Support

If you find any bugs, please use the issue tracker at:

<https://github.com/johannes-titz/cofad/issues>

If you need answers on how to use the package, drop an e-mail at
johannes at titz.science or johannes.titz at gmail.com

## Contributing

Comments and feedback of any kind are very welcome! We will thoroughly
consider every suggestion on how to improve the code, the documentation,
and the presented examples. Even minor things, such as suggestions for
better wording or improving grammar in any part of the package, are more
than welcome.

If you want to make a pull request, please check that you can still
build the package without any errors, warnings, or notes. Overall,
simply stick to the R packages book: <https://r-pkgs.org/> and follow
the code style described here: <http://r-pkgs.had.co.nz/r.html#style>

## Acknowledgments

We want to thank Thomas Schäfer and Isabell Winkler for testing cofad
and giving helpful feedback.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
line-spacing="2">

<div id="ref-furr2004" class="csl-entry">

Furr, R. M. (2004). Interpreting effect sizes in contrast analysis.
*Understanding Statistics*, *3*, 1–25.
<https://doi.org/10.1207/s15328031us0301_1>

</div>

<div id="ref-rosenthal1985" class="csl-entry">

Rosenthal, R., & Rosnow, R. L. (1985). *Contrast analysis: Focused
comparisons in the analysis of variance*. Cambridge, England: Cambridge
University Press.

</div>

<div id="ref-rosenthal2000" class="csl-entry">

Rosenthal, R., Rosnow, R. L., & Rubin, D. B. (2000). *Contrasts and
Effect Sizes in Behavioral Research: A Correlational Approach*.
Cambridge University Press.

</div>

<div id="ref-sedlmeier2018" class="csl-entry">

Sedlmeier, P., & Renkewitz, F. (2018). *Forschungsmethoden und Statistik
für Psychologen und Sozialwissenschaftler* (3rd ed.). Hallbergmoos,
Germany: Pearson Studium.

</div>

</div>
#

# cofad 0.2.1

* small improvements in documentation, references and paper for the official
publication at journal of open source software

# cofad 0.2.0

* Added a `NEWS.md` file to track changes to the package.
* Added Shiny GUI.
* Improved structure of the package.
* Fixed Bug with 0-variance conditions.
* Improved README.
* Improved examples.
* Improved documentation.
* Added and documented data sets.
* Added function for aggregated data.
---
title: "cran-comments"
output: html_document
---
## Test environments
* local Win 7, R 3.6.0
* local Scientific Linux 7.6, R 3.6.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs or NOTES. 


## Note in Win builder
checking CRAN incoming feasibility ... NOTE

## Additional info

Resubmission after Fixing.

Analysis  --> analysis
Designs   --> designs
,and      --> , and

add ISBN to reference Description text in the form
---
title: 'cofad: An R package and shiny app for contrast analysis'
tags:
  - contrast analysis
  - factorial design
  - shiny
  - R
authors:
  - name: Johannes Titz
    orcid: 0000-0002-1102-5719
    affiliation: 1
  - name: Markus Burkhardt
    orcid: 0000-0002-2498-2999
    affiliation: 1
affiliations:
 - name: Department of Psychology, TU Chemnitz, Germany
   index: 1
date: 15 September 2021
bibliography: library_paper.bib
---

# Summary

Contrast analysis offers a powerful method to test a hypothesis about means in factorial designs. Despite many advantages, the method is unknown to most researchers. For instance, even if researchers have a specific hypothesis for group means in an experimental design, they usually run an omnibus test or multiple post-hoc tests. The first approach reduces power while the second increases the type I error because of multiple testing. Contrast analysis is always the better alternative as it provides a single specific test with high power. The method is efficient, informative and shines with its simplicity [@furr2004]. To revive and popularize contrast analysis, a software is needed that is flexible but also intuitive to use.

# Statement of Need

`cofad` is an R package and shiny app for conducting **co**ntrast analysis in **fa**ctorial **d**esigns. Although it is possible to use contrasts in most statistical software tools, these are ad-hoc solutions that lack a consistent framework. In contrast, `cofad` provides a single function for between, within and mixed designs in the tradition and language of @rosenthal1985, @Rosenthal2000 and @sedlmeier2018. In addition, a graphical user interface in the form of a shiny app can either be used locally or online at https://cofad.titz.science. 

`cofad` is targeted at researchers looking for better ways to evaluate their theories. In particular contrast analysis has some tradition in psychology and the social sciences. But other fields using factorial designs will also benefit from higher test power. Practical research examples of using `cofad` can be found in a german textbook on data anlysis with R [@sedlmeier2021]. Furthermore, the package has been used in undergraduate courses at TU Chemnitz and is part of a free online course for R (https://rlernen.de).

# Alternatives

The R packages `lsmeans` [@lenth2016] and `multcomp` [@hothorn2008] provide contrast coding from a multiple comparison perspective. These tools are used to test typical (multiple) contrasts such as *all possible pairs* or *treatment versus control*. Although it is possible to use a more specific contrast it involves multiple steps and is less convenient than in `cofad`. For instance, the order of contrasts has to strictly adhere to the order of the independent variable, which can quickly lead to mistakes that are difficult to spot. Furthermore, typical effect sizes in contrast analysis such as $r_\mathrm{effectsize}$ are not reported. `cofad` is easier to use because the model and the contrast can be set in a single step. Furthermore, errors in specifying contrasts are unlikely because a named vector has to be provided with all conditions of the independent variable.

SPSS [@ibm2020] offers a function for planned contrasts in an ANOVA. But to our knowledge it cannot handle within designs or mixed designs. A minor annoyance is that adding contrasts involves many steps and the interface does not show which contrasts are mapped to which group. Again, `cofad` is more intuitive to use and prevents input errors.

# Acknowledgments

We want to thank Thomas Schäfer and Isabell Winkler for testing `cofad` and giving helpful feedback.

# References
