---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# c212: Methods for Detecting Safety Signals in Clinical Trials Using Body-Systems (System Organ Classes) 

<!-- badges: start -->
<!-- badges: end -->

The goal of c212 is to provide a self-contained set of methods, which use groupings of adverse events, to aid clinical trial safety investigators, statisticians and researchers, in the early detection of adverse events.

## Installation

You can install the released version of c212 from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("c212")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("rcarragh/c212")
```

## Example

This is a basic example which shows how to apply the Double False Discovery Rate to a set of multiple hypotheses:


```r
library(c212)
data(c212.FDR.data)
c212.err.cntrl(c212.FDR.data, method="DFDR", alpha = 0.05)
#>           B j         AE        p
#> 1 Bdy-sys_3 1   Adv-Ev_6 0.000000
#> 2 Bdy-sys_3 2   Adv-Ev_7 0.000011
#> 3 Bdy-sys_3 3   Adv-Ev_8 0.000021
#> 4 Bdy-sys_3 5 Adv-Ev_910 0.000039
#> 5 Bdy-sys_3 6 Adv-Ev_911 0.000079
#> 6 Bdy-sys_3 7 Adv-Ev_912 0.003554
#> 7 Bdy-sys_3 4   Adv-Ev_9 0.010411
```

This is an example of how to apply the Berry and Berry model:


```r
library(c212)
data(c212.trial.data)
mod.BB <- c212.BB(c212.trial.data, burnin = 100, iter = 200)
#> Global Simulation Parameters:
#> 	Simulation Type: 1
#> 	w_alpha (width): 1.000000
#> 	m alpha (control): 6.000000
#> 	w_beta (width): 1.000000
#> 	m beta (control): 6.000000
#> 	w_gamma (width): 1.000000
#> 	m gamma (control): 6.000000
#> 	sigma_MH_alpha: 3.000000
#> 	sigma_MH_beta: 3.000000
#> 	sigma_MH_gamma: 0.200000
#> 	sigma_MH_theta: 0.200000
#> 	default weight: 0.500000
#> MCMC chain fitting complete.
#> [1] "MCMC fitting complete."
```

# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

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
* Focusing on what is best not just for us as individuals, but for the
  overall community

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

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
this [email](rcarragh@gmail.com).
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

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

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
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.

# Contributing to `c212`

The author is interested in extending the software to include new methods, particularly the area of safety analysis, and would welcome collaborations in this area. 

## Issues and Extensions

For any detected bugs, issues or possible extensions please feel free to create an issue in the c212 github repository or contact the author directly.

When reporting issues, if possible, please include a minimal example with a small data set which reproduces the issue.

---
title: 'c212: An R Package for the Detection of Safety Signals in Clinical Trials Using Body-Systems (System Organ Classes)'
tags:
  - Adverse event
  - Safety
  - Pharmacovigilance
  - Clinical trials
  - False Discovery Rate
  - Body-systems
  - System organ classes
  - Bayesian hierarchy
  - Multiple comparisons.
authors:
  - name: Raymond Carragher
    orcid: 0000-0002-0120-625X
    affiliation: "1, 2, 4"
  - name: Chris Robertson
    affiliation: "2, 3"
affiliations:
  - name: Strathclyde Institute of Pharmacy and Biomedical Sciences, University of Strathclyde, Glasgow, UK
    index: 1
  - name: Department of Mathematics and Statistics, University of Strathclyde, Glasgow, UK
    index: 2
  - name: Health Protection Scotland, NHS National Services Scotland, Glasgow, UK
    index: 3
  - name: Health Data Research (UK), University of Strathclyde, Glasgow, UK
    index: 4
bibliography: paper.bib
date: \today
---

# Summary
Safety in clinical trials may be characterised by the incidence or occurrence of adverse events. The statistical analysis of this data is 
complicated by the large number of adverse events recorded, with low event rates, small effect sizes and low power all contributing to the 
difficulty in determining a robust safety profile for a treatment during the trial process. 

In addition to end of trial analyses, a number of interim analyses may take place at different time points during the trial lifecycle. These 
offer the additional statistical challenge of testing accumulating data, with possibly differing recruitment rates on trial arms contributing 
to a lack of balance in the data.

Adverse events are typically defined by medical dictionaries, which provide a common reference terminology for use in and between clinical trials.
There are a number of medical dictionaries in current use, all of which provide similar services. One such dictionary is [MedDRA](http://www.meddra.org/) 
(Medical Dictionary for Regulatory Activities), which was developed by the [ICH](http://www.ich.org/) (International Council for Harmonisation of 
Technical Requirements for Pharmaceuticals for Human Use) and is widely used by regulatory bodies, clinical research 
organisations (CROs), and pharmaceutical companies. WHO-ART (World Health Organisation Adverse Reaction Terminology) is a similar dictionary 
maintained by the [Uppsala Monitoring Centre](http://www.umc-products.com) for the [World Health Organisation Collaborating Centre for 
International Drug Monitoring](http://bioportal.bioontology.org/ontologies/WHO-ART). MedDRA and WHO-ART have a similar hierarchical structure 
consisting of System Organ Classes (SOC) and various grouping and descriptor terms. 

The MedDRA hierarchical structure consists of five levels: System Organ Class (SOC), High Level Group Terms (HLGT), High Level Terms (HLT), 
Preferred Terms (PT), and Lower Level Terms (LLT). The PT is a single medical description of a symptom or observation while the LLT is how a 
patient or data recorder would describe a symptom or observation. Each LLT belongs to one PT and, in general, data will be recorded at the LLT 
level but reported at the PT level (the adverse event). As of 2020 there are 27 SOCs and over 80,000 LLTs.

The grouping of adverse events by SOC (or body-system) provides for possible relationships between the adverse events within a SOC.
One consequence of this is the possibility that, for treatments which may affect a particular SOC, there may be raised rates for a number of adverse events
within that SOC.
A number of methods have recently been proposed to address the statistical issues in adverse event analysis by using these groupings of adverse
events by body-system or SOC, taking into account the additional information provided by these relationships to increase the power of detecting real 
adverse event effects.
These methods, which include both error controlling procedures for multiple hypothesis 
testing [@BH1995; @HZZ2010; @MA2012; @Y2008; @M2006],
and Bayesian modelling approaches [@BB2004; @XMC2011; @phdthesis], are implemented in the `R` package `c212` (Table \ref{table:1}). 

| Method                                                | Description                           |
| ----------------------------------------------------- | ------------------------------------- |
| Benjamini-Hochberg procedure [@BH1995]                 | Control of False Discovery Rate       |
| Group Benjamini-Hochberg procedure [@HZZ2010]          | Control of False Discovery Rate       |
| Double false discovery rate [@MA2012]                  | Control of False Discovery Rate       |
| Subset BH-procedure [@Y2008]                           | Control of False Discovery Rate       |
| Bonferroni correction [@M2006]                         | Control of Familywise Error Rate      |
| Berry and Berry model [@BB2004]                        | Bayesian model for end of trial data  |
| Berry and Berry model without point-mass [@XMC2011]    | Bayesian model for end of trial data  |
| Interim analysis model [@phdthesis]                    | Bayesian model for interim trial data |
| Interim analysis model without point-mass [@phdthesis] | Bayesian model for interim trial data |
: Methods in the c212 package. \label{table:1}

# Statement of Need 
The detection of safety issues in the post-marketing phase of a treatment's life cycle, as opposed to the trial phase, can have a serious effect 
on the health of patients and also a financial impact both for the companies developing the treatments, and the regulatory bodies responsible for 
overseeing them. 

The `R` package `c212` provides a self-contained set of methods for clinical trial safety investigators, statisticians and researchers, to aid in
the early detection of adverse events. It is designed to be easy to use with a simple data input format and interface. 

The primary use case for the software is in the statistical analysis of adverse event incidence and occurrence data during clinical trials. 
A second goal of the package is to provide reference implementations of the methods in Table \ref{table:1} for use by researchers, both in the 
area of safety in clinical trials, as well those developing or testing methods for 
handling error rates when testing multiple hypotheses.
Beyond safety in clinical trials, the package will be useful to any project which deals with multiple hypothesis testing, or projects where two 
groups of comparative data may be modelled by hierarchical Bayesian binomial or Poisson models, with recent extensions of the Bayesian models to 
observational data being developed [@RC2020].


The `c212` package is currently being used both for clinical trial safety analysis [@FKL2016; @Munsaka2018; @WWM2018] and as a research tool in the
investigation and development of new safety methods [@Tan2020; @Tan2019].
It has also been [implemented](https://visual-analytics.shinyapps.io/index/) as part of an [ASA Biopharm Safety Working Group](https://community.amstat.org/biop/workinggroups/safety/safety-home) workstream.



# Overview

The Bayesian models, under assumptions of conditional independence, are fitted using a Gibbs sampling Markov Chain Monte-Carlo (MCMC) 
method [@RCR1999].
The posterior distributions of the model parameters are used to assess which adverse events may have increased rates on the treatment arm. 
In the case of the Berry and Berry model, which is binomial, the `theta` model parameter, representing the increase in the log-odds of an event 
occurring on the treatment arm, is used for this purpose [@BB2004].
For the interim analysis models, which are Poisson based, the increase in the log rate of an event on the treatment arm is used for
adverse event assessment. As in the Berry and Berry model this is represented by the parameter `theta` [@phdthesis].
Functions for generating summary statistics and highest probability intervals are provided using the services of the `coda` package [@coda2006].
The main convergence diagnostics available directly within the package are the Gelman-Rubin and Geweke statistics [@GCSR2004], again from the 
`coda` package. Access to the raw samples is available for further processing should that be required.
The error controlling procedures included in the package follow exactly the method definitions in the papers which introduced 
them [@BH1995; @HZZ2010; @MA2012; @Y2008; @M2006].
The following sections contain examples which cover the main uses of the software. The data sets and functions used are fully documented in the 
package.

## Multiple Hypothesis Testing

The data set `c212.FDR.data` contains the results from a set of multiple hypothesis tests for adverse events grouped by 
body-system. The Group Benjamini-Hochberg procedure is applied to the data with control of the False Discovery Rate set at the 5% level.
The `c212.err.cntrl` function, which provides the interface to the error controlling methods, returns the set of p-values/hypotheses deemed 
significant.

```r
library(c212)
data(c212.FDR.data, package="c212")
head(c212.FDR.data, 2)
B j AE p
1 Bdy-sys_5 7 Adv-Ev_24 1.000000
2 Bdy-sys_6 9 Adv-Ev_34 0.358949

p.sig <- c212.err.cntrl(c212.FDR.data, method = "GBH", alpha = 0.05)
print(p.sig)
          B j         AE        p p_weighted
1 Bdy-sys_3 6 Adv-Ev_911 0.000079   0.000000
2 Bdy-sys_3 3   Adv-Ev_8 0.000021   0.000000
3 Bdy-sys_3 1   Adv-Ev_6 0.000000   0.000000
4 Bdy-sys_3 5 Adv-Ev_910 0.000039   0.000000
5 Bdy-sys_3 2   Adv-Ev_7 0.000011   0.000000
6 Bdy-sys_3 4   Adv-Ev_9 0.010411   0.000000
7 Bdy-sys_3 7 Adv-Ev_912 0.003554   0.000000
8 Bdy-sys_2 1   Adv-Ev_2 0.005333   0.005333
9 Bdy-sys_2 2   Adv-Ev_3 0.016013   0.016013
```


## Berry and Berry End of Trial Analysis

The data set `c212.trial.data` contains sample end of trial adverse event incidence counts.
The data is modelled using the Berry and Berry model as follows:

```r
library(c212)
data(c212.trial.data, package="c212")
head(c212.trial.data, 2)
B j AE Group Count Total
1 Bdy-sys_2 1 Adv-Ev_2 1 20 450
2 Bdy-sys_2 4 Adv-Ev_5 2 21 450

mod.BB <- c212.BB(c212.trial.data)
```

`mod.BB` contains the raw samples generated from the model fitting procedure. To perform a convergence check:

```r
conv.BB = c212.convergence.diag(mod.BB)
c212.print.convergence.summary(conv.BB)
```

In order to assess which adverse events may be associated with treatment the function `c212.ptheta` is used. This calculates the posterior 
probability of an increase in log-odds of an event occurring on the treatment arm. 
A threshold may be used to view the adverse events which exceed some defined level, for example: 0.90:

```r
theta.pos_BB <- c212.ptheta(mode.BB)

theta_pos.BB[theta_pos.BB$ptheta > 0.9,]
B AE ptheta
2 Bdy-sys_2 Adv-Ev_2 0.9410417
6 Bdy-sys_3 Adv-Ev_6 1.0000000
7 Bdy-sys_3 Adv-Ev_7 1.0000000
8 Bdy-sys_3 Adv-Ev_8 1.0000000
9 Bdy-sys_3 Adv-Ev_9 0.9873500
10 Bdy-sys_3 Adv-Ev_910 0.9999333
11 Bdy-sys_3 Adv-Ev_911 0.9999917
12 Bdy-sys_3 Adv-Ev_912 0.9964917
```

The high posterior probabilities may indicate a possible association with treatment for these adverse events.

## Interim Analysis
Apart from the function used to fit the model, the procedure for fitting and accessing interim analysis data is exactly the same as for the 
Berry and Berry model. 

```r
library(c212)
data(c212.trial.interval.data1)
head(c212.trial.interval.data1, 2)
Interval I_index         B       AE Group Count Exposure
0.0-180.0       0 Bdy-sys_1 Adv-Ev_1     1    87 160133.7
0.0-180.0       0 Bdy-sys_1 Adv-Ev_1     2   103 163224.6
mod.BB.interim <- c212.BB.interim(c212.trial.interval.data1)
```

# Software Details and Availability

The `c212` package was initially released to CRAN in 2017 and has been through a number of release cycles.
Before each release a full set of unit and functional tests are performed on the package development system, including memory checks with 
valgrind [@valgrind] and Google address sanitizer [@asan].
The package documentation also contains tests and examples based on data included in the package.

The `c212` package is most easily downloaded and installed directly from CRAN [@cran] or, alternatively, 
from the corresponding GitHub repository [@github].

The authors are interested in extending the software to include new methods, particularly in the area of safety analysis, and would welcome 
collaborations in this area. Any support issues or questions can be addressed directly to the corresponding author, through the associated CRAN 
maintainer email address, or through the Github repository.

## Performance
The Bayesian models are expensive to fit in terms of both computation and memory. The main issue is the number of parameters in the model. For the
Berry and Berry model with $N$ adverse events, and $B$ SOCs there are $(2 \times N + 5 \times B + 6)$ parameters in total. With $C$ parallel MCMC 
chains and $I$ iterations of the MCMC sampler, this will require space for $C \times I \times (2 \times N + 5 \times B + 6)$ double precision 
numbers to store the samples. For larger datasets, the calculation of the convergence diagnostics and summary statistics may be time consuming 
as the samples are passed to the `coda` package. 

**Example:** A trial with 23 SOCs and 497 AEs, 5 parallel chains and 40,000 iterations after burn-in, will require space for 
133,800,000 doubles, equating to approximately 1GB of memory storage for the samples alone (assuming 8 bytes per double).

It is recommended that at least twice the storage required for the samples be available for fitting any of the Bayesian models.

# Acknowledgements

This work was funded by the Engineering and Physical Sciences Research Council (EPSRC) UK (award reference 1521741) and Frontier 
Science (Scotland) Ltd.

# References
