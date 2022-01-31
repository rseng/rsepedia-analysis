---
title: 'StatAid: An R package with a graphical user interface for data analysis'
authors:
- affiliation: 1, 2
  name: Vincent Alcazer
  orcid: 0000-0003-1843-6286
date: "20 October 2020"
bibliography: paper.bib
tags:
- R
- Data analysis
- Medicine
- Science
- Survival analysis
affiliations:
- index: 1
  name: Cancer Research Center of Lyon, INSERM U1052, Lyon, FRANCE
- index: 2
  name: Hospices Civils de Lyon, Lyon, FRANCE

---

# Summary

Data analysis is a crucial step in every research project in life science. Every clinician or researcher is one day faced with the need to perform statistical analysis. However, few free accessible solutions exist to date and most of the relevant software programs need a paid license. R is a free language which allows one to perform statistical analysis [@RCoreTeam:2017].
While the R environment is very powerful, its learning curve can be very steep in the beginning, especially for people with no previous coding skills or those with less time to learn them. A graphical user interface has already been provided as an independent package, but its features are limited for medical and applied life science studies and its usage remains difficult and unintuitive for new users [@Fox:2005]. Other free software programs exist, such as iNZight or Jamovi. However, while providing solutions with multiple features such as variable recoding, they do not guide the user through the analysis and can lack some key features such as time-dependent outcome analysis.

`StatAid` is a free open-source software provided as an R package which allows clinicians and researchers to perform statistical analysis through an intuitive graphical interface. It has been developed using the Shiny package [@Chang:2020], while Golem was used for package compilation and deployment [@Guyader:2020].

This software guides the user through all the steps of data analysis and includes multiple features such as:

- Exploratory data analysis: distribution, count, missing-values and outliers check

- Descriptive analysis, simple comparative analysis and publication ready 'table 1' output 

- Publication-ready graph customization 

- Paired data analysis (e.g. for repeated measures and matched case-control studies)

- Univariate analysis and models for continuous and categorical outcome: correlation, linear and logistic regression 

- Univariate analysis and models for time-dependent outcome: Kaplan-Meier curves and Cox regression 

- Multivariate analysis and models for continuous, categorical and time-dependent outcomes

Its user-friendly interface can guide clinicians or researchers, even those without previous software experience, through statistical analysis. In addition to a local version, a ready-to-use online version [(https://vincentalcazer.shinyapps.io/StatAid/)](https://vincentalcazer.shinyapps.io/StatAid/) is also available, providing access to `StatAid` everywhere, even on hospital/research center computers where external software installation is not allowed.
 

# Statement of need 

`StatAid` is an R package designed to fit the needs of every-day statistical analysis in science. This package provides all the tools necessary to perform data analysis in an intuitive, ready-to-use graphical interface. Users without coding skills or the availability of softwares with paid-licenses can easily perform all the steps of a good statistical analysis, from data-exploration/quality check to multivariate modeling. 

Compared with other similar free software, `StatAid` has been designed to quickly produce publication-ready graphs and tables by guiding the user through their data analysis and providing multiple graph customization options. By limiting the number of choices and integrating different checks and variable controls, `StatAid` helps the user prevent bad test use or bad graph choice. Besides, as an evolving software, `StatAid` also provides the possibility for users to request the implementation of additional features or to contribute to software development. Its open-source aspect can also be seen as a security for people working with sensitive data (e.g. data from clinical trials/patients). The online version of `StatAid` enables it to be accessable everywhere, even on computers with restrictive policies for software installation.

`StatAid` was designed to be used by clinicians, researchers, students, and any person wanting to perform statistical analysis with no prior coding skills. Primarily designed for medical/life science data analysis, `StatAid` can also easily be extended to other fields.



# References

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Welcome to StatAid

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

`StatAid` is a free open-source software provided as an R package
allowing clinicians and researchers to perform statistical analysis
through an intuitive graphical interface. It has been developed with the
R software, using the [Shiny package](https://shiny.rstudio.com/).
[Golem](https://github.com/ThinkR-open/golem) has been used for package
compilation and deployment.

The software guides the users through the steps of a good data analysis,
including multiple features such as:

<ul>

<li>

Exploratory data analysis: distribution, count, missing-values and
outliers check

</li>

<li>

Descriptive analysis, simple comparative analysis and publication ready
‘table 1’ output

</li>

<li>

Publication-ready graph customization

</li>

<li>

Paired data analysis (matched case-control studies, repeated measures)

</li>

<li>

Univariate analysis and models for continuous and categorical outcome:
Correlation, linear and logistic regression

</li>

<li>

Univariate analysis and models for time-dependent outcome: Kaplan-Meier
curves and cox regression

</li>

<li>

Multivariate analysis and models for continuous, categorical and
time-dependent outcomes

</li>

</ul>

# Getting started

## Online version

StatAid has a ready-to-use online version available at
<https://vincentalcazer.shinyapps.io/StatAid/>.

## Local version

You can install the development version from
[GitHub](https://github.com/VincentAlcazer/StatAid) else by cloning the
repository or directly by downloading the package in R:

``` r

install.packages("remotes")
remotes::install_github("VincentAlcazer/StatAid")

StatAid::run_app()
```

## Quick-start user guide

If you are not familiar with StatAid or just want to have an overview of
the different possibilities, you can check the [StatAid’s quick-start
user
guide](https://github.com/VincentAlcazer/StatAid/blob/master/STATAID_QUICK_START_USER_GUIDE.pdf)

## Citing StatAid

If you found StatAid useful and used it for your research, please cite
the [paper published in the Journal of Open Source
Software.](https://joss.theoj.org/papers/10.21105/joss.02630)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02630/status.svg)](https://doi.org/10.21105/joss.02630)

# Troubleshooting and contribution

All troubleshooting and contributions can be found on the [Github
page.](https://github.com/VincentAlcazer/StatAid/issues)

## Bug report

If you encounter any problem with the software or find a bug, please
report it on GitHub:

  - Create a [new
    issue](https://github.com/VincentAlcazer/StatAid/issues) on the
    Github page
  - Try to describe the problem/bug with reproductible steps

## Feature request

To ask for new feature implementation/current feature enhancemenet:

  - Create a [new
    issue](https://github.com/VincentAlcazer/StatAid/issues) on the
    Github page
  - Briefly describe the research question you want to answer and the
    type of data you have
  - If possible: provide pictures of the graph you would like to make or
    references from the paper you saw the analysis in.

## Contribution proposal

Contributions to new features or code enhancement are welcomed by
creating a new [pull
request.](https://github.com/VincentAlcazer/StatAid/pulls)
# StatAid 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(https://www.contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/.
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

# Welcome to StatAid

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

`StatAid` is a free open-source software provided as an R package allowing clinicians and researchers to perform statistical analysis through an intuitive graphical interface. It has been developed with the R software, using the [Shiny package](https://shiny.rstudio.com/). [Golem](https://github.com/ThinkR-open/golem) has been used for package compilation and deployment.  

 The software guides the users through the steps of a good data analysis, including multiple features such as:
<ul><li> Exploratory data analysis: distribution, count, missing-values and outliers check</li>
<li> Descriptive analysis, simple comparative analysis and publication ready 'table 1' output </li>
<li> Publication-ready graph customization</li>
<li> Paired data analysis (matched case-control studies, repeated measures) </li>
<li> Univariate analysis and models for continuous and categorical outcome: Correlation, linear and logistic regression</li>
<li> Univariate analysis and models for time-dependent outcome: Kaplan-Meier curves and cox regression </li>
<li> Multivariate analysis and models for continuous, categorical and time-dependent outcomes </li></ul>
<li> ROC Curves </li></ul>


# Getting started

## Online version 

StatAid has a ready-to-use online version available at [https://vincentalcazer.shinyapps.io/StatAid/](https://vincentalcazer.shinyapps.io/StatAid/).

## Local version

You can install the development version from [GitHub](https://github.com/VincentAlcazer/StatAid) else by cloning the repository or directly by downloading the package in R:

```{r Github install }

install.packages("remotes")
remotes::install_github("VincentAlcazer/StatAid")

StatAid::run_app()
```

## Quick-start user guide

If you are not familiar with StatAid or just want to have an overview of the different possibilities, you can check the [StatAid's quick-start user guide](https://github.com/VincentAlcazer/StatAid/blob/master/STATAID_QUICK_START_USER_GUIDE.pdf)

## Citing StatAid

If you found StatAid useful and used it for your research, please cite the [paper published in the Journal of Open Source Software.](https://joss.theoj.org/papers/10.21105/joss.02630) [![DOI](https://joss.theoj.org/papers/10.21105/joss.02630/status.svg)](https://doi.org/10.21105/joss.02630)


# Troubleshooting and contribution

All troubleshooting and contributions can be found on the [Github page.](https://github.com/VincentAlcazer/StatAid/issues)

## Bug report

If you encounter any problem with the software or find a bug, please report it on GitHub:

- Create a [new issue](https://github.com/VincentAlcazer/StatAid/issues) on the Github page
- Try to describe the problem/bug with reproductible steps

## Feature request

To ask for new feature implementation/current feature enhancemenet:

- Create a [new issue](https://github.com/VincentAlcazer/StatAid/issues) on the Github page
- Briefly describe the research question you want to answer and the type of data you have
- If possible: provide pictures of the graph you would like to make or references from the paper you saw the analysis in.

## Contribution proposal

Contributions to new features or code enhancement are welcomed by creating a new [pull request.](https://github.com/VincentAlcazer/StatAid/pulls)

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{regression_table}
\alias{regression_table}
\title{All-in-one publication ready regression table for glm models
This function output a publication-ready regression table from a dataset. Each column of the dataset is independantly tested as a x variable in a univariate model.}
\usage{
regression_table(data, y_var, family = "gaussian")
}
\arguments{
\item{data}{A dataframe with row corresponding to samples/patients and columns to variables.}

\item{y_var}{A character string corresponding to the y variable.}

\item{family}{Generalized linear model family (gaussian or binomial)}
}
\description{
All-in-one publication ready regression table for glm models
This function output a publication-ready regression table from a dataset. Each column of the dataset is independantly tested as a x variable in a univariate model.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_app.R
\name{run_app}
\alias{run_app}
\title{Run the Shiny Application}
\usage{
run_app(...)
}
\arguments{
\item{...}{A series of options to be used inside the app.}
}
\description{
Run the Shiny Application
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{aov_pipe_pval}
\alias{aov_pipe_pval}
\title{Anova p value pipe
Allows the extraction of an ANOVA p value in a pipe}
\usage{
aov_pipe_pval(data, variable, group)
}
\arguments{
\item{data}{A dataframe with row corresponding to samples/patients and columns to variables.}

\item{variable}{A character string corresponding to the studied variable.}

\item{group}{A character string corresponding to the comparative groups.}
}
\description{
Anova p value pipe
Allows the extraction of an ANOVA p value in a pipe
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{regression_dataframes}
\alias{regression_dataframes}
\title{Regression base and diagnosis dataframe
Provides augmented data for regression diagnosis}
\usage{
regression_dataframes(y_var, x_var, data, model = "lm", cor_type = "pearson")
}
\arguments{
\item{y_var}{A character string corresponding to the y variable.}

\item{x_var}{A character string corresponding to the x variable(s).}

\item{data}{A dataframe with row corresponding to samples/patients and columns to variables.}

\item{model}{Type of generalised linear model (lm or gam).}

\item{cor_type}{Type of correlation to show on the graph (pearson or spearman).}
}
\description{
Regression base and diagnosis dataframe
Provides augmented data for regression diagnosis
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{regression_table_cox}
\alias{regression_table_cox}
\title{All-in-one publication ready regression table for cox models
This function output a publication-ready regression table from a dataset. Each column of the dataset is independantly tested as a x variable in a univariate model.}
\usage{
regression_table_cox(data, y_var, time_var)
}
\arguments{
\item{data}{A dataframe with row corresponding to samples/patients and columns to variables.}

\item{y_var}{A character string corresponding to the y variable = time-dependant outcome (0-1 or dead-alive for example).}

\item{time_var}{A numeric variable corresponding to the time variable.}
}
\description{
All-in-one publication ready regression table for cox models
This function output a publication-ready regression table from a dataset. Each column of the dataset is independantly tested as a x variable in a univariate model.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{regression_table_multi}
\alias{regression_table_multi}
\title{All-in-one publication ready regression table for glm models
This function output a publication-ready regression table from a dataset. Each column of the dataset is included as a x variable in a multivariate model.}
\usage{
regression_table_multi(data, y_var, family = "gaussian")
}
\arguments{
\item{data}{A dataframe with row corresponding to samples/patients and columns to variables.}

\item{y_var}{A character string corresponding to the y variable.}

\item{family}{Generalized linear model family (gaussian or binomial)}
}
\description{
All-in-one publication ready regression table for glm models
This function output a publication-ready regression table from a dataset. Each column of the dataset is included as a x variable in a multivariate model.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{autoplot}
\alias{autoplot}
\title{Automatised plot
This is a wrapper function to set a plot with different parameters from ggplot.}
\usage{
autoplot(
  data,
  variable,
  group,
  group_filter_vector = NULL,
  na_exclude_group = T,
  plot_type = "Boxplot",
  add_points = T,
  error_bar = "IC95",
  stat = "param"
)
}
\arguments{
\item{data}{A dataframe with row corresponding to samples/patients and columns to variables.}

\item{variable}{A character string corresponding to the studied variable.}

\item{group}{A character string corresponding to the comparative groups.}

\item{group_filter_vector}{A character string/vector defining groups sublevels to be included in the analysis.}

\item{na_exclude_group}{Should missing values be excluded from groups?}

\item{plot_type}{Set the plot type (Boxplot, Barchart_mean or  Barchart_count).}

\item{add_points}{Add data points as a scatter plot to the graph.}

\item{error_bar}{Type of error bar to be shown on the barcharts (IC95 or hide).}

\item{stat}{Type of statistics to show (param, non_param or no).}
}
\description{
Automatised plot
This is a wrapper function to set a plot with different parameters from ggplot.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{regression_table_multi_cox}
\alias{regression_table_multi_cox}
\title{All-in-one publication ready regression table for cox models
This function output a publication-ready regression table from a dataset. Each column of the dataset is included as a x variable in a multivariate model.}
\usage{
regression_table_multi_cox(data, y_var, time_var)
}
\arguments{
\item{data}{A dataframe with row corresponding to samples/patients and columns to variables.}

\item{y_var}{A character string corresponding to the y variable.}

\item{time_var}{A numeric variable corresponding to the time variable.}
}
\description{
All-in-one publication ready regression table for cox models
This function output a publication-ready regression table from a dataset. Each column of the dataset is included as a x variable in a multivariate model.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{autoplot_paired}
\alias{autoplot_paired}
\title{Automatised plot for paired data
This is a wrapper function to set a plot with different parameters from ggplot for paired data.}
\usage{
autoplot_paired(
  data,
  timepoints,
  group,
  plot_type = "Mean_lines",
  add_points = F,
  add_lines = T,
  add_individual_lines = F,
  error_bar = "IC95",
  stat = "param",
  alpha_line = 0.5
)
}
\arguments{
\item{data}{A dataframe with row corresponding to samples/patients and columns to variables.}

\item{timepoints}{A character vector corresponding to the studied paired variables.}

\item{group}{A character string corresponding to the comparative groups.}

\item{plot_type}{Set the plot type (Boxplot, Barchart_mean or  Barchart_count).}

\item{add_points}{Add data points as a scatter plot to the graph.}

\item{add_lines}{Add line to paired data}

\item{add_individual_lines}{Add each individual (samples) lines to the plot.}

\item{error_bar}{Type of error bar to be shown on the barcharts (IC95 or se or hide).}

\item{stat}{Type of statistics to show (param, non_param or no).}

\item{alpha_line}{Transparency of the paired lines}
}
\description{
Automatised plot for paired data
This is a wrapper function to set a plot with different parameters from ggplot for paired data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{descriptive_table}
\alias{descriptive_table}
\title{All-in-one publication ready descriptive table
This function output a publication-ready descriptive table from a dataset. Each columns of the dataset is compared between the defined groups.}
\usage{
descriptive_table(
  data,
  group,
  na.include = F,
  percent_type = 1,
  padj_method = "none",
  show_methods = F,
  exclude_vector = c("Patient_id", "patient_id", "Sample_ID", "Whole_cohort")
)
}
\arguments{
\item{data}{A dataframe with row corresponding to samples/patients and columns to variables.}

\item{group}{A character string corresponding to the comparative groups.}

\item{na.include}{Should missing values be included in groups?}

\item{percent_type}{For categorical variables, should percentage be calculated on the row (1) or columns (2)?}

\item{padj_method}{Select the p.value adjustment method (none, fdr, holm or bonferonni).}

\item{show_methods}{Should statistical tests names be precised in a supplementary column?}

\item{exclude_vector}{columns to exclude from the analysis.}
}
\description{
All-in-one publication ready descriptive table
This function output a publication-ready descriptive table from a dataset. Each columns of the dataset is compared between the defined groups.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{ROC_curves}
\alias{ROC_curves}
\title{Dedicated function to plot ROC Curves
This function output individual ROC Curves for each x-variable provided.}
\usage{
ROC_curves(data, y_var, x_var)
}
\arguments{
\item{data}{A dataframe with row corresponding to samples/patients and columns to variables.}

\item{y_var}{A character string corresponding to the y variable. !!Must be a 2-level factor (TRUE/FALSE, 1/0)}

\item{time_var}{A vector of numeric variable to test}
}
\description{
Dedicated function to plot ROC Curves
This function output individual ROC Curves for each x-variable provided.
}
