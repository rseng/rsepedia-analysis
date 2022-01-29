[![linux](https://github.com/JGCRI/ambrosia/actions/workflows/build_linux.yml/badge.svg)](https://github.com/JGCRI/ambrosia/actions/workflows/build_linux.yml) [![osx](https://github.com/JGCRI/ambrosia/actions/workflows/build_osx.yml/badge.svg)](https://github.com/JGCRI/ambrosia/actions/workflows/build_osx.yml) [![windows](https://github.com/JGCRI/ambrosia/actions/workflows/build_windows.yml/badge.svg)](https://github.com/JGCRI/ambrosia/actions/workflows/build_windows.yml) [![codecov](https://codecov.io/gh/JGCRI/ambrosia/branch/master/graph/badge.svg)](https://codecov.io/gh/JGCRI/ambrosia)
[![DOI](https://zenodo.org/badge/69679416.svg)](https://zenodo.org/badge/latestdoi/69679416)
[![status](https://joss.theoj.org/papers/b1af4bb026f674b31d1175075607deef/status.svg)](https://joss.theoj.org/papers/b1af4bb026f674b31d1175075607deef)


# `ambrosia`: An R package for calculating and analyzing food demand that is responsive to changing incomes and prices

## Summary
The `ambrosia` R package was developed to calculate food demand for staples and non-staple commodities that is responsive to changing levels of incomes and prices. `ambrosia` implements the framework to quantify food demand as established by Edmonds et al. (2017) and allows the user to explore and estimate different variables related to the food demand system. Currently `ambrosia` provides three main functions:
1. calculation of food demand for any given set of parameters including income levels and prices,
2. estimation of calibration parameters within a given a dataset.  Note:  `ambrosia` is used to calculate parameters for the food demand model implemented in the [Global Change Analysis Model](http://www.globalchange.umd.edu/gcam/).
3. exploration and preparation of raw data before starting a parameter estimation.

# Statement of need and audience
An important motivation to develop `ambrosia` is functionalizing and separating out the different components of the sophisticated food demand framework into usable R functions that can be easily parameterized and customized by the user.Thus, the tool not only enables easy use and future development, but also enables easy modularization of the code within other systems. 

`ambrosia` has been developed to help researchers explore questions related to trends in food demand empirically. Since the equations of the model are grounded in peer reviewed research while the code itself is written in R (which increases usability), the tool is useful to researchers interested in,

1)	analyzing and exploring trends in food demand with a computational model that is responsive to changes in incomes and prices that can easily be implemented on any time series (dataset);
2)  re-estimating calibration parameters of the food demand model using custom data, thus effectively allowing the user to calibrate the model to custom data;
3)	incorporating a detailed food demand model in their own earth system and economic models.


## Getting Started with `ambrosia`

`ambrosia` can be directly installed from its GitHub repository using the R `devtools` package. From an R prompt, run the command,

```r
devtools::install_github('JGCRI/ambrosia', build_vignettes = TRUE)
```

## User tutorial and examples

A list of examples along with a user tutorial describing the different features in `ambrosia` are available in the [`ambrosia_vignette`](https://jgcri.github.io/ambrosia/articles/ambrosia_vignette.html). To load the vignette with examples for the major functions within ambrosia, run the following command

```r
vignette("ambrosia_vignette")
```
The example below shows how a user can get an estimate of demand using some sample parameters (See the table below for description of parameters).

```r
#Get a sample data set
Test_Data <- data.frame(Y=seq(0.1,30, by=0.1))

#Add sample values of Ps and Pn
Test_Data %>% mutate(Ps=0.1,Pn=0.2) -> Test_Data

#Now calculate food demand
Food_Demand <- food.dmnd(Test_Data$Ps, Test_Data$Pn, Test_Data$Y)

```
The code from the example can be used to visualize the food demand for staples and non-staples as follows,

![A simple plot of food demand for staples and non-staples for changing incomes and constant prices.](vignettes/example_3.png)

#### Description of calibration parameters

The 11 calibration parameters are described in table below with values from the latest version of ambrosia. The table also contains an acceptable range for each of the parameters. The original parameters were calculated using a Markov Chain Monte Carlo (MCMC) approach. The range is calculated as the 95% Joint Confidence Interval of the range of the parameters apperaing in all Monte Carlo samples with likelihood values above the 5th percentile.Parameters in the table below are calibrated on the basis of national level data on food consumption and food consumption prices for staples and non-staple products.

| Parameter name | Description                                                                                                                                                                                                                                                           | Units        | Value | Range of parameter values |
|----------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------|-------|---------------------------|
| `A_s`            | A scale term used to derive expenditure share for staple demand                                                                                                                                                                                                       | Unitless     | 1.28  | 1.25 -1.40                |
| `A_n`            | A scale term used to derive expenditure share for non-staple demand.                                                                                                                                                                                                  | Unitless     | 1.14  | 0.9 - 1.16               |
| `xi_ss`          | Price elasticity of staple goods. Unit change in per capita demand for staples as a result of unit increase in price (in $ per person per day).                                                                                                                       | Elasticity   | -0.19 | -0.27 - -0.07         |
| `xi_cross`          | Cross price elasticity between staples and non-staples which in combination with the other price elasticities is used to derive substitution elasticity.                                                                                                              | Elasticity   | 0.21  | 0.09 - 0.27          |
| `xi_nn`          | Price elasticity of non-staple goods. Unit change in per capita demand for non-staples as a result of unit increase in price (in $ per person per day).                                                                                                                       | Elasticity   | -0.3 | -0.46 - -0.10
| `nu1_n`          | Income elasticity for non-staple goods. Unit change in per capita demand for non-staples for unit change in income (in thousand USD).                                                                                                                                 | Elasticity   | 0.5   | 0.46 - 0.61               |
| `lambda_s`       | Income elasticity for staple goods. Unit change in per capita demand for staples for unit change in income (in thousand USD)                                                                                                                                          | Elasticity   | 0.1   | 0.075 - 0.16                |
| `k_s`            | Log of Income level at which staple demand is anticipated to be at its highest                                                                                                                                                                                   | Thousand USD | 16    | 10 -17                    |
| `psscl`          | Additional scaling term used to derive the expenditure shares for staples. This is applied to price of staples (`Ps`/`Pm` * `psscl`), where Ps is the price of staples and Pm is the price of materials and to the expenditure shares of staples (`alpha_s`).                 | Unitless     | 100   | 80 - 120                  |
| `pnscl`          | Additional scaling term used to derive the expenditure shares for non-staples. This is applied to price of non-staples (`Pn`/`Pm` * `pnscl`), where Pn is the price of non-staples and Pm is the price of materials and to the expenditure shares of non-staples (`alpha_n`). | Unitless     | 20    | 18 - 25                   |
| `Pm`            | Price of materials (everything else in the economy other than food products)                                                                                                                                                                                   | $ per day | 5    | 2 -6                    |

#### Simple example of equations to derive demand using calibration parameters described above.

Below is a simple example of how quantities of demand for staples (`Q_s`) and non-staples (`Q_n`) in thousand calories are calculated using the above 11 calibration parameters above for an income level of `Y` in thousand USD for a staple price of `Ps` in $ per calorie per day and non-staple price of `Pn` in $ per calorie per day.  

```r
# 1) Staple demand

Q_s <- A_s * Ps ^ xi_ss * Pn ^ xi_cross * Income_Term_staples   

where,

Income_term_staples <- (k_s * Y) ^ (lambda_s / Y)

Pn <- Pn/Pm

Ps <- Ps/Pm

# 2) Non-staple demand

Q_n <- A_n * Ps ^ xi_cross * Pn ^ xi_nn * Income_Term_Non_staples

where,

Income_Term_Non_staples <- Y ^ (2 * nu1_n)

Pn <- Pn/Pm

Ps <- Ps/Pm

 
```


## Contributing to `ambrosia`
We welcome contributions to `ambrosia` from the development community. Please contact us at the email IDs below if you want to collaborate! The `ambrosia` GitHub repository is accessible here: [GitHub Repository](https://github.com/JGCRI/ambrosia). In order to report issues with `ambrosia`, please open an issue in the above mentioned Github Repository.

For more information about contributing, please contact Kanishka Narayan at kanishka.narayan@pnnl.gov or Chris Vernon at chris.vernon@pnnl.gov


# Availability

## Operating system
Mac OS X; Linux; Windows 10

## Programming language
R (>= 3.5.0)

## Dependencies
dplyr (>= 0.7)

nleqslv (>= 3.2)

reshape2 (>= 1.4.3)

ggplot2 (>= 2.2.1)

cluster (>= 2.0)

tidyr  (>= 0.7.1)

## Code repository

Name- GitHub; `JGCRI/ambrosia`

Identifier- https://github.com/JGCRI/ambrosia/tree/v1.3.5

License- BSD 2-Clause
---
title: 'ambrosia: An R package for calculating and analyzing food demand that is responsive to changing incomes and prices'
authors:
- affiliation: 1
  name: Kanishka Narayan
  orcid: 0000-0001-8483-6216
- affiliation: 1
  name: Chris R. Vernon
  orcid: 0000-0002-3406-6214
- affiliation: 1
  name: Stephanie Waldhoff
  orcid: 0000-0002-8073-0868
- affiliation: 1
  name: James A. Edmonds
  orcid: 0000-0002-3210-9209
- affiliation: 2
  name: Ryna Cui
  orchid: 0000-0002-8531-4126
date: "14 October 2020"
output:
  word_document: default
  pdf_document:
    fig_caption: yes
  html_document:
    df_print: paged
bibliography: paper.bib
tags:
- R
- GCAM
- food
affiliations:
- index: 1
  name: Joint Global Change Research Institute, Pacific Northwest National Laboratory, College Park, MD, USA
- index: 2
  name: University of Maryland
---

# Summary
The `ambrosia` R package was developed to calculate food demand for staples and non-staple commodities that is responsive to changing levels of incomes and prices. `ambrosia` implements the framework to quantify food demand as established by @edmonds2017global and allows the user to explore and estimate different variables related to the food demand system. Currently `ambrosia` provides three main functions:

(1)	calculation of food demand for any given set of income levels and prices;
(2)	estimation of calibration parameters within a given a dataset.  Note:  `ambrosia` is used to calculate the calibration parameters for the food demand model implemented in the Global Change Analysis Model [GCAM; @calvin2019gcam];
(3)	exploration and preparation of raw data before starting a calibration parameter estimation.

# Statement of need
An important motivation to develop `ambrosia` is functionalizing and separating out the different components of the sophisticated food demand framework from @edmonds2017global (summarized below) into usable R functions that can be easily parameterized and customized by the user. Thus, `ambrosia` has been developed to help researchers explore questions related to trends in food demand empirically. Since the equations of the model are grounded in peer reviewed research while the code itself is written in R (which increases usability), the tool is useful to researchers interested in,

1)	analyzing and exploring trends in food demand with a computational model that is responsive to changes in incomes and prices that can easily be implemented on any time series (dataset);
2)  re-estimating calibration parameters of the food demand model using custom data, thus effectively allowing the user to calibrate the model to custom data;
3)	incorporating a detailed food demand model in their own earth system and economic models.

`ambrosia` is part of an ecosystem of tools within the Global Change Intersectoral Modeling System (GCIMS) that help users computationally explore science and policy questions related to different dimensions of human-Earth systems [@pnnl_2020]. The parameters calculated from `ambrosia` are utilized directly in GCAM [@calvin2019gcam] to represent forecasts of food demand. `ambrosia` ensures that the parameters that are used within GCAM are scientifically and empirically sound and also ensures reproducibility of the parameters for validation to comply with the commitment of GCIMS to FAIR guiding principles for scientific data management [@wilkinson2016fair]. The code is structured to ensure that the parameters can be updated and tested effectively with changes to the underlying data.    

Thus, the tool not only enables easy use and future development, but also enables easy modularization of the code within other systems. The sections below contain a detailed discussion of the different functions and customization options available within the tool.

# Summary of the Edmonds et al. framework
The @edmonds2017global model represents a food demand model for staples and non-staple commodities at different levels of prices and incomes. Demand for staples is described as increasing when income is lower, eventually peaks at under 1000$ per person per capita,  and then begins to decline as higher income ranges are approached. Demand for non-staples increases with income over all income ranges; however, total (staple + non-staple) demand saturates at high income level.

The @edmonds2017global approach uses 11 calibration parameters where the parameters are fit using pooled cross-sectional-timeseries observations and a Bayesian Markov Chain Monte Carlo method [MCMC; @hastings1970]. The framework represents demand for three categories of goods: staples (s), non-staples (n) and materials (m) where materials represent everything in the economy other than staples and non-staple food commodities. The demand for these three categories changes with changes in income (Y) and prices (P), with the response to price changes varying with income. Expenditures on these three goods are assumed to exhaust income.  

Demand for these three categories can be represented mathematically as,

(1) Staple demand: $q_{s} = A_{s}(x^{h_{s}(x)})(w_{s}^{e_{ss}(x)})(w_{n}^{e_{sn}(x)})$

(2) Non-staple food demand: $q_{n} = A_{n}(x^{h_{n}(x)})(w_{s}^{e_{ns}(x)})(w_{n}^{e_{nn}(x)})$

(3) Materials demand : $q_{m} = x - w_{s}q_{s} - w_{n}q_{n}$

where $w_{i}$ is $P_{i}/P_{m}$, $x$ is $Y/P_{m}$ and $A_{i}$ are constants.

$e_{ij}$ is defined in a general way,

(4) $e_{ij}(x)=g_{ij}*f_{i}(x){\alpha}_{j}$

where $g_{i,j}$ are constants, $i >= j$ and $f_{i}(x) = ({\delta}ln(x^{h_{i}(x)}))/({\delta}ln(x))$.

If $h$ and $e$ were constants, $h$ would be an income elasticity as $x= Y/P_{m}$ and $e_{ij}$ would be own and cross price elasticity as $w_{i}= P_{i}/P_{m}$.

The following functional forms are chosen for $h_{s}$ and $h_{n}$,

(5) $h_{s}(x) = ({\lambda}/x)(1+({\kappa}/ln(x)))$

(6) $h_{n}(x) = {\nu}/(1-x)$.

In addition to the above, two other scaling parameters are applied when normalizing the demand values to that of materials. These are $psscl$ for staples and $pnscl$ for non-staples.

The parameters are fit using a weighted least square log likelihood function [@carrol1988] described below.

(7) $ln(L) = {\Sigma}_{i=1}^{N}(w_{i}(y_{i}-\hat{y_{i}})^{2})/2{\sigma^{2}}$

where, $y_{i}$ is the $i$th data value and $\hat{y_{i}}$ is the corresponding model output and $w_{i}$ is the weight assigned to the data point. Since the parameters were fit based on regional data, the regional population was used as the weight.

By applying the 11 parameters to the equations described above, the user can generate estimates of demand for staples and non-staple commodities in thousand calories  across different income levels and prices.     

# Main functions and customization

The [```ambrosia_vignette``](https://jgcri.github.io/ambrosia/articles/ambrosia_vignette.html) provides usable examples for all the major functions within the code. 

The `ambrosia` package can be easily loaded as a standard R package after installation from GitHub. The user can calculate demand for staples and non-staples using the  ```food.dmnd()``` function. The user will have to pass in a dataset with the price of staples ($Ps$), price of non-staples ($Pn$), incomes ($Y$) (Current income proxy used in `ambrosia` is GDP per capita in thousand USD). In addition to the dataset, the user must pass a vector of 11 calibration parameters. In order to functionalize the parameters, the code contains a function called ```vec2param()``` that will generate a parameter structure that can be used by the food demand function. The food demand function is implemented using equations (1), (2), (3) described above. The user can also calculate and analyze price elasticities using the function ```calc1eps()```. These elasticities are calculated in accordance with equation (4) described above.

An interactive version of the food demand model can be launched by the user through the ```runapp()``` function to explore the impact of different parameters.

One of the benefits of using ```ambrosia``` is that a user can estimate their own calibration parameters with a custom data set using the log-likelihood maximization approach. To enable this, ```ambrosia``` is equipped with a function ```create.dataset.for.parameter.fit()``` that will help a user generate a dataset that is appropriate for parameter estimation. The user can re-create the training data used to calculate the parameters for GCAM using the ```Process_Demand_Data.R``` under the ```scripts``` directory.

Users can complete the parameter estimation on the dataset returned by ```create.dataset.for.parameter.fit()``` with a call to the ```calculate.ambrosia.params()``` function. This function builds on the @edmonds2017global approach by maximizing the log-likelihood score using the ```optim()``` function. Note that the user can also choose to use a different method (for example , a  MCMC) to maximize the log-likelihood function by first setting up the function using the ```mc.setup()``` function. The code contains an example of a MCMC implementation in C++ under ```scripts/cpp ```.

# Acknowledgements

This research was supported by the U.S. Department of Energy, Office of Science, as part of research in MultiSector Dynamics, Earth and Environmental System Modeling Program. The Pacific Northwest National Laboratory is operated for DOE by Battelle Memorial Institute under contract DE-AC05-76RL01830. The views and opinions expressed in this paper are those of the authors alone.


# References
The paper can be compiled locally using the following command:

`pandoc paper.md --bibliography paper.bib -o paper_local.pdf`

See https://pandoc.org/installing.html to install `pandoc`.
Spaceholder for outputs from the vignettes.
This folder contains functions that can be used to calculate the training data used for generating parameter values for GCAM.    

-`Process_Demand_Data.R` was used to create the training dataset to calculate the current parameters for GCAM.

-`functions.R` and `FAOStat_commod_matching.R` contains various helper functions used for cleaning up FAO data

These scripts should be used for data processing by users intending to re-create or update the current parameters in GCAM.

# Example of fitting parameters using an MCMC implemented in C++

This code shows how we used to use RInside to run the `ambrosia` from a C++ program.  At one time this was how we used our C++
Markov Chain Monte Carlo code with the model, but these days we've
switched over to R-native MCMC packages. If someone really wants to run the model from C++, it
should be easy to rework this example to run with the new
version.

Using this code, the user can fit a log-likelihood function in ambrosia and then use the same to run a Markov Chain Monte Carlo parameter fit.
