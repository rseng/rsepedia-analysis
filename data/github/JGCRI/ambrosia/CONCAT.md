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
---
title: "`ambrosia`: User Tutorial"
author: "Kanishka Narayan"
date: "2021-03-09"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{`ambrosia`: User Tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Introduction
```r
#load libraries
library(ambrosia)
library(dplyr)
library(ggplot2)
library(knitr)
```
# Part 1: Calculating food demand, exploring demand side variables for a given set of parameters 

This section explains how a user can set up the parameter structure for the `ambrosia` and explore the basic demand side variables.


## Example 1.1: Calculate food demand

Demand in thousand calories can be calculated for different levels of income (`Y`) in thousand USD per capita for different prices of staples (`Ps`) and non-staples (`Pn`) both in $ per capita per day using the parameters calculated in example 1.1. 

```r
#Get a sample data set
Test_Data <- data.frame(Y=seq(0.1,30, by=0.1))

#Add sample values of Ps and Pn
Test_Data <- Test_Data %>% mutate(Ps=0.1,Pn=0.2)

#Calculate food demand
Food_Demand <- food.dmnd(Test_Data$Ps,Test_Data$Pn,Test_Data$Y)

```
## Example 1.2: Visualize food demand

The demand along with the price and income elasticities calculated using the `food.dmnd()` function  can be visualized easily.   

```r
#Add income and total demand to the data frame created in example 1.1 and create a plot.  
Food_Demand$Total_Demand <- Food_Demand$Qs+Food_Demand$Qn
Food_Demand$Y <- seq(0.1,30, by=0.1)

#Create the plot
g <- ggplot()+
    geom_line(data=Food_Demand,aes(x= Y,y= Qs,color= "Staple Demand"))+
    geom_line(data=Food_Demand,aes(x= Y,y= Qn,color= "Non Staple Demand"))+
    xlab("GDP per capita in thousand USD" )+
    ylab("Thousand calories")

#Plot
plot(g)    
```
![](example_3.png){width=50% height=50%} 

## Example 1.3: Visualize budget shares

Using the function above will create a dataframe with estimates of demand for each level of price and incomes and also the budget shares (shares of incomes spent) for staples and non-staples.

```r
#Rename the budget shares columns
budget_shares <- Food_Demand %>% 
                 mutate(share_of_staples= alpha.s*100,
                        share_of_non_staples= alpha.n*100,
                        share_of_materials= alpha.m*100)

#Add some income or Y values
budget_shares$Y <- seq(0.1,30, by=0.1)

#Create the plot
g<-ggplot()+
  geom_line(data=budget_shares,aes(x=Y,y=share_of_staples,color="Staple share"))+
  geom_line(data=budget_shares,aes(x=Y,y=share_of_non_staples,color="Non Staple share"))+
    xlab("GDP per capita in thousand USD" )+
    ylab("Percent of Income")

#Plot
plot(g)
```
![](example_4.png){width=50% height=50%} 

## Example 1.4: Calculate income elasticities

The demand code iteratively solves for the budget shares using a Broyden solver with changing incomes and re-calculates income and price elasticities for changes in budget shares.The user can separately analyze the income elasticities by using two functions (one for staples and other for non-staples) from within the parameter structure (a temporary parameter structure is generated in the example below)

For income elasticities, the functions generated within the parameter structure stored within the object `yfunc` return the elasticities themselves for a given level of income. Setting the optional parameter within these functions to TRUE would return the income term itself (i.e. `Y`^ `elas`). 

```r
#Note that setting the second argument to "TRUE" in the functions below would return the Y term (Y^elas) as opposed to the elasticity.

#Get the parameter structure

tmp_param <- vec2param()

# Income elasticity for staples
Food_Demand$eta.s <- tmp_param$yfunc[[1]](Y=Food_Demand$Y,FALSE)

# Income elasticity for non-staples
Food_Demand$eta.n <- tmp_param$yfunc[[2]](Y=Food_Demand$Y,FALSE)

#Create data for the plot
g<-ggplot()+
  geom_line(data=Food_Demand,aes(x=Y,y=eta.s,color="Income elasticities for staples"))+
  geom_line(data=Food_Demand,aes(x=Y,y=eta.n,color="Income elasticities for non-staples"))+
    xlab("GDP per capita in thousand USD" )+
    ylab("elasticity value")

#Plot
plot(g)
```
![](example_5.png){width=50% height=50%} 

## Example 1.5: Calculate price elasticities

Similar to the income elasticities, the user can also calculate and analyze price elasticities. These elasticities. The function `calc1eps()` can be used to recalculate price elasticities by passing the following parameters:

1)	different budget shares (alphas)

2)	income elasticities (calculated in example 1.4)

3)	A matrix of values for elasticities (this is a part of the parameter structure set up in example 1.4)

The price elasticities are calculated as an epsilion matrix with price elasticities for staples, non-staples and cross price elasticities. 

```r
#calc1eps returns a matrix of elasticities where the matrices are elasticities for staples(matrix 1), cross price #elasticities(matrix 2,3),non_staple elasticities(matrix 4))

#Get staple price and income elasticities
Food_Demand$staple_price_elasticity <- calc1eps(Food_Demand$alpha.s, Food_Demand$alpha.n , Food_Demand$eta.s , Food_Demand$eta.n, tmp_param$xi)[1:300]
Food_Demand$eta.n <- tmp_param$yfunc[[2]](Y=Food_Demand$Y,FALSE)

#Get non-staple price and income elasticities
Food_Demand$non_staple_price_elasticity <- calc1eps(Food_Demand$alpha.s,Food_Demand$alpha.n,Food_Demand$eta.s,Food_Demand$eta.n,tmp_param$xi)[901:1200]
Food_Demand$eta.n <- tmp_param$yfunc[[2]](Y=Food_Demand$Y,FALSE)

g<-ggplot()+
  geom_line(data=Food_Demand,aes(x=Y,y=staple_price_elasticity,color="Staple price elasticity"))+
  geom_line(data=Food_Demand,aes(x=Y,y=non_staple_price_elasticity,color="Non staple price elasticity"))+
    xlab("GDP per capita in thousand USD" )+
    ylab("elasticity value")

#Plot
plot(g)
```
![](example_6.png){width=50% height=50%} 

# Part 2: Understanding calibration parameters

## Example 2.1: Get parameters

As mentioned in the documentation, `ambrosia` can generate estimates of food demand for staples and non-staples in thousand calories for given set of prices and income levels.

In addition to the price and income parameters, `ambrosia` requires 11 calibration parameters for the model. These parameters can be set using the function `vec2param()`. In order to understand the 11 parameters, their units, ranges and how they can be calculated, the user can call,

```r
?vec2param()
```
To use `vec2param()` to introduce new calibration parameters in the food demand function, the user would have to call,

```r
Food_Demand <- food.dmnd(Test_Data$Ps,Test_Data$Pn,Test_Data$Y, params= vec2param()) 
```
The user can set the parameters by passing a 11 parameter vector directly to `vec2param()` or by setting the 11 parameters individually. The user can also directly call `vec2param()` to get the default parameter structure with default values.

Note that the current calibration parameters are calculated on the basis of national level data on food consumption and prices. These can effectively be calculated using new inputs thus allowing the user to calibrate the model to custom data.  

```r
#Use a sample vector of parameters. 
original_param_vector <- c(1.28,1.14,-0.19,0.21,-0.33,0.5,0.1,16,5.06,100,20) 

#These are the names of the parameters
parameter_names <- c('A_s', 'A_n', 'xi_ss', 'xi_cross', 'xi_nn', 'nu1_n',
                     'lambda_s', 'k_s', 'Pm', 'psscl','pnscl')

#Create a dataframe of parameters
parameter_data <- data.frame(parameter_names,original_param_vector)

#Return a parameter structure
tmp_param <- vec2param(original_param_vector)

#Return a table of parameters
kable(parameter_data ,col.names = c("parameter_name","value"),format = "pandoc")
```

This section now explains how a user can calculate model parameters from raw data. This way, `ambrosia` can be calibrated to custom data sets.  

# Part 3: Running the interactive version of the model

The user can also launch an interactive version of ambrosia to visually understand how the parameters are used to calculate the outputs and also to test new parameter values. Note that using the interactive version requires the `shiny` library to be installed.

The user can explore results by changing the main income and price parameters or test the impact of new values for the calibration parameters on the results.

## Example 3.1: Run the interactive version of the model

### Note: `shiny`is required to be installed before running the app using the function below.

```r
runapp()
```

# Part 4: Recalculating calibration parameters on the basis of new input data  

## Example 4.1: Get a dataset for parameter fitting from a sample training dataset

One of the benefits of using `ambrosia` is that a user can estimate their own parameters with a custom data set using the log-likelihood maximization approach. To enable this, `ambrosia` is equipped with a function `create.dataset.for.parameter.fit()` that will help a user generate a dataset that is appropriate for parameter estimation. Note: The user can re-create the training data used to calculate the parameters for GCAM using the `Process_Demand_Data.R` under the scripts folder.

There are a few steps that the function will perform on a sample dataset.

1)	It will ensure that the user’s dataset contains all columns required for parameter estimation

2)	It will filter out anomalies and outliers using parameterized cutoff values selected by the user. This step is necessary since data on food consumption and prices are often incomplete which may lead to unrealistically high or low values of consumption or prices in the dataset.

3)	After this, the function will create clusters of observations from the dataset based on income levels, and prices of staples and non-staples. This step is necessary because this being economic data, the variance (required for the likelihood function) can only be calculated within different clusters. The code will also check for a user specified minimum number of clusters(If there are anomalies within the dataset, the clustering can be incorrect leading in a small number of clusters). The clustering is implemented using the Divisive Analysis Clustering Algorithm (DIANA)

4) Once the clustering is completed, the code will calculate the variance which is the variance in food demand for staples (${\sigma^{2}}Q_{s}$) and non-staples ${\sigma^{2}}Q_{n}$ .Note that the user can chose a lower limit on the variance calculated. The default value of the lower limit is 0.01. 

In addition to a data frame, the function will return a CSV file output called “Processed_Data_for_MC.csv” that is stored in the outputs folder that will be used for the parameter estimation. The example below illustrates how to use the function on a raw dataset.

The example below illustrates how to use the function on a sample dataset 

```r
#Load the training data
data("training_data") 

#Create a sample dataset
sample_data<- training_data %>% 
   filter(year %in% c(2010:2015)) 

#Create a dataset for parameter fit
temp_data <- create.dataset.for.parameter.fit(data=sample_data,min_clusters = 20,min_price_pd = 20,min_cal_fd = 1700,outdir=tempdir()) 

```
## Example 4.2: Analyze distribution of variance for non-staples

The dataset returned by this function can now be used for parameter estimation. The user can also plot the variances for staples and non-staples to ensure there is a valid distribution and the data is not skewed. 

The below plots show the variance for pre-generated data. 

```r
#Load example processed_data 
data("processed_data_example")

#Create a plot using example loaded above or data generated in example 2.1
g<-ggplot()+
   geom_histogram(data=processed_data_example,aes(x=sig2Qn),bins = 40)+
   xlab("Variance within clusters for non-staples")+
   ggtitle("Distribution of variance for non-staples for selected sample dataset")

#Plot
plot(g)   
```
![](example_7.png){width=50% height=50%} 

## Example 4.3: Analyze distribution of variance for staples
```r
#Load example processed_data 
data("processed_data_example")

#Create a plot using example loaded above or data generated in example 2.1
g<-ggplot()+
   geom_histogram(data=processed_data_example,aes(x=sig2Qs),bins = 40)+
   xlab("Variance within clusters for staples")+
   ggtitle("Distribution of variances for staples for selected sample dataset")

#Plot   
plot(g)   
```
![](example_8.png){width=50% height=50%} 

## Example 4.4: Calculate actual parameters (option 1- Optimization)


Finally, the user can complete the parameter estimation on the dataset returned by `create.dataset.for.parameter.fit()` with a call to the `calculate.ambrosia.params()` function. ambrosia builds on the Edmonds et al.(2017) approach by maximizing the log-likelihood score using the `optim()` function. Note that the user can also choose to use a different method (for example, the original MCMC) to maximize the log-likelihood function by first setting up the log likelihood function using the `mc.setup()` function.  

The following steps are involved in the parameter estimation function,

1)	First a log-likelihood function is set up with the data returned by the function above. This function calculates a log likelihood score using a weighted least squares approach (where the population of the region is the weight) 

2)	Now, the value returned by this function will be maximized using `optim()`. The user can provide a seed of initial parameters to begin the optimization process (the lowest possible seed would be the lowest values of all 11 parameters). The default seed is set to the original parameters from Edmonds et al. The user can now specify the optimization method to be used. The default is set to the “BFGS” method, but the user can also run the optimization using methods such as “Neldor-Mead” etc.

3)	The function will now return a vector of parameters that can be used to derive estimates of food demand (Similar to Example 1 above). The function also prints out the maximized value of the log-likelihood function, so that the user can verify the efficiency and effectiveness of the parameter estimation.


### Warning: Calculating parameters takes a long time (around 30 minutes depending upon machine). Max iterations in the below example are set to 2 just for speed. 

```r
#load sample data
data("processed_data_example")

#write it to a temporary csv
write.csv(processed_data_example,"Processed_Data_for_MC.csv")

#calculate new parameters using this saved dataset
new_parameters <- calculate.ambrosia.params(datadir="Processed_Data_for_MC.csv",
                                            optim_method = "BFGS",
                                            original_param_vector= c(1.28,1.14,-0.19,0.21,-0.33,0.5,0.1,16,5.06,100,20),
                                            outdir = tempdir(),
                                            max_iterations = 2)
```

## Example 4.5: Calculate actual parameters (option 2- Set up a likelihood function and use any other method)

While `calculate.ambrosia.params()` is useful to quickly calculate a sample of parameters, the user can also use their own approaches to calculate parameter samples. For example, if a user wanted to create samples for the parameters using a MCMC as opposed to an optimization approach, `ambrosia` has inbuilt functions that allow a user to set up a likelihood function which can be maximized using any other method chosen by the user. The example below shows how the user can set up a likelihood function

```r

#The log likelihood function set up below can be used to fit parameters either using the maximization function or any #other approach like a Markov Chain Monte Carlo

#load sample data
data("processed_data_example")

#write it to a temporary csv
write.csv(processed_data_example,"Processed_Data_for_MC.csv")

#Now, set up a log likelihood function using the saved data. Note that the trace_param parameter is used 
#to trace parameter values during a fitting process like a MCMC calculation. See Part 4 below for more #details.

likelihood_function <- mc.setup("Processed_Data_for_MC.csv", trace_param= FALSE)


```
## Example 4.6: Calculate the log likelihood score based on the parameters (vector) calculated in Example 4.4 and the log likelihood function calculated in Example 4.5

```r
if(!(exists("new_parameters"))){

print("New parameters not calculated as shown in example 4.4. Loading pre-saved parameters")

new_parameters <- c(1.28,1.14,-0.19,0.21,-0.33,0.5,0.1,16,5.06,100,20)

}

#calculate the likelihood score using the function calculated in example 2.5
likelihood_score <- likelihood_function(c(new_parameters))

```

# Part 5: Exploring data from a parameter fitting process (MCMC or maximization or any other process) 

Now if a user calculated new parameters (either using `calculate.ambrosia.params()` or by setting up the log likelihood function using `mc.setup()` and maximizing the log likelihood using another method), `ambrosia` also allows the user to explore the results of a parameter fit.

## Example 5.1: Generating parameter data

Parameter data can be generated in two ways,

1) If `calculate.ambrosia.params()` is used to calculate parameter values, the user would have to set the `trace_param` to TRUE to generate the data (.dat) file

### Warning: Calculating parameters takes a long time (around 30 minutes depending upon machine). Max iterations in the below example are set to 2 just for speed. 

```r

new_parameters <- calculate.ambrosia.params(datadir="Processed_Data_for_MC.csv",
                                            optim_method = "BFGS",
                                            original_param_vector= c(1.28,1.14,-0.19,0.21,-0.33,0.5,0.1,16,5.06,100,20),
                                            outdir = tempdir(),
                                            trace_param = TRUE,
                                            max_iterations =2)

```


2) If the log-likelihood function was set up using `mc.setup()` directly and this was maximized using some other method, then the same `trace_param` parameter needs to be set to TRUE before the fitting process to ensure the .dat file is generated. 

```r

likelihood_function <- mc.setup("Processed_Data_for_MC.csv", trace_param = TRUE) 

```

Using the above code, the user can generate a .dat file which will track parameter samples during the fitting process along with the likelihood value itself. 

The user can now read in the file, assign parameter names and generate density plots for the same.  

## Example 5.2: Read in mc.data and create density plots
```r
#Load example data
data("mc_data_example")

#write it to a temporary file as an example
write.table(mc_data_example,"mc_data_example.dat")

#Read the table
#nparam is set to 11 since we are using the 11 
mc_example <- read.mc.data("mc_data_example.dat",varnames = namemc(nparam = 11))

#Create the density plot
mcparam.density(mc_example)

```
![](example_9.png){width=50% height=50%}

## Example 5.3: Re-create maximum probability density

The user can also extract the parameters that yield the highest log likehood score from amongst the samples.

```r
#Load example data
data("mc_data_example")

#write it to a temporary file as an example
write.table(mc_data_example,"mc_data_example.dat")

#Read the table
mc_example <- read.mc.data("mc_data_example.dat",varnames = namemc(nparam = 11))

#The function  below will return a vector of parameters which yield the highest log likelihood score
ml_parameters <- mcparam.ML(mc_example)

#Set up the log-likelihood function
func_MC <- mc.setup("Processed_Data_for_MC.csv")

#Now calculate the probability density.The value should be -962
probability_density <- func_MC(c(ml_parameters))
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand-plots.R
\name{make.demand.plot}
\alias{make.demand.plot}
\title{Plot staple, nonstaple, and total demand output from the model}
\usage{
make.demand.plot(alldata, xdata, xlabel, max.yval = NULL)
}
\arguments{
\item{alldata}{Data frame of output returned from \code{\link{food.dmnd}}.}

\item{xdata}{Vector of values for the x-axis.}

\item{xlabel}{Character string to use for the x-axis label}

\item{max.yval}{Maximum value to display on the y-axis.  See details.}
}
\value{
Plot of food demand by staples and non-staples relative to income
}
\description{
The plot will have quantity on the y-axis.  The x-axis data is passed in as a
separate vector, so it can be any of the input values used in the model (or
even in principle a variable that is only indirectly related to the input
values).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand.R
\name{lamks2nu1y0}
\alias{lamks2nu1y0}
\title{Convert the lambda and ks parameters to nu1 and y0}
\usage{
lamks2nu1y0(df)
}
\arguments{
\item{df}{Data frame of lambda and ks parameters.}
}
\value{
data frame with converted values
}
\description{
The staple income elasticity has two formulations of its parameters.  The one
that uses nu1 (income elasticity at Y==1) and y0 (value of Y for which the
elasticity is zero) is more intuitive, but it is incomplete.  That is, theere
are some valid models that simply cannot be expressed this way.  The model
can alternatively be expressed in terms of the lambda and ks parameters in
its equations.  This function allows an easy conversion from the latter
representation to the former.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand.R
\name{prepare.obs}
\alias{prepare.obs}
\title{Prepare observations for use in the model.}
\usage{
prepare.obs(obs)
}
\arguments{
\item{obs}{Data frame of observed data.}
}
\description{
Convert units on prices and rename columns with long, wordy names to make
them easier to work with.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcpar-analysis.R
\name{read.mc.data}
\alias{read.mc.data}
\title{Read the Monte Carlo results file}
\usage{
read.mc.data(filename, varnames = namemc())
}
\arguments{
\item{filename}{Name of the file where the results are stored}

\item{varnames}{Names of the variables stored in the file (the MC output
format didn't have column headers). If omitted, a default set of names is
supplied.}
}
\value{
A table with parameter names in accordance with the food demand model
}
\description{
Read the Monte Carlo results file
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand.R
\name{calc.hicks.actual}
\alias{calc.hicks.actual}
\title{Calculate the actual Hicks elasticities using the Slutsky equation.}
\usage{
calc.hicks.actual(eps, alpha.s, alpha.n, alpha.m)
}
\arguments{
\item{eps}{Elasticity values calculated by \code{\link{calc.elas.actual}}}

\item{alpha.s}{Budget fraction for staples}

\item{alpha.n}{Budget fraction for nonstaples}

\item{alpha.m}{Budget fraction for materials}
}
\value{
dataframe with hicks elasticities
}
\description{
Given actual price and income elasticities (i.e., calculated using finite
difference derivatives), compute the corresponding Hicks elasticities.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand.R
\name{compute.bias.corrections}
\alias{compute.bias.corrections}
\title{Compute regional bias corrections.}
\usage{
compute.bias.corrections(params, obs.trn)
}
\arguments{
\item{params}{Model parameter structure (described in
\code{\link{food.dmnd}})}

\item{obs.trn}{Data frame of training observations.}
}
\value{
calculated bias corrected values
}
\description{
Compute regional bias correction for a set of parameters and a training set of observations.
Bias correction factors are \emph{multiplied} by model data to get bias
corrected model output.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand-mc.R
\name{mc.setup}
\alias{mc.setup}
\title{Create a posterior log-pdf function for monte carlo sampling the model parameters.}
\usage{
mc.setup(
  obsdata_filename,
  logprior = NULL,
  logfile = NULL,
  chunksize = 10,
  trace_param = FALSE
)
}
\arguments{
\item{obsdata_filename}{File name for the oberved data.}

\item{logprior}{Optional function that takes a parameter vector and returns a
log prior probability density.  It need not (and should not) support the
log-posterior function's optional argument described in the details section.}

\item{logfile}{Optional file name for logging.}

\item{chunksize}{See details for description of chunking and why it is
needed.}

\item{trace_param}{Set this to TRUE to generate the data file (.dat file) during the statistical fitting procedure.}
}
\value{
A function that computes the log-posterior probability density for an
input vector of parameters.
}
\description{
The function returned takes a vector of parameter values as an argument and
returns a the log of the posterior pdf calculated by comparing this data to
the supplied observed data.  There is an optional second argument that is
described further below.  The log-post function should be usable in any MC
sampler.  It will also work in an optimizer, provided you configure it to
look for a maximum rather than a minimum.
}
\details{
The observed data should have the following columns: Ps, Pn, Y, Qs, Qn,
sigQs, sigQn.  There is a special provision to process the dataset produced
for GCAM

The demand system employs a nonlinear equation solver, which is used to
solve for the budget fractions that appear in the demand equations.  Although
structured as a system of equations, each time is actually independent of all
the others, but the solver has no way of knowing this.  Because the
complexity of the solver scales nonlinearly with the number of equations, we
break the times in the dataset up into chunks of modest size.

The optional second argument to the log-posterior function produced by this
function gives the number of parameter sets that are concatenated into a
single vector.  This capability is supplied to make life easier back when we
were calling these functions from a C code that used SSE/AVX instructions to
vectorize the sampler.  Samplers written in R can and should pretend that
this argument doesn't exist.
}
\section{TODO}{


\itemize{
\item{Provide an example dataset for the observed data.}
\item{Dispense with this chunk business and just run a 1-d solver in a loop.
What was I thinking?}
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand.R
\docType{data}
\name{testdata}
\alias{testdata}
\alias{y.vals}
\alias{Ps.vals}
\alias{Pn.vals}
\alias{Pm.vals}
\alias{samp.params}
\alias{x1}
\alias{x0}
\title{Sample values for the demand model.}
\format{
An object of class \code{numeric} of length 20.

An object of class \code{numeric} of length 20.

An object of class \code{numeric} of length 20.

An object of class \code{numeric} of length 20.

An object of class \code{list} of length 4.

An object of class \code{numeric} of length 9.

An object of class \code{numeric} of length 9.
}
\usage{
y.vals

Ps.vals

Pn.vals

Pm.vals

samp.params

x1

x0
}
\description{
These variables give examples of the data structures used in the model, as
well as baseline price, GDP and parameter values that should produce
reasonable results.

y.vals: Logarithmically-spaced pcGDP values

Ps.vals: Evenly spaced Ps values

Pn.vals: Evenly spaced Pn values

Pm.vals: Evenly spaced Pm values

samp.params: Example model parameter structure

x1: Monte carlo parameter vector.  This vector encodes the same model as
\code{samp.params}, but using the Monte Carlo version of the parameter
formulation for the staple income elasticity model.

x0: Parameters used to generate the test data
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcpar-analysis.R
\name{namemc}
\alias{namemc}
\title{Return the list of names for the parameters in the model.}
\usage{
namemc(nparam = 11)
}
\arguments{
\item{nparam}{Number of parameters.  Legal values are either 10 or 11.}
}
\value{
A vector of parameter names for the selected version of the model.
}
\description{
Supports both the 10 and 11 parameter version.  Adds the "LL" tag to the end to
cover the log likelihood column appened by the monte carlo code.
}
\details{
The 11 parameter version of the model is the canonical version (and the one
published in the paper)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcpar-analysis.R
\name{mcparam.density}
\alias{mcparam.density}
\title{Create a density plot for all of the MC variables}
\usage{
mcparam.density(mc.data)
}
\arguments{
\item{mc.data}{Data frame of Monte Carlo output.}
}
\value{
A density plot
}
\description{
Create a density plot for all of the MC variables
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand-mc.R
\name{mc.food.dmnd.byyear}
\alias{mc.food.dmnd.byyear}
\title{Compute food demand by year}
\usage{
mc.food.dmnd.byyear(obsdata, x, regions = NULL)
}
\arguments{
\item{obsdata}{Data frame of observed food demand}

\item{x}{Vector of model parameters in lambda-ks format.}

\item{regions}{Vector of regions to include in the computation.  If NULL,
include all.}
}
\value{
Data frame with results of food demand by year
}
\description{
The input parameters must be in the Monte Carlo formulation (i.e., lambda-ks,
\emph{not} nu1-y0 format
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand.R
\name{food.dmnd}
\alias{food.dmnd}
\title{Calculate food demand using the Edmonds, et al. model.}
\usage{
food.dmnd(Ps, Pn, Y, params = vec2param(), rgn = NULL)
}
\arguments{
\item{Ps}{Vector of staple food prices.}

\item{Pn}{Vector of nonstaple food prices.}

\item{Y}{Vector of per capita income.}

\item{params}{Model parameters structure (see details)}

\item{rgn}{Optional name for this calculation. If provided, the region will
be included as an extra column in the output data frame.}
}
\description{
The Edmonds model divides food consumption into two categories,
\emph{staples}, which represent basic foodstuffs, and \emph{nonstaples},
which represent higher-quality foods.  Demand for staples increases at low
income, but eventually peaks and begins to decline with higher income.
Demand for nonstaples increases with income over all income ranges; however,
total (staple + nonstaple) demand saturates asymptotically at high income.
}
\section{Arguments}{


Ps and Pn are food prices (staple food price and nonstaple price) in
international dollars per 1000 (dietary) calories.  Y is per-capita GDP in
international dollars. Ps, Pn, Y may be vectors but must all be the same
length.

Params is a structure:
\describe{
\item{xi}{2x2 array of the xi elasticities}
\item{A}{Leading coefficients in the quantity calculations}
\item{yfunc}{Length-2 list of functions giving Y^eta(Y) (see note below)}
\item{Pm}{Price of ``materials''. Loosely speaking, this parameter controls
how valuable food is relative to everything else in the economy.}
}

Note that we don't need elasticity parameters for the materials component
because we calculate Qm as a residual.  That is, whatever portion of
household budgets is not spent on food is by definition spent on materials.
}

\section{Output}{


The return value is a data frame with the following elements.
\describe{
 \item{Qs}{Quantity for staple foods (S)}
 \item{Qn}{Quantity for nonstaple foods (N)}
 \item{Qm}{Quantity for ``materials'' (M).  This quantity
represents an aggregate of everything else besides food that consumers buy.}
 \item{alpha.s}{Budget fraction for S}
 \item{alpha.n}{Budget fraction for N}
 \item{alpha.m}{Budget fraction for M}
}
Demand for staple and nonstaple foods are given in thousands of dietary
calories per person per day.  Units for materials demand are unspecified.
}

\section{Income Elasticity Functions}{


For one of the functional forms used for the income behavior, eta(Y), eta
blows up, but Y^(eta(Y)) is well behaved.  Therefore, the eta functions need
to be able to calculate not just eta(Y), but Y^(eta(Y)), so they can handle
the limiting cases.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcpar-analysis.R
\name{mcparam.itercount}
\alias{mcparam.itercount}
\title{Get the iteration count for a monte carlo dataset.}
\usage{
mcparam.itercount(niter, nproc, npset)
}
\arguments{
\item{niter}{Total number of iterations in the mcpar loop.}

\item{nproc}{Number of processors allocated to the calculation.}

\item{npset}{Number of parameter sets per processor}
}
\value{
iteration count for MCMC calculation
}
\description{
This is only necessary for the C++ MC code, due to some limitations in the
way it records its output.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand.R
\name{food.dmnd.byyear}
\alias{food.dmnd.byyear}
\title{Tabulate food demand by year for a model.}
\usage{
food.dmnd.byyear(obsdata, params, bc = NULL, region = NULL)
}
\arguments{
\item{obsdata}{Table of observed prices and incomes.}

\item{params}{Model parameter structure, as described in
\code{\link{food.dmnd}}}

\item{bc}{Vector of regional bias correction factors.  Optional. If
omitted, no bias correction will be applied.}

\item{region}{Region to apply the analysis to.}
}
\value{
food demand by year as a data frame
}
\description{
Tabulate food demand by year for in input model using the observed prices and
incomes for a given region.  If \code{region} is \code{NULL}, do it for all
regions and concatenate the result.
}
\details{
This function is intended to generate model predictions that can be compared
to observed consumption.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand.R
\name{calc.elas.actual}
\alias{calc.elas.actual}
\title{Calculate actual elasticities using numerical derivatives.}
\usage{
calc.elas.actual(Ps, Pn, Y, params, basedata = NULL)
}
\arguments{
\item{Ps}{Staple food prices}

\item{Pn}{Nonstaple food prices}

\item{Y}{Per-capita income}

\item{params}{Model parameters.  See description in \code{\link{food.dmnd}}.}

\item{basedata}{Model results for the base values of Ps, Pn, and Y.}
}
\value{
data frame with elasticities
}
\description{
Given a set of prices and incomes, and model parameters,calculate the
elasticities using numerical derivatives.  Optionally, you can pass the model
results for the base values, if you've already calculated them.
}
\details{
The inputs Ps, Pn, and Y can be vectors, but if they are, they must all be
the same length
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcpar-analysis.R
\name{mcparam.ML}
\alias{mcparam.ML}
\title{Get the maximum a-posteriori (MAP) parameters}
\usage{
mcparam.ML(mc.data)
}
\arguments{
\item{mc.data}{Data frame of Monte Carlo output}
}
\value{
A vector of parameters which yield the highest log likelihood.
}
\description{
The MAP parameters are those that produced the largest posterior probability
density.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_parameters.R
\name{calculate.ambrosia.params}
\alias{calculate.ambrosia.params}
\title{Calculates the 11 parameters for ambrosia using data calculated by \code{\link{create.dataset.for.parameter.fit}} by maximizing log likelihood}
\usage{
calculate.ambrosia.params(
  optim_method = "BFGS",
  original_param_vector = c(1.28, 1.14, -0.19, 0.21, -0.33, 0.5, 0.1, 16, 5.06, 100,
    20),
  datadir = "outputs/Processed_Data_for_MC.csv",
  outdir = "tests/testthat/test_outputs/",
  max_iterations = 100,
  print_progress = FALSE,
  trace_param = FALSE
)
}
\arguments{
\item{optim_method}{The optimization method to be used for maximization of the log likelihood. The default is set to BFGS}

\item{original_param_vector}{Original parameter vector to be used as the starting point for the optimization.These parameters are taken from Edmonds et al 2017}

\item{datadir}{Directory to the data calculated by \code{\link{create.dataset.for.parameter.fit}}}

\item{outdir}{Directory to store output csv. Default is set to test_output folder.}

\item{max_iterations}{A maximum number of iterations that can be passed to optim. This is largely meant for testing purposes.Default is set to 100 for BFGS.}

\item{print_progress}{A parameter that allows the user to track progress of function.}

\item{trace_param}{Setting this to TRUE will generate a .dat file with all parameter samples and log likelihhod value during a parameter fitting exercise (eg. MCMC)}
}
\value{
A vector with the 11 parameters which are as follows,

A_s : A scaling term for staple commodities

A_n : A scaling perm for nonstaples commodities

xi_ss : Price elasticity for staples

xi_cross : Cross price elasticity

xi_nn : Price elasticity for non-staples

nu_1n : Income elasticity for non-staples

lambda_s : Income elasticity for staples

k_s : Income level at which total demand begins saturating

Pm : Scaling parameter for price of materials (everything other than staples and non-staples)

psscl : Price scaling parameter for staples that is applied to the Price of staples (Ps) and the alpha of staples

pnscl : Price scaling parameter for staples that is applied to the Price of non-staples (Pn) and the alpha of non-staples
}
\description{
Calculates the 11 parameters for ambrosia using data calculated by \code{\link{create.dataset.for.parameter.fit}} by maximizing log likelihood
}
\details{
The following steps are involved in the parameter estimation function.

1) First a log-likelihood function is set up with the data.

2) Next, the value returned by the log-likelihood function will be maximized using optim(). The user can provide
 a seed of initial parameters to begin the optimization process (the lowest possible seed would be the lowest
 values of all 11 parameters). The default seed is set to the original parameters from Edmonds et al (2017). The user
 can now specify the optimization method to be used. The default is set to the "BFGS" method, but the user can
 also run the optimization using methods such as "Neldor-Mead".

3) Finally, the function will now return a vector of parameters that can be used to derive estimates of food demand (using the \code{\link{food.dmnd}} function). The function also prints out the maximized value of the log-likelihood function, so that the user can verify the efficiency and effectiveness of the parameter estimation.
}
\author{
KBN 2020
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand-mc.R
\docType{data}
\name{plohi}
\alias{plohi}
\title{Recommended parameter limits for the model.}
\format{
An object of class \code{matrix} (inherits from \code{array}) with 2 rows and 11 columns.
}
\usage{
plohi
}
\description{
Model parameters outside of these ranges are likely to produce models that
are numerically ill-behaved.  The first row of the matrix has recommended
minimum values; the second has recommended maximum values.
}
\details{
In other functions we talk about the "10-parameter" and "11-parameter" versions
of the model.  The 10-parameter version was used in early tests and is now
deprecated; therefore, these recommendations apply only to the 11-parameter
version.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui-fcns.R
\name{guiwidgets}
\alias{guiwidgets}
\alias{xi.matrix.input}
\alias{y0.input.box}
\alias{column.input.table}
\alias{eta.selector}
\title{Input widgets for food demand GUI app}
\usage{
xi.matrix.input()

y0.input.box()

column.input.table(
  inputids,
  defvals,
  min,
  max,
  step,
  labels = c("Staple (A_s)", "Nonstaple (A_n)")
)

eta.selector(id, label2 = "\\\\(\\\\eta=f(Y)\\\\)", sel = 1)
}
\arguments{
\item{inputids}{Character vector of parameter identifiers}

\item{defvals}{Vector of default values}

\item{min}{Minimum value}

\item{max}{Maximum value}

\item{step}{Slider step size}

\item{labels}{Character vector of labels}
}
\description{
Input widgets for food demand GUI app
}
\section{Functions}{
\itemize{
\item \code{xi.matrix.input}: Input boxes for xi matrix

\item \code{y0.input.box}: Input boxes for y0

\item \code{column.input.table}: Draw input grid for other inputs

\item \code{eta.selector}: Draw selector widget for eta
}}

\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand-plots.R
\name{mc.make.byyear.plot}
\alias{mc.make.byyear.plot}
\title{Make the by-year plot for a set of monte carlo results by sampling the distribution}
\usage{
mc.make.byyear.plot(
  mc.data,
  obsdata,
  bias.correction = NULL,
  region = NULL,
  nsamp = 30,
  pltrgn = NULL
)
}
\arguments{
\item{mc.data}{Monte Carlo results data}

\item{obsdata}{Data frame of observed food demand data}

\item{bias.correction}{Regional bias correction factors (default = none)}

\item{region}{Vector of regions to plot (default = all)}

\item{nsamp}{Number of samples to draw from the Monte Carlo distribution}

\item{pltrgn}{Regions to include in the plot.  If \code{NULL}, include them
all.}
}
\value{
by-year plot for a set of monte carlo results by sampling the distribution
}
\description{
Make the by-year plot for a set of monte carlo results by sampling the distribution
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand.R
\name{calc1q}
\alias{calc1q}
\title{Calculate demand quantities for a single set of inputs}
\usage{
calc1q(Ps, Pn, Y, eps, Ysterm, Ynterm, Acoef, psscl, pnscl)
}
\arguments{
\item{Ps}{Price of staple foods}

\item{Pn}{Price of nonstaple foods}

\item{Y}{Per-capita income}

\item{eps}{Matrix of Marshall elasticities}

\item{Ysterm}{Income term in the demand equation for staples}

\item{Ynterm}{Income term in the demand equation for nonstaples}

\item{Acoef}{Leading multiplier parameter.}

\item{psscl}{Price scaling parameter for Staple food products}

\item{pnscl}{Price scaling parameter for non-staple food products}
}
\value{
Quantities of staples and non-staples (Qs and Qn)
}
\description{
This is a helper function for the \code{mapply} calls in
\code{\link{food.dmnd}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand-plots.R
\name{make.byyear.plot}
\alias{make.byyear.plot}
\title{Plot model results by year}
\usage{
make.byyear.plot(byyear.data, pltrgn = NULL)
}
\arguments{
\item{byyear.data}{Output from \code{\link{food.dmnd.byyear}}}

\item{pltrgn}{Region to plot}
}
\value{
Plot of model results by year
}
\description{
Plot model results by year
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand-plots.R
\name{make.byincome.plot}
\alias{make.byincome.plot}
\title{Plot model results by per-capita income}
\usage{
make.byincome.plot(obsdata, params, region = NULL)
}
\arguments{
\item{obsdata}{Data frame of observed food demand data}

\item{params}{Model parameter structure}

\item{region}{Regions to include in the plot.}
}
\value{
Plot of model results by income levels
}
\description{
Plot model results by per-capita income
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand.R
\name{eta.constant}
\alias{eta.constant}
\alias{eta.s}
\alias{eta.n}
\title{Generate an income elasticity function with constant income elasticity.}
\usage{
eta.constant(eta0)

eta.s(nu1, y0, mc.mode = FALSE)

eta.n(nu1)
}
\arguments{
\item{eta0}{Value of the constant elasticity.}

\item{nu1}{Income elasticity at Y=1}

\item{y0}{Value of Y for which elasticity = 0 (this is generally the peak of
the curve).}

\item{mc.mode}{If true, then treat the first two parameters not as nu1 and
y0, but as direct specifications of k and lambda.  This flag is necessary
becaue nu1 and y0, while more intuitive to work with, are not a complete
specification of the parameter space. (That is, there are valid models that can only
be specified in terms of k and lambda.)}
}
\value{
A function suitable for use as either \code{eta.s} or \code{eta.n} in
the food demand model.
}
\description{
These functions generate a function that can be used as the income elasticity
functions for staple or nonstaple foods.  Staple and nonstaple foods each
have their own functional forms, given by \code{eta.s} and \code{eta.n}.
There is also a constant elasticity function that can be used for testing,
though it is not formally part of the model design.
}
\details{
#' Income elasticity functions have the following signature:
\itemize{
 \item{\code{function(Y, calcQ)}}
}
where Y is the per-capita GDP, and calcQ is a flag indicating whether the
function should calcluate the income quantity term or the income elasticity.
}
\section{Functions}{
\itemize{
\item \code{eta.s}: Generate an income elasticity function for staple foods.

\item \code{eta.n}: Generate an income elasticity function for nonstaple foods.
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{runapp}
\alias{runapp}
\title{Launch the interactive GCAM food demand model}
\usage{
runapp()
}
\description{
The interactive model allows users to adjust model parameters and see how
they affect the model output.  To run the interactive model you must have the
R "shiny" package installed.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand-mc.R
\name{vec2param}
\alias{vec2param}
\title{Convert a vector of parameters into a params structure.}
\usage{
vec2param(
  x = NULL,
  A_s = 1.28,
  A_n = 1.14,
  xi_ss = -0.19,
  xi_cross = 0.21,
  xi_nn = -0.33,
  nu1_n = 0.5,
  lambda_s = 0.1,
  k_s = 16,
  Pm = 5.06,
  psscl = 100,
  pnscl = 20
)
}
\arguments{
\item{x}{Vector of 11 model parameters using the statistical formulation. If this is NULL, a vector is generated based on default parameter values for 11 parameters.}

\item{A_s}{A scale term used to derive expenditure share for staple demand.}

\item{A_n}{A scale term used to derive expenditure share for non-staple demand.}

\item{xi_ss}{Price elasticity of staple goods. Unit change in per capita demand for staples as a result of unit increase in price (in $ per person per day).}

\item{xi_cross}{Cross price elasticity between staples and non-staples which in combination with the other price elasticities is used to derive substitution elasticity.}

\item{xi_nn}{Price elasticity of non-staple goods. Unit change in per capita demand for non-staples as a result of unit increase in price (in $ per person per day).}

\item{nu1_n}{Income elasticity for non-staple goods. Unit change in per capita demand for non-staples for unit change in income (in thousand USD).}

\item{lambda_s}{Income elasticity for staple goods. Unit change in per capita demand for staples for unit change in income (in thousand USD).}

\item{k_s}{Exponent of Income level at which staple demand is anticipated to be at its highest. Log of this term (log(k_s)) is the income level in Thousand USD.k_s and lamda_s are used in conjunction to derive the income elasticity for staples.}

\item{Pm}{Price of materials in $ per capita per day. This is basically the price consumers pay for goods in the economy others than food. The range based on the statistical formulation is 1.94 - 5.91.}

\item{psscl}{Additional scaling term used to derive the expenditure shares for staples. This is applied to price of staples (Ps/Pm * psscl), where Ps is the price of staples and Pm is the price of materials and to the expenditure shares of staples (alpha_s).}

\item{pnscl}{Additional scaling term used to derive the expenditure shares for non-staples. This is applied to price of non-staples (Pn/Pm * pnscl), where Pn is the price of non-staples and Pm is the price of materials and to the expenditure shares of non-staples (alpha_n).}
}
\value{
Parameter structure suitable for use in \code{\link{food.dmnd}}.
}
\description{
This function allows the user to pass in a vector of each of the 11 parameters for ambrosia.
The user can set individual parameter values described below or can pass in a direct vector of 11 parameters.
We also look at the number of
parameters passed in.  If it is 10, we assume you want a constant income elasticity for staples (etas = constant).  If
it's 11, we assume you want to use a dynamic income elasticity for staples (etas = eta.s(lambda, k)).  If it's anything else,
we throw an error. The table in the details section provides a description of the parameters, the units and the acceptable range.
}
\details{
The documentation provides an explanation of the parameters along with acceptable ranges for the parameters.
The default parameters are the parameter values that yielded the maximum log likelihood values duing the statistical fitting procedure.
The range of the parameters is derived by filtering the range of parameters in the MCMC for the 95th percentile confidence interval.
The parameters are not independent of each other.
Note that xi_cross is used for both xi_sn and xi_ns, forcing them to be equal.

If there are only 10 parameters, then the first 8 are as above,
and the next to last is eta_s.

The parameters in the vector are described in the table below:\tabular{lllll}{
   Parameter name \tab Description \tab Units \tab Value \tab Range of parameter values (Based on 95th percentile confidence interval) \cr
   \code{A_s} \tab A scale term used to derive expenditure share for staple demand \tab Unitless \tab 1.28 \tab 1.25 -1.40 \cr
   \code{A_n} \tab A scale term used to derive expenditure share for non-staple demand. \tab Unitless \tab 1.14 \tab 0.9 - 1.16 \cr
   \code{xi_ss} \tab Price elasticity of staple goods. Unit change in per capita demand for staples as a result of unit increase in price \code{Ps} (in $ per person per day). \tab Elasticity \tab -0.19 \tab -0.27 - -0.07 \cr
   \code{xi_cross} \tab Cross price elasticity between staples and non-staples which in combination with the other price elasticities is used to derive substitution elasticity. \tab Elasticity \tab 0.21 \tab 0.09 - 0.27 \cr
   \code{xi_nn} \tab Price elasticity of non-staple goods. Unit change in per capita demand for non-staples as a result of unit increase in price \code{Pn} (in $ per person per day). \tab Elasticity \tab -0.3 \tab -0.46 - -0.10 \cr
   \code{nu1_n} \tab Income elasticity for non-staple goods. Unit change in per capita demand for non-staples for unit change in income \code{Y} (in thousand USD). \tab Elasticity \tab 0.5 \tab 0.46 - 0.61 \cr
   \code{lambda_s} \tab Income elasticity for staple goods. Unit change in per capita demand for staples for unit change in income \code{Y} (in thousand USD) \tab Elasticity \tab 0.1 \tab 0.075 - 0.16 \cr
   \code{k_s} \tab Exponent of Income level at which staple demand is anticipated to be at its highest \tab Thousand USD \tab 16 \tab 10 -17 \cr
   \code{Pm} \tab Price of materials (everything else in the economy other than food products) \tab $ per day \tab 5 \tab 2 -6 \cr
   \code{psscl} \tab Additional scaling term used to derive the expenditure shares for staples. This is applied to price of staples (\code{Ps}/\code{Pm} * \code{psscl}), where Ps is the price of staples and Pm is the price of materials and to the expenditure shares of staples (\code{alpha_s}). \tab Unitless \tab 100 \tab 80 - 120 \cr
   \code{pnscl} \tab Additional scaling term used to derive the expenditure shares for non-staples. This is applied to price of non-staples (\code{Pn}/\code{Pm} * \code{pnscl}), where Pn is the price of non-staples and Pm is the price of materials and to the expenditure shares of non-staples (\code{alpha_n}). \tab Unitless \tab 20 \tab 18 - 25 \cr
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand.R
\name{calc1eps}
\alias{calc1eps}
\title{Calculate the exponents in the demand equation.}
\usage{
calc1eps(alpha.s, alpha.n, eta.s, eta.n, xi)
}
\arguments{
\item{alpha.s}{Budget fraction for staple foods.}

\item{alpha.n}{Budget fraction for nonstaple foods.}

\item{eta.s}{Income elasticity for staple foods.}

\item{eta.n}{Income elasticity for nonstaple foods.}

\item{xi}{Matrix of Hicks elasticities}
}
\value{
matrix of price elasticities, where xi[1],xi[4] are staple and non-staple price elasticities, rest are cross elasticities
}
\description{
This is an approximation to the Slutsky equations, inasmuch as we use these
as the exponents directly, instead of solving for exponents that produce
these values as the elasticities.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui-fcns.R
\docType{data}
\name{guiconst}
\alias{guiconst}
\alias{xidefault}
\alias{elasmin}
\alias{elasmax}
\alias{elasstep}
\alias{etastep}
\alias{spacer}
\title{Constants used in the GUI}
\format{
An object of class \code{numeric} of length 3.

An object of class \code{numeric} of length 1.

An object of class \code{numeric} of length 1.

An object of class \code{numeric} of length 1.

An object of class \code{numeric} of length 1.

An object of class \code{html} (inherits from \code{character}) of length 1.
}
\usage{
xidefault

elasmin

elasmax

elasstep

etastep

spacer
}
\description{
Constants used in the GUI

xidefault : Default xi values

elasmin : Minimum elasticity for the slider bar

elasmax : Maximum elasticity for the slider bar

elasstep : Elasticity bar step size

etastep : Eta slider bar step size

spacer : Spacer for labels
}
\keyword{datasets}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{calc.pop.weight}
\alias{calc.pop.weight}
\title{Calculate a weight factor based on the population.}
\usage{
calc.pop.weight(input.data)
}
\arguments{
\item{input.data}{Data frame of food demand input data from FAO}
}
\value{
Data frame of input data with a population weight column added.
}
\description{
These weight factors can be used to give high population regions more
influence in the model fit relative to low population regions.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand.R
\name{food.dmnd.byincome}
\alias{food.dmnd.byincome}
\title{Tabulate food demand by per-capita-income}
\usage{
food.dmnd.byincome(obsdata, params, region = NULL)
}
\arguments{
\item{obsdata}{Data frame of observed prices and incomes.}

\item{params}{Model parameter structure.  See notes in
\code{\link{food.dmnd}}.}

\item{region}{Name of a single region to tabulate.  If \code{NULL}, tabulate all
regions and concatenate the tables.}
}
\value{
food demand by per capita income returned as a data frame.
}
\description{
Tabulate staple and nonstaple demand as a function of per-capita GDP.  This
function uses observed prices for the comparison
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_parameters.R
\name{create.dataset.for.parameter.fit}
\alias{create.dataset.for.parameter.fit}
\title{Create dataset with observational error for log likelihood calculation using clustering.}
\usage{
create.dataset.for.parameter.fit(
  min_price_pd = 20,
  min_cal_fd = 1000,
  min_clusters = 300,
  lower_limit_sigma = 0.01,
  data = NULL,
  outdir = "tests/testthat/test_outputs/",
  print_progress = FALSE
)
}
\arguments{
\item{min_price_pd}{Minimum price paid for non-staples.}

\item{min_cal_fd}{Minimum calories for the food demand model}

\item{min_clusters}{Minimum number of clusters to be generated by clustering algoritm. It is recommended to not lowewr this parameter below 20.}

\item{lower_limit_sigma}{Lower limit for sigma values calculated}

\item{data}{A data.frame or data.table with the raw data. Data should contain following names,

{s_cal_pcap_day_thous} (Containing 1000 calories per capita per day for staples)

{ns_cal_pcap_day_thous} (Containing 1000 calories per capita per day for non-staples)

{gdp_pcap_thous} (Containing GDP per capita)

{s_usd_p1000cal} (Price of 1000 calories for staples per person per day)

{ns_usd_p1000cal} (Price of 1000 calories for non-staples per person per day)}

\item{outdir}{Directory to store output csv. Default is set to test_output folder.}

\item{print_progress}{A parameter that allows the user to track progress of function.}
}
\value{
A dataframe  called Processed_Data_for_MC with the following columns

sig2Qn- which is the observational error for non-staples

sig2Qs- which is the observational error for staples
}
\description{
There are a few steps that the function will perform on a sample dataset (see details).
}
\details{
The steps that will be performed by the function are :

1) It will ensure that the user's dataset contains all columns required for parameter estimation.

2) It will filter out anomalies and outliers using parameterized cutoff values selected by the user.
  This step is necessary since data on food consumption and prices are often incomplete which may lead to
  unrealistically high or low values of consumption or prices in the dataset.

3) After this, the function will create clusters of observations from the dataset based on income levels,
  and prices of staples and non-staples. This step is necessary because this being economic data,
  the observational error can only be calculated within different clusters. The code will also check for a  user
  specified minimum number of clusters (if there are anomalies within the dataset, the clustering can be incorrect leading in a small number of clusters).
  The clustering is implemented using the Divisive Analysis Clustering Algorithm (DIANA).

4) Once the clustering is completed, the code will calculate the observational error
  which is the variance in food demand for staples and non-staples .Note that the user can chose a
  lower limit on the observational error calculated. The default value of the lower limit is 0.01.
}
\author{
KBN 2020
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{assign.sigma.Q}
\alias{assign.sigma.Q}
\title{Assign observational errors to observed demand quantities}
\usage{
assign.sigma.Q(input.data, min.group = 5)
}
\arguments{
\item{input.data}{Data frame of observational input.}

\item{min.group}{Minimum group size for clustering}
}
\value{
A dataframe updated with the calculated observational error (sig2Qs and sig2Qn columns)
}
\description{
Assign sigma (observational error) values for Qs and Qn in an input data set.
We do this by clustering the input on Ps, Pn, and Y and then taking the
variance of the Qs and Qn in each cluster.
}
\details{
Observational errors are estimated by clustering observed data by the demand
model input values (i.e., prices for staple and nonstaple foods) and
calculating the variance of observations in each cluster.  In order for this
to work, you have to ensure that there are enough observations in each
cluster to produce a variance that is at least somewhat reliable.  The
tradeoff here is that the larger you make the clusters, the better the
variance estimate is, but less alike the observations in the cluster actually
are (meaning some of the variance is not observational error, but actual
difference in demand.  This tradeoff is controlled by the \code{min.group}
argument.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{recursive.partition}
\alias{recursive.partition}
\title{Partition input data into clusters with a minimum number of members}
\usage{
recursive.partition(input.data, cluster.vars, min.members = 5)
}
\arguments{
\item{input.data}{Data frame of data to be clustered.}

\item{cluster.vars}{Vector of names of variables to use in the clustering.}

\item{min.members}{Desired minimum number of members per cluster.}
}
\value{
List of data frames, one data frame for each cluster.
}
\description{
Partition the input data into clusters with a given minimum number of members
per cluster.  This turns out to be kind of hard to do because the structure
produced by the clustering algorithm isn't conducive to finding the cluster
that a too-small cluster was split from.  Instead, we find a dendrogram cut
that produces all clusters larger than the minimum.  Then we recursively
partition all of the clusters that are large enough that they could be split
in two.
}
\details{
The purpose of this function is to allow us to estimate the observational
error in historical food demand observations.  By clustering observations
that have similar input values (prices, GDP) we can get something
approximating repeated measurements of similar situations.  Setting a minimum
number of members per cluster allows us to have enough measurements per
grouping to get a resasonable estimate of the variance.

This function returns a list of data frames, with each list item being a
single cluster.  This irretrievably scrambles the order of the rows, so if
recovering the original order is important, include an ID column.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand.R
\name{apply.bias.corrections}
\alias{apply.bias.corrections}
\title{Apply bias corrections to model outputs}
\usage{
apply.bias.corrections(mod, bc)
}
\arguments{
\item{mod}{Data frame of model outputs}

\item{bc}{Vector of regional bias correction factors}
}
\value{
A data frame with bias corrected values
}
\description{
Apply bias corrections to model outputs
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcpar-analysis.R
\name{mcparam.clip.tails}
\alias{mcparam.clip.tails}
\title{Filter a Monte Carlo distribution by quantiles}
\usage{
mcparam.clip.tails(mc.data, qlo = 0.01, qhi = 0.99)
}
\arguments{
\item{mc.data}{Data frame of Monte Carlo output}

\item{qlo}{Lower quantile for filtering}

\item{qhi}{Upper quantile for filtering.}
}
\value{
Logical vector with \code{TRUE} for rows that are in the main body of
the distribution, \code{FAlSE} for those that aren't.
}
\description{
This filters on all parameter simultaneously, so any sample that is outside
the requested quantile range in the marginal distribution of any parameter
will be excluded.
}
\details{
The log-posterior column is handled a bit differently; it isn't filtered on
the high end, just the low end.

The purpose of this function was to trim values far out on the tails that
made the plot scales impossible to read.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/food-demand.R
\name{merge.trn.tst}
\alias{merge.trn.tst}
\title{Create a merged dataset with training and test data, each labeled accordingly}
\usage{
\method{merge}{trn.tst}(obs.trn, obs.tst)
}
\arguments{
\item{obs.trn}{Data frame of observations in the training set.}

\item{obs.tst}{Data frame of observations in the testing set.}
}
\value{
merged dataframe
}
\description{
Create a merged dataset with training and test data, each labeled accordingly
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcpar-analysis.R
\name{mcparam.sample}
\alias{mcparam.sample}
\title{Sample the MC results using bootstrap sampling.}
\usage{
mcparam.sample(mc.data, nsamp = 100, func = NULL)
}
\arguments{
\item{mc.data}{Data frame of Monte Carlo output}

\item{nsamp}{Number of samples to draw}

\item{func}{Optional function to apply to the samples drawn from the MC
distribution.}
}
\value{
Data frame or list (see details)
}
\description{
Optionally, apply a function to the sampled values.
}
\details{
If \code{func} is \code{NULL}, the return value will be a data frame of
sampled parameter values.  Otherwise the return value will be a list with the
data frame just described in the first element and the output of \code{func}
in the second.
}
