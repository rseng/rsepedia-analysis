[![License](https://img.shields.io/github/license/RETURN-project/BenchmarkRecovery)](https://choosealicense.com/licenses/apache-2.0/)
[![Build Status](https://github.com/RETURN-project/BenchmarkRecovery/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/RETURN-project/BenchmarkRecovery/actions)
[![codecov](https://codecov.io/gh/RETURN-project/BenchmarkRecovery/graph/badge.svg)](https://codecov.io/gh/RETURN-project/BenchmarkRecovery)
[![codecov](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4320502.svg)](https://doi.org/10.5281/zenodo.4320502)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8B-orange)](https://fair-software.eu)

# Benchmarking recovery indicators derived from remote sensing time series

This project simulates Landsat data and evaluates the performance of recovery indicators with respect to data and disturbance characteristics.

## Background

The context of this project is the study of the recovery of tropical forests after an abrupt disturbance (typically a forest fire) using satellite images as a data source.

The speed of recovery after a disturbance is known to be correlated with the concept of resilience. This is true not only for forests, but for many dynamical systems. To put it simply: forests that recover fast are more resilient. Forests that recover slowly may be in danger of permanent disappearance.

The specialized literature proposes different metrics for measuring the recovery speed. The performance of these metrics depends on many factors. Some of them are natural, such as the intensity of the perturbation or the seasonality. Others are technical, such as the sampling frequency or the spatial resolution.

## Purpose


The purpose of this project is to **efficiently** **compare** the  **reliability** of different post-disturbance recovery **metrics**.

## Mechanics

1. **Infers** time series' **parameters** and characteristics from optical satellite image data​
2. Uses those parameters to **create** a large collection of synthetic (but realistic) **time series**
3. **Calculates** several state-of-the-art recovery **metrics**

## Simplified workflow

![Simplified workflow](./img/flow.png)

## Install

You can install the master version from R via:

```r
library(devtools)
install_github("RETURN-project/BenchmarkRecovery")
```

## Citation info

Please cite as:
```
Wanda De Keersmaecker, & Pablo Rodríguez-Sánchez. (2020, December 14). RETURN-project/BenchmarkRecovery: Benchmarking recovery metrics derived from remote sensing time series (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.432050
```
Citation can be exported to different formats (BibTeX, JSON, ...) with our [Zenodo link](https://zenodo.org/record/4320503#.X9dFkFOYWhc).

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

[Apache](https://choosealicense.com/licenses/apache/)
---
title: "Plot recovery indicators"
author: "Wanda De Keersmaecker"
date: "6/10/2020"
output: html_document
vignette: >
   %\VignetteIndexEntry{Plot recovery indicators}
   %\VignetteEngine{knitr::rmarkdown}
   %\usepackage[utf8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

# Load libraries
library(BenchmarkRecovery)
```

## Illustrate recovery indicators for yearly time series

```{r}
# Let's generate yearly time series 
ts1 <- -2:5
ys1 <- c(6,5,1,1.5,2,3,3,4)
ys2 <- c(7,6,1,1.5,2,3,3,4)
ys3 <- c(7,6,1,1.5,2,3,3,4)+2

tpert <- 0
ts_pre <- c(-2,-1)
ts_post <- c(4,5)

# calculate recovery indicators
rri_y1 <- rri(ts1, ys1, tpert = tpert, ts_pre = ts_pre, ts_post = ts_post)
r80p_y1 <- r80p(ts1, ys1, ts_pre = ts_pre, ts_post = ts_post)
yryr_y1 <- yryr(ts1,ys1,tpert,deltat = 5)

rri_y2 <- rri(ts1, ys2, tpert = tpert, ts_pre = ts_pre, ts_post = ts_post)
r80p_y2 <- r80p(ts1, ys2, ts_pre = ts_pre, ts_post = ts_post)
yryr_y2 <- yryr(ts1,ys2,tpert,deltat = 5)

rri_y3 <- rri(ts1, ys3, tpert = tpert, ts_pre = ts_pre, ts_post = ts_post)
r80p_y3 <- r80p(ts1, ys3, ts_pre = ts_pre, ts_post = ts_post)
yryr_y3 <- yryr(ts1,ys3,tpert,deltat = 5)


# plot the results
rec1 <- paste0('RRI = ', round(rri_y1, digits = 3), ', R80p = ', round(r80p_y1, digits = 3), ', YrYr = ', round(yryr_y1, digits = 3))
rec2 <- paste0('RRI = ', round(rri_y2, digits = 3), ', R80p = ', round(r80p_y2, digits = 3), ', YrYr = ', round(yryr_y2, digits = 3))
rec3 <- paste0('RRI = ', round(rri_y3, digits = 3), ', R80p = ', round(r80p_y3, digits = 3), ', YrYr = ', round(yryr_y3, digits = 3))

plot(ts1, ys1, 'o', xlab = 'Year', ylab = 'Response', ylim = c(min(c(ys1,ys2,ys3)), max(c(ys1,ys2,ys3))))
lines(ts1, ys2, 'o', col = 'red')
lines(ts1, ys3, 'o', col = 'orange')
legend('topright', legend = c(rec1, rec2, rec3),
       col = c('black', 'red', 'orange'), lty = 1, cex = 1)
```

---
title: "Plotting decay functions"
author: "Pablo Rodriguez-Sanchez"
date: ""
output: html_document
vignette: >
   %\VignetteIndexEntry{Plotting decay functions}
   %\VignetteEngine{knitr::rmarkdown}
   %\usepackage[utf8]{inputenc}
---

```{r setup, include=FALSE}
# Default chunk configuration
knitr::opts_chunk$set(echo = TRUE, 
                      eval = TRUE, 
                      fig.width=8, 
                      fig.height=6) 

# Load libraries
library(BenchmarkRecovery)
```

## Comparing decay functions

```{r simulate}
# Set the problem
## Set the times
ts <- seq(10, 100, by = 1)

## Set the parameters (the same for all)
offset <- 0
pert <- 5
tpert <- 15
thalf <- 3
noise <- 0.1

# Generate the time series
ys_p <- piecewise(ts, offset, pert, tpert, thalf, noise)
ys_e <- exponential(ts, offset, pert, tpert, thalf, noise)
ys_r <- realistic(ts, offset, pert, tpert, thalf, noise)
```

### Plots

#### Complete time series

```{r plot, echo=FALSE}
# Plot the time series
plot(ts, ys_p, 'l', col = 'black', xlab = 'time', ylab = 'state')
lines(ts, ys_e, col = 'red')
lines(ts, ys_r, col = 'blue')
legend('topright', legend = c('piecewise', 'exponential', 'realistic'),
       col = c('black', 'red', 'blue'), lty = 1, cex = 1)
```

#### Detail (perturbation)

```{r plot-pert, echo=FALSE}
# Plot the time series
plot(ts, ys_p, 'l', col = 'black', xlab = 'time', ylab = 'state', xlim = c(tpert - thalf, tpert + thalf))
lines(ts, ys_e, col = 'red')
lines(ts, ys_r, col = 'blue')
legend('topright', legend = c('piecewise', 'exponential', 'realistic'),
       col = c('black', 'red', 'blue'), lty = 1, cex = 1)
```

#### Detail (transient)

```{r plot-transient, echo=FALSE}
# Plot the time series
plot(ts, ys_p, 'l', col = 'black', xlab = 'time', ylab = 'state', xlim = c(tpert, tpert + 7*thalf))
lines(ts, ys_e, col = 'red')
lines(ts, ys_r, col = 'blue')
legend('topright', legend = c('piecewise', 'exponential', 'realistic'),
       col = c('black', 'red', 'blue'), lty = 1, cex = 1)
```

#### Detail (asymptotics / long term dynamics)

```{r plot-asymptotic, echo=FALSE}
# Plot the time series
plot(ts, ys_p, 'l', col = 'black', xlab = 'time', ylab = 'state', xlim = c(0.5*max(ts), max(ts)), ylim = c(-4*noise, 4*noise))
lines(ts, ys_e, col = 'red')
lines(ts, ys_r, col = 'blue')
legend('topright', legend = c('piecewise', 'exponential', 'realistic'),
       col = c('black', 'red', 'blue'), lty = 1, cex = 1)
```

### Compare standard deviations

```{r}
sds <- c(sd(ys_p), sd(ys_e), sd(ys_r))
print(sds)
```
---
title: "Plot benchmark results"
author: "Wanda De Keersmaecker"
date: "10/1/2020"
output: html_document
vignette: >
   %\VignetteIndexEntry{Prototype of benchmark study}
   %\VignetteEngine{knitr::rmarkdown}
   %\usepackage[utf8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, 
                      eval = FALSE,
                      fig.width=18,
                      fig.height=10)
library(BenchmarkRecovery)
library(reshape2)
library(plyr)
```

```{r load-data}
# inputs
# input folder
ifolder <- '../data/'#'/home/wanda/Documents/data/benchmarkRecovery/Run_20201016_expRecPeriod/'#ifolder <- '/home/wanda/Documents/data/benchmarkRecovery/Run_20201005_4steps/'#'../data/'#ifolder <- '/home/wanda/Documents/data/benchmarkRecovery/Run_20201005_4steps/'
# Folder where outputs will be written
ofolder <- '../data/'#'/home/wanda/Documents/data/benchmarkRecovery/Run_20201016_expRecPeriod/Figures/'#ofolder <-'/home/wanda/Documents/data/benchmarkRecovery/Run_20201005_4steps/Figures'#'../data/' #ofolder <- '/home/wanda/Documents/data/benchmarkRecovery/Run_20201005_4steps/Figures'
# Name of the input dataset
basename <- 'LSTS_RndmSample_NoFire_5_Tree_80_scl_30_npnt_20000_VI' 
caseList <- c('seasAmp','remSd','distT','distRec','missVal','distMag')# 'distMag',evaluated time series characteristics for 

```

```{r general-settings,  echo = F, include=F}
simFullName  <- list('Disturbance magnitude',
                     'Disturbance timing',
                     'Recovery period',
                     'Seasonal amplitude',
                     'Noise level',
                     'Missing values')
names(simFullName) <- c('distMag',
                        'distT',
                        'distRec',
                        'seasAmp',
                        'remSd',
                        'missVal')
metric_list <- c('MAPE', 'R2', 'nTS', 'RMSE')
metric_names <- c('MAPE', 'R²', 'nTS', 'RMSE'); names(metric_names) <- c('MAPE', 'R2', 'nTS', 'RMSE')
tempRes_list <- c('quarterly', 'annual', 'dense', 'all')

```

```{r prepare-data,  echo = F, include=F}
dat_list <- list()
for(mm in 1:length(metric_list)){
  metric <- metric_list[mm]
  for(vr in 1:length(caseList)){
  evr <- caseList[vr]# name of parameter 
  
  RRI_dat <- loadRData(file.path(ifolder, paste0(basename, '_RRI_', metric, '_' , evr, '.rda')))
  R80p_dat <- loadRData(file.path(ifolder, paste0(basename, '_R80p_', metric, '_' , evr, '.rda')))
  YrYr_dat <- loadRData(file.path(ifolder, paste0(basename, '_YrYr_', metric, '_' , evr, '.rda')))
  
  tot_dat <- melt(rbind(RRI_dat, R80p_dat, YrYr_dat))
  
  tot_dat$Period <- revalue(factor(tot_dat$nPostMin), c("1"="Short", "4"="Long"))
  
   if((evr == 'remSd') || (evr == 'seasAmp') || (evr == 'missVal')) {
    tot_dat$variable <-mapvalues(tot_dat$variable, from = levels(tot_dat$variable), to = c("very low", "low", 'medium', 'high'))
    }  else{
    tot_dat$variable <-mapvalues(tot_dat$variable, levels(tot_dat$variable), to = c("very low", "low", 'medium', 'high'))
    }
  
  tot_dat$param = simFullName[[caseList[[vr]]]]
  if(vr == 1){totp_dat <- tot_dat}else{totp_dat <- rbind(totp_dat,tot_dat)}
  }
  totp_dat$param[totp_dat$param == 'SD remainder'] <- 'Noise level'
  totp_dat$paramType <- 'Environmental parameter'
  totp_dat[(totp_dat$param == 'Disturbance magnitude' | totp_dat$param == 'Recovery period' | totp_dat$param == 'Disturbance timing'), ]$paramType <- 'Disturbance parameter'
  totp_dat$param <- factor(totp_dat$param, levels = simFullName)
  
  dat_list[[metric]] <- totp_dat
}
```

Compare the performance of each recovery indicator
```{r plot-compare-indicators,  echo = F, include=F}
for(mm in 1:length(metric_list)){
  metric <- metric_list[mm]
  # plot 
  plt <- plotMet(dat_list[[metric]],  'Metric', metric_names[[metric]])
  png(file.path(ofolder, paste0(basename, '_', metric, 'Met.png')),width = 1311,height =628 )
  print(plt)
  dev.off()
  print(plt)
}


```

Which characteristics influence the performance the most?

```{r plot-sensitivity-overall, echo = F, include=F}
for(tr in 1:length(tempRes_list)){
  tempRes <- tempRes_list[tr]
  for(mm in 1:length(metric_list)){
    metric <- metric_list[mm]
    data <- dat_list[[metric]]
    if(tempRes != 'all'){data <- data[data$Dense == tempRes, ]}
    xlbl <- 'Parameter value'
    ylbl <- metric_names[[metric]]
    scales = 'free_y'
    # plot 
    plt <- plotSensBar(data, xlbl, ylbl, scales)
    png(file.path(ofolder, paste0(basename, '_', tempRes, '_',metric,'_Env.png')),width = 1911,height =1828 )
    print(plt)
    dev.off()
  }
  
}

```

```{r plot-sensitivity-separate, echo = F, include=F}

for(mm in 1: length(metric_list)){
  metric <- metric_list[mm]
  for(pp in 1:length(simFullName)){
    data <- dat_list[[metric]]
    data <- data[data$param == simFullName[pp],]
    data$Dense <- revalue(factor(data$Dense, levels = c('dense', 'quarterly', 'annual')), c("dense" = "no", "annual"="annual", "quarterly" = "quarterly"))
  data$Smooth <- revalue(factor(data$Smooth, levels = c('raw', 'smoothed', 'segmented')), c("raw"="no", "smoothed"="rolling mean", "segmented" = "segmentation"))
  xlbl <- simFullName[pp]
  ylbl <- metric_names[[metric]]
  scales = 'free_y'
  plt <- plotSens(data, xlbl, ylbl, scales)
  png(file.path(ofolder, paste0(basename, '_',metric,'_', names(simFullName[pp]),'_Sens.png')),width = 1911,height =1828 )
    print(plt)
    dev.off()
  }
}
```

How can we improve the performance?
```{r plot-preprocessing, echo = F, include=F}

for(mm in 1:length(metric_list)){
  metric <- metric_list[mm]
  data <- dat_list[[metric]]
  data$Dense <- revalue(factor(data$Dense, levels = c('dense', 'quarterly', 'annual')), c("dense" = "no", "annual"="annual", "quarterly" = "quarterly"))
  data$Smooth <- revalue(factor(data$Smooth, levels = c('raw', 'smoothed', 'segmented')), c("raw"="no", "smoothed"="rolling mean", "segmented" = "segmentation"))

  xlbl <- 'Temporal aggregation'
  ylbl <- metric_names[[metric]]
  scales = 'free_y'
  # plot 
  plt <-pltPrepBox(data, xlbl, ylbl, scales) 
  png(file.path(ofolder, paste0(basename, '_', metric,'_Prep.png')),width = 1911,height =1828 )
  print(plt)
  dev.off()
}

```

---
title: "Prototype of benchmark study"
author: "Wanda De Keersmaecker, Pablo Rodriguez-Sanchez"
date: ""
output: html_document
vignette: >
   %\VignetteIndexEntry{Prototype of benchmark study}
   %\VignetteEngine{knitr::rmarkdown}
   %\usepackage[utf8]{inputenc}
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, 
                      eval = T,
                      fig.width=18,
                      fig.height=10)

# Load libraries
library(devtools)
library(BenchmarkRecovery)
library(zoo)
library(plyr)
library(reshape2)
library(ggplot2)
library(profvis)
library(parallel)
library(pbapply)
library(bfast)
library(strucchange)
library(colorspace)

# Calculate the number of cores
no_cores <- detectCores()
```

# Extract time series characteristics from sampled time series
The following characteristics are derived from the decomposed time series:

* fraction of missing values
* seasonal amplitude (one value for input each time series)
* average seasonal pattern
* offset 
* standard deviation of the remainder component (one value per input time series)
* fitted ARMA model in the remainder component (one model per input time series)

```{r load-data}
# inputs
ifolder <- '../data/'
ofolder <- '../data/' # Folder where outputs will be written
basenames <- c('toyset', 'otherset') # Name of the input dataset (contains sampled satellite time series)
nyr <- 18 # Number of years in observation period
nobsYr <- 365 # Number of observations per year
caseList <- c('seasAmp','remSd','distMag','distT','distRec','missVal')# parameters to be evaluated
```

```{r aux-function, echo=FALSE}
# NOTE: consider refactoring this chunk
# For instance, moving extract_settings to ./R and/or creating a settings generator
# and working through load/save sttngs into/from a file

# TODO: document this function
extract_settings <- function(basename, ifolder, nyr, nobsYr, ofolder = '') {
  ifileVI <- paste0(basename, '.rda') # For instance, toyset.rda
  dfVi <- loadRData(file = file.path(ifolder, ifileVI)) # For instance /data/toyset.rda 
  
  # First decompose time series into seasonality, trend and remainder:
  tmp <- decompTSbfast(dfVi, nyr, nobsYr)
  dataVISeasbf <-  tmp[[1]] # Sesonality (fitted harmonic function)
  dataVIRembf <- tmp[[2]] # Remainder
  dataVITrbf <- tmp[[3]] # Trend (linear trend without break)
  dataVISeasCoef <- tmp[[4]] # Coefficients of fitted harmonic functions 
  
  # Derive characteristics
  tsVIMissVal <- rowSums(is.na(dfVi))/(dim(dfVi)[2]-2) # Fraction missing values
  
  seasVImax <- apply(dataVISeasbf[,-c(1,2)], 1, max) # Seasonal amplitude for each pixel
  seasS <- dataVISeasbf[dataVISeasbf[,1]<0,] # Average seasonal pattern, only southern hemisphere to avoid interference of seasonal cycles
  seasVImean <- colMeans(as.matrix(seasS[,-c(1,2)]))
  
  TrVImean <- mean(rowMeans(as.matrix(dataVITrbf[,-c(1,2)])), na.rm=T) # Offset 
  
  Rem_VIsd <- apply(dataVIRembf[,-c(1,2)], 1, sd, na.rm=T) # SD of remainder per pixel
 
    # Settings  simulation
  STnobsYr <- 365
  Vqntl <- c( .05, .25, .4, .6, .75, .95)#c( .05, .25, .5, .75, .95)#c( .5)# c( .05, .25, .5, .75)# set of quantiles used to derive realistic values (for number of 
  
  sttngs <- list()
  sttngs$'seasAmp' <- list(type = 'dist', vals = c(0,0,quantile(seasVImax, Vqntl)), obs = seasVImax, fix = quantile(seasVImax, c(.4, .6)))
  sttngs$'missVal' <- list(type = 'dist',  vals = c(1-1/16,1-1/16,quantile(tsVIMissVal, Vqntl)), obs = tsVIMissVal, fix = quantile(tsVIMissVal, c(.4, .6)))
  sttngs$'remSd' <- list(type = 'dist',  vals = c(0,0,quantile(Rem_VIsd, Vqntl)), obs = Rem_VIsd, fix = quantile(Rem_VIsd, c(.4, .6)))
  sttngs$'nyr' <- list(type = 'range', vals = c(25), fix = 25)#seq(6,36, by = 6)
  sttngs$'distMag' <- list(type = 'range', vals = -c(0.1,0.2,0.25,0.35,0.4,0.5), fix = -c(0.25,0.35))
  sttngs$'distT' <- list(type = 'range', vals = c(3,6,9,12,15,18), fix = c(9,12))#seq(3,33, by = 6)
  sttngs$'distRec' <- list(type = 'range', vals = c(0.5,2,2.25,3.75,4,6.5)*STnobsYr, fix = c(2.25,3.75)*STnobsYr) #seq(0.5,6.5,by=0.5)
  sttngs$'nDr' <- list(type = 'range',  vals = c(0), fix = 0)
  sttngs$'distType' <- list(type = 'cat',  vals = c('piecewise'), fix = 'piecewise')#piecewise exponential diffEq
  sttngs$'DistMissVal' <- list(type = 'cat', vals = 'random', fix = 'random')
  sttngs$'trAv' <- list(type = 'range', vals = TrVImean, fix = TrVImean)
  sttngs$'general' <- list(
    eval = caseList,#parameters to be evaluated, can be  'distT','distRec', 'missVal' 'distMag', 'seasAmp', 'remSd'
    nTS = 100,
    nobsYr = STnobsYr,
    seasAv = seasVImean,
    # remcoef = Rem_VIcoef,
    parSetUp = 'int') # Parameter set-up: can be avg dist, comb, or int
  
  # remove redundant variables
  rm(list=setdiff(ls(), c("sttngs",'ifolder', 'ofolder', 'basename', 'ncores', 'pars')))
  
  # Save if desired
  save_results = (ofolder != '')
  if(save_results) {
      save(sttngs, file = file.path(ofolder, paste0(basename,  '_simTS_settings.rda')))
  }
  
  return(sttngs)
}

```

```{r extract-settings}
# Start clock
start_time <- Sys.time()

# Extract settings from data
sttngs_list <- mclapply(basenames, FUN = extract_settings, ifolder = ifolder, nyr = nyr, nobsYr = nobsYr, ofolder = ofolder, mc.cores = no_cores - 1)

# Stop clock
end_time <- Sys.time()
print(end_time - start_time)
```

# Simulate time series, measure recovery and evaluate performance
Based on the characteristics of the measured time series, time series are simulated.

## Define simulation settings 
First, the simulation settings are defined in a settings list: 

* __seasAmp__: seasonal amplitude
* __missVal__: fraction of missing values in time series
* __remSd__: standard deviation of the remainder
* __nyr__: time series length (number of years)
* __distMag__: disturbance magnitude, should be a negative value to simulate a drop
* __distT__: timing disturbance (disturbance year)
* __distRec__: recovery half time after disturbance (number of observations)
* __nDr__: number of simulated droughts 
* __distType__: Type of recovery (piecewise, exponential or realistic)
* __DistMissVal__: defines how the introduced missing values should be distributed: random or at an equal interval
* __trAv__: offset
* __eval__: the parameters that should be evaluated
* __nTS__: number of time series simulated per value of the evaluated parameter
* __nobsYr__: number of observations simulated per year (this should equal the frequency of the sampled time series)
* __seasAv__: represents the seasonal pattern
* __parSetUp__: the parameter values selection approach (avg, dist, comb, or int)

For each of these parameters, the type is defined. There are three main parameter types: *dist*, *range*, and *cat*. For the *dist* parameters, parameter values that are observed from sampled time series are available. For *range* parameters no observed values are available, but an expected range of their values can be set. Finally, the *cat* parameters are categoric.

Next to the type of the parameter, a set of parameter values (*vals*) need to be defined. If the sensitivity of the recovery indicators to the parameter is being evaluated, the performance of the recovery indicators is evaluated with respect to each of these parameter values. For *dist* parameters, the observed parameter values (*obs*) need to be additionally provided. If the parameter set-up follows the *int* approach (see next section for more details), the values for the evaluated *dist* parameter are not set to predefined fixed values, but are randomly sampled over an interval. These intervals are defined in the *vals* setting: the first two values refer to the first interval, the next two values to the second interval, and so on.

While evaluating the sensitivity to one specific parameter, the values of all other parameters also need to be set. Four approaches can be used to achieve this. First, the *avg* and *int* approaches set the parameters to their average value. This equals the mean value of the distribution of observed values of *dist* parameters, the mean value of the range of values (given by *vals*) for *range* parameters and a randomly selected value (selected from the values give by *vals*) for *cat* parameters. Second, the *dist* approach selects values for the parameters given by the likelihood of their occurrence. The likelihood is defined by the histogram of observed values for *dist* parameters and a random selection of values is made for *range* and *cat* parameters. Third, the *comb* approach defines all combinations of evaluated values for each parameter (as defined in *vals*). 

The following inputs are needed to calculate the recovery indicators:

__funSet__: list of settings for the computation of the recovery indicators. More than one value for each setting is allowed (yet an equal number of values for each parameter is required). The recovery indicators are then derived for each set of values of the setting parameters.

+ *freq*:  'dense', 'annual', or 'quarterly'. Defines the observation frequency. For 'dense' the original frequency is used. For 'annual' and 'quarterly', the time series are converted to annual or quarterly frequency, respectively.
+ *input*: 'smoothed', 'raw', 'segmented'. Defines the type of time series that is used for the recovery indicators. For 'raw', the simulated time series are directly used to calculate recovery, for 'smooth' a time series smoothing algorithm (rolling mean) is used before recovery calculation, for 'segmented' the trend component of the piecewise regression (BFAST0n) is used.
+ *nPre*: the number of years before the disturbance used to derive the pre-disturbance values
+ *nDist*: the number of years after the disturbance used to derive the value during the disturbance 
+ *nPostMin* and *nPostMax*: the post-disturbance values are derived between nPostMin and nPostMax years after the disturbance
+ *h*: in case *input* equals 'segmented', the *h* value is used in the segmentation algorithm to define the minimal segment size either given as fraction relative to the sample size or as an integer giving the minimal number of observations in each segment
+ *breaks*: in case *input* equals 'segmented', the criterium given by *breaks* is used in the segmentation algorithm to define the optimal number of segments. Can be set to 'BIC' or 'LWZ' (but has been deactivated)
+ *seas*: in case *input* equals 'segmented', *seas* denotes whether a seasonal term needs to be used in the segmentation algorithm 


```{r set-funSet}
# recovery settings
funSet <- list('freq' = c(rep('annual', 6),  rep('dense',6), rep('quarterly',6)),
               'input' = rep(c('raw', 'smoothed','segmented'),6),# settings for the recovery indicators
               'nPre' = rep(2,18),
               'nDist' = c(rep(0,6),rep(1,12)),
               'nPostMin' = rep(c(4,4,4,1,1,1), 3),
               'nPostMax' = c(rep(5,3), rep(1,3), rep(6,3), rep(2,3), rep(6,3), rep(2,3)),
               'h' = rep(0.15,18),
               'seas' = c(rep(F,6),rep(T,12)))

# Save all configurations (one per basename)
# although all of them are now identical
# NOTE: decide if we need this
mclapply(basenames, FUN = function(basename) { save(funSet, file = file.path(ofolder, paste0(basename, '_recSettings.rda')))  })


```

The specified settings are then used to simulate time series, calculate recovery indicators and evaluate their performance:
```{r set-input-table}
# Create all the tests to be performed and store them on a table
# This is the 'homework' list for the HPC
sttngs$general$eval
testsTable <- expand.grid(basename = basenames, case = caseList, stringsAsFactors = FALSE) # All combinations of cases and filenames

# Assign a column with the settings
# NOTE: in the future, perhaps use only basename (as sttngs is redundant)
testsTable$settings <- rep(NA, nrow(testsTable)) # Initialize settings column
testsTable$settings <- sttngs_list
testsTable <- tibble::tibble(testsTable) # Tibblify for increased readability

# Print for inspection
print(testsTable)
```

```{r simulate-and-extract, eval=TRUE}
# Start clock
start_time <- Sys.time()
set_fast_options()

# Run the sensitivity analysis
results_list <- mcmapply(FUN = evalParam, evr = testsTable$case, basename = testsTable$basename, sttngs = testsTable$settings, # Iterate along these parameters
                         MoreArgs = list(funSet = funSet, ofolder = ofolder), # These parameters remain constant
                         mc.cores = no_cores - 1)

# stop clock
end_time <- Sys.time()
print(end_time - start_time)

```

# Plot the performance indicators

```{r plot2, include = F}
# general settings
# characteristics for which a plot needs to be made


simFullName  <- list('Disturbance magnitude',
                     'Number of droughts',
                     'Time series length',
                     'Seasonal amplitude',
                     'SD remainder',
                     'Disturbance timing',
                     'Recovery period [years]',
                     'Missing values')
names(simFullName) <- c('distMag',
                        'nDr', 
                        'len',
                        'seasAmp',
                        'remSd',
                        'distT',
                        'distRec',
                        'missVal')
```

Compare the performance of each recovery indicator
```{r compare-indicators, include=F}
for(vr in 1:length(caseList)){
  evr <- caseList[vr]# name of parameter that will be evaluated in the simulation
  # setvr <- sttngs[[evr]]# settings of simulation
  
  # NOTE: retrieve data from results_list instead of from saved files
  RRI_rsq <- loadRData(file.path(ofolder, paste0(basename, '_RRI_R2_' , evr, '.rda')))
  R80p_rsq <- loadRData(file.path(ofolder, paste0(basename, '_R80p_R2_' , evr, '.rda')))
  YrYr_rsq <- loadRData(file.path(ofolder, paste0(basename, '_YrYr_R2_' , evr, '.rda')))
  
  RRI_mape <- loadRData(file.path(ofolder, paste0(basename, '_RRI_MAPE_' , evr, '.rda')))
  R80p_mape <- loadRData(file.path(ofolder, paste0(basename, '_R80p_MAPE_' , evr, '.rda')))
  YrYr_mape <- loadRData(file.path(ofolder, paste0(basename, '_YrYr_MAPE_' , evr, '.rda')))
  
  RRI_nTS <- loadRData(file.path(ofolder, paste0(basename, '_RRI_nTS_' , evr, '.rda')))
  R80p_nTS <- loadRData(file.path(ofolder, paste0(basename, '_R80p_nTS_' , evr, '.rda')))
  YrYr_nTS <- loadRData(file.path(ofolder, paste0(basename, '_YrYr_nTS_' , evr, '.rda')))
  
  tot_rsq <- melt(rbind(RRI_rsq, R80p_rsq, YrYr_rsq))
  tot_mape <- melt(rbind(RRI_mape, R80p_mape, YrYr_mape))
  tot_nTS <- melt(rbind(RRI_nTS, R80p_nTS, YrYr_nTS))
  
  
  tot_rsq$Period <- revalue(factor(tot_rsq$nPostMin), c("1"="Short", "4"="Long"))
  tot_mape$Period <- revalue(factor(tot_mape$nPostMin), c("1"="Short", "4"="Long"))
  tot_nTS$Period <- revalue(factor(tot_nTS$nPostMin), c("1"="Short", "4"="Long"))
  
   if((evr == 'remSd') || (evr == 'seasAmp') || (evr == 'missVal')) {
    tot_rsq$variable <- mapvalues(tot_rsq$variable, from = levels(tot_rsq$variable), to = c("no", "low", 'medium', 'high'))
    tot_mape$variable <-mapvalues(tot_mape$variable, from = levels(tot_mape$variable), to = c("no", "low", 'medium', 'high'))
    tot_nTS$variable <-mapvalues(tot_nTS$variable, from = levels(tot_nTS$variable), to = c("no", "low", 'medium', 'high'))
  }  else{
    tot_rsq$variable <-mapvalues(tot_rsq$variable, levels(tot_rsq$variable), to = c("low", 'medium', 'high'))
    tot_mape$variable <-mapvalues(tot_mape$variable, levels(tot_mape$variable), to = c("low", 'medium', 'high'))
    tot_nTS$variable <-mapvalues(tot_nTS$variable, levels(tot_nTS$variable), to = c("low", 'medium', 'high'))
  }
  # tot_rsq <- tot_rsq[(tot_rsq$Dense == 'dense' & tot_rsq$Smooth == 'raw' & tot_rsq$Period == 'Long'),]
  tot_rsq$param = simFullName[[caseList[[vr]]]]
  if(vr == 1){totp_rsq <- tot_rsq}else{totp_rsq <- rbind(totp_rsq,tot_rsq)}
  
  # tot_mape <- tot_mape[(tot_mape$Dense == 'dense' & tot_mape$Smooth == 'raw' & tot_mape$Period == 'Long'),]
  tot_mape$param = simFullName[[caseList[[vr]]]]
  if(vr == 1){totp_mape <- tot_mape}else{totp_mape <- rbind(totp_mape,tot_mape)}
  
  # tot_nTS <- tot_nTS[(tot_nTS$Dense == 'dense' & tot_nTS$Smooth == 'raw' & tot_nTS$Period == 'Long'),]
  tot_nTS$param = simFullName[[caseList[[vr]]]]
  if(vr == 1){totp_nTS <- tot_nTS}else{totp_nTS <- rbind(totp_nTS,tot_nTS)}
}

xlbl <- 'Metric'
  
  # plot R2
  data <- totp_rsq
  ylbl <- 'R²'
  pltR2 <- plotMet(data, xlbl, ylbl)
  png(file.path(ofolder, paste0(basename, '_RsqMet.png')),width = 1311,height =628 )
  print(pltR2)
  dev.off()
  
  # plot MAPE
  data <- totp_mape
  ylbl <- 'MAPE'
  pltMAPE <- plotMet(data,  xlbl, ylbl)
  png(file.path(ofolder, paste0(basename, '_MAPEMet.png')),width = 1311,height =628 )
  print(pltMAPE)
  dev.off()
   # plot fraction of time series processed
  data <- totp_nTS
  ylbl <- 'Fraction'
  pltnTS <- plotMet(data,  xlbl, ylbl)
  png(file.path(ofolder, paste0(basename, '_nTSMet.png')),width = 1311,height =628 )
  print(pltnTS)
  dev.off()
  
  print(pltR2)
  print(pltMAPE)
  print(pltnTS)

```

Which characteristics influence the performance the most?

```{r compare-performance, include=F}
# compare effect of each parameter on the 
caseList <- c( 'seasAmp', 'remSd', 'missVal','distRec','distT','distMag')# evaluated time series characteristics for which a plot needs to be made

simFullName  <- list('Disturbance magnitude',
                     'Disturbance timing',
                     'Recovery period',
                     'Number of droughts',
                     'Time series length',
                     'Seasonal amplitude',
                     'SD remainder',
                     'Missing values')
names(simFullName) <- c('distMag',
                        'distT',
                        'distRec',
                        'nDr', 
                        'len',
                        'seasAmp',
                        'remSd',
                        'missVal')

for(vr in 1:length(caseList)){
  evr <- caseList[vr]# name of parameter that will be evaluated in the simulation
  # setvr <- sttngs[[evr]]# settings of simulation
  
  RRI_rsq <- loadRData(file.path(ofolder, paste0(basename, '_RRI_R2_' , evr, '.rda')))
  R80p_rsq <- loadRData(file.path(ofolder, paste0(basename, '_R80p_R2_' , evr, '.rda')))
  YrYr_rsq <- loadRData(file.path(ofolder, paste0(basename, '_YrYr_R2_' , evr, '.rda')))
  
  RRI_mape <- loadRData(file.path(ofolder, paste0(basename, '_RRI_MAPE_' , evr, '.rda')))
  R80p_mape <- loadRData(file.path(ofolder, paste0(basename, '_R80p_MAPE_' , evr, '.rda')))
  YrYr_mape <- loadRData(file.path(ofolder, paste0(basename, '_YrYr_MAPE_' , evr, '.rda')))
  
  RRI_nTS <- loadRData(file.path(ofolder, paste0(basename, '_RRI_nTS_' , evr, '.rda')))
  R80p_nTS <- loadRData(file.path(ofolder, paste0(basename, '_R80p_nTS_' , evr, '.rda')))
  YrYr_nTS <- loadRData(file.path(ofolder, paste0(basename, '_YrYr_nTS_' , evr, '.rda')))
  
  tot_rsq <- melt(rbind(RRI_rsq, R80p_rsq, YrYr_rsq))
  tot_mape <- melt(rbind(RRI_mape, R80p_mape, YrYr_mape))
  tot_nTS <- melt(rbind(RRI_nTS, R80p_nTS, YrYr_nTS))
  
  tot_rsq$Period <- revalue(factor(tot_rsq$nPostMin), c("1"="Short", "4"="Long"))
  tot_mape$Period <- revalue(factor(tot_mape$nPostMin), c("1"="Short", "4"="Long"))
  tot_nTS$Period <- revalue(factor(tot_nTS$nPostMin), c("1"="Short", "4"="Long"))
  
   if((evr == 'remSd') || (evr == 'seasAmp') || (evr == 'missVal')) {
    tot_rsq$variable <- mapvalues(tot_rsq$variable, from = levels(tot_rsq$variable), to = c("no", "low", 'medium', 'high'))
    tot_mape$variable <-mapvalues(tot_mape$variable, from = levels(tot_mape$variable), to = c("no", "low", 'medium', 'high'))
    tot_nTS$variable <-mapvalues(tot_nTS$variable, from = levels(tot_nTS$variable), to = c("no", "low", 'medium', 'high'))
  }  else{
    tot_rsq$variable <-mapvalues(tot_rsq$variable, levels(tot_rsq$variable), to = c("low", 'medium', 'high'))
    tot_mape$variable <-mapvalues(tot_mape$variable, levels(tot_mape$variable), to = c("low", 'medium', 'high'))
    tot_nTS$variable <-mapvalues(tot_nTS$variable, levels(tot_nTS$variable), to = c("low", 'medium', 'high'))
  }
  tot_rsq <- tot_rsq[(tot_rsq$Dense == 'dense' & tot_rsq$Smooth == 'raw' & tot_rsq$Period == 'Long'),]
  tot_rsq$param = simFullName[[caseList[[vr]]]]
  if(vr == 1){totp_rsq <- tot_rsq}else{totp_rsq <- rbind(totp_rsq,tot_rsq)}
  
  tot_mape <- tot_mape[(tot_mape$Dense == 'dense' & tot_mape$Smooth == 'raw' & tot_mape$Period == 'Long'),]
  tot_mape$param = simFullName[[caseList[[vr]]]]
  if(vr == 1){totp_mape <- tot_mape}else{totp_mape <- rbind(totp_mape,tot_mape)}
  
  tot_nTS <- tot_nTS[(tot_nTS$Dense == 'dense' & tot_nTS$Smooth == 'raw' & tot_nTS$Period == 'Long'),]
  tot_nTS$param = simFullName[[caseList[[vr]]]]
  if(vr == 1){totp_nTS <- tot_nTS}else{totp_nTS <- rbind(totp_nTS,tot_nTS)}
}


  totp_rsq$paramType <- 'Environmental parameter'
  totp_rsq[(totp_rsq$param == 'Disturbance magnitude' | totp_rsq$param == 'Recovery period' | totp_rsq$param == 'Disturbance timing'), ]$paramType <- 'Disturbance parameter'
  totp_mape$paramType <- 'Environmental parameter'
  totp_mape[(totp_mape$param == 'Disturbance magnitude' | totp_mape$param == 'Recovery period' | totp_mape$param == 'Disturbance timing' ),]$paramType  <- 'Disturbance parameter'
  totp_nTS$paramType <- 'Environmental parameter'
  totp_nTS[(totp_rsq$param == 'Disturbance magnitude' | totp_nTS$param == 'Recovery period' | totp_nTS$param == 'Disturbance timing' ),]$paramType  <- 'Disturbance parameter'

data <- totp_rsq
  data$param <- factor(data$param, levels = rev(unlist(simFullName[caseList])))
  xlbl <- 'Parameter value'
  ylbl <- 'R²'
  pltR2 <- plotEnv(data, xlbl, ylbl, scales = 'free_y')
  png(file.path(ofolder, paste0(basename, '_Rsq_Env.png')),width = 1311,height =628 )
  print(pltR2)
  dev.off()
  
  data <- totp_mape
  data$param <- factor(data$param, levels = rev(unlist(simFullName[caseList])))
  data$value[is.infinite(data$value)] <- NA
  xlbl <- 'Parameter value'
  ylbl <- 'MAPE'
  pltMAPE <- plotEnv(data, xlbl, ylbl, scales = 'free_y')
  png(file.path(ofolder, paste0(basename, '_MAPE_Env.png')),width = 1311,height =628 )
  print(pltMAPE)
  dev.off()
  
  data <- totp_nTS
  data$param <- factor(data$param, levels = rev(unlist(simFullName[caseList])))
  xlbl <- 'Parameter value'
  ylbl <- 'Fraction'
  pltnTS <- plotEnv(data, xlbl, ylbl, scales = 'free_y')
  png(file.path(ofolder, paste0(basename, '_nTS_Env.png')),width = 1311,height =628 )
  print(pltnTS)
  dev.off()
  
  print(pltR2)
  print(pltMAPE)
  print(pltnTS)
```

How can we improve the performance?

```{r improve-performance, include=F}
for(vr in 1:length(caseList)){
  evr <- caseList[vr]# name of parameter that will be evaluated in the simulation
  # setvr <- sttngs[[evr]]# settings of simulation
  
  RRI_rsq <- loadRData(file.path(ofolder, paste0(basename, '_RRI_R2_' , evr, '.rda')))
  R80p_rsq <- loadRData(file.path(ofolder, paste0(basename, '_R80p_R2_' , evr, '.rda')))
  YrYr_rsq <- loadRData(file.path(ofolder, paste0(basename, '_YrYr_R2_' , evr, '.rda')))
  
  RRI_rmse <- loadRData(file.path(ofolder, paste0(basename, '_RRI_RMSE_' , evr, '.rda')))
  R80p_rmse <- loadRData(file.path(ofolder, paste0(basename, '_R80p_RMSE_' , evr, '.rda')))
  YrYr_rmse <- loadRData(file.path(ofolder, paste0(basename, '_YrYr_RMSE_' , evr, '.rda')))
  
  RRI_mape <- loadRData(file.path(ofolder, paste0(basename, '_RRI_MAPE_' , evr, '.rda')))
  R80p_mape <- loadRData(file.path(ofolder, paste0(basename, '_R80p_MAPE_' , evr, '.rda')))
  YrYr_mape <- loadRData(file.path(ofolder, paste0(basename, '_YrYr_MAPE_' , evr, '.rda')))
  
  RRI_nTS <- loadRData(file.path(ofolder, paste0(basename, '_RRI_nTS_' , evr, '.rda')))
  R80p_nTS <- loadRData(file.path(ofolder, paste0(basename, '_R80p_nTS_' , evr, '.rda')))
  YrYr_nTS <- loadRData(file.path(ofolder, paste0(basename, '_YrYr_nTS_' , evr, '.rda')))
  
  tot_rsq <- melt(rbind(RRI_rsq, R80p_rsq, YrYr_rsq))
  tot_rmse <- melt(rbind(RRI_rmse, R80p_rmse, YrYr_rmse))
  tot_mape <- melt(rbind(RRI_mape, R80p_mape, YrYr_mape))
  tot_nTS <- melt(rbind(RRI_nTS, R80p_nTS, YrYr_nTS))
  
  tot_rsq$Period <- revalue(factor(tot_rsq$nPostMin), c("1"="Short", "4"="Long"))
  tot_rmse$Period <- revalue(factor(tot_rmse$nPostMin), c("1"="Short", "4"="Long"))
  tot_mape$Period <- revalue(factor(tot_mape$nPostMin), c("1"="Short", "4"="Long"))
  tot_nTS$Period <- revalue(factor(tot_nTS$nPostMin), c("1"="Short", "4"="Long"))
  
   if((evr == 'remSd') || (evr == 'seasAmp') || (evr == 'missVal')) {
    tot_rsq$variable <- mapvalues(tot_rsq$variable, from = levels(tot_rsq$variable), to = c("no", "low", 'medium', 'high'))
    tot_rmse$variable <-mapvalues(tot_rmse$variable, from = levels(tot_rmse$variable), to = c("no", "low", 'medium', 'high'))
    tot_mape$variable <-mapvalues(tot_mape$variable, from = levels(tot_mape$variable), to = c("no", "low", 'medium', 'high'))
    tot_nTS$variable <-mapvalues(tot_nTS$variable, from = levels(tot_nTS$variable), to = c("no", "low", 'medium', 'high'))
  } else{ 
    tot_rsq$variable <- mapvalues(tot_rsq$variable, from = levels(tot_rsq$variable), to = c( "low", 'medium', 'high'))
    tot_rmse$variable <-mapvalues(tot_rmse$variable, from = levels(tot_rmse$variable), to = c( "low", 'medium', 'high'))
    tot_mape$variable <-mapvalues(tot_mape$variable, from = levels(tot_mape$variable), to = c( "low", 'medium', 'high'))
    tot_nTS$variable <-mapvalues(tot_nTS$variable, from = levels(tot_nTS$variable), to = c( "low", 'medium', 'high'))
    }
  
  lbls <- c("raw, BAP", 'piecewise, BAP', 'smoothed, BAP', "raw, dense", 'piecewise, dense', 'smoothed, dense', "raw, quarterly", 'piecewise, quarterly', 'smoothed, quarterly')
  sname <- caseList[[vr]]
  xlbl <- simFullName[[sname]]
  
  # plot R2
  data <- tot_rsq
  ylbl <- 'R²'
  pltR2 <- plotSens(data, lbls, xlbl, ylbl, scales = 'fixed')
  png(file.path(ofolder, paste0(basename, '_Rsq_',evr,'.png')),width = 1311,height =628 )
  print(pltR2)
  dev.off()
  
  # plot RMSE
  data <- tot_rmse
  ylbl <- 'RMSE'
  pltRMSE <- plotSens(data, lbls, xlbl, ylbl, scales = 'free_y')
  png(file.path(ofolder, paste0(basename, '_RMSE_',evr,'.png')),width = 1311,height =628 )
  print(pltRMSE)
  dev.off()
   # plot MAPE
  data <- tot_mape
  ylbl <- 'MAPE'
  pltMAPE <- plotSens(data, lbls, xlbl, ylbl, scales = 'free_y')
  png(file.path(ofolder, paste0(basename, '_MAPE_',evr,'.png')),width = 1311,height =628 )
  print(pltMAPE)
  dev.off()
   # plot fraction of time series processed
  data <- tot_nTS
  ylbl <- 'Fraction'
  pltnTS <- plotSens(data, lbls, xlbl, ylbl, scales = 'free_y')
  png(file.path(ofolder, paste0(basename, '_nTS_',evr,'.png')),width = 1311,height =628 )
  print(pltnTS)
  dev.off()
  
  print(pltR2)
  print(pltRMSE)
  print(pltMAPE)
  print(pltnTS)
}
```



---
title: "Avoiding nested loops with apply functions"
author: "Pablo Rodriguez-Sanchez"
date: "9/3/2020"
output: html_document
vignette: >
   %\VignetteIndexEntry{Using apply functions}
   %\VignetteEngine{knitr::rmarkdown}
   %\usepackage[utf8]{inputenc}
---

## Problem description
Currently, our vignette `sensitivityAnalysis.Rmd` is using `evalParam` for creating outputs (in the form of saved files) corresponding to each one of the 96 permutations of the following parameters:

- `par1`, whose possible values are `(distMag, distRec, distT, missVal, remSd, seasAmp)`
- `par2`, whose possible values are `(R80p, RRI, simTS, YrYr)`
- `par3`, whose possible values are `(MAPE, nTS, R2, RMSE)`

The loop corresponding to `par1` happens at vignette level via an `lapply` function. The other two happen inside `evalParam` via nested `for` loops.

For the sake of code clarity, speed and parallelization, it could be a good idea to refactor `evalParam` in a more atomic way. By atomic I mean accepting one or more parameters identifying which one of the 96 possible cases we want to calculate, and calculating that and only that one. Depending on the amount of identification parameters, the whole process can be efficiently looped via vectorization, `lapply` (for a single id per row) or `mapply` (for multiple ids per row).

The `evalParam` function is already quite complicated, so we'll assume vectorization is not feasible. The two remaining possibilities are thus:

1. Use a single identifier (i.e.: `evalParam(..., id = "R80p_nTS_distT")`, to later be called via `lapply`).
2. Use a multiple identifier (i.e: `evalParam(..., par1 = "R80p", par2 = "nTS", par3 = "distT")`, to later be called via `mapply`).

###  Graphical summary
Currently `evalParam` is doing too much:

![](../img/diagram_now.png)

It will be advisable to split like:

![](../img/diagram_future.png)

In the present vignette I show a simple case study of using apply functions in this kind of problems.

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Generate the input data frame
We will generate a simple data set with all the 12 permutations corresponding to these identifiers:

- `type`, whose possible values are `("A", "B", "C")`.
- `gender`, whose possible values are `(1, 2)`.
- `country`, whose possible values are `("NL", "BE")`.

plus a single column containing some measurement (just a random number in this case).

```{r cases-definition}
# Identifiers
v1 <- c("A", "B", "C")
v2 <- c(1, 2)
v3 <- c("NL", "BE")

# Measurements
v4 <- runif(n = length(v1) * length(v2) * length(v3))
```

```{r permute-all-cases}
idf <- expand.grid(type = v1, gender = v2, country = v3, stringsAsFactors = FALSE)
idf <- cbind(idf, measurement = v4)
```

```{r see-dataset-nonames, echo=FALSE}
print(idf)
```

It is usually a good idea to assign meaningful names to the rows.

```{r create-single-ids}
# Auxiliary function.
# Creates a single row identifier by combining type, gender and country
# For instance: "B_2_NL"
create_id <- function(type, gender, country) {
  id = paste(type, gender, country, sep = "_")
}

rownames(idf) <- create_id(idf$type, idf$gender, idf$country) # Rownames checks that the names are not duplicated
```

```{r see-dataset, echo=FALSE}
print(idf)
```

If the name is well chosen (and ours is) it can even be redundant with the other three id columns. But for this tutorial we'll keep all of them, in order to investigate different ways of applying functions to a given row.

## Analyze the data
We want to perform the following analysis: 

- for each row
  - multiply the measurement by 2 if the country is Belgium.
  - multiply the measurement by -2 if the country is the Netherlands. 
  
The three functions below do the same, and only differ in the way the input row is specified:

```{r functions}
# Analyze (brute-force)
# This method is expected to be called row by row (for instance, inside a loop). It will crash otherwise
analyze <- function(row) {
  
  # Auxiliary function. 
  # Returns 2 for Belgium and -2 for The Netherlands
  country_to_number <- function(country) {
    if(country == "BE") {
      return(2)
    } else {
      return(-2)
    }
  }
  
  number <- country_to_number(row$country)
  return(number * row$measurement)
}

# Analyze (prepared for lapply)
# Same as analyze, but a single id (the row name) has to be provided
lanalyze <- function(id, data) {
  row <- data[id, ]
  analyze(row)
}

# Analyze (prepared for mapply)
# Same as analyze, but type, gender and country have to be provided
manalyze <- function(type, gender, country, data) {
  id <- create_id(type, gender, country)
  lanalyze(id, data)
}
```

The analysis itself is deliberately silly. It is just an example of an action to be:

1. Performed on an input data set.
2. Controlled by an input data set (in this case, the same one).

### Apply to a desired subset
```{r apply-sub}
lresults <- lapply(c("A_1_NL", "B_2_BE"), lanalyze, data = idf) # Using lapply (single row identifier)
mresults <- mapply(manalyze, c("A", "B"), c(1, 2), c("NL", "BE"), MoreArgs = list(data = idf)) # Using mapply (multiple row identifiers)
```

```{r print-results-sub}
print(lresults)
print(mresults)
```

### Apply to the whole dataset

```{r apply}
lresults <- lapply(create_id(idf$type, idf$gender, idf$country), lanalyze, data = idf) # Using lapply (single row identifier)
mresults <- mapply(manalyze, idf$type, idf$gender, idf$country, MoreArgs = list(data = idf)) # Using mapply (multiple row identifiers)
```

```{r print-results}
print(lresults)
print(mresults)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_eval_performance.R
\name{rsq}
\alias{rsq}
\title{R squared}
\usage{
rsq(x, y)
}
\arguments{
\item{x}{vector of x values}

\item{y}{vector of y values}
}
\value{
the R squared between the x and y variables
}
\description{
R squared
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{toyset}
\alias{toyset}
\title{Sample of 100 Landat NBR time series over tropical forest}
\format{
A list with 100 observations and 6576 variables (longitude(lon), latitude(lat) and NBR value for a set of dates):
\describe{
  \item{lon}{longitude}
  \item{lat}{latitude}
  ...
}
}
\source{
\url{https://www.usgs.gov/land-resources/nli/landsat}
}
\usage{
toyset
}
\description{
A dataset containing 100 sampled NBR time series over tropical forest
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_general.R
\name{TSdecompress}
\alias{TSdecompress}
\title{The TSdecompress recovers a vector that has been compressed by the TScompress function in its original format.}
\usage{
TSdecompress(ts)
}
\arguments{
\item{ts}{a vector compressed be the TScompress function}
}
\value{
the vector restored in its original format
}
\description{
The TSdecompress recovers a vector that has been compressed by the TScompress function in its original format.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_general.R
\name{loadRData}
\alias{loadRData}
\title{Load RData file and returns it. This function is a substitute for the load function, allowing to assign a user defined variable name when loading a RData file.}
\usage{
loadRData(fileName)
}
\arguments{
\item{fileName}{the path to the Rdata file that needs to be loaded}
}
\value{
R object
}
\description{
Load RData file and returns it. This function is a substitute for the load function, allowing to assign a user defined variable name when loading a RData file.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_characterize.R
\name{getARMAcoef}
\alias{getARMAcoef}
\title{Fit ARMA model: This function automatically fits an ARMA model without seasonal component.}
\usage{
getARMAcoef(tsx)
}
\arguments{
\item{tsx}{a time series object for which an ARMA model needs to be fitted.}
}
\value{
ARMA model coefficients. A list containing the ARMA coefficients and their order (number of AR coefficients, non-seasonal differences and MA coefficients
}
\description{
Fit ARMA model: This function automatically fits an ARMA model without seasonal component.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_general.R
\name{TScompress}
\alias{TScompress}
\title{The TScompress function compresses a time series vector. The purpose is to avoid storing many NA values.}
\usage{
TScompress(ts)
}
\arguments{
\item{ts}{a vector that will be compressed}
}
\value{
a vector containing, in order, the length of the input vector, the number of observations without NA, the observations without NA, and the values of the observations that are no NA
}
\description{
The TScompress function compresses a time series vector. The purpose is to avoid storing many NA values.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_eval_performance.R
\name{rmse}
\alias{rmse}
\title{#' Intercept, slope of the linear model fit and Shapiro-Wilk normality test of the residuals
#'
#' @param x vector of x values
#' @param y vector of y values
#'
#' @return the intercept and slope of the linear fit and p value of Shapiro-Wilk normality test
#' @export
linFit <- function(x, y) {
  ind <- which(is.infinite(x) | is.infinite(y))
  if (length(ind)>0){
    x <- x[-ind]
    y <- y[-ind]
  }
  if((sum(is.na(x)==F) > 3) &&(sum(is.na(y)==F)>3) && (length(unique(x[is.na(x) ==F]))>1)){
    mod <- lm(y~x)
    coef <- summary(mod)$coefficients
    tst <- shapiro.test(mod$residuals)#p-value > 0.05 implies that the distribution of the data is not significantly different from normal distribution
    out <- c(coef[1,1], coef[2,1], tst$p.value)# intercept and slope
  }else(out <- c(NA,NA,NA))
  out
}
RMSE}
\usage{
rmse(val, meas)
}
\arguments{
\item{val}{vector of x values}

\item{meas}{vector of y values}
}
\value{
the RMSE of the two vectors
}
\description{
#' Intercept, slope of the linear model fit and Shapiro-Wilk normality test of the residuals
#'
#' @param x vector of x values
#' @param y vector of y values
#'
#' @return the intercept and slope of the linear fit and p value of Shapiro-Wilk normality test
#' @export
linFit <- function(x, y) {
  ind <- which(is.infinite(x) | is.infinite(y))
  if (length(ind)>0){
    x <- x[-ind]
    y <- y[-ind]
  }
  if((sum(is.na(x)==F) > 3) &&(sum(is.na(y)==F)>3) && (length(unique(x[is.na(x) ==F]))>1)){
    mod <- lm(y~x)
    coef <- summary(mod)$coefficients
    tst <- shapiro.test(mod$residuals)#p-value > 0.05 implies that the distribution of the data is not significantly different from normal distribution
    out <- c(coef[1,1], coef[2,1], tst$p.value)# intercept and slope
  }else(out <- c(NA,NA,NA))
  out
}
RMSE
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_general.R
\name{toAnnualTS}
\alias{toAnnualTS}
\title{Convert time series to annual frequency: The toAnnualTS function converts a time series with n observations per year to an annual time series (one observation per year). The main concept is to select observations per year closest to a given day of year that have no missing value (NA). Here, the day of year for which the seasonality is maximum is being used.}
\usage{
toAnnualTS(tsseas, tsi, obspyr, dtmax = 1/12)
}
\arguments{
\item{tsseas}{vector of observations (time series) representing the seasonal component of the time series to be converted}

\item{tsi}{vector of observations (time series) that needs to be converted to an annual time series}

\item{obspyr}{number of observations per year of the time series to be converted}

\item{dtmax}{maximum time (expressed in year, so 1/12 equals one month) between selected observation and the seasonal maximum}
}
\value{
vector of observations (time series) with annual observation frequency
}
\description{
Convert time series to annual frequency: The toAnnualTS function converts a time series with n observations per year to an annual time series (one observation per year). The main concept is to select observations per year closest to a given day of year that have no missing value (NA). Here, the day of year for which the seasonality is maximum is being used.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_characterize.R
\name{decompTSbfast}
\alias{decompTSbfast}
\title{Decompose time series into trend, seasonality and remainder: This function decomposes time series into three components using BFAST01 functionality: trend, seasonality and remainder. Trends are fitted using linear regression without breaks, seasonality is fitted using a first order harmonic function and the remainder equals the anomalies (i.e. time series - trend - seasonality).}
\usage{
decompTSbfast(df, nyr, nobsYr)
}
\arguments{
\item{df}{a dataframe with time series that need to be decomposed. The dataframe needs to be structured as follows: each row represents a sampled pixel. The first two columns contain the latitude and longitude of the pixel. The next columns contain the time series values for each observation date.}

\item{nyr}{number of years of the input time series}

\item{nobsYr}{number of observations per year of the input time series}
}
\value{
a list containing the estimated seasonality, remainder, trend and seasonality coefficients. The seasonality is a dataframe with the seasonality of each pixel. Each row represents a sampled pixel. The first two columns contain the latitude and longitude of the pixel. The next columns contain the seasonality values for each observation date. The trend and remainder are dataframes with the trend and remainder of each pixel (dataframe is structured in the same way as the seasonality). Seasonality_coefficients is a dataframe with the coeficients of the fitted harmonic function. Each row represents a sampled pixel. The first two columns contain the latitude and longitude of the pixel. The next columns contain the coefficients of the fitted harmonic function.
}
\description{
Decompose time series into trend, seasonality and remainder: This function decomposes time series into three components using BFAST01 functionality: trend, seasonality and remainder. Trends are fitted using linear regression without breaks, seasonality is fitted using a first order harmonic function and the remainder equals the anomalies (i.e. time series - trend - seasonality).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts_generator.R
\name{simulCase}
\alias{simulCase}
\title{Simulation of nrep disturbance time series.}
\usage{
simulCase(
  nrep,
  nyr,
  nobsYr,
  nDr,
  seasAv,
  seasAmp,
  trAv,
  remSd,
  distMaglim,
  distTy,
  distReclim,
  mval,
  mvaldist,
  distType
)
}
\arguments{
\item{nrep}{number of time series to simulate}

\item{nyr}{number of years that need to be simulated}

\item{nobsYr}{number of observations per year that will be simulated}

\item{nDr}{number of drought years that are introduced [i.e. setting seasonality of a year equal to its minimum value]. These drought years are randomly chosen for each of the simulated time series.}

\item{seasAv}{average seasonality profile}

\item{seasAmp}{seasonality amplitude}

\item{trAv}{offset value of time series}

\item{remSd}{standard deviation of the remainder}

\item{distMaglim}{limits of the disturbance magnitude, should be a vector with the minimum and maximum value. If the minimum equals the maximum value, the disturbance magnitude is fixed for each simulated time series (and equal to the minimum value). When the minimum value does not equal the maximum value, a disturbance magnitude is randomly chosen in the given interval for each simulated time series.}

\item{distTy}{year of the disturbance. If distTy equals one, the disturbance will take place in the first year. The exact disturbance date (day or year) is randomly chosen per time series.}

\item{distReclim}{limits of the halftime period of the recovery [number of observations], should be a vector with the minimum and maximum value. If the minimum equals the maximum value, the recovery period is fixed for each simulated time series (and equal to the minimum value). When the minimum value does not equal the maximum value, a recovery period is randomly chosen in the given interval for each simulated time series.}

\item{mval}{number of missing values to be introduced. If mval equals NA, no missing values are introduced. For missing values with a random interval (see mvaldist), this should equal the fraction of missing values (mval equal to 0.1 will result in an NA value for 10 percent of the time series). For missing values having a regular interval, every mval observations one value is kept (eg for a daily time series, a mval equal to 5 will result in one observation every 5 days).}

\item{mvaldist}{the distribution of the missing values. Should equal 'random'or 'interval'.}

\item{distType}{the type of disturbance. piecewise refers to a linear decay function, while exponential refers to an exponential decay}
}
\value{
a list with the simulated time series, offest, seasonality, remainder, disturbance component and parameters used for the simulation. The time series (components) are stored as matrix where each row is a time series and the columns are associated with the observation numbers.
}
\description{
Simulation of nrep disturbance time series.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_recovery.R
\name{calcFrazier}
\alias{calcFrazier}
\title{Calculate recovery metrics from a time series with known disturbance date. The calcFrazier function derives the RRI, R80P and YrYr recovery indicators,
defined by Frazier et al. (2018). The indicators are originally developped for annual long-term time series of optical vegetation indices.
Yet, in order to be able to derive the indicators as well for dense and/or short time series, a modified version is suggested.
Here, the user can define the time period before, during and after the disturbance that is used to derive the indicators.
To reduce the interference of the seasonal pattern of dense time series, the chosen time period should cover blocks of n years.
(Frazier, R. J., Coops, N. C., Wulder, M. A., Hermosilla, T., & White, J. C. (2018). Analyzing spatial and temporal variability in short-term rates
of post-fire vegetation return from Landsat time series. Remote Sensing of Environment, 205, 32-45.)}
\usage{
calcFrazier(tsio, tdist, obspyr, nPre, nDist, nPostMin, nPostMax)
}
\arguments{
\item{tsio}{vector of observations (time series with a fixed observation frequency)}

\item{tdist}{observation number of disturbance, indicating the timing of the disturbance}

\item{obspyr}{number of observations per year}

\item{nPre}{number of years prior to the disturbance used to calculate the pre-disturbance value}

\item{nDist}{number of years used to quantify the time series value during the disturbance}

\item{nPostMin}{the post-disturbance condition is quantified starting from nPostMin years after the disturbance}

\item{nPostMax}{max number of years after the disturbance used to quantify the post-disturbance condition}
}
\value{
a list containing the RRI recovery indicator, R80p recovery indicator and YrYr recovery indicator
}
\description{
Calculate recovery metrics from a time series with known disturbance date. The calcFrazier function derives the RRI, R80P and YrYr recovery indicators,
defined by Frazier et al. (2018). The indicators are originally developped for annual long-term time series of optical vegetation indices.
Yet, in order to be able to derive the indicators as well for dense and/or short time series, a modified version is suggested.
Here, the user can define the time period before, during and after the disturbance that is used to derive the indicators.
To reduce the interference of the seasonal pattern of dense time series, the chosen time period should cover blocks of n years.
(Frazier, R. J., Coops, N. C., Wulder, M. A., Hermosilla, T., & White, J. C. (2018). Analyzing spatial and temporal variability in short-term rates
of post-fire vegetation return from Landsat time series. Remote Sensing of Environment, 205, 32-45.)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_eval_performance.R
\name{mape}
\alias{mape}
\title{MAPE}
\usage{
mape(val, meas)
}
\arguments{
\item{val}{vector of x values}

\item{meas}{vector of y values}
}
\value{
the MAPE of the two vectors
}
\description{
MAPE
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_recovery.R
\name{calcBFASTrec}
\alias{calcBFASTrec}
\title{Post-disturbance slope and recovery metrics derived from BFAST0n trend segments. The calcBFASTrec function derives a set of recovery indicators after fitting a segmented trend in the time series. Using the breakpoints function of the strucchange package, a segmented trend is fitted (hereafter called BFAST0n trend segments). The detected break showing the largest change (in absolute values) is assumed to represent the disturbance. Using the segmented trend and detected disturbance date, the RRI, R80p, YrYr and the slope of the post-disturbance trend segment are derived as recovery indicators.}
\usage{
calcBFASTrec(tsio, obspyr, h, nPre, nDist, nPostMin, nPostMax, seas = F)
}
\arguments{
\item{tsio}{vector of observations (time series)}

\item{obspyr}{number of observations in one year}

\item{h}{This parameter defines the minimal segment size either given as fraction relative to the sample size or as an integer giving the minimal number of observations in each segment.}

\item{nPre}{number of years prior to the disturbance used to calculate the pre-disturbance value}

\item{nDist}{number of months used to quantify the time series value during the disturbance}

\item{nPostMin}{min number of years after the disturbance used to quantify the recovery}

\item{nPostMax}{max number of years after the disturbance used to quantify the recovery}

\item{seas}{TRUE or FALSE, include seasonal term when detecting breaks?}

\item{breaks}{'BIC' or 'LWZ': criteria used to define the optimal number of breaks (desactivated)}
}
\value{
a list containing  the RRI, R80p, YrYr recovery indicator derived from the BFAST0n trend segments and slope of the trend segment after the disturbance (sl).
}
\description{
Post-disturbance slope and recovery metrics derived from BFAST0n trend segments. The calcBFASTrec function derives a set of recovery indicators after fitting a segmented trend in the time series. Using the breakpoints function of the strucchange package, a segmented trend is fitted (hereafter called BFAST0n trend segments). The detected break showing the largest change (in absolute values) is assumed to represent the disturbance. Using the segmented trend and detected disturbance date, the RRI, R80p, YrYr and the slope of the post-disturbance trend segment are derived as recovery indicators.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts_generator.R
\name{simulTS}
\alias{simulTS}
\title{Simulate one time series with disturbance}
\usage{
simulTS(
  nyr,
  nobsyr,
  tMiss,
  nDr,
  seasAv,
  seasAmp,
  trAv,
  remSd,
  distMag,
  distT,
  distRec,
  distType
)
}
\arguments{
\item{nyr}{number of years that need to be simulated}

\item{nobsyr}{number of observations per year that will be simulated}

\item{tMiss}{timing of missing values [observation number]. If tMiss equals NA, no missing values are introduced.}

\item{nDr}{number of drought years that are introduced [i.e. setting seasonality of a year equal to its minimum value]. These drought years are randomly chosen.}

\item{seasAv}{average seasonality profile}

\item{seasAmp}{seasonality amplitude}

\item{trAv}{offset value of time series}

\item{remSd}{standard deviation of the remainder}

\item{distMag}{magnitude of the disturbance}

\item{distT}{timing of the disturbance [observation number]}

\item{distRec}{duration of the recovery [number of observations]}

\item{distType}{type of disturbance-recovery process: 'piecewise' represents a step function with linear recovery, 'exponential' represents an exponential decay}

\item{remMod}{ARMA model of remainder}
}
\value{
a list containign the years for which a drought was introduced and a time series object, containing the simulated seasonality, trend, remainder, disturbance, and the sum of these components.
}
\description{
Simulate one time series with disturbance
}
