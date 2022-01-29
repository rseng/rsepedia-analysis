# CategoricalTimeSeries.jl


| **Documentation**| **Appveyor** |
|:---------------:|:---------------:|
|[![Documentation Status](https://readthedocs.org/projects/categoricaltimeseriesjl/badge/?version=latest)](https://categoricaltimeseriesjl.readthedocs.io/en/latest/?badge=latest)| [![Build status](https://ci.appveyor.com/api/projects/status/ik7hvhu73kpvpr0r?svg=true)](https://ci.appveyor.com/project/johncwok/categoricaltimeseries-jl)| 

Toolbox helpful for the study of *categorical* time-series.
It contains methods used in [spectral envelope](https://categoricaltimeseriesjl.readthedocs.io/en/latest/Spectral_properties/)  analysis, [serial dependence](https://categoricaltimeseriesjl.readthedocs.io/en/latest/Correlations/) analysis, [clustering](https://categoricaltimeseriesjl.readthedocs.io/en/latest/Data_clustering/), and [motif recognition](https://categoricaltimeseriesjl.readthedocs.io/en/latest/Motif_recognition/). 

## Installation
The package can be installed via:
```Julia
using Pkg
Pkg.add("CategoricalTimeSeries")
```

## Documentation
The [documentation](https://categoricaltimeseriesjl.readthedocs.io/en/latest/) is available at https://categoricaltimeseriesjl.readthedocs.io/en/latest/ and provides comprehensive explanations, examples and a descriptive list of all usefull functions. 
## Issues
If you are experiencing troubles/bugs with some functionallities of the package, please open an issue and I will try to resolve it.
Alternatively, you can also open an issue to give me feedback or suggest new features.

## List of improvements
Despite all my efforts to make this package as complete as possible, I haven't had the time to implement certain improvements.
If you want to help, your contribution will be appreciated!

Here is a list of potential improvements to the existing methods:
- Implement common windowing functions for the **spectral envelope** method: The original paper post processes the results with a smoothing triangular function (of width ```m```), however the core estimation is still based on the periodogram. This effectively passes the data through a normalized boxcar function which is known to introduce bias.
As is the case in power-spectral density estimation, (e.g. Welch's method) a carefully choosen window such as the Hanning, Hamming or Blackman function could help reducing the bias and improve the accuracy of the results. An even better improvement would be to use the window functions of the slepian sequence, to allow averaging without frequency range restriction. However, introducing such windows will change the mathematical starting point of the spectral envelope and I do not know if closed-form solutions can be found for all windows.
- Implement simulated annealing for the **Information bottleneck** method: due to the probabilistic nature of this algorithm, results do not always converge to the absolute minima of the cost function. 
At the moment, the package samples n optimizations and selects the one that has the lowest cost. This is a working approach, but can be computationally costly depending on the amount of datapoints and categories. Simulated annealing could be a way to reduce computation time.
- Find a way to display the results of the **Information bottleneck** clusterings in a more visually comprehensive manner. At the moment, a DataFrame filled with 1 and 0s indicating appartenance to different clusters is used,
which can be hard to interpret when many categories are involved.
- add doc entry for ```apply_mapping```` function and correct export from apply_mappings to apply_mapping.

---
title: 'CategoricalTimeSeries.jl: A toolbox for categorical time-series analysis'
tags:
  - Julia
  - Categorical Time-series
  - Spectral analysis
  - Association measurement
  - Clustering

authors:
 - name: Corentin Nelias
   orcid: 0000-0001-6266-5575 
   affiliation: "1, 2"
affiliations:
 - name: Max Planck Institute for Dynamics and Self-Organization
   index: 1
 - name: Department of Physics, Georg-August-Universität Göttingen
   index: 2
date: 8 septembre 2021
bibliography: paper.bib
---

# Introduction

**CategoricalTimeSeries.jl** is a [Julia](https://github.com/JuliaLang/julia) toolbox made for analysing categorical time-series. 

Categorical time-series are time-sequenced data in which the values at each time point are categories rather than measurements.
The common approach to deal with categorical time-series consists in transforming the data via a mapping to obtain a real-valued sequence.
This enables the use of traditional time-series analysis methods. However, most of these methods (power-spectral density estimation, correlation coefficients, dimensionality reduction, etc.)
are not invariant under general transformations and will produce different results based on the choice of mapping. 
Therefore, depending on the type of categorical data and the problem at hand, it is desirable to have methods that work with the direct categorical values themselves. 

The purpose of **CategoricalTimeSeries.jl** is to provide such tools. The package comes with extensive documentation available online: https://categoricaltimeseriesjl.readthedocs.io/en/latest/

# Statement of need

While several implementations of categorical time-series analysis methods are already available, they are written in different languages, some of which are not free (e.g., Matlab). Additionally, no implementations for methods such as the *spectral envelope* or the *random projection* (see *Overview of functionality* below) are available online. This package centralizes and implements most of the standard methods of categorical time-series analysis in a single toolbox fully written in the Julia language. 


# Overview of functionality

This toolbox was designed to be easy to use and to produce results that are simple to plot. 
Consequently, the methods implemented in the package take the inputs as 1-D arrays of any type. 
Type conversion and pre-processing (when needed) are done automatically within the methods without the need for additional coding by the user.
The results are either formatted in a way that can be plotted directly with the  ```Plots.jl``` library, or a helper function is provided for visualization and interpretation.

The main areas of functionality are:

**Spectral analysis**:
The spectral envelope method [@Stoffer:1998] is used to study the power-spectrum of categorical time-series. 
As stated in the Introduction section, the power-spectrum of a time-series is not invariant under a generic transformation. 
A wrong choice of mapping can potentially flatten certain peaks and render them unnoticeable.
For each frequency, the spectral envelope seeks the mapping that maximizes the value of the power-spectrum normalized by the total variance.
The ```spectral_envelope``` function takes a time-series (1-D array) as the input and returns all the frequencies of the spectrum and the values of the intensity associated with the optimal mappings. It also returns the mappings.
For a finer study of the mappings themselves, the ```get_mappings``` function can be used, instead.

**Association analysis**:
The notion of auto-correlation function is not formally defined for a categorical time-series [@Weiss:2018].
Yet it might be of interest to know how interdependent the values of the time-series are. 
We implemented several coefficients generalizing the concept of linear correlations to categorical time-series.
Cramer's coefficient, Cohen's coefficient, and Theil's U can be computed via the ```cramer_coefficient```, ```cohen_coefficient```, and ```theils_u``` functions, respectively.
They take a 1-D array representing the time-series to study and an array of lags storing the lag values at which the coefficients are evaluated as the inputs. 

**Motif recognition**:
Time-series can present repeating motifs that are worthwhile identifying. However, simple line-search algorithms are not adapted for all motifs [@Pevzner:2000].
Moreover, the lack of proper distance measurement complicates the search in the context of categorical time-series.  
An implementation using the *random projection* method [@Buhler:2002] is used here.
The identification of potential motifs is performed by the ```detect_motifs``` function.
It takes a time-series (1-D array), the length of the motifs to look for, and the number of allowed errors as input arguments. 
It returns an instance of the ```pattern``` structure which stores properties of the identified motif such as shapes, repetition number, and positions.


**Data clustering**:
If certain categories in a time-series present functional similarities, one might wish to cluster them together into a single equivalent representation.
This reduces the total number of categories and can simplify the analysis of the time-series. For this purpose, we use an implementation based on the *Information bottleneck* concept [@Tishby:2000; @Strouse:2016].
After an initial bottleneck model of the structure ```IB``` is instantiated, it can be optimized with the ```IB_optimize!``` function to reveal potential clusters of categories. An overview of the results can be obtained with the ```print_results``` function. 

# Acknowledgements

The author thanks Nori Jacoby for discussing and providing insight on the *Information bottleneck* concept.

# References
# Motif recognition

Time-series sometimes present **repeating motifs** (or patterns) that are worthwhile identifying. The detection of such motifs can be difficult depending on the amount of noise in the time-series. <br/>

In the case of categorical time-series, the lack of a proper metric to measure distance between motifs can make their detection tricky. Improper distances like the number of differences between the two motifs is commonly used.

This package proposes a detection algorithm based on JEREMY BUHLER and MARTIN TOMPA's paper "[Finding Motifs Using Random Projections](https://pubmed.ncbi.nlm.nih.gov/12015879/)". This algorithm although very precise is not exact. Therefore, when you are done detecting potential motifs with the `detect_motifs` function, you can refine your results with `find_motifs` for an exact search.
<br/> The main functions return instances of a class called **pattern**:
- - -
**pattern — Class**
- - -
A class storing useful information about found motifs in a time-series. An array of `pattern` instances is returned when the searching algorithm is done running.
>**Attributes**:

- **shape** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): Array containing the shape (or contour) of the first found repetition of the motif.
- **instances** ([Array{Array{Any,1},1}](https://docs.julialang.org/en/v1/base/arrays/)): all the different shapes from the motif's repetitions, they can vary a bit from one to the next.
- **positions** ([Array{Int,1}](https://docs.julialang.org/en/v1/base/arrays/)): the positions at which the different repetitions of the motif were found.

## Main functions
- - -
**detect_motifs — Function**
- - -
```Julia
detect_motifs(ts, w, d, t = w - d; iters = 1000, tolerance = 0.95)
```
Detects all motifs of length 'w' occuring more often than chance, being identical to each other up to 'd' differences inside of input time-series 'ts'.
Returns an array of `pattern`, inside of which the patterns are classified by how frequently they are observed. The first elements is therefore the most frequently observed motif, and so on.
> **Parameters**:

>>* **ts** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): input time-series in which motifs are searched for.
>>* **w** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): length of motifs to look for.
>>* **d** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): allowed errors (differences) between motifs repetitions.
>>* **t = w - d** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): size of the masks to use for random projection in the detection (defaults to w - d).
>>* **iters = 1000** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): the numbers of iterations for the random projection process (defaults to 1000)
>>* **tolerance = 0.95** ([Float](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): threshold of motif identification. If set to 1, only matrix entries that are strictly superior to the (probabilistic) threshold are taken into account. Defaults to 0.7, meaning that matrix entries need to be bigger than 0.7*threshold.

> **Returns** :
>>* **motifs** : list of `pattern` instances sorted by frequency of occurence. motifs[1] is therefore the most frequent motif, motifs[2] the second most observed and so on.

- - -
**find_motifs — Function**
- - -
```Julia
find_motifs(ts, shape, d)
```
Given a motif of shape 'shape' (array{any,1}), looks for all the repetitions of it which differ only up to 'd' differences inside of the input time-series 'ts'.
Input:

>**Parameters**:
>>* **ts** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)) : time-series in which to look for motifs
>>* **shape** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): shape (aray{any,1}) of the motif to look for.
>>* **d** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): allowed errors (differences) between motifs

>* **Returns** :

>>* **motif** : an instance of `pattern` containing the found repetition of the input 'shape'.

## Plotting
To help visualize results, two simple plotting functions are provided.
- - -
**plot_motif — Function**
- - -
```Julia
plot_motif(m::pattern)
```
Plots all repetitions of an input `pattern` instance on top of each other to see how similar they are to each other.
> **Parameters**:

>>* **m** : Instance of the `pattern` class

- - -
**plot_motif — Function**
- - -
```Julia
plot_motif(m::pattern, ts)
```
Plots all repetitions of an input `pattern` instance on top of the input time-series 'ts' to better visualize their repartition in time.
> **Parameters**:

>>* **m** : Instance of the `pattern` class
>>* **ts** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): Input time-series


## Example
From Michael Brecker's improvisation over the piece ["confirmation"](https://github.com/johncwok/CategoricalTimeSeries.jl/tree/main/test), we extract a time-series of pitch intervals (difference from one note to the next).
A spectral envelope analysis reveals a peak at period 6~7, so we look for motifs of length 7 and allow for 1 error between them.
After detection, we visualize the most frequent motif:
```
using DelimitedFiles
using CategoricalTimeSeries

data_path = joinpath(dirname(dirname(pathof(CategoricalTimeSeries))), "test", "confirmation")
data = readdlm(data_path)
pitch = mod.(data, 12) #Removing octave position: not needed
intervals = pitch[2:end] .- pitch[1:end-1] #getting interval time-series.
m = detect_motifs(intervals, 7, 1; iters = 700, tolerance = 0.7)
plot_motif(m[1]) #plotting most frequent motif
```

<img src=https://user-images.githubusercontent.com/34754896/104308882-9c2c9e80-54d1-11eb-8882-cc31b7b2af8b.PNG width = "500">

We notice that the motif `[-1, -2, 10, -10, 2, 3, 5]` seems to be the underlying (consensus) shape. In musical notation, this motif would look like this (written in C major):
<img src=https://user-images.githubusercontent.com/34754896/104315350-1ca3cd00-54db-11eb-864d-3a1da9d5efeb.PNG width = "500">

We do an exact search with 1 error allowed to check if our previous detection missed any repetitions, and plot the found motif on top of each other:
```
consensus_shape = [-1, -2, 10, -10, 2, 3, 5]
motif = find_motifs(intervals, consensus_shape, 1)
plot_motif(motif)
```
<img src=https://user-images.githubusercontent.com/34754896/104308882-9c2c9e80-54d1-11eb-8882-cc31b7b2af8b.PNG width = "500">

Here, we obtain the same plot as before but this is not necessarily always the case. Knowing the consensus motif usually allows to find its repetitions more efficiently.

Now, we visualize the repetitions of the motif in the time-series:
```
plot_motif(motif, data)
```
<img src=https://user-images.githubusercontent.com/34754896/136664810-50cb437a-6924-4ba8-a562-c6e6784affcd.PNG width = "800">
# Miscalleneous
Here we regroup all additional useful functions that do not necessarily deserve a section of their own.

- - -
**rate_evolution — Function**
- - -
```
rate_evolution(series)
```
The rate of evolution is a way to test the stationarity of a categorical time-series.
If the rate evolves more or less linearly, then the time-series can reasonably be considered stationary.
It is most informative to plot the rate of evolution of each categories on the same graph for a direct visual inspection.
> **Parameters**:

>>* **series** ([Array{any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array of categorical time-series.

> **Returns**: `RATE`, Array containing, for each category, an array representing it's rate of evolution.

- - -
**LaggedBivariateProbability — Function**
- - -
```
LaggedBivariateProbability(serie, Lags::Array{Int64,1}, Category1, Category2)
```
Returns the lagged bivariate probability of two given categories, Pij.
Given i and j two categories, and l a lag (or array of lags),
Pij is the probability to have the category j at time t + l, if we have i at time t.
> **Parameters**:

>>* **serie** ([Array{any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array of categorical time-series.
>>* **category1**
>>* **category2**

> **Returns**: `pij`, Array containing, for each value in `lags`, the lagged bivariate probability.

- - -
**varcov — Function**
- - -
```Julia
varcov(ts::Array{Float64,2})
```
Computes the covariance-variance matrix of a given multivariate time-series. This can also be used for a univariate time-series but the input should still be 2-D.
> **Parameters**:

>>* **ts** ([Array{Float,2}](https://docs.julialang.org/en/v1/base/arrays/)): 2-D input array of multivariate time-series.

> **Returns**: `cov_matrix` the correpsonding covariance matrix.

- - -
**power_spectrum — Function**
- - -
```Julia
power_spectrum(x::Array{Float64,1}, window::Int, step::Int)
```
Computes an estimation of the power-spectrum of the input time-series `x`.
> **Parameters**:

>>* **x** ([Array{Float,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array of real-valued time-series.
>>* **window** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): Integer specifying the size of the window for averaging Must be shorter than length(x). Recommended value is 1/10th of length(x).
>>* **step** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): Parameters controlling the overlap between the windows. Shouldn't be biggger than div(window,2).  

> **Returns**: `pxx`, the estimated power-spectrum.
# Categorical Time-Series Analysis

## Introduction
**CategoricalTimeSeries.jl** is a Julia package regrouping methods of categorical time-series analysis. The term *categorical* commonly refers to two types of data: *nominal* and *ordinal*. <br/>
**Nominal**, or labeled values represent *discrete* units that have no intrinsic *order*, like common types of pet:

- Cat
- Dog
- Bird

**Ordinal** values on the other hand represent *discrete* and *ordered* units, like the size of a coffee cup:

1. Small
2. Medium
3. Large

When categorical data is layed out in function of time, one speaks of *categorical time-series*. <br/>
<br/>
Often, especially when dealing with ordinal time-series, it is enough to map the different values to a set of integers to carry the analysis. However, when it is not sufficient **CategoricalTimeSeries.jl** is here to help.
## Overview
The package lets you carry four main kind of analysis: **Spectral analysis**, **Data clustering**, **correlations analysis** and **motif recognition**. Other functionnalities are also avalaible (see misc.). These methods are agnostic to the type of data used (ordinal or nominal) as they do not rely on a pre-established ordering. Here is a quick overview of these methods, for more details go to the specific sections.
#### Spectral analysis
The standard approach to study spectral properties in categorical time-series is to map the different values to a set of numbers. While the overall shape of the spectrum is usually unaffected by this operation, peaks representing cyclic behaviors in the time-series can completely disappear depending on the choice of mapping. To tackle this issue, one can use the spectral envelope method. It was developed by *David S. Stoffer*  in order to identify optimal mappings.
#### Data clustering
A categorical time-series might have many different possible values, but these values are rarely independant. Some categories can be loosely equivalent to one-another. For example, among the values `"male"`, `"female"`, `"house"` and `"car"`, it is evident that `"male"` and `"female"` refer to similar concepts and could be grouped in one single category. Data clustering is the task of identifying such relationships in the time-series in order to reduce it to a more efficient representation. Here, an implementation of the *information bottleneck* concept is used to this end.
#### Motif recognition
Time-series sometimes present repeating motifs (or patterns) that are worthwhile identifying. In categorical time-series, the lack of proper distance measurement complicates this task, nonetheless several methods have been developed to this end. An implementation using the concept of *random projection* is used here.
#### Correlation analysis
The notion of autocorrelation function is formally not defined for a categorical time-series. Yet, it might be of interest to know how inter-dependent the values of the time-series are. Efforts have been made to generalize the concept of linear correlations to categorical time-series. This package implements several of these methods.

##Installation
The source code is available on [GitHub](https://github.com/johncwok/CategoricalTimeSeries.jl), otherwise,
**CategoricalTimeSeries.jl** can be installed with
```Julia
using Pkg
Pkg.add("CategoricalTimeSeries")
```
To use it, you need to import it:
```Julia
using CategoricalTimeSeries
```
# Correlations
The study of categorical data prevents the usage of standard tools like the autocorrelation function, as they are often not defined. The following functions provide ways to study categorical serial dependences.  
Most of these methods are described in C. Weiss's book "[An Introduction to Discrete-Valued Time Series](https://onlinelibrary.wiley.com/doi/book/10.1002/9781119097013)" (2018).
## Main functions
- - -
**cramer_coefficient — Function**
- - -
```Julia
cramer_coefficient(series, lags)
```
Measures average association between elements of ```series``` at time t and time t + ```lags```. Cramer's V is an unsigned measurement : its values lies in [0,1], 0 being perfect independence and 1 perfect dependence. k can be biased, for more informations, refer to [1].
> **Parameters**:

>>* **series** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **lags** ([Array{Int,1}](https://docs.julialang.org/en/v1/base/arrays/)): lag values at which cramer's coefficient is computed. Alternatively, `lags` can be an integer, a single integer value will then be returned.  

> **Returns**: `V`, the value of cramer's coefficient for each value in `lags`.

- - -
**cohen_coefficient — Function**
- - -
```Julia
cohen_coefficient(series, lags)
```
Measures average association between elements of ```series``` at time t and time t + ```lags```.
Cohen's k is a signed measurement : its values lie in [-pe/(1 -pe), 1], with positive (negative) values indicating positive (negative) serial dependence at `lags`. pe is probability of agreement by chance.

> **Parameters**:

>>* **series** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **lags** ([Array{Int,1}](https://docs.julialang.org/en/v1/base/arrays/)): lag values at which Cohen's coefficient is computed. Alternatively, `lags` can be an integer, a single integer value will then be returned.  

> **Returns**: `K`, the value of Cohen's coefficient for each value in `lags`.

- - -
**theils_u — Function**
- - -
```Julia
theils_u(series, Lags)
```
Measures average portion of information known about `series` at t + `lags` given that `series` is known at time t. Theil's U makes use of concepts borrowed from *information theory*
U is an unsigned measurement: its values lies in [0,1], 0 meaning no information shared and 1 complete knowledge (determinism).

> **Parameters**:

>>* **series** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **lags** ([Array{Int,1}](https://docs.julialang.org/en/v1/base/arrays/)): lag values at which Theil's U is computed. Alternatively, `lags` can be an integer, a single integer value will then be returned.  

> **Returns**: `U`, the value of Theil's U for each value in `lags`.

## Confidence interval
Depending on the length of the time-series and the method used, the estimated value of serial dependence might fluctuate a lot around its true value.
It is therefore useful to relate estimations to a corresponding confidence interval to know how significant given results are. The following function provides a confidence interval via bootstrap:
- - -
**bootstrap_CI — Function**
- - -
```Julia
bootstrap_CI(series, lags, coef_func, n_iter = 1000, interval_size = 0.95)
```
Returns a top and bottom limit of the confidence interval at values of `lags`. The width of the confidence interval can be choosen (defaults to 95%). The returned confidence interval corresponds to the null hypothesis (no serial dependence), if the estimated serial dependence lies in this interval, no significant correlations can be claimed.

> **Parameters**:

>>* **series** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **lags** ([Array{Int,1}](https://docs.julialang.org/en/v1/base/arrays/)): lag values at which the CI is computed.
>>* **coef_func** ([function](https://docs.julialang.org/en/v1/manual/functions/)): the function for which the CI needs to be computed.
            `coef_func` can be one of the following **functions** : `cramer_coefficient`, `cohen_coefficient` or `theils_U`.
>>* **n_iter** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): number of iterations for the bootstrap procedure. The higher, the more precise but more computationaly demanding. Defaults to 1000.
>>* **interval_size** ([Float](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): Desired size of the confidence interval. Defaults to 0.95, for a 95% confidence interval.

> **Returns**: `(top_values, bottom_values)`, the top and bottom limit for confidence interval, for each point in `lags`.

## Example
Using the Pewee [birdsong data](https://github.com/johncwok/CategoricalTimeSeries.jl/tree/main/test) (1943) one can do a serial dependence plot using Cohen's cofficient as follow :
```
using DelimitedFiles, Plots
using CategoricalTimeSeries

#reading 'pewee' time-series test folder.
data_path = joinpath(dirname(dirname(pathof(CategoricalTimeSeries))), "test", "pewee.txt")
series = readdlm(data_path,',')[1,:]
lags = collect(1:25)
v = cohen_coefficient(series, lags)
t, b = bootstrap_CI(series, lags, cohen_coefficient)
a = plot(lags, v, xlabel = "Lags", ylabel = "K", label = "Cohen's k")
plot!(a, lags, t, color = "red", label = "Limits of 95% CI"); plot!(a, lags, b, color = "red", label = "", dpi = 600)
```
<img src=https://user-images.githubusercontent.com/34754896/139043402-7321dfda-c741-473d-bcb2-57a5ee217946.png width = "600">
# Spectral Envelope

The **spectral envelope** is a tool to study cyclic behaviors in categorical data. It is more informative than the traditional approach of attributing a different number to each category for power-spectral density estimation. <br/>

For each frequency in the spectrum, the **spectral envelope** finds an optimal real-numbered mapping that maximizes the normed power-spectral density at this point. Therefore, no matter what mapping is choosen for the different categories, the power-spectral density will always be bounded by the spectral envelope.

The spectral envelope was defined by David S. Stoffer in *DAVID S. STOFFER, DAVID E. TYLER, ANDREW J. MCDOUGALL*, [Spectral analysis for categorical time series: Scaling and the spectral envelope](https://www.jstor.org/stable/2337182).

## Main functions
- - -
**spectral_envelope — Function**
- - -
```Julia
spectral_envelope(ts; m = 3)
```
Computes the spectral envelope of an input categorical time-series.  
The degree of smoothing can be chosen by the user.

> **Parameters**:

>>* **ts** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **m** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): Smoothing parameter. corresponds to how many neighboring points
        are to be involved in the smoothing (weighted average). Defaults to 3.  

> **Returns**: `(freq, se, eigvecs)`, with `freq` the frequencies of the power-spectrum, `se` <br/> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the values of the spectral envelope for each frequency in 'freq'.
    `eigvecs` contains <br/> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the optimal real-valued mapping for each frequency point.

- - -
**get_mappings — Function**
- - -
```
get_mappings(data, freq; m = 3)
```

Computes, for a given frequency `freq`, the optimal mappings for the categories in `data`. Scans the vincinity of `freq` to find the maximum of the spectral envelope, prints a sum up and returns the obtained mappings.
> **Parameters**:

>>* **data** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **freq** ([Float](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): Frequency for which the mappings are wanted. The vincinity of 'freq' will be scanned to find maximal value of the spectral envelope.  
>>* **m** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): Smoothing parameter. corresponds to how many neighboring points
        are to be involved in the smoothing (weighted average). Defaults to 3.  

> **Returns**: `mappings`, the optimal mappings for the found maxima around 'freq'.



## Example
Applying the spectral envelope to study a [segment of DNA](https://github.com/johncwok/CategoricalTimeSeries.jl/tree/main/test) from the Epstein-Barr virus and plotting the results:
```
using DelimitedFiles, Plots
using CategoricalTimeSeries

data_path = joinpath(dirname(dirname(pathof(CategoricalTimeSeries))), "test", "DNA_data.txt")
data = readdlm(data_path, ',')
f, se, eigvecs = spectral_envelope(data; m = 0)

plot(f, se, xlabel = "Frequency", ylabel = "Intensity", title = "test data: extract of Epstein virus DNA", label = "spectral envelope")
```
<img src=https://user-images.githubusercontent.com/34754896/136663948-a1ada6b7-691e-4e75-9fea-f905240c261e.PNG width = "500">

To get the associated optimal mapping for the peak at frequency 0.33:
```
mappings = get_mappings(data, 0.33; m = 0)
>> position of peak: 0.33 strengh of peak: 0.02
print(mappings)
>> Dict{SubString{String}, Float64} with 4 entries:
  "A" => -0.59
  "T" => 0.55
  "C" => 0.0
  "G" => 0.6
```
# Information bottleneck

The information bottleneck (IB) concept can be used in the context of categorical data analysis to do **clustering**, or in other words, to look for categories which have equivalent functions.  
Given a time-series, the IB looks for a concise representation of the data that preserves as much meaningful information as possible. In a sense, it is a lossy compression algorithm. The information to preserve can be seen as the ability to make predictions: given a specific context, how much of what is coming next can we predict ?
The goal of this algorithm is to cluster categorical data while preserving predictive power.  
To learn more about the information bottleneck you can look at [[1](https://arxiv.org/abs/1604.00268)] or [[2](https://doi.org/10.1080/09298215.2015.1036888)]

## Quick start
To do a simple IB clustering of a categorical time series, the first step is to instantiate an ```IB``` model. Then optimize it via the ```IB_optimize!``` function to obtain to obtain the optimal parameters.
```
data = readdlm("/path/to/data/")
model = IB(data) #you can call IB(x, beta). beta is a real number that controls the amount of compression.
IB_optimize!(model)
```
The data needs to be presented as a 1-D array, otherwise IB interprets it as a probability distribution (see below).

To see the results, you can use:
```
print_results(model)
```
Rows are clusters and columns correspond to the input categories. The result is the probability **p(t|x)** of a category belonging to a given cluster. Since most of the probabilities are very low, ```print_results``` **sets every p(t|x) > 0.1 to 1**. **p(t|x) < 0.1** are set to **0 otherwise** for ease of readability (see further usage for more options).
The optimized parameters may vary from one optimization to the next as the algorithm is not deterministic, to obtain global optima, use the ```search_optima``` function (see below).

## Further usage
To have a better grasp of the results produced by IB clustering, it is important to understand the parameters influencing the algorithm of **IB** model structures.
The two most important parameters are the amount **compression** and the definition of the **context**. They are provided upon instanciation:
- - -
**IB — Type**
- - -
```
IB(x, y, β = 100, algorithm = "IB")
IB(x, β = 100, algorithm = "IB")
IB(pxy::Array{Float64,2}, β = 100, algorithm = "IB")
```

> **Parameters**:

>>* **x** ([Array{Int,1} or Array{Float,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **y** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): Context used for data compression. If not provided, defaults to "next element", meaning for each element of x, y represent the next element in the series. This means that the IB model will try to preserve as much information between 'x' and it's next element. (see `get_y` function)
>>* **β** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): parameter controlling the degree of compression. The smaller `β` is, the more compression. The higher `β`, the bigger the mutual information I(X;T) between the final clusters and original categories is.
There are two undesirable situations: if `β` is too small, maximal compression is achieved and all information is lost. If `β` is too high, there is no compression.<br/>  with "IB" algorithm, a high `β` value (~200) is a good starting point. With "DIB" algorithm, `β` > ~5 can already be too high to achieve any compression. <br/> `β` values > ~1000 break optimization because all metrics are effectively 0.
>>* **algorithm** ([String](https://docs.julialang.org/en/v1/manual/strings/)): The kind of compression algorithm to use. "IB" choses the original IB algorithm (Tishby, 1999) which does *soft* clustering, "DIB" choses the *deterministic* IB algorithm (DJ Strouse, 2016) doing *hard* clustering. The former seems to produce more meaningfull clustering. Defaults to "IB".
>>* **pxy** ([Array{Float,2}](https://docs.julialang.org/en/v1/base/arrays/)): joint probability of element 'x' to occur with context 'y'. If not provided, is computed automatically. From `x` and `y`.


> **Returns**: instance of the `IB` mutable struct.


- - -
**get_y — Function**
- - -
```
get_y(data, type = "nn")

```
Defines and return the **context** associated with the input time-series `data`.
> **Parameters**:

>>* **data** ([Array{Any,1}](https://docs.julialang.org/en/v1/base/arrays/)): 1-D Array containing input categorical time-series.
>>* **type** ([String](https://docs.julialang.org/en/v1/manual/strings/)): type of context to use. Possible values are "nn" or "an". Defaults to "nn" (for *next neighbor*). This means, if data = ["a","b","c","a","b"], the "nn" context vector y is ["b","c","a","b"]. Chosing "an" (for adjacent neighbors) not only includes the next neighbor but also the previous neighbor, every element of y is then a tuple of previous and next neighbor.

> **Returns**: `y`, associated context to `data`.


## Additional functions
- - -
**calc_metrics — Function**
- - -
```
calc_metrics(model::IB)
```
Computes the different metrics (*H(T), I(X;T), I(Y;T)* and *L*) of an IB model based on its internal probability distributions.
> **Parameters**:
>>* **model**: an IB model

> **Returns**: (ht, ixt, iyt, L), metrics. ht is the entropy of the clustered representation. ixt is the mutual information between input data and clustered representation. iyt is the mutual information between context and clustered representation. L is the loss function.

- - -
**search_optima — Function**
- - -
```
search_optima!(model::IB, n_iter = 10000)
```
Optimization is not 100% guaranteed to converge to a **global maxima**. this function initializes and optimizes the provided `IB` model `n_iter` times, then, the optimization with the lowest `L` value is selected. The provided `IB` is updated in place. <br/>
> **Parameters**:

>>* **model**: an IB model
>>* **n_iter** ([Int](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): defined how many initialization/optimization are performed for the optima search.

> **Returns**: `nothing`. The update is done in place.

- - -
**print_results — Function**
- - -
```
print_results(m::IB, disp_thres = 0.1)
```
Displays the results of an optimized IB model.
> **Parameters**:

>>* **m**: an IB optimized model
>>* **disp_thres** ([Float](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): The probability threshold to consider that a category belongs to a given threshold. This makes reading the results more easy. Defaults to 0.1.

> **Returns**: `nothing`. Print the results.

If you want to get the **raw probabilities** `p(t|x)` after optimization (`print_results` filters it for ease of readability), you can access them with :
```
pt_x = model.qt_x
```
Similarly, you can also get p(y|t) or p(t) with `model.qy_t` and `model.qt`.<br/>


- - -
**get_IB_curve — Function**
- - -
```
`get_IB_curve(m::IB, start = 0.1, stop = 400, step = 0.05; glob = false)`
```
Scans the IB plane with various values of beta to get the optimal curve in the IB plane.
> **Parameters**:

>>* **m**: an IB optimized model
>>* **start** ([Float](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): The start β value.
>>* **stop** ([Float](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): The ending β value
>>* **step** ([Float](https://docs.julialang.org/en/v1/manual/integers-and-floating-point-numbers/)): The steps in β values that the function takes upon optimizing the provided model.
>>* **glob** ([Bool](https://docs.julialang.org/en/v1/manual/types/): if True, each optimization is done with the help of `search_optima` (more computationally demanding). Default to False.

> **Returns**: (ixt, iyt) the values of mutual information between data and clusters and context and clusters for each β value used by the function.


## Examples
- - -
Here is a concrete example with data from [Bach chorales](https://github.com/johncwok/CategoricalTimeSeries.jl/tree/main/test). The input categories are the 7 types of diatonic chords described in classical music theory. In this case, the data (input series and context) have already been compiled into a co-occurence table, so we instantiate the IB model with a probability distribution:

```
using CategoricalTimeSeries
using CSV, DataFrames

data_path = joinpath(dirname(dirname(pathof(CategoricalTimeSeries))), "test", "bach_histogram")
bach = DataFrame(CSV.File(data_path))
pxy = Matrix(bach)./sum(Matrix(bach)) #normalizing the co-occurence table to have probabilities.
model = IB(pxy, 20) #instantiating the model with co-occurence probabilities.
search_optima!(model, 1000)
print_results(model)
```

The output is in accordance with western music theory, but reveals interesting features. It tells us that we can group category 1, 3: this corresponds to the *tonic* function in classical harmony. Category 5 and 7 are joined: this is the *dominant* function.
Interestingly, category 2 and 4 (*subdominant* function) have not been clustered together, revealing that they are not used in a totally equivalent fashion in Bach's chorales.

<img src=https://user-images.githubusercontent.com/34754896/139694753-c4f97601-7ebe-4b70-8781-1e51e4deb842.PNG width = "400">

- - -
In the next example, we instantiate the model with a time-series ([saxophone solo](https://github.com/johncwok/IntegerIB.jl/tree/master/data)) and define our own context.

```
using CategoricalTimeSeries
using CSV, DataFrames

data_path = joinpath(dirname(dirname(pathof(CategoricalTimeSeries))), "test", "coltrane_afro_blue")
data = DataFrame(CSV.File(data_path))[!,1]  #time-series of notes from saxophone solo (John Coltrane).
context = get_y(data, "an") # "an" stands for adjacent neighbors.
model = IB(data, context, 500) # giving the context as input during instantiation.
IB_optimize!(model)
```
- - -
Now, we show how to plot the IB curve:

```
using Plots, CSV, DataFrames
using CategoricalTimeSeries

data_path = joinpath(dirname(dirname(pathof(CategoricalTimeSeries))), "test", "bach_histogram")
bach = DataFrame(CSV.File(data_path))
pxy = Matrix(bach)./sum(Matrix(bach)) #normalizing the co-occurence table to have probabilities.
model = IB(pxy, 1000) #instantiating the model with co-occurence probabilities.
x, y = get_IB_curve(model)
a = plot(x, y, color = "black", linewidth = 2, label = "Optimal IB curve", title = "Optimal IB curve \n Bach's chorale dataset")
scatter!(a, x, y, color = "black", markersize = 1.7, xlabel = "I(X;T) \n", ylabel = "- \n I(Y;T)", label = "", legend = :topleft)
```

<img src=https://user-images.githubusercontent.com/34754896/90395817-72438d00-e095-11ea-8872-3030db40539c.PNG width = "600">

## Acknowledgments
Special thanks to Nori Jacoby from whom I learned a lot on the subject. The IB part of this code was tested with his data and reproduces his results. <br/>
The present implementation is adapted from DJ Strouse's paper https://arxiv.org/abs/1604.00268 and his python implementation.


[1]: https://arxiv.org/abs/1604.00268
[2]: https://doi.org/10.1080/09298215.2015.1036888
