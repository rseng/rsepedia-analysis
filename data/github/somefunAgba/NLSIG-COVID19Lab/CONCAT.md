# NLSIG-COVID19Lab


![Actions Status](https://github.com/somefunAgba/NLSIG-COVID19Lab/actions/workflows/paper.yml/badge.svg)
[![GitHub issues](https://img.shields.io/github/issues/somefunAgba/NLSIG-COVID19Lab)](https://github.com/somefunAgba/NLSIG-COVID19Lab/issues)

[![License: BSD-3-License](https://img.shields.io/badge/License-BSD%203--Clause-success.svg)](https://github.com/somefunAgba/NLSIG-COVID19Lab/blob/main/LICENSE)
[![View NLSIG-COVID19Lab on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/84043-nlsig_covid19lab)
[![PWC](https://img.shields.io/endpoint.svg?url=https://paperswithcode.com/badge/the-nlogistic-sigmoid-function/covid-19-modelling-on-who)](https://paperswithcode.com/sota/covid-19-modelling-on-who?p=the-nlogistic-sigmoid-function)

A playground for descriptive modelling and monitoring of the time-series COVID-19 pandemic growth with the nlogistic-sigmoid function.

<img alt="NLSIG_COVID19LAB" src="nlsig_avatar.png"/>

nlogistic-sigmoid function (NLSIG) is a modern logistic-sigmoid function definition for modelling growth (or decay) processes. It features two logistic metrics (YIR and XIR) for monitoring growth from a two-dimensional (x-y axis) perspective.

## Papers (Links)
* [![DOI](https://joss.theoj.org/papers/10.21105/joss.03002/status.svg)](https://doi.org/10.21105/joss.03002)

* [NLSIG Conference Presentation Slides](nlsigcv19_confslide.pdf) *Best Student Paper* at the **2nd African Symposium on Big Data, Analytics and Machine Intelligence and 6th TYAN International Thematic Workshop December 3-4, 2020**.
 
* [NLSIG Preprint](https://arxiv.org/abs/2008.04210)

## Data Source
World Health Organization

## Getting Started : MATLAB

**Matlab Toolbox Requirements**
- Optimization Toolbox
- Statistics and Machine Learning Toolbox

**MATLAB Release Compatibility**
- Compatible with R2020a and later releases

**External Dependencies**
- None

### MATLAB App: Installation
Provided is a MATLAB App to allow for easy use. 

1. In the *applet folder*, double-click or right-click on the App Installer: *NLSIG-COVID19Lab.mlappinstall*

2. A confirmation dialog **Install into My Apps** pops up. Select *Install*

3. The App is then Installed, and can be accessed from the **MY APPS** section in the **APPS tab** on MATLAB's top panel or toolstrip.
**Hover** on the NLSIG-COVID19Lab App icon to see the **App details** and **install location**.

The default install location is in the Add-Ons location and is specific to the platform.
This can be viewed from the Home tab's Environment section, click **Preferences > MATLAB > Add-Ons**.

The install location then can be somewhere like this:

Windows - C:\Users\username\Documents\MATLAB\Add-Ons\Apps\NLSIGCOVID19Lab

Linux - ~/MATLAB/Add-Ons/Apps/NLSIGCOVID19Lab.

Mac - ~/Library/Application Support/MathWorks/MATLAB/Add-Ons/Apps/NLSIGCOVID19Lab

### MATLAB App: Using

Click the NLSIG-COVID19Lab App icon to start the App.

See the screenshots below

To **Start Modelling**: Click ![model](osspaper/play_24.png)

To **Update the Local Database**: Click ![model](osspaper/import_24.png)

To View Available Country-Codes: Switch to the **List:Country-codes** Tab

To Set options for Modelling: Click on the **Options** Tab
 
![GUI Layout showing the Total COVID-19 Infections of the World. \label{fig:iwdcigui}](osspaper/inf_wd_ci_gui.png)

![GUI Layout showing the Total COVID-19 Deaths of the World. \label{fig:dwdcigui}](osspaper/dth_wd_ci_gui.png)

#### Metrics: Interpretation
As at the 21st of February 2021:

*For infections*: 

**YIR = 0.822 [0.82, 0.843]** indicates that the numbers are past the peak; 

**XIR = 2.95 [2.88, 2.96]** indicates that this time is clearly a post-peak period. 

*For deaths*: 

**YIR = 0.664 [0.61, 0.718]** indicates that the numbers are no longer increasing and are recently past the peak for this phase; 

**XIR = 1.66 [1.41, 1.97]** indicates that this time is an early post-peak period. 

## Metrics

**R2GoF** R2 Goodness of Fit

**KSGoF** Kolmogorov-Smirnov Goodness of Fit 

**KSdist** Kolmogorov-Smirnov Distance Statistic

**YIR** Y-to-Inflection Ratio (Here Y = Infections or Deaths)

`YIR < 0.5` indicates generally increasing motion of growth

`YIR ~= 0.5` indicates generally that the increase has peaked. 

`YIR > 0.5` indicates generally reducing motion of growth

`YIR ~= 0` indicates either that the growth is flattening or could be increasing. 

**XIR** X-to-Inflection Ratio (Here X = Time in Days)

`XIR < 1` indicates a pre-peak period

`XIR ~= 1` indicates a peak-period. 

`XIR > 1` indicates a post-peak period.

`XIR ~= 0` indicates either a post-peak period or an early pre-peak. 

<!-- **Toy Example**

For infections: the YIR = 0.4916 [0.4908, 0.5063] indicates that the numbers are peaking and may start to decrease soon; the XIR = 0.9843 [0.9826, 1.0146] indicates that this time is close to a peak period. 

For deaths: the YIR = 0.4584 [0.4241, 0.5079] indicates that the numbers are still increasing but may likely peak soon; the XIR = 0.9266 [0.8634, 1.0245] indicates that this time is most-likely a peak period, close to a post-peak period. -->

## Frontend API 
Examples of Frontend APIs available for this software package can be found in the:

1. ``examples_m_api`` folder and 

2. ``examples_mlx_api`` folder   

<!-- 	You should see:
	'view_ccode.m'

	'upd_all.m'

	'query_single.m'

	'query_batch.m'

	'query_all.m'

	First, it is recommended to start with 'query_single.m'. 
	The country code for the world here is ``WD``.

	### 'view_ccode.m'
	View all country codes.
	Example: type ``view_ccode`` in the command window.

	### 'upd_all.m'
	Update data on the COVID-19 pandemic for all country codes. This needs
	a good internet connection.
	Example: type ``upd_all`` in the command window.

	### 'query_single.m'
	Query COVID-19 pandemic for selected country code.

	### 'query_batch.m'
	Query COVID-19 pandemic for a batch of selected country codes.

	### 'query_all.m'
	Query COVID-19 pandemic for all country codes. -->

## Saved Results
Saved model fit results and logistic metrics for infections and deaths can be found in the *assets* folder and *measures* folder

### *assets* folder
Stores all graphics for the model fit of infections and deaths in a folder named by the last date time-stamp in the data. 
Graphics are individually saved using the country code. 

For example: ``WDi.pdf`` and ``WDd.pdf`` respectively indicates the
saved graphics of the COVID-19 infections and deaths model fit for the World to the last date time-stamp in the data in pdf format.

### *measures* folder
Stores all estimated logistic metrics for infections and deaths till 
the last date time-stamp in the data in the *infs* and *dths* 
subfolders respectively.
	
## Automated Tests
Automated Tests for the app and package functionalties can be found in the *tests* folder.
	
<!-- #### Example
 --><!-- Running 'query_single.m' with the search_code as ``WD``
gave the following model fit for the ongoing COVID-19 pandemic with respect to the last updated date of the data. -->

<!-- **WORLD COVID-19 Infections**
<p align="center">
 <img alt="WDi" src="landing/WDi.png" width=500px/>
</p>

**WORLD COVID-19 Deaths**
<p align="center">
<img alt="WDd" src="landing/WDd.png" width=500px/>
</p>

**UK COVID-19 Infections**
<p align="center">
 <img alt="USi" src="landing/GBi.png" width=500px/>
</p>

**UK COVID-19 Deaths**
<p align="center">
<img alt="USd" src="landing/GBd.png" width=500px/>
</p>


**USA COVID-19 Infections**
<p align="center">
 <img alt="USi" src="landing/USi.png" width=500px/>
</p>

**USA COVID-19 Deaths**
<p align="center">
<img alt="USd" src="landing/USd.png" width=500px/>
</p>

**CHINESE COVID-19 Infections**
<p align="center">
 <img alt="CNi" src="landing/CNi.png" width=500px/>
</p>

**CHINESE COVID-19 Deaths**
<p align="center">
<img alt="CNd" src="landing/CNd.png" width=500px/>
</p>
 -->

<!--#### Recovered-->

 
## Contributions, Issues or Support
If you are:

- interested in dedicating the time to port to other languages or contribute to the software

- want to report issues or problems with the software

- seek miscellanous support

Then, you may contact me, by creating a [new issue](https://github.com/somefunAgba/NLSIG-COVID19Lab/issues/new/choose).

## License
This work is free software under the [BSD 3-Clause "New" or "Revised" License](https://github.com/somefunAgba/NLSIG-COVID19Lab/blob/main/LICENSE) 

## Citation Details
Somefun et al., (2021). NLSIG-COVID19Lab: A modern logistic-growth tool (nlogistic-sigmoid) for descriptively modelling the dynamics of the COVID-19 pandemic process. Journal of Open Source Software, 6(60), 3002, https://doi.org/10.21105/joss.03002
NLSIG-COVID19Lab: A nlogistic-sigmoid modelling laboratory 
                  for the COVID-19 pandemic growth
				  
Copyright (c) 2020, Oluwasegun Somefun

All rights reserved.

This program is free software; you can redistribute it and/or modify
it under the terms of the `ISC License` or `BSD 3-Clause License`.

See the LICENSE.txt file for the `BSD 3-Clause License`.

`ISC License` is stated below:

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

---
title: '`NLSIG-COVID19Lab`: A modern logistic-growth tool (nlogistic-sigmoid) for descriptively modelling the dynamics of the COVID-19 pandemic process'
tags:
  - Matlab
  - COVID-19
  - logistic function
  - machine learning
  - neural networks
  - optimization
  - regression
  - epidemiology
authors:
  - name: Oluwasegun A. Somefun^[Corresponding author.]
    orcid: 0000-0002-5171-8026
    affiliation: 1
  - name: Kayode F. Akingbade
    affiliation: 1
  - name: Folasade M. Dahunsi
    affiliation: 1
affiliations:
 - name: Federal University of Technology Akure, Nigeria
   index: 1
date: 15 December 2020
bibliography: paper.bib
---

# Summary

The growth (flow or trend) dynamics in any direction for most natural phenomena such as epidemic spread, population growth, 
adoption of new ideas, and many more, can be approximately modeled by the logistic-sigmoid curve. 
In particular, the logistic-sigmoid function with time-varying parameters is the core trend 
model in Facebook's Prophet model for time-series growth-forecasting 
at scale [@taylorForecastingScale2018] on big data. 
The scientific basis for this prevalence is given in [@bejanConstructalLawOrigin2011]. 
Such growth processes can be viewed as complex input--output systems that involve 
multiple peak inflection phases with respect to time, an idea that 
can be traced back in the crudest sense to [@reedSummationLogisticCurves1927]. A modern definition for the logistic-sigmoid growth, which considers restricted growth from  a two-dimensional perspective, is the nlogistic-sigmoid function (`NLSIG`) [@somefunLogisticsigmoidNlogisticsigmoidModelling2020] 
or logistic neural-network (`LNN`) pipeline. 

In this context, `NLSIG-COVID19Lab` functions as a `NLSIG` playground for descriptive modelling 
of the COVID-19 epidemic growth in each affected country of the world and in the world as a whole. 


# Statement of need

Epidemiological models such as the SEIRD variants 
[@leeEstimationCOVID19Spread2020;@okabeMathematicalModelEpidemics2020] are just another form of representing sigmoidal growth [@xsRichardsModelRevisited2012]. However, it has been noted 
[@christopoulosNovelApproachEstimating2020] that the SEIRD-variant models yield largely exaggerated forecasts. 
Observing the current state of the COVID-19 pandemic, this is concern is borne out, in 
the results of various applications of logistic modelling [@batistaEstimationStateCorona2020;@wuGeneralizedLogisticGrowth2020]
that have largely led to erroneous assessments of the epidemic's progress and its future projection, leading policymakers astray [@matthewWhyModelingSpread2020]. 

Notably, two recurring limitations of the logistic definitions in the literature and other software packages exist. These two limitations are trends that have persisted since the first introduction of the logistic-sigmoid function [@bacaerVerhulstLogisticEquation2011]. 

The first is that the co-domain of logistic function is assumed to be infinite. This assumption violates the natural principle of finite growth. 
The second is that during optimization, estimation of the logistic hyperparameters  for the individual logistic-sigmoids that make the multiple logistic-sigmoid sum are computed separately, instead of as a unified function. The effect of this is that as the number of logistic-sigmoids 
considered in the sum increases, regression analysis becomes more cumbersome and complicated, as can be observed in a number of works [@leeEstimationCOVID19Spread2020;@batistaEstimationStateCorona2020;
@hsiehRealtimeForecastMultiphase2006;@wuGeneralizedLogisticGrowth2020;
@chowellNovelSubepidemicModeling2019;@taylorForecastingScale2018]. 

These limitations are efficiently overcome by the nlogistic-sigmoid function `NLSIG` (or logistic neural-network pipeline) for describing logistic growth. We note that the `NLSIG` is a logistic neural-network machine-learning tool under active development. The benefits it provides at a functional level are:
	
 - unified function definition,
 
 - functional simplicity and efficient computation,
	
 - improved nonlinear modelling power.
		
Ultimately, the development of the `NLSIG-COVID19Lab` was motivated by research needs, in that it 
illustrates the power of the nlogistic-sigmoid neural pipeline. \linebreak `NLSIG-COVID19Lab` provides an optimization workflow with functions to make modelling and monitoring the COVID-19 pandemic easier and reliable. Notably, instead of engaging in false prophecy 
or predictions on the cumulative growth of an ongoing growth phenomena, whose source is both uncertain and 
complex to encode in current mathematical models [@christopoulosEfficientIdentificationInflection2016;@matthewWhyModelingSpread2020], this software package makes projections by means of:

- two-dimensional perspective metrics: Y-to-Inflection Ratio (YIR, here Y = Infections or Deaths); X-to-Inflection Ratio (XIR, here X = Time in Days) for robust monitoring of the growth-process being modelled in an area or locale of interest, and

- an adaptation of the Dvoretzky–Kiefer–Wolfowitz (DKW) inequality for the Kolmogorov–Smirnov (KS) test to construct a non-parametric confidence interval of uncertainty on the nlogistic-sigmoid model with a 99% probability ($\alpha=0.01$) by default. 

`NLSIG-COVID19Lab` is useful as a quick real-time monitoring tool for the COVID-19 pandemic. It was designed to be used by humans: both researchers and non-researchers. 

`NLSIG-COVID19Lab` is currently written in MATLAB but will be implemented in other programming languages in the future. 
 
The user-client end (both user application scripts and graphical user interface) of the `NLSIG-COVID19Lab` 
is designed to provide a friendly interface demonstrating the `NLSIG` modelling power for time-series growth processes from data. 
In this case, the growth-process is the time-series COVID-19 pandemic growth from official datasets (see \autoref{fig:iwdcigui} and \autoref{fig:dwdcigui}).

![GUI Layout showing the Total COVID-19 Infections in the World. \label{fig:iwdcigui}](inf_wd_ci_gui.png){ width=70% } 

![GUI Layout showing the Total COVID-19 Deaths in the World. \label{fig:dwdcigui}](dth_wd_ci_gui.png){ width=70% } 


### Core Data Source
At the time of writing, the COVID-19 database of `NLSIG-COVID19Lab` is sourced from the:

* World Health Organization

* Johns Hopkins University Center for Systems Science and Engineering


<!-- # Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x]$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array]{l]
0\textrm{ if ] x < 0\cr
1\textrm{ else]
\end{array]\right.$$

You can also use plain \LaTeX for equations
\begin{equation]\label{eq:fourier]
\hat f(\omega) = \int_{-\infty]^{\infty] f(x) e^{i\omega x] dx
\end{equation]
and refer to \autoref{eq:fourier] from text.
 -->

# Related research and software

To the best of our knowledge, we are unaware of any other software packages or tool providing a similar purpose or functionality for describing the logistic growth of the COVID-19 pandemic from a realistic finite two-dimensional perspective of natural growth.

This application of the `NLSIG` to modelling the COVID-19 pandemic was selected as the best paper at the *2nd African Symposium on Big Data, Analytics and Machine Intelligence and 6th TYAN International Thematic Workshop, December 3-4, 2020*.


# Acknowledgements

This work received no funding. 

# References
