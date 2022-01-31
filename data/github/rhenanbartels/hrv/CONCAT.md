[![Build Status](https://travis-ci.org/rhenanbartels/hrv.svg?branch=develop)](https://travis-ci.org/rhenanbartels/hrv)
[![codecov](https://codecov.io/gh/rhenanbartels/hrv/branch/develop/graph/badge.svg)](https://codecov.io/gh/rhenanbartels/hrv)
[![Current version at PyPI](https://img.shields.io/pypi/v/hrv.svg)](https://pypi.python.org/pypi/hrv)
[![Downloads per month on PyPI](https://img.shields.io/pypi/dm/hrv.svg)](https://pypi.python.org/pypi/hrv)
![Supported Python Versions](https://img.shields.io/pypi/pyversions/hrv.svg)
![Software status](https://img.shields.io/pypi/status/rows.svg)
[![License: LGPLv3](https://img.shields.io/pypi/l/hrv.svg)](https://github.com/rhenanbartels/hrv/blob/develop/LICENSE.md)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3960216.svg)](https://doi.org/10.5281/zenodo.3960216)


Pythonic package for Heart Rate Variability Analysis
===============================

version number: 0.2.10

author: Rhenan Bartels


![logo](https://raw.githubusercontent.com/rhenanbartels/hrv/develop/docs/logo-small.png)

Overview
--------

**hrv** is a simple Python module that brings the most widely used
techniques to extract information about cardiac autonomic functions through RRi series and Heart Rate Variability (HRV) analyses without losing the **Power** and **Flexibility**
of a native Python object.

In other words, the **hrv** module eases the manipulation, inspection, pre-processing, visualization, and analyses of HRV-related information. Additionally, it is written with idiomatic code and tries to implement the API of a built-in object, which might make it intuitive for Python users.


For a more in-depth explanation of how hrv works, please [Read the documentation](https://hrv.readthedocs.io/en/latest/)
---
title: 'HRV: a Pythonic package for Heart Rate Variability Analysis'
tags:
  - Python
  - cardiac autonomic system
  - heart rate variability 
authors:
  - name: Rhenan Bartels
    orcid: 0000-0002-6088-2186
    affiliation: "1"
  - name: Tiago Peçanha
    orcid: 0000-0003-4968-5525
    affiliation: "2"
affiliations:
  - name: Pulmonary Engineering Laboratory, Biomedical Engineering Program, Federal University of Rio de Janeiro, Rio de Janeiro, Brazil
    index: 1
  - name: Faculty of Medicine, University of São Paulo, Brazil
    index: 2
date: 25 October 2019
bibliography: paper.bib
---

# Summary
Heart rate variability (HRV) is a non-invasive tool to assess the cardiac autonomic integrity and cardiovascular homeostasis [@electrophysiology1996heart]. HRV quantifies the instantaneous variations in the RR intervals (RRi), which is produced by the balanced action of the parasympathetic and sympathetic branches of the autonomic nervous system (ANS) over the sinoatrial (SA) node, and modulated by different physiological inputs (e.g., respiration, blood pressure, temperature, emotions, etc) [@malik1990heart]. Specifically, the parasympathetic activation produces fast and short-lasting bradycardia, which results in high-frequency oscillations (i.e., 0.15 - 0.40 Hz) in heart rate. On the other hand, sympathetic activation produces slower and longer-lasting variations in heart rate, which results in low-frequency (i.e., ~ 0.1 Hz) oscillations in heart rate.  

Increased HRV indicates a predominance of the parasympathetic over the sympathetic activation in the SA node, and indicates enhanced cardiac autonomic flexibility and improved overall health [@malik1990heart]. Conversely, a reduced HRV is usually accompanied by sympathetic dominance and suggests increased rigidity and loss of control of the ANS to the cardiovascular system, which is a sign of disease vulnerability [@malik1990heart]. For instance, a reduced SDNN (i.e., standard deviation of RR intervals, a simple statistical-based time domain  index of HRV) has been shown to overperform traditional cardiovascular risk parameters in predicting mortality in a cohort of heart failure patients [@nolan1998prospective]. A reduced HRV has also been linked with metabolic dysfunction [@weissman2006power], increased inflammation [@sajadieh2004increased], depression [@sgoifo2015autonomic],  psychiatric disorders [@degiorgio2010rmssd], sleep disturbance [@burton2010reduced], among others. Based on this wide prognostic utility, the interest in approaches to evaluate HRV has shown an exponential growth in different medicine specialties and research fields in the recent years.  

HRV is routinely assessed using linear methods, through the calculation of different indices either in time- or frequency-domain. Time-domain consists of a collection of statistical metrics, such as the average value of RRi (mRRi), the standard deviation of RRi (SDNN; the NN stands for natural or sinusal intervals), the standard deviation of the successive differences (SDSD), the number or percentage of RRi longer than 50ms (NN50 and pNN50) and the root mean squared of successive difference in adjacent RRi (RMSSD - equation 1) [@electrophysiology1996heart].  Each of these indices quantifies different facets of the HRV, which are promoted by different autonomic sources. SDNN quantifies overall variability behind HRV, which is produced by both parasympathetic and sympathetic branches. NN50, pNN50, and RMSSD quantify beat-to-beat HRV, which is produced predominantly by the parasympathetic action in the heart. 

$$RMSSD = \sqrt{\frac{1}{N-1}\sum_{j-1}^{N-1}(RRi_{j} - RRi_{j+1})^2}$$ - Equation 1

where $N$ is the count of RRi values and $RRi_{j}$ is the jth RRi value.

The frequency-domain analysis quantifies the extent of contribution of each frequency component to the overall heart rate fluctuation (Figure 2). The main frequency components are the VLF (i.e., very low frequency; $<$ 0.04 Hz); LF (low frequency; 0.04-0.15 Hz) and the HF (i.e., high frequency; 0.15-0.40 Hz). The HF component is coupled with the respiratory fluctuation (i.e., respiratory sinus arrhythmia) and is produced by the parasympathetic modulation on the heart. The LF is mainly coupled with variations in the blood pressure (i.e.,  Mayer waves), and is thought to represent the modulation of both parasympathetic and sympathetic branches on the heart. The VLF does not have a defined physiological source, but it may involve alterations in heart rate produced by hormones and body temperature [@electrophysiology1996heart]. 

Roughly, frequency domain analysis involves the calculation of the spectral energy content of each frequency component through a power spectral density (PSD) estimation. Several methods have been developed to perform the PSD estimation and they are generally divided into two categories that provide comparable results: non-parametric and parametric methods, each with respective pros and cons [@electrophysiology1996heart]. The Welch periodogram [@welch1967use] is a non-parametric approach based on the Fourier Transform and consists of the average of several PSD estimations on different segments of the same RRi series, which is an important approach to reduce the spectral estimation variability [@welch1967use]. On the other hand, the autoregressive technique is the most widely used parametric method to estimate the spectral components of the HRV signal [@berntson1997heart]. The PSD estimation with the autoregressive method consists of a parametric representation of the RRi series and the frequency response of the estimated model. From the estimated PSD, generally, the following indices presented in Table 1 are calculated.

| Variable | Units | Frequency Band |
| -------- | ----- | -------------- |
| Total Power | ms<sup>2</sup> | 0 - 0.4 Hz |
| VLF | ms<sup>2</sup> | $<$ 0.04 Hz |
| LF | ms<sup>2</sup> | 0.04 - 0.15 Hz |
| HF | ms<sup>2</sup> | 0.15 - 0.4 Hz |
| LF/HF | | |
| LF<sub>n.u</sub> | normalized units $\frac{LF}{Total Power - VLF}$ | |
| HF<sub>n.u</sub> | normalized units $\frac{HF}{Total Power - VLF}$ | |

Non-linear indices are also frequently used to extract information from the ANS based on the heart rate fluctuations patterns. The Poincaré ellipse plot belongs to the non-linear methods and consists of a diagram in which each RRi is plotted as a function of the previous RRi value [@berntson1997heart]. In addition to the visual information about the RRi scatter given by the plot, two indices are extracted from this diagram: SD1 and SD2. The former reflects the short term fluctuations of the heart rate, and for this reason is highly correlated with the RMSSD, pNN50 and HF indices, while the latter reflects both short and long terms of the fluctuation of the heart rate and correlates with SDNN and LF indices. Additionaly, the SD1 index represents the standard deviation spread orthogonally to the identity line (y=x) and it is the ellipse width, whereas the SD2 index represents the standard deviation spread along the identity line and specifies the length of the ellipse. At the end of the following section, the Poincaré plot of a given RRi series is depicted using the module presented in the current article. 

The calculation of the SD1 and SD2 can be derived from SDSD and SDNN values as shown by equations 2 and 3 below:

$$SD1 = \sqrt{2SDNN^{2} - 2SD2^{2}}$$ - Equation 2

$$SD2 = \sqrt{2SDNN^{2} - \frac{1}{2}SDSD^{2}}$$ - Equation 3

There are several software packages written in many different programming languages that offer functions to work with RRi signals. Some of them have a command-line interface [@rodriguez2008rhrv] and others offer a user's interface to improve the interaction with the RRi series and the analyses  [@tarvainen2014kubios,@bartels2017sinuscor]. 
Specifically for Python, there is also open-source packages available and ready to work on HRV analisys, such as [hrvanalysis](https://github.com/Aura-healthcare/hrvanalysis), pyhrv [@Gomes2019], and heartpy [@van2019analysing]. Although these modules do a great work offering many of the most widely used techniques to deal with tachograms and to extract relevant information from HRV signals, their functions interface (API) relies on RRi signals stored as Python iterable or numpy arrays and is based mostly on the procedural programming paradigm. 

The `hrv` is a simple and open-source Python module that comes with the most common techniques for filtering, detrending and extracting information about the ANS from the RRi signals without losing the power and flexibility of a native Python object and a numpy arrays [@numpy]. It brings the necessary methods to work with a tachogram encapsulated in a Python class. In other words, once an RRi class is instantiated there are several methods available for visualization, descriptive statistics, slicing the signal in shorter segments, and displaying the metadata of the series.

With many software available to work with HRV analysis, the main reason why the `hrv` module is being developed is to improve and simplify the interaction with an RRi series with idiomatic Python code, closer to the native objects of this language. The object-oriented approach offered by the present module allows a strong relation between the RRi series and its methods, especially regarding time information. With a class representing and encapsulating the RRi object, each RRi value is bound to its respective time information, and therefore, after actions like slicing and filtering, the RRi series still keeps track of its information. Additionally, the instance's properties help to keep the state of the RRi series, informing, for instance, if it is already detrended and/or resampled.

The following sections present the basic workflow with an RRi series and gives a better overview of the functionalities available in the `hrv` module, starting with reading a file containing a tachogram, visualizing the given RRi series, dealing with noise filtering and detrending and, finally, calculating the time/frequency domain and non-linear HRV indices.

# Basic Usage

This section presents a non-exhaustive walkthrough of the features offered by the ```hrv``` module. To have access to the source code and more usage examples, please refer to the [software repository](https://github.com/rhenanbartels/hrv), the complete [documentation](https://hrv.readthedocs.io/en/latest/index.html) or the [notebooks](https://github.com/rhenanbartels/hrv/tree/master/notebooks) with some use cases of the ```hrv``` module.

Once the RRi series is created in Python using the `hrv.io` submodule, which supports text, CSV and hrm (Polar<sup>TM</sup>) files, or from any Python iterable (i.e lists, tuples, etc), an RRi instance with the necessary methods to implement the Python iterable pattern is created. With the RRi object it is possible to iterate (i.e `[r for r in rri_series]`), search for a value at a given index (i.e `rri_series[0]`), and slice the tachogram (i.e `rri_series[5:10]`). As the RRi class also implements some of the behaviors of the numpy array [@numpy], it is possible to perform math operations with the tachogram, i.e: `rri_series / 1000`.

The RRi class also has methods for basic statistical metrics calculation, such as average, standard deviation, min and max, and others. In order to access a complete Python dictionary containing all available statistical metrics of an RRi instance, it is possible to call the `describe()` method. Features for visualization are also present in the RRi class. In order to visualize the time series represented by the RRi series, the ``plot()`` method can be called. The visualization of the histograms showing the distributions of RRi or heart rate time series is also possible with the method `hist()`.


## Read a file containing RRi values and visualising it 

The following code snippet shows how to read a RRi series from a single column CSV file and plot the respective series with black lines.

```python
from hrv.io import read_from_csv  

rri = read_from_csv('path/to/file.csv')

fig, ax = rri.plot(color='k')

```

![RR intervals of a young subject at rest condition produced with the `plot()` method from the RRi class.](rri_series.png)

To retrieve statistical properties of a RRi series the method `describe()` can be invoked:

```python
desc = rri.describe()

desc
----------------------------------------
                   rri          hr
----------------------------------------
min             750.00       66.30
max             905.00       80.00
mean            805.50       74.78
var            2646.25       20.85
std              51.44        4.57
median          805.00       74.54
amplitude       155.00       13.70


print(desc['std'])
{'rri': 51.44171459039833, 'hr': 4.5662272355549725}
```

## Pre-processing
### Filtering the RRi series

In some cases and for many different reasons, the tachogram may present with movement artifacts or undesired RRi values, which may jeopardize the analysis results. One way to deal with this scenario is to apply filters to the RRi series. For this reason, the `hrv` package offers four lowpass filters for noise removal: moving average, which given an order value $N$, replaces every RRi value by the average of its $N$ neighbors values;  the moving median, which works similarly to the moving average filter, but apply the median function; the quotient filter [@piskorski2005filtering], that removes the RRi values which the ratio with its adjacent RRi is greater than 1.2 or smaller than 0.8; and finally, the threshold filter, which is inspired in Kubios [@tarvainen2014kubios] threshold-based artifact correction algorithm: each RRi is compared to a local value consisting of the median of adjacent RRi. If the difference between a given RRi and the local median is greater than the threshold in milliseconds this RRi is considered an ectopic beat. Ectopic RRi values are replaced with cubic spline interpolation of the entire tachogram.

```python
from hrv.filters import moving_median, quotient

filt_rri_median = moving_median(rri, order=3)
filt_rri_quotient = quotient(rri)

filt_rri_median.plot(ax=ax)
filt_rri_quotient.plot(ax=ax)
```

![The left panel shows the original RRi (blue line) and after filtering with a moving median filter with order equal to 3 (orange line). The left panel depicts the original RRi (blue line) and after filtering with a quotient filter (orange line). This picture was created using the `plot()` method from the `RRi` instance.](rri_filtered.png)

### Detrending the RRi series

Although the very-low-frequency components of the PSD function might have useful information, they are generally removed from the RRi signals before the frequency-domain analysis is performed. This pre-processing step before the frequency-domain analysis is important to remove intrinsic slow trends that are present in the HR fluctuation. This non-stationary behavior may contaminate the overall dynamic of the RRi series and influence the results, especially the VLF and LF measures [@tarvainen2002advanced]. For this reason, several methods have been developed to extract the frequency components responsible for the non-stationary behavior of the RRi series. 

Among the methods available in the literature for detrending the RRi series, the `hrv` module offers the polynomial detrend, which consists of the subtraction of a $Nth$ degree polynomial from the RRi signal, where $N$ is smaller than the length of the tachogram. It also offers the Smoothness Priors method [@tarvainen2002advanced], which is widely used in HRV analyses and acts as a lowpass filter to remove complex trends from the RRi series. Finally, the `hrv` module also offers a detrending method that uses the Savitsky-Golay lowpass filter to remove low-frequency trends from the RRi series. The following code fragment applies the polynomial detrend with a degree equal to 1 and the Savitsky-Golay filter to remove the slow frequency components from an RRi recorded during rest.

```python
from hrv.detrend import polynomial_detrend, sg_detrend
from hrv.sampledata import load_rest_rri

rri = load_rest_rri()

detrended_rri_poly = polynomial_detrend(rri, degree=1)
detrended_rri_sg = sg_detrend(rri, window_length=51, polyorder=3)

detrended_rri_poly.plot()
detrended_rri_sg.plot()
```

![The left panel shows the original RRi (blue line) and after detrending with polynomial function with degree equal to 1 (black line). The left panel depicts the original RRi (blue line) and after detreding with a Savitsky-Golay lowpass filter (black line).](rri_detrended.png)


## Analyses
### Time Domain 

In order to calculate the time-domain indices, the function `time_domain` can be imported from the submodule `hrv.classical` and applied to any Python iterable containing the RRi series including the RRi instance from the module presented in this article.

```python
from hrv.classical import time_domain

results = time_domain(rri)
print(results)

{'mhr': 66.528130159638053,
 'mrri': 912.50302419354841,
 'nn50': 337,
 'pnn50': 33.971774193548384,
 'rmssd': 72.849900286450023,
 'sdnn': 96.990569261440797
 'sdsd': 46.233829821038042}
```

### Frequency Domain

Similarly to the `time_domain` function, to calculate the frequency-domain indices, the `frequecy_domain`, which is also placed in the `hrv.classical` submodule, can be used. The `frequency_domain` function present in the `hrv` module takes care of the pre-processing steps:  the detrending of the RRi series (which the default is a linear function, but can be any custom Python function), interpolation using cubic splines (also accepts linear interpolation) and resampling at a given frequency, the default is 4Hz.

When Welch’s method is selected, a window function (default: hanning), the number of RRi values per segment and the length of superposition between adjacent segments can be chosen. When the AR method is selected, the order of the model (default 16) can be set.

The area under the curve of each frequency range in the estimated PSD is calculated using the trapezoidal method. As a default, the hrv module uses the frequencies cutoffs shown in Table 1 to limit the integration range of each frequency domain indices, however, it is possible to set the frequency range of VLF, LF, and HF  in the `frequency_domain` function call.

```python
from hrv.classical import frequency_domain

results = frequency_domain(
    rri=rri,
    fs=4.0,
    method='welch',
    interp_method='cubic',
    detrend='linear'
)
print(results)

{'hf': 1874.6342520920668,
 'hfnu': 27.692517001462079,
 'lf': 4894.8271587038234,
 'lf_hf': 2.6110838171452708,
 'lfnu': 72.307482998537921,
 'total_power': 7396.0879278950533,
 'vlf': 626.62651709916258}
```

![Power Spectral Density of a RRi series estimated with the Welch's method.](psd.png)

### Non-linear
Finally, among the non linear metrics, `hrv` module offers SD1 and SD2, which can be calculated with the `non_linear` function from the `hrv.classical` submodule.

```python
from hrv.classical import non_linear

results = non_linear(rri)
print(results)

{'sd1': 51.538501037146382,
 'sd2': 127.11460955437322}
```

The respective Poincaré plot of a given RRi series can be depicted with the `poincare_plot()` method, as follows:

```python
rri.poincare_plot()
```

![Poincaré plot of a given RRi series.](poincare.png)

# Dependencies

The `hrv` package depends on the following modules: numpy [@numpy], matplotlib [@hunter2007matplotlib], scipy [@scipy] and spectrum [@cokelaer2017spectrum].

# Acknowledgements

The authors are grateful to CAPES (Coordenadoria de Aperfeiçoamento de Pessoal de Nível Superior), FAPERJ (Fundação de Amparo à Pesquisa do Estado do Rio de Janeiro), FAPESP (2016/23319-0 - Fundação de Amparo à Pesquisa do Estado do São Paulo) and CNPq (Conselho Nacional de Desenvolvimento Científico e Tecnológico) for financial support. We are also grateful to Wilson Mello a.k.a Bakudas for creating the nice logo of the `hrv` module.

# References
First Steps
===========

Installation
############

To install use pip:

.. code-block:: bash

   pip install hrv

Or clone the repo:

.. code-block:: bash

   git clone https://github.com/rhenanbartels/hrv.git
   python setup.py install


Basic Usage
######################

**Create an RRi instance**

Once you create an RRi object you will have the power of a native Python iterable object.
This means, that you can loop through it using a **for loop**, get a just a part of the series using native
**slicing** and much more. Let us try it:


.. code-block:: python

    from hrv.rri import RRi

    rri_list = [800, 810, 815, 750, 753, 905]
    rri = RRi(rri_list)

    print(rri)
    RRi array([800., 810., 815., 750., 753., 905.])

**Slicing**

.. code-block:: python

    print(rri[0])
    800.0

    print(type(rri[0]))
    numpy.float64

    print(rri[::2])
    RRi array([800., 815., 753.])

**Logical Indexing**

.. code-block:: python

    from hrv.rri import RRi

    rri = RRi([800, 810, 815, 750, 753, 905])
    rri_ge = rri[rri >= 800]

    rri_ge
    RRi array([800., 810., 815., 905.])

**Loop**

.. code-block:: python

    for rri_value in rri:
        print(rri_value)

    800.0
    810.0
    815.0
    750.0
    753.0
    905.0

**Note:**
When time information is not provided, time array will be created using the cumulative sum of successive RRi. After cumulative sum, the time array is subtracted from the value at `t[0]` to make it start from `0s`

**RRi object and time information**

.. code-block:: python

    from hrv.rri import RRi

    rri_list = [800, 810, 815, 750, 753, 905]
    rri = RRi(rri_list)

    print(rri.time)
    array([0.   , 0.81 , 1.625, 2.375, 3.128, 4.033]) # Cumsum of rri values minus t[0]

    rri = RRi(rri_list, time=[0, 1, 2, 3, 4, 5])
    print(rri.time)
    [0. 1. 2. 3. 4. 5.]

**Note:**
Some validations are made in the time list/array provided to the RRi class, for instance: 

 * RRi and time list/array must have the same length;
 * Time list/array can not have negative values;
 * Time list/array must be monotonic increasing.

**Basic math operations**

With RRi objects you can make math operatins just like a numpy array:

.. code-block:: python

    rri
    RRi array([800., 810., 815., 750., 753., 905.])

    rri * 10
    RRi array([8000., 8100., 8150., 7500., 7530., 9050.])

    rri + 200
    RRi array([1000., 1010., 1015.,  950.,  953., 1105.])

**Works with Numpy functions**

.. code-block:: python

    import numpy as np

    rri = RRi([800, 810, 815, 750, 753, 905])

    sum_rri = np.sum(rri)
    print(sum_rri)
    4833.0

    mean_rri = np.mean(rri)
    print(mean_rri)
    805.5

    std_rri = np.std(rri)
    print(std_rri)
    51.44171459039833
RRi Visualization
=================

The RRi class brings a very easy way to visualize your series:

Plot RRi Series
###############

.. code-block:: python

    from hrv.io import read_from_text

    rri = read_from_text('path/to/file.txt')
    fig, ax = rri.plot(color='k')

.. image:: ../figures/rri_fig.png
   :width: 500 px

RRi histogram and Heart Rate Histogram
######################################

.. code-block:: python

    rri.hist()

    rri.hist(hr=True)

.. image:: ../figures/rri_hist.png
   :width: 500 px

.. image:: ../figures/hr_hist.png
   :width: 500 px
Pre-Processing
==============

Filters
#######

**Moving Average**

.. code-block:: python

    from hrv.filters import moving_average
    filt_rri = moving_average(rri, order=3)

    fig, ax = rri.plot()
    filt_rri.plot(ax=ax)

.. image:: ../figures/mov_avg.png
    :width: 500 px

**Moving Median**

.. code-block:: python

    from hrv.filters import moving_median
    filt_rri = moving_median(rri, order=3)

    fig, ax = rri.plot()
    filt_rri.plot(ax=ax)


.. image:: ../figures/mov_median.png
    :width: 500 px

**Quotient**

`Read more`_

.. _Read more: https://www.ncbi.nlm.nih.gov/pubmed/17322593

.. code-block:: python

    from hrv.filters import quotient
    filt_rri = quotient(rri)

    fig, ax = rri.plot()
    filt_rri.plot(ax=ax)

.. image:: ../figures/quotient.png
    :width: 500 px

**Threshold Filter**

This filter is inspired by the threshold-based artifact correction algorithm offered by kubios_ <sup>&reg;</sup> .
To elect outliers in the tachogram series, each RRi is compared to the median value of local RRi (default N=5).
All the RRi which the difference is greater than the local median value plus a threshold is replaced by
cubic_ interpolated RRi.

.. _kubios: https://www.kubios.com/
.. _cubic: https://en.wikiversity.org/wiki/Cubic_Spline_Interpolation

The threshold filter has five pre-defined strength values:

    * Very Low: 450ms
    * Low: 350ms
    * Medium: 250ms
    * Strong: 150ms
    * Very Strong: 50ms

It also accepts custom threshold values (in milliseconds).
The following snippet shows the ectopic RRi removal:

.. code-block:: python

    from hrv.filters import threshold_filter
    filt_rri = threshold_filter(rri, threshold='medium', local_median_size=5)

    fig, ax = rri.plot()
    filt_rri.plot(ax=ax)

.. image:: ../figures/threshold_filter.png
    :width: 500 px

Detrending
##########

The **hrv** module also offers functions to remove the non-stationary trends from the RRi series.
It allows the removal of slow linear or more complex trends using the following approaches:

**Polynomial models**

Given a degree a polynomial filter is applied to the RRi series and subtracted from the tachogram

.. code-block:: python

    from hrv.detrend import polynomial_detrend

    rri_detrended = polynomial_detrend(rri, degree=1)

    fig, ax = rri.plot()
    rri_detrended.plot(ax, color='k')

.. image:: ../figures/polynomial_detrend.png
    :width: 500 px


**Smoothness priors**

Developed by Tarvainen *et al*, allow the removal of complex trends. Visit here_ for more information.
It worth noticing that the detrended RRi with the Smoothness priors approach is also interpolated
and resampled using frequency equals to ```fs```.

.. _here: https://ieeexplore.ieee.org/document/979357

.. code-block:: python

    from hrv.detrend import smoothness_priors

    rri_detrended = smoothness_priors(rri, l=500, fs=4.0)

    fig, ax = rri.plot()
    rri_detrended.plot(ax, color='k')

.. image:: ../figures/smoothness_priors.png
    :width: 500 px

**Note:**
this approach depends on a numpy matrix inversion and due to floating-point precision it might
present round-off errors in the trend calculation

**Savitzky-Golay**

Uses the lowpass filter known as  Savitzky-Golay filter to smooth the RRi series and remove slow components from the tachogram

.. code-block:: python

    from hrv.detrend import sg_detrend

    rri_detrended = sg_detrend(rri, window_size=51, polyorder=3)

    fig, ax = rri.plot()
    rri_detrended.plot(ax, color='k')

.. image:: ../figures/savitzky_golay_detrend.png
    :width: 500 px
RRi statistics
==============

Basic Statistics
################
The RRi object implements some basic statistics information about its values:

* mean
* median
* standard deviation
* variance
* minimum
* maximum
* amplitude

Some examples:

.. code-block:: python

    from hrv.rri import RRi

    rri = RRi([800, 810, 815, 750, 753, 905])

    # mean
    rri.mean()
    805.5

    # median
    rri.median()
    805.0

You can also have a complete overview of its statistical charactheristic

.. code-block:: python

    desc = rri.describe()
    desc

    ----------------------------------------
                       rri          hr
    ----------------------------------------
    min             750.00       66.30
    max             905.00       80.00
    mean            805.50       74.78
    var            2646.25       20.85
    std              51.44        4.57
    median          805.00       74.54
    amplitude       155.00       13.70

    print(desc['std'])
    {'rri': 51.44171459039833, 'hr': 4.5662272355549725}

RRi Basic Information
#####################

.. code-block:: python

    rri = RRi([800, 810, 815, 750, 753, 905])
    rri.info()

    N Points: 6
    Duration: 4.03s
    Interpolated: False
    Detrended: False
    Memory Usage: 0.05Kb


Analysis
========

Time Domain Analysis
####################

.. code-block:: python

    from hrv.classical import time_domain
    from hrv.io import read_from_text

    rri = read_from_text('path/to/file.txt')
    results = time_domain(rri)
    print(results)

    {'mhr': 66.528130159638053,
     'mrri': 912.50302419354841,
     'nn50': 337,
     'pnn50': 33.971774193548384,
     'rmssd': 72.849900286450023,
     'sdnn': 96.990569261440797,
     'sdsd': 46.233829821038042}

Frequency Domain Analysis
#########################
.. code-block:: python

    from hrv.classical import frequency_domain
    from hrv.io import read_from_text

    rri = read_from_text('path/to/file.txt')
    results = frequency_domain(
        rri=rri,
        fs=4.0,
        method='welch',
        interp_method='cubic',
        detrend='linear'
    )
    print(results)

    {'hf': 1874.6342520920668,
     'hfnu': 27.692517001462079,
     'lf': 4894.8271587038234,
     'lf_hf': 2.6110838171452708,
     'lfnu': 72.307482998537921,
     'total_power': 7396.0879278950533,
     'vlf': 626.62651709916258}

Non-linear Analysis
###################

.. code-block:: python

    from hrv.classical import non_linear
    from hrv.io import read_from_text

    rri = read_from_text('path/to/file.txt')
    results = non_linear(rri)
    print(results)

    {'sd1': 51.538501037146382,
     'sd2': 127.11460955437322}

It is also possible to depict the Poincaré Plot, from which SD1 and SD2 are derived:

.. code-block:: python

    rri.poincare_plot()

.. image:: ../figures/poincare.png
   :width: 500 px
Read RRi files
==============

Read .txt files
####################

Text files contains a single column with all RRi values.
Example of RRi text file

.. code-block:: bash

    800
    810
    815
    750

.. code-block:: python

    from hrv.io import read_from_text

    rri = read_from_text('path/to/file.txt')

    print(rri)
    RRi array([800., 810., 815., 750.])

Read Polar® (.hrm) files
########################################

The .hrm files contain the RRi acquired with Polar®

A complete guide for .hrm files can be found here_

.. _here: https://www.polar.com/files/Polar_HRM_file%20format.pdf

.. code-block:: python

    from hrv.io import read_from_hrm

    rri = read_from_hrm('path/to/file.hrm')

    print(rri)
    RRi array([800., 810., 815., 750.])

Read .csv files
#####################################
Example of csv file:

.. code-block:: bash

    800,
    810,
    815,
    750,

.. code-block:: python

    from hrv.io import read_from_csv

    rri = read_from_csv('path/to/file.csv')

    print(rri)
    RRi array([800., 810., 815., 750.])

**Note**:
When using **read_from_csv** you can also provide a column containing time information. Let's check it.

.. code-block:: bash

    800,1
    810,2
    815,3
    750,4

.. code-block:: python

    rri = read_from_csv('path/to/file.csv', time_col_index=1)

    print(rri)
    RRi array([800., 810., 815., 750.])

    print(rri.time)
    array([0., 1., 2., 3., 4.])

RRi Sample Data
###############

The hrv module comes with some sample data. It contains:

* RRi collected during rest
* RRi collected during exercise
* RRi containing ectopic beats

**Rest RRi**

.. code-block:: python

    from hrv.sampledata import load_rest_rri

    rri = load_rest_rri()
    rri.plot()


.. image:: ../figures/rest_rri.png
    :width: 500 px

**Exercise RRi**

.. code-block:: python

    from hrv.sampledata import load_exercise_rri

    rri = load_exercise_rri()
    rri.plot()

.. image:: ../figures/exercise_rri.png
    :width: 500 px

**Noisy RRi**

.. code-block:: python

    from hrv.sampledata import load_noisy_rri

    rri = load_noisy_rri()
    rri.plot()

.. image:: ../figures/noisy_rri.png
    :width: 500 px
.. hrv documentation master file, created by
   sphinx-quickstart on Wed Dec 18 12:07:31 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to hrv's documentation!
===============================

**hrv** is a simple Python module that brings the most widely used
techniques to extract information about cardiac autonomic functions through RRi series and Heart Rate Variability (HRV) analyses without losing the **Power** and **Flexibility**
of a native Python object.

In other words, the **hrv** module eases the manipulation, inspection, pre-processing, visualization, and analyses of HRV-related information. Additionally, it is written with idiomatic code and tries to implement the API of a built-in object, which might make it intuitive for Python users.


.. _numpy: http://numpy.org


Contents
******************

.. toctree::
   :maxdepth: 3

   firststeps
   readfiles
   statistics
   visualization
   timeslicing
   preprocessing
   analysis
   nonstationary
   contribution
Non-stationary RRi series
============

In some situations like physical exercise, the RRi series might present
a non-stationary behavior. In cases like these, classical approaches are not
recommended once the statistical properties of the signal vary over time.

The following figure depicts the RRi series recorded on a subject riding a bicycle.
Without running analysis and only visually inspecting the time series, is possible
to tell that the average and the standard deviation of the RRi are not constant
as a function of time.

.. image:: ../figures/exercise_rri.png
    :width: 500 px

In order to extract useful information about the dynamics of non-stationary RRi series,
the following methods applies the classical metrics in shorter running adjacent segments, 
as illustrated in the following image:

.. image:: ../figures/sliding_segments.png
    :width: 500 px

For example, for a segment size of **30s** (S) and **15s** (O) overlap a signal with **300s** (D) will have P segments:

P = int((D - S) / (S - O)) + 1

P = int((300 - 30) / (30 - 15)) + 1 = 19 segments

Time Varying
############

Time domain indices applied to shorter segments

.. code-block:: python

    from hrv.sampledata import load_exercise_rri
    from hrv.nonstationary import time_varying

    rri = load_exercise_rri()
    results = time_varying(rri, seg_size=30, overlap=0)
    results.plot(index="rmssd", marker="o", color="r")


.. image:: ../figures/tv_rmssd.png
    :width: 500 px


Plot the results from **time varying** together with its respective RRi series


.. code-block:: python

    from hrv.sampledata import load_exercise_rri
    from hrv.nonstationary import time_varying

    rri = load_exercise_rri()
    results = time_varying(rri, seg_size=30, overlap=0)
    results.plot_together(index="rmssd", marker="o", color="k")

.. image:: ../figures/tv_together.png
    :width: 500 px

Short Time Fourier Transform
############################

To be implemented.
Time Slicing
============

It is also possible to slice RRi series with time range information
(in **seconds**).

In the following example, we are taking a slice that starts at `100s ` and ends at `200s`.

.. code-block:: python

    from hrv.io import read_from_text

    rri = read_from_text('path/to/file.txt')
    rri_range = rri.time_range(start=100, end=200)

    fig, ax = rri_range.plot(marker='.')

.. image:: ../figures/rri_range.png
    :width: 500 px

Time offset can be reset from the RRi series range:

.. code-block:: python

    rri_range.reset_time(inplace=True)

.. image:: ../figures/rri_range_reset.png
    :width: 500 px
Contribution start guide
========================

The preferred way to start contributing for the project is creating a virtualenv (you can do by using virtualenv,
virtualenvwrapper, pyenv or whatever tool you'd like). Only **Python 3.x** are supported


Preparing the eviroment
#######################

Create the virtualenv:

.. code-block:: bash

    mkvirtualenv hrv

Install all dependencies:

.. code-block:: bash

    pip install -r requirements.txt

Install development dependencies:

.. code-block:: bash

    pip install -r dev-requirements.txt

Running the tests
#################

In order to run the tests, activate the virtualenv and execute pytest:

.. code-block:: bash

    workon <virtualenv>
    pytest -v
    # or
    make test


Coding and Docstring styles
##########################


Generally, we try to use Python common styles conventions as described
in `PEP 8`_ and `PEP 257`_, which are also followed by the `numpy`_ project.

.. _PEP 8: https://www.python.org/dev/peps/pep-0008/ 
.. _PEP 257: https://www.python.org/dev/peps/pep-0257/
.. _numpy: https://numpydoc.readthedocs.io/en/latest/format.html

Example
********
.. code-block:: python

    def moving_average(rri, order=3):
        """
        Low-pass filter. Replace each RRi value by the average of its ⌊N/2⌋
        neighbors. The first and the last ⌊N/2⌋ RRi values are not filtered

        Parameters
        ----------
        rri : array_like
        sequence containing the RRi series
        order : int, optional
        Strength of the filter. Number of adjacent RRi values used to calculate
        the average value to replace the current RRi. Defaults to 3.

        .. math::
        considering movinge average of order equal to 3:
        RRi[j] = sum(RRi[j-2] + RRi[j-1] + RRi[j+1] + RRi[j+2]) / 3

        Returns
        -------
        results : RRi array
        instance of the RRi class containing the filtered RRi values

        See Also
        -------
        moving_median, threshold_filter, quotient

        Examples
        --------
        >>> from hrv.filters import moving_average
        >>> from hrv.sampledata import load_noisy_rri
        >>> noisy_rri = load_noisy_rri()
        >>> moving_average(noisy_rri)
        RRi array([904., 918., 941.66666667, ..., 732.66666667, 772.33333, 808.])
        """

We also encourage the  use of code linters, such `isort`_ , `black`_ and `autoflake`_.

.. _isort: https://github.com/timothycrosley/isort
.. _black: https://github.com/psf/black
.. _autoflake: https://github.com/myint/autoflake


.. code-block:: bash

    autoflake --in-place --recursive --remove-unused-variables --remove-all-unused-imports .
    sort -rc .
    black .
    
