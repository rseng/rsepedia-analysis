---
title: '```beamshapes```: a Python package to generate directivity patterns for various sound source models'
authors:
- affiliation: 1, 2, 3
  name: Thejasvi Beleyur
  orcid: 0000-0001-5360-4383
date: "May 2021"
output: pdf_document
bibliography: paper.bib
tags:
- Python
- acoustics
- bioacoustics
affiliations:
- index: 1
  name: Centre for the Advanced Study of Collective Behaviour, University of Konstanz,
    Konstanz
- index: 2
  name: Acoustic and Functional Ecology, Max Planck Institute for Ornithology, Seewiesen
- index: 3
  name: Max Planck Institute for Animal Behavior, Radolfzell
---

# Summary

Sound sources such as human beings or loudspeakers often exhibit a 'directionality' in how loud they sound at different angles.
A listener or microphone placed at various angles at a fixed radius will often pick up sometimes drastically different sound levels. The same sound source may however sometimes produce omnidirectional sound fields. The directionality of a sound source can be modelled as a combination of the frequency of the emitted
sound and the geometry of the vibrating and non-vibrating parts of the sound source itself. The resulting pattern of sound radiation with angular location is called the *directivity* of the source [@beranek2012acoustics]. 

Directivity ($D$) describes the relative sound levels at a given angle $D(\theta)$
with relation to the level on-axis ($D(0)$, and thus $directivity = \frac{D_{\theta}}{D_{0}}$). Directivity
functions exist for a wide variety of sound sources that can be modelled analytically. A well-known example 
of a directivity function is of the piston in an infinite baffle. The piston is a circular surface of radius
$a$, vibrating back and forth about a hole in an infinite wall or baffle. The directivity is described by $\frac{2J_{1}(ka \times sin \theta)}{ka \times sin \theta}$, where $J_{1}$ is the Bessel's function of the first kind, $k$ is the wavenumber, where $k = \frac{2\pi f}{c}$, where $f$ is the frequency of the sound, and $c$ is the speed of sound. The angle of the receiver is $\theta$, which varies from 0-2$\pi$ radians in the azimuth. The ```beamshapes``` package aims to provide an easy interface to generate direcitivities of various sound sources.

Directivity functions can be used in a two-fold manner: 1) to deliberately engineer devices to suit particular specifications, e.g., loudspeaker sound fields [@beranek2012acoustics] and 2) to infer parameters of a sound source itself having assumed a relevant model, e.g., estimating the direction of call emission [@guarato2011method] and mouth aperture of bat echolocation calls [@jakobsen2013convergent;@kounitsky2015bats].


![Directivity patterns ($\frac{D_{\theta}}{D_{0}}$) of the currently implemented sound sources for a common set of $ka$ values for comparison, where $k$ is the wavenumber ($\frac{2\pi\:f}{c}$) and $a$ is the piston radius - or its equivalent. The directivity pattern shows the ratio between the off-axis sound level at ($\theta^{\circ}$) to the on-axis level at ($0^{\circ}$) in decibels. A) piston in an infinite baffle B) piston in a sphere with the half-angle aperture $\alpha = 30^{\circ}$, and where the piston radius $a=R\:sin \:\alpha$ C) oscillating cap of a sphere with  $\alpha=30^{\circ}$, and the equivalent piston radius is $a=R\:sin \:\alpha$ D) vibrating point on a sphere, here the $kR=ka$ for comparison with the other models. \label{fig:pistonsphereinf}](paper_related/piston_sphere_baffle.png)


# Statement of need

A host of published directivity functions exist in the literature, but it is my experience that their computational implementations remain as in-house scripts that are often in proprietary language platforms. To my knowledge, there are no computational implementations of multi-parameter sound source directivities. For instance, the ```levitate``` package [@levitate] implements directivities of 'simple' sound sources that can be defined with one parameter (e.g., piston in a baffle, circular ring). Since analytical solutions are available, the directivities of simple source models can be rapidly implemented and computed. In contrast to simple source models, the directivities of other models (as in this package) involve more parameters (e.g., piston in a sphere, oscillating cap of a sphere) and do not have analytical solutions. Their directivity calculations require numerical routines and run-time optimisation. The advantage of more involved source models is however their ability to capture aspects of experimental sound sources. In this paper, I present ```beamshapes```, a Python package that currently implements directivity patterns for two 'simple' and two 'involved' sound source models. As of this publication, ```beamshapes``` version 0.2.0 implements four sound sources: 1) piston in an infinite baffle, 2) point source on a sphere 3) oscillating cap of a sphere, and 4) piston in a sphere. 

Computational implementations of directivity functions often require long run-times due to the intensive numerical routines and arbitrary precision math involved. Long run-times hinder scientific projects in reducing the number of models and parameter space that can be explored. ```beamshapes``` boasts parallelised code to generate significant speed-ups in run-times. 

The availability of openly-available directivity functions will hopefully stir the acoustics, and specifically the bio-acoustics community to rigorously test and compare sound radiation with models. The availability of multiple implementations allows comparison of data with multiple models. Until recently, computational power has perhaps been a limiting factor in calculating directivity functions. Models with easily calculable outputs have thus been favoured (e.g., piston in an infinite baffle), especially in the field of bioacoustics [@strother1970acoustical;@mogensen1979sound]. 


Despite the recent availability of computational power, older, simpler models with limited biological relevance continue to dominate the field of bio-acoustics. For instance, the piston in an infinite baffle only predicts the beam-shape for a $\pm90^{\circ}$ range off-axis, and assumes front-back symmetry (\autoref{fig:pistonsphereinf} A). This is unrealistic for most vocalising animals, especially echolocating bats and odontocetes. However, the piston in an infinite baffle continues to be a standard reference model for multiple studies ranging over the past decade [@jakobsen2013convergent;@kounitsky2015bats;@macaulay2020high]. The piston in a sphere, for instance, recreates many of the highly directional central and side lobes and that is seen in bats, while also predicting sound radiation behind the source (\autoref{fig:pistonsphereinf} B). In contrast to the echolocation literature, the bird [@witkin1977importance;@larsen1990directionality;@brumm2002sound;@patricelli2007differences;@patricelli2008acoustic;@yorzinski2010birds]  and frog call literature [@gerhardt1975sound;@rodriguez2020sound] for instance has been dominated by quantitative characterisations of sound radiation with no attempts at directly comparing measurements to model predictions. Using model-based directivity patterns to infer source parameters allows the discovery of common parameter spaces that birds occupy and facilitates cross-species comparisons. The oscillating cap of a sphere and vibrating point on a sphere (\autoref{fig:pistonsphereinf} C,D) are two other models of potential relevance to bioacousticians attempting to describe the sound radiation of bird and frog calls for instance. 

Future releases of ```beamshapes``` are scheduled to include directivity patterns for additional models of interest such as rectangular cap of a sphere or piston in a closed finite baffle.


# Software packages used in this work
```beamshapes``` relies on the Python open-source ecosystem and is built on the ```numpy``` [@2020NumPy], ```scipy``` [@2020SciPy], ```sympy``` [@meurer2017sympy], ```mpmath``` [@mpmath] and ```flint``` [@hart2011flint] libraries. 

# Acknowledgements
TB thanks Gaurav Dhariwal for his continual math advice and inputs, and Tim Mellow for 
strong support through sharing his Mathematica code and timely clarifications. TB thanks Lasse Jakobsen and Holger Goerlitz for productive discussions leading to this project. This work was executed through a combination of TB's private time, a DAAD stipend, the IMPRS for Organismal Biology,  and then finally by a CASCB Medium grant.

# References
[![Build Status](https://app.travis-ci.com/thejasvibr/bat_beamshapes.svg?branch=dev)](https://app.travis-ci.com/thejasvibr/bat_beamshapes)
[![status](https://joss.theoj.org/papers/820fd37f8255a8c533d6cc4c9475ecb5/status.svg)](https://joss.theoj.org/papers/820fd37f8255a8c533d6cc4c9475ecb5)
[![Documentation Status](https://readthedocs.org/projects/beamshapes/badge/?version=latest)](https://beamshapes.readthedocs.io/en/latest/?badge=latest)


# Beamshapes

Calculate directivity of various sound source models and get their beamshapes.

This package is released under an MIT license. 

## Getting started with *beamshapes*

```
>>> import matplotlib.pyplot as plt 
>>> import numpy as np 
>>> import beamshapes
>>> from beamshapes import piston_in_infinite_baffle_directivity as PIB # short alias
>>> input_parameters = {'k':50, 'a':0.1}
>>> angles = np.linspace(-np.pi/2,np.pi/2,50)
>>> _, directionality = PIB(angles, input_parameters) # output the dB(D(theta)/D(on-axis))
>>> plt.figure();a0 = plt.subplot(111, projection='polar')
>>> plt.plot(angles, directionality) # plot the beamshape !!
```

For more detailed use-cases, check out the [example gallery online](https://beamshapes.readthedocs.io/en/latest/gallery_examples/index.html)!

## Installation 

*PyPi installation (>=version 0.2.1)*

```pip install beamshapes```

*Local installation instructions*
For the steps below to work you need to have a working Python installation that you can access from the command line. It is recommended to do the installation in a  [virtual environment](https://realpython.com/effective-python-environment/#virtual-environments). 

1. Clone the Github repository ```git clone https://github.com/thejasvibr/bat_beamshapes.git```
1. Change directories to the downloaded repo, and switch to the *dev* branch: ```git checkout dev``` 
1. Install the dependencies with ```pip install -r beamshapes/tests/requirements_test.txt```
1. Install *beamshapes* with ```python setup.py install```


## Detailed documentation 
For more details on the concepts and source documentation - please check out the [online docs](https://beamshapes.rtfd.io).

## Citation information 
If you use this package - please cite the paper: 

APA-style format

*Beleyur, T., (2022). beamshapes: a Python package to generate directivity patterns for various sound source models. Journal of Open Source Software, 7(69), 3740, https://doi.org/10.21105/joss.03740*


Bibtext format: 

```
@article{Beleyur2022,
  doi = {10.21105/joss.03740},
  url = {https://doi.org/10.21105/joss.03740},
  year = {2022},
  publisher = {The Open Journal},
  volume = {7},
  number = {69},
  pages = {3740},
  author = {Thejasvi Beleyur},
  title = {```beamshapes```: a Python package to generate directivity patterns for various sound source models},
  journal = {Journal of Open Source Software}
}
```
If possible, and space allows also do mention the package version (e.g. beamshapes 0.2.X) .
You can access the version number of the package being used with 
```
>>> import beamshapes
>>> print(beamshapes.__version__)
```



## Future implementations
* Piston on a cylinder
* Rectangular piston on a prolate spheroid ([paper](https://asa.scitation.org/doi/pdf/10.1121/1.1778840?casa_token=wDAHTxJBISUAAAAA:MW-OSeGIkft-mces_mJgFBuyOhzI1qpPbc_7Xuu9EhDDD8CF8vnCIYaGyVivUb2qOpFda4GkPWto))
#  Contributing 

## Be nice, patient and respectful
Like it says above, while interacting with each other, contributing or raising issues - be nice, patient and respectful. We are all in it for the joy of cool science and software - let's make it about that. 

## Issues 
If you are facing issues running the software, or suspect bugs please report the following in the issue raised:

* ```beamshapes``` version 
* the smallest possible reproducible example that consistently generates the bug in your system
* OS and version
* Your Python type (direct installation, anaconda, etc.)

## Directivity implementations

If you're in doubt about whether an implementation is of relevance to the ```beamshapes``` package, raise an issue with the proposed model and let's discuss. For now, ```beamshapes``` is specifically looking to add source models that have a 'rich' front-back directivity ie. that are not front-back symmetric. 

If you wish to contribute a directivity implementation for a source-model see the generalised workflow page [here](https://beamshapes.readthedocs.io/en/latest/general_workflow.html). 
For reasons of code homogeneity, I request that you keep the template structure suggested in the [how-to](https://beamshapes.readthedocs.io/en/latest/general_workflow.html) page. All implementations must have tests that at least check the coded implementation matches the published results of the original paper/book that described the model. (If you are planning to implement a directivity that isn't published or have calculated yourself - we still need tests, but I'd be curious to hear how to check the validity of the implementation's results!)

##  Unit-tests
All unit-tests in the ```beamshapes``` package are written using the ```unittest``` package. Be aware that running all the tests will take a decent amount of time on an average PC/laptop (anywhere upwards of 20 minutes!).

### Running unit-tests
To run the unittests, first clone and install the ```beamshapes``` package (see the [README](README.md)).To run the tests you will need to install all the packages in the ```tests/requirements_test.txt``` file. From the root of the directory, open your Python command line tool of choice and run:

```python -m unittest```  

To run specific unit tests, use the commande ```python -m unittest beamshapes/tests/{name-of-test-here.py}```, where *{name-of-test-here.py}* is the 
test you want to run.

### Writing unit-tests
The unit-tests in ```beamshapes``` should test the computational implementation for correctness based on previously published results. See [here](beamshapes/tests/testing.md) for more on the testing approach. 

### Test coverage
If you'd like to estimate test coverage:

* Install the [```coverage```](https://coverage.readthedocs.io/en/6.0.2/) package: ```pip install coverage``` 
* Change directories and enter the repo ```cd bat_beamshapes```
* Run the tests with ```coverage``` and ```unittest``` with: ```coverage run -m unittest discover```
* To generate an overall report on coverage: ```coverage report```

## Documentation

### Building docs locally
To build the documentation locally, you must first install all the required packages:

```pip install -r docs/docs_requirements.txt``` 

And then:

* Move into the ```docs``` folder with ```cd docs```
* Make an empty ```build``` folder with ```mkdir build```
* Generate the docs with ```sphinx-build source build```

All the formatted web-pages are in the ```build/``` directory. Start with the ```index.html``` to begin browsing from the landing page. 



# Testing approach

All the source models in ```beamshapes``` are tested for their correctness by comparing
with published results. As of version 0.2, the models and their outputs are tested by
comparing the ```beamshapes``` output with plots in Beranek & Mellow 2012.

## Checking directivity patterns of source models 
Bernanek & Mellow 2012 calculate two types of directivity patterns: 1) the 'on-axis' pattern which describes the on-axis sound intensity across different *ka* values for instance, and 2) the 'full' directivity pattern, looking at the relative sound intensity (re. on-axis level) across the azimuth, and for different *ka* values. 

Tests pass if they are within the threshold of the published pattern. ## Speeding up symbolic calculations

Right now I'm stuck with SymPy calcultions that are taking for ever. 
I have 3 options right now, 

## What are other people using for the same task?


1. Mellow, Tim: *Mathematica*

1. Aarts Ronald, Janssen, Augustus J. E. M.: ronald.m.aarts@philips.com : *NO reply*
    * Comparing Sound Radiation from a Loudspeaker with that from a Flexible Spherical Cap on a Rigid Sphere 2011

1. Wojciech P.Rdzanek: wprdzank@ur.edu.pl :*replied, uses Mathematica*
    * The acoustic power of a vibrating clamped circularplate revisited in the wide low frequency rangeusing expansion into the radial polynomials 2016
1. Mark Poletti, NZ : *replied, uses MATLAB*


## Trying out FriCAS

* I've been getting a hang of FriCAS - it seems to be a very powerful, yet low level language. 
* Write and read text from/to files using the writeFile syntax http://fricas-wiki.math.uni.wroc.pl/SandBoxTextFiles

## I've been messing up the plain Hankel function of the 2nd kind nd the *spherical Hankel function of the 2nd kind!!!*
## Not having a formal math education means, when I encounter the derivative of a legendre function (Pn'(x)) -- I have no idea whether the derivative is with respect to n or x :P

* FriCAS is pretty cool as a system, though the learning curve made me stick to Python in the end. The confusion between multiple derivatives. 


## Troubleshooting 

*copy from ```ts_pistonsphere.py```*
### Troubleshooting piston in a sphere

> For ka>3 the calculations show values that deviate from the 
textbook groundtruth by 2-5 dB -- which is not ignorable!!
> I suspected the problem may come from:
    * low mpmath.mp.dps (decimal places) -- changing from 50-400 had no
    effects, which is odd
    * low N (matrix size/number of terms calculcated) -- changing from 
    baseline of 12+f(ka) --> 15+f(ka) had no effect
    * the directivity calculations were done with numpy pre 5th may, 
    and then changed to mpmath backend --> no effect. 

> 2021-05-65: I now suspect the problem lies perhaps with the quadrature 
terms. What if the quadrature is not 'accurate' enough? Here I'll test this idea
    * Some points to support this idea. The default quadrature method behind
    mpmath.quad is the 'tanh-sinh' algorithm. Instead of directly lambdifying 
    the Imn term into a standard mpmath.quad function I 'manually' made a 
    quadrature function for it to manipulate the options. 
    * For dps 200. Using the default 'tanh-sinh' leads to an estimated error 
    of e-203, while using 'gauss-legendre' leads to an estimated error of e-382. 
    Perhaps this is where the error is arising from. I noticed the integration 
    error increases with increasing m,n values. Perhaps this is why for bigger ka's, 
    (ka>3), the predictions get messier than for small ka's? There is at least 
    a connection here. 

The exact quadrature algorithm used seems to play a big difference. 

*The Imn term quadrature algorithm is NOT the problem* -- though choosing gauss-legendre
reduced execution time by 1/2!!!!

2021-05-06: Despite troubleshooting for so long on the piston in a sphere, I"m not having any luck understanding where the issue is coming from. 

At this point I have two options:
* stick with the piston in a sphere, but run all calculations with ka <= 3 
OR
* continue to troubleshoot with SAGE. 
* The challenge here is basically that 


2021-05-12:
I'm realising there is some effect of *N*, the number of terms used to get the beamshape. 
While the overall effect of increasing N for ka=3 is weak (an decrease of ~0.2 dB max|error|,
and decrease of 0.7 dB of mean|error|), the LUsolve residuals dramatically increase 2-3 orders. 

This tells me there may be issues with the numerical integration in general. Perhaps the numerical integration gets more messy with N terms.


Also -- WHAT ABOUT (PY)-FLINT. Check out  this [page](https://fredrikj.net/blog/2018/11/announcing-python-flint-0-2/), and this [one](https://fredrikj.net/python-flint/index.html)?


Issues faced while installing Flint
> first need to install MPIR. Downloaded source from Github.
    * the ```./configure``` command didn't work because the previous steps weren't complete -- needed to run ```autoreconf -i``` and then the ./configure 
> Also don't forget you have to install GMP+ MPFR too. 

One thing I realised is the installation flow for 'source' installations is 
> ./configure
> make 
> make check 
> make install --> this step typically needs permissions, and so `sudo make install` with 
the password prompt is what is actually needed. 

Installing python-flint with ```pip``` failed. Installing python-flint from Git source failed. Installing with Conda -- This worked!!! 

### What it's like to use py-flint

> Incredibly fast !!! eg. show the example of the Lm integral (eqn. 12.018) for a dps=500, it happens nearly instantaneously in py-flint, while with mpmath it takes 10's of seconds. 
2021-05-12
16:21 - I'm getting different values for P'm(cos(alpha)) between my Flint implementation using eqn. 12.98 and the native legendre_p(n,z).diff(z) version from sympy. 

16:44 -- and that's because there's a typo in eqn. 12.98 --> it should be sin(theta)^2 instead of sin(theta)

22:52 : FLINT is really, really fast now that I think of it. For dps=150, and ka=10 the Flint version of the code runs in ~9.5mins while the parallelised mpmath version takes about that long!!




# bat_beamshapes
 - Thejasvi Beleyur & Gaurav Dhariwal
Modelling bat beamshapes using various acoustic models, and trying to fit them with data. 

The beamshape models being implemented are:

1. Vibrating cap of a sphere: Aarts & Janssen 2010, Comparing sound radiation from a loudspeaker with that from a flexible spherical cap on a rigid sphere, confererence paper 2010

2. Piston in a closed finite baffle: Mellow and Kärkkäinen 2005, On the sound field of an oscillating disk in a finite open and closed circular baffle, JASA

3. Piston in an infinite baffle
---
title: "```beamshapes```: a Python package to generate directivity patterns for various sound source models"
authors:
  - name: Thejasvi Beleyur
    thanks:  1) Centre for the Advanced Study of Collective Behaviour, Konstanz 2) Acoustic and Functional Ecology, Max Planck Institute for Ornithology, Seewiesen 3) Max Planck Institute for Animal Behavior, Radolfzell
    email: thejasvib@gmail.com

abstract: |
  ```beamshapes``` is an open-source Python package that implements various directivity patterns for sound sources. While there is an abundance of published directivity patterns in the literature - their computational implementations often remain as in-house scripts in proprietary languages. ```beamshapes``` overcomes this gap, and provides acousticians and bioacousticians easily accessible implementations of sound source directivities.

keywords:
  - acoustics
  - bioacoustics
  - directivity 
  - sound radiation
  - Python 
  - echolocation

bibliography: ../../paper.bib
biblio-style: unsrt
output: 
  rticles::arxiv_article:
    keep_tex: true
    latex_engine: xelatex
---

# Summary

Sound sources such as human beings or loudspeakers often exhibit a 'directionality' in how loud they sound at different angles.
A listener or microphone placed at various angles at a fixed radius will often pick up sometimes drastically different sound levels. The same sound source may however sometimes produce omnidirectional sound fields. The directionality of a sound source can be modelled as a combination of the frequency of the emitted
sound and the geometry of the vibrating and non-vibrating parts of the sound source itself. The resulting pattern of sound radiation with angular location is called the *directivity* of the source [@beranek2012acoustics]. 

Directivity ($D$) describes the relative sound levels at a given angle $D(\theta)$
with relation to the level on-axis ($D(0)$, and thus $directivity = \frac{D_{\theta}}{D_{0}}$). Directivity
functions exist for a wide variety of sound sources that can be modelled analytically. A well known example 
of a directivity function is of the piston in an infinite baffle. The piston is a circular surface of radius
$a$, vibrating back and forth about a hole in an infinite wall or baffle. The directivity is described by $\frac{2J_{1}(ka \times sin \theta)}{ka \times sin \theta}$, where $J_{1}$ is the Bessel's function of the first kind, $k$ is the wavenumber, where $k = \frac{2\pi f}{c}$, where $f$ is the frequency of the sound, and $c$ is the speed of sound. The angle of the receiver is $\theta$, which varies from from 0-2$\pi$ radians in the azimuth. 

Directivity functions can be used in a two-fold manner 1) to deliberately engineer devices to suit particular specifications, eg. loudspeaker sound fields [@beranek2012acoustics] and 2) to infer parameters of a sound source itself having assumed a relevant model, eg. estimating direction of call emission [@guarato2011method] and mouth aperture of bat echolocation calls [@jakobsen2013convergent;@kounitsky2015bats].


![Directivity patterns ($\frac{D_{\theta}}{D_{0}}$) of the currently implemented sound sources for a common set of $ka$ values for comparison, where $k$ is the wavenumber ($\frac{2\pi\:f}{c}$) and $a$ is the piston radius - or its equivalent. The directivity pattern shows the ratio between the off-axis sound level at ($\theta^{\circ}$) to the on-axis level at ($0^{\circ}$) in decibels. A) piston in an infinite baffle B) piston in a sphere with the half-angle aperture $\alpha = 30^{\circ}$, and where the piston radius $a=R\:sin \:\alpha$ C) oscillating cap of a sphere with  $\alpha=30^{\circ}$, and the equivalent piston radius is $a=R\:sin \:\alpha$ D) vibrating point on a sphere, here the $kR=ka$ for comparison with the other models. \label{fig:pistonsphereinf}](../piston_sphere_baffle.png)


# Statement of need

A host of published directivity functions exist in the literature, but it is my experience that their computational implementations remain as inhouse scripts that are often in proprietary language platforms. To my knowledge there are no computational implementations of sound source directivities that are open-source and developed using modern software practices such as version-control and unit-testing. In this paper I present ```beamshapes```, a Python package that implements directivity patterns for sound source models. As of this publication, ```beamshapes``` version 0.2.0 implements four sound sources 1) piston in an infinite baffle, 2) point source of a sphere 3) oscillating cap of a sphere and 4) piston in a sphere. 

Computational implementations of directivity functions often require long run-times due to the intensive numerical routines and arbitrary precision math involved. Long run-times hinder scientific projects in reducing the number of models and parameter space that can be explored. ```beamshapes``` boasts parallelised code to generate significant speed-ups in run-times. 

The availability of openly-available directivity functions will hopefully stir the acoustics, and specifically the bio-acoustics community to rigorously test and compare sound radiation with models. The availability of multiple implementations allows comparison of data with multiple models. Until recently, computational power has perhaps been a limiting factor to calculating directivity functions. Models with easily calculable outputs have thus been favoured (eg. piston in an infinite baffle), especially in the field of bioacoustics [@strother1970acoustical;@mogensen1979sound]. 


Despite the recent availability of computational power, older, simpler models with limited biological relevance continue to dominate the field of bio-acoustics. For instance, the piston in an infinite baffle only predicts the beam-shape for a $\pm90^{\circ}$ range off-axis, and assumes front-back symmetry (\autoref{fig:pistonsphereinf} A). This is unrealistic for most vocalising animals, especially echolocating bats and odontocetes. However, the piston in an infinite baffle continues to be a standard reference model for multiple studies ranging over the past decade [@jakobsen2013convergent;@kounitsky2015bats;@macaulay2020high]. The piston in a sphere, for instance recreates many of the highly directional central and side lobes and that are seen in bats, while also predicting sound radiation behind the source (\autoref{fig:pistonsphereinf} B). In contrast to the echolocation literature, the bird [@witkin1977importance;@larsen1990directionality;@brumm2002sound;@patricelli2007differences;@patricelli2008acoustic;@yorzinski2010birds]  and frog call literature [@gerhardt1975sound;@rodriguez2020sound] for instance has been dominated by quantitative characterisations of sound radiation with no attempts at directly comparing measurements to model predictions. Using model-based directivity patterns to infer source parameters allows the discovery of common parameter spaces that birds occupy, and facilitates cross-species comparisons. The oscillating cap of a sphere and vibrating point on a sphere (\autoref{fig:pistonsphereinf} C,D) are two other models of potential relevance to bioacousticians attempting to describe the sound radiation of bird and frog calls for instance. 

Future releases of ```beamshapes``` are scheduled to include directivity patterns for additional models of interest such as rectangular cap of a sphere or piston in a closed finite baffle.


# Software packages used in this work
```beamshapes``` relies on the Python open-source ecosystem and is built on the numpy, scipy, sympy, mpmath and flint libraries [@2020NumPy;@2020SciPy;@meurer2017sympy;@mpmath;@hart2011flint]. 

# Package repository
```beamshapes``` can be currently accessed at https://github.com/thejasvibr/bat_beamshapes.git and the documentation with examples are hosted at https://beamshapes.readthedocs.io/en/latest/ .

# Acknowledgements
TB thanks Gaurav Dhariwal for his continual math advice and inputs, and Tim Mellow for 
strong support through sharing his Mathematica code and timely clarifications. TB thanks Lasse Jakobsen and Holger Goerlitz for productive discussions leading to this project. This work was executed through a combination of TB's private time, a DAAD stipend, the IMPRS for Organismal Biology,  and then finally by a CASCB Medium grant.

# References
Examples
========
A set of examples showing the capabilities  of the :code:`beamshapes` package. 

To get an overview of how to use `all` the currently implemented packages, start with `4 models at a time`. 
Misc 
====

.. automodule:: beamshapes.utilities 
       :members:


.. automodule:: beamshapes.flint_parallelisation 
       :members:
Source Models API
=================

.. automodule:: beamshapes.piston_in_infinite_baffle
       :members: d_theta_func, piston_in_infinite_baffle_directivity

.. automodule:: beamshapes.point_source_on_a_sphere
       :members: point_source_on_a_sphere_directivity, d_zero_func, d_theta_func

.. automodule:: beamshapes.cap_in_sphere
       :members: cap_in_sphere_directivity, d_zero, d_theta

.. automodule:: beamshapes.piston_in_sphere
       :members: d_theta, d_zero, piston_in_sphere_directivity

.. automodule:: beamshapes.piston_in_sphere_flint
       :members: make_Mmn_pll, dtheta, dzero, piston_in_sphere_directivity
Introduction
============
The :code:`beamshapes` package implements directivity functions (:math:`\frac{D_{\theta}}{D_{0}}`) of various published 
sound radiation models.

What is a directivity function? It's a function that quantifies how sound level changes as you change
the frequency of the emitted sound, and location of the receiver. 

Check out a general introduction to the concepts of the package :doc:`here <general_intro>`.

Why `beamshapes`?
~~~~~~~~~~~~~~~~~
While there are many sound radiation models described in the literature, there aren't
that many (also see `levitate <https://github.com/AppliedAcousticsChalmers/levitate/blob/master/levitate/transducers.py>`_ ) openly available computational implementations of their beamshapes. Existing packages focus on implementing directivities with analytical solutions - which can be calculated directly and quickly. `beamshapes` aims to increase the breadth of implemented directivities beyond those with analytical solutions. 

Who is this useful for? 
~~~~~~~~~~~~~~~~~~~~~~~
Acousticians and bio-acousticians looking to assess model-fits or perform 
parameter estimation on their sound sources. Check out more on the how
to use this package in the :doc:`examples <gallery_examples/index>`. 


Package installation
~~~~~~~~~~~~~~~~~~~~

`pip installation` : Install the latest stable version with :code:`pip install beamshapes` 

`Local installation` : Install from the `GitHub repo <https://github.com/thejasvibr/bat_beamshapes>`_ directly by cloning and 
following the instructions in the README.


Source models implemented
~~~~~~~~~~~~~~~~~~~~~~~~~
* Point source on a sphere
* Piston in an infinite baffle 
* Oscillating cap of a sphere 
* Piston in a sphere
Concepts behind :code:`beamshapes`
=============================
Sound sources don't radiate sound uniformly most of the time. Even us humans for instance, we emit more energy while talking to the front than to the back. The exact (non-uniform) pattern in which sound is radiated defines the 'beamshape' of a source. The beamshape is typically a combination of the frequency of emitted sound and the geometry of the vibrating surface and its associated (non-vibrating) surfaces. 

Let's go through the main concepts required to use the `beamshapes` package.

Source models
-------------
The source model refers to the source of sound and its geometric properties, and assumptions of the vibrations etc. According
to the situation in hand, different source models may be physically/biologically relevant! 


There are two main ways to predict how sound will radiate from a source - analytical or numerical modelling. Analytical models
start with equations defining the physics of sound radiation and move on to produce mathematical solutions. Numerical
methods use various numerical algorithms (eg. finite-element method) to computationally simulate sound radiation. :code:`beamshapes`
specifically implements analytical source models with pre-defined solutions. The advantage of using such analytical models is the lower
number of parameters needed to describe and understand the resulting sound radation. 

Below is a brief description of the source models and the parameters relevant to them. For more information on each of the model's please refer to 
the source references. 

Piston in an infinite baffle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
.. image:: ./_static/piston_in_inf_baffle.png.png
	:width: 200

`Piston in an infinite baffle schematic`

Here a rigid circular disk (the 'piston') vibrates back and forth across a hole (of matching size) set in a huge baffle .
The parameters needed to define this model are the wavenumber (`k` - see below for a list of all common abbreviations) and piston radius (`a`).  

Reference : Chp 13, Beranek, L. L., & Mellow, T. (2019). Acoustics: sound fields and transducers. Academic Press.

Point on a sphere
~~~~~~~~~~~~~~~~~
.. image:: ./_static/point_on_sphere.png.png
	:width: 200

`Point on a sphere schematic`

An infinitesimally small portion of a sphere's surface (the 'point') is considered to vibrate. The rest of the sphere does not vibrate.
The parameters needed to define this model are the wavenumber (`k`) and sphere radius (`R`). 

Reference: Chp 12, Beranek, L. L., & Mellow, T. (2019). Acoustics: sound fields and transducers.
Academic Press.

Oscillating cap of a sphere
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: ./_static/oscillating_cap_sphere.png.png
	:width: 200
	
`Oscillating cap of a sphere schematic`

The 'sliced' part of a sphere is the 'cap' in this case. The cap moves with an axial velocity :math:`u_{0}`.
Portions of the cap closer to the periphery vibrate less than portions close to the center, this is summarised by the 
relation :math:`u(R,\theta) = u_{0}cos \theta`, where :math:`\theta` is the distended angle from the cap's centre. 

The parameters needed to define this model are wavenumber (`k`), sphere radius (`R`), and the aperture angle of the cap (:math:`\alpha`). 

Reference: Chp 12, Beranek, L. L., & Mellow, T. (2019). Acoustics: sound fields and transducers.
Academic Press.

Piston in a sphere
~~~~~~~~~~~~~~~~~~

.. image:: ./_static/piston_in_a_sphere.png.png
	:width: 200

`Piston in a sphere schematic`

A sphere is sliced, and the 'cap' is discarded. The 'open' portion of the sliced sphere is now replaced with a piston. 
This piston in the sphere vibrates to produce sound. The parameters needed to define this model are wavenumber (`k`), sphere radius (`R`), 
aperture angle of the piston (:math:`\alpha`), and piston radius (`a`). 

Reference: Chp 12, Beranek, L. L., & Mellow, T. (2019). Acoustics: sound fields and transducers.
Academic Press.

Common parameters and abbreviations
-----------------------------------
The inputs will tend to be model-specific, but the common input parameters
to keep in mind are:

    #. `k`, wavenumber. This is :math:`\frac{2\pi}{\lambda}` - this is another way of defining the frequency of the vibration. :math:`\lambda` is the wavelength of the sound, also defined as :math:`\frac{v_\text{sound}}{\text{frequency}}`. 
    #. `a` : piston radius, wherever applicable
    #. `R` : sphere radius, wherever applicable. 
    #. :math:`\alpha, \theta, \phi`: various angles describing the size of the oscillating portion. These angles are in radians - not degrees!



Implementing a directivity function
===================================
This document describes the general workflow behind implementing the directivity 
function for a published sound radiation model. 

Pre-coding
~~~~~~~~~~
#. Find a suitable model with verifiable results (in the form of plots/numerical results)
#. Try and request underlying code for later comparison in case of deviations in results

Coding
~~~~~~
#. Implement as much of the model's components as SymPy objects. 
#. Always add the equation numbers in a comment above the variable/equation. 
#. Name the objects and intermediate variables to be as similar as possible to the names used in the publication. There will be some terms that are extremely long, split them where and when possible. Splitting long terms helps in later verification, and makes for pretty code. 
#. Convert the objects into functions using :code:`lambdify` and the appropriate backend (sympy, scipy, mpmath)
#. Put all the relevant code for the source model into one module. If necessary implement additional convenience functions in a separate model. 
#. Try switching between backends (:code:`numpy, scipy, sympy, mpmath`) if the code throws unexpected errors or is taking too long. 
#. The final directivity function should follow the `{name-of-the-model}_directivity`, and accept 2 inputs. The first input :code:`angles` should be an array/list-like object describing the location/s at which :math:`\frac{D_{\theta}}{D_{0}}` are to be calculated, and a :code:`params` dictionary object - which parametrises the source model. 

Result verification 
~~~~~~~~~~~~~~~~~~~

#. Replicate the key plots from the publication to ensure the correctness of your implementation. 

#. In case you can't directly compare the results by running the original publication's code - use a data digitisation tool like `WebPlotDigitizer <https://apps.automeris.io/wpd/>`_ to extract data from plots directly. The digitisation itself can have its own errors - so make sure to account for it while looking at discrepancies between your output and the published data.

#. If the results don't match - check the source of discrepancy. There is a high chance of a typo over the course of entering 10's of equations and variables! 

#. In case the code itself reflects the underlying equations correctly: check for 'patterns' in matching or discrepancy to try and isolate the problem. Does the match improve with the angle, does it get better at lower/higher frequencies, is it worse when the emitter is larger - or with increasing decimal precision (when using a :code:`mpmath` backend implementation). 


Optimisation
~~~~~~~~~~~~

* Certain models involve numerical methods that can take a long time to run (integration, summations, etc.). Check how to reduce run times without affecting output correctness eg. parallelising the code
* SymPy based implementations and compatible backends are strongly preferred. However, if a correct implementation takes too long (5-10 minutes each time) - consider switching to an alternative library like :code:`flint`.

To Do
~~~~~

* What when results don't match -- possible typo in the original publication/code -- how to proceed then. 



Notes for Piston in a Sphere
============================

updated 2021-07-18

The previous post on this page highlighted what seemed to be two discrepancies (check commit 46c11ec..), the first being a potential typo in equation 12.98, which described
:math:`\frac{\partial}{\partial \theta} P_n(cos \theta)` . Upon closer inspection I realised there was no typo, and it was an interpretational error on my part. 

However, discrepancy 2) related to a difference in `m` and `n` index order between the equations in the book and the Mathematica code implementation.
This post dives into more detail. 

The `m` and `n` index order discrepancy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

I. Expectations from substitutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The solution for :math:`K_{mn} = \int^{\pi}_{\alpha} P_{n}(\cos\theta) P_{m}(\cos\theta)\sin\theta\:d\theta` (eqn. 12.107) is given in App.II, eqn. 70. In the case where  :math:`m \neq n`, the solution is:

.. math:: 

    \frac{sin\:\alpha( P_{m}(cos\:\alpha)P^{\prime}_{n}(cos\:\alpha) - P_{n}(cos\:\alpha)P^{\prime}_{m}(cos\:\alpha))}{m(m+1) - n(n+1)}

Which we'll visually re-arrange for better comparison after substitution:

.. math::

    \frac{sin\:\alpha}{m(m+1) - n(n+1)}\bigg( P_{m}(cos\:\alpha)P^{\prime}_{n}(cos\:\alpha) - P_{n}(cos\:\alpha)P^{\prime}_{m}(cos\:\alpha) \bigg)



Where :math:`P_{n}(cos \:\theta)` is the Legendre polynomial of order `n` , and :math:`P^{\prime}_{n}(cos\:\theta)` (eqn. 12.98) is:

.. math::

    P^{\prime}_{n}(cos \theta) = \frac{\partial}{\partial \theta}P_{n}(cos \theta) = - \frac{n(n+1)}{(2n+1)sin \theta}(P_{n-1}(cos \theta) - P_{n+1}(cos \theta))
    
    = \frac{n(n+1)}{(2n+1)sin \theta}(P_{n+1}(cos \theta) - P_{n-1}(cos \theta))

When we do the substitutions for :math:`P^{\prime}_{n}(cos\:\alpha)`, :math:`P^{\prime}_{m}(cos\:\alpha)` (Appendix II, eqn.70) and :math:`\theta = \alpha`,
the full term is expected to be:

.. math::

    \frac{sin\:\alpha}{m(m+1) - n(n+1)} \\
    \left( P_{m}(cos\:\alpha)\frac{n(n+1)}{(2n+1)sin\:\alpha}(P_{n+1}(cos\:\alpha) - P_{n-1}(cos\:\alpha)) \\
     - P_{n}(cos\:\alpha)\frac{m(m+1)}{(2m+1)sin\:\alpha}(P_{m+1}(cos\:\alpha) - P_{m-1}(cos\:\alpha)) \right)

As of now the :code:`beamshapes`  piston in a sphere implementation follows the above equation. 
    
II. The textbook code implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The code implementation used to generate Fig. 12.23 has the equivalent of:

.. math::

    \frac{sin\:\alpha}{m(m+1) - n(n+1)} \\
    \left( P_{n}(cos\:\alpha)\frac{m(m+1)}{(2m+1)(sin\:\alpha)}(P_{m+1}(cos\:\alpha)-P_{m-1}(cos\:\alpha)) \\
     - P_{m}(cos\:\alpha)\frac{n(n+1)}{(2n+1)(sin\:\alpha)}(P_{n+1}(cos\:\alpha)-P_{n-1}(cos\:\alpha)) \right)


The `m` and `n` indices have been switched in the :math:`P_{m/n}(cos\:\alpha)` and the :math:`P^{\prime}_{m/n}(cos\:\alpha)` terms -- but the 
:math:`m(m+1) - n(n+1)` denominator term remains the same order as in section `I` .

III. Comparing directivity patterns from sections I and II
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The switch in `m` and `n` indices between II and III for leads to different  directivity patterns - here `ka=1` and `ka=3` is shown. 
 
.. image:: _static/pistoninsphere_deviation_2021-05-30.png
    :width: 400

.. image:: _static/pistoninsphere_deviation_2021-05-30_ka=3.png
    :width: 400


Perhaps the current :code:`beamshapes` implementation is the result of a coding error? Unlikely, as switching the order of the :math:`P_{m/n}` and :math:`P^{\prime}_{m/n}` terms in the :code:`beamshapes` implementations recreates Fig. 12.23 (not shown here). 

IV. The correct solution is the Mathematica code implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The code implementation is the correct solution - Tim Mellow has checked its correctness also through numerical integration.
The current implementations of the piston in a sphere are based on the Mathematica code implementation.

Acknowledgements
~~~~~~~~~~~~~~~~
Thanks to Tim Mellow for clarifying the discrepancy between published and implemented models, and Gaurav Dhariwal for re-checking the math once more. 

References
~~~~~~~~~~
* Chp 12, Beranek, L. L., & Mellow, T. (2012). Acoustics: sound fields and transducers. Academic Press. (also see the online Errata)
* To see code implementations check out the :code:`piston_in_sphere` documentation
Notes for 'Acoustic directivity of rectangular pistons on prolate spheroids'
============================================================================

* A 'prolate spheroid' is obtained when you rotate an ellipse on its major axis - a symmetric egg basically. 
* The coordinate system is defined by :math:`\xi, \eta, \phi` , just like how in spherical coordinates there's :math:`r,\theta,\phi`, except with the ranges (eqn 2):

.. math::

    1 \leq \xi < \infty

   -1 \leq \eta \leq 1

    0 \leq \phi \leq 2 \pi


Defined variables of interest
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* The `Directivity function` is given by eqn. 15
.. math::

    f(\theta,\phi) = \sum_{m=0}^\infty \sum_{l=m}^\infty \frac{\epsilon_{m}S^{(1)}_{ml}(h,cos\theta)}{R^{(4)\prime}_{ml}(h, \xi) N_{ml}}i^{l+1} \times [\tilde I^{s}_{ml} cos m\phi + \tilde I^{c}_{ml} sin m\phi]


* :math:`\epsilon_{m}` has a piecewise nature, where:

.. math::

    \epsilon_{m}=\begin{cases}
          1 \quad &\text{if} \, m = 0 \\
          1 \quad &\text{if} \, m \neq 0 \\
     \end{cases}

* `h`, 'size parameter' :math:`h=kd/2`, where `k` is the wavenumber and `d` is the 'interfocal distance of the generating ellipse'
* :math:`R^{(4)}_{ml}(h, \xi)` : prolate spheroidal radial function of the 4th kind  (eqn.5), where :math:`R^{(4)}_{ml}(h, \xi) = R^{(1)}_{ml}(h, \xi) - iR^{(2)}_{ml}(h, \xi)`. Defined in [2]
* :math:`S^{(1)}_{ml}(h, \eta)`: prolate spheroidal angle function of the 1st kind (eqn. 4). Defined in [2].
* :math:`N_{ml}` : prolate spheroidal angle normalization factor (eqn. 11)
* :math:`A_{ml}, B_{ml}` : unknown coefficients to be estimated -- related to the boundary condition of 'particle velocity at the spheroid-fluid interface'. Defined in [3].
* :math:`\tilde I^{s}_{ml} cos m\phi` (eqn. 12) and :math:`\tilde I^{c}_{ml} sin m\phi` '..define the size, shape, and location of the piston in the spheroisal baffle..'. In detail, 
.. math::

    \tilde I^{s}_{ml} = \int \int_{S_{i}} (\xi^{2}_{0} - \eta^{2})^{1/2}S^{(1)}_{ml}(h, \eta) cos m\phi \:d\eta d\phi

    \tilde I^{c}_{ml} = \int \int_{S_{i}} (\xi^{2}_{0} - \eta^{2})^{1/2}S^{(1)}_{ml}(h, \eta) sin m\phi \:d\eta d\phi

* The `directivity` itself is given by eqn. 16: 

.. math::

    F(\theta, \phi) = f(\theta,\phi)/f(\theta_{0},\phi_{0})

where :math:`\theta_{0},\phi_{0}` 'define the direction of maximum response'


Notes 
~~~~~
* :math:`R^{(1,2,4)}_{ml}(h, \xi)` and :math:`S^{(1)}_{ml}(h, \eta)` are 
* There seem to be something relevant `here <https://docs.scipy.org/doc/scipy/reference/special.html>`_ (Scipy implementations).


References
~~~~~~~~~~
#. Boisvert & Buren 2004, Acoustic directivity of rectangular pistons on prolate spheroids, JASA, 116, 1932 (2004); doi: 10.1121/1.1778840
#. C. Flammer, Spheroidal Wave Functions, Stanford University Press, Stanford, CA, 1957
#. J. E. Boisvert and A. L. Van Buren, ‘‘Acoustic radiation impedance of rectangular pistons on prolate spheroids,’’ J. Acoust. Soc. Am. 111, 867–874 (2002)
.. beamshapes documentation master file, created by
   sphinx-quickstart on Tue Apr 20 19:50:14 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Beamshapes: directivity patterns for various sound sources
==========================================================

.. include:: ./intro_why_who.rst

.. toctree::
      :maxdepth: 4
      :caption: Concepts

      ./general_intro.rst

.. toctree::
      :maxdepth: 4
      :caption: Use Cases

      gallery_examples/index.rst

Wishlist for future releases 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Piston in a finite closed baffle
* One-sided piston radiator (baffle=cylinder width)
* Rectangular cap of a sphere 
* Piston on a prolate spheroid

Interested?? Contribute, check out the general :doc:`workflow tips here <general_workflow>`.
All of the sound-source models and directivity calculations are from Leo Beranek & Tim Mellow's 
`Acoustics: sound fields and transducers. (Academic Press)` -  check it out for models
that may be of interest or to get an idea of how the wishlist models look like!



Contributors
~~~~~~~~~~~~
Thejasvi Beleyur (maintainer, thejasvib@gmail.com)

Gaurav Dhariwal


Acknowledgements
~~~~~~~~~~~~~~~~

Many thanks to Tim Mellow for sharing Mathematica code to help with porting to Python.
Also thanks to Holger R. Goerlitz for layout feedback (still in progress!) and 
Neetash MR for inspiring the package logo!


.. toctree::
   :maxdepth:1
   :caption: API reference:

   source_models.rst
   misc.rst

.. toctree::
    :maxdepth: 1
    :caption: dev notes

    developer_notes.rst


    
    
    
    



Model notes
~~~~~~~~~~~
Here are notes for various source model implementations - either under-way or those with some kinks that need to be sorted out. 

* :doc:`Implementing a directivity function <general_workflow>`

* :doc:`Piston on a prolate spheroid <notes_piston_on_prolate_spheroid>`

* :doc:`Piston in a sphere - deviations from fig. 12.23 <notes_piston_in_sphere>`

* :doc:`Rectangular cap of a sphere - implementation notes <notes_rectcap_of_sphere>`
Rectangular cap of a sphere
===========================

.. |Amn| replace:: :math:`A_{mn}`
.. |A0n| replace:: :math:`A_{0n}`
.. |Imn| replace:: :math:`I_{mn}`
.. |I0n| replace:: :math:`I_{0n}`
.. |dthetaphi| replace:: :math:`D(\theta,\pi)`
.. |d00| replace:: :math:`D(0,0)`
.. |date| date::

`Last updated`: |date|

The rectangular cap of a sphere models a square/rectangular part of a sphere's surface
vibrating (in contrast to the cap of  sphere, which models a circular cap). 

Here we'll try a top-down approach to understand the relevant variables and their corresponding code implementations in further detail. 

Directivity functions |dthetaphi| , |d00|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The off-axis directivity function |dthetaphi| is given by: 

.. math::
         
    D(\theta,\pi) = -\frac{4\pi}{k^{2}S}\sum^{N}_{n=0}\sum^{n/2}_{m=0}A_{mn}j^{n}P^{2m}_{n}(cos \:\theta)cos\:2m\phi \:\:\: (12.83)

The on-axis directivity function |d00| is given by: 

.. math::

    D(0,0) = -\frac{4\pi}{k^{2}S}\sum^{N}_{n=0}A_{0n}j^{n} \:\:\: (12.85)

The variables :math:`S, A_{mn}, A_{0n}` are defined below. 

.. math::

    S = 4R^{2}\Bigg(arctan\bigg(\frac{tan\:\alpha\:tan\:\beta}{sec^2\:\alpha + \sqrt{sec^{2}\:\alpha + tan^{2}\:\beta}}\bigg) \\ + arctan\bigg(\frac{tan\:\alpha\:tan\:\beta}{sec^2\:\beta + \sqrt{sec^{2}\:\beta + tan^{2}\:\alpha}}\bigg)\Bigg) \quad (12.69) \\
    \\ 

    A_{mn} = \frac{(2n+1)^2(n-2m)!I_{mn}}{j2\pi(n+2m)!\bigg(nh^{(2)}_{n-1}(kR) - (n+1)h^{(2)}_{n+1}(kR)\bigg)} \quad (12.75)

And |A0n| is all |Amn| where :math:`m=0`. 

Definition of |Imn|:
~~~~~~~~~~~~~~~~~~~~
|Amn| has the variable |Imn|, defined as: 

.. math::
    
    I_{mn} = \\
    \int^{arctan\frac{tan\:\beta}{tan\:\alpha}}_{0} cos\:2m\phi\:\int^{arctan\frac{tan\:\alpha}{cos\:\phi}}_{0} P^{2m}_{n}(cos\:\theta)sin\:\theta\:d\theta d\phi \\ 
    + \int^{\frac{\pi}{2}+arctan\frac{tan\:\alpha}{tan\:\beta}}_{\frac{\pi}{2}-arctan\frac{tan\:\alpha}{tan\:\beta}} cos\:2m\phi \:\int^{arctan\frac{tan\:\beta}{sin\:\phi}}_{0} P^{2m}_{n}(cos\:\theta)sin\:\theta\:d\theta d\phi  \\
    + \int^{\pi}_{\pi-arctan \frac{tan\:\beta}{tan\:\alpha}} cos\:2m\phi\:\int^{arctan\:\frac{tan\:\alpha}{-cos\phi}}_{0} P^{2m}_{n}(cos\:\theta)sin\:\theta\:d\theta d\phi \\
    \quad (12.76)

And |I0n| is all |Imn| where :math:`m=0`, given by: 

.. math::
    
    I_{0n} = \\
    \int^{arctan \frac{tan\:\beta}{tan\:\alpha}}_{0} \frac{tan\:\alpha}{\sqrt{cos^{2}\:\phi + tan^{2}\:\alpha}}
    P^{-1}_{n}\bigg(\frac{\cos\:\phi}{\sqrt{cos^{2}\:\phi + tan^{2}\:\alpha}}\bigg)d\phi \\
    + \int^{\pi/2+arctan\frac{tan\:\alpha}{tan\:\beta}}_{\pi/2-arctan\frac{tan\\:alpha}{tan\:\beta}} \frac{tan\:\beta}{\sqrt{sin^{2}\:\phi + tan^{2}\:\beta}}
    P^{-1}_{n}\bigg(\frac{sin\:\phi}{\sqrt{sin^{2}\:\phi + tan^{2}\:\beta}}\bigg)d\phi \\
    + \int^{\pi}_{\pi-arctan\frac{tan\:\beta}{tan\:\alpha}} \frac{tan\:\alpha}{\sqrt{cos^{2}\:\phi + tan^{2}\:\alpha}}
    P^{-1}_{n}\bigg(\frac{-cos\:\phi}{\sqrt{cos^{2}\:\phi + tan^{2}\:\alpha}}\bigg)d\phi \\
    \quad (12.77)


Computational implementation of |Imn|
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The code implementation looks different from these equations, and that is because the solutions to some of the integrals seem to 
have been included already. 

Let's examine them one-by-one. 


|Imn| is coded with the equivalent of :code:`Int[m_,n_,beta_]`:

.. math::

    I_{mn} = \\
    \frac{1}{PQ}arctan\frac{sin\:\beta}{sin\:\alpha}\sum^{P}_{p=1}\sum^{Q}_{q=1}arctan\:\frac{tan\:\alpha}{cos\:\phi_{1}}cos\:2m\phi_{1}\:P^{2m}_{n}(cos\theta_{1})\:sin\:\theta_{1} \\
    + \frac{2}{PQ}arctan\frac{sin\:\alpha}{sin\:\beta}\sum^{P}_{p=1}\sum^{Q}_{q=1}arctan\:\frac{tan\:\beta}{sin\:\phi_{2}}cos\:2m\phi_{2}\:P^{2m}_{n}(cos\theta_{2})\:sin\:\theta_{2} \\
    + \frac{1}{PQ}arctan\frac{sin\:\beta}{sin\:\alpha}\sum^{P}_{p=1}\sum^{Q}_{q=1}arctan\:\frac{tan\:\alpha}{-cos\:\phi_{3}}cos\:2m\phi_{3}\:P^{2m}_{n}(cos\theta_{3})\:sin\:\theta_{3} \\

where,

.. math::
    P=100, Q=100 \\
    \phi_{1}(p,\beta) = \frac{p-1/2}{P}arctan\:\frac{sin\:\beta}{sin\:\alpha} \\
    \phi_{2}(p,\beta) = \frac{\pi}{2} - arctan\:\frac{sin\:\alpha}{sin\:\beta} + 2\frac{p-1/2}{P}arctan\:\frac{sin\:\alpha}{sin\:\beta} \\
    \phi_{3}(p.\beta) = \pi - arctan\:\frac{sin\:\beta}{sin\:\alpha} + \frac{p-1/2}{P}arctan\:\frac{sin\:\beta}{sin\:\alpha} \\
    \theta_{1}(p,q,\beta) = \frac{q-1/2}{Q}\bigg(arctan\:\frac{tan\:\alpha}{cos\:\phi_{1}} \bigg) \\
    \theta_{2}(p,q,\beta) = \frac{q-1/2}{Q}\bigg(arctan\:\frac{tan\:\beta}{sin\:\phi_{2}} \bigg) \\
    \theta_{3}(p,q,\beta) = \frac{q-1/2}{Q}\bigg(arctan\:\frac{tan\:\alpha}{-cos\:\phi_{3}} \bigg) \\

Comparing the math |Imn| definitions and code implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let's begin to compare the code implementation of the solution with the original |Imn| terms:

Math |Imn| term1:

.. math::

    \int^{arctan\frac{tan\:\beta}{tan\:\alpha}}_{0} cos\:2m\phi\:\int^{arctan\frac{tan\:\alpha}{cos\:\phi}}_{0} P^{2m}_{n}(cos\:\theta)sin\:\theta\:d\theta d\phi \\ 
    
    implemented\:as\: :

    \frac{1}{PQ}arctan\frac{sin\:\beta}{sin\:\alpha}\sum^{P}_{p=1}\sum^{Q}_{q=1}arctan\:\frac{tan\:\alpha}{cos\:\phi_{1}}cos\:2m\phi_{1}\:P^{2m}_{n}(cos\theta_{1})\:sin\:\theta_{1} \\

|Imn| term2:

.. math::

    \int^{\frac{\pi}{2}+arctan\frac{tan\:\alpha}{tan\:\beta}}_{\frac{\pi}{2}-arctan\frac{tan\:\alpha}{tan\:\beta}} cos\:2m\phi \:\int^{arctan\frac{tan\:\beta}{sin\:\phi}}_{0} P^{2m}_{n}(cos\:\theta)sin\:\theta\:d\theta d\phi  \\    
    
    implemented\:as\: :

    \frac{2}{PQ}arctan\frac{sin\:\alpha}{sin\:\beta}\sum^{P}_{p=1}\sum^{Q}_{q=1}arctan\:\frac{tan\:\beta}{sin\:\phi_{2}}cos\:2m\phi_{2}\:P^{2m}_{n}(cos\theta_{2})\:sin\:\theta_{2} \\

|Imn| term 3:

.. math::

    \int^{\pi}_{\pi-arctan \frac{tan\:\beta}{tan\:\alpha}} cos\:2m\phi\:\int^{arctan\:\frac{tan\:\alpha}{-cos\phi}}_{0} P^{2m}_{n}(cos\:\theta)sin\:\theta\:d\theta d\phi \\
    
    implemented\:as\: :

     \frac{1}{PQ}arctan\frac{sin\:\beta}{sin\:\alpha}\sum^{P}_{p=1}\sum^{Q}_{q=1}arctan\:\frac{tan\:\alpha}{-cos\:\phi_{3}}cos\:2m\phi_{3}\:P^{2m}_{n}(cos\theta_{3})\:sin\:\theta_{3} \\


Definition of |I0n|:
~~~~~~~~~~~~~~~~~~~~

|I0n| (all `m=0`) terms:

.. math:: 

    \int^{arctan\:\frac{tan\:\beta}{tan\:\alpha}}_{0}\:\frac{tan\:\alpha}{\sqrt{cos^{2}\:\phi+tan^{2}\:\alpha}}
    P^{-1}_{n}\bigg( \frac{cos\:\phi}{\sqrt{cos^{2}\:\phi+tan^{2}\:\alpha}}\bigg)d\phi \\
    + \int^{\pi/2 + arctan\:\frac{tan\:\alpha}{tan\:\beta}}_{\pi/2-arctan\:\frac{tan\:\alpha}{tan\:\beta}}\:\frac{tan\:\beta}{\sqrt{sin^{2}\:\phi+tan^{2}\:\beta}}
    P^{-1}_{n}\bigg( \frac{sin\:\phi}{\sqrt{sin^{2}\:\phi+tan^{2}\:\beta}}\bigg)d\phi \\
    + \int^{\pi}_{\pi-arctan\:\frac{tan\:\beta}{tan\:\alpha}}\:\frac{tan\:\alpha}{\sqrt{cos^{2}\:\phi+tan^{2}\:\alpha}}
    P^{-1}_{n}\bigg( \frac{-cos\:\phi}{\sqrt{cos^{2}\:\phi+tan^{2}\:\alpha}}\bigg)d\phi \\

Computational implementation of |I0n|:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

    \frac{1}{P} arctan\:\frac{sin\:\beta}{sin\:\alpha} \sum^{P}_{p=1} \frac{tan\:\alpha}{\sqrt{cos^{2}\:\phi_{1}+tan^{2}\alpha}}
    P^{-1}_{n}\bigg(\frac{cos\:\phi_{1}}{\sqrt{cos^{2}\:\phi_{1}+tan^{2}\:\alpha}}\bigg) \\
    + \frac{2}{P} arctan\:\frac{sin\:\alpha}{sin\:\beta}\sum^{P}_{p=1}\frac{tan\:\beta}{\sqrt{sin^{2}\:\phi_{2}+tan^{2}\beta}}
    P^{-1}_{n}\bigg(\frac{sin\:\phi_{2}}{\sqrt{sin^{2}\:\phi_{2}+tan^{2}\:\beta}}\bigg) \\
    \frac{1}{P} arctan\:\frac{sin\:\beta}{sin\:\alpha} \sum^{P}_{p=1} \frac{tan\:\alpha}{\sqrt{cos^{2}\:\phi_{3}+tan^{2}\alpha}}
    P^{-1}_{n}\bigg(\frac{-cos\:\phi_{3}}{\sqrt{cos^{2}\:\phi_{3}+tan^{2}\:\alpha}}\bigg) \\

Comparing |I0n| in code and definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The |I0n| term doesn't have any substitutions/solutions in the code implementation it. The main difference is that the numerical integration is built into the 
implementation (an alternative approach to use a inbuilt numerical integration routine). 


Acknowledgements
----------------
Thanks to Tim Mellow for sharing the `Mathematica` code behind the textbook figures. 
