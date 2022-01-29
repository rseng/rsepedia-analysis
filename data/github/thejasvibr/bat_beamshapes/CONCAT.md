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
