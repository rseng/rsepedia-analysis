---
title: "PyGRB: A pure Python gamma-ray burst analysis package."
tags:
  - python
  - astrophysics
  - cosmology
  - gamma-ray bursts
  - gravitational lensing
authors:
  - name: James R. Paynter
    orcid: 0000-0002-9672-4008
    affiliation: 1
affiliations:
- name: School of Physics, University of Melbourne, Parkville, Victoria, 3010, Australia
  index: 1
date: 16 July 2020
bibliography: paper.bib
---


# Introduction


Gamma-ray bursts (GRBs) are short and intense bursts of low energy (keV -- MeV) gamma radiation.
Their cosmological origin and transient nature makes them important probes of the universe and its structure.
Since their discovery astronomers have sought to model the high energy emission.
A popular and enduring model, although phenomenological, is the fast-rise exponential-decay (FRED) pulse [@norris2005; @norris1996],

$$
S(t|A,\Delta,\tau,\xi) = A \exp \left[ - \xi \left(  \frac{t - \Delta}{\tau} + \frac{\tau}{t-\Delta}  \right)   \right].
$$

# Statement of Need

The analysis of large amounts of light-curves requires downloading FITS files from the relevant server.
To do this by hand is a tiring prospect for the 2,704 GRBs observed by the Burst and Transient Source Explorer (BATSE) experiment [@batse].
Having downloaded a FITS file, a scientist would then need to unpack the data from the file, and extract the relevant tables to construct a light-curve.
They may then want to plot the light-curve for publication, requiring them to write more software to appropriately represent the data.
Ultimately, they may want to look at population statistics, or compare different GRB pulses.
There is a gap in the market for software to get researchers off the ground quickly.
This is where `PyGRB` comes in.


# PyGRB

`PyGRB` is a pure Python, open source pulse-fitting package which aims to bring gamma-ray burst light-curve fitting and analysis into the 21st century.
`PyGRB` is able to download the pre-binned BATSE light curves (``bfits`` files), in addition to ``tte_list`` time-tagged photon arrival times.
FITS I/O functionality is provided by `Astropy` [@astropy]. `PyGRB` is built on top of the `Bilby` Bayesian inference library [@bilby], through which `PyGRB` utilises the `Dynesty` [@dynesty] and `Nestle` [@nestle] nested sampling packages [@skilling; @feroz; @multinest].

`PyGRB` makes visually appealing and scientifically instructive light-curves from the four broadband energy channel BATSE data.
The main feature of `PyGRB` is its ability to fit analytic light-curve models to data.
In particular, Bayesian model selection allows the user to determine the most appropriate pulse parameterisation for a particular burst.
Available pulse parameterisations are Gaussian pulses, FRED pulses and FRED variations.
Residual fitting is additionally possible, for which we implement a sine-Gaussian function.

![BATSE trigger 7475, GRB 990316 with FRED fit by `PyGRB`](../docs/source/images/B_7475__d_NL200__rates_F.png)


Ultimately, the model selection of `PyGRB` is used to determine if two light-curves are statistically identical, which would be indicative of a gravitational lensing event [@paczynski; @blaes; @mao].
It is quite difficult to compare GRB light-curves occurring at different times due to the variability of the gamma-ray background.
Comparing GRBs observed by different satellites is another matter altogether, owing to the different energy sensitivities, time resolution, and detector geometry.
`PyGRB` creates a unified, abstracted framework allowing for the comparison of gamma-ray bursts based on their fitted pulse parameters, rather than visual or bin-wise statistical comparisons of their light-curves, which is inherently fraught with opportunities for mishap.

In the present release only BATSE functionality is available.
However, due to the modular nature of the code, including additional GRB catalogues is as simple as creating the relevant `fetch`, and `preprocess` modules.
Future releases will allow for the easy comparison of gamma-ray bursts observed by different satellites.
As `PyGRB` is open source, the community is actively encouraged to contribute functionality for the many available GRB catalogues.

`PyGRB` is released under the BSD 3-Clause license.
The source code may be found at https://github.com/JamesPaynter/PyGRB, or alternatively the package may be installed from PyPi via ``pip install PyGRB``.
The online documentation, tutorials and examples are hosted at https://pygrb.readthedocs.io.

# Acknowledgements

I would like to thank Rachel Webster for introducing me to the problem of gravitationally lensed gamma-ray bursts, which spawned this project.
I would like to thank Eric Thrane for introducing me to modern computational statistics, particularly through `Bilby`, on whose shoulders this package stands.
Additionally I would like to thank Julian Carlin, Aman Chokshi for supporting me on my journey from novice programmer to published developer.

# References
|JOSS| |DOI| |pypi| |version| |Travis| |Coverage| |Docs| |AstroPy|

.. image:: https://github.com/JamesPaynter/PyGRB/blob/master/docs/source/images/logo.png
    :align: center
    :alt: PyGRB logo

.. inclusion-marker-one-liner-start

A gamma-ray burst (GRB) light-curve analysis package.

.. inclusion-marker-one-liner-end



.. inclusion-marker-what-it-does-start

Introduction
------------
*PyGRB* is a package to download gamma-ray burst (GRB) .FITS files from the relevant data archives (eg. NASA HEARSAC).
At the moment only `BATSE <https://heasarc.gsfc.nasa.gov/FTP/compton/data/batse/>`__ data can be downloaded and analysed with the software, although with only slight tweaks GRBs from other satellites can be easily analysed.
The code is then able to create light-curves from either pre-binned data or time-tagged photon-event data.
Light-curves may then be fitted with with pulse models, for further analysis.
Model fitting is done with nested sampling, powered by `Bilby <https://lscsoft.docs.ligo.org/bilby/index.html>`__, and `Dynesty <https://dynesty.readthedocs.io/>`__ and/or `Nestle <https://github.com/kbarbary/nestle>`__.


Installation
^^^^^^^^^^^^
*PyGRB* may be installed manually through cloning the repository

.. code-block:: console

  $ git clone https://github.com/JamesPaynter/PyGRB
  $ cd PyGRB
  $ pip install -r requirements.txt
  $ pip install .

or by downloading the compiled version from `PyPI <https://pypi.org/project/PyGRB/>`__

.. code-block:: console

  $ pip install pygrb


Installation of *PyGRB* and its dependencies should take no longer than a couple of minutes.

Then import *PyGRB* through ``import PyGRB``.

.. inclusion-marker-what-it-does-end


.. inclusion-marker-pulse-types-start

Pulse types
------------
Description of GRB pulse phenomenology.

.. image:: https://github.com/JamesPaynter/PyGRB/blob/master/docs/source/images/equations/FRED.gif
    :align: center
    :alt: FRED eqn: $I(t) = A \exp{ - \xi \left( \frac{t - \Delta}{\tau} + \frac{\tau}{t-\Delta} \right)}$


.. inclusion-marker-pulse-types-end

`See documentation for more <https://pygrb.readthedocs.io/en/latest/user/pulses.html>`__



.. role:: python(code)
   :language: python

.. image:: https://github.com/JamesPaynter/PyGRB/blob/master/docs/source/images/BATSE_trigger_7475_rates_rates.png
    :align: center
    :alt: BATSE trigger 7475


Usage
------

.. inclusion-marker-usage-start

Say we would like to fit a GRB light-curve such as the above, and determine its pulse parameters.
First we must load the relevant modules.

.. code-block:: python

  from PyGRB.main.fitpulse import PulseFitter
  from PyGRB.backend.makemodels import create_model_from_key


The :python:`PulseFitter` class is the main workhorse of the software.

.. code-block:: python

  GRB = PulseFitter(7475, times = (-2, 60),
            datatype = 'discsc', nSamples = 200, sampler = 'nestle',
            priors_pulse_start = -5, priors_pulse_end = 30)


The first argument specifies the BATSE trigger to be analysed, in this case trigger 7475.
Times can either be specified as :python:`'T90'`, :python:`'full'`, or a tuple of start and end times.
In the case of trigger 7475, most of the action happens over about (-2, 60), so we choose this interval for our times.
The :python:`nSamples` parameter determines how many live points the nested sampler is initiated with.
The :python:`sampler` parameter is used to choose between samplers.
The :python:`priors_pulse_start` and :python:`priors_pulse_end` parameters are used to set the (uniform) interval over which the program will allow the pulse start times.
The :python:`datatype` parameter specifies which kind of data we would like to download and analyse.
Typically :python:`'discsc'` is the most useful.
:python:`'tte'` is better for short GRBs.
The data will be downloaded and stored in :code:`data/`.



:python:`create_model_from_key` allows us to specify pulse models based on a simple key. The simple pulse type, a fast-rise exponential-decay (FRED) pulse, is utilised by

.. code-block:: python

  key = 'F'
  model = create_model_from_key(key)


Finally, we run the model through the sampler

.. code-block:: python

  GRB.main_multi_channel(channels = [0, 1, 2, 3], model = model)


The data products are stored in :code:`products/`.


.. inclusion-marker-usage-end


We should be left with a light-curve that looks like this:

.. image:: https://github.com/JamesPaynter/PyGRB/blob/master/docs/source/images/B_7475__d_NL200__rates_F.png
    :align: center
    :alt: BATSE trigger 7475


`See documentation for more <https://pygrb.readthedocs.io/en/latest/user/usage.html>`__


Under the Hood
---------------


.. image:: https://github.com/JamesPaynter/PyGRB/blob/master/docs/source/images/pulse_fit_animation.gif
    :align: center
    :alt: a GRB light-curve fit animation

There is a typo in this animation, the two fractions should take the same sign (+ve).
The -2 is an amplitude normalisation factor.


`See documentation for more <https://pygrb.readthedocs.io/en/latest/user/sampling.html>`__


Contribute
----------

'PyGRB' is an open-source software package freely available under the BSD 3-Clause License.
Users may request new features by opening a `GitHub Issue`_, or may contribute their own additions and improvements via a pull request.
Similarly, if you run into problems while using `PyGRB`, or require technical support, do not hesitate to request support through a `GitHub Issue`_.
If you use `PyGRB` in your work and would like to further collaborate on GRBs or gravitational lensing, I would be more than willing to discuss it over email or `GitHub Issue`_.

An incomplete list of possible improvements:

- Include support for uneven bin sizes and data gaps.

- Include compatibility with other GRB catalogues that are publicly available.

  - `Swift BAT <https://swift.gsfc.nasa.gov/results/batgrbcat/>`__

  - `Fermi GBM <https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/bursts/>`__

  - `Konus Wind <https://gcn.gsfc.nasa.gov/konus_grbs.html>`__

- Include capability to download and plot GRB spectra in addition to light-curves.

- Increase coverage to 100%



.. _GitHub Issue: https://github.com/JamesPaynter/PyGRB/issues

.. |AstroPy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org/
    :alt: astropy

.. |Travis| image:: https://travis-ci.com/JamesPaynter/PyGRB.svg?branch=master
    :alt: Travis Badge
    :target: https://travis-ci.com/JamesPaynter/PyGRB

.. |Coverage| image:: https://codecov.io/gh/JamesPaynter/PyGRB/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/JamesPaynter/PyGRB
    :alt: CodeCov - Coverage Status

.. |Docs| image:: https://readthedocs.org/projects/pygrb/badge/?version=latest
    :target: https://pygrb.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |pypi| image:: https://badge.fury.io/py/PyGRB.svg
    :target: https://badge.fury.io/py/PyGRB

.. |version| image:: https://img.shields.io/pypi/pyversions/pygrb.svg
    :target: https://pypi.org/project/pygrb/

.. |JOSS| image:: https://joss.theoj.org/papers/8aff0347e6993ec23b060052a80aaaa0/status.svg
    :target: https://joss.theoj.org/papers/8aff0347e6993ec23b060052a80aaaa0
    
.. |DOI| image:: https://zenodo.org/badge/213533787.svg
   :target: https://zenodo.org/badge/latestdoi/213533787
Welcome to PyGRB's documentation!
=================================

Instructions on how to use the package are provided in User Documentation.
Higher level API information can be find in API Reference.


* :ref:`user-docs`
* :ref:`api-docs`

.. _user-docs:

.. toctree::
   :maxdepth: 3
   :caption: User Documentation

   user/readme
   user/usage
   user/pulses
   user/sampling
   user/analysis


.. _api-docs:

.. toctree::
   :maxdepth: 5
   :caption: API Reference

   api/PyGRB.backend
   api/PyGRB.fetch
   api/PyGRB.main
   api/PyGRB.postprocess
   api/PyGRB.preprocess

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Usage
=====

.. role:: python(code)
   :language: python

.. figure:: ../images/BATSE_trigger_7475_rates_rates.png
   :figwidth: 80%
   :width: 80%
   :align: center
   :alt: BATSE trigger 7475

   BATSE trigger 7475


.. include::
  ../../../README.rst
  :start-after: inclusion-marker-usage-start
  :end-before: inclusion-marker-usage-end


----

The complete script for the above tutorial is here:

.. literalinclude:: ../../../examples/basic.py
    :name: basic.py
    :caption: basic.py

----

Something like the following will be printed to the console:

.. code-block:: console

  18:55 bilby INFO    : Single likelihood evaluation took 4.089e-04 s
  18:55 bilby INFO    : Using sampler Nestle with kwargs {'method': 'multi', 'npoints': 200, 'update_interval': None, 'npdim': None, 'maxiter': None, 'maxcall': None, 'dlogz': None, 'decline_factor': None, 'rstate': None, 'callback': <function print_progress at 0x00000185F8598AE8>, 'steps': 20, 'enlarge': 1.2}
  it=  6048 logz=-4174.731348
  18:55 bilby INFO    : Sampling time: 0:00:15.464772
  18:55 bilby INFO    : Summary of results:
  nsamples: 6249
  log_noise_evidence:    nan
  log_evidence: -4174.477 +/-  0.368
  log_bayes_factor:    nan +/-  0.368

The total number of posterior samples, `nsamples`, in this instance is 6249.
This is the length of posterior chains accessed through the code block at the bottom of this page.
The `log_noise_evidence` is used for `Bilby`'s gravitational wave inference.
It is irrelevant for `PyGRB`, as `PyGRB`'s Poisson likelihood does not specify an additional noise model.
The `log_evidence` is the Bayesian evidence in support of this data given the specified model.
The Bayes factor is the comparison of evidence between two models.
Here, `log_bayes_factor` is `NaN` as there is only one model specified, and the comparison cannot be made.
`PyGRB` does it's model comparison external to the inbuilt `Bilby` functionality.

Errors like the following can be ignored.

.. code-block:: console

  RuntimeWarning: overflow encountered in multiply ...

They are bad evaluations of the likelihood and have no real effect on the results.
In the unlikely event that they happen consistently, it may mean that the sampler is stuck in a bad region of parameter space.
Restarting the sampler for that particular evaluation should fix the problem.


This will create three sets of files in the products directory.

1.  Nested sampling posterior chains, which are stored as JSON files.

2.  Corner plots of the posteriors.

3.  Higher data products, the light-curves of each of the individual channels.
    These include by default analysis of the light-curve fit residuals to test for goodness-of-fit.


These are put into subdirectories based on the trigger number and the number of live points. Files for this example can be found under `/products/7475_model_comparison_200/`.


.. figure:: ../images/B_7475__d_NL200__rates_F.png
    :figwidth: 80%
    :width: 100%
    :align: center
    :alt: BATSE trigger 7475

    BATSE trigger 7475 with FRED fit


.. figure:: ../images/trigger7475chan3.PNG
    :figwidth: 80%
    :width: 100%
    :align: center
    :alt: BATSE trigger 7475 channel 3

    BATSE trigger 7475 with FRED pulse fit on a single channel (channel 3).
    The green shaded regions are the 1-sigma statistical (Poisson) uncertainty.
    The correlogram is a visualisation of the autocorrelation of the residuals in the second panel.
    For a "good" fit, one expects 95% (99%) of the points to lie within the 95% (99%) confidence intervals.
    It is interesting to note in this case there is a sinusoidal structure all the way to the end of the pulse, even though the deviations from the fit are not large.
    The probability plot tests the divergence of the residuals from zero for normality.
    Again we see that the residuals are normally distributed to a good approximation.


.. figure:: ../images/trigger7475post3.PNG
    :figwidth: 80%
    :width: 100%
    :align: center
    :alt: BATSE trigger 7475 channel 3 corner plot

    BATSE trigger 7475 with FRED pulse fit  on a single channel (channel 3) posterior corner plot.
    The parameter space is undersampled, (the posterior histograms are patchy), but that does not matter since this tutorial is for illustrative purposes only.



In answer to our initial question, pulse parameters may be read straight of the posterior distribution (or accessed through the posterior chain JSON file).

The following code block could be appended to the end of *basic.py*

.. code-block:: python

  # pick a channel to analyse
  channel = 0

  # create the right filename to open
  result_label = f'{GRB.fstring}{GRB.clabels[channel]}'
  open_result  = f'{GRB.outdir}/{result_label}_result.json'
  result = bilby.result.read_in_result(filename=open_result)

  # channel 0 parameters are subscripted with _a
  parameter = 'background_a'

  # the posterior samples are accessed through the result object
  posterior_samples = result.posterior[parameter].values

  # to find the median
  import numpy as np
  posterior_samples_median = np.median(posterior_samples)
.. figure:: ../images/logo.png
    :figwidth: 100%
    :width: 100%
    :align: center
    :alt: PyGRB logo


.. include::
  ../../../README.rst
  :start-after: inclusion-marker-what-it-does-start
  :end-before: inclusion-marker-what-it-does-end

Note
^^^^

Running `PyGRB` will implicitly import `Bilby`, which is used for Bayesian inference.
Since `Bilby` is primarily used for gravitational wave inference, importing `Bilby` will check for installations of `gwpy` and `lalsuite`, resulting in the following warnings if they are not installed.


.. code-block:: console

  13:23 bilby WARNING : You do not have gwpy installed currently. You will  not be able to use some of the prebuilt functions.
  13:23 bilby WARNING : You do not have lalsuite installed currently. You will not be able to use some of the prebuilt functions.
  13:23 bilby WARNING : You do not have gwpy installed currently. You will  not be able to use some of the prebuilt functions.
  13:23 bilby WARNING : You do not have lalsuite installed currently. You will not be able to use some of the prebuilt functions.
  13:23 bilby WARNING : You do not have gwpy installed currently. You will  not be able to use some of the prebuilt functions.
  13:23 bilby WARNING : You do not have lalsuite installed currently. You will not be able to use some of the prebuilt functions.
  13:23 bilby WARNING : You do not have lalsuite installed currently. You will not be able to use some of the prebuilt functions.


This has no effect on running of `PyGRB`.
.. _sampling:

Sampling
========

.. figure:: ../images/pulse_fit_animation.gif
    :figwidth: 80%
    :width: 100%
    :align: center
    :alt: a GRB light-curve fit animation

    Animation of nested sampling converging on channel 3 of BATSE trigger 8099 with a FRED pulse fit.


The equation in the nested sampling animation should read

.. math::

    S(t|A,\Delta,\tau,\xi) = A \exp \left[ - \xi \left(  \frac{t - \Delta}{\tau} + \frac{\tau}{t-\Delta}  \right)  -2 \right]


Photon Counting
---------------

Photon emission is a stochastic process.
The mean number of photons emitted during time interval :math:`t` is

.. math::

  \left< N(t) \right> = \lambda t,

where :math:`\left< N(t) \right>` is the expectation value for :math:`N` over the integration time, with Poisson rate :math:`\lambda`.
The probability distribution for the random variable :math:`N` is

.. math::

	\Pr(N(t) = n) = \frac{(\lambda t)^n e^{-\lambda t}}{n!}

For photons :math:`\vec{\gamma}=\{\gamma_1, \gamma_2, \ldots, \gamma_n\}` emitted via a homogeneous Poisson process (:math:`\lambda(t)=\lambda`), then the inter-arrival time between photon :math:`\gamma_i` and :math:`\gamma_{i+1}` is exponentially distributed.
The inter-arrival times between photon :math:`\gamma_i` and :math:`\gamma_{i+k}` follow a gamma distribution with

.. math::

	X \sim \Gamma(k, \lambda^{-1}),

where :math:`k=1` recovers the exponential distribution

.. math::

	X \sim \exp (\lambda^{-1}).

High energy detectors like BATSE accumulate photons at discrete multiples of their clock cycles (sampling frequencies).
BATSE has a sampling frequency of 500 kHz (a clock cycle of 2 :math:`{\mu s}`).
These detectors record each photon arrival time (time-tagged event, `tte`) to the nearest integer multiple of the clock cycle.
For BATSE, hardware limitations restricted this to the first 32,768 Large Area Detector (LAD) events, inclusive of all detectors.
Only the shortest, moderately bright :math:`{\gamma}`-ray bursts are contained completely within the `tte` data.
After the `tte` period has been exhausted, BATSE counts are collected in 64ms intervals (Discriminator Science Data, `discsc`) for :math:`{\sim 240}` seconds after the trigger.

.. figure:: ../images/T2611.png
   :figwidth: 70%
   :width: 100%
   :align: center
   :alt: BATSE trigger 2611

   BATSE trigger 2611. This is a short bright :math:`{\gamma}`-ray burst for which the `tte` data cuts out before the event is over.

Notes
^^^^^

Pre-trigger `tte` data is available for all 8 LADs.
Post-trigger `tte` data is available for the triggered detectors only.


BATSE sums the counts of the triggered detectors aboard the satellite.
These are then converted to rates in units of counts per second by dividing by the integration time (bin width).
To use Poisson statistics, `PyGRB` converts rates back to counts by multiplying by the bin widths, and forcing them to be integer.
This does not take into account the detector deadtime (see below).
Using rates rather than counts would underestimate the uncertainty in the radiation field, since for a Poisson distribution

.. math::

  \sigma_\text{Poisson}\sim \sqrt{N}.


A true Poisson process occurs at each detector.
The photon arrival times are convolved with the response of the detector, resulting in a Poisson distributed count spectrum.
Due to the onboard summing across triggered detectors we are forced to apply Poisson statistics to the summed counts, rather than the counts at each detector.


For :math:`{\gamma}`-ray bursts which are completely resolved in `tte` data, it is possible to analyse the counts at each detector.
However, since `PyGRB`'s main focus is the analysis of `discsc` data (prebinned), which is summed over the triggered detectors, this has not yet been implemented.
It is entirely possible to run the program independently over each triggered detector for `tte` data.
There is not yet a unified joint-likelihood framework which will consider a single radiation field across these detectors.


Dead Time
^^^^^^^^^

BATSE has a small dead time after each photon count of approximately one clock cycle.
This dead time is proportional to the energy of the incident photon (or particle event) which triggered the count :cite:`2010JGRAGjesteland`.

.. math::
  \tau \sim \alpha \ln \frac{E_{\gamma}}{E_0}

Where :math:`E_{\gamma}` is the energy of the incident photon, :math:`E_0= 5.5` keV is the reset level of the detector, and :math:`{\alpha}=0.75` :math:`{\mu s}` is the signal decay time :cite:`2008GeoRLGrefenstette`.
This means that photon counting is not a true Poisson process when the count rate approaches the sampling frequency.
Rather, it follows a truncated Poisson distribution.


Poisson rate
------------

The rate passed to the likelihood function is the sum of the individual pulses, specified in :ref:`pulses`.
An example rate for two FRED pulses would be:

.. math::

    \begin{split}
    S(t|A_1,\Delta_1,\tau_1,\xi_1,A_2,\Delta_2,\tau_2,\xi_2) =
                      &A_1 \exp \left[ - \xi_1 \left(  \frac{t - \Delta_1}{\tau_1}
                              + \frac{\tau_1}{t-\Delta_1}  \right)  -2 \right] \\
                    + &A_2 \exp \left[ - \xi_2 \left(  \frac{t - \Delta_2}{\tau_2}
                              + \frac{\tau_2}{t-\Delta_2}  \right)  -2 \right]
    \end{split}


Background
^^^^^^^^^^

The background is by default modelled as a constant.
This is sufficient for all but the longest gamma-ray bursts, excepting periods of unusually high background variability.
More complex background models, such as a polynomial, can be included by specifying them as a rate function and including relevant priors.

The complete rate is then:

.. math::

  \begin{split}
  S(t|A_1,\Delta_1,\tau_1,\xi_1,A_2,\Delta_2,\tau_2,\xi_2) = B +
                    &A_1 \exp \left[ - \xi_1 \left(  \frac{t - \Delta_1}{\tau_1}
                            + \frac{\tau_1}{t-\Delta_1}  \right)  -2 \right] \\
                  + &A_2 \exp \left[ - \xi_2 \left(  \frac{t - \Delta_2}{\tau_2}
                            + \frac{\tau_2}{t-\Delta_2}  \right)  -2 \right]
  \end{split}

Likelihood
----------

The rate function is then passed into the Poisson likelihood, which is a sum of the specified rates.


Priors
------

The default priors are

+---------------------------------+------------------------------------+--------------------+--------------+--------------+
| parameter                       | minimum                            | maximum            | type         | units        |
+---------------------------------+------------------------------------+--------------------+--------------+--------------+
| :math:`\Delta_i`                | \-\-                               | \-\-               | uniform      | seconds      |
+---------------------------------+------------------------------------+--------------------+--------------+--------------+
| :math:`\Delta_{i+1}`            | :math:`\Delta_i`                   | \-\-               | uniform      | seconds      |
+---------------------------------+------------------------------------+--------------------+--------------+--------------+
| :math:`B`                       | :math:`10^{-1}`                    | :math:`10^{3}`     | log\-uniform | counts / bin |
+---------------------------------+------------------------------------+--------------------+--------------+--------------+
| :math:`A`                       | :math:`10^{0}`                     | :math:`10^{5}`     | log\-uniform | counts / bin |
+---------------------------------+------------------------------------+--------------------+--------------+--------------+
| :math:`\tau`                    | :math:`10^{-3}`                    | :math:`10^{3}`     | log\-uniform | seconds      |
+---------------------------------+------------------------------------+--------------------+--------------+--------------+
| :math:`\xi`                     | :math:`10^{-3}`                    | :math:`10^{3}`     | log\-uniform | \-\-         |
+---------------------------------+------------------------------------+--------------------+--------------+--------------+
| :math:`\gamma`                  | :math:`10^{-1}`                    | :math:`10^{1}`     | log\-uniform | \-\-         |
+---------------------------------+------------------------------------+--------------------+--------------+--------------+
| :math:`\nu`                     | :math:`10^{-1}`                    | :math:`10^{1}`     | log\-uniform | \-\-         |
+---------------------------------+------------------------------------+--------------------+--------------+--------------+
| :math:`\Delta_\text{res}`       | \-\-                               | \-\-               | uniform      | seconds      |
+---------------------------------+------------------------------------+--------------------+--------------+--------------+
| :math:`A_\text{res}`            | :math:`10^{0}`                     | :math:`10^{3}`     | log\-uniform | counts / bin |
+---------------------------------+------------------------------------+--------------------+--------------+--------------+
| :math:`\tau_\text{res}`         | :math:`10^{-3}`                    | :math:`10^{3}`     | log\-uniform | counts / bin |
+---------------------------------+------------------------------------+--------------------+--------------+--------------+
| :math:`\omega`                  | :math:`10^{-3}`                    | :math:`10^{3}`     | log\-uniform | \-\-         |
+---------------------------------+------------------------------------+--------------------+--------------+--------------+
| :math:`\varphi`                 | :math:`-\pi`                       | :math:`\pi`        | uniform      | radians      |
+---------------------------------+------------------------------------+--------------------+--------------+--------------+

The priors on :math:`\Delta_i` and :math:`\Delta_\text{res}` are determined by the length of the light-curve passed to `PulseFitter`.


Further reading
---------------

For more on nested sampling, the reader is referred to the following links.

https://lscsoft.docs.ligo.org/bilby/basics-of-parameter-estimation.html

https://dynesty.readthedocs.io/en/latest/overview.html

https://dynesty.readthedocs.io/en/latest/dynamic.html




.. BATSE Data Types
.. ----------------
..
.. , and in (Time-to-Spill, `tts`),


References
----------

`BATSE Appendix G <https://heasarc.gsfc.nasa.gov/docs/cgro/nra/appendix_g.html#V.%20BATSE%20GUEST%20INVESTIGATOR%20PROGRAM>`_

.. bibliography:: ../refs.bib
  :style: unsrt
Application to Gravitational Lensing
====================================

This software package was created with the intent to discern similar looking gamma-ray bursts from statistically identical gamma-ray bursts.
A pair of statistically identical gamma-ray bursts would be indicative of a gravitational lensing event.

The material presented here is supplementary to a research paper currently out to review.
These scripts, using the *PyGRB* package, allow the reader to reproduce the results of the paper (with enough computational power).


.. toctree::
    :caption: Application to candidate lenses

    grb/973
    grb/3770
    grb/3770walk
    grb/3770ext
    grb/3770priors
.. _pulses:

Pulse Types
===========

.. role:: python(code)
   :language: python

A standard gamma-ray burst looks like this


.. figure:: ../images/BATSE_trigger_7475_rates_rates.png
    :figwidth: 80%
    :width: 100%
    :align: center
    :alt: BATSE trigger 7475

    BATSE trigger 7475


Gaussian pulse
--------------

A simple pulse parameterisation one might imagine is a Gaussian, like those used to model emission lines in eg. quasar spectra.
The equation for a Gaussian pulse is:

.. math::

    S(t|A,\Delta,\sigma) = A \exp \left[ \frac{\left( t - \Delta \right)^2}{2\sigma^2} \right]

.. figure:: ../images/B_7475__d_NL200__rates_G.png
    :figwidth: 80%
    :width: 100%
    :align: center
    :alt: BATSE trigger 7475 with gaussian fit

    BATSE trigger 7475 with gaussian fit


However, we can immediately see that such a pulse parameterisation does not catch the fast rise of the GRB pulse, nor the slow decay.
The residuals have consistent structure.

FRED pulse
----------

The standard pulse parameterisation used to model gamma-ray bursts is a fast-rise exponential-decay (FRED) curve.

.. math::

    S(t|A,\Delta,\tau,\xi) = A \exp \left[ - \xi \left(  \frac{t - \Delta}{\tau} + \frac{\tau}{t-\Delta}  \right)   \right]

.. figure:: ../images/B_7475__d_NL200__rates_F.png
    :figwidth: 80%
    :width: 100%
    :align: center
    :alt: BATSE trigger 7475 with FRED fit

    BATSE trigger 7475 with FRED fit


The fit is better than a Gaussian, but again there is structure un-accounted for in the residuals.

FRED-X pulse
------------

We try again with an extended fast-rise exponential-decay model, FRED-X.

.. math::

    S(t|A,\Delta,\tau,\xi,\gamma,\nu) = A \exp \left[ -\xi^\gamma \left(\frac{t - \Delta}{\tau}\right)^\gamma - \xi^\nu \left(\frac{\tau}{t-\Delta}\right)^\nu\right]

.. figure:: ../images/B_7475__d_NL200__rates_X.png
    :figwidth: 80%
    :width: 100%
    :align: center
    :alt: BATSE trigger 7475 with FRED-X fit

    BATSE trigger 7475 with FRED-X fit


Again, an element of structure in the residual persists.

Sine-Gaussian residual
----------------------

We use a sine-gaussian residual function to account for these residuals.

.. math::

    \text{res}(t)= A_\text{res} \exp \left[ - \left(\frac{t-\Delta_\text{res}} {\lambda_\text{res}}\right)^2 \right] \cos\left(\omega t + \varphi \right)

.. figure:: ../images/B_7475__d_NL200__rates_X.png
    :figwidth: 80%
    :width: 100%
    :align: center
    :alt: BATSE trigger 7475 with FRED fit and sine-Gaussian residual

    BATSE trigger 7475 with FRED fit and sine-Gaussian residual (not implemented yet)


A model selection script
------------------------

Say, then, that we would like to fit on of each of these models in turn to the light-curve.
To create a model, we specify a list of keys.
For a Gaussian, FRED, and FRED-X pulse, the keys are as follows:

.. code-block:: python

  keys = ['G', 'F', 'X']


Residuals are included by placing a lower case :python:`'s'` following the pulse to which it should be applied.
Residuals must be applied following a pulse, and cannot be used as a standalone.
Let's add an additional three pulse fits, all with a residual, to our list of keys.

.. code-block:: python

  keys += ['Gs', 'Fs', 'Xs']


So our complete list would look like

.. code-block:: python

  keys = ['G', 'F', 'X', 'Gs', 'Fs', 'Xs']



Now we need to convert our keys into models for the nested sampling analysis.

.. code-block:: python

  model_dict = {}
  for key in keys:
      model_dict[key] = create_model_from_key(key)
  models = [model for key, model in model_dict.items()]


Finally, we feed each of the models in turn to the sampler through :python:`main_multi_channel`.
This function further splits the model (in this case into 4), and tests each of the channels individually.

.. code-block:: python

  for model in models:
      GRB.main_multi_channel(channels = [0, 1, 2, 3], model = model)

----

The complete script for the above tutorial is here:

.. literalinclude:: ../../../examples/intermediate.py
    :name: intermediate.py
    :caption: intermediate.py

----

The model selection script will output a set of tables that looks something like the following (residual models not included here for computational brevity).

Channel 1

+-------+----------+-------+--------+
| Model |   ln Z   | error | ln BF  |
+-------+----------+-------+--------+
| G     | -4610.80 |  0.34 |  0.00  |
+-------+----------+-------+--------+
| F     | -4181.16 |  0.38 | 429.64 |
+-------+----------+-------+--------+
| X     | -4175.97 |  0.39 | 434.83 |
+-------+----------+-------+--------+

Channel 2

+-------+----------+-------+--------+
| Model |   ln Z   | error | ln BF  |
+-------+----------+-------+--------+
| G     | -4912.74 |  0.34 |  0.00  |
+-------+----------+-------+--------+
| F     | -4173.26 |  0.38 | 739.48 |
+-------+----------+-------+--------+
| X     | -4164.73 |  0.39 | 748.01 |
+-------+----------+-------+--------+

Channel 3

+-------+----------+-------+--------+
| Model |   ln Z   | error | ln BF  |
+-------+----------+-------+--------+
| G     | -4551.10 |  0.34 |  0.00  |
+-------+----------+-------+--------+
| F     | -4043.61 |  0.37 | 507.49 |
+-------+----------+-------+--------+
| X     | -4036.39 |  0.39 | 514.71 |
+-------+----------+-------+--------+

Channel 4

+-------+----------+-------+-------+
| Model |   ln Z   | error | ln BF |
+-------+----------+-------+-------+
| G     | -3714.97 |  0.28 |  0.00 |
+-------+----------+-------+-------+
| F     | -3709.29 |  0.31 |  5.68 |
+-------+----------+-------+-------+
| X     | -3709.44 |  0.31 |  5.53 |
+-------+----------+-------+-------+

Which tells us that a FRED-X model is preferred in this case for all channels.
.. _3770priors:

.. role:: python(code)
   :language: python


GRB 950830 priors
-----------------

Investigation of the effects of priors on model selection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the following script we test the effect of the priors on :math:`\gamma` and :math:`\nu` on the model selection.

To compare the same model with different prior ranges, such as over a longer time period, then save the two models in the dict under different names, eg.

.. code-block:: python

  keys = [key_1, key_2]
  model_dict['model 1'] = create_model_from_key(key_1)
  model_dict['model 2'] = create_model_from_key(key_2)
  models = [model for key, model in model_dict.items()]


This is utlised in the prior test script below.

----

.. literalinclude:: ../../../../examples/gravlensnaturepriors.py
    :name: gravlensnaturepriors.py
    :caption: gravlensnaturepriors.py

----


.. We find that
.. _3770:

.. role:: python(code)
   :language: python


GRB 950830 basic analysis
-------------------------

.. figure:: ../../images/trigger3770.PNG
   :figwidth: 80%
   :width: 100%
   :align: center
   :alt: BATSE trigger 3770 with FRED lens fit.

   BATSE trigger 3770 with FRED lens fit.


This script compares the two simplest pulse models for a gravitational lensing event for GRB 950830.


----

.. literalinclude:: ../../../../examples/gravlensnaturebasic.py
    :name: gravlensnaturebasic.py
    :caption: gravlensnaturebasic.py

----


Channel 1

+-------+---------+-------+-------+
| Model |   ln Z  | error | ln BF |
+-------+---------+-------+-------+
| FL    | -692.95 |  0.29 |  3.75 |
+-------+---------+-------+-------+
| FF    | -696.70 |  0.31 |  0.00 |
+-------+---------+-------+-------+

Channel 2

+-------+---------+-------+-------+
| Model |   ln Z  | error | ln BF |
+-------+---------+-------+-------+
| FL    | -663.86 |  0.32 |  6.67 |
+-------+---------+-------+-------+
| FF    | -670.53 |  0.36 |  0.00 |
+-------+---------+-------+-------+

Channel 3

+-------+---------+-------+-------+
| Model |   ln Z  | error | ln BF |
+-------+---------+-------+-------+
| FL    | -719.66 |  0.34 |  5.16 |
+-------+---------+-------+-------+
| FF    | -724.82 |  0.39 |  0.00 |
+-------+---------+-------+-------+

Channel 4

+-------+---------+-------+-------+
| Model |   ln Z  | error | ln BF |
+-------+---------+-------+-------+
| FL    | -629.50 |  0.30 |  0.86 |
+-------+---------+-------+-------+
| FF    | -630.36 |  0.35 |  0.00 |
+-------+---------+-------+-------+


The total ln evidence in favour of lensing is 16.38, which is considered strong evidence.
Thus this gamma-ray burst warrants further investigation.


.. figure:: ../../images/B_3770_YL2000__delmu_delt.png
    :figwidth: 80%
    :width: 100%
    :align: center
    :alt: BATSE trigger 3770 magnification ratio and time delay posterior histogram.

    BATSE trigger 3770 magnification ratio and time delay posterior histogram.


The magnification ratios and time delays are concordant, as one would expect from a gravitational lensing event.


.. figure:: ../../images/B_3770_YL2000__mass.png
    :figwidth: 80%
    :width: 100%
    :align: center
    :alt: BATSE trigger 3770 lens mass posterior histogram.

    BATSE trigger 3770 lens mass posterior histogram.

The inferred lens mass suggests an intermediate mass black hole as the gravitational deflector.
.. _3770ext:

.. role:: python(code)
   :language: python


GRB 950830 Nature script
------------------------


.. figure:: ../../images/trigger3770.PNG
   :figwidth: 80%
   :width: 100%
   :align: center
   :alt: BATSE trigger 3770 with FRED lens fit.

   BATSE trigger 3770 with FRED lens fit.

----

.. literalinclude:: ../../../../examples/gravlensnature.py
    :name: gravlensnature.py
    :caption: gravlensnature.py

----
.. _3770walk:

.. role:: python(code)
   :language: python


GRB 950830 code walkthrough
---------------------------

.. figure:: ../../images/trigger3770.PNG
   :figwidth: 80%
   :width: 100%
   :align: center
   :alt: BATSE trigger 3770 with FRED lens fit.

   BATSE trigger 3770 with FRED lens fit.


Walkthrough
^^^^^^^^^^^

We first attempt to fit all possible combinations of pulse types.
We find that in the strong majority of cases, a lensing model is preferred over a model of two individual pulses.



Imports:

.. code-block:: python

  import numpy as np

  from PyGRB.main.fitpulse import PulseFitter
  from PyGRB.backend.makemodels import create_model_from_key


Set up the :python:`PulseFitter` object:

.. code-block:: python

  sampler  = 'dynesty'
  nSamples = 2000
  GRB = PulseFitter(3770, times = (-.1, 1),
                datatype = 'tte', nSamples = nSamples, sampler = sampler,
                priors_pulse_start = -.1, priors_pulse_end = 0.6,
                priors_td_lo = 0,  priors_td_hi = 0.5)


Note that in some rare instances the code will choose to fit a small residual quite late to the noise, and the 'gravitationally lens' that noise.
So in the case where one is fitting a lens model, a better choice of start time priors might be :python:`priors_pulse_end = 0.3`.

Set up the models to fit to the light-curve.
:python:`lens_keys` are a single pulse, which is then 'gravitationally lensed' -- the light-curve is duplicated at the time delay and rescaled by the magnification ratio.
:python:`null_keys` are two pulses.

.. code-block:: python

  lens_keys = ['FL', 'FsL', 'XL', 'XsL']
  null_keys = ['FF', 'FsFs', 'XX', 'XsXs']
  keys = lens_keys + null_keys

  model_dict = {}
  for key in keys:
      model_dict[key] = create_model_from_key(key)
  models = [model for key, model in model_dict.items()]


Then run the model fitting code with

.. code-block:: python

  for model in models:
      GRB.main_multi_channel(channels = [0, 1, 2, 3], model = model)


Parallelisation
^^^^^^^^^^^^^^^

Parallelisation is done through the use of the indices parameter, which can be used to send one channel of one model to a different CPU in the cluster.
The difference would be changing the last two lines of code:

.. code-block:: python
  :emphasize-lines: 1-2,10-11

  # 8 models x 4 channels = 32 separate nested sampling runs.
  # This can be cheaply parallelised by sending each run to a different CPU
  lens_keys = ['FL', 'FsL', 'XL', 'XsL']
  null_keys = ['FF', 'FsFs', 'XX', 'XsXs']
  keys = lens_keys + null_keys

  model_dict = {}
  for key in keys:
      model_dict[key] = create_model_from_key(key)
  models = [model for key, model in model_dict.items()]

  indices = np.arange(32)
  GRB._split_array_job_to_4_channels(models = models, indices = indices, channels = [0, 1, 2, 3])


.. Example slurm (batch submission) scripts can be found in the dev-james branch of the repository.

After several days on a HPC, hopefully all the jobs will have finished.
The more parameters, the longer the job will take.
The :python:`'FsFs'` and :python:`'XsXs'` models have 18 and 20 parameters each, and may take up to a week depending on the chosen :python:`nSamples` (and other sampler keyword arguments).
Higher level analysis can be done with the :py:meth:`~PyGRB.postprocess.make_evidence_tables.get_evidence_from_models` method:

.. code-block:: python

  GRB.get_evidence_from_models(model_dict = model_dict)


To create plots of the gravitational lensing parameters, use the :python:`lens_calc` method:

.. code-block:: python

  for model in models:
      lens_bounds = [(0.37, 0.42), (0.60, 1.8)]
      GRB.lens_calc(model = model, lens_bounds = lens_bounds)


which will create the following two plots

.. figure:: ../../images/B_3770_YL2000__delmu_delt.png
    :figwidth: 80%
    :width: 100%
    :align: center
    :alt: BATSE trigger 3770 magnification ratio and time delay posterior histogram.

    BATSE trigger 3770 magnification ratio and time delay posterior histogram.


The magnification ratios and time delays are concordant, as one would expect from a gravitational lensing event.


.. figure:: ../../images/B_3770_YL2000__mass.png
    :figwidth: 80%
    :width: 100%
    :align: center
    :alt: BATSE trigger 3770 lens mass posterior histogram.

    BATSE trigger 3770 lens mass posterior histogram.

The inferred lens mass suggests an intermediate mass black hole as the gravitational deflector.
.. _973:

.. role:: python(code)
   :language: python


GRB 911031
----------

.. figure:: ../../images/BATSE_trigger_973_rates_rates.png
   :figwidth: 80%
   :width: 100%
   :align: center
   :alt: BATSE trigger 973.

   BATSE trigger 973.

----

.. literalinclude:: ../../../../examples/gravlens.py
   :name: gravlens.py
   :caption: gravlens.py

----

Channel 1 (red)

+-------+----------+-------+--------+
| Model |   ln Z   | error | ln BF  |
+-------+----------+-------+--------+
| FL    | -3852.93 |  0.60 |  0.00  |
+-------+----------+-------+--------+
| FF    | -3623.66 |  0.67 | 229.27 |
+-------+----------+-------+--------+


Channel 2 (orange)

+-------+----------+-------+--------+
| Model |   ln Z   | error | ln BF  |
+-------+----------+-------+--------+
| FL    | -3967.64 |  0.61 |  0.00  |
+-------+----------+-------+--------+
| FF    | -3665.85 |  0.68 | 301.79 |
+-------+----------+-------+--------+


Channel 3 (green)

+-------+----------+-------+--------+
| Model |   ln Z   | error | ln BF  |
+-------+----------+-------+--------+
| FL    | -4014.63 |  0.62 |  0.00  |
+-------+----------+-------+--------+
| FF    | -3640.90 |  0.70 | 373.73 |
+-------+----------+-------+--------+


Channel 4 (blue)

+-------+----------+-------+-------+
| Model |   ln Z   | error | ln BF |
+-------+----------+-------+-------+
| FL    | -2986.56 |  0.57 |  0.00 |
+-------+----------+-------+-------+
| FF    | -2970.48 |  0.62 | 16.08 |
+-------+----------+-------+-------+


We see that the Bayesian evidence strongly prefers a two pulse model.


Notes
^^^^^
If you want to fit two pulses, one FRED and one FRED-X, but you do not know the order in which you would like to fit them in (eg. the two pulses are overlapping), then you will need to do both a ['FX'] model and a ['XF'] model.
The best fitting model is then the one with the higher Bayesian evidence.
This is because the order of the pulses is automatically encoded in the priors.
We do this so that the results are unimodal for degenerate cases, eg ['FF'].
PyGRB.postprocess package
=========================

PyGRB.postprocess.abp module
----------------------------

.. automodule:: PyGRB.postprocess.abp
    :members:
    :undoc-members:
    :show-inheritance:

PyGRB.postprocess.make\_evidence\_tables module
-----------------------------------------------

.. automodule:: PyGRB.postprocess.make_evidence_tables
    :members:
    :undoc-members:
    :show-inheritance:

PyGRB.postprocess.plot\_analysis module
---------------------------------------

.. automodule:: PyGRB.postprocess.plot_analysis
    :members:
    :undoc-members:
    :show-inheritance:

PyGRB.postprocess.plot\_gl\_posteriors module
---------------------------------------------

.. automodule:: PyGRB.postprocess.plot_gl_posteriors
    :members:
    :undoc-members:
    :show-inheritance:
PyGRB.fetch package
===================

PyGRB.fetch.get\_BATSE module
-----------------------------

.. automodule:: PyGRB.fetch.get_BATSE
    :members:
    :undoc-members:
    :show-inheritance:
PyGRB.backend package
=====================

PyGRB.backend.admin module
--------------------------

.. automodule:: PyGRB.backend.admin
    :members:
    :undoc-members:
    :show-inheritance:

PyGRB.backend.makekeys module
-----------------------------

.. automodule:: PyGRB.backend.makekeys
    :members:
    :undoc-members:
    :show-inheritance:

PyGRB.backend.makemodels module
-------------------------------

.. automodule:: PyGRB.backend.makemodels
    :members:
    :undoc-members:
    :show-inheritance:

PyGRB.backend.makepriors module
-------------------------------

.. automodule:: PyGRB.backend.makepriors
    :members:
    :undoc-members:
    :show-inheritance:

PyGRB.backend.multipriors module
--------------------------------

.. automodule:: PyGRB.backend.multipriors
    :members:
    :undoc-members:
    :show-inheritance:

PyGRB.backend.rate\_functions module
------------------------------------

.. automodule:: PyGRB.backend.rate_functions
    :members:
    :undoc-members:
    :show-inheritance:

PyGRB.backend.rateclass module
------------------------------

.. automodule:: PyGRB.backend.rateclass
    :members:
    :undoc-members:
    :show-inheritance:
PyGRB.preprocess package
========================

PyGRB.preprocess.BATSEpreprocess module
---------------------------------------

.. automodule:: PyGRB.preprocess.BATSEpreprocess
    :members:
    :undoc-members:
    :show-inheritance:

PyGRB.preprocess.GRB\_class module
----------------------------------

.. automodule:: PyGRB.preprocess.GRB_class
    :members:
    :undoc-members:
    :show-inheritance:

PyGRB.preprocess.abstract module
--------------------------------

.. automodule:: PyGRB.preprocess.abstract
    :members:
    :undoc-members:
    :show-inheritance:

PyGRB.preprocess.grb module
---------------------------

.. automodule:: PyGRB.preprocess.grb
    :members:
    :undoc-members:
    :show-inheritance:

PyGRB.preprocess.simulated\_grb module
--------------------------------------

.. automodule:: PyGRB.preprocess.simulated_grb
    :members:
    :undoc-members:
    :show-inheritance:
PyGRB.main package
==================

PyGRB.main.fitpulse module
--------------------------

.. automodule:: PyGRB.main.fitpulse
    :members:
    :undoc-members:
    :show-inheritance:
