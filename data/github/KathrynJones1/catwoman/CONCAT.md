---
title: 'catwoman: A transit modelling Python package for asymmetric light curves'
tags:
  - Python
  - astronomy
  - exoplanets
  - transit 
authors:
  - name: Kathryn Jones
    orcid: 0000-0002-2316-6850 
    affiliation: 1
  - name: Néstor Espinoza
    orcid: 0000-0001-9513-1449
    affiliation: 2
affiliations:
 - name: University of Bern, Center for Space and Habitability, Gesellschaftsstrasse 6, CH-3012, Bern, Switzerland
   index: 1
 - name: Space Telescope Science Institute, 3700 San Martin Drive, Baltimore, MD 21218, USA
   index: 2 
date: 11 June 2020
bibliography: paper.bib

aas-doi:
aas-journal: Astronomical Journal
---

# Summary

When exoplanets pass in front of their stars from our point of view on Earth, they imprint a transit signature on the stellar light curve which, to date, has been assumed to be symmetric in time, owing to the planet being modelled as a circular area occulting the stellar surface [see, e.g., @Mandel02; @Kreidberg15; @Luger19]. However this signature might be asymmetric due to several possible effects, one of which is the different temperature/pressure and/or chemical compositions the different terminator regions a transiting planet could have [see, e.g., @Powell19]. Being able to model these asymmetric signatures directly from transit light curves could give us an unprecedented glimpse into planetary 3-dimensional structure, helping constrain models of atmospheric evolution, structure and composition.

``catwoman`` is a Python package that models these asymmetric transit light curves, calculating light curves for any radially symmetric stellar limb darkening law and where planets are modelled as two semi-circles, of different radii, using the integration algorithm developed in [@Kreidberg15] and implemented in the ``batman`` library, from which ``catwoman`` builds upon. It is fast and efficient and open source with full documentation available to view at https://catwoman.readthedocs.io .
     
The light curves are modelled as follows: The decrease in flux, $\delta$, as a planet transits its star can be approximated by the sum 

\begin{eqnarray}
\label{eq:theproblem}
\delta = \sum_{i=1}^{N} I\left(x_m\right)\Delta A(x_m,R_{p,1},R_{p,2},\varphi,d),
\end{eqnarray}

splitting the semi-circles into iso-intensity bands centred on the star and for each intersectional segment (see \autoref{fig:strips}) you multiply its area, $\Delta A$, by the intensity of the star and then sum these strips to generate the full $\delta$ for a specific separation between the centre of the star and planet, $d$. The code then increments $d$ by a small pre-determined amount (based on the time array given by the user) and recalculates $\delta$.

![Diagram of the geometric configuration during transit of two stacked semi-circles (one of radius $R_{p,1}$, and another of radius $R_{p,2}$) that model the different limbs of an exoplanet transiting in front of a star. The area of the star has been divided in different sections of radius $x_i$ (dashed circles) --- between each subsequent section, the star is assumed to have a radially symmetric intensity profile (e.g., blue band between $x_{i-1}$ and $x_i$ above). In order to obtain the light curve, the challenge is to calculate the sum of the intersectional areas between a given iso-intensity band and the semi-circles, $\Delta A$ (blue band with dashed grey lines). Note the stacked semi-circles are inclined by an angle $\varphi$ with respect to the planetary orbital motion.\label{fig:strips}](strips.png)

The width of the iso-intensity bands determines the truncation error of the model. The model is first initialised with parameters including a maximum truncation error either set by the user or taken as the pre-set value as 1ppm. As in ``batman``, ``catwoman`` first calculates many models, with varying widths and geometrically searches for a width that produces an error less than 1% away (and always less than) the specified level. The model then uses this width value to calculate the desired light curves. A lower specified error, and therefore thinner iso-intensity bands, produces more accurate light curves, however more steps are needed to calculate $\delta$ which takes more time.  

``catwoman`` also allows for $\varphi$, the angle of rotation of the semi-circles, to vary as a free parameter, which is something no other model has tried to implement, accounting for the possibility of spin-orbit misalignments of the planet. The two semi-circle radii, $R_{p,1}$ and $R_{p,2}$, and other orbital variables are also completely free parameters.

``catwoman`` was designed to be used by astronomical researchers. For a realistic light curve with 100 in-transit data points, ``catwoman`` takes around 340 seconds to produce 1 million quadratic-limb-darkened light curves on a single 1.3 GHz Intel Core i5 processor. It is used in Espinoza & Jones (in prep.).

# Acknowledgements
We would like to thank the Max Plank Institute of Astronomy, Heidelberg, for providing the funding for this project and hosting Kathryn Jones as a summer student at the Institute. 

# References
.. :changelog:
1.0.12 (24-08-20)
~~~~~~~~~~~~~~~~~~~
-fixed a bug with the phi angle being corrected for orbital motion
-fixed a bug when, due to orbital motion, phi increases above 90 degrees or below -90 degrees

1.0.10 + 1.0.11 (12-08-20)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
-made sure numpy>=1.16.2 is installed when installing catwoman

1.0.9 (10-08-20)
~~~~~~~~~~~~~~~~~
-actually removed unused files (see below)

1.0.8 (10-08-20)
~~~~~~~~~~~~~~~~~
-removed unused files from directory so they aren't include in installation

1.0.7 (29-06-20)
~~~~~~~~~~~~~~~~~
-fixed a bug with installing catwoman straight from source archive, involved removing unused python files from import

1.0.6 (27-06-20)
~~~~~~~~~~~~~~~~~~
- fixed setup.py file so catwoman specifies and downloads dependencies automatically when pip installed in a new environment

catwoman: A transit modelling Python package for asymmetric light curves
==========================================================================
.. image:: https://github.com/KathrynJones1/catwoman/raw/master/docs/cw.png

``catwoman`` is a Python package that models asymmetric transit lightcurves where planets are modelled as two semi-circles with different radii in any orientation, for any radially symmetric stellar limb darkening law. 

``catwoman`` uses the integration algorithm developed for the ``batman`` package (Kreidberg 2015), from which ``catwoman`` builds upon. 

For a detailed introduction and more information please visit https://catwoman.readthedocs.io/.

``catwoman`` was peer reviewed by JOSS (Journal of Open Source Software) and Jones & Espinoza 2020 can be viewed here

.. image:: https://joss.theoj.org/papers/10.21105/joss.02382/status.svg
   :target: https://doi.org/10.21105/joss.02382

``catwoman`` has been demonstrated and used in `Espinoza & Jones 2021 <https://ui.adsabs.harvard.edu/abs/2021arXiv210615687E/abstract>`_ .


Installation
=============
You can install ``catwoman`` with pip (recommended):

::

	$ pip install catwoman


.. _quickstart:

Quickstart
============
This explains how to quickly and easily plot a catwoman transit using the quadratic limb darkening law. For a more detailed explanation of the parameters, inputs and possible outputs, see the :ref:`tutorial` tab.

``catwoman`` is a Python package that models asymmetric transit lightcurves where planets are modelled as two semi-circles. The key parameters involved in the asymmetry include ``params.rp`` and ``params.rp2`` which define the radius of each semi-circle and ``params.phi`` which is the angle of rotation of the top semi-circle defined from -90° to 90° like so: 

.. image:: phidiagram.png 


The first step is to import ``catwoman`` and the packages needed for it to run and to plot the results:

::
	
	import catwoman
	import numpy as np
	import matplotlib.pyplot as plt

Next, following a similar procedure as to that in ``batman``, initialise a ``TransitParams`` object to store the input parameters of the transit:

:: 

	params  = catwoman.TransitParams()
	params.t0 = 0. 				#time of inferior conjuction (in days)
	params.per = 1.				#orbital period (in days)
	params.rp = 0.1 			#top semi-circle radius (in units of stellar radii)
	params.rp2 = 0.1			#bottom semi-circle radius (in units of stellar radii)
	params.a = 15.				#semi-major axis (in units of stellar radii)
	params.inc = 90.			#orbital inclination (in degrees)
	params.ecc = 0. 			#eccentricity
	params.w = 90.				#longitude of periastron (in degrees)
	params.u = [0.1, 0.3]			#limb darkening coefficients [u1, u2]
	params.limb_dark = "quadratic" 		#limbs darkening model
	params.phi = 0.				#angle of rotation of top semi-circle (in degrees) 

Next make the time array to specify the times we want to calculate the model for:

::

	t = np.linspace(-0.05, 0.05, 1000)

Then, to initialise the model and calculate a light curve:

::
	
	model = catwoman.TransitModel(params,t) 	#initalises model
	flux = model.light_curve(params) 		#calculates light curve

To view the light curve: 

::
	
	plt.plot(t, flux)
	plt.xlabel("Time from central transit/days")
	plt.ylabel("Relative flux")
	plt.show()


.. image:: Simplesymmetric.png

To model an asymmetric planet, simply change ``params.rp`` and/or ``params.rp2`` and ``params.phi`` to change the orientation of the system.


Let's try this by re-initialising the parameters we want to change so that one of the semi-circles is 0.5% larger than the other and they are orientated with φ = 90°. There is no need to initialise the full model again here, whenever the light_curve function is run, it updates the parameters:

::

	params.rp = 0.1
	params.rp2 = 0.1005
	params.phi = 90.

Now we calculate the flux again for this new system:

:: 	

	flux2 = model.light_curve(params)

To view this new light curve:

:: 	

	plt.plot(t, flux2)
	plt.xlabel("Time from central transit/days")
	plt.ylabel("Relative flux")
	plt.show()

.. image:: Asymmetric.png

To clearly see the difference between this and the symmetric planet, we can plot the residuals as so:

:: 
	
	res = (flux2 - flux)*10**6
	plt.plot(t, res)
	plt.xlabel("Time from central transit/days")
	plt.ylabel("Difference in relative flux/ppm")
	plt.show()

.. image:: Asymmetric_diff.png


.. _installation:

Installation
============
pip
---
You can install ``catwoman`` with pip (recommended):

::

	$ pip install catwoman
.. _tutorial:
  
Tutorial
============

This tutorial explains in detail the different features of ``catwoman``.

Initialising the model
----------------------

As shown in the :ref:`quickstart`, to start setting up the model, one has to initialise a variety of parameters:
::

	import catwoman
	import numpy as np
	import matplotlib.pyplot as plt
	
	params = catwoman.TransitParams() 	#object to store transit parameters
	params.t0 = 0.				#time of inferior conjuction (in days)
	params.per = 1.				#orbital period (in days)
	params.rp = 0.1				#top semi-circle radius (in units of stellar radii)
	params.rp2 = 0.1005			#bottom semi-circle radius (in units of stellar radii)
	params.a = 15.				#semi-major axis (in units of stellar radii)
	params.inc = 90.			#orbital inclination (in degrees)
	params.ecc = 0.				#eccentricity
	params.w = 90.				#longitude of periastron (in degrees)
	params.u = [0.1, 0.3]                   #limb darkening coefficients [u1, u2]
	params.limb_dark = "quadratic"		#limb darkening model
	params.phi = 0.				#angle of rotation of top semi-circle

	time = np.linspace(-0.04, 0.04, 1000)	#array of times to calculate the light curves for
	model = catwoman.TransitModel(params, time, max_err = 0.1)	#initialises the model

As in ``batman``, the initialisation step automatically calculates the array of separation of centres between the star and the planet and also pre-runs the light_curve function numerous times in order to find the approriate integration step size for a given ``max_err``. 

``catwoman`` does this for all the supported limb darkening laws ("quadratic", "logarithmic", "exponential", "nonlinear", "linear", "power2", "uniform" and "custom").

*Note*: The default for ``max_err`` is 1ppm and describes the allowed error (in ppm) between the smallest integration step size and the selected integration step size. The lower the specified ``max_err``, the smaller the step size and the longer this initialisation step will take to run.

Calculating light curves
-----------------------------  

To calculate a light curve we run the ``light_curve`` function like so:
::
	
	flux = model.light_curve(params) 		#calculates light curve

This flux can now be plotted:
:: 
	
	plt.plot(time, flux)
	plt.xlabel("Time from central transit/days")
	plt.ylabel("Relative flux")
	plt.show()

.. image:: tutorialbasic.png
				  
Alternatively, if you wanted to change a parameter, you can do this by simply redefining the parameter of interest, say it is the ``params.rp``:
::

	params.rp = 0.09 			#top semi-circle radius (in units of stellar radii)

Now the new flux can be quickly calculated without having to re-initialise the model:
::

	flux2 = model.light_curve(params) 	#calculates light curve

To plot the two fluxes:
::

	plt.plot(time, flux)
	plt.plot(time, flux2)
        plt.xlabel("Time from central transit/days")
        plt.ylabel("Relative flux")
        plt.show()

.. image:: tutorial_newparam.png

This can be repeated for any ``params`` change. However if you want to change the ``time`` or ``max_err``, the model will need to be reinitialised as a new integration step size will need to be calculated.

This can make it easy to loop over certain parameter inputs and plot many light curves quickly. For example, we can make the light curves for a range of ``phi`` values like so:
::

	flux = np.zeros((7,len(time)))
	params.rp = 0.1
	params.rp2 = 0.15
	
	for i in range(0,7):
		params.phi = (i-3)*30			#updates angle of rotation
		flux[i] = model.light_curve(params)	#calculates light curve
		plt.plot(time,flux[i],label=str((i-3)*30)+'°')
	plt.xlabel("Time from central transit/days")
	plt.ylabel("Relative flux")
	plt.xlim(-0.015, 0.015)
	plt.legend()
	plt.show()

.. image:: tutorial_changephi.png

The residuals can also be easily plotted:
::

	for i in range(1,7):
        	plt.plot(time,flux[i]-flux[0],label='flux('+(str((i-3)*30)+'°) - flux(-90°)'))
		plt.xlabel("Time from central transit/days")
	plt.ylabel("Relative flux")
	plt.legend()
	plt.show()

.. image:: tutorial_phires2.png


Limb darkening functions
------------------------- 

As for ``batman``, ``catwoman`` allows you to choose one of the following limb darkening functions for the star:

.. math::

	\begin{align}
	  I(\mu) &= I_0                            						& &\text{(uniform)} 		\\
	  I(\mu) &= I_0[1 - c_1(1-\mu)]								& &\text{(linear)}		\\
	  I(\mu) &= I_0[1 - c_1(1 - \mu) - c_2(1-\mu)^2]	 				& &\text{(quadratic)}		\\
  	  I(\mu) &= I_0[1 - c_1(1 - \mu) - c_2(1-\sqrt{\mu})]                                   & &\text{(square-root)}         \\
  	  I(\mu) &= I_0[1 - c_1(1 - \mu) - c_2\mu\ln{\mu}]                                      & &\text{(logarithmic)}         \\
  	  I(\mu) &= I_0\left[1 - c_1(1 - \mu) - c_2/(1-\exp{\mu})\right]                  	& &\text{(exponential)}         \\
  	  I(\mu) &= I_0\left[1 - c_1(1 - \mu^{c_2})\right]                  	& &\text{(power2)}         \\
	  I(\mu) &= I_0[1 - c_1(1-\mu^{1/2}) - c_2(1- \mu) - c_3(1-\mu^{3/2}) - c_4(1-\mu^2)]  	& &\text{(nonlinear)}				
	\end{align}

where :math:`\mu = \sqrt{1-x^2}` where x is the normalised stellar radial coordinate defined between :math:`0 \leq x \leq 1` and :math:`I_O` is the normalisation constant for these laws so that integrated over the whole star, the total intensity is unity.
For each limb-darkening law you will need to provide the correct number of coefficients in order for the package to run.

Error tolerance
----------------
As mentioned in *Initialising the model*, the model contains a parameter called ``max_err``. If this is not specified, it will be set to the default ``max_err = 1.0``.

Whenever the model calculates a light curve from the parameters given it essentially splits up the planet into a series of very small strips of area in order to calculate the intensity of light that is blocked by the planet moving in front of the star at a particular time (see figure below). 

The width of these strips determines the accuracy of the light curve model and this is set by a scaling factor (``fac``). 

Once the model is initialised, internally, the program will calculate the light curve using an extremely small ``fac = 5e-4`` and an extremely large ``fac = 1`` and then find the error (or the largest difference) between their values.
If this is not equal to the ``max_err`` then the ``fac`` that produces an error within 1% of the ``max_err`` is found using a geometric search between the smallest and largest ``fac`` values.

As multiple light curves are being calculated during this step, this is the most time-intensive part of the package. However once the model has been initialised (and the appropriate ``fac`` value has been determined), as previously explained, this doesn't need to be repeated if some of the parameters are changed.

.. image:: strips.png     

Supersampling
---------------
As in ``batman``, for long exposure times there is the option of calculating the average value of the light curve model over the time of exposure of the samples. Set up the model including the additional parameters ``supersample_factor`` and `exp_time` (in days) like so:
::
	model = catwoman.TransitModel(params, time, supersample_factor = 5, exp_time = 0.001)

This will produce a model calculated by splitting up the samples into 5 sub-samples over the duration of the 0.001 day exposure. When a light curve is calculated, it will keep these sub-samples separate until the end where it will calculate the mean of these and reshape the light curve back to the original intended size, as specified by the ``time`` array.   

Parallelisation
----------------
As ``catwoman`` is built upon ``batman``, the library also inherits its support for OpenMP and OpenACC for CPU parallelisation and GPU acceleration, respectively. The former is active by default, but the latter is usually not. We refer users to the ``batman`` documentation to understand how to enable OpenACC on, e.g., NVIDIA GPUs.

 


.. catwoman documentation master file, created by
   sphinx-quickstart on Tue Sep 17 13:30:53 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: cw.png

Welcome to catwoman's documentation!
====================================

When exoplanets pass in front of their stars, they imprint a transit signature on the stellar light curve which to date has been assumed to be symmetric in time, owing to the planet being modelled as a circular area occulting the stellar surface. However, this signature might be asymmetric due to different temperature/pressure and/or chemical compositions in the different terminator regions of the transiting planet.

``catwoman`` is a Python package that allows to model these asymmetric transit lightcurves, calculating light curves for any radially symmetric stellar limb darkening law, and where planets are modelled as two semi-circles, of different radii, using the integration algorithm developed in Kreidberg (2015) and implemented in the ``batman`` library, from which ``catwoman`` builds upon.

Please cite `Jones & Espinoza 2020 <https://joss.theoj.org/papers/10.21105/joss.02382>`_ and `Espinoza & Jones 2021 <https://ui.adsabs.harvard.edu/abs/2021arXiv210615687E/abstract>`_ if you use ``catwoman`` in your research.   

If you find a bug or have any problems with catwoman, please `opening an issue <https://github.com/KathrynJones1/catwoman/issues>`_ on the project's GitHub and we will try to get back to you as soon as possible.  

Table of Contents
==================

.. toctree::
   :maxdepth: 2
	
   installation
   quickstart
   tutorial
   API
	

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
.. _api:
  
API
============
catwoman has two main classes: `TransitModel` and `TransitParams`.


.. module:: catwoman

.. autoclass:: catwoman.TransitParams
   :inherited-members:


.. autoclass:: catwoman.TransitModel
   :members: calc_err, light_curve

