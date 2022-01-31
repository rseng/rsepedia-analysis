<!-- Fill in the information below before opening an issue. -->

#### What needs to be created or improved?
<!-- Provide a clear and concise description of the issue. -->

#### Can you provide an example?
<!-- Provide a link or minimal code snippet that demonstrates the issue. -->
```python
import lightkurve
# optional: insert code here ...
```

#### What is the goal / expected behaviour?
<!-- Describe the behavior you expected and how it differs from the behavior observed in the example. -->



0.0.1 (2018-09-24)
++++++++++++++++++

- minimum viable product

1.0 (2020-10-16)
++++++++++++++++++

- first release including
    - query via lightkurve
    - light curve de-trending (with lightkurve and K2SC)
    - flare search
    - flare characterization
    - injection/recovery diagnostics
    - statistical analysis tools

1.0.1 (2021-02-20)
++++++++++++++++++

- fixed the lightcurve version to <2 
- altaipony.__version__ now returns a version number


2.0 (2021-05-11)
++++++++++++++++++

- now available
    - compatibility with lightcurve >2
    - docs subpage added for de-trending
    - experimental detrending function provided
    
2.0.1 (2021-06-29)
++++++++++++++++++    

- fixed links in the README
- JOSS version
|ci-badge| |docs-badge| |license-badge| |requirements-badge| |joss-badge| |zenodo-badge|


.. |joss-badge| image:: https://joss.theoj.org/papers/10.21105/joss.02845/status.svg
   :target: https://doi.org/10.21105/joss.02845

..  |zenodo-badge| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.5040830.svg
                    :target: https://doi.org/10.5281/zenodo.5040830
   
.. |ci-badge| image:: https://travis-ci.org/ekaterinailin/AltaiPony.svg?branch=master
              :target: https://travis-ci.org/ekaterinailin/AltaiPony

.. |docs-badge| image:: https://readthedocs.org/projects/altaipony/badge/?version=latest
	      :target: https://altaipony.readthedocs.io/en/latest/?badge=latest
	      :alt: Documentation Status
	      
.. |license-badge|  image:: https://img.shields.io/github/license/mashape/apistatus.svg   
		    :target: https://github.com/ekaterinailin/AltaiPony/blob/master/LICENSE 
		    :alt: GitHub	
.. |requirements-badge| image:: https://requires.io/github/ekaterinailin/AltaiPony/requirements.svg?branch=master
                       :target: https://requires.io/github/ekaterinailin/AltaiPony/requirements/?branch=master
                       :alt: Requirements Status


.. image:: logo.png
   :height: 100px
   :width: 100px
   :alt: Logo credit: Elizaveta Ilin, 2018

AltaiPony
=========

De-trend light curves from Kepler, K2, and TESS missions, and search them for flares. Inject and recover synthetic flares to account for de-trending and noise loss in flare energy and determine energy-dependent recovery probability for every flare candidate. Uses the ``K2SC`` and ``lightkurve`` under the cover, as well as ``pandas``, ``numpy``, ``pytest``, ``astropy`` and more.

Find the documentation at altaipony.readthedocs.io_

Installation
^^^^^^^^^^^^^

Use pip to install AltaiPony

>>> pip install altaipony


Or install directly from the repository:

>>> git clone https://github.com/ekaterinailin/AltaiPony.git
>>> cd AltaiPony
>>> python setup.py install



Getting Started
^^^^^^^^^^^^^^^^

See this notebook_ for an easy introduction, also docs_.


Problems?
^^^^^^^^^

 Often, when something does not work in **AltaiPony**, and this documentation is useless, troubleshooting can be done by diving into the extensive **lightkurve** docs_. Otherwise, you can always shoot Ekaterina an email_ or directly open an issue on GitHub_. Many foreseeable problems will be due to bugs in **AltaiPony** or bad instructions on this website.


Contribute to AltaiPony
^^^^^^^^^^^^^^^^^^^^^^^

**AltaiPony** is under active development on Github_. If you use **AltaiPony** in your research and find yourself missing a functionality, I recommend opening an issue on GitHub_ or shooting Ekaterina an email_. Please do either of the two before you open a pull request. This may save you a lot of development time.

How to cite this work
^^^^^^^^^^^^^^^^^^^^^

If you end up using this package for your science, please cite Ilin et al. (2021) [a]_ and Davenport (2016) [b]_.

Please also cite `lightkurve` as indicated in their docs [1]_. 

Depending on the methods you use, you may also want to cite 

  - Maschberger and Kroupa (2009) [2]_ (MMLE power law fit)
  - Wheatland (2004) [3]_ (MCMC power law fit)
  - Aigrain et al. (2016) [4]_ and their softwar [5]_ (K2SC de-trending)


.. [a] Ekaterina Ilin, Sarah J. Schmidt, Katja Poppenhäger, James R. A. Davenport, Martti H. Kristiansen, Mark Omohundro (2021). "Flares in Open Clusters with K2. II. Pleiades, Hyades, Praesepe, Ruprecht 147, and M67" Astronomy & Astrophysics, Volume 645, id.A42, 25 pp.  	https://doi.org/10.1051/0004-6361/202039198 

.. [b] James R. A. Davenport "The Kepler Catalog of Stellar Flares" The Astrophysical Journal, Volume 829, Issue 1, article id. 23, 12 pp. (2016). https://doi.org/10.3847/0004-637X/829/1/23

.. [1] https://docs.lightkurve.org/about/citing.html

.. [2] Thomas Maschberger, Pavel Kroupa, "Estimators for the exponent and upper limit, and goodness-of-fit tests for (truncated) power-law distributions" Monthly Notices of the Royal Astronomical Society, Volume 395, Issue 2, May 2009, Pages 931–942, https://doi.org/10.1111/j.1365-2966.2009.14577.x

.. [3] Wheatland, Michael S. "A Bayesian approach to solar flare prediction." The Astrophysical Journal 609.2 (2004): 1134. https://doi.org/10.1086/421261

.. [4] Aigrain, Suzanne; Parviainen, Hannu; Pope, Benjamin "K2SC: flexible systematics correction and detrending of K2 light curves using Gaussian process regression" Monthly Notices of the Royal Astronomical Society, Volume 459, Issue 3, p.2408-2419 https://doi.org/10.1093/mnras/stw706

.. [5] Aigrain, Suzanne; Parviainen, Hannu; Pope, Benjamin "K2SC: K2 Systematics Correction." Astrophysics Source Code Library, record ascl:1605.012 https://ui.adsabs.harvard.edu/abs/2016ascl.soft05012A/abstract


.. _Appaloosa: https://github.com/jradavenport/appaloosa/
.. _altaipony.readthedocs.io: https://altaipony.readthedocs.io/en/latest/
.. _notebook: https://github.com/ekaterinailin/AltaiPony/blob/master/notebooks/Getting_Started.ipynb
.. _docs: https://altaipony.readthedocs.io/en/latest/
.. _Github: https://github.com/ekaterinailin/AltaiPony/issues/new
.. _email: eilin@aip.de
Frequently Asked Questions
=======================================


What is detected as a flare and not a flare? 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A flare in a Kepler or TESS light curve is a series of data points that fullfils the flare definition criteria. In AltaiPony, flares are positive excursions from the de-trended light curve above a certain noise threshold. You can play around with a number of parameters - details are explained in the section on `Defining Flare Candidates`_.

What are the default settings for flare detection? 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Default settings for flare detection are explained in the section on `Defining Flare Candidates`_.

How does **AltaiPony** handle the ramp ups at the end of TESS orbits? 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**AltaiPony** currently does not explicitly handle these ramp ups, so they can cause false positive detections. Using PDCSAP_FLUX instead of SAP_FLUX flux avoid these problems for the most part. However, if your algorithm can deal with these, you can pass a custom de-trending function to 

::

    FlareLightCurve.detrend("custom", func=<your custom function>) 


``func`` should *take* a ``FlareLightCurve`` as an argument, and can have arbitrary numbers of keyword arguments. It should also *return* a ``FlareLightCurve``.


What is defined as a the start and stop time of a flare? 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The start and stop times of a flare candidate mark the first and last+1 timestamp that fullfil all flare definition criteria. Accordingly the start and stop cadences and indices can be used to mask flares. An example:

::

    import numpy as np
    [...]
    # Take table of flare candidates in a FlareLightCurve (``flc``) to find all indices:
    flareindices = [list(np.arange(row.istart, row.istop)) for i, row in flc.flares.iterrows()]
    # Flatten the list:
    l = [i for sublist in l for i in sublist]
    # Mask the flares:
    flc.detrended_flux[l] = np.nan

.. _Defining Flare Candidates: https://altaipony.readthedocs.io/en/latest/tutorials/altai.html#defining-flare-candidates

.. AltaiPony documentation master file, created by
   sphinx-quickstart on Wed Sep 26 15:22:31 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

======================================================
Flare science in Kepler, K2 and TESS light curves
======================================================

**AltaiPony** is a toolbox for statistical flare studies in photometric time series from Kepler, K2, and TESS, including flare search and characterization, injection/recovery diagnostics, and statistical analysis of flare frequency distributions along with extensive documentation and tutorials.

Jump to `Getting Started`_ to get an overview!

.. K2 light curves are best de-trended using **K2SC** [5]_. For TESS and Kepler we employ a Savitzky-Golay filter. You can use most methods from lightkurve_'s ``KeplerLightCurve`` and ``TessLightCurve`` classes. **AltaiPony** can be used to inject and recover synthetic flare signatures. Jump to the core class here_. **AltaiPony** features the flare candidate definition from **Appaloosa** [b]_. 

.. _user-docs:

.. toctree::
	:caption: AltaiPony
	:maxdepth: 2

	install
        tutorials/index
   	api/index
        faq
	
	
Problems?
^^^^^^^^^

 Often, when something does not work in **AltaiPony**, and this documentation is useless, troubleshooting can be done by diving into the extensive **lightkurve** docs_. Otherwise, you can always shoot Ekaterina an email_ or directly open an issue on GitHub_. Many foreseeable problems will be due to bugs in **AltaiPony** or bad instructions on this website.


Contribute to AltaiPony
^^^^^^^^^^^^^^^^^^^^^^^

**AltaiPony** is under active development on Github_. If you use **AltaiPony** in your research and find yourself missing a functionality, I recommend opening an issue on GitHub_ or shooting Ekaterina an email_. Please do either of the two before you open a pull request. This may save you a lot of development time.

How to cite this work
^^^^^^^^^^^^^^^^^^^^^

If you end up using this package for your science, please cite Ilin et al. (2021) [a]_ and Davenport (2016) [b]_.

Please also cite **lightkurve** as indicated in their docs [1]_. 

Depending on the methods you use, you may also want to cite 

  - Maschberger and Kroupa (2009) [2]_ (MMLE power law fit)
  - Wheatland (2004) [3]_ (MCMC power law fit)
  - Aigrain et al. (2016) [4]_ and their software [5]_ (**K2SC** de-trending)


.. [a] Ekaterina Ilin, Sarah J. Schmidt, Katja Poppenhäger, James R. A. Davenport, Martti H. Kristiansen, Mark Omohundro (2021). "Flares in Open Clusters with K2. II. Pleiades, Hyades, Praesepe, Ruprecht 147, and M67" Astronomy & Astrophysics, Volume 645, id.A42, 25 pp.  	https://doi.org/10.1051/0004-6361/202039198 

.. [b] James R. A. Davenport "The Kepler Catalog of Stellar Flares" The Astrophysical Journal, Volume 829, Issue 1, article id. 23, 12 pp. (2016). https://doi.org/10.3847/0004-637X/829/1/23

.. [1] https://docs.lightkurve.org/about/citing.html

.. [2] Thomas Maschberger, Pavel Kroupa, "Estimators for the exponent and upper limit, and goodness-of-fit tests for (truncated) power-law distributions" Monthly Notices of the Royal Astronomical Society, Volume 395, Issue 2, May 2009, Pages 931–942, https://doi.org/10.1111/j.1365-2966.2009.14577.x

.. [3] Wheatland, Michael S. "A Bayesian approach to solar flare prediction." The Astrophysical Journal 609.2 (2004): 1134. https://doi.org/10.1086/421261

.. [4] Aigrain, Suzanne; Parviainen, Hannu; Pope, Benjamin "K2SC: flexible systematics correction and detrending of K2 light curves using Gaussian process regression" Monthly Notices of the Royal Astronomical Society, Volume 459, Issue 3, p.2408-2419 https://doi.org/10.1093/mnras/stw706

.. [5] Aigrain, Suzanne; Parviainen, Hannu; Pope, Benjamin "K2SC: K2 Systematics Correction." Astrophysics Source Code Library, record ascl:1605.012 https://ui.adsabs.harvard.edu/abs/2016ascl.soft05012A/abstract


.. _k2sc: https://github.com/OxES/k2sc 
.. _Appaloosa: https://github.com/jradavenport/appaloosa
.. _lightkurve: https://github.com/KeplerGO/lightkurve
.. _docs: http://docs.lightkurve.org/index.html
.. _email: eilin@aip.de
.. _GitHub: https://github.com/ekaterinailin/AltaiPony
.. _here: https://altaipony.readthedocs.io/en/latest/api/altaipony.flarelc.FlareLightCurve.html
.. _Quickstart: https://altaipony.readthedocs.io/en/latest/install.html
.. _Getting Started: https://altaipony.readthedocs.io/en/latest/install.html#getting-started
Quickstart
=======================================

Installation
^^^^^^^^^^^^


Use pip to install AltaiPony

::
	
    pip install altaipony


Or install directly from the repository:

::
    
    git clone https://github.com/ekaterinailin/AltaiPony.git
    cd AltaiPony
    python setup.py install

This package depends, on `lightkurve`, `k2sc`, `numpy`, `pandas` and some other packages, most of which will be installed automatically. Have a look at `requirements.txt` in the repository to see a more extensive list.
   

Getting Started
^^^^^^^^^^^^^^^^

Working with Kepler and TESS light curves is similar, K2 light curves need extra attention, so we treat them separately. For each mission, there is a notebook that will guide you through the basic applications. We recommend taking a look at the `Finding Data`_ tutorial.

Kepler Light Curves
...................

If you are working with ``KeplerLightCurve`` objects, i.e. light curves from the original Kepler mission, try this_ notebook. This shows how to fetch a light curve from MAST, de-trend it with a Savitzky-Golay_ filter from scipy_, and find some flares.

TESS Light Curves
...................

If you are working with ``TessLightCurve`` objects, i.e. light curves from the TESS mission, you may want to try this other_ notebook instead. In this one, you can also test the injection recovery feature of **AltaiPony**, and obtain recovery probability and a corrected equivalent duration (*ED*, aka luminosity independent flare energy).

K2 Light Curves
...................

If you are working with ``KeplerLightCurve`` objects, i.e. light curves from the K2 mission, this last notebook_ is for you. It will run **k2sc**'s GP de-trending and find flares in a typical K2 long cadence light curve despite its heavy instrumental artifacts.


Next Steps
^^^^^^^^^^^

Define your own flare finding
.............................

Once you have tried basic **AltaiPony** on your light curves, you can start to adjust the flare finding parameters to your application, as explained in the `Finding Flares`_ tutorial.


Test the performance of your flare finding algorithm
.....................................................

You may then want to test the perfomance of your chosen flare finding setup by injecting and recoving synthetic flares into your light curves. **AltaiPony** provides a framework to do so, explained in the `Synthetic Flare Injection and Recovery`_ tutorial. Check out the `visualization`_ notebook for nice plots.

Analyze flare frequency distributions
......................................

For a statistical analysis of your flares, **AltaiPony** also features a set of tools for the analysis of flare frequency distributions, including visualization, and different methods for power law fitting. For starters, check out the tutorial on `Flare Frequency Distributions and Power Laws`_. If you want to go hands on, start with the `beginner`_ notebook. For more advanced applications, like working with samples of multiple stars and their flares, go to the `advanced`_ notebook. 


.. _Aigrain et al. 2016: http://ascl.net/1605.012
.. _fork: https://github.com/ekaterinailin/k2sc
.. _notebook: https://github.com/ekaterinailin/AltaiPony/blob/master/notebooks/01_Getting_Started.ipynb
.. _this: https://github.com/ekaterinailin/AltaiPony/blob/master/notebooks/03_Kepler_Light_Curves_With_Flares.ipynb
.. _Savitzky-Golay: http://www.statistics4u.info/fundstat_eng/cc_filter_savgolay.html
.. _scipy: https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.signal.savgol_filter.html
.. _other: https://github.com/ekaterinailin/AltaiPony/blob/master/notebooks/04_TESS_Light_Curves_With_Flares.ipynb
.. _in this tutorial: https://altaipony.readthedocs.io/en/latest/tutorials/altai.html
.. _Finding Flares: https://altaipony.readthedocs.io/en/latest/tutorials/altai.html
.. _Finding Data: https://altaipony.readthedocs.io/en/latest/tutorials/lcio.html
.. _Synthetic Flare Injection and Recovery: https://altaipony.readthedocs.io/en/latest/tutorials/fakeflares.html
.. _visualization: https://github.com/ekaterinailin/AltaiPony/blob/master/notebooks/05_Visualize_Injection_Recovery.ipynb
.. _beginner: https://github.com/ekaterinailin/AltaiPony/blob/master/notebooks/02_Beginner_Flare_Frequency_Distributions_and_Power_Laws.ipynb
.. _advanced: https://github.com/ekaterinailin/AltaiPony/blob/master/notebooks/06_Advanced_Flare_Frequency_Distributions_and_Power_Laws.ipynb
.. _Flare Frequency Distributions and Power Laws: https://altaipony.readthedocs.io/en/latest/tutorials/ffds.html
Finding Data
=====

There are several ways to create a ``FlareLightCurve`` in **AltaiPony**. For a quickstart, all you need is your target's TIC, EPIC or KIC (and campaign, quarter, or sector). Local TESS, Kepler and K2 light curves can be read with ``from_path``, but you can can also query MAST with ``fom_mast``. You can read in a light curve directly or convert a target pixel file. The latter is required if you want to use **k2sc** for de-trending. Familiarity with lightkurve_ is advantageous in any case.

Fetch a TESS light curve from MAST:

>>> from altaipony.lcio import from_mast
>>> flc = from_mast("TIC 29780677", mode="LC", mission="TESS", c=1)

Then, calling

>>> flc.plot()

will show

.. image:: ticplot.png
  :width: 600
  :alt: TESS flare light curve

Or a K2 target pixel file, and convert it to a ``FlareLightCurve`` using some extra arguments like ``flux_type``, ``cadence``, and ``aperture_mask`` :

>>> flc = from_mast("EPIC 211119999", mode="TPF", mission="K2", c=4, flux_type="PDCSAP_FLUX", cadence="long", aperture_mask="default")
>>> flc.plot()

.. image:: epicplot.png
  :width: 600
  :alt: K2 flare light curve

You can also fetch some remote file, de-trend it, store it as an AltaiPony light curve, and read it back in again, like so:

>>> flc = from_mast("TIC 358108509", mode="LC", mission="TESS", c=1)
>>> dflc = flc.detrend("savgol")
>>> dflc.to_fits("ponylc.fits")
>>> ...
>>> from altaipony.lcio import from_path
>>> rflc = from_path("ponylc.fits", mode="AltaiPony", mission="TESS")

.. module:: altaipony.lcio

.. automodapi:: altaipony.lcio
   :no-heading:
   :no-inherited-members:
   :no-inheritance-diagram:


.. _lightkurve: https://github.com/KeplerGO/lightkurve
 
   
Synthetic Flare Injection and Recovery
=====

To characterize how well the de-trending and flare finding procedures actually find and characterized flares in a light curve you can inject synthetic flares into it and run the procedures to compare the recovered events to the injected ones. In **AltaiPony**, the shapes of the flares are generated using the emprical flare model from `Davenport et al. (2014)`_.

To run a series of injections, call the ``sample_flare_recovery()`` method on a ``FlareLightCurve``.

::

    from altaipony.lcio import from_mast
    flc = from_mast("TIC 29780677", mode="LC", c=2, mission="TESS")
    flc = flc.detrend("savgol")
    flc, fake_flc = flc.sample_flare_recovery(inject_before_detrending=True, mode="savgol", 
                                              iterations=50, fakefreq=1, ampl=[1e-4, 0.5], 
                                              dur=[.001/6., 0.1/6.])

``flc`` is the original light curve with a new attribute ``fake_flares``, which is a DataFrame_ that includes the following columns:

* ``amplitude``: the synthetic flare's relative amplitude
* ``duration_d`` : the synthetic flare's duration in days [1]_
* ``ed_inj``: injected equivalent duration in seconds
* ``peak_time``: time at which the synthetic flare flux peaks 	
* all columns that appear in the `flares` attribute of `FlareLightCurve`, see here_. If the respective row has a value, the synthetic flare was recovered with some results, otherwise **AltaiPony** could not re-discover this flare at all.

``fake_flc`` is just like the original one, but without the ``fake_flares`` attribute. Instead, its flux contains synthetic flares from the last iteration run by ``sample_flare_recovery``. We return it because it is often useful to see what one actually injects:

::  

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12,5))
    fakeflc.plot(ax=ax, c="r", label="TIC 29780677 with injected flares")
    flc.plot(ax=ax, c="k")

.. image:: ticplotinjected.png
  :width: 700
  :alt: TESS flare light curve with injected flares


Visualizing results
--------------------

Let's pick GJ 1243, a famous flare star, inject and recover some flares and look at the results. We follow the same step as before, but inject about 1000 synthetic flares. The setup below produces the minimum number of 1 flare per light curve chunk per iteration, and the TESS light curve of GJ 1243 is split into two uniterrupted segments, therefore, we run 500 iterations. You can get faster results by increasing ``fakefreq``. 

::

    from altaipony.lcio import from_mast
    flc = from_mast("GJ 1243", mode="LC", c=15, mission="TESS")
    flc = flc.detrend("savgol")
    flc, fake_flc = flc.sample_flare_recovery(inject_before_detrending=True, mode="savgol", 
                                              iterations=500, fakefreq=.01, ampl=[1e-4, 0.5], 
                                              dur=[.001/6., 0.1/6.])


We can now look at what fraction of the injected equivalent duration of flares with different recovered amplitudes and durations is recovered:

::

    fig = flc.plot_ed_ratio_heatmap(flares_per_bin=.3)
    plt.title("GJ 1243")


.. image:: edratio.png
  :width: 500
  :alt: ED ratio of recovered flares in GJ 1243

Similarly, we can illustrate what fraction of flares with different injected amplitudes and full-width-at-half-maximum values (:math:`t_{1/2}` in `Davenport et al. (2014)`_) is recovered:

::

    fig = flc.plot_recovery_probability_heatmap(flares_per_bin=.3)
    plt.title("GJ 1243");


.. image:: recprob.png
  :width: 500
  :alt: recovery probability of synthetic flares in GJ 1243


Flare characterization
-----------------------

What can we do with all these synthetic flares? We can use them to characterize the flare candidates in the original light curve. To do this, call the ``characterize_flares`` method on your ``FlareLightCurve``:

::
  
   flc = flc.characterize_flares(ampl_bins=10, dur_bins=10)


This method will tile up your sample of fake flares into amplitude and duration bins twice. First, it will tile up the sample into a matrix based on the *recovered* amplitude and durations. Second, it will do the same with the *injected* properties, and so include also those injected flares that were not recovered. 

The first matrix can be used to map each flare candidate's recovered equivalent duration to a value that accounts for losses dealt to the ED by photometric noise, and introduced by the de-trending procedure (if you chose ``inject_before_detrending=True`` above). The typical injected amplitude and duration of flares in that tile of the matrix can then be used by the second matrix to derive the candidate's recovery probability from the ratio of lost to recovered injected flares.

The results from this mapping are stored in the ``flares`` attribute, which now contains the following additional columns in the table:


* ``dur``: ``= tstop - tstart``


* ``ed_ratio``: ratio of recovered ED to injected ED in the synthetic flares in the matrix tile that contains flares with measured properties that are most similar to the candidate flare.
* ``ed_ratio_count``: number of synthetic flares in the tile
* ``ed_ratio_std``: standard deviation of ED ratios in the tile
* ``ed_corr``: ``= rec_err / ed_ratio``
* ``ed_corr_err``: quadratically propagated uncertainties, including ``ed_rec_err`` and ``ed_ratio_std``


As in ``ed_ratio`` but with amplitude:


* ``amplitude_ratio``
* ``amplitude_ratio_count``
* ``amplitude_ratio_std``
* ``amplitude_corr``
* ``amplitude_corr_err`` : uncertainty propagated from ``amplitude_ratio_std``


As in ``amplitude_ratio`` but with duration in days:


* ``duration_ratio``
* ``duration_ratio_count``
* ``duration_ratio_std``
* ``duration_corr``
* ``duration_corr_err``


As in the columns but now for recovery probability:


* ``recovery_probability``: float between 0 and 1
* ``recovery_probability_count``
* ``recovery_probability_std``


"Properties" always refers to amplitude and duration or FWHM.

For a subset of these parameters, ``flc.flares`` could look like this:

.. image:: characterized.png
  :width: 700
  :alt: characterized flares


.. rubric:: Footnotes

.. [1] At the moment this is not a very meaningful quantity because the decay of the flare goes on to infitiny! We may define full width at 1% of the fluxe or something as an approximation but that is for later and I am getting distracted. But we need it to map between injected and recovered flares, that is why it's hanging around in that table.


.. _DataFrame: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html

.. _here: https://altaipony.readthedocs.io/en/latest/api/altaipony.flarelc.FlareLightCurve.html#altaipony.flarelc.FlareLightCurve
.. _Davenport et al. (2014): https://ui.adsabs.harvard.edu/abs/2014ApJ...797..122D/abstract
Flare Frequency Distributions and Power Laws
==============================================

*You can try out everything in this section in the tutorial notebook in the [Github] repository.*

Once you have found all the flares, you can compute statistical measures using your flare table. 

The `FFD` module allows you to compute Flare Frequency Distributions. You can use it to

- convert the flare table into a cumulative flare frequency distribution
- fit the power law exponent :math:`\alpha` and intercept :math:`\beta`, 
- plot the resulting function in the cumulative form,
- test if the power law assumption must be rejected, 
- test if the distribution is truncated at the high energy end,
- apply statistical corrections to the flare properties using the ``ed_corr``, ``recovery_probability`` attributes of the flares in the flare table that you may obtain from performing *injection and recovery of synthetic flares* with ``FlareLightCurve.characterize_flares()``.

Finally, if your flare table contains contributions from multiple stars that you think generate flares that can be described by the same power law but with different detection thresholds, you can use the `multiple_stars` keyword to account for this to a first order approximation. 

*Note that results from samples with less than 100-200 flares should be interpreted with caution.*

A simple flare sample
-----------------------------

In the simplest of all cases, there is one star that was observed for a certain time with high cadence and very low noise. For the resulting light curve we obtained a table of flare candidates, for instance, using ``FlareLightCurve.find_flares()``.

Now, you can directly use the ``FlareLightCurve.flares`` table, or any ``pandas.DataFrame`` where the recovered flare energies column is named ``ed_rec``.

Assume we have such a ``FlareLightCurve`` called ``flc`` with the required attribute ``flc.flares``, we can create a ``FFD`` object 

::

    from altaipony.ffd import FFD
    simple_ffd = FFD(f=flc.flares)

``simple_ffd.f`` contains the table with the energies of the flares. We can also specify the observing time it took to detect the event listed in the table:

::

    simple_ffd.tot_obs_time = 20.
    
The unit is up to you, and you should know which one you are using. If you do not specify ``tot_obs_time``, the FFD frequencies will instead be the number counts, i.e. ``simple_ffd.tot_obs_time=1.``.

Convert the flare table into a cumulative flare frequency distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first method you apply before doing anything else is ``FFD.ed_and_freq()``. It gives you the sorted array of energies, their corresponding frequencies, and number counts for each event with a certain energy, that is the cumulative flare frequency distribution:

::

    import matplotlib.pyplot as plt
    ed, freq, counts = simple_ffd.ed_and_freq()
    plt.figure(figsize=(8, 6))
    plt.scatter(ed, freq, c="k")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("ED [s]")
    plt.ylabel("cumulative number of flares per time");
    
    
.. image:: FFD.jpg
  :width: 400
  :alt: a simple FFD


Fit a power law to the FFD
-----------------------------
  
Next, let's fit a power law to this distribution. We can use a Modified Maximum Likelihood Estimator (MMLE) approach detailed in Maschberger and Kroupa (2009) [1]_ to find the slope :math:`\alpha` and then do a simple least squares fit to estimate the intercept :math:`\beta`:

::

    simple_ffd.fit_powerlaw("mmle")    

The results can be accessed with `simple_ffd.alpha`, `simple_ffd.alpha_err`, `simple_ffd.beta`, and `simple_ffd.beta_err`, respectively. The uncertainties on :math:`\beta` are bootstrapped.

Alternatively, we can the Bayesian flare prediction approach explained in Wheatland (2004) [2]_ to find :math:`\alpha` and :math:`\beta` using the MCMC method:

::

    simple_ffd.fit_powerlaw("mcmc")    

The results can be accessed with `simple_ffd.alpha`, `simple_ffd.alpha_up_err`, and `simple_ffd.alpha_low_err`; and `simple_ffd.beta`, `simple_ffd.beta_up_err`, and `simple_ffd.beta_low_err`, respectively. Upper and lower uncertainties represent the 16th and 84th percentiles of the marginalized posterior distributions for  :math:`\alpha` and :math:`\beta`.

The MCMC method is preferable because it fit both :math:`\beta` and :math:`\alpha` simultaneously, but the sampling may be too costly for a quick estimate of the parameters.

Plot the resulting function in the cumulative form
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use `plot_powerlaw` to plot the result on top of the FFD with the code snippet below:

::

    fig, ax = plt.subplots(1, figsize=(8,6))
    ax.scatter(ed, freq, c="k")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("ED [s]")
    ax.set_ylabel("cumulative number of flares per time")
    simple_ffd.plot_powerlaw(ax, c="r", label=fr'$\alpha=$-{simple_ffd.alpha:.1f}')
    plt.legend();


.. image:: powerlaw.jpg
  :width: 400
  :alt: a simple FFD with powerlaw

Statistical tests
------------------

Test if the power law assumption must be rejected
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The stabilised Kolmogorov-Smirnov statistic, suggested by Maschberger and Kroupa (2009) [1]_, tests if we must reject the power law hypothesis for our FFD. It is not meaningful in absolute terms. But whenever we compare FFDs and/or their power law fits with each other it gives us a better sense of the statistical robustness of a sample at different significance levels. 

For this hypothesis test, we must define a significance level, which is 5% per default. Above this limit we must reject the null-hypothesis. In our context, this is the hypothesis that the distribution follows the power law with the parameters we calculated.

::

    ffd.is_powerlaw(sig_level=0.05)


Test if the distribution is truncated at the high energy end
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An interesting question in flare statistics is whether or not there is a high energy limit seen in the FFD of any given star. It is hard to tell by eye, because the high-energy tail is sparsly populated with events, and log-log plots are deceptive. We may, however, ask, how small the highest observed energy can be to be consistent with an infinite power law distribution. ``FFD.is_powerlaw_truncated()`` performs this exceedance test, a left-sided hypothesis test suggested by Maschberger and Kroupa (2009) [1]_

For this, we generate a random sample of power law distributions and determine their maximum energies. These power law distributions have the same power law exponent, the same minimum detected energy and the same total number of events each. If a large fraction of the maximum energies in the random sample above the maximum detected energy it is more likely that the power law distribution is in fact truncated. As a default value we use percentile :math:`=2.5\%`:

.. image:: truncation.png
  :width: 550
  :alt: exceedance test FFD
  
Dealing with multi-star samples
--------------------------------
  
The above example and the more involved case of when your flare sample 

- stems from multiple light curves with different detection limits and/or
- was characterized using ``FlareLightCurve.characterize_flares``

is demonstrated in this_ notebook on Github.
  
.. rubric:: Footnotes

.. [Github] https://github.com/ekaterinailin/AltaiPony/blob/master/notebooks/Beginner_Flare_Frequency_Distributions_and_Power_Laws.ipynb 

.. [1] Thomas Maschberger, Pavel Kroupa, Estimators for the exponent and upper limit, and goodness-of-fit tests for (truncated) power-law distributions, Monthly Notices of the Royal Astronomical Society, Volume 395, Issue 2, May 2009, Pages 931–942, https://doi.org/10.1111/j.1365-2966.2009.14577.x

.. [2] Wheatland, M. S. "A Bayesian approach to solar flare prediction." The Astrophysical Journal 609.2 (2004): 1134. https://doi.org/10.1086/421261
  
  
  .. _this: https://github.com/ekaterinailin/AltaiPony/blob/master/notebooks/Advanced_Flare_Frequency_Distributions_and_Power_Laws.ipynb
De-trending Light Curves
========================

AltaiPony currently features three different de-trending aproaches that aim at removing astrophysical and instrumental trends in the time series while preserving the flares.

Call the de-trender on a ``FlareLightCurve`` with 

::
   
     detrended_flc = flc.detrend(<detrending method>, **kwargs)


For <detrending method> you can pass one the following:

- "k2sc" - if you wish to de-trend a K2 light curve, this is the best way to deal with, among other issues in K2, the Kepler spacecraft roll (see [1]_, [2]_ and [3]_ for details)
- "savgol" - applies a Savitky-Golay [4]_ filter to the light curve. This method is quick and gives good results in both Kepler and TESS light curves, but you may have to play with the "window_length" keyword to get optimal results
- "custom" - you can use your own de-trending pipeline. Pass a function that takes a ``FlareLightCurve`` and returns a FlareLightCurve to `flc.detrend("custom", func=<your custom de-trending function>)` [5]_. 




.. [1] Ekaterina Ilin, Sarah J. Schmidt, Katja Poppenhäger, James R. A. Davenport, Martti H. Kristiansen, Mark Omohundro (2021). "Flares in Open Clusters with K2. II. Pleiades, Hyades, Praesepe, Ruprecht 147, and M67" Astronomy & Astrophysics, Volume 645, id.A42, 25 pp.  	https://doi.org/10.1051/0004-6361/202039198 

.. [2] Aigrain, Suzanne; Parviainen, Hannu; Pope, Benjamin (2016). "K2SC: flexible systematics correction and detrending of K2 light curves using Gaussian process regression" Monthly Notices of the Royal Astronomical Society, Volume 459, Issue 3, p.2408-2419 https://doi.org/10.1093/mnras/stw706

.. [3] Aigrain, Suzanne; Parviainen, Hannu; Pope, Benjamin "K2SC: K2 Systematics Correction." Astrophysics Source Code Library, record ascl:1605.012 https://ui.adsabs.harvard.edu/abs/2016ascl.soft05012A/abstract

.. [4] Savitzky, Abraham. and Golay, M. J. E. (1964). "Smoothing and Differentiation of Data by Simplified Least Squares Procedures.". Analytical Chemistry, Volume 36, Issue 8 https://doi.org/10.1021/ac60214a047

.. [5] One such function, used to de-trend Kepler and TESS short (1 min, 2 min and 20s) cadence light curve is currently accessible via ``altaipony.customdetrend.custom_detrending``, and will be documented and tested in detail soon (Ilin+2021 in prep.). 
Finding Flares
=====

First you'll need a de-trended light curve (more details on de-trending here_.). Let's pick a TESS light curve:

::
   
     rawflc = from_mast("TIC 29780677", mode="LC", c=2, mission="TESS")

If you have a raw ``FlareLightCurve`` we can call ``rawflc``, for Kepler and TESS light curves use:

::

    flc = rawflc.detrend("savgol")

The following snippet shows the difference: The raw flux is still stored in ``flc.flux``, the de-trended flux is in ``flc.detrended_flux``

::
   
     import matplotlib.pyplot as plt
    plt.figure(figsize=(12,5))
    plt.plot(flc.time, flc.flux / np.nanmedian(flc.flux)+0.1, c="r", label="PDCSAP_FLUX")
    plt.plot(flc.time, flc.detrended_flux / np.nanmedian(flc.detrended_flux), "b", label="detrended flux")
    plt.xlabel("Time - 2457000 [BKJD days]")
    plt.ylabel(r"Flux [e$^-$s$^{-1}$]")
    plt.xlim(flc.time[0], flc.time[-1])
    plt.ylim(.95,1.30)
    plt.legend(loc=2,fontsize=13);

.. image:: ticplotdetrend.png
  :width: 600
  :alt: de-trended TESS flare light curve

**Note:** K2 is more difficult, and computationally intense, but doable with:

::
    
    flc = rawflc.detrend("k2sc")

Now you have a de-trended light curve ``flc``, and you can search it for flares:

::
    
    flc = flc.find_flares()

This will return the initial light curve with a new attribute - ``flares``, which is a DataFrame_ with the following columns:

* ``ampl_rec`` - recovered amplitude measured relative to the quiescent stellar flux
* ``ed_rec`` - recovered equivalent duration of the flare, that is, the are under the light curve with quiescent flux subtracted.
* ``ed_rec_err`` - the minimum uncertainty on equivalent duration derived from the uncertainty on the flare flux values (see `Davenport (2016)`_ for details, Eq. (2)).
* ``cstart, cstop, istart, istop, tstart, tstop`` - start and end of flare candidates in units of cadence, array index and actual time in days.
* ``dur`` - which is ``tstop-tstart``
* ``total_n_valid_data_points`` -  number of data points search for flare epochs in the de-trended light curve.

In our case, it should look like this:

.. image:: flaretable.png
  :width: 700
  :alt: de-trended TESS flare light curve flares

Basic flare definition
^^^^^^^^^^^^^^^^^^^^^^^^^

In ``FlareLightCurve.find_flares()``, the flare candidate definition follows the criteria in `Chang et al. (2015)`_ Eqn. (3) a-d. 

* Flare candidate data points must be positive excursions from the median quiescent flux value.
* The positive offset must be at least :math:`N_1` :math:`\sigma` above the local scatter of the light curve. If the local scatter is not given explicitly by the ``sigma`` keyword, ``FlareLightCurve.detrended_flux_err`` will be used instead, which is equal to PDCSAP_FLUX_ERR in Kepler and TESS light curves.
* The positive offset + ``FlareLightCurve.detrended_flux_err`` must be at least :math:`N_2` :math:`\sigma` above the local scatter.
* The number of consecutive data points fulfilling the above criteria must be at least :math:`N_3`.

You can pass :math:`N_{1,2,3}` and ``sigma`` explicitly like 

::

    FlareLightCurve.find_flares(N1=3, N2=2, N3=3, sigma=<local_scatter_array>)


The default settings are: ``N1=3``, ``N2=2``, ``N3=3``. ``sigma`` defaults to ``FlareLightCurve.detrended_flux_err``.  So, if you do not want to pass an array of local scatter values with the keyword argument ``sigma`` to ``find_flares()``, the :math:`N_2` specification  automatically becomes the more restrictive criterion. In this scenario, choosing ``N1=3`` and ``N2=2`` check for the same criterion.

**Note 1:** You can only apply the ``find_flares()`` method once to each de-trended light curve to avoid accidently listing flare candidates obtained with diffent sets of criteria in a single table. If you want to try different sets you have to create a copy of your ``FlareLightCurve`` for each set.

**Note 2:** Another argument to tinker with is the ``minsep`` keyword. The default is ``minsep=3``, meaning that candidate flare events within 3 data points of each other are combined into one.

Extended flare definition
^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to the basic flare definition above, which is sufficient for flare candidate detection, you may want to add the decay phase of a flare candidate to the flagged data points. This is useful if you are looking for a more precise energy estimate or a better mask that will cover the gradual tails of the flares. To set this up, you can extend teh previous command like

::

    FlareLightCurve.find_flares(N1=3, N2=2, N3=3, sigma=<local_scatter_array>, addtails=True, tailthreshdiff=<decrease in N1 and N2>)

If the `addtails` flag is set, datapoints will be added after the detected stop times of flare candidates if 

* they are positive outliers, 
* if they fulfill the N1 criterion but with N1 reduced by `tailthreshdiff`, and if
* if they fulfill the N2 criterion but with N2 reduced by `tailthreshdiff`.

Data points are added successively until a point no longer meets either one of these criteria.

.. _here: https://altaipony.readthedocs.io/en/latest/tutorials/lcio.html
.. _DataFrame: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html
.. _Davenport (2016): https://iopscience.iop.org/article/10.3847/0004-637X/829/1/23
.. _Chang et al. (2015): https://ui.adsabs.harvard.edu/abs/2015ApJ...814...35C/abstract
Tutorials
==================

.. toctree::
	:caption: Tutorials for a detailed flare study
	:maxdepth: 1

	lcio
        detrend
	altai
	fakeflares
	ffds

FFD
====

`FFD` stands for the Flare Frequency Distribution class in **AltaiPony**. It allows you to analyse the statistical properties of flare samples.

.. module:: altaipony.ffd


API Documentation
==================

.. toctree::
	:caption: Core Functions
	:maxdepth: 1

	flarelc
	ffd

FlareLightCurve
================

`FlareLightCurve`_ is the core class in **AltaiPony** for light curve de-trending, flare finding, and characterizing.

.. module:: altaipony.ffd

.. automodapi:: altaipony.flarelc
   :no-heading:
   :inherited-members:

