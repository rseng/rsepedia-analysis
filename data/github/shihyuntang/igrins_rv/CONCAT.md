[![status](https://joss.theoj.org/papers/37282917527e6c195d9dff80107388fd/status.svg)](https://joss.theoj.org/papers/37282917527e6c195d9dff80107388fd)
[![DOI](https://zenodo.org/badge/266670787.svg)](https://zenodo.org/badge/latestdoi/266670787)

# [IGRINS RV](https://github.com/shihyuntang/igrins_rv): A Radial Velocity Pipeline for IGRINS


``IGRINS RV`` is a ``python`` open source pipeline for extracting radial velocities (RVs) from spectra taken with the Immersion GRating INfrared Spectrometer (IGRINS). It uses a modified forward modeling technique that leverages telluric absorption lines as a common-path wavelength calibrator. ``IGRINS RV`` achieves an RV precision in the H and K bands of around 25-30 m/s for narrow-line stars, and it has successfully recovered the planet-induced RV signals of both HD 189733 and &tau; Boo A. Visit [Stahl et al. 2021](https://ui.adsabs.harvard.edu/abs/2021AJ....161..283S/abstract) to see the published paper.

If you have any questions, suggestions, or wish to report a bug, please let us know by either opening an issue or contacting us (asa.stahl@rice.edu or sytang@lowell.edu).
More on how to contribute can be found in the [FAQ](https://github.com/shihyuntang/igrins_rv/wiki/FAQ#q-how-do-i-report-bugs) page.

### Installation
A detailed installation guide can be found on the [GitHub wiki page](https://github.com/shihyuntang/igrins_rv/wiki/Installation).

### How to run
The best way to learn how IGRINS RV works is to first play with the example data, following steps in [Demo Run With Example Data of GJ281](https://github.com/shihyuntang/igrins_rv/wiki/Demo-Run-With-Example-Data-of-GJ281). While you are waiting for the example code to finish, read the detailed documentation on the [GitHub wiki page](https://github.com/shihyuntang/igrins_rv/wiki) so you know how to alter the example commands to fit your OWN science targets.

***
If you use ``IGRINS RV`` in your work, please go to the [FAQ page](https://github.com/shihyuntang/igrins_rv/wiki/FAQ#q-how-do-i-properly-cite-igrins-rv) to see how to properly cite ``IGRINS RV``.

Acknowledgements:\
Many thanks to Dr. Gregory Mace for helping improve the user experence with `IGRINS RV`!
---
title: "IGRINS RV: A Python Package for Precision Radial Velocities with Near-Infrared Spectra"
tags:
  - Python
  - astronomy
  - radial velocity
  - IGRINS
authors:
 - name: Shih-Yun Tang
   orcid: 0000-0003-4247-1401
   affiliation: "1, 2"
 - name: Asa G. Stahl
   orcid: 0000-0002-0848-6960
   affiliation: "3"
 - name: Christopher M. Johns-Krull
   affiliation: "3"
 - name: L. Prato
   orcid: 0000-0001-7998-226X
   affiliation: "1, 2"
 - name: Joe Llama
   orcid: 0000-0003-4450-0368
   affiliation: "1"
affiliations:
 - name: Lowell Observatory, 1400 W. Mars Hill Road, Flagstaff, AZ 86001, USA
   index: 1
 - name: Department of Astronomy and Planetary Sciences, Northern Arizona University, Flagstaff, AZ 86011, USA
   index: 2
 - name: Department of Physics and Astronomy, Rice University, 6100 Main Street, Houston, TX 77005, USA
   index: 3
aas-doi: 10.3847/1538-3881/abf5e7
aas-journal: Astronomical Journal
date:
bibliography: master.bib
---

# Summary

The relative radial velocity of a star with respect to the Sun can be calculated from its electromagnetic spectrum using the Doppler Effect. This line-of-sight motion, called the Radial Velocity (RV), is an essential tool for astrophysicists. RVs are not only used to detect and characterize exoplanets, but also play a key role in studies of binary stars, star clusters, and moving group member identification.

In the past decade, RVs have primarily been measured from spectra in the optical wavelength regime. This is partly because of advancements in detector technology, but also because of the paucity of Earth's atmospheric absorption features (telluric lines) in the optical. Yet for fainter, cooler, smaller stellar object like M-type stars (stars with mass less than half of the Sun), which emit more energy in the Near-Infrared (NIR), observations in the NIR can save a considerable amount of exposure time. Also, M-type stars are the most common type of star. This along with its size increases the detectability of Earth-like planets around them. Moreover, the stellar activity that can drive false positive exoplanet detections, e.g., star spots carried into view by stellar rotation, is shown to be less severe in the NIR compared to optical.


# Statement of need

Current RV pipelines and techniques that can deliver RV precision in tens of m/s (or better) in the NIR, e.g., ``PySHELL`` [@cale19], ``wobble`` [@bede19], ``SERVAL`` [@zech18] or the PCA-based cross-correlating method used for the SPIRou spectrograph [@mout20], all require instruments that are highly stabilized and have well-characterized wavelength solutions. For example, the iSHELL spectrograph can be equipped with the methane isotopologue gas cell, and the SPIRou and the CARMENES (NIR channel) spectrographs come with uranium-neon hollow-cathode lamps and stabilized Fabry-Perot etalons. The Immersion GRating INfrared Spectrometer (IGRINS) spectrograph [@yuk10; @park14; @mace16; @mace18], on the other hand, was not designed to be RV-stable and comes with no means of wavelength calibration accurate enough to achieve RVs precise to tens of m/s using existing techniques. A new approach to extract precision RVs is needed for an echelle spectrograph like IGRINS, which offers fertile ground for RV science with its high resolution (R ~ 45,000) and broad spectral grasp (the full H and K bands).

``IGRINS RV`` is a pipeline tailored for extracting precision RVs from spectra taken with IGRINS on different facilities. This pipeline is built on the forward-modeling methodology that was successfully applied to CSHELL and PHOENIX spectra [@croc12] that utilized telluric lines as a common-path wavelength calibrator. Compared to RVs obtained by cross-correlation with stellar templates adopted by past studies, ``IGRINS RV`` gives three times higher RV precision, about 25--30 m/s, around narrow-line stars in both H and K bands, shown by years of monitoring on two RV standard stars, GJ\ 281 and HD\ 26257. ``IGRINS RV`` also pushes this technique, using telluric lines as wavelength calibrator for RV calculations, to its limits as studies found the stability of the telluric lines is about 10--20 m/s [@seif08; @fig10]. Moreover, ``IGRINS RV`` is also tailored to take into account specific aspects of the IGRINS instrument, like the variations in spectral resolution across the detector and the year-long K band detector defocus.

``IGRINS RV`` has demonstrated its effectiveness in identifying orbiting companions by successfully recovering the planet-induced RV signals of HD\ 189733 and Tau\ Boo\ A. ``IGRINS RV`` lets users choose to obtain absolute RVs or relative RVs, depending on whether their priority is coarse RV characterization or more precise RV monitoring. The code extends the science capabilities of an already powerful spectrograph, which lacked a publicly available RV pipeline until now. It facilitates the detection and/or characterization of exoplanets, binary stars, star clusters, and moving group members, and it enables such studies to be done in a more precise and uniform way.

``IGRINS RV`` makes use of the ``astropy`` [@astropy2013; @astropy2018] on handing sky coordinates and barycentric velocity correction, ``scipy`` [@2020SciPy] and ``numpy`` [@2020NumPy] on mathmatical calculation, ``nlopt`` [@john08; @box1965new] on the optimization process, ``pandas`` [@reback2020pandas; @mckinney2010] on data management, and ``matplotlib`` [@Hunter2007] on plotting. We also used a part of code from `BMC` [@marcos2021] for peak detection. ``IGRINS RV`` requires that the ``igrins plp v2.2.0`` [@leegul17] and ``Telfit`` [@gull14] packages be pre-installed. Detailed documentation and tutorials can be found on the GitHub wiki page.


# Acknowledgements

Partial support for this work was provided by NASA Exoplanet Research Program grant 80-NSSC19K-0289 to L. Prato. CMJ would like to acknowledge partial support for this work from grants to Rice University provided by NASA (award 80-NSSC18K-0828) and the NSF (awards AST-2009197 and AST-1461918). We are grateful for the generous donations of John and Ginger Giovale, the BF Foundation, and others which made the IGRINS-DCT program possible. Additional funding was provided by the Mt. Cuba Astronomical Foundation and the Orr Family Foundation. IGRINS was developed under a collaboration between the University of Texas at Austin and the Korea Astronomy and Space Science Institute (KASI) with the financial support of the US National Science Foundation under grants AST-1229522 and AST-1702267 to the University of Texas at Austin, and of the Korean GMT Project at KASI.


# References
# IGRINS Radial Relocity Pipeline [IGRINS RV](https://github.com/shihyuntang/igrins_rv)

**This package is NOT yet ready for public use.**

``IGRINS RV`` if a ``python`` open source pipeline for extracting radial velocity (RV) for the IGRINS instrument. It's core technique use a modified forward modeling technique that utilize the telluric absorption lines as the common-path wavelength calibrator. With RV precision at K band around 25 m/s and H band 50 m/s, ``IGRINS RV`` had successfully retrieval of planet-induced RV signals for both HD 189733 and &tau; Boo A. Visit [arXiv]() to see the published paper.

``IGRINS RV`` is now under intense development. \
If you wish to join the internal review and testing, please contact us (asa.stahl@rice.edu or sytang@lowell.edu).

News:\
**2020/08/04: igrins_rv v0.9-beta.1 public beta is under internal review and testing!!**\
2020/06/23: igrins_rv v0.85 is under internal review.. will come to public soon!!

***
Please visit the [GitHub wiki page](https://github.com/shihyuntang/igrins_rv/wiki) for a more detailed documentation and tutorials.
Here only contains a short XXX

Setting up and running igrins_rv is easy if you have Conda installed.
NOTE: ``IGRINS RV`` has only be tested on Linux machine...

## Packages Installation (part 1)
**Installation with conda**\
Single command to setup environment with all pkg needed (using the with the environment.yml file):
```
conda env create
```

This will create an environment called ``igrins_rv``. You can use
```
conda info --envs
```
to check all available environments.

Use
```
conda (source) activate igrins_rv
```
to enter the environment.\
If you do it correctly, your command line will now looks like:
```
(igrins_rv) -->
```
**If you do not have conda**, the basic requirement for running ``IGRINS RV``\
is python.3.7, and the following packages/versions:
>    - python=3.7
>    - matplotlib=3.1.3
>    - numpy=1.18.1
>    - pandas=1.0.3
>    - scipy=1.4.1
>    - astropy=4.0
>    - multiprocess
>    - cython=0.29.15
>    - requests=2.23.0
>    - nlopt=2.6.1
>    - pip:
>      - pysynphot=0.9.14

## Packages installation (part 2)
**``IGRINS RV`` comes with the files for Telfit, but they need to be installed.**

Note that the original Telfit pkg is available on https://github.com/kgullikson88/Telluric-Fitter, by @kgullikson88.\
But the files that come with ``IGRINS RV`` are **modified**, such that one line in ``./src/TelluricFitter.py`` is different:
```
line 672: resolution=None, vac2air=self.air_wave)
           --> resolution=None, vac2air=False)
```
To install Telfit, First, check your mechine's ``gcc`` version by
```
gcc --version
```
**Telfit will not work with version later than v9.X.X**\
Not sure for v8.X.X, but v5.5.0 & v7.5.0 work fine.

Also, Telfit only recognize "gcc", so please make sure the new installed gcc **is callable** by
```
gcc
```
If you intend to use more the 6 threads to run this pipeline, we recommended you go into the ``setup.py`` file, and\
modify ``num_rundirs = `` to a number that is **THREE** times the threads you intend to use.

Now, enter the ``igrins_rv`` environment (within which Telfit must be installed)\
and cd to Telluric-Fitter-master, then run:
```
(igrins_rv) --> python setup.py build
(igrins_rv) --> python setup.py install
```

## Running ``IGRINS RV``

For a detailed understanding of how ``IGRINS RV`` works and how it should be run, it is HIGHLY recommended that the user reads the associated [paper]().
``IGRINS RV`` is divided into four main steps: 
[Setup](https://github.com/shihyuntang/igrins_rv#step-0---setup-), 
[Telluric Modelling](https://github.com/shihyuntang/igrins_rv#step-1---telluric-modelling-), 
[Initial Convergence](https://github.com/shihyuntang/igrins_rv#step-2---initial-convergence-), and 
[Analysis](https://github.com/shihyuntang/igrins_rv#step-3---analysis-).
Each is provided in the package as a separate file and # can be run from the command line with keywords specifying all relevant information.

**All steps are run as**
```
(igrins_rv) --> python main_stepX.py [target name] [keywords]
```
where X is the step number (0-4) and information on the keywords for the step can be found with
```
(igrins_rv) --> python main_stepX.py [target name] -h
```
(The -h command here means --help.)


### Step 0 - Setup: 
Collects and organizes the available data for the target spectra and their associated telluric (A0) standards. 

**If the user have** a copy of the IGRINS Observation MASTERLOG, Step 0 will automatically consult a master observing log file\
to retrieve the relevant information on filenames, coordinates, observing times, and conditions.
```
(igrins_rv) --> python main_step0.py [target name] -
```
**If the user do not have** a copy of the IGRINS Observation MASTERLOG, the output files it generates must be compiled by the user themselves. 


### Step 1 - Telluric Modelling: 
Defines the wavelength regions to be analyzed; generates a synthetic, high-resolution telluric template for use in later model fits on a night by night basis. 

### Step 2 - Initial Convergence: 
Required if the average RV of the target star is unknown to $>$ 5 km/s precision. Performs an abbreviated analysis of the target star observations in order to converge to coarsely accurate RVs, which will be used as starting points for the more precise analysis in the next step; simultaneously does the same for target star's vsini. Only the single most precise echelle region is used, and all separate exposures for a given observation are combined into one higher S/N spectrum before fitting occurs. 

### Step 3 - Analysis: 
Performs a full analysis of each target star observation to produce accurate and precise RVs. All the wavelength regions defined in Step 1 are used, and the code analyzes each exposure that is part of a given observation separately (this allows estimates of the RV uncertainties. In other words, if a star is observed with the ABBA beam pattern, Step 2 analyzes one spectrum made from the combination of the four exposures, while Step 3 analyzes each A and B separately. 

Unless the target vsini is already known to high accuracy, an initial run of Step 3 in which vsini is allowed to vary is required. This will provide an estimate of vsini, which will then be plugged in as a fixed value in the second run of Step 3. 

Even if vsini is already known, it is recommended to run Step 3 twice, as it allows the user to confirm that RVs between consecutive runs have converged within the calculated uncertainties. Additional runs may be necessary depending on the quality of the observed spectra and the vsini of the target star. 

All above three steps can be ran by:
```
(igrins_rv) --> python main_step1.py [-h][keywords]
```
The `-h` command will show the basic information of main_step1.py, and tell you what parameters you can put in.


