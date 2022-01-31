# About macrobot

Macrobot is an image analysis software for studying plant-pathogen interactions on macroscopic level. Currently the macrobot software can detect and quantify the following plant-pathogen interactions:
- Barley powdery mildew (Blumeria graminis f. sp hordei)
- Wheat powdery mildew (Blumeria graminis f. sp tritici)
- Wheat yellow (stripe) rust (Puccinia striiformis f.sp. tritici)
- Wheat brown (leaf) rust (P. graminis f. sp. tritici)

<img src="https://github.com/snowformatics/macrobot/blob/master/docs/images/Slide1.png" width="50%" height="50%"><br>
Figure 1: Powdery mildew on barley plants

The hardware system is based on a custom fully automated multispectral 2D imaging station (Figure 2).

<img src="https://github.com/snowformatics/macrobot/blob/master/docs/images/Bild8.png" width="50%" height="50%"><br>
Figure 2: Macrobot Module

See the macrobot hardware in action:
https://www.youtube.com/watch?v=SmoKQ_uMp34&t=56s

The entire pipline from image aquisition to image analysis is shown in Figure 3.

<img src="https://github.com/snowformatics/macrobot/blob/master/paper/figure.png" width="70%" height="70%"><br>
Figure 3: Software pipeline

# Citation

Lueck et al., (2020). BluVision Macro - a software for automated powdery mildew and rust disease quantification on detached leaves.. Journal of Open Source Software, 5(51), 2259, https://doi.org/10.21105/joss.02259

# Documentation
https://macrobot.readthedocs.io/en/latest/index.html


# Installation
Macrobot software was build and successfully tested on Windows operation system (Windows 7 and 10).

->Install Anaconda (https://www.anaconda.com/distribution/)

`conda create --name macrobot python=3.8`

`conda activate macrobot`

`conda install pip`

`pip install macrobot`


# Usage

1. Create a folder for the result. We will create a new folder on the desktop called mb_results.
2. Open the Ananconda prompt and activate your macrobot environment if you are not already there.<br/>`conda activate macrobot`<br/>
3. Macrobot is a command line program which requires the following arguments:
* source path (-s) - the path with the images coming from the Macrobot hardware system
* destination path (-d) - the path to store the results
* pathogen (-p) - which pathogen to predict ("mildew" or "rust")
4. For a test case we will use a test image set which will be automatically downloaded by the start of the software.
To tell the software to use the test images, we will enter "test_images" for the source path -s argument
5. Start the software with the following command for mildew (adapt the destination path):<br/>`mb -s test_images -d C:\Users\name\Desktop\mb_results\ -p mildew`<br/> or rust
<br/>`mb -s test_images -d C:\Users\name\Desktop\mb_results\ -p rust`<br/>
6. In your destination folder should appear all results:
* A csv file with the predicted values per leaf
* A report html file in folder report which allows and easy control over the pipeline.
* Images created by the software (white=pathogen, red=leaf detection, black=background)

If you want to use a real world experiments, make sure to provide the following folder structure with five images per plate (see documentation)

# Tests
cd to installation path and test folder e.g. d:\Anaconda\envs\mb_test\Lib\site-packages\macrobot\tests

Run pytest:

`pytest`


# Contributions:
We are strongly looking for contributions, some ideas how to support our software could be found here:
https://github.com/snowformatics/macrobot/wiki/Contributions

# References:
https://github.com/snowformatics/macrobot/wiki/References
---
title: 'BluVision Macro - a software for automated powdery mildew and rust disease quantification on detached leaves.'
tags:
  - plant phenotyping
  - powdery mildew
  - barley
  - wheat
  - rust
  - pucchinia
  - python
  - pathogen
authors:
 - name: Stefanie Lueck
   orcid: 0000-0003-0536-835X
   affiliation: 1
 - name: Ulrike Beukert
   orcid: 0000-0002-9482-3512
   affiliation: 2
 - name: Dimitar Douchkov
   orcid: 0000-0001-6603-4930
   affiliation: 1
   
 
affiliations:
 - name: Leibniz-Institut für Pflanzengenetik und Kulturpflanzenforschung Gatersleben, Stadt Seeland, Sachsen-Anhalt
   index: 1
 - name: Julius Kühn-Institut Quedlinburg, Sachsen-Anhalt
   index: 2
   
date: 20 April 2020
bibliography: paper.bib
---
 
# Summary

Powdery mildews and rusts are in the Top 10 of the major fungal pathogens in plant pathology [@2012_2; @2011_2; @2011_4]. The effect of powdery mildews on crop yields can amount to 40% of harvested grain [@2001]. Besides being an agriculturally important pathogen, the powdery mildews of wheat and barley are important models for studying the plant-pathogen interactions [@2014_2]. Wheat leaf rust and stripe rust are two other fungal pathogens, frequently causing epidemics with up to 70% yield losses [@2005; @2011_2] in combination with a decreased grain quality [@2012]. 

Crop protection against pathogens is mostly provided by the application of chemical agents [@2013], however, many of them have detrimental effects on non-target species [@2016].  Therefore, the trend is towards reducing the pesticide application and development of alternative and integrated protection methods [@2011_3]. The most sustainable method for crop protection appears to be the use of natural genetic resources and breeding for disease resistance [@2009_3]. Plant breeding is as old as the domestication of the first agricultural plants (>10000 years) [@2011_4]. The “Green revolution” [@2009] of the 1950-1960 and the more recent “omics” revolution (genomics, proteomics or metabolomics, etc.) [@2007] introduced many new approaches and technologies but the observation methods of the breeders remained mostly unchanged. Recently, new high-throughput observation methods for Phenomics [@2009_2] were introduced, thus providing the fundament for another breeding revolution. However, the powdery mildews and rusts are, as the majority of the plant pathogens, microscopic organisms in the initial and most critical stages of the infection process. Phenotyping of these early stages was significantly held back by the lack of technology for high-throughput phenotyping on a macroscopic and microscopic scale. 

To solve this problem, we have developed the BluVision Macro framework aimed to allow strictly quantitative assessment of disease and host responses on a macroscopic level. The system consists of a hardware part – the Macrobot [@2020] - a multimodal imaging station and robotized sample magazine/loader, and the BluVision Macro software, described in this article. The system is designed to work with samples placed in microtiter plates (MTP), which are well-established standard in biology and medicine. The loading of the MTPs to the imaging station and the image acquisition is fully automated. The system uses a 14-bit monochrome camera (Thorlabs 8050M-GE-TE) at a resolution of 3296×2472 px. The illumination is based on narrow bandwidth isotropic LED light sources (Metaphase Exolight-ISO-14-XXX-U) with 365nm (UV), 470nm (blue), 530nm (green) and 625nm (red) peak wavelength. For each plate monochrome images in all illumination  wavelengths are acquired separately and stored in 16-bit TIFF image files. Additionally, a background illumination image is taken and used for the separation of the foreground and background.

![Macrobot image acquisition and analysis workflow.The 14-bit monochrome camera of the Macrobot (A), acquires images in 5 different illumination scenarios (channels) for each plate (B). The red (625 nm), green (530 nm), and blue (470 nm) channels are combined into an RGB composite image (C), which is further segmented for detection of the infected area of the leaves (F). The UV (365 nm) and the backlight channel, are used for the leaf segmentation (D). Combined prediction for infected area and leaves serves to calculate the infected area as a percent of the leaf surface (E).\label{fig:example}](Macrobot_Figure1.png)




The image analysis pipeline currently contains software modules for powdery mildew of barley and wheat (*Blumeria graminis* f. sp. *hordei* resp. *tritici*), wheat stripe rust (*Puccinia striiformis* f.sp. *tritici*), wheat leaf rust (*P. triticina*), and will work without major modification for barley leaf, stripe rusts and probably other leaf diseases with a similar appearance. Phenotyping of other disease and non-disease related phenotypes will require development of dedicated modules. The system is running in production mode and generates phenotyping data for powdery mildew and rust disease resistance screens at the Leibniz Institute of Plant Genetics and Crop Plant Research (IPK) in Gatersleben, Germany, and the Julius Kühn-Institute (JKI) in Quedlinburg, Germany.


 
# Installation
The macrobot software can be installed with *Ananconda* and *pip* and was built and tested in the Microsoft Windows environment. <br>
It requires Python 3.8 or higher, numpy [@2011], scikit-image [@2014], opencv-python [@2000], pytest and jinja2. 

# Acknowledgment
This work was supported by grants from the German Federal Ministry of Research and Education (BMBF) – DPPN (FKZ 031A05) and GeneBank 2.0  (FKZ 031B0184)

# Additional information
This paper is dedicated in memory of Patrick Schweizer.

# References
macrobot package
================

Submodules
----------

macrobot.bgt module
-------------------

.. automodule:: macrobot.bgt
   :members:
   :undoc-members:
   :show-inheritance:

macrobot.helpers module
-----------------------

.. automodule:: macrobot.helpers
   :members:
   :undoc-members:
   :show-inheritance:

macrobot.mb\_pipeline module
----------------------------

.. automodule:: macrobot.mb_pipeline
   :members:
   :undoc-members:
   :show-inheritance:

macrobot.prediction module
--------------------------

.. automodule:: macrobot.prediction
   :members:
   :undoc-members:
   :show-inheritance:

macrobot.puccinia module
------------------------

.. automodule:: macrobot.puccinia
   :members:
   :undoc-members:
   :show-inheritance:

macrobot.cli module
-----------------------------

.. automodule:: macrobot.cli
   :members:
   :undoc-members:
   :show-inheritance:

macrobot.segmentation module
----------------------------

.. automodule:: macrobot.segmentation
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: macrobot
   :members:
   :undoc-members:
   :show-inheritance:
=======
Authors
=======

* Stefanie Lueck <lueck@ipk-gatersleben.de>
* Dimitar Douchkov <douchkov@ipk-gatersleben.de>
============
Installation
============

The macrobot software could be easily installed with pip. We recommend to use Anaconda for an easy installation process.

Macrobot software was build and successfully tested on Windows operation system (Windows 7 and 10).

Via Anaconda (recommended):
===========================

1.) Install Anaconda (https://www.anaconda.com/distribution/)

2.) Open Anaconda prompt and create a new environment:

``conda create --name macrobot python=3.8``

3.) Activate environment:

``conda activate macrobot``

4.) Install macrobot software:

``conda install pip``

``pip install macrobot``
================================
Macrobot image analysis pipeline
================================

Overview
========
Image analysis pipeline consists of 4 major steps:

1.) Creating a true 3-channel RGB image from the blue, green and red channel image

.. image:: images/plate_ex.png

2.) Frame segmentation

.. image:: images/frame.png

3.) Leaf segmentation

.. image:: images/leaf.png

4.) Pathogen prediction (white is pathogen, black background)

.. image:: images/predict.png


Analysis time
=============

Analysis time is approximately 5 seconds per plate (Intel(R) Core(TM) i7-9700 CPU @ 3.00GHz; 32 GB RAM; Windows 10 x64).





.. macrobot documentation master file, created by
   sphinx-quickstart on Fri Apr 24 13:33:04 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to macrobot's documentation!
====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   introduction.rst
   installation.rst
   use_case.rst
   macrobot_pipeline.rst
   pathogens.rst
   authors.rst
   license.rst
   changelog.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
============
Introduction
============

Macrobot Software Overview
==========================
Macrobot is a software and hardware framework for high-throughput image acquisition and analysis of microscopic and macroscopic images in plant pathology.

Deep investigation of the intimate details of the plant-pathogen interactions is essential to truly understand the defense mechanism of the plants and the evading strategies of the pathogens. By using this knowledge in plant breeding we may significantly diminish the enormous disease-related losses in agriculture by simultaneous reduction of application of potentially hazardous pesticides.

Currently the macrobot software can detect and quantify the following plant-pathogen interactions:

- Barley powdery mildew (Blumeria graminis f. sp hordei) on barely leaves
- Wheat powdery mildew (Blumeria graminis f. sp tritici) on wheat leaves (Figure 1)
- Wheat yellow rust (Puccinia graminis f.sp. tritici) on wheat leaves
- Wheat brown rust (Puccinia dispersa f. sp. tritici) on wheat leaves


.. image:: images/bgt.png

Figure 1: Powdery mildew on barley plants

Macrobot Image Acquisition
==========================

The hardware system is based on a custom fully automated multispectral 2D imaging station (Figure 2).

.. image:: images/macrobot.png

Figure 2: Macrobot Module

See the macrobot hardware in action: https://www.youtube.com/watch?v=SmoKQ_uMp34&t=56s

Macrobot Analysis Pipeline
==========================
Our image analysis pipeline (Figure 3) is aimed to detect macroscopic disease symptoms for barley and wheat powdery mildew as well as yellow and brown rust.

.. image:: images/figure.png

Figure 3: Software pipeline

===================
Supported pathogens
===================

Barley powdery mildew (Blumeria graminis f. sp hordei)

.. list-table::

    * - .. figure:: images/bgh_leaf.png

           Fig 1. Barley powdery mildew leaf segmentation

      - .. figure:: images/bgh_pred.png

           Fig 2. Barley powdery mildew pathogen prediction (7 days after inoculation)

Wheat powdery mildew (Blumeria graminis f. sp tritici)

.. list-table::

    * - .. figure:: images/bgt_leaf.png

           Fig 3. Wheat powdery mildew leaf segmentation

      - .. figure:: images/bgt_pred.png

           Fig 4. Wheat powdery mildew pathogen prediction (6 days after inoculation)

Wheat yellow (stripe) rust (Puccinia striiformis f.sp. tritici)

.. list-table::

    * - .. figure:: images/yell_leaf.png

           Fig 5. Wheat yellow rust leaf segmentation

      - .. figure:: images/yell_pred.png

           Fig 6. Wheat yellow rust pathogen prediction (15 days after inoculation)

Wheat brown (leaf) rust (P. graminis f. sp. tritici)

.. list-table::

    * - .. figure:: images/brown_leaf.png

           Fig 7. Wheat brown rust leaf segmentation

      - .. figure:: images/brown_pred.png

           Fig 8. Wheat brown rust pathogen prediction (6 days after inoculation)
================
Use case example
================


1. Create a folder for the result. We will create a new folder on the desktop called mb_results.
2. Open the Ananconda prompt and activate your macrobot environment if you are not already there.

``conda activate macrobot``

3. Macrobot is a command line program which requires the following arguments:

* source path (-s) - the path with the images coming from the Macrobot hardware system
* destination path (-d) - the path to store the results
* pathogen (-p) - which pathogen to predict ("mildew" or "rust")

4. For a test case we will use a test image set which will be automatically downloaded by the start of the software. To tell the software to use the test images, we will enter "test_images" for the source path -s argument

5. Start the software with the following command (adapt the destination path):

``mb -s test_images -d C:\Users\name\Desktop\mb_results\ -p mildew``

6. In your destination folder should appear all results:

* A csv file with the predicted values per leaf
* A report html file in folder report which allows and easy control over the pipeline.
* Images created by the software (white=pathogen, red=leaf detection, black=background)

If you want to use a real world experiments, make sure to provide the following folder structure with five images per plate:


| my_folder
| ├── experiment1
| │   ├── dai
|         └── plateID1
|             └── plateID1_backlight.tif
|             └── plateID1_green.tif
|             └── plateID1_blue.tif
|             └── plateID1_red.tif
|             └── plateID1_uvs.tif
| ├── experiment2
| │   ├── dai
|         └── plateID2
|             └── plateID2_backlight.tif
|             └── plateID2_green.tif
|             └── plateID2_blue.tif
|             └── plateID2_red.tif
|             └── plateID2_uvs.tif


Test image link:
https://doi.ipk-gatersleben.de/DOI/d92b5ec6-a83c-4ce6-99ab-ee4243764024/a9b2216b-4bb4-4300-8fc9-a259c485b236/2=========
Changelog
=========

Version 0.1
===========

- RGB image for frame detection
- Bgt and Bgh support

Version 0.2
===========

- Use UV image for frame detection

Version 0.3
===========
- Switched Functional programming -> OOP which allows easy integration of other pipelines
- Storing only final images
- 50x faster
- Added leaf and stripe rust pathogens

Version 0.4
===========
- Added a test image set feature
- Improve documentation
- Added integration test