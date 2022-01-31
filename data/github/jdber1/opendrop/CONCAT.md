# Contributing to OpenDrop

Thank you for taking the time to contribute.

We welcome all contributions such as:

* Bug reports
* Feature proposals
* Code patches

OpenDrop is licensed under [GNU GPLv3](https://github.com/jdber1/opendrop/blob/master/LICENSE).

## Bug reports and feature proposals

Please submit any bug reports and feature proposals by opening a new [GitHub issue](https://github.com/jdber1/opendrop/issues). Try to do a brief search of existing issues to see if the problem has already been raised to avoid duplicates. Include any information you think is relevant for replicating a bug, we will follow up for more details if needed.

Any other queries or help can be asked by creating a new issue as well.

## Contributing code

Code contributions are accepted via pull requests. Before making large changes however, please create a new issue so that we can discuss any proposed changes. Try to keep modifications focused and avoid correcting formatting changes in irrelevant code, this will make it easier to see what has actually changed.

We currently aim to support Python versions 3.6 and above.

### Preparing a development environment

Make sure your system has GTK+ 3 installed.

1. Clone the repository.
2. Create a new Python virtual environment (e.g. in the '.venv/' subdirectory of the project root).
3. Install the app's Python dependencies with `pip install -r requirements.txt` and manually install the OpenCV Python bindings to the virtual environment. Also install the build tool [Scons](https://pypi.org/project/SCons/).
4. Add the project root to the Python path using a .pth [path configuration file](https://docs.python.org/3/library/site.html). This will let you test the app as you develop.

### Developing

OpenDrop makes use of some GLib compiled resource features. Whenever a .ui file or a file in the 'opendrop/assets/' directory is modified, the GLib resource bundle needs to be rebuilt. This is done by running `scons opendrop/data.gresource` in the project root.

Versions are generated from git tags.

### Code style

There is no stringent coding style in place. Mainly just follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) conventions and maintain a column width of around 110 (not strict).

Please include type annotations that are "as descriptive as possible" unless exact type-safety is overly arduous.
---
title: 'OpenDrop: Open-source software for pendant drop tensiometry & contact angle measurements'
tags:
  - Python
  - pendant drop
  - drop shape analysis
  - tensiometry
  - surface tension
  - interfacial tension
  - contact angle
authors:
  - name: Eugene Huang
    affiliation: 1
  - name: Adam Skoufis
    affiliation: 2
  - name: Terence Denning
    affiliation: 3
  - name: Jianzhong Qi
    orcid: 0000-0001-6501-9050
    affiliation: 3
  - name: Raymond R. Dagastine
    orcid: 0000-0002-2154-4846
    affiliation: 1
  - name: Rico F. Tabor
    orcid: 0000-0003-2926-0095
    affiliation: 2
  - name: Joseph D. Berry
    orcid: 0000-0002-0961-7782
    affiliation: 1
affiliations:
 - name: Department of Chemical Engineering, University of Melbourne, Parkville 3010, Australia
   index: 1
 - name: School of Chemistry, Monash University, Clayton 3800, Australia
   index: 2
 - name: Computing and Information Systems, University of Melbourne, Parkville 3010, Australia
   index: 3
date: 3 February 2020
bibliography: paper.bib
---

# Summary
 Systems where two or more fluids exist in discrete phases
are ubiquitous in nature and in many manufacturing processes. The
common surface (or interface) between two fluids that do not mix
exists in a state of tension, an intrinsic property known as
interfacial tension. The contact angle is another fundamental property of interest when the
interface between two fluids is also in contact with a surface, for
example a water drop resting on a leaf. The contact angle is dependent on the surface energy of the solid and 
describes how liquids spread on a surface – vital information for
dynamic liquid-solid processes such as coating and painting.

Accurate measurements of interfacial tension allow researchers in
industry and academia to make deductions regarding the chemical
composition and crucially, the behavior of the interfaces, enabling
optimal design of devices and processes. In many real formulations or
applied systems, this basic but critical parameter can be quite
challenging to accurately measure. In addition, precise measurements
of the contact angle between a fluid-fluid interface and a solid
surface are critical in order to deduce wetting and spreading
characteristics of liquids on surfaces, and to calculate the surface energy of a solid by measuring the contact angle of a series of liquids on one type of surface. These surface properties are important when considering, to
name two examples, the application of paints to surfaces and
pesticides to plants. It is therefore clear that accurate, rapid and
reproducible measurements of interfacial tension and contact
angle are imperative for effective design, implementation and
optimization of processes involving multiphase systems.

The experimental apparatus required for measurements of interfacial tension and contact angle is conceptually extremely simple, requiring only a needle, a camera, and a light source. The complexity (and associated cost of commercial instruments) comes from the image processing and the complicated numerical algorithm required to calculate these quantities from the acquired experimental image. In 2015, we released OpenDrop, which enables interfacial tension measurements more rapidly, cheaply and accurately than commercial options [@Berry2015]. The only cost to the user is the camera required (approx. $20 - $1K depending upon application), whereas commercial instruments are much more expensive (~$50K). 

Here we present the latest version of OpenDrop. The new version, Barracuda, is able to measure interfacial tension and also contact angle in a variety of configurations with field-leading accuracy and reproducibility. The performance of OpenDrop compared to currently available commercial instrumentation is shown in \autoref{fig:ift} and \autoref{fig:ca}.
The contact angle measurement capability is new for this release, but has been used successfully in previous studies [@Prathapan2017]. OpenDrop has been written in Python because it is open-source, free, runs on multiple operating systems (including Linux, Mac OSX and Windows), and is easily integrable with a large number of mature, 3rd party open source libraries. In particular, OpenDrop utilises the sophisticated image processing capabilities of the OpenCV library in order to extract drop profiles from experimental images for input into the requisite numerical algorithm. Further, the ease of readability and modular nature of Python encourages and supports collaboration, and gives OpenDrop significant pedagogic value. Python can also be easily integrated with other languages, of particular importance to pendant drop tensiometry and contact angle measurements where integration of code needed to control cameras and associated software is a critical requirement. The previous version is in use in many research groups around the world, and is also used in teaching laboratories including Monash University. 

The availability of the software allows the interested user to
effectively implement, explore and further develop the techniques for
both research and teaching at a small fraction of the cost of
commercial options. 

![Comparison of the surface or interfacial tension of different systems calculated with OpenDrop against values reported in the literature.\label{fig:ift}](iftFigure.pdf){ width=60% }

<!-- ![Comparison of the surface or interfacial tension of different systems calculated with OpenDrop against values reported in the literature.\label{fig:ift}](iftFigure.pdf)-->
![Comparison of contact angles calculated in OpenDrop from experimental images in the literature against values calculated with commercial instrumentation.The images used are taken from [@Nie2017], [@Stacy2009] and [@Brown2016]. \label{fig:ca}](conAnFigure.pdf){ width=60% }


<!-- Consequently, OpenDrop will make significant impact
in both research and education by providing inexpensive access to
high-fidelity information on the stability, function, and behaviour of
interfaces, via a simple and user-friendly interface, with open-source
software that will enable users to implement their own functionality. -->




# Acknowledgements



# References
|logo|
======

.. |logo| raw:: html

    <img src="https://opendrop.readthedocs.io/en/latest/_images/opendrop_logo_wide.png" width="185px" alt="Logo">

.. START 

.. image:: https://joss.theoj.org/papers/10.21105/joss.02604/status.svg
    :target: https://doi.org/10.21105/joss.02604

OpenDrop is a fully-featured image analysis software for performing pendant drop tensiometry and contact angle measurements. Images can be loaded from the file system or acquired directly from USB webcams or GenICam (GigE |nbsp| Vision, USB3 |nbsp| Vision) compliant industrial cameras.

.. raw:: html

    <img src="https://raw.githubusercontent.com/jdber1/opendrop/master/screenshots/ift-prepare.png" width="700px" alt="OpenDrop screenshot">

The software is released under the **GNU GPL** open source license, and available for free.

.. |nbsp| unicode:: 0xA0
    :trim:

*********
Installation & User Guide
*********
OpenDrop is currently distributed as a Python package and can be installed on most operating systems (Windows, macOS, Linux).

For installation instructions and user guides, visit: https://opendrop.readthedocs.io/

Example images have been provided in the 'example_images' folder.

*********
Support
*********
For any questions, issues, or feedback feel free to `open an issue on the GitHub repo <https://github.com/jdber1/opendrop/issues>`_.

We can also be contacted by email `here <mailto:opendrop.dev@gmail.com>`_.
.. title:: Overview

.. toctree::
    :hidden:
    :maxdepth: 2

    installation/index
    usage/index
    developers/index

.. image:: images/opendrop_logo_wide.png
    :width: 370
    :align: center

|

.. include:: ../README.rst
    :start-after: .. START
    :end-before: *********

For installation instructions, see ":doc:`installation/index`".

----

| **Git repo:**
| https://github.com/jdber1/opendrop/

| **Questions, issues, or feedback:**
| https://github.com/jdber1/opendrop/issues
###############
Developer notes
###############

Stub.
############
Installation
############

**************
Release builds
**************

Stand-alone builds for Windows are provided for certain major releases and do not require the installation of
additional software: https://github.com/jdber1/opendrop/releases/.

Releases for Linux and macOS don't exist yet and OpenDrop should instead be installed as a Python package. See next section.

****************************
Building package from source
****************************

OpenDrop requires Python 3.6 or higher, the GTK 3 library, OpenCV Python bindings, and the following build dependencies:

    * Boost.Math
    * SUNDIALS ARKODE

Other required Python packages will be automatically installed by pip.

Platform specific build instructions follow.


Ubuntu
======

#. Install OpenCV.

   * If on Ubuntu 17.10 (or later)::

       sudo apt install python3-opencv

   * Alternatively there is an unofficial opencv-python_ package that can be installed using pip::
       
       pip3 install opencv-python


#. Install SUNDIALS. Unfortunately ``libsundials-dev`` from the Ubuntu repositories are too old, we require at least version 4.0.0 and above. Here are brief instructions for installing SUNDIALS from source.

   #. Download the latest version from the `releases page <https://computing.llnl.gov/projects/sundials/sundials-software>`_. (Note: the latest version requires a CMake version newer than available in Ubuntu < 20.04. If this affects you, try an older version of SUNDIALS like 4.0.0 instead.)

   #. Extract and change into the source directory, e.g.::

       tar -xvf sundials-5.7.0.tar.gz
       cd sundials-5.7.0/
   
   #. Create a build directory::

       mkdir build
       cd build/

   #. Configure, build, and install (make sure ``cmake`` and ``build-essential`` are installed from the Ubuntu repos)::

       cmake \
         -DEXAMPLES_INSTALL=OFF \
         -DBUILD_ARKODE=ON \
         -DBUILD_CVODE=OFF \
         -DBUILD_CVODES=OFF \
         -DBUILD_IDA=OFF \
         -DBUILD_IDAS=OFF \
         -DBUILD_KINSOL=OFF \
         -DBUILD_STATIC_LIBS=OFF \
         -DCMAKE_BUILD_TYPE=Release \
         ..
       make
       sudo make install

#. Install Boost.Math. If on Ubuntu 20.04 or newer, run::

       sudo apt install libboost-dev
   
   The ``libboost-dev`` package on older versions of Ubuntu is not recent enough and Boost will need to be
   installed from source. We need at least Boost 1.71.0.

#. Follow the `installation instructions here <https://pygobject.readthedocs.io/en/latest/getting_started.html#ubuntu-logo-ubuntu-debian-logo-debian>`_ for installing PyGObject and GTK.

#. Use pip to install OpenDrop from the repo::

       pip3 install git+https://github.com/jdber1/opendrop.git

   Run ``pip3 uninstall opendrop`` to uninstall.

#. Run ``python3 -m opendrop`` to launch the app.


macOS
=====

1. Install the latest version of Python 3 and pip. You can do so using a third-party package manager like MacPorts_ or Homebrew_.

2. - Install the unofficial opencv-python_ package by running::

         pip install opencv-python

     (Make sure ``pip`` refers to your Python 3's pip installation.)
   - Alternatively, OpenCV and its python bindings can also be installed using the `opencv Homebrew formula <https://formulae.brew.sh/formula/opencv>`_ or `opencv MacPorts port <https://www.macports.org/ports.php?by=library&substr=opencv>`_.

3. - If Homebrew was used to install Python 3, PyGObject and GTK can also be installed by running::

         brew install pygobject3 gtk+3

   - or if MacPorts was used, run::

         sudo port install py36-gobject3 gtk3

     (Instead of the ``py36-`` prefix, use ``py37-`` or ``py38-`` if Python 3.7/3.8 is the version installed.)

#. Install Boost.Math and SUNDIALS. (todo: Add MacPorts and Homebrew example).

4. Use pip to install OpenDrop from the repo::

       pip install git+https://github.com/jdber1/opendrop.git

   Run ``pip uninstall opendrop`` to uninstall.

5. Run ``python3 -m opendrop`` to launch the app.


Windows
=======

Installing OpenDrop as a Python package is possible on Windows using platforms like MSYS2 or Anaconda.  
The process is not very straightforward so your mileage may vary.


.. _opencv-python: https://pypi.org/project/opencv-python/
.. _MacPorts: https://www.macports.org/
.. _Homebrew: https://brew.sh/
Notes
=====

User input validation is not yet implemented, invalid user input may cause OpenDrop to crash or print errors to the console.
GenICam integration
===================

Install a GenTL producer, (e.g. see `harvesters README <https://github.com/genicam/harvesters#installing-a-gentl-producer>`_).

OpenDrop checks the environment variable GENICAM_GENTL64_PATH (specified by the GenTL standard) for GenTL producers. To verify that a GenTL producer is installed correctly, you can run::

    $ echo $GENICAM_GENTL64_PATH
    /opt/mvIMPACT_Acquire/lib/x86_64

(todo: Add details.)
Contact Angle
=============

A wizard-style window will guide you through the process of performing a contact angle analysis.

Image acquisition
-----------------

The contact angle image acquisition page is the same as the one for interfacial tension analyses.

Image processing
----------------

.. image:: images/conan_image_processing.png
    :alt: Contact Angle Image processing
    :width: 914
    :align: center

The image processing window requires you to define the 'drop region' and 'surface line' of the image. Click on the 'Drop region' button in the 'Tools' panel then drag over the image preview to define the region. Similarly, click on the 'Surface line' button and drag a line to define the surface that the drop is sitting on. With the 'Surface line' button depressed and the preview widget focused, use the arrow keys for finer adjustments of the surface line.

.. image:: images/conan_image_processing_defined.png
    :alt: Contact Angle Image processing, regions defined
    :width: 914
    :align: center

Once the drop region is defined, a blue outline will be drawn over the preview showing the drop profile that has been extracted.

The intersection angle between the drop profile and the surface line will be the contact angle measured.

In a contact angle analysis, OpenDrop uses image thresholding to separate the foreground from the background. Click on the 'Foreground detection' button to open a dialog bubble which will allow you to adjust the threshold value. A blue overlay is painted over parts of the image deemed to be in the foreground.

Click on 'Start analysis' to begin analysing the input images, or begin capturing and analysing images if using a camera.

Results
-------

.. image:: images/conan_results.png
    :alt: Contact Angle Results
    :width: 914
    :align: center

The results page for a contact angle analysis is quite simple.

A summary table is shown on the bottom half with a results visualizer on the top half. Graphs of the left and right contact angles are also available if more than one image is analysed.


Saving
------

.. image:: images/conan_save_dialog.png
    :alt: Contact Angle Save dialog
    :width: 828
    :align: center

Once an analysis is finished, click on the 'Save' button in the footer to open the save dialog. All data will be saved in a folder with name determined by the 'Name' entry, and in a parent directory determined by the 'Parent' selection. 

As a convenience, you may choose to save some pre-made plots.

.. image:: images/conan_save_output.png
    :alt: Contact Angle Example save output
    :width: 619
    :align: center

An example save output is shown above, and screenshots of the contents of some files are shown below. (All coordinates are with respect to the origin being on the top-left corner of the image with increasing x and y in the right and down directions respectively.)

.. figure:: images/conan_timeline_csv.png
    :alt: Contact Angle timeline.csv screenshot
    :width: 927
    :align: center

    timeline.csv

.. figure:: images/conan_profile_extracted_csv.png
    :alt: Contact Angle profile_extracted.csv screenshot
    :width: 159
    :align: center

    drop1/profile_extracted.csv (each row is an (x, y) coordinate pair)

.. figure:: images/conan_surface_csv.png
    :alt: Contact Angle surface.csv screenshot
    :width: 159
    :align: center

    drop1/surface.csv (The coefficients of the surface line; first column is gradient, second column is y-intercept)

.. figure:: images/conan_tangents_csv.png
    :alt: Contact Angle tangents.csv screenshot
    :width: 159
    :align: center

    drop1/tangents.csv (The coefficients of the tangent lines at the contact point. First row is left tangent, second row is right tangent. First column is gradient, second column is y-intercept)Interfacial Tension
===================

A wizard-style window will guide you through the process of performing an interfacial tension analysis.

Image acquisition
-----------------

First, choose an image input method. OpenDrop currently supports opening images from the local filesystem or capturing images with a USB camera.

Local filesystem
^^^^^^^^^^^^^^^^

.. image:: images/ift_local_filesystem.png
    :alt: IFT Image acquisition, local filesystem
    :width: 914
    :align: center

Click on 'Choose files' to open the file chooser dialog and select an individual image or a sequence of images. When analysing a sequence of images, 'Frame interval' refers to the time interval (in seconds) between each image. Sequences of images are ordered in lexicographic order.

USB camera
^^^^^^^^^^

.. image:: images/ift_usb_camera.png
    :alt: IFT Image acquisition, USB camera
    :width: 914
    :align: center

Click on 'Connect camera' to open the camera chooser dialog.

.. image:: images/ift_usb_camera_camera_chooser.png
    :alt: IFT Image acquisition, USB camera chooser dialog
    :width: 828
    :align: center

OpenDrop uses OpenCV to capture images from a connected camera. 'Camera index' refers to the device index argument passed to the OpenCV function ``cv2.VideoCapture()``. An index of 0 refers to the first connected camera (usually a laptop's in-built webcam if present), an index of 1 refers to the second camera, and so on. Currently, there does not appear to be a way in OpenCV to query a list of valid device indices and associated device names, so in a multi-camera setup, some trial-and-error is required.

'Frame interval' refers to the time interval (in seconds) between capturing images.


Physical parameters
-------------------

.. image:: images/ift_physical_parameters.png
    :alt: IFT Physical parameters
    :width: 914
    :align: center


'Inner density' refers to the density of the drop.

'Outer density' refers to the density of the surrounding medium.

'Needle diameter' refers to the diameter of the needle the drop is suspended from.

'Gravity' refers to the gravitational acceleration.


Image processing
----------------

.. image:: images/ift_image_processing.png
    :alt: IFT Image processing
    :width: 914
    :align: center

The image processing window requires you to define the 'drop region' and 'needle region' of the image. Click on the 'Drop region' or 'Needle region' buttons in the 'Tools' panel, then drag over the image preview to define the associated region.

.. image:: images/ift_image_processing_regions_defined.png
    :alt: IFT Image processing, regions defined
    :width: 914
    :align: center

Once each region is defined, a blue outline will be drawn over the preview showing the drop or needle profile that has been extracted.

OpenDrop uses OpenCV's Canny edge detector to detect edges in the image, click on the 'Edge detection' button in the 'Tools' panel to open a dialog bubble which will allow you to adjust the lower and upper threshold parameters of the Canny edge detector. Thin blue lines are drawn over the preview to show detected edges.

The extracted needle profile is used to determine the diameter in pixels of the needle in the image. Along with the needle diameter in millimetres given in the 'Physical parameters' page, a metres-per-pixel scale can be determined, which is then used to derive other physical properties of the drop after the image is analysed.

Click on 'Start analysis' to begin analysing the input images, or begin capturing and analysing images if using a camera.


Results
-------

.. image:: images/ift_results.png
    :alt: IFT Results
    :width: 914
    :align: center

The results page shows the current status of the analysis. Data shown in the window is updated as the analysis progresses.

There are two main views, the 'Individual Fit' view and the 'Graphs' view. The 'Graphs' view is not available when analysing a single image.

Individual Fit
^^^^^^^^^^^^^^

The 'Individual Fit' view shows analysis details for an individual image. Pick an analysis in the lower panel to preview its details in the upper panel.

The 'Drop profile' tab on the right of the upper panel shows the fitted drop profile (drawn in magenta) over the extracted drop profile (drawn in blue).

.. image:: images/ift_results_drop_profile.png
    :alt: IFT Results, drop profile
    :width: 591
    :align: center

The 'Fit residuals' tab shows a plot of the fit residuals. The horizontal axis is the 'drop profile parameter', ranging from 0 to 1, with 0 corresponding to one end of the drop edge outline, and 1 corresponding to the other end. The vertical axis is some dimensionless quantity indicating the deviation of the extracted profile from the fitted profile.

.. image:: images/ift_results_fit_residuals.png
    :alt: IFT Results, fit residuals
    :width: 591
    :align: center

The 'Log' tab shows the history of any messages logged by the fitting routine.

.. image:: images/ift_results_log.png
    :alt: IFT Results, log
    :width: 591
    :align: center

Graphs
^^^^^^

.. image:: images/ift_results_graphs.png
    :alt: IFT Results, graphs
    :width: 914
    :align: center

The 'Graphs' view shows plots of interfacial tension, volume, and surface area over time.

Cancel or discard analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^

You may cancel an in progress analysis by clicking on the 'Cancel' button in the footer (not shown in the screenshots above). To discard the results of a finished analysis, click the 'Back' button, which will return you to the 'Image processing' page, or close the window to return to the Main Menu.


Saving
------

.. image:: images/ift_save_dialog.png
    :alt: IFT Save dialog
    :width: 828
    :align: center

Once an analysis is finished, click on the 'Save' button in the footer to open the save dialog. All data will be saved in a folder with name determined by the 'Name' entry, and in a parent directory determined by the 'Parent' selection. 

As a convenience, you may choose to save some pre-made plots.

.. image:: images/ift_save_output.png
    :alt: IFT Example save output
    :width: 662
    :align: center

An example save output is shown above, and screenshots of the contents of some files are shown below.

.. figure:: images/ift_timeline_csv.png
    :alt: IFT timeline.csv screenshot
    :width: 1000
    :align: center

    timeline.csv

.. figure:: images/ift_profile_fit_csv.png
    :alt: IFT profile_fit.csv screenshot
    :width: 159
    :align: center

    water_in_air1/profile_fit.csv (each row is an (x, y) coordinate pair)

.. figure:: images/ift_profile_extracted_csv.png
    :alt: IFT profile_extracted.csv screenshot
    :width: 159
    :align: center

    water_in_air1/profile_extracted.csv (each row is an (x, y) coordinate pair)

.. figure:: images/ift_profile_fit_residuals_csv.png
    :alt: IFT profile_fit_residuals.csv screenshot
    :width: 159
    :align: center

    water_in_air1/profile_fit_residuals.csv (first column is 'drop profile parameter', second column is residual)

.. figure:: images/ift_params_ini.png
    :alt: IFT params.ini screenshot
    :width: 485
    :align: center

    water_in_air1/params.ini
(todo: This page is out of date and should be updated.)

*****
Usage
*****

When OpenDrop is launched, the Main Menu window will first appear.

.. image:: images/main_menu.png
    :alt: Main Menu window
    :width: 465
    :align: center

Click on either of the 'Interfacial Tension' or 'Contact Angle' buttons to begin a new analysis of the respective type.

.. include:: ift.rst
.. include:: conan.rst
.. include:: genicam.rst
.. include:: notes.rst
