
---
title: 'pyFBS: A Python package for Frequency Based Substructuring'
tags:
  - Python
  - Structural dynamics
  - Frequency Based Substructuring
  - System Equivalent Model Mixing
  - Transfer Path Analysis
authors:
  - name: Tomaž Bregar
    affiliation: 1 
  - name: Ahmed El Mahmoudi
    affiliation: 2
  - name: Miha Kodrič
    affiliation: 3
  - name: Domen Ocepek
    affiliation: 3
  - name: Francesco Trainotti
    affiliation: 2
  - name: Miha Pogačar
    affiliation: 3
  - name: Mert Göldeli
    affiliation: 2	
  - name: Gregor Čepon
    affiliation: 3
  - name: Miha Boltežar
    affiliation: 3
  - name: Daniel J. Rixen
    affiliation: 2
affiliations:
 - name: Gorenje d.o.o., Partizanska 12, 3503 Velenje, Slovenia
   index: 1
 - name: Technical University of Munich, Institute of Applied Mechanics, Boltzmannstr.  15, 85748 Garching, Germany
   index: 2
 - name: Faculty of Mechanical Engineering, University of Ljubljana, Aškerčeva 6, 1000 Ljubljana, Slovenia
   index: 3
date: 13 August 2017
bibliography: paper.bib
---

# Summary

In science, engineering and technology complex problems are often decomposed into smaller, simpler subsystems. 
Each subsystem can then be analyzed and evaluated separately. 
This approach can often reduce the complexity of the overall problem and provide invaluable insight into the optimization and troubleshooting of each individual component. 
The subsystems can also be assembled back together and with that the system can be analyzed as a whole.

Dynamic Substructuring (DS) is an engineering concept where dynamic systems are modeled and analyzed in terms of their components or so-called substructures. 
There are several ways of formulating the dynamics of substructures. One of them is with Frequency Response Functions (FRFs), which describe the response as the result of a unit harmonic force. 
The method is well suited for experimental approaches where FRFs are obtained from measurement of components. Such approaches were already investigated in the 1970s [@KLOSTERMAN_1971_PHD]  and 1980s [@MARTINEZ_1984_COMBINEDEXPANALYTICAL; @KLOSTERMAN_1984_SMURF; @JETMUNDSEN_1988_FBS; @URGUEIRA_1989_DYNAMIC]. 
Due to complicated formulations and difficulties in obtaining good measurements, the method was hardly applicable. 
Thanks to better measurement hardware and proper formulation of the problem, Frequency Based Substructuring (FBS) has gained popularity in recent years [@deKlerk2008; @vanderSeijs2016; @RIXEN_2006_GUITAR]. 
With this approach, it is also possible to build hybrid models in which experimentally characterized and numerically modelled parts are combined.

pyFBS is a Python package for Frequency Based Substructuring. The package implements an object-oriented approach for dynamic substructuring. 
Current state-of-the-art methodologies in frequency based substructuring are available in pyFBS. Each method can be used as a standalone or interchangeably with others. 
Also a 3D display is available so a user can simply place and orient associated input/outputs used with each method.
Furthermore, basic and application examples are provided with the package together with real experimental and numerical data [@ahmed]. 


# Statement of need

Evaluating structural dynamics is a necessary step in the development of any complex mechanical system. 
Vibro-acoustic character together with the visual design contributes to the customer's perception of a premium product.
With DS the vibro-acoustic performance of a product can be analysed in terms of its subcomponents. This approach is highly beneficial, as almost all complex products are designed modularly.
pyFBS helps the user to perform experimental modelling in DS. It enables an intuitive way to position sensors and impact locations on the analysed structures.
The position and orientation of sensors/impacts can be then be used in each implemented DS method. 
Furthermore, experimental measurements can be simulated with ease from a numerical model, where the same positional information can be used.   

To the best of authors knowledge there is currently no open source software alternative, which would enable the user to use dynamic substructuring methodologies. 
pyFBS has been designed to be used for scientific research in the field of dynamic substructuring. 
It is currently being used by a number of undergraduate students and postgraduate researchers. 


# Features

pyFBS enables the user to use state-of-the-art dynamic substructuring methodologies in an intuitive manner.
Currently implemented features are listed below. 

## 3D display

Structures and positions of impacts, sensors and channels can be visualized in 3D in \autoref{fig:3D}. 
The 3D display is built on top of PyVista [@sullivan2019pyvista] and enables an intuitive way to display relevant data. 
Sensors and impacts can be interactively positioned on structures and the updated positions can be directly used within the pyFBS.
With this feature the experimental setup can be prepared in advance, to avoid possible mistakes in experimental modelling.
Furthermore, various animations can be performed directly in the 3D display, such as the animation of mode shapes or operational deflection shapes.

![An example of a simple structure depicted in the pyFBS 3D display.\label{fig:3D}](./images/figure.png)

## FRF synthetization

Frequency Response Functions can be synthetized for predefined positions of channels and impacts in a numerical model. 
Currently, mode superposition FRF synthetization is supported, where mass and stiffness matrices are imported from FEM software. 
Damping can be introduced as modal damping for each mode shape. 
Additionally, noise can be introduced to the response so a realistic set of FRFs, representing experimental measurements, can be obtained.

## Virtual Point Transformation (VPT) 

VPT projects measured dynamics on the predefined interface displacement modes (IDMs) [@solvingRDOF]. 
The interface is usually considered to be rigid; therefore, only 6 rigid IDMs are used in the transformation. 
After applying the transformation, a collocated set of FRFs is obtained, which can afterwards directly be used in DS. 
Expanded VPT is also supported, where directly measured rotational response is included in the transformation [@Bregar2020].

## System Equivalent Model Mixing (SEMM)

SEMM enables mixing of two equivalent frequency-based models into a hybrid model [@Klaassen2018; @semm_svd]. 
The models used can either be of numerical or experimental nature. One of the models provides the dynamic properties (overlay model) and the second model provides a set of degrees of freedom. 
A numerical model is commonly used as a parent model and an experimental model is used as an overlay model. 

# Discussion

The development of the pyFBS is an ongoing effort. Currently the package is used as a research tool. 
In the future, more examples on the topic of DS and applications of Transfer Path Analysis (TPA) are going to be introduced in the documentation. 
Furthermore, implementation of Operational Source Identification (OSI) is going to be integrated into the pyFBS in the near future.

# Acknowledgments

The pyFBS package was developed as a part of collaboration between the Laboratory for Dynamics of Machines and Structures (LADISK), Faculty of Mechanical Engineering, University of Ljubljana (UL FME) 
and the Chair of Applied Mechanics (AM), Technical University of Munich (TUM).

# References
pyFBS is an ongoing open-source project and everybody is welcome to contribute to:
 * asking questions (https://gitlab.com/pyFBS/pyFBS_support/-/issues),
 * reporting bugs (https://gitlab.com/pyFBS/pyFBS/-/issues),
 * feature requests (https://gitlab.com/pyFBS/pyFBS/-/issues),
 * adding new code (by creating Merge Request).

Asking questions
----------------
Asking questions and making them public is immensely helpful for anyone 
who will ever encounter a similar problem or dilemma. 
The community of developers and other ``pyFBS`` users will try to answer your question 
which will benefit both developers and users.

To ask a question, please create an issue on our support page https://gitlab.com/pyFBS/pyFBS_support/-/issues. 
You can also write to us via email at info.pyfbs@gmail.com.

Reporting bugs
--------------
``pyFBS`` is a still-evolving library and that's why you might encounter some unexpected and unintuitive code behavior. 
If you encounter such a bug, please report it by opening an issue at https://gitlab.com/pyFBS/pyFBS/-/issues and mark it with the Bug label.  
To make it easier for the developers to reproduce the bug in order to fix it, please submit a minimal working example
(e.g. support the problem with a screenshot or sample files).

Feature requests
----------------
We will welcome suggestions for improving existing and introducing new functionalities to the ``pyFBS``. 
Please post them by opening the issue at https://gitlab.com/pyFBS/pyFBS/-/issues and mark it with the New feature label. 
Add a brief description of the proposed feature and outline its benefits. 
In addition, you can equip your suggestions with pictures or links to relevant references.

Adding new code
---------------
Contribution to ``pyFBS`` can also be made by adding new code or writing documentation on the existing one.
Before adding a new code, please open the issue with the appropriate label (New feature, Documentation). 
This allows us to determine through the discussion whether the proposed subject fits ``pyFBS``, 
if the proposed topic is not already under development, 
and also to assign the users who will work on the topic.

Cloning
^^^^^^^
Before starting writing your own code you have to download the latest version of the ``pyFBS`` library by running:

.. code-block:: 

    git clone https://gitlab.com/pyFBS/pyFBS.git
    cd pyFBS
    python -m pip install -e .

You can also fork the repository and clone it from your account.

Creating new branch
^^^^^^^^^^^^^^^^^^^
New code is always added via a new branch. 
Please use an informative and descriptive branch name. 
As a branch name, you can also use the number of the issue you are trying to resolve (like iss5).

.. code-block:: 

    git branch iss5

Coding
^^^^^^
Once you create your own branch, you can start making changes to the repository. 
When adding new functionalities, try to follow the current code structure. 
You should always add new code to the most content-related file. 
If the new functionality does not relate to any of the existing files, start a discussion in the open issue 
on how it would make sense to implement this functionality.

Code style
**********
Code should follow the philosophy of the `Python programming language <https://en.wikipedia.org/wiki/Python_(programming_language)#Design_philosophy_and_features>`_:
 - Beautiful is better than ugly.
 - Explicit is better than implicit.
 - Simple is better than complex.
 - Complex is better than complicated.
 - Readability counts.

The naming and code layout convention should follow the PEP 8. 
The exception is line widths which are permitted to exceed 79 characters.

Adding additional comments between lines of code is most welcome, as it greatly increases the intelligibility of the code.

Documentation
^^^^^^^^^^^^^

Consistent documentation is a crucial feature of any great library. 
A brief description of the function's content, input and output parameters, and an example of its use 
allows other users to implement this function correctly.

Documentation style
*******************

Every function must have a docstring as demonstarted in the following example.  

.. code-block:: python

    def my_function(my_param_1, my_param_2):
        """
        Returns sum of my_param_1 and my_param_2.
        
        :parm my_param_1: the first parameter of summation
        :type my_param_1: array
        :parm my_param_2: the second parameter of summation
        :type my_param_2: array
        
        :return: sum of ``my_param_1`` + ``my_param_2``
        :rtype: array
        
        Example: 
        >>>a = my_function_1(1, 1)
        >>>print(a)
        2
        """
        result = my_param_1 + my_param_2
        return result

Docstrings are defined inside ``""" """``. 
First, provide a brief introduction of the function. 
Then, the input parameters are listed using ``:parm parameter_name:`` command, followed by a short description of the parameters.
It is advisable to define the type of the input parameter using ``:parm parameter_name:`` command, followed by the variable type.
At the end of the list of all the input parameters, the output of the function is defined using the ``:return:`` command.
Type of the output is defined by ``:rtype:``
The docstring can end with an example demonstrating the implementation of the function and expected result.

Notebook examples
*****************
``pyFBS`` library has a collection of notebooks in which most of the functionalities are depicted using simple examples. 
If you want to add the example of new functionalities, please find the appropriate existing notebook or create a new 
one inside the ``.\examples\`` folder. 
These notebooks are also a part of the testing procedure, so make sure that all notebooks run without errors.

Online documentation
********************
Documentation, displayed at https://pyfbs.readthedocs.io/en/latest/, is located at ``.\docs\``.
These pages include the theoretical background of methods included in ``pyFBS``. 
All files are in form of `Restructured Text (reST) <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_.

At the beginning of each topic there is an introduction of the method followed by the fundamental equations and relevant references. 
In the end, an extensive description of how to use the introduced topic using the ``pyFBS`` functionalities is provided. 

Testing
^^^^^^^
After making the changes, please test changes locally first before creating a merge request.
The code testing is fully automated. To test the code, you have to install the ``tox`` library:

.. code-block:: 

    pip install tox

In addition, it is necessary to install requirements for developers, listed in ``requirements_dev.txt``, using the command:

.. code-block:: 

    pip install -r requirements_dev.txt 

Testing code
************
Once the ``tox`` is installed, you just have to run the ``tox.ini`` script using the command:

.. code-block:: 

    cd pyFBS
    tox

The ``tox`` script will create a virtual environment and will test all notebook examples and all tests that are defined inside ``.\test\`` folder.

Testing documentation
*********************

Before building the documentation, execute the following command:

.. code-block:: 

    cd doc
    pip install -r requirements_dev.txt 

Documentation is tested separately by running commands:

.. code-block:: 

    make clean
    make html

Generated documentation can be found in ``.\_build\html``. 
Here you can open HTML pages in your browser and see your changes. 

Creating Merge Request
^^^^^^^^^^^^^^^^^^^^^^
When the changes pass all local tests, it's' time to create a merge request.
When creating a merge request, add a short description and assign code reviewers that will check the changes and accept the merge.
Creating a merge request will automatically run the continuous integration (CI) testing. 
If a merge request resolves one or more issues, mention this in the description of the merge request using ``Closes #4, #6``.
This will automatically close mentioned issues once the branch will be merged. 
More useful commands can be found here: https://docs.gitlab.com/ee/user/project/issues/managing_issues.html
Development Lead
----------------

* Tomaž Bregar 

* Ahmed El Mahmoudi

* Miha Kodrič

* Domen Ocepek

* Francesco Trainotti

* Miha Pogačar

* Mert Göldeli

Contributors
------------

* Gregor Čepon

* Miha Boltežar

* Daniel J. Rixen

For a full list of `contributors`_ check the repository.

.. _contributors: https://gitlab.com/pyFBS/pyFBS/-/graphs/master
.. image:: https://gitlab.com/pyFBS/pyFBS/-/raw/master/docs/logo/logo-big.png
	:align: right
	:width: 300

pyFBS
-----
	
pyFBS is a Python package for Frequency Based Substructuring and Transfer Path Analysis. It enables the user to use state-of-the-art dynamic substructuring methodologies in an intuitive manner. 

With the package also basic and application examples are provided, together with real datasets so you can directly try out the capabilities of the pyFBS.

|pypi| |docs| |codequality| |MIT|

Features
--------

* 3D display

* FRF synthetization

* Virtual Point Transformation

* System Equivalent Model Mixing

* Singular Vector Transformation

For more information on features, basic and application examples check out the `documentation`_. 

Citation
--------
A paper about the pyFBS will be submitted to the Journal of Open Source Software journal. If you will be using pyFBS in your scientific research, please consider citing the paper.

License
-------
Licensed under the MIT license.

.. _documentation: https://pyfbs.readthedocs.io/en/latest/intro.html

.. |pypi| image:: https://img.shields.io/pypi/v/pyfbs?style=flat-square
   :target: https://pypi.org/project/pyfbs/

.. |docs| image:: https://readthedocs.org/projects/pyfbs/badge/?version=latest
   :target: https://pyfbs.readthedocs.io/en/latest/?badge=latest

.. |MIT| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
   
.. |codecov| image:: https://codecov.io/gl/pyFBS/pyFBS/branch/\x6d6173746572/graph/badge.svg?token=XSGM89JGMF
   :target: https://codecov.io/gl/pyFBS/pyFBS

.. |codequality| image:: https://app.codacy.com/project/badge/Grade/dbb59e10c07543b6b61c083a09eac500    
   :target: https://www.codacy.com/gl/pyFBS/pyFBS/dashboard?utm_source=gitlab.com&amp;utm_medium=referral&amp;utm_content=pyFBS/pyFBS&amp;utm_campaign=Badge_Grade=====
Usage
=====

To use :mod:`pyFBS` within a project simply import the package:

.. code-block:: python

	import pyFBS

************
Example data
************
	
To test out the capabilities of :mod:`pyFBS` also two example datasets are available directly with the package. 
For each testbench structure a dictionary is available containing relative path to the predefined datasets.

Predefined datasets are also used directly in basic and application examples.

Academic testbench
==================
The first testbench is an academic example :func:`pyFBS.download_lab_testbench`. 
The testbench is used to evaluate and compare different dynamic substructuring methodologies. 

.. figure:: ./examples/data/structure_t.png
   :width: 800px
   
   An example of a academic substructuring testbench depicted in the pyFBS 3D display.


Example datasets for academic testbench contain:

* STL files of the testbench (e.g. ``./lab_testbench/STL/A.stl``),

* FEM of each substructure (e.g. ``./lab_testbench/FEM/A.full``),

* Excel files of positional data for sensors and impacts (e.g. ``./lab_testbench/Measurements/coupling_example.xlsx``),

* Experimental FRF measurements (e.g. ``./lab_testbench/Measurements/Y_A.p``).

Automotive testbench
====================

The second testbench is an automotive example, which can be downloaded with :func:`pyFBS.download_automotive_testbench`. 
The automotive testbench was designed to represent an engine-transmission unit’s suspension from a real car.


.. figure:: ./examples/data/structure.png
   :width: 800px
   
   An example of a automotive testbench depicted in the pyFBS 3D display.
   
Example datasets for automotive testbench contain:

* STL files of the testbench (e.g. ``./automotive_testbench/STL/receiver.stl``),

* Excel files of positional data for sensors and impacts (e.g. ``./automotive_testbench/Measurements/A.xlsx``),

* Experimental FRF measurements (e.g. ``./automotive_testbench/Measurements/A.p``).

********
Features
********

3D display
==========
With the pyFBS substructures and positions of impacts, sensors and channels can be visualized in 3D display :mod:`pyFBS.view3D`. 
The 3D display uses PyVista [1]_ for the visualization and enables an intuitive way to display relevant data. 
Sensors and impacts can be interactively positioned on the substructures and the updated positions can be directly used within pyFBS. 
Furthermore, various animations can be performed directly in the 3D display, such as the animation of mode shapes or operational deflection shapes.

One of the main features of the pyFBS is also the ability to synthetize FRFs directly from the predefined positions of channels and impacts. 
Currently, mode superposition FRF synthetization is supported, where mass and stiffness matrices are imported from FEM software. 
Damping can be introduced as modal damping for each mode shape. Additionally, noise can be introduced to the response so a realistic set of FRFs, representing experimental measurements, can be obtained.


FRF synthetization
==================
One of the main features of the pyFBS is also the ability to synthetize FRFs directly from the predefined positions of channels and impacts :mod:`pyFBS.MK_model`. 
Currently, mode superposition FRF synthetization is supported, where mass and stiffness matrices are imported from FEM software. 
Damping can be introduced as modal damping for each mode shape. 
Additionally, noise can be introduced to the response so a realistic set of FRFs, representing experimental measurements, can be obtained.


Virtual Point Transformation
============================
Within the pyFBS Virtual Point Transformation (VPT) :mod:`pyFBS.VPT` is implemented [2]_. 
VPT projects measured dynamics on the predefined interface displacement modes (IDMs). 
The interface is usually considered to be rigid; therefore, only 6 rigid IDMs are used in the transformation. 
After applying the transformation, a collocated set of FRFs is obtained, which can afterwards directly be used in DS. 
Expanded VPT is also supported, where directly measured rotational response is included in the transformation [3]_.


System Equivalent Model Mixing
==============================
The pyFBS supports System Equivalent Model Mixing (SEMM) [4]_. 
SEMM enables mixing of two equivalent frequency-based models into a hybrid model. 
The models used can either be of numerical or experimental nature. 
One of the models provides the dynamic properties (overlay model) and the second model provides a set of degrees of freedom. 
A numerical model is commonly used as a parent model and an experimental model is used as an overlay model. 


Singular Vector Transformation
==============================
Singular Vector Transformation (SVT) consists in projecting the acquired data into subspaces composed by dominant singular vectors, 
which are extracted directly from the available FRF datasets [5]_. 
No geometrical and/or analytical model is required. 
If some basic requirements are met, the reduced orthonormal frequency dependent basis would be able to control and 
observe most of the rigid and flexible vibration modes of interest over a broad frequency range. 
The SVT can tackle challenging scenarios with flexible behaving interfaces and lightly damped systems.

.. rubric:: References

.. [1] Sullivan C, Kaszynski A. PyVista: 3D plotting and mesh analysis through a streamlined interface for the Visualization Toolkit (VTK). Journal of Open Source Software. 2019 May 19;4(37):1450.
.. [2] de Klerk D, Rixen DJ, Voormeeren SN, Pasteuning P. Solving the RDoF problem in experimental dynamic sybstructuring. InInternational modal analysis conference IMAC-XXVI 2008 (pp. 1-9).
.. [3] Tomaž Bregar, Nikola Holeček, Gregor Čepon, Daniel J. Rixen, and Miha Boltežar. Including directly measured rotations in the virtual point transformation. Mechanical Systems and Signal Processing, 141:106440, July 2020.
.. [4] Steven WB Klaassen, Maarten V. van der Seijs, and Dennis de Klerk. System equivalent model mixing. Mechanical Systems and Signal Processing, 105:90–112, 2018.
.. [5] Trainotti F, Bregar T, Klaassen SW, Rixen DJ. Experimental decoupling of substructures by singular vector transformation. Mechanical Systems and Signal Processing. 2022 Jan 15;163:108092.


=======
Credits
=======

.. include:: ../AUTHORS.rst==========
Contribute
==========

.. include:: ../CONTRIBUTING.rst============
Installation
============
:mod:`pyFBS` is supported on Python versions 3.5+. You can install the pyFBS with following the instructions.


****
PyPI
****

:mod:`pyFBS` can be installed from PyPI using ``pip``:

.. code-block:: python

	pip install pyFBS

******
Source
******
The source code of the pyFBS can be also downloaded (cloned) by running:

.. code-block:: python

	git clone https://gitlab.com/pyFBS/pyFBS.git
	cd pyFBS
	pip install -e .Welcome to pyFBS's documentation!
=================================

.. include:: ../README.rst

Acknowledgements
----------------


The :mod:`pyFBS` was developed as a part of collaboration between the Laboratory for Dynamics of Machines and Structures (`LADISK <http://ladisk.si/>`_), University of Ljubljana, Faculty of Mechanical Engineering (UL FME) and the Chair of Applied Mechanics (`AM <https://www.mw.tum.de/am/home/>`_), Technical University of Munich (TUM).

.. raw:: html

	<a href="http://ladisk.si/">
		<img src="./_static/ladisk_logo.svg" alt="LADISK logo" title="LADISK" align="right" width=200 />
	</a>

	
Laboratory for Dynamics of Machines and Structures (LADISK), UL FME
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The LADISK research group of the Faculty of Mechanical Engineering at the University of Ljubljana is the leading research laboratory in the field of vibroacoustic and structural dynamics in Slovenia. 
The research activities include the simulation and measurement of structural vibrations as well as the development of methods for modeling rigid-flexible multibody systems with one-sided contacts. 
In the field of structural dynamics, research activities focus on theoretical and experimental modal analysis as well as the management of vibration fatigue. 
Different optical and acoustic methods are researched to obtain full-field dynamic parameters of the structure.
In recent years the group has been actively involved in the research of substructuring approaches in the frequency domain and their application for efficient noise, vibration and harshness engineering. 
	
	
.. raw:: html

	<a href="https://www.mw.tum.de/">
		<img src="./_static/tum_logo.svg" alt="TUM logo" title="TUM" align="right" width=200 />
	</a>


Chair of Applied Mechanics (AM), TUM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Chair of Applied Mechanics (Technical University of Munich) covers a wide range of research fields, that directly impact the future of industry and society, but are always related to Structural Dynamics and Mechatronics.
These research fields are divided into three areas: Experimental Dynamics, Numerical Methods, and Robotics. They range from fundamental research, such as development of new numerical methods, to applied research in cooperation with industry. 
This research can help on the path to a sustainable future by providing tools and technologies for a more advanced product design, that reduces the amount of prototypes and, thus, saves development costs and minimizes material consumption.

..
   .. toctree::
      :maxdepth: 1
      :caption: Why pyFBS?

      intro

.. toctree::
   :maxdepth: 2
   :caption: About

   intro
   about


.. toctree::
   :maxdepth: 1
   :caption: Examples
   
   ./examples/basic_examples
   ./examples/frequency_based_substructuring
   ./examples/transfer_path_analysis
   ./examples/case_studies
   


.. toctree::
   :maxdepth: 2
   :caption: Code Documentation

   ./code_documentation/core
***************
About
***************

.. toctree::
    :hidden:
 
    ./credits.rst 
    ./license.rst 
    ./contribute.rst
    ./installation.rst 
    ./usage.rst 

.. panels::
    :column: col-lg-6 col-md-6 col-sm-12 col-xs-12 p-3

    Credits
    ^^^^^^^^^^^^

    See the authors of the pyFBS project.

    .. link-button:: credits
        :type: ref
        :text: Credits
        :classes: btn-outline-primary btn-block stretched-link
    
    ---
    License
    ^^^^^^^^^^^^

    The project is under the MIT license.

    .. link-button:: license
        :type: ref
        :text: License
        :classes: btn-outline-primary btn-block stretched-link
    
    ---
    Contribute
    ^^^^^^^^^^^^

    Everyone can contribute to the pyFBS library, take a look at you can you do it.

    .. link-button:: contribute
        :type: ref
        :text: Contribute
        :classes: btn-outline-primary btn-block stretched-link
 
    ---
    Installation
    ^^^^^^^^^^^^

    Installation process for Python version 3.8.

    .. link-button:: installation
        :type: ref
        :text: Installation
        :classes: btn-outline-primary btn-block stretched-link
    ---
    Usage
    ^^^^^^^^^^^^

    Get familiar with the example data and features that pyFBS offers.

    .. link-button:: usage
        :type: ref
        :text: Usage
        :classes: btn-outline-primary btn-block stretched-link
==========
Why pyFBS?
==========
:mod:`pyFBS` is a Python package for Frequency Based Substructuring and Transfer Path Analysis. 
The package implements an object-oriented approach for dynamic substructuring. 
Current state-of-the-art methodologies in frequency based substructuring are available in pyFBS. 
Each method can be used as a standalone or interchangeably with others. 
Furthermore, basic and application examples are provided with the package together with real experimental and numerical data. 
The pyFBS has been designed to be used for scientific research in the field of dynamic substructuring. 
It is currently being used by a number of undergraduate students and postgraduate researchers. 

.. figure:: ./logo/pyFBS_logo_presenttion.gif
   :width: 800px

**********************
Dynamic Substructuring
**********************
In science, engineering and technology complex problems are often decomposed into smaller, simpler subsystems. 
Each subsystem can then be analyzed and evaluated separately. 
This approach can often reduce the complexity of the overall problem and provide invaluable insight into the optimization and troubleshooting of each individual component. 
The subsystems can also be assembled back together and with that the system can be analyzed as a whole.

Dynamic Substructuring (DS) is an engineering concept where dynamic systems are modeled and analyzed in terms of their components or so-called substructures. 
There are several ways of formulating the dynamics of substructures. One of them is with Frequency Response Functions (FRFs), which describe the response as the result of a unit harmonic force. 
The method is well suited for experimental approaches where FRFs are obtained from measurement of components. Such approaches were already investigated in the 70s [1]_  
and 80s [2]_ [3]_ [4]_  [5]_. 
Due to complicated formulations and difficulties in obtaining good measurements, the method was hardly applicable. 
Thanks to better measurement hardware and proper formulation of the problem,  Frequency Based Substructuring (FBS) has gained popularity in recent years [6]_ [7]_ [8]_.  
With this approach, it is also possible to build hybrid models in which experimentally characterized and numerically modelled parts are combined.

.. rubric:: References

.. [1] Albert L. Klosterman. A combined experimental and analytical procedure for improving automotive system dynamics. PhD thesis, University of Cincinnati, Department of Mechanical Engineering, 1971.
.. [2] David R. Martinez, Thomas G. Carrie, Dan L. Gregory, and A. Keith Miller. Combined Experimental/Analytical Modelling using component modes synthesis. In 25th Structures, Structural Dynamics and Materials Conference, 140–152. Palm Springs, CA, USA, 1984.
.. [3] John R. Crowley, Albert L. Klosterman, G. Thomas Rocklin, and Havard Vold. Direct structural modification using frequency response functions. Proceedings of IMAC II, February 1984.
.. [4] Bjorn Jetmundsen, Richard L. Bielawa, and William G. Flannelly. Generalized frequency domain substructure synthesis. Jnl. American Helicopter Society, 33(1):55–64, 1988.
.. [5] Antonio Paulo Vale Urgueira. Dynamic analysis of coupled structures using experimental data. PhD thesis, Imperial College, London, 1989.
.. [6] D.J. Rixen, T. Godeby, and E. Pagnacco. Dual assembly of substructures and the fbs method: application to the dynamic testing of a guitar. International Conference on Noise and Vibration Engineering, ISMA, September 18-20 2006.
.. [7] de Klerk D, Rixen DJ, Voormeeren SN. General framework for dynamic substructuring: history, review and classification of techniques. AIAA journal. 2008 May;46(5):1169-81.
.. [8] Maarten V. van der Seijs, Dennis de Klerk, and Daniel J. Rixen. General framework for transfer path analysis: history, theory and classification of techniques. Mechanical Systems and Signal Processing, 68-69:217–244, February 2016.==========
pyFBS.SEMM
==========

.. autofunction:: pyFBS.SEMM
.. autofunction:: pyFBS.identification_algorithm
==============
pyFBS.MK_model
==============

.. autoclass:: pyFBS.MK_model
    :members:
============
pyFBS.view3D
============

.. autoclass:: pyFBS.view3D
    :members:

========
pyFBS.IO
========

.. autofunction:: pyFBS.load_uff_file_PAK
=========
pyFBS.SVT
=========

.. autoclass:: pyFBS.SVT
    :members:

=============
pyFBS.utility
=============

.. autofunction:: pyFBS.modeshape_sync_lstsq
.. autofunction:: pyFBS.modeshape_scaling_DP
.. autofunction:: pyFBS.MCF
.. autofunction:: pyFBS.flatten_FRFs
.. autofunction:: pyFBS.unflatten_modes
.. autofunction:: pyFBS.complex_plot
.. autofunction:: pyFBS.complex_plot_3D
.. autofunction:: pyFBS.mode_animation
.. autofunction:: pyFBS.coh_frf
.. autofunction:: pyFBS.dict_animation
.. autofunction:: pyFBS.CMIF
.. autofunction:: pyFBS.TSVD
.. autofunction:: pyFBS.M
.. autofunction:: pyFBS.angle
.. autofunction:: pyFBS.rotation_matrix_from_vectors
.. autofunction:: pyFBS.unit_vector
.. autofunction:: pyFBS.angle_between
.. autofunction:: pyFBS.generate_channels_from_sensors
.. autofunction:: pyFBS.generate_sensors_from_channels
.. autofunction:: pyFBS.coh_on_FRF
.. autofunction:: pyFBS.orient_in_global
.. autofunction:: pyFBS.orient_in_global_2

=========
pyFBS.VPT
=========

.. autoclass:: pyFBS.VPT
    :members:

========
Core API
========

.. toctree::
	:maxdepth: 2

	view3D
	MCK
	VPT
	SEMM
	SVT

=========
Utilities
=========

.. toctree::
	:maxdepth: 2

	IO
	utility======================
Transfer Path Analysis
======================
Transfer-path analysis (TPA) is a reliable and effective diagnostic tool for the characterization of actively vibrating components and the propagation of noise and vibrations to the connected passive substructures [1]_. 
TPA offers the ability to analyse the vibration transfer between the individual components of the assembly, distinguish the partial transfer-path contribution and predict the receiver's response.

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Clasiccal TPA">

.. only:: html

    .. figure:: ./data/classical_tpa.svg   
       :target: ./tpa/classical_tpa.html

       Clasiccal TPA

.. raw:: html

    </div>

.. toctree::
   :hidden:

   ./tpa/classical_tpa





.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Component-based TPA">

.. only:: html

    .. figure:: ./data/component-based_tpa.svg   
       :target: ./tpa/component-based_tpa.html

       Component-based TPA

.. raw:: html

    </div>

.. toctree::
   :hidden:

   ./tpa/component-based_tpa





.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Transmssibility-based TPA">

.. only:: html

    .. figure:: ./data/otpa.svg   
       :target: ./tpa/transmissibility-based_tpa.html

       Transmissibility-based TPA

.. raw:: html

    </div>

.. toctree::
   :hidden:

   ./tpa/transmissibility-based_tpa


.. rubric:: References

.. [1] Maarten V. van der Seijs, Dennis de Klerk, and Daniel J. Rixen. General framework for transfer path analysis: history, theory and classification of techniques. Mechanical Systems and Signal Processing, 68-69:217–244, February 2016.==============================
Case studies
==============================
These examples show applications of the pyFBS on more complex problems. Explore this application examples to see how pyFBS can be used on more complex dynamic problems.

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Operational Deflection Shapes">

.. only:: html

    .. figure:: ./data/ods.gif	   
       :target: ./case_studies/06_ODS.html

       Operational Deflection Shapes

.. raw:: html

    </div>

.. toctree::
   :hidden:

   ./case_studies/06_ODS


  
 
.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Experimental Modal Analysis">

.. only:: html

    .. figure:: ./data/modal_2_min.gif
       :target: ./case_studies/09_EMA.html

       Experimental Modal Analysis

.. raw:: html

    </div>

.. toctree::
   :hidden:

   ./case_studies/09_EMA
   

  
 
.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Transmission Simulator">

.. only:: html

    .. figure:: ./data/ten_display_four.png   
       :target: ./case_studies/10_TS.html

       Transmission Simulator

.. raw:: html

    </div>

.. toctree::
   :hidden:

   ./case_studies/10_TS==============================
Frequency Based Substructuring
==============================
The methodology to divide large and complex systems into several subsystems is a common practice in the field of structural dynamics. 
Structural dynamic analyses can be carried out more efficiently if complex systems are divided into smaller subsystems, analysed separately, and later coupled using dynamic substructuring (DS) methods.
In terms of the modeling domain, a frequency-based substructuring (FBS) is often preferred by experimentalists due to its ease of use and implementation with directly measured Frequency Response Functions (FRFs). 
In this context, datasets of measured transfer functions constitute the dynamic models of the substructures involved in the assembly/disassembly process.
These examples show state-of-the-art techniques to successfully couple or decouple substructures using FBS framework.



.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Virtual Point Transformation">

.. only:: html

    .. figure:: ./data/vpt_scheme.svg   
       :target: ./fbs/virtual_point_transformation.html

       Virtual Point Transformation

.. raw:: html

    </div>

.. toctree::
   :hidden:

   ./fbs/virtual_point_transformation

   


.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Singular Vector Transformation">

.. only:: html

    .. figure:: ./data/SVT.png   
       :target: ./fbs/singular_vector_transformation.html

       Singular Vector Transformation

.. raw:: html

    </div>

.. toctree::
   :hidden:

   ./fbs/singular_vector_transformation






.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="System Equivalent Model Mixing">

.. only:: html

    .. figure:: ./data/semm_scheme.svg   
       :target: ./fbs/system_equivalent_model_mixing.html

       System Equivalent Model Mixing

.. raw:: html

    </div>

.. toctree::
   :hidden:

   ./fbs/system_equivalent_model_mixing


|
|
|
|
|
|
|
|
|

One can distinguish between coupling and decoupling of dynamic systems as follows:

* **Coupling** is the process of assembly sub-systems by imposing physical boundary conditions to the common interface.
* **Decoupling** aims at identifying a standalone sub-system from the assembly by removing the influence of the other subsystem exerted through the interface connection.

In any case, the dynamic interaction between sub-systems (or substructures) is confined to a set of interface DoFs.
Let's consider the linearized equations of motion of a system composed by :math:`n` substructures in the frequency domain [1]_:

.. math::
    \mathbf{Z}(\omega)\,\boldsymbol{u}(\omega)=\boldsymbol{f}(\omega)+\boldsymbol{g}(\omega)

where :math:`\mathbf{Z}(\omega)` represent the block-diagonal frequency-dependent dynamic stiffness matrix of individual subsystems impedances,
the vector :math:`\boldsymbol{u}(\omega)` represents the displacements to the external force vector :math:`\boldsymbol{f}(\omega)`,
applied to the assembly, and :math:`\boldsymbol{g}(\omega)` is the vector of reaction forces
between the substructures:

.. math::
    \mathbf{Z}(\omega)=\begin{bmatrix}\mathbf{Z}^1(\omega) & & \mathbf{0}\\ & \ddots &  \\ \mathbf{0} & &\mathbf{Z}^n(\omega) \end{bmatrix}, \quad
    \boldsymbol{u}(\omega)=\begin{bmatrix}\boldsymbol{u}^1(\omega) \\ \vdots \\ \boldsymbol{u}^n(\omega) \end{bmatrix}, \quad
    \boldsymbol{f}(\omega)=\begin{bmatrix}\boldsymbol{f}^1(\omega) \\ \vdots \\ \boldsymbol{f}^n(\omega) \end{bmatrix}, \quad
    \boldsymbol{g}(\omega)=\begin{bmatrix}\boldsymbol{g}^1(\omega) \\ \vdots \\ \boldsymbol{g}^n(\omega) \end{bmatrix}

Let's now assume two interacting subsystems to be assembled:

.. figure:: ./data/assembly.svg
   :width: 300px
   :align: center

Considering a partition between internal :math:`(\star)_1,(\star)_3` and interface :math:`(\star)_2` DoFs, the vector of displacements, forces and reaction forces can be written as:

.. math::
    \boldsymbol{u}=\begin{bmatrix}
    \boldsymbol{u}_1^{\mathrm A} \\ \boldsymbol{u}_2^{\mathrm A}  \\ \boldsymbol{u}_2^{\mathrm B} \\ \boldsymbol{u}_3^{\mathrm B} 
    \end{bmatrix}, \quad 
    \boldsymbol{f}=\begin{bmatrix}
    \boldsymbol{f}_1^{\mathrm A} \\ \boldsymbol{f}_2^{\mathrm A}  \\ \boldsymbol{f}_2^{\mathrm B} \\ \boldsymbol{f}_3^{\mathrm B} 
    \end{bmatrix}, \quad
    \boldsymbol{g}=\begin{bmatrix}
    \boldsymbol{0} \\ \boldsymbol{g}_2^{\mathrm A}  \\ \boldsymbol{g}_2^{\mathrm B} \\ \mathbf{0}
    \end{bmatrix}.   

In order to assembly/disassembly the subsystems' dynamics, two physical conditions must be enforced:

* compatibility of displacements,
* equilibrium of forces.

Compatibility of displacements
******************************

The first interface condition to be fulfilled is the compatibility of displacements at the matching interface DoFs of the two subsystems:

.. math::
    \boldsymbol{u}^{\mathrm A}_{2}=\boldsymbol{u}^{\mathrm B}_{2}

This can be re-formulated by operating on the full set of physical DoFs :math:`\boldsymbol{u}` as:

.. math::
    \mathbf{B}\,\boldsymbol{u}=\mathbf{0}; \quad \mathbf{B}=
    \begin{bmatrix}
    \mathbf{0} & -\mathbf{I} &  \mathbf{I} & \mathbf{0} 
    \end{bmatrix}.

The signed Boolean matrix :math:`\mathbf{B}` maps the corresponding matching interface DoFs. Each row of the matrix identifies a single pair of interface DoFs to be connected.
Alternatively, by substituting the physical coordinates :math:`\boldsymbol{u}` with a set of unique generalized coordinates :math:`\boldsymbol{q}`:

.. math::
    \boldsymbol{u}=\mathbf{L}\boldsymbol{q} \implies \begin{cases}\boldsymbol{u}^{\mathrm A}_{1}=\boldsymbol{q}_{1} \\ \boldsymbol{u}^{\mathrm A}_{2}=\boldsymbol{q}_{2} \\ \boldsymbol{u}^{\mathrm B}_{2}=\boldsymbol{q}_{2} \\ \boldsymbol{u}^{\mathrm B}_{3}=\boldsymbol{q}_{3}  \end{cases}; \quad \mathbf{L}=
    \begin{bmatrix}
    \mathbf{I} & \mathbf{0} &  \mathbf{0} \\ 
    \mathbf{0} & \mathbf{I} & \mathbf{0} \\
    \mathbf{0} & \mathbf{I} &  \mathbf{0} \\
    \mathbf{0} & \mathbf{0} &  \mathbf{I} 
    \end{bmatrix}.

The localization Boolean matrix :math:`\mathbf{L}` maps the physical DoFs of all subsystems to the generalized global set :math:`\boldsymbol{q}`.

.. tip::
    The use of this coordinate transformation should remind you of a common finite element assembly procedure.

By using a unique set of coordinates :math:`\boldsymbol{q}`, it is made implicit that the compatibility of displacements for :math:`\boldsymbol{u}` is automatically satisfied:

.. math::
    \mathbf{B}\,\boldsymbol{u}=\mathbf{B}\,\mathbf{L}\,\boldsymbol{q}=\mathbf{0} \quad \forall \boldsymbol{q}

This means that :math:`\mathbf{B}` and :math:`\mathbf{L}` are each other's nullspaces:

.. math::
    \begin{array}{lcc}
    \mathbf{L}=\text{null}(\mathbf{B}), \\ \mathbf{B}^{\mathrm T}=\text{null}(\mathbf{L}^{\mathrm T}).\end{array}

Equilibrium of forces
*********************

The second conditions requires the force equilbrium at matching interface DoFs to be satisfied according to the *actio et reactio* principle:

.. math::
    \boldsymbol{g}^{\mathrm A}_{2}=-\boldsymbol{g}^{\mathrm B}_{2}

By back-projecting the vector of reaction forces :math:`\boldsymbol{g}` to the Boolean localization space :math:`\mathbf{L}`, the interface forces are directly paired:

.. math::
    \begin{array}{lcc}
    \mathbf{L}^{\mathrm T}\boldsymbol{g}=\mathbf{0} \implies \begin{cases} \boldsymbol{g}^{\mathrm A}_{1}=\mathbf{0} \\ \boldsymbol{g}^{\mathrm A}_{2}+\boldsymbol{g}^{\mathrm B}_{2}=\mathbf{0} \\ \boldsymbol{g}^{\mathrm B}_{3}=\mathbf{0} \end{cases}\end{array}.

Alternatively, by using the signed Boolean matrix :math:`\mathbf{B}` the reaction forces :math:`\boldsymbol{g}` are replaced by a set of Lagrange multipliers :math:`\boldsymbol{\lambda}`, 
which represent the intensity of the interface forces:

.. math::
    \begin{array}{lcc}
    \boldsymbol{g}=-\mathbf{B}^{\mathrm T}\boldsymbol{\lambda} \implies \begin{cases} \boldsymbol{g}^{\mathrm A}_{1}=\mathbf{0} \\ \boldsymbol{g}^{\mathrm A}_{2}= \boldsymbol{\lambda} \\\boldsymbol{g}^{\mathrm B}_{2}=-\boldsymbol{\lambda} \\ \boldsymbol{g}^{\mathrm B}_{3}=\mathbf{0} \end{cases}
    \end{array}

Using the definition of Lagrange multipliers for the interface forces automatically satisfies the equilibrium condition. 
This can be verified by exploiting the mathematical relationship between :math:`\mathbf{L}` and :math:`\mathbf{B}`:

.. math::
    \mathbf{L}^{\mathrm T}\boldsymbol{g}=-\mathbf{L}^{\mathrm T}\mathbf{B}^{\mathrm T}\boldsymbol{\lambda}=\mathbf{0} \quad \forall \,\boldsymbol{g}

Combining the equation of motion with the introduced interface conditions, the frequency-based formulation of the substructuring problem becomes:

.. math::
    \begin{cases}
    \mathbf{Z}(\omega)\,\boldsymbol{u}(\omega)=\boldsymbol{f}(\omega)+\boldsymbol{g}(\omega) \\ \mathbf{B}\,\boldsymbol{u}(\omega)=\mathbf{0} \\ \mathbf{L}^{\mathrm T}\boldsymbol{g}(\omega)=\mathbf{0}
    \end{cases}

From here on the frequency-dependence will be omitted for simplicity.
Solving the above equations of motion could be expensive due to the interface unknown to be resolved, i.e. :math:`\boldsymbol{u}` and :math:`\boldsymbol{g}`. 
Hence, the primal and dual formulations:

* **Primal**: satisfying a priori compatibility and solving for a unique set of interface displacements.
* **Dual**: satisfying a priori the equilibrium condition and solving for a new set of interface forces.

.. note::
    The concept of primal vs dual has roots in the optimization field. According to the duality principle, 2 mirrored formulation of the same optimization problem exist: 
    each variable in the primal form becomes a constraint in the dual and viceversa; 
    therefore, the objective direction is inversed (maximization and minimization for primal and dual respectively).

A primal formulation
********************

Primal assembly
===============

The primal assembly starts by defining a unique set of generalized coordinates :math:`\boldsymbol{q}`. 
The physical DoFs of all subsystems are mapped to :math:`\boldsymbol{q}` by applying the appropriate localization Boolean matrix :math:`\mathbf{L}`  
as :math:`\boldsymbol{u} = \mathbf{L}\,\boldsymbol{q}`. The compatibility is thus satisfied a priori. The equations of motion of the substructuring problem become:

.. math::
    \begin{cases}
    \mathbf{Z}^{\mathrm{A|B}}\mathbf{L}\boldsymbol{q}=\boldsymbol{f}+\boldsymbol{g} \\ \mathbf{L}^{\mathrm T}\boldsymbol{g}=\mathbf{0}
    \end{cases}, \qquad \mathbf{Z}^{\mathrm{A|B}}=\begin{bmatrix}
    \mathbf{Z}^\mathrm A & \mathbf{0} \\ 
    \mathbf{0} & \mathbf{Z}^\mathrm B
    \end{bmatrix}

The system is solved by premultiplying the first row by :math:`\mathbf{L}^\text{T}`, thus eliminating the interface forces and solving for the generalized interface displacement :math:`\boldsymbol{q}`:

.. math::
    \mathbf{\tilde{Z}}^{\mathrm{AB}}\boldsymbol{q}=\boldsymbol{p}, \qquad \mathbf{\tilde{Z}}^{\mathrm{AB}}=\mathbf{L}^T\mathbf{Z}^{\mathrm{A|B}}\mathbf{L} \quad \text{and}  \quad 
    \boldsymbol{p}=\mathbf{L}^\mathrm T\boldsymbol{f}

where :math:`\mathbf{\tilde{Z}}^\mathrm{AB}` is the primally assembled impedance for the generalized DoFs. 
Note that the size of the matrix is reduced according to the number of generalized DoFs being considered.
This assembly procedure is analogous to a finite element assembly. 
Here, the subsystems act like the super-elements and their dynamic properties are 'added' through the matching interface DoFs.

Consider the simple system depicted bellow. The primally assembled impedance can be written as follows:

.. figure:: ./data/assembly.svg
   :width: 300px
   :align: center

.. math::
    \mathbf{\tilde{Z}}^\mathrm{AB}=\mathbf{L}^\mathrm T\mathbf{Z}^\mathrm{A|B}\mathbf{L} \implies 
    \begin{bmatrix}\mathbf{Z}^\mathrm A_{11} & \mathbf{Z}^\mathrm A_{12} & \mathbf{0} \\ 
    \mathbf{Z}^\mathrm A_{21} & \mathbf{Z}^\mathrm A_{22}+\mathbf{Z}^\mathrm B_{22}& \mathbf{Z}^\mathrm B_{23} \\
    \mathbf{0} & \mathbf{Z}^\mathrm B_{32} &\mathbf{Z}^\mathrm B_{33}\end{bmatrix}= \mathbf{L}^\mathrm T\begin{bmatrix}\mathbf{Z}^\mathrm A_{11} & \mathbf{Z}^\mathrm A_{12} & \mathbf{0}& \mathbf{0} \\ 
    \mathbf{Z}^\mathrm A_{21} & \mathbf{Z}^\mathrm A_{22}  & \mathbf{0}& \mathbf{0} \\ 
    \mathbf{0} & \mathbf{0} & \mathbf{Z}^\mathrm B_{22} &\mathbf{Z}^\mathrm B_{23} \\
    \mathbf{0} & \mathbf{0} & \mathbf{Z}^\mathrm B_{32} &\mathbf{Z}^\mathrm B_{33} 
    \end{bmatrix} \mathbf{L}

with the localization matrix :math:`\mathbf{L}` defined as follows:

.. math::
    \qquad \mathbf{L}= \begin{bmatrix}
    \mathbf{I} & \mathbf{0} &  \mathbf{0} \\ 
    \mathbf{0} & \mathbf{I} & \mathbf{0} \\
    \mathbf{0} & \mathbf{I} &  \mathbf{0} \\
    \mathbf{0} & \mathbf{0} &  \mathbf{I} 
    \end{bmatrix}.

Primal disassembly
==================

A substructure decoupling procedure consists in the removal of the dynamic influence of a subsystem from the assembly in order to retrieve the remaining subsystem. 
In that sense, it can be considered as the 'reverse' operation of substructure coupling.
The reference representation of a disassembly procedure with corresponding DoFs is depicted:

.. figure:: ./data/disassembly.svg
   :width: 400px
   :align: center

Considering a partition between internal :math:`(\star)_1,(\star)_3` and interface :math:`(\star)_2` DoFs, the vector of displacements, forces and reaction forces can be written as:

.. math::
    \boldsymbol{u}=\begin{bmatrix}
    \boldsymbol{u}^\mathrm {AB}_{1} \\ \boldsymbol{u}^\mathrm {AB}_{2} \\ \boldsymbol{u}^\mathrm {AB}_{3} \\ \boldsymbol{u}^\mathrm {A}_{1} \\ \boldsymbol{u}^\mathrm {A}_{2}
    \end{bmatrix}, \quad 
    \boldsymbol{f}=\begin{bmatrix}
    \boldsymbol{0} \\ \boldsymbol{f}^\mathrm {AB}_{2} \\ \boldsymbol{f}^\mathrm {AB}_{3} \\ \boldsymbol{0} \\ \mathbf{0}
    \end{bmatrix}, \quad
    \boldsymbol{g}=\begin{bmatrix}
    \boldsymbol{0} \\ \boldsymbol{g}_2^\mathrm {AB} \\ \boldsymbol{0} \\ \mathbf{0} \\ -\boldsymbol{g}_2^\mathrm {A} 
    \end{bmatrix}.

The decoupling operation in the primal domain can be formulated as a coupling between the assembled system dynamic 
stiffness :math:`\mathbf{Z}^\text{AB}` and the negative dynamic stiffness of :math:`\mathbf{Z}^\text{A}`. The minus sign represents the subtraction operation:

.. math::

    \mathbf{\tilde{Z}}^\mathrm{B}=\mathbf{L}^\mathrm T\mathbf{Z}^\mathrm{AB|A}\mathbf{L} \implies 
    \begin{bmatrix}
    \cdot & \cdot & \cdot & \cdot\\ 
    \cdot & \mathbf{Z}^\mathrm B_{22} & \mathbf{Z}^\mathrm B_{23} & \cdot\\
    \cdot & \mathbf{Z}^\mathrm B_{32} &\mathbf{Z}^\mathrm B_{33} & \cdot \\
    \cdot & \cdot & \cdot & \cdot
    \end{bmatrix}=
    \mathbf{L}^T\begin{bmatrix}\cdot & \cdot & \cdot & \cdot& \cdot \\ 
    \cdot & \mathbf{Z}^\mathrm{AB}_{22} & \mathbf{Z}^\mathrm{AB}_{23} & \cdot& \mathbf{0} \\ 
    \cdot & \mathbf{Z}^\mathrm{AB}_{32} & \mathbf{Z}^\mathrm{AB}_{33} & \cdot& \mathbf{0} \\ 
    \cdot & \cdot & \cdot & \cdot & \cdot \\ 
    \cdot & \mathbf{0} & \mathbf{0} & \cdot  & -\mathbf{Z}^\mathrm A_{22} \\ 
    \end{bmatrix}\mathbf{L}

with the localization matrix :math:`\mathbf{L}` defined as:

.. math::
    \qquad \mathbf{L}= \begin{bmatrix}
    \cdot & \cdot &  \cdot& \cdot  \\ 
    \cdot & \mathbf{I} & \mathbf{0}& \cdot  \\
    \cdot & \mathbf{0} &  \mathbf{I}& \cdot  \\
    \cdot & \cdot &  \cdot&\cdot  \\
    \cdot & \mathbf{I} &  \mathbf{0}& \cdot 
    \end{bmatrix}

A dual formulation
******************

Dual assembly
===============

Starting from the general formulation of the substructuring problem, 
the dual approach chooses Lagrange multipliers :math:`\boldsymbol{\lambda}` as set of coupling forces according to the relation :math:`\boldsymbol{g} = - \mathbf{B}^\text{T} \boldsymbol{\lambda}` [2]_. 
The equilibrium is thus satisfied a priori. The equations of motion of the substructuring problem become:

.. math::
    \begin{cases}
    \mathbf{Z}^\mathrm{A|B}\boldsymbol{u}=\boldsymbol{f}-\mathbf{B}^\mathrm{T}\boldsymbol{\lambda} \\ \mathbf{B}\,\boldsymbol{u}=\mathbf{0} 
    \end{cases}

This is often written in a symmetrical form as:

.. math::
    \begin{bmatrix}
    \mathbf{Z}^\mathrm{A|B} & \mathbf{B}^\mathrm T \\ 
    \mathbf{B} & \mathbf{0}
    \end{bmatrix}\begin{bmatrix}
    \boldsymbol{u} \\ \boldsymbol{\lambda} 
    \end{bmatrix}=\begin{bmatrix}
    \boldsymbol{f} \\ \mathbf{0}
    \end{bmatrix}

To solve the system equations, let's first write them in the admittance notation:

.. math::
    \begin{cases}
    \boldsymbol{u}=\mathbf{Y}^\mathrm{A|B}(\boldsymbol{f}-\mathbf{B}^\mathrm{T}\boldsymbol{\lambda}) \\ \mathbf{B}\,\boldsymbol{u}=\mathbf{0} 
    \end{cases}, \qquad \mathbf{Y}^\mathrm{A|B}=\begin{bmatrix}
    \mathbf{Y}^\mathrm A & \mathbf{0} \\ 
    \mathbf{0} & \mathbf{Y}^\mathrm B
    \end{bmatrix}

By substituting the first line in the second line (compatibility constraint) and solving for :math:`\boldsymbol{\lambda}`:

.. math::
    \boldsymbol{\lambda}=\left(\mathbf{B}\mathbf{Y}^\mathrm{A|B}\mathbf{B}^\mathrm{T}\right)^{-1}\mathbf{B}\mathbf{Y}^\mathrm{A|B}\boldsymbol{f}

.. note::
    Interpretation [3]_ [4]_:

    * A displacement gap :math:`\boldsymbol{\Delta{u}}=\mathbf{B}\mathbf{Y}^\mathrm{A|B}\boldsymbol{f}` is formed between the still uncoupled subsystems' interface as a result of the applied excitation :math:`\boldsymbol{f}`.         
    * The interface forces :math:`\boldsymbol{\lambda}`, defined by the Lagrange multipliers, are applied in order to close the gap and keep the subsystems together.         
    * The stiffness operator between the force and the gap is named interface dynamic stiffness and is obtained by inverting the so called interface flexibility matrix :math:`\mathbf{Y}_\mathrm{int}=\mathbf{B}\mathbf{Y}^\mathrm{A|B}\mathbf{B}^\mathrm{T}`.

By substituting back the :math:`\boldsymbol{\lambda}` in the first line of the governing equation of motion:

.. math::
    \boldsymbol{u}=\mathbf{Y}^\mathrm{A|B}(\boldsymbol{f}-\mathbf{B}^\mathrm{T}\boldsymbol{\lambda}) \implies \boldsymbol{u}=\mathbf{Y}^\mathrm{A|B}\boldsymbol{f}-\mathbf{Y}^\mathrm{A|B}\mathbf{B}^\mathrm{T}\left(\mathbf{B}\mathbf{Y}^\mathrm{A|B}\mathbf{B}^\mathrm{T}\right)^{-1}\mathbf{B}\mathbf{Y}^\mathrm{A|B}\boldsymbol{f}

Due to its derivation, this formulation is referred to as **Lagrange multipliers - frequency based substructuring (LM-FBS)**.

.. note::
    This final result can be interpreted as follows:

    * If an external force :math:`\boldsymbol{f}` acts on the uncoupled system, a response :math:`\boldsymbol{u}_\text{uncoupled} = \mathbf{Y}^\mathrm{A|B} \boldsymbol{f}` is produced and creates an interface incompatibility gap :math:`\Delta \boldsymbol{u}`.
    * Interface forces :math:`\boldsymbol{\lambda}` close this gap, and thus cause response of the system :math:`\boldsymbol{u}_\text{coupled} = \mathbf{Y}^\mathrm{A|B}\mathbf{B}^\mathrm{T}\left(\mathbf{B}\mathbf{Y}^\mathrm{A|B}\mathbf{B}^\mathrm{T}\right)^{-1}\mathbf{B}\mathbf{Y}^\mathrm{A|B}\boldsymbol{f}`
    * The dually assembled response :math:`\boldsymbol{u}` is the combination of the response to external forces exciting the uncoupled system and to the interface forces needed to keep the substructures together while being excited by :math:`\boldsymbol{f}`.

The dually assembled admittance is written as:

.. math::
    \mathbf{\hat{Y}}^\mathrm{AB}=\left[\mathbf{I}-\mathbf{Y}^\mathrm{A|B}\mathbf{B}^\mathrm{T}\left(\mathbf{B}\mathbf{Y}^\mathrm{A|B}\mathbf{B}^\mathrm{T}\right)^{-1}\mathbf{B}\right]\mathbf{Y}^\mathrm{A|B}

The dually assembled admittance :math:`\mathbf{\hat{Y}}^\mathrm{AB}` contains twice the interface DoFs and has the same size of the uncoupled admittance :math:`\mathbf{Y}^\mathrm{A|B}`. 
The redundant DoFs may be removed when deemed necessary. 
This can be performed by exploiting the primal relations :math:`\boldsymbol{u}=\mathbf{L}\,\boldsymbol{q}` and :math:`\boldsymbol{p}=\mathbf{L}^\mathrm{T}\,\boldsymbol{f}` 
to obtain the primally assembled admittance:

.. math::
    \mathbf{\tilde{Y}}^\mathrm{AB}=\mathbf{L}^+\mathbf{\hat{Y}}^\mathrm{AB}(\mathbf{L}^\mathrm T)^+

where :math:`(\star)^+` denotes a pseudo-inversion.

.. note::

    Note how the interface problem boils down to solving the linear problem :math:`\mathbf{Y}_\mathrm{int}\boldsymbol{\lambda}=\boldsymbol{\Delta{u}}`, 
    where the :math:`\mathbf{Y}_\mathrm{int}` contains the core of the assembly operation:

    .. math::
        \mathbf{Y}_\mathrm{int}=\mathbf{B}\mathbf{Y}^\mathrm{A|B}\mathbf{B}^\mathrm{T}=\mathbf{Y}^\mathrm{A}_{22}+\mathbf{Y}^\mathrm{B}_{22}
    
    Therefore, according to the LM-FBS definition, the zeros of the sum of the interface FRFs of the subsystems to be coupled, correspond to the resonances of the assembled system.

Dual disassembly
================

The dual decoupling problem consists of finding the interface forces that suppress the influence of A on AB, thus isolating the uncoupled response of subsystem B [5]_.
Following the definition of decoupling, the equilibrium condition states that the interface forces that ensure the compatibility act in opposite direction on the assembled system AB.
The decoupling can finally be formulated as a standard coupling procedure with a negative admittance for the system to be disassembled:

.. math::

    \begin{cases}
    \boldsymbol{u}=\mathbf{Y}^\mathrm{AB|A}(\boldsymbol{f}-\mathbf{B}^\mathrm{T}\boldsymbol{\lambda}) \\ \mathbf{B}\,\boldsymbol{u}=\mathbf{0} 
    \end{cases}, \qquad \mathbf{Y}^\mathrm{AB|A}=\begin{bmatrix}
    \mathbf{Y}^\mathrm{AB} & \mathbf{0} \\ 
    \mathbf{0} & -\mathbf{Y}^\mathrm A
    \end{bmatrix}

Analogously to the coupling formulation, the interface problem is solved as follows:

.. math::

    \boldsymbol{\lambda}=\left(\mathbf{B}\mathbf{Y}^\mathrm{AB|A}\mathbf{B}^\mathrm{T}\right)^{-1}\mathbf{B}\mathbf{Y}^\mathrm{AB|A}\boldsymbol{f}

.. note::
    Interpretation: 

    The interface flexibility matrix :math:`\mathbf{Y}_\mathrm{int}=\mathbf{B}\mathbf{Y}^\mathrm{AB|A}\mathbf{B}^\mathrm{T}` relates the unknown 
    multipliers :math:`\boldsymbol{\lambda}` with the interface displacements :math:`\boldsymbol{u}_\mathrm{int}=\mathbf{B}\mathbf{Y}^\mathrm{AB|A}\boldsymbol{f}`. 
    Let's consider the separate contributions of the dynamics of AB and A. 
    Given an external excitation :math:`\boldsymbol{f}^\mathrm{AB}`, a displacement 
    :math:`\boldsymbol{u}_\mathrm{int}=\mathbf{B}^\mathrm{AB}\mathbf{Y}^\mathrm{AB}\boldsymbol{f}^\mathrm{AB}` occurs at the interface as a result of the combined 
    dynamics of A and B. 
    The connection forces :math:`\boldsymbol{\lambda}=\left(\mathbf{B}^\mathrm{AB}\mathbf{Y}^\mathrm{AB}\mathbf{B}^\mathrm{{AB}^{T}}-\mathbf{B}^\mathrm{A}\mathbf{Y}^\mathrm{A}\mathbf{B}^\mathrm{{A}^{T}}\right)^{-1}\boldsymbol{u}_\mathrm{int}`
    are then estimated to compensate for the dynamic contribution of the subsystem A, thus removing its influence on :math:`\mathbf{B}`.

By solving according to the LM-FBS:

.. math::
    \mathbf{\hat{Y}}^\mathrm{B}=\left[\mathbf{I}-\mathbf{Y}^\mathrm{AB|A}\mathbf{B}^\mathrm{T}\left(\mathbf{B}\mathbf{Y}^\mathrm{AB|A}\mathbf{B}^\mathrm{T}\right)^{-1}\mathbf{B}\right]\mathbf{Y}^\mathrm{AB|A}

.. tip::
    By extending the decoupling interface from :math:`n_2` interface DoFs to the :math:`n_1` internal DoFs, the interface observability and controllability is
    improved, which is beneficial for the efficiency of the decoupling procedure.

Decoupling offers a broader amount of potentially matching DoFs with respect to coupling. 
However, while the interface measurements play the core role in the substructuring process, the internal DoFs, theoretically, do not bring anything new to the game. 
Afterall, the interface decoupling problem to be solved remained :math:`\mathbf{Y}_\mathrm{int}\boldsymbol{\lambda}=\boldsymbol{u}_\mathrm{int}`. 

In real life, however, erroneous modeling of the interface dynamics and measurement errors are the daily bread for experimentalists. 
The use of additional (internal) information between AB and A can help improving the observability, controllability and conditioning of the interface problem.
The interface strategies commonly used are [6]_ [7]_:

- Standard interface: :math:`n_c=n_e=n_2`. The interface matrix is square and full rank. Compatibility and equilibrium are enforced at the interface only,
- Extended interface: :math:`n_c=n_e=n_2+n_1`. The interface matrix is square and (without measurement errors) singular. In practice, modeling and measurement errors affect the measurements. The interface matrix is ill-conditioned and a singular value truncation is often performed. The additional internal compatibility and equilibrium constraints contributes to increase the observability and controllability of the interface dynamics.
- Non-collocated overdetermined: :math:`n_c=n_2+n_1,\,n_e=n_2`. The compatibility condition is extended to the internal DoFs and the linear problem :math:`\mathbf{Y}_\text{int}\boldsymbol{\lambda}=\boldsymbol{u}_\text{int}` is overdetermined. Since measurements are not perfect, there is no exact solution for the unknown multipliers :math:`\boldsymbol{\lambda}`. The optimal solution is found via pseudo-inverse in a least squares sense.

Following these considerations, a generalized version of the LM-FBS is written as follows:

.. math::
    \mathbf{\hat{Y}}^\mathrm{B}=\left[\mathbf{I}-\mathbf{Y}^\mathrm{AB|A}\mathbf{B}_\mathrm f^\mathrm{T}\left(\mathbf{B}_\mathrm u\mathbf{Y}^\mathrm{AB|A}\mathbf{B}_\mathrm f^\mathrm{T}\right)^{+}\mathbf{B}_\mathrm u\right]\mathbf{Y}^\mathrm{AB|A}

where potentially different Boolean matrices can be used for the compatibility :math:`\mathbf{B}_\mathrm u` and equilibrium :math:`\mathbf{B}_\mathrm f` conditions 
and :math:`(\star)^+` denotes a pseudo-inversion.

.. rubric:: References

.. [1] de Klerk, D., Rixen, D. J., Voormeeren,S. (2008). General Framework for Dynamic Substructuring: History, Review and Classification of Techniques. In: AIAA Journal 46.5.
.. [2] De Klerk, D., Rixen, D. J., De Jong, J. (2006) The frequency based substructuring (FBS) method reformulated according to the dual domain decomposition method. In: 24th International Modal Analysis Conference. St.Louis, MO .
.. [3] Tiso, P., Allen, M. S., Rixen, D., Abrahamsson, T., Van der Seijs, M., Mayes, R. L. (2020) Substructuring in Engineering Dynamics - Emerging Numerical and Experimental Techniques. Springer .
.. [4] Rixen, D., Godeby, T., Pagnacco, E. (2006) Dual Assembly of substructures and the FBS Method: Application to the Dynamic Testing of a Guitar. In: International Conference on Noise and Vibration Engineering, ISMA. KUL. Leuven, Belgium, Sept.
.. [5] van der Seijs, M. V. (2016) Experimental dynamic substructuring: Analysis and design strategies for vehicle development. Delft University of Technology.
.. [6] Voormeeren, S., Rixen, D. (2012) A family of substructure decoupling techniques based on a dual assembly approach. In: Mechanical Systems and Signal Processing 27, pp. 379– 396 .
.. [7] D’Ambrogio, W. and Fregolent, A. (2011) Direct decoupling of substructures using primal and dual formulation. In: Proceedings of 29th IMAC, a Conference on Structural Dynamics, pp. 47–76 .Basic examples
==============

This examples show how to use basic features of pyFBS. Explore this basic examples to get familiar with the pyFBS workflow.

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Static display">

.. only:: html

    .. figure:: ./data/interaction.gif	   
       :target: ./basic_examples/01_static_display.html

       Static display

.. raw:: html

    </div>

.. toctree::
   :hidden:

   ./basic_examples/01_static_display


  
 
.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Interactive positioning">

.. only:: html

    .. figure:: ./data/snapping.gif	   
       :target: ./basic_examples/02_interactive_display.html

       Interactive positioning

.. raw:: html

    </div>

.. toctree::
   :hidden:

   ./basic_examples/02_interactive_display
   


   
.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="FRF synthetization">

.. only:: html

    .. figure:: ./data/FRF_syn-FRF-visualization.svg 
       :target: ./basic_examples/03_FRF_synthetization.html

       FRF synthetization

.. raw:: html

    </div>

.. toctree::
   :hidden:

   ./basic_examples/03_FRF_synthetization============================
Virtual Point Transformation
============================
Virtual point transformation (VPT) projects measured dynamics (input and output signals) into a subspace composed by the predefined interface deformation modes (IDMs) [1]_. 
By default, only 6 rigid IDMs are used in the transformation, thus retaining only the dynamics loading the surrounded interface in a purely rigid manner. 
Rigid IDMs can also be extended by the user with flexible interface modes.
Current implementation of the :class:`pyFBS.VPT` additionaly supports the expansion where directly measured rotational response is included in the transformation [2]_. 

.. note:: 
   Download example showing the basic use of the VPT: :download:`04_VPT.ipynb <../../../examples/04_VPT.ipynb>`
   
Consider an example for the VPT where 9 impacts and 9 channels (3 tri-axial accelerometers) are positioned around the interface:

.. figure:: ./../data/pic_vpt.png
   :width: 700px

.. tip::
	For the VP with 6 DoFs it is common practice to use 9 impacts and 9 channels around the interface. 
	Larger number of impacts/channels is encouraged but the required experimental effort is increased.
	By increasing the number of measurements over the VP DoFs the benefits are two-hand: 
	observability and controlability of the interface is better, and the measurement errors 
	(bias errors due to sensor misalignement and uncorrelated measurement noise) are filtered our from the VP to some extend. 

   
***************************
Interface deformation modes
***************************

.. tip::
	It is advised impacts and channels are evenly distributed across the interface for better controlability and observability of the VP DoFs.
	Even distribution should be ensured in terms of positions and directions as well.
	Also ensure that impacts and channels do not point straight to the VP in order to generate rotational response and moment around the VP.

For the VPT, positional data is required for channels (``df_chn_up``), impacts (``df_imp_up``) and for virtual points (``df_vp`` and ``df_vpref``). 
Positioning can be performed using interactive positioning or imported from Excel datasheets.

.. tip::
	Placing sensors closer to VP reduces the effects of flexible IDMs. 
	However, with a decreased distance the uncertainties associated with the position and orientation are increased.
	Therefore, sensor should be mounted on the strcture with great care.
	Using accelerometer faces as impact locations is discouraged as sensor overload is quickly reached in this manner.

While performing VPT, special care should be put into determining ``Grouping`` number and DoFs ``Description``. 
The ``Description`` tells which DoFs you would like to reconstruct using VPT. 
``Grouping`` number organizes which impacts and sensors belong to which virtual point. 
Arbitrary number of VPs can be reconstructed at once if grouping numbers are properly defined.

.. figure:: ./../data/df_vp.PNG
   :width: 750px

After the positions are defined a class instance of :class:`pyFBS.VPT` can be created.

.. code-block:: python

	vpt = pyFBS.VPT(df_chn_up,df_imp_up,df_vp,df_vpref)


Interface displacement reduction :func:`pyFBS.VPT.define_IDM_U` and interface force reduction :func:`pyFBS.VPT.define_IDM_F` are directly defined. 
Both reduction matrices are avialable as class variables ``vpt.Tu`` and ``vpt.Tf``.


**********************
Application of the VPT
**********************
After the reduction matrices are defined the VPT can be applied directly on an FRF matrix.

.. code-block:: python

		vpt.apply_VPT(freq,FRF)

Transformed FRF matrix is then available as a class variable ``vpt.vptData``.

.. raw:: html

   <iframe src="../../_static/FRF_VP.html" height="500px" width="750px" frameborder="0"></iframe>

.. tip::
	Check driving-point FRFs if passivity criterion is met by simply clicking on corresponding circles.

Reciprocity check is also already implemented in pyFBS:

.. code-block:: python

	reciprocity = pyFBS.coh_on_FRF(vpt.vptData[:,:6,:6])

.. raw:: html

   <iframe src="../../_static/VP_rec.html" height="270px" width="350px" frameborder="0"></iframe>

******************************
Measurement quality indicators
******************************

One of the primary advantages of the VPT is also the ability to evaluate consistency of the performed measurements. 
Measurement consistency is evaluated by expanding the reduced virtual DoFs back to the original DoFs and comparing them with the measured ones.

.. code-block:: python

		vpt.consistency([1],[1])
		
If the interface would be perfectly rigid, the filtered response would be equal to the measured one. 
However, if the interface is not completely rigid or if predetermined positions and orientations are not perfect, the filtered response will vary from the measured response.

Both channel/sensor (``vpt.specific_sensor`` and ``vpt.overall_sensor``) and impact (``vpt.specific_impact`` and ``vpt.overall_impact``) consistency can be evaluated after the transformation. 

.. raw:: html

   <iframe src="../../_static/VP_spec_cons.html" height="290px" width="600px" frameborder="0"></iframe>

.. warning::
	A low value of quality indicators indicates inconsistent measurements. 
	Always check if sensor/impacts are correctly positioned, orientated, sensors' sensitivity is correct, responses do not exceed sensors range, etc. 
	Low values could also indicate flexible interface behaviour; in that case, check CMIF value [3]_ if rigid IDMs are in fact dominant.

.. tip::
	Take A LOT of pictures of your experiment in case you need to correct some sensor positions or orientations later.

.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!

.. rubric:: References

.. [1] de Klerk D, Rixen DJ, Voormeeren SN, Pasteuning P. Solving the RDoF problem in experimental dynamic sybstructuring. InInternational modal analysis conference IMAC-XXVI 2008 (pp. 1-9).
.. [2] Tomaž Bregar, Nikola Holeček, Gregor Čepon, Daniel J. Rixen, and Miha Boltežar. Including directly measured rotations in the virtual point transformation. Mechanical Systems and Signal Processing, 141:106440, July 2020.
.. [3] Allemang RJ, Brown DL. A complete review of the complex mode indicator function (CMIF) with applications. InProceedings of ISMA International Conference on Noise and Vibration Engineering, Katholieke Universiteit Leuven, Belgium 2006 Sep (pp. 3209-3246).==============================
Singular Vector Transformation
==============================
Singular Vector Transformation (SVT) consists in projecting the acquired data into subspaces composed by dominant singular vectors, which are extracted directly from the available FRF datasets. 
No geometrical and/or analytical model is required. 
If some basic requirements are met, the reduced orthonormal frequency dependent basis would be able to control and observe most of the rigid and flexible vibration modes of interest over a broad frequency range. 
The SVT can tackle challenging scenarios with flexible behaving interfaces and lightly damped systems.



.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Singular Vector Transformation">

.. only:: html

    .. figure:: ./../data/SVT.png 
       :target: 12_SVT.html

       Singular Vector Transformation

.. raw:: html

    </div>

.. toctree::
   :hidden:

   12_SVT





.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="FBS Decoupling with SVT">

.. only:: html

    .. figure:: ./../data/seven_three.png   
       :target: 13_SVT_decoupling.html

       FBS Decoupling with SVT

.. raw:: html

    </div>

.. toctree::
   :hidden:

   13_SVT_decoupling


|
|
|
|
|
|
|
|

Singular value decomposition
****************************

The singular value decomposition can be used to decompose and sort the dynamical information of the measured FRF dataset:

.. math::
    \mathbf{Y}(\omega)=\mathbf{U}(\omega)\mathbf{\Sigma}(\omega)\mathbf{V}(\omega)^\text{H}=\mathbf{U}(\omega)\textbf{diag}(\sigma_i(\omega))\mathbf{V}(\omega)^\text{H}

where :math:`\mathbf{U}(\omega)` and :math:`\mathbf{V}(\omega)` are the orthonormal frequency-dependent left and right singular vectors and :math:`\mathbf{\Sigma}(\omega)` 
is the diagonal matrix containing the non-negative singular values of the matrix.

The column vectors of :math:`\mathbf{U}(\omega)` and :math:`\mathbf{V}(\omega)` can be considered as approximate mode shapes and approximate modal participation 
factors at the frequency :math:`\omega`. The associated singular values are the scaling factors representing how dominant the contribution (**dynamically**) 
of a singular state is at that frequency.

.. note::
    Note how the singular 'motions' are frequency dependent and span an orthogonal space (mathematical orthogonality), 
    while the characteristic modal vectors are frequency independent and are :math:`\mathbf{K}`-orthogonal and :math:`\mathbf{M}`-orthogonal (physical orthogonality).

.. tip::
    **Idea!** Why don't we use the singular vectors as reduction bases to weaken the interface problem?

    A subset of dominant singular vectors in :math:`\mathbf{U}(\omega)` and :math:`\mathbf{V}(\omega)` could be used to construct meaningful reduction 
    spaces that naturally contain rigid as well as flexible local interface deformations.

Methodology
***********

An extended (direct) decoupling problem is considered, where the goal is finding the interface forces that suppress the influence of A on AB, thus isolating the uncoupled response of subsystem B.

.. figure:: ./../data/disassembly.svg
   :width: 400px
   :align: center

Let's assume collocation between the inputs between AB and A and outputs between AB and A. 
Let's also assume that inputs and outputs do not share the same location and are distributed all over the common area between AB and A (interface and internal) 
such that the interface deformation of interest is controlled and observed in a similar manner.

The proposed procedure is summarized [1]_:

1.  Place inputs and outputs in a collocated manner between AB and A. Distribute them homogeneously over the full area (interface and internal), 
    such that all the vibration modes of interest are controlled and observed through similar dynamic spaces.
2.  Perform a SVD on :math:`\mathbf{Y}^\text{A}`, such that :math:`\mathbf{Y}^\text{A}=\mathbf{U}^\text{A}\mathbf{\Sigma}^\text{A}(\mathbf{V}^\text{A})^\text{H}`.
3.  Select truncated basis of left- and right SVs :math:`\mathbf{U}_r^\text{A}`, :math:`\mathbf{V}_r^\text{A}` for the reduction 
    of measured displacements and forces. The subspaces' vectors are complex, orthonormal and frequency-dependent.
4.  Use the same reduction bases for :math:`\mathbf{Y}^\text{A}` and :math:`\mathbf{Y}^\text{AB}`. This is necessary for imposing collocated compatibility and equilibrium 
    between the substructures.
5.  Decouple the reduced components according to the LM-FBS algorithm.

.. note::
    The SVD can be applied analogously on :math:`\mathbf{Y}^\text{AB}`. For simplicity, the derivation will be performed with the singular bases of :math:`\mathbf{Y}^\text{A}`.

Mathematical derivation
***********************

Let's apply an SVD on the measured dataset :math:`\mathbf{Y}^\text{A}`:

.. math::
    \mathbf{Y}^\text{A}(\omega)=\mathbf{U}^\text{A}(\omega)\mathbf{\Sigma}^\text{A}(\omega)(\mathbf{V}^\text{A}(\omega))^\text{H}

The left and right singular subspaces are complex orthonormal frequency-dependent bases.
Proper reduced bases :math:`\mathbf{U}_r^\text{A}` and :math:`\mathbf{V}_r^\text{A}` are chosen per frequency line according to a user-defined criteria. 
The frequency dependence will be omitted for simplicity from here on.

Displacement reduction
======================

The same reduction basis :math:`\mathbf{U}_r^\text{A}` is chosen for A and AB to preserve a collocated compatibility in the reduced domain:

.. math::
    \begin{array}{lcc}
    \boldsymbol{u}^\text{A}=\mathbf{U}_r^\text{A}\boldsymbol{\zeta}^\text{A}+\boldsymbol{\mu}^\text{A} \\
    \boldsymbol{u}^\text{AB}=\mathbf{U}_r^\text{A}\boldsymbol{\zeta}^\text{AB}+\boldsymbol{\mu}^\text{AB}
    \end{array}

The residual :math:`\boldsymbol{\mu}^\text{A}` represent the low-value singular states not included in :math:`\mathbf{U}_r^\text{A}`. 
The residual :math:`\boldsymbol{\mu}^\text{AB}` indicates the motion of AB that could not be described by the reduced singular subspace of A.
The transformation from the measured to the reduced space is obtained as:

.. math::
    \boldsymbol{\zeta}=\mathbf{U}_r^+\boldsymbol{u}=\mathbf{U}_r^\text{H}\boldsymbol{u}

where

.. math::
    \boldsymbol{u}=\begin{bmatrix}
    \boldsymbol{u}^\text{AB}  \\ \boldsymbol{u}^\text{A}  
    \end{bmatrix}, \quad \boldsymbol{\zeta}=\begin{bmatrix}
    \boldsymbol{\zeta}^\text{AB}  \\ \boldsymbol{\zeta}^\text{A}  
    \end{bmatrix} \quad \text{and}\quad \mathbf{U}_r=\begin{bmatrix}
    \mathbf{U}_r^\text{A} & \mathbf{0} \\ \mathbf{0} & \mathbf{U}_r^\text{A} 
    \end{bmatrix}

The superscripts :math:`(\star)^+` and :math:`(\star)^\text{H}` denote the Moore-Penrose pseudo-inverse and Hermitian operators respectively.

.. note::
    Note that, due to the orthogonality of :math:`\mathbf{U}_r`, no mathematical inversion has to be computed. 

To check the quality of the transformation, the filtered and measured displacements can be compared:

.. math::
    \tilde{\boldsymbol{u}}=\mathbf{U}_r\mathbf{U}_r^\text{H}\boldsymbol{u}=\mathbf{F}_\text{u}\boldsymbol{u}.

Force reduction
===============

The same reduction basis :math:`\mathbf{V}_r^\text{A}` is chosen for A and AB to preserve a collocated equilibrium in the reduced domain:

.. math::
    \begin{array}{lcc}
    \boldsymbol{\eta}^\text{A}=(\mathbf{V}_r^\text{A})^\text{H}\boldsymbol{g}^\text{A}\\
    \boldsymbol{\eta}^\text{AB}=(\mathbf{V}_r^\text{A})^\text{H}\boldsymbol{g}^\text{AB}
    \end{array}

The transformation from the reduced to the measured space is obtained as:

.. math::
    \boldsymbol{g}=\begin{bmatrix}
    \boldsymbol{g}^\text{AB}  \\ \boldsymbol{g}^\text{A}  
    \end{bmatrix}, \quad \boldsymbol{\eta}=\begin{bmatrix}
    \boldsymbol{\eta}^\text{AB}  \\ \boldsymbol{\eta}^\text{A}  
    \end{bmatrix}\quad \text{and}\quad \mathbf{V}_r=\begin{bmatrix}
    \mathbf{V}_r^\text{A} & \mathbf{0} \\ \mathbf{0} & \mathbf{V}_r^\text{A} 
    \end{bmatrix}

The superscripts :math:`(\star)^+` and :math:`(\star)^\text{H}` denote the Moore-Penrose pseudo-inverse and Hermitian operators respectively.
The filtered forces can be compared to the measured ones:

.. math::
    \tilde{\boldsymbol{g}}=\mathbf{V}_r\mathbf{V}_r^\text{H}\boldsymbol{g}=\mathbf{F}_g\boldsymbol{g}

Finally, the compatibility and equilibrium conditions are imposed in the reduced space observed by :math:`\mathbf{U}_r` and controlled by :math:`\mathbf{V}_r` respectively:

.. math::
    \begin{cases} 
    \boldsymbol{u}=\mathbf{Y}^\text{AB|A}(\mathbf{f}-\mathbf{V}_r\mathbf{B}^\text{T}_{\eta}\boldsymbol{\lambda}_{\eta}) \\
    \mathbf{B}_{\zeta}\mathbf{U}_r^\text{H}\boldsymbol{u}=\mathbf{0}  
    \end{cases}

Thus, a weakening of the interface problem is obtained. The solution of the equations follows the LM-FBS formulation:

.. math::
    \boldsymbol{\zeta}=\left[\mathbf{I}-\mathbf{Y}^\text{AB|A}_{\zeta\eta}\mathbf{B}_{\eta}^\text{T}\left(\mathbf{B}_{\zeta}\mathbf{Y}^\text{AB|A}_{\zeta\eta}\mathbf{B}_{\eta}^\text{T}\right)^{-1}\mathbf{B}_{\zeta}\right]\mathbf{Y}^\text{AB|A}_{\zeta\eta}\boldsymbol{\eta}=\mathbf{Y}^\text{B}_{\zeta\eta}\boldsymbol{\eta}

where the reduced uncoupled admittance :math:`\mathbf{Y}^\text{AB|A}_{\zeta\eta}` can be written as:

.. math::
    \mathbf{Y}^\text{AB|A}_{\zeta\eta}=\mathbf{U}_r^\text{H}\mathbf{Y}^\text{AB|A}\mathbf{V}_r=\begin{bmatrix} \mathbf{Y}^\text{AB}_{\zeta\eta}  & \mathbf{0} \\ \mathbf{0} & -\mathbf{Y}^\text{A}_{\zeta\eta}
    \end{bmatrix}=\begin{bmatrix} \mathbf{Y}^\text{AB}_{\zeta\eta}  & \mathbf{0} \\ \mathbf{0} & -\mathbf{\Sigma}_r^\text{A}
    \end{bmatrix}

Filtering :math:`\mathbf{Y}^\text{AB}` by reduction and back-projection (see filter matrices :math:`\mathbf{F}_\text{u}` and :math:`\mathbf{F}_\text{f}`) using the singular subspaces of :math:`\mathbf{Y}^\text{A}` 
can be useful to determine how well the local dynamics of AB is described through the reduced dominant controllability and observability spaces of A.

Requirements
************

Let's remind ourselves of the main assumptions for a successful SVT:

- Inputs in A and in the A part of AB must be collocated and so must the output in A and in the A part of AB. This is necessary to since no geometrical reduction is applied afterwards. In experimental practice, this requirement should not be too challenging,

- Inputs and output spaces should contain the same dynamical information. The balance between controllability and observability spaces is needed to minimize the physical inconsistency of the reduced model. In experimental practice, an homogeneous distribution of inputs and outputs throughout the available common area (A and the A part of AB) should be ensured,

- The reduced orthonormal subspaces of A should be able to map the observed and controlled dynamics in AB.

Benefits
********

The SVT relies on a simple concept: the measured dynamics is projected into a subspace composed by dominant singular states extracted from the measurement itself.
What are the main advantages and features of the SVT?

- No need for geometrical/analytical information,
- Flexible behaviour potentially included (if properly observed/controlled),
- Frequency-dependent transformation,
- Efficient least square smoothing of random error and outliers,
- Better conditioning of the interface problem,
- Reduce sensitivity to measurement location/direction bias.

.. rubric:: References

.. [1] F.Trainotti, T.Bregar, S.W.B.Klaassen and D.J.Rixen. Experimental Decoupling of Substructures by Singular Vector Transformation. in: Mechanical System and Signal Processing ('Under Review'), 2021############
VPT Coupling
############

The frequency-based substructure coupling is embedded in ``pyFBS``. 
In particular, the admittance-based dual formulation named Lagrange-Multiplier Frequency-Based Substructuring (LM-FBS) is implemented [1]_. 
In the following, a basic coupling of two numerically-generated substructures is presented. 
The virtual point transformation [2]_ is applied to impose collocated matching DoFs at the interface. 
This can also be performed analogously with experimentally acquired data.

.. note:: 
   Download example showing a substructure coupling application: :download:`07_coupling.ipynb <../../../examples/07_FBS_coupling.ipynb>`

.. tip::
    Why use virtual point when coupling substructures?
    
    * Complex interfaces are very often inaccessible for the measurement equipment. Therefore measurements are often performed away from the actual interface. Also, the collocation of DoFs at the neighboring interfaces is almost impossible to ensure. Virtual point, on the other hand, can be defined in an arbitrary location at the interface that coincides for all substructures.
    * Virtual point accounts for rotational degrees of freedom, which are mandatory for the successful coupling of the substructures. 
    * By the reduction of measurements to the virtual point the interface problem is weakened. That means that compatibility and equilibrium conditions on the measured DoFs are a bit more relaxed. In this manner, unwanted stiffening effects can be avoided. 
    * Due to the reduction, measurement errors such as bias from sensor positioning or uncorrelated measurement noise are filtered out to some extend.
         
    
Example Datasets and 3D view
****************************
Load the required predefined datasets and open the 3D viewer in the background as already shown in `3D Display <../../../html/examples/basic_examples/01_static_display.html>`_. 
Especially for the illustration of different substructures and the assembly, the 3D viewer subplot capabilities of `PyVista <https://docs.pyvista.org/index.html>`_ can be used.

.. code-block:: python

    view3D = pyFBS.view3D(show_origin = False, show_axes = False,shape =  (1,3),title = "Overview")
    
Add an STL file of substructure A to the 1-1 subplot and show the corresponding accelerometers, channels and impacts.

.. code-block:: python

    view3D.plot.subplot(0,0)
    view3D.plot.isometric_view()
    view3D.plot.add_text("A structure", position='upper_left', font_size=10, color="k", font="times", name="A_structure")

    view3D.add_stl(stl_dir_A,color = "#83afd2",name = "A");
    view3D.show_acc(df_acc_A)
    view3D.show_imp(df_imp_A)
    view3D.show_chn(df_chn_A)
    
.. figure:: ./../data/seven_one.png
   :width: 500px
    
    
Add an STL file of substructure B to the 1-2 subplot and show the corresponding accelerometers, channels and impacts.

.. code-block:: python

    view3D.plot.subplot(0,1)
    view3D.plot.isometric_view()
    view3D.plot.add_text("B structure", position='upper_left', font_size=10, color="k", font="times", name="B_structure")

    view3D.add_stl(stl_dir_B,color = "#83afd2",name = "B");
    view3D.show_acc(df_acc_B,overwrite = False)
    view3D.show_imp(df_imp_B,overwrite = False)
    view3D.show_chn(df_chn_B,overwrite = False)
 
.. figure:: ./../data/seven_two.png
   :width: 500px
   
    
Add an STL file of the assembly AB to the 1-2 subplot and show the corresponding reference accelerometers, channels and impacts.
 
.. code-block:: python

    view3D.plot.subplot(0,2)
    view3D.plot.isometric_view()
    view3D.plot.add_text("AB structure", position='upper_left', font_size=10, color="k", font="times", name="AB_structure");

    view3D.add_stl(stl_dir_AB,color = "#83afd2",name = "AB");
    view3D.show_acc(df_acc_AB,overwrite = False)
    view3D.show_imp(df_imp_AB,overwrite = False)
    view3D.show_chn(df_chn_AB,overwrite = False)
    
.. figure:: ./../data/seven_three.png
   :width: 500px
        
Each separate subplot view can also be linked or unlinked:

.. code-block:: python

    view3D.plot.link_views()
    #view3D.plot.unlink_views()

.. tip::
    With ``pyFBS`` you can simply prepare the experiments before hand! 
    Position your virtual accelerometers and sensors on 3D model, generate numerical FRFs and make sure that everything is in order. 
    Then, perform your experiment, following the sensor setup you prepared in your virtual example. 
    Simply replace numerical FRFs with experimental ones, and results are only few clicks away!

.. 
    Numerical model
    ***************
    Load the corresponding .full and .rst file from the example datasets. For more information on .full and .rst files refer to the :download:`03_FRF_synthetization.ipynb <../../../examples/03_FRF_synthetization.ipynb>` example.

    .. code-block:: python

        full_file_AB = r"./lab_testbench/FEM/AB.full"
        ress_file_AB = r"./lab_testbench/FEM/AB.rst"

        full_file_B = r"./lab_testbench/FEM/B.full"
        ress_file_B = r"./lab_testbench/FEM/B.rst"

        full_file_A = r"./lab_testbench/FEM/A.full"
        ress_file_A = r"./lab_testbench/FEM/A.rst"
        
    Create an MK model for each component:

    .. code-block:: python

        MK_A = pyFBS.MK_model(ress_file_A,full_file_A,no_modes = 100,allow_pickle= True,recalculate = False)
        MK_B = pyFBS.MK_model(ress_file_B,full_file_B,no_modes = 100,allow_pickle= True,recalculate = False)
        MK_AB = pyFBS.MK_model(ress_file_AB,full_file_AB,no_modes = 100,allow_pickle= True,recalculate = False)
        
    Update locations of channels and impacts to snap to the nearest FE node.

    .. code-block:: python

        df_chn_A_up = MK_A.update_locations_df(df_chn_A)
        df_imp_A_up = MK_A.update_locations_df(df_imp_A)

        df_chn_B_up = MK_B.update_locations_df(df_chn_B)
        df_imp_B_up = MK_B.update_locations_df(df_imp_B)

        df_chn_AB_up = MK_AB.update_locations_df(df_chn_AB)
        df_imp_AB_up = MK_AB.update_locations_df(df_imp_AB)
        
    Perform the FRF sythetization for each component based on the updated locations.

    .. code-block:: python

        MK_A.FRF_synth(df_chn_A_up,df_imp_A_up,f_start = 0,modal_damping = 0.003)
        MK_B.FRF_synth(df_chn_B_up,df_imp_B_up,f_start = 0,modal_damping = 0.003)
        MK_AB.FRF_synth(df_chn_AB_up,df_imp_AB_up,f_start = 0,modal_damping = 0.003)
    
    
Virtual point transformation
****************************

.. tip::
    It would be impractical to measure interface admittance for both substructures in multiple DoFs at the interface and furthermore ensure, 
    that these DoFs are perfectly collocated. Therefore we adopt VPT in order to obtain a collocated full-DoF interface admittance matrix for each substructure.

The VPT can be performed directly on the measured/generated FRFs. See the :download:`04_VPT.ipynb <../../../examples/04_VPT.ipynb>` example for more options and details.

.. code-block:: python

    df_vp = pd.read_excel(pos_xlsx, sheet_name='VP_Channels')
    df_vpref = pd.read_excel(pos_xlsx, sheet_name='VP_RefChannels')

    vpt_A = pyFBS.VPT(df_chn_A_up,df_imp_A_up,df_vp,df_vpref)
    vpt_B = pyFBS.VPT(df_chn_B_up,df_imp_B_up,df_vp,df_vpref)
    
Apply the defined VP transformation on the FRFs:

.. code-block:: python

    vpt_A.apply_VPT(MK_A.freq,MK_A.FRF)
    vpt_B.apply_VPT(MK_B.freq,MK_B.FRF)
    
Extract the requried FRFs and the frequency vector:

.. code-block:: python

    freq = MK_A.freq
    Y_A = vpt_A.vptData
    Y_B = vpt_B.vptData
    
LM-FBS Coupling
***************

First, construct an admittance matrix for the uncoupled system, containing substructure admittances:

.. math::
    \mathbf{Y}^\text{A|B} = \begin{bmatrix} 
    \mathbf{Y}^\text{A} & \mathbf{0} \\
    \mathbf{0} & \mathbf{Y}^\text{B}
    \end{bmatrix}

.. code-block:: python

    Y_AnB = np.zeros((Y_A.shape[0],Y_A.shape[1]+Y_B.shape[1],Y_A.shape[2]+Y_B.shape[2]), dtype=complex)

    Y_AnB[:,:Y_A.shape[1],:Y_A.shape[2]] = Y_A
    Y_AnB[:,Y_A.shape[1]:,Y_A.shape[2]:] = Y_B

Next the compatibility and the equilibrium conditions has to be defined through the signed Boolean matrices ``Bu`` and ``Bf``. 

.. math::
    \mathbf{B}_\text{u}\,\boldsymbol{u} = \mathbf{0}

.. math::
    \boldsymbol{g} = - \mathbf{B}_\text{f}^\text{T} \boldsymbol{\lambda}

Make sure that the correct DoFs are selected for the coupling. In the following example the 6 virtual/generalized DoFs at the interface are matched, 
so the size of the Boolean matrix should be 6  ×  30 (30 is the sum of all DoFs from both substructures A and B).

.. code-block:: python

    k = 6

    Bu = np.zeros((k,Y_A.shape[1]+Y_B.shape[1]))
    Bu[:k,6:6+k] = 1*np.eye(k)
    Bu[:k,12:12+k] = -1*np.eye(k)

    Bf = np.zeros((k,Y_A.shape[2]+Y_B.shape[2]))
    Bf[:k,6:6+k] = 1*np.eye(k)
    Bf[:k,12:12+k] = -1*np.eye(k)

.. figure:: ./../data/Bu.png
   :width: 500px

.. figure:: ./../data/Bf.png
   :width: 500px
    
    
For the LM FBS method, having defined :math:`\mathbf{Y^{\text{A|B}}}`, :math:`\mathbf{B}_\text{u}` and :math:`\mathbf{B}_\text{f}` is already sufficient to perform coupling:

.. math::
    \mathbf Y^{\text{AB}} = \mathbf Y^{\text{A|B}} - \mathbf Y^{\text{A|B}}\,\mathbf B^\mathrm{T} \left( \mathbf B \mathbf Y^{\mathrm{A|B}} \mathbf{B}^\mathrm{T} \right)^{-1} \mathbf B \mathbf Y^\text{A|B}

.. code-block:: python

    Y_ABn = np.zeros_like(Y_AnB,dtype = complex)

    Y_int = Bu @ Y_AnB @ Bf.T
    Y_ABn = Y_AnB - Y_AnB @ Bf.T @ np.linalg.pinv(Y_int) @ Bu @ Y_AnB
    
Results
*************

First extract the FRFs at the reference DoFs:

.. code-block:: python

    arr_coup = [0,1,2,3,4,5,18,19,20,21,22,23,24,25,26,27,28,29]
    Y_AB_coupled = Y_ABn[:,arr_coup,:][:,:,arr_coup]
    Y_AB_ref = MK_AB.FRF
    
The coupled and the reference results can then be compared and evaluated:
    
.. raw:: html

   <iframe src="../../_static/VPT_coupling.html" height="500px" width="750px" frameborder="0"></iframe>

.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!

.. rubric:: References

.. [1] de Klerk D, Rixen DJ, Voormeeren SN. General framework for dynamic substructuring: history, review and classification of techniques. AIAA journal. 2008 May;46(5):1169-81.
.. [2] van der Seijs MV, van den Bosch DD, Rixen DJ, de Klerk D. An improved methodology for the virtual point transformation of measured frequency response functions in dynamic substructuring. In4th ECCOMAS thematic conference on computational methods in structural dynamics and earthquake engineering 2013 Jun (No. 4).

    ==============================
System Equivalent Model Mixing
==============================
With System Equivalent Model Mixing (SEMM), different dynamic models of the same system can be mixed into a single hybrid model using the Lagrange-multiplier frequency-based substructuring method. 
Hybrid model follows the dynamic behavior of a precise overlay model that is expanded to the all degrees of freedom of an equivalent parent model. 
Application of SEMM comprises the expansion of the experimental dynamics to the unmeasurable DoFs for use in coupling and decoupling processes, identification of inconsistent measurements in FBS, and improving the accuracy of the experimental response models.



.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="System Equivalent Model Mixing">

.. only:: html

    .. figure:: ./../data/SEMM_result.png 
       :target: 05_SEMM.html

       System Equivalent Model Mixing

.. raw:: html

    </div>

.. toctree::
   :hidden:

   05_SEMM

|
|
|
|
|
|
|
|
|

System Equivalent Model Mixing (SEMM) was introduced by Klaassen et al [1]_ [2]_. 
The method forms a hybrid dynamic model by mixing equivalent (typically numerical and experimental) models of the same structure. 
The main goal of the method is to use the substructuring approach 
to expand the dynamics contained in a sparse overlay model :math:`\textbf{Y}^{\text{ov}}` 
onto a denser DoF space of a parent model :math:`\textbf{Y}^{\text{par}}`.

Conceptually the procedure is a two-step substructuring process, for which in the first step the experimental model 
is expanded to the denser DoF space by coupling the experimental and numerical model and then removing the excessive 
numerical contribution by introducing a decoupling of a removed model (which is a (sub)-model of the parent model).

.. figure:: ./../data/idea_of_SEMM.svg
   :width: 650px
   :align: center

Basic SEMM method
*****************

Typically, the interface in this process is considered as a set of all collocated DoFs on the parent and overlay model. 
Accordingly, the entire set of parent DoFs can be divided into internal (i) and boundary (b) subset of DoFs, 
for which the former is unique to the parent model and the latter is shared with the overlay (and removed) model.

In order to increase clarity, the same distribution will be taken into account when writing admittance matrices:

.. math::
    \mathbf{Y}^{\text{par}}\triangleq
		\begin{bmatrix}
			\mathbf{Y}_{\text{ii}}&\mathbf{Y}_{\text{ib}}\\
			\mathbf{Y}_{\text{bi}}&\mathbf{Y}_{\text{bb}}
		\end{bmatrix}^{\text{par}},\quad
		\mathbf{Y}^{\text{ov}}\triangleq
		\begin{bmatrix}
			\mathbf{Y}_{\text{bb}}
		\end{bmatrix}^{\text{ov}},\quad
		\mathbf{Y}^{\text{rem}}\triangleq
		\begin{bmatrix}
			\mathbf{Y}_{\text{bb}}
		\end{bmatrix}^{\text{par}}.

.. figure:: ./../data/num_to_par.svg
   :width: 700px
   :align: center

The LM-FBS equation is used to generate the hybrid model in the SEMM method. 
The block-diagonal equation of motion (incorporating the coupling and decoupling step) can be formulated as:

.. math::
    \begin{bmatrix}
			\boldsymbol{u}^{\text{par}}\\
			\boldsymbol{u}^{\text{rem}}\\
			\boldsymbol{u}^{\text{ov}}
		\end{bmatrix}
		=
		\begin{bmatrix}
			\mathbf{Y}^{\text{par}}& & \\
			&-\mathbf{Y}^{\text{rem}}& \\
			& &\mathbf{Y}^{\text{ov}}
		\end{bmatrix}
		\begin{bmatrix}
			\boldsymbol{f}^{\text{par}}\\
			\boldsymbol{0}\\
			\boldsymbol{0}
		\end{bmatrix}
		-
		\begin{bmatrix}
			\boldsymbol{g}^{\text{par}}\\
			\boldsymbol{g}^{\text{rem}}\\
			\boldsymbol{g}^{\text{ov}}
		\end{bmatrix}

.. figure:: ./../data/Y_matrix.svg
   :width: 250px
   :align: center

Following the LM-FBS methodology, :math:`\boldsymbol{u}` represents the displacement vector and :math:`\boldsymbol{f}` represents 
the external force vector (with the latter only acting on the parent model). 
Vector :math:`\boldsymbol{g}` represents the interface forces between the equivalent models. 
The compatibility and equilibrium conditions are defined by the following equations:

.. math::
    \mathbf{B}\boldsymbol{u}=\mathbf{0},
		\\
		\boldsymbol{g}=-\mathbf{B}^{\text{T}}\,\boldsymbol{\lambda}

where the signed Boolean matrix is defined as:

.. math::
		\mathbf{B}=
		\begin{bmatrix}
			\mathbf{B}^{\text{par}}&\mathbf{B}^{\text{rem}}&\mathbf{B}^{\text{ov}}
		\end{bmatrix}
		=
		\left[
		\begin{array}{rr|r|r}
			\mathbf{0} & -\mathbf{I} & \mathbf{I} & \mathbf{0}\\
			\mathbf{0} & \mathbf{0} & -\mathbf{I} & \mathbf{I}	
		\end{array}
		\right].

.. figure:: ./../data/B_matrix.svg
   :width: 400px
   :align: center

The substructuring procedure can be performed using the well known LM-FBS formulation:

.. math::
    \overline{\mathbf{Y}}
		=
		\mathbf{Y}
		-
		\mathbf{Y}\,\mathbf{B}^{\text{T}}\,
		\left(\mathbf{B}\, \mathbf{Y}\,\mathbf{B}^{\text{T}}\right)^{-1} \,
		\mathbf{B}\, \mathbf{Y},

where

.. math::
	\mathbf{Y}
		\triangleq
		\begin{bmatrix}
			\mathbf{Y}^{\text{par}}& & \\
			&-\mathbf{Y}^{\text{rem}}& \\
			& &\mathbf{Y}^{\text{ov}}
		\end{bmatrix}.

In order to retain the primal DoFs, the dual formulation can be reformulated using the localization matrix, 
which results in a single-line form of the basic SEMM method:

.. math::
    \mathbf{Y}^{\text{SEMM}}=
		\begin{bmatrix}
			\mathbf{Y}
		\end{bmatrix}^{\text{par}}
		-
		\begin{bmatrix}
			\mathbf{Y}_{\text{ib}}\\
			\mathbf{Y}_{\text{bb}}
		\end{bmatrix}^{\text{par}}
		%\,
		\left( \mathbf{Y}^{\text{rem}}\right) ^{-1}
		%\,
		\left( \mathbf{Y}^{\text{rem}}-\mathbf{Y}^{\text{ov}}\right)
		%\,
		\left( \mathbf{Y}^{\text{rem}}\right)^{-1}
		%\,
		\begin{bmatrix}
			\mathbf{Y}_{\text{bi}}&\mathbf{Y}_{\text{bb}}
		\end{bmatrix}^{\text{par}}.

.. tip::
	Note how the overlay admittance matrix :math:`\textbf{Y}^{\text{ov}}` is not inverted in the single-line formulation of basic SEMM? 
	In cases when overlay model is obtained through experimental testing this proves to be especially convinient, as, unlike in 
	classic FBS coupling, ill-conditioning due to experimental errors does not amplify in the inversion.

.. warning::
	When mixing equivalent models using basic SEMM, a common occurrence is the presence of spurious peaks at the internal DoFs of the hybrid model. 
	These peaks correspond neither the eigenfrequencies of the overlay nor the parent model. 


Why spurious peaks appear in the hybrid model?
==============================================

In order to be able to answer to this question, let's have a 
look at the physical model of a simple 4-DoF system. Three equivalent models are formed, but the 
dynamics of the overlay and the removed model are only obtained at the boundary DoF (thus its internal DoFs are presented transparent).

.. figure:: ./../data/formation_of_hybrid_model_basic_semm.svg
   :width: 700px
   :align: center

.. note::
	It turns out out that the hybrid model generated with SEMM is actually an 8-DoF system, even though it is collected at the 4 parent-related DoFs. 
	Thus, it is reasonable to expect the hybrid model has 8 natural frequencies and corresponding modes. 
	Compared to the overlay dynamics, this is a source of the spurious modes.

Extension of SEMM method - extended interface
*********************************************

As already established in decoupling applications, extending the interface can prove beneficial 
for the substructuring prediction. 
Therefore, Klaassen [1]_ has proposed an extension to the method, where the compatibility and equilibrium conditions 
in the decoupling step (parent - removed) are extended to the full (internal + boundary) set of DoFs.

Such a formulation has been proven to be able to eliminate the problem with spurious peaks. 
This ability is essential to improve the method's applicability. 
If the removed interface is extended to all the internal DoFs, then the removed model has the following form:

.. math::
    \mathbf{Y}^{\text{rem}}=
		\begin{bmatrix}
			\mathbf{Y}_{\text{ii}}&\mathbf{Y}_{\text{ib}}\\
			\mathbf{Y}_{\text{bi}}&\mathbf{Y}_{\text{bb}}
		\end{bmatrix}^{\text{par}}.

The equation of motion remains the same as in the basic formulation, 
only the removed model is extended to the internal DoFs:

.. math::
    \overline{\mathbf{Y}}
		=
		\mathbf{Y}
		-
		\mathbf{Y}\,\mathbf{B}^{\text{T}}\,
		\left(\mathbf{B}\, \mathbf{Y}\,\mathbf{B}^{\text{T}}\right)^{-1} \,
		\mathbf{B}\, \mathbf{Y},
		\text{ where: }
		\mathbf{Y}
		\triangleq
		\begin{bmatrix}
			\mathbf{Y}^{\text{par}}& & \\
			&-\mathbf{Y}^{\text{rem}}& \\
			& &\mathbf{Y}^{\text{ov}}
		\end{bmatrix}.

The final version of the fully extended SEMM method can also be written in a single-line notation 
where relations defined in matrix :math:`\textbf{B}` are directly applied:

.. math::
    \mathbf{Y}^{\text{SEMM}}=
		\mathbf{Y}^{\text{par}}
		-
		\mathbf{Y}^{\text{par}}
		%\,
		\left( \begin{bmatrix} \mathbf{Y}_{\text{bi}}&\mathbf{Y}_{\text{bb}}\end{bmatrix}^{\text{rem}}\right ) ^{+}
		%\,
		\left( \mathbf{Y}^{\text{rem}}_{\text{bb}}-\mathbf{Y}^{\text{ov}}\right)
		%\,
		\left( \begin{bmatrix} \mathbf{Y}_{\text{ib}}\\ \mathbf{Y}_{\text{bb}}\end{bmatrix}^{\text{rem}}\right ) ^{+}
		%\,
		\mathbf{Y}^{\text{par}}.

.. note::
    In the case of a fully extended formulation of the SEMM method, special attention must be paid to the manipulation 
    of the interface matrix that is inverted during decoupling. 
    It turns out that this matrix is always singular, so the inverse does not exist. 
    The problem is solved by using a pseudoinverse instead which enables the inversion of the matrix and hence the calculation 
    of the response of the hybrid model.

.. warning::
    Due to the use of pseudoinverse, the connections between equivalent models are no longer perfect. Therefore, in the hybrid model, 
    duplicated DoFs do not have the same response. It turns out that the only real response that reflects the 
    experimental model most accurately is the response at the DoFs of the parent model.	This is already considered in the single line notation of the fully extended SEMM.

.. tip::
	SEMM has been successfully employed in a variety of applications, such as dynamic coupling
	within frequency-based substructuring [3]_ and substructuring-based joint identification [4]_. 
	Another possible application of the expansion method is the detection of inconsistent measurements,
	which was implemented in [5]_ [6]_.
	Although the framework is typically considered for experimental-numerical hybrid modeling, the same procedure can be
	applied for purely experimental applications, such as improving full-field high-speed [7]_ or acoustic [8]_ camera
	measurements. 

.. rubric:: References

.. [1] Klaassen SW, van der Seijs MV, de Klerk D. System equivalent model mixing. Mechanical Systems and Signal Processing. 2018 May 15;105:90-112.
.. [2] Klaassen SW, van der Seijs MV. Introducing semm: A novel method for hybrid modelling. InDynamics of Coupled Structures, Volume 4 2018 (pp. 117-125). Springer, Cham.
.. [3] Pasma E, Klaassen S, Nieuwenhuijse L, Van Der Seijs M, Lennström D. Application of system equivalent model mixing (SEMM) to model the structural dynamic properties of a complex vehicle component using numerical and experimental data. Proceedings of ISMA2018. 2018.
.. [4] Klaassen SW, Rixen DJ. Using SEMM to Identify the Joint Dynamics in Multiple Degrees of Freedom Without Measuring Interfaces. InDynamic Substructures, Volume 4 2020 (pp. 87-99). Springer, Cham.
.. [5] Kodrič M, Čepon G, Boltežar M. Experimental framework for identifying inconsistent measurements in frequency-based substructuring. Mechanical Systems and Signal Processing. 2021 Jun 1;154:107562.
.. [6] Saeed Z, Firrone CM, Berruti TM. Joint identification through hybrid models improved by correlations. Journal of Sound and Vibration. 2021 Mar 3;494:115889.
.. [7] Bregar T, Zaletelj K, Čepon G, Slavič J, Boltežar M. Full-field FRF estimation from noisy high-speed-camera data using a dynamic substructuring approach. Mechanical Systems and Signal Processing. 2021 Mar 1;150:107263.
.. [8] Ocepek D, Kodrič M, Čepon G, Boltežar M. On the estimation of structural admittances from acoustic measurement using a dynamic substructuring approach. Applied Acoustics. 2021 Sep 1;180:108115.##############
VPT Decoupling
##############

The decoupling of susbtructures is performed in a similar fashion as the coupling, with the only difference that a minus sign 
must be applied on the subsystem to be decoupled. The operation is also based on the Lagrange-Multiplier Frequency-Based Substructuring (LM-FBS) [1]_
formulation. In the following, a basic decoupling of two numerically-generated substructures is presented. 
The virtual point transformation [2]_ is applied to impose collocated matching DoFs at the interface. 
This can also be performed analogously with experimentally acquired data.

.. note:: 
   Download example showing a substructure decoupling application: :download:`08_decoupling_VPT.ipynb <../../../examples/08_FBS_decoupling_VPT.ipynb>`
    
.. tip::
    Why use virtual point when decoupling substructures?

    * Complex interfaces are very often inaccessible for the measurement equipment. Therefore measurements are often performed away from the actual interface. Also, the collocation of DoFs at the neighboring interfaces is almost impossible to ensure. Virtual point, on the other hand, can be defined in an arbitrary location at the interface that coincides for all substructures.
    * Virtual point accounts for rotational degrees of freedom, which are mandatory for the successful decoupling of the substructures. 
    * By the reduction of measurements to the virtual point the interface problem is weakened. That means that compatibility and equilibrium conditions on the measured DoFs are a bit more relaxed. In this manner, unwanted stiffening effects can be avoided. 
    * Due to the reduction, measurement errors such as bias from sensor positioning or uncorrelated measurement noise are filtered out to some extend.
 

Example Datasets and 3D view
****************************

Load the required predefined datasets and open the 3D viewer in the background as already shown in `3D Display <../../../html/examples/basic_examples/01_static_display.html>`_. Also for decoupling, a subplot representation, as already presented in `Coupling <../../../html/examples/fbs/07_coupling.html>`_, can be used.
    
.. figure:: ./../data/eight_three.png
   :width: 500px
   
.. tip::
    With ``pyFBS`` you can simply prepare the experiments before hand! 
    Position your virtual accelerometers and sensors on 3D model, generate numerical FRFs and make sure that everything is in order. 
    Then, perform your experiment, following the sensor setup you prepared in your virtual example. 
    Simply replace numerical FRFs with experimental ones, and results are only few clicks away!
..    
    Numerical model
    ***************
    Load the corresponding .full and .ress file from the example datasets. For more information on .full and .ress files refer to the :download:`03_FRF_synthetization.ipynb <../../../examples/03_FRF_synthetization.ipynb>` example

    .. code-block:: python

        full_file_AB = r"./lab_testbench/FEM/AB.full"
        ress_file_AB = r"./lab_testbench/FEM/AB.rst"

        full_file_B = r"./lab_testbench/FEM/B.full"
        ress_file_B = r"./lab_testbench/FEM/B.rst"

        full_file_A = r"./lab_testbench/FEM/A.full"
        ress_file_A = r"./lab_testbench/FEM/A.rst"
        
    Create an MK model for each component:

    .. code-block:: python

        MK_A = pyFBS.MK_model(ress_file_A,full_file_A,no_modes = 100,allow_pickle= True,recalculate = False)
        MK_B = pyFBS.MK_model(ress_file_B,full_file_B,no_modes = 100,allow_pickle= True,recalculate = False)
        MK_AB = pyFBS.MK_model(ress_file_AB,full_file_AB,no_modes = 100,allow_pickle= True,recalculate = False)
        
    Update locations of channels and impacts to snap to the nearest FE node.

    .. code-block:: python

        df_chn_A_up = MK_A.update_locations_df(df_chn_A)
        df_imp_A_up = MK_A.update_locations_df(df_imp_A)

        df_chn_B_up = MK_B.update_locations_df(df_chn_B)
        df_imp_B_up = MK_B.update_locations_df(df_imp_B)

        df_chn_AB_up = MK_AB.update_locations_df(df_chn_AB)
        df_imp_AB_up = MK_AB.update_locations_df(df_imp_AB)
        
    Perform the FRF sythetization for each component based on the updated locations.

    .. code-block:: python

        MK_A.FRF_synth(df_chn_A_up,df_imp_A_up,f_start = 0,modal_damping = 0.003)
        MK_B.FRF_synth(df_chn_B_up,df_imp_B_up,f_start = 0,modal_damping = 0.003)
        MK_AB.FRF_synth(df_chn_AB_up,df_imp_AB_up,f_start = 0,modal_damping = 0.003)
    
Virtual point transformation
****************************
.. tip::
    It would be impractical to measure interface admittance for both substructures in multiple DoFs at the interface and furthermore ensure, 
    that these DoFs are perfectly collocated. Therefore we adopt VPT in order to obtain a collocated full-DoF interface admittance matrix for each substructure.


The VPT can be performed directly on the generated data. See the :download:`04_VPT.ipynb <../../../examples/04_VPT.ipynb>` example for more options and details.

.. code-block:: python

    df_vp = pd.read_excel(pos_xlsx, sheet_name='VP_Channels')
    df_vpref = pd.read_excel(pos_xlsx, sheet_name='VP_RefChannels')

    vpt_AB = pyFBS.VPT(df_chn_AB_up,df_imp_AB_up,df_vp,df_vpref)
    vpt_B = pyFBS.VPT(df_chn_B_up,df_imp_B_up,df_vp,df_vpref)
    
Apply the defined VP transformation on the FRFs:

.. code-block:: python

    vpt_AB.apply_VPT(MK_AB.freq,MK_AB.FRF)
    vpt_B.apply_VPT(MK_B.freq,MK_B.FRF)
    
Extract the requried FRFs and the frequency vector:

.. code-block:: python

    freq = MK_AB.freq
    Y_AB = vpt_AB.vptData
    Y_B = vpt_B.vptData
    
LM-FBS Decoupling
*****************
First, construct an admittance matrix for the uncoupled system, containing substructure admittances. 
Note that the operation is equivalent to the one performed for the coupling case with the difference that a minus sign here is applied on the subsystem to be decoupled.

.. math::
    \mathbf{Y}^\text{AB|B} = \begin{bmatrix} 
    \mathbf{Y}^\text{AB} & \mathbf{0} \\
    \mathbf{0} & -\mathbf{Y}^\text{B}
    \end{bmatrix}

.. code-block:: python

    Y_ABnB = np.zeros((Y_AB.shape[0],Y_AB.shape[1]+Y_B.shape[1],Y_AB.shape[2]+Y_B.shape[2]), dtype=complex)

    Y_ABnB[:,:Y_AB.shape[1],:Y_AB.shape[2]] = Y_AB
    Y_ABnB[:,Y_AB.shape[1]:,Y_AB.shape[2]:] = -1*Y_B


Next the compatibility and the equilibrium conditions has to be defined through the signed Boolean matrices ``Bu`` and ``Bf``. 

.. math::
    \mathbf{B}_\text{u}\,\boldsymbol{u} = \mathbf{0}

.. math::
    \boldsymbol{g} = - \mathbf{B}_\text{f}^\text{T} \boldsymbol{\lambda}

Make sure that the correct DoFs are selected for the decoupling. For this case, the interface is extended to the internal DoFs common to both AB and B, 
making in total 6 + 12 compatibility and equilibrium conditions. 
Adding compatibility and equilibrium at the internal DoFs contributes to increase the observability and controllability 
of the interface dynamics [3]_.

.. code-block:: python

    k = 6 + 12 

    Bu = np.zeros((k,Y_AB.shape[1]+Y_B.shape[1]))
    Bu[:k,6:6+k] = 1*np.eye(k)
    Bu[:k,24:24+k] = -1*np.eye(k)

    Bf = np.zeros((k,Y_AB.shape[1]+Y_B.shape[1]))
    Bf[:k,6:6+k] = 1*np.eye(k)
    Bf[:k,24:24+k] = -1*np.eye(k)
    
.. figure:: ./../data/Bu_decoupling.png
   :width: 600px

.. figure:: ./../data/Bf_decoupling.png
   :width: 600px
    
For the LM FBS method, having defined :math:`\mathbf{Y^{\text{AB|B}}}`, :math:`\mathbf{B}_\text{u}` and :math:`\mathbf{B}_\text{f}` is already sufficient to perform coupling:

.. math::
    \mathbf Y^{\text{A}} = \mathbf Y^{\text{AB|B}} - \mathbf Y^{\text{AB|B}}\,\mathbf B^\mathrm{T} \left( \mathbf B \mathbf Y^{\mathrm{AB|B}} \mathbf{B}^\mathrm{T} \right)^{-1} \mathbf B \mathbf Y^\text{AB|B}

.. code-block:: python

    Y_An = np.zeros_like(Y_ABnB,dtype = complex)

    Y_int = Bu @ Y_ABnB @ Bf.T
    Y_An =Y_ABnB - Y_ABnB @ Bf.T @ np.linalg.pinv(Y_int) @ Bu @ Y_ABnB
    
Results
*************
First extract the FRFs at the reference DoFs:

.. code-block:: python

    arr_coup = [0,1,2,3,4,5]
    Y_A_coupled = Y_An[:,arr_coup,:][:,:,arr_coup]
    Y_A_ref = MK_A.FRF
    
The decoupled and the reference results can then be compared and evaluated:
   
.. raw:: html

   <iframe src="../../_static/VP_decoupling.html" height="500px" width="750px" frameborder="0"></iframe>
   
.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!

.. rubric:: References

.. [1] de Klerk D, Rixen DJ, Voormeeren SN. General framework for dynamic substructuring: history, review and classification of techniques. AIAA journal. 2008 May;46(5):1169-81.
.. [2] van der Seijs MV, van den Bosch DD, Rixen DJ, de Klerk D. An improved methodology for the virtual point transformation of measured frequency response functions in dynamic substructuring. In4th ECCOMAS thematic conference on computational methods in structural dynamics and earthquake engineering 2013 Jun (No. 4).
.. [3] Voormeeren SN, Rixen DJ. A family of substructure decoupling techniques based on a dual assembly approach. Mechanical Systems and Signal Processing. 2012 Feb 1;27:379-96.##############
SVT Decoupling
##############

pyFBS has implemented the novel SVD-based approach for interface reduction in LM-FBS. The SVT is the first engineering tool to tackle the issue of flexible interfaces in lightly damped system for experimental frequency-based substructuring. The methodology can be applied without any knowledge of system geometry and treat efficiently measurement error by combining reduction, filtering and regularization in a single transformation step.

.. note:: 
   Download example showing a substructure decoupling application with SVT: :download:`19_FBS_decoupling_SVT.ipynb <../../../examples/19_FBS_decoupling_SVT.ipynb>`
    
Example Datasets and 3D view
****************************

Load the required predefined datasets and open the 3D viewer in the background as already shown in `3D Display <../../../html/examples/basic_examples/01_static_display.html>`_. Also for decoupling, a subplot representation, as already presented in `Coupling <../../../html/examples/fbs/07_coupling.html>`_, can be used.
    
.. figure:: ./../data/eight_three.png
   :width: 500px
   
    
Experimental model
******************
Load the experimental measurements file from the example datasets.

.. code-block:: python

	exp_A = r"./lab_testbench/Measurements/Y_A.p"
	exp_B = r"./lab_testbench/Measurements/Y_B.p"
	exp_AB = r"./lab_testbench/Measurements/Y_AB.p"

	freq, _Y_A_exp = np.load(exp_A, allow_pickle = True)
	_, _Y_B_exp = np.load(exp_B, allow_pickle = True)
	_, _Y_AB_exp = np.load(exp_AB, allow_pickle = True)

	Y_A_exp = np.transpose(_Y_A_exp, (2, 0, 1))
	Y_B_exp = np.transpose(_Y_B_exp, (2, 0, 1))
	Y_AB_exp = np.transpose(_Y_AB_exp, (2, 0, 1))
	

Singular vector transformation
******************************
The SVT can be performed directly on the measured data. Make sure that input and output DoFs involved in the transformation are the same in the subsystem to be decoupled (B) and the assembled system (AB). Furthermore, the same reduction spaces must be used for both the systems (B and AB) in order to guarantee a proper compatibility and equilibrium. The reduced singular subspaces are extracted from the subsystem B and a number of 6 DoFs is used.

.. code-block:: python

	k = 6
	svt = pyFBS.SVT(df_chn_B,df_imp_B,freq,Y_B_exp,[1,10],k)
	
    
Apply the defined SVT to systems B and AB:

.. code-block:: python

	_,_,FRF_B_sv= svt.apply_SVT(df_chn_B,df_imp_B,freq,Y_B_exp)
	_,_,FRF_AB_sv= svt.apply_SVT(df_chn_AB,df_imp_AB,freq,Y_AB_exp)


LM-FBS Decoupling
*****************
The uncoupled global admittance is constructed using the transformed FRF datasets and imposing the decoupling operation 
through the minus sign. 

.. code-block:: python

  	Y_AB_un = np.zeros((len(freq),2*k+6,2*k+6),dtype = complex)

	Y_AB_un[:,0:2*k,0:2*k] = FRF_AB_sv
	Y_AB_un[:,2*k:,2*k:] = -1*FRF_B_sv

	plt.spy(np.abs(Y_AB_un[100])) # display at arbitrary frequency to check for shape

.. figure:: ./../data/Y_SVT.png
   :width: 300px

The compatibility and the equilibrium conditions are defined via the signed Boolean matrices.

.. code-block:: python

	Bu = np.zeros((k,2*k+6))
	Bu[:k,0:k] = 1*np.eye(k)
	Bu[:k,2*k:2*k+6] = -1*np.eye(k)

	Bf = Bu
    
.. figure:: ./../data/Bu_SVT.png
   :width: 400px
.. figure:: ./../data/Bf_SVT.png
   :width: 400px
   
Apply the LM-FBS based on the defined compatibility and equilibrium conditions.

.. code-block:: python

	Y_int = Bu @ Y_AB_un @ Bf.T
	Y_A_dec  = Y_AB_un - Y_AB_un @ Bf.T @ np.linalg.pinv(Y_int) @ Bu @ Y_AB_un

Results
*******

First extract the FRFs at the reference DoFs:

.. code-block:: python

	arr_ = [6,7,8,9,10,11]
	Y_A_LMFBS = Y_A_dec[:,arr_,:][:,:,arr_]
	Y_A_ref = Y_A_exp
    
The decoupled and the reference results for A can be compared:
   
.. raw:: html

   <iframe src="../../_static/decoupling_SVT.html" height="500px" width="750px" frameborder="0"></iframe>

.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!System Equivalent Model Mixing
==============================

System Equivalent Model Mixing (SEMM) [1]_ enables the mixing of equivalent models into a hybrid model in the frequency domain. 
The models used can either be of numerical or experimental nature. 
The overlay model provides the dynamic properties which are expanded to the DoFs of the parent model. 
Therefore the overlay model is usually represented by the experimental model and parent model with the numerical model.

.. note:: 
   Download example showing the basic use of SEMM: :download:`05_SEMM.ipynb <../../../examples/05_SEMM.ipynb>`

..
   DoF-set of parent model is contained from internal (i) and boundary (b) DoFs. 
   Boundary DoFs must overlap with the overlay model so the dynamic coupling can be performed, while the internal DoFs of the parent model can be unique to its own. 
   The equivalent models, appearing in the SEMM method, are arranged by separating internal and boundary DoFs in the admittance matrices:

   .. math::

      \begin{equation}\label{parent_overlay}
      \mathbf{Y}^{\text{par}}=
      \begin{bmatrix}
      \mathbf{Y}_{\text{ii}}&\mathbf{Y}_{\text{ib}}\\
      \mathbf{Y}_{\text{bi}}&\mathbf{Y}_{\text{bb}}
      \end{bmatrix}^{\text{par}},\quad
      \mathbf{Y}^{\text{ov}}=
      \begin{bmatrix}
      \mathbf{Y}_{\text{bb}}
      \end{bmatrix}^{\text{ov}},\quad
      \mathbf{Y}^{\text{rem}}=
      \begin{bmatrix}
      \mathbf{Y}_{\text{bb}}
      \end{bmatrix}^{\text{rem}}.
      \end{equation}

   After satisfying compatibility and equilibrium conditions between equivalent models, the basic form of the SEMM method is defined using the equation:

   .. math::

      \mathbf{Y}^{\text{SEMM}}=
      \begin{bmatrix}
      \mathbf{Y}
      \end{bmatrix}^{\text{par}}
      -
      \begin{bmatrix}
      \mathbf{Y}_{\text{ib}}\\
      \mathbf{Y}_{\text{bb}}
      \end{bmatrix}^{\text{par}}
      %\,
      \left( \mathbf{Y}^{\text{rem}}\right) ^{-1}
      %\,
      \left( \mathbf{Y}^{\text{rem}}-\mathbf{Y}^{\text{ov}}\right)
      %\,
      \left( \mathbf{Y}^{\text{rem}}\right)^{-1}
      %\,
      \begin{bmatrix}
      \mathbf{Y}_{\text{bi}}&\mathbf{Y}_{\text{bb}}
      \end{bmatrix}^{\text{par}}

   By extending the removed model to all DoFs of numerical model, the fully extend formulation of SEMM method follows equation:

   .. math::

      \mathbf{Y}^{\text{SEMM}}=
      \mathbf{Y}^{\text{par}}
      -
      \mathbf{Y}^{\text{par}}
      %\,
      \left( \begin{bmatrix} \mathbf{Y}_{\text{bi}}&\mathbf{Y}_{\text{bb}}\end{bmatrix}^{\text{rem}}\right ) ^{+}
      %\,
      \left( \mathbf{Y}^{\text{rem}}_{\text{bb}}-\mathbf{Y}^{\text{ov}}\right)
      %\,
      \left( \begin{bmatrix} \mathbf{Y}_{\text{ib}}\\ \mathbf{Y}_{\text{bb}}\end{bmatrix}^{\text{rem}}\right ) ^{+}
      %\,
      \mathbf{Y}^{\text{par}}

Example data import
*******************

In the beginning, it necessary to define numerical and experimental data (response models). Datasets used in this example are from a laboratory testbench and are available directly within the :mod:`pyFBS`.

Experimental model
------------------
The experimental model must be properly aranged so that the FRFs are of correct shape. 
The first dimension represents the frequency depth, second the response points, and the third excitation points. 
The experimental model is used as the overlay model.

.. code-block:: python

   exp_file = r"./lab_testbench/Measurements/Y_AB.p"

   freq, Y_exp = np.load(exp_file, allow_pickle = True)
   Y_exp = np.transpose(Y_exp, (2, 0, 1))


Numerical model
---------------
FRFs used for the parent model numerical model can be imported to Python or generated with the mass and stiffness matrix using :mod:`pyFBS.MK_model.FRF_synth`.
Locations and directions for which FRFs are generated are defined in an Excel file and can be parsed into a :mod:`pandas.DataFrame`.

.. code-block:: python

   stl = r"./lab_testbench/STL/AB.stl"
   xlsx = r"./lab_testbench/Measurements/AM_measurements.xlsx"

   full_file = r"./lab_testbench/FEM/AB.full"
   rst_file = r"./lab_testbench/FEM/AB.rst"

   MK = pyFBS.MK_model(rst_file, full_file, no_modes = 100, recalculate = False)

   df_chn = pd.read_excel(xlsx, sheet_name='Channels_AB')
   df_imp = pd.read_excel(xlsx, sheet_name='Impacts_AB')

   MK.FRF_synth(df_chn,df_imp, 
                f_start=0,
                f_end=2002.5,
                f_resolution=2.5,
                modal_damping = 0.003,
                frf_type = "accelerance")

As experimental, also the numerical model must be properly arranged. The first dimension represents the frequency depth, the second the response points, and the third excitation points. 
Additionally, the frequency resolution of the numerical and experimental model has to be the same.

The numerical model is not necessarily a square matrix, as it can also be rectangular. 
But the numerical model must contain at least all the DoFs contained in the experimental model.

Application of SEMM
*******************

The function enables the implementation of three SEMM method formulations: ``basic``, ``fully-extend`` and ``fully-extend-svd``, the choice of which is defined with the ``SEMM_type`` parameter.
The ``red_comp`` and ``red_eq`` parameters can be used to influence the number of eigenvalues used to ensure equilibrium and compatibility conditions when the ``fully-extend-svd`` formulation is used [2]_.

The result is a hybrid model that contains the DoFs represented in the numerical model.

In the example below, only a part of the experimental response matrix is included in SEMM. The remaining DoFs are used to evaluate SEMM.
It is essential that the order of measurements in the experimental model ``Y_exp`` coincides with the order of measurements in the parameters ``df_chn_exp`` and ``df_imp_exp`` and the same must be valid for the numerical model. 
Function :mod:`pyFBS.SEMM` will automatically match corresponding DoFs.

.. code-block:: python

   Y_AB_SEMM = pyFBS.SEMM(MK.FRF, Y_exp[:, 0:15, 5:20],
                          df_chn_num = df_chn, 
                          df_imp_num = df_imp, 
                          df_chn_exp = df_chn[0:15], 
                          df_imp_exp = df_imp[5:20], 
                          SEMM_type='fully-extend-svd', red_comp=10, red_eq=10)

Finally, the results of the hybrid model can be compared with the reference experimental and the numerical model.

..
   .. code-block:: python

      s1 = 24
      s2 = 24

      display(df_chn.iloc[[s1]])
      display(df_imp.iloc[[s2]])

      plt.figure(figsize = (12,8))

      plt.subplot(211)
      plt.semilogy(MK.freq,np.abs(MK.FRF[:,s1,s2]), label = "Num.")
      plt.semilogy(freq,np.abs(Y_exp[:, s1,s2]), label = "Exp.")
      plt.semilogy(freq,np.abs(Y_AB_SEMM[:, s1,s2]), label = "SEMM")
      plt.ylabel("Accelerance [m/s$^2$/N]")
      plt.legend()

      plt.subplot(413)
      plt.plot(MK.freq,np.angle(MK.FRF[:,s1,s2]))
      plt.plot(freq,np.angle(Y_exp[:, s1,s2]))
      plt.plot(MK.freq,np.angle(Y_AB_SEMM[:,s1,s2]))
      plt.xlabel("f [Hz]")
      plt.ylabel("Angle [rad]")
   
.. raw:: html

   <iframe src="../../_static/SEMM_plot.html" height="500px" width="750px" frameborder="0"></iframe>

.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!
   
.. rubric:: References

.. [1] Steven WB Klaassen, Maarten V. van der Seijs, and Dennis de Klerk. System equivalent model mixing. Mechanical Systems and Signal Processing, 105:90–112, 2018.
.. [2] Steven WB Klaassen and D. J. Rixen. The Inclusion of a Singular-value Based Filter in SEMM. in: Proceedings of the 38th International Modal Analysis Conference, A Conference on Structural Dynamics, (2020), 2020.
==============================
Virtual Point Transformation
==============================
So... You want to couple two substructures and you already know you need to measure frequency response functions for all translational and rotational 
degrees of freedom directly at the interface. But you can't even reach the interface with standard tri-axial accelerometers due to the geometry of the interface? 
And what about the rotational degrees of freedom? 

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Virtual Point Transformation">

.. only:: html

    .. figure:: ./../data/pic_vpt.png 
       :target: 04_VPT.html

       Virtual Point Transformation

.. raw:: html

    </div>

.. toctree::
   :hidden:

   04_VPT




.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Coupling">

.. only:: html

    .. figure:: ./../data/seven_three.png   
       :target: 07_coupling.html

       FBS Coupling with VPT

.. raw:: html

    </div>

.. toctree::
   :hidden:

   07_coupling


  
 
.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Decoupling">

.. only:: html

    .. figure:: ./../data/seven_three.png   
       :target: 08_decoupling.html

       FBS Decoupling with VPT

.. raw:: html

    </div>

.. toctree::
   :hidden:

   08_decoupling

|
|
|
|
|
|
|
|
|

In general, the main challenges when coupling two substructures can be summarized with the following points: 

* Non-collocated degrees of freedom on neighboring substructures.
* Lack of rotational degrees of freedom due to limitations in measurement equipment and excitation capabilities.
* Random and systematic measurement errors.

Above-mentioned difficulties can be elegantly overcome using virtual point transformation (VPT). Let us consider a simple interface where the response 
is captured using tri-axial accelerometers and excited with unidirectional forces around the interface as depicted below. The virtual point is placed 
in the middle of the hole.

.. figure:: ./../data/vp_a_b_ab.png
   :width: 700px

Virtual point displacements
***************************

If we observe only the interface, the modes that dominantly represent its response are the translational motions in :math:`x`, :math:`y` and :math:`z` axis, 
as well as rotations about :math:`x`, :math:`y` and :math:`z` axis [1]_. These are named interface deformation modes (IDMs) [2]_. 

.. |IDM1| image:: ./../data/IDM_1.gif
    :width: 120px
.. |IDM2| image:: ./../data/IDM_2.gif
    :width: 120px
.. |IDM3| image:: ./../data/IDM_3.gif
    :width: 120px
.. |IDM4| image:: ./../data/IDM_4.gif
    :width: 120px
.. |IDM5| image:: ./../data/IDM_5.gif
    :width: 120px
.. |IDM6| image:: ./../data/IDM_6.gif
    :width: 120px

.. table::         

  +------+------+------+------+------+------+
  ||IDM1|||IDM2|||IDM3|||IDM4|||IDM5|||IDM6|| 
  +------+------+------+------+------+------+


In the following, we assume rigid IDMs 
are dominant while flexible IDMs can be neglected. 
If we would select one point at the interface (virtual point), an arbitrary channel at :math:`i`-th sensor is then rigidly connected to that point.

.. figure:: ./../data/vp_acc.svg
   :width: 200px

If we treat interface as perfectly rigid, the VP has six rigid displacements :math:`\boldsymbol{q}=[q_X,q_Y,q_Z,q_{\theta_X},q_{\theta_Y},q_{\theta_Z}]^\text{T}`. 
Displacement :math:`u_x^i` can be expressed from :math:`\boldsymbol{q}`; one musk ask himself, which movements of the VP contribute to the :math:`u_x^i` response 
(in the VP coordinate system :math:`XYZ`):

.. figure:: ./../data/disp.svg
   :width: 600px

.. math::
    u_x^i = \underbrace{1 \cdot q_X}_\mathrm{contribution\,from\,\textit{X}\,translation} + 0 \cdot q_Y + 0 \cdot q_Z + 0 \cdot q_{\theta_X} + 
    \underbrace{r_Z \cdot q_{\theta_Y}}_\mathrm{contribution\,from\,\textit{Y}\,rotation} - 
    \underbrace{r_Y \cdot q_{\theta_Z}}_\mathrm{contribution\,from\,\textit{Z}\,rotation}

Similarly, responses :math:`u_y^i` and :math:`u_z^i` are obtained from :math:`\boldsymbol{q}` 
and the relation between one tri-axial sensor and a virtual point can be written as follows:

.. math::
    \begin{bmatrix}
    u_x^i \\ u_y^i \\ u_z ^i
    \end{bmatrix}
    = 
    \begin{bmatrix}
    1 & 0 & 0 & 0 & r_Z & -r_Y \\
    0 & 1 & 0 & -r_Z & 0 & r_X \\
    \underbrace{0}_{\text{translation in }X} & \underbrace{0}_{\text{translation in }Y} & \underbrace{1}_{\text{translation in }Z} & \underbrace{r_Y}_{\text{rotation around }x} & \underbrace{-r_X}_{\text{rotation around }y} & \underbrace{0}_{\text{rotation around }Z} 
    \end{bmatrix}
    \begin{bmatrix}
    q_X \\ q_Y \\ q_Z \\ q_{\theta_X} \\ q_{\theta_Y} \\ q_{\theta_Z}
    \end{bmatrix}

Columns of the :math:`\mathbf{R}` matrix are rigid IDMs, which are assembled from relative sensor position with regard to the VP.
For cases where the orientation of sensor channels and VP mismatch, IDMs are transformed in the direction of the sensor channels [2]_:

.. math::
    \begin{bmatrix}
    u_x^i \\ u_y^i \\ u_z ^i
    \end{bmatrix}
    =
    \begin{bmatrix}
    e_{x,\,X} & e_{x,\,Y} & e_{x,\,Z} \\
    e_{y,\,X} & e_{y,\,Y} & e_{y,\,Z} \\
    e_{z,\,X} & e_{z,\,Y} & e_{z,\,Z}
    \end{bmatrix}
    \begin{bmatrix}
    1 & 0 & 0 & 0 & r_Z & -r_Y \\
    0 & 1 & 0 & -r_Z & 0 & r_X \\
    0 & 0 & 1 & r_Y & -r_X & 0 
    \end{bmatrix}\, \boldsymbol{q}
    \quad \Rightarrow \quad \boldsymbol{u}^i = \mathbf{R}_i \, \boldsymbol{q}

where :math:`[e_{x,X}, e_{x,Y}, e_{x,Z}]^\text{T}` is orientation of :math:`x` sensor channel, :math:`[e_{y,X}, e_{y,Y}, e_{y,Z}]^\text{T}` is orientation of :math:`y`
sensor channel and :math:`[e_{z,X}, e_{z,Y}, e_{z,Z}]^\text{T}` is orientation of :math:`z` sensor channel, all in :math:`XYZ` coordinate system.
Two tri-axial accelerometers are not sufficient to properly describe :math:`\boldsymbol{q}` as it is not possible to describe all three rotations with 
only two tri-axial sensors only. Therefore, it is suggested to use three or more tri-axial accelerometers and the reduction of measurements is performed [2]_. 
The use of laser vibrometers [3]_ or rotational accelerometers in VPT is also encouraged [4]_ as both yield promising results [5]_. 
All relations between measured sensors displacements :math:`\boldsymbol{u}` and :math:`\boldsymbol{q}` can be assembled into:

.. math::
    \begin{bmatrix}
    \boldsymbol{u}^1 \\ \boldsymbol{u}^2 \\ \boldsymbol{u}^3 \\ \vdots
    \end{bmatrix}
    =
    \begin{bmatrix}
    \mathbf{R}_1 \\
    \mathbf{R}_2 \\
    \mathbf{R}_3 \\
    \vdots \\
    \end{bmatrix}
    \, \boldsymbol{q}
    \quad \Rightarrow \quad \boldsymbol{u} = \mathbf{R}\,\boldsymbol{q}

Since the interface rigidity assumption is not always fully satisfied, a residual term :math:`\boldsymbol{\mu}` is added 
which represents the flexible motion not comprised within rigid IDMs:

.. math::
    \boldsymbol{u} = \mathbf{R}\,\boldsymbol{q} + \boldsymbol{\mu}

To find the solution :math:`\boldsymbol{q}` that best approximates the measured response :math:`\boldsymbol{u}`, a residual cost function 
:math:`\boldsymbol{\mu}^\text{T}\boldsymbol{\mu}` has to be minimized [6]_. A solution is found in a least-square sense:

.. math::
    \boldsymbol{q} = \Big( \mathbf{R}^{\text{T}} \mathbf{R} \Big)^{-1} \mathbf{R}^{\text{T}} \, \boldsymbol{u} = \mathbf{T}_{\text{u}} \,\boldsymbol{u} \quad \Rightarrow 
    \quad \mathbf{T}_{\text{u}} = \Big( \mathbf{R}^{\text{T}} \mathbf{R} \Big)^{-1} \mathbf{R}^{\text{T}}

To gain more flexibility over the transformation, symmetric weighting matrix :math:`\mathbf{W}` is introduced [2]_. 
If necessary, one could add more weight or even exclude specific displacements from the VPT. The :math:`\mathbf{W}` can also be defined per frequency line, 
so various displamenets can be managed across the entire frequency range.

.. math::
    \boldsymbol{q} = \Big( \mathbf{R}^{\text{T}} \mathbf{W} \mathbf{R} \Big)^{-1} \mathbf{R}^{\text{T}} \mathbf{W} \, \boldsymbol{u} = \mathbf{T}_{\text{u}} \,
    \boldsymbol{u} \quad \Rightarrow \quad \mathbf{T}_{\text{u}} = \Big( \mathbf{R}^{\text{T}} \mathbf{W} \mathbf{R} \Big)^{-1} \mathbf{R}^{\text{T}} \mathbf{W}

The matrix :math:`\mathbf{T}_u` is a displacement transformation matrix that projects 
translational displacements into a subspace composed by six rigid IDMs, retaining only the dynamics that load the interface 
in a purely rigid manner. 
The residual flexibility is filtered out from :math:`\boldsymbol{q}`. 
For cases where the interface exhibits non-neglectable flexible behaviour, additional DoFs can be added to the VP to describe this flexible motion [7]_.

Virtual point forces
********************

Similar procedure is applied to transform interaface loads onto the VP. Virtual forces and moments 
:math:`\boldsymbol{m}=[m_X,m_Y,m_Z,m_{\theta_X},m_{\theta_Y},m_{\theta_Z}]^\text{T}`
as a consequence of a single unidirectional (:math:`i`-th) impact are expressed as (note that forces :math:`\boldsymbol{f}` are not uniquely defined by 
:math:`\boldsymbol{m}` as there are infinitely many possible solutions for :math:`\boldsymbol{f}`) [2]_:

.. math::
    \begin{bmatrix} m_X \\ m_Y \\ m_Z \\ m_{\theta_X} \\ m_{\theta_Y} \\ m_{\theta_Z} \end{bmatrix} =
    \begin{bmatrix}
        1 & 0 & 0 \\
        0 & 1 & 0 \\
        0 & 0 & 1 \\
        0 & -r_Z & r_Y \\
        r_Z & 0 & -r_X \\
        -r_Y & r_X & 0
    \end{bmatrix}\begin{bmatrix}e_X \\ e_Y \\ e_Z \end{bmatrix}f^i \quad \Rightarrow \quad \boldsymbol{m} = \mathbf{R}_i^\text{T} f^i

where vector :math:`[e_X,\,e_Y,\,e_Z]^\text{T}` is the direction of the impact.
If we assemble all relations between virtual loads :math:`\boldsymbol{m}` and (preferably 9 or more) forces :math:`f` we obtain:

.. math::
    \boldsymbol{m} = \begin{bmatrix}\mathbf{R}_1 \\ \mathbf{R}_2 \\ \mathbf{R}_3 \\ \vdots\end{bmatrix} \boldsymbol{f} 
    \quad \Rightarrow \quad \boldsymbol{m} = \mathbf{R}^\text{T} \boldsymbol{f} 

The problem now is underdetermined and forces :math:`\boldsymbol{f}` are found by standard optimization. Optimal solution is the one that minimizes cost function 
:math:`\boldsymbol{f}^\text{T}\boldsymbol{f}` while subjected to :math:`\mathbf{R}^\text{T} \boldsymbol{f} - \boldsymbol{m} = \mathbf{0}` [6]_:

.. math::
    \boldsymbol{f} = \mathbf{R}\left(\mathbf{R}^\text{T} \mathbf{R} \right)^{-1} \boldsymbol{m} = \mathbf{T}^\text{T}_\text{f} \,\boldsymbol{m} \quad \Rightarrow \quad
    \mathbf{T}^\text{T}_\text{f} = \mathbf{R}\left(\mathbf{R}^\text{T} \mathbf{R} \right)^{-1}

Preferences in measurement data may be given using weightening matrix:

.. math::
    \boldsymbol{f} = \mathbf{W}\mathbf{R}\left(\mathbf{R}^\text{T} \mathbf{W} \mathbf{R} \right)^{-1} \boldsymbol{m} = 
    \mathbf{T}^\text{T}_\text{f}\, \boldsymbol{m} \quad \Rightarrow \quad
    \mathbf{T}^\text{T}_\text{f} = \mathbf{W}\mathbf{R}\left(\mathbf{R}^\text{T} \mathbf{W} \mathbf{R} \right)^{-1}

Virtual point admittance
************************

With both transformation matrices, virtual point FRFs can be computed using [2]_:

.. math::
    \mathbf{Y}_{\text{qm}} = \mathbf{T}_{\text{u}} \, \mathbf{Y}_{\text{uf}} \, \mathbf{T}^{\text{T}}_{\text{f}}

where :math:`\mathbf{Y}_{\text{uf}}` is the measured FRF matrix and :math:`\mathbf{Y}_{\text{qm}}` 
is the full-DoF VP FRF matrix with perfectly collocated motions and loads. That means that the VP FRF matrix should be reciprocal. 

.. figure:: ./../data/reciprocity.svg
   :width: 400px

For the driving-point FRFs at the diagonal of the matrix, passivity can be evaluated since the driving-point FRF should always be minimum-phase function:

.. math::
    \angle \mathbf{Y}_{ii} = 
    \begin{cases}
        \in [-180^{\circ},\,0^{\circ}] \quad &\text{for receptance FRFs}\\
        \in [-90^{\circ},\,90^{\circ}] \quad &\text{for mobility FRFs}\\
        \in [0^{\circ},\,180^{\circ}] \quad &\text{for accelerance FRFs}
    \end{cases}

Measurement quality indicators
******************************

Measurement quality indicators compare original responses with responses transformed into VP and then projected back to the initial location [8]_. 
A good agreement between responses indicates, that the transformation is performed correctly in a physical sense.
If vector :math:`\boldsymbol{q}` is known, expansion of responses to the :math:`\boldsymbol{u}` is also possible. However, as there is no flexible motion 
in :math:`\boldsymbol{q}`, there is also no flexible motion in reconstructed :math:`\boldsymbol{u}` for that is filtered out. 
Therefore, we introduce :math:`\boldsymbol{\tilde{u}}`, reconstructed filtered displacements:

.. math::
    \tilde{\boldsymbol{u}} = \mathbf{R}\, \boldsymbol{q} = \mathbf{R} \Big( \mathbf{R}^{\text{T}} \mathbf{W} \mathbf{R} \Big)^{-1} 
    \mathbf{R}^{\text{T}} \mathbf{W}\, \boldsymbol{u} = \mathbf{F}_{\text{u}}\,\boldsymbol{u}\quad \Rightarrow \quad \mathbf{F}_{\text{u}} 
    = \mathbf{R}\,\mathbf{T}_{\text{u}} = \mathbf{R}\,\Big( \mathbf{R}^{\text{T}} \mathbf{W} \mathbf{R} \Big)^{-1} \mathbf{R}^{\text{T}} \mathbf{W}

.. math::
    \boldsymbol{u} = \mathbf{Y} \boldsymbol{f}

.. math::
    \tilde{\boldsymbol{u}} = \mathbf{F}_\text{u} \mathbf{Y} \boldsymbol{f}
	
Overall sensor consistency is calculated from norms of filtered and non-filtered responses:

.. math::
	\boldsymbol{\rho}_\boldsymbol{u} = \frac{||\tilde{\boldsymbol{u}}||}{||\boldsymbol{u}||}
	
Values of :math:`\rho_\boldsymbol{u}` close to 1 indicate that all sensors are correctly positioned, aligned, and calibrated. 
Low values at higher frequencies indicate the presence of flexible interface modes.
Specific sensor consistency is evaluated per measurement channel with the use of coherence criterion:

.. math::
	\rho_{{u}_i} = \text{coh}(\tilde{u}_i,\,u_i) \quad \tilde{u}_i\,\in\,\tilde{\boldsymbol{u}},\,u_i\,\in\,\boldsymbol{u}
	
Using specific sensor consistency a sensor with incorrect position, direction or calibration can be identified.
A similar procedure is also used to evaluate the overall and specific quality of individual impacts. Sum of responses for each impact is defined:

.. math::
    \boldsymbol{y} = \boldsymbol{w}^\text{T}\mathbf{Y}

Sum of responses is calculated again, but this time for forces first projected onto the VP and then back-projected to the initial location:

.. math::
    \tilde{\boldsymbol{y}} = \boldsymbol{w}^\text{T}\mathbf{Y}\mathbf{F}_\text{f} \quad \Rightarrow \quad \mathbf{F}_\text{f} = \mathbf{R} \mathbf{T}_\text{f}

Overall impact consistency is calculated as follows:

.. math::
	\boldsymbol{\rho}_\boldsymbol{f} = \frac{||\tilde{\boldsymbol{y}}^\text{T}||}{||\boldsymbol{y}^\text{T}||}

A high value of the overall impact consistency indicator implies that the impact 
forces can be fully represented by the rigid IDMs. Specific impact consistency is evaluated as:

.. math::
	\rho_{{f}_i} = \text{coh}(\tilde{y}_i,\,y_i) \quad \tilde{y}_i\,\in\,\tilde{\boldsymbol{y}},\,y_i\,\in\,\boldsymbol{y}

A low specifc impact consistency values can indicate an incorrect position, direction, 
double impact or that the impact resulted in signal overload at any channel.

.. rubric:: References

.. [1] de Klerk D, Rixen DJ, Voormeeren SN, Pasteuning P. Solving the RDoF problem in experimental dynamic sybstructuring. InInternational modal analysis conference IMAC-XXVI 2008 (pp. 1-9).
.. [2] van der Seijs MV, van den Bosch DD, Rixen DJ, de Klerk D. An improved methodology for the virtual point transformation of measured frequency response functions in dynamic substructuring. In4th ECCOMAS thematic conference on computational methods in structural dynamics and earthquake engineering 2013 Jun (No. 4).
.. [3] Trainotti, Francesco, Tobias FC Berninger, and Daniel J. Rixen. Using Laser Vibrometry for Precise FRF Measurements in Experimental Substructuring. Dynamic Substructures, Volume 4. Springer, Cham, 2020. 1-11.
.. [4] Tomaž Bregar, Nikola Holeček, Gregor Čepon, Daniel J. Rixen, and Miha Boltežar. Including directly measured rotations in the virtual point transformation. Mechanical Systems and Signal Processing, 141:106440, July 2020.
.. [5] Bregar T, El Mahmoudi A, Čepon G, Rixen DJ, Boltežar M. Performance of the Expanded Virtual Point Transformation on a Complex Test Structure. Experimental Techniques. 2021 Feb;45(1):83-93.
.. [6] Häußler, Michael, and Daniel Jean Rixen. Optimal transformation of frequency response functions on interface deformation modes. Dynamics of Coupled Structures, Volume 4. Springer, Cham, 2017. 225-237.
.. [7] Pasma EA, van der Seijs MV, Klaassen SW, van der Kooij MW. Frequency based substructuring with the virtual point transformation, flexible interface modes and a transmission simulator. InDynamics of Coupled Structures, Volume 4 2018 (pp. 205-213). Springer, Cham.
.. [8] van der Seijs MV. Experimental dynamic substructuring: Analysis and design strategies for vehicle development (Doctoral dissertation, Delft University of Technology).
==============================
Singular Vector Transformation
==============================

Singular Vector Transformation (SVT) consists in projecting the acquired data into subspaces composed by dominant singular vectors, which are extracted directly from the available FRF datasets. No geometrical and/or analytical model is required. If some basic requirements are met, the reduced orthonormal frequency dependent basis would be able to control and observe most of the rigid and flexible vibration
modes of interest over a broad frequency range. The SVT can tackle challenging scenarios with flexible behaving interfaces and lightly damped systems. The method combines reduction with filtering and regularization. It shows an overall low sensitivity to measurement error and significantly reduces the condition number of the interface problem. [1]_

.. note:: 
   Download example showing the basic use of the SVT: :download:`18_SVT.ipynb <../../../examples/18_SVT.ipynb>`
   
Consider an example for the SVT where 21 impacts and 21 sensor channels are shared between the subsystem B and the B part of system AB within a decoupling application:
   
.. figure:: ./../data/SVT.png
   :width: 800px
   
******************************
Singular Vector Transformation 
******************************

The sensor channels (``df_chn_B``) and impacts (``df_imp_B``) involved in the transformation process must be assigned. The grouping number (``group``) is given to simplify the DoFs selection process.
Then, the frequency vector (``freq_B``) and the FRF matrix (``FRF_B``) from which the reduced singular subspaces are extracted need to be defined. The integer (``n_svs``) defines the chosen amount of retained singular component along the frequency range of interest.
As a result, a class instance of :class:`pyFBS.SVT` can be created.

.. code-block:: python

    svt = pyFBS.SVT(df_chn_B,df_imp_B,freq,FRF_B,group,n_svs)


The reduction matrices corresponding to displacement and force transformation from the measured to the singular domain are available as class variables ``svt.Tu`` and ``svt.Tf``.


**********************
Application of the SVT
**********************
After the reduction matrices are defined, the SVT can be applied directly on the FRF datasets involved in the substructuring scheme. 

.. code-block:: python

    _,_,FRF_B_sv= svt.apply_SVT(df_chn_B,df_imp_B,freq_B,FRF_B)
    _,_,FRF_AB_sv= svt.apply_SVT(df_chn_AB,df_imp_AB,freq_AB,FRF_AB)


*****************
Consistency check
*****************

The consistency of the applied transformation can be evaluated by comparing the measured FRFs with the SVT-filtered (reduced and back-transformed) ones.
The filtered FRFs can be obtained by using the class variables ``svt.Fu`` and ``svt.Ff``.

.. code-block:: python

    FRF_B_filt=np.zeros_like(FRF_B,dtype = complex)
    FRF_AB_filt=np.zeros_like(FRF_AB,dtype = complex)
    for i in np.arange(len(svt.freq)):
        FRF_B_filt[i,:,:]=svt.Fu[i,:,:] @ FRF_B[i,:,:]@ svt.Ff[i,:,:]
        FRF_AB_filt[i,:,:]=svt.Fu[i,:,:] @ FRF_AB[i,:,:]@ svt.Ff[i,:,:]

.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!

.. rubric:: References

.. [1] F.Trainotti, T.Bregar, S.W.B.Klaassen and D.J.Rixen. Experimental Decoupling of Substructures by Singular Vector Transformation. in: Mechanical System and Signal Processing ('Under Review'), 2021###################################
Operational Transfer Path Annalysis
###################################

Operational TPA (OPTA) coming soon!

.. note:: 
   Download example showing a numerical example of the operational TPA: :download:`15_TPA_operational.ipynb <../../../examples/17_TPA_operational.ipynb>`

..
   What is OPTA?
   ******************************

   Consider a system of substructures A and B, coupled at the interface, as depicted below.
   Substructure A is treated as an active component with operational excitation acting in :math:`\boldsymbol{u}_1`. 
   Meanwhile, no excitation force is acting on passive substructure B. 
   Responses in :math:`\boldsymbol{u}_3`, :math:`\boldsymbol{u}_4`, and also in interface DoFs :math:`\boldsymbol{u}_2` are hence a consequence of active force :math:`\boldsymbol{f}_1` only. 

   .. figure:: ./../data/in_situ.png
      :width: 250px
      :align: center

   Using LM-FBS notation, responses at the indicator sensors :math:`\boldsymbol{u}_4` can be expressed in terms of subsystem admittances [1]_:

   .. math::

      \boldsymbol{u}_4 = \textbf{Y}_{41}^{\text{AB}} \boldsymbol{f}_1 = \textbf{Y}_{42}^{\text{B}} \underbrace{ \Big(\textbf{Y}_{22}^{\text{A}} + \textbf{Y}_{22}^{\text{B}}\Big)^{-1} \textbf{Y}_{21}^{\text{A}} \boldsymbol{f}_1 }_{\boldsymbol{g}_2^{\text{B}}}.

   Expressing :math:`\boldsymbol{g}_2^{\mathrm{B}}` yields:

   .. math::

      \boldsymbol{g}_2^{\text{B}} = \Big( \textbf{Y}_{42}^{\text{B}} \Big)^+ \boldsymbol{u}_4.

   Number of indicator responses :math:`\boldsymbol{u}_4` should preferably exceed number of interface forces :math:`\boldsymbol{g}_2^{\mathrm{B}}` (or be at least equal).

   How to calculate equivalent forces?
   ***********************************

   In order to determine equivalent forces, the following steps should be performed:

   1. Measurement of admittance matrices :math:`\textbf{Y}_{42}^{\text{B}}` and :math:`\textbf{Y}_{32}^{\text{B}}` (note that for this assembly must be taken apart and only passive side is considered).
   2. Measurement of responses :math:`\boldsymbol{u}_4` on an assembly subjected to the operational excitation.

   Virtual Point Transformation
   ============================

   To simplify the measurement of the :math:`\textbf{Y}_{42}^{\text{B}}` and :math:`\textbf{Y}_{32}^{\text{B}}` the VPT can be applied on the interface excitation to transform forces at the interface into virtual DoFs (from :math:`\textbf{Y}_{\mathrm{uf}}` to :math:`\textbf{Y}_{\mathrm{um}}`): 

   .. math::

      \textbf{Y}_{\text{um}} = \textbf{Y}_{\text{uf}} \, \textbf{T}_{\text{f}}.

   For the VPT, positional data is required for channels (``df_chn_up``), impacts (``df_imp_up``) and for virtual points  (``df_vp`` and ``df_vpref``):

   .. code-block:: python

      df_imp_B = pd.read_excel(xlsx_pos, sheet_name='Impacts_B')
      df_chn_B = pd.read_excel(xlsx_pos, sheet_name='Channels_B')
      df_vp = pd.read_excel(xlsx_pos, sheet_name='VP_Channels')
      df_vpref = pd.read_excel(xlsx_pos, sheet_name='VP_RefChannels')

      vpt_B = pyFBS.VPT(df_chn_B, df_imp_B, df_vp, df_vpref)

   Defined force transformation is then applied on the FRFs and requried admittance matrices :math:`\textbf{Y}_{42}^{\text{B}}` and :math:`\textbf{Y}_{32}^{\text{B}}` are extracted as follows:

   .. code-block:: python

      Y42_B = MK_B.FRF[:,:9,:9] @ vpt_B.Tf
      Y32_B = MK_B.FRF[:,9:12,:9] @ vpt_B.Tf

   For more options and details about :mod:`pyFBS.VPT` see the :download:`04_VPT.ipynb <../../../examples/04_VPT.ipynb>` example.

   Calculation of interface forces
   ================================

   Interface forces are calculated in the following manner:

   .. code-block:: python

      g2_B = np.linalg.pinv(Y42_B) @ u4

   On-board validation
   ===================

   Finally, interface forces are applied to build up predicted response at 
   Completeness of the interface forces is then evaluated via comparison of predicted and actual response using on-board validation:

   .. code-block:: python

      u3_tpa = Y32_B @ g2_B

      o = 0

      u3 = plot_frequency_response(freq, np.hstack((u3_tpa[:,o:o+1], u3_op[:,o:o+1])))

   .. raw:: html

      <iframe src="../../_static/on_board_matrix_inverse.html" height="460px" width="100%" frameborder="0"></iframe>

   .. [1] Van der Seijs, M. V. "Experimental dynamic substructuring: Analysis and design strategies for vehicle development." (2016).#############
Pseudo-forces
#############

Pseudo-forces respresent the source excitations for responses at the passive side of the assembly. 
The same methodology as for the equivalent forces also applies for the pseudo-forces. 
They surpress the responses at the passive substructure, caused by the source running at operating conditions. 

.. note:: 
   Download example showing a numerical example of pseudo-forces: :download:`17_pseudo-forces.ipynb <../../../examples/16_TPA_pseudo_forces.ipynb>`

A set of negative signed pseudo-forces is acting at arbitrary locations on the active side of the assembly. 
They counteract the operational vibrations transmitted through the interface to the passive side [1]_.
That means that, if a set of pseudo-forces and operational excitations are applied at the system simmultaniously, 
the response at the interface DoFs :math:`\boldsymbol{u}_2` and consequently at the passive side in general equals zero.

.. figure:: ./../data/pseudo_forces_1.svg
   :width: 300px
   :align: center

Using a set of indicator response DoFs :math:`\boldsymbol{u}_4` at the passive side this can be formulated as:

.. math::

   \textbf{0} = \underbrace{\textbf{Y}_{41}^{\text{AB}} \boldsymbol{f}_1}_{\boldsymbol{u}_4} 
   + \textbf{Y}_{4\text{ps}}^{\text{AB}} \big( - \boldsymbol{f}_{\text{ps}} \big).

Expressing :math:`\boldsymbol{f}_{\mathrm{ps}}` yields:

.. math::

   \boldsymbol{f}_{\text{ps}} = \Big( \textbf{Y}_{4\text{ps}}^{\text{AB}} \Big)^+ \boldsymbol{u}_4

or a set of pseudo-forces, that are valid source descriptions for any receiver B [2]_.

.. tip::
   Pseudo-forces are property of the active component and are transferable to any assembly with modified passive side.

Applying the pseudo-forces while the source excitation is turned off 
yields the same responses as if the assembly was subjected to the :math:`\boldsymbol{f}_1`:

.. figure:: ./../data/pseudo_forces_2.svg
   :width: 300px
   :align: center

.. math::

   \boldsymbol{u}_3^{\text{TPA}} = \textbf{Y}_{3\text{ps}}^{\text{AB}} \boldsymbol{f}_{\text{ps}}.

The responses :math:`\boldsymbol{u}_3` remain independent of :math:`\boldsymbol{f}_{\mathrm{ps}}`, 
as they are not considered in the calculation of the latter.
By comparing predicted :math:`\boldsymbol{u}_3^{\mathrm{TPA}}` and measured :math:`\boldsymbol{u}_3` 
it is possible to evaluate if transfer paths through the interface are sufficiently described by :math:`\boldsymbol{f}_{\mathrm{ps}}`.

.. tip ::
   This approach can be useful for an on-board validation, when the prediction is performed on the assembly AB, 
   or a cross validation, when applied to the assembly with an modified passive side.

.. note::
   In-situ TPA is considered as a special case of pseudo-force method, where pseudo-forces are located directly at the interface DoFs.

How to calculate pseudo-forces?
*******************************

In order to determine pseudo-forces, the following steps should be performed:

1. Selection of a set of pseudo-forces :math:`\boldsymbol{f}_{\mathrm{ps}}` on the source. The locations of the pseudo-forces can be 
   chosen arbitrary, so it is reasonable to select locations that can be easily reached with impact hammer of shaker. 
   Care should be taken that the sufficient number of :math:`\boldsymbol{f}_{\mathrm{ps}}` is selected to fully replicate the source excitation.

.. tip::
   A minimum number of 6 pseudo-forces is required to fully describe rigid interface behaviour. If flexible interface behavior is not to be 
   neglected, this number should be increased reasonabely. Mind that by increasing the number of pseudo-forces, number of indicator DoFs is 
   also growing, resulting in larger experimental effort.

2. Measurement of admittance matrix :math:`\textbf{Y}_{4\text{ps}}^{\text{AB}}`.
   Often, measurement campaign is carried out on non-operating system 
   using impact hammer due to rapid FRF aquisition for each impact location.
3. Measurement of operational responses :math:`\boldsymbol{u}_4` while the assembly 
   is subjected to the operational excitation.

.. tip::
   Number of indicator responses :math:`\boldsymbol{u}_4` should preferably exceed the number of pseudo-forces :math:`\boldsymbol{f}_2^{\mathrm{eq}}`.
   An over-determination of at least a factor of 1.5 improves the results of the inverse force identification. If the amount of channels is not a limitation, a factor of 2 is suggested.
   Positions of the :math:`\boldsymbol{u}_4` should be located in the proximity of the interface and must be carefully considered to maximize the observability of the pseudo-forces. 
   If all recommendations are met, sufficient rank and low condition number of :math:`\textbf{Y}_{4\text{ps}}^{\text{AB}}` matrix prevents amplification of measurement errors in the inversion.

.. note::
   Note that no dismounting of the assembly is required in the described measurement campaign.

Once the measurement campaing is complete, pseudo-forces are simply calculated in the following manner:

.. code-block:: python

   f_ps = np.linalg.pinv(Y_4ps) @ u4_op

.. tip::

   In cases when the excitation source exhibits tonal excitation behavior, responses outside the excitation orders may fall below the noise floor of the measurement equipment. 
   The use of regularisation techniques is advisable in such cases to prevent the measurement noise from building up the pseudo-forces 
   (Singular Value Truncation or Tikhonov regularisation, for more info see [2]_ [3]_).

Finally, pseudo-forces are evaluated through on-board validation:

.. code-block:: python

   u3_tpa = Y_3ps @ f_ps

   o = 0

   u3 = plot_frequency_response(freq, np.hstack((u3_tpa[:,o:o+1], u3_op[:,o:o+1])))

.. raw:: html

   <iframe src="../../_static/u3_comparison_ps.html" height="460px" width="750px" frameborder="0"></iframe>

.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!

.. rubric:: References

.. [1] Janssens MH, Verheij JW. A pseudo-forces methodology to be used in characterization of structure-borne sound sources. Applied Acoustics. 2000 Nov 1;61(3):285-308.
.. [2] Van der Seijs, M. V. "Experimental dynamic substructuring: Analysis and design strategies for vehicle development." (2016).
.. [3] Haeussler, M. (2021). Modular sound & vibration engineering by substructuring. Technische Universität München.
====================
Classical TPA
====================
Classic TPA methods describe source excitations in terms of the interface forces between active and passive side. 
It is mainly used to troubleshoot NVH problems in existing products.

   
.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Matrix Inverse">

.. only:: html

    .. figure:: ./../data/matrix_inverse.png   
       :target: 14_matrix_inverse_TPA.html

       Matrix Inverse

.. raw:: html

    </div>

.. toctree::
   :hidden:

   14_matrix_inverse_TPA##############################
Matrix Inverse
##############################

The Matrix Inverse method determines operational interface forces between active and passive side based on structural admittance and responses at the passive side.

.. note:: 
   Download example showing a numerical example of the matrix inverse method: :download:`14_TPA_matrix_inverse.ipynb <../../../examples/13_TPA_matrix_inverse.ipynb>`

What is Matrix Inverse method?
******************************

Consider a system of substructures A and B, coupled at the interface, as depicted below.
Substructure A is treated as an active component with operational excitation acting in :math:`\boldsymbol{u}_1`. 
Meanwhile, no excitation force is acting on passive substructure B. 
Responses in :math:`\boldsymbol{u}_3`, :math:`\boldsymbol{u}_4`, and also in interface DoFs :math:`\boldsymbol{u}_2` are hence a consequence of active force :math:`\boldsymbol{f}_1` only. 

.. figure:: ./../data/in-situ.svg
   :width: 250px
   :align: center

Using LM-FBS notation, responses at the indicator sensors :math:`\boldsymbol{u}_4` can be expressed in terms of subsystem admittances [1]_:

.. math::

   \boldsymbol{u}_4 = \textbf{Y}_{41}^{\text{AB}} \boldsymbol{f}_1 = \textbf{Y}_{42}^{\text{B}} \underbrace{ \Big(\textbf{Y}_{22}^{\text{A}} + \textbf{Y}_{22}^{\text{B}}\Big)^{-1} \textbf{Y}_{21}^{\text{A}} \boldsymbol{f}_1 }_{\boldsymbol{g}_2^{\text{B}}}.

It can be seen that responses at B arise due to application of interface forces :math:`\boldsymbol{g}_2^{\mathrm{B}}` at the interface on the passive side. Expressing :math:`\boldsymbol{g}_2^{\mathrm{B}}` yields:

.. math::

   \boldsymbol{g}_2^{\text{B}} = \Big( \textbf{Y}_{42}^{\text{B}} \Big)^+ \boldsymbol{u}_4.

Responses at the passive side can than be predicted based on :math:`\boldsymbol{g}_2^{\mathrm{B}}`:

.. math::

   \tilde{\boldsymbol{u}}_3 = \Big( \textbf{Y}_{32}^{\text{B}} \Big)^+ \boldsymbol{g}_2^{\text{B}}.

.. tip::

   By comparing predicted :math:`\tilde{\boldsymbol{u}}_3` and measured :math:`\boldsymbol{u}_3` it is possible to evaluate if transfer paths through the 
   interface are sufficiently described by :math:`\boldsymbol{g}_2^{\mathrm{B}}`.

How to calculate interface forces?
***********************************

In order to determine interface forces :math:`\boldsymbol{g}_2^{\mathrm{B}}`, the following steps should be performed:

1. Measurement of admittance matrices :math:`\textbf{Y}_{42}^{\text{B}}` and :math:`\textbf{Y}_{32}^{\text{B}}` (note that for this assembly must be taken apart and only passive side is considered).
2. Measurement of responses :math:`\boldsymbol{u}_4` on an assembly subjected to the operational excitation.

.. tip::
   Number of indicator responses :math:`\boldsymbol{u}_4` should preferably exceed the number of interface forces :math:`\boldsymbol{g}_2^{\mathrm{B}}` (or be at least equal but this is not recommended).
   An over-determination of at least a factor of 1.5 improves the results of the inverse force identification. If the amount of channels is not a limitation, a factor of 2 is suggested.
   Positions of the :math:`\boldsymbol{u}_4` should be located in the proximity of the interface and must be carefully considered to maximize the observability of the interface forces. 
   If all recommendations are met, sufficient rank and low condition number of :math:`\textbf{Y}_{42}^{\text{B}}` matrix prevents amplification of measurement errors in the inversion.

Virtual Point Transformation
============================

.. tip::
   The virtual point [2]_, typically used in frequency based substructuring (FBS) applications, has the advantage of taking into account moments in the transfer paths that are otherwise not measurable with
   conventional force transducers. Hence the description of the interface is more complete.

To simplify the measurement of the :math:`\textbf{Y}_{42}^{\text{B}}` and :math:`\textbf{Y}_{32}^{\text{B}}` the VPT can be applied on the interface excitation to transform forces at the interface into virtual DoFs (from :math:`\textbf{Y}_{\mathrm{uf}}` to :math:`\textbf{Y}_{\mathrm{um}}`): 

.. math::

   \textbf{Y}_{\text{um}} = \textbf{Y}_{\text{uf}} \, \textbf{T}_{\text{f}}.

For the VPT, positional data is required for channels (``df_chn_up``), impacts (``df_imp_up``) and for virtual points  (``df_vp`` and ``df_vpref``):

.. code-block:: python

   df_imp_B = pd.read_excel(xlsx_pos, sheet_name='Impacts_B')
   df_chn_B = pd.read_excel(xlsx_pos, sheet_name='Channels_B')
   df_vp = pd.read_excel(xlsx_pos, sheet_name='VP_Channels')
   df_vpref = pd.read_excel(xlsx_pos, sheet_name='VP_RefChannels')

   vpt_B = pyFBS.VPT(df_chn_B, df_imp_B, df_vp, df_vpref)

Defined force transformation is then applied on the FRFs and requried admittance matrices :math:`\textbf{Y}_{42}^{\text{B}}` and :math:`\textbf{Y}_{32}^{\text{B}}` are extracted as follows:

.. code-block:: python

   Y42_B = MK_B.FRF[:,:9,:9] @ vpt_B.Tf
   Y32_B = MK_B.FRF[:,9:12,:9] @ vpt_B.Tf

For more options and details about :mod:`pyFBS.VPT` see the :download:`04_VPT.ipynb <../../../examples/04_VPT.ipynb>` example.

Calculation of interface forces
================================

Interface forces are calculated in the following manner:

.. code-block:: python

   g2_B = np.linalg.pinv(Y42_B) @ u4

On-board validation
===================

Finally, interface forces are applied to build up predicted response at the passive side. 
Completeness of the interface forces is then evaluated via comparison of predicted and actual response using on-board validation:

.. code-block:: python

   u3_tpa = Y32_B @ g2_B

   o = 0

   u3 = plot_frequency_response(freq, np.hstack((u3_tpa[:,o:o+1], u3_op[:,o:o+1])))

.. raw:: html

   <iframe src="../../_static/on_board_matrix_inverse.html" height="460px" width="100%" frameborder="0"></iframe>

.. warning::

   Poor agreement between :math:`\boldsymbol{u}_3` and :math:`\tilde{\boldsymbol{u}}_3`  indicates that there might be other sources that significantly contribute to the target output, 
   the predefined forces are not correct or the on-board validation sensor is too far from the interface.

.. tip::

   In cases when the excitation source exhibits tonal excitation behavior, responses outside the excitation orders may fall below the noise floor of the measurement equipment. 
   The use of regularisation techniques is advisable in such cases to prevent the measurement noise from building up the interface forces 
   (Singular Value Truncation or Tikhonov regularisation, for more info see [3]_ [4]_).

.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!

.. rubric:: References

.. [1] van der Seijs MV. Experimental dynamic substructuring: Analysis and design strategies for vehicle development (Doctoral dissertation, Delft University of Technology).
.. [2] van der Seijs MV, van den Bosch DD, Rixen DJ, de Klerk D. An improved methodology for the virtual point transformation of measured frequency response functions in dynamic substructuring. In4th ECCOMAS thematic conference on computational methods in structural dynamics and earthquake engineering 2013 Jun (No. 4).
.. [3] Thite AN, Thompson DJ. The quantification of structure-borne transmission paths by inverse methods. Part 1: Improved singular value rejection methods. Journal of Sound and Vibration. 2003 Jul 3;264(2):411-31.
.. [4] Thite AN, Thompson DJ. The quantification of structure-borne transmission paths by inverse methods. Part 2: Use of regularization techniques. Journal of Sound and Vibration. 2003 Jul 3;264(2):433-51.##############################
in-situ Transfer Path Analysis
##############################

The in-situ Transfer Path Analysis is a method that utilizes equivalent forces to describe operational excitations [1]_.
With the possibility to perform operational measurements on the target assembly, dismounting of any part can be avoided.

.. note:: 
   Download example showing a numerical example of the in-situ TPA: :download:`11_insitu_TPA.ipynb <../../../examples/11_TPA_in-situ.ipynb>`

What is in-situ TPA?
********************

Consider a system of substructures A and B, coupled at the interface, as depicted below.
Substructure A is treated as an active component with operational excitation acting in :math:`\boldsymbol{u}_1`. 
Meanwhile, no excitation force is acting on passive substructure B. 
Responses in :math:`\boldsymbol{u}_3`, :math:`\boldsymbol{u}_4`, and also in interface DoFs :math:`\boldsymbol{u}_2` are hence a consequence of active force :math:`\boldsymbol{f}_1` only. 

.. figure:: ./../data/in-situ.svg
   :width: 250px
   :align: center

Source internal structure-borne excitations :math:`\boldsymbol{f}_1` are often unmeasurable in practice.
In-situ TPA introduces the set of equivalent forces, acting on interface DoFs, that cause the same displacements on B as :math:`\boldsymbol{f}_1`.
Therefore, application of forces :math:`\boldsymbol{f}_1` and reaction forces :math:`-\boldsymbol{f}_2^{\mathrm{eq}}` should annul any response on the passive side, e.q. for :math:`\boldsymbol{u}_4`:

.. figure:: ./../data/in-situ_blocked.svg
   :width: 250px
   :align: center

.. math::

   \textbf{0} = \underbrace{\textbf{Y}_{41}^{\text{AB}} \boldsymbol{f}_1}_{\boldsymbol{u}_4} + \textbf{Y}_{42}^{\text{AB}} \big( - \boldsymbol{f}_2^{\text{eq}} \big).

Expressing :math:`\boldsymbol{f}_2^{\mathrm{eq}}` yields:

.. math::

   \boldsymbol{f}_2^{\text{eq}} = \Big( \textbf{Y}_{42}^{\text{AB}} \Big)^+ \boldsymbol{u}_4

or a set of equivalent forces, that are valid source descriptions for any receiver B [2]_.

.. tip::
   Equivalent forces from in-situ TPA are property of the active component and are transferable to any assembly with modified passive side.

TPA methods offer a useful tool to assess the completeness of the interface description in a form of on-board validation [2]_.
The responses :math:`\boldsymbol{u}_3` remain independent of :math:`\boldsymbol{f}_2^{\mathrm{eq}}`, as they are not considered in the calculation of the latter.
Therefore, response in :math:`\boldsymbol{u}_3` can be predicted based on :math:`\boldsymbol{f}_2^{\mathrm{eq}}`:

.. figure:: ./../data/in-situ_eq_source.svg
   :width: 250px
   :align: center

.. math::

   \boldsymbol{u}_3^{\text{TPA}} = \textbf{Y}_{32}^{\text{AB}} \boldsymbol{f}_2^{\text{eq}}.

By comparing predicted :math:`\boldsymbol{u}_3^{\mathrm{TPA}}` and measured :math:`\boldsymbol{u}_3` 
it is possible to evaluate if transfer paths through the interface are sufficiently described by :math:`\boldsymbol{f}_2^{\mathrm{eq}}`.

.. tip ::
   This approach can be useful for an on-board validation, when the prediction is performed on the assembly AB, 
   or a cross validation, when applied to the assembly with an modified passive side.

How to calculate equivalent forces?
***********************************

In order to determine equivalent forces, the following steps should be performed:

1. Measurement of admittance matrices :math:`\textbf{Y}_{42}^{\text{AB}}` and :math:`\textbf{Y}_{32}^{\text{AB}}`.
   Often, measurement campaign is carried out on non-operating system 
   using impact hammer due to rapid FRF aquisition for each impact location.
2. Measurement of operational responses :math:`\boldsymbol{u}_4` while the assembly 
   is subjected to the operational excitation.

.. warning::
   To ensure that the equivalent forces are independent of the receiver structure, the operating excitation must originate solely from the source structure.

.. tip::
   Number of indicator responses :math:`\boldsymbol{u}_4` should preferably exceed the number of equivalent forces :math:`\boldsymbol{f}_2^{\mathrm{eq}}`.
   An over-determination of at least a factor of 1.5 improves the results of the inverse force identification. If the amount of channels is not a limitation, a factor of 2 is suggested.
   Positions of the :math:`\boldsymbol{u}_4` should be located in the proximity of the interface and must be carefully considered to maximize the observability of the equivalent forces. 
   If all recommendations are met, sufficient rank and low condition number of :math:`\textbf{Y}_{42}^{\text{AB}}` matrix prevents amplification of measurement errors in the inversion.

.. warning::
   One should keep in mind that we assume considered systems are linear. Especially with respect to systems with nonlinear interface properties, 
   the limitations of the method should be kept in mind. Because the FRF determination through impact hammer or shaker is done during non-operation, operation FRFs might differ as 
   non-linear components are changing their dynamic behaviour depending on the operational speed. 
   Consequently, errors arise in equivalent forces and their transferabillity is limited.

Virtual Point Transformation
============================

.. tip::
   The virtual point [3]_, typically used in frequency based substructuring (FBS) applications, has the advantage of taking into account moments in the transfer paths that are otherwise not measurable with
   conventional force transducers. Hence the description of the interface is more complete [4]_.

To simplify the measurement of the :math:`\textbf{Y}_{42}^{\text{AB}}` and :math:`\textbf{Y}_{32}^{\text{AB}}` the VPT can be applied on the interface excitation to transform forces at the interface into virtual DoFs (from :math:`\textbf{Y}_{\mathrm{uf}}` to :math:`\textbf{Y}_{\mathrm{um}}`): 

.. math::

   \textbf{Y}_{\text{um}} = \textbf{Y}_{\text{uf}} \, \textbf{T}_{\text{f}}^\text{T}.

For the VPT, positional data is required for channels (``df_chn_up``), impacts (``df_imp_up``) and for virtual points  (``df_vp`` and ``df_vpref``). 
Note that only transformation of forces is required for succesful identification of equivalent forces using in-situ TPA, 
but in order to define the VPT object in pyFBS, VP for displacements must also be defined (but neglected latter in the trasnformation):

.. code-block:: python

   df_chn_AB = pd.read_excel(pos_xlsx, sheet_name='Channels_AB')
   df_imp_AB = pd.read_excel(pos_xlsx, sheet_name='Impacts_AB')
   df_vp = pd.read_excel(pos_xlsx, sheet_name='VP_Channels')
   df_vpref = pd.read_excel(pos_xlsx, sheet_name='VP_RefChannels')

   vpt_AB = pyFBS.VPT(df_chn_AB_up,df_imp_AB_up,df_vp,df_vpref)

Defined force transformation is then applied on the FRFs:

.. code-block:: python

   Y_um = Y_uf @ vpt_AB.Tf

and requried admittance matrices :math:`\textbf{Y}_{42}^{\text{AB}}` and :math:`\textbf{Y}_{32}^{\text{AB}}` are extracted as follows:

.. code-block:: python

   Y_42 = Y_um[:,:9,:6]
   Y_32 = Y_um[:,9:,:6]

Therefore, the interface is loaded with three forces (:math:`f_x,\,f_y,\,f_z`) and three moments (:math:`m_x,\,m_y,\,m_z`). 
Consistency of the VPT can be additionally evaluated using specific and overall impact consistency:

.. code-block:: python

   vpt_AB.consistency([1],[1])
   barchart(np.arange(1,10,1), vpt_AB.specific_impact, title='Specific Impact Consistency')
   plot_coh(freq, vpt_AB.overall_impact, title='Overall Impact Consistency')

.. raw:: html

   <iframe src="../../_static/specific_impact_consistency.html" height="300px" width="295px" frameborder="0"></iframe>
   <iframe src="../../_static/overall_impact_consistency.html" height="300px" width="595px" frameborder="0"></iframe>

For more options and details about :mod:`pyFBS.VPT` see the :download:`04_VPT.ipynb <../../../examples/04_VPT.ipynb>` example.

Calculation of equivalent forces
================================

Equivalent forces at the interface are calculated in the following manner:

.. code-block:: python

   f_eq = np.linalg.pinv(Y_42) @ u4_op

.. tip::

   In cases when the excitation source exhibits tonal excitation behavior, responses outside the excitation orders may fall below the noise floor of the measurement equipment [5]_. 
   The use of regularisation techniques is advisable in such cases to prevent the measurement noise from building up the equivalent forces 
   (Singular Value Truncation or Tikhonov regularisation, for more info see [6]_ [7]_ [8]_).

On-board validation
===================

Finally, equivalent forces are evaluated through on-board validation:

.. code-block:: python

   u3_tpa = Y_32 @ f_eq

   o = 0

   u3 = plot_frequency_response(freq, np.hstack((u3_tpa[:,o:o+1], u3_op[:,o:o+1])))

.. raw:: html

   <iframe src="../../_static/u3_comparison.html" height="460px" width="750px" frameborder="0"></iframe>

Additionally, a coherence criterion can be used to objectively evaluate interface completeness:

.. code-block:: python 

   coh_data = coh(u3_tpa, u3_op)
   plot_coh_group(freq, coh_data)

.. raw:: html

   <iframe src="../../_static/on_board_coherence.html" height="340px" width="750px" frameborder="0"></iframe>

See also cross-validation for further evaluation of the equivalent forces completeness [9]_.

Partial response contribution
=============================

If we now focus on only one response at the passive side at :math:`\boldsymbol{u}_3`, 
we can build this responses solely by the equivalent forces. 
Partial response is calculated for each equivalent force using equation:

.. math::

   u_{i,j} = Y_{ij}^\mathrm{AB} f_j^\text{eq}

Then if we sum all partial :math:`u_{i,j}` we obtain total response of the passive side [2]_. 
But it is also interesting to examine all partial responses individually. 
Based on them, we can determine, which equivalent force contributes the most to the passive side responses, 
or in other words, which transfer path is most critical on our product.

.. code-block:: python

   sel_i = 1
   u_partial = []

   for j in range(6):
      gg = _Y_temp[:,sel_i:sel_i+1,j:j+1] @ f_eq[:,j:j+1,0:1]
      u_partial.append(gg[:,0,0])
      
   u_partial = np.asarray(u_partial).T

The partial responses can then be displayed as a heatmap:

.. raw:: html

   <iframe src="../../_static/TP_contribution.html" height="230px" width="750px" frameborder="0"></iframe>

.. tip::
   Using the graphical presentation above, the most dominant transfer path can be pinpointed.

.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!

.. rubric:: References

.. [1] Moorhouse AT. On the characteristic power of structure-borne sound sources. Journal of sound and vibration. 2001 Nov 29;248(3):441-59.
.. [2] van der Seijs MV. Experimental dynamic substructuring: Analysis and design strategies for vehicle development (Doctoral dissertation, Delft University of Technology).
.. [3] van der Seijs MV, van den Bosch DD, Rixen DJ, de Klerk D. An improved methodology for the virtual point transformation of measured frequency response functions in dynamic substructuring. In4th ECCOMAS thematic conference on computational methods in structural dynamics and earthquake engineering 2013 Jun (No. 4).
.. [4] Haeussler, M., Mueller, T., Pasma, E. A., Freund, J., Westphal, O., Voehringer, T., & AG, Z. F. (2020). Component TPA: benefit of including rotational degrees of freedom and over-determination. In ISMA 2020-International Conference on Noise and Vibration Engineering (pp. 1135-1148).
.. [5] Haeussler, M., Kobus, D. C., & Rixen, D. J. (2021). Parametric design optimization of e-compressor NVH using blocked forces and substructuring. Mechanical Systems and Signal Processing, 150, 107217.
.. [6] Wernsen, M. W. F., van der Seijs, M. V., & de Klerk, D. (2017). An indicator sensor criterion for in-situ characterisation of source vibrations. In Sensors and Instrumentation, Volume 5 (pp. 55-69). Springer, Cham.
.. [7] Haeussler, M. (2021). Modular sound & vibration engineering by substructuring. Technische Universität München.
.. [8] Wernsen MW. Observability and transferability of in-situ blocked force characterisation.
.. [9] El Mahmoudi, A., Trainotti, F., Park, K., & Rixen, D. J. (2019). In-situ TPA for NVH analysis of powertrains: an evaluation on an experimental test setup. In AAC 2019: Aachen acoustics colloquium/aachener akustik kolloquium.#############
Free velocity
#############

Using free velocity concept equivalent forces can be expressed from the interface motion of the source component 
while operating in free-free conditions.

.. note:: 
   Download example showing a numerical example of the free velocity: :download:`16_free_velocity.ipynb <../../../examples/15_TPA_free_velocity.ipynb>`

Free velocity concept
*********************

When the interface of the source structure is left free, vibrations at the interface DoFs can be treated as a 
free displacements :math:`\boldsymbol{u}_2^{\mathrm{free}}`.

.. figure:: ./../data/free_velocity.svg
   :width: 160px
   :align: center

By definition [1]_, equivalent forces, applied in the opposite direction must cancel out responses at the interface 
caused by the source in operation:

.. math::
   \mathbf{0} = \underbrace{\mathbf{Y}_{21}^\text{A} \boldsymbol{f}_1}_{\boldsymbol{u}_2^\text{free}} + \mathbf{Y}_{22}^\text{A} (- \boldsymbol{f}_2^\text{eq})

To derive the equivalent forces :math:`\boldsymbol{f}_2^{\mathrm{eq}}` from the free velocities the 
free admittance matrix of the uncoupled source component :math:`\mathbf{Y}_{22}^{\text{A}}` is needed:

.. math::

   \boldsymbol{f}_2^{\mathrm{eq}} = \left( \mathbf{Y}_{22}^\text{A} \right)^{-1} \boldsymbol{u}_2^\text{free}

.. tip::
   Equivalent forces are property of the active component and are transferable to any assembly with modified passive side [1]_.

How to calculate equivalent forces?
***********************************

In order to determine equivalent forces, the following steps should be performed:

1. Measurement of admittance matrix :math:`\textbf{Y}_{22}^{\text{A}}`.
   Often, measurement campaign is carried out on non-operating system
   using impact hammer due to rapid FRF aquisition for each impact location.
2. Measurement of free velocities :math:`\boldsymbol{u}_2^\text{free}` while the source structure 
   is subjected to the operational excitation.

.. tip::
   At first glance measurement campaign looks fairly simple. Both required measurements have to be performed on a sorce structure only,
   so the experimental effort is quite low. However...

.. warning::
   ...running source objects at free-free conditions is often challenging as some of them require some sort of support to be able to run in 
   operation. The active components often needs to be connected to a certain load or mount for operating.

.. tip::
   In practice, source description using free velocities is limited for lower frequency range due to 
   unreal running conditions if the interface is left free. 
   Hence free velocity concept is more suited for frequency range well above rigid body modes of the source. 
   The method is expected to perform best when differences between operational and free accelerations are small. 
   See [2]_ for some practical considerations.

Virtual Point Transformation
============================

.. tip::
   The virtual point [3]_, typically used in frequency based substructuring (FBS) applications, has the advantage of taking into account moments in the transfer paths that are otherwise not measurable with
   conventional force transducers. Hence the description of the interface is more complete.

To simplify the measurement of the :math:`\textbf{Y}_{22}^{\text{A}}` the VPT can be applied on the interface excitation to 
transform displacements and forces at the interface into virtual DoF: 

.. math::

   \textbf{Y}_{\text{qm}} = \textbf{T}_{\text{u}} \, \textbf{Y}_{\text{uf}} \, \textbf{T}_{\text{f}}^\text{T}.

For the VPT, positional data is required for channels (``df_chn_up``), impacts (``df_imp_up``) and for virtual points  (``df_vp`` and ``df_vpref``):

.. code-block:: python

   df_acc_A = pd.read_excel(xlsx_pos, sheet_name='Sensors_A')
   df_chn_A = pd.read_excel(xlsx_pos, sheet_name='Channels_A')
   df_imp_A = pd.read_excel(xlsx_pos, sheet_name='Impacts_A')

   df_vp = pd.read_excel(xlsx_pos, sheet_name='VP_Channels')
   df_vpref = pd.read_excel(xlsx_pos, sheet_name='VP_RefChannels')

After the reduction matrices are defined the VPT can be applied directly on an FRF matrix:

.. code-block:: python

   vpt = pyFBS.VPT(df_chn_A_up, df_imp_A_up, df_vp, df_vpref)
   vpt.apply_VPT(MK_A.freq, MK_A.FRF)

For more options and details about :mod:`pyFBS.VPT` see the :download:`04_VPT.ipynb <../../../examples/04_VPT.ipynb>` example.

Calculation of equivalent forces
================================

Equivalent forces at the interface are calculated in the following manner:

.. code-block:: python

   f_eq = np.linalg.pinv(Y_42) @ u2_free

.. tip::

   In cases when the excitation source exhibits tonal excitation behavior, responses outside the excitation orders may fall below the noise floor of the measurement equipment. 
   The use of regularisation techniques is advisable in such cases to prevent the measurement noise from building up the equivalent forces 
   (Singular Value Truncation or Tikhonov regularisation, for more info see [1]_ [4]_).

Cross validation
================

.. tip::

   TPA methods offer a useful tool to assess the completeness of the source description in a form of cross validation. 
   
As stated previously, equivalent forces are a property of the source only and are thus transferable to any assembly with modified passive side. 
Response of the new assembly when subjected to the operational excitation of the source (:math:`\boldsymbol{u}_3`) 
can be predicted based on :math:`\boldsymbol{f}_2^{\mathrm{eq}}` identified from source in free conditions
and admittance of the new assembly (:math:`\mathbf{Y}_{32}^{\mathrm{AB}}`):   

.. math::
   \boldsymbol{u}_3 = \mathbf{Y}_{32}^{\mathrm{AB}}\, \boldsymbol{f}_2^{\mathrm{eq}}

.. code-block:: python

   u3 = Y32_AB @ f2_eq

.. raw:: html

   <iframe src="../../_static/free_velocity_cross.html" height="460px" width="750px" frameborder="0"></iframe>

.. tip::
   We can see that equivalent forces are indeed independent of the passive substructure.

.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!

.. rubric:: References

.. [1] Van der Seijs, M. V. "Experimental dynamic substructuring: Analysis and design strategies for vehicle development." (2016).
.. [2] Wagner P, Bianciardi F, Corbeels P, Hülsmann A. High frequency source characterization of an e-motor using component-based TPA.
.. [3] van der Seijs MV, van den Bosch DD, Rixen DJ, de Klerk D. An improved methodology for the virtual point transformation of measured frequency response functions in dynamic substructuring. In4th ECCOMAS thematic conference on computational methods in structural dynamics and earthquake engineering 2013 Jun (No. 4).
.. [4] Haeussler, M. (2021). Modular sound & vibration engineering by substructuring. Technische Universität München.
==========================
Transmissibility-based TPA
==========================
The family of transmissibility-based TPA determines sound and vibration transfer from transmissibilities between sensors while assembly is subjected to operational excitation. 
This approach is convinient when TPA is used solely to identify the dominant path contributions in existing products. 

   
.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Operational TPA">

.. only:: html

    .. figure:: ./../data/labels.png   
       :target: 15_operational_TPA.html

       Operational TPA

.. raw:: html

    </div>

.. toctree::
   :hidden:

   15_operational_TPA====================
Component-based TPA
====================
For an independent characterization of the source structure, component-based TPA adopts a set of equivalent forces applied at the interface between active and passive side. 
These equivalent or blocked forces are valid for any assembly with a modified passive side.

   
.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="in-situ TPA">

.. only:: html

    .. figure:: ./../data/tpa.png   
       :target: 11_insitu_TPA.html

       in-situ TPA

.. raw:: html

    </div>

.. toctree::
   :hidden:

   11_insitu_TPA





.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Free velocity">

.. only:: html

    .. figure:: ./../data/free_velocity.PNG   
       :target: 16_free_velocity.html

       Free velocity

.. raw:: html

    </div>

.. toctree::
   :hidden:

   16_free_velocity




.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="Pseudo-forces">

.. only:: html

    .. figure:: ./../data/pseudo_forces.PNG   
       :target: 17_pseudo_forces.html

       Pseudo-forces

.. raw:: html

    </div>

.. toctree::
   :hidden:

   17_pseudo_forces#################################
Experimental modal analysis (EMA)
#################################

The pyFBS can be also used with other Python packages for structural dynamics. 
One of those packages is `pyEMA <https://pypi.org/project/pyEMA/>`_ which can be used to perform an Experimental Modal Analysis (EMA). 
In this example an integration of the two packages is shown on a frame of the automotive testbench example. 
Real experimental data is used in this example and it is also available directly within the pyFBS.

.. note:: 
   Download example showing an Experimental Modal Analysis (EMA) application: :download:`09_EMA.ipynb <../../../examples/09_experimental_modal_analysis_EMA.ipynb>`

Example Datasests and 3D display
********************************

Load the required predefined datasets. Open a 3Dviewer in the background. Add the STL file of the assembly to the 3D display:

.. figure:: ./../data/nine_one.png
   :width: 800px
   
Experimental Modal analysis - Assembly
======================================
Experimental modal analysis in Python can be performed with pyEMA package. Currently single reference modal identification with LSCD/LSFD is supported. For more details check the pyEMA documentation.

pyEMA
*****

Perform the LSCF/LSFD experimental identification of modal parameters:

.. code-block:: python

    Y =  Y_m1[:,:,0].T

    modal_1 = pyEMA.Model(Y,freq,pol_order_high=60,lower = 0,upper = 1000)
    modal_1.get_poles()
    modal_1.select_poles()

    H_acc, modes_1 = modal_1.get_constants(whose_poles=modal_1,least_squares_type="old")
    pd.DataFrame({"Nat. freq [Hz]": modal_1.nat_freq,"Damping [/]": modal_1.nat_xi},index = np.arange(len(modal_1.nat_freq))+1)

Stable poles are selected using stability chart:

.. figure:: ./../data/nine_two.png
   :width: 800px
   
The following eigenfrequencies and corresponding damping ratios are obtained:
   
.. figure:: ./../data/nine_three.png
   :width: 200px
   
For the animation a mesh can be manually created. In this example line connections between the points are made. For more details on :mod:`pv.PolyData` check an example from `PyVista <https://docs.pyvista.org/index.html>`_.

.. code-block:: python

	pos_array = df_acc[["Position_1","Position_2","Position_3"]].to_numpy()*1000
	faces = np.hstack([[2,1,2],[2,3,2],[2,3,13],[2,0,13],[2,0,1],[2,0,1],[2,5,1],[2,5,12],[2,4,12],[2,9,4],
	                   [2,9,8],[2,8,7],[2,0,4],[2,5,11],[2,10,11],[2,10,6],[2,2,6],[2,3,7]]).astype(np.int8)
	point_cloud_1 = pv.PolyData(pos_array,faces)
	pts_1 = point_cloud_1.points.copy()

	view3D_1.plot.add_mesh(point_cloud_1,name = "mesh",render_lines_as_tubes = True, line_width=10, color = "k",clim = [-1,1], cmap="coolwarm",scalars = np.zeros(pts_1.shape[0]),style = "wireframe");
    
Mode shape animation
********************
The identified mode shape can be animated directly in the 3D view:

.. code-block:: python

    mode_select_1 = 0

    emp_1 = pyFBS.orient_in_global(modes_1[:,mode_select_1],df_chn,df_acc)
    mode_dict = pyFBS.dict_animation(emp_1,"modeshape",pts = pts_1,mesh = point_cloud_1,r_scale = 50)

    mode_dict["freq"] = modal_1.nat_freq[mode_select_1]
    mode_dict["damp"] = modal_1.nat_xi[mode_select_1]
    mode_dict["mcf"] = pyFBS.MCF(modes_1[:,mode_select_1])

    view3D_1.add_modeshape(mode_dict,run_animation = True,add_note = True)

.. figure:: ./../data/modal_1_min.gif
   :width: 800px
   
EMA - Assembly without the source
=================================

3D View
*******
You can open multiple displays at once and performs simoultenous analyses. Open a second 3D display:

.. code-block:: python

    view3D_2 = pyFBS.view3D()
    
Add the STL files of the assembly without the source structure to the 3D view:

.. code-block:: python

    view3D_2.add_stl(stl_rec,name = "receiver_0",color = "#e0e0e0",opacity = .1)
    view3D_2.add_stl(stl_tm,name = "transmission_mount_0",color = "#83afd2",opacity = .1)
    view3D_2.add_stl(stl_rm,name = "roll_mount_0",color = "#83afd2",opacity = .1)
    view3D_2.add_stl(stl_em,name = "engine_mount_0",color = "#83afd2",opacity = .1);
 
.. figure:: ./../data/nine_one_1.png
   :width: 800px

pyEMA
*****
Perform the LSCF/LSFD experimental identification of modal parameters:

.. code-block:: python

    # select the reference DoF
    Y =  Y_m2[:,:,0].T

    modal_2 = pyEMA.Model(Y,freq,pol_order_high=60,lower = 0,upper = 1000)
    modal_2.get_poles()
    modal_2.select_poles()
    H_acc, modes_2 = modal_2.get_constants(whose_poles=modal_2,least_squares_type="old")
    pd.DataFrame({"Nat. freq [Hz]": modal_2.nat_freq,"Damping [/]": modal_2.nat_xi},index = np.arange(len(modal_2.nat_freq))+1)

Stable poles are selected using stability chart:

.. figure:: ./../data/nine_four.png
   :width: 800px

The following eigenfrequencies and corresponding damping ratios are obtained:
   
.. figure:: ./../data/nine_five.png
   :width: 200px
   
Create a mesh and add it to the 3D view:
****************************************

.. code-block:: python

    pos_array = df_acc[["Position_1","Position_2","Position_3"]].to_numpy()*1000
    
    point_cloud_2 = pv.PolyData(pos_array,faces)
    pts_2 = point_cloud_2.points.copy()
    view3D_2.plot.add_mesh(point_cloud_2,name = "mesh",render_lines_as_tubes = True, line_width=10, color = "k",clim = [-1,1], cmap="coolwarm",scalars = np.zeros(pts_2.shape[0]),style = "wireframe")
    
Mode shape animation
********************
The second set of identified mode shapes can be animated directly in the 3D view. As two instances of the pyFBS.view3D can be created in the same analysis a side-by-side comparison is possible:

.. code-block:: python

    mode_select_2 = 0

    emp_2 = pyFBS.orient_in_global(modes_2[:,mode_select_2],df_chn,df_acc)
    mode_dict = pyFBS.dict_animation(emp_2,"modeshape",pts = pts_2,mesh = point_cloud_2,r_scale = 50)

    mode_dict["freq"] = modal_2.nat_freq[mode_select_2]
    mode_dict["damp"] = modal_2.nat_xi[mode_select_2]
    mode_dict["mcf"] = pyFBS.MCF(modes_2[:,mode_select_2])

    view3D_2.add_modeshape(mode_dict,run_animation = True,add_note = True)
    
.. figure:: ./../data/modal_2_min.gif
   :width: 800px

.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!#############################
Transmission simulator in FBS
#############################

The method of Transmission Simulator (TS) is primarly used direclty in the modal domain. Nevertheless, the concept of attaching a transmission simulator at the interface can be used also in the frequency domain. In this example a frequency based subtructuring is performed on a complex structure. First a transmission simulator is decoupled from the receiver structure and aftewards a source structure is coupled to the receiver structure. At the interface the VPT is used to obtain collocated interface DoFs.

.. note:: 
   Download example showing a Transmission Simulator (TS) application: :download:`10_TS.ipynb <../../../examples/10_FBS_transmission_simulator.ipynb>`.

Example Datasests and 3D view
*****************************
As already shown in the `3D Display <../../../html/examples/basic_examples/01_static_display.html>`_ one can load the predefined datasets from an example and add a structure from STL file to the 3D display. This allows both the acceleration sensors and excitation points to be visualized. Also for the transmission simulator application, a subplot representation, as already presented in `Coupling <../../../html/examples/fbs/07_coupling.html>`_, can be used.
	
.. figure:: ./../data/ten_display_four.png
   :width: 500px
   
Virtual point transformation
****************************
The VPT can be performed on experimental data. See the :download:`04_VPT.ipynb <../../../examples/04_VPT.ipynb>` example for more options and details.

.. code-block:: python

	df_vp = pd.read_excel(xlsx_A, sheet_name='VP Channels')
	df_vpref = pd.read_excel(xlsx_A, sheet_name='VP RefChannels')

	vpt_TS = pyFBS.VPT(df_chn_TS,df_imp_TS,df_vp,df_vpref,sort_matrix = False)
	vpt_BTS = pyFBS.VPT(df_chn_BTS,df_imp_BTS,df_vp,df_vpref,sort_matrix = False)
	vpt_A = pyFBS.VPT(df_chn_A,df_imp_A,df_vp,df_vpref,sort_matrix = False)

	vpt_TS.apply_VPT(freq,Y_TS)
	vpt_BTS.apply_VPT(freq,Y_BTS)
	vpt_A.apply_VPT(freq,Y_A)

	Y_TS_tran = vpt_TS.vptData
	Y_BTS_tran = vpt_BTS.vptData
	Y_A_tran = vpt_A.vptData
	
LM-FBS coupling and decoupling
==============================
First, construct an admittance matrix for the uncoupled system, containing substructure admittances:

.. code-block:: python

	Y_AnB = np.zeros((800,30+36+18,30+36+18),dtype = complex)

	Y_AnB[:,0:36,0:36] = Y_BTS_tran
	Y_AnB[:,36:36+30,36:36+30] = -1*Y_TS_tran
	Y_AnB[:,30+36:,30+36:] = Y_A_tran


Next the compatibility and the equilibrium conditions has to be defined through the signed Boolean matrices ``Bu`` and ``Bf``. 

.. code-block:: python

	k = 18
	Bu = np.zeros((2*k,36+30+18))
	Bu[:k,0:k] = 1*np.eye(k)
	Bu[:k,36:36+k] = -1*np.eye(k)

	Bu[k:,0:k] = 1*np.eye(k)
	Bu[k:,36+30:36+30+k] = -1*np.eye(k)

	Bf = np.zeros((2*k,36+30+18))
	Bf[:k,0:k] = 1*np.eye(k)
	Bf[:k,36:36+k] = -1*np.eye(k)

	Bf[k:,0:k] = 1*np.eye(k)
	Bf[k:,36+30:36+30+k] = -1*np.eye(k)
	
.. figure:: ./../data/Bu_TS.png
   :width: 500px

.. figure:: ./../data/Bf_TS.png
   :width: 500px
   
Apply the LM-FBS based on the defined coompatibility and equilibrium conditions.

.. code-block:: python

	Y_ABn = np.zeros_like(Y_AnB,dtype = complex)

	Y_int = Bu @ Y_AnB @ Bf.T
	Y_ABn = Y_AnB - Y_AnB @ Bf.T @ np.linalg.pinv(pyFBS.TSVD(Y_int,reduction = 22)) @ Bu @ Y_AnB
	
First extract the FRFs at the reference DoFs:

.. code-block:: python

	arr_out = [30,31,32,33,34,35]
	arr_in = [30,31,32,33,34,35]

	Y_AB_coupled = Y_ABn[:,arr_out,:][:,:,arr_in]
	Y_AB_ref = Y_AB_ref[:,:,:]
	
Finnaly, the coupled and the reference results can be compared and evaluated:

.. raw:: html

   <iframe src="../../_static/FRF_TS.html" height="500px" width="750px" frameborder="0"></iframe>

.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!#############################
Operational Deflection Shapes
#############################

The 3D display of the pyFBS can also be used to animate any objects. Animation can either be performed on meshes or on the predefined objects (such as accelerometers). 
In this example an Operational Deflection Shape (ODS) of an automotive testbench is animated.

.. note:: 
   Download example showing an application of the ODS: :download:`06_ODS.ipynb <../../../examples/06_operational_deflection_shapes_ODS.ipynb>`


Example Datasets and 3D view
----------------------------

As already shown in the `3D Display <../../html/examples/basic_examples/01_static_display.html>`_ one can load predefined datasets from the available examples and add a structure from STL file to the 3D view. This allows both the sensors and excitation points (impacts) to be visualized.

    
.. figure:: ./../data/six_one.png
   :width: 800px

Experimental example
********************
Load the experimental data, which will be used for the operational deflection shape animation

.. code-block:: python

	_file = r"./automotive_testbench/Measurements/ODS.p"
	freq, Y_ODS = np.load(_file,allow_pickle = True)
	
Checkout a single FRF:

.. code-block:: python

	select_out = 5
	select_in = 1

	Y_ODS = pyFBS.plot_frequency_response(freq, Y_ODS[:,select_out:select_out+1,select_in:select_in+1])

	
.. raw:: html

   <iframe src="../../_static/Y_ODS.html" height="500px" width="750px" frameborder="0"></iframe>
	
Accelerometer animation and GIF export
--------------------------------------
The objects placed in the 3D view can be simply animated.
In this example an operational deflection shape at a certain impact position is animated.
The pyFBS supports also an export to a GIF file. 
Before running the animation just set the output directory ``view3D.gif_dir`` and set the variable ``view3D.take_gif = True``. 
When the GIF is exporting the animation can lag within the 3D display. 

.. code-block:: python

	freq_sel = -1
	select_in = 6

	emp_2 = pyFBS.orient_in_global(Y_ODS[freq_sel,:,select_in],df_chn,df_acc)

	mode_dict = pyFBS.dict_animation(emp_2,"object",object_list = view3D.global_acc,r_scale=30)
	mode_dict["freq"] = freq[freq_sel]
	
	view3D.take_gif = True
	view3D.gif_dir = "..\\pyFBS\\output.gif"

	view3D.add_objects_animation(mode_dict,run_animation = True,add_note= True)
	

.. figure:: ./../data/ods.gif
   :width: 400px

.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!##########
3D Display
##########

The :mod:`pyFBS` can be used for 3D visualization. The 3D display enables depiction of structures, sensors, impacts, channels and virtual points in a simple and intuitive, Pythonic manner. 
Furthermore, the 3D display supports motion animation, where objects or mode shapes can be animated with ease. For the 3D visualization a python package `PyVista <https://docs.pyvista.org/index.html>`_ is used.

.. note:: 
   Download example showing the basic use of the 3D display: :download:`01_static_display.ipynb <../../../examples/01_static_display.ipynb>`


To open a blank 3D display simply make an instance of the :class:`pyFBS.view3D`.

.. code-block:: python

	view3D = pyFBS.view3D()

A rendering window will open in the background, which will not pause the code execution (for more details refer to the :class:`pyvista.BackgroundPlotter`). By default a coordinate system is placed in the origin and an orientation marker is placed in the bottom-left corner. 

*****************
Geometric objects
*****************
Geometric objects can be added to the 3D display in a simple manner. For simple objects (cylinders, spheres, boxes, ...) `PyVista methods <https://docs.pyvista.org/examples/00-load/create-geometric-objects.html#sphx-glr-examples-00-load-create-geometric-objects-py>`_ 
can be used for the geometry generation. For displaying a more complex geometric objects a STL file can be loaded in the 3D display (example datasets can be downloaded with :func:`pyFBS.download_lab_testbench`).

.. code-block:: python
 
	path_to_stl = r"./lab_testbench/STL/A.stl"
	view3D.add_stl(path_to_stl,name = "AB")

After the code execution the geometric object will apear in the rendering window.

.. figure:: ./../data/3D_view.png
   :width: 800px


Multiple geometric objects can be added to the display with different colors and even opacity (checkout :func:`pyvista.BackgroundPlotter.add_mesh` for all options). 
When adding multiple geometric objects, care should be taken that different ``name`` variable is provided, otherwise the object with the same name will be overwritten (discarded from the 3D display). 
   
***********
I/O Objects
***********
In the 3D display accelerometers, channels, impacts and virtual points can be shown. The positional and orientation information for each separate degree of freedom is defined in a :mod:`pandas.DataFrame`. 

Accelerometers
==============
Accelerometers can be added to 3D display directly from the :mod:`pd.DataFrame`.

.. code-block:: python

	path_to_xlsx = r"./lab_testbench/Measurements/AM_measurements.xlsx"
	
	df_acc = pd.read_excel(path_to_xlsx, sheetname='Sensors_AB')
	view3D.show_acc(df_acc)

After the code execution accelerometers will be shown in the 3D display.

.. figure:: ./../data/acc.png
   :width: 800px


Channels
========
Channels associated with accelerometers can be added to the 3D display directly from the :mod:`pd.DataFrame`.

.. code-block:: python

	df_chn = pd.read_excel(path_to_xlsx, sheetname='Channels_AB')
	view3D.show_chn(df_chn)

After the code execution channels will be shown in the 3D display.	

.. figure:: ./../data/chn.png
   :width: 800px


Impacts
=======
Impacts can be added to the 3D display directly from the :mod:`pd.DataFrame`.

.. code-block:: python

	df_imp = pd.read_excel(path_xlsx, sheetname='Impacts_AB')
	view3D.show_imp(df_imp)

After the code execution impacts will be shown in the 3D display.	


.. figure:: ./../data/imp.png
   :width: 800px



Virtual points
==============
Virtual points can also be added to the 3D display directly from the :mod:`pd.DataFrame`.

.. code-block:: python

	df_vps = pd.read_excel(path_xlsx, sheetname='VP_Channels')
	view3D.show_vp(df_vps)

After the code execution virtual points will be shown in the 3D display.	
	
.. figure:: ./../data/VP.png
   :width: 800px


Labels
======
Accelerometer, channels, impacts and virtual points can also be labeled or enumerated based on the information from the corresponding :mod:`pd.DataFrame`.

.. code-block:: python

	view3D.label_acc(df_vps)
	view3D.label_chn(df_chn)
	view3D.label_imp(df_imp)
	view3D.label_vp(df_vps)

Corresponding labels will appear in the 3D display after the code execution.

.. figure:: ./../data/labels.png
   :width: 800px

   
*******************************
Interaction with the 3D display
*******************************
Basic interaction with the 3D display is relatively simple. Mouse ``left-click`` can be used to rotate the rendering scene and ``middle-click`` to pan the rendering scene. 
For more information refer to the `PyVista plotting shortcuts <https://docs.pyvista.org/plotting/plotting.html>`_.

.. figure:: ./../data/interaction.gif
   :width: 800px

.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!=======================
Interactive positioning
=======================
With the :mod:`pyFBS` accelerometers, impacts and virtual points can be added and positioned interactively within the 3D display. 

.. note:: 
   Example showing the interactive positioning: :download:`02_interactive_display.ipynb <../../../examples/02_interactive_display.ipynb>`.

Dynamic interaction with objects in the :mod:`pyFBS` is achieved with the sphere widgets from PyVista. 
Altogether four sphere widgets are used (one for translation and three for rotations).
   
***********
Translation
***********
Translation of an object in the 3D display can be performed by moving a black sphere widget. An example of translation in the 3D display is depicted on a GIF bellow.

.. figure:: ./../data/translation.gif
   :width: 800px


********
Rotation
********
To rotate an object in the 3D display, three sphere widgets are available, for rotation around each axis. 
Red sphere widget is used for rotation around `X`-axis, green for rotation around `Y`-axis and blue for rotation around `Z`-axis. 
An example of rotation in the 3D display is depicted on a GIF bellow.

.. figure:: ./../data/rotation.gif
   :width: 800px


********
Snapping
********
If a mesh from an STL file is available there is possible to snap the position of the object to the surface of the geometric object.
Furthermore, the orientation of the object is also aligned with the surface normal at the snapping point. 
The snapping option can be disabled by holding down the letter ``T`` when moving the object in the 3D display.


.. figure:: ./../data/snapping.gif
   :width: 800px


***********
I/O Objects
***********
If a predefined dataset :mod:`pandas.DataFrame` is available for accelerometers, 
impacts and virtual points, it can be used to place interactive objects already on the predefined positions.

Accelerometers
==============   
To enable the snapping to mesh option, first load an STL file in the 3D display:

.. code-block:: python

	stl = r"./lab_testbench/STL/A.stl"
	mesh = view3D.add_stl(stl,name = "ts")

Accelerometers can then be placed on the predefined positions, which can then be moved around and rotated in the 3D display. 
If there is no predefined data, new accelerometers can be added to the 3D display by pressing the letter ``P``.
	
.. code-block:: python

	view3D.add_acc_dynamic(mesh,predefined = df_sensors)

The new updated positions and orientations can be obtained directly from the :class:`pyFBS.view3D`:

.. code-block:: python

	df_acc_updated = view3D.get_acc_data()


Channels
========   

Channels can be defined based on the positions and orientations of accelerometers.

.. code-block:: python

	df_chn_updated = pyFBS.utility.generate_channels_from_sensors(df_acc_updated)

Currently, all the accelerometers are considered to be tri-axial. However, possible redundant channels can simply be discarded from the ``df_chn_updated``.

Impacts
=======

Interactive impacts can be added also from the predefined positions.

.. code-block:: python

	view3D.add_imp_dynamic(mesh,predefined = df_impacts)
	
The updated positions and orientations can be obtained directly.	

.. code-block:: python

	df_imp_updated = view3D.get_imp_data()


Virtual points
==============

In a simmilar manner also interactive virtual points can be added to the 3D display.

.. code-block:: python

	view3D.add_vp_dynamic(mesh,predefined = df_vp)

The updated positions and orientations can be obtained directly.

.. code-block:: python

	df_vp_updated = view3D.get_vp_data()

Virtual point channels and forces can be defined based on the positions and orientations of the added VP. 
In this manner, 6-DoF VP for rigid IDMs is obtained. Additional DoFs can be added to the dataframe manually.

.. code-block:: python

	df_vp, df_vp_ref = pyFBS.generate_VP_from_position(df_vp_updated)


************************
Export updated positions
************************
The updated datasets can be exported to Excel file in a simple manner with the ``pd.ExcelWriter``.

.. code-block:: python
 
	with pd.ExcelWriter('./output_file.xlsx') as writer:  
		df_acc_updated.to_excel(writer, sheet_name='Sensors',index = False)
		df_imp_updated.to_excel(writer, sheet_name='Impacts',index = False)
		df_chn_updated.to_excel(writer, sheet_name='Channels',index = False)


.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!==================
FRF synthetization
==================
The :mod:`pyFBS` package enables user-friendly modal analysis and FRF synthetization based on the mass and stiffness matrices imported from the FEM software. 
Currently, only data import from Ansys is supported. 

.. note:: 
   Download example showing the basic use of the FRF synthetization: :download:`03_FRF_synthetization.ipynb <../../../examples/03_FRF_synthetization.ipynb>`

Numerical analysis of continuous systems requires their discretization by division into finite elements [1]_. 
The dynamic properties of the system are described by the equilibrium equation, 
where the external forces are equal to the internal forces resulting from inertia, damping and stiffness. 
The basic equation for a linear dynamical system with viscous damping has the following form:

.. math::

   \mathbf{M} \, \boldsymbol{\ddot{x}}(t) + \mathbf{C} \, \boldsymbol{\dot{x}}(t) + \mathbf{K} \, \boldsymbol{x}(t) = \boldsymbol{f}(t)

where :math:`\mathbf{M}` represents the mass matrix, :math:`\mathbf{C}` the damping matrix, :math:`\mathbf{K}` the stiffness matrix, 
:math:`\boldsymbol{x}` vector of responses at all :math:`n` degrees of freedom (DoFs), 
and :math:`\boldsymbol{f}` vector of externally applied forces. 

For consistent modeling of dynamic systems within the frequency domain the following assumptions have to be valid:

* Linearity - response amplitude is linearly proportional to the excitation amplitude.

* Time invariance - mass, stiffness, and damping characteristics are time-independent.

* Passivity - the energy flow in the system is always positive or equal to zero.

* Reciprocity - the response of the structure remains the same if the excitation and response location are switched.

* Stability - the response of the system is bounded if the excitation of the system is bounded.

If we apply the Fourier transform, a time function of response :math:`\boldsymbol{x}(t)` can be written in 
the frequency domain :math:`(\boldsymbol{x}(\omega))` and basic equation takes the following form: 

.. math::
    \mathbf{M} \, \boldsymbol{\ddot{x}}(\omega) + \mathbf{C} \, \boldsymbol{\dot{x}}(\omega) + \mathbf{K} \, \boldsymbol{x}(\omega) = \boldsymbol{f}(\omega).

By following the fact that :math:`\boldsymbol{\ddot{x}}(\omega) = - \omega^2 \boldsymbol{x}(\omega)` 
we can express :math:`\boldsymbol{x}(\omega)` on the right-hand side: 

.. math::
    (- \omega^2\,\mathbf{M} + \text{j}\, \omega\, \mathbf{C} + \mathbf{K}) \boldsymbol{x}(\omega) = \boldsymbol{f}(\omega).

The left-hand side can be collected into a single frequency-dependent matrix :math:`\mathbf{Z}(\omega)`:

.. math::
    \mathbf{Z}(\omega) \boldsymbol{x}(\omega) = \boldsymbol{f}(\omega).

The impedance matrix :math:`\mathbf{Z}(\omega)` is also called the dynamic stiffness matrix (often referred to as mechanical impedance or apparent mass). 
By inverting the mechanical impedance an admittance notation :math:`\mathbf{Y}(\omega)` can be obtained:

.. math::
    \boldsymbol{u}(\omega) = \mathbf{Y}(\omega) \boldsymbol{f}(\omega). \qquad \mathbf{Y}(\omega) = 
    \left(\mathbf{Z}(\omega)\right)^{-1} = \left(- \omega^2\,\mathbf{M} + \text{j}\, \omega\, \mathbf{C} + \mathbf{K} \right)^{-1}.

The last equation presents the direct harmonic method for FRF synthetization. The :math:`\mathbf{Y}(\omega)` denotes the frequency response function matrix and is often referred
to as admittance or dynamic flexibility. The admittance matrix is sometimes also denoted with the letter :math:`\mathbf{H}`. 

An admittance :math:`\mathbf{Y}_{ij}` is defined as the response at :math:`i`-th DoF when a unit excitation force is applied at :math:`j`-th DoF. 
A whole :math:`j`-th column from the admittance matrix :math:`\mathbf{Y}` can be determined from a single excitation point; 
therefore, the admittance matrix is fully coupled. A single admittance FRF is a global observation of the system dynamics.

.. figure:: ./../data/YY_matrix.svg
   :width: 250px
   :align: center

The impedance :math:`\mathbf{Z}_{ji}` is defined as the force at :math:`j`-th DoF when a unitary displacement is imposed at :math:`i`-th DoF, 
while all remaining DoFs are fixed. 
Therefore, the impedance matrix :math:`\mathbf{Z}` is commonly sparse in the sense that a single impedance term is a local observation of the system dynamics.

.. figure:: ./../data/ZZ_matrix.svg
   :width: 250px
   :align: center

Eigenvalue problem
******************

To determine the eigenfrequencies and eigenvectors, the considered system responds with free vibration while the damping is neglected. 
The equilibrium equation takes the form of a homogeneous second-order differential equation:

.. math::

    \mathbf{M}\,\boldsymbol{\ddot{x}}(t)+\mathbf{K}\,\boldsymbol{x}(t)=\boldsymbol{0}

Euler's identity represents the solution of the differential equation: :math:`\boldsymbol{x}(t)=\boldsymbol{X}\,e^{\text{i}\,\omega\,t}` and
:math:`\boldsymbol{\ddot{x}}(t)=-\omega^2\boldsymbol{X}\,e^{\text{i}\,\omega\,t}`.
By transformation to the modal domain the equation of motion takes the form 
:math:`\mathbf{M}(-\omega^2\boldsymbol{X}\,e^{\text{i}\,\omega\,t})+\mathbf{K}\,(\boldsymbol{X}\,e^{\text{i}\,\omega\,t})=\boldsymbol{0}`.
By knowing :math:`e^{\text{i}\,\omega\,t}\neq 0` for any time :math:`t` we obtain: 

.. math::

   (\mathbf{K} - \omega^2\,\mathbf{M}) \boldsymbol{X} = \boldsymbol{0}

To get non-trivial solution it is necessary to satisfy :math:`\text{det}(\mathbf{K}-\omega_r^2\,\mathbf{M})=\boldsymbol{0}`.
By solving the determinant, the eigenvalues :math:`\omega_1^2, \omega_2^2, \dots` are determined which represent undamped natural frequencies. 
Each eigenvalue result in corresponding eigenvector :math:`\boldsymbol{\psi}_1, \boldsymbol{\psi}_2, \dots` which prepresents mode shape.
Eigenvalues and eigenvectors can be organised into a matrix form:

.. math::

   [^{\nwarrow}{\pmb{\omega}_{r}^2}_{\searrow}]=\begin{bmatrix} 
   \omega_1^2 & 0 & \cdots & 0 \\
   0 & \omega_2^2 & \cdots & 0 \\
   \vdots & \vdots & \ddots & \vdots \\
   0 & 0 & \cdots & \omega_N^2
   \end{bmatrix}, 
   \qquad
   [\pmb{\Psi}] = \begin{bmatrix} \boldsymbol{\psi}_1 & \boldsymbol{\psi}_2 & \cdots & \boldsymbol{\psi}_N\end{bmatrix}

Modal mass and modal stifness are calculated using following equiations:

.. math::

   \pmb{\Psi}^{\text{T}}\,\mathbf{M}\,\pmb{\Psi} = [^{\nwarrow}{m_{r}}_{\searrow}], 
   \qquad
   \pmb{\Psi}^{\text{T}}\,\mathbf{K}\,\pmb{\Psi} = [^{\nwarrow}{k_{r}}_{\searrow}]

Mass normalised mode is calculated using equation: :math:`\boldsymbol{\phi}_r=\boldsymbol{\psi}_r\,\frac{1}{\sqrt{m_{r}}}`
and has following properies:

.. math::

   \pmb{\Phi}^{\text{T}}\,\mathbf{M}\,\pmb{\Phi} = [\mathbf{I}], 
   \qquad
   \pmb{\Phi}^{\text{T}}\,\mathbf{K}\,\pmb{\Phi} = [^{\nwarrow}{\pmb{\omega}_{r}^2}_{\searrow}]

For FRF generation mode superposition method can be used, where contributions of modes are superimposed at each frequency line using equation: 

.. math::

   \alpha_{i, j} = \sum_{r=1}^{m}\frac{\boldsymbol{\phi}_{i,r}\,\boldsymbol{\phi}_{j,r}}{\omega_r^2-\omega+i\eta_r\,\omega_r^2}

where :math:`\eta_r` represents modal damping at :math:`r`-th natural frequency and can be neglected for lightly damped systems. 
Index :math:`i` represents the location of response and index :math:`j` stands for the location of excitation.
The number of modes used for reconstruction is equal to :math:`m` and is usually much lower than the number of DoFs (:math:`m \ll n`). 
Therefore, modal truncation occurs.

Using the :mod:`pyFBS` package, it is easy to calculate eigenfrequencies and modal shapes from an imported finite element model 
and visualize them using an animated 3D display.
Also, the FRF synthetization is user friendly and is supported with the mode superposition and direct harmonic method.

MK model initialization
***********************

First, the MK model is initialized (with class :class:`pyFBS.MK_model`) by importing ``.rst`` and ``.full`` files, 
which contain the information on the locations of finite element nodes, their DoFs, the connection between the nodes, 
the mass and stiffness matrix of the system.

.. code-block:: python

    import pyFBS
    from pyFBS.utility import *

    full_file = r"./lab_testbench/FEM/B.full"
    rst_file = r"./lab_testbench/FEM/B.rst"

    MK = pyFBS.MK_model(rst_file, full_file, no_modes = 100, allow_pickle = False, recalculate = False)

In this step also the eigenfrequencies and eigenvectors of the system are simultaneously calculated. Eigenvectors are mass normalised. The number of calculated eigenvalues is limited by the ``no_modes`` parameter.

.. tip::

   In the case of models with a huge number of DoFs, the process of solving the eigenproblem can take quite some time (depending on the complexity of the model and the computational power of the computer).
   By setting ``read_rst = True`` in :class:`pyFBS.MK_model` initialization modal parameters will be imported directly from the ``.rst`` file and not calculated again inside Python.

.. warning::
    The current version uses pickle module to store the modal parameters. The pickle module is not secure. Only unpickle data you trust.
	
.. tip::
    The imported model will only contain as many eigenvalues and eigenforms as there were calculated in Ansys. 
    If you want to use more, you need to re-solve the problem in Python (by setting: ``allow_pickle = True, read_rst == False``) or Ansys.



Mode shape visualization
************************

After the MK model is defined, the calculated mode shapes can be animated. For nicer representation, the STL file can be added to visualize the undeformed shape of the structure.
The 3D display is opened in a new window and allows the user to interact with the added model.

.. code-block:: python

    stl = r"./lab_testbench/STL/B.stl"
    view3D = pyFBS.view3D(show_origin= True)
    view3D.add_stl(stl,name = "engine_mount",color = "#8FB1CC",opacity = .1)   

To animate mode shape, a mesh of finite elements must be added to display. 
Colormap of the model can be changed using ``cmap`` parameter, which supports all `PyVista colormap choices <https://docs.pyvista.org/examples/02-plot/cmap.html>`_.

.. code-block:: python

    view3D.plot.add_mesh(MK.mesh, scalars = np.ones(MK.mesh.points.shape[0]), cmap = "coolwarm", show_edges = True)

Mode shape can be selected with the method ``get_modeshape``. 
Animation parameters are defined with the function ``dict_animation``, which is imported from :mod:`pyFBS.utility`. 
Here you can set the frame rate (``fps``), the relative scale of deformation (``r_scale``) and the number of points in the animation sequence (``no_points``).

.. code-block:: python

    select_mode = 6
    _modeshape = MK.get_modeshape(select_mode)

    mode_dict = dict_animation(_modeshape,"modeshape",pts = MK.pts, mesh = MK.mesh, fps=30, r_scale=10, no_points=60)
    view3D.add_modeshape(mode_dict,run_animation = True)

Animation is visible in the previously defined pop-up window. The following figure shows animated 7th mode shape. 

.. figure:: ./../data/mode_shape_animation3.gif
   :width: 800px
   
To show undeformed mesh you can simply click the button in the pop-up window or call a method :func:`clear_modeshape()`:

.. code-block:: python

    view3D.clear_modeshape()

Visualization of impacts and responses
======================================

Locations and directions of impacts and responses must be passed with a :mod:`pd.DataFrame`. 
They can be either read from an Excel file (as shown below) or generated directly with pyFBS (see `Interactive display example <https://pyfbs.readthedocs.io/en/latest/examples/02_interactive_display.html>`_).
The parameter ``df_acc`` must include the following columns header: ``Position_1``, ``Position_2``,  ``Position_3``, ``Orientation_1``, ``Orientation_2``, ``Orientation_3``. 
The position parameters describe the location of accelerometers in the global coordinate system. 
Orientation of coordinate systems of accelerometer regarding the global coordinate system is defined with orientation parameters which are defined with Euler angles in degrees.
The parameter ``df_chn`` and ``df_imp`` must include the following columns header: ``Position_1``, ``Position_2``,  ``Position_3``, ``Direction_1``, ``Direction_2``, ``Direction_3``. 
The directions presents unit vector directions. 

.. code-block:: python

    # Path to .xslx file
    xlsx = r"./lab_testbench/Measurements/AM_measurements.xlsx"

    # Import and show locations of accelereometers
    df_acc = pd.read_excel(xlsx, sheet_name='Sensors_B')
    view3D.show_acc(df_acc,overwrite = True)

    # Import and show directions of accelereometers channels
    df_chn = pd.read_excel(xlsx, sheet_name='Channels_B')
    view3D.show_chn(df_chn)

    # Import and show locations and directions of impacts
    df_imp = pd.read_excel(xlsx, sheet_name='Impacts_B')
    view3D.show_imp(df_imp,overwrite = True)

.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
            <iframe src="https://kitware.github.io/vtk-js/examples/SceneExplorer/index.html?fileURL=https://dl.dropbox.com/s/n5nhk4f9wsd8l9t/FRF_synthetization_imp_chn.vtkjs?dl=0" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>

Defining DoFs of synthetized FRFs
*********************************

FRFs can currently only be synthetized on the nodes of the numerical model. 
Therefore, it is necessary to find the nodes closest to the desired locations in the numerical model and update them. 
The orientation of the generated FRFs is independent of the direction in the numerical model and will not change with the updated location.

Locations of impacts and channels can be updated to the nodes of the numerical model with the :func:`pyFBS.MK_model.update_locations_df`:

.. code-block:: python

    df_chn_up = MK.update_locations_df(df_chn)
    df_imp_up = MK.update_locations_df(df_imp)

Updated locations can also be displayed in the 3D display. 
By setting ``overwrite`` to ``False`` added locations of impacts and responses won't override previously added features, thus everything is shown simultaneously.

.. code-block:: python

    view3D.show_chn(df_chn_up, color = "y", overwrite = False)
    view3D.show_imp(df_imp_up, color = "y", overwrite = False)

.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
            <iframe src="https://kitware.github.io/vtk-js/examples/SceneExplorer/index.html?fileURL=https://dl.dropbox.com/s/n5nhk4f9wsd8l9t/FRF_synthetization_imp_chn.vtkjs?dl=0" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>

FRF synthetization
******************

FRFs for relativelly small systems can be efficiently computed using the full harmonic method:

.. code-block:: python

    MK.FRF_synth_full(f_start = 0, f_end = 2000, f_resolution = 1, frf_type = "accelerance")

.. warning::

    Using this formulation, a full admitance matrix is generated. 
    For large systems, this can be very consimung in terms of computational times and memory storage.

Another possibility to generate FRFs is by applying the mode superposition method.
FRFs are synthetized at given locations and directions in ``df_channel`` and ``df_impact`` parameters. 
Even if we forget to define updated response and excitation locations, the function will automatically find 
the nearest nodes in the numerical model from which the FRFs are generated. 
Frequency properties are defined in parameters ``f_start``, ``f_end`` and ``f_resolution``. 
The number of modes used for modal superposition FRF generation is defined in the ``no_modes`` parameter 
and coefficient of modal damping is defined in parameter ``modal_damping``.
The resulting FRFs can be in the form of ``accelerance``, ``mobility`` or ``receptance``, 
which is defined in the ``frf_type`` parameter.

.. code-block:: python

    MK.FRF_synth(df_channel = df_chn, df_impact = df_imp, 
                 f_start = 0, f_end = 2000, f_resolution = 1, 
                 limit_modes = 50, modal_damping = 0.003, 
                 frf_type = "accelerance")

.. tip::
    Compared to the direct harmonic, this method can be applied only to a reduced (arbitratily selected) set of DoFs and is 
    therefore preferable in terms of computational and memory cost. Note, that typically also modal truncation is applied in the process.
    Nevertheless, even if you are interested in a small frequency region only, you should always include eigenvalues from 
    outside this region to account for upper- and lower-residuals.

The DoFs in the FRF matrix row follows the order of responses in the ``df_channel`` parameter, 
and the DoFs column matches the order of excitations in ``df_impact``.

Adding noise
============

To analyze various real-life experiments, numerically obtained FRFs are often intentionally contaminated with random noise to follow experimental data. 
Noise can be added to FRFs by the ``add_noise`` method.

.. code-block:: python

    MK.add_noise(n1 = 2e-1, n2 = 2e-1, n3 = 5e-2 ,n4 = 5e-2)

FRF visualization
=================

An experimental measurement is imported to compare all FRFs.

.. code-block:: python

    exp_file = pyFBS.example_lab_testbench["meas"]["Y_B"]

    freq, Y_B_exp = np.load(exp_file,allow_pickle = True)

Comparison of different FRFs can be performed visually:

.. code-block:: python

    o = 3
    i = 0

    pyFBS.plot_frequency_response(freq, 
                np.hstack((MK.FRF_noise[:,o:o+1,i:i+1], MK.FRF[:,o:o+1,i:i+1], Y_B_exp[:,o:o+1,i:i+1])),
                labels=('Num. FRF + noise', 'Num. FRF', 'Exp. FRF'))
	
.. raw:: html

   <iframe src="../../_static/FRF_synth.html" height="460px" width="100%" frameborder="0"></iframe>




.. panels::
    :column: col-12 p-3

    **That's a wrap!**
    ^^^^^^^^^^^^

    Want to know more, see a potential application? Contact us at info.pyfbs@gmail.com!

.. rubric:: References

.. [1] e Silva, Júlio M. Montalvão, and Nuno MM Maia, eds. Modal analysis and testing. Vol. 363. Springer Science & Business Media, 2012.