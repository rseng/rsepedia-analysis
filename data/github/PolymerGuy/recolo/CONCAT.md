---
title: 'RECOLO: A Python package for the reconstruction of surface pressure loads from kinematic fields using the virtual fields method'
tags:
  - Python
  - Virtual fields method
  - Load reconstruction
  - Parameter identification
authors:
  - name: Sindre Nordmark Olufsen
    affiliation: "1, 2"
  - name: Rene Kaufmann
    affiliation: 1
  - name: Egil Fagerholt
    affiliation: 1
  - name: Vegard Aune
    affiliation: "1, 2"
affiliations:
 - name: Structural Impact Laboratory (SIMLab), Department of Structural Engineering, NTNU - Norwegian University of Science and Technology, Trondheim, Norway
   index: 1
 - name: Centre for Advanced Structural Analysis (CASA), Department of Structural Engineering, NTNU - Norwegian University of Science and Technology, Trondheim, Norway
   index: 2
date: 30 October 2021
bibliography: paper.bib
---

# Summary
In experimental mechanics,  conducting non-intrusive measurements of surface pressure distributions acting on blast-loaded structures remains a challenge even in controlled, laboratory environments (see e.g., [@Pannell2021]). Still, for the design of tomorrow's sustainable and material-efficient structures, detailed knowledge of how pressure loads from extreme loading events interact with deformable structures is essential. When pressure loads are imposed on a deformable structure, fluid-structure interaction (FSI) effects are known to cause non-trivial loading scenarios which are difficult to quantify (see e.g., [@Aune2021]).
This project aims at reconstructing the full-field surface pressure loads acting on a deforming structure employing the virtual fields method (VFM) on full-field kinematic measurements [@Kaufmann2019].
Even though the current framework is limited to reconstructions of full-field pressure information from deformation data of thin plates in pure bending, it also allows for future extensions to other loading and deformation scenarios.
Provided that the properties of the structure are known,
the transient pressure distribution on the plate can be reconstructed. To understand the capabilities and accuracy
associated with the pressure reconstruction methodology, the package provides tools for performing virtual experiments based on analytical data or data from finite element simulations. The current implementation is based on the deflectometry technique, using the grid method to obtain the deformation measurements and corresponding kinematics of the structure.

This Python package is made for RECOnstructing surface pressure LOads, ``RECOLO``, acting on plated structures based on deformation measurements using the VFM [@Pierron2012].
The current implementation determines the surface pressure acting on a thin plate undergoing small deformations, assuming linear, elastic material behaviour. However, the framework will be extended to large plastic deformations, allowing the two-way interaction between the pressure loading and the deformation of the plate to be studied.
Other VFM toolkits such as PeriPyVFM are readily available but typically aimed at determining material properties from deformation and load measurements. Hence, as opposed to other VFM toolkits, ``RECOLO`` assumes that the material properties are known and use the full-field deformation measurements to reconstruct the pressure loading.

``RECOLO`` contains a collection of tools enabling the user to perform virtual experiments on synthetically generated data as well
 as performing pressure reconstruction on experimental datasets. The pressure reconstruction algorithm is based on the work by [@Kaufmann2019].
The implementation is based on numerical operations provided by NumPy [@Numpy] and SciPy [@SciPy] as well as visualization by Matplotlib [@Matplotlib].



# Statement of need
``RECOLO`` was established to quantify the blast loading acting on plated structures in a purpose-built shock tube apparatus at SIMLab, NTNU [@Aune2016]. To the authors' best knowledge there a no open-source software providing the functionality necessary to perform pressure load reconstruction based on the kinematics of plated structures during fast transient dynamics, motivating the ``RECOLO`` project.

The methodology developed in this project is directly applicable to obtain new, unique insight into surface pressure distributions on plated structures subjected to blast loading. This project is part of the ongoing research within the SIMLab research group at NTNU.

# Acknowledgements
The authors gratefully appreciate the financial support from the Research Council of Norway (RCN) through the Centre for Advanced Structural Analysis (SFI-CASA RCN Project No. 237885) and the SLADE KPN project (RCN Project No. 294748). The financial support by the Norwegian Ministry of Justice and Public Security is also greatly appreciated.


# References![](docs/logo.png)
=============
[![codecov](https://codecov.io/gh/PolymerGuy/recolo/branch/master/graph/badge.svg?token=7J4EH3C399)](https://codecov.io/gh/PolymerGuy/recolo)
[![CircleCI](https://circleci.com/gh/PolymerGuy/recolo.svg?style=svg&circle-token=3403eba7b905e1a626d1c797ed5ca4e3daba76df)](https://circleci.com/gh/PolymerGuy/recolo)
[![MIT License][license-shield]][license-url]
[![Documentation Status](https://readthedocs.org/projects/recolo/badge/?version=latest)](https://recolo.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03980/status.svg)](https://doi.org/10.21105/joss.03980)




About this project
------------------
In experimental mechanics, measuring the pressure load acting on a surface in a non-intrusive manner is of high interest for several applications. However, techniques allowing for such non-intrusive measurements have not been freely available to the experimental mechanics community.

**RECOLO** provides tools for reconstructing distributed pressure loads acting on thin, elastic plates. The Python package implements the virtual fields method (VFM), where the input is the kinematic fields governing the plate dynamics.

A virtual lab is also provided, allowing synthetic data to be generated based on input from finite element simulations.

Example kinematic fields pressure is shown below:
![alt text](docs/kinematics.gif)

which gives the following pressure field:
![alt text](docs/pressure.gif)

The documentation is hosted at https://recolo.readthedocs.io/en/latest/


Getting Started
---------------
Clone the repo by using `git`:

```bash
git clone https://github.com/PolymerGuy/recolo.git
```

when in the folder with the repo, make a virtual environment and install all dependencies:

```bash
# Make virtual environment using venv
python -m venv env
# Activate the virtual environment
source ./env/bin/activate
# Install dependencies
pip install -r requirements.txt
```

To check that everything is working, run all tests:
```bash
python -m pytest --pyargs recolo --cov=./
```


Building the documentation from source
--------------------------------------
```bash
# Enter the documentation folder
cd docs
# Rebuild docs
make html
```

The documentation is now found in ./_build_/html


How to contribute
-----------------
The RECOLO project welcomes your expertise and enthusiasm!

You are more than welcome to contribute to this project, e.g., by:
* Report a bug using the Github issue tracker.
* Fix an already reported bug via a pull-request.
* Add new functionality via a pull-request.
* Help revise pull-requests.

When you wish to submit new or revised code, we encourage the following procedure:

* Fork the repo
* Implement your changes/additions
* Make a pull-request to the -dev branch

How to cite us
--------------
If you use this toolkit in your scientific work, consider citing one or more of the following:

- S. N. Olufsen, R. Kaufmann, E. Fagerholt, V. Aune (2022). RECOLO: A Python package for the reconstruction of surface pressure loads from kinematic fields using the virtual fields method. Journal of Open Source Software, 7(71), 3980, https://doi.org/10.21105/joss.03980


[license-shield]: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square
[license-url]: https://choosealicense.com/licenses/mit
Overview
========
**RECOLO** provides tools for reconstructing the pressure field acting on a plate based on
the corresponding deformation measurements. The tools allow the user to perform virtual experiments on synthetically generated data as well as performing pressure reconstruction on experimental datasets of full-field deformation measurements (See figure below).

.. image:: ./figures/recolo.png
   :scale: 100 %
   :alt: Illustration of RECOLO
   :align: center

Currently, this tool is limited to small elastic deformations, but will be extended to
large plastic deformations in the future. In order to assess the capabilities of this toolkit, a suite of tools
for performing virtual experiments are provided. The tools allow the user to generate synthetic input for the load
reconstruction based on the results from finite element simulations, providing strict control over all aspects of the imposed loading
(spatial and temporal variations) as well as experimental aspects such as signal-to-noise ratio. The synthetic data is also a powerful tool for testing of the toolkit, facilitating adoption, development and validation.

The toolkit can be seen as two complementary parts with the following functionality:

* Tools for determining the dynamic pressure field acting on a thin plate based on deformation measurements
    * In-plane deformation measurements using "the grid method"
    * Slope and deflection measurements using deflectometry
    * Controlling the signal-to-noise ratio using a variety of filters
    * Pressure reconstruction using the virtual fields method
* Tools for generating synthetic input based on dynamic finite element simulations of plates subjected to arbitrary pressure fields.
    * Parser for the finite element suite Abaqus
    * Artificial grid deformation for generating images for deflectometry

Getting started
---------------
To install the toolkit, please see :doc:`install`

For a quick start showing just the pressure reconstruction, see :doc:`quickstart`

For a full examples including generation of synthetic data, please see :doc:`deflectometryAbaqus`.

For a full examples on an experimental dataset, please see :doc:`impactHammer`.

For examples showing how to set up other synthetic experiments with various level of detail,
please see the examples available on Github  https://github.com/PolymerGuy/recolo/tree/master/examples

The individual components of the toolkit can be browsed in the API-documentation.Experiment with impact hammer
=============================

This example presents the force reconstruction from an experiment with an instrumented hammer impacting a steel plate.
The impacting force from the hammer is reconstructed from the kinematic fields measured by deflectometry.

Due to the large size of the dataset, the data is hosted on https://dataverse.no/ and can either be downloaded manually
or using this toolkit.

Let's now go through the necessary steps for force reconstruction based on deflectometry.
First, we need to import the Recolo toolkit::

    import recolo

The experimental data is conveniently downloaded and accessed via the ImpactHammerExperiment class::

    exp_data = recolo.demoData.ImpactHammer()

After the download has completed, the force measurements can be accessed as::

    hammer_force, hammer_time = exp_data.hammer_data()

The experiment was performed on a 300 mm x 300 mm rectangular plate with the following properties::

     mat_E = 190.e9  # Young's modulus [Pa]
     mat_nu = 0.3  # Poisson's ratio []
     density = 7934
     plate_thick = 4.95e-3

The stiffness of the plate is calculated as::

     plate = recolo.make_plate(mat_E, mat_nu, density, plate_thick)

The experimental configuration is given by::

    grid_pitch = 7.0  # pixels
    grid_pitch_len = 2.5 / 1000.  # m
    mirror_grid_distance = 1.63  # m
    pixel_size_on_mirror = grid_pitch_len / grid_pitch * 0.5

We now set the pressure reconstruction window size::

     win_size = 30

as well as filter and downsampling settings::

     downsampling_factor = 5
     filter_time_sigma = 6
     filter_space_sigma = 2

The grid images are used as input to deflectomerty and the slope fields of the plate are determined::

    slopes_y, slopes_x = recolo.deflectomerty.slopes_from_images(exp_data.path_to_data_folder, grid_pitch,
                                                             mirror_grid_distance, pixel_size_on_grid_plane,
                                                             ref_img_ids=ref_img_ids,
                                                             only_img_ids=use_imgs,
                                                             crop=(76, -35, 10, -10), window="triangular",
                                                             correct_phase=False)

The slope fields are then integrated to determine the deflection fields::

     # Integrate slopes to get deflection fields
     disp_fields = recolo.slope_integration.disp_from_slopes(slopes_x, slopes_y, pixel_size,
                                                            zero_at="bottom corners", zero_at_size=5,
                                                            extrapolate_edge=0, downsample=downsampling_factor)

Based on these fields, the kinematic fields (slopes and curvatures) are calculated::

     kin_fields = recolo.kinematic_fields_from_deflections(disp_fields, downsampling_factor * pixel_size_on_mirror, sampling_rate,
                                                       filter_time_sigma=filter_time_sigma)


Now, the pressure reconstuction can be initiated. First we define the Hermite16 virtual fields::

     virtual_field = recolo.virtual_fields.Hermite16(win_size, pixel_size_on_mirror)

and initialize the pressure reconstruction::

     times = []
     presses = []

     for i, field in enumerate(kin_fields):
          recon_press = recolo.solver_VFM.calc_pressure_thin_elastic_plate(field, plate, virtual_field)
          presses.append(recon_press)
          times.append(field.time)

     presses = np.array(presses)

The results can then be visualized::

     center = int(presses.shape[1] / 2)

     # Plot the results
     plt.figure(figsize=(7,5))
     plt.plot(times, np.sum(presses, axis=(1, 2)) * ((pixel_size_on_mirror * downsampling_factor) ** 2.), label="VFM force from whole plate")
     plt.plot(times, np.sum(presses[:,20:50,20:50], axis=(1, 2)) * ((pixel_size_on_mirror * downsampling_factor) ** 2.), label="VFM force from subsection of plate")
     plt.plot(hammer_time, hammer_force, label="Impact hammer")
     plt.xlim(left=0.0008, right=0.003)
     plt.ylim(top=500, bottom=-100)
     plt.xlabel("Time [ms]")
     plt.ylabel(r"Force [N]")

     plt.legend(frameon=False)
     plt.tight_layout()
     plt.show()


The resulting plot looks like this:

.. image:: ./figures/hammer_force.png
   :scale: 80 %
   :alt: The results
   :align: center

A few things should be noted:
     * The force level is highly sensitive to the area over which the pressure is integrated.
     * The deviations are believed to be caused by interaction between deflection at the position of the hammer and the boundary conditions.
     * Filtering influences the force amplitude, but relatively large filter kernels can be used without decreasing the force level.

Abaqus experiment with grid deformation and deflectometry
=========================================================

Let's now go through the necessary steps for doing pressure reconstruction.
First, we need to import the tools::

     import recolo
     import numpy as np

The example data can be downloaded from the recolo/examples/AbaqusExamples/AbaqusRPTs folder.
The dataset corresponds to a 300 x 300 mm  thin plate exposed to a sinusoidal pressure distribution in space and a saw-tooth shaped history in time.
::

     mat_E = 210.e9  # Young's modulus [Pa]
     mat_nu = 0.33  # Poisson's ratio []
     density = 7700
     plate_thick = 5e-3
     plate = recolo.make_plate(mat_E, mat_nu, density, plate_thick)
     

We now set the pressure reconstuction window size. 
Note that as we here use noise free data on a relatively coarse dataset, a very small window size is used::

     win_size = 6

Additive gaussian noise is added to the deformed grid images to mimic real world noise levels::
     
     noise_std = 0.009

We now load Abaqus data::

     abq_sim_fields = recolo.load_abaqus_rpts("path_to_abaqus_data"))


In this case, the deflection fields from Abaqus are used to generate grid images with the corresponding distortion.
The grid images are then used as input to deflectomerty and the slope fields of the plate are determined::

     slopes_x = []
     slopes_y = []
     undeformed_grid = recolo.artificial_grid_deformation.deform_grid_from_deflection(abq_sim_fields.disp_fields[0, :, :],
                                                                                     abq_sim_fields.pixel_size_x,
                                                                                     mirror_grid_dist,
                                                                                     grid_pitch,
                                                                                     img_upscale=upscale,
                                                                                     img_noise_std=noise_std)
     for disp_field in abq_sim_fields.disp_fields:
          deformed_grid = recolo.artificial_grid_deformation.deform_grid_from_deflection(disp_field,
                                                                                          abq_sim_fields.pixel_size_x,
                                                                                          mirror_grid_dist,
                                                                                          grid_pitch,
                                                                                          img_upscale=upscale,
                                                                                          img_noise_std=noise_std)

          disp_x, disp_y = recolo.deflectomerty.disp_from_grids(undeformed_grid, deformed_grid, grid_pitch)
          slope_x = recolo.deflectomerty.angle_from_disp(disp_x, mirror_grid_dist)
          slope_y = recolo.deflectomerty.angle_from_disp(disp_y, mirror_grid_dist)
          slopes_x.append(slope_x)
          slopes_y.append(slope_y)

     slopes_x = np.array(slopes_x)
     slopes_y = np.array(slopes_y)
     pixel_size = abq_sim_fields.pixel_size_x / upscale

The slope fields are then integrated to determine the defletion fields::

     # Integrate slopes to get deflection fields
     disp_fields = recolo.slope_integration.disp_from_slopes(slopes_x, slopes_y, pixel_size,
                                                            zero_at="bottom corners", zero_at_size=5,
                                                            extrapolate_edge=0, downsample=1)
     
Based on these fields, the kinematic fields (slopes and curvatures) are calculated. 
::

     kin_fields = recolo.kinematic_fields_from_deflections(disp_fields, pixel_size,
                                                     abq_sim_fields.sampling_rate,filter_space_sigma=10)

Now, the pressure reconstuction can be initiated. First we define the Hermite16 virtual fields::

     virtual_field = recolo.virtual_fields.Hermite16(win_size, abq_sim_fields.pixel_size_x)

and initialize the pressure reconstruction::

     pressure_fields = np.array(
     [recolo.solver_VFM.pressure_elastic_thin_plate(field, plate, virtual_field)
                                                      for field in kin_fields])


The results can then be visualized::

     import matplotlib.pyplot as plt
     # Plot the correct pressure in the center of the plate
     times = np.array([0.0, 0.00005, 0.00010, 0.0003, 0.001]) * 1000
     pressures = np.array([0.0, 0.0, 1.0, 0.0, 0.0]) * 1e5
     plt.plot(times, pressures, '-', label="Correct pressure")

     # Plot the coreconstructed pressure in the center of the plate
     center = int(pressure_fields.shape[1] / 2)
     plt.plot(abq_sim_fields.times * 1000., pressure_fields[:, center, center], "-o",label="Reconstructed pressure")

     plt.xlim(left=0.000, right=0.3)
     plt.ylim(top=110000, bottom=-15)
     plt.xlabel("Time [ms]")
     plt.ylabel(r"Overpressure [kPa]")

     plt.legend(frameon=False)
     plt.tight_layout()
     plt.show()

The resulting plot looks like this:

.. image:: ./figures/minimalExamplePressure.png
   :scale: 80 %
   :alt: The results
   :align: center

Minimal example
===============

Let's now go through the necessary steps for doing pressure reconstruction.
First, we need to import the tools::

     import recon
     import numpy as np

The example data can be downloaded from the recon/examples/AbaqusExamples/AbaqusRPTs folder. 
The dataset corresponds to a 300 x 300 mm  thin plate exposed to a sinusoidal pressure distribution in space and a saw-tooth shaped history in time.
::

     mat_E = 210.e9  # Young's modulus [Pa]
     mat_nu = 0.33  # Poisson's ratio []
     density = 7700
     plate_thick = 5e-3
     plate = recon.make_plate(mat_E, mat_nu, density, plate_thick)
     

We now set the pressure reconstuction window size. 
Note that as we here use noise free data on a relatively coarse dataset, a very small window size is used::

     win_size = 6

We now load Abaqus data::

     abq_sim_fields = recon.load_abaqus_rpts("path_to_abaqus_data"))

The Abaqus data contains the out-of-plane deflection and acceleration at each node of the plate for every time step.
Based on these fields, the kinematic fields (slopes and curvatures) are calculated. By default the accelerations are determined from the 
deflection fields directly, but we here choose to use the acceleration fields from Abaqus. This is done by setting the keyword argument "acceleration_field".
::

     kin_fields = recon.kinematic_fields_from_deflections(abq_sim_fields.disp_fields, 
                                                            abq_sim_fields.pixel_size_x,
                                                            abq_sim_fields.sampling_rate,
                                                            acceleration_field=abq_sim_fields.accel_fields)

Now, the pressure reconstuction can be initiated. First we define the Hermite16 virtual fields::

     virtual_field = recon.virtual_fields.Hermite16(win_size, abq_sim_fields.pixel_size_x)

and initialize the pressure reconstruction::

     pressure_fields = np.array(
     [recon.solver_VFM.pressure_elastic_thin_plate(field, plate, virtual_field) 
                                                      for field in kin_fields])


The results can then be visualized::

     import matplotlib.pyplot as plt
     # Plot the correct pressure in the center of the plate
     times = np.array([0.0, 0.00005, 0.00010, 0.0003, 0.001]) * 1000
     pressures = np.array([0.0, 0.0, 1.0, 0.0, 0.0]) * 1e5
     plt.plot(times, pressures, '-', label="Correct pressure")

     # Plot the coreconstructed pressure in the center of the plate
     center = int(pressure_fields.shape[1] / 2)
     plt.plot(abq_sim_fields.times * 1000., pressure_fields[:, center, center], "-o",label="Reconstructed pressure")

     plt.xlim(left=0.000, right=0.3)
     plt.ylim(top=110000, bottom=-15)
     plt.xlabel("Time [ms]")
     plt.ylabel(r"Overpressure [kPa]")

     plt.legend(frameon=False)
     plt.tight_layout()
     plt.show()

The resulting plot looks like this:

.. image:: ./figures/minimalExamplePressure.png
   :scale: 80 %
   :alt: The results
   :align: center
.. Recolo documentation master file, created by
   sphinx-quickstart on Mon May  3 13:35:33 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Recolo's documentation!
==================================

This python package provides tools for reconstuction of pressure loads using the Virtual Fields Method (VFM).




.. toctree::
   :maxdepth: 2
   :caption: Getting started:

   overview
   install
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Virtual experiments:

   deflectometryAbaqus

.. toctree::
   :maxdepth: 2
   :caption: Real experiments:

   impactHammer

.. toctree::
   :maxdepth: 2
   :caption: Theory:

   VFM

.. autosummary::
   :toctree: _autosummary
   :caption: API Documentation:
   :template: custom-module-template.rst
   :recursive:

   recolo
   

Citing us:
----------
If you use this toolkit as part of your research, please cite the following:


S. N. Olufsen, R. Kaufmann, E. Fagerholt, V. Aune (2022). RECOLO: A Python package for the reconstruction of surface pressure loads from kinematic fields using the virtual fields method. Journal of Open Source Software, 7(71), 3980, https://doi.org/10.21105/joss.03980
VFM
===
Pressure Reconstruction
-----------------------
The Virtual fields method can be applied to solve a broad range of problems in solid mechanics [:cite:t:`Pierron2012`]. This python package is limited to the particular case of load reconstruction during fast transient dynamics.
The procedure for surface pressure reconstruction is based on the work of :cite:t:`Pierron2012` and :cite:t:`Kaufmann2019,Kaufmann2019`.

The dynamic equilibrium equations for a thin plate can be written written on weak form using the principal of virtual work as :cite:`dym1973solid`:

.. math::
   :nowrap:

    \begin{equation}
    \underbrace{\int\limits_{V} \rho ~ \boldsymbol{a} ~ \boldsymbol{u^*} dV}_{W_{inertial}^*} ~=~ \underbrace{ -\int\limits_{V} \boldsymbol{\sigma} : \boldsymbol{\varepsilon ^*} dV}_{W_{int}^*} + \underbrace{ \int\limits_{S} \overline{\boldsymbol{T}} \boldsymbol{u^*} ~dS + \int\limits_{V} \rho ~ \boldsymbol{F_{Vol}} ~ \boldsymbol{u^*} ~dV}_{W_{ext}^*} ~~~,
    \end{equation}

where :math:`W_{inertial}^*`, :math:`W_{int}^*` and :math:`W_{ext}^*` denotes the inertial virtual work, the internal virtual work and the external virtual work, respectively.

For the particular case of a thin plate represented by an isotropic linear elastic material, the principal of virtual work can be written using the Kirchoff-Love theory as:

.. math::
   :nowrap:

    \begin{equation}
	    \begin{aligned}
	    \int\limits_{S} p w^{*} dS ~ =
	    ~& D_{xx} \int\limits_{S} \left( \kappa _{xx} \kappa ^* _{xx} +\kappa _{yy} \kappa ^* _{yy} + 2 \kappa _{xy} \kappa ^* _{xy} \right) dS \\
	    +~&D_{xy} \int\limits_{S} \left( \kappa _{xx} \kappa ^* _{yy} +\kappa _{yy} \kappa ^*  _{xx} -2 \kappa _{xy} \kappa ^* _{xy} \right) dS
	    +~\rho ~ t_S \int\limits_{S} a ~w^{*} dS ~~~.
	    \end{aligned}
    \end{equation}

where :math:`S` is the surface of the plate, :math:`p` is the pressure acting on the surface of the plate.
The deformation of the plate is given by the curvatures :math:`\kappa` and the acceleration :math:`a`. The density of the plate material is denoted :math:`\rho`, the plate thickness is denoted :math:`t_S`, and :math:`D_{xx}` and :math:`D_{xy}` are the plate bending stiffness matrix components. Virtual quantities are marked with :math:`^*`.

As local pressure values are of interest, the surface is divided into subdomains. By assuming a constant pressure distribution within each subdomain, the integrals in the above equation is reformulated as discrete sums and the pressure :math:`p` is solved for:

.. math::
   :nowrap:

    \begin{equation}
    \begin{aligned}
    p ~ = ~
    \Biggl( ~&D_{xx} \sum\limits_{i = 1} ^{N}  \kappa ^{i} _{xx} \kappa ^{*i} _{xx} +\kappa _{yy}^{i} \kappa ^{*i}  _{yy} +2\kappa ^{i}_{xy} \kappa ^{*i} _{xy} \\
    +~&D_{xy} \sum\limits_{i = 1} ^{N} \kappa ^{i}_{xx} \kappa ^{*i} _{yy} +\kappa ^{i}_{yy} \kappa ^{*i} _{xx}
    -2\kappa ^{i}_{xy} \kappa ^{*i} _{xy}
    +~\rho ~ t_S \sum\limits_{i = 1} ^{N} a^{i} ~w^{*i} \Biggr) ~ \left( \sum\limits_{i = 1} ^{N} w^{*i} \right) ^{-1} ~~~,
    \end{aligned}
    \end{equation}

where :math:`N` is number of discrete surface elements.

The virtual fields based on 4-node Hermite 16 element shape functions :cite:t:`zienkiewicz1977` are available for pressure reconstruction, see :cite:t:`Pierron2012` for more details.

Bibliography
------------

.. bibliography::
   :style: plain
Installation
=============
In order to get started with Recolo, you need to install it on your computer.

By cloning the repo:
---------------------

These instructions will get you a copy of the project up and running on your
local machine for development and testing purposes.

Prerequisites::

    This toolkit is tested on Python 3.7
    We recommend the use of virtualenv

Start to clone this repo to your preferred location::

   git clone https://github.com/PolymerGuy/recolo.git


We recommend that you always use virtual environments, either by virtualenv or by Conda env

Virtual env::

    python -m venv env
    source ./env/bin/activate #On Linux and Mac OS
    env\Scripts\activate.bat #On Windows
    pip install -r requirements.txt


You can now run an example::

    $ python ./examples/AbaqusExperiments/Reconstruction_minimal.py

Running the tests
------------------
The tests should always be launched to check your installation.
These tests are integration and unit tests

If you cloned the repo, you have to call pytest from within the folder::

    python -m pytest --pyargs recolo --cov=./


{{ fullname | escape | underline}}

.. automodule:: {{ fullname }}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Module Attributes') }}

   .. autosummary::
      :toctree:
   {% for item in attributes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block functions %}
   {% if functions %}
   .. rubric:: {{ _('Functions') }}

   .. autosummary::
      :toctree:
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: {{ _('Classes') }}

   .. autosummary::
      :toctree:
      :template: custom-class-template.rst
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   .. rubric:: {{ _('Exceptions') }}

   .. autosummary::
      :toctree:
   {% for item in exceptions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

{% block modules %}
{% if modules %}
.. rubric:: Modules

.. autosummary::
   :toctree:
   :template: custom-module-template.rst
   :recursive:
{% for item in modules %}
   {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :members:
   :show-inheritance:
   :inherited-members:

   {% block methods %}
   .. automethod:: __init__

   {% if methods %}
   .. rubric:: {{ _('Methods') }}

   .. autosummary::
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Attributes') }}

   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}