# **Quantification of lidar uncertainties - Qlunc**



[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4579842.svg)](https://doi.org/10.5281/zenodo.4579842)



## What is `Qlunc`?
`Qlunc` is a python-based, open, freely available software that aims to quantify errors when measuring with a lidar device. The code has an objected-oriented structure; by using python objects and simulating real lidar components the code puts all together in modules to eventually build up a lidar [digital twin](https://en.wikipedia.org/wiki/Digital_twin). The code is meant to be as modular as possible and offers the possibility of creating different lidar objects on parallel (see [Tutorial2.ipynb](https://github.com/SWE-UniStuttgart/Qlunc/blob/main/Tutorials/Tutorial2.ipynb)), with different components at the same time. This allows to easily combine different modules with different characteristics simulating different lidar devices.

![Qlunc basic structure image](https://github.com/SWE-UniStuttgart/Qlunc/blob/main/Pictures_repo_/Qlunc_GralStructure.JPG)

Currently, the code can calculate uncertainties coming from photonics, including photodetector (with or without trans-impedance amplifier) and optical amplifier uncertainties, as well as optics module uncertainty including scanner pointing accuracy, probe volume assessment (both continuous wave and pulsed lidars) and optical circulator uncertainties. For each module the Guide to the Expression of Uncertainty in Measurement ([GUM](http://www.bipm.org/en/publications/guides/gum.html)) is applied to calculate uncertainty expansion, taking into account that components are considered uncorrelated. 

### Creating a lidar device
The user creates the different lidar components by instantiating a python `class`, including its functional parameters and defining the function that is used to obtain the specific component uncertainty. Then, each module (also python objects) is "filled" with the corresponding components and their uncertainties are computed following uncertainty expansion method according to the GUM model. Once each component is 'ensembled' building up the different modules, the lidar object is created and the modules included. As a result, the desired lidar digital twin is created, the uncertainty of which is computed again by following [GUM](http://www.bipm.org/en/publications/guides/gum.html) suggestions about uncertainty expansion.

### Creating atmospheric conditions
The user creates also atmospheric scenarios to account for the different atmospheric conditions the lidar has to deal with. Atmospheric inputs, basically temperature 
and humidity, either single values or time series coming from peripherals are both accepted.

### `Qlunc` available capabilities

#### Uncertainties
The next step is to ask for the uncertainty we are interested in, either coming from a component, module or lidar object. Indeed, the flexibility of the code allows the user not just to assess global lidar uncertainty, but also to query uncertainties coming from specific modules or even single components.

#### Plots
 - Can draw photodetector uncertainties comparison including shot noise, thermal noise, dark current noise and, if needed, trans-impedance amplifier noise.
 - Scanning points and their uncertainty in meters (VAD and Scanning lidar).
 - Optical signal to noise ration from the optical amplifier

## How to use `Qlunc`

:warning: **Please downolad the latest release (V0.91).**

### Create an environment and install dependencies

1) Having [Anaconda](https://docs.anaconda.com) installed is a prerequisite if we want to work in a different environment than `base`, and it is recommended. Then, based on the requirements added to the ``environment.yaml`` file on the repository, where are included the name of the environment and the tools/packages we want to install, we build the new environment. 

2) In the Anaconda prompt, go to the directory where you have clone/download `Qlunc` and type:

```
conda env create -f environment.yml 
conda activate <envname>
```

3) Your environment is ready to rumble. You have now a new environment, called `Qlunc_Env` by default, with all the packages needed to run `Qlunc`.

4) In case you don't want to create a new environment, just install the requirements listed in the *Requirements* section below.

### Download or [clone](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository) the repository to a local directory

By downloading or cloning the repository you will get several folders within which `Qlunc` is organized. The most importants to know are:

### Main
This is the core of `Qlunc`. Here the user creates the classes describing the components, modules and general inputs of the lidar device and instantiate the classes.
 - `Template_yaml_inputs_file.yml` and `Qlunc_inputs.yml`. The former is a yaml template where user introduces the lidar components values, modules and general lidar features as well as atmospheric scenarios; the latter can be taken as an example showing how to fill in the template.
 - `Qlunc_Classes.py` contains the code which _creates_ all the lidar digital twins. Each lidar module/component is assigned to a python `class`.
 - `Qlunc_Instantiate.py` instantiate the lidar classes taking the values from `Qlunc_inputs.yml`.
### UQ_Functions
 - Contains the functions that compute the uncertainties coming from different devices, calculting also the uncertainty propagation corresponding to the different      modules and lidar uncertainty as well. Users can define their own functions to calculate specific module uncertainties, and combined/expanded uncertainties as well. 
### Utils
 - Contains scripts meant to do different tasks. Importing packages and stand alone funtions which don´t interface directly with `Qlunc` but are necessary to compute calculations. Also contains a `Qlunc_Plotting.py` script to automate plots and `Scanning_patterns.py` to introduce users' pre-defined scanning patterns.
###  TestFile_Qlunc
 - A working example is provided to show how the process looks like. In this test case, a lidar is built up with its modules and components, puting all together to set up a lidar device. User can find more information on how to run this test file in the `readme.md` file dropped in this folder.
### Tutorials
- Containing 3 [Jupyter Notebook-based tutorials](https://github.com/SWE-UniStuttgart/Qlunc/tree/Qlunc-V0.9/Tutorials); `Turial0.ipynb`, `Tutorial1.ipynb` and `Tutorial2.ipynb` with their corresponding yaml files. The tutorials are also available through the Binder service to ease accessibility and reproducibility. Users can find more information about these tutorials in the corresponding `readme.md` file dropped in this folder.
## Requirements
The following python libraries and tools should be installed beforehand and are included in the `environment.yml` file:

- matplotlib==3.2.1
- numpy==1.18.5 
- pandas==1.2.1
- pyyaml==5.4.1
- scipy==1.6.0
- sympy==1.7.1
- xarray==0.15.1
- xarray-extras==0.4.2
- python==3.7.9
- spyder==4.2.1
- netcdf4==1.5.3
- notebook
- jupyterlab
- termcolor==1.1.0

## Author
[Francisco Costa](https://www.ifb.uni-stuttgart.de/en/institute/team/Costa-Garcia/)

## Contributions
Contributions are very welcome!
If you are wishing to:
- Colaborate: please contact the author
- Report problems or enhance the code: post an [issue](https://docs.github.com/en/issues/tracking-your-work-with-issues/quickstart) or make a [pull request](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request)
- Seek for support: please contact the author

## License
`Qlunc` is licensed under **[SD 3-Clause License](https://github.com/SWE-UniStuttgart/Qlunc/blob/main/LICENSE)**

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Citing and Contact

<div itemscope itemtype="https://schema.org/Person"><a itemprop="sameAs" content="https://orcid.org/0000-0003-1318-9677" href="https://orcid.org/0000-0003-1318-9677" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">Francisco Costa García</a></div>

University of Stuttgart - Stuttgart Wind Energy
 
email: costa@ifb.uni-stuttgart.de
 
## Tutorial0:

Qlunc's presentation and basics. What's Qlunc and how does it work.

#### Try it yourself:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/SWE-UniStuttgart/Qlunc.git/HEAD?filepath=Tutorials%2FTutorial0.ipynb)

## Tutorial1:
This tutorial aims to facilitate the introduction to Qlunc. 
Will go through the code and create a lidar digital twin, with its modules and components. Will ask for uncertainties either lidar general one or component specific uncertainty. We will see some graphical interesting results. Will see how to access design lidar data.

#### Try it yourself:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/SWE-UniStuttgart/Qlunc.git/HEAD?filepath=Tutorials%2FTutorial1.ipynb)

## Tutorial2:
Along `Tutorial1` we learn how to virtually create a lidar device, being the input to the uncertainty estimation model. In `Tutorial2` we design two lidar devices and compare against each other. These devices differ from each other just in the scanner head, so we will "build up" two optic modules by changing the scanner component.

#### Try it yourself:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/SWE-UniStuttgart/Qlunc.git/HEAD?filepath=Tutorials%2FTutorial2.ipynb)
---
title: 'Qlunc: Quantification of lidar uncertainty'
tags:
  - wind lidar
  - lidar hardware uncertainty
  - OpenScience
  - OpenLidar
  - digital twin
authors:
  - name: Francisco Costa
    orcid: 0000-0003-1318-9677
    affiliation: 1
  - name: Andrew Clifton
    orcid: 0000-0001-9698-5083
    affiliation: 1
  - name: Nikola Vasiljevic
    orcid: 0000-0002-9381-9693
    affiliation: 2
  - name: Ines Würth
    orcid: 0000-0002-1365-0243
    affiliation: 1
affiliations:
 - name: Stuttgart Wind Energy (SWE), Allmandring 5b, 70569 Stuttgart, Germany
   index: 1
 - name: DTU Wind Energy, Frederiksborgvej 399, 4000 Roskilde Denmark 
   index: 2
date: xx March 2021
bibliography: paper.bib
---
# Summary

Wind lidar is a flexible and versatile remote sensing device for wind energy applications [@Hauke] that measures the wind vector remotely using laser light backscattered from aerosols. It is a key tool for wind energy and meteorology. As with any measurement method, it is essential to estimate its uncertainty.
Qlunc, which stands for **Q**uantification of **l**idar **unc**ertainty, is an open-source Python-based tool to create a digital twin of lidar hardware, and to estimate the uncertainty of wind lidar wind speed measurements.

Qlunc contains models of the uncertainty contributed by individual lidar components and modules (represented by Python objects, which in turn represent physical lidar  objects), that then are combined, considering their different natures, to estimate the uncertainties in wind lidar measurements. The modules are based on the OpenLidar architecture [@OpenLidar] and can be easily adapted for particular use cases thanks to the modularity of the code (see \autoref{fig:QluncStructure}). The terminology for the components and modules defined within Qlunc has also been aligned with a community-driven wind lidar ontology, which is in development [@OntoStack;@sheet2rdf]. 

![Qlunc basic structure.\label{fig:QluncStructure}](Qlunc_BasicStructure_diagram.png)
 
The first release is focused on velocity azimuth display (VAD)[@Browning] scans and forward-looking nacelle-mounted measuring modes, which are common wind-energy-industry applications. Besides uncertainty estimations, Qlunc’s functions could be extended for other applications, for example to compare different  wind velocity vector calculation methods. This, combined with the underlying open-source code, defines an attractive scenario for sharing knowledge and fostering collaboration on wind lidars. 

# Statement of Need

Wind lidars are measuring devices, and as for any other measuring systems, their measurements have uncertainties [@Borraccino_2016]. Therefore, as already stated, it is crucial to assess their measurement uncertainty in order to increase confidence in lidar technology.

Measurement uncertainty means doubt about the validity of the result of a measurement [@GUM]. It represents the dispersion of the values attributed to a measurand. The ability to simulate uncertainty through a model such as Qlunc is important for judging measurement data but can also be useful for designing and setting up experiments and optimizing lidar design. Because wind lidar is important for wind energy applications [@Clifton_2018], better models for wind lidar hardware (e.g., Qlunc) and measurement processes (e.g., through MOCALUM [@mocalum] or YADDUM [@yaddum], with which Qlunc can feasibly combine) will directly contribute to the adoption of wind lidar for wind energy applications. 

This project is influenced by fundamental open science principles [@OpenScience]. The scope is to create an open, standardized and collaborative framework to describe both generic and specific lidar architectures, characterize lidar uncertainties, and provide the tools for others to contribute within this framework. 
 
# Future development roadmap

Over the next year, we plan to implement further lidar hardware modules in the model and compute their combined uncertainties. In addition, we will identify main data processing methods and include those that we consider the highest contributors to uncertainty. 

We also plan to further align the terminology used in Qlunc with the IEA Wind Task 32 controlled vocabulary for wind lidar [@task32ControlledVocabulary]. This will make it easier for users to understand what each of the modules and components do, and promotes interoperability.

All documentation from the project, tutorials, and raw code will be published through a website, to enable users to dive into the numerical framework and get used to the Qlunc routines.
We welcome contributions from the wind lidar community.

# Acknowledgements

This work is part of the LIKE ([Lidar Knowledge Europe](https://www.msca-like.eu/)) project. The project LIKE H2020-MSCA-ITN-2019, Grant number 858358 is funded by the European Union.
 
# References
 


## Working example to create a lidar digital twin
:warning: To do this exercise interactively try [Tutorial0](https://github.com/SWE-UniStuttgart/Qlunc/blob/main/Tutorials/Tutorial0.ipynb)


In the following lines it is shown the `Qlunc`'s workflow and the code structure it is based on; the `WorkingExample` files are a use case example.
Now, imagine we want to create a lidar object made up with one single module. This module will contain one single component with properties Property_1 and Property_2. The steps we have to follow are: 

 1) Fill up yaml file with inputs
 2) Create a component with its propertie(s)
 3) Create a module containing the component(s)
 4) Create a lidar containing the module(s)
 5) Ask for the uncertainties we are interested in by using _dot notation_

### 1. Fill up the inputs yaml file
Before creating the classes for the different components we will fill in the yaml file with the corresponding values for the components and decide the components and the modules that we want to include in the lidar device for uncertainty calculations. Users have a [yaml template for Qlunc inputs](https://github.com/SWE-UniStuttgart/Qlunc/blob/main/Main/Template_yaml_inputs_file.yml) in the repository.

The minimum information defining your component might be:
 - Name: Provide an ID to our object
 - Property: As many as the component has (e.g. for the photodetector could be `wavelength`, `load resistor` and  `gain`).
 - Uncertainty function: Function developed by the user decribing the uncertainty of the lidar module.
 
**Warning!** When introducing the component in each module, the name should be the same as in the component instance (e.g. if the name of your module instance is _Module_A_ the name to use in the yaml file should be the same). 

  ```
   YAML file:
    >>## Components:
    >> 
    >>  Component_A:
    >>  
    >>    Name: ComponentA
    >>   
    >>    Property_A: property_1_value 
    >>    
    >>    Property_B: property_2_value 
    >>   
    >>    Uncertainty function: Uncertainty_ComponentA  # Function describing the module uncertainty in _Module_A_ due to their components.
    >>---
    >>## Modules:
    >> 
    >>  Module_A: 
    >>  
    >>    Name: ModuleA
    >>   
    >>    Component: _Component_A_                   # It has to be the same name as the instance (*).
    >>   
    >>    Uncertainty function: Uncertainty_ModuleA  # Function describing the module uncertainty in _Module_A_ due to their components.
    >>---
    >>## Lidar:
    >> 
    >>  Lidar_A:
    >>  
    >>    Name: LidarA
    >>   
    >>    Module: _Module_A_                         # It has to be the same name as the instance (**).
    >>   
    >>    Uncertainty function: Uncertainty_ModuleA  # Function describing the module uncertainty in _Module_A_ due to their components.
```

### 2. Creating the component digital twin
The components are included as python classes, for example a component, _Component_A_, is created by instantiating class _Comp_A_:

- Creating a class for the component _Component_A_:
```
  >> class Comp_A():
  >> 
  >>   def __init__(self, property_1, property_2, unc_func)
  >>   
  >>      self.property_1  = property_1
  >>      
  >>      self.property_2  = property_2
  >>      
  >>      self.uncertainty = unc_func 
``` 
- Then we instantiate class _Comp_A_ to create the object representing the lidar component digital twin:
```
  >> Component_A (*) = Comp_A (name        = C_A,
  >>   
  >>                           property_1  = property_1_value,         # picked from the yaml file
  >>                       
  >>                           property_2  = property_2_value,         # picked from the yaml file
  >>                       
  >>                           uncertainty = Component_A_uncertainty)  # parameter describing uncertainty in _Comp_A_.
```
The uncertainty function is a function either found in literature or developed by the user that discribes the uncertatinty of the component.

### 3. Creating the module digital twin
For the modules we create a class and include the 

- Creating a class for the _Module_A_:
 ``` 
  >> class Mod_A():
  >> 
  >>   def __init__(self, name, component, unc_func)
  >>   
  >>      self.name        = name
  >>      
  >>      self.component   = component   
  >>      
  >>      self.uncertainty = unc_func  
``` 
- Then we instantiate class _Mod_A_ to create the Module object:
```
  >> Module_A (**) = Mod_A (name        = ModuleA, 
  >> 
                            component   = Component_A,                      # Including _Component_A_ in the module. 
                                                                            # It has to be the same as in the instance (*).
                            uncertainty = Module_A_uncertainty_function.py) # Uncertainty describing uncertainty in _Module_A_. Following GUM.                      
```
### 4. Creating the lidar

Then once we have created the module(s), we can made up a lidar object just in the same way:

- Creating a class for the _Lidar_A_ device:
```
  >> ## Lidar:

  >> class lidar():
  >> 
  >>   def __init__(self, name, module, unc_func)
  >>   
  >>      self.name        = name
  >>      
  >>      self.module      = module
  >>             
  >>      self.uncertainty = unc_func  
```  
- Then we instantiate class _lidar_ to create the lidar object:
```
  >> Lidar_A = lidar (name        = LidarA, 
  >> 
  >>                  module      = Module_A,                         # Actually, picked from the yaml file. 
  >>                                                                  # It has to be the same name as the instance (**).                
  >>                  uncertainty = Lidar_A_uncertainty_function.py)  # Uncertainty describing uncertainty in _Lidar_A_. Following GUM.
```
Then, we have created a lidar (python-based) object called _Lidar_A_, made up of one module, _Module_A_, which contains one single component, _Component_A_, with properties _Property_1_ and _Property_2_.

### 5. Asking for uncertainties
The modularity of the code  allows user either to ask for _Photodetector1_ uncertainty (component uncertainty), _Photonics_ uncertainty (module unceratinty) or global lidar uncertainty. using the dot notation we can write:
```
>> Lidar_A.module.component.uncertainty(Lidar_A, AtmosphericScenario,cts,Qlunc_yaml_inputs) # for the component uncertainty included in Module
>> Lidar_A.module.uncertainty(Lidar_A, AtmosphericScenario,cts,Qlunc_yaml_inputs) # for the Module uncertainty
>> Lidar_A.uncertainty(Lidar_A, AtmosphericScenario,cts,Qlunc_yaml_inputs) # for the lidar global uncertainty
```
![Uncertainty_WF](https://github.com/SWE-UniStuttgart/Qlunc/blob/main/Pictures_repo_/FlowChartUnc.JPG)
**Qlunc/Main/.**
- `Qlunc_Classes.py`: Script to define the classes representing the modules and components the lidar is made of
- `Qlunc_Instantiate.py`: Script where classes are instantited to *create* the lidar objects. The modules'/components' values are taken from the `Qlunc_inputs.yml` file
- `Qlunc_inputs.yml`: .yml file where user introduces the modules/components values, define the lidar characteristics, create the atmospheric scenarios and choose the      plot options
- `Template_yaml_inputs_file.yml`: Template for the users to test

