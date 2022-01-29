# PyAstroPol
Instrumental Polarization Analysis of Astronomical Optics

[![DOI](https://zenodo.org/badge/288481693.svg)](https://zenodo.org/badge/latestdoi/288481693)

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02693/status.svg)](https://doi.org/10.21105/joss.02693)

## Overview
The package has one simple goal : compute 4x4 Mueller matrix for the given optical system, and it is developed keeping astronomical optics in view.
It uses geometric optics approach i.e., all the analysis uses strictly ray treatment. Users should keep in mind that this is NOT for optical design i.e., it is presumed that the user already knows the optical system that is to be analyzed.

The package imports following external libraries, all of which are ubiquitous. They accompany any decent scientific `Python` distribution hence this package should function with virtually no dependency issues however, it should be noted that it is developed with `Python3.6`.
```python
numpy
matplotlib
```
**Documentation on the `Classes` is hosted at [ReadTheDocs](https://pyastropol.readthedocs.io/index.html)**.




## Getting Started

### Installation

The package has following distinct components :
1. *Code Base*: The directory `PyAstroPol/PyAstroPol`. It containes source code of the packge that is only accessed by the application and not by the user.
2. *User Data*: The directories `PyAstroPol/Examples` and `PyAstroPol/Materials`. They contain data pertaining to the package that user should be able to access, which may also be used by the application in the runtime.

The present *installation* scheme is outlined below. The users are requested to provide their valuable feedback about the preferences regarding the installation (e.g., having two directories for root and user data, installing in development mode etc.). This will be of great help in devising a better installation scheme for the next versions. 

Follow these steps to start using the package.

1. Go to the user directory where the package is to be installed.  
`cd <User directory>`   
Download the package from the Github.   
`git clone https://github.com/hemanthpruthvi/PyAstroPol.git`  
Rename the top directory from `PyAstroPol.git` to `PyAstroPol`

2. Go to `PyAstroPol` root directory  
`cd <User directory/PyAstroPol>`  
Install the required dependencies by running    
`pip install -r requirements.txt`

3. Add `<User directory>/PyAstroPol` to the `PYTHONPATH` environment variable.  
In Windows systems, this option can be found at `Control Panel > All Control Panel Items > System > Advanced system settings > Environment Variables`   
In Linux systems, this can be done with the command line  
`export PYTHONPATH=<User directory>/PyAstroPol`

4. Import this package to your `Python` script using   
`import PyAstroPol as pap`

### Examples

[PyAstroPol/Examples/](https://github.com/hemanthpruthvi/PyAstroPol/tree/master/Examples) contains several examples files to demonstrate the applications of the package. They are provided in the form of `IPython` notebooks, and running them is a good way to quick-start using the package. **They also function as the test cases**. 

### Analysis on own

As previously mentioned, this is not a design software. Hence, one needs to know the optical system they wish to analyze. As per the framework of the `PyAstroPol` optical system, there are thee types of objects :
1. Source  
2. Components  
3. Detector   

Following steps illustrate how to devise a simple optical system.  
1. import the required modules.
```python
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import PyAstroPol as pap
```  
2. Create a source, and optionally create a source for display. For analysis one can define a source with a lot of rays (say 10000), and for display one can define a source with fewer rays (say 10).  
```python
S_analysis = pap.Source(10000, Clear=20)        # Source for analysis, with 10k rays and 20 mm size
S_display = pap.Source(10, Clear=20)            # Source for disply, with 10 rays and 20 mm size
```  
3. Create a component such as surface, lens etc., and position it. 
```python
L = pap.UncoatedLens(50, Thick=10, R1=200, R2=-200)     # Simple bi-convex lens of 50 mm size
L.translateOrigin(z=100.0)                              # Move the lens from default position (origin)
```  
4. Create a detector and position it.
```python
D = pap.Detector(50)                    # Detector of size 50 mm
D.translateOrigin(z=200.0)              # Move the detector from default position (origin)
```  
5. Put them together to create the optical system.
```python
O_system = pap.System(S_analysis, [L], D, dRays=S_display)
O_system.propagateRays()                                # Propagate rays in the optical system
```  
6. Display the optical system using matplotlib 3d axis.
```python
Fig = plt.figure()
Ax = Fig.add_subplot(111, projection='3d')
O_system.draw(Ax)
plt.show()
```  
7. Compute the Mueller matrix and print it.
```python
MM, T = O_system.getSystemMuellerMatrix()       # Compute Mueller matrix for the system
print(MM)
```

## Directories

[PyAstroPol/PyAstroPol/](https://github.com/hemanthpruthvi/PyAstroPol/tree/master/PyAstroPol)  
It is the main directory containing all the source files.

[PyAstroPol/Materials/](https://github.com/hemanthpruthvi/PyAstroPol/tree/master/Materials)  
It has the refractive index data for different materials in a formatted manner. These files are loaded by the code to look-up the refractive index information of the given material. Users can easily create such files using following steps.
1. Download wavelength vs refractive index file as `.csv` from popular refractive index database [RefractiveIndexInfo](https://refractiveindex.info/).
2. Rename the file to an appropriate material name e.g., for Aluminium the file name is `Al.csv`.
3. Copy the `.csv` file into `PyAstroPol/Materials/` directory.
4. Format the material file using provided function i.e., `formatMaterialFile(MaterialName)` (without file extensions).
5. The material is ready to be used by the code e.g., `M1 = Surface(50, n2='Al', Mirror=True)`.
6. An [example](https://github.com/hemanthpruthvi/PyAstroPol/blob/master/Examples/09_FormatMaterialFile.ipynb) is also provided.

[PyAstroPol/Docs/](https://github.com/hemanthpruthvi/PyAstroPol/tree/master/Docs)  
It contains documentation related codes and files.
`Theory_and_Implementation_Notes.ipynb` details the formulation behind the codes. **Users interested in development are encouraged to refer this document**.

## Conventions used in this package  
The most important aspect to remember while using the package is **the convention**, which is described below. 
### For astronomy : 
Positive X-axis : West  
Positive Y-axis : Zenith  
Positive Z-axis : North  
Positive Latitude : North  
Positive Hour Angle : West  
Positive Declination : North  
### For optics : 
Complex refractive index is **n-*i*k** where **n** and **k** are positive real numbers.    
Jones vector corresponding to positive Stokes-V is <img src="https://render.githubusercontent.com/render/math?math=\frac{1}{\sqrt 2} \begin{bmatrix} 1 \\ -i \end{bmatrix}">.

## Contributing
Any mode of contribution is highly encouraged.
1. Bug reporting : Open an issue in github with the following details.
    - Description of the bug
    - Python, numpy and matplotlib versions
    - Operating system details
    - Snippet of the code causing the issue
2. Feature request : Open an issue in github with the following details.
    - Description of the feature
    - Description of the application
    - If possible, an example
3. Example request : Open an issue in github with following details.
    - Description of the application
    - Expected output from the example
4. Bug fixes : Open a pull request in github with following details.
    - Description of the bug corresponding to the fix
5. Feature addition : Open a pull request in github with following details.
    - Description of the feature
    - Description of the application
    - At least one Example using the particular feature
6. Other : Open an issue on github with a description.

Kindly use appropriate Tags as well.

## TODO
1. Add polarizing elements such as birefringent elements, waveplates and polarizers.
2. Add feature to create and save coatings as files.
3. Add rectandular and elliptical apertures.## References
Although the materials are downloaded from [https://refractiveindex.info/](https://refractiveindex.info/), the original source of the data is very much different. They are listed here for the sake of reference and reproducability.

**Ag** and **Al** : 
McPeak et.al., 2015 ([DOI](https://doi.org/10.1021/ph5004237)).

**FusedSilica** : 
Malitson, 1965 ([DOI](https://doi.org/10.1364/JOSA.55.001205)).

**SodaLime** : 
Rubin, 1985 ([DOI](https://doi.org/10.1016/0165-1633(85)90052-8)).

**Zerodur&reg;*** : 
SCHOTT, TIE-43 ([Source](https://refractiveindex.info/download/data/2007/schott_tie-43_optical_properties_of_zerodur_november_2007_us.pdf)).

**N-XXXXX*** glasses : 
SCHOTT catalog ([Source](https://www.schott.com/d/advanced_optics/ac85c64c-60a0-4113-a9df-23ee1be20428/1.17/schott-optical-glass-collection-datasheets-english-may-2019.pdf)).


\* The refractive index values are computed by using Sellmeier dispersion model with coefficents stated by the manufacturer.## References
Although the materials are downloaded from [https://refractiveindex.info/](https://refractiveindex.info/), the original source of the data is very much different. They are listed here for the sake of reference and reproducability.

**Ag** and **Al** : 
McPeak et.al., 2015 ([DOI](https://doi.org/10.1021/ph5004237)).

**FusedSilica** : 
Malitson, 1965 ([DOI](https://doi.org/10.1364/JOSA.55.001205)).

**SodaLime** : 
Rubin, 1985 ([DOI](https://doi.org/10.1016/0165-1633(85)90052-8)).

**Zerodur&reg;*** : 
SCHOTT, TIE-43 ([Source](https://refractiveindex.info/download/data/2007/schott_tie-43_optical_properties_of_zerodur_november_2007_us.pdf)).

**N-XXXXX*** glasses : 
SCHOTT catalog ([Source](https://www.schott.com/d/advanced_optics/ac85c64c-60a0-4113-a9df-23ee1be20428/1.17/schott-optical-glass-collection-datasheets-english-may-2019.pdf)).


\* The refractive index values are computed by using Sellmeier dispersion model with coefficents stated by the manufacturer.---
title: 'PyAstroPol: A Python package for the instrumental polarization analysis of the astronomical optics.'

tags:
  - Python
  - Astronomy
  - Polarization
  - Optics

authors:
  - name: Hemanth Pruthvi.
    orcid: 0000-0002-4892-6561
    affiliation: "1"

affiliations:
 - name: Leibniz-institut für Sonnenphysik, Freiburg, Germany.
   index: 1

date: 25 August 2020

bibliography: paper.bib
---

# Statement of Need

The instrumental polarization analysis is one of the key aspects of the optical system analysis of the astronomical telescopes and instruments. The majority of the optical surfaces in the astronomical instruments induce a change in the state of polarization of the incoming light, which could lead to inaccurate measurements of the state of polarization or _polarimetry_. Hence, polarimetric instruments inevitably require the information about the polarization properties of the optical system or the _instrumental polarization_. This can be determined through experiments, analysis or a combination of both. As the telescopes become larger and instruments become more complex, the instrumental polarization analysis has become all the more crucial [e.g., @TMT_2015; @DKIST_2016].  

# Summary

The Python package `PyAstroPol` provides a means to analyze the polarization properties of a given optical system with relative ease and minimal dependencies. The simple end goal is to calculate the Mueller matrix [e.g., @Gil_2016] of the optical system.  
Data pipelines of the astronomical instruments can easily distribute this package alongside their polarimetric calibration routines.

In the polarization analysis of the astronomical telescopes, various approaches have been adopted depending on the complexity of the system. A significant part of the complexity is because the instrumental polarization is often time-dependent. For the solar telescopes with Coelostats, Mueller matrices were analytically derived as a function of time [e.g., @KTT_1985; @VTT_2005]. For Thirty Meter Telescope (TMT), a combination of analytical and numerical methods is used [@TMT_2015]. However, for Daniel K. Inouye Solar Telescope (DKIST), Zemax – a commercial software – is used, along with the team's in-house tools [@DKIST_2016].

`PyAstroPol` aims to combine the better parts of the aforementioned programs: it is an open-source tool written in the `Python` programming language and offers a range of functions to model the astronomical optics. There is no open-source software which has polarization ray propagation features as per my knowledge. Hence, Zemax OpticStudio&reg; modelling is used for comparison. The salient features and important limitations are listed below.

Salient features :   

* `PyAstroPol` calculates the Mueller matrix of a given imaging-type optical system, by coherently adding the electric field vectors after propagation.
* Astronomical sources can be directly placed in the model using relevant coordinates, namely, declination, hour angle and latitude of the telescope site.
* Off-axis components are included, as they have a significant effect on polarization.
* Effect of multi-layered coatings, such as oxide layers and protective coatings, on the state of polarization is included. These are also significant in the polarization analysis of the astronomical optics [e.g., @VanHarten_2009].
* All the data, such as points of incidence, polarization directions, complex electric field values and more, are readily available to the user for any further analysis.
* Material refractive index information can be downloaded from the popular online source [https://refractiveindex.info/](https://refractiveindex.info/) as `.csv`. It can be formatted and used with this software.
* Spot diagram is possible at any instance of the beam path, as a by-product of ray tracing.

Important Limitations:   

* All the analysis uses strictly the ray treatment. Hence, all the limitations of the rays optics shall be applicable.
* Only circular optics (apertures) can be modelled at the moment.
* Birefringent components are not included yet. The justification for this choice is that the behaviour of the birefringent components is fairly straightforward as they strongly polarize the light.
* The visualization features are limited.
* `PyAstroPol` is neither a design software nor interactive. That is, the user must know the optical system that is to be analyzed, and the system must be updated every time a component is changed.     

Results of the code have been verified against previous works [@Pruthvi_2018], and Zemax OpticStudio&reg;. A variety of examples have been provided, and they should facilitate the quick-start. One of the examples also provides the aforementioned comparison with commercial software. 

# Acknowledgements

This work has been carried out as a part of the ongoing project _Jets in the solar atmosphere_, funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) – Projektnummer 407727365. I thank DVS Phanindra (IIA, Bengaluru) and V. Sreekanth Reddy (CHESS, Hyderabad) for the discussion, and Mathias Waidele (KIS, Freiburg) for his feedback. I also thank the PI of the DFG project Markus Roth (KIS, Freiburg) for the opportunity to carry out this work. I thank reviewers [`@caldarolamartin`](https://github.com/caldarolamartin), [`@aquilesC`](https://github.com/aquilesC) and [`@mwcraig`](https://github.com/mwcraig), and editor [`@pibion`](https://github.com/pibion) for their valuable feedback in improving this package. 

# References
