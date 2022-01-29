<img src="images/logo.PNG" alt="logo" height="250"/>

# AstroPaint
_A python package for painting the sky_ 

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/syasini/AstroPaint/master?filepath=tutorial.ipynb)
[![Documentation Status](https://readthedocs.org/projects/astropaint/badge/?version=master)](https://astropaint.readthedocs.io/en/master/?badge=master)
![Python package](https://github.com/syasini/AstroPaint/workflows/Python%20package/badge.svg?branch=develop&event=push)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4243176.svg)](https://doi.org/10.5281/zenodo.4243176)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02608/status.svg)](https://doi.org/10.21105/joss.02608)

You can install **AstroPaint** by running the following in the command line:

`git clone https://github.com/syasini/AstroPaint.git`

`cd AstroPaint`

`pip install [-e] .` 

the `-e` argument will install the package in editable mode which is suitable for developement. If you want to modify the code use this option.

**Important Note**:
If you want the sample catalogs to be cloned automatically
 along with the
 rest of the repository, make sure you have [Git Large File Storage (`git lfs`)](https://git-lfs.github.com/) installed. 

If you are a conda user, please consider creating a new environment before
 installation:
 
 `conda create -n astropaint python=3.7`
 
 `conda activate astropaint`


# Workflow

Converting catalogs to mock maps with AstroPaint is extremely simple. Here is what an example session looks like:

```python
from astropaint import Catalog, Canvas, Painter

catalog = Catalog(data=your_input_data)

canvas = Canvas(catalog, nside)

painter = Painter(template=your_radial_profile)

painter.spray(canvas)
```

That's it! Now you can check out your masterpiece using

`canvas.show_map()`

![BG](images/BG_websky_cover.png)

# What is AstroPaint?

AstroPaint is a python package for generating and visualizing sky maps of a wide range of astrophysical signals 
originating from dark matter halos or the gas that they host. AstroPaint creates a whole-sky mock map of the 
target signal/observable, at a desired resolution, by combining an input halo catalog and the radial/angular 
profile of the astrophysical effect. The package also provides a suite of tools that can facilitate analysis
 routines such as catalog filtering, map manipulation, and cutout stacking. The simulation suite has an 
 Object-Oriented design and runs in parallel, making it both easy to use and readily scalable for production 
 of high resolution maps with large underlying catalogs. Although the package has been primarily developed 
 to simulate signals pertinent to galaxy clusters, its application extends to halos of arbitrary size or 
 even point sources.

# Package Structure 

See our [documentation](https://astropaint.readthedocs.io/) and [this
 chart](https://www.mindmeister.com/1417665103/astropaint-astropaint-py?fullscreen=1)
 to understand the package structure and see what methods are available so
  far. 


# Examples

## Nonsense Template

Here's an example script that paints a nonsense template on a 10 x 10 [sqr deg]
 patch of the `Sehgal` catalog: 


```python
import numpy as np
from astropaint import Catalog, Canvas, Painter

# Load the Sehgal catalog
catalog = Catalog("Sehgal")

# cutout a 10x10 sqr degree patch of the catalog
catalog.cut_lon_lat(lon_range=[0,10], lat_range=[0,10])

# pass the catalog to canvas
canvas = Canvas(catalog, nside=4096, R_times=5)

# define a nonsense template and plot it
def a_nonsense_template(R, R_200c, x, y, z):
    
    return np.exp(-(R/R_200c/3)**2)*(x+y+z)

# pass the template to the painter
painter = Painter(template=a_nonsense_template)

# plot the template for halos #0, #10, and #100 for R between 0 to 5 Mpc 
R = np.linspace(0,5,100)
painter.plot_template(R, catalog, halo_list=[0,10,100])
```
<p align="center">
<img src="images/a_random_template.png" alt="template" height="300"/>
</p>
The painter automatically extracts the parameters `R_200c` and `x,y,z
` coordinates of the halo from the catalog that the canvas was initialized
 with. Let's spray ths canvas now:
 
```python
# spray the template over the canvas
painter.spray(canvas)

# show the results
canvas.show_map("cartview", lonra=[0,10], latra=[0,10])
```
<p align="center">
<img src="images/a_random_map.png" alt="map" height="400"/>
</p>

_Voila!_

You can use the `n_cpus` argument in the spray function to paint in parallel and speed things up! 
Setting `n_cpus=-1` uses all the available cpus.   

<p align="center">
<img src="images/parallel.gif" alt="parallel" width="450"/>
</p>
   
## Stacking
You can easily stack cutouts of the map using the following:

```python
deg_range = [-0.2, 0.2] # deg
halo_list = np.arange(5000) # stack the first 5000 halos

# stack the halos and save the results in canvas.stack
stack = canvas.stack_cutouts(halo_list=halo_list, lon_range=deg_range, lat_range=deg_range)

plt.imshow(canvas.stack)
```
<p align="center">
<img src="images/a_random_stack.png" alt="stack" height="300"/>
</p>
 If this is taking too long, use `parallel=True` for *parallel stacking*. 

## Line-Of-Sight integration of 3D profiles

AstroPaint only allows you to paint 2D (line-of-sight integrated) profiles on
 your catalog halos, so if you already have the analytical expression of
  the projected profile you want to paint, we are in business. However, not
   all 3D profiles can be LOS integrated analytically (e.g. generalized NFW
    or Einasto, etc), and integrating profiles numerically along every
     single LOS is generally expensive. In order to alleviate this problem, AstroPaint offers two python decorators
 `@LOS_integrate` and `@interpolate` which make 3D -> 2D projections effortless.
 
 To convert a 3D profile into a 2D LOS integrated profile, all you need to do
  is add the `@LOS_integrate` to the definition.
  
 For example, here's how you can turn a 3D top hat profile 
 
 ```python
def tophat_3D(r, R_200c):
    """Equals 1 inside R_200c and 0 outside"""
    
    tophat = np.ones_like(r)
    tophat[r > R_200c]=0 
    
    return tophat
```

into a 2D projected one:

```python  
from astropaint.lib.utilities import LOS_integrate

@LOS_integrate
def tophat_2D(R, R_200c):
    """project tophat_3D along the line of sight"""

    return tophat_3D(R, R_200c)
``` 
This function integrates the `tophat_3D` function along every single line of
 sight. If you have many halos in a high resolution map, this can take
  forever. The trick to make this faster would be to integrate along a
   several LOSs and interpolate the values in between. This is what the
    `@interpolate` decorator does. So, a faster version of the `tophat_2D
    ` function can be constructed as the following:
    

```python  
from astropaint.lib.utilities import interpolate

@interpolate(n_samples=20)
@LOS_integrate
def tophat_2D_interp(R, R_200c):
    """project and interpolate tophat_3D along the line of sight"""
 
    return tophat_3D(R, R_200c)
```   
This is much faster, but the speed comes at a small price. If your 3D profile
 is not smooth, the interpolated 2D projection will slightly deviate from the
  exact integration. 
 <p align="center">
 <img src="images/tophat_interp.png" alt="interp" height="300"/>
 </p>
You can minimize this deviation by increasing the `n_samples` argument of the
 `@interpolate` decorator, but that will obviously decrease the painting speed.
 
 Does this plot agree with what you would expect a LOS integrated top hat
  profile (a.k.a. a solid sphere) to look like? 

## Painting Optical Depth and kSZ Profiles on the WebSky Catalog  

Let's use the `Battaglia16` gas profiles to paint tau (optical depth) and
 kinetic Sunyaev-Zeldovich (kSZ) on the WebSky catalog halos. 
 
 ```python
from astropaint.profiles import Battaglia16
 
 tau_painter = Painter(Battaglia16.tau_2D_interp)
```
 
 Since the shape of the profile is smooth, we won't lose accuracy by using the
  interpolator. 
<p align="center">
<img src="images/battaglia16_tau.png" alt="tau" height="300"/>
</p> 

Let's paint this on a 5x5 sqr deg patch of the WebSky catalog with a mass
 cut of 8E13 M_sun. 
 
 ```python
catalog = Catalog("WebSky_lite")
catalog.cut_lon_lat(lon_range=[5,10], lat_range=[5,10])
catalog.cut_M_200c(8E13)

canvas = Canvas(catalog, nside=8192, R_times=3)

tau_painter.spray(canvas)
``` 
<p align="center">
<img src="images/tau_map_battaglia.png" alt="tau_map" height="300"/>
</p>
The `Battaglia16.kSZ_T` function uses this tau and multiplies it by the
 dimensionless velocity of the halos to get the kSZ signal. 
 
```python 
kSZ_painter = Painter(Battaglia16.kSZ_T)
kSZ_painter.spray(canvas)
```
And here is what it looks like:
<p align="center">
<img src="images/ksz_map_battaglia.png" alt="ksz_map" height="300"/>
</p>


# Art Gallery 

Just because AstroPaint is developed for probing new science and doing
 serious stuff, it doesn't mean you can't have fun with it! Check out our
  [cool web app](https://astropaint-art-gallery.herokuapp.com/) to get your hands dirty with some paint. 

**Made with AstroPaint**

<img src="images/blue_drops.png" height="250">  <img src="images/spongy_terror.png" height="250">  <img src="images/burning_twilight.png" height="250">


# How to contribute

If you would like to contribute to AstroPaint, take the following steps:

1) Fork this repository
2) Clone it on your local machine
3) Create a new branch (be as explicit as possible with the branch name)
4) Add and Commit your changes to the local branch
5) Push the branch to your forked repository
6) Submit a pull request on this repository

See [this repository](https://github.com/firstcontributions/first-contributions) or [Kevin Markham's step-by-step guide](https://www.dataschool.io/how-to-contribute-on-github/) for more detailed
 instructions. 

Developement happens on the `develop` branch, so make sure you are always in sync with the latest version and submit your pull requests to this branch. 

---
title: 'AstroPaint: A Python Package for Painting Halo Catalogs into
 Celestial Maps'
tags:
  - python
  - astrophysics
  - simulation
  - visualization
  - extragalactic foregrounds
authors:
  - name: Siavash Yasini^[corresponding author]
    orcid: 0000-0003-1978-6325
    affiliation: 1 
  - name: Marcelo Alvarez 
    affiliation: "2, 3"
  - name: Emmanuel Schaan 
    orcid: 0000-0002-4619-8927
    affiliation: "2, 3"
  - name: Karime Maamari
    affiliation: "1, 5"
  - name: Shobeir K. S. Mazinani
    affiliation: 4
  - name: Nareg Mirzatuny
    affiliation: 1
  - name: Elena Pierpaoli
    affiliation: 1
affiliations:
 - name: University of Southern California 
   index: 1
 - name: Lawrence Berkeley National Laboratory 
   index: 2
 - name: University of California, Berkeley 
   index: 3
 - name: Aetna Inc.
   index: 4
 - name: Argonne National Lab 
   index: 5
date:  31 July 2020
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Overview 

`AstroPaint` is a python package for generating and visualizing
 sky maps of a wide range of astrophysical signals originating from dark
  matter halos or the gas that they host. `AstroPaint` creates a whole-sky mock map of 
 the target signal/observable, at a desired resolution, by combining an input
  halo catalog and the radial/angular profile of the astrophysical effect
   (see the
   workflow
   section for details). 
  The package also provides a suite of tools that
     can facilitate analysis routines such as catalog filtering, map manipulation, 
     and cutout stacking. The simulation suite has an Object-Oriented design and
      runs in parallel, making it both easy to use and readily scalable for
       production of high resolution maps with large underlying catalogs. Although the package has been
        primarily developed to simulate signals pertinent to galaxy clusters, its application extends to halos of arbitrary size or even point
         sources. 
             
          

![Map of the Birkinshaw-Gull effect painted with AstroPaint using the
 WebSky catalog \label{fig:BG}](../images/BG_websky_cover.png)

# Statement of Need 

Studying the large scale structure of the universe heavily relies on
 observations of astrophysical signals at various frequencies. Examples of such
  studies include detection or characterization of objects such as galaxies, clusters, or voids
   through either gravitational lensing, electromagnetic scattering, absorption or emission events in the optical, radio, or x-ray
    frequency bands. Such studies typically require simulated high resolution
     maps of various astrophysical effects to emulate both the signal and
      noise (foregrounds) components. For example, in a study that aims	
	to evaluate the detection significance of the Birkinshaw-Gull (BG)
	effect – a probe of the transverse
	velocities of halos [@Birkinshaw:1983; @Yasini:2018] – using the Simons
	Observatory [@SO:2019] or CMB-S4 [@CMB-S4:2019], one needs a mock
	map of the BG effect
	\autoref{fig:BG}
	as well as maps of potential contaminants such as kinetic and
	thermal Sunyaev-Zeldovich effects (kSZ and tSZ) [@Sunyaev:1970] for the
	 same set of objects. 

     
While it is possible to create realistic maps of astrophysical effects through
 hydrodynamical simulations [@Dolag:2015], these methods are numerically
  expensive for large numbers of objects and reproducing them under different
   cosmologies and initial conditions can be prohibitive. An alternative
    strategy for creating mock observations of extended objects
  such as galaxies and galaxy cluster halos is to simulate the
   positions of these objects (either semi-analytically or through N-body
    simulations [@Stein:2020; @Stein:2018; @Sehgal:2010]) and then synthetically
     paint the desired signal at the location of the halos. `AstroPaint` is
      developed to help researchers in creating mock maps using the latter
       strategy. 

AstroPaint can also be used to create templates for detecting astrophysical
 effects in image data. For example, to detect kSZ for an ensemble of
  galaxies in a CMB map, one needs a
  template of this effect for the observed patch of the sky. Such a template can
   be generated by taking the catalog of the target galaxies along with their
    velocities and painting kSZ profiles around them on a map using
     `AstroPaint`. 
 
# Package Structure and Workflow 


`AstroPaint` consists of three main objects that interact with each other: `Catalog`, `Canvas`, and `Painter`. 


`Catalog` contains the locations, velocities, and masses of the objects. 
`Canvas` contains the map of the astrophysical signal in HEALPix format
 [@Healpy:2019]. 
`Painter` contains the template for the radial profile of the signal to be
 painetd on the `Canvas` in circular discs centered at the location of the
  halos in the
  `Catalog`.   

 These objects are sequentially passed into each other according to the
  following workflow: 

```python
from astropaint import Catalog, Canvas, Painter

catalog = Catalog(data=input_data)
canvas = Canvas(catalog, nside)
painter = Painter(template=radial_profile)

painter.spray(canvas)
```

The output map array can be accessed via `canvas.pixels` or directly
 visualized using `canvas.show_map()`. Here `input_data` is the dataframe that
  hold the locations, velocities, and
 masses of the halos. `nside` is a parameter in `healpy` that determines the
  total number of pixels (`npix = 12 * nside ** 2)` and
   consequently the resolution of the map . Finally, `radial_profile` is a one-dimensional function that determines the shape
   of the profile. A mind map visualization of the package structure can be
    found in [here](https://www.mindmeister.com/1417665103/astropaint-astropaint-py?fullscreen=1).   


# Acknowledgements

We would like to thank Simone Ferraro, Tim Morton, Vera Gluscevic, George Stein, Mathew Madhavacheril, 
Zack Li, and Alex van Engelen for their incredibly helpful comments and
 feedback. SY is grateful to the BCCP
  group at UC Berkeley for their
 hospitality during
 Summer 2019 where this project was inaugurated. EP is supported by NASA
  80NSSC18K0403 and the Simons Foundation award number 615662; EP and SY are supported by NSF AST-1910678.

# References
