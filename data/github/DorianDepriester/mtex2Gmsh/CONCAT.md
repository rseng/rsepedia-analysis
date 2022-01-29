# MTEX2Gmsh [![DOI](https://joss.theoj.org/papers/10.21105/joss.02094/status.svg)](https://doi.org/10.21105/joss.02094) [![View MTEX2Gmsh on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://fr.mathworks.com/matlabcentral/fileexchange/71469-mtex2gmsh) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



This toolbox for Matlab allows to generate meshes from EBSD data. It is intended to perform Finite Element Analysis (FEA) at grain scale on polycrystal imaged by EBSD. It is based on [MTEX](http://mtex-toolbox.github.io/) and [Gmsh](http://gmsh.info/).

## :thinking: How it works
This toolbox defines the class named `gmshGeo`. Once the grains are computed using MTEX, an instance of `gmshGeo` can be constructed. This object can be used to generate a Gmsh-readable file, in order to mesh it and perform FEA.

## :construction_worker: Requirements
This toolbox has been designed for MATLAB R2013b, but it may work on newer versions. In addition, the following are required:
- The [MTEX toolbox](https://mtex-toolbox.github.io/) (v 5.3.1 or newer) should be installed in your MATLAB session;
- The [Gmsh software](http://gmsh.info/) (v 4.7.1 or newer) should be installed on your computer (at least its binary should accessible).

It works on both Windows and Unix-like plateform (Linux and Mac OS).

<details><summary><b>:penguin: Linux users</b></summary>
When running the ``mesh`` command, you may stumble on the error below:

    /MATLAB/sys/os/glnxa64/libstdc++.so.6: version `GLIBCXX_3.4.21' not found (required by gmsh)
    
If so, instead of running 

    matlab
    
run

    LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libstdc++.so.6" matlab

</details>

## :mag: Documentation and examples
Here is an example of mesh obtained from the EBSD map called ``aachen`` in MTEX:
![aachen example](https://doriandepriester.github.io/MTEX2Gmsh/Examples/aachen.png)


Visit the corresponding [site](https://doriandepriester.github.io/MTEX2Gmsh/) to see other examples and [full documentation](https://doriandepriester.github.io/MTEX2Gmsh/html/index.html). Alternatively, you can check out the [``docs/Examples``](https://github.com/DorianDepriester/MTEX2Gmsh/tree/master/docs/Examples) folder.

## :gear: Unit test
The aforementioned examples can be easily reproduced. In addition, the reader can check out the reproductibility of minimal example on [Code Ocean](https://codeocean.com/capsule/8758800/tree/v2).

## :books: Reference
If you use this work, please cite the following paper:

> Depriester et al., (2020). MTEX2Gmsh: a tool for generating 2D meshes from EBSD data. *Journal of Open Source Software*, 5(52), 2094, https://doi.org/10.21105/joss.02094

In BibTeX, use the following entry:
````
@article{MTEX2Gmsh,
  doi = {10.21105/joss.02094},
  url = {https://doi.org/10.21105/joss.02094},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {52},
  pages = {2094},
  author = {Dorian Depriester and R\'egis Kubler},
  title = {{MTEX2Gmsh}: a tool for generating {2D} meshes from {EBSD} data},
  journal = {Journal of Open Source Software}
}
````

## :ambulance: Bug report
Please, use the [Issue](https://github.com/DorianDepriester/MTEX2Gmsh/issues) tab to report any bug or whish for new feature.

## :handshake: Contribute
You can easily edit the present code so that it fits your needs (as long as this edit complies with the MIT licence). You are also welcome to contribute. In this case, please read [``CONTRIBUTING.md``](CONTRIBUTING.md).
# Contribute to this work
MTEX2Gmsh is developped under the [MIT licence](LICENSE). Thus, you are welcome to contribute in any way.
This file gives the guidelines to whom may be interested in contributing. 

## Ressources
Before requesting help or trying to contribute, please first check at the following ressources:
- [Online documentation](https://doriandepriester.github.io/MTEX2Gmsh/html/index.html),
- [Issue tracker](/issues).

## Bugs
Use the Issue tracker for reporting a bug or requesting help. Please, keep in mind that I ([@DorianDepriester](https://github.com/DorianDepriester)) am currently the only 
maintainer of this toolbox. 
Thus, I may take a while for considering your issues, depending on external factors such as other works, vacations etc.

## Testing
If you make any change in the core code, run the examples scripts (in [``/docs/Examples``](/docs/Examples)) to check that everything runs fine with a large variety of geometries.

## Submit changes
If your contribution passes the aforementioned tests and if you want to share it, create a pull request (see [here](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request) for details).
Again, do not feel offended if I do not respond quickly.

## Code of conduct
The Contributor Covenant V2 applies here. This code of conduct is available at https://www.contributor-covenant.org/version/2/0/code_of_conduct/
---
title: 'MTEX2Gmsh: a tool for generating 2D meshes from EBSD data'
tags:
  - Polycrystals
  - FEM
  - Grain Boundary
  - EBSD
  - mesh
  - Matlab
authors:
  - name: Dorian Depriester
    orcid: 0000-0002-2881-8942
    affiliation: 1
  - name: Régis Kubler
    orcid: 0000-0001-7781-5855
    affiliation: 1
affiliations:
  - name: MSMP laboratory (EA 7350), Ecole Nationale Supérieure d'Arts et Métiers, 2 cours des Arts et Métiers - 13617 Aix-en-Provence, France
    index: 1
date: 13 december 2019
bibliography: paper.bib

---

# Summary
In material sciences applied to crystalline materials, such as metals or ceramics, the grain morphology (size and shape) and the crystallographic texture are of great importance for understanding the macroscopic behaviour of the materials.
Micromechanics of polycrystalline aggregates consists in evaluating the thermo-mechanical behaviour of the aggregates at their grain scale. If the investigated material is subjected to macroscopic deformation, the local strain can be obtained either experimentally, thanks to full-field measurement methods such as microgrid technique [@Allais:1994] or Digital Image Correlation (DIC) [@Hild:2002], or thanks to numerical simulation of the microstructure. The latter needs to take into account the mechanical heterogeneities (due to the different constituents) and the anisotropy of each phase, depending on its crystalline orientation.

Orientation Imaging Microscopy (OIM), usually made from Electron Backscatter Diffraction (EBSD), is now widely used as a characterization technique. Indeed, it is in great interest for investigating the grain morphology and local crystal orientations in crystalline materials. Raw EBSD data can be considered as matrices of measurements of crystallographic data: each dot contains information about the phase and its orientation at the corresponding position.

In order to perform Finite Element Analysis (FEA) on a polycrystal, one needs to first generate a mesh based on either EBSD or reconstructed grains. In this mesh, the Grain Boundaries (GBs) must be accurately described since they play an important role in the overall behaviour of aggregates. Indeed, it is known that GBs increase the energy of the materials. The interfacial energy between two adjacent grains due to their boundary depends, among other parameters, on their misorientation and on the surface normal of the boundary [@Priester:2012]. In addition, @Zhong:2017 mentioned that the GB curvature is one of the most important properties of a microstructure. For instance, the driving force for grain growth depends on the local curvature of the GBs.

@Latypov:2016 proposed a program to generate regular pseudo-3D mesh, consisting in brick elements with only one element in thickness. Nevertheless, this program results in serrated descriptions of the GBs because of the regular structure of EBSD data. In addition, the element size must be constant, possibly resulting in a huge number of elements, depending on the size and the spatial resolution of the orientation map. @Dancette:2016 proposed the following method to generate a conforming mesh with smooth GBs:

* computation of the GBs based on a proper criterion;
* grain reconstruction using a graph theory-based method;
* spline interpolation of the GBs;
* meshing.

The criterion used by the previous authors for defining the GBs, called weight in the context of graph theory, was specially designed for cubic phases. The geometry was meshed using the Gmsh software [@Geuzaine:2009].

As a conclusion, it appears that no existing tool for generating meshes from EBSD data is able to provide a robust grain description (e.g. suitable for any kind of phase and geometry) together with customizable features (e.g. variable element sizes). The proposed software, named [`MTEX2Gmsh`](https://github.com/DorianDepriester/mtex2Gmsh) works regardless the number of phases and the symmetries of those phases. In addition, it provides a smooth and accurate definition of the GBs. It is based on the MTEX toolbox for Matlab [@Bachmann:2011] and the Gmsh software. Figure 1 schematically illustrates the proposed algorithm. [`MTEX2Gmsh`](https://github.com/DorianDepriester/mtex2Gmsh) allows to mesh the volume with a couple of options, such as:

* increasing element size with increasing distance from the grains boundaries;
* element type (tetrahedron, wedge or brick elements);
* nesting the Region of Interest (ROI) into a larger medium.

This sofware comes with an Abaqus plugin for importing the mesh and allocating the phase and Euler Angles of each grain.


![Schematic representation of the algorithm used in `MTEX2Gmsh`: 1) once the grains are reconstructed using to `MTEX` [@Bachmann:2011], the algorithm fetches all triple junctions (TJ) in the whole map; 2) each grain boundary is divided into TJ-to-TJ segments; 3) all those segments are smoothed using Bspline approximation; 3) this decriptions of the grains can be converted into Gmsh-readable files [@Geuzaine:2009], allowing to mesh the whole region efficiently. The Bspline approximation results in very accurate definitions of the GBs, with limited serration (usually introduced by the EBSD resolution) and limited number of elements.](GraphicalAbstract.png)


# References
# Running PRISMS-plasticity simulations from EBSD: an example
The MATLAB script named [Copper.m](Copper.m) illustrates the steps by step procedure used for converting EBSD into data suitable for PRISMS-Plasticity.

## Generate input data for PRISMS-Plasticity
All you have to do is running this script file to generate the files needed for running the crystal plasticity simulation, namely the mesh file and the orienation file.

## Running PRISMS-Plasticity
The PRISMS-Plasticity parameters are already set in [``prm.prm``](prm.prm) for simulating a tensile test along the x direction, up to 3% elongation. Just run
   
    path/to/prisms_binary prm.prm
 
 It is higly advised to run this simulation on multiple threads in order to speed it up. This can be done with MPI, e.g.:
 
     mpirun -np 8 path/to/prisms_binary prm.prm
     
 ---
 **NOTE**
 
*This simulation takes about 6 hours on 32 threads on Intel Xeon Gold 6242 CPU @ 2.80GHz.*
 
 ---
 
## Results
Here is an image of the grains we get from MTEX:

<img src="Grains.png" width="400">

and the figure below illustrates the equivalent strain at the end of tensile test:

<img src="Eq_strain.jpeg" width="500">
## Purpose of this toolbox
In order to evaluate the thermo-mechanical behaviour of crystalline materials (such as metals or ceramics) at microscopic scale, one usually perform numerical simulation at grain scale using the Finite Element Method. In order to proceed, one must first create a mesh which is representative of the real material.

The microstructure of crystalline materials is usually made from Electron Backscattered Diffraction (EBSD) technique. Thus, this toolbox is designed to generate meshes from EBSD in a robust and accurate way.

## Examples
[![Example: aachen.m](./Examples/aachen.png)](./Examples/aachen.png)
[``aachen.m``](Examples/aachen.m)


[![Example: titanium_medium.m](./Examples/titanium_medium.png)](./Examples/titanium_medium.png)
[``titanium_medium.m``](Examples/titanium_medium.m)

[![Example: twins.m](./Examples/twins.png)](./Examples/twins.png)
[``twins.m``](Examples/twins.m)

## Documentation
### Online
You can navigate the documentation [here](html/index.html).

### From MATLAB
#### Full documentation
Once the toolbox is installed on your Matlab session, open the documentation of the present toolbox by typing:

    doc
    
Then, click on "MTEX2Gmsh toolbox", under the _Supplemental Software_ section (bottom right).

#### Help functions
The ``gmshGeo`` class is the core of this toolbox. For comprehensive details about it, just type

    help gmshGeo
    
The following command will print all the ``GmshGeo`` methods:

    methods gmshGeo
    
For details about a given method (let say ``plot``):

    help gmshGeo/plot

## Reference
For further details, check out the corresponding paper [[1]](#1). If you use this project, please cite it as follows:

````
@article{MTEX2Gmsh,
  doi = {10.21105/joss.02094},
  url = {https://doi.org/10.21105/joss.02094},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {52},
  pages = {2094},
  author = {Dorian Depriester and R\'egis Kubler},
  title = {MTEX2Gmsh: a tool for generating 2D meshes from EBSD data},
  journal = {Journal of Open Source Software}
}
````

<a id="1">[1]</a> Depriester et al., (2020). MTEX2Gmsh: a tool for generating 2D meshes from EBSD data. *Journal of Open Source Software*, 5(52), 2094, https://doi.org/10.21105/joss.02094
# MTEX2Abaqus
Abaqus Plugin to import a geometry generated from MTEXGmsh and assign local properties (phase and mean orientation of each grain).

## Installation
In your folder containing Abaqus plugins (eg: \USER\user\abaqus_plugin\ by default on Windows), just drag and drop the parent folder (named MTEX2Abaqus)
