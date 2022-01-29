![baniere.png](http://www.normalesup.org/~benmiquel/misc/baniere.png)

# Coral v1.1


Coral is a spectral PDEs solver for the plane-layer geometry, written in modern Fortran and highly scalable. 
The name comes from the initial motivation, which was **Co**nvection in **Ra**pidly rotating **L**ayers (eventhough the code has now a much more general scope).

See current version and history [here](https://github.com/BenMql/coral/wiki/versions)


## Tutorials
+ [Installation](https://github.com/BenMql/coral/wiki/1_Installation)
+ [First Run](https://github.com/BenMql/coral/wiki/2_First_run)
+ [Restarting a run](https://github.com/BenMql/coral/wiki/3_Restarting_a_run)
+ [Data vizualisation](https://github.com/BenMql/coral/wiki/4_Data_visualisation)

## Reference
+ [Input files--General overview](https://github.com/BenMql/coral/wiki/5_Input_files). 
+ [Input files--Equations](https://github.com/BenMql/coral/tree/master/doc): python notebooks with TeX formulas
+ [Common issues](https://github.com/BenMql/coral/wiki/6_Common_issues)

## Benchmarks
+ [Rayleigh-Benard](https://github.com/BenMql/coral/wiki/benchmark01_rayleigh_benard)
+ [Rotating Rayleigh-Benard](https://github.com/BenMql/coral/wiki/benchmark02_RRBC_julien_JFM96)
+ [Rotating convective dynamo](https://github.com/BenMql/coral/wiki/benchmark03_rotating_convective_MHD)

![coral.jpg](http://www.normalesup.org/~benmiquel/misc/coral.jpg)
# Contributing to Coral

Greetings, potential contributor!

We are excited to find you here, with perhaps the intent to contribute to this open-source project. We summarily describe in this file how you can make a very much appreciated contribution to the Fortran solver itself, the documentation, data-processing routines, etc.



## How to contribute

A good starting point is perhaps [Git's guide to making open-source contributions](https://opensource.guide/how-to-contribute/).

At the time of writing, Coral is still in its relative infancy, by which we mean that there is a long list of planned modifications
that range from unifying coding style and source files templates, to implementing new features. Thus before you put some work into 
any form of contribution, please create an issue that motivates (found a bug; thought of a new feature; etc.) and describes your proposed contribution.
It is my hope that in a near future contributions can be initiated in a more independent way. However, for the time being, this is how 
conflicts with planned development would be avoided.

After making sure that the proposed contribution fits in the planned roadmap, it will be included after the pull request has been reviewed.

#### Guidelines
Disclaimer: some of the sources of coral do not follow the guidelines below. The reason why is that the development started a while ago, 
when I did not reckon that the project would be shared with a community. Once in a while, I clean one of these ancient sources so that it
follows a better coding hygiene, but there remains some files with poor syntax style. Nonetheless, please follow a few guidelines for new contributions:
+ Fortran code: Use modern fortran standard (f90 to f2008) and syntax (mostly lower case; long, meaningful variable names; line breaks and indents for legibility). Avoid deprecated constructions. 
 Comment and document (using any syntax, or Doxygen if you know it). Debug and put your code through a grinder (e.g. use `gfortran -pedantic -std=f2003`, use multiple compilers, etc.).
+ Python: please use exclusively `python3`. I suggest following [these guidelines](https://www.python.org/dev/peps/pep-0008/)



## Copyrights
While you retain copyrights on your contributions, these contributions have to be licensed under the GNU Public Licence v3, similarly to the rest of Coral. 

---
title: 'Coral: a parallel spectral solver for fluid dynamics and partial differential equations'
tags:
  - fluid dynamics
  - PDE
  - fortran
  - chebyshev
  - plane layer
authors:
  - name: Benjamin Miquel
    orcid: 0000-0001-6283-0382
    affiliation: 1
affiliations:
 - name: Université Paris-Saclay, CEA, CNRS, Service de Physique de l’Etat Condensé, 91191 Gif-sur-Yvette,France
   index: 1
date: 07 January 2021
bibliography: paper.bib
---

# Summary

`Coral` is a fast, flexible, and efficient time-stepper for solving a large class of partial differential equations, at the core of which are the Navier-Stokes equations that govern fluid motions. Written in Fortran and employing the MPI standard for parallelization, the scalability of `Coral` allows the code to leverage the resources of high-performance computing infrastructures (up to hundreds of thousands of core, see @decomp2d), while running efficiently on laptops and workstations. Equations are entered by the user in the form of a plain text file following a simple and legible syntax. No coding proficiency in Fortran is required. This flexibility makes `Coral` suitable for both students and researchers with no coding experience.


# Statement of need

Natural and industrial flows exist in numerous different flavours, including homogeneous incompressible flows, shear flow, stably or unstably stratified flows, rotating flows, and flows of an electrically conducting fluid. These flows, however, have in common that they can be modelled by sets of (quadratic) advection-diffusion equations for the velocity, and possibly for the density, the temperature, the salinity, the magnetic field, etc. Hard-coding the sets of equations corresponding to each of these flow configurations is complex, time-consuming, and error-prone. These difficulties impede the development of new models. While `Coral` was initially motivated by the study of **Co**nvection in **Ra**pidly rotating **L**ayers, its scope has broadened and now encompasses solving homogeneous quadratic partial differential equations in a plane-layer geometry, i.e., a 3D domain with periodic boundary conditions along the two horizontal directions $x$ and $y$. Internally, `Coral` expands the variables along Fourier basis (horizontal directions) and Chebyshev polynomials (vertical direction). Transforms from physical to spectral space and domain decomposition are handled by the 2decomp&fft library [@decomp2d]. The quasi-inverse technique permits employing an arbitrarily large numbers of Chebyshev polynomials, resulting in the ability to resolve thin boundary layers characteristic of turbulent flows without suffering from loss of accuracy. Early versions of Coral have been used for studies concerning the turbulent motion of convective flows in presence of internal heat sources and sinks [@miquelPRF19; @miquelJFM20].

# Validation and examples

Coral has been validated on a variety of test cases (gathered in `etc/benchmarks`) found in the literature: Rayleigh-Bénard convection [@chandra], rotating convection [@julienJFM96], and convective dynamos [@stellmachPRE04; @cooperPRF20]. Those accuracy benchmarks, bound to grow in number, also constitute a library of examples for defining PDEs in Coral.

# State of the field

Among the existing flexible spectral solvers for marching in time PDEs in Cartesian geometries, alternatives to Coral include Dedalus [@dedalus], spectralDNS [@spectralDNS], FluidDyn [@fluiddyn], and FluidSim [@fluidsim]. For more complex geometries, options include `nek5000` [@nek5000], Nektar++ [@nektar], Freefem++ [@freefem], and Fenics [@fenics].

# Acknowledgements

The author warmly thanks Basile Gallet, Keith Julien, and Nick Featherstone for discussions and encouragement during the genesis of this project. This work was granted access to the HPC resources of TGCC and CINES under the allocation 2020-A0082A10803 attributed by GENCI (Grand Equipement National de Calcul Intensif). This work was supported in part by the National Science Foundation under Grants No. DMS-1317666, and No. NASA-NNX17AM01G.

# References
