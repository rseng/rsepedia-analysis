# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
 advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
 address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
 professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at matador@ml-evs.science. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
---
title: 'matador: a Python library for analysing, curating and performing high-throughput density-functional theory calculations'
tags:
  - Python
  - density-functional theory
  - ab initio
  - crystal structure prediction
  - materials discovery
  - databases
  - castep
  - quantum espresso
  - mongodb
authors:
  - name: Matthew L. Evans
    orcid:  0000-0002-1182-9098 
    affiliation: 1
  - name: Andrew J. Morris
    orcid: 0000-0001-7453-5698
    affiliation: 2
affiliations:
 - name: Theory of Condensed Matter Group, Cavendish Laboratory, University of Cambridge, J. J. Thomson Avenue, Cambridge, CB3 0HE, U.K.
   index: 1
 - name: School of Metallurgy and Materials, University of Birmingham, Edgbaston, Birmingham, B15 2TT, U.K.
   index: 2

date: July 2020
bibliography: paper.bib
---

# Summary

The properties of materials depend heavily on their atomistic structure; knowledge of the possible stable atomic configurations that define a material is required to understand the performance of many technologically and ecologically relevant devices, such as those used for energy storage [@NaP; @CuP]. First-principles crystal structure prediction (CSP) is the art of finding these stable configurations using only quantum mechanics [@JM]. Density-functional theory (DFT) is a ubiquitous theoretical framework for finding approximate solutions to quantum mechanics; calculations using a modern DFT package are sufficiently robust and accurate that insight into real materials can be readily obtained. The computationally intensive work is performed by well-established, low-level software packages, such as CASTEP [@castep] or Quantum Espresso [@qe], which are able to make use of modern high-performance computers. In order to use these codes easily, reliably and reproducibly, many high-level libraries have been developed to create, curate and manipulate the calculations from these low-level workhorses; `matador` is one such framework.

# Statement of need

The purpose of `matador` is fourfold:

- to promote the use of local databases and high-throughput workflows to increase the reproducibility of the computational results, 
- to perform reliable analysis of the stability, structure and properties of materials derived from calculations, 
- to provide tools to create customisable, publication-quality plots of phase diagrams, spectral properties and electrochemistry,
- to make the above functionality available to those with limited programming experience.

# `matador`

`matador` is a Python 3.6+ library and set of command-line tools for performing and analysing high-throughput DFT calculations using the CASTEP [@castep] and Quantum Espresso [@qe] packages. It is well-tested and fully-documented at [ReadTheDocs](https://matador-db.readthedocs.io), and comes with several tutorials and examples. The package is available on PyPI under the name [`matador-db`](https://pypi.org/project/matador-db). As with many projects, `matador` is built on top of the scientific Python ecosystem of NumPy [@numpy], SciPy [@scipy] and matplotlib [@matplotlib].

`matador` has been developed with high-throughput CSP in mind and has found use in the application of CSP to energy storage materials [@NaP; @CuP]; in this use case, a single compositional phase diagram can consist of tens of thousands of structural relaxation calculations. This package is aimed at users of CASTEP or Quantum Espresso who are comfortable with the command-line, yet maybe lack the Python knowledge required to start from scratch with more sophisticated packages. There are many mature packages that provide overlapping functionality with `matador`, the most widespread of which being the Atomic Simulation Environment (ASE) [@ase] and pymatgen [@pymatgen]. A translation layer to and from the structure representation of both of these packages is provided, such that analysis can be reused and combined. 

![Li-Zn-P ternary phase diagram created with matador, plot generated with matplotlib [@matplotlib] and python-ternary [@ternary].\label{fig:hull}](./LiZnP_hull.pdf)

# Overview of functionality

There are two ways of working with `matador`, either from the command-line interface (CLI) or through the Python library directly, with some features that are unique to each. The functionality of `matador` can be broadly split into three categories: 

1. *Creation and curation of databases of the results of first-principles calculations.*

`matador` allows for the creation of [MongoDB](https://mongodb.com) databases of CASTEP (6.0+) geometry optimisations from the command-line, using `matador import`. Only calculations that are deemed "successful", and that have fully-specified inputs are stored, with errors displayed for the rest. The resulting database can be queried with `matador query`, either with Python or through the powerful CLI. The results can be filtered for structural "uniqueness" and written to one of several supported file types and exported for use in other frameworks, such as ASE or pymatgen. Prototyping of structures directly from the database is achieved using `matador swaps`, which uses the same interface as `matador query` to return structure files with "swapped" elements [@NaP].

2. *High-throughput calculations and automated workflows.*

The `run3` executable bundled with `matador` allows for high-throughput calculations to be performed with little setup and no programming knowledge. Specialised support for CASTEP and the post-processing tool OptaDOS [@optados; @optados2] is provided to perform high-throughput geometry optimisations, orbital-projected band structures and densities of states, phonon calculations and elastic properties, however `run3` can also be used to run generic MPI programs concurrently on a set of structures. Sensible defaults for these workflows are provided by leveraging the open-source SeeK-path [@seekpath] and spglib [@spglib] libraries. The bundled `dispersion` script and associated library functionality allows for the creation of publication-quality spectral and vibrational property plots, in a similar fashion to the `sumo` package [@sumo]. The `matador.compute` module behind `run3` also powers the `ilustrado` genetic algorithm code [@ilustrado].

3. *Stability and structural analysis (with an emphasis on battery materials).*

The construction of reliable compositional phase diagrams requires several independent calculations to be performed on different atomic configurations with a compatible set of external parameters. These can be generated from a database query using `matador hull`, which allows the user to filter between different sets of calculations, and, where relevant, `matador voltage` can provide the electrochemical properties of that same phase diagram. Structural fingerprints implemented include pair distribution functions, powder X-ray diffraction patterns, and periodic crystal bond graphs. As more calculations are performed, changes to phase diagrams stored in the local database can be tracked with `matador hulldiff`. Phase diagrams can also be constructed from multiple energy values per structure, for example to show the effects of finite temperature [@CuP], or in the specific case of ensemble-based exchange-correlation functionals like the Bayesian Error Estimate Functional (BEEF) [@beef]. An example of a ternary phase diagram is shown \autoref{fig:hull}.

# Acknowledgements

We acknowledge all the contributors, users and testers of this package, primarily Angela Harper, James Darby, Jordan Dorrell and Matthew Cliffe. M.E. would like to acknowledge the EPSRC Centre for Doctoral Training in Computational Methods for Materials Science for funding under grant number EP/L015552/1. A.J.M. acknowledges funding from EPSRC (EP/P003532/1). The authors acknowledge networking support via the EPSRC Collaborative Computational Projects, CCP9 (EP/M022595/1) and CCP-NC (EP/T026642/1). Much of the development and testing was performed on the Cambridge Service for Data Driven Discovery (CSD3) operated by the University of Cambridge Research Computing Service (http://www.csd3.cam.ac.uk/), provided by Dell EMC and Intel using Tier-2 funding from the Engineering and Physical Sciences Research Council, and DiRAC funding from the Science and Technology Facilities Council (www.dirac.ac.uk).
Date: 13:41 16/04/2020  
Command: matador voltage -int -c NaSnP --db me388_NaPSn --pathways -hc 0.0 --res  
Version: 0.8b1-772-g94db115-voltages-and-volumes  

--------------------------------------------------------------------------------------------------------------------------
```
                  Root                   Pressure  Volume/fu  Hull dist.   Space group     Formula      # fu   Prov.  
                                          (GPa)     (Ang^3)   (meV/atom)  
--------------------------------------------------------------------------------------------------------------------------
Na-OQMD_8535-CollCode44758                   0.04       36.1          0.0  P6_3/mmc         Na          2      ICSD  
Na15Sn4-OQMD_90363-CollCode105167           -0.00      560.7          0.0    I-43d        Na15Sn4       2      ICSD  
NaP-CollCode182164                           0.02       93.9          0.0   P6_3cm         Na3P         6      ICSD  
NaSn-CollCode167672                         -0.01      283.7          0.0   P3_212        Na7Sn3        6      ICSD  
NaPSn-GA-dw2swh-1x103                       -0.03      320.0          0.0     P-1         Na8P4Sn       2       GA   
NaP-CollCode56530                           -0.02      196.1          0.0    C2/m          Na5P4        1      ICSD  
NaPSn-GA-dw2swh-1x79                        -0.06      208.5          0.0     P-1         Na5P3Sn       4       GA   
NaP-CollCode421420                          -0.04       41.9          0.0   P2_1/c          NaP         8      ICSD  
NaSn-CollCode641294                         -0.01       59.5          0.0  I4_1/acd        NaSn         16     ICSD  
NaPSn-GA-npf6il-1x3                         -0.04       78.2          0.0     P1           NaPSn        2       GA   
NaPSn-GA-npf6il-2x47                        -0.04       78.2          0.0     P1           NaPSn        2       GA   
NaP-CollCode60774                           -0.03      251.5          0.0 P2_12_12_1       Na3P7        4      ICSD  
NaSn-NaSn2-Baggeto                           0.02       80.1          0.0    Cmmm          NaSn2        2     SWAPS  
NaP-na3p11collo                              0.01      326.1          0.0    Pbcn         Na3P11        4      ICSD  
NaSn-146907-4243-58                         -0.02      105.6          0.0    R-3m          NaSn3        4     SWAPS  
NaP-NaP7                                     0.03      192.1          0.0  I4_1/acd        NaP7         8     SWAPS  
P-OQMD_5708-CollCode29273                    0.02       23.6          0.0    P2/c            P          84     ICSD  
P3Sn-OQMD_3387-CollCode16293                 0.00       85.9          0.0    R-3m          P3Sn         2      ICSD  
P3Sn4-OQMD_2953-CollCode15014                0.01      165.8          0.0    R-3m          P3Sn4        1      ICSD  
Sn-Collo                                     0.03       28.2          0.0  I4_1/amd         Sn          2      ICSD  
```Date: 11:58 23/04/2020  
Command: matador voltage -c LiCoP -int --db oqmd --spin any --volume --res  
Version: 0.8b1-772-g94db115-voltages-and-volumes  

--------------------------------------------------------------------------------------------------------------------------
```
                  Root                   Pressure  Volume/fu  Hull dist.   Space group     Formula      # fu   Prov.  
                                          (GPa)     (Ang^3)   (meV/atom)  
--------------------------------------------------------------------------------------------------------------------------
OQMD 5154                                    0.00       55.8          0.0   P63/mmc        Li3P         2      ICSD  
OQMD 12657                                   0.00      122.3          0.0    Amm2         Co6LiP4       1      ICSD  
OQMD 8021                                    0.00       22.8          0.0    Pnma           CoP         4      ICSD  
OQMD 6978                                    0.00       38.9          0.0    P21/c         CoP2         4      ICSD  
OQMD 4397                                    0.00       56.9          0.0    Im-3          CoP3         4      ICSD  
OQMD 17007                                   0.00       30.3          0.0    P21/c          LiP         8      ICSD  
OQMD 4845                                    0.00       31.7          0.0    Pnma          Co2P         4      ICSD  
OQMD 11352                                   0.00      192.0          0.0   P212121        Li3P7        4      ICSD  
OQMD 4371                                    0.00      164.2          0.0   I41/acd        LiP7         8      ICSD  
OQMD 24584                                   0.00       23.8          0.0     P-1            P          42     ICSD  
OQMD 594128                                  0.00       10.8          0.0    Fm-3m          Co          1      ICSD  
OQMD 30676                                   0.00       18.2          0.0    R-3m           Li          3      ICSD  
```---
name: Issue template
about: Default issue template
title: ''
labels: ''
assignees: ''

---

If you think you have found a bug, please raise an issue on GitHub, providing information about what you were trying to do, the function/script you ran, the error message/output and your matador version. If you are able to, please try to replicate the problem on the master branch before posting.

If you have a feature request, or you would like to contribute a feature, please raise an issue so that its suitability for this package can be discussed. Ideally, you could do this before writing your implementation, as the feature may already exist! Code contributions must include tests for any new functionality and must pass existing tests and flake8 linting run by the CI. Any new modules must use the black auto-formatter.
# Examples

This folder contains an assortment of examples, some as Jupyter notebooks and
some as tutorials to be followed according to the text files within. 

A full explanation, with links to the appropriate Binder can be found in the
[Examples section](https://docs.matador.science/en/stable/examples_index.html)
of the online documentation.

Examples in the `./interactive/` folder can be run without any external data,
either locally or with Binder. Those in the `./non-interactive/` folder require
external data but provide useful code snippets. 
