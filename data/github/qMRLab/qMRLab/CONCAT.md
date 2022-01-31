---
title: 'qMRLab: Quantitative MRI analysis, under one umbrella'
tags:
  - Matlab
  - Octave
  - quantitative magnetic resonance imaging
  - mri
  - neuroimaging
authors:
  - name: Agah Karakuzu
    orcid: 0000-0001-7283-271X
    affiliation: "1, 4"
  - name: Mathieu Boudreau
    orcid: 0000-0002-7726-4456
    affiliation: "1, 4"
  - name: Tanguy Duval
    orcid: 0000-0002-1228-5192
    affiliation: 1
  - name: Tommy Boshkovski
    orcid: 0000-0002-5243-5204
    affiliation: 1
  - name: Ilana R. Leppert
    affiliation: 2
  - name: Jean-François Cabana
    orcid: 0000-0003-0579-5378
    affiliation: 6
  - name: Ian Gagnon
    orcid: 0000-0001-6815-504X
    affiliation: 1
  - name: Pascale Beliveau
    orcid: 0000-0002-1971-4877
    affiliation: 7
  - name: G. Bruce Pike
    orcid: 0000-0001-8924-683X
    affiliation: "2, 5"
  - name: Julien Cohen-Adad
    orcid: 0000-0003-3662-9532
    affiliation: "1, 3"
  - name: Nikola Stikov
    orcid: 0000-0002-8480-5230
    affiliation: "1, 4"
affiliations:
 - name: NeuroPoly Lab, Institute of Biomedical Engineering, Polytechnique Montreal, Montreal, Canada
   index: 1
 - name: McConnell Brain Imaging Center, Montreal Neurological Institute, McGill University, Montreal, Canada
   index: 2
 - name: Functional Neuroimaging Unit, CRIUGM, University of Montreal, Montreal, Canada
   index: 3
 - name: Montreal Heart Institute, University of Montréal, Montréal, Canada
   index: 4
 - name: Departments of Radiology and Clinical Neuroscience,  Hotchkiss Brain Institute, University of Calgary, Calgary, Canada
   index: 5
 - name: Chaudière-Appalaches Integrated Health and Social Services Center, Sainte-Marie, Canada
   index: 6
 - name: Department of radio-oncology, Research Center of the University of Montreal Hospital Center, Montreal, Canada
   index: 7
date: 02 March 2020
bibliography: paper.bib
---

# Summary

Magnetic resonance imaging (MRI) has revolutionized the way we look at the human body. However, conventional MR scanners are not measurement devices. They produce digital images represented by “shades of grey”, and the intensity of the shades depends on the way the images are acquired. This is why it is difficult to compare images acquired at different clinical sites, limiting the diagnostic, prognostic, and scientific potential of the technology.

Quantitative MRI (qMRI) aims to overcome this problem by assigning units to MR images, ensuring that the values represent a measurable quantity that can be reproduced within and across sites. While the vision for quantitative MRI is to overcome site-dependent variations, this is still a challenge due to variability in the hardware and software used by MR vendors to produce quantitative MRI maps.

Although qMRI has yet to enter mainstream clinical use, imaging scientists see great promise in the technique's potential to characterize tissue microstructure. However, most qMRI tools for fundamental research are developed in-house and are difficult to port across sites, which in turn hampers their standardization, reproducibility, and widespread adoption.

To tackle this problem, we developed qMRLab, an open-source software package that provides a wide selection of qMRI methods for data fitting, simulation and protocol optimization \autoref{fig:header}. It not only brings qMRI under one umbrella, but also facilitates its use through documentation that features online executable notebooks, a user friendly graphical user interface (GUI), interactive tutorials and blog posts.

![qMRLab is an open-source software for quantitative MRI analysis It provides a myriad of methods to characterize microstructural tissue properties, from relaxometry to magnetization transfer.\label{fig:header}](https://github.com/qMRLab/qMRLab/raw/master/docs/logo/page_header.png)

MATLAB is the native development language of qMRLab, primarily because it is by far the most common choice among MRI methods developers. However, we have made a strong effort to lower licensing and accessibility barriers by supporting Octave compatibility and Docker containerization.

qMRLab started as a spin-off project of qMTLab [@cabana:2015]. In the meantime, a few other open-source software packages were developed, addressing  the lack of qMRI consistency  from different angles. QUIT [@wood:2018] implemented an array of qMRI methods in C++, which is highly favorable as an on-site solution because of its speed. The hMRI toolbox [@tabelow:2019] was developed as an SPM [@ashburner:2012] module that expands on the multi-parametric mapping method [@weiskopf:2008]. Other tools such as mrQ [@mezer:2016] and QMAP [@samsonov:2011] are also primarily designed for brain imaging. In addition to the arrival of community-developed tools for multi-modal qMRI processing, the field of diffusion MRI (dMRI) has recently witnessed an increase in the development of open-source software, bringing more transparency to the diffusion-based microstructural characterization of the brain. For example, the Dmipy Toolbox implemented an array of multi-compartment diffusion models in Python [@fick:2019], another Python package DIPY brought together dMRI processing methods at multiple levels [@garyfallidis:2014], Camino was developed in Java to provide users with yet another dMRI pipeline [@cook:2006] and TractoFlow [@theaud:2020] introduced container-based reproducible dMRI workflows developed in Nextflow pipeline orchestration tool [@di:2017]. Yet, brain imaging is not the only qMRI area slowed down by lack of consistency. Recently we published a preprint demonstrating notable disagreements between cardiac qMRI methods [@hafyane:2018]. Open-source software can go a long way in explaining these discrepancies, and the cardiac imaging community was recently introduced to TOMATO [@werys:2020], an open C++ framework for parametric cardiac MRI.

As open-source practices in the realm of qMRI become more popular, the need for effective communication of these tools also increases. This is important not only because we need consistency and transparency in the implementations, but also because non-specialist qMRI users would benefit from better understanding of the methodology. To this end, we envision qMRLab as a powerful tool with which users can easily interact with various techniques, perform simulations, design their experiments and fit their data. We reinforce this vision through our web portal (https://qmrlab.org) that includes interactive tutorials, blog posts and Jupyter Notebooks running on BinderHub, all tailored to a wide range of qMRI methods. The qMRLab portal is open for community contributions.

Currently, qMRLab is used by dozens of research labs around the world, mostly, but not limited to, application in brain and spinal cord imaging. A list of published studies using qMRLab is available on our GitHub repository.

While closed solutions may be sufficient for qualitative MRI (shades of grey lack standardized units), quantitative MRI will not realize its potential if we cannot peek inside the black box that generates the numbers. With qMRLab we want to open the black boxes developed in-house and reach a critical mass of users across all MR vendor platforms, while also encouraging developers to contribute to a central repository where all features and bugs are in the open. We hope that this concept will level the field for MR quantification and open the door to vendor-neutrality. We've been sitting in our MR cathedrals long enough. It is now time to join the MR bazaar [@raymond:1999]!

# Acknowledgements

This research was undertaken thanks, in part, to funding from the Canada First Research Excellence Fund through the TransMedTech Institute. The work is also funded in part by the Montreal Heart Institute Foundation, Canadian Open Neuroscience Platform (Brain Canada PSG), Quebec Bio-imaging Network (NS, 8436-0501 and JCA, 5886, 35450), Natural Sciences and Engineering Research Council of Canada (NS, 2016-06774 and JCA, RGPIN-2019-07244), Fonds de Recherche du Québec (JCA, 2015-PR-182754), Fonds de Recherche du Québec - Santé (NS, FRSQ-36759, FRSQ-35250 and JCA, 28826), Canadian Institute of Health Research (JCA, FDN-143263 and GBP, FDN-332796), Canada Research Chair in Quantitative Magnetic Resonance Imaging (950-230815), CAIP Chair in Health Brain Aging, Courtois NeuroMod project and International Society for Magnetic Resonance in Medicine (ISMRM Research Exchange Grant).

# References
# Changelog
All notable changes to this project will be documented in this file.

## Release [2.4.1] - 2020-09-02

## New ✨
- 🆕 model: `inversion_recovery` 
    - Add general equation fitting in addition to Barral's model.

### Improvements 🚀
- GUI (JOSS review by @mfroeling)
    - Please see changes [here](https://github.com/qMRLab/qMRLab/pull/400).
- Documentation (JOSS review by @grlee77)
    - Please see changes [here](https://github.com/qMRLab/qMRLab/pull/399)

### Bug Fixes🐛
- `FilterClas` bug [fix](https://github.com/qMRLab/qMRLab/pull/385).

### Other
- Change citation reference to JOSS paper
    - Karakuzu A., Boudreau M., Duval T.,Boshkovski T., Leppert I.R., Cabana J.F., 
    Gagnon I., Beliveau P., Pike G.B., Cohen-Adad J., Stikov N. (2020), qMRLab: 
    Quantitative MRI analysis, under one umbrella doi: 10.21105/joss.02343

## Release [2.4.0] - 2020-02-14

### New ✨
- 🆕 model: `mp2rage` 
    - Fit MP2RAGE data to create a T1map.
    - The original codebase is [here](https://github.com/JosePMarques/MP2RAGE-related-scripts).
    - Check out [qMRLab's MP2RAGE blog post](https://qmrlab.org/2019/04/08/T1-mapping-mp2rage.html) by @mathieuboudreau!
- 🆕 model: `mono_t2`
    - Fit MESE data to create a T2map.
- 🆕 simulator: `Monte-Carlo Diffusion`
    - Monte Carlo simulator for 2D diffusion is able to generate synthetic 
    diffusion signal from any 2D axon packing.
    - An MRathon project by @Yasuhik, @TomMingasson and @tanguyduval. 
- 🆕 Changelog ❤️

### Improvements 🚀
- Model: `qsm_sb` 
    - With the new echo combination implementation, `qsm_sb` can now take 
      multi-echo GRE data. 
    - An MRathon project by @jeremie-fouquet.
- Get rid of redundant buttons in GUI `Protocol` panel. 

### Bug Fixes🐛
- `qMRgenBatch` account for models w/o fixed required inputs (e.g. `mp2rage`).
- Remove old built packages from `qmrlab/mcrgui`.
- Fix `qmrlab/octjn` dependencies.

### Removed 🧹

## Release [2.3.1] - 2020-01-07

### New ✨
- 🆕 static member function: getProvenance 
    - Scrape details and add more (optional) to save sidecar `*.json` files for maps.
    - See an example use [here](https://github.com/qMRLab/qMRWrappers/blob/master/mt_sat/mt_sat_wrapper.m).
- 🆕 Docker image: `qmrlab/minimal`
    - qMRLab + Octave - Jupyter for [qMRFlow](https://github.com/qMRLab/qMRflow) pipelines.    

### Improvements 🚀
- New MATLAB/Octave env: `ISNEXTFLOW` 
    - Deals with the `load_nii` case for symlinked inputs.
    - Enforces `gzip -d --force` if `ISNEXTFLOW` 
    - Commonly used by `qMRWrappers` 

### Bug Fixes🐛
- N/A

### Removed 🧹
- N/A 

## Release [2.3.0] - 2019-05-08

### New ✨

- 🆕 model: `Processing/filtermap` 
    - Apply 2D/3D spatial filtering, primarily intended for fieldmaps. 
        - `Polynomial`
        - `Gaussian` 
        - `Median` 
        - `Spline` 
- 🆕 model: `qsm_sb` 
    - Fast quantitative susceptibility mapping:
        - `Split-Bregman` 
        - `L1 Regularization`
        - `L2 Regulatization` 
        - `No Regularization` 
        - `SHARP background filtering` 
- 🆕 model: `mt_ratio` 
    - Semi-quantitative MTR. 
- 🆕 GUI 3D toolbox:
    - An array of UI tools for the visualization and brief statistical
      inspection of the data using ROI tools. 
- 🆕 functionality `qMRgenJNB`:
    - Create a Jupyter Notebook for any model. 
    - Insert Binder Badge to the documentation. 
- 🆕 Azure release pipelines and deployment protocols:
    - Set self-hosted Azure agent to compile qMRLab and ship in a Docker image
    - `qmrlab/mcrgui`: Use qMRLab GUI in a Docker image. 
    - `qmrlab/octjn`: Use qMRLab in Octave in Jupyter Env. 
    - See `/Deploy` folder for furhter details. 
    - [qMRLab DockerHub page.](https://hub.docker.com/orgs/qmrlab)

### Improvements 🚀
- Model: `vfa_t1`:
    - Bloch simulations are added 
    - Performance improvement 
- Model: `ir_t1` 
    - Parameter descriptions are improved. 
- Model: `b1_dam`
    - Protocol descriptions has been updated. 
- `FitTempResults`:
    - Is now saved every 5 minutes instead of every 20 voxels. 
    
### Bug Fixes🐛
- GUI fixes. 

### Removed 🧹
- N/A # [qMRLab](https://qmrlab.org)
[![Build Status](https://dev.azure.com/neuropoly/qMRLab/_apis/build/status/qMRLab?branchName=master)](https://dev.azure.com/neuropoly/qMRLab/_build/latest?definitionId=15&branchName=master) [![OCTAVE_CI_TESTS](https://github.com/qMRLab/qMRLab/actions/workflows/octaveTests.yml/badge.svg)](https://github.com/qMRLab/qMRLab/actions/workflows/octaveTests.yml) [![codecov](https://codecov.io/gh/qMRLab/qMRLab/branch/master/graph/badge.svg?token=mfvAuaJ8P2)](https://codecov.io/gh/qMRLab/qMRLab) [![Documentation Status](https://readthedocs.org/projects/pip/badge/?version=latest)](https://qmrlab.readthedocs.io/en/master/?version=latest) [![License](https://img.shields.io/github/license/mashape/apistatus.svg)](https://opensource.org/licenses/MIT) [![Website](https://img.shields.io/badge/Website-qmrlab.org-red.svg)](https://qmrlab.org) [![Compile](https://vsrm.dev.azure.com/neuropoly/_apis/public/Release/badge/46799180-12b7-4319-b8e7-94e1d39047e7/2/6)](https://dev.azure.com/neuropoly/qMRLab/_release?_a=releases&view=mine&definitionId=2) [![Twitter Follow](https://img.shields.io/twitter/follow/qmrlab.svg?style=social&label=Follow)](https://twitter.com/qmrlab)

![hdr](https://github.com/qMRLab/documentation/blob/master/logo/header_new.png?raw=true)

<b><h2 align="center">
qMRLab is an open-source software for quantitative MR image analysis.</h2></b> 

<p align="center">Our main goal is to provide the community with an intuitive tool for <b>data fitting</b>, <b>plotting</b>, <b>simulation</b> and <b>protocol optimization</b> for a myriad of different quantitative models.
The modularity of the qMRLab framework makes it easy to add any additional modules and we encourage everyone to contribute their favorite recipe for qMR!</p>

<p align="center">
For <b>documentation 📖</b>, visit the <a href="http://qmrlab.readthedocs.io/">Documentation website</a>.</p>

<p align="center">
If you are a <b>developer 🛠</b>, please visit the <a href="https://github.com/neuropoly/qMRLab/wiki">Wiki page<a>.</p>

<p align="center">
Please report any <b>bug 🐛</b> or <b>suggestions 💭</b> in <a href="https://github.com/neuropoly/qMRLab/issues">GitHub</a>.</p>

<p align="center">
For <b>interactive tutorials 🎚</b>, <b>blog posts 🖋</b> and more, you can visit <a href="https://qmrlab.org">qMRLab portal</a>.</p>

<p align="center">
qMRLab is a fork from the initial project <a href="https://github.com/neuropoly/qMTLab">qMTLab</a>.</p>  

***

<b><p align="center" style="font-size:24px">
FROM SCANNER</p></b> 

<p align="center">
Version controlled, fully transparent & vendor-neutral  <a href="https://github.com/qMRLab/pulse_sequences"> pulse sequences 🌀</a>.</p>  

<p align="center">
⬇️</p>  

<p align="center">
Data-driven, container-mediated, platform-agnostic & BIDS compatible  <a href="https://github.com/qMRLab/pulse_sequences"> workflows 🔀</a>.</p>

<p align="center">
⬇️</p>  

<p align="center">
 Reproducible & modern <a href="https://qmrlab.org/blog.html"> publication 📚</a> objects.</p>

<b><p align="center" style="font-size:24px">
TO PUBLICATION</p></b> 

***

## References

**qMRLab**

* [Karakuzu et al. *The qMRLab workflow: From acquisition to publication.* ISMRM 2019](https://www.ismrm.org/19/program_files/DP23.htm#005)
* [Duval et al. *Quantitative MRI made easy with qMRLab.* ISMRM 2018](http://archive.ismrm.org/2018/2288.html)
* [Cabana et al. *Quantitative magnetization transfer imaging made easy with qMTLab: Software for data simulation, analysis, and visualization.* Concepts in Magn Reson 2015](https://onlinelibrary.wiley.com/doi/abs/10.1002/cmr.a.21357)

**Applications**

* [Soustelle et al. *Correlations of quantitative MRI metrics with myelin basic protein (MBP) staining in a murine model of demyelination.* NMR in Biomed 2019](https://onlinelibrary.wiley.com/doi/abs/10.1002/nbm.4116)
* [Kim et al. *Rapid framework for quantitative magnetization transfer imaging with interslice magnetization transfer and dictionary‐driven fitting approaches.* Mag Res Med 2019](https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.27850)
* [Boudreau et al. *All you need is a browser: eliminating barriers to MRI education with open-source interactive tutorials.* Junior Fellows Symposium: Africa Challenge, ISMRM 2019](https://www.ismrm.org/19/program_files/W06.htm)
* [Romero and Sinha.  *Magnetization Transfer Saturation Imaging of Human Calf Muscle: Reproducibility and Sensitivity to Regional and Sex Differences.* J Magn Reson Imaging 2019](https://onlinelibrary.wiley.com/doi/abs/10.1002/jmri.26694)
* [Michálek et al.  *Fast and accurate compensation of signal offset for T<sub>2</sub> mapping.* Magn Reson Mater Phy 2019](https://link.springer.com/article/10.1007/s10334-019-00737-3)
* [Barbieri et al. *Circumventing the Curse of Dimensionality in Magnetic Resonance Fingerprinting through a Deep Learning Approach.* arXiv:1811.11477](https://arxiv.org/abs/1811.11477)
* [Varma et al. *Low duty-cycle pulsed irradiation reduces magnetization transfer and increases the inhomogeneous magnetization transfer effect.* J of Mag Res 2018](https://www.sciencedirect.com/science/article/abs/pii/S1090780718302088)
* [Campbell et al. *Promise and pitfalls of g-ratio estimation with MRI.* NeuroImage 2018](https://www.sciencedirect.com/science/article/pii/S1053811917306857)
* [Boudreau and Pike. *Sensitivity regularization of the Cramér‐Rao lower bound to minimize B1 nonuniformity effects in quantitative magnetization transfer imaging.* Mag Res Med 2018](https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.27337)
* [Mingasson et al. *AxonPacking: an open-source software to simulate arrangements of axons in white matter.* Frontiers in neuroinformatics 2017](https://www.frontiersin.org/articles/10.3389/fninf.2017.00005/full)

**Interactive Tutorials**

* [Boudreau M. *Relaxometry Series: MP2RAGE T<sub>1</sub> Mapping.* qMRLab.org 2019](https://qmrlab.org/2019/04/08/T1-mapping-mp2rage.html)
* [Boudreau M. *Relaxometry Series: Variable Flip Angle T<sub>1</sub> Mapping.* qMRLab.org 2018](https://qmrlab.org/jekyll/2018/12/11/T1-mapping-variable-flip-angle.html)
* [Boudreau M. *Relaxometry Series: Inversion Recovery T<sub>1</sub> Mapping.* qMRLab.org 2018](https://qmrlab.org/jekyll/2018/10/23/T1-mapping-inversion-recovery.html)

**Awards**

* Karakuzu A. ISMRM Research Exchange Grant 2020
* Karakuzu A. et al. Quantitative MR Study Group Competition, second place, ISMRM 2019 
* Boudreau M. et al. Junior Fellows Symposium Challenge, Africa challenge winner, ISMRM 2019
* Karakuzu A. et al. [Magnetic Moments](https://www.youtube.com/watch?v=67GKiK3iFr0), People's Choice Award, ISMRM 2018 

**Sponsored by**

* [Institut TransMedTech Montréal (ITMT)](https://www.polymtl.ca/transmedtech/)
* [Canadian Open Neuroscience Platform (CONP)](https://conp.ca/)
* [International Society of Magnetic Resonance in Medicine (ISMRM)](http://ismrm.org)

## Citation

If you use qMRLab in you work, please cite:

Karakuzu A., Boudreau M., Duval T.,Boshkovski T., Leppert I.R., Cabana J.F., Gagnon I., Beliveau P., Pike G.B., Cohen-Adad J., Stikov N. (2020), qMRLab: Quantitative MRI analysis, under one umbrella doi: 10.21105/joss.02343

[![OSF](https://img.shields.io/badge/DOI-10.21105%2Fjoss.02343-green.svg?style=for-the-badge)](https://doi.org/10.21105/joss.02343)

***

<p align="center">
The MIT License (MIT)</p>

<p align="center">Copyright (c) 2016 NeuroPoly, Ecole Polytechnique, Universite de Montreal</p>

<p align="center">Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:</p>

<p align="center">The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.</p>

<p align="center">THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
</p># Contributor Covenant Code of Conduct

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

We do not tolerate harassment or other, inappropriate behavior in our community.
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

This Code of Conduct applies within all project spaces, and it also applies when
an individual is representing the project or its community in public spaces.
Examples of representing a project or community include using an official
project e-mail address, posting via an official social media account, or acting
as an appointed representative at an online or offline event. Representation of
a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at agah.karakuzu@polymtl.ca. 

As the first interim [Benevolent Dictator for Life (BDFL)](https://en.wikipedia.org/wiki/Benevolent_dictator_for_life),
Agah Karakuzu can take any action she deems appropriate
for the safety of the `qMRLab` community, including but not limited to:

* facilitating a conversation between the two parties involved in the violation of the code of conduct
* requesting a contributor apologize for their behaviour
* asking a contributor or multiple contributors to enter a cooling off period that puts a
  time-limited pause on a particular discussion topic
* asking a contributor to no longer participate in the development of `qMRLab`.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html and from the Code of Conduct by [`ME-ICA/tedana`](https://github.com/ME-ICA/tedana/edit/master/CODE_OF_CONDUCT.md)

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
# Contributing to `qMRLab`
We love your input! We want to make contributing to this project as easy and transparent as possible, whether it's:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features
- Becoming a maintainer

## We Develop with Github
We use GitHub to host code, to track issues and feature requests, as well as accept pull requests.

## We Use [Github Flow](https://guides.github.com/introduction/flow/index.html), So All Code Changes Happen Through Pull Requests
Pull requests are the best way to propose changes to the codebase (we use [Github Flow](https://guides.github.com/introduction/flow/index.html)). We actively welcome your pull requests:

1. Fork the repo and create your branch from `master`.
2. If you've added code that should be tested, add tests.
3. If you've changed APIs, update the documentation.
4. Ensure the test suite passes.
5. Make sure your code lints.
6. Issue that pull request!

## Any contributions you make will be under the MIT Software License
In short, when you submit code changes, your submissions are understood to be under the same [MIT License](http://choosealicense.com/licenses/mit/) that covers the project. Feel free to contact the maintainers if that's a concern.

## Report bugs using Github's [issues](https://github.com/briandk/transcriptase-atom/issues)
We use GitHub issues to track public bugs. Report a bug by [opening a new issue](); it's that easy!

## Write bug reports with detail, background, and sample code
[This is an example](http://stackoverflow.com/q/12488905/180626) of a bug report written by Brian A. Danielak, and I think it's not a bad model. Here's [another example from Craig Hockenberry](http://www.openradar.me/11905408), an app developer whom I greatly respect.

**Great Bug Reports** tend to have:

- A quick summary and/or background
- Steps to reproduce
  - Be specific!
  - Give sample code if you can. [Brian A. Danielak's stackoverflow question](http://stackoverflow.com/q/12488905/180626) includes sample code that *anyone* with base requirements can run to reproduce the problem
- What you expected would happen
- What actually happens
- Notes (possibly including why you think this might be happening, or stuff you tried that didn't work)

People *love* thorough bug reports.

## Recognizing contributors
We welcome and recognize all contributions from documentation to testing to code development.
You can see a list of current contributors in the [contributors tab](https://github.com/qMRLab/qMRLab/graphs/contributors).

## License
By contributing, you agree that your contributions will be licensed under its MIT License.

## References
This document was adapted from the open-source contribution guidelines for [Facebook's Draft](https://github.com/facebook/draft-js/blob/a9316a723f9e918afde44dea68b5f9f39b7d9b00/CONTRIBUTING.md)## Purpose
_Describe the problem or feature in addition to a link to the issues._

## Approach
_How does this change address the problem?_

#### Open Questions and Pre-Merge TODOs
- [ ] Use github checklists. When solved, check the box and explain the answer.

- [ ] Review that changed source files/lines are related to the pull request/issue
_If any files/commits were accidentally included, cherry-pick them into another branch._

- [ ] Review that changed source files/lines were not accidentally deleted
_Fix appropriately if so._

- [ ] Test new features or bug fix
_If not implemented/resolved adequately, solve it or inform the developer by requesting changes in your review._
_Preferably, set breakpoints in the locations that the code was changed and follow allong line by line to see if the code behaves as intended._

##### Manual GUI tests (general)

- [ ] Does the qMRLab GUI open?
- [ ] Can you change models?
- [ ] Can you load a data folder for a model?
- [ ] Can you view data?
- [ ] Can you zoom in the image?
- [ ] Can you pan out of the image?
- [ ] Can you view the histogram of the data?
- [ ] Can you change the color map?
- [ ] Can you fit dataset (Fit data)?
- [ ] Can you save/load the results?
- [ ] Can you open the options panel?
- [ ] Can you change option parameters?
- [ ] Can you save/load option paramters?
- [ ] Can you select a voxel?
- [ ] Can you fit the data of that voxel ("View data fit")?
- [ ] Can you simulate and fit a voxel ("Single Voxel Curve")?
- [ ] Can you run a Sensitivity Analysis?
- [ ] Can you simulate a Multi Voxel Distribution?

##### Specifications

  - Version:
  - Platform:
  - Subsystem:
%%
%
% script_Laplacian_unwrap_Sharp_Fast_TV_gre3D.m : example reconstruction
% demonstrating closed-form L2-regularized and Split-Bregman L1-regularized
% QSM algorithms
%
%
% script_Laplacian_unwrap_Sharp_Magnitude_TV_gre3D.m : example
% reconstruction demonstrating magnitude weighted L2- and L1-regularized
% QSM algorithms

% MB Notes
% 
% Data and code generously provided by Dr. Berkin Bilgic
%
% References : 
% 
% Bilgic, B. , Fan, A. P., Polimeni, J. R., Cauley, S. F., Bianciardi, M. , Adalsteinsson, E. , Wald, L. L. and Setsompop, K. (2014), Fast quantitative susceptibility mapping with L1‐regularization and automatic parameter selection. Magn. Reson. Med., 72: 1444-1459. doi:10.1002/mrm.25029
% 
% https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.25029
% 
# AxonPacking: Simulator of White Matter Axons Arrangement

author : Tom Mingasson    
contact : mingasson.tom@gmail.com       
institution : University Polytechnique Montreal, NeuroPoly   
date : 2016 

<img src="https://github.com/neuropoly/axonpacking/blob/master/img1.jpeg" width="800px" align="center" />

## Description 

Here is  a new random disks packing algorithm for numerical simulation of the white matter. 

White matter tissue is divided in three compartments:  axons, myelin sheath and extra-axonal space. Axons are assumed to be parallel cylinders, therefore the invariance along the fiber axis makes it possible to consider this problem in 2D. The dense packing of axons is thus equivalent to the generation of random 2-dimensional packing of N perfectly round and non-compressible disks. Axon diameter distributions follow a Gamma distribution (defined by its mean µ and variance σ2). Interestingly the g-ratio is fairly constant across species and white matter regions (31,32) and is dependent mostly on the diameter of the axon according to the relationship presented in (Ikeda M, Oka Y. Brain Behav. 2012):  gratio= 0.220 * log(DIAMETER_unmyelinated) +0.508. 

The different steps to process packing are the following: first, the diameters of the disks are randomly chosen using a gamma or lognormal distribution parameterized with the mean (d_mean), variance (d_var) and number of axons (N).  Then, the positions of disks are initialized on a grid, and they migrate toward the center of the packing area until the maximum disk density is achieved. 


The software packing provides microstructure features (Fiber Volume Fraction FVF, Myelin Volume Fraction MVF, Axon Volume Fraction AVF, fraction restricted FR).

<img src="https://github.com/neuropoly/axonpacking/blob/master/img2.png" width="1000px" align="middle" />

## Scripts

- main.m
- axons_setup.m
- process_packing.m
- compute_gratio.m
- compute_statistics.m
- progressBar.m

## How to use it ?

### INPUTS
In ‘main.m’ change the inputs

- N : the number of disks i.e axons to include in the simulation
- d_mean and d_var : the diameter distribution parameters : mean µ and variance σ² of disk diameters 
- Delta : the fixed gap between the edge of disks Δ 
- iter_max : the number of iterations i.e disk migrations performed by the algorithm before computing the outputs 

#### Help 	
The disk density increases over the migrations and tends toward a limit value. It is necessary to first launch the algorithm with a high number of iterations iter_max. The disk density i.e FVF is calculated every 'iter_fvf' iterations to assess the sufficient number of iterations to reach convergence. 'iter_fvf' is a user defined integer: iter_fvf = iter_max/10 by default. 

When d_mean closed to 3 um, d_var  closed to 1 um : 
 - if N about 1000, iter_max = 30000 is sufficient. 
 - if N about 100, iter_max = 10000 is sufficient. 

#### Example  	
N = 100;            
d_mean = 3;         
d_var  = 1;        
Delta  = 0; 
iter_max = 10000;                            

### OUPUTS
The function ‘computeStatistics.m’ provides MVF, AVF, FVF and FR for each packing image defined by the input combinations. To do that it creates a binary mask.

Outputs are stored in 3 matlab structures. 

- in 'axons.mat' :  the axon features (N, d_mean, d_var, Delta, g_ratio and the drawn diameters d)
- in 'packing.mat' : the packing results (initial positions of disks (initial_positions) and final positions of disks (final_positions))
- in 'stats.mat' : the statistics results with the values for each metric computed in the packing (FVF, FR, MVF, AVF) 

A png image of the final packing with three different labels (intra-axonal, myelin and extra-axonal) is saved by default in the current script folder.
---
name: Basic template
about: Describe your problem using a basic template
title: ''
labels: ''
assignees: ''

---
## Expected Behavior


## Actual Behavior


## Steps to Reproduce the Problem

  1.
  1.
  1.

## Specifications

  - Version:
  - Platform:
  - Subsystem:---
name: Bug report - Q&A
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---
## *Who* is the bug affecting?
<!-- Ex. All supervisors, Sally Supervisor, Level 1 CCs -->

## *What* is affected by this bug?
<!-- Ex. supervision, sending messages, texter profiles -->

## *When* does this occur?
<!-- Ex. After ending a conversation, every night at 3pm, when I sign off -->

## *Where* on the platform does it happen?
<!-- Ex. In the a Supervisor chat box, on the conversation profile page, on the two-factor screen -->


## *How* do we replicate the issue?
<!-- Please be specific as possible. Use dashes (-) or numbers (1.) to create a list of steps -->


## Expected behavior (i.e. solution)
<!-- What should have happened? -->


## Other Comments---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

<!--- Provide a general summary of the issue in the Title above -->

## Expected Behavior
<!--- Tell us what should happen -->

## Current Behavior
<!--- Tell us what happens instead of the expected behavior -->

## Possible Solution
<!--- Not obligatory, but suggest a fix/reason for the bug, -->

## Steps to Reproduce
<!--- Provide a link to a live example, or an unambiguous set of steps to -->
<!--- reproduce this bug. Include code to reproduce, if relevant -->
1.
1.
1.
1.

## Context (Environment)
<!--- How has this issue affected you? What are you trying to accomplish? -->
<!--- Providing context helps us come up with a solution that is most useful in the real world -->

<!--- Provide a general summary of the issue in the Title above -->

## Detailed Description
<!--- Provide a detailed description of the change or addition you are proposing -->

## Possible Implementation
<!--- Not obligatory, but suggest an idea for implementing addition or change -->

# DICOM to NIfTI conversion, DICOM and NIfTI tools, NIfTI visualization (version 2019.06.24)

# dicm2nii
Convert DICOM into NIfTI. It can also convert PAR/XML/REC, HEAD/BRIK, MGZ and BrainVoyager files into NIfTI.

# nii_tool
Create, load, save NIfTI file. Support both version 1 and 2 NIfTI, and variety of data type.

# nii_viewer
Visualize NIfTI. Can also visualize any file convertible to NIfTI by dicm2nii.

# nii_moco
Perform motion correction on a NIfTI.

# nii_stc
Perform slice timing correction on a NIfTI.

# nii_xform
Transform a NIfTI into different resolution, or into a template space.

# dicm_hdr, dicm_img, dicm_dict
Read DICOM header and image, independent of Matlab Image Processing Toolbox. 

# rename_dicm, sort_dicm, anonymize_dicm
DICOM tools performing the tasks as indicated by the name.
# imtool3D
This is an image viewer designed to view a 3D stack of image slices. For example, if you load into matlab a DICOM series of CT or MRI images, you can visualize the images easily using this tool. It lets you scroll through slices, adjust the window and level, make ROI measurements, and export images into standard image formats (e.g., .png, .jpg, or .tif) or the 3D mask as NIFTI file (.nii). Use imtool3D_nii to load NIFTI files.
This tool is written using the object-oriented features of matlab. This means that you can treat the tool like any graphics object and it can easily be embedded into any figure. So if you're designing a GUI in which you need the user to visualize and scroll through image slices, you don't need to write all the code for that! Its already done in this tool! Just create an imtool3D object and put it in your GUI figure.

<p align="center">
  <img src="Capture.PNG" width="600">
</p>
  
imtool3D is used heavily by several other projects:
* [qMRLab](https://github.com/qMRLab/qMRLab)
* [imquest](https://gitlab.oit.duke.edu/railabs/SameiResearchGroup/imquest)
* [lesionTool](https://gitlab.oit.duke.edu/railabs/SameiResearchGroup/lesionTool)

# Dependencies
* Matlab's image processing toolbox (ROI tools are disabled otherwise)
* [dicm2nii](https://github.com/xiangruili/dicm2nii) (if NIFTI images are used)

# Tuto
## include in a GUI
````matlab
% Add viewer in a panel in the middle of the GUI
GUI = figure;
tool = imtool3D([],[.1 .1 .8 .8],GUI)

% set MRI image
load mri % example mri image provided by MATLAB
D = squeeze(D);
D = permute(D(end:-1:1,:,:),[2 1 3]); % LPI orientation
tool.setImage(D)
tool.setAspectRatio([1 1 2.5]) % set voxel size to 1mm x 1mm x 2.5mm

````

# what is new in this fork? 
* Support for 5D volumes (scroll through time and volumeS with arrows)
* Keyboard shortcut
* Multi-label mask
* Save mask
* NIFTI files (.nii) support (double click on a. nii file in Matlab filebrowser) 
* New tools for mask (interpolate slices, active contour...)
* Convert Mask2poly and poly2mask
* splines in polygons (double click a circle)
* 3 planes view

# Authors
Justin Solomon (Original release)  
Tanguy Duval (4D (time) and 5D (different contrast); multi-label mask; active_contour, undo button, mask2poly, poly2mask, shortcuts)  

# Original release
https://fr.mathworks.com/matlabcentral/fileexchange/40753-imtool3d
# BIDS for MATLAB / Octave

This repository aims at centralising MATLAB/Octave tools to interact with datasets conforming to the BIDS (Brain Imaging Data Structure) format.

For more information about BIDS, visit https://bids.neuroimaging.io/.

See also [PyBIDS](https://github.com/bids-standard/pybids) for Python and the [BIDS Starter Kit](https://github.com/bids-standard/bids-starter-kit).

## Implementation

Starting point was `spm_BIDS.m` from [SPM12](https://github.com/spm/spm12) ([documentation](https://en.wikibooks.org/wiki/SPM/BIDS#BIDS_parser_and_queries)) reformatted in a `+bids` package with dependencies to other SPM functions removed.

## Example

```Matlab
BIDS = bids.layout('/home/data/ds000117');
bids.query(BIDS, 'subjects')
```# AMICO

Implementation of the linear framework for Accelerated Microstructure Imaging via Convex Optimization (AMICO) described here:

> **Accelerated Microstructure Imaging via Convex Optimization (AMICO) from diffusion MRI data**  
> *Alessandro Daducci, Erick Canales-Rodriguez, Hui Zhang, Tim Dyrby, Daniel Alexander, Jean-Philippe Thiran*  
> NeuroImage 105, pp. 32-44 (2015)

## Code implementation

This is the first/original version implementation of the AMICO framework and is written in MATLAB. **THIS CODE IS NO LONGER MAINTAINED**.

NB: the official version of AMICO is written in **Python** can be found [`here`](https://github.com/daducci/AMICO).

# Installation

## Download and install external software

- **NODDI MATLAB toolbox**. [Download](http://mig.cs.ucl.ac.uk/index.php?n=Download.NODDI) the software and follow the instructions provided [here](http://mig.cs.ucl.ac.uk/index.php?n=Tutorial.NODDImatlab) to install it.  

- **CAMINO toolkit**. [Download](http://cmic.cs.ucl.ac.uk/camino//index.php?n=Main.Download) the software and follow the instructions provided [here](http://cmic.cs.ucl.ac.uk/camino//index.php?n=Main.Installation) to install it.  

- **SPArse Modeling Software (SPAMS)**. [Download](http://spams-devel.gforge.inria.fr/downloads.html) the software and follow the instructions provided [here](http://spams-devel.gforge.inria.fr/doc/html/doc_spams003.html) to install it.  

## Setup paths/variables in MATLAB

Add the folder containing the source code of AMICO to your `MATLAB PATH`.

Copy the file `AMICO_Setup.txt` and rename it to `AMICO_Setup.m`. Modify its content to set the paths to your specific needs, as follows:

- `AMICO_code_path` : path to the folder containing the *MATLAB source code* of AMICO (this repository). E.g. `/home/user/AMICO/code/matlab`.

- `NODDI_path` : path to the folder containing the *source code* of the NODDI toolbox (in case you want to use NODDI, not needed for ActiveAx). E.g. `/home/user/NODDI_toolbox_v0.9`.

- `CAMINO_path` : path to the `bin` folder containing the *executables* of the Camino toolkit (in case you want to use ActiveAx, not needed for NODDI). E.g. `/home/user/camino/bin`.

- `SPAMS_path` : path to the folder containing the *source code* of the SPAMS Library. E.g. `/home/user/spams`.

- `AMICO_data_path` : path to the folder where you store all your datasets. E.g. `/home/user/AMICO/data`. Then, the software assumes the folder structure is the following:

    ```
    ├── data
        ├── Study_01                 --> all subjects acquired with protocol "Study_01"
            ├── Subject_01
            ├── Subject_02
            ├── ...
        ├── Study_02                 --> all subjects acquired with protocol "Study_02"
            ├── Subject_01
            ├── Subject_02
            ├── ...
        ├── ...
    ```
  This way, the kernels need to be computed only *once per each study*, i.e. same protocol (number of shells, b-values etc), and subsequently adapted to each subject (specific gradient directions) very efficiently.


# Getting started

Tutorials/demos are provided in the folder [`doc/demos/`](doc/demos/) to help you get started with the AMICO framework.
# Tutorials/demos

This folder contains a series of tutorials/demos to show how to use the AMICO framework.# NODDI fitting tutorial

This tutorial shows how to use the AMICO framework to **fit the NODDI model**, using the example dataset distributed with the [NODDI Matlab Toolbox](http://mig.cs.ucl.ac.uk/index.php?n=Tutorial.NODDImatlab).

## Download data for this tutorial

1. Download the original DWI data from [here](http://www.nitrc.org/frs/download.php/5508/NODDI_example_dataset.zip).
2. Create the folder `NoddiTutorial/Tutorial` in your data folder and extract into it the content of the downloaded archive `NODDI_example_dataset.zip`.
3. Copy the scheme file `NODDI_DWI.scheme` distributed with this tutorial into the folder.

## Setup AMICO

Setup the AMICO environment:

```matlab
clearvars, clearvars -global, clc

% Setup AMICO
AMICO_Setup

% Pre-compute auxiliary matrices to speed-up the computations
AMICO_PrecomputeRotationMatrices(); % NB: this needs to be done only once and for all
```

## Load the data

Load the data:

```matlab
% Set the folder containing the data (relative to the data folder).
% This will create a CONFIG structure to keep all the parameters.
AMICO_SetSubject( 'NoddiTutorial', 'Tutorial' );

% Override default file names
CONFIG.dwiFilename    = fullfile( CONFIG.DATA_path, 'NODDI_DWI.hdr' );
CONFIG.maskFilename   = fullfile( CONFIG.DATA_path, 'roi_mask.hdr' );
CONFIG.schemeFilename = fullfile( CONFIG.DATA_path, 'NODDI_DWI.scheme' );

% Load the dataset in memory
AMICO_LoadData
```

The output will look like:

```
-> Loading and setup:
	* Loading DWI...
		- dim    = 128 x 128 x 50 x 81
		- pixdim = 1.875 x 1.875 x 2.500
	* Loading SCHEME...
		- 81 measurements divided in 2 shells (9 b=0)
	* Loading MASK...
		- dim    = 128 x 128 x 50
		- voxels = 5478
   [ DONE ]
```

## Generate the kernels

Generate the kernels corresponding to the different compartments of the NODDI model:

```matlab
% Setup AMICO to use the 'NODDI' model
AMICO_SetModel( 'NODDI' );

% Generate the kernels corresponding to the protocol
AMICO_GenerateKernels( false );

% Resample the kernels to match the specific subject's scheme
AMICO_ResampleKernels();
```

The output will look something like:

```
-> Generating kernels for protocol "NoddiTutorial":
	* Creating high-resolution scheme:
	  [ DONE ] 
	* Simulating "NODDI" kernels:
		- A_001... [2.4 seconds]
		- A_002... [2.4 seconds]	

        ...

	    - A_144... [2.5 seconds]
	    - A_145... [0.0 seconds]
      [ 362.5 seconds ]
   [ DONE 

-> Resampling rotated kernels for subject "Tutorial":
	- A_001...  [0.5 seconds]
	- A_002...  [0.4 seconds]

    ...

	- A_144...  [0.4 seconds]
	- A_145...  [0.0 seconds]
	- saving... [5.2 seconds]
   [ 67.1 seconds ]
```


## Fit the model

Actually **fit** the NODDI model using the AMICO framework:

```matlab
AMICO_Fit()
```

The output will look something like:

```
-> Fitting NODDI model to data:
   [ 0h 0m 15s ]

-> Saving output maps:
   [ AMICO/FIT_*.nii ]
```

![NRMSE for COMMIT](https://github.com/daducci/AMICO/blob/master/matlab/doc/demos/NODDI/RESULTS_Fig1.png)

The results will be saved as NIFTI/ANALYZE files in `NoddiTutorial/Tutorial/AMICO/`.


# ActiveAx fitting tutorial

This tutorial shows how to use the AMICO framework to **fit the ActiveAx model**, using the example dataset distributed with the [ActiveAx tutorial](http://cmic.cs.ucl.ac.uk/camino/index.php?n=Tutorials.ActiveAx) in the Camino toolkit.

## Download data for this tutorial

1. Download the original DWI data from [here](http://dig.drcmr.dk/activeax-dataset/).

2. Create the folder `ActiveAxTutorial/Tutorial` in your data directory and merge the downloaded datasets into one:
```bash
export FSLOUTPUTTYPE=NIFTI
fslmerge -t DWI.nii \\
DRCMR_ActiveAx4CCfit_E2503_Mbrain1_B13183_3B0_ELEC_N90_Scan1_DIG.nii \\
DRCMR_ActiveAx4CCfit_E2503_Mbrain1_B1925_3B0_ELEC_N90_Scan1_DIG.nii \\
DRCMR_ActiveAx4CCfit_E2503_Mbrain1_B1931_3B0_ELEC_N90_Scan1_DIG.nii \\
DRCMR_ActiveAx4CCfit_E2503_Mbrain1_B3091_3B0_ELEC_N90_Scan1_DIG.nii
```

3. Download the scheme file from [here](http://cmic.cs.ucl.ac.uk/camino/uploads/Tutorials/ActiveAxG140_PM.scheme1) and save it into the same folder.
4. Download the binary mask of the corpus callosum from [here](http://hardi.epfl.ch/static/data/AMICO_demos/ActiveAx_Tutorial_MidSagCC.nii) to the same folder.

## Setup AMICO

Setup the AMICO environment:

```matlab
clearvars, clearvars -global, clc

% Setup AMICO
AMICO_Setup

% Pre-compute auxiliary matrices to speed-up the computations
AMICO_PrecomputeRotationMatrices(); % NB: this needs to be done only once and for all
```

## Load the data

Load the data:

```matlab
% Set the folder containing the data (relative to the data folder).
% This will create a CONFIG structure to keep all the parameters.
AMICO_SetSubject( 'ActiveAxTutorial', 'Tutorial' );

% Override default file names
CONFIG.dwiFilename    = fullfile( CONFIG.DATA_path, 'DWI.nii' );
CONFIG.maskFilename   = fullfile( CONFIG.DATA_path, 'ActiveAx_Tutorial_MidSagCC.nii' );
CONFIG.schemeFilename = fullfile( CONFIG.DATA_path, 'ActiveAxG140_PM.scheme1' );

% Load the dataset in memory
AMICO_LoadData
```

The output will look like:

```
-> Loading and setup:
	* Loading DWI...
		- dim    = 128 x 256 x 3 x 372
		- pixdim = 0.400 x 0.400 x 0.500
	* Loading SCHEME...
		- 372 measurements divided in 4 shells (12 b=0)
	* Loading MASK...
		- dim    = 128 x 256 x 3
		- voxels = 338
   [ DONE ]
```

## Generate the kernels

Generate the kernels corresponding to the different compartments of the ActiveAx model:

```matlab
% Setup AMICO to use the 'ActiveAx' model
AMICO_SetModel( 'ActiveAx' );

% Generate the kernels corresponding to the protocol
AMICO_GenerateKernels( false );

% Resample the kernels to match the specific subject's scheme
AMICO_ResampleKernels();
```

The output will look something like:

```
-> Generating kernels for protocol "ActiveAxTutorial":
	* Creating high-resolution scheme:
	  [ DONE ]
	* Simulating "ActiveAx" kernels:
		- A_001... [5.3 seconds]
		- A_002... [5.1 seconds]

        ...

		- A_028... [5.7 seconds]
		- A_029... [0.7 seconds]
	  [ 149.2 seconds ]
   [ DONE ]
   
-> Resampling rotated kernels for subject "Tutorial":
	- A_001...  [0.9 seconds]
	- A_002...  [0.9 seconds]

    ...

	- A_028...  [0.9 seconds]
	- A_029...  [0.0 seconds]
	- saving... [5.2 seconds]
   [ 33.4 seconds ]
```


## Fit the model

Actually **fit** the ActiveAx model using the AMICO framework:

```matlab
AMICO_Fit()
```

The output will look something like:

```
-> Fitting ACTIVEAX model to data:
   [ 0h 0m 0s ]

-> Saving output maps:
   [ AMICO/FIT_*.nii ]
```

![NRMSE for COMMIT](https://github.com/daducci/AMICO/blob/master/matlab/doc/demos/ActiveAx/RESULTS_Fig1.png)

The results will be saved as NIFTI/ANALYZE files in `ActiveAxTutorial/Tutorial/AMICO`.


# MOcov [![Build Status](https://travis-ci.org/MOcov/MOcov.svg?branch=master)](https://travis-ci.org/MOcov/MOcov)

MOcov is a coverage report generator for Matlab and GNU Octave.


### Features

- Runs on both the [Matlab] and [GNU Octave] platforms.
- Can be used directly with continuous integration services, such as [coveralls.io] and [Shippable].
- Integrates with [MOxUnit], a unit test framework for Matlab and GNU Octave.
- Supports the Matlab profiler.
- Writes coverage reports in HTML, JSON and XML formats.
- Distributed under the MIT license, a permissive free software license.


### Installation

- Using the shell (requires a Unix-like operating system such as GNU/Linux or Apple OSX):

    ```bash
    git clone https://github.com/MOcov/MOcov.git
    cd MOcov
    make install
    ```
    This will add the MOcov directory to the Matlab and/or GNU Octave search path. If both Matlab and GNU Octave are available on your machine, it will install MOcov for both.

- Manual installation:

    + Download the zip archive from the [MOcov] website.
    + Start Matlab or GNU Octave.
    + On the Matlab or GNU Octave prompt, `cd` to the `MOcov` root directory, then run:
    
        ```matlab
        cd MOcov            % cd to MOcov subdirectory
        addpath(pwd)        % add the current directory to the Matlab/GNU Octave path
        savepath            % save the path
        ```


### Determining coverage

Coverage can be determined for evaluating a single expression or evaluation of a single function handle; for typical use cases this invokes running a test suite. 

There are two methods to generate coverage while evaluating such an expression or function handle:

1. the 'file' method (default)

    - Coverage information is stored internally by the function `mocov_line_covered`, which keeps this information through the use of persistent variables. Initially the coverage information is reset to being empty.
    - This method considers all files in a directory (and its subdirectories).
    - A temporary directory is created where modified versions of each file is stored.
    - Prior to evaluting the expression or function handle, for each file, MOcov determines which of its lines can be executed. Each line that can be executed is prefixed by a call to `mocov_line_covered`, which cause it to update internal state to record the filename and line number that was executed, and the result stored in the temporary directory.
    - The search path is updated to include the new temporary directory.
    
    After evaluating the expression or function handle, the temporary directory is deleted and the search path restored. Line coverage information is then extracted from the internal state of `mocov_line_covered`.
    
    This method runs on both GNU Octave and Matlab, but is typically slow.

2. the 'profile' method
    - It uses the Matlab profiler. 
    - This method runs on Matlab only (not on GNU Octave), but is generally faster.


### Use cases

Typical use cases for MOcov are:

-   Locally run code with coverage for code in a unit test framework on GNU Octave or Matlab. Use

    ```matlab    
        mocov('-cover','path/with/code',...
                '-expression','run_test_command',...
                '-cover_json_file','coverage.json',...
                '-cover_xml_file','coverage.xml',...
                '-cover_html_dir','coverage_html',
                '-method','file');
    ```

    to generate coverage reports for all files in the `'path/with/code'` directory when `running eval('run_test_command')`. Results are stored in JSON, XML and HTML formats.

-   As a specific example of the use case above, when using the [MOxUnit] unit test platform such tests can be run as

    ```matlab
        success=moxunit_runtests('path/with/tests',...
                                    '-with_coverage',...
                                    '-cover','/path/with/code',...
                                    '-cover_xml_file','coverage.xml',...
                                    '-cover_html_dir','coverage_html');
    ```

    where `'path/with/tests'` contains unit tests. In this case, `moxunit_runtests` will call the `mocov` function to generate coverage reports.

-   On the Matlab platform, results from `profile('info')` can be stored in JSON, XML or HTML formats directly. In the following:

    ```matlab
        % enable profiler
        profile on;

        % run code for which coverage is to be determined
        <your code here>

        % write coverage based on profile('info')
        mocov('-cover','path/with/code',...
                '-profile_info',...
                '-cover_json_file','coverage.json',...
                '-cover_xml_file','coverage.xml',...
                '-cover_html_dir','coverage_html');
    ```

    coverage results are stored in JSON, XML and HTML formats.

-   Use with continuous integration service, such as [Shippable] or [travis-ci] combined with [coveralls.io]. See the   [travis.yml configuration file] in the [MOxUnit] project for an example.


### Use with travis-ci and Shippable
MOcov can be used with the [Travis-ci] and [Shippable] services for continuous integration testing. This is achieved by setting up a `travis.yml` file. Due to recursiveness issues, MOcov cannot use these services to generate coverage reports for itself; for an example in the related [MOxUnit] project, see the [travis.yml configuration file] file.


### Compatibility notes
- Because GNU Octave 3.8 and 4.0 do not support `classdef` syntax, 'old-style' object-oriented syntax is used for the class definitions. 


### Limitations
- The 'file' method uses a very simple parser, which may not work as expected in all cases.
- Currently there is only support to generate coverage reports for files in a single directory (and its subdirectory).


### Contact
Nikolaas N. Oosterhof, nikolaas dot oosterhof at unitn dot it


### Contributions
- Thanks to Scott Lowe and Anderson Bravalheri for their contributions.


### License

(The MIT License)

Copyright (c) 2015-2017 Nikolaas N. Oosterhof

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



[GNU Octave]: http://www.gnu.org/software/octave/
[Matlab]: http://www.mathworks.com/products/matlab/
[MOxUnit]: https://github.com/MOxUnit/MOxUnit
[MOcov]: https://github.com/MOcov/MOcov
[MOxUnit .travis.yml]: https://github.com/MOxUnit/MOxUnit/blob/master/.travis.yml
[Travis-ci]: https://travis-ci.org
[coveralls.io]: https://coveralls.io/
[travis.yml configuration file]: https://github.com/MOxUnit/MOxUnit/blob/master/.travis.yml
[Shippable]: https://shippable.com

# MOxUnit [![Build Status](https://travis-ci.org/nno/MOxUnit.svg?branch=master)](https://travis-ci.org/MOxUnit/MOxUnit) ![Test](https://github.com/MOxUnit/MOxUnit/workflows/CI/badge.svg) [![Coverage Status](https://coveralls.io/repos/github/MOxUnit/MOxUnit/badge.svg?branch=master)](https://coveralls.io/github/MOxUnit/MOxUnit?branch=master) <!-- omit in toc -->

MOxUnit is a lightweight unit test framework for Matlab and GNU Octave.

- [Features](#features)
- [Installation](#installation)
- [Defining MOxUnit tests](#defining-moxunit-tests)
- [Running MOxUnit tests](#running-moxunit-tests)
- [Use with CI](#use-with-ci)
  - [Octave](#octave)
  - [Matlab](#matlab)
- [Compatibility notes](#compatibility-notes)
- [Limitations](#limitations)

## Features

- Runs on both the [Matlab] and [GNU Octave] platforms.
- Uses object-oriented TestCase, TestSuite and TestResult classes, allowing for user-defined extensions.
- Can be used directly with continuous integration services, such as [Travis-ci] and [Shippable].
- Supports JUnit-like XML output for use with Shippable and other test results visualization approaches.
- Supports the generation of code coverage reports using [MOCov]
- Provides compatibility with the (now unsupported) Steve Eddin's [Matlab xUnit test framework], and with recent Matlab test functionality.
- Distributed under the MIT license, a permissive free software license.

## Installation

- Using the shell (requires a Unix-like operating system such as GNU/Linux or Apple OSX):

    ```bash
    git clone https://github.com/MOxUnit/MOxUnit.git
    cd MOxUnit
    make install
    ```
    This will add the MOxUnit directory to the Matlab and/or GNU Octave searchpath. If both Matlab and GNU Octave are available on your machine, it will install MOxUnit for both.

- Manual installation:

    + Download the [[MOxUnit zip archive] from the [MOxUnit] website, and extract it. This should
      result in a directory called ``MOxUnit-master``.
    + Start Matlab or GNU Octave.
    + On the Matlab or GNU Octave prompt, go to the directory that contains the new ``MOxUnit-master`` directory, then run:

        ```matlab
        % change to the MOxUnit subdirectory
        %
        % Note: if MOxUnit was retrieved using 'git', then the name of
        %       top-level directory is 'MOxUnit', not 'MOxUnit-master'
        cd MOxUnit-master/MoxUnit

        % add the current directory to the Matlab/GNU Octave path
        moxunit_set_path()

        % save the path
        savepath
        ```

## Defining MOxUnit tests

To define unit tests, write a function with the following header:

```matlab
function test_suite=test_of_abs
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
```

### Important <!-- omit in toc -->

- It is crucial that the output of the main function is a variable named `test_suite`, and that the output of `localfunctions` is assigned to a variable named `test_functions`.
- As of Matlab 2016b, Matlab scripts (such as `initTestSuite.m`) do not have access to subfunctions in a function if called from that function. Therefore it requires using localfunctions to obtain function handles to local functions. The "try-catch-end" statements are necessary for compatibility with older versions of GNU Octave, which do not provide the `localfunctions` function.
- Alas, the call to `localfunctions` **cannot** be incorporated into `initTestSuite` so the entire code snippet above has to be the header of each test file

Then, define subfunctions whose name start with `test` or end with `test` (case-insensitive). These functions can use the following `assert*` functions:

- `assertTrue(a)`: assert that `a` is true.
- `assertFalse(a)`: assert that `a` is false.
- `assertEqual(a,b)`: assert that `a` and `b` are equal.
- `assertElementsAlmostEqual(a,b)`: assert that the floating point arrays `a` and `b` have the same size, and that corresponding elements are equal within some numeric tolerance.
- `assertVectorsAlmostEqual(a,b)`: assert that floating point vectors `a` and `b` have the same size, and are equal within some numeric tolerance based on their vector norm.
- `assertExceptionThrown(f,id)`: assert that calling `f()` throws an exception with identifier `id`. (To deal with cases where Matlab and GNU Octave throw errors with different identifiers, use `moxunit_util_platform_is_octave`. Or use `id='*'` to match any identifier).

As a special case, `moxunit_throw_test_skipped_exception('reason')` throws an exception that is caught when running the test; `moxunit_run_tests` will report that the test is skipped for reason `reason`.

For example, the following function defines three unit tests that tests some possible inputs from the builtin `abs` function:

```matlab
function test_suite=test_of_abs
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_abs_scalar
    assertTrue(abs(-1)==1)
    assertEqual(abs(-NaN),NaN);
    assertEqual(abs(-Inf),Inf);
    assertEqual(abs(0),0)
    assertElementsAlmostEqual(abs(-1e-13),0)

function test_abs_vector
    assertEqual(abs([-1 1 -3]),[1 1 3]);

function test_abs_exceptions
    % GNU Octave and Matlab use different error identifiers
    if moxunit_util_platform_is_octave()
        assertExceptionThrown(@()abs(struct),'');
    else
        assertExceptionThrown(@()abs(struct),...
                             'MATLAB:UndefinedFunction');
    end
```

Examples of unit tests are in MOxUnit's `tests` directory, which test some of MOxUnit's functions itself.

## Running MOxUnit tests

- `cd` to the directory where the unit tests reside. For MOxUnit itself, the unit tests are in the directory `tests`.
- run the tests using `moxunit_runtests`. For example, running `moxunit_runtests` from MOxUnit's `tests` directory runs tests for MOxUnit itself, and should give the following output:

  ```matlab
  >> moxunit_runtests
  suite: 98 tests
  ............................................................
  ......................................
  --------------------------------------------------

  OK (passed=98)
  ans =
    logical
    1
  ```

- `moxunit_runtests`, by default, gives non-verbose output and runs all tests in the current directory. This can be changed using the following arguments:
  - `-verbose`: show verbose output.
  - `-quiet`: supress all output
  - `directory`: run unit tests in directory `directory`.
  - `file.m`: run unit tests in file `file.m`.
  - `-recursive`: add files from directories recursively.
  - `-logfile logfile.txt`: store the output in file `logfile.txt`.
  - `-junit_xml_file xmlfile`: store JUnit-like XML output in file `xmlfile`.

- To test MOxUnit itself from a terminal, run:

  ```
  make test
  ```

## Use with CI

MOxUnit can be used with the [Travis-ci] service for continuous integration (CI) testing. This is achieved by setting up a [.travis.yml configuration file](.travis.yml). This file is also used by [Shippable]. As a result, the test suite is run automatically on both [Travis-ci] and [Shippable] every time it is pushed to the github repository, or when a pull request is made. If a test fails, or if all tests pass after a test failed before, the developers are notified by email.

### Octave

The easiest test to set up on Travis and/or Shippable is with [GNU Octave]. Make sure your code is Octave compatible. Note that many Matlab projects tend to use functionality not present in Octave (such as particular functions), whereasand writing code that is both Matlab- and Octave-compatible may require some additional efforts.

A simple `.travis.yml` file for a project could look like that:

```yaml
language: generic
os: linux
      
before_install:
  - sudo apt-get install octave

before_script:
  - git clone https://github.com/MOxUnit/MOxUnit.git
  - make -C MOxUnit install

script:        
  - make test
```

In this case `make test` is used to run the tests. To avoid a Makefile and run tests directly through Octave, the script has to call Octave directly to run the tests:

  ```yaml
  # ...
  before_script:
  - git clone https://github.com/MOxUnit/MOxUnit.git

  script:
    - octave --no-gui --eval "addpath('~/git/MOxUnit/MOxUnit');moxunit_set_path;moxunit_runtests('tests')"
  ```

Note that MOxUnit tests **itself** on travis, with [this](https://github.com/MOxUnit/MOxUnit/blob/master/.travis.yml) travis file.

### Matlab

Travis [now supports Matlab](https://docs.travis-ci.com/user/languages/matlab/) directly. You can use MOxUnit with it, but its tricky because:
  1) Travis only supports Matlab 2020a and, presumably, higher (at the time of writing 2020a is the newest version).
  2) Makefile installation does not work with Matlab on travis.
  3) Nor does calling Matlab from the command line in a usual way - with ` matlab -nodesktop -nosplash ...` . Instead it has to be called with the `-batch` flag.
  
  Therefore, `.travis.yml` file looks as follows:
  ```yml
  language: matlab
  matlab: R2020a
  os: linux

  # Just clone MOxUnit, `don't make install` it (!)
  before_script:
    - git clone https://github.com/MOxUnit/MOxUnit.git
      
  script: 
    - matlab -batch 'back=cd("./MOxUnit/MOxUnit/"); moxunit_set_path(); cd(back); moxunit_runtests tests -verbose; exit(double(~ans))'
  ```

  `exit(double(~ans))` ensures that the build fails if MOxUnit tests fail.

## Compatibility notes

- Because GNU Octave 3.8 does not support `classdef` syntax, 'old-style' object-oriented syntax is used for the class definitions. For similar reasons, MOxUnit uses the `lasterror` function, even though its use in Matlab is discouraged.
- Recent versions of Matlab (2016 and later) do not support tests defined just using "initTestSuite", that is without the use of `localfunctions` (see above). To ease the transition, consider using the Python script `tools/fix_mfile_test_init.py`, which can update existing .m files that do not use `localfunctions`.

  For example, the following command was used on a Unix-like shell to preview changes to MOxUnit's tests:

  ```bash
    find tests -iname 'test*.m' | xargs -L1 tools/fix_mfile_test_init.py
  ```

  and adding the `--apply` option applies these changes, meaning that found files are rewritten:

  ```bash
    find tests -iname 'test*.m' | xargs -L1 tools/fix_mfile_test_init.py --apply
  ```
- Recent versions of Matlab define a `matlab.unittest.Test` class for unit tests. An instance `t` can be used with MOxUnit using the `MOxUnitMatlabUnitWrapperTestCase(t)`, which is a `MOxUnitTestCase` instance. Tests that are defined through

  ```matlab
  function tests=foo()
      tests=functiontests(localfunctions)

  function test_funcA(param)

  function test_funcA(param)
  ```

  can be run using MOxUnit as well (and included in an ``MOxUnitTestSuite`` instance using its with ``addFile``) instance, with the exception that currently setup and teardown functions are currently ignored.

## Limitations

Currently MOxUnit does not support:
- Documentation tests require [MOdox].
- Support for setup and teardown functions in `TestCase` classes.
- Subclasses of MOxUnit's classes (`MOxUnitTestCase`, `MOxUnitTestSuite`, `MOxUnitTestReport`) have to be defined using "old-style" object-oriented syntax.
- Subtests

## Acknowledgements <!-- omit in toc -->

- The object-oriented class structure was inspired by the [Python unit test] framework.
- The `assert*` function signatures are aimed to be compatible with Steve Eddin's [Matlab xUnit test framework].

## Contact <!-- omit in toc -->

Nikolaas N. Oosterhof, n dot n dot oosterhof at googlemail dot com.

## Contributions <!-- omit in toc -->

- Thanks to Scott Lowe, Thomas Feher, Joel LeBlanc, Anderson Bravalheri, Sven Baars, 'jdbancal', Marcin Konowalczyk for contributions.

## License <!-- omit in toc -->

(The MIT License)

Copyright (c) 2015-2020 Nikolaas N. Oosterhof

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge,
publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


[GNU Octave]: http://www.gnu.org/software/octave/
[Matlab]: http://www.mathworks.com/products/matlab/
[Matlab xUnit test framework]: http://it.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework
[MOdox]: https://github.com/MOdox/MOdox
[MOxUnit]: https://github.com/MOxUnit/MOxUnit
[MOxUnit zip archive]: https://github.com/MOxUnit/MOxUnit/archive/master.zip
[MOcov]: https://github.com/MOcov/MOcov
[Python unit test]: https://docs.python.org/2.6/library/unittest.html
[Travis-ci]: https://travis-ci.org
[Shippable]: https://app.shippable.com/


 "MPPCA": 4d image denoising and noise map estimation by exploiting  data redundancy in the PCA domain using universal properties of the eigenspectrum of
     random covariance matrices, i.e. Marchenko Pastur distribution
    
      [Signal, Sigma] = MPdenoising(data, mask, kernel, sampling)
           output:
               - Signal: [x, y, z, M] denoised data matrix
               - Sigma: [x, y, z] noise map
           input:
               - data: [x, y, z, M] data matrix
               - mask:   (optional)  region-of-interest [boolean]
               - kernel: (optional)  window size, typically in order of [5 x 5 x 5]
               - sampling: 
                        1. full: sliding window (default for noise map estimation, i.e. [Signal, Sigma] = MPdenoising(...) )
                        2. fast: block processing (default for denoising, i.e. [Signal] = MPdenoising(...))
     
      Authors: Jelle Veraart (jelle.veraart@nyumc.org)
     Copyright (c) 2016 New York Universit and University of Antwerp
           
          Permission is hereby granted, free of charge, to any non-commercial entity
          ('Recipient') obtaining a copy of this software and associated
          documentation files (the 'Software'), to the Software solely for
          non-commercial research, including the rights to use, copy and modify the
          Software, subject to the following conditions: 
           
            1. The above copyright notice and this permission notice shall be
          included by Recipient in all copies or substantial portions of the
          Software. 
           
            2. THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
          EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIESOF
          MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
          NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BELIABLE FOR ANY CLAIM,
          DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
          OTHERWISE, ARISING FROM, OUT OF ORIN CONNECTION WITH THE SOFTWARE OR THE
          USE OR OTHER DEALINGS IN THE SOFTWARE. 
           
            3. In no event shall NYU be liable for direct, indirect, special,
          incidental or consequential damages in connection with the Software.
          Recipient will defend, indemnify and hold NYU harmless from any claims or
          liability resulting from the use of the Software by recipient. 
           
          4. Neither anything contained herein nor the delivery of the Software to
          recipient shall be deemed to grant the Recipient any right or licenses
           under any patents or patent application owned by NYU. 
           
            5. The Software may only be used for non-commercial research and may not
          be used for clinical care. 
           
            6. Any publication by Recipient of research involving the Software shall
          cite the references listed below.
     
     REFERENCES
          Veraart, J.; Fieremans, E. & Novikov, D.S. Diffusion MRI noise mapping
          using random matrix theory Magn. Res. Med., 2016, early view, doi:
          10.1002/mrm.26059
##############################################################################                                                      
      JSONLab: An open-source MATLAB/Octave JSON encoder and decoder             
##############################################################################

* Copyright (C) 2011-2017  Qianqian Fang <q.fang at neu.edu>
* License: BSD or GNU General Public License version 3 (GPL v3), see License*.txt
* Version: 1.5 (Nominus - alpha)


#################
Table of Contents
#################
.. contents::
  :local:
  :depth: 3

============
Introduction
============

JSON (`JavaScript Object Notation <http://www.json.org/>`_) is a highly portable, 
human-readable and " `fat-free <http://en.wikipedia.org/wiki/JSON>`_" text format 
to represent complex and hierarchical data. It is as powerful as `XML <http://en.wikipedia.org/wiki/XML>`_, but less verbose. JSON format is widely used for data-exchange in applications, and is essential for the wild success 
of (programming) `Ajax <http://en.wikipedia.org/wiki/Ajax_>`_ and `Web2.0 <http://en.wikipedia.org/wiki/Web_2.0>`_.

UBJSON (Universal Binary JSON) is a binary JSON format, specifically 
optimized for compact file size and better performance while keeping
the semantics as simple as the text-based JSON format. Using the UBJSON
format allows to wrap complex binary data in a flexible and extensible
structure, making it possible to process complex and large dataset 
without accuracy loss due to text conversions.

We envision that both JSON and its binary version will serve as part of 
the mainstream data-exchange formats for scientific research in the future. 
It will provide the flexibility and generality achieved by other popular 
general-purpose file specifications, such as  `HDF5 <http://www.hdfgroup.org/HDF5/whatishdf5.html>`_, with significantly 
reduced complexity and enhanced performance.

JSONLab is a free and open-source implementation of a JSON/UBJSON encoder 
and a decoder in the native MATLAB language. It can be used to convert a MATLAB 
data structure (array, struct, cell, struct array and cell array) into 
JSON/UBJSON formatted strings, or to decode a JSON/UBJSON file into MATLAB 
data structure. JSONLab supports both MATLAB and `GNU Octave <http://www.gnu.org/software/octave/>`_ (a free MATLAB clone).

================
Installation
================

The installation of JSONLab is no different than any other simple
MATLAB toolbox. You only need to download/unzip the JSONLab package
to a folder, and add the folder's path to MATLAB/Octave's path list
by using the following command:

.. code:: shell

    addpath('/path/to/jsonlab');

If you want to add this path permanently, you need to type "pathtool", 
browse to the jsonlab root folder and add to the list, then click "Save".
Then, run "rehash" in MATLAB, and type "which loadjson", if you see an 
output, that means JSONLab is installed for MATLAB/Octave.


================
Using JSONLab
================

JSONLab provides two functions, loadjson.m -- a MATLAB->JSON decoder, 
and savejson.m -- a MATLAB->JSON encoder, for the text-based JSON, and 
two equivalent functions -- loadubjson and saveubjson for the binary 
JSON. The detailed help info for the four functions can be found below:

----------
loadjson.m
----------

.. code-block:: matlab

      <pre>
        data=loadjson(fname,opt)
           or
        data=loadjson(fname,'param1',value1,'param2',value2,...)

        parse a JSON (JavaScript Object Notation) file or string

        authors:Qianqian Fang (q.fang <at> neu.edu)
        created on 2011/09/09, including previous works from 

                Nedialko Krouchev: http://www.mathworks.com/matlabcentral/fileexchange/25713
                   created on 2009/11/02
                François Glineur: http://www.mathworks.com/matlabcentral/fileexchange/23393
                   created on  2009/03/22
                Joel Feenstra:
                http://www.mathworks.com/matlabcentral/fileexchange/20565
                   created on 2008/07/03

        $Id$

        input:
             fname: input file name, if fname contains "{}" or "[]", fname
                    will be interpreted as a JSON string
             opt: a struct to store parsing options, opt can be replaced by 
                  a list of ('param',value) pairs - the param string is equivalent
                  to a field in opt. opt can have the following 
                  fields (first in [.|.] is the default)

                  opt.SimplifyCell [0|1]: if set to 1, loadjson will call cell2mat
                                for each element of the JSON data, and group 
                                arrays based on the cell2mat rules.
                  opt.FastArrayParser [1|0 or integer]: if set to 1, use a
                                speed-optimized array parser when loading an 
                                array object. The fast array parser may 
                                collapse block arrays into a single large
                                array similar to rules defined in cell2mat; 0 to 
                                use a legacy parser; if set to a larger-than-1
                                value, this option will specify the minimum
                                dimension to enable the fast array parser. For
                                example, if the input is a 3D array, setting
                                FastArrayParser to 1 will return a 3D array;
                                setting to 2 will return a cell array of 2D
                                arrays; setting to 3 will return to a 2D cell
                                array of 1D vectors; setting to 4 will return a
                                3D cell array.
                  opt.ShowProgress [0|1]: if set to 1, loadjson displays a progress bar.

        output:
             dat: a cell array, where {...} blocks are converted into cell arrays,
                  and [...] are converted to arrays

        examples:
             dat=loadjson('{"obj":{"string":"value","array":[1,2,3]}}')
             dat=loadjson(['examples' filesep 'example1.json'])
             dat=loadjson(['examples' filesep 'example1.json'],'SimplifyCell',1)

        license:
            BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details 
       </pre>

----------
savejson.m
----------

.. code-block:: matlab

      <pre>
        json=savejson(rootname,obj,filename)
           or
        json=savejson(rootname,obj,opt)
        json=savejson(rootname,obj,'param1',value1,'param2',value2,...)

        convert a MATLAB object (cell, struct or array) into a JSON (JavaScript
        Object Notation) string

        author: Qianqian Fang (q.fang <at> neu.edu)
        created on 2011/09/09

        $Id$

        input:
             rootname: the name of the root-object, when set to '', the root name
               is ignored, however, when opt.ForceRootName is set to 1 (see below),
               the MATLAB variable name will be used as the root name.
             obj: a MATLAB object (array, cell, cell array, struct, struct array,
             class instance).
             filename: a string for the file name to save the output JSON data.
             opt: a struct for additional options, ignore to use default values.
               opt can have the following fields (first in [.|.] is the default)

               opt.FileName [''|string]: a file name to save the output JSON data
               opt.FloatFormat ['%.10g'|string]: format to show each numeric element
                                of a 1D/2D array;
               opt.ArrayIndent [1|0]: if 1, output explicit data array with
                                precedent indentation; if 0, no indentation
               opt.ArrayToStruct[0|1]: when set to 0, savejson outputs 1D/2D
                                array in JSON array format; if sets to 1, an
                                array will be shown as a struct with fields
                                "_ArrayType_", "_ArraySize_" and "_ArrayData_"; for
                                sparse arrays, the non-zero elements will be
                                saved to _ArrayData_ field in triplet-format i.e.
                                (ix,iy,val) and "_ArrayIsSparse_" will be added
                                with a value of 1; for a complex array, the 
                                _ArrayData_ array will include two columns 
                                (4 for sparse) to record the real and imaginary 
                                parts, and also "_ArrayIsComplex_":1 is added. 
               opt.ParseLogical [0|1]: if this is set to 1, logical array elem
                                will use true/false rather than 1/0.
               opt.SingletArray [0|1]: if this is set to 1, arrays with a single
                                numerical element will be shown without a square
                                bracket, unless it is the root object; if 0, square
                                brackets are forced for any numerical arrays.
               opt.SingletCell  [1|0]: if 1, always enclose a cell with "[]" 
                                even it has only one element; if 0, brackets
                                are ignored when a cell has only 1 element.
               opt.ForceRootName [0|1]: when set to 1 and rootname is empty, savejson
                                will use the name of the passed obj variable as the 
                                root object name; if obj is an expression and 
                                does not have a name, 'root' will be used; if this 
                                is set to 0 and rootname is empty, the root level 
                                will be merged down to the lower level.
               opt.Inf ['"$1_Inf_"'|string]: a customized regular expression pattern
                                to represent +/-Inf. The matched pattern is '([-+]*)Inf'
                                and $1 represents the sign. For those who want to use
                                1e999 to represent Inf, they can set opt.Inf to '$11e999'
               opt.NaN ['"_NaN_"'|string]: a customized regular expression pattern
                                to represent NaN
               opt.JSONP [''|string]: to generate a JSONP output (JSON with padding),
                                for example, if opt.JSONP='foo', the JSON data is
                                wrapped inside a function call as 'foo(...);'
               opt.UnpackHex [1|0]: conver the 0x[hex code] output by loadjson 
                                back to the string form
               opt.SaveBinary [0|1]: 1 - save the JSON file in binary mode; 0 - text mode.
               opt.Compact [0|1]: 1- out compact JSON format (remove all newlines and tabs)

               opt can be replaced by a list of ('param',value) pairs. The param 
               string is equivalent to a field in opt and is case sensitive.
        output:
             json: a string in the JSON format (see http://json.org)

        examples:
             jsonmesh=struct('MeshNode',[0 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1],... 
                      'MeshTetra',[1 2 4 8;1 3 4 8;1 2 6 8;1 5 6 8;1 5 7 8;1 3 7 8],...
                      'MeshTri',[1 2 4;1 2 6;1 3 4;1 3 7;1 5 6;1 5 7;...
                                 2 8 4;2 8 6;3 8 4;3 8 7;5 8 6;5 8 7],...
                      'MeshCreator','FangQ','MeshTitle','T6 Cube',...
                      'SpecialData',[nan, inf, -inf]);
             savejson('jmesh',jsonmesh)
             savejson('',jsonmesh,'ArrayIndent',0,'FloatFormat','\t%.5g')

        license:
            BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details
       </pre>

-------------
loadubjson.m
-------------

.. code-block:: matlab

      <pre>
        data=loadubjson(fname,opt)
           or
        data=loadubjson(fname,'param1',value1,'param2',value2,...)

        parse a JSON (JavaScript Object Notation) file or string

        authors:Qianqian Fang (q.fang <at> neu.edu)
        created on 2013/08/01

        $Id$

        input:
             fname: input file name, if fname contains "{}" or "[]", fname
                    will be interpreted as a UBJSON string
             opt: a struct to store parsing options, opt can be replaced by 
                  a list of ('param',value) pairs - the param string is equivalent
                  to a field in opt. opt can have the following 
                  fields (first in [.|.] is the default)

                  opt.SimplifyCell [0|1]: if set to 1, loadubjson will call cell2mat
                                for each element of the JSON data, and group 
                                arrays based on the cell2mat rules.
                  opt.IntEndian [B|L]: specify the endianness of the integer fields
                                in the UBJSON input data. B - Big-Endian format for 
                                integers (as required in the UBJSON specification); 
                                L - input integer fields are in Little-Endian order.
                  opt.NameIsString [0|1]: for UBJSON Specification Draft 8 or 
                                earlier versions (JSONLab 1.0 final or earlier), 
                                the "name" tag is treated as a string. To load 
                                these UBJSON data, you need to manually set this 
                                flag to 1.

        output:
             dat: a cell array, where {...} blocks are converted into cell arrays,
                  and [...] are converted to arrays

        examples:
             obj=struct('string','value','array',[1 2 3]);
             ubjdata=saveubjson('obj',obj);
             dat=loadubjson(ubjdata)
             dat=loadubjson(['examples' filesep 'example1.ubj'])
             dat=loadubjson(['examples' filesep 'example1.ubj'],'SimplifyCell',1)

        license:
            BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details 
      </pre>

-------------
saveubjson.m
-------------


.. code-block:: matlab


      <pre>
        json=saveubjson(rootname,obj,filename)
           or
        json=saveubjson(rootname,obj,opt)
        json=saveubjson(rootname,obj,'param1',value1,'param2',value2,...)

        convert a MATLAB object (cell, struct or array) into a Universal 
        Binary JSON (UBJSON) binary string

        author: Qianqian Fang (q.fang <at> neu.edu)
        created on 2013/08/17

        $Id$

        input:
             rootname: the name of the root-object, when set to '', the root name
               is ignored, however, when opt.ForceRootName is set to 1 (see below),
               the MATLAB variable name will be used as the root name.
             obj: a MATLAB object (array, cell, cell array, struct, struct array,
             class instance)
             filename: a string for the file name to save the output UBJSON data
             opt: a struct for additional options, ignore to use default values.
               opt can have the following fields (first in [.|.] is the default)

               opt.FileName [''|string]: a file name to save the output JSON data
               opt.ArrayToStruct[0|1]: when set to 0, saveubjson outputs 1D/2D
                                array in JSON array format; if sets to 1, an
                                array will be shown as a struct with fields
                                "_ArrayType_", "_ArraySize_" and "_ArrayData_"; for
                                sparse arrays, the non-zero elements will be
                                saved to _ArrayData_ field in triplet-format i.e.
                                (ix,iy,val) and "_ArrayIsSparse_" will be added
                                with a value of 1; for a complex array, the 
                                _ArrayData_ array will include two columns 
                                (4 for sparse) to record the real and imaginary 
                                parts, and also "_ArrayIsComplex_":1 is added. 
               opt.ParseLogical [1|0]: if this is set to 1, logical array elem
                                will use true/false rather than 1/0.
               opt.SingletArray [0|1]: if this is set to 1, arrays with a single
                                numerical element will be shown without a square
                                bracket, unless it is the root object; if 0, square
                                brackets are forced for any numerical arrays.
               opt.SingletCell  [1|0]: if 1, always enclose a cell with "[]" 
                                even it has only one element; if 0, brackets
                                are ignored when a cell has only 1 element.
               opt.ForceRootName [0|1]: when set to 1 and rootname is empty, saveubjson
                                will use the name of the passed obj variable as the 
                                root object name; if obj is an expression and 
                                does not have a name, 'root' will be used; if this 
                                is set to 0 and rootname is empty, the root level 
                                will be merged down to the lower level.
               opt.JSONP [''|string]: to generate a JSONP output (JSON with padding),
                                for example, if opt.JSON='foo', the JSON data is
                                wrapped inside a function call as 'foo(...);'
               opt.UnpackHex [1|0]: conver the 0x[hex code] output by loadjson 
                                back to the string form

               opt can be replaced by a list of ('param',value) pairs. The param 
               string is equivalent to a field in opt and is case sensitive.
        output:
             json: a binary string in the UBJSON format (see http://ubjson.org)

        examples:
             jsonmesh=struct('MeshNode',[0 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1],... 
                      'MeshTetra',[1 2 4 8;1 3 4 8;1 2 6 8;1 5 6 8;1 5 7 8;1 3 7 8],...
                      'MeshTri',[1 2 4;1 2 6;1 3 4;1 3 7;1 5 6;1 5 7;...
                                 2 8 4;2 8 6;3 8 4;3 8 7;5 8 6;5 8 7],...
                      'MeshCreator','FangQ','MeshTitle','T6 Cube',...
                      'SpecialData',[nan, inf, -inf]);
             saveubjson('jsonmesh',jsonmesh)
             saveubjson('jsonmesh',jsonmesh,'meshdata.ubj')

        license:
            BSD or GPL version 3, see LICENSE_{BSD,GPLv3}.txt files for details
      </pre>

---------
examples
---------

Under the ``"examples"`` folder, you can find several scripts to demonstrate the
basic utilities of JSONLab. Running the ``"demo_jsonlab_basic.m"`` script, you 
will see the conversions from MATLAB data structure to JSON text and backward.
In ``"jsonlab_selftest.m"``, we load complex JSON files downloaded from the Internet
and validate the ``loadjson/savejson`` functions for regression testing purposes.
Similarly, a ``"demo_ubjson_basic.m"`` script is provided to test the saveubjson
and loadubjson functions for various matlab data structures.

Please run these examples and understand how JSONLab works before you use
it to process your data.

=======================
Known Issues and TODOs
=======================

JSONLab has several known limitations. We are striving to make it more general
and robust. Hopefully in a few future releases, the limitations become less.

Here are the known issues:

  * 3D or higher dimensional cell/struct-arrays will be converted to 2D arrays
  
  * When processing names containing multi-byte characters, Octave and MATLAB can give different field-names; you can use feature('DefaultCharacterSet','latin1') in MATLAB to get consistent results
  
  * savejson can not handle class and dataset.
  
  * saveubjson converts a logical array into a uint8 ([U]) array
  
  * an unofficial N-D array count syntax is implemented in saveubjson. We are actively communicating with the UBJSON spec maintainer to investigate the possibility of making it upstream 
  
  * loadubjson can not parse all UBJSON Specification (Draft 9) compliant files, however, it can parse all UBJSON files produced by saveubjson.

==========================
Contribution and feedback
==========================

JSONLab is an open-source project. This means you can not only use it and modify
it as you wish, but also you can contribute your changes back to JSONLab so
that everyone else can enjoy the improvement. For anyone who want to contribute,
please download JSONLab source code from its source code repositories by using the
following command:


.. code:: shell

      git clone https://github.com/fangq/jsonlab.git jsonlab

or browsing the github site at

.. code:: shell

      https://github.com/fangq/jsonlab
 

alternatively, if you prefer svn, you can checkout the latest code by using

.. code:: shell

       svn checkout svn://svn.code.sf.net/p/iso2mesh/code/trunk/jsonlab jsonlab

You can make changes to the files as needed. Once you are satisfied with your
changes, and ready to share it with others, please cd the root directory of 
JSONLab, and type

.. code:: shell

      git diff --no-prefix > yourname_featurename.patch
 

or

.. code:: shell

      svn diff > yourname_featurename.patch

You then email the .patch file to JSONLab's maintainer, Qianqian Fang, at
the email address shown in the beginning of this file. Qianqian will review 
the changes and commit it to the subversion if they are satisfactory.

We appreciate any suggestions and feedbacks from you. Please use the following
mailing list to report any questions you may have regarding JSONLab:

`forum/jsonlab-users <https://groups.google.com/forum/?hl=en#>`_

(Subscription to the mailing list is needed in order to post messages).
