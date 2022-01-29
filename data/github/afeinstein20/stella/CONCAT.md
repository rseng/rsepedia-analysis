---
title: '\texttt{stella}: Convolutional Neural Networks for Flare Identification in \textit{TESS}'
tags:
  - Python
  - astronomy
  - PMS stars
  - stellar activity
  - stellar rotation
authors:
  - name: Adina D. Feinstein
    orcid: 0000-0002-9464-8101 
    affiliation: "1, 2"
  - name: Benjamin T. Montet
    orcid: 0000-0001-7516-8308
    affiliation: 3
  - name: Megan Ansdell
    affiliation: 4
 
affiliations:
 - name: Department of Astronomy and Astrophysics, University of Chicago, 5640 S. Ellis Ave, Chicago, IL 60637, USA
   index: 1
 - name: NSF Graduate Research Fellow
   index: 2
 - name: School of Physics, University of New South Wales, Sydney, NSW 2052, Australia
   index: 3
 - name: Flatiron Institute, Simons Foundation, 162 Fifth Ave, New York, NY 10010, USA
   index: 4

date: 4 May 2020
bibliography: paper.bib

aas-doi: 10.3847/1538-3881/abac0a
aas-journal: Astrophysical Journal
---

# Summary

Nearby young moving groups are kinematically bound systems of stars that are believed to have formed at the same time.
With all member stars having the same age, they provide snapshots of stellar and planetary evolution. 
In particular, young ($<$ 800 Myr) stars have increased levels of activity, seen in both fast rotation periods, large spot modulation, and increased flare rates [@zuckerman:2004; @ilin:2019].
Flare rates and energies can yield consequences for the early stages of planet formation, particularly with regards to their atmospheres.
Models have demonstrated that the introduction of superflares ($> 5\%$ flux increase) are able to irreparably alter the chemistry of an atmosphere [@venot:2016] and expedite atmospheric photoevaporation [@lammer:2007]. 
Thus, understanding flare rates and energies at young ages provides crucial keys for understanding the exoplanet population we see today.

Previous methods of flare detection with both *Kepler* [@borucki:2010] and *Transiting Exooplanet Survey Satellite* (*TESS*; @ricker:2014) data have relied on detrending a light curve and using outlier detection heuristics for identifying flare events [@Davenport:2016; @allesfitter]. 
More complex methods, such as a RANdom SAmple Consensus (RANSAC) algorithm has been tested as well [@vida18]. RANSAC algorithms identfy and subtract inliers (the underlying light curve) before searching for outliers above a given detection threshold.
Low-amplitude flares can easily be removed with aggressive detrending techniques (e.g. using a small window-length to remove spot modulation). 
Additionally, low energy flares likely fall below the outlier threshold, biasing the overall flare sample towards higher energy flares.
As flares exhibit similar temporal evolution (a sharp rise followed by an exponential decay, with the exception of complex flare groups), machine learning algorithms may prove suitable for identifying such features without light curve detrending.

`stella` is an open-source Python package for identifying flares in the *TESS* two-minute data with convolutional neural networks (CNNs).
Users have the option to use the models created in @Feinstein:2020 or build their own cutomized networks.
The training, validation, and test sets for our CNNs use the flare catalog presented in @guenther:2020. These light curves are publicly available through the Mikulski Archive for Space Telescopes and can be downloaded through \texttt{stella} as a wrapper around the \texttt{lightkurve} package [@lightkurve]; they are not, by default, included in the package.
It takes approximately twenty minutes to train a \texttt{stella} model from scratch and $<1$ minute to predict flares on a single sector light curve.
The package also allows users to measure rotation periods and fit flares to extract underlying flare parameters. Further documentation and tutorials can be found at \url{adina.feinste.in/stella}.

# Acknowledgements

We acknowledge contributions from Travis expert Rodrigo Luger.
This material is based upon work supported by the National Science Foundation Graduate Research Fellowship Program under Grant No. (DGE-1746045).
This work was funded in part through the NASA *TESS* Guest Investigator Program, as a part of Program G011237 (PI Montet).
<p align="center">
  <img width = "500" src="./figures/stella_logo.png"/>
</p>

<p align="center">
  <a href="https://github.com/afeinstein20/stella/actions?query=workflow%3Astella-tests"><img src="https://github.com/afeinstein20/stella/workflows/stella-tests/badge.svg"?color=D35968/></a>
  <a href="https://ui.adsabs.harvard.edu/abs/2020AJ....160..219F/abstract"><img src="https://img.shields.io/badge/read-the_paper-3C1370.svg?style=flat"/></a>
  <a href="https://afeinstein20.github.io/stella/"><img src="https://img.shields.io/badge/read-the_docs-3C1370.svg?style=flat"/></a>
  <a href="https://pypi.org/project/stella"><img alt="PyPI - Downloads" src="https://img.shields.io/pypi/dm/stella?color=D35968"></a>
  <a href="https://doi.org/10.21105/joss.02347">   <img src="https://joss.theoj.org/papers/10.21105/joss.02347/status.svg?color=D35968"></a>
</p>


</p>
stella is a Python package to create and train a neural network to identify stellar flares.
Within stella, users can simulate flares as a training set, run a neural network, and feed
in their own data to the neural network model. stella returns a probability at each data point
that that data point is part of a flare or not. stella can also characterize the flares identified.
</p>


To install stella with pip:

	pip install stella

Alternatively you can install the current development version of stella:

        git clone https://github.com/afeinstein20/stella
        cd stella
        python setup.py install

<p>
If your work uses the stella software, please cite <a href="https://ui.adsabs.harvard.edu/abs/2020JOSS....5.2347F/abstract">Feinstein, Montet, & Ansdell (2020)</a>.
</p>
<p>
If your work discusses the flare rate of young stars in the TESS Southern Hemisphere or the details of the CNNs, please cite <a href="https://ui.adsabs.harvard.edu/abs/2020arXiv200507710F/abstract">Feinstein et al. (AJ, 2020)</a>.
</p>

<p>
<b><u>Bug Reports, Questions, & Contributions</u></b>
</p>
<p>
stella is an open source project under the MIT license. 
The source code is available on GitHub. In case of any questions or problems, please contact us via the Git Issues. 
Pull requests are also welcome through the GitHub page.
</p># Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at afeinstein@uchicago.edu. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
