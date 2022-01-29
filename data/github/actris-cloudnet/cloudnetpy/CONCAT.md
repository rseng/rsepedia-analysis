# CloudnetPy

![](https://github.com/actris-cloudnet/cloudnetpy/workflows/CloudnetPy%20CI/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/cloudnetpy/badge/?version=latest)](https://cloudnetpy.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/cloudnetpy.svg)](https://badge.fury.io/py/cloudnetpy)
[![DOI](https://zenodo.org/badge/233602651.svg)](https://zenodo.org/badge/latestdoi/233602651)
[![status](https://joss.theoj.org/papers/959971f196f617dddc0e7d8333ff22b7/status.svg)](https://joss.theoj.org/papers/959971f196f617dddc0e7d8333ff22b7)

CloudnetPy is a Python software for producing vertical profiles of cloud properties from ground-based 
remote sensing measurements. The Cloudnet processing combines cloud radar, optical lidar, microwave 
radiometer and model data. Measurements and model data are brought into common grid and 
classified as ice, liquid, aerosol, insects, and so on. 
Then, geophysical products such as ice water content can be 
retrieved in the further processing steps. See [Illingworth et. al. (2007)](https://doi.org/10.1175/BAMS-88-6-883) for more details about the concept.

CloudnetPy is a rewritten version of the original Cloudnet Matlab code. CloudnetPy features several revised methods, extensive documentation, and more.

* CloudnetPy documentation: https://cloudnetpy.readthedocs.io/en/latest/
* Cloudnet data portal: https://cloudnet.fmi.fi

<img src="docs/source/_static/20190423_mace-head_classification.png">

## Installation

### From PyPI
```
$ python3 -m pip install cloudnetpy
```

### From the source
```
$ git clone https://github.com/actris-cloudnet/cloudnetpy
$ cd cloudnetpy/
$ python3 -m venv venv
$ source venv/bin/activate
(venv) $ python3 -m pip install .
```
## Citing
If you wish to acknowledge CloudnetPy in your publication, please cite:
>Tukiainen et al., (2020). CloudnetPy: A Python package for processing cloud remote sensing data. Journal of Open Source Software, 5(53), 2123, https://doi.org/10.21105/joss.02123

## Contributing

We encourage you to contribute to CloudnetPy! Please check out the [contribution guidelines](CONTRIBUTING.md) about how to proceed.

## License
MIT
# Acknowledgements

CloudnetPy contains code and ideas from several people and would
not exist without their prior contribution.

Below is a partial list. Please add your name if it has been left off.

- Hannes Griesche
- Robin Hogan
- Anthony Illingworth
- Anniina Korpinen
- Ewan O'Connor
- Willi Schimmel
- Simo Tukiainen
- Minttu Tuononen
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
reported by contacting the project team at cloudnet@fmi.fi. All
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
# How to Contribute to CloudnetPy

First off, thanks for taking the time to contribute! :+1:

Reporting bugs
--------------

* Create a new [issue](https://github.com/actris-cloudnet/cloudnetpy/issues)
* Describe the problem and steps to reproduce it
* Suggest a fix if you can but this is not mandatory

Proposing code improvements
---------------------------

If you have a good idea how to improve the code without changing the actual outcome
of the processing, you can directly create a [pull request](https://github.com/actris-cloudnet/cloudnetpy/pulls) (PR). 
Examples of these kind of changes are, for example, better 

* unit tests
* implementations of the helper functions
* docstrings

Try to keep your pull requests small. This will increase the chance to be accepted.

Proposing new features
----------------------

Changes that would alter the actual outcome of the processing need to be carefully
reviewed and tested before accepting. Examples include, for example, modifications in the

* methodology (e.g. replacing the wet bulb method with a newer version)
* threshold values used in different functions / methods
* constant values

Also for these sort of proposals, you can open a new [issue](https://github.com/actris-cloudnet/cloudnetpy/issues)
where the idea can be discussed before (possible) implementation.
---
title: 'CloudnetPy: A Python package for processing cloud remote sensing data'
tags:
  - Python
  - cloud radar
  - lidar
  - microwave radiometer
  - remote sensing
authors:
  - name: Simo Tukiainen
    orcid: 0000-0002-0651-4622
    affiliation: 1
  - name: Ewan O'Connor
    affiliation: 1
  - name: Anniina Korpinen
    affiliation: 1
affiliations:
 - name: Finnish Meteorological Institute, Helsinki, Finland
   index: 1
date: 13 February 2020
bibliography: paper.bib
---

# Summary

Active ground-based remote sensing instruments such as cloud radars and lidars 
provide vertical profiles of clouds and aerosols with high vertical and 
temporal resolution. Cloud radars typically operate in the sub-millimeter 
wavelength region, around 35 or 94 GHz, 
and are sensitive to clouds, particularly ice clouds, rain and insects. Lidars operating 
at visible and near-infrared wavelengths 
on the other hand, are more sensitive to liquid clouds and aerosols. 
Combining these two complementary data sources with temperature and humidity profiles 
from a numerical weather prediction model or radiosonde makes it possible to accurately classify 
the various scattering hydrometeors in the atmosphere, diagnosing them as: rain drops, 
ice particles, melting ice particles, liquid droplets, supercooled liquid droplets, 
drizzle drops, insects and aerosol particles. 
Furthermore, adding a passive microwave radiometer, an instrument measuring 
liquid water path, attenuation corrections and quantitative retrievals of geophysical 
products such as ice water content, liquid water content 
and drizzle properties become feasible [@OConnorEtAl05; @HoganEtAl06].

Methodology and prototype software to combine these different data sources, 
and to retrieve target classification and other products, were developed within 
the EU-funded Cloudnet project [@IllingworthEtAl07]. Since Cloudnet started in 
2002, the network has expanded from 3 stations to a coordinated
and continuously operated network of around 15 stations across Europe. 
The network routinely collects, processes and distributes Cloudnet data (http://cloudnet.fmi.fi). 
While the current methodology has been validated, it is important to develop the Cloudnet software 
so that it can efficiently handle large amounts of data and reliably perform 
continuous data processing. In the forthcoming years, Cloudnet will be one of 
the key components in ACTRIS (Aerosol, Clouds and Trace Gases Research 
Infrastructure) [@ACTRIS_handbook], where the Cloudnet framework 
will process gigabytes of cloud remote sensing data per day 
in near real time. The ACTRIS RI is now in its implementation phase and 
aims to be fully operational in 2025. 

CloudnetPy is a Python implementation of the Cloudnet processing scheme. 
CloudnetPy covers the full Cloudnet processing chain starting from the raw 
measurements and providing similar functionality to the original, 
proprietary Cloudnet software written in Matlab and C. The output from CloudnetPy
is no longer identical to the original scheme because several methods have been 
revised and improved during the refactoring process. For example, as most modern cloud 
radars are polarimetric, CloudnetPy uses the linear depolarization ratio 
to improve the detection of the melting layer and insects. Liquid layer detection is 
now based on the lidar attenuated backscatter profile shape instead of relying only on
threshold values [@TuononenEtAl19]. Detailed 
verification of the updated methods is a subject of future studies. 
The CloudnetPy API is designed to serve the operational cloud remote sensing data 
processing in ACTRIS, but it will be straightforward for site operators and the 
scientific community with access to the raw data to run the software, improve 
existing methods and develop new products.


# Acknowledgements

The research was supported by the European Union Framework Programme for Research and 
Innovation, Horizon 2020 (ACTRIS-2,  grant  no.  654109). The authors would like to 
thank the Academy of Finland for supporting ACTRIS activities in Finland
and Lauri Kangassalo for providing comments on the manuscript.

# References
