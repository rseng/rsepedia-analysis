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
File format description
=======================

All Cloudnet files use ``NETCDF4_CLASSIC`` data model, i.e., ``HDF5`` file format.


Level 1b files
--------------

Lidar file
..........

**Dimensions:**

+------------------------------+
|**Name**                      |
+------------------------------+
|time                          |
+------------------------------+
|range                         |
+------------------------------+


**Variables (all lidars):**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|time                          |float32        |time                     |Time UTC                                                              |hours since 2021-06-21 00:00:00         |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|range                         |float32        |range                    |Range from instrument                                                 |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|altitude                      |int32          |                         |Altitude of site                                                      |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|wavelength                    |int32          |                         |Laser wavelength                                                      |nm                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|calibration_factor            |float32        |                         |Backscatter calibration factor                                        |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|beta_smooth                   |float32        |time, range              |Smoothed attenuated backscatter coefficient                           |sr-1 m-1                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|height                        |float32        |range                    |Height above mean sea level                                           |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|beta                          |float32        |time, range              |Attenuated backscatter coefficient                                    |sr-1 m-1                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|tilt_angle                    |float32        |                         |Tilt angle from vertical                                              |degrees                                 |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|beta_raw                      |float32        |time, range              |Raw attenuated backscatter coefficient                                |sr-1 m-1                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

**Variables (CHM15K specific):**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

**Variables (CL51 specific):**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|laser_energy                  |float32        |time                     |Laser pulse energy                                                    |%                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|detection_status              |float32        |time                     |Detection status                                                      |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|range_resolution              |float32        |                         |Range resolution                                                      |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|background_light              |float32        |time                     |Background light                                                      |mV                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|message_number                |float32        |                         |Message number                                                        |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|warning_flags                 |float32        |time                     |Warning flags                                                         |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|scale                         |float32        |                         |Scale                                                                 |%                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|unit_id                       |float32        |                         |Ceilometer unit number                                                |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|window_transmission           |float32        |                         |Window transmission estimate                                          |%                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|message_subclass              |float32        |                         |Message subclass number                                               |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|backscatter_sum               |float32        |time                     |Sum of detected and normalized backscatter                            |sr-1                                    |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|software_level                |float32        |                         |Software level ID                                                     |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|laser_temperature             |float32        |time                     |Laser temperature                                                     |C                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|number_of_gates               |float32        |                         |Number of range gates in profile                                      |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

Model file
..........

**Dimensions:**

+------------------------------+
|**Name**                      |
+------------------------------+
|time                          |
+------------------------------+
|level                         |
+------------------------------+
|flux_level                    |
+------------------------------+
|frequency                     |
+------------------------------+
|soil_level                    |
+------------------------------+


**Variables (all models):**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|specific_liquid_atten         |float32        |frequency, time, level   |Specific one-way attenuation due to liquid water, per unit liquid wat |(dB km-1)/(g m-3)                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_pressure                  |float32        |time                     |Surface pressure                                                      |Pa                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|qi                            |float32        |time, level              |Gridbox-mean ice water mixing ratio                                   |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|q                             |float32        |time, level              |Specific humidity                                                     |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_wind_u_10m                |float32        |time                     |Zonal wind at 10m                                                     |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_wind_v_10m                |float32        |time                     |Meridional wind at 10m                                                |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|ql                            |float32        |time, level              |Gridbox-mean liquid water mixing ratio                                |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_net_lw                    |float32        |time                     |Surface net downward longwave flux                                    |W m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|K2                            |float32        |frequency, time, level   |Dielectric parameter (K^2) of liquid water                            |dB km-1                                 |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|time                          |float32        |time                     |Hours UTC                                                             |hours since 2021-06-21 00:00:00 +00:00  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|uwind                         |float32        |time, level              |Zonal wind                                                            |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|specific_gas_atten            |float32        |frequency, time, level   |Specific one-way attenuation due to atmospheric gases                 |dB km-1                                 |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|flx_height                    |float32        |time, flux_level         |Height above ground                                                   |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|gas_atten                     |float32        |frequency, time, level   |Two-way attenuation from the ground due to atmospheric gases          |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_down_sens_heat_flx        |float32        |time                     |Sensible heat flux                                                    |W m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|horizontal_resolution         |float32        |                         |Horizontal resolution of model                                        |km                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|rh                            |float32        |time, level              |Relative humidity                                                     |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|specific_saturated_gas_atten  |float32        |frequency, time, level   |Specific one-way attenuation due to atmospheric gases for saturated a |dB km-1                                 |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|wwind                         |float32        |time, level              |Vertical wind                                                         |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_ls_rain                   |float32        |time                     |Large-scale rainfall amount                                           |kg m-2                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|specific_dry_gas_atten        |float32        |frequency, time, level   |Specific one-way attenuation due to atmospheric gases for dry air (no |dB km-1                                 |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|frequency                     |float32        |frequency                |Microwave frequency                                                   |GHz                                     |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_ls_snow                   |float32        |time                     |Large-scale snowfall amount                                           |kg m-2                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|cloud_fraction                |float32        |time, level              |Cloud fraction                                                        |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|pressure                      |float32        |time, level              |Pressure                                                              |Pa                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_net_sw                    |float32        |time                     |Surface net downward shortwave flux                                   |W m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|vwind                         |float32        |time, level              |Meridional wind                                                       |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|temperature                   |float32        |time, level              |Temperature                                                           |K                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|latitude                      |float32        |                         |Latitude of model gridpoint                                           |degrees_N                               |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_down_lat_heat_flx         |float32        |time                     |Latent heat flux                                                      |W m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|longitude                     |float32        |                         |Longitude of model gridpoint                                          |degrees_E                               |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|height                        |float32        |time, level              |Height above ground                                                   |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|forecast_time                 |float32        |time                     |Time since initialization of forecast                                 |hours                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

**Variables (HARMONIE-FMI-6-11 specific):**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_down_sw_direct            |float32        |time                     |Direct downwelling shortwave flux                                     |W m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|qs                            |float32        |time, level              |Gridbox-mean snow mixing ratio                                        |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|qg                            |float32        |time, level              |Gridbox-mean graupel mixing ratio                                     |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|conv_cloud_fraction           |float32        |time, level              |Convective cloud fraction                                             |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|ls_cloud_fraction             |float32        |time, level              |Large scale cloud fraction                                            |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_total_cloud_fraction      |float32        |time                     |Surface total cloud fraction                                          |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_turb_mom_u                |float32        |time                     |Surface zonal turbulent momentum flux                                 |kg m-2 s-1                              |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|soil_depth                    |float32        |time, soil_level         |Depth below ground                                                    |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|qr                            |float32        |time, level              |Gridbox-mean rain mixing ratio                                        |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_temp                      |float32        |time                     |Surface temperature                                                   |K                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|tke                           |float32        |time, level              |Turbulent kinetic energy                                              |J m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_ls_graupel                |float32        |time                     |Large-scale graupel amount                                            |kg m-2                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_down_lw                   |float32        |time                     |Surface downwelling longwave flux                                     |W m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|omega                         |float32        |time, level              |Vertical wind in pressure coordinates                                 |Pa s-1                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_turb_mom_v                |float32        |time                     |Surface meridional turbulent momentum flux                            |kg m-2 s-1                              |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_down_sw                   |float32        |time                     |Surface downwelling shortwave flux                                    |W m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_pressure_amsl             |float32        |time                     |Surface pressure at mean sea level                                    |Pa                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_down_sw_direct_normal     |float32        |time                     |Direct normal downwelling shortwave flux                              |W m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_rh_2m                     |float32        |time                     |Relative humidity at 2m                                               |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

**Variables (ICON-IGLO-12-23 specific):**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|ql_diag                       |float32        |time, level              |Total specific liquid water (diagnostic)                              |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_down_sw_direct            |float32        |time                     |Direct downwelling shortwave flux                                     |W m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|conv_cloud_fraction           |float32        |time, level              |Convective cloud fraction                                             |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|ls_cloud_fraction             |float32        |time, level              |Large scale cloud fraction                                            |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|toa_net_sw                    |float32        |time                     |Top of atmosphere net downward shortwave flux                         |W m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_cloud_fraction            |float32        |time                     |Surface total cloud fraction                                          |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|soil_depth                    |float32        |time, soil_level         |Depth below ground                                                    |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|qr                            |float32        |time, level              |Gridbox-mean rain mixing ratio                                        |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_turb_mom_v                |float32        |time                     |Surface meridional turbulent momentum flux                            |kg m-2 s-1                              |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|altitude                      |float32        |                         |Height of station above mean sea level                                |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|q_diag                        |float32        |time, level              |Total specific humidity (diagnostic)                                  |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|turb_heat_coeff               |float32        |time, flux_level         |Turbulent diffusion coefficients for heat                             |m2 s-1                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|soil_temperature              |float32        |time, soil_level         |Soil temperature                                                      |K                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_wind_gust_10m             |float32        |time                     |Wind gust at 10m                                                      |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|qs                            |float32        |time, level              |Gridbox-mean snow mixing ratio                                        |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_turb_mom_u                |float32        |time                     |Surface zonal turbulent momentum flux                                 |kg m-2 s-1                              |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_conv_snow                 |float32        |time                     |Convective snowfall amount                                            |kg m-2                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|turb_mom_coeff                |float32        |time, flux_level         |Turbulent diffusion coefficients for momentum                         |m2 s-1                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_temp                      |float32        |time                     |Surface temperature                                                   |K                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_q_2m                      |float32        |time                     |Specific humidity at 2m                                               |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_up_sw_diffuse             |float32        |time                     |Diffuse upwelling shortwave flux                                      |W m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_roughness_length          |float32        |time                     |Surface roughness length                                              |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_dewpoint_temp_2m          |float32        |time                     |Dew point temperature at 2m                                           |K                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_conv_rain                 |float32        |time                     |Convective rainfall amount                                            |kg m-2                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_temp_2m                   |float32        |time                     |Temperature at 2m                                                     |K                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|toa_net_lw                    |float32        |time                     |Top of atmosphere net downward longwave flux                          |W m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|qi_diag                       |float32        |time, level              |Total specific ice water (diagnostic)                                 |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_albedo                    |float32        |time                     |Surface albedo                                                        |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_land_cover                |float32        |time                     |Land cover                                                            |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_height_amsl               |float32        |time                     |Surface height above mean sea level                                   |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_down_sw_diffuse           |float32        |time                     |Diffuse downwelling shortwave flux                                    |W m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

**Variables (ECMWF specific):**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_cloud_fraction            |float32        |time                     |Surface total cloud fraction                                          |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_conv_snow                 |float32        |time                     |Convective snowfall amount                                            |kg m-2                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_bl_height                 |float32        |time                     |Boundary layer height                                                 |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_geopotential              |float32        |time                     |Geopotential                                                          |m2 s-2                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|omega                         |float32        |time, level              |Vertical wind in pressure coordinates                                 |Pa s-1                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_conv_rain                 |float32        |time                     |Convective rainfall amount                                            |kg m-2                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_temp_2m                   |float32        |time                     |Temperature at 2m                                                     |K                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_down_sw                   |float32        |time                     |Surface downwelling shortwave flux                                    |W m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sfc_ls_precip_fraction        |float32        |time                     |Large-scale precipitation fraction                                    |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

Mwr file
........

**Dimensions:**

+------------------------------+
|**Name**                      |
+------------------------------+
|time                          |
+------------------------------+


**Variables (all mwrs):**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|LWP                           |float32        |time                     |Liquid water path                                                     |g m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|time                          |int32          |time                     |Time UTC                                                              |seconds since 2001-01-01 00:00:00       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

Radar file
..........

**Dimensions:**

+------------------------------+
|**Name**                      |
+------------------------------+
|time                          |
+------------------------------+
|range                         |
+------------------------------+
|chirp_sequence                |
+------------------------------+


**Variables (all radars):**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|v                             |float32        |time, range              |Doppler velocity                                                      |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|time                          |float32        |time                     |Time UTC                                                              |hours since 2021-06-21 00:00:00         |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|latitude                      |float32        |                         |Latitude of site                                                      |degrees_north                           |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|altitude                      |int32          |                         |Altitude of site                                                      |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|Ze                            |float32        |time, range              |Radar reflectivity factor.                                            |dBZ                                     |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|longitude                     |float32        |                         |Longitude of site                                                     |degrees_east                            |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|height                        |float32        |range                    |Height above mean sea level                                           |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|nyquist_velocity              |float32        |chirp_sequence           |Nyquist velocity                                                      |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|radar_frequency               |float32        |                         |Radar transmit frequency                                              |GHz                                     |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|range                         |float32        |range                    |Range from instrument                                                 |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

**Variables (RPG-FMCW-94 specific):**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|voltage                       |float32        |time                     |Voltage                                                               |V                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|time_ms                       |int32          |time                     |Time ms                                                               |ms                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|pc_temperature                |float32        |time                     |PC temperature                                                        |K                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|azimuth                       |float32        |time                     |Azimuth angle                                                         |degrees                                 |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|status_flag                   |float32        |time                     |Status flag for heater and blower                                     |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|brightness_temperature        |float32        |time                     |Brightness temperature                                                |K                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|if_power                      |float32        |time                     |IF power at ACD                                                       |uW                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|quality_flag                  |int32          |time                     |Quality flag                                                          |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|antenna_separation            |float32        |                         |Antenna separation                                                    |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|antenna_gain                  |float32        |                         |Antenna gain                                                          |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|program_number                |int32          |                         |Program number                                                        |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|model_number                  |int32          |                         |Model number                                                          |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|sample_duration               |float32        |                         |Sample duration                                                       |s                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|range_resolution              |float32        |chirp_sequence           |Vertical resolution of range                                          |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|dual_polarization             |int32          |                         |Dual polarisation type                                                |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|chirp_start_indices           |int32          |chirp_sequence           |Chirp sequences start indices                                         |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|elevation                     |float32        |time                     |Elevation angle above horizon                                         |degrees                                 |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|noise_threshold               |float32        |                         |Noise filter threshold factor                                         |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|wind_direction                |float32        |time                     |Wind direction                                                        |degrees                                 |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|FFT_window                    |int32          |                         |FFT window type                                                       |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|width                         |float32        |time, range              |Spectral width                                                        |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|transmitted_power             |float32        |time                     |Transmitted power                                                     |W                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|transmitter_temperature       |float32        |time                     |Transmitter temperature                                               |K                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|number_of_spectral_samples    |int32          |chirp_sequence           |Number of spectral samples in each chirp sequence                     |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|file_code                     |int32          |                         |File code                                                             |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|pressure                      |float32        |time                     |Pressure                                                              |Pa                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|receiver_temperature          |float32        |time                     |Receiver temperature                                                  |K                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|input_voltage_range           |int32          |                         |ADC input voltage range (+/-)                                         |mV                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|lwp                           |float32        |time                     |Liquid water path                                                     |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|temperature                   |float32        |time                     |Temperature                                                           |K                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|wind_speed                    |float32        |time                     |Wind speed                                                            |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|calibration_interval          |int32          |                         |Calibration interval in samples                                       |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|integration_time              |float32        |chirp_sequence           |Integration time                                                      |s                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|antenna_diameter              |float32        |                         |Antenna diameter                                                      |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|skewness                      |float32        |time, range              |Skewness of spectra                                                   |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|half_power_beam_width         |float32        |                         |Half power beam width                                                 |degrees                                 |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|number_of_averaged_chirps     |int32          |chirp_sequence           |Number of averaged chirps in sequence                                 |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|rain_rate                     |float32        |time                     |Rain rate                                                             |mm h-1                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

**Variables (BASTA specific):**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

**Variables (MIRA specific):**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|nfft                          |int32          |                         |Number of FFT Points                                                  |count                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|zrg                           |int32          |                         |Number of Range Gates                                                 |count                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|width                         |float32        |time, range              |Spectral width                                                        |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|SNR                           |float32        |time, range              |Signal-to-noise ratio                                                 |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|ldr                           |float32        |time, range              |Linear depolarisation ratio                                           |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|rg0                           |int32          |                         |Number of Lowest Range Gates                                          |count                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|prf                           |int32          |                         |Pulse Repetition Frequency                                            |Hz                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|nave                          |int32          |                         |Number of Spectral Avreages                                           |count                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

Level 1c files
--------------

Categorize file
...............

**Dimensions:**

+------------------------------+
|**Name**                      |
+------------------------------+
|time                          |
+------------------------------+
|height                        |
+------------------------------+
|model_time                    |
+------------------------------+
|model_height                  |
+------------------------------+


**Variables:**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|category_bits                 |int32          |time, height             |Target categorization bits                                            |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|Z_error                       |float32        |time, height             |Error in radar reflectivity factor                                    |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|q                             |float32        |model_time, model_height |Specific humidity                                                     |1                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|Tw                            |float32        |time, height             |Wet-bulb temperature                                                  |K                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|model_height                  |float32        |model_height             |Height of model variables above mean sea level                        |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|insect_prob                   |float32        |time, height             |Insect probability                                                    |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|is_undetected_melting         |int32          |time                     |Presence of undetected melting layer                                  |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|time                          |float32        |time                     |Time UTC                                                              |hours since 2021-06-21 00:00:00         |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|uwind                         |float32        |model_time, model_height |Zonal wind                                                            |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|Z                             |float32        |time, height             |Radar reflectivity factor                                             |dBZ                                     |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|beta_error                    |float32        |                         |Error in attenuated backscatter coefficient                           |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|is_rain                       |int32          |time                     |Presence of rain                                                      |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|beta                          |float32        |time, height             |Attenuated backscatter coefficient                                    |sr-1 m-1                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|lwp_error                     |float32        |time                     |Error in liquid water path                                            |g m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|Z_sensitivity                 |float32        |height                   |Minimum detectable radar reflectivity                                 |dBZ                                     |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|beta_bias                     |int32          |                         |Bias in attenuated backscatter coefficient                            |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|radar_liquid_atten            |float32        |time, height             |Approximate two-way radar attenuation due to liquid water             |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|lidar_wavelength              |float32        |                         |Laser wavelength                                                      |nm                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|pressure                      |float32        |model_time, model_height |Pressure                                                              |Pa                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|lwp                           |float32        |time                     |Liquid water path                                                     |g m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|v                             |float32        |time, height             |Doppler velocity                                                      |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|v_sigma                       |float32        |time, height             |Standard deviation of mean Doppler velocity                           |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|vwind                         |float32        |model_time, model_height |Meridional wind                                                       |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|temperature                   |float32        |model_time, model_height |Temperature                                                           |K                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|latitude                      |float32        |                         |Latitude of site                                                      |degrees_north                           |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|altitude                      |int32          |                         |Altitude of site                                                      |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|model_time                    |float32        |model_time               |Model time UTC                                                        |decimal hours since midnight            |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|Z_bias                        |int32          |                         |Bias in radar reflectivity factor                                     |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|longitude                     |float32        |                         |Longitude of site                                                     |degrees_east                            |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|radar_gas_atten               |float32        |time, height             |Two-way radar attenuation due to atmospheric gases                    |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|quality_bits                  |int32          |time, height             |Data quality bits                                                     |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|height                        |float32        |height                   |Height above mean sea level                                           |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|radar_frequency               |float32        |                         |Radar transmit frequency                                              |GHz                                     |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

Level 2 files
--------------

Classification file
...................

**Dimensions:**

+------------------------------+
|**Name**                      |
+------------------------------+
|time                          |
+------------------------------+
|height                        |
+------------------------------+


**Variables:**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|cloud_top_height_amsl         |float32        |time                     |Height of cloud top above mean sea level                              |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|time                          |float32        |time                     |Time UTC                                                              |hours since 2021-06-21 00:00:00         |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|detection_status              |int32          |time, height             |Radar and lidar detection status                                      |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|latitude                      |float32        |                         |Latitude of site                                                      |degrees_north                           |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|altitude                      |int32          |                         |Altitude of site                                                      |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|cloud_base_height_amsl        |float32        |time                     |Height of cloud base above mean sea level                             |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|longitude                     |float32        |                         |Longitude of site                                                     |degrees_east                            |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|target_classification         |int32          |time, height             |Target classification                                                 |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|height                        |float32        |height                   |Height above mean sea level                                           |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|cloud_base_height_agl         |float32        |time                     |Height of cloud base above ground level                               |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|cloud_top_height_agl          |float32        |time                     |Height of cloud top above ground level                                |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

Drizzle file
............

**Dimensions:**

+------------------------------+
|**Name**                      |
+------------------------------+
|time                          |
+------------------------------+
|height                        |
+------------------------------+


**Variables:**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|drizzle_lwf_error             |float32        |time, height             |Random error in drizzle liquid water flux                             |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|drizzle_lwf_bias              |float32        |                         |Possible bias in drizzle liquid water flux                            |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|drizzle_lwc_bias              |float32        |                         |Possible bias in drizzle liquid water content                         |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|Do_error                      |float32        |time, height             |Random error in drizzle median diameter                               |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|time                          |float32        |time                     |Time UTC                                                              |hours since 2021-06-21 00:00:00         |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|drizzle_N                     |float32        |time, height             |Drizzle number concentration                                          |m-3                                     |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|Do                            |float32        |time, height             |Drizzle median diameter                                               |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|drizzle_N_error               |float32        |time, height             |Random error in drizzle number concentration                          |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|beta_corr                     |float32        |time, height             |Lidar backscatter correction factor                                   |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|drizzle_lwc_error             |float32        |time, height             |Random error in drizzle liquid water content                          |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|v_drizzle                     |float32        |time, height             |Drizzle droplet fall velocity                                         |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|mu_error                      |float32        |                         |Random error in drizzle droplet size distribution shape parameter     |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|S                             |float32        |time, height             |Lidar backscatter-to-extinction ratio                                 |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|drizzle_lwc                   |float32        |time, height             |Drizzle liquid water content                                          |kg m-3                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|Do_bias                       |float32        |                         |Possible bias in drizzle median diameter                              |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|v_air                         |float32        |time, height             |Vertical air velocity                                                 |m s-1                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|mu                            |float32        |time, height             |Drizzle droplet size distribution shape parameter                     |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|latitude                      |float32        |                         |Latitude of site                                                      |degrees_north                           |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|drizzle_lwf                   |float32        |time, height             |Drizzle liquid water flux                                             |kg m-2 s-1                              |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|altitude                      |int32          |                         |Altitude of site                                                      |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|drizzle_retrieval_status      |int32          |time, height             |Drizzle parameter retrieval status                                    |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|longitude                     |float32        |                         |Longitude of site                                                     |degrees_east                            |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|S_error                       |float32        |time, height             |Random error in lidar backscatter-to-extinction ratio                 |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|v_drizzle_error               |float32        |time, height             |Random error in drizzle droplet fall velocity                         |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|height                        |float32        |height                   |Height above mean sea level                                           |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

Iwc file
........

**Dimensions:**

+------------------------------+
|**Name**                      |
+------------------------------+
|time                          |
+------------------------------+
|height                        |
+------------------------------+


**Variables:**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|iwc_inc_rain                  |float32        |time, height             |Ice water content including rain                                      |kg m-3                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|time                          |float32        |time                     |Time UTC                                                              |hours since 2021-06-21 00:00:00         |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|latitude                      |float32        |                         |Latitude of site                                                      |degrees_north                           |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|altitude                      |int32          |                         |Altitude of site                                                      |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|iwc                           |float32        |time, height             |Ice water content                                                     |kg m-3                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|iwc_bias                      |float32        |                         |Possible bias in ice water content, one standard deviation            |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|longitude                     |float32        |                         |Longitude of site                                                     |degrees_east                            |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|iwc_retrieval_status          |int32          |time, height             |Ice water content retrieval status                                    |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|iwc_sensitivity               |float32        |height                   |Minimum detectable ice water content                                  |kg m-3                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|height                        |float32        |height                   |Height above mean sea level                                           |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|iwc_error                     |float32        |time, height             |Random error in ice water content, one standard deviation             |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+

Lwc file
........

**Dimensions:**

+------------------------------+
|**Name**                      |
+------------------------------+
|time                          |
+------------------------------+
|height                        |
+------------------------------+


**Variables:**

+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|**Name**                      |**Data type**  |**Dimension**            |**Long name**                                                         |**Unit**                                |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|lwp_error                     |float32        |time                     |Error in liquid water path                                            |g m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|time                          |float32        |time                     |Time UTC                                                              |hours since 2021-06-21 00:00:00         |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|latitude                      |float32        |                         |Latitude of site                                                      |degrees_north                           |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|altitude                      |int32          |                         |Altitude of site                                                      |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|longitude                     |float32        |                         |Longitude of site                                                     |degrees_east                            |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|lwc_retrieval_status          |int32          |time, height             |Liquid water content retrieval status                                 |                                        |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|lwc_error                     |float32        |time, height             |Random error in liquid water content, one standard deviation          |dB                                      |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|lwc                           |float32        |time, height             |Liquid water content                                                  |kg m-3                                  |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|height                        |float32        |height                   |Height above mean sea level                                           |m                                       |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
|lwp                           |float32        |time                     |Liquid water path                                                     |g m-2                                   |
+------------------------------+---------------+-------------------------+----------------------------------------------------------------------+----------------------------------------+
==========
Quickstart
==========

Processing is easy using CloudnetPy's high level APIs. You only need some
measurement data from your instruments. And if you don't have it, you can
always try `these example files <http://lake.fmi.fi/cloudnet-public/cloudnetpy_test_input_files.zip>`_.

Radar processing
----------------

In the first example we convert a raw METEK MIRA-36 cloud radar file into
Cloudnet netCDF file that can be used in further processing steps.

.. code-block:: python

    from cloudnetpy.instruments import mira2nc
    uuid = mira2nc('raw_mira_radar.mmclx', 'radar.nc', {'name': 'Mace-Head'})

where ``uuid`` is an unique identifier for the generated ``radar.nc`` file.
For more information, see `API reference <api.html#instruments.mira2nc>`__ for this function.

Lidar processing
----------------

Next we convert a raw Jenoptik CHM15k ceilometer (lidar) file into Cloudnet netCDF file
and process the signal-to-noise screened backscatter. Also this converted lidar
file will be needed later.

.. code-block:: python

    from cloudnetpy.instruments import ceilo2nc
    uuid = ceilo2nc('raw_chm15k_lidar.nc', 'lidar.nc', {'name':'Mace-Head', 'altitude':5})

where ``uuid`` is an unique identifier for the generated ``lidar.nc`` file.
For more information, see `API reference <api.html#instruments.ceilo2nc>`__ for this function.

MWR processing
--------------

Processing of multi-channel HATPRO microwave radiometer (MWR) data is not part of CloudnetPy.
Thus, site operators need to run custom processing software to retrieve integrated liquid
water path (LWP) from raw HATPRO measurements.

However, with a 94 GHz RPG cloud radar, a separate MWR instrument is not necessarely
required. RPG radars contain single MWR channel providing a rough estimate
of LWP, which can be used in CloudnetPy. Nevertheless, it is always
recommended to equip a measurement site with a dedicated multi-channel
radiometer if possible.

Model data
----------

Model files needed in the next processing step can be downloaded
from the `Cloudnet http API <https://actris-cloudnet.github.io/dataportal/>`_.
Several models may be available depending on the site and date.
The list of different model models can be found `here <https://cloudnet.fmi.fi/api/models/>`_.

Categorize processing
---------------------

After processing the raw radar and raw lidar files, and acquiring
the model and mwr files, a Cloudnet categorize file can be created.

In the next example we create a categorize file starting from the
``radar.nc`` and ``lidar.nc`` files generated above. The required
``ecmwf_model.nc`` and ``hatpro_mwr.nc`` files are
included in the provided `example input files <http://devcloudnet.fmi.fi/files/cloudnetpy_test_input_files.zip>`_.

.. code-block:: python

   from cloudnetpy.categorize import generate_categorize
   input_files = {
       'radar': 'radar.nc',
       'lidar': 'lidar.nc',
       'model': 'ecmwf_model.nc',
       'mwr': 'hatpro_mwr.nc'
       }
   uuid = generate_categorize(input_files, 'categorize.nc')

where ``uuid`` is an unique identifier for the generated ``categorize.nc`` file.
For more information, see `API reference <api.html#categorize.generate_categorize>`__ for this function.
Note that with a 94 GHz RPG cloud radar, the ``radar.nc`` file can be used as input
for both inputs: ``'radar'`` and ``'mwr'``.


Processing products
-------------------

In the last example we create the smallest and simplest Cloudnet
product, the classification product. The product-generating functions always
use a categorize file as an input.

.. code-block:: python

    from cloudnetpy.products import generate_classification
    uuid = generate_classification('categorize.nc', 'classification.nc')

where ``uuid`` is an unique identifier for the generated ``classification.nc`` file.
Corresponding functions are available for other products
(see :ref:`Product generation`).
Developer's Guide
=================

CloudnetPy is hosted by Finnish Meteorological Institute (FMI) and
will be used to process cloud remote sensing data in the
ACTRIS research infrastructure. We are happy to welcome the cloud remote sensing community
to provide improvements in the methods and their implementations, writing
tests and fixing bugs.

How to contribute
-----------------

Instructions can be found from `CloudnetPy's Github page <https://github.com/actris-cloudnet/cloudnetpy/blob/master/CONTRIBUTING.md>`_.

Testing
-------

To run the CloudnetPy test suite, first
clone the whole repository from `GitHub
<https://github.com/actris-cloudnet/cloudnetpy>`_:

.. code-block:: console

	$ git clone https://github.com/actris-cloudnet/cloudnetpy

Testing environment
...................

Now, create a virtual environment and install pytest and CloudnetPy:

.. code-block:: console

    $ cd cloudnetpy
    $ python3 -m venv venv
    $ source venv/bin/activate
    (venv) $ pip3 install pytest .

Unit tests
..........

.. code-block:: console

    (venv) $ pytest

End-to-end test
...............

.. code-block:: console

    (venv) $ python3 tetsts/e2e_test.py

.. note::

   Cloudnetpy performs relatively complicated scientific processing, converting
   noisy measurement data into higher level products. Most of the
   Cloudnetpy's low-level functions are unit tested, but it is
   difficult to write unambiguous tests for the high-level API calls.
   However, the quality of the processed files can be at least roughly
   checked using CloudnetPy's quality control functions.


Coding guidelines
-----------------

- Use `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ standard.

- Write `Google-style docstrings <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`_.

- Check your code using, e.g., `Pylint <https://www.pylint.org/>`_.
=========================
Installation Instructions
=========================

CloudnetPy can be installed on any computer supporting Python3.6 (or higher).
The actual installation procedure depends on the operating system. The
instructions below are for Ubuntu.

Python Installation
-------------------

.. code-block:: console
		
   $ sudo apt update && sudo apt upgrade
   $ sudo apt install python3 python3-venv python3-pip python3-tk

Virtual Environment
-------------------

Create a new virtual environment and activate it:

.. code-block:: console
		
   $ python3 -m venv venv
   $ source venv/bin/activate


Pip-based Installation
----------------------

CloudnetPy is available from Python Package Index, `PyPI
<https://pypi.org/project/cloudnetpy/>`_.
Use Python's package manager, `pip <https://pypi.org/project/pip/>`_,
to install CloudnetPy package into the virtual environment:

.. code-block:: console
		
   (venv)$ pip3 install cloudnetpy

CloudnetPy is now ready for use from that virtual environment.

.. note::

   CloudnetPy codebase is rapidly developing and the PyPI package does not
   necessarily contain all the latest features and modifications. To get an up-to-date
   version of CloudnetPy, download it directly from `GitHub
   <https://github.com/actris-cloudnet/cloudnetpy>`_.


========
Overview
========

Cloudnet processing scheme
--------------------------

CloudnetPy is a Python package implementing the so-called Cloudnet processing scheme
(`Tukiainen 2020`_). The Cloudnet processing produces vertical profiles of cloud properties
from the ground-based remote sensing measurements. Cloud radar, optical lidar, microwave radiometer
and thermodynamical (model or radiosonde) data are combined to accurately characterize
clouds up to 15 km with high temporal and vertical resolution.

.. figure:: _static/example_data.png
	   :width: 500 px
	   :align: center

           Example of input data used in Cloudnet processing: Radar reflectivity factor (top), mean
           doppler velocity (2nd), lidar backscatter coefficient (3rd),
           and liquid water path from microwave radiometer (bottom).

The measurement and model data are brought into common grid and classified as ice,
liquid, aerosol, insects, and so on. Then, geophysical products such as ice water content
can be retrieved in further processing steps. A more detailed description can be
found in `Illingworth 2007`_ and references in it.

.. note::

    Near real-time Cloudnet data can be accessed at https://cloudnet.fmi.fi.

Statement of need
-----------------

In the forthcoming years, Cloudnet will be one of the key components in `ACTRIS`_ (Aerosol,
Clouds and Trace Gases Research Infrastructure), where the Cloudnet framework will be used
to process gigabytes of cloud remote sensing data per day in near real time. The ACTRIS
RI is now in its implementation phase and aims to be fully operational in 2025.

To fulfill requirements from ACTRIS, a robust, open software that can reliably process
large amounts of data is needed. The CloudnetPy software package is aimed to perform operational
ACTRIS cloud remote sensing processing, providing quality controlled data products from
around 15 measurement sites in Europe. Unlike the original proprietary Cloudnet software,
CloudnetPy is open source and includes tests, documentation and user-friendly API that
the research community can use to further develop the existing methods and to create
new products.

.. _Tukiainen 2020: https://doi.org/10.21105/joss.02123
.. _Illingworth 2007: https://journals.ametsoc.org/doi/abs/10.1175/BAMS-88-6-883
.. _ACTRIS: http://actris.eu/

.. important::

   CloudnetPy is a rewritten Python implementation of the original Matlab processing prototype code.
   CloudnetPy features several revised methods and bug fixes, open source codebase,
   netCDF4 file format and extensive documentation.

Operational Cloudnet processing will have CloudnetPy in its core but additionally include a
calibration database and comprehensive quality control / assurance procedures:

.. figure:: _static/CLU_workflow.png
	   :width: 650 px
	   :align: center

           Workflow of operational Cloudnet processing in ACTRIS.


See also:

- Cloudnet data portal: https://cloudnet.fmi.fi/
- CloudnetPy source: https://github.com/actris-cloudnet/cloudnetpy
- ACTRIS home: http://actris.eu/
- ACTRIS data portal: http://actris.nilu.no/
.. cloudnet documentation master file, created by
   sphinx-quickstart on Wed Dec  5 21:38:46 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

========================
CloudnetPy documentation
========================

Welcome! This is the documentation for CloudnetPy, the Python
implementation of the Cloudnet processing scheme.


.. toctree::
   :maxdepth: 4

   overview
   installation
   quickstart
   api
   fileformat
   guide



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
API reference
=============


High-level functions
--------------------

CloudnetPy's high-level functions provide a simple mechanism to process
cloud remote sensing measurements into Cloudnet products. A full processing
goes in steps. Each step produces a file which used as an input for the
next step.

Raw data conversion
...................

Different Cloudnet instruments provide raw data in various formats (netCDF, binary, text)
that first need to be converted into homogeneous Cloudnet netCDF files
containing harmonized units and other metadata. This initial processing step
is necessary to ensure that the subsequent processing steps work with
all supported instrument combinations.

.. autofunction:: instruments.mira2nc

.. autofunction:: instruments.rpg2nc

.. autofunction:: instruments.basta2nc

.. autofunction:: instruments.ceilo2nc

.. autofunction:: instruments.pollyxt2nc

.. autofunction:: instruments.hatpro2nc

.. autofunction:: instruments.disdrometer2nc


The categorize file
...................

The categorize file concatenates all input data into common
time / height grid.

.. autofunction:: categorize.generate_categorize


Product generation
..................

Starting from the categorize file, several geophysical products can be
generated.

.. autofunction:: products.generate_classification

.. autofunction:: products.generate_iwc

.. autofunction:: products.generate_lwc

.. autofunction:: products.generate_drizzle


Visualizing results
...................

CloudnetPy offers an easy-to-use plotting interface:

.. autofunction:: plotting.generate_figure

There is also possibility to compare CloundetPy files with the
Matlab-processed legacy files
(tagged "legacy" in the `Cloudnet data portal <https://cloudnet.fmi.fi>`_):

.. autofunction:: plotting.compare_files


Categorize modules
------------------

Categorize is CloudnetPy's subpackage. It contains
several modules that are used when creating the Cloudnet
categorize file.


datasource
..........

.. automodule:: categorize.datasource
   :members:

radar
.....

.. automodule:: categorize.radar
   :members:

lidar
.....

.. automodule:: categorize.lidar
   :members:

mwr
...

.. automodule:: categorize.mwr
   :members:

model
.....

.. automodule:: categorize.model
   :members:

classify
........

.. automodule:: categorize.classify
   :members:

melting
.......

.. automodule:: categorize.melting
   :members:

freezing
........

.. automodule:: categorize.freezing
   :members:


falling
.......

.. automodule:: categorize.falling
   :members:


insects
.......

.. automodule:: categorize.insects
   :members:


atmos
.....

.. automodule:: categorize.atmos
   :members:


droplet
.......

.. automodule:: categorize.droplet
   :members:


Products modules
----------------

Products is CloudnetPy's subpackage. It contains
several modules that correspond to different Cloudnet
products.

classification
..............

.. automodule:: products.classification
   :members:


iwc
...

.. automodule:: products.iwc
   :members:


lwc
...

.. automodule:: products.lwc
   :members:


drizzle
.......

.. automodule:: products.drizzle
   :members:


product_tools
.............

.. automodule:: products.product_tools
   :members:


Misc
----

Documentation for various modules with low-level
functionality.

concat_lib
..........

.. automodule:: concat_lib
   :members:


utils
.....

.. automodule:: utils
   :members:


cloudnetarray
.............

.. automodule:: cloudnetarray
   :members:


output
......

.. automodule:: output
   :members:

