---
title: 'Virgo: A Versatile Spectrometer for Radio Astronomy'
tags:
  - Python
  - astronomy
  - radio astronomy
  - digital signal processing
  - astrophysics
authors:
  - name: Apostolos Spanakis-Misirlis
    orcid: 0000-0001-6928-6877
    affiliation: 1
  - name: Cameron L. Van Eck
    affiliation: 2
    orcid: 0000-0002-7641-9946
  - name: E.P. Boven
    affiliation: 3
affiliations:
 - name: Department of Informatics, University of Piraeus, Greece
   index: 1
 - name: Dunlap Institute for Astronomy and Astrophysics, University of Toronto, 50 St. George Street, Toronto, ON M5S 3H4, Canada
   index: 2
 - name: C.A. Muller Radioastronomie Station (CAMRAS), Dwingeloo, the Netherlands
   index: 3
date: 20 January 2021
bibliography: paper.bib

---

# Introduction

For the past few decades, radio astronomy has been a rapidly developing area of
observational astronomy. This is due to the fact that a variety of celestial objects emit
electromagnetic radiation at radio wavelengths, which has led to the
development of radio telescopes capable of revealing the otherwise-hidden
astrophysical properties of the universe. An important requirement that makes radio
astronomy observations and analysis possible is an appropriate software pipeline
compatible with the spectrometers with which radio observatories are equipped. In
this work, we present `Virgo`: a versatile software solution for radio telescopes.

# Statement of Need

`Virgo` is a Python package for the acquisition, processing and analysis of
data from radio telescopes. It is an easy-to-use open-source spectrometer and
radiometer based on the GNU Radio framework (https://www.gnuradio.org), and is conveniently
applicable to any radio telescope working with a GNU Radio-supported software-defined
radio (SDR; a radio receiver architecture where some conventional hardware-based
steps are replicated in software, @Dillinger:2003). Although in-house acquisition
solutions have been independently created by many, composing their own software
for simple observations, these are not comparable in terms of installation count,
ease of use, cross-platform support, functionality, reproducibility, maintenance
and level of documentation that this package has to offer.

Designed to be used by students, educators and amateurs in the field of radio
astronomy, `Virgo` has already been adopted by a number of small and
large-aperture radio telescopes, permitting both spectral and continuum
observations with great success. These instruments include the ISEC TLM-18 (18m),
the ACRO RT-320 (3.2m), the JRT (1.9m), and the PICTOR Telescope (1.5m), among others.

Although the hardware aspect of a radio telescope is generally handled by newcomers
with relative ease, the skill set needed to integrate a complete software pipeline to
support observations is not something most users are equipped with. `Virgo` tackles
this problem by providing non-experts with a tool to collect and interpret data from
radio telescopes, without requiring expertise in digital signal processing and
software engineering. An example use case is classroom experiments in which students
build a small-aperture antenna connected to a low-noise amplifier followed by an SDR,
and with the help of `Virgo`, obtain data to map out the galactic distribution of
neutral hydrogen and/or derive the rotation curve of the Milky Way. An example
observation of the 21-cm hydrogen line acquired and processed with `Virgo` is shown
in \autoref{fig:example}. The package's versatility also provides a convenient solution for
researchers wishing to rapidly deploy low-cost radio telescopes with commercial hardware.


![Observation of galactic clouds of neutral hydrogen toward the constellation of Cygnus ($\alpha = 20^{\mathrm{h}}$, $\delta = 40^{\circ}$, $l = 77^{\circ}$, $b = 3^{\circ}$), observed by the TLM-18 Telescope in New Jersey, U.S. with `Virgo`. The average spectrum (top left), the calibrated spectrum (top center), the dynamic spectrum (top right) and the time series along with the total power distribution (bottom) are all plotted by the software automatically.\label{fig:example}](example.pdf)

# Features

One of the key features of `Virgo` is its implementation of a four-tap weighted overlap-add (WOLA) Fourier transform (FT)
spectrometer, offering a significant reduction in spectral leakage compared to
a simple FT filterbank spectrometer that does not make use of the WOLA method, with a minimal
increase in computational requirements [@Crochiere1996]\footnote{The package also supports
a plain FT filterbank pipeline for observatories with limited computational resources.}. In addition to
its data-acquisition functionality that performs data reduction by time-averaging
spectra in real time, `Virgo` also carries out automated analysis of the recorded
samples. The time-averaged spectrum, the calibrated spectrum, the dynamic spectrum
(waterfall), the time series (power vs time) and the total power distribution of the
observation are all automatically computed and plotted with the help of the `Numpy`
[@Harris2020] and `Matplotlib` [@Hunter:2007] packages.

Because of the nature of the RF instrumentation that radio telescopes are equipped with, the spectra acquired by SDRs
have an unwanted frequency-dependant sensitivity, also known as the bandpass shape.
In general, this frequency response makes it difficult to distinguish true signals,
originating from the sky, from instrumentation artifacts. For that reason, `Virgo`
performs bandpass calibration by taking the ratio of the observed spectrum over
the calibration spectrum. However, because this ratio is arbitrarily scaled (due to
the difference in the noise floor levels), the power axis is automatically rescaled
to units of signal-to-noise ratio. In case the resulted spectrum has an unwanted
slope due to e.g. inconsistent conditions between the calibration and the observation,
the software can also automatically correct poorly-calibrated spectra using linear
regression.

Furthermore, `Virgo` supports optional median operations, both
in the frequency and time domain, for the suppression of narrowband and/or
short-duration radio frequency interference (RFI), while allowing the user to export
the raw observation data as a FITS [@Pence2010] or csv-formatted file. Frequencies contaminated
with RFI may also be stamped out with the software's built-in channel masking
capability.

Lastly, the package may also be used for the detection of giant pulses (irregularly intense
bursts of radio emission by pulsars). Due to the frequency-dependent dispersion
caused by the plasma distribution in the interstellar medium (ISM), observed pulses
naturally appear smeared in time, depending on the dispersion measure (integrated column density of
free electrons from the observer to the source). To prevent implied degradations of
the signal-to-noise ratio, incoherent dedispersion is optionally applied to dynamic spectra
of pulsar observations, compensating for the unwanted dispersive effects of the ISM.

By additionally providing the observer with an important set of utilities, `Virgo` also
makes for a great tool for planning (radio) observations. This includes the ability to
compute the position of astronomical sources in the sky for a given date (see \autoref{fig:predict}),
and conversely, to estimate the right ascencion and declination given the observer's coordinates along with
the altitude and azimuth the telescope is pointing to, with the help of the `Astropy` package [@astropy:2013; @astropy:2018].
On top of that, the package comes with a basic set of calculators for quickly carrying out a variety of
computations involving the theoretical sensitivity and performance of a given instrument.

![Example prediction of the location of the Cygnus A radio galaxy (3C 405) in the celestial sphere of the observer.\label{fig:predict}](predict.pdf)

Likewise, the software provides a handy tool for retrieving HI profiles based on the
Leiden/Argentine/Bonn (LAB) Survey of Galactic HI [@Kalberla2005]. These spectra (see \autoref{fig:profile}
for an example) can be correlated with the integrated 21 cm all-sky map previewer shown
in \autoref{fig:map}.

![Sample HI profile ($\alpha = 20^{\mathrm{h}}30^{\mathrm{m}}$, $\delta = 45^{\circ}$) obtained with the package's `virgo.simulate()` function.\label{fig:profile}](profile.pdf)

![21 cm all-sky map rendered by the software. The red dot indicates the position of the telescope's beam in the sky, provided by the user.\label{fig:map}](map.pdf)

Moreover, the package comes with an integrated frequency-domain RFI measurement pipeline,
allowing observers to rapidly carry out a survey outlining the compatibility of the
telescope's environment with radio observation standards.

Last but not least, the software's modularity allows users to effortlessly integrate
`Virgo`'s functionalities into other pipelines, permitting a variety of automation
applications invloving the acquisition and/or processing of telescope data.

# Example Usage

`Virgo` can either be called directly from the command line using e.g.,

`virgo -rf 10 -if 20 -bb 20 -f 1420e6 -b 5e6 -c 2048 -t 1 -d 60 -s 10 -n 20 -m 35 -r 1420.4057517667e6 -C calibration.dat -W obs.fits`,

or imported and used as a package:
```python
import virgo

# Define observation parameters
obs = {
    'dev_args': '',
    'rf_gain': 10,
    'if_gain': 20,
    'bb_gain': 20,
    'frequency': 1420e6,
    'bandwidth': 5e6,
    'channels': 2048,
    't_sample': 1,
    'duration': 60
}

# Check source position
virgo.predict(lat=39.83, lon=-74.87, source='Cas A', date='2020-12-26')

# Begin data acquisition in 10 sec
virgo.observe(obs_parameters=obs, obs_file='observation.dat', start_in=10)

# Analyze data, mitigate RFI and export the data as a FITS file
virgo.plot(obs_parameters=obs, n=20, m=35, f_rest=1420.4057517667e6,
           obs_file='observation.dat', cal_file='calibration.dat',
           rfi=[1419.2e6, 1419.3e6], waterfall_fits='obs.fits',
           slope_correction=True, plot_file='plot.png')
```

# Future Work

As `Virgo` gets adopted by more and more radio telescopes, the need for expanding the
software's capabilities grow. Additional features that have been proposed include but
are not limited to a more robust and intelligent system for the detection and
mitigation of RFI, and the support for a data acquisition/analysis pipeline for
pulsar astronomy, both of which offer the potential of appealing to a broader user
audience.

# References
# Virgo: A Versatile Spectrometer for Radio Astronomy
<p align="center">
  <img src="https://i.imgur.com/lH2OOTd.png?raw=true" alt="Virgo Spectrometer"/>
</p>

<p align="center">
  <a href="https://joss.theoj.org/papers/612c2634c1dd83749e95a93449740861"><img src="https://joss.theoj.org/papers/612c2634c1dd83749e95a93449740861/status.svg"></a>
  <img src="https://img.shields.io/badge/python-2.7%20%7C%203.x-green"/>
  <img src="https://img.shields.io/pypi/v/astro-virgo"/>
  <img src="http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat"/>
  <img src="https://img.shields.io/github/license/0xCoto/Virgo?color=yellow"/>
</p>


## About Virgo
**Virgo** is an easy-to-use **open-source** spectrometer and radiometer based on [Python](https://www.python.org) and [GNU Radio](https://wiki.gnuradio.org) (GR) that is conveniently applicable to any radio telescope working with a GR-supported software-defined radio (SDR). In addition to data acquisition, Virgo also carries out automated analysis of the recorded samples, producing an **averaged spectrum**, a **calibrated spectrum**, a **dynamic spectrum (waterfall)**, a **time series (power vs time)** and a **total power distribution** plot.

Lastly, an important set of utilities is provided to observers, making the package for a great tool for planning (radio) observations, estimating the system sensitivity of an instrument, and many more.

## Statement of Need

Designed to be used by students, educators and amateurs in the field of radio astronomy, Virgo has already been adopted by a number of small and large-aperture radio telescopes, permitting both spectral and continuum observations with great success. These instruments include the ISEC TLM-18 (18m), the ACRO RT-320 (3.2m), the JRT (1.9m), and the PICTOR Telescope (1.5m), among others.

Although the hardware aspect of a radio telescope is generally handled by newcomers with relative ease, the skill set needed to integrate a complete software pipeline to support observations is not something most users are equipped with. Virgo tackles this problem by providing non-experts with a tool to collect and interpret data from radio telescopes, without requiring expertise in digital signal processing and software engineering. An example use case is classroom experiments in which students build a small-aperture antenna connected to a low-noise amplifier followed by an SDR, and with the help of Virgo, obtain data to map out the galactic distribution of neutral hydrogen and/or derive the rotation curve of the Milky Way.

### Key Features

- 4-tap weighted overlap-add (WOLA) Fourier transform spectrometer
- - Reduced FFT sidelobes
- - Plain FT filterbank pipeline also supported for observatories with limited computational resources
- Adjustable SDR parameters
- - Device arguments
- - RF/IF/BB Gain
- Header file
- - Observation parameters automatically passed to corresponding `.header` file
- - Includes logged MJD (at observation *t<sub>0</sub>*)
- Spectral line support
- - Spectrum calibration
- - - *y* axis is automatically rescaled to S:N units with line masking
- - - Optional automatic slope correction (based on linear regression) for poorly-calibrated spectra
- - Supports median operation for RFI mitigation on the frequency-domain (adjustable *n*-factor)
- - RFI channel masking
- - Adjustable *f*<sub>rest</sub> for observation of any spectral line (not just HI)
- - Secondary axes for relative velocity automatically adjusted accordingly
- - Prevention against strong narrowband RFI rescaling subplot
- - The average spectra, calibration spectra and calibrated spectra are optionally saved as a `csv` file for further analysis
- Continuum support
- - Supports median operation for time-varying RFI mitigation (adjustable *n*-factor)
- - Total power distribution (histogram) displayed, both for raw and clean data
- - - Best Gaussian fits computed automatically
- - Prevention against strong short-duration RFI rescaling subplot
- - Time series optionally saved as a `csv` file for further analysis
- Pulsars
- - Incoherent dedispersion support for giant pulse search (and FRB follow-up, assuming DM is known)
- Dynamic spectrum (waterfall)
- - Optionally saved as a `FITS` file for further advanced/custom analysis
- Decibel support
- - Power units optionally displayed in dB
- Observation planning toolkit
- - Predict source altitude & azimuth vs time
- - Quickly convert galactic to equatorial and Alt/Az to RA/Dec
- - Plot telescope position on the 21 cm all-sky survey
- - Simulate 21 cm profiles based on the LAB HI survey
- Basic calculation toolkit for system sensitivity & performance. Computes:
- - Antenna gain (in dBi, linear or K/Jy)
- - Effective aperture
- - Half-power beamwidth
- - Noise figure to noise temperature and vice versa
- - Antenna gain-to-noise-temperature (G/T)
- - System equivalent flux density (SEFD)
- - Radiometer equation (S:N estimation)
- Built-in tool for conducting rapid RFI surveys
- Argument-parsing support
- Works directly from the command line (`virgo -h`), or as a Python module

---

## Telescopes based on the Virgo Spectrometer
- ISEC TLM-18 Telescope (18m)
- ACRO RT-320 (3.2m)
- SALSA Vale Telescope (2.3m) [potentially soon, but already tested]
- SALSA Brage Telescope (2.3m) [potentially soon, but already tested]
- JRT (1.9m)
- PICTOR Telescope (1.5m)
- NanoRT Telescope (15cm)
- and more!

### Example Observation
<p align="center">
  <img src="https://i.imgur.com/ROPPWza.png" alt="Example Observation"/>
</p>
Observation of galactic clouds of neutral hydrogen toward the constellation of Cygnus (α = 20h, δ = 40° , l = 77° , b = 3°), observed by the TLM-18 Telescope in New Jersey, U.S. with Virgo. The average spectrum (top left), the calibrated spectrum (top center), the dynamic spectrum (top right) and the time series along with the total power distribution (bottom) are all plotted by the software automatically.

### Example Source Location Prediction
<p align="center">
  <img src="https://i.imgur.com/jnGJEvQ.png" alt=""/>
</p>

### Example HI Profile Retrieval
<p align="center">
  <img src="https://i.imgur.com/HHSkDJM.png" alt="Example HI Profile Retrieval"/>
</p>

### Example HI Map
<p align="center">
  <img src="https://i.imgur.com/bvg4r4c.png" alt="Example HI map plot"/>
</p>
The red dot indicates the position of the telescope's beam in the sky.

## Data Acquisition Flowgraph
**Virgo** is a **four-tap WOLA Fourier transform** spectrometer. The raw I/Q samples are processed in real time using GNU Radio, with the amount of data stored to file being drastically reduced for further analysis. The following flowgraph handles the acquisition and early-stage processing of the data:

![alt text](https://i.imgur.com/5tR7WjL.png "Data Acquisition Flowgraph")

### Example radio map acquired and processed with the help of Virgo (PICTOR Northern HI Survey)
![alt text](https://i.imgur.com/pYgMAhW.png "PICTOR HI Survey")

## Installation
To use **Virgo**, make sure **[Python](https://www.python.org/)** and **[GNU Radio](https://wiki.gnuradio.org/index.php/InstallingGR)** (with **[gr-osmosdr](https://osmocom.org/projects/gr-osmosdr/wiki)**) are installed on your machine.

For Debian/Ubuntu and derivates, the installation is straightforward:

```
sudo apt install gnuradio gr-osmosdr
```

**Note:** The `GNU Radio` and `gr-osmosdr` dependencies are only required for **acquiring data with the necessary hardware** (software-defined radio). They are not required for planning observations, analyzing data, running calculations or any other functionalities provided by the package. For more information, please refer to the [Dependencies section](https://virgo.readthedocs.io/en/latest/installation.html#dependencies).

Once Python and GNU Radio are installed on your system, run

```
pip install astro-virgo
```

## Documentation
To learn how to use Virgo, please read through the documentation **[here](https://virgo.readthedocs.io/en/latest/)**.

If you believe something is not clarified in the documentation page, you are encouraged to **[create an issue](https://github.com/0xCoto/Virgo/issues/new)** (or send me an [e-mail](mailto:0xcoto@protonmail.com) and I'll be happy to help.

## Contributing
If you wish to contribute to the package (either with code, ideas or further documentation), please read through the **[Contributor Guidelines](https://github.com/0xCoto/Virgo/blob/master/docs/contributing.md)**.

## Credits
**Virgo** was created by **[Apostolos Spanakis-Misirlis](https://www.github.com/0xCoto/)**.

**Contact:** [0xcoto@protonmail.com](mailto:0xcoto@protonmail.com)

---

Special thanks to **Dr. Cameron Van Eck**, **Paul Boven** and **Dr. Cees Bassa** for their valuable contributions.
# Contributing

When contributing to Virgo, please first discuss the change you wish to make via an issue,
email, or any other method with the owner of this repository before making a change.

Please note the Code of Conduct. We hope you follow it in all your interactions with the project.

## Pull Requests

1. Update the README.md with details of changes to the interface, this includes new environment 
   variables, exposed ports, useful file locations and container parameters.
2. Increase the version numbers in any examples files and the README.md to the new version that this
   Pull Request would represent. The versioning scheme Virgo uses is [SemVer](http://semver.org/).

## Code of Conduct

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

### Standards

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

### Responsibilities

The maintainer of the project is responsible for clarifying the standards of acceptable
behavior and is expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

The project maintainer has the right and responsibility to remove, edit, or
reject comments, commits, code, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

### Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project maintainer at 0xcoto@protonmail.com. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project maintainer is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

### Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
