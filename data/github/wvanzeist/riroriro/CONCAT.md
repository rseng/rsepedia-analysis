# Contributing Guidelines

Contributions to Riroriro are welcome. If you have an issue, a code contribution or a documentation contribution, thanks for helping to improve Riroriro!

## Issues

When creating an issue report, please be as specific as possible when describing how to reproduce an issue, and include both the intended/expected result and what you are actually getting.

## Code contributions

If you are performing a bugfix or adding a feature, please fork the repository, create a dedicated branch for your contribution and create a pull request when you are done.

## Documentation contributions

Riroriro has documentation written in reStructuredText, which can be found in the docs/ directory. If you have an improvement for the documentation, please put it on a dedicated branch as with code contributions. If you are adding a new feature in code contributions, you are encouraged to add documentation about it in the docs/ directory.# riroriro

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4588070.svg)](https://doi.org/10.5281/zenodo.4588070)

**Riroriro** is a set of Python modules containing functions to simulate the gravitational waveforms of mergers of black holes and/or neutron stars, and calculate several properties of these mergers and waveforms, specifically relating to their observability by gravitational wave detectors. Riroriro combines areas covered by previous gravitational wave models (such as gravitational wave simulation, SNR calculation, horizon distance calculation) into a single package with broader scope and versatility in Python, a programming language that is ubiquitous in astronomy. Aside from being a research tool, Riroriro is also designed to be easy to use and modify, and it can also be used as an educational tool for students learning about gravitational waves.

The modules “inspiralfuns”, “mergerfirstfuns”, “matchingfuns”, “mergersecondfuns” and “gwexporter”, in that order, can be used to simulate the strain amplitude and frequency of a merger gravitational waveform. The module “snrcalculatorfuns” can compare such a simulated waveform to a detector noise spectrum to calculate a signal-to-noise ratio (SNR) for that signal for that detector. The module “horizondistfuns” calculates the horizon distance of a merger given its waveform, and the module “detectabilityfuns” evaluates the detectability of a merger given its SNR.

Riroriro is installable via pip:

    pip install riroriro

More information on the pip installation can be found here: https://pypi.org/project/riroriro/

Tutorials for Riroriro can be found here: https://github.com/wvanzeist/riroriro_tutorials

Full documentation of each of the functions of Riroriro can be found here: https://wvanzeist.github.io/

Riroriro is one of several Python packages associated with **BPASS** (Binary Population And Spectral Synthesis), a suite of programs that simulates the evolution of a population of binary and single-star systems from a wide range of initial conditions. Each of these associated packages are named after native animals of New Zealand. The riroriro (*Gerygone igata*, also known as the grey warbler) is a small bird that can be recognised by its distinctive melodious call but is rarely seen, similarly to how black hole binary mergers are detected by their gravitational wave signals rather than visually.

The central website of BPASS, which also contains links to related programs, can be found here: https://bpass.auckland.ac.nz

## Paper

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02968/status.svg)](https://doi.org/10.21105/joss.02968)

A paper describing Riroriro has been published in the Journal of Open Source Software. If you use Riroriro in your work, please cite this paper! https://doi.org/10.21105/joss.02968

    @ARTICLE{2021JOSS....6.2968V,
           author = {{van Zeist}, Wouter G.~J. and {Stevance}, H{\'e}lo{\"i}se F. and {Eldridge}, J.~J.},
            title = "{Riroriro: Simulating gravitational waves and evaluating their detectability in Python}",
          journal = {The Journal of Open Source Software},
         keywords = {Python, neutron stars, astronomy, gravitational waves, black holes, General Relativity and Quantum Cosmology, Astrophysics - High Energy Astrophysical Phenomena},
             year = 2021,
            month = mar,
           volume = {6},
           number = {59},
              eid = {2968},
            pages = {2968},
              doi = {10.21105/joss.02968},
    archivePrefix = {arXiv},
           eprint = {2103.06943},
     primaryClass = {gr-qc},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2021JOSS....6.2968V},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
---
title: 'Riroriro: Simulating gravitational waves and evaluating their detectability in Python'
tags:
  - Python
  - astronomy
  - gravitational waves
  - black holes
  - neutron stars
authors:
  - name: Wouter G. J. van Zeist
    affiliation: 1
  - name: Héloïse F. Stevance
    affiliation: 1
  - name: J. J. Eldridge
    affiliation: 1
affiliations:
  - name: Department of Physics, University of Auckland, New Zealand
    index: 1
date: 2021
bibliography: paper.bib
---

# Summary

`Riroriro` is a Python package to simulate the gravitational waveforms of binary mergers of black holes and/or neutron stars, and calculate several properties of these mergers and waveforms, specifically relating to their observability by gravitational wave detectors.

The gravitational waveform simulation of `Riroriro` is based upon the methods of @buskirk2019, a paper which describes a computational implementation of an earlier theoretical gravitational waveform model by @huerta2017, using post-Newtonian expansions and an approximation called the implicit rotating source to simplify the Einstein field equations and simulate gravitational waves. `Riroriro`'s calculation of signal-to-noise ratios (SNR) of gravitational wave events is based on the methods of @barrett2018, with the simpler gravitational wave model `Findchirp` [@findchirp] being used for comparison and calibration in these calculations.

# Statement of Need

The field of gravitational wave astronomy has been especially active since the first observation of gravitational waves was announced in 2016 [@gw150914discovery]. Observations of gravitational waves from binary mergers can provide unique information about their progenitors and stellar populations, especially when combined with electromagnetic observations in the field called multi-messenger astronomy. A major factor in the successful detection and analysis of gravitational wave signals is the creation of simulations of such signals which observed data can be compared to. Because of this, multiple gravitational wave models have been created over the years. In particular, the gravitational wave observatories LIGO and Virgo have created their own models to use as templates in gravitational wave searches, with the main software for this being `LALSuite` [@lalsuite]. Various research groups have created waveform models, with some examples of recent sophisticated waveform models being `IMRPhenomXPHM` [@imrphenom] and `SEOBNRv4PHM` [@seobnr], which are both also available in `LALSuite`.

We have not tested if the waveform model of `Riroriro` that is based on @huerta2017 and @buskirk2019 is accurate enough to use for parameter estimation of detected gravitational wave transients as this was not within the scope of our project; it is likely that accurate parameter estimation requires careful modelling of the ringdown phase, especially for the most massive mergers. However, we use a level of accuracy adequate for our aim of modelling the detectability of gravitational wave transients predicted by stellar population syntheses. Furthermore, the code of `Riroriro` is structured and commented in such a way that each step in the process of the simulation is individually identifiable and modifiable by users. Users could also substitute in functions from other algorithms or even use the detectability modules on waveforms from other sources, as long as the user puts these in Riroriro’s format.

`Riroriro` combines areas covered by previous models (such as gravitational wave simulation, SNR calculation, horizon distance calculation) into a single package with broader scope and versatility in Python, a programming language that is ubiquitous in astronomy. Aside from being a research tool, `Riroriro` is also designed to be easy to use and modify, and it can also be used as an educational tool for students learning about gravitational waves.

# Features

Features of `Riroriro` include:

- Simulating the gravitational waveform signal from a binary merger of two black holes, two neutron stars or a black hole and a neutron star and outputting the data of this signal in terms of frequency and strain amplitude.
- Using a gravitational wave output and given a detector noise spectrum (such spectra are made publicly available by LIGO), calculating the signal-to-noise ratio (SNR) of the signal at a given distance assuming optimal alignment.
- Calculating the horizon distance (maximum distance at which an event could be observed) for a gravitational wave model and a given detector.
- Given the optimal-alignment SNR of an event, evaluating its detectability, the probability that the event would be detected with a SNR above the commonly used threshold of 8, if the alignment would be arbitrary. These results could then be combined with population synthesis calculations to estimate how many of the predicted mergers would be detected.

In addition, to help users get started with `Riroriro`, we have created Jupyter Notebook tutorials ([see here](https://github.com/wvanzeist/riroriro_tutorials)).

# Research

`Riroriro` has been used for research in conjunction with `BPASS`, a suite of computer programs that simulates the evolution of a population of binary and single-star systems from a wide range of initial conditions and predicts their electromagnetic spectral emission [@bpass1; @bpass2]. There is also a Python interface for `BPASS` called `Hoki` [@hoki]. This research took rates of formation of merging systems from `BPASS` and then evaluated the detectability of the gravitational wave signals from those systems using `Riroriro`. This was done to obtain predictions of the rates at which gravitational waves of different types would be expected to be observed, which can then be directly compared to those events found by the LIGO/Virgo gravitational wave observatories [@massdistribution].

# Acknowledgments

HFS and JJE acknowledge support from the University of Auckland and also the Royal Society of New Zealand Te Apārangi under the Marsden Fund​.

# References
Documentation
=============

This is the documentation for the functions of riroriro. Tutorials for the usage of these functions can be found at https://github.com/wvanzeist/riroriro_tutorials

.. toctree::
  :maxdepth: 2
  :caption: Contents:

  riroriro/index.rst
  riroriro/inspiralfuns.rst
  riroriro/mergerfirstfuns.rst
  riroriro/matchingfuns.rst
  riroriro/mergersecondfuns.rst
  riroriro/gwexporter.rst
  riroriro/snrcalculatorfuns.rst
  riroriro/horizondistfuns.rst
  riroriro/detectabilityfuns.rst
  riroriro/tests.rst
***************
mergerfirstfuns
***************

This is the documentation for the mergerfirstfuns module, which consists of the pre-matching parts of the procedure for simulating the merger/ringdown portions of gravitational waves from binary black holes, collected into modular functions.

quasi_normal_modes
==================

``quasi_normal_modes(eta)``

Calculation of the final spin and quasi-normal mode factor used in the
calculation of angular frequency for the merger/ringdown waveform, based on
Buskirk et al. (2019) equations 19 and 20.

Parameters
----------
eta: float
    Symmetric mass ratio of the binary, can be obtained from
    get_M_and_eta() in inspiralfuns.

Returns
-------
(sfin,wqnm): tuple of floats
    The first constant is the final spin, the second is the quasi-normal
    mode factor used in subsequent calculations.

gIRS_coefficients
=================

``gIRS_coefficients(eta,sfin)``

Calculation of several gIRS (generic implicit rotating source)-related
coefficients used in the calculation of angular frequency for the merger/
ringdown waveform, based on Buskirk et al. (2019) Appendix C.

Parameters
----------
eta: float
    Symmetric mass ratio of the binary, can be obtained from
    get_M_and_eta() in inspiralfuns.
sfin: float
    Final spin value, from quasi_normal_modes().
    
Returns
-------
(alpha,b,C,kappa): tuple of floats
    Four gIRS-related constants used in subsequent calculations.
    (NOTE: alpha is not used by anything in mergerfirstfuns but *is* used
    in mergersecondfuns.)
    
merger_freq_calculation
=======================

``merger_freq_calculation(wqnm,b,C,kappa)``

Calculation of orbital angular frequency for the merger/ringdown portion,
based on Buskirk et al. (2019) equations 17 and 18.

Parameters
----------
wqnm: float
    Quasi-normal mode factor, from quasi_normal_modes().
b: float
    A gIRS coefficient, from gIRS_coefficients().
C: float
    A gIRS coefficient, from gIRS_coefficients().
kappa: float
    A gIRS coefficient, from gIRS_coefficients().

Returns
-------
[fhat,m_omega]: list of lists of floats
    First list is the values over time of a sort of frequency parameter
    called fhat (f^), second list is the angular frequency.
    
fhat_differentiation
====================

``fhat_differentiation(fhat)``

Calculation of derivative of fhat used by amplitude calculations in
mergersecondfuns.

Parameters
----------
fhat: list of floats
    Values of a sort of frequency parameter called fhat (f^) over time,
    from merger_freq_calculation().
    
Returns
-------
fhatdot: list of floats
    Values of the time-derivative of fhat over time.
    
merger_time_conversion
======================

``merger_time_conversion(M)``

Calculating times in real units corresponding to the times in geometric
units used by other merger/ringdown functions.

Parameters
----------
M: float
    Total mass of the binary, can be obtained from get_M_and_eta() in
    inspiralfuns.
    
Returns
-------
m_time: list of floats
    The list of timesteps used by other merger/ringdown functions, but in
    seconds instead of geometric units.
****************
mergersecondfuns
****************

This is the documentation for the mergersecondfuns module, which consists of the post-matching parts of the procedure for simulating the merger/ringdown portions of gravitational waves from binary black holes, collected into modular functions.

merger_phase_calculation
========================

``merger_phase_calculation(min_switch_ind,final_i_index,i_phase,m_omega)``

Calculation of the orbital phase for the merger/ringdown portion, based
on Buskirk et al. (2019) equation 21.

Parameters
----------
min_switch_ind: int
    The index in the merger/ringdown data where the switch from inspiral to
    merger/ringdown should occur, from min_switch_ind_finder() in
    matchingfuns.
final_i_index: int
    The last index in the inspiral data before the switch to the merger/
    ringdown, from final_i_index_finder() in matchingfuns.
i_phase: list of floats
    Values of orbital phase at each timestep for the inspiral portion, from
    inspiral_phase_freq_integration() in inspiralfuns.
m_omega: list of floats
    Values of angular frequency over time for the merger/ringdown portion,
    from merger_freq_calculation() in mergerfirstfuns.
    
Returns
-------
m_phase: list of floats
    Values of orbital phase over time for the merger/ringdown portion.

phase_stitching
===============

``phase_stitching(final_i_index,i_phase,m_phase)``

Stitching together the inspiral and merger/ringdown portions of the phase
lists to give a combined list with the correct matching.

Parameters
----------
final_i_index: int
    The last index in the inspiral data before the switch to the merger/
    ringdown, from final_i_index_finder() in matchingfuns.
i_phase: list of floats
    Values of orbital phase at each timestep for the inspiral portion, from
    inspiral_phase_freq_integration() in inspiralfuns.
m_phase: list of floats
    Values of orbital phase over time for the merger/ringdown portion, from
    merger_phase_calculation().
    
Returns
-------
i_m_phase: list of floats
    Values of orbital phase over time for the entire duration of the
    gravitational waveform.

merger_strain_amplitude
=======================

``merger_strain_amplitude(min_switch_ind,final_i_index,alpha,i_amp,m_omega,fhat,fhatdot)``

Calculating the amplitude of strain for the merger/ringdown portion, based
on Buskirk et al. (2019) equation 16.

Parameters
----------
min_switch_ind: int
    The index in the merger/ringdown data where the switch from inspiral to
    merger/ringdown should occur, from min_switch_ind_finder() in
    matchingfuns.
final_i_index: int
    The last index in the inspiral data before the switch to the merger/
    ringdown, from final_i_index_finder() in matchingfuns.
alpha: float
    A gIRS coefficient, from gIRS_coefficients() in mergerfirstfuns.
i_amp: list of floats
    The values of the amplitude of the GW strain over time for the inspiral
    portion, from inspiral_strain_amplitude() in inspiralfuns.
m_omega: list of floats
    Values of angular frequency over time for the merger/ringdown portion,
    from merger_freq_calculation() in mergerfirstfuns.
fhat: list of floats
    Values of a sort of frequency parameter called fhat (f^) over time,
    from merger_freq_calculation() in mergerfirstfuns.
fhatdot: list of floats
    Values of the time-derivative of fhat over time, from
    fhat_differentiation() in mergerfirstfuns.
    
Returns
-------
m_amp: list of floats
    The values of the amplitude of the GW strain over time for the
    merger/ringdown portion.

amplitude_stitching
===================

``amplitude_stitching(final_i_index,i_amp,m_amp)``

Stitching together the inspiral and merger/ringdown portions of the
amplitude lists to give a combined list with the correct matching.

Parameters
----------
final_i_index: int
    The last index in the inspiral data before the switch to the merger/
    ringdown, from final_i_index_finder() in matchingfuns.
i_amp: list of floats
    The values of the amplitude of the GW strain over time for the inspiral
    portion, from inspiral_strain_amplitude() in inspiralfuns.
m_amp: list of floats
    The values of the amplitude of the GW strain over time for the
    merger/ringdown portion, from merger_strain_amplitude().
    
Returns
-------
i_m_amp: list of floats
    The values of the amplitude of the GW strain over time for the entire
    duration of the gravitational waveform.

merger_polarisations
====================

``merger_polarisations(final_i_index,m_amp,m_phase,i_Aorth)``

Calculating the values of the two polarisations of strain for the merger.

Parameters
----------
final_i_index: int
    The last index in the inspiral data before the switch to the merger/
    ringdown, from final_i_index_finder() in matchingfuns.
m_amp: list of floats
    The values of the amplitude of the GW strain over time for the
    merger/ringdown portion, from merger_strain_amplitude().
m_phase: list of floats
    Values of orbital phase over time for the merger/ringdown portion, from
    merger_phase_calculation().
i_Aorth: list of floats
    The values of the orthogonal/plus polarisation of strain over time for
    the inspiral portion, from inspiral_strain_polarisations() in
    inspiralfuns.
    
Returns
-------
[m_Aorth,m_Adiag]: list of lists of floats
    The first list is the values of the orthogonal/plus polarisation of
    strain over time, the second list is the diagonal/cross polarisation.

polarisation_stitching
======================

``polarisation_stitching(final_i_index,i_Aorth,i_Adiag,m_Aorth,m_Adiag)``

Stitching together the inspiral and merger/ringdown portions of the
polarisation lists to give combined lists with the correct matching.

Parameters
----------
final_i_index: int
    The last index in the inspiral data before the switch to the merger/
    ringdown, from final_i_index_finder() in matchingfuns.
i_Aorth: list of floats
    The values of the orthogonal/plus polarisation of strain over time for
    the inspiral portion, from inspiral_strain_polarisations() in
    inspiralfuns.
i_Adiag: list of floats
    The values of the diagonal/cross polarisation of strain over time for
    the inspiral portion, from inspiral_strain_polarisations() in
    inspiralfuns.
m_Aorth: list of floats
    The values of the orthogonal/plus polarisation of strain over time for
    the merger/ringdown portion, from merger_polarisations().
m_Adiag: list of floats
    The values of the diagonal/cross polarisation of strain over time for
    the merger/ringdown portion, from merger_polarisations().
    
Returns
-------
[i_m_Aorth,i_m_Adiag]: list of lists of floats
    The first list is the combined orthogonal/plus polarisation values, the
    second list is the combined diagonal/cross polarisation values.
*****
tests
*****

The riroriro/tests directory contains unit tests for the main Riroriro modules, utilising pytest. The tests can be performed with Tox by using the standard command ``tox -e test`` in the install directory.
*****************
snrcalculatorfuns
*****************

This is the documentation for the snrcalculatorfuns module, which consists of parts of the procedure for calculating the SNR of a gravitational waveform (from gwexporter or otherwise), collected into modular functions.

polynomial_redshift
===================

``polynomial_redshift(d)``

Polynomial approximation of calculating redshift corresponding to a given
distance.

Parameters
----------
d: float
    A luminosity distance, in Mpc.
    
Returns
-------
z: float
    The redshift corresponding to the input distance.
    
redshift_distance_adjustment
============================

``redshift_distance_adjustment(inputarray,d,z)``

Adjusts the frequencies and amplitudes in the input gravitational waveform
to account for the effects of distance/redshift.

Parameters
----------
inputarray: numpy.ndarray
    The time, frequency and amplitude data of the gravitational waveform,
    in the format used by waveform_exporter() in gwexporter.
d: float
    The luminosity distance to the merging binary, in Mpc.
z: float
    The redshift corresponding to the input distance.
    
Returns
-------
adjustedarray: numpy.ndarray
    inputarray, but with the frequency and amplitudes adjusted.
    
frequency_limits
================

``frequency_limits(inputarray)``

Calculates the upper and lower limits of the frequency of the gravitational
waveform in inputarray, which are used by amplitude_interpolation().

Parameters
----------
inputarray: numpy.ndarray
    The time, frequency and amplitude data of the gravitational waveform;
    should have been adjusted by redshift_distance_adjustment().
    
Returns
-------
(freqmax,freqmin): tuple of floats
    The upper and lower limits of the waveform signal frequency,
    respectively.
    
findchirp_fourier
=================

``findchirp_fourier(inputarray,findchirp_array,d,z)``

Approximation of a Fourier transform on the gravitational waveform data,
using the frequency spectrum output by the simpler model FINDCHIRP (Allen
et al., 2012) for calibration.
NOTE: May in the future be replaced by something fft-based.

Parameters
----------
inputarray: numpy.ndarray
    The time, frequency and amplitude data of the gravitational waveform;
    should have been adjusted by redshift_distance_adjustment().
findchirp_array: numpy.ndarray
    The array output by FINDCHIRP. The second column is frequency, the
    fourth is (Fourier-transformed) strain amplitude, the other columns
    are irrelevant. A grid of sample findchirp_arrays can be found at
    https://drive.google.com/drive/folders/12TYxYKtBL1iuFHG_ySFhS12Aqv4JHGOr
d: float
    The luminosity distance to the merging binary, in Mpc.
z: float
    The redshift corresponding to the input distance.
    
Returns
-------
fourieramp: list
    Fourier-transformed/calibrated amplitudes at each frequency value in
    inputarray.
    
amplitude_interpolation
=======================

``amplitude_interpolation(inputarray,fourieramp,noisearray,freqmax,freqmin)``

The simulated gravitational waveform data and the detector noise spectrum
are assumed to have amplitude data at different sets of frequencies, so
this function uses scipy's interp1d to calculate the waveform amplitude
values at the frequencies used by the detector data.

Parameters
----------
inputarray: numpy.ndarray
    The time, frequency and amplitude data of the gravitational waveform;
    should have been adjusted by redshift_distance_adjustment().
fourieramp: list
    Fourier-transformed/calibrated amplitudes at each frequency value in
    inputarray, from findchirp_fourier().
noisearray: numpy.ndarray
    Data on the noise spectrum of the detector; it is assumed that
    frequency values are in the first column and ASD noise levels in the
    second.
freqmax: float
    The upper limit of the waveform signal frequency, from
    frequency_limits().
freqmin: float
    The lower limit of the waveform signal frequency, from
    frequency_limits().

Returns
-------
noise_freq_amp: list
    Waveform amplitudes as in fourieramp, but over the set of frequencies
    in noisearray rather than those in inputarray.
    
individual_detector_SNR
=======================

``individual_detector_SNR(noisearray,noise_freq_amp)``

Calculates the single-detector optimal-alignment SNR by comparing the
waveform frequency spectrum and detector noise spectrum using the method of
Barrett et al. (2018).

Parameters
----------
noisearray: numpy.ndarray
    Data on the noise spectrum of the detector; it is assumed that
    frequency values are in the first column and ASD noise levels in the
    second.
noise_freq_amp: list
    Amplitudes of the simulated gravitational waveform, over the set of
    frequencies of noisearray, from amplitude_interpolation().
    
Returns
-------
ind_SNR: float
    The SNR of the simulated gravitational waveform, for the detector in
    noisearray and assuming optimal alignment.
**********
gwexporter
**********

This is the documentation for the gwexporter module, which consists of functions to output the gravitational waveforms from inspiralfuns, mergerfirstfuns, matchingfuns and mergersecondfuns into a file.

waveform_exporter
=================

``waveform_exporter(time,freq,amp,path)``

Function to export a simulated gravitational waveform into a file.

Parameters
----------
time: list of floats
    The time at each data point. For a BH-BH merger, use i_m_time from
    time_frequency_stitching in matchingfuns. For a BH-NS or NS-NS merger,
    use i_time (realtimes) from inspiral_time_conversion in inspiralfuns.
freq: list of floats
    The frequency of the GW signal at each data point. For a BH-BH merger,
    use i_m_freq from frequency_SI_units in matchingfuns. For a BH-NS or
    NS-NS merger, use i_freq (freq) from inspiral_phase_freq_integration
    in inspiralfuns.
amp: list of floats
    The amplitude of the GW strain at each data point. For a BH-BH merger,
    use i_m_amp from amplitude_stitching in mergersecondfuns. For a BH-NS
    or NS-NS merger, use i_amp from inspiral_strain_amplitude in
    inspiralfuns.
path: str
    The file path to the location/document where you want to save the
    simulated gravitational waveform data.
    
Returns
-------
An output file containing an array wherein the first column is the time,
the second is the frequency and the third is the amplitude.

waveform_arrayer
================

``waveform_arrayer(time,freq,amp)``

Function to collate important data of the simulated gravitational waveform
(for SNR calculation etc.) into a single array for ease of storage. Similar
to waveform_exporter(), but instead of outputting the area into a file,
this function outputs the data inline in Python, in a numpy.ndarray.

Parameters
----------
time: list of floats
    The time at each data point. For a BH-BH merger, use i_m_time from
    time_frequency_stitching in matchingfuns. For a BH-NS or NS-NS merger,
    use i_time (realtimes) from inspiral_time_conversion in inspiralfuns.
freq: list of floats
    The frequency of the GW signal at each data point. For a BH-BH merger,
    use i_m_freq from frequency_SI_units in matchingfuns. For a BH-NS or
    NS-NS merger, use i_freq (freq) from inspiral_phase_freq_integration
    in inspiralfuns.
amp: list of floats
    The amplitude of the GW strain at each data point. For a BH-BH merger,
    use i_m_amp from amplitude_stitching in mergersecondfuns. For a BH-NS
    or NS-NS merger, use i_amp from inspiral_strain_amplitude in
    inspiralfuns.
    
Returns
-------
exportarray: numpy.ndarray
    An array containing the important data of the simulated gravitational
    waveform: the first column is the time, the second is the frequency and
    the third is the amplitude.
*****************
detectabilityfuns
*****************

This is the documentation for the detectabilityfuns module, which consists of parts of the procedure for calculating the detectability fraction of a merger given its optimal-alignment SNR, collected into modular functions.

cdf_generator
=============

``cdf_generator(N=10**6)``

Generates the cumulative distribution function (CDF) of the projection
function Theta, for use with detectability_calculator(), based on Finn
(1996), Belczynski et al. (2013), Belczynski et al. (2014).

Parameters
----------
N: int
    The number of random samples of Theta you want to take to build the
    CDF. Default: 10**6.
    
Returns
-------
Theta_CDF: function
    The CDF of the projection function Theta.
min_CDF: float
    The lower boundary of the range over which Theta_CDF is defined.
max_CDF: float
    The upper boundary of the range over which Theta_CDF is defined.

detectability_calculator
========================

``detectability_calculator(Theta_CDF,min_CDF,max_CDF,SNR_in)``

Given the optimal-alignment SNR of a merger, this function returns the
fraction of arbitrary orientations in which the merger would be expected to
be observable (i.e. have a SNR above 8).

Parameters
----------
Theta_CDF: function
    The CDF of the projection function Theta, from cdf_generator().
min_CDF: float
    The lower boundary of the range over which Theta_CDF is defined, from
    cdf_generator().
max_CDF: float
    The upper boundary of the range over which Theta_CDF is defined, from
    cdf_generator().
SNR_in: float
    The optimal-alignment SNR of the merger in question, can be obtained
    from snrcalculatorfuns.
    
Returns
-------
det: float
    The detectability fraction of the merger.

specific_orientation_SNR
========================

``specific_orientation_SNR(theta,phi,iota,psi,SNR_in,angle_unit='rad')``

Given the optimal-alignment SNR of a merger, this function returns the SNR
that would result if the detector and binary had a specific orientation/
alignment, specified by four angles as in Finn (1996), Belczynski et al.
(2013), Belczynski et al. (2014).

Parameters
----------
theta: float
    The relative latitude, one of the angles describing the direction of
    the line of sight to the gravitational wave source relative to the axes
    of the detector’s arms (sky-location coordinates of the binary). Ranges
    from 0 to π rad (180 deg).
phi: float
    The relative longitude, one of the angles describing the direction of
    the line of sight to the gravitational wave source relative to the axes
    of the detector’s arms (sky-location coordinates of the binary). Ranges
    from 0 to 2π rad (360 deg).
iota: float
    The inclination angle of the binary. Ranges from 0 to π rad (180 deg).
psi: float
    The polarisation angle of the binary. Ranges from 0 to π (180 deg).
SNR_in: float
    The optimal-alignment SNR of the merger in question, can be obtained
    from snrcalculatorfuns.
angle_unit: str
    Specifies whether the input angles are given in 'rad' or 'deg'; the
    default is 'rad'.

Returns
-------
SNR_out: float
    The SNR of the merger in question at the specific orientation given by
    the input angles.
************
inspiralfuns
************

This is the documentation for the inspiralfuns module, which consists of parts of the procedure for simulating the inspiral portions of gravitational waves, collected into modular functions.

get_M_and_eta
=============

``get_M_and_eta(**kwargs)``

Gives total mass (M) and symmetric mass ratio (eta) from either m1 and m2
OR logMc and q; M and eta are used by many functions in the GW synthesis.

Parameters
----------
First method

m1: float
    Mass of one object in binary in solar masses.
m2: float
    Mass of other object in binary in solar masses.

Second method

logMc: float
    log10(the chirp mass of the binary in solar masses).
q: float
    The mass ratio of the objects in the binary.

Returns
-------
(M,eta): tuple of floats
    The total mass of the binary, followed by the symmetric mass ratio.

startx
======

``startx(M,flow)``

Gives starting value/lower boundary for integration of post-Newtonian
parameter, based on Buskirk et al. (2019) equation 22.

Parameters
----------
M: float
    Total mass of the binary, can be obtained from get_M_and_eta().
flow: float
    Lower cutoff frequency for the gravitational waveform, which we usually
    set to be 10 Hz.
    
Returns
-------
value: float
    The starting value for the post-Newtonian integration.

endx
====

``endx(eta,merger_type)``

Gives ending value/upper boundary for integration of post-Newtonian
parameter, based on Buskirk et al. (2019) equation 23.

Parameters
----------
eta: float
    Symmetric mass ratio of the binary, can be obtained from
    get_M_and_eta().
merger_type: string
    'BH' for a BH-BH merger, 'NS' for a BH-NS or NS-NS merger
    
Returns
-------
value: float
    The ending value for the post-Newtonian integration.

PNderiv
=======

``PNderiv(x,M,eta)``

Encodes the differential equation for the post-Newtonian parameter that is
integrated by x_integration(), based on Huerta et al. (2017) and Buskirk et
al. (2019). This function should usually not be called directly, but rather
by x_integration().

Parameters
----------
x: float
    The post-Newtonian parameter, the variable being integrated over.
M: float
    Total mass of the binary, can be obtained from get_M_and_eta().
eta: float
    Symmetric mass ratio of the binary, can be obtained from
    get_M_and_eta().

Returns
-------
Mdxdt: float
    The value of M * (dx/dt) for the input x, as given by the differential
    equation.

PN_parameter_integration
========================

``PN_parameter_integration(start,end,M,eta)``

Integrates the PNderiv() differential equation for the post-Newtonian
parameter, x.

Parameters
----------
start: float
    The starting value/lower boundary of the integration, from startx().
end: float
    The ending value/upper boundary of the integration, from endx().
M: float
    Total mass of the binary, can be obtained from get_M_and_eta().
eta: float
    Symmetric mass ratio of the binary, can be obtained from
    get_M_and_eta().
    
Returns
-------
[x,xtimes,dt]: list of lists of floats
    First list is the series of values of the post-Newtonian parameter x
    that has been integrated, second list is the time corresponding to each
    value of x (data point), third list is the timestep between each pair
    of data points.
    
inspiral_time_conversion
========================

``inspiral_time_conversion(xtimes,M)``

Converting times in geometric units from x_integration() to times in real
units.

Parameters
----------
xtimes: list of floats
    Times in geometric units of data points in the integration of the post-
    Newtonian parameter, from PN_parameter_integration().
M: float
    Total mass of the binary, can be obtained from get_M_and_eta().
    
Returns
-------
realtimes: list of floats
    xtimes, but in seconds instead of geometric units.
    
inspiral_phase_freq_integration
===============================

``inspiral_phase_freq_integration(x,dt,M)``

Integration of orbital phase and angular frequency for the inspiral, using
the post-Newtonian parameter, based on Buskirk et al. (2019) equation 7.

Parameters
----------
x: list of floats
    Values of the post-Newtonian parameter over time, from
    PN_parameter_integration().
dt: list of floats
    Timesteps in geometric units between each value of xtimes, from
    PN_parameter_integration().
M: float
    Total mass of the binary, can be obtained from get_M_and_eta().
    
Returns
-------
[i_phase,omega,freq]: list of lists of floats
    First list is the values of orbital phase at each timestep, second list
    is the angular frequency, third list is the frequency of the GW signal.
    
radius_calculation
==================

``radius_calculation(x,M,eta)``

Calculation of orbital radius (and time-derivative of radius) for the
binary for each timestep during the inspiral, based on Buskirk et al.
(2019).

Parameters
----------
x: list of floats
    Values of the post-Newtonian parameter over time, from
    PN_parameter_integration().
M: float
    Total mass of the binary, can be obtained from get_M_and_eta().
eta: float
    Symmetric mass ratio of the binary, can be obtained from
    get_M_and_eta().
    
Returns
-------
[r,rdot]: list of lists of floats
    First list is the values of the orbital radius (in geometric units) at
    each timestep, second list is the time-derivative of the radius (used
    by strain calculations).

a1_a2_calculation
=================

``a1_a2_calculation(r,rdot,omega,D,M,eta)``

Calculation of A1 and A2, two coefficients used in the calculation of
strain polarisations, based on Buskirk et al. (2019) equation 9.

Parameters
----------
r: list of floats
    Values of the orbital radius over time, from radius_calculation().
rdot: list of floats
    Values of the time-derivative of the radius, from radius_calculation().
omega: list of floats
    Values of the angular frequency over time, from
    inspiral_phase_freq_integration().
D: float
    Distance from the detector to the binary, in Mpc. IMPORTANT: if you
    want to feed the strain values into the SNR calculator, use the default
    distance of 100 Mpc here and instead set the distance when using the
    SNR functions.
M: float
    Total mass of the binary, can be obtained from get_M_and_eta().
eta: float
    Symmetric mass ratio of the binary, can be obtained from
    get_M_and_eta().
    
Returns
-------
[A1,A2]: list of lists of floats
    The first list is the values  of the A1 parameter used in strain
    calculation over time, the second list is the A2 parameter.

inspiral_strain_polarisations
=============================

``inspiral_strain_polarisations(A1,A2,i_phase)``

Calculating the values of the two polarisations of strain for the inspiral,
using the coefficients from a1_a2_calculation().

Parameters
----------
A1: list of floats
    Values of the first strain coefficient over time, from
    a1_a2_calculation().
A2: list of floats
    Values of the second strain coefficient over time, from
    a1_a2_calculation().
i_phase: list of floats
    Values of the orbital phase at each timestep, from
    inspiral_phase_freq_integration().
    
Returns
-------
[Aorth,Adiag]: list of lists of floats
    The first list is the values of the orthogonal/plus polarisation of
    strain over time, the second list is the diagonal/cross polarisation.
    
inspiral_strain_amplitude
=========================

``inspiral_strain_amplitude(Aorth,Adiag)``

Calculating the amplitude of the strain from the polarisations.

Parameters
----------
Aorth: list of floats
    The values of the orthogonal/plus polarisation of strain over time,
    from inspiral_strain_polarisations().
Adiag: list of floats
    The values of the diagonal/cross polarisation of strain over time, from
    inspiral_strain_polarisations().
    
Returns
-------
i_amp: list of floats
    The values of the amplitude of the GW strain over time (unitless).

list_size_reducer
=================

``list_size_reducer(reduction_factor,your_list)``

Optional function to reduce the size of the lists output by the inspiral
functions (not the merger lists, as those are much shorter), in order to
reduce filesize to conserve storage space.
NOTES:
The typical reduction factor we have used in our research using this code
is 100.
The inspiral lists used by the matching/merger portions are realtimes,
omega, i_phase and i_amp so if you reduce one of these you should reduce
all of them.

Parameters
----------
reduction_factor: int
    The factor you want to reduce the list length by.
your_list: list
    The list you want to reduce.
    
Returns
-------
reduced_list: list
    your_list, in reduced form.
**********************
riroriro Documentation
**********************

This is the documentation for riroriro.

Reference/API
=============

.. automodapi:: riroriro
************
matchingfuns
************

This is the documentation for the matchingfuns module, which consists of parts of the procedure for matching together the inspiral and merger/ringdown portions of the gravitational waveforms of binary black holes, collected into modular functions.

MQdiff
======

``MQdiff(i,i_time,i_omega,m_time,m_omega)``

A function that calculates a "matching quantity" that is subsequently used
to determine the best offset and switching point for matching together the
inspiral and merger/ringdown portions (the precise method was revised
several times). This function should usually not be called directly, but
rather by min_switch_ind_finder().

Parameters
----------
i: int
    An index in the range of merger/ringdown data points.
i_time: list of floats
    Real time values for the inspiral portion, from
    inspiral_time_conversion() in inspiralfuns.
i_omega: list of floats
    Values of angular frequency over time for the inspiral portion, from
    inspiral_phase_freq_integration() in inspiralfuns.
m_time: list of floats
    Real time values for the merger/ringdown portion, from
    merger_time_conversion() in mergerfirstfuns.
m_omega: list of floats
    Values of angular frequency over time for the merger/ringdown portion,
    from merger_freq_calculation() in mergerfirstfuns.
    
Returns
-------
df_diff: float (or nan)
    A value describing the difference in gradient of frequency between the
    inspiral and merger/ringdown portions for points of matching frequency.

min_switch_ind_finder
=====================

``min_switch_ind_finder(i_time,i_omega,m_time,m_omega)``

Finds the index in the merger/ringdown data where the switch from inspiral
to merger/ringdown should occur, as part of the matching process.

Parameters
----------
i_time: list of floats
    Real time values for the inspiral portion, from
    inspiral_time_conversion() in inspiralfuns.
i_omega: list of floats
    Values of angular frequency over time for the inspiral portion, from
    inspiral_phase_freq_integration() in inspiralfuns.
m_time: list of floats
    Real time values for the merger/ringdown portion, from
    merger_time_conversion() in mergerfirstfuns.
m_omega: list of floats
    Values of angular frequency over time for the merger/ringdown portion,
    from merger_freq_calculation() in mergerfirstfuns.
    
Returns
-------
min_switch_ind: int
    The index in the merger/ringdown data where the switch from inspiral to
    merger/ringdown should occur.
    
final_i_index_finder
====================

``final_i_index_finder(min_switch_ind,i_omega,m_omega)``

Finds what the last index in the inspiral data before the switch to the
merger/ringdown should be, as part of the matching process.

Parameters
----------
min_switch_ind: int
    The index in the merger/ringdown data where the switch from inspiral to
    merger/ringdown should occur, from min_switch_ind_finder().
i_omega: list of floats
    Values of angular frequency over time for the inspiral portion, from
    inspiral_phase_freq_integration() in inspiralfuns.
m_omega: list of floats
    Values of angular frequency over time for the merger/ringdown portion,
    from merger_freq_calculation() in mergerfirstfuns.
    
Returns
-------
final_i_index: int
    The last index in the inspiral data before the switch to the merger/
    ringdown.

time_offset_finder
==================

``time_offset_finder(min_switch_ind,final_i_index,i_time,m_time)``

Calculates what the offset between the time values of the inspiral and
merger/ringdown portions should be to match them together.

Parameters
----------
min_switch_ind: int
    The index in the merger/ringdown data where the switch from inspiral to
    merger/ringdown should occur, from min_switch_ind_finder().
final_i_index: int
    The last index in the inspiral data before the switch to the merger/
    ringdown, from final_i_index_finder().
i_time: list of floats
    Real time values for the inspiral portion, from
    inspiral_time_conversion() in inspiralfuns.
m_time: list of floats
    Real time values for the merger/ringdown portion, from
    merger_time_conversion() in mergerfirstfuns.
    
Returns
-------
time_offset: float
    The offset between the time values of the inspiral and merger/ringdown
    portions.
    
time_frequency_stitching
========================

``time_frequency_stitching(min_switch_ind,final_i_index,time_offset,i_time,i_omega,m_time,m_omega)``

Stitches together the inspiral and merger/ringdown portions of the time and
angular frequency lists to give combined lists for these with the correct
matching.

Parameters
----------
min_switch_ind: int
    The index in the merger/ringdown data where the switch from inspiral to
    merger/ringdown should occur, from min_switch_ind_finder().
final_i_index: int
    The last index in the inspiral data before the switch to the merger/
    ringdown, from final_i_index_finder().
time_offset: float
    The offset between the time values of the inspiral and merger/ringdown
    portions, from time_offset_finder().
i_time: list of floats
    Real time values for the inspiral portion, from
    inspiral_time_conversion() in inspiralfuns.
i_omega: list of floats
    Values of angular frequency over time for the inspiral portion, from
    inspiral_phase_freq_integration() in inspiralfuns.
m_time: list of floats
    Real time values for the merger/ringdown portion, from
    merger_time_conversion() in mergerfirstfuns.
m_omega: list of floats
    Values of angular frequency over time for the merger/ringdown portion,
    from merger_freq_calculation() in mergerfirstfuns.
    
Returns
-------
[i_m_time,i_m_omega]: list of lists of floats
    The first list is the combined time values, the second list is the
    combined angular frequency values.
    
frequency_SI_units
==================

``frequency_SI_units(i_m_omega,M)``

The angular frequency in geometric units translated to ordinary/temporal
frequency in SI units (Hz). Useful for plotting and also required for the
SNR calculator.

Parameters
----------
i_m_omega: list of floats
    Values of angular frequency over time for the entire duration of the
    gravitational waveform, from time_frequency_stitching().
M: float
    Total mass of the binary, can be obtained from get_M_and_eta() in
    inspiralfuns.

Returns
-------
i_m_freq: list of floats
    Values of frequency in Hz for the entire duration of the gravitational
    waveform.
***************
horizondistfuns
***************

This is the documentation for the horizondistfuns module, which consists of parts of the procedure for calculating the horizon distance of a gravitational waveform (from gwexporter or otherwise), collected into modular functions.

compact_SNR_calculation
=======================

``compact_SNR_calculation(inputarray,findchirp_array,noisearray_list,method,d)``

Runs through all of the functions of snrcalculatorfuns to obtain a SNR from
an individual detector. This function is mainly included not to be called
directly, but rather by horizon_distance_calculation().

Parameters
----------
inputarray: numpy.ndarray
    The time, frequency and amplitude data of the gravitational waveform,
    in the format used by waveform_exporter() in gwexporter.
findchirp_array: numpy.ndarray
    The array output by FINDCHIRP. The second column is frequency, the
    fourth is (Fourier-transformed) strain amplitude, the other columns
    are irrelevant. A grid of sample findchirp_arrays can be found at
    https://drive.google.com/drive/folders/12TYxYKtBL1iuFHG_ySFhS12Aqv4JHGOr
noisearray_list: list of numpy.ndarrays
    Each item in this list should be an array describing the noise spectrum
    of a detector; in each noise spectrum, it is assumed that frequency
    values are in the first column and ASD noise levels in the second.
method: str
    If 'quad', returns the quadrature SNR across the detectors in
    noisearray_list. If 'mean', returns the mean of the SNRs with each
    individual detector (simulating one random detector in operation). If
    only one detector is included in noisearray_list, these methods are
    equivalent.
d: float
    The luminosity distance to the merging binary, in Mpc.
    
Returns
-------
final_SNR: float
    The SNR of the simulated gravitational waveform, for the detectors in
    noisearray and assuming optimal alignment.
    
horizon_distance_calculation
============================

``horizon_distance_calculation(inputarray,findchirp_array,noisearray_list,method)``

Calculates the horizon distance (maximum distance at which something can
be observed) given optimal alignment for a given merger.

Parameters
----------
inputarray: numpy.ndarray
    The time, frequency and amplitude data of the gravitational waveform,
    in the format used by waveform_exporter() in gwexporter.
findchirp_array: numpy.ndarray
    The array output by FINDCHIRP. The second column is frequency, the
    fourth is (Fourier-transformed) strain amplitude, the other columns
    are irrelevant. A grid of sample findchirp_arrays can be found at
    https://drive.google.com/drive/folders/12TYxYKtBL1iuFHG_ySFhS12Aqv4JHGOr
noisearray_list: list of numpy.ndarrays
    Each item in this list should be an array describing the noise spectrum
    of a detector; in each noise spectrum, it is assumed that frequency
    values are in the first column and ASD noise levels in the second.
method: str
    If 'quad', uses the quadrature SNR across the detectors in
    noisearray_list. If 'mean', uses the mean of the SNRs with each
    individual detector (simulating one random detector in operation). If
    only one detector is included in noisearray_list, these methods are
    equivalent.
    
Returns
-------
horizon_dist: float
    The horizon distance of the given merger, for the given detector(s).
