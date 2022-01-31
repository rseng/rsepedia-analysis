# Contributing to WODEN

If you have a cool feature request, a burning passion to add some functionality
to `WODEN`, or have found an annoying bug, please reach out. I welcome the help!
Here is a quick guide on how to go about it.

## I have an idea for a feature
Please raise an issue suggesting the new feature, and chuck the `enhancement`
label on the issue. Try and describe the feature as completely as possible and
why it would be useful for you. We can then start a dialogue and work out the
best way forward / where in the priority list of new features to implement this
it might lie.

## I found a bug
`#sadface`, these things happen. Please raise an issue and add the `bug` label on the issue (please check the issue doesn't already exist). Include a description of what's going wrong, including an example command that causes the error for you.

## I have an idea/found a bug AND I want to fix it myself
Awesome, go for it! First off, make sure you've submitted an issue so we know
the bug/feature isn't being fixed or implemented in another branch somewhere. After
a little discussion to clarify the goals, make your own branch or fork, and
then start working on that branch/fork. For any piece of work, I'd ask that:

 - If fixing code, ensure the fix still passes any existing unit test as detailed in [the online documentation here](https://woden.readthedocs.io/en/joss_review/testing/cmake_testing.html#what-do-the-tests-actually-do). If the fix changes the expected outcomes, discuss it on the issue, and if appropriate, you can update the existing unit tests.
 - If making new code, write new unit tests that cover as much of the new code as possible.
 - Document whatever you've done so we all know what good stuff you've added in.

If you've not written test code before, or made documentation via `sphinx`, just let me know on whatever git issue you've started (maybe use the "help wanted" label) and we can work together on it.

Once all that's done, submit a pull request and we'll get it in the main branch.

## Response time
At the moment, it's just me (Jack Line) writing the code, so please
be patient with me and I'll get back to you as fast as my schedule allows.
# WODEN

[![status](https://joss.theoj.org/papers/bbc90ec4cd925ade93ed0781e571d247/status.svg)](https://joss.theoj.org/papers/bbc90ec4cd925ade93ed0781e571d247)

[![Documentation Status](https://readthedocs.org/projects/woden/badge/?version=latest)](https://woden.readthedocs.io/en/latest/?badge=latest) [![codecov](https://codecov.io/gh/JLBLine/WODEN/branch/joss_review/graph/badge.svg?token=Q3JFCI5GOC)](https://codecov.io/gh/JLBLine/WODEN) _*note code coverage only applies to `python3` and `C` code, `CUDA` code is not currently supported by code coverage software_

> The `WODEN` documentation lives [here on readthedocs](https://woden.readthedocs.io/en/latest/). If your internet has broken and you have already installed `WODEN`, you can build a local copy by navigating into `WODEN/docs/sphinx` and running `make html` (you'll need to have `doxygen` installed).

`WODEN` is C / CUDA code designed to be able to simulate low-frequency radio interferometric data. It is written to be simplistic and *fast* to allow all-sky simulations. Although `WODEN` was primarily written to simulate Murchinson Widefield Array ([MWA, Tingay et al. 2013](https://doi.org/10.1017/pasa.2012.007)) visibilities, it is becoming less instrument-specific as time goes on. `WODEN` outputs `uvfits` files.

The unique part of `WODEN` is that it can simulate shapelet model sources (along with point and Gaussian) that are compatible with the `RTS` ([Mitchell et al. 2008](https://ieeexplore.ieee.org/document/4703504?arnumber=4703504 "IEEExplorer")). These models are generated with SHApelet Modelling For Interferometers ([SHAMFI](https://github.com/JLBLine/SHAMFI)), specified with the `--woden_srclist` SHAMFI option. It also includes a script to convert a multi-scale CLEAN component list out of [WSClean](https://sourceforge.net/projects/wsclean/) into a `WODEN`-style srclist (when running `WSClean` use the `-save-source-list` option). `WODEN` can also produce visibilities that can be fed directly into the `RTS` to allow testing of calibration and modelling methodologies.

If you have feature requests or want to contribute to `WODEN`, have a read of the
[guide to contributing](CONTRIBUTION_GUIDE.md) to get started. I welcome the feedback and/or help!

Jack Line \
January 2022


## 1. Installation
Read the comprehensive [installation guide on readthedocs](https://woden.readthedocs.io/en/latest/installation/installation.html#dependencies). In short, you will need the dependencies:

- CMake - https://cmake.org version >= 3.10
- NVIDIA CUDA - https://developer.nvidia.com/cuda-downloads
- json-c - https://github.com/json-c/json-c
- ERFA - https://github.com/liberfa/erfa/releases
- HDF5 - https://www.hdfgroup.org/downloads/hdf5/
- PAL - https://github.com/Starlink/pal/releases
- python >= 3.6

Once you have those, installation is done via `CMake`. Ideally, this will be enough:
```bash
$ cd WODEN
$ mkdir build && cd build
$ cmake ..
$ make -j 4
$ sudo make install #(this is optional)
```
with a couple of post-compilation environment variables needed (if you don't want to `make install`). Checkout the [installation guide on readthedocs](https://woden.readthedocs.io/en/latest/installation/installation.html#dependencies) for full details.

## 2. Testing
There are two routes to test WODEN:
- via the unit/integration tests run by `ctest`, which will test your dependencies (see [cmake testing on readthedocs](https://woden.readthedocs.io/en/latest/testing/cmake_testing.html))
- via simple example scripts, which will test your installation (see [testing via scripts on readthedocs](https://woden.readthedocs.io/en/latest/testing/script_testing.html))

## 3. Usage

### 3.1 Python wrapper

The recommended way to run `WODEN` is via the wrapper script `run_woden.py`, which will launch the `woden` executable. All [ `run_woden.py` arguments are explained here on readthedocs](https://woden.readthedocs.io/en/latest/API_reference/python_code/run_woden.html) (or you can type `run_woden.py --help`).

As `WODEN` is written primarily for the MWA, the simplest way to feed in observational properties is via a `metafits` file. You can [obtain metafits files here](https://asvo.mwatelescope.org/). A minimalistic example command using a metafits file looks like

```bash
run_woden.py \
    --ra0=50.67 --dec0=-37.2 \
    --cat_filename=srclist_msclean_fornaxA_phase1+2.txt \
    --metafits_filename=1202815152_metafits_ppds.fits \
    --primary_beam=MWA_FEE
```
where the array layout, observational settings, frequency and time resolution are all read in from the metafits file. All you have to specify is a phase centre (`ra0`, `dec0`), a sky model (`cat_filename`), and a primary beam (in this case the MWA Fully Embedded Element beam - the delays are also read from the metafits). For full control, you can also specify all observational settings explicitly:

```bash
run_woden.py \
    --ra0=50.67 --dec0=-37.2 \
    --cat_filename=srclist_msclean_fornaxA_phase1+2.txt \
    --primary_beam=MWA_FEE \
    --MWA_FEE_delays=[6,4,2,0,8,6,4,2,10,8,6,4,12,10,8,6] \
    --lowest_channel_freq=169.6e+6 \
    --freq_res=10e+3 \
    --num_time_steps=240 \
    --time_res=0.5 \
    --date=2018-02-28T08:47:06 \
    --array_layout=MWA_phase2_extended.txt
```

`WODEN` can read in an array layout as specified in local east, north, and height, and can simulate for any location on Earth (defaults to the MWA but you can specify Earth location with `--longitude` / `--latitude`).

### 3.2 Direct command line

Alternatively, one can run WODEN directly via the command line, using a `.json` file such as
```sh
woden input_file.json
```
where an `input_file.json` reads as
```json
{
  "ra0": 50.6700000000,
  "dec0": -37.2000000000,
  "num_freqs": 128,
  "num_time_steps": 240,
  "cat_filename": "srclist_msclean_fornaxA_phase1+2.txt",
  "time_res": 0.50000,
  "frequency_resolution": 10000.000,
  "chunking_size": 0,
  "jd_date": 2458165.9714583335444331,
  "LST": 72.79744734,
  "array_layout": "WODEN_array_layout.txt",
  "lowest_channel_freq": 1.6959500000e+08,
  "latitude": -26.70331944,
  "coarse_band_width": 1.2800000000e+06,
  "band_nums": [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
}
```
with the input values matching the options listed in the table above. This will only create binary output files however, and not `uvfits` files.

### 3.3 Using multi-scale CLEAN components from WSCLEAN
To generate a `WODEN`-style source catalogue from a `WSClean` multi-scale CLEAN, simply use the script `convert_WSClean_list_to_WODEN.py`, which will be linked to the build directory upon compilation. An exmaple command is
```sh
convert_WSClean_list_to_WODEN.py \
 --file=msclean_output_from_WSClean-sources.txt \
 --outname=srclist-woden_wsclean-model.txt
```

## 4. WODEN sky model
The WODEN sky model uses point sources, elliptical Gaussians, and shapelet models. All sources are currently given a power-law spectral behaviour. A full [breakdown of the sky model lives here on readthedocs](https://woden.readthedocs.io/en/latest/operating_principles/skymodel.html). As an idea, here is a simple sky model of thee point sources:

```
SOURCE multi_point P 3 G 0 S 0 0
COMPONENT POINT 4.0 -27.0
LINEAR 1.8e+08 10.0 0 0 0 -0.4
ENDCOMPONENT
COMPONENT POINT 3.0 -37.0
LINEAR 1.3e+08 1.0.0 0 0 0 -0.786
ENDCOMPONENT
COMPONENT POINT 5.0 -47.0
LINEAR 3.9e+08 0.04 0 0 0 .02
ENDCOMPONENT
ENDSOURCE
```
---
title: '`WODEN`: A CUDA-enabled package to simulate low-frequency radio interferometric data'
tags:
  - Python
  - C
  - CUDA
  - radio astronomy
  - interferometers
authors:
  - name: Jack L. B. Line
    orcid: 0000-0002-9130-5920
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: International Centre for Radio Astronomy Research, Curtin Institute of Radio Astronomy, Perth, WA 6102, Australia
   index: 1
 - name: ARC Centre of Excellence for All Sky Astrophysics in 3 Dimensions (ASTRO 3D)
   index: 2
date: XXX
bibliography: paper.bib

---
# Summary

`WODEN` is designed to simulate the response of a class of telescope known as an interferometer, producing output "visibilities" for a given astrophysical sky model. Simulated observations allow us to test other software packages that are designed to calibrate and analyse real interferometric data, including verifying expected behaviour with known inputs, and testing new sky modelling techniques. The `WODEN` sky model can be specified in Dirac-delta like functions on the sky (known in the field as "point sources"), elliptical Gaussian models, or built out of "shapelet" basis functions, allowing complicated morphologies to be created. Users are able to input a bespoke layout for the interferometer, vary a number of observational parameters including time of day, length of observation and frequency coverage, and select from a number of predefined primary beams which encode the response of the receiving elements of an interferometer. This allows simulations of a number of telescopes to be undertaken. `WODEN` works with input Stokes $I,Q,U,V$ polarisations as a sky model, simulating telescopes with dual linear polarisations, and outputting linear Stokes polarisations.

The core functionality of `WODEN` is written in CUDA as interferometric simulations are computationally intensive but embarrassingly parallel. The performance of CUDA allows for large-scale simulations to be run including emission from all directions in the sky. This is paramount for interferometers with a wide field of view such as the Murchison Widefield Array [MWA, @Tingay2013]. A Python wrapper is used to take advantage of community packages such as [astropy](https://www.astropy.org/) [@astropy2013; @astropy2018] and [pyerfa](https://pypi.org/project/pyerfa/) [@pyerfa] and to present a user-friendly interface to `WODEN`. Those simulating MWA observations can use the MWA `metafits` file to quickly feed in observational parameters to `WODEN` to match real data.

`WODEN` can be run to two levels of precision: a `woden_float` precision (which uses a mix of 32 and 64 bit floating precision), and a `woden_double` (which uses nearly entirely 64 bit precision). In the section titled "Estimation of accuracy and computational speed" below, `WODEN` is shown to produce visibilities to within 0.2% of the expected values when running in `woden_float` mode, and 0.000002% in `woden_double` mode, for baselines of length $\le 10$km.

# Underlying methodolgy

An interferometer creates visibilities $V$ by cross-correlating signals detected between pairs of antennas or dishes (baselines), described by coordinates $u,v,w$. Each visibility is sensitive to the entire sky, directions of which we describe by the direction cosines $l,m,n$. Ignoring the antenna response, the full integral over the sky can be discretised as

\begin{equation} \label{eq:RIME}
V_s(u_i,v_i,w_i) = \\ \sum_j \mathcal{S}_s(l_j,m_j) \exp[-2\pi i(u_il_j + v_im_j + w_i(n_j-1))],
\end{equation}

where $u_i,v_i,w_i$ are the visibility coordinates of the $i^{\mathrm{th}}$ baseline, $l_j$, $m_j$, $n_j$ is the sky position of the $j^{\mathrm{th}}$ component in the sky model, and $\mathcal{S}(l_j,m_j)$ is the flux density of that component in a given Stokes polarisation $s$. `WODEN` simulates dual-linear polarisation antennas, with each
antenna/station having its own primary beam shape. I can define the response of a dual polarisation antenna to direction $l,m$ as

$$
\mathbf{J}(l,m) =
\begin{bmatrix}
g_{\mathrm{ns}}(l,m) & D_{\mathrm{ns}}(l,m) \\
D_{\mathrm{ew}}(l,m) & g_{\mathrm{ew}}(l,m)
\end{bmatrix},
$$

where $g$ are gain terms, $D$ are leakage terms, and $\mathrm{ns}$ refers to north-south and $\mathrm{ew}$ east-west aligned antennas. When calculating the cross-correlation responses from antennas 1 and 2 towards direction $l,m$ to produce linear polarisation visibilities, these gains and leakages interact with the four Stokes polarisations $I,Q,U,V$ as

\begin{equation}\label{eq:RIME_full}
\begin{bmatrix}
V_{12\,XX}(l,m) \\
V_{12\,XY}(l,m) \\
V_{12\,YX}(l,m) \\
V_{12\,YY}(l,m)
\end{bmatrix} =
\mathbf{J}_1(l,m) \otimes \mathbf{J}_2^*(l,m)
\begin{bmatrix}
1 & 1 & 0 & 0 \\
0 & 0 & 1 & i \\
0 & 0 & 1 & -i \\
1 & -1 & 0 & 0
\end{bmatrix}
\begin{bmatrix}
V_{12\,I}(l,m) \\
V_{12\,Q}(l,m) \\
V_{12\,U}(l,m) \\
V_{12\,V}(l,m)
\end{bmatrix}
\end{equation}



where $*$ denotes a complex conjugate, and $\otimes$ an outer product (the result
of this outer product is written explicitly in the `WODEN` documentation [here](https://woden.readthedocs.io/en/joss_review/operating_principles/visibility_calcs.html)). For each baseline, frequency, and time step, `WODEN` calculates all four linear Stokes polarisations ($V_{XX}, V_{XY}, V_{YX}, V_{YY}$) as defined above for all $l_j,m_j$ in the sky model, and then sums over $j$, to produce four full-sky linear Stokes polarisation visibilities per baseline/frequency/time.

For a telescope like the MWA, the primary beam $\mathbf{J}(l,m)$ is a complicated pattern on the sky, which is sensitive to emission from directly overhead to all the way down to the horizon. To truly capture the effects of astrophysical foregrounds we therefore have to simulate the entire sky. The MWA Fully Embedded Element [FEE, @Sokolowski2017] model is currently the most accurate representation of the MWA primary beam, and is incorporated into `WODEN`.

As the sky model of `WODEN` is a list of Right Ascension and Declinations with associated flux densities, the user has full control over the projection of the sky into visibilities. To simulate discrete foregrounds, one can simply input any sky catalogue specified in RA/Dec. For diffuse sky models, one could for example input a list of point source/elliptical Gaussians following the HEALPix projection [@HEALPix2005], or employ a TAN or SIN FITS [@FITS2002] projection. ``WODEN`` will simply calculate the measurement equation for all directions in the sky model.

# Statement of need

Under this discrete sky formalism, upwards of $j\ge25\times10^6$ components can be required to achieve the angular resolution required. Furthermore, $u,v,w$ are time and frequency dependent, so to sample in frequency of order 500 times and 100 samples in time, there are of order $10^{12}$ visibility calculations to make. This makes CUDA acceleration paramount.

Alternative approaches to interferometric simulations exist, such as [pyuvsim](https://github.com/RadioAstronomySoftwareGroup/pyuvsim) [@Lanman2019], which sacrifices speed for excellent precision, and [RIMEz](https://github.com/upenneor/rimez), which decomposes the sky into spherical harmonics rather than discrete points. `WODEN` was designed with the Australian MWA Epoch of Reionisation (EoR) processing pipeline in mind, which uses a calibration and foreground removal software called the `RTS` [@Mitchell2008] in search of signals from the very first stars [see @Yoshiura2021 for a recent use of this pipeline]. The `RTS` creates a sky model using the same formalism above, however the code is not optimised enough to handle the volume of sources to simulate the entire sky. To test the `RTS` method of sky generation, we therefore needed a fast and discretised method. Another excellent CUDA accelerated simulation package, [OSKAR](https://github.com/OxfordSKA/OSKAR) [@OSKAR], addresses these two points. However, the `RTS` also generates parts of the sky model via shapelets [see @Line2020 for an overview], which `OSKAR` cannot. Furthermore, in real data, the precession/nutation of the Earth's rotational axis causes sources to move from the sky coordinates as specified in the RA, DEC J2000 coordinate system. The `RTS` is designed to undo this precession/nutation, and so a simulation fed into the `RTS` should *contain* precession. `WODEN` adds in this precession using the same method as the `RTS` to be consistent. This unique combination of CUDA, shapelet foregrounds, the MWA FEE primary beam, along with source precession, created the need for `WODEN`. These effects should not preclude other calibration packages from using `WODEN` outputs however, meaning `WODEN` is not limited to feeding data into the `RTS` alone.

# Estimation of accuracy and computational speed

The goal of this section is to test the accuracy of the functionality of `WODEN`, including reading of inputs, the array coordinate calculations, the precession/nutation correction, $l,m,n$ and $u,v,w$ calculations, flux density frequency extrapolation via spectral index, calculation of Equation \ref{eq:RIME_full}, and writing out of the data to `uvfits` files.

To test the absolute accuracy of `WODEN`, we first need a set of input parameters that have an analytically predictable outcome. If we ignore the beam response and polarisation, set the flux density of a source to one, and consider a single baseline and sky direction, the measurement equation (Equation \ref{eq:RIME}) becomes[^1]

[^1]: Note there is no negative at the front inside the exponential for $V(u,v,w)$. After numerous comparisons to other simulation packages, and imaging to check the input source positions match, I find dropping the negative gives the correct outputs.

\begin{equation} \label{eq:RIME_simple}
V(u,v,w) = \exp[2\pi i(ul + vm + w(n-1))].
\end{equation}

We can use Euler's formula to split $V$ into real and imaginary components. If
I label the phase for a particular source and baseline as

$$
  \phi = 2\pi \left( ul + vm + w(n-1)\right)
$$

then the real and imaginary parts of the visibility $V_{re}$, $V_{im}$ are

$$
  V_{re} = \cos(\phi), \quad V_{im} = \sin(\phi).
$$

If we can therefore set $\phi$ to a number of values which produce known sine and cosine outputs, by selecting specific combinations of $u,v,w$ and $l,m,n$, we can simulate visibilities with predictable outputs. First of all, consider the simplified case $\phi_{\mathrm{simple}}$ when $u,v,w = 1,1,1$. In that case,

$$
  \frac{\phi_{\mathrm{simple}}}{2\pi} = l + m + (n-1).
$$

If we further set $l = m$, we end up with

$$
\begin{aligned}
  \frac{\phi_{\mathrm{simple}}}{2\pi} = 2l + (n-1), \\
  l = \sqrt{\left( \frac{1 - n^2}{2} \right)}
\end{aligned}
$$

It can be shown (via [Wolfram Alpha](https://www.wolframalpha.com/widgets/view.jsp?id=c07cc70f1e81887dfd0971d3fe17cfcd)) that a solution for $n$ is

$$
  n = \frac{\sqrt{2}\sqrt{-\phi_{\mathrm{simple}}^2 - 4\pi\phi_{\mathrm{simple}} + 8\pi^2} + \phi_{\mathrm{simple}} + 2\pi}{6\pi}
$$

which we can then use to calculate values for $l,m$ through

$$
  l = m = \sqrt{\frac{1 - n^2}{2}}.
$$

Practically then, if we input the following combinations of $l,m,n$ into Equation \ref{eq:RIME_simple} our output visibilities should exactly match the $\cos(\phi)$, $\sin(\phi)$ values.

\begin{table}[h]
\begin{center}
\begin{tabular}{ c c c c c }
\hline
$\phi_{\mathrm{simple}}$ & $l,m$ & $n$ & $\cos(\phi)$ & $\sin(\phi)$ \\
\hline
\hline
$0$ & 0.0 & 1.0 & $1.0$ & $0$ \\
$\pi/6$ & 0.0425737516338956 & 0.9981858300655398 & $\sqrt{3}/2$ & $0.5$ \\
$\pi/4$ & 0.0645903244635131 & 0.9958193510729726 & $\sqrt{2}/2$ & $\sqrt{2}/2$ \\
$\pi/3$ & 0.0871449863555500 & 0.9923766939555675 & $0.5$ & $\sqrt{3}/2$ \\
$\pi/2$ & 0.1340695840364469 & 0.9818608319271057 & $0.0$ & $1.0$ \\
$2\pi/3$ & 0.1838657911209207 & 0.9656017510914922 & $-0.5$ & $\sqrt{3}/2$ \\
$3\pi/4$ & 0.2100755148372292 & 0.9548489703255412 & $-\sqrt{2}/2$ & $\sqrt{2}/2$ \\
$5\pi/6$ & 0.2373397982598921 & 0.9419870701468823 & $-\sqrt{3}/2$ & $0.5$ \\
$\pi$ & 0.2958758547680685 & 0.9082482904638630 & $-1.0$ & $0.0$ \\
$7\pi/6$ & 0.3622725654470420 & 0.8587882024392495 & $-\sqrt{3}/2$ & $-0.5$ \\
$5\pi/4$ & 0.4003681253515569 & 0.8242637492968862 & $-\sqrt{2}/2$ & $-\sqrt{2}/2$ \\
\hline
\end{tabular}
\caption{$l,m,n$ combinations used in accuracy test}
\label{tab:lmn_combos}
\end{center}
\end{table}

To test for a range of baseline lengths, we can make a simplification where we set all baseline coordinates to be equal, i.e. $u = v = w = b$ where $b$ is some length in units of wavelength. In this form, the phase including the baseline length $\phi_{b}$ is

$$
  \phi_{b} = 2\pi b\left( l + m + n - 1 \right) = b\phi_{\mathrm{simple}}.
$$

As sine/cosine are periodic functions, the following is true:

$$
  \phi_{\mathrm{simple}} = \phi_{\mathrm{simple}} + 2\pi \mathrm{n}
$$

where $\mathrm{n}$ is some integer. This means for a given $\phi_{\mathrm{simple}}$, we can find an appropriate $b$ that should still result in the expected sine and cosine outputs by setting

\begin{gather*}
  b\phi_{\mathrm{simple}} = \phi_{\mathrm{simple}} + 2\pi \mathrm{n}, \\
  b = \frac{\phi_{\mathrm{simple}} + 2\pi \mathrm{n}}{\phi_{\mathrm{simple}}}
\end{gather*}

for a range of $\mathrm{n}$ values. The values of $\mathrm{n}$ and the resultant size of b that I use in testing are shown in Table \ref{tab:b_values}.

\begin{table}[h]
\begin{center}
\begin{tabular}{c c c c c c }
\hline
$\phi_{\mathrm{simple}}$ & $b(\mathrm{n=1})$ & $b(\mathrm{n=10})$ & $b(\mathrm{n=100})$ & $b(\mathrm{n=1000})$ & $b(\mathrm{n=10000})$ \\
\hline
\hline
$0$ & 6.3 & 62.8 & 628.3 & 6283.2 & 62831.9 \\
$\pi/6$ & 13.0 & 121.0 & 1201.0 & 12001.0 & 120001.0 \\
$\pi/4$ & 9.0 & 81.0 & 801.0 & 8001.0 & 80001.0 \\
$\pi/3$ & 7.0 & 61.0 & 601.0 & 6001.0 & 60001.0 \\
$\pi/2$ & 5.0 & 41.0 & 401.0 & 4001.0 & 40001.0 \\
$2\pi/3$ & 4.0 & 31.0 & 301.0 & 3001.0 & 30001.0 \\
$3\pi/4$ & 3.7 & 27.7 & 267.7 & 2667.7 & 26667.7 \\
$5\pi/6$ & 3.4 & 25.0 & 241.0 & 2401.0 & 24001.0 \\
$\pi$ & 3.0 & 21.0 & 201.0 & 2001.0 & 20001.0 \\
$7\pi/6$ & 2.7 & 18.1 & 172.4 & 1715.3 & 17143.9 \\
$5\pi/4$ & 2.6 & 17.0 & 161.0 & 1601.0 & 16001.0 \\
\hline
\end{tabular}
\caption{Range of baseline lengths used in conjunction with the $l,m,n$ coordinates in Table~\ref{tab:lmn_combos}.}
\label{tab:b_values}
\end{center}
\end{table}

`WODEN` reads in an input array layout specified in local east, north, height $E,N,H$ coordinates. It then converts those into local $X,Y,Z$ coordinates via the equations

\begin{gather}
  X = -\sin(\phi_{\mathrm{lat}})N + \cos(\phi_{\mathrm{lat}})H \label{eq:xyz_calc1} \\
  Y = E \label{eq:xyz_calc2} \\
  Z = \cos(\phi_{\mathrm{lat}})N + \sin(\phi_{\mathrm{lat}})H \label{eq:xyz_calc3}
\end{gather}

where $\phi_{\mathrm{lat}}$ is the latitude of the array. $X,Y,Z$ are used to calculate the $u,v,w$ coodinates (c.f. Chapter 4 in @TMSthird). If we place our interferometer at a $\phi_{\mathrm{lat}} = 0.0^\circ$ and set the local sidereal time (LST) to zero, the calculation of $u,v,w$ becomes

\begin{equation}\label{eq:uvw_simple}
u = E; \, v = N; \, w = H;
\end{equation}

allowing us to set $E, N, H = b$ for our values on $b$ in Table \ref{tab:b_values}. Furthermore, we can convert our $l,m$ values from Table \ref{tab:lmn_combos} into RA,Dec ($\alpha, \delta$) via:

\begin{gather}
  \delta = \arcsin(l)  \label{eq:dec_simple} \\
  \alpha = \arcsin \left( \frac{l}{\cos(\arcsin(l))} \right) \label{eq:ra_simple}
\end{gather}

Following the `RTS`, `WODEN` first of all calculates *X,Y,Z* using the array latitude at the time of the observation. It then uses the `PAL` [@PAL2013] [palPrenut](https://github.com/Starlink/pal/blob/master/palPrenut.c) function to generate a rotation matrix to rotate the local *X,Y,Z* coordinates back to the J2000 epoch, as well as the LST and latitude of the array. This accounts for the precession/nutation of the Earth with respect to the J2000 RA/Dec coordinates that the sky model is specified in. To manifest the outcomes of Equations \ref{eq:uvw_simple}, \ref{eq:dec_simple}, and \ref{eq:ra_simple}, we have to apply the opposite rotation about $\phi_{\mathrm{lat}}$ as defined by Equations \ref{eq:xyz_calc1}, \ref{eq:xyz_calc2}, and \ref{eq:xyz_calc3}, as well as the rotations applied via [palPrenut](https://github.com/Starlink/pal/blob/master/palPrenut.c) to account for precession/nutation, to our input $E,N,H$ coordinates.

Figure \ref{fig:WODEN_accuracy} shows the result of running multiple simulations, each with: an array layout with a single baseline; a single time step and frequency channel; a single point source sky model; a primary beam model with gains of one and zero leakage. All possible combinations of $l,m,n$ and $b$ as listed in Tables \ref{tab:lmn_combos} and \ref{tab:b_values} are run. Each simulation is run with the parameters specified in Table \ref{tab:acc_sim_settings}.

\renewcommand{\arraystretch}{1.4}
\begin{table}[h]
\begin{center}
\begin{tabular}{p{0.25\linewidth} p{0.15\linewidth} p{0.5\linewidth}}
\hline
Parameter & Value & Manifestation in simulation \\
\hline
\hline
Date (UTC) & 2020-01-01 12:00:00.0 & \texttt{WODEN} must correct for precession and nutation \\
Latitude (deg) & 0.1095074 & After precess/nut correction, latitude is 0.0$^\circ$ \\
Longitude (deg) & 79.6423588 & After precess/nut correction, LST is 0.0$^\circ$ \\
Frequency (MHz) & 299.792458 & Means $\lambda = 1$, so wavelength scaled $u,v,w = E,N,H$ \\
Reference frequency for sky model (MHz) & 150 & \texttt{WODEN} has to extrapolate the flux density \\
Spectral Index & -0.8 & Needed to extrapolate flux density \\
Reference Stokes I flux density (Jy) & 1.7401375 & Should be extrapolated to a flux of 1.0 at the simulation frequency \\
\hline
\hline
\end{tabular}
\caption{Common settings for the simulations run to produce the results in Figure \ref{fig:WODEN_accuracy}}
\label{tab:acc_sim_settings}
\end{center}
\end{table}

![The absolute fractional difference (in percent) of visibilities calculated by `WODEN`, compared to their expected values, with the real component shown on the left, and the imaginary shown on the right. The green triangles show an older version of `WODEN` which used only 32 bit precision; the orange square show the v1.1 `woden_float` version which uses a mixture of 32 and 64 bit precision; the blue crosses show the `woden_double` mode which uses nearly entirely 64 bit precision.\label{fig:WODEN_accuracy}](quantify_woden_accuracy.png)

All array layouts, sky models, and simulations are run by `WODEN/test_installation/absolute_accuracy/run_the_absolute_accuracy_test.sh`, which can be run as part of a test suite bundled with `WODEN`. This script reads the values out of the output `uvfits` files, and produces the plot in Figure \ref{fig:WODEN_accuracy}.

[Version 1.0](https://github.com/JLBLine/WODEN/releases/tag/v1.0.0) of `WODEN` was fully 32 bit, which produced the green triangles in Figure \ref{fig:WODEN_accuracy}, with longer baselines consistently a few percent off expectations. A two time processing slow down by moving to a combined 32 and 64 bit `woden_float` mode (orange squares) improves the accuracy to $\le 0.2$% on the longer baselines. The entirely 64 bit `woden_double` precision mode is consistent in precision across baseline length, sitting at < 2e-6% accuracy. The `woden_float` and `woden_double` executables are available in [Version 1.1](https://github.com/JLBLine/WODEN/releases/tag/v1.1.0), and can be switched between via a command line option in `run_woden.py`. It should be noted that these offset errors are deterministic, meaning comparison between different simulations out of [Version 1.0](https://github.com/JLBLine/WODEN/releases/tag/v1.0.0) `WODEN` are consistent; these errors matter most when comparing to real data.

As 32 and 64 bit precision calculations are performed in physically different parts of an NVIDIA GPU, with cards typically having less double precision hardware that single, the `woden_double` version is slower that the `woden_float`. Each card will show a different slow-down between the two modes. As a test, I ran a simulation using a catalogue of over 300,000 sources. The number of sources above the horizon and the simulation settings used are listed in Table \ref{tab:benchmark_sim}, along with the speed difference between the `woden_float` and `woden_double` executables for two different NVIDIA GPU cards.

\renewcommand{\arraystretch}{1}
\begin{table}[h]
\begin{center}
\begin{tabular}{l l}
\hline
Parameters & Value \\
\hline
\hline
Time steps & 14 \\
Frequency channels & 80 \\
Point sources components & 207673 \\
Gaussian components & 1182 \\
Shapelet components (basis functions) & 62 (10400) \\
Primary beam model & MWA FEE \\
GTX 1080 Ti \texttt{woden\_{}float} simulation time & 10min 39sec \\
GTX 1080 Ti \texttt{woden\_{}double} simulation time & 55min 46sec \\
V100 \texttt{woden\_{}float} simulation time & 4min 35sec \\
V100 \texttt{woden\_{}double} simulation time & 5min 55sec \\
\end{tabular}
\caption{Benchmark simulation to compare \texttt{woden\_{}float} and \texttt{woden\_{}double} speeds. Each shapelet component can have several basis function calculations, each more expensive that a point source component calculation. The MWA FEE is the most computationally expensive beam model included with \texttt{WODEN}.}
\label{tab:benchmark_sim}
\end{center}
\end{table}

Given this > 5 times slow down on a desktop card, having the option to toggle between `woden_float` and `woden_double` allows quick experimentation using `woden_float` and longer science-quality runs with `woden_double`. Luckily, for cards like the V100, the slowdown is around 1.3. Note that these simulations can easily be broken up and run across multiple GPUs if available, reducing the real time taken to complete the simulations.

# Example application

In @Line2020, we compared two methods to model Fornax A: a combination of point and elliptical Gaussians, compared to shapelets (see Figure \ref{fig:ForA}). We were able to quickly compare the computational efficiency of the methods using a desktop, and comment on their respective strengths and weaknesses in regard to foreground removal for EoR purposes. Furthermore, as we could control the simulations, we could compare the methods in the absence of other processing systematics that are present in the real data from the MWA, which dominated the comparison when using the `RTS` alone.

![Two methods to simulate Fornax A visibilities are compared here [both imaged using `WSClean` @Offringa2014; @Offringa2017], with point and elliptical Gaussians on the left, and shapelets on the right.\label{fig:ForA}](FornaxA_model_comparison.png)

# Documentation

The documentation for `WODEN` can be found on Read the Docs at [woden.readthedocs.io](https://woden.readthedocs.io/en/latest/), including a detailed installation guide, ways to test a local installation, details of the calculations `WODEN` makes under the hood, and worked examples, which are also included in the `GitHub` repo.

# Acknowledgements

I acknowledge direct contributions from Tony Farlie (who taught me how pointer arithmetic works in `C`) and contributions from Bart Pindor and Daniel Mitchell (through their work in the `RTS` and through advising me on `CUDA`). I would like to thank Chris Jordan who acted as a sounding board as I learned `C` and `CUDA`. Finally, I would like to thank both Matthew Kolopanis and Paul La Plante for reviewing the code and giving useful suggestions on how to improve the code.

This research was supported by the Australian Research Council Centre of Excellence for All Sky Astrophysics in 3 Dimensions (ASTRO 3D), through project number CE170100013. The International Centre for Radio Astronomy Research (ICRAR) is a Joint Venture of Curtin University and The University of Western Australia, funded by the Western Australian State government. This work was supported by resources provided by the Pawsey Supercomputing Centre with funding from the Australian Government and the Government of Western Australia.

# References
# Generating Code coverage
At the moment this is slightly stone-henge, due to the fact that CUDA is not supported for code coverage, and you can't automate a GPU test for free. For now, users must run the unit tests locally (using `ctest`, see [ReadTheDocs](https://woden.readthedocs.io/en/latest/testing/cmake_testing.html) for instructions). They can then run `source create_cov_reports.sh`, which will grab the outputs of the `C` tests from `ctest` and covert them into an appropriate formats, as well as running the `python` tests using `coverage`, which also generates appropriate outputs. To get `coverage`, do something like:

```bash
pip install coverage
```

Once those outputs have been created, you can then run `source send_reports_to_codecov.sh` to update the [codecov](https://about.codecov.io/) hosted coverage report. You will need an environment variable `WODEN_CODECOV_TOKEN` to be able to do this. Read the Community Guidelines (this doesn't exist yet) on how to get hold of that token (with great power...)

It's possible in the future that we can separate out the CUDA tests from the C/python tests, and automate that whole part, to automagically generate the whole "codecov" report.
## Example Simulations

For the `MWA_EoR1` and `EDA2_haslam` simulations, you need to download two extra sky models that are larger (total of about 100 MB). Follow the instructions [online here](https://woden.readthedocs.io/en/latest/examples/example_simulations.html) on how to run these simulations.
.. WODEN documentation master file, created by
   sphinx-quickstart on Mon Mar 22 14:37:03 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. _WSClean: https://sourceforge.net/projects/wsclean/
.. _Mitchell et al. 2008: https://ieeexplore.ieee.org/document/4703504?arnumber=4703504
.. _SHAMFI: https://github.com/JLBLine/SHAMFI
.. _Tingay et al. 2013: https://doi.org/10.1017/pasa.2012.007

The WODEN visibility simulator
=================================

``WODEN`` is C / CUDA code designed to be able to simulate low-frequency radio interferometric data. It is written to be simplistic and *fast* to allow all-sky simulations. Although ``WODEN`` was primarily written to simulate Murchinson Widefield Array (MWA, `Tingay et al. 2013`_) visibilities, it is becoming less instrument-specific as time goes on. `WODEN` outputs `uvfits` files.

``WODEN`` has been written with Stokes polarisations in mind, and as such, users can input a fully Stokes `I,Q,U,V` model, which is then propagated fully through the polarised instrumental response (depending on which primary beam you select), and output into Stokes `XX,XY,YX,YY` polarisations. Note however that currently only a power-law SED is supported, a rotation measure is currently unsupported.

The unique part of ``WODEN`` is that it can simulate shapelet model sources (along with point and Gaussian) that are compatible with the ``RTS`` (`Mitchell et al. 2008`_). These models are generated with SHApelet Modelling For Interferometers (`SHAMFI`_), specified with the ``--woden_srclist`` option. It also includes a script to convert a multi-scale CLEAN component list out of `WSClean`_ into a ``WODEN``-style srclist (when running ``WSClean`` use the ``-save-source-list`` option). ``WODEN`` can also produce visibilities that can be fed directly into the ``RTS`` to allow testing of calibration and modelling methodologies.

Documentation
-----------------

.. toctree::
   :maxdepth: 2

   installation/installation
   testing/cmake_testing
   testing/script_testing
   operating_principles/operating_principles
   examples/example_simulations
   API_reference/API_index

.. Indices and tables
.. --------------------
..
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
Testing installation via scripts
=================================

This is a straight-forward way to check your installation is working; just run some simple small simulations to check various functionality. Easiest way to check is to create images of the results. I've included a second set of scripts to convert the outputs to measurement sets and image them using ``WSClean``.

.. note:: I've tried to make these tests computationally low, so they need < 1.5 GB of GPU RAM, and < 8 GB system RAM. To do this, I've had to set ``--precision=float`` for some of the tests, to keep the memory requirements down. Running all the simulations will need about 600 MB of storage, with the imaging adding a further 800 MB (for a total of < 1.5 GB storage). The simulations should take < 2 minutes on most GPUs, and imaging less that 10 minutes for most CPUs (far less for fancier CPUs). I ran these tests fine on my laptop which has an Intel i7 2.8 GHz CPU, 16 GB system RAM, and an NVIDIA GeForce 940MX card with 2 GB RAM.

Running the simulations
------------------------

Once you've built ``WODEN``, and you have variables defined by calling ``WODEN/build/init_WODEN.sh``,
navigate to ``WODEN/test_installation``. To run all the tests immediately, you'll need to have obtained the defined the environment variable ``MWA_FEE_HDF5`` (see :ref:`Post compilation (optional)` for instructions on how to define that). To run all scripted test (including MWA FEE simulations)::

  $ cd WODEN/test_installation
  $ ./run_all_simulations.sh

This should rattle through a number of tests, with various combinations of component types (point, Gaussian, or shapelet) and different primary beams. If you don't want to run MWA FEE simulations, you can run::

  $ ./run_all_but_MWAFEE_simulations.sh

If you then want to run MWA FEE tests alone at a later date, you can run::

  $ ./run_only_MWAFEE_simulations.sh

If you want to incrementally run through tests, you can navigate through the ``single_component_models``, ``grid_component_models``, and ``different_beam_models`` directories to run each set individually.

Imaging the simulations
------------------------

Dependencies
^^^^^^^^^^^^^

To image the tests, you'll need a way to convert the uvfits files to measurement sets, and then image them. Which means dependencies (again).

+ **CASA** - https://casa.nrao.edu/casa_obtaining.shtml. Download an appropriate tarball (at the time 6.2 is the best version as it's ``python3``) and decompress it::

  $ wget https://casa.nrao.edu/download/distro/casa/release/rhel/casa-6.2.0-124.tar.xz
  $ tar -xvf casa-6.2.0-124.tar.xz

  That's it, it doesn't need installation
+ **WSClean** - https://wsclean.readthedocs.io/en/latest/installation.html. Head to this link to find out how to install ``WSClean``. You can of course use any other CLEANing software you want, but ``WSClean`` is excellent.

I'm assuming if you want to simulate interferometric data, you'll have some kind of FITS file imager already, but if not, ``DS9`` is a good place to start - https://sites.google.com/cfa.harvard.edu/saoimageds9.

Imaging scripts
^^^^^^^^^^^^^^^^

To use the imaging scripts, you **must** set an environment variable ``CASA_DIR`` to point towards the ``casa/bin`` of where you installed your ``casa``. For example, I did this::

  $ export CASA_DIR="/usr/local/casa-6.2.0-124/bin"

This will enable the scripts to convert the uvfits files to measurement sets. Once that variable is set, you can either image all the test outputs with::

  $ ./run_all_imaging.sh

or run the following as required::

  $ ./run_all_but_MWAFEE_imaging.sh
  $ ./run_only_MWAFEE_imaging.sh

Expected outcomes
------------------------
Nearly all of these simulations use the MWA phase 1 array layout, with a maximum baseline length of about 3 km, giving a resolution of around 3 arcmin.

Absolute Accuracy
^^^^^^^^^^^^^^^^^^^^^^^^
There is no imaging here, but runs an end-to-end simulation with a set of array layouts and sky models that should yield exact visibilities. The exact method is described in the JOSS paper.
The scripts run here are the exact scripts I used to create the plot in the JOSS paper.

.. todo:: Link the JOSS paper here once published

Single Component Models
^^^^^^^^^^^^^^^^^^^^^^^^

You should end up with three very basic images, each of just a single component of type point, Gaussian, and shapelet::

  $ cd WODEN/test_installation/single_component_models/images
  $ ls *image.fits
     single_gauss-image.fits  single_point-image.fits  single_shapelet-image.fits

which should look like

.. image:: single_component_plots.png
   :width: 600pt

For these simulations, I've switched off the primary beam, and set the spectral index to zero. I've also intentionally set the Gaussian and shapelet models to produce the same output, as the very first shapelet basis function is a Gaussian. All sources should have an integrated flux density of 1 Jy. If you're a sadist like me and still use ``kvis`` (https://www.atnf.csiro.au/computing/software/karma/) to look at FITS files, you can zoom into the source, and press 's' which will measure the integrated flux for you on the command line. This is quick and dirty, but gives us a good indication that the flux scale for all source types is working::

  points  mean mJy/Beam     std dev      min          max          sum
  2601     +16.917           +90.4008     -0.00196195  +999.997     +44001
  Total flux: 1000.00 mJy
  npoints  mean mJy/Beam     std dev      min          max          sum
  2601     +16.9164          +44.104      -0.110186    +264.247     +43999.5
  Total flux: 999.97 mJy
  npoints  mean mJy/Beam     std dev      min          max          sum
  2601     +16.916           +44.1038     -0.104652    +264.247     +43998.6
  Total flux: 999.95 mJy

This shows that we are within 50 micro Jy of the expected 1 Jy (taking into account that this is a CLEANed image with pixelisation effects).

Grid Component Models
^^^^^^^^^^^^^^^^^^^^^^^^

This should end up with three 5 by 5 grids, of the three component types::

  $ cd WODEN/test_installation/grid_component_models/images
  $ ls *image.fits
     grid_gauss-image.fits  grid_point-image.fits  grid_shapelet-image.fits

which should look like

.. image:: grid_component_plots.png
   :width: 600pt

The CLEAN isn't fantastic here as I've intentionally simulated a small amount of data to keep the size of the outputs down. But this tests that we can have multiple components and they are located at the requested positions (at every degree marker). I've included a very low-res model of PicA for the shapelet components, testing that we can have multiple shapelets with multiple basis functions. I've thrown in random position angles for the Gaussian and shapelets for a bit of variety.

Different Beam Models
^^^^^^^^^^^^^^^^^^^^^^^^

This should end up with a larger grid of a mix of components, with 4 different beam types (None, Gaussian, EDA2 (analytic dipole), and MWA FEE)::

  $ cd WODEN/test_installation/different_beam_models/images
  $ ls *image.fits
     multi-comp_grid_EDA2-image.fits      multi-comp_grid_MWA_FEE-image.fits
     multi-comp_grid_Gaussian-image.fits  multi-comp_grid_None-image.fits

The images with no beam, the Gaussian beam, and MWA FEE beam should look like this:

.. image:: different_beam_plots.png
   :width: 600pt

In the sky model, the top half are point sources, bottom left are shapelets, and bottom right are Gaussians. Again, limited data, so the CLEAN has some residuals. But we've successfully run a simulation with all three component types. We should also see different results for the Gaussian and MWA FEE beam plots, which we do, as we've used different primary beams. In particular I've made the Gaussian small enough of the sky to chop off the top left corner. The MWA FEE beam has a larger foot print.

For the EDA2 image, I've called the EDA2 array layout to override the settings in the metafits. The EDA2 has very short baselines, maximum of around 30 metres. If you compare the MWA phase 1 psf and the EDA psf we should be able to see the difference:

.. image:: MWA-vs-EDA2_psf.png
   :width: 600pt

This tests that we can override the array layout with a specified text file. Unsurprisingly, this turns our EDA2 image of the same sky model into a bunch of blobs:

.. image:: EDA2_layout_image.png
   :width: 300pt

but this is what we expect. That's it for the simple installation tests. If you want to really test out the simulation capabilities of ``WODEN``, check out the :ref:`WODEN demonstrated via examples`  section, which has bigger and better simulations.

Deleting test outputs
------------------------
If you don't want a bunch of files hanging around on your system for no reason, just run::

  $ ./delete_sim_outputs.sh
  $ ./delete_images.sh

which will nuke the outputs for you.
Testing via ``ctest``
======================

This is totally optional, but you can run the unit / integration tests I use for
developing ``WODEN`` to check your system / the code runs as expected. The tests
are compiled using the same ``cmake`` script as the main code.

ALL tests are compiled in both FLOAT and DOUBLE precision, meaning the same
test code is used to test the two different precision versions, with the compiler
dropping in the necessary precision. Unless explicitly noted in the details
in :ref:`What do the tests actually do?`, the tests require the same level
of accuracy from the FLOAT and DOUBLE precision versions. Most of the time
the FLOAT version is accurate to an absolute tolerance of 1e-7, and the DOUBLE
to 1e-15. Read the tolerances for specific functions in the sections listed in
:ref:`What do the tests actually do?`.

Dependencies
-------------

The tests use the following C library:

* **Unity** - https://github.com/ThrowTheSwitch/Unity

I installed ``unity`` via::

  $ git clone https://github.com/ThrowTheSwitch/Unity.git
  $ cd Unity
  $ mkdir build && cd build
  $ cmake ..
  $ make -j 4
  $ sudo make install

That way, ``cmake`` can find ``unity``. However, you don't need to install Unity anywhere, as ``WODEN`` uses the ``C`` code directly. You just need to tell ``WODEN`` where Unity lives in your system (for example, you could download a release version e.g. version 2.5.2 - see example below for how to link without installation).

You'll also need to initiate the ``git submodule`` that runs code coverage. Simply navigate to the ``WODEN`` directory and run::

  $ git submodule init
  $ git submodule update

which will pull in the relevant ``CMake-codecov`` dependencies. This allows us to track code coverage for the ``python`` and ``C`` code (no free tools exist for ``CUDA`` at the time of writing, boooo).

To tell ``cmake`` to build tests, you add ``TARGET_GROUP=test`` to your command to tell CMake to build tests instead of the main code::

  $ cd $WODEN_DIR
  $ cmake .. -DTARGET_GROUP=test -DUNITY_ROOT=/usr/local/Unity-2.5.2
  $ make -j 4

(where ``-DUNITY_ROOT=`` is needed if you didn't install ``unity``).

.. warning:: once you have done this, to go back to compiling the main ``woden`` executable, you need to run::

    $ cmake .. -DTARGET_GROUP=production

    otherwise you'll just keep building the tests.

Running tests
--------------

Once that compiles, you can run the tests by running::

  $ ctest

You should see something like the following if successful::

  $ ctest
  Test project /home/jline/software/WODEN/build
         Start  1: C_test_RTS_ENH2XYZ_local_float
    1/87 Test  #1: C_test_RTS_ENH2XYZ_local_float .......................   Passed    0.00 sec
         Start  2: C_test_calc_XYZ_diffs_float
    2/87 Test  #2: C_test_calc_XYZ_diffs_float ..........................   Passed    0.00 sec
         Start  3: C_test_RTS_PrecessXYZtoJ2000_float
    3/87 Test  #3: C_test_RTS_PrecessXYZtoJ2000_float ...................   Passed    0.00 sec
         Start  4: C_test_RTS_ENH2XYZ_local_double
    4/87 Test  #4: C_test_RTS_ENH2XYZ_local_double ......................   Passed    0.00 sec
         Start  5: C_test_calc_XYZ_diffs_double
    5/87 Test  #5: C_test_calc_XYZ_diffs_double .........................   Passed    0.00 sec
         Start  6: C_test_RTS_PrecessXYZtoJ2000_double
    6/87 Test  #6: C_test_RTS_PrecessXYZtoJ2000_double ..................   Passed    0.00 sec
         Start  7: C_test_null_comps_float
    7/87 Test  #7: C_test_null_comps_float ..............................   Passed    0.00 sec
         Start  8: C_test_fill_chunk_src_with_pointgauss_float
    8/87 Test  #8: C_test_fill_chunk_src_with_pointgauss_float ..........   Passed    0.03 sec
         Start  9: C_test_fill_chunk_src_with_shapelets_float
   ...etc etc

.. note:: To test MWA Fully Embedded Element beam code, you must have the environment variable::

    MWA_FEE_HDF5=/path/to/mwa_full_embedded_element_pattern.h5

  declared. If you don't have that, the MWA FEE beam code tests will just be skipped.

If you want more detail of what the tests are doing, run::

  $ ctest --verbose

What do the tests actually do?
---------------------------------

The tests are all located in ``WODEN/cmake_testing``, and each directory within contains tests
for a different file from ``WODEN/src``. Within each test directory, there are separate files for testing different functions, which include the function name. As an example, the directory ``WODEN/cmake_testing/array_layout`` contains tests for the file ``WODEN/src/array_layout.c``, and contains test files that test the following functions::

  cmake_testing/array_layout/test_calc_XYZ_diffs.c -> src/array_layout.c::calc_XYZ_diffs
  cmake_testing/array_layout/test_RTS_ENH2XYZ_local.c -> src/array_layout.c::RTS_ENH2XYZ_local
  cmake_testing/array_layout/test_RTS_PrecessXYZtoJ2000.c -> src/array_layout.c::RTS_PrecessXYZtoJ2000

The ``C`` and ``CUDA`` functions are tested using the `Unity`_ library, which has useful functions like::

  TEST_ASSERT_FLOAT_EQUAL();
  TEST_ASSERT_DOUBLE_WITHIN();
  TEST_ASSERT_NULL();

allowing a simple testing of values. If a test says outputs are tested to be
equal, it refers to the ``TEST_ASSERT_FLOAT_EQUAL`` or ``TEST_ASSERT_DOUBLE_EQUAL``
function.

.. _`Unity`: https://github.com/ThrowTheSwitch/Unity

.. note:: For those unfamiliar with ``CMake`` testing, even though the tests are located in ``WODEN/cmake_testing/``, when you run ``ctest``, the test files are copied and run in ``WODEN/build/cmake_testing``, so any output from the tests will be located there.

The sections below give an outline of the tests performed in each directory.

``C`` code tests:

.. toctree::
   :maxdepth: 1

   cmake_testing/array_layout
   cmake_testing/chunk_sky_model
   cmake_testing/create_sky_model
   cmake_testing/FEE_primary_beam
   cmake_testing/primary_beam
   cmake_testing/shapelet_basis
   cmake_testing/visibility_set
   cmake_testing/woden_settings

``CUDA`` code tests:

.. toctree::
   :maxdepth: 1

   cmake_testing/calculate_visibilities
   cmake_testing/FEE_primary_beam_cuda
   cmake_testing/fundamental_coords
   cmake_testing/primary_beam_cuda
   cmake_testing/source_components

``python`` code tests:

.. toctree::
   :maxdepth: 1

   cmake_testing/run_woden

.. note:: To be able to test ``CUDA`` functions that are designed to work solely in GPU memory, it's necessary to write wrapper functions that allocate GPU memory, pass the data into the ``CUDA`` code to be tested, and then copy the results back into host memory. I've kept these 'intermediate' test functions inside the ``*.cu`` files that contain the code being tested, as it's not straight forward / performance degrading to have them in separate files. On casual inspection it looks like there are many functions in the ``*.cu`` files I haven't written tests for, but the extra functions are there *because* of testing. Sigh.
``create_sky_model``
=========================
Tests for the functions in ``WODEN/src/create_sky_model.c``. The functions
read in the sky model from a text file, and crop anything below the horizon.

test_read_source_catalogue.c
*********************************
``create_sky_model::read_source_catalogue`` reads in the sky model from a text
file. Each test reads in a sky model, uses ``read_source_catalogue`` to
create a sky model, and tests the correct information has been read in.
``read_source_catalogue`` returns an integer as an error message (0 good, 1 bad),
so some text files below are written to cause failure. See the table below
for each test sky model and the expected result. The way each SOURCE and it's
associated COMPONENTs are stored can affect the way the sky model is cropped,
so all tests check that generated ``source_catalogue_t`` struct is structured
correctly. For all floating point values, when compiling in FLOAT mode, test
asserts that values are within an absolute tolerance of 1e-7, and 1e-15 when
compiling in DOUBLE mode.

.. list-table::
   :widths: 25 50 25
   :header-rows: 1

   * - Sky model
     - Test case
     - Test outcomes
   * - srclist_no-comp_numbers.txt
     - Is missing the line that contains number of COMPONENTs in the SOURCE.
     - Check fails
   * - srclist_badcoeff.txt
     - Has a bad SCOEFF line where one number is sdfasdfasdfasdfasdfasf
     - Check fails
   * - srclist_badspell.txt
     - Contains an incorrect spelling of COMPONENT
     - Check fails
   * - srclist_singlegauss.txt
     - Contains a single GAUSSIAN SOURCE
     - Check sky model values
   * - srclist_singlepoint.txt
     - Contains a single POINT SOURCE
     - Check sky model values
   * - srclist_comment.txt
     - Contains a commented line (and a single POINT SOURCE)
     - Check sky model values
   * - srclist_singleshape.txt
     - Contains a single SHAPELET SOURCE
     - Check sky model values
   * - srclist_empty_line.txt
     - Contains an empty line  (and a single POINT SOURCE)
     - Check sky model values
   * - srclist_threecomponents.txt
     - Contains one SOURCE with three COMPONENTs
     - Check sky model values
   * - srclist_mulitple_source-components.txt
     - Contains multiple SOURCEs each with multiple COMPONENTs
     - Check sky model values
   * - srclist_threesources.txt
     - Contains multiple SOURCEs each with a single COMPONENT
     - Check sky model values
   * - srclist_linear.txt
     - Contains a single POINT COMPONENT with the LINEAR keyword specifying the SED
     - Check sky model values


test_horizon_test.c
*********************************
``create_sky_mode::horizon_test`` takes information on a single COMPONENT of a
SOURCE and tests whether it is above or below the horizon. Depending on
whether we are cropping the sky model by SOURCE or by COMPONENT, it updates
various counters that effect the sky model cropping. (cropping by SOURCE
throws away the whole SOURCE if one COMPONENT is below the horizon, cropping
by COMPONENT only throws away the COMPONENTs below the horizon). SHAPELETs are complicating
factors as a single position can match multiple basis function parameters (of
any length) so ``horizon_test`` does some logic to count how many SHAPELET
parameters are being retained.

First set of tests check that cropping on POINT/GAUSSIAN type COMPONENTs work for
both cropping by SOURCE and cropping by COMPONENT. The second set of tests check
that when cropping by COMPONENT, and we have multiple SHAPELET COMPONENTs, that
the correct number of SHAPELET basis functions parameters are retained.

test_crop_sky_model.c
*********************************
``create_sky_mode::crop_sky_model`` calculates the azimuth / zenith angle of
all COMPONENTs in a sky model (for the first LST step of a simulation), and then
crops out either SOURCEs or COMPONENTs that are below the horizon (by using
``horizon_test``). Once cropped, it then calculates the az/za for all time
steps in the simulation for the surviving COMPONENTs.

The tests here are split across 4 sky models - just POINTs, just GAUSSIANs,
just SHAPELETs, and a mix of all three. For each sky model type, the tests are
run for both the cropping by SOURCE and cropping by COMPONENT case. All tests
are run for four different LSTs - a total of 32 tests. The RA/Dec of the
COMPONENTs and the LSTs are chosen in such a way that different combinations
should result in different COMPONENT/SOURCEs being discarded. All tests check
that the correct COMPONENT/SOURCEs are retained (including the SHAPELET basis
function parameters), and that the correct az/za are calculated for each time
step. The sky model setup and expected results are stored in ``test_crop_sky_model.h``.

The azimuth and zenith angle outputs are tested to match expectations to within
an absolute tolerance of 1e-6 for the FLOAT compiled code, and 1e-12 for the
DOUBLE.
``primary_beam_cuda``
=========================
Tests for the functions in ``WODEN/src/primary_beam_cuda.cu``. These functions
calculate the beam responses for the EDA2 and Gaussian beam models.

test_gaussian_beam.c
*********************************
This calls ``primary_beam_cuda::test_kern_gaussian_beam``, which in turn
tests ``primary_beam_cuda::kern_gaussian_beam``, the kernel that calculates
the Gaussian primary beam response. As a Gaussian is an easy function to
calculate, I've setup tests that calculate a north-south and east-west strip
of the beam response, and then compare that to a 1D Gaussian calculation.

As ``kern_gaussian_beam`` just takes in *l,m* coords, these tests just generate
100 *l,m* coords that span from -1 to +1. The tests check whether the kernel
produces the expected coordinates in the *l* and *m* strips, as well as changing
with frequency as expected, by testing 5 input frequencies with a given
reference frequency. For each input frequency :math:`\nu`, the output is
checked against the following calculations:

 - When setting *m* = 0, assert gain = :math:`\exp\left[-\frac{1}{2} \left( \frac{l}{\sigma} \frac{\nu}{\nu_0} \right)^2 \right]`
 - When setting *l* = 0, assert gain = :math:`\exp\left[-\frac{1}{2} \left( \frac{m}{\sigma} \frac{\nu}{\nu_0} \right)^2 \right]`

where :math:`\nu_0` is the reference frequency, and :math:`\sigma_0` the std of
the Gaussian in terms of *l,m* coords. These calculations are made using ``C``
with 64 bit precision.  The beam responses are tested to be within an absolute
tolerance of 1e-10 from expectations for the FLOAT compiled code, and 1e-16 for
the DOUBLE compiled code.

test_analytic_dipole_beam.c
***********************************
This calls ``primary_beam_cuda::test_analytic_dipole_beam``, which in turn
tests ``primary_beam_cuda::calculate_analytic_dipole_beam``, code that copies
az/za angles into GPU memory, calculates an analytic dipole response toward
those directions, and then frees the az/za coords from GPU memory.

Nothing exiting in this test, just call the function for 25 directions on
the sky, for two time steps and two frequencies (a total of 100 beam calculations),
and check that the real beam gains match stored expected values, and the imaginary
values equal zero. The expected values have been generated using the DOUBLE
precision compiled code, and so the absolute tolerance of within 1e-12 is set
by how many decimal places I've stored in the lookup table. The FLOAT precision
must match within 1e-6 of these stored values.
``shapelet_basis``
=========================
Tests for the functions in ``WODEN/src/shapelet_basis.c``.

``test_create_sbf.c``
****************************
``shapelet_basis::create_sbf`` just creates a massive look-up table of shapelet
basis function values. Here we call the function and test 20 different
array indexes and assert that they are equal to expected values. Exciting stuff.

For FLOAT compiled code, the absolute tolerance threshold on values is set to
1e-7, and 1e-15 for DOUBLE compiled code.
``calculate_visibilities``
===========================
Tests for the functions in ``WODEN/src/calculate_visibilities.cu`` (there is only one function). ``calculate_visibilities::calculate_visibilities`` is the gateway function
to all ``CUDA`` functionality in ``WODEN``. It takes in simulations settings and
a sky model, and performs the necessary coorindate and measurement equation calculations, as well as summations over sky model components to generate visibilities.

The tests below are fairly limited in the number of variables tested, as they
are integration tests to make sure all the functionality expected is called by ``calculate_visibilities::calculate_visibilities``, and is able to talk to the
GPU correctly. More rigorous testing of this functionality is included in other
test suites in ``cmake_testing``.

The tests below all use ``test_calculate_visibilities_common.c`` to setup a
``source_catalogue_t`` sky model struct with the following sky model setups.
These differing sky models should cover all COMPONENT combinations possible, to
make sure ``calculate_visibilities`` is calling the appropriate functions:

.. list-table::
   :widths: 30 30 30 30 30
   :header-rows: 1

   * - Num SOURCEs
     - Num POINT in each SOURCE
     - Num GAUSS in each SOURCE
     - Num SHAPELET in each SOURCE
     - Total COMPONENTS
   * - 1
     - 1
     - 0
     - 0
     - 1
   * - 1
     - 0
     - 1
     - 0
     - 1
   * - 1
     - 0
     - 0
     - 1
     - 1
   * - 1
     - 1
     - 1
     - 1
     - 3
   * - 3
     - 1
     - 0
     - 0
     - 3
   * - 3
     - 0
     - 1
     - 0
     - 3
   * - 3
     - 0
     - 0
     - 1
     - 3
   * - 3
     - 1
     - 1
     - 1
     - 9
   * - 3
     - 3
     - 0
     - 0
     - 9
   * - 3
     - 0
     - 3
     - 0
     - 9
   * - 3
     - 0
     - 0
     - 3
     - 9
   * - 3
     - 3
     - 3
     - 3
     - 27

All sources are given an *RA,Dec* equal to the phase centre, a spectral index
of zero, and a flux density of 0.3333333333333333 Jy. This way, the output visibilities
(in the absence of a primary beam) should be fully real, and equal to the sum of the number of
COMPONENTs in the sky model multiplied by 0.3333333333333333. This numerical 1/3
flux is a good test of the precision of the FLOAT and DOUBLE compiled codes.

GAUSSIAN and SHAPELET components with any size cause a reduction of the real part
of the visibility due to the extended size on the sky. For these tests, I've set
their major and minor axes to 1e-10 to make them behave similarly to point sources.

For all tests below, I setup a simulation with three baselines, three frequencies,
and two timesteps. For all sky models, frequencies, time steps, and baselines, I check:

 - The *u,v,w* are as expected
 - The real part of the visibilities are equal to number of COMPONENTs times 0.3333333333333333 Jy (modulu the beam repsonse - see below)

To keep testing the *u,v,w* straight forward, I've set the baseline lengths in :math:`X` and :math:`Y` equal, (i.e. :math:`X_{\mathrm{diff}} = X_{\mathrm{ant1}} - X_{\mathrm{ant2}} = Y_{\mathrm{diff}}`), and the length in :math:`Z` to zero. With this configuration, the
following is true:

 - :math:`u = X_{\mathrm{diff}}(\cos(ha_0) + \sin(ha_0))`
 - :math:`v = X_{\mathrm{diff}}\sin(\phi_{\mathrm{lat}})(-\cos(ha_0) + \sin(ha_0))`
 - :math:`w = X_{\mathrm{diff}}\cos(\phi_{\mathrm{lat}})(\cos(ha_0) - \sin(ha_0))`

where :math:`ha_0` is the hour angle of the phase centre, and :math:`\phi_{\mathrm{lat}}`
the latitude of the array. The allows us to check the *u,v,w* are changing with time.

For FLOAT compiled code, the absolute tolerance threshold on the *u,v,w*
values is set to 1e-5, and 1e-12 for DOUBLE compiled code.

``test_calculate_visibilities_nobeam.c``
*********************************************
This runs the tests explained above, whilst switching the primary beam off. This
really does check that real visibilities are equal to the number of COMPONENTs
times 0.3333333333333333 for all 12 sky model configurations, and the imaginary
equal to zero.

For FLOAT compiled code, the absolute tolerance threshold on
values is set to 1e-6, and 1e-9 for DOUBLE compiled code.

``test_calculate_visibilities_gaussbeam.c``
*********************************************
This runs the same tests as ``test_calculate_visibilities_nobeam.c``, but applies
a Gaussian primary beam model. As all sky model COMPONENTs are set at the same location,
only one beam gain per time step should be applied to visibilities, so for each time
step, we can check whether the visibilities equal this gain times the number of
COMPONENTs. The Gaussian beam is fully real and has no cross-pol values, so only
check that the real XX and YY visibilities have value, and ensure all other
visibility information is zero.

For FLOAT compiled code, the absolute tolerance threshold on
values is set to 1e-5, and 1e-8 for DOUBLE compiled code.

``test_calculate_visibilities_edabeam.c``
*********************************************
Runs exactly the same tests as ``test_calculate_visibilities_gaussbeam.c``, but
using the analytic single dipole beam (the EDA2 beam). Obviously tests the
visiblities match the gain values that are expected for the EDA2 beam and not
the Gaussian Beam test.

For FLOAT compiled code, the absolute tolerance threshold on
values is set to 1e-5, and 1e-8 for DOUBLE compiled code.

``test_calculate_visibilities_mwafeebeam.c``
*********************************************
Again, runs the same tests as ``test_calculate_visibilities_gaussbeam.c``, but
this time for the coarse resolution MWA FEE primary beam. As this model is
complex and includes mutual coupling, both the real and imaginary values
for all XX, XY, YX, and YY polarisations are tested.

For FLOAT compiled code, the absolute tolerance threshold on
values is set to 4e-3, and 1e-7 for DOUBLE compiled code. The expected gains
of the MWA FEE beam are taken from the DOUBLE compiled code, and so the large
threshold for the FLOAT here is mostly due to the inaccuracy of the FLOAT
MWA FEE beam code (see :ref:`FEE_primary_beam_cuda_cmake` for more discussion on this).
``primary_beam``
=========================
Tests for the functions in ``WODEN/src/primary_beam.c``. These functions
setup primary beam settings, ready to calculate beam responses on the GPU.

test_calc_para_angle.c
*********************************
``primary_beam::calc_para_angle`` calculates the parallactic angle for
all COMPONENTs, for all time steps. This test calls ``calc_para_angle`` for
three POINT, three GAUSSIAN, and three SHAPELET COMPONENTs, for three different
LSTs. The different COMPONENT types are given different RA/Decs. The
results are compared to expected values, which are stored in
``expected_para_angles.h``. For FLOAT compiled code, the absolute tolerance
threshold is set to 1e-7, and 1e-15 for DOUBLE compiled code.


test_fill_primary_beam_settings.c
***********************************
``primary_beam::fill_primary_beam_settings`` prepares a ``beam_settings_t``
struct to be used by ``calculate_visibilities::calculate_visibilities``. The
az,za coords have already been calculated by
``chunk_sky_model::create_chunked_sky_models``, which is sufficient for some
beam models. There are two beam models that need further inputs:

   - ``GAUSS_BEAM``: the Gaussian beam function uses *l,m,n* coords to incorporate projection effects, so ``fill_primary_beam_settings`` calculates the hour angle of all COMPONENTs for all time steps, to feed into ``fundamental_coords::kern_calc_lmn`` later down the line
   - ``FEE_BEAM``: needs the parallactic angle to rotate the telescope-based coords into the Stokes frame

The tests here call ``fill_primary_beam_settings`` for the four primary
beam types, and perform the following checks:

 - ``GAUSS_BEAM``:
    - Assert that ``beam_settings->beamtype == GAUSS_BEAM``
    - Assert that a number of constants are copied from ``woden_settings`` into ``beam_settings``
    - Assert that::

        src->point_gaussbeam_decs
        src->point_gaussbeam_has
        src->gauss_gaussbeam_decs
        src->gauss_gaussbeam_has
        src->shape_gaussbeam_decs
        src->shape_gaussbeam_has

      have been set to the correct values for all COMPONENTs and time steps
 - ``FEE_BEAM``:
    - Assert that ``beam_settings->beamtype == FEE_BEAM``
    - Assert that::

        src->sin_point_para_angs
        src->cos_point_para_angs
        src->sin_gauss_para_angs
        src->cos_gauss_para_angs
        src->sin_shape_para_angs
        src->cos_shape_para_angs

      have been set to the correct values for all COMPONENTs and time steps
 - ``ANALY_DIPOLE``:
    - Assert that ``beam_settings->beamtype == ANALY_DIPOLE``
 - ``NO_BEAM``:
    - Assert that ``beam_settings->beamtype == NO_BEAM``

For FLOAT compiled code, the absolute tolerance threshold on values is set to
1e-7, and 1e-15 for DOUBLE compiled code. For values that are 64 bit in both the
FLOAT and DOUBLE versions the values are tested using ``TEST_ASSERT_EQUAL_DOUBLE``.
``FEE_primary_beam``
=========================
Tests for the functions in ``WODEN/src/FEE_primary_beam.c``.

test_RTS_MWAFEEInit.c
*********************************
``create_sky_model::RTS_MWAFEEInit`` reads in stored spherical harmonic
coefficients from ``mwa_full_embedded_element_pattern.h5`` and stores them in an ``RTS_MWA_FEE_beam_t`` struct, to be used later in MWA FEE beam calculations.
There are four generated arrays that matter, which are::

  RTS_MWA_FEE_beam_t->Q1 (double _Complex)
  RTS_MWA_FEE_beam_t->Q2 (double _Complex)
  RTS_MWA_FEE_beam_t->M (double)
  RTS_MWA_FEE_beam_t->N (double)

The MWA beam pointing direction on the sky is controlled by a set of 16 delays.
A different delay setting is stored in a different table in the ``hdf5`` file.
In these tests, three different delays settings are tested at 50MHz, 150MHz, and
250MHz (a total of nine tests). For each combination of settings, the values
of the four arrays are tested at 6 locations against stored known values (
which can be found in ``test_RTS_MWAFEEInit.h``), and must be within an absolute
tolerance of 1e-7 to pass. Only six array entries each are tested purely
to keep the stored values to compare to down a reasonable number (but suffice
as a test that things are working correctly).
``visibility_set``
=========================
Tests for the functions in ``WODEN/src/visibility_set.c``. These functions handle
a ``visibility_set_t`` struct, which holds the output visibilities. Functions
here include mallocing, filling, and freeing attributes.

``test_fill_timefreq_visibility_set.c``
*****************************************
This calls ``visibility_set::fill_timefreq_visibility_set``, which uses a
populated ``woden_settings_t`` struct to fill in the following attributes in
a ``visibility_set_t`` struct::

  visibility_set->allsteps_sha0s
  visibility_set->allsteps_cha0s
  visibility_set->allsteps_lsts
  visibility_set->allsteps_wavelengths

Where ``sha0`` and ``cha0`` are sine and cosine of the hour angle of the phase
centre, respectively. This test runs with three time steps with LSTs of
:math:`0`, :math:`\pi/6`, :math:`\pi/4`, so that ``allsteps_sha0s`` and
``allsteps_cha0s`` can be analytically predicted. ``allsteps_wavelengths`` are
calculated using input frequencies, which are set to :math:`c/2`, :math:`3c/4`,
and :math:`c`, meaning the expected output wavelengths should be :math:`2`, :math:`4/3`,
and :math:`1`. The third let's us test the accuracy of the wavelength calculation.

All angles and LSTs are stored at 64bit precision, so both FLOAT and DOUBLE
code versions are tested to within an absolute tolerance of 1e-15. The precision
of ``allsteps_wavelengths`` is set by the user, and is tested to within 1e-7
for FLOAT and 1e-15 for DOUBLE.

``test_malloc_and_free.c``
*****************************************
Very basic test of ``visibility_set::setup_visibility_set``,
``visibility_set::free_visi_set_inputs``, and ``visibility_set::free_visi_set_outputs``,
which are functions that either ``malloc`` or ``free`` specific attributes in a
``visibility_set_t`` struct. Tests by calling each function and checking that
the following attributes are NOT a NULL if a ``malloc`` was called, and the correct
attributes are NULL if ``free`` was called.

``test_write_visi_set_binary.c``
*****************************************
Tests ``visibility_set::write_visi_set_binary``, which writes out the contents
of a ``visibility_set_t`` struct to a binary file. Tests by filling a
``visibility_set_t`` struct with some simple non-repeating values, and then
uses ``write_visi_set_binary``. Then reads that binary file in and checks the
contents matched what was in the ``visibility_set_t`` struct.

``test_write_visi_set_text.c``
*****************************************
Tests ``visibility_set::write_visi_set_text``, which writes out a subset of
the contents of a ``visibility_set_t`` struct to a text file. This test
works similarly to :ref:`test_write_visi_set_binary.c`, by calling
``write_visi_set_text`` with a known set of inputs, and checking the text file
it writes out contains the known inputs.
``run_woden``
=========================
Tests for the functions in ``WODEN/src/run_woden.py``. These functions handle:
parsing user arguments; calculating astronomical constants;
reading in variables from an MWA metafits if requested; writing out a ``.json``
file to input into the ``woden`` executable; calling the ``woden`` executable;
reading the binary file written out by ``woden``; creating a ``uvfits`` file;
tidying up after the simulator.

In the following, ``rw`` is short for ``run_woden.py``, and so ``rw.calc_jdcal``
means the function ``calc_jdcal`` in ``run_woden.py``.

test_command.py
*******************************************************
Tests the ``rw.command`` function, which should call things on the command line
from within ``python``. Test by running the command::

   $ echo cheese > example.txt

and then reading the word "cheese" out of the file "example.txt" that should
have been created.

test_calc_jdcal.py
*******************************************************
Tests the ``rw.calc_jdcal`` function, which should calculate the Julian Date and
split into a integer day and fractional day values. Just test by calling
``rw.calc_jdcal`` with two known date strings, and checking the output values
match expectations.

test_get_uvfits_date_and_position_constants.py
*******************************************************
Tests the ``rw.get_uvfits_date_and_position_constants`` function,
which should calculate the LST, GST0 (greenwich sidereal time at 0 hours
of the given date), DEGPDY (rotational speed of the Earth) and UT1UTC (
difference between UT1 and UTC) for a given Long/Lat/Height array location and
UTC date. Test with two combinations of different UTC/Long/Lat/Height and
check the returned values are as expected.

test_RTS_encoding.py
*******************************************************
Tests the ``rw.RTS_encode_baseline`` function, which should take two antenna
numbers and create the ``BASELINE`` number as per AIPS uvfits file convention.
Tests by running with four antenna pairs, and ensuring the output values match
expectation.

Secondly, tests the function ``rw.RTS_decode_baseline``, which should separate
the encoded BASELINE number back into two antenna numbers. Test by decoding the
same four BASELINE numbers and ensuring the correct antenna numbers are found.

test_make_antenna_table.py
*******************************************************
Tests the ``rw.make_antenna_table`` function, which should create
the antenna table that goes into a uvfits file. Test by giving it a
known set of input parameters and checking those values end up in
the correct location and format of output antenna table. The following parameters
are tested as correct:

.. code-block:: python

   ant_table.data['ANNAME']
   ant_table.data['STABXYZ']
   ant_table.data['ORBPARM']
   ant_table.data['NOSTA']
   ant_table.data['MNTSTA']
   ant_table.data['STAXOF']
   ant_table.data['POLTYA']
   ant_table.data['POLAA']
   ant_table.data['POLCALA']
   ant_table.data['POLTYB']
   ant_table.data['POLAB']
   ant_table.data['POLCALB']
   ant_table.header['ARRAYX']
   ant_table.header['ARRAYY']
   ant_table.header['ARRAYZ']
   ant_table.header['FREQ']
   ant_table.header['GSTIA0']
   ant_table.header['DEGPDY']
   ant_table.header['UT1UTC']
   ant_table.header['XYZHAND']
   ant_table.header['FRAME']
   ant_table.header['RDATE']
   ant_table.header['TIMSYS']
   ant_table.header['ARRNAM']
   ant_table.header['NUMORB']
   ant_table.header['NOPCAL']
   ant_table.header['POLTYPE']
   ant_table.header['CREATOR']

test_create_uvfits.py
*******************************************************
Tests the ``rw.create_uvfits`` function, which should take a whole
heap of inputs and write out a uvfits. Test by running function with a known
set of inputs, reading in the created file, and checking contents match
expectations. Along checking all the same parameters in the antenna table as
checked by :ref:`test_make_antenna_table.py`, the following parameters are
checked against the inputs:

.. code-block:: python

    data_table.data['UU']
    data_table.data['VV']
    data_table.data['WW']
    data_table.data['BASELINE']
    ##Astropy automatically adds the header value to the DATE array,
    ##so need to subtract before comparison
    data_table.data['DATE'] - data_table.header['PZERO4']
    ##Check the actual visisbility values are correct
    data_table.data.data
    data_table.header['CTYPE2']
    data_table.header['CRVAL2']
    data_table.header['CRPIX2']
    data_table.header['CDELT2']
    data_table.header['CTYPE3']
    data_table.header['CRVAL3']
    data_table.header['CRPIX3']
    data_table.header['CDELT3']
    data_table.header['CTYPE4']
    data_table.header['CRVAL4']
    data_table.header['CRPIX4']
    data_table.header['CDELT4']
    data_table.header['CTYPE5']
    data_table.header['CRVAL5']
    data_table.header['CRPIX5']
    data_table.header['CDELT5']
    data_table.header['CTYPE6']
    data_table.header['CRVAL6']
    data_table.header['CRPIX6']
    data_table.header['CDELT6']
    data_table.header['PSCAL1']
    data_table.header['PZERO1']
    data_table.header['PSCAL2']
    data_table.header['PZERO2']
    data_table.header['PSCAL3']
    data_table.header['PZERO3']
    data_table.header['PSCAL4']
    data_table.header['PZERO4']
    data_table.header['PSCAL5']
    data_table.header['PZERO5']
    data_table.header['OBJECT']
    data_table.header['OBSRA']
    data_table.header['OBSDEC']
    data_table.header['GITLABEL']
    data_table.header['TELESCOP']
    data_table.header['LAT']
    data_table.header['LON']
    data_table.header['ALT']
    data_table.header['INSTRUME']

test_enh2xyz.py
*******************************************************
Tests the ``rw.enh2xyz`` function, which should calculate the local X,Y,Z coords
using the local east, north, height. Test using the cases where latitude is 0
and -30 deg, which have analytically predictable outcomes. This runs the same
test as descibred in :ref:`test_RTS_ENH2XYZ_local.c`.

test_load_data.py
*******************************************************
Tests the ``rw.load_data`` function, which should read in a binary
file as output by ``woden_float`` or ``woden_double``, into various arrays.
Test by writing out a binary file with known input params, reading in that
binary file using ``rw.load_data``, and comparing the inputs to outputs. The
test is run in both 32 and 64 bit precision.

test_write_json.py
*******************************************************
Test the ``rw.write_json`` function, which writes an input file to feed into
either ``woden_float`` or ``woden_double``. A number of tests are run, all of
which call ``rw.write_json`` using a minimum set of example input arguments.
The resulting ``.json`` is then read back in, and the following parameters are
checked as correct:

.. code-block:: python

   json_data['ra0']
   json_data['dec0']
   json_data['num_freqs']
   json_data['num_time_steps']
   json_data['cat_filename']
   json_data['time_res']
   json_data['frequency_resolution']
   json_data['chunking_size']
   json_data['jd_date']
   json_data['LST']
   json_data['array_layout']
   json_data['lowest_channel_freq']
   json_data['latitude']
   json_data['coarse_band_width']
   json_data['band_nums']

The following tests run with the following optional arguments:

 - ``test_write_gaussian_beam``: checks that extra arguments that control the Gaussian primary beam are written correctly
 - ``test_write_MWA_FEE_beam``: checks that extra arguments that control the MWA FEE beam are written correctly
 - ``test_write_EDA2_beam``: checks that extra arguments that control the EDA2 beam are written correctly
 - ``test_write_no_precession``: checks that the option to turn off precession is added when asked for

test_make_baseline_date_arrays.py
*******************************************************
Tests the ``rw.make_baseline_date_arrays`` function, which should make
the DATE and BASELINE arrays that are needed to populate a uvfits file. Test
by giving the function a known date string, number of antennas, number of time
steps, and time resolution, and checking the output arrays match expectations.

test_remove_phase_tracking.py
*******************************************************
Tests the ``rw.remove_phase_tracking`` function, which should remove
the phase tracking applied to visibilities. The original MWA correlator
did not phase track, so the ``RTS`` expects no phase tracking on the data, so
to input ``WODEN`` simulations into the ``RTS``, have to undo the phase-tracking.
The ``RTS`` calculates it's own ``u,v,w``, so I only fiddle the visibilities
here so be warned.

This test starts by creating a random array layout via:

.. code-block:: python

  num_antennas = 50
  ##Make a random array layout
  east = np.random.uniform(-1000, 1000, num_antennas)
  north = np.random.uniform(-1000, 1000, num_antennas)
  height = np.random.uniform(0, 10, num_antennas)

These coordinates can then be used the calculate *u,v,w* coodinates for a given
array location (I'm using the MWA site) and phase-centre.

First of all, for 10 frequency channels (100MHz to 190MHz at 10MHz resolution),
and for 10 time steps (at a 2s resolution), calculate the phase-tracked
measurement equation:

.. math::

    V_{\textrm{phased}} = \exp\left[2\pi i \left(ul + vm + w(n-1) \right) \right]

where the :math:`u,v,w` and :math:`l,m,n` are calculated with a phase centre of RA, Dec =
:math:`40^\circ, -50^\circ`, and I calculate a single :math:`l,m,n` for a source at
RA, Dec = :math:`10^\circ, -15^\circ` (so in this setup, :math:`u,v,w` change with
time, and :math:`l,m,n` are constant).

I also calculate the  measurement equation without phase tracking, where I calculate
:math:`u_{\mathrm{zen}},v_{\mathrm{zen}},w_{\mathrm{zen}}` and
:math:`l_{\mathrm{zen}},m_{\mathrm{zen}},n_{\mathrm{zen}}`, using the zenith of
the instrument as a coordinate system centre, and use the following
equation:

.. math::

    V_{\textrm{unphased}} = \exp\left[2\pi i \left(u_{\mathrm{zen}}l_{\mathrm{zen}} + v_{\mathrm{zen}}m_{\mathrm{zen}} + w_{\mathrm{zen}}n_{\mathrm{zen}} \right) \right]

(in this setup, :math:`u_{\mathrm{zen}},v_{\mathrm{zen}},w_{\mathrm{zen}}`
are constant with time, and :math:`l_{\mathrm{zen}},m_{\mathrm{zen}},n_{\mathrm{zen}}`
change with time).

I then use :math:`V_{\textrm{phased}}` as an input to ``rw.remove_phase_tracking``
along with :math:`w`, and use that to unwrap the phase tracking. I then assert
that the output of ``rw.remove_phase_tracking`` matches :math:`V_{\textrm{unphased}}`.

test_argument_inputs.py
*******************************************************
These tests run ``rw.get_parser``, which runs the command line parser, and
``rw.check_args``, which checks the ``args`` collected by ``rw.get_parser``
are parsed correctly. It also does sanity checks on certain combinations of args
such that we don't feed WODEN arguments that won't work. The following tests are run
with the expected outcomes:

 - ``test_parser_fails``: There are three required arguments, ``--ra0``, ``--dec0``, and ``--cat_filename``. Check the parser errors if missing.
 - ``test_missing_args_without_metafits_fails``: If the user doesn't supply the ``--metafits`` arg, there are a minimum set of arguments that must be input. Check ``rw.check_args`` errors if they are missing
 - ``test_metafits_read_fails``: Check ``rw.check_args`` errors if there is a bad path to a metafits file
 - ``test_read_metafits_succeeds``: Check the correct values are read in from a known metafits file
 - ``test_EDA2_args_work``: Check the correct arguments are selected for an EDA2 beam simulation
 - ``test_GaussBeam_args_work``: Check that Gaussian beam related arguments work as expected. Iteratively check that if arguments with defaults are not given (e.g. ``--gauss_ra_point``) that they are set to their defaults, and if they *are* supplied, that they match the given value.
 - ``test_MWAFEEBeam_args_work``: Checks that the MWA FEE primary beam is handled by ``ra.check_args`` correctly. The function should error out if certain paths to the hdf5 file that holds the spherical harmonic information is missing, and if the delays have been specified incorrectly. Check that things work when the correct arguments are given.


test_read_uvfits_into_pyuvdata.py
*******************************************************
This tests the absolute minimal compliance with `pyuvdata`_. The test calls
``rw.create_uvfits`` with a set of dummy input variables to make a file
called ``unittest_example.uvfits``. It then simply checks that the following
lines don't throw an error:

.. code-block:: python

  from pyuvdata import UVData
  UV = UVData()
  UV.read('unittest_example.uvfits')

This really only tests that the correct keywords and arrays are present in the
output ``unittest_example.uvfits`` to a level that appeases ``pyuvdata``.
The test is setup to skip if the user has not installed ``pyuvdata``.

.. _`pyuvdata`: https://pyuvdata.readthedocs.io/en/latest/index.html
``chunk_sky_model``
=========================
Tests for the functions in ``WODEN/src/chunk_sky_model.c``. These functions
handle splitting the sky model up into chunks that can fit into GPU memory
during simulation. There is a variable (which can be controlled by the user
via ``--chunking_size``) that controls the maximum number of measurement
equations that can be calculated simultaneously, so the some of the functions
below need to know the number of baselines, frequencies etc to calculate
the maximum number of COMPONENTs to put into each chunk.


``test_null_comps.c``
****************************
Tests the functions ``chunk_sky_model::null_point_comps``,
``chunk_sky_model::null_gauss_comps``, and ``chunk_sky_model::null_shapelet_comps``.
These functions set either POINT, GAUSSIAN, or SHAPELET COMPONENT attributes of
the sky model to ``NULL``. Tested here by setting up a dummy populated sky model,
and making sure each function sets the correct attributes to ``NULL``.
Also tests that other COMPONENT attributes are left as they are.


``test_fill_chunk_src_with_pointgauss.c``
***********************************************
Tests ``chunk_sky_model::fill_chunk_src_with_pointgauss``, which calculates
how many POINT and GAUSSIAN COMPONENTs should go into each chunk. Tests by
trying 12 different combinations of number of POINTS, number of GAUSSIANs,
chunk size, number time steps, number baselines, and number of frequencies.
Iterates over the function and checks that each chunked sky model has the
correct number of output COMPONENTs of each type. Checks by setting attributes
of the sky model to the array index from the original full sky model, and checks
that every attribute is being split correctly into the resultant chunked sky
models. Furthermore, as azimuth and zenith angles are stored for all time steps
in the sky model, the dummy sky model is setup with repeating/tiled arrays
to make sure the correct az/za coords are being copied from the full sky
model into the cropped sky models.

``test_fill_chunk_src_with_shapelets.c``
***********************************************
Tests ``chunk_sky_model::fill_chunk_src_with_shapelets``, which handles
chunking the SHAPELET COMPONENTs of the sky model. Works the same as pointgauss
in that it sets sky model attributes to array indexes, to check whether the
correct values are being chunked into the smaller chunked sky models. On top
of that, due to the way SHAPELETs can have multiple basis functions per
COMPONENT, has an extra layer of sky model trickery to have repeating indexes
inside the basis function attributes, which can then be traced and checked
once the chunking has been performed. Tests here vary the number of SHAPELETs,
the number of basis functions per SHAPELET, and the number of time steps.

``test_create_chunked_sky_models.c``
***********************************************
``chunk_sky_model::create_chunked_sky_models`` uses all functions above to take
in a full sky model and created an array of chunked sky models. This test
runs the same testing for both functions above with a sky model containing
various combinations of POINT, GAUSSIAN, and SHAPELET COMPONENTs, as well
as different chunking sizes and time settings.
``fundamental_coords``
=========================
Tests for the functions in ``WODEN/src/fundamental_coords.cu``.

test_lmn_coords.c
*********************************
This runs ``fundamental_coords::test_kern_calc_lmn``, which tests
``fundamental_coords::kern_calc_lmn``, which calculates *l,m,n* coords.

This runs two control tests, both that generate analytically predictable
outcomes. Both set the phase centre to *RA*:math:`_{\textrm{phase}}`, *Dec*:math:`_{\textrm{phase}}` = :math:`0^\circ, 0^\circ`. One
test holds *Dec* = :math:`0^\circ`, and varies *RA*, the other holds
*RA* = :math:`0^\circ`, and varies *Dec*.  Under these settings the following
should be true:

.. list-table:: Outcomes when *Dec* = :math:`0^\circ`
   :widths: 25 25 25 25
   :header-rows: 1

   * - *RA*
     - *l*
     - *m*
     - *n*
   * - :math:`\frac{3\pi}{2}`
     - :math:`-1`
     - :math:`0`
     - :math:`0`
   * - :math:`\frac{5\pi}{3}`
     - :math:`-\frac{\sqrt{3}}{2}`
     - :math:`0`
     - :math:`0.5`
   * - :math:`\frac{7\pi}{4}`
     - :math:`-\frac{\sqrt{2}}{2}`
     - :math:`0`
     - :math:`\frac{\sqrt{2}}{2}`
   * - :math:`\frac{11\pi}{6}`
     - :math:`-0.5`
     - :math:`0`
     - :math:`\frac{\sqrt{3}}{2}`
   * - :math:`0`
     - :math:`0`
     - :math:`0`
     - :math:`1`
   * - :math:`\frac{\pi}{6}`
     - :math:`0.5`
     - :math:`0`
     - :math:`\frac{\sqrt{3}}{2}`
   * - :math:`\frac{\pi}{4}`
     - :math:`\frac{\sqrt{2}}{2}`
     - :math:`0`
     - :math:`\frac{\sqrt{2}}{2}`
   * - :math:`\frac{\pi}{3}`
     - :math:`\frac{\sqrt{3}}{2}`
     - :math:`0`
     - :math:`0.5`
   * - :math:`\frac{\pi}{2}`
     - :math:`1.0`
     - :math:`0`
     - :math:`0`

.. list-table:: Outcomes when *RA* = :math:`0^\circ`
   :widths: 25 25 25 25
   :header-rows: 1

   * - *Dec*
     - *l*
     - *m*
     - *n*
   * - :math:`-\frac{\pi}{2}`
     - :math:`0`
     - :math:`-1`
     - :math:`0`
   * - :math:`-\frac{\pi}{3}`
     - :math:`0`
     - :math:`-\frac{\sqrt{3}}{4}`
     - :math:`0.5`
   * - :math:`-\frac{\pi}{4}`
     - :math:`0`
     - :math:`-\frac{\sqrt{2}}{2}`
     - :math:`\frac{\sqrt{2}}{2}`
   * - :math:`-\frac{\pi}{6}`
     - :math:`0`
     - :math:`-0.5`
     - :math:`\frac{\sqrt{3}}{2}`
   * - :math:`0`
     - :math:`0`
     - :math:`0`
     - :math:`1`
   * - :math:`\frac{\pi}{6}`
     - :math:`0`
     - :math:`0.5`
     - :math:`\frac{\sqrt{3}}{2}`
   * - :math:`\frac{\pi}{4}`
     - :math:`0`
     - :math:`\frac{\sqrt{2}}{2}`
     - :math:`\frac{\sqrt{2}}{2}`
   * - :math:`\frac{\pi}{3}`
     - :math:`0`
     - :math:`\frac{\sqrt{3}}{2}`
     - :math:`0.5`
   * - :math:`\frac{\pi}{2}`
     - :math:`0`
     - :math:`1.0`
     - :math:`0`

The tests ensure that the inputs yield the outputs as expected, and in
the process check that the execution of the kernel yields the correct number
of outputs. Note that this function is entirely 64-bit whether in FLOAT or
DOUBLE compile mode. The absolute tolerance for outputs vs expectation tabled
above is 1e-15 (this function is a good 'un).

test_uvw_coords.c
*********************************
This runs both ``fundamental_coords::test_kern_calc_uvw`` as well as
``fundamental_coords::test_kern_calc_uvw_shapelet``, which in turn test
``fundamental_coords::kern_calc_uvw`,
``fundamental_coords::kern_calc_uvw_shapelet`` respectively.

Both kernels calculate *u,v,w* coords (in wavelenths) in slightly different circumstances.
``kern_calc_uvw`` calculates *u,v,w* coords towards a specified *RA,Dec* phase centre,
for a given set of baseline lengths :math:`X_{\mathrm{diff}}, Y_{\mathrm{diff}}, Z_{\mathrm{diff}}`,
for a number of LSTs and frequencies (meaning *u,v,w* change with time and frequency).
``kern_calc_uvw_shapelet`` does the above for a number of *u,v,w* coordinate systems,
each centred on a different SHAPELET component.

Both kernels are tested for scaling by wavelength, and changing by time. To generate
analytically predictable outcomes, the phase centre is again set to
*RA*:math:`_{\textrm{phase}}`, *Dec*:math:`_{\textrm{phase}}` = :math:`0^\circ, 0^\circ`.
Under these conditions, the following is true:

.. math::

   \begin{eqnarray}
   u & = & \left[\sin(H_{\textrm{phase}}) X_{\mathrm{diff}} + \cos(H_{\textrm{phase}}) Y_{\mathrm{diff}} \right] / \lambda \\
   v & = & Z_{\mathrm{diff}} / \lambda \\
   w & = & \left[\cos(H_{\textrm{phase}}) X_{\mathrm{diff}} - \sin(H_{\textrm{phase}}) Y_{\mathrm{diff}} \right] / \lambda
   \end{eqnarray}

where :math:`H_{\textrm{phase}}` is the hour angle of the phase centre. These tests
check that this holds true over multiple time and frequency steps. In the case
of ``kern_calc_uvw_shapelet``, this is checked for each SHAPELET component,
making sure that the outputs are ordered as expected. Both the FLOAT and DOUBLE
versions are tested within a tolerance of 1e-16 (this test compares the ``CUDA``
code to the ``C`` code calculation of the above equations, so they agree
very nicely).
.. _mwa_hyperbeam: https://pypi.org/project/mwa-hyperbeam/

.. _FEE_primary_beam_cuda_cmake:

``FEE_primary_beam_cuda``
===========================
Tests for the functions in ``WODEN/src/FEE_primary_beam_cuda.cu``. This is more
of an integration test rather than a suite of individual function tests.

test_RTS_FEE_beam.c
*********************************
This runs ``create_sky_model::test_RTS_CUDA_FEE_beam``, which in turn calls
the following functions:

 - ``create_sky_model::get_HDFBeam_normalisation`` - get values to normalise to zenith
 - ``create_sky_model::copy_FEE_primary_beam_to_GPU`` - move values from host to device
 - ``create_sky_model::calc_CUDA_FEE_beam`` - calculate the MWA FEE beam response
 - ``create_sky_model::free_FEE_primary_beam_from_GPU`` - free things from the device

after copying the MWA FEE beam values from the device back to the host for testing.

The MWA beam pointing direction on the sky is controlled by a set of 16 delays.
In these tests, three different delays settings are tested at 50MHz, 150MHz, and
250MHz (a total of nine tests). Each test is run with ~5000 sky directions,
spanning the whole sky. For each combination of settings, the beam gains
output by ``test_RTS_CUDA_FEE_beam`` are compared to those stored in the header
``test_RTS_FEE_beam.h``.

That header ``test_RTS_FEE_beam.h`` is stitched together from values stored
in text files like ``hyperbeam_zenith_200.txt`` and ``hyperbeam_zenith_200_rot.txt``,
which are generated using the script ``compare_to_hyperdrive.py`` using the
python package `mwa_hyperbeam`_. Each text file stores the az/za, real and imaginary
values for the gain and leakage for both the north-south and east-west dipoles,
either with the parallactic angle rotation applied or not.

All delay and frequency combinations are run with both parallactic angle rotation
applied and not. For the FLOAT precision, the real and imaginary must match
the ``hyperbeam`` values to within a absolute tolerance of 3e-2. For DOUBLE,
they must match to within 1e-13.

.. note:: Given that the accuracy of the FLOAT precision is <= 3%, I suggest if you are comparing simulated visibilities to real data, that DOUBLE is the only way to go (which is accurate to <= 0.00000000001%). However, if you are running some comparison simulations, the inaccuracy is a constant bias, and is still a good representation of the beam model, so is probably fine for internal comparison.

To convince yourself sensible values are stored in those test files, a very rough
plotting script is included in as ``WODEN/cmake_testing/FEE_primary_beam_cuda/plot_beam_results.py``,
which converts the beam gains and leakages into Stokes
XX and YY polarisations, assuming a fully Stokes I sky. The script can be used
as::

  python plot_beam_results.py hyperbeam_zenith_200_rot.txt

which will produce a plot like the below (this is log10(gain) for a parallactic
rotated zenith pointing at 200 MHz).

.. image:: hyperbeam_zenith_200_rot.png
  :width: 400

(Plot looks a little warped purely because I've just done a scatter plot which
is a quick and dirty way of showing the info).

If you are *really* interested in the differences, you can run::

  $ source plot_all_beam_diffs.sh

which will produce a bunch of plots in a directory
``WODEN/cmake_testing/FEE_primary_beam_cuda/beam_plots``. Included are difference
plots, showing the offset from the ``WODEN`` output to ``hyperbeam`` for real
and imaginary in the gain and leakage terms of the Jones matrix. An example at
FLOAT precision for zenith at 100MHz is (note this is the difference
in gain, NOT log10(gain)):

.. image:: hyperbeam_zenith_100_rotzenith_100_rot_float_diff.png
  :width: 400


The equivalent plot for DOUBLE is shown below, showing the vast improvement in
accuracy:

.. image:: hyperbeam_zenith_100_rotzenith_100_rot_double_diff.png
  :width: 400
``array_layout``
=========================
Tests for the functions in ``WODEN/src/array_layout.c``. These functions handle
reading in array coords of *e,n,h* and converting them into *X,Y,Z*, used later
to calculate *u,v,w*. Also handles precession of the array layout and LST from the
simulation date back to J2000, so the array is in the same coord frame as the
sky model. Should be the equivalent to rotating the sky model forward to current
day and ensures the output visibilities are in the J2000 frame.

.. _test_RTS_ENH2XYZ_local.c:

``test_RTS_ENH2XYZ_local.c``
*****************************
``array_layout::RTS_ENH2XYZ_local`` transforms a local east, north, height coord
into *X,Y,Z* coords, which can be used to calculate *u,v,w*. Tests here
generate some test *e,n,h* coords, and then tests for two latitude test cases:

 - latitude = :math:`0^\circ`, should result in:
    - :math:`X = h`
    - :math:`Y = e`
    - :math:`Z = n`
 - latitude = :math:`-30^\circ`, should result in:
    - :math:`X = \frac{n}{2} + \frac{\sqrt(3)}{2}h`
    - :math:`Y = e`
    - :math:`Z = \frac{\sqrt(3)}{2}n - \frac{h}{2}`

Both the FLOAT and DOUBLE compiled versions are 64 bit precision. Both are
tested to match a 64 bit ``C`` calculation of the above equations to within an
absolute tolerance of 1e-13.

``test_RTS_PrecessXYZtoJ2000.c``
*********************************
``array_layout::RTS_PrecessXYZtoJ2000`` takes observation date ``X,Y,Z`` coords
and precesses them back to J2000, and updates the LST accordingly.
Test by giving a known set of *X,Y,Z* coordinates and julian date, and
checks output precessed *X,Y,Z* and LST match the expected values stored in
``test_RTS_XYZ_common.h``.

``test_calc_XYX_diffs.c``
****************************
Tests the function ``array_layout::calc_XYZ_diffs``, which reads in an array
text file in *e,n,h* coords and transforms to *X,Y,Z*, precesses the locations
back to J2000 if requested, and then calculates the baseline length in *X,Y,Z*.
It also rotates the latitude of the array back to J2000 and in doing so
updates the LST.

These tests simply provide a set of known input coords (read in from
``example_array_layout.txt``), runs ``array_layout::calc_XYZ_diffs`` for both
the precession and no precession cases, and check the output values match the
expected values stored in ``test_RTS_XYZ_common.h``.
``woden_settings``
=========================
Tests for the functions in ``WODEN/src/woden_settings.c``. These functions handle
reading input simulation settings from a ``.json`` file, and filling a
``woden_settings_t`` struct based on them. ``woden_settings_t`` is passed
around by most functions internally to ``WODEN`` to tell them what to do.

``test_read_json_settings.c``
******************************
Tests the function ``woden_settings::read_json_settings``. This function
reads in the simulation settings held in a ``.json`` file (usually written by
``run_woden.py``). This test runs reads multiple ``.json`` files in, each
with different settings. All test ``.json`` files share these common core
settings:

.. code-block:: json

   {
   "ra0": 0.0000000000,
   "dec0": -27.0000000000,
   "num_freqs": 16,
   "num_time_steps": 4,
   "cat_filename": "srclist_singlepoint.txt",
   "time_res": 2.0,
   "frequency_resolution": 40000.0,
   "chunking_size": 5000,
   "jd_date": 2457278.2010995,
   "LST": 0.44312771,
   "array_layout": "example_array_layout.txt",
   "lowest_channel_freq": 1.6703500000e+08,
   "latitude": -26.70331944,
   "coarse_band_width": 1.2800000000e+06,
   "sky_crop_components": "True",
   "band_nums": [1,4,9]
   }

For all tests, these core attributes are tested as being read in correctly into
``woden_settings_t`` after calling ``read_json_settings``. The test
files listed in this table are all run to test different simulation setups and
eventualities.

.. list-table::
   :widths: 25 50
   :header-rows: 1

   * - run_woden_nobeam.json
     - The NO_BEAM primary beam is selected
   * - run_woden_EDA2.json
     - The ANALY_DIPOLE primary beam is selected
   * - run_woden_gaussian_bespoke.json
     - The GAUSSIAN primary beam is selected, using input values to set the FWHM and reference frequency
   * - run_woden_gaussian_default.json
     - The GAUSSIAN primary beam is selected, using default values for the FWHM and reference frequency
   * - run_woden_MWAFEE.json
     - The MWA_FEE primary beam is selected, and associated delays and path to hdf5 file are read in correctly
   * - run_woden_MWAFEE_baddelay.json
     - Contains a bad set of MWA FEE delays and should throw an error
   * - run_woden_MWAFEE_nopath.json
     - Has no path to the MWA FEE hdf5 file so should throw an error
   * - run_woden_multiple_beams.json
     - Contains multiple primary beam selections so should throw an error
   * - run_woden_noprecession.json
     - Checks that precession is switched off if requested

``test_setup_lsts_and_phase_centre.c``
*****************************************
Tests ``woden_settings::setup_lsts_and_phase_centre`` which uses the input
simulation parameters to calculate the sine and cosine of the declination of
the phase centre, and the LST at the centre of every time integration. Runs
three tests, two with :math:`\mathrm{Dec}_{\mathrm{phase}} = \pi/2`, where one has a
zero initial LST, another with non-zero initial LST, and third test with
:math:`\mathrm{Dec}_{\mathrm{phase}}` set to the latitude of the MWA and a
non-zero initial LST. All tests have multiple time steps, and are tested against
equivalent calculations made in ``C`` with 64 bit precision. The FLOAT code
outputs are tested with an absolute tolerance of 1e-7 and the DOUBLE a tolerance
of 1e-15.
``source_components``
=========================
Tests for the functions in ``WODEN/src/source_components.cu``. These functions
calculate visibilities for the different COMPONENT types, as well as calling
the various beam models to be applied to the visibilities.


test_apply_beam_gains.c
************************************
This calls ``source_components::test_kern_apply_beam_gains``, which tests
``source_components::kern_apply_beam_gains``. This kernel applies
beam gain and leakage terms to Stokes visibilities to create linear Stokes
polarisation visibilities via:

.. math::

   \begin{eqnarray}
   \mathrm{V}^{XX}_{12} = (g_{1x}g_{2x}^{\ast} + D_{1x}D_{2x}^{\ast})\mathrm{V}^{I}_{12}
    +  (g_{1x}g_{2x}^{\ast} - D_{1x}D_{2x}^{\ast})\mathrm{V}^{Q}_{12} \\
    +  (g_{1x}D_{2x}^{\ast} + D_{1x}g_{2x}^{\ast})\mathrm{V}^{U}_{12}
    +  i(g_{1x}D_{2x}^{\ast} - D_{1x}g_{2x}^{\ast})\mathrm{V}^{V}_{12}
   \end{eqnarray}

.. math::

   \begin{eqnarray}
   \mathrm{V}^{XY}_{12} =
        (g_{1x}D_{2y}^{\ast} + D_{1x}g_{2y}^{\ast})\mathrm{V}^{I}_{12}
     +  (g_{1x}D_{2y}^{\ast} - D_{1x}g_{2y}^{\ast})\mathrm{V}^{Q}_{12} \\
     +  (g_{1x}g_{2y}^{\ast} + D_{1x}D_{2y}^{\ast})\mathrm{V}^{U}_{12}
     +  i(g_{1x}g_{2y}^{\ast} - D_{1x}D_{2y}^{\ast})\mathrm{V}^{V}_{12}
   \end{eqnarray}

.. math::

   \begin{eqnarray}
   \mathrm{V}^{XY}_{12} =
        (g_{1x}D_{2y}^{\ast} + D_{1x}g_{2y}^{\ast})\mathrm{V}^{I}_{12}
     +  (g_{1x}D_{2y}^{\ast} - D_{1x}g_{2y}^{\ast})\mathrm{V}^{Q}_{12} \\
     +  (g_{1x}g_{2y}^{\ast} + D_{1x}D_{2y}^{\ast})\mathrm{V}^{U}_{12}
     +  i(g_{1x}g_{2y}^{\ast} - D_{1x}D_{2y}^{\ast})\mathrm{V}^{V}_{12}
   \end{eqnarray}

.. math::

   \begin{eqnarray}
   \mathrm{V}^{YY}_{12} =
        (D_{1y}D_{2y}^{\ast} + g_{1y}g_{2y}^{\ast})\mathrm{V}^{I}_{12}
     +  (D_{1y}D_{2y}^{\ast} - g_{1y}g_{2y}^{\ast})\mathrm{V}^{Q}_{12} \\
     +  (D_{1y}g_{2y}^{\ast} + g_{1y}D_{2y}^{\ast})\mathrm{V}^{U}_{12}
     +  i(D_{1y}g_{2y}^{\ast} - g_{1y}D_{2y}^{\ast})\mathrm{V}^{V}_{12}
   \end{eqnarray}

where :math:`\ast` means complex conjugate, :math:`g_x, D_x, D_y, g_y` are beam
gain and leakage terms, with subscript 1 and 2 meaning antenna 1 and 2, and
:math:`\mathrm{V}^I, \mathrm{V}^Q, \mathrm{V}^U, \mathrm{V}^V` the Stokes
visibilities. These tests try a number of combinations of values that have
a simple outcome, and tests that the function returns the expected values. The
combinations are shown in the table below. For all combinations, the beam gain
and leakage is used for both antennas in the above equations. Each entry is a
complex values and should be read as *real,imag*. Both FLOAT and DOUBLE code
are just tested to a 32 bit accuracy here as these a simple numbers that require
little precision.

.. list-table::
   :widths: 25 25 25 25 25 25 25 25 25 25 25 25
   :header-rows: 1

   * - :math:`g_x`
     - :math:`D_x`
     - :math:`D_y`
     - :math:`g_y`
     - :math:`\mathrm{V}^I`
     - :math:`\mathrm{V}^Q`
     - :math:`\mathrm{V}^U`
     - :math:`\mathrm{V}^V`
     - :math:`\mathrm{V}^{XX}`
     - :math:`\mathrm{V}^{XY}`
     - :math:`\mathrm{V}^{YX}`
     - :math:`\mathrm{V}^{YY}`
   * - 1,0
     - 0,0
     - 0,0
     - 1,0
     - 1,0
     - 0,0
     - 0,0
     - 0,0
     - 1,0
     - 0,0
     - 0,0
     - 1,0
   * - 1,0
     - 0,0
     - 0,0
     - 1,0
     - 0,0
     - 1,0
     - 0,0
     - 0,0
     - 1,0
     - 0,0
     - 0,0
     - -1,0
   * - 1,0
     - 0,0
     - 0,0
     - 1,0
     - 0,0
     - 0,0
     - 1,0
     - 0,0
     - 0,0
     - 1,0
     - 1,0
     - 0,0
   * - 1,0
     - 0,0
     - 0,0
     - 1,0
     - 0,0
     - 0,0
     - 0,0
     - 1,0
     - 0,0
     - 0,1
     - 0,-1
     - 0,0
   * - 0,0
     - 1,0
     - 1,0
     - 0,0
     - 1,0
     - 0,0
     - 0,0
     - 0,0
     - 1,0
     - 0,0
     - 0,0
     - 1,0
   * - 0,0
     - 1,0
     - 1,0
     - 0,0
     - 0,0
     - 1,0
     - 0,0
     - 0,0
     - -1,0
     - 0,0
     - 0,0
     - 1,0
   * - 0,0
     - 1,0
     - 1,0
     - 0,0
     - 0,0
     - 0,0
     - 1,0
     - 0,0
     - 0,0
     - 1,0
     - 1,0
     - 0,0
   * - 0,0
     - 1,0
     - 1,0
     - 0,0
     - 0,0
     - 0,0
     - 0,0
     - 1,0
     - 0,0
     - 0,-1
     - 0,1
     - 0,0
   * - 2,0
     - 2,0
     - 2,0
     - 2,0
     - 1,0
     - 0,0
     - 0,0
     - 0,0
     - 8,0
     - 8,0
     - 8,0
     - 8,0
   * - 2,0
     - 2,0
     - 2,0
     - 2,0
     - 0,0
     - 1,0
     - 0,0
     - 0,0
     - 0,0
     - 0,0
     - 0,0
     - 0,0
   * - 2,0
     - 2,0
     - 2,0
     - 2,0
     - 0,0
     - 0,0
     - 1,0
     - 0,0
     - 8,0
     - 8,0
     - 8,0
     - 8,0
   * - 2,0
     - 2,0
     - 2,0
     - 2,0
     - 0,0
     - 0,0
     - 0,0
     - 1,0
     - 0,0
     - 0,0
     - 0,0
     - 0,0
   * - 1,2
     - 3,4
     - 5,6
     - 7,8
     - 1,0
     - 0,0
     - 0,0
     - 0,0
     - 30,0
     - 70,8
     - 70,-8
     - 174,0
   * - 1,2
     - 3,4
     - 5,6
     - 7,8
     - 0,0
     - 1,0
     - 0,0
     - 0,0
     - -20,0
     - -36,0
     - -36,0
     - -52,0
   * - 1,2
     - 3,4
     - 5,6
     - 7,8
     - 0,0
     - 0,0
     - 1,0
     - 0,0
     - 22,0
     - 62,8
     - 62,-8
     - 166,0
   * - 1,2
     - 3,4
     - 5,6
     - 7,8
     - 0,0
     - 0,0
     - 0,0
     - 1,0
     - -4,0
     - -4,-16
     - -4,16
     - -4,0

test_calc_measurement_equation.c
************************************
This calls ``source_components::test_kern_calc_measurement_equation``, which
calls ``source_components::kern_calc_measurement_equation``, which in turn
is testing the device code ``source_components::calc_measurement_equation``
(which is used internally in ``WODEN`` in another kernel that calls multiple
device functions, so I had to write a new kernel to test it alone.)

The following methodology is also written up the JOSS paper (TODO link JOSS
paper once it's accepted), where it's used in a sky model + array layout through
to visibility end-to-end version. Here we directly test the fuction above.

``calc_measurement_equation`` calculates the phase-tracking measurement equation:

.. math::

  V(u,v,w) =  \exp \left[ 2\pi i\left( ul + vm + w(n-1) \right) \right]

We can use Euler's formula to split this into real and imaginary components. If
I label the phase for a particular source and baseline as

.. math::

  \phi = 2\pi \left( ul + vm + w(n-1)\right)

then the real and imaginary parts of the visibility :math:`V_{re}`, :math:`V_{im}` are

.. math::

  V_{re} = \cos(\phi) \\
  V_{im} = \sin(\phi)

A definitive test of the ``calc_measurement_equation`` function then is to then set
:math:`\phi` to a number of values which produce known sine and cosine outputs, by
selecting specific combinations of *u,v,w* and *l,m,n*. First of all, consider the case when
*u,v,w = 1,1,1*. In that case,

.. math::

  \frac{\phi_{\mathrm{simple}}}{2\pi} = l + m + (n-1).

if we further set *l == m*, we end up with

.. math::

  \frac{\phi_{\mathrm{simple}}}{2\pi} = 2l + (n-1), \\
  l = \sqrt{\left( \frac{1 - n^2}{2} \right)}

I shoved those two equations into `Wolfram Alpha`_ who assured me that a solution
for *n* here is

.. _Wolfram Alpha: https://www.wolframalpha.com/widgets/view.jsp?id=c07cc70f1e81887dfd0971d3fe17cfcd

.. math::

  n = \frac{\sqrt{2}\sqrt{-\phi_{\mathrm{simple}}^2 - 4\pi\phi_{\mathrm{simple}} + 8\pi^2} + \phi_{\mathrm{simple}} + 2\pi}{6\pi}

which we can then use to calculate values for *l,m* through

.. math::

  l = m = \sqrt{\frac{1 - n^2}{2}}.

By selecting the following values for :math:`\phi`, we can create the following
set of *l,m,n* coords, which have the a known set of outcomes:

.. list-table::
   :widths: 25 25 25 25 25
   :header-rows: 1

   * - :math:`\phi_{\mathrm{simple}}`
     - *l,m*
     - *n*
     - :math:`\cos(\phi)`
     - :math:`\sin(\phi)`
   * - :math:`0`
     - 0.0
     - 1.0
     - :math:`1.0`
     - :math:`0`
   * - :math:`\pi/6`
     - 0.0425737516338956
     - 0.9981858300655398
     - :math:`\sqrt{3}/2`
     - :math:`0.5`
   * - :math:`\pi/4`
     - 0.0645903244635131
     - 0.9958193510729726
     - :math:`\sqrt{2}/2`
     - :math:`\sqrt{2}/2`
   * - :math:`\pi/3`
     - 0.0871449863555500
     - 0.9923766939555675
     - :math:`0.5`
     - :math:`\sqrt{3}/2`
   * - :math:`\pi/2`
     - 0.1340695840364469
     - 0.9818608319271057
     - :math:`0.0`
     - :math:`1.0`
   * - :math:`2\pi/3`
     - 0.1838657911209207
     - 0.9656017510914922
     - :math:`-0.5`
     - :math:`\sqrt{3}/2`
   * - :math:`3\pi/4`
     - 0.2100755148372292
     - 0.9548489703255412
     - :math:`-\sqrt{2}/2`
     - :math:`\sqrt{2}/2`
   * - :math:`5\pi/6`
     - 0.2373397982598921
     - 0.9419870701468823
     - :math:`-\sqrt{3}/2`
     - :math:`0.5`
   * - :math:`\pi`
     - 0.2958758547680685
     - 0.9082482904638630
     - :math:`-1.0`
     - :math:`0.0`
   * - :math:`7\pi/6`
     - 0.3622725654470420
     - 0.8587882024392495
     - :math:`-\sqrt{3}/2`
     - :math:`-0.5`
   * - :math:`5\pi/4`
     - 0.4003681253515569
     - 0.8242637492968862
     - :math:`-\sqrt{2}/2`
     - :math:`-\sqrt{2}/2`

.. note:: If you try and go higher in :math:`\phi` then because I set :math:`l == m` you no longer honour :math:`\sqrt{l^2 + m^2 + n^2} <= 1.0` I think this range of angles is good enough coverage though.

This is a great test for when :math:`u,v,w = 1`, but we want to test a range of
baseline lengths to check our function is consistent for short and long baselines.
We can play another trick, and set all baseline coords to be equal, i.e. :math:`u = v = w = b` where :math:`b` is baseline length. In this form, the phase including
the baseline length :math:`\phi_{b}` is

.. math::

  \phi_{b} = 2\pi b\left( l + m + n - 1 \right) = b\phi_{\mathrm{simple}}.

As sine/cosine are periodic functions, the following is true:

.. math::

  \phi_{\mathrm{simple}} = \phi_{\mathrm{simple}} + 2\pi \mathrm{n}

where :math:`\mathrm{n}` is some integer. This means for a given :math:`\phi_{\mathrm{simple}}` from
the table above, we can find an appropriate :math:`b` that should still result in the
expected sine and cosine outputs by setting

.. math::

  b\phi_{\mathrm{simple}} = \phi_{\mathrm{simple}} + 2\pi \mathrm{n} \\
  b = \frac{\phi_{\mathrm{simple}} + 2\pi \mathrm{n}}{\phi_{\mathrm{simple}}}

for a range of :math:`\mathrm{n}` values. The values of :math:`\mathrm{n}` and the
resultant size of :math:`b` that I use in testing are shown in the table below (note for :math:`\phi_{\mathrm{simple}} = 0` I just set :math:`b = 2\pi \mathrm{n}` as the effects of :math:`l,m,n` should set everything to zero regardless of baseline coords).

.. list-table::
   :widths: 25 25 25 25 25 25 25
   :header-rows: 1

   * - :math:`\phi_{\mathrm{simple}}`
     - :math:`b(\mathrm{n=0})`
     - :math:`b(\mathrm{n=1})`
     - :math:`b(\mathrm{n=10})`
     - :math:`b(\mathrm{n=100})`
     - :math:`b(\mathrm{n=1000})`
     - :math:`b(\mathrm{n=10000})`
   * - :math:`0`
     - 0.0
     - 6.3
     - 62.8
     - 628.3
     - 6283.2
     - 62831.9
   * - :math:`\pi/6`
     - 1.0
     - 13.0
     - 121.0
     - 1201.0
     - 12001.0
     - 120001.0
   * - :math:`\pi/4`
     - 1.0
     - 9.0
     - 81.0
     - 801.0
     - 8001.0
     - 80001.0
   * - :math:`\pi/3`
     - 1.0
     - 7.0
     - 61.0
     - 601.0
     - 6001.0
     - 60001.0
   * - :math:`\pi/2`
     - 1.0
     - 5.0
     - 41.0
     - 401.0
     - 4001.0
     - 40001.0
   * - :math:`2\pi/3`
     - 1.0
     - 4.0
     - 31.0
     - 301.0
     - 3001.0
     - 30001.0
   * - :math:`3\pi/4`
     - 1.0
     - 3.7
     - 27.7
     - 267.7
     - 2667.7
     - 26667.7
   * - :math:`5\pi/6`
     - 1.0
     - 3.4
     - 25.0
     - 241.0
     - 2401.0
     - 24001.0
   * - :math:`\pi`
     - 1.0
     - 3.0
     - 21.0
     - 201.0
     - 2001.0
     - 20001.0
   * - :math:`7\pi/6`
     - 1.0
     - 2.7
     - 18.1
     - 172.4
     - 1715.3
     - 17143.9
   * - :math:`5\pi/4`
     - 1.0
     - 2.6
     - 17.0
     - 161.0
     - 1601.0
     - 16001.0

This gives a range of baseline lengths from 1 to :math:`> 10^4` wavelengths.

In this test, I run every combination of :math:`l,m,n` and :math:`u,v,w = b` for each
:math:`\phi_{\mathrm{simple}}` from the tables above, and assert that the real and
imaginary of every output visibility match the expected values of
:math:`\sin(\phi_{\mathrm{simple}})` and :math:`\cos(\phi_{\mathrm{simple}})`.
When compiling with FLOAT precision, I assert the outputs must be within
an absolute tolerance of 2e-3, and for DOUBLE a tolerance of 2e-9.

The error scales with the length of baseline, as shown in this plot below. Here,
I have plotted the fractional offset of the recovered value of
:math:`\sin(\phi_{\mathrm{simple}})` (imaginary part of the visibility) and
:math:`\cos(\phi_{\mathrm{simple}})` (real part of the visibility), compared
to their analytically expected outcome. I've plotted each :math:`\phi_{\mathrm{simple}}`
as a different symbol, with the FLOAT in blue and DOUBLE in orange.

.. image:: measure_eq_results.png
  :width: 800

You can see as you increase baseline length, a general trend of increasing error
is seen. Note these are the absolute differences plotted here, to work on a log10
scale. For the DOUBLE results, it's close to 50% a negative or positive offset
from expected. The FLOAT results that perform the worst (where
:math:`\phi_{\mathrm{simple}} = 3\pi/4,\, 7\pi/6`) correspond to values of :math:`b`
with large fractional values (see the table above), showing how the 32 bit
precision fails to truthfully report large fractional numbers.

As a second test, I setup 10,000 *u,v,w* ranging from -1000 to 1000 wavelengths,
and 3600 *l,m,n* coordinates that span the entire sky, and run them through
``source_components::test_kern_calc_measurement_equation`` and check they
equal the equivalent measurement equation as calculated by ``C`` in 64 bit precisions.
This checks the kernel works for a range of input coordinates.
I assert the ``CUDA`` outputs must be within an absolute tolerance of 1e-7 for
the FLOAT code, and 1e-15 for the DOUBLE code.

   .. TODO One day, could add in estimation of effect on calibration

test_extrap_stokes.c
************************************
This calls ``source_components::test_kern_extrap_stokes``, which
calls ``source_components::kern_extrap_stokes``, which handles extrapolating
a reference flux density to a number of frequencies, given a spectral index.

Five test cases are used, with the following parameters:

.. list-table::
   :widths: 25 25 25 25 25 25
   :header-rows: 1

   * - Reference Freq (MHz)
     - Spectral Index
     - Stokes *I*
     - Stokes *Q*
     - Stokes *U*
     - Stokes *V*
   * - 50
     - 0.0
     - 1.0
     - 0.0
     - 0.0
     - 0.0
   * - 100
     - -0.8
     - 1.0
     - 0.0
     - 0.0
     - 0.0
   * - 150
     - 0.5
     - 1.0
     - 1.0
     - 0.0
     - 0.0
   * - 200
     - -0.5
     - 1.0
     - 0.0
     - 1.0
     - 0.0
   * - 250
     - 1.0
     - 1.0
     - 0.0
     - 0.0
     - 1.0

Each of these test cases is extrapolated to 50, 100, 150, 200, 250 MHz. The
``CUDA`` outputs are testing as being equal to

.. math::

   S_{\mathrm{extrap}} = S_{\mathrm{ref}} \left( \frac{\nu_{\mathrm{ref}}}{\nu_{\mathrm{extrap}}} \right)^{\alpha}

as calculated in 64 bit precision in ``C`` code in ``test_extrap_stokes.c``.
The FLOAT complied code must match the ``C`` estimate to within an absolute
tolerance of 1e-7, and a tolerance of 1e-15 for the DOUBLE compiled code.

test_get_beam_gains.c
************************************
This calls ``source_components::test_kern_get_beam_gains``, which
calls ``source_components::kern_get_beam_gains``, which in turn is testing the
device code ``source_components::get_beam_gains``. This function handles grabbing
the pre-calculated beam gains for a specific beam model, time, and frequency
(assuming the beam gains have already been calculated). ``kern_get_beam_gains``
is setup to call ``get_beam_gains`` for multiple inputs and recover them into a
set of output arrays.

Beam gain calculations are stored in ``primay_beam_J*`` arrays, including
all frequency and time steps, as well as all directions on the sky. This test
sets all real entries in the four ``primay_beam_J*`` beam gain arrays to the
value of their index. In this way, we can easily predict the expected value
in the outputs as being the index of the beam gain we wanted to select. I've
set the imaginary to zero.

This test runs with two time steps, two frequency steps, three baselines,
and four beam models. Three different outcomes are expected given the beam model:

 - ANALY_DIPOLE, GAUSS_BEAM: The values of the gains are testing to match the expected index. The leakage terms are tested to be zero as the models have no leakage terms
 - FEE_BEAM: Both the beam gain and leakage terms are tested as this model includes leakage terms
 - NO_BEAM: The gain terms are tested to be 1.0, and leakage to be 0.0

Both FLOAT and DOUBLE code are tested to a 32 bit accuracy here as these are
simple numbers that require little precision.

test_source_component_common.c
************************************
This calls ``source_components::test_source_component_common``, which
calls ``source_components::source_component_common``. ``source_component_common``
is run by all visibility calculation functions (the functions
``kern_calc_visi_point``, ``kern_calc_visi_gauss``, ``kern_calc_visi_shape``).
It handles calculating the *l,m,n* coordinates and beam response for all
COMPONENTs in a sky model, regardless of the type of COMPONENT.

Similarly to the tests in :ref:`test_lmn_coords.c`, I setup a slice of 9 *RA*
coordinates, and hold the *Dec* constant, set the phase centre to
*RA*:math:`_{\textrm{phase}}`, *Dec*:math:`_{\textrm{phase}}` = :math:`0^\circ, 0^\circ`.
This way I can analytically predict what the *l,m,n* calculated coordinates
should be (which are tested to be within 1e-15 of expected values).

In these tests I run with three time steps, two frequency steps (100 and 200 MHz),
and five baselines (the coordinates of which don't matter, but change the size
of the outputs, so good to have a non-one value). I input a set of *az,za* coords
that match the *RA,Dec* coords for an *LST* = 0.0. As I have other tests that check
the sky moves with time, I just set the sky to be stationary with time here, to
keep the test clean.

For each primary beam type, I run the 9 COMPONENTs through the test, and check
the calcualted *l,m,n* are correct, and check that the calculated beam values
match a set of expected values, which are stored in ``test_source_component_common.c``. As with previous tests varying the primary beam, I check that leakage terms
should be zero when the model doesn't include them.

The absolute tolerance values used for the different beam models, for the two
different precisions are shown in the table below. Note I've only stored the
expected values for the ANALY_DIPOLE and FEE_BEAM to 1e-7 accuracy, as the
accuracy of these functions beam functions is tested elsewhere. The
GAUSS_BEAM values are calculated analytically in the same test as
described in :ref:`test_gaussian_beam.c`.

.. list-table::
   :widths: 25 25 25
   :header-rows: 1

   * - Beam type
     - FLOAT tolerance
     - DOUBLE tolerance
   * - GAUSS_BEAM
     - 1e-7
     - 1e-12
   * - ANALY_DIPOLE
     - 1e-6
     - 1e-7
   * - FEE_BEAM
     - 3e-2
     - 1e-7

test_kern_calc_visi_point.c
************************************
This calls ``source_components::test_kern_calc_visi_point``, which
calls ``source_components::kern_calc_visi_point``. This kernel calculates
the visibility response for POINT COMPONENTs for a number of sky directions, for
all time and frequency steps, and all baselines.

I set up a grid of 25 *l,m* coords with *l,m* coords ranging over -0.5, -0.25,
0.0, 0.25, 0.5. I run a simulation with 10 baselines, where I set *u,v,w* to:

.. math::

   u,v = 100(b_{\mathrm{ind}} + 1) \\
   w = 10(b_{\mathrm{ind}} + 1)

where :math:`b_{\mathrm{ind}}` is the baseline index, meaning the test covers
the baseline range :math:`100 < u,v <= 1000` and :math:`10 < w <= 100`. The test
also runs three frequncies, 150, 175, 200 MHz, and two time steps. As I am providing
predefined *u,v,w*, I don't need to worry about LST effects, but I simulate with
two time steps to make sure the resultant visibilities end up in the right order.

Overall, I run three groups of tests here:

 - Keeping the beam gains and flux densities constant at 1.0
 - Varying the flux densities with frequency and keeping the beam gains constant at 1.0. When varying the flux, I set the Stokes I flux of each component to it's index + 1, so we end up with a range of fluxes between 1 and 25. I set the spectral index to -0.8.
 - Varying the beam gains with frequency and keeping the flux densities constant at 1.0. As the beam can vary with time, frequency, and direction on sky, I assign each beam gain a different value. As *num_freqs*num_times*num_components* = 375, I set the real of all beam gains to :math:`\frac{1}{375}(B_{\mathrm{ind}} + 1)`, where :math:`B_{\mathrm{ind}}` is the beam value index. This way we get a unique value between 0 and 1 for all beam gains, allowing us to test time/frequency is handled correctly by the function under test

Each set of tests is run for all four primary beam types, so a total of 12 tests
are called. Each test calls ``kern_calc_visi_point``, which should calculate
the measurement equation for all baselines, time steps, frequency steps, and COMPONENTs.
It should also sum over COMPONENTs to get the resultant visibility for each
baseline, time, and freq. To test the outputs, I have created equivalent ``C``
functions at 64 bit precisions in ``test_kern_calc_visi_common.c`` to calculate
the measurement equation for the given inputs. For all visibilities, for the FLOAT version
I assert the ``CUDA`` code output must match the ``C`` code output to
within an fractional tolerance of 1e-5 to the ``C`` value, for both the real and
imaginary parts. For the DOUBLE code, the fractional tolerance is 1e-13. I've
switched to fractional tolerance here as the range of magnitudes covered by
these visibilities means a small absolute tolernace will fail a large magnitude
visibility when it reports a value that is correct to 1e-11%.

test_kern_calc_visi_gauss.c
************************************
This calls ``source_components::test_kern_calc_visi_gauss``, which
calls ``source_components::kern_calc_visi_gauss``. This kernel calculates
the visibility response for GAUSSIAN COMPONENTs for a number of sky directions, for
all time and frequency steps, and all baselines.

This runs all tests as described by :ref:`test_kern_calc_visi_point.c`, plus a
fourth set of tests that varies the position angle, major, and minor axis of the
input GAUSSIAN components, for a total of 16 tests. Again, I have ``C`` code to
test the ``CUDA`` code against. I assert the ``CUDA`` code output must match the
``C`` code output to within an fractional tolerance of 1e-5 to the ``C`` value,
for both the real and imaginary parts. For the DOUBLE code, the fractional
tolerance is 1e-13.

test_kern_calc_visi_shape.c
************************************
This calls ``source_components::test_kern_calc_visi_gauss``, which
calls ``source_components::kern_calc_visi_gauss``. This kernel calculates
the visibility response for GAUSSIAN COMPONENTs for a number of sky directions, for
all time and frequency steps, and all baselines.

This runs all tests as described by :ref:`test_kern_calc_visi_gauss.c`, plus a
fifth set of tests that gives multiple shapelet basis function parameters to the
input SHAPELET components, for a total of 20 tests. Again, I have ``C`` code
to test the ``CUDA`` code against.

The final 5th test really pushes the FLOAT code hard, as the range of magnitudes
of the visibilities is large. As a result, FLOAT code is tested to to within a
fractional tolerance of 1e-2 to the ``C`` values (which happens mostly when
the expected value is around 1e-5 Jy, so a fractional offset of 1e-2 is an
absolute offset of 1e-7 Jy), for both the real and imaginary parts.
For the DOUBLE code, the fractional tolerance is 1e-12.

test_update_sum_visis.c
************************************
This calls ``source_components::test_kern_update_sum_visis``, which in turn calls
``source_components::kern_update_sum_visis``. This kernel gathers pre-calculated
primary beam values, unpolarised visibilities, Stokes flux densities, and
combines them all into linear polarisation Stokes visibilities. It then
sums all COMPONENTs together onto the final set of linear polarisation
visibilities.

This code runs three sets of tests, each with three baselines, four time steps, three frequencies, and ten COMPONENTs. For each set of tests, the input
arrays that are being tested have their values set to the index of the value,
making the summation able to test whether the correct time, frequency, and
beam indexes are being summed into the resultant visibilities. The
three sets of tests consist of:

   - Varying the beam gains, while keeping the flux densities and unpolarised measurement equations constant.
   - Varying the flux densities, while keeping the beam gains and unpolarised measurement equations constant.
   - Varying the unpolarised measurement equations, while keeping the beam gains and flux densities constant.

Each set of tests is run for all primary beam types, for a total of 12 tests.
The different beam models have different expected values depending on whether
they include leakage terms or not.
*************
Installation
*************

WODEN is built on CUDA so you will need an NVIDIA GPU to run it. Currently, WODEN has only been tested and run on linux, specifically Ubuntu 16.04 up to 20.04, the OzStar super cluster of Swinburne University, and Garrawarla cluster of Pawsey. If you're mad keen to run on Windows or Mac, please contact Jack at jack.l.b.line@gmail.com and we can give it a go.

Dependencies
##############

``WODEN`` has a number of dependencies so it doesn't reinvent the wheel. A brief list of them here is followed by detailed instructions on how I installed them in the following subsection. Note that the explicit installation instructions I have included for ``json-c``, ``erfa``, and ``pal`` are the only way I have reliably managed to install these packages - the package installation manager sometimes does whacky things for them.

- **CMake** - https://cmake.org version >= 3.10
- **NVIDIA CUDA** - https://developer.nvidia.com/cuda-downloads
- **json-c** - https://github.com/json-c/json-c
- **ERFA** - https://github.com/liberfa/erfa/releases
- **HDF5** - https://www.hdfgroup.org/downloads/hdf5/
- **PAL** - https://github.com/Starlink/pal/releases
- **python >= 3.6**

How to install dependencies
****************************

These instructions are for Ubuntu 20.04, but can be used as a guide for other
linux-like systems.

+ **CMake** - https://cmake.org version >= 3.10::

   $ sudo snap install cmake

+ **NVIDIA CUDA** - https://developer.nvidia.com/cuda-downloads. I typically download the runfile option, which you run as::

  $ sudo sh cuda_11.2.2_460.32.03_linux.run

  but I do NOT install the drivers at this point, as I'll already have drivers. Up to you and how your system works. Also, don't ignore the step of adding something like ``export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-11.2/lib64`` to your ``~/.bashrc``, or your system won't find ``CUDA``.
+ **json-c** - https://github.com/json-c/json-c. This is a typical ``cmake`` installation::

  $ git clone https://github.com/json-c/json-c.git
  $ cd json-c
  $ mkdir build && cd build
  $ cmake ..
  $ make -j 4
  $ sudo make install

  When you run ``cmake ..`` you should find out what dependencies you are missing and can install them as needed.
+ **ERFA** - https://github.com/liberfa/erfa/releases. I think it's best to install a release version of ``ERFA``. Comes with a ``configure`` file, while the ``git`` repo doesn't. An installation route would look like::

  $ wget https://github.com/liberfa/erfa/releases/download/v2.0.0/erfa-2.0.0.tar.gz
  $ tar -xvf erfa-2.0.0.tar.gz
  $ cd erfa-2.0.0
  $ ./configure
  $ make -j 4
  $ sudo make install
+ **HDF5** - https://www.hdfgroup.org/downloads/hdf5/ - just do::

  $ sudo apt install libhdf5-serial-dev
+ **PAL** - https://github.com/Starlink/pal/releases - ``PAL`` is a little mental with it's default installation paths. I *HIGHLY* recommend downloading a release version, and then using the ``--without-starlink`` option::

  $ wget https://github.com/Starlink/pal/releases/download/v0.9.8/pal-0.9.8.tar.gz
  $ tar -xvf pal-0.9.8.tar.gz
  $ cd pal-0.9.8
  $ ./configure --prefix=/usr/local --without-starlink
  $ make
  $ sudo make install

  Doing it this way installs things in normal locations, making life easier during linking.
+ **python >= 3.6** - the best way to run ``WODEN`` is through the script ``run_woden.py``, which has a number of package dependencies. One of these is ``pyerfa``, which uses f-strings during installation, so you have to use a python version >= 3.6. Sorry. The requirements can be found in ``WODEN/docs/sphinx/sphinx/requirements_testing.txt``, which you can install via something like::

  $ pip3 install -r requirements_testing.txt

For completeness, those packages are::

  sphinx_argparse
  breathe
  astropy
  numpy
  pyerfa
  palpy
  matplotlib

The ``sphinx_argparse, breathe`` packages are used for the documentation. Further packages of ``palpy, matplotlib`` are only used in the ``test_installation/absolute_accuracy`` test, so if you're aiming for a minimal installation, you only need ``numpy, astropy, and pyerfa``.

Phew! That's it for now.

Compiling ``WODEN``
######################

In an ideal world, if the installation of your dependencies went perfectly and
you have a newer NVIDIA GPU, you should be able to simply run::

  $ git clone https://github.com/JLBLine/WODEN.git
  $ cd WODEN
  $ mkdir build && cd build
  $ cmake ..
  $ make -j 4

et voila, your code is compiled. If this worked, and you're happy to install ``WODEN`` into the system default location, just run::

  $ sudo make install

(usually the default is something like ``/usr/local`` hence you need admin privileges). If complilation fails or you're not used to ``cmake``, check out the 'Machine specifics' for help. If you don't want to install or don't have admin rights, head to the 'Post Compilation' section below to finish off your installation.

.. warning:: Even if the code compiled, if your GPU has a compute capability < 5.1, newer versions of ``nvcc`` won't compile code that will work. You'll get error messages like "No kernel image available". Check out how to fix that in 'Machine specifics' below.

Machine specifics
######################
``cmake`` is pretty good at trying to find all the necessary libraries, but every machine is unique, so often you'll need to point ``cmake`` in the correct direction. To that end, I've include 4 keywords: ``JSONC_ROOT``, ``ERFA_ROOT``, ``HDF5_ROOT``, ``PAL_ROOT`` that you can pass to ``cmake``. When passing an option to ``cmake``, you add ``-D`` to the front. For example, on ``OzStar``, I used the command::

  $ cmake ..  -DJSONC_ROOT=/fred/oz048/jline/software/json-c/install/

which tells ``cmake`` to look for ``libjson-c.so`` in paths like ``${JSONC_ROOT}/lib`` or ``${JSONC_ROOT}/lib64``, and ``json.h`` in paths like ``${JSONC_ROOT}/include`` and ``${JSONC_ROOT}/include/json-c``. Read the errors out of ``cmake`` to see which libraries it can't find and add whatever you need to your ``cmake`` command to point to the correct libraries.

.. note:: If you install a dependency in an unusual place on you machine, you have to make sure ``woden`` can find it at run time. So if you compiled with the ``json-c`` library in the ``cmake`` example above, you'd need to call ``export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/fred/oz048/jline/software/json-c/install/lib64`` before you call ``woden`` (or put that line in your ``~/.bashrc`` or equivalent).

All NVIDIA GPUs have a specific compute capability, which relates to their internal architecture. You can tell the compiler which architecture to compile for, which in theory should make compilation quicker, and ensure the code runs correctly on your GPU. You can find out the compute value here (https://developer.nvidia.com/cuda-gpus), and pass it to CMake via::

  $ cmake .. -DCUDA_ARCH=6.0

(for a compute capability of 6.0, for example).

.. warning:: For newer ``CUDA`` versions, some compute capabilities are deprecated, so the compiler leaves them out by default. For example, using ``CUDA`` version 11.2, compute capabilities 3.5 to 5.0 are ignored. If you card has a compute capability of 5.0, you **must** include the flag ``-DCUDA_ARCH=5.0``, otherwise the `nvcc` compiler will not create an executable capable of running on your device.

If you need to pass extra flags to your CUDA compiler, you can do so by adding something like the following::

  -DCMAKE_CUDA_FLAGS="-Dsomeflag"


Post compilation (required if you don't run ``make install``)
###############################################################

If you don't run ``make install``, ``run_woden.py`` won't be able to find the ``woden`` executable. Default installation locations often need admin privileges. If you can't install to them (or just want to keep ``WODEN`` contained inside a single directory), you can instead just add::

  source /path/to/your/location/WODEN/build/init_WODEN.sh

to your ``~/.bash_rc`` (where you replace ``/path/to/your/location`` to wherever you installed ``WODEN``). This will create the variable ``$WODEN_DIR``, and add it to your ``$PATH``. Furthermore, ``init_WODEN.sh`` is generated by the script ``src/update_init_WODEN.py``, which looks through ``CMakeCache.txt`` for the locations of ``ERFA``, ``HDF5``, ``JSONC``, ``PAL``. It then appends lines to ``init_WODEN.sh`` to add these locations to ``LD_LIBRARY_PATH``, so ``woden`` can find these libraries at run time. For example, on my machine, ``init_WODEN.sh`` ends up looking like::

  ##This line finds the current directory at sets the env variable WODEN_DIR
  export WODEN_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
  ##This adds the line to PATH
  export PATH=$WODEN_DIR:$PATH
  ##Add library paths to LD_LIBRARY_PATH so the can be found at runtime
  export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial/:$LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=/usr/local/lib/:$LD_LIBRARY_PATH

.. note:: Every time you run ``make``, ``init_WODEN.sh`` is regenerated, so any edits you make will be overwritten. I suggest any other customisation of you ``LD_LIBRARY_PATH`` happens in your ``~/.bashrc`` or equivalent.

Post compilation (optional)
###############################

If you want to use the MWA FEE primary beam model, you must have the stored spherical harmonic coefficients hdf5 file ``mwa_full_embedded_element_pattern.h5``. You can then define this environment variable in your ``~/.bash_rc``::

  export MWA_FEE_HDF5=/path/to/your/location/mwa_full_embedded_element_pattern.h5

so again ``run_woden.py`` can find it. There is a command line option ``--hdf5_beam_path`` in ``run_woden.py`` which you can use instead of this environment variable if you want.

If you don't have the spherical harmonic file you can obtain it via the command::

  $ wget http://ws.mwatelescope.org/static/mwa_full_embedded_element_pattern.h5
API Reference
==============

Here you can find documentation of all the functions grouped by coding language.
``python`` script documentation also includes a summary of all arguments, which
can be reproduced on the command line using ``script.py --help``.

.. note::  All functions that begin with ``RTS`` are either borrowed or adapted
   from the RTS (`Mitchell et al. 2008`_) calibration package, with permission
   from the original authors. All credit to the original authors. The RTS code
   can currently be found at the `RTS github`_.

.. _Mitchell et al. 2008: https://doi.org/10.1109/JSTSP.2008.2005327
.. _RTS github: https://github.com/ICRAR/mwa-RTS.git

``python`` code
-----------------

.. toctree::
  :maxdepth: 1

  python_code/run_woden
  python_code/convert_WSClean_list_to_WODEN
  python_code/uv2ms

``C`` code
----------------------

The precision of most functions is determined at compilation time, via the
following ``#ifdef`` statement found in ``WODEN/include/woden_precision_defs.h``:

.. code-block:: C

    #ifdef DOUBLE_PRECISION
    /*! If -DDOUBLE_PRECISION flag is added at compilation,
    then user_precision_t is set to double */
    typedef double user_precision_t;
    /*! If -DDOUBLE_PRECISION flag is added at compilation,
    then user_precision_complex_t is set to double _Complex */
    typedef double _Complex user_precision_complex_t;
    #else
    /*! If -DDOUBLE_PRECISION flag is NOT added at compilation,
    then user_precision_t defaults to float */
    typedef float user_precision_t;
    /*! If -DDOUBLE_PRECISION flag is NOT added at compilation,
    then user_precision_complex_t defaults to float _Complex */
    typedef float _Complex user_precision_complex_t;
    #endif

So where you see ``user_precision_t`` in the API, this is either a ``float`` or
``double``, and similarly ``user_precision_complex_t`` is either a ``float _Complex``
or ``double _Complex``.

.. toctree::
  :maxdepth: 1

  C_code/array_layout
  C_code/chunk_sky_model
  C_code/constants
  C_code/create_sky_model
  C_code/FEE_primary_beam
  C_code/primary_beam
  C_code/print_help
  C_code/shapelet_basis
  C_code/visibility_set
  C_code/woden_precision_defs
  C_code/woden_settings
  C_code/woden_struct_defs

``CUDA`` code
-------------------------

Similarly to ``C`` code, the precision of most functions is determined at
compilation time, via the following ``#ifdef`` statement found in
``WODEN/include/cudacomplex.h``:

.. code-block:: C

   #ifdef DOUBLE_PRECISION
   /*! If -DDOUBLE_PRECISION flag is added at compilation,
   then cuUserComplex is set to cuDoubleComplex */
   typedef cuDoubleComplex cuUserComplex;
   #else
   /*! If -DDOUBLE_PRECISION flag is NOTE added at compilation,
   then cuUserComplex is set to cuFloatComplex */
   typedef cuFloatComplex cuUserComplex;
   #endif

meaning that ``cuUserComplex`` in the API either means ``cuFloatComplex`` or
``cuDoubleComplex`` depending on compilation flags.

.. toctree::
  :maxdepth: 1

  CUDA_code/calculate_visibilities
  CUDA_code/cudacheck
  CUDA_code/cudacomplex
  CUDA_code/FEE_primary_beam_cuda
  CUDA_code/fundamental_coords
  CUDA_code/primary_beam_cuda
  CUDA_code/source_components
``create_sky_model``
=====================

API documentation for ``create_sky_model.c``.

TODO Should stick a bunch of information about sky model formats in here

.. doxygenfile:: create_sky_model.h
   :project: WODEN
``shapelet_basis``
===================
API documentation for ``shapelet_basis.c``

.. doxygenfile:: shapelet_basis.h
   :project: WODEN
``primary_beam``
================

API documentation for ``primary_beam.c``.

.. doxygenfile:: primary_beam.h
   :project: WODEN
``FEE_primary_beam``
====================
The MWA Fully Embedded Element primary beam model (`Sokolowski et al. 2017`_) is
at the time of writing the most accurate and advanced MWA primary beam model.
The model is stored in spherical harmonic coefficients, at a frequency
resolution of 1.28MHz, in an ``hdf5`` file called
``mwa_full_embedded_element_pattern.h5``. I have done my best to document the RTS
code I have used. The RTS code makes use of the Legendre Polynomial code
`available here`_ by John Burkardt, and is included in the ``WODEN`` distribution
as ``legendre_polynomial.c``.

.. _Sokolowski et al. 2017: https://doi.org/10.1017/pasa.2017.54
.. _available here: https://people.sc.fsu.edu/~jburkardt/c_src/laguerre_polynomial/laguerre_polynomial.html

.. doxygenfile:: FEE_primary_beam.h
   :project: WODEN
``visibility_set``
===================

API documentation for ``visibility_set.c``.

.. doxygenfile:: visibility_set.h
   :project: WODEN
``woden_precision_defs``
=========================

A header to set the precision required internally to ``WODEN``.

.. note:: I haven't worked out how to get doxygen/sphinx to show the whole defintion here, which includes an `#ifdef`, so only shows the definition for float below.

.. doxygenfile:: woden_precision_defs.h
   :project: WODEN
``chunk_sky_model``
=====================

Functions to take a WODEN sky model and break it up (chunk) into a number of
smaller sky models, to distribute the processing across the GPU.

.. doxygenfile:: chunk_sky_model.h
   :project: WODEN
``constants``
===============

API documentation for ``constants.h``.

.. doxygenfile:: constants.h
   :project: WODEN
``woden_struct_defs``
=======================

API documentation for ``woden_struct_defs.h``.

.. doxygenfile:: woden_struct_defs.h
   :project: WODEN
``array_layout``
================

API documentation for ``array_layout.c``.

.. doxygenfile:: array_layout.h
   :project: WODEN
``print_help``
===================

API documentation for ``print_help.h``.

.. doxygenfile:: print_help.h
   :project: WODEN
``woden_settings``
===================

API documentation for ``woden_settings.c``.

.. doxygenfile:: woden_settings.h
   :project: WODEN
``uv2ms.py``
=============

Helper script to convert ``uvfits`` files to ``measurement sets`` using ``casa``.

.. _uv2ms command line running options:

*Command line running options*
-------------------------------

.. argparse::
   :filename: ../../src/uv2ms.py
   :func: get_parser
   :prog: uv2ms.py
``convert_WSClean_list_to_WODEN.py``
=====================================

Convenience script included to convert ``WSClean`` CLEAN components into a
sky model compatible with ``WODEN``.

.. _convert_wsclean command line running options:

*Command line running options*
-------------------------------

.. argparse::
   :filename: ../../src/convert_WSClean_list_to_WODEN.py
   :func: get_parser
   :prog: convert_WSClean_list_to_WODEN.py
``run_woden.py``
=================

This is the main way to control the ``WODEN`` executable, and ensure good
arguments are supplied. I suggest you only run ``WODEN`` through this script.

.. _run_woden command line running options:

*Command line running options*
-------------------------------

.. argparse::
   :filename: ../../src/run_woden.py
   :func: get_parser
   :prog: run_woden.py

*Function documentation*
------------------------

.. automodule:: run_woden
   :members:
``primary_beam_cuda``
======================
API documentation for ``primary_beam_cuda.cu``.

.. doxygenfile:: primary_beam_cuda.h
   :project: WODEN
``cudacheck``
==============

API documentation for ``cudacheck.h``.

.. doxygenfile:: cudacheck.h
   :project: WODEN
``calculate_visibilities``
===========================

API documentation for ``calculate_visibilities.cu``

.. doxygenfile:: calculate_visibilities.h
   :project: WODEN
``fundamental_coords``
=======================
API documentation for ``fundamental_coords.cu``

.. doxygenfile:: fundamental_coords.h
   :project: WODEN
``cudacomplex``
================

This header is in the ``RTS``, and contains useful CUDA operators. All credit to
the original author, R. G. Edgar. I've added in ``double`` definitions, as well
as my ``cuUserComplex`` def, which allows ``float`` or ``double`` to be
selected during compilation using the flag ``--DDOUBLE_PRECISION`` (so even
though in the below it says ``typedef cuFloatComplex cuUserComplex``, this
depends on compilation).

.. doxygenfile:: cudacomplex.h
   :project: WODEN
``FEE_primary_beam_cuda``
=========================

The MWA Fully Embedded Element primary beam model (`Sokolowski et al. 2017`_) is
at the time of writing the most accurate and advanced MWA primary beam model.
The model is stored in spherical harmonic coefficients, at a frequency
resolution of 1.28MHz, in an ``hdf5`` file called
``mwa_full_embedded_element_pattern.h5``. I have done my best to document the RTS
code I have used. The RTS code makes use of the Legendre Polynomial code
`available here`_ by John Burkardt, and is included in the ``WODEN`` distribution
as ``legendre_polynomial.c``.

.. _Sokolowski et al. 2017: https://doi.org/10.1017/pasa.2017.54
.. _available here: https://people.sc.fsu.edu/~jburkardt/c_src/laguerre_polynomial/laguerre_polynomial.html

.. doxygenfile:: FEE_primary_beam_cuda.h
   :project: WODEN
``source_components``
======================
API documentation for ``source_components.c``

.. doxygenfile:: source_components.h
   :project: WODEN
.. _`this googledoc link to woden-srclist_pumav3.txt`: https://drive.google.com/file/d/1GFnQPXVGsS_7eE5EKTuRp6naNO6IHNFI/view?usp=sharing

MWA EoR1 simulation
====================

.. note:: Running the simulation and making all the images will take up around 1.8 GB storage.

You'll need to download the skymodel from `this googledoc link to woden-srclist_pumav3.txt`_ and put it in the correct directory. If you're comfortable with ``wget`` you can do::

  $ cd WODEN/examples/MWA_EoR1
  $ wget 'https://docs.google.com/uc?export=download&id=1GFnQPXVGsS_7eE5EKTuRp6naNO6IHNFI' -O woden-srclist_pumav3.txt

This skymodel contains over 300,000 sources, and is based on the GLEAM catalogue, with some embellishment. The simulation we're going to run is of the MWA 'EoR1' field, centred at RA, Dec = :math:`60^\circ, -30^\circ`. This field contains Fornax A, and is a good way to demonstrate a sky with point, Gaussian, and shapelet models, as well as the effect of the MWA FEE primary beam. To run the simulation, simply run::

  $ ./MWA_EoR1_simulation.sh

Which contains the command::

  time run_woden.py \
    --ra0=60.0 --dec0=-27.0 \
    --num_freq_channels=16 --num_time_steps=14 \
    --freq_res=80e+3 --time_res=8.0 \
    --cat_filename=woden-srclist_pumav3.txt \
    --metafits_filename=../metafits/1136380296_metafits_ppds.fits \
    --band_nums=1,2,3,4,5 \
    --output_uvfits_prepend=./data/MWA_EoR1 \
    --primary_beam=MWA_FEE \
    --sky_crop_components

Running this took 55 mins 46 secs seconds on my GPU (running with the ``--precision=float`` flag runs in 10 min 39 sec). I've reduced the time and frequency resolution as specified in the ``metafits`` file to keep the size of the outputs smaller on your machine. If you wanted to run the full resolution data of this observation, (2s, 40kHz), you can just remove the ``--num_freq_channels, --num_time_steps, --freq_res, --time_res`` arguments.

I've included two imaging commands::

  $ ./MWA_EoR1_smaller_imaging.sh
  $ ./MWA_EoR1_larger_imaging.sh

The 'larger' image is 6000 by 6000 pixels and has large *w*-terms, so it took about an hour on my desktop. If you don't want to hang around, just run the smaller imaging. The smaller image looks like this:

.. image:: MWA_EoR1_plot_smaller.svg
   :width: 600px

Here I've blown up the colour scale just to highlight how many sources there are (especially FornaxA booming in the corner). The CLEAN isn't great here as I haven't made the image big enough (you can see some aliasing of Fornax A around the image).

If you have the patience to make the bigger image, it looks like this:

.. image:: MWA_EoR1_plot_larger.svg
   :width: 600px

On the left we see the full image, which clearly shows the main lobe of the MWA primary beam. On the right I have zoomed into north of the main lobe, and you can see sources sat in the northern beam sidelobe of the primary at the top, and the edge of the primary lobe at the bottom. This is as expected due to the MWA primary beam shape (see the :ref:`MWA Fully Embedded Element` section if you are unfamiliar with the beam).
.. _WODEN demonstrated via examples:

``WODEN`` demonstrated via examples
=====================================

Below are a number of example simulations with various settings. You can find and run these examples in the ``WODEN/examples`` directory. I'll let you know how much storage you'll need for each set of simulation and images (things can get large with radio data).

Two of the sky models are large (a total of about 100 MB), so instead of bundling them into the github, I've added links into the relevant instructions to download them.

.. note:: For all simulation times reported in the below, I used a single NVIDIA GeForce GTX 1080 Ti with 12 GB of RAM.

The examples are:

.. toctree::

   fornaxA_sim
   MWA_EoR1_sim
   eda2_haslam_sim

**Fornax A simulation** - two examples that compare a point/Gaussian model to a shapelet model. This serves as a general introduction on how to simulate an MWA observation using a ``metafits`` file, with a few extra commands. It also serves as a comparison of running with ``woden_float`` and ``woden_double``.

**MWA EoR1 simulation** - demonstrates using a larger (>300,000) source catalogue

**EDA2 Haslam Map simulation** - this demonstrates using ``WODEN`` without a ``metafits`` file, using a text file to describe the array layout, and using the ``EDA2`` beam.

.. warning:: If you have a GPU with small amounts of RAM (say 2GB) some of these simulations won't work to DOUBLE precision, you won't have enough memory. You can add the ``--precision=float`` argument to switch to a lower memory requirement (for the loss of accuracy).
.. _`this googledoc link to pygsm_woden-list_100MHz_n256.txt`: https://drive.google.com/file/d/1TEELux33UClRTiZBFOzGJHF-XbLnZjUV/view?usp=sharing
.. _`pygdsm`: https://github.com/telegraphic/pygdsm

EDA2 Haslam Map simulation
===========================

.. note:: Running the simulation and making all the images will take up around 1.8 GB storage.

In this simulation, we'll use an all-sky healpix image with nside 256 (generated using `pygdsm`_). For the sky model, I have converted every healpixel into a point source. You'll need to download the skymodel from `this googledoc link to pygsm_woden-list_100MHz_n256.txt`_ and put it in the correct directory. If you're comfortable with ``wget`` you can do::

  $ cd WODEN/examples/EDA2_haslam
  $ wget 'https://docs.google.com/uc?export=download&id=1TEELux33UClRTiZBFOzGJHF-XbLnZjUV' -O pygsm_woden-list_100MHz_n256.txt

To run the command, do::

  $ ./EDA2_haslam_simulation.sh

which took 2 hours 11 mins on my machine (this is running in DOUBLE precision, it takes 1 hour 15 mins at FLOAT precision). This simulates 393,216 point sources for an array of 255 antennas. The command run is::

  run_woden.py \
    --ra0=74.79589467 --dec0=-27.0 \
    --time_res=10.0 --num_time_steps=10 \
    --freq_res=10e+3 --coarse_band_width=10e+4 \
    --lowest_channel_freq=100e+6 \
    --cat_filename=pygsm_woden-list_100MHz_n256.txt \
    --array_layout=../../test_installation/array_layouts/EDA2_layout_255.txt \
    --date=2020-02-01T12:27:45.900 \
    --output_uvfits_prepend=./data/EDA2_haslam \
    --primary_beam=EDA2 \
    --sky_crop_components \
    --band_nums=1,2,3,4,5

Here is a line by line explanation of the command.

::

  --ra0=74.79589467 --dec0=-27.0

sets the phase centre of the simulation.

::

  --time_res=10.0 --num_time_steps=10

means there will be 10 time samples with 10 seconds between each sample.

::

  --freq_res=10e+3 --coarse_band_width=10e+4 \
  --lowest_channel_freq=100e+6 --band_nums=1,2,3,4,5

this combination of arguments will create 5 uvfits file outputs, each containing 10 frequency channels of width 10 kHz. The lowest band will start at 100 MHz, giving a total frequency coverage from 100 MHz to 100.5 MHz.

::

  --cat_filename=pygsm_woden-list_100MHz_n256.txt

points towards the sky model.

::

  --array_layout=../../test_installation/array_layouts/EDA2_layout_255.txt

points towrads an array file that contains local east, north, height coordinates (in metres). This is used in conjunction with latitude to generate baseline coordinates. The default ``--latitude`` is set to the MWA which is right next to the EDA2 so good enough for the example.

::

  -date=2020-02-01T12:27:45.900

sets a UTC date which is used in conjunction with ``--longitude`` to calculate the LST (again, defaults to MWA which is good for purpose here).


::

  --output_uvfits_prepend=./data/EDA2_haslam

sets the naming convention for the outputs, in conjunction with ``--band_nums=1,2,3,4,5`` will produce the outputs::

  ./data/EDA2_haslam_band01.uvfits
  ./data/EDA2_haslam_band02.uvfits
  ./data/EDA2_haslam_band03.uvfits
  ./data/EDA2_haslam_band04.uvfits
  ./data/EDA2_haslam_band05.uvfits

\

::

  --primary_beam=EDA2

selects the EDA2 primary beam.

::

  --sky_crop_components

this means that sky model is cropped by ``COMPONENT`` and not by ``SOURCE``. This model has the diffuse sky as a single ``SOURCE``, so some ``COMPONENT`` s are always below the horizon so need this flag to not crop the whole sky out.

.. note:: The real EDA2 instrument has 256 antennas. ``CASA`` only allows a maximum of 255 elements in an array table, so imaging becomes a nightmare. For this example, to make an image, I've just left out an antenna to make my life easier.

Once you've run that, you can make an image via::

  $ ./EDA2_haslam_imaging.sh

and you'll see this:

.. image:: EDA2_all_sky_image.png
  :width: 600px

where we can see that the EDA2 can see essentially the whole sky, albeit at poor resolution.
.. _`Line et al. 2020`: https://doi.org/10.1017/pasa.2020.18

Fornax A simulation
=========================================

This example not only compares two sky model types, but compares the speed of the
``float`` vs the ``double`` precision of ``WODEN``.

.. note:: Running and imaging both Fornax A simulations will need 3.6 GB of storage.

If you're itching to immediately run something, you can run and image the first simulation with the following commands. Note you'll have to have followed :ref:`Post compilation (optional)` to get the MWA FEE beam working.

Read on to find out what you've done or read the section before running the commands.

::

  $ cd WODEN/examples/FornaxA
  $ ./FornaxA_msclean_simulation.sh
  $ ./FornaxA_msclean_imaging.sh

In this simulation we'll compare two different models of Fornax A which were used in `Line et al. 2020`_. First up we'll look at a model made from ``WSClean`` multi-scale CLEAN outputs. The command looks like this::

  run_woden.py \
    --ra0=50.67 --dec0=-37.2 \
    --freq_res=80e+3 --num_freq_channels=16 \
    --time_res=8.0 --num_time_steps=14 \
    --metafits_filename=../metafits/1202815152_metafits_ppds.fits \
    --cat_filename=srclist_msclean_fornaxA_phase1+2.txt \
    --band_nums=1,2,3,4,5 \
    --output_uvfits_prepend=./data/FornaxA_msclean \
    --primary_beam=MWA_FEE \
    --precision=${precision}

where ``${precision}`` is either "float" or "double" to choose the precision of the simulation. Here we have set the phase centre to RA, Dec = 50.67, -37.2, the fine channel frequency to 80 kHz with 16 fine channels per coarse band, the time sampling to every 8 seconds for 14 time samples, and used a ``metafits`` file for all other observational settings. I've selected to run the first 5 coarse bands via the ``--band_nums`` parameters, which combined with the ``--output_uvfits_prepend`` argument should create 10 ``uvfits`` files::

  data/FornaxA_msclean_${precision}_band01.uvfits
  data/FornaxA_msclean_${precision}_band02.uvfits
  data/FornaxA_msclean_${precision}_band03.uvfits
  data/FornaxA_msclean_${precision}_band04.uvfits
  data/FornaxA_msclean_${precision}_band05.uvfits

each of which will contain 16 frequency and 14 time steps. I've selected to use the ``MWA_FEE`` primary beam, which will use the MWA fully embedded element (FEE) primary beam pattern (using the delays specified in the ``metafits`` to point the beam). As described in :ref:`Post compilation (optional)`, you'll need to grab an hdf5 file and set an environment variable to point to it for this to work.

The sky model is specified using ``--cat_filename``, where we have used ``convert_WSClean_list_to_WODEN.py`` to convert outputs from WSClean into a ``WODEN`` sky model via::

  convert_WSClean_list_to_WODEN.py \
    --file=msclean_output_from_WSClean-sources.txt \
    --outname=srclist_msclean_fornaxA_phase1+2.txt

The sky model contains 4544 point and 1736 Gaussian components. The "float" precision version took about 48 seconds on my card, with the "double" taking about 144 seconds. If you run the imaging, you should get something that looks like this:

.. image:: FornaxA_msclean-image.png
   :width: 400pt

This is an MWA phase II extended array simulation, hence the ~ arcmin resolution. Next, we can compare this to an equivalent shapelet simulation, where the only change to the command is to change the sky model::

  --cat_filename=srclist_shapelets_fornaxA_phase1+2.txt

You can run and image the shapelet simulation via::

  $ ./FornaxA_shapelet_simulation.sh
  $ ./FornaxA_shapelet_imaging.sh

with the "float" simulation taking 59 seconds on my GPU, the "double" taking 144 seconds, and the image looking like:

.. image:: FornaxA_shapelets-image.png
   :width: 400pt
``WODEN`` sky model format
============================

.. _Line et al. 2020: https://doi.org/10.1017/pasa.2020.18
.. _SHAMFI readthedocs: https://shamfi.readthedocs.io/en/latest/


The ``WODEN`` source catalogue is a modified version of the ``RTS`` srclist. In the current version of ``WODEN``, you create one single SOURCE which can include as many COMPONENTS as desired, each of type ``POINT``, ``GAUSSIAN`` or ``SHAPELET``. A ``POINT`` is a dirac delta point source model, a GAUSSIAN is a 2D Gaussian model (with a major, minor, and position angle), and a ``SHAPELET`` model uses multiple 'shapelet' basis functions to build a model. For details on the model types, see `Line et al. 2020`_. If you want to build a shapelet model, you can use the software ``SHAMFI``, which you can read about on the `SHAMFI readthedocs`_.

Currently, every source is given a simple power-law frequency behaviour as:

.. math::
  S = S_0 \left( \frac{\nu_0}{\nu} \right)^\alpha

where :math:`S` is the flux density at frequency :math:`\nu`, with a reference flux density :math:`S_0`, reference frequency :math:`\nu_0`, and spectral index  :math:`\alpha`.

Point sources
^^^^^^^^^^^^^^^^^^^^

An example of a single SOURCE with a single point source COMPONENT is::

  SOURCE source_name P 1 G 0 S 0 0
  COMPONENT POINT 4.0 -27.0
  LINEAR 1.8e+08 10.0 0 0 0 -0.8
  ENDCOMPONENT
  ENDSOURCE

An explanation of each line and value follows.

::

  SOURCE source_name P 1 G 0 S 0 0

Initialises the SOURCE, giving it the name ``source_name``, and specifying the number and type of components (P = point, G = gaussian, S = shapelet). For shapelet, the two numbers are total number of coefficients and total number of components. Read on further for more explanation of shapelets.

::

  COMPONENT POINT 4.0 -27.0

Initialises a component, specifying the type (either POINT, GAUSSIAN, SHAPELET) and the RA and DEC (hours, deg). So this line means a point source at RA,DEC = 4h, -27deg.

::

  LINEAR 1.8e+08 10.0 0 0 0 -0.8

Specifies a reference Stokes flux density as *LINEAR Freq I Q U V SI*, where the Freq is in Hz, Stokes params *I,Q,U,V* are all in units of Jy, and SI is the spectral index. It's labelled ``LINEAR`` as a power-law is linear in log-log space. This example line specifies we have a source that has a flux density of purely Stokes I of 10 Jy at 180 MHz, with a spectral index if -0.8.

::

  ENDCOMPONENT

This line ends the component.

::

  ENDSOURCE

This line ends the source.


To add multiple point sources, simply repeat the ``COMPONENT`` / ``ENDCOMPONENT`` sections with new details, i.e.

::

  SOURCE multi_point P 3 G 0 S 0 0
  COMPONENT POINT 4.0 -27.0
  LINEAR 1.8e+08 10.0 0 0 0 -0.4
  ENDCOMPONENT
  COMPONENT POINT 3.0 -37.0
  LINEAR 1.3e+08 1.0.0 0 0 0 -0.786
  ENDCOMPONENT
  COMPONENT POINT 5.0 -47.0
  LINEAR 3.9e+08 0.04 0 0 0 .02
  ENDCOMPONENT
  ENDSOURCE

noting that at the very top line, I have updated ``P 3`` to reflect there are now three point sources. These numbers are used to quickly allocate memory, that's why they re included.

.. note:: ``WODEN`` crops everything below the horizon out of the sky model. It can do this one of two ways - either by ``COMPONENT`` or by ``SOURCE``. In the example above, we have three COMPONENT in one SOURCE. If you ask ``WODEN`` to crop by ``SOURCE``, if just one of the ``COMPONENTS`` is below the horizon, it'll crop the *entire* source.

Gaussian sources
^^^^^^^^^^^^^^^^^^^^

An example srclist containing a single gaussian::

  SOURCE gaussian_source P 0 G 1 S 0 0
  COMPONENT GAUSSIAN 3.378 -37.2
  LINEAR 1.8e+08 10.0 0 0 0 -0.8
  GPARAMS 45.0000000000 6.0 3.0
  ENDCOMPONENT
  ENDSOURCE

where all lines have the same meaning as the point source, and the meaning of the extra line::

  GPARAMS 45.0000000000 6.0 3.0

which specifies the Gaussian parameters as ``GPARAMS pa(deg) major_axis(arcmin) minor_axis(arcmin)``. The major and minor axes are specified as FWHM. Note this line needs to sit in between the lines starting with ```COMPONENT GAUSSIAN`` and ```ENDCOMPONENT``.

Shapelet sources
^^^^^^^^^^^^^^^^^^^^

To generate shapelet models compatible with WODEN, simply use ``SHAMFI`` to fit an image with the ``--woden_srclist`` option (again see `SHAMFI readthedocs`_. for more detail). This will ensure all normalisations are correct. An example sky model (made by hand so the normalisations *won't* be correct) is::

  SOURCE shapelet_source P 0 G 0 S 1 3
  COMPONENT SHAPELET 3.378 -37.2
  FREQ 1.8e+08 10.0 0 0 0
  SPARAMS 45.0000000000 6.0 3.0
  SCOEFF 0 0 0.92342
  SCOEFF 1 10 0.0002354
  SCOEFF 4 5 0.004567
  ENDCOMPONENT
  ENDSOURCE

which generates a single shapelet component, including 3 shapelet basis functions, hence ``S 1 3`` in the first line. The ``SPARAMS`` line is similar to the ``GAUSSIAN`` line with ``SPARAMS pa(deg) major_axis(arcmin) minor_axis(arcmin)``. The extra lines like::

  SCOEFF 0 0 0.92342

encode the order of the shapelet basis function (see `Line et al. 2020`_ for details) and fitted coefficient as ``SCOEFF p1 p2 coeff_value``. You can add as many ``SCOEFF`` lines as necessary, with a maximum order < 100. If you use ``SHAMFI``, the coefficients will be scaled such that the Stokes I flux density of the full source will be 10 Jy at 180 MHz for this example. You may have noticed the SED information here is different::

  FREQ 1.8e+08 10.0 0 0 0

This line will still assume a power-law frequency behaviour, with a reference flux of 10 Jy at 180 MHz, but use a default SI = -0.8.

Putting it all together
^^^^^^^^^^^^^^^^^^^^^^^^^

An example skymodel with four sources, the first with all component types, the next three with a single component of each type,  would look something like this::

  SOURCE multi_sources P 3 G 1 S 2 7
  COMPONENT SHAPELET 3.378 -37.2
  FREQ 1.8e+08 10.0 0 0 0
  SPARAMS 45.0000000000 6.0 3.0
  SCOEFF 0 0 0.92342
  SCOEFF 1 10 0.0002354
  SCOEFF 4 5 0.004567
  ENDCOMPONENT
  COMPONENT SHAPELET 3.12 -32.2
  FREQ 1.8e+08 3.1 0 0 0
  SPARAMS 56.0000000000 9.0 3.0
  SCOEFF 0 0 0.02345
  SCOEFF 3 0 -0.234234
  SCOEFF 21 34 0.82342
  SCOEFF 31 5 -0.00876234
  ENDCOMPONENT
  COMPONENT GAUSSIAN 3.378 -37.2
  LINEAR 1.8e+08 10.0 0 0 0 -0.8
  GPARAMS 45.0000000000 6.0 3.0
  ENDCOMPONENT
  COMPONENT POINT 4.0 -27.0
  LINEAR 1.8e+08 10.0 0 0 0 -0.8
  ENDCOMPONENT
  COMPONENT POINT 3.0 -37.0
  LINEAR 1.8e+08 0.6 0 0.2 0 -0.8
  ENDCOMPONENT
  COMPONENT POINT 5.0 -47.0
  LINEAR 70E+6 87.0 0 0 0 -0.8
  ENDCOMPONENT
  ENDSOURCE
  SOURCE source_name P 1 G 0 S 0 0
  COMPONENT POINT 4.0 -27.0
  LINEAR 1.8e+08 10.0 0 0 0 -0.8
  ENDCOMPONENT
  ENDSOURCE
  SOURCE gaussian_source P 0 G 1 S 0 0
  COMPONENT GAUSSIAN 3.378 -37.2
  LINEAR 1.8e+08 10.0 0 0 0 -0.8
  GPARAMS 45.0000000000 6.0 3.0
  ENDCOMPONENT
  ENDSOURCE
  SOURCE shapelet_source P 0 G 0 S 1 3
  COMPONENT SHAPELET 3.378 -37.2
  LINEAR 1.1e+08 10.0 2.0 0 0.8 -0.7
  SPARAMS 45.0000000000 6.0 3.0
  SCOEFF 0 0 0.92342
  SCOEFF 1 10 0.0002354
  SCOEFF 4 5 0.004567
  ENDCOMPONENT
  ENDSOURCE
.. _Sokolowski et al. 2017: https://doi.org/10.1017/pasa.2017.54
.. _polarised_source_and_FEE_beam.ipynb which lives here: https://github.com/JLBLine/polarisation_tests_for_FEE
.. _Tingay et al. 2013: https://doi.org/10.1017/pasa.2012.007
.. _Wayth et al. 2017: https://doi.org/10.1017/pasa.2017.27

Primary Beams
================
``WODEN`` has been written to include stationary primary beams. That means the beam is pointed at a constant azimuth / zenith angle during an observation. There are currently three primary beams available, which are detailed below.


MWA Fully Embedded Element
----------------------------

The Murchison Widefield Array (MWA, `Tingay et al. 2013`_) has 16 bow-tie dipoles arranged in a 4 by 4 grid as recieving elements, yielding a grating-lobe style primary beam.

``WODEN`` incudes a GPU-implementation of the MWA Fully Embedded Element (FEE) Beam pattern (`Sokolowski et al. 2017`_), which to date is the most accurate model of the MWA primary beam. This model is defined in a spherical harmonic coordinate system, which is polarisation-locked to instrumental azimuth / elevation coordinates. ``WODEN`` however uses Stokes parameters to define it's visibilities, and so a rotation of the beam about parallactic angle (as calculated using ``erfa``) is applied to align the FEE beam to move it into the Stokes frame.

Due to convention issues with whether 'X' means East-West or North-South, and whether azimuth starts towards North and increase towards East, we also find is necessary to reorder outputs and apply a sign flip to two of the outputs of the MWA FEE code. For an *exhaustive* investigation into why this is necessary to obtain the expected Stokes parameters, see `polarised_source_and_FEE_beam.ipynb which lives here`_

I can define the Jones matrix of the primary beam as:

.. math::

  \mathbf{J_\mathrm{linear}} =
    \begin{bmatrix}
    g_{x} & D_{x} \\
    D_{y} & g_{y} \\
    \end{bmatrix}.

Here, the subscript :math:`x` means a polarisation angle of :math:`0^\circ` and :math:`y` an angle of :math:`90^\circ`, :math:`g` means a gain term, and :math:`D` means a leakage term (so :math:`x` means North-South and :math:`y` is East-West). Under this definition, a typical zenith-pointing looks like this:

.. image:: MWAFEE_jones.png
  :width: 400pt

These plots are all sky, with northward at the top. If we assume the sky is totally Stokes I, this will yield instrumental polarisations (where 'XX' is North-South and 'YY' is East-West) like this:

.. image:: MWAFEE_instrumental_pols.png
  :width: 400pt

The MWA beam is electronically steered, which can be defined via integer delays and supplied to the MWA FEE beam. ``run_woden.py`` can read these directly from an MWA metafits file, or can be directly supplied using the ``--MWA_FEE_delays`` argument.


.. warning:: The frequency resolution of the MWA FEE model is 1.28 MHz. I have NOT yet coded up a frequency interpolation, so the frequency response for a given direction looks something like the below. This is coming in the future.

.. image:: MWAFEE_beam_vs_freq.svg
  :width: 400pt

In fact, when running using the MWA FEE band, I only calculate the beam response once per coarse band. If you set your ``--coarse_band_width`` to greater than 1.28 MHz you'll make this effect even worse. If you stick to normal MWA observational params (with the default 1.28 MHz) all will be fine.

EDA2
------

The 2nd version of the Engineering Development Array (EDA2, `Wayth et al. 2017`_), is an SKA_LOW test station, which swaps the planned logarithmic 'christmas tree' dipoles for MWA bow-tie dipoles. Currently, ``WODEN`` just assumes a perfect dipole with an infinite ground screen as a beam model. This makes the primary beam entirely real, with no leakage terms. Explicitly, the beam model is

.. math::

  \mathcal{G} = 2\sin\left(\pi \frac{2h}{\lambda} \cos(\theta) \right) \\
  g_x = \mathcal{G}\arccos\left(\sin(\theta)\cos(\phi)\right) \\
  g_y = \mathcal{G}\arccos\left(\sin(\theta)\sin(\phi)\right)


where :math:`h` is the height of the dipole, :math:`\lambda` is the wavelength, :math:`\theta` is the zenith angle, :math:`\phi` is the azimuth angle. I've set :math:`h=0.3` m.

The beams basically see the whole sky (this image shows some :math:`\mathbf{J_\mathrm{linear}}` values at 70 MHz):

.. image:: EDA2_jones.png
  :width: 400pt

.. note:: The EDA2 beam is neither physically nor electronically steered, so it always points towards zenith.

Gaussian
----------

This is a toy case of a symmetric (major = minor) Gaussian primary beam. The beam gets smaller on the sky with increasing frequency, but both polarisations are identical. You can control the pointing of the beam (which remains constant in az/za for a single observation) via an initial RA/Dec pointing (``--gauss_ra_point``, ``--gauss_dec_point``), and the FWHM of the beam (``--gauss_beam_FWHM``) at a reference frequency (``--gauss_beam_ref_freq``).

I've implemented this beam by creating a cosine angle coordinate system locked to the initial hour angle and declination of the specified RA,Dec pointing :math:`l_\mathrm{beam}, m_\mathrm{beam}, n_\mathrm{beam}`. The beam is then calculated as

.. math::

  G(l_\mathrm{beam}, m_\mathrm{beam}) = \exp \left( -\left( al_\mathrm{beam}^2 + 2bl_\mathrm{beam}m_\mathrm{beam} + cm_\mathrm{beam}^2 \right)  \right)


where

.. math::

  a  =  \frac{\cos(\phi_{\mathrm{PA}})^2}{2\sigma_l^2} + \frac{\sin(\phi_{\mathrm{PA}})^2}{2\sigma_m^2} \\
  b  =  -\frac{\sin(2\phi_{\mathrm{PA}})}{4\sigma_l^2} + \frac{\sin(2\phi_{\mathrm{PA}})}{4\sigma_m^2} \\
  c  =  \frac{\sin(\phi_{\mathrm{PA}})^2}{2\sigma_l^2} + \frac{\cos(\phi_{\mathrm{PA}})^2}{2\sigma_m^2}.

Currently, I have set the position angle of the beam :math:`\phi_{\mathrm{PA}}=0` the std :math:`\sigma_l = \sigma_m` to be equal, as:

.. math::

  \sigma_l = \sigma_m = \frac{\sin(\varphi_0)}{ 2\sqrt{2\ln(2)} }\frac{\nu_0}{\nu}

where :math:`\varphi_0` is the desired FWHM at reference frequency :math:`\nu_0`, and :math:`\nu` is the frequency to calculate the beam at.

An example of a zenith pointing, with :math:`\varphi_0 = 10^\circ, \nu_0=100` MHz looks like:

.. image:: Gaussian_jones_zenith.png
  :width: 400pt

Using the same settings with an off-zenith pointing yields:

.. image:: Gaussian_jones_offzenith.png
  :width: 400pt

which at least visually looks like we are getting realistic-ish projection effects of the beam towards the horizon.

.. note:: The machinery is there to have different major / minor axes and a position angle if this is desired. Just open an `issue on the github`_ if you want this implemented.

.. _`issue on the github`: https://github.com/JLBLine/WODEN/issues
Frequency specifications and uvfits outputs
============================================

``WODEN`` writes visibilities into ``uvfits`` files, which are output in linear polarisations of ``XX, YY, XY, YX``, where ``X`` refers to a North-South aligned receiver and ``Y`` a East-West aligned receiver. These outputs are split up across frequency, into 'coarse bands'. Each coarse band is split into 'fine channels' like so:

.. image:: frequency_wording.svg
   :width: 600pt

When running MWA simulations using a ``metafits`` file, these frequency options are filled in automatically for the user, and can be overridden / alternatively supplied by keywords to ``run_woden.py``. The arguments that match the diagram above are listed here.

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - ``run_woden.py`` argument
     - Diagram label
   * - ``--coarse_band_width``
     - Coarse band width
   * - ``--freq_res``
     - Fine frequency channel width
   * - ``--lowest_channel_freq``
     - Lowest frequency channel
   * - ``--band_nums=1,2,3``
     - Band 01, Band 02, Band 03

As many coarse bands as needed can be run, allowing for straight-forward splitting of the simulation across multiple GPUs. We'll look at an example of that later.  If you run with following arguments::

  run_woden.py \
  --freq_res=10e+3 --coarse_band_width=1e+6 \
  --time_res=1.0 --num_time_steps=10 \
  --lowest_channel_freq=100e+6 \
  --band_nums=1,2,3 \
  --output_uvfits_prepend=epic_output

this will produce three uvfits files named with the following properties:

.. list-table::
   :widths: 20 10 10
   :header-rows: 1

   * - Output name
     - Lowest frequency channel (Hz)
     - Number frequency channels
   * - ``epic_output_band01.uvfits``
     - 100e+6
     - 100
   * - ``epic_output_band02.uvfits``
     - 101e+6
     - 100
   * - ``epic_output_band03.uvfits``
     - 102e+6
     - 100


.. note:: The command above won't work as many arguments are missing; I've left them out here to concentrate on the arguments that define the uvfits outputs.

You can run with whatever band numbers you want::

    run_woden.py \
    --freq_res=10e+3 --coarse_band_width=1e+6 \
    --time_res=1.0 --num_time_steps=10 \
    --lowest_channel_freq=100e+6 \
    --band_nums=4,7,24 \
    --output_uvfits_prepend=epic_output

which will create uvfits files like:

.. list-table::
   :widths: 20 10 10
   :header-rows: 1

   * - Output name
     - Lowest frequency channel (Hz)
     - Number frequency channels
   * - ``epic_output_band04.uvfits``
     - 104e+6
     - 100
   * - ``epic_output_band07.uvfits``
     - 107e+6
     - 100
   * - ``epic_output_band24.uvfits``
     - 124e+6
     - 100
``WODEN`` operating principles
===============================

Read on below for ``WODEN`` specifics, including the analytic equations used to transform the sky models into visibilities.

.. toctree::
   :maxdepth: 1

   visibility_calcs
   skymodel
   frequency_wording
   primary_beams
.. _`Thompson, Moran, & Swenson 2017`: https://link.springer.com/book/10.1007/978-3-319-44431-4
.. _`Line et al. 2020`: https://doi.org/10.1017/pasa.2020.18

Visibility Calculations
========================

This section assumes a basic understanding on radio interferometry, assuming you know what visibilities and baselines are, and are familiar with the :math:`u,v,w` and :math:`l,m,n` coordinate systems. I can recommend `Thompson, Moran, & Swenson 2017`_ if you are looking to learn / refresh these concepts. Some of this section is basically a copy/paste from `Line et al. 2020`_.

.. note:: In `Line et al. 2020`_, I detailed that I was using the ``atomicAdd`` functionality in CUDA. This is no longer true, as I've found a loop inside my kernels is actually faster than being fully parallel and using ``atomicAdd``. The calculations being made have remained the same.

Measurement Equation and Point Sources
----------------------------------------

``WODEN`` analytically generates a sky model directly in visibility space via the measurement equation (c.f. `Thompson, Moran, & Swenson 2017`_). Ignoring the effects of instrumental polarisation and the primary beam, this can be expressed as:

.. math::

  V_s(u,v,w) =   \int \mathcal{S}_s(l,m) \exp[-2\pi i(ul + vm + w(n-1))] \dfrac{dldm}{n},

where :math:`V_s(u,v,w)` is the measured visibility in some Stokes polarisation :math:`s` at baseline coordinates :math:`u,v,w`, given the sky intensity :math:`\mathcal{S}`, which is a function of the direction cosines :math:`l,m`, and :math:`n=\sqrt{1-l^2-m^2}`. This can be discretised for point sources such that

.. math::

    V_s(u_i,v_i,w_i) = \sum_j \mathcal{S}_s(l_j,m_j) \exp[-2\pi i(u_il_j + v_im_j + w_i(n_j-1))],

where :math:`u_i,v_i,w_i` are the visibility co-ordinates of the :math:`i^{\mathrm{th}}` baseline, and :math:`l_j`, :math:`m_j`, :math:`n_j` is the sky position of the :math:`j^{\mathrm{th}}` point source.

Stokes parameters :math:`\mathcal{S}_I, \mathcal{S}_Q, \mathcal{S}_U, \mathcal{S}_V` are all extrapolated from an input catalogue, along with the position on the sky. :math:`u,v,w` are set by a supplied array layout, phase centre, and location on the Earth.

.. note:: :math:`u_i,v_i,w_i` and :math:`\mathcal{S}` are also functions of frequency, so must be calculated for each frequency steps as required.

Apply Linear Stokes Polarisations and the Primary Beam
---------------------------------------------------------
``WODEN`` simulates dual-linear polarisation antennas, with each antenna/station having it's own primary beam shape. I can define the response of a dual polarisation antenna to direction :math:`l,m` as

.. math::
   \mathbf{J}(l,m) =
   \begin{bmatrix}
   g_x(l,m) & D_x(l,m) \\
   D_y(l,m) & g_y(l,m)
   \end{bmatrix},

where :math:`g` are gain terms, :math:`D` are leakage terms, and :math:`x` refers to a north-south aligned antenna, and :math:`y` an east-west aligned antenna. When calculating the cross-correlation responses from antennas 1 and 2 towards direction :math:`l,m` to produce linear polarisation visibilities, these gains and leakages interact with the four Stokes polarisations :math:`I,Q,U,V` as

.. math::
   \begin{bmatrix}
   V_{12\,XX}(l,m) \\
   V_{12\,XY}(l,m) \\
   V_{12\,YX}(l,m) \\
   V_{12\,YY}(l,m)
   \end{bmatrix} =
   \mathbf{J}_1(l,m) \otimes \mathbf{J}_2^*(l,m)
   \begin{bmatrix}
   1 & 1 & 0 & 0 \\
   0 & 0 & 1 & i \\
   0 & 0 & 1 & -i \\
   1 & -1 & 0 & 0
   \end{bmatrix}
   \begin{bmatrix}
   V_{12\,I}(l,m) \\
   V_{12\,Q}(l,m) \\
   V_{12\,U}(l,m) \\
   V_{12\,V}(l,m)
   \end{bmatrix}


where :math:`*` denotes a complex conjugate, and :math:`\otimes` an outer product. Explicitly, each visibility is

.. math::
   \begin{eqnarray*}
   V_{12\,XX} = (g_{1x}g_{2x}^{\ast} + D_{1x}D_{2x}^{\ast})\mathrm{V}^{I}_{12}
     +  (g_{1x}g_{2x}^{\ast} - D_{1x}D_{2x}^{\ast})\mathrm{V}^{Q}_{12} \\
     +  (g_{1x}D_{2x}^{\ast} + D_{1x}g_{2x}^{\ast})\mathrm{V}^{U}_{12}
     +  i(g_{1x}D_{2x}^{\ast} - D_{1x}g_{2x}^{\ast})\mathrm{V}^{V}_{12}
   \end{eqnarray*}
.. math::
   \begin{eqnarray*}
   V_{12\,XY} =
        (g_{1x}D_{2y}^{\ast} + D_{1x}g_{2y}^{\ast})\mathrm{V}^{I}_{12}
     +  (g_{1x}D_{2y}^{\ast} - D_{1x}g_{2y}^{\ast})\mathrm{V}^{Q}_{12} \\
     +  (g_{1x}g_{2y}^{\ast} + D_{1x}D_{2y}^{\ast})\mathrm{V}^{U}_{12}
     +  i(g_{1x}g_{2y}^{\ast} - D_{1x}D_{2y}^{\ast})\mathrm{V}^{V}_{12}
   \end{eqnarray*}
.. math::
   \begin{eqnarray*}
   V_{12\,YX} =
        (D_{1y}g_{2x}^{\ast} + g_{1y}D_{2x}^{\ast})\mathrm{V}^{I}_{12}
     +  (D_{1y}g_{2x}^{\ast} - g_{1y}D_{2x}^{\ast})\mathrm{V}^{Q}_{12} \\
     +  (D_{1y}D_{2x}^{\ast} + g_{1y}g_{2x}^{\ast})\mathrm{V}^{U}_{12}
     +  i(D_{1y}D_{2x}^{\ast} - g_{1y}g_{2x}^{\ast})\mathrm{V}^{V}_{12}
   \end{eqnarray*}
.. math::
   \begin{eqnarray*}
   V_{12\,YY} =
        (D_{1y}D_{2y}^{\ast} + g_{1y}g_{2y}^{\ast})\mathrm{V}^{I}_{12}
     +  (D_{1y}D_{2y}^{\ast} - g_{1y}g_{2y}^{\ast})\mathrm{V}^{Q}_{12} \\
     +  (D_{1y}g_{2y}^{\ast} + g_{1y}D_{2y}^{\ast})\mathrm{V}^{U}_{12}
     +  i(D_{1y}g_{2y}^{\ast} - g_{1y}D_{2y}^{\ast})\mathrm{V}^{V}_{12}
   \end{eqnarray*}

For each baseline, frequency, and time step, ``WODEN`` calculates all four linear polarisations as defined above for all directions :math:`l_j,m_j` in the sky model, and then sums over :math:`j`, to produce four full-sky linear Stokes polarisation visibilities per baseline/frequency/time.


Gaussian and Shapelet sources
------------------------------
You can inject morphology into your sources analytically by tranforming a visibility into a Gaussian or Shapelet source. We utilise the ``RTS`` methodology of inserting a visibility "envelope" :math:`\xi` into the visibility equation:

.. math::

  V(u_i,v_i,w_i) = \sum_j \xi_j(u_i,v_i)\mathcal{S}(l_j,m_j) \exp[-2\pi i(u_il_j + v_im_j + w_i(n_j-1))],

For a Gaussian, this envelope looks like

.. math::

    \begin{align}
    &\xi_j = \exp\left( -\dfrac{\pi^2}{4\ln(2)} \left( k_x^2\theta_\mathrm{maj}^2 + k_y^2\theta_\mathrm{min}^2\right) \right); \\
    &k_x =  \cos(\phi_{\textrm{PA}})v_i + \sin(\phi_{\textrm{PA}})u_i; \\
    &k_y = -\sin(\phi_{\textrm{PA}})v_i + \cos(\phi_{\textrm{PA}})u_i;
    \end{align}

where :math:`\theta_\mathrm{maj}` and :math:`\theta_\mathrm{min}` are the major and minor axes and :math:`\phi_{\textrm{PA}}` the position angle of an elliptical Gaussian.

For a shapelet model, the envelope looks like:

.. math::

    \begin{align}
    &\xi_j = \sum^{p_k +p_l < p_\mathrm{max}}_{k,l} C_{p_k,p_l} \tilde{B}_{p_k,p_l}(k_x,k_y); \label{eq:shape-env} \\
    &k_x =  \dfrac{\pi}{\sqrt{2\ln(2)}} \left[\cos(\phi_{PA})v_{i,j} + \sin(\phi_{PA})u_{i,j} \right]; \label{eq:scale-shape-x} \\
    &k_y = \dfrac{\pi}{\sqrt{2\ln(2)}} \left[-\sin(\phi_{PA})v_{i,j} + \cos(\phi_{PA})u_{i,j} \right], \label{eq:scale-shape-y}
    \end{align}


where :math:`u_{i,j},v_{i,j}` are visibility co-ordinates for baseline :math:`i`, calculated with a phase-centre :math:`RA_j,\delta_j`, which corresponds to the central position :math:`x_0,y_0` used to fit the shapelet model in image-space. The shapelet basis function values :math:`\tilde{B}_{p_k,p_l}(u,v)` can be calculated by interpolating from one dimensional look-up tables of :math:`\tilde{B}(k_x;1)`, and scaling by the appropriate :math:`\beta` (c.f. Equation 1 in `Line et al. 2020`_ - see for a introduction and breakdown of shapelets bais functions).

You can see the difference between the three types of sky model component below. You can generate this plot yourself, checkout the section :ref:`Grid Component Models`.

.. image:: ../testing/grid_component_plots.png
   :width: 800px
