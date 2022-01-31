
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Spectroscopy Analysis Tools (spant) <img src="man/figures/logo.png" align="right" width=130/>

[![R build
status](https://github.com/martin3141/spant/workflows/R-CMD-check/badge.svg)](https://github.com/martin3141/spant/actions)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03646/status.svg)](https://doi.org/10.21105/joss.03646)
[![](http://cranlogs.r-pkg.org/badges/spant)](http://cran.rstudio.com/web/packages/spant/index.html)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/spant)](https://cran.r-project.org/package=spant)
[![Coverage
Status](https://coveralls.io/repos/github/martin3141/spant/badge.svg?branch=master)](https://coveralls.io/github/martin3141/spant?branch=master)

## Overview

spant provides a full suite of tools to build automated analysis
pipelines for Magnetic Resonance Spectroscopy (MRS) data. The following
features and algorithms are included:

-   Advanced fully-automated metabolite fitting algorithm - ABfit
    <https://onlinelibrary.wiley.com/doi/10.1002/mrm.28385>.
-   Robust retrospective frequency and phase correction - RATS
    <https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.27605>.
-   Flexible data types to support single voxel, dynamic and
    spectroscopic imaging data types.
-   Raw data import from individual coils and dynamic measurements, eg
    support for importing individual FIDs from Siemens TWIX formatted
    data.
-   Publication quality plotting.
-   Extensive set of pre-processing steps (phasing, coil-combination,
    zero-filling, HSVD filtering…)
-   Quantum mechanical based simulation for experimental design and
    basis-set generation.
-   Set of metabolite, macromolecule and lipid parameters for typical
    brain analyses.
-   Voxel registration to anatomical images for partial volume
    concentration corrections.

## Basic installation

Download and install the latest version of R
(<https://cloud.r-project.org/>), or with your package manager if using
a recent Linux distribution, eg `sudo apt install r-base`.

It is also strongly recommended to install RStudio Desktop
(<https://rstudio.com/products/rstudio/download>) to provide a modern
environment for interactive data analysis.

Once R and RStudio have been installed, open the RStudio application and
type the following in the Console (lower left panel) to install the
latest stable version of spant:

``` r
install.packages("spant", dependencies = TRUE)
```

Or the the development version from GitHub (requires the `devtools`
package):

``` r
install.packages("devtools")
devtools::install_github("martin3141/spant", ref = "devel", dependencies = TRUE)
```

## Documentation

Quick introduction to the basic analysis workflow :
<https://martin3141.github.io/spant/articles/spant-intro.html>

Short tutorials : <https://martin3141.github.io/spant/articles/>

Function reference : <https://martin3141.github.io/spant/reference/>

Once the spant library has been loaded with `library(spant)`, type
`?spant` on the console for instructions on how to access the offline
documentation. Note that offline help on the available functions can be
quickly shown in RStudio using `?function_name`, eg `?read_mrs`.

## Ubuntu 20.04 installation

CRAN packages need to be compiled on Linux, and therefore you may need
to ensure some additional system libraries are installed. spant may be
installed from a clean installation of Ubuntu 20.04 with the following
commands pasted into the terminal:

``` ubuntu
sudo apt install -y r-base libcurl4-openssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libtiff5-dev
mkdir -p ~/R/x86_64-pc-linux-gnu-library/3.6
Rscript -e 'install.packages("spant", dependencies = TRUE)'
```

## Ubuntu 21.10 installation

``` ubuntu
sudo apt install -y r-base libcurl4-openssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libtiff5-dev
mkdir -p ~/R/x86_64-pc-linux-gnu-library/4.0
Rscript -e 'install.packages("spant", dependencies = TRUE)'
```

## Anaconda environment installation

Firstly install Anaconda in the standard way for your platform:
<https://docs.anaconda.com/anaconda/install/index.html>.

Create a text file, called `spant_requirements.yml`, containing the
following:

``` conda
name: spant
channels:
    - conda-forge
    - r
dependencies:
    - r-base
    - r-essentials
    - r-abind
    - r-plyr
    - r-foreach
    - r-pracma
    - r-stringr
    - r-signal
    - r-matrixcalc
    - r-minpack.lm
    - r-nnls
    - r-r.utils
    - r-graphicsqc
    - r-expm
    - r-smoother
    - r-readr
    - r-magrittr
    - r-ptw
    - r-mmand
    - r-RNifti
    - r-RNiftyReg
    - r-fields
    - r-MASS
    - r-numDeriv
    - r-nloptr
    - r-irlba
    - r-tibble
    - r-jsonlite
    - r-viridisLite
    - r-shiny
    - r-miniUI
    - r-knitr
    - r-rmarkdown
    - r-testthat
    - r-ragg
    - r-doParallel
```

Create and activate the environment:

``` conda
conda env create -f spant_requirements.yml
conda activate spant
```

Start R and install spant in the usual way:

``` r
install.packages("spant", dependencies = TRUE)
```

Big thanks to [João M.C. Teixeira](https://github.com/joaomcteixeira)
for figuring out this method of install.
# spant 1.17.0
* New function added (recon_twix_2d) for reconstructing basic phase encoded twix
2D MRSI data.
* File name argument for read_mrs now permits globbing, eg read_mrs("*.dcm").
* Improved plotting of metabolite maps containing infinite values.
* Improvements to GE p-file reader.
* Bug fix for TWIX MRSI voxel dimensions.
* get_mrsi2d_seg now returns partial volume maps as well as a data frame.
* Improved error handling for spec_op function.
* apodise_xy can now be applied to multi-coil and dynamic datasets.

# spant 1.16.0
* The package has been published in the Journal for Open Source Software :
"spant: An R package for magnetic resonance spectroscopy analysis. JOSS 2021,
6(67), 3646".
* The basis result of the HSVD function is now an mrs_data object.
* The complete model signal is now output by the HSVD function.
* image function x-axis updated to be consistent with other plotting methods.
* Minor refactor of the simulation code and a dependency swap from complexplus
to exmp packages.

# spant 1.15.0
* Added a unit test for reading and writing LCM .basis formatted files.
* Added FWHM estimates for tCr and tCho in ABfit.
* Added water suppression efficiency and water FWHM measures to the diagnostic
output of svs_1h_brain_analysis.
* Added get_head_dyns and get_tail_dyns to return the first and last dynamic
scans within a dataset.
* Fixed CI errors.
* Improved installation instructions.
* Removed comb_csv_results function and reduced the number of required packages.

# spant 1.14.0
* Added glycerol simulation parameters, e.g. get_mol_paras("glyc").
* Bug fix for read_ima_* functions.
* Improved y = 0 baseline for stackplot when setting bl_lty parameter.
* Removed norm_mrs function and replaced with scale_spec for simple data
scaling tasks, eg scaling based on the integration of a spectral region.
* Added spec_op function for performing simple summary operations on spectral
segments.
* Changed the name of scale_mrs function to scale_mrs_amp.
* Added mean_mrs_list function.
* Improved LCM RAW and BASIS readers (contribution from Alex Craven).

# spant 1.13.0
* ABfit frequency shifts limits are now specified in ppm rather than Hz to
improve consistency between field strengths.
* NAA linewidth is now estimated and output by ABfit when NAAG is absent - 
useful for BRAINO phantom scans.
* Warning now given when spectra are mathematically combined and are not both
in the same time/frequency domain.
* sum_mrs function added to simplify combining spectra within a pipe.
* scale_mrs function added to simplify scaling a spectrum within a pipe.
* add_noise function added to simplify generating simulated data with a pipe.
* hsvd_filt function now accepts a frequency range in ppm units and gives the
option to return the model rather than the filtered data.
* Bug fix for hsvd_filt function where the max_damp argument was ignored.
* SVS reference scans, found in some TWIX files, are now removed by default.
* Bug fix for reading list/data when only one coil element is used.

# spant 1.12.0
* Added 2HG and citrate simulation parameters. e.g. get_mol_paras("2hg").
* Better print output for molecular definitions.
* Added metabolite and basis simulation vignettes.
* Bug fix for setting the ppm reference when reading LCModel RAW files.
* Bug fix for ABfit CRLB calculation of combined signals, eg tNAA, tCr.

# spant 1.11.0
* Options added to allow extra information to attached to mrs_data and
fit_result objects as a data frame.
* New functions (combine_fit_XXX) for working with multiple fit results
contained within a list structure.
* Basis set and noise region are now checked for validity in ABfit.
* Added spant_abfit_benchmark function.
* Tentative functions for performing "standard" 1H brain analyses: 
svs_1h_brain_analysis and svs_1h_brain_batch_analysis.
* Improved support for LCModel analyses.

# spant 1.10.0
* Fix for NIfTI MRS reader/writer.
* ortho3 now shows correct labels for orientations other than RAS.
* ortho3_int function renamed to ortho3_inter.
* Argument order change to plot_voi_overlay and plot_voi_overlay_seg to be more
consistent with ortho3.
* Regression fix for partial volume segmentation plotting.
* Echo time parameter is now stored in the meta structure.
* ABfit now performs a 1D phase parameter search before the prefit stage to
improve reliability. May be disabled with the prefit_phase_search fit option.

# spant 1.9.0
* NIfTI MRS reader and writer now uses the header extension for metadata. Thanks
to Jon Clayden for adding extension read/write support to the RNifti package.
* Default plots now have gridlines in the y-direction and the plot line is now
thicker and coloured blue.
* Opacity option added to the plotting functions (alpha).
* Bug fix for comb_coils with SVS data.
* Bug fixes for Siemens geometry information.
* Changed the ordering of arguments to write_mrs and write_mrs_nifti to improve 
consistency with other functions.
* Internal function (ortho3) now used to plot voxel locations on MRI.
* Internal dicom reader function added.
* Tentative support for Siemens and Philips DICOM MRS format added.
* Internal changes to the way orientation information is handled.

# spant 1.8.0
* Added gridplot function.
* New functions added for down sampling.
* Added signal space projection method for MRSI.
* Geometry information is now read from Siemens twix files.
* GitHub actions are now used for continuous integration instead of Travis and
AppVeyor.
* Added precomp function to avoid repeated computation.
* mrs_data2mat function now collapses all dimensions to dynamics by default.
* mrs_data objects now store the nucleus.
* Added a reader for old Varian format (fid/procpar) data.

# spant 1.7.0
* Added write_mrs function which guesses the output format from the file
extension or can be specified as an argument. write_mrs_XXX functions have been
depreciated.
* read_mrs function now tries to guess the format from the file extension.
* Added json sidecar to NIFTI MRS export function.
* Added the option to read MRS data from a NIFTI file and json sidecar using the
read_mrs function.
* Changed default crop_spec region to between 4.0 and 0.2 ppm.

# spant 1.6.0
* Bug fix for GE P file reader.
* Added downsample_mrs function.
* Tentative function (write_mrs_nifti) to write MRS data as a NIFTI format file 
- for evaluation purposes only.

# spant 1.5.0
* Added an option to ABfit to allow the metabolite amplitudes to be negative 
(ahat_calc_method).
* Removed lsei package dependency.

# spant 1.4.0
* Added the option to plot a y = 0 baseline trace for stackplot.mrs_data.
* Added convenience functions to read and write mrs_data to rds format.

# spant 1.3.0
* Added get_fit_table function to combine all fit tables in a fit_result object
into a single dataframe.
* ppm function can now be applied to fit result objects.
* Bug fix for plot_slice_fit_inter and added the option to specify a
denominator.
* Added the option to specify a denominator to plot_slice_fit.

# spant 1.2.1
* Added the option to display a progress bar in fit_mrs function for better
conformance to "Writing R Extensions" in non-interactive use.
* Changed test tolerance to accommodate differences with OSX.

# spant 1.2.0
* Performance improvement for HSVD water filter.
* ABfit unit tests are now run on simulated data to improve consistency between
different platforms.
* Improvements to fit amplitude scaling code.
* Improved checking for mrs_data processing functions.
* Added preprocessing steps vignette.

# spant 1.1.0
* Improved ppm labels for ABfit plot results.
* Bug fix for plot_slice_fit when using fits from masked data.
* Updated unit tests.
* Added vignette on manually adjusting ABfit baseline smoothness.
* Updated the introduction vignette.
* Bug fixes for image and stackplot functions for masked MRS data.
* Bug fix for using RATS with masked MRS data.
* Changes to prepare for for R 4.0.0.

# spant 1.0.0
* ABfit analysis method has been added, and is now default for the mrs_fit
function.
* Added reader for LCModel RAW format data.
* Added read_mri_dyn_dir function for reading dynamic MRS exported from Siemens
scanners.
* Bug fixes for 2D MRSI voxel segmentation calculation.
* Bug fixes for Siemens IMA format reader for SVS data.
* Optional colourbar added to ortho3 function.

# spant 0.19.0
* Added Asc, BHB, Cho, PEth and Ser simulation parameters.
* Added ker option to get_mrsi_voi function.
* Added append_basis function to combine two basis sets objects.
* The align function now accepts more than one reference frequency.

# spant 0.17.0
* Added a function to grid shift 2D MRSI data in the x/y direction.
* Better plotting/fitting support for masking data by setting data points to NA.
* Bug fix for interactive voxel selection position indicator.
* Added mask_xy to mask voxels in a centred rectangular region.
* Minor changes to improve parallel processing support.

# spant 0.16.0
* SNR is now defined as the max data point divided by the standard deviation of
the noise (n.b. factor of two has been removed in-line with upcoming terminology
paper).
* Default rats method improved to work with multidimensional datasets.
* Added norm_mrs function to normalise the intensity of spectral data.
* Added bc_constant function to correct spectral baselines by a constant offset.
* Added re_weighting function to apply a resolution enhancement weighting to the
FID.
* Performance improvement for apodise_xy function.
* sd function now works for mrs_data.
* Added 2D MRSI support for Siemens IMA format.

# spant 0.15.0
* Bug fix for using auto_phase function with a single spectrum.
* Bug fix for comb_coils not returning unaveraged data when requested.
* Added options to combine metabolite signals from the stackplot of a fit object
and adjust the plot margins.
* Added comb_fits function.
* Added collapse_to_dyns function for mrs_fit objects.
* RDA reader now extracts geometry information.

# spant 0.14.0
* Added options to omit basis signals, change label names and combine lipid and
MM signals from the stackplot of a fit object.
* Added auto_phase function for zeroth order phase-correction of simple spectra.
* Added get_subset function to aid MRSI slicing.
* Added decimate_mrs function.
* Added fit_amps function to quickly extract amplitude estimates from a fit
object.
* Bug fix for int_spec function.
* sim_basis function arguments updated to accept acq_par objects.

# spant 0.13.0
* Various bug fixes for Siemens TWIX reader.
* rats and tdsr functions now use the mean spectrum as the default reference.
* Added the option to remove the x axis in an mrs_data plot.
* Added ylim and y_scale options to fit plotting.
* Added %$% operator from magrittr package.
* Added an interpolation option to calc_spec_snr.
* Added hline and vline options to image.mrs_data.

# spant 0.12.0
* Fit results stackplot now has the option to display labels.
* Added the option to reverse eddy current correction.
* Improved GE p-file reader.
* diff function can now be applied to mrs_data objects.
* Complex functions: Re, Im, Mod, Arg and Conj can now be applied to mrs_data 
objects.
* Default simulations for Glc and NAAG have been improved.

# spant 0.11.0
* Added mar argument to plot command.
* td2fd and fd2td now give warnings when used with data already in the target
domain.
* Improved documentation formatting consistency and fixed some spelling errors.
* Added rats method.

# spant 0.10.0
* The names of in-built pulse sequence functions now all start with seq_* to
make them easier to find.
* Added new functions to simulate the following MRS sequences: CPMG, MEGA-PRESS, 
STEAM, sLASER. sLASER sequence kindly contributed by Pierre-Gilles Henry.
* Bug fix for get_mol_names function.
* stackplot function now accepts labels argument and time-domain plotting.
* def_acq_paras function now accepts arguments to override the defaults.
* Added a source field to mol.paras object to cite the origin of the values.
* Option to restore plotting par defaults.
* The magrittr pipe operator is now exported.

# spant 0.9.0
* Updated plotting modes to be one of : "re", "im", "mod" or "arg".
* Updated int_spec function to use "re", "im", or "mod".
* Added a function to replicate data across a particular dimension.
* Added a convenience function to simulate normal looking 1H brain MRS data.
* phase and shift functions now accept vector inputs.

# spant 0.7.0
* Added new function for frequency drift correction.
* Added support for Siemens ima and TWIX SVS data.
* Added support for GE p-file SVS data.
* Added apply_axes fn.
* Support for reading SPM style segmentation results (spm_pve2categorical).

# spant 0.6.0
* Interactive plotting function added for fit results - plot_fit_slice_inter.
* Bug fix for appending dynamic results.
* Bug fix for reading list data files without reference data.
* Bug fix for append_dyns function.
* basis2mrs_data function has been extended to allow the summation of basis
elements and setting of individual amplitudes.
* Added a shift function for manual frequency shift adjustment.
* Added initial unit tests and automated coveralls checking.

# spant 0.5.0
* A default brain PRESS basis is now simulated by the fit_mrs function when the
basis argument isn't specified.
* Added calc_peak_info function for simple singlet analyses.
* crop_spec function now maintains the original frequency scale.
* The basis set used for analyses has now been added to the fit result object.
* Bug fix for simulating basis sets with one element.
* lb function can now be used with basis-set objects.
* Bug fix for spar_sdat reader for non-localised MRS.
* AppVeyor now being used to test Windows compatibility - John Muschelli.

# spant 0.4.0
* Bug fix for SPAR/SDAT SVS voxel dimensions.
* MRSI support added for Philips SPAR/SDAT data.
* Fit plots now default to the full spectral range unless xlim is specified.
* Fit plots allow the x, y, z, coil, dynamic indices to be specified.
* Added the option to subtract the baseline from fit plots.

# spant 0.3.0
* Added stackplot method for fit objects.
* Added functions for registering and visualising SVS volumes on images and 
performing partial volume correction.
* Philips "list data" also now reads noise scans.
* calc_coil_noise_cor, calc_coil_noise_sd functions added to aid coil 
combination.
* Documentation updates for plotting methods.
* Added some simulation methods to userland.

# spant 0.2.0
* Added Siemens RDA format reader.
* Added Philips "list data" format reader.
* Added Bruker paravision format reader.
* Added PROPACK option for HSVD based filtering.
* Added a coil combination function.
* Bug fix for incorrect ppm scale on fit plots when fs != 2000Hz.
* Bug fix for VARPRO analytical jacobian calculation.

# spant 0.1.0
* First public release.
# Contributing to spant development

Contributions to spant are most welcome and generally take the form of bug reports/feature requests, documentation and suggestions for code changes via a pull request.

## Bug reports/feature requests

Please report bugs and make feature requests via <https://github.com/martin3141/spant/issues/>.

When filing an issue, the most important thing is to include a minimal reproducible example so that we can quickly verify the problem, and then figure out how to fix it.

## Documentation

If you would like to help new users with a particular apect of MRS data processing or analysis please consider contributing a short document in R markdown format to be added to the package. For examples please see the doc folder in the main GitHub repository.

Function reference documentation is derived from roxygen2 code comments and are best contributed via a pull request, as described below.

## Pull requests

To contribute a change to spant, follow these steps:

1. Create a branch in git and make your changes.
1. Push branch to github and issue pull request (PR).
1. Discuss the pull request.
1. Iterate until either we accept the PR or decide that it's not a good fit for spant.

If you're not familiar with git or github, please start by reading <http://r-pkgs.had.co.nz/git.html>. Please use the following style <http://adv-r.had.co.nz/Style.html>.

Please ensure all checks pass with `devtools::check(args = c('--as-cran'))` and consider adding one or more unit tests to the repository tests directory to confirm expected behavior with `devtools::test()`. CI checks are set to run following pushes to the master branch using GitHub Actions, see <https://orchid00.github.io/actions_sandbox> for details.## Test environments

* Windows 7, R 4.0.3
* Linux Ubuntu 20.04, R 4.0.3

## R CMD check results

0 errors | 0 warnings | 1 note

Note - Imports includes 25 non-default packages.

## Downstream dependancies

There are currently no downstream dependencies for this package.---
title: 'spant: An R package for magnetic resonance spectroscopy analysis'
tags:
  - R
  - spectroscopy
  - MRS
  - NMR
  - medical imaging
  - neuroimaging
authors:
  - name: Martin Wilson
    orcid: 0000-0002-2089-3956
    affiliation: 1
affiliations:
 - name: Centre for Human Brain Health and School of Psychology, University of Birmingham, Birmingham, UK
   index: 1
date: 4 August 2021
bibliography: paper.bib
---

# Summary

Magnetic Resonance Spectroscopy (MRS) allows the measurement of small molecules (metabolites) in the body without the use of harmful radiation. Based on the same basic principles and technology behind Magnetic Resonance Imaging (MRI), most modern MRI scanners are also capable of acquiring MRS — making the technique highly suited to a number of clinical applications (@oz:2014). Despite the success of MRS in the research environment, clinical translation has proven slow due to a number of technical and practical reasons, with challenges associated with reliable data processing and analysis having particular importance (@wilson:2019a). The `spant` (SPectroscopy ANalysis Tools) package has been developed to: (1) provide open-source implementations of traditional and modern MRS processing and analysis techniques for routine analysis (@near:2021) and (2) aid the development, validation and comparison of new algorithms and analysis pipelines.

# Statement of need

Traditional MRS analysis was dominated by the use of proprietary software, either supplied by scanner manufactures or offline tools such as LCModel (@provencher:1993) and jMRUI (@naressi:2001). In more recent years there has been a steadily increasing trend toward the use of open-source methods — with some early examples including TARQUIN (@reynolds:2006; @wilson:2011) and AQSES (@poullet:2007). This trend is set to continue with the recent transition of LCModel to an open-source license, and an acceleration in the development of new open-source methods and packages such as Vespa (@soher:2011), Gannet (@edden:2014), FID-A (@simpson:2017), Osprey (@oeltzschner:2020), suspect (@rowland:2021) and FSL-MRS (@clarke:2021a). The availability of the MRSHub (<https://mrshub.org/>), a new community orientated software sharing and support platform, and the recent development of the NIfTI MRS file format (@clarke:2021b), to aid data sharing and interoperability, are set to further enhance the ecosystem of open-source MRS analysis tools.

The vast majority of recently developed open-source MRS analysis tools have been written in either MATLAB or Python. Whilst all languages have strengths and weaknesses, R is particularly suited to the interactive exploration and batch processing of large and complex datasets — typical of MRS and neuroimaging studies. For example, the acquisition and storage of high dimensional datasets, including three spatial axes, chemical shift and coil axes are becoming more common for MRS. For a typical study, these scans may be acquired at multiple time-points for multiple participants split across one or more groups (e.g. control and treatment) - requiring both single subject and group level analyses.

The `spant` package was developed to combine traditional and modern MRS data processing techniques with strengths of R, including: plotting/visualization, statistics, machine learning and data wrangling. Furthermore, `spant` may be used to conveniently combine MRS results with other imaging modalities, due to the availability of a wide range of R packages focused on image processing (@muschelli:2019) and support for the NIfTI data format (@whitcher:2011; @clayden:2021). `spant` also supports the majority of common MR vendor data formats allowing complete pipelines to be developed, from raw time-domain samples to metabolite quantities derived from spectral fitting.

`spant` is part of the MRSHub, demonstrating its acceptance and interest from the MRS community. At the time of writing, `spant` has been used to develop and validate two new MRS spectroscopy analysis algorithms: RATS (@wilson:2019b) and ABfit (@wilson:2021), and has also been used to study cancer (@franco:2021), Alzheimer’s Disease (@montal:2021) and psychosis (@fisher:2020) — confirming its suitability for both MRS methods research and clinical studies.

# Acknowledgements

Particular thanks go to Dr Jonathan D. Clayden and Dr Robert W. Cox for their work on the `RNifti` package (@clayden:2021) and NIfTI standard (<https://nifti.nimh.nih.gov/>) which have substantially expanded the capabilities of `spant`.

# References---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-",
  fig.width = 6,
  fig.height = 5,
  dev = "ragg_png"
)
```

# Spectroscopy Analysis Tools (spant) <img src="man/figures/logo.png" align="right" width=130/>
[![R build status](https://github.com/martin3141/spant/workflows/R-CMD-check/badge.svg)](https://github.com/martin3141/spant/actions)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03646/status.svg)](https://doi.org/10.21105/joss.03646)
[![](http://cranlogs.r-pkg.org/badges/spant)](http://cran.rstudio.com/web/packages/spant/index.html)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/spant)](https://cran.r-project.org/package=spant)
[![Coverage Status](https://coveralls.io/repos/github/martin3141/spant/badge.svg?branch=master)](https://coveralls.io/github/martin3141/spant?branch=master)

## Overview
spant provides a full suite of tools to build automated analysis pipelines for
Magnetic Resonance Spectroscopy (MRS) data. The following features and
algorithms are included:

* Advanced fully-automated metabolite fitting algorithm - ABfit https://onlinelibrary.wiley.com/doi/10.1002/mrm.28385.
* Robust retrospective frequency and phase correction - RATS https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.27605.
* Flexible data types to support single voxel, dynamic and spectroscopic imaging data types.
* Raw data import from individual coils and dynamic measurements, eg support for importing individual FIDs from Siemens TWIX formatted data.
* Publication quality plotting.
* Extensive set of pre-processing steps (phasing, coil-combination, zero-filling, HSVD filtering...)
* Quantum mechanical based simulation for experimental design and basis-set generation.
* Set of metabolite, macromolecule and lipid parameters for typical brain analyses.
* Voxel registration to anatomical images for partial volume concentration corrections.

## Basic installation
Download and install the latest version of R (https://cloud.r-project.org/), or with your package manager if using a recent Linux distribution, eg `sudo apt install r-base`.

It is also strongly recommended to install RStudio Desktop (https://rstudio.com/products/rstudio/download) to provide a modern environment for interactive data analysis.

Once R and RStudio have been installed, open the RStudio application and type the following in the Console (lower left panel) to install the latest stable version of spant:
```{r cran, eval = FALSE}
install.packages("spant", dependencies = TRUE)
```

Or the the development version from GitHub (requires the `devtools` package):
```{r github, eval = FALSE}
install.packages("devtools")
devtools::install_github("martin3141/spant", ref = "devel", dependencies = TRUE)
```

## Documentation
Quick introduction to the basic analysis workflow : https://martin3141.github.io/spant/articles/spant-intro.html

Short tutorials : https://martin3141.github.io/spant/articles/

Function reference : https://martin3141.github.io/spant/reference/

Once the spant library has been loaded with `library(spant)`, type `?spant` on the console for instructions on how to access the offline documentation. Note that offline help on the available functions can be quickly shown in RStudio using `?function_name`, eg `?read_mrs`.

## Ubuntu 20.04 installation
CRAN packages need to be compiled on Linux, and therefore you may need to ensure some additional system libraries are installed. spant may be installed from a clean installation of Ubuntu 20.04 with the following commands pasted into the terminal:

```{ubuntu install, eval = FALSE}
sudo apt install -y r-base libcurl4-openssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libtiff5-dev
mkdir -p ~/R/x86_64-pc-linux-gnu-library/3.6
Rscript -e 'install.packages("spant", dependencies = TRUE)'
```

## Ubuntu 21.10 installation
```{ubuntu 2110 install, eval = FALSE}
sudo apt install -y r-base libcurl4-openssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libtiff5-dev
mkdir -p ~/R/x86_64-pc-linux-gnu-library/4.0
Rscript -e 'install.packages("spant", dependencies = TRUE)'
```

## Anaconda environment installation

Firstly install Anaconda in the standard way for your platform: https://docs.anaconda.com/anaconda/install/index.html.

Create a text file, called `spant_requirements.yml`, containing the following:

```{conda yml, eval = FALSE}
name: spant
channels:
    - conda-forge
    - r
dependencies:
    - r-base
    - r-essentials
    - r-abind
    - r-plyr
    - r-foreach
    - r-pracma
    - r-stringr
    - r-signal
    - r-matrixcalc
    - r-minpack.lm
    - r-nnls
    - r-r.utils
    - r-graphicsqc
    - r-expm
    - r-smoother
    - r-readr
    - r-magrittr
    - r-ptw
    - r-mmand
    - r-RNifti
    - r-RNiftyReg
    - r-fields
    - r-MASS
    - r-numDeriv
    - r-nloptr
    - r-irlba
    - r-tibble
    - r-jsonlite
    - r-viridisLite
    - r-shiny
    - r-miniUI
    - r-knitr
    - r-rmarkdown
    - r-testthat
    - r-ragg
    - r-doParallel
```

Create and activate the environment:

```{conda, eval = FALSE}
conda env create -f spant_requirements.yml
conda activate spant
```

Start R and install spant in the usual way:

```{r conda, eval = FALSE}
install.packages("spant", dependencies = TRUE)
```

Big thanks to [João M.C. Teixeira](https://github.com/joaomcteixeira) for figuring out this method of install.---
output: github_document
---

<!-- index.md is generated from index.Rmd. Please edit that file -->

```{r setup, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-",
  fig.width = 5,
  fig.height = 4
)
```
Spectroscopy Analysis Tools (spant)
=====

[![Travis Build Status](https://travis-ci.org/martin3141/spant.svg?branch=master)](https://travis-ci.org/martin3141/spant) 
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/martin3141/spant?branch=master&svg=true)](https://ci.appveyor.com/project/martin3141/spant)
[![](http://cranlogs.r-pkg.org/badges/spant)](http://cran.rstudio.com/web/packages/spant/index.html)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/spant)](https://cran.r-project.org/package=spant)

## Overview
spant provides a full suite of tools to build automated analysis pipelines for
Magnetic Resonance Spectroscopy (MRS) data. The following features are included:

* Raw data import/export.
* Flexible data types to support single voxel, dynamic and spectroscopic imaging data types.
* Publication quality plotting.
* Extensive set of pre-processing steps (phasing, coil-combination, zero-filling, HSVD filtering...)
* Quantum mechanical based simulation for experimental design and basis-set generation.
* Set of metabolite, macromolecule and lipid parameters for typical brain analyses.
* VARPRO based fitting and interfaces for TARQUIN and LCModel for metabolite quantitation.
* Voxel registration to anatomical images for partial volume concentration corrections.

## Installation
You can install the stable version of spant from CRAN:
```{r cran, eval = FALSE}
install.packages("spant", dependencies = TRUE)
```

Or the the development version from GitHub (requires `devtools` package):
```{r github, eval = FALSE}
install.packages("devtools")
devtools::install_github("martin3141/spant")
```---
title: "spant MEGA-PRESS GABA analysis report"
output: html_document
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(spant)
```

## Edited spectrum
```{r ed_spec, echo = FALSE}
plot(zf(anal_results$ed_gaba), xlim = c(4,0.5))
```

## Frequency drift
```{r freq_drift, echo = FALSE}
plot(res$shifts, xlab = "Dynamic", ylab = "Frequency shift (Hz)")
```

## Edit-off spectrum
```{r ed_off_spec, echo = FALSE}
plot(zf(anal_results$ed_off), xlim = c(4,0.5))
```

## Edit-on spectrum
```{r ed_on_spec, echo = FALSE}
plot(zf(anal_results$ed_on), xlim = c(4,0.5))
```

## Edit-off fit
```{r ed_off_fit, echo = FALSE}
plot(anal_results$ed_off_fit, xlim = c(4,0.5))
```

<style type="text/css">
.table {
    width: 30%;
}
</style>

```{r table, echo = FALSE}
knitr::kable(col.names=c("Amplitude"),t(anal_results$ed_off_fit$results[6:20]),format="markdown")
```---
title: "Common preprocessing steps"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Common preprocessing steps}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
library(ragg)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 5,
  dev = "ragg_png"
)
```

## Reading raw data and plotting
Load the spant package:
```{r, message = FALSE}
library(spant)
```

Load some example data for preprocessing:
```{r}
fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT", package = "spant")
mrs_data <- read_mrs(fname, format = "spar_sdat")
```


Plot the spectral region between 4 and 0.5 ppm:
```{r}
plot(mrs_data, xlim = c(4, 0.5))
```

Apply a 180 degree phase adjustment and plot:
```{r}
mrs_data_p180 <- phase(mrs_data, 180)
plot(mrs_data_p180, xlim = c(4, 0.5))
```

Apply 3 Hz Guassian line broadening:
```{r}
mrs_data_lb <- lb(mrs_data, 3)
plot(mrs_data_lb, xlim = c(4, 0.5))
```

Zero fill the data to twice the original length and plot:
```{r}
mrs_data_zf <- zf(mrs_data, 2)
plot(mrs_data_zf, xlim = c(4, 0.5))
```

Apply a HSVD filter to the residual water region and plot together with the
original data:
```{r}
mrs_data_filt <- hsvd_filt(mrs_data)
stackplot(list(mrs_data, mrs_data_filt), xlim = c(5, 0.5), y_offset = 10,
          col = c("black", "red"), labels = c("original", "filtered"))
```

Apply a 0.1 ppm frequency shift and plot together with the original data:
```{r}
mrs_data_shift <- shift(mrs_data, 0.1, "ppm")
stackplot(list(mrs_data, mrs_data_shift), xlim = c(4, 0.5), y_offset = 10,
          col = c("black", "red"), labels = c("original", "shifted"))
```

Multiple processing commands may be conveniently combined with the pipe operator
"%>%" :
```{r}
mrs_data_proc <- mrs_data %>% hsvd_filt %>% lb(2) %>% zf
plot(mrs_data_proc, xlim = c(5, 0.5))
```---
title: "Basis simulation"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Basis simulation}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
library(ragg)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 5,
  dev = "ragg_png"
)
```

## Basis simulation

Basis simulation is necessary step for modern MRS analysis and the this vignette will explain how to achieve this with spant. It is advisable to follow the examples given in the [metabolite simulation vignette](spant-metabolite-simulation.html) before following this guide.

Load the spant package:

```{r, message = FALSE}
library(spant)
```

A basis set is a collection of signals to be fit to the MRS data. In spant we start with  a list of molecular definitions containing the relevant information for each signal - such as chemical shifts and j-coupling values:

```{r}
mol_list <- list(get_mol_paras("lac"),
                 get_mol_paras("naa"),
                 get_mol_paras("cr"),
                 get_mol_paras("gpc"))
```

In the next step we convert these chemical properties into a collection of signals (a spant `basis_set` object) with the `sim_basis` function. When fitting, the signal parameters (e.g. sampling frequency) and pulse sequence (e.g. echo-time) must match the MRS data acquisition protocol.

```{r}
basis <- sim_basis(mol_list, pul_seq = seq_slaser_ideal,
                   acq_paras = def_acq_paras(N = 2048, fs = 2000, ft = 127.8e6),
                   TE1 = 0.008, TE2 = 0.011, TE3 = 0.009)

stackplot(basis, xlim = c(4, 0.5), y_offset = 50, labels = basis$names)
```

In 1H MRS broad resonances from lipids and macromolecules are often included in addition to metabolites:

```{r}
mol_list_mm <- append(mol_list, list(get_mol_paras("MM09", ft = 127.8e6)))

basis_mm <- sim_basis(mol_list_mm, pul_seq = seq_slaser_ideal,
                   acq_paras = def_acq_paras(N = 2048, fs = 2000, ft = 127.8e6),
                   TE1 = 0.008, TE2 = 0.011, TE3 = 0.009)

stackplot(basis_mm, xlim = c(4, 0.5), y_offset = 50, labels = basis_mm$names)
```

Note the field strength is often required to simulate these broad resonances as their linewidth is usually specified in ppm. spant also includes the functions `sim_basis_1h_brain` and
`sim_basis_1h_brain_press` to produce commonly used sets of basis signals:

```{r}
sim_basis_1h_brain() %>% stackplot(xlim = c(4, 0.5), y_offset = 20, labels = .$names)
```

Basis sets can be exported for use with LCModel with the `write_basis` function, and sim_basis_1h_brain has the option `lcm_compat` to remove signals that are usually generated within the LCModel package:

```{r}
sim_basis_1h_brain(lcm_compat = TRUE) %>%
  stackplot(xlim = c(4, 0.5), y_offset = 20, labels = .$names)
```---
title: "ABfit baseline options"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{ABfit baseline options}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
library(ragg)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 5,
  dev = "ragg_png"
)
```

## Introduction

A good baseline estimate is a prerequisite for accurate metabolite quantitation. The default MRS analysis algorithm in spant (ABfit) is designed to find accurate baseline estimates by automatically adapting the level of baseline flexibility to match the data complexity. This process is inevitably a balancing act between a baseline that is too smooth (resulting in greater bias), and too flexible (resulting in greater variance). Therefore, ABfit has a number of fitting options to adjust the baseline according to user preference. In this vignette the most common adjustments are demonstrated.

## Default analysis

Load the spant analysis package:
```{r setup, message = FALSE}
library(spant)
```
Read an example dataset from file and simulate matching basis set:
```{r read_data}
fname    <- system.file("extdata", "philips_spar_sdat_WS.SDAT", package = "spant")
mrs_data <- read_mrs(fname, format = "spar_sdat")
basis    <- sim_basis_1h_brain_press(mrs_data)
```
Run a default ABfit analysis and plot the result:
```{r abfit_default, results = "hide"}
fit_res <- fit_mrs(mrs_data, basis)
plot(fit_res)
```

The above fit looks good, with a smooth baseline and no significant signals in the residual (top trace) above the noise level. We can find the automatically determined level of baseline smoothness by inspecting the results table in the ``fit_res`` object:
```{r abfit_bl_value}
fit_res$res_tab$bl_ed_pppm
```
The baseline flexibility was found to be `r round(fit_res$res_tab$bl_ed_pppm, 1)` ED per ppm, where ED is the effective dimension -- analogous to the number of spline functions required per ppm. Whilst the automated fit looks reasonable at first glance, let's try and convince ourselves we can't do better with manual adjustments to the algorithm.

## Custom analyses

Changing the default behaviour of ABfit is achieved by supplying an options structure to the ``fit_mrs`` function. The ``abfit_opts`` function generates the default fitting options, which may be modified by supplying arguments. To manually specify the baseline flexibility we set the ``auto_bl_flex`` option to ``FALSE`` and set the ``bl_ed_pppm`` option to the desired level. A greater value results in more baseline flexibility, let's try a value of 8 ED ppm:

```{r abfit_flex, results = "hide"}
opts    <- abfit_opts(auto_bl_flex = FALSE, bl_ed_pppm = 8)
fit_res <- fit_mrs(mrs_data, basis, opts = opts)
plot(fit_res)
```

The baseline is clearly more flexible, resulting in a slightly improved residual, however some baseline features are likely to be due to instability from noise, rather than true spectral features. For the next analysis let's investigate 1 ED pppm:

```{r abfit_stiff, results = "hide"}
opts    <- abfit_opts(auto_bl_flex = FALSE, bl_ed_pppm = 1)
fit_res <- fit_mrs(mrs_data, basis, opts = opts)
plot(fit_res)
```

Now we have a much smoother (almost linear) baseline, which comes at the cost of having broad unmodelled signals in the residual -- ultimately resulting in biased metabolite levels. An alternative to manually specifying a fixed level of baseline flexibility is to adjust the criterion used for automated estimation. The ``aic_smoothing_factor`` can be set to a smaller value (default = 5) to encourage more flexible baselines, whilst still being adaptive to any broad spectral features:

```{r abfit_aic, results = "hide"}
opts    <- abfit_opts(aic_smoothing_factor = 1)
fit_res <- fit_mrs(mrs_data, basis, opts = opts)
plot(fit_res)
```

It can be informative to visualise the individual spline components used for baseline modelling by saving these in the results object:
```{r abfit_bspline, results = "hide"}
opts    <- abfit_opts(export_sp_fit = TRUE)
fit_res <- fit_mrs(mrs_data, basis, opts = opts)
stackplot(fit_res, omit_signals = basis$names)
```

The default number of spline functions for ABfit is 15 per PPM which may be verified from the above plot. Let's try increasing to 25:
```{r abfit_bspline_more, results = "hide"}
opts    <- abfit_opts(export_sp_fit = TRUE, bl_comps_pppm = 25)
fit_res <- fit_mrs(mrs_data, basis, opts = opts)
stackplot(fit_res, omit_signals = basis$names)
```

Clearly the density of spline functions has increased, however the baseline smoothness remains very close to the default. The general principle for ABfit (derived from P-splines) is to over specify the number of baseline modelling spline functions, and rely on a penalty factor to encourage smoothness.
```{r abfit_bl_value_more_splines}
fit_res$res_tab$bl_ed_pppm
```
Inspecting the automatically determined level of baseline flexibility shows the ED per ppm value remains very close to the default analysis despite the change in spline function density -- precisely the desired behaviour.---
title: "Metabolite simulation"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Metabolite simulation}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
library(ragg)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 5,
  dev = "ragg_png"
)
```

## Simple simulation
Load the spant package:
```{r, message = FALSE}
library(spant)
```

Output a list of pre-defined molecules available for simulation:

```{r, message = FALSE}
get_mol_names()
```

Get and print the spin system for myo-inositol:

```{r, message = FALSE}
ins <- get_mol_paras("ins")
print(ins)
```

Simulate and plot the simulation at 7 Tesla for a pulse acquire sequence
(seq_pulse_acquire), apply 2 Hz line-broadening and plot.
```{r, message = FALSE}
sim_mol(ins, ft = 300e6, N = 4096) %>% lb(2) %>% plot(xlim = c(3.8, 3.1))
```

Other pulse sequences may be simulated including: seq_cpmg_ideal,
seq_mega_press_ideal, seq_press_ideal, seq_slaser_ideal, seq_spin_echo_ideal,
seq_steam_ideal. Note all these sequences assume chemical shift displacement is
negligible. Next we simulate a 30 ms spin-echo sequence and plot:

```{r, message = FALSE}
ins_sim <- sim_mol(ins, seq_spin_echo_ideal, ft = 300e6, N = 4086, TE = 0.03)
ins_sim %>% lb(2) %>% plot(xlim = c(3.8, 3.1))
```

Finally we simulate a range of echo-times and plot all results together to see
the phase evolution:

```{r, message = FALSE, fig.height = 8}
sim_fn <- function(TE) {
  te_sim <- sim_mol(ins, seq_spin_echo_ideal, ft = 300e6, N = 4086, TE = TE)
  lb(te_sim, 2)
}

te_vals <- seq(0, 2, 0.4)

lapply(te_vals, sim_fn) %>% stackplot(y_offset = 150, xlim = c(3.8, 3.1),
                                      labels = paste(te_vals * 100, "ms"))
```

See the [basis simulation](spant-basis-simulation.html) vignette for how to combine these simulations into a basis set for MRS analysis.

## Custom molecules

For simple signals that do not require j-coupling evolution, for example
singlets or approximations to macromolecule or lipid resonances, the
`get_uncoupled_mol` function may be used. In this example we simulated two broad
Gaussian resonances at 1.3 and 1.4 ppm with differing amplitudes:

```{r, message = FALSE}
get_uncoupled_mol("Lip13", c(1.3, 1.4), c("1H", "1H"), c(2, 1), c(10, 10),
                  c(1, 1)) %>% sim_mol %>% plot(xlim = c(2, 0.8))
```

Molecules that aren't defined within spant, or need adjusting to match a particular scan, may be manually defined by constructing a `mol_parameters` object. In the following code we define an imaginary molecule based on Lactate, with the addition of a second spin group containing a singlet at 2.5 ppm. Whilst this molecule could be defined as a single group, it is more computationally efficient to split non j-coupled spin systems up in this way. Note the lineshape is set to a Lorentzian (Lorentz-Gauss factor lg = 0) with a width of 2 Hz. It is generally a good idea to simulate resonances with narrower lineshapes that you expect to see in experimental data, as it is far easier to make a resonance broader than narrower.

```{r, message = FALSE}
nucleus_a <- rep("1H", 4)

chem_shift_a <- c(4.0974, 1.3142, 1.3142, 1.3142)

j_coupling_mat_a <- matrix(0, 4, 4)
j_coupling_mat_a[2,1] <- 6.933
j_coupling_mat_a[3,1] <- 6.933
j_coupling_mat_a[4,1] <- 6.933

spin_group_a <- list(nucleus = nucleus_a, chem_shift = chem_shift_a, 
                     j_coupling_mat = j_coupling_mat_a, scale_factor = 1,
                     lw = 2, lg = 0)

nucleus_b <- c("1H")
chem_shift_b <- c(2.5)
j_coupling_mat_b <- matrix(0, 1, 1)

spin_group_b <- list(nucleus = nucleus_b, chem_shift = chem_shift_b, 
                     j_coupling_mat = j_coupling_mat_b, scale_factor = 3,
                     lw = 2, lg = 0)

source <- "This text should include a reference on the origin of the chemical shift and j-coupling values."

custom_mol <- list(spin_groups = list(spin_group_a, spin_group_b), name = "Cus",
              source = source, full_name = "Custom molecule")

class(custom_mol) <- "mol_parameters"
```

In the next step we output the molecule definition as formatted text and plot it.

```{r, message = FALSE}
print(custom_mol)
custom_mol %>% sim_mol %>% lb(2) %>% zf %>% plot(xlim = c(4.4, 0.5))
```


Once your happy the new molecule is correct, please consider contributing it to the package if you think others would benefit.---
title: "Introduction to spant"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to spant}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
library(ragg)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 5,
  dev = "ragg_png"
)
```

## Reading raw data and plotting
Load the spant package:
```{r, message = FALSE}
library(spant)
```

Get the path to a data file included with spant:
```{r}
fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT", package = "spant")
```

Read the file and save to the workspace as ``mrs_data``:
```{r}
mrs_data <- read_mrs(fname, format = "spar_sdat")
```

Output some basic information about the data:
```{r}
print(mrs_data)
```

Plot the spectral region between 5 and 0.5 ppm:
```{r}
plot(mrs_data, xlim = c(5, 0.5))
```

## Basic preprocessing
Apply a HSVD filter to the residual water region and align the spectrum to the tNAA resonance at 2.01 ppm:
```{r}
mrs_proc <- hsvd_filt(mrs_data)
mrs_proc <- align(mrs_proc, 2.01)
plot(mrs_proc, xlim = c(5, 0.5))
```

## Basis simulation
Simulate a typical basis set for short TE brain analysis, print some basic information and plot:
```{r, fig.height=9}
basis <- sim_basis_1h_brain_press(mrs_proc)
print(basis)
stackplot(basis, xlim = c(4, 0.5), labels = basis$names, y_offset = 5)
```

Perform ABfit analysis of the processed data (``mrs_proc``):
```{r, results = "hide"}
fit_res <- fit_mrs(mrs_proc, basis)
```

Plot the fit result: 
```{r}
plot(fit_res)
```

Extract the estimated amplitudes from ``fit_res`` and print as a ratio to total-creatine in column format:
```{r}
amps <- fit_amps(fit_res)
print(t(amps / amps$tCr))
```

Unscaled amplitudes, CRLB error estimates and other fitting diagnostics, such as SNR, are given in the results table:

```{r}
fit_res$res_tab
```

Spectral SNR:
```{r}
fit_res$res_tab$SNR
```

Linewidth of the tNAA resonance in PPM:
```{r}
fit_res$res_tab$tNAA_lw
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_display.R
\name{stackplot.fit_result}
\alias{stackplot.fit_result}
\title{Plot the fitting results of an object of class \code{fit_result} with
individual basis set components shown.}
\usage{
\method{stackplot}{fit_result}(
  x,
  xlim = NULL,
  y_offset = 0,
  dyn = 1,
  x_pos = 1,
  y_pos = 1,
  z_pos = 1,
  coil = 1,
  n = NULL,
  sub_bl = FALSE,
  labels = FALSE,
  label_names = NULL,
  sig_col = "black",
  restore_def_par = TRUE,
  omit_signals = NULL,
  combine_lipmm = FALSE,
  combine_metab = FALSE,
  mar = NULL,
  show_grid = TRUE,
  grid_nx = NULL,
  grid_ny = NA,
  ...
)
}
\arguments{
\item{x}{fit_result object.}

\item{xlim}{the range of values to display on the x-axis, eg xlim = c(4,1).}

\item{y_offset}{separate basis signals in the y-axis direction by this value.}

\item{dyn}{the dynamic index to plot.}

\item{x_pos}{the x index to plot.}

\item{y_pos}{the y index to plot.}

\item{z_pos}{the z index to plot.}

\item{coil}{the coil element number to plot.}

\item{n}{single index element to plot (overrides other indices when given).}

\item{sub_bl}{subtract the baseline from the data and fit (logical).}

\item{labels}{print signal labels at the right side of the plot.}

\item{label_names}{provide a character vector of signal names to replace the
defaults determined from the basis set.}

\item{sig_col}{colour of individual signal components.}

\item{restore_def_par}{restore default plotting par values after the plot has
been made.}

\item{omit_signals}{a character vector of basis signal names to be removed
from the plot.}

\item{combine_lipmm}{combine all basis signals with names starting with "Lip"
or "MM".}

\item{combine_metab}{combine all basis signals with names not starting with
"Lip" or "MM".}

\item{mar}{option to adjust the plot margins. See ?par.}

\item{show_grid}{plot gridlines behind the data (logical). Defaults to TRUE.}

\item{grid_nx}{number of cells of the grid in x and y direction. When NULL
the grid aligns with the tick marks on the corresponding default axis (i.e.,
tickmarks as computed by axTicks). When NA, no grid lines are drawn in the
corresponding direction.}

\item{grid_ny}{as above.}

\item{...}{further arguments to plot method.}
}
\description{
Plot the fitting results of an object of class \code{fit_result} with
individual basis set components shown.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_subset}
\alias{get_subset}
\title{Extract a subset of MRS data.}
\usage{
get_subset(
  mrs_data,
  x_set = NULL,
  y_set = NULL,
  z_set = NULL,
  dyn_set = NULL,
  coil_set = NULL,
  fd_set = NULL,
  td_set = NULL
)
}
\arguments{
\item{mrs_data}{MRS data object.}

\item{x_set}{x indices to include in the output (default all).}

\item{y_set}{y indices to include in the output (default all).}

\item{z_set}{z indices to include in the output (default all).}

\item{dyn_set}{dynamic indices to include in the output (default all).}

\item{coil_set}{coil indices to include in the output (default all).}

\item{fd_set}{frequency domain data indices to include in the output (default
all).}

\item{td_set}{time-domain indices to include in the output (default all).}
}
\value{
selected subset of MRS data.
}
\description{
Extract a subset of MRS data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis_set.R
\name{shift_basis}
\alias{shift_basis}
\title{Apply frequency shifts to basis set signals.}
\usage{
shift_basis(basis, shifts)
}
\arguments{
\item{basis}{the basis to apply the shift to.}

\item{shifts}{a vector of frequency shifts to apply in ppm units. Must be the
same length as there are basis elements.}
}
\value{
modified basis set object.
}
\description{
Apply frequency shifts to basis set signals.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{fs}
\alias{fs}
\title{Return the sampling frequency in Hz of an MRS dataset.}
\usage{
fs(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\value{
sampling frequency in Hz.
}
\description{
Return the sampling frequency in Hz of an MRS dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{acquire}
\alias{acquire}
\title{Simulate pulse sequence acquisition.}
\usage{
acquire(sys, rec_phase = 180, tol = 1e-04, detect = NULL)
}
\arguments{
\item{sys}{spin system object.}

\item{rec_phase}{receiver phase in degrees.}

\item{tol}{ignore resonance amplitudes below this threshold.}

\item{detect}{detection nuclei.}
}
\value{
a list of resonance amplitudes and frequencies.
}
\description{
Simulate pulse sequence acquisition.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{get_1h_brain_basis_paras_v3}
\alias{get_1h_brain_basis_paras_v3}
\title{Return a list of \code{mol_parameter} objects suitable for 1H brain MRS
analyses.}
\usage{
get_1h_brain_basis_paras_v3(ft, metab_lw = NULL, lcm_compat = FALSE)
}
\arguments{
\item{ft}{transmitter frequency in Hz.}

\item{metab_lw}{linewidth of metabolite signals (Hz).}

\item{lcm_compat}{when TRUE, lipid, MM and -CrCH molecules will be excluded
from the output.}
}
\value{
list of \code{mol_parameter} objects.
}
\description{
Return a list of \code{mol_parameter} objects suitable for 1H brain MRS
analyses.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spant.R
\name{def_fs}
\alias{def_fs}
\title{Return the default sampling frequency in Hz.}
\usage{
def_fs()
}
\value{
sampling frequency in Hz.
}
\description{
Return the default sampling frequency in Hz.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{max_mrs}
\alias{max_mrs}
\title{Apply the max operator to an MRS dataset.}
\usage{
max_mrs(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\value{
MRS data following max operator.
}
\description{
Apply the max operator to an MRS dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rats.R
\name{rats}
\alias{rats}
\title{Robust Alignment to a Target Spectrum (RATS).}
\usage{
rats(
  mrs_data,
  ref = NULL,
  xlim = c(4, 0.5),
  max_shift = 20,
  p_deg = 2,
  sp_N = 2,
  sp_deg = 3,
  max_t = 0.2,
  basis_type = "poly",
  rescale_output = TRUE,
  phase_corr = TRUE
)
}
\arguments{
\item{mrs_data}{MRS data to be corrected.}

\item{ref}{optional MRS data to use as a reference, the mean of all dynamics
is used if this argument is not supplied.}

\item{xlim}{optional frequency range to perform optimisation, set to NULL
to use the full range.}

\item{max_shift}{maximum allowable frequency shift in Hz.}

\item{p_deg}{polynomial degree used for baseline modelling. Negative values
disable baseline modelling.}

\item{sp_N}{number of spline functions, note the true number will be sp_N +
sp_deg.}

\item{sp_deg}{degree of spline functions.}

\item{max_t}{truncate the FID when longer than max_t to reduce time taken,
set to NULL to use the entire FID.}

\item{basis_type}{may be one of "poly" or "spline".}

\item{rescale_output}{rescale the bl_matched_spec and bl output to improve
consistency between dynamic scans.}

\item{phase_corr}{apply phase correction (in addition to frequency). TRUE by
default.}
}
\value{
a list containing the corrected data; phase and shift values in units
of degrees and Hz respectively.
}
\description{
Robust Alignment to a Target Spectrum (RATS).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_io.R
\name{read_mrs}
\alias{read_mrs}
\title{Read MRS data from a file.}
\usage{
read_mrs(
  fname,
  format = NULL,
  ft = NULL,
  fs = NULL,
  ref = NULL,
  n_ref_scans = NULL,
  full_fid = FALSE,
  omit_svs_ref_scans = TRUE,
  verbose = FALSE,
  extra = NULL
)
}
\arguments{
\item{fname}{filename of the dpt format MRS data.}

\item{format}{string describing the data format. Must be one of the
following : "spar_sdat", "rda", "dicom", "twix", "pfile", "list_data",
"paravis", "dpt", "lcm_raw", "rds", "nifti", "varian". If not specified,
the format will be guessed from the filename extension.}

\item{ft}{transmitter frequency in Hz (required for list_data format).}

\item{fs}{sampling frequency in Hz (required for list_data format).}

\item{ref}{reference value for ppm scale (required for list_data format).}

\item{n_ref_scans}{override the number of water reference scans detected in
the file header (GE p-file only).}

\item{full_fid}{export all data points, including those before the start
of the FID (default = FALSE), TWIX format only.}

\item{omit_svs_ref_scans}{remove any reference scans sometimes saved in
SVS twix data (default = TRUE).}

\item{verbose}{print data file information (default = FALSE).}

\item{extra}{an optional data frame to provide additional variables for use
in subsequent analysis steps, eg id or grouping variables.}
}
\value{
MRS data object.
}
\description{
Read MRS data from a file.
}
\examples{
fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT", package = "spant")
mrs_data <- read_mrs(fname)
print(mrs_data)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{sim_mol}
\alias{sim_mol}
\title{Simulate a \code{mol_parameter} object.}
\usage{
sim_mol(
  mol,
  pul_seq = seq_pulse_acquire,
  ft = def_ft(),
  ref = def_ref(),
  fs = def_fs(),
  N = def_N(),
  xlim = NULL,
  ...
)
}
\arguments{
\item{mol}{\code{mol_parameter} object.}

\item{pul_seq}{pulse sequence function to use.}

\item{ft}{transmitter frequency in Hz.}

\item{ref}{reference value for ppm scale.}

\item{fs}{sampling frequency in Hz.}

\item{N}{number of data points in the spectral dimension.}

\item{xlim}{ppm range limiting signals to be simulated.}

\item{...}{extra parameters to pass to the pulse sequence function.}
}
\value{
\code{mrs_data} object.
}
\description{
Simulate a \code{mol_parameter} object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{Ncoils}
\alias{Ncoils}
\title{Return the total number of coil elements in an MRS dataset.}
\usage{
Ncoils(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\description{
Return the total number of coil elements in an MRS dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mol_parameters.R
\name{get_uncoupled_mol}
\alias{get_uncoupled_mol}
\title{Generate a \code{mol_parameters} object for a simple spin system with one resonance.}
\usage{
get_uncoupled_mol(
  name,
  chem_shift,
  nucleus,
  scale_factor,
  lw,
  lg,
  full_name = NULL
)
}
\arguments{
\item{name}{abbreviated name of the molecule.}

\item{chem_shift}{chemical shift of the resonance (PPM).}

\item{nucleus}{nucleus (1H, 31P...).}

\item{scale_factor}{multiplicative scaling factor.}

\item{lw}{linewidth in Hz.}

\item{lg}{Lorentz-Gauss lineshape parameter (between 0 and 1).}

\item{full_name}{long name of the molecule (optional).}
}
\value{
mol_parameters object.
}
\description{
Generate a \code{mol_parameters} object for a simple spin system with one resonance.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abfit.R
\name{bbase}
\alias{bbase}
\title{Generate a spline basis, slightly adapted from : "Splines, knots, and
penalties", Eilers 2010.}
\usage{
bbase(N, number, deg = 3)
}
\arguments{
\item{N}{number of data points.}

\item{number}{number of spline functions.}

\item{deg}{spline degree : deg = 1 linear, deg = 2 quadratic, deg = 3 cubic.}
}
\value{
spline basis as a matrix.
}
\description{
Generate a spline basis, slightly adapted from : "Splines, knots, and
penalties", Eilers 2010.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svs_batch_fit.R
\name{svs_1h_brain_batch_analysis}
\alias{svs_1h_brain_batch_analysis}
\title{Batch interface to the standard SVS 1H brain analysis pipeline.}
\usage{
svs_1h_brain_batch_analysis(
  metab_list,
  w_ref_list = NULL,
  mri_seg_list = NULL,
  mri_list = NULL,
  output_dir_list = NULL,
  extra = NULL,
  ...
)
}
\arguments{
\item{metab_list}{list of file paths or mrs_data objects containing MRS
metabolite data.}

\item{w_ref_list}{list of file paths or mrs_data objects containing MRS
water reference data.}

\item{mri_seg_list}{list of file paths or nifti objects containing segmented
MRI data.}

\item{mri_list}{list of file paths or nifti objects containing anatomical
MRI data.}

\item{output_dir_list}{list of directory paths to output fitting results.}

\item{extra}{a data frame with the same number of rows as metab_list,
containing additional information to be attached to the fit results table.}

\item{...}{additional options to be passed to the svs_1h_brain_analysis
function.}
}
\value{
a list of fit_result objects.
}
\description{
Batch interface to the standard SVS 1H brain analysis pipeline.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{est_noise_sd}
\alias{est_noise_sd}
\title{Estimate the standard deviation of the noise from a segment of an mrs_data
object.}
\usage{
est_noise_sd(mrs_data, n = 100, offset = 100, p_order = 2)
}
\arguments{
\item{mrs_data}{MRS data object.}

\item{n}{number of data points (taken from the end of array) to use in the
estimation.}

\item{offset}{number of final points to exclude from the calculation.}

\item{p_order}{polynomial order to fit to the data before estimating the
standard deviation.}
}
\value{
standard deviation array.
}
\description{
Estimate the standard deviation of the noise from a segment of an mrs_data
object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{downsample_mrs_td}
\alias{downsample_mrs_td}
\title{Downsample an MRS signal by a factor of 2 by removing every other data point
in the time-domain. Note, signals outside the new sampling frequency will be
aliased.}
\usage{
downsample_mrs_td(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data object.}
}
\value{
downsampled data.
}
\description{
Downsample an MRS signal by a factor of 2 by removing every other data point
in the time-domain. Note, signals outside the new sampling frequency will be
aliased.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{hsvd_vec}
\alias{hsvd_vec}
\title{HSVD of a complex vector.}
\usage{
hsvd_vec(y, fs, comps = 40, irlba = TRUE, max_damp = 0)
}
\arguments{
\item{y}{time domain signal to be filtered as a vector.}

\item{fs}{sampling frequency of y.}

\item{comps}{number of Lorentzian components to use for modelling.}

\item{irlba}{option to use irlba SVD (logical).}

\item{max_damp}{maximum allowable damping factor. Default value of 0 ensures
resultant model is damped.}
}
\value{
basis matrix and signal table.
}
\description{
HSVD method as described in:
Barkhuijsen H, de Beer R, van Ormondt D. Improved algorithm for noniterative
and timedomain model fitting to exponentially damped magnetic resonance
signals. J Magn Reson 1987;73:553-557.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/precomp.R
\name{set_precomp_verbose}
\alias{set_precomp_verbose}
\title{Set the verbosity of the precompute function.}
\usage{
set_precomp_verbose(verbose = NA)
}
\arguments{
\item{verbose}{can be TRUE or FALSE.}
}
\description{
Set the verbosity of the precompute function.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{spm_pve2categorical}
\alias{spm_pve2categorical}
\title{Convert SPM style segmentation files to a single categorical image where
the numerical values map as: 0) Other, 1) CSF, 2) GM and 3) WM.}
\usage{
spm_pve2categorical(fname)
}
\arguments{
\item{fname}{any of the segmentation files (eg c1_MY_T1.nii).}
}
\value{
nifti object.
}
\description{
Convert SPM style segmentation files to a single categorical image where
the numerical values map as: 0) Other, 1) CSF, 2) GM and 3) WM.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/amp_scaling.R
\name{apply_pvc}
\alias{apply_pvc}
\title{Convert default LCM/TARQUIN concentration scaling to molal units with partial
volume correction.}
\usage{
apply_pvc(fit_result, p_vols, te, tr)
}
\arguments{
\item{fit_result}{a \code{fit_result} object to apply partial volume
correction.}

\item{p_vols}{a numeric vector of partial volumes.}

\item{te}{the MRS TE.}

\item{tr}{the MRS TR.}
}
\value{
a \code{fit_result} object with a rescaled results table.
}
\description{
Convert default LCM/TARQUIN concentration scaling to molal units with partial
volume correction.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{pg_extrap_xy}
\alias{pg_extrap_xy}
\title{Papoulis-Gerchberg (PG) algorithm method for k-space extrapolation.}
\usage{
pg_extrap_xy(
  mrs_data,
  img_mask = NULL,
  kspace_mask = NULL,
  intensity_thresh = 0.15,
  iters = 50
)
}
\arguments{
\item{mrs_data}{MRS data object.}

\item{img_mask}{a boolean matrix of voxels with strong signals to be
extrapolated. Must be twice the dimensions of the input data.}

\item{kspace_mask}{a boolean matrix of kspace points that have been sampled.
Typically a circle for MRSI, but defaults to the full rectangular area of
k-space covered by the input data. Must match the x-y dimensions of the input
data.}

\item{intensity_thresh}{used to define img_mask based on the strength of the
signal in each voxel. Defaults to intensities greater than 15\% of the
maximum. Ignored if img_mask is specified as argument.}

\item{iters}{number of iterations to perform.}
}
\value{
extrapolated \code{mrs_data} object.
}
\description{
PG method as described in: Haupt CI, Schuff N, Weiner MW, Maudsley AA.
Removal of lipid artifacts in 1H spectroscopic imaging by data extrapolation.
Magn Reson Med. 1996 May;35(5):678-87. Extrapolation is performed to expand
k-space coverage by a factor of 2, with the aim to reduce Gibbs ringing.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{Re.mrs_data}
\alias{Re.mrs_data}
\title{Apply Re operator to an MRS dataset.}
\usage{
\method{Re}{mrs_data}(z)
}
\arguments{
\item{z}{MRS data.}
}
\value{
MRS data following Re operator.
}
\description{
Apply Re operator to an MRS dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{sub_mean_dyns}
\alias{sub_mean_dyns}
\title{Subtract the mean dynamic spectrum from a dynamic series.}
\usage{
sub_mean_dyns(mrs_data)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}
}
\value{
subtracted data.
}
\description{
Subtract the mean dynamic spectrum from a dynamic series.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{sim_basis_1h_brain}
\alias{sim_basis_1h_brain}
\title{Simulate a basis-set suitable for 1H brain MRS analysis acquired with a PRESS
sequence. Note, ideal pulses are assumed.}
\usage{
sim_basis_1h_brain(
  pul_seq = seq_press_ideal,
  acq_paras = def_acq_paras(),
  xlim = c(0.5, 4.2),
  lcm_compat = FALSE,
  ...
)
}
\arguments{
\item{pul_seq}{pulse sequence function to use.}

\item{acq_paras}{list of acquisition parameters or an mrs_data object. See
\code{\link{def_acq_paras}}.}

\item{xlim}{range of frequencies to simulate in ppm.}

\item{lcm_compat}{exclude lipid and MM signals for use with default LCModel
options.}

\item{...}{extra parameters to pass to the pulse sequence function.}
}
\value{
basis object.
}
\description{
Simulate a basis-set suitable for 1H brain MRS analysis acquired with a PRESS
sequence. Note, ideal pulses are assumed.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ift_shift_mat}
\alias{ift_shift_mat}
\title{Perform an ifft and ifftshift on a matrix with each column replaced by its
shifted ifft.}
\usage{
ift_shift_mat(mat_in)
}
\arguments{
\item{mat_in}{matrix input.}
}
\value{
output matrix.
}
\description{
Perform an ifft and ifftshift on a matrix with each column replaced by its
shifted ifft.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis_set.R
\name{sort_basis}
\alias{sort_basis}
\title{Sort the basis-set elements alphabetically.}
\usage{
sort_basis(basis)
}
\arguments{
\item{basis}{input basis.}
}
\value{
sorted basis.
}
\description{
Sort the basis-set elements alphabetically.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_io.R
\name{write_mrs}
\alias{write_mrs}
\title{Write MRS data object to file.}
\usage{
write_mrs(mrs_data, fname, format = NULL)
}
\arguments{
\item{mrs_data}{object to be written to file.}

\item{fname}{the filename of the output.}

\item{format}{string describing the data format. Must be one of the
following : "nifti", "dpt", "lcm_raw", "rds". If not specified, the format
will be guessed from the filename extension.}
}
\description{
Write MRS data object to file.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spant.R
\name{def_ref}
\alias{def_ref}
\title{Return the default reference value for ppm scale.}
\usage{
def_ref()
}
\value{
reference value for ppm scale.
}
\description{
Return the default reference value for ppm scale.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abfit.R
\name{abfit_opts}
\alias{abfit_opts}
\title{Return a list of options for an ABfit analysis.}
\usage{
abfit_opts(
  init_damping = 5,
  maxiters = 1024,
  max_shift = 0.078,
  max_damping = 15,
  max_phase = 360,
  lambda = NULL,
  ppm_left = 4,
  ppm_right = 0.2,
  zp = TRUE,
  bl_ed_pppm = 2,
  auto_bl_flex = TRUE,
  bl_comps_pppm = 15,
  export_sp_fit = FALSE,
  max_asym = 0.25,
  max_basis_shift = 0.0078,
  max_basis_damping = 2,
  maxiters_pre = 1000,
  algo_pre = "NLOPT_LN_NELDERMEAD",
  min_bl_ed_pppm = NULL,
  max_bl_ed_pppm = 7,
  auto_bl_flex_n = 20,
  pre_fit_bl_ed_pppm = 1,
  remove_lip_mm_prefit = FALSE,
  pre_align = TRUE,
  max_pre_align_shift = 0.1,
  pre_align_ref_freqs = c(2.01, 3.03, 3.22),
  noise_region = c(-0.5, -2.5),
  optimal_smooth_criterion = "maic",
  aic_smoothing_factor = 5,
  anal_jac = TRUE,
  pre_fit_ppm_left = 4,
  pre_fit_ppm_right = 1.8,
  phi1_optim = FALSE,
  phi1_init = 0,
  max_dphi1 = 0.2,
  max_basis_shift_broad = 0.0078,
  max_basis_damping_broad = 2,
  ahat_calc_method = "lh_pnnls",
  prefit_phase_search = TRUE,
  freq_reg = NULL,
  output_all_paras = FALSE
)
}
\arguments{
\item{init_damping}{initial value of the Gaussian global damping parameter
(Hz). Very poorly shimmed or high field data may benefit from a larger value.}

\item{maxiters}{The maximum number of iterations to run for the detailed fit.}

\item{max_shift}{The maximum allowable shift to be applied in the
optimisation phase of fitting (ppm).}

\item{max_damping}{maximum permitted value of the global damping parameter
(Hz).}

\item{max_phase}{the maximum absolute permitted value of the global
zero-order phase term (degrees). Note, the prefit_phase_search option is not
constrained by this term.}

\item{lambda}{manually set the the baseline smoothness parameter.}

\item{ppm_left}{downfield frequency limit for the fitting range (ppm).}

\item{ppm_right}{upfield frequency limit for the fitting range (ppm).}

\item{zp}{zero pad the data to twice the original length before fitting.}

\item{bl_ed_pppm}{manually set the the baseline smoothness parameter (ED per
ppm).}

\item{auto_bl_flex}{automatically determine the level of baseline smoothness.}

\item{bl_comps_pppm}{spline basis density (signals per ppm).}

\item{export_sp_fit}{add the fitted spline functions to the fit result.}

\item{max_asym}{maximum allowable value of the asymmetry parameter.}

\item{max_basis_shift}{maximum allowable frequency shift for individual basis
signals (ppm).}

\item{max_basis_damping}{maximum allowable Lorentzian damping factor for
individual basis signals (Hz).}

\item{maxiters_pre}{maximum iterations for the coarse (pre-)fit.}

\item{algo_pre}{optimisation method for the coarse (pre-)fit.}

\item{min_bl_ed_pppm}{minimum value for the candidate baseline flexibility
analyses (ED per ppm).}

\item{max_bl_ed_pppm}{minimum value for the candidate baseline flexibility
analyses (ED per ppm).}

\item{auto_bl_flex_n}{number of candidate baseline analyses to perform.}

\item{pre_fit_bl_ed_pppm}{level of baseline flexibility to use in the coarse
fitting stage of the algorithm (ED per ppm).}

\item{remove_lip_mm_prefit}{remove broad signals in the coarse fitting stage
of the algorithm.}

\item{pre_align}{perform a pre-alignment step before coarse fitting.}

\item{max_pre_align_shift}{maximum allowable shift in the pre-alignment step
(ppm).}

\item{pre_align_ref_freqs}{a vector of prominent spectral frequencies used in
the pre-alignment step (ppm).}

\item{noise_region}{spectral region to estimate the noise level (ppm).}

\item{optimal_smooth_criterion}{method to determine the optimal smoothness.}

\item{aic_smoothing_factor}{modification factor for the AIC calculation.}

\item{anal_jac}{use a analytical approximation to the jacobian in the
detailed fitting stage.}

\item{pre_fit_ppm_left}{downfield frequency limit for the fitting range in
the coarse fitting stage of the algorithm (ppm).}

\item{pre_fit_ppm_right}{upfield frequency limit for the fitting range in the
coarse fitting stage of the algorithm (ppm).}

\item{phi1_optim}{apply and optimise a frequency dependant phase term.}

\item{phi1_init}{initial value for the frequency dependant phase term (ms).}

\item{max_dphi1}{maximum allowable change from the initial frequency
dependant phase term (ms).}

\item{max_basis_shift_broad}{maximum allowable shift for broad signals in the
basis (ppm). Determined based on their name beginning with Lip or MM.}

\item{max_basis_damping_broad}{maximum allowable Lorentzian damping for broad
signals in the basis (Hz). Determined based on their name beginning with Lip
or MM.}

\item{ahat_calc_method}{method to calculate the metabolite amplitudes. May be
one of: "lh_pnnls" or "ls".}

\item{prefit_phase_search}{perform a 1D search for the optimal phase in the
prefit stage of the algorithm.}

\item{freq_reg}{frequency shift parameter.}

\item{output_all_paras}{include more fitting parameters in the fit table,
e.g. individual shift and damping factors for each basis set element.}
}
\value{
full list of options.
}
\description{
Return a list of options for an ABfit analysis.
}
\examples{
opts <- abfit_opts(ppm_left = 4.2, noise_region = c(-1, -3))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis_set.R
\name{write_basis_tqn}
\alias{write_basis_tqn}
\title{Generate a basis file using TARQUIN.}
\usage{
write_basis_tqn(basis_file, metab_data, opts = NULL)
}
\arguments{
\item{basis_file}{filename of the basis file to be generated.}

\item{metab_data}{MRS data object to match the generated basis parameters.}

\item{opts}{list of options to pass to TARQUIN.}
}
\description{
Generate a basis file using TARQUIN.
}
\examples{
\dontrun{
write_basis_tqn('test.basis',mrs_data,c("--echo","0.04"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spant.R
\name{set_lcm_cmd}
\alias{set_lcm_cmd}
\title{Set the command to run the LCModel command-line program.}
\usage{
set_lcm_cmd(cmd)
}
\arguments{
\item{cmd}{path to binary.}
}
\description{
Set the command to run the LCModel command-line program.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interactive_plotting.R
\name{ortho3_inter}
\alias{ortho3_inter}
\title{Display an interactive orthographic projection plot of a nifti object.}
\usage{
ortho3_inter(
  underlay,
  overlay = NULL,
  xyz = NULL,
  zlim = NULL,
  zlim_ol = NULL,
  alpha = 0.7,
  ...
)
}
\arguments{
\item{underlay}{underlay image to be shown in grayscale.}

\item{overlay}{optional overlay image.}

\item{xyz}{x, y, z slice coordinates to display.}

\item{zlim}{underlay intensity limits.}

\item{zlim_ol}{overlay intensity limits.}

\item{alpha}{transparency of overlay.}

\item{...}{other options to be passed to the ortho3 function.}
}
\description{
Display an interactive orthographic projection plot of a nifti object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{fd2td}
\alias{fd2td}
\title{Transform frequency-domain data to the time-domain.}
\usage{
fd2td(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data in frequency-domain representation.}
}
\value{
MRS data in time-domain representation.
}
\description{
Transform frequency-domain data to the time-domain.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{append_coils}
\alias{append_coils}
\title{Append MRS data across the coil dimension, assumes they matched across the
other dimensions.}
\usage{
append_coils(...)
}
\arguments{
\item{...}{MRS data objects as arguments, or a list of MRS data objects.}
}
\value{
a single MRS data object with the input objects concatenated together.
}
\description{
Append MRS data across the coil dimension, assumes they matched across the
other dimensions.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{Ndyns}
\alias{Ndyns}
\title{Return the total number of dynamic scans in an MRS dataset.}
\usage{
Ndyns(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\description{
Return the total number of dynamic scans in an MRS dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{apply_axes}
\alias{apply_axes}
\title{Apply a function over specified array axes.}
\usage{
apply_axes(x, axes, fun, ...)
}
\arguments{
\item{x}{an array.}

\item{axes}{a vector of axes to apply fun over.}

\item{fun}{function to be applied.}

\item{...}{optional arguments to fun.}
}
\value{
array.
}
\description{
Apply a function over specified array axes.
}
\examples{
z <- array(1:1000, dim = c(10, 10, 10))
a <- apply_axes(z, 3, fft)
a[1,1,] == fft(z[1,1,])
a <- apply_axes(z, 3, sum)
a[1,1,] == sum(z[1,1,])
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{read_tqn_result}
\alias{read_tqn_result}
\title{Reader for csv results generated by TARQUIN.}
\usage{
read_tqn_result(result_f, remove_rcs = TRUE)
}
\arguments{
\item{result_f}{TARQUIN result file.}

\item{remove_rcs}{omit row, column and slice ids from output.}
}
\value{
list of amplitudes, crlbs and diagnostics.
}
\description{
Reader for csv results generated by TARQUIN.
}
\examples{
\dontrun{
result <- read_tqn_result(system.file("extdata","result.csv",package="spant"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{seconds}
\alias{seconds}
\title{Return a time scale vector to match the FID of an MRS data object.}
\usage{
seconds(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\value{
time scale vector in units of seconds.
}
\description{
Return a time scale vector to match the FID of an MRS data object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_fp}
\alias{get_fp}
\title{Return the first time-domain data point.}
\usage{
get_fp(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\value{
first time-domain data point.
}
\description{
Return the first time-domain data point.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_display.R
\name{get_fit_map}
\alias{get_fit_map}
\title{Get a data array from a fit result.}
\usage{
get_fit_map(fit_res, name)
}
\arguments{
\item{fit_res}{\code{fit_result} object.}

\item{name}{name of the quantity to plot, eg "tNAA".}
}
\description{
Get a data array from a fit result.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{Nx}
\alias{Nx}
\title{Return the total number of x locations in an MRS dataset.}
\usage{
Nx(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\description{
Return the total number of x locations in an MRS dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_write_nifti.R
\name{write_mrs_nifti}
\alias{write_mrs_nifti}
\title{Write MRS data object to file in NIFTI format.}
\usage{
write_mrs_nifti(mrs_data, fname)
}
\arguments{
\item{mrs_data}{object to be written to file.}

\item{fname}{the filename of the output NIFTI MRS data.}
}
\description{
Write MRS data object to file in NIFTI format.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spant.R
\name{set_def_acq_paras}
\alias{set_def_acq_paras}
\title{Set the default acquisition parameters.}
\usage{
set_def_acq_paras(
  ft = getOption("spant.def_ft"),
  fs = getOption("spant.def_fs"),
  N = getOption("spant.def_N"),
  ref = getOption("spant.def_ref"),
  nuc = getOption("spant.nuc")
)
}
\arguments{
\item{ft}{transmitter frequency in Hz.}

\item{fs}{sampling frequency in Hz.}

\item{N}{number of data points in the spectral dimension.}

\item{ref}{reference value for ppm scale.}

\item{nuc}{resonant nucleus.}
}
\description{
Set the default acquisition parameters.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{int_spec}
\alias{int_spec}
\title{Integrate a spectral region.}
\usage{
int_spec(mrs_data, xlim = NULL, freq_scale = "ppm", mode = "re")
}
\arguments{
\item{mrs_data}{MRS data.}

\item{xlim}{spectral range to be integrated (defaults to full range).}

\item{freq_scale}{units of xlim, can be : "ppm", "hz" or "points".}

\item{mode}{spectral mode, can be : "re", "im", "mod" or "cplx".}
}
\value{
an array of integral values.
}
\description{
See spec_op function for a more complete set of spectral operations.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pulse_sequences.R
\name{seq_pulse_acquire}
\alias{seq_pulse_acquire}
\title{Simple pulse and acquire sequence with ideal pulses.}
\usage{
seq_pulse_acquire(spin_params, ft, ref)
}
\arguments{
\item{spin_params}{spin system definition.}

\item{ft}{transmitter frequency in Hz.}

\item{ref}{reference value for ppm scale.}
}
\value{
list of resonance amplitudes and frequencies.
}
\description{
Simple pulse and acquire sequence with ideal pulses.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{rm_dyns}
\alias{rm_dyns}
\title{Remove a subset of dynamic scans.}
\usage{
rm_dyns(mrs_data, subset)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}

\item{subset}{vector containing indices to the dynamic scans to be
removed.}
}
\value{
MRS data without the specified dynamic scans.
}
\description{
Remove a subset of dynamic scans.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{td2fd}
\alias{td2fd}
\title{Transform time-domain data to the frequency-domain.}
\usage{
td2fd(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data in time-domain representation.}
}
\value{
MRS data in frequency-domain representation.
}
\description{
Transform time-domain data to the frequency-domain.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{sim_brain_1h}
\alias{sim_brain_1h}
\title{Simulate MRS data with a similar appearance to normal brain (by default).}
\usage{
sim_brain_1h(
  acq_paras = def_acq_paras(),
  type = "normal_v1",
  pul_seq = seq_press_ideal,
  xlim = c(0.5, 4.2),
  full_output = FALSE,
  amps = NULL,
  ...
)
}
\arguments{
\item{acq_paras}{list of acquisition parameters or an mrs_data object. See
\code{\link{def_acq_paras}}.}

\item{type}{type of spectrum, only "normal" is implemented currently.}

\item{pul_seq}{pulse sequence function to use.}

\item{xlim}{range of frequencies to simulate in ppm.}

\item{full_output}{when FALSE (default) only output the simulated MRS data.
When TRUE output a list containing the MRS data, basis set object and
corresponding amplitudes.}

\item{amps}{a vector of basis amplitudes may be specified to modify the
output spectrum.}

\item{...}{extra parameters to pass to the pulse sequence function.}
}
\value{
see full_output option.
}
\description{
Simulate MRS data with a similar appearance to normal brain (by default).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{gen_F_xy}
\alias{gen_F_xy}
\title{Generate the Fxy product operator with a specified phase.}
\usage{
gen_F_xy(sys, phase, detect = NULL)
}
\arguments{
\item{sys}{spin system object.}

\item{phase}{phase angle in degrees.}

\item{detect}{detection nuclei.}
}
\value{
product operator matrix.
}
\description{
Generate the Fxy product operator with a specified phase.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{lb}
\alias{lb}
\alias{lb.mrs_data}
\alias{lb.basis_set}
\title{Apply line-broadening (apodisation) to MRS data or basis object.}
\usage{
lb(x, lb, lg = 1)

\method{lb}{mrs_data}(x, lb, lg = 1)

\method{lb}{basis_set}(x, lb, lg = 1)
}
\arguments{
\item{x}{input mrs_data or basis_set object.}

\item{lb}{amount of line-broadening in Hz.}

\item{lg}{Lorentz-Gauss lineshape parameter (between 0 and 1).}
}
\value{
line-broadened data.
}
\description{
Apply line-broadening (apodisation) to MRS data or basis object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{kspace2img_xy}
\alias{kspace2img_xy}
\title{Transform 2D MRSI data from k-space to image space in the x-y direction.}
\usage{
kspace2img_xy(mrs_data)
}
\arguments{
\item{mrs_data}{2D MRSI data.}
}
\value{
MRSI data in image space.
}
\description{
Transform 2D MRSI data from k-space to image space in the x-y direction.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abfit.R
\name{calc_ed_from_lambda}
\alias{calc_ed_from_lambda}
\title{Calculate the effective dimensions of a spline smoother from lambda.}
\usage{
calc_ed_from_lambda(spline_basis, deriv_mat, lambda)
}
\arguments{
\item{spline_basis}{spline basis.}

\item{deriv_mat}{derivative matrix.}

\item{lambda}{smoothing parameter.}
}
\value{
the effective dimension value.
}
\description{
Calculate the effective dimensions of a spline smoother from lambda.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{calc_spec_snr}
\alias{calc_spec_snr}
\title{Calculate the spectral SNR.}
\usage{
calc_spec_snr(
  mrs_data,
  sig_region = c(4, 0.5),
  noise_region = c(-0.5, -2.5),
  p_order = 2,
  interp_f = 4,
  full_output = FALSE
)
}
\arguments{
\item{mrs_data}{an object of class \code{mrs_data}.}

\item{sig_region}{a ppm region to define where the maximum signal value
should be estimated.}

\item{noise_region}{a ppm region to defined where the noise level should be
estimated.}

\item{p_order}{polynomial order to fit to the noise region before estimating
the standard deviation.}

\item{interp_f}{interpolation factor to improve detection of the highest
signal value.}

\item{full_output}{output signal, noise and SNR values separately.}
}
\value{
an array of SNR values.
}
\description{
SNR is defined as the maximum signal value divided by the standard deviation
of the noise.
}
\details{
The mean noise value is subtracted from the maximum signal value to reduce DC
offset bias. A polynomial detrending fit (second order by default) is applied
to the noise region before the noise standard deviation is estimated.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{ppm}
\alias{ppm}
\alias{ppm.mrs_data}
\alias{ppm.fit_result}
\title{Return the ppm scale of an MRS dataset or fit result.}
\usage{
ppm(x, ft = NULL, ref = NULL, fs = NULL, N = NULL)

\method{ppm}{mrs_data}(x, ft = NULL, ref = NULL, fs = NULL, N = NULL)

\method{ppm}{fit_result}(x, ft = NULL, ref = NULL, fs = NULL, N = NULL)
}
\arguments{
\item{x}{MRS dataset of fit result.}

\item{ft}{transmitter frequency in Hz, does not apply when the object is a
fit result.}

\item{ref}{reference value for ppm scale, does not apply when the object is a
fit result.}

\item{fs}{sampling frequency in Hz, does not apply when the object is a
fit result.}

\item{N}{number of data points in the spectral dimension, does not apply when the object is a
fit result.}
}
\value{
ppm scale.
}
\description{
Return the ppm scale of an MRS dataset or fit result.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{rep_dyn}
\alias{rep_dyn}
\title{Replicate a scan in the dynamic dimension.}
\usage{
rep_dyn(mrs_data, times)
}
\arguments{
\item{mrs_data}{MRS data to be replicated.}

\item{times}{number of times to replicate.}
}
\value{
replicated data object.
}
\description{
Replicate a scan in the dynamic dimension.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{inv_even_dyns}
\alias{inv_even_dyns}
\title{Invert even numbered dynamic scans starting from 1 (2,4,6...).}
\usage{
inv_even_dyns(mrs_data)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}
}
\value{
dynamic MRS data with inverted even numbered scans.
}
\description{
Invert even numbered dynamic scans starting from 1 (2,4,6...).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{mean_dyns}
\alias{mean_dyns}
\title{Calculate the mean dynamic data.}
\usage{
mean_dyns(mrs_data)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}
}
\value{
mean dynamic data.
}
\description{
Calculate the mean dynamic data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{fit_mrs}
\alias{fit_mrs}
\title{Perform a fit based analysis of MRS data.}
\usage{
fit_mrs(
  metab,
  basis = NULL,
  method = "ABFIT",
  w_ref = NULL,
  opts = NULL,
  parallel = FALSE,
  time = TRUE,
  progress = "text",
  extra = metab$extra
)
}
\arguments{
\item{metab}{metabolite data.}

\item{basis}{basis class object or character vector to basis file in
LCModel .basis format.}

\item{method}{'ABFIT' (default), 'VARPRO', 'VARPRO_3P', 'TARQUIN' or
'LCMODEL'.}

\item{w_ref}{water reference data for concentration scaling (optional).}

\item{opts}{options to pass to the analysis method.}

\item{parallel}{perform analyses in parallel (TRUE or FALSE).}

\item{time}{measure the time taken for the analysis to complete
(TRUE or FALSE).}

\item{progress}{option is passed to plyr::alply function to display a
progress bar during fitting. Default value is "text", set to "none" to
disable.}

\item{extra}{an optional data frame to provide additional variables for use
in subsequent analysis steps, eg id or grouping variables.}
}
\value{
MRS analysis object.
}
\description{
Note that TARQUIN and LCModel require these packages to be installed, and
the functions set_tqn_cmd and set_lcm_cmd (respectively) need to be used to
specify the location of these software packages.
}
\details{
Fitting approaches described in the following references:
ABfit
Wilson, M. Adaptive baseline fitting for 1H MR spectroscopy analysis. Magn
Reson Med 2012;85:13-29.

VARPRO
van der Veen JW, de Beer R, Luyten PR, van Ormondt D. Accurate quantification
of in vivo 31P NMR signals using the variable projection method and prior
knowledge. Magn Reson Med 1988;6:92-98.

TARQUIN
Wilson, M., Reynolds, G., Kauppinen, R. A., Arvanitis, T. N. & Peet, A. C.
A constrained least-squares approach to the automated quantitation of in vivo
1H magnetic resonance spectroscopy data. Magn Reson Med 2011;65:1-12.

LCModel
Provencher SW. Estimation of metabolite concentrations from localized in vivo
proton NMR spectra. Magn Reson Med 1993;30:672-679.
}
\examples{
fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT", package =
"spant")
svs <- read_mrs(fname)
\dontrun{
basis <- sim_basis_1h_brain_press(svs)
fit_result <- fit_mrs(svs, basis)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/amp_scaling.R
\name{scale_amp_ratio}
\alias{scale_amp_ratio}
\title{Scale fitted amplitudes to a ratio of signal amplitude.}
\usage{
scale_amp_ratio(fit_result, name)
}
\arguments{
\item{fit_result}{a result object generated from fitting.}

\item{name}{the signal name to use as a denominator (usually, "tCr" or
"tNAA").}
}
\value{
a \code{fit_result} object with a rescaled results table.
}
\description{
Scale fitted amplitudes to a ratio of signal amplitude.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{conv_mrs}
\alias{conv_mrs}
\title{Convolve two MRS data objects.}
\usage{
conv_mrs(mrs_data, conv)
}
\arguments{
\item{mrs_data}{MRS data to be convolved.}

\item{conv}{convolution data stored as an mrs_data object.}
}
\value{
convolved data.
}
\description{
Convolve two MRS data objects.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_tail_dyns}
\alias{get_tail_dyns}
\title{Return the last scans of a dynamic series.}
\usage{
get_tail_dyns(mrs_data, n = 1)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}

\item{n}{the number of dynamic scans to return.}
}
\value{
last scans of a dynamic series.
}
\description{
Return the last scans of a dynamic series.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_display.R
\name{plot_bc}
\alias{plot_bc}
\title{Convenience function to plot a baseline estimate with the original data.}
\usage{
plot_bc(orig_data, bc_data, ...)
}
\arguments{
\item{orig_data}{the original data.}

\item{bc_data}{the baseline corrected data.}

\item{...}{other arguments to pass to the stackplot function.}
}
\description{
Convenience function to plot a baseline estimate with the original data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{back_extrap_ar}
\alias{back_extrap_ar}
\title{Back extrapolate time-domain data points using an autoregressive model.}
\usage{
back_extrap_ar(
  mrs_data,
  extrap_pts,
  pred_pts = NULL,
  method = "burg",
  rem_add = TRUE,
  ...
)
}
\arguments{
\item{mrs_data}{mrs_data object.}

\item{extrap_pts}{number of points to extrapolate.}

\item{pred_pts}{number of points to base the extrapolation on.}

\item{method}{character string specifying the method to fit the model. Must
be one of the strings in the default argument (the first few characters are
sufficient). Defaults to "burg".}

\item{rem_add}{remove additional points from the end of the FID to maintain
the original length of the dataset. Default to TRUE.}

\item{...}{additional arguments to specific methods, see ?ar.}
}
\value{
back extrapolated data.
}
\description{
Back extrapolate time-domain data points using an autoregressive model.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/amp_scaling.R
\name{scale_amp_water_ratio}
\alias{scale_amp_water_ratio}
\title{Scale metabolite amplitudes as a ratio to the unsuppressed water amplitude.}
\usage{
scale_amp_water_ratio(fit_result, ref_data, ...)
}
\arguments{
\item{fit_result}{a result object generated from fitting.}

\item{ref_data}{a water reference MRS data object.}

\item{...}{additional arguments to get_td_amp function.}
}
\value{
a \code{fit_result} object with a rescaled results table.
}
\description{
Scale metabolite amplitudes as a ratio to the unsuppressed water amplitude.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ft_shift}
\alias{ft_shift}
\title{Perform a fft and ffshift on a vector.}
\usage{
ft_shift(vec_in)
}
\arguments{
\item{vec_in}{vector input.}
}
\value{
output vector.
}
\description{
Perform a fft and ffshift on a vector.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_display.R
\name{print.mrs_data}
\alias{print.mrs_data}
\title{Print a summary of mrs_data parameters.}
\usage{
\method{print}{mrs_data}(x, full = FALSE, ...)
}
\arguments{
\item{x}{mrs_data object.}

\item{full}{print all parameters (default FALSE).}

\item{...}{further arguments.}
}
\description{
Print a summary of mrs_data parameters.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spant.R
\name{check_tqn}
\alias{check_tqn}
\title{Check the TARQUIN binary can be run}
\usage{
check_tqn()
}
\description{
Check the TARQUIN binary can be run
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{align}
\alias{align}
\title{Align spectra to a reference frequency using a convolution based method.}
\usage{
align(
  mrs_data,
  ref_freq = 4.65,
  zf_factor = 2,
  lb = 2,
  max_shift = 20,
  ret_df = FALSE,
  mean_dyns = FALSE
)
}
\arguments{
\item{mrs_data}{data to be aligned.}

\item{ref_freq}{reference frequency in ppm units. More than one frequency
may be specified.}

\item{zf_factor}{zero filling factor to increase alignment resolution.}

\item{lb}{line broadening to apply to the reference signal.}

\item{max_shift}{maximum allowable shift in Hz.}

\item{ret_df}{return frequency shifts in addition to aligned data (logical).}

\item{mean_dyns}{align the mean spectrum and apply the same shift to each
dynamic.}
}
\value{
aligned data object.
}
\description{
Align spectra to a reference frequency using a convolution based method.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis_set.R
\name{read_basis_ac}
\alias{read_basis_ac}
\title{Read a basis file in LCModel .basis format (for testing only).}
\usage{
read_basis_ac(basis_file, ref = def_ref(), sort_basis = TRUE)
}
\arguments{
\item{basis_file}{path to basis file.}

\item{ref}{assumed ppm reference value.}

\item{sort_basis}{sort the basis set based on signal names.}
}
\value{
basis object.
}
\description{
Read a basis file in LCModel .basis format (for testing only).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{decimate_mrs_fd}
\alias{decimate_mrs_fd}
\title{Decimate an MRS signal to half the original sampling frequency by filtering
in the frequency domain before down sampling.}
\usage{
decimate_mrs_fd(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data object.}
}
\value{
decimated data at half the original sampling frequency.
}
\description{
Decimate an MRS signal to half the original sampling frequency by filtering
in the frequency domain before down sampling.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spant.R
\name{def_nuc}
\alias{def_nuc}
\title{Return the default nucleus.}
\usage{
def_nuc()
}
\value{
number of data points in the spectral dimension.
}
\description{
Return the default nucleus.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{spant_mpress_drift}
\alias{spant_mpress_drift}
\title{Example MEGA-PRESS data with significant B0 drift.}
\format{
An object of class \code{mrs_data} of length 13.
}
\usage{
spant_mpress_drift
}
\description{
Example MEGA-PRESS data with significant B0 drift.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{zero_nzoc}
\alias{zero_nzoc}
\title{Zero all non-zero-order coherences.}
\usage{
zero_nzoc(sys, rho)
}
\arguments{
\item{sys}{spin system object.}

\item{rho}{density matrix.}
}
\value{
density matrix.
}
\description{
Zero all non-zero-order coherences.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spant.R
\name{get_tqn_cmd}
\alias{get_tqn_cmd}
\title{Print the command to run the TARQUIN command-line program.}
\usage{
get_tqn_cmd()
}
\description{
Print the command to run the TARQUIN command-line program.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{fd_conv_filt}
\alias{fd_conv_filt}
\title{Frequency-domain convolution based filter.}
\usage{
fd_conv_filt(mrs_data, K = 25, ext = 1)
}
\arguments{
\item{mrs_data}{MRS data to be filtered.}

\item{K}{window width in data points.}

\item{ext}{point separation for linear extrapolation.}
}
\description{
Frequency-domain convolution based filter.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis_set.R
\name{sim_basis_tqn}
\alias{sim_basis_tqn}
\title{Simulate a basis file using TARQUIN.}
\usage{
sim_basis_tqn(
  fs = def_fs(),
  ft = def_ft(),
  N = def_N(),
  ref = def_ref(),
  opts = NULL
)
}
\arguments{
\item{fs}{sampling frequency}

\item{ft}{transmitter frequency}

\item{N}{number of data points}

\item{ref}{chemical shift reference}

\item{opts}{list of options to pass to TARQUIN.}
}
\description{
Simulate a basis file using TARQUIN.
}
\examples{
\dontrun{
write_basis_tqn('test.basis',mrs_data,c("--echo","0.04"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_display.R
\name{fit_res2csv}
\alias{fit_res2csv}
\title{Write fit results table to a csv file.}
\usage{
fit_res2csv(fit_res, fname, unscaled = FALSE)
}
\arguments{
\item{fit_res}{fit result object.}

\item{fname}{filename of csv file.}

\item{unscaled}{output the unscaled result table (default = FALSE).}
}
\description{
Write fit results table to a csv file.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{get_1h_brain_basis_paras_v2}
\alias{get_1h_brain_basis_paras_v2}
\title{Return a list of \code{mol_parameter} objects suitable for 1H brain MRS
analyses.}
\usage{
get_1h_brain_basis_paras_v2(ft, metab_lw = NULL, lcm_compat = FALSE)
}
\arguments{
\item{ft}{transmitter frequency in Hz.}

\item{metab_lw}{linewidth of metabolite signals (Hz).}

\item{lcm_compat}{when TRUE, lipid, MM and -CrCH molecules will be excluded
from the output.}
}
\value{
list of \code{mol_parameter} objects.
}
\description{
Return a list of \code{mol_parameter} objects suitable for 1H brain MRS
analyses.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{sum_dyns}
\alias{sum_dyns}
\title{Calculate the sum of data dynamics.}
\usage{
sum_dyns(mrs_data)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}
}
\value{
sum of data dynamics.
}
\description{
Calculate the sum of data dynamics.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{crop_xy}
\alias{crop_xy}
\title{Crop an MRSI dataset in the x-y direction}
\usage{
crop_xy(mrs_data, x_dim, y_dim)
}
\arguments{
\item{mrs_data}{MRS data object.}

\item{x_dim}{x dimension output length.}

\item{y_dim}{y dimension output length.}
}
\value{
selected subset of MRS data.
}
\description{
Crop an MRSI dataset in the x-y direction
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{rep_mrs}
\alias{rep_mrs}
\title{Replicate a scan over a given dimension.}
\usage{
rep_mrs(
  mrs_data,
  x_rep = 1,
  y_rep = 1,
  z_rep = 1,
  dyn_rep = 1,
  coil_rep = 1,
  warn = TRUE
)
}
\arguments{
\item{mrs_data}{MRS data to be replicated.}

\item{x_rep}{number of x replications.}

\item{y_rep}{number of y replications.}

\item{z_rep}{number of z replications.}

\item{dyn_rep}{number of dynamic replications.}

\item{coil_rep}{number of coil replications.}

\item{warn}{print a warning when the data dimensions do not change.}
}
\value{
replicated data object.
}
\description{
Replicate a scan over a given dimension.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{auto_phase}
\alias{auto_phase}
\title{Perform zeroth-order phase correction based on the minimisation of the
squared difference between the real and magnitude components of the
spectrum.}
\usage{
auto_phase(mrs_data, xlim = NULL, ret_phase = FALSE)
}
\arguments{
\item{mrs_data}{an object of class \code{mrs_data}.}

\item{xlim}{frequency range (default units of PPM) to including in the phase.}

\item{ret_phase}{return phase values (logical).}
}
\value{
MRS data object and phase values (optional).
}
\description{
Perform zeroth-order phase correction based on the minimisation of the
squared difference between the real and magnitude components of the
spectrum.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{calc_sd_poly}
\alias{calc_sd_poly}
\title{Perform a polynomial fit, subtract and return the standard deviation of the
residuals.}
\usage{
calc_sd_poly(y, degree = 1)
}
\arguments{
\item{y}{array.}

\item{degree}{polynomial degree.}
}
\value{
standard deviation of the fit residuals.
}
\description{
Perform a polynomial fit, subtract and return the standard deviation of the
residuals.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{Arg.mrs_data}
\alias{Arg.mrs_data}
\title{Apply Arg operator to an MRS dataset.}
\usage{
\method{Arg}{mrs_data}(z)
}
\arguments{
\item{z}{MRS data.}
}
\value{
MRS data following Arg operator.
}
\description{
Apply Arg operator to an MRS dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{spin_sys}
\alias{spin_sys}
\title{Create a spin system object for pulse sequence simulation.}
\usage{
spin_sys(spin_params, ft, ref)
}
\arguments{
\item{spin_params}{an object describing the spin system properties.}

\item{ft}{transmitter frequency in Hz.}

\item{ref}{reference value for ppm scale.}
}
\value{
spin system object.
}
\description{
Create a spin system object for pulse sequence simulation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{grid_shift_xy}
\alias{grid_shift_xy}
\title{Grid shift MRSI data in the x/y dimension.}
\usage{
grid_shift_xy(mrs_data, x_shift, y_shift)
}
\arguments{
\item{mrs_data}{MRSI data in the spatial domain.}

\item{x_shift}{shift to apply in the x-direction in units of voxels.}

\item{y_shift}{shift to apply in the y-direction in units of voxels.}
}
\value{
shifted data.
}
\description{
Grid shift MRSI data in the x/y dimension.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varpro_3_para.R
\name{varpro_3_para_opts}
\alias{varpro_3_para_opts}
\title{Return a list of options for VARPRO based fitting with 3 free parameters:
\itemize{
\item zero'th order phase correction
\item global damping
\item global frequency shift.
}}
\usage{
varpro_3_para_opts(
  nstart = 20,
  init_damping = 2,
  maxiters = 200,
  max_shift = 5,
  max_damping = 5,
  anal_jac = FALSE,
  bl_smth_pts = 80
)
}
\arguments{
\item{nstart}{position in the time-domain to start fitting, units of data
points.}

\item{init_damping}{starting value for the global Gaussian line-broadening
term - measured in Hz.}

\item{maxiters}{maximum number of levmar iterations to perform.}

\item{max_shift}{maximum global shift allowed, measured in Hz.}

\item{max_damping}{maximum damping allowed, FWHM measured in Hz.}

\item{anal_jac}{option to use the analytic or numerical Jacobian (logical).}

\item{bl_smth_pts}{number of data points to use in the baseline smoothing
calculation.}
}
\value{
list of options.
}
\description{
Return a list of options for VARPRO based fitting with 3 free parameters:
\itemize{
\item zero'th order phase correction
\item global damping
\item global frequency shift.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spant.R
\name{def_acq_paras}
\alias{def_acq_paras}
\title{Return (and optionally modify using the input arguments) a list of the
default acquisition parameters.}
\usage{
def_acq_paras(
  ft = getOption("spant.def_ft"),
  fs = getOption("spant.def_fs"),
  N = getOption("spant.def_N"),
  ref = getOption("spant.def_ref"),
  nuc = getOption("spant.def_nuc")
)
}
\arguments{
\item{ft}{specify the transmitter frequency in Hz.}

\item{fs}{specify the sampling frequency in Hz.}

\item{N}{specify the number of data points in the spectral dimension.}

\item{ref}{specify the reference value for ppm scale.}

\item{nuc}{specify the resonant nucleus.}
}
\value{
A list containing the following elements:
\itemize{
\item ft transmitter frequency in Hz.
\item fs sampling frequency in Hz.
\item N number of data points in the spectral dimension.
\item ref reference value for ppm scale.
\item nuc resonant nucleus.
}
}
\description{
Return (and optionally modify using the input arguments) a list of the
default acquisition parameters.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{mask_xy_mat}
\alias{mask_xy_mat}
\title{Mask a 2D MRSI dataset in the x-y dimension.}
\usage{
mask_xy_mat(mrs_data, mask, value = NA)
}
\arguments{
\item{mrs_data}{MRS data object.}

\item{mask}{matrix of boolean values specifying the voxels to mask, where a
value of TRUE indicates the voxel should be removed.}

\item{value}{the value to set masked data to (usually NA or 0).}
}
\value{
masked dataset.
}
\description{
Mask a 2D MRSI dataset in the x-y dimension.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mol_parameters.R
\name{get_mol_names}
\alias{get_mol_names}
\title{Return a character array of names that may be used with the
\code{get_mol_paras} function.}
\usage{
get_mol_names()
}
\value{
a character array of names.
}
\description{
Return a character array of names that may be used with the
\code{get_mol_paras} function.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{ecc}
\alias{ecc}
\title{Eddy current correction.}
\usage{
ecc(metab, ref, rev = FALSE)
}
\arguments{
\item{metab}{MRS data to be corrected.}

\item{ref}{reference dataset.}

\item{rev}{reverse the correction.}
}
\value{
corrected data in the time domain.
}
\description{
Apply eddy current correction using the Klose method.
}
\details{
In vivo proton spectroscopy in presence of eddy currents.
Klose U.
Magn Reson Med. 1990 Apr;14(1):26-30.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{get_1h_brain_basis_paras_v1}
\alias{get_1h_brain_basis_paras_v1}
\title{Return a list of \code{mol_parameter} objects suitable for 1H brain MRS
analyses.}
\usage{
get_1h_brain_basis_paras_v1(ft, metab_lw = NULL, lcm_compat = FALSE)
}
\arguments{
\item{ft}{transmitter frequency in Hz.}

\item{metab_lw}{linewidth of metabolite signals (Hz).}

\item{lcm_compat}{when TRUE, lipid, MM and -CrCH molecules will be excluded
from the output.}
}
\value{
list of \code{mol_parameter} objects.
}
\description{
Return a list of \code{mol_parameter} objects suitable for 1H brain MRS
analyses.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{phase}
\alias{phase}
\title{Apply phasing parameters to MRS data.}
\usage{
phase(mrs_data, zero_order, first_order = 0)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{zero_order}{zero'th order phase term in degrees.}

\item{first_order}{first order (frequency dependent) phase term in ms.}
}
\value{
MRS data with applied phase parameters.
}
\description{
Apply phasing parameters to MRS data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{plot_voi_overlay}
\alias{plot_voi_overlay}
\title{Plot a volume as an image overlay.}
\usage{
plot_voi_overlay(mri, voi, export_path = NULL, zlim = NULL, ...)
}
\arguments{
\item{mri}{image data as a nifti object or path to data file.}

\item{voi}{volume data as a nifti object or path to data file.}

\item{export_path}{optional path to save the image in png format.}

\item{zlim}{underlay intensity limits.}

\item{...}{additional arguments to the ortho3 function.}
}
\description{
Plot a volume as an image overlay.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmark.R
\name{spant_abfit_benchmark}
\alias{spant_abfit_benchmark}
\title{Simulate and fit some spectra with ABfit for benchmarking purposes. Basic
timing and performance metrics will be printed.}
\usage{
spant_abfit_benchmark(noise_reps = 10, return_res = FALSE, opts = abfit_opts())
}
\arguments{
\item{noise_reps}{number of spectra to fit with differing noise samples.}

\item{return_res}{return a list of fit_result objects.}

\item{opts}{ABfit options structure.}
}
\description{
Simulate and fit some spectra with ABfit for benchmarking purposes. Basic
timing and performance metrics will be printed.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/amp_scaling.R
\name{scale_amp_molar}
\alias{scale_amp_molar}
\title{Apply water reference scaling to a fitting results object to yield metabolite
quantities in millimolar (mM) units (mol/litre).}
\usage{
scale_amp_molar(fit_result, ref_data, w_att = 0.7, w_conc = 35880, ...)
}
\arguments{
\item{fit_result}{a result object generated from fitting.}

\item{ref_data}{water reference MRS data object.}

\item{w_att}{water attenuation factor (default = 0.7).}

\item{w_conc}{assumed water concentration (default = 35880). Default value
corresponds to typical white matter. Set to 43300 for gray matter, and 55556
for phantom measurements.}

\item{...}{additional arguments to get_td_amp function.}
}
\value{
a \code{fit_result} object with a rescaled results table.
}
\description{
See the LCModel manual section on water-scaling for details on the
assumptions and relevant references.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{apodise_xy}
\alias{apodise_xy}
\title{Apodise MRSI data in the x-y direction with a k-space filter.}
\usage{
apodise_xy(mrs_data, func = "hamming", w = 2.5)
}
\arguments{
\item{mrs_data}{MRSI data.}

\item{func}{must be "hamming" or "gaussian".}

\item{w}{the reciprocal of the standard deviation for the Gaussian function.}
}
\value{
apodised data.
}
\description{
Apodise MRSI data in the x-y direction with a k-space filter.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{comb_metab_ref}
\alias{comb_metab_ref}
\title{Combine a reference and metabolite mrs_data object.}
\usage{
comb_metab_ref(metab, ref)
}
\arguments{
\item{metab}{metabolite mrs_data object.}

\item{ref}{reference mrs_data object.}
}
\value{
combined metabolite and reference mrs_data object.
}
\description{
Combine a reference and metabolite mrs_data object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{sim_zero}
\alias{sim_zero}
\title{Simulate an mrs_data object containing complex zero valued samples.}
\usage{
sim_zero(fs = def_fs(), ft = def_ft(), N = def_N(), ref = def_ref(), dyns = 1)
}
\arguments{
\item{fs}{sampling frequency in Hz.}

\item{ft}{transmitter frequency in Hz.}

\item{N}{number of data points in the spectral dimension.}

\item{ref}{reference value for ppm scale.}

\item{dyns}{number of dynamic scans to generate.}
}
\value{
mrs_data object.
}
\description{
Simulate an mrs_data object containing complex zero valued samples.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{mean_mrs_list}
\alias{mean_mrs_list}
\title{Return the mean of a list of mrs_data objects.}
\usage{
mean_mrs_list(mrs_list)
}
\arguments{
\item{mrs_list}{list of mrs_data objects.}
}
\value{
mean \code{mrs_data} object.
}
\description{
Return the mean of a list of mrs_data objects.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{Mod.mrs_data}
\alias{Mod.mrs_data}
\title{Apply Mod operator to an MRS dataset.}
\usage{
\method{Mod}{mrs_data}(z)
}
\arguments{
\item{z}{MRS data.}
}
\value{
MRS data following Mod operator.
}
\description{
Apply Mod operator to an MRS dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{circ_mask}
\alias{circ_mask}
\title{Create a logical circular mask spanning the full extent of an n x n matrix.}
\usage{
circ_mask(d, n, offset = 1)
}
\arguments{
\item{d}{diameter of the mask.}

\item{n}{number of matrix rows and columns.}

\item{offset}{offset the mask centre in matrix dimension units.}
}
\value{
logical n x n mask matrix.
}
\description{
Create a logical circular mask spanning the full extent of an n x n matrix.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{td_conv_filt}
\alias{td_conv_filt}
\title{Time-domain convolution based filter.}
\usage{
td_conv_filt(mrs_data, K = 25, ext = 1)
}
\arguments{
\item{mrs_data}{MRS data to be filtered.}

\item{K}{window width in data points.}

\item{ext}{point separation for linear extrapolation.}
}
\description{
Time-domain convolution based filter described by:
Marion D, Ikura M, Bax A. Improved solvent suppression in one-dimensional and
twodimensional NMR spectra by convolution of time-domain data. J Magn Reson
1989;84:425-430.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_seg_ind}
\alias{get_seg_ind}
\title{Get the indices of data points lying between two values (end > x > start).}
\usage{
get_seg_ind(scale, start, end)
}
\arguments{
\item{scale}{full list of values.}

\item{start}{smallest value in the subset.}

\item{end}{largest value in the subset.}
}
\value{
set of indices.
}
\description{
Get the indices of data points lying between two values (end > x > start).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_display.R
\name{n2coord}
\alias{n2coord}
\title{Print fit coordinates from a single index.}
\usage{
n2coord(n, fit_res)
}
\arguments{
\item{n}{fit index.}

\item{fit_res}{\code{fit_result} object.}
}
\description{
Print fit coordinates from a single index.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{scale_spec}
\alias{scale_spec}
\title{Scale mrs_data to a spectral region.}
\usage{
scale_spec(
  mrs_data,
  xlim = NULL,
  operator = "sum",
  freq_scale = "ppm",
  mode = "re",
  mean_dyns = TRUE
)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{xlim}{spectral range to be integrated (defaults to full range).}

\item{operator}{can be "sum" (default), "mean", "l2", "max", "min" or
"max-min".}

\item{freq_scale}{units of xlim, can be : "ppm", "Hz" or "points".}

\item{mode}{spectral mode, can be : "re", "im", "mod" or "cplx".}

\item{mean_dyns}{mean the dynamic scans before applying the operator. The
same scaling value will be applied to each individual dynamic.}
}
\value{
normalised data.
}
\description{
Scale mrs_data to a spectral region.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_read_ima.R
\name{read_ima_coil_dir}
\alias{read_ima_coil_dir}
\title{Read a directory containing Siemens MRS IMA files and combine along the coil
dimension. Note that the coil ID is inferred from the sorted file name and
should be checked when consistency is required between two directories.}
\usage{
read_ima_coil_dir(dir, extra = NULL)
}
\arguments{
\item{dir}{data directory path.}

\item{extra}{an optional data frame to provide additional variables for use
in subsequent analysis steps, eg id or grouping variables.}
}
\value{
mrs_data object.
}
\description{
Read a directory containing Siemens MRS IMA files and combine along the coil
dimension. Note that the coil ID is inferred from the sorted file name and
should be checked when consistency is required between two directories.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{get_mrs_affine}
\alias{get_mrs_affine}
\title{Generate an affine for nifti generation.}
\usage{
get_mrs_affine(mrs_data, x_pos = 1, y_pos = 1, z_pos = 1)
}
\arguments{
\item{mrs_data}{input data.}

\item{x_pos}{x_position coordinate.}

\item{y_pos}{y_position coordinate.}

\item{z_pos}{z_position coordinate.}
}
\value{
affine matrix.
}
\description{
Generate an affine for nifti generation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{mean_dyn_blocks}
\alias{mean_dyn_blocks}
\title{Calculate the mean of adjacent dynamic scans.}
\usage{
mean_dyn_blocks(mrs_data, block_size)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}

\item{block_size}{number of adjacent dynamics scans to average over.}
}
\value{
dynamic data averaged in blocks.
}
\description{
Calculate the mean of adjacent dynamic scans.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{lw2beta}
\alias{lw2beta}
\title{Covert a linewidth in Hz to an equivalent beta value in the time-domain ie:
x * exp(-t * t * beta).}
\usage{
lw2beta(lw)
}
\arguments{
\item{lw}{linewidth in Hz.}
}
\value{
beta damping value.
}
\description{
Covert a linewidth in Hz to an equivalent beta value in the time-domain ie:
x * exp(-t * t * beta).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spant.R
\name{check_lcm}
\alias{check_lcm}
\title{Check LCModel can be run}
\usage{
check_lcm()
}
\description{
Check LCModel can be run
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{elliptical_mask}
\alias{elliptical_mask}
\title{Create an elliptical mask stored as a matrix of logical values.}
\usage{
elliptical_mask(xN, yN, x0, y0, xr, yr, angle)
}
\arguments{
\item{xN}{number of pixels in the x dimension.}

\item{yN}{number of pixels in the y dimension.}

\item{x0}{centre of ellipse in the x direction in units of pixels.}

\item{y0}{centre of ellipse in the y direction in units of pixels.}

\item{xr}{radius in the x direction in units of pixels.}

\item{yr}{radius in the y direction in units of pixels.}

\item{angle}{angle of rotation in degrees.}
}
\value{
logical mask matrix with dimensions fov_yN x fov_xN.
}
\description{
Create an elliptical mask stored as a matrix of logical values.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{get_voi_cog}
\alias{get_voi_cog}
\title{Calculate the centre of gravity for an image containing 0 and 1's.}
\usage{
get_voi_cog(voi)
}
\arguments{
\item{voi}{nifti object.}
}
\value{
triplet of x,y,z coordinates.
}
\description{
Calculate the centre of gravity for an image containing 0 and 1's.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_display.R
\name{plot_slice_map}
\alias{plot_slice_map}
\title{Plot a slice from a 7 dimensional array.}
\usage{
plot_slice_map(
  data,
  zlim = NULL,
  mask_map = NULL,
  mask_cutoff = 20,
  interp = 1,
  slice = 1,
  dyn = 1,
  coil = 1,
  ref = 1,
  denom = NULL,
  horizontal = FALSE
)
}
\arguments{
\item{data}{7d array of values to be plotted.}

\item{zlim}{smallest and largest values to be plotted.}

\item{mask_map}{matching map with logical values to indicate if the
corresponding values should be plotted.}

\item{mask_cutoff}{minimum values to plot (as a percentage of the maximum).}

\item{interp}{map interpolation factor.}

\item{slice}{the slice index to plot.}

\item{dyn}{the dynamic index to plot.}

\item{coil}{the coil element number to plot.}

\item{ref}{reference index to plot.}

\item{denom}{map to use as a denominator.}

\item{horizontal}{display the colourbar horizontally (logical).}
}
\description{
Plot a slice from a 7 dimensional array.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{sim_resonances}
\alias{sim_resonances}
\title{Simulate a MRS data object containing a set of simulated resonances.}
\usage{
sim_resonances(
  freq = 0,
  amp = 1,
  lw = 0,
  lg = 0,
  phase = 0,
  freq_ppm = TRUE,
  acq_paras = def_acq_paras(),
  fp_scale = TRUE,
  back_extrap_pts = 0,
  sum_resonances = TRUE
)
}
\arguments{
\item{freq}{resonance frequency.}

\item{amp}{resonance amplitude.}

\item{lw}{line width in Hz.}

\item{lg}{Lorentz-Gauss lineshape parameter (between 0 and 1).}

\item{phase}{phase in degrees.}

\item{freq_ppm}{frequencies are given in ppm units if set to TRUE, otherwise
Hz are assumed.}

\item{acq_paras}{list of acquisition parameters. See
\code{\link{def_acq_paras}}}

\item{fp_scale}{multiply the first data point by 0.5.}

\item{back_extrap_pts}{number of data points to back extrapolate.}

\item{sum_resonances}{sum all resonances (default is TRUE), otherwise return
a dynamic mrs_data object.}
}
\value{
MRS data object.
}
\description{
Simulate a MRS data object containing a set of simulated resonances.
}
\examples{
sim_data <- sim_resonances(freq = 2, lw = 5)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_display.R
\name{gridplot}
\alias{gridplot}
\title{Arrange spectral plots in a grid.}
\usage{
gridplot(x, ...)
}
\arguments{
\item{x}{object for plotting.}

\item{...}{arguments to be passed to methods.}
}
\description{
Arrange spectral plots in a grid.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_read_twix.R
\name{read_siemens_txt_hdr}
\alias{read_siemens_txt_hdr}
\title{Read the text format header found in Siemens IMA and TWIX data files.}
\usage{
read_siemens_txt_hdr(input, version = "vd")
}
\arguments{
\item{input}{file name to read or raw data.}

\item{version}{software version, can be "vb" or "vd".}
}
\value{
a list of parameter values
}
\description{
Read the text format header found in Siemens IMA and TWIX data files.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{fit_diags}
\alias{fit_diags}
\title{Calculate diagnostic information for object of class \code{fit_result}.}
\usage{
fit_diags(x, amps = NULL)
}
\arguments{
\item{x}{\code{fit_result} object.}

\item{amps}{known metabolite amplitudes.}
}
\value{
a dataframe of diagnostic information.
}
\description{
Calculate diagnostic information for object of class \code{fit_result}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{Nz}
\alias{Nz}
\title{Return the total number of z locations in an MRS dataset.}
\usage{
Nz(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\description{
Return the total number of z locations in an MRS dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{fp_phase_correct}
\alias{fp_phase_correct}
\title{Perform a zeroth order phase correction based on the phase of the first data
point in the time-domain.}
\usage{
fp_phase_correct(mrs_data, ret_phase = FALSE)
}
\arguments{
\item{mrs_data}{MRS data to be corrected.}

\item{ret_phase}{return phase values (logical).}
}
\value{
corrected data or a list with corrected data and optional phase
values.
}
\description{
Perform a zeroth order phase correction based on the phase of the first data
point in the time-domain.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spant.R
\name{def_ft}
\alias{def_ft}
\title{Return the default transmitter frequency in Hz.}
\usage{
def_ft()
}
\value{
transmitter frequency in Hz.
}
\description{
Return the default transmitter frequency in Hz.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dicom_reader.R
\name{dicom_reader}
\alias{dicom_reader}
\title{A very simple DICOM reader.}
\usage{
dicom_reader(
  input,
  tags = list(sop_class_uid = "0008,0016"),
  endian = "little",
  debug = FALSE
)
}
\arguments{
\item{input}{either a file path or raw binary object.}

\item{tags}{a named list of tags to be extracted from the file.
eg tags <- list(spec_data = "7FE1,1010", pat_name = "0010,0010")}

\item{endian}{can be "little" or "big".}

\item{debug}{print out some debugging information, can be "little" or "big".}
}
\value{
a list with the same structure as the input, but with tag codes
replaced with the corresponding data in a raw format.
}
\description{
Note this reader is very basic and does not use a DICOM dictionary or try to
convert the data to the correct datatype. For a more robust and sophisticated
reader use the oro.dicom package.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{set_td_pts}
\alias{set_td_pts}
\title{Set the number of time-domain data points, truncating or zero-filling as
appropriate.}
\usage{
set_td_pts(mrs_data, pts)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{pts}{number of data points.}
}
\value{
MRS data with pts data points.
}
\description{
Set the number of time-domain data points, truncating or zero-filling as
appropriate.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{reson_table2mrs_data}
\alias{reson_table2mrs_data}
\title{Generate mrs_data from a table of single Lorentzian resonances.}
\usage{
reson_table2mrs_data(
  reson_table,
  acq_paras = def_acq_paras(),
  back_extrap_pts = 0
)
}
\arguments{
\item{reson_table}{as produced by the hsvd function.}

\item{acq_paras}{list of acquisition parameters. See}

\item{back_extrap_pts}{number of data points to back extrapolate
\code{\link{def_acq_paras}}}
}
\value{
mrs_data object.
}
\description{
Generate mrs_data from a table of single Lorentzian resonances.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pulse_sequences.R
\name{seq_pulse_acquire_31p}
\alias{seq_pulse_acquire_31p}
\title{Simple pulse and acquire sequence with ideal pulses.}
\usage{
seq_pulse_acquire_31p(spin_params, ft, ref)
}
\arguments{
\item{spin_params}{spin system definition.}

\item{ft}{transmitter frequency in Hz.}

\item{ref}{reference value for ppm scale.}
}
\value{
list of resonance amplitudes and frequencies.
}
\description{
Simple pulse and acquire sequence with ideal pulses.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{fp_phase}
\alias{fp_phase}
\title{Return the phase of the first data point in the time-domain.}
\usage{
fp_phase(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\value{
phase values in degrees.
}
\description{
Return the phase of the first data point in the time-domain.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{comb_fit_tables}
\alias{comb_fit_tables}
\title{Combine all fitting data points into a single data frame.}
\usage{
comb_fit_tables(fit_res, inc_basis_sigs = FALSE, inc_indices = TRUE)
}
\arguments{
\item{fit_res}{a single fit_result object.}

\item{inc_basis_sigs}{include the individual fitting basis signals in the
output table, defaults to FALSE.}

\item{inc_indices}{include indices such as X, Y and coil in the output,
defaults to TRUE. These are generally not useful for SVS analysis.}
}
\value{
a data frame containing the fit data points.
}
\description{
Combine all fitting data points into a single data frame.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{Im.mrs_data}
\alias{Im.mrs_data}
\title{Apply Im operator to an MRS dataset.}
\usage{
\method{Im}{mrs_data}(z)
}
\arguments{
\item{z}{MRS data.}
}
\value{
MRS data following Im operator.
}
\description{
Apply Im operator to an MRS dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmark.R
\name{spant_simulation_benchmark}
\alias{spant_simulation_benchmark}
\title{Simulate a typical metabolite basis set for benchmarking. Timing metrics will
be printed on completion.}
\usage{
spant_simulation_benchmark(sim_reps = 10, N = 1024)
}
\arguments{
\item{sim_reps}{number of times to simulate the basis set.}

\item{N}{number of FID data points to simulate.}
}
\description{
Simulate a typical metabolite basis set for benchmarking. Timing metrics will
be printed on completion.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{interleave_dyns}
\alias{interleave_dyns}
\title{Interleave the first and second half of a dynamic series.}
\usage{
interleave_dyns(mrs_data)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}
}
\value{
interleaved data.
}
\description{
Interleave the first and second half of a dynamic series.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{sum_mrs}
\alias{sum_mrs}
\title{Sum two mrs_data objects.}
\usage{
sum_mrs(a, b, force = FALSE)
}
\arguments{
\item{a}{first mrs_data object to be summed.}

\item{b}{second mrs_data object to be summed.}

\item{force}{set to TRUE to force mrs_data objects to be summed, even if they
are in different time/frequency domains.}
}
\value{
a + b
}
\description{
Sum two mrs_data objects.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{bc_als}
\alias{bc_als}
\title{Baseline correction using the ALS method.}
\usage{
bc_als(mrs_data, lambda = 10000, p = 0.001)
}
\arguments{
\item{mrs_data}{mrs_data object.}

\item{lambda}{lambda parameter.}

\item{p}{p parameter.}
}
\value{
baseline corrected data.
}
\description{
Baseline correction using the ALS method.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{get_mrsi2d_seg}
\alias{get_mrsi2d_seg}
\title{Calculate the partial volume estimates for each voxel in a 2D MRSI dataset.}
\usage{
get_mrsi2d_seg(mrs_data, mri_seg, ker)
}
\arguments{
\item{mrs_data}{2D MRSI data with multiple voxels in the x-y dimension.}

\item{mri_seg}{MRI data with values corresponding to the segmentation class.
Must be 1mm isotropic resolution.}

\item{ker}{MRSI PSF kernel in the x-y direction compatible with the mmand
package, eg: mmand::shapeKernel(c(10, 10), type = "box").}
}
\value{
a data frame of partial volume estimates and individual segmentation
maps.
}
\description{
Localisation is assumed to be perfect in the z direction and determined by
the ker input in the x-y direction.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{append_dyns}
\alias{append_dyns}
\title{Append MRS data across the dynamic dimension, assumes they matched across the
other dimensions.}
\usage{
append_dyns(...)
}
\arguments{
\item{...}{MRS data objects as arguments, or a list of MRS data objects.}
}
\value{
a single MRS data object with the input objects concatenated together.
}
\description{
Append MRS data across the dynamic dimension, assumes they matched across the
other dimensions.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tdsr.R
\name{tdsr}
\alias{tdsr}
\title{Time-domain spectral registration.}
\usage{
tdsr(mrs_data, ref = NULL, xlim = c(4, 0.5), max_t = 0.2)
}
\arguments{
\item{mrs_data}{MRS data to be corrected.}

\item{ref}{optional MRS data to use as a reference, the mean of all dynamics
is used if this argument is not supplied.}

\item{xlim}{optional frequency range to perform optimisation, set to NULL
to use the full range.}

\item{max_t}{truncate the FID when longer than max_t to reduce time taken.}
}
\value{
a list containing the corrected data; phase and shift values in units
of degrees and Hz respectively.
}
\description{
An implementation of the method published by Near et al MRM 73:44-50 (2015).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{sum_coils}
\alias{sum_coils}
\title{Calculate the sum across receiver coil elements.}
\usage{
sum_coils(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data split across receiver coil elements.}
}
\value{
sum across coil elements.
}
\description{
Calculate the sum across receiver coil elements.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spant.R
\name{def_N}
\alias{def_N}
\title{Return the default number of data points in the spectral dimension.}
\usage{
def_N()
}
\value{
number of data points in the spectral dimension.
}
\description{
Return the default number of data points in the spectral dimension.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spant.R
\name{get_lcm_cmd}
\alias{get_lcm_cmd}
\title{Print the command to run the LCModel command-line program.}
\usage{
get_lcm_cmd()
}
\description{
Print the command to run the LCModel command-line program.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{recon_twix_2d_mrsi}
\alias{recon_twix_2d_mrsi}
\title{Reconstruct 2D MRSI data from a twix file loaded with read_mrs.}
\usage{
recon_twix_2d_mrsi(twix_mrs)
}
\arguments{
\item{twix_mrs}{raw dynamic data.}
}
\value{
reconstructed data.
}
\description{
Reconstruct 2D MRSI data from a twix file loaded with read_mrs.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{get_1h_brain_basis_paras}
\alias{get_1h_brain_basis_paras}
\title{Return a list of \code{mol_parameter} objects suitable for 1H brain MRS
analyses.}
\usage{
get_1h_brain_basis_paras(ft, metab_lw = NULL, lcm_compat = FALSE)
}
\arguments{
\item{ft}{transmitter frequency in Hz.}

\item{metab_lw}{linewidth of metabolite signals (Hz).}

\item{lcm_compat}{when TRUE, lipid, MM and -CrCH molecules will be excluded
from the output.}
}
\value{
list of \code{mol_parameter} objects.
}
\description{
Return a list of \code{mol_parameter} objects suitable for 1H brain MRS
analyses.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_slice}
\alias{get_slice}
\title{Return a single slice from a larger MRSI dataset.}
\usage{
get_slice(mrs_data, z_pos)
}
\arguments{
\item{mrs_data}{MRSI data.}

\item{z_pos}{the z index to extract.}
}
\value{
MRS data.
}
\description{
Return a single slice from a larger MRSI dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_td_amp}
\alias{get_td_amp}
\title{Return an array of amplitudes derived from fitting the initial points in the
time domain and extrapolating back to t=0.}
\usage{
get_td_amp(mrs_data, nstart = 10, nend = 50, method = "spline")
}
\arguments{
\item{mrs_data}{MRS data.}

\item{nstart}{first data point to fit.}

\item{nend}{last data point to fit.}

\item{method}{method for measuring the amplitude, one of "spline" or "exp".}
}
\value{
array of amplitudes.
}
\description{
Return an array of amplitudes derived from fitting the initial points in the
time domain and extrapolating back to t=0.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{calc_coil_noise_cor}
\alias{calc_coil_noise_cor}
\title{Calculate the noise correlation between coil elements.}
\usage{
calc_coil_noise_cor(noise_data)
}
\arguments{
\item{noise_data}{\code{mrs_data} object with one FID for each coil element.}
}
\value{
correlation matrix.
}
\description{
Calculate the noise correlation between coil elements.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis_set.R
\name{mrs_data2basis}
\alias{mrs_data2basis}
\title{Convert an mrs_data object to basis object - where basis signals are spread
across the dynamic dimension in the MRS data.}
\usage{
mrs_data2basis(mrs_data, names)
}
\arguments{
\item{mrs_data}{mrs_data object with basis signals spread across the dynamic dimension.}

\item{names}{list of names corresponding to basis signals.}
}
\value{
basis set object.
}
\description{
Convert an mrs_data object to basis object - where basis signals are spread
across the dynamic dimension in the MRS data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{Ny}
\alias{Ny}
\title{Return the total number of y locations in an MRS dataset.}
\usage{
Ny(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\description{
Return the total number of y locations in an MRS dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{gausswin_2d}
\alias{gausswin_2d}
\title{Create a two dimensional Gaussian window function stored as a matrix.}
\usage{
gausswin_2d(xN, yN, x0, y0, xw, yw)
}
\arguments{
\item{xN}{number of pixels in the x dimension.}

\item{yN}{number of pixels in the y dimension.}

\item{x0}{centre of window function in the x direction in units of pixels.
Note, only integer values are applied.}

\item{y0}{centre of window function in the y direction in units of pixels.
Note, only integer values are applied.}

\item{xw}{the reciprocal of the standard deviation of the Gaussian window in
x direction.}

\item{yw}{the reciprocal of the standard deviation of the Gaussian window in
y direction.}
}
\value{
matrix with dimensions fov_yN x fov_xN.
}
\description{
Create a two dimensional Gaussian window function stored as a matrix.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mol_parameters.R
\name{get_mol_paras}
\alias{get_mol_paras}
\title{Get a \code{mol_parameters} object for a named molecule.}
\usage{
get_mol_paras(name, ...)
}
\arguments{
\item{name}{the name of the molecule.}

\item{...}{arguments to pass to molecule definition function.}
}
\description{
Get a \code{mol_parameters} object for a named molecule.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{calc_coil_noise_sd}
\alias{calc_coil_noise_sd}
\title{Calculate the noise standard deviation for each coil element.}
\usage{
calc_coil_noise_sd(noise_data)
}
\arguments{
\item{noise_data}{\code{mrs_data} object with one FID for each coil element.}
}
\value{
array of standard deviations.
}
\description{
Calculate the noise standard deviation for each coil element.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{hsvd}
\alias{hsvd}
\title{HSVD of an mrs_data object.}
\usage{
hsvd(mrs_data, comps = 40, irlba = TRUE, max_damp = 10)
}
\arguments{
\item{mrs_data}{mrs_data object to be decomposed.}

\item{comps}{number of Lorentzian components to use for modelling.}

\item{irlba}{option to use irlba SVD (logical).}

\item{max_damp}{maximum allowable damping factor.}
}
\value{
basis matrix and signal table.
}
\description{
HSVD method as described in:
Barkhuijsen H, de Beer R, van Ormondt D. Improved algorithm for noniterative
and timedomain model fitting to exponentially damped magnetic resonance
signals. J Magn Reson 1987;73:553-557.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pulse_sequences.R
\name{seq_press_ideal}
\alias{seq_press_ideal}
\title{PRESS sequence with ideal pulses.}
\usage{
seq_press_ideal(spin_params, ft, ref, TE1 = 0.01, TE2 = 0.02)
}
\arguments{
\item{spin_params}{spin system definition.}

\item{ft}{transmitter frequency in Hz.}

\item{ref}{reference value for ppm scale.}

\item{TE1}{TE1 sequence parameter in seconds (TE=TE1+TE2).}

\item{TE2}{TE2 sequence parameter in seconds.}
}
\value{
list of resonance amplitudes and frequencies.
}
\description{
PRESS sequence with ideal pulses.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pulse_sequences.R
\name{seq_slaser_ideal}
\alias{seq_slaser_ideal}
\title{sLASER sequence with ideal pulses.}
\usage{
seq_slaser_ideal(spin_params, ft, ref, TE1 = 0.008, TE2 = 0.011, TE3 = 0.009)
}
\arguments{
\item{spin_params}{spin system definition.}

\item{ft}{transmitter frequency in Hz.}

\item{ref}{reference value for ppm scale.}

\item{TE1}{first echo time (between exc. and 1st echo) in seconds.}

\item{TE2}{second echo time (between 2nd echo and 4th echo) in seconds.}

\item{TE3}{third echo time (between 4th echo and 5th echo) in seconds.}
}
\value{
list of resonance amplitudes and frequencies.
}
\description{
sLASER sequence with ideal pulses.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{resample_img}
\alias{resample_img}
\title{Resample an image to match a target image space.}
\usage{
resample_img(source, target, interp = 3L)
}
\arguments{
\item{source}{image data as a nifti object.}

\item{target}{image data as a nifti object.}

\item{interp}{interpolation parameter, see nifyreg.linear definition.}
}
\value{
resampled image data as a nifti object.
}
\description{
Resample an image to match a target image space.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{get_mrsi_voxel_xy_psf}
\alias{get_mrsi_voxel_xy_psf}
\title{Generate a MRSI voxel PSF from an \code{mrs_data} object.}
\usage{
get_mrsi_voxel_xy_psf(mrs_data, target_mri, x_pos, y_pos, z_pos)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{target_mri}{optional image data to match the intended volume space.}

\item{x_pos}{x voxel coordinate.}

\item{y_pos}{y voxel coordinate.}

\item{z_pos}{z voxel coordinate.}
}
\value{
volume data as a nifti object.
}
\description{
Generate a MRSI voxel PSF from an \code{mrs_data} object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{fp_scale}
\alias{fp_scale}
\title{Scale the first time-domain data point in an mrs_data object.}
\usage{
fp_scale(mrs_data, scale = 0.5)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{scale}{scaling value, defaults to 0.5.}
}
\value{
scaled mrs_data object.
}
\description{
Scale the first time-domain data point in an mrs_data object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{mean.mrs_data}
\alias{mean.mrs_data}
\title{Calculate the mean spectrum from an mrs_data object.}
\usage{
\method{mean}{mrs_data}(x, ...)
}
\arguments{
\item{x}{object of class mrs_data.}

\item{...}{other arguments to pass to the colMeans function.}
}
\value{
mean mrs_data object.
}
\description{
Calculate the mean spectrum from an mrs_data object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{Conj.mrs_data}
\alias{Conj.mrs_data}
\title{Apply Conj operator to an MRS dataset.}
\usage{
\method{Conj}{mrs_data}(z)
}
\arguments{
\item{z}{MRS data.}
}
\value{
MRS data following Conj operator.
}
\description{
Apply Conj operator to an MRS dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_even_dyns}
\alias{get_even_dyns}
\title{Return even numbered dynamic scans starting from 1 (2,4,6...).}
\usage{
get_even_dyns(mrs_data)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}
}
\value{
dynamic MRS data containing even numbered scans.
}
\description{
Return even numbered dynamic scans starting from 1 (2,4,6...).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pulse_sequences.R
\name{seq_spin_echo_ideal_31p}
\alias{seq_spin_echo_ideal_31p}
\title{Spin echo sequence with ideal pulses.}
\usage{
seq_spin_echo_ideal_31p(spin_params, ft, ref, TE = 0.03)
}
\arguments{
\item{spin_params}{spin system definition.}

\item{ft}{transmitter frequency in Hz.}

\item{ref}{reference value for ppm scale.}

\item{TE}{echo time in seconds.}
}
\value{
list of resonance amplitudes and frequencies.
}
\description{
Spin echo sequence with ideal pulses.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{l2_reg}
\alias{l2_reg}
\title{Perform l2 regularisation artefact suppression.}
\usage{
l2_reg(
  mrs_data,
  thresh = 0.05,
  b = 1e-11,
  A = NA,
  xlim = NA,
  thresh_xlim = NULL
)
}
\arguments{
\item{mrs_data}{input data for artefact suppression.}

\item{thresh}{threshold parameter to extract lipid signals from mrs_data
based on the spectral integration of the thresh_xlim region in magnitude
mode.}

\item{b}{regularisation parameter.}

\item{A}{set of spectra containing the artefact basis signals. The thresh
parameter is ignored when A is specified.}

\item{xlim}{spectral limits in ppm to restrict the reconstruction range.
Defaults to the full spectral width.}

\item{thresh_xlim}{spectral limits in ppm to integrate for the threshold map.}
}
\value{
l2 reconstructed mrs_data object.
}
\description{
Perform l2 regularisation artefact suppression using the method proposed by
Bilgic et al. JMRI 40(1):181-91 2014.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_display.R
\name{print.fit_result}
\alias{print.fit_result}
\title{Print a summary of an object of class \code{fit_result}.}
\usage{
\method{print}{fit_result}(x, ...)
}
\arguments{
\item{x}{\code{fit_result} object.}

\item{...}{further arguments.}
}
\description{
Print a summary of an object of class \code{fit_result}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{gen_F}
\alias{gen_F}
\title{Generate the F product operator.}
\usage{
gen_F(sys, op, detect = NULL)
}
\arguments{
\item{sys}{spin system object.}

\item{op}{operator, one of "x", "y", "z", "p", "m".}

\item{detect}{detection nuclei.}
}
\value{
F product operator matrix.
}
\description{
Generate the F product operator.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{is.def}
\alias{is.def}
\title{Check if an object is defined, which is the same as being not NULL.}
\usage{
is.def(x)
}
\arguments{
\item{x}{object to test for being NULL.}
}
\value{
logical value.
}
\description{
Check if an object is defined, which is the same as being not NULL.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{get_voi_seg_psf}
\alias{get_voi_seg_psf}
\title{Return the white matter, gray matter and CSF composition of a volume.}
\usage{
get_voi_seg_psf(psf, mri_seg)
}
\arguments{
\item{psf}{volume data as a nifti object.}

\item{mri_seg}{segmented brain volume as a nifti object.}
}
\value{
a vector of partial volumes expressed as percentages.
}
\description{
Return the white matter, gray matter and CSF composition of a volume.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{downsample_mrs_fd}
\alias{downsample_mrs_fd}
\title{Downsample an MRS signal by a factor of 2 using an FFT "brick-wall" filter.}
\usage{
downsample_mrs_fd(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data object.}
}
\value{
downsampled data.
}
\description{
Downsample an MRS signal by a factor of 2 using an FFT "brick-wall" filter.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_odd_dyns}
\alias{get_odd_dyns}
\title{Return odd numbered dynamic scans starting from 1 (1,3,5...).}
\usage{
get_odd_dyns(mrs_data)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}
}
\value{
dynamic MRS data containing odd numbered scans.
}
\description{
Return odd numbered dynamic scans starting from 1 (1,3,5...).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{crossprod_3d}
\alias{crossprod_3d}
\title{Compute the vector cross product between vectors x and y. Adapted from
http://stackoverflow.com/questions/15162741/what-is-rs-crossproduct-function}
\usage{
crossprod_3d(x, y)
}
\arguments{
\item{x}{vector of length 3.}

\item{y}{vector of length 3.}
}
\value{
vector cross product of x and y.
}
\description{
Compute the vector cross product between vectors x and y. Adapted from
http://stackoverflow.com/questions/15162741/what-is-rs-crossproduct-function
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_sh_dyns}
\alias{get_sh_dyns}
\title{Return the second half of a dynamic series.}
\usage{
get_sh_dyns(mrs_data)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}
}
\value{
second half of the dynamic series.
}
\description{
Return the second half of a dynamic series.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{set_mask_xy_mat}
\alias{set_mask_xy_mat}
\title{Set the masked voxels in a 2D MRSI dataset to given spectrum.}
\usage{
set_mask_xy_mat(mrs_data, mask, mask_mrs_data)
}
\arguments{
\item{mrs_data}{MRSI data object.}

\item{mask}{matrix of boolean values specifying the voxels to set, where a
value of TRUE indicates the voxel should be set to mask_mrs_data.}

\item{mask_mrs_data}{the spectral data to be assigned to the masked voxels.}
}
\value{
updated dataset.
}
\description{
Set the masked voxels in a 2D MRSI dataset to given spectrum.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interactive_plotting.R
\name{ortho3}
\alias{ortho3}
\title{Display an orthographic projection plot of a nifti object.}
\usage{
ortho3(
  underlay,
  overlay = NULL,
  xyz = NULL,
  zlim = NULL,
  zlim_ol = NULL,
  alpha = 0.7,
  col_ol = viridisLite::viridis(64),
  orient_lab = TRUE,
  rescale = 1,
  crosshairs = TRUE,
  ch_lwd = 1,
  colourbar = TRUE,
  bg = "black",
  mar = c(0, 0, 0, 0),
  smallplot = c(0.63, 0.65, 0.07, 0.42)
)
}
\arguments{
\item{underlay}{underlay image to be shown in grayscale.}

\item{overlay}{optional overlay image.}

\item{xyz}{x, y, z slice coordinates to display.}

\item{zlim}{underlay intensity limits.}

\item{zlim_ol}{overlay intensity limits.}

\item{alpha}{transparency of overlay.}

\item{col_ol}{colour palette of overlay.}

\item{orient_lab}{display orientation labels (default TRUE).}

\item{rescale}{rescale factor for the underlay and overlay images.}

\item{crosshairs}{display the crosshairs (default TRUE).}

\item{ch_lwd}{crosshair linewidth.}

\item{colourbar}{display a colourbar for the overlay (default TRUE).}

\item{bg}{plot background colour.}

\item{mar}{plot margins.}

\item{smallplot}{smallplot option for positioning the colourbar.}
}
\description{
Display an orthographic projection plot of a nifti object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_display.R
\name{plot.mrs_data}
\alias{plot.mrs_data}
\title{Plotting method for objects of class mrs_data.}
\usage{
\method{plot}{mrs_data}(
  x,
  dyn = 1,
  x_pos = 1,
  y_pos = 1,
  z_pos = 1,
  coil = 1,
  fd = TRUE,
  x_units = NULL,
  xlim = NULL,
  y_scale = FALSE,
  x_ax = TRUE,
  mode = "re",
  lwd = NULL,
  bty = NULL,
  label = "",
  restore_def_par = TRUE,
  mar = NULL,
  xaxis_lab = NULL,
  xat = NULL,
  xlabs = TRUE,
  yat = NULL,
  ylabs = TRUE,
  show_grid = TRUE,
  grid_nx = NULL,
  grid_ny = NA,
  col = NULL,
  alpha = NULL,
  ...
)
}
\arguments{
\item{x}{object of class mrs_data.}

\item{dyn}{the dynamic index to plot.}

\item{x_pos}{the x index to plot.}

\item{y_pos}{the y index to plot.}

\item{z_pos}{the z index to plot.}

\item{coil}{the coil element number to plot.}

\item{fd}{display data in the frequency-domain (default), or time-domain
(logical).}

\item{x_units}{the units to use for the x-axis, can be one of: "ppm", "hz",
"points" or "seconds".}

\item{xlim}{the range of values to display on the x-axis, eg xlim = c(4,1).}

\item{y_scale}{option to display the y-axis values (logical).}

\item{x_ax}{option to display the x-axis values (logical).}

\item{mode}{representation of the complex numbers to be plotted, can be one
of: "re", "im", "mod" or "arg".}

\item{lwd}{plot linewidth.}

\item{bty}{option to draw a box around the plot. See ?par.}

\item{label}{character string to add to the top left of the plot window.}

\item{restore_def_par}{restore default plotting par values after the plot has
been made.}

\item{mar}{option to adjust the plot margins. See ?par.}

\item{xaxis_lab}{x-axis label.}

\item{xat}{x-axis tick label values.}

\item{xlabs}{x-axis tick labels.}

\item{yat}{y-axis tick label values.}

\item{ylabs}{y-axis tick labels.}

\item{show_grid}{plot gridlines behind the data (logical). Defaults to TRUE.}

\item{grid_nx}{number of cells of the grid in x and y direction. When NULL
the grid aligns with the tick marks on the corresponding default axis (i.e.,
tickmarks as computed by axTicks). When NA, no grid lines are drawn in the
corresponding direction.}

\item{grid_ny}{as above.}

\item{col}{set the line colour, eg col = rgb(0.5, 0.5, 0.5).}

\item{alpha}{set the line transparency, eg alpha = 0.5 is 50\% transparency.
Overrides any transparency levels set by col.}

\item{...}{other arguments to pass to the plot method.}
}
\description{
Plotting method for objects of class mrs_data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_display.R
\name{stackplot}
\alias{stackplot}
\title{Produce a plot with multiple traces.}
\usage{
stackplot(x, ...)
}
\arguments{
\item{x}{object for plotting.}

\item{...}{arguments to be passed to methods.}
}
\description{
Produce a plot with multiple traces.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pulse_shapes.R
\name{get_guassian_pulse}
\alias{get_guassian_pulse}
\title{Generate a gaussian pulse shape.}
\usage{
get_guassian_pulse(angle, n, trunc = 1)
}
\arguments{
\item{angle}{pulse angle in degrees.}

\item{n}{number of points to generate.}

\item{trunc}{percentage truncation factor.}
}
\description{
Generate a gaussian pulse shape.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interactive_plotting.R
\name{plot_slice_map_inter}
\alias{plot_slice_map_inter}
\title{Plot an interactive slice map from a data array where voxels can be selected
to display a corresponding spectrum.}
\usage{
plot_slice_map_inter(
  mrs_data,
  map = NULL,
  xlim = NULL,
  slice = 1,
  zlim = NULL,
  mask_map = NULL,
  denom = NULL,
  mask_cutoff = 20,
  interp = 1,
  mode = "re",
  y_scale = FALSE,
  ylim = NULL,
  coil = 1,
  fd = TRUE
)
}
\arguments{
\item{mrs_data}{spectral data.}

\item{map}{array of values to be plotted, defaults to the integration of the
modulus of the full spectral width.}

\item{xlim}{spectral region to plot.}

\item{slice}{the slice index to plot.}

\item{zlim}{smallest and largest values to be plotted.}

\item{mask_map}{matching map with logical values to indicate if the
corresponding values should be plotted.}

\item{denom}{map to use as a denominator.}

\item{mask_cutoff}{minimum values to plot (as a percentage of the maximum).}

\item{interp}{map interpolation factor.}

\item{mode}{representation of the complex spectrum to be plotted, can be one
of: "re", "im", "mod" or "arg".}

\item{y_scale}{option to display the y-axis values (logical).}

\item{ylim}{intensity range to plot.}

\item{coil}{coil element to plot.}

\item{fd}{display data in the frequency-domain (default), or time-domain
(logical).}
}
\description{
Plot an interactive slice map from a data array where voxels can be selected
to display a corresponding spectrum.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_display.R
\name{image.mrs_data}
\alias{image.mrs_data}
\title{Image plot method for objects of class mrs_data.}
\usage{
\method{image}{mrs_data}(
  x,
  xlim = NULL,
  mode = "re",
  col = NULL,
  plot_dim = NULL,
  x_pos = NULL,
  y_pos = NULL,
  z_pos = NULL,
  dyn = 1,
  coil = 1,
  restore_def_par = TRUE,
  y_ticks = NULL,
  vline = NULL,
  hline = NULL,
  ...
)
}
\arguments{
\item{x}{object of class mrs_data.}

\item{xlim}{the range of values to display on the x-axis, eg xlim = c(4,1).}

\item{mode}{representation of the complex numbers to be plotted, can be one
of: "re", "im", "mod" or "arg".}

\item{col}{Colour map to use, defaults to viridis.}

\item{plot_dim}{the dimension to display on the y-axis, can be one of: "dyn",
"x", "y", "z", "coil" or NULL. If NULL (the default) all spectra will be
collapsed into the dynamic dimension and displayed.}

\item{x_pos}{the x index to plot.}

\item{y_pos}{the y index to plot.}

\item{z_pos}{the z index to plot.}

\item{dyn}{the dynamic index to plot.}

\item{coil}{the coil element number to plot.}

\item{restore_def_par}{restore default plotting par values after the plot has
been made.}

\item{y_ticks}{a vector of indices specifying where to place tick marks.}

\item{vline}{draw a vertical line at the value of vline.}

\item{hline}{draw a horizontal line at the value of hline.}

\item{...}{other arguments to pass to the plot method.}
}
\description{
Image plot method for objects of class mrs_data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{decimate_mrs_td}
\alias{decimate_mrs_td}
\title{Decimate an MRS signal by filtering in the time domain before downsampling.}
\usage{
decimate_mrs_td(mrs_data, q = 2, n = 4, ftype = "iir")
}
\arguments{
\item{mrs_data}{MRS data object.}

\item{q}{integer factor to downsample by (default = 2).}

\item{n}{filter order used in the downsampling.}

\item{ftype}{filter type, "iir" or "fir".}
}
\value{
decimated data.
}
\description{
Decimate an MRS signal by filtering in the time domain before downsampling.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{plot_voi_overlay_seg}
\alias{plot_voi_overlay_seg}
\title{Plot a volume as an overlay on a segmented brain volume.}
\usage{
plot_voi_overlay_seg(mri_seg, voi, export_path = NULL, ...)
}
\arguments{
\item{mri_seg}{segmented brain volume as a nifti object.}

\item{voi}{volume data as a nifti object.}

\item{export_path}{optional path to save the image in png format.}

\item{...}{additional arguments to the ortho3 function.}
}
\description{
Plot a volume as an overlay on a segmented brain volume.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis_set.R
\name{read_basis}
\alias{read_basis}
\title{Read a basis file in LCModel .basis format.}
\usage{
read_basis(basis_file, ref = def_ref(), sort_basis = TRUE)
}
\arguments{
\item{basis_file}{path to basis file.}

\item{ref}{assumed ppm reference value.}

\item{sort_basis}{sort the basis set based on signal names.}
}
\value{
basis object.
}
\description{
Read a basis file in LCModel .basis format.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis_set.R
\name{write_basis}
\alias{write_basis}
\title{Write a basis object to an LCModel .basis formatted file.}
\usage{
write_basis(basis, basis_file, fwhmba = 0.1)
}
\arguments{
\item{basis}{basis object to be exported.}

\item{basis_file}{path to basis file to be generated.}

\item{fwhmba}{parameter used by LCModel.}
}
\description{
Write a basis object to an LCModel .basis formatted file.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\alias{\%$\%}
\alias{readNifti}
\alias{writeNifti}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{magrittr}{\code{\link[magrittr:exposition]{\%$\%}}, \code{\link[magrittr:pipe]{\%>\%}}}

  \item{RNifti}{\code{\link[RNifti]{readNifti}}, \code{\link[RNifti]{writeNifti}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_2d_psf}
\alias{get_2d_psf}
\title{Get the point spread function (PSF) for a 2D phase encoded MRSI scan.}
\usage{
get_2d_psf(FOV = 160, mat_size = 16, sampling = "circ", hamming = FALSE)
}
\arguments{
\item{FOV}{field of view in mm.}

\item{mat_size}{acquisition matrix size (not interpolated).}

\item{sampling}{can be either "circ" for circular or "rect" for rectangular.}

\item{hamming}{should Hamming k-space weighting be applied (default FALSE).}
}
\value{
A matrix of the PSF with 1mm resolution.
}
\description{
Get the point spread function (PSF) for a 2D phase encoded MRSI scan.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis_set.R
\name{append_basis}
\alias{append_basis}
\title{Combine a pair of basis set objects.}
\usage{
append_basis(basis_a, basis_b)
}
\arguments{
\item{basis_a}{first basis.}

\item{basis_b}{second basis.}
}
\value{
combined basis set object.
}
\description{
Combine a pair of basis set objects.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_voxel}
\alias{get_voxel}
\title{Return a single voxel from a larger mrs dataset.}
\usage{
get_voxel(mrs_data, x_pos = 1, y_pos = 1, z_pos = 1, dyn = 1, coil = 1)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{x_pos}{the x index to plot.}

\item{y_pos}{the y index to plot.}

\item{z_pos}{the z index to plot.}

\item{dyn}{the dynamic index to plot.}

\item{coil}{the coil element number to plot.}
}
\value{
MRS data.
}
\description{
Return a single voxel from a larger mrs dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{inv_odd_dyns}
\alias{inv_odd_dyns}
\title{Invert odd numbered dynamic scans starting from 1 (1,3,5...).}
\usage{
inv_odd_dyns(mrs_data)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}
}
\value{
dynamic MRS data with inverted odd numbered scans.
}
\description{
Invert odd numbered dynamic scans starting from 1 (1,3,5...).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{lw2alpha}
\alias{lw2alpha}
\title{Covert a linewidth in Hz to an equivalent alpha value in the time-domain ie:
x * exp(-t * alpha).}
\usage{
lw2alpha(lw)
}
\arguments{
\item{lw}{linewidth in Hz.}
}
\value{
beta damping value.
}
\description{
Covert a linewidth in Hz to an equivalent alpha value in the time-domain ie:
x * exp(-t * alpha).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{ft_dyn}
\alias{ft_dyn}
\title{Apply the Fourier transform over the dynamic dimension.}
\usage{
ft_dyn(mrs_data, ft_shift = FALSE)
}
\arguments{
\item{mrs_data}{MRS data where the dynamic dimension is in the time-domain.}

\item{ft_shift}{apply FT shift to the output, default is FALSE.}
}
\value{
transformed MRS data.
}
\description{
Apply the Fourier transform over the dynamic dimension.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{comb_fit_list_fit_tables}
\alias{comb_fit_list_fit_tables}
\title{Combine all fitting data points from a list of fits into a single data frame.}
\usage{
comb_fit_list_fit_tables(
  fit_list,
  add_extra = TRUE,
  harmonise_ppm = TRUE,
  inc_basis_sigs = FALSE,
  inc_indices = TRUE,
  add_res_id = TRUE
)
}
\arguments{
\item{fit_list}{list of fit_result objects.}

\item{add_extra}{add variables in the extra data frame to the output (TRUE).}

\item{harmonise_ppm}{ensure the ppm scale for each fit is identical to the
first.}

\item{inc_basis_sigs}{include the individual fitting basis signals in the
output table, defaults to FALSE.}

\item{inc_indices}{include indices such as X, Y and coil in the output,
defaults to TRUE. These are generally not useful for SVS analysis.}

\item{add_res_id}{add a res_id column to the output to distinguish between
datasets.}
}
\value{
a data frame containing the fit data points.
}
\description{
Combine all fitting data points from a list of fits into a single data frame.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{sim_basis}
\alias{sim_basis}
\title{Simulate a basis set object.}
\usage{
sim_basis(
  mol_list,
  pul_seq = seq_pulse_acquire,
  acq_paras = def_acq_paras(),
  xlim = NULL,
  ...
)
}
\arguments{
\item{mol_list}{list of \code{mol_parameter} objects.}

\item{pul_seq}{pulse sequence function to use.}

\item{acq_paras}{list of acquisition parameters or an mrs_data object. See
\code{\link{def_acq_paras}}}

\item{xlim}{ppm range limiting signals to be simulated.}

\item{...}{extra parameters to pass to the pulse sequence function.}
}
\value{
basis object.
}
\description{
Simulate a basis set object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pulse_sequences.R
\name{seq_spin_echo_ideal}
\alias{seq_spin_echo_ideal}
\title{Spin echo sequence with ideal pulses.}
\usage{
seq_spin_echo_ideal(spin_params, ft, ref, TE = 0.03)
}
\arguments{
\item{spin_params}{spin system definition.}

\item{ft}{transmitter frequency in Hz.}

\item{ref}{reference value for ppm scale.}

\item{TE}{echo time in seconds.}
}
\value{
list of resonance amplitudes and frequencies.
}
\description{
Spin echo sequence with ideal pulses.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{rectangular_mask}
\alias{rectangular_mask}
\title{Create a rectangular mask stored as a matrix of logical values.}
\usage{
rectangular_mask(xN, yN, x0, y0, xw, yw, angle)
}
\arguments{
\item{xN}{number of pixels in the x dimension.}

\item{yN}{number of pixels in the y dimension.}

\item{x0}{centre of rectangle in the x direction in units of pixels.}

\item{y0}{centre of rectangle in the y direction in units of pixels.}

\item{xw}{width in the x direction in units of pixels.}

\item{yw}{width in the y direction in units of pixels.}

\item{angle}{angle of rotation in degrees.}
}
\value{
logical mask matrix with dimensions fov_yN x fov_xN.
}
\description{
Create a rectangular mask stored as a matrix of logical values.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spant.R
\docType{package}
\name{spant-package}
\alias{spant}
\alias{spant-package}
\title{spant: spectroscopy analysis tools.}
\description{
spant provides a set of tools for reading, visualising and processing
Magnetic Resonance Spectroscopy (MRS) data.
}
\details{
To get started with spant, take a look at the introduction vignette:

\code{vignette("spant-intro", package="spant")}

Full list of vignettes:

\code{browseVignettes(package = "spant")}

Full list of functions:

\code{help(package = spant, help_type = "html")}

An online version of the documentation is available from:

\url{https://martin3141.github.io/spant/}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://martin3141.github.io/spant/}
  \item \url{https://github.com/martin3141/spant/}
  \item Report bugs at \url{https://github.com/martin3141/spant/issues/}
}

}
\author{
\strong{Maintainer}: Martin Wilson \email{martin@pipegrep.co.uk} (\href{https://orcid.org/0000-0002-2089-3956}{ORCID})

Other contributors:
\itemize{
  \item Yong Wang [contributor]
  \item John Muschelli [contributor]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{nifti_flip_lr}
\alias{nifti_flip_lr}
\title{Flip the x data dimension order of a nifti image. This corresponds to
flipping MRI data in the left-right direction, assuming the data in save in
neurological format (can check with fslorient program).}
\usage{
nifti_flip_lr(x)
}
\arguments{
\item{x}{nifti object to be processed.}
}
\value{
nifti object with reversed x data direction.
}
\description{
Flip the x data dimension order of a nifti image. This corresponds to
flipping MRI data in the left-right direction, assuming the data in save in
neurological format (can check with fslorient program).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{apply_mrs}
\alias{apply_mrs}
\title{Apply a function across given dimensions of a MRS data object.}
\usage{
apply_mrs(mrs_data, dims, fun, ..., data_only = FALSE)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{dims}{dimensions to apply the function.}

\item{fun}{name of the function.}

\item{...}{arguments to the function.}

\item{data_only}{return an array rather than an MRS data object.}
}
\description{
Apply a function across given dimensions of a MRS data object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{mean_dyn_pairs}
\alias{mean_dyn_pairs}
\title{Calculate the pairwise means across a dynamic data set.}
\usage{
mean_dyn_pairs(mrs_data)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}
}
\value{
mean dynamic data of adjacent dynamic pairs.
}
\description{
Calculate the pairwise means across a dynamic data set.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interactive_plotting.R
\name{plot_slice_fit_inter}
\alias{plot_slice_fit_inter}
\title{Plot a 2D slice from an MRSI fit result object.}
\usage{
plot_slice_fit_inter(
  fit_res,
  map = NULL,
  map_denom = NULL,
  slice = 1,
  zlim = NULL,
  interp = 1,
  xlim = NULL
)
}
\arguments{
\item{fit_res}{\code{fit_result} object.}

\item{map}{fit result values to display as a colour map. Can be specified as
a character string or array of numeric values. Defaults to "tNAA".}

\item{map_denom}{fit result values to divide the map argument by. Can be
specified as a character string (eg "tCr") or array of numeric values.}

\item{slice}{slice to plot in the z direction.}

\item{zlim}{range of values to plot.}

\item{interp}{interpolation factor.}

\item{xlim}{spectral plot limits for the x axis.}
}
\description{
Plot a 2D slice from an MRSI fit result object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pulse_sequences.R
\name{seq_steam_ideal}
\alias{seq_steam_ideal}
\title{STEAM sequence with ideal pulses.}
\usage{
seq_steam_ideal(spin_params, ft, ref, TE = 0.03, TM = 0.02)
}
\arguments{
\item{spin_params}{spin system definition.}

\item{ft}{transmitter frequency in Hz.}

\item{ref}{reference value for ppm scale.}

\item{TE}{sequence parameter in seconds.}

\item{TM}{sequence parameter in seconds.}
}
\value{
list of resonance amplitudes and frequencies.
}
\description{
STEAM sequence with ideal pulses.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{sd}
\alias{sd}
\title{Calculate the standard deviation spectrum from an mrs_data object.}
\usage{
sd(x, na.rm)
}
\arguments{
\item{x}{object of class mrs_data.}

\item{na.rm}{remove NA values.}
}
\value{
sd mrs_data object.
}
\description{
Calculate the standard deviation spectrum from an mrs_data object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ft_shift_mat}
\alias{ft_shift_mat}
\title{Perform a fft and fftshift on a matrix with each column replaced by its
shifted fft.}
\usage{
ft_shift_mat(mat_in)
}
\arguments{
\item{mat_in}{matrix input.}
}
\value{
output matrix.
}
\description{
Perform a fft and fftshift on a matrix with each column replaced by its
shifted fft.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{peak_info}
\alias{peak_info}
\title{Search for the highest peak in a spectral region and return the frequency,
height and FWHM.}
\usage{
peak_info(
  mrs_data,
  xlim = c(4, 0.5),
  interp_f = 4,
  scale = "ppm",
  mode = "real"
)
}
\arguments{
\item{mrs_data}{an object of class \code{mrs_data}.}

\item{xlim}{frequency range (default units of PPM) to search for the highest
peak.}

\item{interp_f}{interpolation factor, defaults to 4x.}

\item{scale}{the units to use for the frequency scale, can be one of: "ppm",
"hz" or "points".}

\item{mode}{spectral mode, can be : "real", "imag" or "mod".}
}
\value{
list of arrays containing the highest peak frequency, height and FWHM
in units of PPM and Hz.
}
\description{
Search for the highest peak in a spectral region and return the frequency,
height and FWHM.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{sim_basis_1h_brain_press}
\alias{sim_basis_1h_brain_press}
\title{Simulate a basis-set suitable for 1H brain MRS analysis acquired with a PRESS
sequence. Note, ideal pulses are assumed.}
\usage{
sim_basis_1h_brain_press(
  acq_paras = def_acq_paras(),
  xlim = c(0.5, 4.2),
  lcm_compat = FALSE,
  TE1 = 0.01,
  TE2 = 0.02
)
}
\arguments{
\item{acq_paras}{list of acquisition parameters or an mrs_data object. See
\code{\link{def_acq_paras}}}

\item{xlim}{range of frequencies to simulate in ppm.}

\item{lcm_compat}{exclude lipid and MM signals for use with default LCModel
options.}

\item{TE1}{TE1 of PRESS sequence (TE = TE1 + TE2).}

\item{TE2}{TE2 of PRESS sequence.}
}
\value{
basis object.
}
\description{
Simulate a basis-set suitable for 1H brain MRS analysis acquired with a PRESS
sequence. Note, ideal pulses are assumed.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{img2kspace_xy}
\alias{img2kspace_xy}
\title{Transform 2D MRSI data to k-space in the x-y direction.}
\usage{
img2kspace_xy(mrs_data)
}
\arguments{
\item{mrs_data}{2D MRSI data.}
}
\value{
k-space data.
}
\description{
Transform 2D MRSI data to k-space in the x-y direction.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pulse_sequences.R
\name{seq_mega_press_ideal}
\alias{seq_mega_press_ideal}
\title{MEGA-PRESS sequence with ideal localisation pulses and Gaussian shaped
editing pulse.}
\usage{
seq_mega_press_ideal(
  spin_params,
  ft,
  ref,
  ed_freq = 1.89,
  TE1 = 0.015,
  TE2 = 0.053,
  BW = 110,
  steps = 50
)
}
\arguments{
\item{spin_params}{spin system definition.}

\item{ft}{transmitter frequency in Hz.}

\item{ref}{reference value for ppm scale.}

\item{ed_freq}{editing pulse frequency in ppm.}

\item{TE1}{TE1 sequence parameter in seconds (TE=TE1+TE2).}

\item{TE2}{TE2 sequence parameter in seconds.}

\item{BW}{editing pulse bandwidth in Hz.}

\item{steps}{number of hard pulses used to approximate the editing pulse.}
}
\value{
list of resonance amplitudes and frequencies.
}
\description{
MEGA-PRESS sequence with ideal localisation pulses and Gaussian shaped
editing pulse.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{mrs_data2vec}
\alias{mrs_data2vec}
\title{Convert mrs_data object to a vector.}
\usage{
mrs_data2vec(mrs_data, dyn = 1, x_pos = 1, y_pos = 1, z_pos = 1, coil = 1)
}
\arguments{
\item{mrs_data}{MRS data object.}

\item{dyn}{dynamic index.}

\item{x_pos}{x index.}

\item{y_pos}{y index.}

\item{z_pos}{z index.}

\item{coil}{coil element index.}
}
\value{
MRS data vector.
}
\description{
Convert mrs_data object to a vector.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{sd.mrs_data}
\alias{sd.mrs_data}
\title{Calculate the standard deviation spectrum from an mrs_data object.}
\usage{
\method{sd}{mrs_data}(x, na.rm = FALSE)
}
\arguments{
\item{x}{object of class mrs_data.}

\item{na.rm}{remove NA values.}
}
\value{
sd mrs_data object.
}
\description{
Calculate the standard deviation spectrum from an mrs_data object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{get_mrsi_voi}
\alias{get_mrsi_voi}
\title{Generate a MRSI VOI from an \code{mrs_data} object.}
\usage{
get_mrsi_voi(mrs_data, target_mri = NULL, map = NULL, ker = mmand::boxKernel())
}
\arguments{
\item{mrs_data}{MRS data.}

\item{target_mri}{optional image data to match the intended volume space.}

\item{map}{optional voi intensity map.}

\item{ker}{kernel to rescale the map data to the target_mri.}
}
\value{
volume data as a nifti object.
}
\description{
Generate a MRSI VOI from an \code{mrs_data} object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{spec_op}
\alias{spec_op}
\title{Perform a mathematical operation on a spectral region.}
\usage{
spec_op(
  mrs_data,
  xlim = NULL,
  operator = "sum",
  freq_scale = "ppm",
  mode = "re"
)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{xlim}{spectral range to be integrated (defaults to full range).}

\item{operator}{can be "sum" (default), "mean", "l2", "max", "min" or
"max-min".}

\item{freq_scale}{units of xlim, can be : "ppm", "hz" or "points".}

\item{mode}{spectral mode, can be : "re", "im", "mod" or "cplx".}
}
\value{
an array of integral values.
}
\description{
Perform a mathematical operation on a spectral region.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{get_mrsi_voxel}
\alias{get_mrsi_voxel}
\title{Generate a MRSI voxel from an \code{mrs_data} object.}
\usage{
get_mrsi_voxel(mrs_data, target_mri, x_pos, y_pos, z_pos)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{target_mri}{optional image data to match the intended volume space.}

\item{x_pos}{x voxel coordinate.}

\item{y_pos}{y voxel coordinate.}

\item{z_pos}{z voxel coordinate.}
}
\value{
volume data as a nifti object.
}
\description{
Generate a MRSI voxel from an \code{mrs_data} object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varpro_basic.R
\name{varpro_basic_opts}
\alias{varpro_basic_opts}
\title{Return a list of options for a basic VARPRO analysis.}
\usage{
varpro_basic_opts(method = "fd_re", nnls = TRUE)
}
\arguments{
\item{method}{one of "td", "fd", "fd_re".}

\item{nnls}{restrict basis amplitudes to non-negative values.}
}
\value{
full list of options.
}
\description{
Return a list of options for a basic VARPRO analysis.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{set_ref}
\alias{set_ref}
\title{Set the ppm reference value (eg ppm value at 0Hz).}
\usage{
set_ref(mrs_data, ref)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{ref}{reference value for ppm scale.}
}
\description{
Set the ppm reference value (eg ppm value at 0Hz).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{mat2mrs_data}
\alias{mat2mrs_data}
\title{Convert a matrix (with spectral points in the column dimension and dynamics
in the row dimensions) into a mrs_data object.}
\usage{
mat2mrs_data(
  mat,
  fs = def_fs(),
  ft = def_ft(),
  ref = def_ref(),
  nuc = def_nuc(),
  fd = FALSE
)
}
\arguments{
\item{mat}{data matrix.}

\item{fs}{sampling frequency in Hz.}

\item{ft}{transmitter frequency in Hz.}

\item{ref}{reference value for ppm scale.}

\item{nuc}{resonant nucleus.}

\item{fd}{flag to indicate if the matrix is in the frequency domain (logical).}
}
\value{
mrs_data object.
}
\description{
Convert a matrix (with spectral points in the column dimension and dynamics
in the row dimensions) into a mrs_data object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{re_weighting}
\alias{re_weighting}
\title{Apply a weighting to the FID to enhance spectral resolution.}
\usage{
re_weighting(mrs_data, re, alpha)
}
\arguments{
\item{mrs_data}{data to be enhanced.}

\item{re}{resolution enhancement factor (rising exponential factor).}

\item{alpha}{alpha factor (Guassian decay)}
}
\value{
resolution enhanced mrs_data.
}
\description{
Apply a weighting to the FID to enhance spectral resolution.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{add_noise}
\alias{add_noise}
\title{Add noise to an mrs_data object.}
\usage{
add_noise(mrs_data, sd = 0.1, fd = TRUE)
}
\arguments{
\item{mrs_data}{data to add noise to.}

\item{sd}{standard deviation of the noise.}

\item{fd}{generate the noise samples in the frequency-domain (TRUE) or
time-domain (FALSE). This is required since the absolute value of the
standard deviation of noise samples changes when data is Fourier transformed.}
}
\value{
mrs_data object with additive normally distributed noise.
}
\description{
Add noise to an mrs_data object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{bc_constant}
\alias{bc_constant}
\title{Remove a constant baseline offset based on a reference spectral region.}
\usage{
bc_constant(mrs_data, xlim)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{xlim}{spectral range containing a flat baseline region to measure the
offset.}
}
\value{
baseline corrected data.
}
\description{
Remove a constant baseline offset based on a reference spectral region.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{zero_fade_spec}
\alias{zero_fade_spec}
\title{Fade a spectrum to zero by frequency domain multiplication with a tanh
function. Note this operation distorts data points at the end of the FID.}
\usage{
zero_fade_spec(mrs_data, start_ppm, end_ppm)
}
\arguments{
\item{mrs_data}{data to be faded.}

\item{start_ppm}{start point of the fade in ppm units.}

\item{end_ppm}{end point of the fade in ppm units.}
}
\value{
modified mrs_data object.
}
\description{
Fade a spectrum to zero by frequency domain multiplication with a tanh
function. Note this operation distorts data points at the end of the FID.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{scale_mrs_amp}
\alias{scale_mrs_amp}
\title{Scale an mrs_data object by a scalar or vector or amplitudes.}
\usage{
scale_mrs_amp(mrs_data, amp)
}
\arguments{
\item{mrs_data}{data to be scaled.}

\item{amp}{multiplicative factor, must have length equal to 1 or
Nspec(mrs_data).}
}
\value{
mrs_data object multiplied by the amplitude scale factor.
}
\description{
Scale an mrs_data object by a scalar or vector or amplitudes.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/precomp.R
\name{set_precomp_mode}
\alias{set_precomp_mode}
\title{Set the precompute mode.}
\usage{
set_precomp_mode(mode = NA)
}
\arguments{
\item{mode}{can be one of: "default", "overwrite", "clean" or "disabled".}
}
\description{
Set the precompute mode.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_display.R
\name{gridplot.mrs_data}
\alias{gridplot.mrs_data}
\title{Arrange spectral plots in a grid.}
\usage{
\method{gridplot}{mrs_data}(
  x,
  rows = NA,
  cols = NA,
  mar = c(0, 0, 0, 0),
  oma = c(3.5, 1, 1, 1),
  bty = "o",
  restore_def_par = TRUE,
  ...
)
}
\arguments{
\item{x}{object of class mrs_data.}

\item{rows}{number of grid rows.}

\item{cols}{number of grid columns.}

\item{mar}{option to adjust the plot margins. See ?par.}

\item{oma}{outer margin area.}

\item{bty}{option to draw a box around the plot. See ?par.}

\item{restore_def_par}{restore default plotting par values after the plot has
been made.}

\item{...}{other arguments to pass to the plot method.}
}
\description{
Arrange spectral plots in a grid.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{mask_dyns}
\alias{mask_dyns}
\title{Mask an MRS dataset in the dynamic dimension.}
\usage{
mask_dyns(mrs_data, mask)
}
\arguments{
\item{mrs_data}{MRS data object.}

\item{mask}{vector of boolean values specifying the dynamics to mask, where a
value of TRUE indicates the spectrum should be removed.}
}
\value{
masked dataset.
}
\description{
Mask an MRS dataset in the dynamic dimension.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{comb_coils}
\alias{comb_coils}
\title{Combine coil data based on the first data point of a reference signal.}
\usage{
comb_coils(
  metab,
  ref = NULL,
  noise = NULL,
  scale = TRUE,
  scale_method = "sig_noise_sq",
  sum_coils = TRUE,
  noise_region = c(-0.5, -2.5),
  average_ref_dyns = TRUE,
  ref_pt_index = 1
)
}
\arguments{
\item{metab}{MRS data containing metabolite data.}

\item{ref}{MRS data containing reference data (optional).}

\item{noise}{MRS data from a noise scan (optional).}

\item{scale}{option to rescale coil elements based on the first data point
(logical).}

\item{scale_method}{one of "sig_noise_sq", "sig_noise" or "sig".}

\item{sum_coils}{sum the coil elements as a final step (logical).}

\item{noise_region}{the spectral region (in ppm) to estimate the noise.}

\item{average_ref_dyns}{take the mean of the reference scans in the dynamic
dimension before use.}

\item{ref_pt_index}{time-domain point to use for estimating phase and scaling
values.}
}
\value{
MRS data.
}
\description{
By default, elements are phased and scaled prior to summation. Where a
reference signal is not given, the mean dynamic signal will be used
instead.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{hz}
\alias{hz}
\title{Return the frequency scale of an MRS dataset in Hz.}
\usage{
hz(mrs_data, fs = NULL, N = NULL)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{fs}{sampling frequency in Hz.}

\item{N}{number of data points in the spectral dimension.}
}
\value{
frequency scale.
}
\description{
Return the frequency scale of an MRS dataset in Hz.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_io.R
\name{read_mrs_tqn}
\alias{read_mrs_tqn}
\title{Read MRS data using the TARQUIN software package.}
\usage{
read_mrs_tqn(fname, fname_ref = NA, format, id = NA, group = NA)
}
\arguments{
\item{fname}{the filename containing the MRS data.}

\item{fname_ref}{a second filename containing reference MRS data.}

\item{format}{format of the MRS data. Can be one of the following:
siemens, philips, ge, dcm, dpt, rda, lcm, varian, bruker, jmrui_txt.}

\item{id}{optional ID string.}

\item{group}{optional group string.}
}
\value{
MRS data object.
}
\description{
Read MRS data using the TARQUIN software package.
}
\examples{
fname <- system.file("extdata","philips_spar_sdat_WS.SDAT",package="spant")
\dontrun{
mrs_data <- read_mrs_tqn(fname, format="philips")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{ssp}
\alias{ssp}
\title{Signal space projection method for lipid suppression.}
\usage{
ssp(mrs_data, comps = 5, xlim = c(1.5, 0.8))
}
\arguments{
\item{mrs_data}{MRS data object.}

\item{comps}{the number of spatial components to use.}

\item{xlim}{spectral range (in ppm) covering the lipid signals.}
}
\value{
lipid suppressed \code{mrs_data} object.
}
\description{
Signal space projection method as described in:
Tsai SY, Lin YR, Lin HY, Lin FH. Reduction of lipid contamination in MR
spectroscopy imaging using signal space projection. Magn Reson Med 2019
Mar;81(3):1486-1498.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{vec2mrs_data}
\alias{vec2mrs_data}
\title{Convert a vector into a mrs_data object.}
\usage{
vec2mrs_data(
  vec,
  fs = def_fs(),
  ft = def_ft(),
  ref = def_ref(),
  nuc = def_nuc(),
  dyns = 1,
  fd = FALSE
)
}
\arguments{
\item{vec}{the data vector.}

\item{fs}{sampling frequency in Hz.}

\item{ft}{transmitter frequency in Hz.}

\item{ref}{reference value for ppm scale.}

\item{nuc}{resonant nucleus.}

\item{dyns}{replicate the data across the dynamic dimension.}

\item{fd}{flag to indicate if the matrix is in the frequency domain (logical).}
}
\value{
mrs_data object.
}
\description{
Convert a vector into a mrs_data object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_read_ima.R
\name{read_ima_dyn_dir}
\alias{read_ima_dyn_dir}
\title{Read a directory containing Siemens MRS IMA files and combine along the
dynamic dimension. Note that the coil ID is inferred from the sorted file
name and should be checked when consistency is required.}
\usage{
read_ima_dyn_dir(dir, extra = NULL)
}
\arguments{
\item{dir}{data directory path.}

\item{extra}{an optional data frame to provide additional variables for use
in subsequent analysis steps, eg id or grouping variables.}
}
\value{
mrs_data object.
}
\description{
Read a directory containing Siemens MRS IMA files and combine along the
dynamic dimension. Note that the coil ID is inferred from the sorted file
name and should be checked when consistency is required.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/precomp.R
\name{precomp}
\alias{precomp}
\title{Save function results to file and load on subsequent calls to avoid repeat
computation.}
\usage{
precomp(file, fun, ...)
}
\arguments{
\item{file}{file name to write the results.}

\item{fun}{function to run.}

\item{...}{arguments to be passed to fun.}
}
\description{
Save function results to file and load on subsequent calls to avoid repeat
computation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_display.R
\name{plot_slice_fit}
\alias{plot_slice_fit}
\title{Plot a 2D slice from an MRSI fit result object.}
\usage{
plot_slice_fit(
  fit_res,
  map,
  map_denom = NULL,
  slice = 1,
  zlim = NULL,
  interp = 1
)
}
\arguments{
\item{fit_res}{\code{fit_result} object.}

\item{map}{fit result values to display as a colour map. Can be specified as
a character string or array of numeric values. Defaults to "tNAA".}

\item{map_denom}{fit result values to divide the map argument by. Can be
specified as a character string (eg "tCr") or array of numeric values.}

\item{slice}{slice to plot in the z direction.}

\item{zlim}{range of values to plot.}

\item{interp}{interpolation factor.}
}
\description{
Plot a 2D slice from an MRSI fit result object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_acq_paras}
\alias{get_acq_paras}
\title{Return acquisition parameters from a MRS data object.}
\usage{
get_acq_paras(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\value{
list of acquisition parameters.
}
\description{
Return acquisition parameters from a MRS data object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{beta2lw}
\alias{beta2lw}
\title{Covert a beta value in the time-domain to an equivalent linewidth in Hz:
x * exp(-i * t * t * beta).}
\usage{
beta2lw(beta)
}
\arguments{
\item{beta}{beta damping value.}
}
\value{
linewidth value in Hz.
}
\description{
Covert a beta value in the time-domain to an equivalent linewidth in Hz:
x * exp(-i * t * t * beta).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{array2mrs_data}
\alias{array2mrs_data}
\title{Convert a 7 dimensional array in into a mrs_data object. The array dimensions
should be ordered as : dummy, X, Y, Z, dynamic, coil, FID.}
\usage{
array2mrs_data(
  data_array,
  fs = def_fs(),
  ft = def_ft(),
  ref = def_ref(),
  nuc = def_nuc(),
  fd = FALSE
)
}
\arguments{
\item{data_array}{7d data array.}

\item{fs}{sampling frequency in Hz.}

\item{ft}{transmitter frequency in Hz.}

\item{ref}{reference value for ppm scale.}

\item{nuc}{nucleus that is resonant at the transmitter frequency.}

\item{fd}{flag to indicate if the matrix is in the frequency domain (logical).}
}
\value{
mrs_data object.
}
\description{
Convert a 7 dimensional array in into a mrs_data object. The array dimensions
should be ordered as : dummy, X, Y, Z, dynamic, coil, FID.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{is_fd}
\alias{is_fd}
\title{Check if the chemical shift dimension of an MRS data object is in the
frequency domain.}
\usage{
is_fd(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\value{
logical value.
}
\description{
Check if the chemical shift dimension of an MRS data object is in the
frequency domain.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{mvifftshift}
\alias{mvifftshift}
\title{Perform an ifftshift on a matrix, with each column replaced by its shifted
result.}
\usage{
mvifftshift(x)
}
\arguments{
\item{x}{matrix input.}
}
\value{
output matrix.
}
\description{
Perform an ifftshift on a matrix, with each column replaced by its shifted
result.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{rep_array_dim}
\alias{rep_array_dim}
\title{Repeat an array over a given dimension.}
\usage{
rep_array_dim(x, rep_dim, n)
}
\arguments{
\item{x}{array.}

\item{rep_dim}{dimension to extend.}

\item{n}{number of times to repeat.}
}
\value{
extended array.
}
\description{
Repeat an array over a given dimension.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{calc_peak_info_vec}
\alias{calc_peak_info_vec}
\title{Calculate the FWHM of a peak from a vector of intensity values.}
\usage{
calc_peak_info_vec(data_pts, interp_f)
}
\arguments{
\item{data_pts}{input vector.}

\item{interp_f}{interpolation factor to improve the FWHM estimate.}
}
\value{
a vector of: x position of the highest data point, maximum peak
value in the y axis, FWHM in the units of data points.
}
\description{
Calculate the FWHM of a peak from a vector of intensity values.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{diff_mrs}
\alias{diff_mrs}
\title{Apply the diff operator to an MRS dataset in the FID/spectral dimension.}
\usage{
diff_mrs(mrs_data, ...)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{...}{additional arguments to the diff function.}
}
\value{
MRS data following diff operator.
}
\description{
Apply the diff operator to an MRS dataset in the FID/spectral dimension.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{collapse_to_dyns}
\alias{collapse_to_dyns}
\alias{collapse_to_dyns.mrs_data}
\alias{collapse_to_dyns.fit_result}
\title{Collapse MRS data by concatenating spectra along the dynamic dimension.}
\usage{
collapse_to_dyns(x, rm_masked = FALSE)

\method{collapse_to_dyns}{mrs_data}(x, rm_masked = FALSE)

\method{collapse_to_dyns}{fit_result}(x, rm_masked = FALSE)
}
\arguments{
\item{x}{data object to be collapsed (mrs_data or fit_result object).}

\item{rm_masked}{remove masked dynamics from the output.}
}
\value{
collapsed data with spectra or fits concatenated along the dynamic
dimension.
}
\description{
Collapse MRS data by concatenating spectra along the dynamic dimension.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{comb_fit_list_result_tables}
\alias{comb_fit_list_result_tables}
\title{Combine the fit result tables from a list of fit results.}
\usage{
comb_fit_list_result_tables(fit_list, add_extra = TRUE, add_res_id = TRUE)
}
\arguments{
\item{fit_list}{a list of fit_result objects.}

\item{add_extra}{add variables in the extra data frame to the output (TRUE).}

\item{add_res_id}{add a res_id column to the output to distinguish between
datasets.}
}
\value{
a data frame combine all fit result tables with an additional id
column to differentiate between data sets. Any variables in the extra data
frame may be optionally added to the result.
}
\description{
Combine the fit result tables from a list of fit results.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{mrs_data2mat}
\alias{mrs_data2mat}
\title{Convert mrs_data object to a matrix, with spectral points in the column
dimension and dynamics in the row dimension.}
\usage{
mrs_data2mat(mrs_data, collapse = TRUE)
}
\arguments{
\item{mrs_data}{MRS data object or list of MRS data objects.}

\item{collapse}{collapse all other dimensions along the dynamic dimension, eg
a 16x16 MRSI grid would be first collapsed across 256 dynamic scans.}
}
\value{
MRS data matrix.
}
\description{
Convert mrs_data object to a matrix, with spectral points in the column
dimension and dynamics in the row dimension.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qm_simulation.R
\name{qn_states}
\alias{qn_states}
\title{Get the quantum coherence matrix for a spin system.}
\usage{
qn_states(sys)
}
\arguments{
\item{sys}{spin system object.}
}
\value{
quantum coherence number matrix.
}
\description{
Get the quantum coherence matrix for a spin system.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varpro.R
\name{varpro_opts}
\alias{varpro_opts}
\title{Return a list of options for VARPRO based fitting.}
\usage{
varpro_opts(
  nstart = 20,
  init_g_damping = 2,
  maxiters = 200,
  max_shift = 5,
  max_g_damping = 5,
  max_ind_damping = 5,
  anal_jac = TRUE,
  bl_smth_pts = 80
)
}
\arguments{
\item{nstart}{position in the time-domain to start fitting, units of data
points.}

\item{init_g_damping}{starting value for the global Gaussian line-broadening
term - measured in Hz.}

\item{maxiters}{maximum number of levmar iterations to perform.}

\item{max_shift}{maximum shift allowed to each element in the basis set,
measured in Hz.}

\item{max_g_damping}{maximum permitted global Gaussian line-broadening.}

\item{max_ind_damping}{maximum permitted Lorentzian line-broadening for each
element in the basis set, measured in Hz.}

\item{anal_jac}{option to use the analytic or numerical Jacobian (logical).}

\item{bl_smth_pts}{number of data points to use in the baseline smoothing
calculation.}
}
\value{
list of options.
}
\description{
Return a list of options for VARPRO based fitting.
}
\examples{
varpro_opts(nstart = 10)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_head_dyns}
\alias{get_head_dyns}
\title{Return the first scans of a dynamic series.}
\usage{
get_head_dyns(mrs_data, n = 1)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}

\item{n}{the number of dynamic scans to return.}
}
\value{
first scans of a dynamic series.
}
\description{
Return the first scans of a dynamic series.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{read_lcm_coord}
\alias{read_lcm_coord}
\title{Read an LCModel formatted coord file containing fit information.}
\usage{
read_lcm_coord(coord_f)
}
\arguments{
\item{coord_f}{path to the coord file.}
}
\value{
list containing a table of fit point and results structure containing
signal amplitudes, errors and fitting diagnostics.
}
\description{
Read an LCModel formatted coord file containing fit information.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{calc_spec_diff}
\alias{calc_spec_diff}
\title{Calculate the sum of squares differences between two mrs_data objects.}
\usage{
calc_spec_diff(mrs_data, ref = NULL, xlim = c(4, 0.5))
}
\arguments{
\item{mrs_data}{mrs_data object.}

\item{ref}{reference mrs_data object to calculate differences.}

\item{xlim}{spectral limits to perform calculation.}
}
\value{
an array of the sum of squared difference values.
}
\description{
Calculate the sum of squares differences between two mrs_data objects.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{max_mrs_interp}
\alias{max_mrs_interp}
\title{Apply the max operator to an interpolated MRS dataset.}
\usage{
max_mrs_interp(mrs_data, interp_f = 4)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{interp_f}{interpolation factor.}
}
\value{
Array of maximum values (real only).
}
\description{
Apply the max operator to an interpolated MRS dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{resample_voi}
\alias{resample_voi}
\title{Resample a VOI to match a target image space using nearest-neighbour
interpolation.}
\usage{
resample_voi(voi, mri)
}
\arguments{
\item{voi}{volume data as a nifti object.}

\item{mri}{image data as a nifti object.}
}
\value{
volume data as a nifti object.
}
\description{
Resample a VOI to match a target image space using nearest-neighbour
interpolation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{get_voi_seg}
\alias{get_voi_seg}
\title{Return the white matter, gray matter and CSF composition of a volume.}
\usage{
get_voi_seg(voi, mri_seg)
}
\arguments{
\item{voi}{volume data as a nifti object.}

\item{mri_seg}{segmented brain volume as a nifti object.}
}
\value{
a vector of partial volumes expressed as percentages.
}
\description{
Return the white matter, gray matter and CSF composition of a volume.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{zf_xy}
\alias{zf_xy}
\title{Zero-fill MRSI data in the k-space x-y direction.}
\usage{
zf_xy(mrs_data, factor = 2)
}
\arguments{
\item{mrs_data}{MRSI data.}

\item{factor}{zero-filling factor, a factor of 2 returns a dataset with
twice the original points in the x-y directions. Factors smaller than one
are permitted, such that a factor of 0.5 returns half the k-space points in
the x-y directions.}
}
\value{
zero-filled data.
}
\description{
Zero-fill MRSI data in the k-space x-y direction.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svs_batch_fit.R
\name{svs_1h_brain_analysis}
\alias{svs_1h_brain_analysis}
\title{Standard SVS 1H brain analysis pipeline.}
\usage{
svs_1h_brain_analysis(
  metab,
  basis = NULL,
  w_ref = NULL,
  mri_seg = NULL,
  mri = NULL,
  output_dir = NULL,
  extra = NULL,
  decimate = NULL,
  rats_corr = TRUE,
  ecc = FALSE,
  comb_dyns = TRUE,
  hsvd_filt = FALSE,
  scale_amps = TRUE,
  te = NULL,
  tr = NULL,
  preproc_only = FALSE,
  method = "ABFIT",
  opts = NULL
)
}
\arguments{
\item{metab}{filepath or mrs_data object containing MRS metabolite data.}

\item{basis}{basis set object to use for analysis.}

\item{w_ref}{filepath or mrs_data object containing MRS water reference data.}

\item{mri_seg}{filepath or nifti object containing segmented MRI data.}

\item{mri}{filepath or nifti object containing anatomical MRI data.}

\item{output_dir}{directory path to output fitting results.}

\item{extra}{data.frame with one row containing additional information to be
attached to the fit results table.}

\item{decimate}{option to decimate the input data by a factor of two. The
default value of NULL does not perform decimation unless the spectral width
is greater than 20 PPM.}

\item{rats_corr}{option to perform rats correction, defaults to TRUE.}

\item{ecc}{option to perform water reference based eddy current correction,
defaults to FALSE.}

\item{comb_dyns}{option to combine dynamic scans, defaults to TRUE.}

\item{hsvd_filt}{option to apply hsvd water removal, defaults to FALSE.}

\item{scale_amps}{option to scale metabolite amplitude estimates, defaults to
TRUE.}

\item{te}{metabolite mrs data echo time in seconds.}

\item{tr}{metabolite mrs data repetition time in seconds.}

\item{preproc_only}{only perform the preprocessing steps and omit fitting.
The preprocessed metabolite data will be returned in this case.}

\item{method}{analysis method to use, see fit_mrs help.}

\item{opts}{options to pass to the analysis method.}
}
\value{
a fit_result or mrs_data object depending on the preproc_only option.
}
\description{
Standard SVS 1H brain analysis pipeline.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{read_tqn_fit}
\alias{read_tqn_fit}
\title{Reader for csv fit results generated by TARQUIN.}
\usage{
read_tqn_fit(fit_f)
}
\arguments{
\item{fit_f}{TARQUIN fit file.}
}
\value{
A data frame of the fit data points.
}
\description{
Reader for csv fit results generated by TARQUIN.
}
\examples{
\dontrun{
fit <- read_tqn_fit(system.file("extdata","fit.csv",package="spant"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{sim_noise}
\alias{sim_noise}
\title{Simulate an mrs_data object containing simulated Gaussian noise.}
\usage{
sim_noise(
  sd = 0.1,
  fs = def_fs(),
  ft = def_ft(),
  N = def_N(),
  ref = def_ref(),
  dyns = 1,
  fd = TRUE
)
}
\arguments{
\item{sd}{standard deviation of the noise.}

\item{fs}{sampling frequency in Hz.}

\item{ft}{transmitter frequency in Hz.}

\item{N}{number of data points in the spectral dimension.}

\item{ref}{reference value for ppm scale.}

\item{dyns}{number of dynamic scans to generate.}

\item{fd}{return data in the frequency-domain (TRUE) or time-domain (FALSE)}
}
\value{
mrs_data object.
}
\description{
Simulate an mrs_data object containing simulated Gaussian noise.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_display.R
\name{plot.fit_result}
\alias{plot.fit_result}
\title{Plot the fitting results of an object of class \code{fit_result}.}
\usage{
\method{plot}{fit_result}(
  x,
  dyn = 1,
  x_pos = 1,
  y_pos = 1,
  z_pos = 1,
  coil = 1,
  xlim = NULL,
  data_only = FALSE,
  label = NULL,
  plot_sigs = NULL,
  n = NULL,
  sub_bl = FALSE,
  mar = NULL,
  restore_def_par = TRUE,
  ylim = NULL,
  y_scale = FALSE,
  show_grid = TRUE,
  grid_nx = NULL,
  grid_ny = NA,
  ...
)
}
\arguments{
\item{x}{fit_result object.}

\item{dyn}{the dynamic index to plot.}

\item{x_pos}{the x index to plot.}

\item{y_pos}{the y index to plot.}

\item{z_pos}{the z index to plot.}

\item{coil}{the coil element number to plot.}

\item{xlim}{the range of values to display on the x-axis, eg xlim = c(4,1).}

\item{data_only}{display only the processed data (logical).}

\item{label}{character string to add to the top left of the plot window.}

\item{plot_sigs}{a character vector of signal names to add to the plot.}

\item{n}{single index element to plot (overrides other indices when given).}

\item{sub_bl}{subtract the baseline from the data and fit (logical).}

\item{mar}{option to adjust the plot margins. See ?par.}

\item{restore_def_par}{restore default plotting par values after the plot has
been made.}

\item{ylim}{range of values to display on the y-axis, eg ylim = c(0,10).}

\item{y_scale}{option to display the y-axis values (logical).}

\item{show_grid}{plot gridlines behind the data (logical). Defaults to TRUE.}

\item{grid_nx}{number of cells of the grid in x and y direction. When NULL
the grid aligns with the tick marks on the corresponding default axis (i.e.,
tickmarks as computed by axTicks). When NA, no grid lines are drawn in the
corresponding direction.}

\item{grid_ny}{as above.}

\item{...}{further arguments to plot method.}
}
\description{
Plot the fitting results of an object of class \code{fit_result}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abfit.R
\name{abfit_opts_v1_9_0}
\alias{abfit_opts_v1_9_0}
\title{Return a list of options for an ABfit analysis to maintain comparability with
analyses performed with version 1.9.0 (and earlier) of spant.}
\usage{
abfit_opts_v1_9_0(...)
}
\arguments{
\item{...}{arguments passed to \link{abfit_opts}.}
}
\value{
full list of options.
}
\description{
Return a list of options for an ABfit analysis to maintain comparability with
analyses performed with version 1.9.0 (and earlier) of spant.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{crop_td_pts}
\alias{crop_td_pts}
\title{Crop \code{mrs_data} object data points in the time-domain.}
\usage{
crop_td_pts(mrs_data, start = NULL, end = NULL)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{start}{starting data point (defaults to 1).}

\item{end}{ending data point (defaults to the last saved point).}
}
\value{
cropped \code{mrs_data} object.
}
\description{
Crop \code{mrs_data} object data points in the time-domain.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{zf}
\alias{zf}
\alias{zf.mrs_data}
\alias{zf.basis_set}
\title{Zero-fill MRS data in the time domain.}
\usage{
zf(x, factor = 2)

\method{zf}{mrs_data}(x, factor = 2)

\method{zf}{basis_set}(x, factor = 2)
}
\arguments{
\item{x}{input mrs_data or basis_set object.}

\item{factor}{zero-filling factor, factor of 2 returns a dataset with
twice the original data points.}
}
\value{
zero-filled data.
}
\description{
Zero-fill MRS data in the time domain.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{median_dyns}
\alias{median_dyns}
\title{Calculate the median dynamic data.}
\usage{
median_dyns(mrs_data)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}
}
\value{
median dynamic data.
}
\description{
Calculate the median dynamic data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pulse_sequences.R
\name{seq_cpmg_ideal}
\alias{seq_cpmg_ideal}
\title{CPMG style sequence with ideal pulses.}
\usage{
seq_cpmg_ideal(spin_params, ft, ref, TE = 0.03, echoes = 4)
}
\arguments{
\item{spin_params}{spin system definition.}

\item{ft}{transmitter frequency in Hz.}

\item{ref}{reference value for ppm scale.}

\item{TE}{echo time in seconds.}

\item{echoes}{number of echoes.}
}
\value{
list of resonance amplitudes and frequencies.
}
\description{
CPMG style sequence with ideal pulses.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{Nspec}
\alias{Nspec}
\title{Return the total number of spectra in an MRS dataset.}
\usage{
Nspec(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\description{
Return the total number of spectra in an MRS dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{set_lw}
\alias{set_lw}
\title{Apply line-broadening to an mrs_data object to achieve a specified linewidth.}
\usage{
set_lw(mrs_data, lw, xlim = c(4, 0.5))
}
\arguments{
\item{mrs_data}{data in.}

\item{lw}{target linewidth in units of ppm.}

\item{xlim}{region to search for peaks to obtain a linewidth estimate.}
}
\value{
line-broadened data.
}
\description{
Apply line-broadening to an mrs_data object to achieve a specified linewidth.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_display.R
\name{stackplot.mrs_data}
\alias{stackplot.mrs_data}
\title{Stackplot plotting method for objects of class mrs_data.}
\usage{
\method{stackplot}{mrs_data}(
  x,
  xlim = NULL,
  mode = "re",
  x_units = NULL,
  fd = TRUE,
  col = NULL,
  alpha = NULL,
  x_offset = 0,
  y_offset = 0,
  plot_dim = NULL,
  x_pos = NULL,
  y_pos = NULL,
  z_pos = NULL,
  dyn = 1,
  coil = 1,
  bty = NULL,
  labels = NULL,
  lab_cex = 1,
  right_marg = NULL,
  bl_lty = NULL,
  restore_def_par = TRUE,
  show_grid = NULL,
  grid_nx = NULL,
  grid_ny = NA,
  lwd = NULL,
  ...
)
}
\arguments{
\item{x}{object of class mrs_data.}

\item{xlim}{the range of values to display on the x-axis, eg xlim = c(4,1).}

\item{mode}{representation of the complex numbers to be plotted, can be one
of: "re", "im", "mod" or "arg".}

\item{x_units}{the units to use for the x-axis, can be one of: "ppm", "hz",
"points" or "seconds".}

\item{fd}{display data in the frequency-domain (default), or time-domain
(logical).}

\item{col}{set the colour of the line, eg col = rgb(1, 0, 0, 0.5).}

\item{alpha}{set the line transparency, eg alpha = 0.5 is 50\% transparency.
Overrides any transparency levels set by col.}

\item{x_offset}{separate plots in the x-axis direction by this value.
Default value is 0.}

\item{y_offset}{separate plots in the y-axis direction by this value.}

\item{plot_dim}{the dimension to display on the y-axis, can be one of: "dyn",
"x", "y", "z", "coil" or NULL. If NULL (the default) all spectra will be
collapsed into the dynamic dimension and displayed.}

\item{x_pos}{the x index to plot.}

\item{y_pos}{the y index to plot.}

\item{z_pos}{the z index to plot.}

\item{dyn}{the dynamic index to plot.}

\item{coil}{the coil element number to plot.}

\item{bty}{option to draw a box around the plot. See ?par.}

\item{labels}{add labels to each data item.}

\item{lab_cex}{label size.}

\item{right_marg}{change the size of the right plot margin.}

\item{bl_lty}{linetype for the y = 0 baseline trace. A default value NULL
results in no baseline being plotted.}

\item{restore_def_par}{restore default plotting par values after the plot has
been made.}

\item{show_grid}{plot gridlines behind the data (logical). Defaults to TRUE.}

\item{grid_nx}{number of cells of the grid in x and y direction. When NULL
the grid aligns with the tick marks on the corresponding default axis (i.e.,
tickmarks as computed by axTicks). When NA, no grid lines are drawn in the
corresponding direction.}

\item{grid_ny}{as above.}

\item{lwd}{plot linewidth.}

\item{...}{other arguments to pass to the matplot method.}
}
\description{
Stackplot plotting method for objects of class mrs_data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{mvfftshift}
\alias{mvfftshift}
\title{Perform a fftshift on a matrix, with each column replaced by its shifted
result.}
\usage{
mvfftshift(x)
}
\arguments{
\item{x}{matrix input.}
}
\value{
output matrix.
}
\description{
Perform a fftshift on a matrix, with each column replaced by its shifted
result.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis_set.R
\name{basis2mrs_data}
\alias{basis2mrs_data}
\title{Convert a basis object to an mrs_data object - where basis signals are spread
across the dynamic dimension.}
\usage{
basis2mrs_data(basis, sum_elements = FALSE, amps = NULL, shifts = NULL)
}
\arguments{
\item{basis}{basis set object.}

\item{sum_elements}{return the sum of basis elements (logical)}

\item{amps}{a vector of scaling factors to apply to each basis element.}

\item{shifts}{a vector of frequency shifts (in ppm) to apply to each basis
element.}
}
\value{
an mrs_data object with basis signals spread across the dynamic
dimension or summed.
}
\description{
Convert a basis object to an mrs_data object - where basis signals are spread
across the dynamic dimension.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ift_shift}
\alias{ift_shift}
\title{Perform an iffshift and ifft on a vector.}
\usage{
ift_shift(vec_in)
}
\arguments{
\item{vec_in}{vector input.}
}
\value{
output vector.
}
\description{
Perform an iffshift and ifft on a vector.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_metab}
\alias{get_metab}
\title{Extract the metabolite component from an mrs_data object.}
\usage{
get_metab(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\value{
metabolite component.
}
\description{
Extract the metabolite component from an mrs_data object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{mask_xy}
\alias{mask_xy}
\title{Mask an MRSI dataset in the x-y direction}
\usage{
mask_xy(mrs_data, x_dim, y_dim)
}
\arguments{
\item{mrs_data}{MRS data object.}

\item{x_dim}{x dimension output length.}

\item{y_dim}{y dimension output length.}
}
\value{
masked MRS data.
}
\description{
Mask an MRSI dataset in the x-y direction
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/amp_scaling.R
\name{scale_amp_molal_pvc}
\alias{scale_amp_molal_pvc}
\title{Apply partial volume correction to a fitting result object.}
\usage{
scale_amp_molal_pvc(fit_result, ref_data, p_vols, te, tr, ...)
}
\arguments{
\item{fit_result}{result object generated from fitting.}

\item{ref_data}{water reference MRS data object.}

\item{p_vols}{a numeric vector of partial volumes.}

\item{te}{the MRS TE in seconds.}

\item{tr}{the MRS TR in seconds.}

\item{...}{additional arguments to get_td_amp function.}
}
\value{
A \code{fit_result} object with a rescaled results table.
}
\description{
Apply partial volume correction to a fitting result object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_dyns}
\alias{get_dyns}
\title{Extract a subset of dynamic scans.}
\usage{
get_dyns(mrs_data, subset)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}

\item{subset}{vector containing indices to the dynamic scans to be
returned.}
}
\value{
MRS data containing the subset of requested dynamics.
}
\description{
Extract a subset of dynamic scans.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_ref}
\alias{get_ref}
\title{Extract the reference component from an mrs_data object.}
\usage{
get_ref(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\value{
reference component.
}
\description{
Extract the reference component from an mrs_data object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{get_fh_dyns}
\alias{get_fh_dyns}
\title{Return the first half of a dynamic series.}
\usage{
get_fh_dyns(mrs_data)
}
\arguments{
\item{mrs_data}{dynamic MRS data.}
}
\value{
first half of the dynamic series.
}
\description{
Return the first half of a dynamic series.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{crop_spec}
\alias{crop_spec}
\title{Crop \code{mrs_data} object based on a frequency range.}
\usage{
crop_spec(mrs_data, xlim = c(4, 0.2), scale = "ppm")
}
\arguments{
\item{mrs_data}{MRS data.}

\item{xlim}{range of values to crop in the spectral dimension eg
xlim = c(4, 0.2).}

\item{scale}{the units to use for the frequency scale, can be one of: "ppm",
"hz" or "points".}
}
\value{
cropped \code{mrs_data} object.
}
\description{
Crop \code{mrs_data} object based on a frequency range.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{get_svs_voi}
\alias{get_svs_voi}
\title{Generate a SVS acquisition volume from an \code{mrs_data} object.}
\usage{
get_svs_voi(mrs_data, target_mri)
}
\arguments{
\item{mrs_data}{MRS data.}

\item{target_mri}{optional image data to match the intended volume space.}
}
\value{
volume data as a nifti object.
}
\description{
Generate a SVS acquisition volume from an \code{mrs_data} object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{Npts}
\alias{Npts}
\title{Return the number of data points in an MRS dataset.}
\usage{
Npts(mrs_data)
}
\arguments{
\item{mrs_data}{MRS data.}
}
\value{
number of data points.
}
\description{
Return the number of data points in an MRS dataset.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_reg.R
\name{reslice_to_mrs}
\alias{reslice_to_mrs}
\title{Reslice a nifti object to match the orientation of mrs data.}
\usage{
reslice_to_mrs(mri, mrs)
}
\arguments{
\item{mri}{nifti object to be resliced.}

\item{mrs}{mrs_data object for the target orientation.}
}
\value{
resliced imaging data.
}
\description{
Reslice a nifti object to match the orientation of mrs data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{hsvd_filt}
\alias{hsvd_filt}
\title{HSVD based signal filter.}
\usage{
hsvd_filt(
  mrs_data,
  xlim = c(-30, 30),
  comps = 40,
  irlba = TRUE,
  max_damp = 10,
  scale = "hz",
  return_model = FALSE
)
}
\arguments{
\item{mrs_data}{MRS data to be filtered.}

\item{xlim}{frequency range to filter, default units are Hz which can be
changed to ppm using the "scale" argument.}

\item{comps}{number of Lorentzian components to use for modelling.}

\item{irlba}{option to use irlba SVD (logical).}

\item{max_damp}{maximum allowable damping factor.}

\item{scale}{either "hz" or "ppm" to set the frequency units of xlim.}

\item{return_model}{by default the filtered spectrum is returned. Set
return_model to TRUE to return the HSVD model of the data.}
}
\value{
filtered data or model depending on the return_model argument.
}
\description{
HSVD based signal filter described in:
Barkhuijsen H, de Beer R, van Ormondt D. Improved algorithm for noniterative
and timedomain model fitting to exponentially damped magnetic resonance
signals. J Magn Reson 1987;73:553-557.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrs_data_proc.R
\name{shift}
\alias{shift}
\title{Apply a frequency shift to MRS data.}
\usage{
shift(mrs_data, shift, units = "ppm")
}
\arguments{
\item{mrs_data}{MRS data.}

\item{shift}{frequency shift (in ppm by default).}

\item{units}{of the shift ("ppm" or "hz").}
}
\value{
frequency shifted MRS data.
}
\description{
Apply a frequency shift to MRS data.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitting.R
\name{fit_amps}
\alias{fit_amps}
\title{Extract the fit amplitudes from an object of class \code{fit_result}.}
\usage{
fit_amps(
  x,
  inc_index = FALSE,
  sort_names = FALSE,
  append_common_1h_comb = TRUE
)
}
\arguments{
\item{x}{\code{fit_result} object.}

\item{inc_index}{include columns for the voxel index.}

\item{sort_names}{sort the basis set names alphabetically.}

\item{append_common_1h_comb}{append commonly used 1H metabolite combinations
eg tNAA = NAA + NAAG.}
}
\value{
a dataframe of amplitudes.
}
\description{
Extract the fit amplitudes from an object of class \code{fit_result}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spant.R
\name{set_tqn_cmd}
\alias{set_tqn_cmd}
\title{Set the command to run the TARQUIN command-line program.}
\usage{
set_tqn_cmd(cmd)
}
\arguments{
\item{cmd}{path to binary.}
}
\description{
Set the command to run the TARQUIN command-line program.
}
