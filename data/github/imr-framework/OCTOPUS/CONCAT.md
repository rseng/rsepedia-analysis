---
title: 'Off-resonance CorrecTion OPen soUrce Software (OCTOPUS)'
tags:
  - Python
  - MRI
  - off-resonance correction
authors:
  - name: Marina Manso-Jimeno
    orcid: 0000-0002-1141-2049
    affiliation: "1, 2"
  - name: John Thomas Vaughan Jr.
    orcid: 0000-0002-6933-3757
    affiliation: "1, 2"
  - name: Sairam Geethanath
    orcid: 0000-0002-3776-4114
    affiliation: 2

affiliations:
 - name: Department of Biomedical Engineering, Columbia University in the City of New York, USA
   index: 1
 - name: Columbia Magnetic Resonance Research Center, Columbia University in the City of New York, USA
   index: 2
date: 15 OCTOBER 2020
bibliography: paper.bib
---

# Summary

`OCTOPUS` is a Python-based software for correction of off-resonance
artifacts in Magnetic Resonance (MR) images. It implements three different
methods for correction of both Cartesian and non-Cartesian data: Conjugate Phase Reconstruction (CPR), 
frequency-segmented CPR and Multi-Frequency Interpolation(MFI). `OCTOPUS` is easy to integrate into other two and three-dimensional reconstruction pipelines, which makes the tool highly flexible 
and customizable.

# Statement of need

Off-resonance is an MR artifact which occurs due to field inhomogeneities, differences in tissue 
susceptibilities and chemical shift [@Noll1991]. These phenomena can cause the phase of off-resonant spins to accumulate along the
read-out direction, which can turn into blurring, geometrical distortion
and degradation in the reconstructed image [@LukPat2001]. Images
acquired using long readout trajectories and/or at high fields where the
field homogeneity is lower are more prone to this problem. However,
such acquisition scenarios also deliver desirable properties, such as
short scanning times, gradient efficiency, motion tolerance, and better
signal-to-noise ratio [@Chen2008].

Multiple successful off-resonance correction methods have been reported
in the literature [@Schomberg1999]. Most of them are based on Conjugate
Phase Reconstruction (CPR), a method that counteracts the accumulated
phase by demodulating k-space data with its conjugate [@Maeda1988].
Faster and more efficient implementations that the original CPR 
have been developed, such as frequency-segmented CPR [@Noll1992] and
Multi-Frequency Interpolation (MFI) [@Man1997]. Frequency-segmented CPR reconstructs 
the corrected image by combining the pixels of "L" base images according to each pixel value on a field map. Each base image corresponds to the data demodulated at a fixed frequency, with 
the frequency values for each base image equally spaced within the field map frequency range.
MFI  works in a similar way as frequency-segmented CPR, with main differences being that it 
requires a smaller number of base images (L) and that these images are added together into the corrected image using a set of
linear coefficients derived from the field map. 

One can find optimised off-resonance correction capabilities within
existing packages. Examples are: SPIRiT [@Lustig2010], a MATLAB-based
approach for auto-calibrated parallel imaging reconstruction; Ostenson's
MFI implementation for Magnetic Resonance Fingerprinting (MRF)
[@Ostenson2017]; FUGUE, a tool for Echo-Planar Imaging (EPI) distortion
correction part of the FSL library [@Jenkinson2012]; and the MIRT
toolbox, a MATLAB-based MRI reconstruction package that offers field
inhomogeneity correction using iterative reconstruction
methods [@Sutton2003; @Fessler2005]. Nylund's thesis [@Nylund2014] also
contains source MATLAB code for fs-CPR and MFI correction of spiral
images.

All of these implementations are highly specific, defined for a
particular k-space trajectory or application, and/or include a single
correction method. SPIRiT is devoted to correct data acquired using 
parallel imaging methods; Ostenson's package only corrects MRF spiral data and implements 
only one correction method; and FUGUE corrects distortion solely on EPI images. These limitations typically lead researchers to
adapt their data in an attempt to fit them into the available pipelines
or to write their own version of the methods. Either approach results in
a significant investment of time and effort and can generate isolated
implementations and inconsistent results. Furthermore, most of the
available packages are also MATLAB-based, which unlike Python, requires users to pay a license fee.

`OCTOPUS` is aimed at filling this gap in MR off-resonance correction packages. It provides
Python open-source code for three fundamental methods (CPR, fs-CPR, and
MFI). The implementation is independent of the application and the image
acquisition scheme, easing its integration into any reconstruction
pipeline. `OCTOPUS` can also run in the browser through Google Colab, a freely hosted Jupyter notebook environment that allows one to execute Python code in the browser.
Given this feature, `OCTOPUS` is the first zero-footprint off-resonance
correction software, meaning it doesn't require software download, installation, or configuration on a user's local machine.

# Functionality and limitations
`OCTOPUS` is aimed at MR researchers working with long-readout or field-inhomogeneity sensitive k-space trajectories or 
MR acquisition methods. A short demo is provided in the next section. `OCTOPUS` corrects or reduces geometric distortion and/or blurring present in the images due to off-resonance effects by 
leveraging other Python libraries, specifically NumPy [@2020NumPy], SciPy [@2020SciPy], scikit-image [@scikit-image], 
NiBabel[@Nibabel], Matplotlib [@Matplotlib], OpenCV [@itseez2015opencv], Pydicom [@darcy_mason_2020_3891702], and PyNUFFT[@pynufft]. 
The expected output is an image with recovered, sharper edges and undistorted shape.

Also, `OCTOPUS` corrects off-resonance independently of whether the trajectory used to acquire the data was Cartesian or non-Cartesian. 
The input of the correction methods can be either image or raw data. However, using raw data as input is more efficient
and may avoid non Cartesian trajectory-dependent artifacts. `OCTOPUS` is also able to correct 3D multi-slice and multi-channel data by feeding it to the tool in a slice- and channel-wise manner and then applying channel combination with the user's method of choice.

Presently, the software limitations include correction restricted to data acquired in the absence of 
acceleration techniques, long correction times for large datasets, and degraded correction quality in the presence of highly-inhomogeneous
fields. Additionally, the tool has been only tested on Cartesian, EPI, and spiral data.

# Short demo
To illustrate the usage of the package, we performed in silico numerical
simulations using a single-shot EPI trajectory, a single-shot spiral trajectory and a
simulated field map. For these experiments we used a Shepp-Logan head phantom, which simulates a section of the skull and is widely used
to test reconstruction algorithms [@Shepp1974]. Figure 1 shows all inputs and outputs of the experiment. The steps were:

1. Forward model simulation of off-resonance effect on a 128x128
   Shepp-Logan phantom and 256 mm<sup>2</sup> FOV.

   + Using single-shot EPI and spiral trajectories. Figure 1 shows simplified versions of both trajectories for visualization purposes.

   + Using a  simulated field map based on a blurred version of the phantom image with frequency ranges of -/+ 100, -/+150 and -/+200 Hz.

2. Correction of the results of the forward model with CPR, fs-CPR and MFI .

![Top row (left-right): Shepp-Logan phantom image (128x128), Simplified single-shot EPI k-space trajectory, Simplified single-shot spiral k-space trajectory, and simulated field map (128x128). Bottom row (left-right): EPI experiment results and Spiral experiment results.](JOSS_figs/simfig.png)

In both experiments, 'OCTOPUS' successfully corrected the
off-resonance induced blurring and/or geometrical distortion. Note how the EPI-corrupted images show geometric distortion in the phase-encode direction while spiral corrupted images show blurred and distorted edges.

To test the effect of noise on the correction performance we introduced different levels of noise to a single-shot EPI trajectory-based simulation and measured the peak signal-to-noise ratio (pSNR) and Structural Similarity Index (SSIM). 

![Effect of different noise leves on OCTOPUS correction performance measured using pSNR and SSIM.](JOSS_figs/noise_sim_epi.png)

As expected, PSNR and SSIM are reduced as the off-resonance range widens and the noise level in the original image increases. Nevertheless, in all cases, the three implemented methods improve the metrics with respect to the off-resonance corrupted image.

Finally, to demonstrate the correction capabilities in 3D multi-slice and multi-channel data, we corrected phantom images of a Stack-of-Spirals acquisition with matrix size of 72x72, FOV=240 mm<sup>2</sup> and 54 slices. The images were acquired on a Siemens 3T Prisma scanner using a 20-channel head coil. Figure 3 shows three representative slices and their off-resonance corrected versions. The regions of the images highlighted in red show improved image quality and enhaced edges.

![Off-resonance correction of three slices of a Stack-of-Spirals 3D acquisition.](JOSS_figs/SoS_ORC_v2.png)

# Acknowledgements

This study was funded (in part) by the 'MR Technology Development Grant'
and the 'Seed Grant Program for MR Studies' of the Zuckerman Mind Brain
Behavior Institute at Columbia University (PI: Geethanath) and the 'Fast
Functional MRI with sparse sampling and model-based reconstruction' of
the National Institute of Biomedical Imaging and Bioengineering (PI:
Fessler and, supplement, sub-award to Geethanath).

# References

<p align="center">
<img src="OCTOPUSLogo_Rect.png"/>
</p>

[![DOI](https://zenodo.org/badge/214007114.svg)](https://zenodo.org/badge/latestdoi/214007114)

# OCTOPUS: Off-resonance CorrecTion OPen-soUrce Software
`OCTOPUS` is an open-source tool that provides off-resonance correction  methods for  Magnetic Resonance (MR) images. In particular, the implemented techniques are Conjugate Phase Reconstruction (CPR)[[1]](#references), frequency-segmented CPR [[2]](#references) and Multi-Frequency Interpolation (MFI) [[3]](#references).

Off-resonance is a type of MR image artifact. It originates as an accumulation of phase from off-resonant spins along the read-out direction due to field inhomogeneities, tissue susceptibilities and chemical shift among other possible sources [[4]](#references). Long read-out k-space trajectories are therefore more prone to this artifact and its consequences on the image. The image effects are tipycally blurring and/or geometrical distortion, and consequently, quality deterioration [[5]](#references).

`OCTOPUS` leverages existing techniques and outputs artifact-corrected or mitigated image reconstruction given the raw data from the scanner, k-space trajectory and field map. It is targeted to MR scientists, researchers, engineers and students who work with off-resonance-prone trajectories, such as spirals.

To learn more about the used methods and their implementation visit the [API docs][api-docs].

## Installation
1. Install Python (>=Python 3.6)
2. Create and activate a virtual environment (optional but recommended)
3. Copy and paste this command in your terminal
```pip install MR-OCTOPUS```

**Otherwise, [skip the installation!]** <mark>Run `OCTOPUS` in your browser instead.</mark>

## Quick start
The [Examples folder] contains scripts and data to run off-resonance correction on numerical simulations and phantom images for different k-space trajectories and field maps.

After the [installation] is completed, download the [example data]. Now you can run two types of demos. More information about these and other experiments can be found in the [Wiki page].

### 1. Numerical simulations

`numsim_cartesian.py` and `numsim_spiral.py` run a forward model on a 192x192 Shepp-Logan phantom image. They simulate the off-resonance effect of a cartesian and spiral k-space trajectory, respectively, given a simulated field map.

With `OCTOPUS.fieldmap.simulate` you can experiment the effect of the type of field map and its frequency range on the output corrupted image.

The corrupted image is then corrected using CPR, fs-CPR and MFI and the results are displayed.

### 2. In vitro experiment
If you want to use `OCTOPUS` to correct real data, you can use `ORC_main.py` as a template.
1. Fill the `settings.ini` file with the paths for your inputs and outputs. NOTE: the default settings are configured to run the script using the sample data provided.
2. Input your field of view (FOV), gradient raster time (dt), and echo time (TE).
```python
FOV =   # meters
dt =    # seconds
TE =    # seconds
```
3. Check that the dimensions of your inputs agree.
	`rawdata` dims = `ktraj` dims
5. Specify the number of frequency segments for the fs-CPR and MFI methods
```python
Lx =    # L=Lmin * Lx
```
6. Run the script.
The program will display an image panel with the original image and the corrected versions.

### 3. Command line implementation - **NEW!**
Now you can easily run OCTOPUS using commands on your terminal. After installation, type:
```python
OCTOPUS_cli path/to/container/folder rawdata_file kspace_trajectory_file field_map_file correction_method
```
For more information about the command line implemetation and its required arguments, type:
```python
OCTOPUS_cli -h
```
For more information about how to structure your data, visit the [API docs][api-docs]. 

## Skip the installation! - `OCTOPUS` in your browser

There's no need to go through the installation process. Using this [template][colab-template] you can now run off-resonance correction in your browser!

As a demo, you can use the [example data] provided for the [in vitro experiment].

## Contributing and Community guidelines
`OCTOPUS` adheres to a code of conduct adapted from the [Contributor Covenant] code of conduct.
Contributing guidelines can be found [here][contrib-guidelines].

## References
1. Maeda, A., Sano, K. and Yokoyama, T. (1988), Reconstruction by weighted correlation for MRI with time-varying gradients. IEEE Transactions on Medical Imaging, 7(1): 26-31. doi: 10.1109/42.3926
2. Noll, D. C., Pauly, J. M., Meyer, C. H., Nishimura, D. G. and Macovskj, A. (1992), Deblurring for non‐2D fourier transform magnetic resonance imaging. Magn. Reson. Med., 25: 319-333. doi:10.1002/mrm.1910250210
3. Man, L., Pauly, J. M. and Macovski, A. (1997), Multifrequency interpolation for fast off‐resonance correction. Magn. Reson. Med., 37: 785-792. doi:10.1002/mrm.1910370523
4. Noll, D. C., Meyer, C. H., Pauly, J. M., Nishimura, D. G. and Macovski, A. (1991), A homogeneity correction method for magnetic resonance imaging with time-varying gradients. IEEE Transactions on Medical Imaging, 10(4): 629-637. doi: 10.1109/42.108599
5. Schomberg, H. (1999), Off-resonance correction of MR images. IEEE Transactions on Medical Imaging, 18( 6): 481-495. doi: 10.1109/42.781014

[api-docs]: https://mr-octopus.readthedocs.io/en/latest/
[Contributor Covenant]: http://contributor-covenant.org
[contrib-guidelines]: https://github.com/imr-framework/OCTOPUS/blob/master/CONTRIBUTING.md
[installation]: #installation
[in vitro experiment]: #2-in-vitro-experiment
[Examples folder]: https://github.com/imr-framework/OCTOPUS/tree/master/OCTOPUS/Examples
[example data]: https://github.com/imr-framework/OCTOPUS/blob/master/OCTOPUS/Examples/examples_zip.zip
[colab-template]: https://colab.research.google.com/drive/1hEIj5LaF19yOaWkSqi2uWXyy3u6UgKoP?usp=sharing
[skip the installation!]: #skip-the-installation---octopus-in-your-browser
[Wiki page]: https://github.com/imr-framework/OCTOPUS/wiki/Welcome-to-the-OCTOPUS-wiki!
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

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
reported by contacting the team at <imr-framework2018@gmail.com>. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/# Contributing to `OCTOPUS`
:thumbsup: :tada: Thanks for taking time to contribute! :thumbsup: :tada:

Here are guidelines (not rules!) for contributing to `OCTOPUS`. Use your best judgment, and feel free to propose
changes to this document in a pull request.

## Code of Conduct
This project and everyone participating in it is governed by the
[`OCTOPUS` Code of Conduct](https://github.com/imr-framework/pypulseq/blob/master/CODE_OF_CONDUCT.md).
By participating, you are expected to uphold this code. Please report unacceptable behavior to
[imr.framework2018@gmail.com](mailto:imr.framework2018@gmail.com).

## Pull requests
Follow the coding conventions laid out in the [Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/). Ensure source code is
documented as per the Numpy convention [[numpy1]], [[numpy2]]. If you notice any `OCTOPUS` code not adhering to
[PEP8](https://www.python.org/dev/peps/pep-0008/), submit a pull request or open an issue.

## Issues
Please adhere to the appropriate templates when reporting bugs or requesting features. The templates are automatically
presented via Github's 'New Issue' feature.

[email]: mailto:imr.framework2018@gmail.com
[code_of_conduct]: https://github.com/imr-framework/pypulseq/blob/master/CODE_OF_CONDUCT.md
[style_guide]: https://www.python.org/dev/peps/pep-0008/
[numpy1]: https://numpydoc.readthedocs.io/en/latest/format.html
[numpy2]: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
