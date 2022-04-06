<a href="https://zenodo.org/badge/latestdoi/360498327"><img src="https://zenodo.org/badge/360498327.svg" alt="DOI"></a>

<a style="border-width:0" href="https://doi.org/10.21105/joss.04024"><img src="https://joss.theoj.org/papers/10.21105/joss.04024/status.svg" alt="DOI badge" ></a>

![](https://github.com/SandoghdarLab/PiSCAT/blob/28eafb06ea4f6b468dde70e33fa970d1699974cf/docs/Fig/PiSCAT_logo_bg.png)

# PiSCAT: An open source package in Python for interferometric Scattering Microscopy ([Homepage](https://piscat.readthedocs.io))

iSCAT microscopy was introduced about two decades ago ([1](https://link.aps.org/doi/10.1103/PhysRevLett.93.037401)) and demonstrated to be the method of choice for label-free imaging and tracking of matter at nanometric scale ([2](https://doi.org/10.1021/acs.nanolett.9b01822)), with a wide range of applications such as detection of gold nanoparticles, single dye molecules, viruses, and small proteins ([3](https://en.wikipedia.org/wiki/Interferometric_scattering_microscopy)).
The image of a nanoparticle in iSCAT microscopy is formed via the interference between the light scattered from the particle and a reference field which is a part of the incident laser light. The photostable scattering signal from nanoparticles allows for very long measurements at high speeds, all the way up to megahertz, limited only by the available technology, e.g. of cameras or scanners. Recording fast and long videos however, produces a large volume of data which needs to undergo several stages of computationally demanding analysis. We introduce **PiSCAT** as a python-based package for the analysis of variuos iSCAT measurements and related experiments. PiSCAT aims to facilitate high-performance quantitative analysis of big data and provide a generally open-access platform to enable and speed up the research in iSCAT and related communities. To facilitate the use of PiSCAT, we offer tutorials with live-code features in which we present state-of-the-art algorithms for iSCAT microscopy. These cover important educative materials in [jupyter notebooks](https://jupyter.org/), supported with a web-based [documentation page](https://piscat.readthedocs.io).

In this first release, we provide analysis tools for the sensitive detection of single unlabelled proteins via wide-field iSCAT microscopy. Proteins are only a few nanometers in size with a molecular weight of a few to several hundred kDa. They were detected via iSCAT already in 2014 for small proteins down to the Bovines Serumalbumin (BSA) protein with a mass of 65 kDa ([4](https://doi.org/10.1038/ncomms5495)). iSCAT microscopy is since employed in several more advanced applications such as real-time investigation of cellular secretion ([5](https://doi.org/10.3791/58486),[6](https://doi.org/10.1021/acs.nanolett.7b04494)) and quantitative mass spectrometry of single proteins ([7](https://doi.org/10.1126/science.aar5839)).

## Documentation

The documentation webpage of PiSCAT modules can be found
[here](https://piscat.readthedocs.io).

The outputs from most of the PiSCAT localization and tracking methods are of [Panda data frame type](https://pandas.pydata.org/pandas-docs/stable/reference/frame.html). This data structure has the ability to be easily appended/extended with more information based on different levels of analysis. The data structures containing the results of localization and tracking routines can be saved as csv, mat and HDF5 files. This helps users to work with the analyzed information using different softwares namely, MATLAB and Microsoft Excel. HDF5 is a well-known format that is readable in different programming languages and supports large, complex, heterogeneous data. HDF5 uses a "file directory" like structure that allows users to organize data within the file in structured ways and to embed metadata as well, making it self-describing. 


## Installation

### From PyPi

To install PiSCAT using PyPi, enter the following command in the console:

```
pip install piscat
```

### Local installation of PiSCAT

Clone/download this repository and unzip it. In the project directory enter the following command

```
pip install -e .
```

## Running PiSCAT GUI

Once the installation is done and the python environment is activated, enter the following command in the console:

```
python -m piscat
```

## Running PiSCAT Tutorials

Once the installation is done and the python environment is activated, enter the following command in the console:

```
python -m piscat.Tutorials
```

## Citing PiSCAT

<a style="border-width:0" href="https://doi.org/10.21105/joss.04024"><img src="https://joss.theoj.org/papers/10.21105/joss.04024/status.svg" alt="DOI badge" ></a>


## Testing

To run the tests, please activate the PiSCAT virtual environment. In the project directory, in the console, enter the following command:

```
python setup.py test
```

## Installation of PiSCAT virtual environment in the PyCharm IDE:

1.	Follow the hyper links and the install [ Python 3.9](https://www.python.org/downloads/) and [PyCharm](https://www.jetbrains.com/pycharm/download/#section=windows).
2.	Create a virtual environment based on the instructions provided [here](https://www.jetbrains.com/help/pycharm/creating-virtual-environment.html).
3.  Follow [this link](https://www.jetbrains.com/help/pycharm/creating-and-running-setup-py.html) to select PiSCAT venv as the interpreter, to install the setup.py file and then to run a setup.py task. 

# Contributing

Contributions to PiSCAT are always welcome, and they are greatly appreciated! Our contribution policy can be found [here](https://github.com/SandoghdarLab/PiSCAT/blob/main/CONTRIBUTING.md).

# Contributing

Contributions to PiSCAT are always welcome, and they are greatly appreciated!
A list of open problems can be found [here]( https://github.com/SandoghdarLab/PiSCAT/issues).
Of course, it is also always appreciated to bring own ideas and problems to the community!


Please submit all contributions to the official [Github repository](https://github.com/SandoghdarLab/PiSCAT/) in the form of a Merge Request. Please do not submit git diffs or files containing the changes.

`PiSCAT` is an open-source python package under the license of [GNUv3](https://github.com/SandoghdarLab/PiSCAT/blob/main/LICENSE). Thus we consider the act of contributing to the code by submitting a Merge Request as the "Sign off" or agreement to the GNUv3 license.

You can contribute in many different ways:

## Types of Contributions

### Report Bugs

Report bugs at https://github.com/SandoghdarLab/PiSCAT/issues.

### Fix Issues

Look through the Github issues. Different tags are indicating the status of the issues.
The "bug" tag indicates problems with PiSCAT, while the "enhancement" tag shows ideas that should be added in the future.

### Write Documentation

The documentation of PiSCAT can be found [here](https://piscat.readthedocs.io/). [Jupyter notebooks](https://github.com/SandoghdarLab/PiSCAT/tree/main/piscat/Tutorials/JupyterFiles) and [GUI](https://github.com/SandoghdarLab/PiSCAT/tree/main/piscat/GUI) are used to provide an
interactive start to PiSCAT. It is always appreciated if new document notebooks are provided
since this helps others a lot.

## Get Started!

Ready to contribute? Here is how to set up `PiSCAT` for local development.

1. Fork the `PiSCAT` repo on GitHub.
2. Clone your fork locally:
```bash
    $ git clone https://github.com/USERNAME/PiSCAT.git
    $ cd PiSCAT
```
3. Install your local copy into a virtualenv.
```bash
    $ pip install virtualenvwrapper-win (windows)/ pip install virtualenvwrapper (linux)
    $ mkvirtualenv PiSCAT
    $ pip install -e .
```

By following the comments, you can also install all dependencies with the specific versions that we tested for Python 3.9.7:
```bash
    $ pip install -r requirements.txt
```

4. Create a branch for local development:
```bash
    $ git checkout -b name-of-your-bugfix-or-feature
```
   Now you can make your changes locally.

   To get all packages needed for development, a requirements list can be found [here](https://github.com/SandoghdarLab/PiSCAT/blob/main/setup.py).

5. Commit your changes and push your branch to GitHub::
```bash
    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature
```
6. Submit a Merge Request on Github.

## Merge Request Guidelines

Before you submit a Merge Request, check that it meets these guidelines:

1. All functionality that is implemented through this Merge Request should be covered by unit tests. These are implemented in `PiSCAT\tests`. It would be necessary to add your unit test if the new implementation has some features that are not covered by our unit tests.
2. If the Merge Request adds functionality, the docs should be updated. Put your new functionality into a function with a docstring.
3. If you have a maintainer status for `PiSCAT`, you can merge Merge Requests to the master branch. However, every Merge Request needs to be reviewed by another developer. Thus it is not allowed to merge a Merge Request, which is submitted by oneself.

# Sample data information

For testing different classes and functions in PiSCAT, you can download the following sample video.

Name | size      | shape             | type | subject                      | link
----|-----------|-------------------|------|------------------------------|----
control_4999_128_128_uint16_2.33FPS.raw | 156 (MB)  | (4999, 128, 128)  |uint16| Tutorial [1](https://piscat.readthedocs.io/Tutorial1/Tutorial1.html), [2](https://piscat.readthedocs.io/Tutorial2/Tutorial2.html), [3](https://piscat.readthedocs.io/Tutorial3/Tutorial3.html) |[Download](https://owncloud.gwdg.de/index.php/s/tzRZ7ytBd1weNDl)
5nm_GNPs_128x128_uint16_3333fps_10Acc | 625 (MB)  | (20000, 128, 128) |uint16| Tutorial [4](https://piscat.readthedocs.io/Tutorial4/Tutorial4.html)               |[Download](https://owncloud.gwdg.de/index.php/s/Cq5vU8qIAFIWwEh)
00_darkframes_fullFPS | 40.2 (MB) | (1287, 128, 128)  |uint16| Tutorial [4](https://piscat.readthedocs.io/Tutorial4/Tutorial4.html)               |[Download](https://owncloud.gwdg.de/index.php/s/Cq5vU8qIAFIWwEh)  # How to build the documentation

The documentation uses `Sphinx` including some extensions. 

## Install requirements

```
pip install -r requirements.txt
```

## Build

```
Linux/Mac os --> make clean html
Windows 10 --> 1. make clean, 2. make html
```

Starting point for output is found at 

```
./_build/html/index.html
```

# Tutorial for the correction of the fixed pattern noise in iSCAT images recorded using a sCMOS camera

The weak residual stripes [[1](https://iopscience.iop.org/article/10.1088/1361-6463/ac2f68)] remaining after
Differential Rolling Average ([DRA](https://piscat.readthedocs.io/code_reference.html#piscat.BackgroundCorrection.DifferentialRollingAverage)) are called fixed pattern noise (FPN). A modern CMOS camera typically uses multiple ADCs to improve the imaging speed. The mismatch between the gain and bias of the ADCs [[2](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4014637)] however leaves a FPN on the recorded images. The problem is that the gain and bias parameters temporally fluctuate and therefore FPN is visible even after the subtraction of two consecutive images from each other. The periodicity of such a pattern in our [imaging condition](https://piscat.readthedocs.io/tutorials.html) is about ten pixels which is in the order of the size of a diffraction-limited spot (DLS) of a nano-scatterer. We therefore discuss the effect of FPN on the noise floor behaviour as well as PSF detection sensitivity. We then use a variety of PiSCAT algorithmic approaches based on 
[wavelet transform](https://piscat.readthedocs.io/code_reference.html#piscat.Preproccessing.FrequencyFPNc.update_wFPN) 
[[3](https://www.sciencedirect.com/science/article/abs/pii/S0923596517301522)], [Fourier transform](https://piscat.readthedocs.io/code_reference.html#piscat.Preproccessing.FrequencyFPNc.update_fFPN)[[4](https://www.mdpi.com/1424-8220/18/12/4299)] and [column projection FPN filtering](https://piscat.readthedocs.io/code_reference.html#piscat.Preproccessing.MedianProjectionFPNc) in order to correct for the FPN in the DRA frames.

### Previously ...

In the previous tutorials, we demonstrated how to use PiSCAT's APIs for [setting up the PiSCAT modules and downloading a demo iSCAT video](https://piscat.readthedocs.io/Tutorial1/Tutorial1.html#setting-up-the-piscat-modules-and-downloading-a-demo-iscat-video), [performing basic checks on the acquisition process](https://piscat.readthedocs.io/Tutorial1/Tutorial1.html#examining-the-status-line-removing-it), [suppressing the temporal instability of the laser light](https://piscat.readthedocs.io/Tutorial1/Tutorial1.html#normalization-of-the-power-in-the-frames-of-a-video) and [basic data visualization](https://piscat.readthedocs.io/Tutorial1/Tutorial1.html#display-and-inspect-a-loaded-video). Furthermore, finally we investigated [the noise floor trend](https://piscat.readthedocs.io/Tutorial2/Tutorial2.html#the-effect-of-power-normalization-on-the-detection-limit) along with computing the Differential Rolling Average [DRA](https://piscat.readthedocs.io/Tutorial2/Tutorial2.html#frame-averaging-to-boost-snr-of-imaged-proteins-followed-by-visualization-of-their-signal-via-differential-imaging) of the frames.



```python
# Only to ignore warnings 
import warnings
warnings.filterwarnings('ignore')

# Setting up the path to the PiSCAT modules
import os
import sys
current_path = os.path.abspath(os.path.join('..'))
dir_path = os.path.dirname(current_path)
module_path = os.path.join(dir_path)
data_path = os.path.join(dir_path, 'Tutorials', 'Demo data')#The path to the demo data
if module_path not in sys.path:
    sys.path.append(module_path)
    
# Downloading a control video for this tutorial 
from piscat.InputOutput import download_tutorial_data
download_tutorial_data('control_video')    

# Examining the status line in a loaded/downloaded video and removing the status line, PN+DRA
from piscat.InputOutput import reading_videos
from piscat.Visualization import JupyterDisplay
from piscat.InputOutput import read_status_line
from piscat.Preproccessing import normalization
from piscat.BackgroundCorrection import DifferentialRollingAverage
import numpy as np

df_video = reading_videos.DirectoryType(data_path, type_file='raw').return_df()
paths = df_video['Directory'].tolist()
video_names = df_video['File'].tolist()
demo_video_path = os.path.join(paths[0], video_names[0])#Selecting the first entry in the list
video = reading_videos.video_reader(file_name=demo_video_path, type='binary', img_width=128, img_height=128, 
                                    image_type=np.dtype('<u2'), s_frame=0, e_frame=-1)#Loading the video
status_ = read_status_line.StatusLine(video)#Reading the status line
video_remove_status, status_information  = status_.find_status_line()#Examining the status line & removing it
video_pn, _ = normalization.Normalization(video=video_remove_status).power_normalized()
DRA_PN = DifferentialRollingAverage(video=video_pn, batchSize=120)
RVideo_PN_, gainMap1D_DRA = DRA_PN.differential_rolling(FPN_flag=False, FFT_flag=False)
```

```lang-none    
    ---Status line detected in column---
    
    start power_normalized without parallel loop---> Done
    
    --- start DRA ---
    100%|#########| 4758/4758 [00:00<?, ?it/s]
```


### Median Projection FPN Correction (mFPNc)
To obtain the simplest and easiest correction approach for lowering the additive element of FPN, subtract the median of each column from the corresponding ones in a differential image formed from the difference of two successive batches B 1 and B 2 of raw frames [[1](https://iopscience.iop.org/article/10.1088/1361-6463/ac2f68)].

To optimize the computational performance, FPNc algorithms and [DifferentialRollingAverage](https://piscat.readthedocs.io/code_reference.html#piscat.BackgroundCorrection.DifferentialRollingAverage) have been integrated inside the same class.

![](../Fig/mFPN.png)

```python
DRA_PN_mFPNc = DifferentialRollingAverage(video=video_pn, batchSize=120, mode_FPN='mFPN')
RVideo_PN_mFPNc, gainMap1D_mFPN = DRA_PN_mFPNc.differential_rolling(FPN_flag=True, 
                                                                  select_correction_axis='Both', 
                                                                  FFT_flag=False)
```

```lang-none    
    --- start DRA + mFPN_axis: Both---
    100%|#########| 4758/4758 [00:00<?, ?it/s]

    median FPN correction without parallel loop --->
    100%|#########| 4759/4759 [00:00<?, ?it/s]
    Done
    
    median FPN correction without parallel loop ---> 
    100%|#########| 4759/4759 [00:00<?, ?it/s]

    Done
```    

### mFPN mean-signature

For each correction axis (column or row), the `DifferentialRollingAverage` class returns the corrected video (RVideo_PN_mFPNc) and a set of 1D mean-signatures. These projections assist users in visualizing the impact of FPN. One of these signatures for mFPN correction is plotted in the following cell.


```python
# plotting 1D projection of the first frame for the corrected axis
import matplotlib.pyplot as plt

# For Jupyter notebooks only:
%matplotlib inline
plt.plot(gainMap1D_mFPN[0][500])
plt.xlabel('Pixels')
plt.ylabel('cpFPN projection of DRA')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

```

![](output_7_0.png)

### Column Projection FPN Correction (cpFPNc)

This is a heuristic algorithm in which we first extract the median signature of the column noise by just calculating the median of the image columns as shown in the figure below (subfigure (a)). The image that we work on is a differential image made from two batches B_1 and B_2 of raw frames using the [Differential Rolling Average](https://piscat.readthedocs.io/code_reference.html#piscat.BackgroundCorrection.DifferentialRollingAverage) processing. We note that the median value of pixels in a column will not contain the signal of a particle present in the field of view. This is demonstrated using the blue circles in the figure which encircle the particle's signal in all the steps. The 1D projection plot shows the pixel values with blue dots, the median profile with a red trace and the two upper and lower Median Absolute Deviation (MAD) profiles with black traces. The pixel values corresponding to the blue circle region are above the median and the MAD bands. We then build a binary mask per column, with a threshold value similar to the MAD value. Multiplying the binary mask by the differential image results in a Masked Difference (MD) image. In this manner, we end up partitioning a group of pixels in a column which has values close to their median value, thus omitting the particle's signal. The region marked with the blue circle in the binary mask has zero values and therefore the particle signal is missing in the MD image.

As illustrated in subfigure (b), we compute the column-wise mean of the MD image to obtain the FPN mean-signature with sub-gray level precision. Replicating the 1D mean signature gives us a 2D-FPN map. Finally as shown in subfigure (c), we remove 2D-FPN map from the difference image to get the FPN corrected differential image. 

![](../Fig/cpFPN.png)


```python
DRA_PN_cpFPNc = DifferentialRollingAverage(video=video_pn, batchSize=120, mode_FPN='cpFPN')
RVideo_PN_cpFPNc, gainMap1D_cpFPN = DRA_PN_cpFPNc.differential_rolling(FPN_flag=True, 
                                                                      select_correction_axis='Both', 
                                                                      FFT_flag=False)
```

```lang-none        
    --- start DRA + cpFPN_axis: Both---
    100%|#########| 4758/4758 [00:00<?, ?it/s]

    cpFPN correction without parallel loop ---> 
    100%|#########| 4759/4759 [00:00<?, ?it/s]
    Done
    
    cpFPN correction without parallel loop ---> 
    100%|#########| 4759/4759 [00:00<?, ?it/s]
    Done
```    

### FFT fixed pattern noise correction (fFPNc):


The new method introduced by Zeng, Qingjie, et al.[[3](https://www.mdpi.com/1424-8220/18/12/4299)] separates image structures from column-wise pattern by windowing in the spectral domain. By using an iterative two-stage filtering technique, the algorithm was able to eliminate stripe nonuniformity from coarse to fine. By combining spectral and spatial filtering, this approach preserves image information. The details of this method are depicted in the diagram below.

![](../Fig/fFPN.png)


```python
DRA_PN_wFPN = DifferentialRollingAverage(video=video_pn, batchSize=120, mode_FPN='fFPN')
RVideo_PN_fFPN, gainMap1D_fFPN = DRA_PN_wFPN.differential_rolling(FPN_flag=True, select_correction_axis='Both', 
                                                                 FFT_flag=False, inter_flag_parallel_active=False,
                                                                    max_iterations=30)
```

```lang-none            
    --- start DRA + fFPN_axis: Both---
    100%|#########| 4758/4758 [00:00<?, ?it/s]

    ---start fFPNc without Parallel---
    100%|#########| 4759/4759 [00:00<?, ?it/s]

    ---start fFPNc without Parallel---
    100%|#########| 4759/4759 [00:00<?, ?it/s]
```

### Display and compare different FPNc:

The [JupyterFPNcDisplay](https://piscat.readthedocs.io/code_reference.html#piscat.Visualization.JupyterFPNcDisplay) class in ``display_jupyter``, is developed for visualization of different FPNc videos. For each FPN correction method, images are projected along the columns/rows and plotted with blue dots. For each column, we obtain the mean of those pixel values which is plotted in red.


```python
from piscat.Visualization import JupyterFPNcDisplay
list_videos = [RVideo_PN_, RVideo_PN_mFPNc, RVideo_PN_cpFPNc, RVideo_PN_fFPN]
list_titles = ['DRA_PN', 'mFPNc', 'cpFPNc', 'fFPNc']
%matplotlib inline
JupyterFPNcDisplay(list_videos=list_videos, list_titles=list_titles, correction_axis=1, 
                     numRows=1, numColumns=4, imgSizex=15, imgSizey=15, median_filter_flag=False, color='gray')
```

![](../Fig/tu3_vid1.png)

### FPNc spatial benchmarking

The lateral extent of the microscope PSF is a key spatial feature for detecting nanoparticles, as we discuss it thoroughly in [the protein localization section](https://piscat.readthedocs.io/Tutorial4/Tutorial4.html). The spatial periodicity of the FPN is sometimes in the order of the PSF size. Therefore, in the presence of FPN some PSF-like features are added to the detected signals. Therefore we can also compare the quality of FPN correction simply by considering the number of false particle detections in the blank videos with FPNc being done with different methods but with the same hyperparameters of localization algorithm. The following interactive displays depict mFPNc and cpFPNc that have 8 and 7 false detections while fFPN video has 2 false detection for the same threshold of 6e-5 on the frame number 1000. 


```python
from piscat.Localization import PSFsExtraction
%matplotlib inline
PSF_1 = PSFsExtraction(video=RVideo_PN_mFPNc)
PSFs = PSF_1.psf_detection_preview(function='dog', 
                            min_sigma=1.6, max_sigma=1.7, sigma_ratio=1.1, threshold=6e-5,
                            overlap=0, mode='BOTH', frame_number=[1000], IntSlider_width='400px', 
                                   title='Localization threshold on mFPNc')

PSF_l = PSFsExtraction(video=RVideo_PN_cpFPNc)
PSFs = PSF_l.psf_detection_preview(function='dog', 
                            min_sigma=1.6, max_sigma=1.7, sigma_ratio=1.1, threshold=6e-5,
                            overlap=0, mode='BOTH', frame_number=[1000], IntSlider_width='400px', 
                                   title='Localization threshold on cpFPNc')

PSF_l = PSFsExtraction(video=RVideo_PN_fFPN)
PSFs = PSF_l.psf_detection_preview(function='dog',
                           min_sigma=1.6, max_sigma=1.7, sigma_ratio=1.1, threshold=6e-5,
                           overlap=0, mode='BOTH', frame_number=[1000], IntSlider_width='400px',
                                  title='Localization threshold on fFPNc')

```

![](../Fig/tu3_vid2.png)

### FPNc temporal benchmarking


The noise floor curve illustrates the comparison between mFPNc, cpFPNc and fFPNc for the difference in their performance regarding background temporal fluctuations. FPN comprises both additive and multiplicative terms; mFPN and cpFPN correct the additive component, whereas fFPN corrects FPN in the frequency domain by converting the multiplicative terms to additive terms. The multiplicative part of the FPN is thus adjusted by edge preserving filtering, resulting in an increase in the feasible minimum temporal noise floor. To summarize, fFPN is most appropriate when FPN is caused by a mismatch in ADC gain (multiplicative), but with our camera, FPN is caused by a disparity in ADC offset (additive). In this instance, mFPN or cpFPN are more effective.


```python
# Noise floor analysis

from piscat.BackgroundCorrection import NoiseFloor
l_range = list(range(30, 200, 30))
noise_floor_DRA = NoiseFloor(video_pn, list_range=l_range, FPN_flag=False)

noise_floor_mFPN = NoiseFloor(video_pn, list_range=l_range, select_correction_axis='Both',
                           FPN_flag=True, mode_FPN='mFPN')

noise_floor_cpFPN = NoiseFloor(video_pn, list_range=l_range, select_correction_axis='Both',
                           FPN_flag=True, mode_FPN='cpFPN')

noise_floor_fFPN = NoiseFloor(video_pn, list_range=l_range, select_correction_axis='Both',
                           FPN_flag=True, mode_FPN='fFPN',  max_iterations=10)
```

```python
import matplotlib.pyplot as plt
%matplotlib inline
fig = plt.figure(figsize=(15, 10))
plt.plot(l_range, noise_floor_DRA.mean, label='PN+DRA')
plt.plot(l_range, noise_floor_mFPN.mean, label='PN+DRA+mFPN')
plt.plot(l_range, noise_floor_cpFPN.mean, label='PN+DRA+cpFPN')
plt.plot(l_range, noise_floor_fFPN.mean, label='PN+DRA+fFPN')
plt.xlabel("Batch size", fontsize=18)
plt.ylabel("Noise floor", fontsize=18)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.legend()
``` 

![png](output_21_1.png)

### Bibliography
1. [Mirzaalian Dastjerdi, Houman, et al. "Optimized analysis for sensitive detection and analysis of single proteins via interferometric scattering microscopy." Journal of Physics D: Applied Physics (2021).](http://iopscience.iop.org/article/10.1088/1361-6463/ac2f68)

2. [Snoeij, Martijn F., et al. "A CMOS imager with column-level ADC using dynamic column fixed-pattern noise reduction." IEEE Journal of Solid-State Circuits 41.12 (2006): 3007-3015.](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4014637)

3. [Cao, Yanlong, et al. "A multi-scale non-uniformity correction method based on wavelet decomposition and guided filtering for     uncooled long wave infrared camera." Signal Processing: Image Communication 60 (2018): 13-21.](https://www.sciencedirect.com/science/article/abs/pii/S0923596517301522)

4. [Zeng, Qingjie, et al. "Single infrared image-based stripe nonuniformity correction via a two-stage filtering method." Sensors 18.12 (2018): 4299.](https://www.mdpi.com/1424-8220/18/12/4299)
# Detection & contrast estimation of the proteins in iSCAT videos 
The static version of tutorial documents are presented here. Once the installation of PiSCAT on your local computer is completed, the dynamic version of the tutorial files can be found in the local PiSCAT directory located at `"./Tutorials/JupyterFiles/"`. Based on the number of available CPU cores for parallel processing, this tutorial needs 12-20 GB of computer memory (RAM) to run.
## Previously on PiSCAT tutorials...
Previously, we demonstrated how to use PiSCAT's APIs for 
[setting up the PiSCAT modules and downloading a demo iSCAT video](
https://piscat.readthedocs.io/Tutorial1/Tutorial1.html#setting-up-the-piscat-modules-and-downloading-a-demo-iscat-video), 
[performing basic checks on the acquisition process](https://piscat.readthedocs.io/Tutorial1/Tutorial1.html#examining-the-status-line-removing-it) and 
[basic data visualization](https://piscat.readthedocs.io/Tutorial1/Tutorial1.html#display-and-inspect-a-loaded-video).


```python
# Only to ignore warnings
import warnings
warnings.filterwarnings('ignore')

# Setting up the path to the PiSCAT modules
import os
import sys
current_path = os.path.abspath(os.path.join('..'))
dir_path = os.path.dirname(current_path)
module_path = os.path.join(dir_path)
if module_path not in sys.path:
    sys.path.append(module_path)
     
# Downloading a measurement video for this tutorial 
from piscat.InputOutput import download_tutorial_data
download_tutorial_data('Tutorial3_video')

# Examining the status line in a loaded/downloaded video and removing the line
from piscat.InputOutput import reading_videos
from piscat.Visualization import JupyterDisplay,JupyterSubplotDisplay 
from piscat.InputOutput import read_status_line
from piscat.Preproccessing import normalization
from piscat.BackgroundCorrection import DifferentialRollingAverage
import numpy as np

data_path = os.path.join(dir_path, 'Tutorials', 'Demo data', 'Tutorial3', 'Tutorial3_1')#The path to the measurement data
df_video = reading_videos.DirectoryType(data_path, type_file='raw').return_df()
paths = df_video['Directory'].tolist()
video_names = df_video['File'].tolist()
demo_video_path = os.path.join(paths[0], video_names[0])#Selecting the first entry in the list
video = reading_videos.video_reader(file_name=demo_video_path, type='binary', img_width=128, img_height=128, 
                                    image_type=np.dtype('<u2'), s_frame=0, e_frame=-1)#Loading the video
status_ = read_status_line.StatusLine(video)#Reading the status line
video_remove_status, status_information  = status_.find_status_line()#Examining the status line & removing it
```


```lang-none 
    The directory with the name  Demo data  already exists in the following path: PiSCAT\Tutorials
    Directory  Tutorial3  already exists!
    
    Directory  Histogram  already exists
    ---Status line detected in column---
```


## Dark frame correction
The gray values of an iSCAT image can be directly used to quantitatively estimate the mass of the detected proteins.
 The read-out digital value of a pixel, however, is not entirely built from the collected photons. Part of the recorded signal stems from an extra offset voltage that exists even in imaging under no light condition. It might be helpful to perform an additional preprocessing step on iSCAT videos to subtract a calibrating dark frame from every recorded frame. The dark count in the read-out value of a pixel is a function of acquisition parameters such as exposure time, frame rate and etc. The dark frame correction of iSCAT videos would allow us to obtain the accurate contrasts of the proteins i.e. the true contrast value in
an absolute sense which is independent of the acquisition parameters. In the following cell, we first compute the pixel-wise average of a calibrating dark video to form the mean dark frame and then subtract this frame from the protein measurement videos. 


```python
data_path = os.path.join(dir_path, 'Tutorials', 'Demo data', 'Tutorial3', 'Tutorial3_2')#The path to the dark video data
df_video = reading_videos.DirectoryType(data_path, type_file='raw').return_df()
paths = df_video['Directory'].tolist()
video_names = df_video['File'].tolist()
demo_video_path = os.path.join(paths[0], video_names[0])#Selecting the first entry in the list
video_darkFrames = reading_videos.video_reader(file_name=demo_video_path, type='binary', img_width=128, img_height=128, 
                                    image_type=np.dtype('<u2'), s_frame=0, e_frame=-1)#Loading the video
status_ = read_status_line.StatusLine(video_darkFrames)#Reading the status line
video_darkFrames_remove_status, status_information  = status_.find_status_line()#Examining the status line & removing it
#Computing the mean frame of the dark video
mean_dark_frame = np.mean(video_darkFrames_remove_status, axis=0)
#Alternatively, the mean dark count could also be a good measure of the global offset due to dark counts, given as below,
#mean_dark_frame = np.mean(video_darkFrames_remove_status)
video_remove_status_dc = np.subtract(video_remove_status, mean_dark_frame)#Subtracting the mean dark frame from the measurement 

#Visualization of iSCAT frames before and after correction of dark counts
list_titles=['Raw video', 
             'Video after \ndark frame correction',
            'Difference']

from piscat.Visualization.display_jupyter import JupyterSubplotDisplay
# For Jupyter notebooks only:
%matplotlib inline
JupyterSubplotDisplay(list_videos=[video_remove_status, video_remove_status_dc, 
                                       video_remove_status - video_remove_status_dc], 
                    numRows=1, numColumns=3, list_titles=list_titles, imgSizex=15, imgSizey=5, IntSlider_width='500px',
                    median_filter_flag=False, color='gray')
```


```lang-none 
    ---Status line detected in column---
```


![](../Fig/tu4_vid1.png)


Next, we perform 
[power normalization to suppress the temporal instability of the laser light](
https://piscat.readthedocs.io/Tutorial1/Tutorial1.html#normalization-of-the-power-in-the-frames-of-a-video) and 
[DRA](https://piscat.readthedocs.io/Tutorial2/Tutorial2.html#frame-averaging-to-boost-snr-of-imaged-proteins-followed-by-visualization-of-their-signal-via-differential-imaging) with a batch size of 500 frames.  


```python
#From previous tutorials: power normalization, DRA
video_pn, _ = normalization.Normalization(video=video_remove_status_dc).power_normalized()
video_pn = video_remove_status_dc
batchSize = 500
DRA_PN = DifferentialRollingAverage(video=video_pn, batchSize=batchSize)
RVideo_PN = DRA_PN.differential_rolling(FFT_flag=False)
```

 ```lang-none    
    start power_normalized without parallel loop---> Done
    
    --- start DRA ---
    100%|#########| 18999/18999 [00:00<?, ?it/s]
```

## Localization of proteins [[1](https://iopscience.iop.org/article/10.1088/1361-6463/ac2f68)]:
In this section, we directly work with the dynamic features that are remained in the DRA videos. As mentioned earlier 
the system response of a wide-field microscope for weakly scattering objects can be well approximated with a 2D Gaussian 
function. There exist a variety of localization algorithms available in the localization toolbox of PiSCAT. 
Difference of Gaussian ([DoG](https://piscat.readthedocs.io/code_reference.html#piscat.Localization.PSFsExtraction.psf_detection)) 
algorithm, for example, is suitable to perform a very efficient localization of proteins with pixel precision particle 
localization. [[1](https://iopscience.iop.org/article/10.1088/1361-6463/ac2f68)]

In the following cell, a DRA video is being processed with a suitable set of parameters. The minima and maxima of the sigma values for the DoG kernels are lower and upper limits of the PSF size (in pixels) that one expects once the microscope response function is approximated with a 2D Gaussian function. The sigma ratio and the threshold values in this cell are set with respect to the contrast of the particles we are seeking to detect. We begin this analysis by presenting an interactive PiSCAT class that enables us to tune the DoG detection parameter and visualizes the localized particles dynamically.   


```python
from piscat.Localization import particle_localization

PSF_l = particle_localization.PSFsExtraction(video=RVideo_PN)
PSFs = PSF_l.psf_detection_preview(function='dog',  
                            min_sigma=1.6, max_sigma=1.7, sigma_ratio=1.1, threshold=1.5e-4,
                            overlap=0, mode='BOTH', frame_number=[500, 7385])
```


![](../Fig/tu4_vid2.png)


Once we get to a working set of parameters for our localization algorithm we run the detection algorithm for all the frames of the video. 


```python
PSFs_dog = PSF_l.psf_detection(function='dog', 
                            min_sigma=1.6, max_sigma=1.7, sigma_ratio=1.1, threshold=1.5e-4,
                            overlap=0, mode='BOTH')
```


```lang-none     
    ---start PSF detection with parallel loop---
    100%|#########| 19000/19000 [00:00<?, ?it/s]
```


Detected particles are listed in a Panda data frame named as `PSFs_dog`. The information stored in this data structure is printed in the following,


```python
PSFs_dog
```


<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>y</th>
      <th>x</th>
      <th>frame</th>
      <th>center_intensity</th>
      <th>sigma</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>15.0</td>
      <td>6.0</td>
      <td>0</td>
      <td>-0.001362</td>
      <td>1.6</td>
    </tr>
    <tr>
      <th>1</th>
      <td>30.0</td>
      <td>104.0</td>
      <td>0</td>
      <td>-0.002936</td>
      <td>1.6</td>
    </tr>
    <tr>
      <th>2</th>
      <td>99.0</td>
      <td>11.0</td>
      <td>0</td>
      <td>-0.004625</td>
      <td>1.6</td>
    </tr>
    <tr>
      <th>3</th>
      <td>102.0</td>
      <td>8.0</td>
      <td>0</td>
      <td>-0.001364</td>
      <td>1.6</td>
    </tr>
    <tr>
      <th>4</th>
      <td>116.0</td>
      <td>8.0</td>
      <td>0</td>
      <td>-0.001927</td>
      <td>1.6</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>131851</th>
      <td>46.0</td>
      <td>7.0</td>
      <td>18999</td>
      <td>-0.001506</td>
      <td>1.6</td>
    </tr>
    <tr>
      <th>131852</th>
      <td>105.0</td>
      <td>39.0</td>
      <td>18999</td>
      <td>-0.002903</td>
      <td>1.6</td>
    </tr>
    <tr>
      <th>131853</th>
      <td>109.0</td>
      <td>3.0</td>
      <td>18999</td>
      <td>-0.000172</td>
      <td>1.6</td>
    </tr>
    <tr>
      <th>131854</th>
      <td>114.0</td>
      <td>71.0</td>
      <td>18999</td>
      <td>-0.002101</td>
      <td>1.6</td>
    </tr>
    <tr>
      <th>131855</th>
      <td>118.0</td>
      <td>84.0</td>
      <td>18999</td>
      <td>-0.002118</td>
      <td>1.6</td>
    </tr>
  </tbody>
</table>
<p>131856 rows × 5 columns</p>
</div>


### Filtering of DRA frames prior to the localization of particles with band-pass filters and/or Radial Variance Transform ([2](https://doi.org/10.1364/OE.420670))
Further filtering of DRA frames would facilitate even more robust localization of particles.
This could be, for example, a simple image conditioning routine such as a high-pass Fourier filter which can easily remove large features in the image due to lateral instability of the illumination profile.
In some other cases, while imaging relatively large proteins (bigger than 120KDa), we would have prominent sidelobes in the recorded PSFs. Such proteins have a very strong radially symmetric signature 
that can be exploited in order to tell the particles apart from the other features in the image. In the following cell, we demonstrate the application of [Radial Variance Treansform (RVT) filtering functionality](https://piscat.readthedocs.io/code_reference.html#piscat.Preproccessing.RadialVarianceTransform) on the iSCAT images together with the DoG particle localizer.
As demonstrated in the following cell, with the application of such filters to the images one would need to re-tune the DoG localization parameters.  


```python
# Step1- Localize particles one more time with RVT. 
from piscat.Preproccessing import filtering
rvt_ = filtering.RadialVarianceTransform(inter_flag_parallel_active=False)
filtered_video = rvt_.rvt_video(video=RVideo_PN, rmin=2, rmax=3, kind="basic", highpass_size=None,
            upsample=1, rweights=None, coarse_factor=1, coarse_mode='add',
            pad_mode='constant')
    
# Step2- Just as before we deploy DoG localization algorithm but this time on the filtered DRA video 
PSF_localized_rvt = particle_localization.PSFsExtraction(video=filtered_video, flag_transform=True)
PSFs_RVT_dog = PSF_localized_rvt.psf_detection(function='dog',  
                            min_sigma=1.6, max_sigma=1.7, sigma_ratio=1.1, threshold=1.2e-8,
                            overlap=0, mode='BOTH')

# Step3- Visualization of Localized proteins in the DRA videos with DoG with and without RVT filtering on the left and right correspondingly  
from piscat.Visualization import JupyterPSFs_subplotLocalizationDisplay
JupyterPSFs_subplotLocalizationDisplay(list_videos=[RVideo_PN, filtered_video], list_df_PSFs=[PSFs_dog, PSFs_RVT_dog], 
                                    numRows=1, numColumns=2, list_titles=['DoG', 'RVT'], median_filter_flag=False, 
                                       color='gray', imgSizex=15, imgSizey=15, IntSlider_width='400px', step=1, value=11931)

# Alternatively: We can perform both the processes of steps 1 and 2 together efficiently:
PSF_l = particle_localization.PSFsExtraction(video=RVideo_PN, flag_transform=False)
PSFs_RVT = PSF_l.psf_detection(function='RVT', 
                            min_radial=2, max_radial=3,  rvt_kind="basic", 
                            highpass_size=None, upsample=1, rweights=None, coarse_factor=1, coarse_mode="add",
                            pad_mode="constant", threshold=1.5e-7)
```


```lang-none     
    ---start RVT without Parallel---
    100%|#########| 19000/19000 [00:00<?, ?it/s]

    ---start PSF detection with parallel loop---
    100%|#########| 19000/19000 [00:00<?, ?it/s]

    ---start PSF detection with parallel loop---
    100%|#########| 19000/19000 [00:00<?, ?it/s]
```


![](../Fig/tu4_vid3_1.png)


### Deploying several 2D localization algorithms
Often times it is beneficial to localize the proteins with several localization algorithms as each have their own advantageous. Here is an example where we improve the PSF localization accuracy to sub-pixel level by using the [2D Gaussian fitting](https://piscat.readthedocs.io/code_reference.html#piscat.Localization.PSFsExtraction.fit_Gaussian2D_wrapper) method and append the estimated Gaussian parameters to the previous data structure. Such nonlinear fitting routines are much slower than localization kernels such as DoG but instead, they provide more information on the lateral shape of the PSFs of the proteins.


```python
PSFs = PSF_l.fit_Gaussian2D_wrapper(PSF_List=PSFs_dog, scale=5, internal_parallel_flag=True)  
PSFs.info()
```


```lang-none   
    ---Fitting 2D gaussian with parallel loop---
    100%|#########| 131856/131856 [00:00<?, ?it/s]

    <class 'pandas.core.frame.DataFrame'>
    RangeIndex: 131856 entries, 0 to 131855
    Data columns (total 18 columns):
     #   Column                Non-Null Count   Dtype  
    ---  ------                --------------   -----  
     0   y                     131856 non-null  float64
     1   x                     131856 non-null  float64
     2   frame                 131856 non-null  float64
     3   center_intensity      131856 non-null  float64
     4   sigma                 131856 non-null  float64
     5   Sigma_ratio           131856 non-null  float64
     6   Fit_Amplitude         124098 non-null  float64
     7   Fit_X-Center          124098 non-null  float64
     8   Fit_Y-Center          124098 non-null  float64
     9   Fit_X-Sigma           124098 non-null  float64
     10  Fit_Y-Sigma           124098 non-null  float64
     11  Fit_Bias              124098 non-null  float64
     12  Fit_errors_Amplitude  123550 non-null  float64
     13  Fit_errors_X-Center   123965 non-null  float64
     14  Fit_errors_Y-Center   123929 non-null  float64
     15  Fit_errors_X-Sigma    123734 non-null  float64
     16  Fit_errors_Y-Sigma    123759 non-null  float64
     17  Fit_errors_Bias       123813 non-null  float64
    dtypes: float64(18)
    memory usage: 18.1 MB
```    


## Tracking proteins
Since we perform rolling average analysis before forming differential images, a landing event of a protein would span over multiple frames. So far the proteins are individually localized per frame. In this step, the localization events will be linked together in order to obtain protein trajectories by employing the particle linking routines from [trackpy](http://soft-matter.github.io/trackpy/v0.4.2/index.html). [Linking](https://piscat.readthedocs.io/code_reference.html#piscat.Trajectory.Linking) between the localized proteins is formed when two PSFs are not further than 2 pixels apart (`search_range`) and not more than 10 frames distant temporally (`Memory`).


```python
from piscat.Trajectory.particle_linking import Linking

linking_ = Linking()
linked_PSFs = linking_.create_link(psf_position=PSFs, search_range=2, memory=10)

print("Number of Particles {}".format(linking_.trajectory_counter(linked_PSFs)))
```


```lang-none   
    Frame 18999: 10 trajectories present.
    Number of Particles 422
```
    

## Spatio-temporal filtering [[1](https://iopscience.iop.org/article/10.1088/1361-6463/ac2f68)]
Transient perturbations in the experimental setup may lead to a quasi-speckle like background fluctuation in 
DRA videos. Some of the features from these vibrations could be wrongly identified to be proteins. To identify the 
true particles and keep track of them reliably, we need to have a closer look at the spatial and temporal behaviour of 
each of the protein trajectory candidates.

**Spatial filtering**: In PiSCAT we have a [SpatialFilter](https://piscat.readthedocs.io/code_reference.html#piscat.Localization.SpatialFilter) class that enables users to filter [outlier_frames](https://piscat.readthedocs.io/code_reference.html#piscat.Localization.SpatialFilter.outlier_frames) which suffer from a sudden strong vibration or particle flying by, [dense_PSFs ](https://piscat.readthedocs.io/code_reference.html#piscat.Localization.SpatialFilter.dense_PSFs) and [non-symmetric_PSFs](https://piscat.readthedocs.io/code_reference.html#piscat.Localization.SpatialFilter.symmetric_PSFs) that may not properly resemble the iPSF that one expects from the experimental setup. All of these filters have the threshold parameter that defines the sensitivity of each filter.

**Temporal filtering**:  Since for the detection of very small proteins (ones with few tens of KDa mass), we perform DRA with certain batch size, the signal of the particle is smeared over the frames in which averaging has been done.  In other words, a candidate protein trajectory should have a length comparable to the length of the batches in DRA. In the following, we plot the histogram of protein trajectory lengths in the first cell. We also calculate the median of the trajectory length in order to find out the correct temporal length to be used as the thresholding value in the second cell where we deploy functionalities from the [TemporalFilter](https://piscat.readthedocs.io/code_reference.html#piscat.Trajectory.TemporalFilter) class. 


```python
#Histogram of temporal length of particle trajectories
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

his_all_particles = linked_PSFs['particle'].value_counts()
# For Jupyter notebooks only:
%matplotlib inline
fig = plt.figure(figsize=(12, 4))

plt.hist(his_all_particles, bins=50, fc='C1', ec='k', density=False)
plt.ylabel('#Counts', fontsize=18)
plt.xlabel('Length of linking', fontsize=18)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))  
plt.title('Linking', fontsize=18)
print('Median of linking length is {}'.format(np.median(his_all_particles)))
plt.show()
```

```lang-none   
    Median of linking length is 280.0
```


![](output_20_1.png)


Spatial filtering of PSFs and temporal filtering of linked trajectories are done and the resultant data structure is printed below. 


```python
from piscat.Localization import localization_filtering

# Spatial filters
spatial_filters = localization_filtering.SpatialFilter()
PSFs_filtered = spatial_filters.outlier_frames(linked_PSFs, threshold=20)
PSFs_filtered = spatial_filters.dense_PSFs(PSFs_filtered, threshold=0)
PSFs_filtered = spatial_filters.symmetric_PSFs(PSFs_filtered, threshold=0.7)

# Temporal filters
from piscat.Trajectory import TemporalFilter

t_filters = TemporalFilter(video=RVideo_PN, batchSize=batchSize)
all_trajectories, linked_PSFs_filter, his_all_particles = t_filters.v_trajectory(df_PSFs=PSFs_filtered, threshold=270)
    
# Printing results
PSFs.info()
```


```lang-none   
   
    start removing crappy frames ---> Done!
    
    ---Cleaning the df_PSFs that have drift without parallel loop---
    100%|#########| 18359/18359 [00:00<?, ?it/s]

    Number of PSFs before filters = 115112
    
    Number of PSFs after filters = 105370
    
    start V_trajectories without parallel loop---> 

    100%|#########| 140/140 [00:00<?, ?it/s]

    Done
    
    <class 'pandas.core.frame.DataFrame'>
    RangeIndex: 131856 entries, 0 to 131855
    Data columns (total 18 columns):
     #   Column                Non-Null Count   Dtype  
    ---  ------                --------------   -----  
     0   y                     131856 non-null  float64
     1   x                     131856 non-null  float64
     2   frame                 131856 non-null  float64
     3   center_intensity      131856 non-null  float64
     4   sigma                 131856 non-null  float64
     5   Sigma_ratio           131856 non-null  float64
     6   Fit_Amplitude         124098 non-null  float64
     7   Fit_X-Center          124098 non-null  float64
     8   Fit_Y-Center          124098 non-null  float64
     9   Fit_X-Sigma           124098 non-null  float64
     10  Fit_Y-Sigma           124098 non-null  float64
     11  Fit_Bias              124098 non-null  float64
     12  Fit_errors_Amplitude  123550 non-null  float64
     13  Fit_errors_X-Center   123965 non-null  float64
     14  Fit_errors_Y-Center   123929 non-null  float64
     15  Fit_errors_X-Sigma    123734 non-null  float64
     16  Fit_errors_Y-Sigma    123759 non-null  float64
     17  Fit_errors_Bias       123813 non-null  float64
    dtypes: float64(18)
    memory usage: 18.1 MB
```


### Visualization of detected proteins before and after spatiotemporal filtering


```python
from piscat.Visualization import JupyterPSFs_subplotLocalizationDisplay
JupyterPSFs_subplotLocalizationDisplay(list_videos=[RVideo_PN, RVideo_PN], list_df_PSFs=[PSFs_dog, linked_PSFs_filter], 
                                        numRows=1, numColumns=2, 
                                        list_titles=['Before Spatiotemporal filtering', 'After Spatiotemporal filtering'], 
                                        median_filter_flag=False, color='gray', imgSizex=15, imgSizey=15, 
                                        IntSlider_width='400px', step=1, value=0)
```


![](../Fig/tu4_vid4.png)


## Saving analysis results


```python
from piscat.InputOutput import read_write_data
saving_directory = os.path.join(dir_path, 'Tutorials', 'Demo data','Histogram')
# Saving the filtered particles as an csv file
read_write_data.save_df2csv(linked_PSFs_filter, path=saving_directory, name='filtered_particles')
# Saving all the analysis done for the detected proteins as mat file to be used in matlab for further analysis
read_write_data.save_mat(all_trajectories, path=saving_directory, name='all_trajectories')
# Saving all the analysis results as HDF5 file to be used for example at a later stage for further analysis in PiSCAT
read_write_data.save_list_to_hdf5(all_trajectories, path=saving_directory, name='histogram_data')
```

## Loading results
In the following cell, we provide an example for loading previously saved results of data type of hdf5.


```python
loading_directory = os.path.join(dir_path, 'Tutorials', 'Demo data','Histogram')
# Listing all the files which have HDF5 in the above directory
df_video = reading_videos.DirectoryType(loading_directory, type_file='h5').return_df()
paths = df_video['Directory'].tolist()
HDF5_names = df_video['File'].tolist()
# Choosing the first entry in the list to load 
hist_data_path = os.path.join(paths[0], HDF5_names[0])

from piscat.InputOutput import read_write_data
all_trajectories = read_write_data.load_dict_from_hdf5(hist_data_path)
```


## Estimation of the protein contrast [[1](https://iopscience.iop.org/article/10.1088/1361-6463/ac2f68)]
Once the trajectory of a protein landing or take-off event is built from the individual localization events, the central intensity value of the protein signal in each frame of the DRA video is extracted. This temporal trace can be used to estimate the contrast of the proteins. In the following, we provide two cells to demonstrate this analysis. The contrast estimation of a protein is illustrated in the first cell using the temporal intensity trace of the protein. The same protein is marked and visualized in the second cell where this overlay is done on the DRA videos. 

In the [contrast estimation class](https://piscat.readthedocs.io/code_reference.html#piscat.Analysis.PlotProteinHistogram), we first extend the trace from the detected region (sandwiched between two vertical dashed lines) to grow to twice the batch size if possible. This extended trace is then smoothed using a windowing average which corresponds to 5% of the length of the trace. The mean of the signal at the detected region is computed. A positive mean value would mean the extremum is a maximum and vice versa. We then make sure that the extremum is always a peak from which we separate left and right arms of the signal. Each side of the arms is then separately fitted with a line. The intersection point of these lines can be used to read off the contrast value of the protein in addition to simply reading the peak value of the signal. In case the baseline of the intensity profiles are not symmetric or not zero-valued, one can take Prominence of the extrema as a measure of the contrast which is shown here with a vertical orange line.


```python
# The contrast estimation cell
from piscat.Analysis import PlotProteinHistogram
his_ = PlotProteinHistogram(intersection_display_flag=True, imgSizex=10, imgSizey=5)
his_.plot_contrast_extraction(particles=all_trajectories, batch_size=batchSize, video_frame_num=RVideo_PN.shape[0], MinPeakWidth=100,
                              MinPeakProminence=0, pixel_size=0.66, particles_num='#16')
```


![](output_30_0.png)
    

```python
# Marking and visualizing the detected proteins in the DRA videos, here the particle with ID number 16
from piscat.Visualization import JupyterSelectedPSFs_localizationDisplay

# For Jupyter notebooks only:
%matplotlib inline
JupyterSelectedPSFs_localizationDisplay(video=RVideo_PN, particles=all_trajectories, particles_num='#16', 
                                          frame_extend=0, median_filter_flag=True, 
                                          flag_fit2D=False, color='gray', imgSizex=5, imgSizey=5)
```


![](../Fig/tu4_vid5.png)


## Histogram of the protein contrasts [[1](https://iopscience.iop.org/article/10.1088/1361-6463/ac2f68)]
In the following cell, the distribution of the contrasts of the proteins 
(dark, bright and total) which were previously estimated using 
[three different methods of fitting, peak and prominence](https://piscat.readthedocs.io/Tutorial3/Tutorial3.html#estimation-of-the-protein-contrast) are visualized in the histograms using the ([PlotProteinHistogram](https://piscat.readthedocs.io/code_reference.html#piscat.Analysis.PlotProteinHistogram)) module. Here, we employ the Gaussian Mixture Model (GMM) as a well-established method for identifying the modes or components in a population as well as their features [[3]](https://www.annualreviews.org/doi/abs/10.1146/annurev-statistics-031017-100325).


```python
from piscat.Analysis import PlotProteinHistogram
# For Jupyter notebooks only:
%matplotlib inline

his_ = PlotProteinHistogram(intersection_display_flag=False)
his_(folder_name='', particles=all_trajectories, batch_size=batchSize, video_frame_num=RVideo_PN.shape[0], 
     MinPeakWidth=200, MinPeakProminence=0, pixel_size=0.66)
his_.plot_histogram(bins=6, upper_limitation=6e-3, lower_limitation=-6e-3, step_range=1e-6, face='g', 
                    edge='k', Flag_GMM_fit=True, max_n_components=3, scale=1e1, external_GMM=False)
```


![](output_33_0.png)
    

### Bibliography 
1. [Mirzaalian Dastjerdi, Houman, et al. "Optimized analysis for sensitive detection and analysis of single proteins via interferometric scattering microscopy." Journal of Physics D: Applied Physics (2021).](http://iopscience.iop.org/article/10.1088/1361-6463/ac2f68)
2. [Kashkanova, Anna D., et al. “Precision single-particle localization using radial variance transform.” Optics Express 29.7 (2021): 11070-11083.](https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-29-7-11070&id=449504)
3. [McLachlan, Geoffrey J., Sharon X. Lee, and Suren I. Rathnayake. "Finite mixture models." Annual review of statistics and its application 6 (2019): 355-378.](https://www.annualreviews.org/doi/abs/10.1146/annurev-statistics-031017-100325)# Differential imaging of averaged iSCAT frames    
The static version of tutorial documents are presented here. Once the installation of PiSCAT on your local computer is completed, the dynamic version of the tutorial files can be found in the local PiSCAT directory located at `"./Tutorials/JupyterFiles/"`.  Based on the number of available CPU cores for parallel 
processing, this tutorial needs 5-7 GB of computer memory (RAM) to run.
## Previously on PiSCAT tutorials...
In the last tutorial, 
we [set up the PiSCAT modules and downloaded a demo iSCAT video](Tutorial1.ipynb#Setting-up-the-PiSCAT-modules-and-downloading-a-demo-iSCAT-video), 
[did some basic checks on the acquisition process](Tutorial1.ipynb#Examining-the-status-line-&-removing-it), 
[suppressed the temporal instability of the laser light](Tutorial1.ipynb#Normalization-of-the-power-in-the-frames-of-a-video) 
and used some of the [basic data visualization](Tutorial1.ipynb#Display-and-inspect-a-loaded-video) 
tools provided in PiSCAT for inspection of the iSCAT videos.

```python
# Only to ignore warnings
import warnings
warnings.filterwarnings('ignore')

# Setting up the path to the PiSCAT modules
import os
import sys
current_path = os.path.abspath(os.path.join('..'))
dir_path = os.path.dirname(current_path)
module_path = os.path.join(dir_path)
if module_path not in sys.path:
    sys.path.append(module_path)

# Downloading a control video for this tutorial 
from piscat.InputOutput import download_tutorial_data
download_tutorial_data('control_video')

# Examining the status line in a loaded/downloaded video and removing the line
from piscat.InputOutput import reading_videos
from piscat.Visualization import JupyterDisplay
from piscat.InputOutput import read_status_line
from piscat.Preproccessing import normalization
import numpy as np

data_path = os.path.join(dir_path, 'Tutorials', 'Demo data')#The path to the demo data
df_video = reading_videos.DirectoryType(data_path, type_file='raw').return_df()
paths = df_video['Directory'].tolist()
video_names = df_video['File'].tolist()
demo_video_path = os.path.join(paths[0], video_names[0])#Selecting the first entry in the list
video = reading_videos.video_reader(file_name=demo_video_path, type='binary', img_width=128, img_height=128, 
                                    image_type=np.dtype('<u2'), s_frame=0, e_frame=-1)#Loading the video
status_ = read_status_line.StatusLine(video)#Reading the status line
video_remove_status, status_information  = status_.find_status_line()#Examining the status line & removing it

# Normalization of the power in the frames of a video
video_pn, _ = normalization.Normalization(video=video_remove_status).power_normalized()
```
```lang-none 
    The directory with the name  Demo data  already exists in the following path: PiSCAT\Tutorials
    
    The data file named  Control  already exists in the following path: PiSCAT\Tutorials\Demo data

    ---Status line detected in column---
    
    start power_normalized without parallel loop---> Done
```

## Frame averaging to boost SNR of imaged proteins, followed by visualization of their signal via differential imaging
The illumination profile and imaged speckles from the coverglass are among static features in iSCAT videos that can 
be removed by subtracting two subsequent frames to obtain a differential image which will only include dynamic 
features. As illustrated in the figure below, these features are new relative to the reference image, 
which is itself being rolled forward. In the calculation of the differential image, each image is the mean 
frame of a batch of $L$ number of camera frames. In order to apply Differential Rolling Average (DRA), an object of the 
class [Differential_Rolling_Average](https://piscat.readthedocs.io/code_reference.html#piscat.BackgroundCorrection.DifferentialRollingAverage) is 
instantiated and deployed [[1](http://iopscience.iop.org/article/10.1088/1361-6463/ac2f68)].

![](../Fig/DRA.png)

```python
#For Jupyter notebooks only:
%matplotlib inline

from piscat.BackgroundCorrection import DifferentialRollingAverage
DRA_PN = DifferentialRollingAverage(video=video_pn, batchSize=200)
RVideo_PN, _ = DRA_PN.differential_rolling(FFT_flag=False)

from piscat.Visualization import JupyterDisplay
JupyterDisplay(RVideo_PN, median_filter_flag=False, color='gray', imgSizex=5, imgSizey=5, IntSlider_width='500px', step=100)
```

```lang-none    
    --- start DRA ---
    100%|#########| 4598/4598 [00:00<?, ?it/s]
```

![](../Fig/tu2_vid1.png)

## The effect of power normalization on the detection limit 
Here, we perform a quantitative analysis of the influence of the laser power fluctuations on the sensitivity limit of our scheme using [noise_floor class](https://piscat.readthedocs.io/code_reference.html#piscat.BackgroundCorrection.NoiseFloor) to analyze the noise floor trend as a function of the batch size [[1](http://iopscience.iop.org/article/10.1088/1361-6463/ac2f68)].


```python
# Noise floor analysis

from piscat.BackgroundCorrection import NoiseFloor
l_range = list(range(30, 200, 30))
noise_floor_DRA = NoiseFloor(video_remove_status, list_range=l_range)
noise_floor_DRA_pn = NoiseFloor(video_pn, list_range=l_range)

import matplotlib.pyplot as plt
%matplotlib inline
plt.plot(l_range, noise_floor_DRA.mean, label='DRA')
plt.plot(l_range, noise_floor_DRA_pn.mean, label='PN+DRA')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.xlabel("Batch size", fontsize=18)
plt.ylabel("Noise floor", fontsize=18)
plt.legend()
plt.show()
```

![](output_7_0.png)

We see about 10% improvement in our detection limit with performing power normalization on top of the differential rolling averaging with the best results obtained when the batch size corresponds to 120 frames.

### Bibliography
1. [Mirzaalian Dastjerdi, Houman, et al. "Optimized analysis for sensitive detection and analysis of single proteins via interferometric scattering microscopy." Journal of Physics D: Applied Physics (2021).](http://iopscience.iop.org/article/10.1088/1361-6463/ac2f68)
# Suppression of the laser light fluctuations in a wide-field iSCAT measurement 
In this tutorial, we normalize the recorded power in each video frame in order to suppress 
the temporal instability of the laser light. By doing so we also touch upon the usage of some 
of the very basic PiSCAT packages such as the  [InputOutput module](https://piscat.readthedocs.io/code_reference.html#piscat-inputoutput) 
which provides functionalities for loading iSCAT videos and performing some basic checks on the acquisition process 
through the recorded meta-data. The [Visualization module](https://piscat.readthedocs.io/code_reference.html#piscat-visualization) 
provides a variety of data visualization tools for inspecting iSCAT videos and presentating the analysis results. The normalization 
of laser light fluctuations is one of the early-stage analysis tools that 
are in the [Preprocessing module](https://piscat.readthedocs.io/code_reference.html#piscat-preproccessing). Based on 
the number of available CPU cores for parallel processing, this tutorial needs 2-3 GB of computer memory (RAM) to run.

## Setting up the PiSCAT modules and downloading a demo iSCAT video


```python
# Only to ignore warnings 
import warnings
warnings.filterwarnings('ignore')

# Setting up the path to the PiSCAT modules
import os
import sys
current_path = os.path.abspath(os.path.join('..'))
dir_path = os.path.dirname(current_path)
module_path = os.path.join(dir_path)
if module_path not in sys.path:
    sys.path.append(module_path)
    
# Downloading a blank video for this tutorial
from piscat.InputOutput import download_tutorial_data
download_tutorial_data('control_video')
```

    The directory with the name  Demo data  already exists in the following path: PiSCAT\Tutorials
    
    The data file named  Control  already exists in the following path: PiSCAT\Tutorials\Demo data

## The binary iSCAT videos in a path and loading videos 
In this section, we tabulate all the image and video files of a certain data type which are available in a particular path and demonstrate how to read one exemplary video. We mainly use PhotonFocus cameras and store its recordings as binary files. Here we work on the video that was downloaded earlier. 


```python
import numpy as np
from piscat.InputOutput import reading_videos

#Setting up the path to a data set of the type 'raw' in a particular path 'data_path'
data_path = os.path.join(dir_path, 'Tutorials', 'Demo data', 'Control')
df_video = reading_videos.DirectoryType(data_path, type_file='raw').return_df()
paths = df_video['Directory'].tolist()
video_names = df_video['File'].tolist()

#Choosing the first entry in the video list and loading it
demo_video_path = os.path.join(paths[0], video_names[0])
video = reading_videos.video_reader(file_name=demo_video_path, type='binary', img_width=128, img_height=128, 
                                    image_type=np.dtype('<u2'), s_frame=0, e_frame=-1)

help(reading_videos.video_reader)#Calling help on an imported module/class to know more about it.

```

    Help on function video_reader in module piscat.InputOutput.reading_videos:
    
    video_reader(file_name, type='binary', img_width=128, img_height=128, image_type=dtype('float64'), s_frame=0, e_frame=-1)
        This is a wrapper that can be used to call various video/image readers.
        
        Parameters
        ----------
        file_name: str
            Path of video and file name, e.g. test.jpg.
        
        type: str
            Define the video/image format to be loaded.
        
                * 'binary': use this flag to load binary
                * 'tif': use this flag to load tif
                * 'avi': use this flag to load avi
                * 'png': use this flag to load png
        
        optional_parameters:
            These parameters are used when video 'bin_type' define as binary.
        
            img_width: int
                 For binary images, it specifies the image width.
        
            img_height: int
                For binary images, it specifies the image height.
        
            image_type: str
                Numpy.dtype('<u2') --> video with uint16 pixels data type
        
                * "i"  (signed) integer, "u" unsigned integer, "f" floating-point
                * "<" active little-endian
                * "1" 8-bit, "2" 16-bit, "4" 32-bit, "8" 64-bit
        
            s_frame: int
                Video reads from this frame. This is used for cropping a video.
        
            e_frame: int
                Video reads until this frame. This is used for cropping a video.
        
        Returns
        -------
        @returns: NDArray
            The video/image

## Display and inspect a loaded video
As mentioned earlier, the [Visualization module](https://piscat.readthedocs.io/code_reference.html#piscat-visualization) consists 
of several classes which provide display functionalities. Some of these classes may have the word `jupyter` in their name, 
for example, `display_jupyter`. The reason behind this is that such a class has functionalities similar to its twin 
class namely `display`, but adjusted to be used in Jupyter notebooks. The median filter flag passed as an argument to 
the display classes can be used to achieve a proper visualization of a video albeit having hot or dead pixels. In order 
to scroll left/right through the video frames, you can use the mouse wheel as well as the keyboard arrows button. The 
last line in these images is the meta-data of the measurement that the PhotonFocus camera records in each frame as the status-line.

```python
#For Jupyter notebooks only:
%matplotlib inline

from piscat.Visualization import JupyterDisplay_StatusLine
JupyterDisplay_StatusLine(video, median_filter_flag=False, color='gray', imgSizex=5, imgSizey=5, IntSlider_width='500px', 
                          step=1)
```

    ---Status line detected in column---

![](../Fig/tu1_vid1.png)

## Examining the status line & removing it
The status of the frame-acquisition process is encoded in the status line at the last row of an image. 
We check out the status line to make sure that all the images are recorded properly in high frame rate measurements. 
In such measurements, the acquisition buffer can overflow and some frames could be missed in the recording process.

```python
from piscat.InputOutput import read_status_line
from piscat.Visualization import JupyterDisplay
from IPython.display import display

status_ = read_status_line.StatusLine(video)
video_remove_status, status_information  = status_.find_status_line()

JupyterDisplay(video_remove_status, median_filter_flag=False, color='gray', imgSizex=5, imgSizey=5, IntSlider_width='500px', 
               step=10)
```

    ---Status line detected in column---

![](../Fig/tu1_vid2.png)

## Normalization of the power in the frames of a video
The [Preprocessing module](https://piscat.readthedocs.io/code_reference.html#piscat-preproccessing) provides 
several normalization techniques. In the following step, we correct for the fluctuations in the laser light 
intensity. The summation of all the pixels in an image is the recorded power $P$ in that frame which is then used to 
form the average frame power in a video through $\overline{P}$. The corresponding normalization subroutine returns 
both the power normalized video and the fluctuations in power given by $P/\overline{P} -1$ [[1](http://iopscience.iop.org/article/10.1088/1361-6463/ac2f68)]. The fluctuating 
trend as demonstrated below is in the order of 1E-3 which is in and above the range of contrasts that we expect 
from single proteins of mass few tens of kDa.


```python
# Normalization of the power in the frames of a video
from piscat.Preproccessing import Normalization
video_pn, power_fluctuation = Normalization(video=video_remove_status).power_normalized()

import matplotlib.pyplot as plt
plt.plot(power_fluctuation, 'b', linewidth=1, markersize=0.5)
plt.xlabel('Frame #', fontsize=18)
plt.ylabel(r"$p / \bar p - 1$", fontsize=18)
plt.title('Intensity fluctuations in the laser beam', fontsize=13)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
```

    Directory  piscat_configuration  Created 

    start power_normalized without parallel loop---> Done

![](output_10_1.png)

## Save a video
Finally, we save an analyzed video again as a binary file in order to demonstrate video writing functionalities of the [InputOutput module](https://piscat.readthedocs.io/code_reference.html#piscat-inputoutput).


```python
from piscat.InputOutput import write_video
write_video.write_binary(dir_path=data_path, file_name='demo_1_output.bin', data=video_remove_status, type='original')
```

    Directory  20211015-133647  Created 

### Bibliography
1. [Mirzaalian Dastjerdi, Houman, et al. "Optimized analysis for sensitive detection and analysis of single proteins via interferometric scattering microscopy." Journal of Physics D: Applied Physics (2021).](http://iopscience.iop.org/article/10.1088/1361-6463/ac2f68)
# Introduction to the tutorials for protein detection analysis

In this set of tutorials, we pick out few exemplary iSCAT videos and demonstrate how PiSCAT functionalities can be used in order to detect single unlabelled proteins via wide-field iSCAT microscopy to begin with and continue all the way to performing quantitative analysis on the detected particles. In doing so, we also present a number of necessary preprocessing tools together with various insightful visualization modules in PiSCAT. 

The demo iSCAT videos that we use here, have been recorded with a [Photonfocus](https://www.photonfocus.com/products/camerafinder/camera/mv1-d1024e-160-cl/)  CMOS camera running at 4 kHz and imaging a Region Of Interest (ROI) of 128x128 pixels which is about 5.6x5.6 µm^2 given a magnification of 242X, i.e., the pixel size on the object side is 43.7 nm. The laser operates at a wavelength of 420 nm and the Numerical Aperture (NA) of the objective lens is 1.4. At the acquisition time, every recorded frame is the average of 10 consecutive frames, so effectively we record 400 fps. As a case study, we detect gold nanoparticles (GNP) as small as 5nm with a concentration of 25 nM which were injected into 0.1 M SA buffer onto an uncoated coverslip. Prior to such measurements, we also record a few calibration videos. For example, before injecting nanoparticles, we record iSCAT images of the empty cuvette filled only with the medium. These videos are called blank or control measurements. We also record camera frames while its shutter is closed, just to see how much signal we get in such dark frames under no incident light condition. These demo videos are however quite large and therefore we give the option to the users to download them when they are needed in the tutorials.

The first processing step tackles Intensity fluctuations as they play a role in the formation of very tiny signals in iSCAT. Mechanical vibrations or thermal drifts, for example, can cause temporal instability in the power of a laser beam. In sensitive iSCAT applications as in single protein detection experiments, a small change in the incident light is picked up by the microscope and can limit the extent to which we could average frames. The well depth of pixels in CMOS cameras is finite and this puts an upper limit to the number of photons one can integrate in a single frame. Thus, averaging high number of frames is necessary in order to achieve a detectable signal for very weak scattering objects with contrast less than 1%. 

In the first tutorial, we work with an iSCAT demo video and normalize the pixel values in each frame to the sum of all pixels in the same frame to create a power normalized video in which the temporal instability of the laser light is suppressed. We import PiSCAT modules, run some basic checks on the acquisition process, suppress the temporal instability of the laser light and use some of the basic data visualization tools in PiSCAT.

We discuss about the fluctuations in the recorded signal (e.g. shot noise) whose origins lies in optics which can be smoothed for example through frame averaging. PiSCAT contains a variety of efficient functionalities to tackle noises of such nature and in the second tutorial these methods will be introduced and applied to iSCAT videos.

In the last tutorial we go through the entire analysis pipeline for the detection of single proteins. After performing the above mentioned preprocessing analysis on the 5nm GNP measurements, DLSs are localized in each frame of the video and then linked together to build up trajectories. We later examine the spatio-temporal characteristics of the candidate particles and filter out the outliers. In the next step, the landing profile of a protein is extracted from the intensity trace of the central pixel of the protein in its entire trajectory. The particle contrast can then be estimated via several methods available in PiSCAT. In case the injected sample consists of more than one type of protein, there could exist a number of sub-populations in the distribution of the detected particles. We finally show how to use clustering routines to distinguish multiple modes in the contrast histograms.
.. PiSCAT documentation master file, created by
   sphinx-quickstart on Mon Jan 25 16:37:05 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: https://zenodo.org/badge/360498327.svg
   :target: https://zenodo.org/badge/latestdoi/360498327

PiSCAT: An open source package in Python for interferometric Scattering Microscopy
==================================================================================

iSCAT microscopy was introduced about two decades ago
(`[1] <https://link.aps.org/doi/10.1103/PhysRevLett.93.037401>`_) and demonstrated to be the method of choice for label-free
imaging and tracking of matter at nanometric scale (`[2] <https://doi.org/10.1021/acs.nanolett.9b01822>`_), with a wide range
of applications such as detection of gold nanoparticles, single dye molecules, viruses, and small proteins
(`[3] <https://en.wikipedia.org/wiki/Interferometric_scattering_microscopy>`_).
The image of a nanoparticle in iSCAT microscopy is formed via the interference between the light scattered from the
particle and a reference field which is a part of the incident laser light. The photostable scattering signal from
nanoparticles allows for very long measurements at high speeds, all the way up to megahertz, limited only by the available
technology, e.g. of cameras or scanners. Recording fast and long videos however, produces a large volume of data
which needs to undergo several stages of computationally demanding analysis. We introduce **PiSCAT** as a python-based
package for the analysis of variuos iSCAT measurements and related experiments.
PiSCAT aims to facilitate high-performance quantitative analysis of big data and provide a generally open-access platform
to enable and speed up the research in iSCAT and related communities. To facilitate the use of PiSCAT, we offer tutorials
with live-code features in which we present state-of-the-art algorithms for iSCAT microscopy. These cover important educative materials
in `jupyter notebooks <https://jupyter.org/>`_, supported with a web-based
`documentation page <https://piscat.readthedocs.io>`_.

In this first release, we provide analysis tools for the sensitive detection of single unlabelled proteins via .
wide-field iSCAT microscopy. Proteins are only a few nanometers in size with a molecular weight of a few to several hundred
kDa. They were detected via iSCAT already in 2014 for small proteins down to the Bovines Serumalbumin (BSA) protein with
a mass of 65 kDa (`[4] <https://doi.org/10.1038/ncomms5495>`_). iSCAT microscopy is since employed in several more
advanced applications such as real-time investigation of cellular secretion (`[5] <https://doi.org/10.3791/58486>`_,
`[6] <https://doi.org/10.1021/acs.nanolett.7b04494>`_) and quantitative mass spectrometry of single proteins (`[7] <https://doi.org/10.1126/science.aar5839>`_).


Documentation
-------------

The documentation webpage of PiSCAT modules can be found
`here <https://piscat.readthedocs.io>`_.

The outputs from most of the PiSCAT localization and tracking methods are of `Panda data frame type <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`_.
This data structure has the ability to be easily appended/extended with more information based on different levels of analysis.
The data structures containing the results of localization and tracking routines can be saved as csv,
mat and HDF5 files. This helps users to work with the analyzed information using different softwares namely,
MATLAB and Microsoft Excel. HDF5 is a well-known format that is readable in different programming languages and supports
large, complex, heterogeneous data. HDF5 uses a "file directory" like structure that allows users to organize data within
the file in structured ways and to embed metadata as well, making it self-describing.


PiSCAT:
-------

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   tutorials
   code_reference
   appendix
   bibliography


Bibliography
------------
1. Lindfors, Kalkbrenner, et al. "Detection and spectroscopy of gold nanoparticles using supercontinuum white light confocal microscopy." Physical review letters 93.3 (2004): 037401.

2. Taylor, Richard W., and Vahid Sandoghdar. "Interferometric scattering microscopy: seeing single nanoparticles and molecules via rayleigh scattering." Nano letters 19.8 (2019): 4827-4835.

3. https://en.wikipedia.org/wiki/Interferometric_scattering_microscopy

4. Piliarik, Marek, and Vahid Sandoghdar. "Direct optical sensing of single unlabelled proteins and super-resolution imaging of their binding sites." Nature communications 5.1 (2014): 1-8.

5. Gemeinhardt, André, et al. "Label-free imaging of single proteins secreted from living cells via iSCAT microscopy." JoVE (Journal of Visualized Experiments) 141 (2018): e58486.

6. McDonald, Matthew P., et al. "Visualizing single-cell secretion dynamics with single-protein sensitivity." Nano letters 18.1 (2018): 513-519.

7. Young, Gavin, et al. "Quantitative mass imaging of single biological macromolecules." Science 360.6387 (2018): 423-427.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Appendix
========

Summary of the PiSCAT pipeline for efficient nanoparticle detection
-------------------------------------------------------------------

The implementation flowchart for detecting single unlabeled proteins from iSCAT videos is shown in the following Fig [1].

.. image:: ./Fig/flowchart.png
  :width: 400
  :alt:  PiSCAT pipeline [1]
  :align: center


PiSCAT publications
-------------------
[1] Dastjerdi, Houman Mirzaalian, et al. "Optimized analysis for sensitive detection and analysis of single proteins via interferometric scattering microscopy." Journal of Physics D: Applied Physics 55.5 (2021): 054002. (`Journal <https://iopscience.iop.org/article/10.1088/1361-6463/ac2f68>`_)













------------
Bibliography
------------

.. bibliography:: references.bib
   :all:
   :style: unsrt=============
API Reference
=============


piscat.Analysis
---------------

.. autoclass:: piscat.Analysis.ReadProteinAnalysis
    :members:

.. autoclass:: piscat.Analysis.PlotProteinHistogram
    :members:

.. autofunction:: piscat.Analysis.protein_analysis


piscat.BackgroundCorrection
---------------------------

.. autoclass:: piscat.BackgroundCorrection.DifferentialRollingAverage
    :members:

.. autoclass:: piscat.BackgroundCorrection.NoiseFloor
    :members:


piscat.InputOutput
------------------

.. autoclass:: piscat.InputOutput.CameraParameters
    :members:

.. autoclass:: piscat.InputOutput.CPUConfigurations
    :members:

.. autoclass:: piscat.InputOutput.Image2Video
    :members:

.. autoclass:: piscat.InputOutput.StatusLine
    :members:

.. autofunction:: piscat.InputOutput.save_mat

.. autofunction:: piscat.InputOutput.read_mat

.. autofunction:: piscat.InputOutput.save_dic_to_hdf5

.. autofunction:: piscat.InputOutput.save_list_to_hdf5

.. autofunction:: piscat.InputOutput.load_dict_from_hdf5

.. autofunction:: piscat.InputOutput.save_df2csv

.. autofunction:: piscat.InputOutput.save_dic2json

.. autofunction:: piscat.InputOutput.read_json2dic

.. autofunction:: piscat.InputOutput.video_reader

.. autofunction:: piscat.InputOutput.read_binary

.. autofunction:: piscat.InputOutput.read_tif

.. autofunction:: piscat.InputOutput.read_avi

.. autofunction:: piscat.InputOutput.read_png

.. autofunction:: piscat.InputOutput.read_fits

.. autoclass:: piscat.InputOutput.DirectoryType
    :members:

.. autofunction:: piscat.InputOutput.write_binary

.. autofunction:: piscat.InputOutput.write_MP4

.. autofunction:: piscat.InputOutput.write_GIF


piscat.Localization
-------------------

.. autoclass:: piscat.Localization.RadialCenter
    :members:

.. autoclass:: piscat.Localization.PSFsExtraction
    :members:

.. autoclass:: piscat.Localization.SpatialFilter
    :members:

.. autoclass:: piscat.Localization.DirectionalIntensity
    :members:

.. autofunction:: piscat.Localization.gaussian_2d

.. autofunction:: piscat.Localization.fit_2D_Gaussian_varAmp

.. autofunction:: piscat.Localization.blob_frst

.. autofunction:: piscat.Localization.feature2df

.. autofunction:: piscat.Localization.list2dataframe

piscat.Preproccessing
---------------------

.. autoclass:: piscat.Preproccessing.FFT2D
    :members:

.. autoclass:: piscat.Preproccessing.Filters
    :members:

.. autoclass:: piscat.Preproccessing.RadialVarianceTransform
    :members:

.. autoclass:: piscat.Preproccessing.FastRadialSymmetryTransform
    :members:

.. autoclass:: piscat.Preproccessing.GuidedFilter
    :members:

.. autoclass:: piscat.Preproccessing.GrayGuidedFilter
    :members:

.. autoclass:: piscat.Preproccessing.MultiDimGuidedFilter
    :members:

.. autoclass:: piscat.Preproccessing.MedianProjectionFPNc
    :members:

.. autoclass:: piscat.Preproccessing.ColumnProjectionFPNc
    :members:

.. autoclass:: piscat.Preproccessing.FrequencyFPNc
    :members:

piscat.Trajectory
-----------------

.. autoclass:: piscat.Trajectory.Linking
    :members:

.. autoclass:: piscat.Trajectory.TemporalFilter
    :members:

.. autofunction:: piscat.Trajectory.protein_trajectories_list2dic


piscat.Visualization
--------------------

.. autoclass:: piscat.Visualization.ContrastAdjustment
    :members:

.. autoclass:: piscat.Visualization.Display
    :members:

.. autoclass:: piscat.Visualization.DisplayDataFramePSFsLocalization
    :members:

.. autoclass:: piscat.Visualization.DisplayPSFs_subplotLocalizationDisplay
    :members:

.. autoclass:: piscat.Visualization.DisplaySubplot
    :members:

.. autoclass:: piscat.Visualization.JupyterDisplay
    :members:

.. autoclass:: piscat.Visualization.JupyterDisplay_StatusLine
    :members:

.. autoclass:: piscat.Visualization.JupyterPSFs_localizationDisplay
    :members:

.. autoclass:: piscat.Visualization.JupyterPSFs_localizationPreviewDisplay
    :members:

.. autoclass:: piscat.Visualization.JupyterPSFs_subplotLocalizationDisplay
    :members:

.. autoclass:: piscat.Visualization.JupyterPSFs_TrackingDisplay
    :members:

.. autoclass:: piscat.Visualization.JupyterSelectedPSFs_localizationDisplay
    :members:

.. autoclass:: piscat.Visualization.JupyterPSFs_2_modality_subplotLocalizationDisplay
    :members:

.. autoclass:: piscat.Visualization.JupyterFPNcDisplay
    :members:

.. autofunction:: piscat.Visualization.plot2df

.. autofunction:: piscat.Visualization.plot3

.. autofunction:: piscat.Visualization.plot_histogram
Introduction to the tutorials for protein detection analysis
============================================================

In this set of tutorials, we pick out few exemplary iSCAT videos and demonstrate how PiSCAT functionalities can be used in order to detect single unlabelled proteins via wide-field iSCAT microscopy to begin with and continue all the way to performing quantitative analysis on the detected particles. In doing so, we also present a number of necessary preprocessing tools together with various insightful visualization modules in PiSCAT.

The demo iSCAT videos that we use here, have been recorded with a `Photonfocus <https://www.photonfocus.com/products/camerafinder/camera/mv1-d1024e-160-cl/>`_  CMOS camera running at 4 kHz and imaging a Region Of Interest (ROI) of 128x128 pixels which is about 5.6x5.6 µm^2 given a magnification of 242X, i.e., the pixel size on the object side is 43.7 nm. The laser operates at a wavelength of 420 nm and the Numerical Aperture (NA) of the objective lens is 1.4. At the acquisition time, every recorded frame is the average of 10 consecutive frames, so effectively we record 400 fps. As a case study, we detect gold nanoparticles (GNP) as small as 5nm with a concentration of 25 nM which were injected into 0.1 M SA buffer onto an uncoated coverslip. Prior to such measurements, we also record a few calibration videos. For example, before injecting nanoparticles, we record iSCAT images of the empty cuvette filled only with the medium. These videos are called blank or control measurements. We also record camera frames while its shutter is closed, just to see how much signal we get in such dark frames under no incident light condition. These demo videos are however quite large and therefore we give the option to the users to download them when they are needed in the tutorials.

The first processing step tackles Intensity fluctuations as they play a role in the formation of very tiny signals in iSCAT. Mechanical vibrations or thermal drifts, for example, can cause temporal instability in the power of a laser beam. In sensitive iSCAT applications as in single protein detection experiments, a small change in the incident light is picked up by the microscope and can limit the extent to which we could average frames. The well depth of pixels in CMOS cameras is finite and this puts an upper limit to the number of photons one can integrate in a single frame. Thus, averaging high number of frames is necessary in order to achieve a detectable signal for very weak scattering objects with contrast less than 1%.

In the first tutorial, we work with an iSCAT demo video and normalize the pixel values in each frame to the sum of all pixels in the same frame to create a power normalized video in which the temporal instability of the laser light is suppressed. We import PiSCAT modules, run some basic checks on the acquisition process, suppress the temporal instability of the laser light and use some of the basic data visualization tools in PiSCAT.

We discuss about the fluctuations in the recorded signal (e.g. shot noise) whose origins lies in optics which can be smoothed for example through frame averaging. PiSCAT contains a variety of efficient functionalities to tackle noises of such nature and in the second tutorial these methods will be introduced and applied to iSCAT videos.

In the last tutorial we go through the entire analysis pipeline for the detection of single proteins. After performing the above mentioned preprocessing analysis on the 5nm GNP measurements, DLSs are localized in each frame of the video and then linked together to build up trajectories. We later examine the spatio-temporal characteristics of the candidate particles and filter out the outliers. In the next step, the landing profile of a protein is extracted from the intensity trace of the central pixel of the protein in its entire trajectory. The particle contrast can then be estimated via several methods available in PiSCAT. In case the injected sample consists of more than one type of protein, there could exist a number of sub-populations in the distribution of the detected particles. We finally show how to use clustering routines to distinguish multiple modes in the contrast histograms. The static version of tutorial documents are presented here.
Once the installation of PiSCAT on your local computer is completed, the dynamic version of the tutorial files can be found in the local PiSCAT directory located at `"./Tutorials/JupyterFiles/"`.





---------
Tutorials
---------

.. toctree::
    :maxdepth: 1

    download_tutorial_data.md
    Tutorial1/Tutorial1.md
    Tutorial2/Tutorial2.md
    Tutorial3/Tutorial3.md
    Tutorial4/Tutorial4.md
Installation
============

From PyPi
---------

To install PiSCAT using PyPi, enter the following command in the console:

``pip install piscat``


Local installation of PiSCAT
----------------------------
Clone/download this repository and unzip it. In the project directory enter the following command:

``pip install -e .``


Running PiSCAT GUI
------------------
Once the installation is done and the python environment is activated, enter the following command in the
console:

``python -m piscat``


Running PiSCAT Tutorials
------------------------
Once the installation is done and the python environment is activated, enter the following command in the console:

``python -m piscat.Tutorials``


Testing
-------
To run the tests, please activate the PiSCAT virtual environment. In the project directory,
in the console, enter the following command:

``python setup.py test``


Installation of PiSCAT virtual environment in the PyCharm IDE:
--------------------------------------------------------------

1.	Follow the hyper links and the install `Python 3.8 <https://www.python.org/downloads/>`_
and `PyCharm <https://www.jetbrains.com/pycharm/download/#section=windows>`_.

2.	Create a virtual environment based on the instructions provided
`here <https://www.jetbrains.com/help/pycharm/creating-virtual-environment.html>`_.

3.  Follow `this link <https://www.jetbrains.com/help/pycharm/creating-and-running-setup-py.html>`_
to select PiSCAT venv as the interpreter, to install the setup.py file and then to run a setup.py task.





