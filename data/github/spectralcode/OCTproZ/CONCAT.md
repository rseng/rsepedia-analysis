 # <img style="vertical-align:middle" img src="images/octproz_icon.png" width="64"> OCTproZ - Vision and long-term goals


General vision
----------
OCTproZ should:
- be an easy to use OCT signal processing and visualization software for researchers that utilize OCT.
- enable everyone to use custom build OCT systems without the hassle of having to write the processing software
- __provide open OCT processing algorithms to make scientific research more reproducible__



Long-term goals
----------
- Implement OCT processing pipeline in other frameworks besides CUDA (e.g. OpenCL, OpenMP, C ++ AMP) to support various computer hardware for processing
- Plugin database for most common OCT system hardware configurations
- Advanced visualization for live view of OCT data (e.g. speckle variance, spectroscopic)
- Python or scripting interface either as plugin or integrated within OCTproZ
 # <img style="vertical-align:middle" img src="images/octproz_icon.png" width="64"> OCTproZ 

OCTproZ is an open source software for optical coherence tomography (OCT) processing and visualization. A plug-in system enables the integration of custom OCT systems and software modules.

<p align="center">
  <img src="images/octproz_screenshot_ubuntu.png" width="640">
</p>

The output windows in the screenshot above show OCT images of a strawberry. 

Video
--------
[Live OCT visualization: Burning an acorn with a laser](https://www.youtube.com/watch?v=d_O83kFXiRA&feature=youtu.be)


Features
--------

* **Real-time OCT processing and visualization with single GPU**  </br>
The full [OCT processing pipeline](https://spectralcode.github.io/OCTproZ/#processing-section) is implemented in [CUDA](https://developer.nvidia.com/cuda-zone) and visualization is performed with [OpenGL](https://www.opengl.org). Depending on the GPU used, OCTproZ can be used for MHz-OCT. 

* **Plug-in system** </br>
Plug-ins enable the integration of custom OCT systems and software modules. There are two kinds of plug-ins for OCTproZ: _Acquisition Systems_ and _Extensions_. An Acquisition System controls the OCT hardware and provides raw data to OCTproZ. Extensions have access to processed OCT data and can be used to extend the functionality of OCTproZ. 

* **Cross platform** </br>
OCTproZ runs on Windows and Linux. </br>
It has been successfully tested on Windows 10, Ubuntu 16.04, Ubuntu 18.04 and JetPack 4.4.1 (Jetson Nano 4 GB)


Processing Pipeline
--------
<p align="center">
  <img src="images/processing_pipeline_linear_v1_1_0.png" >
</p>

A detailed overview of the OCTproZ processing pipeline can be found [here](https://spectralcode.github.io/OCTproZ/#processing-section).

Performance
----------
Performance highly depends on the used computer hardware and the size of the of the OCT data. A test data set with 12 bit per sample, 1024 samples per raw A-scan, 512 A-scans per B-scan and 256 B-scans per volume was used to measure the performance on different systems:

GPU           | A-scan rate without live 3D view | A-scan rate with live 3D view
------------- | ------------- | -------------
NVIDIA Quadro K620  | ~ 300 kHz ( ~2.2 volumes/s) | ~ 250 kHz ( ~1.9 volumes/s)
NVIDIA GeForce GTX 1080 Ti  | ~ 4.8 MHz (~ 36 volumes/s) | ~ 4.0 MHz (~ 30 volumes/s)

You can find more information [here](performance.md).


Plug-ins
----------
The following plug-ins are currently available:
</br></br>
__Acquisition Systems:__
|Name | Description |
|-----|-----|
|[Virtual OCT System](octproz_project/octproz_plugins/octproz_virtual_oct_system)| Can be used to load already acquired OCT raw data from the disk|


__Extensions:__
|Name | Description |
|-----|-----|
|[Demo Extension](octproz_project/octproz_plugins/octproz_demo_extension)| This demo extension is for developers. It has no useful functionality, but the code can be used as a template for developing custom extensions.|
|[Image Statistics](https://github.com/spectralcode/ImageStatisticsExtension)| Displays useful image statistics, such as a histogram, in real time of currently acquired B-scans |
|[Socket Stream](https://github.com/spectralcode/SocketStreamExtension)| Streaming of OCT data via TCP/IP. Just for slow OCT acquisitions.|
|[Phase Extraction](https://github.com/spectralcode/PhaseExtractionExtension)| Can be used to determine a suitable resampling curve for k-linearization.|

The easiest way to develop custom plug-ins is to clone/download the entire OCTproZ project, compile the DevKit and OCTproZ and use the existing examples as templates. Have a look at the [plugin developer guide](https://spectralcode.github.io/OCTproZ/developer.html). 


Download and Installation
----------
To run OCTproZ a cuda-compatible graphics card with current drivers is required.

A precompiled package for Windows (64bit) can be downloaded from:
[GitHub release section](https://github.com/spectralcode/OCTproZ/releases). Extract the zip archive and execute OCTproZ, installation is not necessary.

If you need OCTproZ for a different operating system, the easiest way is to compile it yourself. See the compiling section.

Test Dataset
----------
A test dataset that can be used with the Virtual OCT System can be downloaded from [here](https://figshare.com/articles/SSOCT_test_dataset_for_OCTproZ/12356705). 

User Manual
----------
An online version of the user manual can be found [here](https://spectralcode.github.io/OCTproZ/index.html). 

Developer Guide
----------
The plugin developer guide can be found [here](https://spectralcode.github.io/OCTproZ/developer.html). 

Compiling
---------
Compiling instructions can be found [here](BUILD.md).

Contributing
----------
Contribution guidelines can be found [here](CONTRIBUTING.md).

Long-term goals
----------
Vision and long-term goals can be found [here](vision.md).

Known issues
----------
On some Linux distributions floating dock widgets lose mouse focus when dragged. See: [Qt bug](https://bugreports.qt.io/browse/QTBUG-65640)


Publication
----------
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02580/status.svg)](https://doi.org/10.21105/joss.02580)


BibTeX:
```
@article{Zabic2020,
  doi = {10.21105/joss.02580},
  url = {https://doi.org/10.21105/joss.02580},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {54},
  pages = {2580},
  author = {Miroslav Zabic and Ben Matthias and Alexander Heisterkamp and Tammo Ripken},
  title = {Open Source Optical Coherence Tomography Software},
  journal = {Journal of Open Source Software}
}
```


License
----------
OCTproZ is licensed under [GPLv3](LICENSE).</br>
The DevKit is licensed under [MIT license](octproz_project/octproz_devkit/LICENSE). # <img style="vertical-align:middle" img src="images/octproz_icon.png" width="64"> Building OCTproZ

OCTproZ can be build on Windows and Linux. 

# Compiling
Building OCTproZ from source requires: 
- Installation of [Qt 5](https://www.qt.io/offline-installers) (version 5.10.1 or later)
- Installation of [CUDA Toolkit](https://developer.nvidia.com/cuda-downloads) (version 8 or later)
- __Windows:__ MSVC compiler that is compatible with your CUDA version (see [CUDA installation guide for Windows](https://docs.nvidia.com/cuda/cuda-installation-guide-microsoft-windows/index.html#system-requirements)) To get the MSVC compiler it is the easiest to search online where/how to get it as this changes from time to time. Pay attention that you get the right version of the MSVC compiler as described in the CUDA guide. <br>
__Linux:__ Development environment that is compatible with your CUDA version (see [CUDA installation guide for Linux](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#system-requirements)) and the third-party libraries mentioned in the [CUDA installation guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#install-libraries)

A more detailed Linux guide for the installation steps above can be found in the section "Installing development tools to build OCTproZ on Linux" below. 

How to compile:
1. Clone/Download the OCTproZ project. The destination path should not contain any spaces!
2. Start Qt Creator and open [octproz_project.pro](octproz_project/octproz_project.pro)
3. Configure project by selectig appropriate kit in Qt Creator (on Windows you need the MSVC compiler)
4. Change the CUDA architecture flags in [cuda.pri](octproz_project/octproz/pri/cuda.pri) if necessary for your hardware ([more info](https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/))
5. Build _octproz_project_ (_right click on _octproz_project_  -> Build_)
6. Run OCTproZ (right click on _octproz_project_ __or__ right click on _octproz_ -> Run)

The project file _octproz_project.pro_ ensures that the devkit is compiled first as it generates a folder with files that are used by OCTproZ and any plug-in at compile time.  </br>



# Installing development tools to build OCTproZ on Linux

## Debian based systems
The following instructions have been tested with Ubuntu 18.04.

### 1. Install Qt 5:
Open a terminal (ctrl + alt + t) and type the following commands:
```
sudo apt-get install build-essential
sudo apt-get install qtcreator
sudo apt-get install qt5-default
```

Qt documentation and examples which are not not required, but recommended can be installed with these commands:
```
sudo apt-get install qt5-doc
sudo apt-get install qt5-doc-html qtbase5-doc-html
sudo apt-get install qtbase5-examples
```


### 2. Install CUDA
Follow the [CUDA installation guide for Linux](https://docs.nvidia.com/cuda/cuda-installation-guide-linux) carefully.

When you get to the post-installation actions, you need to know that adding paths to the PATH variable can be done with the terminal-based text editor Nano by modifying the file ".bashrc".

Open .bashrc with Nano:
```
nano /home/$USER/.bashrc
```
Now insert the cuda relevant paths as statet in the cuda installation guide at the end of the file.

To save the file press on your keyboard
```
ctrl + o
```
And to close the file press
```
ctrl + x
```

After this you should verify [that the CUDA toolkit can find and communicate correctly with the CUDA-capable hardware.](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/#verify-installation)

Finally you need to [install some third-party libraries](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/#install-libraries):
```
sudo apt-get install g++ freeglut3-dev build-essential libx11-dev \
    libxmu-dev libxi-dev libglu1-mesa libglu1-mesa-dev
```


That is all! Now you are able to compile OCTproZ by opening the OCTproZ project files with Qt Creator. 


# OCTproZ on the NVIDIA Jetson Nano with JetPack

To compile and run OCTproZ on the Jetson Nano you need Qt with enabled _OpenGL desktop_ option.

For this Qt can be built from source on the Jetson Nano:

### 1. Install build dependencies
Enable the "Source code" option in Software and Updates > Ubuntu Software under the "Downloadable from the Internet" section.

Open a terminal (ctrl + alt + t) an install the build dependencies like this
```
sudo apt-get build-dep qt5-default
```

### 2. Get the Qt source
```
git clone https://code.qt.io/qt/qt5.git
cd qt5
git checkout 5.12.9
```
then
```
git submodule update --init --recursive
cd ~
mkdir qt5-build
cd qt5-build
```

### 3. Configure and build
In this step you configure Qt for _OpenGL desktop_ (this is necessary!) and remove some packages with _-skip_ (this is optional. Removing those packages slightly reduces the build time)

```
../qt5/configure -qt-xcb -opengl desktop -nomake examples -nomake tests -skip qtwebengine -skip qtandroidextras -skip qtcanvas3d -skip qtcharts -skip qtconnectivity -skip qtdatavis3d -skip qtdeclarative -skip qtpurchasing -skip qtquickcontrols -skip qtquickcontrols2 -skip qtwinextras
```

After the configuration was done a _Configure summary_ will be displayed. Please verify that there is a _yes_ in the line with _Desktop OpenGL_. Now you can start the build process:

```
make
sudo make install
```
Be aware  that _make_ will take about 5 hours on the Jetson Nano. 

When everything has been successfully completed you can start Qt Creator and build OCTproZ!

References:
- [wiki.qt.io Building Qt 5 from Git](https://wiki.qt.io/Building_Qt_5_from_Git)
- [stackoverflow.com Qt-default version issue on migration from RPi4 to NVIDIA Jetson Nano](https://stackoverflow.com/questions/62190967/qt-default-version-issue-on-migration-from-rpi4-to-nvidia-jetson-nano)


# Troubleshooting

After installing Qt 5.12.11 with the offline installer, you may get the error message:

```
NMAKE : fatal error U1077: "C:\Program": Rückgabe-Code "0x1""
```

One way to solve this issue is to close Qt Creator and (re-)install __Qt Creator 4.14.2__ with the offline installer form here:
https://download.qt.io/official_releases/qtcreator/4.14/4.14.2/

Then deleted the file _toolchains.xml_. You can find the file here:
```
C:\Users\%USERNAME%\AppData\Roaming\QtProject\qtcreator\toolchains.xml
```

After these steps reopen Qt Creator and everything should work fine. # <img style="vertical-align:middle" img src="images/octproz_icon.png" width="64"> OCTproZ - Contribution Guidelines

There are four main ways to contribute: issue reports, feature requests, code contributions and plug-in development.  


Issue Reports
----------
- Check if the issue you want to report or a similar one was already reported on the [OCTproZ issue page](https://github.com/spectralcode/OCTproZ/issues). Feel free to comment on already existing issue reports to provide additional information. 
- Report new issues by using the "Bug report"-template if applicable.


Feature Requests
----------
- Check if the feature you want to request or a similar one was already requested on the [OCTproZ issue page](https://github.com/spectralcode/OCTproZ/issues). Feel free to comment on already existing feature requests to provide additional information. 
- Request new features by using the "Feature request"-template if applicable.


Code Contributions
----------
- Before creating a [pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests), please announce your plans in the [OCTproZ issue page](https://github.com/spectralcode/OCTproZ/issues). This enables an open discussion which will get the best out of your efforts. 
- Currently there is no style guide. New code should be more or less consistent with the style of already existing code.


Plug-ins
----------
- If you have developed a plug-in for OCTproZ and want it to be mentioned in the [OCTproZ README.md](https://github.com/spectralcode/OCTproZ/blob/master/README.md), please contact us.
- Contact information can be found at the bottom of the [user manual](https://spectralcode.github.io/OCTproZ/index.html)
 # <img style="vertical-align:middle" img src="images/octproz_icon.png" width="64"> OCTproZ - Performance Information

Processing rate highly depends on the size of the raw data, the used computer hardware and resource usage by background or system processes. With modern computer hardware and typical data dimensions for OCT, OCTproZ achieves A-scan rates in the MHz range.

A test data set with 12 bit per sample, 1024 samples per raw A-scan, 512 A-scans per B-scan and 256 B-scans per volume was used to measure the performance on different systems:

| |**Office Computer**|**Lab Computer**|**Gaming Computer**|
|:-----|:-----|:-----|:-----|
|CPU|Intel® Core i5-7500|AMD Ryzen™ Threadripper 1900X|AMD Ryzen™ 5 1600|
|RAM|16 GB|32 GB|16 GB|
|GPU|NVIDIA Quadro K620|NVIDIA GeForce GTX 1080 Ti| NVIDIA GeForce GTX 1080
|Operating system|Windows 10|Ubuntu 16.04| Windows 10|
|A-scan rate with 3D view| ~ 250 kHz (~ 1.9 volumes/s)|~ 4.0 MHz (~ 30 volumes/s)|~ 1.9 MHz (~ 15 volumes/s)|
|A-scan rate without 3D view| ~ 300 kHz (~ 2.2 volumes/s)|~ 4.8 MHz (~ 36 volumes/s)|~ 2.4 MHz (~ 18 volumes/s)|

<br>

| Embedded System |**NVIDIA Jetson Nano**|
|:-----|:-----|
|CPU|ARMv8 Processor rev 1(v8l) x 4|
|RAM|4 GB|
|GPU|NVIDIA Tegra X1 (128-core Maxwell)|
|Operating system|Ubuntu 18.04 (JetPack 4.4.1)|
|A-scan rate with 3D view| ~ 27 kHz (~ 0.2 volumes/s)|
|A-scan rate without 3D view| ~ 116 kHz (~ 0.89 volumes/s)|

Office Computer, Lab Computer: <br>
The performance was measured with the full processing pipeline of OCTproZ v1.0.0. The same performance is expected with OCTproZ v1.2.0 if live sinusoidal scan distortion correction is disabled.

Gaming Computer:<br>
The performance was measured with OCTproZ v1.2.0 with disabled live sinusoidal scan distortion correction.

Here are the relevant parameters that were used with Virtual OCT System and OCTproZ to determine the performance:

| |**Office Computer**|**Lab Computer**|**Gaming Computer**|**Jetson Nano**|
|:-----|:-----|:-----|:-----|:-----|
|**Virtual OCT System Settings**| | | |
|bit depth [bits]|12|12|12|12|
|Samples per raw A-scan|1024|1024|1024|1024|
|A-scan per B-scan|512|512|512|512|
|B-scans per buffer|32|256|256|32|
|Buffers per volume|8|1|1|8|
|Buffers to read from file|16|2|2|16|
|Wait after file read [us]|100|100|100|100|
|**OCTproZ Settings**| | | | |
|Bit shift sample values by 4|enabled|enabled|enabled|enabled|
|Flip every second B-scan|enabled|enabled|enabled|enabled|
|k-linearization|enabled|enabled|enabled|enabled|
|Dispersion Compensation|enabled|enabled|enabled|enabled|
|Windowing|enabled|enabled|enabled|enabled|
|Fixed-Pattern Noise Removal|enabled|enabled|enabled|enabled|
|B-scans for noise determination:|1|26|1|1|
|&emsp;once at start of measurement|enabled|enabled|enabled|enabled|
|&emsp;continuously|disabled|disabled|disabled|disabled|
|Sinusoidal scan correction|disabled|disabled|disabled|disabled|
|Log scaling|enabled|enabled|enabled|enabled|
|Stream Processed Data to Ram|enabled|disabled|disabled|enabled|



How to Determine Performance
--------
 OCTproZ provides live performance information within the sidebar in the "Processing"-tab. Live performance estimation is performed and updated every 5 seconds:
<p align="center">
  <img src="images/performance_sidebar.png" >
</p>

 It is also possible to use the [NVIDIA Visual Profiler](https://developer.nvidia.com/nvidia-visual-profiler) to analyze performance in more detail.

 For example, the following screenshot from the NVIDIA Visual Profiler shows the performance analysis of the measurement (without 3D live view) from the table at the beginning of this document with the lab computer:

 <p align="center">
  <img src="images/visualprofilerLabPC.png" >
</p>

The individual kernels are marked alphanumerically: <br>
&emsp;a) data conversion <br>
&emsp;b) kernel that combines k-linearization, windowing and dispersion compensation<br>
&emsp;c) IFFT<br>
&emsp;d) subtraction step of fixed pattern noise removal<br>
&emsp;e) truncate and logarithm<br>
&emsp;f) backward scan correction<br>
&emsp;g) copy B-scan frame to display buffer<br> 
&emsp;h) copy en face view to display buffer<br>


Additional Information
--------
- Processing happens in batches. One batch is equal to one buffer and the size of the buffer has impact on processing performance. If it is too small the processing may be slower than possible. If it is too large the application may crash as a larger buffer size results in higher GPU memory usage, which can exceed the available memory on the used GPU 
- The optimal buffer size for a specific GPU needs to be determined experimentally 
- In Virtual OCT System the buffer size can be changed by changing _bit depth_, _Samples per raw A-scan_, _A-scans per B-scan_ and _B-scans per buffer_.
- Buffer size in bytes = 	ceil(bitDepth/8) * SamplesPerRawAscan * AscansPerBscan * BscansPerBuffer
- When _B-scans per buffer_ is changed in Virtual OCT System, you should also change _Buffers per Volume_ and _Buffers to read from file_ accordingly 
- If OCTproZ crashes after setting the parameters in Virtual OCT System and starting the processing, try reducing the buffer size (for example instead of _B-scans per buffer_: 256, _Buffers per volume_: 1, _Buffers to read from file_: 2, you could try: _B-scans per buffer_: 128, _Buffers per volume_: 2, _Buffers to read from file_: 4)
- In Virtual OCT System a value greater than 2 for _Buffers to read from file_ will result in a slower processing rate displayed by OCTproZ. The reason for that is that Virtual OCT System takes more time to provide the raw data if more than two buffers should be read from a file. The processing itself is not slowed down just the time between two batches is increased. 

For performance measurement, you can use the provided [test data set](https://figshare.com/articles/SSOCT_test_dataset_for_OCTproZ/12356705). To replicate the measurements from above you need to set the value for _Samples per raw A-scan_ to 1024. This will cause the resulting OCT images to look distorted as the test data set was recorded with 1664 samples per raw A-scan. This is expected behavior that does not invalidate the performance measurement.


Performance and buffer size
--------
The following bar graph shows the A-scan rate for different buffer sizes. The __Gaming Computer__ setup without 3D live view described above was used.
To change the buffer size _Buffers to read from file_ was kept at a value of 2 and only _B-scans per buffer_ and  _Buffers per volume_ were changed. 

 <p align="center">
  <img src="images/gamingPCbufferSizePerformance.png" >
</p>---
title: 'Open Source Optical Coherence Tomography Software'
tags:
  - C++
  - CUDA
  - Optical Coherence Tomography
authors:
  - name: Miroslav Zabic
    orcid: 0000-0002-4494-6127
    affiliation: "1, 2"
  - name: Ben Matthias
    affiliation: 2 
  - name: Alexander Heisterkamp
    affiliation: 1
  - name: Tammo Ripken
    affiliation: 2        
affiliations:
 - name: Institute of Quantum Optics, Leibniz University Hannover, Welfengarten 1, 30167 Hannover, Germany
   index: 1
 - name: Industrial and Biomedical Optics Department, Laser Zentrum Hannover e.V., Hollerithallee 8, 30419 Hannover, Germany
   index: 2
date: 22 May 2020
bibliography: paper.bib
---

# Summary

Optical coherence tomography (OCT) is a non-invasive imaging technique that is often described as the optical equivalent to ultrasound imaging. 
The basic building block of OCT acquisitions is an optical interference pattern that can be processed into a depth profile, which is also called A-scan. Several adjacent A-scans can be merged into a cross-sectional image. Most research that incorporates OCT requires a software solution for processing of the acquired raw data.

Here we present an open source software package for OCT processing with an easy to use graphical user interface. The implemented OCT processing pipeline enables A-scan processing rates in the MHz range. Custom OCT systems, or any other source of Fourier Domain OCT raw data, can be integrated via a developed plug-in system, which also allows the development of custom post processing modules.

# 1. Introduction

Optical coherence tomography (OCT) is a non-invasive imaging technique used primarily in the medical field, especially in ophthalmology. The core element of any OCT system is an optical interferometer that generates a spectral fringe pattern by combining a reference beam and the backscattered light from a sample. To obtain an interpretable image from this acquired raw OCT signal several processing steps are necessary, whereby the inverse Fourier transform represents an essential step. As the possible acquisition speed for raw OCT data has increased constantly, more sophisticated methods were needed for processing and live visualization of the acquired OCT data. A particularly impressive setup was presented by Choi et al. [@choi2012spectral] that utilizes twenty FPGA-modules for real-time OCT signal processing and a graphics processing unit (GPU) for volume rendering. Nowadays, processing is typically done on graphics cards [@zhang2010real; @rasakanthan2011processing; @sylwestrzak2012four; @jian2013graphics; @wieser2014high], not FPGAs, because implementing algorithms on GPUs is more flexible and takes less time [@li2011scalable]. Most of the publications that describe OCT GPU processing do not provide the actual software implementation.
A commendable exemption is the GPU accelerated OCT processing pipeline published by Jian et al. The associated source code, which demonstrates an implementation of OCT data processing and visualization and does not include any advanced features such as a graphical user interface (GUI), already consists of several thousand lines. Thus, the most time consuming task of Fourier Domain OCT (FD-OCT) system development is not the optical setup, but the software development. The software can be separated into hardware control and signal processing, whereby the former being a highly individual, hardware-dependent software module and the latter being a generic software module, which is almost identical for many systems. To drastically reduce OCT system development time, we present OCTproZ, an open source OCT processing software package that can easily be extended, via a plug-in system, for many different hardware setups. In this paper we give a brief overview of the key functionality and structure of the software.

# 2. Basic overview of OCTproZ

OCTproZ performs live signal processing and visualization of OCT data. It is written in C++, uses the cross-platform application framework Qt [@qt] for the GUI and utilizes Nvidia’s computer unified device architecture (CUDA) [@cuda] for GPU parallel computing. A screenshot of the application can be seen in Figure \ref{fig:screenshot}.

 ![Screenshot of OCTproZ v1.0.0 Processing settings visible in the left panel can be changed before processing is started or while processing is in progress. Processed data is live visualized in 2D as cross sectional images (B-scan and en face view) and in 3D as interactive volume rendering. The live view shows a piece of wood with a couple layers of tape and a laser burned hole. \label{fig:screenshot}](figures/20191122_screenshot3d.png)


The software can be separated into three parts: main application, development kit (DevKit) and plug-ins. The main application, OCTproZ itself, contains the logic for the GUI, processing and visualization. The DevKit, which is implemented as static library, provides the necessary interface for plug-in development. Plug-ins can be one of two kinds: “Acquisition Systems” or “Extensions”. The former represent OCT systems and provide raw data to OCTproZ, the later are software modules that can extend the functionality of an OCT system (e.g., software control of a liquid lens) or provide additional custom defined post processing steps. All plug-ins are dynamic libraries that can be loaded into OCTproZ during runtime. 

# 3. Processing Pipeline

Raw data, i.e. acquired spectral fringe pattern, from the OCT system is transferred to RAM until enough data for a user-defined amount of cross-sectional images, so-called B-scans, is acquired. Via direct memory access (DMA) this raw data batch is then copied asynchronously to GPU memory where OCT signal processing is executed. If the processed data needs to be stored or post processing steps are desired the processed OCT data can be transferred back to RAM with the use of DMA. An overview of the processing steps is depicted in Figure \ref{fig:processing}.

 ![Processing pipeline of OCTproZ v1.2.0. Each box inside "OCTproZ GPU Processing" represents a CUDA kernel. Some processing steps are combinend into a single kernel (e.g. k-linearization, dispersion compensation and windowing) to enhance processing performance. \label{fig:processing}](figures/processing_pipeline_linear_v2.png) 

The processing pipeline consists of data conversion, k-linearization, numerical dispersion compensation, windowing, fast Fourier transform, fixed-pattern noise removal, truncate, logarithm, backward scan correction, sinusoidal scan correction and visualization. A detailed description of each processing step can be found in the user manual. Here we just want to mention that the implementation of the fixed-pattern noise removal is based on a publication by Moon et al. [@moon2010reference] and the volume viewer is based on source code from an open source raycaster [@pilia2018raycaster].
In order to avoid unnecessary data transfer to host memory, CUDA-OpenGL interoperability is used which allows the processed data to remain in GPU memory for visualization. 

# 4. Processing Performance 
Processing rate highly depends on the size of the raw data, the used computer hardware and resource usage by background or system processes. With common modern computer systems and typical data dimensions for OCT, OCTproZ achieves A-scan rates in the MHz range. Exemplary, Table 1 shows two computer systems and their respective processing rates for the full processing pipeline. However, since the 3D live view is computationally intensive the processing rate changes noticeably depending on whether the volume viewer is activated or not. The used raw data set consists of samples with a bit depth of 12, 1024 samples per raw A-scan, 512 A-scans per B-scan and 256 B-scans per volume. As the volume is processed in batches, the batch size was set for each system to a reasonable number of B-scans per buffer to avoid GPU memory overflow. It should be noted that this performance evaluation was done with OCTproZ v1.0.0 but is also valid for v1.2.0 if the newly introduced processing step for sinusoidal scan distortion correction is disabled.

<p style="text-align: center;"><small><b>Table 1</b>: Comparison of two computer systems and their respective processing rates for raw data sets with 12 bit per sample, 1024 samples per raw A-scan, 512 A-scans per B-scan and 256 B-scans per volume.</small></p>

| |**Office Computer**|**Lab Computer**
:-----|:-----|:-----
CPU|Intel® Core i5-7500|AMD Ryzen Threadripper 1900X
RAM|16 GB|32 GB
GPU|NVIDIA Quadro K620|NVIDIA GeForce GTX 1080 Ti
Operating system|Windows 10|Ubuntu 16.04
B-scans per buffer|32|256
With 3D live view:| | 
  A-scans per second|~**$250 \cdot 10^{3}$**|~**$4.0 \cdot 10^{6}$**
  Volumes per second|~**$1.9$**|~**$30$**
Without 3D live view:| | 
   A-scans per second|~**$300 \cdot 10^{3}$**|~**$4.8 \cdot 10^{6}$**
   Volumes per second|~**$2.2$**|~**$36$**


# 5. Plug-in System

To develop custom plug-ins for OCTproZ and thus extend its functionality, a development kit is provided. It consists of a static library and a collection of C++ header files that specify which classes and methods have to be implemented to create custom plug-ins. Currently two kinds of plug-ins exist: Acquisition Systems and Extensions. For both we made examples including source code publicly available which may be used together with the open source and cross-platform integrated development environment Qt Creator as starting point for custom plug-in development.

The main task of an Acquisition System is to provide raw data to OCTproZ. In most cases, this means that the implementation of an Acquisition System contains the software control of a data acquisition unit.

Extensions have a wide area of use cases. As they are able to receive raw data and processed data via the Qt signals and slots mechanism, they are suitable for custom post-processing routines. The exact implementation of an Extension is mainly up to the developer and can also include hardware control. Therefore, Extensions are ideal for hardware control algorithms that rely on feedback from live OCT images. The best example of this is wavefront sensorless adaptive optics with a wavefront modulator such as a deformable mirror. Particular care must be taken if high speed OCT systems are used, as the acquisition speed of OCT data may exceed the processing speed of the custom Extension. In this case, a routine within the Extension should be implemented that discards incoming OCT data if previous data is still being processed. 

# 6. Conclusion
In this paper, we introduced OCTproZ, an open source software package for live OCT signal processing. With the presented plug-in system, it is possible to develop software modules to use OCTproZ with custom OCT systems, thus reducing the OCT system development time significantly. OCTproZ is meant to be a collaborative project, where everyone involved in the field of OCT is invited to improve the software and share the changes within the community. 
We especially hope for more open source publications within the OCT community to reduce the time necessary for the replication of OCT processing algorithms and thereby accelerate scientific progress.

# Funding
This work was partially funded by the European Regional Development Fund and the state of Lower Saxony as part of the project OPhonLas (EFRE-SER 2014-2020, 85007492). 

![](figures/efre.png) 


# References

---
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

**Used hardware**
 - OS: [e.g. Win 10 64bit]
 - GPU: [e.g. Nvidia GTX 1060]

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
