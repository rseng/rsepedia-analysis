<p align="center">
  <img src="logo.png" alt=""/>
</p>

### A python based [QGIS](https://qgis.org/en/site/index.html) plugin 
[![status](https://joss.theoj.org/papers/aba2f441ab3c99e7694c97345e1255a0/status.svg)](https://joss.theoj.org/papers/aba2f441ab3c99e7694c97345e1255a0)
[![Documentation Status](https://readthedocs.org/projects/sar-tools/badge/?version=latest)](https://sar-tools.readthedocs.io/en/latest/?badge=latest)
[![License: GPL 3.0](https://img.shields.io/badge/License-GPL_3.0-green.svg)](https://opensource.org/licenses/gpl-license)
[![Open Source Love svg1](https://badges.frapsoft.com/os/v1/open-source.svg?v=103)](https://github.com/ellerbrock/open-source-badges/)
[![DOI](https://zenodo.org/badge/238603440.svg)](https://zenodo.org/badge/latestdoi/238603440)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-ffd040.svg)](https://www.python.org/)
[![GitHub release](https://img.shields.io/github/release/Narayana-Rao/PolSAR-tools.svg)](https://github.com/Narayana-Rao/PolSAR-tools/releases)
[![GitHub commits](https://img.shields.io/github/commits-since/Narayana-Rao/PolSAR-tools/v0.7.svg)](https://GitHub.com/Narayana-Rao/PolSAR-tools/commit/)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Narayana-Rao/PolSAR-tools/graphs/commit-activity)
[![Website http://www.mrslab.in/qgisplugin/](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](http://www.mrslab.in/qgisplugin/)
<p align="center">
<a href="https://hits.seeyoufarm.com"><img src="https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https://github.com/Narayana-Rao/PolSAR-tools&count_bg=%2379C83D&title_bg=%23555555&icon=go.svg&icon_color=%2300ADD8&title=hits&edge_flat=false"/></a>
</p>


## General Information
-------------------
This plugin generates derived SAR parameters (viz. vegetation indices, polarimetric decomposition parameters) from input polarimetric matrix (C3, T3, C2, T2). The input data needs to be in [PolSARpro](https://earth.esa.int/web/polsarpro/home)/[ENVI](https://www.l3harrisgeospatial.com/Software-Technology/ENVI) format (\*.bin and \*.hdr). It requires [numpy](https://numpy.org/), [matplotlib](https://matplotlib.org/) python libraries pre-installed.

## Installation

> **__Note:__** PolSAR tools requires QGIS version >=3.0.

* The easiest way (requires internet connection) : 
    - Open QGIS -> Plugins -> Manage and Install Plugins... -> select ```All``` tab -> search for ```PolSAR tools``` --> select and install plugin
* Alternative way (offline installation) : 
    - Go to [releases](https://github.com/Narayana-Rao/SAR-tools/releases) of this repository -> select desired version -> download the ```.zip``` file.
    - Open QGIS -> Plugins -> Manage and Install Plugins... -> ```install from ZIP``` tab --> select the downloaded zip --> install plugin (ignore warnings, if any).
 
## Up and running

After successful installation, find the plugin by opening **QGIS** --> Plugins --> ``PolSAR tools`` --> Process. As shown in the following figure.

<p align="center">
  <img src="help/source/files/figures/open_ui.png" alt="Opening the plugin"/>
  <p align="center"> <em>Opening the plugin</em> </p>
</p>

<p align="center">
  <img height=80% width=80% src="help/source/files/figures/main_ui.png" alt="GUI-Main window layout"/>
  <p align="center"> <em>GUI-Main window layout</em> </p>
</p>



**Layout**:

1.  Data type tabs: Functions are arranged according to the data type (full-, compact- and dual-pol).
2.  Function details viewer: Contains a list of functions for respective data tab. 
3. Derived parameter selection, required input variables and constraints.
4. Input data folder
5. Logger: displays the log of processing parameters
6. progressbar: shows the progress of the current task.
7. Credits and quick help.


Additional ``reset`` button to clear the environment, ``view data`` button to import the data into **QGIS** environment and ``Process`` button to start processing after selecting valid input data variables. 


## Available functionalities:
-----------------------------
  * Full-pol :
    * Model free 4-Component decomposition for full-pol data (MF4CF)[[11]](#11)
    * Model free 3-Component decomposition for full-pol data (MF3CF)[[4]](#4)
	* Radar Vegetation Index (RVI) [[8]](#8) 
    * Generalized volume Radar Vegetation Index (GRVI) [[2]](#2)
    * Polarimetric Radar Vegetation Index (PRVI) [[1]](#1)
    * Degree of Polarization (DOP) [[10]](#10) 

  * Compact-pol : 
    * Model free 3-Component decomposition for compact-pol data (MF3CC) [[4]](#4)
    * Improved S-Omega decomposition for compact-pol data (iS-Omega) [[7]](#7)
    * Compact-pol Radar Vegetation Index (CpRVI)  [[6]](#6)
    * Degree of Polarization (DOP)  [[10]](#10) 

  * Dual-pol:
	* Dual-pol Radar Vegetation Index (DpRVI) [[5]](#5)
	* Radar Vegetation Index (RVI) [[9]](#9)
    * Degree of Polarization (DOP) [[10]](#10) 
    * Polarimetric Radar Vegetation Index (PRVI) [[1]](#1)

## Example usage
> Note: All the following processing steps should be done in sequential manner. Sample data for all the polarization modes is provided in [sample_data](/sample_data/) folder.


**STEP 1**: Open the plugin as explained in [Up and Running section](#Up-and-running).

**STEP 2**: Select the polarimetric data type (Full/compact/dual).

<p align="center">
  <img src="help/source/files/figures/step2.png" alt="Opening the plugin" height=50% width=50%/>
  <p align="center"> <em>Selecting the polarimetric mode</em> </p>
</p>

**STEP 3**: Select the parameter/descriptor from the dropdown menu.

<p align="center">
  <img src="help/source/files/figures/step3.png" alt="Opening the plugin" height=50% width=50%/>
  <p align="center"> <em>Selecting the polarimetric descriptor</em> </p>
</p>

**STEP 4**: Provide the required input variables.
<p align="center">
  <img src="help/source/files/figures/step4.png" alt="Opening the plugin" height=50% width=50%/>
  <p align="center"> <em>Selecting the input variables</em> </p>
</p>

**STEP 5**: Select the input matrix folder.

<p align="center">
  <img src="help/source/files/figures/step5.png" alt="Opening the plugin" height=90% width=90%/>
  <p align="center"> <em>Selecting the input folder</em> </p>
</p>

**STEP 6**: Wait for the logger to prompt ```->> Ready to process.``` --> click process
> **__Note:__** Do not click process button more than once while it is processing. It may crash the QGIS and the plugin.
It is possible that the plugin may show not responding for larger datasets but please wait for the process to complete.

<p align="center">
  <img src="help/source/files/figures/step6.png" alt="Opening the plugin" height=90% width=90%/>
  <p align="center"> <em>Processing the data for selected descriptor</em> </p>
</p>

**STEP 7** (optional): Click view data to import the data into QGIS for vizualisation of the generated descriptors.

<p align="center">
  <img src="help/source/files/figures/step7a.png" alt="Opening the plugin" height=90% width=90%/>
  <p align="center"> <em>Importing the data into QGIS for visualization</em> </p>
</p>
<p align="center">
  <img src="help/source/files/figures/step7b.png" alt="Opening the plugin" height=90% width=90%/>
  <p align="center"> <em>Imported data in QGIS</em> </p>
</p>

## Functions description

Description and the details of all the core functions of this plugin are available here: [Functions_description](help/Functions_description.md)

## Contributions
1) Contribute to the software

    [Contribution guidelines for this project](help/CONTRIBUTING.md)


2) Report issues or problems with the software
	
	Please raise your issues here : <https://github.com/Narayana-Rao/SAR-tools/issues>

3) Seek support

	Please write to us: <bnarayanarao@iitb.ac.in> 

## References
-------------
<a id="1">[1]</a> 
Chang, J.G., Shoshany, M. and Oh, Y., 2018. Polarimetric Radar Vegetation Index for Biomass Estimation in Desert Fringe Ecosystems. IEEE Transactions on Geoscience and Remote Sensing, 56(12), pp.7102-7108.

<a id="2">[2]</a> 
Ratha, D., Mandal, D., Kumar, V., McNairn, H., Bhattacharya, A. and Frery, A.C., 2019. A generalized volume scattering model-based vegetation index from polarimetric SAR data. IEEE Geoscience and Remote Sensing Letters, 16(11), pp.1791-1795.

<a id="3">[3]</a> 
Mandal, D., Kumar, V., Ratha, D., J. M. Lopez-Sanchez, A. Bhattacharya, H. McNairn, Y. S. Rao, and K. V. Ramana, 2020. Assessment of rice growth conditions in a semi-arid region of India using the Generalized Radar Vegetation Index derived from RADARSAT-2 polarimetric SAR data, Remote Sensing of Environment, 237: 111561.

<a id="4">[4]</a> 
Dey, S., Bhattacharya, A., Ratha, D., Mandal, D. and Frery, A.C., 2020. Target Characterization and Scattering Power Decomposition for Full and Compact Polarimetric SAR Data. IEEE Transactions on Geoscience and Remote Sensing.

<a id="5">[5]</a> 
Mandal, D., Kumar, V., Ratha, D., Dey, S., Bhattacharya, A., Lopez-Sanchez, J.M., McNairn, H. and Rao, Y.S., 2020. Dual polarimetric radar vegetation index for crop growth monitoring using sentinel-1 SAR data. Remote Sensing of Environment, 247, p.111954.

<a id="6">[6]</a> 
Mandal, D., Ratha, D., Bhattacharya, A., Kumar, V., McNairn, H., Rao, Y.S. and Frery, A.C., 2020. A Radar Vegetation Index for Crop Monitoring Using Compact Polarimetric SAR Data. IEEE Transactions on Geoscience and Remote Sensing, 58 (9), pp. 6321-6335.

<a id="7">[7]</a> 
V. Kumar, D. Mandal, A. Bhattacharya, and Y. S. Rao, 2020. Crop Characterization Using an Improved Scattering Power Decomposition Technique for Compact Polarimetric SAR Data. International Journal of Applied Earth Observations and Geoinformation, 88: 102052.

<a id="8">[8]</a> 
Kim, Y. and van Zyl, J.J., 2009. A time-series approach to estimate soil moisture using polarimetric radar data. IEEE Transactions on Geoscience and Remote Sensing, 47(8), pp.2519-2527.

<a id="9">[9]</a> 
Trudel, M., Charbonneau, F. and Leconte, R., 2012. Using RADARSAT-2 polarimetric and ENVISAT-ASAR dual-polarization data for estimating soil moisture over agricultural fields. Canadian Journal of Remote Sensing, 38(4), pp.514-527.

<a id="10">[10]</a> 
Barakat, R., 1977. Degree of polarization and the principal idempotents of the coherency matrix. Optics Communications, 23(2), pp.147-150.

<a id="11">[11]</a> S. Dey, A. Bhattacharya, A. C. Frery, C. Lopez-Martinez and Y. S. Rao, "A Model-free Four Component Scattering Power Decomposition for Polarimetric SAR Data," in IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 2021. doi: [10.1109/JSTARS.2021.3069299](https://doi.org/10.1109/JSTARS.2021.3069299). 
---
title: 'PolSAR tools: A QGIS plugin for generating SAR descriptors'
tags:
  - SAR
  - QGIS
  - Vegetation indices
  - Polarimetric decompositions
authors:
  - name: Narayanarao Bhogapurapu
    orcid: 0000-0002-6496-7283
    affiliation: 1
  - name: Subhadip Dey
    orcid: 0000-0002-4979-0192
    affiliation: 1
  - name: Dipankar Mandal
    orcid: 0000-0001-8407-7125
    affiliation: 1
  - name: Avik Bhattacharya
    orcid: 0000-0001-6720-6108
    affiliation: 1
  - name: Y. S. Rao
    orcid: 0000-0002-6351-2391
    affiliation: 1
affiliations:
 - name: Microwave Remote Sensing Lab, Centre of Studies in Resources Engineering, Indian Institute of Technology Bombay, Mumbai-400076, India
   index: 1

date: 2 December 2020
bibliography: paper.bib

---

# Statement of need

The demand for processing tools increases with the increasing number of ***Synthetic Aperture Radar (SAR)*** satellite missions and datasets. However, to process SAR data, a minimal number of free tools are available ([PolSARpro](https://earth.esa.int/web/polsarpro/home), [SNAP](https://step.esa.int/main/toolboxes/snap/)) that consolidate all necessary pre-processing steps. Bearing this in mind, there is a need to develop specific tools for the remote sensing user community to derive polarimetric descriptors like vegetation indices and decomposition parameters. Besides, to the best of our knowledge, there are no such free tools available on a GIS platform, which is often necessary for the processing of SAR remote sensing datasets. 

Hence we have developed a plugin for ```QGIS``` that supports data for the three available polarimetric modes (i.e., full-, compact, and dual). The ```PolSAR tools``` plugin generates polarimetric descriptors (viz., vegetation indices, polarimetric decomposition parameters) from the 3x3 (C3/T3) or the 2x2 (C2/T2) covariance (coherency) matrices obtained from the ESA's [PolSARpro](https://earth.esa.int/web/polsarpro/home) software. The input data needs to be in PolSARpro format (```*.bin``` and ```*.hdr```). The plugin is coded in Python and is dependant on the QGIS framework. It uses the following Python libraries (bundled with QGIS): [numpy](https://numpy.org/), and [gdal](https://gdal.org/).

# Background
Conventional model-based decomposition methods utilize diverse scattering models and typical hierarchical rule to enumerate power components leading to numerous limitations. The ***polarimetric decomposition techniques*** incorporated in this QGIS plugin are model-free, i.e., no prior scattering models are assumed to compute the powers. The proposed decomposition techniques utilize certain novel roll-invariant target characterization parameters to decompose the total power into even bounce, odd bounce, and diffused power components. It is guaranteed that the proposed technique's powers are non-negative, which is seldom true with the existing methodologies.

In terms of target descriptors, we often use ***vegetation indices*** as plant growth proxies. While appreciating the potential of vegetation indices derived from optical remote sensing sensors, regional to global products have been supported for operational uses. These days, the Earth Observation (EO) community relies on SAR imaging technology due to its all-weather imaging capability among its numerous advantages. The SAR images are presently processed by several downstream users and are more frequently interpreted by non-radar specialists. This paradigm shift allows the utility of radar-derived vegetation indices towards a quintessential goal of Analysis Ready Data (ARD) products.

Recently, we proposed three vegetation indices: ```GRVI``` (Generalized Radar Vegetation Index) [@ratha2019generalized], ```CpRVI``` (Compact-pol Radar Vegetation Index) [@mandal2020radar], and Dual-pol Radar Vegetation Index (```DpRVI```) [@mandal2020dual] for distinct acquisition modes. These vegetation indices have provided a better opportunity to estimate biophysical parameters with fitted models directly. The retrieval of biophysical parameters from SAR observations is of vital importance for in-season monitoring of crop growth.

# PolSAR tools Audience
```PolSAR tools``` are intended for students, researchers, and polarimetry experts to derive different SAR descriptors, utilizing the QGIS and Python ecosystem of diverse tools. Especially for non-domain and application users, the plugin interface provides an easy way to process the pre-processed ***polarimetric SAR data***. 


# PolSAR tools Functionality

The functionalities of the ```PolSAR tools``` are organized into three modules to handle the data from three different SAR polarization modes. The following is the list of the available functions in the ```PolSAR tools```:

* Full-pol : 
    * Radar Vegetation Index (RVI) [@Kim_2009]
    * Generalized volume Radar Vegetation Index (GRVI) [@ratha2019generalized]
    * Polarimetric Radar Vegetation Index (PRVI) [@chang2018polarimetric] 
    * Model free 3-Component decomposition for full-pol data (MF3CF) [@dey2020target]
    * Model free 4-Component decomposition for full-pol data (MF4CF) [@dey2021mf4cf]
    * Degree of Polarization (DOP) [@barakat1977degree]
* Compact-pol :
    * Model free 3-Component decomposition for compact-pol data (MF3CC) [@dey2020target]
    * Improved S-Omega decomposition for compact-pol data (iS-Omega) [@kumar2020crop]
    * Compact-pol Radar Vegetation Index (CpRVI) [@mandal2020radar]
    * Degree of Polarization (DOP) 
 * Dual-pol :
    * Radar Vegetation Index (RVI) 
    * Dual-pol Radar Vegetation Index (DpRVI) [@mandal2020dual], 
    * Polarimetric Radar Vegetation Index (PRVI) 
    * Degree of Polarization (DOP) [@barakat1977degree]

# Acknowledgements
The authors would like to thank the developers of [QGIS Plugin Builder](https://github.com/g-sherman/Qgis-Plugin-Builder). Authors acknowledge the [GEO-AWS Earth Observation Cloud Credits Program](https://www.earthobservations.org/aws.php), which supported the computation, development, and testing of ```PolSAR tools``` on AWS cloud platform through the project: *AWS4AgriSAR-Crop inventory mapping from SAR data on cloud computing platform*.
	
# References
## Functions description

**Full-pol functions**
----------------------
Full-pol functionalities require the SAR data in the form of covariance (C3) or coherency matrix (T3). A typical file structures of T3 and C3 matrices are as follows:

<center>

<table>
<thead>
  <tr>
    <th colspan="2">C3 matrix files</th>
    <th colspan="2">T3 matrix files</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>C11.bin</td>
    <td>C11.hdr</td>
    <td>T11.bin</td>
    <td>T11.hdr</td>
  </tr>
  <tr>
    <td>C12_real.bin</td>
    <td>C12_real.hdr</td>
    <td>T12_real.bin</td>
    <td>T12_real.hdr</td>
  </tr>
  <tr>
    <td>C12_imag.bin</td>
    <td>C12_imag.hdr</td>
    <td>T12_imag.bin</td>
    <td>T12_imag.hdr</td>
  </tr>
  <tr>
    <td>C13_real.bin</td>
    <td>C13_real.hdr</td>
    <td>T13_real.bin</td>
    <td>T13_real.hdr</td>
  </tr>
  <tr>
    <td>C13_imag.bin</td>
    <td>C13_imag.hdr</td>
    <td>T13_imag.bin</td>
    <td>T13_imag.hdr</td>
  </tr>
  <tr>
    <td>C22.bin</td>
    <td>C22.hdr</td>
    <td>T22.bin</td>
    <td>T22.hdr</td>
  </tr>
  <tr>
    <td>C23_real.bin</td>
    <td>C23_real.hdr</td>
    <td>T23_real.bin</td>
    <td>T23_real.hdr</td>
  </tr>
  <tr>
    <td>C23_imag.bin</td>
    <td>C23_imag.hdr</td>
    <td>T23_imag.bin</td>
    <td>T23_imag.hdr</td>
  </tr>
  <tr>
    <td>C33.bin</td>
    <td>C33.hdr</td>
    <td>T33.bin</td>
    <td>T33.hdr</td>
  </tr>
</tbody>
</table>

</center>
<br>

 * ```GRVI``` (Generalized volume based Radar Vegetation Index): This functionality computes the generalized volume based radar vegetation index for full polarimetric SAR data. The required input and the computed output are as follows:

    ````python
        input : input_T3/C3_folder, window_size
        output: GRVI.bin
    ````
    
    The formulation of GRVI is as follows:

    <center>

    ![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\centering\text{GRVI}=\left(1-\text{GD}_{\text{GV}}\right)\Big(\frac{p}{q}\Big)^{2\,\text{GD}_{\text{GV}}},\quad0\le\text{GRVI}\le1)

    </center> 

    where, GD<sub>GV</sub> is the geodesic distance between Kennaugh matrices (K) of the observed and the generalized volume scattering model, p,q are minimum and maximum value of distances between K matrices of the observed and elementary targets respectively. A detailed explanation of GRVI is available in [[2]](#2).


 * ```MF3CF``` (Model Free 3-Component decomposition for Full-pol data): This functionality computes the model free 3 component scattering power decomposition for full polarimetric SAR data. The required input and the computed output are as follows:

    ````python
        input : input_T3/C3_folder, window_size
        output: Ps_FP.bin, Pd_FP.bin, Pv_FP.bin, Theta_FP.bin
    ````
    
    The formulation of the scattering powers (P<sub>s</sub>: Surface, P<sub>d</sub>: Double bounce, P<sub>v</sub>: volume) is as follows:
    <center>

    ![tfp](https://latex.codecogs.com/svg.latex?\Large&space;\noindent\\\P_{d}^{\text{FP}}=\frac{m_{\text{FP}}{\text{Span}}}{2}{\left(1-\sin2\theta_{\text{FP}}\right)}\\\P_{v}^{\text{FP}}={\text{Span}}\left(1-m_{\text{FP}}\right)\\\P_{s}^{\text{FP}}=\frac{m_{\text{FP}}{\text{Span}}}{2}\left(1+\sin2\theta_{\text{FP}}\right))
    
    </center>

    where m<sub>FP</sub> is degree of polarization, &theta;<sub>FP</sub> scattering type parameter, Span is the sum of the diagonal elements os coherence matrix (T3).  The derivation of these parameters in-terms of coherancey matrix (T3) elements is as shown below. Further details can be obtained from [[4]](#4)

    <center>

    ![mfp](https://latex.codecogs.com/svg.latex?\Large&space;m_{\text{FP}}=\sqrt{1-\frac{27|\mathbf{T3}|}{\big(\mathrm{Trace}(\mathbf{T3})\big)^3}};\qquad{}\tan\theta_{\text{FP}}=\frac{m_{\text{FP}}{\text{Span}}\left(T_{11}-T_{22}-T_{33}\right)}{T_{11}\left(T_{22}+T_{33}\right)+m_{\text{FP}}^{2}{\text{Span}}^{2}})
    
    </center> 

    <center>

    ![spanfp](https://latex.codecogs.com/svg.latex?\Large&space;\text{Span}=T_{11}+T_{22}+T_{33})
    
    </center>


 * ```RVI``` (Radar Vegetation Index): This functionality computes the Radar vegetation index for full polarimetric SAR data. The required input and the computed output are as follows:
    ````python
        input : input_T3/C3_folder, window_size
        output: RVI_FP.bin
    ````
    
    The formulation of RVI is as follows:

    <center>

    ![rvifp](https://latex.codecogs.com/svg.latex?\Large&space;\text{RVI}_{fp}=\frac{4\lambda_1}{\lambda_1+\lambda_2+\lambda_3})
        
    </center>

    where, &lambda;<sub>1</sub>, &lambda;<sub>2</sub> and &lambda;<sub>3</sub> are the eigen values of coherency matrix (T3) in descending order (&lambda;<sub>1</sub>> &lambda;<sub>2</sub>>&lambda;<sub>3</sub>). Further details can be found in [[8]](#8)

<!-- <center>

![rvifp](https://latex.codecogs.com/svg.latex?\Large&space;\text{RVI}_{fp}=\frac{8\sigma^\circ_{\text{HV}}}{\sigma^\circ_{\text{HH}}+\sigma^\circ_{\text{VV}}+2\sigma^\circ_{\text{HV}}})
    
</center>  -->


 * ```PRVI``` (Polarimetric Radar Vegetation Index) : This functionality computes the polarimetric Radar vegetation index for full polarimetric SAR data. The required input and the computed output are as follows:
    
    ````python
        input : input_T3/C3_folder, window_size
        output: PRVI_FP.bin
    ````
    The formlation of PRVI interms of degree of polarization and cross-pol backscatter intensity can be expressed as follows: 
    <center>

    ![grvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{PRVI}_{fp}=(1-\text{DOP}_{fp})\sigma^\circ_{\text{XY}})
        
    </center> 

    where DOP<sub>fp</sub> 3D Barakt degree of polarization and can be expressed as shown below. Further details on the PRVI can be found in [[1]](#1)

 * ```DOP``` (Degree of Polarization):  This functionality computes the 3D Barakat degree of polarization for full polarimetric SAR data. The required input and the computed output are as follows:

    ````python
        input : input_T3/C3_folder, window_size
        output: DOP_FP.bin
    ````
    <center>

    ![dopfp](https://latex.codecogs.com/svg.latex?\Large&space;\text{DOP}_{fp}=\sqrt{1-\frac{27\times\text{det([T3])}}{\text{(Trace[T3])}^3}})
        
    </center> 

    Further details on the Barakat Degree of polarization can be found in [[10]](#10)





**Compact-pol functions**
-------------------------
Compact-pol functionalities require the SAR data in the form of 2x2 covariance matrix (C2). A typical file structures of C2 matrix is as follows:

<center>

<table>
<thead>
  <tr>
    <th colspan="2">C2 matrix files</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>C11.bin</td>
    <td>C11.hdr</td>
  </tr>
  <tr>
    <td>C12_real.bin</td>
    <td>C12_real.hdr</td>
  </tr>
  <tr>
    <td>C12_imag.bin</td>
    <td>C12_imag.hdr</td>
  </tr>
  <tr>
    <td>C22.bin</td>
    <td>C22.hdr</td>
  </tr>
</tbody>
</table>

</center>
<br>

 * ```MF3CC``` (Model Free 3-Component decomposition for Compact-pol data): This functionality computes the model free 3 component scattering power decomposition for compact polarimetric SAR data. The required input and the computed output are as follows:
  
    ````python
        input : input_C2_folder, window_size, tau
        output: Ps_CP.bin, Pd_CP.bin, Pv_CP.bin, Theta_CP.bin
    ````  

    The formulation of the scattering powers (P<sub>s</sub>: Surface, P<sub>d</sub>: Double bounce, P<sub>v</sub>: volume) is as follows:
    
    <center>

    ![cppowers](https://latex.codecogs.com/svg.latex?\Large&space;\\\P_{d}^{\text{CP}}=\frac{m_{\text{FP}}{S_0}}{2}{\left(1-\sin2\theta_{\text{CP}}\right)};\\\P_{v}^{\text{CP}}={S_0}\left(1-m_{\text{CP}}\right);\\\P_{s}^{\text{CP}}=\frac{m_{\text{CP}}{S_0}}{2}\left(1+\sin2\theta_{\text{CP}}\right))
    
    </center>    

    where m<sub>CP</sub> is degree of polarization, &theta;<sub>CP</sub>: scattering type parameter, S<sub>0</sub>,  - S<sub>3</sub>, are Stokes parameters. The derivation of these parameters in-terms of covariance matrix (C2) elements is as shown below. Further details can be obtained from [[4]](#4)

    <center>

    ![mcp](https://latex.codecogs.com/svg.latex?\Large&space;m_{\text{CP}}=\sqrt{1-\frac{4|\mathbf{C2}|}{\big(\mathrm{Trace}(\mathbf{C2})\big)^2}};\qquad{}\tan\theta_{\text{CP}}=\frac{m_{\text{CP}}{S_0}\left(\text{OC}-\text{SC}\right)}{\text{OC}\times\text{SC}+m_{\text{CP}}^{2}{S_0}^{2}})
    
    </center> 



    <center>

    ![sparm](https://latex.codecogs.com/svg.latex?\Large&space;\\\S_0=\text{C11+C22};\qquad{}S_1=\text{C11-C22};\\\S_2=\text{C12+C21};\qquad{}S_3=\pm\text{j(C12-C21)})
    
    ![cparm](https://latex.codecogs.com/svg.latex?\Large&space;\text{SC}=\frac{S_0-S_3}{2};\qquad{}\text{OC}=\frac{S_0+S_3}{2};)
    
    </center>



 * ```CpRVI``` (Compact-pol Radar Vegetation Index): This functionality computes the compact-pol radar vegetation index for compact polarimetric SAR data. The required input and the computed output are as follows:

    ````python
        input : input_C2_folder, window_size
        output: CpRVI.bin
    ````

    The formulation of the CpRVI is as follows:
    <center>
    
    ![cprvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{CpRVI}=\left(1-\dfrac{3}{2}\text{GD}_{\text{ID}}\right)\Big(\frac{p}{q}\Big)^{2(\frac{3}{2}\text{GD}_{\text{ID}})})

    ![pqcp](https://latex.codecogs.com/svg.latex?\Large&space;p=\text{min\\{SC,OC\\}},q=\text{max\\{SC,OC\\}})
    
    ![cparm](https://latex.codecogs.com/svg.latex?\Large&space;\text{SC}=\frac{S_0-S_3}{2};\qquad{}\text{OC}=\frac{S_0+S_3}{2};)
    

    ![sparm](https://latex.codecogs.com/svg.latex?\Large&space;\\\S_0=\text{C11+C22};\qquad{}S_1=\text{C11-C22};\\\S_2=\text{C12+C21};\qquad{}S_3=\pm\text{j(C12-C21)})


    </center> 

    where, GD<sub>ID</sub> is the geodesic distance between Kennaugh matrices (K) of the observed and the ideal depolarizer, p,q are minimum and maximum values of SC and OC which are functions of stocks parameters (S<sub>0</sub>, S<sub>1</sub>, S<sub>2</sub>, and S<sub>3</sub>). A detailed explanation of CpRVI is available in [[6]](#6).




 * ```iS-Omega``` (improved S-&Omega; decomposition): 
    This functionality computes the scattering powers for compact polarimetric SAR data. This is an improved decomposition technique based on Stokes vector(S) and the polarized power fraction (&Omega;). The required input and the computed output are as follows:
    
    ````python
        input : input_C2_folder, window_size, tau, psi, chi
        output: Ps_iSOmega.bin, Pd_iSOmega.bin,Pv_iSOmega.bin
    ````

    The stokes paramters can be written in terms of the covariance matrx (C2) elements as follows:
    
    <center>
    
    ![sparm](https://latex.codecogs.com/svg.latex?\Large&space;\\\S_0=\text{C11+C22};\qquad{}S_1=\text{C11-C22};\\\S_2=\text{C12+C21};\qquad{}S_3=\pm\text{j(C12-C21)})    

    </center> 

    Then, the parameters Same-sense Circular (SC) and Opposite-sense Circular (OC) can be expressed as follows:

    <center>

    ![cparm](https://latex.codecogs.com/svg.latex?\Large&space;\text{SC}=\frac{S_0-S_3}{2};\qquad{}\text{OC}=\frac{S_0+S_3}{2};)
    </center>    
    <!-- <center>

    ![cparm](https://latex.codecogs.com/svg.latex?\Large&space;\Vec{\mathbf{S}}=\begin{bmatrix}S_{0}\\\S_{1}\\\S_{2}\\\\S_{3}\end{bmatrix}=\begin{bmatrix}C_{11}+C_{22}\\\C_{11}-C_{22}\\\C_{12}+C_{21}\\\pm\left(C_{12}-C_{21}\right)\end{bmatrix})

    </center> -->

    Now, based on the ratio of SC and OC the decomposition powers can be derived as given below. Further details can be found in [[7]](#7)

    <br>
    <br>
    <center>

    ![cparm](https://latex.codecogs.com/svg.latex?\Large&space;\text{SC/OC}<1;\qquad{}\qquad{}\qquad{}\text{SC/OC}>1\\\P_s=\Omega\left(S_{0}-\left(1-\Omega\right)\text{SC}\right);\qquad{}P_s=\Omega\left(1-\Omega\right)\text{OC}\\\P_d=\Omega\left(1-\Omega\right)\text{SC};\qquad{}P_d=\Omega\left(S_{0}-\left(1-\Omega\right)\text{OC}\right))

    ![cparm](https://latex.codecogs.com/svg.latex?\Large&space;P_v=S_{r0}\left(1-\Omega\right))

    
    </center>    



 * ```DOP``` (Degree of Polarization): This functionality computes the degree of polarization for compact polarimetric SAR data. The required input and the computed output are as follows:
    

    ````python
        input : input_c2_folder, window_size, tau
        output: DOP_CP.bin
    ````  
    
    The conventional degree of polarization in terms of stokes paramters can be written as follows:

    <center>

    ![dopcp](https://latex.codecogs.com/svg.latex?\Large&space;\text{DOP}_{cp}=\frac{\sqrt{S^2_1+S^2_2+S^2_3}}{S_0})
        
    </center> 
    
    where, 
    <br>

    <center>
    
    ![sparm](https://latex.codecogs.com/svg.latex?\Large&space;\\\S_0=\text{C11+C22};\qquad{}S_1=\text{C11-C22};\\\S_2=\text{C12+C21};\qquad{}S_3=\pm\text{j(C12-C21)})    

    </center>  







**Dual-pol**
------------
Dual-pol functionalities require the SAR data in the form of 2x2 covariance matrix (C2). A typical file structures of C2 matrix is as follows:

<center>

<table>
<thead>
  <tr>
    <th colspan="2">C2 matrix files</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>C11.bin</td>
    <td>C11.hdr</td>
  </tr>
  <tr>
    <td>C12_real.bin</td>
    <td>C12_real.hdr</td>
  </tr>
  <tr>
    <td>C12_imag.bin</td>
    <td>C12_imag.hdr</td>
  </tr>
  <tr>
    <td>C22.bin</td>
    <td>C22.hdr</td>
  </tr>
</tbody>
</table>

</center>
<br>

 * ```DpRVI``` (Dual-pol Radar Vegetation Index): This functionality computes the dual polarimetric radar vegetation index for dual polarimetric (HH | HV), (VV | VH) SAR data. The required input and the computed output are as follows:

    ````python
        input : input_C2_folder, window_size
        output: DpRVI.bin
    ````
    
    The formulation of DpRVI is as follows:
    <br>

    <center>

    ![dprvi](https://latex.codecogs.com/svg.latex?\Large&space;\text{DpRVI}=1-\Big(\frac{\lambda_1}{\lambda_1+\lambda_2}\Big)\sqrt{1-\frac{4\times\text{det([C2])}}{\text{(Trace[C2])}^2}})
    
    </center> 

    where, C2 is co-variance matrix,  and  &lambda;<sub>1</sub> and &lambda;<sub>2</sub> are the eigen values of C2 matrix in descending order (&lambda;<sub>1</sub> > &lambda;<sub>2</sub>). Further details on DpRVI can be obtained from [[5]](#5)




 * ```RVI``` (Radar Vegetation Index): This functionality computes the radar vegetation index for dual polarimetric (HH | HV), (VV | VH) SAR data. The required input and the computed output are as follows:

    ````python
        input : input_c2_folder, window_size
        output: RVI_dp.bin
    ````
    The formulation of RVI is as follows:

    <center>

    ![rvidp](https://latex.codecogs.com/svg.latex?\Large&space;\text{RVI}_{dp}=\frac{4\sigma^\circ_{\text{XY}}}{\sigma^\circ_{\text{XX}}+\sigma^\circ_{\text{XY}}})
        
    </center> 

    where, &sigma;<sup>o</sup><sub>XX</sub> and &sigma;<sup>o</sup><sub>XY</sub> co- and cross-pol backscatter intensities.

 * ```PRVI``` (Polarimetric Radar Vegetation Index): This functionality computes the polarimetric radar vegetation index for dual polarimetric (HH | HV), (VV | VH) SAR data. The required input and the computed output are as follows:

    ````python
        input : input_c2_folder, window_size
        output: PRVI_dp.bin
    ````
    
    The formulation of PRVI is as follows: 
    <br>

    <center>

    ![prvidp](https://latex.codecogs.com/svg.latex?\Large&space;\text{PRVI}_{dp}=\Big(1-\sqrt{1-\frac{4\times\text{det([C2])}}{\text{(Trace[C2])}^2}}\Big)\sigma^\circ_{\text{XY}})
        
    </center> 




 * ```DOP``` (Degree of Polarization):  This functionality computes the 2D Barakat degree of polarization for dual polarimetric (HH | HV), (VV | VH) SAR data. The required input and the computed output are as follows:
    
    ````python
        input : input_c2_folder, window_size
        output: dop_dp.bin
    ````

    <center>

    ![dopdp](https://latex.codecogs.com/svg.latex?\Large&space;\text{DOP}_{dp}=\sqrt{1-\frac{4\times\text{det([C2])}}{\text{(Trace[C2])}^2}})
        
    </center> 

    Further details on the Barakat Degree of polarization can be found in [[10]](#10)




## References
-------------
<a id="1">[1]</a> 
Chang, J.G., Shoshany, M. and Oh, Y., 2018. Polarimetric Radar Vegetation Index for Biomass Estimation in Desert Fringe Ecosystems. IEEE Transactions on Geoscience and Remote Sensing, 56(12), pp.7102-7108.

<a id="2">[2]</a> 
Ratha, D., Mandal, D., Kumar, V., McNairn, H., Bhattacharya, A. and Frery, A.C., 2019. A generalized volume scattering model-based vegetation index from polarimetric SAR data. IEEE Geoscience and Remote Sensing Letters, 16(11), pp.1791-1795.

<a id="3">[3]</a> 
Mandal, D., Kumar, V., Ratha, D., J. M. Lopez-Sanchez, A. Bhattacharya, H. McNairn, Y. S. Rao, and K. V. Ramana, 2020. Assessment of rice growth conditions in a semi-arid region of India using the Generalized Radar Vegetation Index derived from RADARSAT-2 polarimetric SAR data, Remote Sensing of Environment, 237: 111561.

<a id="4">[4]</a> 
Dey, S., Bhattacharya, A., Ratha, D., Mandal, D. and Frery, A.C., 2020. Target Characterization and Scattering Power Decomposition for Full and Compact Polarimetric SAR Data. IEEE Transactions on Geoscience and Remote Sensing.

<a id="5">[5]</a> 
Mandal, D., Kumar, V., Ratha, D., Dey, S., Bhattacharya, A., Lopez-Sanchez, J.M., McNairn, H. and Rao, Y.S., 2020. Dual polarimetric radar vegetation index for crop growth monitoring using sentinel-1 SAR data. Remote Sensing of Environment, 247, p.111954.

<a id="6">[6]</a> 
Mandal, D., Ratha, D., Bhattacharya, A., Kumar, V., McNairn, H., Rao, Y.S. and Frery, A.C., 2020. A Radar Vegetation Index for Crop Monitoring Using Compact Polarimetric SAR Data. IEEE Transactions on Geoscience and Remote Sensing, 58 (9), pp. 6321-6335.

<a id="7">[7]</a> 
V. Kumar, D. Mandal, A. Bhattacharya, and Y. S. Rao, 2020. Crop Characterization Using an Improved Scattering Power Decomposition Technique for Compact Polarimetric SAR Data. International Journal of Applied Earth Observations and Geoinformation, 88: 102052.

<a id="8">[8]</a> 
Kim, Y. and van Zyl, J.J., 2009. A time-series approach to estimate soil moisture using polarimetric radar data. IEEE Transactions on Geoscience and Remote Sensing, 47(8), pp.2519-2527.

<a id="9">[9]</a> 
Trudel, M., Charbonneau, F. and Leconte, R., 2012. Using RADARSAT-2 polarimetric and ENVISAT-ASAR dual-polarization data for estimating soil moisture over agricultural fields. Canadian Journal of Remote Sensing, 38(4), pp.514-527.

<a id="10">[10]</a> 
Barakat, R., 1977. Degree of polarization and the principal idempotents of the coherency matrix. Optics Communications, 23(2), pp.147-150.

<a id="11">[11]</a> 
Lee, J.S. and Pottier, E., 2009. Polarimetric radar imaging: from basics to applications. CRC press.## Contributing

Contribute to the software

  * Setting up environment
    - Download and install [anaconda](https://www.anaconda.com/products/individual) (python version: >3.0)
    - Download and install [Qt Designer](https://build-system.fman.io/qt-designer-download) (Qt version: 5.0)

  * Preparing your own descriptor function:

    - All the core functions are arranged in separate modules. The generic structure of a module is as follows:
        
        ````python
            def your_function_name (data_stack, **vars):
                ...
                code
                ...
                return your_descriptor
                  
        ````
      data_stack : 3-D array of the polarimetric matrix (N x N x 9 (T3/C3); N x N x 4 (C2/T2)).

      \**vars :  list of required variables(E.g. **window_size**, **ellipticity** etc.)


> A template module is provided in the functions folder [mod_template](../functions/mod_template.py)

  * Updating the GUI 
    - Open the **mainWindow.ui** file in the Qt Designer and update the elements.
    - link the module with the ui elements in [SAR_Tools.py](../SAR_Tools.py)

> Please contribute by [forking](http://help.github.com/forking/) and sending a [pull request](https://docs.github.com/en/github/getting-started-with-github/github-glossary#pull).


.. SAR tools documentation master file, created by
   sphinx-quickstart on Mon Dec 28 22:30:50 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PolSAR tools's documentation!
=========================================
PolSAR tools is a QGIS plugin, which generates derived SAR parameters ( vegetation indices, polarimetric decomposition parameters) from input polarimetric matrix (C3, T3, C2, T2).

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   files/General_Information
   files/functions_description
   files/References



.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

About
=====
The current version of PolSAR tools is v0.6.3 and is licensed under the GPL-3.0 license. 

General Information
===================

This plugin generates derived SAR parameters (viz. vegetation indices, polarimetric decomposition parameters) from input polarimetric matrix (C3, T3, C2, T2). The input data needs to be in `PolSARpro`_/`ENVI`_ format (\*.bin and \*.hdr). It requires `numpy`_, `matplotlib`_ python libraries pre-installed.

Installation
-------------------

.. note::

    PolSAR tools requires QGIS version >=3.0.


* The easiest way (requires internet connection) : 
	- Open QGIS -> Plugins -> Manage and Install Plugins... --> select ``All`` tab -> search for ``PolSAR tools`` --> select and install plugin
* Alternative way (offline installation) : 
	- Go to `releases`_ of PolSAR tools -> select desired version -> download the ``.zip`` file.
	- Open QGIS -> Plugins -> Manage and Install Plugins... --> ``install from ZIP`` tab --> select the downloaded zip --> install plugin (ignore warnings, if any).

.. _Up and running:

Up and running
--------------------

After successful installation, find the plugin by opening **QGIS** --> Plugins --> ``PolSAR tools`` --> Process. As shown in the following figure.

.. figure:: figures/open_ui.png
    :scale: 80%
    :align: center
    
    Opening the plugin 

.. figure:: figures/main_ui.png
    :scale: 60%
    :align: center
    
    GUI-Main window layout

**Layout**:

1. Data type tabs: Functions are arranged according to the data dype (full-, compact- and dual-pol).
2. Function detials viewer: Contains list of functions for respective data tab. 
3. Derived arameter selection, required input variables and constraints.
4. Input data folder
5. Logger: displays the log of procesing parameters
6. progressbar: displays the progress of current task.
7. Credits and quick help.


Additional ``reset`` button to clear the envirinment, ``view data`` button to import the data into **QGIS** environment and ``Process`` button to start processing after selecting valid input data variables. 

 

Available functionalities
--------------------------
  1. **Full-pol** 

    * Model free 3-Component decomposition for full-pol data (MF3CF) 
    * Radar Vegetation Index (RVI) 
    * Generalized volume Radar Vegetation Index (GRVI) 
    * Polarimetric Radar Vegetation Index (PRVI) 
    * Degree of Polarization (DOP) 

  2. **Compact-pol**

    * Model free 3-Component decomposition for compact-pol data (MF3CC) 
    * Improved S-Omega decomposition for compact-pol data (iS-Omega) 
    * Compact-pol Radar Vegetation Index (CpRVI) 
    * Degree of Polarization (DOP)  

  3. **Dual-pol**

    * Dual-pol Radar Vegetation Index (`DpRVI <functions/dual_pol/DpRVI.html>`_) 
    * Radar Vegetation Index (`RVI <functions/dual_pol/RVI_dp.html>`_) 
    * Degree of Polarization (`DOP <functions/dual_pol/DOP_dp.html>`_)
    * Polarimetric Radar Vegetation Index (`PRVI <functions/dual_pol/PRVI_dp.html>`_) 

Example usage
--------------
.. note::

    All the following processing steps should be done in sequential manner. Sample data for all the polarization modes is provided in [sample_data](/sample_data/) folder.


**STEP 1**: Open the plugin as explained in :ref:`Up and running` section.

**STEP 2**: Select the polarimetric data type (Full/compact/dual).

.. figure:: figures/step2.png
    :scale: 50%
    :align: center
    
    Selecting the polarimetric mode

**STEP 3**: Select the parameter/descriptor from the dropdown menu.

.. figure:: figures/step3.png
    :scale: 50%
    :align: center
    
    Selecting the polarimetric descriptor

**STEP 4**: Provide the required input variables.

.. figure:: figures/step4.png
    :scale: 50%
    :align: center
    
    Selecting the input variables

**STEP 5**: Select the input matrix folder.

.. figure:: figures/step5.png
    :scale: 45%
    :align: center
    
    Selecting the input folder

**STEP 6**: Wait for the logger to prompt ```->> Ready to process.``` --> click process

.. note::
    Do not click process button more than once while it is processing. It may crash the QGIS and the plugin.
    It is possible that the plugin may show not responding for larger datasets but please wait for the process to complete.

.. figure:: figures/step6.png
    :scale: 45%
    :align: center
    
    Processing the data for selected descriptor

**STEP 7** (optional): Click view data to import the data into QGIS for vizualisation of the generated descriptors.

.. figure:: figures/step7a.png
    :scale: 45%
    :align: center
    
    Importing the data into QGIS for visualization

.. figure:: figures/step7b.png
    :scale: 45%
    :align: center
    
    Imported data in QGIS

Functions description
---------------------

Description and the details of all the core functions of this plugin are available here: (`Functions description <functions_description.html>`_)

Contributions
-------------

1) Contribute to the software

    `Contribution guidelines for this project  <https://github.com/Narayana-Rao/PolSAR-tools/blob/master/help/CONTRIBUTING.md>`_


2) Report issues or problems with the software
  
  Please raise your issues here : https://github.com/Narayana-Rao/PolSAR-tools/issues

3) Seek support

  Please write to us: bnarayanarao@iitb.ac.in


.. _PolSARpro: https://earth.esa.int/web/polsarpro/home
.. _ENVI: https://www.l3harrisgeospatial.com/Software-Technology/ENVI
.. _numpy: https://numpy.org/
.. _matplotlib: https://matplotlib.org/
.. _releases: https://github.com/Narayana-Rao/SAR-tools/releasesReferences
----------

References of the research work used in this plugin.



Functions Description 
======================


**Full-pol functions**
----------------------
Full-pol functionalities require the SAR data in the form of covariance (C3) or coherency matrix (T3). A typical file structures of T3 and C3 matrices are as follows: 

+-----------------------------+-----------------------------+
|       C3 matrix files       |       T3 matrix files       |
+==============+==============+==============+==============+
| C11.bin      | C11.hdr      | T11.bin      | T11.hdr      |
+--------------+--------------+--------------+--------------+
| C12_real.bin | C12_real.hdr | T12_real.bin | T12_real.hdr |
+--------------+--------------+--------------+--------------+
| C12_imag.bin | C12_imag.hdr | T12_imag.bin | T12_imag.hdr |
+--------------+--------------+--------------+--------------+
| C13_real.bin | C13_real.hdr | T13_real.bin | T13_real.hdr |
+--------------+--------------+--------------+--------------+
| C13_imag.bin | C13_imag.hdr | T13_imag.bin | T13_imag.hdr |
+--------------+--------------+--------------+--------------+
| C22.bin      | C22.hdr      | T22.bin      | T22.hdr      |
+--------------+--------------+--------------+--------------+
| C23_real.bin | C23_real.hdr | T23_real.bin | T23_real.hdr |
+--------------+--------------+--------------+--------------+
| C23_imag.bin | C23_imag.hdr | T23_imag.bin | T23_imag.hdr |
+--------------+--------------+--------------+--------------+
| C33.bin      | C33.hdr      | T33.bin      | T33.hdr      |
+--------------+--------------+--------------+--------------+

Following are the avaialble functions for full-pol data:

.. toctree::
    :maxdepth: 3

    functions/full_pol/RVI_fp
    functions/full_pol/GRVI
    functions/full_pol/MF3CF
    functions/full_pol/PRVI_fp
    functions/full_pol/DOP_fp
 


**Compact-pol functions**
-------------------------
Compact-pol functionalities require the SAR data in the form of 2x2 covariance matrix (C2). A typical file structures of C2 matrix is as follows:

+-----------------------------+
|       C2 matrix files       |
+==============+==============+
| C11.bin      | C11.hdr      |
+--------------+--------------+
| C12_real.bin | C12_real.hdr |
+--------------+--------------+
| C12_imag.bin | C12_imag.hdr |
+--------------+--------------+
| C22.bin      | C22.hdr      |
+--------------+--------------+

.. toctree::
    :maxdepth: 3

    functions/compact_pol/cprvi
    functions/compact_pol/iS_Omega
    functions/compact_pol/MF3CC
    functions/compact_pol/DOP_cp



**Dual-pol**
------------
Dual-pol functionalities require the SAR data in the form of 2x2 covariance matrix (C2). A typical file structures of C2 matrix is as follows:

+-----------------------------+
|       C2 matrix files       |
+==============+==============+
| C11.bin      | C11.hdr      |
+--------------+--------------+
| C12_real.bin | C12_real.hdr |
+--------------+--------------+
| C12_imag.bin | C12_imag.hdr |
+--------------+--------------+
| C22.bin      | C22.hdr      |
+--------------+--------------+

.. toctree::
    :maxdepth: 3

    functions/dual_pol/RVI_dp
    functions/dual_pol/DpRVI
    functions/dual_pol/PRVI_dp
    functions/dual_pol/DOP_dp
    

``MF3CF`` (Model Free 3-Component decomposition for Full-pol data)
=====================================================================
This functionality computes the model free 3 component scattering power decomposition for full polarimetric SAR data. The required input and the computed output are as follows:

.. code-block:: python
        
        input : input_T3/C3_folder, window_size
        output: Ps_FP.bin, Pd_FP.bin, Pv_FP.bin, Theta_FP.bin

The formulation of the scattering powers (:math:`P_s` : Surface, :math:`P_d` : Double bounce, :math:`P_v` : volume) is as follows:

.. math::
    
    P_{d}^{\text{FP}}=\frac{m_{\text{FP}}{\text{Span}}}{2}{\left(1-\sin2\theta_{\text{FP}}\right)}\\P_{v}^{\text{FP}}={\text{Span}}\left(1-m_{\text{FP}}\right)\\P_{s}^{\text{FP}}=\frac{m_{\text{FP}}{\text{Span}}}{2}\left(1+\sin2\theta_{\text{FP}}\right)

where :math:`m_\text{FP}` is degree of polarization, :math:`\theta_\text{FP}` scattering type parameter, Span is the sum of the diagonal elements os coherence matrix (T3).  The derivation of these parameters in-terms of coherancey matrix (T3) elements is as shown below. Further details can be obtained from [[4]](#4)

.. math::

    m_{\text{FP}}=\sqrt{1-\frac{27|\mathbf{T3}|}{\big(\mathrm{Trace}(\mathbf{T3})\big)^3}};\qquad{}\tan\theta_{\text{FP}}=\frac{m_{\text{FP}}{\text{Span}}\left(T_{11}-T_{22}-T_{33}\right)}{T_{11}\left(T_{22}+T_{33}\right)+m_{\text{FP}}^{2}{\text{Span}}^{2}}
    
    \text{Span}=T_{11}+T_{22}+T_{33}
    
``PRVI`` (Polarimetric Radar Vegetation Index)
===============================================
This functionality computes the polarimetric Radar vegetation index for full polarimetric SAR data. The required input and the computed output are as follows:

.. code-block:: python

        input : input_T3/C3_folder, window_size
        output: PRVI_FP.bin

The formlation of PRVI interms of degree of polarization and cross-pol backscatter intensity can be expressed as follows: 

.. math::

    \text{PRVI}_{fp}=(1-\text{DOP}_{fp})\sigma^\circ_{\text{XY}}


where, :math:`\text{DOP}_{fp}` 3D Barakt degree of polarization and can be expressed as shown below. Further details on the PRVI can be found in [[1]](#1)

.. math::

    \text{DOP}_{fp}=\sqrt{1-\frac{27\times\text{det([T3])}}{\text{(Trace[T3])}^3}}

``DOP`` (Degree of Polarization)
=================================
This functionality computes the 3D Barakat degree of polarization for full polarimetric SAR data. The required input and the computed output are as follows:

.. code-block:: python

    input : input_T3/C3_folder, window_size
    output: DOP_FP.bin

.. math::

    \text{DOP}_{fp}=\sqrt{1-\frac{27\times\text{det([T3])}}{\text{(Trace[T3])}^3}}

Further details on the Barakat Degree of polarization can be found in [[10]](#10)

``GRVI`` (Generalized volume based Radar Vegetation Index)
==========================================================
This functionality computes the generalized volume based radar vegetation index for full polarimetric SAR data. The required input and the computed output are as follows:

.. code-block:: python

    input : input_T3/C3_folder, window_size
    output : GRVI.bin

    
The formulation of GRVI is as follows:

.. math::

    \text{GRVI} = \left(1 - \text{GD}_{\text{GV}}\right)\Big(\frac{p}{q}\Big)^{2\,\text{GD}_{\text{GV}}}, \quad 0\le \text{GRVI} \le 1

where, :math:`\text{GD}_{\text{GV}}` is the geodesic distance between Kennaugh :math:`(\mathbf{K})` matrices of the observed and the generalized volume scattering model, :math:`p,q` are minimum and maximum value of distances between :math:`\mathbf{K}` matrices of the observed and elementary targets respectively. A detailed explanation of GRVI is available in.




``RVI`` (Radar Vegetation Index)
=================================
This functionality computes the Radar vegetation index for full polarimetric SAR data. The required input and the computed output are as follows:
    
.. code-block:: python

    input : input_T3/C3_folder, window_size
    output : RVI_FP.bin
   
The formulation of RVI is as follows:

.. math::

    \text{RVI}_{fp} = \frac{4 \times \lambda_3}{\lambda_1+\lambda_2+\lambda_3} 

where, :math:`\lambda_1, \lambda_2` and :math:`\lambda_3` are the eigen values of coherency matrix (T3) in descending order (:math:`\lambda_1 > \lambda_2 > \lambda_3`). Further details can be found in [[8]](#8)``MF3CC`` (Model Free 3-Component decomposition for Compact-pol data)
=======================================================================
 This functionality computes the model free 3 component scattering power decomposition for compact polarimetric SAR data. The required input and the computed output are as follows:

.. code-block:: python

    input : input_C2_folder, window_size, tau
    output: Ps_CP.bin, Pd_CP.bin, Pv_CP.bin, Theta_CP.bin

The formulation of the scattering powers (:math:`P_s` : Surface, :math:`P_d`: Double bounce, :math:`P_v`: volume) is as follows:

.. math::

    P_{d}^{\text{CP}}=\frac{m_{\text{FP}}{S_0}}{2}{\left(1-\sin2\theta_{\text{CP}}\right)};\\P_{v}^{\text{CP}}={S_0}\left(1-m_{\text{CP}}\right);\\P_{s}^{\text{CP}}=\frac{m_{\text{CP}}{S_0}}{2}\left(1+\sin2\theta_{\text{CP}}\right)

where :math:`m_\text{CP}` is degree of polarization; :math:`\theta_\text{CP}` : scattering type parameter; :math:`S_0, S_3`, are Stokes parameters. The derivation of these parameters in-terms of covariance matrix (C2) elements is as shown below. Further details can be obtained from [[4]](#4)

.. math::

    m_{\text{CP}}=\sqrt{1-\frac{4|\mathbf{C2}|}{\big(\mathrm{Trace}(\mathbf{C2})\big)^2}};\qquad{}\tan\theta_{\text{CP}}=\frac{m_{\text{CP}}{S_0}\left(\text{OC}-\text{SC}\right)}{\text{OC}\times\text{SC}+m_{\text{CP}}^{2}{S_0}^{2}}\\
    S_0=\text{C11+C22};\qquad{}\qquad{}\qquad{}\qquad{}\qquad{}S_1=\text{C11-C22};\\
    S_2=\text{C12+C21};\qquad{}\qquad{}\qquad{}\qquad{}S_3=\pm\text{j(C12-C21)};\\
    \text{SC}=\frac{S_0-S_3}{2};\qquad{}\qquad{}\qquad{}\qquad{}\qquad{}\text{OC}=\frac{S_0+S_3}{2};``iS-Omega`` (improved S-:math:`\Omega` decomposition)
========================================================
This functionality computes the scattering powers for compact polarimetric SAR data. This is an improved decomposition technique based on Stokes vector(S) and the polarized power fraction (:math:`\Omega`). The required input and the computed output are as follows:

.. code-block:: python

    input : input_C2_folder, window_size, tau, psi, chi
    output: Ps_iSOmega.bin, Pd_iSOmega.bin,Pv_iSOmega.bin

The stokes paramters can be written in terms of the covariance matrx (C2) elements as follows:

.. math::

    S_0=\text{C11+C22};\qquad{}S_1=\text{C11-C22};\\
    S_2=\text{C12+C21};\qquad{}S_3=\pm\text{j(C12-C21)}

Then, the parameters Same-sense Circular (:math:`\text{SC}`) and Opposite-sense Circular (:math:`\text{OC}`) can be expressed as follows:

.. math::

    \text{SC}=\frac{S_0-S_3}{2};\qquad{}\text{OC}=\frac{S_0+S_3}{2};

Now, based on the ratio of :math:`\text{SC}` and :math:`\text{OC}` the decomposition powers can be derived as given below. Further details can be found in [[7]](#7)

.. math::

    \text{SC/OC}<1;\qquad{}\qquad{}\qquad{}\qquad{}\qquad{}\qquad{}\qquad{}\text{SC/OC}>1\\P_s=\Omega\left(S_{0}-\left(1-\Omega\right)\text{SC}\right);\qquad{}\qquad{}\qquad{}P_s=\Omega\left(1-\Omega\right)\text{OC}\\P_d=\Omega\left(1-\Omega\right)\text{SC};\qquad{}\qquad{}\qquad{}P_d=\Omega\left(S_{r0}-\left(1-\Omega\right)\text{OC}\right)
    \\P_v=S_{0}\left(1-\Omega\right)\qquad{}\qquad{}\qquad{}\qquad{}\qquad{}P_v=S_{0}\left(1-\Omega\right)``CpRVI`` (Compact-pol Radar Vegetation Index)
================================================
This functionality computes the compact-pol radar vegetation index for compact polarimetric SAR data. The required input and the computed output are as follows:

.. code-block:: python

        input : input_C2_folder, window_size
        output: CpRVI.bin

The formulation of the CpRVI is as follows:

.. math::
    
    \text{CpRVI}=\left(1-\dfrac{3}{2}\text{GD}_{\text{ID}}\right)\Big(\frac{p}{q}\Big)^{2(\frac{3}{2}\text{GD}_{\text{ID}})}\\
    p=\text{min\{SC,OC\}},q=\text{max\{SC,OC\}}\\
    \text{SC}=\frac{S_0-S_3}{2};\qquad{}\text{OC}=\frac{S_0+S_3}{2};\\
    S_0=\text{C11+C22};\qquad{}S_1=\text{C11-C22};\\
    S_2=\text{C12+C21};\qquad{}S_3=\pm\text{j(C12-C21)}\\


where, :math:`\text{GD}_\text{ID}` is the geodesic distance between Kennaugh matrices (:math:`\mathbf{K}`) of the observed and the ideal depolarizer, :math:`p, q` are minimum and maximum values of :math:`\text{SC}` and :math:`\text{OC}` which are functions of stocks parameters (:math:`S_0`, :math:`S_1`, :math:`S_2`, and :math:`S_3`). A detailed explanation of CpRVI is available in [[6]](#6).``DOP`` (Degree of Polarization)
=================================
This functionality computes the degree of polarization for compact polarimetric SAR data. The required input and the computed output are as follows:

.. code-block:: python

        input : input_c2_folder, window_size, tau
        output: DOP_CP.bin

The conventional degree of polarization in terms of stokes paramters can be written as follows:

.. math::

    \text{DOP}_{cp}=\frac{\sqrt{S^2_1+S^2_2+S^2_3}}{S_0}

where, 

.. math::
    
    S_0=\text{C11+C22};\qquad{}S_1=\text{C11-C22};\\
    S_2=\text{C12+C21};\qquad{}S_3=\pm\text{j(C12-C21)}



``PRVI`` (Polarimetric Radar Vegetation Index)
===============================================
This functionality computes the polarimetric radar vegetation index for dual polarimetric (HH | HV), (VV | VH) SAR data. The required input and the computed output are as follows:

.. code-block:: python

        input : input_c2_folder, window_size
        output: PRVI_dp.bin

The formulation of PRVI is as follows: 

.. math::

    \text{PRVI}_{dp}=\Big(1-\sqrt{1-\frac{4\times\text{det([C2])}}{\text{(Trace[C2])}^2}}\Big)\sigma^\circ_{\text{XY}}

where, :math:`\text{[C2]}` is co-variance matrix and :math:`\sigma^\circ_{\text{XY}}` is corss-pol backscatter intensity.
``DOP`` (Degree of Polarization)
=================================
This functionality computes the 2D Barakat degree of polarization for dual polarimetric (HH | HV), (VV | VH) SAR data. The required input and the computed output are as follows:

.. code-block:: python

    input : input_c2_folder, window_size
    output: dop_dp.bin

.. math::

    \text{DOP}_{dp}=\sqrt{1-\frac{4\times\text{det([C2])}}{\text{(Trace[C2])}^2}}

where, :math:`\text{[C2]}` is co-variance matrix. Further details on the Barakat Degree of polarization can be found in [[10]](#10)


``RVI`` (Radar Vegetation Index)
================================
This functionality computes the radar vegetation index for dual polarimetric (HH | HV), (VV | VH) SAR data. The required input and the computed output are as follows:

.. code-block:: python

        input : input_c2_folder, window_size
        output: RVI_dp.bin

The formulation of RVI is as follows:

.. math::

    \text{RVI}_{dp} = \frac{4 \times \sigma^\circ_{\text{XY}}}{\sigma^\circ_{\text{XX}}+\sigma^\circ_{\text{XY}}}  

where, :math:`\sigma^\circ_{\text{XX}}` is co-pol backscatter intensity and :math:`\sigma^\circ_{\text{XY}}` is corss-pol backscatter intensity``DpRVI`` (Dual-pol Radar Vegetation Index)
============================================
This functionality computes the dual polarimetric radar vegetation index for dual polarimetric (HH | HV), (VV | VH) SAR data. The required input and the computed output are as follows:

.. code-block:: python

        input : input_C2_folder, window_size
        output: DpRVI.bin

The formulation of DpRVI is as follows:

.. math::

    \text{DpRVI} = 1- \text{DOP}_{dp}\Big(\frac{\lambda_1}{\lambda_1+\lambda_2}\Big)

where,

.. math::
	\text{DOP}_{dp} = \sqrt{1-\frac{4\times \text{det ([C2])}}{\text{(Trace [C2])}^2}}

:math:`\text{[C2]}` is co-variance matrix,  and :math:`\lambda_1, \lambda_2` are the eigen values of :math:`\langle\mathbf{[C2]}\rangle` matrix in descending order. Further details on DpRVI can be obtained from [[5]](#5)


