[![Website monip.org](https://img.shields.io/website-up-down-green-red/http/monip.org.svg)](https://sevimcengiz.github.io/)
[![Documentation Status](https://readthedocs.org/projects/ansicolortags/badge/?version=latest)](https://sevimcengiz.github.io/oryx/)
![Open Source? Yes!](https://badgen.net/badge/Open%20Source%20%3F/Yes%21/blue?icon=github)
[![Ask Me Anything !](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)](https://github.com/sevimcengiz)
# Oryx-MRSI
<img src="https://user-images.githubusercontent.com/5468765/108315274-9c0c7680-71d4-11eb-9040-7e6248ea55b8.png" width="100" height="100">
Oryx-MRSI is a fully automated and complementary software for a comprehensive multi-slice proton magnetic resonance spectroscopic imaging (1H-MRSI) data analysis. It includes multi-slice MRSI raw data and LCModel .coord file output visualizations, tissue fraction calculation, chemical shift correction, metabolite maps generation, registration onto MNI152 brain atlas, and atlas-based ROI analysis.

# Features
- ```Main Page``` The user needs to provide parameters for the multi-slice 1H-MRSI data.
  Required parameters: 
 
   a. Exclusion criteria for fCSF, SNR, FWHM, CRLB 
 
   b. RF bandwidhth of the system for chemical shift correction 
 
   c. Cut-off value for the probabilistic binary map after registration
   
   d. Chemical shift correction is ```On``` or ```Off```
   
   e. RFOV dir is ```RL``` or ```AP```
   
   f. Chemical shift dir (AP) is ```A``` or ```P```
   
   g. Chemical shift dir (LR) is ```L``` or ```R```
   
   h. Chemical shift dir (FH) is ```F``` or ```H```
   
   i. Reference metabolite

- ```Load Data``` Reads the raw 1H-MRSI data and LCModel .coord otput files for raw data and Coord file visualization of spectra.

- ```Co-registration``` Coregisters FOV,Press-Box(VOI), all voxels of spectra considering chemical shift correction if chemical shift correction is 'On'.

- ```Segmentation``` Calculation of WM, GM, CSF fractions in each voxel of all metabolites considering chemical shift correction. 

- ```FWHM-SNR``` Visualization of FWHM and SNR values for all voxels.

- ```Spectral Quality``` Visualization of included voxels into the 1H-MRSI data analysis after exclusion criteria values considering FWHM, SNR, CRLB, and fCSF. 

- ```Metabolite Map``` All metabolite results are used to create metabolite maps including: 
  - concentration map,
  - concentration map to Ins ratio, 
  - concentration map to Cr+PCr ratio,
  - CSF corrected concentration map, 
  - CSF corrected concentration map to Ins ratio,
  - CSF corrected concentration map to Cr+PCr ratio.

- ```Registration``` Generates MNI152 brain atlas-Registered metabolite maps including the outputs created in the previous module.

- ```ROI Analysis``` Region of interest (ROI) analysis at multiple brain atlases like [Schafer2018 100/400 Parcels on 7 resting-state (rs-fMRI) networks](https://pubmed.ncbi.nlm.nih.gov/28981612/) or [MNI thr 0/25/50 brain regions](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases). 


# Prerequirements
- [MATLAB R2020b](https://www.fil.ion.ucl.ac.uk/spm/software/download/)
- [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSL) for FLIRT 
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/download/)
- [GUI Layout Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox)
- [Widgets Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/66235-widgets-toolbox-compatibility-support)
- Oryx-MRSI is tested on MAC (2.9 GHz Quad-Core Intel Core i7, 16 GB 2133 MHz LPDDR3, Radeon Pro 560 4 GB
Intel HD Graphics 630 1536 MB ) and Ubuntu 18.04.4 LTS (Memory 32GIB, Processor Intel Core i7-9800X CPU @3.8GHzx16, Graphics GeForce RTX 2070/PCle/SSE2)

# Installation
Oryx-MRSI uses FSL-Flirt function so using FSL from MATLAB should be ready. 

If you want to install FSL into your computer, check [this link](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)

If you use ```MAC```, check [this link](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/MacOsX) (Advance Usage part-Using FSL from MATLAB)

If you use ```LINUX```, check [this link](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/Linux) (Using FSL from MATLAB)

Plase download SPM12 using [this link](https://www.fil.ion.ucl.ac.uk/spm/software/download/)

Download Oryx-MRSI from Github repository,

Addpath Oryx-MRSI with subfolders.

Addpath SPM12 with subfolders.

Please make sure that FSL usage from Matlab command window installation is completed properly.
Before running a data analysis using Oryx-MRSI, let's check that FSL usage is from Matlab is done.

Please open matlab and run ```check_fsl_usage_from_matlab.m``` script which is given under Oryx-MRSI Github repo.

If there is no error, FSL usage from Matlab is completely installed.

If you get an error, plese check these:

   If you use ```MAC```, check [this link](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/MacOsX) (Advance Usage part-Using FSL from MATLAB)
    
   If you use ```LINUX```, check [this link](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation/Linux) (Using FSL from MATLAB)

# How to get started and User Documentation
[Oryx-MRSI Documentation](https://sevimcengiz.github.io/oryx/)

# Developers

Sevim Cengiz

Muhammed Yildirim

Abdullah Bas

Esin Ozturk Isik

Should you publish material that made use of Oryx-MRSI, please cite the following publication:

**Cengiz S, Yildirim M, Bas A, Ozturk-Isik E. ORYX-MRSI: A data analysis software for multi-slice 1H-MRSI. International Society for Magnetic Resonance in Medicine. Virtual Meeting, May 15-20, 2021. (Digital Poster)**

# Release
- Version 1.0

# Help and Support
- There isn't known any bug or issue up to now. 
- If you see any bug or issue, please  submit a topic in issues, or contact: sevim_cengiz@icloud.com
- If you support or contribute the code, most welcome to Oryx-MRSI Github Repository.

# License
- MIT License

# Acknowledgement
- This project was funded by TUBITAK 115S219. We thank all open-source MR and MRS tools. 
- Oryx-MRSI uses some functions of [FID-A](https://github.com/CIC-methods/FID-A), check [this link](https://github.com/CIC-methods/FID-A/blob/master/LICENSE.txt) for license.
- Oryx-MRSI uses some functions of [Gannet](https://github.com/richardedden/Gannet3.1).
    - Edden RAE, Puts NAJ, Harris AD, Barker PB, Evans CJ. Gannet: A batch-processing tool for the quantitative analysis of gamma-aminobutyric acid-edited MR      spectroscopy spectra. J. Magn. Reson. Imaging 2014;40:1445–1452. doi: 10.1002/jmri.24478)
- Oryx-MRSI uses some functions of [Osprey](https://github.com/schorschinho/osprey), check [this link](https://github.com/schorschinho/osprey/blob/develop/LICENSE.md) for license.
    - G Oeltzschner, HJ Zöllner, SCN Hui, M Mikkelsen, MG Saleh, S Tapper, RAE Edden. Osprey: Open-Source Processing, Reconstruction & Estimation of Magnetic Resonance Spectroscopy Data. J Neurosci Meth 343:108827 (2020).
- Oryx-MRSI uses some functions of [MRS_MRI_libs](https://github.com/chenkonturek/MRS_MRI_libs), check [this link](https://github.com/chenkonturek/MRS_MRI_libs/blob/master/LICENSE) for license.
- Oryx-MRSI uses some functionf of [NIFTI-Matlab](https://github.com/NIFTI-Imaging/nifti_matlab), check [this link](https://github.com/NIFTI-Imaging/nifti_matlab/blob/master/license.txt) for license.
- Oryx-MRSI uses some functions written by [Jamie Near](https://www.mcgill.ca/psychiatry/jamie-near) (McGill University)
- Oryx-MRSI uses some functions written by [H.Ratiney](https://www.creatis.insa-lyon.fr/site7/en/users/ratiney) (CREATIS-LRMN) 
- Oryx-MRSI uses Schaefer2018_100/400Parcels_7Networks_order_FSLMNI152_2mm.nii. See Github [link](https://github.com/ThomasYeoLab/Standalone_CBIG_fMRI_Preproc2016), for [license](https://github.com/ThomasYeoLab/Standalone_CBIG_fMRI_Preproc2016/blob/master/LICENSE.md).
- Oryx-MRSI uses MNI152_T1_2mm_brain.nii.gz, MNI-maxprob-thr0/25/50-2mm.nii.gz acquired from [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Atlases), for [license](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/License). 
- If there is any function that I forget to mention about name/link/citation, please let me know.

MR_libs
=======
### Description

These are MATLAB libraries for post-processing, analysing and simulating Magnetic Resonance Spectroscopy and Imaging (MRS &amp; MRI) data.

### Configuration

Data files & header files are required to be added to MATLAB path.  

### File List
**MRS_lib:** contains functions for post-processing, analysing and simulating Magnetic Resonance Spectroscopy data
* io/ 
  * **mrs_readSPAR.m**  : reads .SPAR Philips MRS header file
  * **mrs_readSDAT.m**  : reads .SDAT Philips MRS data file
  * **mrs_readLIST.m**  : reads .list Philips MRS header file produced by delayed reconstruction
  * **mrs_readDATA.m**  : reads .data Philips MRS data file produced by delayed reconstruction
  * **mrs_readSIN.m**   : reads .sin Philips raw MRS header file
  * **mrs_readRAW.m**   : reads .raw Philips raw MRS data file
  * **mrs_readRDA.m**   : reads .RDA Siemens MRS file (header information and data)
  * **mrs_readGEpfile.m** : reads Pxxxx.7 GE MRS raw data file (header information and data)
  * **mrs_readLcmodelBASIS.m**   : reads the LCModel input .basis file, which contains the basis set of model metabolite spectra 
  * **mrs_readLcmodelBasisRAW.m** : reads th LCModel input .RAW file, which contains time domain data of one metabolite spectrum 
  * **mrs_readLcmodelTABLE.m**   : reads the metabolite absolute and relative concentration and their SDs from the LCModel output .table file 
  * **mrs_readLcmodelCOORD.m**   : reads the LCmodel output .coord file, which contains the coordinates of all curves on the one-page output
  * **mrs_readLcmodelRAW.m** : reads LCModel output .RAW file which contains time domain data of each metabolite spectrum
  * **mrs_readJmruiTXT.m**   : reads .txt MRS data file from jMRUI
  * **mrs_writeSDAT.m** : writes MRS data to Philips .SDAT file
  * **mrs_writeSPAR.m** : writes MRS header information to Philips .SPAR file
  * **mrs_editSPAR.m**  : edits Philips .SPAR file
  * **mrs_writeLcmodelIN.m** : writes .in file for creating a LCmodel baisis set 
  * **mrs_writeLcmodelBasisRAW.m** : writes .RAW file, which contains time domain data of one metabolite spectrum    
* postprocess/
  * **mrs_truncate.m**  : truncates points from the end of spectra
  * **mrs_zerofill.m**  : fills zeros to the end of spectra
  * **mrs_apod.m**      : applies line-broadening filter 
  * **mrs_fft.m**       : applies fourier transformation 
  * **mrs_ifft.m**      : applies inverse fourier transformation 
  * **mrs_rephase.m**   : rephases spectra with specified phase value  
  * **mrs_zeroPHC.m**   : applies automatic zero-order phase correction to a spectrum
  * **mrs_manualzeroPHC.m**   : allows users to manually apply zero-order phase correction of a spectrum
  * **mrs_firstPHC.m**  : applies automatic first-order phase correction to a spectrum
  * **mrs_realign.m**   : aligns the target peaks in the spectra
  * **mrs_average.m**   : calculates averaged data
  * **mrs_lowpass.m**   : filters out high frequency components above specified frequency
  * **mrs_highpass.m**  : filters out low frequency components below specified frequency     
* simulation/  
  * **mrs_PRESS.m**      : simulates signal acquired using Position Resolved Spectroscopy (PRESS)
  * **mrs_STEAM.m**      : simulates signal acquired using Stimulated Echo Acquisition Mode (STEAM) 
  * **mrs_sLASER.m**     : simulates the signal acquired using semi-localised by adiabatic selective refocusing(sLASER) sequence 
  * **mrs_ISISscheme.m** : demonstrates how Image Selective in vivo Spectroscopy (ISIS) works
  * **mrs_simulateFID.m**: simulates a Free Induction Decay (FID) or Half-Echo. 
* utils/
  * **mrs_selectPeakrange.m** : allows users to manually pick the peak range interactively
  * **mrs_findPeak.m**        : locates the highest positive peak or lowest negative peak
  * **mrs_fitPeak.m**         : fits a peak in the given range of a spectrum with a lorenztian curve  
  * **mrs_lorentzFit.m**      : fits data with a Lorenztian function by minimising the squared error
  * **mrs_lorentzFun.m**      : defines the Lorentzian function
  * **mrs_gaussianFit.m**     : fits data with a Gaussian function by minimising the squared error  
  * **mrs_gaussianFun.m**     : defines Gaussian function function   
  * **mrs_points2Hz.m**       : converts unit from points to Hz
  * **mrs_points2ppm.m**      : converts unit from points to ppm
  * **mrs_ppm2Hz.m**          : converts unit from ppm to Hz
  * **mrs_plotSpectra.m**     : displays spectra 
  * **mrs_plotBASISspectra.m**: displays spectra in LCmodel .basis file  
  * **mrs_viewCSI.m**         : displays a spectrum from a selected voxel      
  * **mrs_rot90.m**           : rotates the spectra images 90 degree clockwise
  * **mrs_T1corr.m**          : applies T1 correction 
  * **mrs_T2corr.m**          : applies T2 correction 
  * **mrs_calFWHM.m**         : calculates the full width at haflf maximum height of the peak of interest  
  * **mrs_calTemp.m**         : calculates the temperature based on chemical shift difference of the water resonance and the temperature-independent reference resonance
* voxel_planner/
  * **create_stdMaskvoxel.m&.fig** : allows you to plan MRS voxel size & location with a guid of MR structural images of a standard brain (MNI152_T1_1mm_brain.nii.gz) and functional region masks (HarvardOxford-cort-maxprob-thr0-1mm.nii.gz) 
  
  ![Alt text](https://raw.github.com/chenkonturek/MR_libs/master/Images/MRS_voxel_planner.PNG)  

**MRI_lib** contains functions for post-processing, analysing and simulating Magnetic Resonance Imaging data
* io/
  * **mri_readHDR.m**  : reads .hdr MRI ANALYZE 7.5 header file 
  * **mri_readIMG.m**  : reads .img MRI ANALYZE 7.5 data file 
  * **mri_writeIMG.m** : writes data to .img file 
  * **mri_readPAR.m**  : reads .PAR Philips MRI header file
  * **mri_readREC.m**  : reads .REC Philips MRI data file
* utils/
  * **mri_IRcurve.m**        : calculates the inversion recovery (IR) curve
  * **mri_absIRcurveFit.m**  : fits the data to absolute inversion recovery curve 
  * **mri_createSVOI.m**     : creates a mask for locating the spectroscopic VOI in MR images based on given information of MR images and MRS voxel.
  * **mri_locateMask.m**     : creates a mask for locating the spectroscopic VOI in MR images based on Philips .PAR and .SPAR header files.
  * **mri_dispSVOI.m**       : displays the spectroscopic VOI on top of MR images

**NMR_lib:** contains functions for NMR simulation 
  * **nmr_bloch.m** : defines full Bloch equations in rotating frame
  * **nmr_calT2star.m** : calculates T2*
  * **nmr_getGamma** : returns gyromagnetic ratio value for different nuclei

**Examples:** contains example scripts. (Please email me if you want the data files) 
  * **example1.m** : demonstrates how to use MR_libs to post-process MRS data.
  * **example2.m** : demonstrates how to use MR_libs to locate spectroscopic voxel in MR images. 
  * **example3.m** : demonstrates how to use MR_libs to to simulate tye half-echo signal acquired using PRESS and STEAM sequences.respectively.
  * **example4.m** : demonstrates how to use MR_libs to do Bloch Equation simulation. 
  * **example5.m** : demonstrates how to use MR_libs to simulate magnetisation profile produced by an RF pulse.
  * **example6.m** : demonstrates how to use MR_libs to estimate T1.
  * **example7.m** : demonstrates how to use MR_libs to displays a spectrum from a selected voxel (CSI data).
  * **example8.m** : demonstrates how to use MR_libs to create .RAW for a LCModel model spectrum, and create .in file for making a LCmodel basisset.

### Acknowledgements

I would like to thank Professor Penny Gowland and Dr. Susan Francis for their supervision. 
I would also like to express my enormous thanks to Dr. Mary Stephensons and Dr. Emma Halls for their help and contributions.      


# MultiVoxelReaderWriter

for multivoxel spectroscopic data reader,you follow the following steps:

 1.) run multivoxel_readPhilips_Sdat.m file
 2.) a=fftshift(fft(sig,[],4),4);
 3.) overlay(a,number of slice, number of row, number of column
 4.) jMRUI_Metabolite_Map

rootname='/Users/sevimcengiz/Documents/CIL_LAB/PD_Analyze/PD_Patients/K_04/exam_1/spectra/K_04_press_1_raw_act';
[sig,Fe,Fref,B0]=Multivoxel_readPhilips_Sdat(rootname);
 a=fftshift(fft(sig,[],4),4);
%overlay(a,3,14,14);

An example usage:

[sig,Fe,Fref,B0]=Multivoxel_readPhilips_Sdat('your directory'); % It doesn't consis of .sdat/spar file extension at the end of the your directory.
a=fftshift(fft(sig,[],4),4);
overlay(a,3,14,14) % slice=3, row=14, column=14

Note 1.) If you want to look at a specific spectra with a determined row or column, you run the following code:
figure,plot(squeeze(real(a(row,column,:)))) %figure,plot(squeeze(real(a(2,2,:))))
 
Note 2.) If you want to change the line color or linewidth of the spectra figured, you got to overlay.m and change the color and linewidth.

Note 3.) If you want to see the metabolite map using LCModel concentration results, you go to LCModel_Metabolite_Map.m file.