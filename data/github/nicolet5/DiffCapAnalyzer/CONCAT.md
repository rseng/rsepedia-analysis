# DiffCapAnalyzer [![Build Status](https://travis-ci.com/nicolet5/DiffCapAnalyzer.svg?branch=master)](https://travis-ci.com/nicolet5/DiffCapAnalyzer) ![Python status](https://img.shields.io/badge/python-3.6%20%7C%203.7%20%7C%203.8-blue) [![status](https://joss.theoj.org/papers/543dd52dffe6a9551613a8533b25247a/status.svg)](https://joss.theoj.org/papers/543dd52dffe6a9551613a8533b25247a)

## Package for Battery Cycling Data Visualization and Analysis
This package is intended to quantitatively analyze raw cycling data by identifying and parameterizing peaks found in the differential capacity plots. Differential capacity plots (dQ/dV) can be very powerful for uncovering battery performance characteristics, as the peaks that appear in these plots correspond to the various electrochemical events. However, because of the large amount of data gathered during cycing experiments, many researchers report subsets of cycles and purely qualitative conclusions. This package adds the ability to quantify this type of analysis by cleaning battery cycling datasets and obtaining peak locations, widths, amplitudes, and other descriptors for every charge/discharge cycle in the data. To this end, this tool develops individualized fitted models for each charge and discharge cycle in the data set, comprised of a gaussian baseline and pseudo-voigt distributions at each peak location. 

Additionally, there is a DASH based visualization app that can be used as the user interface. Users can upload raw cycling data, either collected via a MACCOR or an Arbin cycler. The app will then process the data and add a few files to the database: the raw data, the cleaned data, and the peak descriptors for every cycle. The app also allows users to scroll through cycles and better understand the differential capacity curves. Additionally, there is a section to evaluate the fit of the gaussian baseline, and tailor the peak finding process. The user can also download the peak descriptors using the "Download CSV" file button in the app. 

Additionally, some machine learning was done to classify between two different cathode chemistries, LiCoO2 and LiFePO4. Data sets for these chemistries were obtained from the [CALCE website](https://web.calce.umd.edu/batteries/data.htm). Once this data was cleaned and labelled, a 20-80 test-train split was done and a support vector classifier was utilized, with a final test set accuracy of 77%. 

### Software Dependencies 
- Python3 
- For python packages see requirements.txt

### How to Install
To run the app and fully utilize DiffCapAnalyzer and the corresponding examples, simply clone this repo and from the root directory run: 
```
pip install -Ur requirements.txt
```
This will install all packages necessary for DiffCapAnalyzer. 

To use the DiffCapAnalyzer functions without the app, you can pip install instead of cloning the entire repo: 
```
pip install diffcapanalyzer 
```
This will install the DiffCapAnalyzer modules for use in the example notebooks, or for using the package outside of the Dash app. For example usage of the functions outside of the Dash app, see `examples/ProcessData_PlotDescriptors_Examples.ipynb`, which serves as a guide for processing one file and plotting the results through a jupyter notebook instead of the Dash App. Additionally, `exmaples/Detailed_Steps_Example.ipynb` provides additional insight into data cleaning, peak finding, and model fitting. 


## Dash App Tutorial
After cloning this repo annd installing the requirements, the Dash app can be used to interact with the underlying DiffCapAnalyer functions to process cycling data and visualize results. To start the app run the following command in terminal from the root directory:
```
python app.py
```
Which should return
```
 * Running on http://someurl/
```
Type or copy that URL in browser to launch the app locally. From the browser you should be able to see something like this: 

<img src="docs/images/app_initial_screen.png" alt="App Preview Image" width="1000"/>

The app initially loads with example data, as shown in the above image. You can explore that example data or you can upload new data. Example data sets to explore are found in this repo under `data/ARBIN/` and under `data/MACCOR/`. The Arbin data is obtained from the CALCE website (https://web.calce.umd.edu/batteries/data.htm), and are good examples of real dQ/dV data. The data under `MACCOR/` is a fabricated example, mainly meant to document the format of `MACCOR` data expected. Using one of the `ARBIN` data files, we can upload one by choosing the `ARBIN` option from the dropdown, and selecting our file for upload. Once the file is uploaded, the data is cleaned, each cycle has a model fit, and descriptors of each cycle are extracted. This process may take a few minutes: 

<img src="docs/images/app_loading_after_upload.png" alt="App Loading Image" width="1000"/>

NOTE: Any data you upload via this app remains on your local machine. The app stores the uploaded data in a local database found here: `data/databases/dQdV.db` 

Once the data is uploaded, you can explore the cycles by clicking through the slider. You can also zoom in and out on the plots and hover over the plot to see the data values: 

<img src="docs/images/app_explore_cycles.png" alt="App Explore Image" width="1000"/>

Scrolling down, you can explore the model a bit further. You can plot the peak descriptors, including peak locations, areas, and heights for both the positive and negative side of the differential capacity curve (Left Plot). 
Additionally, you can adjust the model fit, see a preview of model fit on a representative cycle, and update the model found in the database (Right Plot). The peak threshold can be altered in the app, which is a measure of which relative heights of peaks are identified as peaks. If the peak threshold is decreased, more peaks will be found (e.g. smaller peaks will be labelled as peaks), if it is increased, less peaks will be found (e.g. only the largest peaks will be identified as peaks). Currently, only the peak threshold is adjustable from the app, but other model fitting parameters could be surfaced in future versions. 

<img src="docs/images/app_descriptors_model.png" alt="App Model Image" width="1000"/>

The intended workflow in the above section of the app would be to try out various peak thresholds, update the preview of the model, and then once that preview looks good, update the model in the database for all the cycles. The resulting individual models can then be explored by reselecting the data from the dropdown and clicking through the cycles with the slider.

At the very bottom of the page, you can download the descriptors, cleaned cycles, and model points as CSVs, for exploration on your own. In the cleaned cycles CSV, the column titled `Smoothed_dQ/dV` is the processed (clean) cycling data. each column in the descriptors CSV is formatted as `sorted` descriptor type (e.g. FWHM, area, location, height) - `c` or `d` for charge/discharge - peak number. For example `sortedloc-c-2` would indicate the location of the second peak in the charge cyle. The `sorted` prefix is an indication these peaks were sorted, such that, e.g., if the first peak disappears in one of the cycles, the second peak does not become `peak 1`, it remains as `peak 2` in that cycle even though there is no `peak 1`. This is so this second peak is still associated with all of the other second peaks from other cycles. 

## Organization of the project
```
|   app.py
|   LICENSE
|   README.md
|   requirements.txt
|   runTests
|   setup.py
|   __init__.py
|
+---data
|   +---ARBIN
|   |   |   README.md
|   |   |
|   |   +---CS2_33
|   |   |
|   |   \---K2_016
|   |
|   +---databases
|   |       dQdV.db
|   |       init_database.db
|   |
|   +---MACCOR
|   |       example_data.csv
|   |
|   \---ML_data
|           c_descriptors.xlsx
|           descriptors_without_heights.xlsx
|           final_descriptors.xlsx
|           k_descriptors.xlsx
|           svc_model.sav
|           svc_results.png
|
+---diffcapanalyzer
|       app_helper_functions.py
|       chachifuncs.py
|       databasefuncs.py
|       databasewrappers.py
|       descriptors.py
|       __init__.py
|
+---docs
|   |   Poster.pdf
|   |
|   +---images
|   |       diagram.png
|   |
|   \---paper
|       |   paper.md
|       |
|       \---images
|               cleaning_dqdv.png
|               fitting_dqdv.png
|
+---examples
|   |   ProcessData_PlotDescriptors_Examples.ipynb
|   |
|   \---ML
|           SVC_Model.ipynb
|
\---tests
    |   test_app_helper_functions.py
    |   test_chachifuncs.py
    |   test_databasefuncs.py
    |   test_databasewrappers.py
    |   test_descriptors.py
    |   __init__.py
    |
    \---test_data
            test_data.csv
            test_data_mac.csv
```

## Running Tests
From the root directory: 
```
./runTests
```
If the executable file does not work, you can run the tests with the command: 

```
pytest --cov-report term --cov=diffcapanalyzer tests/
```


## Data Requirements
At the moment, the package can only process CSV files and relies on specific column headers for each type of file (Arbin vs. Maccor). Please reference the `data` directory for example files. The column headers for each data type must include and appear exactly as the following: 
* Arbin: 
    * Cycle_Index
    * Data_Point
    * Voltage(V)
    * Current(A)
    * Discharge_Capacity(Ah)
    * Charge_Capacity(Ah)
    * Step_Index
* MACCOR: 
    * Rec
    * Cycle C Step
    * TestTime
    * StepTime
    * Cap. [Ah]
    * Voltage [V]
    * Md
    * Current [A]

## Have Any Issues? Want to Contribute?
Please do so! Input and contributions from users is invaluable in shaping DiffCapAnalyzer to be the most useful for the community. If you use the software and notice a bug, issue, or think of a feature that would be great to have, report it via a [Github issue](https://github.com/nicolet5/DiffCapAnalyzer/issues). If you have modifications or new feature contributions to the DiffCapAnalyzer code, awesome! Submit a [pull request](https://github.com/nicolet5/DiffCapAnalyzer/pulls) with sufficient documentation on the changes, unit tests, and instructions on how to validate your changes/additions. For support with the tool, feel free to open a [Github issue](https://github.com/nicolet5/DiffCapAnalyzer/issues) or contact me via email at <nicole.thompson140@gmail.com>. ## Example ARBIN Data

The ARBIN data contained here was obtained by the Center for Advanced Life Cycle Engineering Battery Research Group. The data can be found online [here](https://web.calce.umd.edu/batteries/data.htm).
---
title: 'DiffCapAnalyzer: A Python Package for Quantitative Analysis of Total Differential Capacity Data'
tags:
  - Python
  - differential capacity
  - dQ/dV
  - cycling data
authors:
  - name: Nicole L. Thompson
    orcid: 0000-0003-3411-7373
    affiliation: 1
  - name: Theodore A. Cohen
    orcid: 0000-0001-7170-3211
    affiliation: "2, 3, 4"
  - name: Sarah Alamdari
    orcid: 0000-0002-4515-1262
    affiliation: 1
  - name: Chih-Wei Hsu
    orcid: 0000-0002-7328-5452
    affiliation: 1
  - name: Grant A. Williamson
    orcid: 0000-0003-2778-8304
    affiliation: 2
  - name: David A. C. Beck
    orcid: 0000-0002-5371-7035
    affiliation: "1, 5" 
  - name: Vincent C. Holmberg
    orcid: 0000-0002-9591-8951
    affiliation: "1, 2, 6" 
affiliations:
  - name: Dept. of Chemical Engineering, University of Washington
    index: 1
  - name: Molecular Engineering and Sciences Institute, University of Washington
    index: 2
  - name: Dept. of Materials Science and Engineering, University of Washington
    index: 3
  - name: Dept. of Chemistry, University of Washington
    index: 4
  - name: eScience Institute, University of Washington
    index: 5
  - name: Clean Energy Institute, University of Washington
    index: 6
date: 3 July 2020
bibliography: paper.bib
---

## Summary 
In order to study long-term degradation and charge storage mechanisms in batteries, researchers often cycle these electrochemical cells for hundreds or even thousands of charge and discharge cycles. The raw data produced during cycling can be interpreted via a variety of techniques that each highlight specific aspects of how the battery is functioning. Differential capacity (dQ/dV) analysis, one such technique, results in plots of the differential capacity – the charge introduced into the battery during a small change in voltage – _vs._ the voltage. Electrochemical reactions result in significant charge introduced into the cell across a small voltage window. In the differential capacity plot, this behavior results in a peak for each electrochemical reaction. Therefore, differential capacity plots are particularly useful for highlighting the various electrochemical events occurring within the cell, specific to each cycle [@marzocca_atwater; @torai_nakagomi_yoshitake_yamaguchi_oyama_2016; @aihara_ito_omoda_yamada_fujiki_watanabe_park_doo_2016; @christophersen_shaw_2010; @christophersen_bloom_thomas_gering_henriksen_battaglia_howell_2006; @weng_cui_sun_peng_2013]. In turn, these peaks carry important characteristics of the electrochemical reaction. For example, the location of the peak indicates at what voltage the reaction occurs, and the area of the peak is linked to the amount of charge exchanged in the reaction.

We present DiffCapAnalyzer, a Python package for extracting and tracking differential capacity curve features through multiple charge and discharge cycles. DiffCapAnalyzer provides cleaned dQ/dV curves, peak locations, peak heights, peak areas, and other characteristics specific to each cycle from raw battery cycling data.

## Statement of Need
Traditionally, when using differential capacity plots, researchers have drawn conclusions based on an arbitrarily chosen subset of cycles and reported mainly qualitative claims on how peaks shift during cycling, due to the difficulties in analyzing the full amount of data produced in the differential capacity plots. Additionally, although it is known that peak shapes and areas correlate to important electrochemical events, only a few papers report using peak deconvolution as a method to interpret dQ/dV plots [@torai_nakagomi_yoshitake_yamaguchi_oyama_2016; @aihara_ito_omoda_yamada_fujiki_watanabe_park_doo_2016; @bian_model_2019; @huang_incremental_2019; @he_comparative_2020]. Further, there does not exist any standardized method for peak deconvolution of differential capacity plots. These issues can largely be attributed to the lack of software designed for investigating sets of dQ/dV curves. Prior to DiffCapAnalyzer, no open source software has been available to researchers for the analysis of differential capacity curves through peak deconvolution.

## Description
The software described herein, DiffCapAnalyzer, has been developed to address the drawbacks associated with differential capacity analysis by processing cycling data in a chemistry-agnostic manner. This is done by calculating differential capacity from the given raw cycling data using Equation 1, cleaning and smoothing the dQ/dV plots, and performing automatic peak locating and deconvolution for every cycle within the dataset. 

$$(dQ/dV)_i=(Q_i-Q_{(i-1)})/(V_i-V_{(i-1)})\tag*{(1)}$$

In differential capacity curves without any cleaning or smoothing, there is significant noise and large step-wise changes present. This is a common problem when the denominator of Equation 1 approaches zero [@bloom_jansen_abraham_knuth_jones_battaglia_henriksen_2005; @huang_incremental_2019]. Therefore, in order to accurately identify peaks, the data is cleaned by removing points such that the voltage difference between datapoints is at least 0.001 V. Subsequently, the curve is smoothed using a Savitzky-Golay filter, which is a moving polynomial of specified order fit over a specified number of data points. At the current state of the software, the polynomial order of the Savitzky-Golay filter is set at 3 with a window length of 9 data points, as these seemed the best parameters on the data tested to preserve important features while removing noise. This cleaning process is summarized in \autoref{Figure 1}:
 
![Cleaning process on an example differential capacity curve.\label{Figure 1}](images/cleaning_dqdv.png)

Once the data is clean, the software automatically finds peaks in the dQ/dV curves utilizing the PeakUtils Python package [@peakutils], and returns the peak heights and the peak locations, as shown by an example cycle in \autoref{Figure 2}a. These peak heights and locations are then used to inform the model build, which is individualized to each cycle contained in the dataset. The model consists of Pseudo-Voigt distributions positioned at the identified peak locations and a baseline Gaussian distribution that captures the area that is not part of the Pseudo-Voigt distributions. The Pseudo-Voigt distribution described by Equations 2 and 3 is simply the linear combination of a Gaussian and a Lorentzian. This distribution is often used in fitting experimental spectral data due to it being a highly generalizable function able to fit a large variety of peak shapes [@wertheim_butler_west_buchanan_1974].

$$f_v(x,A,\mu,\sigma,\alpha)=\frac{(1- \alpha)A}{\sigma_g \sqrt{2 \pi}}\exp{[-{(x- \mu)}^2/2 {\sigma_g}^2]}+\frac{\alpha A}{\pi}[\frac{\sigma}{{(x-\mu)}^2 + \sigma^2}]\tag*{(2)}$$


$$\sigma_g = \sigma/\sqrt{2 \ln{2}}\tag*{(3)}$$


Once the model is generated, an optimized fit is found by allowing all parameters to vary except the center position of the Pseudo-Voigt peaks, which are assigned via the previously identified peak locations. \autoref{Figure 2}b presents an example of an initial model fit and the model fit once optimized specifically for that charge cycle.

![Fitting process on an example differential capacity curve.\label{Figure 2}](images/fitting_dqdv.png)

Further example model fits can be found on [GitHub](https://github.com/nicolet5/DiffCapAnalyzer/). From this model, peak areas, widths, and shapes can be extracted and examined to give further insight into the electrochemical processes occurring. The software also utilizes an SQLite database backend to store raw data, cleaned data, model parameters, and peak descriptors for each cycle. In addition to the data processing abilities of the software, a Dash-based web application has been developed where users can upload their own raw data to be processed and visualize the resulting dQ/dV plots and peak descriptors. Users can also evaluate the model fit, alter the threshold for peak detection, and update the model and descriptors in the database. From this application users can also download the cycle descriptors data as a CSV file for their own uses. Further instructions and descriptions of the software functionality can be found in the [DiffCapAnalyzer Github repo](https://github.com/nicolet5/DiffCapAnalyzer/).

In summary, DiffCapAnalyzer provides the ability to quantitatively analyze battery cycling data. The peak descriptors obtained from DiffCapAnalyzer could be used to classify electrochemical events, visualize battery degradation over time, or as features in state of health analyses. This package also lays the groundwork for a standardized method for cleaning and analyzing this type of data. We hope that this Python package advances the field of electrochemistry and enables researchers to better analyze, interpret, and present their battery cycling data.

## Acknowledgments
This project was supported by: Data Intensive Research Enabling Clean Technology (DIRECT) National Science Foundation (NSF) National Research Traineeship (DGE-1633216), the State of Washington through the University of Washington (UW) Clean Energy Institute and the UW eScience Institute, in part upon the work of S.A. and G.A.W. supported by the NSF Graduate Research Fellowship under Grants No. DGE-1762114 and DGE-1256082, respectively, and via funding from the Washington Research Foundation.

## References



