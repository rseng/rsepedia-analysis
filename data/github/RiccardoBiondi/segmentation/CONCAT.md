# Authors #

----------
- Riccardo Biondi - Dept. of Experimental, Diagnostic and Specialty Medicine of Bologna University ([riccardo.biondi4@studio.unibo.it](mailto:riccardo.biondi4@studio.unibo.it))
- Nico Curti - eDIMESLab, Dept. of Experimental, Diagnostic and Specialty Medicine of Bologna University ([nico.curti2@unibo.it](mailto:nico.curti2@unibo.it))
- Enrico Giampieri - eDIMESLab, Dept. of Experimental, Diagnostic and Specialty Medicine of Bologna University ([enrico.giampieri@unibo.it](mailto:enrico.giampieri@unibo.it))
- Gastone Castellani - Dept. of Experimental, Diagnostic and Specialty Medicine of Bologna University ([gastone.castellani@unibo.it](mailto:gastone.castellani@unibo.it))
| **Authors**  | **Project** |  **Build Status** | **License** | **Code Quality** | **Coverage** |
|:------------:|:-----------:|:-----------------:|:-----------:|:----------------:|:------------:|
| [**R. Biondi**](https://github.com/RiccardoBiondi) <br/> [**N. Curti**](https://github.com/Nico-Curti) | **COVID-19 Lung Segmentation** [![status](https://joss.theoj.org/papers/31abd09499e0535e2d65cd40f4cb1766/status.svg)](https://joss.theoj.org/papers/31abd09499e0535e2d65cd40f4cb1766)| **Linux** : [![Build Status](https://travis-ci.com/RiccardoBiondi/segmentation.svg?branch=master)](https://travis-ci.com/RiccardoBiondi/segmentation) <br/>  **Windows** : [![Build status](https://ci.appveyor.com/api/projects/status/om6elsnkoi22xii3?svg=true)](https://ci.appveyor.com/project/RiccardoBiondi/segmentation) | [![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://github.com/RiccardoBiondi/segmentation/blob/master/LICENSE.md) | **Codacy** : [![Codacy Badge](https://app.codacy.com/project/badge/Grade/cc0fd47ae8e44ab1943b1f74c2a3d7e2)](https://www.codacy.com/manual/RiccardoBiondi/segmentation?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=RiccardoBiondi/segmentation&amp;utm_campaign=Badge_Grade) <br/> **Codebeat** : [![CODEBEAT](https://codebeat.co/badges/927db14b-36fc-42ed-88f1-09b2a9e1b9c0)](https://codebeat.co/projects/github-com-riccardobiondi-segmentation-master) | [![codecov](https://codecov.io/gh/RiccardoBiondi/segmentation/branch/master/graph/badge.svg)](https://codecov.io/gh/RiccardoBiondi/segmentation) |

[![Project CI](https://github.com/RiccardoBiondi/segmentation/workflows/CTLungSeg%20CI/badge.svg)](https://github.com/RiccardoBiondi/segmentation/actions/workflows/python.yml)
[![Docs CI](https://github.com/RiccardoBiondi/segmentation/workflows/CTLungSeg%20Docs%20CI/badge.svg)](https://github.com/RiccardoBiondi/segmentation/actions/workflows/docs.yml)

[![docs](https://readthedocs.org/projects/covid-19-ggo-segmentation/badge/?version=latest)](https://covid-19-ggo-segmentation.readthedocs.io/en/latest/?badge=latest)
[![GitHub pull-requests](https://img.shields.io/github/issues-pr/RiccardoBiondi/segmentation.svg?style=plastic)](https://github.com/RiccardoBiondi/segmentation/pulls)
[![GitHub issues](https://img.shields.io/github/issues/RiccardoBiondi/segmentation.svg?style=plastic)](https://github.com/RiccardoBiondi/segmentation/issues)

[![GitHub stars](https://img.shields.io/github/stars/RiccardoBiondi/segmentation.svg?label=Stars&style=social)](https://github.com/RiccardoBiondi/segmentation/stargazers)
[![GitHub watchers](https://img.shields.io/github/watchers/RiccardoBiondi/segmentation.svg?label=Watch&style=social)](https://github.com/RiccardoBiondi/segmentation/watchers)

# COVID-19 Lung Segmentation

This package allows to isolate the lung region and identify ground glass lesions
on chest CT scans of patients affected by COVID-19.
The segmentation approach is based on color quantization, performed by K-means
clustering.
This package provides a series of scripts to isolate lung regions, pre-process
the images, estimate K-means centroids and labels of the lung regions.

1. [Overview](#Overview)
2. [Contents](#Contents)
3. [Prerequisites](#Prerequisites)
4. [Installation](#Installation)
5. [Usage](#Usage)
6. [License](#License)
7. [Contribution](#Contribution)
8. [References](#References)
9. [Authors](#Authors)
10. [Acknowledgments](#Acknowledgments)
11. [Citation](#Citation)

## Overview

COronaVirus Disease (COVID-19) has widely spread all over the world since the
beginning of 2020. It is acute, highly contagious, viral infection mainly
involving the respiratory system. Chest CT scans of patients affected by this
condition have shown peculiar patterns of Ground Glass Opacities (GGO) and Consolidation (CS) related to the severity and the stage of the disease.

In this scenario, the correct and fast identification of these patterns is a
fundamental task. Up to now this task is performed mainly using manual or
semi-automatic techniques, which are time-consuming (hours or days) and
subjected to the operator expertise.

This project provides an automatic pipeline for the segmentation of
GGO areas on chest CT scans of patient affected by COVID-19.
The segmentation is achieved with a color quantization algorithm, based on
k-means clustering, grouping voxel by color and texture similarity.

**Example of segmentation**. **Left:** Original image: **Right** original image with identified ground-glass areas.

<a href="https://github.com/RiccardoBiondi/segmentation/blob/master/images/results.png">
  <div class="image">
    <img src="https://github.com/RiccardoBiondi/segmentation/blob/master/images/results.png" width="800" height="400">
  </div>
</a>

The pipeline was tested on 15 labeled chest CT scans, manually segmented by
expert radiologist.
The goodness of the segmentation was estimated using Dice(0.67 ± 0.12),
Sensitivity(0.666 ± 0.15), Specificity(0.9993 ± 0.0005) and
Precision(0.75± 0.20) scores.

These results make the pipeline suitable as initialization for more accurate
methods

## Contents

COVID-19 Lung segmentation is composed of scripts and modules:
- scripts allows to isolate lung regions, find the centroids for colour quantization and segment the images.
- modules allows to load and save the images from and to different extensions and perform operations on image series.

To refer to script documentation:

| **Script** | **Description** |
|:----------:|:---------------:|
| [lung_extraction](https://covid-19-ggo-segmentation.readthedocs.io/en/latest/script.html#lung-extraction) | Extract lung from CT scans 										 																																				|
| [train](https://covid-19-ggo-segmentation.readthedocs.io/en/latest/script.html#train) | Apply colour quantization on a series of stacks to estimate the centroid to use for segmentation																																													|
| [labeling](https://covid-19-ggo-segmentation.readthedocs.io/en/latest/script.html#labeling) |Segment the input image by using pre-estimated centroids or user-provided set|

To refer to modules documentation:

| **Module**| **Description**|
|:---------:|:--------------:|
| [utils](https://covid-19-ggo-segmentation.readthedocs.io/en/latest/modules.html#utils) | method to load, save and preprocess stack																																										|
| [method](https://covid-19-ggo-segmentation.readthedocs.io/en/latest/modules.html#method) | method to filter the image tensor |
| [segmentation](https://covid-19-ggo-segmentation.readthedocs.io/en/latest/modules.html#segmentation) | contains useful function to segment stack of images and select ROI																										|

For each script described below, there are a PowerShell and a shell script that
allows their execution on multiple patients scans. Moreover it also provide a
snakemake pipeline.

## Prerequisites

Supported python version: ![Python version](https://img.shields.io/badge/python-3.5|3.6|3.7|3.8-blue.svg)

First of all ensure to have the right python version installed.

This script use opencv-python, numpy and SimpleITK: see
[requirements](https://github.com/RiccardoBiondi/segmentation/blob/master/requirements.txt)
for more informations.

The lung extraction is performed by using a pre-trained UNet, so please ensure to
have installed the [lungmask](https://github.com/JoHof/lungmask) package.
For more information about how the network is trained, please refer
to https://doi.org/10.1186/s41747-020-00173-2.

> :warning: The OpenCV requirement binds the minimum Python version of this project
> to Python 3.5!

To run the tests you need to install ```PyTest``` and ```Hypothesis```.
Installation instructions are available at: [PyTest](https://docs.pytest.org/en/6.2.x/getting-started.html), [Hypothesis](https://docs.pytest.org/en/6.2.x/getting-started.html)

## Installation

Download the project or the latest release:

```bash
git clone https://github.com/RiccardoBiondi/segmentation
```

Now you can install the package using pip:

```bash
pip install segmentation/
```

### Testing

Testing routines use ```PyTest``` and ```Hypothesis``` packages. please install
these packages to perform the test. o install the package in development mode you need to add also this requirement:

- pytest >= 3.0.7

- hypothesis >= 4.13.0

> :warning: pytest versions above 6.1.2 are not supported by python 3.5

A full set of test is provided in [testing](https://github.com/RiccardoBiondi/segmentation/blob/master/testing) directory.
You can run the full list of test with:

```bash
python -m pytest
```

## Usage

This modules provides some script to segment a single scan, to automate the segmentation for multiple patients and to train your centroid set.
In the following paragraph, we will see how to use all the features. To achieve this purpose,
we will use, as example, the public dataset *COVID-19 CT Lung and Infection Segmentation Dataset*, published by Zenodo[5].

### Download Data

Firstly, we have to download and prepare the data.
All the data will be stored and organized in a folder named *Example*.

Download data into the Examples folder

using Bash:

```bash
  $ mkdir Examples
  $ wget https://zenodo.org/record/3757476/files/COVID-19-CT-Seg_20cases.zip -P ./Examples
  $ unzip ./Examples/COVID-19-CT-Seg_20cases.zip -d ./Examples/COVID-19-CT
```

Or PowerShell:

```PowerShell

  PS \> New-Item  -Path . -Name "Examples" -ItemType "directory"
  PS \> Start-BitsTransfer -Source https://zenodo.org/record/3757476/files/COVID-19-CT-Seg_20cases.zip -Destination .\Examples\
  PS \> Expand-Archive -LiteralPath .\Examples\COVID-19-CT-Seg_20cases.zip -DestinationPath .\Examples\COVID-19-CT -Force
```

### Single Scan
Once you have download the data and installed the module, you can start to segment the images.
Input CT scans must be in Hounsfield units(HU) since grey-scale
images are not allowed.
The input allowed formats are the ones supported by SimpleITK.
If the input is a DICOM series, pass the path to the directory containing
the series files.
Please ensure that the folder contains only one series.
As output will save the segmentation as *nrrd*.

To segment a single CT scan run the following from the bash or PowerShell:

```bash
   python -m CTLungSeg --input='./Examples/COVID-19-CT/coronacases_003.nii.gz'  --output='./Examples/coronacases_003_label.nrrd'
```

### Multiple Scans

In the case of multiple patients segmentation, you have to repeat the segmentation process many times:  We have automated this process using bash(for Linux) and PowerShell(for Windows) scripts.
We have also provided a snakemake pipeline for the whole segmentation procedure in a multi-processing environment.
In the following paragraph, we will explain how to organize your data to benefits from this automation.

#### Script

To run the scripts,, you have to organize the data into three folders:

- input folder: contains all and only the CT scans to segment
- temporary folder: empty folder. Will contain the scans after the lung segmentation
- output folder: empty folder, will contain the labels files.

As examples we will segmenta the *coronacases_002* and the *coronacases_005* patients.

From bash:

```bash
  $ mkdir ./Examples/INPUT
  $ mkdir ./Examples/LUNG
  $ mkdir ./Examples/OUTPUT
  $ mv ./Examples/COVID-19-CT/coronacases_002.nii.gz ./Examples/COVID-19-CT/coronacases_005.nii.gz ./Examples/INPUT
```
or from PowerShell

```PowerShell
  PS \> New-Item -Path "Examples" -Name "INPUT" -ItemType "directory"
  PS \> New-Item -Path "Examples" -Name "LUNG" -ItemType "directory"
  PS \> New-Item -Path "Examples" -Name "OUTPUT" -ItemType "directory"
  PS \> Move-Item -Path "Examples\COVID-19-CT\coronacases_002.nii.gz" -Destination "Examples\INPUT"
  PS \> Move-Item -Path "Examples\COVID-19-CT\coronacases_005.nii.gz" -Destination "Examples\INPUT"
```

Now you can proceed with the **lung segmentation**. To achieve this purpose run
from PowerShell the  script:

 ```PowerShell
  PS \> ./lung_extraction.ps1 ./Examples/INPUT ./Examples/LUNG
 ```

Or its equivalent bash version:

```bash
  $ ./lung_extraction.sh./Examples/INPUT ./Examples/LUNG
```

Once you have successfully isolated the lung, you are ready to perform the GGO
segmentation. Run the labelling scrip from PowerShell :

```PowerShell
  PS /> ./labeling.ps1 ./Examples/LUNG ./Examples/OUTPUT
```

Or its corresponding bash version:

```bash
$ ./labeling.sh ./Examples/LUNG ./Examples/OUTPUT
```

##### Train your own centroid set

It is possible to train your centroid set instead of using the pre-trained one.

In this case you have to prepare these folders :
  - TRAIN : will contain the scans in the training set
  - TLUNG : will stores the scans after lung extraction

We will use *coronaceses_003* and *coronaceses_008* as training set.

From bash:

```bash
  $ mkdir ./Examples/TRAIN
  $ mkdir ./Examples/TLUNG
  $ mv ./Examples/COVID-19-CT/coronacases_003.nii.gz ./Examples/COVID-19-CT/coronacases_008.nii.gz ./Examples/TRAIN
```

or Powershell:

```PowerShell
  PS \> New-Item -Path ".\Examples" -Name "TRAIN" -ItemType "directory"
  PS \> New-Item -Path ".\Examples" -Name "TLUNG" -ItemType "directory"
  PS \> Move-Item -Path ".\Examples\COVID-19-CT\coronacases_003.nii.gz" -Destination "Examples\TRAIN"
  PS \> Move-Item -Path ".\Examples\COVID-19-CT\coronacases_008.nii.gz" -Destination "Examples\TRAIN"
```

First of all, you have to perform the lung extraction on the train scans,
as before run:

```bash
  $ ./lung_extraction.sh ./Examples/TRAIN/ ./Examples/TLUNG/
```

or its corresponding PowerShell version. Now, to estimate the centroid set, run:

```bash
  $ ./train.sh ./Examples/TLUNG/ ./centroid.pkl.npy
```

or its corresponding PowerShell version.

#### Snakemake

If you have not installed snakemake, you can find the instruction [here](https://snakemake.readthedocs.io/en/stable/).
To use the snakemake pipeline, you have to create two folders:

  - INPUT : contains all and only the CT scans to segment
  - OUTPUT : empty folder, will contain the segmented scans as *nrrd*.

As before we will use as examples *coronacases_002* and *coronacases_005* patients

> :notes: If you already run the script version, these folder are ready

Execute from bash

```bash
  $ mkdir ./Examples/INPUT
  $ mkdir ./Examples/OUTPUT
  $ mv ./Examples/COVID-19-CT/coronacases_002.nii.gz ./Examples/COVID-19-CT/coronacases_005.nii.gz ./Examples/INPUT
```

or PowerShell

```PowerShell
  PS \> New-Item -Path "Examples" -Name "INPUT" -ItemType "directory"
  PS \> New-Item -Path "Examples" -Name "OUTPUT" -ItemType "directory"
  PS \> Move-Item -Path ".\Examples\COVID-19-CT\coronacases_002.nii.gz" -Destination "Examples\INPUT"
  PS \> Move-Item -Path ".\Examples\COVID-19-CT\coronacases_005.nii.gz" -Destination "Examples\INPUT"

```

Now, from command line, execute:

```bash
  snakemake --cores 1 --config input_path='./Examples/INPUT/'
  output_path='./Examples/OUTPUT/'
```

> :notes: This command works both for Bash and Powershell

> :warning: It will create a folder named **LUNG** inside the INPUT,
> which contains the results of the lung extraction step.

#### Train Your Centroids

As before, you can decide to train your centroid set. To achieve this purpose, using the snakemake pipeline, you have to prepare three folders :

  - INPUT: will contains all the scans to segment
  - OUTPUT: will contain the segmented scans
  - TRAIN: will contain all the scans of the training set. (**NOTE** Cannot be the INPUT folder)

> :warning: INPUT and TRAIN folder cannot be the same

> :notes: This will train the centroid set, and after that perform the segmentation on the scans in the input folder.
> So the INPUT folder is organized as before.

Now run Snakemake with the following configuration parameters :

```bash
  snakemake --cores 1 --config input_path='./Examples/INPUT/'
  output_path='.Examples/OUTPUT/' train_path='./Examples/TRAIN/' centroid_path='./Examples/centorids.pkl.npy'
```

## License

The `COVID-19 Lung Segmentation` package is licensed under the MIT "Expat" License.
[![License](https://img.shields.io/github/license/mashape/apistatus.svg)](LICENSE.md)

## Contribution

Any contribution is more than welcome.
Just fill an [issue](./.github/ISSUE_TEMPLATE/ISSUE_TEMPLATE.md) or a
[pull request](./.github/PULL_REQUEST_TEMPLATE/PULL_REQUEST_TEMPLATE.md)
and we will check ASAP!

See [here](https://github.com/RiccardoBiondi/segmentation/blob/master/CONTRIBUTING.md)
for further informations about how to contribute with this project.

## References

<blockquote>1- Hofmanninger, J., Prayer, F., Pan, J. et al. Automatic lung segmentation in routine imaging is primarily a data diversity problem, not a methodology problem. Eur Radiol Exp 4, 50 (2020). https://doi.org/10.1186/s41747-020-00173-2. </blockquote>

<blockquote>2- Bradski, G. (2000). The OpenCV Library. Dr. Dobb&#x27;s Journal of Software Tools.</blockquote>

<blockquote>3- Yaniv, Z., Lowekamp, B.C., Johnson, H.J. et al. SimpleITK Image-Analysis Notebooks: a Collaborative Environment for Education and Reproducible Research. J Digit Imaging 31, 290–303 (2018). https://doi.org/10.1007/s10278-017-0037-8.</blockquote>

<blockquote>4- Lowekamp Bradley, Chen David, Ibanez Luis, Blezek Daniel The Design of SimpleITK  Frontiers in Neuroinformatics 7, 45 (2013) https://www.frontiersin.org/article/10.3389/fninf.2013.00045.</blockquote>

<blockquote>5- Ma Jun, Ge Cheng, Wang Yixin, An Xingle, Gao Jiantao, Yu Ziqi, Zhang Minqing, Liu Xin, Deng Xueyuan, Cao Shucheng, Wei Hao, Mei Sen, Yang Xiaoyu, Nie Ziwei, Li Chen, Tian Lu, Zhu Yuntao, Zhu Qiongjie, Dong Guoqiang, & He Jian. (2020). COVID-19 CT Lung and Infection Segmentation Dataset (Verson 1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.3757476.</blockquote>

## Authors

* <img src="https://avatars3.githubusercontent.com/u/48323959?s=400&v=4" width="25px"> **Riccardo Biondi** [git](https://github.com/RiccardoBiondi)

* <img src="https://avatars0.githubusercontent.com/u/24650975?s=400&v=4" width="25px"> **Nico Curti** [git](https://github.com/Nico-Curti), [unibo](https://www.unibo.it/sitoweb/nico.curti2)

* <img src="https://avatars2.githubusercontent.com/u/1419337?s=400&v=4" width="25px;"/> **Enrico Giampieri** [git](https://github.com/EnricoGiampieri), [unibo](https://www.unibo.it/sitoweb/enrico.giampieri)

* <img src="https://www.unibo.it/uniboweb/utils/UserImage.aspx?IdAnagrafica=236217&IdFoto=bf094429" width="25px;"/> **Gastone Castellani** [unibo](https://www.unibo.it/sitoweb/gastone.castellani)

See also the list of [contributors](https://github.com/RiccardoBiondi/segmentation/contributors) [![GitHub contributors](https://img.shields.io/github/contributors/RiccardoBiondi/segmentation.svg?style=plastic)](https://github.com/RiccardoBiondi/segmentation/graphs/contributors/) who participated to this project.

## Acknowledgments

The authors acknowledge all the members of the Department of Radiology, IRCCS Azienda
Ospedaliero-Universitaria di Bologna and the SIRM foundation, Italian Society of Medical and
Interventional Radiology for the support in the development of the project and analysis of the data.

## Citation

If you have found `COVID-19 Lung Segmentation` helpful in your research, please
consider citing the original paper

```BibTeX
@article{app11125438,
  author = {Biondi, Riccardo and Curti, Nico and Coppola, Francesca and Giampieri, Enrico and Vara, Giulio and Bartoletti, Michele and Cattabriga, Arrigo and Cocozza, Maria Adriana and Ciccarese, Federica and De Benedittis, Caterina and Cercenelli, Laura and Bortolani, Barbara and Marcelli, Emanuela and Pierotti, Luisa and Strigari, Lidia and Viale, Pierluigi and Golfieri, Rita and Castellani, Gastone},
  title = {Classification Performance for COVID Patient Prognosis from Automatic AI Segmentation—A Single-Center Study},
  journal = {Applied Sciences},
  volume = {11},
  year = {2021},
  number = {12},
  article-number = {5438},
  url = {https://www.mdpi.com/2076-3417/11/12/5438},
  issn = {2076-3417},
  doi = {10.3390/app11125438}
}
```

or just this project

```BibTeX
@misc{COVID-19 Lung Segmentation,
  author = {Biondi, Riccardo and Curti, Nico and Giampieri, Enrico and Castellani, Gastone},
  title = {COVID-19 Lung Segmentation},
  year = {2020},
  publisher = {GitHub},
  howpublished = {\url{https://github.com/RiccardoBiondi/segmentation}},
}

```
---
title: 'COVID-19 Lung Segmentation'
tags:
  - radiomics
  - artificial-intelligence
  - machine-learning
  - deep-learning
  - medical-imaging
  - chest-CT
  - python3
authors:
  - name: Riccardo Biondi^[co-first author]
    orcid: 0000-0003-0872-7098
    affiliation: 1
  - name: Nico Curti^[co-first author]
    orcid: 0000-0001-5802-1195
    affiliation: 2
  - name: Enrico Giampieri
    orcid: 0000-0003-2269-2338
    affiliation: 2
  - name: Gastone Castellani
    orcid: 0000-0003-4892-925X
    affiliation: 1

affiliations:
  - name : Department of Experimental, Diagnostic and Specialty Medicine of Bologna University
    index: 1
  - name: eDIMESLab, Department of Experimental, Diagnostic and Specialty Medicine of Bologna University
    index: 2
date: 16/06/2021
bibliography: paper.bib
---

# Summary

The `COVID-19 Lung Segmentation` project provides a novel, unsupervised and
fully automated pipeline for the semantic segmentation of ground-glass opacity
(GGO) areas in chest Computer Tomography (CT) scans of patients affected by COVID-19.
In the project we provide a series of scripts and functions for the automated
segmentation of lungs 3D areas, segmentation of GGO areas, and estimation of
radiomic features.

Both PowerShell and bash scripts are provided for the scripts management.
A possible Snakemake pipeline for the whole segmentation procedure applied
to several CT scans (in a multi-processing environment) is included into
the project.

A detailed description of the whole pipeline of processing has been already discussed
in @app11125438, where we have showed also the results obtained on public
datasets [@zenodo].
In that work we proved the efficiency of the proposed unsupervised method for the
identification of GGO areas and extraction of informative radiomic features.
Radiomic features were collected and used to predict clinically relevant
scores, with particular focus on mortality and the PREDI-CO score
[@Bartoletti2020].

# Statement of Need

COronaVirus Disease (COVID-19) has widely spread all over the world since the
beginning of 2020.
It is an acute, highly contagious, viral infection mainly involving the respiratory system.
Chest CT scans of patients affected by this condition have shown peculiar patterns
of Ground Glass Opacities (GGO) and Consolidation (CS) related to the severity
and the stage of the disease.

The correct and fast identification of these patterns is a fundamental task.
Up to now, this task has mainly been performed using manual or semi-automatic techniques,
which are time-consuming (hours or days), with results dependent on the operator's expertise.

This project provides an automated pipeline for the segmentation of
GGO areas on chest CT scans of patient affected by COVID-19.
The segmentation is achieved with a color quantization algorithm, based on k-means
clustering, which groups the voxels by color and texture similarity. This
approach is preceeded by the lung segmentation, achieved by a public available
U-Net model [@Hofmanninger2020;@lungmask].

The pipeline's performance has been tested on a dataset of 15 labeled chest CT scans.
These scans were segmented and validated by an expert radiologist.
Ten of these scans were extracted from the public dataset
*COVID-19 CT Lung and Infection Segmentation Dataset* [@zenodo]
published on Zenodo.
The Department of Diagnostic and Preventive Medicine of the IRCCS Policlinic Sant'Orsola-Malpighi
provided another 82 scans, with the 5 labeled scans used for the evaluation.

We tested the segmentation performances using the dice coefficient and specificity,
sensitivity, and precision scores.
The average value and the corresponding standard deviation at $1\sigma$ are reported in
the following table.

|  Dice Score  |  Sensitivity |    Specificity   |   Precision  |
|:------------:|:------------:|:----------------:|:------------:|
|$0.67\pm 0.12$|$0.66\pm 0.15$|$0.9992\pm 0.0005$|$0.75\pm 0.20$|

The proposed unsupervised segmentation pipeline is able to approximate the gold
standard with satisfactory results.
Given that the amount of information required for the k-means method training is considerably lower than for CNN methods, while still retaining good results, this segmentation can be implemented with in-patient training [@app11125438];
as a reference, a 3D U-Net-based method [@yan2020covid19] required two order of magnitude more training samples to achieve comparable results.
With this work we aimed to prove that semi-supervised approaches to segmentation are promising,
as they would combine the best effort of highly trained physicians to develop true gold standard
segmentation and the expertise of data analysts to augment those segmentation in full blown models.
While the proposed pipeline is not yet at the accuracy level necessary for assisted diagnostics,
we surmise that our pipeline can be successfully used as a first segmentation method to be used as training for other, more specific methods.

# Acknowledgments

The authors acknowledge all the members of the Department of Radiology, IRCCS Azienda
Ospedaliero-Universitaria di Bologna and the SIRM foundation, Italian Society of Medical and
Interventional Radiology for the support in the development of the project and analysis of the data.

# References
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to make participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

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

This Code of Conduct applies within all project spaces, and it also applies when
an individual is representing the project or its community in public spaces.
Examples of representing a project or community include using an official
project e-mail address, posting via an official social media account, or acting
as an appointed representative at an online or offline event. Representation of
a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at riccardo.biondi4@studio.unibo.it. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
<!-- Please use this line to close one or multiple issues when this pullrequest gets merged
You can add another line right under the first one:
resolves #1234
resolves #1235
-->

#### This PR changes :grey_question:
<!-- Explain your changes -->

#### Any other comments?
<!--
This is a template helping you to create an issue which can be processed as quickly as possible. This is the bug reporting section for the NumPyNet library.
-->

#### Detailed Description :grey_question:
<!-- your description -->

#### Steps/Code to Reproduce :grey_question:
<!-- to add code example fence it with triple backticks and optional file extension
    ```.py
    // Python code example
    ```
 or attach as .txt or .zip file
-->

#### Expected Behavior :grey_question:
<!-- Description of the expected result(s) -->

#### Actual Behavior :grey_question:
<!-- Description (possibly with some shell reports) of the actual result(s) -->

#### Operating System / Platform :grey_question:
<!-- Example
- OS: Windows 10 Pro
- System type: x64
- Processor: i7-6500U
- RAM: 8 GB
-->

#### Python Version :grey_question:
<!--
$ python --version
Python 3.8.2
-->

#### COVID-19 Lung Segmentation Version (`segmentation.__version__`) :grey_question:
<!--
$ python -c "import segmentation; print(segmentation.__version__)"
'1.0.0'
-->
Script
======
Since the segmentation of several patients is time-consuming, some scripts
to automatize this process are provided. The scripts are in two versions: bash and
PowerShell. The segmentation approach is different, instead of to perform the
entire segmentation of all the patient all at once, the segmentation steps are
divided into two different steps :

- lung extraction
- labeling

The following examples will use the data previously downloaded from the public dataset.

Lung Extraction
---------------

This script allows to run the lung segmentation on the whole set of patients,
that is the preliminary step of the GGO identification. Its implemented both for
bash and PowerShell.

In this example we will segment *coronacases_002* and *coronacases_005* patients.
To perform the segmentation, you have to create two folers:

  - input: contains only the scan to segment
  - lung: contains the lung extraction results.



Ensure that in the input folder there are only the files corresponding to the scan
to segment. The supported input format are all the ones supported by SimpleITK_.
you can provide as input also DICOM series, simply arrange them into one folder
for each scan. Please ensure that all the subsamples contain only one series.

Now you can simply run the following command from bash :

.. code-block:: bash

  mkdir ./Examples/INPUT
  mkdir ./Examples/LUNG
  mv ./Examples/COVID-19-CT/coronacases_002.nii.gz ./Examples/COVID-19-CT/coronacases_005.nii.gz ./Examples/INPUT
  lung_extraction.sh ./Examples/INPUT ./Examples/LUNG

or its equivalent for powershell

.. code-block:: powershell

  New-Item -Path "Examples" -Name "INPUT" -ItemType "directory"
  New-Item -Path "Examples" -Name "LUNG" -ItemType "directory"
  Move-Item -Path "Examples\COVID-19-CT\coronacases_002.nii.gz" -Destination "Examples\INPUT"
  Move-Item -Path "Examples\COVID-19-CT\coronacases_005.nii.gz" -Destination "Examples\INPUT"
  lung_extraction.ps1 .\Examples\INPUT .\Examples\LUNG

For lung extraction, a pre-trained UNet model was used. The model and the
code used to apply it belong to this_ repository. For more details, please
refers here_.

Labeling
--------

Once you have isolated the lung, you can run the actual segmentation.
In this case you have to create an other empty folder, wich will contains all the results.
As input the script requires the results of the lung extraction(in this case stored in LUNG folder).

Run from bash:

.. code-block:: bash

  mkdir ./Examples/OUTPUT
  labeling.sh ./Examples/LUNG ./Examples/OUTPUT

or its equivalent for powershell

  .. code-block:: powershell

    New-Item -Path "Examples" -Name "OUTPUT" -ItemType "directory"
    labeling.ps1 .\Examples\LUNG .\Examples\OUTPUT

This will run the segmentation by using the already estimated centroids. If you
want to use another set of centroids, simply provide as third arguments the path
of the file in which the set of centroids is saved

Train
-----

Even if with the script a set of pre-estimated centroids is provided, we also provide
a script to train another set of centroids. To perform the training simply organize
the scans resulting from the lung extraction into the same folder, this will be the
training set. Now simply the training script:

.. code-block:: bash

  python -m CTLungSeg.train --input='./Examples/LUNG' --output='./Examples/new_centroids.pkl.npy'

Once you have run this script, a brief recap of the training parameter will be
displayed :

.. code-block:: bash

  I m Loading...
  Loaded 20 files from ./Examples/LUNG
  *****Starting clustering*****
  Number of subsamples--> 100
  Total images --> 4000
  Centroid initialization technique-->KMEANS_RANDOM_CENTERS
  I m clustering...
  100%|█████████████████████████████████████████████████████████████████████████████████████| 100/100 [00:14<00:00,  2.86s/it]
  I m saving...
  [DONE]

All the images will be divided into N subsamples, and a K-means clustering is
performed for each subsample, after that a second clustering is performed in order
to refine the clustering and provide the set of centroids.
To control the parameters simply provides the following arguments when the script
is execute:

* init : centroid initialization algorithm: if 0 the centroids will be initialized randomly, if 1 the K-means++ center will be used.

* n : number of subsamples, as default as 100.

Once the training is complete, the centroid file will be stored in `.pkl.npy`
format.

.. note::

  please notice that this process may be time consuming and computational expansive

.. _SimpleITK: https://simpleitk.readthedocs.io/en/master/IO.html
.. _this: https://github.com/JoHof/lungmask
.. _here: https://eurradiolexp.springeropen.com/articles/10.1186/s41747-020-00173-2
COVID-19 Lung Segmentation
==========================

The SARS-CoV-2 virus has widely spread all over the world since the beginning of 2020.
This virus affects lung areas and causes respiratory illness. In this scenario is
highly desirable a method to identify in CT images the lung injuries caused by COVID-19.
The approach proposed here is based on colour quantization to identify the infection
regions inside the lung(Ground Glass Opacities, Consolidation and Pleural Effusion).

To achieve this purpose we have used the colour quantization approach to segment the
chest CT scans of patients affected by COVID-19. Use this technique as medical
image segmentation means to reduce the number of colours in the image to the number
of anatomical structures and tissue present in the anatomical region; in this
way we  assign to each kind of tissue a characteristic colour: so must exist a
relationship between the kind of tissue and the colour used to represent it.

For CT scan which is in greyscale, each colour is represented by a single value
given by the Hounsfield Units(HU): voxels colours are proportional to HU, which
are defined as a linear transformation of the linear attenuation coefficient.
HU normalize the coefficient of a particular tissue according to a reference one,
usually, water, as we can see in the equation below :

.. math::

  	HU = 1000\times\frac{\mu - \mu_{H_2 O}}{\mu_{H_2 O}}

In the end, each colour results proportional to the linear attenuation coefficient,
different from each tissue, so exist a relation between the GL and the tissue type
that makes these techniques available.

Colour quantization and the properties of digital images allow us to consider also
other properties of the image beside the single voxel intensity.
This purpose can be achieved by building a suitable colour space:

In digital image processing, images are represented with a 3D tensor, in which the
first two dimensions represent the height and width of the image and the last one
the number of channels. Grayscale images require only one channel, so each pixel
has a numeric value whose range may change according to the image format.
On the other hand colour images requires 3 channels, and the value of each channel
represent the level of the primary colour stored in this particular channel, so each
colour is represented by 3 different values, according to the Young model.
In this work, the different channels are used to takes into account different properties,
exploited by the application of different filters. This allows us to consider also
neighbouring voxels, which is really suitable for the segmentation since the
lesions areas involve many closest voxels, not only a single one. We have also
used these features to discriminate between other lung regions like bronchi by
exploit shape information.
The used image features are displayed in the figure below:

.. image:: images/Multi_Channel.png
   :height: 500px
   :width: 500px
   :scale: 100 %
   :alt:
   :align: left

Once we have built the colour space, we have to found the characteristic colour of
each tissue under study, which is represented by centroids in the colour space.
To perform this task and achieve the centroids estimation a simple -means
clustering was used.
K-means clustering requires prior knowledge about the number of clusters, which
in our case is given by the anatomical structure of the lung, so we can consider
a different cluster for each anatomical structure.
Once we have estimated the centroids for each tissue, we use that for the actual
segmentation by assign each voxel to the cluster of the closest centroids: in this
way the estimation step, that we will call "train", needs to be performed only once,
so can be time expansive since is not involved in the actual segmentation.

This package provides a set of already estimates centroids, together with scripts
to perform the actual segmentation. Also, a script to train your own set of centroids
is provided.
Snakemake
=========

Since the segmentation of several patients is time-consuming, we have provided a
snakemake pipeline to automate the process. This pipeline also allows to train
other set of centroids and use it for the segmentation. This file allows to
customize the usage of the hardware resources, like the number of threads and the
amount of memory.

As before, this examples will use he data previously downloaded from the public dataset

Segment Multiple Scan
---------------------

First of all, you have to create two folders:

  - INPUT : contains all and only the CT scans to segment
  - OUTPUT : empty folder, will contain the segmented scans as *nrrd*.

Now simply execute from command line

.. code-block:: bash

  snakemake --cores 1 --config input_path='./Examples/INPUT/' --output_path='./Examples/OUTPUT/'

.. note::

  It will create a folder named **LUNG** inside the INPUT, which
  contains the results of the lung extraction step.

Train a Centroid Set
--------------------

Prepare three folders:
  - INPUT: will contains all the scans to segment
  - OUTPUT: will contain the segmented scans
  - TRAIN: will contain all the scans of the training set.

Now run Snakemake with the following configuration parameters :

.. code-block:: bash

  snakemake --cores 1 --config input_path='./Examples/INPUT/' --output_path='./Examples/OUTPUT/'
  --train_path='./Examples/TRAIN/' --centroid_path='.Examples/centorid_set.pkl.npy'

This will train the centroid set and use them to segment the input scans.

.. note::

  This will create a folder named LUNG inside INPUT and TRAIN which
  contains the scans after lung extraction.

.. warning::

  The `TRAIN` folder cannot be the same of `INPUT`!

Configuration
-------------

We have provided a configuration file (config.yaml) which allows to manage the
resources and the path, which we usually provide from command-line.

**Threads**:

  - *threads_labelling* : Set the number of threads to use for the labelling process (default = 8);

  - *threads_lung_extraction* : Set the number of threads to use for the lung_extraction (default = 8);

  - *threads_train* : Set the number of threads to use for the training process (default = 8).

**Memory**:

  - memory_labelling : 8
  - memory_lung_extraction : 8
  - memory_train : 8

**Training Parameters**:

It is possible to specify the parameters for the training step:

  - n_subsamples : number of subsamples in which the slice of the training set  will be divided during the training;

  - centroid_initialization : technique to use for the initialization of the centroids during k-means (0 for random initialization, 1 for k-means++)
Installation
=================

Supported python versions :
|python version|

The full list of prerequisites is the following:

- numpy>=1.17
- opencv-python
- tdqm
- SimpleITK

And, for testing:

- PyTest>=3.0.7
- Hypothesis>=4.13.0

The lung extraction is performed by using pre-trained UNet, so please ensure to
have installed the lungmask_ package. For more information about how the network
is trained, please refers here_

Installation
------------

First of all, ensure to have the right python version and the package for the
lung extraction correctly installed

To install this package first of all you have to clone the repositories from GitHub:

.. code-block:: bash

  git clone https://github.com/RiccardoBiondi/segmentation

Now you can install tha package using pip

.. code-block:: bash

  pip install segmentation/

Testing
-------

Testing routines use pytest_ and hypothesis_ packages. please install
these packages to perform the test:

.. code-block:: bash

  pip install pytest>=3.0.7
  pip install hypothesis>=4.13.0

.. warning::
  Pytest versions above 6.1.2 are not available for python 3.5


All the full set of test is provided in the testing_ directory.
You can run the full list of test from the segmentation folder:

.. code-block:: bash

  cd segmentation
  python -m pytest


.. |python version| image:: https://img.shields.io/badge/python-3.5|3.6|3.7|3.8-blue.svg
.. _pytest: https://pypi.org/project/pytest/6.0.2/
.. _hypothesis: https://hypothesis.readthedocs.io/en/latest/
.. _testing: https://github.com/RiccardoBiondi/segmentation/tree/master/testing
.. _lungmask: https://github.com/JoHof/lungmask
.. _here: https://eurradiolexp.springeropen.com/articles/10.1186/s41747-020-00173-2
Modules
=======

Together with the scripts, a series of modules are provided. Each module
contains a series of functions for image processing which are used during the
script developing. The modules are the following, each of them provides a different
kind of functions.

Utils
-----

This modules provides all the functions to read and write images in a medical image
format like '.nrrd' or '.nifti'. All the formats supported by SimpleITK_ are allowed.

.. _SimpleITK: https://simpleitk.readthedocs.io/en/master/IO.html

.. automodule:: utils
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members: _read_dicom_series, _read_image
   :special-members:

Method
------

This module contains the implementation of all the filter used for the
processing of images inside the script. The functions are based on SimpleITK_
methods

.. _OpenCV: https://opencv.org/

.. automodule:: method
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:
   :special-members:

Segmentation
------------

This module contains the implementation of the functions used to perform the
tasks on each script.

.. automodule:: CTLungSeg.segmentation
   :members:
   :undoc-members:
   :show-inheritance:
   :private-members:
   :special-members:
.. CTLungSeg documentation master file, created by
   sphinx-quickstart on Wed Oct 21 12:55:39 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CTLungSeg's documentation!
=====================================

Automated Pipeline for the segmentation of Ground Glass Opacities on chest CT
scans of COVID-19 affected patients.

This package provides a fast way to isolate the lung region and identify ground glass
lesions on CT images of patients affected by COVID-19.
The segmentation approach is based on colour quantization,
performed using K-means clustering. This package provides a series of scripts to
isolate lung regions, pre-process the images, estimate K-means centroids and
labels the lung regions; together with methods to perform thresholding,
morphological and statistical operations on image series.

.. image :: ../../images/results.png

Usage Example
=============

Once you have installed you can directly start to segment the images.
Input CT scans must be in Hounsfield units(HU), gray-scale images are not allowed.
The input allowed formats are the ones supported by SimpleITK_ .

Example Data
------------

As Example data, we will use the ones of the public dataset  *COVID-19 CT Lung and Infection Segmentation Dataset*, published by Zenodo.
How to organize them depends on the purpose and will be explained for each tutorial.

Firstly, create the Examples folder, which will contain the dataset and the results.
After, you will download the .zip containing the data and unzip it.

So, run from bash:

.. code-block:: bash

  mkdir Examples
  wget https://zenodo.org/record/3757476/files/COVID-19-CT-Seg_20cases.zip -P ./Examples
  unzip ./Examples/COVID-19-CT-Seg_20cases.zip -d ./Examples/COVID-19-CT

or PowerShell:

.. code-block:: powershell

    New-Item  -Path . -Name "Examples" -ItemType "directory"
    Start-BitsTransfer -Source https://zenodo.org/record/3757476/files/COVID-19-CT-Seg_20cases.zip -Destination .\Examples\
    Expand-Archive -LiteralPath .\Examples\COVID-19-CT-Seg_20cases.zip -DestinationPath .\Examples\COVID-19-CT -Force

Single Patient Example
----------------------

To segment a single CT scan, run the following command from the bash or
PowerShell :

.. code-block:: bash

   python -m CTLungSeg --input='.Examples/COVID-19-CT/coronacases_002.nii.gz'  --output='./Examples/coronacases_002_label.nrrd'

Which takes as input the CT scan in each format supported by SimpleITK_. If the
input is a Dicom series, simply pass the path to the directory which contains
the series files, please ensure that in the folder there is only one series.

The output label will be saved as '.nrrd'.

Multiple Patient Example
------------------------

The segmentation of multiple patients can be time-consuming and tedious, so we have provided a series of scripts to automate this procedure.  Moreover, we have provided a snakemake pipeline.
To see their usage, please refer to script and snakemake contents.


.. _SimpleITK: https://simpleitk.org/


.. toctree::
  :maxdepth: 2
  :caption: Contents:

  installation
  theory
  modules
  script
  snakemake
  ./examples/examples
  references



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
References
----------

- Hofmanninger, J., Prayer, F., Pan, J. et al. Automatic lung segmentation in routine imaging is primarily a data diversity problem, not a methodology problem. Eur Radiol Exp 4, 50 (2020). https://doi.org/10.1186/s41747-020-00173-2

- Bradski, G. (2000). The OpenCV Library. Dr. Dobb&#x27;s Journal of Software Tools.

- Yaniv, Z., Lowekamp, B.C., Johnson, H.J. et al. SimpleITK Image-Analysis Notebooks: a Collaborative Environment for Education and Reproducible Research. J Digit Imaging 31, 290–303 (2018). https://doi.org/10.1007/s10278-017-0037-8

- Lowekamp Bradley, Chen David, Ibanez Luis, Blezek Daniel The Design of SimpleITK  Frontiers in Neuroinformatics 7, 45 (2013) https://www.frontiersin.org/article/10.3389/fninf.2013.00045

- Ma Jun, Ge Cheng, Wang Yixin, An Xingle, Gao Jiantao, Yu Ziqi, Zhang Minqing, Liu Xin, Deng Xueyuan, Cao Shucheng, Wei Hao, Mei Sen, Yang Xiaoyu, Nie Ziwei, Li Chen, Tian Lu, Zhu Yuntao, Zhu Qiongjie, Dong Guoqiang, & He Jian. (2020). COVID-19 CT Lung and Infection Segmentation Dataset (Verson 1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.3757476.
Examples
========

.. toctree::
  :maxdepth: 1

  ./pipeline_workflow.ipynb
  ./body_segmentation.ipynb
