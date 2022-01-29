# Orange hackathon

[![DOI](https://zenodo.org/badge/205136200.svg)](https://zenodo.org/badge/latestdoi/205136200)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu)
[![Research Software Directory](https://img.shields.io/badge/rsd-Research%20Software%20Directory-00a3e3.svg)](https://www.research-software.nl/software/orangehackathon)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=e-mental-health_orange-hackathon&metric=alert_status)](https://sonarcloud.io/dashboard?id=e-mental-health_orange-hackathon)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/71fd89c125d4482ab7b5d69f9c547073)](https://www.codacy.com/manual/eriktks/orange-hackathon?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=e-mental-health/orange-hackathon&amp;utm_campaign=Badge_Grade)
[![Maintainability](https://api.codeclimate.com/v1/badges/bba74d89a390004d9ed0/maintainability)](https://codeclimate.com/github/e-mental-health/orange-hackathon/maintainability)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/4691/badge)](https://bestpractices.coreinfrastructure.org/projects/4691)

The Orange Hackathon was a software development event in which researchers and students expanded the platform system [Orange](http://orange.biolab.si) with modules which can automatically analyze online therapeutic therapies

The hackathon was a part of the project [What Works When for Whom](https://www.esciencecenter.nl/project/what-works-when-for-whom) in which three partners collaborate: [Tactus](https://tactus.nl), [Psychology, Health and Technology](www.utwente.nl/en/bms/pht/) of the University of Twente and the [Netherlands eScience Center](esciencecenter.nl).

The contact person for this hackathon is Erik Tjong Kim Sang e.tjongkimsang@esciencecenter.nl


## Installation

In order to use the software developed in the hackathon, you need to have [Orange](http://orange.biolab.si) running at your computer. Next you need to install the hackathon software by following the instructions below or by watching our installation videos ([video 1](https://vimeo.com/485242301), [video 2](https://vimeo.com/487662686))

### Orange installation

We install Orange via Anaconda:

1. Download Anaconda from http://www.anaconda.com/distribution (choose the operation system of your computer: Windows, macOS or Linux)
2. When the download is complete: start the Anaconda prompt: under Windows: Anaconda prompt under Anaconda3
3. In the Anaconda3 window, type:
  * conda config --add channels conda-forge (followed by Enter)
  * conda install orange3
  * conda install orange3-text
  * on Windows: mkdir %userprofile%\orange
  * cd %userprofile%\\orange (or cd $HOME/orange on macOS or Linux)
  * git clone https://github.com/e-mental-health/orange-hackathon

**Starting Orange**

1. start the Anaconda prompt: under Windows: Anaconda prompt under Anaconda3
2. In the Anaconda3 window, type:
  * cd %userprofile%\\orange\\orange-hackathon (or cd ~/orange/orange-hackathon)
  * pip install .
  * python -m Orange.canvas

## Widgets
**Mail2Tsv**: 
Convert email files to .tsv files, for easy importing as a corpus.

_Instructions:_ select your input folder, and fill in the filter. All files from the subdirectories will be attempted to convert. By default the filter is set to all_documents/ and inbox/, resulting in only the conversion of files inside those subdirectories. 

Once a batch of emails is converted, it will be automatically send them out as a corpus, if any other node is connected. If emails have previously been converted, and the output path/file is still a valid .tsv file, it is enough to press the "Output file > corpus" button. This button only needs to be pressed in a fresh session, if the convert button has been pressed, the output file will automatically be converted.

* this widget requires Tkinter to be installed

**EmailSorter**:
Sorts emails converted by the Mail2Tsv widget in chronological order, based on the "date" column. 

_Instructions:_ No interaction is required, it will start automatically once it receives an input. Clicking on the widget will give the option to sort in ascending, or descending order.

## External data

The WRAD.Wt dictionary used by the module DAAP can be downloaded from

[https://github.com/DAAP/WRAD](https://github.com/DAAP/WRAD)

Note that we also use the [text module](https://github.com/biolab/orange3-text) of Orange3

## Examples

### Small pipeline

<img src="https://raw.githubusercontent.com/e-mental-health/orange-hackathon/master/images/orange-small.jpg" width="60%">

### Medium pipeline

<img src="https://raw.githubusercontent.com/e-mental-health/orange-hackathon/master/images/orange-medium.jpg" width="100%">

### Tactus pipeline

<img src="https://raw.githubusercontent.com/e-mental-health/orange-hackathon/master/images/tactus-pipeline.jpg" width="100%">

The Tactus pipeline is used for filtering and visualizing the therapy mails from the organization [Tactus](https://tactus.nl) in the project [What Works When for Whom](https://www.esciencecenter.nl/project/what-works-when-for-whom). It contains six modules:

1. Tactus Mail Loader: read all mails of one therapy session (xml data)
2. Sort Emails: sort mails chronologically
3. Mark duplicates: mark all phrases of 20 tokens or longer that appeared in an earlier mail
4. Remove Marked Text: remove all marked text from the mails
5. LIWC: perform LIWC analysis of the mail text: count relevant tokens
6. Line Plot: make a graph of the counted tokens
