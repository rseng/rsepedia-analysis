[![Build Status](https://app.travis-ci.com/EISy-as-Py/hardy.svg?branch=master)](https://app.travis-ci.com/EISy-as-Py/hardy)
[![Coverage Status](https://coveralls.io/repos/github/EISy-as-Py/hardy/badge.svg?branch=master)](https://coveralls.io/github/EISy-as-Py/hardy?branch=master)
[![Documentation Status](https://readthedocs.org/projects/hardy/badge/?version=latest)](https://hardy.readthedocs.io/en/latest/?badge=latest)
[![Anaconda-Server Badge](https://anaconda.org/pozzorg/hardy/badges/platforms.svg)](https://anaconda.org/pozzorg/hardy)
[![Anaconda-Server Badge](https://anaconda.org/pozzorg/hardy/badges/installer/conda.svg)](https://conda.anaconda.org/pozzorg)
[![Anaconda-Server Badge](https://anaconda.org/pozzorg/hardy/badges/license.svg)](https://anaconda.org/pozzorg/hardy)


[![DOI](https://joss.theoj.org/papers/10.21105/joss.03829/status.svg)](https://doi.org/10.21105/joss.03829)

<img src=https://github.com/EISy-as-Py/hardy/blob/master/doc/images/EIS_Formats.PNG width=400 p align="right">

# Project HARDy

 _"HARDy: Handling Arbitrary Recognition of Data in python"_
A package to assist in discovery, research, and classification of YOUR data, no matter who you are!

## Project Objective

Numerical and visual transformation of experimental data to improve its classification and cataloging

This project was part of DIRECT Capstone Project at University of Washington and was presented at the showcase, follow this
<a href=https://prezi.com/view/5ugf5HyDxZevQlOHmuyO/>link </a>  for the presentation

## Requirements:
Package HARDy has following main dependencies:
1. Python = 3.7
2. Tensorflow = 2.0

The detailed list of dependencies is reflected in the <a href=https://github.com/EISy-as-Py/hardy/blob/master/environment.yml><code>environment.yml</code></a> file

## Installation:
The package HARDy can be installed using following command:

<code>conda install -c pozzorg hardy </code>

Alternatively, you can also install it using the GitHub repository in following steps:

*Please note that currently v1.0 is the most stable release

1. In your terminal, run <code>git clone https://github.com/EISy-as-Py/hardy.git</code>
2. Change the directory to hardy root directory, by running <code>cd hardy</code> 
3. Run <code>git checkout v1.0</code>
4. Run <code>python setup.py install</code>
5. To check installation run, <code>python -c "import hardy"</code> in your terminal

For other methods of installation like using environment file and installation using pip, please visit <a href=https://hardy.readthedocs.io/en/latest/installation.html>Installation<a> page.

## Usage:

HARDy uses Keras for training Convolutional Neural Network & Keras-tuner for the hyperparameter optimization. The flow of information is shown in image below:

<p align="center"><img src="https://github.com/EISy-as-Py/hardy/blob/master/doc/images/HARDy_diagram.png" width=700 alt="information flow of how the package works"/></p>

An example jupyter notebook to run HARDy using single script is available at this link <a href=https://github.com/EISy-as-Py/hardy/blob/master/doc/examples/example_HARDy_script.ipynb>Example Notebook</code></a>

To perform various transformations, training Neural Network and Hyperparameter Optimization, Hardy utilizes following <code>.yaml</code> configuration files:

* <a href=https://github.com/EISy-as-Py/hardy/blob/master/hardy/arbitrage/README.md>tform_config.yaml</a>
* <a href=https://github.com/EISy-as-Py/hardy/blob/master/hardy/recognition/README.md>cnn_config.yaml</a>
* <a href=https://github.com/EISy-as-Py/hardy/blob/master/hardy/recognition/README.md>tuner_config.yaml</a>

The instructions for modifying or writing your own configuration file can be accessed by clicking on the configuration files listed above.

The notebooks and documentations can also be accessed at this link <a href=https://hardy.readthedocs.io/en/latest/index.html>Documentations</a>

## Visualization
 In order to increase the density of data presented to the convolutional neural network and add a visual transformation of the data, we adopted a new plotting technique that takes advantage of how images are read by computers. Using color images, we were able to encode the experimental data in the pixel value, using different series per each image channel. The results are data- dense images, which are also pretty to look at.

 <p align="center"><img src="https://github.com/EISy-as-Py/hardy/blob/master/doc/images/data_visualization.PNG" width=700 alt=" details on the proposed visual transformation to increased the images data density"/></p>


## Mission:
We have been commissioned by Professor Lilo Pozzo to create a new tool for research and discovery, For her lab and for high throughput researchers everywhere.
Our vision of the final product:
 * A package which can approach any large, labeled dataset (such as those familiar to High Throughput Screening (HTS) researchers).
 * Perform a (procedurally generated and data-guided) wide array of transformations on the data to produce completely novel ways of examining the data, maybe not Human-Readable but in a certainly machine-readable format.
 * Train "A Machine Learning Algorithm" (We currently focus on Visual-Processing CNNs but are open to anything!) to classify the existing labled data based on each of the aforementioned transformations.
 * Report back to the user:
    * Which versions of the Model/Algorithm worked best?
    * Which transformations appeared the most useful? (AKA were used across many of the most successful models)
    * What Data "Fingerprints" should we pay the most attention to?
 * Present a User Interface, to allow non-programmers to interact with and use the chosen classifier(s?) in their work.

 ## Use Cases:
 The package is designed to deal with a diverse set of labeled data. These are some of the use cases we see benefitting from using the _HARDy_ package.

 <p align="center"><img src="https://github.com/EISy-as-Py/hardy/blob/master/doc/images/use_cases.PNG" width=500 alt="possible use cases for the HARDy package"/></p>


 ## Modules Overview:
 * __handling.py__         :  Functions related to configuration, importing/exporting, and other sorts of back-end useful tasks.
 * __arbitrage.py__        :  Data Pre-Analysis, Transformations, and other preparation to be fed into the learning algorithm.
 * __recognition.py__      :  Setup, training and testing of single convolutional neural network (CNN) or hyperparameters optimization for CNNs.
 * __data_reporting.py__   :  Output and reporting of any/all results. Tabular summary of runs, visual performance comparison, as well as parallel coordinate plots and feature maps

 ## Community Guidlines:
 We welcome the members of open-source community to extend the functionalities of HARDy, submit feature requests and report bugs.
 
 ### Feature Request:
 If you would like to suggest a feature or start a discussion on possible extension of HARDy, please feel free to <a href="https://github.com/EISy-as-Py/hardy/issues/new/choose">raise an issue</a>
 
 ### Bug Report:
 If you would like to report a bug, please follow <a href="https://github.com/EISy-as-Py/hardy/issues/new/choose">this link</a>
 
 ### Contributions:
 If you would to contribute to HARDy, you can fork the repository, add your contribution and generate a pull request. The complete guide to make contributions can be found at this <a href="https://github.com/EISy-as-Py/hardy/blob/master/CONTRIBUTIONS.md">link</a>

 ## Acknowledgment

 Maria Politi acknowledges support from the National Science Foundation through NSF-CBET grant 1917340
# Welcome to the guide for making contributions to HARDy

Please report all the bugs in HARDy modules and mistakes in documentation by creating an <a href=https://github.com/EISy-as-Py/hardy/issues/new/choose>issue</a>. You can also request features using the same link. Alternatively, if you want to make contribution to hardy, please feel free to create a pull request. The procedure to make contributions is as follows:

### Step 1. Setting up the repository

The first step requires you to fork the repository on your GitHub using the button on top right corner of HARDy page. The button is highlighted in image below:


### Step 2. Cloning the repository

The next step is to clone the repository by running the following command in your temrinal/command prompt:
```
git clone https://github.com/your-usename/hardy
```

### Step 3. Making changes

The repository is then ready for the changes to be made. After making the changes, make sure you add, commit and push to your forked branch by running following commands:

```
git add changed_file
git commit -m "A decriptive message outlining the changes made into the repository"
git push
```

### Step 4. Generating pull request

The final step is to generate a pull request. The instructions for opening a pull request are outlined in this <a href="https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request">tutorial</a>


## Things to note
To achieve the purpose of continous integration, <code>hardy</code> utilizes <code>travis CI</code>. travis CI is configured to perform PEP8 compliance test (flake8) along with running unittests. Any failure in compliance would result in build error or build failure on travis CI.

To run these test on your local machine, following code snippets can be used in your terminal/command prompt:

For flake8,
```
conda install flake8
cd hardy
flake8
```

For unittest,
```
cd hardy
python -m unittest discover
```### Project HARDy: Main Package
Congrats, you're in the codebase for project HARDy. 
You can see that so far (edited 2020-05-05), we don't have a "Main" function for the package, everything is contianed in these submodule folders.


### Module Descriptions:
In general, the folders here are submodules that each have their own __init__.py file and can be called independantly, but they also follow a workflow into eachother!

The __Handling__ module should call no others, but is available to be called on as needed to load, save, and pre-process data.

The __Arbitrage__ module has the arbitrage.py function in it (which is mostly the class wrapper function to store and execute data transformations, and then to produce images as needed!
 * Inside this module is also transformations.py, which contains plans and code detailing all of the 1D, 2D, and more complex transformations that the arbitrage.py class will call!
 
The __classifier__ module (To be renamed __Recognition???__) has all of our CNN construction code using the Keras package, as well as wrapping functions to call the above modules and execute data operations

## Instructions for using the tform_config file

An example transformation configuration file is shown below:

<img src="https://github.com/EISy-as-Py/hardy/blob/master/doc/images/Quickstart_TransformConfig.PNG" width=700 />

The transformation configuration relies on list of transformation a user intends to perform. The name of each transformation follows the variable naming rules for python. The names for transformation must be listed under <code>tform_command_list</code>.

The rules for plotting must be defined in <code>tform_command_dict</code>. The header of each entry in the dictionary corresponds to the transformation name which should be same as entered in the <code>tform_command_list</code>. The operations performed on the data are defined as arrays under this entry.

As many as <b>six</b> definitions can be entered under transformation command of dictionary. Each command follows the structure of \[column_number, mathematical operation, plotting_value].

The <code>column_number</code> corresponds to the column number according to the data in csv file. <code>Mathematical operation</code> is the operation that needs to be performed on this column and <code>plotting_value</code> corresponds to the color and orientation of the plot in final image that is to be read by machine learning algorithm.

The scheme for <code>plotting_values</code> is as follow:

```
- 0: Red on x-axis
- 1: Green on x-axis
- 2: Blue on x-axis
- 3: Red on y-axis
- 4: Green on y-axis
- 5: Blue on y-axis
```

Currently supported mathematical operations are as follows:

```
- raw: returns raw data without performing any operation
- exp: exponential
- nlog: natural log
- log10: logarithm tranformation with base 10
- reciprocal: reciprocal
- cumsum: cumulative sum
- derivative_1d: Differential with respect to 1 dimension
- derivative_2d: 2-D differentiation
- power: can be used for array multiplication or to take user defined power for array 
```
## Instructions for using Convolutional Neural Network Configuration files

Each configuration file in the recognition folder provides the example Hyperparameter space over which Neural Network model is built. The configuration file can be placed anywhere in the system and relative path must be passed as argument in <code>run_hardy.hardy_multi_transform</code> module

### Configuration file for CNN

This file provides the input information along with the Hyperparameter space for Neural Network.

<img src="https://github.com/EISy-as-Py/hardy/blob/master/doc/images/Quickstart_cnn_config.PNG" width=500 p align="right" />

* __cnn_config.yaml__

_A configuration file which contains the hyperparameters to use in the single convolutional neural network.
The configuration file is easy to fill out and interact with._


__Note__: Make sure that the hyperparameters found in the config. file are also used in the cnn model


Currently supported keys for the <code>cnn_config.yaml</code> includes:

```
-> kernel_size

-> epochs

-> activation

-> input_shape

-> filter_size

-> num_classes

-> learning_rate

-> patience
```

__Note__: All this information must be entered to successfully execute the Machine Learning Step. The detailed information about the options for keys can be found in the config 
<a href=https://github.com/EISy-as-Py/hardy/blob/master/hardy/recognition/cnn_config.yaml>file</a> itself.

<hr>

### Configuration file for Hyperparameter Optimization

* __tuner_config.yaml__
    
A configuration file containing the hyperparamter search space for the tuning step. This should substitute the single cnn model. 

The first part deals with defining the tuner run:   

For definition of tuner run, following keys are currently supported:

<img src="https://github.com/EISy-as-Py/hardy/blob/master/doc/images/Quickstart__tuner_config_run.PNG" width=450 p align="right" />


```
-> num_classes

-> epochs

-> patience

-> input_shape

-> max_trials

-> exec_per_trial

-> search_function
```

<img src="https://github.com/EISy-as-Py/hardy/blob/master/doc/images/Quickstart__tuner_config_space.PNG" width=400 p align="right" />

The second section, deals with the actual hyperparameter search space to use in the tuning operation:

For hyperparameters search space, following keys are currently supported:
```
-> layers

-> filters

-> kernel_size

-> activation

-> pooling

-> optimizer

-> learning_rate
```

__Note__: All this information must be entered to successfully execute the Hyperparameter Step. The detailed information about the options for keys can be found in the config <a href=https://github.com/EISy-as-Py/hardy/blob/master/hardy/recognition/tuner_config.yaml>file</a> itself.
---
name: Discussion Starter
about: Use this to reach out to us for discussion
title: "[Discussion] : "
labels: question
assignees: ''

---

Hello!
We're always glad to start a discussion for your project or implementation,
Just let us know a bit about your Use-Case and we'll see what advice we can use to get you started.

** Your Use-Case **
Say Hi! Who are you, and what are you trying to do? We love to hear about new projects.


** Your Data **
<u> File Type and organization:</u> [csv? excel?], [header? uniform shape?]
<u> Data quantity: </u> [csv? excel?], [header? uniform shape?]
<u> Data amount: </u>[Existing amount?], [Frequency of new data?]
<u> Data Organization: </u> [Labeled by Filenames? Sorted by folder? Other?]

** Your Computing **
<u> RAM and Processing:</u> How fast can your computer run?
<u> Time constratints:</u> How patient will you be for a run?
<u> Access to other Computing:</u> If we can support Cluster/ Supercomputer script options, would that help?
---
name: Feature request
about: Suggest an idea or improvement
title: ''
labels: enhancement
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
---
name: Bug report
about: Create a report to help us improve
title: "[BUG] : "
labels: bug
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
 - OS: [e.g. iOS, Windows, Linux]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
# Functional Specifications
# And User-Case Stories
### __Updated 2020-04-10__
--------------------------------------------------------

## Meet the Team:
* __Alice (she/her)__, a Data Scientist post-doc who has years of data and research. She wants to set up the project, but then focus on other work like writing papers and grants. That's what PostDocs do, right?

* __Bonnie (she/her)__, a Quality Control Manager for a mid-sized startup. She's setting up a QC monitoring lab at the new factory, and will have a massive influx of data that she needs to manage in order to ensure constant quality of the Batteries we produce. That data includes incoming Material Assays (digestion UV-Vis, XRD, ShortAngle Scattering, Particle Size), as well as electrochemical testing (CV cycling, EIS, and even performance cycling). So much Data! And it all needs to be tracked and tested for pass/fail.

* __Claude (they/them)__, a ChemE Undergrad working in the lab three days a week. Pretty competent, but doesn't have a background in either programming or  the research done in its lab. They doesn't want to screw up, because he needs Grad School recs from Alice and Dave, and is afraid of Bob.

* __Dave (Your Highness)__, a stern old Professor who's lab they all work in. Your Highness is not a micro-manager, but simply expects results to be delivered on time. Dave is skeptical of the error and quality of high-throughput data, so needs to be convinced that all the results have well defined error explainations. The last programming Dave did was in Visual Basic in 1998.

---------------------------------------------------------

## USER STORIES:

### Alice:
 * She will use her old data and __Import & Classify__ tools to make a classifyer and help categorize incoming data.
 * She can then use our __CNN Modeling__ tools to bootstrap a Classifyer model, to confidently identify and classify her data into the expected models.
 * She will hand off this model to the rest of the lab, and get back to papers and research. If she gets her hands on more data, the model will already be set up to accept it and improve.

### Bob: 
 * He will set up his lab experiment to export data to a single folder on his computer or a server he has access to. [Future Work: Including binary _Export Tools_ to read directly from the machine's format such as biologic .mpt, etc]
 * Using a few files as examples, he can use __SQL Configuration__ tools to set up a SQL database for his project, and also save __Import Configuration__ files so that data importing is always as easy and reliable as possible!
 * On a daily basis, he can use our __Import Tools__ to grab all related files and process their raw data (F-series, not yet T-series).
 * The __SQL Export__ tools will put all of the raw data, experimental metadata, and processed results to the Lab's SQL Database. 
 * He will use our __Image Export__ tools to save images to be processed, AND/OR he can directly send these images (perhaps without actually making the image files???) to the __Classification Model__.
 * The __Classification Model__ will make several decisions, and report back each of the following statistics AS WELL AS the Confidence and Model-Training statistics that accompany the decision!
	* __Smooth, Noisy-But-Ok, Noisy-And-Fail__: Bob needs to know whether there was something wrong in the data, because some experiments may need to be repeated if the data is bad. There may even be important trends in why some data is noisier than others.
	* __Frequency Settings (Max too low? Min too High? Min too low?)__: Depending on the set-up, Bob maybe needs to change the settings on his equipment to make sure the results are reasonable. 
	* __Dominant Circuit Behavior__: What type of circuit-fit does this data represent? 
	* __MANY OTHER POSISBILITIES__: Bob will work with Alice to continue training the model, so that eventually it can do almost all the processing for them, and they can focus the NEXT Project: What to do with all that Wonderful SQL Meta-Data! But that's a project for another Semester...

### Claude: 
 * They will be the most hands-on with the equiment and data collection, so need to have a reliable and user friendly interface to make sure that everything is running smoothly. Our __HTML User Interface__ will eventually be JUST what they need!
 * When they get to the lab, Bob will hand them a list of what should be run today. 
 * First, they will run iniital checks on the equipment using two or three well known test circuits (an RC, an LRC, and a Randle's setup). When they put those data sets into the __User Interface__, it will use the __Classification Model__ behind the scenes to confirm that these are confidently the correct results, that noise is low, and that no frequency parameters are unusual. 
 * They will save these Calibration results in a separate "Calibration Tracking" __SQL Database__ using the __GUI__, and then switch the program to load/store data in the main __Experiment SQL Database__. 
 * They will run experiments as normal, and do their homework like a good little student... If it's not automated yet, they will export the data and feed it into the __GUI__ every 30 minutes or so, and it will quickly give feedback on whether things are going according to plan. If not, they can make adjustments as recommended and/or text Bob for advice. 

### Dave: 
 * Dave's new best friend is the __HTML User Interface__, because unlike Grad-Students and Post-Docs it doesn't worry about his feelings and so he can see for himself all of the data with proper, accurate confidence intervals. 
 * Of course, Dave doesn't actually do that very often unless there is a paper coming soon. Dave really uses the __GUI Exporter__ to make beautiful plots and to check the results of Dave's own hypotheses about the experiments. 
 * This way, when that crazy black-box-machine-learning-voodoo-witchcraft tells Bob or Alice that they have discovered a fantastic breakthrough, Dave can explain that he knew that would happen all along.
# Component Specification:
List of MODULES, each containing Component Spec "Cards". 
If there's a clever way of formatting this, we should consider it!

Each "Module" contains:

__module.py__
 * Purpose and Description
 * Timeline + Milestones (?)
 * Current Status
 * LIST of CARDS / FUNCTIONS


Each "Card" contains:

__func_name():__
 * Purpose and general description
 * INPUT: description of imputs- types and sources
 * OUTPUT: description of returns
 * ACTIONS: Any non-obvious actions the function can take other than IN/OUT (example - option to saving files)
 * NOTES:
 -------------------------------------------------------------------------------
 
 ## handling.py
 The zeroth-level functions, related to importing and setting up data. The twin goals are to be fast, but also as broad and thoughtful about the data we may import as possible. 
 We may crib from work done in eisy/data_managment.py since we contributed to both project.
 __Note:__ Some of this may have Options to User-Interface, which is great but should DEFAULT be off, so we don't deal with it during Automated Processes!
 (Similarly, travis doesn't like User-Interfacing...)
 __Timeline + Milestones__:
  * 2020-04-21: Basic 2D, 3D importer that works with the simple files we're using 
  * 2020-04-28: Importer that works with each of the Raw Data types that we're using as part of our part-1 program. __HAND OFF!__
  * 2020-05-12: Debugging and Tests that cover some basic and logical use-cases (header, Too-Many-Columns, etc)
  * 2020-06-09: "Smart-ish" handling of the use cases mentioned. This should be the final __HAND OFF__
  
 __Current Status__:
  * (2020-04-14)
  * Just creating files and setup, no progress yet
  * Planning Functions, in compontent spec document (here!)
  * Considering how much to Frankenstien from prior work.
  
 __Module List__:
 #### import_2d():
  * Basic file to import a csv that has two columns (with or without column names?) - Essentially pd.read_csv, but with error handling to work better in our program
  * INPUT: filepath/name, error-handling instructions (interact=False, return2=True)
  * OUTPUT: x_data, y_data  <-- Two (Tuples? Series? Discuss) 1D datasets 
  * ACTIONS: Checks that the file is actualy 2D data - if not, either asks for guidance, errors, or returns the first two columns depending on instructions
 
 #### import_3d():
  * Basic file to import a csv that has three columns (with or without column names?) - Essentially pd.read_csv, but with error handling to work better in our program
  * INPUT: filepath/name, error-handling instructions (interact=False, return3=True)
  * OUTPUT: x_data, y_data, z_data  <-- Three (Tuples? Series? Discuss) 1D datasets 
  * ACTIONS: Checks that the file is actualy 3D data - if not, either asks for guidance, errors, or returns the first two columns depending on instructions
  
 #### get_classification_files():
  * Our project will expect a specific though simple folder structure: main __"/train/"__ data folder will have (2? more?) classification folders with the classifier names as the folder names - inside will be all the raw data in (csv?) format. 
  * INPUT: Base folder containing classified data folders
  * OUTPUT: Dictionary, where keys are the classifications (folder names) and each have a (tuple/list) of file names to load
  * ACTIONS: Checks the file and data format. May try and load one/more files? 
  * NOTE: Should we estimate size of all files and eventually try and estimate program-time? that may be useful...
  
 
 ## arbitrage.py
 The first-level functions, which will take input data (either 2D, 3D, or eventually nD...), and perform transformations to generate the full set of data-columns that we will test against!
 This package will contain several 'sections':
  * __Transformation Functions:__ This is the mathematical side, and to start out we will be able to perform a variety of 1D or 2D transformations such as Log, Inverse, accumulate, Integrate, derrive, etc.
  * __Complex Transforms__: Some data transformations are combinations of the ones above (you can integrate AFTER you log-ify, for instance)
  * __Perform Transformation__: Wrapper function(s) to actually /do/ the transformations. This MAY be smart-ish and perform only specified ones, perhaps from a CONFIG or LIST? (That List and Config is an idea we'll talk about in __yNot__ probably!)
  * __Association Functions:__ This is optional, but once you've done the transformations (or before?) you can check how well the data CORRELATES/ASSOCIATES, both to itself and to a "Standard" dataset AKA starting with linear data and transforming it! This may give us "SCORES" for each transformation, which we can use to prioritize!
 
 __Note (again):__ Some of this may have Options to User-Interface, which is great but should DEFAULT be off, so we don't deal with it during Automated Processes!
 (Similarly, travis doesn't like User-Interfacing...)
 
 __Timeline + Milestones__:
  * 2020-04-21: List of the high-priority functions and Simple-Transformations, with progress and timeline to get them all done soon.
  * 2020-04-28: Passing Tests and can __HAND OFF__ to the classifier - a DataFrame of "all" the transformed data columns. Recieve Handoff from handling, and begin to Integrate
  * 2020-05-12: Complex Transforms - consider what other things we may want, and discuss feedback with group
  * 2020-06-09: Make Decision on Association functions and __HAND OFF__ if so. Otherwise, simply focus on new group priorities
  * 2020-06-23: IF yNot function is doing Configuration ideas, make Stretch-Goal learning gameplan... TBD...
 
 __Current Status__:
  * (2020-04-14)
  * Just creating files and setup, no progress yet
  * Planning Functions, in compontent spec document (here!)
  * Considering how much to Frankenstien from prior work. Configuration?
 
 __Module List__:
 
 ### *SECTION: Basic 1D Transformations*
 #### transform_log():
  * INPUT: 1D data array with NO NEGATIVE VALUES
  * OUTPUT: Logrythmic transform of that data 
  * Note: Consider shifting or abs() for negative data? no? simply don't call log transform for negative data?
 #### transform_reciprocal():
  * INPUT: 1D data array - Limits tbd?
  * OUTPUT: all values inverted (1/x)
  * Note: twice should return itself!
 #### transform_cumsum():
  * INPUT: 1D data array - 
  * OUTPUT: cumulative sum of data (aka integrated with unit-steps)
  * Note: Not in high priority list?
 #### transform_1d_derrivative():
  * INPUT: 1D data array
  * OUTPUT: the step-by-step delta (Note: copy last delta to retain length?)
  * Note: Also not in high-priority list? Should also be able to complete the loop w/ cumsum.
 #### transform_exp():
  * INPUT: 1D data array - Limits?
  * OUTPUT: e^x of each datapoint. 
  * Note: may be redundant in general with log? should be able to complete that loop!
 #### transform_0to1():
  * INPUT: 1D data array, data handling case instructions
  * OUTPUT: that array shifted and scaled to the 0-to-1 basis (by FIRST shifting to min=0, THEN scaling to max=1)
  * Note: option to leave data alone if min is already 0-to-1, or if max After Shift is already 0-to-1 (case: data begins 0.2-0.4, can either scale 0-1 or leave as is!)

 ### *SECTION: Basic 2D Transformations*
 #### transform_2D_int():
  * INPUT: 2 equal size 1D arrays Y, X, to be integrated (Y)dx - [Optional offset value? to use as the Plus-C]
  * OUTPUT: The integral of Y dx (BOX? Trapz?) - offset if instructed to.
  * Note: Error handling? What if not Sorted/Linear in X? Should we sort by X first? (Or, what if reciprocal data ie CV Sweeps?)
 #### transform_2D_der():
  * INPUT: 2 equal size 1D arrays Y, X - Sorted? in either? 
  * OUTPUT: the Single-point derivatives dY/dX, also the offset so you /could/ integrate it back again!
  * Note: could use some sort of average or smoothing to reduce noise? However that would be LOSSY DATA PRACTICE!
 #### transform_prod():
  * INPUT: 2 equal size 1D arrays X, Y , [Optional Power arguments? or do those in the 1D cases and use as inputs?]
  * OUTPUT: Product of each x*y, maybe with power-math included (options)
  
### *SECTION: More Complex Transformations*
 #### transform_fourier_wavelets():
  * Ok so this is the only High-Priority one that I'm genuinely concerned with... while you "CAN" try to do a transform on a whole dataset, that gets noisy and lossy. 
    What I want to investigate is "Wavelet Filtering" Fourier transform, which we learned about at a Data Sci seminar last quarter?
    (Or otherwise, there's a whole realm of Signal-transforming science, I can research that...)
  * INPUT:  2 equal size 1D arrays X, Y - Sorted in X?? (Frequency range parameters? or is that the X-size?)
  * OUTPUT: 2D? output matrix or Meshgrid - in X-Freq space (for each wavelet size, return the match(-1 to 1?)*amplitude at each X?)
  * NOTE: This will have to be a group discussion- we need TEST DATA that should work in this space, and then we can report that back!
 #### multi_transform():
  * Wrapping function, to perform multiple transformations all together... Not sure which of these may be useful but I can see possible value in knowing the integral of a log function, for example. 
  * INPUT: X, [Y if 2D], Multiple transforms to perform... Is this what classes are for??
  * OUTPUT: Data output from the final transform listed. 
  * NOTE: This is low-priority, and should only be done if we convince ourselves that it's useful... RELATED, if we get the "Smart" learning functionality, maybe we can combine things this way
 #### transform___():
  *
  * INPUT: 
  * OUTPUT:
  
### *SECTION: Perform Transformation*
 #### get_xy():
  * Uses handling.py to load a file, check that we have xy data, and do a quick analysis on it! 
  * INPUT: file name
  * OUTPUT: each of 1D arrays X and Y, plus messages OR list of "Approved" transformations??
 #### perform_transform():
  * Wrapping function, to take some input data and return a dataframe with every transform that we want to use:
  * INPUT: 1D arrays X, Y, some sort of list or guiance for what transforms to do
  * OUTPUT: Pandas dataframe of X and Y with all of their transforms as requested.
  * NOTE: There might be a better way to do this - 
 #### generate_linear_transforms():
  * Creates a "sample" linear 2D dataset, possibly following range instructions, and performs all(?) transforms on the data
  * INPUT: [ALL OPTIONAL? Default 0-1 and X=Y], [List of transforms to do?? default True to perform "ALL"]
  * OUTPUT: pandas dataframe with all X-transformationa and all Y-transformations - (Maybe in Standard names? Maybe Not?)
  * NOTE: This could get messy - And maybe use the same wrapping function above to perform the listed tranforms?
  
 ### *SECTION: Association Functions*
  * NOTE: Not really planning these yet, that will be scoped or descoped based on Team Update by __2020-04-28__
 #### setup_correlation_matrix():
  * Sets up the 2D matrix of "scores" to judge the correlations 
 #### correlate_to_transforms():
  * For a given transform, determine (?) how correlated the data is to all other columns in the dataframe
 #### correlate_to_linear():
  * Compares the given transform to the linear_transform() function transforms. Any that correlate are probably good ideas?!?
 #### correlate_to_null():
  * Maybe compares the given transform to a flat line of low but nonzero values (all value = 0.1)? 
  * Not sure what's the best way to do this... 
  * What I'm TRYING to do is to identify/flag "BORING" data, which are probably BAD transforms to use?
 #### grade_all_transforms():
  * Wrapping function. Given a "fully" transformed dataset (or generate it here?), run all correlations and use some fancy math or grading (we generate?) to give each column a __SCORE__
  * The __SCORE__ should reflect how "Interesting" we think the data is (which is a topic for discussion, but all zeros is not interesting)
  * INPUT: dataframe ready to be "graded", OR give an XY dataset and we'll call the transform functions on it
  * OUTPUT: dictionary of Key,Values where each Key is a transform (or column key), and each Value is a "grade" to estimate how interesting we think the data may be
  * NOTE: This is SUPER arbitrary and is the 'creative' part of the STRETCH-GOALS of the project.s!
 #### grade_all_files(): 
  * Load all the files in a list, perform transformations, and grade. Hopefully this is a fast function so you can do a large list of files.
  * Then combine all the grades to get an average idea of what transforms we consider "good" 
  * INPUT: list of files
  * OUTPUT: dictionary of key, values as before, where results are averaged across dataset
  * ACTIONS: Optionally, save results as a report csv (or append to existing csv report??), to track grades over time
 #### load_transform_results():
  * IF we've run this program before, we should have a file that has previous "grades", this time based on the model training
  * If one transform shows up in a lot of the best-performing models, we should bump it to the top of the Transform To-Do List!
 #### combine_grades_results():
  * Somehow we should combine the 'grades' with the 'prior results' to get our new guess for what are the 'best' transformations to try
  * this will allow our model to do the highest-profile transfomations first and hopefully get good results in fewer attempts.
  
## recognition.py
This is the CNN, Baby! Finally we get to have the Machine do the Learning, instead of doing it ourselves!
Setup and Configure a CNN model, run optimization, and stash the result! 
Since we're stepping through this document and the package in a rough workflow order, this may also call functions from the previous two functions. 
It occurs to me that we haven't turned things into images yet? That can go wherever we want, but it may as well go here because the Image size and data might be an important factor in running the CNN that we need control over!

 __Timeline + Milestones__:
  * 2020-04-21: Initialize CNN, possibly with setup from preveious project EISy... Also set up some sort of image creator to feed it. __NEED: Raw Data input to use__
  * 2020-04-28: Manually able to generate an image set and run the CNN with a given list of metrics. __DEMONSTRATE For Group__
  * 2020-05-12: Integrate with hand-offs from above, given files and list of Transforms, can loop through multiple models and generate Report.
  * 2020-06-09: Improve "strategy" for looping, and discuss calculation time, duration, and next-steps. __HAND OFF__ to data_reporting.py

__Current Status__:
  * (2020-04-14)
  * Just creating files and setup, no progress yet
  * Planning Functions, in compontent spec document (here!)
  * Considering how much to Frankenstien from prior work.

__Module List__:
 #### func_name():
 * Purpose and general description
 * INPUT: description of imputs- types and sources
 * OUTPUT: description of returns
 * ACTIONS: Any non-obvious actions the function can take other than IN/OUT (example - option to saving files)
 * NOTES:
 
 
## data_reporting.py
  This will handle the INPUT/OUTPUT and Reporting for the results of the model... Group vision here is TBD, but I have a sketch of the various outputs that we might want to report:
   * __Run_list.csv:__ A record with all of the manual parameters of each time we run the program. Should record meta info about how the run was set up, how many/how large files, and how long the run took (useful so we can get a feel for the computing power needed). Of course, also report some sort of 'success' metrics like the best performing result (and settings?) or the performance curve (how good were the best ones all together?). Append every new run with the run ID to point to it's result file...
   * __Report_yymmdd_ID.csv:__ A report generated for each run. */maybe/* only generate if we pass a certian performance standard? (aka the run "passes"). Maybe a large CSV including all of the parameters of the best 10-or-so model runs, (or all above a certain %-spec?)
   * __Best_Transformations.csv:__ A Meta-report that keeps records of any Transformations that occur often across the best performing models for the given dataset. Maybe simply append list with every run (or every *good* run), or perhaps keep this as a "LeaderBoard", which can be Replaced as new tests change the algorythms' opinion of the "best" transformations...
   
 __Timeline + Milestones__:
  * 2020-04-21: 
  * 2020-04-28: 
  * 2020-05-12: 
  * 2020-06-09: 

__Current Status__:
  * (2020-04-14)
  * Just creating files and setup, no progress yet
  * Planning Functions, in compontent spec document (here!)
  * Considering how much to Frankenstien from prior work.

__Module List__:
 #### func_name():
 * Purpose and general description
 * INPUT: description of imputs- types and sources
 * OUTPUT: description of returns
 * ACTIONS: Any non-obvious actions the function can take other than IN/OUT (example - option to saving files)
 * NOTES:
 
 
 ## y_not.py *...or yNot.py ?*
  Mostly TBD, we included this as a catchall or pun, but have no specific plans.:
  
 __Timeline + Milestones__:
  * 2020-04-21: 
  * 2020-04-28: 
  * 2020-05-12: 
  * 2020-06-09: 

__Current Status__:
  * (2020-04-14)
  * Just creating files and setup, no progress yet
  * Planning Functions, in compontent spec document (here!)
  * Considering how much to Frankenstien from prior work.

__Module List__:
 #### func_name():
 * Purpose and general description
 * INPUT: description of imputs- types and sources
 * OUTPUT: description of returns
 * ACTIONS: Any non-obvious actions the function can take other than IN/OUT (example - option to saving files)
 * NOTES:
 
  ## hardy.py *...or yNot.py ?*
  Core Python Notebook? Do we need/want a Core notebook to wrap all of the functions in?
  
 __Timeline + Milestones__:
  * 2020-04-21: 
  * 2020-04-28: 
  * 2020-05-12: 
  * 2020-06-09: 

__Current Status__:
  * (2020-04-14)
  * Just creating files and setup, no progress yet
  * Planning Functions, in compontent spec document (here!)
  * Considering how much to Frankenstien from prior work.

__Module List__:
 #### func_name():
 * Purpose and general description
 * INPUT: description of imputs- types and sources
 * OUTPUT: description of returns
 * ACTIONS: Any non-obvious actions the function can take other than IN/OUT (example - option to saving files)
 * NOTES:
 
---
title: 'HARDy: Handling Arbitrary Recognition of Data in Python'
tags:
    - Feature Engineering
    - Kernel methods
    - Machine Learning
    - Python
authors:
    - name: Maria Politi^[co-first author][^1]
      email: politim@uw.edu
      orcid: 0000-0002-5815-3371
      affiliation: 1
    - name: Abdul Moeez^[co-first author]
      email: amoeez@uw.edu
      orcid: 0000-0002-9582-0372
      affiliation: 2
    - name: David Beck
      email: dacb@uw.edu
      orcid: 0000-0002-5371-7035
      affiliation: "1,3"
    - name: Stuart Adler
      email: stuadler@uw.edu
      affiliation: 1
    - name: Lilo Pozzo
      email: dpozzo@uw.edu
      orcid: 0000-0001-7104-9061
      affiliation: 1
affiliations:
    - name: University of Washington, Department of Chemical Engineering, Seattle, WA, USA
      index: 1
    - name: University of Washington, Department of Materials Science and Engineering, Seattle, WA, USA
      index: 2
    - name: eScience Institute, University of Washington, Seattle, WA, USA
      index: 3


date: 10 April December 2021
bibliography: paper.bib

---

`HARDy` is a Python-based package that helps evaluate differences in data through feature engineering coupled with kernel methods. The package provides an extension to machine learning by adding layers of feature transformation and representation. The workflow of the package is as follows:

- _Configuration_: Sets attribute for user-defined transformations, machine learning hyperparameters or hyperparameter space
- _Handling_: Imports pre-labelled data from `.csv` files and loads into the catalogue. Later the data will be split into training and testing sets
- _Arbitrage_: Applies user defined numerical and visual transformations to all the data loaded.
- _Recognition_: Machine Learning module that applies user defined hyperparameter search space for training and evaluation of model
- _Data-Reporting_: Imports result of machine learning models and reports it into dataframes and plots

# Statement of Need

High Throughput Experimentation (HTE) and High Throughput Testing (HTT) have exponentially increased the volume of experimental data available to scientists. One of the major bottlenecks in their implementation is the data analysis. The need for autonomous binning and classification has seen an increase in the employment of machine learning approaches in discovery of catalysts, energy materials and process parameters for design of experiment [@williams2019enabling; @becker2019low]. However, these solutions rely on specific sets of hyperparameters for their machine learning models to achieve the desired purpose. Furthermore, numerical data from experimental characterization of materials carries diversity in both features and magnitude. These features are traditionally extracted using deterministic models based on empirical relationships between variables of the process under investigation. As an example, X-ray diffraction (XRD) data is easier to characterize in linear form as compared to small angle X-ray scattering data, which requires transformation of axis to log-log scale.

One of the most widely applied strategy to enhance the performance of machine learning model is Combined Automatic Machine Learning (AutoML) for CASH (Combined Alogrithm Selection and Hyperparameter Optimization) [@hutter2019automated]. However, these packages are only limited to hyper-parameter tuning and data features remain untouched. To improve the effectiveness of machine learning models, some of the popular feature engineering strategies used for simple numerical data include binning, binarization, normalization, Box-Cox Transformations and Quantile Sketch Array (QSA) [@zheng2018feature; @nargesian2017learning]. Moreover, Deep Feature Synthesis has also shown promising results. Here features are generated from relational databases by performing multi-layer mathematical transformation operations [@kanter2015deep].

`HARDy` presents an infrastructure which aids in the identification of the best combination of numerical and visual transformations to improve data classification through Convolutional Neural Networks (CNN). `HARDy` exploits the difference between human-readable images of experimental data (i.e., Cartesian representation) and computer-readable plots, which maximizes the data density presented to an algorithm and reduce superfluous information. `HARDy` uses configuration files, fed to the open-source package `KerasTuner`, removing the need for the user to manually generate unique parameters combinations for each neural network model to be investigated.



# Description and Use Case

The Python-based package `HARDy` is a modularly structured package which classifies data using 2D convolutional neural networks. A schematic for the package can be found in figure 1.

![Workflow schematics for HARDy. The data files, left-most column, are subject to numerical and visual transformations, according to the rules outlined in the user defined configuration files. The images are then fed into either a CNN or a tuner, for which the hyperparameter space is controlled through another configuration file. Finally, each transformation produces a report comprising of the best trained model file, the log of training session and the model validation result.](./images/HARDy_diagram.png)

The package was tested on a set of simulated Small Angle Scattering (SAS) data to be classified into four different particle models: spherical, ellipsoidal, cylindrical and core-shell spherical. A total of ten thousand files were generated for each model. The data was generated using \textit{sasmodels}. The geometrical and physical parameters used to obtain each spectrum were taken from a published work discussing a similar classification task [@ArchibaldRichardK2020Caas]. The name of each SAS model was used as label for the data, allowing for further validation of the test set results. These models were selected as they present similar parameters and data features, which at times make it challenging to distinguish between them.

![Summary table of few transformations visualized using cartesian coordinate and RGB representation along with the respective fittings. The left-most column shows the transformations applied to the data: the scattering vector _q_ and scattering intensity _I(q)_. The final model accuracies are also provided. The original data label corresponds to the icon on the left-most graph, whereas the icon under each fit correspond to the model prediction.](./images/transformation_run_example.png)

First, the pre-labelled data was loaded. A subset of the files, three thousand files in total, was identified as the testing set. All the ML models initialized in the same code run were validated using the same testing set. A user-provided list of transformations, inputted through a configuration file, was then applied to the data. Different trials can be specified, so that multiple sets of transformations can be investigated. Both Cartesian and RGB plots representations were compared. The latter visualization option was obtained by encoding the data into the pixel values of each channel composing a color image, for a total of six-channels available (i.e., 3 RGB channels in horizontal/vertical orthogonal directions).

The data was then fed into a convolutional neural network, whose hyperparameters and structure were defined using another configuration file. Alternatively, it is also possible to train multiple classifiers for a single transformation trial through the use of a tuner, by instead providing a hyperparameter space and a search method. The classification results, as well as the best performing trained model were saved for each transformation run. The package also allows to visually compare, through parallel coordinates plots (see documentation), the performance of each transformation. Figure 2 shows a summary of few runs comparing the two visualization strategies and their best performing model accuracies.

Comprehensive results for all transformations tested are available in the documentation. It can be noticed that data representation using Cartesian coordinate plots yielded a higher number of instances in which the accuracy of the trained machine learning model was ~25\%. This value corresponds to machine learning model's inability to recognize differences in a four-class classification task. On the other hand, the RGB plots show, on average, higher accuracy for the same combinations of numerical transformations. To further validate the results, mathematical fitting was performed on a test set using the SASmodels package. The fitting was based on probabilities determined by the ML model for each label. In scenarios where the output probability was below 70\%, the data was also fitted using the second highest possible SAS model.

![Model fitting examples for various validation files used under the following transformation: _log(q)_ & _derivative I(q)_. The model was trained using RGB images and yielded the highest classification accuracy of 90.1%. The bottom labels used correspond to the predicted model that was assigned to the data by the model, whereas the top label correspond to the original model used to generate the data. ](./images/fitting_example.png)

The average chi-square parameter of the fitted data was determined to be 7.5. Approximately 11 \% of the data had a probability lower than 70\%. In all cases,as seen in figure 3, if the neural network was not able to correctly label the data with the highest probability label, the second highest probability label was the correct one.

In conclusion, `HARDy` can significantly improve data classification so that automatic data fitting and modeling can be executed without human intervention and without compromising on reliability. We also note that data representation for computer-classification tasks may not follow human intuition and/or standard conventions. `HARDy` serves a key role in the optimization of visual data representations for CNN classification tasks. Finally, the flexibility of `HARDy` allows for deployment of the task on a supercomputing cluster system, possibly removing the limitations given by the high computational power required to run these ML algorithms. All configuration files and scripts used to run the example presented in this paper can be found in the package documentation.


# Acknowledgements
This project was supported by: National Science Foundation through NSF-CBET grant no. 1917340, the Data Intensive Research Enabling Clean Technology (DIRECT) National Science Foundation (NSF) National Research Traineeship (DGE-1633216), the State of Washington through the University of Washington (UW) Clean Energy Institute and the UW eScience Institute.
This work was facilitated through the use of advanced computational, storage, and networking infrastructure provided by the Hyak supercomputer system and funded by the STF at the University of Washington.

# References

[^1]: corresponding author
Guide to numerical transformations
==================================

Following numerical transformations are currently defined in the
:code:`HARDy`

.. automodule:: hardy.arbitrage.transformations
   :members:

Defining numerical transformations
----------------------------------

User defined numerical transformations can be integrated in the
:code:`HARDy` by integrating the transformation definition
inside the :code:`hardy.arbitrage.transformations` module. The
example transformation definition is shown below::

    def transformation_function(args):
        y = f(x)

Most of the transformations defined in :code:`HARDy` are one
dimensional transformations i-e they require only one arguments.
If the user want to include more than one argument and metadata
(exponents), refer to function definition of :code:`power` or
:code:`derivative_2d` defined in :code:`hardy.arbitrage.transformations`
module. If the metadata is different than the length of metadata in
:code:`power`, function :code:`apply_tform` in
:code:`hardy.arbitrage.arbitrage` needs to be modified
accordingly as well.


Data\ Reporting
=============================

Submodules
----------

reporting module
--------------------------------------

.. automodule:: hardy.data_reporting.reporting
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: hardy.data_reporting
   :members:
   :undoc-members:
   :show-inheritance:
.. Example of running HARDy

Examples
========

.. toctree::
    :maxdepth: 1
    :glob:

    examples/example_sas_data.ipynb
    examples/example_HARDy_script.ipynb
    examples/How_to_write_Configuration_files.ipynb
    examples/How_to_make_predictions_using_trained_model.ipynb.. Hardy Information

hardy package
=============

The ``HARDy`` package is composed of following modules::

   - Handling: Module to transform data
   - Arbitrage: Module to create data set ready to be fed in Machine Learning Model
   - Recognition: Machine Learning Engine
   - Data Reporting: Post Classification analysis tool

Subpackages
-----------

.. toctree::

   hardy.handling
   hardy.arbitrage
   hardy.recognition
   hardy.data_reporting


Module contents
---------------

.. automodule:: hardy
   :members:
   :undoc-members:
   :show-inheritance:
.. Hardy Information

Module Information
==================

The ``HARDy`` package is composed of following modules::

   - Handling: Module to transform data
   - Arbitrage: Module to create data set ready to be fed in Machine Learning Model
   - Recognition: Machine Learning Engine
   - Data Reporting: Post Classification analysis tool

Modules
-------

.. toctree::

   hardy.handling
   hardy.arbitrage
   hardy.recognition
   hardy.data_reporting


Arbitrage
=======================

Submodules
----------

arbitrage module
--------------------------------

.. automodule:: hardy.arbitrage.arbitrage
   :members:
   :undoc-members:
   :show-inheritance:

transformations module
--------------------------------------

.. automodule:: hardy.arbitrage.transformations
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: hardy.arbitrage
   :members:
   :undoc-members:
   :show-inheritance:
.. hardy documentation master file, created by
   sphinx-quickstart on Fri Apr 17 16:41:00 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HARDy
=====
:code:`HARDy` Handling, Arbitrage, Recognition of Data using python is a
package to assist in discovery, research and classification of data.

:code:`HARDy` utilizes feature engineering and kernal methods simultaneously
to improve the accuracy of convolutional neural network. :code:`HARDy` enables
the user to search over given range of hyperparameters to achieve the best
prediction score for a machine learning model. The summary of results can be
printed out using single argument custom built plotly function.

Installation
------------
The easiest way to install :code:`HARDy` is using :code:`conda`::
 
   conda install -c pozzorg hardy 

For detailed installation instructions and other methods of installation, visit
`Installation Page <https://hardy.readthedocs.io/en/latest/installation.html>`_

Dependencies
~~~~~~~~~~~~

:code:`HARDy` requires following modules::

- Python (3.7)
- Tensorflow (2.0)
- Keras Tuner (1.0.0)

Complete list of dependencies can be found in the :code:`environment` file

Examples and Documentations!
=================================

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   gstarted
   advanced
   examples
   modules
   faqs

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Handling
======================

Submodules
----------

handling module
------------------------------

.. automodule:: hardy.handling.handling
   :members:
   :undoc-members:
   :show-inheritance:

pre\_processing module
-------------------------------------

.. automodule:: hardy.handling.pre_processing
   :members:
   :undoc-members:
   :show-inheritance:

to\_catalogue module
-----------------------------------

.. automodule:: hardy.handling.to_catalogue
   :members:
   :undoc-members:
   :show-inheritance:


visualization module
-----------------------------------

.. automodule:: hardy.handling.visualization
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: hardy.handling
   :members:
   :undoc-members:
   :show-inheritance:
Advanced Functionalities
========================
The guide provided on `Getting Started` page uses the wrapper
function :code:`hardy_main` which sequentially executes
:code:`data_wrapper`, :code:`classifier_wrapper` and report generation.
The structure of :code:`HARDy` package is shown in the image below:

.. image:: ../joss-paper/images/HARDy_diagram.png
    :width: 400
    :align: center
    :alt: image explaining the structure of HARDy

:code:`data_wrapper` parses through the transformation configuration
file and loads it into the environment. The :code:`data wrapper` also 
loads the :code:`.csv` files into the environment. The data wrapper
then applies the transformations, outlined in the configuration file,
to the data files. The :code:`data_wrapper` stores the tranformed
information into the :code:`.pkl` file or into :code:`.png` files
depending on the user defined arguments.

To save time :code:`data_wrapper` is capable of parallelizing the 
transformations process. The parallelization is controlled through the 
:code:`n_threads` parameter. By default, :code:`n_threads` is set a 1.

:code:`classifier_wrapper` loads the data generated by :code:`data_wrapper`. 
By default, :code:`classifier_wrapper` deletes the :code:`.pkl` generated
by the :code:`data_wrapper`. The loaded data is then ran through the 
convolutional neural network (CNN) or hyperparameter tuning sessions
depending upon the specified inputs from the user.

Each tuning session is then reported into :code:`data_reporting` module.
The tuned/trained model is then validated againt :code:`num_test_files_class`.
Then the report folder having :code:`project_name` in :code:`data_path`
is created. The report folder contains folders for each transformation run.
These transformation folders contains model validation results, best tuned model
for a particular transformation and model hyperparameter details.

Using Advanced Tools & Customization
------------------------------------

.. toctree::
    :maxdepth: 1
    :glob:

    transformations
    examples/rgb_cart.ipynb
    examples/using_data_reporting.ipynb
    examples/numpy_directory.ipynb
    examples/kfoldvalidation.ipynb
    examples/feature_mapping.ipynbGetting Started
===============
:code:`HARDy` works by leveraging the information density in data
representation which makes it easier for machine to make inference
as compared to the human readable representation. This is achieved
through numerical and visual transformation. :code:`HARDy` in the
first stage achieves numerical transformation and then it transforms
it into RGB representation as shown in the image below:

.. image:: images/hardy_gstarted.png
    :width: 400
    :align: center
    :alt: image explaining the numerical and visual transformation

:code:`Hardy`, by default, is configured to take minimal inputs
from the user and perform numerical and visual transformations 
on its own. The numerical transformations follow rules defined
by the user in a :code:`.yaml` configuration files. The user can
perform either hyperparamter search to evaluate best hyperparameters
or run a simple convolutional neural network (CNN).
The hyperparameter space for both tuning session
or :code:`CNN` can be defined in a :code:`.yaml` configuration
file. The guide to write configuration files is available at
`Guide to write configuration files 
<https://hardy.readthedocs.io/en/latest/examples/How_to_write_Configuration_files.html>`_

Data Preparation
----------------
:code:`HARDy` is configured to input :code:`.csv` files only. Before
starting your :code:`HARDy` run, make sure following conditions are met:
    * Data files are only in :code:`.csv` format.
    * The :code:`.csv` files must have a header of same length.
    * The training data files must also include the labels in  their filenames. The labels should be unique and must not overlap. For example, label `core-shell sphere` for scattering model overlaps with the label `sphere`. This overlap should be avoided.
    * The data files must have same number of data rows.

The wrapper function, :code:`run_hardy`, takes care of all the numerical
and visual transformations along with hyperparameter tuning and CNN runs.
The example script to run HARDy is as follows:

Importing HARDy library
-----------------------

The following code snippet imports :code:`hardy` into your respective environment::

    import hardy.run_hardy as run

Defining path variables
-----------------------
Defining the path to :code:`.csv` files::

    raw_data_path = 'path/to/raw/data/'

Defining the path to numerical configuration file::

    tform_config_path = './hardy/arbitrage/tform_config.yaml'

Defining the path to tuner or CNN configuration::

    classifier_config_path = './hardy/recognition/'

Setting up parameters
---------------------
The wrapper function of :code:`HARDy (hardy-main)` is capable to intake a large number
of user defined inputs. However, some important parameters that should be defined before
the runs are described here:

    * :code:`scale`: The RGB images generated by :code:`HARDy` are NxN dimensions, where N
      is the length of data rows in the csv files. It was observed that if the data
      has 500 rows, the RGB image would have 500x500 dimension. This exponentially
      increases the demand for RAM. To overcome this issue, :code:`scale` argument
      has been introduced to scale down the images. It uses :code:`scipy.misc.imresize`
      to scale down the image. :code:`scale` represents a fractional value for example
      0.1, 0.2.. 1.0
    * :code:`target_size`: target size is the product of :code:`scale` and :code:`N`.
      It is tuple of shape 1x2. For example, (100, 100)
    * :code:`num_test_files_class`: this is an integer value representing the number
      of files per category/class. This number of test of files will be seperated out
      from the training and only be utilized during final validation of model. The files
      are randomly selected. To control which data files are used, main function is capable
      of taking :code:`seed` as the input.
    * :code:`classes`: list of classes/categories in which data is categorized. Remembering
      the order of classes is important since this order is utilized for classfiying the
      new data inputs for a trained.
    * :code:`iterator_mode`: Since :code:`Keras` is capable of using either numpy array
      representation of data as well as image data, this string value controls which
      iterator to be used. If the value is :code:`arrays`, the RGB data representation
      will be in numpy array format and :code:`Keras` will use numpy iterator. If value
      other than the :code:`arrays`, PNG images will be generated, saved on disk and will
      be fed into the :code:`Keras` using a directory iterator.
    * :code:`n_threads`: This parameter controls the parallel processing of transformations.
      Since transformation configuration file may include saveral transformations, this parameter
      can assign each transformation to each n_thread available. A value of 6 would means
      each of the 6 threads will be utilized to run 6 transformations in parallel. The
      :code:`n_thread` parameter doesn't control the threads used for traning and tuning of
      neural networks. The neural network multiprocessing is either controlled through CUDA
      or default tensforflow operations.

Other arguments that can be given as input of :code:`hardy_main` function are described below

Executing hardy run
-------------------
The following code starts the numerical and visual transformations along with the
hyperparameter tuning session::

    run.hardy_main(raw_data_path, tform_config_path, classifier_config_path, batch_size=64,
    scale=0.2, num_test_files_class=750, target_size=(100, 100), iterator_mode='arrays',
    classifier='tuner', n_threads=1, classes=['class_1', 'class_2', 'class_3'],
    project_name='my_project_name')

The following arguments are acceptable in the :code:`hardy_main()` function:

    * ``raw_data_path``: data_path for the .csv files or images
    * ``tform_config_path``: path for transformation configuration files (.yaml)
    * ``classifier_config_path``: path for hyperparameter search (.yaml)
    * ``batch_size``: batch size for splitting of training and testing of data in machine learning model
    * ``scale``: the scale to which plots are reduce
    * ``num_test_files_class``: The number of test files per class. These files would be reserved for final testing of machine learning model
    * ``target_size``: number of data points in the csv files or dimension of images
    * ``iterator_mode``: if "arrays", the data is fed into machine learning model in array structure. For other values, images files are saved first in .png format and then fed into machine learning model through directory iterators.
    * ``classifier``: tuner or cnn model. Tuner means hyperparameter search while other options execute pre-defined convolutional neural network.
    * ``n_thread``: number of threads used for parallel transformation of data
    * ``classes``: labels or categories in data. If .csv files are used, the label must be present in the filename. If images are used, the images must be contained in respective folders
    * ``project_name``: name for the project. Folder with same name will be created in the raw_data_path containing all the results for the run
    * ``plot_format``: format of the plot to be used for training and testing of data. RGBrgb corresponds to usage of RGB images while any other argument will use cartesian coordinate system.
    * ``skiprows``: Used to skip the metadata contained in the csv files. It must be of same length for all classes.
    * ``split``: The fraction of data used for training and testing of machine learning model. This is different from num_test_files_class since the later one is never fed into machine learning model until the best hyperparameter search is done.
    * ``seed``: the seed used for random-selection of num_test_files_class
    * ``k_fold``: Boolean value indicating whether k-fold validation need to be performed or not
    * ``k``: value indicating how many k-folds need to be performed

Evaluating Results
------------------

After the :code:`HARDy` run is complete, the results for each transformation can be
found under the path::

    raw_data_path/project_name/transformation_name

The results include best trained model, evaluation result for best model and hyperparameter
configuration for best model. The reports can be analyzed through :code:`data_reporting`
module in `HARDy`. Its usage is described in :code:`Advanced Functionalities` section.


.. Installation doc

Installation
============

Installation using Conda
------------------------
The easiest way to install :code:`HARDy` is using :code:`conda`::
 
   conda install -c pozzorg hardy

Installation using Git
----------------------
:code:`HARDy` can also be installed using Git. Currently version 1.0
is the most stable version. To install version 1.0, follow the following
steps:

* In your terminal, run::

    git clone https://github.com/EISy-as-Py/hardy.git

* Change the directory to hardy root directory, by running::

    cd hardy
    
* Run::

    git checkout v1.0
    
* Run::

    python setup.py install

To check installation run the following command in your terminal::
    
    python -c "import hardy"


Installation using ``evironment.yml`` (Recommended)
---------------------------------------------------
To avoid installing each dependency one by one, we recommend using
environment.yml provided in the github repository. To install the
environment run the following code in your terminal::

    conda install --name hardy --file environment.yml

To proceed with the installation, the hardy environment needs to be
acitivated through::

    conda activate hardy
    
The final step is run the installation command for :code:`HARDy`::

    conda install -c pozzorg hardy

To check installation run the following command in your terminal::
    
    python -c "import hardy"

Installation using pip
----------------------
:code:`HARDy` is also configured to be installed using pip. Currently
version 1.0 is the most stable version. To install version 1.0, run the
following commands in the terminal::

    git clone https://github.com/EISy-as-Py/hardy.git
    cd HARDy
    git checkout v1.0
    pip install .

To check installation run the following command in your terminal::
    
    python -c "import hardy"


.. toctree::
    :maxdepth: 1
    :glob:Recognition
=========================

Submodules
----------

cnn module
----------------------------

.. automodule:: hardy.recognition.cnn
   :members:
   :undoc-members:
   :show-inheritance:

tuner module
------------------------------

.. automodule:: hardy.recognition.tuner
   :members:
   :undoc-members:
   :show-inheritance:


Module contents
---------------

.. automodule:: hardy.recognition
   :members:
   :undoc-members:
   :show-inheritance:
