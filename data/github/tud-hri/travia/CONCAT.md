# Travia: a TRAffic VIsualization and Annotation tool
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03607/status.svg)](https://doi.org/10.21105/joss.03607)

In recent years many open datasets have been published that contain vehicle trajectories of human drivers on public roads. Such datasets can be useful for all
kinds of research targeting traffic, e.g. traffic-flow studies and the development of autonomous vehicles. These datasets contain very large amounts of valuable
data, however, they often focus on collecting data in one country and one situation, e.g. highway traffic in Germany. To broaden the scope of a research
project, it can be valuable to include data from different datasets. This is often difficult because traffic data comes in *.csv files, often without a way to
visualize and annotate the data for a quick start of your project. This is why Travia was created to combine the potential of different datasets by bringing 
them
together in one tool for visualization and annotation.

Datasets from four different projects can currently be visualized using Travia. The NGSIM project collected data using cameras on tall buildings next to roads
in the USA. The PNeuma project captured all traffic during the morning rush hour in the business district of Athens, for 5 days. To maximize precision and
coverage this was done using drones. Finally, the HighD and ExiD projects also used drones to record traffic at different locations on German highways.

TraViA was built as a broad basis for visualization and manual annotation, it can easily be extended to incorporate the specific needs of your research project.
Three examples of such specific implementations are included in Travia to illustrate the possibilities and provide a starting point for other developers. The
first example considers automatically detecting and annotating specific situations in the HighD datasets. The second example is the plotting of heatmap
overlays, to aid with the design of reward functions for autonomous driving. Finally, it includes an example of plotting vehicle signals on an annotated
selection of data.

Travia was tested with Python 3.8 on Windows and Ubuntu.

## Installation

To install TraViA, simply clone the repository and install all requirements from the `requirements.txt` file using pip (run `pip install -r requirements.txt` in
the terminal from the travia folder). That's all!

To use TraViA, you'll need some traffic data, this data is not included in the repository. Follow the steps below to get started and test if TraViA is
functioning properly.

## Getting started, visualize data in 3 steps

To get started, you only need to execute the next 3 steps.

1. acquire data
2. extract in the correct folder
4. run

### Acquire data

The data to visualize is not included in Travia, you need to acquire the data and a valid license to use it yourself. All implemented datasets are free for use
for scientific purposes. In this example we will acquire the NGSim data set, but if you already have a copy of one of the supported datasets you can skip this
step. To acquire a copy of NGSim data, open a browser and navigate
to [this](https://data.transportation.gov/Automobiles/Next-Generation-Simulation-NGSIM-Vehicle-Trajector/8ect-6jqj)
page. Scroll down to find the attachments section under "about this dataset". Locate the attachment called "US-101-LosAngeles-CA.zip" and click it to start the
download. This contains data recorded on the US-101 highway.

### Extract in the right folder

To use the data, it has to be placed in the `travia/data` folder to be found automatically when loading the dataset. The data folder has three sub-folders, one
for each data source. For our example, we need the NGSim sub-folder, but if you have data from another project please see the `README.txt` file in the
corresponding sub-folder for instructions. For now, navigate to your freshly downloaded data and extract the contents of `US-101-LosAngeles-CA.zip`
to `travia/data/NGSim`. This should create a new folder, `travia/data/NGSim/US-101-LosAngeles-CA`, which contains more zip files. Extract the contents of
`aerial-ortho-photos.zip` and `vehicle-trajectory-data.zip` to this folder. This will create more sub-folders, a representation of the correct end-result can be
found in the `README.txt` file located in `travia/data/NGSim`, please verify that your folder structure is correct.

### Run

Now you're ready to run `visualize.py`, this is the main run script in TraViA. A dialog will show asking you to select which dataset you want to load, 
please select the NGSim dataset that you just downloaded. After this, the data will be loaded from the *.csv file and converted into a dataset object. 
Loading a new NGSim dataset for the first time will take a while since the data will be smoothed first. This is done to calculate headings for all vehicles 
because these are not included in the dataset. For more information, check out the section on "Particularities with the data". This smoothing is also needed 
for PNeuma datasets.

Once done, the TraViA user interface should be displayed. A screenshot of how it should look can be found below. This screenshot shows a portion of the
HighD dataset, so don't worry if the traffic looks different. To verify TraViA is working correctly, a video is included with the
repository (`./paper/image/US101_0805_0820.mp4`). You can compare the visualization of the NGSim data from the example with this video to verify TraViA is
working correctly.

![TraViA user interface](./paper/images/screenshot.png "the Travia user interface")

## How to use TraViA

TraViA was built with two use cases in mind, the first being to provide a quick visualization and annotation opportunity for researchers working with open
traffic data. For this use case, all interaction between the user and TraViA can take place in the GUI. An example of such use would be a researcher that
performed numerical analysis on the data from the HighD dataset in a custom-made script. If this analysis shows something happening at frame x for vehicle y
that is hard to explain numerically, TraViA can be used to a) investigate what is happening in a visualization and b) annotate this event for future use. For
such use, the user interface of TraViA is explained below in the subsection user interface.

The second potential use-case is for researchers to extend TraViA for their use, for example, to visualise other (closed access) datasets or to implement
Autonomous vehicle controllers that operate in the traffic data environment. To aid in this use case, a high-level overview of the architecture of TraViA is
provided in the subsection High-Level Software Architecture.

TraViA was **not** designed to be imported as a module or library for use in other software projects. For that reason, TraViA has no `setup.py` file, is not
available on PyPi and there is no API provided for software interfacing with TraViA.

### User Interface

The Travia user interface consists of four parts (see the figure above for more information). From top to bottom, you'll see a vehicle information pane,
displaying all data of the selected vehicle, a view of the road and traffic, the annotation overview and controls, and finally the timeline controls.

The traffic view can be panned by dragging with the mouse, zoomed by using the scroll wheel, and rotated using the knob on the top left. To display information
on a single vehicle you can select it by clicking it. The same holds for annotations, they are represented by bars just above the timeline. By clicking a bar,
the annotation is selected and its information is displayed. It is not possible to selected multiple annotations or vehicles at once.

The "start annotation" button will open a new annotation that starts at the current frame. After this, it is possible to edit the start and end frames, ego
vehicle ID (the vehicle of interest), and add a note. The "stop annotation" button will then save the annotation.

The time can be controlled by using the play/pause button in the lower part of the interface. The time can be fast-forwarded using the >> button and reversed 
with the << button. The > and < buttons will execute a single frame step. The rec button will start a video recording of the visualization, this recording 
can be stopped by clicking the pause button and will then be exported to the user's video folder. A single frame can be saved as an image through the view 
menu. The buttons "Create Plots" and "Create Heatmap" are used for example function and only work with HighD datasets. Please see the section HighD example 
tools below for more information.

### High-Level Software Architecture

![A conceptual UML class diagram showing the main classes in TraViA.
The design with base classes and child classes for specific data sources make it easy to add data from other sources. The separate dataset and visualization
master objects allow for maintaining the source format in low-level storage while translating to a generic high-level implementation in the visualization
master.](./paper/images/UML_class_diagram.png)

There are three main objects in TraViA used for visualizing the data, a conceptual UML class diagram can be found in the figure above. The highest level 
object is the Graphical User Interface (GUI). This inherits from the Qt QMainWindow object and handles all user interaction. The state of the GUI is 
periodically updated by a visualizationMaster object. This is the object that translates the specific format of the dataset to store the data for every 
vehicle in a generic vehicle object. This way of working ensures that the dataset itself is used in its original format such that all original documentation 
can be used and other scripts based on that documentation will still work. But TraViA uses the same generic `Vehicle` object for all datasets such that tools 
developed in TraViA can easily be used on multiple datasets.
  
The visualization master object also keeps track of the clock to allow real-time visualization. For every data source, there is a specific visualizationMaster 
that inherits from a visualizationMaster base class. The data itself is stored in a Dataset object, again there is a specific class for each source, inheriting 
from a base class. As said, the Dataset objects preserve the format of the data source by storing the contents of `*.csv` data files in pandas DataFrame 
objects. Only for the pNEUMA data, the original format of storing a vehicle per row is reshaped to storing a frame per row to obtain a DataFrame with equal 
length rows. All data stored in the Dataset objects are converted to SI units, global coordinates are converted to suitable local coordinate systems in meters. 
Metadata on the dataset is stored as `Dataset` object attributes.

## Datasets and how to get them

As said before, Travia supports data from four different projects:

- [HighD](https://www.highd-dataset.com/)
- [ExiD](https://www.exid-dataset.com/)
- [NGSim](https://ops.fhwa.dot.gov/trafficanalysistools/ngsim.htm)
- [PNeuma](https://open-traffic.epfl.ch/)

Follow these links to the websites of the projects. There, you can find more information on how to obtain a copy of the data, and a license to use it.

### Particularities with the data

There are some particularities with these datasets that had to be accounted for in Travia. These particularities and the solutions implements in Travia are
listed here. Please have a look at the code to get a better idea of the specific implementation. If you find any other issues or extensions of these issues
that are not accounted for. Please submit a bug report and/or a pull request containing a fix on the GitHub repository.

**Heading detection with Kalman smoother**

Only the Exid dataset provides headings for the vehicles, all of the other datasets do not. For the HighD data, this is no problem since the highways are 
straight, and the headings can be assumed to be either 0 or pi at these high velocities. But for the NGSim and PNeuma datasets, headings have to be 
estimated to visualize the data. This is done automatically using an unscented Kalman smoother in `processing/kalmansmoothing.py`. The smoother uses a 
bicycle model to estimate the dynamics of vehicles. These smoothers were manually tuned to provide results that are good enough for visualization purposes, 
but only that. If you want to use the smoothed data for anything else than visualization, please verify that the smoother did not introduce problems for 
your approach.

**missing vehicle sizes in PNeuma**

The PNeuma dataset does not include vehicle sizes for the individual vehicles, only a very specific vehicle type. To account for this, generic vehicle
dimensions were estimated for all vehicle types. These dimensions can be found in `dataobjects/enums/vehicletype.py`.

**Mixed files in NGSim**

At two locations of the NGSim project, Peachtree and Lankershim, the collected data was not split into files per time slot. It was all stored in a single file,
but the vehicle IDs were reset at the second timeslot. To display the data properly, it has to be split into different datasets. This is done automatically when
loading one of these files. The functions that do this can be found in `processing/NGSIMsplitting.py`. An added issue with the Peachtree dataset is the fact
that the timestamps are also corrupted. The first timestamp is correct but from that point on 100 seconds are added every time step instead of 100 ms. This is
also fixed in the same function.

**Typos in NGSim**

Some NGSim files contain typos in the header. This is taken care of when loading the NGSim *.csv file in `dataobjects/ngsimdataset.py` (`read_ngsim_csv`
function).

**NGSim data in feet**

All NGSim data is provided in feet. All values are converted to SI-units in `dataobjects/ngsimdataset.py` (`read_ngsim_csv` function).

**Velocity in km/h in Pneuma data**
The PNeuma data is provided with velocities in km/h, this is converted to m/s in `dataobjects/pneumadataset.py` (`read_pneuma_csv` function).

**Background for ExiD dataset 00-18**
The background image for ExiD datasets 00 to 18 does not seem to match the recorded data. This issue also exists in the visualizer provided with the ExiD 
data. Because there is no fix possible without access to the raw data and specific location of the recording, this issue remains open. 

### Where to place the data

If you want to load the data automatically, it should be stored in the data folder. This folder contains sub-folders for every project, these sub-folders
contain README.txt files explaining the folder structure in this sub-folder. In general, the procedure is the same as in the example above, download the 
data and
extract the files in the data folder as specified in the README.txt file.

### How data is stored

All data comes in *.csv files, these are great for general use but take some time to convert to usable Python objects. To shorten the loading time, Travia will
load the raw data from the *.csv files once and then store the created Python objects using pickle. However, opening pickle objects from an unknown source poses
a security risk. To minimize this risk, Travia encrypts the pickle files using a key that will be generated and saved in your home folder at first use.

This ensures that you will only open files that were created on your computer (or with access to your key). However, it also means that you cannot share *
.pkl files with other people or between computers. To enable your annotations to be shared, they will also be saved in a *.csv file. These annotation files will
be found and loaded automatically when first creating the python objects.

### How to select a dataset for visualization
There are three ways to specify which dataset is loaded for visualization. You can specify which dataset to load in the file `visualize.py` directly. To do 
this look for the code below `if __name__ == '__main__':` and uncomment the part needed to load a specific dataset, comment the other parts. As you can see, 
each dataset is identified with a DatasetID enum object. The NGSim and PNeuma dataset are split into parts based on recording location and time, the HighD 
data is split into indexed sub-sets. The enum objects correspond to the divisions made in the data source. For the example used in the getting started 
part, the dataset_id should be a dataset ID corresponding to a dataset from the `NGSIM-US101`-location, in this specific case, it could be `NGSimDatasetID.
US101_0805_0820`. Here, `US101` refers to the location of the recording and `0805_0820` to the time of the recording (between 08:05 and 08:20 AM).

Alternatively, you can pass arguments from the command line to specify which dataset to load. This will override any `dataset_id` specified in `visualize.py`.
To specify which dataset to load a source and a dataset ID need to be specified. The source argument (`-s`) should be one of {highd, ngsim, pneuma}. The 
dataset argument (`-d`) should be the name of the enum value, e.g. to load HighD dataset 1 run `python visualize.py -s highd -d DATASET_01`.

If no dataset ID is specified in the file or arguments, TraViA will ask the user what dataset to use with a dialog. This dialog will only show the dataset 
IDs corresponding to data files that can be found in the data folder.  

The script to start TraViA is a file called `visualize.py` and it can be found in the main project folder. This file is where you indicate which dataset you
want to load. Open the file now and look for the code below `if __name__ == '__main__':`. Uncomment the part needed to load an NGSim dataset, and comment the
other parts. 

### HighD example tools

To illustrate the possibilities for extension of Travia, three example functionalities have been implemented on the HighD dataset

- Rendering an overlay heatmap
- Plotting annotation data
- Automatic annotation of situations

The rendering of an overlay and plotting of annotation data are both implemented in `gui/gui.py` (`_toggle_overlay` and `_create_plots` respectively). These
function have comments and should be enough to get you started on implementing your overlay of plots.

The automatic annotation of situations can be performed using the `annotate_highd.py`. It will automatically detect interesting scenarios and add annotations
accordingly. please have a look at `processing/automaticannotationhighd.py` for more information.

For more help with these example implementations or if you have specific questions about them, please submit an issue on GitHub.

## Contributing to TraViA

Contributions to TraViA are very welcome, if you want to make a contribution you can open a pull request on GitHub. Please stick to pythons pep-8 guidelines and
clean-coding standards in general. The only deviation from pep-8 (that was made on purpose) is the maximum line width, this was increased from 120 characters to
160 character because we all have wide-screen monitors.

Also, please keep in mind the design choices that were made when creating TraViA. To summarize them: TraViA is created to serve as a broad basis for research on
traffic data in many disciplines. It was created such that it can easily be extended for a specific use but was not meant to be imported in other projects.
For these reasons, contributions that are too specific and are build for only one research project will be rejected (but can be made open access in separate
repositories). Contributions that are general enough for use in different projects but that are only useful for a sub-group of traffic data users might be added
to TraViA in a separate branch.

If you have questions about using or extending TraViA, or if you found a bug, please open an issue on GitHub.
--- 
title: 'TraViA: a Traffic data Visualization and Annotation tool in Python'
tags: 
  - naturalistic traffic data
  - visualization
  - annotation
  - HighD
  - NGSim
  - pNEUMA
  - Python 
authors: 
  - name: Olger Siebinga
    orcid: 0000-0002-5614-1262 
    affiliation: 1 
affiliations: 
 - name: Human-Robot Interaction group, Department of Cognitive Robotics, Faculty 3mE, Delft University of Technology, Mekelweg 2, 2628 CD Delft, the Netherlands
   index: 1 
date: 24 June 2021
bibliography: paper.bib
--- 

# Summary

In recent years, multiple datasets containing traffic recorded in the real world and containing human-driven trajectories have been made available to researchers.
Among these datasets are the HighD, pNEUMA, and NGSIM datasets. TraViA, an open-source Traffic data Visualization and Annotation tool
was created to provide a single environment for working with data from these three datasets. Combining the data in a single visualization tool enables
researchers to easily study data from all sources. TraViA was designed in such a way that it can easily be extended to visualize data from other datasets and
that specific needs for research projects are easily implemented.

# Statement of need

The combination of drones, cameras, and image recognition techniques might sound like a recipe for a spy movie. But actually, this combination allows for the
collection of rich traffic datasets. The recipe is straightforward: hover a drone above a location with traffic, record a video, and use image
recognition to generate bounding boxes for all vehicles. The result is a dataset containing human-driven trajectories at the location of interest that can be
used for many scientific purposes, e.g., to study traffic flow, model human behavior, or design autonomous vehicle controllers.

Because the required ingredients are easily accessed all over the world, multiple such datasets have been published in recent years. In Germany, the highD
project [@Krajewski2018] recorded all traffic at 6 different high-way locations; in Athens, Greece, all traffic in the city's business district was recorded
using 10 drones for 5 days in the pNEUMA project [@Barmpounakis2020]; and American highway traffic was recorded using fixed base cameras in the NGSIM
project [@NGSIM2016]. Combined, these datasets span different countries, types of vehicles, and environments, a combination valuable for researchers with
different backgrounds. Example usages of these datasets are validating human behavior models (e.g., by @Talebpour2015a and @Treiber2008) or testing autonomous vehicle controllers (e.g., by @Schwarting2019). 

Currently, it is difficult to leverage the powerful combination of multiple datasets because all the datasets come in different formats, and it is often
difficult to get a good and real-time visualization of the data. Some visualization tools exist (one is provided with the highD data [@Krajewski2018] and another 
example for NGSIM data can be found in @Sazara2017) but they are specifically made for only one of these datasets and are very basic in the sense that they 
provide little control over simulation time and no insight in raw values per vehicle per frame. In addition to difficulties with
visualization, finding, and annotating situations of interest in these massive datasets is a time-consuming task and keeping track of the annotations for the
different datasets requires some bookkeeping skills.

TraViA was developed to provide a solution for these problems. TraViA can be used to visualize and annotate data from highD, pNEUMA, and NGSIM and uses
generic vehicle objects to store the state of vehicles at a specific time. This makes it possible to validate and test models or controllers on multiple 
datasets in parallel, without having to cope with the different dataset formats.

# Software Functionality 

TraViA is written in Python 3 and has a graphical user interface developed in PyQt5. A screenshot of TraViA is provided
in \autoref{fig:screenshot}. This screenshot shows the capabilities of TraViA in a single image. The main features of TraViA are:
 
* Advanced information display based on raw data for every vehicle in every dataset by leveraging generic vehicle objects
* Dynamic visualization of the traffic scene with possibilities to zoom, pan, and rotate for an optimal view
* Exporting the visualization to a video or single image
* An interactive timeline that shows dataset annotations, which are saved as python objects for easy manipulation 

![A screenshot of the TraViA software visualizing a frame of the highD
dataset. The main features of TraViA are highlighted in this image. \label{fig:screenshot}](images/screenshot.png)

TraViA was designed for use as a stand-alone program. It uses abstract classes as a basis for all dataset-specific objects to enable easy implementation of 
new datasets (for a class diagram and more information on how to do this, please see the readme file in the repository). It was specifically created to serve 
as a tool for generic visualization and annotation such that it can be used by researchers from different 
fields. To show the capabilities of TraViA and to provide a starting point for other researchers that want to use TraViA for their work, three example 
implementations of tools for specific purposes are included with TraVia. The first example is the functionality to automatically detect and annotate 
specific scenarios (e.g., lane changes), the second is functionality to plot specific vehicle signals over the course of an 
annotation, and the third is a function to plot a heatmap overlay for use in autonomous vehicle reward function development. All of these example tools are only 
implemented for use with the highD dataset.

# Usage of TraViA in Science and Education
Currently, TraVia is being used by the author for model validation of an inverse-reinforcement-learning-based driver model. A publication on this validation 
is currently being prepared for submission. Besides that, TraViA is used for educational purposes, allowing students at TU Delft to explore big naturalistic 
datasets by providing them with an accessible, GUI-based starting point.

# Acknowledgements

I thank Nissan Motor Co. Ltd. for funding this work, and I also thank my supervisors David Abbink and
 Arkady Zgonnikov for their valuable help and advice. Finally, I thank Joris Giltay for his help with testing the instructions in the readme file. 

# References
