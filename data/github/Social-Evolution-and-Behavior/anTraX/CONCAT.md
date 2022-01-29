[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/Social-Evolution-and-Behavior/anTraX)](https://github.com/Social-Evolution-and-Behavior/anTraX/releases)
![GitHub commits since latest release](https://img.shields.io/github/commits-since/Social-Evolution-and-Behavior/anTraX/latest?color=red)
![GitHub last commit](https://img.shields.io/github/last-commit/Social-Evolution-and-Behavior/anTraX?color=red)
[![Read the Docs](https://img.shields.io/readthedocs/antrax)](http://antrax.readthedocs.io)

![trails](https://github.com/Social-Evolution-and-Behavior/anTraX/blob/master/docs/images/trails.png)

# **anTraX**: high throughput tracking of color-tagged insects

anTraX is a software for video tracking of ants and other small animals tagged with a unique pattern of color dots. It was designed for behavioral experiment using the Clonal Raider Ant [*Ooceraea biroi*](https://en.m.wikipedia.org/wiki/Ooceraea_biroi), but can be used for any other model system. anTraX is a **brute force** type tracking algorithm, which was designed to handle high throuput long duration experiments (many colonies over many days). Therefore, it will require considerable computational resources. 

The software was designed and written by Jonathan Saragosti and Asaf Gal of the [Laboratory of Social Evolution and Behavior](https://www.rockefeller.edu/research/2280-kronauer-laboratory/) in the Rockefeller University, and is distributed under the [GPLv3](https://github.com/Social-Evolution-and-Behavior/CATT/blob/master/LICENSE) licence.


## Requirements

anTraX works natively on machines running Linux or OSX operating system. It benefits significantly from a multicore system. It is recommended to have at least 2GB of RAM per used core, and a similar sized swap. Computational GPU will speedup the classification phase considerabley. 


## Installation and usage

See online documentation at https://antrax.readthedocs.io/. An example dataset to be used in the tutorial can be found [here](https://github.com/Social-Evolution-and-Behavior/anTraX-data). All datasets used to benchmark anTraX are [available from the Zenodeo repository](https://zenodo.org/record/3740547).

## References

The anTraX preprint is now [available on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.04.29.068478v1).

### Requirements

anTraX works natively on machines running Linux or OSX operating system. It will benefit significantly from a multicore system. It is recommended to have at least 2GB of RAM per used core, and a similar sized swap. Computational GPU will speedup the classification phase considerably. 

** Python:**

anTraX requires Python version 3.6 or above. It is highly recommended to install and use anTraX inside a virtual environment, using conda or any other environment manager. Full list of dependenices and version numbers can be found [here](antrax_dependencies.txt). 

** MATLAB:** 

anTraX comes with binraries compiled with MATLAB 2019a. If you have MATLAB 2019a or above installed on your machine, you do not need to install the Runtime. In case you don't, please follow the instructions to install the freely available [MATLAB Runtime](https://www.mathworks.com/products/compiler/matlab-runtime.html) for version 2019a.

** FFmpeg:**

anTraX uses [FFmpeg](https://www.ffmpeg.org/) to read video files.

To install on Ubuntu, open a terminal and type:

```console
sudo apt install ffmpeg
```

To install on OSX using [homebrew](https://brew.sh/), open a terminal and type:

```console
brew install ffmpeg
```
### Get anTraX

Change into your favorite place to install code packages, then clone the anTraX repository:

```console
git clone http://github.com/Social-Evolution-and-Behavior/anTraX.git
```

To install a specific version, use the git checkout command:

```console
cd anTraX
git checkout <version>
```

Alternatively, anTraX the latest version of anTraX can be directly dowloaded [here](https://github.com/Social-Evolution-and-Behavior/anTraX/archive/master.zip). 


### Install anTraX

**Install the python package: **

To install the python package, run in the anTraX folder:

```console
pip install .
```
** Set envoronment variables: **

If you are using a full MATLAB installation, add these lines to your bash profile file (usually `~/.bash_profile` on OSX,   `.profile` on Linux):

```bash
export ANTRAX_PATH=<full path to anTraX repository>
export ANTRAX_USE_MCR=False
```

Otherwise, if you are using MCR, add these lines:

```bash
export ANTRAX_MCR=<full path to MCR installation>
export ANTRAX_PATH=<full path to anTraX repository>
export ANTRAX_USE_MCR=True
```

For the changes to take effect, run:
```console
source ~/.bash_profile
```

** Setup MATLAB (full MATLAB mode):**

**Note**: You do not need to do these steps if you are using the compiled binaries!

In case you have MATLAB installed (i.e. you are not using the compiled binaries), you will need to install the [MATLAB engine for python](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html). Take care to install this in the same python environment you installed anTraX to.

You will also need to compile a mex file. In the matlab console, type:

```console
cd <path-to-anTraX>/matlab/external/popenmatlab
mex popenr.c
```
If you intend to use anTraX within an interactive MATLAB session (requires an active license), add the repository to your search path. In the MATLAB console, type:

```console
addpath(genpath(<path-to-anTraX>/matlab));
```
Save the path if you want it to hold for future sessions.

**Developer mode:**

If you use a full MATLAB installation, changes to the MATLAB source code of anTraX will take effect immediately, as it is always run from the repository directory. If you use MCR, you will need to compile your changes for them to take effect.

Changes made to the python source code will take effect only after reinstalling the package. Alternatively, for changes to take effect immediately when you make them, you can install the package in editable mode:

```console
cd anTraX
pip install -e .
```

**Installation on computer clusters:**

Please refer to the [HPC section](hpc.md#installation).
### Classification workflow

In this step we will classify each tracklet that was marked as a possible single-ant tracklet. Each of these tracklets will be assigned as either:

* An ID from the list of possible IDs
* A non-ant tracklet, either as a general category or a specific one if such exist in the classifier.
* A multi-ant tracklet
* An ambiguous tracklet ('Unknown') in case the classifier couldn't make a decision.  

Each tracklet is classified by first classifying all blob images belonging to that tracklet, using the *blob classifier* and then weighting these classifications to produce a whole tracklet classification. The blob classifier is a trained deep convolutional network (CNN), that needs to be trained on a trainset of pre-classified blob images. 

### Creating a training set

anTraX includes an interactive GUI application to prepare such a training set from a tracked experiment. To launch the GUI, type and run in the terminal:

```console
antrax extract-trainset <expdir> [--session <session>]
```

The app will display blobs images from a randomly selected tracklet, as well as info about the tracklet in the text-box below. Using the buttons on the right, you can select the label appropriate for the tracklet, and either export all images to the train set (using the ***export all*** button), or select a subset using the 'select frames' button (a frame selection window will open). You can move between tracklets using the ***Next*** and ***Back*** buttons.
See the [tips and best practices](tips.md#training-and-classification) page regarding what are good training set examples.

By default, the app will load tracklets from the first tracked movie. Tracklets from other movies can be loaded using the ***Tracklets*** menu.

![trainset extraction](images/extract-trainset1.png)

![frame selection](images/extract-trainset-frame-selector.png)


### Merging training sets

The exported examples are saved as images in the experimental directory, under `/session/classifier/examples/id/`. When creating a classifier specific to that experiment, training can be done from that directory. However, in order to create a general classifier to be used on many experiments, it is necessary to create a training set that contains examples from several experiments. This can be done with the command:

```console
antrax merge-trainset <source-classdir> <dest-classdir>
```

This will merge all the examples from the source classifier directory (usually `expdir/session/classifier`) into the destination classifier directory. It is recommended to keep multi-experiment classifiers separate from any specific experimental directory to avoid confusion.

**Important:** the user is responsible for making sure the lists of labels match when merging trainsets. Otherwise, problems might occur.

### Training the classifier

To train the classifier, run: 

```console
antrax train <classdir> [OPTIONS]
```

The `classdir` argument can be either a full path to the experimental directory containing the classifier, or to the directory containing the `examples` directory, and where the classifier will be saved.

The train command accepts the following options:

`--scratch`

By default, anTraX will load a pre-trained classifier if it exists in the `classdir` directory, and run incremental training. If the `--scratch` flag is used, it will initialize a new classifier instead. 

`--ne <ne>`

Number of training epochs to run. Default is 5.

`--hsymmetry`

Use this flag if your tagging is symmetrical to horizontal flips (this is used in the dataset augmentation process). This option can be used either for a new classifier, or together with the `--scratch` options.

`--target-size <size>`

The side length (in pixels) of the input image to the classifier (anTraX always use square images for classification). This option will only take effect if a new classifier is trained. By default, the size will be that of the first image read from the trainset. All other images will be resized to the target-size. This option can be used either for a new classifier, or together with the `--scratch` options.

`--name <name>`

Use a custom name for the classifier. 

`--arch <arch-name>`

If you want to use a non-default architecture for the CNN (see next section for options). Acceptable values are `small|large|wide|MobileNetV2|custom`. In case `custom` is chosen, you will need to also include the `modelfile` option. It is highly recommended to start with the default architecture, and explore different ones only for complex problems. This option can be used either for a new classifier, or together with the `--scratch` options.

`--modelfile <jsonfile>`

A path to a json file containing a CNN model [serialized by keras](https://www.tensorflow.org/guide/keras/save_and_serialize). Note that only the model architecture is loaded, not a trained model. The number of classes in the model must match the number of classes in your tracking problem. This option can be used only together with `--arch custom'.

`--aug-options <options>`
  
A comma separated list of augmentation parameters. This option can be used to change some of the augmentation parameters used by TensorFlow to train the classifier. Currently supports the options `width_shift_range` (default value 10), `height_shift_range` (default value 10), `shear_range` (default value 5), `rotation_range` (default value 15), `zoom_range` (default value 0.1). Details of these parameters can be found in the [Keras documentation of  ImageDataGenerator](https://keras.io/api/preprocessing/image/). Example `--aug-options shear_range=10,rotation_range=20`.


### anTraX-supplied CNN architectures

anTraX defines a few possible CNN architectures. To see a complete specification of the architectures, look at the [models.py](https://github.com/Social-Evolution-and-Behavior/anTraX/blob/master/antrax/models.py) file in the anTraX repository.

* **small**: This is the default architecture, which we have found to give the best trade-off between accuracy, training time, and number of examples needed. It has 3 convolutional layers.  

* **wide**: Also a 3-layered model, but with wider layers and hence slower to train.

* **large**: A 4-layered model.

* **MobileNetV2**: A pre-trained [MobileNetV2](https://arxiv.org/abs/1801.04381) model. 


### Classifying tracklets

To run classification in batch mode:

```console
antrax classify <experiments> [OPTIONS]
```

The `experiments` argument can be either a full path to an experimental directory, a full path to a text file with a list of experimental directories (all of which will run in parallel), or a full path to a folder that contains one or more experimental directories (all of which will run in parallel).

The classify command accepts the following options:

`--classifier <path-to-classifier file>`

Explicit path to a classifier (`.h5` file created by the train process). By default, anTraX will use the classifier file that exists in the default location in the experimental directory `expdir/session/classifier/classifier.h5`. If it doesn't exist, an error will be raised.

`--movlist <list of movie indices>`

By default, anTraX will track all movies in the experiment. This can be changed by using this option. Example for valid inputs incluse: `4`, `3,5,6`, `1-5,7`.

`--session <session name>`

If your experiment contains more than one configured session, anTraX will run on the last configured one. Use this option to choose a session explicitly.


### Validating classification and retraining 

Once the tracklets in the experiment are classified, the *extract-trainset* app can be used to validate the results and export new examples if needed. For a classified experiment, the app will display the assigned label for each displayed tracklet. If it is incorrect, or if it is unclassified (labeled 'Unknown' although it is identifiable), you can add its images to the training set the same way as before. Don't forget to choose the correct label before exporting!

The ***Filter by autoID*** option in the ***Tracklet*** menu can be used to show only tracklets assigned with a specific label. 

Optionally, you can formally evaluate the performance of the classifier by marking each classification as ***Correct***, ***Wrong*** or ***Should not have an ID*** (if the tracklet is not identifiable). Doing that for a set of tracklets will give an estimate of the classifier performance per-tracklet and per-frame.

Once the training set has been expanded, train the classifier and re-run classification in the same way as before.

![trainset extraction](images/extract-trainset2.png)
![trails](images/trails.png)

# **anTraX** -  high throughput video tracking of color tagged insects

anTraX is a software for video tracking ants and other small animals that are marked with a unique pattern of color dots. It was originally designed for behavioral experiments using the Clonal Raider Ant [*Ooceraea biroi*](https://en.m.wikipedia.org/wiki/Ooceraea_biroi), but can be used for any other model system. anTraX is a **brute force** type tracking algorithm, prioritizing accuracy over speed, therefore requiring considerable computational resources. It was designed to handle high throughput, long duration experiments (many colonies over many days), and benefits from running on computer clusters.

The software was designed and written by Jonathan Saragosti and Asaf Gal of the [Laboratory of Social Evolution and Behavior](https://www.rockefeller.edu/research/2280-kronauer-laboratory/) in the Rockefeller University, and is distributed under the [GPLv3](https://github.com/Social-Evolution-and-Behavior/anTraX/blob/master/LICENSE) licence.

###References

Gal A, Saragosti J and Kronauer DJ. **anTraX: high throughput video tracking of color-tagged insects**. *bioRxiv* (2020). [doi:10.1101/2020.04.29.068478](https://www.biorxiv.org/content/10.1101/2020.04.29.068478v1)

### Marking insects

* Choose colors that are easily separated in RGB space. As a rule of thumb, colors that are easily separated by a human eye in the video will result in better tracking (surprise!). The color set Blue, Green, Orange, Pink, and Yellow usually results in good accuracy.

* It is important to make sure that all tags are visible to the camera for a large part of the frames, at least when an ant is on the move. For that reason, placing tags one behind the other is preferable over side-by-side tagging . 

* It is not advised to use a variable number of tags per ant (for example, using "no tag" as an additional color). As tags are often not visible (e.g., when an ant is grooming or on its side), this can lead to a high error rate (i.e. obscured tag is identified as a no-tag ID).

* It is possible to use symmetric color combinations (e.g., Blue-Green and Green-Blue) in the same experiment. However, note that in that case, an asymmetry in the ant appearance in the image is necessary for good classification. For example, for *O. biroi*, a thorax tag and abdomen tag will work nicely, assuming the head and the antennae are breaking the image symmetry (don't paint the head!).    

### Recording videos

* The program is designed to identify dark ants on a light background. The better the contrast, the easier it is to segment the ants from their background. 

* It is recommended to tune the camera, especially the brightness, contrast, and saturation levels, to make sure tag colors are maximally distinguishable. A bit of over saturation is usually advised. 

* As the software segments the image using background subtraction, it is important to make sure the image is stable throughout the experiment. Particularly, try to minimize camera movement, arena movements, and illumination fluctuations. 

* Sometimes, it is necessary to pause the experiment and to take care of the experimental colonies. It is recommended to start a new video subdirectory to hold movies if the arenas cannot be returned to their exact position relative to the camera, or if there is a change in the arena that can result in a changing background. Backgrounds can be created separately for each subdirectory. 

* Video resolution should be high enough to identify the color tags. As a rule of thumb, a tag size of at least 5x5 pixels is recommended. 

* Video frame rate should be high enough to make sure ant blobs in consecutive frames overlap considerably. For experiments with the clonal raider ant *O. biroi*, a framerate of 10 fps is enough.

### Computer configuration and parallel execution 

* A tracking job will typically use around 125-150% CPU. Therefore, it is recommended to set the number of MATLAB workers to about 0.5-0.75 of the number of threads supported by your machine. It is also recommended to have at least 2GB of RAM and 4GB of swap per worker. 

* Tracking jobs are usually long. It is advised to use a desktop computer or workstation and not a laptop computer. Make sure to turn off auto sleep on the computer. 

### Tracking settings and parameters

* What is a good background

* median vs max

* What is a good ROI mask

* Good segmentation 

* Why single ant threshold is important


### Training and classification

* Choosing good examples for training a classifier is an art. On the one hand, it is a good practice to span the space of possible presentations of the ants, and choose "atypical" images. On the other hand, as this space is huge, we are practically guaranteed to be in an under-sampled regime, so it is better not to use examples in ambiguous areas of the image space. As a rule of thumb, don't use images that are not easily identifiable for a human. 
  
* If you see a systematic classification error between two IDs, especially if there is no clear visual reason for this error, it is likely that you have some contamination of the training set. Try checking the example directory of the classified label and look for examples that belong to the actual ant ID. Delete them and retrain the classifier.

* Neural networks are best trained using negative examples. Otherwise, they will produce high rate of false positives in cases where none of the known labels exist. In our case, this is represented by the "Unknown" class (and the "Multi" class if it is included). It is important to add as many examples as possible to this class. Good images are those that a human cannot identify at all. Tracklets labeled as "Unknown" are left for the propagation algorithm to assign.

* Why are some tracklets identified as "Unknown" although they are very obvious to recognize? This might be due to a few possible reasons. It might be just a sampling issue: as short tracklets are, well, short, maybe their images happen to be in a region of image space not covered by the training set (remember that assignments requires high confidence, so even though the blob classifier assigns it correctly, the algorithm removes the assignment due to low confidence). In that case, adding these images to the training set will help future rounds of classification.  Alternatively, the training set is contaminated, so a conflicting example reduces the confidence of the classifier (it is always recommended to check the set every now and then). 

### ID propagation 


### Debugging and fixing the tracking

### Working with xy trajectories

* ### The results files organization

### Loading tracking results in Python


```python
from antrax import *

ex = axExperiment(<expdir>, session=None)
antdata = axAntData(ex, movlist=None, antlist=None, colony=None)

```

### Loading tracking results in MATLAB

```matlab
ex = trhandles.load(<expdir>, session=None)
antdata = loadxy(ex, 'movlist', movlist, 'colony', colony)
```### Release notes for version 1.0.2

Some organizational changes. This is the version submitted for review.

### Release notes for version 1.0.1

Bug fix with running graph-explorer in MCR mode.

Documentation updates.

### Release notes for version 1.0.0

First public release

### Release notes for version 0.1.0

First software release 


### Validating tracking results

Generally speaking, the performance of the tracking algorithm can be captured using two separate measures. The first is the rate of assignment, defined as the ratio of assigned locations in the experiments to the total possible assignments (i.e. the number of IDs times the number of frames). The second measure is the assignment error, defined as the ratio of wrong assignments to the total number of assignments made. While the assignment rate can be computed directly and precisely from the tracking results, the error rate in assigning IDs for a given data set needs to be tested against human annotation of the same dataset. Often, it is impossible to do because the recording duration of these datasets is typically long (many hours); therefore it is impractical to manually annotate them in full. Moreover, the performance of tracking will usually depend on many factors that are specific to a given dataset, so for every new experiment we will need to redo the error estimation procedure.

As an alternative, we use a validation procedure in which a sequence of randomly selected test points is presented to a human observer, where each test point corresponded to a location assignment made by the software to a specific ID in a specific frame. The user is then asked to classify the assignment as either ‘correct’ or ‘incorrect’. If the user is unsure of the correctness of the assignment, they could skip to the next one. The number of times this process is repeated depends on the required power of the error estimation. 

### Using the validation app

anTraX includes a graphical interface to easily validate tracking results. To launch:

```console
antrax validate <expdir> [--session <session>]
```

The app will randomly select an assigment of an ID to a blob. The user should respond by pressing either ***Correct***, ***Incorrect***, or ***Can't say***. The text box will display the accumulated error and its confidence interval. Use the ***Assigment*** menu to estimate the error for a specific ID, time range, or a colony (for multi-colony experiments).

![Validation app](images/validation1.png)

### Run a batch job 

Once the session is configured, the next step will be to run the blob-tracking step on all the videos in the experiment (see the anTraX publication for details). To start tracking, simply execute the following command in a terminal:

```console
antrax track <experiments> [OPTIONS]
```

The `experiments` argument can be a full path to an experimental directory, a full path to a text file with a list of experimental directories (all of which will run in parallel), or a full path to a folder that contains one or more experimental directories (all of which will run in parallel). Note that each of the experiments needs to be configured separately before running the batch job. 

The track command accepts the following options:

`--nw <number of workers>`

anTraX will parallelize the tracking jobs by video. By default, it will use two MATLAB workers. Depending on your machine power, this can be changed by using this option. 

`--movlist <list of movie indices>`

By default, anTraX will track all movies in the experiment. This can be changed by using this option. Example for valid inputs incluse: `4`, `3,5,6`, `1-5,7`.

`--session <session name>`

If your experiment contains more than one configured session, anTraX will run on the last configured one. Use this option to choose a session explicitly.

### Checking job progress

Depending on the number and length of video files, tracking jobs can be very long. anTraX will print a report to terminal when a task (tracking of single video) starts/ends. Logs for each task can be found in the experimental directory, under `session/logs/`. Note that depending on your machine settings, the log file might not be updated in real time.


### Propagating IDs on tracklet graphs

To run propagation in batch mode:

```console
antrax solve <experiments> [OPTIONS]
```

The `experiments` argument can be either a full path to an experimental directory, a full path to a text file with a list of experimental directories (all of which will run in parallel), or a full path to a folder that contains one or more experimental directories (all of which will run in parallel).


The solve command accepts the following options:

`--nw <number of workers>`

anTraX will parallelize the tracking by video. By default, it will use two MATLAB workers. Depending on your machine power, this can be changed by using this option. 

`--glist <list of graph indices>`

By default, anTraX will track all graphs in the experiment. This can be changed by using this option. Graph indices are an enumeration of movie groups according to the options set in the session configuration. Example for valid inputs include: `4`, `3,5,6`, `1-5,7`.

`--clist <list of colony indices>`

If the experiment is multi-colony, by default, anTraX will run all colonies in the experiment. This can be changed by using this option. Colony indices are 1 to the number of colonies. Example for valid inputs include: `4`, `3,5,6`, `1-5,7`.


`--session <session name>`

If your experiment contains more than one configured session, anTraX will run on the last configured one. Use this option to choose a session explicitly.


### Using the graph explorer to view and debug ID assignments

To launch:

```console
antrax graph-explorer <expdir> [--session]
```

![Graph Explorer](images/graph-explorer1.png)

* To load a graph, choose a movie (and a colony if applicable) under the 'Graph' menu.
* To highlight a specific ID, choose it in the 'View' menu. The subgraph of that ID will highlight in the graph plot, and its trajectory will highlight in the xy plots on the left.
* To see details of a specific tracklet, select its node in the graph plot. The details will appear in the text box on the right. The part of the trajectory corresponding to that tracklet will show as thickened lines on the xy plots on the left.
* When selecting a single-animal tracklet, the first cropped image of that tracklet will appear in the top left panel. You can move between frames by using the drop down menu above the tracklet text box, or press 'Image all frames' to open a separate window with a montage view of all the cropped images of the tracklet. When a multi-animal tracklet is selected, you can load its cropped image by pressing the 'Image' button.
* Use the 'Show frame' button to show an image with the blob marked on its corresponding video frame.
* You can assign an ID to a tracklet to override or augment the results of the propagation algorithm. You can also mark a tracklet as single-animal. Note that you will need to rerun the propagation step to observe the effects of these fixes.
* You can export the cropped images of a specific tracklets using the buttons on the bottom row. This is useful in case you tumble upon a misclassified tracklet.### DeepLabCut

[DeepLabCut](http://www.mousemotorlab.org/deeplabcut) is a popular software that uses deep neural networks  for pose-tracking animals in videos, developed by the [Mathis lab](http://www.mousemotorlab.org/) at Harvard University. 
anTraX includes an interface to export cropped single-animal examples from tracked experiments to be labeled in the DeepLabCut interface, as well as options to run trained DLC models in the anTraX interface without the need to export cropped videos from entire experiment. This integration allows to create efficient pipelines that pose-track marked individual animals in large groups.

### The anTraX/DLC workflow

* Track an experiment using anTraX.
* Export a subset of single-animal images from anTraX to DeepLabCut.
* Train a DLC model
* Use the trained DLC model to track all single ant tracklets in the anTraX session. This is done using a custom anTraX function that feeds cropped single ant videos to DLC.
* The estimated bodypart positions from DLC is saved and loaded together with the centroid position of the ants.  

### Install DeepLabCut

Install the DeepLabCut package into your python environment:

```console
pip install deeplabcut
```

### Export single ant videos for training 

To export training images from anTraX to DLC project, run:

```console
antrax export-dlc-trainset <expdir> <dlcdir> [OPTIONS]
```

This will export example single ant frames from the experiment in `expdir` to the DeepLabCut project directory `dlcdir`.  If `dlcdir` doesnt exist, a new DLC project will be created (a date string will be appended to that directory name per the DeepLabCut convention).

`--nimages <nimages>`

By default, 100 randomly selected images will be exported. Change this using this option.

`--movlist <movlist>`

By default, anTraX will select frames from all movies in the experiment. This can be changed by using this option. Example for valid inputs incluse: `4`, `3,5,6`, `1-5,7`.

`--antlist <antlist>`

By default, anTraX will select frames from all ants in the experiment. This can be changed by using this option, providing a comma separated list of IDs. Example: `--antlist BB,GP,PO`.

`--video`

By default, the extracted frames will be saved as images in the DeepLabCut project directory. You can select to save them as videos using this flag.

`--session <session>`

If your experiment includes more than one configured session, anTraX will export examples  from the last configured one. Use this option to choose a session explicitly.

### Train a DeepLabCut model

At this point, you are ready to train your DeepLabCut model. Refer to the package [webpage](http://www.mousemotorlab.org/deeplabcut) for a tutorial.

### Run

Once the model is trained and ready, you can use the anTraX interface to easily run it on the tracked experiment:

```console
antrax dlc <experiments> --cfg <path-to-dlc-cfg-file> [OPTIONS]
```

The `experiments` argument can be either a full path to an experimental directory, a full path to a text file with a list of experimental directories (all of which will run in parallel), or a full path to a folder that contains one or more experimental directories (all of which will run in parallel).

The required argument `cfg` is the full path to the project file of the trained DeepLabCut model.

The classify command accepts the following options:

`--movlist <list of movie indices>`

By default, anTraX will process all movies in the experiment. This can be changed by using this option. Example for valid inputs incluse: `4`, `3,5,6`, `1-5,7`.

`--session <session name>`

If your experiment contains more than one configured session, anTraX will run on the last configured one. Use this option to choose a session explicitly.

### Loading and analyzing postural data

The DeepLabCut pose tracking results can be loaded using the anTraX python interface:

```python
from antrax import *

ex = axExperiment(<expdir>, session=None)
antdata = axAntData(ex, movlist=None, antlist=None, colony=None)
antdata.set_dlc()
```

Note that loading pose tracking data for a full experiment might take some time if the experiment is long. Consider loading a partial dataset using the `movlist` argument for initial exploration before doing a full analysis.

Refer to the example jupyter notebook for an elaborated example.






### The experiment directory

The input and output of anTraX are organized in an  "experiment directory" (a.k.a. `expdir`). Under this directory, the program expects to find a subdirectory named `videos`, containing the input files. While running anTraX, a subdirectory for each tracking session will be created under the experimental directory and will store the parameters and results for that session.

### Video files

The input files are the raw videos of the same animal group, recorded sequentially under the same conditions. anTraX can process any [file format and video encoding readable by ffmpeg](http://www.ffmpeg.org/general.html#Supported-File-Formats_002c-Codecs-or-Features); however, all videos are expected to have the same file format, codec, frame size, and framerate. Moreover, all videos must have the same base name, followed by a file index suffix separated by a `_` character. 
Optionally, for convenience, the videos can be organized into subdirectories. This can be a useful way to separate  some meaningful partition of the experiment such as periods of consecutive recording, change in experimental conditions or feeding events, or just to partition a very long experiment.  The subdirectories should indicate the indexes of the videos stored in them. For example, subdirectory named "1_24" will contain videos 1 to 24.

![expdir structure](images/expdir_structure.png "structure of the experimental directory")

### The frame data files

Optionally, each video will be accompanied by a frame information file, with the same name of the video and with a `.dat` extension. This file should contain a header row with the variable names (e.g. a timestamp, sensor data, etc.) and a value row for each frame in the video. If the file contains a variable named `dt`, it will be interpreted as the precise inter-frame interval, and will be used for tracking instead of the video frame rate parameter.

![dat file example](images/dat_file_example.png)

To see an example for an experimental directory, download one of the [example datasets](datasets.md).





### JAABA

[JAABA (the Janelia Automatic Animal Behavior Annotator)](http://jaaba.sourceforge.net/index.html) is a machine learning-based system created by the [Branson lab](https://www.janelia.org/lab/branson-lab) at HHMI Janelia Farm. It enables users to automatically compute interpretable, quantitative statistics describing video of behaving animals. In a nutshell, it uses a set of user-labeled examples to train a classifier that can spot more occurences of that behavior in new (unseen) data. JAABA works by projecting trajectory data into a high dimentional space of so-called "per-frame features", in which the underlying machine learning algorithm searches for regularities. 
anTraX includes an option to generate these perframe features for experiments tracked in anTraX, as well as functions to run JAABA classifier from the anTraX interface, thus greatly simplifying the use of JAABA for classifying behavioral features in these experiments. 

### The anTraX/JAABA workflow

* Track an experiment using anTraX. 
* Write tracks and per-frame data in JAABA-readable format.
* Train a classifier using the JAABA interface
* Classify the entire dataset using either the JAABA interface or anTraX interface.
* The JAABA-generated scores for each behavioral classifier will be imported together with the spatial coordinates of each animal.

### Install JAABA

Install the JAABA package from GitHub:

```console
git clone https://github.com/kristinbranson/JAABA.git
```

Copy anTraX configuration files into the JAABA directory:

```console
cp $ANTRAX_PATH/matlab/jaaba/*.xml $ANTRAX_JAABA_PATH/perframe/params/
```

Add this variable to your bash profile file:

```bash
export ANTRAX_JAABA_PATH=<full path to JAABA repository>
```
Don't forget to source!

### Write tracks and perframe data for JAABA


```console
antrax export-jaaba <expdir>  

```

The export-jaaba command accepts the following options:


`--nw <number of workers>`

anTraX will parallelize the tracking by video. By default, it will use two MATLAB workers. Depending on your machine, this can be changed by using this option.

`--movlist <movlist>`

By default, all movies will be processed, which might take some time. Change this using this option.

`--session <session name>`

If your experiment contains more than one configured session, anTraX will run on the last configured one. Use this option to choose a session explicitly.

*** Note: *** The `export-jaaba` command does not currently support the `--mcr` or the `--hpc` options.

### The JAABA directory structure

anTraX will create a directory called `jaaba` under the session directory. In that directory, a subdirectory for each movie in the experiment will be created (JAABA considers each movie to be a separate experiment, and will therefore call each of these subdirectories an experimental directory). In each of these directories, there will be a soft link to the movie, named `movie.mp4` (or `movie.avi` etc. if your video files extension is different), and a mat file called `trx.mat` containing the trajectories in JAABA-compatible format. A subdirectory called `perframe` will also be generated, and will hold the per-frame data.

![expdir structure](images/jaaba_directory_structure.png "jaaba directory structure")

### JAABA perframe features

JAABA computed a long list of perframe features (see the [JAABA original publication](https://www.nature.com/articles/nmeth.2281), supplementary material section 13  for full details). Many of these features are not valid for multi-animal tracklets (e.g. kinematic features and appearance features). Therefore, anTraX will write NaN values for these features in frames that corresponds to multi-animal tracklets.

### anTraX-specific perframe features

In addition to JAABA's list of perframe features, anTraX will also include an additional set of features:

* *antrax_blob_area*: The real blob area as reported by the segmentation algorithm. 
* *antrax_dblob_area*: The derivative of the real blob area.
* *antrax_dist_to_wall*: Distance to closest point on the ROI perimeter.
* *antrax_angle_to_wall*: Angle between blob orientation and the closest point on the ROI perimeter.
* *antrax_ddist_to_wall*: Derivative of the distance from ROI perimeter.
* *antrax_dangle_to_wall*: Derivative of the angle between blob orientation and ROI perimeter.
* *antrax_dist_to_center*: Distance to the arena's center (defined as the ROI center of mass).
* *antrax_angle_to_center*: Angle between blob orientation and the arena's center.
* *antrax_ddist_to_center*: Derivative of the distance to center.
* *antrax_dangle_to_center*: Derivative of the angle to the arena's center.
* *antrax_dist_to_openwall*: Distance to the closest open point in the ROI perimeter (NaN if ROI is fully closed).
* *antrax_angle_to_openwall*: Angle between blob orientation and the closest open point in the ROI perimeter (NaN if ROI is fully closed).
* *antrax_ddist_to_openwall*: Derivative of the distance to closest open point in the ROI perimeter (NaN if ROI is fully closed).
* *antrax_dangle_to_openwall*: Derivative of the angle between blob orientation and the closest open point in the ROI perimeter (NaN if ROI is fully closed).
* *antrax_dist_to_median*: Distance between the animal's centroid to the median location of all other animals.
* *antrax_angle_to_median*: Angle between the blob orientation and the median location of all other animals.
* *antrax_ddist_to_median*: Derivative of the distance between the animal's centroid to the median location of all other animals.
* *antrax_dangle_to_median*: Derivative of the angle between the blob orientation and the median location of all other animals.
* *antrax_nants_in_blob*: The number of individual animals in the blob.
* *antrax_frac_in_blob*: The normalized number of individual animals in the blob. 

### Training a classifier

Once the export step is finished, you are ready to use the JAABA interface to load the data and train the classifier. Refer to the [JAABA documentation page](http://jaaba.sourceforge.net/index.html) for information on this step.

### Applying the classifier to an experiment 


To run a trained JAABA classifier (defined a `.jab` file), run the command:

```console
antrax run-jaaba <expdir>  --jab <jabfile>
```

The run-jaaba command accepts the following options:

`--nw <number of workers>`

anTraX will parallelize the tracking by video. By default, it will use two MATLAB workers. Depending on your machine power, this can be changed by using this option.

`--movlist <movlist>`

By default, all movies will be processed, which might take some time. Change this using this option.

`--session <session name>`

If your experiment contains more than one configured session, anTraX will run on the last configured one. Use this option to choose a session explicitly.

*** Note: *** The `run-jaaba` command does not currently support the `--mcr` or the `--hpc` options.

### Loading and analyzing JAABA scores

For each ant in each frame, JAABA will assign a classification score. A positive score will imply a positive classification, and a negative score will imply a negative classification. The larger the absolute value of the score is, the stronger is the confidence in the classification. 

To load the results using the anTraX python interface:

```python
from antrax import *

ex = axExperiment(<expdir>, session=None)
antdata = axAntData(ex, movlist=None, antlist=None, colony=None)
antdata.set_jaaba()
```

Consider loading a partial dataset using the `movlist` argument for initial exploration before doing a full analysis.

Refer to the example jupyter notebook for a more elaborate example.


As described in the anTraX publication, the software was benchmarked using 9 different datasets, representing a variety of possible use cases. All datasets, including the anTraX configuration files used to track them are [available for download](https://zenodo.org/record/3740547).





![datasets table 1](images/datasets_table1.png)


![datasets table 2](images/datasets_table2.png)


### Dataset J16

A colony of 16 clonal raider ants (*O. biroi*) tracked for 24 hours. This videos in this example are taken with a simple webcam, using low resolution (10 pix/mm, representing around 5x5 pixels per tag).

![J16](images/J16.png)

### Dataset A36

A colony of 36 clonal raider ants (*O. biroi*) tracked for 24 hours. The videos in this example are taken with FLIR-Flea3 12MP camera at a relatively high resolution (25 pix/mm), enabling the program to distinguish between more colors, and to augment the basic tracking with pose tracking using DeepLabCut (see anTraX publication).

![A36](images/A36.png)


### Dataset V25

A colony of 25 clonal raider ants (*O. biroi*) tracked for 6 hours. The videos in this example are taken with a simple webcam, but using higher resolution than in dataset J16. Although the resolution is high, the image quality in this example is reduced by an acrylic cover placed between the ants and the camera. This dataset is an example for tracking an open boundry arena, where ants can leave and enter through a specific part of the boundry.

![V25](images/V25.png)


### Dataset G6X16

This dataset is an example for tracking multiple colonies within the same video (6 colonies of 16 *O. biroi* ants). In addition, the image quality in this example is reduced by the low contrast between the ants and the background, the petri dish covers, and light reflections from those covers.

![G6X16](images/G6X16.png)


### Dataset T10

A colony of 10 *Temnothorax nylanderi* ants recorded for 6 hours using a webcam. The ants are marked with 4 tags each. 

![T10](images/T10.png)


### Dataset C12

A colony of 12 *Camponotous fellah* ants recorded for 6 hours using a webcam. The ants are marked with 3 tags each. Although the ants are big in this example, their high velocity compared to the frame rate poses a challange and leads to an increased rate of linking errors. In addition, two ants are marked with the same color combination, but the algorithm was able to individually track them by using slight appearance variations. 

![C12](images/C12a.png)


### Dataset C32

A colony of 12 *Camponotous* ants (unidentified species) recorded for 24 hours. The ants are marked with 3 tags each. One ant (a virgin queen) is unmarked and identified by the classifier by other appearance markers (size, wings).

![C32](images/C32.png)


### Dataset D7

A group of 7 fruit flies marked by one tag each and recorded for 1 hour.

![D7](images/D7.png)


### Dataset D16

A group of 16 fruit flies marked by two tags each and recorded for 5 hours.

![D16](images/D16.png)


In this tutorial, we will track a small example dataset, included with anTraX. The dataset consists of a thirty minute recording of a colony of 16  *Ooceraea biroi* ants, split into 6 video files, each of  5 minutes duration. This dataset is a short segment of the longer [J16 benchmark dataset](datasets.md#dataset-j16), which is also available for dowload together with all the other benchmark datasets. Unlike these larger datasets, the JS16 dataset is appropriate to track on a laptop/desktop computer in a reasonable time.

All commands in the tutorial are to be entered in the bash terminal of your system, in the same environment anTraX was installed into.

### Download the test dataset

```console
git clone http://github.com/Social-Evolution-and-Behavior/anTraX-data.git
```

Inside the repository there is a directory called 'JS16' (the experimental directory). For a full  explanation of the structure of the experimental directory, refer to the [Preparing data for anTraX](data_organization.md) section.

### Open the antrax app

The dataset includes a pre-configured anTraX session. To explore and change the parameters, open the anTraX configuration app:

```console
antrax configure <path-to-JS16>
```

For a full  explanation of the configuration process and the tracking parameters, refer to the [Configuring a tracking session](configuration.md) section. 

### Track

The first step is the tracking. To run it, enter in the terminal:

```console
antrax track <path-to-JS16> --nw 3
```

The `--nw 3` option tells anTraX to run 3 parallel tracking threads. For a full listing of the options for the track command, refer to the [Runnning the tracking](tracking.md) section.  

### Train a classifier

Once tracking is complete, the next step is to train a blob classifier. The dataset includes an already trained classifier and a set of examples. But for the purpose of this example, let's run an additional round of training with 3  epochs:

```console
antrax train <path-to-JS16>/antrax_demo/classifier --ne 3
```

For a full explanation of the training step, including how to generate a training dataset, see the [Classifying tracklets](classification.md) section.  

### Classify tracklets

Once the classifier is trained, we can classify all the tracklets in the experiment:

```console
antrax classify <path-to-JS16>
```
For a full  listing of the options for the classify command, refer to the [Classifying tracklets](classification.md#classifying-tracklets)  section.  

### Run graph propagation 

The final step of the algorithm is running the graph propagation, or the 'solve' step:

```console
antrax solve <path-to-JS16>
```

For a full explanation of this step and the all the command options, refer to the [Graph propagation](propagation.md) section. 

#### Validate tracking

Now that tracking is complete, we can verify its accuracy and estimate the tracking error:

```console
antrax validate <path-to-JS16>
```
See the [Validating tracking results](validation.md) section and the anTraX publication for explanation about the validation process.

### Open graph-explorer

To debug the tracking, manually fix an important point in the experiment, or just view the tracklet graph, use the graph-explorer app:

```console
antrax graph-explorer <path-to-JS16> 
```

See [Using the graph explorer to view and debug ID assignments](propagation.md#using-the-graph-explorer-to-view-and-debug-id-assignments) for details about using the app.

### Loading and analyzing tracks

To load and analyze the tracking results, see the  [Working with tracking results - python](analysis_nb.ipynb) and  [Working with tracking results - matlab](analysis_matlab.ipynb) pages.

### Installation

In principle, anTraX installation on an HPC environment is the same as installation on any other Linux machine. The main difference is that typically, you will not have administrator privileges to intall system-wide packages on the HPC. Luckily, there are not many of those required by anTraX. You will also need to install MATLAB Runtime for version 2019a. You **do not** need to install MATLAB engine for python.

We recommend using a conda environemnt to setup anTraX on HPC, as it enables installation of required system packages such as [ffmpeg](https://anaconda.org/conda-forge/ffmpeg). If some packages are still missing and are not available in conda, work with your system administrator to find a solution.

If you plan on using DeepLabCut, install it into the python environment, and set it to 'light mode'  in  the bash profile file:

```bash
export DLClight=True 
```

### anTraX workflow on HPC

As computer clusters do not typically support interactive work, this will need to be done on a PC. An example for a tracking workflow using a computer cluster will be as follows:

1. Prepare the experimental directory/ies on a PC.
2. Configure a tracking session for each experimental directory on the PC. 
3. Sync the experimental directories into the HPC environment.
4. Run batch tracking on the HPC (see below).
5. If you have a trained blob classifier, jump to step 7. Otherwise, sync tracking results back to the PC, and create a training set using the interactive interface. 
6. Sync again to the HPC, and train the classifier. 
7. Run  the `classify` and `solve` commands in batch mode on the HPC.
8. Sync the results back to your PC.

It is recommended to use an incremental tool such as `rsync` to speed up data transfer. 

### Batch run on HPC environment 

If you are on an HPC environment, using the SLURM workload manager, you can run each of the batch commands (`track`, `classify`, `solve` and `dlc`) using the `--hpc` flag:

```console
antrax <command> <experiments> --hpc [--dry] [--hpc-options <opts>]
```

anTraX will then submit a SLURM job for each experiment, each containing a task for each video in the experiment. 

The optional `--dry` flag will create a SLURM jobfile, but will not submit it. It is useful to make changes to the sbatch options not currently supported by the anTraX interface.

The optional `--hpc-options` argument can control some of the SLURM options and accepts a comma separated list of some of the following options:

`throttle=<throttle>`

Number of tasks to run in parallel for each job.

`partition=<partition>`

The partition to run in (otherwise, use the system default).

`email=<email>`

Send an email for start/end of each job.

`cpus=<n>`

Allocate a specific number of cpus per task. The default value will vary according to the command (2 for tracking, 4 for classification/propagation/dlc, 12 for training).

`time=<time>`

Allocate time for task. Argument is a time string in the format supported by the `sbatch --time <time>` command (see [here](https://slurm.schedmd.com/sbatch.html)).

`mem-per-cpu=<mem>`

Allocate memory per task CPU. Argument in the format supported by the `sbatch --mem-per-cpu <mem>` command (see [here](https://slurm.schedmd.com/sbatch.html)).


**Note:** the commands `export-jaaba` and `run-jaaba`  do not currently support hpc mode.### What is a session?

A tracking *session* is a run of the algorithm with a set of settings and parameters. In the typical case, you will only create one session per experiment. However, it is sometimes useful to play around with a different parameter set without overwriting existing results, or track different parts of the experiment with different parameter sets. In these cases, multiple sessions should be created. The session, together with its parameters and results, is stored as a subdirectory of the experimental directory and is named by the session identifier name.

### Launching the anTraX app

To create and configure a tracking session, simply launch the anTraX app by entering the command into a bash terminal (don't forget to activate your virtual/conda environment if are using one):

```console
antrax configure [expdir]
```

The optional argument `expdir` is a full path to the experimental directory to be configured. If omitted, a directory selection dialog box will appear to select the experiment. You can move between experiments by using the options in the `Experiment` menu.

Any configuration changes are saved on the fly. When finished, just exit the app and the session will be saved. 


### Creating/loading a tracking session

If the experiment contains a previously  configured session, it will automatically load. Otherwise, you will be prompted to create a new one. Once a session is loaded/created, the configuration workflow will appear as tabs in the application window (see images below).

You can move between sessions, or create new ones, by using the options in the `Session` menu.

### Displaying video frames
The anTraX application window is divided into two main parts: configuration panel on the left, and the frame viewer on the right. The configuration panel contains multiple tabs corresponding to the algorithm step. The displayed image will be augmented according to the configuration tab currently active.
The frames in the experiment can be browsed using the selectors on the top part of the configuration panel (outlined in purple), which will appear in most of the configuration tabs. A frame in the experiment can be defined either by its video index (the ***Movie number*** selector)  and the frame index in that video (the ***Movie frame*** selector), or by its total index in the experiment (the ***Absolute frame*** selector).

![Frame display selection](images/frame_selection.png)

### Creating a background image

The first step in the configuration process is to generate a background image.

Use the ***method*** dropdown to select between the possible background computation methods. The ***median*** method computes the background as the per-pixel per-channel median of a set of randomly selected frames. The ***max*** method computes the background as the per-pixel per-channel max value. Select the number of frames for generating a background image. Obviously the more frames that are used, the better the background frame is, especially when median method is used. However, 20 will usually give a good trade-off between computation time and quality. Frames are randomly selected from the frame range.

***One BG*** option will generate a single background to be used throughout the tracking. With this option, you can select a frame range for selecting frames (useful for cases where some parts of the experiment are more suitable for background calculation).

***BG per subdir*** will generate a separate background frame for each subdirectory of videos. This option is useful for cases where there are movements between the location of the arena in the frame between recording sessions, or some other change in filming conditions.

The ***Create BG*** button will start the background creation process. Depending on the parameters, this might take several minutes. After the computation is done, the new background will be displayed. If several backgrounds are created, you can choose which one is displayed from the BG file dropdown menu.

The background images are saved as png files in the directory `expdir/session/parameters/background/`.


![Create background tab](images/background_creation.png)

### Setting the spatial scale
In the second tab, the spatial scale of the videos will be defined. This is required to have all the parameters and results in real world units, which is essential for  parameters to be generalized between experiments and tracking results comparable between experiments.

To set the scale, choose a feature in the image of which the dimensions are known. Choose the appropriate tool from the drop down menu (either ***Circle*** or ***Line***), press the ***Draw*** button, and adjust the tool to fit the feature.

When done, enter the Length/Diameter of the feature in *mm* in the box, and finish by pressing the ***Done*** button.

![Scale tab](images/scale.png)

### Creating an ROI mask

The ***ROI Mask*** (Region Of Interest) is used to define the regions of the image in which tracking is performed. 

To set a mask, start by either a white mask ("track everywhere") by pressing the ***Reset to White*** or a black mask ("track nowhere") by pressing the ***Reset to Black***. Then, add and remove regions by selecting a tool from the dropdown and drawing on the image. Adjust by dragging the anchor points. When done, double click the tool. You can repeat this process untill the ROI is ready.

The ROI mask is saved as png files in the directory: `expdir/session/parameters/masks/`.

![ROI mask](images/roi-mask.png)

### Multi-colony experiments

The ***Multi-Colony*** option is used to control how a mask with several disconnected ROIs should be treated if these regions correspond to separate ant colonies. If checked, each of these regions will be treated as a separate colony, containing a full and fixed set of identified ants, and will be saved separately. Use the dropdown to control the numbering order of the ROIs, and the Assign colony labels buttons to manually assign labels to each numbered colony (avoid white spaces in the labels).

The colony masks are saved as png files in the directory: `expdir/session/parameters/masks/`.

![Multi-colony experiment](images/multi-colony.png)

### Open boundary ROIs

The ***Open Boundary*** option is used to mark parts of the ROI perimeter that are "open" to ants getting in and out of the ROI. This is used to optimize tracking in these regions. Otherwise, the ROI is assumed to be completely closed. To mark a segment of the boundary as 'open', click ***Add*** and adjust the shape, so its intersection with the ROI boundary will be the open region. Double-click to finish. The 'open' segment will be marked with thick blue line.

![Open boundry ROI](images/open-boundry.png)

### Tuning the segmentation

After subtracting the image from the background, anTrax segments the image into foreground and background, with the foreground being composed of several connected components ('blobs'). This is a multi-parameter process that should be tuned for each experiment. 

On the ***Segmentation*** tab, several of the segmentation parameters can be tuned, while displaying the results. The control parameters include:

* ***Segmentation threshold:*** in units of gray value difference between image and background.
	
* ***Adaptive threshold:*** if checked, the threshold will be adjusted locally as a function of the background brightness, causing the segmentation to be more sensitive in darker areas.
	
* ***Min area (pixels):*** Blobs below this threshold are discarded.
	
* ***Closing (pixels):*** Optional morphological closing, that merge blobs separated by few pixels.
	
* ***Opening (pixels):*** Optional morphological opening, eroding thin pixel lines.
Min intensity (gray level). Blobs with maximum intensity lower than this value will be discarded.

* ***Convex hull:*** Fill in the blob convex hull. Useful in cases where there is bad contrast between parts of the ant and the background.
	
* ***Fill holes:*** Useful when very bright tags are used, that do not have good contrast with the background and appear as 'holes' in the blob.
	
* ***Min intensity:*** Optional blob filter, which discard blobs with maximal intensity lower than the threshold value.

The display of the segmented frame can be configured usng the checkboxes below the frame selectors: ROI mask can be turned on/off, blobs can be shown as convex hull curves (default) or as colored segmented regions (better to check the fine details of the segmentation). Text showing the blob area in pixels and maximum intensity value can be displayed. 

![Image segmentation](images/segmentation.png)

### Tuning single individual size range

anTrax uses the size of an individual ant for filtering possible single ant tracklets for classification and for calibrating the linking algorithm. The single ant size is defined by the possible size range, which is adjusted in the ***single ant*** tab. For tuning these range parameters, the blobs detected in the displayed frames are marked with green outlines if they are in the single ant range, with red if they are larger, and with pink if they are smaller. It is recommended to scan a decent number of frames throughout the experiment to look for near-threshold cases. Note, the range doesn't need to perfectly classify blobs; rather, it captures the possible size range for single ants. If your experiment contains individuals with variable size, choose the range to capture all the individuals, even at the 'price' of classifying some multi-animal blobs as single-animal.

![Single ant size range](images/single_ants.png)


### Tuning the linking step

Linking is the process of connecting blobs from consecutive frames into tracklets. Linked blobs are assumed to represent a case where some or all of the ants that are included in the blob from the first frame are also included in the blob in the second frame. A blob can be linked to zero, one, or multiple blobs in the the other frame. 

As described in the paper, anTraX uses optical flow to link blobs. This is used whenever a blob has more than one possible blob to link to inside its "linking cluster".

The displayed image can be selected as the "previous frame", "current frame", or "blend" (an overlay of the two frames in different color channels). It is also possible to display the linking clusters by selecting the checkbox. Clusters with more than one blob in one of the frames (which will undergo optical flow) will be marked with bright blue contours, while other clusters will be marked with dim blue contours.

* **max speed**: This is the typical maximum velocity possible by the tracked individual animals. It is used to calculate the max possible distance for an animal to move between frames (taking into account the specific time interval). 

* **linking cluster coefficient**: The coefficient multiplies the max possible distance to get a conservative radius for the linking cluster (this parameter should not be changed regularly). 

* **optical flow cutoff coefficient**: When optical flow is used, a flow index will be computed between each blob pair. The pair will be linked if the index is above a cutoff threshold. The cutoff threshold is set by the middle of the single animal range set in the previous tab, times this coefficient. Generally, the linking process is relatively robust to the precise value of this parameter. However, in some cases, increasing it might lead to an increase in false negative linking errors (missed true links), and decreasing it might lead to increase in false positive linking errors (wrong links). 

![Linking tunning tab](images/linking.png)

### Classification options

The **Classification** tab is used to configure how blob images are processed for classification. Note that the processing controlled by the options in this tab does not effect blob segmentation and linking.

* The ***Color correction*** option applies a white-reference correction to the frame before extracting the images for classification (see the anTraX publication, appendix section 2.1). Applying color correction can help in cases of inhomogeneous lighting, that creates variation in tag colors in different regions of the image. It also help to reduce inter-experiment color variability.. However, in cases of poor contrast, color correction might actually degrade the image quality, so it is recommended to look closely at the effect of these option using on your images using the GUI before applying it.
* Blob images are passed to the classifier after masking, using the segmentation mask generated during tracking. Sometimes, however, the optimal segmentation for tracking is not the optimal segmentation for classification. For example, in cases where some tag colors have poor contrast with the background, they are cut from the blob. The ***Blob mask dilation*** option apply a dilation operation to enlarge the blob before extracting its image for classification. Using the **Show segmentation** check box will visualize the added blob area.



![Classification tab](images/BlobMaskDilation1.png)

### Entering individual tags information 

The **IDs** tab is used to configure the list of ant IDs used in the experiment. First, set the tagging type as either *untagged* (experiments without color tags; no classification will be done), *group-tagged* (non individual tags; tracklet classification will be done, but no graph propagation) or *individual-tagged* (ants marked with unique IDs). 

Before classification, you will need to provide the program a list of the ants in the experiment (identified by their color tags). In case your classifier is trained to identify other types of objects (food, brood, prey insect, etc.) you will need to provide these as well.

The list of labels must match the the one the classifier is trained with (read more about [classifiers](classification.md)).

If your experiment is a multi-colony one, it is assumed the ID list is the same for all colonies in the experiment. If it is not the case, give a list that include all possible IDs, and adjust it using a config file as described [below](configuration.md#modifying-the-id-list-using-temporal-config-file).

The label list is defined by the file `expdir/session/parameters/labels.csv`. Each row in the file contains two entries. The first is the label ID, and the second is the category. Three categories exist: ant_labels, noant_labels, and other_labels. The list must include the label 'Unknown' in the 'other_labels' category. If the animal size variability in the experiment is considerable, and many multi-individual blobs are classified as single-animal (see discussion [above](configuration.md#tuning-single-individual-size-range)), it is recommended to include also a 'Multi' class in the 'other_labels' category.  

![Labels file](images/labels.png)

The ***IDs*** tab is an easy way to configure the list of labels for cases where animals are marked with two color tags. First, check the boxes of the color tags used. A label list containing all possible combinations will be created. Next, trim the list to only include the combinations actually used. Finally, add no-ant labels as needed. 

### Tuning the graph propagation step

The **Propagation** tab is used to configure the graph propagation step. 

Propation is done in parallel on movie groups. You can select how movies are grouped using the ***Group by*** drop down list. The ***movie*** option will processed separately. This is the fastest option, but will be less optimal near the start/end of the movie. It is appropriate for either long movies or when movies are not continuous in time. The ***subdir*** option will group movies according to the subdirectory organization (see the [data organization page ](data_organization.md)). The ***experiment*** option will group all movies in the experiment together. This is the slowest option, and is recommended only for short experiements (less than 24 hours). The ***custom*** option allows the user to define custom movie groups.

The groups are enumerated in sequence from 1 to the number of groups. This enumeration is used in the [batch run](propagation.md#propagating-ids-on-tracklet-graphs) to identify the group.

The **Pairs search depth** parameter controls the graph radius over which the algorithm searches for topological propagation opportunities (see the [anTraX paper]). A high number might give better results, but will overhead the algorithm. A value between 5-10 usually gives the best tradeoff.

The **Max iteration** parameter imposes a limit on the number of total iterations the propagation algorithm performs. In most cases, the algorithm converges after a few iterations. For the rare cases where it does not, a value of 10 usually represents a good tradeoff.

The **Impossible speed** parameter is used to filter out cases where the algorithm assigns IDs to tracklets that represent an impossible traveling speed for an ant. Use a very conservative value here, as many times there is a lag in the movie (as a results of skipping frames in the recording) that is not captured correctly by the interval frame data. A value of at least twice the max speed is recommended. 

The **temporal config** option tells the program to use the commands in the temporal config file (see below).

The **manual config** option tells the program to use manual ID assigments (see [later in the documentation](propagation.md#using-the-graph-explorer-to-view-and-debug-id-assigments)).

### Modifying the ID list using temporal config file

Optionally, a configuration file can be written for adjusting the ID list per colony or per time. Currently, the config supports ***remove** commands. The file should be text file located in `expdir/session/parameters/ids.cfg`. Each line in the file is interpreted as a command in the format:

```console
command colony id from to
```

The time arguments `from` and `to` can be either `start` for the first frame in the experiment, `end` for the last frame in the experiment, `m` followed by a number for the first frame in a movie (e.g. `m3`), 'f' followed by a number for a specific frame in the experiment (e.g. `f4000`), or a combination of a movie and frame for a specific frame in a specific movie (e.g. `m9f1000`).

To remove the id GP in colony C1 for the entire experiment:

```console
remove C1 GP start end
```

To remove YY in colony C5 from movie 22 to the end:

```console
remove C5 YY m22 end
```

To remove BG in all colonies for frames 20000 to 30000:

```console
remove all BG f20000 f30000
```

To remove PP in colony A from frame 100 in movie 4 to frame 2000 in movie 7:

```console
remove A PP m4f100 m7f2000
```

### The 'other options' tab


TBA