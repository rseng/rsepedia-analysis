# pymicmac

[![Build Status](https://travis-ci.org/ImproPhoto/pymicmac.svg?branch=master)](https://travis-ci.org/ImproPhoto/pymicmac)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/cedd804840ca4de0af4c6bae6939b28d)](https://www.codacy.com/app/ImproPhoto/pymicmac?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ImproPhoto/pymicmac&amp;utm_campaign=Badge_Grade)
[![Anaconda-Server Badge](https://anaconda.org/improphoto/pymicmac/badges/installer/conda.svg)](https://conda.anaconda.org/improphoto)
[![DOI](https://zenodo.org/badge/57877195.svg)](https://zenodo.org/badge/latestdoi/57877195)

This software is a result of the Netherlands eScience Center project [Improving Open-Source Photogrammetric Workflows for Processing Big Datasets](https://www.esciencecenter.nl/project/improving-open-source-photogrammetric-workflows-for-processing-big-datasets).

`pymicmac` provides a Python interface for `MicMac` workflows execution and distributed computing tools for `MicMac`.
[`MicMac`](http://micmac.ensg.eu) is a photogrammetric suite which contains many different tools to execute photogrammetric workflows.

In short, a photogrammetric workflow contains at least:

 - **Tie-points detection.** First,  key  features  in  the  images  are  extracted.   This  is  done
for example with the [SIFT algorithm](https://www.cs.ubc.ca/~lowe/papers/ijcv04.pdf).  Second, the key features are cross-matched
between different images to detect points corresponding to the same physical
locations and that are visible in different images.  The detected points are called tie-points
(they are also referred to as homologous points in related literature).

 - **Bundle block adjustment.** The camera positions and orientations are estimated and the parameters calibrated.

 - **Dense image matching.** The detected tie-points are matched and 3D projected to produce the dense point cloud. The 3D points are back projected in the images to correct for projective deformation. This creates a metrically correct True-orthophoto.

The MicMac suite contains several tools dedicated to each of the steps of the programmatic worflow. The tie-point detection is done with `Tapioca`, the bundle block adjustment is done with `Tapas` and the dense matching and point cloud generation is done with `Malt`, `Tawny` and `Nuage2Ply`.

 `pymicmac` addresses two main issues with `MicMac`:

1 `pymicmac` helps you when running a sequence of `MicMac` commands. The sequence of commands is defined in an XML. During execution `pymicmac` creates isolated execution folders separating the input data from the intermediates and the output data. `pymicmac` also adds CPU/MEM monitoring for the commands.

2 `pymicmac` contains distributed computing versions of `Tapioca` and of the matching pipeline (`Malt`, `Tawny`, `Nuage2Ply`). This allows running `Tapioca` and the matching pipeline in distributed systems such as SGE clusters or a bunch of ssh-reachable machines.

The `micmac-run-workflow` tool uses the sequential commands execution tool of [`pycoeman`](https://github.com/NLeSC/pycoeman). Detailed information is provided in the [Instructions](#instructions) section.

In section [Large image sets](#large-image-sets) we provide some tips on how to use `MicMac` and `pymicmac` for processing large image sets using distributed computing for the tie-points extraction and the dense image matching, and tie-points reduction for the bundle block adjustment.

A step-by-step tutorial is also available in [Tutorial](https://github.com/ImproPhoto/pymicmac/tree/master/docs/TUTORIAL.md).

## Installation

The easiest way to install pymicmac is by using the [Anaconda package](https://anaconda.org/ImproPhoto/pymicmac).

```bash
conda install -c improphoto -c conda-forge pymicmac
```

## Development installation

Clone this repository and install it with pip (using a virtualenv is recommended):

```
git clone https://github.com/ImproPhoto/pymicmac
cd pymicmac
pip install -e .
```

Python dependencies: `pycoeman` and `noodles` (see https://github.com/NLeSC/pycoeman and https://github.com/NLeSC/noodles for installation instructions)

Other Python  dependencies (numpy, tabulate, matplotlib, lxml) are automatically installed by `pip install -e .` but some system libraries have to be installed (for example freetype is required by matplotlib and may need to be installed by the system admin)

For now `pymicmac` works only in Linux systems. Requires Python 3.5.

## Instructions

The tool `micmac-run-workflow` is used to execute entire photogrammetric workflows with MicMac or portions of it. We recommend splitting the workflow in three pieces: (1) tie-points extraction, (2) bundle block adjustment and (3) dense image matching. Each time the tool is executed, it creates an independent execution folder to isolate the processing from the input data. The tool can be executed as a python script (see example in `tests/run_workflow_test.sh`) or can be imported as a python module. Which MicMac commands are executed is specified with a XML configuration file.

### Workflow XML configuration file

The Workflow XML configuration file format is the sequential commands XML configuration file format used by [`pycoeman`](https://github.com/NLeSC/pycoeman). For `pymicmac`, usually the first tool in any Workflow XML configuration file links to the list of images. So, we can use `<requirelist>` to specify a file with the list of images. Next, some XML examples:

- tie-points extraction:
```
<SeqCommands>
  <Component>
    <id> Tapioca </id>
    <requirelist>list.txt</requirelist>
    <command> mm3d Tapioca All ".*jpg" -1 </command>
  </Component>
</SeqCommands>
```

- Bundle block adjustment
```
<SeqCommands>
  <Component>
    <id> Tapas </id>
    <requirelist>list.txt</requirelist>
    <require>tie-point-detection/Homol</require>
    <command> mm3d Tapas Fraser ".*jpg" Out=TapasOut </command>
  </Component>
</SeqCommands>
```

- Dense image matching:
```
<SeqCommands>
  <Component>
    <id> Malt </id>
    <requirelist>list.txt</requirelist>
    <require> param-estimation/Ori-TapasOut</require>
    <command> mm3d Malt GeomImage ".*jpg" TapasOut "Master=1.jpg" "DirMEC=Results" UseTA=1 ZoomF=1 ZoomI=32 Purge=true </command>
  </Component>
  <Component>
    <id> Nuage2Ply </id>
    <command> mm3d Nuage2Ply "./Results/NuageImProf_STD-MALT_Etape_8.xml" Attr="1.jpg" Out=1.ply</command>
  </Component>
</SeqCommands>
```

Following the examples above, we could execute a whole photogrammetric workflow with:
```
micmac-run-workflow -d /path/to/data -e tie-point-detection -c tie-point-detection.xml
micmac-run-workflow -d /path/to/data -e param-estimation -c param-estimation.xml
micmac-run-workflow -d /path/to/data -e matching -c matching.xml
```

*NOTE*: all file and folder names, specified in `<require>` and `<requirelist>` must be provided relative to the folder where all the data is - `/path/to/data`.

### Monitoring

The tool used by `pymicmac` to run the commands, `pycoeman` keeps the log produced by each command in a .log file. Additionally, it stores a .mon and a .mon.disk with the CPU/MEM/disk monitoring. There are `pycoeman` tools for calculating statistics and creating plots for monitoring (see https://github.com/NLeSC/pycoeman). There are some tools in `pymicmac` to analyze the logs of several `MicMac` commands: `micmac-tapas-log-anal`, `micmac-redtiep-log-anal`, `micmac-redtiep-log-anal`, `micmac-campari-log-anal`, `micmac-gcpbascule-log-anal`, `micmac-gcpbascule-log-plot`, `micmac-campari-log-plot` (the MicMac command if part of the tool name).

## Large image sets

For the tie-points extraction and dense image matching the processing can be easily enhanced by using distributed computing (clusters or clouds). The reason is that the processes involved can be easily split in independent chunks (in each chunk one or more images are processed). For the bundle block adjustment, this is not the case since the involved processes usually require having data from all the images simultaneously in memory. In this case, we propose to use tie-points reduction to deal with large image sets.

For more information about distributed computing and tie-points reduction, see our paper: Martinez-Rubi, Oscar, Francesco Nex, Marc Pierrot-Deseilligny, and Ewelina Rupnik. “Improving FOSS Photogrammetric Workflows for Processing Large Image Datasets.” Open Geospatial Data, Software and Standards 2 (May 15, 2017): 12. [https://doi.org/10.1186/s40965-017-0024-5.](https://doi.org/10.1186/s40965-017-0024-5).

### Distributed computing

Some steps of the photogrammetric workflow, namely the tie-points extraction and the dense image matching, can be executed more efficiently on distributed computing systems exploiting the innate data parallelism in photogrametry.

For example, the `Tapioca` tool (tie-points extraction) first extracts the features for each image and then cross-matches the features between image pairs. The proposed distributed computing solution divides the list of all image pairs in chunks where each chunk can be processed mostly independently. At the end of the results from each chunk processing need to be combined.

We use the parallel commands execution tools of `pycoeman`. The various parallel/distributed commands are specified in an XML configuration file which is similar to the Workflow XML configuration file. An example XML configuration, where the `Tapioca` processes two data chunks, each containing half of the image pairs:

```
<ParCommands>
  <Component>
    <id>0_Tapioca</id>
    <requirelist>DistributedTapioca/0_ImagePairs.xml.list</requirelist>
    <require>DistributedTapioca/0_ImagePairs.xml</require>
    <command>mm3d Tapioca File 0_ImagePairs.xml -1</command>
    <output>Homol</output>
  </Component>
  <Component>
    <id>1_Tapioca</id>
    <requirelist>DistributedTapioca/1_ImagePairs.xml.list</requirelist>
    <require>DistributedTapioca/1_ImagePairs.xml</require>
    <command>mm3d Tapioca File 1_ImagePairs.xml -1</command>
    <output>Homol</output>
  </Component>
</ParCommands>
```

The `pycoeman` tool `coeman-mon-plot-cpu-mem` can be used to get a plot of the aggregated CPU and MEM usage.

#### Distributed MicMac Tools

##### Tapioca

The tools `micmac-disttapioca-create-pairs`, `micmac-disttapioca-create-config`, `micmac-disttapioca-combine` together with the parallel commands execution tools of `pycoeman` are used to run `Tapioca` on distributed computing systems.

In order to run `Tapioca` on a distributed system, a list containing image pairs is needed. The `micmac-disttapioca-create-pairs` tool creates such a list. The `micmac-disttapioca-create-config` tool splits the image pairs XML file in multiple chunks and creates an XML configuration file compatible with `pycoeman`:

```
micmac-disttapioca-create-config -i [input XML image pairs] -f [folder for output XMLs and file lists, one for each chunk] -n [number of image pairs per output XML, must be even number] -o [XML configuration file]
```

Next, the distributed tool can be executed on any hardware system supporting `pycoeman` (see https://github.com/NLeSC/pycoeman) using the `coeman-par-local`, `coeman-par-ssh` or `coeman-par-sge` tools.

After the distributed `Tapioca` has finished, the outputs from the different chunks need to be combined. The `micmac-disttapioca-combine` tool combines all outputs in a final Homol folder:

```
micmac-disttapioca-combine -i [folder with subfolders, each subfolder with the results of the processing of a chunk] -o [output combined folder]
```

##### Matching (Malt, Tawny, Nuage2Ply)

To generate the final point cloud on a distributed computing system, the dense point matching is parallelized by `micmac-distmatching-create-config` and the parallel commands execution tool of `pycoeman`.

The algorithm in `micmac-distmatching-create-config` is restricted to aerial images when the camera orientation obtained in the parameter estimation step of the photogrammetric workflow is in a cartographic reference system. From the estimated camera positions and assuming that the Z direction (along which the pictures were taken) is always pointing to the ground, the tool computes the XY bounding box that includes all the XY camera positions. The bounding box is divided in tiles like shown in the figure below:

![exampledistmatching](docs/distmatching_example.png)

Each tile can be then processed by an independent process. For each tile the images intersecting their XY position with the tile are used. If needed this set of images is extended by the nearest neighbours to guarantee a minimum of 6 images per tile.

The `micmac-distmatching-create-config` generates the tiles from the initial large list of imagesand creates an XML configuration file suitable for `pycoeman`:

```
micmac-distmatching-create-config -i [orientation folder] -t [Homol folder] -e [images format] -o [XML configuration file] -f [folder for extra tiles information] -n [number of tiles in x and y]
```

Next, the distributed tool can be executed on any hardware systems supported by `pycoeman` (see https://github.com/NLeSC/pycoeman). The  `coeman-par-local`, `coeman-par-ssh` or `coeman-par-sge` tools are used for this purpose.

### Tie-points reduction

For more efficient workflow execution, we propose to perform a tie-points reduction step before the bundle adjustment. This extra step is added to the processing chain of the parameters estimation step. For detailed explanation, please refer to the Tie-point reduction section in this [report](http://knowledge.esciencecenter.nl/content/report-tie-points.pdf).

There are two tools for this purpose: `RedTieP` and `OriRedTieP`. The former tool should be preceeded by `NO_AllOri2Im` and  `Martini` should preceed the latter. For examples, see `tests/param-estimation_reduction.xml` and  `tests/param-estimation_orireduction.xml`.

Note that after running the tie-points reduction tools, the Homol folder has to be changed (see the examples).
Also note that when running `RedTieP`, it is possible to use parallel execution mode together with the tool `micmac-noodles`. See the example in `tests/param-estimation_reduction.xml`.

The`micmac-homol-compare` tool can be used to compute the reduction factors.

Contributions are welcome!
Please create a pull request, and make sure that:

1. we follow the [GitHub Flow](https://guides.github.com/introduction/flow/) branching model.
2. For other development and coding style conventions, see the [NLeSC Style Guide](https://nlesc.gitbooks.io/guide/content/).
4. Don't include extra dependencies without a good reason. Only use licenses compattible with the license of this project- Apache v2.0.
5. Please, document your code, and provide unit tests.

Make sure that after your contribution at least the [`run_workflow_test`](https://github.com/ImproPhoto/pymicmac/blob/master/tests/run_workflow_test.sh) passes.
# Tutorial

This tutorial is intended to drive you in the process to convert a bunch of images in some folder of your (Linux) computer into a colored dense point cloud using MicMac and pymicmac.

## Installation

Before anything else, we need to install all the required software: MicMac and pymicmac and their dependencies.
In this tutorial, we do the installation on Ubuntu 16.04. Note that some steps and libraries names may be different in other Linux distributions.

First, we install MicMac (we get MicMac from its SVN repository):
```
# Install mercurial to be able to download MicMac
sudo apt-get install mercurial

# We install it in the next location (please change accordingly to your system)
cd /home/oscar/sw

# Clone the repo (this requires user/password)
hg clone https://geoportail.forge.ign.fr/hg/culture3d

# Install MicMac
cd culture3d
mkdir build
cd build
cmake ..
make -j24
make install

# Assuming that we installed micmac in /home/oscar/sw/culture3d, add the next lines to your .bashrc (in your case replace accordingly)
export LD_LIBRARY_PATH="/home/oscar/sw/culture3d/lib:$LD_LIBRARY_PATH"
export PATH="/home/oscar/sw/culture3d/bin:$PATH"

# Source .bashrc to activate the installation
source ~/.bashrc
```

Second, pymicmac is a Python 3.5 package so we need to have a Python 3.5 installation. We recommend using Anaconda:
```
# Get the latest Anaconda installer (in 32 or 64 bits depending on your system)
wget https://repo.continuum.io/archive/Anaconda3-4.2.0-Linux-x86_64.sh

# Install it
bash Anaconda3-4.2.0-Linux-x86_64.sh

# Create a anaconda environment for all the installation of pymicmac
conda create --name python35 python=3.5

# Add this line in your .bashrc
source activate python35

# Source .bashrc to activate the installation
source ~/.bashrc
```

Third, pymicmac has some system library requirements (freetype, ssl, ffi) and also requires pycoeman and noodles:
```
# Install pycoeman dependencies
sudo apt-get install libfreetype6-dev libssl-dev libffi-dev

# Install pycoeman
pip install git+https://github.com/NLeSC/pycoeman

# Install noodles
pip install git+https://github.com/NLeSC/noodles
```

Finally, we install pymicmac:
```
# Install pymicmac
pip install git+https://github.com/ImproPhoto/pymicmac
```

We can test the installation:
```
mm3d -help
micmac-run-workflow -h
```

## Processing a dataset

Now that we have all performed the installation of all tne required software, we will generate a colored dense point cloud.
We assume that the data is in `/home/oscar/data/GRONAU/4ms_60m_1000`. Concretely the folder looks like:

```
ls /home/oscar/data/GRONAU/4ms_60m_1000
coord_List2D.xml gcp_List3D.xml GrapheHom.xml Ori-IniCal
IMG_0990.JPG ...
```

In addition to the set of JPG images, we also have GCPs files, a `GrapheHom.xml` file and a `Ori-IniCal` folder.
The GCPs file `gcp_List3D.xml` have the 3D positions of the Ground Control Points (GCPs) and Check Points (CPs), and their 2D positions in the images are registered in the `coord_List2D.xml` file. The `GrapheHom.xml` contains the list of valid image pairs (extracted from geotag info). The `Ori-IniCal` folder contains a XML file with the initial calibration information.

For pymicmac it will make our live easier if we create a file with the list of images:
```
cd /home/oscar/data/GRONAU/4ms_60m_1000
ls *JPG > images.list
```

Our photogrammetric pipeline with MicMac consists of running Tapioca to extract tie-points, followed by Tapas to perform the bundle block adjustment, and finally Malt, Tawny and Nuage2Ply to get the colored dense point cloud.

The commands to be executed are configured in pymicmac with XML. During execution, the commands are executed in a folder than the one where the original data is stored. For the required data links are created. pymicmac monitors the CPU/MEM and disk usage. With pymicmac we can configure any photogrammetric workflow and we can split it in parts in any way we want. In the next subsections, we present a couple of examples of different strategies to execute workflows.

### Single XML with entire workflow
In this example we define a single XML to execute the entire workflow. The pymicmac XML configuration file is called `Workflow.xml` and its content is:
```
<SeqCommands>
  <Component>
    <id> Tapioca </id>
    <require> GrapheHom.xml </require>
    <requirelist> images.list </requirelist>
    <command> mm3d Tapioca File GrapheHom.xml -1 </command>
  </Component>
  <Component>
    <id> Tapas </id>
    <command> mm3d Tapas Fraser ".*JPG" InCal=IniCal Out=TapasOut </command>
    <require> Ori-IniCal </require>
  </Component>
  <Component>
    <id> GCPBascule </id>
    <command> mm3d GCPBascule ".*JPG" TapasOut GCPBOut gcp_List3D.xml coord_List2D.xml </command>
    <require> gcp_List3D.xml coord_List2D.xml </require>
  </Component>
  <Component>
    <id> Malt </id>
    <command> mm3d Malt Ortho ".*JPG" GCPBOut </command>
  </Component>
  <Component>
    <id> Tawny </id>
    <command> mm3d Tawny Ortho-MEC-Malt </command>
  </Component>
  <Component>
    <id> Nuage2Ply </id>
    <command> mm3d Nuage2Ply MEC-Malt/NuageImProf_STD-MALT_Etape_8.xml Attr=Ortho-MEC-Malt/Orthophotomosaic.tif Out=pointcloud.ply </command>
  </Component>
</SeqCommands>
```

The workflow is executed with pymicmac with:
```
cd /home/oscar/data/GRONAU/4ms_60m_1000/
micmac-run-workflow -d . -c Workflow.xml -e WorkflowOutput
```

Note that in each command the required files/folders are specified with the tags `<require>` and `<requirelist>`. The specified locations are relative paths to the folder specified in the `-d` option of the `micmac-run-workflow` command, i.e. `.` which is `/home/oscar/data/GRONAU/4ms_60m_1000/`. Also note that there is no need to duplicate required items. For example, all the commands need the images (provided with `<requirelist>`) but it is enough to specify this in the first command.

Executing `micmac-run-workflow` will first create the `WorkflowOutput` folder, then it will create links inside this folder to all the required data (specified with `<require>` and `<requirelist>`) and finally will run all the commands. After the execution is finished, for each of the commands we will find in `WorkflowOutput` a log file, a mon file and and mon.disk file. These contain respectively the log of the command execution, the CPU/MEM usage and the disk usage.
pycoeman has tools to obtain statistics of the CPU/MEM usage:
```
coeman-mon-stats -t Tapioca,Tapas,GCPBascule,Malt,Tawny,Nuage2Ply -f WorkflowOutput
```


### Various XMLs to (re-)execute parts of the workflow
While the previous example is useful, one can argue if only to have CPU/MEM/disk usage and a clean environment it is worthy the hassle to install and learn how to operate pymicmac, and that is a completely valid point. However, it is not for cases like the previous one that pymicmac was made.

In this second example we will see when pymicmac is really beneficial. Now we want to run a workflow similar to the previous one but with a tie-points reduction step before Tapas. This will make Tapas execution faster. We also want to compare the workflow with tie-points reduction with the case when no reduction is used. More concretely we will look at the residuals of Tapas and the the errors of GCPBascule. We would also like to see the impact of Tapas if reduction is used (less memory and faster execution)

Running Tapioca is common in the two workflows we want to test. Thus, we divide the workflows in two parts. The first part is Tapioca and the second part is Tapas and GCPBascule with and without tie-points reduction. In this example, since we are only interested in Tapas and GCPBascule we will not run Malt, Tawny and Nuage2Ply.

The XML to run Tapioca is called `Tapioca.xml` and its content is:
```
<SeqCommands>
  <Component>
    <id> Tapioca </id>
    <require> GrapheHom.xml </require>
    <requirelist> images.list </requirelist>
    <command> mm3d Tapioca File GrapheHom.xml -1 </command>
  </Component>
</SeqCommands>
```

We execute it with:
```
cd /home/oscar/data/GRONAU/4ms_60m_1000/
micmac-run-workflow -d . -c Tapioca.xml -e TapiocaOutput
```

This will create the `TapiocaOutput` folder, make links to the images and to the `GrapheHom.xml` file and will execute Tapioca inside `TapiocaOutput`. When execution is finished, in addition to the `Homol` created by Tapioca, we will also have a log file, a mon file and a mon.disk. We can use `coeman-mon-stats` to get a overview on the CPU/MEM usage of Tapioca:
```
coeman-mon-stats -t Tapioca -f TapiocaOutput
```

Next, we define the rest of the workflow when no tie-points reduction is done. The `WorkflowNoTPR.xml` is:
```
<SeqCommands>
  <Component>
    <id> Tapas </id>
    <require> Ori-IniCal TapiocaOutput/Homol </require>
    <requirelist> images.list </requirelist>
    <command> mm3d Tapas Fraser ".*JPG" InCal=IniCal Out=TapasOut </command>
  </Component>
  <Component>
    <id> GCPBascule </id>
    <command> mm3d GCPBascule ".*JPG" TapasOut GCPBOut gcp_List3D.xml coord_List2D.xml </command>
    <require> gcp_List3D.xml coord_List2D.xml </require>
  </Component>
```

Note that in this case in Tapas we need to specify that we require `images.list` (we always do this in the first command in the XML) and that we also require the `Homol` folder. This will be inside the `TapiocaOutput` folder created before (which will be inside `/home/oscar/data/GRONAU/4ms_60m_1000/` so it is fine to use the relative path). We can now run the workflow with:
```
cd /home/oscar/data/GRONAU/4ms_60m_1000
micmac-run-workflow -d . -c WorkflowNoTPR.xml -e NoTPROutput
```

This will create the `NoTPROutput` folder, then it will create links to the images, the `Homol` folder and the rest of files, and finally will run Tapas and GCPBascule. Like before, after the execution is finished we will find log files, mon files and mon.disk files also in `NoTPROutput` folder.

Following we define the workflow with tie-points reduction. We can reuse the same `Homol`folder than before so there is not need to rerun Tapioca. The `WorkflowTPR.xml` is:
```
<SeqCommands>
  <Component>
    <id> NO_AllOri2Im </id>
    <require> TapiocaOutput/Homol </require>
    <requirelist> images.list </requirelist>
    <command> mm3d TestLib NO_AllOri2Im ".*JPG" Quick=1 </command>
  </Component>
  <Component>
    <id> RedTieP </id>
    <command> mm3d RedTiep ".*JPG" NumPointsX=12 NumPointsY=12  WeightAccGain=0.00; rm Homol; mv Homol-Red Homol </command>
  </Component>
  <Component>
    <id> Tapas </id>
    <require> Ori-IniCal </require>
    <command> mm3d Tapas Fraser ".*JPG" InCal=IniCal Out=TapasOut </command>
  </Component>
  <Component>
    <id> GCPBascule </id>
    <command> mm3d GCPBascule ".*JPG" TapasOut GCPBOut gcp_List3D.xml coord_List2D.xml </command>
    <require> gcp_List3D.xml coord_List2D.xml </require>
  </Component>
```

First NO_AllOri2Im will run (which is required by RedTiep, the tie-points reduction tool) and since this is the first command we add the requires for the images and the `Homol` folder. Next, the actual reduction with RedTieP will be done. Note that after the reduction we will replace the `Homol` folder with the one output by the tool. But do no panic!, the `rm Homol` only deletes the link to the `Homol` folder. The full set of tie-points is still there. This is one of the benefits of separating the execution of the several parts of the workflow and of using links. Finally, we will run Tapas and GCPBascule exactly as before but in this case they will use a reduced set of tie-points. We execute it with:
```
cd /home/oscar/data/GRONAU/4ms_60m_1000
micmac-run-workflow -d . -c WorkflowTPR.xml -e TPROutput
```

After the execution is done, we will find log files, mon files and mon.disk files in the `TPROutput` folder.

#### Comparison of workflows

We want to compare the results in the two workflows. First, we see the different CPU/MEM usage of the executed commands in both cases:
```
coeman-mon-stats -t Tapioca,NO_AllOri2Im,RedTieP,Tapas,GCPBascule -f TapiocaOutput,NoTPROutput,TPROutput
```
The output will look something like:
```
##########################
Time/CPU/MEM tools monitor
##########################
#Command     ExeFolder      Time[s]    Avail. CPU    Max. CPU    Mean CPU    Avail. MEM[GB]    Max. MEM[GB]    Mean MEM[GB]
----------   -----------  ---------  ------------  ----------  ----------  ----------------  --------------  --------------
Tapioca      TapiocaOutput   508433           400       400        384.42             13.13            3.44            1.78
NO_AllOri2Im TPROutput          607           400       182.8       38.64             13.13            3.11            2.86
RedTieP      TPROutput           91           400        87.1       14.76             13.13            3.05            2.8
Tapas        NoTPROutput      65199           400       266.4      100.28             13.13            9.18            8.82
Tapas        TPROutput         3154           400       321.9      111.54             13.13            4.73            3.46
GCPBascule   NoTPROutput          7           400       101.1       72.01             13.13            0.7             0.68
GCPBascule   TPROutput            6           400       104         88.84             13.13            1.61            1.59
```
Tapas when reduction is done is 30x faster and uses much less RAM. We can also see the actual reduction that has been done in tie-points set in the second workflow:
```
micmac-homol-compare -o TapiocaOutput/Homol -c NoTPROutput/Homol,TPROutput/Homol
```
The output will look like:
```
###########
Ratio Homol
###########
#Name                      Homol dec
-----------------------  -----------
NoTPROutput                   1.0000
TPROutput                     0.0669
```
Only 6.7% of the tie-points were used! But did the images correctly oriented with only 6.7% of the tie-points?.
Well, let's look at the Tapas residuals.
The tool `micmac-tapas-log-anal` opens the Tapas log files, counts the number of iterations and for the last one it shows the residuals:
```
micmac-tapas-log-anal -f NoTPROutput,TPROutput
```
will report something like:
```
##########################
Tapas last residuals/worts
##########################
#Name                NumIter       Res      Wor
-----------------  ---------  --------  -------
NoTPROutput              160  0.642464  1.83815
TPROutput                132  0.727202  1.02918
```
The residuals are a bit higher if tie-points reduction is applied. We could expect this. Note that less iterations were required with less tie-points. Next, we check what really happened with the GCPs points. The tool `micmac-gcpbascule-log-anal` reads the GCPBascule logs and computes the errors of the orientation in the GCPs.
```
micmac-gcpbascule-log-anal -x gcp_List3D.xml -f NoTPROutput,TPROutput
```
will report something like:
```
###########################
GCPBascule Dists statistics
###########################
KOs
#Name
-----------------  -----------------------
NoTPROutput        -
TPR_12_12_0.00_N   -

GCPs
#Name                 Min     Max    Mean     Std    Median
-----------------  ------  ------  ------  ------  --------
NoTPROutput        0.0419  0.1191  0.0799  0.0278    0.0681
TPROutput          0.0583  0.1182  0.096   0.0242    0.1101

CPs
#Name                 Min     Max    Mean     Std    Median
-----------------  ------  ------  ------  ------  --------
NoTPROutput        0.0231  0.2374  0.0888  0.051     0.0734
TPROutput          0.0287  0.2662  0.1048  0.0616    0.0785
```

The errors (in meters) of using a reduced set of tie-points increase less than 2 centimeters.

The previous example is the sort of case in which pymicmac will make your life much easier, i.e. when you have to rerun parts of the workflow with different parameters (or commands) and then you want to compare between the different workflows.



### What else can do pymicmac for you? Distributed computing

In the previous example we ran two workflows and we compared them. We saw that by using pymicmac the whole process is a bit less tedious.
In addition to the benefits highlighted in the previous example, there are tools in pymicmac that are crucial when processing large image sets.

We saw that we can use tie-points reduction to decrease the memory usage and to speed-up Tapas (the bundle block adjustment). But what happens with Tapioca and with the dense image matching (Malt, Tawny and Nuage2Ply)? We saw in one of the tables before that while Tapas without tie-points reduction took around 60,000 seconds, Tapioca took around 500,000 seconds. That is almost a week! And the dataset was not even that large, just a few hundred images. Tie-points reduction is useful to run Tapas with large sets faster and with less memory. However, with large image sets issues also arise in Tapioca and later also in the dense image matching. In these cases, a feasible choice is to use distributed computing facilities such as clouds or clusters. pymicmac has tools to run a distributed version of Tapioca and a distributed version of the dense image matching pipeline (Malt, Tawny, Nuage2Ply) in clusters (with SGE queuing system) and in a bunch of ssh-reachable machines. In order to port these tasks to distributed computing systems, a small modification in the algorithms that perform these tasks is required. The key idea is to divide the processing in chunks that can be processed independently (and in different machines), and combine the results in the end.

Next we show how to use these tools in a SGE cluster.

#### Distributed Tapioca

In this example, we have transfered the dataset used in the previous examples to a SGE cluster. We have stored the data in the following location `/var/scratch/orubi/data/medium/4ms_60m_1000`. This is a shared location, so all the nodes of the cluster can access it. The folder contains:
```
coord_List2D.xml gcp_List3D.xml GrapheHom.xml Ori-IniCal images.list
IMG_0990.JPG ...
```

In order to run the distributed Tapioca, first we need to divide the processing in chunks. The tool `micmac-disttapioca-create-config` can be used to define the different chunks and the processing that needs to be done for each chunk.
```
micmac-disttapioca-create-config -i GrapheHom.xml -o DistributedTapioca.xml -f ChunksData -n 50
```

The previous command divides the image pairs defined in `GrapheHom.xml` in chunks of 50 image pairs (if you do not have a file like `GrapheHom.xml` you can use `micmac-disttapioca-create-pairs` to create one with all possible image pairs). Which image pairs are used in each chunk is stored in files in the `ChunksData` folder. For each chunk it defines the commands to execute and all the commands are stored in `DistributedTapioca.xml`.

Now we use the tool `coeman-par-sge` in pycoeman to run the list of commands in the cluster and in parallel:
```
cd /var/scratch/orubi/data/medium/4ms_60m_1000
coeman-par-sge -d . -c DistributedTapioca.xml -s /home/orubi/sw/export_paths.sh -r /local/orubi/DistributedTapioca -o DistributedTapiocaAllOutputs
```
`/var/scratch/orubi/data/medium/4ms_60m_1000` is our shared data location, so all paths in the XML files are relative to this location. The various commands to run in the cluster are specified in `DistributedTapioca.xml`. The file `/home/orubi/sw/export_paths.sh` is a file that sets the environment in the nodes. Once a job is going to be executed in a certain node, the node needs to have the software available (MicMac, pymicmac and pycoeman). The execution of the jobs in the nodes will be done in the location specified, i.e. `/local/orubi/DistributedTapioca`. For each job the required data for the chunk will copied from the shared data location (`/var/scratch/orubi/data/medium/4ms_60m_1000`) to the local disk (`/local/orubi/DistributedTapioca`) so each job has a local copy of the required data for faster access. The chunk will be processed locally in the node and the output data (the partial `Homol` folders) will be copied back to the shared location.

Once all the jobs have finished (check qsub) we need to combined the partial `Homol` folders:
```
micmac-disttapioca-combine -i DistributedTapiocaAllOutputs -o Homol
```
Now we have a `Homol` folder that was created in a distributed manner. We can plot the combined CPU/MEM usage for all the nodes of the cluster:
```
coeman-mon-plot-cpu-mem -i DistributedTapiocaAllOutputs -r 20
```

#### Distributed dense image matching

In the previous subsection we ran Tapioca in a distributed system. After Tapioca we need to run Tapas (and maybe GCPBascule or other processes) in order to obtain the images orientation. We saw before that we can use tie-points reduction to speed-up Tapas. After we have th eimage orientation, we are ready to generate the dense colored point cloud. In MicMac this can be done with the dense image matching pipeline that consists of Malt, Tawny and Nuage2Ply. These processes are very time consuming and will generate a large amount of intermediate data that for large datasets will fill your disk storage easily.

We have developed a distributed tool to run the dense image matching pipeline (Malt, Tawny and Nuage2Ply). Right now the solution only works for aerial images oriented in a cartographic reference system. In this example we show how to run it in a SGE cluster. We assume that we have the data (images) again in `/var/scratch/orubi/data/medium/4ms_60m_1000` and in this folder we also have the image orientation in the folder `Ori-Final`.

In order to run the distributed dense image matching, first we need to divide the processing in chunks. The tool `micmac-distmatching-create-config` can be used to define the different chunks and the processing that needs to be done for each chunk:
```
micmac-distmatching-create-config -i TPROutput/Ori-GCPBOut -e JPG -o DistributedMatching.xml -f DistributedMatchingConfigFolder -n 60,60
```
The previous command divides the area in 3600 tiles. Which images are processed in each tile is defined in the `DistributedMatchingConfigFolder` folder. For each tile it defines the commands to execute and all the commands are stored in `DistributedMatching.xml`.

Now we use the tool `coeman-par-sge` in pycoeman to run the list of commands in the cluster and in parallel:
```
cd /var/scratch/orubi/data/medium/4ms_60m_1000
coeman-par-sge -d . -c DistributedMatching.xml -s /home/orubi/sw/export_paths.sh -r /local/orubi/DistributedMatching -o DistributedMatchingAllOutputs
```

`/var/scratch/orubi/data/medium/4ms_60m_1000` is our shared data location, so all paths in the XML files are relative to this location. The various commands to run in the cluster are specified in `DistributedMatching.xml`. The file `/home/orubi/sw/export_paths.sh` is a file that sets the environment in the nodes. Once a job is going to be executed in a certain node, the node needs to have the software available (MicMac, pymicmac and pycoeman). The execution of the jobs in the nodes will be done in the location specified, i.e. `/local/orubi/DistributedMatching`. For each job the required data for the tile will copied from the shared data location (`/var/scratch/orubi/data/medium/4ms_60m_1000`) to the local disk (`/local/orubi/DistributedMatching`) so each job has a local copy of the required data for faster access. The tile will be processed locally in the node and the output data (the pointlcoud of the tile) will be copied back to the shared location.

After the exectuion is finished, we have a `DistributedMatchingAllOutputs` folder that contains subfolders and each subfolder contains a ply file. We can plot the combined CPU/MEM usage for all the nodes of the cluster:
```
coeman-mon-plot-cpu-mem -i DistributedMatchingAllOutputs -r 20
```
