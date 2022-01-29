Structure From Motion Pipeline
------------------------------

[![Build Status](https://travis-ci.org/NLeSC/structure-from-motion.svg?branch=develop)](https://travis-ci.org/NLeSC/structure-from-motion)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/13ba2c747cde4bc4ba5809873aa40e7d)](https://www.codacy.com/app/NLeSC/structure-from-motion?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=NLeSC/structure-from-motion&amp;utm_campaign=Badge_Grade)
[![Coverage Status](https://coveralls.io/repos/github/NLeSC/structure-from-motion/badge.svg?branch=)](https://coveralls.io/github/NLeSC/structure-from-motion?branch=)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.45937.svg)](http://dx.doi.org/10.5281/zenodo.45937)

Please cite the tool with its DOI if you are using it in your scientific publication.



This repo contains a complete _Structure from Motion_ pipeline. Structure from Motion is a technique to construct a 3-D point cloud from a set of images (or a video) of an object. The software in this repository relies heavily on a number of third party libaries, notably Bundler, CMVS, PMVS, and SIFT.


* Go [here](docs/install-ubuntu-14.10.md) for the installation instructions;
* A conceptual overview of the pipeline is documented [here](docs/structure_from_motion.md);
* The current pipeline has many options that can be configured. [This document](/docs/tuning_guide.md) describes which option does what and how it affects the characteristics of the resulting point cloud;
* [This document](docs/related_work.md) lists a couple of key people, their websites, and tools;
* [Here](docs/ideas.md) we describe some ideas we never found time to look into;
* You can run the pipeline with [docker](https://www.docker.com/) using [this docker image](https://hub.docker.com/r/nlesc/structure-from-motion/). Find the instructions [here](docs/docker.md).




Example
--------

![example-output](docs/images/example-output.png "Example Output")

This software includes a small example, in this case [a rock on the parking lot outside of our building](https://www.google.com/maps/place/52%C2%B021'24.6%22N+4%C2%B057'15.1%22E/@52.356789,4.9542065,49m/data=!3m1!1e3!4m2!3m1!1s0x0:0x0). See [here](docs/example.md) for some info on how to test the pipeline on the example.




Copyrights & Disclaimers
------------------------

The software is copyrighted by the Netherlands eScience Center and 
releases under the GNU general public license (GPL), Version 2.0.

See <http://www.esciencecenter.nl> for more information on the 
Netherlands eScience Center.



See the "LICENSE" and "NOTICE" files for more information. 

Contributions are welcome! Please create a pull request and adhere to the following list:

1. The code follows the standard style, check by fixing all errors
2. Use the GitHub Flow branching model
3. For other development and coding style conventions, see the NLeSC Style [Guide](https://guide.esciencecenter.nl/index.html)
4. Don't include extra dependencies without a good reason. Only use licenses compattible with the license of this project 
5. Please document your code, and provide unit tests
    # extract all frames from the video file:
    ./frames-extractor.py 00017.MTS ~/tmp/frames/in
    
    # add the same camera exif data (focal length, camera make, camera model) to each frame
    ./add-exif-data.py ~/tmp/frames/in 2.8mm Panasonic HC-X900 1920 1080
    
    # make a new directory with links to a subset of all frames, for example, each 50th frame:
    ./frame-subsetter.py ~/tmp/frames/in ~/tmp/frames/out-49 49
test dir

In some cases, it can be advantageous to preprocess the images. The repository includes scripts to:
* [crop](cropper)  the edges;
* [calculate](masker) a mask;
* [resize](resizer) image to make them more managable.

Both these rely on imagemagick for the heavy lifting. Imagemgick can be installed with:

```
sudo apt-get install imagemagick
```


Usage example:
```
./resize-images.sh --in ./testimage --out resized
```


'resize-images.sh' resizes all ``*.jpg|*.JPG`` (but not ``*.jpeg|*.JPEG`` ) in the ``--in`` directory , such that the smallest dimension will be 500 pixels long after conversion. The output is written in the ``--out`` directroy, which is created if it doesn't already exist.

'resize-images.sh' needs imagemagick's ``convert``. Install imagemagick using:

```
sudo apt-get install imagemagick
```


This script tries to guess the characteristics of the object of interest by assuming it is located in the middle of a photo. It then uses selection growing to calculate a mask for each image, cutting out the background. We found that this script works well for some sets of images, but not others.  


Example usage:
```
./crop-image-sides.sh --in testimage --out cropped --cropSides 5 --cropTop 3
```
crops one fifth from the sides of the images in the 'testimage' directory, and one third from the top (nothing from the bottom).


'crop-image-sides.sh' crops all ``*.jpg|*.JPG`` (but not ``*.jpeg|*.JPEG``) in the ``--in`` directory. The amount of cropping is controlled by the ``--cropSides`` and ``--cropTop`` arguments. 

'crop-image-sides.sh' needs imagemagick's ``convert``. Install imagemagick using

```
sudo apt-get install imagemagick
```

This directory contains the pilot job framework file used for commiting jobs to a cluster. 
The example folder (examples/rock) contains an example input for the pipeline, in this case a rock.

After [installing](install-ubuntu-14.10.md), the pipeline can be started by cd'ing into the data directory, and starting the 'run-sfm.py' script from there:

```
$ cd ${HOME}/structure-from-motion/examples/rock
$ python ../../run-sfm.py
```

Alternatively, the docker image can be used, see [here](docker.md).

```
$ cd examples/rock
$ sudo docker run -u $UID -v $PWD:/data nlesc/structure-from-motion
```

Viewing the resulting sparse pointcloud (in bundle/bundle.out) and dense pointcloud (in pmvs/models/optio-0000.ply) with for example [meshlab](http://meshlab.sourceforge.net/).

To test if the point cloud was correctly generated, you can use a test script which prints the number of points in the generated cloud:

```
$ cd ${HOME}/structure-from-motion
$ test/number_of_points.py examples/rock

# Using 'bundle.out' file from here: examples/rock/bundle/bundle.out
# Using 'option-0000.ply' from here: examples/rock/pmvs/models/option-0000.ply
# The results are:
#    nPointsSparse = 9753
#    nPointsDense = 2497111

9753
2497111
```
Using Docker to run the Stucture From Motion Pipeline
=====================================================

To facilitate running the pipeline with as little effort as possible we have created a docker image.

Docker is a system for fast deployment of applications using virtual machines. See here for more information: https://www.docker.com/



Quick HOWTO:

1. Install Docker: ```sudo apt-get install docker.io```,
1. retrieve the docker image by running ```docker pull nlesc/structure-from-motion```,
1. go to the examples directory 'examples/rock',
1. start the image using ```sudo docker run -u $UID -v "$PWD:/data" nlesc/structure-from-motion```. To process your own set of images, run this command inside the directory containing your image files.

By default the docker image will run the entire structure-from-motion pipeline on all pictures in the current working directory. If instead you would like a terminal session to play around with the image try this command:
````
sudo docker run -u $UID -v "$PWD:/data" -i -t nlesc/structure-from-motion /bin/bash
````
The main script to run the pipeline is called run-sfm.py'.
The image can also be built from source. To do this yourself, you need to checkout the submodules:
````
git submodule update --init --recursive
````
Then, build the image using the Dockerfile in the repository root directory:
````
sudo docker build -t sfm_image .
sudo docker run -u $UID -v $PWD:/data sfm_image
````
Structure From Motion
=====================

Structure from motion is a technique where a collection of images of a single object is transformed into a pointcloud.

See the Wikipedia page on Structure from Motion: http://en.m.wikipedia.org/wiki/Structure_from_motion

Basic Workflow
--------------

The process consists of 6 basic steps shown in the workflow below:

![pipeline](images/sfm.png "SFM Pipeline")

- focal point extraction -- extract the focal point and sensor size from the exif information in each image.
- keypoint detection -- detects "point of interest" in each image.
- keypoint matching -- compare keypoints of each image pair to see if and how they overlap.
- bundle adjustment -- determine the camera positions in each image, using multiple overlapping images as input. This also produces an initial sparse pointcloud.
- undistort images -- fix any distortion in the images caused by the camera.
- reconstruction of 3D structure -- combine the images into a dense pointcloud.

There are many different implementations for each step in this workflow. In our workflow we use the following combination of tools:

- [SIFT](http://www.cs.ubc.ca/~lowe/keypoints/) for keypoint detection. Note that SIFT is patented and can only be used for research purposes. 
- [Bundler](http://www.cs.cornell.edu/~snavely/bundler/) for keypoint matching, bundle adjustment and undistort images. 
- [CMVS/PMVS2](http://www.di.ens.fr/cmvs/) to reconstruct the 3D structure using multi-view stereo.

Note that we don't use the original versions provides in the links above, but instead use more up-to-date versions from github. The details can be found in the 
[install guide](./install-ubuntu-14.10.md).

This document collects some of the ideas that we never had time to look into.

Here's the list:
* [Brute force calculation of point clouds](#brute-force-calculation-of-point-clouds)
* [Quick feedback system](#quick-feedback-system)
* [Camera parameters sensitivity analysis](#camera-parameters-sensitivity-analysis)
* Improve visual quality of objects
* [Alternative keypoint detectors](#alternative-keypoint-detectors)
* [Improve accuracy of key matching by adding easily identifiable objects](#improve-accuracy-of-key-matching-by-adding-easily-identifiable-objects)

<label name="brute-force-calculation-of-point-clouds" />

### Brute force calculation of point clouds
* **context:** In our experience it is difficult to know which optimal settings to use when constructing a point cloud. There are many knobs to turn, and it's often not clear how the settings interact in terms of performance, memory requirements, quality of the result point cloud, etc.
* **proposed solution:** start construction of the point cloud using different settings, and then either combine the results, or select a good one (automatically or by asking the user for visual inspection).


<label name="quick-feedback-system" />

### Quick feedback system
* **context:** It turns out that it is quite difficult to take 'good' pictures during data acquisition. In the Via Appia data set, we generally see at least a few images not being used for the point cloud. Furthermore, some photos generate many keypoints, while others have few, and additionally, some keypoints are really informative while others aren't. Photos can be good or bad due to various reasons, such as:

    * angle between adjacent photos is too small (in particular for 'photos' derived from video frames)
    * photos that are blurry
    * photos with wrong aperture (software assumes pinhole camera)
    * photos with much background
    * photos with low contrast (dark areas such as shadows, often of the photographer; light areas such as walls and other man-made structures.

    Add to this the many settings of modern cameras, and it becomes a multidimensional nonlinear optimization problem.
 
* **proposed solution:** All in all, we think the most robust way of dealing with these factors is to come up with a system that is capable of providing quick feedback to the user. We already did a lot of work on increasing the performance of the pipeline, but this is still offline. It would be good to have quick feedback on how much the pointcloud improved as a result of the photo you _just_ took. Such a setup would probably involve wireless cameras that upload their photos to a cluster/cloud, with almost immediate feedback on the number, location, coverage of keypoints; further diagnostics on the resulting sparse/dense point clouds may be provided (albeit with a small delay, perhaps in the order of minutes). This way, archeologists can quickly get a feel for what makes a good photo for their purposes, given the prevalent lighting conditions, camera settings, photograph positions, etc., ultimately resulting in higher-quality datasets.

<label name="camera-parameters-sensitivity-analysis" />

### Camera parameters sensitivity analysis
* **context:** To get a better feel for the optimal camera model and camera settings, we could do a sensitivity analysis of different cameras, and vary the settings used on each camera
* **proposed solution:** Vary:
   * camera model
   * flash settings
   * aperture
   * ISO

<!--- I tried adding a label here, but I couldn't for the life of me figure out why this particular link does not work while the others do -->

### Improve visual quality of point clouds/objects

* **context:** The current pipeline spits out a point cloud, with colored points. The visual representation can still be improved by calculating meshes (surfaces between points) and then by calculating textures on top of the outer surface of objects. This is good for making visually attractive representations of the objects, and calculating meshes has the added bonus of being able to calculate certains metrics (e.g. volume, surface area) that may be of interest to archeologists.
* **proposed solution:** there are a couple of tools available which calculate meshes. We experimented with using meshlab. This works OK, but scripting the mesh calculating was a bit ugly (though not impossible). Perhaps other tools can be used as well, for instance Blender&Python can do mesh and texture calculations.

<label name="alternative-keypoint-detectors" />

### Alternative keypoint detectors
* **context:** We currently use SIFT to do the keypoint detection. This works OK in principle, but has the potential drawback of being patented. 
* **proposed solution:** Other keypoint detectors are available, some of which are supposedly quicker (although that is not really where most time is spent, so maybe it's not worth optimizing). Avoiding license issues may be a reason to switch from SIFT to something else though. Also, it's worth investigating whether the results from different keypoint identifiers can be concatenated for a better result.

<label name="improve-accuracy-of-key-matching-by-adding-easily-identifiable-objects" />

### Improve accuracy of key matching by adding easily identifiable objects
* **context:** Key matching is sometimes difficult, in particular when the object has symmetry or repeating shapes (e.g. standard windows, pillars, tiles, etc).
* **proposed solution:** Adding small objects to the scene before the photographs are taken can help correctly stitch together the photographs. Ideally, the objects are rigid, high contrast, and uniquely identifiable from any angle and distance. Perhaps [QR codes](http://en.wikipedia.org/wiki/QR_code) could be used.
























































Tuning guide for the structure from motion pipeline
===================================================

Many of the tools in the structure from motion pipeline require tuning
to improve the quality of the output and/or the performance. For many
tools, the best configuration to use may also depend on the image
resolution or number of images that are used. In this document we 
describe what setting we have tried so far.

SIFT
----

Sift generates the keypoints for each image. Each keypoint describes a
_distinctive feature_ in the image in a scale an rotation independent 
way. By matching the keypoints in each image with all other images, the
SfM pipeline can determine which images overlap (partly). 

More information of sift can be found 
[here](http://en.wikipedia.org/wiki/Scale-invariant_feature_transform).

The version of sift we use can be found in the ``bundler_sfm/src/Sift.cpp`` 
file. In this sift implementation, the following settings are important:

- ``SIFT::DoubleImSize`` this setting determines if the the image 
  should be doubled in size before the sift algorithm is run. Sift internally
  downscales the image repeatedly to detect features of different sizes. It 
  always downscales once before detecting the first features. Therefore, very 
  small features cannot be detected unless the image is doubled in size first.
  It is typically good to set this to ``true`` for low resolution images 
  (e.g. 1024x768) and ``false`` for high resolution images (as produced by 
  modern cameras). In our pipeline the default is ``false``.
  
- ``SIFT::PeakThreshInit`` this setting determines the minimum contrast required 
  for a point to be considered as a keypoint. Dark areas in images typically result 
  in _noisy_ keypoints which are easily confused with others. In our pipeline the
  default is ``0.08``. Using lower values will include more keypoints from darker 
  areas.
  
KeyPoint matching
-----------------

After sift, the keypoints generated for the images are compared using a keypoint matcher. 
The matchers we use, ``KeyMatchFull`` or ``KeyMatchPart`` are part of the bundler tool set. 
These matchers compare the keypoints using a _approximate nearest neighbor KD tree_. More
information on the implementation of these trees can be found [here](https://www.cs.umd.edu/~mount/ANN/)

The matcher we of use can be found in the ``bundler_sfm/src/KeyMatchPart.cpp`` 
file. In this matcher implementation, the following settings are important:

- ``ratio`` is the fifth (optional) parameter to KeyMatchPart. During matching, each keypoint in one image is compared to all keypoints in another image by computing the euclidean distance between the feature vectors of the two keypoints. The ratio of the distance of the two best matches (the best and the runner up) are then computed. If this ratio is close to 1, the match is considered to be bad, since there are multiple 'potential' matches for the keypoint. If the ratio is closer to 0, the match is considered good, since the difference (and thus the distance) between the best match and the runner up is large. The ratio parameter determines the threshold above which matches are discarded. The default in our pipeline is `0.6`. Higher values will make the matching less strict and thus produce more matches of lower quality. 


Bundler
-------

After matching, bundler takes the keypoints and attempts to reconstruct the camera positions in 3D space for each images. In addition, a sparse point cloud is created that contains a relatively small number of object points in 3D space. 

Bundler reads it configuration from a text file `options.txt` which is generated by the `run-sfm.py` script we us in our pipeline. In this configuration file, the following settings are important:

- ``--use_ceres`` is used to switch between the ceres solver and the internal solver of bundler. Using ceres significantly improves the performance of bundler, as it is capable of using all cores in a machine. 

- ``--construct_max_connectivity`` is used to instruct bundler to add images in the order of how _connected_ they are to other images. That is, images containing features than can be matched to many other images are added first. The alternative is to add images based on the number of matches, which tends to add images based on how similar the are to the current set. 

- ``--projection_estimation_threshold 1.1`` is the RANSAC threshold used when doing pose estimation to add in a new image. Lower values will result in a stricter selection of which estimates are valid. The default value used in bundler is 4, which seemed to result in too noisy estimates of the camera positions. A lower value resulted higher quality result. 

CMVS/PMVS2
----------

After estimating the camera positions, PMVS2 is used to estimate the 3D positions of object points. Before 
PMVS2 is run, CMVS is used to read the output of bundler and create a configuration file for PMVS2 in `pmvs/option-0000`. In this configuration file specifies which images should be used for the reconstruction, and at what resolution the input images should be used. In this configuration file, the following settings are important:

- ``timages`` the images actually used in the reconstruction of the 3D object. Usually only a subset of the input images is used. 

- ``level`` this setting determines how much the input images are down sampled before 3D reconstruction. Level 0 means full resolution, 1 uses half the resolution, etc. We use 0 for the best result. Reducing this to 1 will significantly reduce the computation time and the number of points in the result.

- ``threshold`` this setting determines which patch reconstructions are accepted. Lower values will accept more patches but will produce a in noisier result. Higher values will accept less patches and will produce a result with less errors, but more missing points. We use the default of `0.7`.

- ``maxAngle`` determines the minimal angle between cameras before they are considered for 3D reconstruction. If the baseline between the cameras is too small, the 3D reconstruction tend to have higher errors. We use the default angle of 10 degrees.  

- ``CPU`` determined the number of cores used by PMVS2. Our `run-sfm.py` script detects the number of 
cores automatically and generated the correct configuration.


 





Install guide
=============
This install guide explains how to install the structure from motion pipeline in Ubuntu 14.10. 

* [Setting up a virtual machine (optional)](#set-up-a-virtual-machine)
* [Installing packages from Ubuntu repositories](#install-packages-from-ubuntu-repositories)
* [Downloading and installing other tools](#download-and-install-other-tools)
* [Running an example](#run-an-example)

<a name="set-up-a-virtual-machine"></a>
## Setting up a virtual machine (optional)


Creating a virtual machine is an optional step, for example, if you are using Windows. You can skip this step and go directly to [Installing packages from Ubuntu repositories](#install-packages-from-ubuntu-repositories) instead. 

Download and install the latest [VirtualBox](https://www.virtualbox.org/wiki/Downloads) if you have't already. 

Download the latest [Ubuntu](https://www.ubuntu.com/download/desktop) iso. 

Create an image in VirtualBox and install Ubuntu. 

I configured virtualbox to use:

  * 5000 MB memory
  * virtual harddisk 
      * type: VDI
      * dynamically allocated storage
      * 64 GB diskspace
      
After the VM has been created, you can set the following properties:

  * 2 cores
  * 128 MB video memory
  * 3D acceleration enabled
  
These options are mainly determined by the limitations of the machine you run on (this is about as much as my laptop can handle). Generally, using more cores and more memory is a good idea.      

Start the virtual machine, install Ubuntu from the iso we just downloaded. Tick the box about downloading any available updates when installing.

We also installed the following optional packages:

    sudo apt-get install virtualbox-guest-utils 
    sudo apt-get install virtualbox-guest-x11
    
Doing so allows you to share the clipboard between the host and the guest.


<a name="install-packages-from-ubuntu-repositories"></a>
## Installing packages from the Ubuntu repositories

Once you have Ubuntu up and running we need to install the necessary tools and libraries. Open a terminal  and install the following packages:

```
# git
# You will need git to clone the lastest versions of the structure from motion
# software from github:
sudo apt-get install git 

# cmake
# You need cmake to generate the Makefiles needed to build the
# structure from motion software from github:
sudo apt-get install cmake

# gfortran
# You need a Fortran compiler to compile (parts of) the structure from motion software:
sudo apt-get install gfortran

# Glog
# a logging library from google (https://github.com/google/glog):
sudo apt-get install libgoogle-glog-dev

# Atlas
# The "Automatically Tuned Linear Algebra Software" provides C and Fortran77 
# interfaces to a portably efficient BLAS implementation, as well as a few routines 
# from LAPACK (http://math-atlas.sourceforge.net/):
sudo apt-get install libatlas-base-dev

# Eigen3
# C++ template library for linear algebra: matrices, vectors, numerical solvers,
# and related algorithms (http://eigen.tuxfamily.org).
sudo apt-get install libeigen3-dev
    
# SuiteSparse
# suite of sparse matrix algorithms (http://faculty.cse.tamu.edu/davis/suitesparse.html):
sudo apt-get install libsuitesparse-dev

# zlib
# a library implementing the deflate compression method found in gzip and PKZIP:
sudo apt-get install zlib1g-dev

# libjpeg
# a library implementing the loading of jpeg images:
sudo apt-get install libjpeg-dev

# libboost
# library with 'all the features you wanted in C++ but weren't there'
sudo apt-get install libboost-dev

# Python imaging library may not be installed by default on the 
# lighter flavors of Ubuntu (e.g. Lubuntu 14.10 or Ubuntu 14.10 server) 
sudo apt-get install python-pil

# for viewing the point clouds afterwards
sudo apt-get install meshlab

```

<a name="download-and-install-other-tools"></a>

## Downloading and installing other tools


### Cloning NLeSC's structure-from-motion repository

Our repository includes two other repositories as submodules: 

* [bundler_sfm](http://www.cs.cornell.edu/~snavely/bundler/)
* [cmvs](http://www.di.ens.fr/cmvs/)/[pmvs](http://www.di.ens.fr/pmvs/)

To make sure you get the contents of the submodules when checking out the structure-from-motion repository, use the ``--recursive`` option to ``git clone``:


```
cd ${HOME}
git clone --recursive https://github.com/NLeSC/structure-from-motion.git
  
```




### Installing Ceres


The [Ceres Solver](http://ceres-solver.org) is _"an open source C++ library for modeling and solving large, 
complicated optimization problems. It is a feature rich, mature and performant library which has been used
in production at Google since 2010."_ This solver is needed by the _bundle adjustment_ step of the structure from motion (SfM) pipeline. We already installed the Ceres dependencies (originally described 
[here](http://ceres-solver.org/building.html)) in the previous section, so now we can proceed to download and install Ceres:

```
cd ${HOME}/structure-from-motion
wget http://ceres-solver.org/ceres-solver-1.10.0.tar.gz
tar zxf ceres-solver-1.10.0.tar.gz
mkdir ceres-bin
cd ceres-bin
cmake ../ceres-solver-1.10.0
make -j3
make test
sudo make install
  
```


### Compiling Bundler


[Bundler](http://www.cs.cornell.edu/~snavely/bundler/) is a structure-from-motion (SfM) system for unordered
image collections. Bundler takes a set of images, image features, and image matches as input, and produces a 
3D reconstruction of camera and (sparse) scene geometry as output.

Next, compile bundler_sfm: 

```
cd ${HOME}/structure-from-motion/bundler_sfm
make
  
```


### Compiling CMVS/PMVS2

[PMVS2](http://www.di.ens.fr/pmvs/) is multi-view stereo software that takes a set of images and camera 
parameters (generated by bundler), and then reconstructs 3D structure of an object or a scene visible in the images. The software outputs a _dense point cloud_, that is, a set of oriented points where both the 3D coordinate and the surface normal are estimated for each point. 

[CMVS](http://www.di.ens.fr/cmvs/) is software for _clustering views for multi-view stereo_. It is basically a 
pre-processor for PMVS2 that takes the output of bundler and generates one or more (optimized) configuration files for PMVS2. CMVS is normally used to split the PMVS2 processing in multiple independent parts, for example when creating a 3D reconstruction on the basis of thousands of images (which would be too much for PMVS2 to handle all at once). However, even when the number of images used is small, there is an advantage in using CMVS as it also removes unused images from the data set, and provides the order in which PMVS2 should process the images. This significantly reduces the processing time needed by PMVS2. The version we included in the structure-from-motion repository is a fork of [pmoulon](https://github.com/pmoulon/CMVS-PMVS). It contains both CMVS and PVMS2, adds a cmake configuration, and contains several bug and performance fixes. 

Compile CMVS/PMVS like this:

```
cd ${HOME}/structure-from-motion/cmvs-pmvs/program
mkdir build
cd build
cmake ..
make
  
```






<a name="run-an-example"></a>
## Running an example

The pipeline can be started by ``cd``'ing into a data directory, and starting the 'run-sfm.py' script from there:

```
cd ${HOME}/structure-from-motion/examples/rock
python ../../run-sfm.py
  
```
On my laptop the example finished in 1 hour. After finishing, you can open a .ply file from ```${HOME}/structure-from-motion/examples/rock/bundle``` with MeshLab to see the outcome. 
This documents gives a small overview of the various developments on structure-from-motion that we found and/or have experience with

List of SFM Tools
=================

Integrated tools for SFM
------------------------

https://github.com/dddExperiments/SFMToolkit  
http://www.visual-experiments.com/demos/sfmtoolkit/  
Theia: http://cs.ucsb.edu/~cmsweeney/theia/sfm.html  

Keypoint Detection
------------------------

Image libraries that contain SIFT/SURF/BRISK

http://opencv.org/about.html  
(An example:  http://stackoverflow.com/questions/5461148/sift-implementation-with-opencv-2-2)

http://www.vlfeat.org/

### Sift

Good explanation of what sift does:

http://www.aishack.in/2010/05/sift-scale-invariant-feature-transform/  

Original sift:

http://www.cs.ubc.ca/~lowe/keypoints/  

Some alternative implementations of sift:

* http://www.robots.ox.ac.uk/~vedaldi/code/siftpp.html
* http://robwhess.github.io/opensift/
* http://www.cs.unc.edu/~ccwu/siftgpu/

### Surf

http://www.vision.ee.ethz.ch/~surf/

### Brisk

https://github.com/rghunter/BRISK

Bundle adjustment
-----------------

Bundler Tool  

http://www.cs.cornell.edu/~snavely/  
https://github.com/snavely  

http://grail.cs.washington.edu/projects/mcba/  


Theia 

https://github.com/sweeneychris/TheiaSfM

Theia is an alternative to bundler (and the processing pipeline proceeding bundler). 
It consists of a library containing all the elements needed to do bundle adjustment
(keypoint detection, keypoint matching, etc.) and contains several example applications 
implementing the entire pipeline. Theia contains several state-of-the-art algorithms, 
such as a cascade hashing based keypoint matching, and a global SfM approach that 
considers the entire view graph at the same time instead of incrementally adding 
more and more images to the reconstruction. Late 2014, Theia was still in active 
development and not completely stable, but it is likely to become an efficient 
replacement for bundler. 

Clustering
----------

http://www.di.ens.fr/pmvs/  
http://www.di.ens.fr/cmvs/  

Misc
----

These libraries are used by some of the steps:

http://www.cs.utexas.edu/users/dml/Software/graclus.html  
http://ceres-solver.org/  


Relevant Papers
===============

http://foto.hut.fi/seura/julkaisut/pjf/pjf_e/2014/PJF2014_Lehtola_et_al.pdf
