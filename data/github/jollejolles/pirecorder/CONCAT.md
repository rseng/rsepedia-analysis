[![Downloads](https://pepy.tech/badge/pirecorder)](https://pepy.tech/project/pirecorder)

# pirecorder
**A Python package for controlled and automated image and video recording with the raspberry pi**

*pirecorder* is a Python package, built on the [picamera](http://picamera.readthedocs.io/) and [OpenCV](https://opencv.org/) libraries, that provides a flexible solution for the collection of consistent image and video data with the raspberry pi. It was developed to overcome the need for a complete solution to help researchers, especially those with limited coding skills, to easily set up and configure their raspberry pi to run large numbers of controlled and automated image and video recordings using optimal settings.

A paper accompanying this package is published in the Journal of Open Source Software:

*Jolles, J.W. (2020). pirecorder: controlled and automated image and video recording with the raspberry pi. J. Open Source Softw. 5, 2584. doi: [10.21105/joss.02584](http://doi.org/10.21105/joss.02584).*

<p align="center"><img src="https://github.com/jollejolles/pirecorder/blob/master/images/pirecorder-logo-large.jpg"></p>

## Key Features
* **Controlled recording using custom, easy-to-edit configuration files**
* **Record single images and videos, timelapses, and sequences of videos**
* **Configure camera settings interactively via a live camera stream**
* **Dynamically draw the region of interest for your recordings**
* **Automatic naming of files and folders with relevant and custom labels**
* **Easy scheduling and automating recordings in the future**
* **Direct control of all modules via simple terminal commands**
* **Convert (folders of) images and videos with resize, monitor, and label options**
* **Dedicated documentation website with detailed guides and tutorials**
* **Jupyter notebook tutorial files**

## An overview
<p align="center"><a href="https://www.youtube.com/watch?v=pcVHpijd6wc" title="A quick overview of the pirecorder package" target="_blank"><img src="https://github.com/jollejolles/pirecorder/blob/master/images/pirecorder-video.jpg" width="400" height="298"></a></p>

## Modules
*pirecorder* consists of a main `PiRecorder` module to run image and video recordings, `stream` and `camconfig` modules with interactive user interfaces for help setting up, calibrating, and configuring the camera, a `schedule` module for scheduling future recordings, and a `convert` module for the easy converting of (folders of) recorded images and videos.

## Install
**Note:** pirecorder relies on picamera, which is not properly integrated in the latest Raspberry Pi OS (Bullseye). Therefore use the previous OS Buster to use pirecorder and picamera. You can download the OS [here](https://downloads.raspberrypi.org/raspios_armhf/images/raspios_armhf-2021-05-28/) and find more information from the Raspberry Pi foundation [here](https://www.raspberrypi.com/news/bullseye-camera-system/). I am frustrated about this remarkable change too.

To install the latest release, simply open a terminal window and enter:

```
pip install pirecorder
```

To install the latest development version, enter:

```
pip install git+https://github.com/jollejolles/pirecorder.git --upgrade
```

## Dependencies
*pirecorder* builds strongly on the [picamera](http://picamera.readthedocs.io/) package, uses [numpy](http://www.numpy.org/), [pyyaml](https://pyyaml.org), and [opencv](http://opencv.org) for some of its core functionality, and relies on various utility functions of the accompanying [pythutils](https://github.com/jolle/pythutils) package. Scheduling functionality is based on *CronTab* and the associated [python-crontab](https://pypi.org/project/python-crontab/) package. All dependencies are automatically installed with *pirecorder* except for:
* *OpenCV*: has to be manually installed due to its various dependencies on raspberry pi. Click [here](https://github.com/JolleJolles/pirecorder/tree/master/docs/other/install-opencv.md) for a quick install guide.
* *FFmpeg*: is only needed for the convert functionality of *pirecorder*. Click [here](https://github.com/JolleJolles/pirecorder/tree/master/docs/other/install-ffmpeg-raspberry-pi.md) and [here](https://github.com/JolleJolles/pirecorder/tree/master/docs/other/install-ffmpeg-osx.md) for guides to install on raspberry pi and OS X respectively.

## Documentation
For detailed documentation and tutorials about *pirecorder* and all its functionalities, see the dedicated website [jollejolles.github.io/pirecorder/](http://jollejolles.github.io/pirecorder/).
1. [quick guide ](https://jollejolles.github.io/pirecorder/quick-guide.html)
2. [the pirecorder package](https://jollejolles.github.io/pirecorder/pirecorder-package.html)
3. [setting up your raspberry pi](https://jollejolles.github.io/pirecorder/1-setting-up-raspberry-pi.html)
4. [installing pirecorder](https://jollejolles.github.io/pirecorder/2-installing-pirecorder.html)
5. [position and calibrate the camera](https://jollejolles.github.io/pirecorder/3-position-and-calibrate-camera.html)
6. [configure recording settings](https://jollejolles.github.io/pirecorder/4-configure-recording-settings.html)
7. [configure camera settings](https://jollejolles.github.io/pirecorder/5-configure-camera-settings.html)
8. [record and schedule recordings](https://jollejolles.github.io/pirecorder/6-recording-and-scheduling.html)
9. [converting media](https://jollejolles.github.io/pirecorder/7-convert-media.html)
10. [run from the command line](https://jollejolles.github.io/pirecorder/8-run-from-commandline.html)

## Tests
To test all functionalities of the pirecorder package, run the `tests/test.py` file ([here](https://github.com/JolleJolles/pirecorder/tree/master/tests/test.py)), or alternatively run commands manually using the documented jupyter files [here](https://github.com/JolleJolles/pirecorder/tree/master/notebooks). Note that running the tests will require user input as some of the functionalities are interactive.

## Development
*pirecorder* is developed by Dr Jolle Jolles, a research fellow at the Max Planck Institute of Animal Behavior, and at the Zukunftskolleg, Institute of Advanced Study at the University of Konstanz.

For an overview of version changes see the [CHANGELOG](https://github.com/jollejolles/pirecorder/blob/master/CHANGELOG) and for detailed changes see the [commits page](https://github.com/jollejolles/pirecorder/commits/). Please submit bugs or feature requests to the GitHub issue tracker [here](https://github.com/jollejolles/pirecorder/issues).

Contributions to this package are welcomed via the usual pull request mechanism.

## Citing
If you use pirecorder in your research, please cite the accompanying paper:

```
@misc{Jolles2020,
      title = {pirecorder: controlled and automated image and video recording with the raspberry pi},
      author = {Jolles, Jolle W.},
      year = {2020},
      volume = {5},
      number = {54},
      pages = {2584},
      doi = {https://doi.org/10.21105/joss.02584}
}
```

## License
Released under a Apache 2.0 License. See [LICENSE](https://github.com/JolleJolles/pirecorder/blob/master/LICENSE) for details.
---
title: 'pirecorder: Controlled and automated image and video recording with the raspberry pi'
tags:
  - Python
  - raspberry pi
  - camera
  - recording
  - video
  - automation
authors:
  - name: Jolle W. Jolles
    orcid: 0000-0003-0872-7098
    affiliation: "1, 2"
affiliations:
 - name: Department of Collective Behaviour, Max Planck Institute of Animal Behaviour, Konstanz, Germany
   index: 1
 - name: Zukunftskolleg, Institute of Advanced Study, University of Konstanz, Germany
   index: 2
date: 6 Jul 2020
bibliography: paper.bib
---

# Summary
A fundamental component of empirical research is the acquisition of accurate, consistent, and often significant amounts of data. Specifically, researchers often require large numbers of controlled and often parallel image and video recordings. For this the raspberry pi, a small, single-board computer that brings together open-source principles with sensor and controller interfaces, and highly customisable programming capabilities, provides a great, low cost solution. Indeed, in recent years, the raspberry pi has been increasingly taken up by the scientific community [@Fletcher2019] and used in a wide range of projects that required the collection of high quality image data, from sub-micron resolution microscopy [@Aidukas2019], and deep sea video recordings [@Phillips2019], to motion-triggered camera trapping [@Nazir2017; @Prinz2016], high-throughput behavioural assessments [@Geissmann2017; @Todd2017; @Jolles2019], long-term home cage monitoring [@Singh2019], and the automated tracking of animal groups [@Alarcon-Nieto2018; @Jolles2018; @Jolles2020].

# Statement of need
So far, researchers have often relied on writing their own recordings scripts to take still photographs and videos from the command line (using `raspistill` and `raspivid`), control the camera module with `picamera` in Python [@Jones2017], or trigger recordings by motion-detection  ([Motion](https://motion-project.github.io)). Also some specific solutions exist, such as a web-based interface to run recordings [@Singh2019] and advanced software that converts the raspberry pi in a dedicated behavioural profiling machine [@Geissmann2017]. However, there is still a need for a complete solution that helps researchers, especially those with limited coding skills, to easily set up and configure their raspberry pi and run large numbers of controlled and automated image and video recordings. Here I present `pirecorder` to overcome this need.

# Functionality
`pirecorder` is a Python package, built on the picamera [@Jones2017] and OpenCV [@Bradski2000] libraries, that provides a flexible solution for the collection of consistent image and video data. It consists of a number of interconnected modules to facilitate key aspects of media recording: 1) setting-up and configuring the camera, 2) recording images, videos, time-lapses, and standardised video sequences with automatic file-naming, 3) easy scheduling of future recordings, and 4) converting of recorded media with resize, timestamp, and monitoring options. All functionalities are designed to make it very straightforward, even for users with limited coding experience, to configure, initiate, schedule, and convert recordings. In particular, `pirecorder` offers interactive streaming functionalities to facilitate users in positioning and focusing the camera, selecting the desired white-balance and other image parameters using trackbars, and set the ideal camera shutter speed. Furthermore, `pirecorder` comes with a dedicated documentation website with detailed information and tutorials ([jollejolles.github.io/pirecorder](https://jollejolles.github.io/pirecorder/)) as well as a set of annotated [Jupyter Notebooks](https://github.com/JolleJolles/pirecorder/tree/master/notebooks) to help users integrate the raspberry pi and `pirecorder` in their work.

![Screenshots of pirecorder in action, from configuring the camera with the interactive video stream, running recordings, testing and scheduling future recordings, and converting recorded media.](Figure1.jpg)

A core functionality of `pirecorder` is that it works with configuration files. These files make it possible to store a wide range of camera and recording settings that are then automatically used for recordings without requiring further user input. Furthermore, multiple configuration files can be stored and used, such as to easily start recordings for different experimental contexts or treatments. Configuration files can be edited directly, or parameters can be set in python or using the interactive video stream functionalities. Recordings can be easily initiated remotely, such as via an SSH connection, and scheduled to automatically start and stop at specific times in the future. By its use of configuration files and the automatic naming of files, `pirecorder` also makes it possible to start controlled recordings on multiple raspberry pi's simultaneously, such as with [csshX](https://github.com/brockgr/csshx), which sends a command to multiple computers at once.  

# Use cases
`pirecorder` has already been used successfully in a number of studies, such as to facilitate the high-throughput recording of large numbers of individuals and shoals of fish [@Jolles2018; @Jolles2019; @Jolles2020] and more recently, the autonomous long-term recording of fish each day, every day, for the first four-month of their life (in prep). By facilitating and streamlining controlled and automated image and video recordings, I hope `pirecorder` will be used by scientists to help simplify and improve their collection of high quality data and thereby ultimately enhance their research.

# Acknowledgements
I would like to thank Lucas Koerner and Jan Heuschele for helpful feedback on the package and paper. This work was made possible by a Postdoctoral Fellowship from the Alexander von Humboldt Foundation, a  Postdoctoral Fellowship from the Zukunftskolleg, Institute of Advanced Study at the University of Konstanz, and a research grant from the Dr. J.L. Dobberke Foundation.

# References
---
layout: page
title: 7 Converting media
nav_order: 9
---
# Convert image and video files
{: .no_toc }

The convert module of the *pirecorder* package facilitates the converting of recorded media, both of individual images and videos as well as directories of videos, with the ability to resize, add timestamps on each frame, and monitor folders for automatic conversion. It can be set to optimally use the number of computer cores available, can be run directly from the command line, and works across different operating systems.
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta .fs-4 .fw-300 }

1. TOC
{:toc}
---

## Install dependencies

As dependencies both `FFmpeg` and `OpenCV` are needed to optimally convert the different media types. To help install FFmpeg on raspberry pi follow [this guide](other/install-ffmpeg-raspberry-pi.md), to install ffmpeg on os X follow [this guide](other/install-ffmpeg-osx.md), and to install OpenCV, follow [this guide](other/install-opencv.md).

## Convert a directory of videos
To convert a directory of videos, simply import and run the `Convert` module and provide the directory where the media files are located (`indir` parameter), and where you want to store them (`outdir` parameter). This conversion is very fast. By default it will look for files with filetype "h264", but this can be changed with the `type` parameter. For example, to convert all videos in a folder called `videos` and store them in a nested folder called `converted` run:

```
from pirecorder import Convert
Convert(indir = "videos", outdir = "videos/converted")
```

If the `outdir` does not exist it will be automatically created. By default, videos that are converted will not be overwritten (`overwrite=False`). Therefore, if some or all videos in the `indir` are already converted it will skip those automatically.

## Continuously monitor a folder for new files
With the `Convert` module it is also possible to continuously monitor a folder for new files with a set delay to wait between subsequent checks. This makes it easy to integrate with an automated media recording workflow. Simply add `sleeptime=XX` where `XX` is the time in seconds between subsequent checks.

## Convert images to video

Using the `Convert` module it is also easy to convert a directory of (timelapse) images to video. For this you need to set the `type` parameter to the image format you use, and set the `imgfps` parameter to the desired framerate of the video. For example, to create a video of 10fps:

```
Convert(indir = "media/vidimages", outdir = "media", type = ".png", imgfps = 10)
```

To then for example change the framerate to 30fps, simple run:

```
Convert(indir = "media/vidimages", outdir = "media", type = ".png", imgfps = 30, overwrite = True)
```

To convert a folder consisting of multiple image folders, you can use the `listfiles` function from my [pythutils package](https://github.com/jollejolles/pythutils), which is automatically installed with `pirecorder`, as follows:

```
imagedirs = listfiles(dir = "parentdir", type = "dir", keepdir = True)
for imagedir in imagedirs:
    Convert(indir = imagedir, outdir = "converted", type = ".png")
```

## Set number of converting pools
The `Convert` module can run multiple conversions at the same time. This works optimally when linked to the number of processing cores of the computer being used. By default it will presume a minimum of 4 cores are available. To change this, use the `pools` parameter, e.g. `pools = 6`. When running the convert functionality from python you can stop the converting by entering `ctrl+c`, and when running it in a jupyter notebook simply press the stop button in the menu bar.

## Delete originals
By default the original videos are not deleted. If you want this to be done automatically, such as when you incorporate the `Convert` module in your own automation functions, then add `delete = True`. Be careful to test your desired functionality first to not loose any unwanted data!

## Display frame number on each frame
To display the frame number on the top-left corner of each frame, simply add `withframe = True`. The conversion will now be a bit slower as it will go through each frame to draw the frame number, but due to the pooling should still work fine.

## Resize the video
By default the generated media will have the same dimensions as the originals. However, these dimensions can be decreased (as well as increased if needed), such as when wanting to reduce the file size. To do so use the `resizeval` parameter, which defaults to 1 to keep the same size. For example, to create a video with half the dimensions of the originals use `resizeval = 0.5`.

## Convert directly from the command line
When pirecorder is installed it is also possible to use the pirecorder convert functionality straight from the command line using the `convert` command, just like any other native command. You can use the same parameters as when running the `Convert` command in python, e.g.:

```
convert --indir VIDEOS --outdir CONVERTED --type ".h264" --withframe True /
--pools 4 --resizeval 0.5 --sleeptime 5 --delete False
```

---
Convert module documentation
{: .text-delta .fs-5}

```
Module to convert a directory of media files with potential to resize the
media, write the unique frame number on each frame, and continuously monitor
a folder for updated files. Multiple files can be converted simultaneously
with the pools parameter that optimally uses the computer's processing cores.

Parameters
-----------
indir : str, default = ""
    Directory containing the videos.
outdir : str, default = ""
    Directory where the converted videos should be stored. If the Directory
    does not exist yet it will be newly created.
type : str, default = ".h264"
    The filetype of the media to convert.
withframe : bool, default = False
    Type of conversion, either very fast conversion using FFmpeg or
    using OpenCV to draw the frame number on each video frame.
delete : bool, default = False
    If the original videos should be delete or not.
pools : int, default = 4
    Number of simultaneous converting processing that should be allowed.
    Works optimally when equal to the number of computer processing cores.
resizeval : float, default = 1
    Float value to which video should be resized.
imgfps : int, default = 25
    Framerate for conversion of images to video.
sleeptime : 2, default = None
    Time in seconds between subsequent checks of file folder. To not
    continuously monitor a folder set to None.
```
---
layout: page
title: 2 Installing pirecorder
nav_order: 4
---

# Installing pirecorder
{: .no_toc }

Okay, so you got your raspberry pi up and running and fully set up. Now it is time to install the pirecorder package.
{: .fs-6 .fw-300 }

The short explanation is very simple: just open a terminal window and enter `pip install pirecorder`. The more detailed explanation can be found below. Note also that pirecorder can be installed on non-raspberry pi systems, such as to use the `stream` and `convert` modules.

## Table of contents
{: .no_toc .text-delta .fs-4 .fw-300 }

1. TOC
{:toc}
---

## Creating a virtual environment

When you start working with Python, it is great practice to create isolated Python environments to work on your specific projects. The standard python environment is used by a large number of system scripts and therefore best to leave alone. I therefore strongly suggest to start by creating a virtual Python environment and install pirecorder there.  

To create virtual Python environments on your raspberry pi, first we need to install the virtual environment modules. I recommend to use `virtualenv`. To install this module and the helpful wrapper module type in:

```
sudo pip3 install virtualenv virtualenvwrapper
```

To get it to work easily on the command line, edit the file .bashrc (`nano ~/.bashrc`) and append the following lines to the bottom of the file:

```
#Virtualenvwrapper settings:
export VIRTUALENVWRAPPER_PYTHON=/usr/bin/python3
export WORKON_HOME=$HOME/.virtualenvs
export VIRTUALENVWRAPPER_VIRTUALENV=/usr/local/bin/virtualenv
source /usr/local/bin/virtualenvwrapper.sh
export VIRTUALENVWRAPPER_ENV_BIN_DIR=bin
```

Exit the file (`ctrl+x`, then `y`, and `Enter`). Now reload the file to make the changes come into effect:

```
source ~/.bashrc
```

Now we can simply create a new virtual environment with the command:

```
mkvirtualenv NAME
```

where NAME is the name you want to give to your virtual environment. Each virtual environment will have their own installed python packages, which can be checked with the command `pip freeze`.

## Installing pirecorder and dependencies

It is very easy to install pirecorder with pip. First make sure you are in the desired virtual environment if you have one (e.g. `workon MYENV`) and simply enter:

```
pip install pirecorder
```

This command will also install the majority of dependencies, which include `numpy` and my `pythutils` package [link](https://github.com/jollejolles/pythutils). One of the main dependencies, `picamera` should already be installed by default. If it is not, follow the simple steps [here](https://picamera.readthedocs.io/en/release-1.13/install.html).

`OpenCV` is a dependency that needs to be manually installed. Follow my 5min guide for Mac, ubuntu and raspberry pi [here](other/install-opencv.md). And if you plan on using the converter functionality then you will additionally need to install `FFmpeg`. Click [here](other/install-ffmpeg-raspberry-pi.md) for my guide to install it on raspberry pi and [here](other/install-ffmpeg-osx.md) for my guide to install it on OS X.

You should now be fully set up and have pirecorder working. To test it, open a terminal window and type in `python3` to enter python, and then import the pirecorder module:

```
import pirecorder
```

If you don't get any message, pirecorder is installed succesfully!

## Running the PiRecorder module for the first time

The main functionality of the pirecorder package is the `PiRecorder` module. To use it, simply create a recorder instance:

```
rec = pirecorder.PiRecorder()
```

As the PiRecorder functionality is a class instance it needs to be stored as a variable. Above we used the variable name `Rec`, but any variable is fine as long as you are consistent in using it.

The first time the PiRecorder instance is run, automatically a `pirecorder` directory will be created in the user's home directory with a default configuration file (`pirecorder.conf`) that will be used as the basis for future recordings. To also make it is easy to see what recordings you did back in time, all commands and output created with the PiRecorder instance will be stored in a log file in the setup directory with a date and time stamp.

## Testing
It is  possible to test all the functionalities of the pirecorder package. Simply run the `tests/test.py` file ([here](https://github.com/JolleJolles/pirecorder/tree/master/tests/test.py)), or alternatively run the main functionalities manually using the documented jupyter notebook files [here](https://github.com/JolleJolles/pirecorder/tree/master/notebooks). Note that running the tests will require user input as some of the functionalities are interactive.

---
PiRecorder module documentation
{: .text-delta .fs-5}

```
Sets up the rpi with a pirecorder folder with configuration and log files
and initiates a recorder instance for controlled image and video recording

Parameters
----------
configfile : str, default = "pirecorder.conf"
    The name of the configuration file to be used for recordings. If the
    file does not exist yet, automatically a new file with default
    configuration values will be created.

Returns
-------
self : class
    PiRecorder class instance that can be used to set the configuration,
    start a video stream to calibrate and configure the camera, to set the
    shutterspeed and white balance automatically, to start recordings, and
    to schedule future recordings.
```
---
layout: page
title: The pirecorder package
nav_order: 2
---
# The pirecorder package
{: .no_toc }

*pirecorder* is a Python package with a number of inter-connected modules, developed with the aim to facilitate running controlled and automatic image and video recordings using optimal settings with the raspberry pi.
{: .fs-6 .fw-300 }

So far, researchers have often relied on writing their own recordings scripts to take still photographs and videos from the command line. Although so some specific software solutions exist, what was missing is a complete solution that helps researchers, especially those with limited coding skills, to easily set up and configure their raspberry pi to run large numbers of controlled and automated image and video recordings. `pirecorder` was developed to over come this need.

The package consists of a main `PiRecorder` module to run recordings, `stream` and `camconfig` modules for help setting up, calibrating, and configuring the camera, a `schedule` module for scheduling future recordings, and a `convert` module for the easy converting of (folders of) recorded images and videos.

A core component of *pirecorder* is that it uses configuration files and timeplans that can be easily called, modified, and stored by the user, with automatic naming of files and folders. *pirecorder* also works directly from the terminal without the need to code in Python (see [this page](8-run-from-commandline.md) and comes with detailed documentation and tutorials (this website). This also has the aim to further help people with limited coding experience to set up their raspberry pi and make controlled and automated recordings.


## Table of contents
{: .no_toc .text-delta .fs-4 .fw-300 }

1. TOC
{:toc}
---

## Quick install

pirecorder can be easily installed with pip, which should already be automatically installed on your system. To install the latest release, simply open a terminal window and enter:

```
pip install pirecorder
```

To install the latest development version, enter:

```
pip install git+https://github.com/jollejolles/pirecorder.git --upgrade
```

See the [quick usage guide](quick-guide.md) for quickly getting you up and running or the [setting-up your raspberry pi](1-setting-up-raspberry-pi.md) and [installing pirecorder](2-installing-pirecorder.md) pages for more in-depth documentation and tutorials.

## Dependencies
*pirecorder* is both Python 2.7 and 3 compatible. It builds strongly on the [picamera](http://picamera.readthedocs.io/) package, uses [numpy](http://www.numpy.org/), [pyyaml](https://pyyaml.org), and [opencv](http://opencv.org) for some of its core functionality, and relies on various utility functions of my [pythutils](https://github.com/JolleJolles/pythutils) package. The scheduling functionality is based on *CronTab* and the associated [python-crontab](https://pypi.org/project/python-crontab/) package.

All dependencies are automatically installed with *pirecorder* except for:
* *OpenCV*: has to be manually installed due to various dependencies on the raspberry pi. Click [here](other/install-opencv.md) for a quick install guide.
* *FFmpeg*: is only needed for the convert functionality of *pirecorder*, so if you plan on using that functionality you should make sure it is installed. Click [here](other/install-ffmpeg-raspberry-pi.md) for my guide to install it on raspberry pi and [here](other/install-ffmpeg-osx.md) for my guide to install it on OS X.

## PiRecorder module
The main functionality of *pirecorder* is the `PiRecorder` module. This class initiates a PiRecorder instance that sets up the raspberry pi to record either A) a single image, b) a sequence of images, C) a single video, or D) a session of video recordings. When `PiRecorder` is run for the first time, it creates a "pirecorder" setup directory in the user's home folder to store all relevant setup files. This includes the default configuration file (*pirecorder.conf*) with all recording and camera settings as well as a log file (*pirecorder.log*) that will store all terminal output when `PiRecorder` is used to help keep a history log of your recordings.

### Configuration file
*PiRecorder* is set up in such a way that it is very easy to set and save custom recording and camera settings that are then automatically used for future recordings without further user input. Multiple configuration files can be created and called for specific recording settings and the configuration file(s) can be easily edited with any text editor as well as updated from the command line with the `settings` function.

A large number of custom recording parameters can be set, divided into 1) general user recording parameters, 2) camera settings, 3) video recording settings, 4) image recording settings, and 5) custom settings. A detailed overview and description of all the configuration settings can be found by calling `print(pirecorder.PiRecorder.settings.__doc__)` in Python and is explained in-depth in the [configure recording settings](4-configure-recording-settings.md) and the [configure camera settings](5-configure-camera-settings.md) guides.

### Recording modes
There are four recording modes, which can be set with the `rectype` parameter:

1. `img`: Records a single image with the custom settings.
2. `imgseq`: Creates a controlled sequence of images (i.e. timelapse) based on A) the set duration (`imgtime`) and B) the set total number of images to be recorded (`imgnr`) and the provided time delay between images (`imgwait`).
3. `vid`: Records a single video. Specific settings that can be set for this mode are `vidfps` (the framerate of the video), `vidduration` (the duration of the video), and `viddelay` (extra recording time in seconds that will be added to vidduration).
4. `vidseq`: Starts a series of standardized videos using the custom settings whereby, after each recording has ended, the user is asked if a new recording should be started or the program should exit.

### Automatic file naming
Files are automatically stored in the configured directory (`recdir`), by default a directory called `recordings` in the pirecorder directory, and named according to the provided `label`, the computer name, the date and time, and potentially the session number or image sequence number.

## Other modules
In addition to the main recording module, *pirecorder* contains a number of modules to facilitate setting-up and configuring the raspberry pi camera, schedule future recordings, and convert recorded media:

- `stream`: Opens a live video stream with user interface to calibrate the raspberry pi camera in terms of its position, focus, and region of interest (roi). For more detail, see the [calibrate camera page](3-position-and-calibrate-camera.md).
- `camconfig`: Opens a live video stream with user interface to dynamically, both manually and automatically, set the camera settings, including camera rotation, shutterspeed, whitebalance, iso, exposure compensation, brightness, contrast, saturation, and sharpness. For more detail, see the [configure camera settings page](5-configure-camera-settings.md)
- `schedule`: Automatically start image and video recording in the future according to custom recording schedules. For more detail, see the [schedule recordings page](6-recording-and-scheduling.md).
- `convert`: Convert (folders of) images or videos to videos with the option to resize, add timestamps on each frame, and monitor folders for automatic conversion. For more detail, see the [convert media page](7-convert-media.md).

## Development
*pirecorder* is developed by Dr Jolle Jolles, a research fellow at the Max Planck Institute of Animal Behavior, and at the Zukunftskolleg, Institute of Advanced Study at the University of Konstanz. For more information about his work, see his [academic website](http://jollejolles.com) or his [google scholar profile](https://scholar.google.nl/citations?user=VCZqbK4AAAAJ).

For an overview of version changes see the [CHANGELOG](https://github.com/jollejolles/pirecorder/blob/master/CHANGELOG) and for detailed changes see the [commits page](https://github.com/jollejolles/pirecorder/commits/).

Please submit bugs or feature requests to the GitHub issue tracker [here](https://github.com/jollejolles/pirecorder/issues).

## Citing
If you use pirecorder in your research, please cite it as follows:

```
@misc{Jolles2019,
      title = {pirecorder: controlled and automated image and video recording with the raspberry pi},
      author = {Jolles, Jolle W.},
      year = {2019}
      url = {http://doi.org/10.5281/zenodo.2529515},
      doi = {10.5281/zenodo.2529515}
}
```

## License
Released under a Apache 2.0 License. See [LICENSE](https://github.com/jollejolles/pirecorder/blob/master/LICENSE) for details.
---
layout: default
title: Home
nav_order: 0
permalink: /
---

# pirecorder documentation
{: .fs-9 }

*pirecorder* is a package to facilitate controlled and automated image and video recording for the raspberry pi, specifically developed with the biological sciences in mind.
{: .fs-6 .fw-300 }

![pirecorder logo](/pirecorder/assets/images/pirecorder-logo-large.jpg)

This website provides detailed documentation for setting-up your raspberry pi and installing and working with *pirecorder*. The guides provided here should also enable you to start working with your raspberry pi and making automated recordings, even if you have very limited knowledge of the raspberry pi and/or Python.

In addition to the documentation on this website, you can find a couple Jupyter notebooks in the [pirecorder repository](https://github.com/JolleJolles/pirecorder/tree/master/notebooks/) with all the relevant commands to start, configure, and run pirecorder, schedule recordings, and convert media. Jupyter is great to run Python interactively in your browser. To install it, simply open a terminal window and enter: `python -m pip install jupyter`, and then type `jupyter notebook` in Terminal to start an instance of Jupyter.

If you have issues with installing or working with *pirecorder* or have any feature requests, please submit them to the GitHub issue tracker [here](https://github.com/JolleJolles/pirecorder/issues).
---
layout: page
title: 6 Record and schedule recordings
nav_order: 8
---
# Record and schedule recordings=
{: .no_toc }

With the pirecorder package it is very easy to make recordings as well as schedule recordings in the future. You can follow the below guide directly or first configure your [recording](4-configure-recording-settings.md) and [camera](5-configure-camera-settings.md).
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta .fs-4 .fw-300 }

1. TOC
{:toc}
---

## Recording
To start your recordings, first import pirecorder, start your instance with the right configuration file, and use the `record` function:

```
import pirecorder
rec = pirecorder.PiRecorder()
rec.record()
```

This will use all recording and camera settings as detailed in your configuration file. It is possible to run recordings immideatealy without any configuration. Then it will just use the default (`configfile = "pirecorder.conf"`) and use the automatic mode to dynamically set the right shutterspeed and white balance. You may wish to set the `rectype` and related parameters though to get the recording that you want. Follow the documentation [here](4-configure-recording-settings.md).

It is also possible to run recordings straight from the terminal without requiring any further user input using the `record` command, which makes it very easy to run controlled recordings without requiring any user input. It will use the custom settings as provided in the configuration file, which you can change with the `--configfile` parameter:

```
record --configfile "custom.conf"
```

## Schedule recordings
Besides starting recordings directly, it is possible to schedule recordings to start recordings (repeatedly) in the future. For this there is the `schedule` function, which creates unique recording jobs (`jobname`) with specific `timeplan`s. An overview with a concise description of all parameters can be found at the bottom of this page.

When using the schedule function make sure you don't have multiple recording jobs at the same time as the camera can ofcourse only deal with a single recording at a time. Also scheduling is not possible with the "vidseq" rectype as that option waits for user input between videos.

### Set the timeplan for your future recordings
The `timeplan` parameter expects a code string that is build on [CRON](https://en.wikipedia.org/wiki/Cron) and should consists of the following parts:

```
* * * * *
- - - - -
| | | | |
| | | | +----- day of week (0 - 7) (sunday = 0 or 7)
| | | +---------- month (1 - 12)
| | +--------------- day of month (1 - 31)
| +-------------------- hour (0 - 23)
+------------------------- min (0 - 59)
```

Each of these parts supports wildcards (\*), ranges (2-5), and lists (2,5,6,11). For example, if you want to schedule a recording at 22:00, every workday of the week, enter the code `0 22 * * 1-5`. The minimum possible time between subsequent scheduled recordings is one minute. Smaller intervals between recordings is of course possible for images with the imgseq command with the `Record` method.

It is important to make sure that the PiRecorder configuration timing settings are within the timespan between subsequent scheduled recordings based on the provided timeplan. For example, a video duration of 20 min and a scheduled recording every 15 min between 13:00-16:00 (`*/15 13-16 * * *`) will fail. This will be checked automatically.

### Test your timeplan
To test a timeplan before linking it to a job, simply set the `test` parameter to `True`. For example,

```
rec.schedule(timeplan = "*/10 */2 10-15 * *", test = True)
```

will state "Your timeplan will run Every 10 minutes, every 2 hours, between day 10 and 15 of the month".

### Plan a recording job
When you have a correct timeplan that works as desired, you can link it to a new or existing job with the `jobname` parameter. Note that by default when creating a new job it won't be enabled immediately. For example, to create a new job named "rec1" :

```python
rec.schedule(timeplan = "*/10 */2 10-15 * *", jobname = "rec1")
```

### Enable/disable a job
To enable/disable a job, enter the `jobname` parameter and set the `enable` parameter to either `True` or `False` (the default). For example, to enable an existing job named "rec1":

```
rec.schedule(jobname = "rec1", enable = True)
```

### Show all jobs
It is easy to see all existing jobs and if they are disabled or when the next recording will start. To show this information, simply set the `showjobs` parameter to `True`:

```
rec.schedule(showjobs = True)
```

would for example show:

```
Job       Timeplan             Next recording
---------------------------------------------
rec1      */10 */2 10-15 *     disabled
rec2      * 3 * *              2020-06-01 03:00:00
```


### Remove jobs
Using the `schedule` function it is also easy to remove jobs. For example, to remove a specific recording job named "rec1" enter:

```
rec.schedule(jobname = "rec1", clear = "job")
```

And to remove all jobs:

```
rec.schedule(jobname = "rec1", clear = "all")
```

### Schedule recordings from the command line
Like for the other modules it is also possible to directly schedule recordings from the command line! You can enter all parameters like explained above by adding to the `schedule` command in terminal. For example:

```
schedule --timeplan "20 8 * * *" --test True
schedule --jobname "rec1" --timeplan "20 8 * * *" --enable False
schedule --jobname "rec1" --enable True
schedule --showjobs True
```

Entering just `schedule` without any additional parameters will just show an overview of all current jobs and their status.

---
Schedule module documentation
{: .text-delta .fs-5}
```
Parameters
----------
jobname : str, default = None
    Name for the scheduled recorder task to create, modify or remove.
timeplan : string, default = None
    Code string representing the time planning for the recorder to run
    with current configuration set. Build on CRON, the time plan should
    consist of the following parts:
    * * * * *
    - - - - -
    | | | | |
    | | | | +----- day of week (0 - 7) (sunday = 0 or 7)
    | | | +---------- month (1 - 12)
    | | +--------------- day of month (1 - 31)
    | +-------------------- hour (0 - 23)
    +------------------------- min (0 - 59)
    Each of the parts supports wildcards (*), ranges (2-5), and lists
    (2,5,6,11). For example, if you want to schedule a recording at
    22:00, every workday of the week, enter the code '0 22 * * 1-5' If
    uncertain, crontab.guru is a great website for checking your CRON
    code. Note that the minimum time between subsequent scheduled
    recordings is 1 minute. Smaller intervals between recordings is
    possible for images with the imgseq command with the Record method.
enable : bool, default = None
    If the scheduled job should be enabled or not.
showjobs : bool, default = False
    If the differently timed tasks should be shown or not.
delete : [None, "job", "all"], default = None
    If a specific job ('job'), all jobs ('all') or no jobs (None)
    should be cleared from the scheduler.
test : bool; default = False
    Determine if the timeplan is valid and how often it will run the
    record command.
configfile : str, default = "pirecorder.conf"
    The name of the configuration file to be used for the scheduled
    recordings. Make sure the file exists, otherwise the default
    configuration settings will be used.

Note: Make sure Recorder configuration timing settings are within the
timespan between subsequent scheduled recordings based on the provided
timeplan. For example, a video duration of 20 min and a scheduled
recording every 15 min between 13:00-16:00 (*/15 13-16 * * *) will fail.
This will be checked automatically.
```
---
layout: page
title: 8 Run from command line
nav_order: 10
---
# Run pirecorder from command line

When pirecorder is installed, it comes with a added functionality to run of all modules straight from the command line for ease of use.
{: .fs-6 .fw-300 }

All the same parameters can be set as when using pirecorder in Python, so see the respective documentation pages for a detailed explanation. Below just the bash commands are given with all parameters for the different modules.

### Camera stream
```
stream --cameratype "v2" --framerate 8 --vidsize 0.2 --rotation 0 \
       --imgoverlay "/home/pi/overlay.jpg"
```

### Camera configuration Stream
```
camconfig --auto True --framerate 20 --iso 200 --res (1640,1232) --vidsize 0.3
```

### Recording
```
record --configfile "pirecorder.conf"
```

### Scheduling
```
schedule --jobname None --timeplan "* * * * *" --enable True --showjobs False \
         --delete "job" --test True --configfile "pirecorder.conf"
```

### Converting
```
convert --indir VIDEOS --outdir CONVERTED --type ".h264" --withframe True \
        --pools 4 --resizeval 0.5 --sleeptime 5 --delete False
```
---
layout: page
title: Contact
nav_order: 12
---
# Contact

Hi, I am Jolle Jolles, the developer of *pirecorder*. I am a research fellow at the Max Planck Institute of Animal Behavior, Konstanz and at the Zukunftskolleg, Institute of Advanced Study at the University of Konstanz. For more information about my research, see my [academic website](http://jollejolles.com) and my [google scholar profile](https://scholar.google.nl/citations?user=VCZqbK4AAAAJ).
{: .fs-5 .fw-300 }

If you need further help running *pirecorder* you can email me [j.w.jolles@gmail.com](mailto:j.w.jolles@gmail.com). I also always like hearing from you when *pirecorder* has been useful for you! For bugs or feature requests it is easiest if you post an issue on the GitHub issue tracker [here](https://github.com/jollejolles/pirecorder/issues).
{: .fs-5 .fw-300 }
---
layout: page
title: 3 Position and calibrate the camera
nav_order: 5
---

# Calibrate the camera
{: .no_toc }

If pirecorder is successfully installed on your raspberry pi, the next thing you may want to do is properly position and calibrate your camera. This page will explain the various ways in which pirecorder can help with that. Note that the stream functionality also works on non-raspberry pi systems.
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta .fs-4 .fw-300 }

1. TOC
{:toc}
---

## Overview of the pirecorder stream functionality
pirecorder comes with a `stream` module that displays an interactive live video stream on the desktop to help position the raspberry pi camera and objects in the camera view, draw a region of interest to be used for recordings, and zoom in on part of the video. It records the clicks and movements of the mouse and responds to the following keypresses:

- `f`-key: Display the stream fullscreen/display the stream in a window
- `c`-key: Display/hide a diagonal cross across the screen
- `s`-key: Save the coordinates of the rectangular area when drawn
- `e`-key: Erase the rectangular area when drawn
- `z`-key: Show a zoomed-in section of the video inside the rectangular area in maximum resolution
- `n`-key: Refresh the zoom-in image
- `o`-key: If the potential overlay image should be shown or not
- `[`- and `]`-keys: Decrease or increase the relative opacity of the potential overlay image with 5%
- `esc`-key: Exit the the zoom window as well as the calibrate function completely

The `stream` module can be run independently:
```
import pirecorder
pirecorder.Stream()
```

as part of a `PiRecorder` instance:
```
import pirecorder
rec = pirecorder.PiRecorder()
rec.stream()
```

and directly from the command line:
```
stream
```

When running stream directly from the command line it is also possible to add a configuration file parameter such that when a rectangle is drawn, its coordinates will be stored as the ROI in the configfile (e.g. `stream --configfile pirecorder.conf`).

## Adjust the stream settings
There are a 5 parameters you can set for the video stream:
- `cameratype`: By providing this parameter the maximum resolution of the raspberry pi camera will be automatically determined. The default is "v2". Other possible values are "v1" and "hq".
- `vidsize`: Sets the size of the stream display window proportional to the maximum resolution of the provided camera type. The smaller the vidsize parameter, the more likely the video stream will run smoothly.
- `framerate`: By default the video stream will have a framerate of 8fps. A lower framerate may be desirable on older raspberry pi models. The maximum possible framerate is 5 when using an image overlay.
- `rotation`: This parameter enables showing the video stream in normal orientation (`0`) or up-side-down (`180`).
- `imgoverlay`: This parameter makes it possible to show an image overlay on the live video stream, which is further detailed below.

For example, to show a video stream in a relatively small window at a framerate of 15fps up-side down for a raspberry pi with a "v1" camera connected:

```
pirecorder.Stream(cameratype="v1", vidsize=0.2, framerate=15, rotation=180)
```

## Overlay an image on the video stream
It is also possible to overlay an image on the video stream. This can be highly beneficial to help position the camera exactly like in a previous recording or arrange objects in the field of view of the camera exactly like before, such as a custom maze for behavioural experiments.

To overlay an image, simply add the path to the image file to the `imgoverlay` parameter. The image will be stretched to have the same dimensions as the video stream. By default it will be shown with an opacity of 0.5 (range 0-1). The opacity can be reduced and increased interactively with the `[` and `]` keys.

## Positioning the camera
The stream function comes with the option to display a simple diagonal white cross, which can help to accurately position the camera or objects of interest below the camera. Simple toggle the cross with the `c`-key.

## Setting the region of interest
Clicking and drawing the mouse will draw a rectangular area. This can be used directly to store the coordinates of the region of interest that should be used for recordings, i.e. to only record the region inside the rectangular area.

Simply draw and redraw the rectangular area until you are happy and press the `s`-key. Now the coordinates of the region of interest will be displayed and, if running the `stream` functionality from a `PiRecorder` instance, stored automatically in the configuration file (e.g. `rec.config.cus.roi`). If you stored the region of interest accidentally or want to remove the drawn rectangle simple enter the `e`-key. To exit press the `esc`-key.

## Show a zoomed-in region at full resolution
Besides drawing the rectangular area for creating a region of interest for recordings, it can also be used to zoom-in on part of the video at the maximum image resolution. This can be very helpful to help you improve the focus of the raspberry pi camera. To do so, simply draw a rectangular area around the region you want to zoom-in to with the mouse and when satisfied press the `z`-key. Now a second window will open that will show the zoomed-in region of the video as an image.

As this is at the maximum image resolution, which is a lot higher than the maximum video resolution, it will only show a static image, but thus with much more detail than the video stream. You can refresh the image by pressing the `n`-key. To exit, press the `esc`-key.

---
Stream module documentation
{: .text-delta .fs-5}

```
Opens a video stream with user interface to help position and
adjust the camera

parameters
-----------
system : str, default = "auto"
    If the system should be automatically determined. Should detect if
    the computer is a raspberry pi or not. If this somehow fails, set
    to "rpi" manually.
framerate : int, default = 8
    The framerate of the displayed video stream. Lower framerates take
    longer to start up. When using an image overlay, maximum possible
    framerate is 5, to avoid getting shutter effects
vidsize : float, default = 0.2
    The relative size of the video window to the maximum resolution of
    the raspberry pi camera type.
rotation : int, default = 180
    If the camera should be rotated or not. 0 and 180 are valid values.
cameratype : str, default = "v2"
    The raspberry camera type used. Should be either "v1", "v2", or "hq"
imgoverlay : str, default = None
    The path to an image that will be overlaid on the video stream. This
    can be helpful for accurate positioning of the camera in line with
    previous recordings or setups.

interface
-----------
This interactive module stores mouse position and clicks and responds to
the following keypresses:
f-key : display the stream fullscreen/display the stream in a window
c-key : display/hide a diagonal cross across the screen
s-key : save the coordinates of the rectangular area when drawn
e-key : erase the rectangular area when drawn
z-key : show a zoomed-in section of the video inside the rectangular
    area in maximum resolution
n-key : refresh the zoom-in image
o-key : if the potential overlay image should be shown or not
[- and ]-keys : decrease or increase the relative opacity of the
    potential overlay image with 5%
esc-key : exit the the zoom window; exit the calibrate function
```
---
layout: page
title: 5 Configure camera settings
nav_order: 7
---

# Configure camera settings
{: .no_toc }

When you have pirecorder running and are happy with the recording settings, you may want to further configure the camera settings to get your optimal recordings. Below I explain the various ways in which you can do that very easily with pirecorder.
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta .fs-4 .fw-300 }

1. TOC
{:toc}
---

## Camera settings
It is possible to set a large number of camera settings with pirecorder: `rotation`, `gains`, `brightness`, `contrast`, `saturation`, `iso`, `sharpness`, `compensation`, and `shutterspeed`. A concise description all these parameters can be found at the bottom of this page. Below I go through the details of how to use the different parameters. As explained on the [configure recording settings page](4-configure-recording-settings.md), you can set them manually by editing your configuration file or when in Python and having a PiRecorder instance with the `settings` function.

To change the orientation of the camera set the `rotation` parameter either to `0` (normal) or `180` (upside down); to change the white balance set the `gains` parameter. This should be a tuple of blue and red gains values between 0 and 8, e.g. `(1.5, 1.85)`; to film in black and white set `saturation` to `-100`; to change the differences in luminance and color and get a more contrasting image increase the `contrast` parameter; to increase the sharpness of your recordings increase the `sharpness` parameter from a value from `-100` to `100`.

To get the optimal brightness level for your recordings you can tweak a number of parameters. For example, to get a brighter image, you could increase the `compensation`, increase the `brightness`, increase the `iso` or set a longer `shutterspeed`. Changing the compensation and iso values will change the amount of light that is allowed to enter the camera before recordings, which will affect the noise in the image. Brightness on the other hand is done in the post processing by stretching/compressing the range of values.

It is best to set `iso` to a low value and `compensation` to `0` to reduce image noise and keep `brightness` at its default intermediate value. Then if the image is too dark and image noise is not a problem for you, then first increase the `iso` and `compensation` before increasing the `brightness` parameter.

The `shutterspeed` parameter is important to get consistent images and can therefore be manually set. In many cases it may be desirable to not get motion blur in the image. In that case it is important to first determine the required shutterspeed, and then adjust the amount of light available, such as for experiments. If that is not possible, then start tweaking the other parameters as explained above. In some cases a very long shutterspeed may be needed (e.g. >1s), which may result in very blurry images but still allows objects to be detected and tracked for example.

## Automatic shutterspeed and whitebalance
By default, PiRecorder automatically sets the shutterspeed and white balance dynamically for and during each recording (`automode = True`). However, to get consistent recordings, it may be preferable to set the shutterspeed and white balance at a fixed value. This can be done both automatically and interactively (explained below) with pirecorder.

To get the optimal shutterspeed and white balance for the current conditions automatically you can use the `autoconfig` function, which will directly update the configuration file. This function will use the framerate provided in the configuration file for calibration, so make sure that is set properly. Then to run autoconfig, simply enter:

```
rec.autoconfig()
```

## Change the camera settings interactively
pirecorder also comes with a very handy interactive tool (`camconfig`) that enables you to set the camera settings dynamically. `camconfig` opens a live video stream and a separate window with a trackbar for each of the camera settings. You can slide your parameters of interest between the possible values and see live how the resulting recording will look like. To run camconfig and store the values automatically in your configuration file, use the function linked to your PiRecorder instance:

```
rec.camconfig()
```

You can exit the stream without saving with the `esc`-key and with saving with the `s`-key.

An important setting is the `automatic` mode. By default this is set to `True` such that it automatically gets the optimal shutterspeed and white balance (blue and red gains), which is visible by the respective trackbars sliding automatically to their optimal values. When you are relatively happy with these values it is a good time to use the non-automatic mode as then you are able to further tweak these values to your wishes.

`camconfig` will use the framerate as provided in the configuration file, but you can also dynamically update the framerate while the video stream is open, which in turn will influence the range of shutterspeeds possible (see above). It is also possible to use camconfig stand alone without a PiRecorder instance. Simply import pirecorder and run the function:

```
import pirecorder
pirecorder.Camconfig()
```

Then when the `s`-key is pressed it will output a dictionary of all values.

## Correct for raspberry-pi specific brightness
In some cases you may want to use multiple raspberry pi's to record the same scene. Due to small differences in the hardware or camera positioning there may be slight differences in the brightness of the recordings of the raspberry pi's. To correct for those you can use the `brighttune` parameter. Each raspberry pi will use the `brightness` parameter as default, so it can be easily set consistently across all raspberry pi's, but then for each specific raspberry pi you can tweak this value from -10 to +10.

---
Camera settings documentation
{: .text-delta .fs-5}
```
Parameters
---------------
rotation : int, default = 0
    Custom rotation specific to the Raspberry Pi, should be either 0 or 180.
gains : tuple, default = (1.0, 2.5)
    Sets the blue and red gains to acquire the desired white balance.
    Expects a tuple of floating values (e.g. "(1.5, 1.85)"). Can be
    automatically set with the autoconfig() function and interactively with the
    camconfig() function using a live video stream.
brightness : int, default = 45
    Sets the brightness level of the camera. Expects an integer value
    between 0 and 100. Higher values result in brighter images.
contrast : int, default = 20
    Sets the contrast for the recording. Expects an integer value between 0
    and 100. Higher values result in images with higher contrast.
saturation : int, default 0
    Sets the saturation level for the recording. Expects an integer value
    between -100 and 100.
iso : int, default = 200
    Sets the camera ISO value. Should be one of the following values:
    [100, 200, 320, 400, 500, 640, 800]. Higher values result in brighter
    images but with higher gain.
sharpness : int, default = 50
    Sets the sharpness of the camera. Expects an integer value between -100
    and 100. Higher values result in sharper images.
compensation : int, default = 0
    Adjusts the camera’s exposure compensation level before recording.
    Expects a value between -25 and 25, with each increment representing
    1/6th of a stop and thereby a brighter image.
shutterspeed : int, detault = 10000
    Sets the shutter speed of the camera in microseconds, i.e. a value of
    10000 would indicate a shutterspeed of 1/100th of a second. A longer
    shutterspeed will result in a brighter image but more motion blur.
    Important to consider is that the framerate of the camera will be
    adjusted based on the shutterspeed. At low shutterspeeds (i.e. above
    ~ 0.2s) the required waiting time between images increases considerably
    due to the raspberry pi hardware. To control for this, automatically a
    standard `imgwait` time should be chosen that is at least 6x the
    shutterspeed. For example, for a shutterspeed of 300000 imgwait should
    be > 1.8s.
brighttune : int, default = 0
    A rpi-specific brightness compensation factor to standardize light
    levels across multiple rpi"s, an integer between -10 and 10.    
```
---
layout: page
title: 1 Setting up your RPi
nav_order: 3
---

# Setting up your raspberry pi
{: .no_toc }

The below steps will guide you through some basic settings to get your raspberry pi set up for working with *pirecorder*. Most users may not need to follow these steps, but they are shown for completeness.
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta .fs-4 .fw-300 }

1. TOC
{:toc}
---

## Start up your raspberry pi
This guide presumes you have a raspberry pi with working version of Raspbian that is powered on with direct access (screen+keyboard+mouse) or VNC. SSH is also an option but then no direct view of the camera is possible and so is only recommended after initial setting up and calibration of the camera.

## Update your raspberry pi

Make sure that your raspberry pi is fully up to date. Open a terminal window and enter:

```
sudo apt-get update && sudo apt-get upgrade -y
sudo apt-get autoremove -y
```

Also make sure your firmware is fully up to date to get the latest drivers to work with the picamera
```
sudo apt-get dist-upgrade -y
```

## Enable and test the camera

To be able to use the camera we need to enable it in the configuration menu. Enter:

```
sudo raspi-config
```

go to `5 Interfacing options`, then `P1 Camera`, and click `yes`. Alternatively you can go to the main menu and use the Raspberry Pi Configuration tool there. Now reboot your pi if you have already set up your raspberry pi camera. If not, turn of your raspberry pi and now connect one part of the ribbon cable to the camera and the other end to the raspberry pi by pulling up the edges of the plastic clip of the camera module port and sliding in the ribbon cable, making sure the cable is the right way around.

You can test the camera quickly by entering the command `raspistill -t 0 -k` in a terminal window. To exit again, press `ctrl+c`. If you get an error message, then double check the cable is properly connected to the raspberry pi and the camera and try restarting.

## Set the screen resolution

For configuring and calibrating the camera and using the dynamic interfaces that come with pirecorder, a screen resolution of at least 800x600 is strongly recommended, and at least 1024x1024 is preferable.

The resolution can be set with `raspi-config`: option `7 Advanced options`, then `A5 Resolution`. When changing the resolution a restart will be required.

## Setup python for working with the camera

*pirecorder* uses Python, which comes pre-installed with Raspbian on any raspberry pi. However, we need to update some python development tools for the raspberry pi camera to work properly. To do this simply open a terminal window and enter:

```
sudo apt-get install python-setuptools python-dev build-essential libpq-dev
```

## Optional: Enable filesharing
For easy transfering of files with your raspberry pi it may be good to enable filesharing. There are various ways to share files, depending on your system. I recommend to use netatalk:

```
sudo apt-get install netatalk -y
```

Add the home directory to be shared:
```
sudo nano /etc/netatalk/afp.conf
```
And add the following text to the bottom of the file

```
[Homes]
  basedir regex = /home
```

Now restart the service:

```
sudo systemctl restart netatalk
```
---
layout: page
title: Quick usage guide
nav_order: 1
---

# Quick usage guide
{: .no_toc }

This guide is meant to get you up and running and using pirecorder in no time. We will go through all the steps from setting up your raspberry pi, installing pirecorder, configuring pirecorder, recording, scheduling, and finally converting your recorded media.
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta .fs-4 .fw-300 }

1. TOC
{:toc}

## Setting up your raspberry pi
1. Connect a camera with ribbon cable to the raspberry pi and roughly position it where you need it
2. Connect your raspberry pi to power to turn it on and either have direct control or VNC access
3. Open a terminal window and make sure your raspberry pi is fully up to date:
```
sudo apt-get update && sudo apt-get upgrade -y
```
4. Also make sure your firmware is fully up to date to get the latest drivers to work with the picamera
```
sudo apt-get dist-upgrade -y
```
5. Make sure the camera is enabled
```
sudo raspi-config
```
go to `5 Interfacing options`, then `P1 Camera`, and click `yes`. Now reboot your raspberry pi.

6. Make sure the camera is properly connected:
```
raspistill -t 0 -k
```
To exit again, press `ctrl+c`

## Install pirecorder
1. install pirecorder (and most dependencies) with pip:
```
pip install pirecorder
```
2. install dependencies for OpenCV:
```
sudo apt-get install libhdf5-dev libhdf5-serial-dev -y
sudo apt-get install libqtwebkit4 libqt4-test -y
sudo apt-get install libatlas-base-dev libatlas3-base libjasper-dev libqtgui4 python3-pyqt5 -y
```
3. Install OpenCV with pip:
```
pip install opencv-contrib-python==4.1.0.25
```
Note: we use a specific version as the latest version may not always work properly on raspberry pi yet.

## Use pirecorder for the first time
1. Start python (`python`) and import the pirecorder package
```
import pirecorder
```
2. Initiate a recording instance, which will automatically create the `/home/pi/pirecorder` folder and default configuration file:
```
rec = pirecorder.PiRecorder()
```
Note: the variable to store the `PiRecorder` instance can be any name, here `rec` is chosen as an example.

## Optionally: Correctly position the camera and adjust its focus
1. Open a stream instance:
```
rec.stream()
```
2. Change the position of the camera until satisfied. Optionally press the `c`-key to show a diagonal cross, which can help with positioning.
3. Drag with the mouse to create a rectangular area of a region you want to use to adjust the camera focus. Press the `z`-key to show this region zoomed in. Now adjust the focus by slightly turning the camera lens, press the `n`-key to refresh the image, and continue with these steps until you are satisfied. Press the `esc`-key to exit the zoomed in window.

## Optionally: Store the region of interest
Still in the stream instance, use the mouse again to draw a rectangle that encompasses the region of the camera stream that should be recorded. When satisfied with the region press the `s`-key to store the coordinates in the configuration file. Press the `esc`-key to exit the stream.

## Set the recording settings
Configure your recording settings with the `settings` function. Key is the `rectype` parameter, which can be `img` (a single image), `vid` (a single video), `imgseq` (a timelape), `vidseq` (a sequence of videos). By default it will record with the maximal resolution, which can be altered with `imgdims` and `viddims`. Videos are recorded with 24fps, which can be changed with the `vidfps` parameter. The `label` parameter will help identify your recorded files. Other relevant parameters are `vidduration`, `viddelay`, `imgnr`, `imgtime`, `imgwait`, `imgquality`, and `vidquality`. For example, to set the `rectype`, `vidduration`, and `label` for your recordings:
```
rec.settings(rectype = "vid", vidduration = 60, label = "test")
```

## Optionally: set the camera settings
When you have pirecorder running and are happy with the recording settings, you may want to further configure the camera settings to get your optimal recordings. A large number of parameters can be set. Have a look at the detailed documentation [here](6-configure-camera-settings.md). Very handy is the `camconfig` function with which it is possible to open an interactive video stream and adjust the camera settings with trackbars:
```
rec.camconfig()
```

## Start a recording
Now you are happy with your settings you can simply start a recording with the `record` function:
```
rec.record()
```

## Schedule a recording
To schedule recordings in the future you need to set a `timeplan` (build on CRON) and create jobs with unique `jobname`. These jobs can then be enabled or disable with the `enable` parameter and removed completely with the `delete` parameter (either enter the jobname or "all"). It is also possible to test timeplans with the `test` parameter.

An example. If your recording configuration is set to record a 15min video and you want to record a new video every half an hour between 07.00 and 17.00 on each week day:

```
rec.schedule(timeplan = "0/30 7-17 * * 1-5", jobname = "weekjob", enable = True)
```

## Convert media
Finally, when you have your recorded media you may want to convert it. Videos are automatically recorded in the `.h264` format and may need converting to `.mp4` to be easily viewable. If you have recorded a timelapse of images you may also want to convert that to a video. Or maybe you want to add framenumbers to the frames of a video to later accurately determine certain events. You can use the `convert` module for all these commands, which also works on non-raspberry pi systems and can be run straight from the terminal. For example, to convert a folder of videos:

```
convert --indir VIDEOS --outdir VIDEOS/CONVERTED
```
---
layout: page
title: 4 Configure recording settings
nav_order: 6
---
# Configure recording settings
{: .no_toc }

When you have pirecorder running and have positioned your camera, the next thing you may want to do is change the default recording settings. Below I guide you through the various options.
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta .fs-4 .fw-300 }

1. TOC
{:toc}
---

## Working with the configuration file
As is explained in the [pirecorder package page](pirecorder-package.md), one of the main features of PiRecorder is that it works with a simple-to-use configuration file to set all your camera and recording settings to easily run repeated and automated recordings.

The configuration file can be set with the `configfile` parameter when initiating the PiRecorder instance, which defaults to "pirecorder.conf". For example, to work with a special config file for infrared recordings, you can run:

```
import pirecorder
rec = pirecorder.PiRecorder(configfile = "irsettings.conf")
```

The configuration files can be simply changed manually with any text editor or with your favorite terminal editor (e.g. `nano pirecorder/pirecorder.conf`). You can also easily update the configuration via the `settings` function and only add the parameters you want to change. For example:

```
rec.settings(label="test", rectype="img")
```

An overview of all the recording parameters with concise description that are possible to set is provided at the bottom of this page and can also be called directly in Python:

```
print(rec.settings.__doc__)
```

To see your current configuration settings, simply type in:

```
print(rec.config)
```

## Where and what to record
Two important configurations to set are the `recdir` and `rectype` parameters. By default PiRecorder will store all media in a recordings folder inside the pirecorder folder (`pirecorder/recordings`). You can change this to any folder name that you like. If no name is provided, the files will be directly stored in the users' home directory. If you want to store media on a mounted media, name the mounted folder "NAS" and use `recdir = "NAS"`, then it will automatically check if it is mounted and otherwise record in the default location.

With PiRecorder you can record single images (`rectype = img`), a sequence/timelapse of images (`rectype = imgseq`), a single video (`rectype = vid`), or multiple sessions of videos (`rectype = vidseq`). The "vidseq" recording type will record multiple videos with the same recording settings but wait after each finished recording for user input to continue with the next recording or exit. Each new video will be treated as a new "session" and have a corresponding session number in its filename (e.g. "S01", "S02" etc). The benefit of this recording option is that it is even quicker to record multiple videos one after the other with the same parameters and to have a simple automatic filenaming system that keeps those videos together.

## Automatic naming of files
The naming of the (folders of) recorded images and videos is all done automatically. Each filename will include a user-provided label that can be set with the `label` parameter, the host computer name, the date, the time and potentially the session number or count nr. For example, single images will have a filename like `test_200601_pi13_101300.jpg`, image sequences `test_200601_pi13_img00231_101300.jpg` and videos `test_180312_pi3_S01_101300.h264`.

The `subdirs` parameter makes it possible to automatically store the files of each separate recording in its own folder with a unique filename. This is especially helpful when recording image sequences of video sequences, such that those files are all stored together in their own folder.

## Settings for image recording
To set the resolution for images you can use the `imgdims` parameter, which defaults to the maximum resolution for the v1.5 camera model (2592, 1944). The v2 model has a max resolution of 3280 x 2464 pixels, and the hq camera 4056 x 3040 pixels. The `imgquality` parameter specifies the quality that the jpeg encoder should attempt to maintain. Use values between 1 and 100, where higher values are higher quality. Playing with this setting can help to considerably reduce the file size of your recordings while keeping the same quality.

To control your image sequences you can set three parameters: `imgnr`, `imgtime`, and `imgwait`. PiRecorder will use the minimum of `imgnr` and the nr of images based on `imgwait` and `imgtime`. When the value provided for imgwait is too low relative to the provided shutterspeed it will be automatically set to the minimum value of 0.45s. With a fast enough shutterspeed it is possible to record multiple images per second, but depends on the model of raspberry pi you use. Also, when a delay is provided that is less than ~x5 the shutterspeed, the camera processing time will take more time than the provided imgwait parameter and so images are taken immediately one after the other. To take a sequence of images at the exact right delay interval the imgwait parameter should be at least 5x the shutterspeed (e.g. shutterspeed of 400ms needs imgwait of 2s.

For example, to record a sequence of 10 images at very high resolution at 1 image a minute:

```
rec.settings(imgwait = 1, imgnr = 10, imgtime = 15, imgwait = 60, imgquality = 90 \
             imgdims = (3280, 2464))
rec.record()
```

## Settings for video recording
For recording video you can set the `vidduration` and `viddelay` to get the right recording duration. The viddelay is extra recording time in seconds that will be added to vidduration. Its use is to add a standard amount of time to the video that can be easily cropped or skipped, such as for tracking, but still provides useful information.

To set the resolution for video recording use the `viddims` parameter. Note that the maximum video resolution for all currently existing raspberry pi camera's ("v1","v2" and "hq") is 1080p, but that it is also possible to record in a different format as long as the total number of pixels does not exceed that, such as 1640 x 1232.

To set the framerate of the video use the `vidfps` parameter. With smaller resolutions higher framerates are possible, see [this page](https://picamera.readthedocs.io/en/release-1.13/fov.html#camera-modes) for more information. 40fps with the max resolution of 1640 x 1232, and 90fps with 1280 x 720 is possible but may result in dropped frame, so it is safer to stay just slightly below that.

The `vidquality` parameter specifies the quality that the h264 encoder should attempt to maintain. Use values between 10 and 40, where 10 is extremely high quality, and 40 is extremely low.

For example, to take a single video for 10 minutes with 20s extra time, with a 1640x1232 resolution at 24fps, with a relatively low quality and thus file size:

```
rec.set_config(rectype = "vid", vidduration = 600, viddelay = 20, vidquality = 30 \
               viddims = (1640, 1232), vidfps = 24)
rec.record()
```

---
Recording settings documentation
{: .text-delta .fs-5}
```
Parameters
---------------
recdir : str, default = "pirecorder/recordings"
    The directory where media will be stored. Default is "recordings". If
    different, a folder with name corresponding to location will be created
    inside the home directory. If no name is provided (""), the files are
    stored in the home directory. If "NAS" is provided it will additionally
    check if the folder links to a mounted drive.
subdirs : bool, default = False
    If files of individual recordings should be stored in subdirectories
    or not, to keep all files of a single recording session together.
label : str, default = "test"
    Label that will be associated with the specific recording and stored in
    the filenames.
rectype : ["img", "imgseq", "vid", "vidseq"], default = "img"
    Recording type, either a single image or video or a sequence of images
    or videos.
cameratype : str, default = None
    The raspberry cameratype used. Can be either None, "v1", "v2", or "hq"
    to indicate the different models and will help set the maximum recording
    resolution.
imgdims : tuple, default = (2592, 1944)
    The resolution of the images to be taken in pixels. The default is the
    max resolution for the v1.5 model, the v2 model has a max resolution of
    3280 x 2464 pixels, and the hq camera 4056 x 3040 pixels.
viddims : tuple, default = (1640, 1232)
    The resolution of the videos to be taken in pixels. The default is the
    max resolution that does not return an error for this mode.
imgfps : int, default = 1
    The framerate for recording images. Will be set automatically based on
    the imgwait setting so should not be set by user.
vidfps : int, default = 24
    The framerate for recording video.
imgwait : float, default = 5.0
  The delay between subsequent images in seconds. When a delay is provided
    that is less than ~x5 the shutterspeed, the camera processing time will
    take more time than the provided imgwait parameter and so images are
    taken immideately one after the other. To take a sequence of images at
    the exact right delay interval the imgwait parameter should be at least
    5x the shutterspeed (e.g. shutterspeed of 400ms needs imgwait of 2s).
imgnr : int, default = 12
    The number of images that should be taken. When this number is reached,
    the recorder will automatically terminate.
imgtime : integer, default = 60
    The time in seconds during which images should be taken. The minimum of
    a) imgnr and b) nr of images based on imgwait and imgtime will be used.
imgquality : int, default = 50
    Specifies the quality that the jpeg encoder should attempt to maintain.
    Use values between 1 and 100, where higher values are higher quality.
vidduration : int, default = 10
    Duration of video recording in seconds.
viddelay : int, default = 0
    Extra recording time in seconds that will be added to vidduration. Its
    use is to add a standard amount of time to the video that can be easily
    cropped or skipped, such as for tracking, but still provides useful
    information, such as behaviour during acclimation.
vidquality : int, default = 11
    Specifies the quality that the h264 encoder should attempt to maintain.
    Use values between 10 and 40, where 10 is extremely high quality, and
    40 is extremely low.
```
---
layout: default
title: Other guides
nav_order: 11
has_children: true
permalink: /docs/other-guides
---

# Other guides

Here you can find additional guides that may be of help when installing and running pirecorder.
{: .fs-6 .fw-300 }
---
layout: page
title: Install opencv
parent: Other guides
nav_order: 10
---

# Installing opencv on raspberry pi, ubuntu, and Mac OS X
{: .no_toc }

Installing [OpenCV](https://opencv.org) has never been very simple, especially when it could only be build from source. This was especially painful when wanting to run it on a Raspberry Pi as building and installing OpenCV took such a lot of time, especially on the older models. Luckily as of last year it is possible to [install OpenCV with pip](https://pypi.org/project/opencv-python)!

Below I guide you through the basic steps necessary to get OpenCV to wor on Mac, Ubuntu and Raspberry Pi. If you want more background information, see the excellent article [here](https://www.pyimagesearch.com/2018/09/19/pip-install-opencv/) by Adrian Rosebrock from [pyimagesearch.com](http://PyImageSearch.com).

## Table of contents
{: .no_toc .text-delta .fs-4 .fw-300 }

1. TOC
{:toc}
---

## Install pip
Pip is the main package manager for python that we will also use to install OpenCV. Pip should already be installed on your system (see [here](https://pip.pypa.io/en/stable/installing/)), but if it's not, we can install it with wget. Open a Terminal window and enter:

```
wget https://bootstrap.pypa.io/get-pip.py
```

Now to install pip for Python 3 enter:

```
sudo python3 get-pip.py
```

## Install prerequisites
Next, for some versions of Raspbian we may need to install some additional packages. First make sure `apt-get` is fully up-to-date by entering the following in Terminal:

```
sudo apt-get update
```

Now install the prerequisites:

```
sudo apt-get install libhdf5-dev libhdf5-serial-dev -y
sudo apt-get install libqtwebkit4 libqt4-test -y
sudo apt-get install libatlas-base-dev libatlas3-base libjasper-dev libqtgui4 python3-pyqt5 -y
```

## Install OpenCV with pip
Finally, we can enter install OpenCV very simply with the command:

```
pip install opencv-contrib-python
```

However, before running above command, it is important to note that the latest version of opencv may not always be fully functional on the Raspberry Pi. Therefore I recommend to run the below command that installs the latest known working version:

```
pip install opencv-contrib-python==4.1.0.25
```

If you still get an error message such as <em>Could not find a version that satisfies the requirement opencv-contrib-python (from versions: ) No matching distribution found for opencv-contrib-python</em>, try the alternative to use apt-get instead of pip:

```
sudo apt-get install python-opencv
```

## Testing
Now let's just make sure that OpenCV is working. Open a terminal window and enter `python3` to start Python. Now to make sure you have installed OpenCV correctly enter:

```
import cv2
cv2.__version__
```

Your terminal window should look like:

```
$ python3
Python 3.7.3 (default, Dec 20 2019, 18:57:59)
[GCC 8.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import cv2
>>> cv2.__version__
'4.1.0'
```

You are done!

Note: an alternative to getting to run the latest version of opencv on Raspberry Pi is to run the following command rather than `python3`: `LD_PRELOAD=/usr/lib/arm-linux-gnueabihf/libatomic.so.1 python3`.
---
layout: page
title: Install ffmpeg on RPi
parent: Other guides
nav_order: 11
---

# Installing ffmpeg on raspberry pi with .h264 support
{: .no_toc }

The raspberry pi is great for recording video. One issue however is that the `.h264` container it records in is hard to work with. It is therefore often desirable to convert videos to widely applicable formats like `.mp4` to be able to view them properly and get the right meta information. For this I recommend the program `FFmpeg`.

The `pirecorder` package also comes with a special `Convert` class with a number of helpful functionalities to make it very easy to convert videos, folders of videos as well as images to video. See the guide [here](7-convert-media.md).

Installing [FFmpeg](https://www.ffmpeg.org) on a Raspberry Pi is not as simple as downloading an executable from the command line, but it is also not too difficult. Follow these steps:

## Table of contents
{: .no_toc .text-delta .fs-4 .fw-300 }

1. TOC
{:toc}
---

## Install h264 library
Open a terminal window on the raspberrypi (or via SSH connection) and type in the following commands:
* Download h264 library:
```
git clone --depth 1 https://code.videolan.org/videolan/x264
```
* Change directory to the x264 folder:
```
cd x264
```
* Configure installation:
```
./configure --host=arm-unknown-linux-gnueabi --enable-static --disable-opencl
```
* Create the installation:
```
make -j4
```

* Install h264 library on your system:
```
sudo make install
```

## Install ffmpeg with h264
* Change to home directory:
```
cd ~
```
* Download ffmpeg:
```
git clone git://source.ffmpeg.org/ffmpeg --depth=1
```
* Change to ffmpeg directory:
```
cd ffmpeg
```
* Configure installation:
```
./configure --extra-ldflags="-latomic" --arch=armel --target-os=linux --enable-gpl --enable-omx --enable-omx-rpi --enable-libx264 --enable-nonfree
```
* Make the installation:
```
make -j4
```
*Note this step may take a long time!*
* Now finally run the installation:
```
sudo make install
```

There are many options available and many other ways to convert h264 videos with ffmpeg, but this command is the quickest of all methods that I tested.

Note: If you are working with an older model of the raspberrypi (&lt; 3 B+) then you may not have 4 cores available. You will then have to change `make -j4` to `make -j`.

## Convert (h264) video
Now you are ready to convert (h264) videos on your Raspberry Pi. To convert a single video with ffmpeg:

```
ffmpeg -i USER_VIDEO.h264 -vcodec copy USER_VIDEO.mp4
```

The `Convert` class of the pirecorder package builds on ffmpeg and added a number of functionalities to make it easier to help you convert your media, especially as recorded with `pirecorder`. Read its documentation [here](https://github.com/JolleJolles/pirecorder/wiki/pirecorder-convert/). For example, to convert a folder of videos, add frame numbers to the topleft corner of each video frame, and resize the video by half:

```
convert --indir VIDEOS --outdir CONVERTED --withframe True --resizeval 0.5
```
---
layout: page
title: Install ffmpeg on OSX
parent: Other guides
nav_order: 12
---

# Installing ffmpeg on Mac OS X
FFmpeg is a great program that you can run from the command line to help convert more or less any media format. This short guide will help you install ffmpeg on Mac.

The easiest way to install ffmpeg is to use [HomeBrew](https://brew.sh) a package manager for Mac. If you don’t have homebrew installed on your mac already, run the following command using terminal:

```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

Once you have Homebrew installed, you can now simply install ffmpeg from terminal with the following command:

```
brew install ffmpeg
```

To install ffmpeg with specifical modules, instead of running the above command run below command or remove those modules you do not need:

```
brew install ffmpeg --with-chromaprint --with-fdk-aac --with-fontconfig --with-freetype --with-frei0r --with-game-music-emu --with-libass --with-libbluray --with-libbs2b --with-libcaca --with-libgsm --with-libmodplug --with-librsvg --with-libsoxr --with-libssh --with-libvidstab --with-libvorbis --with-libvpx --with-opencore-amr --with-openh264 --with-openjpeg --with-openssl --with-opus --with-rtmpdump --with-rubberband --with-sdl2 --with-snappy --with-speex --with-tesseract --with-theora --with-tools --with-two-lame --with-wavpack --with-webp --with-x265 --with-xz --with-zeromq --with-zim
```
