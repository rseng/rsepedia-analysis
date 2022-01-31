![Raspberry Shake event logo](https://raw.githubusercontent.com/raspishake/rsudp/master/docs/_static/logo.png)
# rsudp
### Continuous sudden motion and visual monitoring of Raspberry Shake data
*Written by Ian Nesbitt (@iannesbitt), Richard Boaz, and Justin Long (@crockpotveggies)*

[![PyPI](https://img.shields.io/pypi/v/rsudp)](https://pypi.org/project/rsudp/)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/rsudp)](https://pypi.org/project/rsudp/)
[![GitHub](https://img.shields.io/github/license/raspishake/rsudp)](https://github.com/raspishake/rsudp/blob/master/LICENSE)
[![Documentation](https://img.shields.io/badge/docs-passed-brightgreen)](https://raspishake.github.io/rsudp/)
[![Build Status](https://scrutinizer-ci.com/g/raspishake/rsudp/badges/build.png?b=master)](https://scrutinizer-ci.com/g/raspishake/rsudp/build-status/master)
[![Code Coverage](https://scrutinizer-ci.com/g/raspishake/rsudp/badges/coverage.png?b=master)](https://scrutinizer-ci.com/g/raspishake/rsudp/?branch=master)
[![Scrutinizer Code Quality](https://scrutinizer-ci.com/g/raspishake/rsudp/badges/quality-score.png?b=master)](https://scrutinizer-ci.com/g/raspishake/rsudp/?branch=master)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02565/status.svg)](https://doi.org/10.21105/joss.02565)

`rsudp` is a tool for receiving and interacting with data casts from [Raspberry Shake](https://raspberryshake.org) personal seismographs and Raspberry Boom pressure transducer instruments.

`rsudp` has [full documentation here](https://raspishake.github.io/rsudp/). We also have [tutorial instructions](https://raspishake.github.io/rsudp/index.html#tutorial) to install, set up, and run `rsudp` there. Additionally, our documentation features [YouTube walkthroughs](https://raspishake.github.io/rsudp/youtube.html), [notes for contributors](https://raspishake.github.io/rsudp/contributing.html), a brief [Developer's guide](https://raspishake.github.io/rsudp/theory.html), and [module documentation](https://raspishake.github.io/rsudp/#code-documentation).

We now have [a paper](https://doi.org/10.21105/joss.02565) published in The Journal of Open Source Software! You can reference `rsudp` using the following citation:

> Nesbitt et al., (2021). rsudp: A Python package for real-time seismic monitoring with Raspberry Shake instruments. Journal of Open Source Software, 6(68), 2565, https://doi.org/10.21105/joss.02565


`rsudp` contains ten main features:
1. **Alert** - an earthquake/sudden motion alert trigger, complete with a bandpass filter and stream deconvolution capabilities
2. **AlertSound** - a thread that plays a MP3 audio file in the event of the alert module signalling an alarm state
3. **Plot** - a live-plotting routine to display data as it arrives on the port, with an option to save plots some time after an alarm
4. **Tweeter** - a thread that broadcasts a Twitter message when the alert module is triggered, and optionally can tweet saved plots from the plot module
5. **Telegrammer** - a thread similar to the Tweeter module that sends a [Telegram](https://telegram.org) message when an alarm is triggered, which can also broadcast saved images
6. **Writer** - a simple miniSEED writer
7. **Forward** - forward a data cast to one or several IP/port destinations
8. **RSAM** - computes RSAM (Real-time Seismic AMplitude) and either prints or forwards it to an IP/port destination
9. **Custom** - run custom code when an `ALARM` message is received
10. **Print** - a debugging tool to output raw data to the command line

`rsudp` is written in Python but requires no coding knowledge to run. Simply follow the [instructions to install the software](https://raspishake.github.io/rsudp/installing.html), go to your Shake's web front end, [configure a UDP datacast](https://manual.raspberryshake.org/udp.html#configuring-a-data-stream-the-easy-way) to your computer's local IP address, [start the software](https://raspishake.github.io/rsudp/running.html) from the command line, and watch the data roll in.

![Earthquake plot recorded on a Raspberry Shake 4D](https://raw.githubusercontent.com/raspishake/rsudp/master/docs/_static/4d-event.png)

(Above) a plot of an earthquake on the four channels of a Raspberry Shake 4D (EHZ---the geophone channel, and EHE, EHN, and ENZ---the accelerometer east, north, and vertical channels).
# Changelog
## changes in 1.1.0
- changing version tag to reflect peer review status, and creating Zenodo record

## changes in 1.0.3
- changed matplotlib pin to be `<3.2` rather than `==3.1.1` to address [#21](https://github.com/raspishake/rsudp/issues/21)
- modified logos
- fixed unicode output error (emoji caused error on Windows machines)
- added version printout to rs lib initialization sequence
- simplified dependency imports
- fixed an error with the packetization script that caused it to break on files where the first trace in a given stream was not the shortest
- adjusted default alert settings to be more in line with what testing says is optimal
- adding entrypoint convenence functions `rs-settings`, `rs-log`, and `rs-tailf` to make editing settings and monitoring log output easier
- minor changes to paper manuscript and bib file

## changes in 1.0.2
- corrected install script to fix [#14](https://github.com/raspishake/rsudp/issues/14)
- corrected social media URL destination [#15](https://github.com/raspishake/rsudp/issues/15)
- adding feature requested in [#9](https://github.com/raspishake/rsudp/issues/9) (additional text for twitter posts [documented here](https://raspishake.github.io/rsudp/settings.html#tweets-twitter-notification-module))
- edited language as requested in [#13](https://github.com/raspishake/rsudp/issues/13)
- added feature requested in [#17](https://github.com/raspishake/rsudp/issues/17) to control forwarding of `ALARM` and `RESET` messages as well as data
- added feature requested in [#18](https://github.com/raspishake/rsudp/issues/18) to forward messages to multiple destinations (this changes the syntax of the `"address"` and `"port"` fields of the `"forward"` settings section to lists)
- changed logging structure to be more downstream-friendly. downstream software can now initialize logging to `/tmp/rsudp/XYZ.log` by calling `rsudp.start_logging(logname='XYZ.log')`
- addressed ([#23](https://github.com/raspishake/rsudp/issues/23)) which prevented data from being written to the output directory
- added tests for alert, alertsound, consumer, custom, forward, plot, printraw, rsam, producer, packetize, Telegram, Twitter, and write modules as suggested [in review](https://github.com/openjournals/joss-reviews/issues/2565) ([#22](https://github.com/raspishake/rsudp/issues/22))
- added `.coveragerc` and code coverage basics

## changes in 1.0.1
- added `rsudp.c_rsam` Real-time Seismic Amplitude Measurement consumer
- modified install scripts for clarity
- added Windows batch scripts for installation, updates, and running, to match Unix ones

## changes in 1.0.0
- settings changed to deconvolve plot channels by default
- added the ability to post to multiple Telegram chats (by spinning up multiple independent threads)
- moved several functions to a new `helpers.py` module
- simplified several functions to make them more readable
- changed doc structure to github pages compatible

## changes in 0.4.3
- added ability to run tests with any data file containing at least one of `SHZ, E[H,N][E,N,Z], HDF` channels (even miniSEED, which gets converted to text first then read by the pre-producer)
- cut whitespace from the beginning of included MP3s
- added standardized queue message constructors to `rsudp.raspberryshake`
- removed warning filters
- fixed plot trace offset issue
- fixed a problem where UTC would appear after link in telegram and tweet messages
- fixed problem with precision in event `UTCDateTime` objects
- fixed unit capitalization in plot y-label
- added an exit code to the test function
- added a custom thread class (`rsudp.raspberryshake.ConsumerThread`) for consumers to inherit which contains all internal flags that the Producer needs to function
- added additional trove classifiers
- alarm time in plot, telegrams, and tweets now has 0.01 second precision
- alarm time now reports directly from `rsudp.c_alert.Alert` instead of Producer
- fixed a circular import issue which manifest on RPi
- added earth gravity fraction deconvolution option ("GRAV", which is basically "ACC"/9.81)
- added testing capabilities using `rs-test`
- added a script to translate seismic data to Raspberry Shake UDP packet format for testing
- changed warning and error message colors in terminal stdout
- alert module stdout STA/LTA messages now colorized
- added `rsudp.c_custom` as an independent thread to run custom code
- added and expanded explicit docstrings and comments, as well as Sphinx `conf.py` file
- turned off alert module STA/LTA live printed output when `settings['settings']['debug']` is `False` in order to keep systemd log file size down
- streamlined alert sound module operation; no longer writes temporary sound file on every alert when using `ffplay`
- added `rsudp.__version__` linked to version in `setup.py`

## changes in 0.4.2
- the station's [Flinn-Engdahl region](https://en.wikipedia.org/wiki/Flinn%E2%80%93Engdahl_regions) is added to tweets when the station inventory is available through FDSN
- changed message format to include live link to StationView
- added redundancy to tweets
- added Telegram bot module (https://core.telegram.org/bots)

## changes in 0.4.1
- fixed [#1](https://github.com/raspishake/rsudp/issues/1) which caused the writing module to crash on some machines
- added a module that posts to twitter in the event of an alert, and can also post screenshots saved by the plot module
- fixed [#2](https://github.com/raspishake/rsudp/issues/2) which caused an error at 28>1 30>1 and 31>1 month rollovers
    - should also address issues arising from leap second addition
- added Twitter bot module (https://developer.twitter.com/en/apps)
---
title: 'rsudp: A Python package for real-time seismic monitoring with Raspberry Shake instruments'
tags:
  - Python
  - seismology
  - monitoring
  - earthquake alerts
  - Raspberry Shake
authors:
  - name: Ian M. Nesbitt^[Corresponding author]
    orcid: 0000-0001-5828-6070
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Richard I. Boaz
    affiliation: 1
  - name: Justin Long
    affiliation: 3
affiliations:
  - name: Raspberry Shake, S.A.
    index: 1
  - name: School of Earth and Climate Sciences, University of Maine
    index: 2
  - name: Waterbear Energy
    index: 3
date: 11 June 2020
bibliography: paper.bib

---

# Statement of Need

The uses of low-cost seismographs in science and education are becoming more
widely known as these devices become more popular
[@Anthony:2018; @Walter:2019; @Diaz:2020; @Lecocq:2020; @Subedi:2020; @Winter:2021].
Raspberry Shake seismographs are commonly used in schools, by Shake community
members, and other individuals having no formal training in seismology. The
existence of this class of instruments highlighted the need for easy-to-use
visualization and notification software to complement these devices. Because
all Raspberry Shake instruments are able to forward data as user datagram
protocol (UDP) packets, taking the opportunity to exploit the existence of this
streaming data was obvious.

While the plotting may be the centerpiece of the program, perhaps the most
useful aspect of `rsudp` for researchers is its ability to monitor sudden motion
and trigger various actions when events are detected. This software's ability
to monitor data and trigger alerts with little processing overhead could be
critical to monitoring units in the field. Additionally, `rsudp` was designed for
extensibility, meaning that it leaves room for users to add their own code to
be run when events are detected. The demands of real-time seismic processing
require that calculations must be made quickly and remain stable for weeks or
months without user intervention. `rsudp` aims to achieve both of these things,
maintaining a codebase lean enough to run on Raspberry Pi but intuitive enough
that users can learn the theory of real time continuous data processing and
contribute code of their own. Programs that do similar tasks are usually not as
fully-featured, cost money, are unmaintained, are difficult to fork and
customize, or are complex to set up and run. We have tried to keep dependencies
to a minimum, the code base understandable, and installation simple across
multiple platforms.

Similar JAVA programs, including Swarm [@USGS:2020], jAmaSeis
([http://www.iris.edu/hq/jamaseis/](http://www.iris.edu/hq/jamaseis/)), and
SeisGram2K ([http://alomax.free.fr/seisgram/SeisGram2K.html](http://alomax.free.fr/seisgram/SeisGram2K.html))
have broader scope but less extensibility, and while they can all be set up to
run with the Rasbperry Shake, they can not read Raspberry Shake UDP format.
Therefore, accessing near-realtime data will necessarily use more bandwidth
and place processing load on the Shake itself. More powerful network processing
suites like Earthworm
([http://www.earthwormcentral.org/](http://www.earthwormcentral.org/)) are
difficult to set up and do not easily produce kiosk-ready live visualizations.
SeisComP4 ([https://www.seiscomp.de](https://www.seiscomp.de)), while arguably
the industry standard for network processing, requires a license for full
functionality, and is typically meant for high-level seismological
institutions.

![Chart of producer and consumer threads and the organization of data flow in `rsudp`. In order to maximize computational efficiency, features are broken into modules—each module constituting a thread—and data is passed to each module through an asynchronous queue. Inset: thread hierarchy and ownership chart, color-coded by function. Note that the Plot module is owned by the main thread, since `matplotlib` objects can only be created and destroyed by the main thread.\label{fig:flow}](flow.png)

# Summary

`rsudp` is a multi-featured, continuous monitoring tool for both Raspberry
Shake seismographs⁠, used to record both weak and strong ground motion⁠—and
Raspberry Boom pressure transducer instruments, used to record infrasound
waves. To encourage hands-on community involvement, `rsudp` is open-source,
written in Python, and utilizes easy-to-use tools common to the seismology
community, including `matplotlib` visualizations [@Hunter:2007] and the `obspy`
seismic framework for Python [@Beyreuther:2007; @Megies:2011; @Krischer:2015].
`rsudp` is multi-threaded and architected according to a modular
producer-consumer data-flow paradigm (\autoref{fig:flow}). The detection
algorithm employs a recursive short-term/long-term average ratio (STA/LTA)
computation threshold function from `obspy`, executed repeatedly within a loop
over the incoming data.

![An earthquake trace plotted with a spectrogram on multiple data channels in `rsudp`. The spectrograms are a representation of the fraction of maximum frequency power of the signal on each channel over the duration of the plot. Note that the first channel is data recorded with a geophone (EHZ), and the next three are accelerometers (ENE, ENN, ENZ).\label{fig:event}](event.png)

`rsudp` can be used by seismologists as a data analysis tool operating in real
time, and as a way for students, citizen scientists, and other end-users to
easily visualize and conceptualize live-streaming seismic data
(\autoref{fig:event}). Using the application’s simple and straightforward
framework, power-users can run their own custom code in the case of detected
strong motion. The distribution already contains many useful data-processing
modules, including: sound alerts, automated and instantaneous social media
notifications, data-forwarding, real-time seismic amplitude (RSAM) forwarding,
integrated logging, a miniSEED data archiver, and external script execution
(for example, to control input/output pins or some other programmable action).
The combination of speed, easy-to-interpret visualization, and ease of
customization makes `rsudp` a valuable and instructive companion to the
Raspberry Shake family of instruments for researchers, students, and amateur
seismologists alike.


# Acknowledgements

Financial support for this project comes from Raspberry Shake S.A. We are
grateful to Trinh Tuan Vu for his help authoring Windows setup scripts, Fabian
Walter and Calum Chamberlain for helpful reviews, and to Leif Lobinsky for
design input.


# References
:py:data:`rsudp.entry_points` (convenience fx)
=====================================================

.. versionadded:: 1.0.3

These are some convenience functions for editing rsudp's settings and
monitoring log output.

.. automodule:: rsudp.entry_points
    :members:


`Back to top ↑ <#top>`_
Contributing to this project
#####################################

.. |newissue| raw:: html

   <a href="https://github.com/raspishake/rsudp/issues/new" target="_blank">new issue</a>

.. |trigger_on_off| raw:: html

   <a href="https://docs.obspy.org/tutorial/code_snippets/trigger_tutorial.html#advanced-example" target="_blank">trigger on-off events</a>

Code contributions
*********************************

Contributions to this project are always welcome.
If you have questions or comments about how this software works,
we want to hear from you.
Even if coding isn't your thing,
we want to make it easier for you to get involved.
We monitor both our forums at https://community.raspberryshake.org, and our GitHub
issues page at https://github.com/raspishake/rsudp/issues.

See our resources at :ref:`add_your_own` to learn more about creating your own
consumer modules like the ones already used to plot,
send social media messages, play sounds, and more.

Since the Producer function passes an ``ALARM`` queue message when it sees
:py:class:`rsudp.c_alert.Alert` indicate an alarm state,
other modules can be easily added and programmed to do something when they
see this message as well.

The :py:class:`rsudp.c_custom.Custom` class makes running custom code easy.
If you have suggestions for feature addition of a new module, please open a
|newissue| with the "enhancement" tag.

If you're a developer or feeling adventurous,
here are some fun potential projects:

- Windows batch scripts similar to the provided UNIX ones
- GPIO pin interaction module (lights, motor control, buzzers, etc.)
- IFTTT integration
- Integration into other social media apps beyond Telegram and Twitter
- plot |trigger_on_off| for times set in ``ALARM`` and ``RESET`` messages using :py:func:`matplotlib.pyplot.axvline`::

    on_events = [UTCDateTime1, UTCDateTime3]
    for time in on_events:
        plt.axvline(x=time, color='r', linestyle=(0, (14,14)), linewidth=1, alpha=0.7)
    off_events = [UTCDateTime2, UTCDateTime4]
    for time in off_events:
        plt.axvline(x=time, color='g', linestyle=(0, (14,14)), linewidth=1, alpha=0.7)

- a more efficient plotting routine (I'm kidding, that's actually not a fun one)
- a way to run the plot module with the Agg backend in matplotlib, allowing for the creation of screenshots without the need for a plot window


Bugs
***********************

This software, like most, contains bugs and errors.
If you find a bug, please create a GitHub issue.
Be sure to describe the problem clearly, attach your logs
(:code:`/tmp/rsudp/rsudp.log`) and/or copy/paste command line output
in triple backticks \`\`\` like this \`\`\` to format it as code.

`Back to top ↑ <#top>`_
:py:data:`rsudp.raspberryshake` (main library)
=====================================================

.. note::
    If you are starting this program from a command line by using
    ``rs-client`` or a start script, the functions in this
    library will be executed by the :mod:`rsudp.client` automatically.
    See :ref:`running` for details.

    If you are a developer looking for helper functions to create a new module,
    you have come to the right place.

This is the main library powering rsudp.
It contains common functions to open a port listener,
get data packets, parse the information in those packets,
and help consumer modules update their data streams.

Prior to working with data from the port,
this library must be initialized using :func:`rsudp.raspberryshake.initRSlib`:

.. code-block:: python

    >>> import rsudp.raspberryshake as rs
    >>> rs.initRSlib(dport=8888, rsstn='R940D')

.. note:: This request will time out if there is no data being sent to the port.

After initializing the library, the :func:`rsudp.raspberryshake.getDATA`
function and its derivatives will be available for use.


.. automodule:: rsudp.raspberryshake
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
.. _install:

Installing rsudp
#####################################

.. role:: bash(code)
   :language: bash

.. |installation_video| raw:: html

   <a href="https://youtu.be/e-kyg55GZyA" target="_blank">installation tutorial video</a>

Installation is covered in our |installation_video|.


Installing on Linux & MacOS
*********************************

A UNIX installer script is available at :bash:`unix-install-rsudp.sh`.
This script checks whether or not you have Anaconda or Miniconda installed,
then downloads and installs it if need be.
This script has been tested on both :bash:`x86_64` and :bash:`armv7l`
architectures (meaning that it can run on your home computer or a Raspberry Pi)
and will download the appropriate Anaconda distribution, set up a virtual Python environment,
and leave you ready to run the program. To install using this method, open a Terminal and
enter the following command:

.. code-block:: bash

    bash unix-install-rsudp.sh

.. warning::
    In order for this installer to work correctly,
    your Anaconda/Miniconda base directory must be in your home folder,
    for example: :bash:`~/anaconda3` or :bash:`~/miniconda3`.
    If this is not the case, you could end up with a second conda installation overriding your first.

.. note::
    **The installer script will pause partway through to ask if you would like to make the**
    :bash:`conda` **command executable by default.**
    **This is done by appending the line below to your** :bash:`~/.bashrc` **file.**
    This is generally harmless, but if you have a specific objection to it,
    hitting any key other than "y" will cause the script to skip this step.
    You will have to manually run the :bash:`conda` executable in this case, however.
    If you choose to do it manually later, follow the instructions in the section below.

You are now ready to proceed to the next section, :ref:`settings`.


.. _source:

Run the sourcing line
-----------------------------------------------------------------

If you are running UNIX and having trouble running ``conda activate``
commands, your operating system probably doesn't know where to look
for Anaconda. To fix this, we need to tell it to read Anaconda's
definitions telling it where to look.

If you are running an x86 (desktop or laptop) or AMD type processor,
run the following command:

.. code-block:: bash

    . $HOME/miniconda3/etc/profile.d/conda.sh

or on ARMv7 (Raspberry Pi) architecture with Raspbian OS:

.. code-block:: bash

    . $HOME/berryconda3/etc/profile.d/conda.sh

where :bash:`$HOME` is the home directory of the current user.

.. note::

    You can run :bash:`uname -m` to check your computer's architecture.


Add the sourcing line to your :py:data:`~/.bashrc`
-----------------------------------------------------------------

The UNIX installer script *should* do this step automatically,
but if you have this problem consistently, you may need to add this
line to your ``~/.bashrc`` file.
The following step will append the sourcing line to
the end of your :bash:`~/.bashrc` is the following (architecture-dependent):

On x86/AMD systems:

.. code-block:: bash

    echo ". $HOME/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc

or on ARMv7:

.. code-block:: bash

    echo ". $HOME/berryconda3/etc/profile.d/conda.sh" >> ~/.bashrc


Updating
---------------------------------

Unix users can update the repository to the latest development version by running the following commands:

.. code-block:: bash

    cd /rsudp/location
    git pull
    bash unix-install-rsudp.sh

The update script will replace the previous default settings file
(:bash:`~/.config/rsudp/rsudp_settings.json`) with a new settings file.
If you use the default settings file, you will need to copy some old values over to the new file.
The reason for this is that the default settings file may change (i.e. add or modify sections of values)
and thus must be rewritten when updating. On Linux, backed up settings files will be named
:bash:`~/.config/rsudp/rsudp_settings.json.~x~`, where :bash:`x` is an integer.
On Mac, the backed up file will simply be named :bash:`~/.config/rsudp/rsudp_settings.json~`.
To back up the settings file yourself to a location that will not be overwritten,
you can do a command similar to the following:

.. code-block:: bash

    cp ~/.config/rsudp/rsudp_settings.json ~/.config/rsudp/rsudp_settings.json.bak


Installing on Windows
*********************************

The Easy Way
---------------------------------

You can follow these steps to both install and update rsudp.

.. |github_download| raw:: html

   <a href="https://github.com/raspishake/rsudp/releases/latest/download/rsudp.zip" target="_blank">Download</a>

.. |github_latest| raw:: html

   <a href="https://github.com/raspishake/rsudp/releases/latest/" target="_blank">latest release</a>

1. |github_download| and unzip the software from the |github_latest| in the GitHub repository.
2. Double click the file named ``win-install-rsudp.bat`` in the unzipped folder. You may need administrator privileges for this step.

The install will take several minutes. When it is done, you will have a new settings file at
``~/.config/rsudp/rsudp_settings.json``. Edit this file to change how rsudp runs.

For explanations of the various settings fields and values, head to :ref:`settings`.


Advanced Users
---------------------------------

.. |miniconda3| raw:: html

   <a href="https://docs.conda.io/en/latest/miniconda.html" target="_blank">Miniconda3</a>


1. Download and install Anaconda3 or |miniconda3|.
2. Open an Anaconda Prompt.
3. Execute the following lines of code:

.. code-block:: bash

    conda config --append channels conda-forge
    conda create -n rsudp python=3 matplotlib=3.1.1 numpy=1.16.4 future scipy lxml sqlalchemy cryptography obspy
    conda activate rsudp
    pip install rsudp

.. |windows_tutorial| raw:: html

   <a href="https://windowsloop.com/install-ffmpeg-windows-10/" target="_blank">this tutorial</a>

If you wish to play sounds on Windows, please follow steps 1-8 in |windows_tutorial|
in order to install :code:`ffmpeg` and add it to your system's path variable.


You are now ready to proceed to the next section, :ref:`settings`.


`Back to top ↑ <#top>`_
.. _youtube:

YouTube walkthroughs
########################################

.. aspect ratio from YouTube embed code is 560x400, which is 16x9
    this also equates to 400x225 or 800x450

1. Installation

    .. raw:: html

        <iframe width="560" height="400" src="https://www.youtube-nocookie.com/embed/e-kyg55GZyA" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

2. Adjust settings and run

    .. raw:: html

        <iframe width="560" height="400" src="https://www.youtube-nocookie.com/embed/HA9k3CzmgLI" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

3. Set up YouTube streaming (coming soon)


Live demo
#########################################

Live streaming Raspberry Shake & Boom data from our office in Panama (when available).

    .. raw:: html

        <iframe width="560" height="400" src="https://www.youtube-nocookie.com/embed/live_stream?channel=UCqxETqiBOCMH7fy-d0XgIWw" frameborder="0" allowfullscreen></iframe>


Detected earthquake
#########################################

.. |earthquake| raw:: html

   <a href="https://www.emsc-csem.org/Earthquake/earthquake.php?id=806235" target="_blank">an earthquake</a>

rsudp detects |earthquake|!

    .. raw:: html

        <iframe width="560" height="400" src="https://www.youtube-nocookie.com/embed/pT_PkKKxFeM" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
:py:data:`rsudp.c_rsam` (RSAM analysis)
=====================================================

.. automodule:: rsudp.c_rsam
    :members:

.. note: This module is courtesy of Shake community member Justin Long (@crockpotveggies).

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
:py:data:`rsudp.packetloss` (track packets)
=====================================================

.. automodule:: rsudp.packetloss
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
:py:data:`rsudp.t_testdata` (test pre-producer)
=====================================================

.. automodule:: rsudp.t_testdata
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
:py:data:`rsudp.c_consumer` (master consumer)
=====================================================

.. automodule:: rsudp.c_consumer
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
.. _flow:

Program architecture and flow
#################################################

rsudp is laid out in a way that is standard for many continuous monitoring
softwares, but which may be unfamiliar to some developers.

The program's organization relies on two hierarchies: thread hierarchy and
data hierarchy. This allows the program to distribute data and messages
between different modules efficiently while maintaining programmatic
integrity required by rsudp's dependencies (:py:mod:`matplotlib.pyplot`
objects, for example, can only be owned by the master thread).

.. _flow_diagram:
.. figure::  _static/rsudp-flow.png
    :align:   center

    Data and (inset) thread hierarchies for rsudp.
    This figure shows a flow chart of data as it makes its way through
    rsudp's producer-consumer architecture,
    and eventually to the end-consumer processing the data.
    Thread hierarchies are shown in the inset chart.
    If you are looking to add a module to rsudp,
    it will be a worker class (red) and will receive data from
    the master consumer like the rest of the workers.


Thread layout
*************************************************

The program relies on the :py:mod:`threading` and :py:mod:`queue` modules
to spin up multiple threads which all receive and work with data
independently of one another.

First, the :py:mod:`rsudp.client` (the "main" or "parent" thread) gathers
and parses settings. The client then instantiates the relevant
:py:class:`threading.Thread` objects and passes settings and queues to
them. Next, it starts each of these consumer threads as "child" threads,
then finally the producer, also as a child thread.

If the :py:class:`rsudp.c_plot.Plot` thread is enabled, the program will
start that last, since plotting must be run as a loop controlled by the
main thread (:py:mod:`rsudp.client`) which can only be done once it has
started all of its child threads.


.. _producer-consumer:

Producer-consumer message passing
*************************************************

Data is read off the port by the data producer
(:py:class:`rsudp.p_producer.Producer`) and passed directly to the
(:py:class:`rsudp.c_consumer.Consumer`) via a first-in first-out (FIFO)
queue object.

This master consumer then duplicates these messages and
passes them to each sub-consumer queue destination in a list. At the
other end of each of these queues is a sub-consumer module, denoted
with a :py:data:`c_` before its module name.

Sub-consumers read messages from their queues and process data in
their logic loops. Some build :py:class:`obspy.core.stream.Stream` with
the data passed to them, while some ignore the data and watch for
status messages.

.. _message-types:

Message types
=================================================

.. versionadded:: 0.4.3 the :code:`RESET` message type was added.

Currently, the message types are as follows.

========= ==========================================
 Message              Format example
========= ==========================================
 data      ``b"{'EHZ', 1582315130.292, 14168, 14927, 16112, 17537, 18052, 17477, 15418, 13716, 15604, 17825, 19637, 20985, 17325, 10439, 11510, 17678, 20027, 20207, 18481, 15916, 13836, 13073, 14462, 17628, 19388}"``
 ALARM     ``b'ALARM 2020-02-23T06:56:40.598944Z'``
 RESET     ``b'RESET 2020-02-23T06:56:55.435214Z'``
 IMGPATH   ``b'IMGPATH 2020-02-23T06:59:19.211704Z /home/pi/rsudp/screenshots/R24FA-2020-02-23-065919.png'``
 TERM      ``b'TERM'``
========= ==========================================

.. note::

    The above message formats are Python bytes objects, not traditional
    strings. The difference between a bytes object and a string is
    outlined briefly
    `here <https://www.geeksforgeeks.org/byte-objects-vs-string-python/>`_.

**ALARM** messages are sent by :py:class:`rsudp.p_producer.Producer`
when it sees the :py:data:`rsudp.c_consumer.Alert.alarm` flag set to
``True``. This can trigger all sorts of actions. For example, when the
:py:class:`rsudp.c_alertsound.AlertSound` module is enabled and sees
this message, it uses ffmpeg or libav to play a sound. The social media
classes :py:class:`rsudp.c_tweet.Tweeter` and
:py:class:`rsudp.c_telegram.Telegrammer` both use this message to
instantly broadcast to their respective platforms.

**RESET** messages are sent by :py:class:`rsudp.p_producer.Producer`
when it sees the :py:data:`rsudp.c_consumer.Alert.alarm` flag set to
``True``. Similar to ALARM messages, consumers can be programmed for
an essentially infinite number of things upon seeing this message.

**IMGPATH** messages are placed on the master queue by the
:py:func:`rsudp.c_plot.Plot.savefig` function, if and when a screenshot
figure is saved to disk. This is currently only used by the social media
modules, :py:class:`rsudp.c_tweet.Tweeter` and
:py:class:`rsudp.c_telegram.Telegrammer` which then send the saved image
to their respective social media platforms' APIs for broadcast.

**TERM** messages are the universal signal for rsudp to quit.
They generally start at the Producer and are passed through the
data hierarchy as normal data would.


.. _add_your_own:

Adding your own consumer modules
*************************************************

Adding consumer modules is easy in theory, when you understand the
workings of rsudp's layout. Using the existing modules' code architecture
is likely useful and should be encouraged, so feel free to follow along
with what we have already laid out in the code base.

There are three main things that need to happen in order to add a consumer.

1. Create a new module, named ``c_mymodule.py``, in the ``rsudp`` directory.
2. Add a section to your settings file which will tell rsudp what settings to pass to your module.
3. Add code to :py:func:`rsudp.client.run` to pass settings and a queue to your module, and start it.

And some optional things to do in case you plan on submitting a pull request:

4. Add documentation in the form of reStructuredText-formatted docstrings (see examples below)
5. Add testing capability to your module.


Sample consumer
=================================================

Below is a sample consumer construction to modify for your own purposes.
It receives all queue messages (outlined in :ref:`producer-consumer`)
and can be made to do pretty much whatever you wish,
until it receives a ``TERM`` queue message.

This consumer is created from a
:py:class:`rsudp.raspberryshake.ConsumerThread` object,
which in turn modifies the :py:class:`threading.Thread` class.

.. code-block:: python

    import sys
    from rsudp.raspberryshake import ConsumerThread
    from rsudp import printM

    class MyModule(ConsumerThread): # this means MyModule will be based on the ConsumerThread class
        '''
        Documentation of your new module class goes here.
        Below is the format of two types of *param* string, which tell the
        documentation parser to inform users that this object needs the user to
        pass it a queue in order to work correctly.

        The first string, for the ``q`` parameter, has the type as the
        middle object and the caption after. The second one, ``thing1``
        could either be a string or a boolean value,
        so we move the type for it to its own row with the types listed after.
        Sphinx, the documentation formatter, will be able to combine these into
        one object describing the parameter.

        :param queue.Queue q: queue of data and messages sent by :class:`rsudp.c_consumer.Consumer`
        :param thing1: a passed parameter that's either a string or a boolean (True/False)
        :type thing1: bool or str
        '''
        def __init__(self, q, thing1    # ... probably some more parameters to pass to the class
                    )
            super().__init__()
            self.sender = 'MyModule'
            self.alive = True
            self.queue = q
            self.thing1 = thing1
            # ... lots of other stuff to initialize your module
            printM(self.thing1, sender=self.sender)

        def getq(self):
            '''
            Reads data from the queue and returns the queue object.

            Since this function returns something, it has return
            strings (*rtype* stands for return type) so that the
            user reading the documentation knows what they'll get
            back if they call it.

            :rtype: bytes
            :return: The queue object.
            '''
            d = self.queue.get()
            self.queue.task_done()
            return d

        def run(self):
            '''
            Documenting how my cool module runs!

            Right now, its only function is to get and read queue messages
            to see if one of them has the ``TERM`` message in it,
            at which point it quits.
            '''
            printM('Starting.', sender=self.sender)
            # some stuff to execute here at runtime before looping

            while self.alive:
                # main loop, do something until self.alive == False
                d = self.getq()
                if 'TERM' in str(d):
                    self.alive = False

            # now exit
            printM('Exiting.', sender=self.sender)
            sys.exit()


Adding your module to the settings file
=================================================

An example settings section is given here.
As a reminder, each settings section except the last one
is required to have a comma after its closing brace to conform
to JSON standards.
Here let's assume this is not the last JSON section,
so we include the comma:

.. code-block::

    "mymodule": {
        "enabled": true,
        "thing1": "first thing"},


Adding your module to ``client.py``
=================================================

Since all modules are started from the client's :py:func:`rsudp.client.run`
function, you will need to add a section of code to the client to tell it
how to start your module.
An example based on the JSON section above is given here.

.. code-block:: python

    from c_mymodule import MyModule

    # ... lots of other stuff in client.py

    def run(settings, debug):

        # ... setting up other modules

        if settings['mymodule']['enabled']:
            # first, gather settings
            thing1 = settings['mymodule']['thing1']
            # then, set up queue
            q = mk_q()
            # then, start a MyModule instance with the settings you got earlier
            mymod = MyModule(q=q, thing1=thing1)
            # now, pass this instance to the process list to be started below
            mk_p(mymod)

        # ...

        # this part already exists, but just to show you where in sequence your code should be:
        start()

        # ...


.. _add_testing:

Testing your module
=================================================

Formal testing of new modules is easy in rsudp.

The :py:func:`rsudp.client.test` function is set to run any enabled
module by default. If the module is not enabled in the default
settings, you can add a line to the
:py:func:`rsudp.test.make_test_settings` that specifies

.. code-block:: python

    settings['your_module']['enabled'] = True

The second step is to add your test to the dictionary of tests in
:py:mod:`rsudp.test`, so that it gets reported. For example:

.. code-block:: python

    TEST = {
            # other tests
            # ...
            'c_mytest':             ['something I am testing for  ', False],
            'c_anotherone':         ['some other thing I test     ', False],
    }

Each dictionary item is constructed as a two-item list,
where the first item is the description string,
and the second is the status of the test
(False is failure and True is passing).

Then, in your module, you can import the test dictionary and modify
the status of your tests like so:

.. code-block:: python

    from rsudp.raspberryshake import ConsumerThread
    from rsudp.test import TEST

    class MyModule(ConsumerThread):
        def __init__(self, q    # ...
                    )
            super().__init__()
            # ... stuff to initialize your module
            if abc:
                # this test occurs during initialization
                TEST['c_mytest'][1] = True

        def run(self):
            # some stuff here also
            if xyz:
                # this test is done at runtime
                TEST['c_anotherone'][1] = True
            while self.alive:
                # main loop, do something until self.alive == False
                # or you receive the TERM message
            # now exit
            printM('Exiting.', self.sender)
            sys.exit()


Suggesting features
*************************************************

As with other issues, if you have an idea for a feature addition but have
questions about how to implement it, we encourage you to post to our
forums at https://community.raspberryshake.org.

Thanks for supporting open source!


`Back to top ↑ <#top>`_
:py:data:`rsudp.packetize` (data formatter)
=====================================================

The ``packetize`` module is a utility that turns miniSEED data into text files
containing the Raspberry Shake UDP data format (see :ref:`producer-consumer`).
It can be run either from another python module using the
:py:func:`rsudp.packetize.packetize` function, or from the command line.

Python usage:

.. code-block:: python

    from rsudp.packetize import packetize
    packetize('test.mseed', 'output.txt')

Command line usage:

.. code-block:: bash

    conda activate rsudp
    python packetize.py -i test.mseed -o output.txt

.. note::

    Command line usage must be done from within an environment in which
    ``obspy`` is installed as a python3 package.

.. automodule:: rsudp.packetize
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
.. _troubleshooting:

Troubleshooting rsudp
#################################################

In general, troubleshooting should be fairly straightforward.
Most output goes either to the command line or to ``/tmp/rsudp/rsudp.log``.
If you have a recurring issue, please see if the logs give any indication
as to what is going on, and check relevant sections below.

By far the most common problem you'll run into if you are a first-time user
is the inability to see data coming in.

.. _remote:

Remote (rsudp) side
*************************************************

Can't run the ``conda`` command
=================================================

See :ref:`source`.

No data received
=================================================

If you are getting an ``IOError('No data received')``, most likely one of three
things is wrong:

#. Your Shake is not forwarding to the correct address and port — see :ref:`local`
#. rsudp is configured to listen on the wrong port
#. There is a firewall between your computer and the Shake — see :ref:`middle`

In an error like this, your rsudp output will display some helpful info just
above the error text. Scroll up to where you see something like this:

.. code-block:: bash

    2020-02-22 21:49:01 [Init] ERROR: No data received in 10 seconds; aborting.
    2020-02-22 21:49:01 [Init]        Check that the Shake is forwarding data to:
    2020-02-22 21:49:01 [Init]        IP address: 192.168.1.118    Port: 8888
    2020-02-22 21:49:01 [Init]        and that no firewall exists between the Shake and this computer.

Make note of the values after **IP address** and **Port**.

First, check to make sure the address and port have been configured on
the Shake side to cast data to the address ``192.168.1.118`` and port ``8888``
— see :ref:`local`.

If you are sending data to a different port than 8888, check the ``"port"``
value in your settings file (``~/.config/rsudp/rsudp_settings.json``) reflects
the port to which data is being sent.

Finally, if you are for example sending data from a location outside your home
network to somewhere inside your home network, your router may not be letting
data come through the port you specified, and you may need to specify rsudp to
send data to the router's public (externally facing) IP address.

Data stops flowing or is inconsistent
=================================================

This may be due to network problems. rsudp is designed to be able to ingest all
data sent to it by the Shake. However, since the Shake uses UDP, which is not a
guaranteed delivery protocol, some packets may be dropped.

There are several reasons why this might happen.

#. You are using WiFi and there is an unstable connection
#. The router nodes between the Shake and your computer may restart or be overloaded
#. The Shake may be a great distance from your computer (across the globe)
#. The Shake or your computer may have a slow connection

Typically, if the Shake is connected via Ethernet, and sending to a computer that
also uses Ethernet, you will experience on average zero dropped packets in a given
24 hour period. However, if you are sending data across an unstable connection,
you could experience 40% or more dropped packets.

To monitor packet loss over time, you can run our :py:mod:`rsudp.packetloss` script.
For example, to report dropped packets on port 8888 in periods of 1 hour at a time:

.. code-block:: bash

    conda activate rsudp
    rs-packetloss -p 8888 -f 3600

where ``-p 8888`` specifies port 8888 and ``-f 3600`` specifies 3600 seconds between
reports.

This will run indefinitely until the CTRL+C keys are pressed.


.. _local:

Shake (local) side
*************************************************

To set up or change Datacasting (also known as data forwarding or UDP forwarding)
navigate to your Shake's home page, then click on ``Settings > DATACAST``.

In the above example, you should configure your Shake to send data to:

================= ================
Label              Value
================= ================
Target Host IP     192.168.1.118
Target Port        8888
================= ================

Then press the blue plus button on the right side of the row.

.. _middle:

Middle (firewalls)
*************************************************

Home network
=================================================

Almost every home router in existence has a firewall between the outside of the
network it resides on and the "inside", i.e. the local in-home network it is
responsible for. (If you're working on a :ref:`school-net`, this works slightly
differently)

Most home routers also have a feature called "Port Forwarding" which will forward
data through the firewall from an external port to an internal port at a specific
IP address.

In rsudp's case: if we assume your Shake is somewhere else (i.e. not on your home
network) then it will be forwarding data to the external side of your router, and
you will need to tell your router to let that data through and where to send it.

First of all, you will need to know your router's IP address. There are many
online services that will do this. One of the safer ways to figure it out is just
`searching "what is my IP" on DuckDuckGo
<https://duckduckgo.com/?q=what+is+my+IP&t=canonical&ia=answer>`_
(DuckDuckGo will not store your information, while many other sites will).
Your IP should appear right under the search bar.

Let's say DuckDuckGo tells you that your IP address is ``28.14.122.178``.

Let's look at the following configuration:

============== ================ ======================
Device          IP               Public or Private IP
============== ================ ======================
Your Shake      130.112.21.12    Public
Your router     28.14.122.178    Public (external)
Your router     192.168.1.1      Private (internal)
Your computer   192.168.1.118    Private
============== ================ ======================

In this case, you must configure your Shake to forward UDP data to address
``28.14.122.178`` at, for example, port ``8888`` (i.e. port 8888 on the external side
of your router). Then, configure your router to forward data on external UDP port
``8888`` to internal address ``192.168.1.118`` and port ``8888``.

You should then be able to receive data on your computer.

.. note::

    Some internet service providers (ISPs) do not let you change your router's
    settings yourself. In this case, you will need to call them and ask them to
    configure port forwarding for external port ``8888`` to forward data to the same
    port at the internal IP address ``192.168.1.118``.

.. _school-net:

School or university network
=================================================

If you are on a school or university network, often security is much more strict.
In your home network, data is usually free to move around internally on the
network. On school networks, individual devices are usually not allowed to talk
much to each other. So even if your Shake is on the internal network, you may
still need to notify the school's IT team to give your Shake permission to send
data to another computer on the network.

They may be able to help with configuration of the setup as well, although they
usually have difficult jobs, so don't be too hard on them!


Other issues
*************************************************

For more information about Raspberry Shake UDP output, see the |rs-manual-udp| of the
|rs-manual|.

If you are having a technical support issue other than one described above,
please post the issue you are having to our forum at
https://community.raspberryshake.org. We would be glad to help you solve your
issue there.

.. |newissue| raw:: html

   <a href="https://github.com/raspishake/rsudp/issues/new" target="_blank">new issue</a>

.. |rs-manual-udp| raw:: html

   <a href="https://manual.raspberryshake.org/udp.html#udp" target="_blank">UDP section</a>

.. |rs-manual| raw:: html

   <a href="https://manual.raspberryshake.org/" target="_blank">Raspberry Shake User Manual</a>

If it turns out that we cannot solve it without a bug fix in the code, please
submit a |newissue| on GitHub.
Be sure to describe the problem clearly, attach your logs
(:code:`/tmp/rsudp/rsudp.log`) and/or copy/paste command line output
in triple backticks \`\`\` like this \`\`\` to format it as code.

Our small team thanks you for your patience and cooperation!


`Back to top ↑ <#top>`_
:py:data:`rsudp.p_producer` (data producer)
=====================================================

.. automodule:: rsudp.p_producer
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
:py:data:`rsudp.client` (central module)
=====================================================

.. automodule:: rsudp.client
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
:py:data:`rsudp.c_plot` (plot data)
=====================================================

.. _event_plot:
.. figure::  _static/event.png
    :align:   center

    An event recorded on a vertical geophone by rsudp and saved by the plotting module.


.. automodule:: rsudp.c_plot
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
:py:data:`rsudp.c_telegram` (Telegram alerts)
=====================================================

.. automodule:: rsudp.c_telegram
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
:py:data:`rsudp.c_testing` (test consumer)
=====================================================

.. automodule:: rsudp.c_testing
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
:py:data:`rsudp.c_write` (write data)
=====================================================

.. automodule:: rsudp.c_write
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
:py:data:`rsudp.c_forward` (forward data)
=====================================================

.. automodule:: rsudp.c_forward
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
.. _running:

Running rsudp
#################################################

After modifying the settings file to your liking, you are ready to start rsudp.

.. role:: bash(code)
   :language: bash

Starting rsudp on Unix
*************************************************

Unix users may prefer to start with the easy-to-use start script available in the git repository:

.. code-block:: bash

    bash unix-start-rsudp.sh


Starting rsudp on Windows
*************************************************

After modifying the settings file, Windows users can double click
on the batch file named ``win-start-rsudp.bat`` to start.


.. _running-manually:

Starting manually on Windows or Unix
*************************************************

.. |start_run_tutorial| raw:: html

   <a href="https://youtu.be/HA9k3CzmgLI" target="_blank">rsudp start/run tutorial video</a>

This start method is covered in our |start_run_tutorial|.

1. First, to activate the conda environment, type :bash:`conda activate rsudp`.
If you can't do this, and you're on Unix, you may need to :ref:`source`.

2. Next, configure a datacast stream (formerly known as a UDP stream)
to forward data to an open port on the computer where this program is running.
By default this port is :code:`8888`.

3. The UNIX installer will create a settings file in :bash:`$HOME/.config/rsudp/rsudp_settings.json`.
Change the settings in this file to control how the client operates.

.. note::

    Windows users will need to type :bash:`rs-client -d default` to dump the settings to a file
    the first time they run this program.

.. note::

    To dump the default settings to a different location of your choosing, type
    :bash:`rs-client -d /path/to/settings.json`.

.. note::

    As stated above, to rebuild and overwrite the default settings file in
    :bash:`$HOME/.config/rsudp/rsudp_settings.json`, type :bash:`rs-client -d default`

4. After modifying the settings file to your liking,
type :bash:`rs-client` to use the settings file at :bash:`$HOME/.config/rsudp/rsudp_settings.json`,
or :bash:`rs-client -s /path/to/settings.json` to run with a settings file other than the default one.

.. note::

    This library can only handle incoming data from one Shake per port.
    If for some reason more than one Shake is sending data to the port,
    the software will only process data coming from the IP of the first Shake it sees sending data.
    All data coming from any other Shake(s) will be ignored.

.. _run-test:

Running in demonstration/testing mode
*************************************************

See more about this functionality in :ref:`test`.

To start, open a Terminal or Anaconda Prompt window.

1. Activate the conda environment by typing :bash:`conda activate rsudp`.
2. Type :bash:`rs-test` and press enter.

Test data will begin flowing through the program.
Several features will be tested, including the
earthquake detection functionality, the alarm sound,
and the plot.


Quitting
*************************************************

You can force-stop rsudp in all operating systems by either closing the plot window (if open)
or by pressing :kbd:`Ctrl`\ +\ :kbd:`C` with the terminal window in focus.


`Back to top ↑ <#top>`_
:py:data:`rsudp.__init__` (initialization functions)
=====================================================

These are the initialization functions in rsudp.
The useful things here are likely ``printM``, ``printW`` and ``printE``
which interface with the logging utility to print messages, warnings,
and errors in color and save them to the logs.

.. automodule:: rsudp.__init__
    :members:


`Back to top ↑ <#top>`_
:py:data:`rsudp.c_alertsound` (play sound)
=====================================================

.. automodule:: rsudp.c_alertsound
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
`rsudp` |version|
#####################################
.. image:: https://img.shields.io/github/stars/raspishake/rsudp?style=social
    :target: https://github.com/raspishake/rsudp

| |github|
| **Continuous sudden motion and visual monitoring of Raspberry Shake data**
| *by Ian M. Nesbitt and Richard I. Boaz*

.. |github| raw:: html

   <a href="https://github.com/raspishake/rsudp" target="_blank">https://github.com/raspishake/rsudp</a>

.. |raspberryshake| raw:: html

   <a href="https://raspberryshake.org" target="_blank">Raspberry Shake</a>

Welcome to rsudp's documentation.
This program was written to parse and process live UDP data streams from
|raspberryshake| personal seismographs and
Raspberry Boom pressure transducer instruments.
rsudp allows users the options to see their data in real time, create alert parameters,
and be notified in various ways when their instrument detects sudden motion.
It is written in Python and is therefore highly customizable.

| In order to get a feel for what rsudp can do, check out our :ref:`youtube` page.
| If you prefer to read in-depth documentation, follow our written :ref:`tutorial`.
| Or, if you know what you're looking for, find it in :ref:`modules`.



.. _tutorial:

.. toctree::
    :numbered:
    :maxdepth: 2
    :caption: Tutorial guide

    about
    installing
    settings
    running
    daemon
    troubleshooting
    testing

.. toctree::
    :caption: Video resources

    youtube

.. toctree::
    :numbered:
    :maxdepth: 2
    :caption: Developers' guide

    theory
    contributing


.. _modules:

Code documentation
========================

The modules available in rsudp are organized by type below.

------------

.. toctree::
    :maxdepth: 2
    :caption: Library

    init
    raspberryshake
    helpers
    entry_points

.. toctree::
    :maxdepth: 2
    :caption: Clients

    client
    packetloss

.. toctree::
    :maxdepth: 2
    :caption: Producer and Consumer

    p_producer
    c_consumer

.. toctree::
    :maxdepth: 2
    :caption: Sub-Consumers

    c_alert
    c_rsam
    c_alertsound
    c_plot
    c_tweet
    c_telegram
    c_forward
    c_write
    c_custom

.. toctree::
    :maxdepth: 2
    :caption: Testing modules

    t_testdata
    c_testing
    test
    packetize



Function index
==================

Need to look something up?

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
:py:data:`rsudp.test` (test helpers)
=====================================================

.. versionadded:: 0.4.3

The test module. Here are the currently available tests, descriptions,
and their initial state:

.. code-block:: python

    TEST = {
        # permissions
        'p_log_dir':            ['log directory               ', False],
        'p_log_std':            ['stdout logging              ', False],
        'p_log_file':           ['logging to file             ', False],
        'p_output_dirs':        ['output directory structure  ', False],
        'p_screenshot_dir':     ['screenshot directory        ', False],
        'p_data_dir':           ['data directory              ', False],

        # network
        'n_port':               ['port                        ', False],
        'n_internet':           ['internet                    ', False],
        'n_inventory':          ['inventory (RS FDSN server)  ', False],

        # core
        'x_packetize':          ['packetizing data            ', False],
        'x_send':               ['sending data                ', False],
        'x_data':               ['receiving data              ', False],
        'x_masterqueue':        ['master queue                ', False],
        'x_processing':         ['processing data             ', False],
        'x_ALARM':              ['ALARM message               ', False],
        'x_RESET':              ['RESET message               ', False],
        'x_IMGPATH':            ['IMGPATH message             ', False],
        'x_TERM':               ['TERM message                ', False],

        # dependencies
        'd_pydub':              ['pydub dependencies          ', False],
        'd_matplotlib':         ['matplotlib backend          ', False],

        # consumers
        'c_plot':               ['plot                        ', False],
        'c_write':              ['miniSEED write              ', False],
        'c_miniseed':           ['miniSEED data               ', False],
        'c_print':              ['print data                  ', False],
        'c_alerton':            ['alert trigger on            ', False],
        'c_alertoff':           ['alert trigger off           ', False],
        'c_play':               ['play sound                  ', False],
        'c_img':                ['screenshot exists           ', False],
        'c_tweet':              ['Twitter text message        ', False],
        'c_tweetimg':           ['Twitter image message       ', False],
        'c_telegram':           ['Telegram text message       ', False],
        'c_telegramimg':        ['Telegram image              ', False],
        'c_forward':            ['forwarding                  ', False],
        'c_rsam':               ['RSAM transmission           ', False],
        'c_custom':             ['custom code execution       ', False],
    }

.. note::

    If you wish to add your own consumer module, the easiest way to test
    its functionality is to follow the instructions in
    :ref:`add_testing`, then add the relevant test to this dictionary.
    Then, you would import the :py:data:`rsudp.test.TEST` variable
    and modify the test result (:py:data:`TEST['your_test'][1] = True`)
    if the test passed.

    If your module is set not to start by default, and you are using the
    default settings file for testing, you will need to set
    ``settings['your_module']['enabled'] = True`` in
    :py:func:`rsudp.test.make_test_settings` prior to running the tests.




.. automodule:: rsudp.test
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
About rsudp
#####################################

.. |raspberryshake| raw:: html

   <a href="https://raspberryshake.org" target="_blank">Raspberry Shake S.A.</a>

.. |github| raw:: html

   <a href="https://github.com/raspishake/rsudp" target="_blank">here</a>

rsudp is an open source program developed for the community by |raspberryshake|
to actively monitor and plot UDP data cast output from Raspberry Shake instruments.

Broad overview
*************************************

rsudp is a collection of python monitoring and notification tools
to accompany Raspberry Shake seismographs and pressure transducers.
It contains code for continuous plotting, continuous monitoring of sudden motion,
and distribution of alert notifications via multiple media including audible sound
and the Telegram and Twitter social media platforms.

rsudp is able to run as a stable headless daemon or GUI application continuously
on most systems with sufficent RAM, and can be used to power both a display kiosk
and notification system simultaneously.

Why it's special
*************************************

rsudp's uses as an educational tool in both seismology and computer science are broad.

The demands of real-time seismic processing
require that calculations must be made quickly and
remain stable for weeks or months without user intervention.
rsudp aims to achieve both of these things,
while maintaining a codebase lean enough to run on Raspberry Pi
but intuitive enough that users can learn the theory of
real time continuous data processing and contribute code of their own.
The project's source repository is |github|.

Programs that do similar tasks are usually not as fully-featured, cost money,
are unmaintained, or are complex to set up and run.
We have tried to keep dependencies to a minimum and installation simple
across multiple platforms.

In addition, we spent many hours making rsudp's plotting routines a beautiful
and informative way to explore the vibrations that move through Earth's crust
(see 4d-event_ figure below).
Whether looking at an earthquake trace from far away or a car going by on a street,
the plots are designed to show the user the character of the vibration in an easy
to grasp yet informative format.


.. _4d-event:
.. figure::  _static/4d-event.png
    :align:   center

    An earthquake event recorded on the accelerometer and geophone channels of a
    Raspberry Shake 4D, recorded by rsudp and saved by its plotting module.


While the plotting may be the centerpiece of the program,
perhaps the most useful aspect of rsudp is its ability to monitor sudden motion
and trigger multiple outcomes and actions when events are detected.
Additionally, the way it was designed leaves room for developers
to add their own code to be run when events are detected.

rsudp is a special and unique piece of software in the seismological community
which brings easy-to-use open-source monitoring to low cost instrumentation.

What can't rsudp do?
*************************************

.. |license| raw:: html

   <a href="https://github.com/raspishake/rsudp/blob/master/LICENSE" target="_blank">license</a>

.. warning::

    **Standard performance disclaimer**

    It is extremely important that you do not rely on this code to save life or property.
    It is not a substitute for earthquake early warning (EEW), or state or local official
    communication and alert systems.

    Although this software can detect earthquakes and sudden motion events,
    Raspberry Shake makes no guarantee and provides no warranty in any way,
    implied or explicit, for the performance of this software in earthquake detection.

    Raspberry Shake assumes no liability for false positives, false negatives,
    errors running the Alert module, or any other part of this software;
    it is meant for hobby and non-professional notification use only.

    If you need professional-grade software to provide a warning intended to save life
    or property, please contact Raspberry Shake directly or look elsewhere.
    See sections 16 and 16b of the |license| for further details.


As noted in the |license|,
Raspberry Shake S.A. does not assume any liability or make any guarantee regarding
the performance of this software in detecting earthquakes.
Due to the unpredictable nature of earthquakes and the fact that this is not professional
monitoring software, this software should not be used to protect life or property.

Due to the limitations of the Raspberry Pi 3B's RAM modules, rsudp will run but occasionally
suffer from memory errors if the plotting module is enabled.
It can be programmed as a daemon in order to restart in the event of one of these errors,
which means that its monitoring capability may have brief periods of non-availability
on these devices (after :ref:`install`, see :ref:`daemon`).
The Raspberry Pi 4B does not seem to suffer nearly as much from this issue,
but this has not been tested extensively.

`Back to top ↑ <#top>`_
.. _settings:

Modules and Settings
#################################################

.. role:: bash(code)
    :language: bash

.. role:: json(code)
    :language: json

.. role:: pycode(code)
    :language: python



You will need to adjust the settings file before running :bash:`rsudp` in order to
both receive data and suit your earthquake detection and notification needs.
Additionally after changing settings, you will need to restart :bash:`rsudp` to load new values.

By default, the rsudp settings file will live in :bash:`$HOME/.config/rsudp/rsudp_settings.json`,
where :bash:`$HOME` is your home directory (often shortened to :bash:`~/`).
To output default settings to a different location, see :ref:`running-manually`.

For convenience, you can use rsudp's built-in :bash:`rs-settings` command
(:func:`rsudp.entry_points.ep_edit_settings`) to call your system's default editor to edit the settings
file.


:code:`settings` (general settings)
*************************************************

The :json:`"settings"` portion of the settings file contains some basic items:
:json:`"port"`, :json:`"station"`, :json:`"output_dir"`, and :json:`"debug"`.
Change :json:`"port"` if you are receiving the data at a different port than :json:`8888`.
To set your station name, change the value set for :json:`"station"`.
:json:`"output_dir"` will contain folders for miniSEED data and plot screenshots,
which are explained in the relevant sections (write and plot) below.
The directory specified here will be created if it doesn't already exist.
:json:`"debug"` controls how much text is sent to the command line STDOUT
(even if this is false, output will always be sent to a log at :code:`/tmp/rsudp/rsudp.log`).


:code:`plot` (live data plot)
*************************************************

:json:`"plot"` controls :class:`rsudp.c_plot.Plot`, the thread containing the GUI plotting
algorithm.
This module can plot seismogram data from a list of 1-4 Shake channels, and calculate and
display a spectrogram beneath each.

By default the plotted :json:`"duration"` in seconds is :json:`30`.
The plot will refresh at most once per second, but slower processors may take longer.
The longer the duration, the more processor power it will take to refresh the plot,
especially when the spectrogram is enabled.
To disable the spectrogram, set :json:`"spectrogram"` to :json:`false` in the settings file.
To put the plot into fullscreen window mode, set :json:`"fullscreen"` to :json:`true`.
To put the plot into kiosk mode, set :json:`"kiosk"` to :json:`true`.

.. note::

    Kiosk mode will force the plot to fill the entire screen.
    To exit, press Ctrl+W or Alt+Tab (Command+Tab on Mac OS) to bring up a window switcher).

.. note::

    On a Raspberry Pi 3B+, plotting 600 seconds of data and a spectrogram from one channel,
    the update frequency is approximately once every 5 seconds,
    but more powerful processors will be able to accommodate a higher refresh speed.

.. note::

    Because the plot module is queue-based, it will not drop any packets received, no matter the processor.
    Dropped packets (if you experience them) are most likely a sign of network issues
    where the missing data never actually arrives at the receiving machine.

By default, the :json:`"channels"` field is :json:`["HZ", "HDF"]`.
This will resolve to at least one channel of any Shake input.
:json:`"HZ"` will match either :json:`"SHZ"` or :json:`"EHZ"` depending on your Shake digitizer model,
and :json:`"HDF"` will match the pressure transducer channel on a Raspberry Boom or Shake & Boom.
If one of the channels in the list doesn't exist in the data sent to the port, it will be ignored.

The program will use the Raspberry Shake FDSN service to search for an inventory response file
for the Shake you specify in the :json:`"station"` field.
If it successfully finds an inventory,
setting "deconvolve" to :json:`true` will deconvolve the channels plotted to either :json:`"ACC"` (acceleration in m/s^2),
:json:`"VEL"` (velocity in m/s), or :json:`"DISP"` (displacement in m).
The default is :json:`"CHAN"` which lets the program deconvolve the channel
to its native units (acceleration for accelerometers, and velocity for geophones).
This means that the Shake must both have the 4.5 Hz geophone distributed by RS,
and be forwarding data to the Shake server, in order to deconvolve successfully.
For the time being, the Raspberry Boom will display in counts of Voltage, i.e., not a deconvolved unit.

If the :ref:`alert` module is enabled, setting :json:`"eq_screenshots"` to :json:`true`
will result in screenshots being saved whenever there is an :code:`ALARM`
is internally forwarded for further processing (see Alert section below).
The script will save one PNG figure per alert to the :code:`screenshots` directory
inside of :json:`"output_dir"` when the leading edge of the quake is about 70% of the way across the plot window.
This will only occur when the alarm gets triggered, however, so make sure to test your alert settings thoroughly.

`Back to top ↑ <#top>`_

.. _alert:

:code:`alert` (STA/LTA earthquake detection trigger)
*********************************************************************************

.. |license| raw:: html

   <a href="https://github.com/raspishake/rsudp/blob/master/LICENSE" target="_blank">license</a>

.. warning::

    **Standard performance disclaimer**

    It is extremely important that you do not rely on this code to save life or property.
    It is not a substitute for earthquake early warning (EEW), or state or local official
    communication and alert systems.

    Although this software can detect earthquakes and sudden motion events,
    Raspberry Shake makes no guarantee and provides no warranty in any way,
    implied or explicit, for the performance of this software in earthquake detection.

    Raspberry Shake assumes no liability for false positives, false negatives,
    errors running the Alert module, or any other part of this software;
    it is meant for hobby and non-professional notification use only.

    If you need professional-grade software to provide a warning intended to save life
    or property, please contact Raspberry Shake directly or look elsewhere.
    See sections 16 and 16b of the |license| for further details.

.. |obspy_stalta| raw:: html

   <a href="https://docs.obspy.org/tutorial/code_snippets/trigger_tutorial.html#recursive-sta-lta" target="_blank">here</a>

:json:`"alert"` controls the :class:`rsudp.c_alert.Alert` module (please see Warning above).
The alert module is a fast recursive STA/LTA sudden motion detector that utilizes obspy's
:py:func:`obspy.signal.trigger.recursive_sta_lta` function
(more detailed information on how to use that function |obspy_stalta|).
STA/LTA algorithms calculate a ratio of the short term average of station noise to the long term average.
The data can be highpass, lowpass, or bandpass filtered by changing the :json:`"highpass"`
and :json:`"lowpass"` parameters from their defaults (:json:`0` and :json:`50` respectively).
By default, the alert will be calculated on raw count data
from the vertical geophone channel (either :json:`"SHZ"` or :json:`"EHZ"`).
It will throw an error if there is no Z channel available (i.e. if you have a Raspberry Boom with no geophone).
If you have a Boom and still would like to run this module, change the default channel :json:`"HZ"` to :json:`"HDF"`.

Like in the plot module, the alert module deconvolves the instrument response if a response file exists
for your :json:`"station"` on the Raspberry Shake FDSN server.
Same as above, if the response file exists,
setting :json:`"deconvolve"` to :json:`true` will cause the alert function to
calculate the STA/LTA ratio on deconvolved data (again :json:`"ACC"`, :json:`"VEL"`, or :json:`"DISP"`).

If the STA/LTA ratio goes above a certain value (defined by :json:`"threshold"`),
then the :py:class:`rsudp.p_producer.Producer` thread will generate an :code:`ALARM` "event packet",
to be distributed to every consumer module.
This tells all consumers listening for :code:`ALARM` messages to do something.

When the ratio goes back below the :json:`"reset"` value, the alarm is reset.
The Producer will then send a :code:`RESET` message to the queues.

For more information on the packets generated by the Producer, see :ref:`producer-consumer`.

Recommendations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The detection and filtering settings that we've found work well are below for different scenarios.

General use
"""""""""""""""""""""""""""""""""""

For a station with sudden motion (footsteps nearby occasionally),
or one atop unconsolidated sediment:

.. code-block::

    "alert": {
        "enabled": true,
        "channel": "HZ",
        "sta": 6,
        "lta": 30,
        "threshold": 4.5,
        "reset": 0.5,
        "highpass": 0.8,
        "lowpass": 9,
        "deconvolve": false,
        "units": "VEL"},

Quiet vault
"""""""""""""""""""""""""""""""""""

For a very quiet station placed atop bedrock:

.. code-block::

    "alert": {
        "enabled": true,
        "channel": "HZ",
        "sta": 6,
        "lta": 30,
        "threshold": 1,
        "reset": 0.2,
        "highpass": 0.8,
        "lowpass": 9,
        "deconvolve": false,
        "units": "VEL"},

Classroom demonstrations
"""""""""""""""""""""""""""""""""""

For a classroom looking to detect jumps but not necessarily earthquakes, start with
the settings below. The main difference here is that there is no bandpass filter
applied to the signal before it is put into the STA/LTA algorithm, which changes
the calculation needed for exceedence of the threshold. Adjust the
:json:`"threshold"` downward, closer to :json:`1.7` if :json:`1.8` is too high.

.. code-block::

    "alert": {
        "enabled": true,
        "channel": "HZ",
        "sta": 6,
        "lta": 30,
        "threshold": 1.8,
        "reset": 1.6,
        "highpass": 0,
        "lowpass": 50,
        "deconvolve": false,
        "units": "VEL"},

Using :code:`"exec"`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    .. deprecated:: 0.4.3

        You can change the :json:`"exec"` field and supply a path to
        executable Python code to run with the :py:func:`exec` function.
        :py:func:`exec` functionality will move to its own module in version 0.4.3
        (see :ref:`customcode` and the :py:class:`rsudp.c_custom.Custom` class),
        and this part of the alert module will be fully removed in a future release.


`Back to top ↑ <#top>`_


:code:`RSAM` (Real-time Seismic AMplitude)
*************************************************

.. versionadded:: 1.0.1

This module calculates the Real-time Seismic Amplitude Measurement (RSAM) of the data stream every few seconds
and can forward this data to another location on the network.

:json:`"interval"` is a float that specifies the number of seconds to wait between each RSAM analysis.

:json:`"quiet"` controls the amount of data printed to the console in debug mode.
When :json:`"quiet"` is :json:`true`, the module will not print any RSAM analysis,
If debug mode is on and :json:`"quiet"` is :json:`false`, then the module will
print the analysis to the console every :json:`"interval"` seconds.

:json:`"fwaddr"` and :json:`"fwport"` specify the forwarding address and port to which to
optionally send RSAM data. If one of these fields is :json:`false` then no data will be
forwarded. If these fields are populated with valid IP and port, data will be forwarded every
:json:`"interval"` seconds.

:json:`"fwformat"` specifies the format of data to be forwarded. There are three formats,
:json:`"LITE"`, :json:`"JSON"`, and :json:`"CSV"`, which can be used depending on the
endpoint processing method and size constraints.

:json:`"channel"` specifies the channel to use for RSAM analysis (only one can be chosen).

:json:`"deconvolve"` specifies whether the instrument response should be removed from the data stream
prior to RSAM calculations.

To run the RSAM module, set :json:`"enabled"` to :json:`true`.


:code:`alarmsound` (play sounds upon alerts)
*************************************************

.. |pydub_deps| raw:: html

   <a href="https://github.com/jiaaro/pydub#dependencies" target="_blank">this page</a>

If alarmsound's :json:`"enabled"` is :json:`true` and you have either :bash:`ffmpeg` or :bash:`libav` installed,
:class:`rsudp.c_alertsound.AlertSound` plays an MP3 sound every time it receives an :code:`ALARM` queue message.
For details on installation of these dependencies, see |pydub_deps|.

The rsudp software will install several small MP3 files.
The :json:`"mp3file"` is :json:`"doorbell"` (two doorbell chimes) by default,
but there are a few more aggressive alert sounds, including: a three-beep sound :json:`"beeps"`,
a sequence of sonar pings :json:`"sonar"`,
and a continuous alarm beeping for 5 seconds, :json:`"alarm"`.
You can also point the :json:`"mp3file"` field to an MP3 file somewhere in your filesystem.
For example, if your username was :code:`pi` and you had a file called `earthquake.mp3` in your Downloads folder,
you would specify :json:`"mp3file": "/home/pi/Downloads/earthquake.mp3"`.
The program will throw an error if it can't find (or load) the specified MP3 file.
It will also alert you if the software dependencies for playback are not installed.

To test the sound output, ensure you have the correct dependencies installed (see below),
change :json:`"enabled"` to :json:`true`, start rsudp,
wait for the trigger to warm up, then stomp, jump, or Shake to trigger the sound.

Installing :code:`pydub` dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you would like to play sounds when the STA/LTA trigger activates,
you will need to take the following installation steps beforehand:

On Linux
"""""""""""""""""""""""""""""""""""""""""""""""""""""

.. |ffmpeg| raw:: html

   <a href="http://ffmpeg.org/" target="_blank">ffmpeg</a>

.. |ffmpeg_dl| raw:: html

   <a href="http://ffmpeg.org/download.html#build-mac" target="_blank">from the ffmpeg website</a>

|ffmpeg| comes installed by default on some OS flavors
and is available on most Linux package managers.

Debian and Raspbian users can simply type :bash:`sudo apt update; sudo apt install ffmpeg`

On MacOS
"""""""""""""""""""""""""""""""""""""""""""""""""""""

Users with Homebrew can install by doing :bash:`brew install ffmpeg`

Users without Homebrew will need to install using a binary build |ffmpeg_dl|.

On Windows
"""""""""""""""""""""""""""""""""""""""""""""""""""""

.. |ffmpeg_win| raw:: html

   <a href="https://windowsloop.com/install-ffmpeg-windows-10/" target="_blank">this installation guide</a>

Windows users will need to do a couple of extra steps to get :code:`ffmpeg` installed.
Following steps 1-8 in |ffmpeg_win| should be sufficient to get things working.

`Back to top ↑ <#top>`_


:code:`telegram` (Telegram notification module)
*************************************************

.. |telegram| raw:: html

    <a href="https://t.me/" target="_blank">Telegram</a>

.. |sasmex| raw:: html

    <a href="https://sasmex.net/" target="_blank">SASMEX</a>

.. |sasmex_telegram| raw:: html

    <a href="https://t.me/sasmex" target="_blank">Telegram channel here</a>


|telegram| is a free and open source messaging and notification system,
used by several earthquake notification agencies including the
Mexican national early warning system (|sasmex|, |sasmex_telegram|).
It has the bonus of being much, much easier to set up than Twitter,
and will not as readily lock your account if there happen to be many posts in a short time period
(in comparison to Twitter).

If :json:`"enabled"` is :json:`true`, and bot :json:`"token"` key is correctly entered,
:class:`rsudp.c_telegram.Telegrammer` will use the Telegram bot API to create alerts when an
:code:`ALARM` message arrives on the queue.
If :json:`"send_images"` is :json:`true`, then the module will also send a saved image of the event,
if :json:`"eq_screenshots"` is set to :json:`true` in the :json:`"plot"` module.

If any text is put in the :json:`"extra_text"` string, then the software will insert that text
(no longer than 3900 characters) into the message after the UTC designation and prior to the
stationview hyperlink.
This works similarly to the :json:`"extra_text"` field in the Twitter module below.
(See :ref:`examples`.)

.. warning::

    Starting the software with an :json:`"extra_text"` string in excess of 3900 characters
    will yield a warning and the :json:`"extra_text"` string will be truncated
    in order to avoid the message being rejected for exceeding the 4096 character limit.

.. _setting-up-telegram:

Setting up a Telegram Bot
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here is a brief overview of the steps to set up a Telegram bot in order to make and distribute
Telegram alerts from rsudp.

.. |so_answer| raw:: html

    <a href="https://stackoverflow.com/a/32572159" target="_blank">this stackoverflow answer</a>


#. Download |telegram|, create a profile, and sign in.
#. Create a Telegram bot by sending the :code:`/start` message to the :code:`@BotFather` account.
#. Follow the instructions. Your messages to :code:`@BotFather` should look something like the following:

    #. :code:`/start`

    #. :code:`/newbot`

    #. :code:`Your Shake Bot Name`

    #. :code:`your_shake_bot_id`

    #. :code:`@BotFather` will then give you an access token for your new bot.

#. Enter your bot's access token in the :json:`"token"` field of the settings file.
#. Enter a user or group ID into the :json:`"chat_id"` field (or multiple separated by commas),
    which you can find by following the instructions in |so_answer|.

If you wish to broadcast telegrams to a group or a channel, first add the bot to the group using your
user account, then follow the instructions in the previous link,
where you will see the group chat ID appear as a field in the last JSON entry.
This chat ID may be negative, in which case you must enter the negative sign into :json:`"chat_id"`
as well.


`Back to top ↑ <#top>`_


:code:`tweets` (Twitter notification module)
*************************************************

If :json:`"enabled"` is :json:`true`, and all API keys have been generated and are correctly entered,
then the :class:`rsudp.c_tweet.Tweeter` class will use the Twitter API to
create tweets when an ALARM message arrives on the queue.
If :json:`"tweet_images"` is :json:`true`, then the module will also tweet a saved image of the event,
if :json:`"eq_screenshots"` is set to :json:`true` in the "plot" module. If any text is put in the
:json:`"extra_text"` string, then the software will insert that text (no longer than 103 characters)
into the tweets after a single space. See examples below.

.. _eq-tweet-examples:

Examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As a comparison point, an unmodified tweet with :code:`"extra_text": ""` might look like
the following:

.. _eq-tweet:
.. figure::  _static/eq_tweet.png
    :align:   center

    An example tweet sent with the "extra_text" parameter empty (this is the default).


Changing the :json:`"extra_text"` parameter to :code:`"extra_text": "from #Williamstown #MA"`
would render something like this:

.. _eq-tweet-extra:
.. figure::  _static/eq_tweet_extra.png
    :align:   center

    An example tweet sent with the "extra_text" parameter filled.

.. warning::

    Starting the software with an :json:`"extra_text"` string in excess of 103 characters
    will yield a warning and the :json:`"extra_text"` string will be truncated
    in order to avoid the tweet being rejected for exceeding the 280 character limit.


.. _setting-up-twitter:

Setting up Twitter Apps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here is a brief overview of the steps to set up a Twitter app (also known as an API bot)
in order to make and distribute tweets from rsudp.

.. |tw_signup| raw:: html

    <a href="https://twitter.com/signup" target="_blank">Create a twitter profile</a>

.. |tw_dev| raw:: html

    <a href="https://developer.twitter.com/en.html" target="_blank">Twitter developer account</a>

.. |tw_api_app| raw:: html

    <a href="https://opensource.com/article/17/8/raspberry-pi-twitter-bot" target="_blank">Twitter API app</a>

#. |tw_signup| (or use an existing account).
#. Register this account as a |tw_dev|.
#. Create a |tw_api_app| inside said developer account.
#. Generate consumer keys and API keys for that app.

Once you have generated the four API keys required for authentication
(consumer API key, consumer API secret, access token, and access token secret),
you may enter them into your settings file in the appropriate fields:
:json:`"api_key"`, :json:`"api_secret"`, :json:`"access_token"`, and :json:`"access_secret"`.

`Back to top ↑ <#top>`_


:code:`write` (miniSEED writer)
*************************************************

:json:`"write"` controls :class:`rsudp.c_write.Write`, a very simple STEIM2 miniSEED writer class.
If :json:`"enabled"` is :json:`true`, seismic data is appended to a miniSEED file with a
descriptive name in the data directory inside of :json:`"output_dir"` every 10 seconds.
By default, :json:`"all"` channels will be written to their own files.
You can change which channels are written by changing this to, for example, :json:`["EHZ", "ENZ"]`,
which will write the vertical geophone and accelerometer channels from RS4D output.

`Back to top ↑ <#top>`_


.. _datacast-forwarding:

:code:`forward` (datacast forwarding)
*************************************************

The :json:`"forward"` module controls :class:`rsudp.c_forward.Forward`, a UDP datacast forwarding module.
You can forward UDP packets containing data and/or alarm state messages to a list of destinations specified
in :json:`"address"` and :json:`"port"`, just like you would from the Shake's web front end.

By default, :json:`["all"]` channels are forwarded. To forward only data from EHZ and ENZ
channels, set this field to a list, e.g. :json:`["EHZ", "ENZ"]`.

To change the types of messages that are forwarded, change the boolean fields :json:`"fwd_data"` and
:json:`"fwd_alarms"` accordingly. Setting :code:`"fwd_data": true` will forward data from the specified
channels, while :code:`"fwd_alarms": true` will forward :code:`ALARM` and :code:`RESET` messages. These can
both be set to true simultaneously.

To take advantage of this forwarding capability in another piece of software (such as NodeRED), it may help
to consult the :ref:`message-types`.

Forwarding to multiple destinations (such as in a classroom setting) is easy. Say you want to send alarm
messages to several Raspberry Pis running NodeRED in a classroom. Simply create equal-length lists of
addresses and ports in the forward settings like so::

    "forward": {
        "enabled": false,
        "address": ["192.168.1.250","192.168.1.251","192.168.1.252","192.168.1.253"],
        "port": [8888,8888,8888,8888],
        "channels": ["all"],
        "fwd_data": false,
        "fwd_alarms": true},

This will create one Forward thread per destination and distribute :code:`ALARM` and :code:`RESET`
messages to each simultaneously. Each Pi node can then be configured to listen to its own port 8888
(127.0.0.1:8888) to read these messages.

`Back to top ↑ <#top>`_


.. _customcode:

:code:`custom` (run custom code)
*************************************************

.. versionadded:: 0.4.3

.. warning:: Do not use this module unless you understand the implications of running unchecked code.

:json:`"custom"` controls the execution of a custom python code file specified by the :json:`"codefile"` field.
If :json:`"enabled"` is :json:`true` and a python file is found at the path specified,
this thread will run the specified file using python's :py:func:`exec` function.

Be very careful when using this module, as the :py:func:`exec` function is known to have problems.
Notably, :py:func:`exec` does not check the passed file for errors prior to running.
Also, the passed file cannot have Windows line endings (see warning below).
Additionally, if the code takes too long to execute,
you could end up losing data packets from the queue, so keep it simple.
Sending a message or a tweet, which should either succeed or time out in a few seconds,
is really the intended purpose, and this can typically be achieved by setting up a different module anyway
(see Twitter and Telegram modules).

In testing, we were able to run scripts with execution times of 30 seconds without losing any data packets.
Theoretically you could run code that takes longer to process than that,
but the issue is that the longer it takes the function to process code,
the longer the module will go without processing data from the queue
(the queue can hold up to 2048 packets, which for a RS4D works out to ~128 seconds of data).
Another way of saying this is: you could miss whatever subsequent earthquakes occur while :pycode:`exec()` is running.
A better way to run your own code would be to fork this repository
and create a new thread that does the thing you want when it sees an ALARM data packet on the queue.
That way, the code will be checked for errors prior to running.

.. |lineendings_howto| raw:: html

   <a href="https://stackoverflow.com/questions/17579553/windows-command-to-convert-unix-line-endings" target="_blank">this stackoverflow question</a>

.. |lineendings_wiki| raw:: html

   <a href="https://en.wikipedia.org/wiki/Newline" target="_blank">here</a>

.. warning::

    If you are running Windows and have code you want to pass to the :py:func:`exec` function,
    Python requires that your newline characters are in the UNIX style (:code:`\n`), not the standard Windows style (:code:`\r\n`).
    To convert, follow the instructions in one of the answers to |lineendings_howto|.
    If you're not sure what this means, please read about newline/line ending characters |lineendings_wiki|.
    If you are certain that your code file has no Windows newlines, you can set :json:`"win_override"` to :json:`true`.

    Code will not execute on Windows unless this field is set to :json:`true`.

`Back to top ↑ <#top>`_


:code:`printdata` (print data to console)
*************************************************

:json:`"printdata"` controls the data output module :class:`rsudp.c_printraw.PrintRaw`,
which simply prints Shake data packets to stdout as it receives them.
Change :json:`"enabled"` to :json:`true` to activate.

`Back to top ↑ <#top>`_


You are now ready to proceed to the next section, :ref:`running`.


.. _defaults:

Default settings
*************************************************

By default, the settings are as follows:

.. code-block:: json

    {
    "settings": {
        "port": 8888,
        "station": "Z0000",
        "output_dir": "@@DIR@@",
        "debug": true},
    "printdata": {
        "enabled": false},
    "write": {
        "enabled": false,
        "channels": ["all"]},
    "plot": {
        "enabled": true,
        "duration": 90,
        "spectrogram": true,
        "fullscreen": false,
        "kiosk": false,
        "eq_screenshots": false,
        "channels": ["all"],
        "deconvolve": true,
        "units": "CHAN"},
    "forward": {
        "enabled": false,
        "address": ["192.168.1.254"],
        "port": [8888],
        "channels": ["all"],
        "fwd_data": true,
        "fwd_alarms": false},
    "alert": {
        "enabled": true,
        "channel": "HZ",
        "sta": 6,
        "lta": 30,
        "threshold": 3.95,
        "reset": 0.9,
        "highpass": 0.8,
        "lowpass": 9,
        "deconvolve": false,
        "units": "VEL"},
    "alertsound": {
        "enabled": false,
        "mp3file": "doorbell"},
    "custom": {
        "enabled": false,
        "codefile": "n/a",
        "win_override": false},
    "tweets": {
        "enabled": false,
        "tweet_images": true,
        "api_key": "n/a",
        "api_secret": "n/a",
        "access_token": "n/a",
        "access_secret": "n/a",
        "extra_text": ""},
    "telegram": {
        "enabled": false,
        "send_images": true,
        "token": "n/a",
        "chat_id": "n/a",
        "extra_text": ""},
    "rsam": {
        "enabled": false,
        "quiet": true,
        "fwaddr": "192.168.1.254",
        "fwport": 8887,
        "fwformat": "LITE",
        "channel": "HZ",
        "interval": 10,
        "deconvolve": false,
        "units": "VEL"}
    }


................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
:py:data:`rsudp.c_custom` (run code)
=====================================================

.. automodule:: rsudp.c_custom
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
.. _daemon:

rsudp as a ``systemctl`` daemon
################################################

If you need rsudp to restart every time your Debain/Ubuntu machine boots,
this tutorial will work well for you.
This includes both raspbian (Raspberry Pi) and Ubuntu-like systems.
If you have an OS other than Linux, you will need to look elsewhere for
documentation regarding the creation and activation of daemon programs.


Setup instructions
=====================

First of all, clone rsudp to ``~/bin/rsudp``:

.. code-block:: bash

    mkdir -p ~/bin
    cd ~/bin
    git clone https://github.com/raspishake/rsudp

If you have not already done so, install rsudp using the provided
installer script.

.. code-block:: bash

    bash ~/bin/rsudp/unix-install-rsudp.sh

Next, enter the following commands in order to set up a user systemd
directory structure.

.. code-block:: bash

    mkdir -p ~/.config/systemd/user/default.target.wants

Now create a new service file

.. code-block:: bash

    nano /home/pi/.config/systemd/user/rsudp.service

and paste the following code in it::

    [Unit]
    Description=rsudp daemon
    Documentation=https://github.com/iannesbitt/rsudp
    After=graphical.target

    [Service]
    ExecStartPre=/bin/sleep 10
    ExecStart=/bin/bash /home/pi/bin/rsudp/unix-start-rsudp.sh
    ExecStop=kill $MAINPID
    ExecReload=kill $MAINPID; /bin/bash /home/pi/bin/rsudp/unix-start-rsudp.sh
    Restart=always
    RestartSec=10s

    [Install]
    WantedBy=default.target

Next, execute the following lines:

.. code-block:: bash

    systemctl --user daemon-reload
    sudo loginctl enable-linger "$USER"
    systemctl --user start rsudp.service

rsudp should start within about 30 seconds.
You can monitor the program's output in real time by entering the following command
(source function: :func:`rsudp.entry_points.ep_tailf_log`):

.. code-block:: bash

    rs-tailf

Which is the equivalent of ``tail -f /tmp/rsudp/rsudp.log`` on Linux/MacOS
and ``Get-Content -Path "C:/tmp/rsudp/rsudp.log" -Wait`` on Windows.

If it does start correctly, you can enable the daemon to run permanently with this command:

.. code-block:: bash

    systemctl --user enable rsudp.service

You can test its enablement by restarting the entire system:

.. code-block:: bash

    sudo reboot

Restarting the daemon
==================================

Finally, if you need to restart the rsudp daemon service
(this may be necessary if your Shake changes IP or the network connection
is interrupted, or if rsudp freezes for some reason): 

.. code-block:: bash

    systemctl --user restart rsudp.service


Troubleshooting the daemon
=================================

If rsudp fails to start, you can run ``tail -n 30 -f /tmp/rsudp/rsudp.log``
to see what the error might be, or ``systemctl --user status rsudp.service``
to check whether the service file is misconfigured somehow.

A running daemon will show its status with green text saying "active (running)",
whereas a failed start will show red or grey text that will say
something like "inactive (failed)" or "inactive (dead)"
and will have some diagnostic text with which you can troubleshoot.


`Back to top ↑ <#top>`_

.. _test:

Testing and demoing rsudp
#################################################

rsudp includes a small piece of software meant to test the
ability to function and process seismic data.
This could be useful if you are looking to demonstrate rsudp's
functionality to someone, perhaps a classroom of students,
or if you need to test the core functionality like alerts
and sounds.

Testing can also be useful for discovering whether or not a specific
piece of data will trigger the alarm using settings from a custom
settings file. For instructions on how to do that, see
:ref:`test_settings` and :ref:`custom_data` below.

.. note::

    The testing functions are useful for figuring out local problems.
    That is, the testing capabilities of rsudp are meant to discover
    whether or not the software can feed data to `itself` and
    process it.

    If you can run this testing program without any problems
    but you are having issues getting the software to see data from
    your own Raspberry Shake or boom, see the :ref:`troubleshooting`
    page.


Using the testing functionality
=================================================

The testing modules of this software are designed to read a small
data file from disk and send its contents to the test port one
line at a time. The program functions essentially as it would if
it were receiving remote data, but instead it is feeding data
to itself.

This means that you can demo the software even if you don't have
a Raspberry Shake, or even use the testing functionality to check
whether or not an arbitrary piece of archival Raspberry Shake
data will trigger an alarm in the software.

To run this software, make sure you have installed this software
using the instructions in :ref:`install`, and that you can enter
your conda environment (see :ref:`run-test`).

Once you have done that, the test command ``rs-test`` will become
available.

Type ``rs-test`` to watch earthquake detection in
action. The test will last about 120 seconds, over which time
various bits of functionality will be tested, including ports,
directory permissions, internet, processing routines,
alert functionality, sound-playing capability, and more.

It does test whether it can see the internet at large,
and whether it can send data to its own port
(we've chosen 18888 as a test port).
However, it does not test the ability to receive data from a
remote shake. If you are having trouble with that, please see the
:ref:`troubleshooting` page.


.. _test_settings:

Settings during testing
=================================================

Default settings are slightly different during testing than they would
ordinarily be. Find a summary of what gets changed at
:py:func:`rsudp.test.make_test_settings`.

To specify a settings file to use, use the ``-s`` flag. This is the same
as it would be if you were telling the ``rs-client`` to start with a
specific settings file. Usage looks like this:

.. code-block:: bash

    rs-test -s custom_settings.json

.. note::

    If you need to dump and edit a custom settings file to test with, you can
    use the client's settings dump:

    .. code-block:: bash

        rs-client -d custom_settings.json


.. _testing_flow:

Data flow during testing
=================================================

During testing, the typical data flow as depicted in
:ref:`flow` must be created artificially.
So, instead of getting data from the Raspberry Shake as usual,
the :py:class:`rsudp.t_testdata.TestData` thread reads a file and
sends the individual lines in that file to the data port.
The Producer then reads that data and data flow through the rest
of the architecture continues as normal.

.. _test_diagram:
.. figure::  _static/test-flow.png
    :align:   center

    Flow chart of test data hierarchy,
    based on the :ref:`flow` diagram, showing how data
    makes its way through the program during testing.
    Note that there is no Raspberry Shake in the hierarchy
    as there would be in ordinary operation, but instead
    data is generated from a text file at
    ``rsudp/test/testdata``.


Testing your own modules
=================================================

Read about adding testing capabilities to new modules in
:ref:`add_testing`.


.. _custom_data:

Using your own data
=================================================

.. |canread| raw:: html

   <a href="https://docs.obspy.org/packages/autogen/obspy.core.stream.read.html#supported-formats" target="_blank">can read</a>


Included in this software is a function that will convert
small seismic data files (basically anything that obspy |canread|
that was recorded at the same sample rate as the Raspberry Shake)
to the UDP packet format required by rsudp.

This function is documented at :py:func:`rsudp.packetize.packetize`
and it is integrated into the testing script. You can tell the testing
script to convert and use a miniSEED file on disk by doing the following:

.. code-block:: bash

    rs-test -f test.mseed

This will create a text file named ``test.mseed.txt`` in the same directory
which will be used to feed data to the producer during testing.

`Back to top ↑ <#top>`_
:py:data:`rsudp.helpers` (helper functions)
=====================================================

These are some helper functions in rsudp.

.. automodule:: rsudp.helpers
    :members:


`Back to top ↑ <#top>`_
:py:data:`rsudp.c_alert` (STA/LTA alarm)
=====================================================

.. automodule:: rsudp.c_alert
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
:py:data:`rsudp.c_tweet` (Twitter alerts)
=====================================================

.. automodule:: rsudp.c_tweet
    :members:

................

* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

`Back to top ↑ <#top>`_
