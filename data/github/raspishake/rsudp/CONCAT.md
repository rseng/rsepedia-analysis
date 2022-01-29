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
