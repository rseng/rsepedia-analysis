<h2 align="center">Octopus Sensing</h2>
<p align="center">
  <img src="https://octopus-sensing.nastaran-saffar.me/_static/octopus-sensing-logo-small.png" alt="Octopus Sensing Logo">
</p>
<p align="center">
  <img src="https://img.shields.io/github/workflow/status/octopus-sensing/octopus-sensing/Python%20Check?label=checks" alt="GitHub Workflow Status">
  <img src="https://img.shields.io/codecov/c/gh/octopus-sensing/octopus-sensing" alt="Codecov">
  <img src="https://img.shields.io/pypi/v/octopus-sensing" alt="PyPI">
  <img src="https://img.shields.io/pypi/l/octopus-sensing" alt="PyPI - License">
</p>

Octopus Sensing is a tool to help you run scientific experiments that involve recording data synchronously from
multiple sources in human-computer interaction studies. You write steps of an experiment scenario, for example showing a stimulus and then a questionnaire. The tool takes care of the rest.

It can collect data from multiple devices such as OpenBCI EEG headset, Shimmer sensor (GSR and PPG),
Video and Audio and so forth simultaneously. Data collection can be started and stopped synchronously across all devices.
Collected data will be tagged with the timestamp of the start and stop of the experiment, the ID of
the experiment, etc.

The aim is to make the scripting interface so simple that people with minimum or no software
development skills can define experiment scenarios with no effort.
Also, this tool can be used as the base structure for creating real-time data processing systems like systems with capabilities of recognizing emotions, stress, cognitive load, or analyzing human behaviors.


**To see the full documentation visit the [Octopus Sensing website](https://octopus-sensing.nastaran-saffar.me/).**

Main features
--------------

* Controls data recording from multiple sources using a simple unified interface
* Tags an event on collected data, such as the start of an experiment, and events during the experiment, etc.
* Can show stimuli (images and videos) and questionnaires
* Monitoring interface that visualizes collected data in real-time
* Offline visualization of data from multiple sources simultanously

Copyright
---------
Copyright © 2020-2022 Nastaran Saffaryazdi, Aidin Gharibnavaz

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

See [License file](https://github.com/nastaran62/octopus-sensing/blob/master/LICENSE) for full terms.
# Octupus Sensing - Unity Example
In this example, We're going to send predefined messages (START, STOP and, TERMINATE) from Unity to Octopus sensing, to control the recording of data. A python server running the Octopus sensing library is used and Unity (C#) is used as the client.

Each button has the ability to send a specific message when pressed. See the UI below.

![](ui-unity.jpg)

## References
- [HTTP Endpoint](https://octopus-sensing.nastaran-saffar.me/api/octopus_sensing.device_message_endpoint.html)
- [Predefined messages](https://octopus-sensing.nastaran-saffar.me/api/octopus_sensing.device_message_endpoint.html)---
title: 'Octopus Sensing: A Python library for human behavior studies'
tags:
 - Python
 - Javascript
 - Human-Computer-Interaction(HCI)
 - Human behavior research
 - Physiological Signals
 - Electroencephalography
 - Multimodal Sensors
 - Synchronous data acquisition
 - Data visuaization
 - Affective Computing
 - Experimental design

authors:
 - name: Nastaran Saffaryazdi
   orcid: 0000-0002-6082-9772
   affiliation: 1
 - name: Aidin Gharibnavaz
   orcid: 0000-0001-6482-3944
   affiliation: 3
 - name: Mark Billinghurst
   orcid: 0000-0003-4172-6759
   affiliation: 1, 2
affiliations:
- name: Empathic Computing Laboratory, Auckland Bioengineering Institute, University of Auckland
  index: 1
- name: Empathic Computing Laboratory, University of South Australia
  index: 2
- name: Independent Researcher
  index: 3
date: 1 December 2021
bibliography: paper.bib

---

# Summary
Designing user studies and collecting data is critical to exploring and automatically recognizing human behavior. It is currently possible to use a range of sensors to capture heart rate, brain activity, skin conductance, and a variety of different physiological cues [@seneviratne2017survey]. These data can be combined to provide information about a user's emotional state [@egger2019emotion; @dzedzickis2020human], cognitive load [@vanneste2021towards; @mangaroska2022exploring], or other factors. However, even when data are collected correctly, synchronizing data from multiple sensors is time-consuming and prone to errors. Failure to record and synchronize data is likely to result in errors in analysis and results, as well as the need to repeat the time-consuming experiments several times.
To overcome these challenges, `Octopus Sensing` facilitates synchronous data acquisition from various sources and provides some utilities for designing user studies, real-time monitoring, and offline data visualization. The primary aim of `Octopus Sensing` is to provide a simple scripting interface so that people with basic or no software development skills can define sensor-based experiment scenarios with less effort.

# Statement of need
Several changes occur in the body and mind due to various internal and external stimuli. Nowadays, researchers use various sensors to measure and monitor these responses to determine an individual's state [@kreibig2010autonomic; @chen2021physiological; @sun2020multimodal] and to assist patients [@hassouneh2020development] or monitor mental health [@jiang2020snapshot]. Monitoring and analyzing human responses can be used to improve social interactions [@verschuere2006psychopathy; @hossain2019observers] and improve quality of life by creating intelligent devices such as Intelligent Tutoring Systems [@dewan2019engagement], creating adaptive systems [@aranha2019adapting], or creating interactive robots and virtual characters [@val2020affective; @hong2021multimodal].

Researchers have recently attempted to gain a deeper understanding of humans by simultaneously studying physiological and behavioral changes in the human body [@shu2018review; @koelstra2011deap]. Acquiring and analyzing data from different sources with various hardware and software is complex, time-consuming, and challenging. Additionally, human error can easily affect synchronously recording data in multiple formats. These tasks decrease the pace of progress in human-computer interaction and human behavior research.

There are only a few frameworks that support synchronous data acquisition and design. [iMotions](https://imotions.com/) has developed software for integrating and synchronizing data recording through a wide range of various sensors and devices. Despite having many great features, iMotions is commercial software and not open-source. In contrast, there are a few open-source programs for conducting human studies. [Psychopy](https://www.psychopy.org/) [@peirce2019psychopy2] is a powerful open-source, cross-platform software that is mainly used for building experiments' pipelines in behavioral science with visual and auditory stimuli. It can also record data from a few devices and send triggers to them. Another effort in this area is [LabStreamingLayer (LSL) LabRecorder](http://labstreaminglayer.org/). Although LSL LabRecorder provides synchronized, multimodal streaming through a wide range of devices, an extra application still needs to be run for acquiring data from each sensor separately.

`Octopus Sensing` is a lightweight open-source multi-platform library that facilitates synchronous data acquisition from various sources through a unified interface and could be easily extended to process and analyze data in real-time. We designed the `Octopus Sensing` library to minimize the effect of network failure in synchronous data streaming and reduce the number of applications that we should run for data streaming through different devices. Rather than creating a standalone software or framework, we created a library that could be easily integrated with other applications. `Octopus Sensing` provides a real-time monitoring system for illustrating and monitoring signals remotely using a web-based platform. The system also offers offline data visualization to see various human responses simultaneously.

# Overview

`Octopus Sensing` is a tool to help in running scientific experiments that involve recording data synchronously from multiple sources. It can simultaneously collect data from various devices such as [OpenBCI EEG headset](https://openbci.com/), [Shimmer3 sensor](https://shimmersensing.com), camera, and audio-recorder without running another software for data recording. Data recording can be started, stopped, and triggered synchronously across all devices through a unified interface.

The main features of `Octopus Sensing` are that it:

* manages data recording from multiple sources using a simple unified interface,
* minimizes human errors from manipulating data in synchronous data collection,
* provides some utilities for designing studies like showing different stimuli or designing questionnaires,
* offers a monitoring interface that prepares and visualizes collected data in real-time, and
* provides offline visualization of data from multiple sources simultaneously.

# Methodology
`Octopus Sensing` synchronizes data recording by using `multiprocessing` in Python. By instantiating the `Device` class, `Octopus Sensing` creates a process for the device. Each device's process has three threads: `data acquisition` thread, `trigger handling and data recording` thread, and `real-time data` thread. The `data acquisition` thread is responsible for acquiring data from a sensor. The `trigger handling and data recording` thread handles trigger messages through a message queue for synchronous data recording. It also records data in a file or files. The `real-time data` thread listens on a queue for requests and returns the last three seconds of the recorded data in the same queue. This data is being used in real-time monitoring and can be used for real-time processing and creating real-time feedback in the future. The `Device Coordinator` sends different triggers such as the start of recording or end of recording to different devices by putting the message in all devices' trigger queues at the same time. The `Device Coordinator` can also send the trigger over the network for devices that are not embedded in the `Octopus Sensing`. The following image shows the overall view of the main components of the `Octopus Sensing` and their relations.
-![Ovrall view of Octopus Sensing](OCS-diagram.png)

# Research perspective
We used `Octopus Sensing` to design several human emotion recognition user studies. We developed a user study using `Octopus Sensing` for recording facial video, brain activity, and physiological signals during a watching video task. The recorded data was used to make multimodal emotion recognition models [@Saffaryazdi2022emotion]. This scenario which is common in physiological emotion recognition studies has been included in the repository as an example and explained in the tutorial. In another study, we collected multimodal data during face-to-face conversations to make models for emotion recognition during interactive tasks [@Saffaryazdi2022conv]. We developed this user study's scenario in Python using the `Octopus Sensing` library and conducted the study only by running our developed program, without running any other software or any supervision for data recording or data synchronization.

This tool can be used to build real-time data processing systems to recognize emotions, stress, cognitive load, or analyze human behavior. Our final goal is to extend its capabilities to provide real-time emotion recognition using multimodal data. Furthermore, we plan to integrate it with `Psychopy` in the future and combine multimodal data collection and monitoring with `Psychopy` features when designing scenarios. Additionally, we plan to support LSL in the future. By supporting LSL, other applications that already support LSL could work with `Octopus Sensing`.


# Acknowledgement
We acknowledge the [Empatic Computing Laboratory (ECL)](http://empathiccomputing.org/) for financial support and for providing feedback, and Professor Suranga Nanayakkara for encouragement and feedback.

# References
.. _octopus_sensing_monitoring:

***************************
Octopus Sensing Monitoring
***************************

A web-based real-time monitoring for `Octopus Sensing <https://octopus-sensing.nastaran-saffar.me/>`_. You can
monitor your data from any machine in the same network.

`Octopus Sensing monitoring <https://github.com/octopus-sensing/octopus-sensing-monitoring>`_ is
a separated project and can be installed for Octopus Sensing if we need monitoring features.


Installation
------------

It required Python 3.7 or later. And it needs to be installed on the same machine where `Octopus
Sensing` is running.

You can use `pip` to install it:

.. code-block:: bash

   $ pip install octopus-sensing-monitoring


Then simply run it by invoking `octopus-sensing-monitoring` from the command line.

You can also use one of the Python package managers like `pipenv <https://pipenv.pypa.io/en/latest/>`_
or `poetry <https://python-poetry.org/>`_ to prevent package conflict.

.. code-block:: bash

   $ pipenv install octopus-sensing-monitoring
   $ pipenv run octopus-sensing-monitoring


The monitoring will listen on `8080` port. Open a web page and point to the machine's IP. For
example, in the same machine, open http://localhost:8080 . Or replace `localhost` with the machine's
IP and open it from any other machine.

Starting endpoint in your code
------------------------------

In the code that running Octopus Sensing, you need to start the monitoring endpoint as well. To do so, add this to your code:

>>> from octopus_sensing.device_coordinator import DeviceCoordinator
>>> from octopus_sensing.monitoring_endpoint import MonitoringEndpoint
>>> # Create coordinator instance
>>> coordinator = DeviceCoordinator()
>>> # Add your devices
>>> ...
>>> # Creating the endpoint instance and start it.
>>> monitoring_endpoint = MonitoringEndpoint(coordinator)
>>> monitoring_endpoint.start()
>>> ...
>>> # It's a good idea to stop it after your software terminated.
>>> monitoring_endpoint.stop()


Testing with fake data
----------------------

For testing purposes, you can ask the server to generate fake data instead of fetching data from
`Octopus Sensing`. To do so, add `--fake` flag when running the script:

.. code-block:: bash

   $ octopus-sensing-monitoring --fake

Naming your devices
-------------------

In `Octopus Sensing`, when you're creating an instance of devices, you need to provide a `name`. At the
moment, device names are hard coded in this monitoring app. So you need to use these names for your
devices in order for them to appear on the web page.

* For OpenBCIStreaming or BrainFlowOpenBCIStreaming use `eeg` (i.e. `OpenBCIStreaming(name="eeg", ...)` )
* For Shimmer3Streaming use `shimmer`
* For the camera, you need to create instance of `CameraStreaming` and name it `webcam`

Security notice
---------------

Note that the webserver accepts requests from any machine, and it uses `http` protocol, which
is not encrypted. Don't run it on a network that you don't trust.

**Copyright**

Copyright © 2020,2021 `Aidin Gharibnavaz <https://aidinhut.com>`

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

See `License file <LICENSE>` for full terms.
.. _visualizer_config_guide:

*******************************************
How to prepare a config file for Visualizer
*******************************************

A config file includes one section for defining specific configs for each kind of data.

EEG
--------------------
To display EEG data, first of all we should specify the path of recorded EEG file by giving value to the `path` option. 
Also, we should identify the `sampling_rate` of recorded data. 
Octopus-sensing-visualizer has provided several options for displaying raw or processed EEG data.
By giving `True` value to each option, we specify which kind of graphs to be displayed.

Octopus-sensing-visualizer supports the following graphs for displaying the EEG signals:
    - `display_signal`: If True, displays raw signals for all channels
    - `display_power_band_bars`: If True, displays a bar chart of average of each power band signal
    - `display_alpha_signal`: If True, displays alpha band signal extracted from all channels.
    - `display_beta_signal`: If True, displays beta band signal extracted from all channels.
    - `display_gamma_signal`: If True, displays gamma band signal extracted from all channels.
    - `display_delta_signal`: If True, displays delta band signal extracted from all channels.
    - `display_theta_signal`: If True, displays theta band signal extracted from all channels.
    - `window_size`: It specifies the size of window for measuring power bands in seconds
    - `overlap`: Shows the overlap between consequences windows in measuring power bands

Example
"""""""

>>> [EEG]
>>> path=data/OpenBCI_01_01.csv
>>> sampling_rate=128
>>> display_signal=False
>>> display_power_band_bars=True
>>> display_alpha_signal=False
>>> display_beta_signal=False
>>> display_gamma_signal=False
>>> display_delta_signal=False
>>> display_theta_signal=False
>>> window_size=3
>>> overlap=2

GSR
----
All options related to GSR signal. 
    - `path`: Path to GSR data (csv file)
    - `sampling_rate`: recording sampling rate
    - `display_signal`: If True, displays GSR raw signal
    - `display_phasic`: If True, extracts phasic component and displays it
    - `display_tonic`: If True, extracts tonic component and displays it

Octopus-sensing-visualizer uses `neurokit library <https://neurokit.readthedocs.io/en/latest/>`_ 
for extracting GSR components.

Example
"""""""

>>> [GSR]
>>> path=data/gsr-01-01.csv
>>> sampling_rate=128
>>> display_signal=True
>>> display_phasic=True
>>> display_tonic=True
```

PPG
----
All options related to PPG signal. 

    - `path`: Path to PPG data (csv file)
    - `sampling_rate`: recording sampling rate
    - `display_signal`: If True, displays PPG raw signal
    - `display_hr`: If True, extracts heart rate (hr) and displays it
    - `display_hrv`: If True, extracts heart rate variability (hrv) and displays it
    - `display_breathing_rate`: If True, extracts breathing rate (br) and displays it
    - `window_size`: The window size for extracting hr, hrv and br.
    - `overlap`: Shows the overlap between consequences windows

Octopus-sensing-visualizer uses `heartpy library <https://github.com/paulvangentcom/heartrate_analysis_python>`_ 
for extracting hr, hrv and breathing rate.

Example
"""""""

>>> path=octopus_sensing_visualizer/test_data/ppg_video-43-00-08.csv
>>> sampling_rate=128
>>> display_signal=True
>>> display_hr=True
>>> display_hrv=True
>>> display_breathing_rate=True
>>> window_size=20
>>> overlap=19





.. _contribution:

***********************
Contribution Guideline
***********************

Your contribution in Octopus Sensing could be by reporting bugs, requesting new features and code development.

Bug report and requesting new features
======================================
To report bugs and request new features, please create a new issue. To make sure the bug has not already been reported, please search under `GitHub Issues <https://github.com/octopus-sensing/octopus-sensing/issues>`_.
When you find no open issue that addresses your problem, create a new one. Make sure to include a title, clear description, as much relevant information as possible, as well as a code sample or an executable test case showing the expected bug happens.

Code Development
=================
To learn about the code convention, code documentation and how to contribute in code development
see :ref:`development`. 

Submit changes
==============
Please send a pull request with a list of what you have done and an appropriate comment about your changes.
To know more about pull requests, read `pull requests document <https://docs.github.com/en/pull-requests>`_ 
                        
support
========
Please contact Nastaran (nsaffar@gmail.com) or Aidin (aidin@aidinhut.com) if you have any questions or seek support.
         
     
Thanks!

Octopus Sensing Team.. _tutorial:

*************
Tutorial
*************


What Are We Building?
----------------------

This tutorial will show how to design a simple scenario with octopus-sensing step by step.

The example scenario is the most common in emotion recognition research in affective computing. In this scenario, we learn how to record data from different sources synchronously when an event happens and stop data recording by finishing the event.

**By following these examples, we learn how to:**

    1. Record data from various sources synchronously.
    2. Being synchronized with other software like Matlab and unity.
    3. Running the scenario and creating triggers in another application and recording data
        synchronously using Octopus Sensing
    4. Use various kinds of stimuli in octopus-sensing.
    5. Providing some utilities for designing experiments.
    6. Monitor and data in real-time.
    7. Reading recorded data in real-time
    8. Preprocess and visualize data offline.
    9. Watching video scenario

**Prerequisites**

Create a project and install the `octopus-sensing` package by following the instructions on :ref:`quick_start`. We recommend using `pipenv` to do so.
And then, copy the source of examples from `examples` package in octopus-sensing repository to your project directory and run them.

1- Record data from various sources synchronously
-------------------------------------------------
The most crucial feature of octopus-sensing is synchronous data recording from different sensors.
Octopus-sensing supports a set of sensors with a python library for data streaming.
Also, it supports synchronous data recording using other software like Matlab and Unity.
In this section, we learn how to record data from different sensors with internal drivers.
(devices with a python driver for data acquisition).

Adding a sensor
""""""""""""""""
Imagine you want to record video using your built-in webcam by pressing a key on the keyboard
and stop recording after 5 seconds.

Firstly we should create a CameraStreaming object with a specific name and an output path for recording data.

>>> from octopus_sensing.devices import CameraStreaming
>>> from octopus_sensing.device_coordinator import DeviceCoordinator
>>> from octopus_sensing.common.message_creators import start_message, stop_message
>>> my_camera = CameraStreaming(camera_no=0,
...                             name="camera",
...                             output_path="./output")

Then we should add the created object to the `DeviceCoordinator`. As the name suggests, the device coordinator is responsible for coordination, like starting to record data in all devices at once, stopping data recording, triggering (marking data at a point), and terminating devices. When a device is added to the device coordinator, it will be initialized and prepared for recording.

>>> device_coordinator = DeviceCoordinator()
>>> device_coordinator.add_devices([my_camera])

We are now developing a simple code to start data recording by pressing a key and stopping recording after 5 seconds.

>>> input("Press a key to start data recording")
>>> device_coordinator.dispatch(start_message(experiment_id, stimuli_id))
>>> time.sleep(5)
>>> device_coordinator.dispatch(stop_message(experiment_id, stimuli_id))
>>> time.sleep(0.5)
>>> device_coordinator.terminate()

Octopus-sensing provides a set of default messages for handling different actions like
starting and stopping the recording or terminating the program.
To identify recorded files, Octopus-sensing needs an experiment ID and stimulus ID.
They are two strings and can be anything you want.
For example, we use the id of the recorded subject as experiment ID.
Defining stimulus ID is essential for identifying the recorded data related to each stimulus
when we have different stimuli.

To see the completed example see `add_sensors example <https://github.com/octopus-sensing/octopus-sensing/blob/master/examples/add_sensors.py>`_.
The name of the recorded file will be `camera-{experiment_id}.avi` and will be saved in `output/camera` path.

Adding several sensors
""""""""""""""""""""""

To add each sensor, we should first create an instance of it and then add it to the device coordinator device list.
The device coordinator will manage synchronous data recording by sending some markers to all devices in its device_list.

>>> from octopus_sensing.devices import Shimmer3Streaming
>>> from octopus_sensing.devices import CameraStreaming
>>> from octopus_sensing.devices import BrainFlowOpenBCIStreaming
>>> from octopus_sensing.device_coordinator import DeviceCoordinator
>>> from octopus_sensing.common.message_creators import start_message, stop_message
>>> my_shimmer = Shimmer3Streaming(name="shimmer",
...                                saving_mode=SavingModeEnum.CONTINIOUS_SAVING_MODE,
...                                output_path="./output")
>>> my_camera = CameraStreaming(camera_no=0,
...                             name="camera",
...                             output_path="./output")
>>> my_openbci =
...     BrainFlowOpenBCIStreaming(name="OpenBCI",
...                               output_path="./output",
...                               board_type="cyton-daisy",
...                               saving_mode=SavingModeEnum.CONTINIOUS_SAVING_MODE,
...                               channels_order=["Fp1", "Fp2", "F7", "F3",
...                                               "F4", "F8", "T3", "C3",
...                                               "C4", "T4", "T5", "P3",
...                                               "P4", "T6", "O1", "O2"])
>>> device_coordinator.add_device(my_shimmer)
>>> device_coordinator.add_devices([my_openbci, my_shimmer, my_camera])
>>> input("Press a button to start data recording")
>>> device_coordinator.dispatch(start_message(experiment_id, stimuli_id))
>>> time.sleep(5)
>>> device_coordinator.dispatch(stop_message(experiment_id, stimuli_id))
>>> device_coordinator.terminate()

By running this example, according to the `saving_mode` option that we passed to Shimmer3Streaming and  BrainFlowOpenBCIStreaming,
the recorded file/s will be different. The default value of saving mode is continuous.
It means if we have several stimuli, all data will be recorded in one file and only some markers indicate where the event happened. In the SEPARATED_SAVING_MODE the data recorded during each stimulus will be recorded in a separate file.
In the recorded file for Shimmer3 and OpenBCI, data samples have been recorded from when the sensor initialized to when it received the terminate message.
The last column of data is the trigger column, which shows in what sample and time the device has received the start and stop triggers 
(pressing the button and 5 seconds after that). If we change the saving mode to separate (`SavingModeEnum.SEPARATED_SAVING_MODE`), it will record one file for each stimulus (For this example, one file), and the name of stimuli will appear in the file name.

Octopus Sensing can simultaneously record data from several cameras, an audio recorder, and several Shimmer3 OpenBCI sensors.
To learn more about supported sensors, see :ref:`devices`.

**Troubleshooting**

Keep in your mind, before running the code, connect the OpenBCI USB dongle, turn on the OpenBCI board. Also, turn on the Shimmer3 sensor and pair Bluetooth and the serial port for Shimmer3 streaming.
(Shimmer password: 1234)

For example, in Linux, you can do it as follow:
    1. hcitool scan   //It shows the mac-address of the device. for shimmer it is 00:06:66:F0:95:95
    2. vim /etc/bluetooth/rfcomm.conf write the below line in it: rfcomm0{ bind no; device 00:06:66:F0:95:95; channel 1; comment "serial port" }
    3. sudo rfcomm connect rfcomm0 00:06:66:F0:95:95 // This is for reading Bluetooth data from a serial port


2- Synchronization with other software
---------------------------------------
Octopus Sensing also can send synchronization markers to external devices which record data through other
software like `Matlab <https://au.mathworks.com/products/matlab.html>`_.

First, we should create an instance of `SocketNetworkDevice` and allocate an IP address and port.
Then add it to the `DeviceCoordinator` like other devices. By adding it to the `DeviceCoordinator`, it will start
listening on specified IP address and port.

>>> from octopus_sensing.devices.socket_device import SocketNetworkDevice
>>> socket_device = SocketNetworkDevice("0.0.0.0", 5002)
>>> device_coordinator.add_devices([socket_device])

Then a client can connect to this server to receive triggers. In the following code, we created a simple scenario
that sends several triggers to a simple data recorder in Matlab.

**Server Code in python**

By running the server code, it starts listening. Before to begin sending markers, make sure
that client code is running, and it has connected to the server.
See the complete example in `send_trigger_to_remote_device example <https://github.com/octopus-sensing/octopus-sensing/blob/master/examples/remote_device_example/send_trigger_to_remote_device.py>`_.

>>> from octopus_sensing.device_coordinator import DeviceCoordinator
>>> from octopus_sensing.devices import SocketNetworkDevice
>>> from octopus_sensing.common.message_creators import start_message, stop_message
>>> device_coordinator = DeviceCoordinator()
>>> socket_device = SocketNetworkDevice("0.0.0.0", 5002)
>>> device_coordinator.add_devices([socket_device])
>>> time.sleep(2)
>>> input("If a client has connected successfully, press enter to start sending marker")
>>> message = start_message("test", "00")
>>> device_coordinator.dispatch(message)
>>> time.sleep(2)
>>> message = stop_message("test", "00")
>>> device_coordinator.dispatch(message)
>>> time.sleep(2)
>>> message = start_message("test", "01")
>>> device_coordinator.dispatch(message)
>>> time.sleep(2)
>>> message = stop_message("test", "01")
>>> device_coordinator.dispatch(message)
>>> time.sleep(3)
>>> device_coordinator.terminate()

**Client Code in Matlab**

We created a simple data recorder in this example which, in parallel, listens to the network.
By running matlabRecorder in Matlab, firstly, it tries to connect to the specified server.
Then it starts listening to specified port asynchronously. Parallel to this, it is recording some numbers in a file.
As soon as it receives a marker, it will add it to the recorded line in the file.
See this example in `matlabRecorder example <https://github.com/octopus-sensing/octopus-sensing/blob/master/examples/remote_device_example/matlabRecorder.m>`_.


>>> function matlabRecorder()
>>>     global marker
>>>     marker = "";
>>>     tcpipClient = tcpip('localhost',5002,'NetworkRole','Client');
>>>     tcpipClient.ReadAsyncMode = 'continuous';
>>>     tcpipClient.Terminator = 10;
>>>     tcpipClient.BytesAvailableFcn = @setMarker;
>>>     tcpipClient.BytesAvailableFcnMode = 'terminator';
>>>     fopen(tcpipClient);
>>>     file_out = fopen("file_out.csv", 'w');
>>>     i = double(0);
>>>     while(1)
>>>         if marker == "terminate"
>>>             break
>>>         elseif marker == ""
>>>             fprintf(file_out, "%d, %s\n", i, "");
>>>         else
>>>             fprintf(file_out, "%d,%s\n", i, marker);
>>>             marker = "";
>>>         end
>>>         i =  i + 1;
>>>         pause(0.1);
>>>     end
>>>     fclose(file_out);
>>>     fclose(tcpipClient)
>>>
>>> end
>>>
>>> function setMarker(obj, event)
>>>     global marker;
>>>     data = fscanf(obj);
>>>     marker = erase(data, char(10));
>>> end


3- Receiving Messages over Network
-----------------------------------
Octopus Sensing provides an endpoint that listens for incoming Message requests by starting it.
It passes the message to the Device Coordinator to dispatch them to the devices.
It accepts HTTP POST requests. The Body can be serialized in one of 'json', 'msgpack'
or 'pickle'.
This feature can be used when we have designed the overall scenario with other programming languages or the scenario
is running in other software like Unity or Matlab. In this cases, we should write a simple code in python that uses
Octopus Sensing for data recording and our scenario will just send triggers as an http request.

On the server-side first of all, we should create the device_coordinator and add the desired devices to it. Then we should
create an endpoint as follows, pass the DeviceCoordinator instance to it and start it.

>>> from octopus_sensing.device_message_endpoint import DeviceMessageHTTPEndpoint
>>> message_endpoint = DeviceMessageHTTPEndpoint(device_coordinator, port=9331)
>>> message_endpoint.start()

An HTTP server will be started by running this code, which is listening on port 9331.
When it receives a trigger, it passes it to the DeviceCoordinator, and DeviceCoordinator
dispatches it to all the added devices.

On the client-side, if the language is python, we should first connect to the server
by giving the machine's address and the specified port of the server. In this example, we provide the
address of the local machine because both client and server is running on the same machine

>>> import msgpack
>>> import http.client
>>> http_client = http.client.HTTPConnection("127.0.0.1:9331", timeout=3)

Then we can send a message as follows:

>>> http_client.request("POST", "/",
...                     body=msgpack.packb({'type': 'START',
...                                         'experiment_id': experiment_id,
...                                         'stimulus_id': stimuli_id}),
...                     headers={'Accept': 'application/msgpack'})
>>> response = http_client.getresponse()
>>> assert response.status == 200

See the full example in `endpoint_example <https://github.com/octopus-sensing/octopus-sensing/tree/master/examples/endpoint_example>`_.


4- Use various kinds of stimuli in octopus-sensing
--------------------------------------------------
In this example, we learn how to record data in parallel with displaying image stimuli.

To display stimuli, Octopus-Sensing provides a set of predefined stimuli, including video and image.
To display image stimuli, we used `GTK <https://athenajc.gitbooks.io/python-gtk-3-api/content/>`_. We should specify the path of the image stimulus and the duration time
for displaying it. See :ref:`stimuli` for stimuli API documentation.

>>> from octopus_sensing.stimuli import ImageStimulus
>>> stimulus = ImageStimulus(stimuli_id, os.path.join(stimuli_path, stmulus_name), 5)
>>> stimulus.show_standalone()

Similarly, we can create a video stimulus. Octopus Sensing uses
`VLC media player <https://www.videolan.org/vlc/>`_ to display video stimuli.
You should have VLC installed on your system.

>>> from octopus_sensing.stimuli import VideoStimulus
>>> stimulus = VideoStimulus(stimuli_id, os.path.join(stimuli_path, stmulus_name))
>>> stimulus.show()

The following code is the complete example of recording physiological data using Shimmer3
sensor while a set of images are displaying. See `simple_scenario example <https://github.com/octopus-sensing/octopus-sensing/blob/master/examples/simple_scenario.py>`_. 
In this example, you can have video stimuli with uncommenting video stimuli lines and commenting image stimuli lines.

>>> import time
>>> import os
>>> from octopus_sensing.devices import Shimmer3Streaming
>>> from octopus_sensing.device_coordinator import DeviceCoordinator
>>> from octopus_sensing.common.message_creators import start_message, stop_message
>>> from octopus_sensing.stimuli import ImageStimulus
>>>
>>>
>>> def simple_scenario(stimuli_path):
>>>     # Reading image stimuli and assigning an ID to them based on their alphabetical order
>>>     stimuli_list = os.listdir(stimuli_path)
>>>     stimuli_list.sort()
>>>     stimuli = {}
>>>     i = 0
>>>     for item in stimuli_list:
>>>         stimuli[i] = item
>>>         i += 1
>>>
>>>     print("initializing")
>>>     # Creating an instance of sensor
>>>     my_shimmer = Shimmer3Streaming(name="Shimmer3_sensor",
>>>                                    output_path="./output")
>>>
>>>     # Creating an instance of device coordinator
>>>     device_coordinator = DeviceCoordinator()
>>>
>>>     # Adding sensor to device coordinator
>>>     device_coordinator.add_devices([my_shimmer])
>>>
>>>     experiment_id = "p01"
>>>
>>>     # A delay to be sure initialing devices have finished
>>>     time.sleep(3)
>>>
>>>     input("\nPress a key to run the scenario")
>>>
>>>     for stimuli_id, stmulus_name in stimuli.items():
>>>         # Starts data recording by displaying the image
>>>         device_coordinator.dispatch(start_message(experiment_id, stimuli_id))
>>>
>>>         # Displaying an image may start with some milliseconds delay after data recording because of GTK
>>>         # initialization in show_image_standalone. If this delay is important to you, use other tools for displaying image stimuli
>>>         # Since image is displaying in another thread we have to manually create the same delay in current
>>>         # thread to record data for 10 seconds
>>>         stimulus = ImageStimulus(stimuli_id, os.path.join(stimuli_path, stmulus_name), 5)
>>>         stimulus.show_standalone()
>>>         time.sleep(5)
>>>
>>>         # Stops data recording by closing image
>>>         device_coordinator.dispatch(stop_message(experiment_id, stimuli_id))
>>>         input("\nPress a key to continue")
>>>
>>>     # Terminate, This step is necessary to close the connection with added devices
>>>     device_coordinator.terminate()


Since the default saving mode is continuous, Shimmer3 will record all data in one file.
For each stimulus, the device records two triggers in the file, one for the start of the stimulus and one for the end of the stimulus.


5- Utilities for designing experiments
--------------------------------------
Octopus Sensing provides some utilities using `GTK <https://athenajc.gitbooks.io/python-gtk-3-api/content/>`_ for
designing a questionnaire, displaying images, and some widgets like creating a timer. We used all of these utilities in
the `full_scenario example <https://github.com/octopus-sensing/octopus-sensing/blob/master/examples/full_scenario>`_. Look at this example to find a simple scenario by
displaying a fixation cross image, displaying a video clip and data recording, and then creating and showing a questionnaire
after each stimulus.
Also, go to the API section and look at the :ref:`questionnaire` and :ref:`windows` documentation to know more about utilities.

6- Monitoring
--------------
See :ref:`octopus_sensing_monitoring` to know more about monitoring and how to use it.
See the example in `full_scenario example <https://github.com/octopus-sensing/octopus-sensing/blob/master/examples/full_scenario>`_ as an example to know more about how to monitor data.

7- Reading recorded data in real-time
---------------------------------------

You can read the data that Octopus Sensing is recording, in real-time, through an HTTP endpoint. To
do so, you can use the same endpoint that Monitoring is using: `MonitoringEndpoint`.

To do so, start the Monitoring Endpoint in the usual way:

>>> from octopus_sensing.device_coordinator import DeviceCoordinator
>>> from octopus_sensing.monitoring_endpoint import MonitoringEndpoint
>>> # Create coordinator instance
>>> coordinator = DeviceCoordinator()
>>> # Add your devices
>>> ...
>>> # Creating the endpoint instance and start it.
>>> monitoring_endpoint = MonitoringEndpoint(coordinator)
>>> monitoring_endpoint.start()
>>> ...

On the client-side (a separate application), simply send a GET request:

>>> import json
>>> import http.client
>>> http_client = http.client.HTTPConnection("127.0.0.1:9330", timeout=3)
>>> http_client.request("GET", "/",
...                     headers={"Accept": "application/json"})
>>> response = http_client.getresponse()
>>> assert response.status == 200
>>> recorded_data = json.loads(response.read())


8- Preprocess and visualize data offline
----------------------------------------

If you used continuous `saving_mode` and want to split them into several files for processing,
Octopus Sensing provides this feature by adding only one line to the end of the previous example.

>>> from octopus_sensing.preprocessing.preprocess_devices import preprocess_devices
>>> preprocess_devices(device_coordinator,
...                    output_path,
...                    shimmer3_sampling_rate=128,
...                    signal_preprocess=True)

By passing the instance of `DeviceCoordinator` as a parameter to `preprocess_devices` function,
it will apply the preprocessing step on all added devices that implemented preprocessing.
For audio and video, we don't need any general preparation.
But, the OpenBCI and Shimmer3 sensor will apply three or two preprocessing steps according to the passed parameters.
It will resample the recorded data for Shimmer3 in this example to a sampling rate of 128 Hz.
Then it will split data based on start and stop triggers.
Then, since `signal_preprocess` is True, it will apply bandpass filtering and cleaning noises.
Finally, this data will be recorded in the specified output path and ready to be used for analysis.

See :ref:`octopus_sensing_visualizer` to know more about visualizer and how to use it.

9- Watching video scenario
---------------------------

Octopus Sensing provides the common scenario in emotion recognition studies. 
In this scenario, the data is recorded during a watching video task, and the user can report emotions using a questionnaire.
Every step in the code is fully commented. By reading and running this example, you can learn how to
do every step in the scenario, monitor data in real-time, and visualize data after finishing the scenario.
See the example in `full_scenario example <https://github.com/octopus-sensing/octopus-sensing/blob/master/examples/full_scenario>`_.    

.. image:: _static/octopus-sensing-logo-small.png
   :alt: Octopus Sensing Logo

Octopus Sensing
====================================================================================

**Welcome to Octopus Sensing's documentation!**

Octopus Sensing is a tool to help you run scientific experiments that involve recording data synchronously from
multiple sources in human-computer interaction studies. You write steps of an experiment scenario, for example showing a stimulus and then
a questionnaire. The tool takes care of the rest.

It can collect data from multiple devices such as OpenBCI EEG headset, Shimmer sensor (GSR and PPG),
Video and Audio and so forth simultaneously. Data collection can be started and stopped synchronously across all devices.
Collected data will be tagged with the timestamp of the start and stop of the experiment, the ID of
the experiment, etc.

The aim is to make the scripting interface so simple that people with minimum or no software
development skills can define experiment scenarios with no effort. 
Also, this tool can be used as the base structure for creating real-time data processing systems like systems with capabilities of recognizing emotions, stress, cognitive load, or analyzing human behaviors.

Main features
~~~~~~~~~~~~~

* Controls data recording from multiple sources using a simple unified interface
* Tags an event on collected data, such as the start of an experiment, and events during the experiment, etc.
* Can show stimuli (images and videos) and questionnaires
* Monitoring interface that visualizes collected data in real-time
* Offline visualization of data from multiple sources simultaneously 

.. image:: _static/OCS-diagram.png
   :alt: Octopus Sensing Components

Quick Start
=========== 
See :ref:`quick_start` to learn how to install Octopus Sensing.

Tutorial
========
See :ref:`tutorial` to learn how to use Octopus Sensing.

API Reference
=============
See :ref:`api_reference` for the complete documentation of all the public classes and methods available in Octopus Sensing.

Contribution Guideline
=======================
See :ref:`contribution` to learn how to contribute to the Octopus Sensing.

Support
========
Please contact Nastaran (nsaffar@gmail.com) or Aidin (aidin@aidinhut.com) if you have any questions or seek support.

Copyright
=========
Copyright © 2020-2022 Nastaran Saffaryazdi, Aidin Gharibnavaz

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

See `License file <https://github.com/octopus-sensing/octopus-sensing/blob/master/LICENSE>`_ for full terms.



Content
========

.. toctree::
   :maxdepth: 3
   :caption: .

   quick_start
   tutorial
   api_reference
   contribution
   development
   visualizer
   monitoring.. _octopus_sensing_visualizer:

***************************
Octopus Sensing Visualizer
***************************

Octopus Sensing Visualizer is a web-based real-time visualizer for `Octopus Sensing <https://octopus-sensing.nastaran-saffar.me/>`_.
It can be used for offline data visualization. You can visualize the recorded multimodal data as the raw data. Also, it can extract
some essential features or components of data and display them in a single window. Using this tool, you can observe the effect of an event on recorded data simultaneously.

`Octopus Sensing Visualizer <https://github.com/octopus-sensing/octopus-sensing-visualizer>`_ is
a separated project and can be installed if we need to visualize data.
It can be used for displaying recorded data with
the same format as we recorded through Octopus Sensing.

**Note**

If we want to display data that is recorded by Octopus Sensing,
we should apply :ref: `api_reference/preprocessing` module to prepare data for the visualizer while we record data or later.



Installation
------------
It requires Python 3.7 or later.

You can use `pip` to install it:

.. code-block:: bash

   $ pip install octopus-sensing-visualizer

You can also use one of the Python package managers like `pipenv <https://pipenv.pypa.io/en/latest/>`_ or
`poetry <https://python-poetry.org/>`_ to prevent package conflict.

.. code-block:: bash

   $ pipenv install octopus-sensing-visualizer

How to use it
--------------
At first, you should create an `octopus_sensing_visualizer_config.conf` in the current directory.
This config file includes the path to the data and the type of graphs that we want to visualize.
See :ref:`visualizer_config_guide` to know how to prepare this file.

Then simply run the server by invoking `octopus-sensing-visualizer` from the command line.
For example, if you use `pipenv <https://pipenv.pypa.io/en/latest/>`_ as the package manager, run it as follows:

.. code-block:: bash

   $ pipenv run octopus-sensing-visualizer


The visualizer will listen on `8080` port. Open a web page and point to the machine's IP. For
example, in the same machine, open http://localhost:8080 . Or replace `localhost` with the machine's
IP and open it from any other machine.


Security notice
---------------
Note that the webserver accepts requests from any machine, and it uses `http` protocol, which
is not encrypted. Don't run it on a network that you don't trust.


How to prepare a config file
----------------------------
.. toctree::
   :maxdepth: 1

   visualizer_config_guide

Which data can be visualized
-----------------------------
It can be used for displaying recorded data using any software if the data is prepared with
the same format as we recorded through Octopus Sensing.

**EEG** :
A CSV file with 16 columns for 16 channels. The number of rows shows the number of samples.

**GSR**:
A CSV file with one column. The number of rows shows the number of samples.

**PPG**:
A CSV file with one column. The number of rows shows the number of samples.


User interface
--------------
The user interface includes graphs with True value in the config file.
Also, it has a slide bar that allows us to go forward and backward in data samples.
There is a text box for setting window size. This window size shows the length of the horizontal axis in time.
Setting the bigger values leads to displaying a wider window of data in each moment and so fewer details.
The sliding bar is moving each second, and all graphs will be updated in each second.

**Copyright**

Copyright © 2021 [Nastaran Saffaryazdi]

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

See `License file <https://github.com/nastaran62/octopus-sensing/blob/master/LICENSE>`_  for full terms.
.. _quick_start:

***********
Quick start
***********

requirements
============

You need `Python <https://python.org>`_ installed on your computer (version 3.7 or higher). Refer to
`this guide <https://realpython.com/installing-python/>`_ if you need help.

Quickstart Using init script (Linux & Mac)
==========================================

Octopus Sensing comes with a script that helps you quickly start a project. It uses
`Pipenv <https://pipenv.pypa.io/>`_ to create a `virtual
environment <https://docs.python.org/3/tutorial/venv.html>`_ in order to keep everything clean. It
will also create a sample application.


.. code-block:: bash

    mkdir my-awesome-project
    cd my-awesome-project
    curl --output init.sh https://raw.githubusercontent.com/nastaran62/octopus-sensing/master/init_script/init.sh
    # It's a good idea to read any script before executing it.
    bash ./init.sh
    rm ./init.sh


The created `main.py` file is an example application. To run it:

:code:`pipenv run python main.py`


If you don't want to use the script, you can use the following methods instead.

Installation using Pipenv (All Platforms)
=========================================

We recommend using a package manager like `Pipenv <https://pipenv.pypa.io/>`_ instead of globally
installing Octopus Sensing using `pip` to prevent package conflicts. To do so, follow these
commands. (This is same as what the above script does.)


.. code-block:: bash

    mkdir my-awesome-project
    cd my-awesome-project
    # Or replace it with your python version
    pipenv --python python3.8
    pipenv install octopus-sensing



It installs Octopus Sensing inside the virtual environment created by Pipenv. You need to use
`pipenv` to run your code. For example:

:code:`pipenv run python main.py`


Refer to `Pipenv website <https://pipenv.pypa.io/>`_ for more info.

Installation using pip (All Platforms)
======================================

You can use `pip` to install `octopus-sensing` as simple as:

:code:`pip3 install octopus-sensing`

(You might need to replace `pip3` with `pip` depending on your system.)

Then it can be imported like:

:code:`import octopus_sensing`


Installation from source (All Platforms)
========================================

If you want to compile it from source for development purposes or to have the un-released features,
please refer to :ref:`development`.

Troubleshooting
===============

- Pip cannot install PyGObject on Windows. If users want to use `octopus-sensing.stimuli` or `octopus-sensing.windows` packages, they need to install it manually themselves. See `PyGObject documentation <https://pygobject.readthedocs.io/en/latest/getting_started.html#windows-getting-started>`_ to know how to install PyGObject on Windows.
.. _api_reference:

*************
API Reference
*************

octopus_sensing
=============================

.. toctree::
   :maxdepth: 4

   api/octopus_sensing.device_coordinator
   api/octopus_sensing.device_message_endpoint
   api/octopus_sensing.devices
   api/octopus_sensing.common
   api/octopus_sensing.stimuli
   api/octopus_sensing.questionnaire
   api/octopus_sensing.windows
   api/octopus_sensing.preprocessing
   


.. _development:

***********
Development
***********

Installing from source
======================

you need `poetry`:

>>> pip3 install poetry

(Refer to `Poetry installation guide <https://python-poetry.org/docs/#installation>`_
for alternative ways of installing Poetry.)

Then `cd` to where the source code located and run:

>>> poetry install
>>> poetry build

It will create a virtual environment and installs `Octopus Sensing` with its dependencies in it.

Coding Style
==============

We're following Python PEP 8 for our coding style.

For formatting the code, `autopep8 <https://github.com/hhatto/autopep8>`_ is a good tool.
Many editors support it. You can use it to automate the code formatting.

In-code Comments
~~~~~~~~~~~~~~~~~~~
Add comments to clarify *why* you did things this way. Usually, it's easy to figure out *what* a piece
of code does, but *why* is harder or impossible to figure out.

Also, add a comment for complex algorithms even if they are written very clearly.

Doc Strings
~~~~~~~~~~~~
Every public method or function should have doc string. We also generate our API Reference document
from these doc strings. So ensure they are clear and addresses all the functionality and exceptions
of a method. We are using `NumPy style <https://numpydoc.readthedocs.io/>` for creating doc strings.

Static Type Checking
======================

We use `mypy <http://www.mypy-lang.org/>`_ to type-check our sources. Every variable or parameter
in the source code should have a type, unless MyPy can automatically determine the type.

To run MyPy use our make file:

.. code-block:: bash

   $ make mypy


Tests
======
We're using `pytest <https://docs.pytest.org>`_ to run the tests. You can simply invoke the tests using:

.. code-block:: bash

   $ make test

There are two sets of tests: An integration test that checks the overall health of the library by running
a full scenario, and a lot of small unit tests that are testing functionalities individually.

All tests are located in `octopus-sensing/tests` directory.

Adding Support for a New Device
===============================

To add support for a new device, you need to create a subclass of `octopus-sensing.devices.device.Device`.

The device's code will be run in a separate process. Because Python is not good with threading, and also
because it minimizes the effect of other parts of the application on the data collection.

The Device process and the parent process (the Device Coordinator) talk with each other using Message Queues.

When implementing the `Device` class, you need to override `__init__` and `_run` methods.

In `__init__`, you will receive the required parameters. For example, if the device needs configuration options.
You also need to receive the same parameters as in the `Device` class (your base class) and pass them to your
base class using `super()`.

The `_run` method will be run in the separate process. You need to initialize the device here, and start
recording the data. At the same time, you need to check the messages in `self.message_queue`. So, usually,
you need to do your data recording in a separate thread, and check the messages in the main thread.

The following code can be used as a starting point for adding a device. And also have a look at the devices
currently implemented in Octopus Sensing for some sample codes.

.. code-block:: python

    import threading

    from octopus_sensing.common.message_creators import MessageType
    from octopus_sensing.devices.device import Device

    # We inherit from Device class
    class SampleDevice(Device):

        # 'name' and 'output_path' are from our parent, the Device class
        def __init__(self, config_flag, output_path, name=None):
            # Parent should always be called. It does some initialization of itself.
            # We're passing the parameters we received to it.
            super().__init__(name=name, output_path=output_path)

            # Keeping the config parameter
            self._config_flag = config_flag

            # Note that we don't do anything with the device here.
            # Everything should be done after the process is created,
            # in the _run method.

        # Note that this is '_run' and not 'run'!
        # You should never override 'run'.
        def _run(self):
            # Initialize your device here.
            self._device_handle = ...
            # Then we start a thread for recording the data.
            # We will use this flag to tell the thread to finish recording.
            self._record = True
            threading.Thread(target=self._record_data).start()

            # We're checking messages in the main thread.
            while True:
                # This will block until a message receives from the parent (the deivce coordinator)
                message = self.message_queue.get()
                if message.type == MessageType.TERMINATE:
                    # This will cause the recording thread to exit. (see its code)
                    self._record = False
                    # Exiting the main loop. It will cause the process to finish and terminate.
                    # (since there's nothing after this.)
                    break


        def _record_data(self):
            # This is running in another thread (see _run)
            # Do the actual data recording here.
            while self._record:
                data = self._device_handle.read()
                # Write it to a file for example.

            # Depending on the device, you might want to start recording data
            # when you received the START message in the message_queue, and
            # stop recording when you received the STOP message.
******************
Device Coordinator
******************

.. automodule:: octopus_sensing.device_coordinator
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: get_monitoring_data, MonitoringCache
.. _preprocessing:

Preprocessing
======================================


OpenBCI
---------------------------------------------

.. automodule:: octopus_sensing.preprocessing.openbci
   :members:
   :undoc-members:
   :show-inheritance:

Shimmer3
----------------------------------------------

.. automodule:: octopus_sensing.preprocessing.shimmer3
   :members:
   :undoc-members:
   :show-inheritance:


Preprocess Devices
---------------------------------------------------------

.. automodule:: octopus_sensing.preprocessing.preprocess_devices
   :members:
   :undoc-members:
   :show-inheritance:
   :Exclude-members:


.. _stimuli:

Stimuli
================================

.. automodule:: octopus_sensing.stimuli
   :members:
   :undoc-members:
   :show-inheritance:

Image Stimulus
--------------------------------------------

.. automodule:: octopus_sensing.stimuli.image_stimulus
   :members:
   :undoc-members:
   :show-inheritance:


Video Stimulus
----------------------------------------------

.. automodule:: octopus_sensing.stimuli.video_stimulus
   :members:
   :undoc-members:
   :show-inheritance:
.. _questionnaire:

Questionnaire
======================================

Questionnaire
--------------------------------------------

.. automodule:: octopus_sensing.questionnaire.questionnaire
   :members:
   :undoc-members:
   :show-inheritance:


Base question
----------------------------------------------

.. automodule:: octopus_sensing.questionnaire.question
   :members:
   :undoc-members:
   :show-inheritance:

Text Question
----------------------------------------------------

.. automodule:: octopus_sensing.questionnaire.text_question
   :members:
   :undoc-members:
   :show-inheritance:


Opinion question
-------------------------------------------------------

.. automodule:: octopus_sensing.questionnaire.opinion_question
   :members:
   :undoc-members:
   :show-inheritance:




Common
=======


Message
--------

.. automodule:: octopus_sensing.common.message
   :members:
   :undoc-members:
   :show-inheritance:

Predefined Messages
-------------------

.. automodule:: octopus_sensing.common.message_creators
   :members:
   :undoc-members:
   :show-inheritance:

HTTP Endpoint
-------------------------------------------------

.. automodule:: octopus_sensing.device_message_endpoint
   :members:
   :undoc-members:
   :show-inheritance:.. _windows:

Windows
================================
Creates some common dialogs/windows using Gtk, which can be used in designing Experiment scenarios


Message Window
-----------------------------------------------

.. automodule:: octopus_sensing.windows.message_window
   :members:
   :undoc-members:
   :show-inheritance:


Image Window
---------------------------------------------

.. automodule:: octopus_sensing.windows.image_window
   :members:
   :undoc-members:
   :show-inheritance:


Timer Window
---------------------------------------------

.. automodule:: octopus_sensing.windows.timer_window
   :members:
   :undoc-members:
   :show-inheritance:

.. _devices:

*******
Devices
*******

Device (Base class)
-------------------

.. automodule:: octopus_sensing.devices.device
   :members:
   :undoc-members:
   :show-inheritance:


MonitoredDevice
---------------

.. automodule:: octopus_sensing.devices.monitored_device
   :members:
   :undoc-members:
   :show-inheritance:

Camera
------

.. automodule:: octopus_sensing.devices.camera_streaming
   :members:
   :undoc-members:
   :show-inheritance:


Audio
------

.. automodule:: octopus_sensing.devices.audio_streaming
   :members:
   :undoc-members:
   :show-inheritance:


Shimmer3
--------

.. automodule:: octopus_sensing.devices.shimmer3_streaming
   :members:
   :undoc-members:
   :show-inheritance:


OpenBCI (brainflow)
-------------------

.. automodule:: octopus_sensing.devices.brainflow_openbci_streaming
   :members:
   :undoc-members:
   :show-inheritance:



OpenBCI (pyOpenBCI)
--------------------

.. automodule:: octopus_sensing.devices.openbci_streaming
   :members:
   :undoc-members:
   :show-inheritance:



Brainflow
---------

.. automodule:: octopus_sensing.devices.brainflow_streaming
   :members:
   :undoc-members:
   :show-inheritance:



Open Vibe
---------

.. automodule:: octopus_sensing.devices.open_vibe_streaming
   :members:
   :undoc-members:
   :show-inheritance:



Network
--------

.. automodule:: octopus_sensing.devices.network_devices.socket_device
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: octopus_sensing.devices.network_devices.http_device.HttpNetworkDevice
   :members:
   :undoc-members:
   :show-inheritance:



Common
-------

.. automodule:: octopus_sensing.devices.common
   :members:
   :undoc-members:
   :show-inheritance:
