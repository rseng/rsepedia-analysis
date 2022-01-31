PiVR has been developed by David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

.. _FAQ:

FAQ
===================

.. _FAQ_undistort:

Distorted Images
-----------------

*Question:*
  Why do my images/videos look distorted?

*Answer:*
  Every lens will introduce radial distortions to the image.
  Since the Raspbbery Pi Camera lens is not a high quality lens the
  radial distortion can become very obvious. The distortion is a
  function of the lens you are using meaning the distortion is
  identical for all videos/images taken by the same camera.

  To get remove camera distortions in videos, :ref:`please see
  here<undistort_h264_GUI>`

Progress Bar
-------------

*Question:*
  When running a tracking/virtual reality/video experiment on PiVR, why
  is there no progress bar?

*Answer:*
  PiVR was designed to process each frame as quickly as it can. This is
  necessary to produce realistic virtual realities. Having a progress
  bar undermines this goal as the act of updating the progress bar
  increase latency.

  @ Developers: If you think this problem can be solved, please
  `let us know <https://gitlab.com/LouisLab/pivr/issues>`_ - if there
  is a way to update a progress bar in less than a
  millisecond, this could and should be implemented.

Raspberry Pi does not turn on
-----------------------------

*Question:*
  I am unable to start the Raspberry Pi. It seems to start, but it
  never progresses to the Desktop.

*Answer:*
  We have observed this when the SD card was completely full! This
  can happen if you run experiments until you encounter an error
  saying that the experimental data can not be saved. If you then
  turn off the Raspberry Pi, the OS is unable to start.

  We have been able to recover the data by removing the SD card from
  the PiVR setup and putting into a Ubunut workstation card reader.
  The files on the SD card can then be rescued. After removing the
  files, the Raspberry Pi was able to start again.

  Alternatively, you can also just format the SD card and re-install
  Raspbian and PiVR, of course.

  To avoid this error in the future, please always remove data before
  the SD card fills up!

PiVR won't start a video/experiment
-----------------------------------

*Question:*
  I want to record a video/start and experiment, but nothing happens.

*Answer:*
  There can be several causes for this. Usually, the terminal will be
  very helpful in troubleshooting why you are having a particular
  problem. Here, we list common problems and how to fix them. If you
  do not find your error message listed here, please take a
  picture/write it down and open an `issue <https://gitlab
  .com/LouisLab/pivr/issues>`_

    a) *PermissionError: [Errno 13] Permission denied: ...*
        Please make sure to select a folder where the experiment
        should be saved. That is the button next to "Save in:"
    b) *The camera is already using port %d ' % splitter_port) ...*
        We have seen this error after incorrectly closing the PiVR
        software. Restart the Raspberry Pi to fix this.

PiVR software won't start
--------------------------

*Question:*
  I can't start the PiVR software on the Raspberry Pi

*Answer:*
  There could be several causes for this. To get to the bottom of the
  this bug please do the following:

    a) Restart the Raspberry Pi.
    b) Doubleclick on the desktop shortcut as if starting the software
       (this step is necessary!).
    c) Open the terminal.
    d) Change directory to the PiVR folder and start PiVR manually.
       If you installed PiVR with the provided script you would type::

          cd PiVR/PiVR
          python3 start_GUI.py

    e) See below if your error is explained below. If not, please take a
       picture/write it down and open an
       `issue <https://gitlab.com/LouisLab/pivr/issues>`_.

1) Camera error

.. figure:: Figures/FAQ/1_cam_not_enabled.png
   :alt: 1_cam_not_enabled.png

This could be due to several issues:

  - The camera might not have been enabled. The provided installation
    script normally enables the camera but it worth
    `double-checking <https://projects.raspberrypi.org/en/projects/getting-started-with-picamera/2>`_.

  - If this does not solve the problem, double-check whether the cable
    connecting the Raspberry Pi and the camera is plugged in with the
    correct orientation.

  - If this does not solve the problem it is likely that either the
    camera cable or the camera itself is faulty. It is more likely that
    the cable is faulty as we have experienced this more than once
    already. The solution is to buy a new camera cable and test it.
    If this does not fix the problem, try replacing the camera.

    Do note that we have seen batches of camera cables that seem to have
    a high (~50%) failure rates. We recommend using reputable sellers
    such as Mouser or DigiKey or better specialist sellers such as
    adafruit.

I have camera glitches/interference
-----------------------------------

*Question:*
  I sometimes have lines appearing in my image.

   .. figure:: Figures/FAQ/2_camera_interference.jpg
      :alt: 2_camera_interference.jpg

*Answer:*
  It seems that this happens due to interference of the signal in the
  `camera cable <https://forums.raspberrypi.com/viewtopic.php?t=300009>`_.

  Make sure the camera cable is a short as possible. Ideally, the camera
  cable should be completely straight. If that doesn't help you could
  try to shield the cable or use a
  `completely different solution <https://www.thinesolutions.com/cable-extension-kit>`_.
  Note that we haven *not* tested this but it claims to solve the
  problem.PiVR has been developed by David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

Community Building
******************

.. _Labs_using_PiVR:

Labs using PiVR
================

- Carrillo Lab (https://www.carrillolab.com/), University of Chicago, United States

- Louis Lab (https://louis.mcdb.ucsb.edu/), University of California, Santa Barbara, United States

- Schilder Lab (http://www.personal.psu.edu/rjs360/index.html), Pennsylvania State University, United States

- Sprecher Lab (http://www.sprecherlab.com), University of Fribourg, Switzerland

.. Note::

    If you are using PiVR, please :ref:`contact us <contact>` so that we
    can list you here.

.. _PiVR_in_press:

PiVR in the press
==================

- `Scientific American <https://www.scientificamerican.com/article/fruit-flies-plug-into-the-matrix/>`__

- `Eureka Alert <https://www.eurekalert.org/pub_releases/2020-07/p-arp070720.php>`__

- `Laboratory News <https://www.labnews.co.uk/article/2030709/serving-virtual-reality-raspberry-pi-to-flies>`__

- `Nano Werk <https://www.nanowerk.com/news2/biotech/newsid=55649.php>`__

- `Science Daily <https://www.sciencedaily.com/releases/2020/07/200714143044.htm>`__

- `SciTechDaily <https://scitechdaily.com/virtual-reality-system-for-small-animals-based-on-raspberry-pi/>`__

- `UCSB News <https://www.news.ucsb.edu/2020/019955/affordable-alternative>`__

PiVR on blogs
==============

- `LabOnTheCheap.com <https://www.labonthecheap.com/pivr-415-to-create-a-vr-platform-for-small-creatures/>`__

- `Open-Neuroscience.com <https://open-neuroscience.com/post/pivr/>`__
PiVR has been developed by David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

.. _Benchmark:

PiVR Benchmark
***************

.. _Concepts:

Concepts
========

When creating a virtual reality two different measurements of time
are relevant:

#. The **loop time**: How often per second does the system take input
   from the subject? This is usually expressed in **Hz**.
#. The **latency**: How long between taking the input and presenting the
   subject with the updated stimulus that makes the virtual reality?
   This is usually expressed in milliseconds (**ms**).

Next, how can one measure time on a computer? Most operating systems
(Windows, MacOS and Linux) are *non-real time operating systems*.
This means that the OS can **not** guarantee that a given line of
code is executed at a fixed interval. If one asks the system to
simply record the time to measure *loop time* and *latency* on
introduces uncertainty, possibly in the order of milliseconds.

Luckily PiVR has a real-time operating system - on the GPU:
Digital cameras work by reading out lines on their sensor and the
camera used by PiVR is no different. (Read `here
<https://picamera.readthedocs.io/en/release-1.13/fov.html#>`__ for a
fantastic introduction).

For example at 90 frames per second at 640x480 the camera must read
out 480 lines in ~11ms. It must therefore spend no more than 2.3^-5s
or 23.1us per line. This means every 23.1us the camera is pushing one
frame line to the GPU (see `here
<https://picamera.readthedocs.io/en/release-1.13/fov.html#division-of-labor>`__
for detailed information).

The GPU records the timestamp of when the first line of a new frame
is received based on its real-time hardware clock. We
take advantage of this while doing the PiVR measurements. This timestamp
is called presentation timestamps (**PTS**). (More info `here
<https://picamera.readthedocs.io/en/release-1.13/_modules/picamera/frames.html?highlight=timestamp##>`__).

.. _Benchmark_loop_time:

PiVR loop time
==============

PiVR loop time has been tested by recording the PTS while doing
real-time tracking at 640x480 resolution and a variety of framerates
using a Raspberry Pi 3.

.. figure:: Figures/benchmarking/1_loop_time_640x480.jpg
  :width: 100 %
  :alt: 1_loop_time_640x480.jpg

For example at 30fps each camera timestamp should be 1/30 = 0.033s
which of course is 33ms.
The reported timestamps are perfect for 30 and 40fps. At 50, 60 and
70fps we start loosing a few frames (6/5,598, 25/7,198 and 16/8,398,
respectively) meaning PiVR could not keep up with the images coming
in therefore the loop time was larger than 1/fps.

For more detail please see `Timing Performance of PiVR
<https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000712#sec025>`__

.. _Benchmark_latency:

PiVR latency
=============

To measure latency we turned the LED on at a defined frame and asked
how long it takes for PiVR to detect the change in light intensity.

This was based on the excellent suggestion of `Reviewer #1
<https://journals.plos.org/plosbiology/article/peerReview?id=10.1371/journal.pbio.3000712>`__:

#. Run closed loop tracking at 70fps at 640x480 with a facsimile larva.
#. At each 50th frame of a trial (and multiples thereof), turn the red
   LED on.
#. Take the median of a set of 10 pixels in the search box and compare
   the average current intensity to the intensity of the previous frame.
#. If this difference is larger than 2, consider the LED “ON” and
   immediately turn the LED “OFF”.
#. If the LED has been turned “OFF”in the previous frame, turn the LED
   “ON”again.

The code used for this can be found on the Gitlab repository of PiVR
`Gitlab repository of PiVR <https://gitlab.com/LouisLab/pivr>`__ in
branch ‘LED_flash_test’ in the file ‘fast_tracking.py’ with
`relevant code going from line 1457 to 1584
<https://gitlab.com/LouisLab/pivr/-/blob/LED_flash_test/PiVR/fast_tracking.py>`__.

Using the methodology described above, we found that the detection of a
‘LED ON’ event was either detected during the earliest possible frame or
one frame later. This difference in latency was dependent on the region
of the image (region of interest, ROI) used to compare the pixel
intensity. The detection latency was systematically longer when the
ROI was in the top (green box) as compared to the bottom of the image
(magenta box).

.. figure:: Figures/benchmarking/2_latency_70fps_640x480.jpg
  :width: 100 %
  :alt: 2_latency_70fps_640x480.jpg


As illustrated in the figure below, this can be explained by the fact
that turning the LED “ON” while the camera is grabbing the frame through
the rolling shutter cannot lead to a successful detection of the
change in light intensity. During frame #50 the LED is still “OFF” in
the top part of the image whereas it is “ON” in the lower part of
the image.

.. figure:: Figures/benchmarking/3_rolling_shutter_effect.png
  :width: 100 %
  :alt: 3_rolling_shutter_effect.png

As shown above the number of frames between two consecutive events
corresponding to the LED being turned “ON” was 3 frames if the ROI was
located in the top part of the image (green box). The schematic
diagram below outlines our interpretation of this 3-frame
delay: The LED is being turned “ON” during frame (j) while the frame is
being recorded (vertical arrow in frame j + 1). If the PiVR algorithm
is set to monitor a ROI in the upper part of the image, this ROI is
read before the LED intensity has switched to the “ON” state. Therefore,
the LED “ON”event is not observed during frame j. Instead, it is
detected during processing of frame j + 1. This delay happens again
when the second LED “ON” event is detected in frame j + 4.

.. figure:: Figures/benchmarking/4_top_diagram.jpg
  :width: 100 %
  :alt: 4_top_diagram.jpg

By contrast, if PiVR is set to monitor a ROI in the lower part of the
image, the LED turn “ON” event can be detected while processing frame j
+ 1
(see below).

.. figure:: Figures/benchmarking/5_bottom_diagram.jpg
  :width: 100 %
  :alt: 5_bottom_diagram.jpg

Consequently, the time that elapses between the two LED “ON” events is
either 5 frames (or 71ms at 70 Hz) for a ROI in the top part of the
image or 2 frames for a ROI in the bottom part of the image.

Based on these observations, we conclude that the LED “ON” event must
be happening **while** the image is being captured, which corresponds
to a duration of **1.34ms**. The longest possible latency between the
action of an animal (in this example during frame j) and the actuation
of the update of the LED is therefore 2 frames + 1.34ms
(time window associated with frame j + 2), which is equivalent to less
than 30 ms (frame duration at 70 Hz is **14.29 ms**.

.. _PiVR_loop_time_high_res:

PiVR loop time (High Res)
==========================

PiVR v1.5.0  was released on 27th of March, 2021. It allows users to
perform :ref:`online tracking <TrackingLabel>` with resolutions other
than 640x480.

Higher resolutions impact the :ref:`closed loop times <Concepts>` as
a larger image is being analyzed to define the animal in each frame.

We tested the performance using both a
`Raspberry Pi 3 Model B Plus Rev 1.3
<https://www.raspberrypi.org/products/raspberry-pi-3-model-b-plus/>`__
and the newer and more powerful
`Raspberry Pi 4 Model B Rev 1.2
<https://www.raspberrypi.org/products/raspberry-pi-4-model-b/>`__

The **Raspberry Pi 3** was able to handle 30fps without dropping frames
at 1024x768. At 1296x972 it was only able to handle 15fps without
dropping frames.

.. figure:: Figures/benchmarking/6_RPi3_Image_processing_time.jpg
  :width: 100 %
  :alt: 6_RPi3_Image_processing_time.jpg

The `50% more powerful
<https://www.raspberrypi.org/documentation/hardware/raspberrypi/bcm2711/README.md>`__
**Raspberry Pi 4** on the other hand was able to handle 40fps at
1024x768 and 30fps at 1296x972:

.. figure:: Figures/benchmarking/7_RPi4_Image_processing_time.jpg
  :width: 100 %
  :alt: 7_RPi4_Image_processing_time.jpg
PiVR has been developed by David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

.. _contact:

Contact
===================

Questions?
----------

If something does not work as you expect, one of the best ways to get
information is to open an "issue" on the `PiVR gitlab repository
<https://gitlab.com/LouisLab/PiVR/issues>`__.

Found a Bug?
------------

Please open an "issue" on the `PiVR gitlab repository
<https://gitlab.com/LouisLab/PiVR/issues>`__. It is imperative to
give as much information in when you encountered the bug. Ideally,
you are able to let us know exactly how to reproduce the bug. At
least you need to let us know:

#. The OS (on the Raspberry Pi it is probably Raspbian)
#. What exactly you were doing when the error occured.

Other questions (scientific, technological, copyright...)
------------------------------------------------------------------------

Please contact the corresponding author of this work, `Matthieu Louis
<https://www.mcdb.ucsb.edu/people/faculty/matthieu-louis>`__.
PiVR has been developed by David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

.. _PIVR Hardware Manual:

Build your own PiVR
===================

It doesn't take long and it isn't too hard to build PiVR. See this
timelapse video where a PiVR setup is being 3D printed, the Printed Circuit Board (PCB) soldered
and then built.

.. raw:: html

    <iframe width="75%"
    src="https://www.youtube.com/embed/w5tIG6B6FWo" frameborder="0"
    allow="accelerometer; autoplay; encrypted-media; gyroscope;
    picture-in-picture" allowfullscreen></iframe>

Standard and High Powered version
---------------------------------

PiVR comes in two versions:

#. :ref:`The Standard version<PiVR Standard Version>`:
   LED strips are used for both the backlight and stimulation light.
   We have measured around :math:`2{\mu}W/mm^2` for red light
   intensity and around :math:`8000 Lux` when measuring white light
   intensity. Please see the :ref:`Bill of Materials <BOM>` for
   components you must buy to build the setup.

#. :ref:`The High Powered Version<PiVR HP Version>`: The
   stronger LEDs require not only a different arena but also
   dedicated LED drivers which provide a constant current. We have
   measured light intensities in the order of :math:`40 -
   50{\mu}W/mm^2` for the red LEDs.

.. _PiVR Standard Version:

Building the Standard version
-----------------------------

.. Warning::

    You will be handling sensitive electronic equipment. As you might
    be electrically charged, always touch something metallic that is
    connected to the ground before unpacking anything from an
    electrostatical shielding bag. A radiator will do fine.

.. note::

    You can of course change the order of how to build the different
    parts. I build them in this order as there are necessary delays
    such as installing the operation system on the Raspberry Pi.

.. note::

    We used an Ultimaker 3 with a Ultimaker Print Core AA (0.8mm) to
    print the PLA and a Ultimaker Print Core BB (0.8) to print the
    PVA support material. STL files were converted using Ultimaker
    CURA software.

#. 3D print each part in the folder PiVR/Hardware/3D printer files.
   We used an Ultimaker 3 and printed with Standard PLA. We used PVA as
   support material. For best results print Casing, CameraHolder,
   TowerBottom and TowerTop with support material.

    .. figure:: Figures/manual_hardware/standard/1_3D_printed_parts.jpg
       :width: 100 %
       :alt: 1_3D_printed_parts.jpg

#. Obtain the Printed Circuit Board (PCB). Find the blueprint in
   PiVR\Hardware\PCB. I have used the software
   `Fritzing <http://fritzing.org/home/>`_ to design the board. I
   have used the company  `Aisler.net <https://aisler.net/>`_ to print
   the circuit boards.

    .. figure:: Figures/manual_hardware/standard/2_0_PCB.jpg
       :width: 100 %
       :alt: 2_0_PCB.jpg

#. Get the following items to solder the PCB board:

    .. important::

        Before touching electrical components always touch something
        metallic that is connected to the ground. For example, a radiator.

    #. Four (4) Transistors: 30N06L
    #. One (1) GPIO Header
    #. One (1) Barrel connector, 5.5mm Jack with 2.1mm center pole
       diameter
    #. Four (4) Barrel connectors, 3.5mm Jack with 1.35mm center pole
    #. Four (4) :math:`10k{\Omega}` resistors
    #. Break off a 5 pin stretch from the Breakaway Headers. This can
       be done using scissors or pliers.

   .. figure:: Figures/manual_hardware/standard/2_7_PCB_components.jpg
      :width: 100 %
      :alt: 2_7_PCB_components.jpg

#. Solder the components on the PCB. See
   :ref:`this section<detailed_soldering_instructions>` for
   detailed soldering instructions.

    .. important::

        The correct **orientation** of the **GPIO header** and the
        **Transistors** is crucial for PiVR to work correctly.

   .. figure:: Figures/manual_hardware/standard/3_solderedPCB.jpg
      :width: 100 %
      :alt: 3_solderedPCB.jpg

#. Cut off the excess wire of the resistors and the transistors, e.g.
   with scissors.

   .. figure:: Figures/manual_hardware/standard/4_PCB_removed_wire.jpg
      :width: 100 %
      :alt: 4_PCB_removed_wire.jpg

#. Unpack the Touchscreen, remove the stand-off screws. Attach the
   monitor cable that came with the Touchscreen and the 4"
   5 pin cable.

    .. important::

        The monitor cable must be inserted in the correct orientation.
        When you look into the receptacle you'll see that only one
        side has connectors. Make sure that you insert the cable's
        connectors on the same side.

    .. important::

        Note the orientation of the 4" 5pin cable! Left is (+) while
        Right is (-).

   .. figure:: Figures/manual_hardware/standard/6_Touchscreen.jpg
      :width: 100 %
      :alt: 6_Touchscreen.jpg

#. Place the Casing on top of the Touchscreen (it will only fit in
   the shown orientation). Organize the 4" 5pin cable and the
   monitor cable as shown in the picture. Use the M2.5x10mm screws to
   fix the casing to the touchscreen.


   .. figure:: Figures/manual_hardware/standard/7_2_Casing_on_Touchscreen.jpg
      :width: 100 %
      :alt: 7_2_Casing_on_Touchscreen.jpg

#. Prepare the SD card: Format the SD card using SD Formatter and load with
   NOOBs installation files as instructed
   `here <https://www.raspberrypi.org/help/noobs-setup/2/>`_:

#. Connect monitor cable with the Raspberry Pi (with inserted SD
   card). Again, make sure you insert the cable in the correct
   orientation. Use M2.5x10 screws to attach the Raspberry Pi to the
   Casing.

   .. figure:: Figures/manual_hardware/standard/8_0_Raspberry_monitor_cable.jpg
      :width: 100 %
      :alt: 8_0_Raspberry_monitor_cable.jpg

#. Attach the PCB on the right side of the casing using M2.5x10mm
   screws. Plug the 4" 5pin cable into the PCB in the **correct**
   orientation as shown!

   .. figure:: Figures/manual_hardware/standard/10_PCB_RPi_on_casing.jpg
      :width: 100 %
      :alt: 10_PCB_RPi_on_casing.jpg

#. Use the GPIO Ribbon cable to connect the PCB with the
   Raspberry Pi. Thread the long camera cable through the slit as
   shown in the image below. Connect it to the Raspberry Pi camera port.

   .. figure:: Figures/manual_hardware/standard/12_Casing_stuffed.jpg
      :width: 100 %
      :alt: 12_Casing_stuffed.jpg

#. Attach the Backside to the CasingPedestal with M2.5x10 screws.

   .. figure:: Figures/manual_hardware/standard/13_Backside_attached_to_Pedestal.jpg
      :width: 100 %
      :alt: 13_Backside_attached_to_Pedestal.jpg

#. Slide the Backside into the Casing

#. Drop a 2.5mm nut in each hole in the TowerPedestal. Use the M2.5x10
   screws to attach the TowerBottom to the Tower Pedestal

   .. figure:: Figures/manual_hardware/standard/14_TowerPedestal_attached.jpg
      :width: 100 %
      :alt: 14_TowerPedestal_attached.jpg

#. Using a hammer, drive the dowel pins into the TowerBottom. Then,
   attach the TowerTop to it. In principle you can stack more TowerTops
   on top.

   .. figure:: Figures/manual_hardware/standard/15_TowerBottom_Dowel.jpg
      :width: 100 %
      :alt: 15_TowerBottom_Dowel.jpg

#. Attach the 800nm Longpass Filter to the Camera using Parafilm. It
   is best to wear gloves for this step.

   .. Important::

      If you are **not** using the camera or lens in the image below,
      please :ref:`provide<create_own_undistort_files_label>`
      your own undistort files before running experiments.

   .. figure:: Figures/manual_hardware/standard/16_Camera_with_LP_filter.jpg
      :width: 100 %
      :alt: 16_Camera_with_LP_filter.jpg

#. Thread the camera cable from the Casing through the slit in the
   TowerBottom and through the slit of the Camera Holder.

    .. important::

        Note the orientation to avoid having to curl the camera cable
        in the camera holder

   .. figure:: Figures/manual_hardware/standard/17_Threading_of_Cam_cable.jpg
      :width: 100 %
      :alt: 17_Threading_of_Cam_cable.jpg

#. Attach the Camera Cable to the Camera in the **correct**
   orientation. Then, screw the camera to the Camera Holder using the
   M2.5x10 screws. It is **not** necessary to fixate the screws with
   nuts!

   .. figure:: Figures/manual_hardware/standard/18_Camera_attached_to_CamHolder.jpg
      :width: 100 %
      :alt: 18_Camera_attached_to_CamHolder.jpg

#. Drop a 2.5mm nut in the hole in the Camera holder and use it to
   fasten the M2.5x10 screw. Then attach the CameraHolder to the
   Tower.

   .. figure:: Figures/manual_hardware/standard/19_PiVR_complete.jpg
      :width: 100 %
      :alt: 19_PiVR_complete.jpg

#. Plug the 5V power source into the micro USB slot of the Raspberry
   Pi (right side). After a couple of seconds the monitor should
   display a colorful image. Then the operating system installation
   will commence. Select the *Recommended* OS.

   .. figure:: Figures/manual_hardware/standard/20_Install_Raspian.jpg
      :width: 100 %
      :alt: 20_Install_Raspian.jpg

#. On the first startup the OS asks a couple of questions. The
   most important one is the language - make sure you choose the
   correct Keyboard layout. Make sure the Raspberry Pi is
   connected to the internet and download the
   :download:`PiVR installation file <../install_PiVR.sh>`

#. Open the terminal. Then, change directory to the 'Downloads' folder
   (or wherever you downloaded the file) and type::

      bash install_PiVR.sh

   .. figure:: Figures/manual_hardware/standard/21_install_PiVR_software.jpg
      :width: 100 %
      :alt: 21_install_PiVR_software.jpg

#. Now the arena will be built. In the folder PiVR/Lasercutter_Files/
   you can find two vector graphic files that can be used to Lasercut
   a 20cm or 30cm arena, circular holes for M8 screws and small lines
   indicating the distance of 1cm on each side. For one arena you
   will need two acrylic plates.

   .. figure:: Figures/manual_hardware/standard/22_Arena_overview.jpg
      :width: 100 %
      :alt: 22_Arena_overview.jpg

#. Cut the 850nm (infrared) LED strips to the desired length (e.g. 20cm
   on a 20cm Arena) and attach them to the arena. You can choose the
   horizontal distance yourself. I usually use a distance of 3cm.

    .. important::

        It will make soldering much easier if you make sure the
        (+) and the (-) between the LED strips is consistent!

   Solder (+) to (+) from one side of the arena to the other. Then
   repeat for (-).

   .. figure:: Figures/manual_hardware/standard/23_soldered_arena.jpg
      :width: 100 %
      :alt: 23_soldered_arena.jpg

#. Attach the female Barrel Jack to a convenient copper dot on the
   LED strip. Then fix the Female Barrel Jack using super
   glue/hot glue gun. Make sure you are leaving space for the M8
   screw to pass through.

    .. important::

        Usually the red wire of the Barrel Jack indicates (+)!

   .. figure:: Figures/manual_hardware/standard/24_arena_barrel_jack_attached.jpg
      :width: 100 %
      :alt: 24_arena_barrel_jack_attached.jpg

#. If you want to add a Stimulation LED strip (e.g. Red Stimulation
   light), just attach it in between the infrared LED strips, solder
   it as you did the 850nm LED strips and attach the female Barrel
   connector at a convenient location and fix it using the hot glue
   gun.

   .. figure:: Figures/manual_hardware/standard/25_arena_second_strip.jpg
      :width: 100 %
      :alt: 25_arena_second_strip.jpg

#. After inserting the M8 screws into the holes, thread a M8 nut on
   each of the screws about 2cm in. Put the second plate on top of
   the first and fasten it by threading a second M8 nut on top of the
   plate. Make sure the top plate is completely level by using a spirit
   level!

   .. figure:: Figures/manual_hardware/standard/26_second_plate_fastened.jpg
      :width: 100 %
      :alt: 26_second_plate_fastened.jpg

#. To connect PiVR with the arena a cable needs to be constructed.
   You will need two 5.5mm Male Barrel Jack, two 3.5 male solderable
   Barrel Jacks and around 20 Gauge wire.

   .. figure:: Figures/manual_hardware/standard/27_cables_overview.jpg
      :width: 100 %
      :alt: 27_cables_overview.jpg

#. See :ref:`here<detailed_soldering_instructions_cable>` how to
   construct a cable in detail. Briefly, start by cutting a
   reasonable long piece of the wire, e.g. 50cm, but this depends on
   your application. Attach one side of the cable to the 3.5mm barrel
   jack. You may solder it, but be careful to only use minute amounts of
   solder. Then solder on the 5.5mm barrel jack on the other side,
   fixing it using the shrinking tubes.

   .. figure:: Figures/manual_hardware/standard/28_cable_done.jpg
      :width: 100 %
      :alt: 28_cable_done.jpg

#. Start the PiVR software by double clicking on the shortcut on the
   Desktop. Under 'Options' Select
   :ref:`'Optimize Image' <AdjustImageLabel>`.

   .. figure:: Figures/manual_hardware/standard/29_Starting_Software.jpg
      :width: 100 %
      :alt: 29_Starting_Software.jpg

#. Connect the 12V power source (make sure you have an appropriate
   Ampere rating for the amount of LEDs you use!) to the 5.5mm Input
   on the setup. Do not plug it into the wall socket just yet!

   .. Warning::

       Do not plug the 12V power source into the wall socket while you
       are handling the arena wires.

   Then connect the 3.5mm cable with the appropriate receptacle
   closest to the 5.5mm plug. Then plug the other side into the IR
   LEDs on the arena.

   Now you can plug in the 12V power source into the wall socket.

   .. figure:: Figures/manual_hardware/standard/30_connecting_IR_cables.jpg
      :width: 100 %
      :alt: 30_connecting_IR_cables.jpg

#. Turn the camera on ('Cam On'). Then move the 'Backlight Intensity'
   slider to something like 400'000. You should see how the image on
   the top left of the screen lights up.

   .. note::

      Since the camera has a 800nm Longpass filter you shouldn't see
      anything in the camera preview as long as the infrared light of
      the arena is off, **unless** you have a strong source of
      infrared radiation around, e.g. the Sun.

   .. figure:: Figures/manual_hardware/standard/30_1_testing_IR.jpg
      :width: 100 %
      :alt: 30_1_testing_IR.jpg


#. Connect a second 3.5mm cable just below the first. The other side
   goes into the first Stimulation Light in the arena.

   .. figure:: Figures/manual_hardware/standard/31_connecting_redlight_cable.jpg
      :width: 100 %
      :alt: 31_connecting_redlight_cable.jpg

#. When moving the slider labelled 'Channel 1' the stimulation LED
   should light up.

   .. figure:: Figures/manual_hardware/standard/31_1_testing_channel1.jpg
      :width: 100 %
      :alt: 31_1_testing_channel1.jpg

#. If these tests have been successful, congratulations, you've built
   your own PiVR.

.. _detailed_soldering_instructions:

Detailed PCB soldering instructions (Standard Version)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Warning::

   Important! Make sure that the pins are not connected due to
   imprecise soldering!

#. I prefer to solder the components on the PCB board in this
   particular sequence as I find it easiest to keep the components in
   place. Otherwise there is no reason to not solder components in
   any sequence you prefer!

#. To solder the PCB board you will need the following elements:

    #. Four (4) Transistors: 30N06L

       .. figure:: Figures/manual_hardware/standard/2_1_transistors.jpg
          :width: 100 %
          :alt: 2_1_transistors.jpg

    #. One (1) GPIO Header

       .. figure:: Figures/manual_hardware/standard/2_2_GPIO-Header.jpg
          :width: 100 %
          :alt: 2_2_GPIO-Header.jpg

    #. One (1) Barrel connector, 5.5mm Jack with 2.1mm center pole
       diameter

       .. figure:: Figures/manual_hardware/standard/2_3_Jack5_5.jpg
          :width: 100 %
          :alt: 2_3_Jack5_5.jpg

    #. Four (4) Barrel connectors, 3.5mm Jack with 1.35mm center pole

       .. figure:: Figures/manual_hardware/standard/2_4_Jack3_5.jpg
          :width: 100 %
          :alt: 2_4_Jack3_5.jpg

    #. Four (4) :math:`10k{\Omega}` resistors

       .. figure:: Figures/manual_hardware/standard/2_5_resistors.jpg
          :width: 100 %
          :alt: 2_5_resistors.jpg

    #. Break off a 5 pin stretch from the Breakaway Headers. This can
       be done using scissors or pliers.

       .. figure:: Figures/manual_hardware/standard/2_6_Headers.jpg
          :width: 100 %
          :alt: 2_6_Headers.jpg

#. Take one of  the small barrel plug and place it into the
   leftmost possible spot on the PCB board as shown.

   .. figure:: Figures/manual_hardware/standard/S_1.jpg
      :width: 100%
      :alt: S_1.jpg

#. Flip the PCB board while holding the small barrel plug in
   place. By placing it on the table, it should not move and allow
   you to easily solder the three pins of the barrel plug to the PCB
   as shown.

   .. figure:: Figures/manual_hardware/standard/S_2.jpg
      :width: 100%
      :alt: S_2.jpg

#. Continue to solder the other three small barrel plugs, one by one,
   onto the the PCB board.

   .. figure:: Figures/manual_hardware/standard/S_3.jpg
      :width: 100%
      :alt: S_3.jpg

#. Next, place the GPIO header in **exactly** the orientation shown in
   the image below onto the PCB board.

   .. figure:: Figures/manual_hardware/standard/S_4.jpg
      :width: 100%
      :alt: S_4.jpg

#. Flip the PCB board with the GPIO header around. As it now stands
   on the table it should be easy to solder. You do not have to
   solder every single pin to the PCB (minimum is shown on top
   picture) but it is recommended to solder more, ideally all. **Be
   sure the solder between the pins does not touch**

   .. figure:: Figures/manual_hardware/standard/S_5.jpg
      :width: 100%
      :alt: S_5.jpg

#. Place the 5-pin stretch of breakaway headers into the holes on
   the far right on the PCB. Make sure to place them in the correct
   orientation as shown in the picture (short side on the 'back' of
   the PCB).

   .. figure:: Figures/manual_hardware/standard/S_6.jpg
      :width: 100%
      :alt: S_6.jpg

#. Flip the PCB with the 5-pin stretch of breakaway headers around
   and solder the header to the PCB board.

   .. figure:: Figures/manual_hardware/standard/S_7.jpg
      :width: 100%
      :alt: S_7.jpg

#. Now to the resistors. Place a resistor in the indicated position:

   .. figure:: Figures/manual_hardware/standard/S_8.jpg
      :width: 100%
      :alt: S_8.jpg

#. Flip the PCB board around. If the resistor falls out, just fixate
   it by bending the wire as indicated here. Then solder it the the
   PCB board.

   .. figure:: Figures/manual_hardware/standard/S_9.jpg
      :width: 100%
      :alt: S_9.jpg

#. Do the same for the other three resistors.

   .. figure:: Figures/manual_hardware/standard/S_10.jpg
      :width: 100%
      :alt: S_10.jpg

#. Now take the large barrel connector and place it on the PCB at the
   indicated position

   .. figure:: Figures/manual_hardware/standard/S_11.jpg
      :width: 100%
      :alt: S_11.jpg

#. Flip the PCB board around and solder the large barrel connector to
   the PCB board.

   .. figure:: Figures/manual_hardware/standard/S_12.jpg
      :width: 100%
      :alt: S_12.jpg

#. Next, take one of the transistors and place it **exactly** as
   shown onto the PCB board.

   .. figure:: Figures/manual_hardware/standard/S_13.jpg
      :width: 100%
      :alt: S_13.jpg

#. Flip the PCB board around and solder the transistor to the board.
   Make sure the solder of the different pins does not touch the
   contact of one of the other pins! **Warning: Transistors are more
   heat sensitive compared to the other components you have used so
   far. Make sure to not let them heat up too much!**

   .. figure:: Figures/manual_hardware/standard/S_14.jpg
      :width: 100%
      :alt: S_14.jpg

#. Do the same for the other three transistors.

   .. figure:: Figures/manual_hardware/standard/S_15.jpg
      :width: 100%
      :alt: S_15.jpg

#. You must get rid of the elongated wiring of the transistors and
   especially the resistors as not doing so will 1) increase risk of
   shorting components and 2) it will physically be very hard to put
   the PCB board into the casing. While it is probably best to use
   the shown wire clipper, it is also possible to do that using
   normal scissor.

   .. figure:: Figures/manual_hardware/standard/S_16.jpg
      :width: 100%
      :alt: S_16.jpg

#. Congratulations - You have finished soldering your very own PCB
   board to build PiVR!

.. _PiVR HP Version:

Building the High Powered Version
---------------------------------

.. Warning::

   **PLEASE** avoid looking directly into the LED light. It might
   damage your eyes.

   If you are in a situation where you might accidentally look into
   the LED light at 100% power, please consider purchasing eye saftey
   googles for the wavelength you are using, for example from `here
   <https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=762>`__.

.. Warning::

    You will be handling sensitive electronic equipment. As you might
    be electrically charged, always touch something metallic that is
    connected to the ground before unpacking anything from an
    electrostatical shielding bag. A radiator will do fine.

.. note::

    You can of course change the order of how to build the different
    parts. I build them in this order as there are necessary delays
    such as installing the operation system on the Raspberry Pi.

.. note::

    We used an Ultimaker 3 with a Ultimaker Print Core AA (0.8mm) to
    print the PLA and a Ultimaker Print Core BB (0.8) to print the
    PVA support material. STL files were converted using Ultimaker
    CURA software.

#. 3D print each part in the folder PiVR/Hardware/3D printer files.
   We used an Ultimaker 3 and printed with Standard PLA. We used PVA as
   support material. For best results print Casing, CameraHolder,
   TowerBottom and TowerTop with support material.

    .. figure:: Figures/manual_hardware/standard/1_3D_printed_parts.jpg
       :width: 100 %
       :alt: 1_3D_printed_parts.jpg

#. Obtain the Printed Circuit Board (PCB). Find the blueprint in
   PiVR\Hardware\PCB. I have used the company `Aisler.net
   <https://aisler.net/>`_ to print the circuit boards.

    .. figure:: Figures/manual_hardware/high_power/2_1_PCB_overview.jpg
       :width: 100 %
       :alt: 2_1_PCB_overview.jpg

#. Get the following items to solder the PCB board:

    .. important::

        Before touching electrical components always touch something
        metalic that is connected to the ground. For example a radiator.

    #. Three (3) LED driver "MiniPuck", 700mA (or whatever strength
       you need)
    #. One (1) Transistors: 30N06L
    #. One (1) GPIO Header
    #. One (1) Barrel connector, 5.5mm Jack with 2.1mm center pole
       diameter
    #. Four (4) Barrel connectors, 3.5mm Jack with 1.35mm center pole
    #. One (1) :math:`10k{\Omega}` resistors
    #. Break off a 5 pin stretch from the Breakaway Headers. This can
       be done using scissors or pliers.

   .. figure:: Figures/manual_hardware/high_power/2_2_PCB_components_HP.jpg
      :width: 100 %
      :alt: 2_7_PCB_components_HP.jpg

#. Solder the components on the PCB. See
   :ref:`this section<detailed_soldering_instructions_HP>` for
   detailed soldering instructions.

    .. important::

        The correct **orientation** of the **GPIO header** and the
        **Transistor** is crucial for PiVR to work correctly.

   .. figure:: Figures/manual_hardware/high_power/3_solderedPCB_HP.jpg
      :width: 100 %
      :alt: 3_solderedPCB_HP.jpg

#. Cut off the excess wire of the resistors and the Transistors, e.g.
   with scissors.

   .. figure:: Figures/manual_hardware/high_power/4_PCB_removed_wire_HP.jpg
      :width: 100 %
      :alt: 4_PCB_removed_wire_HP.jpg

#. Unpack the Touchscreen, remove the stand-off screws. Attach the
   monitor cable that came with the Touchscreen and the 4"
   5 pin cable.

    .. important::

        The monitor cable must be inserted in the correct orientation.
        When you look into the receptacle you'll see that only one
        side has connectors. Make sure that you insert the cable's
        connectors on the same side.

    .. important::

        Note the orientation of the 4" 5pin cable! Left is (+) while
        Right is (-).

   .. figure:: Figures/manual_hardware/standard/6_Touchscreen.jpg
      :width: 100 %
      :alt: 6_Touchscreen.jpg

#. Place the Casing on top of the Touchscreen (it will only fit in
   the shown orientation). Organize the 4" 5pin cable and the
   monitor cable as shown in the picture. Use the M2.5x10mm screws to
   fix the casing to the touchscreen.


   .. figure:: Figures/manual_hardware/standard/7_2_Casing_on_Touchscreen.jpg
      :width: 100 %
      :alt: 7_2_Casing_on_Touchscreen.jpg

#. Prepare the SD card: Format the SD card using SD Formatter and load with
   NOOBs installation files as instructed
   `here <https://www.raspberrypi.org/help/noobs-setup/2/>`_:

#. Connect monitor cable with the Raspberry Pi (with inserted SD
   card). Again, make sure you insert the cable in the correct
   orientation. Use M2.5x10 screws to attach the Raspberry Pi to the
   Casing.

   .. figure:: Figures/manual_hardware/standard/8_0_Raspberry_monitor_cable.jpg
      :width: 100 %
      :alt: 8_0_Raspberry_monitor_cable.jpg

#. Attach the PCB board on the right side of the casing using M2.5x10mm
   screws. Plug the 4" 5pin cable into the PCB  in the **correct**
   orientation as shown!

   .. figure:: Figures/manual_hardware/high_power/9_PCB_RPi_on_casing_HP.jpg
      :width: 100 %
      :alt: 9_PCB_RPi_on_casing_HP.jpg

#. Use the GPIO Ribbon cable to connect the PCB board with the
   Raspberry Pi. Thread the long camera cable through the slit as
   shown in the image below. Connect it to the Raspberry Pi Camera port.

   .. figure:: Figures/manual_hardware/high_power/10_Casing_stuffed.JPG
      :width: 100 %
      :alt: 10_Casing_stuffed.JPG

#. Attach the Backside to the CasingPedestal with M2.5x10 screws.

   .. figure:: Figures/manual_hardware/standard/13_Backside_attached_to_Pedestal.jpg
      :width: 100 %
      :alt: 13_Backside_attached_to_Pedestal.jpg

#. Now just slide the Backside (with attached CasingPedestal) into the
   Casing.

   .. figure:: Figures/manual_hardware/high_power/11_Backside_into_casing.JPG
      :width: 100 %
      :alt: 11_Backside_into_casing.JPG


#. Drop a 2.5mm nut in each hole in the TowerPedestal. Use the M2.5x10
   screws to attach the TowerBottom to the Tower Pedestal

   .. figure:: Figures/manual_hardware/standard/14_TowerPedestal_attached.jpg
      :width: 100 %
      :alt: 14_TowerPedestal_attached.jpg

#. Using a hammer, drive the dowel pins into the TowerBottom. Then
   attach the TowerTop to it. In principle you can stack more TowerTops
   on top.

   .. figure:: Figures/manual_hardware/standard/15_TowerBottom_Dowel.jpg
      :width: 100 %
      :alt: 15_TowerBottom_Dowel.jpg

#. Attach the 800nm Longpass Filter to the Camera using Parafilm. It
   is best to wear gloves for this step.

   .. Important::

      If you are **not** using the camera or lens in the image below,
      please :ref:`providee<create_own_undistort_files_label>`
      your own undistort files before running experiments.

   .. figure:: Figures/manual_hardware/standard/16_Camera_with_LP_filter.jpg
      :width: 100 %
      :alt: 16_Camera_with_LP_filter.jpg

#. Thread the camera cable from the Casing through the slit in the
   TowerBottom and through the slit of the Camera Holder.

    .. important::

        Note the orientation to avoid having to curl the camera cable
        in the camera holder

   .. figure:: Figures/manual_hardware/standard/17_Threading_of_Cam_cable.jpg
      :width: 100 %
      :alt: 17_Threading_of_Cam_cable.jpg

#. Attach the Camera Cable to the Camera in the **correct**
   orientation. Then screw the camera to the Camera Holder using the
   M2.5x10 screws. It is **not** necessary to fixate the screws with
   nuts!

   .. figure:: Figures/manual_hardware/standard/18_Camera_attached_to_CamHolder.jpg
      :width: 100 %
      :alt: 18_Camera_attached_to_CamHolder.jpg

#. Drop a 2.5mm nut in the hole in the Camera holder and use it to
   fasten the M2.5x10 screw. Then attach the CameraHolder to the
   Tower.

   .. figure:: Figures/manual_hardware/standard/19_PiVR_complete.jpg
      :width: 100 %
      :alt: 19_PiVR_complete.jpg

#. Plug the 5V power source into the micro USB slot of the Raspberry
   Pi (right side). After a couple of seconds the monitor should
   display a colorful image. Then the operating system installation
   will commence. Select the *Recommened* OS.

   .. figure:: Figures/manual_hardware/standard/20_Install_Raspian.jpg
      :width: 100 %
      :alt: 20_Install_Raspian.jpg

#. On the first startup the OS asks a couple of questions. The
   most important one is the language - make sure you choose the
   correct Keyboard layout. Make sure the Raspberry Pi is
   connected to the internet and download the
   :download:`PiVR installation file <../install_PiVR.sh>`

#. Open the terminal. Then change directory to the 'Downloads' folder
   (or wherever you downloaded the file) and type::

      bash install_PiVR.sh

   .. figure:: Figures/manual_hardware/standard/21_install_PiVR_software.jpg
      :width: 100 %
      :alt: 21_install_PiVR_software.jpg

#. While the PiVR software is being installed on the Raspberry Pi (it
   will take > 10 minutes) the light pad can be built. Take the
   diffuser. If you have access to a lasercutter, you can use
   the *svg* file in 'Hardware\Lasercutter_Files' to prepare it.
   Alternatively, you can use a drill with a M8 or M8.5 drill bit to
   drill 4 holes as shown below.

   .. figure:: Figures/manual_hardware/high_power/12_diffuser.jpg
      :width: 100 %
      :alt: 12_diffuser.jpg

#. Put in the M8 screws with a nut below and above the diffuser as
   shown below.

   .. figure:: Figures/manual_hardware/high_power/12_1_diffuser_screw_example.jpg
      :width: 100 %
      :alt: 12_1_diffuser_screw_example.jpg

#. Next, take the Aluminum heat sink and use a Sharpie and a scale to
   draw a grid with 1cm squares on the heat sink.

   .. figure:: Figures/manual_hardware/high_power/12_2_LED_heat_sink.jpg
      :width: 100 %
      :alt: 12_2_LED_heat_sink.jpg

#. Take the thermally conductive double sided tape. I usually first
   attach it to the LED before sticking them on the aluminum heat
   sink.

   .. Note::

      While you can chose LEDs of any wavelength, there is only space
      (and power) for up to 9 LEDs. In the interest of homogeneous
      illumination only a single wavelength can be used per light pad.

   Do make sure that the LEDs are space 5 cm apart using the
   previously drawn grid.

   .. figure:: Figures/manual_hardware/high_power/12_3_LEDs_on_heat_sink.jpg
      :width: 100 %
      :alt: 12_3_LEDs_on_heat_sink.jpg

#. Now take the 850 nm LED strip for background illumination and
   place it in between the high powered LED lights.

   .. figure:: Figures/manual_hardware/high_power/12_4_background_LEDs.jpg
      :width: 100 %
      :alt: 12_4_background_LEDs.jpg

#. Now, for the first row of LEDs, do the following: connect the (+) of
   the LED on the right with the (-) with the LED on its left.

   .. note::

      As heat transfers very quickly from the LED to the aluminum
      board it can be hard to liquify the solder on the LEDs. If you
      run into this problem, put some solder on the wire as well
      and keep it liquid until touching the contact on the LED.
      (Thank you Ruud Schilder for this suggestion)

   .. figure:: Figures/manual_hardware/high_power/12_6_LED_connected.jpg
      :width: 100 %
      :alt: 12_6_LED_connected.jpg

#. For the LED on the right, keep a longer, wire from (-). For the LED
   on the left, keep a longer wire from the (+).

   .. figure:: Figures/manual_hardware/high_power/12_7_one_row_done.jpg
      :width: 100 %
      :alt: 12_7_one_row_done.jpg

#. Take one of the smaller heat shrinking tubes and put it over the
   wire coming from the (+) of the LED. Take a female barrel plug and
   solder the red wire with the wire coming from the (+) of the LED
   row. When the solder has cooled down, push the heat shrink wire
   over the just soldered wire and heat it for stability of the cable.

   .. figure:: Figures/manual_hardware/high_power/12_8_connect_positive.jpg
      :width: 100 %
      :alt: 12_8_connect_positive.jpg

#. Repeat this for the wire coming from the (-) of the LED and
   connect it to the black wire of the barrel plug.

   .. figure:: Figures/manual_hardware/high_power/12_9_row_done.jpg
      :width: 100 %
      :alt: 12_9_row_done.jpg

#. Now is a good time to test whether everything works as expected so
   far. First, head over to your PiVR setup. The installation should
   have finished by now. Before plugging in the power supply you need
   to change a setting in the PiVR software:

   Go to "Options->High Power LEDs".

   In the popup window select "High Powered LED configuration" and hit
   "Exit and save changes".

   To restart the PiVR software, go to "File->Save and Exit".

   .. figure:: Figures/manual_hardware/high_power/13_PiVR_HP_settings.jpg
      :width: 100 %
      :alt: 13_PiVR_HP_settings.jpg

#. You also have to manually adapt the Pulse Width Modulation (PWM)
   frequency to the manufacturers recommendation. As you are using a
   MiniPuck, the
   `manual <https://www.ledsupply.com/content/pdf/MiniPuck_F004.pdf>`__
   states that the PWM frequency must be between 100-2kHz. See the
   :ref:`software manual<DefineGPIOsLabel>` on how to modify the PWM
   frequency.

#. Next, you will have to prepare a cable as described
   :ref:`here<detailed_soldering_instructions_cable>`.

#. Connect an appropriate power supply (make sure it is 12V and has
   enough ampere for your particular setup) to the wall and to the
   PiVR setup.

   In this particular example we are using 3 MiniPuck that have an
   output current of 700mA. The 850nm IR LEDs are rated at 24W/5M=2A of
   which we are using 3 x 20cm which is 2A x 0.06 = 120mA. The power
   supply needs to provide at least 3 x 700mA + 120mA = 2220mA. It is
   recommended to never be to close to the limit, hence I chose a 5A
   wall adapter.

   .. figure:: Figures/manual_hardware/high_power/14_power_supply.jpg
      :width: 100 %
      :alt: 14_power_supply.jpg

#. Connect the cable to the *second small connector from the top*.
   The large barrel connector needs to be connected to the female
   barrel plug you previously soldered to the LEDs.

   .. Warning::

        Important: Only connect PiVR with the LED after changing the
        settings in "Options->High Power LEDs" to "High Powered LED
        configuration". Otherwise the (very bright) LEDs will be
        completely turned on immediately.

   .. figure:: Figures/manual_hardware/high_power/15_connect_Channel_1.jpg
      :width: 100 %
      :alt: 15_connect_Channel_1.jpg

#. To test the LED, go to "Options->Optimize Image" and use the
   slider below "Channel 1". You should see the LED light up. If it
   does not, see :ref:`here<troubleshooting_LEDs>` for some
   possible problems that might cause the failure.

   .. figure:: Figures/manual_hardware/high_power/16_testing.jpg
      :width: 100 %
      :alt: 16_testing.jpg

#. After confirming that the first row of LEDs works as expected,
   continue by soldering the other rows up exactly in the same
   fashion. You will quickly notice that there are too many cables
   around. I find that using super  glue to organize the wires really
   does help. I recommend using gloves while handling super glue.

   Do make sure to keep the edges of the aluminum heat plate free as
   indicated with the red crosses below.

   It's a good idea to test each row individually for functionality
   to quickly detect soldering errors.

   .. figure:: Figures/manual_hardware/high_power/17_superglueing.jpg
      :width: 100 %
      :alt: 17_superglueing.jpg

#. Finally, connect the 850nm infrared LED strips. From each row
   connect one of the (+) with the (+) on the next row and the repeat
   for the (-). Connect the red wire of the barrel plug with the (+)
   of the LED strip.

   .. figure:: Figures/manual_hardware/high_power/18_850nm_LED_strip.jpg
      :width: 100 %
      :alt: 18_850nm_LED_strip.jpg

#. To test whether these LEDs work, take one of your cables and
   connect it to the top plug as shown below. Since these are infrared
   LEDs you will not be able to see with the bare eye, you need to use
   the PiVR camera. Then, go to "Options->Optimize Image" and drag the
   slider below 'Background' to the right. The LEDs should light up
   in the PiVR camera as shown in the image below.

   .. figure:: Figures/manual_hardware/high_power/19_850nm_LED_test.jpg
      :width: 100 %
      :alt: 19_850nm_LED_test.jpg

#. Now you are ready to test the complete setup: Make sure that
   the cable powering the infrared LEDs is connected to the top
   outlet at PiVR (top arrow).

   The high-power LEDs can be connect in any order to the 3 remaining
   outlets (the three bottom arrows).

   The diffuser should now be placed on top of the aluminum heat sink.

   .. figure:: Figures/manual_hardware/high_power/20_example_working.jpg
      :width: 100 %
      :alt: 20_example_working.jpg

#. If these tests have been successful, congratulations, you've built
   your own high-power PiVR.

.. _detailed_soldering_instructions_HP:

Detailed PCB soldering instructions (High Power Version)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Warning::

   Important! Make sure that the pins are not connected due to
   imprecise soldering!

#. I prefer to solder the components on the PCB board in this
   particular sequence as I find it easiest to keep the components in
   place. Otherwise there is no reason to not solder components in
   any sequence you prefer!

#. To solder the PCB board you will need the following elements:

    #. One (1) Transistor: 30N06L

       .. figure:: Figures/manual_hardware/high_power/2_3_transistor_HP.jpg
          :width: 100 %
          :alt: 2_3_transistor_HP.jpg

    #. Three (3) LED Drivers "MiniPuck"

       .. figure:: Figures/manual_hardware/high_power/2_4_LEDDrivers_HP.jpg
          :width: 100 %
          :alt: 2_4_LEDDrivers_HP.jpg

    #. One (1) GPIO Header

       .. figure:: Figures/manual_hardware/standard/2_2_GPIO-Header.jpg
          :width: 100 %
          :alt: 2_2_GPIO-Header.jpg

    #. One (1) Barrel connector, 5.5mm Jack with 2.1mm center pole
       diameter

       .. figure:: Figures/manual_hardware/standard/2_3_Jack5_5.jpg
          :width: 100 %
          :alt: 2_3_Jack5_5.jpg

    #. Four (4) Barrel connectors, 3.5mm Jack with 1.35mm center pole

       .. figure:: Figures/manual_hardware/standard/2_4_Jack3_5.jpg
          :width: 100 %
          :alt: 2_4_Jack3_5.jpg

    #. One (1) :math:`10k{\Omega}` resistors

       .. figure:: Figures/manual_hardware/high_power/2_5_Resistance.jpg
          :width: 100 %
          :alt: 2_5_Resistance.jpg

    #. Break off a 5 pin stretch from the Breakaway Headers. This can
       be done using scissors or pliers.

       .. figure:: Figures/manual_hardware/standard/2_6_Headers.jpg
          :width: 100 %
          :alt: 2_6_Headers.jpg

#. Take one of  the small barrel plug and place it into the
   leftmost possible spot on the PCB board as shown.

   .. figure:: Figures/manual_hardware/high_power/S_1.jpg
      :width: 100%
      :alt: S_1.jpg

#. Flip the PCB board while holding the small barrel plug in
   place. By placing it on the table, it should not move and allow
   you to easily solder the three pins of the barrel plug to the PCB
   as shown.

   .. figure:: Figures/manual_hardware/high_power/S_2.jpg
      :width: 100%
      :alt: S_2.jpg

#. Continue to solder the other three small barrel plugs, one by one,
   onto the the PCB board.

   .. figure:: Figures/manual_hardware/high_power/S_3.jpg
      :width: 100%
      :alt: S_3.jpg

#. Next, place the GPIO header in **exactly** the orientation shown in
   the image below onto the PCB board.

   .. figure:: Figures/manual_hardware/high_power/S_4.jpg
      :width: 100%
      :alt: S_4.jpg

#. Flip the PCB board with the GPIO header around. As it now stands
   on the table it should be easy to solder. **Be
   sure the solder between the pins does not touch**

   .. figure:: Figures/manual_hardware/high_power/S_5.jpg
      :width: 100%
      :alt: S_5.jpg

#. Place the 5-pin stretch of breakaway headers into the holes on
   the far right on the PCB. Make sure to place them in the correct
   orientation as shown in the picture (short side on the 'back' of
   the PCB).

   .. figure:: Figures/manual_hardware/standard/S_6.jpg
      :width: 100%
      :alt: S_6.jpg

#. Flip the PCB with the 5-pin stretch of breakaway headers around
   and solder the header to the PCB board.

   .. figure:: Figures/manual_hardware/high_power/S_7.jpg
      :width: 100%
      :alt: S_7.jpg

#. Now take the large barrel connector and place it on the PCB at the
   indicated position

   .. figure:: Figures/manual_hardware/high_power/S_8.jpg
      :width: 100%
      :alt: S_8.jpg

#. Flip the PCB board around and solder the large barrel connector to
   the PCB board.

   .. figure:: Figures/manual_hardware/high_power/S_9.jpg
      :width: 100%
      :alt: S_9.jpg

#. Take the MiniPuck LED drivers and put them into the labelled position
   on the PCB.

   .. figure:: Figures/manual_hardware/high_power/S_10.jpg
      :width: 100%
      :alt: S_10.jpg

#. Turn the PCB around and solder each MiniPuck LED drivers to the PCB.

   .. figure:: Figures/manual_hardware/high_power/S_11.jpg
      :width: 100%
      :alt: S_11.jpg

#. After soldering the MiniPuck LED drivers, flip the PCB again. You
   will notice that they have relatively long pins (top part of image).
   You need to remove them using a wire cutter.

   .. figure:: Figures/manual_hardware/high_power/S_12.jpg
      :width: 100%
      :alt: S_12.jpg

#. Next, the :math:`10k{\Omega}` resistor will be soldered to the PCB.

   .. figure:: Figures/manual_hardware/high_power/S_13.jpg
      :width: 100%
      :alt: S_13.jpg

#. Flip the PCB board around. If the resistor falls out, just fixate
   it by bending the wire as indicated here. Then solder it the the
   PCB board.

   .. figure:: Figures/manual_hardware/high_power/S_14.jpg
      :width: 100%
      :alt: S_14.jpg

#. Similar to the MiniPuck LED drivers above, you need to remove the
   leftover wire from the resistor after soldering it to the PCB.

   .. figure:: Figures/manual_hardware/high_power/S_15.jpg
      :width: 100%
      :alt: S_15.jpg

#. Next, take transistor and place it **exactly** as
   shown onto the PCB board.

   .. figure:: Figures/manual_hardware/high_power/S_16.jpg
      :width: 100%
      :alt: S_16.jpg

#. Flip the PCB board around and solder the transistor to the board.
   Make sure the solder of the different pins does not touch the
   contact of one of the other pins! **Warning: Transistors are more
   heat sensitive compared to the other components you have used so
   far. Make sure to not let them heat up too much!**

   .. figure:: Figures/manual_hardware/high_power/S_17.jpg
      :width: 100%
      :alt: S_17.jpg

#. You must get rid of the elongated wiring of the transistors
   as not doing so will make it physically be very hard to put
   the PCB board into the casing. While it is probably best to use
   the shown wire clipper, it is also possible to do that using
   normal scissor.

   .. figure:: Figures/manual_hardware/high_power/S_18.jpg
      :width: 100%
      :alt: S_18.jpg

#. Congratulations - You have finished soldering your very own PCB
   board to build the high-powered version of PiVR!

.. _LED_measurements:

Testing homogeneity of the stimulation light
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The principle of creating virtual realities depends on the homogeneity
of the stimulation light. Below I will lay out how to measure the
light intensity and strategies on how to handle inhomogeneity.

You can use any light/power meter. We have used the PM100D with the
S130VC sensor (both from Thorlabs).

.. figure:: Figures/manual_hardware/homogeneity/1_tools.jpg
   :width: 100%
   :alt: 1_tools.jpg

Draw a grid on the diffuser. I like to use numbers for rows and the
alphabet for columns for clear separation while writing down the
light intensity.

Wear eye protection and **set the LED to 100%**.

.. Warning::

   If you are using a unusual wavelength, e.g. deep blue at ~400nm or
   below you must also protect your skin! Please consult the LED
   manufacturer safety recommendation sheet.

.. Note::

   Why not using a low light intensity, e.g. 5%? PiVR is using Pulse
   Width Modulation (PWM) to control light intensity. The light meter
   is integrating the number of photons it reads over some time t
   (which can be hard to find in the light meter manual). PWM turns
   the LED on and off rapidly (as defined
   :ref:`here<DefineGPIOsLabel>`). The easiest workaround is to just
   measure at 100% light intensity. Since PWM is guaranteed to scale
   linearly this will lead to accurate predictions of light
   intensities at e.g. 10% power.

Measure the light intensity for each square in the grid and note the
light intensity. This will give you a good idea of the light
intensity provided by the light pad and the homogeneity.

.. figure:: Figures/manual_hardware/homogeneity/2_measurement_example.jpg
   :width: 100%
   :alt: 2_measurement_example.jpg

As you can quickly see, the light intensity is not perfectly
homogeneous.

.. figure:: Figures/manual_hardware/homogeneity/3_Single_diffuser.jpg
   :width: 100%
   :alt: 3_Single_diffuser.jpg

There are several ways to optimize for homogeneity: The easiest
option is to add another diffuser plate. This will usually improve
homogeneity but knock off quite a bit of the light intensity. If you
compare the light intensity you can see that peak light intensity
decreased from :math:`40{\mu}W/mm^2` to :math:`10{\mu}W/mm^2`

.. figure:: Figures/manual_hardware/homogeneity/4_double_diffuser.jpg
   :width: 100%
   :alt: 4_double_diffuser.jpg

Finally, the reason why the light intensity falls of at the edges is
due to the fact that quite a bit of light can escape by the sides
instead of going through the diffuser. To counter this, the a
reflective border can be used - just cardboard and aluminum foil.

.. figure:: Figures/manual_hardware/homogeneity/5_double_diffuser_border.jpg
   :width: 100%
   :alt: 5_double_diffuser_border.jpg

Conclusion: While the homogeneity in the center (indicated with a
circle of diameter 90mm) is reasonable, it falls off at the edges. By
trying different heights of the second diffuser plate and by
introducing a border you can try to mitigate this.

There is an obvious software solution to this problem: For each
plate, do these measurements save on the PiVR setup. When creating
the virtual arena, tie stimulus intensity could take these
measurments into account and normalize to the maximum intensity. This
might be implemented in the future. Please check the gitlab
repository for any open issues if you need this functionality.

Additional Information
----------------------

.. _detailed_soldering_instructions_cable:

Detailed soldering instructions (cable)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. To construct 4 cables that connects PiVR with the light pad you
   will need the following components:

   .. figure:: Figures/manual_hardware/cable/1_overview.jpg
      :width: 100%
      :alt: 1_overview.jpg

   #. For each cable you will need one large barrel plug,

      .. figure:: Figures/manual_hardware/cable/1_3_large_barrel_connectors.jpg
        :width: 100%
        :alt: 1_3_large_barrel_connectors.jpg

   #. one small barrel plug, consisting of the plug (on the right side) and the holder (left side),

      .. figure:: Figures/manual_hardware/cable/1_2_small_barrel_connectors.jpg
        :width: 100%
        :alt: 1_2_small_barrel_connectors.jpg

   #. one large and two smaller heat shrinking tubes,

      .. figure:: Figures/manual_hardware/cable/1_4_heat_shrink_tubing.jpg
        :width: 100%
        :alt: 1_4_heat_shrink_tubing.jpg

   #. and some wire:

      .. figure:: Figures/manual_hardware/cable/1_5_wire.jpg
        :width: 100%
        :alt: 1_5_wire.jpg

#. Cut ~50 cm (or however long you want the wire to be). Take the
   3.5mm barrel connector and place the red wire into the center hole
   of the connector. Place the black wire into the hole on the side
   of the connector.

   .. figure:: Figures/manual_hardware/cable/2_setting_it_up.jpg
      :width: 100%
      :alt: 2_setting_it_up.jpg

#. Solder the red wire to the 3.5mm barrel connector. If there is a
   lot of leftover wire as in the image below, it is a good idea to
   cut it away.

   .. figure:: Figures/manual_hardware/cable/3_solder_pos_cable.jpg
      :width: 100%
      :alt: 3_solder_pos_cable.jpg

#. Now solder the black wire to the 3.5mm barrel connector. Again,
   make sure to remove leftover wire.

   Finally, thread the cylindrical black plastic screw on top of the
   just soldered 3.5mm barrel connector.

   .. figure:: Figures/manual_hardware/cable/4_solder_neg.jpg
      :width: 100%
      :alt: 4_solder_neg.jpg

#. After screwing the back onto the barrel connector, it should look
   like this:

   .. figure:: Figures/manual_hardware/cable/5_small_side_finished.jpg
      :width: 100%
      :alt: 5_small_side_finished.jpg

#. On the opposite side of the cable, start by placing the heat
   shrink tubes as indicated below: One small set of tubes on each of
   the thinner wires (left arrows) and one larger one on the male
   large barrel plug (right arrow).

   .. figure:: Figures/manual_hardware/cable/6_heat_shrink_tubes.jpg
      :width: 100%
      :alt: 6_heat_shrink_tubes.jpg

#. Solder the red wire of the large barrel connector to the (+) wire.
   In our example that would be red wire.

   .. figure:: Figures/manual_hardware/cable/7_large_side_soldered.jpg
      :width: 100%
      :alt: 7_large_side_soldered.jpg

#. Using a heat source, use the heat shrink tubing to stablize the
   cable: Start with the thinner wires and finally push the
   larger heat shrink tube over the other two.

   .. figure:: Figures/manual_hardware/cable/8_shrinking_the_tubes.jpg
      :width: 100%
      :alt: 8_shrinking_the_tubes.jpg

.. _troubleshooting_LEDs:

Troubleshooting LEDs
^^^^^^^^^^^^^^^^^^^^

If you can not turn the LED on using PiVR, the first test should be
whether there is an error in soldering the LEDs. To test this, plug
the 12V wall adaptor directly into the LEDs on the illumination plate.

If the LEDs do *not* turn on, it is very likely you have a loose
connection on the illumination pad. Carefully check each solder patch
. It is also possible you mixed up (+) and (-) somehwere. Finally, it
is possible that you wall adaptor is broken, however, this never
happened to me.

If the LEDs do turn on when you connect the 12V wall adaptor directly
to the LED, there are several options what could be going wrong: (1)
The cable might have a faulty connection. (2) The PCB board might
have a faulty soldering connection.

To address (1), you can use the 'resistance' option on your voltmeter
to check if the two ends of the cable are electrically connected. If
you read '0' (or a very small number), the cable is good.

To address (2), you will have to take out the PCB of the PiVR setup
again. Carefully inspect the soldering. Make sure everything that can
be soldered is soldered and that no two solder patches touch each
other. If that checks out good, double check if you soldered the GPIO
header onto the board the correct way.

If you are still not getting any light out of the LEDs, feel free to
submit at ticket on
`the PiVR repository <https://gitlab.com/LouisLab/pivr/-/issues>`__


PiVR has been developed by David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

.. _BOM:

Bill of Materials (BOM)
========================

All components you need to build your

#. :ref:`Standard PiVR <BOM_std>`

#. :ref:`High Powered PiVR <BOM_HP>`

.. _BOM_std:

Standard PiVR setup
--------------------

For detailed information, scroll to the right.

This table is identical to the one found on
`Gitlab <https://gitlab.com/LouisLab/pivr/blob/master/PiVR/Hardware/
Inventory>`__ which might be more convenient when ordering parts.

.. csv-table:: Table Title
   :file: ../../Hardware/Inventory/PiVR_inventory_standard.csv
   :header-rows: 1

.. _BOM_HP:

High Powered PiVR setup
-----------------------

This table is identical to the one found on
`Gitlab <https://gitlab.com/LouisLab/pivr/blob/master/PiVR/Hardware/
Inventory>`__ which might be more convenient when ordering parts.

.. csv-table:: Table Title
   :file: ../../Hardware/Inventory/PiVR_inventory_HighPower.csv
   :header-rows: 1PiVR has been developed by David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

.. _tools:

Tools
*****

.. important::

    The software has different functionality if run on a Raspberry Pi
    as compared to any other PC. This manual is for the PC (Windows,
    MacOS and Linux if *not* run on a Raspberry Pi) version of the
    software

The Menubar
===========

To select a different window use the Menu Bar at the top of the window

.. figure:: Figures/tools/1_menubar.png
   :alt: 1_menubar.png

The Analysis Menu
=================

The Analysis Menu lets you choose a folder (either with one
experiment or a folder containing several experiments) and run
different analyses.

.. figure:: Figures/tools/2_AnalysisFrame.png
   :alt: 2_menubar.png


To analyze an experiment (or several) first press the "Press to select
data to analyze" button on the top left. Select a folder. You may:

#. Select a single experiment. In that case you **must** uncheck the
   button on the right next to "More than one folder".

#. Select a folder containing several experiments. If this folder
   only contains folders you want to apply the same analysis, you do
   not have to do anything. If there are other files or folders,
   please indicate a commonality between all the folders in the entry
   field below. This can be the date (e.g. 2019) a genotype (e.g.
   Or42a) or..

   The "Number of files" below the entry field indicates how many
   files and folders in the folder you selected is taken into account
   if the current input is used.

Next, select what analysis you want to do. Currently there are several
options:

Distance to source
-------------------

If you are interested in the distance the animal has to a point in
its environment try this option.

We used it in our publication to produce `Figure 2D
<https://journals.plos.org/plosbiology/article/figure?id=10.1371/journal.pbio.3000712.g002>`__:
A fruit fly larva behaving in an environment with a single odor
source. We were interested in the distance to this odor source during
the experiment.

To use this analysis program, you must have used the
:ref:`Tracking<TrackingLabel>` tab on PiVR **or** you must have
recorded a :ref:`Video<VideoLabel>` or taken
:ref:`Full Frame Images<FullFrameLabel>` with subsequent
:ref:`Post-Hoc Single Animal Analysis<tools_single_animal_tracking>`.

.. important::

   If you want to use the multiple experiment options, please ensure
   that you have **identical** Framerates and **identical** Recording
   length

After selecting a folder and pressing "Press to start analysis" a
popup will emerge that will ask you to "Select a source". Do the
following:

#. Press on the button that says "Select Source"

#. Click on the image where the "Source" is. In our case this would
   be an odor source roughly in the center of the dish. In your case
   it could be any *single* point

#. Once you are satisfied with the arrow placement, press "Analyze"

.. figure:: Figures/tools/3_DistanceToSource.png
   :alt: 3_DistanceToSource.png

Once the analysis has finished you will find:

#. "distance_to_source.csv" in each experimental folder. This
   csv file contains the calculated distance to source in mm for each
   frame.

#. "Distance_to_source.png" in each experimental folder. This plot is
   intended to give a quick overview about the distance to source of
   this experiment.

If you have analyzed multiple experiments you will also find:

#. "all_distance_to_source.csv" in the parental folder. This file
   contains the calculated distance to source in mm for each frame
   (rows) for all the analyzed experiments (columns)

#. "Median_Distance_to_source.png" in the parental folder. This plot
   is intended to give a quick overview about the distance to source
   of all the analyzed experiments. Individual experiments are
   plotted in light grey and the median in red.

See here for actual code:
:class:`analysis_scripts.AnalysisDistanceToSource`

Distance to source, VR
-----------------------

If you ran a virtual reality experiment with a single maximum point
and you wish to calculate the distance to that point, try this
option.

We used this analysis in our publication to produce `Figure 2E
<https://journals.plos.org/plosbiology/article/figure?id=10.1371/journal.pbio.3000712.g002>`__:
A fruit fly larva behaving in an virtual
odor reality. We were interested in the distance to this virtual odor
source during the experiment.

To use this analysis program, you must have used the
:ref:`Virtual Reality <VRLabel>` tab on PiVR.

.. important::

   If you want to use the multiple experiment options, please ensure
   that you have **identical** Framerates and **identical** Recording
   length

After selecting a folder and pressing "Press to start analysis" the
script will automatically read the virtual reality arena presented
during the experiment and calculate the distance of the animal to the
maximum point of stimulus intensity in the virtual reality.

Once the analysis has finished you will find:

#. "distance_to_VR_max.csv" in each experimental folder. This
   csv file contains the calculated distance to the maximum point of
   stimulus intensity in the virtual reality in mm for each frame.

#. "Distance_to_VR_max.png" in each experimental folder. This plot is
   intended to give a quick overview about the distance to source of
   this experiment.

If you have analyzed multiple experiments you will also find:

#. "all_distance_to_VR_max.csv" in the parental folder. This file
   contains the calculated distance to the maximum point of stimulus
   intensity in the virtual reality in mm for each frame
   (rows) for all the analyzed experiments (columns)

#. "Median_Distance_to_VR_max.png" in the parental folder. This plot
   is intended to give a quick overview about the distance to source
   of all the analyzed experiments. Individual experiments are
   plotted in light grey and the median in red.

See here for actual code:
:class:`analysis_scripts.AnalysisVRDistanceToSource`

.. _tools_single_animal_tracking:

Single animal tracking (post-hoc)
---------------------------------

If you have recorded image sequences or videos of single animals
behaving and you wish to track their position, try this option.

While it works best if the PiVR
:ref:`Full Frame Recording<FullFrameLabel>` or the PiVR
:ref:`Video<VideoLabel>` recording option was used, the analysis
program should be able to handle other image sequences and video files
as well.

.. important::
    Make sure that **one experiment/trial is in one folder**. For
    example, if you have recorded three videos, *'Video_1.mp4'*,
    *'Video_2.mp4'* and *'Video_3.mp4'* each video **must** be in its own
    folder (e.g. *'Video_1.mp4'* goes into *'Folder_1'*, *'Video_2.mp4'*
    goes into *'Folder_2'* etc.),for the analysis software to work.

Once you press the 'Press to start analysis' button, the software
will check whether it can find metadata in order to perform the
tracking. If it does not find it, it will ask for user input.
Specifically, it needs:

#. The frame rate the video/image sequence was recorded in.

#. How many pixels are one mm. The software will use the
   :ref:`identical popup as on PiVR<PixMMLabel>` with the only
   difference being that you have to select a file with a known
   distance.

#. An estimate of the maximum animal speed. If unsure, it is better
   to overestimate the speed. If the tracking algorithm does not
   produce the desired result the maximum animal speed should be
   lowered.

#. An estimate of the maximum length of the animal. If unsure, it is
   better to overestimate the length. If the tracking algorithm does
   not produce the desired result, the maximum length should be lowered.

This tracking algorithm is identical to the tracking algorithm used for
live tracking of the animal in PiVR.

The output is therefore almost identical to a tracking experiment run on
PiVR. See :ref:`here<OutputTrackingLabel>` for an explanation.

The only two differences are:

   #. You will find a file called "DATE_Time_heuristics.csv" in the
      analyzed folder. This file contains useful information in order
      for you to define a new animal in "list_of_available_organisms
      .json".

   #. The "DATE_Time" for both the heuristics.csv and the data.csv
      indicate the time of the analysis **NOT** the time of recording
      the data.

The Image Data Handling Menu
=============================

The Image Data Handling Menu lets you choose a folder (either with
one experiment or a folder containing several experiments) and
convert the image data.

.. figure:: Figures/tools/4_ImageDataHandlingFrame.png
   :alt: 4_ImageDataHandlingFrame.png

To convert image data (either a series of full frame images or
videos) first press the "Press to select data to modify" button on
the top left. Select a folder. You may:

#. Select a single experiment/folder. In that case you **must**
   uncheck the button on the right next to "More than one folder".

#. Select a folder containing several experiments. If this folder
   only contains you want to apply the same image data conversion,
   you do not have to do anything. If there are other files or
   folders, please indicate a commonality between all the folder in
   the entry field below. This can be the date (e.g. 2019) a genotype
   (e.g. Or42a) or...

    The "Number of files" below the entry field indicates how many
    files and folders in the folder you selected will be taken into
    account if the current input is used.

Next, you have to choose the image conversion you want to perform
using the dropdown menu under "What modifications do you want to do?"
. Currently there are the following options:

Image Seq. to matrix
--------------------

This is intended to be used after recording a series of images
with the :ref:`Full Frame Recording Option<FullFrameLabel>` option of
PiVR. The disadvantage of having a large number of single image files
instead of one large file (with the identical size) is the time it
takes the PC to read and write (e.g. copy, manipulate etc.) many
single image files. This program will help you to "pack up" your
image files. There are a variety of options you can choose from on
the right side of the GUI:

#. Zip Images: If this checkbox is marked, the images will be zipped.
   You will find a zip file called 'images.zip' in the folder where
   the original images were located. The zip file is uncompressed -
   the goal here is to have one file instead of many single files,
   not to save disk space!

#. Delete Original: If this checkbox is marked, the script will
   delete all images that are considered to be part of the recorded
   image series. After zipping the images, it is useful to delete the
   original images as they will only make data handling slow, but:
   **Make sure to only have this option on if at least one of the
   other options is selected, otherwise your data will be lost!**

#. Greyscale/Color: This dropdown menu lets you choose whether you
   want the images to be saved in greyscale (standard, especially if
   a standard PiVR version is being used) or if the input colors are
   color images and you want the output to be saved as color images.
   *Be careful with the color option, this has not been fully tested
   yet-.*

#. Save \*.npy: If this checkbox is marked, the script will save the
   images in a
   `Numpy array <https://docs.scipy.org/doc/numpy/reference/arrays.html>`_.
   This can be very handy if you want to use python to run the
   downstream analysis of the data.

#. Save \*.mat: If this checkbox is marked, the script will save the
   images in a
   `Matlab like array <https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.savemat.html>`_.
   This can be very handy if you want to use Matlab to run downstream
   analysis of the data.

.. _h264_to_avi_GUI:

Video conversion
----------------

This is intended to be used after recording a video with the
:ref:`Video<VideoLabel>` option of PiVR. Many users will find the h264
video file straight from PiVR inconvenient to work with, since: (1) Many
video players (e.g. VLC) have problems decoding the video
and (2) the metadata of the video seem to not always be correct.

.. note::
   I found that the metadata of the videos recorded on the Raspberry
   Pi are not completely correct. For example, when reading a video
   file with the imageio module (using ffmpeg) the number of frames
   is given as "inf" and the framerate seems to be always at 25, even
   though the video was recorded at a different framerate. This
   script takes care of this bug by using the 'experiment_settings
   .json' file created when using PiVR to record a video.

There are several options available to define the desired output:

#. avi/mp4/None: If you want to watch the video, it is probably a
   good idea to convert the video either into avi or mp4 as your
   video player will be better able to handle these formats (and the
   metadata of this file will be correct, see above). If "None" is
   chosen, the video will *not* be converted.

#. h264/raw video: Besides the format, the codec can be a problem for
   some programs you want to open a video. Here, you can choose
   between the efficient h264 codec which will lead to very small
   file sizes and the raw video codec. The latter will lead to
   significantly larger video files! **h264 is recommended!**

#. Greyscale/Color: This dropdown menu, lets you choose whether you
   want the video to be in color (only works if original video is in
   color, of course) or converted to greyscale.

#. Save \*.npy: If this checkbox is marked, the video will save each
   frame in a `Numpy array <https://docs.scipy.org/doc/numpy/reference/arrays.html>`_.
   This can be very handy if you want to use python to run the
   downstream analysis of the data.

    .. Note::
       Video encoding has been perfected over the years.
       Decompressing a video often leads to surprisingly large files,
       especially for long videos, or videos with a high framerate. If
       the uncompressed video is larger than your computer has RAM this
       script will most likely fail.

#. Save \*.mat: If this checkbox is marked, the script will save the
   images in a
   `Matlab like array <https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.savemat.html>`_.
   This can be very handy if you want to use Matlab to run downstream
   analysis of the data.

    .. Note::
       Video encoding has been perfected over the years.
       Decompressing a video often leads to surprisingly large files,
       especially for long videos, or videos with a high framerate. If
       the uncompressed video is larger than your computer has RAM this
       script will most likely fail.

.. _undistort_h264_GUI:

Undistort Video
---------------

Videos recorded with PiVR have obvious lens distortions which introduce
a mild fisheye effect.

.. figure:: Figures/tools/4_1distorted.png
   :alt: 4_1distorted.png

Using this option corrects this artifact using openCV functions. An
in-depth explanation can be found `on their website
<https://docs.opencv.org/master/dc/dbb/tutorial_py_calibration.html>`_.

After running the algorithm on the image shown above, the lines
appear significantly straighter now:

.. figure:: Figures/tools/4_2undistorted.png
   :alt: 4_2undistorted.png

.. Important::
   The algorithm needs two matrices, *dist.npy* and *mtx.npy* which are
   provided for the standard Raspberry Pi camera lens. If you are
   using a different lens, you must redefine these matrices. Please
   see `here <https://docs.opencv.org/4.5.2/dc/dbb/tutorial_py_calibration.html>`__
   or `here <https://gitlab.com/davidtadres/cameracalibrations>`__ to
   learn how to update them.

There are several options available to define the desired output:

#. avi/mp4/None: If you want to watch the video, it is probably a
   good idea to convert the video either into avi or mp4 as your
   video player will be better able to handle these formats (and the
   metadata of this file will be correct, see above). If "None" is
   chosen, the video will *not* be converted.

#. h264/raw video: Besides the format, the codec can be a problem for
   some programs you want to open a video. Here, you can choose
   between the efficient h264 codec which will lead to very small
   file sizes and the raw video codec. The latter will lead to
   significantly larger video files! **h264 is recommended!**

#. Save \*.npy: If this checkbox is marked, the video will save each
   frame in a `Numpy array <https://docs.scipy.org/doc/numpy/reference/arrays.html>`_.
   This can be very handy if you want to use python to run the
   downstream analysis of the data.

    .. Note::
       Video encoding has been perfected over the years.
       Decompressing a video often leads to surprisingly large files,
       especially for long videos, or videos with a high framerate. If
       the uncompressed video is larger than your computer has RAM this
       script will most likely fail.

#. Save \*.mat: If this checkbox is marked, the script will save the
   images in a
   `Matlab like array <https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.savemat.html>`_.
   This can be very handy if you want to use Matlab to run downstream
   analysis of the data.

    .. Note::
       Video encoding has been perfected over the years.
       Decompressing a video often leads to surprisingly large files,
       especially for long videos, or videos with a high framerate. If
       the uncompressed video is larger than your computer has RAM this
       script will most likely fail.

Display tracked experiment
==========================

This tool allows you to display the tracked animal in its arena
similar to a video player.

The experiment must have been generated using PiVR: either on on the
Raspberry Pi using the :ref:`Tracking<TrackingLabel>` or
:ref:`Virtual reality<VRLabel>` or by using the
:ref:`post-hoc analysis<tools_single_animal_tracking>`
option in the 'Tools' Menu.

.. figure:: Figures/tools/5_DisplayTrackedExperiment.png
   :alt: 5_DisplayTrackedExperiment.png

Besides enabling you to conveniently see where the animal was in each
frame, this tool allows you to:

#. Correct false Head/Tail assignments by swapping them

#. Save a video (in mp4 format) of the experiment.

To display an experiment, press on the "Select data to analyze"
button and select a folder containing a single experiment.

.. figure:: Figures/tools/5_1_ExprimentSelected.png
   :alt: 5_1_ExprimentSelected.png

The software will automatically read the metadata of the experiment
and will be displayed on the left side of the GUI.

In the center, the 'Overview of tracking.png' file is shown.

On the right you may choose a colormap before pressing "Press to show
behavior of the animal".

.. figure:: Figures/tools/5_2_ExperimentDisplay.png
   :alt: 5_2_ExperimentDisplay.png

The window above will emerge as a popup. While it is open the main
GUI is unavailable!

This window has the following functionality:

#. Updating Overview: This button will allow you to turn the main
   window off when playing back the experiment. This can be useful if
   you are only interested in the small images on the right.

#. The main figure in the center of the window is created by placing
   the raw small image (sm_raw.npy) into the reconstructed background
   image (Background.jpg) using the bounding box coordinates
   (bounding_boxes.npy). It also displays the detected Centroid, Tail
   and Head (\*_data.csv) directly on the animal.

#. The three buttons below the main window on the left, "Showing
   Centroid", "Showing Head" and "Showing Tail" allow you to turn off
   the different parts shown in the main figure.

#. The toolbar below allows you to interact with the main window
   (zoom etc.).

#. The slider below lets you scroll through the experiment. Pressing
   the "Start Playing" button will display the experiment.

#. The Dropdown menu below the slider lets you play back the
   experiment at a variety of speeds you can choose from.

#. If there is a particular frame you want to go to, you can enter
   that number in the field below and press "Jump to frame".

#. If the head/tail assignment has been made incorrectly, you have two
   options to fix it here.

   #. *Manual Head Tail correction*: Define the timepoints in which
      the head and tail have been incorrectly identified in 'Start
      frame' and 'End frame. Press 'Swap Head & Tail'

   #. *Automatic Head Tail correction*: You need to slide the slider
      to a frame with incorrect head tail classification. Then press
      'Automatically swap Head Tail'. The algorithm will swap head
      and tail for previous and following frames until it finds a
      frame where the head and tail could not be assigned, for
      example during a U-Turn.

#. On the top right, the small raw image (sm_raw.npy) is being
   displayed.

#. Below, the binary small image (sm_thresh.npy) is being displayed.

Creation of a video file
------------------------

In the "Save as Video" box at the bottom right you may create a
video of the experiment. There are several options to customize the
output file:

.. note::
   If you get an error message it is possible that you have not
   installed the correct version of ffmpeg. Usually the error message
   should indicate what you have to install. If you need help, please
   open a `ticket <https://gitlab.com/LouisLab/pivr/-/issues>`__.

#. If you want to create a video displaying the the animal at the
   correct position, all you need to do is define the start and end
   frame of the video and hit 'Save Video'. If you do not want the
   head, tail and/or centroid to be displayed in the video, turn it
   off as you would for the interactive visualization.

#. If you have been stimulating your animal with a time dependent
   stimulus and you want your video to indicate when the stimulus was
   turned on, press the button titled 'select stim file' and select
   your stimulus file. The stimulus file must adhere to the
   :ref:`PiVR guidlines <PrepTimeStimLabel>`. Then hit 'Save
   Video'

   .. figure:: Figures/tools/5_3_time_dependent_stimulus.png
      :alt: 5_3_time_dependent_stimulus.png

#. If you have been stimulating your animal with a static Virtual
   Arena and you want your video to indicate where the animal was
   relative to the position of the Virtual Arena. Press the button
   titled 'select VR arena' and select your VR file. Then hit 'Save
   Video'

   .. figure:: Figures/tools/5_4_static_VR_stimulus.png
      :alt: 5_4_static_VR_stimulus.png

#. If you have been presenting your animal with a dynamic virtual
   reality (i.e. one that changes over time) you may also include
   this in the video. Press the button titled 'select VR arena' and
   select the file you used when presenting the animal with the dynamic
   virtual reality. Finally, indicate in 'Update arena every nth
   frame' after how many frames the arena needs to be updated. For
   example, if you have been presenting an arena that updates at 15Hz
   while recording at 30Hz you need to update every 2nd frame
   (30/15=2).

   .. figure:: Figures/tools/5_5_dynamic_VR_stimulus.png
      :alt: 5_5_dynamic_VR_stimulus.png

When pressing "Save video" the video will be saved as "Video.mp4" in
the experimental folder. **This usually takes a significant amount of
time, even on a fast computer**.

.. _ToolsMultianimalTracking:

Multi-Animal Tracking
======================

This tool allows you to identify more than one animal in a video or
in a series of images.

The identification of the animals is achieved via background
subtraction of the mean image. Each animal is given an arbitrary
number. In the next frame, the animal closest to that number in the
previous frame is assigned.

This guide is intended to get you started quickly. If you are
interested in the more technical aspects, please see
:py:class:`multi_animal_tracking.MultiAnimalTracking`

.. note::
   Identifying more than one animal in an experiment is
   computationally challenging. There are several specialized tools
   such as the `multiworm tracker <https://omictools.com/mwt-tool>`_,
   `Ctrax <http://ctrax.sourceforge.net/>`_,
   `idtracker <http://www.idtracker.es/>`_,
   `idtracker.ai <http://idtracker.ai/en/latest/>`_ and
   `MARGO <https://github.com/de-Bivort-Lab/margo>`_. The PiVR
   multi-animal tracking software has **not** been benchmarked
   against these tools. This software has several limitations. It is
   probably ok to use in cases where you are interested in counting
   how many animals are in a general area. It is not recommend to use
   the tracker for other parameters, such as calculating animal speed
   (due to loss of identity after collision and 'jumps' in the
   trajectory) and similar parameters.

.. figure:: Figures/tools/6_MultiAnimalTracker.png
   :alt: 6_MultiAnimalTracker.png

After selecting a folder containing a video of image sequence of an
experiment containing multiple animals, press the "Press to show the
behavior of the animals" button.

The software will now load the image data (which can take a
considerable amount of time) and display a new window which will help
you to optimize the tracking parameters.

.. figure:: Figures/tools/6_1_MultiAnimalTrackerWindow.png
   :alt: 6_1_MultiAnimalTrackerWindow.png

.. note::
   You might notice that the image is distorted. This is due to the
   Raspberry Pi camera lens. See :ref:`here <FAQ_undistort>`.

First, have a look at the all-important "# of animals in the current
frame" on the bottom right of the popup. In this video there are only
6 animals, but the algorithm detects 10 objects as "animals". The
blobs identified as animals are indicated using small rectangles in
the main figure.

.. important::
   The goal is to have have as many frames as possible contain only the
   expected amount of animals.

To achieve this goal, start by defining the area in which animals can
move into. In this particular experiment, the outline of the petri
dish can be clearly seen:

#. Press the "Select Rectangular ROI". The button will turn red.

#. Using the mouse, create a rectangle in the main window. The result
   is indicated. The "# of animals in the current frame" is
   immediately updated.

.. figure:: Figures/tools/6_2_MultiAnimalTrackerSelectRectangle.png
   :alt: 6_2_MultiAnimalTrackerSelectRectangle.png

While this is already better, there are still two blobs wrongly
identified as animals, one at (y=300,x=100) and the other at (y=270,
x=480).

By increasing the "Threshold (STD from Mean) number these
mis-identified blobs are not mistaken for animals animal:

.. figure:: Figures/tools/6_3_MultiAnimalTrackerSelectTreshold.png
   :alt: 6_3_MultiAnimalTrackerSelectTreshold.png

There are now no mis-identified animals. While there are a total of 6
animals in the image, the algorithm can only detect 5. To
understand why that is the case, press the magnifying class symbol
below the main figure and draw a rectangle around the region you want
to enlarge.

.. figure:: Figures/tools/6_4_MultiAnimalTrackerSelectZoom.png
   :alt: 6_4_MultiAnimalTrackerSelectZoom.png

It is now obvious that 4 out of 6 animals are properly identified.
Two of the animals, at (y=260, x=360) are very close to each other
and are identified as one animals, however. *There is no way for this
tracking algorithm to separate animals that are touching each other!*

Now you need to make sure that the image parameters are such that you
get the expected animal number in all frames. To not have to go
through each of the frames manually, you can just press "Auto-detect
blobs" on the bottom right. This will run a script that will take the
current image detection parameters into account and just count how
many blobs are counted as animals. The result is plotted in the
figure on the top right.

.. figure:: Figures/tools/6_5_MultiAnimalTrackerSelectDetectBlobs.png
   :alt: 6_5_MultiAnimalTrackerSelectDetectBlobs.png

This result indicates that at the beginning of the video there are
several frames where only 5 animals can be detected. By visually
inspecting these frames it becomes obvious that this is due to the
two animals touching each other as they already do in the first frame.

Next, find frames that have the wrong amount of animals. Try to fix
them using the image parameter settings.

.. figure:: Figures/tools/6_6_MultiAnimalTrackerSelectSearchProblematicFrames1.png
   :alt: 6_6_MultiAnimalTrackerSelectSearchProblematicFrames1.png

There can be situations where the animal number will just be wrong
and it can not be fixed. The algorithm can handle this if the number
of those frames is low.

Once you have optimized the image parameters, go to the frame where the
correct amount of animals is detected. **This is crucial as it tell
the algorithm how many animals to expect!**. Then press "Track animals".

A new popup will open. It indicates how many animals (and where) are
detected in this frame. Each animal gets an arbitrary number. If the
number of animals is correct, press "Looks good, start tracking".
Else press "Not good, take me back".

.. figure:: Figures/tools/6_7_MultiAnimalTrackerStartTracking.png
   :alt: 6_7_MultiAnimalTrackerStartTracking.png

The tracking itself is computationally quite expensive and therefore
usually takes a while to complete. To speed up tracking, you can
press the "Updating Overview" button above the main window.

Once tracking has concluded the result will be displayed in the main
window as shown below.

.. figure:: Figures/tools/6_8_MultiAnimalTrackerTrackingResult.png
   :alt: 6_8_MultiAnimalTrackerTrackingResult.png

If there are huge gaps in the trajectory, for example because an
animal could not be detected for a while, you can press the
"interpolate centroids" button. This will calculate realistic
possible distances travelled (based on 'Max Speed Animal[mm/s])
between frames and try to connect trajectories. This is a
**untested** feature - use at your own risk. Ideally you should not
have to use this option.

Multi Animal Tracking Output
-----------------------------

After using the Multi-Animal Tracker, you will find two new files in
the experimental folder:

#. The "Background.npy" file which is just the mean of all images,
   resulting in the background image used during tracking.

#. The "XY_positions.csv" file contains the X and Y coordinates for
   each identified animal for each frame. For frames with not enough
   animals, the corresponding row will be empty.

.. _create_VR_arena:

Creating a new Virtual Arena
=============================

How to create a new virtual arena, the essential tool that gives
Pi\ **VR** its name?

First, a brief overview of what the different elements of an arena
file mean:

You can find several example virtual arena files in the
PiVR/VR_arenas folder: As an example lets analyze
"640x480_checkerboard.csv"

The file itself is a 2D matrix with 640 columns and 480 rows. Each
value in the matrix defines what happens on that pixel on the camera.

For example, if you define the value 100 at position column = 75
and row = 90 here in the virtual arena and then present this virtual
arena to the animal using the :ref:`Closed Loop stimulation<VRLabel>`
tab on PiVR, if the animal is at pixel column = 75 and row = 90 the
intensity of 100 will be played back.

A big advantage of virtual realities over real environments is the
control the experimenter has over the experimental conditions. When
running an experiment, one often has to repeat trials many times.
In real environments, the experimenter can never have identical
initial conditions, e.g. the animals is placed at a slightly
different position, the animal moves in different conditions before
the tracking even starts etc. With virtual reality, this factor
(which often introduces variability into data) can be alleviated. The
experimenter can define a virtual reality arena and the animal will
always be presented with the identical initial conditions.

To do this, lets examine another example virtual arena file you can
find in PiVR/VR_arenas: "640x480_gaussian_centred_animal_pos[250,
240,0.0].csv"

This file has a string at the end of its filename: *animal_pos[250,
240,0.0]*. This string indicates where the animal must start in
relation to the virtual reality. In this example, wherever the animal
is in the real image when the experiment starts, the virtual reality
will be translated so that it is at x-coordinate 250 and y coordinate
240.

In addition, if the third value (here 0.0) is defined, the movement
of the animal during detection is taken into account: The virtual
arena is rotated so that the animal always starts going into the same
direction relative to the virtual arena. The angle you may use goes
from -pi to +pi (see `atan2 <https://en.wikipedia.org/wiki/Atan2>`__).

.. figure:: Figures/tools/7_explain_names.png
   :alt: 7_explain_names.png

If you want to create a virtual reality arena from scratch, for
example in python, all you need to do is create a matrix with the
correct dimension, fill it with values between 0 and 100 as you see
fit and export the file as csv.

.. code-block:: python

   import numpy as np

   virtual_arena = np.zeros((480,640),dtype=np.float64)
   # Define parts of the arena where the animal is supposed to be
   # stimulated e.g. by typing
   virtual_arena[:,0:100] = 100
   virtual_arena[:,101:200] = 75
   virtual_arena[:,201:300] = 50
   # this will give you a very coarse grained virtual arena that will
   # stimulate strongly if the animal is on the left side, stimulate
   # 75% if the animal is a bit more on the right and 50% if it is
   # still on the left but almost in the middle. The rest is
   # unstimulated as of yet.

   # Now you need to save the arena. Let say you want to have the animal
   # start ascending the gradient from the middle (essentially animal
   # has to move to the left in the virtual arena)
   np.savetext("Path/On/Your/PC/640x480_my_awesome_arena_animal_pos[300,240,0.0].csv")

Alternatively, you can use the "Tools->Draw VR Arena" option on the
:ref:`PC version of PiVR<software_installation_PC>`.

.. figure:: Figures/tools/8_DrawVRArena_menubar.png
   :alt: 8_DrawVRArena_menubar.png

You will find the following empty canvas. You can open a previously
defined virtual arena or work on the blank canvas. Either way, you
have the the option to create *gaussian shaped 2D circles* and
*rectangles*. To "draw" such a gaussian shaped 2D circle, you can
either press on "Draw Gaussian Circle with Mouse" (and then click
somewhere on the canvas) or you can press on "Draw Gaussian Circle at
defined coordinates".

.. figure:: Figures/tools/9_DrawVRArena_overview.png
   :alt: 9_DrawVRArena_overview.png

You can change the Gaussian shaped 2D circle by changing its Sigma,
its size and the intensity.

Analogous, you can define the size and the intensity of the rectangle
by entering the desired value.

In the example below, I used the standard settings for the gaussian
circle size but changed the "Coordinates" to the values shown. Then I
pressed on "Draw Gaussian Circle at defined coordinates"(red squares).
Then I modified the "Coordinates" of the rectangle on the right to
x=100 and y=100 and pressed on "Draw Rectangle at defined
coordinates" (green squares).

.. figure:: Figures/tools/9_1_DrawVRArena_example_1.png
   :alt: 9_1_DrawVRArena_example_1.png

There are 3 additional buttons that expand the possibilities of drawing
virtual arenas:

#. Invert: If this is on, whenever you draw a circle or rectangle, it
   subtracts the values from the intensity that is already present.
   This option was used to create the following virtual arena:
   "640x480_volcano_animal_pos[320,240,0.0].csv"

#. Overwrite: will just overwrite the previous pixel values with the
   new values with no regard to the previous value (as opposed to
   "Invert" and "Additive")

#. Additive: If you place a 50% rectangle somewhere and then place
   another on top of it, usually nothing will happen as the absolute
   value is being drawn. If this is on, the values are "added".

   .. note::
      If you go above 100%, everything will be normalized.

On the right side of the canvas you can define the starting position
of the animal in the virtual arena. Besides just x and y position,
you can define the direction from which the animal is coming from.

You can use the mouse, either by clicking (just x and y position) or by
clicking and dragging (x, y and angle). Alternatively, you can use
precise coordinates.

.. note::
   Angle can go from -pi to +pi. See
   `atan2 <https://en.wikipedia.org/wiki/Atan2>`__ for visualization.

.. figure:: Figures/tools/9_2_DrawVRArena_example_2.png
   :alt: 9_2_DrawVRArena_example_2.png

Once you are done with the arena, make sure to save it. Then you can
just quit the window.

.. _create_dynamic_VR_arena:

Creating a new Dynamic Virtual Reality
======================================

Before attempting to create a new dynamic virtual reality please make
sure you understand exactly how a
:ref:`static virtual reality is created <create_VR_arena>`.

The main differences between the static and the virtual reality are:

#. Usage of `numpy arrays <https://docs.scipy
   .org/doc/numpy/reference/arrays.ndarray.html>`__ instead of csv
   files.

#. The filename includes the speed at which the virtual reality is
   "played back" as a "Hz[#]" term.

   .. figure:: Figures/tools/10_1_dynamic_arena_name.png
      :alt: 10_1_dynamic_arena_name.png

#. Instead of going from 0 - 100, stimulation in defined as 8 bit (0
   - 255)

#. For details on how exactly this is implemented, please see
   :ref:`here <Dynamic Virtual Realities explained>`

There are several limitations compared to static virtual realities:

#. Dynamic virtual realities can **not** orient themselves
   relative to the the position of the animal at the start of the
   experiment.

#. RAM of PiVR is limited (but this is getting better with the RPi4
   which has up to 4Gb). As the arena that is presented must be
   loaded into RAM at the beginning of the experiment, it can not be
   too large. On our RPi3+ (1Gb of RAM) we successfully used dynamic
   virtual realities of 135Mb for 5 minute experiments.

#. It is not possible to use the "Adjust Power" field!

How does a dynamic virtual reality file look like? Essentially it is
just the 3D version of the static virtual arena with the difference
that an **integer value of 255 indicates that the LED is completely
on** (as opposed to 100 in the static case). The third axis contains
all the 'frames' of the virtual reality 'video'.

As an example, follow the instructions below to convert a real odor
plume measurement to a dynamic virtual reality that PiVR can use:

#. Get the odor plume data from `Álvarez-Salvado et al
   ., 2018 <https://datadryad.org/stash/dataset/doi:10.5061/dryad
   .g27mq71>`__ and unzip the file. The relevant file is
   "10302017_10cms_bounded_2.h5".

#. While downloading, you can make sure you have the correct conda
   environment to handle the data. Below is the code to create new
   environment and install the necessary packages (Windows)::

      conda create -n dynamic_odor_env python=3.7 -y
      activate dynamic_odor_env
      conda install -c anaconda h5py -y
      conda install matplotlib -y

#. Now start python and adapt the following code::

      import numpy as np
      import h5py

      # Adapt this to the location of your file
      data_path = 'C:\\10302017_10cms_bounded_2.h5'

      # Read the h5 datafile using the h5py library
      data = h5py.File(data_path, 'r')

      # The h5 datafile is a dictionary. As the keys are unkown,
      # cycle through to learn what files are stored in this fie.
      for key in data.keys():
         print(key) #Names of the groups in HDF5 file.
      # ok, the relevant (only) key is 'dataset2'
      # Using this group key we can now access the data
      actual_data_temp = data['dataset2']
      actual_data = actual_data[()]
      # make sure you do have the data:
      print(actual_data.shape)
      # (3600, 406, 216)
      # This information, together with the metadate indicates that
      # the first dimension holds the measurements of the timepoints
      # (Frequency of 15Hz for a duration of 4 minutes = 3600 frames).
      # The other two dimensions are spatial dimensions.

      # Below the reason it is necessary to use uint8 (8bit integer
      # precision instead of float (just remove the #)
      # import sys
      # print('Using 64bit float precision: ' +
      #    repr(sys.getsizeof(np.zeros((actual_data.shape[0],640,480),
      #    dtype=np.float64))/1e6) + 'Megabytes' )
      # print('Better to use uint8 as I ll only need to use: ' +
      #    repr(sys.getsizeof(np.zeros((actual_data.shape[0],640, 480),
      #    dtype=np.uint8))/1e6) + 'Megabytes in RAM' )

      # Now to convert the 64bit floating point to 8 bit
      # find maximum value:
      max_value = np.amax(actual_data)
      # find minimum value:
      min_value = np.amin(actual_data)

      # convert data to uint8 to save a ton of space
      convert_factor = 255/(max_value + (- min_value))
      downsampled = (actual_data+(-min_value)) * convert_factor

       # We need to create an empty arena with 640x480 pixels using the
       # correct datatype..
       correct_size_arena = np.zeros((480,640,actual_data.shape[0]),
          dtype=np.uint8)
       # and place the moving arena in to the center.
       start_y = int(correct_size_arena.shape[0]/2 -
          downsampled.T.shape[0]/2) # 132.0
       start_x = int((correct_size_arena.shape[1]/2 -
           downsampled.T.shape[1]/2)) # 117
       correct_size_arena[start_y:int(start_y+downsampled.T.shape[0]),
          start_x:int(start_x+downsampled.T.shape[1]),:] = downsampled.T

       # In order to see how the just created arena looks like, run
       # the code below (uncomment by removeing the #)
       # import matplotlib.pyplot as plt
       # fig,ax = plt.subplots()
       # ax.imshow(correct_size_arena[:, :, 0])
       # plt.show()

       # Even after downsampling, file is way to large (>1Gb). To
       # solve this, only save 1/8 of the file (~30seconds).
       np.save('C:\\dynamic_odor_plume_eigth_Hz[15].npy',
          correct_size_arena[:,:,0:int(correct_size_arena.shape[2]/8)])

#. Now you must transfer the file containing the new dynamic virtual
   reality to your PiVR setup and select it when running your
   experiment, identical to what you would do to
   :ref:`select a static virtual arena <VRLabel>`.





PiVR has been developed by David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

.. _Output:

Explanation of PiVR output
**************************

.. _OutputTrackingLabel:

Tracking
========

After running a :ref:`tracking <TrackingLabel>` experiment you will
find a folder with the "DATE_TIME_EXP.GROUP" as its name. An example
would be "2019.01.11_14-00-05_CantonS". This is an experiment
conducted on the 11th of January 2019. "CantonS" is the value that
was entered in the field "Exp. Group".

This folder will contain the following files:

"DATE_TIME_data.csv"
--------------------
is probably the most important file. It contains the following data
for each frame of the experiment:

   #. The frame (=image) number into the experiment
   #. The time in seconds since the experiment started
   #. The X (column) coordinate of the Centroid (Check
      :ref:`here <OutputPointExp>` for comparison with midpoint)
   #. The Y (row) coordinate of the Centroid
   #. The X (column) coordinate of the head
   #. The Y (row) coordinate of the head
   #. The X (column) coordinate of the tail
   #. The Y (row) coordinate of the tail
   #. The X (column) coordinate of the midpoint (Check
      :ref:`here <OutputPointExp>` for comparison with centroid)
   #. The Y (row) coordinate of the midpoint
   #. The Y-min (row) coordinate of the bounding box (See
      :ref:`here <bbox_label>` for explanation
   #. The Y-max (row) coordinate of the bounding box.
   #. The X-min (row) coordinate of the bounding box.
   #. The X-max (row) coordinate of the bounding box.
   #. The local threshold used to extract the binary image during tracking.

"Background.jpg"
-----------------
contains the reconstructed background image. See
:ref:`here <AnimalDetectionExplanationLabel>` for explanation
where it is coming from and what it means.


"experiment_settings.json"
--------------------------
is a `json file <https://en.wikipedia.org/wiki/JSON>`_ and contains a
lot of useful experimental information:

#. Camera Shutter Speed [us]: Shutter speed in microseconds

#. Exp. Group: The string that was entered by the user during the
   experiment

#. Experiment Date and Time: exactly that

#. Framerate: The frequency at which PiVR tracked the animal

#. Model Organism: While tracking, PiVR used the parameters of this
   animal to optimize tracking. See :ref:`here<define_new_animal_label>`
   for how to modify this parameter.

#. PiVR info (recording): version number, git branch and git hash of
   the PiVR software that was used to record the experiment.

#. PiVR info (tracking): version number, git branch and git hash of
   the PiVR software that was used to track the experiment. If online
   tracking is being done, this is identical to the info above.

#. Pixel per mm: For PiVR to be able to track the animal, it needs
   to know how many pixels indicate one mm. This has been set by the
   user as described :ref:`here<PixMMLabel>`.

#. Recording time: The time in seconds that PiVR was tracking the
   animal

#. Resolution: The camera resolution in pixels that PiVR used while
   tracking.

#. Time delay due to Animal Detection[s]: For the
   :ref:`autodetection<AnimalDetectionExplanationLabel>` the animal
   must move. The time it took between pressing "start" and successful
   animal detection is saved here.

#. Virtual Reality arena name: If no virtual arena was presented,
   it will say 'None'

#. backlight 2 channel: If Backlight 2 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 18)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. backlight channel: If Backlight 1 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 18)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].
   This would normally be defined as [18, 40000].

#. output channel 1: If Channel 1 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 17)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. output channel 2: If Channel 2 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 27)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. output channel 3: If Channel 3 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 13)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. output channel 4: If Channel 4 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 13)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

"first_frame_data.json"
-----------------------
is a `json file <https://en.wikipedia.org/wiki/JSON>`__ and contains
information that collected during
:ref:`animal detection<AnimalDetectionExplanationLabel>` (Source
code :class:`pre_experiment.FindAnimal`.)

#. bounding box col max: The X_max value of the bounding box of the
   animal detected in the first frame during animal detection.

#. bounding box col min: The X_min value of the bounding box of the
   animal detected in the first frame during animal detection.

#. bounding box row max: The Y_min value of the bounding box of the
   animal detected in the first frame during animal detection.

#. bounding box row min: The Y_max value of the bounding box of the
   animal detected in the first frame during animal detection.

#. centroid col: The X value of the centroid of the animal detected
   in the first frame during animal detection.

#. centroid row: The Y value of the centroid of the animal detected
   in the first frame during animal detection.

#. filled area: The filled area in pixels of the blob defined as
   the animal in the first frame during animal detection

Optional files
--------------

The rest of this subsection describes files that need to be explicitly
saved as described here. TODO LINK

.. _bbox_label:

"bounding_boxes.npy"
~~~~~~~~~~~~~~~~~~~~~~
is a `Numpy <https://numpy.org/>`_ file. It contains the coordinates
of the bounding box of the small image. The bounding box defines the
Y/X coordinates of the small image.

This file comes in
`shape <https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.shape.html>`_
[4, # of frames] with:

======  =========================
[0, :]  contains the Y_min values
------  -------------------------
[1, :]  contains the Y_max values
------  -------------------------
[2, :]  contains the X_min values
------  -------------------------
[3, :]  contains the X_max values
======  =========================

These values are necessary to describe where in the full image
frame the small image that has been saved during the experiment is
located. The bounding box is the rectangle that contains all image
information used during this frame. Below an illustration on how
the different values are used to construct the bounding box.

.. figure:: Figures/ouput/bbox_illustration.png
  :width: 100 %
  :alt: bbox_illustration.png

.. note::

  Why Y/X and not X/Y? In image processing the convention is to
  reference points in (Rows, Columns) which translates to Y/X. The
  underlying image processing libraries work with the (Rows,
  Columns) convention. See for example `here
  <https://scikit-image.org/docs/dev/api/skimage.measure.html#skimage
  .measure.regionprops>`__.
  PiVR therefore follows this convention.

"centroids.npy"
~~~~~~~~~~~~~~~~~~~~~~
is a `Numpy <https://numpy.org/>`__ file. It contains the coordinates
of the centroid of the blob identified during the experiment.
See :ref:`here<OutputPointExp>` to see the centroid compared to the
midpoint.

The file comes in `shape <https://docs.scipy
.org/doc/numpy/reference/generated/numpy.ndarray.shape.html>`__
[# of frames, 2] with:

======  ==============================
[:, 0]  contains the centroid Y values
------  ------------------------------
[:, 1]  contains the centroid X values
======  ==============================

These values are identical to what you will find in the
"DATE_TIME_data.csv" file

"midpoints.npy"
~~~~~~~~~~~~~~~~~~~~~~
is a `Numpy <https://numpy.org/>`__ file. It contains the coordinates
of the midpoint extracted from the skeleton during the experiment.
See :ref:`here<OutputPointExp>` to see the midpoint compared to the
centroid.

The file comes in `shape <https://docs.scipy
.org/doc/numpy/reference/generated/numpy.ndarray.shape.html>`__
[# of frames, 2] with:

======  ==============================
[:, 0]  contains the midpoint Y values
------  ------------------------------
[:, 1]  contains the midpoint X values
======  ==============================

These values are identical to what you will find in the
"DATE_TIME_data.csv" file


"heads.npy"
~~~~~~~~~~~~~~~~~~~~~~
is a `Numpy <https://numpy.org/>`_ file. It contains the coordinates
of the head position :ref:`assigned during tracking<HeadTailLabel>`.

The file comes in `shape <https://docs.scipy
.org/doc/numpy/reference/generated/numpy.ndarray.shape.html>`__
[# of frames, 2] with:

======  ==========================
[:, 0]  contains the head Y values
------  --------------------------
[:, 1]  contains the head X values
======  ==========================

These values are identical to what you will find in the
"DATE_TIME_data.csv" file

"tails.npy"
~~~~~~~~~~~~~~~~~~~~~~
is a `Numpy <https://numpy.org/>`_ file. It contains the coordinates
of the tail position :ref:`assigned during tracking<HeadTailLabel>`.

The file comes in `shape <https://docs.scipy
.org/doc/numpy/reference/generated/numpy.ndarray.shape.html>`__
[# of frames, 2] with:

======  ==========================
[:, 0]  contains the tail Y values
------  --------------------------
[:, 1]  contains the tail X values
======  ==========================

These values are identical to what you will find in the
"DATE_TIME_data.csv" file

.. _OutputVRArenaLabel:

VR Arena and Dynamic VR Arena
=============================

After running a :ref:`VR Arena<VRLabel>` experiment you will
find a folder with the "DATE_TIME_EXP.GROUP" as its name. An example
would be "2019.01.11_14-00-05_CantonS". This is an experiment
conducted on the 11th of January 2019. "CantonS" is the value that
was entered in the field "Exp. Group".


This folder will contain the following files:

"DATE_TIME_data.csv"
--------------------
is probably the most important file. It contains the following data
for each frame of the experiment:

   #. The frame (=image) number into the experiment
   #. The time since the experiment started
   #. The X (column) coordinate of the Centroid (Check
      :ref:`here <OutputPointExp>` for comparison with midpoint)
   #. The Y (row) coordinate of the Centroid
   #. The X (column) coordinate of the head
   #. The Y (row) coordinate of the head
   #. The X (column) coordinate of the tail
   #. The Y (row) coordinate of the tail
   #. The X (column) coordinate of the midpoint (Check
      :ref:`here <OutputPointExp>` for comparison with centroid)
   #. The Y (row) coordinate of the midpoint
   #. The Y (row) coordinate of the midpoint
   #. The Y-min (row) coordinate of the bounding box (See
      :ref:`here <bbox_label>` for explanation
   #. The Y-max (row) coordinate of the bounding box.
   #. The X-min (row) coordinate of the bounding box.
   #. The X-max (row) coordinate of the bounding box.
   #. The local threshold used to extract the binary image during tracking.
   #. The stimulus delivered stimulus in %

"RESOLUTION_NAME.csv"
---------------------
for example "640x480_checkerboard.csv". This is the virtual arena
presented to the animal. In case the virtual arena is positioned
relative to the starting position and the movement of the animal
(such as the "640x480_gaussian_centred_animal_pos[250,240,0.0].csv"
arena), this file will *final* translated and rotated arena as it was
presented to the animal.

.. note::
   If a dynamic virtual reality has been presented, this file will
   not be present - it would simply take too long and take up too
   much space. This is one reason why dynamic virtual realities can
   not be translated and rotated at the moment.

"Background.jpg"
-----------------
contains the reconstructed background image. See
:ref:`here <AnimalDetectionExplanationLabel>` for explanation
where it is coming from and what it means.

"experiment_settings.json"
--------------------------
is a `json file <https://en.wikipedia.org/wiki/JSON>`_ and contains a
lot of useful experimental information:

#. Search box size: The Search box used to locate the animal during
   the experiment

#. Exp. Group: The string that was entered by the user during the
   experiment

#. Experiment Date and Time: exactly that

#. Framerate: The frequency at which PiVR tracked the animal

#. Model Organism: While tracking, PiVR used the parameters of this
   animal to optimize tracking. See **Todo** here for how to modify
   this parameter.

#. Pixel per mm: For PiVR to be able to track the animal, it needs
   to know how many pixels indicate one mm. This has been set by the
   user as described :ref:`here<PixMMLabel>`.

#. Recording time: The time in seconds that PiVR was tracking the
   animal

#. Resolution: The camera resolution in pixels that PiVR used while
   tracking. Currently only 640x480 is possible.

#. Time delay due to Animal Detection[s]: For the
   :ref:`autodetection<AnimalDetectionExplanationLabel>` the animal
   must move. The time it took between pressing "start" and successful
   animal detection is saved here.

#. Virtual Reality arena name: As no virtual arena was presented,
   it will say 'None'

#. backlight 2 channel: If Backlight 2 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 18)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. backlight channel: If Backlight 1 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 18)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].
   This would normally be defined as [18, 40000].

#. output channel 1: If Channel 1 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 17)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. output channel 2: If Channel 2 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 27)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. output channel 3: If Channel 3 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 13)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. output channel 4: If Channel 4 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 13)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

"first_frame_data.json"
-----------------------
is a `json file <https://en.wikipedia.org/wiki/JSON>`_ and contains
information that collected during
:ref:`animal detection<AnimalDetectionExplanationLabel>` (Source
code :class:`pre_experiment.FindAnimal`.)

#. bounding box col max: The X_max value of the bounding box of the
   animal detected in the first frame during animal detection.

#. bounding box col min: The X_min value of the bounding box of the
   animal detected in the first frame during animal detection.

#. bounding box row max: The Y_min value of the bounding box of the
   animal detected in the first frame during animal detection.

#. bounding box row min: The Y_max value of the bounding box of the
   animal detected in the first frame during animal detection.

#. centroid col: The X value of the centroid of the animal detected
   in the first frame during animal detection.

#. centroid row: The Y value of the centroid of the animal detected
   in the first frame during animal detection.

#. filled area: The filled area in pixels of the blob defined as
   the animal in the first frame during animal detection

Optional files
--------------

The rest of this subsection describes files that need to be explicitly
saved as described here. TODO LINK

"bounding_boxes.npy"
~~~~~~~~~~~~~~~~~~~~~~
is a `Numpy <https://numpy.org/>`_ file. It contains the coordinates
of the bounding box of the small image. The bounding box defines the
Y/X coordinates of the small image

This file comes in
`shape <https://docs.scipy.org/doc/numpy/reference/generated/numpy
.ndarray.shape.html>`__
[4, # of frames] with:

======  =========================
[0, :]  contains the Y_min values
------  -------------------------
[1, :]  contains the Y_max values
------  -------------------------
[2, :]  contains the X_min values
------  -------------------------
[3, :]  contains the X_max values
======  =========================

These values are necessary to describe where in the full image
frame the small image that has been saved during the experiment is
located. The bounding box is the rectangle that contains all image
information used during this frame. Below an illustration on how
the different values are used to construct the bounding box.

.. figure:: Figures/ouput/bbox_illustration.png
  :width: 100 %
  :alt: bbox_illustration.png

.. note::

  Why Y/X and not X/Y? In image processing the convention is to
  reference points in (Rows, Columns) which translates to Y/X. The
  underlying image processing libraries work with the (Rows,
  Columns) convention. See for example `here
  <https://scikit-image.org/docs/dev/api/skimage.measure.html#skimage
  .measure.regionprops>`__.
  PiVR therefore follows this convention.

"centroids.npy"
~~~~~~~~~~~~~~~~~~~~~~
is a `Numpy <https://numpy.org/>`_ file. It contains the coordinates
of the centroid of the blob identified during the experiment.
See :ref:`here<OutputPointExp>` to see the centroid compared to the
midpoint.

The file comes in `shape <https://docs.scipy
.org/doc/numpy/reference/generated/numpy.ndarray.shape.html>`__
[# of frames, 2] with:

======  ==============================
[:, 0]  contains the centroid Y values
------  ------------------------------
[:, 1]  contains the centroid X values
======  ==============================

These values are identical to what you will find in the
"DATE_TIME_data.csv" file

"midpoints.npy"
~~~~~~~~~~~~~~~~~~~~~~
is a `Numpy <https://numpy.org/>`_ file. It contains the coordinates
of the midpoint extracted from the skeleton during the experiment.
See :ref:`here<OutputPointExp>` to see the midpoint compared to the
centroid.

The file comes in `shape <https://docs.scipy
.org/doc/numpy/reference/generated/numpy.ndarray.shape.html>`__
[# of frames, 2] with:

======  ==============================
[:, 0]  contains the midpoint Y values
------  ------------------------------
[:, 1]  contains the midpoint X values
======  ==============================

These values are identical to what you will find in the
"DATE_TIME_data.csv" file

"heads.npy"
~~~~~~~~~~~~~~~~~~~~~~
is a `Numpy <https://numpy.org/>`_ file. It contains the coordinates
of the head position :ref:`assigned during tracking<HeadTailLabel>`.

The file comes in `shape <https://docs.scipy
.org/doc/numpy/reference/generated/numpy.ndarray.shape.html>`__
[# of frames, 2] with:

======  ==========================
[:, 0]  contains the head Y values
------  --------------------------
[:, 1]  contains the head X values
======  ==========================

These values are identical to what you will find in the
"DATE_TIME_data.csv" file

"tails.npy"
~~~~~~~~~~~~~~~~~~~~~~
is a `Numpy <https://numpy.org/>`_ file. It contains the coordinates
of the tail position :ref:`assigned during tracking<HeadTailLabel>`.

The file comes in `shape <https://docs.scipy
.org/doc/numpy/reference/generated/numpy.ndarray.shape.html>`__
[# of frames, 2] with:

======  ==========================
[:, 0]  contains the tail Y values
------  --------------------------
[:, 1]  contains the tail X values
======  ==========================

These values are identical to what you will find in the
"DATE_TIME_data.csv" file

"stimulation.npy"
~~~~~~~~~~~~~~~~~
is a `Numpy <https://numpy.org/>`_ file. It contains the stimulus
delivered to the animal during the experiment.

These values are identical to what you will find in the
"DATE_TIME_data.csv" file

Full Frame Recording
====================

After taking a lot of images with
:ref:`Full Frame Recording<FullFrameLabel>`, find a folder with the
"DATE_TIME_EXP.GROUP" as its name. An example would be
"2019.01.11_14-00-05_CantonS". This is an experiment conducted on the
11th of January 2019. "CantonS" is the value that was entered in the field
"Exp. Group".

"DATE_TIME_data.csv"
--------------------
contains the following data for each frame of the video:

   #. Frame (=image) number into the experiment
   #. Time in seconds since the experiment started
   #. Channel 1 stimulus delivered
   #. Channel 2 stimulus delivered
   #. Channel 3 stimulus delivered
   #. Channel 4 stimulus delivered

Image files
------------

Usually lots upon lots of them. Each image is saved separately directly
into this folder.

"experiment_settings.json"
--------------------------
is a `json file <https://en.wikipedia.org/wiki/JSON>`_ and contains a
lot of useful experimental information:

#. Experiment Date and Time: Exactly as advertised

#. Framerate: The framerate the video was recorded in

#. Exp. Group: The string that was entered by the user during the
   experiment

#. Model Organism: If selected, what animal has been indicated
   during the experiment.

#. Pixel per mm: If defined (see :ref:`here<PixMMLabel>`) a useful
   parameter for analysis.

#. Recording time: The time in seconds that PiVR was recording
   this video.

#. Resolution: The camera resolution in pixels that PiVR used while
   recording the video.

#. Virtual Reality arena name: As no virtual arena was presented,
   it will say 'None'

#. backlight 2 channel: If Backlight 2 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 18)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. backlight channel: If Backlight 1 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 18)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].
   This would normally be defined as [18, 40000].

#. output channel 1: If Channel 1 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 17)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. output channel 2: If Channel 2 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 27)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. output channel 3: If Channel 3 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 13)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. output channel 4: If Channel 4 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 13)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

Video
=====

After recording a :ref:`video<VideoLabel>`, you will
find a folder with the "DATE_TIME_EXP.GROUP" as its name. An example
would be "2019.01.11_14-00-05_CantonS". This is an experiment
conducted on the 11th of January 2019. "CantonS" is the value that
was entered in the field "Exp. Group".

"EXPGRP_VIDEO.h264"
-------------------
is the video file. This video file on its own is not perfectly useful
(at least in my hands) as h264 seems to be a bit of an exotic file
format that many video players can not handle without problems.

In order to directly convert this file, see the :ref:`Image Data
handling instructions <h264_to_avi_GUI>`. If you want to a
GUI-free version of these modules, check out the
"convert_h264_to_AVI.py" at
https://gitlab.com/davidtadres/pivr_bonus

.. note::
   I have tried to directly convert the image using ffmpeg. I believe
   there is a bug somewhere in the encoder of the camera as ffmpeg
   reads that the video is "inf" long. The scripts above take the
   video metadata from "experiment_settings.json" to properly convert
   the video.

The standard lens introduces a lot of radial aberrations at the
edges! To fix them see :ref:`undistort h264 video <undistort_h264_GUI>`.

.. figure:: Figures/ouput/undistortExample.png
  :width: 100 %
  :alt: undistortExample.png

"DATE_TIME_data.csv"
--------------------
contains the following data for each frame of the video:

   #. Frame (=image) number into the experiment
   #. Time in seconds since the experiment started
   #. Channel 1 stimulus delivered
   #. Channel 2 stimulus delivered
   #. Channel 3 stimulus delivered
   #. Channel 4 stimulus delivered

"experiment_settings.json"
--------------------------
is a `json file <https://en.wikipedia.org/wiki/JSON>`_ and contains a
lot of useful experimental information:

#. Experiment Date and Time: Exactly as advertised

#. Framerate: The framerate the video was recorded in

#. Exp. Group: The string that was entered by the user during the
   experiment

#. Model Organism: If selected, what animal has been indicated
   during the experiment.

#. Pixel per mm: If defined (see :ref:`here<PixMMLabel>`) a useful
   parameter for analysis.

#. Recording time: The time in seconds that PiVR was recording
   this video.

#. Resolution: The camera resolution in pixels that PiVR used while
   recording the video.

#. Virtual Reality arena name: As no virtual arena was presented,
   it will say 'None'

#. backlight 2 channel: If Backlight 2 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 18)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. backlight channel: If Backlight 1 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 18)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].
   This would normally be defined as [18, 40000].

#. output channel 1: If Channel 1 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 17)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. output channel 2: If Channel 2 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 27)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. output channel 3: If Channel 3 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 13)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].

#. output channel 4: If Channel 4 has been defined (as
   described :ref:`here<DefineGPIOsLabel>`) the chosen GPIO (e.g. 13)
   and the maximal PWM frequency (e.g. 40000) is saved as a [list].


Get started
===========

Different experiments necessitate different analysis. In the original
PiVR publication
`PiVR publication <https://doi.org/10.1371/journal.pbio.3000712>`__
a number of different experiments were run
and the analysis
`analysis <https://gitlab.com/LouisLab/pivr_publication>`__ and data
of those has been
made `public <https://doi.org/10.25349/D9ZK50>`__.
These scripts are all annotated and you
should be able to run them on your computer with the original data
to understand what is happening in them. Then you can adapt them/use
them with your data.

In addition, the PiVR software has a couple of built in analysis
tools when run on a PC (i.e. not on a Raspberry Pi):

.. _OutputPointExp:

Visualization of different points on the animal
===============================================

What exactly do the terms "Centroid" and "Midpoint" mean? I will try
to illustrate the difference so that you may choose the appropriate
parameter for your experiment:

#. To identify the animal the tracking algorithm identifies a "blob"
   that has significantly different pixel intensity values compared
   to the background.

#. The centroid is the center of mass (in 2D) of these pixels.

#. The midpoint is the center of the skeletonized blob.

   .. figure:: Figures/ouput/CM_illustration.png
      :width: 100 %
      :alt: CM_illustration.png.. Intended Audience: User who wants to better understand how the
   software works without looking at the code.

PiVR has been developed by David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

.. _CodeExplanationLabel:

Code explanation
****************

Read this guide if you:

#. Want to gain a high-level understanding of how the detection
   and tracking algorithm work.

#. Do not (yet) want to look into the source code

.. _AnimalDetectionExplanationLabel:

Animal detection
================

The :ref:`tracking algorithm <AnimalTrackingLabel>`  depends on the
approximate location of the animal and a background image in order to
identify the animal. The *Animal Detection* modes
(:ref:`guide <AnimalDetectionGuideLabel>`, and
:ref:`selection <AnimalDetectionLabel>`) allow the user to choose
between one of three Modes depending on the experimental
setup and question asked. (:ref:`Source Code <AnimalDetectionClass>`)

.. _CodeExplanationMode1Label:

Standard - Mode 1
-----------------

**User perspective:**

#. Place animal in the final dish on the arena, taking the
   :ref:`optimal image parameters into account <SettingUpArenaCamera>`
#. Press *Start*

**Behind the scenes:**

#. When the user presses *Start* the camera starts sending pictures
   to the tracking software.
#. For the first frame, the image is just filtered using a gaussian
   filter with sigma depending on the properties saved in
   'list_of_available_organisms.json'. In short, the sigma will be
   half the minimal expected cross-section of the animal. Then it
   takes the next frame.

    .. figure:: /Figures/code_explanation/FirstImage.png
        :width: 100 %
        :alt: FirstImage.png

#. Starting with the second frame, the mean of the current and all
   previously taken images is taken.

    .. figure:: /Figures/code_explanation/SecondMeanImage.png
        :width: 100 %
        :alt: SecondMeanImage.png

#. The current frame is then subtracted from the mean of all previous
   frames

    .. figure:: /Figures/code_explanation/SecondSubtractedImage.png
        :width: 100 %
        :alt: SecondSubtractedImage.png

#. The subtracted image has a trimodal distribution of pixels - The
   main peak is the background (grey arrows). The intensity values
   for pixels where the animal was before but isn't anymore has
   positive values (cyan arrow). The intensity values for the pixels
   where the animal only recently moved in has negative values
   (magenta arrow)

    .. figure:: /Figures/code_explanation/HistogramMode1.png
        :width: 100 %
        :alt: HistogramMode1.png

#. The threshold used to define the animal is subtracting 2 *
   the standard deviation of the smoothed image from the mean value
   of pixel intensities

    .. figure:: /Figures/code_explanation/HistogramMode1Tresh.png
        :width: 100 %
        :alt: HistogramMode1Tresh.png

#. The current image gets a threshold using the threshold defined above

    .. figure:: /Figures/code_explanation/Mode1Binary.png
        :width: 100 %
        :alt: Mode1Binary.png

#. Using this binary image, the function
   `label <http://scikit-image.org/docs/dev/api/skimage.measure
   .html#skimage.measure.label>`_ and then the function
   `regionprops <http://scikit-image.org/docs/dev/api/skimage.measure
   .html#skimage.measure.regionprops>`_ of the scikit-image library
   is applied.
#. Using the *filled_area* parameter of the regionprops function, the
   blobs are checked for minimal size defined in the
   'list_of_available_organisms.json' file for the animal in
   question. As soon as an image is found where a single blob is
   larger than the minimal size, the animal counts as being identified.
#. The first picture is then filtered again with a gaussian filter
   with a sigma of 1. This is defined as the background image. This
   image of course contains the animal as it was at the original
   position. For many experiments it is unlikely that the animal can
   be at exactly the same position as in the first frame. If your
   experiment makes it likely that the animal is in the exact same
   position more than once, you might want to look at Animal
   Detection :ref:`Mode#2 <CodeExplanationMode2Label>` and
   :ref:`Mode#3 <CodeExplanationMode3Label>`.

    .. figure:: /Figures/code_explanation/Mode1Background.png
        :width: 100 %
        :alt: Mode1Background.png

#. All of this is done in the function
   :meth:`pre_experiment.FindAnimal.find_roi_mode_one_and_three`
#. Finally, another frame is taken from the camera. *Importantly, the
   threshold is now calculated locally!* This allows for a better
   separation of the animal and the background by defining a local
   treshold. The image gets filtered using the local threshold and
   then subtracted from the first (filtered) image. After applying
   the `label <http://scikit-image.org/docs/dev/api/skimage.measure
   .html#skimage.measure.label>`_ and then the  `regionprops
   <http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage
   .measure.regionprops>`_ functions of the scikit-image library the
   largest blob is defined as the animal. The centroid, the filled
   area and the bounding box are saved for the tracking algorithm.

#. Example data collected for this particular example with Animal
   Detection Mode#1:

    .. figure:: /Figures/code_explanation/Mode1Output.png
        :width: 100 %
        :alt: Mode1Output.png

    **Result of Animal Detection Mode 1:** *Top Left:* Background
    image saved for the rest of the experiment. Using Mode 1 the
    animal will always be present in the background image used during
    the actual experiment. If this is a problem for your experiment,
    please check :ref:`Mode#2<CodeExplanationMode2Label>` and
    :ref:`Mode#3<CodeExplanationMode3Label>`.
    *Top Right:* Animal identified during Animal Detection
    :ref:`Mode#1<CodeExplanationMode1Label>`. The red box indicates
    the bounding box of the animal.
    *Bottom Left:* Close-up of the detected animal. The bounding box
    is defined with 4 coordinates: The smallest row (Row Min, or
    Y-min), the smallest column (Col Min, or X-min), the largest row
    (Row Max, or Y-max) and the largest column (Col Max, or X-Max).
    Also see the `regionprops <http://scikit-image
    .org/docs/dev/api/skimage.measure.html#skimage.measure.regionprops>`_
    function documentation (*bbox*). *Bottom Right:* The binary image
    used by the algorithm to define the animal. The  centroid is
    indicated as well as the filled area. Both parameters are defined
    using the `regionprops <http://scikit-image
    .org/docs/dev/api/skimage.measure.html#skimage.measure
    .regionprops>`_ function.

#. This is done using:
   :func:`pre_experiment.FindAnimal.define_animal_mode_one`
#. Next, the :ref:`actual tracking algorithm
   <AnimalTrackingLabel>` will be called and
   run until the experiment is over.

.. _CodeExplanationMode2Label:

Pre-Define Background - Mode 2
------------------------------

**User perspective:**

#. Prepare PiVR for the experiment by putting everything *except* the
   animal **exactly** at the position where it will be during the
   experiment. Take the :ref:`optimal image parameters into account
   <SettingUpArenaCamera>`
#. Press *Start*
#. User will be asked to take a picture by clicking 'OK'.
#. Then user will be asked to put the animal **without changing
   anything else in the Field of View of the camera**
#. Press *Start*

**Behind the scenes:**

#. The image taken without the animal is being used as the background
   image.

    .. figure:: /Figures/code_explanation/Mode2Background.png
        :width: 100 %
        :alt: Mode2Background.png

#. After the user puts the animal and hits 'Ok' the same algorithm as
   in Mode 1 and 3 is searching for the animal: First a new picture
   is taken and filtered using a gaussian filter with Sigma = 1:

    .. figure:: /Figures/code_explanation/Mode2FirstImage.png
        :width: 100 %
        :alt: Mode2FirstImage.png

#. This image is subtracted from the background image:

    .. figure:: /Figures/code_explanation/Mode2SubtractedImage.png
        :width: 100 %
        :alt: Mode2SubtractedImage.png

#. The pixel intensity values are a bimodal distribution.

    .. figure:: /Figures/code_explanation/Mode2Histogram.png
        :width: 100 %
        :alt: Mode2Histogram.png

#. The threshold is defined as all values larger 2 times the Standard
   deviation of the the image + the mean pixel intensity values of
   the image (usually 0)

    .. figure:: /Figures/code_explanation/Mode2HistogramThresh.png
        :width: 100 %
        :alt: Mode2HistogramThresh.png

    **Threshold to identify animal:** On the left the background
    pixels with a approximate value of zero are seen. Note the
    logarithmic scale. On the right the yellow rectangle indicates
    all the pixel intensity values that will count as the identified
    animal

#. The threshold is used to binarize the image.

    .. figure:: /Figures/code_explanation/Mode2Binary.png
        :width: 100 %
        :alt: Mode2Binary.png

#. Using this binary image, the function `label <http://scikit-image
   .org/docs/dev/api/skimage.measure.html#skimage.measure.label>`_ and
   then the function `regionprops <http://scikit-image
   .org/docs/dev/api/skimage.measure.html#skimage.measure
   .regionprops>`_ of the scikit-image library is applied.
#. Using the *filled_area* parameter of the regionprops function, the
   blobs are checked for minimal size defined in the
   'list_of_available_organisms.json' file for the animal in question.
   As soon as an image is found where a single blob is larger than
   the minimal size, the animal counts as being identified.
#. This is done in :meth:`pre_experiment.FindAnimal.find_roi_mode_two`
#. Unlike :ref:`Mode1 <CodeExplanationMode1Label>` it is *not*
   necessary to take local thresholding as the animal should be
   clearly visible compared to the background. The global threshold
   is used to create a binary image. The largest blob is defined as
   the animal.
#. Example data collected for this particular example with Animal
   Detection Mode#2:

    .. figure:: /Figures/code_explanation/Mode2Output.png
        :width: 100 %
        :alt: Mode2Output.png

    **Result of Animal Detection Mode 1:** *Top Left:* Background
    image saved for the rest of the experiment. *Top Right:* Animal
    identified during :ref:`Animal Detection Mode 2
    <CodeExplanationMode2Label>`. The red box  indicates the bounding
    box of the animal. *Bottom Left:* Close-up of  the detected
    animal. The bounding box is defined with 4 coordinates: The
    smallest row (Row Min, or Y-min), the smallest column (Col Min,
    or X-min), the largest row (Row Max, or Y-max) and the largest
    column (Col Max, or X-Max). Also see the `regionprops
    <http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.regionprops>`_
    function documentation (*bbox*). *Bottom Right:* The binary image
    used by the algorithm to define the animal. The centroid is
    indicated as well as the filled area. Both parameters are defined
    using the `regionprops <http://scikit-image
    .org/docs/dev/api/skimage.measure.html#skimage.measure
    .regionprops>`_ function.

#. This is done in
   :func:`pre_experiment.FindAnimal.define_animal_mode_two`
#. Next, the :ref:`actual tracking algorithm <AnimalTrackingLabel>` will
   be called and run until the experiment is over.

.. _CodeExplanationMode3Label:

Reconstruct Background by Stitching - Mode 3
--------------------------------------------

**User perspective:**

#. Place animal in the final dish on the arena, taking the
   :ref:`optimal image parameters into account <SettingUpArenaCamera>`
#. Press *Start*

Behind the scenes:

#. As soon as the user hits "Start", the camera will start streaming
   images from the camera. Here a fruit fly larva can be seen in the
   center of the image

    .. figure:: /Figures/code_explanation/Mode3FirstImage.png
        :width: 100 %
        :alt: Mode3FirstImage.png

#. A second image then taken.

   .. figure:: /Figures/code_explanation/Mode3SecondImage.png
        :width: 100 %
        :alt: Mode3SecondImage.png

#. The new image is now subtracted from previous image(s).

   .. figure:: /Figures/code_explanation/Mode3SecondSubtractedImage.png
        :width: 100 %
        :alt: Mode3SecondSubtractedImage.png

#. The histogram of this image exhibits a trimodal distribution - a
   very large peak at around 0 indicating all the background pixels,
   positive pixel intensities indicating coordinates that the
   animal occupied in the past and has left and negative pixel values
   indicating coordinates that the animal has not occupied in the
   first frame but does now.

   .. figure:: /Figures/code_explanation/Mode3HistogramSubtractedImage.png
        :width: 100 %
        :alt: Mode3HistogramSubtractedImage.png

#. The threshold is defined as 4 minus the mean pixel intensity of
   the subtracted image.

   .. figure:: /Figures/code_explanation/Mode3HistogramMode1Tresh.png
        :width: 100 %
        :alt: Mode3HistogramMode1Tresh.png

#. The threshold is then used to binarize the first image:

   .. figure:: /Figures/code_explanation/Mode3Binary.png
        :width: 100 %
        :alt: Mode3Binary.png

#. Using this binary image, the function
   `label <http://scikit-image.org/docs/dev/api/skimage.measure
   .html#skimage.measure.label>`_ and then the function
   `regionprops <http://scikit-image.org/docs/dev/api/skimage.measure
   .html#skimage.measure.regionprops>`_ of the scikit-image library
   is applied.
#. Using the *filled_area* parameter of the regionprops function, the
   blobs are checked for minimal size defined in the
   'list_of_available_organisms.json' file for the animal in
   question. As soon as an image is found where a single blob is the
   region of interest in calculated depending on the maximum speed,
   maximum size and pixel/mm.

   .. figure:: /Figures/code_explanation/Mode3ROIthresh.png
        :width: 100 %
        :alt: Mode3ROIthresh.png

#. So far this has been identical to
   :ref:`Mode 1<CodeExplanationMode1Label>` - both Modes use the
   identical function up to this point:
   :func:`pre_experiment.FindAnimal.find_roi_mode_one_and_three`

#. Now comes the trickiest part of Mode 3: The animal must be
   identified as complete as possible. For this the region of
   interest is used:

   .. figure:: /Figures/code_explanation/Mode3ROI_Raw.png
        :width: 100 %
        :alt: Mode3ROI_Raw.png

#. The histogram of this region of interest looks very different
   compared to the histogram of the whole frame. Importantly, the
   larva now clearly stands out from the background as can be seen in
   the histogram.

   .. figure:: /Figures/code_explanation/Mode3ROI_histogram.png
        :width: 100 %
        :alt: Mode3ROI_histogram.png

#. The threshold is now adjusted until only **a single object with
   animal-like properties** is left in the thresholded region of
   interest. ^

   ^animal like properties are defined in the
   "list_of_available_organisms.json" file

   .. figure:: /Figures/code_explanation/Mode3MovingTheSTD.png
        :width: 100 %
        :alt: Mode3MovingTheSTD.png

#. **It is important to understand the limitations of this approach!**
   This animal identification will only work if the animal is
   clearly separated from the background. It will not work, for
   example if the animal is close to the edge of the petri dish, if
   there are other structures in the arena etc...!

   All this is done using the function:
   :func:`pre_experiment.FindAnimal.define_animal_mode_three`

#. After identifying the animal, the algorithm waits until the animal
   has left the inital position.

   To do this it will continue capturing images, binarizing them
   and finally subtracting them from the identfied animal. Only when
   the original animal is 99.9% reconstructed has the animal left the
   original position.

   This is done using the function
   :func:`pre_experiment.FindAnimal.animal_left_mode_three`

   .. figure:: /Figures/code_explanation/Mode3Waiting.png
        :width: 100 %
        :alt: Mode3Waiting.png

#. When the animal has left, the region that was occupied by the
   animal in the first frame is replaced by the pixels at the same
   coordinates of the image where the animal has left the original
   position. This should lead to a "clean" background image, meaning
   the animal shouldn't be present (as it would be in Mode 1).

   This is done using the function:
   :func:`pre_experiment.FindAnimal.background_reconstruction_mode_three`

   .. figure:: /Figures/code_explanation/Mode3reconstruction.png
        :width: 100 %
        :alt: Mode3reconstruction.png

#. Before the tracking can start, the location of the animal
   after background reconstruction must be saved. This is done using
   the function:
   :func:`pre_experiment.FindAnimal.animal_after_box_mode_three`


.. _AnimalTrackingLabel:

Animal Tracking
===============

After the detection of the animal with any of the described animal
detection modes (:ref:`Mode 1<CodeExplanationMode1Label>`,
:ref:`Mode 2<CodeExplanationMode2Label>` or
:ref:`Mode 3<CodeExplanationMode3Label>`) the tracking algorithm
(:func:`fast_tracking.FastTrackingVidAlg.animal_tracking`)
starts.

.. note::

     The camera will start recording at the pre-set frame rate. If the
     frame rate exceed the time it takes for PiVR to process the frame,
     the next frame will be dropped. For example, if you are recording
     at 50 frames per second, each frame has to be processed in 20ms
     (1/50=0.02 seconds). If the frame processing takes 21ms, the
     next frame is dropped and PiVR will be idle for the next 19ms
     until the next frame arrives.

#. Once the tracking algorithm starts, the camera and GPU start sending
   images at the defined frame rate to the CPU, i.e. at 30fps the CPU
   will receive one new image ever 33ms. Source code here:
   :func:`fast_tracking.FastTrackingControl.run_experiment`

#. The images are then prepared for tracking. Source code here:
   :func:`fast_tracking.FastTrackingVidAlg.write`

#. Next the animal_tracking function which does all the heavy lifting
   described below is called:
   :func:`fast_tracking.FastTrackingVidAlg.animal_tracking`.

#. The Search Box is defined depending on the previous animal
   position, the selected organism and that organisms specific
   parameters. For details check source code here:
   :func:`start_GUI.TrackingFrame.start_experiment_function`

   .. figure:: /Figures/code_explanation/FigS4a_FlowchartIndicateSearchBox.png
        :width: 100 %
        :alt: FigS4a_FlowchartIndicateSearchBox.png

#. The content of the Search Box of the current frame is then
   subtracted from the background frame (generated during
   :ref:`animal detection<AnimalDetectionExplanationLabel>`).

   .. figure:: /Figures/code_explanation/FigS4Subtract.png
        :width: 100 %
        :alt: FigS4Subtract.png

#. Upon inspection of the histogram of the subtracted image it becomes
   clear that the animal has clearly different pixel intensity values
   compared to the background.

   .. figure:: /Figures/code_explanation/FigS4_Histogram.png
        :width: 100 %
        :alt: FigS4_Histogram.png

#. The threshold is determined by calculating the mean pixel intensity
   of the subtracted image and subtracting 3 times the standard
   deviation.

   .. figure:: /Figures/code_explanation/FigS4_Histogram2.png
        :width: 100 %
        :alt: FigS4_Histogram2.png

#. This threshold is then used to binarize the current small image:

   .. figure:: /Figures/code_explanation/FigS4_Binary.png
        :width: 100 %
        :alt: FigS4_Binary.png

#. Using that binary image, the function
   `label <http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.label>`_
   and then the function
   `regionprops <http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.regionprops>`_
   of the scikit-image library is applied. This way all 'blobs' are
   identified. Using the animal parameters defined before, the blob
   that looks the most like the sought after animal is assigned to
   being the animal.

   After detection of the animal, the image of the animal is saved
   (blue rectangle) and the Search Box for the next frame prepared
   (red rectangle)

   .. figure:: /Figures/code_explanation/FigS4_Animal_detected.png
        :width: 100 %
        :alt: FigS4_Animal_detected.png

.. _HeadTailLabel:

Head Tail Classification
========================

The Head Tail classification is based upon a Hungarian algorithm.

#. First, the binary image is skeletonized (either with `thin
   function
   <https://scikit-image.org/docs/dev/api/skimage.morphology.html#skimage.morphology.thin>`_
   or the `skeletonize function
   <https://scikit-image.org/docs/dev/api/skimage.morphology.html#skimage.morphology.skeletonize>`_)

   .. figure:: /Figures/code_explanation/FigS5_Skeletonize.png
        :width: 100 %
        :alt: FigS5_Skeletonize.png

#. Using the rule that the endpoints of that skeleton must only have
   one neighbour, the endpoints are defined.

   .. figure:: /Figures/code_explanation/FigS5_Endpoints.png
        :width: 100 %
        :alt: FigS5_Endpoints.png

#. To define the head and tail the following conditions must
   be met:

   A. The aspect ratio of the long axis over the short axis must be
      at least 1.25

   B. The skeleton must have exactly 2 endpoints

   C. The length of the skeleton must be larger than half the
      mean length of the previous 3 frames.

#. Next, the distance of each endpoint to a reference point is
   calculated:

   A. In case the tail has not yet been assigned (happens in the first
      frame) use the centroid of the previous frame as the reference
      point.

   B. In case of not having been able to assign a tail in the
      previous frame, e.g. due to the violation of any of the rules
      shown above, also use the centroid of the previous frame as the
      reference point.

   C. Otherwise (in most cases) the endpoint that has been
      assigned the tail in the previous frame is used as the
      reference point.

   .. figure:: /Figures/code_explanation/FigS5_DistanceEndpoints.png
        :width: 100 %
        :alt: FigS5_DistanceEndpoints.png

#. Whichever endpoint has the shorter distance the previous reference
   point is assigned the label 'Tail'.

   .. figure:: /Figures/code_explanation/FigS5_TailClassified.png
        :width: 100 %
        :alt: FigS5_TailClassified.png

.. _Dynamic Virtual Realities explained:

Dynamic Virtual Realities
=========================

PiVR is capable of presenting virtual realities that change over
time, or 'dynamic virtual realities'. The difference to the
'standard' virtual reality is that the dynamic virtual reality will
change its stimulus profile over time even if the animal is not
moving. See the following move for an application example where a
fruit fly larva is tested in a virtual odor plume.

.. raw:: html

    <iframe width="75%"
    src="https://www.youtube.com/embed/RGAta8VqIlw" frameborder="0"
    allow="accelerometer; autoplay; encrypted-media; gyroscope;
    picture-in-picture" allowfullscreen></iframe>

|

To create such a dynamic virtual reality from scratch,
:ref:`see here<create_dynamic_VR_arena>`.

How does PiVR present such a dynamic virtual reality? The stimulation
file is a 3D numpy array with the 3rd dimension coding for time. I
will call the virtual arena at a given time *slice*.

.. figure:: /Figures/code_explanation/dVR_Numpy3D_array.png
     :width: 100 %
     :alt: dVR_Numpy3D_array.png

For each timepoint during the experiment, one of the x/y slices is
presented to the animal as a virtual arena. How does PiVR know when to
update the arena? You could, for example update the arena 15
times per second (Hz), or maybe only once every minute (1Hz).

The relevant parameter is saved directly in the filename of the arena
as Hz[xx].npy

.. figure:: Figures/tools/10_1_dynamic_arena_name.png
   :alt: 10_1_dynamic_arena_name.png

For example, you are running your experiment at 30 frames per second.
You want to update the arena 15 times per second. The following
pseudo-code will explain how PiVR updates the arena:

.. code-block:: python

    # User defined update parameter
    recording_framerate = 30 # Hz
    dynamic_virtual_reality_update_rate = 15 # Hz
    update_VR_every_x_frame = recording_framerate/dynamic_virtual_reality_update_rate
    recording_time = 5 # minutes
    recording_length = recording_framerate * recording_time * 60
    current_VR_time_index = 0

    # The for loop indicates a running experiment where PiVR processes
    # one frame per iteration
    for current_frame in range(recording_length):
        detect_animal_function()
        detect_head_function()

        # Here we update VR according to the update_VR_every_x_frame
        # Note, the modulo operator % returns the remainder of a
        # division. For example 0%2=0, 1%2=1, 2%2=0, 3%2=1 etc.
        if current_frame% update_VR_every_x_frame == 0:
            current_VR_time_index += 1

        # Now, we present the current slice of the virtual reality
        present_VR(stimulation_file[:, :, current_VR_time_index])

        # As the size of the stimulus file in time can be smaller than
        # the length of recording_time the index needs to be reset to
        # zero.
        if current_VR_time_index > size_of_stimulus_file:
            current_VR_time_index = 0

As an example for how this would play out, consider an experiment
with the following parameters:

    #. Experiment time = 5 minutes (300 seconds)
    #. The recording framerate = 30Hz, hence there will be 9000 frames
    #. dynamic VR update rate = 15Hz
    #. dynamic VR file has 450 time slices.

You can immediately see that there are less *slices* in the provided
stimulus file compared to the other experimental parameters: We want
to record for 300 seconds while updating 15 times per second with
only 450 slices. We will run out of slices after only 450/15=30 seconds!

PiVR solves this by just starting over at the beginning when it
encounters the end of the stimulus file. The figure shown below is
visualizing this:

.. figure:: Figures/code_explanation/dVR_counter.png
   :alt: dVR_counter.png
PiVR has been developed by David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

.. _software_installation:

PiVR software installation
==========================

.. _software_installation_RPi:

Install PiVR on the Raspberry Pi
--------------------------------

The easiest way to install PiVR on your Raspberry Pi is to just
follow the :ref:`instructions during hardware construction<PIVR
Hardware Manual>`.

If you just want to download the installation script to your
Raspberry Pi, :download:`press
here <../install_PiVR.sh>`

Open the terminal. Then change directory to the 'Downloads' folder
(or wherever you downloaded the file) and type::

      bash install_PiVR.sh

.. _software_installation_PC:

Install PiVR on a PC
---------------------

#. Install `miniconda <https://docs.conda.io/en/latest/miniconda.html>`__
   on your computer.

#. Install `git <https://git-scm.com/downloads>`__ on your computer

.. note::
   If you have Windows, you may try :ref:`this
   guide<software_install_Win10>` which will install the software
   more or less automatically.

.. note::
   If you have Ubuntu, you may try
   :ref:`this guide<software_install_Ubuntu>` which will install the
   software more or less automatically.

#. Now, create an empty conda environment:

   .. code-block:: python

      conda create --name PiVR_environment

#. Activate the environment you just created by typing:

   Linux/Mac:

   .. code-block:: python

      source activate PiVR_environment

   Windows:

   .. code-block:: python

      activate PiVR_environment

#. Install the a number of packages which are necessary to run the PiVR
   software by copying each line of code into the Terminal

   .. code-block:: python

      conda install -y python=3.9

      conda install -y matplotlib

      conda install -y pandas

      conda install -y scipy

      conda install -y natsort

      conda install -y scikit-image

      conda install -c conda-forge opencv

      Windows/Linux:
      conda install -c conda-forge imageio-ffmpeg

      MacOS:
      conda install -c conda-forge ffmpeg

#. You have now prepared the
   `virtual environment <https://towardsdatascience.com/getting-started-with-python-environments-using-conda-32e9f2779307>`__ PiVR will be running in.

#. Using the anaconda terminal, change the working directory to a
   folder where you want to store the actual PiVR software.

    .. code-block:: python

       cd C:\Users\UserA\Documents>

    .. note::

       You might want to write down the exact path so that you will
       find it again in the future!

#. Download the software by typing:

    .. code-block:: python

       git clone https://gitlab.com/louislab/PiVR

#. Now navigate into the folder you have just downloaded by typing:

   .. code-block:: python

      cd PiVR

#. To start the PiVR software type:

   .. code-block:: python

      python start_GUI.py

.. _software_install_Win10:

Install PiVR on a Windows 10 PC
-------------------------------

.. important::

   If you are having trouble with this installation procedure, do the
   :ref:`manual install<software_installation_PC>`.

.. warning::

   Only Win10, 64bit tested!

#. Open the Anaconda prompt

#. Navigate into a folder where you want to store the PiVR software,
   for example:

   .. code-block:: python

      cd C:\Users\UserA\Documents>

#. Download the software by typing:

   .. code-block:: python

      git clone https://gitlab.com/louislab/PiVR

#. Navigate into the installation folder by typing:

   .. code-block:: python

      cd PiVR\Installation_update

#. Create the Windows 10 virtual environment for the PiVR software to
   run using the provided package list by typing:

   .. code-block:: python

      conda create --name PiVR_environment --file PiVR_Win64.txt

#. Once done, activate the virtual environment by typing:

   .. code-block:: python

      activate PiVR_environment

   You know you successfully activated the virtual enviroment if it
   says '(PiVR)' at the beginnig of the line in the terminal.

#. Start the software by going into the folder where the file
   "start_GUI.py" can be found, which is the parent folder of the
   installation folder you should be in now. So just type:

   .. code-block:: python

      cd ..

#. And to finally start PiVR, type:

   .. code-block:: python

      python start_GUI.py

.. _software_install_Ubuntu:

Install PiVR on a Linux PC
---------------------------

.. important::

   If you are having trouble with this installation procedure, do the
   :ref:`manual install<software_installation_PC>`.

.. warning::

   Only Ubuntu, 64bit tested)

#. Open the Terminal

#. Navigate into a folder where you want to store the PiVR software,
   for example:

   .. code-block:: python

      cd /home/UserA

#. Clone the repository by typing:

   .. code-block:: python

      git clone https://gitlab.com/louislab/PiVR

#. Navigate to the "Installation_update" folder of the repository you
   just cloned:

   .. code-block:: python

      cd /home/UserA/PiVR/PiVR/Installation_update

#. Create the Linux virtual environment for the PiVR software to
   run using the provided package list by typing:

   .. code-block:: python

      conda create --name PiVR_environment --file PiVR_Linux64.txt

#. Once done, activate the virtual environment by typing:

   .. code-block:: python

      source activate PiVR_environment


   You know you successfully activated the virtual enviroment if it
   says '(PiVR)' at the beginnig of the line in the terminal.

#. Start the software by going into the folder where the file
   "start_GUI.py" can be found, which is the parent folder of the
   installation folder you should be in now. So just type:

   .. code-block:: python

      cd ..

#. Start the program by typing:

   .. code-block:: python

      python start_GUI.py


.. _software_start_PC:

Start PiVR on a PC
-------------------

.. note::

   To run PiVR, you of course need to first
   :ref:`install<software_installation_PC>` the software.

#. Open the Anaconda terminal (Windows) or Terminal (MacOS/Linux)

#. Activate the virtual environment you have created during the
   installation. If you followed these instructions type:

   Windows:

   .. code-block:: python

      activate PiVR_environment

   Linux/MacOS:

   .. code-block:: python

      source activate PiVR_environment

#. Change directory to the folder where you downloaded the PiVR
   software into. In the example here we used:

   .. code-block:: python

      cd C:\Users\UserA\Documents\PiVR\PiVR

#. Start PiVR software by typing:

   .. code-block:: python

      python start_GUI.py
PiVR has been developed by David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

.. _advanced_topics:

Advanced topics
***************

.. _debug_tracking:

Simulate Real-Time Tracking
===========================

Imagine setting up your experiment: preparing the animals, booking
the setup/room for a whole afternoon...and then the tracker does not
track the animal half of the time!

It is quite frustrating sitting in a dark room trying to troubleshoot
these kind of problems.

To alleviate situations like these, there is the option to *simulate*
real time tracking after installing :ref:`the PiVR software on a
PC <software_installation_PC>` (=not on the Raspberry Pi).

#. At the PiVR setup, double check that:

   #. you have :ref:`set the resolution <AdjustImageLabel>` to
      resolution you want to use,
   #. that the :ref:`pixel/mm <PixMMLabel>` is set correctly,
   #. that the :ref:`framerate <AdjustImageLabel>` is identical to
      the framerate you are trying to achieve with real-time tracking.
   #. that you have :ref:`selected the correct animal <SelectOrganismLabel>`

#. Then, record a :ref:`video<VideoLabel>` with these settings. Then
   record some more.

#. Transfer the video data to your PC where you have installed the
   PiVR software and select the Debug->Simulate Online Tracking

   .. figure:: Figures/advanced_options/11_simulate_online.png
     :alt: 11_simulate_online.png
     :width: 100 %

#. Make sure the :ref:`Animal Detection Method
   <AnimalDetectionLabel>` is the same as the one you want to use
   during Real-Time tracking.

   .. note::
      This has not been tested with Mode 2

#. Select a single folder. You will now see the metadata created
   while the video was taken. Carefully inspect it to see if the
   settings are as you expect them.

   .. figure:: Figures/advanced_options/12_read_metadata.png
     :alt: 12_read_metadata.png
     :width: 100 %

#. Press the 'Track Single Animal' button - you will get a popup as
   soon as the animal detection algorithm detects a moving object.

   .. figure:: Figures/advanced_options/13_simulate_tracking_detection.png
     :alt: 13_simulate_tracking_detection.png
     :width: 100 %

#. After pressing ok, you will see what the animal detection
   algorithm has defined as the animal.

   It is obvious something has gone wrong here as the image on the
   right (the binary image) has a lot of spots where the image is
   white (=areas which are considered to be the animal)

   .. figure:: Figures/advanced_options/14_simulate_tracking_detection_broken.png
     :alt: 14_simulate_tracking_detection_broken.png
     :width: 100 %

#. After pressing Ok, the tracking algorithm starts - as the animal
   has not been properly identified in the first frame, the tracking
   algorithm is unable to identify the animal during tracking as well:

   .. figure:: Figures/advanced_options/15_simulated_tracking_broken.png
     :width: 100 %
     :alt: 15_simulated_tracking_broken.png

#. After going through the simulated tracking, the potential source
   of the problem has been identified: The animal can not be detected
   correctly. There are many reasons why this could be:

   #. The edge of the dish seems to have moved a bit during the first
      couple of frames (red rectangle). If you are able to stabilize
      the setup to ensure no movement while doing experiments, this
      problem should be solved.

      .. figure:: Figures/advanced_options/16_detection_problems1.png
         :width: 100 %
         :alt: 16_detection_problems1.png

   #. The fact that several spots in the middle of the dish are
      wrongly binarized as the potential animal, indicates that the
      detection algorithm has trouble setting the threshold correctly.
      This problem arises because the threshold is calculated as 2
      standard deviations from the mean of the pixel intensities in
      the subtracted image :func:`pre_experiment.FindAnimal.define_animal_mode_one`.
      While the animal is the darkest spot in the image, the whole
      petri dish is darker than the background which might lead to
      this problem.

      .. figure:: Figures/advanced_options/17_detection_problem2.png
         :width: 100 %
         :alt: 17_detection_problem2.png

      There are two general ways to solve this problem:

         #. Optimize the :ref:`imaging conditions<AdjustImageLabel>`
            so that the animal has a higher contrast to the
            background, which should be as homogenous as possible.
            See :ref:`here<SettingUpArenaCamera>` for an example.

         #. Optimize the animal parameters. You can follow :ref:`this
            guide <define_new_animal_label>` to set stringent animal
            parameters for tracking.

#. There are many ways how tracking can fail. Only a single
   example is described above. I hope the walkthrough will enable
   you to generally get an idea where during tracking the
   algorithm fails.

.. _define_new_animal_label:

Tracking of a new animal
========================

PiVR has been used to track a variety of animals: Walking adult fruit
flies, fruit fly larvae, spiders, fireflies, pillbugs and zebrafish.

For the tracking algorithm to function optimally, it takes several
"animal parameters" into account:

#. The amount of pixels the animal will fill in the image.

#. The "roundness" of the animal.

#. The proportions of the animal.

#. The length of the animal.

#. The speed of the animal.

For each animal you can choose in
:ref:`Options->Select Organism<SelectOrganismLabel>` these parameters
were defined. You can find them in the file
"list_of_available_organisms.json" in source code.

If you want to track an animal that is **not** on the list you can
always try to use the "Not in list" option. However, the tracking
algorithm might not work optimally.

There is a straightforward pipeline to collect the necessary animal
parameters to optimize real-time tracking:

#. Place your (**single!!**) animal in the arena you want to use for
   your experiment.

#. As always, do not forget to :ref:`define the pixel/mm <PixMMLabel>`.

#. Select "Not in List" under :ref:`Options->Select
   Organism<SelectOrganismLabel>`.

#. Record a :ref:`video<VideoLabel>`. If you use a fast animal, make
   sure to select a sufficiently high framerate. **You must use
   640x480 resolution**. As always it is imperative that the camera
   and the arena are stable during recording, i.e. nothing in the
   image should move except the animal!

#. Record for a couple of minutes, i.e. 5 minutes.

#. Make sure you have videos with animals moving as fast as they
   might in your actual experiment.

#. It is also necessary that the animals move for a large fraction of
   the video!

#. Take the videos to your PC on which you have :ref:`installed the
   PiVR software <software_installation_PC>`.

#. To observe what the algorithm is doing, turn the :ref:`Debug mode
   on<DebugModeLabel>`. This is recommend as you will see immediately
   if and where something goes wrong. This can help you to solve
   tracking problems.

#. Analyze each video using the :ref:`Tools->Analysis: Single Animal
   Tracking <tools_single_animal_tracking>` option.

#. If using the Debug mode, you will get informed as soon as the
   algorithm detects an object that is moving. It will also inform
   you how much space (in pixels) the detected animal occupies, its
   eccentricity ('roundness') and a parameter for proportions (Major
   axis over minor axis). If the identified object clearly is **not**
   the animal answer the question with "No" and the algorithm will
   look in the next frame the largest moving object.

   .. figure:: Figures/advanced_options/1_debug_detection1.png
      :alt: 1_debug_detection1.png
      :width: 100 %

#. Next, you will see a side by side comparison of the original
   picture (with a box drawn around the detected animal and the binary
   image you have seen in the previous popup.

   .. figure:: Figures/advanced_options/2_debug_detection2.png
      :alt: 2_debug_detection2.png
      :width: 100 %

#. The algorithm will then start tracking. You will see an overview
   of how the algorithm detects the animal: On the left you can see
   the original image. In the center you can see the binary image:
   The grey area indicates the search box (which depends on defined
   max speed of animal, pixel/mm and framerate) and in white the
   pixels that are below threshold. The black area is not considered
   as it is too fare away from the position of the animal in the
   previous frame. On the right, you can see the result of the
   tracking: A box drawn around the identified animal. In addition,
   you can see the animal parameters you are looking for. These are
   just for your information, read below to see how to comfortably
   get the list of these parameters.

   .. figure:: Figures/advanced_options/3_debug_tracking.png
      :width: 100 %
      :alt: 3_debug_tracking.png

#. After running the Single Animal Tracking algorithm, you will find
   a number of new files in each experimental folder. To get to the
   animal parameters, open the file "DATE_TIME_heuristics.csv", for
   example with excel.

   .. figure:: Figures/advanced_options/4_heuristics.png
      :alt: 4_heuristics.png
      :width: 100 %

#. Each row in the table stands for one frame. The title of the
   column describes the value.

   .. figure:: Figures/advanced_options/5_heuristics_columns.png
      :width: 100 %
      :alt: 5_heuristics_columns.png

#. You need to get the following values to get all animal parameters:

   #. A minimum value for filled area (in mm)

   #. A maximum value for filled area (in mm)

   #. A minimum value for eccentricity

   #. A maximum value for eccentricity

   #. A minimum value for major over minor axis

   #. A maximum value for major over minor axis

   #. Maximum skeleton length

   #. Maximum speed of the animal (mm/s)

#. As the tracking algorithm needs the extreme values to function
   properly, I have found it easiest to plot a Line Plot for each
   experiment for each of the relevant parameters. For example for
   the filled area:

   .. figure:: Figures/advanced_options/6_heuristic_plot_example.png
      :alt: 6_heuristic_plot_example.png
      :width: 100 %

#. Write down the maximum and minimum value for each of relevant
   parameters. In this example, the minimum value for filled area in
   mm would be ~25 and the maximum would be ~90.

#. Do the same for eccentricity, major over minor axis, skeleton
   length and maximum speed (mm/s)

#. Then do the same for a few other videos. The goal is to get
   extreme values without having to put 0 as minimum and infinity as
   maximum.

#. In this example I have found the following parameters:

   #. Minimum value for filled area (in mm): 20

   #. Maximum value for filled area (in mm): 90

   #. Minimum value for eccentricity: 0.4

   #. Maximum value for eccentricity: 1

   #. Minimum value for major over minor axis: 1

   #. Maximum value for major over minor axis: 3.5

   #. Maximum skeleton length: 14

   #. Maximum speed of the animal (mm/s): 350

#. Now go to the PiVR software folder on your PC and find the file
   named: "list_of_available_organisms.json":

   .. figure:: Figures/advanced_options/7_open_organism_json.png
      :alt: 7_open_organism_json.png
      :width: 100 %

#. Open it with an text editor. I often use "Code Writer" that ships
   with Windows. You will see that there are repeating structures: A
   word, defining the name of the animal, then a colon and then some
   image parameters in brackets.

   .. note::
      Json files require correct  `formatting <https://en.wikipedia
      .org/wiki/JSON>`__. Be careful to not accidentally deleting
      commas etc.

   .. figure:: Figures/advanced_options/8_open_organism_json.png
      :alt: 8_open_organism_json.png
      :width: 100 %

#. To enter your animal parameters you have two options: The easiest
   (and safest) option is to choose an animal in the list that you
   are certain to never use and just enter your parameters:

   .. figure:: Figures/advanced_options/9_modified_organism_json.png
      :alt: 9_modified_organism_json.png
      :width: 100 %

   Alternatively, you may also enter a new "cell" at the end of the
   list. There is no limit on the number of different animals that
   can be entered in this list.

#. Now save the file (do **not** rename it - If you want to keep a
   backup, rename the original, i.e. to
   "list_of_available_organisms_original.json".)

#. Restart the PiVR software (so that it reads the newly defined
   animal parameters).

#. If you want to know whether PiVR is able to perform real-time
   tracking, you can open the "experiment_settings.json" files in one
   of the video folders you used to find the animal parameters (or a
   newly created video) and change the "Model Organism" cell name to
   your animal name

   .. figure:: Figures/advanced_options/10_change_experiment_settings.png
      :alt: 10_change_experiment_settings.png
      :width: 100 %

#. Now, select the "Debug->Simulate Online Tracking" window, select a
   video and check :ref:`whether the algorithm can track the animal in
   real-time<debug_tracking>`. If not, you might have to select more
   stringent animal parameters and/or you have to optimize imaging
   conditions.

.. _create_own_undistort_files_label:

Create your own undistort files
================================

The stanard camera introduces a fisheye effect as can be observed below:
In reality, the edges of the dish are straight.

   .. figure:: Figures/advanced_options/18_distorted_image.png
      :alt: 18_distorted_image.png
      :width: 100 %

This can lead to problems when collecting data.

PiVR allows for the correction of these optical effects based on a
function in the opencv library.

   .. Note::

      This option can not be turned on if opencv is not installed. If
      the menu is greyed out make sure to install opencv. In addition,
      you will have 'noCV2' written next to the version number of PiVR.

      If you are on the Raspberry Pi the easiest way to install opencv
      is to wipe the SD card, reinstall the OS and make a clean install
      of the PiVR software using the PiVR installation file.

      On a PC, just install it using conda by first (1) activating the
      PiVR environment and (2) then entering `conda install -c conda-forge opencv`

The image below demonstrates the functionality of the undistort
algorithm.

   .. figure:: Figures/advanced_options/19_undistorted_image.png
      :alt: 19_undistorted_image.png
      :width: 100 %

This option is always turned on after v1.7.0 unless you turned it off
:ref:`as described here <UndistortOptionsLabel>`.

Why is it so important to fix the distorted image? PiVR assigns x and y
position of the animal based on the image it gets. If the input image
is distorted these values will be off. For example, in the trajectory
below the x/y positions of the animal differ visibly between the
distorted original image and the undistorted image

   .. figure:: Figures/advanced_options/24_example_trajectory.png
      :alt: 24_example_trajectory.png
      :width: 100 %

In case of presentation of a virtual reality, the arena gets
presented based on the distorted image. This leads to distortion of the
virtual reality.

For the undistort algorithm to work, it requires lens-specific distortion
coefficients. PiVR comes with these coefficients for the standard lens
you get when you buy from the BOM for all standard resolutions (640x480,
1024x768 and 1296x972).

   .. Important::

      If you are using a **different lens** you must create your own
      undistort files. Read on to learn how to do so.

To create your own undistort files please follow the steps below.

   .. Note::

      You will need to conduct this procedure for every resolution you
      want to employ in your experiments.

#. Print this :download:`chessboard <../Chessboard.jpg>` on a piece of
   paper.

#. Go to your PiVR setup and place the printed chessboard on a well lit
   area (for example the light pad).

#. Make sure you have selected the resolution you want to use in the
   future :ref:`Resolution <AdjustImageLabel>`.

#. Open the :ref:`Timelapse Recording Window <TimelapseLabel>` in the
   recording menu.

#. Select a place to save the files. You might want to keep them just
   in case.

#. Set it to record for 1 minute with one image every 2 seconds.

#. Hit start.

#. Take the camera into your hands and take images of the chessboard
   from different angles. See the collage below for an example of the
   different angles you want to get. Try to get the whole chessboard
   into the Field of view

      .. figure:: Figures/advanced_options/20_grabbing_pictures.png
         :alt: 20_grabbing_pictures.png
         :width: 100 %

#. Once you are done, go to 'Tools' -> 'Undistort, new lens'

      .. figure:: Figures/advanced_options/21_undistortnewLens.png
         :alt: 21_undistortnewLens.png
         :width: 100 %

#. Press the 'Chessboard Images' Button and select the folder where you
   saved the chessboard images you just took.

      .. figure:: Figures/advanced_options/22_undistortnewLensOptions.png
         :alt: 22_undistortnewLensOptions.png
         :width: 100 %

#. Everything will now freeze for a couple of minutes. At one point
   you will start to see parts of the chessboard.
   If the images are not good (e.g. because parts of the chessboard are
   missing from the field of view) you will get an error message.
   Please re-take pictures.

#. Once the algorithm is done, you should have a new set of undistort
   coefficent files on your local setup. If you want to make a copy,
   they are in *PiVR/PiVR/undistort_matrices/user_provided*.

#. To use these matrices in your next experiment, press the 'Options'
   Menu in the Menu Bar. Then select 'Undistort Options'.

      .. figure:: Figures/manual_software/UndistortOptions.png
         :width: 100 %
         :alt: UndistortOptions.png

#. In the popup select 'Use your own undistort files'.

      .. figure:: Figures/advanced_options/23_undisortSelectOwnFiles.png
         :alt: 23_undisortSelectOwnFiles.png
         :width: 100 %

#. Save the settings.

#. From now on the output of online tracking is based on your own lens... PiVR documentation master file, created by
   sphinx-quickstart on Tue Oct 16 19:17:07 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#########################################
PiVR: virtual reality for small animals
#########################################

The Raspberry  **Pi** based **V**\irtual **R**\eality system (PiVR) is
a virtual reality system for small animals. It has been developed by
David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

The tool has been published in PLOS Biology: `PiVR: An affordable and
versatile closed-loop platform to study unrestrained sensorimotor
behavior <https://doi.org/10.1371/journal.pbio.3000712>`__.

In addition, Scientific American featured PiVR in a delightful
`article <https://www.scientificamerican.com/article/fruit-flies-plug-into-the-matrix/>`__
accompanied by a fantastic `video <https://youtu.be/71SI7UScBU4>`__.

You can find a 1 hour presentation of PiVR hosted by
`World Wide Open Source <https://www.world-wide.org/Open-Source/>`__
on `youtube <https://www.youtube.com/watch?v=uG898FL421U>`__.


- The code is open source (`BSD license
  <https://choosealicense.com/licenses/bsd-3-clause/>`__)
- The source code for all the PiVR software can be found on `Gitlab
  <https://gitlab.com/LouisLab/PiVR/>`__
- You can also find a `Bug Tracker <https://gitlab
  .com/LouisLab/pivr/issues>`__ on Gitlab

.. figure:: Figures/other/PiVRLogo.jpg
   :width: 75%
   :alt: larva_banana.jpg

   Leonard the larva - always chasing that virtual banana smell.

What can PiVR do?
=================

PiVR has been used to create virtual odor realities for fruit fly larvae.
--------------------------------------------------------------------------

.. figure:: Figures/example_traces/larval_trajectory.png
   :width: 75%
   :alt: larval_trajectory.png

   Trajectory of a *Drosophila* larva in a virtual odor reality. The
   larva expresses the optogenetic tool Chrimson in the *Or42a*
   expressing olfactory sensory neuron.


   .. raw:: html

       <iframe width="75%"
       src="https://www.youtube.com/embed/C9fYz0iN8w0" frameborder="0"
       allow="accelerometer; autoplay;
       encrypted-media; gyroscope; picture-in-picture"
       allowfullscreen></iframe>

|

PiVR has also been used to create virtual taste realities for adult fruit flies.
---------------------------------------------------------------------------------

.. figure:: Figures/example_traces/fly_trajectory.png
   :width: 75%
   :alt: fly_trajectory.png\

   Trajectory of an adult *Drosophila* fly in a virtual taste reality.
   The fly expresses the optogenetic tool Chrimson in the *Gr66a*
   expressing sensory neurons.

.. raw:: html

    <iframe width="75%"
    src="https://www.youtube.com/embed/aow-CmUcT_o" frameborder="0"
    allow="accelerometer; autoplay; encrypted-media;
    gyroscope; picture-in-picture" allowfullscreen></iframe>

|

PiVR was also used to create a virtual light bulb for a number of animals, including larval zebrafish.
-------------------------------------------------------------------------------------------------------

.. figure:: Figures/example_traces/fish_trajectory.png
   :width: 75%
   :alt: fish_trajectory.png

   Trajectory of a zebrafish (D. rerio) larva exposed to a virtual
   white light source.

.. raw:: html

    <iframe width="75%"
    src="https://www.youtube.com/embed/7J5thhZ7Sro" frameborder="0"
    allow="accelerometer; autoplay; encrypted-media; gyroscope;
    picture-in-picture" allowfullscreen></iframe>

|

PiVR is also able to create dynamic virtual gradients
-----------------------------------------------------

While it is often convenient to present static virtual gradients (see
examples above) animals usually have to navigate an environment that
is changing over time. PiVR is able to present animals with dynamic
virtual realities.

We presented *Drosophila* larva expressing the optogenetic tool
Chrimson in the Or42a olfactory sensory neuron with a dynamic odor
plume based on the measurement of a real odor plume (`Álvarez-Salvado et.
al., <https://elifesciences.org/articles/37815>`__). PiVR thus
enables researchers to study how *Drosophila* larvae are orienting
themselves in a more naturalistic environments.

.. raw:: html

    <iframe width="75%"
    src="https://www.youtube.com/embed/RGAta8VqIlw" frameborder="0"
    allow="accelerometer; autoplay; encrypted-media; gyroscope;
    picture-in-picture" allowfullscreen></iframe>

|

How does it work?
=================

PiVR combines high speed/low latency tracking with temporally defined
light stimuli to create arbitrary virtual realities.

.. raw:: html

    <iframe width="75%" height="75%"
   src="https://www.youtube.com/embed/M23uN5JQRG4"
   title="YouTube video player" frameborder="0" allow="accelerometer;
   autoplay; clipboard-write; encrypted-media; gyroscope;
   picture-in-picture" allowfullscreen></iframe>

PiVR is capable of tracking Drosophila larvae, adult
flies, zebrafish and many other small animals.

The position of the animal in the real world is then mapped onto the
user provided virtual reality and the appropriate stimulus is presented.

|

Sounds great. I want one! How?
==============================

PiVR has been designed to make building one as easy as possible so
that you do not spend a lot of time building the setup and spend more
time running experiments.

Please follow the :ref:`PIVR Hardware Manual` to see step-by-step
instructions on how to build your own PiVR setup.

Don't worry, it's not hard and it won't take too long. Please see the
timelapse video below for an example of one setup being built: from
3D printing to running experiments!

.. raw:: html

    <iframe width="75%"
    src="https://www.youtube.com/embed/w5tIG6B6FWo" frameborder="0"
    allow="accelerometer; autoplay; encrypted-media; gyroscope;
    picture-in-picture" allowfullscreen></iframe>

|

I've got a setup. How do I use it?
===================================

If you are a first time user, check out the :ref:`Step-By-Step
<PiVR_Step_by_step_guide>` Guide which will walk you through each of
the four recording modes:

#. :ref:`Tracking single animal <tracking guide>`

#. :ref:`Virtual Reality experiments <VR_guide>`

#. :ref:`Image Sequence Recording <Full_Frame_guide>`

#. :ref:`Video Recording <Video_guide>`

You have just run an experiment. What to make of the output data?
See :ref:`here<Output>` to understand what each output file means and
what it contains.

To see how PiVR can help you analyse data check out the
:ref:`tools<tools>` available on the PC version of PiVR.

Advanced documentation
----------------------

If you are running into trouble with the closed loop tracking, please
head over to :ref:`How to simulate real time tracking <debug_tracking>`.

If you want to track an animal that is not available under
:ref:`Select Animal <SelectOrganismLabel>`, please read the
:ref:`How to define a new animal <define_new_animal_label>` chapter.

If you want to understand what each button in the GUI is doing,
please see the :ref:`PiVR Software Manual <PiVR Software Manual>`.

If you want to gain a high-level understanding on how the code
identifies the animal and tracks them please read the
:ref:`Code Explanation <CodeExplanationLabel>`

The annotated source code can be found :ref:`here <PiVR source code>`
and on the `Gitlab page <https://gitlab.com/LouisLab/PiVR/>`__.


Content
=======

.. toctree::
   :maxdepth: 1
   :numbered:

   self
   manual_hardware
   BOM
   running_experiments
   code_explanation
   manual_software
   output
   tools
   advanced_topics
   software_installation
   software_documentation
   benchmarking
   FAQ
   contact
   PiVR_mentioned
   legacy

.. Indices and tables
   ==================
   * :ref:`intro`
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
   * :ref:`PiVR source code`
   * :ref:`PiVR Software Manual`
   * :ref:`PIVR Hardware Manual`
   * :ref:`PiVR_Step_by_step_guide`
   * :ref:`Output`
   * :ref:`tools`
   * :ref:`advanced_topics`
   * :ref:`software_installation`
   * :ref:`benchmarking`
   * :ref:`FAQ`
   * :ref:`contact`
   * :ref:`PiVR_mentioned`
   * :ref:`legacy`
   Note: I've deactivated this as it doesn't look great



PiVR has been developed by David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

.. Intended Audience: User who just wants to get an experiment started

.. _PiVR_Step_by_step_guide:

Step-by-step experimental guide
********************************

Read this guide if you:

#. Just want to run an experiment

#. Do not (yet) want to go too much into the details of how the
   tracking algorithm works.

Below are the **essential** steps to run an experiment which are
further elaborated on below:

#. Turn the PiVR camera on by pressing ‘Cam On’.

#. Optimize camera height for desired field of view. Typically,
   this adjustment is only done at the beginning of a series of
   experiments.

#. Focus the camera by turning the lens to get sharpest image possible.

#. Set desired resolution in ‘Options->Optimize Image->Camera Resolution’.
   Please see :ref:`here <PiVR_loop_time_high_res>` for a discussion of
   the limitations when using high resolutions.

#. Set desired framerate in ‘Options->Optimize Image->Recording framerate’.

#. Set the backlight intensity to 400,000 in ‘Options->Optimize Image->Backlight Intensity’.

#. Turn autoexposure off by pressing on ‘Options->Optimize Image->autoexp. on’.

#. Enabling the autoexposure can result in conditions where the larva
   cannot be found by the tracker. Disabling this function ensures
   consistency across trials.

#. Increase the backlight intensity a bit to increase contrast.
   Usually going to 500,000 or 600,000 is sufficient.

#. Define camera distance to the object by selecting ‘Options->Define Pixel/mm’.
   This calibration step is essential for successful tracking as it
   converts pixels into physical distance. This only needs to be adjusted
   after (1) changing the camera height or (2) change of resolution.

#. If you want to just track an animal select the Tracking menu
   in ‘Recording->Tracking’. If you want to present a virtual arena
   select the VR Arena menu in 'Recording->VR Arena'.

#. Define folder to save data.

#. Optional: Enter genotype or experimental condition in ‘Exp. Group’.
   This tag is added to the saved data.

#. Select Animal Detection Mode in ‘Options->Animal Detection Method’.

#. Ensure that the correct animal is selected in ‘Options->Select Organisms’.

#. Make sure time is set to desired value.


.. _tracking guide:

Single animal tracking
======================

Press :ref:`here <TrackingLabel>` to learn how to select the Single
Animal tracking menu. Use this option to track single animals. It is
possible to provide a time-dependent stimulation. See :ref:`here
<VR_guide>` to learn how to present a single animal with a virtual
reality.

.. important::

    Several options will open a pop-up. You must close the pop-up in
    order to interact with the main window.

Select Organism
---------------

For proper identification and tracking you first have to define which
organism you are using. See :ref:`here <SelectOrganismLabel>` to
learn how to select your organism.

.. important::
    If you are using an animal that is **not** in the :ref:`List of
    Organisms <SelectOrganismLabel>`, please see
    :ref:`here <define_new_animal_label>`.

Define pixel per mm
-------------------

As you can put the camera at a variety of distances to the arena it
is **essential** to define the pixel per mm ratio. Please see
:ref:`here <PixMMLabel>` to learn how to do that.

.. _SettingUpArenaCamera:

Setting up the arena and the camera
-----------------------------------


It is important to understand how PiVR is able to detect and track an
animal in order to master this method.

The tracking software can only track what is in the field of view of
the camera. If your animal can leave the field of view of the camera,
tracking will stop and you will see an error.

.. figure:: /Figures/running_experiments/ExampleOutsideFOV.png
    :width: 100 %
    :alt: ExampleOutsideFOV.png

    The animal can run outside of the Field of View (FOV) during the
    experiment as the petri dish is not entirely visible using the
    camera. Adjust the dish that constrains animal movement to be in
    the FOV of the camera.

The algorithm works best if the background has as little contrast as
possible

.. figure:: /Figures/running_experiments/ExampleUnevenBackground.png
    :width: 100 %
    :alt: ExampleUnevenBackground.png

    The background is very uneven as the screw of the arena is
    visible (bottom left) and large portions of the image are just
    black while others are white. Adjust the camera and/or the arena
    so that the image only consists of the
    white background illumination (and the animal).

The algorithm works best if the animal has a **high contrast**
relative to the rest of the structure in the image.

.. figure:: /Figures/running_experiments/ExampleContrast.png
    :width: 100 %
    :alt: ExampleContrast.png

    Left: The fly can clearly be seen relative to the background.
    Right: The fly can be seen, but does not have a lot of contrast
    relative to the background. Try to :ref:`improve the image
    <AdjustImageLabel>` so that it looks more like the one on the left.

To set up PiVR to get an image as shown on the left, please follow
:ref:`these Instructions <SetUpOptimalImageLabel>`.

.. _AnimalDetectionGuideLabel:

Animal detection
----------------

The animal needs to be detected **before** the actual data collection
starts. It is also necessary to save a *Background* image for later
use. Ideally this Background image does not contain the animal that
should be tracked.

There are three different animal detection modes
(:ref:`press here to see how to select them <AnimalDetectionLabel>`), each with
it's  own advantages and drawbacks. See
:ref:`here<AnimalDetectionExplanationLabel>` for a an illustration of
the different methods.


#. **Standard - Mode#1:**  This mode will allow you to track a wide
   variety of animals without a lot of optimization of the
   camera or image.

   We have used this mode to track adult flies

   **Advantages**:

   Should work with pretty much any organism if it has been defined
   before.

   Straightforward to run: Place animal, press start tracking, done

   **Disadvantages**:

   The Background image will have a partially visible animal.

   You can't align any virtual reality arenas relative to the  inital
   movement direction of the animal.

   See details :ref:`here <CodeExplanationMode1Label>`

#. **Pre-define Background - Mode#2:** If you need to define a cleaner
   background image, this mode can be useful. You take an image
   before the animal is placed, then add **only** the animal.

   We have used this mode to track zebrafish larvae

   **Advantages**:

   Can track any animal (probably better than Mode#1) if it has been
   defined before.

   Will give a clean background image

   **Disadvantages**:

   As with Mode#1 you can't align any virtual reality arenas relative
   to the inital movement direction of the animal.

   Extremely sensitive: If you have to move anything (such as holding
   up a lid of a petri dish where the animal is supposed to behave)
   this Mode probably won't work well.

   See details :ref:`here <CodeExplanationMode2Label>`

#. **Reconstruct Background by Stitching - Mode#3:**  Should produce the
   same clean background image as Mode#2. Only works if animal
   clearly stands out (has high contrast) in its local environement.

   We have used this mode with fruit fly larvae.

   **Advantages**:

   Allows the usage of aligned virtual reality arenas as the initial
   movement direction of the animal is detected.

   Will give a clean background image.

   **Disadvantages**:

   High contrast requirement hard to fullfill, therefore this mode
   does **not** work well with fast moving animals. Both because
   animals move quickly to the edge and because the area that must be
   taken into account by the detection algorithm increases with the
   speed of the animal.

   Can be difficult to use.

   See details :ref:`here <CodeExplanationMode3Label>`.

**To Summarize**

If you choose Mode 1 or Mode 3:

#. Place the animal in the arena, taking the :ref:`guide above
   <SettingUpArenaCamera>` into account.
#. Press 'Start Tracking'

If you choose Mode 2:

#. Prepare the arena *without* the animal, taking the :ref:`guide
   above <SettingUpArenaCamera>` into account.
#. Press 'Start Tracking'
#. Place the animal into the arena.
#. Press 'Ok'.

Animal tracking
---------------

After successful detection, the tracking algorithm starts following
the animal automatically without the user having to do anything.

While tracking is in progress, the preview window will overlay most
of the monitor and the GUI is not responsive.

There is no status bar available that could be shown during the
tracking.

It is also not possible to cancel the experiment using a
button.

After Animal Tracking
---------------------

Once tracking is finished (either because the animal was tracked for
the defined time or due to an error), the preview window will become
much smaller again.

After saving all the data (which can take a couple of seconds) the
GUI becomes responsive again.

.. _VR_guide:

Single animal Tracking with Virtual Reality
===========================================

Press :ref:`here <VRLabel>` to learn how to select the VR Arena option.

Besides selecting a virtual arena (as described in the :ref:`link
<VRLabel>`) everything is identical to :ref:`Single Animal Tracking
<tracking guide>`.

.. _Full_Frame_guide:

Taking full frame images
========================

Press :ref:`here <FullFrameLabel>` to learn how to select the single
image recording option.

This option is useful if you have several animals that you want to
record simultaneously at a low frame rate. It is possible to provide a
time-dependent stimulation. Compared to the :ref:`video option
<Video_guide>` the resulting images are uncompressed. The
disadvantage is the lower frame rate and the additional hard drive
space necessary to save all the images.

#. Place your arena with the animals you would like to observe into
   the field of view of the camera.
#. Use :ref:`this guide<SetUpOptimalImageLabel>` to get the optimal
   image for your experiment.
#. Press 'Start Recording Images'

.. important::
    This option uses a lot of hard drive space. Make sure there is
    enough space left on your SD card before doing such an experiment.

.. _Video_guide:

Recording a video
=================

Press :ref:`here <VideoLabel>` to learn how to select the Video
Recording option.

This option is useful if you have several animals that you want to
record simultaneously at a high (determined by your resolution, but at
640x480 you should be able to sustain >80 fps) frame rates. It is
possible to provide a time-dependent stimulation. Compared to the
:ref:`full frame option <Full_Frame_guide>` the video will allow for
much higher frame rates while using a fraction of hard drive space.
Videos are encoded in the h264 format. They can be converted into
any other format using ffmpeg_.

 .. _ffmpeg: https://ffmpeg.org/

#. Place your arena with the animals you would like to observe into
   the field of view of the camera.
#. Use :ref:`this guide<SetUpOptimalImageLabel>` to get the optimal
   image for your experiment.
#. Press 'Start Video Recording'
PiVR has been developed by David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

.. _legacy:

Legacy Options
***************

.. warning::

   This page contains information about options that are not recommended
   to be used.

   They might not be supported in future versions of PiVR. Please use
   the alternative.

.. _legacyTimeStimLabel:

Preparing a Frame Based Time Dependent Stimulus File
====================================================

.. warning::

   Please consider using a time based time dependent stimulus file as
   described :ref:`here <PrepTimeStimLabel>`.

.. figure:: Figures/manual_software/TimeDepStimFileOLD.png
   :width: 100 %
   :alt: TimeDepStimFileOLD.png

The first column (A) is the frame number. E.g. if you are recording at
30 frames per second the row 2-32 will define what’s going on in
that time.

The second column defines what Channel 1 is doing at a given frame. 0
means the light is completely OFF. 100 means the light is completely
ON. A number in between, e.g. 50 means that the light is on at
50/100=50%

The third (Channel 2), the fourth (Channel 3) and the fifth (Channel
4) use the same principle for the other channels.

It is important to notice that the stimulation file needs to be
defined on a very low level: Frame Number. The same stimulus file
will give different stimulations depending on the framerate. Therefore:

    1)	Decide on a framerate for you experiment, as an example we’ll
        say you decide on 30fps
    2)	Decide on a length of your experiment, for example 20 seconds
    3)	Decide on the stimulation pattern, e.g. you want Channel 1 to
        be OFF for the first second and
        Channel 2 to be ON for the first second. Then you want to
        switch, Channel 1 is ON for 1 sec, Channel 2 is OFF for 1 sec
    4)	You will need to set the first 30 (framerate * length of
        stimulus) rows of Channel 1 to 0
    5)	And you will need to set the first 30 (framerate * length of
        stimulus) rows of Channel 2 to 100
    6)	As you don’t care about Channel 3 and 4 you can leave it at zero
    7)	At row # 2 (since you start at row #2 in excel) or frame #
        30 (first column) you set Channel 1 to 100 for 30 rows
        (framerate * length of stimulus) to turn it ON and Channel 2
        to 0 to turn it OFF

Notes:
    A)	If you do not define enough rows for your experiment, e.g. if
        you want to run the 20 seconds experiment at 30frames per
        second but you only define what happens during the first 15
        seconds (by only going to row 15*30=450 instead of row
        20*30=600) the last value for each channel will be
        propagated, e.g. if row 450 is set to 100 and row 451 to
        600 are not defined the value 100 will be used for the rest
        of the experiment.
    B)	If you define more rows than you need for your experiments
        only the stimulation up to the point you record are used
        (this will behave as you probably expect)PiVR has been developed by David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

.. _PiVR Software Manual:

PiVR Software Manual
********************

.. warning::

    If you have the High LED power version of PiVR you **must** take
    care to properly shield yourself and others from the potentially
    very strong LED light to protect eyes and skin!

.. important::

    Several options will open a pop-up. You must close the pop-up in
    order to interact with the main window.

.. important::

    The software has different functionality if run on a Raspberry Pi
    as compared to any other PC. This software manual is for the
    Raspberry Pi version of the software


The Menubar
===========

To select a different window use the Menu Bar at the top of the window.

.. figure:: Figures/manual_software/Menubar.png
   :width: 100 %
   :alt: Menubar.png

The Recording Menu
==================

The Recording Menu lets you choose between different recording
options. There are currently 4 different methods:

    1)	Tracking – Online tracking of a single animal. Possibility of
        delivering a time dependent stimulus.
    2)	VR Arena – Online tracking of a single animal. Present a
        virtual arena that will define how the stimulus is present in
        response to the position of the animal.
    3)  Dynamica VR Arena - Online tracking of a single animal. Present
        a virtual arena as above but which changes over time.
    4)	Full Frame Recording – Record an image sequence. Possibility
        of delivering a time dependent stimulus.
    5)  Timelapse Recording - Record a long image sequence at low
        frequency.
    6)	Video – Record a video (h264 format). Possibility of
        delivering a time dependent stimulus.

.. figure:: Figures/manual_software/RecordingMenu.png
   :width: 100 %
   :alt: RecordingMenu.png

Camera Control Frame
--------------------

In all of the recording options you have access to the Camera control
frame. It can be used to turn the camera preview on (Cam On) and off
(Cam Off). You can also control the size of the preview window.

.. warning::

    The Camera preview is always on top of everything else on the
    screen. Use the Preview Window carefully!

.. figure:: Figures/manual_software/CameraControlFrame.png
   :width: 75 %
   :alt: CameraControlFrame.png

.. _TrackingLabel:

Experiment Control Frame – Tracking
-----------------------------------

The ‘Recording’ Option you choose is printed in bold on top of the
Experiment Control Frame. In this example it is ‘Online Tracking’.

.. figure:: Figures/manual_software/ExperimentControlFrameTracking.png
   :width: 100 %
   :alt: ExperimentControlFrameTracking.png

Online tracking tracks a **single** animal.

You have to select a folder in which the experiment will be saved by
clicking on the button to the right of ‘Save in:’.

You can then give your experiment an identifier. Examples include
genotypes or an experimental treatment. This information will be
saved in your experiment folder.

If you want to present a Time Dependent Stimulus you can press the
button ‘Select Time Dependent Stim File’. Please make sure you follow
the :ref:`guidelines <PrepTimeStimLabel>` to learn how to prepare
the file.

The figure below gives you a quick overview over of the parameters used
by the program:

#. Pixel/mm: **Essential**: This value has to be set by you
   before you run your first experiment! See :ref:`set
   Pixel/mm <PixMMLabel>`. You must change it after changing the
   resolution or adjusting the height of the camera relative to the
   arena!

#. Frame rate: The frame rate you will be using to track the animal. See
   :ref:`adjust image <AdjustImageLabel>` to see how to adjust
   the frame rate.

   .. warning::

      There is a difference between the frame rate the camera
      can deliver and the frame rate the Raspberry Pi can handle. If
      you select a very high frame rate you might get a
      lower frame rate than expected. Always check the
      timestamps in the ‘data.csv’ if you are trying a new,
      higher frame rate than before!

#. VR stim at: N/A
#. Animal Detection Mode: Either Mode 1, Mode 2 or Mode 3. See
   :ref:`Select Animal Detection Mod <AnimalDetectionLabel>`.
#. Cam Resolution: Indicates the resolution you selected. See
   :ref:`adjust image <AdjustImageLabel>` to see how to change
   the resolution.

   .. important::

      :ref:`Online Tracking <TrackingLabel>` has only been tested
      with the following resolutions: 640x480, 1024x768, 1296x972.

#. Animal: **Essential**: for :ref:`Online Tracking
   <TrackingLabel>`. See :ref:`here <SelectOrganismLabel>` for
   how to select an animal. See :ref:`Define new
   animal<define_new_animal_label>` in case you are working with an
   animal which is not listed. If you are having problems detecting
   your animal see :ref:`here <debug_tracking>`.

Next, enter the time you want to track the animal in the field
below ‘Recording Time[s]’. Then hit ‘Start Tracking’.

.. _VRLabel:

Experiment Control Frame – VR Arena
---------------------------------------------
The ‘Recording’ Option you choose is printed in bold on top of the
Experiment Control Frame. In this example it is ‘Closed Loop
Stimulation’.

.. figure:: Figures/manual_software/ExperimentControlFrameVRArena.png
   :width: 100 %
   :alt: ExperimentControlFrameVRArena.png

Closed Loop Stimulation tracks a **single** animal.

You have to select a folder in which the experiment will be saved by
clicking on the button to the right of ‘Save in:’.

You can then give your experiment an identifier. Examples include
genotypes or an experimental treatment. This information will be
saved in your experiment folder.

To present a virtual arena (stimulation depending on the position of
the animal) press the ‘Select VR Arena’ button and select an arena.
Static virtual arenas are csv files. Note that you can present the
virtual arena either at a fixed position and independent of the
starting position of the animal (e.g. file "640x480_checkerboard.csv")
**or** you can have the position of the arena defined by the starting
position of the animal (e.g. file
"640x480_gaussian_centred_animal_pos[250, 240,0.0].csv"). See
:ref:`here <create_VR_arena>` for an in-depth explanation.

To learn how to create a new arena please see
:ref:`Create new VR Arena <create_VR_arena>`.

The figure below gives you a quick overview of the parameters used by
the program:

#. Pixel/mm: **Essential**: This value has to be set by you
   before you run your first experiment! See :ref:`set Pixel/mm
   <PixMMLabel>`. You must change it after changing the resolution
   or adjusting the height of the camera relative to the arena!
#. Frame rate: The frame rate you will be using to track the
   animal. See :ref:`adjust image <AdjustImageLabel>` to see how
   to adjust frame rate.

   .. warning::

       There is a difference between the frame rate the camera can
       deliver and the frame rate the Raspberry Pi can handle. If you
       select a very high frame rate you might get a lower frame rate
       than expected. Always check the timestamps in the ‘data.csv’
       if you are trying a new, higher frame rate than before!

#. VR stim at: Either Head, Centroid, Midpoint or Tail. See
   :ref:`here<SelectBodyPart>` how to turn it on.
#. Animal Detection Mode: Either Mode 1, Mode 2 or Mode 3. See
   :ref:`Select Animal Detection Mod <AnimalDetectionLabel>`.
#. Cam Resolution: Indicates the resolution you selected.  See
   :ref:`adjust image <AdjustImageLabel>` to see how to change
   the resolution.

#. Animal: **Essential**: for :ref:`Closed Loop Experiments
   <VRLabel>`. See :ref:`here <SelectOrganismLabel>` for how to
   select an animal. See :ref:`Define new
   animal<define_new_animal_label>` in case you are working with an
   animal which is not listed. If you are having problems detecting
   your animal see :ref:`here <debug_tracking>`

Next, please enter the time you want to track the animal in the field
below ‘Recording Time[s]’. Then hit ‘Start Tracking VR’.

.. _DynamicVRLabel:

Experiment Control Frame – Dynamic VR
---------------------------------------------
The ‘Recording’ Option you choose is printed in bold on top of the
Experiment Control Frame. In this example it is ‘Dynamic VR’.

.. figure:: Figures/manual_software/ExperimentControlDynamicaVR.png
   :width: 100 %
   :alt: ExperimentControlDynamicaVR.png

Dynamic VR tracks a **single** animal.

You have to select a folder in which the experiment will be saved by
clicking on the button to the right of ‘Save in:’.

You can then give your experiment an identifier. Examples include
genotypes or an experimental treatment. This information will be
saved in your experiment folder.

To present a dynamic virtual arena (stimulation depending on the
position of the animal) press the ‘Select VR Arena’ button and select
an arena.
Dynamic virtual arenas are npy files. See
:ref:`here <create_dynamic_VR_arena>` for an in-depth explanation and
how to create them.

The figure below gives you a quick overview of the parameters used by
the program:

#. Pixel/mm: **Essential**: This value has to be set by you
   before you run your first experiment! See :ref:`set Pixel/mm
   <PixMMLabel>`. You must change it after changing the resolution
   or the adjusting height of the camera relative to the arena!
#. Frame rate: The frame rate you will be using to track the
   animal. See :ref:`adjust image <AdjustImageLabel>` to see how
   to adjust the frame rate.

   .. warning::

       There is a difference between the frame rate the camera can
       deliver and the frame rate the Raspberry Pi can handle. If you
       select a very high frame rate you might get a lower frame rate
       than expected. Always check the timestamps in the ‘data.csv’
       if you are trying a new, higher frame rate than before!

#. VR stim at: Either Head, Centroid, Midpoint or Tail. See
   :ref:`here<SelectBodyPart>` how to turn it on.
#. Animal Detection Mode: Either Mode 1, Mode 2 or Mode 3. See
   :ref:`Select Animal Detection Mod <AnimalDetectionLabel>`.
#. Cam Resolution: Indicates the resolution you selected.  See
   :ref:`adjust image <AdjustImageLabel>` to see how to change
   the resolution.

#. Animal: **Essential**: for :ref:`Closed Loop Experiments
   <VRLabel>`. See :ref:`here <SelectOrganismLabel>` for how to
   select an animal. See :ref:`Define new
   animal<define_new_animal_label>` in case you are working with an
   animal which is not listed. If you are having problems detecting
   your animal see :ref:`here <debug_tracking>`.

Next, enter the time you want to track the animal in the field
below ‘Recording Time[s]’. Then hit ‘Start Tracking, dynamic VR’

.. _FullFrameLabel:

Experiment Control Frame – Full Frame Recording
-----------------------------------------------
The ‘Recording’ Option you choose is printed in bold on top of the
Experiment Control Frame. In this example it is ‘Image Sequence’.

.. figure:: Figures/manual_software/ExperimentControlFrameImageSequence.png
   :width: 100 %
   :alt: ExperimentControlFrameImageSequence.png

Image Sequence just records still images without tracking anything.
The advantage over video is that no compression of the image data is
done. The disadvantage is that it is limited by the time it takes the
Raspberry Pi to write the file on the SD card. If you are using a
higher quality SD card, you will be able to write at a higher
frame rate. However, it will probably always be lower than :ref:`video
<VideoLabel>`.

You have to select a folder in which the experiment will be saved by
clicking on the button to the right of ‘Save in:’.

You can then give your experiment an identifier. Examples include
genotypes or an experimental treatment. This information will be
saved in your experiment folder.

If you want to present a Time Dependent Stimulus you can press the
button ‘Select Time Dependent Stim File’. Please make sure you follow
the :ref:`guidelines <PrepTimeStimLabel>` to learn how to prepare
the file.

The figure below gives you a quick overview of the parameters used by
the program:

#. Pixel/mm: This value indicates how many pixels are in one mm.
   You will need this value to be correct to calculate anything
   with distance afterwards (speed, distance to source etc.) See
   :ref:`set Pixel/mm <PixMMLabel>`. You must change it after
   changing the resolution or adjusting the height of the camera
   relative to the arena!
#. Frame rate: The frame rate at which you will be collecting images. See
   :ref:`adjust image <AdjustImageLabel>` to see how to adjust
   the frame rate.

   .. warning::

      There is a difference between the framerate the camera
      can deliver and the framerate the Raspberry Pi can handle.
      If you select a very high framerate you might get a lower
      framerate than expected. Always check the timestamps in
      the ‘data.csv’ if you are trying a new, higher framerate
      than before!

#. VR stim at: N/A.
#. Animal Detection Mode: N/A.
#. Cam Resolution: Indicates the resolution you selected. See
   :ref:`adjust image <AdjustImageLabel>` to see how to change
   the resolution.
#. Animal: Value that will be saved in ‘experiment_settings.json’.

Select the image format you want your images to be in: jpg, png, rbg,
yuv or rgba. See `here <https://picamera.readthedocs.io/en/release-1
.10/api_camera.html#picamera.camera.PiCamera.capture>`__ for details
on the different formats.

Next, please enter the time you want to track the animal in the field
below ‘Recording Time[s]’.

Then hit ‘Start Recording Images'.

.. _TimelapseLabel:

Experiment Control Frame – Timelapse Recording
-----------------------------------------------
The ‘Recording’ Option you choose is printed in bold on top of the
Experiment Control Frame. In this example it is ‘Timelapse’.

.. figure:: Figures/manual_software/ExperimentControlFrameTimelapse.png
   :width: 100 %
   :alt: ExperimentControlFrameTimelapse.png

Timelapse is similar to 'Image Sequence' (See above) in that it
records still images without tracking anything.
In contrast to 'Image Sequence', it allows the taking of pictures at
less than 2 frames per second, the minimal frame rate for all other
modes.

You have to select a folder in which the experiment will be saved by
clicking on the button to the right of ‘Save in:’.

You can then give your experiment an identifier. Examples include
genotypes or an experimental treatment. This information will be
saved in your experiment folder.

   .. Note::

      Please open a ticket on `gitlab <https://gitlab.com/LouisLab/pivr/-/issue>`__
      if you want to be able to present a time dependent stimulus.

The figure below gives you a quick overview of the parameters used by
the program:

#. Pixel/mm: This value indicates how many pixels are in one mm.
   You will need this value to be correct to calculate anything
   with distance afterwards (speed, distance to source etc.) See
   :ref:`set Pixel/mm <PixMMLabel>`. You should change it after
   changing the resolution or adjusting the height of the camera
   relative to the arena!
#. Frame rate: The frame rate the camera is running.
#. VR stim at: N/A.
#. Animal Detection Mode: N/A.
#. Cam Resolution: Indicates the resolution you selected. See
   :ref:`adjust image <AdjustImageLabel>` to see how to change
   the resolution.
#. Animal: Value that will be saved in ‘experiment_settings.json’.

In `Recording Time` indicate the total time you wish to record.

In 'Time between Images' enter the time between frames.

   .. Warning::

      You must make sure that enough space remains on the Raspberry
      Pi. If you run out of space, the program will most likely throw
      an error and stop recording.

Select the image format you want your images to be in: jpg, png, rbg,
yuv or rgba. See `here <https://picamera.readthedocs.io/en/release-1
.10/api_camera.html#picamera.camera.PiCamera.capture>`__ for details
on the different formats.

Then hit ‘Start Timelapse`.

.. _VideoLabel:

Experiment Control Frame – Video
--------------------------------
The ‘Recording’ Option you choose is printed in bold on top of the
Experiment Control Frame. In this example it is ‘Video’.

.. figure:: Figures/manual_software/ExperimentControlVideo.png
   :width: 100 %
   :alt: ExperimentControlVideo.png

As the name indicates, use this option to record videos. The
advantage of this method over image sequence is its superior speed.
The disadvantage, especially for scientific questions, might be that
it compresses the image file in the temporal domain. See `here
<https://www.vcodex.com/an-overview-of-h264-advanced-video-coding/>`__
for an introduction and the Wikipedia page for more details.

You have to select a folder in which the experiment will be saved by
clicking on the button to the right of ‘Save in:’.

You can then give your experiment an identifier. Examples include
genotypes or an experimental treatment. This information will be
saved in your experiment folder.

If you want to present a Time Dependent Stimulus you can press the
button ‘Select Time Dependent Stim File’. Please make sure you follow
the :ref:`guidelines <PrepTimeStimLabel>` to learn how to prepare
the file.

The box below gives you a quick overview over the parameters used by
the program:

#. Pixel/mm: This value indicates how many pixels are in one mm.
   You will need this value to be correct to calculate anything
   with distance afterwards (speed, distance to source etc.) See
   :ref:`set Pixel/mm <PixMMLabel>`. You must change it after
   changing the resolution or adjusting the height of the camera
   relative to the arena!
#. Frame rate: The frame rate at which you will be recording the video. See
   :ref:`adjust image <AdjustImageLabel>` to see how to adjust
   the framerate.

   .. warning::

      There is a difference between the frame rate the camera
      can deliver and the frame rate the Raspberry Pi can handle.
      If you select a very high frame rate you might get a lower
      frame rate than expected. Always check the timestamps in
      the ‘data.csv’ if you are trying a new, higher frame rate
      than before!

#. VR stim at: N/A.
#. Animal Detection Mode: N/A.
#. Cam Resolution: Indicates the resolution you selected. See
   :ref:`adjust image <AdjustImageLabel>` to see how
   to change the resolution.

   .. important::

      For :ref:`video <VideoLabel>` you cannot use 2592x1944.

#. Animal: Value that will be saved in ‘experiment_settings.json’.

Next, please enter the time you want to track the animal in the field
below ‘Recording Time[s]’. Then hit ‘Start Recording Images'.

.. _PrepTimeStimLabel:

Preparing a Time Dependent Stimulus File
=========================================

In your PiVR folder you can find a folder called
‘time_dependent_stim’. On a fresh install it is supposed to contain a
single file: blueprint_stim_file.csv.

When you open it with, e.g. excel or your csv viewing program of
choice, you'll see that there are 6 columns and many rows:

.. figure:: Figures/manual_software/TimeDepStimFile.png
   :width: 100 %
   :alt: TimeDepStimFile.png

The first column (A) is just an index and not really important. The
second column (B) indicates the time at which the stimulus defined
in the columns labelled 'Channel 1', 'Channel 2', 'Channel 3' and
'Channel 4' is being presented. See :ref:`here <DefineGPIOsLabel>` what
a 'Channel' is.

0 means the light is completely OFF. 100 means the light is completely
ON. A number in between, e.g. 50, means that the light is on at
50/100=50%.

You may use the provided file as a blueprint to create your own
stimulus by adding the stimulus intensity at the desired timepoint.
Note that the stimulus must be between 0 and 100.

Alternatively, you can create another file from
scratch. It is important that the file is a csv file with the identical
column names as provided in the file above.

You can change the time resolution if you wish.

.. Important::

   What is a good time resolution to program into the time dependent
   stimulus file? It depends:

   Internally, PiVR keeps track of time using timestamps from the
   camera.
   It then calls `numpy.searchsorted
   <https://numpy.org/doc/stable/reference/generated/numpy.searchsorted.html>`__
   on the provided 'Time [s]' column.

   The algorithm is fast but at a low time resolution can lead to
   unexpected results as it will always stop whenever it finds a value
   larger than the one it looks for.

   For example, if you provide one timepoint for each second while
   recording at 10 frames per second for the first frame at t=0 it will
   present the stimulus for t=0. At t=0.1s it will already provide
   stimulus defined at 1second.

   A good compromise between precision and file size for e.g. 10 frames
   per second is a resolution of 0.01 seconds (10ms). If you want to
   use higher frame rates AND you need very precise stimuli you should
   increase the resolution to 0.001 seconds (1ms). Anything above is not
   useful considering that PiVR can't run at frequencies above 90 Hz
   (about 10ms per frame).

.. note::

   Before v1.7.0, Time Dependent Stimulus File was defined based on
   frame. The above was implemented to give better control over when
   exactly a stimulus is presented. The previous method could introduce
   incoherence between experiments and it is therefore strongly recommended to
   use the method described above.

   If you must use the frame based Time Dependent Stimulus File you may
   find more information :ref:`here <legacyTimeStimLabel>`.

.. _PixMMLabel:

Set Pixel/mm
======================

In order to set Pixel/mm for your resolution, press the ‘Options’
Menu in the Menu Bar. Then select ‘Define Pixel/mm’.

.. figure:: Figures/manual_software/OptionsDefinePxMm.png
   :width: 100 %
   :alt: OptionsDefinePxMm.png

In the popup window you will see features:

    1)	The resolution you are currently using. The defined value
        will only be valid for this resolution
    2)	The left and right cutoff slider. By moving them you can
        measure the distance.
    3)	A slice of the image taken by the camera. You want to put
        something you can measure horizontally before the camera.
    4)	A text field to enter a length you want to measure.

.. figure:: Figures/manual_software/DistanceConfigurationOverview.png
   :width: 100 %
   :alt: DistanceConfigurationOverview.png

Below is an example of an adjusted distance configuration window.
Once you are satisfied with the adjustments you’ve made, hit the quit
button.

.. figure:: Figures/manual_software/DistanceConfigurationAdjusted.png
   :width: 100 %
   :alt: DistanceConfigurationAdjusted.png

.. _AdjustImageLabel:

Adjust image
============

In order to set any options related to the image, press the ‘Options’
Menu in the Menu Bar. Then select ‘Optimize Image’.

.. figure:: Figures/manual_software/OptionsOptimizeImage.png
   :width: 100 %
   :alt: OptionsOptimizeImage.png

This popup should being used to set up the image in the optimal way:

    1)	Turn the camera on (‘Cam On’) if it’s not on already.
    2)	Adjust the preview size so that you can comfortably see both
        the preview and the popup.
    3)	Set the frame rate as desired.
    4)	Press the ‘Update Preview Framerate’ button.
    5)	Set the resolution you’d like to use for the recording.

        .. important::
            For :ref:`Online Tracking <TrackingLabel>` and
            :ref:`Closed Loop Experiments <VRLabel>` only 640x480,
            1024x764 and 1296x962 have been tested.

    6)	Make sure the autoexposure button says ‘autoexp on’.
    7)	Turn the Backlight Intensity up. It is normal to only see
        something above 150’000. 400’000-500’000
        is often a good value to choose.
    8)	If you have Backlight 2 intensity on one of the GPIOs (see
        :ref:`define GPIO output channels <DefineGPIOsLabel>`) you can
        also adjust Backlight 2 intensity at this point.
    9)	To test your output channels, slide the appropriate slider to
        the right. At the beginning of any experiments, these will be
        turned off again. To keep a stimulus ON for the duration of the
        experiment use the Backlight 2 intensity.

.. figure:: Figures/manual_software/OptimizeImageOverview.png
   :width: 100 %
   :alt: OptimizeImageOverview.png

.. _SetUpOptimalImageLabel:

Set up optimal image
--------------------

In order to set up optimal image parameters I usually do the following:

    #) Turn 'Cam On'.
    #) Set 'autoexp on'.
    #) Pull 'Backlight Intensity' slider all the way to the left
       (Image will be dark).
    #) Now pull the 'Backlight Intensity' slider to the right. As
       soon as I see an image in the camera I go another 100'000 to
       the right - this way I'm not at the lower detection limit of the camera.
    #) Then I turn 'autoexp off'.
    #) Often it can improve the image if I pull the 'Backlight
       Intensity' slider a bit more to the right, effectively
       overexposing the image a bit.

.. _DefineOutputFilesLabel:

Define Output Files
===================

Initial versions of PiVR saved data such as centroid position not only
in the *data.csv* file but also in separate *npy* files.

With version 1.6.9 the goal was to reduce clutter in the experimental
folder. All redundant files are now **not** saved by default.

To keep backward compatibility, this option allows users to save files
explicitly as in the earlier versions of PiVR.

To find the menu, press the 'Options' menu in the Menu Bar. Then select
'Output Files'.

.. figure:: Figures/manual_software/OutputFilesSelection.png
   :width: 100 %
   :alt: OutputFilesSelection.png

The popup will allow you to select any of the previously saved numpy
files:

#. Centroids.npy

#. Heads.npy

#. Tails.npy

#. Midpoints.npy

#. Bounding_boxes.npy

#. Stimulation.npy

.. figure:: Figures/manual_software/OutputFilesSelectionOptions.png
   :width: 100 %
   :alt: OutputFilesSelectionOptions.png

In addition, you have the option to save significant amounts of space by
**not** saving the binary images and the skeletons.

These are still being saved by default as there is currently no way to
create identical files. See
`here <https://gitlab.com/LouisLab/pivr/-/issues/52>`__ for discussion
and examples where it fails.

If you know you won't need the binary images and/or the skeletons you
have the option to select that here.

.. _UndistortOptionsLabel:

Undistort Options
=================

In v1.7.0, the undistort feature was added. See
`here (Gitlab) <https://gitlab.com/LouisLab/pivr/-/issues/64>`__ or
:ref:`here (PiVR.org) <create_own_undistort_files_label>` to see a
detailed explanation of what the problem is and how PiVR is solving it.

To find the menu, press the 'Options' menu in the Menu Bar. Then select
'Undistort Options'.

   .. Note::

      This option cannot be turned on if opencv is not installed. If
      the menu is greyed out make sure to install opencv. In addition,
      you will have 'noCV2' written next to the version number of PiVR.

      If you are on the Raspberry Pi the easiest way to install opencv
      is to wipe the SD card, reinstall the OS and make a clean install
      of the PiVR software using the installation file.

      On a PC, just install it using conda by first (1) activating the
      PiVR environment and (2) entering `conda install -c conda-forge opencv`


.. figure:: Figures/manual_software/UndistortOptions.png
   :width: 100 %
   :alt: UndistortOptions.png

In this menu you can choose to perform undistort during tracking or not.

.. figure:: Figures/manual_software/undistortOptionsPopup.png
   :width: 100 %
   :alt: undistortOptionsPopup.png

If you are **not** using the standard lens that comes with the camera
in the BOM you need to use your own undistort files.

See :ref:`here<create_own_undistort_files_label>` how to create your own
files.

.. _DefineGPIOsLabel:

Define GPIO output channels
===========================

What is a 'Channel'?

There are 4 GPIO’s that can be used to control LEDs: GPIO#18,
GPIO#17, GPIO#27 and GPIO#13. (Side Note: GPIO#18 and
GPIO#13 are special as they are the only ones that are capable of
providing PWM frequencies above 40kHz.)

To give the user maximum flexibility, each of the GPIO's can be
assigned a 'Channel' which can be controlled independently in the
software. This also allows the 'bundling' of GPIO's into Channels.

In order to define GPIO output channels for your resolution, press
the ‘Options’ menu in the Menu Bar. Then select ‘define GPIO output
channels’.

.. figure:: Figures/manual_software/OptionsDefineOutputChannels.png
   :width: 100 %
   :alt: OptionsDefineOutputChannels.png


.. figure:: Figures/manual_software/outputChannelSelection.png
   :width: 100 %
   :alt: outputChannelSelection.png

The images on the far left indicate which of the outputs on the left
of your setups are which GPIO (e.g. the one closest to the LED power
input is GPIO#18).

Channel 1 is always defined as the channel that is used for the
Virtual Arena experiments.

Channel 1, Channel 2, Channel 3 and Channel 4 can be separately
addressed using the time dependent stimulus files.

The standard frequency values are set for the normal PiVR setup
running exclusively with LED strips:

.. list-table:: Standard values
   :header-rows: 1

   * - GPIO #
     - Output Channel
     - PWM Frequency
   * - #18
     - Background
     - 40'000 Hz
   * - #17
     - Channel 1
     - 40'000 Hz
   * - #27
     - Channel 1
     - 40'000 Hz
   * - #13
     - Channel 1
     - 40'000 Hz

If you are building the :ref:`The High Powered Version<PiVR HP Version>`
you have to modify the PWM frequency to match the values in the
`datasheet <https://www.ledsupply.com/content/pdf/MiniPuck_F004.pdf>`__:

.. list-table:: High Powered LED PiVR using MiniPuck
   :header-rows: 1

   * - GPIO #
     - Output Channel
     - PWM Frequency
   * - #18
     - Background
     - 40'000 Hz
   * - #17
     - Channel 1
     - 1'000 Hz
   * - #27
     - Channel 1
     - 1'000 Hz
   * - #13
     - Channel 1
     - 1'000 Hz

.. _DebugModeLabel:

Turn Debug Mode ON/OFF
======================

In order turn debug mode On or Off press ‘Options’ menu in the Menu
Bar. Then go on ‘Turn Debug Mode…’ and select either ‘OFF’ or ‘ON’.

.. figure:: Figures/manual_software/OptionsDebugMode.png
   :width: 100 %
   :alt: OptionsDebugMode.png

.. _AnimalDetectionLabel:

Select Animal Detection Mode
============================

In order define the animal detection method press ‘Options’ menu in
the Menu Bar. Then press ‘Animal Detection Method’.

.. figure:: Figures/manual_software/OptionsAnimalDetectionMethods.png
   :width: 100 %
   :alt: OptionsAnimalDetectionMethods.png

When in either ‘Online Tracking’ or ‘Closed Loop Stimulation' the
animal needs to be detected. There are 3 modes that can be used to
detect the animal. For most cases Mode 1 (Standard) will be fine. If
you need a clear background image consider Mode 2 or Mode 3.

.. figure:: Figures/manual_software/SelectAnimalDetection.png
   :width: 100 %
   :alt: SelectAnimalDetection.png

.. _SelectOrganismLabel:

Select Organism
===============

In order select an organism press ‘Options’ menu in the Menu Bar.
Then go on ‘Select Animal’ and select your animal.

.. figure:: Figures/manual_software/OptionsSelectAnimal.png
   :width: 100 %
   :alt: OptionsSelectAnimal.png


.. _UpdateLabel:

Updating the software
=====================

In order to update the software on your RaspberryPi, press the 'File'
menu in the Menu Bar. Then go on 'Update Software'.

.. note::
    Please make sure you are connected to the Internet when updating.

.. figure:: Figures/manual_software/FileUpdate.png
   :width: 100 %
   :alt: FileUpdate.png

Technicalities:

This will first update our Linux by calling::

    sudo update

Next, it will download the newest version from the `gitlab
<https://gitlab.com/louislab/pivr/>`_ repository by calling::

    git pull

.. _HighLowPowerLabel:

High/Low Power LED switch
==========================
In order to choose between high and low power LED setups press
‘Options’ menu in the Menu Bar. Then go on ‘High Power LEDs’.

.. figure:: Figures/manual_software/OptionsHigLowPowerLEDSwitch.png
   :width: 100 %
   :alt: OptionsHigLowPowerLEDSwitch.png

Select either Standard or High power version depending on the setup
you have.

.. _SelectBodyPart:

Select Body Part for VR stimulation
===================================

When running virtual reality experiments the cells you are interested
in could be at different places of the animal.

PiVR allows you to present the virtual reality depending on
:ref:`different body parts<OutputPointExp>` identified during
tracking.

.. figure:: Figures/manual_software/OptionsVRStimulationPoint.png
   :width: 100 %
   :alt: OptionsVRStimulationPoint.png

You may choose different body parts that are defined during tracking.

.. note::
   As the difference between *centroid* and *midpoint* is not
   straightforward, please see :ref:`here <OutputPointExp>` for an
   explanation.

#. The Head (standard) will probably make a lot of sense in many
   experiments, as a lot of sensory neurons of many animals are
   located there. However, be aware that the Head/Tail classification
   algorithm is not perfect and does make mistakes. There is no
   option to correct for wrong head/tail assignment during the
   experiment!

#. The Centroid is probably the most consistently correct point
   during tracking. Please see
   `here <https://scikit-image.org/docs/dev/api/skimage.measure
   .html#skimage.measure.regionprops>`_ to see how it is defined.

#. The Midpoint is similar to the centroid, but can be different in
   flexible animals such as fruit fly larvae.

#. The tail is the final option to choose from. We have used the
   presentation of the virtual reality based on tail position as a
   control in the past.

.. figure:: Figures/manual_software/VRStimulationPoint_Menu.png
   :width: 100 %
   :alt: VRStimulationPoint_Menu.png

.. _select_animal_color:

Animal Color Selection
=======================

Depending on your experimental setup, the animal can either be *dark
on white* background due to transillumination, or *white on dark*
background due to side illumination.

The standard setting is *dark on white*. If you need to change this
setting, go to Options->Animal Color.

.. figure:: Figures/manual_software/OptionMenuAnimalColor.png
   :width: 100 %
   :alt: OptionMenuAnimalColor.png

Now just press the button above the image that describes your
experiment.

.. figure:: Figures/manual_software/AnimalColorOptions.png
   :width: 100 %
   :alt: AnimalColorOptions.png


































PiVR has been developed by David Tadres and Matthieu Louis (`Louis Lab
<https://labs.mcdb.ucsb.edu/louis/matthieu/>`__).

.. _PiVR source code:

PiVR software documentation
***************************

- :ref:`Graphical User interface<PiVR GUI software_documentation>`
- :ref:`Tracking software<PiVR software_documentation>`
- :ref:`Analysis<PiVR Analysis software_documentation>`
- :ref:`Virtual Arena drawing<VR drawing>`
- :ref:`Image Data Handling<PiVR Image data handling>`

.. _PiVR GUI software_documentation:

PiVR GUI source code
=====================

This page contains the classes used to construct the graphical user
interface (GUI).

.. autoclass:: start_GUI.PiVR
    :members:

.. autoclass:: start_GUI.DynamicVirtualRealityFrame
    :members:

.. autoclass:: start_GUI.TrackingFrame
    :members:

.. _PiVR software_documentation:

PiVR Tracking source code
=========================

Detection
---------
.. autoclass:: pre_experiment.FindAnimal
    :members:

Tracking
--------

.. autoclass:: control_file.ControlTracking
    :members:

.. autoclass:: fast_tracking.FastTrackingControl
    :members:

.. autoclass:: fast_tracking.FastTrackingVidAlg
    :members:

Detection and Tracking Helpers
------------------------------

.. autoclass:: tracking_help_classes.FindROI
    :members:

.. autoclass:: tracking_help_classes.MeanThresh
    :members:

.. autoclass:: tracking_help_classes.CallImageROI
    :members:

.. autoclass:: tracking_help_classes.CallBoundingBox
    :members:

.. autoclass:: tracking_help_classes.DescribeLargestObject
    :members:

.. autoclass:: tracking_help_classes.DrawBoundingBox
    :members:

.. autoclass:: tracking_help_classes.Save
    :members:

DefineOutputChannels
--------------------

.. autoclass:: output_channels.DefineOutputChannels
    :members:

Error Messages
--------------

.. autofunction:: tracking_help_classes.show_vr_arena_update_error

.. _VR drawing:

Virtual Arena drawing
======================

.. autoclass:: VR_drawing_board.VRArena
    :members:

.. autoclass:: VR_drawing_board.GaussianSlope
    :members:

.. autoclass:: VR_drawing_board.Step
    :members:

.. autoclass:: VR_drawing_board.PlaceAnimal
    :members:

.. _PiVR Analysis software_documentation:

PiVR Analysis source code
=========================

.. autoclass:: analysis_scripts.AnalysisDistanceToSource
    :members:

.. autoclass:: analysis_scripts.AnalysisVRDistanceToSource

.. autoclass:: multi_animal_tracking.MultiAnimalTracking
    :members:

.. _PiVR Image data handling:

PiVR Image Data Handling source code
====================================

.. autoclass:: image_data_handling.PackingImages
    :members:

.. autoclass:: image_data_handling.ConvertH264
    :members:

