---
title: 'The DynaGUI package'
tags:
  - Python
  - particle accelerator
  - physics
  - control system
  - dynamic
  - graphical user interface
authors:
  - name: Benjamin Edward Bolling
    orcid: 0000-0002-6650-5365
    affiliation: "1, 2"
date: 1 November 2019
affiliations:
 - name: European Spallation Source ERIC
   index: 1
 - name: MAX IV Laboratory
   index: 2
bibliography: paper.bib
---

# Summary

At large research facilities and industrial complexes, there is a need for control system user interfaces. However, modern facilities also require continuous upgrading, maintenance, and development, which means that also the control systems need to be upgraded. In order to simplify the construction of control systems and diagnostics, the DynaGUI (Dynamic Graphical User Interface) package was created. The main idea of this package is to get rid of the middle-hand coding needed between hardware and the user by supplying the user with a simple GUI toolkit for generating diagnostics- and control-system GUIs in accordance with any user’s need. Initially developed at MAX IV Laboratory [@MAXIVproj], the initial and main users for this package are control rooms at large-scale research facilities, such as particle accelerators.

In order to further enhance the user-friendliness of this package, a simple system of configuration files has been developed, enabling users to configure the applications using any plain text-editor. The DynaGUI package consists of three applications:

- DynaGUI TF, a true/false (boolean) dynamic control system,

- DynaGUI Alarms, a dynamic diagnostics system for continuously monitoring of the values for a set of attributes, and

- DynaGUI NV, a system for observing attributes' numerical-, string-, vector- or waveform values.

Each DynaGUI application's layout is designed for simplicity. At the top of each application is a combobox with the list of attributes defined. Below the combobox is a button called 'Edit DynaGUI', which opens a window for configuring the DynaGUI. Below the 'Edit DynaGUI' button is the dynamic field, in which a dynamic control panel is generated. Below the dynamic control panel field is the DynaGUI status bar, showing the last action carried out, error messages, or if an alarm is active (for DynaGUI Alarms). Below the status bar are the load and save buttons for loading or saving DynaGUI configuration files. These also have a tool-tip function showing the last loaded or saved file (in the current session). The simplest method to launch DynaGUI is via its launcher (Launcher.py), in which the user can select a control system and either directly define the path to the configuration file or browse to it. If the field is left blank, DynaGUI will load with a predefined setup.

DynaGUI TF (True/False) is dynamical in the sense that the user can insert the name of any device's servers and attributes that, in current state, are True or False. The GUI then builds itself up by creating buttons for each device server that will display the boolean of the device in the combobox selected attribute. First, the application will try to connect to the device's server as entered, and then read in the attribute selected in the combobox and paint the button's background colour:  green meaning True, red meaning False, fuchsia meaning attribute not existing or not boolean, and maroon meaning that the device cannot be connected.

The DynaGUI Alarms (Dynamic Alarms GUI) has been designed to monitor numerical values and notify user(s) if a condition is not fulfilled. To edit the alarms GUI, the user has to press Edit DynaGUI to open an edit-panel. In the left window of the edit-panel, the user has to define the list of devices' server domains and signals to monitor. In the right window of the edit-panel, the user has to define descriptions of the different alarm signals as they should appear in the DynaGUI Alarms control panel, as well as the sweep-time (how frequently the system should check the values). The DynaGUI Alarms' dynamic panel can be divided into three columns for each device and attribute it monitors. The first column contains the descriptions of the alarms, the second column contains the numerical values from the last sweep, the fourth column contains the alarm limit, and the third column contains a combobox with two gaps (larger than or smaller than) which points in the way it should between the read value and the alarm limit. For gap criteria that are fulfilled, the background colour of the alarm description is painted lime-green. If a gap criteria is not fulfilled, the computer's speakers will emit an alarm, and the display will show a message in the DynaGUI Alarms' status bar. The  alarm description's background colour becomes red.

DynaGUI NV is the most advanced tool in the package. Each device has two columns. The first column is a button with the server address of the device, whilst the value of the selected attribute is shown in the second column, enabling users to inspect values in a simple and fast manner. DynaGUI NV also allows for launching a controller for the device (using AtkPanel for Tango Controls [@TangoCS]) or for initializing 1D plots or 2D colormaps of the selected attribute for all devices for which the attribute is valid. Each device's control panel button becomes painted in lime-green if the attribute is valid, in fuchsia if the attribute is not existing for the device and maroon if a connection to the device cannot be established. The plot initialization automatically launches for all devices for which the attribute is valid. The plotting tool has been elaborately described in the DynaGUI documentation. Features in the 1D graph tool include  plotting read values in real time, and setting up and plotting equations (or functions) by using NumPy [@numpy] of the read values as well as combining the read values with one another. Examples are shown in Figure 1. The 1D graph tool supports both scalar and vector (or waveform) plotting in real time, and plot data can be saved and loaded.

The package aims to enable users to construct dynamic GUIs for multiple purposes, and the author is open to implementing new functions and control systems on demand. For testing purposes, an artificial control system 'Randomizer' has been created displaying only random values. In this DynaGUI version, the DynaGUI panels are constructed using PyQt [@PyQtReference], with versions 4 and 5 supported. The author is working on fully implementing the PyEPICS [@pyepicsReference] package as well as a new Finance package built on Pandas [@pandasReference] and Matplotlib [@matplotlib] in order to have more sources to monitor live-stream data from and hence to demonstrate the openness that serves as the core of the package.

# Figures

![A dynamic control panel of DynaGUI NV has been configured (top-left), from which a 1D realtime plot has been launched for 4 artificial devices for a made-up attribute (right). Using this tool, two other lines have been set up as functions of two input data streams, with equations defined in the bottom-left figure and then plotted.](figureNV.png)

# Acknowledgements
The author wants to thank Bernhard Meirose at the MAX IV Laboratory for the discussion that inspired and led me into developing this package. The author also recognises and wants to thank Jonas Petersson and Robin Svard at the MAX IV Laboratory for developing the original 2D Spectrogram Application (for monitoring transverse beam position via Beam Position Monitors at the MAX IV Laboratory storage rings). The author also wants to thank Bernhard Meirose at the MAX IV Laboratory for giving the inspiration and awakening the idea to create this package.

# References
[![DOI](https://joss.theoj.org/papers/10.21105/joss.01942/status.svg)](https://doi.org/10.21105/joss.01942)
# DynaGUI README

DynaGUI stands for Dynamic Graphical User Interface and is a method to construct temporary, permanent and/or a set of GUI:s for users in a simple and fast manner. Developed during shift works at a particle accelerator, the initial goal was to fill in some functions that were then missing: Fast dynamic construction of new control system GUI:s for various purposes.

Different devices can have different attributes, depending on what type of device it is. For example, the Beam Position Monitors (BPM:s) of a particle accelerator reveal information about the transverse beam position whilst magnets' power supplies both reads and can set the current set-point of the magnet.

The BPM:s' device servers do, however, not only reveal information about the beam itself but also contain information of how their data should be treated, such as if they may generate an interlock (InterlockEnabled) or if they should have the Automatic Gain Control enabled or disabled. These signals are handled as true or false flags, meaning that a simple GUI can be constructed to read and also control the true or false-flags. This was the initial application in the package, which is the DynaGUI TF. For simplicity reasons, their states are revealed by their colour:
True = Green
False = Red
non-valid attribute = Magenta
Device disconnected = Maroon

For the case of having a GUI that reads numerical values, each device needs two fields: A field with the domain address for the device and another with the numerical value of the defined attribute, which resulted in the second application DynaGUI NV. This application has also evolved to enable 1D and 2D plotting (showing devices along the vertical axis, time along the horizontal axis, and the values using intensity colours).

The third and last application in this package is the DynaGUI Alarms, which allows a user to set up a list of channels returning numerical values continuously and a second list of limit values (or conditions). This application sounds an alarm for the channel with a value not fulfilling the criteria and paints its description background red in the GUI. Since different devices can have different attributes, inside the DynaGUI NV and DynaGUI TF, the devices which do not have the selected attribute obtain a magenta-coloured background.

The package has then evolved to have the ability to analyse data from any file containing plot data and also live-streamed data from other sources, such as the stock market using Pandas. A random data value package has also been implemented for testing and demonstration purposes.

## Installation procedure
In order to use the PyTango package, the TANGO Controls has to be installed, then follow steps b. TANGO Controls can be obtained from https://www.tango-controls.org/downloads/.

If TANGO is not required by the end user, see steps a.

1. In order to setup this package, ensure that Python 3.x (3.7 is recommended) is installed on the computer.
2. Check Python version used with the PIP package manager such that it points to the correct Python version (pip -V).
3. a) Use PIP to install all packages required, see [requirements](requirements.yml), or use conda to create the environment:

    conda env create --file environments.yml

   b) Use PIP to install all packages required, see [requirements with PyTango](requirements_tango.yml), or use conda to create the environment:

    conda env create --file environment_tango.yml

4. If all required Python packages have been successfully installed, the package is ready.

## Getting Started
This section is a tutorial. Begin with Part 1, followed by Parts 2-4.

### Part 1: The Launcher
The package can be launched by executing `python Launcher.py` in a terminal from the package's location. Upon execution, a widget will open that looks like this:
<table>
    <tr>
        <td>
            <img alt="DynaGUI Launcher" src="figureLauncher.png">
        </td>
    </tr>
</table>
Depending on which packages are installed, different buttons will be enabled and disabled. The package that is included with Python is Random and is hence expected to always load properly, and which is hence used in the examples in this section. The ideal procedure would be if this tutorial is completed using the Random package and then repeated using the package DynaGUI is intended to be used with (EPICS, Tango or Finance) with real data acquisition.

### Part 2: The Boolean Controller (also referred to as TF)
Clicking on Boolean Controller beneath the selected package in the launcher will open a widget that looks like this:
<table>
    <tr>
        <td>
            <img alt="DynaGUI TF" src="figureTF.png">
        </td>
    </tr>
</table>
Each button represents an individual device. The status of each device with the attribute selected are indicated in the different colours discussed previously. If the state is 0 (False) or 1 (True), a signal can be sent out to the device such that the state indicated by the attribute is switched to the opposite for the device. Changing attribute can be carried out by selecting a different one from the attributes drop-down menu. Editing the device configuration is possible by clicking on Edit DynaGUI, which opens up a widget with a list of devices and another list of attributes similar to this:
<table>
    <tr>
        <td>
            <img alt="DynaGUI TF" src="figureEditTF.png">
        </td>
    </tr>
</table>

### Part 3: Value Alarms
Clicking on Value Alarms beneath the selected package in the launcher will open a widget that looks like this:
<table>
    <tr>
        <td>
            <img alt="DynaGUI Alarms" src="figureAlarms.png">
        </td>
    </tr>
</table>
Each checkbox represents a device together with an attribute, that is, the full address to a subscriptable object which returns a numerical value. The value next to it is the latest retrieved value of the object. The symbol dropdown menu can be either smaller than or larger than the expected value, which is defined in the numerical inputs to the right. To run, press the bottom button (and then press again to stop). If a row's logics condition is false, the software will notify with sound and mark the row's background red. Similar to the Boolean Controller described in Part 2, the devices can be edited and the alarm timing set by clicking on Edit DynaGUI which opens up a widget similar to this:
<table>
    <tr>
        <td>
            <img alt="DynaGUI TF" src="figureEditAlarms.png">
        </td>
    </tr>
</table>

### Part 4: Plotting (also referred to as NV)
Clicking on Plotting beneath the selected package in the launcher will open a widget that looks like the top-left widget in the figure below:
<table>
    <tr>
        <td>
            <img alt="DynaGUI NV" src="figureNV.png">
        </td>
    </tr>
</table>
This widget is the so-called DynaGUI NV widget and is used for subscribing to different device's various attributes' numerical values (scalars and waveforms/arrays), with its editing widget being equivalent to the one from the Boolean Controller described in Part 2.

From here, pressing the 1D plot button, the selected attribute will be sent for plotting together with the devices that have shown numerical values for the attribute selected. A new widget opens up that prompts the plotting frequency and number of minutes to show in the plot. Pressing Ok will open up a widget that looks like the right widget in the figure above, except for the lines. Press Start Plotting to begin the data acquisition. Clicking on Plot Settings opens up a widget similar to the bottom-left widget. Clicking on Add New Line opens up a new widget prompting for the (legend) name of the new line, followed by a widget prompting for the equation for the new line which can be a function of another or multiple other line(s).

Clicking on the 2D plot button opens up a widget that looks like this:
<table>
    <tr>
        <td>
            <img alt="DynaGUI NV 2D plotting" src="figure2D.png">
        </td>
    </tr>
</table>
This widget shows each device or waveform index along the vertical axis and the time is plotted along the horizontal axis. With this, each pixel in the plot is assigned with its unique numerical value coming from a device and/or waveform at a given point in time. This value is represented in the pixel with the colour defined by the colour-map definition (to the right in the figure above). To see the different values at a point in time, drag the yellow marker line to the position in time and click on Plot Trace, which will open up a widget that looks like this:
<table>
    <tr>
        <td>
            <img alt="DynaGUI NV 2D plotting" src="figure2Dtrace.png">
        </td>
    </tr>
</table>
The vertical axis from the 2D plot was converted to the horizontal axis, and colour intensity converted to the vertical axis.

### Part 5: Getting Started Tutorial Completion
Congratulations, you have now gone through the tutorial of all the applications of the DynaGUI package. Note that all configurations of the different applications can be saved and loaded at a later point in time. More information can be found under the next section User Guide.

## User Guide
The user guide contains more information on how to use the package, see [User Guide](UserGuide.pdf).

## Package Motivation
Package usages include multiple scenarios for simple financial analysis or at research and industrial complexes where a control system is used that can be incorporated within Python. We will address some case studies below for the different applications in this package.

The first package built was DynaGUI TF (TF = True/False) for the Tango control system with the goal to give the user a quick overview of user-defined devices' statuses of various boolean value attributes by giving each device their own button. Using colour coding for each devices' state, DynaGUI TF offers a time-efficient method for checking multiple devices' attribute states. DynaGUI TF was used a lot to check e.g. BPM statuses (e.g. is interlock enabled for them and is automatic gain control enabled). Whether they should be true or false can change over time, and therefore, it was important to check and ensure that all had proper states set.

Since the limitation with the DynaGUI TF was such that it could only read boolean values, the DynaGUI NV was developed such that numerical values could also be read. However, colour-coding cannot be used to show numerical value statuses, meaning that for this GUI, each device gets two objects, with the first being a button showing if the device is connected to and attribute is valid, if device is connected to but the attribute is not valid, or the device is disconnected (using colour-coding green, magenta and maroon, respectively.). The button can be clicked on to launch a control panel for the device or for more information about the device, whilst the second object is a label which shows the actual value (if it has been obtained). This has been used for launching control panels (AtkPanels) which are built on the Tango control system for the devices and for checking if their values are ok or not, such as the water temperatures and water flows for cooling magnets in the accelerator. It can also be used to plot the measured current from Beam Current Monitors and, by setting up user-defined functions, it can e.g. plot the difference in measured current between two (or more) beam charge monitors to monitor the beam current loss between them.

By using DynaGUI NV, the need for an alarm if a numerical value surpasses some limit was realised. Therefore, the DynaGUI Alarms application was realised as an extension to the two former applications, which was used to monitor e.g. transverse beam emittance and beam energy spread at a synchrotron. When the energy spread or emittance became higher than the user-defined limit, an alarm sounded and notified the operators.

As a future development to test the openness of DynaGUI NV and DynaGUI alarms, a "Finance" package was added to them to monitor stock prices: DynaGUI NV which monitors stock prices in a 1D plot with the ability for users to add functions that are applied to realtime stock price data, whilst DynaGUI Alarms sounds an alarm if a stock price becomes higher or lower than a user-defined limit.

## Dependencies
The package depends on multiple Python packages depending on if it is to be used with Tango, EPICS, Finance, Random, or only historical data plotting and browsing.

## History
The package was initially developed in 2019 by Benjamin Bolling during his time as an Accelerator Operator at MAX IV Laboratory. It has since then evolved to its current state as it is today.

## Comparison to other packages
Many graphical user interface (GUI) toolkits exists for building GUI:s with the Python language, e.g. PyQt, Tkinter and wxPython. All three packages, however, requires that the user does now coding and has time to build a GUI. With the idea of this package being its openness and fast-paced construction for GUI:s, a similar package would be PyGTK which, however, requires a little amount of programming. The author is presently unaware of any Python-based dynamic GUI construction package.

## License
This package is intended to be free and open-source. For more information, see [license](LICENCE.txt).
# DynaGUI test procedures
A manual testing procedure for the DynaGUI package.

## Test 0: The DynaGUI Launcher
First test to be done is to ensure that the launcher works as it is supposed to, which can be accomplished by launching a terminal, browsing to the location of the DynaGUI package and then executing | python Launcher.py |. From the launcher window, the three applications DynaGUI TF, DynaGUI Alarms and DynaGUI NV can be launched for the 4 different packages Tango, EPICS, Random and Finance. The Data Viewer package is not yet in a ready stage.


## Test 1: DynaGUI TF
Launch DynaGUI TF from the launcher for Tango (requires PyTango) or EPICS (requires PyEpics). If none of these packages are available and no devices are connectable, the testing package "Random" can be used. To launch DynaGUI TF, select Boolean Controller for the package to test it with.

In the DynaGUI TF, clicking on any item will make it change state if its colour is red or green. Select any attribute in the top combobox to change which attributes to observe (note that for Random it is not possible to view states for different attributes as there are no devices to communicate with).

The DynaGUI TF devices and attributes can be edited by pressing "Edit DynaGUI". Test the GUI editing by adding, editing and/or removing devices and attributes in their respective list. Change maximum number of rows and select or deselect the "Show the Enable All button" checkbox. Press Ok to confirm or Cancel. 

To save a configuration, press "Save" and select location and filename. The file format for DynaGUI TF configurations is '.dg1'. Confirm file configuration saving and loading methods by closing DynaGUI TF, launching DynaGUI TF again and then loading the file that was saved.

Finally, pressing update statuses fetches the current values of all devices having the specified attribute (except for the Random package).


## Test 2: DynaGUI Alarms
Launch DynaGUI Alarms from the launcher for Tango (requires PyTango), EPICS (requires PyEpics), Finance or Random. To launch DynaGUI Alarms, press Value Alarms for the package to test it with. Press the Edit DynaGUI-button to get to the GUI's editing screen, one column has the full device addresses (incl. the attribute) whilst the other one has a description for each respective device address. Hence, ensure both have the same number of rows (otherwise it will not be possible to proceed). Then select number of rows for the GUI and the alarms timer (time between each sweep), and press Ok to proceed or Cancel.

Test selecting and deselecting alarms, and try the Select All button and Unselect All button. Press the bottom button to begin monitoring the alarms with the defined alarm timer defined (green means it is monitoring, red that it is not monitoring). Edit the values next to the alarm descriptions for testing if the alarms are activated or not. The 'larger than'- and 'smaller than'-signs represent how the values should be (for smaller than sign, the read value has to be smaller than the limit, and vice versa). If this condition is not valid, an alarm will sound that reads the alarm description and "in alarm". The last alarm is also shown as a text together with a time stamp. The background colour is coded as following:

- Green: Ok
- Red: Alarm
- Grey: Skipped

The last things to test is saving and loading a DynaGUI Alarms configuration file, which can be done by pressing "Save". Define the filename and its location. Confirm by closing the DynaGUI Alarms application and then loading the saved file by pressing "Load" and browsing to the file.


## Test 3.1: DynaGUI NV
Launch DynaGUI NV from the launcher for Tango (requires PyTango), EPICS (requires PyEpics), Finance or Random. To launch DynaGUI NV, press Plotting for the package to test it with. Begin by repeating the above described test for DynaGUI TF in terms of editing the DynaGUI layout and then saving and loading the configurations. The "Get all attributes" function is currently only available for the Tango and Random packages. For Random, it generates 5 "random" attributes. For Tango, it obtains a list of all attributes for all devices whilst making sure that no duplicates exist. Press Update statuses to fetch the latest values for all devices or stocks.

## Test 3.2: DynaGUI 2D plotting
The next thing to test is the 2D plotting feature, which can be launched by pressing the "2D Plot" button. Define the plotting frequency and number of minutes to show in the spectrogram. There is still a bug with the intensity mapping, so one has to pull e.g. the top bar up or down slightly on the colour map intensity settings on the top-right side. Press start to begin plotting values. To edit the colormap, press Edit CM and select which CM colour to change (CM1-CM3: CM1 = zero-level value (black), CM2 = middle-level value (green), CM3 = high-level value (red)). The intensity setting can be edited by the user (change size of it, where the zero-level value is, and where the top-level value is). Scroll with the mouse hovering above it to zoom in/out the intensity map or the spectrogram. To reset zoom, right-click in the intensity settings or the spectrogram and select all. 

To plot the values at a specific time, press pause, drag the marker (located to the right in the spectrogram) to a point in time, and press plot trace. To plot the change over time rather than the intensity, press "Plotting real values" and it will change its text to "Plotting vs stored". It will then store the current values as reference values and plot all future values in the spectrogram with respect to the reference values (value = read-value - reference value). New reference values can be taken by clicking on "Store current position".

## Test 3.3: DynaGUI 1D plotting
The next thing to test is the 1D plotting feature, which can be launched by pressing the "1D Plot" button. Define the plotting frequency and number of minutes to show in the plot. To start plotting, press "Start Plotting". It will then change text to "Stop Plotting", which has to be clicked in order to pause the plotting and to be able to change plot settings or load plot data. Press Reset Plot to get rid of all plotted data (it will prompt you if you are really sure about this).

The saving and loading plot data needs more work and is not considered as in a ready state yet (and are hence disabled). The 1D/2D/3D-Layout buttons are also disabled and will in the future give the user a possibility to switch directly between a 1D Plot, 2D Plot and 3D Plot layouts.

Under Plot Settings, change the Plotting Frequency, number of minutes for the spectrogram, and the label of the vertical axis. In here, you can add new lines defined by you. Press Add New Line, define the name of it (no spaces and avoid using other characters than underscores),  e.g. myline_1, and press OK. Now define the equation of the new line. To use the values of other lines, type "PV[x]" where x is the number of the line to include (shown in parenthesis, e.g. line 8 means that x should be 8). Example equations to try out:

	2+PV[0]
	PV[1]-PV[2]
	sin(PV[3]-PV[4])

Add as many lines as you want, and press on the "Functions" tab to see all your lines. To remove a line, press "Remove Line" and select the line to remove. To see all mathematical functions available, press "Pre-defined mathematical functions". In a future version, the "user-defined mathematical functions" will also be available (e.g., func1(x) = x^2 + 3*x). Delays can also be defined for the lines (both negative and positive). Note that for user defined functions, the delays of the read-values are not taken into account (they use simply the last value) but will be fixed in a future version. Press Ok to confirm the new plot settings or Cancel. Select if the horizontal axis of the plot should be rescaled or not (e.g. select yes), and if you want the equations/functions to be applied on the previous values in the plot (e.g. select yes). Now start plotting. Press Show legend to see all the lines and also your new lines in the end. Press Hide All to hide all lines (which will change label to Show All that can be pressed to show all the lines). Lines can also be shown or hidden by selecting or unselecting their respective checkboxes. Line colours can be changed by pressing the colour-box next to each line's checkbox (using RGB). Test: Hide all lines. Select to show "myline_1". Press its colour box, type | 0, 255, 255 | and press Ok. The line should now be cyan.





# Welcome!
Thank you for using and/or considering contributing to this open source project. To do a contribution, please follow the following steps:

1. Create your own branch and develop the code and your contribution.
2. If your contribution is useful for the community, proceed with a merge request to the master (community) branch.
3. Document your contribution (API documentation), write a short description of what you have done and its potential use-cases.
4. Start the pull request.

## Issues
If you experience any issue, please submit it as an issue ticket (but please ensure that it will not result in any duplicate) and label it as a 'bug'.

## Feature request
If there is a feature that you would like to see in a later version, please submit it as an issue ticket and label it as an 'enhancement'.

## Contact details
For direct questions or enquiries, you may write an email with the title "DynaGUI" to:
Email: benjaminbolling@icloud.com
