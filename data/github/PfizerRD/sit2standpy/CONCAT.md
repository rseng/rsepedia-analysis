# Sit2StandPy

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02449/status.svg)](https://doi.org/10.21105/joss.02449)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3988351.svg)](https://doi.org/10.5281/zenodo.3988351)

``Sit2StandPy`` is an open source Python package that uses novel algorithms to first detect Sit-to-Stand transitions 
from lumbar-mounted accelerometer data, and then provide quantitative metrics assessing the performance of the 
transitions. A modular framework is employed that would allow for easy modification of parts of the algorithm to suit 
other specific requirements, while still keeping core elements of the algorithm intact. As gyroscopes impose a 
significant detriment to battery life due to power consumption, ``Sit2StandPy``'s use of acceleration only allows for
a single sensor to collect days worth of analyzable data.

## Documentation

[Full documentation](https://sit2standpy.readthedocs.io/en/latest/) is available, containing API references, 
installation instructions, and usage examples.


## Requirements

- Python >= 3.7
- Numpy
- pandas
- Scipy 
- pywavelets
- udatetime

To run the tests, additionally the following are needed

- pytest
- h5py

## Installation

Run in the command line/terminal:

```shell script
pip install sit2standpy
```

pip will automatically collect and install the required packages by default. If you do not want this behavior, run

```shell script
pip install sit2standpy --no-deps
```


## Testing

Automated tests can be run with ``pytest`` through the terminal:

```shell script
pytest --pyargs sit2standpy.tests -v
```

To run the v2 interface tests:
```shell script
pytest --pyargs sit2standpy.v2.tests -v
```

## V2 Interface

Starting with version 1.1.0 a new "v2" interface is available alongside the old interface. Following a sequential
pipeline layout, a basic usage example is:

```python
import sit2standpy as s2s

# transform the data into the appropriate format for H5 or dictionary
# note that "data_transform_function" is your own function to achieve the appropriate format
# if you are looking for a quick example data loader function, you can use the one at
# https://gist.github.com/LukasAdamowicz/b8481ef32e4beeb77c80f29f34c8045e
data = <data_transform/loader_function>(acceleration_data)

sequence = s2s.v2.Sequential()
sequence.add(s2s.v2.WindowDays(hours=[8, 20]))  # window the data into days using only the hours from 8:00 to 20:00
sequence.add(s2s.v2.AccelerationFilter())  # Do the initial filtering and processing required
sequence.add(s2s.v2.Detector(stillness_constraint=True))  # Detect the transitions using the stillness constraint

sequence.predict(data)  # predict and save the results into data

s2s.v2.tabulate_results(data, path_to_csv_output, method='stillness')  # tabulate the results to a csv for easy reading
```


## Old Usage

Basic use is accomplished through the ``Sit2Stand`` object:

```python
import sit2standpy as s2s
import numpy as np  # importing sample data
from sys import version_info
if version_info < (3, 7):
    from pkg_resources import resource_filename
else:
    from importlib import resources

# locate the sample data and load it (depending on python version)
if version_info < (3, 7):
    file_path = resource_filename('sit2standpy', 'data/sample.csv')
    data = np.loadtxt(file_path, delimiter=',')
else:
    with resources.path('sit2standpy.data', 'sample.csv') as file_path:
        data = np.loadtxt(file_path, delimiter=',')

# separate the stored sample data
time = data[:, 0]
accel = data[:, 1:]

# initialize the framework for detection
ths = {'stand displacement': 0.125, 'transition velocity': 0.3, 'accel moving avg': 0.15,
                   'accel moving std': 0.1, 'jerk moving avg': 2.5, 'jerk moving std': 3}
sts = s2s.Sit2Stand(method='stillness', gravity=9.84, thresholds=ths, long_still=0.3, still_window=0.3,
                    duration_factor=4, displacement_factor=0.6, lmin_kwargs={'height': -9.5}, power_band=[0, 0.5],
                    power_peak_kwargs={'distance': 128}, power_stdev_height=True)

# run the sit-to-stand detection
SiSt = sts.apply(accel, time, time_units='us')

# print the list of Transition objects, stored as a dictionary with the time they occurred
print(SiSt)
```

`sit_to_stands` is then a dictionary of `Transition` objects containing information about each of the transitions 
detected


## Contributing

Contributions are welcome.  Please see the [contributions](https://github.com/PfizerRD/sit2standpy/blob/master/CONTRIBUTING.md) document for more information


# Contributing

We encourage any and all contributions to the project. 
If you want to contribute, please [fork the master branch](https://github.com/PfizerRD/sit2standpy/fork) and open a pull request with your edits. 
We will reply to pull requests and are more than happy to discuss any problems or address any questions you may have about the project or process.

## Documentation

Improving documentation is a great way to make first contributions. 
If you found the wording of any of the documentation unclear, or the flow of the documentation difficult to follow, 
please feel free to suggest improvements to make the package more user friendly. All documentation should be in the NumPy standard.

## Code

If you wish to add new features or improve existing code, we are happy to review and approve your suggested changes or
updates. Please just ensure that your addition is well explained, well justified, and clearly commented so that its
purpose and function are evident to the reviewer. 

## Issues

If you find any issues, discover any errors, or can not make the package run feel free to [file an issue](https://github.com/elyiorgos/sleeppy/issues/new). 
When filing your issue please be as detailed and explicit as possible, include the libraries you are using and their versions, the operating system you are using, and any other information you feel is relevant. 
Please also include as much information as necessary to reproduce the problem you encountered.
In the event that you find an issue and are able to resolve it on your own, please submit a pull request with the updated code to help improve the package for others.
---
title: 'Sit2StandPy: An Open-Source Python Package for Detecting and Quantifying Sit-to-Stand Transitions Using an Accelerometer on the Lower Back'
tags:
  - Python
  - Biomechanics
  - Digital Biomarkers
  - Digital Medicine
  - Accelerometer
  - Inertial Sensors
  - Wearable Sensors
  - Human Motion
  - Activity Recognition
authors:
  - name: Lukas Adamowicz
    affiliation: 1
  - name: Shyamal Patel
    affiliation: 1
affiliations:
  - name: Pfizer, Inc. 610 Main Street, Cambridge MA, USA 02139
    index: 1
date: 23 June, 2020
bibliography: paper.bib
---


# Summary

Sit-to-stand transitions have been used to assess mobility for a broad range of conditions, such as Parkinson's Disease [@buckley:2008] and knee osteoarthritis [@bolink:2012]. Assessments are typically performed in the clinic using timed performance tests like the timed-up-and-go [@nguyen:2015; @nguyen:2017] and chair stand tests [@guralnik:1994]. While these assessment tools have demonstrated good psychometric properties, they have two key limitations: (1) assessments are performed episodically as they need to be administered by trained examiners, and (2) assessments performed in the clinic might not provide an adequate measure of real-world mobility. Therefore, there is a growing interest [@pham:2018; @martinez-hernandez:2019] in the use of wearable devices to detect sit-to-stand transitions that occur during daily life and quantify the quality of mobility during these transitions. While some commercial wearable sensor systems have sit-to-stand algorithms (e.g., APDM), and some previous papers have released code [@martinez-hernandez:2019], to the authors' knowledge there are no commonly used, open source, and easy to use and configure, implementations of sit-to-stand algorithms available.

![Location of the wearable device on the lower back.\label{fig:sensor_loc}](sensor_location.png)

``Sit2StandPy`` is an open source Python package that implements novel heuristics-based algorithms to detect Sit-to-Stand transitions from accelerometer data captured using a single wearable device mounted on the lower back (\autoref{fig:sensor_loc}). For each detected transition, the package also calculates objective features to assess quality of the transition. Features include duration, maximum acceleration, minimum acceleration, and spectral arc length (SPARC) [@balasubramanian:2015]. While most previous works have focused on either in-clinic applications [@van_lummel:2013; @nguyen:2015; @nguyen:2017] or the use of multiple sensors [@nguyen:2015; @nguyen:2017], the algorithms in this package can handle data collected under free-living conditions as well as prescribed tasks (e.g., 30-second chair stand task), and it uses only the acceleration data from a single sensor on the lower back with unconstrained device orientation. The practicality of a single-accelerometer approach, which affords a long battery life and improved wearability, makes it well suited for long-term, continuous monitoring at home. 

![A high-level illustration of the input, output, and main processing steps of ``Sit2StandPy``.\label{fig:proc_steps}](alg-overview.png)

As illustrated in \autoref{fig:proc_steps}, ``Sit2StandPy`` takes raw accelerometer data with timestamps as input and returns detected sit-to-stand transitions. Data can be windowed by full days, or parts of days can be selected for each window (e.g., window from 08:00 to 20:00). A high-level interface is provided, which allows the user to access all adjustable parameters. Users provide raw data as input and get detected transitions along with the computed features as output. Additionally, the lower-level methods that are called during the detection are available as well for more fine-grained control. 

With this framework, users can control many parameters of transition detection. In addition, the separation of processing steps and modularity of the sub-processes allows for easy customization if desired. For example, users can easily customize functions for generating additional features to assess quality of a transition.

# Availability

``Sit2StandPy`` is distributed under the MIT License and is published on PyPI, the Python Package Index, and can be installed by running the following in the terminal:

```shell-script
pip install sit2standpy  # install with checking for dependencies
# or
# installation without checking for installed dependencies
pip install sit2standpy --no-deps
```

Sit2StandPy requires the following Python packages:

- Python >=3.7
- NumPy - [@numpy]
- SciPy - [@scipy]
- pandas - [@pandas]
- PyWavelets - [@pywavelets]
- udatetime

``Sit2StandPy`` contains example code with sample data in its GitHub repository. Full documentation with usage examples, installation instructions, and API reference are all available at [https://sit2standpy.readthedocs.io/en/latest/](https://sit2standpy.readthedocs.io/en/latest/)

# Acknowledgements
The Digital Medicine & Translational Imaging group at Pfizer, Inc supported the development of this package.

# References
.. sit2standpy usage

=======================================
Usage examples
=======================================

Basic Use (New/v2 API)
----------------------

Basic usage of ``Sit2StandPy`` to detect transitions in the sample data, using the version 2 API

.. code-block:: python

    >>> import sit2standpy as s2s

    >>> # transform the data into the appropriate format for H5 or dictionary
    >>> # note that "data_transform_function" is your own function to achieve the
    >>> # appropriate format
    >>> # if you are looking for a quick example data loader function, you can
    >>> # use the one at
    >>> # https://gist.github.com/LukasAdamowicz/b8481ef32e4beeb77c80f29f34c8045e
    >>> data = <data_transform/loader_function>(acceleration_data)
    >>>
    >>> sequence = s2s.v2.Sequential()
    >>> # window the data into days using only the hours from 8:00 to 20:00
    >>> sequence.add(s2s.v2.WindowDays(hours=[8, 20]))
    >>> # Do the initial filtering and processing required
    >>> sequence.add(s2s.v2.AccelerationFilter())
    >>> # Detect the transitions using the stillness constraint
    >>> sequence.add(s2s.v2.Detector(stillness_constraint=True))
    >>>
    >>> sequence.predict(data)  # predict and save the results into data
    >>>
    >>> # tabulate the results to a csv for easy reading
    >>> s2s.v2.tabulate_results(data, path_to_csv_output, method='stillness')


Basic Use (Old/v1 API)
----------------------

Basic usage of ``Sit2StandPy`` to detect transitions in sample data:

.. code-block:: python

    >>> import sit2standpy as s2s
    >>> import numpy as np  # importing sample data
    >>> from sys import version_info
    >>> if version_info < (3, 7):
    >>>     from pkg_resources import resource_filename
    >>> else:
    >>>     from importlib import resources
    >>>
    >>> # locate the sample data and load it (depending on python version)
    >>> if version_info < (3, 7):
    >>>     file_path = resource_filename('sit2standpy', 'data/sample.csv')
    >>>     data = np.loadtxt(file_path, delimiter=',')
    >>> else:
    >>>     with resources.path('sit2standpy.data', 'sample.csv') as file_path:
    >>>         data = np.loadtxt(file_path, delimiter=',')
    >>>
    >>> # separate the stored sample data
    >>> time = data[:, 0]
    >>> accel = data[:, 1:]
    >>>
    >>> # initialize the framework for detection
    >>> ths = {'stand displacement': 0.125, 'transition velocity': 0.3, 'accel moving avg': 0.15,
    >>>                    'accel moving std': 0.1, 'jerk moving avg': 2.5, 'jerk moving std': 3}
    >>> sts = s2s.Sit2Stand(method='stillness', gravity=9.84, thresholds=ths, long_still=0.3, still_window=0.3,
    >>>                     duration_factor=4, displacement_factor=0.6, lmin_kwargs={'height': -9.5}, power_band=[0, 0.5],
    >>>                     power_peak_kwargs={'distance': 128}, power_stdev_height=True)
    >>>
    >>> # run the sit-to-stand detection
    >>> SiSt = sts.apply(accel, time, time_units='us')
    >>>
    >>> # print the list of Transition objects, stored as a dictionary with the time they occurred
    >>> print(SiSt)

Advanced Examples (Old/v1 API)
------------------------------

Using :meth:`sit2standpy.Sit2Stand` automatically does all the preprocessing, filtering, and sit-to-stand transition
detection for the user. However, this can be broken up into the constituent parts - preprocessing, filtering, and
detection.

.. code-block:: python

    >>> import sit2standpy as s2s
    >>> from packages_for_importing_data import your_import_data_function
    >>>
    >>> # due to the size of multi day files, no samples are provided with sit2standpy
    >>> accel, time = your_import_data_function()
    >>>
    >>> # PREPROCESSING : conversion of timestamps, and windowing the data
    >>> timestamps, dt, acc_win = s2s.process_timestamps(time, accel,
    >>>                                                  time_units='us',  # time is in microseconds since the epoch
    >>>                                                  window=True,  # window the data
    >>>                                                  hours=('08:00', '20:00'))  # use data from 8am to 8pm each day
    >>>
    >>> # setup filter
    >>> afilt = s2s.AccelerationFilter(power_band=[0, 0.5], power_peak_kw={'distance': 128})
    >>> # filter the windowed data, iterating over the days
    >>> filt_acc, rm_acc, power_peaks = {}, {}, {}
    >>> for day in acc_win.keys():  # dictionary of data for each day, since data was windowed
    >>>     filt_acc[day], rm_acc[day], _, power_peaks[day] = afilt.apply(acc_win[day], 1 / dt)
    >>>
    >>> # setup the detection of the transitions
    >>> still_detect = s2s.detectors.Stillness()  # use the default values
    >>>
    >>> sist = {}
    >>> for day in filt_acc.keys():
    >>>     day_sist = still_detect.apply(acc_win[day], filt_acc[day], rm_acc[day], timestamps[day], dt,
    >>>                                   power_peaks[day])
    >>>     sist.update(day_sist)
.. Sit2StandPy installation file

Installation, Requirements, and Testing
=======================================

Installation
------------

Sit2StandPy can be installed by running any of the following in the terminal:

::

    pip install sit2standpy  # install checking for dependencies with pip
    pip install sit2standpy --no-deps  # install without checking for dependencies
    pip install git+https://github.com/PfizerRD/Sit2StandPy  # alternative to pull from the master branch


Requirements
------------
These requirements will be collected if not already installed, and should require no input from the user.

- Python >= 3.7
- NumPy
- SciPy
- pandas
- pywavelets
- udatetime

These are the versions developed on, and some backwards compatibility may be possible.

To run the tests, additionally the following are needed:

- pytest
- h5py

Testing
-------

Automated tests can be run with ``pytest`` through the terminal:

::

    pytest --pyargs sit2standpy.tests -v

.. Sit2StandPy documentation master file, created by
   sphinx-quickstart on Wed Jul 17 16:13:41 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Sit2StandPy: Automated Sit-to-Stand Transition Detection
========================================================

.. image :: https://joss.theoj.org/papers/10.21105/joss.02449/status.svg
   :target: https://doi.org/10.21105/joss.02449
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3988351.svg
   :target: https://doi.org/10.5281/zenodo.3988351

``Sit2StandPy`` is an open source Python package that uses novel algorithms to first detect Sit-to-Stand transitions
from lumbar-mounted accelerometer data, and then provide quantitative metrics assessing the performance of the
transitions. A modular framework is employed that would allow for easy modification of parts of the algorithm to suit
other specific requirements, while still keeping core elements of the algorithm intact. As gyroscopes impose a
significant detriment to battery life due to power consumption, ``Sit2StandPy``'s use of acceleration only allows for
a single sensor to collect days worth of analyzable data.


Validation
----------
``Sit2StandPy`` was validated, and used to compute the results presented in [1]_. Detailed background and
presentation of the methods can be found there as well.

Capabilities
------------

- Automatic day by day windowing, with option to select a specific window of hours per day
- Optional use of parallel processing to decrease run time
- Quantification of detected transitions
- Modification of functions and methods used in framework
- Access to various parameters controlling the performance and function of the detection

License
-------
Sit2StandPy is open source software distributed under the MIT license.

Papers
------
.. [1] L. Adamowicz et al. "Assessment of Sit-to-Stand Transfers During Daily Life Using an Accelerometer on the Lower Back.'' IEEE Journal of Biomedical and Health Informatics. Under Review.
.. [2] L. Adamowicz, S. Patel. "Sit2StandPy: An Open-Source Python Package for Detecting and Quantifying Sit-to-Stand Transitions Using an Accelerometer on the Lower Back." Journal of Open Source Software. 5(52), 2449. Aug 2020. https://doi.org/10.21105/joss.02449


Contents
--------
.. toctree::
   :maxdepth: 3

   installation
   usage
   refv2/index
   ref/index
sit2standpy.v2.Detector
-----------------------

.. currentmodule:: sit2standpy.v2.detectors

.. autoclass:: Detector
    :members:
sit2standpy.v2.AccelerationFilter
---------------------------------

.. currentmodule:: sit2standpy.v2.filters

.. autoclass:: AccelerationFilter
    :members:
sit2standpy.v2.Sequential
-------------------------

.. currentmodule:: sit2standpy.v2.pipeline

.. autoclass:: Sequential
    :members:
sit2standpy.v2.tabulate_results
-------------------------------

.. currentmodule:: sit2standpy.v2.utility

.. autofunction:: tabulate_results
sit2standpy.v2.get_stillness
----------------------------

.. currentmodule:: sit2standpy.v2.utility

.. autofunction:: get_stillness
sit2standpy.v2.WindowDays
------------------------------------

.. currentmodule:: sit2standpy.v2.day_window

.. autoclass:: WindowDays
    :members:
.. sit2standpy api reference

V2 API Reference
================

.. toctree::
  :maxdepth: 3

  pipeline.Sequential
  day_window.WindowDays
  detectors.Detector
  utility.tabulate_results
  utility.mov_stats
  utility.get_stillness
sit2standpy.v2.mov_stats
------------------------

.. currentmodule:: sit2standpy.v2.utility

.. autofunction:: mov_stats
.. sit2standpy detectors
.. currentmodule:: pysit2stand


sit2standpy.detectors
---------------------

.. toctree::

    detectors.Stillness
    detectors.Displacement


sit2standpy.AccelerationFilter
------------------------------

.. currentmodule:: sit2standpy

.. autoclass:: AccelerationFilter
    :members:
sit2standpy.process_timestamps
------------------------------

.. currentmodule:: sit2standpy

.. autofunction:: process_timestamps
.. Sit2StandPy core functions and classes
.. currentmodule:: Sit2StandPy

Sit2StandPy.Sit2Stand
---------------------

.. autoclass:: Sit2Stand
  :members:sit2standpy.TransitionQuantifier
--------------------------------

.. currentmodule:: sit2standpy

.. autoclass:: TransitionQuantifier
    :members:
sit2standpy.detectors.Stillness
-------------------------------

.. currentmodule:: sit2standpy.detectors

.. autoclass:: Stillness
    :members:
.. Sit2StandPy core functions and classes
.. currentmodule:: Sit2StandPy

Sit2StandPy
-----------

.. toctree::

    core.Sit2Stand
    processing.AccelerationFilter
    processing.process_timestamps
    quantify.TransitionQuantifier
    utility.Transition
    utility.mov_stats
sit2standpy.detectors.Displacement
----------------------------------

.. currentmodule:: sit2standpy.detectors

.. autoclass:: Displacement
    :members:
.. sit2standpy api reference

API Reference
=============

.. toctree::
  :maxdepth: 3

  core
  detectors
sit2standpy.Transition
----------------------

.. currentmodule:: sit2standpy

.. autoclass:: Transition
    :members:
sit2standpy.mov_stats
--------------------------------

.. currentmodule:: sit2standpy

.. autofunction:: mov_stats
