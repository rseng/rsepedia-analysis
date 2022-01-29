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
