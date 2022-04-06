# pvrpm-lcoe
[![Tests](https://github.com/FSEC-Photovoltaics/pvrpm-lcoe/actions/workflows/tests.yml/badge.svg?branch=master&event=push)](https://github.com/FSEC-Photovoltaics/pvrpm-lcoe/actions/workflows/tests.yml) [![black-checker](https://github.com/FSEC-Photovoltaics/pvrpm-lcoe/actions/workflows/black-checker.yml/badge.svg)](https://github.com/FSEC-Photovoltaics/pvrpm-lcoe/actions/workflows/black-checker.yml)

[![DOI](https://joss.theoj.org/papers/10.21105/joss.04093/status.svg)](https://doi.org/10.21105/joss.04093)

---

### Quick Reference
[Documentation](https://pvrpm.readthedocs.io/en/latest/)

[License](LICENSE)

[Issues](https://github.com/FSEC-Photovoltaics/pvrpm-lcoe/issues/)

[Contributing Guide](.github/CONTRIBUTING.md)

---

Photovoltaic simulation tool using [SAM](https://github.com/NREL/SAM)

Currently using SAM version 2021.12.02

Please see the documentation for installation and usage guides.

*Note:* Currently Py-PVRPM has issues running on macOS, pleas see #49 for updates.

If you would like to contribute to PVRPM, have a question, or need support, please see the contributing guide [here](.github/CONTRIBUTING.md) for guides on submitting these issues.
---
title: 'PyPVRPM: Photovoltaic Reliability and Performance Model in Python'
tags:
  - Python
  - LCOE
  - Photovoltaic
  - Solar
  - Energy
  - Cost model
authors:
  - name: Brandon Silva^[co-first author] # note this makes a footnote saying 'co-first author'
    affiliation: 1
  - name: Paul Lunis^[co-first author] # note this makes a footnote saying 'co-first author'
    affiliation: 1
  - name: Marios Theristis^[corresponding author]
    affiliation: 2
  - name: Hubert Seigneur^[corresponding author]
    affiliation: 1

affiliations:
 - name: University of Central Florida
   index: 1
 - name: Sandia National Laboratories
   index: 2

date: 08 December 2021
bibliography: paper.bib

---

# Summary

The ability to perform accurate techno-economic analysis of solar photovoltaic (PV) systems is essential for bankability and investment purposes. Most energy yield models assume an almost flawless operation (i.e., no failures); however, realistically, components fail and get repaired stochastically. This package, PyPVRPM, is a Python translation and improvement of the Language Kit (LK) based PhotoVoltaic Reliability Performance Model (PVRPM), which was first developed at Sandia National Laboratories in Goldsim software [@PVRPM:2011] [@PVRPM:2012]. PyPVRPM allows the user to define a PV system at a specific location and incorporate failure, repair, and detection rates and distributions to calculate energy yield and other financial metrics such as the levelized cost of energy and net present value [@PVRPM:2017]. Our package is a simulation tool that uses NREL's Python interface for System Advisor Model (SAM) [@SAM:2020] [@PYSAM] to evaluate the performance of a PV plant throughout its lifetime by considering component reliability metrics. Besides the numerous benefits from migrating to Python (e.g., speed, libraries, batch analyses), it also expands on the failure and repair processes from the LK version by including the ability to vary monitoring strategies. These failures, repairs, and monitoring processes are based on user-defined distributions and values, enabling a more accurate and realistic representation of cost and availability throughout a PV system's lifetime.  

# Statement of need

As photovoltaic technology becomes cheaper and more widely used, solutions are needed to accurately assess its performance and monitor it to safeguard performance and reduce downtime. Existing PV performance models assume ideal operation without considering any failures or repairs. In contrast, the value of monitoring at different levels (e.g., inverter, combiner, string, module-level monitoring) is unclear since these depend on the quality of components, financial metrics, and climatic conditions. PVRPM was initially developed in Goldsim [@PVRPM:2011] [@PVRPM:2012] and then adapted in LK script [@osti_1761998] to improve the accuracy of energy yield simulations by inclusion of realistic energy losses, including operations and maintenance (O&M) costs [@PVRPM:EVAL]. However, these models lack many needed features such as efficient matrix operations, command-line interfaces, and external libraries to perform computations.

Furthermore, the LK-based PVRPM does not include simulation of monitoring configurations (e.g., inverter-level vs. string-level monitoring), which PyPVRPM implements. PyPVRPM takes the original logic from PVRPM, refines it, and adds new features to make it faster, easily configurable, and more powerful. This is done using established Python libraries such as NumPy and Pandas. PyPVRPM can also batch simulate in an automated manner on a computer cluster. PyPVRPM can simulate any component of a PV plant from the module level up to the grid while also considering different monitoring scenarios.

Researchers and relevant stakeholders will be able (but not limited) to use this package to a) determine realistic energy yields and revenues, b) examine whether new and improved monitoring solutions provide benefits to PV investments, and c) optimize O&M and end-of-life replacement strategies. Repair and failure rates can also be compared to perform risk assessments depending on how many failures are expected in a given PV system's lifetime. This library enables bankability and optimization studies. As open-source software in Python, PyPVRPM opens up new opportunities to the PV community to accelerate the development of new capabilities within the PhotoVoltaic Reliability Performance Model.

# Acknowledgements

This work was supported by the U.S. Department of Energy's Office of Energy Efficiency and Renewable Energy (EERE) under the Solar Energy Technologies Office Award Number DE-EE0008157.
Sandia National Laboratories is a multimission laboratory managed and operated by National Technology & Engineering Solutions of Sandia, LLC, a wholly-owned subsidiary of Honeywell International Inc., for the U.S. Department of Energy's National Nuclear Security Administration under contract DE-NA0003525. This paper describes objective technical results and analysis. Any subjective views or opinions that might be expressed in the paper do not necessarily represent the views of the U.S. Department of Energy or the United States Government.

# References
# Contributing to PVRPM

---

Thank you for taking time to contribute to PVRPM! Please refer to the following guidelines below before you submit an issue or pull request. For making edits to the PVRPM code, please follow how to install the package in edit mode [here](https://pvrpm.readthedocs.io/en/latest/tutorial_1installation.html).

### Quick Reference
[Documentation](https://pvrpm.readthedocs.io/en/latest/)

[License](../LICENSE)

[Issues](https://github.com/FSEC-Photovoltaics/pvrpm-lcoe/issues/)

---

## Asking Questions
Before asking a question, please make sure to read the [documentation](https://pvrpm.readthedocs.io/en/latest/) carefully, as your answer to basic questions will be there. If you can not find the answer to your question, check all of issues, both open and closed, [here](https://github.com/FSEC-Photovoltaics/pvrpm-lcoe/issues?q=) to see if someone else asked a similar question that was resolved. If you still can not find your question, please create a new [issue](https://github.com/FSEC-Photovoltaics/pvrpm-lcoe/issues/) and use the **question issue template** for submitting a question. Questions should be clear and concise, and contain as much information as possible; the more information provided the better an answer can be.  

## Feature Requests
Feature requests be submitted as issues [here](https://github.com/FSEC-Photovoltaics/pvrpm-lcoe/issues) and should follow the **feature request template**.
Feature requests should be concise, specific, and if possible include examples of how the feature should work. It will then be assigned a priority level by one of the members of the repository. For changes to **documentation only**, you can instead use the **documentation change template** proposing changes to the documentation. A PR can then reference this issue to change the documentation. This includes grammar fixes, spelling mistakes, or rewording.

## Reporting Bugs
Bug reports should also be submitted as issues [here](https://github.com/FSEC-Photovoltaics/pvrpm-lcoe/issues) and follow the **bug report template**. Bug reports need to provide the operating system and hardware specifications of the computer that ran into the bug, full stack trace (use the `--trace` flag when running a simulation to get stack traces), and include a `zip` file containing **your PVRPM `YAML` configuration, `JSON` files for the SAM case, and the weather file used with the simulation.** This will allow others to reproduce and confirm the bug on other systems.

## Contributing Code
Contribution to the code base of PVRPM should be done with pull requests from a forked repository of PVRPM (see [GitHub's pull request guide](https://docs.github.com/en/pull-requests)). Contributions **must solve a current bug or feature request in PVRPM!** Pull requests that do not have an associated bug or enhancement in an open issue [here](https://github.com/FSEC-Photovoltaics/pvrpm-lcoe/issues) will not be merged or considered until one is opened and referenced to the pull request.

#### Code style
PVRPM follows [black](https://github.com/psf/black) code formatter for all Python code. **You must format you code properly using black before submitting a PR!** Your code can be properly formatted by running these commands in your Python environment:

```bash
$ pip install black
$ cd /path/to/pvrpm/repo
$ python -m black -l 120 .
```

#### Testing
Before submitting, you can easily test your new additions by running the integrated tests in PVRPM using `pytest`. To run the tests, please follow the installation for testing [here](https://pvrpm.readthedocs.io/en/latest/tutorial_1installation.html) then run:

```bash
$ pytest
```

Depending on your hardware, tests should take around 20 minutes to complete. If your tests include new features in PVRPM that are not configured in the `pytest` configuration file, you may edit the test simulation file at `tests/integration/case/test.yml` to add new the features you implemented. **All tests must pass before submitting a PR!**

Again, thank you for taking the time to make PVRPM better! **Make sure to follow templates when contributing! Issues not following templates may receive delayed responses or be automatically closed.**
---
name: Documentation Change
about: Request a change to documentation
title: "[DOCs] <title of doc change here>"
labels: documentation
assignees: ''

---

**Provide a link to the documentation page that needs edits**
Can be linked to the readthedocs site or documentation file in the GitHub repository. If requesting a new page, list the page name and description of what the page should contain.

**Describe the changes**
Can be generic as "this needs to be clearer" or specific changes like "the grammar here is incorrect, change it to ..."

If wanting to add a new page, describe what the new page is and why it should be added.
---
name: Feature Request
about: Suggest an enhancement to this project
title: "[FEATURE] <replace with feature title>"
labels: enhancement
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]
Please link any relevant issues or commits to this feature request.

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Provide an example of the feature**
Could be an example output result file, example of the configuration file, analogy to a real world phenomena in PV systems, etc.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: Bug Report
about: Submit a bug in the PVRPM package
title: "[BUG] <please replace this with bug title>"
labels: bug
assignees: brandons209

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Run command '..'
2. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Case file**
Please include a zip file with the issue submission, or link to a zip file containing:
  - PVRPM `yaml` configuration file
  - JSON files representing the SAM case (exported from SAM)
  - Weather file used in the SAM case (location can be found by looking at the weather file in the SAM's GUI for the case)

**System Information (please complete the following information):**
 - OS: [e.g. macOS 10.16, Arch Linux, Windows 10, Ubuntu 20.04]
 - CPU: [e.g. AMD Ryzen 7 2700x] 
 - Python Version: [e.g. 3.8, 3.9, 3.10]
 - PVRPM Version: [e.g. 1.7.1, 1.7.2]

**Stack Trace**
If applicable, provide a full stack trace here. Use the `--trace` flag when running a PVRPM simulation to print stack traces.
```
<stack trace>
```

**Additional context**
Add any other context about the problem here.
---
name: Question
about: Ask a question and get support
title: "[QUESTION] <please fill with brief question title>"
labels: question
assignees: brandons209

---

**Ask the question**
Provide a brief and concise question.

**Where have I looked for the answer?**
List out where you looked for the answer, [e.g. I checked the installation documentation, assumptions, readme, other GitHub issues, etc]

**Describe the question**
Please provide as much detail as possible in what you are asking. If applicable include screenshots, references to other issues, an example that relates to the question, etc.

**Provide SAM case and PVRPM configuration**
If your question relates to running a PVRPM simulation, please provide a zip file or a link to a zip file containing:
  - PVRPM `yaml` configuration file
  - JSON files representing the SAM case (exported from SAM)
  - Weather file used in the SAM case (location can be found by looking at the weather file in the SAM's GUI for the case)
Installation
=================================
.. toctree::
  :hidden:

This document covers installation and setup of the tool. The tool requires you to build a valid case in SAM, so if you haven't already download and install SAM from here: https://sam.nrel.gov/download.html

**Currently, the supported SAM version is 2021.12.02!**

SAM can be installed on Windows, MAC, or Linux.


Installation
--------------
Requires python >= 3.8

Works on Windows and Linux x64 OSes.

**Currently, there are issues running on macOS. Please see this issue for updates:** https://github.com/FSEC-Photovoltaics/pvrpm-lcoe/issues/49


**Recommended** using pip:
(Replace `@master` with the branch release name if you want a release version)

.. code-block:: bash
  :linenos:

  # for latest development branch
  pip install git+https://github.com/FSEC-Photovoltaics/pvrpm-lcoe/@master

  # for specific version
  pip install git+https://github.com/FSEC-Photovoltaics/pvrpm-lcoe/@vx.x.x

Using the wheel file downloaded from https://github.com/FSEC-Photovoltaics/pvrpm-lcoe/releases

.. code-block:: bash
  :linenos:

  pip install wheel
  pip install pvrpm-x.x.x-py3-none-any.whl

Manually:

.. code-block:: bash
  :linenos:

  git clone https://github.com/FSEC-Photovoltaics/pvrpm-lcoe
  cd pvrpm-lcoe
  python setup.py install

If you want to build the documentation:

.. code-block:: bash
  :linenos:

  git clone https://github.com/FSEC-Photovoltaics/pvrpm-lcoe
  cd pvrpm-lcoe
  pip install .[docs]
  cd docs
  make html

If you want to run automated tests (will take a while based on compute power):

.. code-block:: bash
  :linenos:

  git clone https://github.com/FSEC-Photovoltaics/pvrpm-lcoe
  cd pvrpm-lcoe
  pip install .[testing]
  pytest

For setting up the package in edit mode to modify PVRPM for fixing bugs or adding features:

.. code-block:: bash
  :linenos:

  git clone https://github.com/FSEC-Photovoltaics/pvrpm-lcoe
  cd pvrpm-lcoe
  pip install -e .

Edits to the code can be made in the `pvrpm-lcoe` folder and be tested by running PVRPM from the command line or custom wrapper script.
Usage
=================================
.. toctree::
  :hidden:

Command line usage for PVRPM:

.. code-block:: bash
  :linenos:

  $ pvrpm --help
    Usage: pvrpm [OPTIONS] COMMAND [ARGS]...

      Perform cost modeling for PV systems using SAM and PVRPM

    Options:
      --help  Show this message and exit.

    Commands:
      run     Run the PVRPM LCOE cost model for the case
      sim     Load the SAM case and test the basic SAM simulation
      verify  Verify the case and pvrpm configuration files


Verify configuration:
-----------------------------

.. code-block:: bash
  :linenos:

  $ pvrpm verify --help
    Usage: pvrpm verify [OPTIONS] CONFIG

      Verify the case and pvrpm configuration files

    Options:
      --case <path>  Path to directory containing json export from SAM for the
                     case
      --help         Show this message and exit.

  $ pvrpm verify --case /path/to/case/jsons/ /path/to/config.yml
    2022-02-23 13:30:06,142--INFO: Configuration verified successfully!

Run SAM simulation:
-----------------------------
Mainly for debugging purposes.

.. code-block:: bash
  :linenos:

  $ pvrpm sim --help
    Usage: pvrpm sim [OPTIONS] CONFIG

      Load the SAM case and test the basic SAM simulation

      The config YAML file should specify module order of simulation

    Options:
      --case <path>  Path to directory containing json export from SAM for the
                     case
      --verbose      Enable verbosity in SAM simulation
      --help         Show this message and exit.

  $ pvrpm verify --case /path/to/case/jsons/ --verbose /path/to/config.yml
    0.08 %  @ -1
    0.16 %  @ -1
    0.24 %  @ -1
    0.32 %  @ -1
    0.40 %  @ -1
    0.48 %  @ -1
    0.56 %  @ -1
    0.64 %  @ -1
    ...
    99.76 %  @ -1
    99.84 %  @ -1
    99.92 %  @ -1

Run PVRPM simulation:
-----------------------------

.. code-block:: bash
  :linenos:

  $ pvrpm run --help
    Usage: pvrpm run [OPTIONS] CONFIG

    Run the PVRPM LCOE cost model for the case

    The config YAML file should specify module order of simulation

    Options:
    --case <path>                   Path to directory containing json export
                                    from SAM for the case
    --threads <num_threads>         Number of threads to use for paralized
                                    simulations, set to -1 to use all CPU
                                    threads
    --realization, --realizations <num_realizations>
                                    Set the number of realizations for this run,
                                    overrides configuration file value
    --results <path/to/results/folder>
                                    Folder to use for saving the results,
                                    overrides configuration file value
    --trace                         Enable debug stack traces
    --debug INTEGER                 Save simulation state every specified number
                                    of days for debugging. Saved to results
                                    folder specified
    --noprogress                    Disable progress bars for realizations
    --help

  $ pvrpm run --case /path/to/case/jsons/ --threads -1 --realizations 10 --trace /path/to/config.yml
    2022-02-23 13:51:36,250--WARNING: Lifetime daily DC and AC losses will be overridden for this run.
    2022-02-23 13:51:36,250--WARNING: There is a non-zero value in the fixed annual O&M costs input. These will be overwritten with the new values.
    2022-02-23 13:51:36,250--WARNING: Degradation is set by the PVRPM script, you have entered a non-zero degradation to the degradation input. This script will set the degradation input to zero.
    2022-02-23 13:51:36,252--INFO: Running base case simulation...
    2022-02-23 13:52:24,146--INFO: Base case simulation took: 47.89 seconds

    2022-02-23 13:56:35,089--INFO: Generating results...
    2022-02-23 13:56:39,469--INFO: Graphs saved to /path/to/results/
    2022-02-23 13:56:43,264--INFO: Results saved to /path/to/results/
PVRPM Assumptions / Limitations
=====================================

PVRPM makes a few assumptions in order to be able to run a realistic simulation efficiently. Along with these assumptions arise limitations to what the model can simulate realistically.

To calculate LCOE, NPV, and other data, PVRPM uses SAMs simulation capabilities for PV systems. With all of the configuration for failures, repairs, and monitoring, it generates availability for DC and AC power and the OM yearly cost for the system. These three parameters are what SAM uses in its simulation for each realization to calculate all of the output statistics.

With this, everything PVRPM does must boil down to these three parameters. In doing so, some assumptions must be made to reduce computation complexity. The primary assumption is that failures, repairs, and monitoring are stochastic and can be modeled using statistical distribution. The distribution should consider the many factors that go into these areas, like weather, type of modules, quality of equipment and monitoring, etc., since PVRPM does not simulate these real-world phenomena. Alongside this PVRPM  does not account for changes in the system size, it is considered static across the lifetime of the system; that is, the number of components stays the same across the simulation lifetime.

As of right now, failures are assumed to be a total component failure, where availability of the component while failed is zero. In a future release, partial failures will be available to simulate reduced functionality because of a partial failure of the component but still provide greater than 0 availability in this degraded state. It also assumes that repair costs for these failures do not rise with inflation; only the labor rates rise with yearly inflation rates. The yearly inflation rate is also set to be the same each year, and new labor costs are calculated only at the beginning of every year. When components are repaired with a warranty, it is assumed that the new component has the full warranty time defined in that component's configuration. Warranties are also only applied on successful repair, so if a warranty runs out before a repair can occur, that component is repaired out of warranty.

Availability is also based upon the daylight hours for the configured location via the weather file in SAM, except for the grid component level, which should be up 24 hours every day. When availability is lost due to a failure, the availability lost is only the daylight hours lost, not total hours considered failed. PVRPM only considers sun-up hours, not sun-rise or sun-set hours. It also only considers whether the entire hour is sun-up or not, so for example, if most of an hour is sun-up with a sun-set towards the end, the entire hour is still considered sun-up.

During the calculation of when repairs, monitoring, and failures take place, it is assumed that as soon as a failure occurs, either monitoring or repair immediately begins the next day, depending on what is configured for that component. This must be accounted for in the user's configuration of the distributions. If multiple monitoring types are defined for a single component, then the quickest monitoring will override the others. For example, if component-level monitoring is defined for modules and cross-level monitoring from the inverter to the module, whatever has the shortest time to detection. Typically, it would be the component level monitoring.

With the configuration and outputs for the realizations, PVRPM also compares the ``base case`` in its outputs. The ``base case`` is what the system would output if there were no failures for any component on all of the component levels. The ``base case`` only has module degradation throughout the system's lifetime, depending on the configuration (which could be none if no module degradation is configured). For modeling small systems with module degradation and failures, you might see that monthly energy production is lower for the ``base case`` than the realizations since modules perform much more poorly in the later years of the simulation. Modules in the realizations will fail and be repaired, so their degradation will reset and operate at peak efficiency. This constitutes an overall higher monthly energy production across the lifetime than the base case. This is expected behavior for this scenario.

Finally, each realization is independent of the other realizations in a simulation run. A realization is simply a single "run" of the system for the duration of its lifetime, with the provided configurations.
Logic Diagram
=================================

Here is a set of diagrams to understand how the simulation runs. Each module has a set of functions: init, reinit, and update. Init describes the logic to handle the initial creation at the beginning of the simulation. Reinit describes how to reinitialize components that are repaired. Update gives the logic of what happens on each day of the simulation.

.. image:: images/pvrpm_logic.drawio.png
Getting Started
=================================
.. toctree::
  :hidden:

**Make sure to follow** :doc:`Installation <tutorial_1installation>` **and install SAM before continuing here!**

PVRPM enhances SAM's models for PV systems to obtain a more accurate LCOE than its base simulation. It also allows studying the effects of monitoring and repair techniques on PV systems throughout their lifetime.

To get started, create a new SAM case in the SAM GUI. From there, you **must choose the Detailed Photovoltaic model** for PVRPM to work. This is required because PVRPM need's specific parameters that only exist in this model.

Then, choose your financial model. It must be a financial model that supports lifetime losses. Any financial model under the `Detailed Photovoltaic Model` will work **except** the `LCOE Calculator (FCR Method)` and `No Financial Model`. Please read SAM's documentation for help in setting up the case as it goes into more detail on these models.

Once that is set up, you can download the example configuration and modify it as needed. Below explains from start to finish of how to run a simulation with PVRPM.

Exporting the SAM Case
--------------------------------
PVRPM works by reading the information in your SAM case and the PVRPM YAML configuration to properly run the simulation. For the first step, you need to export your SAM case to JSON files, which PVRPM uses to read in your SAM case.

To do this:
1. Open the SAM GUI, then open your ``.sam`` case file.
2. Once it is open, click the drop-down menu for the case, which is located next to the case name on the top bar.
3. Click ``Generate Code``, then ``PySAM JSON`` (**not** ``JSON for Inputs``).
4. A file explorer window will open, allowing you to select a place to save the JSON files. Select an empty folder.

Once this is done, the selected folder will contain a few JSON files. The names of these files will follow the convention of ``case-name_module.json`` where ``case-name`` is the name of the case, and ``module`` is the module that JSON represents. Pay attention to what modules you have; you'll need to know that for the next part. You can remove the ``.h`` and ``.so`` files.

Configuring PVRPM
------------------------

This will go over every configuration option to set up the case study step by step. The example configuration is also heavily commented on to help with the parameters. Also, please study the logic diagram as it can help when setting up the configuration here. *Also note all of the values listed in examples are entirely arbitrary and do not represent a realistic case.*

You can download the example configuration file :download:`here <../pvrpm/config/example_config.yml>` or view the example configuration file :doc:`here <example_pvrpm_config>`.

Run Setup
----------------
Here, set the results folder location. On Windows, use only one backslash to separate directories; you do not need to escape them or spaces. Then set the number of realizations you want to run and the confidence interval for calculating results.

The results folder and number of realizations can be overridden from the command line.

.. code-block:: yaml
  :linenos:

  ### Run Setup ###
  # results folder doesn't need to exist, it will be created if it doesnt already
  # For windows use single backslash. Do not need to escape spaces in folder names
  results_folder: /path/to/results/folder
  num_realizations: 2 # number of realizations to run
  conf_interval: 90 #XX % confidence interval around the mean will be calculated for summary results

Case Setup
----------------
Set up the basics of the case. Set the number of trackers, combiners, and transformers. Set ``num_trackers`` to ``0`` if you are not using trackers. The worst-case tracker can be set to true if you are using trackers. This means that failures in tracking components result in them being stuck in the worst way: the module is pointing the opposite of the sun's travel arc.

.. code-block:: yaml
  :linenos:

  ### Case Setup ###
  num_combiners: 2 # total number of DC combiner boxes
  num_transformers: 1 # total number of Transformers
  num_trackers: 2 # total number of trackers

  ### Financial Inputs ###
  present_day_labor_rate: 100 # dollars per hour
  inflation: 2.5 # in percent

  ### Failure Tracker Algorithm ###
  use_worst_case_tracker: false

Component Level Setup
------------------------
Each component level requires a setup of failures, monitoring, and repairs. This is required for module, string, combiner, inverter, disconnect, transformer, and grid. However, if you are not setting a component level to fail (``can_fail: false``), you can remove the rest of the sections below it. For example, if string's ``can_fail: false``, I can remove the ``failure``, ``monitoring``, and ``repair`` sections.

There are many options for types of failures, monitoring, and the component, alongside various combinations to get different behaviors.

Keep in mind the way these operations are as follows:

1. First, using the distributions defined for failures, monitoring, and repairs, a ``time_to_failure``, ``time_to_detection``, and ``time_to_repair`` is generated for each component.

2. ``time_to_failure`` then counts down to 0. Once a failure occurs, ``time_to_detection`` counts down to 0 (if monitoring is defined). Finally, ``time_to_repair`` counts to 0, which repairs the component and resets these values.

Component Behaviors
########################
``can_fail`` can be set to ``true`` or ``false``, which dictates whether the components in the component level (``module``, ``string``, etc.) can fail. If this is ``false``, then nothing will fall for this level. You can remove the ``failures``, ``repairs``, and ``monitoring`` sections.
``can_repair`` dictates if the components can be repaired at this level. Typically, leave this ``true`` if components can fail.
``can_monitor`` turns on or off component-level monitoring. This signifies some type of monitoring in place; if you want to simulate this type of monitoring, set this to ``true``.

Warranty
###############
Components can be set to have a warranty. Components covered under warranty do not incur repair costs when they are repaired. A repaired component resets the warranty to the specified time in this section. For no warranty, remove this section.

Distributions for Failures, Repairs, Monitoring
############################################################
Every failure, repair, and monitoring mode requires a distribution to be defined that dictates how long until the specified action occurs.

PVRPM has some built-in distributions, where only a mean and standard deviation is needed to model the distribution properly. Under the hood, PVRPM uses ``scipy.stats`` distributions to generate samples. However, ``scipy.stats`` documentation for each function is unclear on how to convert the mean and std into usable values for the distribution, which is why PVRPM will do that for you.
However, not every single ``scipy`` distribution is wrapped by PVRPM. These are the distributions wrapped by PVRPM (use these as the distribution option):

  - exponential
  - normal
  - uniform
  - lognormal
  - weibull

Using these distributions as the ``distribution`` parameter for any failure, repair, or monitoring only requires you to provide the mean and standard deviation **in days**. The ``weibull`` distribution also allows you to give the ``shape`` for this distribution instead of the standard deviation. On the `Wikipedia page <https://en.wikipedia.org/wiki/Weibull_distribution>`_ for the weibull distribution is the parameter ``k`` and ``lambda`` is calculated from the mean. Using the ``std`` option, make sure it is large since weibull distributions have large STDs by design.

If these distributions don't properly model your data, you can use any distribution listed in the `scipy.stats <https://docs.scipy.org/doc/scipy/reference/stats.html>`_ module. The ``distribution`` parameter in the configuration should be set to the function name of the distribution in the ``scipy.stats`` module. The ``parameters`` key will then be keyword arguments to the function. Make sure to carefully read scipy's documentation, as each function is different in how you need to define it. Remember, the samples from the distribution represent the number of days before that event occurs for a component.

.. code-block:: yaml
  :linenos:

  distribution: normal
  parameters: # parameters for distribution chosen above, either mean and std for a built in distribution or kwargs to the scipy function
    mean: 1460 # years converted to days
    std: 365 # days


Failure Setup
###############
PVRPM currently has two failure modes: total failures and concurrent failures. Total failure modes use the shortest time to failure taken from the defined failures as the time it takes for a component to completely fail. Every component gets a different time to failure, depending on the samples drawn from the distribution.

For a total failure setup, the failure requires the distribution, its parameters, the labor time to fix a component, and the cost to repair the component.

Optionally, there are two fraction modes for a failure: ``fraction`` and ``decay_fraction``. Setting the ``fraction`` will tell PVRPM to fail that fraction (between 0 and 1) of components in the component level consistently throughout the simulation. This means PVRPM will maintain ``fraction`` of the components with this failure mode throughout the simulation. Remember, PVRPM will always pick the failure mode in this section **with the shortest time to failure**, so if you set two failures mode, where one is always shorter than the other, then the longer failure mode will never occur, **even if the fraction is defined on the longer failure mode**.
The ``decay_fraction`` also selects ``decay_fraction`` of the components to fail; however, it decays with each failure. If you set ``decay_fraction`` to 0.5, then at first, 50 percent of the components will fail with this failure mode, then 25  percent, then 12.5 percent, etc., until it approaches 0, which in reality would mean the number of failures from this mode would be 0 when ``decay_fraction`` is small enough.

A typical setup of failures is to have a long "end of life" failure with a large time, and failures with shorter time to failures with a ``fraction`` or ``decay_fraction``, so some will fail with the shorter failures, and most will fail with the end of life failure.

Concurrent failures work the same way as above, except each failure mode is counted **concurrently**. This means that failure modes defined as concurrent failures **do not have the shortest time picked among the modes; instead, each failure mode will fail the component independent of each other and the total failure mode**. You can view this mode as "partial failures", where failures of this nature happen more often than total failures but cost less and are faster to repair. You can use ``fraction`` and ``decay_fraction`` here as needed.

A typical setup for concurrent failure modes is to list routine failures every year or two to a ``fraction`` of the components.

Total failure mode chooses the quickest time to failure from the different modes, and concurrent failure modes all operate independently of each other; they fail each component independent of other failures. Further note, **when a component is repaired from a *total failure*, all *concurrent failures* get reset** since this is a full replacement, and the partial failures that affected the old component won't affect the new one.

.. code-block:: yaml
  :linenos:

  failures:
    eol_failures: # this key name can be anything you want
      distribution: normal
      parameters: # parameters for distribution chosen above, either mean and std for a built in distribution or kwargs to the scipy function
        mean: 3650 # years converted to days
        std: 365 # days
      labor_time: 2 # in hours
      cost: 322 # in USD
    routine_failures: # this key name can be anything you want
      distribution: normal
      parameters:
        mean: 365 # mean in days, or you can do 1 / (num_failures / year * 365)
        std: 365
      labor_time: 2 # in hours
      cost: 322 # in USD
      fraction: 0.1 # > 0 and < 1, fraction of these components that are normal failures, maintained throughout the simulation
    defective_failures: # this key name can be anything you want
      distribution: exponential
      parameters:
        mean: 100 # mean in days, or you can do 1 / (num_failures / year * 365)
      labor_time: 2 # in hours
      cost: 322 # in USD
      decay_fraction: 0.2 # > 0 and < 1, fraction of these components that are defective

  concurrent_failures: # this happens all in parallel, independent of each other
    cell_failure: # this key name can be anything you want
      distribution: normal
      parameters: # parameters for distribution chosen above, either mean and std for a built in distribution or kwargs to the scipy function
        mean: 365 # years converted to days
        std: 365 # days
      labor_time: 2 # in hours
      cost: 322 # in USD
      decay_fraction: 0.2
    wiring_failure: # this key name can be anything you want
      distribution: normal
      parameters:
        mean: 365 # mean in days, or you can do 1 / (num_failures / year * 365)
        std: 365
      labor_time: 2 # in hours
      cost: 322 # in USD
      fraction: 0.1 # > 0 and < 1, fraction of these components that are normal failures, maintained throughout the simulation

Repair Setup
###############
Repairs are much more straightforward. They only need the distribution and its parameters defined for every repair mode. You can either have one repair mode that applies to all failures or a repair mode for each failure mode. You also must list repairs for total failures and concurrent failures separately.

.. code-block:: yaml
  :linenos:

  repairs:
    all_repairs: # this key name can be anything you want
      distribution: lognormal
      parameters:
        mean: 60 # in days
        std: 20 # in days

  concurrent_repairs:
    cell_repair: # this key name can be anything you want
      distribution: lognormal
      parameters:
        mean: 7 # in days
        std: 3 # in days

    wire_repair: # this key name can be anything you want
      distribution: lognormal
      parameters:
        mean: 3 # in days
        std: 3 # in days


Monitoring
###############
Multiple monitoring modes are available for components. You can remove any section you are not using. It is also optional; you can disable all monitoring, in which components that fail are immediately repaired. The modes available are:

  - Component Level: monitoring at the level of the component, which usually offers quick time to detection.
  - Cross Level: monitoring done at a higher level to lower-level components. Meaning inverter monitoring string, combiner, etc.
  - Independent: monitoring done independently of any component level, such as drone IR imaging.

Component level monitoring is defined under each component level's configuration. It simply requires distribution and parameters that signify the time to detection in days to detect a failed component in this level.

.. code-block:: yaml
  :linenos:

  monitoring:
    normal_monitoring: # this key name can be anything you want
      distribution: exponential
      parameters:
        mean: 5

Cross-level monitoring is a bit more complex. Alongside the distribution and parameters, some thresholds control how the monitoring works. A ``global_threshold`` option defines the fraction of components in the monitored level **must fail** before monitoring can detect failed components. This can be seen as enough modules to fail before monitoring at the inverter can start detecting those failures. In PVRPM, this is replicated by the ``global_threshold`` must be met before ``time to detection`` counts down. There is also a ``failure_per_threshold``, the fraction of **lower level** components that must fail **per upper-level component**. For example, if monitoring at the string with a ``failure_per_threshold`` of 0.1, then 10 percent of modules under a single string must fail before the string monitoring can detect module failures. Both thresholds can be defined simultaneously, but one must be defined for this monitoring to work.

.. code-block:: yaml
  :linenos:

  component_level_monitoring:
  # lists what monitoring each component level has for levels BELOW it
  # this is for cross level monitoring only, for defining monitoring at each level use the monitoring distributions above
  string: # component level that has the monitoring for levels below
    module: # componenet level that is below's key. same keys used above: module, string, combiner, inverter, disconnect, grid, transformer
      global_threshold: 0.2 # fraction on [0, 1] that specifies how many of this component type must fail before detection can occur.
                             # this means that until this threshold is met, component failures can never be detected
                             # In the simulation, the time to detection doesn't count down until this threshold is met, which at that point the compouning function will be used along with the distribution as normal
      failure_per_threshold: 0.1 # the fraction of components that must fail per string for failures to be detected at that specified string, or if the total number of failures reach the global_threshold above
      distribution: normal # distribution that defines how long this monitoring takes to detect a failure at this level (independent)
                           # the value calculated from this distribution will be reduced by the compounding factor for every failure in this level
      parameters:
        mean: 1200
        std: 365

  combiner:
    string:
        failure_per_threshold: 0.2 # this fraction of strings must fail on a specific combiner for detection for those to start
        # failures are not compounded globally in this case, only per each combiner
        distribution: normal # distribution that defines how long this monitoring takes to detect a failure at this level (independent)
        parameters:
          mean: 1825
          std: 365

Independent monitoring works outside of component levels. It represents monitoring that **detects all failures in any component level** instantly. It can happen statically every set number of days, defined by the ``interval``, or at a threshold of failed components. There are a few ways to define this threshold. First, the threshold can be defined as ``global_threshold``, which works differently than cross-level monitoring. This value is based on the **DC availability**, meaning the power reaching the inverter. This is calculated using the operating strings and combiners to determine how many modules reach the inverter. With this, combiners and strings are weighted higher than module failures.

The other way to define a threshold is more similar to cross-level monitoring. Using the ``failure_per_threshold`` sets a threshold of failed components for **each level** that must be reached before monitoring occurs. This uses OR logic, meaning only one level has to drop below this threshold for the independent monitoring for **all levels**.

Finally, you can combine all these arguments together; ``interval``, ``global_threshold``, and ``failure_per_threshold``.

You must specify the labor time for each independent monitoring defined, which is in hours for other parameters. There is also **an optional distribution and parameters** that can be defined as the ``time_to_detection`` for components under the levels when the independent monitoring occurs. Think of it as the time it takes to complete the independent monitoring. Not setting this means that the ``time_to_detection`` gets set to zero when independent monitoring occurs.

.. code-block:: yaml
  :linenos:

  indep_monitoring:
    drone_ir:  # this name can be anything you want!
      interval: 1095 # interval in days when this static monitoring occurs
      cost: 50000 # cost of this monitoring in USD
      labor_time: 1 # in hours
      distribution: normal
      parameters:
        mean: 14
        std: 5
      levels: # the component levels this detects on
        - module
        - string
        - combiner

    drone_ir2: # list as many static monitoring methods as you want
      interval: 365 # this monitoring will happen every 365 days, alongside the threshold.
      # a indep monitoring triggered by a threshold RESETs the countdown to the interval
      global_threshold: 0.1 # if DC availability drops by this threshold amount, then this indep monitoring will occur
      # DC availability is the DC power reaching the inverter(s), which is affected by combiners, strings, and module failures
      failure_per_threshold: 0.2 # this threshold is PER LEVEL, if the availability of ANY of the defined levels drops by this threshold amount, this indep monitoring will occur
      cost: 100000
      labor_time: 1 # in hours
      levels:
        - module
        - combiner
        - string

    drone_ir3: # list as many static monitoring methods as you want
      failure_per_threshold: 0.2 # this threshold is PER LEVEL, people if the availability of ANY of the defined levels drops by this threshold amount, this indep monitoring will occur
      cost: 1500
      labor_time: 1 # in hours
      levels:
        - module

Running the simulation
------------------------
The example configuration provided shows how all these options are defined; please consult it as necessary.

Now that you have your SAM case JSONs, and your PVRPM configuration, you can run the simulation:

.. code-block:: bash
  :linenos:

  pvrpm run --case /path/to/directory/with/jsons /path/to/pvrpm/config.yaml


You can also parallelize realizations to decrease the overall run time. To use all your CPU cores to run PVRPM:

.. code-block:: bash
  :linenos:

  pvrpm run --case /path/to/directory/with/jsons --threads -1 /path/to/pvrpm/config.yaml

  PVRPM will alert you to unknown keys in your configuration if you misspelled something and tell you any incorrect or missing parameters you may have.

Once the simulation is completed, result graphs and CSV files will be saved to the defined results folder.
PVRPM's documentation
=================================
About
--------
In photovoltaics (PV), the ability to perform an accurate techno-economic analysis is essential. Creating economically efficient PV systems is beneficial to both consumers and producers alike. This package: Python PhotoVoltaic Reliability Performance Model (PyPVRPM), fills this need. PyPVRPM is a simulation tool that uses NREL's SAM software to model the performance of a PV plant throughout its lifetime. It simulates failures, repairs, and monitoring across the lifespan based on user-defined distributions and values. This allows a more accurate representation of cost and availability throughout a lifetime than SAM's base simulation done from the GUI. By varying repair rates and monitoring solutions, one can compare different configurations to find the most optimal setup for implementing an actual PV power plant.

A few assumptions are taken in the tool, alongside specific ways calculations are made. Please see the logic diagram to understand how the simulation works. Also, view the example configuration to get an idea of setting up your case.

PVRPM Requires a valid SAM case created in its GUI before use. The SAM case provides the information for PVRPM to operate, alongside simulations using its LCOE calculator defined in the case, weather files, etc. Please see the getting started tutorials to learn how to set up your SAM case and PVRPM YAML configuration.

.. toctree::
  :maxdepth: 2
  :glob:
  :caption: Tutorials:

  tutorial*

.. toctree::
  :maxdepth: 2
  :glob:
  :caption: Examples:

  example*

.. toctree::
  :maxdepth: 2
  :glob:
  :caption: API:

  api/*

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
Example YAML Configuration
===========================
.. toctree::
  :hidden:

Please read all the comments carefully, there are many options for each section of PVRPM.

.. include:: ../pvrpm/config/example_config.yml
  :code: yaml
pvrpm.core.modules package
==========================

Submodules
----------

pvrpm.core.modules.failure module
---------------------------------

.. automodule:: pvrpm.core.modules.failure
   :members:
   :undoc-members:
   :show-inheritance:

pvrpm.core.modules.monitor module
---------------------------------

.. automodule:: pvrpm.core.modules.monitor
   :members:
   :undoc-members:
   :show-inheritance:

pvrpm.core.modules.repair module
--------------------------------

.. automodule:: pvrpm.core.modules.repair
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: pvrpm.core.modules
   :members:
   :undoc-members:
   :show-inheritance:
pvrpm
=====

.. toctree::
   :maxdepth: 4

   pvrpm
pvrpm package
=============

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   pvrpm.core

Module contents
---------------

.. automodule:: pvrpm
   :members:
   :undoc-members:
   :show-inheritance:
pvrpm.core package
==================

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   pvrpm.core.modules

Submodules
----------

pvrpm.core.case module
----------------------

.. automodule:: pvrpm.core.case
   :members:
   :undoc-members:
   :show-inheritance:

pvrpm.core.components module
----------------------------

.. automodule:: pvrpm.core.components
   :members:
   :undoc-members:
   :show-inheritance:

pvrpm.core.enums module
-----------------------

.. automodule:: pvrpm.core.enums
   :members:
   :undoc-members:
   :show-inheritance:

pvrpm.core.exceptions module
----------------------------

.. automodule:: pvrpm.core.exceptions
   :members:
   :undoc-members:
   :show-inheritance:

pvrpm.core.logger module
------------------------

.. automodule:: pvrpm.core.logger
   :members:
   :undoc-members:
   :show-inheritance:

pvrpm.core.simulation module
----------------------------

.. automodule:: pvrpm.core.simulation
   :members:
   :undoc-members:
   :show-inheritance:

pvrpm.core.utils module
-----------------------

.. automodule:: pvrpm.core.utils
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: pvrpm.core
   :members:
   :undoc-members:
   :show-inheritance:
