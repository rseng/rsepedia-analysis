# Guidelines for Contributors

The MSWH software is open to contributions from the coder community around the world. 

The contributors should comply with the 
[University of California Standards of Ethical Conduct](https://policy.ucop.edu/doc/1100172/EthicalValuesandConduct).

The contribution can be a bug fix, an improvement to an existing feature, or a fully new feature. If you wish to contribute any of those please follow this workflow:

1. Fork the MSWH repository.
1. Create a new issue using the [MSWH issue tracker](https://github.com/LBNL-ETA/MSWH/issues) which briefly describes your intended contribution.
1. Await a response from a person responsible for code maintenance.
1. Discuss whether there is a need to address this issue and ways in which the issue could be best addressed.
1. Once an agreement is reached on how to address the issue, develop the code on a new branch on your forked repository. It is highly recommended that the branch includes the MSWH software repo issue number for easier tracking.
1. When you are happy with your code that addresses the issue, and after you confirm that all tests are passing, create a [pull request](https://github.com/LBNL-ETA/MSWH/compare) against the `master` branch of the MSWH repository.
1. Notify the person that was discussing the issue with you using `@githubname` that the feature/fix/improvement is now ready for review.
1. Conduct iterations of receiving review and addressing it until the reviewer approves the PR.
1. The reviewer will then thank you and merge the PR.
1. Now you can start utilizing the updated code!

In addition, if you have questions on how to best apply the existing code for your particular purpose, you may pose those using the [MSWH issue tracker](https://github.com/LBNL-ETA/MSWH/issues).

Thank you for your interest in contributing to our open source MSWH repository!

Milica
## Copyright Notice

Multiscale Solar Water Heating (MSWH) Copyright (c) 2019, The
Regents of the University of California, through Lawrence Berkeley National
Laboratory (subject to receipt of any required approvals from the U.S.
Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative
works, and perform publicly and display publicly, and to permit other to do
so.

## License Agreement

Multiscale Solar Water Heating (MSWH) Copyright (c) 2019, The
Regents of the University of California, through Lawrence Berkeley National
Laboratory (subject to receipt of any required approvals from the U.S.
Dept. of Energy).  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches,
or upgrades to the features, functionality or performance of the source
code ("Enhancements") to anyone; however, if you choose to make your
Enhancements available either publicly, or directly to Lawrence Berkeley
National Laboratory, without imposing a separate written license agreement
for such Enhancements, then you hereby grant the following license: a
non-exclusive, royalty-free perpetual license to install, use, modify,
prepare derivative works, incorporate into other computer software,
distribute, and sublicense such enhancements or derivative works thereof,
in binary and source code form.
# Multiscale Solar Water Heating
**Solar water heating system modeling and simulation for individual and community scale projects**

## Repository Content

Folder | Content
------ | ------
[mswh](mswh) | Python module to calculate solar irradiation on a tilted surface ([mswh/system/source_and_sink.py](mswh/system/source_and_sink.py)). <br><br> Python module with simplified component models ([mswh/system/components.py](mswh/system/components.py)) for Converter (solar collectors, electric resistance heater, gas burner, photovoltaic panels, heat pump), Storage (solar thermal tank, heat pump thermal tank, conventional gas tank water heater), and Distribution (distribution and solar pump, piping losses) components. <br><br> Python module with preconfigured system simulation models ([mswh/system/models.py](mswh/system/models.py)) for: base case gas tank water heaters, solar thermal water heaters (solar collector feeding a storage tank, with a tankless gas water heater backup in a new installation cases and a base case gas tank water heater in a retrofit case) and solar electric water heaters (heat pump storage tank with an electric resistance backup). <br><br> Database with component performance parameters, California specific weather data and domestic hot water end-use load profiles ([mswh/comm/swh_system_input.db](mswh/comm/mswh_system_input.db)). <br><br> Modules to communicate with the database ([mswh/comm/sql.py](mswh/comm/sql.py)), unit conversion and plotting modules in [mswh/tools](mswh/tools).
[scripts](scripts) | Jupyter notebooks with preconfigured models and any side analysis if applicable. Navigate to scripts in [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/LBNL-ETA/MSWH/632fd9860c66e3d5b5cafe0af61a6b42f4c6b4f7) to try them out quickly.
[web](web) | Django web framework to configure project, parametrize components and run simulation from a web browser.
[docs](docs) | API documentation, including a short methodology documentation can be found [here](https://lbnl-eta.github.io/MSWH/). To build HTML or LaTeX use `make html` or `make latex`. A `pdf` version of the Code Documentation can be viewed and downloaded [here](https://github.com/LBNL-ETA/MSWH/blob/v2.0.0/docs/MSWH.pdf).

## Statement of Need

We envision four main groups of users for the MSWH software:

* Researchers and policy developers.
* Solar water heating planners, designers and contractors.
* Homeowners.
* Educators.

The policy developers and researchers could utilize the existing MSWH software by embedding it into some larger analysis framework they construct such that it provides answers to their specific research questions.

The professional planners, designers, and contractors of solar thermal water heating systems might find it useful to have access to a freely available simulation tool such as the MSWH software, that they can use to evaluate alternative system designs.

Homeowners considering transitioning to a solar water heating system may be interested in doing the math before seeking further professional help, or just for their own education and curiosity about both solar water heating systems and system simulation in general.

Educators may wish and find it useful to utilize the MSWH simulation tool in the classroom when teaching the basics of energy simulation.

## Usage

The fastest way to explore the preset simulations is to use the [`MSWH System Tool`](scripts/MSWH&#32;System&#32;Tool.ipynb) notebook. In the notebook the user provides a climate zone for a project, an occupancy for each household and whether any of the occupants stay at home during the day. The notebook can then load a set of example California specific hourly domestic hot water end-use load profiles from a database, size and locate the systems. The user can now simulate the hourly system performance over a period of one representative year, visualize and explore the simulation results using time-series plots for temperature profiles, heat and power rates, or look at annual summaries. Similarly the user can model individual household solar water heating projects and base case conventional gas tank water heater systems, such that the results can be compared between the individual, community and base case systems. All simulation and sizing parameters are exposed in the notebook and the user can easily change them if needed.

If you opt to use the web framework the shortest path to explore the simulaton results after [setting up a local server](#django-web-framework-deployment) is to:

* Click on `Configurations` on the landing page.
* Click on `Simulate` for any of the example preconfigured systems (`Solar Thermal New` or `Solar Electric`). This leads the user to a visualization page with hourly timeseries results for a representative year.
* Play with sizes and performance parameters of preconfigured components.

To configure new system types in the web framework (such as `Solar Thermal Retrofit`) one would need to map it through the backend analogously to the currently preconfigured systems.

An example demonstrating usage of the simulation models for an additional climate outside
of California, that is Banja Luka in Bosnia & Herzegovina, is provided in [this notebook](scripts/MSWH&#32;System&#32;Tool&#32;-&#32;Additional&#32;Climate.ipynb).

## Setup and Installation

1. Make sure that `pip` [is installed](https://pip.pypa.io/en/stable/installing/).

2. Unless you already have [`conda`](https://docs.conda.io/en/latest/) installed, please install the lightweight option [`Miniconda`](https://docs.conda.io/en/latest/miniconda.html) or [`Anaconda`](https://docs.anaconda.com/anaconda/install/) software.

### Simple Installation Using `Conda`

1. If you are familiar with `conda` and experienced with virtual environments
 you can perform the package installation using the following set of commands:

        conda create -n mswh -c conda-forge -c plotly python=3.8 pip git-lfs jupyterlab plotly-orca
        conda activate mswh
        git lfs install
        git clone https://github.com/LBNL-ETA/MSWH.git
        cd MSWH
        pip install -e .

    To ensure functionality of the example notebooks install the following:

        python -m ipykernel install --user --name mswh
        jupyter labextension install jupyterlab-plotly

The examples are best explored using `JupyterLab`. Please check out the
[JupyterLab documentation](https://jupyterlab.readthedocs.io/en/latest/)
for further help as needed.

### Detailed Installation Steps

If for any reason a user encounters difficulties with the simple installation
instructions, the user is encouraged to consult a [more detailed installation guide that is
posted with the code documentation](https://lbnl-eta.github.io/MSWH/source/installation.html).

## Django Web Framework Deployment

### 1. Local

If the installation succeeded, to run the Django application navigate to the `web` folder (there should be a `manage.py` file) and start the development server on your local machine with:

        python manage.py runserver

   Now you can open your browser and type in `localhost:8000` (or `127.0.0.1:8000` if you are on a Windows machine) to start the web interface.

   Note that to build python extensions one needs to have `python3.x-dev` installed.

   Make sure that `DEBUG = True` in `settings.py`, this will ensure that the development server is able to serve local static files.

### 2. Public

#### Override settings locally

To deploy publicly, rename the file `local_settings_TEMPLATE.py` to `local_settings.py` and update the constants.

* `SECRET_KEY = '<random_string>'`

  The random string should be 50 characters long and can created (on Linux) by using the following command as super user:

        </dev/urandom tr -dc '1234567890!#$?*#-.,+qwertyuiopQWERTYUIOPasdfghjklASDFGHJKLzxcvbnmZXCVBNM' | head -c50; echo ""

    An example for a good secret key is this: `SECRET_KEY = 'O&2aYmv%)0B5#U-'9qsLTpfItC9N*V?%3L#fOHxDO,zyUm*S,U'`

* `DEBUG = False`

  Keep the `Debug` constant set to `True` in `settings.py` to get more debugging info for local deployement (using the Django development server). For the public deployement, you should set it to `False` (in `local_settings.py`).

#### Serving static files

> For detailed documentation on how to serve static files, see the official Django documentation:
> * [Managing static files](https://docs.djangoproject.com/en/3.1/howto/static-files/)
> * [Deploying static files](https://docs.djangoproject.com/en/3.1/howto/static-files/deployment/)

As the Django devlopment server is not meant for production and only serves static files if `Debug` is set to `True`, the static files used in the Django project need to be served another way.

At this point, two important aspects regarding how to deploy static files in production will be named:

1. Running this command, will create a folder `static` that will contain a copy of all static files from different directories across the Django project.
    ```
    python manage.py collectstatic
    ```
    > :warning: Run this command every time you update one of the static files in their respective location in the Django project folder.

2. Configure `nginx` to serve static files from the generated `static` folder by adding a `location /static` block to the server block of the `nginx` config file for the domain you serve the Django app with. This is an example  `nginx` server block:
    ```
    server {
      listen 80;
      server_name <domain>;

      access_log  /var/log/nginx/access.log;
      error_log  /var/log/nginx/error.log;

      location / {
        # For testing, using the django development server:
        # proxy_pass http://127.0.0.1:8000/;
        # For production, using gunicorn:
        proxy_pass http://unix:/run/swhweb.sock;
      }

      # Run 'python manage.py collectstatic' command in Django root project folder, so this folder will be created
      location /static {
        root <path_to_repo>/MSWH/web;
        try_files $uri $uri/ =404;
      }
    }
    ```
    Replace `<path_to_repo>` with the actual path to the MSWH repository and `<domain>` with your domain.

## Contributing

All are invited to contribute to the MSWH software through following the [Guidelines for Contributors](contributing.md).

### Automated tests

To run tests, from the `MSWH` folder use the following command modified according to the test module and method you intend to run:

    python -m unittest mswh.{my_module}.tests.{test_my_module}.{MyModuleTests}.{test_my_method}

## Publications

The code was used for the following publications:
* Coughlin, Katie, Milica Grahovac, Mohan Ganeshalingam, Robert Hosbach, and Vagelis Vossos. 2020. Costs and Benefits of Community versus Individual End-use Infrastructure for Solar Water Heating. California Energy Commission. CEC-XXX-2020-XXX. (in press)

* Grahovac, Milica, Katie Coughlin, Mohan Ganeshalingam, Robert Hosbach and Vagelis Vossos. 2020. Costs and Benefits of Community Scale Solar Water Heating. 2020 ACEEE Study on Energy Efficiency in Buildings. Pacific Grove, California. [Link to the paper with a video presentation](https://aceee2020.conferencespot.org/event-data/pdf/catalyst_activity_10923/catalyst_activity_paper_20200812133157248_498ce455_3a9c_4278_9088_6e3fdce5745b)

* Milica Grahovac, Katie Coughlin, Robert Hosbach, Hannes Gerhart, (2020). Multiscale Solar Water Heating. Journal of Open Source Software, 5(56), 2695, [![DOI](https://joss.theoj.org/papers/10.21105/joss.02695/status.svg)](https://doi.org/10.21105/joss.02695)

* Gerhart, H. (2019). Implementation of a Flexible Web Framework for Simulating Python System Models (p. 82). Technical University of Munich; Technical University of Munich. Research performed at LBNL. [Download link](https://gerhart.xyz/thesis.pdf)

* Web deployed version of the Django app is under construction [on this publicly available private website](https://solar.floweragenda.org/).

## About

The software may be distributed under the copyright and a BSD license provided in [legal.md](legal.md).

Milica Grahovac, Robert Hosbach, Katie Coughlin, Mohan Ganeshalingam and Hannes Gerhart created the contents of this repo
in the scope of the CEC "Costs and Benefits of Community vs. Individual End-Use Infrastructure for Solar Water Heating" project.

To cite use format provided at the [DOE CODE](https://www.osti.gov/doecode/biblio/26000) MSWH record.

## Acknowledgements

This work was supported by the California Energy Commission, Public Interest Energy Research Program, under Contract No. PIR-16-022.

We thank the reviewers and the editor of [The Journal of Open Source Software (JOSS)](https://joss.theoj.org/), [Bryn Pickering](https://github.com/brynpickering), [Nithiya Streethran](https://github.com/nmstreethran), and [Stefan Pfenninger](https://github.com/sjpfenninger) for their contributions in improving the code, the examples and the code documentation for the code release 2.0.0.
---
title: 'Multiscale Solar Water Heating'
tags:
  - Python
  - system
  - component
  - simulation
  - solar thermal water heating
  - solar electric water heating
  - heat pump water heating
  - natural gas water heating
  - photovoltaic
  - flat plate solar collector
  - evacuated tubes solar collector
  - hot water demand
  - solar radiation
  - thermal storage
  - solar water heating
  - sizing

authors:
  - name: Milica Grahovac^[corresponding author]
    affiliation: "1"
  - name: Katie Coughlin
    affiliation: "1"
  - name: Hannes Gerhart^[At the time of code creation was at affiliation 1 and 2]
    affiliation: "1, 2"
  - name: Robert Hosbach
    affiliation: "1"

affiliations:
  - name: Lawrence Berkeley National Laboratory, Berkeley, CA, USA
    index: 1
  - name: Technical University of Munich, Munich, Germany
    index: 2

date: September 7, 2020
bibliography: paper.bib
---


# Summary

The Multiscale Solar Water Heating (MSWH) package simulates individual and community scale solar water heating projects and allows for a comparison with the simulation performance of conventional natural gas tank water heaters (WH). The package contains a [Jupyter notebook with examples](https://github.com/LBNL-ETA/MSWH/blob/v2.0.0/scripts/MSWH%20System%20Tool.ipynb), a [graphical user interface (GUI) developed using Django Framework](https://github.com/LBNL-ETA/MSWH/tree/v2.0.0/web) and both functional and unit tests. System performance time series visualizations are available both in example notebooks and through the GUI, either spun off locally or [using a web deployed version](https://solar.floweragenda.org/).

The package was developed in the scope of a California Energy Commission (CEC) funded project looking at costs and benefits of using community versus individual scale solar thermal water heating systems. The database included in the MSWH software focuses primarily on California-specific hot water use profiles and climate data, but can structurally accommodate any further climate zones. The scale refers to the number of households served by a single system. Therefore, one can apply the models to explore the benefits of grouping multiple households to be served by a single solar water heating system in comparison to a system installed in a single household. Another example application of the models is to enable calculation of gas savings when switching from a gas WH to a solar WH in a single household.

The preconfigured system simulation models provided in the package include base-case gas tank WH and the following solar WH configurations with solar storage tanks:

* Solar thermal collector WH with either a tankless or a tank gas WH backup.
* Solar electric photovoltaic WH with a heat pump storage tank and an electric resistance backup.

[This documentation page](https://lbnl-eta.github.io/MSWH/source/models.html#approach-to-component-and-system-modeling-and-simulation) provides more details about the implemented models, the modeling approach, and the references used in some of the model development.

To evaluate a solar water heating project at the design phase by looking at its simulation performance the user should create a system instance for each compared system. This is described in detail in our [example notebooks](https://github.com/LBNL-ETA/MSWH/tree/v2.0.0/scripts). The user needs to specify the following:

* The project location by choosing one of the climate zones for which data is available in our database.
* For each household: count of people supplied by the system and whether there is any daytime household occupancy.

The MSWH software was used to perform the engineering analysis to estimate energy consumption and savings in the @Coughlin:2021 project report as well as in the @Grahovac:2020 research paper. @Gerhart:2019 published a master's thesis that explains the development of [the GUI](https://github.com/LBNL-ETA/MSWH/tree/v2.0.0/web) and how to effectively utilize it to facilitate easy addition of new MSWH software system models.

The MSWH software is both accessible and functionally robust. We performed extensive validation of [component models](https://github.com/LBNL-ETA/MSWH/blob/v2.0.0/mswh/system/tests/test_components.py) and [system models](https://github.com/LBNL-ETA/MSWH/blob/v2.0.0/mswh/system/tests/test_models.py) against performance results obtained using freely available open source tools and certification data generated using commercial tools. The validation models are a part of the test suite and show good agreement in all test comparisons.

The package is structured so that it can be extended with further technologies, applications, and locations as needed. The weather data are currently mostly limited to 16 California climate zones and can be extended to other climate zones. An example climate zone outside of California was added for Banja Luka, Bosnia and Herzegovina, through an [additional example Jupyter notebook](https://github.com/LBNL-ETA/MSWH/blob/v2.0.0/scripts/MSWH&#32;System&#32;Tool&#32;-&#32;Additional&#32;Climate.ipynb).

The code is available as DOE CODE, see @Doecode:2019.

# Statement of Need

A project that prompted the development of this software is described in @Coughlin:2021. The project enquired whether, on the state level in California, there exist any economic benefits from grouping households to be served by one community-level solar water heating installation, in comparison to having a single solar WH installation in each household.

Our primary motivation to develop new software was the combination of the following needs:

* The level of detail sufficient to allow for an investigation of transient effects of thermal storage.
* The simulation time for a single system short enough to allow for over 100,000 simulations to be performed on a personal computer within a reasonable amount of time.
* Simplicity of integration within the larger life-cycle cost framework as presented in @Coughlin:2021 and @Grahovac:2020.

We developed lightweight models that allow for around 120,000 hourly annual system simulation runs, including component auto-sizing and life-cycle cost analysis to be performed on a personal computer with an Intel(R) Core(TM) i7-8700 CPU @ 3.2GHz in about 8 hours. The users can expect an annual solar WH system simulation to complete in less than 0.2 seconds.

Policy developers and researchers could utilize the existing MSWH software by embedding it into a larger analysis framework for their specific research questions. Solar thermal water heating system planners, designers, and contractors may find it useful to have access to a freely available simulation tool that they can use to evaluate various system designs. Homeowners considering transitioning to a solar water heating system may be interested in analyzing a hypothetical system before seeking professional assistance. Further possible use cases are elaborated in [this section of the documentation](https://lbnl-eta.github.io/MSWH/source/models.html#future-applications-statement-of-need).

Simulation tools tend to be inaccessible to non-technical users. The MSWH software provides insight into what happens in a relatively simple simulation model due to the use of readable Python code, while the example notebooks and the GUI allow for instant utilization of the models. These features also make the code suitable for educators.

# Acknowledgements

This work was supported by the California Energy Commission, Public Interest Energy Research Program, under Contract No. PIR-16-022. We thank Vagelis Vossos and Mohan Ganeshalingam for their contributions and support.

We dedicate this paper to Dino Kosić.

# References
## What is in here?

This folder contains the Django source code for the web app.

File/Folder | Content
------ | ------
[data](data) | This folder contains the simulation results stored as pickled files.
[swhweb](swhweb) | This folder contains high-level, Django project related configuration files.
[system](system) | Here, the Django app data lives in.
[templates](templates) | This folder contains html files associated with general, project related content.
[db.sqlite3](db.sqlite3) | This is the Django project database.
[manage.py](manage.py) | This is the main Django utility file used to run the server, add a new Django application, etc.

The basic structure of the web app is split into the following sections available through the tabs in the navigation bar at the top of the page:

* **Home:** Overview of the functionality
* **Configurations:** A list of preconfigured systems. The user can edit, delete and create new configurations and place the system in a climate zone
* **Components:** A list of available system components, grouped in their respective class. The user can edit, delete and create new components
* **Projects:** Any system configuration can be deployed for a specific project. A project can be a single household or a community consisting of multiple households. On the _Project_ page, the user sees all available project, can edit or delete them as well as create new project configurations
* **Visualization:** After having invoked the simulation on the _Configurations_ page, this page provides a list with available simulation results for a respective system configuration. When clicking on an item, the interactive plot as well as a table with annual totals can be accessed for further analyses of the system's performance
* **Admin:** Django provides an administration interface, which can be reached through this tab. A super user has been created as follows:

      User name:  admin
      Password:   $SWH2019
Multiscale Solar Water Heating System Modeling
==============================================



.. toctree::
   :maxdepth: 3

   source/installation
   source/models
   source/modules
   source/acknowledgements
   source/bibliography

   source/legal

.. role:: note
Tools
=====


mswh.tools.plots module
-----------------------

.. automodule:: mswh.tools.plots
    :members:
    :undoc-members:
    :show-inheritance:

mswh.tools.unit\_converters module
----------------------------------

.. automodule:: mswh.tools.unit_converters
    :members:
    :undoc-members:
    :show-inheritance:
References
==========

.. bibliography:: references.bib
   :all:
   :style: alphaSystem and Component Models
===========================


mswh.system.components module
------------------------------

.. automodule:: mswh.system.components
    :members:
    :undoc-members:
    :show-inheritance:


mswh.system.models module
--------------------------

.. automodule:: mswh.system.models
    :members:
    :undoc-members:
    :show-inheritance:

mswh.system.source\_and\_sink module
--------------------------------------

.. automodule:: mswh.system.source_and_sink
    :members:
    :undoc-members:
    :show-inheritance:
Copyright Notice
================

Multiscale Solar Water Heating (MSWH) Copyright (c) 2019, The
Regents of the University of California, through Lawrence Berkeley National
Laboratory (subject to receipt of any required approvals from the U.S.
Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative
works, and perform publicly and display publicly, and to permit other to do
so.

License Agreement
=================

Multiscale Solar Water Heating (MSWH) Copyright (c) 2019, The
Regents of the University of California, through Lawrence Berkeley National
Laboratory (subject to receipt of any required approvals from the U.S.
Dept. of Energy).  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches,
or upgrades to the features, functionality or performance of the source
code ("Enhancements") to anyone; however, if you choose to make your
Enhancements available either publicly, or directly to Lawrence Berkeley
National Laboratory, without imposing a separate written license agreement
for such Enhancements, then you hereby grant the following license: a
non-exclusive, royalty-free perpetual license to install, use, modify,
prepare derivative works, incorporate into other computer software,
distribute, and sublicense such enhancements or derivative works thereof,
in binary and source code form.
Package Installation
====================

Most users will be able to install the MSWH Python package by following the 
`Setup and Installation and Simple Installation Using Conda <https://github.com/LBNL-ETA/MSWH#setup-and-installation>`_ 
section of the README.md file that is displayed at the landing page of the MSWH repository. Please use those instructions
as the primary approach to MSWH package installation.

The set of instructions presented here is intended for technical users that are relatively new to virtual environments 
or Python in general, for users who encountered issues with the 
simple installation instructions available in the `README.md <https://github.com/LBNL-ETA/MSWH>`_ file, or any other 
users looking for a reminder on some of the installation steps. These instructions also show to the users how to 
utilize an alternative Python package management system, `venv <https://docs.python.org/3.8/library/venv.html>`_.

Please make sure to install ``pip`` and, in case you are not using ``venv``, ``conda`` as instructed 
in `Setup and Installation on the readme file <https://github.com/LBNL-ETA/MSWH#setup-and-installation>`_.

Here are the detailed steps to install the `MSWH Python package <https://github.com/LBNL-ETA/MSWH>`_:

#. Since the repo comes with database files, please download, install and see the documentation for `git large file storage <https://git-lfs.github.com/>`_.

#. It is recommended to create a new ``Python`` environment in order to avoid interference with the system-wide Python installation, 
   for example by using `conda <(https://docs.conda.io/en/latest/>`_ or `venv <https://docs.python.org/3.8/library/venv.html>`_. 
   Depending on the approach you take, pick one of the commands below and run it in a terminal to create a new environment named, for instance, ``mswh``.

    If you use ``conda`` from the repo clone folder run:

    .. code-block:: console

        conda create -n mswh python=3.8

    If you use ``venv``, for example on ``Linux``:

    .. code-block:: console

        python3.8 -m venv <path_to_env>/mswh

    With ``<path_to_env>`` as your selected folder path to store virtual
    environments.

#. Now the virtual environment needs to be activated, by running one of the following commands:

    When using ``Anaconda`` or ``Miniconda``:

    .. code-block:: console

        conda activate mswh

    When using ``venv``:

    .. code-block:: console

        source <path_to_env>/mswh/bin/activate

    After having activated the virtual environment, the name of it should appear before the prompt in the terminal.

    For deactivating use:
    
    .. code-block:: console
    
        conda deactivate

#. To make use of example ``Jupyter notebooks`` one should have `JupyterLab <https://jupyter.org/install>`_ installed. 
   To ensure the same Python kernel can be used in a ``Jupyter notebook``, activate the virtual environment and run:

    .. code-block:: console
        
        python -m ipykernel install --user --name mswh

   Users with admin privileges can skip the ``--user`` flag.

   If you have any issues with plots not being displayed when running the example notebooks,
   please install the following:

    .. code-block:: console

        jupyter labextension install jupyterlab-plotly

#. Clone the repository with:
    
    .. code-block:: console
        
        git clone https://github.com/LBNL-ETA/MSWH.git

#.  To install the necessary Python packages navigate to the ``setup.py`` directory and run:

    .. code-block:: console
        
        pip install -e .

    The ``-e`` flag is only necessary if one would like changes to the source code be reflected immediately 
    (without having to rerun the ``setup.py`` script with every change to the source code). 
    If you just want to run the project application, you can omit the ``-e`` flag.

#. To use the plotting capabilities, also required when running tests, please install `orca <https://github.com/plotly/orca>`_.
Multiscale Solar Water Heating (MSWH)
=====================================

Scope
^^^^^

The main purpose of the Multiscale Solar Water Heating (MSWH) software is to model energy use for individual and community scale solar water heating projects in California.

The package contains functional and unit tests and it is structured so that it can be extended with further technologies, applications, and locations.

Usage
^^^^^

The user provides a climate zone for a project, occupancy for each household, and whether any of the occupants stay at home during the day. The software can then load a set of example California specific hourly domestic hot water end-use load profiles from a database, size, and locate the systems. Next, the user can simulate the hourly system performance over a period of one representative year, visualize and explore the simulation results using time-series plots for temperature profiles, heat and power rates, or look at annual summaries. Similarly, the user can model individual household solar water heating projects and base case conventional gas tank water heater systems, such that the results can be compared between the individual, community-scale, and base case systems.

This functionality is readily available through a `Jupyter notebook <https://github.com/LBNL-ETA/MSWH/blob/v2.0.0/scripts/MSWH%20System%20Tool.ipynb>`_ and a `Django web framework graphical user interface (GUI) <https://github.com/LBNL-ETA/MSWH/tree/v2.0.0/web>`_, depending on what level of detail the user would like to access. Please see the README file on the `MSWH repo <https://github.com/LBNL-ETA/MSWH>`_ for further usage and installation instructions.

System performance time series visualizations are available both in example notebooks and through the GUI, either spun off locally or `using a web deployed version <https://solar.floweragenda.org/>`_.

Features
^^^^^^^^

This software package contains the following Python modules:

* Solar irradiation on a tilted surface.

* Simplified component models for:

    * Converters: solar collectors, electric resistance heater, gas burner, photovoltaic panels, heat pump.
    * Storage: solar thermal tank, heat pump thermal tank, conventional gas tank water heater.
    * Distribution: distribution and solar pump, piping losses.

* Preconfigured system simulation models for: 

    * Base case gas tank water heaters.
    * Solar thermal water heaters (solar collector feeding a storage tank, with a tankless gas water heater backup in a new installation cases and a basecase gas tank water heater in a retrofit case). 
    * Solar electric water heaters (grid supported photovoltaic panel powering a heat pump storage tank with an electric resistance backup).

* Database with component performance parameters, California specific weather data, and domestic hot water end-use load profiles.

* Django web framework to configure project, parametrize components and `simulate from a web browser <https://solar.floweragenda.org/>`_.

We also developed component sizing rules and size scaling rules to account for the household occupancy and project scale, respectively. The rules are readily available in the example notebooks and can easily be modified for exploratory purposes that we further describe in the Statement of Need section. For the sizing and scaling rules, we used the following data sources: expert knowledge, web-scraped data with the help of a tool described in :cite:`Gerke:2017`, sizing rules available in :cite:`Csi_thermal:2016`, and certification databases such as :cite:`Ccms:2018` and :cite:`Cec_appliance:2019`.

Approach to Component and System Modeling and Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this section we briefly introduce the characteristics of the underlying models and simulation. 

We performed an extensive literature rewiew prior to developing the models. Modelica buildings library by :cite:`Wetter:2014` exceeds the level of detail but proves too detailed and thus somewhat slow for our particular application. SAM tool (:cite:`Blair:2014`) has a fitting level of detail, provides most of the system models that we needed but for our purposes proves not flexible enough in terms of modifying the system configuration, automating the size scaling, and embedding it into our custom life-cycle cost framework.

Namely, to capture a sufficient level of detail of the California demographics, such as variability in climate zones, household types, and household occupancy, we wanted to be able to simulate a few alternative water heating systems in each of the California sample households. Secondly, to get a more realistic picture of the effect of thermal storage and distribution system losses, we opted to perform a simulation with relatively short time-steps of one hour for a duration of one representative year. We were not able to identify an open source tool that is capable of firstly satisfying the simulation speed requirement combined with the necessary level of detail for our analysis and secondly providing the flexibility for us to customize various integral parts of the analysis such as automate the component and system size scaling, specify hot water load profiles and solar radiation for each household or group of households in the sample.

To satisfy our research need we thus opted to develop lightweight simulation models for all involved systems that would allow for around 120,000 simulation runs together with the component sizing and life-cycle cost analysis to be performed on a computer with a 12-core processor in about 8 hours. The users can expect a single solar water heater simulation model to run in less than one second (the developers were experiencing run times on the order of 0.2 seconds), providing an almost instantaneous experience for a user only seeking to design and investigate a single system.

We developed and implemented simplified fast performing energy balance based component models. We connected the component models into two preconfigured solar water heating systems, that are both provided with the MSWH software. Those models are:

* Solar thermal collector, hot water thermal storage tank, with a selection of backups: gas storage water heater or an instantaneous gas water heater.
* Photovoltaic panel, heat pump tank water heater, with an electric resistance water heater as backup.

We built a simple simulation solver that uses the explicit forward Euler method to solve the balance equations in each simulation time-step.

The component models were either developed from scratch or implemented in Python based on existing models identified in the literature. We implemented the following existing or new models:

* Solar irradiation on a tilted surface model is based on equations found in :cite:`Duffie:2013`.
* Solar collector models and model parameters are based on :cite:`Ashrae:2013` and :cite:`Srcc:2013`.
* We converted the natural gas tank water heater model from :cite:`Lutz:1998` into an hourly time-step model implementation.
* Photovoltaic model is based on a simplified model found in :cite:`Wetter:2014`.
* Heat pump water heater tank is based on :cite:`Sparn:2014`.
* Solar thermal tank is a phenomenological model based on ideas very similar to the model developed for NREL's SAM software (:cite:`Blair:2014`), as described in :cite:`DiOrio:2014`.
* Simplified performance data-based gas burner model was implemented to represent instantaneous gas water heater.
* Simple electric resistance model was implemented to represent instantaneous electric water heater.
* We developed a simplified data based solar and distribution pump model.
* To model the distribution piping network we developed a simplified model that is capable of accounting for thermal losses at stagnation and flows on-demand with correction factors available to help account for the relatively long time-step of one hour.

More details on the hot water demand model used in creating the database of sample hot water use load profiles, as well as extensive detail on the software's solar radiation, component and system models can be found in the project report by :cite:`Coughlin:2021`. :cite:`Gerhart:2019` thesis provides additional details on the solar electric system model development.

Note that the weather data are currently mostly limited to California and can be extended to other climate zones. An example climate zone outside of California was added for Banja Luka, Bosnia and Herzegovina, through an `additional example Jupyter notebook <https://github.com/LBNL-ETA/MSWH/blob/v2.0.0/scripts/MSWH%20System%20Tool%20-%20Additional%20Climate.ipynb>`_. The water consumption profiles can be highly location specific and their development for additional climate zones would require new research efforts. A quick approximation may be made with caution by scaling the California profiles to match the location-specific estimate of the average annual water use. This is possible as the shape of each daily profile can be assumed similar and sufficiently variable to allow for the study of transient and peak load effects at any location. The weather processor is TMY3 enabled and the user may populate the database with additional climates as needed.

The energy sources we consider are solar irradiation, gas, and electricity. The source energy is converted, if needed stored, and distributed to meet the end-use loads for each household.

Upon assembling the components into systems, we perform an annual simulation with hourly timesteps. We solve any differential equations for each time step using an explicit forward Euler method, a first order technique that provides a good approximation given the dynamics of the process observed and the level of detail required in our analysis.

We configure and size each MSWH thermal configuration so that it complies with the CSI-T (California Solar Initiative - Thermal) rebate program sizing requirements. The system model assumes appropriate flow and temperature controls and includes freeze and stagnation protection.

Future Applications - Statement of Need
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When it comes to the future application of the MSWH software, we can envision four main groups of users:

* Researchers and policy developers.
* Solar water heating planners, designers, and contractors.
* Homeowners.
* Educators.

If the features of the existing MSWH software are sufficient for their application, the policy developers and researchers could utilize the existing MSWH software by embedding it into some larger analysis framework they construct such that it provides answers to their specific research questions. Should they require additional system configurations and even additional components, the existing framework should be expanded in line with the structure made available to the user in the MSWH software. When systems are added following the structure of the existing systems, the addition of such a new system to the GUI is made possible by using the flexible web framework.

Solar thermal water heating system planners, designers, and contractors may find it useful to have access to a freely available simulation tool, such as the MSWH software, that they can use to evaluate various system designs. The design parameters that such users can easily modify are household occupancies, climate zone, collector and tank sizes, component performance parameters such as insulation level of any thermal storage tanks, and types of solar collectors. The MSWH software relies on standard collector rating data readily available for most designs found on the market today. For each proposed design the MSWH software will output, among other results, the solar fraction and the backup energy use on an annual level, the two variables allowing for a quick cross-comparison for the proposed designs.

Similarly, homeowners considering transitioning to a solar water heating system may be interested in analyzing a hypothetical system before seeking further professional help. Or, some homeowners may simply be interested in learning about both solar water heating systems and system simulation in general. Another example use case would be to enable the occupants of households that:

* Are retrofitting an existing system due to an increase or decrease in occupancy, or
* Already possess one of the components and are looking to appropriately size the others

to simulate alternatives and compare the obtained energy consumption and solar fraction results for any alternative designs they like to define.

Lastly, simulation tools tend to be inaccessible to non-technical users, both in terms of usage and the chance for the user to understand the underlying codebase just by reading through it. The MSWH software provides a unique insight into what actually happens in a relatively simple mezzo-level simulation model due to the use of readable Python code, while the example notebooks and GUI allow for instant utilization of the models. These features make the code suitable also for educators.

Code Development and Code Contributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We welcome code contributions. The development primarily takes place on the `MSWH GitHub repository <https://github.com/LBNL-ETA/MSWH>`_. Please refer to the `contributing guidelines <https://github.com/LBNL-ETA/MSWH/blob/master/contributing.md>`_ and `README.md <https://github.com/LBNL-ETA/MSWH/blob/master/README.md>`_ for further instructions, including those on running the unit tests.Acknowledgements
================

This work was supported by the California Energy Commission, Public Interest Energy Research Program, under Contract No. PIR-16-022.Python Code Documentation
=========================

.. toctree::
   :maxdepth: 4

Subpackages
-----------

.. toctree::

    mswh.system
    mswh.tools
    mswh.comm
Database Communication
======================


mswh.comm.sql module
--------------------

.. automodule:: mswh.comm.sql
    :members:
    :undoc-members:
    :show-inheritance:
