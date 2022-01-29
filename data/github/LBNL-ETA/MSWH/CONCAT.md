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
