---
title: 'MCALF: Multi-Component Atmospheric Line Fitting'
tags:
  - Python
  - astronomy
  - solar physics
  - spectrum
  - spectra
  - fitting
  - absorption
  - emission
  - voigt
authors:
  - name: Conor D. MacBride
    orcid: 0000-0002-9901-8723
    affiliation: 1
  - name: David B. Jess
    orcid: 0000-0002-9155-8039
    affiliation: "1, 2"
affiliations:
 - name: Astrophysics Research Centre, School of Mathematics and Physics, Queen's University Belfast, Belfast, BT7 1NN, UK
   index: 1
 - name: Department of Physics and Astronomy, California State University Northridge, Northridge, CA 91330, U.S.A.
   index: 2
date: 19 April 2021
bibliography: paper.bib
---

# Summary

Determining accurate velocity measurements from observations of the Sun is of vital importance to solar physicists who are studying the wave dynamics in the solar atmosphere. Weak chromospheric absorption lines, due to dynamic events in the solar atmosphere, often consist of multiple spectral components. Isolating these components allows for the velocity field of the dynamic and quiescent regimes to be studied independently. However, isolating such components is particularly challenging due to the wide variety of spectral shapes present in the same dataset. `MCALF` provides a novel method and infrastructure to determine Doppler velocities in a large dataset. Each spectrum is fitted with a model adapted to its specific spectral shape.

# Statement of need

MCALF is an open-source Python package for accurately constraining velocity information from spectral imaging observations using machine learning techniques. This software package is intended to be used by solar physicists trying to extract line-of-sight (LOS) Doppler velocity information from spectral imaging observations (Stokes $I$ measurements) of the Sun. This `toolkit' can be used to define a spectral model optimised for a particular dataset.

This package is particularly suited for extracting velocity information from spectral imaging observations where the individual spectra can contain multiple spectral components. Such multiple components are typically present when active solar phenomena occur within an isolated region of the solar disk. Spectra within such a region will often have a large emission component superimposed on top of the underlying absorption spectral profile from the quiescent solar atmosphere [@Felipe:2014]. Being able to extract velocity information from such observations would provide solar physicists with a wider range of data products that can be used for science [@Stangalini:2020]. This package implements the novel approach of automated classification of spectral profiles prior to fitting a model.

A sample model is provided for an IBIS Ca $\text{\sc{ii}}$ 8542 Å spectral imaging sunspot dataset. This dataset typically contains spectra with multiple atmospheric components and this package supports the isolation of the individual components such that velocity information can be constrained for each component. The method implemented in this IBIS model has been discussed extensively in @MacBride:2020. There are also several ongoing research projects using this model to extract velocity measurements.

Using this sample model, as well as the separate base (template) model it is built upon, a custom model can easily be built for a specific dataset. The custom model can be designed to take into account the spectral shape of each particular spectrum in the dataset. By training a neural network classifier using a sample of spectra from the dataset labelled with their spectral shapes, the spectral shape of any spectrum in the dataset can be found. The fitting algorithm can then be adjusted for each spectrum based on the particular spectral shape the neural network assigned it. The `toolkit' nature of this package also allows the possibility of utilising existing machine learning classifiers, such as the ``supervised hierarchical $k$-means" classifier introduced in @Panos:2018, which clusters solar flare spectra based on their profile shape.

This package is designed to run in parallel over large data cubes, as well as in serial. As each spectrum is processed in isolation, this package scales very well across many processor cores. Numerous functions are provided to plot the results clearly, some of which are showcased in \autoref{fig:example}. The `MCALF` API also contains many useful functions which have the potential of being integrated into other Python packages. Full documentation as well as examples on how to use `MCALF` are provided at [mcalf.macbride.me](https://mcalf.macbride.me).

![An overview of some of the plotting functions that are included in `MCALF`.\label{fig:example}](figure.pdf)

# Acknowledgements

CDM would like to thank the Northern Ireland Department for the Economy for the award of a PhD studentship. DBJ wishes to thank Invest NI and Randox Laboratories Ltd. for the award of a Research and Development Grant (059RDEN-1) that allowed the computational techniques employed to be developed. DBJ would also like to thank the UK Science and Technology Facilities Council (STFC) for the consolidated grant ST/T00021X/1. The authors wish to acknowledge scientific discussions with the Waves in the Lower Solar Atmosphere (WaLSA; [www.WaLSA.team](https://www.WaLSA.team)) team, which is supported by the Research Council of Norway (project no. 262622) and the Royal Society (award no. Hooke18b/SCTM).

# References


Contributor Covenant Code of Conduct
====================================

Our Pledge
----------

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

Our Standards
-------------

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

Enforcement Responsibilities
----------------------------

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

Scope
-----

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

Enforcement
-----------

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
`mcalf@macbride.me <mailto:mcalf@macbride.me>`_.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

Enforcement Guidelines
----------------------

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

1. Correction
~~~~~~~~~~~~~

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

2. Warning
~~~~~~~~~~

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

3. Temporary Ban
~~~~~~~~~~~~~~~~

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

4. Permanent Ban
~~~~~~~~~~~~~~~~

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

Attribution
-----------

This Code of Conduct is adapted from the
`Contributor Covenant <https://www.contributor-covenant.org>`_,
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by `Mozilla's code of conduct
enforcement ladder <https://github.com/mozilla/diversity>`_.

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
===============================================
MCALF: Multi-Component Atmospheric Line Fitting
===============================================

|Azure Pipelines Status| |Codecov| |PyPI Version| |Zenodo DOI| |Docs Status| |GitHub License|

MCALF is an open-source Python package for accurately constraining velocity
information from spectral imaging observations using machine learning
techniques.

This software package is intended to be used by solar physicists trying
to extract line-of-sight (LOS) Doppler velocity information from
spectral imaging observations (Stokes I measurements) of the Sun.
A ‘toolkit’ is provided that can be used to define a spectral model
optimised for a particular dataset.

This package is particularly suited for extracting velocity information
from spectral imaging observations where the individual spectra can
contain multiple spectral components.
Such multiple components are typically present when active solar phenomenon
occur within an isolated region of the solar disk.
Spectra within such a region will often have a large emission component
superimposed on top of the underlying absorption spectral profile from the
quiescent solar atmosphere.

A sample model is provided for an IBIS Ca II 8542 Å spectral imaging sunspot
dataset.
This dataset typically contains spectra with multiple atmospheric
components and this package supports the isolation of the individual
components such that velocity information can be constrained for each
component.
Using this sample model, as well as the separate base (template) model it is
built upon, a custom model can easily be built for a specific dataset.

The custom model can be designed to take into account the spectral shape of
each particular spectrum in the dataset.
By training a neural network classifier using a sample of spectra from the
dataset labelled with their spectral shapes, the spectral shape of any
spectrum in the dataset can be found.
The fitting algorithm can then be adjusted for each spectrum based on
the particular spectral shape the neural network assigned it.

This package is designed to run in parallel over large data cubes, as well
as in serial.
As each spectrum is processed in isolation, this package scales very well
across many processor cores.
Numerous functions are provided to plot the results in a clearly.
The MCALF API also contains many useful functions which have the potential
of being integrated into other Python packages.

Installation
------------

For easier package management we recommend using `Miniconda`_ (or `Anaconda`_)
and creating a `new conda environment`_ to install MCALF inside.
To install MCALF using `Miniconda`_, run the following commands in your
system's command prompt, or if you are using Windows, in the
'Anaconda Prompt':

.. code:: bash

    $ conda config --add channels conda-forge
    $ conda config --set channel_priority strict
    $ conda install mcalf

MCALF is updated to the latest version by running:

.. code:: bash

    $ conda update mcalf

Alternatively, you can install MCALF using ``pip``:

.. code:: bash

    $ pip install mcalf

Testing
-------

A test suite is included with the package. The package is tested on
multiple platforms, however you may wish to run the tests on your
system also. More details on running our tox/pytest test suite are
available in our `documentation`_.

Getting Started
---------------

Documentation is `available here <https://mcalf.macbride.me/>`_.
Some examples are included `here <examples/>`_.
If you are interested in using this package in your research and would
like advice on how to use this package, please contact `Conor MacBride`_.

Contributing
------------

|Contributor Covenant|

If you find this package useful and have time to make it even better,
you are very welcome to contribute to this package, regardless of how much
prior experience you have.
Types of ways you can contribute include, expanding the documentation with
more use cases and examples, reporting bugs through the GitHub issue tracker,
reviewing pull requests and the existing code, fixing bugs and implementing new
features in the code.

You are encouraged to submit any `bug reports`_ and `pull requests`_ directly
to the `GitHub repository`_.
If you have any questions regarding contributing to this package please
contact `Conor MacBride`_.

Please note that this project is released with a Contributor Code of Conduct.
By participating in this project you agree to abide by its terms.

Citation
--------

If you have used this package in work that leads to a publication, we would
be very grateful if you could acknowledge your use of this package in the
main text of the publication.
Please cite the following publications,

    MacBride CD, Jess DB. 2021
    MCALF: Multi-Component Atmospheric Line Fitting.
    *Journal of Open Source Software*. **6(61)**, 3265.
    (`doi:10.21105/joss.03265 <https://doi.org/10.21105/joss.03265>`_)

..

    MacBride CD, Jess DB, Grant SDT, Khomenko E, Keys PH, Stangalini M. 2020
    Accurately constraining velocity information from spectral imaging
    observations using machine learning techniques.
    *Philosophical Transactions of the Royal Society A*. **379**, 2190.
    (`doi:10.1098/rsta.2020.0171 <https://doi.org/10.1098/rsta.2020.0171>`_)

Please also cite the `Zenodo DOI`_ for the package version you used.
Please also consider integrating your code and examples into the package.

License
-------

MCALF is licensed under the terms of the BSD 2-Clause license.

.. |Azure Pipelines Status| image:: https://dev.azure.com/ConorMacBride/mcalf/_apis/build/status/ConorMacBride.mcalf?repoName=ConorMacBride%2Fmcalf&branchName=main
    :target: https://dev.azure.com/ConorMacBride/mcalf/_build/latest?definitionId=5&repoName=ConorMacBride%2Fmcalf&branchName=main
    :alt: Azure Pipelines
.. |Codecov| image:: https://codecov.io/gh/ConorMacBride/mcalf/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/ConorMacBride/mcalf
    :alt: Codecov
.. |PyPI Version| image:: https://img.shields.io/pypi/v/mcalf
    :target: https://pypi.python.org/pypi/mcalf
    :alt: PyPI
.. |Zenodo DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3924527.svg
    :target: https://doi.org/10.5281/zenodo.3924527
    :alt: DOI
.. |Docs Status| image:: https://readthedocs.org/projects/mcalf/badge/?version=latest&style=flat
    :target: https://mcalf.macbride.me/
    :alt: Documentation
.. |GitHub License| image:: https://img.shields.io/github/license/ConorMacBride/mcalf
    :target: LICENSE.rst
    :alt: License
.. |Contributor Covenant| image:: https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg
    :target: CODE_OF_CONDUCT.rst
    :alt: Code of Conduct

.. _Anaconda: https://www.anaconda.com/products/individual#Downloads
.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _new conda environment: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
.. _documentation: https://mcalf.macbride.me/en/latest/guide/index.html#testing

.. _Conor MacBride: https://macbride.me/

.. _bug reports: https://github.com/ConorMacBride/mcalf/issues
.. _pull requests: https://github.com/ConorMacBride/mcalf/pulls
.. _GitHub repository: https://github.com/ConorMacBride/mcalf

.. _Zenodo DOI: https://doi.org/10.5281/zenodo.3924527
Example Gallery
***************

Here are a collection of examples on how this package can be used.
Models
======

Below are examples of how to use the models included within the models module:
Visualisation
=============

Below are examples of plots produced using functions within the visualisation module:
====================================
example1: Basic usage of the package
====================================

``FittingIBIS.ipynb``
---------------------
This file is an IPython Notebook containing examples of how to use the package
to accomplish typical tasks.

``FittingIBIS.pro``
-------------------
This file is similar to ``FittingIBIS.ipynb`` file, except it written is IDL.
It is not recommended to use the IDL wrapper in production, just use it to
explore the code if you are familiar with IDL and not Python.
If you wish to use this package, please use the Python implementation.
IDL is not fully supported in the current version of the code for reasons
such as, the Python tuple datatype cannot be passed from IDL to Python,
resulting in certain function calls not being possible.

``config.yml``
--------------
This is an example configuration file containing default parameters.
This can be easier than setting the parameters in the code.
The file follows the YAML_ format.

.. _YAML: https://pyyaml.org/wiki/PyYAMLDocumentation
MCALF: Multi-Component Atmospheric Line Fitting
***********************************************

MCALF Documentation
===================

Welcome to MCALF's documentation!

MCALF is an open-source Python package for accurately constraining
velocity information from spectral imaging observations using
machine learning techniques.

These pages document how the package can be interacted with.
Some examples are also provided.
A :any:`Documentation Index <genindex>` and
a :any:`Module Index <modindex>` are available.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   guide/index
   gallery/index
   code_ref/index
   code_of_conduct
   license

Contributor Covenant Code of Conduct
====================================

Our Pledge
----------

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

Our Standards
-------------

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

Enforcement Responsibilities
----------------------------

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

Scope
-----

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

Enforcement
-----------

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
`mcalf@macbride.me <mailto:mcalf@macbride.me>`_.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

Enforcement Guidelines
----------------------

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

1. Correction
~~~~~~~~~~~~~

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

2. Warning
~~~~~~~~~~

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

3. Temporary Ban
~~~~~~~~~~~~~~~~

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

4. Permanent Ban
~~~~~~~~~~~~~~~~

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

Attribution
-----------

This Code of Conduct is adapted from the
`Contributor Covenant <https://www.contributor-covenant.org>`_,
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by `Mozilla's code of conduct
enforcement ladder <https://github.com/mozilla/diversity>`_.

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
==================
User Documentation
==================

|Azure Pipelines Status| |Codecov| |PyPI Version| |Zenodo DOI| |Docs Status| |GitHub License|

MCALF is an open-source Python package for accurately constraining velocity
information from spectral imaging observations using machine learning
techniques.

This software package is intended to be used by solar physicists trying
to extract line-of-sight (LOS) Doppler velocity information from
spectral imaging observations (Stokes I measurements) of the Sun.
A ‘toolkit’ is provided that can be used to define a spectral model
optimised for a particular dataset.

This package is particularly suited for extracting velocity information
from spectral imaging observations where the individual spectra can
contain multiple spectral components.
Such multiple components are typically present when active solar phenomenon
occur within an isolated region of the solar disk.
Spectra within such a region will often have a large emission component
superimposed on top of the underlying absorption spectral profile from the
quiescent solar atmosphere.

A sample model is provided for an IBIS Ca II 8542 Å spectral imaging sunspot
dataset.
This dataset typically contains spectra with multiple atmospheric
components and this package supports the isolation of the individual
components such that velocity information can be constrained for each
component.
Using this sample model, as well as the separate base (template) model it is
built upon, a custom model can easily be built for a specific dataset.

The custom model can be designed to take into account the spectral shape of
each particular spectrum in the dataset.
By training a neural network classifier using a sample of spectra from the
dataset labelled with their spectral shapes, the spectral shape of any
spectrum in the dataset can be found.
The fitting algorithm can then be adjusted for each spectrum based on
the particular spectral shape the neural network assigned it.

This package is designed to run in parallel over large data cubes, as well
as in serial.
As each spectrum is processed in isolation, this package scales very well
across many processor cores.
Numerous functions are provided to plot the results in a clearly.
The MCALF API also contains many useful functions which have the potential
of being integrated into other Python packages.

Installation
------------

For easier package management we recommend using `Miniconda`_ (or `Anaconda`_)
and creating a `new conda environment`_ to install MCALF inside.
To install MCALF using `Miniconda`_, run the following commands in your
system's command prompt, or if you are using Windows, in the
'Anaconda Prompt':

.. code:: bash

    $ conda config --add channels conda-forge
    $ conda config --set channel_priority strict
    $ conda install mcalf

MCALF is updated to the latest version by running:

.. code:: bash

    $ conda update mcalf

Alternatively, you can install MCALF using ``pip``:

.. code:: bash

    $ pip install mcalf

Testing
-------

A test suite is included with the package. The package is tested on
multiple platforms, however you may wish to run the tests on your
system also.

Installing Dependencies
~~~~~~~~~~~~~~~~~~~~~~~

Using MCALF with pip
====================

To run the tests you need a number of extra packages installed. If you
installed MCALF using pip, you can run ``pip install mcalf[tests]`` to
install the additional testing dependencies (and MCALF if it's not
already installed).

Using MCALF with conda
======================

If you want to use MCALF inside a conda environment you should first
follow the conda installation instructions above. Once MCALF is
installed in a conda environment, ask conda to install each of MCALF's
testing dependencies using the following command.
(See `setup.cfg`_ for an up-to-date list of dependencies.)

.. code:: bash

    $ conda install pytest pytest-cov tox

Running Tests
~~~~~~~~~~~~~

Tests should be run within the virtual environment where MCALF and its
testing dependencies were installed. Run the following command to test
your installation,

.. code:: bash

    $ pytest --pyargs mcalf

Editing the Code
~~~~~~~~~~~~~~~~

If you are planning on making changes to your local version of the code,
it is recommended to run the test suite to help ensure that the changes
do not introduce problems elsewhere.

Before making changes, you'll need to set up a development environment.
The SunPy Community have compiled an excellent set of instructions and
is available in their `documentation`_. You can mostly replace
``sunpy`` with ``mcalf``, and install with

.. code:: bash

    $ pip install -e .[tests,docs]

After making changes to the MCALF source, run the MCALF test suite with
the following command (while in the same directory as ``setup.py``),

.. code:: bash

    $ pytest --pyargs mcalf --cov

The tox package has also been configured to run the MCALF test suite.

Getting Started
---------------

The following examples provide the key details on how to use this package.
For more details on how to use the particular classes and function,
please consult the `Code Reference <../code_ref/index.html>`_.
We plan to expand this section with more examples of this package being
used.

.. toctree::

   ../gallery/index
   examples/example1/index
   examples/neural_network/LabellingTutorial

If you are interested in using this package in your research and would
like advice on how to use this package, please contact `Conor MacBride`_.

Contributing
------------

|Contributor Covenant|

If you find this package useful and have time to make it even better,
you are very welcome to contribute to this package, regardless of how much
prior experience you have.
Types of ways you can contribute include, expanding the documentation with
more use cases and examples, reporting bugs through the GitHub issue tracker,
reviewing pull requests and the existing code, fixing bugs and implementing new
features in the code.

You are encouraged to submit any `bug reports`_ and `pull requests`_ directly
to the `GitHub repository`_.
If you have any questions regarding contributing to this package please
contact `Conor MacBride`_.

Please note that this project is released with a Contributor Code of Conduct.
By participating in this project you agree to abide by its terms.

Citation
--------

If you have used this package in work that leads to a publication, we would
be very grateful if you could acknowledge your use of this package in the
main text of the publication.
Please cite the following publications,

    MacBride CD, Jess DB. 2021
    MCALF: Multi-Component Atmospheric Line Fitting.
    *Journal of Open Source Software*. **6(61)**, 3265.
    (`doi:10.21105/joss.03265 <https://doi.org/10.21105/joss.03265>`_)

..

    MacBride CD, Jess DB, Grant SDT, Khomenko E, Keys PH, Stangalini M. 2020
    Accurately constraining velocity information from spectral imaging
    observations using machine learning techniques.
    *Philosophical Transactions of the Royal Society A*. **379**, 2190.
    (`doi:10.1098/rsta.2020.0171 <https://doi.org/10.1098/rsta.2020.0171>`_)

Please also cite the `Zenodo DOI`_ for the package version you used.
Please also consider integrating your code and examples into the package.

License
-------

MCALF is licensed under the terms of the BSD 2-Clause license.

.. |Azure Pipelines Status| image:: https://dev.azure.com/ConorMacBride/mcalf/_apis/build/status/ConorMacBride.mcalf?repoName=ConorMacBride%2Fmcalf&branchName=main
    :target: https://dev.azure.com/ConorMacBride/mcalf/_build/latest?definitionId=5&repoName=ConorMacBride%2Fmcalf&branchName=main
    :alt: Azure Pipelines
.. |Codecov| image:: https://codecov.io/gh/ConorMacBride/mcalf/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/ConorMacBride/mcalf
    :alt: Codecov
.. |PyPI Version| image:: https://img.shields.io/pypi/v/mcalf
    :target: https://pypi.python.org/pypi/mcalf
    :alt: PyPI
.. |Zenodo DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3924527.svg
    :target: https://doi.org/10.5281/zenodo.3924527
    :alt: DOI
.. |Docs Status| image:: https://readthedocs.org/projects/mcalf/badge/?version=latest&style=flat
    :target: https://mcalf.macbride.me/
    :alt: Documentation
.. |GitHub License| image:: https://img.shields.io/github/license/ConorMacBride/mcalf
    :target: ../license.rst
    :alt: License
.. |Contributor Covenant| image:: https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg
    :target: ../code_of_conduct.rst
    :alt: Code of Conduct

.. _Anaconda: https://www.anaconda.com/products/individual#Downloads
.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _new conda environment: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
.. _setup.cfg: https://github.com/ConorMacBride/mcalf/blob/main/setup.cfg
.. _documentation: https://docs.sunpy.org/en/latest/dev_guide/contents/newcomers.html#setting-up-a-development-environment

.. _Conor MacBride: https://macbride.me/

.. _bug reports: https://github.com/ConorMacBride/mcalf/issues
.. _pull requests: https://github.com/ConorMacBride/mcalf/pulls
.. _GitHub repository: https://github.com/ConorMacBride/mcalf

.. _Zenodo DOI: https://doi.org/10.5281/zenodo.3924527
====================================
example1: Basic usage of the package
====================================

``FittingIBIS.ipynb``
---------------------

* `View code <FittingIBIS.ipynb>`_
* :download:`Download FittingIBIS.ipynb <FittingIBIS.ipynb>`

This file is an IPython Notebook containing examples of how to use the package
to accomplish typical tasks.

``FittingIBIS.pro``
-------------------

* :download:`Download FittingIBIS.pro <FittingIBIS.pro>`

This file is similar to ``FittingIBIS.ipynb`` file, except it written is IDL.
It is not recommended to use the IDL wrapper in production, just use it to
explore the code if you are familiar with IDL and not Python.
If you wish to use this package, please use the Python implementation.
IDL is not fully supported in the current version of the code for reasons
such as, the Python tuple datatype cannot be passed from IDL to Python,
resulting in certain function calls not being possible.

``config.yml``
--------------

* :download:`Download config.yml <config.yml>`

This is an example configuration file containing default parameters.
This can be easier than setting the parameters in the code.
The file follows the YAML_ format.

.. _YAML: https://pyyaml.org/wiki/PyYAMLDocumentation
MCALF visualisation
*******************

This sub-package contains:

* Functions to plot the input spectrum and the fitted
  model.
* Functions to plot the spatial distribution and their
  general profile.
* Functions to plot the velocities calculated for a
  spectral imaging scan.

.. automodapi:: mcalf.visualisation
MCALF profiles
**************

This sub-package contains:

* Functions that can be used to model the spectra.
* Voigt profile with a variety of wrappers for different
  applications (`mcalf.profiles.voigt`).
* Gaussian profiles and skew normal distributions
  (`mcalf.profiles.gaussian`).

.. automodapi:: mcalf.profiles

.. automodapi:: mcalf.profiles.voigt

.. automodapi:: mcalf.profiles.gaussian
MCALF utils
***********

This sub-package contains:

* Functions for processing spectra (`mcalf.utils.spec`).
* Functions for smoothing n-dimensional arrays
  (`mcalf.utils.smooth`).
* Functions for masking the input data to limit the
  region computed (`mcalf.utils.mask`).
* Functions for helping with plotting (`mcalf.utils.plot`).
* Classes for managing collections of data (`mcalf.utils.collections`).
* Miscellaneous utility functions (`mcalf.utils.misc`).

.. automodapi:: mcalf.utils

.. automodapi:: mcalf.utils.spec

.. automodapi:: mcalf.utils.smooth

.. automodapi:: mcalf.utils.mask

.. automodapi:: mcalf.utils.plot

.. automodapi:: mcalf.utils.collections

.. automodapi:: mcalf.utils.misc
MCALF
*****

mcalf Package
=============

MCALF: Multi-Component Atmospheric Line Fitting
-----------------------------------------------

MCALF is an open-source Python package for accurately constraining
velocity information from spectral imaging observations using
machine learning techniques.
MCALF models
************

This sub-package contains:

* Base and sample models that can be adapted and fitted to any spectral
  imaging dataset.
* Models optimised for particular data sets that can be used directly.
* Data structures for storing and exporting the fitted parameters, as
  well as simplifying the calculation of velocities.

.. automodapi:: mcalf.models
   :inherited-members:
.. _reference:

**************
Code Reference
**************

.. toctree::
   :maxdepth: 2

   mcalf
   models
   profiles
   visualisation
   utils
