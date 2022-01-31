## Description

## Checklist
- [ ] Follow the [Contributor Guidelines](https://github.com/skypyproject/skypy/blob/main/CONTRIBUTING.rst)
- [ ] Write unit tests
- [ ] Write documentation strings
- [ ] Assign someone from your working team to review this pull request
- [ ] Assign someone from the infrastructure team to review this pull request
---
name: New model proposal
about: Description of new SkyPy model
title: ''
labels: enhancement
assignees: ''

---

## Description

## Inputs
 -

## Outputs
 -

## References
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.

ADR 3: Position sampling and patches of the sky
===============================================

Author: Nicolas Tessore  
Date: 5 February 2021  
Status: Accepted


Context
-------

This ADR addresses two related open problems for SkyPy simulations:

- How can we describe the regions in which to sample the positions of e.g.
  galaxies or supernovae?
- How can we break up the sampling over large fractions of the sky into smaller
  chunks?

The second point is particularly relevant when the simulation becomes so large
that it no longer fits into memory, as well as for parallelisation.


Extension to pipelines
----------------------

This ADR proposes to introduce a new top-level keyword `regions` [alt: `sky`,
`geometry`, `patches`] into pipelines that describes the geometry of the
simulation. For example, to include a "rectangular" and a "circular" region:

```yaml
regions:
- !rectangular [ ... ]
- !circular [ ... ]
```

The `regions` list can be traversed by the pipeline runner to create what are
effectively independent parallel simulations. The list items are objects with
the following interface.


Region interface
----------------

The regions need to support two sets of operations:

- Information queries: For example, there should be a `.area` [alt:
  `.solid_angle`] attribute that returns the solid angle of the region.
- Random sampling: There needs to be at least a `random_point()` [alt:
  `random_position()`, `uniform_in()`] function that can uniformly sample a
  random position from within the region.

When the pipeline runner traverses the list of regions, it can keep track of
the current region in a `$region` reference that can be used where necessary.
For example, to sample from a luminosity function with positions:

```yaml
tables:
  galaxies:
    z, M: !schechter_lf
      ...
      sky_area: $region.area
    ra, dec: !random_point [ $region ]
```


Support for HEALPix maps
------------------------

The above proposal is powerful enough to support advanced features such as
regions that are described by HEALPix maps. There may be a `healpix()` function
that generates a list of regions from HEALPix pixels:

```yaml
regions: !healpix
  nside: 8
```

The resulting list would contain `12 * nside**2 = 768` regions corresponding
to the HEALPix pixels of a map with `nside = 8`.

The function is easily extensible. For example, instead of using all HEALPix
pixels, there might be a footprint that describes a specific survey:

```yaml
regions: !healpix
  mask: survey-footprint-n512.fits
```

The `mask` keyword can be combined with the `nside` parameter to change the
resolution of the mask if requested.

If the HEALPix maps become finely resolved, it may be desirable to combine
several pixels into a single region. There may be a `batch` [alt: `combine`]
keyword for this purpose:

```yaml
regions: !healpix
  nside: 8
  batch: 16
```

The resulting list of regions will contain `768/16 = 48` regions. The `batch`
keyword may also take a quantity of solid angle and automatically choose the
number of pixels to combine accordingly.


Map making
----------

This ADR does not address the problem of how maps will be generated from the
list of regions. For example, a very real use case would be to generate
populations of galaxies and simply count the total number in each HEALPix pixel
to generate a density map. This will be addressed in a separate ADR.


Consequences
------------
The existing `Pipeline` class must be extended to support iterating regions. No
existing interfaces are affected.
# ADR 2: Mpc or Mpc/h
February 6, 2020

## Context
We need to decide on a unit convention as to whether units include the factor /h or not (for instance Mpc or Mpc/h as a unit of distance). For further discussion see e.g. 10.1017/pasa.2013.31

## Decision Drivers
- Flexibility: Mpc/h allows results to be easily propagated across the 'unknown' value of h (0.68 or 0.74 or something else).
- Consistency / least surprise: the default for astropy is Mpc

## Considered Options
- Mpc
- Mpc/h

## Decision Outcome
After [discussion](https://github.com/skypyproject/skypy/issues/23) and offline, Mpc has been chosen to ensure the closest integration and least surprise for astropy.# ADR 1: Considering options for the SkyPy `Model`
January 22, 2020

## Context
Within SkyPy all functions used to create a "simulation" will in practice be taking in some values (either parameters or columns from a table) and creating new column(s) in an output table *or* selecting specific rows from an input table.

The inputs and outputs of these functions are clearly defined so a directed acyclic graph (DAG) can be constructed to determine what order the functions should be run in.

To aid in the creation of the tables and the DAG a helper class or decorator should be used so the person writing the function does not have to worry about the implementation details. This class or decorator is what we are currently referring to as the `Model`.

For clarity in the options below we will assume the following example function:
```python
def redshift_gamma(shape, scale, size):
    """Gamma-distributed redshifts (Smail et al. 1994).

    Sample `size` redshifts from a gamma distribution with the
    given `shape` and `scale` parameters. See `numpy.random.gamma`.
    """

    # redshift distribution
    redshift = np.random.gamma(shape=shape, scale=scale, size=size)

    return redshift
```

## Decision Drivers
- Ease of use: if there is too much boiler plate `Model`s will be annoying to write
- Clarity of implementation: the base `Model` should be easy to read, understand, and debug

## Considered Options

### A base `Model` class
In this implementation all functions must be written inside a class that inherits from the base `Model` class.  A different base class should be used depending on if the function adds a column to a table or selects rows.

The `__init__` method would define all the inputs and outputs and the inherited `__init__` can add this to the DAG.

The `compute` method will contain the custom function.

The `execute` method will call the `compute` method and add the results to the table/mask out rows.

- Ease of use: medium (lots of boiler plate)
- Clarity of implementation: high (Classes are well understood by most developers)

Example:
```python
import BaseModel
import numpy as np

class RedshiftGamma(BaseModel):
    def __init__(self):
        self.inputs = ["shape", "scale", "size"]
        self.outputs = ["redshift"]
        super(RedshiftGamma, self).__init__(self.inputs, self.outputs)
    
    def compute(shape, scale, size):
        """Gamma-distributed redshifts (Smail et al. 1994).

        Sample `size` redshifts from a gamma distribution with the
        given `shape` and `scale` parameters. See `numpy.random.gamma`.
        """

        # redshift distribution
        redshift = np.random.gamma(shape=shape, scale=scale, size=size)

        return redshift
```

### A `Model` decorator
In this implementation all functions must use an `@Model(inputs=[], outputs=[])` decorator. A different decorator should be written for adding columns and selecting rows. The decorator will:

1. Add the `inputs` and `outputs` to the DAG
2. Return a callable function that executes the wrapped function and add the results to the table/mask out rows.

- Ease of use: easy (one line added above a function)
- Clarity of implementation: medium (decorators are functions that return function that return function... This particular implementation will be at least 3 wrappers deep)

Example:
```python
import ModelWrapper
import numpy as np

@ModelWrapper(inputs=["shape", "scale", "size"], outputs=["redshift"])
def redshift_gamma(shape, scale, size):
    """Gamma-distributed redshifts (Smail et al. 1994).

    Sample `size` redshifts from a gamma distribution with the
    given `shape` and `scale` parameters. See `numpy.random.gamma`.
    """

    # redshift distribution
    redshift = np.random.gamma(shape=shape, scale=scale, size=size)

    return redshift
```

### Use the DAG directly
Packages such as [pyungo](https://pypi.org/project/pyungo/) have APIs for most of the functionality we need here with decorators that define `inputs` and `outputs`.  When the compute graph is called all `inputs` and `outputs` are stored in the returned `results` data structure.  Once computed we can write a function that turns this into the final data table.

Also the actual function wrapping only needs to happen for functions contained in the configuration file preventing any un-needed nodes being added to the graph.

We kind of get masking for free here as the DAG does not care if/when the number of rows changes, we just have to be careful when constructing the final table out of the `results`.

- Ease of use: easy (one line added above a function)
- Clarity of implementation: high (we off load this to an existing package that we don't have to maintain)

Example:
```python
def redshift_gamma(shape, scale, size):
    """Gamma-distributed redshifts (Smail et al. 1994).

    Sample `size` redshifts from a gamma distribution with the
    given `shape` and `scale` parameters. See `numpy.random.gamma`.
    """

    # redshift distribution
    redshift = np.random.gamma(shape=shape, scale=scale, size=size)

    return redshift
```

After reading in the config we can wrap all the functions that we need:
```python
from pyungo import Graph

graph = Graph()
graph.register()(redshift_gamma)
res = graph.calculate(data={'shape': 1, 'scale': 1, 'size': 5})
# res is a dict with `shape`, `scale`, `size`, and `redshift`
# this dict can be turned into a table
```

## Decision Outcome
After [discussion](https://github.com/skypyproject/skypy/pull/38) option 3 has been picked.  This will be easiest for developers to write new functions and write clean unit tests.  Within the example given above `pyungo` was just used as an example, other DAG frameworks exist and picking one should be the topic of a different ADR.Citation Guidelines
===================

|JOSS| |Zenodo|


If you use SkyPy for work or research presented in a publication (whether
directly, or as a dependency of another package) we recommend and encourage
the following acknowledgment:

  This research made use of SkyPy, a Python package for forward modeling
  astronomical surveys (Amara et. al., 2021, SkyPy Collaboration, 202x).

where the citations are to our publication in the `Journal of Open Source
Software`_ and the `Zenodo DOI`_ for the specific version of the software that
you used. We also encourage citations within the main text wherever
appropriate. DOIs and BibTeX keys are available through the links above.

.. _Journal of Open Source Software: https://joss.theoj.org/papers/10.21105/joss.03056
.. _Zenodo DOI: https://zenodo.org/record/3755531


.. |JOSS| image:: https://joss.theoj.org/papers/10.21105/joss.03056/status.svg
    :target: https://doi.org/10.21105/joss.03056

.. |Zenodo| image:: https://zenodo.org/badge/doi/10.5281/zenodo.4475347.svg
    :target: https://doi.org/10.5281/zenodo.3755531
Contributor Guidelines
======================

GitHub Workflow
---------------

Fork and Clone the SkyPy Repository
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**You should only need to do this step once**

First *fork* the SkyPy repository. A fork is your own remote copy of the repository on GitHub. To create a fork:

  1. Go to the `SkyPy GitHub Repository <https://github.com/skypyproject/skypy>`_
  2. Click the **Fork** button (in the top-right-hand corner)
  3. Choose where to create the fork, typically your personal GitHub account

Next *clone* your fork. Cloning creates a local copy of the repository on your computer to work with. To clone your fork:

::

   git clone https://github.com/<your-account>/skypy.git


Finally add the ``skypyproject`` repository as a *remote*. This will allow you to fetch changes made to the codebase. To add the ``skypyproject`` remote:

::

  cd skypy
  git remote add skypyproject https://github.com/skypyproject/skypy.git


Create a branch for your new feature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a *branch* off the ``skypyproject`` main branch. Working on unique branches for each new feature simplifies the development, review and merge processes by maintining logical separation. To create a feature branch:

::

  git fetch skypyproject
  git checkout -b <your-branch-name> skypyproject/main


Hack away!
^^^^^^^^^^

Write the new code you would like to contribute and *commit* it to the feature branch on your local repository. Ideally commit small units of work often with clear and descriptive commit messages describing the changes you made. To commit changes to a file:

::

  git add file_containing_your_contribution
  git commit -m 'Your clear and descriptive commit message'


*Push* the contributions in your feature branch to your remote fork on GitHub:

::

  git push origin <your-branch-name>


**Note:** The first time you *push* a feature branch you will probably need to use `--set-upstream origin` to link to your remote fork:

::

  git push --set-upstream origin <your-branch-name>


Open a Pull Request
^^^^^^^^^^^^^^^^^^^

When you feel that work on your new feature is complete, you should create a *Pull Request*. This will propose your work to be merged into the main SkyPy repository.

  1. Go to `SkyPy Pull Requests <https://github.com/skypyproject/skypy/pulls>`_
  2. Click the green **New pull request** button
  3. Click **compare across forks**
  4. Confirm that the base fork is ``skypyproject/skypy`` and the base branch is ``main``
  5. Confirm the head fork is ``<your-account>/skypy`` and the compare branch is ``<your-branch-name>``
  6. Give your pull request a title and fill out the the template for the description
  7. Click the green **Create pull request** button

Status checks
^^^^^^^^^^^^^

A series of automated checks will be run on your pull request, some of which will be required to pass before it can be merged into the main codebase:

  - ``Tests`` (Required) runs the `unit tests`_ in four predefined environments; `latest supported versions`, `oldest supported versions`, `macOS latest supported` and `Windows latest supported`. Click "Details" to view the output including any failures.
  - ``Code Style`` (Required) runs `flake8 <https://flake8.pycqa.org/en/latest/>`__ to check that your code conforms to the `PEP 8 <https://www.python.org/dev/peps/pep-0008/>`_ style guidelines. Click "Details" to view any errors.
  - ``codecov`` reports the test coverage for your pull request; you should aim for `codecov/patch — 100.00%`. Click "Details" to view coverage data.
  - ``docs`` (Required) builds the `docstrings`_ on `readthedocs <https://readthedocs.org/>`_. Click "Details" to view the documentation or the failed build log.

Updating your branch
^^^^^^^^^^^^^^^^^^^^

As you work on your feature, new commits might be made to the ``skypyproject`` main branch. You will need to update your branch with these new commits before your pull request can be accepted. You can achieve this in a few different ways:

  - If your pull request has no conflicts, click **Update branch**
  - If your pull request has conflicts, click **Resolve conflicts**, manually resolve the conflicts and click **Mark as resolved**
  - *merge* the ``skypyproject`` main branch from the command line:

    ::

        git fetch skypyproject
        git merge skypyproject/main

  - *rebase* your feature branch onto the ``skypyproject`` main branch from the command line:
    ::

        git fetch skypyproject
        git rebase skypyproject/main


**Warning**: It is bad practice to *rebase* commits that have already been pushed to a remote such as your fork. Rebasing creates new copies of your commits that can cause the local and remote branches to diverge. ``git push --force`` will **overwrite** the remote branch with your newly rebased local branch. This is strongly discouraged, particularly when working on a shared branch where you could erase a collaborators commits.

For more information about resolving conflicts see the GitHub guides:
  - `Resolving a merge conflict on GitHub <https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/resolving-a-merge-conflict-on-github>`_
  - `Resolving a merge conflict using the command line <https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/resolving-a-merge-conflict-using-the-command-line>`_
  - `About Git rebase <https://help.github.com/en/github/using-git/about-git-rebase>`_

More Information
^^^^^^^^^^^^^^^^

More information regarding the usage of GitHub can be found in the `GitHub Guides <https://guides.github.com/>`_.

Coding Guidelines
-----------------

Before your pull request can be merged into the codebase, it will be reviewed by one of the SkyPy developers and required to pass a number of automated checks. Below are a minimum set of guidelines for developers to follow:

General Guidelines
^^^^^^^^^^^^^^^^^^

- SkyPy is compatible with Python>=3.7 (see `setup.cfg <https://github.com/skypyproject/skypy/blob/main/setup.cfg>`_). SkyPy *does not* support backwards compatibility with Python 2.x; `six`, `__future__` and `2to3` should not be used.
- All contributions should follow the `PEP8 Style Guide for Python Code <https://www.python.org/dev/peps/pep-0008/>`_. We recommend using `flake8 <https://flake8.pycqa.org/>`__ to check your code for PEP8 compliance.
- Importing SkyPy should only depend on having `NumPy <https://www.numpy.org>`_, `SciPy <https://www.scipy.org/>`_ and `Astropy <https://www.astropy.org/>`__ installed.
- Code is grouped into submodules based on broad science areas e.g. `galaxies <https://skypy.readthedocs.io/en/stable/galaxies.html>`_. There is also a `utils <https://skypy.readthedocs.io/en/stable/utils/index.html>`_ submodule for general utility functions.
- For more information see the `Astropy Coding Guidelines <http://docs.astropy.org/en/latest/development/codeguide.html>`_.

Unit Tests
^^^^^^^^^^

Pull requests will require existing unit tests to pass before they can be merged. Additionally, new unit tests should be written for all new public methods and functions. Unit tests for each submodule are contained in subdirectories called ``tests`` and you can run them locally using ``pytest``. For more information see the `Astropy Testing Guidelines <https://docs.astropy.org/en/stable/development/testguide.html>`_.

If your unit tests check the statistical distribution of a random sample, the test outcome itself is a random variable, and the test will fail from time to time. Please mark such tests with the ``@pytest.mark.flaky`` decorator, so that they will be automatically tried again on failure. To prevent non-random test failures from being run multiple times, please isolate random statistical tests and deterministic tests in their own test cases.

Docstrings
^^^^^^^^^^

All public classes, methods and functions require docstrings. You can build documentation locally by installing `sphinx-astropy <https://github.com/astropy/sphinx-astropy>`_ and calling ``make html`` in the ``docs`` subdirectory. Docstrings should include the following sections:

  - Description
  - Parameters
  - Notes
  - References

For more information see the Astropy guide to `Writing Documentation <https://docs.astropy.org/en/stable/development/docguide.html>`_.
Code of Conduct
===============

1. Introduction
###############

The SkyPy Board is committed to creating a friendly and respectful place for research and scientific discussion. All SkyPy Members are expected to show respect and courtesy to others both inside and outside the team.

By being a Member of the SkyPy Project, Members accept to abide by The SkyPy Code of Conduct and accept the procedures by which Code of Conduct incidents are resolved. Any form of behaviour to exclude, intimidate, or harass is a violation of the Code of Conduct.

All Members are required to conform to the Code of Conduct. This Code of Conduct applies to all spaces managed by the Team including, but not limited to, workshops, email lists, and online forums such as GitHub, Slack, and Twitter. The Board may also choose to consider cases that occur outside the Team spaces.

If you believe someone is violating the Code of Conduct, we ask that you report it to the Board who will take the appropriate action to address the situation.

The Board is responsible for enforcing the Code of Conduct.

2. The Code of Conduct
######################

The Board is dedicated to providing a welcoming and supportive environment for all people, regardless of background or identity. As such, we do not tolerate behaviour that is disrespectful or that excludes, intimidates, or harasses  others. We do not tolerate discrimination or harassment based on characteristics that include, but are not limited to, gender identity and expression, sexual orientation, disability, physical appearance, body shape, citizenship, nationality, ethnic or social origin, pregnancy, familial status, veteran status, genetic information, religion or belief (or lack thereof), membership of a national minority, property, age, education, socio-economic status, technical choices, and experience level.

Workshop hosts are expected to assist with the enforcement of the Code of Conduct. Project Members are required to accept the procedures by which the Board resolves any Code of Conduct incidents, which may include storage and processing of their personal information. The Board reserves the right to also investigate possible cases where the Code of Conduct has been violated by Members outside the Project spaces. The Board may decide that the case should be treated as having happened inside a Team space and to take action against individuals in line with the Code of Conduct.

2.1 Expected behaviour
**********************

All team Members are expected to show respect and courtesy to others. All interactions should be professional regardless of platform: either online, by telephone or in-person. To foster a positive and professional learning environment, we encourage the following kinds of behaviours in all Project spaces and all interactions about the Project:

- Use welcoming and inclusive language
- Be respectful of different viewpoints and experiences
- Gracefully accept constructive criticism
- Focus on what is best for the community
- Show courtesy and respect towards other community members

Also, please see the `Four Social Rules`_ for further recommendations.

.. _Four Social Rules: https://www.recurse.com/manual#sub-sec-social-rules

2.2 Unacceptable behaviour
**************************

Examples of unacceptable behaviour by participants in any Project space include:

- Written or verbal comments which have the effect of excluding people on the basis of membership of any specific group
- Causing someone to fear for their safety, such as through stalking, following, or intimidation
- Violent threats or language directed against another person
- The display of sexual or violent images
- Unwelcome sexual attention
- Nonconsensual or unwelcome physical contact
- Sustained disruption of talks, events or communications
- Insults or put-downs
- Sexist, racist, homophobic, transphobic, ableist, or exclusionary jokes
- Excessive swearing
- Incitement to violence, suicide, or self-harm
- Continuing to initiate interaction (including photography or recording) with someone after being asked to stop
- Publication of private communication without consent

3. Complaints procedures
########################

If you feel uncomfortable, don’t wait until it gets worse. `This guide`_ may help provide ways to broach the issues with the SkyPy Board and/or other Members. Please also see articles on Abuse of Power e.g. `48 Ways Managers Abuse Their Power and Destroy Employee Engagement`_.

.. _This guide: https://www.edcc.edu/counseling/documents/Conflict.pdf
.. _48 Ways Managers Abuse Their Power and Destroy Employee Engagement: https://www.linkedin.com/pulse/48-ways-managers-abuse-power-destroy-employee-hanna/

There are three mechanisms for dealing with complaints:

1. Seek advice from the Ombudsperson,
2. Seek advice from the Lead of the Coordination Team,
3. File a complaint directly with all or some of the Board.

If there is any concern about a conflict of interest of any Board member then you may choose which member(s) of the Board to contact. By default the whole Board will discuss the complaint including writing notes on key points in a private SkyPy slack channel visible only to the current members of the Board. Therefore, if you foresee a potential conflict of interest with one (or more) members of the Board, or if would prefer your complaint to remain verbal (not recorded in any way other than in individual current Board members’ memories) then please note this when you register your complaint. When the Board membership changes the discussion of previous complaints will not automatically be shared with new Board members, and will only be shared with new Board members with the consent of the complainee(s).

It is acceptable to contact a member of the Board to say that you would like to have a phone call without giving the reasons for the call, and to ask for the entire contents of the conversation to remain confidential. The only situation in which the Board member would violate your wishes to keep the conversation confidential is if they perceive your safety to be at risk e.g. if you seemed suicidal they would try to contact a local responsible person which may include the emergency services.

Members of the Board and their contact details:

- Brian Nord - Ombudsperson - nord@fnal.gov
- Sarah Bridle - Chair of the Board - pa@sarahbridle.net
- Adam Amara - adam.amara@port.ac.uk

If you feel uncomfortable, talk to someone you trust. If you feel afraid of telling someone you trust, then this is generally a sign that there is a serious problem and you need to get help. Remember, it isn’t always possible to resolve problems on your own. Unfortunately, the spectrum of human nature includes some extremely difficult people who are capable of making your life a misery. If you find yourself in this situation, you have no reason to feel ashamed, or that you should be able to solve the problem on your own.

Your Board is here to help you and will take seriously any violations of the Code of Conduct and take steps to resolve the situation in a timely way.

3.1 Violations of the Code of Conduct
*************************************

If the Board deems that a Member has violated the Code of Conduct then it will first consider whether the offender should be ejected from the Project (i.e. lose their Member status). If not, a note describing the violation will be added to the offending Member’s profile page - which is visible inside the Project. The Board will then review this violation in the context of any previous violations and complaints (including verbal complaints that are not recorded on the offender’s profile page), and may then decide to eject the offender from the Project.

3.2 Ejection from the Project
*****************************

If the Board decides to eject someone from the Project, they will inform the ejected person at the same time as notifying all Members, and remove access to all Project resources, including github, slack and google docs.

The Board understands that taking action against an aggressor could itself cause problems for the victim. The Board will consider this when deciding how to act, and will consult with the victim(s) if this is a concern.

4. Process for updating the Code of Conduct
###########################################

The Code of Conduct is a living document, which is the responsibility of the Board. The Code of Conduct lives on the SkyPy GitHub repository and this is the place where suggestions are welcomed on how it might be updated and improved. Proposed revisions will be considered by the Board. Before they are put into practice they will be emailed to the SkyPy team with a period of 2 weeks for feedback. If no objections are raised, the new Code of Conduct will replace the old one. If there are complaints, the Board will review the issues raised and decide whether (i) to continue with the updated version of the CoC, (ii) keep the old version or (iii) begin the process of drafting a new version.

Credits
#######

This Code of Conduct borrows heavily from the `Carpentries Code of Conduct`_.

.. _Carpentries Code of Conduct: https://docs.carpentries.org/topic_folders/policies/code-of-conduct.html

License
#######

This Code of Conduct is licensed under a Creative Commons Attribution 4.0 International License. We encourage other communities related to ours to use or adapt this code as they see fit. Feedback is welcome.

Update Logs
###########

[update logs will be placed here in future revisions]
===========================================
SkyPy: A package for modelling the Universe
===========================================

|Read the Docs| |GitHub| |Codecov| |Compatibility|

This package contains methods for modelling the Universe, galaxies and the
Milky Way. SkyPy simulates populations of astronomical objects, generating
random realisations of intrinsic and observed properties, with the
intention the simulations can then be compared to data as part of an inference
pipeline.

Currently, SkyPy implements the following modules:

* Galaxies_: morphology, luminosity and redshift distributions
* Pipelines_ to generate populations of astronomical objects

The `full list of features`_ can be found in the `SkyPy Documentation`_.

For more information on the people involved and how SkyPy is developed, please
visit the SkyPy Collaboration website: `http://skypyproject.org`_

.. _Galaxies: https://skypy.readthedocs.io/en/latest/galaxies.html
.. _Pipelines: https://skypy.readthedocs.io/en/latest/pipeline/index.html
.. _full list of features: https://skypy.readthedocs.io/en/latest/feature_list.html
.. _SkyPy Documentation: https://skypy.readthedocs.io/en/latest/
.. _http://skypyproject.org: http://skypyproject.org

Citation
--------

|JOSS| |Zenodo|

If you use SkyPy for work or research presented in a publication please follow
our `Citation Guidelines`_.

.. _Citation Guidelines: CITATION.rst


Installation
------------

|PyPI| |conda-forge|

SkyPy releases are distributed through PyPI_ and conda-forge_. Instructions for
installing SkyPy and its dependencies can be found in the Installation_
section of the documentation.


Examples
--------

SkyPy has a driver script that can run simulation pipelines from the command
line. The documentation contains a description of the Pipeline_ module and
Examples_ that demonstrate how to use it.

.. _PyPI: https://pypi.org/project/skypy/
.. _conda-forge: https://anaconda.org/conda-forge/skypy
.. _Installation: https://skypy.readthedocs.io/en/stable/install.html
.. _Pipeline: https://skypy.readthedocs.io/en/stable/pipeline/index.html
.. _Examples: https://skypy.readthedocs.io/en/stable/examples/index.html


Get in Touch
------------

You are welcome to talk about the SkyPy package and code using our
`Discussions Page`_. For any other questions about the project in general,
please get in touch with the `SkyPy Co-ordinators`_.

 .. _Discussions Page: https://github.com/skypyproject/skypy/discussions
 .. _SkyPy Co-ordinators: mailto:skypy-coordinators@googlegroups.com

Contributing
------------

We love contributions! SkyPy is open source,
built on open source, and we'd love to have you hang out in our community.
For information on how to contribute see our `Contributor Guidelines`_.
All communication relating to The SkyPy Project must meet the standards set out
in the `Code of Conduct`_.

.. _Contributor Guidelines: https://skypy.readthedocs.io/en/latest/developer/contributing.html
.. _Code of Conduct: https://skypy.readthedocs.io/en/stable/project/code_of_conduct.html

.. |PyPI| image:: https://img.shields.io/pypi/v/skypy?label=PyPI&logo=pypi
    :target: https://pypi.python.org/pypi/skypy

.. |conda-forge| image:: https://img.shields.io/conda/vn/conda-forge/skypy?logo=conda-forge
    :target: https://anaconda.org/conda-forge/skypy

.. |Read the Docs| image:: https://img.shields.io/readthedocs/skypy/stable?label=Docs&logo=read%20the%20docs
    :target: https://skypy.readthedocs.io/en/stable

.. |GitHub| image:: https://github.com/skypyproject/skypy/workflows/Tests/badge.svg
    :target: https://github.com/skypyproject/skypy/actions

.. |Compatibility| image:: https://github.com/skypyproject/skypy/actions/workflows/compatibility.yaml/badge.svg
    :target: https://github.com/skypyproject/skypy/actions/workflows/compatibility.yaml

.. |Codecov| image:: https://codecov.io/gh/skypyproject/skypy/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/skypyproject/skypy

.. |Zenodo| image:: https://zenodo.org/badge/221432358.svg
    :target: https://zenodo.org/badge/latestdoi/221432358
    :alt: SkyPy Concept DOI

.. |JOSS| image:: https://joss.theoj.org/papers/d4fac0604318190d6627ab29b568a48d/status.svg
    :target: https://joss.theoj.org/papers/d4fac0604318190d6627ab29b568a48d
Data directory
==============

This directory contains data files included with the package source
code distribution. Note that this is intended only for relatively small files
- large files should be externally hosted and downloaded as needed.
.. _examples-index:

Examples
========

SkyPy implements models for the properties of various astrophysical objects and
a command-line script to run simulation pipelines. In the examples below we
demonstrate the usage of some of these models in both simple python scripts and
configuration files for more complex simulation pipelines.
.. _galaxies_examples:

Galaxies
--------

Galaxy properties including morphology, luminosity and redshift distributions.
###################
Configuration Files
###################

This page outlines how to construct configuration files to run your own routines
with `~skypy.pipeline.Pipeline`.

`SkyPy` is an astrophysical simulation pipeline tool that allows to define any
arbitrary workflow and store data in table format. You may use `SkyPy` `~skypy.pipeline.Pipeline`
to call any function --your own implementation, from any compatible external software or from the `SkyPy library`.
Then `SkyPy` deals with the data dependencies and provides a library of functions to be used with it.

These guidelines start with an example using one of the `SkyPy` functions, and it follows
the concrete YAML syntax necessary for you to write your own configuration files, beyond using `SkyPy`
functions.

SkyPy example
-------------

In this section, we exemplify how you can write a configuration file and use some of the `SkyPy` functions.
In this example, we sample redshifts and magnitudes from the SkyPy luminosity function, `~skypy.galaxies.schechter_lf`.

- `Define variables`:

From the documentation, the parameters for the `~skypy.galaxies.schechter_lf` function are: ``redshift``, the characteristic absolute magnitude ``M_star``, the amplitude ``phi_star``, faint-end slope parameter ``alpha``,
the magnitude limit ``magnitude_limit``, the fraction of sky ``sky_area``, ``cosmology`` and ``noise``.
If you are planning to reuse some of these parameters, you can define them at the top-level of your configuration file.
In our example, we are using ``Astropy`` linear and exponential models for the characteristic absolute magnitude and the amplitude, respectively.
Also, ``noise`` is an optional parameter and you could use its default value by omitting its definition.

  .. code:: yaml

    cosmology: !astropy.cosmology.default_cosmology.get []
    z_range: !numpy.linspace [0, 2, 21]
    M_star: !astropy.modeling.models.Linear1D [-0.9, -20.4]
    phi_star: !astropy.modeling.models.Exponential1D [3e-3, -9.7]
    magnitude_limit: 23
    sky_area: 0.1 deg2

- `Tables and columns`:

You can create a table ``blue_galaxies`` and for now add the columns for redshift and magnitude (note here the ``schechter_lf`` returns a 2D object)

  .. code:: yaml

    tables:
      blue_galaxies:
        redshift, magnitude: !skypy.galaxies.schechter_lf
          redshift: $z_range
      	  M_star: $M_star
      	  phi_star: $phi_star
      	  alpha: -1.3
      	  m_lim: $magnitude_limit
      	  sky_area: $sky_area

`Important:` if cosmology is detected as a parameter but is not set, it automatically uses the cosmology variable defined at the top-level of the file.

This is how the entire configuration file looks like!

.. literalinclude:: luminosity.yml
   :language: yaml

You may now save it as ``luminosity.yml`` and run it using the `SkyPy` `~skypy.pipeline.Pipeline`:

.. plot::
   :include-source: true
   :context: close-figs

    import matplotlib.pyplot as plt
    from skypy.pipeline import Pipeline

    # Execute SkyPy luminosity pipeline
    pipeline = Pipeline.read("luminosity.yml")
    pipeline.execute()

    # Blue population
    skypy_galaxies = pipeline['blue_galaxies']

    # Plot histograms
    fig, axs = plt.subplots(1, 2, figsize=(9, 3))

    axs[0].hist(skypy_galaxies['redshift'], bins=50, histtype='step', color='purple')
    axs[0].set_xlabel(r'$Redshift$')
    axs[0].set_ylabel(r'$\mathrm{N}$')
    axs[0].set_yscale('log')

    axs[1].hist(skypy_galaxies['magnitude'], bins=50, histtype='step', color='green')
    axs[1].set_xlabel(r'$Magnitude$')
    axs[1].set_yscale('log')

    plt.tight_layout()
    plt.show()

You can also run the pipeline directly from the command line and write the outputs to a fits file:

.. code-block:: bash

    $ skypy luminosity.yml luminosity.fits



Don’t forget to check out for more complete examples_!

.. _examples: https://skypy.readthedocs.io/en/stable/examples/index.html


YAML syntax
-----------
YAML_ is a file format designed to be readable by both computers and humans.
Fundamentally, a file written in YAML consists of a set of key-value pairs.
Each pair is written as ``key: value``, where whitespace after the ``:`` is optional.
The hash character ``#`` denotes the start of a comment and all further text on that
line is ignored by the parser.


This guide introduces the main syntax of YAML relevant when writing
a configuration file to use with ``SkyPy``. Essentially, it begins with
definitions of individual variables at the top level, followed by the tables,
and, within the table entries, the features of objects to simulate are included.
Main keywords: ``parameters``, ``cosmology``, ``tables``.


Variables
^^^^^^^^^
* `Variable definition`: a variable is defined as a key-value pair at the top of the file.
  YAML is able to interpret any numeric data with the appropriate type: integer, float, boolean.
  Similarly for lists and dictionaries.
  In addition, SkyPy has added extra functionality to interpret and store Astropy Quantities_.
  Everything else is stored as a string (with or without explicitly using quotes)

  .. code:: yaml

      # YAML interprets
      counter: 100 # An integer
      miles: 1000.0 # A floating point
      name: "Joy" # A string
      planet: Earth # Another string
      mylist: [ 'abc', 789, 2.0e3 ] # A list
      mydict: { 'fruit': 'orange', 'year': 2020 } # A dictionary

      # SkyPy extra functionality
      angle: 10 deg
      distance: 300 kpc


* `Import objects`:
  the SkyPy configuration syntax allows objects to be imported directly from external
  (sub)modules using the ``!`` tag and providing neither a list of arguments or a
  dictionary of keywords. For example, this enables the import and usage of any Astropy cosmology:

  .. code:: yaml

      cosmology: !astropy.cosmology.Planck13 # import the Planck13 object and bind it to the variable named "cosmology"


Parameters
^^^^^^^^^^

* `Parameters definition`: parameters are variables that can be modified at execution.

  For example,

  .. code:: yaml

      parameters:
        hubble_constant: 70
        omega_matter: 0.3


Functions
^^^^^^^^^
* `Function call`: functions are defined as tuples where the first entry is the fully qualified function name tagged with and exclamation mark ``!`` and the second entry is either a list of positional arguments or a dictionary of keyword arguments.

  For example, if you need to call the ``log10()`` and ``linspace()`` NumPy_ functions, then you define the following key-value pairs:

  .. code:: yaml

      log_of_2: !numpy.log10 [2]
      myarray: !numpy.linspace [0, 2.5, 10]

  You can also define parameters of functions with a dictionary of keyword arguments.
  Imagine you want to compute the comoving distance for a range of redshifts and an `Astropy` Planck 2015 cosmology.
  To run it with the `SkyPy` pipeline, call the function and define the parameters as an indented dictionary.

  .. code:: yaml

      comoving_distance: !astropy.cosmology.Planck15.comoving_distance
        z: !numpy.linspace [ 0, 1.3, 10 ]

  Similarly, you can specify the functions arguments as a dictionary:

  .. code:: yaml

      comoving_distance: !astropy.cosmology.Planck15.comoving_distance
        z: !numpy.linspace {start: 0, stop: 1.3, num: 10}

  `N.B.` To call a function with no arguments, you should pass an empty list of
  ``args`` or an empty dictionary of ``kwargs``. For example:

  .. code:: yaml

      cosmo: !astropy.cosmology.default_cosmology.get []


* `Variable reference`: variables can be referenced by their full name tagged with a dollar sign ``$``.
  In the previous example you could also define the variables at the top-level of the file and then reference them:

  .. code:: yaml

      redshift: !numpy.linspace [ 0, 1.3, 10 ]
      comoving_distance: !astropy.cosmology.Planck15.comoving_distance
        z: $redshift

* The `cosmology` to be used by functions within the pipeline only needs to be set up at the top. If a function needs ``cosmology`` as an input, you need not define it again, it is automatically detected.

  For example, calculate the angular size of a galaxy with a given physical size, at a fixed redshift and for a given cosmology:

  .. code:: yaml

      cosmology: !astropy.cosmology.FlatLambdaCDM
        H0: 70
        Om0: 0.3
      size: !skypy.galaxies.morphology.angular_size
        physical_size: 10 kpc
        redshift: 0.2

* `Job completion`: ``.depends`` can be used to force any function call to wait for completion
  of any other job.

  A simple example where, for some reason, the comoving distance needs to be called after
  completion of the angular size function:

  .. code:: yaml

    cosmology: !astropy.cosmology.Planck15
    size: !skypy.galaxies.morphology.angular_size
      physical_size: 10 kpc
      redshift: 0.2
    comoving_distance: !astropy.cosmology.Planck15.comoving_distance
      z: !numpy.linspace [ 0, 1.3, 10 ]
      .depends: size

  By doing so, you force the function call ``redshift`` to be completed before is used to compute the comoving distance.


Tables
^^^^^^

* `Table creation`: a dictionary of table names, each resolving to a dictionary of column names for that table.

  Let us create a table called ``telescope`` with a column to store the width of spectral lines that follow a normal distribution

  .. code:: yaml

      tables:
        telescope:
          spectral_lines: !scipy.stats.norm.rvs
            loc: 550
            scale: 1.6
            size: 100

* `Column addition`: you can add as many columns to a table as you need.
    Imagine you want to add a column for the telescope collecting surface

  .. code:: yaml

      tables:
        telescope:
          spectral_lines: !scipy.stats.norm.rvs
            loc: 550
            scale: 1.6
            size: 100
          collecting_surface: !numpy.random.uniform
            low: 6.9
            high: 7.1
            size: 100

* `Column reference`: columns in the pipeline can be referenced by their full name tagged with a dollar sign ``$``.
  Example: the galaxy mass that follows a lognormal distribution. You can create a table ``galaxies``
  with a column ``mass`` where you sample 10000 object and a second column, ``radius`` which also follows a lognormal distribution
  but the mean depends on how massive the galaxies are:

  .. code:: yaml

    tables:
      galaxies:
        mass: !numpy.random.lognormal
          mean: 5.
          size: 10000
        radius: !numpy.random.lognormal
          mean: $galaxies.mass


* `Multi-column assignment`: multi-column assignment is performed with any 2d-array, where one of the dimensions is interpreted
  as the rows of the table and the second dimension, as separate columns. Or you can do it from a function that returns a tuple.

  We use multi-column assignment in the following example where we sample a two-dimensional array of values from a lognormal distribution and then store them as three columns in a table:

  .. code:: yaml

    tables:
      halos:
        mass, radius, concentration: !numpy.random.lognormal
          size: [10000, 3]


* `Table initialisation`: by default tables are initialised using ``astropy.table.Table()`` however this can be overridden using the ``.init`` keyword to initialise the table with any function call.

  For example, you can stack galaxy properties such as radii and mass:

  .. code:: yaml

    radii: !numpy.logspace [ 1, 2, 100 ]
    mass: !numpy.logspace [ 9, 12, 100 ]
    tables:
      galaxies:
        .init: !astropy.table.vstack [[ $radii, $mass ]]


* `Table reference`: when a function call depends on tables, you need to ensure the referenced table has the necessary content and is not empty.
  You can do that with ``.complete``.

  Example: you want to perform a very simple abundance matching, i.e. painting galaxies within your halos.
  You can create two tables ``halos`` and ``galaxies`` storing the halo mass and galaxy luminosities.
  Then you can stack these two tables and store it in a third table called ``matching``.

  .. code:: yaml

    tables:
      halos:
        halo_mass: !numpy.random.uniform
          low: 1.0e8
          high: 1.0e14
          size: 20
      galaxies:
        luminosity: !numpy.random.uniform
          low: 0.05
          high: 10.0
          size: 20
      matching:
        .init: !astropy.table.hstack
          tables: [ $halos, $galaxies ]
          .depends: [ halos.complete, galaxies.complete ]


.. _YAML: https://yaml.org
.. _NumPy: https://numpy.org
.. _Quantities: https://docs.astropy.org/en/stable/units/
.. _clone(): https://docs.astropy.org/en/stable/api/astropy.cosmology.FLRW.html?highlight=clone#astropy.cosmology.FLRW.clone
#############
Feature List
#############

This page outlines the main features in the `SkyPy` library.

Galaxies
--------

- `Luminosity Distributions`_
- `Morphological Distributions`_
- `Redshift Distributions`_
- `Spectral Energy Distribution Modelling`_
- `Stellar Mass Distributions`_

.. _Luminosity Distributions: https://skypy.readthedocs.io/en/latest/galaxies.html#module-skypy.galaxies.luminosity
.. _Morphological Distributions: https://skypy.readthedocs.io/en/latest/galaxies.html#module-skypy.galaxies.morphology
.. _Redshift Distributions: https://skypy.readthedocs.io/en/latest/galaxies.html#module-skypy.galaxies.redshift
.. _Spectral Energy Distribution Modelling: https://skypy.readthedocs.io/en/latest/galaxies.html#module-skypy.galaxies.spectrum
.. _Stellar Mass Distributions: https://skypy.readthedocs.io/en/latest/galaxies.html#module-skypy.galaxies.stellar_mass
:tocdepth: 3

###################
SkyPy Documentation
###################

This package contains methods for modelling the Universe, galaxies and the
Milky Way. SkyPy simulates populations of astronomical objects, generating
random realisations of intrinsic and observed properties, with the
intention the simulations can then be compared to data as part of an inference
pipeline.

.. Important:: If you use SkyPy for work presented in a publication or talk
   please follow our :doc:`project/citation`.



.. _getting-started:

***************
Getting Started
***************

.. toctree::
   :maxdepth: 1

   install
   feature_list
   configuration_files
   examples/index

.. _user-docs:

******************
User Documentation
******************

Packages
--------

.. toctree::
   :maxdepth: 1

   galaxies
   utils/index

Pipeline
--------

.. toctree::
   :maxdepth: 1

   pipeline/index


.. _developer-docs:

***********************
Developer Documentation
***********************

.. toctree::
   :maxdepth: 1

   developer/contributing


.. _project-details:

***************
Project details
***************

.. toctree::
   :maxdepth: 1

   project/code_of_conduct
   project/citation


*****
Index
*****

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

****************
Acknowledgements
****************

Logo Credit: `Maria Fonseca de la Bella <https://www.linkedin.com/in/mariadelabella>`_
############
Installation
############

This page outlines how to install one of the officially distributed SkyPy
releases and its dependencies, or install and test the latest development
version.

From PyPI
---------

All  SkyPy releases are distributed through the Python Package Index (PyPI_).
To install the latest version use pip_:

.. code:: console

    $ pip install skypy

From conda-forge
----------------

All SkyPy releases are also distributed for conda_ through the `conda-forge`_
channel. To install the latest version for your active conda environment:

.. code:: console

    $ conda install -c conda-forge skypy

From GitHub
-----------

The latest development version of SkyPy can be found on the main branch of
the `skypyproject/skypy`_ GitHub repository. This and any other branch or tag
can be installed directly from GitHub using a recent version of pip:

.. code:: console

    $ pip install skypy@git+https://github.com/skypyproject/skypy.git@main

Dependencies
------------

SkyPy is compatble with Python versions 3.7 or later on Ubuntu, macOS and
Windows operating systems. It has the following core dependencies:

- `astropy <https://www.astropy.org/>`__
- `networkx <https://networkx.github.io/>`_
- `numpy <https://numpy.org/>`_
- `pyyaml <https://pyyaml.org/>`_
- `scipy <https://www.scipy.org/>`_

Installing using pip or conda will automatically install or update these core
dependencies if necessary. SkyPy also has a number of optional dependencies
that enable additional features:

- `h5py <https://www.h5py.org/>`_
- `speclite <https://speclite.readthedocs.io/>`_

To install SkyPy with all optional dependencies using pip:

.. code:: console

    $ pip install skypy[all]

Testing
-------

Once installed, you should be able to import the `skypy` module in python:

.. code:: python

    >>> import skypy

You should also be able to check the installed version number using the `skypy`
command line script:

.. code:: console

    $ skypy --version

You may also want to run the unit tests, for example if you have installed the
development version or you use an unsupported operating system. The unit tests
have the following additional dependencies:

- `pytest-astropy <https://github.com/astropy/pytest-astropy>`_
- `pytest-rerunfailures <https://github.com/pytest-dev/pytest-rerunfailures>`_

The test dependencies can be installed using pip:

.. code:: console

    $ pip install skypy[test]

and the unit tests can then be run using pytest_:

.. code:: console

    $ pytest --pyargs skypy

.. _PyPI: https://pypi.org/project/skypy/
.. _pip: https://pip.pypa.io/
.. _conda: https://docs.conda.io/
.. _conda-forge: https://anaconda.org/conda-forge/skypy
.. _skypyproject/skypy: https://github.com/skypyproject/skypy
.. _pytest: https://docs.pytest.org/
***************************
Galaxies (`skypy.galaxies`)
***************************


Introduction
============


Galaxy Properties
=================

What follows are the physical properties of galaxies simulated by SkyPy, and
the available models for each individual property.


Luminosity
----------

The following models are found in the `skypy.galaxies.luminosity` package.

.. currentmodule:: skypy.galaxies.luminosity
.. autosummary::
   :nosignatures:

   schechter_lf_magnitude


Morphology
----------

The following models are found in the `skypy.galaxies.morphology` package.

.. currentmodule:: skypy.galaxies.morphology
.. autosummary::
   :nosignatures:

   angular_size
   beta_ellipticity
   early_type_lognormal_size
   late_type_lognormal_size
   linear_lognormal_size
   ryden04_ellipticity


Redshift
--------

The following models are found in the `skypy.galaxies.redshift` package.

.. currentmodule:: skypy.galaxies.redshift
.. autosummary::
   :nosignatures:

   redshifts_from_comoving_density
   schechter_lf_redshift
   schechter_smf_redshift
   smail


Spectrum
--------

The following models are found in the `skypy.galaxies.spectrum` package.

SkyPy uses the `speclite <https://speclite.readthedocs.io/>`_ package for
photometric calculations. Some of the following functions take the names of
photometric filters as an input parameter. Users can choose from the available
`Speclite Filters <https://speclite.readthedocs.io/en/latest/filters.html>`_
following the naming syntax described in `speclite.filters.load_filters`, or
create their own named `speclite.filters.FilterResponse`.

.. currentmodule:: skypy.galaxies.spectrum
.. autosummary::
   :nosignatures:

   dirichlet_coefficients
   KCorrectTemplates
   kcorrect


Stellar mass
------------

The following models are found in the `skypy.galaxies.stellar_mass` package.

.. currentmodule:: skypy.galaxies.stellar_mass
.. autosummary::
  :nosignatures:

  schechter_smf_mass


Reference/API
=============

.. automodapi:: skypy.galaxies
.. automodapi:: skypy.galaxies.luminosity
   :include-all-objects:
.. automodapi:: skypy.galaxies.morphology
.. automodapi:: skypy.galaxies.redshift
   :include-all-objects:
.. automodapi:: skypy.galaxies.spectrum
   :include-all-objects:
.. automodapi:: skypy.galaxies.stellar_mass
***************************
Pipeline (`skypy.pipeline`)
***************************

The `~skypy.pipeline` package contains the functionality to run a SkyPy
simulation from end to end. This is implemented in the `~skypy.pipeline.Pipeline`
class and can be called using the :ref:`skypy command line script <skypy-script>`.


.. _skypy-script:

Running ``skypy`` from the command line
=======================================

``skypy`` is a command line script that runs a pipeline of functions defined in
a config file to generate tables of objects and write them to file. For example,
you can use ``skypy`` to run one of the `Examples`_ and write the outputs to
fits files:

.. code-block:: bash

    $ skypy examples/galaxies/sdss_photometry.yml sdss_photometry.fits

To view the progress of the pipeline as it runs you can enable logging using the
`--verbose` flag.

Config files are written in YAML format and read using the
`~skypy.pipeline.load_skypy_yaml` funciton. Each entry in the config specifices
an arbitrary variable, but there are also some particular fields that SkyPy uses:

- `parameters` : Variables that can be modified at execution
- `cosmology` : The cosmology to be used by functions within the pipeline
- `tables` : A dictionary of tables names, each resolving to a dictionary of column names for that table

Every variable can be assigned a fixed value as parsed by pyyaml_.
However, variables and columns can also be evaluated as functions. Fuctions are
defined as tuples where the first entry is the fully qualified function name
tagged with and exclamation mark ``!`` and the second entry is either a list
of positional arguments or a dictionary of keyword arguments. Variables
and columns in the pipeline can also be referenced by their full name tagged
with a dollar sign ``$``. For example:

.. literalinclude:: examples/config.yml
   :language: yaml
   :caption:

.. plot::
  :include-source: false

    import matplotlib.pyplot as plt
    from skypy.pipeline import Pipeline

    pipeline = Pipeline.read('examples/config.yml')
    pipeline.execute()

    z = pipeline['galaxies']['redshift']

    plt.hist(z, histtype='step', density=True, label='redshifts')
    plt.legend()
    plt.xlabel('redshift')

When executing a pipeline, all dependencies are tracked and resolved in order
using a Directed Acylic Graph implemented in networkx_.

.. _Examples: https://skypy.readthedocs.io/en/stable/examples/index.html
.. _pyyaml: https://pyyaml.org/
.. _networkx: https://networkx.github.io/


Using a pipeline from other code
================================

SkyPy pipelines can be executed programmatically from other code. Consider the
following example configuration:

.. literalinclude:: examples/pipeline.yml
   :language: yaml
   :caption:

The `~skypy.pipeline.Pipeline` class can be used to load the configuration file
and run the resulting pipeline. If the configuration defines a `parameters`
section, the definition can be accessed and individual parameter values can be
changed for individual executions of the pipeline:

.. plot::

    import matplotlib.pyplot as plt
    from skypy.pipeline import Pipeline

    # read the example pipeline
    pipeline = Pipeline.read('examples/pipeline.yml')

    # run the pipeline as given
    pipeline.execute()

    # plot the results for the given parameters
    plt.hist(pipeline['galaxy-redshifts'], histtype='step', density=True,
             label='{:.2f}'.format(pipeline.parameters['median-redshift']))

    # change the median redshift parameter in a loop
    for z in [1.2, 1.4, 1.6, 1.8, 2.0]:

        # median redshift parameter
        parameters = {'median-redshift': z}

        # run pipeline with updated parameters
        pipeline.execute(parameters)

        # plot the new results
        plt.hist(pipeline['galaxy-redshifts'], histtype='step', density=True,
                 label='{:.2f}'.format(parameters['median-redshift']))

    # show plot labels
    plt.legend()
    plt.xlabel('redshift')


Reference/API
=============

.. automodapi:: skypy.pipeline
.. include:: ../../CONTRIBUTING.rst
*********************
Utils (`skypy.utils`)
*********************

.. currentmodule:: skypy.utils


Decorators
==========

SkyPy provides a number of convenient decorators to perform common tasks:

- :ref:`broadcast_arguments`
- :ref:`dependent_argument`


.. _broadcast_arguments:

Broadcast arguments to same shape
---------------------------------

The `broadcast_arguments` decorator takes a list of argument names that will
be broadcast to a common shape (using `numpy.broadcast_arrays`) before being
passed to the function.

.. code-block:: python

    import numpy as np

    from skypy.utils import broadcast_arguments

    @broadcast_arguments('a', 'b')
    def a_function_of_arrays(a, b):
        print(np.block([a, b]))

    a = [[1, 2, 3]]      # row vector
    b = [[4], [5], [6]]  # column vector

    a_function_of_arrays(a, b)
    # output: [[1 2 3 4 4 4]
    #          [1 2 3 5 5 5]
    #          [1 2 3 6 6 6]]

Since `broadcast_arguments` requires the final shapes of its arguments, it
should be placed below decorators which modify arguments.


.. _dependent_argument:

Evaluate dependent arguments
----------------------------

Sometimes, the value of one argument (the dependent argument) will be a
function other arguments (the independent arguments). Consider the following
example to plot a parametric surface:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    def plot(x, y, z):
        '''plot the surface given by `z = f(x, y)`'''
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(x, y, z)
        plt.show()

    def f(x, y):
        '''the surface function `z = f(x, y)`'''
        return np.sin(np.hypot(x, y))

    # x and y over a [-5, 5]x[-5, 5] grid
    x, y = np.mgrid[-5.0:5.0:0.1, -5.0:5.0:0.1]

    # compute surface
    z = f(x, y)

    # plot surface
    plot(x, y, z)

The function `plot(x, y, z)` takes three coordinates, where `z` is generated by
a function `f(x, y)`. In this situation, it is often more convenient to pass
the function `f` directly instead of the `z` coordinate, which would require
the receiving function to distinguish whether it is given values or a function.
The `dependent_argument` decorator handles this situation transparently by
taking the name of one dependent argument and the names of one or more
independent arguments, and computing values if a function is provided for the
dependent argument:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    from skypy.utils import dependent_argument

    @dependent_argument('z', 'x', 'y')
    def plot(x, y, z):
        '''plot the surface given by `z = f(x, y)`'''
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(x, y, z)
        plt.show()

    def f(x, y):
        '''the surface function `z = f(x, y)`'''
        return np.sin(np.hypot(x, y))

    # x and y over a [-5, 5]x[-5, 5] grid
    x, y = np.mgrid[-5.0:5.0:0.1, -5.0:5.0:0.1]

    # plot surface
    plot(x, y, f)

Since `dependent_argument` requires the values of all independent arguments, it
should be placed below decorators which modify any of the independent
arguments.


Photometry (`skypy.utils.photometry`)
=====================================

This module contains methods that model spectral energy distributions and
calculate photometric properties.

SkyPy uses the `speclite <https://speclite.readthedocs.io/>`_ package for
photometric calculations. Some of the following functions take the names of
photometric filters as an input parameter. Users can choose from the available
`Speclite Filters <https://speclite.readthedocs.io/en/latest/filters.html>`_
following the naming syntax described in `speclite.filters.load_filters`, or
create their own named `speclite.filters.FilterResponse`.


.. currentmodule:: skypy.utils.photometry
.. autosummary::
   :nosignatures:

   absolute_magnitude_from_luminosity
   luminosity_from_absolute_magnitude
   luminosity_in_band
   mag_ab
   SpectrumTemplates


Random sampling (`skypy.utils.random`)
======================================

.. automodule:: skypy.utils.random


Special functions (`skypy.utils.special`)
=========================================

.. automodule:: skypy.utils.special


Reference/API
=============

.. automodapi:: skypy.utils
.. automodapi:: skypy.utils.photometry
.. include:: ../../CITATION.rst
.. include:: ../../CODE_OF_CONDUCT.rst
