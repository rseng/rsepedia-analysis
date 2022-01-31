# pip requirements files

## Index

- [`default.txt`](default.txt)
  Default requirements
- [`extra.txt`](extra.txt)
  Optional requirements that may require extra steps to install
- [`example.txt`](example.txt)
  Requirements for gallery examples
- [`test.txt`](test.txt)
  Requirements for running test suite
- [`doc.txt`](doc.txt)
  Requirements for building the documentation (see `../doc/`)
- [`developer.txt`](developer.txt)
  Requirements for developers
- [`release.txt`](release.txt)
  Requirements for making releases

## Examples

### Installing requirements

```bash
$ pip install -U -r requirements/default.txt
```

### Running the tests

```bash
$ pip install -U -r requirements/default.txt
$ pip install -U -r requirements/test.txt
```
<!--
Please run black to format your code.
See https://networkx.org/documentation/latest/developer/contribute.html for details.
-->
---
name: Bug report
about: 'Please describe the problem you have encountered'
---

<!-- If you have a general question about NetworkX, please use the discussions tab to create a new discussion -->

<!--- Provide a general summary of the issue in the Title above -->


### Current Behavior
<!--- Tell us what happens instead of the expected behavior -->

### Expected Behavior
<!--- Tell us what should happen -->

### Steps to Reproduce
<!--- Provide a minimal example that reproduces the bug -->

### Environment
<!--- Please provide details about your local environment -->
Python version:
NetworkX version:


### Additional context
<!--- Add any other context about the problem here, screenshots, etc. -->
# Building docs

We use Sphinx for generating the API and reference documentation.

Pre-built versions can be found at

    https://networkx.org/

for both the stable and the latest (i.e., development) releases.

## Instructions

After installing NetworkX and its dependencies, install the Python
packages needed to build the documentation by entering::

    pip install -r requirements/doc.txt

in the root directory.

To build the HTML documentation, enter::

    make html

in the ``doc/`` directory.  This will generate a ``build/html`` subdirectory
containing the built documentation.

To build the PDF documentation, enter::

    make latexpdf

You will need to have LaTeX installed for this.
.. _contributor_guide:

Contributor Guide
=================

.. note::
   This document assumes some familiarity with contributing to open source
   scientific Python projects using GitHub pull requests. If this does not
   describe you, you may first want to see the :ref:`contributing_faq`.

.. _dev_workflow:

Development Workflow
--------------------

1. If you are a first-time contributor:

   * Go to `https://github.com/networkx/networkx
     <https://github.com/networkx/networkx>`_ and click the
     "fork" button to create your own copy of the project.

   * Clone the project to your local computer::

      git clone git@github.com:your-username/networkx.git

   * Navigate to the folder networkx and add the upstream repository::

      git remote add upstream git@github.com:networkx/networkx.git

   * Now, you have remote repositories named:

     - ``upstream``, which refers to the ``networkx`` repository
     - ``origin``, which refers to your personal fork

   * Next, you need to set up your build environment.
     Here are instructions for two popular environment managers:

     * ``venv`` (pip based)

       ::

         # Create a virtualenv named ``networkx-dev`` that lives in the directory of
         # the same name
         python -m venv networkx-dev
         # Activate it
         source networkx-dev/bin/activate
         # Install main development and runtime dependencies of networkx
         pip install -r requirements/default.txt -r requirements/test.txt -r requirements/developer.txt
         #
         # (Optional) Install pygraphviz, pydot, and gdal packages
         # These packages require that you have your system properly configured
         # and what that involves differs on various systems.
         # pip install -r requirements/extra.txt
         #
         # Build and install networkx from source
         pip install -e .
         # Test your installation
         PYTHONPATH=. pytest networkx

     * ``conda`` (Anaconda or Miniconda)

       ::

         # Create a conda environment named ``networkx-dev``
         conda create --name networkx-dev
         # Activate it
         conda activate networkx-dev
         # Install main development and runtime dependencies of networkx
         conda install -c conda-forge --file requirements/default.txt --file requirements/test.txt --file requirements/developer.txt
         #
         # (Optional) Install pygraphviz, pydot, and gdal packages
         # These packages require that you have your system properly configured
         # and what that involves differs on various systems.
         # conda install -c conda-forge --file requirements/extra.txt
         #
         # Install networkx from source
         pip install -e .
         # Test your installation
         PYTHONPATH=. pytest networkx

   * Finally, we recommend you use a pre-commit hook, which runs black when
     you type ``git commit``::

       pre-commit install

2. Develop your contribution:

   * Pull the latest changes from upstream::

      git checkout main
      git pull upstream main

   * Create a branch for the feature you want to work on. Since the
     branch name will appear in the merge message, use a sensible name
     such as 'bugfix-for-issue-1480'::

      git checkout -b bugfix-for-issue-1480

   * Commit locally as you progress (``git add`` and ``git commit``)

3. Test your contribution:

   * Run the test suite locally (see `Testing`_ for details)::

      PYTHONPATH=. pytest networkx

   * Running the tests locally *before* submitting a pull request helps catch
     problems early and reduces the load on the continuous integration
     system.

4. Submit your contribution:

   * Push your changes back to your fork on GitHub::

      git push origin bugfix-for-issue-1480

   * Go to GitHub. The new branch will show up with a green Pull Request
     button---click it.

   * If you want, post on the `mailing list
     <http://groups.google.com/group/networkx-discuss>`_ to explain your changes or
     to ask for review.

5. Review process:

   * Every Pull Request (PR) update triggers a set of `continuous integration
     <https://en.wikipedia.org/wiki/Continuous_integration>`_ services
     that check that the code is up to standards and passes all our tests.
     These checks must pass before your PR can be merged.  If one of the
     checks fails, you can find out why by clicking on the "failed" icon (red
     cross) and inspecting the build and test log.

   * Reviewers (the other developers and interested community members) will
     write inline and/or general comments on your PR to help
     you improve its implementation, documentation, and style.  Every single
     developer working on the project has their code reviewed, and we've come
     to see it as friendly conversation from which we all learn and the
     overall code quality benefits.  Therefore, please don't let the review
     discourage you from contributing: its only aim is to improve the quality
     of project, not to criticize (we are, after all, very grateful for the
     time you're donating!).

   * To update your PR, make your changes on your local repository
     and commit. As soon as those changes are pushed up (to the same branch as
     before) the PR will update automatically.

   .. note::

      If the PR closes an issue, make sure that GitHub knows to automatically
      close the issue when the PR is merged.  For example, if the PR closes
      issue number 1480, you could use the phrase "Fixes #1480" in the PR
      description or commit message.

6. Document changes

   If your change introduces any API modifications, please update
   ``doc/release/release_dev.rst``.

   To set up a function for deprecation:

   - Use a deprecation warning to warn users. For example::

         msg = "curly_hair is deprecated and will be removed in v3.0. Use sum() instead."
         warnings.warn(msg, DeprecationWarning)

   - Add a warning to ``networkx/conftest.py``::

         warnings.filterwarnings(
             "ignore", category=DeprecationWarning, message=<start of message>
         )

   - Add a reminder to ``doc/developer/deprecations.rst`` for the team
     to remove the deprecated functionality in the future. For example:

     .. code-block:: rst

        * In ``utils/misc.py`` remove ``generate_unique_node`` and related tests.

   - Add a note (and a link to the PR) to ``doc/release/release_dev.rst``:

     .. code-block:: rst

        [`#4281 <https://github.com/networkx/networkx/pull/4281>`_]
        Deprecate ``read_yaml`` and ``write_yaml``.


   .. note::

      To reviewers: make sure the merge message has a brief description of the
      change(s) and if the PR closes an issue add, for example, "Closes #123"
      where 123 is the issue number.


Divergence from ``upstream main``
---------------------------------

If GitHub indicates that the branch of your Pull Request can no longer
be merged automatically, merge the main branch into yours::

   git fetch upstream main
   git merge upstream/main

If any conflicts occur, they need to be fixed before continuing.  See
which files are in conflict using::

   git status

Which displays a message like::

   Unmerged paths:
     (use "git add <file>..." to mark resolution)

     both modified:   file_with_conflict.txt

Inside the conflicted file, you'll find sections like these::

   <<<<<<< HEAD
   The way the text looks in your branch
   =======
   The way the text looks in the main branch
   >>>>>>> main

Choose one version of the text that should be kept, and delete the
rest::

   The way the text looks in your branch

Now, add the fixed file::


   git add file_with_conflict.txt

Once you've fixed all merge conflicts, do::

   git commit

.. note::

   Advanced Git users may want to rebase instead of merge,
   but we squash and merge PRs either way.


Guidelines
----------

* All code should have tests.
* All code should be documented, to the same
  `standard <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_
  as NumPy and SciPy.
* All changes are reviewed.  Ask on the
  `mailing list <http://groups.google.com/group/networkx-discuss>`_ if
  you get no response to your pull request.
* Default dependencies are listed in ``requirements/default.txt`` and extra
  (i.e., optional) dependencies are listed in ``requirements/extra.txt``.
  We don't often add new default and extra dependencies.  If you are considering
  adding code that has a dependency, you should first consider adding a gallery
  example.  Typically, new proposed dependencies would first be added as extra
  dependencies.  Extra dependencies should be easy to install on all platforms
  and widely-used.  New default dependencies should be easy to install on all
  platforms, widely-used in the community, and have demonstrated potential for
  wide-spread use in NetworkX.
* Use the following import conventions::

   import numpy as np
   import scipy as sp
   import matplotlib as mpl
   import matplotlib.pyplot as plt
   import pandas as pd
   import networkx as nx

  After importing `sp`` for ``scipy``::

   import scipy as sp

  use the following imports::

   import scipy.linalg  # call as sp.linalg
   import scipy.sparse  # call as sp.sparse
   import scipy.sparse.linalg  # call as sp.sparse.linalg
   import scipy.stats  # call as sp.stats
   import scipy.optimize  # call as sp.optimize

  For example, many libraries have a ``linalg`` subpackage: ``nx.linalg``,
  ``np.linalg``, ``sp.linalg``, ``sp.sparse.linalg``. The above import
  pattern makes the origin of any particular instance of ``linalg`` explicit.

* Use the decorator ``not_implemented_for`` in ``networkx/utils/decorators.py``
  to designate that a function doesn't accept 'directed', 'undirected',
  'multigraph' or 'graph'.  The first argument of the decorated function should
  be the graph object to be checked.

  .. code-block:: python

      @nx.not_implemented_for('directed', 'multigraph')
      def function_not_for_MultiDiGraph(G, others):
          # function not for graphs that are directed *and* multigraph
          pass

      @nx.not_implemented_for('directed')
      @nx.not_implemented_for('multigraph')
      def function_only_for_Graph(G, others):
          # function not for directed graphs *or* for multigraphs
          pass


Testing
-------

``networkx`` has an extensive test suite that ensures correct
execution on your system.  The test suite has to pass before a pull
request can be merged, and tests should be added to cover any
modifications to the code base.
We make use of the `pytest <https://docs.pytest.org/en/latest/>`__
testing framework, with tests located in the various
``networkx/submodule/tests`` folders.

To run all tests::

    $ PYTHONPATH=. pytest networkx

Or the tests for a specific submodule::

    $ PYTHONPATH=. pytest networkx/readwrite

Or tests from a specific file::

    $ PYTHONPATH=. pytest networkx/readwrite/tests/test_yaml.py

Or a single test within that file::

    $ PYTHONPATH=. pytest networkx/readwrite/tests/test_yaml.py::TestYaml::testUndirected

Use ``--doctest-modules`` to run doctests.
For example, run all tests and all doctests using::

    $ PYTHONPATH=. pytest --doctest-modules networkx

Tests for a module should ideally cover all code in that module,
i.e., statement coverage should be at 100%.

To measure the test coverage, run::

  $ PYTHONPATH=. pytest --cov=networkx networkx

This will print a report with one line for each file in `networkx`,
detailing the test coverage::

  Name                                             Stmts   Miss Branch BrPart  Cover
  ----------------------------------------------------------------------------------
  networkx/__init__.py                                33      2      2      1    91%
  networkx/algorithms/__init__.py                    114      0      0      0   100%
  networkx/algorithms/approximation/__init__.py       12      0      0      0   100%
  networkx/algorithms/approximation/clique.py         42      1     18      1    97%
  ...

Adding tests
------------

If you're **new to testing**, see existing test files for examples of things to do.
**Don't let the tests keep you from submitting your contribution!**
If you're not sure how to do this or are having trouble, submit your pull request
anyway.
We will help you create the tests and sort out any kind of problem during code review.

Adding examples
---------------

The gallery examples are managed by
`sphinx-gallery <https://sphinx-gallery.readthedocs.io/>`_.
The source files for the example gallery are ``.py`` scripts in ``examples/`` that
generate one or more figures. They are executed automatically by sphinx-gallery when the
documentation is built. The output is gathered and assembled into the gallery.

You can **add a new** plot by placing a new ``.py`` file in one of the directories inside the
``examples`` directory of the repository. See the other examples to get an idea for the
format.

.. note:: Gallery examples should start with ``plot_``, e.g. ``plot_new_example.py``

General guidelines for making a good gallery plot:

* Examples should highlight a single feature/command.
* Try to make the example as simple as possible.
* Data needed by examples should be included in the same directory and the example script.
* Add comments to explain things are aren't obvious from reading the code.
* Describe the feature that you're showcasing and link to other relevant parts of the
  documentation.

Bugs
----

Please `report bugs on GitHub <https://github.com/networkx/networkx/issues>`_.

Policies
--------

All interactions with the project are subject to the
:doc:`NetworkX code of conduct <code_of_conduct>`.

We also follow these policies:

* :doc:`NetworkX deprecation policy <deprecations>`
* :doc:`Python version support <nep-0029-deprecation_policy>`
Install
=======

NetworkX requires Python 3.8, 3.9, or 3.10.  If you do not already
have a Python environment configured on your computer, please see the
instructions for installing the full `scientific Python stack
<https://scipy.org/install.html>`_.

Below we assume you have the default Python environment already configured on
your computer and you intend to install ``networkx`` inside of it.  If you want
to create and work with Python virtual environments, please follow instructions
on `venv <https://docs.python.org/3/library/venv.html>`_ and `virtual
environments <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_.

First, make sure you have the latest version of ``pip`` (the Python package manager)
installed. If you do not, refer to the `Pip documentation
<https://pip.pypa.io/en/stable/installing/>`_ and install ``pip`` first.

Install the released version
----------------------------

Install the current release of ``networkx`` with ``pip``::

    $ pip install networkx[default]

To upgrade to a newer release use the ``--upgrade`` flag::

    $ pip install --upgrade networkx[default]

If you do not have permission to install software systemwide, you can
install into your user directory using the ``--user`` flag::

    $ pip install --user networkx[default]

If you do not want to install our dependencies (e.g., ``numpy``, ``scipy``, etc.),
you can use::

    $ pip install networkx

This may be helpful if you are using PyPy or you are working on a project that
only needs a limited subset of our functionality and you want to limit the
number of dependencies.

Alternatively, you can manually download ``networkx`` from
`GitHub <https://github.com/networkx/networkx/releases>`_  or
`PyPI <https://pypi.python.org/pypi/networkx>`_.
To install one of these versions, unpack it and run the following from the
top-level source directory using the Terminal::

    $ pip install .[default]

Install the development version
-------------------------------

If you have `Git <https://git-scm.com/>`_ installed on your system, it is also
possible to install the development version of ``networkx``.

Before installing the development version, you may need to uninstall the
standard version of ``networkx`` using ``pip``::

    $ pip uninstall networkx

Then do::

    $ git clone https://github.com/networkx/networkx.git
    $ cd networkx
    $ pip install -e .[default]

The ``pip install -e .[default]`` command allows you to follow the development branch as
it changes by creating links in the right places and installing the command
line scripts to the appropriate locations.

Then, if you want to update ``networkx`` at any time, in the same directory do::

    $ git pull

Extra packages
--------------

.. note::
   Some optional packages (e.g., `gdal`) may require compiling
   C or C++ code.  If you have difficulty installing these packages
   with `pip`, please consult the homepages of those packages.

The following extra packages provide additional functionality. See the
files in the ``requirements/`` directory for information about specific
version requirements.

- `PyGraphviz <http://pygraphviz.github.io/>`_ and
  `pydot <https://github.com/erocarrera/pydot>`_ provide graph drawing
  and graph layout algorithms via `GraphViz <http://graphviz.org/>`_.
- `PyYAML <http://pyyaml.org/>`_ provides YAML format reading and writing.
- `gdal <http://www.gdal.org/>`_ provides shapefile format reading and writing.
- `lxml <http://lxml.de/>`_ used for GraphML XML format.

To install ``networkx`` and extra packages, do::

    $ pip install networkx[default,extra]

To explicitly install all optional packages, do::

    $ pip install pygraphviz pydot pyyaml gdal lxml

Or, install any optional package (e.g., ``pygraphviz``) individually::

    $ pip install pygraphviz

Testing
-------

NetworkX uses the Python ``pytest`` testing package.  You can learn more
about pytest on their `homepage <https://pytest.org>`_.

Test a source distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^

You can test the complete package from the unpacked source directory with::

    pytest networkx

Test an installed package
^^^^^^^^^^^^^^^^^^^^^^^^^

From a shell command prompt you can test the installed package with::

   pytest --pyargs networkx
.. _code_of_conduct:

Code of Conduct
===============


Introduction
------------

This code of conduct applies to all spaces managed by the NetworkX project,
including all public and private mailing lists, issue trackers, wikis, and
any other communication channel used by our community.

This code of conduct should be honored by everyone who participates in
the NetworkX community formally or informally, or claims any affiliation with the
project, in any project-related activities and especially when representing the
project, in any role.

This code is not exhaustive or complete. It serves to distill our common
understanding of a collaborative, shared environment and goals. Please try to
follow this code in spirit as much as in letter, to create a friendly and
productive environment that enriches the surrounding community.


Specific Guidelines
-------------------

We strive to:

1. Be open. We invite anyone to participate in our community. We prefer to use
   public methods of communication for project-related messages, unless
   discussing something sensitive. This applies to messages for help or
   project-related support, too; not only is a public support request much more
   likely to result in an answer to a question, it also ensures that any
   inadvertent mistakes in answering are more easily detected and corrected.

2. Be empathetic, welcoming, friendly, and patient. We work together to resolve
   conflict, and assume good intentions. We may all experience some frustration
   from time to time, but we do not allow frustration to turn into a personal
   attack. A community where people feel uncomfortable or threatened is not a
   productive one.

3. Be collaborative. Our work will be used by other people, and in turn we will
   depend on the work of others. When we make something for the benefit of the
   project, we are willing to explain to others how it works, so that they can
   build on the work to make it even better. Any decision we make will affect
   users and colleagues, and we take those consequences seriously when making
   decisions.

4. Be inquisitive. Nobody knows everything! Asking questions early avoids many
   problems later, so we encourage questions, although we may direct them to
   the appropriate forum. We will try hard to be responsive and helpful.

5. Be careful in the words that we choose.  We are careful and respectful in
   our communication and we take responsibility for our own speech. Be kind to
   others. Do not insult or put down other participants.  We will not accept
   harassment or other exclusionary behaviour, such as:

    - Violent threats or language directed against another person.
    - Sexist, racist, or otherwise discriminatory jokes and language.
    - Posting sexually explicit or violent material.
    - Posting (or threatening to post) other people's personally identifying information ("doxing").
    - Sharing private content, such as emails sent privately or non-publicly,
      or unlogged forums such as IRC channel history, without the sender's consent.
    - Personal insults, especially those using racist or sexist terms.
    - Unwelcome sexual attention.
    - Excessive profanity. Please avoid swearwords; people differ greatly in their sensitivity to swearing.
    - Repeated harassment of others. In general, if someone asks you to stop, then stop.
    - Advocating for, or encouraging, any of the above behaviour.


Diversity Statement
-------------------

The NetworkX project welcomes and encourages participation by everyone. We are
committed to being a community that everyone enjoys being part of. Although
we may not always be able to accommodate each individual's preferences, we try
our best to treat everyone kindly.

No matter how you identify yourself or how others perceive you: we welcome you.
Though no list can hope to be comprehensive, we explicitly honour diversity in:
age, culture, ethnicity, genotype, gender identity or expression, language,
national origin, neurotype, phenotype, political beliefs, profession, race,
religion, sexual orientation, socioeconomic status, subculture and technical
ability.

Though we welcome people fluent in all languages, NetworkX development is
conducted in English.

Standards for behaviour in the NetworkX community are detailed in the Code of
Conduct above. Participants in our community should uphold these standards
in all their interactions and help others to do so as well (see next section).


Reporting Guidelines
--------------------

We know that it is painfully common for internet communication to start at or
devolve into obvious and flagrant abuse.  We also recognize that sometimes
people may have a bad day, or be unaware of some of the guidelines in this Code
of Conduct. Please keep this in mind when deciding on how to respond to a
breach of this Code.

For clearly intentional breaches, report those to the NetworkX Steering Council
(see below). For possibly unintentional breaches, you may reply to the person
and point out this code of conduct (either in public or in private, whatever is
most appropriate). If you would prefer not to do that, please feel free to
report to the NetworkX Steering Council directly, or ask the Council for
advice, in confidence.

You can report issues to the
`NetworkX Steering Council <https://github.com/orgs/networkx/teams/steering-council/members>`__,
at networkx-conduct@groups.io.

If your report involves any members of the Council, or if they feel they have
a conflict of interest in handling it, then they will recuse themselves from
considering your report. Alternatively, if for any reason you feel
uncomfortable making a report to the Council, then you can also contact:

- Senior `NumFOCUS staff <https://numfocus.org/code-of-conduct#persons-responsible>`__: conduct@numfocus.org.


Incident reporting resolution & Code of Conduct enforcement
-----------------------------------------------------------

We will investigate and respond to all complaints. The NetworkX Steering Council
will protect the identity of the reporter, and treat the content of
complaints as confidential (unless the reporter agrees otherwise).

In case of severe and obvious breaches, e.g., personal threat or violent, sexist
or racist language, we will immediately disconnect the originator from NetworkX
communication channels.

In cases not involving clear severe and obvious breaches of this code of
conduct, the process for acting on any received code of conduct violation
report will be:

1. acknowledge report is received
2. reasonable discussion/feedback
3. mediation (if feedback didn't help, and only if both reporter and reportee agree to this)
4. enforcement via transparent decision by the NetworkX Steering Council

The Council will respond to any report as soon as possible, and at most
within 72 hours.


Endnotes
--------

This document is adapted from:

- `SciPy Code of Conduct <http://scipy.github.io/devdocs/dev/conduct/code_of_conduct.html>`_
NetworkX
========

.. image:: https://github.com/networkx/networkx/workflows/test/badge.svg?branch=main
  :target: https://github.com/networkx/networkx/actions?query=workflow%3A%22test%22

.. image:: https://codecov.io/gh/networkx/networkx/branch/main/graph/badge.svg
   :target: https://app.codecov.io/gh/networkx/networkx/branch/main
   
.. image:: https://img.shields.io/github/labels/networkx/networkx/Good%20First%20Issue?color=green&label=Contribute%20&style=flat-square
   :target: https://github.com/networkx/networkx/issues?q=is%3Aopen+is%3Aissue+label%3A%22Good+First+Issue%22
   

NetworkX is a Python package for the creation, manipulation,
and study of the structure, dynamics, and functions
of complex networks.

- **Website (including documentation):** https://networkx.org
- **Mailing list:** https://groups.google.com/forum/#!forum/networkx-discuss
- **Source:** https://github.com/networkx/networkx
- **Bug reports:** https://github.com/networkx/networkx/issues
- **Tutorial:** https://networkx.org/documentation/latest/tutorial.html
- **GitHub Discussions:** https://github.com/networkx/networkx/discussions

Simple example
--------------

Find the shortest path between two nodes in an undirected graph:

.. code:: python

    >>> import networkx as nx
    >>> G = nx.Graph()
    >>> G.add_edge('A', 'B', weight=4)
    >>> G.add_edge('B', 'D', weight=2)
    >>> G.add_edge('A', 'C', weight=3)
    >>> G.add_edge('C', 'D', weight=4)
    >>> nx.shortest_path(G, 'A', 'D', weight='weight')
    ['A', 'B', 'D']

Install
-------

Install the latest version of NetworkX::

    $ pip install networkx

Install with all optional dependencies::

    $ pip install networkx[all]

For additional details, please see `INSTALL.rst`.

Bugs
----

Please report any bugs that you find `here <https://github.com/networkx/networkx/issues>`_.
Or, even better, fork the repository on `GitHub <https://github.com/networkx/networkx>`_
and create a pull request (PR). We welcome all changes, big or small, and we
will help you make the PR if you are new to `git` (just ask on the issue and/or
see `CONTRIBUTING.rst`).

License
-------

Released under the 3-Clause BSD license (see `LICENSE.txt`)::

   Copyright (C) 2004-2022 NetworkX Developers
   Aric Hagberg <hagberg@lanl.gov>
   Dan Schult <dschult@colgate.edu>
   Pieter Swart <swart@lanl.gov>
:orphan:

Geospatial Examples Description
-------------------------------

Functions for reading and writing shapefiles are provided in NetworkX versions <3.0.
However, we recommend that you use the following libraries when working
with geospatial data (including reading and writing shapefiles).

Geospatial Python Libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~

`GeoPandas <https://geopandas.readthedocs.io/>`__ provides
interoperability between geospatial formats and storage mechanisms
(e.g., databases) and Pandas data frames for tabular-oriented processing
of spatial data, as well as a wide array of supporting functionality
including spatial indices, spatial predicates (e.g., test if geometries
intersect each other), spatial operations (e.g., the area of overlap
between intersecting polygons), and more.

See the following examples that use GeoPandas:

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example shows how to build a delaunay graph (plus its dual, the set of Voronoi polygons) f...">

.. only:: html

 .. figure:: /auto_examples/geospatial/images/thumb/sphx_glr_plot_delaunay_thumb.png
     :alt: Delaunay graphs from geographic points

     :ref:`sphx_glr_auto_examples_geospatial_plot_delaunay.py`

.. raw:: html

    </div>

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example shows how to build a graph from a set of geographic lines (sometimes called &quot;lines...">

.. only:: html

 .. figure:: /auto_examples/geospatial/images/thumb/sphx_glr_plot_lines_thumb.png
     :alt: Graphs from a set of lines

     :ref:`sphx_glr_auto_examples_geospatial_plot_lines.py`

.. raw:: html

    </div>

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example shows how to build a graph from a set of polygons using PySAL and geopandas. We&#x27;ll...">

.. only:: html

 .. figure:: /auto_examples/geospatial/images/thumb/sphx_glr_plot_polygons_thumb.png
     :alt: Graphs from Polygons

     :ref:`sphx_glr_auto_examples_geospatial_plot_polygons.py`

.. raw:: html

    </div>

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example shows how to build a graph from a set of points using PySAL and geopandas. In this...">

.. only:: html

 .. figure:: /auto_examples/geospatial/images/thumb/sphx_glr_plot_points_thumb.png
     :alt: Graphs from geographic points

     :ref:`sphx_glr_auto_examples_geospatial_plot_points.py`

.. raw:: html

    </div>

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example shows how to use OSMnx to download and model a street network from OpenStreetMap, ...">

.. only:: html

 .. figure:: /auto_examples/geospatial/images/thumb/sphx_glr_plot_osmnx_thumb.png
     :alt: OpenStreetMap with OSMnx

     :ref:`sphx_glr_auto_examples_geospatial_plot_osmnx.py`

.. raw:: html

    </div>

.. raw:: html

    <div class="sphx-glr-clear"></div>

`PySAL <https://pysal.org/>`__ provides a rich suite of spatial analysis
algorithms. From a network analysis context, `spatial
weights <https://pysal.org/libpysal/api.html#spatial-weights>`__
provide… (Levi please add more here).

See the following examples that use PySAL:

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example shows how to build a delaunay graph (plus its dual, the set of Voronoi polygons) f...">

.. only:: html

 .. figure:: /auto_examples/geospatial/images/thumb/sphx_glr_plot_delaunay_thumb.png
     :alt: Delaunay graphs from geographic points

     :ref:`sphx_glr_auto_examples_geospatial_plot_delaunay.py`

.. raw:: html

    </div>

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example shows how to build a graph from a set of geographic lines (sometimes called &quot;lines...">

.. only:: html

 .. figure:: /auto_examples/geospatial/images/thumb/sphx_glr_plot_lines_thumb.png
     :alt: Graphs from a set of lines

     :ref:`sphx_glr_auto_examples_geospatial_plot_lines.py`

.. raw:: html

    </div>

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example shows how to build a graph from a set of polygons using PySAL and geopandas. We&#x27;ll...">

.. only:: html

 .. figure:: /auto_examples/geospatial/images/thumb/sphx_glr_plot_polygons_thumb.png
     :alt: Graphs from Polygons

     :ref:`sphx_glr_auto_examples_geospatial_plot_polygons.py`

.. raw:: html

    </div>

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example shows how to build a graph from a set of points using PySAL and geopandas. In this...">

.. only:: html

 .. figure:: /auto_examples/geospatial/images/thumb/sphx_glr_plot_points_thumb.png
     :alt: Graphs from geographic points

     :ref:`sphx_glr_auto_examples_geospatial_plot_points.py`

.. raw:: html

    </div>

.. raw:: html

    <div class="sphx-glr-clear"></div>

`momepy <http://docs.momepy.org/en/stable/>`__ builds on top of
GeoPandas and PySAL to provide a suite of algorithms focused on urban
morphology. From a network analysis context, momepy enables you to
convert your line geometry to `networkx.MultiGraph` and back to 
`geopandas.GeoDataFrame` and apply a range of analytical functions aiming at 
morphological description of (street) network configurations.

See the following examples that use momepy:

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example shows how to build a graph from a set of geographic lines (sometimes called &quot;lines...">

.. only:: html

 .. figure:: /auto_examples/geospatial/images/thumb/sphx_glr_plot_lines_thumb.png
     :alt: Graphs from a set of lines

     :ref:`sphx_glr_auto_examples_geospatial_plot_lines.py`

.. raw:: html

    </div>

.. raw:: html

    <div class="sphx-glr-clear"></div>

`OSMnx <https://osmnx.readthedocs.io/>`__ provides a set of tools to retrieve,
model, project, analyze, and visualize OpenStreetMap street networks (and any
other networked infrastructure) as `networkx.MultiDiGraph` objects, and convert
these MultiDiGraphs to/from `geopandas.GeoDataFrame`. It can automatically add
node/edge attributes for: elevation and grade (using the Google Maps Elevation
API), edge travel speed, edge traversal time, and edge bearing. It can also
retrieve any other spatial data from OSM (such as building footprints, public
parks, schools, transit stops, etc) as Geopandas GeoDataFrames.

See the following examples that use OSMnx:

.. raw:: html

    <div class="sphx-glr-thumbcontainer" tooltip="This example shows how to use OSMnx to download and model a street network from OpenStreetMap, ...">

.. only:: html

 .. figure:: /auto_examples/geospatial/images/thumb/sphx_glr_plot_osmnx_thumb.png
     :alt: OpenStreetMap with OSMnx

     :ref:`sphx_glr_auto_examples_geospatial_plot_osmnx.py`

.. raw:: html

    </div>

.. raw:: html

    <div class="sphx-glr-clear"></div>

Key Concepts
~~~~~~~~~~~~

One of the essential tasks in network analysis of geospatial data is
defining the spatial relationships between spatial features (points,
lines, or polygons).

``PySAL`` provides several ways of representing these spatial
relationships between features using the concept of spatial weights.
These include relationships such as ``Queen``, ``Rook``, ...
(Levi please add more here with a brief explanation of each).

``momepy`` allows representation of street networks as both primal
and dual graphs (in a street network analysis sense). The primal approach
turns intersections into Graph nodes and street segments into edges,
a format which is used for a majority of morphological studies. The dual 
approach uses street segments as nodes and intersection topology
as edges, which allows encoding of angular information (i.e an analysis
can be weighted by angles between street segments instead of their length).

``OSMnx`` represents street networks as primal, nonplanar, directed graphs with
possible self-loops and parallel edges to model real-world street network form
and flow. Nodes represent intersections and dead-ends, and edges represent the
street segments linking them. Details of OSMnx's modeling methodology are
available at https://doi.org/10.1016/j.compenvurbsys.2017.05.004

Learn More
~~~~~~~~~~

To learn more see `Geographic Data Science with PySAL and the PyData Stack
<https://geographicdata.science/book/intro.html>`_.
Tutorial
========

.. currentmodule:: networkx

This guide can help you start working with NetworkX.

Creating a graph
----------------

Create an empty graph with no nodes and no edges.

.. nbplot::

    >>> import networkx as nx
    >>> G = nx.Graph()

By definition, a :class:`Graph` is a collection of nodes (vertices) along with
identified pairs of nodes (called edges, links, etc).  In NetworkX, nodes can
be any :py:term:`hashable` object e.g., a text string, an image, an XML object,
another Graph, a customized node object, etc.

.. note:: Python's ``None`` object is not allowed to be used as a node. It
    determines whether optional function arguments have been assigned in many
    functions.

Nodes
-----

The graph ``G`` can be grown in several ways.  NetworkX includes many
:doc:`graph generator functions <reference/generators>` and
:doc:`facilities to read and write graphs in many formats <reference/readwrite/index>`.
To get started though we'll look at simple manipulations.  You can add one node
at a time,

.. nbplot::

    >>> G.add_node(1)

or add nodes from any :py:term:`iterable` container, such as a list

.. nbplot::

    >>> G.add_nodes_from([2, 3])

You can also add nodes along with node
attributes if your container yields 2-tuples of the form 
``(node, node_attribute_dict)``::

    >>> G.add_nodes_from([
    ...     (4, {"color": "red"}),
    ...     (5, {"color": "green"}),
    ... ])

Node attributes are discussed further :ref:`below <attributes>`.

Nodes from one graph can be incorporated into another:

.. nbplot::

    >>> H = nx.path_graph(10)
    >>> G.add_nodes_from(H)

``G`` now contains the nodes of ``H`` as nodes of ``G``.
In contrast, you could use the graph ``H`` as a node in ``G``.

.. nbplot::

    >>> G.add_node(H)

The graph ``G`` now contains ``H`` as a node.  This flexibility is very powerful as
it allows graphs of graphs, graphs of files, graphs of functions and much more.
It is worth thinking about how to structure your application so that the nodes
are useful entities.  Of course you can always use a unique identifier in ``G``
and have a separate dictionary keyed by identifier to the node information if
you prefer.

.. note:: You should not change the node object if the hash depends
   on its contents.

Edges
-----

``G`` can also be grown by adding one edge at a time,

.. nbplot::

    >>> G.add_edge(1, 2)
    >>> e = (2, 3)
    >>> G.add_edge(*e)  # unpack edge tuple*

by adding a list of edges,

.. nbplot::

    >>> G.add_edges_from([(1, 2), (1, 3)])

or by adding any :term:`ebunch` of edges.  An *ebunch* is any iterable
container of edge-tuples.  An edge-tuple can be a 2-tuple of nodes or a 3-tuple
with 2 nodes followed by an edge attribute dictionary, e.g.,
``(2, 3, {'weight': 3.1415})``.  Edge attributes are discussed further
:ref:`below <attributes>`.

.. nbplot::

    >>> G.add_edges_from(H.edges)

There are no complaints when adding existing nodes or edges. For example,
after removing all nodes and edges,

.. nbplot::

    >>> G.clear()

we add new nodes/edges and NetworkX quietly ignores any that are
already present.

.. nbplot::

    >>> G.add_edges_from([(1, 2), (1, 3)])
    >>> G.add_node(1)
    >>> G.add_edge(1, 2)
    >>> G.add_node("spam")        # adds node "spam"
    >>> G.add_nodes_from("spam")  # adds 4 nodes: 's', 'p', 'a', 'm'
    >>> G.add_edge(3, 'm')

At this stage the graph ``G`` consists of 8 nodes and 3 edges, as can be seen by:

.. nbplot::

    >>> G.number_of_nodes()
    8
    >>> G.number_of_edges()
    3

.. note:: 
   
   The order of adjacency reporting (e.g., :meth:`G.adj <networkx.Graph.adj>`,
   :meth:`G.successors <networkx.DiGraph.successors>`,
   :meth:`G.predecessors <networkx.DiGraph.predecessors>`) is the order of
   edge addition. However, the order of G.edges is the order of the adjacencies
   which includes both the order of the nodes and each 
   node's adjacencies. See example below:

.. nbplot::

    >>> DG = nx.DiGraph()
    >>> DG.add_edge(2, 1)   # adds the nodes in order 2, 1
    >>> DG.add_edge(1, 3)
    >>> DG.add_edge(2, 4)
    >>> DG.add_edge(1, 2)
    >>> assert list(DG.successors(2)) == [1, 4]
    >>> assert list(DG.edges) == [(2, 1), (2, 4), (1, 3), (1, 2)]

Examining elements of a graph
-----------------------------

We can examine the nodes and edges. Four basic graph properties facilitate
reporting: ``G.nodes``, ``G.edges``, ``G.adj`` and ``G.degree``.  These
are set-like views of the nodes, edges, neighbors (adjacencies), and degrees
of nodes in a graph. They offer a continually updated read-only view into
the graph structure. They are also dict-like in that you can look up node
and edge data attributes via the views and iterate with data attributes
using methods ``.items()``, ``.data()``.
If you want a specific container type instead of a view, you can specify one.
Here we use lists, though sets, dicts, tuples and other containers may be
better in other contexts.

.. nbplot::

    >>> list(G.nodes)
    [1, 2, 3, 'spam', 's', 'p', 'a', 'm']
    >>> list(G.edges)
    [(1, 2), (1, 3), (3, 'm')]
    >>> list(G.adj[1])  # or list(G.neighbors(1))
    [2, 3]
    >>> G.degree[1]  # the number of edges incident to 1
    2

One can specify to report the edges and degree from a subset of all nodes
using an :term:`nbunch`. An *nbunch* is any of: ``None`` (meaning all nodes),
a node, or an iterable container of nodes that is not itself a node in the
graph.

.. nbplot::

    >>> G.edges([2, 'm'])
    EdgeDataView([(2, 1), ('m', 3)])
    >>> G.degree([2, 3])
    DegreeView({2: 1, 3: 2})

Removing elements from a graph
------------------------------

One can remove nodes and edges from the graph in a similar fashion to adding.
Use methods
:meth:`Graph.remove_node`,
:meth:`Graph.remove_nodes_from`,
:meth:`Graph.remove_edge`
and
:meth:`Graph.remove_edges_from`, e.g.

.. nbplot::

    >>> G.remove_node(2)
    >>> G.remove_nodes_from("spam")
    >>> list(G.nodes)
    [1, 3, 'spam']
    >>> G.remove_edge(1, 3)

Using the graph constructors
----------------------------

Graph objects do not have to be built up incrementally - data specifying
graph structure can be passed directly to the constructors of the various
graph classes.
When creating a graph structure by instantiating one of the graph
classes you can specify data in several formats.

.. nbplot::

    >>> G.add_edge(1, 2)
    >>> H = nx.DiGraph(G)  # create a DiGraph using the connections from G
    >>> list(H.edges())
    [(1, 2), (2, 1)]
    >>> edgelist = [(0, 1), (1, 2), (2, 3)]
    >>> H = nx.Graph(edgelist)  # create a graph from an edge list
    >>> list(H.edges())
    [(0, 1), (1, 2), (2, 3)]
    >>> adjacency_dict = {0: (1, 2), 1: (0, 2), 2: (0, 1)}
    >>> H = nx.Graph(adjacency_dict)  # create a Graph dict mapping nodes to nbrs
    >>> list(H.edges())
    [(0, 1), (0, 2), (1, 2)]

What to use as nodes and edges
------------------------------

You might notice that nodes and edges are not specified as NetworkX
objects.  This leaves you free to use meaningful items as nodes and
edges. The most common choices are numbers or strings, but a node can
be any hashable object (except ``None``), and an edge can be associated
with any object ``x`` using ``G.add_edge(n1, n2, object=x)``.

As an example, ``n1`` and ``n2`` could be protein objects from the RCSB Protein
Data Bank, and ``x`` could refer to an XML record of publications detailing
experimental observations of their interaction.

We have found this power quite useful, but its abuse
can lead to surprising behavior unless one is familiar with Python.
If in doubt, consider using :func:`~relabel.convert_node_labels_to_integers` to obtain
a more traditional graph with integer labels.

Accessing edges and neighbors
-----------------------------

In addition to the views :attr:`Graph.edges`, and :attr:`Graph.adj`,
access to edges and neighbors is possible using subscript notation.

.. nbplot::

    >>> G = nx.Graph([(1, 2, {"color": "yellow"})])
    >>> G[1]  # same as G.adj[1]
    AtlasView({2: {'color': 'yellow'}})
    >>> G[1][2]
    {'color': 'yellow'}
    >>> G.edges[1, 2]
    {'color': 'yellow'}

You can get/set the attributes of an edge using subscript notation
if the edge already exists.

.. nbplot::

    >>> G.add_edge(1, 3)
    >>> G[1][3]['color'] = "blue"
    >>> G.edges[1, 2]['color'] = "red"
    >>> G.edges[1, 2]
    {'color': 'red'}

Fast examination of all (node, adjacency) pairs is achieved using
``G.adjacency()``, or ``G.adj.items()``.
Note that for undirected graphs, adjacency iteration sees each edge twice.

.. nbplot::

    >>> FG = nx.Graph()
    >>> FG.add_weighted_edges_from([(1, 2, 0.125), (1, 3, 0.75), (2, 4, 1.2), (3, 4, 0.375)])
    >>> for n, nbrs in FG.adj.items():
    ...    for nbr, eattr in nbrs.items():
    ...        wt = eattr['weight']
    ...        if wt < 0.5: print(f"({n}, {nbr}, {wt:.3})")
    (1, 2, 0.125)
    (2, 1, 0.125)
    (3, 4, 0.375)
    (4, 3, 0.375)

Convenient access to all edges is achieved with the edges property.

.. nbplot::

    >>> for (u, v, wt) in FG.edges.data('weight'):
    ...     if wt < 0.5:
    ...         print(f"({u}, {v}, {wt:.3})")
    (1, 2, 0.125)
    (3, 4, 0.375)

.. _attributes:

Adding attributes to graphs, nodes, and edges
---------------------------------------------

Attributes such as weights, labels, colors, or whatever Python object you like,
can be attached to graphs, nodes, or edges.

Each graph, node, and edge can hold key/value attribute pairs in an associated
attribute dictionary (the keys must be hashable).  By default these are empty,
but attributes can be added or changed using ``add_edge``, ``add_node`` or direct
manipulation of the attribute dictionaries named ``G.graph``, ``G.nodes``, and
``G.edges`` for a graph ``G``.

Graph attributes
~~~~~~~~~~~~~~~~

Assign graph attributes when creating a new graph

.. nbplot::

    >>> G = nx.Graph(day="Friday")
    >>> G.graph
    {'day': 'Friday'}

Or you can modify attributes later

.. nbplot::

    >>> G.graph['day'] = "Monday"
    >>> G.graph
    {'day': 'Monday'}

Node attributes
~~~~~~~~~~~~~~~

Add node attributes using ``add_node()``, ``add_nodes_from()``, or ``G.nodes``

.. nbplot::

    >>> G.add_node(1, time='5pm')
    >>> G.add_nodes_from([3], time='2pm')
    >>> G.nodes[1]
    {'time': '5pm'}
    >>> G.nodes[1]['room'] = 714
    >>> G.nodes.data()
    NodeDataView({1: {'time': '5pm', 'room': 714}, 3: {'time': '2pm'}})

Note that adding a node to ``G.nodes`` does not add it to the graph, use
``G.add_node()`` to add new nodes. Similarly for edges.

Edge Attributes
~~~~~~~~~~~~~~~

Add/change edge attributes using ``add_edge()``, ``add_edges_from()``,
or subscript notation.

.. nbplot::

    >>> G.add_edge(1, 2, weight=4.7 )
    >>> G.add_edges_from([(3, 4), (4, 5)], color='red')
    >>> G.add_edges_from([(1, 2, {'color': 'blue'}), (2, 3, {'weight': 8})])
    >>> G[1][2]['weight'] = 4.7
    >>> G.edges[3, 4]['weight'] = 4.2

The special attribute ``weight`` should be numeric as it is used by
algorithms requiring weighted edges.

Directed graphs
---------------

The :class:`DiGraph` class provides additional methods and properties specific
to directed edges, e.g.,
:attr:`DiGraph.out_edges`, :attr:`DiGraph.in_degree`,
`DiGraph.predecessors`, `DiGraph.successors` etc.
To allow algorithms to work with both classes easily, the directed versions of
:meth:`neighbors <DiGraph.neighbors>` is equivalent to
`successors <DiGraph.successors>` while `~DiGraph.degree` reports the sum
of `~DiGraph.in_degree` and `~DiGraph.out_degree` even though that may feel inconsistent at times.

.. nbplot::

    >>> DG = nx.DiGraph()
    >>> DG.add_weighted_edges_from([(1, 2, 0.5), (3, 1, 0.75)])
    >>> DG.out_degree(1, weight='weight')
    0.5
    >>> DG.degree(1, weight='weight')
    1.25
    >>> list(DG.successors(1))
    [2]
    >>> list(DG.neighbors(1))
    [2]

Some algorithms work only for directed graphs and others are not well
defined for directed graphs.  Indeed the tendency to lump directed
and undirected graphs together is dangerous.  If you want to treat
a directed graph as undirected for some measurement you should probably
convert it using :meth:`Graph.to_undirected` or with

.. nbplot::

    >>> H = nx.Graph(G)  # create an undirected graph H from a directed graph G

Multigraphs
-----------

NetworkX provides classes for graphs which allow multiple edges
between any pair of nodes.  The :class:`MultiGraph` and
:class:`MultiDiGraph`
classes allow you to add the same edge twice, possibly with different
edge data.  This can be powerful for some applications, but many
algorithms are not well defined on such graphs.
Where results are well defined,
e.g., :meth:`MultiGraph.degree` we provide the function.  Otherwise you
should convert to a standard graph in a way that makes the measurement
well defined.

.. nbplot::

    >>> MG = nx.MultiGraph()
    >>> MG.add_weighted_edges_from([(1, 2, 0.5), (1, 2, 0.75), (2, 3, 0.5)])
    >>> dict(MG.degree(weight='weight'))
    {1: 1.25, 2: 1.75, 3: 0.5}
    >>> GG = nx.Graph()
    >>> for n, nbrs in MG.adjacency():
    ...    for nbr, edict in nbrs.items():
    ...        minvalue = min([d['weight'] for d in edict.values()])
    ...        GG.add_edge(n, nbr, weight = minvalue)
    ...
    >>> nx.shortest_path(GG, 1, 3)
    [1, 2, 3]

Graph generators and graph operations
-------------------------------------

In addition to constructing graphs node-by-node or edge-by-edge, they
can also be generated by

1. Applying classic graph operations, such as:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::

    ~networkx.classes.function.subgraph
    ~networkx.algorithms.operators.binary.union
    ~networkx.algorithms.operators.binary.disjoint_union
    ~networkx.algorithms.operators.product.cartesian_product
    ~networkx.algorithms.operators.binary.compose
    ~networkx.algorithms.operators.unary.complement
    ~networkx.classes.function.create_empty_copy
    ~networkx.classes.function.to_undirected
    ~networkx.classes.function.to_directed

2. Using a call to one of the classic small graphs, e.g.,
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::

    ~networkx.generators.small.petersen_graph
    ~networkx.generators.small.tutte_graph
    ~networkx.generators.small.sedgewick_maze_graph
    ~networkx.generators.small.tetrahedral_graph

3. Using a (constructive) generator for a classic graph, e.g.,
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::

    ~networkx.generators.classic.complete_graph
    ~networkx.algorithms.bipartite.generators.complete_bipartite_graph
    ~networkx.generators.classic.barbell_graph
    ~networkx.generators.classic.lollipop_graph

like so:

.. nbplot::

    >>> K_5 = nx.complete_graph(5)
    >>> K_3_5 = nx.complete_bipartite_graph(3, 5)
    >>> barbell = nx.barbell_graph(10, 10)
    >>> lollipop = nx.lollipop_graph(10, 20)

4. Using a stochastic graph generator, e.g,
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::

    ~networkx.generators.random_graphs.erdos_renyi_graph
    ~networkx.generators.random_graphs.watts_strogatz_graph
    ~networkx.generators.random_graphs.barabasi_albert_graph
    ~networkx.generators.random_graphs.random_lobster

like so:

.. nbplot::

    >>> er = nx.erdos_renyi_graph(100, 0.15)
    >>> ws = nx.watts_strogatz_graph(30, 3, 0.1)
    >>> ba = nx.barabasi_albert_graph(100, 5)
    >>> red = nx.random_lobster(100, 0.9, 0.9)

5. Reading a graph stored in a file using common graph formats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NetworkX supports many popular formats, such as edge lists, adjacency lists,
GML, GraphML, pickle, LEDA and others.

.. nbplot::

    >>> nx.write_gml(red, "path.to.file")
    >>> mygraph = nx.read_gml("path.to.file")

For details on graph formats see :doc:`/reference/readwrite/index`
and for graph generator functions see :doc:`/reference/generators`

Analyzing graphs
----------------

The structure of ``G`` can be analyzed using various graph-theoretic
functions such as:

.. nbplot::

    >>> G = nx.Graph()
    >>> G.add_edges_from([(1, 2), (1, 3)])
    >>> G.add_node("spam")       # adds node "spam"
    >>> list(nx.connected_components(G))
    [{1, 2, 3}, {'spam'}]
    >>> sorted(d for n, d in G.degree())
    [0, 1, 1, 2]
    >>> nx.clustering(G)
    {1: 0, 2: 0, 3: 0, 'spam': 0}

Some functions with large output iterate over (node, value) 2-tuples.
These are easily stored in a `dict` structure if you desire.

.. nbplot::

    >>> sp = dict(nx.all_pairs_shortest_path(G))
    >>> sp[3]
    {3: [3], 1: [3, 1], 2: [3, 1, 2]}

See :doc:`/reference/algorithms/index` for details on graph algorithms
supported.

Drawing graphs
--------------

NetworkX is not primarily a graph drawing package but basic drawing with
Matplotlib as well as an interface to use the open source Graphviz software
package are included.  These are part of the :doc:`networkx.drawing <reference/drawing>`
module and will be imported if possible.

First import Matplotlib's plot interface (pylab works too)

.. nbplot::

    >>> import matplotlib.pyplot as plt

To test if the import of `~networkx.drawing.nx_pylab` was successful draw ``G``
using one of

.. nbplot::

    >>> G = nx.petersen_graph()
    >>> subax1 = plt.subplot(121)
    >>> nx.draw(G, with_labels=True, font_weight='bold')
    >>> subax2 = plt.subplot(122)
    >>> nx.draw_shell(G, nlist=[range(5, 10), range(5)], with_labels=True, font_weight='bold')

when drawing to an interactive display.  Note that you may need to issue a
Matplotlib

>>> plt.show()  # doctest: +SKIP

command if you are not using matplotlib in interactive mode.

.. nbplot::

    >>> options = {
    ...     'node_color': 'black',
    ...     'node_size': 100,
    ...     'width': 3,
    ... }
    >>> subax1 = plt.subplot(221)
    >>> nx.draw_random(G, **options)
    >>> subax2 = plt.subplot(222)
    >>> nx.draw_circular(G, **options)
    >>> subax3 = plt.subplot(223)
    >>> nx.draw_spectral(G, **options)
    >>> subax4 = plt.subplot(224)
    >>> nx.draw_shell(G, nlist=[range(5,10), range(5)], **options)

You can find additional options via :func:`~drawing.nx_pylab.draw_networkx` and
layouts via the :mod:`layout module<networkx.drawing.layout>`.
You can use multiple shells with :func:`~drawing.nx_pylab.draw_shell`.

.. nbplot::

    >>> G = nx.dodecahedral_graph()
    >>> shells = [[2, 3, 4, 5, 6], [8, 1, 0, 19, 18, 17, 16, 15, 14, 7], [9, 10, 11, 12, 13]]
    >>> nx.draw_shell(G, nlist=shells, **options)

To save drawings to a file, use, for example

>>> nx.draw(G)
>>> plt.savefig("path.png")

This function writes to the file ``path.png`` in the local directory. If Graphviz and
PyGraphviz or pydot, are available on your system, you can also use
`networkx.drawing.nx_agraph.graphviz_layout` or
`networkx.drawing.nx_pydot.graphviz_layout` to get the node positions, or write
the graph in dot format for further processing.

>>> from networkx.drawing.nx_pydot import write_dot
>>> pos = nx.nx_agraph.graphviz_layout(G)
>>> nx.draw(G, pos=pos)
>>> write_dot(G, 'file.dot')

See :doc:`/reference/drawing` for additional details.

.. code-links::
.. _contents:

Software for Complex Networks
=============================

.. only:: html

    :Release: |version|
    :Date: |today|

NetworkX is a Python package for the creation, manipulation, and study
of the structure, dynamics, and functions of complex networks.
It provides:

-  tools for the study of the structure and
   dynamics of social, biological, and infrastructure networks;
-  a standard programming interface and graph implementation that is suitable
   for many applications;
-  a rapid development environment for collaborative, multidisciplinary
   projects;
-  an interface to existing numerical algorithms and code written in C,
   C++, and FORTRAN; and
-  the ability to painlessly work with large nonstandard data sets.

With NetworkX you can load and store networks in standard and nonstandard data
formats, generate many types of random and classic networks, analyze network
structure, build network models, design new network algorithms, draw networks,
and much more.

Citing
------

To cite NetworkX please use the following publication:

Aric A. Hagberg, Daniel A. Schult and Pieter J. Swart,
`"Exploring network structure, dynamics, and function using NetworkX"
<http://conference.scipy.org/proceedings/SciPy2008/paper_2/>`_,
in
`Proceedings of the 7th Python in Science Conference (SciPy2008)
<http://conference.scipy.org/proceedings/SciPy2008/index.html>`_, Gäel
Varoquaux, Travis Vaught, and Jarrod Millman (Eds), (Pasadena, CA
USA), pp. 11--15, Aug 2008

.. only:: html

   `PDF <http://conference.scipy.org/proceedings/SciPy2008/paper_2/full_text.pdf>`_
   `BibTeX <http://conference.scipy.org/proceedings/SciPy2008/paper_2/reference.bib>`_

Audience
--------

The audience for NetworkX includes mathematicians, physicists, biologists,
computer scientists, and social scientists. Good reviews of the science of
complex networks are presented in Albert and Barabási [BA02]_, Newman
[Newman03]_, and Dorogovtsev and Mendes [DM03]_. See also the classic texts
[Bollobas01]_, [Diestel97]_ and [West01]_ for graph theoretic results and
terminology. For basic graph algorithms, we recommend the texts of Sedgewick
(e.g., [Sedgewick01]_ and [Sedgewick02]_) and the survey of Brandes and
Erlebach [BE05]_.

Python
------

Python is a powerful programming language that allows simple and flexible
representations of networks as well as clear and concise expressions of network
algorithms.  Python has a vibrant and growing ecosystem of packages that
NetworkX uses to provide more features such as numerical linear algebra and
drawing.  In order to make the most out of NetworkX you will want to know how
to write basic programs in Python.  Among the many guides to Python, we
recommend the `Python documentation <https://docs.python.org/3/>`_ and the text
by Alex Martelli [Martelli03]_.

License
-------

.. include:: ../LICENSE.txt

Bibliography
------------

.. [BA02] R. Albert and A.-L. Barabási, "Statistical mechanics of complex
   networks", Reviews of Modern Physics, 74, pp. 47-97, 2002.
   https://arxiv.org/abs/cond-mat/0106096

.. [Bollobas01] B. Bollobás, "Random Graphs", Second Edition,
   Cambridge University Press, 2001.

.. [BE05] U. Brandes and T. Erlebach, "Network Analysis:
   Methodological Foundations", Lecture Notes in Computer Science,
   Volume 3418, Springer-Verlag, 2005.

.. [Diestel97] R. Diestel, "Graph Theory", Springer-Verlag, 1997.
   http://diestel-graph-theory.com/index.html

.. [DM03] S.N. Dorogovtsev and J.F.F. Mendes, "Evolution of Networks",
   Oxford University Press, 2003.

.. [Martelli03]  A. Martelli, "Python in a Nutshell", O'Reilly Media
   Inc, 2003.

.. [Newman03] M.E.J. Newman, "The Structure and Function of Complex
   Networks", SIAM Review, 45, pp. 167-256, 2003.
   http://epubs.siam.org/doi/abs/10.1137/S003614450342480

.. [Sedgewick02] R. Sedgewick, "Algorithms in C: Parts 1-4:
   Fundamentals, Data Structure, Sorting, Searching", Addison Wesley
   Professional, 3rd ed., 2002.

.. [Sedgewick01] R. Sedgewick, "Algorithms in C, Part 5: Graph Algorithms",
   Addison Wesley Professional, 3rd ed., 2001.

.. [West01] D. B. West, "Introduction to Graph Theory", Prentice Hall,
    2nd ed., 2001.

.. toctree::
   :maxdepth: 1
   :hidden:

   install
   tutorial
   reference/index
   release/index
   developer/index
   auto_examples/index
.. include:: ../INSTALL.rst
{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}
   .. automethod:: __init__

   {% if methods %}
   .. rubric:: {{ _('Methods') }}

   .. autosummary::
      :toctree: generated/

   {% for item in methods %}
       {% if item != "__init__" %}
          ~{{ name }}.{{ item }}
       {% endif %}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Attributes') }}

   .. autosummary::
      :toctree: generated/

   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
.. _mission_and_values:

==================
Mission and Values
==================

Our mission
-----------

NetworkX aims to be the reference library for network science algorithms in
Python. We accomplish this by:

- **being easy to use and install**. We are careful in taking on new
  dependencies, and sometimes cull existing ones, or make them optional. All
  functions in our API have thorough docstrings clarifying expected inputs and
  outputs.
- **providing a consistent API**. Conceptually identical arguments have the
  same name and position in a function signature.
- **ensuring correctness**. Test coverage is close to 100% and code is reviewed by
  at least two core developers before being included in the library.
- **caring for users’ data**. We have a functional API and don't modify
  input data unless explicitly directed to do so.
- **promoting education in network science**, with extensive pedagogical
  documentation.

Our values
----------

- We are inclusive (:ref:`code_of_conduct`). We welcome and mentor newcomers who are
  making their first contribution.
- We are open source and community-driven (:ref:`governance`).
- We focus on graph data structures and algorithms for network science applications.
- We prefer pure Python implementations using native data structures
  (especially dicts) due to their consistent, intuitive interface and amazing
  performance capabilities. We include interfaces to other data structures,
  especially NumPy arrays and SciPy sparse matrices for algorithms that more
  naturally use arrays and matrices or where time or space requirements are
  significantly lower. Sometimes we provide two algorithms for the same result,
  one using each data structure, when pedagogy or space/time trade-offs justify
  such multiplicity.
- We value simple, readable implementations over getting every last ounce of
  performance. Readable code that is easy to understand, for newcomers and
  maintainers alike, makes it easier to contribute new code as well as prevent
  bugs. This means that we will prefer a 20% slowdown if it reduces lines of
  code two-fold, for example.
- We value education and documentation. All functions should have `NumPy-style
  docstrings <https://numpy.org/doc/stable/docs/howto_document.html>`,
  preferably with examples, as well as gallery examples that showcase how that
  function is used in a scientific application.

Acknowledgments
---------------

This document is modified from the `scikit-image` mission and values document.
.. include:: ../../CONTRIBUTING.rst
About Us
========

NetworkX was originally written by Aric Hagberg, Dan Schult, and Pieter Swart,
and has been developed with the help of many others. Thanks to everyone who has
improved NetworkX by contributing code, bug reports (and fixes), documentation,
and input on design, features, and the future of NetworkX.

.. include:: team.rst 

Contributors
------------

If you are a NetworkX contributor, please feel free to
open an `issue <https://github.com/networkx/networkx/issues/new>`_ or
submit a `pull request <https://github.com/networkx/networkx/compare/>`_
to add your name to the bottom of the list.

- Aric Hagberg, GitHub: `hagberg <https://github.com/hagberg>`_
- Dan Schult, GitHub: `dschult <https://github.com/dschult>`_
- Pieter Swart
- Katy Bold
- Hernan Rozenfeld
- Brendt Wohlberg
- Jim Bagrow
- Holly Johnsen
- Arnar Flatberg
- Chris Myers
- Joel Miller
- Keith Briggs
- Ignacio Rozada
- Phillipp Pagel
- Sverre Sundsdal
- Ross M. Richardson
- Eben Kenah
- Sasha Gutfriend
- Udi Weinsberg
- Matteo Dell'Amico
- Andrew Conway
- Raf Guns
- Salim Fadhley
- Fabrice Desclaux
- Arpad Horvath
- Minh Van Nguyen
- Willem Ligtenberg
- Loïc Séguin-C.
- Paul McGuire
- Jesus Cerquides
- Ben Edwards
- Jon Olav Vik
- Hugh Brown
- Ben Reilly
- Leo Lopes
- Jordi Torrents, GitHub: `jtorrents <https://github.com/jtorrents>`_
- Dheeraj M R
- Franck Kalala
- Simon Knight
- Conrad Lee
- Sérgio Nery Simões
- Robert King
- Nick Mancuso
- Brian Cloteaux
- Alejandro Weinstein
- Dustin Smith
- Mathieu Larose
- Romain Fontugne
- Vincent Gauthier
- chebee7i, GitHub: `chebee7i <https://github.com/chebee7i>`_
- Jeffrey Finkelstein
- Jean-Gabriel Young, Github: `jg-you <https://github.com/jgyou>`_
- Andrey Paramonov, http://aparamon.msk.ru
- Mridul Seth, GitHub: `MridulS <https://github.com/MridulS>`_
- Thodoris Sotiropoulos, GitHub: `theosotr <https://github.com/theosotr>`_
- Konstantinos Karakatsanis, GitHub: `k-karakatsanis <https://github.com/k-karakatsanis>`_
- Ryan Nelson, GitHub: `rnelsonchem <https://github.com/rnelsonchem>`_
- Niels van Adrichem, GitHub: `NvanAdrichem <https://github.com/NvanAdrichem>`_
- Michael E. Rose, GitHub: `Michael-E-Rose <https://github.com/Michael-E-Rose>`_
- Jarrod Millman, GitHub: `jarrodmillman <https://github.com/jarrodmillman>`_
- Andre Weltsch
- Lewis Robbins
- Mads Jensen, Github: `atombrella <https://github.com/atombrella>`_
- Edward L. Platt, `elplatt <https://elplatt.com>`_
- James Owen, Github: `leamingrad <https://github.com/leamingrad>`_
- Robert Gmyr, Github: `gmyr <https://github.com/gmyr>`_
- Mike Trenfield
- Jon Crall, Github: `Erotemic <https://github.com/Erotemic>`_
- Issa Moradnejad, Github `<https://github.com/Moradnejad>`_, LinkedIn `<https://linkedin.com/in/moradnejad/>`_
- Brian Kiefer, Github: `bkief <https://github.com/bkief>`_
- Julien Klaus
- Peter C. Kroon, Github: `pckroon <https://github.com/pckroon>`_
- Weisheng Si, Github: `ws4u <https://github.com/ws4u>`_
- Haakon H. Rød, Gitlab: `haakonhr <https://gitlab.com/haakonhr>`_, `<https://haakonhr.gitlab.io>`_
- Efraim Rodrigues, GitHub `<https://github.com/efraimrodrigues>`_, LinkedIn `<https://linkedin.com/in/efraim-rodrigues/>`_
- Erwan Le Merrer
- Søren Fuglede Jørgensen, GitHub: `fuglede <https://github.com/fuglede>`_
- Salim BELHADDAD, LinkedIn `<https://www.linkedin.com/in/salymdotme/>`_
- Jangwon Yie, GitHub `<https://github.com/jangwon-yie>`_, LinkedIn `<https://www.linkedin.com/in/jangwon-yie-a7960065/>`_
- ysitu
- Tomas Gavenciak
- Luca Baldesi
- Yuto Yamaguchi
- James Clough
- Minas Gjoka
- Drew Conway
- Alex Levenson
- Haochen Wu
- Erwan  Le Merrer
- Alex Roper
- P C Kroon
- Christopher Ellison
- D. Eppstein
- Federico Rosato
- Aitor Almeida
- Ferran Parés
- Christian Olsson
- Fredrik Erlandsson
- Nanda H Krishna
- Nicholas Mancuso
- Fred Morstatter
- Ollie Glass
- Rodrigo Dorantes-Gilardi
- Pranay Kanwar
- Balint Tillman
- Diederik van Liere
- Ferdinando Papale
- Miguel Sozinho Ramalho
- Brandon Liu
- Nima Mohammadi
- Jason Grout
- Jan Aagaard Meier
- Henrik Haugbølle
- Piotr Brodka
- Sasha Gutfraind
- Alessandro Luongo
- Huston Hedinger
- Oleguer Sagarra
- Kazimierz Wojciechowski, GitHub `<https://github.com/kazimierz-256>`_, LinkedIn `<https://linkedin.com/in/wojciechowski-kazimierz/>`_
- Gaetano Pietro Paolo Carpinato, Github `<https://github.com/Carghaez>`_, LinkedIn `<https://linkedin.com/in/gaetanocarpinato/>`_
- Arun Nampally, GitHub `<https://github.com/arunwise>`_, LinkedIn `<https://www.linkedin.com/in/arun-nampally-b57845b7/>`_
- Ryan Duve
- Shashi Prakash Tripathi, Github `<https://github.com/itsshavar>`_,LinkedIn `<https://www.linkedin.com/in/itsshashitripathi/>`_
- Danny Niquette
- James Trimble, Github: `jamestrimble <https://github.com/jamestrimble>`_
- Matthias Bruhns, Github `<https://github.com/mbruhns>`_ 
- Philip Boalch
- Matt Schwennesen, Github: `mjschwenne <https://github.com/mjschwenne>`_

A supplementary (but still incomplete) list of contributors is given by the
list of names that have commits in ``networkx``'s
`git <http://git-scm.com>`_ repository. This can be obtained via::

    git log --raw | grep "^Author: " | sort | uniq

A historical, partial listing of contributors and their contributions to some
of the earlier versions of NetworkX can be found
`here <https://github.com/networkx/networkx/blob/886e790437bcf30e9f58368829d483efef7a2acc/doc/source/reference/credits_old.rst>`_.


Support
-------

NetworkX acknowledges support from the following research groups:

- `Center for Nonlinear Studies <http://cnls.lanl.gov>`_, Los Alamos National
  Laboratory, PI: Aric Hagberg

- `Open Source Programs Office <https://developers.google.com/open-source/>`_,
  Google

- `Complexity Sciences Center <http://csc.ucdavis.edu/>`_, Department of
  Physics, University of California-Davis, PI: James P. Crutchfield

- `Center for Complexity and Collective Computation <http://c4.discovery.wisc.edu>`_,
  Wisconsin Institute for Discovery, University of Wisconsin-Madison,
  PIs: Jessica C. Flack and David C. Krakauer

NetworkX acknowledges the following financial support:

- Google Summer of Code via Python Software Foundation

- U.S. Army Research Office grant W911NF-12-1-0288

- DARPA Physical Intelligence Subcontract No. 9060-000709

- NSF Grant No. PHY-0748828

- John Templeton Foundation through a grant to the Santa Fe Institute to
  study complexity

- U.S. Army Research Laboratory and the U.S. Army Research Office under
  contract number W911NF-13-1-0340
.. _contributing_faq:

New Contributor FAQ
*******************

A collection of frequently-asked questions by newcomers to
open-source development and first-time contributors to NetworkX.

Q: I'm new to open source and would like to contribute to NetworkX. How do I get started?
-----------------------------------------------------------------------------------------

To contribute to NetworkX, you will need three things:

  1. The source code
  2. A development environment
  3. An idea of what you'd like to contribute

Steps 1 & 2 are covered extensively in :ref:`Development Workflow <dev_workflow>`.
There is no generic answer for step 3. There are many ways that NetworkX can
be improved, from adding new algorithms, improving existing algorithms,
improving the test suite (e.g. increasing test coverage), and improving the
documentation.
The "best" way to find a place to start is to follow your own personal
interests!
That said, a few places to check for ideas on where to get started:

 - `The issue tracker <https://github.com/networkx/networkx/issues>`_ lists
   known bugs and feature requests. Of particular interest for first-time
   contributors are issues that have been tagged with the `Good First Issue`_
   or `Sprint`_ labels.
 - The `Algorithms discussion`_ includes a listing of algorithms that users
   would like to have but that are not yet included in NetworkX.

.. _Good First Issue: https://github.com/networkx/networkx/issues?q=is%3Aopen+is%3Aissue+label%3A%22Good+First+Issue%22

.. _Sprint: https://github.com/networkx/networkx/issues?q=is%3Aopen+is%3Aissue+label%3ASprint

.. _Algorithms discussion: https://github.com/networkx/networkx/discussions/categories/algorithms

Q: I've found an issue I'm interested in, can I have it assigned to me?
-----------------------------------------------------------------------

NetworkX doesn't typically assign issues to contributors. If you find an issue
or feature request on the issue tracker that you'd like to work on, you should
first check the issue thread to see if there are any linked pull requests.
If not, then feel free to open a new PR to address the issue - no need
to ask for permission - and don't forget to reference the issue number in the PR
comments so that others know you are now working on it!

Q: How do I contribute an example to the Gallery?
-------------------------------------------------

The example gallery is great place to contribute, particularly if you have an
interesting application or visualization that uses NetworkX.
The gallery is generated using :doc:`sphinx-gallery <sphinx-gallery:index>`
from Python scripts stored in the ``examples/`` directory.

For instance, let's say I'd like to contribute an example of visualizing a
`complete graph <networkx.generators.classic.complete_graph>` using a
`circular layout <networkx.drawing.layout.circular_layout>`.
Assuming you have already followed the procedure for
:ref:`setting up a development environment <dev_workflow>`, start by
creating a new branch:

.. code-block:: bash

   git checkout -b complete-graph-circular-layout-example

.. note:: It's generally a good idea to give your branch a descriptive name so
   that it's easy to remember what you are working on.

Now you can begin work on your example. Sticking with the circular layout idea,
you might create a file in ``examples/drawing`` called ``plot_circular_layout.py``
with the following contents::

   import networkx as nx
   import matplotlib.pyplot as plt

   G = nx.complete_graph(10)  # A complete graph with 10 nodes
   nx.draw_networkx(G, pos=nx.circular_layout(G))

.. note:: It may not be clear where exactly an example belongs. Our circular
   layout example is very simple, so perhaps it belongs in ``examples/basic``.
   It would also make sense for it to be in ``examples/drawing`` since it deals
   with visualization. Don't worry if you're not sure: questions like this will
   be resolved during the review process.

At this point, your contribution is ready to be reviewed. You can make the
changes on your ``complete-graph-circular-layout-example`` branch visible to
other NetworkX developers by
`creating a pull request`__. 

.. _PR: https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request

__ PR_

.. seealso:: The :ref:`developer guide <dev_workflow>` has more details on
   creating pull requests.

Q: I want to work on a specific function. How do I find it in the source code?
------------------------------------------------------------------------------

Assuming you have followed the instructions for
:ref:`setting up the development workflow <dev_workflow>`, there are several
ways of determining where the in the **source code** a particular function or
class is defined.

For example, let's say you are interested in making a change to the
`~networkx.drawing.layout.kamada_kawai_layout` function, so you need to know
where it is defined. In an IPython terminal, you can use ``?`` --- the source file is
listed in the ``File:`` field:

.. code-block:: ipython

   In [1]: import networkx as nx
   In [2]: nx.kamada_kawai_layout?

.. code-block:: text

   Signature: <clipped for brevity>
   Docstring: <clipped for brevity>
   File: ~/networkx/networkx/drawing/layout.py
   Type: function

Command line utilities like ``grep`` or ``git grep`` are also very useful.
For example, from the NetworkX source directory:

.. code-block:: bash

   $ grep -r "def kamada_kawai_layout" .
   ./networkx/drawing/layout.py:def kamada_kawai_layout(

Q: What is the policy for deciding whether to include a new algorithm?
----------------------------------------------------------------------

There is no official policy setting explicit inclusion criteria for new
algorithms in NetworkX. New algorithms are more likely to be included if they
have been published and are cited by others. More important than number of
citations is how well proposed additions fit the project :ref:`mission_and_values`.
=======
Roadmap
=======

The roadmap is intended for larger, fundamental changes to
the project that are likely to take months or years of developer time.
Smaller-scoped items will continue to be tracked on our issue tracker.

The scope of these improvements means that these changes may be
controversial, are likely to involve significant discussion
among the core development team, and may require the creation
of one or more :ref:`nxep`.

Installation
------------

We aim to make NetworkX as easy to install as possible.
Some of our dependencies (e.g., graphviz and gdal) can be tricky to install.
Other of our dependencies are easy to install on the CPython platform, but
may be more involved on other platforms such as PyPy.
Addressing these installation issues may involve working with the external projects.

Sustainability
--------------

We aim to reduce barriers to contribution, speed up pull request (PR) review,
onboard new maintainers, and attract new developers to ensure the long-term
sustainability of NetworkX.

This includes:

- improving continuous integration
- making code base more approachable
- creating new pathways beyond volunteer effort
- growing maintainers and leadership
- increasing diversity of developer community

Performance
-----------

Speed improvements, lower memory usage, and the ability to parallelize
algorithms are beneficial to most science domains and use cases.

A first step may include implementing a benchmarking system using something
like airspeed velocity (https://asv.readthedocs.io/en/stable/).
It may also include review existing comparisons between NetworkX
and other packages.

Individual functions can be optimized for performance and memory use.
We are also interested in exploring new technologies to accelerate
code and reduce memory use.  Before adopting any new technologies
we will need to careful consider its impact on code readability
and difficulty of building and installating NetworkX.
For more information, see our :ref:`mission_and_values`.

Many functions can be trivially parallelized.
But, we need to decide on an API and perhaps implement some
helper code to make it consistent.

Documentation
-------------

We’d like to improve the content, structure, and presentation of the NetworkX
documentation. Some specific goals include:

- longer gallery examples
- domain-specific documentation (NetworkX for Geneticists,
  NetworkX for Neuroscientists, etc.)
- examples of how to use NetworkX with other packages

Linear Algebra
--------------

We would like to improve our linear algebra based algorithms.
The code is old and needs review and refactoring.
This would include investigating SciPy's csgraph.
It would also include deciding how to handle algorithms that
have multiple implementations (e.g., some algorithms are implemented in Python,
NumPy, and SciPy).

NumPy has split its API from its execution engine with ``__array_function__`` and
``__array_ufunc__``. This will enable parts of NumPy to accept distributed arrays
(e.g. dask.array.Array) and GPU arrays (e.g. cupy.ndarray) that implement the
ndarray interface. At the moment it is not yet clear which algorithms will work
out of the box, and if there are significant performance gains when they do.

Interoperability
----------------

We'd like to improve interoperability with the rest of the scientific Python
ecosystem.
This includes projects we depend on (e.g., NumPy, SciPy, Pandas, Matplotlib)
as well as ones we don't (e.g., Geopandas).

For example, we would also like to be able to seamlessly exchange graphs with
other network analysis software.
Another way to integrate with other scientific python ecosystem tools is to
take on features from the other tools that are useful. And we should develop
tools to ease use of NetworkX from within these other tools.
Additional examples of interoperability improvements may include providing a more
pandas-like interface for the ```__getitem__``` dunder function of node and
edge views (:ref:`NXEP2`).
Also developing a universal method to represent a graph as a single sequence of
```nodes_and_edges``` objects that allow attribute dicts, nodes and edges as
`discussed for graph generators
<https://github.com/networkx/networkx/issues/3036>`_.

Visualization
-------------

Visualization is not a focus on NetworkX, but it is a major feature for
many users.
We need to enhance the drawing tools for NetworkX.
Mentored Projects
==================

This page maintains a list of mentored project ideas that contributors can work
on if they are interested in contributing to the NetworkX project. Feel free to
suggest any other idea if you are interested on the
`NetworkX GitHub discussions page <https://github.com/networkx/networkx/discussions>`__

These ideas can be used as projects for Google Summer of Code, Outreachy,
NumFOCUS Small Development Grants and university course/project credits (if
your university allows contribution to open source for credit).


Pedagogical Interactive Notebooks for Algorithms Implemented in NetworkX
------------------------------------------------------------------------

- Abstract: NetworkX has a :ref:`wide variety of algorithms <Algorithms>`
  implemented. Even though the algorithms are well documented, explanations of
  the ideas behind the algorithms are often missing and we would like to
  collect these, write Jupyter notebooks to elucidate these ideas and explore
  the algorithms experimentally, and publish the notebooks at
  https://github.com/networkx/notebooks. The goal is to gives readers a
  deeper outlook behind standard network science and graph theory algorithms
  and encourage them to delve further into the topic.

- Recommended Skills: Python, Jupyter notebooks, graph algorithms.

- Expected Outcome: A collection of Interactive Jupyter notebooks which
  explain and explore network algorithms to readers and users of NetworkX.
  For example, see this notebook on
  :doc:`Random Geometric Graphs <content/generators/geometric>`

- Complexity: Depending on the algorithms you are interested to work on.

- Interested Mentors: `@dschult <https://github.com/dschult/>`__,
  `@MridulS <https://github.com/MridulS/>`__,
  `@rossbar <https://github.com/rossbar/>`__

Implement the VF2++ Graph Isomorphism Algorithm
-----------------------------------------------

- Abstract: The `Graph Isomorphism Problem`_ is a famous difficult network problem at
  the boundary between P and NP-Complete. The VF2 algorithm is included with NetworkX
  in a recursive formulation. There is an improved version of this algorithm called
  `VF2++`_ which we intend to implement. We have early attempts at a nonrecursive version
  of the main algorithm that also address subgraph isomorphism and subgraph monomorphism.
  This project involves fully implementing them and extending to directed and multigraph
  settings.

- Recommended Skills: Python, graph algorithms

- Expected Outcome: A new set of functions in NetworkX that implement the VF2++
  algorithm for all problem and graph types in a nonrecursive manner.

- Complexity: Moderate

- Interested Mentors: `@dschult <https://github.com/dschult/>`__,
  `@MridulS <https://github.com/MridulS/>`__, `@boothby <https://github.com/boothby/>`__,

.. _`Graph Isomorphism Problem`: https://en.wikipedia.org/wiki/Graph_isomorphism_problem
.. _VF2++: https://doi.org/10.1016/j.dam.2018.02.018

Completed Projects
==================

- `Louvain community detection algorithm`_ 
    - Program: Google Summer of Code 2021
    - Contributor: `@z3y50n <https://github.com/z3y50n/>`__
    - Link to Proposal:  `GSoC 2021: Community Detection Algorithms <https://github.com/networkx/archive/blob/main/proposals-gsoc/GSoC-2021-Community-Detection-Algorithms.pdf>`__ 

- `Asadpour algorithm for directed travelling salesman problem`_
    - Program: Google Summer of Code 2021
    - Contributor: `@mjschwenne <https://github.com/mjschwenne/>`__
    - Link to Proposal:  `GSoC 2021: Asadpour algorithm <https://github.com/networkx/archive/blob/main/proposals-gsoc/GSoC-2021-Asadpour-Asymmetric-Traveling%20Salesman-Problem.pdf>`__ 

- Pedagogical notebook: `Directed acyclic graphs and topological sort`_
    - Program: Google Summer of Code 2021
    - Contributor:  `@vdshk <https://github.com/vdshk>`__

- Pedagogical notebooks: `Graph assortativity`_ & `Network flow analysis and Dinitz algorithm`_
    - Program: Google Summer of Code 2021
    - Contributor: `@harshal-dupare <https://github.com/harshal-dupare/>`__

- Add On system for NetworkX: `NetworkX-Metis`_
    - Program: Google Summer of Code 2015
    - Contributor: `@OrkoHunter <https://github.com/OrkoHunter/>`__
    - Link to Proposal:  `GSoC 2015: Add On System for NetworkX <https://github.com/networkx/archive/blob/main/proposals-gsoc/GSoC-2015-Add-on-system-for-NetworkX.md>`__

- `NetworkX 2.0 API`_
    - Program: Google Summer of Code 2015
    - Contributor: `@MridulS <https://github.com/MridulS/>`__
    - Link to Proposal: `GSoC 2015: NetworkX 2.0 API <https://github.com/networkx/archive/blob/main/proposals-gsoc/GSoC-2015-NetworkX-2.0-api.md>`__

.. _`Louvain community detection algorithm`: https://github.com/networkx/networkx/pull/4929
.. _`Asadpour algorithm for directed travelling salesman problem`: https://github.com/networkx/networkx/pull/4740
.. _`Directed acyclic graphs and topological sort`: https://github.com/networkx/nx-guides/pull/44
.. _`Graph assortativity`: https://github.com/networkx/nx-guides/pull/42
.. _`Network flow analysis and Dinitz algorithm`: https://github.com/networkx/nx-guides/pull/46
.. _`NetworkX-Metis`: https://github.com/networkx/networkx-metis
.. _`NetworkX 2.0 API`: https://networkx.org/documentation/latest/release/migration_guide_from_1.x_to_2.0.html

..
   Project Idea Template
   ---------------------
   
   - Abstract:
   
   - Recommended Skills:
   
   - Expected Outcome:
   
   - Complexity;
   
   - Interested Mentors:
   
.. _developer:

Developer
*********

.. only:: html

    :Release: |version|
    :Date: |today|

.. toctree::
   :maxdepth: 2

   about_us
   code_of_conduct
   values
   contribute
   projects
   new_contributor_faq
   core_developer
   release
   deprecations
   roadmap
   nxeps/index
.. include:: ../../CODE_OF_CONDUCT.rst
.. _core_dev:

Core Developer Guide
====================

As a core developer, you should continue making pull requests
in accordance with the :ref:`contributor_guide`.
You are responsible for shepherding other contributors through the review process.
You should be familiar with our :ref:`mission_and_values`.
You also have the ability to merge or approve other contributors' pull requests.
Much like nuclear launch keys, it is a shared power: you must merge *only after*
another core developer has approved the pull request, *and* after you yourself have carefully
reviewed it.  (See `Reviewing`_ and especially `Merge Only Changes You
Understand`_ below.) To ensure a clean git history, use GitHub's
`Squash and Merge <https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/merging-a-pull-request#merging-a-pull-request-on-github>`__
feature to merge, unless you have a good reason not to do so.

Reviewing
---------

How to Conduct A Good Review
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Always* be kind to contributors. Nearly all of NetworkX is
volunteer work, for which we are tremendously grateful. Provide
constructive criticism on ideas and implementations, and remind
yourself of how it felt when your own work was being evaluated as a
novice.

NetworkX strongly values mentorship in code review.  New users
often need more handholding, having little to no git
experience. Repeat yourself liberally, and, if you don’t recognize a
contributor, point them to our development guide, or other GitHub
workflow tutorials around the web. Do not assume that they know how
GitHub works (e.g., many don't realize that adding a commit
automatically updates a pull request). Gentle, polite, kind
encouragement can make the difference between a new core developer and
an abandoned pull request.

When reviewing, focus on the following:

1. **API:** The API is what users see when they first use
   NetworkX. APIs are difficult to change once released, so
   should be simple, `functional
   <https://en.wikipedia.org/wiki/Functional_programming>`__ (i.e. not
   carry state), consistent with other parts of the library, and
   should avoid modifying input variables.  Please familiarize
   yourself with the project's :ref:`deprecation_policy`.

2. **Documentation:** Any new feature should have a gallery
   example that not only illustrates but explains it.

3. **The algorithm:** You should understand the code being modified or
   added before approving it.  (See `Merge Only Changes You
   Understand`_ below.) Implementations should do what they claim,
   and be simple, readable, and efficient.

4. **Tests:** All contributions to the library *must* be tested, and
   each added line of code should be covered by at least one test. Good
   tests not only execute the code, but explores corner cases.  It is tempting
   not to review tests, but please do so.

Other changes may be *nitpicky*: spelling mistakes, formatting,
etc. Do not ask contributors to make these changes, and instead
make the changes by `pushing to their branch
<https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/committing-changes-to-a-pull-request-branch-created-from-a-fork>`__,
or using GitHub’s `suggestion
<https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/commenting-on-a-pull-request>`__
`feature
<https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/incorporating-feedback-in-your-pull-request>`__.
(The latter is preferred because it gives the contributor a choice in
whether to accept the changes.)

Our default merge policy is to squash all PR commits into a single
commit. Users who wish to bring the latest changes from ``main``
into their branch should be advised to merge, not to rebase.  Even
when merge conflicts arise, don’t ask for a rebase unless you know
that a contributor is experienced with git. Instead, rebase the branch
yourself, force-push to their branch, and advise the contributor on
how to force-pull.  If the contributor is no longer active, you may
take over their branch by submitting a new pull request and closing
the original. In doing so, ensure you communicate that you are not
throwing the contributor's work away!  You should use GitHub's
``Co-authored-by:`` keyword for commit messages to credit the
original contributor.


Please add a note to a pull request after you push new changes; GitHub
may not send out notifications for these.

Merge Only Changes You Understand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Long-term maintainability* is an important concern.  Code doesn't
merely have to *work*, but should be *understood* by multiple core
developers.  Changes will have to be made in the future, and the
original contributor may have moved on.

Therefore, *do not merge a code change unless you understand it*. Ask
for help freely: we have a long history of consulting community
members, or even external developers, for added insight where needed,
and see this as a great learning opportunity.

While we collectively "own" any patches (and bugs!) that become part
of the code base, you are vouching for changes you merge.  Please take
that responsibility seriously.

Closing issues and pull requests
--------------------------------

Sometimes, an issue must be closed that was not fully resolved. This can be
for a number of reasons:

- the person behind the original post has not responded to calls for
  clarification, and none of the core developers have been able to reproduce
  their issue;
- fixing the issue is difficult, and it is deemed too niche a use case to
  devote sustained effort or prioritize over other issues; or
- the use case or feature request is something that core developers feel
  does not belong in NetworkX,

among others. Similarly, pull requests sometimes need to be closed without
merging, because:

- the pull request implements a niche feature that we consider not worth the
  added maintenance burden;
- the pull request implements a useful feature, but requires significant
  effort to bring up to NetworkX's standards, and the original
  contributor has moved on, and no other developer can be found to make the
  necessary changes; or
- the pull request makes changes that do not align with our values, such as
  increasing the code complexity of a function significantly to implement a
  marginal speedup,

among others.

All these may be valid reasons for closing, but we must be wary not to alienate
contributors by closing an issue or pull request without an explanation. When
closing, your message should:

- explain clearly how the decision was made to close. This is particularly
  important when the decision was made in a community meeting, which does not
  have as visible a record as the comments thread on the issue itself;
- thank the contributor(s) for their work; and
- provide a clear path for the contributor or anyone else to appeal the
  decision.

These points help ensure that all contributors feel welcome and empowered to
keep contributing, regardless of the outcome of past contributions.

Further resources
-----------------

As a core member, you should be familiar with community and developer
resources such as:

-  Our :ref:`contributor_guide`
-  Our :ref:`code_of_conduct`
-  `PEP8 <https://www.python.org/dev/peps/pep-0008/>`__ for Python style
-  `PEP257 <https://www.python.org/dev/peps/pep-0257/>`__ and the `NumPy
   documentation
   guide <https://numpy.org/doc/stable/docs/howto_document.html>`__
   for docstrings. (NumPy docstrings are a superset of PEP257. You
   should read both.)
-  The NetworkX `tag on
   StackOverflow <https://stackoverflow.com/questions/tagged/networkx>`__
-  Our `mailing
   list <http://groups.google.com/group/networkx-discuss/>`__

You are not required to monitor all of the social resources.
Release Process
===============

- Update the release notes:

  1. Review and cleanup ``doc/release/release_dev.rst``,

  2. Fix code in documentation by running
     ``cd doc && make doctest``.

  3. Make a list of merges and contributors by running
     ``doc/release/contribs.py <tag of previous release>``.

  4. Paste this list at the end of the ``release_dev.rst``. Scan the PR titles
     for highlights, deprecations, and API changes, and mention these in the
     relevant sections of the notes.

  5. Rename to ``doc/release/release_<major>.<minor>.rst``.

  6. Copy ``doc/release/release_template.rst`` to
     ``doc/release/release_dev.rst`` for the next release.

  7. Add ``release_<major>.<minor>`` to ``doc/release/index.rst``.

- Delete developer banner on docs::

   git rm doc/_templates/layout.html

- Update ``__version__`` in ``networkx/__init__.py``.

- Commit changes::

   git add networkx/__init__.py
   git commit -m "Designate X.X release"

- Add the version number as a tag in git::

   git tag -s [-u <key-id>] networkx-<major>.<minor> -m 'signed <major>.<minor> tag'

  (If you do not have a gpg key, use -m instead; it is important for
  Debian packaging that the tags are annotated)

- Push the new meta-data to github::

   git push --tags upstream main

  (where ``upstream`` is the name of the
   ``github.com:networkx/networkx`` repository.)

- Review the github release page::

   https://github.com/networkx/networkx/releases

- Pin badges in ``README.rst``::

  - https://github.com/networkx/networkx/workflows/test/badge.svg?tag=networkx-<major>.<minor>
  - https://github.com/networkx/networkx/actions?query=branch%3Anetworkx-<major>.<minor>

- Publish on PyPi::

   git clean -fxd
   pip install -r requirements/release.txt
   python setup.py sdist bdist_wheel
   twine upload -s dist/*

- Unpin badges in ``README.rst``::

   git restore README.rst 

- Update documentation on the web:
  The documentation is kept in a separate repo: networkx/documentation

  - Wait for the CI service to deploy to GitHub Pages
  - Sync your branch with the remote repo: ``git pull``.
  - Copy the documentation built by the CI service.
    Assuming you are at the top-level of the ``documentation`` repo::

      # FIXME - use eol_banner.html
      cp -a latest networkx-<major>.<minor>
      ln -sfn networkx-<major>.<minor> stable
      git add networkx-<major>.<minor> stable
      git commit -m "Add <major>.<minor> docs"
      # maybe squash all the Deploy GitHub Pages commits
      # git rebase -i HEAD~XX where XX is the number of commits back
      # check you didn't break anything
      # diff -r latest networkx-<major>.<minor>
      # you will then need to force the push so be careful!
      git push

- Update ``__version__`` in ``networkx/__init__.py``.

- Create ``doc/_templates/layout.html`` with::

    {% extends "!layout.html" %}

    {% block document %}
      {% include "dev_banner.html" %}
      {{ super() }}
    {% endblock %}

 - Commit and push changes::

    git add networkx/__init__.py doc/_templates/layout.html
    git commit -m "Bump release version"
    git push upstream main

- Update the web frontpage:
  The webpage is kept in a separate repo: networkx/website

  - Sync your branch with the remote repo: ``git pull``.
    If you try to ``make github`` when your branch is out of sync, it
    creates headaches.
  - Update ``build/index.html``.
  - Edit ``build/_static/docversions.js`` and commit
  - Push your changes to the repo.
  - Deploy using ``make github``.

- Post release notes on mailing list.

  - networkx-discuss@googlegroups.com
Deprecations
============

.. _deprecation_policy:

Policy
------

If the behavior of the library has to be changed, a deprecation cycle must be
followed to warn users.

A deprecation cycle is *not* necessary when:

* adding a new function, or
* adding a new keyword argument to the *end* of a function signature, or
* fixing buggy behavior

A deprecation cycle is necessary for *any breaking API change*, meaning a
change where the function, invoked with the same arguments, would return a
different result after the change. This includes:

* changing the order of arguments or keyword arguments, or
* adding arguments or keyword arguments to a function, or
* changing the name of a function, class, method, etc., or
* moving a function, class, etc. to a different module, or
* changing the default value of a function's arguments.

Usually, our policy is to put in place a deprecation cycle over two minor
releases (e.g., if a deprecation warning appears in 2.3, then the functionality
should be removed in 2.5).  For major releases we usually require that all
deprecations have at least a 1-release deprecation cycle (e.g., if 3.0 occurs
after 2.5, then all removed functionality in 3.0 should be deprecated in 2.5).

Note that these 1- and 2-release deprecation cycles for major and minor
releases is not a strict rule and in some cases, the developers can agree on a
different procedure upon justification (like when we can't detect the change,
or it involves moving or deleting an entire function for example).

Todo
----

Make sure to review ``networkx/conftest.py`` after removing deprecated code.

Version 3.0
~~~~~~~~~~~

* In ``readwrite/gml.py`` remove ``literal_stringizer`` and related tests.
* In ``readwrite/gml.py`` remove ``literal_destringizer`` and related tests.
* In ``utils/misc.py`` remove ``is_string_like`` and related tests.
* In ``utils/misc.py`` remove ``make_str`` and related tests.
* In ``utils/misc.py`` remove ``is_iterator``.
* In ``utils/misc.py`` remove ``iterable``.
* In ``utils/misc.py`` remove ``is_list_of_ints``.
* In ``utils/misc.py`` remove ``consume``.
* In ``utils/misc.py`` remove ``default_opener``.
* In ``utils/misc.py`` remove ``empty_generator``.
* Remove ``utils/contextmanagers.py`` and related tests.
* In ``drawing/nx_agraph.py`` remove ``display_pygraphviz`` and related tests.
* In ``algorithms/chordal.py`` replace ``chordal_graph_cliques`` with ``_chordal_graph_cliques``.
* In ``algorithms/centrality/betweenness_centrality_subset.py`` remove ``betweenness_centrality_source``.
* In ``algorithms/centrality/betweenness.py`` remove ``edge_betweeness``.
* In ``algorithms/community_modularity_max.py`` remove old name ``_naive_greedy_modularity_communities``.
* In ``linalg/algebraicconnectivity.py`` remove ``_CholeskySolver`` and related code.
* In ``convert_matrix.py`` remove ``to_numpy_matrix`` and ``from_numpy_matrix``.
* In ``readwrite/json_graph/cytoscape.py``, change function signature for
  ``cytoscape_graph`` and ``cytoscape_data`` to replace the ``attrs`` keyword.
  argument with explicit ``name`` and ``ident`` keyword args.
* In ``readwrite/json_graph/tree.py``, remove ``attrs`` kwarg from ``tree_graph``
  and ``tree_data``.
* Undo changes related to the removal of ``pyyaml``. Remove the
  ``__getattr__`` definitions from ``networkx/__init__.py``,
  ``networkx/readwrite/__init__.py`` and ``networkx/readwrite/nx_yaml.py`` and
  remove ``networkx/readwrite/tests/test_getattr_nxyaml_removal.py``
* Remove ``readwrite/gpickle.py`` and related tests.
* Remove ``readwrite/nx_shp.py`` and related tests (add info in alternatives).
* Remove ``copy`` method in the coreview Filtered-related classes and related tests.
* In ``algorithms/link_analysis/pagerank_alg.py`` replace ``pagerank`` with ``pagerank_scipy``.
* In ``algorithms/link_analysis/pagerank_alg.py`` rename ``pagerank_numpy`` as ``_pagerank_numpy``.
* In ``convert_matrix.py`` remove ``order`` kwarg from ``to_pandas_edgelist`` and docstring
* Remove ``readwrite/json_graph/jit.py`` and related tests.
* In ``utils/misc.py`` remove ``generate_unique_node`` and related tests.
* In ``algorithms/link_analysis/hits_alg.py`` remove ``hub_matrix`` and ``authority_matrix``
* In ``algorithms/link_analysis/hits_alg.py``, remove ``hits_numpy`` and ``hist_scipy``.
* In ``networkx.classes`` remove the ``ordered`` module and the four ``Ordered``
  classes defined therein.
* In ``utils/decorators.py`` remove ``preserve_random_state``.
* In ``algorithms/community/quality.py`` remove ``coverage`` and ``performance``.
* Remove ``testing``.
* In ``linalg/graphmatrix.py`` remove ``adj_matrix``.
* In ``algorithms/similarity.py`` replace ``simrank_similarity`` with ``simrank_similarity_numpy``.
* In ``algorithms/assortativity/mixing.py`` remove ``numeric_mixing_matrix``.
* In ``algorithms/assortativity/connectivity.py`` remove ``k_nearest_neighbors``.
* In ``utils/decorators.py`` remove ``random_state``.
* In ``algorithms/operators/binary.py`` remove ``name`` kwarg from ``union`` and docstring.
* In ``networkx/generators/geometric.py`` remove ``euclidean`` and tests.
* In ``networkx/algorithms/node_classification/`` remove ``hmn.py``, ``lgc.py``,
  and ``utils.py`` after moving the functions defined therein into the newly created
  ``node_classification.py`` module, which will replace the current package.
* In ``networkx/algorithms/link_analysis/pagerank_alg.py``, remove the
  ``np.asmatrix`` wrappers on the return values of ``google_matrix`` and remove
  the associated FutureWarning.
* In ``networkx/convert_matrix.py`` remove ``from_scipy_sparse_matrix`` and
  ``to_scipy_sparse_matrix``.

Core Developers
---------------

NetworkX development is guided by the following core team:



.. raw:: html

   <div class="team-member">
     <a href="https://github.com/boothby" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/569654?u=c29b79275293c22fa3c56a06ed04e004465ef331&v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @boothby"
           />
        </div>
        Kelly Boothby
     </a>
     <div class="team-member-handle">@boothby</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/camillescott" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/2896301?u=bd57c546510c131f4f7f41e3999fb8e6e33a2298&v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @camillescott"
           />
        </div>
        Camille Scott
     </a>
     <div class="team-member-handle">@camillescott</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/dschult" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/915037?u=6a27f396c666c5c2172a1cfc7b0d4bbcd0069eed&v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @dschult"
           />
        </div>
        Dan Schult
     </a>
     <div class="team-member-handle">@dschult</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/ericmjl" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/2631566?u=c5d73d769c251a862d7d4bbf1119297d8085c34c&v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @ericmjl"
           />
        </div>
        Eric Ma
     </a>
     <div class="team-member-handle">@ericmjl</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/hagberg" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/187875?v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @hagberg"
           />
        </div>
        Aric Hagberg
     </a>
     <div class="team-member-handle">@hagberg</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/harshal-dupare" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/52428908?u=cb974ff050563c3610f377b7dbbf4982df6a1b90&v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @harshal-dupare"
           />
        </div>
        Harshal Dupare
     </a>
     <div class="team-member-handle">@harshal-dupare</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/jarrodmillman" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/123428?v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @jarrodmillman"
           />
        </div>
        Jarrod Millman
     </a>
     <div class="team-member-handle">@jarrodmillman</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/mjschwenne" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/19698215?u=96f60d2e90261aa7487fffcc2ebad1158028ecd5&v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @mjschwenne"
           />
        </div>
        Matt Schwennesen
     </a>
     <div class="team-member-handle">@mjschwenne</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/MridulS" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/5363860?u=ce5c6e9388d2fd153ebf8b0bb521c928b0813608&v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @MridulS"
           />
        </div>
        Mridul Seth
     </a>
     <div class="team-member-handle">@MridulS</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/rossbar" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/1268991?u=974707b96081a9705f3a239c0773320f353ee02f&v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @rossbar"
           />
        </div>
        Ross Barnowski
     </a>
     <div class="team-member-handle">@rossbar</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/stefanv" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/45071?u=c779b5e06448fbc638bc987cdfe305c7f9a7175e&v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @stefanv"
           />
        </div>
        Stefan van der Walt
     </a>
     <div class="team-member-handle">@stefanv</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/vdshk" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/43042296?u=01411ddb7d394274117007e8d29019e091a8e00a&v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @vdshk"
           />
        </div>
        Vadim Abzalov
     </a>
     <div class="team-member-handle">@vdshk</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/z3y50n" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/33282622?u=66483bb152faad7fcdb80cb9d1f8b6d391e448bc&v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @z3y50n"
           />
        </div>
        Dimitrios Papageorgiou
     </a>
     <div class="team-member-handle">@z3y50n</div>
   </div>



Emeritus Developers
-------------------

We thank these previously-active core developers for their contributions to NetworkX.



.. raw:: html

   <div class="team-member">
     <a href="https://github.com/bjedwards" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/726274?u=e493f38cb65425f6de7a9568ee3802a183deaa8e&v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @bjedwards"
           />
        </div>
        Benjamin Edwards
     </a>
     <div class="team-member-handle">@bjedwards</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/chebee7i" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/326005?u=a5a33cadf55b2fbdd8b033517f97f763563aa72a&v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @chebee7i"
           />
        </div>
        @chebee7i
     </a>
     <div class="team-member-handle">@chebee7i</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/jfinkels" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/121755?v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @jfinkels"
           />
        </div>
        @jfinkels
     </a>
     <div class="team-member-handle">@jfinkels</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/jtorrents" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/1184374?v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @jtorrents"
           />
        </div>
        Jordi Torrents
     </a>
     <div class="team-member-handle">@jtorrents</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/loicseguin" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/812562?v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @loicseguin"
           />
        </div>
        Loïc Séguin-Charbonneau
     </a>
     <div class="team-member-handle">@loicseguin</div>
   </div>


.. raw:: html

   <div class="team-member">
     <a href="https://github.com/ysitu" class="team-member-name">
        <div class="team-member-photo">
           <img
             src="https://avatars.githubusercontent.com/u/7018196?v=4&s=40"
             loading="lazy"
             alt="Avatar picture of @ysitu"
           />
        </div>
        @ysitu
     </a>
     <div class="team-member-handle">@ysitu</div>
   </div>

.. _NXEP3:

==================================
NXEP 3 — Graph Builders
==================================

:Author: Dan Schult
:Author: Kelly Boothby
:Status: Draft
:Type: Standards Track
:Created: 2020-11-27


Abstract
--------

Graph generators in NetworkX create a Graph starting from an object
specified in the `create_using` argument. Many of these generators
do no more than create edges to add to the graph. Sometimes all we
want the graph for is to generate those edges. It might be better
to allow generators to report either the `edgelist` or any of the
standard graph classes or a custom class. This NXEP proposes a
framework for graph builders which allows a user friendly interface
to these features and decorators to make it easy for developers to
provide these features whether the graph builder algorithm requires
a graph, or just edges.

Motivation and Scope
--------------------

Consider, for example, the function `nx.path_graph(nodes, create_using)`.
It creates the edges for the indicated path and adds them to an empty
graph data structure created using the type `create_using`.
`path_graph` does not use the graph structure to create the edges
being generated and could arguably just yield
the edges without involving the data structure at all.
The parameter `create_using` is used to indicate the type of graph data
structure to use. This could be indicated by passing an edge generator
to the type constructor e.g. `nx.MultiDiGraph(nx.path_edges(n))` instead
of `nx.path_graph(n, create_using=nx.MultiDiGraph)`. The former style
allows the edges to be used without creating any graph data structure if
that is desired. The latter is preferred stylistically because the main
idea of the code phrase is to create a path graph and that comes first.
Details such as what graph type to use are specified later in the phrase.

Separating edge generation from graph data structure creation
arguably makes a cleaner interface where independent tools can be put
together in creative ways. To the extent that users need to generate
edges rather than graphs, having an edge generator that doesn't store
the graph is an advantage. It's not exactly clear how much demand there
is for this feature. But e.g. `nx.utils.pairwise` would no longer be needed
if we had `nx.path_edges(node_iterable)`.

The `create_using` parameter is a mechanism to tell the function what
class of graph data structure to start with. Separating edge generation
from graph construction would mean the edge generator function would
no longer need a type for the graph data structure since there isn't one.
This NXEP proposes one way to provide an interface that separates edge
generation from graph data structure creation when desired, while leaving
a user-friendly mechanism for selection the graph type when desired.

The changes proposed involve any nx function that creates a graph or an
edgelist. The proposal is to make these functions return a graph of
any type or an edgelist based on the choice indicated by the user.
Developers choose whether it is more effective to yield edges or to
return a graph. Decorators are used to construct the surrounding code
to enable the other output styles.

The proposed solution is to provide the user with graph builders that
return either a graph or an edgelist while minimizing the code needed
for developers to support both. The underlying code could choose to
either 1) yield edges, or 2) construct a graph from an input graph
parameter. Two decorators would then add the extra code needed to
construct a single object so users would use the same interface no
matter which style of underlying code was used. The user facing
interface would allow the user to specify a graph data structure
by type, or request an edgelist. One syntax proposal is::

    G = nx.path_graph(9)
    DG = nx.path_graph.DiGraph(9)
    MG = nx.path_graph.MultiGraph(9)
    MDG = nx.path_graph.MultiDiGraph(9)
    CG = nx.path_graph.CustomGraph(9, create_using)
    elist = nx.path_graph.edgelist(9)


Edgelists that only contain pairs of nodes indicating an edge are restrictive.
Some graphs have isolated nodes which would not appear in any node-pair.
Some graphs have node or edge attributes associated with the node or edge.
Multigraphs have edge keys associated with each edge, often as a 3-tuple
(u, v, ekey). This proposal suggests that we adopt the following edgelist
protocol for describing a graph (perhaps there is a better name than edgelist).

An edgelist is a sequence of containers. The length of the container along
with the hashable nature of it's last element determined the type of
information included in the container. All currently used graph information
can be stored in such a sequence. The logic is as follows where C denotes
the container:

+------------------------------+--------+-----------------+
|                              | len(C) | hashable(C[-1]) |
+==============================+========+=================+
|Graph attributes:             |   1    |    False        |
+------------------------------+--------+-----------------+
|Node without attributes:      |   1    |    True         |
+------------------------------+--------+-----------------+
|Node with attributes:         |   2    |    False        |
+------------------------------+--------+-----------------+
|Edge without attributes:      |   2    |    True         |
+------------------------------+--------+-----------------+
|Edge with attributes:         |   3    |    False        |
+------------------------------+--------+-----------------+
|Multiedge without attributes: |   3    |    True         |
+------------------------------+--------+-----------------+
|Multiedge with attributes:    |   4    |    False        |
+------------------------------+--------+-----------------+

Here is some code to process such an edgelist and construct the graph
starting from an empty graph G:

.. code-block:: python

    for C in edgelist:
     if len(C) == 1:
       if not hashable(C[-1]):
         G.graph.update(C[-1])  # C[-1] is a dict of graph attributes
       else:
         G.add_node(C[-1])  # C[-1] is a node
     elif len(C) == 2:
       if not hashable(C[-1]):
         G.add_node(C[0], **C[-1])  # C[-1] is a dict of node attributes
       else:
         G.add_edge(*C)  # C is a node-pair indicating an edge
     elif len(C) == 3:
       if not hashable(C[-1]):
         G.add_edge(*C[:2], **C[-1])  # C -> (u, v, attrdict)
       else:
         G.add_edge(*C)  # C -> (u, v, edge_key)
     elif len(C) == 4:
         assert not hashable(C[-1])
         G.add_edge(*C)  # C -> (u, v, edge_key, attr_dict)
     else:
         raise NetworkXInvalidEdgelist(
             "no container in an edgelist should be larger than 4 objects."
         )

Usage and Impact
----------------

Users will build graphs using similar syntax as before with added flexibility.

Create a wheel graph with 9 spokes (10 nodes):

    >>> G = nx.wheel_graph(9)  # same as current code

Construct a path graph using a MultiDiGraph data structure:

    >>> MDG = nx.path_graph.MultiDiGraph([3, 4, 2, 5, 7, 6])
    >>> # current code:
    >>> MDG = nx.path_graph([3, 4, 2, 5, 7, 6], create_using=MultiDiGraph)

Construct a star graph using a CustomGraph subclass of a NetworkX graph class.

    >>> G = nx.star_graph.CustomGraph(9, MyCustomGraph)
    >>> # current code:
    >>> G = nx.star_graph(9, create_using=MyCustomGraph)

Add a complete graph to an existing graph G:

    >>> G.update(nx.complete_graph.edgelist(range(len(G) - 10, 20))

Iterate over the edges of a randomly generated graph without storing it.

    >>> for u, v in nx.configuration_model_graph.edgelist(deg_sequence):
    >>>     process(u, v)

Developers will use a decorator to indicate whether their graph builder
has underlying code that yields from an edgelist, or returns a graph.

.. code-block:: python

    @graph_builder
    @py_random_state(4)
    def extended_barabasi_albert_graph(n, m, p, q, seed=None)
        # some fancy code that requires we construct G to use graph properties
        # while we decide what edges to add next.
        return G

The `@graph_builder` decorator adds code to enable
e.g. `nx.extended_barabasi_albert_graph.edgelist`.

For most graph builders we simply yield from an edgelist.

.. code-block:: python

    @node_and_edge_builder
    def ladder_graph(n):
        yield from pairwise(range(n))
        yield from pairwise(range(n, 2 * n))
        yield from ((v, v + n) for v in range(n))

The `@node_and_edge_builder` decorator adds code to enable
e.g. `nx.ladder_graph.MultiGraph(6)`. Note that `nx.ladder_graph(6)`
would still return an nx.Graph as it currently does. To make use of the
edgelist functionality yielding edge without graph constructing, the syntax
would be `nx.ladder_graph.edgelist(6)`.


Backward compatibility
----------------------

To reduce backward incompatibility, the base calling structure `nx.path_graph(9)`
works as it currently does. The `create_using` parameter is removed and
replaced by an attribute of the calling function.
So `nx.path_graph(9, nx.DiGraph)` becomes `nx.path_graph.DiGraph(9)`.
The `create_using` parameter could also be retained providing more backward
compatibility at the potential cost of providing at least 2 ways to create
the same graph: `nx.path_graph(9, create_using=nx.DiGraph)`
and `nx.path_graph.DiGraph(9)`. See the Alternatives section.


Due to the renaming of graph generators as graph builders (to avoid confusion
with Python's generator functions) anyone using full-path calling syntax
e.g., `nx.generators.path_graph(9)` will need to change to `nx.path_graph(9)`
or `nx.builders.path_graph(9)` though the latter is discouraged.
This change of name is independent of the main thrust of this proposal.
But it seems a reasonable time to make such a change.

To reduce developer impact, upon inception, we could use all current graph
generators as graph builders by attaching the `@graph_builder` decorator.
Presumably for efficiency many of them should be rewritten to yield
edgelists rather than returning graphs. But this could be done gradually
and when done switch the decorator to `@node_and_edge_builder`. Both should
return equivalent graph builder objects.


Detailed description
--------------------

This can be accomplished through a couple decorators, which could be
adopted gradually -- a big patch initially decorating all existing generators
with `@graph_builder` would immediately support the notation
`nx.complete_graph.edgelist(...)` without impacting existing code.
Later generators could use `@node_and_edge_builder`.

.. code-block:: python

    def node_and_edge_builder(f):
        @wraps(f)
        def graph(*args, **kwargs):
            return nx.Graph(f(*args, **kwargs))
        def digraph(*args, **kwargs):
            return nx.DiGraph(f(*args, **kwargs))
        def multigraph(*args, **kwargs):
            return nx.MultiGraph(f(*args, **kwargs))
        def multidigraph(*args, **kwargs):
            return nx.MultiDiGraph(f(*args, **kwargs))
        def custom_graph(*args, create_using=None, **kwargs):
            g = create_using()
            g.update(f(*args, **kwargs))
            return g
        graph.Graph = graph
        graph.DiGraph = digraph
        graph.MultiGraph = multigraph
        graph.MultiDiGraph = multidigraph
        graph.CustomGraph = custom_graph
        graph.edgelist = f
        return graph

    def graph_builder(f):
        @wraps(f)
        def edgelist(*args, **kwargs):
            g = f(*args, **kwargs)
            return itertools.ichain(
                map(tuple, G.nodes.data()), map(tuple, G.edges.data())
            )
        f.edgelist = edgelist
        f.CustomGraph = f
        def graph(*args, **kwargs):
            return f(*args, create_using=nx.Graph, **kwargs))
        def digraph(*args, **kwargs):
            return f(*args, create_using=nx.DiGraph, **kwargs))
        def multigraph(*args, **kwargs):
            return f(*args, create_using=nx.MultiGraph, **kwargs))
        def multidigraph(*args, **kwargs):
            return f(*args, create_using=nx.MultiDiGraph, **kwargs))
        f.Graph = graph
        f.DiGraph = digraph
        f.MultiGraph = multigraph
        f.MultiDiGraph = multidigraph
        return f

Note: the graph_builder underlying code should accept a create_using
parameter for this implementation to work. We need to think if this is
universally applicable and how to handle builders that shouldn't work
with all four of the major NetworkX graph classes.

Graph.update will need to handle an edgelist input. It currently handles
node-pairs and node-pair with edge key triples for multigraphs. Code like
that shown above in the description of Edgelist should be used.

Example developer usage:

.. code-block:: python

    @node_and_edge_builder
    def path_graph(n):
        """an overly simplified path graph implementation"""
        return pairwise(range(n))

    @graph_builder
    def complete_graph(n, create_using=None):
        """an overly simplified complete graph implementation"""
        if create_using is None:
            create_using = nx.Graph
        g = empty_graph(0, create_using)
        g.update(itertools.combinations(range(n), 2))
        return g


Related Work
------------

This proposal is based on ideas and discussions from #3036 and #1393.


Implementation
--------------

The first major step is to implement the two builder decorators.
Next we need to change the Graph update methods, convert functions, etc.
to process edgelists that contain isolated nodes and data attributes.
Third we should identify any functions that build graphs or edgelists
and decorate them to make them Graph Builders.

Special care should be made to ensure only desired graph types are
accepted and appropriate errors raised when not.

We should rename the generators directory as builders and adjust
documentation where needed appropriately (including old documentation
getting the correct canonical url).

Later steps include going through the existing generator code and switching
that code to yield edgelists instead of returning graphs (where appropriate).


Alternatives
------------

We can just leave the generators as they are and deal with the cost of
creating a graph when one only needs the edgelist. It's not a huge cost
most of the time.

We can split the edge generation from graph creation using
`nx.DiGraph(nx.path_edgelist(9))` and disallowing `create_using`.

We can implement the proposal retaining the `create_using` parameter
for backward compatibility.


Discussion
----------

Most of the ideas here are from
- [`#3036 <https://github.com/networkx/networkx/pull/3036>`]
which built on discussion from
- [`#1393 <https://github.com/networkx/networkx/pull/1393>`]
.. _governance:

=======================================
NXEP 1 — Governance and Decision Making
=======================================

:Author: Jarrod Millman <millman@berkeley.edu>
:Author: Dan Schult <dschult@colgate.edu>
:Status: Draft
:Type: Process
:Created: 2020-06-25

Abstract
========

NetworkX is a consensus-based community project. Anyone with an interest in the
project can join the community, contribute to the project design, and
participate in the decision making process. This document describes how that
participation takes place, how to find consensus, and how deadlocks are
resolved.

Roles And Responsibilities
==========================

The Community
-------------
The NetworkX community consists of anyone using or working with the project
in any way.

Contributors
------------
Any community member can become a contributor by interacting directly with the
project in concrete ways, such as:

- proposing a change to the code or documentation via a GitHub pull request;
- reporting issues on our
  `GitHub issues page <https://github.com/networkx/networkx/issues>`_;
- discussing the design of the library, website, or tutorials on the
  `mailing list <http://groups.google.com/group/networkx-discuss/>`_,
  or in existing issues and pull requests; or
- reviewing
  `open pull requests <https://github.com/networkx/networkx/pulls>`_,

among other possibilities. By contributing to the project, community members
can directly help to shape its future.

Contributors should read the :ref:`contributor_guide` and our :ref:`code_of_conduct`.

Core Developers
---------------
Core developers are community members that have demonstrated continued
commitment to the project through ongoing contributions. They
have shown they can be trusted to maintain NetworkX with care. Becoming a
core developer allows contributors to merge approved pull requests, cast votes
for and against merging a pull request, and be involved in deciding major
changes to the API, and thereby more easily carry on with their project related
activities. Core developers appear as team members on the `NetworkX Core Team page
<https://github.com/orgs/networkx/teams/core-developers/members>`_ and can
be messaged ``@networkx/core-developers``. Core
developers are expected to review code contributions while adhering to the
:ref:`core_dev`.

New core developers can be nominated by any existing core developer.
Discussion about new core developer nominations is one of the few activities
that takes place on the project's private management list. The decision to
invite a new core developer must be made by “lazy consensus”, meaning unanimous
agreement by all responding existing core developers. Invitation must take
place at least one week after initial nomination, to allow existing members
time to voice any objections.

.. _steering_council:

Steering Council
----------------
The Steering Council (SC) members are core developers who have additional
responsibilities to ensure the smooth running of the project. SC members are
expected to participate in strategic planning, approve changes to the
governance model, and make decisions about funding granted to the project
itself. (Funding to community members is theirs to pursue and manage.) The
purpose of the SC is to ensure smooth progress from the big-picture
perspective. Changes that impact the full project require analysis informed by
long experience with both the project and the larger ecosystem. When the core
developer community (including the SC members) fails to reach such a consensus
in a reasonable timeframe, the SC is the entity that resolves the issue.

Steering Council members appear as team members on the `NetworkX Steering
Council Team page
<https://github.com/orgs/networkx/teams/steering-council/members>`_ and
can be messaged ``@networkx/steering-council``. Core

Decision Making Process
=======================

Decisions about the future of the project are made through discussion with all
members of the community. All non-sensitive project management discussion takes
place on the project
`mailing list <http://groups.google.com/group/networkx-discuss/>`_
and the `issue tracker <https://github.com/networkx/networkx/issues>`_.
Occasionally, sensitive discussion may occur on a private list.

Decisions should be made in accordance with our :ref:`mission_and_values`.

NetworkX uses a *consensus seeking* process for making decisions. The group
tries to find a resolution that has no open objections among core developers.
Core developers are expected to distinguish between fundamental objections to a
proposal and minor perceived flaws that they can live with, and not hold up the
decision making process for the latter.  If no option can be found without
an objection, the decision is escalated to the SC, which will itself use
consensus seeking to come to a resolution. In the unlikely event that there is
still a deadlock, the proposal will move forward if it has the support of a
simple majority of the SC. Any proposal must be described by a NetworkX :ref:`nxep`.

Decisions (in addition to adding core developers and SC membership as above)
are made according to the following rules:

- **Minor documentation changes**, such as typo fixes, or addition / correction of a
  sentence (but no change of the NetworkX landing page or the “about”
  page), require approval by a core developer *and* no disagreement or requested
  changes by a core developer on the issue or pull request page (lazy
  consensus). Core developers are expected to give “reasonable time” to others
  to give their opinion on the pull request if they’re not confident others
  would agree.

- **Code changes and major documentation changes** require agreement by *two*
  core developers *and* no disagreement or requested changes by a core developer
  on the issue or pull-request page (lazy consensus).

- **Changes to the API principles** require a :ref:`nxep` and follow the
  decision-making process outlined above.

- **Changes to this governance model or our mission and values**
  require a :ref:`nxep` and follow the decision-making process outlined above,
  *unless* there is unanimous agreement from core developers on the change.

If an objection is raised on a lazy consensus, the proposer can appeal to the
community and core developers and the change can be approved or rejected by
escalating to the SC, and if necessary, a NXEP (see below).

.. _nxep:

Enhancement Proposals (NXEPs)
=============================

Any proposals for enhancements of NetworkX should be written as a formal NXEP
following the template :doc:`nxep-template`. The NXEP must be made public and
discussed before any vote is taken. The discussion must be summarized by a
key advocate of the proposal in the appropriate section of the NXEP.
Once this summary is made public and after sufficient time to allow the
core team to understand it, they vote.
The workflow of a NXEP is detailed in :ref:`nxep0`.

A list of all existing NXEPs is available :ref:`here <nxep_list>`.

Acknowledgments
===============

This document is based on the `scikit-image governance document
<https://scikit-image.org/docs/stable/skips/1-governance.html>`_.
==================================
NXEP X — Template and Instructions
==================================

:Author: <list of authors' real names and optionally, email addresses>
:Status: <Draft | Active | Accepted | Deferred | Rejected | Withdrawn | Final | Superseded>
:Type: <Standards Track | Process>
:Created: <date created on, in yyyy-mm-dd format>
:Resolution: <url> (required for Accepted | Rejected | Withdrawn)


Abstract
--------

The abstract should be a short description of what the NXEP will achieve.

Note that the — in the title is an elongated dash, not -.

Motivation and Scope
--------------------

This section describes the need for the proposed change. It should describe
the existing problem, who it affects, what it is trying to solve, and why.
This section should explicitly address the scope of and key requirements for
the proposed change.

Usage and Impact
----------------

This section describes how users of NetworkX will use features described in this
NXEP. It should be comprised mainly of code examples that wouldn't be possible
without acceptance and implementation of this NXEP, as well as the impact the
proposed changes would have on the ecosystem. This section should be written
from the perspective of the users of NetworkX, and the benefits it will provide
them; and as such, it should include implementation details only if
necessary to explain the functionality.

Backward compatibility
----------------------

This section describes the ways in which the NXEP breaks backward compatibility.

The mailing list post will contain the NXEP up to and including this section.
Its purpose is to provide a high-level summary to users who are not interested
in detailed technical discussion, but may have opinions around, e.g., usage and
impact.

Detailed description
--------------------

This section should provide a detailed description of the proposed change.
It should include examples of how the new functionality would be used,
intended use-cases and pseudo-code illustrating its use.


Related Work
------------

This section should list relevant and/or similar technologies, possibly in other
libraries. It does not need to be comprehensive, just list the major examples of
prior and relevant art.


Implementation
--------------

This section lists the major steps required to implement the NXEP.  Where
possible, it should be noted where one step is dependent on another, and which
steps may be optionally omitted.  Where it makes sense, each step should
include a link to related pull requests as the implementation progresses.

Any pull requests or development branches containing work on this NXEP should
be linked to from here.  (A NXEP does not need to be implemented in a single
pull request if it makes sense to implement it in discrete phases).


Alternatives
------------

If there were any alternative solutions to solving the same problem, they should
be discussed here, along with a justification for the chosen approach.


Discussion
----------

This section may just be a bullet list including links to any discussions
regarding the NXEP:

- This includes links to mailing list threads or relevant GitHub issues.
.. _NXEP0:

============================
NXEP 0 — Purpose and Process
============================

:Author: Jarrod Millman <millman@berkeley.edu>
:Status: Draft
:Type: Process
:Created: 2020-06-25


What is a NXEP?
---------------


NXEP stands for NetworkX Enhancement Proposal.  NXEPs are the primary
mechanisms for proposing major new features, for collecting community input on
an issue, and for documenting the design decisions that have gone into
NetworkX.  A NXEP should provide a concise technical specification of the
feature and a rationale for the feature.  The NXEP author is responsible for
building consensus within the community and documenting dissenting opinions.

Because the NXEPs are maintained as text files in a versioned
repository, their revision history is the historical record of the
feature proposal [1]_.


Types
^^^^^

There are three kinds of NXEPs:

1. A **Standards Track** NXEP describes a new feature or implementation
   for NetworkX.

2. An **Informational** NXEP describes a NetworkX design issue, or provides
   general guidelines or information to the Python community, but does not
   propose a new feature. Informational NXEPs do not necessarily represent a
   NetworkX community consensus or recommendation, so users and implementers are
   free to ignore Informational NXEPs or follow their advice.

3. A **Process** NXEP describes a process surrounding NetworkX, or
   proposes a change to (or an event in) a process.  Process NXEPs are
   like Standards Track NXEPs but apply to areas other than the NetworkX
   language itself.  They may propose an implementation, but not to
   NetworkX's codebase; they require community consensus.  Examples include
   procedures, guidelines, changes to the decision-making process, and
   changes to the tools or environment used in NetworkX development.
   Any meta-NXEP is also considered a Process NXEP.


NXEP Workflow
-------------

The NXEP process begins with a new idea for NetworkX.  It is highly
recommended that a single NXEP contain a single key proposal or new
idea. Small enhancements or patches often don't need
a NXEP and can be injected into the NetworkX development workflow with a
pull request to the NetworkX `repo`_. The more focused the
NXEP, the more successful it tends to be.
If in doubt, split your NXEP into several well-focused ones.

Each NXEP must have a champion---someone who writes the NXEP using the style
and format described below, shepherds the discussions in the appropriate
forums, and attempts to build community consensus around the idea.  The NXEP
champion (a.k.a. Author) should first attempt to ascertain whether the idea is
suitable for a NXEP. Posting to the networkx-discussion `mailing list`_ is the best
way to go about doing this.

The proposal should be submitted as a draft NXEP via a `GitHub pull
request`_ to the ``doc/nxeps`` directory with the name ``nxep-<n>.rst``
where ``<n>`` is an appropriately assigned four-digit number (e.g.,
``nxep-0000.rst``). The draft must use the :doc:`nxep-template` file.

Once the PR for the NXEP is in place, a post should be made to the
mailing list containing the sections up to "Backward compatibility",
with the purpose of limiting discussion there to usage and impact.
Discussion on the pull request will have a broader scope, also including
details of implementation.

At the earliest convenience, the PR should be merged (regardless of
whether it is accepted during discussion).  Additional PRs may be made
by the Author to update or expand the NXEP, or by maintainers to set
its status, discussion URL, etc.

Standards Track NXEPs consist of two parts, a design document and a
reference implementation.  It is generally recommended that at least a
prototype implementation be co-developed with the NXEP, as ideas that sound
good in principle sometimes turn out to be impractical when subjected to the
test of implementation.  Often it makes sense for the prototype implementation
to be made available as PR to the NetworkX repo (making sure to appropriately
mark the PR as a WIP).


Review and Resolution
^^^^^^^^^^^^^^^^^^^^^

NXEPs are discussed on the mailing list.  The possible paths of the
status of NXEPs are as follows:

.. image:: _static/nxep-0000.png

All NXEPs should be created with the ``Draft`` status.

Eventually, after discussion, there may be a consensus that the NXEP
should be accepted – see the next section for details. At this point
the status becomes ``Accepted``.

Once a NXEP has been ``Accepted``, the reference implementation must be
completed.  When the reference implementation is complete and incorporated
into the main source code repository, the status will be changed to ``Final``.

To allow gathering of additional design and interface feedback before
committing to long term stability for a language feature or standard library
API, a NXEP may also be marked as "Provisional". This is short for
"Provisionally Accepted", and indicates that the proposal has been accepted for
inclusion in the reference implementation, but additional user feedback is
needed before the full design can be considered "Final". Unlike regular
accepted NXEPs, provisionally accepted NXEPs may still be Rejected or Withdrawn
even after the related changes have been included in a Python release.

Wherever possible, it is considered preferable to reduce the scope of a
proposal to avoid the need to rely on the "Provisional" status (e.g. by
deferring some features to later NXEPs), as this status can lead to version
compatibility challenges in the wider NetworkX ecosystem.

A NXEP can also be assigned status ``Deferred``.  The NXEP author or a
core developer can assign the NXEP this status when no progress is being made
on the NXEP.

A NXEP can also be ``Rejected``.  Perhaps after all is said and done it
was not a good idea.  It is still important to have a record of this
fact. The ``Withdrawn`` status is similar---it means that the NXEP author
themselves has decided that the NXEP is actually a bad idea, or has
accepted that a competing proposal is a better alternative.

When a NXEP is ``Accepted``, ``Rejected``, or ``Withdrawn``, the NXEP should be
updated accordingly. In addition to updating the status field, at the very
least the ``Resolution`` header should be added with a link to the relevant
thread in the mailing list archives.

NXEPs can also be ``Superseded`` by a different NXEP, rendering the
original obsolete.  The ``Replaced-By`` and ``Replaces`` headers
should be added to the original and new NXEPs respectively.

Process NXEPs may also have a status of ``Active`` if they are never
meant to be completed, e.g. NXEP 0 (this NXEP).


How a NXEP becomes Accepted
^^^^^^^^^^^^^^^^^^^^^^^^^^^

A NXEP is ``Accepted`` by consensus of all interested contributors. We
need a concrete way to tell whether consensus has been reached. When
you think a NXEP is ready to accept, send an email to the
networkx-discussion mailing list with a subject like:

  Proposal to accept NXEP #<number>: <title>

In the body of your email, you should:

* link to the latest version of the NXEP,

* briefly describe any major points of contention and how they were
  resolved,

* include a sentence like: "If there are no substantive objections
  within 7 days from this email, then the NXEP will be accepted; see
  NXEP 0 for more details."

For an example, see: https://mail.python.org/pipermail/networkx-discussion/2018-June/078345.html

After you send the email, you should make sure to link to the email
thread from the ``Discussion`` section of the NXEP, so that people can
find it later.

Generally the NXEP author will be the one to send this email, but
anyone can do it – the important thing is to make sure that everyone
knows when a NXEP is on the verge of acceptance, and give them a final
chance to respond. If there's some special reason to extend this final
comment period beyond 7 days, then that's fine, just say so in the
email. You shouldn't do less than 7 days, because sometimes people are
travelling or similar and need some time to respond.

In general, the goal is to make sure that the community has consensus,
not provide a rigid policy for people to try to game. When in doubt,
err on the side of asking for more feedback and looking for
opportunities to compromise.

If the final comment period passes without any substantive objections,
then the NXEP can officially be marked ``Accepted``. You should send a
followup email notifying the list (celebratory emoji optional but
encouraged 🎉✨), and then update the NXEP by setting its ``:Status:``
to ``Accepted``, and its ``:Resolution:`` header to a link to your
followup email.

If there *are* substantive objections, then the NXEP remains in
``Draft`` state, discussion continues as normal, and it can be
proposed for acceptance again later once the objections are resolved.

In unusual cases, disagreements about the direction or approach may
require escalation to the NetworkX :ref:`steering_council` who
then decide whether a controversial NXEP is ``Accepted``.


Maintenance
^^^^^^^^^^^

In general, Standards track NXEPs are no longer modified after they have
reached the Final state as the code and project documentation are considered
the ultimate reference for the implemented feature.
However, finalized Standards track NXEPs may be updated as needed.

Process NXEPs may be updated over time to reflect changes
to development practices and other details. The precise process followed in
these cases will depend on the nature and purpose of the NXEP being updated.


Format and Template
-------------------

NXEPs are UTF-8 encoded text files using the reStructuredText_ format.  Please
see the :doc:`nxep-template` file and the reStructuredTextPrimer_ for more
information.  We use Sphinx_ to convert NXEPs to HTML for viewing on the web
[2]_.


Header Preamble
^^^^^^^^^^^^^^^

Each NXEP must begin with a header preamble.  The headers
must appear in the following order.  Headers marked with ``*`` are
optional.  All other headers are required. ::

    :Author: <list of authors' real names and optionally, email addresses>
    :Status: <Draft | Active | Accepted | Deferred | Rejected |
             Withdrawn | Final | Superseded>
    :Type: <Standards Track | Process>
    :Created: <date created on, in dd-mmm-yyyy format>
  * :Requires: <nxep numbers>
  * :NetworkX-Version: <version number>
  * :Replaces: <nxep number>
  * :Replaced-By: <nxep number>
  * :Resolution: <url>

The Author header lists the names, and optionally the email addresses
of all the authors of the NXEP.  The format of the Author header
value must be

    Random J. User <address@dom.ain>

if the email address is included, and just

    Random J. User

if the address is not given.  If there are multiple authors, each should be on
a separate line.


References and Footnotes
------------------------

.. [1] This historical record is available by the normal git commands
   for retrieving older revisions, and can also be browsed on
   `GitHub <https://github.com/networkx/networkx/tree/main/doc/developer/nxeps>`_.

.. [2] The URL for viewing NXEPs on the web is
   https://networkx.org/documentation/latest/developer/nxeps/index.html

.. _repo: https://github.com/networkx/networkx

.. _mailing list: https://groups.google.com/group/networkx-discuss/

.. _issue tracker: https://github.com/networkx/networkx/issues

.. _`GitHub pull request`: https://github.com/networkx/networkx/pulls

.. _reStructuredText: http://docutils.sourceforge.net/rst.html

.. _reStructuredTextPrimer: http://www.sphinx-doc.org/en/stable/rest.html

.. _Sphinx: http://www.sphinx-doc.org/en/stable/
.. _NXEP2:

==================================
NXEP 2 — API design of view slices
==================================

:Author: Mridul Seth
:Status: Draft
:Type: Standards Track
:Created: 2020-07-23


Abstract
--------

Iterating over a subset of nodes or edges in a graph is a very common 
operation in networkx analysis.
The graph classes in NetworkX (e.g. :class:`~networkx.Graph`,
:class:`~networkx.DiGraph`, :class:`~networkx.MultiGraph`, etc.) expose the
node and edge data of the graph via :meth:`~networkx.Graph.nodes` and
:meth:`~networkx.Graph.edges`, which return dict view objects, `NodeView`
(or `NodeDataView`) and `EdgeView` (or `EdgeDataView`), respectively.
The node and edge `View` classes have dict-like semantics for item access,
returning the data dict corresponding to a given node or edge.
This NXEP proposes adding support for slicing to the relevant node & edge
`View` classes.

Motivation and Scope
--------------------

While accessing Graph data with `G.nodes` and `G.edges`, the only way of slicing the data
is by casting the view to a list manually and then calling a slice on it.
A slice inherently implies an ordering of the elements. We intend to use the ordering
imposed on the nodes and edges by the iteration order (due to the adjacency data structure).

``G.nodes(data=True)`` returns a NodeDataView of all the nodes, ``G.nodes(data=True)[x]`` returns an attribute dictionary for the node x.
The current way of getting a slice out of the underlying dict view is to cast it to list and then
slice it ``list(G.nodes(data=True))[0:10]``. This bit of code is something that is written a lot of times
by users. For graphs with a lot of nodes and edges, ``G.nodes`` and ``G.edges`` will take a lot of screen space and
when the users try to slice the resulting view (the first instinct) it will error out. Users definitely need to go through
a couple of documentation links before they realise that they need to first cast this NodeDataView to a list and then create
a slice. Updating the documentation to make this more clear would be helpful.
But it also seems good to ease the complexity of this common idiom.

In this NXEP we propose to move the casting as list inside the Node(data)View methods.
Thus ``list(G.nodes(data=True))[0:10]`` either becomes ``G.nodes(data=True)[0:10]``
or it is provided by a new slicing method like ``G.nodes(data=True).slice(10)``
or a new slicing object to allow subscripting like ``G.nodes(data=True).slice[0:10:2]``.
Then users can get a small subset of nodes by creating a slice.

Motivating Use-Case
~~~~~~~~~~~~~~~~~~~

It is common to use :meth:`~networkx.Graph.nodes` and 
:meth:`~networkx.Graph.edges` when using NetworkX interactively, e.g. in a
terminal.
If a graph has very many components (i.e. edges or nodes) then the `repr` of 
`View` object may be very long::

   >>> G = nx.complete_graph(100)   # A graph with 4950 edges
   >>> G.edges                      # Output suppressed
   
In this case, the first instinct of the user is often to inspect only the first
few edges, say 10, via slicing::

   >>> G.edges[0:10]
   Traceback (most recent call last)
      ...
   TypeError: cannot unpack non-iterable slice object

The resulting `TypeError` is opaque and hard to understand in the context of 
what was originally intended.

Usage and Impact
----------------

The main impact and the decision that needs to be taken in this NXEP is with
respect to the user facing API. By implementing this NXEP via subscripting NodeViews,
we may end up adding some ambiguity for users. As for example `G.nodes[x]`
will return an attribute dict but `G.nodes[0:5]` will return a list of first five nodes.
This will be more ambigious with EdgeView as ``G.edges[0, 1]`` will return an
attribute dictionary of the edge between 0 and 1 and ``G.edges[0:1]`` will return the first edge.
We need to find a way to counter this potential confusion.
The alternative proposal of a new slicing method is one possible solution.

For a historical context, in pre 2.0 NetworkX, G.nodes() and G.edges() returned lists.
So, slicing was native behavior like ``G.nodes()[:10]``. One caveat is that the order
of that list could change from one call to the next if the adjacency structure changed
between calls.

In more detail, in pre 2.0 NetworkX, there were 3 ways to access node information:

- ``G.node`` was a dict keyed by node to that node's attribute dict as a value.
- ``G.nodes()`` returned a list.
- ``G.nodes_iter()`` returned an iterator over the nodes.

In line with Python 3's move toward returning dict views and iterators rather than lists,
NetworkX 2.0 introduced a single interface for node information. ``G.nodes`` is a
dict-like object keyed by node to that node's attribute dict.
It also provides set-like operations on the nodes. And it offers a method ``G.nodes.data``
which provides an interface similar to ``dict.items`` but pulling out specific attributes
from the inner attribute dict rather than the entire dict. Functional synonyms
``G.nodes(data="cost", default=1)`` and ``G.nodes.data("cost", 1)`` allow an interface
that looks like a dict keyed by node to a specific node attribute.

Slicing was not provided in NetworkX 2.0 primarily because there was
no inherent order to the nodes or edges as stored in the
dict-of-dict-of-dict data structure. However, in Python 3.6, dicts
became ordered based on insertion order. So, nodes are ordered based
on when they were added to the graph and edges are ordered based on the
adjacency dict-of-dict structure. So, there is now a concept of the "first edge".

With this NXEP we would like to bring the intuitiveness
of slicing behavior back to ``G.edges`` and ``G.nodes`` using the node
add order and edge order based on adjacency storage.

On the computational front, if we create lists to allow slices, we use memory to store the lists.
This is something user would have anyway done with something like ``list(G.nodes(data=True))[0:10]``.
But we can do better with our slicing mechanisms.
We should be able to avoid constucting the entire list simply to get the slices by internally
using code like: ``indx=[n for i, n in enumerate(G.nodes(data=True)) if i in range(x.start, x.stop, s.step)]``
where x is the desired slice object.

Backward compatibility
----------------------

N/A

Detailed description
--------------------

The new implementation will let users slice Node(Data)View and Edge(Data)View.

The following code will be valid::

  >>> G.nodes(data=True)[0:10]
  >>> G.nodes[3:10]
  >>> G.edges[1:10]
  >>> G.edges(data=True)[4:6]

Prelimanary impelementation work is available at https://github.com/networkx/networkx/pull/4086

Alternatively, to get rid of the ambiguity in slicing API with respect to
the dict views we can implement a new
``slice`` method which leads to a less ambigious API.::

  >>> G.nodes(data=True).slice[:10]
  >>> G.nodes.slice[10:30]
  >>> G.edges.slice[10:40]
  >>> G.edges(data=True).slice[5:]


Related Work
------------

N/A


Implementation
--------------

A reference implementation is proposed in 
`#4086 <https://github.com/networkx/networkx/pull/4086/files>`_.

The core of this NXEP is to implement ``slicing`` to Node(Data)View
and Edge(Data)View to allow users to access a subset of nodes and edges without casting them
first to a list. We will do this by adding a check of ``slice`` in the getitem dunder method of
Node(Data)View and Edge(Data)View and returning a list of the sliced values.
For example, the `__getitem__` method for `NodeView` might look something like:

.. code-block:: python

    def __getitem__(self, n):
        if isinstance(n, slice):
            return list(self._nodes).__getitem__(n)
        return self._nodes[n]


We can instead move the check for ``slice`` to an independent ``slice`` method for nodes and edges to
implement this NXEP.

Alternatives
------------

The following list summarizes some alternatives to modifying the `__getitem__`
of the various `View` classes.
The listed alternatives are not mutually exclusive.

- **Improved Documentation** - Add more explicit documentation about the 
  necessity of casting Node(Data)View and Edge(Data)View objects to lists in
  order to be able to use slicing.
- **Improved Exceptions** - Currently, users see the following exception when
  attempting to slice a `View`::

     >>> G.nodes[0:10]
     Traceback (most recent call last)
        ...
     TypeError: unhashable type: 'slice'

  The exception message is not very useful in the context of accessing a subset
  of nodes or edges of a graph. 
  A more specific exception message could be something along the lines of::

     >>> G.nodes[0:10]
     Traceback (most recent call last)
        ...
     NetworkXError: NodeView does not support slicing. Try list(G.nodes)[0:10].

- Instead of changing the behavior of ``__getitem__`` we can impelment a new
  method, something like ``G.nodes.head(x)`` (insipired by pandas) which
  returns the first x nodes.
  This approach could be expanded to using a ``slice`` object directly but
  interfacing it with an independent ``slice`` method of G.nodes and G.edges
  instead of implementing it in getitem dunder method.

  - The nice colon syntax for slices is only available with subscript notation.
    To allow G.nodes.slice to use the nice colon syntax, we could make it a
    property that creates a subscriptable object. Syntax would be ``G.nodes.slice[4:9:2]``.


Discussion
----------

- https://github.com/networkx/networkx/pull/4086

The motivating example for the NXEP is the use-case where users want to
introspect a subset (usually the first few) of the nodes and/or edges.
If we look at the changes proposed by this NXEP and the listed alternatives,
there are several ways that this use-case might be improved.

1. Add a descriptive error message when users try to access ``View`` objects
   with a slice object.
2. Add specialized methods to the slice object (e.g. ``head()`` and ``tail()``
   or ``slice()`` that provide functionality useful for introspection.
3. The approach this NXEP proposes - modify ``View.__getitem__`` to add
   Sequence semantics.

Option 1 (better error messages) changes neither API nor behavior and would
help guide users to the correct solution for the introspection use-case.
The downside is that it does not offer the same level of convenience that
support for slicing does.

Option 2 (``head``, ``tail``, and/or ``slice`` methods) would add new methods
to view a subset of the nodes/edges.
For example::

   >>> G = nx.path_graph(10)
   >>> G.nodes()
   NodeView((0, 1, 2, 3, 4, 5, 6, 7, 8, 9))
   >>> G.nodes().head(3)   # Display the first three nodes
   NodeView((0, 1, 2))

One drawback of the approach is that is introduces new API, which has to be
both discoverable and intuitive in order to make node/edge viewing more
convenient.
For example, is ``G.nodes().head(3)`` or ``G.nodes().slice(0, 10, 2)``
more convenient than ``list(G.nodes())[:3]`` or ``list(G.nodes())[0:10:2]``,
respectively?
Another complication involves choosing the names for the new methods.
``head`` and ``tail`` are intuitive for users coming from `*nix` backgrounds
and have been adopted by other popular libraries like `pandas`.
However, ``head`` and ``tail`` also have meaning in the context of network
science pertaining to e.g. graph edges.
For example, a user might reasonably assume that ``G.edges().tail()`` would
give the set of source nodes in a directed graph, instead of the last `n`
edges.

Option 3 (add sequence semantics to `View` objects) is arguably the most
convenient as it doesn't involve raising any error messages.
However, overriding the behavior of `*View.__getitem__` to mix Mapping and
Sequence semantics is a relatively pervasive change that may have
unforeseen consequences for some use-cases.
Furthermore there is precedent in Python itself for returning un-sliceable view
objects from some mappings, a notable example being the `dict_keys` and
`dict_values` objects returned when accessing components in dictionaries::

   >>> d = {k:v for k, v in zip(range(10), range(10))}
   >>> d.values()[3:6]
   Traceback (most recent call last)
      ...
   TypeError: 'dict_values' object is not subscriptable
   >>> list(d.values())[3:6]
   [3, 4, 5]
    
Since Python dictionaries are now ordered by default (as of 3.6 in CPython),
this behavior may change in the future.

Given the considerations associated with the listed options, the following
course of action is proposed:

- **Adopt option 1** - more informative error messages for the motivating
  use-case (e.g. ``G.edges()[0:10]``) alleviates the need for users to go
  digging through the documentation to find/remember how to get the
  desired behavior.
  Since no new API is introduced nor are there any backwards compatibility
  concerns, this change doesn't require any further design discussion.
  It is possible that this change is enough to resolve the motivating
  use-case satisfactorily - monitor user feedback.
- Option 2 doesn't require any further discussion in a design doc (i.e. NXEP).
  New methods along the lines discussed above can be proposed via PR.
- Defer implementing option 3 for now, but reconsider if:

   - The improved error message is not in itself a sufficient solution
   - Other use-cases are identified for which adding slicing to the `*View`
     objects would be a nice improvement (e.g. improved performance).
.. _nxep_list:

NXEPs
*****

NetworkX Enhancement Proposals (NXEPs) document major changes or proposals.

.. toctree::
   :maxdepth: 1

   nxep-0000
   nxep-0001
   nxep-0002
   nxep-0003

.. toctree::
   :hidden:

   nxep-template
NetworkX 1.0
============

Release date:  8 Jan 2010

Version 1.0 requires Python 2.4 or greater.


New features
------------
This release has significant changes to parts of the graph API
to allow graph, node, and edge attributes.
See http://networkx.lanl.gov/reference/api_changes.html

 - Update Graph, DiGraph, and MultiGraph classes to allow attributes.
 - Default edge data is now an empty dictionary (was the integer 1)
 - Difference and intersection operators
 - Average shortest path
 - A* (A-Star) algorithm
 - PageRank, HITS, and eigenvector centrality
 - Read Pajek files
 - Line graphs
 - Minimum spanning tree (Kruskal's algorithm)
 - Dense and sparse Fruchterman-Reingold layout
 - Random clustered graph generator
 - Directed scale-free graph generator
 - Faster random regular graph generator
 - Improved edge color and label drawing with Matplotlib
 - and much more, see  https://networkx.lanl.gov/trac/query?status=closed&group=milestone&milestone=networkx-1.0

Examples
--------
 - Update to work with networkx-1.0 API
 - Graph subclass example


Version numbering
-----------------

In the future we will use a more standard release numbering system
with major.minor[build] labels where major and minor are numbers and
[build] is a label such as "dev1379" to indicate a development version
or "rc1" to indicate a release candidate.

We plan on sticking closer to a time-based release schedule with smaller
incremental changes released on a roughly quarterly basis.  The graph
classes API will remain fixed, unless we determine there are serious
bugs or other defects in the existing classes, until networkx-2.0 is
released at some time in the future.

Changes in base classes
-----------------------

The most significant changes in are in the graph classes.  All of the
graph classes now allow optional graph, node, and edge attributes.  Those
attributes are stored internally in the graph classes as dictionaries
and can be accessed simply like Python dictionaries in most cases.

Graph attributes
^^^^^^^^^^^^^^^^
Each graph keeps a dictionary of key=value attributes
in the member G.graph.  These attributes can be accessed
directly using G.graph or added at instantiation using
keyword arguments.

>>> G=nx.Graph(region='Africa')
>>> G.graph['color']='green'
>>> G.graph
{'region': 'Africa', 'color': 'green'}

Node attributes
^^^^^^^^^^^^^^^
Each node has a corresponding dictionary of attributes.
Adding attributes to nodes is optional.

Add node attributes using add_node(), add_nodes_from() or G.node

>>> G.add_node(1, time='5pm')
>>> G.add_nodes_from([3], time='2pm')
>>> G.node[1]  # doctest: +SKIP
{'time': '5pm'}
>>> G.node[1]['room'] = 714  # doctest: +SKIP
>>> G.nodes(data=True)  # doctest: +SKIP
[(1, {'room': 714, 'time': '5pm'}), (3, {'time': '2pm'})]

Edge attributes
^^^^^^^^^^^^^^^
Each edge has a corresponding dictionary of attributes.
The default edge data is now an empty dictionary of attributes
and adding attributes to edges is optional.

A common use case is to add a weight attribute to an edge:

>>> G.add_edge(1,2,weight=3.14159)

Add edge attributes using add_edge(), add_edges_from(), subscript
notation, or G.edge.

>>> G.add_edge(1, 2, weight=4.7 )
>>> G.add_edges_from([(3,4),(4,5)], color='red')
>>> G.add_edges_from([(1,2,{'color':'blue'}), (2,3,{'weight':8})])
>>> G[1][2]['weight'] = 4.7
>>> G.edge[1][2]['weight'] = 4  # doctest: +SKIP

Methods changed
---------------

Graph(), DiGraph(), MultiGraph(), MultiDiGraph()
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   Now takes optional keyword=value attributes on initialization.

   >>> G=nx.Graph(year='2009',city='New York')

add_node()
^^^^^^^^^^
   Now takes optional keyword=value attributes or a dictionary of attributes.

   >>> G.add_node(1,room=714)


add_nodes_from()
^^^^^^^^^^^^^^^^
   Now takes optional keyword=value attributes or a dictionary of
   attributes applied to all affected nodes.

   >>> G.add_nodes_from([1,2],time='2pm')  # all nodes have same attribute

add_edge()
^^^^^^^^^^
   Now takes optional keyword=value attributes or a dictionary of attributes.

   >>> G.add_edge(1, 2, weight=4.7 )

add_edges_from()
^^^^^^^^^^^^^^^^
   Now takes optional keyword=value attributes or a dictionary of
   attributes applied to all affected edges.

   >>> G.add_edges_from([(3,4),(4,5)], color='red')
   >>> G.add_edges_from([(1,2,{'color':'blue'}), (2,3,{'weight':8})])


nodes() and nodes_iter()
^^^^^^^^^^^^^^^^^^^^^^^^
   New keyword data=True|False keyword determines whether to return
   two-tuples (n,dict) (True) with node attribution dictionary

   >>> G=nx.Graph([(1,2),(3,4)])
   >>> G.nodes(data=True)  # doctest: +SKIP
   [(1, {}), (2, {}), (3, {}), (4, {})]

copy()
^^^^^^
   Now returns a deep copy of the graph (copies all underlying
   data and attributes for nodes and edges).  Use the class
   initializer to make a shallow copy:

   >>> G=nx.Graph()
   >>> G_shallow=nx.Graph(G) # shallow copy
   >>> G_deep=G.copy() # deep copy

to_directed(), to_undirected()
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   Now returns a deep copy of the graph (copies all underlying
   data and attributes for nodes and edges).  Use the class
   initializer to make a shallow copy:

   >>> G = nx.Graph()
   >>> D_shallow = nx.DiGraph(G) # shallow copy
   >>> D_deep = G.to_directed() # deep copy

subgraph()
^^^^^^^^^^

   With copy=True now returns a deep copy of the graph
   (copies all underlying data and attributes for nodes and edges).

   >>> G = nx.Graph()
   >>> # note: copy keyword deprecated in networkx>1.0
   >>> # H = G.subgraph([],copy=True)  # deep copy of all data

add_cycle(), add_path(), add_star()
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   Now take optional keyword=value attributes or a dictionary of
   attributes which are applied to all edges affected by the method.

   >>> G = nx.Graph()
   >>> G.add_path([0, 1, 2, 3], width=3.2)  # doctest: +SKIP

Methods removed
---------------

delete_node()
^^^^^^^^^^^^^
   The preferred name is now remove_node().

delete_nodes_from()
^^^^^^^^^^^^^^^^^^^
   No longer raises an exception on an attempt to delete a node not in
   the graph.  The preferred name is now remove_nodes_from().

delete_edge()
^^^^^^^^^^^^^
   Now raises an exception on an attempt to delete an edge not in the graph.
   The preferred name is now remove_edge().

delete_edges_from()
^^^^^^^^^^^^^^^^^^^
   The preferred name is now remove_edges_from().

has_neighbor():

   Use has_edge()

get_edge()
^^^^^^^^^^
   Renamed to get_edge_data().	Returns the edge attribute dictionary.

   The fastest way to get edge data for edge (u,v) is to use G[u][v]
   instead of G.get_edge_data(u,v)


Members removed
---------------

directed, multigraph, weighted
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Use methods G.is_directed() and G.is_multigraph().
    All graphs are weighted graphs now if they have numeric
    values in the 'weight' edge attribute.


Methods added
-------------

add_weighted edges_from()
^^^^^^^^^^^^^^^^^^^^^^^^^
   Convenience method to add weighted edges to graph using a list of
   3-tuples (u,v,weight).

get_edge_data()
^^^^^^^^^^^^^^^
   Renamed from get_edge().

   The fastest way to get edge data for edge (u,v) is to use G[u][v]
   instead of G.get_edge_data(u,v)

is_directed()
^^^^^^^^^^^^^
    replaces member G.directed

is_multigraph()
^^^^^^^^^^^^^^^
    replaces member G.multigraph



Classes Removed
---------------

LabeledGraph, LabeledDiGraph
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    These classes have been folded into the regular classes.

UbiGraph
^^^^^^^^
    Removed as the ubigraph platform is no longer being supported.


Additional functions/generators
-------------------------------

ego_graph, stochastic_graph, PageRank algorithm, HITS algorithm,
GraphML writer, freeze, is_frozen, A* algorithm,
directed scale-free generator, random clustered graph.


Converting your existing code to networkx-1.0
---------------------------------------------

Weighted edges
^^^^^^^^^^^^^^

Edge information is now stored in an attribution dictionary
so all edge data must be given a key to identify it.

There is currently only one standard/reserved key, 'weight', which is
used by algorithms and functions that use weighted edges.  The
associated value should be numeric.  All other keys are available for
users to assign as needed.

>>> G=nx.Graph()
>>> G.add_edge(1,2,weight=3.1415) # add the edge 1-2 with a weight
>>> G[1][2]['weight']=2.3 # set the weight to 2.3

Similarly, for direct access the edge data, use
the key of the edge data to retrieve it.

>>> w = G[1][2]['weight']

All NetworkX algorithms that require/use weighted edges now use the
'weight' edge attribute.  If you have existing algorithms that assumed
the edge data was numeric, you should replace G[u][v] and
G.get_edge(u,v) with G[u][v]['weight'].

An idiom for getting a weight for graphs with or without an assigned
weight key is

>>> w= G[1][2].get('weight',1)  # set w to 1 if there is no 'weight' key
.. currentmodule:: networkx

Old Release Log
===============

NetworkX 2.5
------------
Release date: 22 August 2020

Supports Python 3.6, 3.7, and 3.8.

Release notes
~~~~~~~~~~~~~

See :doc:`release_2.5`.


NetworkX 2.4
------------
Release date: 16 October 2019

Supports Python 3.5, 3.6, 3.7, and 3.8.
This is the last release to support Python 3.5.

Release notes
~~~~~~~~~~~~~

See :doc:`release_2.4`.


NetworkX 2.3
------------
Release date: 11 April 2019

Supports Python 3.5, 3.6 and 3.7.
This is our first Python 3 only release.

Release notes
~~~~~~~~~~~~~

See :doc:`release_2.3`.


NetworkX 2.2
------------
Release date: 19 September 2018

Supports Python 2.7, 3.5, 3.6 and 3.7.
This is the last release to support Python 2.

Release notes
~~~~~~~~~~~~~

See :doc:`release_2.2`.


NetworkX 2.1
------------
Release date: 22 January 2018

Supports Python 2.7, 3.4, 3.5, and 3.6.

Release notes
~~~~~~~~~~~~~

See :doc:`release_2.1`.


NetworkX 2.0
------------
Release date: 20 September 2017

Support for Python 3.6 added, drop support for Python 3.3.

See :doc:`migration_guide_from_1.x_to_2.0`.

Release notes
~~~~~~~~~~~~~

See :doc:`release_2.0`.

NetworkX 1.11
-------------
Release date: 30 January 2016

Support for Python 3.5 added, drop support for Python 3.2.

Highlights
~~~~~~~~~~

Pydot features now use pydotplus.
Fixes installation on some machines and test with appveyor.
Restores default center and scale of layout routines.
Fixes various docs including no symbolic links in examples.
Docs can now build using autosummary on readthedocs.org.

NetworkX 1.10
--------------

Release date: 2 August 2015

Support for Python 2.6 is dropped in this release.

Highlights
~~~~~~~~~~

- Connected components now return generators
- new functions including

  + enumerate_all_cliques, greedy_coloring, edge_dfs, find_cycle
    immediate_dominators, harmonic_centrality
  + Hopcraft--Karp algorithm for maximum matchings
  + optimum branchings and arborescences.
  + all_simple_paths

- pyparsing dependence removed from GML reader/parser
- improve flow algorithms
- new generators related to expander graphs.
- new generators for multipartite graphs, nonisomorphic trees,
  circulant graphs
- allow graph subclasses to use dict-like objects in place of dicts
- added ordered graph subclasses
- pandas dataframe read/write added.
- data keyword in G.edges() allows requesting edge attribute directly
- expanded layout flexibility for node subsets
- Kanesky's algorithm for cut sets and k_components
- power function for graphs
- approximation of node connectivity
- transitive closure, triadic census and antichains
- quotient graphs and minors
- longest_path for DAGS
- modularity matrix routines

API changes
~~~~~~~~~~~
See :doc:`api_1.10`.

NetworkX 1.9.1
--------------
Release date: 13 September 2014

Bugfix release for minor installation and documentation issues.

NetworkX 1.9
------------
Release date: 21 June 2014

Support for Python 3.1 is dropped in this release.

Highlights
~~~~~~~~~~
- Completely rewritten maximum flow and flow-based connectivity algorithms with
  backwards incompatible interfaces
- Community graph generators
- Stoer–Wagner minimum cut algorithm
- Linear-time Eulerian circuit algorithm
- Linear algebra package changed to use SciPy sparse matrices
- Algebraic connectivity, Fiedler vector, spectral ordering algorithms
- Link prediction algorithms
- Goldberg–Radzik shortest path algorithm
- Semiconnected graph and tree recognition algorithms

API changes
~~~~~~~~~~~
See :doc:`api_1.9`.

NetworkX 1.8.1
--------------
Release date:  4 August 2013

Bugfix release for missing files in source packaging.


NetworkX 1.8
------------
Release date:  28 July 2013

Highlights
~~~~~~~~~~
- Faster (linear-time) graphicality tests and Havel-Hakimi graph generators
- Directed Laplacian matrix generator
- Katz centrality algorithm
- Functions to generate all simple paths
- Improved shapefile reader
- More flexible weighted projection of bipartite graphs
- Faster topological sort, descendants and ancestors of DAGs
- Scaling parameter for force-directed layout

Bug fixes
~~~~~~~~~
- Error with average weighted connectivity for digraphs, correct normalized laplacian with self-loops, load betweenness for single node graphs, isolated nodes missing from dfs/bfs trees, normalize HITS using l1, handle density of graphs with self loops

- Cleaner handling of current figure status with Matplotlib, Pajek files now don't write troublesome header line, default alpha value for GEXF files, read curved edges from yEd GraphML


For full details of the issues closed for this release (added features and bug fixes) see: https://github.com/networkx/networkx/issues?milestone=1&page=1&state=closed


API changes
~~~~~~~~~~~
See :doc:`api_1.8`


NetworkX 1.7
------------
Release date:  4 July 2012

Highlights
~~~~~~~~~~

- New functions for k-clique community finding, flow hierarchy,
  union, disjoint union, compose, and intersection operators that work on
  lists of graphs, and creating the biadjacency matrix of a bipartite graph.

- New approximation algorithms for dominating set, edge dominating set,
  independent set, max clique, and min-weighted vertex cover.

- Many bug fixes and other improvements.

For full details of the tickets closed for this release (added features and bug fixes) see:
https://networkx.lanl.gov/trac/query?status=closed&group=milestone&milestone=networkx-1.7

API changes
~~~~~~~~~~~
See :doc:`api_1.7`


NetworkX 1.6
------------

Release date:  20 November 2011

Highlights
~~~~~~~~~~

New functions for finding articulation points, generating random bipartite graphs, constructing adjacency matrix representations, forming graph products, computing assortativity coefficients, measuring subgraph centrality and communicability, finding k-clique communities, and writing JSON format output.

New examples for drawing with D3 Javascript library, and ordering matrices with the Cuthill-McKee algorithm.

More memory efficient implementation of current-flow betweenness and new approximation algorithms for current-flow betweenness and shortest-path betweenness.

Simplified handling of "weight" attributes for algorithms that use weights/costs/values.  See :doc:`api_1.6`.

Updated all code to work with the PyPy Python implementation http://pypy.org which produces faster performance on many algorithms.

For full details of the tickets closed for this release (added features and bug fixes) see:
https://networkx.lanl.gov/trac/query?status=closed&group=milestone&milestone=networkx-1.6

API changes
~~~~~~~~~~~
See :doc:`api_1.6`


NetworkX 1.5
------------

Release date:  4 June 2011

For full details of the tickets closed for this release see:
https://networkx.lanl.gov/trac/query?status=closed&group=milestone&milestone=networkx-1.5

Highlights
~~~~~~~~~~

New features
~~~~~~~~~~~~
 - Algorithms for :mod:`generating <networkx.generators.bipartite>`
   and :mod:`analyzing <networkx.algorithms.bipartite>` bipartite graphs
 - :mod:`Maximal independent set <networkx.algorithms.mis>` algorithm
 - :mod:`Erdős-Gallai graphical degree sequence test <networkx.generators.degree_seq>`
 - :mod:`Negative edge cycle test <networkx.algorithms.shortest_paths.weighted>`
 - More memory efficient :mod:`Dijkstra path length <networkx.algorithms.shortest_paths.weighted>` with cutoff parameter
 - :mod:`Weighted clustering coefficient <networkx.algorithms.cluster>`
 - Read and write version 1.2 of :mod:`GEXF reader <networkx.readwrite.gexf>` format
 - :mod:`Neighbor degree correlation <networkx.algorithms.neighbor_degree>`
   that handle subsets of nodes
 - :mod:`In-place node relabeling <networkx.relabel>`
 - Many 'weighted' graph algorithms now take optional parameter to use
   specified edge attribute (default='weight')
   (ticket https://networkx.lanl.gov/trac/ticket/509)

 - Test for :mod:`distance regular <networkx.algorithms.distance_regular>` graphs
 - Fast :mod:`directed Erdős-Renyi graph  <networkx.generators.random_graphs>` generator
 - Fast :mod:`expected degree graph  <networkx.generators.degree_seq>` generator
 - :mod:`Navigable small world  <networkx.generators.geometric>` generator
 - :mod:`Waxman model <networkx.generators.geometric>` generator
 - :mod:`Geographical threshold graph <networkx.generators.geometric>` generator
 - :mod:`Karate Club, Florentine Families, and Davis' Women's Club <networkx.generators.social>` graphs


API changes
~~~~~~~~~~~
See :doc:`api_1.5`


Bug fixes
~~~~~~~~~
 - Fix edge handling for multigraphs in networkx/graphviz interface
   (ticket https://networkx.lanl.gov/trac/ticket/507)
 - Update networkx/pydot interface for new versions of pydot
   (ticket https://networkx.lanl.gov/trac/ticket/506)
   (ticket https://networkx.lanl.gov/trac/ticket/535)
 - Fix negative cycle handling in Bellman-Ford
   (ticket https://networkx.lanl.gov/trac/ticket/502)
 - Write more attributes with GraphML and GML formats
   (ticket https://networkx.lanl.gov/trac/ticket/480)
 - Handle white space better in read_edgelist
   (ticket https://networkx.lanl.gov/trac/ticket/513)
 - Better parsing of Pajek format files
   (ticket https://networkx.lanl.gov/trac/ticket/524)
   (ticket https://networkx.lanl.gov/trac/ticket/542)
 - Isolates functions work with directed graphs
   (ticket https://networkx.lanl.gov/trac/ticket/526)
 - Faster conversion to numpy matrices
   (ticket https://networkx.lanl.gov/trac/ticket/529)
 - Add graph['name'] and use properties to access Graph.name
   (ticket https://networkx.lanl.gov/trac/ticket/544)
 - Topological sort confused None and 0
   (ticket https://networkx.lanl.gov/trac/ticket/546)
 - GEXF writer mishandled weight=0
   (ticket https://networkx.lanl.gov/trac/ticket/550)
 - Speedup in SciPy version of PageRank
   (ticket https://networkx.lanl.gov/trac/ticket/554)
 - Numpy PageRank node order incorrect + speedups
   (ticket https://networkx.lanl.gov/trac/ticket/555)

NetworkX 1.4
------------

Release date:  23 January 2011

New features
~~~~~~~~~~~~
 - :mod:`k-shell,k-crust,k-corona <networkx.algorithms.core>`
 - :mod:`read GraphML files from yEd <networkx.readwrite.graphml>`
 - :mod:`read/write GEXF format files <networkx.readwrite.gexf>`
 - :mod:`find cycles in a directed graph <networkx.algorithms.cycles>`
 - :mod:`DFS <networkx.algorithms.traversal.depth_first_search>` and :mod:`BFS <networkx.algorithms.traversal.breadth_first_search>` algorithms
 - :mod:`chordal graph functions <networkx.algorithms.chordal.chordal_alg>`
 - :mod:`Prim's algorithm for minimum spanning tree <networkx.algorithms.mst>`
 - :mod:`r-ary tree generator <networkx.generators.classic>`
 - :mod:`rich club coefficient <networkx.algorithms.richclub>`
 - NumPy matrix version of :mod:`Floyd's algorithm for all-pairs shortest path  <networkx.algorithms.shortest_paths.dense>`
 - :mod:`read GIS shapefiles <networkx.readwrite.nx_shp>`
 - :mod:`functions to get and set node and edge attributes <networkx.classes.function>`
 - and more, see  https://networkx.lanl.gov/trac/query?status=closed&group=milestone&milestone=networkx-1.4

API changes
~~~~~~~~~~~
 - :mod:`gnp_random_graph() <networkx.generators.random_graphs>` now takes a
   directed=True|False keyword instead of create_using
 - :mod:`gnm_random_graph() <networkx.generators.random_graphs>` now takes a
   directed=True|False keyword instead of create_using

Bug fixes
~~~~~~~~~
  - see  https://networkx.lanl.gov/trac/query?status=closed&group=milestone&milestone=networkx-1.4



NetworkX 1.3
------------

Release date:  28 August 2010

See: https://networkx.lanl.gov/trac/timeline

New features
~~~~~~~~~~~~
 - Works with Python versions 2.6, 2.7, 3.1, and 3.2 (but not 2.4 and 2.5).
 - :mod:`Minimum cost flow algorithms <networkx.algorithms.flow>`
 - :mod:`Bellman-Ford shortest paths <networkx.algorithms.shortest_paths.weighted>`
 - :mod:`GraphML reader and writer <networkx.readwrite.graphml>`
 - :mod:`More exception/error types <networkx.exception>`
 - Updated many tests to unittest style.  Run with: "import networkx; networkx.test()" (requires nose testing package)
 - and more, see  https://networkx.lanl.gov/trac/query?status=closed&group=milestone&milestone=networkx-1.3

API changes
~~~~~~~~~~~
 - :mod:`minimum_spanning_tree() now returns a NetworkX Graph (a tree or forest) <networkx.algorithms.mst>`

Bug fixes
~~~~~~~~~
  - see  https://networkx.lanl.gov/trac/query?status=closed&group=milestone&milestone=networkx-1.3


NetworkX 1.2
------------

Release date:  28 July 2010

See: https://networkx.lanl.gov/trac/timeline

New features
~~~~~~~~~~~~
 - :mod:`Ford-Fulkerson max flow and min cut <networkx.algorithms.flow>`
 - :mod:`Closeness vitality <networkx.algorithms.vitality>`
 - :mod:`Eulerian circuits <networkx.algorithms.euler>`
 - :mod:`Functions for isolates <networkx.algorithms.isolates>`
 - :mod:`Simpler s_max generator <networkx.generators.degree_seq>`
 - Compatible with IronPython-2.6
 - Improved testing functionality: import networkx; networkx.test() tests
   entire package and skips tests with missing optional packages
 - All tests work with Python-2.4
 - and more, see  https://networkx.lanl.gov/trac/query?status=closed&group=milestone&milestone=networkx-1.2


NetworkX 1.1
------------

Release date:  21 April 2010

See: https://networkx.lanl.gov/trac/timeline

New features
~~~~~~~~~~~~
 - :mod:`Algorithm for finding a basis for graph cycles <networkx.algorithms.cycles>`
 - :mod:`Blockmodeling <networkx.algorithms.block>`
 - :mod:`Assortativity and mixing matrices <networkx.algorithms.mixing>`
 - :mod:`in-degree and out-degree centrality <networkx.algorithms.centrality.degree>`
 - :mod:`Attracting components <networkx.algorithms.components.attracting>`
   and  :mod:`condensation <networkx.algorithms.components.strongly_connected>`.
 - :mod:`Weakly connected components <networkx.algorithms.components.weakly_connected>`
 - :mod:`Simpler interface to shortest path algorithms <networkx.algorithms.shortest_paths.generic>`
 - :mod:`Edgelist format to read and write data with attributes <networkx.readwrite.edgelist>`
 - :mod:`Attribute matrices <networkx.linalg.spectrum>`
 - :mod:`GML reader for nested attributes <networkx.readwrite.gml>`
 - Current-flow (random walk)
   :mod:`betweenness <networkx.algorithms.centrality.current_flow_betweenness>`
   and
   :mod:`closeness <networkx.algorithms.centrality.current_flow_closeness>`.
 - :mod:`Directed configuration model <networkx.generators.degree_seq>`,
   and  :mod:`directed random graph model <networkx.generators.random_graphs>`.
 - Improved documentation of drawing, shortest paths, and other algorithms
 - Many more tests, can be run with "import networkx; networkx.test()"
 - and much more, see  https://networkx.lanl.gov/trac/query?status=closed&group=milestone&milestone=networkx-1.1

API changes
~~~~~~~~~~~
Returning dictionaries
**********************
Several of the algorithms and the degree() method now return dictionaries
keyed by node instead of lists.  In some cases there was a with_labels
keyword which is no longer necessary.  For example,

>>> G=nx.Graph()
>>> G.add_edge('a','b')
>>> G.degree()  # doctest: +SKIP
{'a': 1, 'b': 1}

Asking for the degree of a single node still returns a single number

>>> G.degree('a')
1

The following now return dictionaries by default (instead of lists)
and the with_labels keyword has been removed:

 - :meth:`Graph.degree`,
   :meth:`MultiGraph.degree`,
   :meth:`DiGraph.degree`,
   :meth:`DiGraph.in_degree`,
   :meth:`DiGraph.out_degree`,
   :meth:`MultiDiGraph.degree`,
   :meth:`MultiDiGraph.in_degree`,
   :meth:`MultiDiGraph.out_degree`.
 - :func:`clustering`,
   :func:`triangles`
 - :func:`node_clique_number`,
   :func:`number_of_cliques`,
   :func:`cliques_containing_node`
 - :func:`eccentricity`


The following now return dictionaries by default (instead of lists)

 - :func:`pagerank`
 - :func:`hits`


Adding nodes
************
``add_nodes_from`` now accepts ``(node, attrdict)`` two-tuples

>>> G = nx.Graph()
>>> G.add_nodes_from([(1, {'color': 'red'})])

Examples
~~~~~~~~
 - `Mayvi2 drawing <https://networkx.org/documentation/latest/auto_examples/3d_drawing/mayavi2_spring.html>`_
 - `Blockmodel <https://networkx.org/documentation/latest/auto_examples/algorithms/plot_blockmodel.html>`_
 - `Sampson's monastery <https://networkx.org/documentation/latest/auto_examples/drawing/plot_sampson.html>`_
 - `Ego graph <https://networkx.org/documentation/latest/auto_examples/drawing/plot_ego_graph.html>`_

Bug fixes
~~~~~~~~~
 - Support graph attributes with union, intersection, and other graph operations
 - Improve subgraph speed (and related algorithms such as
   connected_components_subgraphs())
 - Handle multigraphs in more operators (e.g. union)
 - Handle double-quoted labels with pydot
 - Normalize betweenness_centrality for undirected graphs correctly
 - Normalize eigenvector_centrality by l2 norm
 - :func:`read_gml` now returns multigraphs

NetworkX 1.0.1
--------------

Release date:  11 Jan 2010

See: https://networkx.lanl.gov/trac/timeline

Bug fix release for missing setup.py in manifest.


NetworkX 1.0
------------

Release date:  8 Jan 2010

See: https://networkx.lanl.gov/trac/timeline

New features
~~~~~~~~~~~~
This release has significant changes to parts of the graph API
to allow graph, node, and edge attributes.
See http://networkx.lanl.gov/reference/api_changes.html

 - Update Graph, DiGraph, and MultiGraph classes to allow attributes.
 - Default edge data is now an empty dictionary (was the integer 1)
 - Difference and intersection operators
 - Average shortest path
 - A* (A-Star) algorithm
 - PageRank, HITS, and eigenvector centrality
 - Read Pajek files
 - Line graphs
 - Minimum spanning tree (Kruskal's algorithm)
 - Dense and sparse Fruchterman-Reingold layout
 - Random clustered graph generator
 - Directed scale-free graph generator
 - Faster random regular graph generator
 - Improved edge color and label drawing with Matplotlib
 - and much more, see  https://networkx.lanl.gov/trac/query?status=closed&group=milestone&milestone=networkx-1.0

Examples
~~~~~~~~
 - Update to work with networkx-1.0 API
 - Graph subclass example


NetworkX 0.99
-------------

Release date:  18 November 2008

See: https://networkx.lanl.gov/trac/timeline

New features
~~~~~~~~~~~~
This release has significant changes to parts of the graph API.
See http://networkx.lanl.gov/reference/api_changes.html

 - Update Graph and DiGraph classes to use weighted graphs as default
   Change in API for performance and code simplicity.
 - New MultiGraph and MultiDiGraph classes (replace XGraph and XDiGraph)
 - Update to use Sphinx documentation system http://networkx.lanl.gov/
 - Developer site at https://networkx.lanl.gov/trac/
 - Experimental LabeledGraph and LabeledDiGraph
 - Moved package and file layout to subdirectories.

Bug fixes
~~~~~~~~~
 - handle root= option to draw_graphviz correctly

Examples
~~~~~~~~
 - Update to work with networkx-0.99 API
 - Drawing examples now use matplotlib.pyplot interface
 - Improved drawings in many examples
 - New examples - see http://networkx.lanl.gov/examples/


NetworkX 0.37
---------------

Release date: 17 August 2008

See: https://networkx.lanl.gov/trac/timeline

NetworkX now requires Python 2.4 or later for full functionality.

New features
~~~~~~~~~~~~
 - Edge coloring and node line widths with Matplotlib drawings
 - Update pydot functions to work with pydot-1.0.2
 - Maximum-weight matching algorithm
 - Ubigraph interface for 3D OpenGL layout and drawing
 - Pajek graph file format reader and writer
 - p2g graph file format reader and writer
 - Secondary sort in topological sort

Bug fixes
~~~~~~~~~
 - Better edge data handling with GML writer
 - Edge betweenness fix for XGraph with default data of None
 - Handle Matplotlib version strings (allow "pre")
 - Interface to PyGraphviz (to_agraph()) now handles parallel edges
 - Fix bug in copy from XGraph to XGraph with multiedges
 - Use SciPy sparse lil matrix format instead of coo format
 - Clear up ambiguous cases for Barabasi-Albert model
 - Better care of color maps with Matplotlib when drawing colored nodes
   and edges
 - Fix error handling in layout.py

Examples
~~~~~~~~
 - Ubigraph examples showing 3D drawing


NetworkX 0.36
---------------

Release date: 13 January 2008

See: https://networkx.lanl.gov/trac/timeline


New features
~~~~~~~~~~~~
  - GML format graph reader, tests, and example (football.py)
  - edge_betweenness() and load_betweenness()

Bug fixes
~~~~~~~~~
  - remove obsolete parts of pygraphviz interface
  - improve handling of Matplotlib version strings
  - write_dot() now writes parallel edges and self loops
  - is_bipartite() and bipartite_color() fixes
  - configuration model speedup using random.shuffle()
  - convert with specified nodelist now works correctly
  - vf2 isomorphism checker updates

NetworkX 0.35.1
---------------

Release date: 27 July 2007

See: https://networkx.lanl.gov/trac/timeline

Small update to fix import readwrite problem and maintain Python2.3
compatibility.


NetworkX 0.35
-------------

Release date: 22 July 2007

See: https://networkx.lanl.gov/trac/timeline

New features
~~~~~~~~~~~~
  - algorithms for strongly connected components.
  - Brandes betweenness centrality algorithm (weighted and unweighted versions)
  - closeness centrality for weighted graphs
  - dfs_preorder, dfs_postorder, dfs_tree, dfs_successor, dfs_predecessor
  - readers for GraphML, LEDA, sparse6, and graph6 formats.
  - allow arguments in graphviz_layout to be passed directly to graphviz

Bug fixes
~~~~~~~~~
  - more detailed installation instructions
  - replaced dfs_preorder,dfs_postorder (see search.py)
  - allow initial node positions in spectral_layout
  - report no error on attempting to draw empty graph
  - report errors correctly when using tuples as nodes #114
  - handle conversions from incomplete dict-of-dict data

NetworkX 0.34
-------------

Release date: 12 April 2007

See: https://networkx.lanl.gov/trac/timeline

New features
~~~~~~~~~~~~
  - benchmarks for graph classes
  - Brandes betweenness centrality algorithm
  - Dijkstra predecessor and distance algorithm
  - xslt to convert DIA graphs to NetworkX
  - number_of_edges(u,v) counts edges between nodes u and v
  - run tests with python setup_egg.py test (needs setuptools)
    else use python -c "import networkx; networkx.test()"
  - is_isomorphic() that uses vf2 algorithm

Bug fixes
~~~~~~~~~
  - speedups of neighbors()
  - simplified Dijkstra's algorithm code
  - better exception handling for shortest paths
  - get_edge(u,v) returns None (instead of exception) if no edge u-v
  - floyd_warshall_array fixes for negative weights
  - bad G467, docs, and unittest fixes for graph atlas
  - don't put nans in numpy or scipy sparse adjacency matrix
  - handle get_edge() exception (return None if no edge)
  - remove extra kwds arguments in many places
  - no multi counting edges in conversion to dict of lists for multigraphs
  - allow passing tuple to get_edge()
  - bad parameter order in node/edge betweenness
  - edge betweenness doesn't fail with XGraph
  - don't throw exceptions for nodes not in graph (silently ignore instead)
    in edges_* and degree_*

NetworkX 0.33
-------------

Release date: 27 November 2006

See: https://networkx.lanl.gov/trac/timeline

New features
~~~~~~~~~~~~
  - draw edges with specified colormap
  - more efficient version of Floyd's algorithm for all pairs shortest path
  - use numpy only, Numeric is deprecated
  - include tests in source package (networkx/tests)
  - include documentation in source package (doc)
  - tests can now be run with
     >>> import networkx
     >>> networkx.test()  # doctest: +SKIP

Bug fixes
~~~~~~~~~
  - read_gpickle now works correctly with Windows
  - refactored large modules into smaller code files
  - degree(nbunch) now returns degrees in same order as nbunch
  - degree() now works for multiedges=True
  - update node_boundary and edge_boundary for efficiency
  - edited documentation for graph classes, now mostly in info.py

Examples
~~~~~~~~
  - Draw edges with colormap



NetworkX 0.32
-------------

Release date: 29 September 2006

See: https://networkx.lanl.gov/trac/timeline

New features
~~~~~~~~~~~~
  - Update to work with numpy-1.0x
  - Make egg usage optional: use python setup_egg.py bdist_egg to build egg
  - Generators and functions for bipartite graphs
  - Experimental classes for trees and forests
  - Support for new pygraphviz update (in nx_agraph.py) , see
    http://networkx.lanl.gov/pygraphviz/ for pygraphviz details

Bug fixes
~~~~~~~~~
  - Handle special cases correctly in triangles function
  - Typos in documentation
  - Handle special cases in shortest_path and shortest_path_length,
    allow cutoff parameter for maximum depth to search
  - Update examples: erdos_renyi.py, miles.py, roget,py, eigenvalues.py


Examples
~~~~~~~~
  - Expected degree sequence
  - New pygraphviz interface

NetworkX 0.31
-------------

Release date: 20 July 2006

See: https://networkx.lanl.gov/trac/timeline

New features
~~~~~~~~~~~~
   - arbitrary node relabeling (use relabel_nodes)
   - conversion of NetworkX graphs to/from Python dict/list types,
     numpy matrix or array types, and scipy_sparse_matrix types
   - generator for random graphs with given expected degree sequence

Bug fixes
~~~~~~~~~
   - Allow drawing graphs with no edges using pylab
   - Use faster heapq in dijkstra
   - Don't complain if X windows is not available

Examples
~~~~~~~~
   - update drawing examples


NetworkX 0.30
-------------


Release date: 23 June 2006

See: https://networkx.lanl.gov/trac/timeline

New features
~~~~~~~~~~~~
   - update to work with Python 2.5
   - bidirectional version of shortest_path and Dijkstra
   - single_source_shortest_path and all_pairs_shortest_path
   - s-metric and experimental code to generate  maximal s-metric graph
   - double_edge_swap and connected_double_edge_swap
   - Floyd's algorithm for all pairs shortest path
   - read and write unicode graph data to text files
   - read and write YAML format text files, http://yaml.org

Bug fixes
~~~~~~~~~
   - speed improvements (faster version of subgraph, is_connected)
   - added cumulative distribution and modified discrete distribution utilities
   - report error if DiGraphs are sent to connected_components routines
   - removed with_labels keywords for many functions where it was
     causing confusion
   - function name changes in shortest_path routines
   - saner internal handling of nbunch (node bunches), raise an
     exception if an nbunch isn't a node or iterable
   - better keyword handling in io.py allows reading multiple graphs
   - don't mix Numeric and numpy arrays in graph layouts and drawing
   - avoid automatically rescaling matplotlib axes when redrawing graph layout

Examples
~~~~~~~~
   - unicode node labels


NetworkX 0.29
-------------

Release date: 28 April 2006

See: https://networkx.lanl.gov/trac/timeline

New features
~~~~~~~~~~~~
   - Algorithms for betweenness, eigenvalues, eigenvectors, and
     spectral projection for threshold graphs
   - Use numpy when available
   - dense_gnm_random_graph generator
   - Generators for some directed graphs: GN, GNR, and GNC by Krapivsky
     and Redner
   - Grid graph generators now label by index tuples.  Helper
     functions for manipulating labels.
   - relabel_nodes_with_function


Bug fixes
~~~~~~~~~
   - Betweenness centrality now correctly uses Brandes definition and
     has normalization option outside main loop
   - Empty graph now labeled as empty_graph(n)
   - shortest_path_length used python2.4 generator feature
   - degree_sequence_tree off by one error caused nonconsecutive labeling
   - periodic_grid_2d_graph removed in favor of grid_2d_graph with
     periodic=True


NetworkX 0.28
-------------

Release date: 13 March 2006

See: https://networkx.lanl.gov/trac/timeline

New features
~~~~~~~~~~~~
  - Option to construct Laplacian with rows and columns in specified order
  - Option in convert_node_labels_to_integers to use sorted order
  - predecessor(G,n) function that returns dictionary of
    nodes with predecessors from breadth-first search of G
    starting at node n.
    https://networkx.lanl.gov/trac/ticket/26

Examples
~~~~~~~~
  - Formation of giant component in binomial_graph:
  - Chess masters matches:
  - Gallery https://networkx.org/documentation/latest/auto_examples/index.html

Bug fixes
~~~~~~~~~
  - Adjusted names for random graphs.
     + erdos_renyi_graph=binomial_graph=gnp_graph: n nodes with
       edge probability p
     + gnm_graph: n nodes and m edges
     + fast_gnp_random_graph: gnp for sparse graphs (small p)
  - Documentation contains correct spelling of Barabási, Bollobás,
    Erdős, and Rényi in UTF-8 encoding
  - Increased speed of connected_components and related functions
    by using faster BFS algorithm in networkx.paths
    https://networkx.lanl.gov/trac/ticket/27
  - XGraph and XDiGraph with multiedges=True produced error on delete_edge
  - Cleaned up docstring errors
  - Normalize names of some graphs to produce strings that represent
    calling sequence

NetworkX 0.27
-------------

Release date: 5 February 2006

See: https://networkx.lanl.gov/trac/timeline

New features
~~~~~~~~~~~~
  - sparse_binomial_graph: faster graph generator for sparse random graphs
  - read/write routines in io.py now handle XGraph() type and
    gzip and bzip2 files
  - optional mapping of type for read/write routine to allow
    on-the-fly conversion of node and edge datatype on read
  - Substantial changes related to digraphs and definitions of
    neighbors() and edges().  For digraphs edges=out_edges.
    Neighbors now returns a list of neighboring nodes with
    possible duplicates for graphs with parallel edges
    See https://networkx.lanl.gov/trac/ticket/24
  - Addition of out_edges, in_edges and corresponding out_neighbors
    and in_neighbors for digraphs.  For digraphs edges=out_edges.

Examples
~~~~~~~~
  - Minard's data for Napoleon's Russian campaign

Bug fixes
~~~~~~~~~
   - XGraph(multiedges=True) returns a copy of the list of edges
     for get_edge()

NetworkX 0.26
-------------

Release date: 6 January 2006

New features
~~~~~~~~~~~~
  - Simpler interface to drawing with pylab
  - G.info(node=None) function returns short information about graph
    or node
  - adj_matrix now takes optional nodelist to force ordering of
    rows/columns in matrix
  - optional pygraphviz and pydot interface to graphviz is now callable as
    "graphviz" with pygraphviz preferred.  Use draw_graphviz(G).

Examples
~~~~~~~~
  - Several new examples showing how draw to graphs with various
    properties of nodes, edges, and labels

Bug fixes
~~~~~~~~~
   - Default data type for all graphs is now None (was the integer 1)
   - add_nodes_from now won't delete edges if nodes added already exist
   - Added missing names to generated graphs
   - Indexes for nodes in graphs start at zero by default (was 1)

NetworkX 0.25
-------------

Release date: 5 December 2005

New features
~~~~~~~~~~~~
  - Uses setuptools for installation http://peak.telecommunity.com/DevCenter/setuptools
  - Improved testing infrastructure, can now run python setup.py test
  - Added interface to draw graphs with pygraphviz
    https://networkx.lanl.gov/pygraphviz/
  - is_directed() function call

Examples
~~~~~~~~
  - Email example shows how to use XDiGraph with Python objects as
    edge data


Documentation
~~~~~~~~~~~~~
  - Reformat menu, minor changes to Readme, better stylesheet

Bug fixes
~~~~~~~~~
   - use create_using= instead of result= keywords for graph types
     in all cases
   - missing weights for degree 0 and 1 nodes in clustering
   - configuration model now uses XGraph, returns graph with identical
     degree sequence as input sequence
   - fixed Dijkstra priority queue
   - fixed non-recursive toposort and is_directed_acyclic graph

NetworkX 0.24
-------------

Release date: 20 August 2005

Bug fixes
~~~~~~~~~
   - Update of Dijkstra algorithm code
   - dfs_successor now calls proper search method
   - Changed to list comprehension in DiGraph.reverse() for python2.3
     compatibility
   - Barabasi-Albert graph generator fixed
   - Attempt to add self loop should add node even if parallel edges not
     allowed

NetworkX 0.23
-------------

Release date: 14 July 2005

The NetworkX web locations have changed:

http://networkx.lanl.gov/     - main documentation site
http://networkx.lanl.gov/svn/  - subversion source code repository
https://networkx.lanl.gov/trac/ - bug tracking and info


Important Change
~~~~~~~~~~~~~~~~
The naming conventions in NetworkX have changed.
The package name "NX" is now "networkx".

The suggested ways to import the NetworkX package are

 - import networkx
 - import networkx as NX
 - from networkx import *

New features
~~~~~~~~~~~~
  - DiGraph reverse
  - Graph generators
     + watts_strogatz_graph now does rewiring method
     + old watts_strogatz_graph->newman_watts_strogatz_graph

Examples
~~~~~~~~

Documentation
~~~~~~~~~~~~~
  - Changed to reflect NX-networkx change
  - main site is now https://networkx.lanl.gov/

Bug fixes
~~~~~~~~~
   - Fixed logic in io.py for reading DiGraphs.
   - Path based centrality measures (betweenness, closeness)
     modified so they work on graphs that are not connected and
     produce the same result as if each connected component were
     considered separately.

NetworkX 0.22
-------------

Release date: 17 June 2005

New features
~~~~~~~~~~~~
  - Topological sort, testing for directed acyclic graphs (DAGs)
  - Dijkstra's algorithm for shortest paths in weighted graphs
  - Multidimensional layout with dim=n for drawing
  - 3d rendering demonstration with vtk
  - Graph generators
     + random_powerlaw_tree
     + dorogovtsev_goltsev_mendes_graph

Examples
~~~~~~~~
  - Kevin Bacon movie actor graph: Examples/kevin_bacon.py
  - Compute eigenvalues of graph Laplacian: Examples/eigenvalues.py
  - Atlas of small graphs: Examples/atlas.py

Documentation
~~~~~~~~~~~~~
  - Rewrite of setup scripts to install documentation and
    tests in documentation directory specified

Bug fixes
~~~~~~~~~
   - Handle calls to edges() with non-node, non-iterable items.
   - truncated_tetrahedral_graph was just plain wrong
   - Speedup of betweenness_centrality code
   - bfs_path_length now returns correct lengths
   - Catch error if target of search not in connected component of source
   - Code cleanup to label internal functions with _name
   - Changed import statement lines to always use "import NX" to
     protect name-spaces
   - Other minor bug-fixes and testing added
:orphan:

*****************************
Preparing for the 3.0 release
*****************************

.. note::
   Much of the work leading to the NetworkX 3.0 release will be included
   in the NetworkX 2.6 and 2.7 releases.  For example, we are deprecating a lot
   of old code in the 2.6 and 2.7 releases.  This guide will discuss this
   ongoing work and will help you understand what changes you can make now
   to minimize the disruption caused by the move to 3.0.

This is a guide for people moving from NetworkX 2.X to NetworkX 3.0

Any issues with these can be discussed on the `mailing list
<https://groups.google.com/forum/#!forum/networkx-discuss>`_.

The focus of 3.0 release is on addressing years of technical debt, modernizing our codebase,
improving performance, and making it easier to contribute.
We plan to release 2.7 near the end of summer and 3.0 near the end of the year.

Default dependencies
--------------------

We no longer depend on the "decorator" library.

Deprecated code
---------------

The 2.6 release deprecates over 30 functions.
See :ref:`networkx_2.6`.
Next Release
============

Release date: TBD

Supports Python ...

NetworkX is a Python package for the creation, manipulation, and study of the
structure, dynamics, and functions of complex networks.

For more information, please visit our `website <https://networkx.org/>`_
and our :ref:`gallery of examples <examples_gallery>`.
Please send comments and questions to the `networkx-discuss mailing list
<http://groups.google.com/group/networkx-discuss>`_.

Highlights
----------

This release is the result of X of work with over X pull requests by
X contributors. Highlights include:


Improvements
------------


API Changes
-----------


Deprecations
------------


Merged PRs
----------

<output of contribs.py>


Contributors
------------

<output of contribs.py>
NetworkX 1.7
============

Release date:  4 July 2012

Highlights
~~~~~~~~~~

- New functions for k-clique community finding, flow hierarchy,
  union, disjoint union, compose, and intersection operators that work on
  lists of graphs, and creating the biadjacency matrix of a bipartite graph.

- New approximation algorithms for dominating set, edge dominating set,
  independent set, max clique, and min-weighted vertex cover.

- Many bug fixes and other improvements.

For full details of the tickets closed for this release (added features and bug fixes) see:
https://networkx.lanl.gov/trac/query?status=closed&group=milestone&milestone=networkx-1.7

Other
-----
* Untested bipartite_random_regular_graph() removed.

:orphan:

*******************************
Migration guide from 1.X to 2.0
*******************************

This is a guide for people moving from NetworkX 1.X to NetworkX 2.0

Any issues with these can be discussed on the `mailing list
<https://groups.google.com/forum/#!forum/networkx-discuss>`_.

At the bottom of this document we discuss how to create code that will
work with both NetworkX v1.x and v2.0.

We have made some major changes to the methods in the Multi/Di/Graph classes.
The methods changed are explained with examples below.

With the release of NetworkX 2.0 we are moving to a view/iterator reporting API.
We have changed many methods from reporting lists or dicts to iterating over
the information. Most of the changes in this regard are in the base classes.
Methods that used to return containers now return views (inspired from
`dictionary views <https://docs.python.org/3/library/stdtypes.html#dict-views>`_
in Python) and methods that returned an iterator have been removed.
The methods which create new graphs have changed in the depth of data copying.
``G.subgraph``/``edge_subgraph``/``reverse``/``to_directed``/``to_undirected``
are affected.  Many now have options for view creation instead of copying data.
The depth of the data copying may have also changed.

One View example is ``G.nodes`` (or ``G.nodes()``) which now returns a
dict-like NodeView while ``G.nodes_iter()`` has been removed. Similarly
for views with ``G.edges`` and removing ``G.edges_iter``.
The Graph attributes ``G.node`` and ``G.edge`` have been removed in favor of
using ``G.nodes[n]`` and ``G.edges[u, v]``.
Finally, the ``selfloop`` methods and ``add_path``/``star``/``cycle`` have
been moved from graph methods to networkx functions.

We expect that these changes will break some code. We have tried to make
them break the code in ways that raise exceptions, so it will be obvious
that code is broken.

There are also a number of improvements to the codebase outside of the base
graph classes. These are too numerous to catalog here, but a couple obvious
ones include:

- centering of nodes in ``drawing/nx_pylab``,
- iterator vs dict output from a few ``shortest_path`` routines

-------

Some demonstrations:

    >>> import networkx as nx
    >>> G = nx.complete_graph(5)
    >>> G.nodes  # for backward compatibility G.nodes() works as well
    NodeView((0, 1, 2, 3, 4))

You can iterate through ``G.nodes`` (or ``G.nodes()``)

    >>> for node in G.nodes:
    ...     print(node)
    0
    1
    2
    3
    4

If you want a list of nodes you can use the Python list function

    >>> list(G.nodes)
    [0, 1, 2, 3, 4]

``G.nodes`` is set-like allowing set operations. It is also dict-like in that you
can look up node data with ``G.nodes[n]['weight']``. You can still use the calling
interface ``G.nodes(data='weight')`` to iterate over node/data pairs. In addition
to the dict-like views ``keys``/``values``/``items``, ``G.nodes`` has a data-view
G.nodes.data('weight').  The new EdgeView ``G.edges`` has similar features for edges.

By adding views NetworkX supports some new features like set operations on
views.

    >>> H = nx.Graph()
    >>> H.add_nodes_from([1, 'networkx', '2.0'])
    >>> G.nodes & H.nodes  # finding common nodes in 2 graphs
    {1}
    >>> # union of nodes in 2 graphs
    >>> G.nodes | H.nodes  # doctest: +SKIP
    {0, 1, 2, 3, 4, 'networkx', '2.0'}

Similarly, ``G.edges`` now returns an EdgeView instead of a list of edges and it
also supports set operations.

    >>> G.edges  # for backward compatibility G.nodes() works as well
    EdgeView([(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)])
    >>> list(G.edges)
    [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]

``G.degree`` now returns a DegreeView. This is less dict-like than the other views
in the sense that it iterates over (node, degree) pairs, does not provide
keys/values/items/get methods. It does provide lookup ``G.degree[n]`` and
``(node, degree)`` iteration. A dict keyed by nodes to degree values can be
easily created if needed as ``dict(G.degree)``.

    >>> G.degree  # for backward compatibility G.degree() works as well
    DegreeView({0: 4, 1: 4, 2: 4, 3: 4, 4: 4})
    >>> G.degree([1, 2, 3])
    DegreeView({1: 4, 2: 4, 3: 4})
    >>> list(G.degree([1, 2, 3]))
    [(1, 4), (2, 4), (3, 4)]
    >>> dict(G.degree([1, 2, 3]))
    {1: 4, 2: 4, 3: 4}
    >>> G.degree
    DegreeView({0: 4, 1: 4, 2: 4, 3: 4, 4: 4})
    >>> list(G.degree)
    [(0, 4), (1, 4), (2, 4), (3, 4), (4, 4)]
    >>> dict(G.degree)
    {0: 4, 1: 4, 2: 4, 3: 4, 4: 4}

The degree of an individual node can be calculated by ``G.degree[node]``.
Similar changes have been made to ``in_degree`` and ``out_degree``
for directed graphs. If you want just the degree values, here are some options.
They are shown for ``in_degree`` of a ``DiGraph``, but similar ideas work
for ``out_degree`` and ``degree``

    >>> DG = nx.DiGraph()
    >>> DG.add_weighted_edges_from([(1, 2, 0.5), (3, 1, 0.75)])
    >>> deg = DG.in_degree   # sets up the view
    >>> [d for n, d in deg]   # gets all nodes' degree values
    [1, 1, 0]
    >>> (d for n, d in deg)    # iterator over degree values
    <generator object <genexpr> ...>
    >>> [deg[n] for n in [1, 3]]   # using lookup for only some nodes
    [1, 0]

    >>> for node, in_deg in dict(DG.in_degree).items():  # works for nx1 and nx2
    ...     print(node, in_deg)
    1 1
    2 1
    3 0
    >>> dict(DG.in_degree([1, 3])).values()    # works for nx1 and nx2
    dict_values([1, 0])
    >>> # DG.in_degree(nlist) creates a restricted view for only nodes in nlist.
    >>> # but see the fourth option above for using lookup instead.
    >>> list(d for n, d in DG.in_degree([1, 3]))
    [1, 0]

    >>> [len(nbrs) for n, nbrs in DG.pred.items()]  # probably slightly fastest for all nodes
    [1, 1, 0]
    >>> [len(DG.pred[n]) for n in [1, 3]]           # probably slightly faster for only some nodes
    [1, 0]

-------

If ``n`` is a node in ``G``, then ``G.neighbors(n)`` returns an iterator.

    >>> n = 1
    >>> G.neighbors(n)
    <dict_keyiterator object at ...>
    >>> list(G.neighbors(n))
    [0, 2, 3, 4]

DiGraphViews behave similar to GraphViews, but have a few more methods.

    >>> D = nx.DiGraph()
    >>> D.add_edges_from([(1, 2), (2, 3), (1, 3), (2, 4)])
    >>> D.nodes
    NodeView((1, 2, 3, 4))
    >>> list(D.nodes)
    [1, 2, 3, 4]
    >>> D.edges
    OutEdgeView([(1, 2), (1, 3), (2, 3), (2, 4)])
    >>> list(D.edges)
    [(1, 2), (1, 3), (2, 3), (2, 4)]
    >>> D.in_degree[2]
    1
    >>> D.out_degree[2]
    2
    >>> D.in_edges
    InEdgeView([(1, 2), (2, 3), (1, 3), (2, 4)])
    >>> list(D.in_edges())
    [(1, 2), (2, 3), (1, 3), (2, 4)]
    >>> D.out_edges(2)
    OutEdgeDataView([(2, 3), (2, 4)])
    >>> list(D.out_edges(2))
    [(2, 3), (2, 4)]
    >>> D.in_degree
    InDegreeView({1: 0, 2: 1, 3: 2, 4: 1})
    >>> list(D.in_degree)
    [(1, 0), (2, 1), (3, 2), (4, 1)]
    >>> D.successors(2)
    <dict_keyiterator object at ...>
    >>> list(D.successors(2))
    [3, 4]
    >>> D.predecessors(2)
    <dict_keyiterator object at ...>
    >>> list(D.predecessors(2))
    [1]

The same changes apply to MultiGraphs and MultiDiGraphs.

-------

The order of arguments to ``set_edge_attributes`` and ``set_node_attributes`` has
changed.  The position of ``name`` and ``values`` has been swapped, and ``name`` now
defaults to ``None``.  The previous call signature of ``(graph, name, value)`` has
been changed to ``(graph, value, name=None)``. The new style allows for ``name`` to
be omitted in favor of passing a dictionary of dictionaries to ``values``.

A simple method for migrating existing code to the new version is to explicitly
specify the keyword argument names. This method is backwards compatible and
ensures the correct arguments are passed, regardless of the order. For example the old code

    >>> G = nx.Graph([(1, 2), (1, 3)])
    >>> nx.set_node_attributes(G, 'label', {1: 'one', 2: 'two', 3: 'three'})  # doctest: +SKIP
    >>> nx.set_edge_attributes(G, 'label', {(1, 2): 'path1', (2, 3): 'path2'})  # doctest: +SKIP

Will cause ``TypeError: unhashable type: 'dict'`` in the new version. The code
can be refactored as

    >>> G = nx.Graph([(1, 2), (1, 3)])
    >>> nx.set_node_attributes(G, name='label', values={1: 'one', 2: 'two', 3: 'three'})
    >>> nx.set_edge_attributes(G, name='label', values={(1, 2): 'path1', (2, 3): 'path2'})

-------

Some methods have been moved from the base graph class into the main namespace.
These are:  ``G.add_path``, ``G.add_star``, ``G.add_cycle``, ``G.number_of_selfloops``,
``G.nodes_with_selfloops``, and ``G.selfloop_edges``.
They are replaced by ``nx.path_graph(G, ...)`` ``nx.add_star(G, ...)``,
``nx.selfloop_edges(G)``, etc.
For backward compatibility, we are leaving them as deprecated methods.

-------

With the new GraphViews (SubGraph, ReversedGraph, etc) you can't assume that
``G.__class__()`` will create a new instance of the same graph type as ``G``.
In fact, the call signature for ``__class__`` differs depending on whether ``G``
is a view or a base class. For v2.x you should use ``G.fresh_copy()`` to
create a null graph of the correct type---ready to fill with nodes and edges.

Graph views can also be views-of-views-of-views-of-graphs. If you want to find the
original graph at the end of this chain use ``G.root_graph``. Be careful though
because it may be a different graph type (directed/undirected) than the view.

-------

``topological_sort``  no longer accepts ``reverse`` or ``nbunch`` arguments.
If ``nbunch`` was a single node source, then the same effect can now be achieved
using the ``subgraph`` operator:

    nx.topological_sort(G.subgraph(nx.descendants(G, nbunch)))

To achieve a reverse topological sort, the output should be converted to a list:

    reversed(list(nx.topological_sort(G)))

-------

Writing code that works for both versions
=========================================

Methods ``set_node_attributes``/``get_node_attributes``/``set_edge_attributes``/``get_edge_attributes``
have changed the order of their keyword arguments ``name`` and ``values``. So, to make it
work with both versions you should use the keywords in your call.

    >>> nx.set_node_attributes(G, values=1.0, name='weight')

-------

Change any method with ``_iter`` in its name to the version without ``_iter``.
In v1 this replaces an iterator by a list, but the code will still work.
In v2 this creates a view (which acts like an iterator).

-------

Replace any use of ``G.edge`` with ``G.adj``. The Graph attribute ``edge``
has been removed. The attribute ``G.adj`` is ``G.edge`` in v1 and will work
with both versions.

-------

If you use ``G.node.items()`` or similar in v1.x, you can replace it with
``G.nodes(data=True)`` which works for v2.x and v1.x.  Iterating over ``G.node```
as in ``for n in G.node:`` can be replaced with ``G``, as in: ``for n in G:``.

-------

The Graph attribute ``node`` has moved its functionality to ``G.nodes``, so code
expected to work with v2.x should use ``G.nodes``.
In fact most uses of ``G.node`` can be replaced by an idiom that works for both
versions. The functionality that can't easily is: ``G.node[n]``.
In v2.x that becomes ``G.nodes[n]`` which doesn't work in v1.x.

Luckily you can still use ``G.node[n]`` in v2.x when you want it to be able to work
with v1.x too. We have left ``G.node`` in v2.x as a transition pointer to ``G.nodes``.
We envision removing ``G.node`` in v3.x (sometime in the future).

-------

Copying node attribute dicts directly from one graph to another can corrupt
the node data structure if not done correctly. Code such as the following:

    >>> # dangerous in v1.x, not allowed in v2.x
    >>> G.node[n] = H.node[n]  # doctest: +SKIP

used to work, even though it could cause errors if ``n`` was not a node in ``G``.
That code will cause an error in v2.x.  Replace it with one of the more safe versions:

    >>> G.nodes[n].update(H.nodes[n])  # works in v2.x

-------

The methods removed from the graph classes and put into the main package namespace
can be used via the associated deprecated methods. If you want to update your code
to the new functions, one hack to make that work with both versions is to write
your code for v2.x and add code to the v1 namespace in an ad hoc manner:

    >>> if nx.__version__[0] == '1':
    ...     nx.add_path = lambda G, nodes: G.add_path(nodes)

Similarly, v2.x code that uses ``G.fresh_copy()`` or ``G.root_graph`` is hard to make
work for v1.x. It may be best in this case to determine the graph type you want
explicitly and call Graph/DiGraph/MultiGraph/MultiDiGraph directly.

Using Pickle with v1 and v2
===========================

The Pickle protocol does not store class methods, only the data. So if you write a
pickle file with v1 you should not expect to read it into a v2 Graph. If this happens
to you, read it in with v1 installed and write a file with the node and edge
information. You can read that into a config with v2 installed and then add those nodes
and edges to a fresh graph. Try something similar to this:

    >>> # in v1.x
    >>> pickle.dump([G.nodes(data=True), G.edges(data=True)], file)  # doctest: +SKIP
    >>> # then in v2.x
    >>> nodes, edges = pickle.load(file)  # doctest: +SKIP
    >>> G = nx.Graph()  # doctest: +SKIP
    >>> G.add_nodes_from(nodes)  # doctest: +SKIP
    >>> G.add_edges_from(edges)  # doctest: +SKIP
NetworkX 1.10
=============

Release date: 2 August 2015

Support for Python 2.6 is dropped in this release.

Highlights
~~~~~~~~~~

- Connected components now return generators
- new functions including

  + enumerate_all_cliques, greedy_coloring, edge_dfs, find_cycle
    immediate_dominators, harmonic_centrality
  + Hopcraft--Karp algorithm for maximum matchings
  + optimum branchings and arborescences.
  + all_simple_paths

- pyparsing dependence removed from GML reader/parser
- improve flow algorithms
- new generators related to expander graphs.
- new generators for multipartite graphs, nonisomorphic trees,
  circulant graphs
- allow graph subclasses to use dict-like objects in place of dicts
- added ordered graph subclasses
- pandas dataframe read/write added.
- data keyword in G.edges() allows requesting edge attribute directly
- expanded layout flexibility for node subsets
- Kanesky's algorithm for cut sets and k_components
- power function for graphs
- approximation of node connectivity
- transitive closure, triadic census and antichains
- quotient graphs and minors
- longest_path for DAGS
- modularity matrix routines

API changes
-----------
* [`#1501 <https://github.com/networkx/networkx/pull/1501>`_]
  ``connected_components``, ``weakly_connected_components``, and
  ``strongly_connected_components`` return now a generator of sets of
  nodes. Previously the generator was of lists of nodes. This PR also
  refactored the ``connected_components`` and ``weakly_connected_components``
  implementations making them faster, especially for large graphs.

* [`#1547 <https://github.com/networkx/networkx/pull/1547>`_]
  The ``func_iter`` functions in Di/Multi/Graphs classes are slated for
  removal in NetworkX 2.0 release. ``func`` will behave like ``func_iter``
  and return an iterator instead of list. These functions are deprecated in
  NetworkX 1.10 release.

New functionalities
-------------------

* [`#823 <https://github.com/networkx/networkx/pull/823>`_]
  A ``enumerate_all_cliques`` function is added in the clique package
  (``networkx.algorithms.clique``) for enumerating all cliques (including
  nonmaximal ones) of undirected graphs.

* [`#1105 <https://github.com/networkx/networkx/pull/1105>`_]
  A coloring package (``networkx.algorithms.coloring``) is created for
  graph coloring algorithms. Initially, a ``greedy_color`` function is
  provided for coloring graphs using various greedy heuristics.

* [`#1193 <https://github.com/networkx/networkx/pull/1193>`_]
  A new generator ``edge_dfs``, added to ``networkx.algorithms.traversal``,
  implements a depth-first traversal of the edges in a graph. This complements
  functionality provided by a depth-first traversal of the nodes in a graph.
  For multigraphs, it allows the user to know precisely which edges were
  followed in a traversal. All NetworkX graph types are supported. A traversal
  can also reverse edge orientations or ignore them.

* [`#1194 <https://github.com/networkx/networkx/pull/1194>`_]
  A ``find_cycle`` function is added to the ``networkx.algorithms.cycles``
  package to find a cycle in a graph. Edge orientations can be optionally
  reversed or ignored.

* [`#1210 <https://github.com/networkx/networkx/pull/1210>`_]
  Add a random generator for the duplication-divergence model.

* [`#1241 <https://github.com/networkx/networkx/pull/1241>`_]
  A new ``networkx.algorithms.dominance`` package is added for
  dominance/dominator algorithms on directed graphs. It contains a
  ``immediate_dominators`` function for computing immediate
  dominators/dominator trees and a ``dominance_frontiers`` function for
  computing dominance frontiers.

* [`#1269 <https://github.com/networkx/networkx/pull/1269>`_]
  The GML reader/parser and writer/generator are rewritten to remove the
  dependence on pyparsing and enable handling of arbitrary graph data.

* [`#1280 <https://github.com/networkx/networkx/pull/1280>`_]
  The network simplex method in the ``networkx.algorithms.flow`` package is
  rewritten to improve its performance and support multi- and disconnected
  networks. For some cases, the new implementation is two or three orders of
  magnitude faster than the old implementation.

* [`#1286 <https://github.com/networkx/networkx/pull/1286>`_]
  Added the Margulis--Gabber--Galil graph to ``networkx.generators``.

* [`#1306 <https://github.com/networkx/networkx/pull/1306>`_]
  Added the chordal p-cycle graph, a mildly explicit algebraic construction
  of a family of 3-regular expander graphs. Also, moves both the existing
  expander graph generator function (for the Margulis-Gabber-Galil
  expander) and the new chordal cycle graph function to a new module,
  ``networkx.generators.expanders``.

* [`#1314 <https://github.com/networkx/networkx/pull/1314>`_]
  Allow overwriting of base class dict with dict-like:
  OrderedGraph, ThinGraph, PrintGraph, etc.

* [`#1321 <https://github.com/networkx/networkx/pull/1321>`_]
  Added ``to_pandas_dataframe`` and ``from_pandas_dataframe``.

* [`#1322 <https://github.com/networkx/networkx/pull/1322>`_]
  Added the Hopcroft--Karp algorithm for finding a maximum cardinality
  matching in bipartite graphs.

* [`#1336 <https://github.com/networkx/networkx/pull/1336>`_]
  Expanded data keyword in G.edges and added default keyword.

* [`#1338 <https://github.com/networkx/networkx/pull/1338>`_]
  Added support for finding optimum branchings and arborescences.

* [`#1340 <https://github.com/networkx/networkx/pull/1340>`_]
  Added a ``from_pandas_dataframe`` function that accepts Pandas DataFrames
  and returns a new graph object. At a minimum, the DataFrame must have two
  columns, which define the nodes that make up an edge. However, the function
  can also process an arbitrary number of additional columns as edge
  attributes, such as 'weight'.

* [`#1354 <https://github.com/networkx/networkx/pull/1354>`_]
  Expanded layout functions to add flexibility for drawing subsets of nodes
  with distinct layouts and for centering each layout around given
  coordinates.

* [`#1356 <https://github.com/networkx/networkx/pull/1356>`_]
  Added ordered variants of default graph class.

* [`#1360 <https://github.com/networkx/networkx/pull/1360>`_]
  Added harmonic centrality to ``network.algorithms.centrality``.

* [`#1390 <https://github.com/networkx/networkx/pull/1390>`_]
  The ``generators.bipartite`` have been moved to
  ``algorithms.bipartite.generators``. The functions are not imported in the
  main  namespace, so to use it, the bipartite package has to be imported.

* [`#1391 <https://github.com/networkx/networkx/pull/1391>`_]
  Added Kanevsky's algorithm for finding all minimum-size separating
  node sets in an undirected graph. It is implemented as a generator
  of node cut sets.

* [`#1399 <https://github.com/networkx/networkx/pull/1399>`_]
  Added power function for simple graphs

* [`#1405 <https://github.com/networkx/networkx/pull/1405>`_]
  Added fast approximation for node connectivity based on White and
  Newman's approximation algorithm for finding node independent paths
  between two nodes.

* [`#1413 <https://github.com/networkx/networkx/pull/1413>`_]
  Added transitive closure and antichains function for directed acyclic
  graphs in ``algorithms.dag``. The antichains function was contributed
  by Peter Jipsen and Franco Saliola and originally developed for the
  SAGE project.

* [`#1425 <https://github.com/networkx/networkx/pull/1425>`_]
  Added generator function for the complete multipartite graph.

* [`#1427 <https://github.com/networkx/networkx/pull/1427>`_]
  Added nonisomorphic trees generator.

* [`#1436 <https://github.com/networkx/networkx/pull/1436>`_]
  Added a generator function for circulant graphs to the
  ``networkx.generators.classic`` module.

* [`#1437 <https://github.com/networkx/networkx/pull/1437>`_]
  Added function for computing quotient graphs; also created a new module,
  ``networkx.algorithms.minors``.

* [`#1438 <https://github.com/networkx/networkx/pull/1438>`_]
  Added longest_path and longest_path_length for DAG.

* [`#1439 <https://github.com/networkx/networkx/pull/1439>`_]
  Added node and edge contraction functions to
  ``networkx.algorithms.minors``.

* [`#1445 <https://github.com/networkx/networkx/pull/1448>`_]
  Added a new modularity matrix module to ``networkx.linalg``,
  and associated spectrum functions to the ``networkx.linalg.spectrum``
  module.

* [`#1447 <https://github.com/networkx/networkx/pull/1447>`_]
  Added function to generate all simple paths starting with the shortest
  ones based on Yen's algorithm for finding k shortest paths at
  ``algorithms.simple_paths``.

* [`#1455 <https://github.com/networkx/networkx/pull/1455>`_]
  Added the directed modularity matrix to the
  ``networkx.linalg.modularity_matrix`` module.

* [`#1474 <https://github.com/networkx/networkx/pull/1474>`_]
  Adds ``triadic_census`` function; also creates a new module,
  ``networkx.algorithms.triads``.

* [`#1476 <https://github.com/networkx/networkx/pull/1476>`_]
  Adds functions for testing if a graph has weighted or negatively weighted
  edges. Also adds a function for testing if a graph is empty. These are
  ``is_weighted``, ``is_negatively_weighted``, and ``is_empty``.

* [`#1481 <https://github.com/networkx/networkx/pull/1481>`_]
  Added Johnson's algorithm; one more algorithm for shortest paths. It
  solves all pairs shortest path problem. This is ``johnson`` at
  ``algorithms.shortest_paths``

* [`#1414 <https://github.com/networkx/networkx/pull/1414>`_]
  Added Moody and White algorithm for identifying ``k_components`` in a
  graph, which is based on Kanevsky's algorithm for finding all minimum-size
  node cut-sets (implemented in ``all_node_cuts`` #1391).

* [`#1415 <https://github.com/networkx/networkx/pull/1415>`_]
  Added fast approximation for ``k_components`` to the
  ``networkx.approximation`` package. This is based on White and Newman
  approximation algorithm for finding node independent paths between two
  nodes (see #1405).

Removed functionalities
-----------------------

* [`#1236 <https://github.com/networkx/networkx/pull/1236>`_]
  The legacy ``ford_fulkerson`` maximum flow function is removed. Use
  ``edmonds_karp`` instead.

Miscellaneous changes
---------------------

* [`#1192 <https://github.com/networkx/networkx/pull/1192>`_]
  Support for Python 2.6 is dropped.
.. _networkx_2.6:

NetworkX 2.6
============

Release date: 08 July 2021

Supports Python 3.7, 3.8, and 3.9.

This release has a larger than normal number of changes in preparation for the upcoming 3.0 release.
The current plan is to release 2.7 near the end of summer and 3.0 in late 2021.
See :doc:`migration_guide_from_2.x_to_3.0` for more details.

NetworkX is a Python package for the creation, manipulation, and study of the
structure, dynamics, and functions of complex networks.

For more information, please visit our `website <https://networkx.org/>`_
and our :ref:`gallery of examples <examples_gallery>`
Please send comments and questions to the `networkx-discuss mailing list
<http://groups.google.com/group/networkx-discuss>`_.

Highlights
----------

This release is the result of 11 months of work with over 363 pull requests by
91 contributors. Highlights include:

- Dropped support for Python 3.6
- Dropped "decorator" library dependency
- Improved example gallery
- Removed code for supporting Jython/IronPython
- The ``__str__`` method for graph objects is more informative and concise.
- Improved import time
- Improved test coverage
- New documentation theme
- Add functionality for drawing self-loop edges
- Add approximation algorithms for Traveling Salesman Problem

New functions:

- Panther algorithm
- maximum cut heuristics
- equivalence_classes
- dedensification
- random_ordered_tree
- forest_str
- snap_aggregation
- networkx.approximation.diameter
- partition_quality
- prominent_group
- prefix_tree_recursive
- topological_generations

NXEPs
-----

**N**\etwork\ **X** **E**\nhancement **P**\roposals capture changes
that are larger in scope than typical pull requests, such as changes to
fundamental data structures.
The following proposals have come under consideration since the previous
release:

- :ref:`NXEP2`
- :ref:`NXEP3`

Improvements
------------

- [`#3886 <https://github.com/networkx/networkx/pull/3886>`_]
  Adds the Panther algorithm for top-k similarity search.
- [`#4138 <https://github.com/networkx/networkx/pull/4138>`_]
  Adds heuristics for approximating solution to the maximum cut problem.
- [`#4183 <https://github.com/networkx/networkx/pull/4183>`_]
  Adds ``equivalence_classes`` to public API.
- [`#4193 <https://github.com/networkx/networkx/pull/4193>`_]
  ``nx.info`` is more concise.
- [`#4198 <https://github.com/networkx/networkx/pull/4198>`_]
  Improve performance of ``transitivity``.
- [`#4206 <https://github.com/networkx/networkx/pull/4206>`_]
  UnionFind.union selects the heaviest root as the new root
- [`#4240 <https://github.com/networkx/networkx/pull/4240>`_]
  Adds ``dedensification`` function in a new ``summarization`` module.
- [`#4294 <https://github.com/networkx/networkx/pull/4294>`_]
  Adds ``forest_str`` for string representation of trees.
- [`#4319 <https://github.com/networkx/networkx/pull/4319>`_]
  pagerank uses scipy by default now.
- [`#4841 <https://github.com/networkx/networkx/pull/4841>`_]
  simrank_similarity uses numpy by default now.
- [`#4317 <https://github.com/networkx/networkx/pull/4317>`_]
  New ``source`` argument to ``has_eulerian_path`` to look for path starting at
  source.
- [`#4356 <https://github.com/networkx/networkx/pull/4356>`_]
  Use ``bidirectional_djikstra`` in ``shortest_path`` for weighted graphs
  to improve performance.
- [`#4361 <https://github.com/networkx/networkx/pull/4361>`_]
  Adds ``nodelist`` argument to ``triadic_census``
- [`#4435 <https://github.com/networkx/networkx/pull/4435>`_]
  Improve ``group_betweenness_centrality``.
- [`#4446 <https://github.com/networkx/networkx/pull/4446>`_]
  Add ``sources`` parameter to allow computing ``harmonic_centrality`` from a
  subset of nodes.
- [`#4463 <https://github.com/networkx/networkx/pull/4463>`_]
  Adds the ``snap`` summarization algorithm.
- [`#4476 <https://github.com/networkx/networkx/pull/4476>`_]
  Adds the ``diameter`` function for approximating the lower bound on the
  diameter of a graph.
- [`#4519 <https://github.com/networkx/networkx/pull/4519>`_]
  Handle negative weights in clustering algorithms.
- [`#4528 <https://github.com/networkx/networkx/pull/4528>`_]
  Improved performance of ``edge_boundary``.
- [`#4560 <https://github.com/networkx/networkx/pull/4560>`_]
  Adds ``prominent_group`` function to find prominent group of size k in
  G according to group_betweenness_centrality.
- [`#4588 <https://github.com/networkx/networkx/pull/4588>`_]
  Graph intersection now works when input graphs don't have the same node sets.
- [`#4607 <https://github.com/networkx/networkx/pull/4607>`_]
  Adds approximation algorithms for solving the traveling salesman problem,
  including ``christofides``, ``greedy_tsp``, ``simulated_annealing_tsp``,
  and ``threshold_accepting_tsp``.
- [`#4640 <https://github.com/networkx/networkx/pull/4640>`_]
  ``prefix_tree`` now uses a non-recursive algorithm. The original recursive
  algorithm is still available via ``prefix_tree_recursive``.
- [`#4659 <https://github.com/networkx/networkx/pull/4659>`_]
  New ``initial_graph`` argument to ``barabasi_albert_graph`` and
  ``dual_barabasi_albert_graph`` to supply an initial graph to the model.
- [`#4690 <https://github.com/networkx/networkx/pull/4690>`_]
  ``modularity_max`` now supports edge weights.
- [`#4727 <https://github.com/networkx/networkx/pull/4727>`_]
  Improved performance of ``scale_free_graph``.
- [`#4739 <https://github.com/networkx/networkx/pull/4739>`_]
  Added `argmap` function to replace the decorator library dependence
- [`#4757 <https://github.com/networkx/networkx/pull/4757>`_]
  Adds ``topological_generations`` function for DAG stratification.
- [`#4768 <https://github.com/networkx/networkx/pull/4768>`_]
  Improved reproducibility of geometric graph generators.
- [`#4769 <https://github.com/networkx/networkx/pull/4769>`_]
  Adds ``margins`` keyword to ``draw_networkx_nodes`` to control node clipping
  in images with large node sizes.
- [`#4812 <https://github.com/networkx/networkx/pull/4812>`_]
  Use ``scipy`` implementation for ``hits`` algorithm to improve performance.
- [`#4847 <https://github.com/networkx/networkx/pull/4847>`_]
  Improve performance of ``scipy`` implementation of ``hits`` algorithm.

API Changes
-----------

- [`#4183 <https://github.com/networkx/networkx/pull/4183>`_]
  ``partition`` argument of `quotient_graph` now accepts dicts
- [`#4190 <https://github.com/networkx/networkx/pull/4190>`_]
  Removed ``tracemin_chol``.  Use ``tracemin_lu`` instead.
- [`#4216 <https://github.com/networkx/networkx/pull/4216>`_]
  In `to_*_array/matrix`, nodes in nodelist but not in G now raise an exception.
  Use G.add_nodes_from(nodelist) to add them to G before converting.
- [`#4360  <https://github.com/networkx/networkx/pull/4360>`_]
  Internally `.nx_pylab.draw_networkx_edges` now always generates a
  list of `matplotlib.patches.FancyArrowPatch` rather than using
  a `matplotlib.collections.LineCollection` for un-directed graphs.  This
  unifies interface for all types of graphs.  In
  addition to the API change this may cause a performance regression for
  large graphs.
- [`#4384 <https://github.com/networkx/networkx/pull/4384>`_]
  Added ``edge_key`` parameter for MultiGraphs in to_pandas_edgelist
- [`#4461 <https://github.com/networkx/networkx/pull/4461>`_]
  Added ``create_using`` parameter to ``binomial_tree``
- [`#4466 <https://github.com/networkx/networkx/pull/4466>`_]
  `relabel_nodes` used to raise a KeyError for a key in `mapping` that is not
  a node in the graph, but it only did this when `copy` was `False`. Now
  any keys in `mapping` which are not in the graph are ignored.
- [`#4502 <https://github.com/networkx/networkx/pull/4502>`_]
  Moves ``maximum_independent_set`` to the ``clique`` module in ``approximation``.
- [`#4536 <https://github.com/networkx/networkx/pull/4536>`_]
  Deprecate ``performance`` and ``coverage`` in favor of ``partition_quality``,
  which computes both metrics simultaneously and is more efficient.
- [`#4573 <https://github.com/networkx/networkx/pull/4573>`_]
  `label_propagation_communities` returns a `dict_values` object of community
  sets of nodes instead of a generator of community sets. It is still iterable,
  so likely will still work in most user code and a simple fix otherwise:
  e.g., add ``iter( ... )`` surrounding the function call.
- [`#4545 <https://github.com/networkx/networkx/pull/4545>`_]
  `prefix_tree` used to return `tree, root` but root is now always 0
  instead of a UUID generate string. So the function returns `tree`.
- [`#4545 <https://github.com/networkx/networkx/pull/4545>`_]
  The variable `NIL` ="NIL" has been removed from `networkx.generators.trees`
- [`#3620 <https://github.com/networkx/networkx/pull/3620>`_]
  The function `naive_greedy_modularity_communities` now returns a
  list of communities (like `greedy_modularity_communities`) instead
  of a generator of communities.
- [`#4786 <https://github.com/networkx/networkx/pull/4786>`_]
  Deprecate the ``attrs`` keyword argument in favor of explicit keyword
  arguments in the ``json_graph`` module.
- [`#4843 <https://github.com/networkx/networkx/pull/4843>`_]
  The unused ``normalized`` parameter has been removed
  from ``communicability_betweeness_centrality``
- [`#4850 <https://github.com/networkx/networkx/pull/4850>`_]
  Added ``dtype`` parameter to adjacency_matrix
- [`#4851 <https://github.com/networkx/networkx/pull/4851>`_]
  Output of `numeric_mixing_matrix` and `degree_mixing_matrix` no longer
  includes rows with all entries zero by default. The functions now accept
  a parameter `mapping` keyed by value to row index to identify each row.
- [`#4867 <https://github.com/networkx/networkx/pull/4867>`_]
  The function ``spring_layout`` now ignores 'fixed' nodes not in the graph

Deprecations
------------

- [`#4238 <https://github.com/networkx/networkx/pull/4238>`_]
  Deprecate ``to_numpy_matrix`` and ``from_numpy_matrix``.
- [`#4279 <https://github.com/networkx/networkx/pull/4279>`_]
  Deprecate ``networkx.utils.misc.is_iterator``.
  Use ``isinstance(obj, collections.abc.Iterator)`` instead.
- [`#4280 <https://github.com/networkx/networkx/pull/4280>`_]
  Deprecate ``networkx.utils.misc.is_list_of_ints`` as it is no longer used.
  See ``networkx.utils.misc.make_list_of_ints`` for related functionality.
- [`#4281 <https://github.com/networkx/networkx/pull/4281>`_]
  Deprecate ``read_yaml`` and ``write_yaml``.
- [`#4282 <https://github.com/networkx/networkx/pull/4282>`_]
  Deprecate ``read_gpickle`` and ``write_gpickle``.
- [`#4298 <https://github.com/networkx/networkx/pull/4298>`_]
  Deprecate ``read_shp``, ``edges_from_line``, and ``write_shp``.
- [`#4319 <https://github.com/networkx/networkx/pull/4319>`_]
  Deprecate ``pagerank_numpy``, ``pagerank_scipy``.
- [`#4355 <https://github.com/networkx/networkx/pull/4355>`_]
  Deprecate ``copy`` method in the coreview Filtered-related classes.
- [`#4384 <https://github.com/networkx/networkx/pull/4384>`_]
  Deprecate unused ``order`` parameter in to_pandas_edgelist.
- [`#4428 <https://github.com/networkx/networkx/pull/4428>`_]
  Deprecate ``jit_data`` and ``jit_graph``.
- [`#4449 <https://github.com/networkx/networkx/pull/4449>`_]
  Deprecate ``consume``.
- [`#4448 <https://github.com/networkx/networkx/pull/4448>`_]
  Deprecate ``iterable``.
- [`#4536 <https://github.com/networkx/networkx/pull/4536>`_]
  Deprecate ``performance`` and ``coverage`` in favor of ``parition_quality``.
- [`#4545 <https://github.com/networkx/networkx/pull/4545>`_]
  Deprecate ``generate_unique_node``.
- [`#4599 <https://github.com/networkx/networkx/pull/4599>`_]
  Deprecate ``empty_generator``.
- [`#4600 <https://github.com/networkx/networkx/pull/4600>`_]
  Deprecate ``default_opener``.
- [`#4617 <https://github.com/networkx/networkx/pull/4617>`_]
  Deprecate ``hub_matrix`` and ``authority_matrix``
- [`#4629 <https://github.com/networkx/networkx/pull/4629>`_]
  Deprecate the ``Ordered`` graph classes.
- [`#4802 <https://github.com/networkx/networkx/pull/4802>`_]
  The ``nx_yaml`` function has been removed along with the dependency on
  ``pyyaml``. Removal implemented via module ``__getattr__`` to patch security
  warnings related to ``pyyaml.Loader``.
- [`#4826 <https://github.com/networkx/networkx/pull/4826>`_]
  Deprecate ``preserve_random_state``.
- [`#4827 <https://github.com/networkx/networkx/pull/4827>`_]
  Deprecate ``almost_equal``.
- [`#4833 <https://github.com/networkx/networkx/pull/4833>`_]
  Deprecate ``run``.
- [`#4829 <https://github.com/networkx/networkx/pull/4829>`_]
  Deprecate ``assert_nodes_equal``, ``assert_edges_equal``, and ``assert_graphs_equal``.
- [`#4850 <https://github.com/networkx/networkx/pull/4850>`_]
  Deprecate ``adj_matrix``.
- [`#4841 <https://github.com/networkx/networkx/pull/4841>`_]
  Deprecate ``simrank_similarity_numpy``.
- [`#4923 <https://github.com/networkx/networkx/pull/4923>`_]
  Deprecate ``numeric_mixing_matrix``.
- [`#4937 <https://github.com/networkx/networkx/pull/4937>`_]
  Deprecate ``k_nearest_neighbors``.

Merged PRs
----------

- Bump release version
- Update release process
- Update website doc
- fix issue #4173: cytoscape_graph(input_data) did modify the original data (#4176)
- Some docstring fixes for draw_networkx_edge_labels() in nx_pylab.py + one typo (#4182)
- TST: add dtype to pandas test (#4185)
- Partitions for quotient graphs (#4183)
- graphml: re-add graph attribute type 'long' after 857aa81 removed it (#4189)
- Test mac osx via actions (#4201)
- DOC: Update docstrings in cytoscape module (#4180)
- rewrite add_nodes_from to relax code meant to allow ironpython pre-2.7.5 (#4200)
- Speed up transitivity, remove redundant call (#4198)
- NXEP 2 — API design of view slices (#4101)
- Cleanup old platforms (#4202)
- Fixed "topolgical_sort" typo (#4211)
- Make optional dependencies default on CPython
- Simplify imports
- Populate setup.py requires from requirements
- Update dependencies
- Remove _CholeskySolver
- to_numpy/scipy array functions should not allow non-nodes in nodelist (#4216)
- fix "see also" links in json_graph.tree (#4222)
- MAINT: changed is_string_like to isinstance (#4223)
- Fix UnionFind.union to select the heaviest root as the new root (#4206)
- CI: Configure circleCI to deploy docs. (#4134)
- MAINT: Update nx.info (#4193)
- Fix indexing in kernighan_lin_bisection (#4177)
- CI: Add GH fingerprint (#4229)
- Create ssh dir for circleci
- CI: update circleci doc deployment. (#4230)
- Revert "CI: Configure circleCI to deploy docs. (#4134)" (#4231)
- DOC: Add discussion to NXEP 2.
- Update format dependencies
- Use black for linting
- Format w/ black==20.8b1
- Check formatting of PRs via black (#4235)
- TST: Modify heuristic for astar path test. (#4237)
- MAINT: Deprecate numpy matrix conversion functions (#4238)
- Add roadmap (#4234)
- Add nx.info to str dunder for graph classes (#4241)
- DOC: Minor reformatting of contract_nodes docstring. (#4245)
- Fix betweenness_centrality doc paper links (#4257)
- Fix bug in has_eulerian_path for directed graphs  (#4246)
- Add PR template (#4258)
- Use seed to make plot fixed (#4260)
- Update giant component example (#4267)
- Update "house with colors" gallery example (#4263)
- Replace degree_histogram and degree_rank with a single example (#4265)
- Update Knuth miles example. (#4251)
- Update "four_grids" gallery example (#4264)
- Improve legibility of labels in plot_labels_and_colors example (#4266)
- Improve readibility of chess_example in gallery (#4252)
- Fix contracted_edge for multiple edges (#4274)
- Add seeds to gallery examples for reproducibility (#4276)
- Add a 3D plotting example with matplotlib to the gallery (#4268)
- Deprecate `utils.is_iterator` (#4279)
- Deprecate utils.is_list_of_ints (#4280)
- Improve axes layout in plot_decomposition example (#4278)
- Update homepage URL (#4285)
- Build docs for deployment on Travis CI (#4286)
- Add simple graph w/ manual layout (#4291)
- Deprecate nx_yaml (#4281)
- Deprecate gpickle (#4282)
- Improve relabel coverage, tweak docstrings (#4299)
- Swith to travis-ci.com
- TST: Increase test coverage of convert_matrix (#4301)
- Add descriptive error message for Node/EdgeView slicing. NEXP2 (#4300)
- Don't import other people's version.py (#4289)
- TST: Refactor to improve coverage. (#4307)
- Improve readwrite test coverage (#4310)
- Fix typo (#4312)
- Update docstring of to_dict_of_dicts.
- Add tests for edge_data param.
- Minor touchups to docstring
- adds dedensification function (#4240)
- TST: improve multigraph test coverage to 100% (#4340)
- Add rainbow coloring example to gallery. (#4330)
- Test on Python 3.9 (#4303)
- Sphinx33 (#4342)
- fix order of yield and seen.update in all cc routines (see #4331 & #3859 & 3823) (#4333)
- Updates to slicing error message for reportviews (#4304)
- Eulerian path fix (#4317)
- Add FutureWarning in preparation for simplifying cytoscape function signatures. (#4284)
- Move a few imports inside functions to improve import speed of the library (#4296)
- Address comments from code review.
- Cleanup algebraicconnectivity (#4287)
- Switch from travis to gh actions (#4320)
- Fix (#4345)
- Fix travis doc deployment
- Fix gdal version on travis
- Update to_dict_of_dict edge_data (#4321)
- Update adjacency_iter to adjacency (#4339)
- Test and document missing nodes/edges in set_{node/edge}_attributes (#4346)
- Update tests and docs for has_eulerian_path (#4344)
- Deprecate nx_shp (#4298)
- Refactor and improve test coverage for restricted_view and selfloop_edges (#4351)
- Enable mayavi in sphinx gallery. (#4297)
- CI: Add mayavi conf to travis and GH for doc deploy (#4354)
- Fix doc build w/ GH actions
- Install vtk before mayavi
- Install vtk before mayavi
- Install vtk before mayavi
- Use bidirection_dijkstra as default in weighted shortest_path (#4356)
- Add unit tests for utils.misc.flatten (#4359)
- Improve test coverage for coreviews.py (#4355)
- Update tutorial.rst - Fixes #4249 (#4358)
- Bugfix for issue 4336, moving try/except and adding else clause (#4365)
- Added nodelist attribute to triadic_census (#4361)
- API: always use list of FancyArrowPatch rather than LineCollection (#4360)
- MNT: make the self-loop appear in all cases (#4370)
- Add additional libraries to intersphinx mapping (#4372)
- Make nx.pagerank a wrapper around different implementations, use scipy one by default (#4319)
- MAINT: remove deprecated numpy type aliases. (#4373)
- DOC: Fix return type for random_tournament and hamiltonian_path (#4376)
- Skip memory leak test for PyPy (#4385)
- add OSMnx example (#4383)
- Update docstring for to_pandas_edgelist and add edgekey parameter (#4384)
- TST: Boost test coverage of nx_pylab module (#4375)
- Fixed issue where edge attributes were being silently overwritten during node contraction (#4273)
- CI: Fix CircleCI doc buiild failure (#4388)
- Improve test coverage of convert module (#4306)
- Add gene-gene network (#4269)
- Ignore expected warnings (#4391)
- Use matrix multiplication operator (#4390)
- code and doc fix for square_clustering algorithm in cluster.py (#4392)
- Remove xml import checks (#4393)
- fix typo in NXEP template (#4396)
- Add Panther algorithm per #3849 (#3886)
- Pagerank followup (#4399)
- Don't import nx from networkx (#4403)
- Modify and document behavior of nodelist param in draw_networkx_edges. (#4378)
- Add circuit plot (#4408)
- Add words graph plot (#4409)
- DOC: Remove repeated words (#4410)
- Add plot for rcm example (#4411)
- Fix small index iteration bug in kernighan_lin algorithm (#4398)
- Use str dunder (#4412)
- Use xetex for uft8 latex backend (#4326)
- Add recommended fonts to travis.yml. (#4414)
- CI: Workaround font naming bug. (#4416)
- DOC: geospatial example using lines (#4407)
- Add plotting examples for geospatial data (#4366)
- Increase coverage in graphviews.py (#4418)
- Refactor gallery (#4422)
- Safer repr format of variables (#4413)
- Updates to docs and imports for classic.py (#4424)
- Remove advanced example section (#4429)
- Add coreview objects to documentation (#4431)
- Add gallery example for drawing self-loops. (#4430)
- Add igraph example (#4404)
- Standard imports (#4401)
- Collect graphviz examples (#4427)
- NXEP 3: Allow generators to yield from edgelists (#4395)
- Update geospatial readme (#4417)
- DOC: Fix broken links in shortest_path docstrings (#4434)
- Improves description bfs_predecessors and bfs_successors. (#4438)
- Deprecate jit (#4428)
- JavaScript example: fix link (#4450)
- Deprecate utils.misc.consume (#4449)
- DOC: Switch from napoleon to numpydoc sphinx extension (#4447)
- Correct networkxsimplex docstring re: multigraph
- Correct networkxsimplex docstring re: multigraph (#4455)
- Maxcut heuristics (#4138)
- binomial_tree() with "create_using parameter (#4461)
- Reorganize tests (#4467)
- Drop Py3.6 support per NEP 29 (#4469)
- Add random_ordered_tree and forest_str (#4294)
- Deprecate iterable (#4448)
- Allow relabel_nodes mapping to have non-node keys that get ignored (#4466)
- Fixed docs + added decorator for k_components approx (#4474)
- Update docs for clustering Fixes #4348 (#4477)
- Handle self-loops for single self-loop (drawing) (#4425)
- Update GH actions links in README (#4482)
- Improve code coverage for cuts.py (#4473)
- Reenable tests (#4488)
- Update Sphinx (#4494)
- Update pre-commit (#4495)
- Simplify example dependencies (#4506)
- Update geospatial readme (#4504)
- Update year (#4509)
- Drop Travis CI (#4510)
- Run pypy tests separately (#4512)
- Simplify version information (#4492)
- Delete old test (#4513)
- Gallery support for pygraphviz examples (#4464)
- TST: An approach to parametrizing read_edgelist tests. (#4292)
- Setup cross-repo doc deploy via actions. (#4480)
- use issue templates to redirect to discussions tab, add a bug report template (#4524)
- Fix performance issue in nx.edge_boundary (#4528)
- clean up list comp (#4499)
- Improve code coverage of swap.py (#4529)
- Clustering for signed weighted graphs (#4519)
- Fix docstrings and remove unused variables (#4501)
- Improving code coverage of chordal.py (#4471)
- Cliques on mutigraph/directed graph types (#4502)
- Approximated Diameter  (#4476)
- `arrows` should be True by default for directed graphs (#4522)
- Remove unnecessary node_list from gallery example (#4505)
- fixing the width argument description of the function draw_networkx (#4479)
- Partially revert #4378 - Modify behavior of nodelist param in draw_networkx_edges. (#4531)
- Replace generate_unique_node internally where not needed (#4537)
- Extend harmonic centrality to include source nodes (#4446)
- improve group betweenness centrality (#4435)
- fixes Github Actions failures (#4548)
- updated cutoff def in weighted.py (#4546)
- Less strict on mayavi constraint for doc building. (#4547)
- Update docstring for ancestor and descendents (#4550)
- TST: Fix error in katz centrality test setup. (#4554)
- Correct mu parameter documentation for LFR (#4557)
- Pin pygeos==0.8 (#4563)
- Unpin pygeos (#4570)
- Test Windows via GH actions (#4567)
- Update documentation and testing of arbitrary_element (#4451)
- added test for max_iter argument
- reformatted test_kernighan_lin.py
- Simplify test pylab (#4577)
- Update README.rst
- Fix search (#4580)
- Add test Kernighan Lin Algorithm (#4575)
- Fix typos (#4581)
- Boiler plate for mentored projects documentation (#4576)
- Deprecate generate_unique_node (#4545)
- Check nodelist input to floyd_warshall (#4589)
- Improve intersection function (#4588)
- Pygraphviz choco (#4583)
- Add prominent group algorithm (#4560)
- Add partition_quality to compute coverage and performance  (coverage and perfor… (#4536)
- Use Pillow for viewing AGraph output and deprecate default_opener (#4600)
- Remove mktemp usage (#4593)
- Add an FAQ to the developer guide for new contributors (#4556)
- Improve test coverage and docs for nonrandomness (#4613)
- Collect label propagation communities in one go (#4573)
- Deprecate networkx.utils.empty_generator. (#4599)
- return earlier from `clique.graph_clique_number` (#4622)
- More for projects page: TSP and Graph Isomorphism (#4620)
- add recommended venv directory to .gitignore (#4619)
- adding weight description to centrality metrices (#4610)
- Add a good first issue badge to README  (#4627)
- add test to regular (#4624)
- Add scipy-1.6.1 to blocklist. (#4628)
- Deprecate hub_matrix and authority_matrix (#4617)
- Fix issue #3153: generalized modularity maximization  (#3260)
- Improve doc example for find_cycle. (#4639)
- Correct and update Atlas example (#4635)
- Remove attr_dict from parameters list in the docstring (#4642)
- Verify edges are valid in is_matching() (#4638)
- Remove old file reference (#4646)
- Deprecate Ordered graph classes (#4629)
- Update CI to use main (#4651)
- Make main default branch (and remove gitwash) (#4649)
- Fix link for Katz centrality definition (#4655)
- fix for negative_edge_cycle weight kwarg to bellman_ford (#4658)
- Refactor bipartite and multipartite layout (#4653)
- Volunteering for mentorship (#4671)
- Adding an iterative version of prefix tree (#4640)
- Increase code coverage tournament (#4665)
- Fix to_vertex_cover (#4667)
- Reorganize minor submodule as subpackage (#4349)
- modularity_max: account for edge weights (#4690)
- Remove instances of random.sample from sets (deprecated in Python 3.9) (#4602)
- Fixing Bug in Transitive Reduction, resulting in loss of node/edge attributes (#4684)
- direct links to the tutorial and discussions in README (#4711)
- Pin upper bound of decorator dep. (#4721)
- fix typo (#4724)
- Updating average_clustering() documentation - Issue #4734 (#4735)
- rm nx import from docstring example. (#4738)
- CI: persist pip cache between circleci runs (#4714)
- Use pydata sphinx theme (#4741)
- O(n^2) -> O(n) implementation for scale_free_graph (#4727)
- TST: be more explicit about instance comparison. (#4748)
- fix typo in docstring (ismorphism -> isomorphism) (#4756)
- CI: Fix cartopy build failure in docs workflow (#4751)
- Add missing __all__'s to utils modules + test. (#4753)
- Add 2 articles for TSP project as references (#4758)
- Improve reproducibilty of geometric graphs (#4768)
- Updated decorator requirement for #4718 (#4773)
- Gallery Example: Drawing custom node icons on network using MPL (#4633)
- Get rid of invalid escape sequences. (#4789)
- imread(url) is deprecated, use pillow + urllib to load image from URL (#4790)
- Add auto-margin scaling in draw_networkx_nodes function (fix for issue 3443) (#4769)
- Update documentation dependencies (#4794)
- Fix sphinx warnings during doc build. (#4795)
- Remove mayavi and cartopy dependencies (#4800)
- make plots less dense, enable plotting for igraph (#4791)
- fix urllib import (#4793)
- Improve documentation look (#4801)
- Add approximation algorithms for traveling salesman problem (#4607)
- adds implementation of SNAP summarization algorithm (#4463)
- Update black (#4814)
- Restructure documentation (#4744)
- Pin upper bound on decorator for 2.6 release. (#4815)
- Use `callable()` to check if the object is calllable (#1) (#4678)
- Remove dictionary from signature of tree_graph and tree_data (#4786)
- Make nx.hits a wrapper around different implementations, use scipy one by default (#4812)
- restructured networksimplex.py and added test_networksimplex.py (#4685)
- Update requirements (#4625)
- Fix Sphinx errors (#4817)
- Add topological_generations function (#4757)
- Add `initial_graph` parameter to simple and dual Barábasi-Albert random graphs (#4659)
- Link to guides (#4818)
- switch alias direction of spring_layout and fruchterman_reingold_layout (#4820)
- Fix to_undirected doc typo (#4821)
- Deprecate preserve_random_state (#4826)
- Fixes read/write_gml with nan/inf attributes (#4497)
- Remove pyyaml dependency via module getattr (#4802)
- Use pytest.approx (#4827)
- DOC: Clarify behaviour of k_crust(G, k) (#4831)
- Limit number of threads used by OMP in circleci. (#4830)
- Deprecate run (#4833)
- Fix a few broken links in the html docs (#4572)
- Refactor testing utilities (#4829)
- Fix edge drawing performance regression (#4825)
- Draft 2.6 release notes (#4828)
- Fix bad import pattern (#4839)
- Add info about testing and examples (#4582)
- Remove unused `normalized` parameter from communicability_betweenness_centrality (#4843)
- add special processing of `multigraph_input` upon graph init (#4823)
- Add dtype argument to adjacency_matrix (#4850)
- Use scipy to compute eigenvalues (#4847)
- Default to NumPy for simrank_similarity (#4841)
- Remove "networkx" from top-level networkx namespace (#4840)
- Designate 2.6rc1 release
- Bump release version
- DOC: point towards web archive link in GML docs (#4864)
- Fix docstring typo (#4871)
- Reformatted  table to address issue #4852 (#4875)
- spring_layout: ignore 'fixed' nodes not in the graph nodes (#4867)
- Deserializing custom default properties graph ml (#4872)
- DOC: Fix links, use DOI links, wayback machine where required (#4868)
- Fix conda instructions (#4884)
- Decode GraphML/yEd shape type (#4694)
- bugfix-for-issue-4353: modify default edge_id format (#4842)
- Raise ValueError if None is added as a node. (#4892)
- Update arrows default value in draw_networkx. (#4883)
- Doc/fix 403 error drawing custom icons (#4906)
- Remove decorator dependency (#4739)
- Update docstrings for dfs and bfs edges and fix cross links (#4900)
- Fix graph_class usage in to_undirected method (#4912)
- Fix assortativity coefficent calculation (#4851)
- Deprecate numeric_mixing_matrix. (#4923)
- Update read_gml docstring with destringizer ex (#4916)
- Update release process (#4866)
- Designate 2.6rc2 release
- Bump release version
- Add 3.0 migration guide (#4927)
- quotient_graph doc fix (#4930)
- Page number for Katz centrality reference (#4932)
- Expand destringizer example in read_gml docstring (#4925)
- move partition checking outside private _quotient_graph function (#4931)
- Fixes #4275 - Add comment to parallel betweenness example (#4926)
- Minor Improvements on Networkx/algorithms/community/quality.py (#4939)
- Fix numeric and degree assortativity coefficient calculation (#4928)
- fix spelling in docstring of conftest.py (#4945)
- fix trouble with init_cycle argument to two TSP functions (#4938)
- split out deprecation. remove all changes to neighbor_degree (#4937)
- Add matrix market to readwrite reference (#4934)
- fix typo for PR number of deprecation (#4949)
- Fix neighbor degree for directed graphs (#4948)
- `descendants_at_distance` also for non-DiGraphs (#4952)
- Changes to rst files to make doctests pass (#4947)
- Fix version pull down (#4954)
- Finalize 2.6 release notes (#4958)

Contributors
------------

- AbhayGoyal
- Suvayu Ali
- Alexandre Amory
- Francesco Andreuzzi
- Salim BELHADDAD
- Ross Barnowski
- Raffaele Basile
- Jeroen Bergmans
- R. Bernstein
- Geoff Boeing
- Kelly Boothby
- Jeff Bradberry
- Erik Brendel
- Justin Cai
- Thomas A Caswell
- Jonas Charfreitag
- Berlin Cho
- ChristopherReinartz
- Jon Crall
- Michael Dorner
- Harshal Dupare
- Andrew Eckart
- Tomohiro Endo
- Douglas Fenstermacher
- Martin Fleischmann
- Martha Frysztacki [frɨʂtat͡skʲ]
- Debargha Ganguly
- CUI Hao
- Floris Hermsen
- Ward Huang
- Elgun Jabrayilzade
- Han Jaeseung
- Mohammed Kashif
- Alex Korbonits
- Mario Kostelac
- Sebastiaan Lokhorst
- Lonnen
- Delille Louis
- Xiaoyan Lu
- Alex Malins
- Oleh Marshev
- Jordan Matelsky
- Fabio Mazza
- Chris McBride
- Abdulelah S. Al Mesfer
- Attila Mester
- Jarrod Millman
- Miroslav Šedivý
- Harsh Mishra
- S Murthy
- Matthias Nagel
- Attila Nagy
- Mehdi Nemati
- Dimitrios Papageorgiou
- Vitaliy Pozdnyakov
- Bharat Raghunathan
- Randy
- Michael Recachinas
- Carlos González Rotger
- Taxo Rubio
- Dan Schult
- Mridul Seth
- Kunal Shah
- Eric Sims
- Ludovic Stephan
- Justin Timmons
- Andrea Tomassilli
- Matthew Treinish
- Milo Trujillo
- Danylo Ulianych
- Alex Walker
- Stefan van der Walt
- Anthony Wilder Wohns
- Levi John Wolf
- Xiangyu Xu
- Shichu Zhu
- alexpsimone
- as1371
- cpurmessur
- dbxnr
- wim glenn
- goncaloasimoes
- happy
- jason-crowley
- jebogaert
- josch
- ldelille
- marcusjcrook
- guy rozenberg
- tom
- walkeralexander
NetworkX 2.0
============

Release date: 20 September 2017

Support for Python 3.6 added, drop support for Python 3.3.

See :doc:`migration_guide_from_1.x_to_2.0`.

NetworkX is a Python package for the creation, manipulation, and study of the
structure, dynamics, and functions of complex networks.

For more information, please visit our `website <https://networkx.org/>`_
and our `gallery of examples
<https://networkx.org/documentation/latest/auto_examples/index.html>`_.
Please send comments and questions to the `networkx-discuss mailing list
<http://groups.google.com/group/networkx-discuss>`_.

Highlights
----------

This release is the result of over two years of work with 1212 commits and
193 merges by 86 contributors. Highlights include:

- We have made major changes to the methods in the Multi/Di/Graph classes.
  There is a `migration guide for people moving from 1.X to 2.0
  <https://networkx.org/documentation/latest/release/migration_guide_from_1.x_to_2.0.html>`_.

- We updated the documentation system.

API Changes
-----------

* Base Graph Class Changes
  With the release of NetworkX 2.0 we are moving towards a view/iterator reporting API.
  We used to have two methods for the same property of the graph, one that returns a
  list and one that returns an iterator. With 2.0 we have replaced them with a view.
  A view is a read-only object that is quick to create, automatically updated, and
  provides basic access like iteration, membership and set operations where appropriate.
  For example, ``G.nodes()`` used to return a list and ``G.nodes_iter()`` an iterator.
  Now ``G.nodes()`` returns a view and ``G.nodes_iter()`` is removed. ``G.degree()``
  returns a view with ``(node, degree)`` iteration, so that ``dict(G.degree())``
  returns a dict keyed by node with degree as value.
  The old behavior

    >>> G = nx.complete_graph(5)
    >>> G.nodes()  # doctest: +SKIP
    [0, 1, 2, 3, 4]
    >>> G.nodes_iter()  # doctest: +SKIP
    <dictionary-keyiterator at ...>

  has changed to

    >>> G = nx.complete_graph(5)
    >>> G.nodes()
    NodeView((0, 1, 2, 3, 4))
    >>> list(G.nodes())
    [0, 1, 2, 3, 4]

  New feature include lookup of node and edge data from the views, property
  access without parentheses, and set operations.

    >>> G.add_node(3, color='blue')
    >>> G.nodes[3]
    {'color': 'blue'}
    >>> G.nodes & {3, 4, 5}
    {3, 4}

  The following methods have changed:

    * Graph/MultiGraph

      * ``G.nodes()``
      * ``G.edges()``
      * ``G.neighbors()``
      * ``G.adjacency_list()`` and ``G.adjacency_iter()`` to ``G.adjacency()``
      * ``G.degree()``
      * ``G.subgraph()``
      * ``G.copy()``
      * ``G.__class__()`` should be replaced with ``G.fresh_copy()``

    * DiGraph/MultiDiGraph

      * ``G.nodes()``
      * ``G.edges()``
      * ``G.in_edges()``
      * ``G.out_edges()``
      * ``G.degree()``
      * ``G.in_degree()``
      * ``G.out_degree()``
      * ``G.reverse()``

  The following deprecated methods will be removed in a future release (3.0?).
      * ``G.node``, ``G.edge`` (replaced by G.nodes, G.edges)
      * ``G.add_path``, ``G.add_cycle``, ``G.add_star`` (Now ``nx.add_path(G,...``)
      * ``G.selfloop_edges``, ``G.nodes_with_selfloops``, ``G.number_of_selfloops``
        (Now ``nx.selfloop_edges(G)``, etc)

  Many subclasses have been changed accordingly such as:
    * AntiGraph
    * OrderedGraph and friends
    * Examples such as ThinGraph that inherit from Graph

* [`#2107 <https://github.com/networkx/networkx/pull/2107>`_]
  The Graph class methods ``add_edge`` and ``add_edges_from`` no longer
  allow the use of the ``attr_dict`` parameter.  Instead use keyword arguments.
  Thus ``G.add_edge(1, 2, {'color': 'red'})`` becomes
  ``G.add_edge(1, 2, color='red')``.
  Note that this only works if the attribute name is a string. For non-string
  attributes you will need to add the edge and then update manually using
  e.g. ``G.edges[1, 2].update({0: "zero"})``

* [`#1577 <https://github.com/networkx/networkx/pull/1577>`_]
  In addition to minimum spanning trees, a new function for calculating maximum
  spanning trees is now provided. The new API consists of four functions:
  ``minimum_spanning_edges``, ``maximum_spanning_edges``,
  ``minimum_spanning_tree``, and ``maximum_spanning_tree``.
  All of these functions accept an ``algorithm`` parameter which specifies the
  algorithm to use when finding the minimum or maximum spanning tree. Currently,
  Kruskal's and Prim's algorithms are implemented, defined as 'kruskal' and
  'prim', respectively. If nothing is specified, Kruskal's algorithm is used.
  For example, to calculate the maximum spanning tree of a graph using Kruskal's
  algorithm, the function ``maximum_spanning_tree`` has to be called like::

      >>> nx.maximum_spanning_tree(G, algorithm='kruskal')

  The ``algorithm`` parameter is new and appears before the existing ``weight``
  parameter. So existing code that did not explicitly name the optional
  ``weight`` parameter will need to be updated::

      >>> nx.minimum_spanning_tree(G, 'mass')  # old
      >>> nx.minimum_spanning_tree(G, weight='mass') # new

  In the above, we are still relying on the functions being imported into the
  top-level  namespace. We do not have immediate plans to deprecate this approach,
  but we recommend the following instead::

       >>> from networkx.algorithms import tree
       # recommended
       >>> tree.minimum_spanning_tree(G, algorithm='kruskal', weight='mass')
       >>> tree.minimum_spanning_edges(G, algorithm='prim', weight='mass')

* [`#1445 <https://github.com/networkx/networkx/pull/1445>`_]
  Most of the ``shortest_path`` algorithms now raise a ``NodeNotFound`` exception
  when a source or a target are not present in the graph.

* [`#2326 <https://github.com/networkx/networkx/pull/2326>`_]
  Centrality algorithms were harmonized with respect to the default behavior of
  the weight parameter. The default value of the ``weight`` keyword argument has
  been changed from ``weight`` to ``None``.  This affects the
  following centrality functions:

  - :func:`approximate_current_flow_betweenness_centrality()`
  - :func:`current_flow_betweenness_centrality()`
  - :func:`current_flow_betweenness_centrality_subset()`
  - :func:`current_flow_closeness_centrality()`
  - :func:`edge_current_flow_betweenness_centrality()`
  - :func:`edge_current_flow_betweenness_centrality_subset()`
  - :func:`eigenvector_centrality()`
  - :func:`eigenvector_centrality_numpy()`
  - :func:`katz_centrality()`
  - :func:`katz_centrality_numpy()`

* [`#2420 <https://github.com/networkx/networkx/pull/2420>`_]
  New community detection algorithm provided. Fluid Communities is an
  asynchronous algorithm based on the simple idea of fluids interacting in an
  environment, expanding and pushing each other. The algorithm is completely
  described in `"Fluid Communities: A Competitive and Highly Scalable Community
  Detection Algorithm" <https://arxiv.org/pdf/1703.09307.pdf>`_.

* [`#2510 <https://github.com/networkx/networkx/pull/2510>`_ and
  `#2508 <https://github.com/networkx/networkx/pull/2508>`_]
  ``single_source_dijkstra``, ``multi_source_dijkstra`` and functions that use
  these now have new behavior when ``target`` is specified. Instead of
  returning dicts for distances and paths a 2-tuple of ``(distance, path)`` is
  returned.  When ``target`` is not specified the return value is still 2
  dicts.

* [`#2553 <https://github.com/networkx/networkx/pull/2553>`_]
  ``set_node_attributes()`` and ``set_edge_attributes()`` now accept
  dict-of-dict input of shape ``{node/edge: {name: value}}`` in addition to
  previous valid inputs: ``{node/edge: value}`` and ``value``. The order of the
  parameters changed also: The second parameter "values" is the value argument
  and the third parameter "name" is the name of the attribute. "name" has
  default value ``None`` in which case "values" must be the newly allowed form
  containing names. Previously "name" came second without default, and "values"
  came third.

* [`#2604 <https://github.com/networkx/networkx/pull/2604>`_] Move selfloop
  methods out of base classes to networkx functions.
  ``G.number_of_selfloops()``, ``G.selfloop_edges()``,
  ``G.nodes_with_selfloops()`` are now ``nx.number_of_selfloops(G)``,
  ``nx.selfloop_edges(G)``, ``nx.nodes_with_selfloops(G)``.

  ``G.node`` and ``G.edge`` are removed. Their functionality are replaced by
  ``G.nodes`` and ``G.edges``.

* [`#2558 <https://github.com/networkx/networkx/pull/2558>`_]
  Previously, the function ``from_pandas_dataframe`` assumed that the dataframe
  has edge-list like structures, but ``to_pandas_dataframe`` generates an
  adjacency matrix.  We now provide four functions ``from_pandas_edgelist``,
  ``to_pandas_edgelist``, ``from_pandas_adjacency``, and ``to_pandas_adjacency``.

* [`#2620 <https://github.com/networkx/networkx/pull/2620>`_]
  Removed ``draw_nx``, please use ``draw`` or ``draw_networkx``.

* [`#1662 <https://github.com/networkx/networkx/pull/1662>`_]
  Rewrote ``topolgical_sort`` as a generator.  It no longer accepts
  ``reverse`` or ``nbunch`` arguments and is slightly faster.
  Added ``lexicographical_topological_sort``, which accepts a key.

Deprecations
------------

The following deprecated functions will be removed in 2.1.

- The function ``bellman_ford`` has been deprecated in favor of
  ``bellman_ford_predecessor_and_distance``.

- The functions ``to_pandas_dataframe`` and ``from_pandas_dataframe`` have been
  deprecated in favor of ``to_pandas_adjacency``, ``from_pandas_adjacency``,
  ``to_pandas_edgelist``, and ``from_pandas_edgelist``.

Contributors
------------

- Niels van Adrichem
- Kevin Arvai
- Ali Baharev
- Moritz Emanuel Beber
- Livio Bioglio
- Jake Bogerd
- Moreno Bonaventura
- Raphaël Bournhonesque
- Matthew Brett
- James Clough
- Marco Cognetta
- Jamie Cox
- Jon Crall
- Robert Davidson
- Nikhil Desai
- DonQuixoteDeLaMancha
- Dosenpfand
- Allen Downey
- Enrico
- Jens Erat
- Jeffrey Finkelstein
- Minas Gjoka
- Aravind Gollakota
- Thomas Grainger
- Aric Hagberg
- Harry
- Yawara ISHIDA
- Bilal AL JAMMAL
- Ryan James
- Omer Katz
- Janis Klaise
- Valentin Lorentz
- Alessandro Luongo
- Francois Malassenet
- Arya McCarthy
- Michael-E-Rose
- Peleg Michaeli
- Jarrod Millman
- Chris Morin
- Sanggyu Nam
- Nishant Nikhil
- Rhile Nova
- Ramil Nugmanov
- Juan Nunez-Iglesias
- Pim Otte
- Ferran Parés
- Richard Penney
- Phobia
- Tristan Poupard
- Sebastian Pucilowski
- Alexander Rodriguez
- Michael E. Rose
- Alex Ryan
- Zachary Sailer
- René Saitenmacher
- Felipe Schneider
- Dan Schult
- Scinawa
- Michael Seifert
- Mohammad Hossein Sekhavat
- Mridul Seth
- SkyTodInfi
- Stacey Smolash
- Jordi Torrents
- Martin Törnwall
- Jannis Vamvas
- Luca Verginer
- Prayag Verma
- Peter Wills
- Ianto Lin Xi
- Heqing Ya
- aryamccarthy
- chebee7i
- definitelyuncertain
- jfinkels
- juliensiebert
- leotrs
- leycec
- mcognetta
- numpde
- root
- salotz
- scott-vsi
- thegreathippo
- vpodpecan
- yash14123
- Neil Girdhar

Merged PRs
----------

- Gml read fix. (#1962)
- Small changes leftover from #1847 (#1966)
- Fix k_core for directed graphs. Add tests (#1963)
- Communicability fix (#1958)
- Allows weight functions in shortest path functions (#1690)
- minor doc changes on weighted.py (#1969)
- Fix minimum_st_edge_cut documentation. (#1977)
- Fix all_node_cuts corner cases: cycle and complete graphs. (#1976)
- Change add_path/star/cycle from methods to functions (#1970)
- branch 'edge-subgraph' from @jfinkels (#1740)
- Corrected eppstein matching (#1955)
- Nose ignore docstrings (#1980)
- Edited Doc Makefile so clean doesn't delete the examples folder (#1967)
- bug fix in convert_matrix.py (#1983)
- Avoid unnecessary eigenval sort in pagerank_numpy (#1986)
- Fix a typo in install.rst (#1991)
- Adds unorderable nodes test for dag_longest_path. (#1999)
- Improve drawing test scripts (typos, newlines, methods) (#1992)
- Improves test coverage for A* shortest path. (#1988)
- Improves test coverage for avg degree connectivity (#1987)
- Fix Graph() docstring to reflect input flexibility (#2006)
-  Fix sphinx autosummary doc generation errors. (#2026)
- Improve gexf.py (#2010)
- Readme.rst should mention Decorator package is required. (#2009)
- fix_duplicate_kwarg: Fix a duplicate kwarg that was causing to_agraph… (#2005)
- Cleans documentation for graph6 and sparse6 I/O. (#2002)
- Remove http server example (#2001)
- Generalize and improve docstrings of node_link.py (#2000)
- fix issue #1948 and PEP8 formatting (#2031)
- Uses weight function for dijkstra_path_length. (#2033)
- Change default role for sphinx to 'obj' (#2027)
- fixed typo s/abritrary/arbitrary/ (#2035)
- Fix bug in dtype-valued matrices (#2038)
- Adds example for using Graph.nodes() with default (#2040)
- Clarifies some examples for relabel_nodes(). (#2041)
- Cleans code and documentation for graph power. (#2042)
- Cleans the classes.function module. (#2043)
- UnboundLocalError if called with an empty graph (#2047)
- Standardized Bellman-Ford function calls (#1910)
- Nobody is in IRC (#2059)
- Uses add_weighted_edges_from function in MST test. (#2061)
- Adds multi-source Dijkstra's algorithm (#2073)
- Adds Voronoi cells algorithm (#2074)
- Fixes several issues with the Girvan-Newman partitioning function. Fixes #1703, #1725, #1799  (#1972)
- Moves is_path from utils to simple_paths. (#1921)
- add max_iter and tol parameter for numpy version (#2013)
- Remove draw_graphviz function. Fixes #1997 (#2077)
- Fixes #1998 edge_load function needs documentation. (#2075)
- Update fixcoverage.py (#2080)
- Support digraphs in approximate min vertex cover (#2039)
- Simplifies code in functions for greedy coloring. (#1680)
- Allows arbitrary metric in geometric generators. (#1679)
- Fix spring_layout for single node graph. (#2081)
- Updates set_{node,edge}_attributes and docs. (#1935)
- Fixes tests for maximal matching. (#1919)
- Adds LFM benchmark graph generator for communities (#1727)
- Adds global and local efficiency functions. (#1521)
- Apply alphas to individual nodes (#1289)
- Code and tests for temporal VF2 (#1653)
- extend convert_bool in gexf.py and graphml.py to all valid boolean  (#1063)
- Remove encoded ... to plain ascii (#2086)
- Use not_implemented_for() for in_degree_centrality() and out_degree_centrality() (#2084)
- Issue 2072 weighted modularity (#2088)
- Simplifies eigenvector centrality implementation. (#1708)
- Fjmalass nodes as tuples (#2089)
- Generator rename (#2090)
- Ensure links in doc ```See also``` sections (#2082)
- Document integer-only numeric mixing (#2085)
- doc sphinx error removal (#2091)
- Correct see also links (#2095)
- Adjust layout.py function signatures, docs, exposure (#2096)
- Adds missing __all__ attributes. (#2098)
- Fixes 2 bugs in dominance frontier code (#2092)
- Created two new files: joint_degree_seq.py and test_joint_degree_seq.… (#2011)
- Adds Borůvka's minimum spanning tree algorithm. (#1873)
- Adds global/local reaching centrality functions. (#2099)
- Remove conflicts from #1894 (Update Exception Classes) (#2100)
- Add Exceptions for missing source in shortest_path (#2102)
- Docs for compose now warn about MultiGraph edgekeys (#2101)
- Improve Notes section on simplex and friends docs. (#2104)
- Add Dinitz' algorithm for maximum flow problems. (#1978)
- Removed duplicated method/doc (add_edges_from) (#1)
- Bugfix for generic_multiedge_match (Issue #2114) (#2124)
- Fix for 2015. (#2)
- add_node, add_edge attr_dict change. (#2132)
- Handle graph name attribute in relabel_nodes (#2136)
- Fix fruchterman reingold bug and add more tests to layouts. (#2141)
- Adds exception: failed power iteration convergence (#2143)
- Tweak iteration logic of HITS (#2142)
- Fix PageRank personalize docstring (#2148)
- Set default source=None for dfs_tree (#2149)
- Fix docs for maximal_matching and tensor_product (#2158)
- Isolate edge key generation in multigraphs (#2150)
- Sort centralities together and outsource dispersion (#2083)
- Changed classic generators to use generators instead of lists (#2167)
- Adds beam search traversal algorithm with example (#2129)
- Turan graph (#2172)
- Removes irrelevant Notes section from docstring (#2178)
- Corrects logarithm base in example (#2179)
- Minor correction in documentation (#2180)
- Add Boykov Kolmogorov algorithm for maximum flow problems. (#2122)
- Remove temporary files after tests are run. (#2202)
- Add support for subgraphs with no edges in convert_matrix.to_scipy_sparse_matrix. (#2199)
- Add support for reading adjacency matrix in readwrite.pajek.parse_pajek. (#2200)
- Moves Graph Atlas to data file. (#2064)
- Refactor Dinitz' algorithm implementation. (#2196)
- Use arrays instead of matrices in scipy.linalg.expm() (#2208)
- Making in_edges equivalent to out_edges (#2206)
- Fix tests failing because of ordering issues. (#2207)
- Fix code escaping. (#2214)
- Add adjlist_outer_dict_factory. (#2222)
- Typo in scale free network generator documentation (#2225)
- Add link to nx.drawing.layout instead of mentionning nx.layout. (#2224)
- Example not working in tutorial (#2230)
- don't assume nodes are sortable when running dag_longest_path (#2228)
- Correct typo (#2236)
- Use ego graph when computing local efficiency (#2246)
- Make harmonic centrality more memory-efficient (#2247)
- have dag_longest_path_length return path length, not edge count (#2237)
- Added transitive_reduction in dag (#2215)
- alpha kwarg not used in pylab label drawing, added it here.   (#2269)
- Make PyDot Support Great Again (#2272)
- Unnecessary array copying in katz_centrality_numpy ? (#2287)
- Switch to faster smallest-last algorithm implementation. (#2268)
- Adds example for getting all simple edge paths. Fixes #718  (#2260)
- Remove obsolete testing tools. (#2303)
- Correct error in minimum_spanning_arborescence (#2285)
- Yield string, not dict, in dfs_labeled_edges. (#2277)
- Removes unnecessary convert_to_(un)directed func (#2259)
- Complete multipartite graph docs (#2221)
- fix LPA bug, see issues/2219 (#2227)
- Generalized degree (#2220)
- Turan docs (#2218)
- Fix broken link to the description of the P2G format. (#2211)
- Test ordering (#2209)
- add example of node weights (#2250)
- added paramether nbunch (#2253)
- Adds unit tests for using dtype with to_numpy_matrix (#2257)
- Adds chain decomposition algorithm. (#2284)
- add the Hoffman-Singleton graph (#2275)
- Allow grid_graph generator to accept tuple dim argument (#2320)
- psuedo -> pseudo (fixing typo) (#2322)
- Corrects navigable small world graph param docs (#2321)
- Fix bug in find_cycle. (#2324)
- flip source target (#2309)
- Simpler version of digitsrep(..) function (#2330)
- change articulation_points so that it only returns every vertex once (#2333)
- Use faster random geometric graph implementation. (#2337)
- Allow community asyn_lpa test to have two answers (#2339)
- Fix broken links and remove pdf files from Makefile (#2344)
- Documents orderable node requirement for isom. (#2302)
- Adds modularity measure for communities. (#1729)
- Simplifies degree sequence graph generators. (#1866)
- Adds tree encoding and decoding functions. (#1874)
- Corrects number_of_edges docs for directed graphs (#2360)
- Adds multigraph keys to Eulerian circuits (#2359)
- Update predecessors/successors in edge subgraph (#2373)
- Fix for #2364 (#2372)
- Raise an Exception for disconnected Graphs in bipartite.sets (#2375)
- fixes typo in NetworkXNotImplemented (#2385)
- Check alternating paths using iterative DFS in to_vertex_cover. (#2386)
- Fix typos in generating NXError in networkx.linalg.graphmatrix.incidence_matrix (#2395)
- [Fixes #2342] remove calls to plt.hold(), deprecated in mpl2.0 (#2397)
- Fix broken links (#2414)
- Fix all tests for 3.6 (#2413)
- Improve bipartite documentation. (#2402)
- correct logic in GEXFWriter (#2399)
- list optional dependencies in setup.py (#2398)
- Gitwash update (#2371)
- Added cytoscape JSON handling (#2351)
- Fix for issues #2328 and #2332 (#2366)
- Workaround for gdal python3.6 at travis and more doctests fixes (#2416)
- Fixed bug on custom attrs usage: unavailable iteritems method for dict. (#2461)
- Fix sphinx errors and class outlines (#2480)
- Note the precondition that graphs are directed and acyclic (#2500)
- Add CONTRIBUTE file (#2501)
- Remove external module (#2521)
- Ensure `make html` doesn't fail build on exit (#2530)
- Cherry pick missing commits (#2535)
- Document release process (#2539)
- Update copyright (#2551)
- Remove deprecated code (#2536)
- Improve docs (#2555)
- WIP: Add note on how to estimate appropriate values for alpha (#2583)
- Travis refactor (#2596)
- Create separate functions for df as edge-lists and adjacency matrices (#2558)
- Use texext for math_dollar (#2609)
- Add drawing tests (#2617)
- Add threshold tests (#2622)
- Update docs (#2623)
- Prep beta release (#2624)
- Refactor travis tests and deploy docs with travis (#2647)
- matplotlib 2.1 deprecated is_string_like (#2659)
- topolgical_sort, lexicographical_topological_sort (#1662)
NetworkX 2.1
============

Release date: 22 January 2018

Supports Python 2.7, 3.4, 3.5, and 3.6.

NetworkX is a Python package for the creation, manipulation, and study of the
structure, dynamics, and functions of complex networks.

For more information, please visit our `website <https://networkx.org/>`_
and our `gallery of examples
<https://networkx.org/documentation/latest/auto_examples/index.html>`_.
Please send comments and questions to the `networkx-discuss mailing list
<http://groups.google.com/group/networkx-discuss>`_.

Highlights
----------

This release is the result of four months of work with 75 pull requests by
37 contributors. Highlights include:

  - Arrows for drawing DiGraph edges are vastly improved!
    And an example to show them.

  - More than 12 new functions for graph generation, manipulation and/or
    new graph algorithms.

    - Add a large clique size heuristic function (#2830)
    - Add rooted product function (#2825)
    - Label Propagation Community Detection (#2821)
    - Minimum cycle basis (#2823)
    - Add Mycielski Operator (#2785)
    - Adds prefix_tree, dag_to_branching, and example. (#2784)
    - Add inverse_line_graph generator from #2241 (#2782)
    - Steiner tree and metric closure. (#2252)
    - Add flow based node and edge disjoint paths. (#2063)
    - Update geometric networks with new models (#2498)
    - Graph edit distance (#2729)
    - Added function for finding a k-edge-augmentation (#2572)

  - G.name is no longer processed by graph operators. It remains as a
    property mechanism to access ``G.graph['name']`` but the user is in
    charge of updating or changing it for copies, subgraphs, unions and
    other graph operations.

Improvements
------------

  - Many bug fixes, documentation changes.
  - Speed improvements especially for subgraphs.
  - Changed input variable names for functions using ``**kwds``
    to avoid name collisions -- especially ``add_node``
  - New examples for arrows and spectral embedding of the grid graph.

API Changes
-----------

* [`#2498 <https://github.com/networkx/networkx/pull/2498>`_]
  In ``geographical_threshold_graph``, starting in NetworkX 2.1 the parameter
  ``alpha`` is deprecated and replaced with the customizable ``p_dist``
  function parameter, which defaults to r^-2
  if ``p_dist`` is not supplied. To reproduce networks of earlier NetworkX
  versions, a custom function needs to be defined and passed as the ``p_dist``
  parameter. For example, if the parameter ``alpha`` = 2 was used in NetworkX 2.0,
  the custom function def custom_dist(r): r**-2 can be passed in versions >=2.1
  as the parameter p_dist = custom_dist to produce an equivalent network.
  Note the change in sign from +2 to -2 in this parameter change.

* [`#2554 <https://github.com/networkx/networkx/issues/2554>`_]
  New algorithms for finding k-edge-connected components and k-edge-connected
  subgraphs in directed and undirected graphs. Efficient implementations are
  provided for the special case of k=1 and k=2. The new functionality is
  provided by:

     :func:`k_edge_components()`

     :func:`k_edge_subgraphs()`

* [`#2572 <https://github.com/networkx/networkx/issues/2572>`_]
  New algorithm finding for finding k-edge-augmentations in undirected graphs.
  Efficient implementations are provided for the special case of k=1 and k=2.
  New functionality is provided by:

   - :func:`k_edge_augmentation()`

* [`#2812 <https://github.com/networkx/networkx/pull/2812>`_]
  Removed ``bellman_ford``, please use
  ``bellman_ford_predecessor_and_distance``.

* [`#2811 <https://github.com/networkx/networkx/pull/2811>`_]
  Removed ``to_pandas_dataframe`` and ``from_pandas_dataframe``, please use
  ``to_pandas_adjacency``, ``from_pandas_adjacency``, ``to_pandas_edgelist``,
  or ``from_pandas_edgelist``.

* [`#2766 <https://github.com/networkx/networkx/pull/2766>`_]
  Add seed keyword argument to random_layout and spring_layout

* [`#2776 <https://github.com/networkx/networkx/pull/2776>`_]
  Add threshold option to spring layout

* [`#2774 <https://github.com/networkx/networkx/pull/2774>`_]
  max_weight_matching returns set of edges

* [`#2753 <https://github.com/networkx/networkx/pull/2753>`_]
  Add directed graphs support for jit_graph reading

* [`#2788 <https://github.com/networkx/networkx/pull/2788>`_]
  Control node-border color in draw_networkx_nodes

Deprecations
------------

* [`#2819 <https://github.com/networkx/networkx/pull/2819>`_]
  Deprecate ``connected_component_subgraphs``, ``biconnected_component_subgraphs``,
  ``attracting_component_subgraphs``, ``strongly_connected_component_subgraphs``,
  ``weakly_connected_component_subgraphs``.
  Instead use: ``[G.subgraph(c) for c in *_components]``

Contributors
------------

- Jack Amadeo
- Boskovits
- Daniel Bradburn
- David Bradway
- Ariel Chinn
- Jon Crall
- Rodrigo Dorantes-Gilardi
- Bradley Ellert
- Adam Erispaha
- Ioannis Filippidis
- ForFer
- Louis Gatin
- Aric Hagberg
- Harry
- Huston Hedinger
- Charles Tapley Hoyt
- James Lamb
- Sanghack Lee
- MD
- Cole MacLean
- Marco
- Jarrod Millman
- Sanggyu Nam
- Viraj Parimi
- Dima Pasechnik
- Richard Penney
- Naresh Peshwe
- Zachary Sailer
- Dan Schult
- Jordi Torrents
- John Wegis
- aparamon
- aweltsch
- gfyoung
- md0000
- mddddd
- talhum


Merged PRs
----------

- Update Release Notes for v2.1 (#2839)
- Update release notes (#2838)
- Update copyright (#2837)
- Add a large clique size heuristic function (#2830)
- Remove automatic processing of G.name attribute (#2829)
- Add rooted product function (#2825)
- Label Propagation Community Detection (#2821)
- change variable names to avoid kwargs clobber (#2824)
- Minimum cycle basis (#2823)
- Deprecate component_subgraphs functions (#2819)
- Temporarily disable sphinx doctests (#2818)
- Adjust docs for graph class edge attrib assignment (#2817)
- Add directed graphs support for jit_graph reading (#2753)
- Arrows as a plot example. (#2801)
- Fix bug in len(edges) for self-loops (#2816)
- MRG: Remove ``to_pandas_dataframe`` and ``from_pandas_dataframe`` (#2811)
- Fix Pydot tests so works with new version 1.2.4 (#2815)
- MRG: Remove ``bellman_ford`` (#2812)
- Combine generator modules and tweak docs (#2814)
- Legacy array printing for NumPy 1.14+ (#2810)
- Fix rare structurally forbidden mappings bug. (#2798)
- Digraph Arrows to fix #2757 (#2760)
- use a generic Integral type for parameters check (#2800)
- Control node-border color in draw_networkx_nodes (#2788)
- Add seed keyword argument to random_layout and spring_layout (#2766)
- Add Mycielski Operator (#2785)
- Adds prefix_tree, dag_to_branching, and example. (#2784)
- Add inverse_line_graph generator from #2241 (#2782)
- Add docs for steiner_tree and metric_closure (#2783)
- Steiner tree and metric closure. (#2252)
- Correct docstring for weight parameter (#2781)
- Switch to xcode 7.3 for osx_image in .travis.yml (#2780)
- Change how sparse6 tests filenames (#2779)
- Add flow based node and edge disjoint paths. (#2063)
- Update geometric networks with new models (#2498)
- [WIP] Graph edit distance 2361 (#2729)
- max_weight_matching returns set of edges (#2774)
- Avoid keyword and attribute clash (#2775)
- Add threshold option to spring layout (#2776)
- Fix bug in expected_degree_graph generator (#2773)
- Add support for incomplete partitions in quotient_graph. (#2771)
- Fix SOURCE_DATE_EPOCH ignored bug (#2735) (#2736)
- Makes write_graph6 less memory-intensive. (#2299)
- all_simple_paths should not return cycles. Fix issue #2762 (#2770)
- Fix typo in write_gml and add test (#2769)
- Fix bug and add checks for non-convergent fiedler_vector (#2681)
- Dictionary comprehensions from #1700 merged conflicts (#2768)
- Fix 2763: Typo `furether` in networkx tutorial documentation (#2764)
- Fix #2726: ensure add_path to add the first node (#2759)
- a minor correction in docs (#2751)
- Speedups for subgraph and copy methods (#2744)
- fix typo in tutorial (#2746)
- Expand documentation regarding strong connectivity (#2732)
- Correct when we raise NetworkXNotImplemented (#2731)
- removed list conversion from _triangles_and_degree_iter (#2725)
- nx_shp fixes (#2721)
- removed reference to create_using from union docs (#2722)
- Copy graph in transitive closure algorithm. (#2718)
- Fix dag_longest_path bug (#2703)
- Fix for inter_community_edges (#2713)
- Fix shortest_simple_paths. Issue #2427 (#2712)
- Update migration_guide_from_1.x_to_2.0.rst (#2694)
- mention `doc.txt` in `requirements/README.md` (#2699)
- docs(centrality/dispersion): updating contributor email address (#2698)
- Fixes bug #2503 by removing arrow labels (#2696)
- Add example of spectral embedding of the grid graph (#2690)
- Fix create_using of nx.from_pandas_adjacency() (#2693)
- Added function for finding a k-edge-augmentation (#2572)
- rm arg `strict` from function `networkx.drawing.nx_pydot.to_pydot` (#2672)
- Fixed problem parsing graphml with nodes in groups (#2644)
- Remove unused imports (#2653)
- Improve subgraph node iteration (#2687)
- Added Kamada-Kawai functions to Sphinx documentation (#2680)
- unpacked dict to provide kwargs when creating nodes from shapefiles (#2678)
- Fix typo in documentation (#2677)
NetworkX 1.11
=============

Release date: 30 January 2016

Support for Python 3.5 added, drop support for Python 3.2.

Highlights
~~~~~~~~~~

Pydot features now use pydotplus.
Fixes installation on some machines and test with appveyor.
Restores default center and scale of layout routines.
Fixes various docs including no symbolic links in examples.
Docs can now build using autosummary on readthedocs.org.

API changes
-----------
* [`#1930 <https://github.com/networkx/networkx/pull/1930>`_]
  No longer import nx_agraph and nx_pydot into the top-level namespace.
  They can be accessed within networkx as e.g. ``nx.nx_agraph.write_dot``
  or imported as ``from networkx.drawing.nx_agraph import write_dot``.

* [`#1750 <https://github.com/networkx/networkx/pull/1750>`_]
  Arguments center and scale are now available for all layout functions.
  The defaul values revert to the v1.9 values (center is the origin
  for circular layouts and domain is [0, scale) for others.

* [`#1924 <https://github.com/networkx/networkx/pull/1924>`_]
  Replace pydot with pydotplus for drawing with the pydot interface.

* [`#1888 <https://github.com/networkx/networkx/pull/1888>`_]
  Replace support for Python3.2 with support for Python 3.5.

Miscellaneous changes
---------------------

* [`#1763 <https://github.com/networkx/networkx/pull/1763>`_]
  Set up appveyor to automatically test installation on Windows machines.
  Remove symbolic links in examples to help such istallation.

Change many doc_string typos to allow sphinx
to build the docs without errors or warnings.

Enable the docs to be automatically built on
readthedocs.org by changing requirements.txt
NetworkX 1.9
============

Release date: 21 June 2014

Support for Python 3.1 is dropped in this release.

Highlights
----------
- Completely rewritten maximum flow and flow-based connectivity algorithms with
  backwards incompatible interfaces
- Community graph generators
- Stoer–Wagner minimum cut algorithm
- Linear-time Eulerian circuit algorithm
- Linear algebra package changed to use SciPy sparse matrices
- Algebraic connectivity, Fiedler vector, spectral ordering algorithms
- Link prediction algorithms
- Goldberg–Radzik shortest path algorithm
- Semiconnected graph and tree recognition algorithms

Flow package
------------

The flow package (:samp:`networkx.algorithms.flow`) is completely rewritten
with backward *incompatible* changes. It introduces a new interface to flow
algorithms. Existing code that uses the flow package will not work unmodified
with NetworkX 1.9.

Main changes
============

1. We added two new maximum flow algorithms (:samp:`preflow_push` and
   :samp:`shortest_augmenting_path`) and rewrote the Edmonds–Karp algorithm in
   :samp:`flow_fulkerson` which is now in :samp:`edmonds_karp`.
   `@ysitu <https://github.com/ysitu>`_ contributed implementations of all new
   maximum flow algorithms. The legacy Edmonds–Karp algorithm implementation in
   :samp:`ford_fulkerson` is still available but will be removed in the next
   release.

2. All maximum flow algorithm implementations (including the legacy
   :samp:`ford_fulkerson`) output now a residual network (i.e., a
   :samp:`DiGraph`) after computing the maximum flow. See :samp:`maximum_flow`
   documentation for the details on the conventions that NetworkX uses for
   defining a residual network.

3. We removed the old :samp:`max_flow` and :samp:`min_cut` functions. The main
   entry points to flow algorithms are now the functions :samp:`maximum_flow`,
   :samp:`maximum_flow_value`, :samp:`minimum_cut` and
   :samp:`minimum_cut_value`, which have new parameters that control maximum
   flow computation: :samp:`flow_func` for specifying the algorithm that will
   do the actual computation (it accepts a function as argument that implements
   a maximum flow algorithm), :samp:`cutoff` for suggesting a maximum flow
   value at which the algorithm stops, :samp:`value_only` for stopping the
   computation as soon as we have the value of the flow, and :samp:`residual`
   that accepts as argument a residual network to be reused in repeated maximum
   flow computation.

4. All flow algorithms are required to accept arguments for these parameters
   but may selectively ignored the inapplicable ones. For instance,
   :samp:`preflow_push` algorithm can stop after the preflow phase without
   computing a maximum flow if we only need the flow value, but both
   :samp:`edmonds_karp` and :samp:`shortest_augmenting_path` always compute a
   maximum flow to obtain the flow value.

5. The new function :samp:`minimum_cut` returns the cut value and a node
   partition that defines the minimum cut. The function
   :samp:`minimum_cut_value` returns only the value of the cut, which is what
   the removed :samp:`min_cut` function used to return before 1.9.

6. The functions that implement flow algorithms (i.e., :samp:`preflow_push`,
   :samp:`edmonds_karp`, :samp:`shortest_augmenting_path` and
   :samp:`ford_fulkerson`) are not imported to the base NetworkX namespace. You
   have to explicitly import them from the flow package:

>>> from networkx.algorithms.flow import (ford_fulkerson, preflow_push,
...        edmonds_karp, shortest_augmenting_path)  # doctest: +SKIP


7. We also added a capacity-scaling minimum cost flow algorithm:
   :samp:`capacity_scaling`. It supports :samp:`MultiDiGraph` and disconnected
   networks.

Examples
========

Below are some small examples illustrating how to obtain the same output than in
NetworkX 1.8.1 using the new interface to flow algorithms introduced in 1.9:

>>> import networkx as nx
>>> G = nx.icosahedral_graph()
>>> nx.set_edge_attributes(G, 'capacity', 1)

With NetworkX 1.8:

>>> flow_value = nx.max_flow(G, 0, 6)  # doctest: +SKIP
>>> cut_value = nx.min_cut(G, 0, 6)  # doctest: +SKIP
>>> flow_value == cut_value  # doctest: +SKIP
True
>>> flow_value, flow_dict = nx.ford_fulkerson(G, 0, 6)  # doctest: +SKIP

With NetworkX 1.9:

>>> from networkx.algorithms.flow import (ford_fulkerson, preflow_push,
...        edmonds_karp, shortest_augmenting_path)  # doctest: +SKIP
>>> flow_value = nx.maximum_flow_value(G, 0, 6)  # doctest: +SKIP
>>> cut_value = nx.minimum_cut_value(G, 0, 6)  # doctest: +SKIP
>>> flow_value == cut_value  # doctest: +SKIP
True
>>> # Legacy: this returns the exact same output than ford_fulkerson in 1.8.1
>>> flow_value, flow_dict = nx.maximum_flow(G, 0, 6, flow_func=ford_fulkerson)  # doctest: +SKIP
>>> # We strongly recommend to use the new algorithms:
>>> flow_value, flow_dict = nx.maximum_flow(G, 0, 6)  # doctest: +SKIP
>>> # If no flow_func is passed as argument, the default flow_func
>>> # (preflow-push) is used. Therefore this is the same than:
>>> flow_value, flow_dict = nx.maximum_flow(G, 0, 6, flow_func=preflow_push)  # doctest: +SKIP
>>> # You can also use alternative maximum flow algorithms:
>>> flow_value, flow_dict = nx.maximum_flow(G, 0, 6, flow_func=shortest_augmenting_path)  # doctest: +SKIP
>>> flow_value, flow_dict = nx.maximum_flow(G, 0, 6, flow_func=edmonds_karp)  # doctest: +SKIP

Connectivity package
--------------------

The flow-based connecitivity and cut algorithms from the connectivity
package (:samp:`networkx.algorithms.connectivity`) are adapted to take
advantage of the new interface to flow algorithms. As a result, flow-based
connectivity algorithms are up to 10x faster than in NetworkX 1.8 for some
problems, such as sparse networks with highly skewed degree distributions.
A few backwards *incompatible* changes were introduced.

* The functions for local connectivity and cuts accept now
  arguments for the new parameters defined for the flow interface:
  :samp:`flow_func` for defining the algorithm that will perform the
  underlying maximum flow computations, :samp:`residual` that accepts
  as argument a residual network to be reused in repeated maximum
  flow computations, and :samp:`cutoff` for defining a maximum flow
  value at which the underlying maximum flow algorithm stops. The big
  speed improvement with respect to 1.8 comes mainly from the reuse
  of the residual network and the use of :samp:`cutoff`.

* We removed the flow-based local connectivity and cut functions from
  the base namespace. Now they have to be explicitly imported from the
  connectivity package. The main entry point to flow-based connectivity
  and cut functions are the functions :samp:`edge_connectivity`,
  :samp:`node_connectivity`, :samp:`minimum_edge_cut`, and
  :samp:`minimum_node_cut`. All these functions accept a couple of nodes
  as optional arguments for computing local connectivity and cuts.

* We improved the auxiliary network for connectivity functions: The node
  mapping dict needed for node connectivity and minimum node cuts is now a
  graph attribute of the auxiliary network. Thus we removed the
  :samp:`mapping` parameter from the local versions of connectivity and cut
  functions. We also changed the parameter name for the auxuliary digraph
  from :samp:`aux_digraph` to :samp:`auxiliary`.

* We changed the name of the function :samp:`all_pairs_node_connectiviy_matrix`
  to :samp:`all_pairs_node_connectivity`. This function now returns a dictionary
  instead of a NumPy 2D array. We added a new parameter :samp:`nbunch` for
  computing node connectivity only among pairs of nodes in :samp:`nbunch`.

* A :samp:`stoer_wagner` function is added to the connectivity package
  for computing the weighted minimum cuts of undirected graphs using
  the Stoer–Wagner algorithm. This algorithm is not based on maximum flows.
  Several heap implementations are also added in the utility package
  (:samp:`networkx.utils`) for use in this function.
  :class:`BinaryHeap` is recommended over :class:`PairingHeap` for Python
  implementations without optimized attribute accesses (e.g., CPython)
  despite a slower asymptotic running time. For Python implementations
  with optimized attribute accesses (e.g., PyPy), :class:`PairingHeap`
  provides better performance.

Other new functionalities
-------------------------

* A :samp:`disperson` function is added in the centrality package
  (:samp:`networkx.algorithms.centrality`) for computing the dispersion of
  graphs.

* A community package (:samp:`networkx.generators.community`) is added for
  generating community graphs.

* An :samp:`is_semiconnected` function is added in the connectivity package
  (:samp:`networkx.algorithms.connectivity`) for recognizing semiconnected
  graphs.

* The :samp:`eulerian_circuit` function in the Euler package
  (:samp:`networkx.algorithm.euler`) is changed to use a linear-time algorithm.

* A :samp:`non_edges` function in added in the function package
  (:samp:`networkx.functions`) for enumerating nonexistent edges between
  existing nodes of graphs.

* The linear algebra package (:samp:`networkx.linalg`) is changed to use SciPy
  sparse matrices.

* Functions :samp:`algebraic_connectivity`, :samp:`fiedler_vector` and
  :samp:`spectral_ordering` are added in the linear algebra package
  (:samp:`networkx.linalg`) for computing the algebraic connectivity, Fiedler
  vectors and spectral orderings of undirected graphs.

* A link prediction package (:samp:`networkx.algorithms.link_prediction`) is
  added to provide link prediction-related functionalities.

* Write Support for the graph6 and sparse6 formats is added in the read/write
  package (:samp:`networkx.readwrite`).

* A :samp:`goldberg_radzik` function is added in the shortest path package
  (:samp:`networkx.algorithms.shortest_paths`) for computing shortest paths
  using the Goldberg–Radzik algorithm.

* A tree package (:samp:`networkx.tree`) is added to provide tree recognition
  functionalities.

* A context manager :samp:`reversed` is added in the utility package
  (:samp:`networkx.utils`) for temporary in-place reversal of graphs.

Miscellaneous changes
---------------------

* The functions in the components package
  (:samp:`networkx.algorithms.components`) such as :samp:`connected_components`,
  :samp:`connected_components_subgraph` now return generators instead of lists.
  To recover the earlier behavior, use :samp:`list(connected_components(G))`.

* JSON helpers in the JSON graph package (:samp:`networkx.readwrite.json_graph`)
  are removed. Use functions from the standard library (e.g.,
  :samp:`json.dumps`) instead.

* Support for Python 3.1 is dropped. Basic support is added for Jython 2.7 and
  IronPython 2.7, although they remain not officially supported.

* Numerous reported issues are fixed.
NetworkX 0.99
=============

Release date:  18 November 2008

See: https://networkx.lanl.gov/trac/timeline

New features
------------
This release has significant changes to parts of the graph API.
See http://networkx.lanl.gov/reference/api_changes.html

 - Update Graph and DiGraph classes to use weighted graphs as default
   Change in API for performance and code simplicity.
 - New MultiGraph and MultiDiGraph classes (replace XGraph and XDiGraph)
 - Update to use Sphinx documentation system http://networkx.lanl.gov/
 - Developer site at https://networkx.lanl.gov/trac/
 - Experimental LabeledGraph and LabeledDiGraph
 - Moved package and file layout to subdirectories.

Bug fixes
---------
 - handle root= option to draw_graphviz correctly

Examples
--------
 - Update to work with networkx-0.99 API
 - Drawing examples now use matplotlib.pyplot interface
 - Improved drawings in many examples
 - New examples - see http://networkx.lanl.gov/examples/

The version networkx-0.99 is the penultimate release before
networkx-1.0.  We have bumped the version from 0.37 to 0.99 to
indicate (in our unusual version number scheme) that this is a major
change to NetworkX.

We have made some significant changes, detailed below, to NetworkX
to improve  performance, functionality, and clarity.

Version 0.99 requires Python 2.4 or greater.

Please send comments and questions to the networkx-discuss mailing list.
http://groups.google.com/group/networkx-discuss

Changes in base classes
-----------------------

The most significant changes are in the graph classes.
We have redesigned the Graph() and DiGraph() classes
to optionally allow edge data.
This change allows Graph and DiGraph to naturally represent
weighted graphs and to hold arbitrary information on edges.

 - Both Graph and DiGraph take an optional argument weighted=True|False.
   When weighted=True the graph is assumed to have numeric edge data
   (with default 1).  The Graph and DiGraph classes in earlier versions
   used the Python None as data (which is still allowed as edge data).

 - The Graph and DiGraph classes now allow self loops.

 - The XGraph and XDiGraph classes are removed and replaced with
   MultiGraph and MultiDiGraph. MultiGraph and MultiDiGraph
   optionally allow parallel (multiple) edges between two nodes.

The mapping from old to new classes is as follows::

 - Graph -> Graph (self loops allowed now, default edge data is 1)
 - DiGraph -> DiGraph (self loops allowed now, default edge data is 1)
 - XGraph(multiedges=False) -> Graph
 - XGraph(multiedges=True) -> MultiGraph
 - XDiGraph(multiedges=False) -> DiGraph
 - XDiGraph(multiedges=True) -> MultiDiGraph


Methods changed
---------------

edges()
^^^^^^^
   New keyword data=True|False keyword determines whether to return
   two-tuples (u,v) (False) or three-tuples (u,v,d) (True)


delete_node()
^^^^^^^^^^^^^
   The preferred name is now remove_node().


delete_nodes_from()
^^^^^^^^^^^^^^^^^^^
   No longer raises an exception on an attempt to delete a node not in
   the graph.  The preferred name is now remove_nodes_from().


delete_edge()
^^^^^^^^^^^^^^
   Now raises an exception on an attempt to delete an edge not in the graph.
   The preferred name is now remove_edge().


delete_edges_from()
^^^^^^^^^^^^^^^^^^^
   The preferred name is now remove_edges_from().


add_edge()
^^^^^^^^^^
   The add_edge() method no longer accepts an edge tuple (u,v)
   directly.  The tuple must be unpacked into individual nodes.

   >>> import networkx as nx
   >>> u='a'
   >>> v='b'
   >>> e=(u,v)
   >>> G=nx.Graph()

   Old

   >>> # G.add_edge((u,v))  # or G.add_edge(e)

   New

   >>> G.add_edge(*e) # or G.add_edge(*(u,v))

   The * operator unpacks the edge tuple in the argument list.

   Add edge now has
   a data keyword parameter for setting the default (data=1) edge
   data.

   >>> # G.add_edge('a','b','foo')  # add edge with string "foo" as data
   >>> # G.add_edge(1,2,5.0)  # add edge with float 5 as data



add_edges_from()
^^^^^^^^^^^^^^^^
   Now can take list or iterator of either 2-tuples (u,v),
   3-tuples (u,v,data) or a mix of both.

   Now has data keyword parameter (default 1) for setting the edge data
   for any edge in the edge list that is a 2-tuple.


has_edge()
^^^^^^^^^^
   The has_edge() method no longer accepts an edge tuple (u,v)
   directly.  The tuple must be unpacked into individual nodes.

   Old:

   >>> # G.has_edge((u,v))  # or has_edge(e)

   New:

   >>> G.has_edge(*e) # or has_edge(*(u,v))
   True

   The * operator unpacks the edge tuple in the argument list.

get_edge()
^^^^^^^^^^
   Now has the keyword argument "default" to specify
   what value to return if no edge is found.  If not specified
   an exception is raised if no edge is found.

   The fastest way to get edge data for edge (u,v) is to use G[u][v]
   instead of G.get_edge(u,v)


degree_iter()
^^^^^^^^^^^^^
   The degree_iter method now returns an iterator over pairs of (node,
   degree).  This was the previous behavior of degree_iter(with_labels=true)
   Also there is a new keyword weighted=False|True for weighted degree.

subgraph()
^^^^^^^^^^
   The argument inplace=False|True has been replaced with copy=True|False.

   Subgraph no longer takes create_using keyword.  To change the graph
   type either make a copy of
   the graph first and then change type or change type and make
   a subgraph.  E.g.

   >>> G=nx.path_graph(5)
   >>> H=nx.DiGraph(G.subgraph([0,1])) # digraph of copy of induced subgraph

__getitem__()
^^^^^^^^^^^^^
   Getting node neighbors from the graph with G[v] now returns
   a dictionary.

   >>> G=nx.path_graph(5)
   >>> # G[0]
   #  {1: 1}

   To get a list of neighbors you can either use the keys of that
   dictionary or use

   >>> G.neighbors(0)  # doctest: +SKIP
   [1]

   This change allows algorithms to use the underlying dict-of-dict
   representation through G[v] for substantial performance gains.
   Warning: The returned dictionary should not be modified as it may
   corrupt the graph data structure.  Make a copy G[v].copy() if you
   wish to modify the dict.


Methods removed
---------------

info()
^^^^^^
   now a function

   >>> G=nx.Graph(name='test me')
   >>> nx.info(G)  # doctest: +SKIP
   Name:                  test me
   Type:                  Graph
   Number of nodes:       0
   Number of edges:       0


node_boundary()
^^^^^^^^^^^^^^^
   now a function

edge_boundary()
^^^^^^^^^^^^^^^
   now a function

is_directed()
^^^^^^^^^^^^^
   use the directed attribute

   >>> G=nx.DiGraph()
   >>> # G.directed
   #  True

G.out_edges()
^^^^^^^^^^^^^
   use G.edges()

G.in_edges()
^^^^^^^^^^^^
   use

   >>> G = nx.DiGraph()
   >>> R = G.reverse()
   >>> R.edges()  # doctest: +SKIP
   []

   or

   >>> [(v,u) for (u,v) in G.edges()]
   []

Methods added
-------------

adjacency_list()
^^^^^^^^^^^^^^^^
Returns a list-of-lists adjacency list representation of the graph.

adjacency_iter()
^^^^^^^^^^^^^^^^
Returns an iterator of (node, adjacency_dict[node]) over all
nodes in the graph.  Intended for fast access to the internal
data structure for use in internal algorithms.


Other possible incompatibilities with existing code
===================================================

Imports
-------
Some of the code modules were moved into subdirectories.

Import statements such as::

  import networkx.centrality
  from networkx.centrality import *

may no longer work (including that example).

Use either

>>> import networkx # e.g. centrality functions available as networkx.fcn()

or

>>> from networkx import * # e.g. centrality functions available as fcn()

Self-loops
----------
For Graph and DiGraph self loops are now allowed.
This might affect code or algorithms that add self loops
which were intended to be ignored.

Use the methods

   - nodes_with_selfloops()
   - selfloop_edges()
   - number_of_selfloops()

to discover any self loops.

Copy
----
Copies of NetworkX graphs including using the copy() method
now return complete copies of the graph.  This means that all
connection information is copied--subsequent changes in the
copy do not change the old graph.  But node keys and edge
data in the original and copy graphs are pointers to the same data.

prepare_nbunch
--------------
Used internally - now called nbunch_iter and returns an iterator.


Converting your old code to Version 0.99
----------------------------------------

Mostly you can just run the code and python will raise an exception
for features that changed.  Common places for changes are

    - Converting XGraph() to either Graph or MultiGraph
    - Converting XGraph.edges()  to  Graph.edges(data=True)
    - Switching some rarely used methods to attributes (e.g. directed)
      or to functions (e.g. node_boundary)
    - If you relied on the old default edge data being None, you will
      have to account for it now being 1.

You may also want to look through your code for places which could
improve speed or readability.  The iterators are helpful with large
graphs and getting edge data via G[u][v] is quite fast.   You may also
want to change G.neighbors(n) to G[n] which returns the dict keyed by
neighbor nodes to the edge data.  It is faster for many purposes but
does not work well when you are changing the graph.

Releases
********

We don't use semantic versioning.  The first number indicates that we have
made a major API break (e.g., 1.x to 2.x), which has happened once and probably
won't happen again for some time.  The point releases are new versions and may
contain minor API breakage.  Usually, this happens after a one cycle deprecation
period.

.. warning::
   Since we don't normally make bug-fix only releases, it may not make sense
   for you to use ``~=`` as a pip version specifier.

.. toctree::
   :maxdepth: 2

   release_dev
   release_2.6
   release_2.5
   release_2.4
   release_2.3
   release_2.2
   release_2.1
   release_2.0
   api_1.11
   api_1.10
   api_1.9
   api_1.8
   api_1.7
   api_1.6
   api_1.5
   api_1.4
   api_1.0
   api_0.99
   old_release_log
NetworkX 1.4
============

Release date:  23 January 2011

New features
------------
 - :mod:`k-shell,k-crust,k-corona <networkx.algorithms.core>`
 - :mod:`read GraphML files from yEd <networkx.readwrite.graphml>`
 - :mod:`read/write GEXF format files <networkx.readwrite.gexf>`
 - :mod:`find cycles in a directed graph <networkx.algorithms.cycles>`
 - :mod:`DFS <networkx.algorithms.traversal.depth_first_search>` and :mod:`BFS <networkx.algorithms.traversal.breadth_first_search>` algorithms
 - :mod:`chordal graph functions <networkx.algorithms.chordal.chordal_alg>`
 - :mod:`Prim's algorithm for minimum spanning tree <networkx.algorithms.mst>`
 - :mod:`r-ary tree generator <networkx.generators.classic>`
 - :mod:`rich club coefficient <networkx.algorithms.richclub>`
 - NumPy matrix version of :mod:`Floyd's algorithm for all-pairs shortest path  <networkx.algorithms.shortest_paths.dense>`
 - :mod:`read GIS shapefiles <networkx.readwrite.nx_shp>`
 - :mod:`functions to get and set node and edge attributes <networkx.classes.function>`
 - and more, see  https://networkx.lanl.gov/trac/query?status=closed&group=milestone&milestone=networkx-1.4

API changes
-----------
 - :mod:`gnp_random_graph() <networkx.generators.random_graphs>` now takes a
   directed=True|False keyword instead of create_using
 - :mod:`gnm_random_graph() <networkx.generators.random_graphs>` now takes a
   directed=True|False keyword instead of create_using


Algorithms changed
==================

Shortest path
-------------

astar_path(), astar_path_length(), shortest_path(), shortest_path_length(),
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
bidirectional_shortest_path(), dijkstra_path(), dijkstra_path_length(),
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
bidirectional_dijkstra()
^^^^^^^^^^^^^^^^^^^^^^^^
   These algorithms now raise an exception when a source and a target are
   specified and no path exist between these two nodes. The exception is
   a NetworkXNoPath exception.

Next Release
============

Release date: TBD

Supports Python 3.8, 3.9, and ...

NetworkX is a Python package for the creation, manipulation, and study of the
structure, dynamics, and functions of complex networks.

For more information, please visit our `website <https://networkx.org/>`_
and our :ref:`gallery of examples <examples_gallery>`.
Please send comments and questions to the `networkx-discuss mailing list
<http://groups.google.com/group/networkx-discuss>`_.

Highlights
----------

This release is the result of X of work with over X pull requests by
X contributors. Highlights include:

.. warning::
   Hash values observed in outputs of 
   `~networkx.algorithms.graph_hashing.weisfeiler_lehman_graph_hash`
   have changed in version 2.7 due to bug fixes. See gh-4946_ for details.
   This means that comparing hashes of the same graph computed with different
   versions of NetworkX (i.e. before and after version 2.7)
   could wrongly fail an isomorphism test (isomorphic graphs always have matching
   Weisfeiler-Lehman hashes). Users are advised to recalculate any stored graph
   hashes they may have on upgrading.

.. _gh-4946: https://github.com/networkx/networkx/pull/4946#issuecomment-914623654

- Dropped support for Python 3.7.

Improvements
------------


API Changes
-----------

- The values in the dictionary returned by
  `~networkx.drawing.layout.rescale_layout_dict` are now `numpy.ndarray` objects
  instead of tuples. This makes the return type of ``rescale_layout_dict``
  consistent with that of all of the other layout functions.
- A ``FutureWarning`` has been added to ``google_matrix`` to indicated that the
  return type will change from a ``numpy.matrix`` object to a ``numpy.ndarray``
  in NetworkX 3.0.

Deprecations
------------

- [`#5055 <https://github.com/networkx/networkx/pull/5055>`_]
  Deprecate the ``random_state`` alias in favor of ``np_random_state``
- [`#5114 <https://github.com/networkx/networkx/pull/5114>`_]
  Deprecate the ``name`` kwarg from ``union`` as it isn't used.
- [`#5143 <https://github.com/networkx/networkx/pull/5143>`_]
  Deprecate ``euclidean`` in favor of ``math.dist``.
- [`#5166 <https://github.com/networkx/networkx/pull/5166>`_]
  Deprecate the ``hmn`` and ``lgc`` modules in ``node_classification``.
- [`#5262 <https://github.com/networkx/networkx/pull/5262>`_]
  Deprecate ``to_scipy_sparse_matrix`` and ``from_scipy_sparse_matrix`` in
  favor of ``to_scipy_sparse_array`` and ``from_scipy_sparse_array``, respectively.


Merged PRs
----------

<output of contribs.py>


Contributors
------------

<output of contribs.py>
NetworkX 2.3
============

Release date: 11 April 2019

Supports Python 3.5, 3.6 and 3.7.
This is our first Python 3 only release.

We're happy to announce the release of NetworkX 2.3!
NetworkX is a Python package for the creation, manipulation, and study of the
structure, dynamics, and functions of complex networks.

For more information, please visit our `website <https://networkx.org/>`_
and our `gallery of examples
<https://networkx.org/documentation/latest/auto_examples/index.html>`_.
Please send comments and questions to the `networkx-discuss mailing list
<http://groups.google.com/group/networkx-discuss>`_.

Highlights
----------

This release is the result of 6 months of work with over 92 pull requests by
30 contributors. Highlights include:

- Dropped support for Python 2. We are no longer supporting Python 2.7 and we will
  start changing code to take advantage of Python 3 features we couldn't before.
- Added some Moral Graph analysis functions.
- Enable matplotlib drawing using curved arrows via connectionstyle parameter.
- Remove ticks and axes labels from matplotlib plots.
- Two new generators of Harary Graphs.
- Added Dual Barabasi-Albert model
- Added VoteRank algorithm
- Added Equitable coloring algorithms
- Added planar layout algorithms
- Les Miserables network example
- Javascript example update

Improvements
------------

- Change default colors to be color-blind friendly
- Many bug fixes and documentation improvements
- Speed up of simple_cycles
- Improvements for reading various formats like GML, GEXF, Graphml
- Allow subclassing to access node_attr_dict_factory


API Changes
-----------
- The G.fresh_copy() mechanism for creating an empty_graph of the same
  type (introduced in v2.0) does not playing nicely with pickle and others.
  So, we have removed the code that caused a need for that. Instead you
  should use the more natural G.__class__() syntax to get an empty_graph
  of the same type as G.

Deprecations
------------
- The Graph.fresh_copy() method should now use Graph.__class__()
- ReverseView class removed in favor of reverse_view() function.

Contributors
------------

- Mike Babst
- Jonathan Barnoud
- Scott Chow
- Jon Crall
- Clayton A Davis
- Michaël Defferrard
- Fredrik Erlandsson
- Eyal
- Tanay Gahlot
- Matthew Gilbert
- Øyvind Heddeland Instefjord
- Hongwei Jin
- Kieran
- Dongkwan Kim
- Julien Klaus
- Warren W. Kretzschmar
- Elias Kuthe
- Eric Ma
- Christoph Martin
- Jarrod Millman
- Issa Moradnejad
- Moradnejad
- Niema Moshiri
- Ramil Nugmanov
- Jens P
- Benjamin Peterson
- Edward L Platt
- Matteo Pozza
- Antoine Prouvost
- Mickaël Schoentgen
- Dan Schult
- Johannes Schulte
- Mridul Seth
- Weisheng Si
- Utkarsh Upadhyay
- damianos
- guidoeco
- jeanfrancois8512
- komo-fr
- last2sword
NetworkX 1.6
============

Release date:  20 November 2011

Highlights
~~~~~~~~~~

New functions for finding articulation points, generating random bipartite graphs, constructing adjacency matrix representations, forming graph products, computing assortativity coefficients, measuring subgraph centrality and communicability, finding k-clique communities, and writing JSON format output.

New examples for drawing with D3 Javascript library, and ordering matrices with the Cuthill-McKee algorithm.

More memory efficient implementation of current-flow betweenness and new approximation algorithms for current-flow betweenness and shortest-path betweenness.

Simplified handling of "weight" attributes for algorithms that use weights/costs/values.

Updated all code to work with the PyPy Python implementation http://pypy.org which produces faster performance on many algorithms.

Graph Classes
-------------

The degree* methods in the graph classes (Graph, DiGraph, MultiGraph,
MultiDiGraph) now take an optional weight= keyword that allows computing
weighted degree with arbitrary (numerical) edge attributes.  Setting
weight=None is equivalent to the previous weighted=False.


Weighted graph algorithms
-------------------------

Many 'weighted' graph algorithms now take optional parameter to
specify which edge attribute should be used for the weight
(default='weight') (ticket https://networkx.lanl.gov/trac/ticket/573)

In some cases the parameter name was changed from weighted, to weight.  Here is
how to specify which edge attribute will be used in the algorithms:

- Use weight=None to consider all weights equally (unweighted case)

- Use weight='weight' to use the 'weight' edge attribute

- Use weight='other' to use the 'other' edge attribute

Algorithms affected are:

to_scipy_sparse_matrix,
clustering,
average_clustering,
bipartite.degree,
spectral_layout,
neighbor_degree,
is_isomorphic,
betweenness_centrality,
betweenness_centrality_subset,
vitality,
load_centrality,
mincost,
shortest_path,
shortest_path_length,
average_shortest_path_length


Isomorphisms
------------

Node and edge attributes are now more easily incorporated into isomorphism
checks via the 'node_match' and 'edge_match' parameters.  As part of this
change, the following classes were removed::

    WeightedGraphMatcher
    WeightedDiGraphMatcher
    WeightedMultiGraphMatcher
    WeightedMultiDiGraphMatcher

The function signature for 'is_isomorphic' is now simply::

    is_isomorphic(g1, g2, node_match=None, edge_match=None)

See its docstring for more details.  To aid in the creation of 'node_match'
and 'edge_match' functions, users are encouraged to work with::

    categorical_node_match
    categorical_edge_match
    categroical_multiedge_match
    numerical_node_match
    numerical_edge_match
    numerical_multiedge_match
    generic_node_match
    generic_edge_match
    generic_multiedge_match

These functions construct functions which can be passed to 'is_isomorphic'.
Finally, note that the above functions are not imported into the top-level
namespace and should be accessed from 'networkx.algorithms.isomorphism'.
A useful import statement that will be repeated throughout documentation is::

    import networkx.algorithms.isomorphism as iso

Other
-----
* attracting_components

  A list of lists is returned instead of a list of tuples.

* condensation

  The condensation algorithm now takes a second argument (scc) and returns a
  graph with nodes labeled as integers instead of node tuples.

* degree connectivity

  average_in_degree_connectivity and average_out_degree_connectivity have
  been replaced with

  average_degree_connectivity(G, source='in', target='in')

  and

  average_degree_connectivity(G, source='out', target='out')

* neighbor degree

  average_neighbor_in_degree and  average_neighbor_out_degreey have
  have been replaced with

  average_neighbor_degree(G, source='in', target='in')

  and

  average_neighbor_degree(G, source='out', target='out')

NetworkX 2.4
============

Release date: 16 October 2019

Supports Python 3.5, 3.6, 3.7, and 3.8.
This is the last release to support Python 3.5.

NetworkX is a Python package for the creation, manipulation, and study of the
structure, dynamics, and functions of complex networks.

For more information, please visit our `website <https://networkx.org/>`_
and our `gallery of examples
<https://networkx.org/documentation/latest/auto_examples/index.html>`_.
Please send comments and questions to the `networkx-discuss mailing list
<http://groups.google.com/group/networkx-discuss>`_.

Highlights
----------

This release is the result of 6 months of work with over 200 commits by
67 contributors. Highlights include:

- Remove deprecated code from 1.x
- Support for Python 3.8
- Switched to pytest for testing
- Last release to support Python 3.5

New Functions:

- barycenter functions
- Bethe Hessian matrix function
- Eulerian Path methods
- group centrality measures
- subgraph monomorphisms
- k-truss algorithms
- onion decomposition
- resistance distance
- asteroidal triples
- non-randomness measures
- linear prufing
- minimum weight bipartite matching
- Incremental closeness centrality
- ISMAGS subgraph isomorphism algorithm
- create chordal graph of a graph

New generators

- Binomial tree generator
- Directed joint degree generator
- Random internet AS graph generator

New for Layouts

- spiral node layout routine
- support for 3d layouts


Improvements
------------
- allow average shortest path to use Floyd-Warshall method
- improve read/write of GML, GEXF, GraphML
- allow string or json object as input to jit_graph
- attempt to allow numpy.array input in place of lists in more places
- faster strongly connected components
- faster Floyd-Warshall Optimization
- faster global efficiency
- faster transitive closure
- fix unionfind; betweenness_subset; lexico-topo-sort; A*;
  inverse_line_graph; async label propagation; edgelist reading;
  Gomory-Hu flow method; label_propagation; partial_duplication;
  shell_layout with 1 node in shell; from_pandas_edgelist
- Documentation improvement and fixes


API Changes
-----------

A utility function is_list_of_ints became is_bunch_of_ints
and now tests int(item)==item instead of isinstance(_, int)
This allows e.g. floats whose values are integer.

Added utility make_list_of_ints to convert containers of
integer values to lists of integers


Deprecations
------------

Removed functions (marked as deprecated in NetworkX 2.1):

- attracting_component_subgraphs
- connected_component_subgraphs
- weakly_connected_component_subgraphs
- strongly_connected_component_subgraphs
- biconnected_component_subgraphs
- See docs for component functions for how to get subgraphs.

Graph Object methods removed (marked as deprecated 2.1)

- G.add_path
- G.add_cycle
- G.add_star
- G.nodes_with_selfloops
- G.number_of_selfloops
- G.selfloop_edges
- These are now NetworkX functions, e.g. nx.add_star(G, 5)
- G.node   --> use G.nodes
- G.fresh_copy   --> use G.__class__

Remove old names for graphview functions.

- ReverseView
- SubGraph
- SubMultiGraph
- SubMultiDiGraph
- SubDiGraph
- GraphView
- DiGraphView
- MultiGraphView
- MultiDiGraphView
- MultiReverseView
- Use reverse_view, subgraph_view and generic_graph_view.

Merged PRs
----------

A total of 205 changes have been committed.

- Bump release version
- algorithms/traversal/edgebfs name fix (#3397)
- Add see also links (#3403)
- Add the reference for the Harary graph generators (#3407)
- typo: swap source and target (#3413)
- Fix spring_layout bug with fixed nodes (#3415)
- Move LFR_benchmark to generators (#3411)
- Add barycenter algorithm (#2939)
- Add bethe hessian matrix (#3401)
- Binomial trees generator (#3409)
- Fix edge_color inconsistency with node_color and description. (#3395)
- Adding module for group centrality measures (#3421)
- Improve edgelist See Also (#3423)
- Typo fix (#3424)
- Add doc warning about self-loops for adamic_adar_index (#3427)
- Fix UnionFind set extraction (#3224)
- add required argument to `write_graphml` example (#3429)
- Fix centrality betweeness subset (#3425)
- Add two versions of Simrank similarity (#3222)
- Fixed typo
- Merge pull request #3436 from nandahkrishna/fix-typo-betweenness-centrality-subset-test
- Reorder and complete doc (#3438)
- added topo_order parameter to functions that rely on topological_sort (#3447)
- Implemented subgraph monomorphism (#3435)
- Set seed in random_degree_sequence_graph docstring test (#3451)
- Replace cb.iterable with np.iterable (#3458)
- don't remove ticks of other pyplot axes (#3476)
- Fix typo in "G>raph Modelling Language" (#3468)
- Naive k-truss algorithm implementation. (#3462)
- Adding onion decomposition (#3461)
- New Feature - Resistance Distance (#3385)
- No multigraphs for betweenness (#3454)
- Wheels are python 3 only
- Fix deprecation warning with Python 3.7 (#3487)
- Fix dfs_preorder_nodes docstring saying "edges" instead of "nodes" (#3484)
- Added group closeness and group degree centralities (#3437)
- Fixed incorrect docs (#3495)
- Fixes Issue #3493 - Bug in lexicographical_topological_sort() (#3494)
- AT-free graph recognition (#3377)
- Update introduction.rst (#3504)
- Full join operation and cograph generator (#3503)
- Optimize the strongly connected components algorithm. (#3516)
- Adding non-randomness measures for graphs (#3515)
- Added safeguards (input graph G) for non-randomness measures  (#3526)
- Optimize the strongly connected components algorithm - Take 2 (#3519)
- Small fix for bug found @ issue #3524 (#3529)
- Restore checking PyPy3 (#3514)
- Linear prufer coding (#3535)
- Fix inverse_line_graph. (#3507)
- Fix A* returning wrong solution (#3508)
- Implement minimum weight full matching of bipartite graphs (#3527)
- Get chordal graph for #1054 (#3353)
- Faster transitive closure computation for DAGs (#3445)
- Write mixed-type attributes correctly in write_graphml_lxml (#3536)
- Fixes some edge cases for inverse_line_graph(). (#3538)
- explicitly stated i.j convention in to_numpy_array
- Incremental Closeness Centrality (undirected, unweighted graphs) (#3444)
- Implement ISMAGS subgraph isomorphism algorithm (#3312)
- Fixes bug in networkx.algorithms.community.label_propagation.asyn_lpa_communities (#3545)
- When exporting to GML, write non 32-bit numbers as strings. (#3540)
- Try to bug Fix #3552 (#3554)
- add Directed Joint Degree Graph generator (#3551)
- typo (#3557)
- Fix a few documentation issues for the bipartite algorithm reference (#3555)
- i,j convention in adj mat i/o in relevant funcs
- Merge pull request #3542 from malch2/doc/update
- Add 3.8-dev to travis
- Fix dict iteration for Py3.8
- Ignore other failures for now
- Fix a typo in docstring for get_edge_data (#3564)
- Fix wrong title (#3566)
- Fix typo in doctring (#3568)
- Fix and Improve docstrings in graph.py (#3569)
- Improved graph class selection table (#3570)
- Add spiral layout for graph drawing (#3534)
- #3575 return coordinates of 3d layouts (#3576)
- Handle k==n within the Watts-Strogatz graph generator (#3579)
- Floyd-Warshall Optimization (#3400)
- Use Sphinx 2.2
- Add missing link to asteroidal docs
- Fix Sphinx warnings
- Fix Sphinx latexpdf build
- Updated Contributor list (#3592)
- Prim from list to set (#3512)
- Fix issue 3491 (#3588)
- Make Travis fail on Python 3.8 failures
- Fix test_gexf to handle default serialisation order of the XML attributes
- Remove future imports needed by Py2
- add internet_as_graph generator (#3574)
- remove cyclical references from OutEdgeDataView (#3598)
- Add minimum source and target margin to draw_networkx_edges. (#3390)
- fix to_directed function (#3599)
- Fixes #3573:GEXF output problem (#3606)
- Global efficiency attempt to speed up (#3604)
- Bugfix: Added flexibility in reading values for label and id (#3603)
- Add method floyd-warshall to average_shortest_path_length (#3267)
- Replaced is with == and minor pycodestyle fixes (#3608)
- Fix many documentation based Issues (#3609)
- Resolve many documentation issues (#3611)
- Fixes #3187  transitive_closure now returns self-loops when cycles present (#3613)
- Add support for initializing pagerank_scipy (#3183)
- Add last 7 lines of Gomory-hu algorithm Fixes #3293 (#3614)
- Implemented Euler Path functions (#3399)
- Fix the direction of edges in label_propagation.py (#3619)
- Removed unused import of random module (#3620)
- Fix operation order in partial_duplication_graph (#3626)
- Keep shells with 1 node away from origin in shell_layout (#3629)
- Allow jit_graph to read json string or json object (#3628)
- Fix typo within incode documentation (#3621)
- pycodestyle and update docs for greedy_coloring.py+tests (#3631)
- Add version badges
- Load long description from README
- Add missing code block (#3630)
- Change is_list_of_ints to make_list_of_ints (#3617)
- Handle edgeattr in from_pandas_edgelist when no columns match request (#3634)
- Make draft of release notes for v2.4
- Shift notes from dev to v2.4 filename.
- Use recent pypy
- Test Py 3.8 on macos
- add check of attr type before converting inf/nan in GEXF (#3636)
- Fix sphinx errors And add links to single_source_dijkstra in docs for dijkstra_path/length (#3638)
- Document subgraph_view (#3627)
- First round of pytest fixes
- Use class methods for class setup/teardown
- Have CIs use pytest
- Use class methods for class setup/teardown, cont.
- Do less testing (until we get it working)
- replace idiom from networkx import * in test files
- Fix assert funcs override
- Fix static methods in link_prediction
- Partially fix v2userfunc tests
- Fix graph/digraph tests
- Fix multigraph checks
- Fix multidigraph checks
- Fix test_function checks
- Fix distance_measures tests
- Fix decorators tests
- Fix some raises in test_mst
- Fix clique tests
- Fix yaml tests
- Fix tests in reportviews
- Fix vf2 tests
- Fix mst tests
- Fix gdal tests
- Convert nose.tools.assert_* functions into asserts
- Remove unused imports
- Fix some warnings
- Update testing instructions
- Reenable all test platforms
- Fix some __init__ warnings
- replace nose yield tests in test_coloring.py
- Add testing, coverage, and dev environment info
- Try pytestimportorskip
- Another pair of variations on pytest.importorskip
- fix typo and try again
- Remove deprecated weakly_connected_component_subgraphs
- replace assert_almost_equal and raises in algorithms/tests
- set places=0 on tests that use old almost_equal
- Update nx.test()
- Have pytest run doctests / not sphinx
- Revert "Remove deprecated weakly_connected_component_subgraphs"
- remove warnings for using deprecated function
- Remove deprecated functions and methods. add to release notes.
- Fix subgraph_view testing
- remove tests of deprecated views and fix use of deprecated G.node
- tracking down use of deprecated functions
- Fix deprecated use of add_path/star/cycle
- reduce warnings for deprecated functions
- skirt issues wih raises in test_harmonic
- reduce the number of warnings by removing deprecated functions
- convert_matrix demo of one way to get doctests to work
- Remove deprecated from examples
- Changes to convert_matrix and others that depend on np.matrix
- clean up doctest deprecated code
- More doctest corrections
- Fix examples
- Remove nose from generators
- Remove nose from utils
- Remove nose from classes
- Replace nose.assert_raises with pytest.raises
- Replace nose.raises with pytest.raises context manager
- Replace `eq_`, `ok_` with assert
- Use pytest for doctest
- Highlight switch to pytest in release notes
- Remove `from nose.tools import *`
- Remove nose.tools.SkipTest
- Finalize transition to pytest
- Merge pull request #3639 from stefanv/pytest-port
- Test Python 3.8 with AppVeyor
- Merge pull request #3648 from jarrodmillman/windows-py3.8
- Remove deprecated weakly_connected_component_subgraphs
- Update release notes
- Update README
- Announce Python 3.8 support
- Designate 2.4rc1 release
- Bump release version
- Remove remaining SkipTests
- fix documentation notes (#3644) (#3645)
- Test Py 3.8.0 on AppVeyor
- Speed up AppVeyor
- Cleanup travis config
- Improve CI caching
- Update Py 3.8 on travis
- Merge pull request #3652 from jarrodmillman/speedup-appveyor
- Finalize release notes

It contained the following 5 merges:

- Fixed typo in betweenness centrality subset test (#3436)
- explicitly stated i.j convention in to_numpy_array (#3542)
- pytest port (#3639)
- Test Python 3.8 with AppVeyor (#3648)
- Cleanup and speedup CI (#3652)

Contributors
------------

- Rajendra Adhikari
- Antoine Allard
- Antoine
- Salim BELHADDAD
- Luca Baldesi
- Tamás Bitai
- Tobias Blass
- Malayaja Chutani
- Peter Cock
- Almog Cohen
- Diogo Cruz
- Martin Darmüntzel
- Elan Ernest
- Jacob Jona Fahlenkamp
- Michael Fedell
- Andy Garfield
- Ramiro Gómez
- Haakon
- Alex Henrie
- Steffen Hirschmann
- Martin James McHugh III
- Jacob
- Søren Fuglede Jørgensen
- Omer Katz
- Julien Klaus
- Matej Klemen
- Nanda H Krishna
- Peter C Kroon
- Anthony Labarre
- Anton Lodder
- MCer4294967296
- Eric Ma
- Fil Menczer
- Erwan Le Merrer
- Alexander Metz
- Jarrod Millman
- Subhendu Ranajn Mishra
- Jamie Morton
- James Myatt
- Kevin Newman
- Aaron Opfer
- Aditya Pal
- Pascal-Ortiz
- Peter
- Jose Pinilla
- Alexios Polyzos
- Michael Recachinas
- Efraim Rodrigues
- Adam Rosenthal
- Dan Schult
- William Schwartz
- Weisheng Si
- Kanishk Tantia
- Ivan Tham
- George Valkanas
- Stefan van der Walt
- Hsi-Hsuan Wu
- Haochen Wu
- Xiangyu Xu
- Jean-Gabriel Young
- bkief
- daniel-karl
- michelb7398
- mikedeltalima
- nandahkrishna
- skhiuk
- tbalint
NetworkX 2.5
============

Release date: 22 August 2020

Supports Python 3.6, 3.7, and 3.8.

NetworkX is a Python package for the creation, manipulation, and study of the
structure, dynamics, and functions of complex networks.

For more information, please visit our `website <https://networkx.org/>`_
and our `gallery of examples
<https://networkx.org/documentation/latest/auto_examples/index.html>`_.
Please send comments and questions to the `networkx-discuss mailing list
<http://groups.google.com/group/networkx-discuss>`_.

Highlights
----------

This release is the result of 10 months of work with over 200 commits by
92 contributors. Highlights include:

- Dropped support for Python 3.5.
- add Pathlib support to work with files.
- improve performance.
- Updated docs and tests.
- Removed code designed to work with Python 2.

New Functions:

- lukes_partitioning
- triadic analysis functions
- functions for trophic levels analysis
- d_separated
- is_regular and other regular graph measures
- graph_hash using Weisfeiler Lehman methods
- common_neighbor_centrality (CCPA link prediction)
- max_weight_clique
- path_weight and is_path
- rescale_layout_dict
- junction_tree

New generators:

- paley_graph
- interval_graph

New layouts:

- multipartite_layout


Improvements
------------

- Add governance documents, developer guide and community structures
- Implement explicit deprecation policy.
- Initiate an NX Enhancement Proposal (NXEP) system
- optimize single_source_shortest_path
- improved consistent "weight" specification in shortest_path routines
- Reduce numpy.matrix usage which is discouraged by numpy.
- improved line color
- better search engine treatment of docs
- lattice and grid_graph and grid_2d_graph can use dim=tuple
- fix initializer of kamada_kawai_layout algorithm
- moral and threshold functions now included in namespace and docs
- scale arrows better when drawing
- more uniform creation of random lobster graphs
- allow editing graph during iteration over connected_components
- better column handling in coversion of pandas DataFrame
- allow simrank_similarity with directed graph input
- ensure VoteRank ability is nonnegative
- speedup kernighan_lin_bisection
- speedup negative weight cycle detection
- tree_isomorphism
- rooted_tree_isomorphism
- Gexf edge attribute "label" is available


API Changes
-----------

- enabled "copy" flag parameter in `contracted_nodes`
- allow partially periodic lattices
- return value for minimum_st_node_cut now always a set
- removed unused "has_numpy" argument from create_py_random_state
- fixed return values when drawing empty nodes and edges
- allow sets and frozensets of edges as input to nx.Graph()
- "weight" can be function for astar, directional_dijksta, all_shortest_path
- allow named key ids for GraphML edge writing
- all keywords are now checked for validity in nx.draw and friends
- EdgeDataView "in" operator checks if nodes are "in nbunch"
- remove completeness condition from minimum weight full matching
- option to sort neighbors in bfs traversal
- draw_networkx accepts numpy array for edgelist
- relabel_nodes with 2 nodes mapped to same node can now create multiedge
- steiner_tree works with MultiGraph
- Add `show` kwarg to view_pygraphviz (#4155)
- Prepare for turning chordal_graph_cliques into a generator (#4162)
- GraphML reader keyword force_multigraph creates MultiGraph even w/o multiedges


Deprecations
------------

- [`#3680 <https://github.com/networkx/networkx/pull/3680>`_]
  Deprecate `make_str(x)` for `str(x)`.
  Deprecate `is_string_like(obj)` for `isinstance(obj, str)`.

- [`#3725 <https://github.com/networkx/networkx/pull/3725>`_]
  Deprecate `literal_stringizer` and `literal_destringizer`.

- [`#3983 <https://github.com/networkx/networkx/pull/3983>`_]
  Deprecate `reversed` context manager.

- [`#4155 <https://github.com/networkx/networkx/pull/4155>`_]
  Deprecate `display_pygraphviz`.

- [`#4162 <https://github.com/networkx/networkx/pull/4162>`_]
  Deprecate `chordal_graph_cliques` returning a set.

- [`#4161 <https://github.com/networkx/networkx/pull/4161>`_]
  Deprecate `betweenness_centrality_source`.

- [`#4161 <https://github.com/networkx/networkx/pull/4161>`_]
  Deprecate `edge_betweeness`.

- [`#4161 <https://github.com/networkx/networkx/pull/4161>`_]
  Rename `_naive_greedy_modularity_communities` as `naive_greedy_modularity_communities`.

Merged PRs
----------

A total of 256 changes have been committed.

- Bump release version
- Update release process
- Drop support for Python 3.5
- fix typo docs
- Remove old Python 2 code
- Enable more doctests
- Fix pydot tests
- Unclear how to test the test helper function
- Pathlib introduced in Py 3.4
- Remove code using sys.version_info to detect Python 2
- Use yield from
- PEP8 fixes to tests
- Remove unused imports
- Use pytest.importorskip
- PEP8 fixes
- Remove unused imports
- Add pep8_speaks conf
- Use itertools accumulate
- Fixes issue 3610: Bug in version attribute of gexf.py
- Ignore W503
- Run doctest without optional dependencies
- Skip doctests when missing dependencies
- Remove sed imports
- Enable tests (#3678)
- `contracted_nodes` copy flag added (#3646)
- Deprecate make_str
- Deprecate is_string_like
- Fix PEP8 issues
- Enable ThinGraph tests (#3681)
- Optimize _single_shortest_path_length (#3647)
- Fix issue 3431: Return error in case of bad input to make_small_graph (#3676)
- avoid duplicate tests due to imports (#3684)
- Fix typo: Laplacion -> Laplacian (#3689)
- Add tests
- Lukes algorithm implementation (#3666)
- Remove shim that worked around using starmap
- Add back to gallery
- Add colormap and color limits to LineCollection (#3698)
- Fix matplotlib deprecation (#3697)
- Adapt SciPy CoC
- Update docs to be more accurate about speed of G.neighbors (#3699)
- Use canonical url to help search engines
- Remove duplicate license parameter (#3710)
- Fix documentation issues for exceptions in a few places
- Fix more documentation issues with exceptions
- Remove old Python 2 code
- Remove boiler plate from top of modules
- Remove superfluous encoding information
- Update examples
- Simplify package docstring
- Remove shebang from non-executables
- Add contributors
- K-truss is defined for edges being in (k-2) triangles and not for k triangles (#3713)
- Enable optional tests on Python 3.8
- Fix test_numpy_type to pass under Python 3.8
- Add links to data files
- Deprecate Python 2/3 compatibility code
- Update style
- Update style
- Separate easy and hard to install optional requirements
- Install optional dependencies by default
- Refactor tests
- Sample code for subgraph copy: add parenthesis to is_multigraph (#3734)
- Fixed typo (#3735)
- fix citation links (#3741)
- remove f strings from setup.py for clear error message < py3.6 (#3738)
- 3511 gml list support (#3649)
- added linestyle as argument (#3747)
- Link to files needed for example (#3752)
- fixed a typo
- Merge pull request #3759 from yohm/patch-1
- remove unused variable so grid_graph supports dim=tuple (#3760)
- Sudoku generator issue 3756 (#3757)
- Fix scaling of single node shells in shall_layout (#3764)
- Adding triadic analysis functions (#3742)
- Improve test coverage
- Update contribs script
- Convert %-format to fstring
- Upgrade to Py36 syntax
- Upgrade to Py36 syntax
- Update string format
- Fix scipy deprecation warnings
- Update year
- Silence known warnings (#3770)
- Fix docstring for asyn_fluidc (#3779)
- Fix #3703 (#3784)
- fix initializer for kamada_kawai_layout (networkx #3658) (#3782)
- Minor comments issue (#3787)
- Adding moral and threshold packages to main namespace (#3788)
- Add weight functions to bidirectional_dijkstra and astar (#3799)
- Shrink the source side of an arrow properly when drawing a directed edge. #3805 (#3806)
- option for partially-periodic lattices (networkx #3586) (#3807)
- Prevent KeyError on subgraph_is_monomorphic (#3798)
- Trophic Levels #3736 (#3804)
- UnionFind's union doesn't accurately track set sizes (#3810)
- Remove whitespace (#3816)
- reconsider the lobster generator (#3822)
- Fix typo (#3838)
- fix typo slightly confusing the meaning (#3840)
- Added fix for issue #3846 (#3848)
- Remove unused variable has_numpy from create_py_random_state (#3852)
- Fix return values when drawing empty nodes and edges  #3833 (#3854)
- Make connected_components safe to component set mutation (#3859)
- Fix example in docstring (#3866)
- Update README.rst website link to https (#3888)
- typo (#3894)
- Made CONTRIBUTING.rst more clearer (#3895)
- Fixing docs for nx.info(), along with necessary tests (#3893)
- added default arg for json dumps for jit_data func (#3891)
- Fixed nx.Digraph to nx.DiGraph (#3909)
- Use Sphinx 3.0.1
- Fix Sphinx deprecation
- Add logo to docs
- allow set of edge nodes (#3907)
- Add extra information when casting 'id' to int() fails. (Resolves #3910) (#3916)
- add paley graph (#3900)
- add paley graph to doc (#3927)
- Update astar.py (#3947)
- use keywords for positional arguments (#3952)
- fix documentation (#3959)
- Add option for named key ids to GraphML writing. (#3960)
- fix documentation (#3958)
- Correct handling of zero-weight edges in all_shortest_paths (#3783)
- Fix documentation typo (#3965)
- Fix: documentation of simrank_similarity_numpy (#3954)
- Fix for #3930 (source & target columns not overwritten when converting to pd.DataFrame) (#3935)
- Add weight function for shortest simple paths for #3948 (#3949)
- Fix defination of communicability (#3973)
- Fix simrank_similarity with directed graph input (#3961)
- Fixed weakening of voting ability (#3970)
- implemented faster sweep algorithm for kernighan_lin_bisection (#3858)
- Fix issue #3926 (#3928)
- Update CONTRIBUTORS.rst (#3982)
- Deprecate context_manager reversed in favor of reversed_view (#3983)
- Update CONTRIBUTORS.rst (#3987)
- Enhancement for voterank (#3972)
- add d-separation algorithm (#3974)
- DOC: added see also section to find_cycle (#3999)
- improve docs for subgraph_view filter_egde (#4010)
- Fix exception causes in dag.py (#4000)
- use raise from for exceptions in to_networkx_graph (#4009)
- Fix exception causes and messages in 12 modules (#4012)
- Fix typo: `np.int` -> `np.int_` (#4013)
- fix a typo (#4017)
- change documentation (#3981)
- algorithms for regular graphs (#3925)
- Typo Hand should be Hans (#4025)
- DOC: Add testing bullet to CONTRIBUTING. (#4035)
- Update Sphinx
- Update optional/test deps
- Add governance/values/nexp/roadmap
- Improve formatting of None in tutorial (#3986)
- Fixes DiGraph spelling in docstring (#3892)
- Update links to Py3 docs (#4042)
- Add method to clear edges only (#3477)
- Fix exception causes and messages all over the codebase (#4015)
- Handle kwds explicitly in draw_networkx (#4033)
- return empty generator instead of empty list (#3967)
- Correctly infer numpy float types (#3919)
- MAINT: Update from_graph6_bytes arg/docs. (#4034)
- Add URLs/banner/titlebar to documentation (#4044)
- Add negative cycle detection heuristic (#3879)
- Remove unused imports (#3855)
- Fixed Bug in generate_gml(G, stringizer=None) (#3841)
- Raise NetworkXError when k < 2 (#3761)
- MAINT: rm np.matrix from alg. conn. module
- MAINT: rm np.matrix from attribute_ac.
- MAINT,TST: Parametrize methods in TestAlgebraicConnectivity.
- MAINT,TST: parametrize buckminsterfullerene test.
- MAINT,TST: Remove unused _methods class attr
- MAINT,TST: Parametrize TestSpectralOrdering.
- excluded self/recursive edges  (#4037)
- WIP: Change EdgeDataView __contains__ feature (2nd attempt) (#3845)
- Index edges for multi graph simple paths (#3358)
- ENH: Add new graph_hashing feature
- Fix pandas deprecation
- Organize removal of deprecated code
- Update sphinx
- ENH: Add roots and timeout to GED (#4026)
- Make gallery more prominent
- Add an implementation for interval_graph and its unit tests (#3705)
- Fixed typo in kamada_kawai_layout docstring (#4059)
- Remove completeness condition from minimum weight full matching (#4057)
- Implemented multipartite_layout (#3815)
- added new Link Prediction algorithm (CCPA) (#4028)
- add the option of sorting node's neighbors during bfs traversal  (#4029)
- TST: remove int64 specification from test. (#4055)
- Ran pyupgrade --py36plus
- Remove trailing spaces
- Tell psf/black to ignore specific np.arrays
- Format w/ black
- Add pre-commit hook to for psf/black
- Merge pull request #4060 from jarrodmillman/black
- Fix a few typos in matching docstrings (#4063)
- fix bug for to_scipy_sparse_matrix function (#3985)
- Update documentation of minimum weight full matching (#4062)
- Add maximum weight clique algorithm (#4016)
- Clear pygraphviz object after creating networkx object (#4070)
- Use newer osx on travis (#4075)
- Install Python after updating brew (#4079)
- Add link to black (#4078)
- Improves docs regarding aliases of erdos-reyni graph generators (#4074)
- MAINT: Remove dependency version info from INSTALL (#4081)
- Simplify top-level directory (#4087)
- DOC: Fix return types in laplacianmatrix. (#4090)
- add modularity to the docs (#4096)
- Allow G.remove_edges_from(nx.selfloops_edges(G)) (#4080)
- MAINT: rm private fn in favor of numpy builtin. (#4094)
- Allow custom keys for multiedges in from_pandas_edgelist (#4076)
- Fix planar_layout docstring (#4097)
- DOC: Rewording re: numpy.matrix
- MAINT: rm to/from_numpy_matrix internally
- Merge pull request #4093 from rossbar/rm_npmatrix
- Remove copyright boilerplate (#4105)
- Update contributor guide (#4088)
- Add function to calculate path cost for a specified path (#4069)
- Update docstring for from_pandas_edgelist (#4108)
- Add max_weight_clique to doc (#4110)
- Update deprecation policyt (#4112)
- Improve modularity calculation (#4103)
- Add team gallery (#4117)
- CI: Setup circle CI for documentation builds (#4119)
- Build pdf (#4123)
- DOC: Suggestions and improvments from tutorial readthrough (#4121)
- Enable 3.9-dev on travis (#4124)
- Fix parse_edgelist behavior with multiple attributes (#4125)
- CI: temporary fix for CI latex installation issues (#4131)
- Updated draw_networkx to accept numpy array for edgelist (#4132)
- Add tree isomorphism (#4067)
- MAINT: Switch to abc-based isinstance checks in to_networkx_graph (#4136)
- Use dict instead of OrderedDict since dict is ordered by default from Python 3.6. (#4145)
- MAINT: fixups to parse_edgelist. (#4128)
- Update apt-get on circleci image (#4147)
- add rescale_layout_dict to change scale of the layout_dicts (#4154)
- Update dependencies
- Remove gdal from requirements
- relabel_nodes now preserves edges in multigraphs (#4066)
- MAINT,TST: Improve coverage of nx_agraph module (#4156)
- Get steiner_tree to work with MultiGraphs by postprocessing (#4160)
- junction_tree for #1012 (#4004)
- API: Add `show` kwarg to view_pygraphviz. (#4155)
- Prepare for turning chordal_graph_cliques into a generator (#4162)
- Docs update (#4161)
- Remove unnecessary nx imports from doctests (#4163)
- MultiGraph from graphml with explicit edge ids #3470 (#3763)
- Update sphinx dep (#4164)
- Add edge label in GEXF writer as an optional attribute (#3347)
- First Draft of Release Notes for v2.5 (#4159)
- Designate 2.5rc1 release
- Bump release version
- Update deprecations in release notes (#4166)
- DOC: Update docstrings for public functions in threshold module (#4167)
- Format python in docstrings (#4168)
- DOC,BLD: Fix doc build warning from markup error. (#4174)

It contained the following 3 merges:

- fixed a typo (#3759)
- Use psf/black (#4060)
- MAINT: Replace internal usage of to_numpy_matrix and from_numpy_matrix (#4093)


Contributors
------------

- Adnan Abdulmuttaleb
- Abhi
- Antoine-H
- Salim BELHADDAD
- Ross Barnowski
- Lukas Bernwald
- Isaac Boates
- Kelly Boothby
- Matthias Bruhns
- Mahmut Bulut
- Rüdiger Busche
- Gaetano Carpinato
- Nikos Chan
- Harold Chan
- Camden Cheek
- Daniel
- Daniel-Davies
- Bastian David
- Christoph Deil
- Tanguy Fardet
- 赵丰 (Zhao Feng)
- Andy Garfield
- Oded Green
- Drew H
- Alex Henrie
- Kang Hong Jin
- Manas Joshi
- Søren Fuglede Jørgensen
- Aabir Abubaker Kar
- Folgert Karsdorp
- Suny Kim
- Don Kirkby
- Katherine Klise
- Steve Kowalik
- Ilia Kurenkov
- Whi Kwon
- Paolo Lammens
- Zachary Lawrence
- Sanghack Lee
- Anton Lodder
- Lukas Lösche
- Eric Ma
- Mackyboy12
- Christoph Martin
- Alex Marvin
- Mattwmaster58
- James McDermott
- Jarrod Millman
- Ibraheem Moosa
- Yohsuke Murase
- Neil
- Harri Nieminen
- Danny Niquette
- Carlos G. Oliver
- Juan Orduz
- Austin Orr
- Pedro Ortale
- Aditya Pal
- PalAditya
- Jose Pinilla
- PranayAnchuri
- Jorge Martín Pérez
- Pradeep Reddy Raamana
- Ram Rachum
- David Radcliffe
- Federico Rosato
- Tom Russell
- Craig Schmidt
- Jonathan Schneider
- Dan Schult
- Mridul Seth
- Karthikeyan Singaravelan
- Songyu-Wang
- Kanishk Tantia
- Jeremias Traub
- James Trimble
- Shashi Tripathi
- Stefan van der Walt
- Jonatan Westholm
- Kazimierz Wojciechowski
- Jangwon Yie
- adnanmuttaleb
- anentropic
- arunwise
- beckedorf
- ernstklrb
- farhanbhoraniya
- fj128
- gseva
- haochenucr
- johnthagen
- kiryph
- muratgu
- ryan-duve
- sauxpa
- tombeek111
- willpeppo
NetworkX 1.5
============

Release date:  4 June 2011

Highlights
~~~~~~~~~~

New features
~~~~~~~~~~~~
 - Algorithms for :mod:`generating <networkx.generators.bipartite>`
   and :mod:`analyzing <networkx.algorithms.bipartite>` bipartite graphs
 - :mod:`Maximal independent set <networkx.algorithms.mis>` algorithm
 - :mod:`Erdős-Gallai graphical degree sequence test <networkx.generators.degree_seq>`
 - :mod:`Negative edge cycle test <networkx.algorithms.shortest_paths.weighted>`
 - More memory efficient :mod:`Dijkstra path length <networkx.algorithms.shortest_paths.weighted>` with cutoff parameter
 - :mod:`Weighted clustering coefficient <networkx.algorithms.cluster>`
 - Read and write version 1.2 of :mod:`GEXF reader <networkx.readwrite.gexf>` format
 - :mod:`Neighbor degree correlation <networkx.algorithms.neighbor_degree>`
   that handle subsets of nodes
 - :mod:`In-place node relabeling <networkx.relabel>`
 - Many 'weighted' graph algorithms now take optional parameter to use
   specified edge attribute (default='weight')
   (ticket https://networkx.lanl.gov/trac/ticket/509)

 - Test for :mod:`distance regular <networkx.algorithms.distance_regular>` graphs
 - Fast :mod:`directed Erdős-Renyi graph  <networkx.generators.random_graphs>` generator
 - Fast :mod:`expected degree graph  <networkx.generators.degree_seq>` generator
 - :mod:`Navigable small world  <networkx.generators.geometric>` generator
 - :mod:`Waxman model <networkx.generators.geometric>` generator
 - :mod:`Geographical threshold graph <networkx.generators.geometric>` generator
 - :mod:`Karate Club, Florentine Families, and Davis' Women's Club <networkx.generators.social>` graphs

Weighted graph algorithms
-------------------------

Many 'weighted' graph algorithms now take optional parameter to
specify which edge attribute should be used for the weight
(default='weight') (ticket https://networkx.lanl.gov/trac/ticket/509)

In some cases the parameter name was changed from weighted_edges,
or weighted, to weight.  Here is how to specify which edge attribute
will be used in the algorithms:

- Use weight=None to consider all weights equally (unweighted case)

- Use weight=True or weight='weight' to use the 'weight' edge attribute

- Use weight='other' to use the 'other' edge attribute

Algorithms affected are:

betweenness_centrality, closeness_centrality, edge_bewteeness_centrality,
betweeness_centrality_subset, edge_betweenness_centrality_subset,
betweenness_centrality_source, load, closness_vitality,
weiner_index, spectral_bipartivity
current_flow_betweenness_centrality,
edge_current_flow_betweenness_centrality,
current_flow_betweenness_centrality_subset,
edge_current_flow_betweenness_centrality_subset,
laplacian, normalized_laplacian, adj_matrix, adjacency_spectrum,
shortest_path, shortest_path_length, average_shortest_path_length,
single_source_dijkstra_path_basic, astar_path, astar_path_length

Random geometric graph
----------------------

The random geometric graph generator has been simplified.
It no longer supports the create_using, repel, or verbose parameters.
An optional pos keyword was added to allow specification of node positions.

Bug fixes
~~~~~~~~~
 - Fix edge handling for multigraphs in networkx/graphviz interface
   (ticket https://networkx.lanl.gov/trac/ticket/507)
 - Update networkx/pydot interface for new versions of pydot
   (ticket https://networkx.lanl.gov/trac/ticket/506)
   (ticket https://networkx.lanl.gov/trac/ticket/535)
 - Fix negative cycle handling in Bellman-Ford
   (ticket https://networkx.lanl.gov/trac/ticket/502)
 - Write more attributes with GraphML and GML formats
   (ticket https://networkx.lanl.gov/trac/ticket/480)
 - Handle white space better in read_edgelist
   (ticket https://networkx.lanl.gov/trac/ticket/513)
 - Better parsing of Pajek format files
   (ticket https://networkx.lanl.gov/trac/ticket/524)
   (ticket https://networkx.lanl.gov/trac/ticket/542)
 - Isolates functions work with directed graphs
   (ticket https://networkx.lanl.gov/trac/ticket/526)
 - Faster conversion to numpy matrices
   (ticket https://networkx.lanl.gov/trac/ticket/529)
 - Add graph['name'] and use properties to access Graph.name
   (ticket https://networkx.lanl.gov/trac/ticket/544)
 - Topological sort confused None and 0
   (ticket https://networkx.lanl.gov/trac/ticket/546)
 - GEXF writer mishandled weight=0
   (ticket https://networkx.lanl.gov/trac/ticket/550)
 - Speedup in SciPy version of PageRank
   (ticket https://networkx.lanl.gov/trac/ticket/554)
 - Numpy PageRank node order incorrect + speedups
   (ticket https://networkx.lanl.gov/trac/ticket/555)
NetworkX 2.2
============

Release date: 19 September 2018

Supports Python 2.7, 3.5, 3.6 and 3.7.
This is the last release to support Python 2.

NetworkX is a Python package for the creation, manipulation, and study of the
structure, dynamics, and functions of complex networks.

For more information, please visit our `website <https://networkx.org/>`_
and our `gallery of examples
<https://networkx.org/documentation/latest/auto_examples/index.html>`_.
Please send comments and questions to the `networkx-discuss mailing list
<http://groups.google.com/group/networkx-discuss>`_.

Highlights
----------

This release is the result of 8 months of work with over 149 commits by
58 contributors. Highlights include:

- Add support for Python 3.7. This is the last release to support Python 2.
- Uniform random number generator (RNG) handling which defaults to global
  RNGs but allows specification of a single RNG for all random numbers in NX.
- Improved GraphViews to ease subclassing and remove cyclic references
  which caused trouble with deepcopy and pickle.
- New Graph method `G.update(H)`

Improvements
------------

Each function that uses random numbers now uses a `seed` argument to control
the random number generation (RNG). By default the global default RNG is
used. More precisely, the `random` package's default RNG or the numpy.random
default RNG. You can also create your own RNG and pass it into the `seed`
argument. Finally, you can use an integer to indicate the state to set for
the RNG. In this case a local RNG is created leaving the global RNG untouched.
Some functions use `random` and some use `numpy.random`, but we have written
a translater so that all functions CAN take a `numpy.random.RandomState`
object. So a single RNG can be used for the entire package.

Cyclic references between graph classes and views have been removed to ease
subclassing without memory leaks. Graphs no longer hold references to views.

Cyclic references between a graph and itself have been removed by eliminating
G.root_graph. It turns out this was an avoidable construct anyway.

GraphViews have been reformulated as functions removing much of the subclass
trouble with the copy/to_directed/subgraph methods. It also simplifies the
graph view code base and API. There are now three function that create
graph views: generic_graph_view(graph, create_using), reverse_view(digraph)
and subgraph_view(graph, node_filter, edge_filter).

GraphML can now be written with attributes using numpy numeric types.
In particular, np.float64 and np.int64 no longer need to convert to Python
float and int to be written. They are still written as generic floats so
reading them back in will not make the numpy values.

A generator following the Stochastic Block Model is now available.

New function `all_topolgical_sort` to generate all possible top_sorts.

New functions for tree width and tree decompositions.

Functions for Clauset-Newman-Moore modularity-max community detection.

Functions for small world analysis, directed clustering and perfect matchings,
eulerizing a graph, depth-limited BFS, percolation centrality,
planarity checking.

The shortest_path generic and convenience functions now have a `method`
parameter to choose between dijkstra and bellmon-ford in the weighted case.
Default is dijkstra (which was the only option before).

API Changes
-----------
empty_graph has taken over the functionality from
nx.convert._prep_create_using which was removed.

The `create_using` argument (used in many functions) should now be a
Graph Constructor like nx.Graph or nx.DiGraph.
It can still be a graph instance which will be cleared before use, but the
preferred use is a constructor.

New Base Class Method: update
H.update(G) adds the nodes, edges and graph attributes of G to H.
H.update(edges=e, nodes=n) add the edges and nodes from containers e and n.
H.update(e), and H.update(nodes=n) are also allowed.
First argument is a graph if it has `edges` and `nodes` attributes.
Otherwise the first argument is treated as a list of edges.

The bellman_ford predecessor dicts had sentinal value `[None]` for
source nodes. That has been changed so source nodes have pred value '[]'


Deprecations
------------

Graph class method `fresh_copy` - simply use `__class__`.
The GraphView classes are deprecated in preference to the function
interface. Specifically, `ReverseView` and `ReverseMultiView` are
replaced by `reverse_view`. `SubGraph`, `SubDiGraph`, `SubMultiGraph`
and `SubMultiDiGraph` are replaced by `subgraph_view`.
And `GraphView`, `DiGraphView`, `MultiGraphView`, `MultiDiGraphView`
are derecated in favor of `generic_graph_view(graph, create_using)`.


Contributors
------------

- Luca Baldesi
- William Bernoudy
- Alexander Condello
- Saurav Das
- Dormir30
- Graham Fetterman
- Robert Gmyr
- Thomas Grainger
- Benjamin M. Gyori
- Ramiro Gómez
- Darío Hereñú
- Mads Jensen
- Michael Johnson
- Pranay Kanwar
- Aabir Abubaker Kar
- Jacek Karwowski
- Mohammed Kashif
- David Kraeutmann
- Winni Kretzschmar
- Ivan Laković
- Daniel Leicht
- Katrin Leinweber
- Alexander Lenail
- Lonnen
- Ji Ma
- Erwan Le Merrer
- Jarrod Millman
- Baurzhan Muftakhidinov
- Neil
- Jens P
- Edward L Platt
- Guillaume Plique
- Miguel Sozinho Ramalho
- Lewis Robbins
- Romain
- Federico Rosato
- Tom Russell
- Dan Schult
- Gabe Schwartz
- Aaron Smith
- Leo Torres
- Martin Váňa
- Ruaridh Williamson
- Huon Wilson
- Haochen Wu
- Yuto Yamaguchi
- Felix Yan
- Jean-Gabriel Young
- aparamon
- armando1793
- aweltsch
- chebee7i
- hongshaoyang
- komo-fr
- leamingrad
- luzpaz
- mtrenfield
- regstrtn
NetworkX 1.8
============

Release date:  28 July 2013

Highlights
~~~~~~~~~~
- Faster (linear-time) graphicality tests and Havel-Hakimi graph generators
- Directed Laplacian matrix generator
- Katz centrality algorithm
- Functions to generate all simple paths
- Improved shapefile reader
- More flexible weighted projection of bipartite graphs
- Faster topological sort, descendants and ancestors of DAGs
- Scaling parameter for force-directed layout

Bug fixes
~~~~~~~~~
- Error with average weighted connectivity for digraphs, correct normalized laplacian with self-loops, load betweenness for single node graphs, isolated nodes missing from dfs/bfs trees, normalize HITS using l1, handle density of graphs with self loops

- Cleaner handling of current figure status with Matplotlib, Pajek files now don't write troublesome header line, default alpha value for GEXF files, read curved edges from yEd GraphML


For full details of the issues closed for this release (added features and bug fixes) see: https://github.com/networkx/networkx/issues?milestone=1&page=1&state=closed

API changes
~~~~~~~~~~~

* Laplacian functions now all return matrices.  To get a numpy array from a matrix use L = nx.laplacian_matrix(G).A

* is_directed_acyclic_graph() now returns false on undirected graphs (instead of raising exception)

* cycles returned from simple_cycles() do not include repeated last node


.. _generators:


Graph generators
****************

.. currentmodule:: networkx


Atlas
-----
.. automodule:: networkx.generators.atlas
.. autosummary::
   :toctree: generated/

   graph_atlas
   graph_atlas_g


Classic
-------
.. automodule:: networkx.generators.classic
.. autosummary::
   :toctree: generated/

   balanced_tree
   barbell_graph
   binomial_tree
   complete_graph
   complete_multipartite_graph
   circular_ladder_graph
   circulant_graph
   cycle_graph
   dorogovtsev_goltsev_mendes_graph
   empty_graph
   full_rary_tree
   ladder_graph
   lollipop_graph
   null_graph
   path_graph
   star_graph
   trivial_graph
   turan_graph
   wheel_graph


Expanders
---------
.. automodule:: networkx.generators.expanders
.. autosummary::
   :toctree: generated/

   margulis_gabber_galil_graph
   chordal_cycle_graph
   paley_graph

Lattice
-------
.. automodule:: networkx.generators.lattice
.. autosummary::
   :toctree: generated/

   grid_2d_graph
   grid_graph
   hexagonal_lattice_graph
   hypercube_graph
   triangular_lattice_graph


Small
-----
.. automodule:: networkx.generators.small
.. autosummary::
   :toctree: generated/

   make_small_graph
   LCF_graph
   bull_graph
   chvatal_graph
   cubical_graph
   desargues_graph
   diamond_graph
   dodecahedral_graph
   frucht_graph
   heawood_graph
   hoffman_singleton_graph
   house_graph
   house_x_graph
   icosahedral_graph
   krackhardt_kite_graph
   moebius_kantor_graph
   octahedral_graph
   pappus_graph
   petersen_graph
   sedgewick_maze_graph
   tetrahedral_graph
   truncated_cube_graph
   truncated_tetrahedron_graph
   tutte_graph


Random Graphs
-------------
.. automodule:: networkx.generators.random_graphs
.. autosummary::
   :toctree: generated/

   fast_gnp_random_graph
   gnp_random_graph
   dense_gnm_random_graph
   gnm_random_graph
   erdos_renyi_graph
   binomial_graph
   newman_watts_strogatz_graph
   watts_strogatz_graph
   connected_watts_strogatz_graph
   random_regular_graph
   barabasi_albert_graph
   dual_barabasi_albert_graph
   extended_barabasi_albert_graph
   powerlaw_cluster_graph
   random_kernel_graph
   random_lobster
   random_shell_graph
   random_powerlaw_tree
   random_powerlaw_tree_sequence
   random_kernel_graph


Duplication Divergence
----------------------
.. automodule:: networkx.generators.duplication
.. autosummary::
   :toctree: generated/

   duplication_divergence_graph
   partial_duplication_graph


Degree Sequence
---------------
.. automodule:: networkx.generators.degree_seq

.. autosummary::
   :toctree: generated/

   configuration_model
   directed_configuration_model
   expected_degree_graph
   havel_hakimi_graph
   directed_havel_hakimi_graph
   degree_sequence_tree
   random_degree_sequence_graph


Random Clustered
----------------
.. automodule:: networkx.generators.random_clustered

.. autosummary::
   :toctree: generated/

   random_clustered_graph


Directed
--------
.. automodule:: networkx.generators.directed
.. autosummary::
   :toctree: generated/

   gn_graph
   gnr_graph
   gnc_graph
   random_k_out_graph
   scale_free_graph


Geometric
---------
.. automodule:: networkx.generators.geometric
.. autosummary::
   :toctree: generated/

   geometric_edges
   geographical_threshold_graph
   navigable_small_world_graph
   random_geometric_graph
   soft_random_geometric_graph
   thresholded_random_geometric_graph
   waxman_graph

Line Graph
----------
.. automodule:: networkx.generators.line
.. autosummary::
   :toctree: generated/

   line_graph
   inverse_line_graph


Ego Graph
---------
.. automodule:: networkx.generators.ego
.. autosummary::
   :toctree: generated/

   ego_graph


Stochastic
----------
.. automodule:: networkx.generators.stochastic
.. autosummary::
   :toctree: generated/

   stochastic_graph


AS graph
--------
.. automodule:: networkx.generators.internet_as_graphs
.. autosummary::
   :toctree: generated/

   random_internet_as_graph


Intersection
------------
.. automodule:: networkx.generators.intersection
.. autosummary::
   :toctree: generated/

   uniform_random_intersection_graph
   k_random_intersection_graph
   general_random_intersection_graph


Social Networks
---------------
.. automodule:: networkx.generators.social
.. autosummary::
   :toctree: generated/

   karate_club_graph
   davis_southern_women_graph
   florentine_families_graph
   les_miserables_graph


Community
---------
.. automodule:: networkx.generators.community
.. autosummary::
   :toctree: generated/

   caveman_graph
   connected_caveman_graph
   gaussian_random_partition_graph
   LFR_benchmark_graph
   planted_partition_graph
   random_partition_graph
   relaxed_caveman_graph
   ring_of_cliques
   stochastic_block_model
   windmill_graph


Spectral
--------
.. automodule:: networkx.generators.spectral_graph_forge
.. autosummary::
   :toctree: generated/

   spectral_graph_forge


Trees
-----
.. automodule:: networkx.generators.trees
.. autosummary::
   :toctree: generated/

   random_tree
   prefix_tree


Non Isomorphic Trees
--------------------
.. automodule:: networkx.generators.nonisomorphic_trees
.. autosummary::
   :toctree: generated/

   nonisomorphic_trees
   number_of_nonisomorphic_trees


Triads
------
.. automodule:: networkx.generators.triads
.. autosummary::
   :toctree: generated/

   triad_graph


Joint Degree Sequence
---------------------
.. automodule:: networkx.generators.joint_degree_seq
.. autosummary::
   :toctree: generated/

   is_valid_joint_degree
   joint_degree_graph
   is_valid_directed_joint_degree
   directed_joint_degree_graph


Mycielski
---------
.. automodule:: networkx.generators.mycielski
.. autosummary::
   :toctree: generated/

   mycielskian
   mycielski_graph


Harary Graph
------------
.. automodule:: networkx.generators.harary_graph
.. autosummary::
   :toctree: generated/

   hnm_harary_graph
   hkn_harary_graph

Cographs
------------
.. automodule:: networkx.generators.cographs
.. autosummary::
   :toctree: generated/

   random_cograph

Interval Graph
---------------
.. automodule:: networkx.generators.interval_graph
.. autosummary::
   :toctree: generated/

   interval_graph

Sudoku
------
.. automodule:: networkx.generators.sudoku
.. autosummary::
   :toctree: generated/

   sudoku_graph
.. _drawing:

*******
Drawing
*******

NetworkX provides basic functionality for visualizing graphs, but its main goal
is to enable graph analysis rather than perform graph visualization. In the
future, graph visualization functionality may be removed from NetworkX or only
available as an add-on package.

Proper graph visualization is hard, and we highly recommend that people
visualize their graphs with tools dedicated to that task. Notable examples of
dedicated and fully-featured graph visualization tools are
`Cytoscape <http://www.cytoscape.org/>`_,
`Gephi <https://gephi.org/>`_,
`Graphviz <http://www.graphviz.org/>`_ and, for
`LaTeX <http://www.latex-project.org/>`_ typesetting,
`PGF/TikZ <https://sourceforge.net/projects/pgf/>`_.
To use these and other such tools, you should export your NetworkX graph into
a format that can be read by those tools. For example, Cytoscape can read the
GraphML format, and so, ``networkx.write_graphml(G, path)`` might be an appropriate
choice.

More information on the features provided here are available at
 - matplotlib:  http://matplotlib.org/
 - pygraphviz:  http://pygraphviz.github.io/


Matplotlib
==========
.. automodule:: networkx.drawing.nx_pylab
.. autosummary::
   :toctree: generated/

   draw
   draw_networkx
   draw_networkx_nodes
   draw_networkx_edges
   draw_networkx_labels
   draw_networkx_edge_labels
   draw_circular
   draw_kamada_kawai
   draw_planar
   draw_random
   draw_spectral
   draw_spring
   draw_shell



Graphviz AGraph (dot)
=====================
.. automodule:: networkx.drawing.nx_agraph
.. autosummary::
   :toctree: generated/

   from_agraph
   to_agraph
   write_dot
   read_dot
   graphviz_layout
   pygraphviz_layout


Graphviz with pydot
===================
.. automodule:: networkx.drawing.nx_pydot
.. autosummary::
   :toctree: generated/

   from_pydot
   to_pydot
   write_dot
   read_dot
   graphviz_layout
   pydot_layout


Graph Layout
============
.. automodule:: networkx.drawing.layout
.. autosummary::
   :toctree: generated/

   bipartite_layout
   circular_layout
   kamada_kawai_layout
   planar_layout
   random_layout
   rescale_layout
   rescale_layout_dict
   shell_layout
   spring_layout
   spectral_layout
   spiral_layout
   multipartite_layout
   
****************
Relabeling nodes
****************
.. currentmodule:: networkx

Relabeling
----------
.. automodule:: networkx.relabel

.. autosummary::
   :toctree: generated/

   convert_node_labels_to_integers
   relabel_nodes


**********
Exceptions
**********

.. automodule:: networkx.exception
.. currentmodule:: networkx

.. autoclass:: networkx.NetworkXException

.. autoclass:: networkx.NetworkXError

.. autoclass:: networkx.NetworkXPointlessConcept

.. autoclass:: networkx.NetworkXAlgorithmError

.. autoclass:: networkx.NetworkXUnfeasible

.. autoclass:: networkx.NetworkXNoPath

.. autoclass:: networkx.NetworkXNoCycle

.. autoclass:: networkx.NodeNotFound

.. autoclass:: networkx.HasACycle

.. autoclass:: networkx.NetworkXUnbounded

.. autoclass:: networkx.NetworkXNotImplemented

.. autoclass:: networkx.AmbiguousSolution

.. autoclass:: networkx.ExceededMaxIterations

.. autoclass:: networkx.PowerIterationFailedConvergence
*********
Utilities
*********

.. automodule:: networkx.utils
.. currentmodule:: networkx.utils

Helper Functions
----------------
.. automodule:: networkx.utils.misc
.. autosummary::
   :toctree: generated/

   arbitrary_element
   is_string_like
   flatten
   iterable
   make_list_of_ints
   make_str
   generate_unique_node
   default_opener
   pairwise
   groups
   create_random_state
   nodes_equal
   edges_equal
   graphs_equal

Data Structures and Algorithms
------------------------------
.. automodule:: networkx.utils.union_find
.. autosummary::
   :toctree: generated/

   UnionFind.union

Random Sequence Generators
--------------------------
.. automodule:: networkx.utils.random_sequence
.. autosummary::
   :toctree: generated/

   powerlaw_sequence
   cumulative_distribution
   discrete_sequence
   zipf_rv
   random_weighted_sample
   weighted_choice

Decorators
----------
.. automodule:: networkx.utils.decorators
.. autosummary::
   :toctree: generated/

   open_file
   not_implemented_for
   nodes_or_number
   np_random_state
   py_random_state
   argmap

Cuthill-Mckee Ordering
----------------------
.. automodule:: networkx.utils.rcm
.. autosummary::
   :toctree: generated/

   cuthill_mckee_ordering
   reverse_cuthill_mckee_ordering
.. _glossary:

Glossary
========

.. glossary::

   dictionary
      A Python dictionary maps keys to values. Also known as "hashes",
      or "associative arrays" in other programming languages.
      See https://docs.python.org/2/tutorial/datastructures.html#dictionaries

   edge
      Edges are either two-tuples of nodes `(u, v)` or three tuples of nodes
      with an edge attribute dictionary `(u, v, dict)`.

   ebunch
      An iteratable container of edge tuples like a list, iterator,
      or file.

   edge attribute
      Edges can have arbitrary Python objects assigned as attributes
      by using keyword/value pairs when adding an edge
      assigning to the `G.edges[u][v]` attribute dictionary for the
      specified edge *u*-*v*.

   nbunch
      An nbunch is a single node, container of nodes or `None` (representing
      all nodes). It can be a list, set, graph, etc.. To filter an nbunch
      so that only nodes actually in `G` appear, use `G.nbunch_iter(nbunch)`.

   node
      A node can be any hashable Python object except None.

   node attribute
     Nodes can have arbitrary Python objects assigned as attributes
     by using keyword/value pairs when adding a node or
     assigning to the `G.nodes[n]` attribute dictionary for the
     specified node `n`.
.. _randomness:

Randomness
==========
.. currentmodule:: networkx

Random Number Generators (RNGs) are often used when generating, drawing
and computing properties or manipulating networks. NetworkX provides
functions which use one of two standard RNGs: NumPy's package `numpy.random`
or Python's built-in package `random`. They each provide the same
algorithm for generating numbers (Mersenne Twister). Their interfaces
are similar (dangerously similar) and yet distinct.
They each provide a global default instance of their generator that
is shared by all programs in a single session.
For the most part you can use the RNGs as NetworkX has them set up and
you'll get reasonable pseudorandom results (results that are statistically
random, but created in a deterministic manner).

Sometimes you want more control over how the numbers are generated.
In particular, you need to set the `seed` of the generator to make
your results reproducible -- either for scientific publication or
for debugging. Both RNG packages have easy functions to set the seed
to any integer, thus determining the subsequent generated values.
Since this package (and many others) use both RNGs you may need to
set the `seed` of both RNGs.  Even if we strictly only used one of the
RNGs, you may find yourself using another package that uses the other.
Setting the state of the two global RNGs is as simple setting the
seed of each RNG to an arbitrary integer:

.. nbplot::

   >>> import random
   >>> random.seed(246)        # or any integer
   >>> import numpy
   >>> numpy.random.seed(4812)

Many users will be satisfied with this level of control.

For people who want even more control, we include an optional argument
to functions that use an RNG.  This argument is called `seed`, but
determines more than the seed of the RNG. It tells the function which
RNG package to use, and whether to use a global or local RNG.

.. nbplot::

    >>> from networkx import path_graph, random_layout
    >>> G = path_graph(9)
    >>> pos = random_layout(G, seed=None)  # use (either) global default RNG
    >>> pos = random_layout(G, seed=42)  # local RNG just for this call
    >>> pos = random_layout(G, seed=numpy.random)  # use numpy global RNG
    >>> random_state = numpy.random.RandomState(42)
    >>> pos = random_layout(G, seed=random_state)  # use/reuse your own RNG

Each NetworkX function that uses an RNG was written with one RNG package
in mind. It either uses `random` or `numpy.random` by default.
But some users want to only use a single RNG for all their code.
This `seed` argument provides a mechanism so that any function
can use a `numpy.random` RNG even if the function is written for `random`.
It works as follows.

The default behavior (when `seed=None`) is to use the global RNG
for the function's preferred package.
If seed is set to an integer value,
a local RNG is created with the indicated seed value and
is used for the duration of that function (including any
calls to other functions) and then discarded.
Alternatively, you can specify `seed=numpy.random` to ensure that
the global numpy RNG is used whether the function expects it or not.
Finally, you can provide a numpy RNG to be used by the function.
The RNG is then available to use in other functions or even other
package like sklearn.
In this way you can use a single RNG for all random numbers
in your project.

While it is possible to assign `seed` a `random`-style RNG for
NetworkX functions written for the `random` package API,
the numpy RNG interface has too
many nice features for us to ensure a `random`-style RNG will work in
all functions. In practice, you can do most things using only `random`
RNGs (useful if numpy is not available). But your experience will be
richer if numpy is available.

To summarize, you can easily ignore the `seed` argument and use the global
RNGs. You can specify to use only the numpy global RNG with
`seed=numpy.random`. You can use a local RNG by providing an integer
seed value. And you can provide your own numpy RNG, reusing it for all
functions. It is easier to use numpy RNGs if you want a single RNG for
your computations.
.. _linalg:


Linear algebra
**************

.. currentmodule:: networkx

Graph Matrix
------------
.. automodule:: networkx.linalg.graphmatrix
.. autosummary::
   :toctree: generated/

   adjacency_matrix
   incidence_matrix

Laplacian Matrix
----------------
.. automodule:: networkx.linalg.laplacianmatrix
.. autosummary::
   :toctree: generated/

   laplacian_matrix
   normalized_laplacian_matrix
   directed_laplacian_matrix
   directed_combinatorial_laplacian_matrix

Bethe Hessian Matrix
--------------------
.. automodule:: networkx.linalg.bethehessianmatrix
.. autosummary::
   :toctree: generated/

   bethe_hessian_matrix

Algebraic Connectivity
----------------------
.. automodule:: networkx.linalg.algebraicconnectivity
.. autosummary::
   :toctree: generated/

   algebraic_connectivity
   fiedler_vector
   spectral_ordering

Attribute Matrices
------------------

.. automodule:: networkx.linalg.attrmatrix
.. autosummary::
   :toctree: generated/

   attr_matrix
   attr_sparse_matrix

Modularity Matrices
-------------------

.. automodule:: networkx.linalg.modularitymatrix
.. autosummary::
   :toctree: generated/

   modularity_matrix
   directed_modularity_matrix

Spectrum
---------
.. automodule:: networkx.linalg.spectrum
.. autosummary::
   :toctree: generated/

   adjacency_spectrum
   laplacian_spectrum
   bethe_hessian_spectrum
   normalized_laplacian_spectrum
   modularity_spectrum
*****************************************
Converting to and from other data formats
*****************************************
.. currentmodule:: networkx

To NetworkX Graph
-----------------
.. automodule:: networkx.convert

.. autosummary::
   :toctree: generated/

   to_networkx_graph


Dictionaries
------------
.. autosummary::
   :toctree: generated/

   to_dict_of_dicts
   from_dict_of_dicts

Lists
-----
.. autosummary::
   :toctree: generated/

   to_dict_of_lists
   from_dict_of_lists
   to_edgelist
   from_edgelist

Numpy
-----
.. automodule:: networkx.convert_matrix
.. autosummary::
   :toctree: generated/

   to_numpy_matrix
   to_numpy_array
   to_numpy_recarray
   from_numpy_matrix
   from_numpy_array

Scipy
-----
.. autosummary::
   :toctree: generated/

   to_scipy_sparse_array
   to_scipy_sparse_matrix
   from_scipy_sparse_matrix

Pandas
------
.. autosummary::
   :toctree: generated/

   to_pandas_adjacency
   from_pandas_adjacency
   to_pandas_edgelist
   from_pandas_edgelist
*********
Functions
*********

.. automodule:: networkx.classes.function

Graph
-----
.. autosummary::
   :toctree: generated/

   degree
   degree_histogram
   density
   info
   create_empty_copy
   is_directed
   to_directed
   to_undirected
   is_empty
   add_star
   add_path
   add_cycle
   subgraph
   subgraph_view
   induced_subgraph
   restricted_view
   reverse_view
   edge_subgraph


Nodes
-----
.. autosummary::
   :toctree: generated/

   nodes
   number_of_nodes
   neighbors
   all_neighbors
   non_neighbors
   common_neighbors


Edges
-----
.. autosummary::
   :toctree: generated/

   edges
   number_of_edges
   density
   non_edges

Self loops
----------
.. autosummary::
   :toctree: generated/

   selfloop_edges
   number_of_selfloops
   nodes_with_selfloops

Attributes
----------
.. autosummary::
   :toctree: generated/

   is_weighted
   is_negatively_weighted
   set_node_attributes
   get_node_attributes
   set_edge_attributes
   get_edge_attributes

Paths
----------
.. autosummary::
   :toctree: generated/

   is_path
   path_weight

Freezing graph structure
------------------------
.. autosummary::
   :toctree: generated/

   freeze
   is_frozen
.. _reference:

Reference
*********

   :Release: |release|
   :Date: |today|

.. toctree::
   :maxdepth: 2

   introduction
   classes/index
   algorithms/index
   functions
   generators
   linalg
   convert
   relabel
   readwrite/index
   drawing
   randomness
   exceptions
   utils
   glossary
Introduction
============

.. currentmodule:: networkx

The structure of NetworkX can be seen by the organization of its source code.
The package provides classes for graph objects, generators to create standard
graphs, IO routines for reading in existing datasets, algorithms to analyze
the resulting networks and some basic drawing tools.

Most of the NetworkX API is provided by functions which take a graph object
as an argument.  Methods of the graph object are limited to basic manipulation
and reporting.  This provides modularity of code and documentation.
It also makes it easier for newcomers to learn about the package in stages.
The source code for each module is meant to be easy to read and reading
this Python code is actually a good way to learn more about network algorithms,
but we have put a lot of effort into making the documentation sufficient and friendly.
If you have suggestions or questions please contact us by joining the
`NetworkX Google group <http://groups.google.com/group/networkx-discuss>`_.

Classes are named using ``CamelCase`` (capital letters at the start of each word).
functions, methods and variable names are ``lower_case_underscore`` (lowercase with
an underscore representing a space between words).


NetworkX Basics
---------------

After starting Python, import the networkx module with (the recommended way)

.. nbplot::

   >>> import networkx as nx

To save repetition, in the documentation we assume that
NetworkX has been imported this way.

If importing networkx fails, it means that Python cannot find the installed
module. Check your installation and your ``PYTHONPATH``.

The following basic graph types are provided as Python classes:

:class:`Graph`
   This class implements an undirected graph. It ignores
   multiple edges between two nodes.  It does allow self-loop
   edges between a node and itself.

:class:`DiGraph`
   Directed graphs, that is, graphs with directed edges.
   Provides operations common to directed graphs,
   (a subclass of Graph).

:class:`MultiGraph`
   A flexible graph class that allows multiple undirected edges between
   pairs of nodes.  The additional flexibility leads to some degradation
   in performance, though usually not significant.

:class:`MultiDiGraph`
   A directed version of a MultiGraph.

Empty graph-like objects are created with

.. nbplot::

   >>> G = nx.Graph()
   >>> G = nx.DiGraph()
   >>> G = nx.MultiGraph()
   >>> G = nx.MultiDiGraph()

All graph classes allow any :term:`hashable` object as a node.
Hashable objects include strings, tuples, integers, and more.
Arbitrary edge attributes such as weights and labels
can be associated with an edge.


The graph internal data structures are based on an
adjacency list representation and implemented using
Python :term:`dictionary` datastructures.
The graph adjacency structure is
implemented as a Python dictionary of
dictionaries; the outer dictionary is keyed by nodes to values that are
themselves dictionaries keyed by neighboring node to the
edge attributes associated with that edge.  This "dict-of-dicts" structure
allows fast addition, deletion, and lookup of nodes and neighbors in
large graphs.  The underlying datastructure is accessed directly
by methods (the programming interface "API") in the class definitions.
All functions, on the other hand, manipulate graph-like objects
solely via those API methods and not by acting directly on the datastructure.
This design allows for possible replacement of the 'dicts-of-dicts'-based
datastructure with an alternative datastructure that implements the
same methods.


Graphs
------

The first choice to be made when using NetworkX is what type of graph
object to use.  A graph (network) is a collection of nodes together
with a collection of edges that are pairs of nodes.  Attributes are
often associated with nodes and/or edges.  NetworkX graph objects come in
different flavors depending on two main properties of the network:

 - Directed: Are the edges **directed**?  Does the order of the edge
   pairs $(u, v)$ matter?  A directed graph is specified by the "Di"
   prefix in the class name, e.g. `DiGraph()`.  We make this distinction
   because many classical graph properties are defined differently for
   directed graphs.

 - Multi-edges: Are multiple edges allowed between each pair of nodes?
   As you might imagine, multiple edges requires a different data
   structure, though clever users could design edge data attributes to
   support this functionality.  We provide a standard data structure
   and interface for this type of graph using the prefix "Multi",
   e.g., `MultiGraph()`.

The basic graph classes are named:
:doc:`Graph </reference/classes/graph>`,
:doc:`DiGraph</reference/classes/digraph>`,
:doc:`MultiGraph </reference/classes/multigraph>`, and
:doc:`MultiDiGraph </reference/classes/multidigraph>`


Nodes and Edges
^^^^^^^^^^^^^^^

The next choice you have to make when specifying a graph is what kinds
of nodes and edges to use.

If the topology of the network is all you
care about then using integers or strings as the nodes makes sense and
you need not worry about edge data.  If you have a data structure
already in place to describe nodes you can simply use that structure
as your nodes provided it is :term:`hashable`.  If it is not hashable you can
use a unique identifier to represent the node and assign the data
as a :term:`node attribute`.

Edges often have data associated with them.  Arbitrary data
can be associated with edges as an :term:`edge attribute`.
If the data is numeric and the intent is to represent
a *weighted* graph then use the 'weight' keyword for the attribute.
Some of the graph algorithms, such as
Dijkstra's shortest path algorithm, use this attribute
name by default to get the weight for each edge.

Attributes can be assigned to an edge by using keyword/value
pairs when adding edges.  You can use any keyword
to name your attribute and can then query the edge
data using that attribute keyword.

Once you've decided how to encode the nodes and edges, and whether you have
an undirected/directed graph with or without multiedges you are ready to build
your network.

Graph Creation
--------------

NetworkX graph objects can be created in one of three ways:

- Graph generators---standard algorithms to create network topologies.
- Importing data from pre-existing (usually file) sources.
- Adding edges and nodes explicitly.

Explicit addition and removal of nodes/edges is the easiest to describe.
Each graph object supplies methods to manipulate the graph.  For example,

.. nbplot::

   >>> import networkx as nx
   >>> G = nx.Graph()
   >>> G.add_edge(1, 2)  # default edge data=1
   >>> G.add_edge(2, 3, weight=0.9)  # specify edge data

Edge attributes can be anything:

.. nbplot::

   >>> import math
   >>> G.add_edge('y', 'x', function=math.cos)
   >>> G.add_node(math.cos)  # any hashable can be a node

You can add many edges at one time:

.. nbplot::

   >>> elist = [(1, 2), (2, 3), (1, 4), (4, 2)]
   >>> G.add_edges_from(elist)
   >>> elist = [('a', 'b', 5.0), ('b', 'c', 3.0), ('a', 'c', 1.0), ('c', 'd', 7.3)]
   >>> G.add_weighted_edges_from(elist)

See the :doc:`/tutorial` for more examples.

Some basic graph operations such as union and intersection
are described in the :ref:`operators module <operators>` documentation.

Graph generators such as :func:`~generators.random_graphs.binomial_graph`
and :func:`~generators.random_graphs.erdos_renyi_graph` are
provided in the :ref:`graph generators <generators>` subpackage.

For importing network data from formats such as GML, GraphML, edge list text files
see the :ref:`reading and writing graphs <readwrite>` subpackage.


Graph Reporting
---------------

Class views provide basic reporting of nodes, neighbors, edges and degree.
These views provide iteration over the properties as well as membership
queries and data attribute lookup. The views refer to the graph data structure
so changes to the graph are reflected in the views. This is analogous to
dictionary views in Python 3. If you want to change the graph while iterating
you will need to use e.g. ``for e in list(G.edges):``. The views provide
set-like operations, e.g. union and intersection, as well as dict-like
lookup and iteration of the data attributes using ``G.edges[u, v]['color']``
and ``for e, datadict in G.edges.items():``. Methods ``G.edges.items()`` and
``G.edges.values()`` are familiar from python dicts. In addition ``G.edges.data()``
provides specific attribute iteration e.g. ``for e, e_color in G.edges.data('color'):``.

The basic graph relationship of an edge can be obtained in two ways.
One can look for neighbors of a node or one can look for edges.
We jokingly refer to people who focus on nodes/neighbors as node-centric
and people who focus on edges as edge-centric.  The designers of NetworkX
tend to be node-centric and view edges as a relationship between nodes.
You can see this by our choice of lookup notation like ``G[u]`` providing neighbors
(adjacency) while edge lookup is ``G.edges[u, v]``.
Most data structures for sparse graphs are essentially adjacency lists and so
fit this perspective. In the end, of course, it doesn't really matter which way
you examine the graph. ``G.edges`` removes duplicate representations of undirected
edges while neighbor reporting across all nodes will naturally report both directions.

Any properties that are more complicated than edges, neighbors and degree are
provided by functions.  For example ``nx.triangles(G, n)`` gives the number of triangles
which include node n as a vertex.  These functions are grouped in the code and
documentation under the term :ref:`algorithms<algorithms>`.


Algorithms
----------

A number of graph algorithms are provided with NetworkX.
These include shortest path, and breadth first search
(see :ref:`traversal<traversal>`),
clustering and isomorphism algorithms and others.  There are
many that we have not developed yet too.  If you implement a
graph algorithm that might be useful for others please let
us know through the
`NetworkX Google group <http://groups.google.com/group/networkx-discuss>`_
or the Github `Developer Zone <https://github.com/networkx/networkx>`_.

As an example here is code to use Dijkstra's algorithm to
find the shortest weighted path:

.. nbplot::

   >>> G = nx.Graph()
   >>> e = [('a', 'b', 0.3), ('b', 'c', 0.9), ('a', 'c', 0.5), ('c', 'd', 1.2)]
   >>> G.add_weighted_edges_from(e)
   >>> print(nx.dijkstra_path(G, 'a', 'd'))
   ['a', 'c', 'd']

Drawing
-------

While NetworkX is not designed as a network drawing tool, we provide
a simple interface to drawing packages and some simple layout algorithms.
We interface to the excellent Graphviz layout tools like dot and neato
with the (suggested) pygraphviz package or the pydot interface.
Drawing can be done using external programs or the Matplotlib Python
package.  Interactive GUI interfaces are possible, though not provided.
The drawing tools are provided in the module :ref:`drawing <drawing>`.

The basic drawing functions essentially place the nodes on a scatterplot
using the positions you provide via a dictionary or the positions are
computed with a layout function. The edges are lines between those dots.

.. nbplot::

   >>> import matplotlib.pyplot as plt
   >>> G = nx.cubical_graph()
   >>> subax1 = plt.subplot(121)
   >>> nx.draw(G)   # default spring_layout
   >>> subax2 = plt.subplot(122)
   >>> nx.draw(G, pos=nx.circular_layout(G), node_color='r', edge_color='b')

See the :doc:`examples </auto_examples/index>` for more ideas.

Data Structure
--------------

NetworkX uses a "dictionary of dictionaries of dictionaries" as the
basic network data structure.  This allows fast lookup with reasonable
storage for large sparse networks.  The keys are nodes so ``G[u]`` returns
an adjacency dictionary keyed by neighbor to the edge attribute
dictionary. A view of the adjacency data structure is provided
by the dict-like object ``G.adj`` as e.g. ``for node, nbrsdict in G.adj.items():``.
The expression ``G[u][v]`` returns the edge attribute dictionary itself.
A dictionary of lists would have also been possible, but not allow
fast edge detection nor convenient storage of edge data.

Advantages of dict-of-dicts-of-dicts data structure:

 - Find edges and remove edges with two dictionary look-ups.
 - Prefer to "lists" because of fast lookup with sparse storage.
 - Prefer to "sets" since data can be attached to edge.
 - ``G[u][v]`` returns the edge attribute dictionary.
 - ``n in G`` tests if node ``n`` is in graph ``G``.
 - ``for n in G:`` iterates through the graph.
 - ``for nbr in G[n]:`` iterates through neighbors.

As an example, here is a representation of an undirected graph with the
edges $(A, B)$ and $(B, C)$.

.. nbplot::

   >>> G = nx.Graph()
   >>> G.add_edge('A', 'B')
   >>> G.add_edge('B', 'C')
   >>> print(G.adj)
   {'A': {'B': {}}, 'B': {'A': {}, 'C': {}}, 'C': {'B': {}}}

The data structure gets morphed slightly for each base graph class.
For DiGraph two dict-of-dicts-of-dicts structures are provided, one
for successors (``G.succ``) and one for predecessors (``G.pred``).
For MultiGraph/MultiDiGraph we use a dict-of-dicts-of-dicts-of-dicts [#turtles]_
where the third dictionary is keyed by an edge key identifier to the fourth
dictionary which contains the edge attributes for that edge between
the two nodes.

Graphs provide two interfaces to the edge data attributes: adjacency
and edges. So ``G[u][v]['width']`` is the same as ``G.edges[u, v]['width']``.

.. nbplot::

   >>> G = nx.Graph()
   >>> G.add_edge(1, 2, color='red', weight=0.84, size=300)
   >>> print(G[1][2]['size'])
   300
   >>> print(G.edges[1, 2]['color'])
   red

.. code-links::

.. rubric:: Footnotes

.. [#turtles] "It's dictionaries all the way down."
SparseGraph6
============

Functions for reading and writing graphs in the *graph6* or *sparse6* file
formats.

According to the author of these formats,

    *graph6* and *sparse6* are formats for storing undirected graphs in a
    compact manner, using only printable ASCII characters. Files in these
    formats have text type and contain one line per graph.

    *graph6* is suitable for small graphs, or large dense graphs. *sparse6* is
    more space-efficient for large sparse graphs.

    -- `graph6 and sparse6 homepage`_

.. _graph6 and sparse6 homepage: http://users.cecs.anu.edu.au/~bdm/data/formats.html

Graph6
------
.. automodule:: networkx.readwrite.graph6
.. autosummary::
   :toctree: generated/

   from_graph6_bytes
   read_graph6
   to_graph6_bytes
   write_graph6

Sparse6
-------
.. automodule:: networkx.readwrite.sparse6
.. autosummary::
   :toctree: generated/

   from_sparse6_bytes
   read_sparse6
   to_sparse6_bytes
   write_sparse6

Adjacency List
==============

.. automodule:: networkx.readwrite.adjlist
.. autosummary::
   :toctree: generated/

   read_adjlist
   write_adjlist
   parse_adjlist
   generate_adjlist
GraphML
=======
.. automodule:: networkx.readwrite.graphml
.. autosummary::
   :toctree: generated/

   read_graphml
   write_graphml
   generate_graphml
   parse_graphml



Multiline Adjacency List
========================

.. automodule:: networkx.readwrite.multiline_adjlist
.. autosummary::
   :toctree: generated/

   read_multiline_adjlist
   write_multiline_adjlist
   parse_multiline_adjlist
   generate_multiline_adjlist
GIS Shapefile
=============
.. automodule:: networkx.readwrite.nx_shp
.. autosummary::
   :toctree: generated/

   read_shp
   write_shp


Pajek
=====
.. automodule:: networkx.readwrite.pajek
.. autosummary::
   :toctree: generated/

   read_pajek
   write_pajek
   parse_pajek
   generate_pajek
JSON
====
.. automodule:: networkx.readwrite.json_graph
.. autosummary::
   :toctree: generated/

   node_link_data
   node_link_graph
   adjacency_data
   adjacency_graph
   cytoscape_data
   cytoscape_graph
   tree_data
   tree_graph
   jit_data
   jit_graph
*************
Matrix Market
*************

The `Matrix Market`_ exchange format is a text-based file format described by
NIST.
Matrix Market supports both a **coordinate format** for sparse matrices and
an **array format** for dense matrices.
The :mod:`scipy.io` module provides the `scipy.io.mmread` and `scipy.io.mmwrite`
functions to read and write data in Matrix Market format, respectively.
These functions work with either `numpy.ndarray` or `scipy.sparse.coo_matrix`
objects depending on whether the data is in **array** or **coordinate** format.
These functions can be combined with those of NetworkX's `~networkx.convert_matrix`
module to read and write Graphs in Matrix Market format.

.. _Matrix Market: https://math.nist.gov/MatrixMarket/formats.html

Examples
========

Reading and writing graphs using Matrix Market's **array format** for dense
matrices::

    >>> import scipy as sp
    >>> import scipy.io  # for mmread() and mmwrite()
    >>> import io  # Use BytesIO as a stand-in for a Python file object
    >>> fh = io.BytesIO()

    >>> G = nx.complete_graph(5)
    >>> a = nx.to_numpy_array(G)
    >>> print(a)
    [[0. 1. 1. 1. 1.]
     [1. 0. 1. 1. 1.]
     [1. 1. 0. 1. 1.]
     [1. 1. 1. 0. 1.]
     [1. 1. 1. 1. 0.]]

    >>> # Write to file in Matrix Market array format
    >>> sp.io.mmwrite(fh, a)
    >>> print(fh.getvalue().decode('utf-8'))  # file contents
    %%MatrixMarket matrix array real symmetric
    %
    5 5
    0.0000000000000000e+00
    1.0000000000000000e+00
    1.0000000000000000e+00
    1.0000000000000000e+00
    1.0000000000000000e+00
    0.0000000000000000e+00
    1.0000000000000000e+00
    1.0000000000000000e+00
    1.0000000000000000e+00
    0.0000000000000000e+00
    1.0000000000000000e+00
    1.0000000000000000e+00
    0.0000000000000000e+00
    1.0000000000000000e+00
    0.0000000000000000e+00

    >>> # Read from file
    >>> fh.seek(0)
    >>> H = nx.from_numpy_array(sp.io.mmread(fh))
    >>> H.edges() == G.edges()
    True

Reading and writing graphs using Matrix Market's **coordinate format** for
sparse matrices::

    >>> import scipy as sp
    >>> import scipy.io  # for mmread() and mmwrite()
    >>> import io  # Use BytesIO as a stand-in for a Python file object
    >>> fh = io.BytesIO()

    >>> G = nx.path_graph(5)
    >>> m = nx.to_scipy_sparse_array(G)
    >>> print(m)
      (0, 1)        1
      (1, 0)        1
      (1, 2)        1
      (2, 1)        1
      (2, 3)        1
      (3, 2)        1
      (3, 4)        1
      (4, 3)        1

    >>> sp.io.mmwrite(fh, m)
    >>> print(fh.getvalue().decode('utf-8'))  # file contents
    %%MatrixMarket matrix coordinate integer symmetric
    %
    5 5 4
    2 1 1
    3 2 1
    4 3 1
    5 4 1

    >>> # Read from file
    >>> fh.seek(0)
    >>> H = nx.from_scipy_sparse_matrix(sp.io.mmread(fh))
    >>> H.edges() == G.edges()
    True
.. _readwrite:

**************************
Reading and writing graphs
**************************

.. toctree::
   :maxdepth: 2

   adjlist
   multiline_adjlist
   edgelist
   gexf
   gml
   gpickle
   graphml
   json_graph
   leda
   sparsegraph6
   pajek
   nx_shp
   matrix_market
LEDA
====
.. automodule:: networkx.readwrite.leda
.. autosummary::
   :toctree: generated/

   read_leda
   parse_leda



GML
===
.. automodule:: networkx.readwrite.gml
.. autosummary::
   :toctree: generated/

   read_gml
   write_gml
   parse_gml
   generate_gml
   literal_destringizer
   literal_stringizer
GEXF
====
.. automodule:: networkx.readwrite.gexf
.. autosummary::
   :toctree: generated/

   read_gexf
   write_gexf
   generate_gexf
   relabel_gexf_graph

Pickle
======
.. automodule:: networkx.readwrite.gpickle
.. autosummary::
   :toctree: generated/

   read_gpickle
   write_gpickle

Edge List
=========
.. automodule:: networkx.readwrite.edgelist
.. autosummary::
   :toctree: generated/

   read_edgelist
   write_edgelist
   read_weighted_edgelist
   write_weighted_edgelist
   generate_edgelist
   parse_edgelist
Shortest Paths
==============
.. automodule:: networkx.algorithms.shortest_paths.generic
.. autosummary::
   :toctree: generated/

   shortest_path
   all_shortest_paths
   shortest_path_length
   average_shortest_path_length
   has_path


Advanced Interface
------------------

.. automodule:: networkx.algorithms.shortest_paths.unweighted
.. autosummary::
   :toctree: generated/

   single_source_shortest_path
   single_source_shortest_path_length
   single_target_shortest_path
   single_target_shortest_path_length
   bidirectional_shortest_path
   all_pairs_shortest_path
   all_pairs_shortest_path_length
   predecessor

.. automodule:: networkx.algorithms.shortest_paths.weighted
.. autosummary::
   :toctree: generated/

   dijkstra_predecessor_and_distance
   dijkstra_path
   dijkstra_path_length
   single_source_dijkstra
   single_source_dijkstra_path
   single_source_dijkstra_path_length
   multi_source_dijkstra
   multi_source_dijkstra_path
   multi_source_dijkstra_path_length
   all_pairs_dijkstra
   all_pairs_dijkstra_path
   all_pairs_dijkstra_path_length
   bidirectional_dijkstra

   bellman_ford_path
   bellman_ford_path_length
   single_source_bellman_ford
   single_source_bellman_ford_path
   single_source_bellman_ford_path_length
   all_pairs_bellman_ford_path
   all_pairs_bellman_ford_path_length
   bellman_ford_predecessor_and_distance

   negative_edge_cycle
   find_negative_cycle
   goldberg_radzik
   johnson


Dense Graphs
------------

.. automodule:: networkx.algorithms.shortest_paths.dense
.. autosummary::
   :toctree: generated/

   floyd_warshall
   floyd_warshall_predecessor_and_distance
   floyd_warshall_numpy
   reconstruct_path


A* Algorithm
------------

.. automodule:: networkx.algorithms.shortest_paths.astar
.. autosummary::
   :toctree: generated/

   astar_path
   astar_path_length

**********
Asteroidal
**********

.. automodule:: networkx.algorithms.asteroidal
.. autosummary::
   :toctree: generated/

   is_at_free
   find_asteroidal_triple
**********
Clustering
**********

.. automodule:: networkx.algorithms.cluster
.. autosummary::
   :toctree: generated/

   triangles
   transitivity
   clustering
   average_clustering
   square_clustering
   generalized_degree
*********
Dominance
*********

.. automodule:: networkx.algorithms.dominance
.. autosummary::
   :toctree: generated/

   immediate_dominators
   dominance_frontiers
******
Clique
******

.. automodule:: networkx.algorithms.clique
.. autosummary::
   :toctree: generated/

   enumerate_all_cliques
   find_cliques
   make_max_clique_graph
   make_clique_bipartite
   graph_clique_number
   graph_number_of_cliques
   node_clique_number
   number_of_cliques
   cliques_containing_node
   max_weight_clique
***********************
Directed Acyclic Graphs
***********************

.. automodule:: networkx.algorithms.dag
.. autosummary::
   :toctree: generated/

   ancestors
   descendants
   topological_sort
   topological_generations
   all_topological_sorts
   lexicographical_topological_sort
   is_directed_acyclic_graph
   is_aperiodic
   transitive_closure
   transitive_closure_dag
   transitive_reduction
   antichains
   dag_longest_path
   dag_longest_path_length
   dag_to_branching
***************
Dominating Sets
***************

.. automodule:: networkx.algorithms.dominating
.. autosummary::
   :toctree: generated/

   dominating_set
   is_dominating_set
*****************************
Approximations and Heuristics
*****************************

.. automodule:: networkx.algorithms.approximation


Connectivity
------------
.. automodule:: networkx.algorithms.approximation.connectivity
.. autosummary::
   :toctree: generated/

   all_pairs_node_connectivity
   local_node_connectivity
   node_connectivity


K-components
------------
.. automodule:: networkx.algorithms.approximation.kcomponents
.. autosummary::
   :toctree: generated/

   k_components


Clique
------
.. automodule:: networkx.algorithms.approximation.clique
.. autosummary::
   :toctree: generated/

   maximum_independent_set
   max_clique
   clique_removal
   large_clique_size


Clustering
----------
.. automodule:: networkx.algorithms.approximation.clustering_coefficient
.. autosummary::
   :toctree: generated/

   average_clustering


Distance Measures
-----------------
.. automodule:: networkx.algorithms.approximation.distance_measures
.. autosummary::
   :toctree: generated/

   diameter


Dominating Set
---------------
.. automodule:: networkx.algorithms.approximation.dominating_set
.. autosummary::
   :toctree: generated/

   min_weighted_dominating_set
   min_edge_dominating_set

Matching
--------
.. automodule:: networkx.algorithms.approximation.matching
.. autosummary::
   :toctree: generated/

   min_maximal_matching

Ramsey
------
.. automodule:: networkx.algorithms.approximation.ramsey
.. autosummary::
   :toctree: generated/

   ramsey_R2


Steiner Tree
------------
.. automodule:: networkx.algorithms.approximation.steinertree
.. autosummary::
   :toctree: generated/

   metric_closure
   steiner_tree


Traveling Salesman
------------------
.. automodule:: networkx.algorithms.approximation.traveling_salesman
.. autosummary::
   :toctree: generated/

   christofides
   traveling_salesman_problem
   greedy_tsp
   simulated_annealing_tsp
   threshold_accepting_tsp
   asadpour_atsp


Treewidth
---------
.. automodule:: networkx.algorithms.approximation.treewidth
.. autosummary::
   :toctree: generated/

   treewidth_min_degree
   treewidth_min_fill_in


Vertex Cover
------------
.. automodule:: networkx.algorithms.approximation.vertex_cover
.. autosummary::
   :toctree: generated/

   min_weighted_vertex_cover


Max Cut
-------
.. automodule:: networkx.algorithms.approximation.maxcut
.. autosummary::
   :toctree: generated/

   randomized_partitioning
   one_exchange
============
D-Separation
============

.. automodule:: networkx.algorithms.d_separation
.. autosummary::
   :toctree: generated/

   d_separated
*************
Voronoi cells
*************

.. automodule:: networkx.algorithms.voronoi
.. autosummary::
   :toctree: generated/

   voronoi_cells
*********
Bipartite
*********

.. automodule:: networkx.algorithms.bipartite


Basic functions
---------------
.. automodule:: networkx.algorithms.bipartite.basic
.. autosummary::
   :toctree: generated/

   is_bipartite
   is_bipartite_node_set
   sets
   color
   density
   degrees

Edgelist
--------
.. automodule:: networkx.algorithms.bipartite.edgelist
.. autosummary::
   :toctree: generated/

   generate_edgelist
   write_edgelist
   parse_edgelist
   read_edgelist


Matching
--------
.. automodule:: networkx.algorithms.bipartite.matching
.. autosummary::
   :toctree: generated/

   eppstein_matching
   hopcroft_karp_matching
   to_vertex_cover
   maximum_matching
   minimum_weight_full_matching


Matrix
------
.. automodule:: networkx.algorithms.bipartite.matrix
.. autosummary::
   :toctree: generated/

   biadjacency_matrix
   from_biadjacency_matrix


Projections
-----------
.. automodule:: networkx.algorithms.bipartite.projection
.. autosummary::
   :toctree: generated/

   projected_graph
   weighted_projected_graph
   collaboration_weighted_projected_graph
   overlap_weighted_projected_graph
   generic_weighted_projected_graph


Spectral
--------
.. automodule:: networkx.algorithms.bipartite.spectral
.. autosummary::
   :toctree: generated/

   spectral_bipartivity


Clustering
----------
.. automodule:: networkx.algorithms.bipartite.cluster
.. autosummary::
   :toctree: generated/

   clustering
   average_clustering
   latapy_clustering
   robins_alexander_clustering


Redundancy
----------
.. automodule:: networkx.algorithms.bipartite.redundancy
.. autosummary::
   :toctree: generated/

   node_redundancy


Centrality
----------
.. automodule:: networkx.algorithms.bipartite.centrality
.. autosummary::
   :toctree: generated/

   closeness_centrality
   degree_centrality
   betweenness_centrality


Generators
----------
.. automodule:: networkx.algorithms.bipartite.generators
.. autosummary::
   :toctree: generated/

   complete_bipartite_graph
   configuration_model
   havel_hakimi_graph
   reverse_havel_hakimi_graph
   alternating_havel_hakimi_graph
   preferential_attachment_graph
   random_graph
   gnmk_random_graph


Covering
--------
.. automodule:: networkx.algorithms.bipartite.covering
.. autosummary::
   :toctree: generated/

   min_edge_cover



******
Minors
******

.. automodule:: networkx.algorithms.minors
.. autosummary::
   :toctree: generated/

   contracted_edge
   contracted_nodes
   identified_nodes
   equivalence_classes
   quotient_graph
*****
Moral
*****

.. automodule:: networkx.algorithms.moral
.. autosummary::
   :toctree: generated/

   moral_graph
************
Wiener index
************

.. automodule:: networkx.algorithms.wiener
.. autosummary::
   :toctree: generated/

   wiener_index
********
s metric
********

.. automodule:: networkx.algorithms.smetric
.. autosummary::
   :toctree: generated/

   s_metric
************
Connectivity
************

.. automodule:: networkx.algorithms.connectivity

Edge-augmentation
-----------------
.. automodule:: networkx.algorithms.connectivity.edge_augmentation
.. autosummary::
   :toctree: generated/

   k_edge_augmentation
   is_k_edge_connected
   is_locally_k_edge_connected

K-edge-components
-----------------
.. automodule:: networkx.algorithms.connectivity.edge_kcomponents
.. autosummary::
   :toctree: generated/

   k_edge_components
   k_edge_subgraphs
   bridge_components
   EdgeComponentAuxGraph

K-node-components
-----------------
.. automodule:: networkx.algorithms.connectivity.kcomponents
.. autosummary::
   :toctree: generated/

   k_components

K-node-cutsets
--------------
.. automodule:: networkx.algorithms.connectivity.kcutsets
.. autosummary::
   :toctree: generated/

   all_node_cuts

Flow-based disjoint paths
-------------------------
.. automodule:: networkx.algorithms.connectivity.disjoint_paths
.. autosummary::
   :toctree: generated/

   edge_disjoint_paths
   node_disjoint_paths

Flow-based Connectivity
-----------------------
.. automodule:: networkx.algorithms.connectivity.connectivity
.. autosummary::
   :toctree: generated/

   average_node_connectivity
   all_pairs_node_connectivity
   edge_connectivity
   local_edge_connectivity
   local_node_connectivity
   node_connectivity

Flow-based Minimum Cuts
-----------------------
.. automodule:: networkx.algorithms.connectivity.cuts
.. autosummary::
   :toctree: generated/

   minimum_edge_cut
   minimum_node_cut
   minimum_st_edge_cut
   minimum_st_node_cut

Stoer-Wagner minimum cut
------------------------
.. automodule:: networkx.algorithms.connectivity.stoerwagner
.. autosummary::
   :toctree: generated/

   stoer_wagner

Utils for flow-based connectivity
---------------------------------
.. automodule:: networkx.algorithms.connectivity.utils
.. autosummary::
   :toctree: generated/

    build_auxiliary_edge_connectivity
    build_auxiliary_node_connectivity
*********
Rich Club
*********

.. automodule:: networkx.algorithms.richclub
.. autosummary::
   :toctree: generated/

   rich_club_coefficient
Node Classification
===================
.. automodule:: networkx.algorithms.node_classification
.. currentmodule:: networkx

Harmonic Function
-----------------
.. automodule:: networkx.algorithms.node_classification.hmn
.. autosummary::
   :toctree: generated/

   harmonic_function


Local and Global Consistency
----------------------------
.. automodule:: networkx.algorithms.node_classification.lgc
.. autosummary::
   :toctree: generated/

   local_and_global_consistency
*****
Flows
*****

.. automodule:: networkx.algorithms.flow


Maximum Flow
------------
.. autosummary::
   :toctree: generated/

   maximum_flow
   maximum_flow_value
   minimum_cut
   minimum_cut_value


Edmonds-Karp
------------
.. autosummary::
   :toctree: generated/

   edmonds_karp


Shortest Augmenting Path
------------------------
.. autosummary::
   :toctree: generated/

   shortest_augmenting_path


Preflow-Push
------------
.. autosummary::
   :toctree: generated/

   preflow_push


Dinitz
------
.. autosummary::
   :toctree: generated/

   dinitz


Boykov-Kolmogorov
-----------------
.. autosummary::
   :toctree: generated/

   boykov_kolmogorov


Gomory-Hu Tree
--------------
.. autosummary::
   :toctree: generated/

   gomory_hu_tree


Utils
-----
.. autosummary::
   :toctree: generated/

   build_residual_network


Network Simplex
---------------
.. autosummary::
   :toctree: generated/

    network_simplex
    min_cost_flow_cost
    min_cost_flow
    cost_of_flow
    max_flow_min_cost


Capacity Scaling Minimum Cost Flow
----------------------------------
.. autosummary::
   :toctree: generated/

    capacity_scaling
**********
Centrality
**********

.. automodule:: networkx.algorithms.centrality

Degree
------
.. autosummary::
   :toctree: generated/

   degree_centrality
   in_degree_centrality
   out_degree_centrality

Eigenvector
-----------
.. autosummary::
   :toctree: generated/

   eigenvector_centrality
   eigenvector_centrality_numpy
   katz_centrality
   katz_centrality_numpy

Closeness
---------
.. autosummary::
   :toctree: generated/

   closeness_centrality
   incremental_closeness_centrality

Current Flow Closeness
----------------------
.. autosummary::
   :toctree: generated/

   current_flow_closeness_centrality
   information_centrality

(Shortest Path) Betweenness
---------------------------
.. autosummary::
   :toctree: generated/

   betweenness_centrality
   betweenness_centrality_source
   betweenness_centrality_subset
   edge_betweenness_centrality
   edge_betweenness_centrality_subset


Current Flow Betweenness
------------------------
.. autosummary::
   :toctree: generated/

   current_flow_betweenness_centrality
   edge_current_flow_betweenness_centrality
   approximate_current_flow_betweenness_centrality
   current_flow_betweenness_centrality_subset
   edge_current_flow_betweenness_centrality_subset

Communicability Betweenness
---------------------------
.. autosummary::
   :toctree: generated/

   communicability_betweenness_centrality

Group Centrality
----------------
.. autosummary::
   :toctree: generated/

   group_betweenness_centrality
   group_closeness_centrality
   group_degree_centrality
   group_in_degree_centrality
   group_out_degree_centrality
   prominent_group

Load
----
.. autosummary::
   :toctree: generated/

   load_centrality
   edge_load_centrality

Subgraph
--------
.. autosummary::
   :toctree: generated/

   subgraph_centrality
   subgraph_centrality_exp
   estrada_index

Harmonic Centrality
-------------------
.. autosummary::
   :toctree: generated/

   harmonic_centrality

Dispersion
----------
.. autosummary::
   :toctree: generated/

   dispersion

Reaching
--------
.. autosummary::
   :toctree: generated/

   local_reaching_centrality
   global_reaching_centrality

Percolation
-----------
.. autosummary::
   :toctree: generated/

   percolation_centrality

Second Order Centrality
-----------------------
.. autosummary::
   :toctree: generated/

   second_order_centrality

Trophic
-------
.. autosummary::
   :toctree: generated/

   trophic_levels
   trophic_differences
   trophic_incoherence_parameter

VoteRank
--------
.. autosummary::
   :toctree: generated/

   voterank
*************
Link Analysis
*************

PageRank
--------

.. automodule:: networkx.algorithms.link_analysis.pagerank_alg
.. autosummary::
   :toctree: generated/

   pagerank
   pagerank_numpy
   pagerank_scipy
   google_matrix

Hits
----

.. automodule:: networkx.algorithms.link_analysis.hits_alg
.. autosummary::
   :toctree: generated/

   hits
   hits_numpy
   hits_scipy
   hub_matrix
   authority_matrix

********
Coloring
********

.. automodule:: networkx.algorithms.coloring
.. autosummary::
   :toctree: generated/

   greedy_color
   equitable_color

Some node ordering strategies are provided for use with :func:`greedy_color`.

.. autosummary::
   :toctree: generated/

    strategy_connected_sequential
    strategy_connected_sequential_dfs
    strategy_connected_sequential_bfs
    strategy_independent_set
    strategy_largest_first
    strategy_random_sequential
    strategy_saturation_largest_first
    strategy_smallest_last
***********************
Distance-Regular Graphs
***********************

.. automodule:: networkx.algorithms.distance_regular
.. autosummary::
   :toctree: generated/

   is_distance_regular
   is_strongly_regular
   intersection_array
   global_parameters
.. _tree:

Tree
====

.. toctree::
   :maxdepth: 2

Recognition
-----------
.. automodule:: networkx.algorithms.tree.recognition
.. autosummary::
   :toctree: generated/

   is_tree
   is_forest
   is_arborescence
   is_branching

Branchings and Spanning Arborescences
-------------------------------------
.. automodule:: networkx.algorithms.tree.branchings
.. autosummary::
   :toctree: generated/

   branching_weight
   greedy_branching
   maximum_branching
   minimum_branching
   maximum_spanning_arborescence
   minimum_spanning_arborescence
   ArborescenceIterator
   Edmonds

Encoding and decoding
---------------------
.. automodule:: networkx.algorithms.tree.coding
.. autosummary::
   :toctree: generated/

   from_nested_tuple
   to_nested_tuple
   from_prufer_sequence
   to_prufer_sequence

Operations
----------
.. automodule:: networkx.algorithms.tree.operations
.. autosummary::
   :toctree: generated/

   join

Spanning Trees
--------------
.. automodule:: networkx.algorithms.tree.mst
.. autosummary::
   :toctree: generated/

   minimum_spanning_tree
   maximum_spanning_tree
   minimum_spanning_edges
   maximum_spanning_edges
   SpanningTreeIterator

Decomposition
-------------
.. automodule:: networkx.algorithms.tree.decomposition
.. autosummary::
   :toctree: generated/

   junction_tree

Exceptions
----------
.. automodule:: networkx.algorithms.tree.coding
   :noindex:
.. autosummary::
   :toctree: generated/

   NotATree
********
Isolates
********

.. automodule:: networkx.algorithms.isolate
.. autosummary::
   :toctree: generated/

   is_isolate
   isolates
   number_of_isolates
.. _isomorphism:

***********
Isomorphism
***********

.. toctree::
   :maxdepth: 2

.. automodule:: networkx.algorithms.isomorphism
.. autosummary::
   :toctree: generated/

   is_isomorphic
   could_be_isomorphic
   fast_could_be_isomorphic
   faster_could_be_isomorphic


Tree Isomorphism
-----------------
.. automodule:: networkx.algorithms.isomorphism.tree_isomorphism
.. autosummary::
   :toctree: generated/

    rooted_tree_isomorphism
    tree_isomorphism


Advanced Interfaces
-------------------
.. toctree::
   :maxdepth: 2

   isomorphism.vf2
   isomorphism.ismags

Bridges
=======

.. automodule:: networkx.algorithms.bridges
.. autosummary::
   :toctree: generated/

   bridges
   has_bridges
   local_bridges
.. _vf2:

*************
VF2 Algorithm
*************

.. automodule:: networkx.algorithms.isomorphism.isomorphvf2

Graph Matcher
-------------
.. currentmodule:: networkx.algorithms.isomorphism

.. autosummary::
   :toctree: generated/

    GraphMatcher.__init__
    GraphMatcher.initialize
    GraphMatcher.is_isomorphic
    GraphMatcher.subgraph_is_isomorphic
    GraphMatcher.isomorphisms_iter
    GraphMatcher.subgraph_isomorphisms_iter
    GraphMatcher.candidate_pairs_iter
    GraphMatcher.match
    GraphMatcher.semantic_feasibility
    GraphMatcher.syntactic_feasibility


DiGraph Matcher
---------------
.. currentmodule:: networkx.algorithms.isomorphism

.. autosummary::
   :toctree: generated/

    DiGraphMatcher.__init__
    DiGraphMatcher.initialize
    DiGraphMatcher.is_isomorphic
    DiGraphMatcher.subgraph_is_isomorphic
    DiGraphMatcher.isomorphisms_iter
    DiGraphMatcher.subgraph_isomorphisms_iter
    DiGraphMatcher.candidate_pairs_iter
    DiGraphMatcher.match
    DiGraphMatcher.semantic_feasibility
    DiGraphMatcher.syntactic_feasibility


Match helpers
-------------
.. currentmodule:: networkx.algorithms.isomorphism

.. autosummary::
   :toctree: generated/

   categorical_node_match
   categorical_edge_match
   categorical_multiedge_match
   numerical_node_match
   numerical_edge_match
   numerical_multiedge_match
   generic_node_match
   generic_edge_match
   generic_multiedge_match

**********************
Lowest Common Ancestor
**********************

.. automodule:: networkx.algorithms.lowest_common_ancestors
.. autosummary::
   :toctree: generated/

   all_pairs_lowest_common_ancestor
   tree_all_pairs_lowest_common_ancestor
   lowest_common_ancestor
***********
Small-world
***********

.. automodule:: networkx.algorithms.smallworld
.. autosummary::
   :toctree: generated/

   random_reference
   lattice_reference
   sigma
   omega
*************
Summarization
*************

.. automodule:: networkx.algorithms.summarization
.. autosummary::
   :toctree: generated/

   dedensify
   snap_aggregation
Chains
======

.. automodule:: networkx.algorithms.chains
.. autosummary::
   :toctree: generated/

   chain_decomposition
***********
Sparsifiers
***********

.. automodule:: networkx.algorithms.sparsifiers
.. autosummary::
   :toctree: generated/

   spanner
********
Matching
********

.. automodule:: networkx.algorithms.matching
.. autosummary::
   :toctree: generated/

   is_matching
   is_maximal_matching
   is_perfect_matching
   maximal_matching
   max_weight_matching
   min_weight_matching
******
Hybrid
******

.. automodule:: networkx.algorithms.hybrid
.. autosummary::
   :toctree: generated/

   kl_connected_subgraph
   is_kl_connected
***************
Link Prediction
***************

.. automodule:: networkx.algorithms.link_prediction
.. autosummary::
   :toctree: generated/

   resource_allocation_index
   jaccard_coefficient
   adamic_adar_index
   preferential_attachment
   cn_soundarajan_hopcroft
   ra_index_soundarajan_hopcroft
   within_inter_cluster
   common_neighbor_centrality
*******
Regular
*******

.. automodule:: networkx.algorithms.regular
.. autosummary::
   :toctree: generated/

   is_regular
   is_k_regular
   k_factor
***********
Communities
***********

.. automodule:: networkx.algorithms.community
.. currentmodule:: networkx


Bipartitions
------------
.. automodule:: networkx.algorithms.community.kernighan_lin
.. autosummary::
   :toctree: generated/

   kernighan_lin_bisection

K-Clique
--------
.. automodule:: networkx.algorithms.community.kclique
.. autosummary::
   :toctree: generated/

   k_clique_communities

Modularity-based communities
----------------------------
.. automodule:: networkx.algorithms.community.modularity_max
.. autosummary::
   :toctree: generated/

   greedy_modularity_communities
   naive_greedy_modularity_communities

Tree partitioning
-----------------
.. automodule:: networkx.algorithms.community.lukes
.. autosummary::
   :toctree: generated/

   lukes_partitioning

Label propagation
-----------------
.. automodule:: networkx.algorithms.community.label_propagation
.. autosummary::
   :toctree: generated/

   asyn_lpa_communities
   label_propagation_communities

Louvain Community Detection
---------------------------
.. automodule:: networkx.algorithms.community.louvain
.. autosummary::
    :toctree: generated/

    louvain_communities
    louvain_partitions

Fluid Communities
-----------------
.. automodule:: networkx.algorithms.community.asyn_fluid
.. autosummary::
   :toctree: generated/

   asyn_fluidc

Measuring partitions
--------------------
.. automodule:: networkx.algorithms.community.quality
.. autosummary::
   :toctree: generated/

   coverage
   modularity
   partition_quality
   performance

Partitions via centrality measures
----------------------------------
.. automodule:: networkx.algorithms.community.centrality
.. autosummary::
   :toctree: generated/

   girvan_newman

Validating partitions
---------------------
.. automodule:: networkx.algorithms.community.community_utils
.. autosummary::
   :toctree: generated/

   is_partition
*************************
Graphical degree sequence
*************************

.. automodule:: networkx.algorithms.graphical
.. autosummary::
   :toctree: generated/

   is_graphical
   is_digraphical
   is_multigraphical
   is_pseudographical
   is_valid_degree_sequence_havel_hakimi
   is_valid_degree_sequence_erdos_gallai


********
Eulerian
********

.. automodule:: networkx.algorithms.euler
.. autosummary::
   :toctree: generated/

   is_eulerian
   eulerian_circuit
   eulerize
   is_semieulerian
   has_eulerian_path
   eulerian_path
*************
Assortativity
*************

.. automodule:: networkx.algorithms.assortativity
.. autosummary::
   :toctree: generated/

Assortativity
-------------
.. autosummary::
   :toctree: generated/

   degree_assortativity_coefficient
   attribute_assortativity_coefficient
   numeric_assortativity_coefficient
   degree_pearson_correlation_coefficient

Average neighbor degree
-----------------------
.. autosummary::
   :toctree: generated/

   average_neighbor_degree


Average degree connectivity
---------------------------
.. autosummary::
   :toctree: generated/

   average_degree_connectivity
   k_nearest_neighbors


Mixing
------
.. autosummary::
   :toctree: generated/

   attribute_mixing_matrix
   degree_mixing_matrix
   numeric_mixing_matrix
   attribute_mixing_dict
   degree_mixing_dict
   mixing_dict

Pairs
-----
.. autosummary::
   :toctree: generated/

   node_attribute_xy
   node_degree_xy

*****
Cores
*****

.. automodule:: networkx.algorithms.core
.. autosummary::
   :toctree: generated/

   core_number
   k_core
   k_shell
   k_crust
   k_corona
   k_truss
   onion_layers
*********
Hierarchy
*********

.. automodule:: networkx.algorithms.hierarchy
.. autosummary::
   :toctree: generated/

   flow_hierarchy
****
Swap
****

.. automodule:: networkx.algorithms.swap
.. autosummary::
   :toctree: generated/

   double_edge_swap
   connected_double_edge_swap

******
Triads
******

.. automodule:: networkx.algorithms.triads
.. autosummary::
   :toctree: generated/

   triadic_census
   random_triad
   triads_by_type
   triad_type
   all_triads
   all_triplets
**********
Efficiency
**********

.. automodule:: networkx.algorithms.efficiency_measures
.. autosummary::
   :toctree: generated/

   efficiency
   local_efficiency
   global_efficiency
**********
Tournament
**********

.. automodule:: networkx.algorithms.tournament
.. autosummary::
   :toctree: generated/

   hamiltonian_path
   is_reachable
   is_strongly_connected
   is_tournament
   random_tournament
   score_sequence
***********
Reciprocity
***********

.. automodule:: networkx.algorithms.reciprocity
.. autosummary::
   :toctree: generated/

   reciprocity
   overall_reciprocity
************
Simple Paths
************

.. automodule:: networkx.algorithms.simple_paths
.. autosummary::
   :toctree: generated/

   all_simple_paths
   all_simple_edge_paths
   is_simple_path
   shortest_simple_paths
.. _ismags:

****************
ISMAGS Algorithm
****************

.. automodule:: networkx.algorithms.isomorphism.ismags

ISMAGS object
-------------
.. currentmodule:: networkx.algorithms.isomorphism

.. autosummary::
   :toctree: generated/

   ISMAGS
********
Covering
********

.. automodule:: networkx.algorithms.covering
.. autosummary::
   :toctree: generated/

   min_edge_cover
   is_edge_cover
.. _algorithms:

**********
Algorithms
**********

.. currentmodule:: networkx

.. toctree::
   :maxdepth: 2

   approximation
   assortativity
   asteroidal
   bipartite
   boundary
   bridges
   centrality
   chains
   chordal
   clique
   clustering
   coloring
   communicability_alg
   community
   component
   connectivity
   core
   covering
   cycles
   cuts
   d_separation
   dag
   distance_measures
   distance_regular
   dominance
   dominating
   efficiency_measures
   euler
   flow
   graph_hashing
   graphical
   hierarchy
   hybrid
   isolates
   isomorphism
   link_analysis
   link_prediction
   lowest_common_ancestors
   matching
   minors
   mis
   non_randomness
   moral
   node_classification
   operators
   planarity
   planar_drawing
   reciprocity
   regular
   rich_club
   shortest_paths
   similarity
   simple_paths
   smallworld
   smetric
   sparsifiers
   structuralholes
   summarization
   swap
   threshold
   tournament
   traversal
   tree
   triads
   vitality
   voronoi
   wiener
*********
Planarity
*********

.. automodule:: networkx.algorithms.planarity
.. autosummary::
   :toctree: generated/

   check_planarity
.. autoclass:: PlanarEmbedding
   :members:****************
Threshold Graphs
****************

.. automodule:: networkx.algorithms.threshold
.. autosummary::
   :toctree: generated/

   find_threshold_graph
   is_threshold_graph
*************
Graph Hashing
*************

.. automodule:: networkx.algorithms.graph_hashing
.. autosummary::
   :toctree: generated/

   weisfeiler_lehman_graph_hash
   weisfeiler_lehman_subgraph_hashes
*******************
Similarity Measures
*******************

.. automodule:: networkx.algorithms.similarity
.. autosummary::
   :toctree: generated/

   graph_edit_distance
   optimal_edit_paths
   optimize_graph_edit_distance
   optimize_edit_paths
   simrank_similarity
   simrank_similarity_numpy
   panther_similarity
   generate_random_paths
***********************
Maximal independent set
***********************

.. automodule:: networkx.algorithms.mis
.. autosummary::
   :toctree: generated/

   maximal_independent_set

**********
Components
**********
.. automodule:: networkx.algorithms.components

Connectivity
------------
.. autosummary::
   :toctree: generated/

   is_connected
   number_connected_components
   connected_components
   node_connected_component

Strong connectivity
-------------------
.. autosummary::
   :toctree: generated/

   is_strongly_connected
   number_strongly_connected_components
   strongly_connected_components
   strongly_connected_components_recursive
   kosaraju_strongly_connected_components
   condensation

Weak connectivity
-----------------
.. autosummary::
   :toctree: generated/

   is_weakly_connected
   number_weakly_connected_components
   weakly_connected_components

Attracting components
---------------------
.. autosummary::
   :toctree: generated/

   is_attracting_component
   number_attracting_components
   attracting_components

Biconnected components
----------------------
.. autosummary::
   :toctree: generated/

   is_biconnected
   biconnected_components
   biconnected_component_edges
   articulation_points

Semiconnectedness
-----------------
.. autosummary::
   :toctree: generated/

   is_semiconnected
**************
non-randomness
**************

.. automodule:: networkx.algorithms.non_randomness
.. autosummary::
   :toctree: generated/

   non_randomness
*****************
Distance Measures
*****************

.. automodule:: networkx.algorithms.distance_measures
.. autosummary::
   :toctree: generated/

   barycenter
   center
   diameter
   eccentricity
   extrema_bounding
   periphery
   radius
   resistance_distance


.. _operators:

Operators
*********

.. automodule:: networkx.algorithms.operators.unary
.. autosummary::
   :toctree: generated/

   complement
   reverse

.. automodule:: networkx.algorithms.operators.binary
.. autosummary::
   :toctree: generated/

   compose
   union
   disjoint_union
   intersection
   difference
   symmetric_difference
   full_join


.. automodule:: networkx.algorithms.operators.all
.. autosummary::
   :toctree: generated/

   compose_all
   union_all
   disjoint_union_all
   intersection_all


.. automodule:: networkx.algorithms.operators.product
.. autosummary::
   :toctree: generated/

   cartesian_product
   lexicographic_product
   rooted_product
   strong_product
   tensor_product
   power
**************
Planar Drawing
**************

.. automodule:: networkx.algorithms.planar_drawing
.. autosummary::
   :toctree: generated/

   combinatorial_embedding_to_pos
******
Cycles
******

.. automodule:: networkx.algorithms.cycles
.. autosummary::
   :toctree: generated/

   cycle_basis
   simple_cycles
   find_cycle
   minimum_cycle_basis
.. _traversal:

Traversal
=========

.. toctree::
   :maxdepth: 2



Depth First Search
------------------
.. automodule:: networkx.algorithms.traversal.depth_first_search
.. autosummary::
   :toctree: generated/

   dfs_edges
   dfs_tree
   dfs_predecessors
   dfs_successors
   dfs_preorder_nodes
   dfs_postorder_nodes
   dfs_labeled_edges

Breadth First Search
--------------------
.. automodule:: networkx.algorithms.traversal.breadth_first_search
.. autosummary::
   :toctree: generated/

   bfs_edges
   bfs_tree
   bfs_predecessors
   bfs_successors
   descendants_at_distance

Beam search
-----------
.. automodule:: networkx.algorithms.traversal.beamsearch
.. autosummary::
   :toctree: generated/

   bfs_beam_edges


Depth First Search on Edges
---------------------------
.. automodule:: networkx.algorithms.traversal.edgedfs
.. autosummary::
   :toctree: generated/

   edge_dfs

Breadth First Search on Edges
-----------------------------
.. automodule:: networkx.algorithms.traversal.edgebfs
.. autosummary::
   :toctree: generated/

   edge_bfs
.. _chordal:

Chordal
=======

.. automodule:: networkx.algorithms.chordal
.. autosummary::
   :toctree: generated/

   is_chordal
   chordal_graph_cliques
   chordal_graph_treewidth
   complete_to_chordal_graph
   find_induced_nodes
****************
Structural holes
****************

.. automodule:: networkx.algorithms.structuralholes
.. autosummary::
   :toctree: generated/

   constraint
   effective_size
   local_constraint
********
Boundary
********

.. automodule:: networkx.algorithms.boundary
.. autosummary::
   :toctree: generated/

   edge_boundary
   node_boundary

****
Cuts
****

.. automodule:: networkx.algorithms.cuts
.. autosummary::
   :toctree: generated/

   boundary_expansion
   conductance
   cut_size
   edge_expansion
   mixing_expansion
   node_expansion
   normalized_cut_size
   volume
***************
Communicability
***************

.. automodule:: networkx.algorithms.communicability_alg
.. autosummary::
   :toctree: generated/

   communicability
   communicability_exp
********
Vitality
********

.. automodule:: networkx.algorithms.vitality
.. autosummary::
   :toctree: generated/

   closeness_vitality
.. _graph:

=========================================
Graph---Undirected graphs with self loops
=========================================

Overview
========
.. currentmodule:: networkx
.. autoclass:: Graph

Methods
=======

Adding and removing nodes and edges
-----------------------------------

.. autosummary::
   :toctree: generated/

   Graph.__init__
   Graph.add_node
   Graph.add_nodes_from
   Graph.remove_node
   Graph.remove_nodes_from
   Graph.add_edge
   Graph.add_edges_from
   Graph.add_weighted_edges_from
   Graph.remove_edge
   Graph.remove_edges_from
   Graph.update
   Graph.clear
   Graph.clear_edges



Reporting nodes edges and neighbors
-----------------------------------
.. autosummary::
   :toctree: generated/

   Graph.nodes
   Graph.__iter__
   Graph.has_node
   Graph.__contains__
   Graph.edges
   Graph.has_edge
   Graph.get_edge_data
   Graph.neighbors
   Graph.adj
   Graph.__getitem__
   Graph.adjacency
   Graph.nbunch_iter



Counting nodes edges and neighbors
----------------------------------
.. autosummary::
   :toctree: generated/

   Graph.order
   Graph.number_of_nodes
   Graph.__len__
   Graph.degree
   Graph.size
   Graph.number_of_edges


Making copies and subgraphs
---------------------------
.. autosummary::
   :toctree: generated/

   Graph.copy
   Graph.to_undirected
   Graph.to_directed
   Graph.subgraph
   Graph.edge_subgraph
.. _digraph:

=========================================
DiGraph---Directed graphs with self loops
=========================================

Overview
========
.. currentmodule:: networkx
.. autoclass:: DiGraph

Methods
=======

Adding and removing nodes and edges
-----------------------------------

.. autosummary::
   :toctree: generated/

   DiGraph.__init__
   DiGraph.add_node
   DiGraph.add_nodes_from
   DiGraph.remove_node
   DiGraph.remove_nodes_from
   DiGraph.add_edge
   DiGraph.add_edges_from
   DiGraph.add_weighted_edges_from
   DiGraph.remove_edge
   DiGraph.remove_edges_from
   DiGraph.update
   DiGraph.clear
   DiGraph.clear_edges



Reporting nodes edges and neighbors
-----------------------------------
.. autosummary::
   :toctree: generated/

   DiGraph.nodes
   DiGraph.__iter__
   DiGraph.has_node
   DiGraph.__contains__
   DiGraph.edges
   DiGraph.out_edges
   DiGraph.in_edges
   DiGraph.has_edge
   DiGraph.get_edge_data
   DiGraph.neighbors
   DiGraph.adj
   DiGraph.__getitem__
   DiGraph.successors
   DiGraph.succ
   DiGraph.predecessors
   DiGraph.pred
   DiGraph.adjacency
   DiGraph.nbunch_iter


Counting nodes edges and neighbors
----------------------------------
.. autosummary::
   :toctree: generated/

   DiGraph.order
   DiGraph.number_of_nodes
   DiGraph.__len__
   DiGraph.degree
   DiGraph.in_degree
   DiGraph.out_degree
   DiGraph.size
   DiGraph.number_of_edges


Making copies and subgraphs
---------------------------
.. autosummary::
   :toctree: generated/

   DiGraph.copy
   DiGraph.to_undirected
   DiGraph.to_directed
   DiGraph.subgraph
   DiGraph.edge_subgraph
   DiGraph.reverse
.. _multigraph:

=================================================================
MultiGraph---Undirected graphs with self loops and parallel edges
=================================================================

Overview
========
.. currentmodule:: networkx
.. autoclass:: MultiGraph

Methods
=======

Adding and removing nodes and edges
-----------------------------------

.. autosummary::
   :toctree: generated/

   MultiGraph.__init__
   MultiGraph.add_node
   MultiGraph.add_nodes_from
   MultiGraph.remove_node
   MultiGraph.remove_nodes_from
   MultiGraph.add_edge
   MultiGraph.add_edges_from
   MultiGraph.add_weighted_edges_from
   MultiGraph.new_edge_key
   MultiGraph.remove_edge
   MultiGraph.remove_edges_from
   MultiGraph.update
   MultiGraph.clear
   MultiGraph.clear_edges



Reporting nodes edges and neighbors
-----------------------------------
.. autosummary::
   :toctree: generated/

   MultiGraph.nodes
   MultiGraph.__iter__
   MultiGraph.has_node
   MultiGraph.__contains__
   MultiGraph.edges
   MultiGraph.has_edge
   MultiGraph.get_edge_data
   MultiGraph.neighbors
   MultiGraph.adj
   MultiGraph.__getitem__
   MultiGraph.adjacency
   MultiGraph.nbunch_iter



Counting nodes edges and neighbors
----------------------------------
.. autosummary::
   :toctree: generated/

   MultiGraph.order
   MultiGraph.number_of_nodes
   MultiGraph.__len__
   MultiGraph.degree
   MultiGraph.size
   MultiGraph.number_of_edges


Making copies and subgraphs
---------------------------
.. autosummary::
   :toctree: generated/

   MultiGraph.copy
   MultiGraph.to_undirected
   MultiGraph.to_directed
   MultiGraph.subgraph
   MultiGraph.edge_subgraph
.. _ordered:

============================================
Ordered Graphs---Consistently ordered graphs
============================================

.. automodule:: networkx.classes.ordered

.. currentmodule:: networkx
.. autoclass:: OrderedGraph
.. autoclass:: OrderedDiGraph
.. autoclass:: OrderedMultiGraph
.. autoclass:: OrderedMultiDiGraph
.. _classes:

***********
Graph types
***********

NetworkX provides data structures and methods for storing graphs.

All NetworkX graph classes allow (hashable) Python objects as nodes
and any Python object can be assigned as an edge attribute.

The choice of graph class depends on the structure of the
graph you want to represent.

Which graph class should I use?
===============================

+----------------+------------+--------------------+------------------------+
| Networkx Class | Type       | Self-loops allowed | Parallel edges allowed |
+================+============+====================+========================+
| Graph          | undirected | Yes                | No                     |
+----------------+------------+--------------------+------------------------+
| DiGraph        | directed   | Yes                | No                     |
+----------------+------------+--------------------+------------------------+
| MultiGraph     | undirected | Yes                | Yes                    |
+----------------+------------+--------------------+------------------------+
| MultiDiGraph   | directed   | Yes                | Yes                    |
+----------------+------------+--------------------+------------------------+

Basic graph types
=================

.. toctree::
   :maxdepth: 2

   graph
   digraph
   multigraph
   multidigraph
   ordered

.. note:: NetworkX uses `dicts` to store the nodes and neighbors in a graph.
   So the reporting of nodes and edges for the base graph classes will not
   necessarily be consistent across versions and platforms.  If you need the
   order of nodes and edges to be consistent (e.g., when writing automated
   tests), please see :class:`~networkx.OrderedGraph`,
   :class:`~networkx.OrderedDiGraph`, :class:`~networkx.OrderedMultiGraph`,
   or :class:`~networkx.OrderedMultiDiGraph`, which behave like the base
   graph classes but give a consistent order for reporting of nodes and edges.

Graph Views
===========

.. automodule:: networkx.classes.graphviews
.. autosummary::
   :toctree: generated/

   generic_graph_view
   subgraph_view
   reverse_view

Core Views
==========

.. automodule:: networkx.classes.coreviews
.. autosummary::
   :toctree: generated/

   AtlasView
   AdjacencyView
   MultiAdjacencyView
   UnionAtlas
   UnionAdjacency
   UnionMultiInner
   UnionMultiAdjacency
   FilterAtlas
   FilterAdjacency
   FilterMultiInner
   FilterMultiAdjacency

Filters
=======

.. note:: Filters can be used with views to restrict the view (or expand it).
   They can filter nodes or filter edges. These examples are intended to help
   you build new ones. They may instead contain all the filters you ever need.

.. automodule:: networkx.classes.filters
.. autosummary::
   :toctree: generated/

   no_filter
   hide_nodes
   hide_edges
   hide_diedges
   hide_multidiedges
   hide_multiedges
   show_nodes
   show_edges
   show_diedges
   show_multidiedges
   show_multiedges
.. _multidigraph:


=================================================================
MultiDiGraph---Directed graphs with self loops and parallel edges
=================================================================

Overview
========
.. currentmodule:: networkx
.. autoclass:: MultiDiGraph

Methods
=======

Adding and Removing Nodes and Edges
-----------------------------------

.. autosummary::
   :toctree: generated/

   MultiDiGraph.__init__
   MultiDiGraph.add_node
   MultiDiGraph.add_nodes_from
   MultiDiGraph.remove_node
   MultiDiGraph.remove_nodes_from
   MultiDiGraph.add_edge
   MultiDiGraph.add_edges_from
   MultiDiGraph.add_weighted_edges_from
   MultiDiGraph.new_edge_key
   MultiDiGraph.remove_edge
   MultiDiGraph.remove_edges_from
   MultiDiGraph.update
   MultiDiGraph.clear
   MultiDiGraph.clear_edges



Reporting nodes edges and neighbors
-----------------------------------
.. autosummary::
   :toctree: generated/

   MultiDiGraph.nodes
   MultiDiGraph.__iter__
   MultiDiGraph.has_node
   MultiDiGraph.__contains__
   MultiDiGraph.edges
   MultiDiGraph.out_edges
   MultiDiGraph.in_edges
   MultiDiGraph.has_edge
   MultiDiGraph.get_edge_data
   MultiDiGraph.neighbors
   MultiDiGraph.adj
   MultiDiGraph.__getitem__
   MultiDiGraph.successors
   MultiDiGraph.succ
   MultiDiGraph.predecessors
   MultiDiGraph.succ
   MultiDiGraph.adjacency
   MultiDiGraph.nbunch_iter


Counting nodes edges and neighbors
----------------------------------
.. autosummary::
   :toctree: generated/

   MultiDiGraph.order
   MultiDiGraph.number_of_nodes
   MultiDiGraph.__len__
   MultiDiGraph.degree
   MultiDiGraph.in_degree
   MultiDiGraph.out_degree
   MultiDiGraph.size
   MultiDiGraph.number_of_edges

Making copies and subgraphs
---------------------------
.. autosummary::
   :toctree: generated/

   MultiDiGraph.copy
   MultiDiGraph.to_undirected
   MultiDiGraph.to_directed
   MultiDiGraph.subgraph
   MultiDiGraph.edge_subgraph
   MultiDiGraph.reverse
