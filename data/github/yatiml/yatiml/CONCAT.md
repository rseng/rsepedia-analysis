############################
Contributing guidelines
############################

We welcome any kind of contribution to our software, from simple comment or question to a full fledged `pull request <https://help.github.com/articles/about-pull-requests/>`_. Please read and follow our `Code of Conduct <CODE_OF_CONDUCT.rst>`_.

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

You have a question
*******************

1. use the search functionality `here <https://github.com/yatiml/yatiml/issues>`_ to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

1. use the search functionality `here <https://github.com/yatiml/yatiml/issues>`_ to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the `SHA hashcode <https://help.github.com/articles/autolinked-references-and-urls/#commit-shas>`_ of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

You want to make some kind of change to the code base
*****************************************************

1. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions `here <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`_ and `here <https://help.github.com/articles/syncing-a-fork/>`_);
1. make sure the existing tests still work by running ``python setup.py test``;
1. add your own tests (if necessary);
1. update or expand the documentation;
1. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the YAtiML repository on GitHub;
1. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`_.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
##########
Change Log
##########

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.

0.8.0
*****

Incompatible changes
--------------------

* Accept explicit tags only if compatible with the recognised type(s)

New functionality
-----------------

* Support for untyped documents and attributes
* Support for Any-typed documents and attributes
* Support for Python dataclasses


Fixes
-----

* Dumping of OrderedDict to a file (but not to a string) produced a stray
  !!omap.
* Various fixes and improvements to development infrastructure


Removed
-------

* Official support for Python 3.5, which is no longer supported upstream. It
  will probably still work, but getting anything to install on 3.5 is getting to
  be pretty difficult so it's probably time to upgrade.


0.7.0
*****

Incompatible changes
--------------------

* Use seasoning functions only on the class they're defined on

New functionality
-----------------

* New yatiml.String to mark string-serialisable classes
* User-defined strings may be used as dictionary keys
* Support for index mappings
* Support for latest ruamel.yaml
* Documentation improvements


0.6.1
*****

Incompatible changes
--------------------

* Use datetime.date instead of datetime.datetime

New functionality
-----------------

* Support for loading and dumping pathlib.Path objects
* Support for Python 3.9


0.6.0
*****

New functionality
-----------------

* New make_loader and make_dumper functions improve ease-of-use
* JSON support
* Support for Mapping and Sequence types
* UnknownNode.require_attribute_value_not() function
* Node.remove_attributes_with_default_values() function
* Recipe for seasoning Enums

Fixes
-----

* Various documentation improvements
* Better error message if constructor raises


0.5.1
*****

Fixes
-----

* Fixed support for Python 3.5.1 (again, sorry)

0.5.0
*****

Incompatible changes
--------------------

* yatiml_* methods should now be called _yatiml_*
* Dropped support for Python 3.4, which is end-of-life

Fixes
-----

* Savourised classes in lists and dicts now load correctly
* Fixed compatibility with the latest versions of ruamel.yaml
* Fixed support for Python 3.5.1

0.4.2
*****

Fixes
-----

* Don't generate cross-references for enum values
* Various small fixes

0.4.1
*****

New functionality
-----------------

* Added fix_union_bool type for fixing Union[int, bool] on Python < 3.7
* Added support for Python 3.7

Fixes
-----

* Return scalar values with the correct type

0.4.0
*****

New functionality
-----------------

* Extended map_to_seq seasoning
* Support for YAML timestamp / Python datetime
* Support for YAML keys with dashes

Fixes
-----

* Much improved error messages

0.3.0
*****

New functionality
-----------------

* Support for classes that are represented by a string in the YAML file
* New unified yatiml.Node interface (API change)

Fixes
-----

* Small improvements to documentation
* Miscellaneous small fixes

0.2.0
*****

New functionality
-----------------

* Support for enumerations
* Support for user-defined string types

Fixes
-----

* Various small tooling fixes
* Some refactoring

0.1.0
*****

* Initial release with basic functionality
###############################################################################
Contributor Covenant Code of Conduct
###############################################################################

Our Pledge
**********

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

Our Standards
*************

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
  advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

Our Responsibilities
********************

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

Scope
*****

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

Enforcement
***********

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at l.veen@esciencecenter.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

Attribution
***********

This Code of Conduct is adapted from the `Contributor Covenant <https://www.contributor-covenant.org>`_, version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html
.. image:: https://readthedocs.org/projects/yatiml/badge/?version=develop
    :target: https://yatiml.readthedocs.io/en/latest/?badge=develop
    :alt: Documentation Build Status

.. image:: https://github.com/yatiml/yatiml/workflows/continuous_integration/badge.svg
    :target: https://github.com/yatiml/yatiml/actions
    :alt: Build Status

.. image:: https://app.codacy.com/project/badge/Grade/bca7a121d9c742d2905eae08a75676c3
    :target: https://www.codacy.com/gh/yatiml/yatiml/dashboard
    :alt: Codacy Grade

.. image:: https://app.codacy.com/project/badge/Coverage/bca7a121d9c742d2905eae08a75676c3
    :target: https://www.codacy.com/gh/yatiml/yatiml/dashboard
    :alt: Code Coverage

.. image:: https://requires.io/github/yatiml/yatiml/requirements.svg?branch=develop
    :target: https://requires.io/github/yatiml/yatiml/requirements/?branch=develop
    :alt: Requirements Status

.. image:: https://zenodo.org/badge/147202299.svg
   :target: https://zenodo.org/badge/latestdoi/147202299

.. image:: https://img.shields.io/badge/rsd-yatiml-00a3e3.svg
   :target: https://www.research-software.nl/software/yatiml

################################################################################
YAtiML
################################################################################

YAML-based file formats can be very handy, as YAML is easy to write by humans,
and parsing support for it is widely available. Just read your YAML file into a
document structure (a tree of nested dicts and lists), and manipulate that in
your code.

As long as that YAML file contains exactly what you expect, that works fine.
But if it contains a mistake, then you're likely to crash the program with a
cryptic error message, or worse (especially if the YAML file was loaded from the
Internet) it may do something unexpected.

To avoid that, you can validate your YAML using various schema checkers. You
write a description of what your YAML file must look like, then feed that to a
library which checks the incoming file against the description. That gives you a
better error message, but it's a lot of work.

YAtiML takes a different approach. Instead of a schema, you write a Python
class. You probably already know how to do that, so no need to learn anything.
YAtiML then generates loading and dumping functions for you, which convert
between YAML and Python objects. If needed, you can add some extra code to make
the YAML look nicer or implement special features.

YAtiML supports Python 3.6 and later.

If you use YAtiML for scientific work, we ask that you cite it. You can
`download a citation in various formats at the Research Software Directory
<https://www.research-software.nl/software/yatiml>`_.

Documentation and Help
**********************

Instructions on how to install and use YAtiML can be found in `the YAtiML
documentation <https://yatiml.readthedocs.io>`_.

Code of Conduct
---------------

Before describing where to ask questions or report bugs, we'd like to point out
that this project is governed by a code of conduct, as described in
CODE_OF_CONDUCT.rst, and we expect you to adhere to it. Please be nice to your
fellow humans.

Questions
---------

If you have a question that the documentation does not answer for you, then you
have found a bug in the documentation. We'd love to fix it, but we need a bit of
help from you to do so. Please do the following:

#. use the `search functionality <https://github.com/yatiml/yatiml/issues>`_
   to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

We'll answer your question, and improve the documentation where necessary.

Bugs
----

Like most software, YAtiML is made by humans, and we make mistakes. If you think
you've found a bug in YAtiML, please let us know! Reporting bugs goes as follows.

#. Use the `search functionality`_ to see if someone already filed the same
   issue.

#. If your issue search did not yield any relevant results, make a new issue.

   Please explain:
    - what you were trying to achieve,
    - what you did to make that happen,
    - what you expected the result to be,
    - what happened instead.

  It really helps to have the actual code for a simple example that demonstrates
  the issue, but excerpts and error messages and a description are welcome too.

#. Finally, apply any relevant labels to the newly created issue.

With that, we should be able to fix the problem.

License
*******

YAtiML is Copyright 2018-2021, Netherlands eScience Center, University of
Amsterdam, and VU University Amsterdam

Distributed under the Apache Software License 2.0.
.. _development:

Development
***********

To get the source code from GitHub, you can do

.. code-block:: console

  git clone git@github.com:yatiml/yatiml.git
  cd yatiml


Run tests (including coverage and type checking) with:

.. code-block:: console

  pip install tox
  tox


A local copy of the documentation can be generated using:

.. code-block:: console

  tox -e docs


Contributing
------------

If you want to contribute some improvements to YAtiML, please use the following
process:

#. (**important**) Announce your plan to the rest of the community *before you
   start working*. This announcement should be in the form of a (new) issue.
#. (**important**) wait until some kind of consensus is reached about your idea
   being a good idea.
#. If needed, fork the repository to your own Github profile and create your
   own feature branch off of the latest ``develop`` commit. While working on
   your feature branch, make sure to stay up to date with the ``develop``
   branch by pulling in changes, possibly from the 'upstream' repository
   (follow the instructions `here
   <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`__ and
   `here <https://help.github.com/articles/syncing-a-fork/>`__).
#. Make sure the existing tests still work by running ``python setup.py test``.

#. Add your own tests (if necessary).

#. Update or expand the documentation.

#. Use ``yapf`` to fix the readability of your code style and ``isort``
   to format and group your imports.

#. `Push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to
   (your fork of) the YAtiML repository on GitHub.

#. create the pull request,
   e.g. following the instructions `here
   <https://help.github.com/articles/creating-a-pull-request/>`_.

In case you feel like you've made a valuable contribution, but you don't know
how to write or run tests for it, or how to generate the documentation: don't
let this discourage you from making the pull request; we can help you! Just go
ahead and submit the pull request, but keep in mind that you might be asked to
append additional commits to your pull request.


Making a release
----------------

YAtiML uses Git on GitHub for version management, using the `Git Flow`_
branching model. Making a release involves quite a few steps, so they're listed
here to help make the process more reliable; this information is really only
useful for the maintainers.

Check metadata
--------------

- Check the metadata in ``setup.py``, and update as necessary.

- Check the copyright date and owners in README.rst and docs/conf.py and update
  as necessary.


Update the changelog
....................

Each release should have an entry in the CHANGELOG.rst describing the new
features and fixed problems. Since we'll want to carry these entries forward,
we'll make them first, on the develop branch. Use the git log to get a list of
the changes, and switch to the development branch:

.. code-block:: bash

  git log <your favourite options>
  git checkout develop

and then edit CHANGELOG.rst and commit.

.. code-block:: bash

  git add CHANGELOG.rst
  git commit -m 'Add version x.y.z to the change log'

Make release branch
...................

To start the release process, make a release branch

.. code-block:: bash

  git checkout -b release-x.y.z develop

YAtiML uses `Semantic Versioning`_, so name the new version accordingly.

Update version
..............

Next, the version should be updated. There is a version tag in ``setup.py`` and
two for the documentation in ``docs/conf.py`` (search for ``version`` and
``release``). There is also an ``__version__`` in ``__init__.py``. On the
development branch, these should be set to ``x.y.z.dev0``, where ``x.y.z`` is
the expected next version. On the release branch, they should be set to
``x.y.z`` (with here the actual number of this release of course).

Check documentation
...................

Since we've just changed the documentation build configuration, the build should
be run locally to test:

.. code-block:: bash

  tox -e docs

It may give some warnings about missing references; they should disappear if
you run the command a second time. Next, point your web browser to
``docs/_build/html/index.html`` and verify that the documentation built
correctly. In particular, the new version number should be in the browser's
title bar as well as in the blue box on the top left of the page.

Run tests
.........

Before we make a commit, the tests should be run, and this is a good idea anyway
if we're making a release. So run ``tox`` and check that everything is in order.

Commit the version update
.........................

This is the usual Git poem:

.. code-block:: bash

  git add setup.py docs/conf.py yatiml/__init__.py
  git commit -m 'Set release version to x.y.z'
  git push --set-upstream origin release-x.y.z

This will trigger the Continuous Integration, so check that that's not giving
any errors while we're at it.

Fix badges
..........

The badges in the README.rst normally point to the development branch versions
of everything. For the master branch, they should point to the master version.
Note that for the ReadTheDocs badge, `develop` should be changed to `latest`,
and that for Codacy there is only one badge, so no change is needed.

.. code-block:: bash

  # edit README.rst
  git add README.rst
  git commit -m 'Update badges to point to master'
  git push

Merge into the master branch
............................

If all seems to be well, then we can merge the release branch into the master
branch and tag it, thus making a release, at least as far as Git Flow is
concerned. We use the ``-X theirs`` option here to resolve the merge conflict
caused by the version update that was done for the previous release, which we
don't have on this branch. The last command is to push the tag, which is
important for GitHub and GitHub integrations.

.. code-block:: bash

  git checkout master
  git merge --no-ff -X theirs release-x.y.z
  git tag -a x.y.z -m 'Release x.y.z'
  git push
  git push origin x.y.z

Build and release to PyPI
.........................

Finally, the new version needs to be built and uploaded to PyPI, so that people
can start using it. To build, use:

.. code-block:: bash

  python3 setup.py sdist bdist_wheel

Then, we can upload to the test instance of PyPI:

.. code-block:: bash

  twine upload --repository-url https://test.pypi.org/legacy/ dist/yatiml-x.y.z*

To test that we can install it, run this in a fresh virtualenv:

.. code-block:: bash

  python3 -m pip install --index-url https://test.pypi.org/simple/ yatiml

And if all seems well, we can upload to the real PyPI:

.. code-block:: bash

  twine upload dist/yatiml-x.y.z*

Make a GitHub Release
.....................

Go to Releases on the GitHub page and make a new release from the tag. For the
release notes, copy-paste from the CHANGELOG and convert from RST to Markdown.

Merge release branch back into develop
......................................

To continue developing, merge the release branch back into develop

.. code-block:: bash

  git checkout develop
  git merge --no-commit release-x.y.z
  git push

Make sure that the badges are set to develop, and that the version number is
set to the next expected version x.y.{z+1}.dev (it's fine if x.{y+1}.0 is what
ends up being released eventually). Then you can commit and continue developing.

.. _`Git Flow`: http://nvie.com/posts/a-successful-git-branching-model/
.. _`Semantic Versioning`: http://www.semver.org
Recipes
=======

Parsed classes
--------------

For some classes, the easiest way to write them in a YAML file is as a string
representing all the values in the class. For example, you may want to have a
namespaced name that looks like ``ns.subns.name`` in the YAML text format, but
on the Python side represent it as a class having one attribute `namespaces`
that is a list of namespaces, and an attribute `name` that is a string.

To do this, you need to override recognition to tell YAML to recognise a string
(because by default it will expect a mapping), and then to add a savorizing
function that parses the string and generates a mapping, attributes from which
will then be fed to your constructor. In order to save objects of your class as
strings, you'll need to add a sweetening function too. The complete solution
looks like this:

.. literalinclude:: examples/parsed_classes.py
  :caption: ``docs/examples/parsed_classes.py``
  :language: python

The tricky part here is in the savorize and sweeten functions. Savorize needs to
build a list of strings, which YAtiML doesn't help with, so it needs to
construct ruamel.yaml objects directly. For each namespace item, it builds a
``yaml.ScalarNode``, which represents a scalar and has a tag to describe the type,
and a value, a string. It also requires a start and end mark, for which we use
dummy values, as this node was generated and is therefore not in the file. The
ruamel.yaml library will raise an Exception if you do not add those. The item
nodes are then added to a ``yaml.SequenceNode``, and the whole thing set as the
value of the ``namespaces`` attribute.

Sweeten does the reverse of course, getting a list of :class:`yatiml.Node`
objects representing the items in the ``namespaces`` attribute, extracting the
string values using :meth:`yatiml.Node.get_value`, then joining them with
periods and finally combining them with the name. Since we're only altering the
top-level node here, we do not need to build a ``yaml.ScalarNode`` ourselves but
can just use :meth:`yatiml.Node.set_value`.


Timestamps and dates
--------------------

YAML has a `timestamp` type, which represents a point in time. The
`ruamel.yaml` library parses this into a python `datetime.date` object, and
will serialise such an object back to a YAML `timestamp`. YAtiML supports this
as well, so all you need to do to use a timestamp or a date is to use
`datetime.date` in your class definition.

Note that the object created by YAtiML may be an instance of `datetime.date` (if
no time is given) or an instance of `datetime.datetime` (if a time is given)
which is a subclass of `datetime.date`. Since Python does not have a
date-without-time type, you cannot currently specify in the type that you want
only a date, without a time attached to it.

If this is an attribute in a class, and date-with-time is not a legal value,
then you should add a check to the __init__ method that raises an exception if
the given value is an instance of `datetime.datetime`. That way, you can't
accidentally make an instance of the class in Python with an incorrect value
either.


Dashed keys
-----------

Some YAML-based formats (like CFF) use dashes in their mapping keys. This is a
problem for YAtiML, because keys get mapped to parameters of ``__init__``,
which are identifiers, and those are not allowed to contain dashes in Python. So
some kind of conversion will have to be made. YAtiML's seasoning mechanism is
just the way to do it: :class:`yatiml.Node` has two methods to convert all
dashes in a mapping's keys to underscores and back:
:meth:`unders_to_dashes_in_keys()` and :meth:`dashes_to_unders_in_keys()`, so
all you need to do is use underscores instead of dashes when defining your
classes, and add seasoning functions. Here's an example:

.. literalinclude:: examples/dashed_keys.py
  :caption: ``docs/examples/dashed_keys.py``
  :language: python

If you've been paying very close attention, then you may be wondering why this
example passes through the recognition stage. After all, the names of the keys
do not match those of the ``__init__`` parameters. YAtiML is a bit flexible in
this regard, and will match a key to a parameter if it is identical after dashes
have been replaced with underscores. The flexibility is only in the recognition
stage, not in the type checking stage, so you do need the seasoning functions.
(The reason to not completely automate this is that YAtiML cannot know if the
YAML side should have dashes or underscores. So you need to specify this
somehow in order to be able to dump correctly, and then it's better to specify
it on loading as well for symmetry.)


.. _seasoning_enumerations:

Seasoning enumerations
----------------------

By default, YAtiML will use an enum member's name to write to the YAML file, and
that's what it will recognise on loading as well. Sometimes, that's not what you
want however. Maybe you want to use the values, or you want to have the names on
the Python side in uppercase (because PEP-8 says so) while you want to use
a lower-case version in the YAML file. In that case, you can apply YAtiML's
seasoning mechanisms like this:

.. literalinclude:: examples/enum_use_values.py
  :caption: ``docs/examples/enum_use_values.py``
  :language: python

or like this:

.. literalinclude:: examples/enum_lowercase.py
  :caption: ``docs/examples/enum_lowercase.py``
  :language: python
Why YAtiML?
==============

YAML-based file formats can be very handy, as YAML is easy to write by humans,
and parsing support for it is widely available. Just read your YAML file into a
document structure (a tree of nested dicts and lists), and manipulate that in
your code.

As long as that YAML file contains exactly what you expect, that works fine.
But if it contains a mistake, then you're likely to crash the program with a
cryptic error message, or worse (especially if the YAML file was loaded from the
Internet) it may do something unexpected.

To avoid that, you can validate your YAML using various schema checkers. You
write a description of what your YAML file must look like, then feed that to a
library which checks the incoming file against the description. That gives you a
better error message, but it's a lot of work.

YAtiML takes a different approach. Instead of a schema, you write a Python
class. You probably already know how to do that, so no need to learn anything.
YAtiML then generates loading and dumping functions for you, which convert
between YAML and Python objects. On loading, the input is checked to ensure that
it matches the intended type, using standard Python type annotations on the
class. If there is an error, a (very!) nice error message is produced. Default
values, specified as usual in the ``__init__`` method, are applied
automatically, saving you another big headache.

If you want to go further, and create a more complex YAML-based file format
like Docker Compose-files or the Common Workflow Language, then YAtiML has you
covered too. It lets you hook into the loading and dumping processes, modifying
the formatting of the YAML file without affecting the Python side, which lets
you implement all sorts of nice formatting features. YAtiML supports class
hierarchies, enumerations, and extension points (parts of the YAML document
where anything goes) as well.
Problem solving
===============

This section collects some useful information about solving problems you may
encounter when working with YAtiML.

Enabling logging
----------------

YAtiML uses the standard Python logging system, using a logger named ``yatiml``.
YAtiML produces log messages at the INFO and DEBUG levels, which makes them
invisible at Python's default log level of WARNING. To be able to see what
YAtiML does, you can either lower the general Python log level (which may cause
other parts of your program to produce (more) log output as well), or you can
lower YAtiML's log level specifically. To do the latter, use this:

.. code-block:: python

  import logging
  import yatiml

  yatiml.logger.setLevel(logging.INFO)

Or you can use ``logging.DEBUG`` for very detailed output.

It seems that you may have to do a ``logging.debug()`` call to get any output at
all, maybe because that causes Python to set up something it needs. There's
probably a good explanation and a better fix for this. If you know, please
contribute.

In order to understand how to interpret the output, it helps to have an idea of
how YAtiML processes a YAML file into Python objects. See `The YAtiML
pipeline`_ below for more on that.

Unions with bool
----------------

While defining your classes, you may want to have an attribute that can be of
multiple types. As described in the tutorial, you would use a ``Union`` for
this. For example something like

.. code-block:: python

  class Setting:
      def __init__(name: str,
                   value: Union[str, int, float, bool]
                   ) -> None:
          self.name = name
          self.value = value

is likely to occur in many YAML-based configuration file formats.

There is a problem with the above code, in that it will give an error message if
you try to read the following input on Python < 3.7, saying that it could not
recognise the bool value:

.. code-block:: yaml

  name: test
  value: true

(Scroll down for the fix, if you don't care for the explanation.)

Arguably, this is a bug in Python's type handling, and the developers of
Python's ``typing`` module seem to agree, because they have fixed this in Python
3.7. What happens here is that in Python, ``bool`` is a subtype of ``int``, in
other words, any ``bool`` value is also an ``int`` value. Furthermore, the
``Union`` class will automatically normalise the types it is passed by removing
redundant types. If you put in a type twice then one copy will be removed for
instance, and also, if you put in a type and also a subtype of that type, then
the subtype will be removed. This makes some sense: if every ``bool`` is an
``int``, then just ``int`` will already match boolean values, and ``bool`` is
redundant.

While this works for Python, it's problematic in YAML, where ``bool`` and
``int`` are unrelated types. Indeed, YAtiML will not accept a boolean value in
the YAML file if you declare the attribute to be an ``int``. And that's where
we get into trouble: Python normalises the above Union to ``Union[str, int,
float]``, and YAtiML reads this and generates an error if you feed it a boolean.

In Python 3.7, the behaviour of Union has changed. While mypy still does the
normalisation internally when checking types, the runtime Union object no longer
normalises. Since the runtime object is what YAtiML reads, this problem does
not occur on Python 3.7 (and hopefully versions after that, the `typing` module
is not entirely stable yet).

A fix for Python < 3.7
''''''''''''''''''''''

So, this is fixed in Python 3.7, but what if you're running on an older version?
In that case you need a work-around, and YAtiML provides one called
``bool_union_fix``. It works like this:

.. code-block:: python

  from yatiml import bool_union_fix

  class Setting:
      def __init__(name: str,
                   value: Union[str, int, float, bool, bool_union_fix]
                   ) -> None:
          self.name = name
          self.value = value

All you need to do is import ``bool_union_fix`` and add it to the list of types
in the ``Union``, and now things will work as expected (also in Python 3.7).

``bool_union_fix`` is essentially a dummy type that is recognised by YAtiML and
treated just like ``bool``. Since it's a separate type, it isn't merged into the
``int``, so it'll still be there for YAtiML to read. Note that you do need the
``bool`` in there as well, to avoid mypy complaining if you try to create a
Setting object in your code with a bool for the value attribute.


The YAtiML pipeline
-------------------

With plain PyYAML or ruamel.yaml, the loading process has two stages. First, the
text is parsed into a parse tree, which consists of nodes. Each node has a tag
and a value. Second, objects are constructed from the nodes, with the type of
the object decided based on the tag, and the contents of the object coming from
the value.

YAtiML inserts three additional stages in between the two existing ones:
recognition, savourising, and type checking.

Recognition determines, for each node, as which type YAtiML will try to process
it. This is mostly based on the object model given to the custom loader. In our
ongoing example, the value corresponding to the ``name`` attribute is expected
to be a string, so YAtiML will try to recognise only a string here. The ``age``
attribute has a union type, and for those YAtiML will look at the value given
and see if it matches one of the types in the union. If it matches exactly one,
it is recognised as that type. If it matches none of them, or multiple, an error
message is given.

When recognising a node that according to the object model is of a custom class
type, YAtiML will try to recognise a mapping node with keys and values according
to the custom class's ``__init__`` method's parameters. If the custom class has
subclasses which are also registered with the loader, then those will be
recognised as well at this point in the document. If both a class and its
subclass are matched, the node is recognised as being of the subclass, i.e.
the recognition process prefers the most derived class. If there are multiple
matching sibling subclasses, the node is declared ambiguous and an error is
raised. Recognition for a custom class can be overridden using a
``_yatiml_recognize()`` method.

Incidentally, a technical term for what the recognition process does is `type
inference`, which explains the name YAtiML: it inserts type inference in the
middle of the YAML processing pipeline.

The second and third stages, savourising and type checking, only apply to custom
classes. To savourise a recognised node, YAtiML calls that node's
``_yatiml_savorize()`` method, after calling those of its base classes, if any.
Savourising is entirely defined by the custom class, the default is to do
nothing. After savourising, the resulting mapping is type checked against the
``__init__`` signature, since Python does not do run-time type checking itself.
This is a safety feature, since the read-in YAML document will often be
untrusted, or if it is, at least a convenience feature, in that what you see in
the ``__init__`` signature is guaranteed to be what you get, thus applying the
principle of least surprise.

Note that no type check is done for built-in types, but for built-in types the
default recognition process is effectively a type check, and it cannot be
overridden. Another way of looking at the type check for custom classes is that
it reduces the requirements on custom recognition functions: they need to merely
disambiguate between derived classes and in unions, rather than performing a
complete type check. That makes it easier to write them.

Error messages
--------------

This section contains some error messages that you may encounter when using
YAtiML, and potential solutions to try if you do. If you run into an error that
you cannot figure out, please make an issue describing the error message and
what you are doing (a short example really helps!). Contributions directly to
the documentation are of course also welcome! See the
:ref:`development` section for information on how to contribute.


_yatiml_recognize missing required argument
'''''''''''''''''''''''''''''''''''''''''''

If you get

.. code-block:: python

  TypeError: _yatiml_recognize() missing 1 required positional argument: 'node'

or

.. code-block:: python

  TypeError: yatiml_savorize() missing 1 required positional argument: 'node'

then you have probably forgotten to add the ``@classmethod`` decorator to your
``_yatiml_recognize()`` or ``_yatiml_savorize()`` function.
.. _installing:

Installing YAtiML
=================

YAtiML is available from PyPI, and you can install it using pip:

.. code-block:: console

  pip install yatiml


or add it to the dependencies of your project in your usual way.

Changes between versions are listed in the file CHANGELOG.rst. YAtiML adheres
to `Semantic Versioning <http://semver.org/>`_. Starting from version 1.0, you
will be able to count on a stable API as long as you pin your dependency to a
major version.
Basic Tutorial
==============

YAtiML is a library for reading and writing YAML from Python.

This tutorial shows how to use YAtiML by example. You can find the example
programs shown below in the ``docs/examples/`` directory in the repository.

A first example
---------------

Let's say that we're organising a drawing contest for kids, and are tracking
submissions in YAML files. We'll need to read the files into a Python program in
order to process them. Here's how to do that with YAtiML:

For the example to run, make sure that you have :ref:`installed YAtiML
<installing>` first.

.. literalinclude:: examples/load_any_yaml.py
  :caption: ``docs/examples/load_any_yaml.py``
  :language: python

If you run this program, it will output

.. code-block:: none

  ordereddict([('name', 'Janice'), ('age', 6)])


Here is the example again one line at a time.

.. code-block:: python

  import yatiml


This loads the YAtiML package, so that we can use it in our Python script.

.. code-block:: python

  yaml_text = (
          'name: Janice\n'
          'age: 6\n')


This makes a string with a YAML document in it. The parentheses are so we can
split the string over multiple lines, and we need a ``\n`` at the end of each
line to explicitly mark it as the end, otherwise everything will be glued
together on a single line and we end up with invalid YAML.

.. code-block:: python

  load = yatiml.load_function()


To load our document from the string, we need a load function. YAtiML doesn't
have a built-in load function. Instead, it makes a custom load function just for
you, if you call :meth:`yatiml.load_function`. We'll call the result ``load``.

This probably looks a bit funny, but it will become clear why we're doing this
in the next example.

.. code-block:: python

  doc = load(yaml_text)

  print(doc)


Here we call our shiny new load function to load the YAML into a Python object.
Then we print the result so that we can see what happened. We will get

.. code-block:: none

  ordereddict([('name', 'Janice'), ('age', 6)])


This again looks a bit funny, but this is almost the same thing as the
dictionary ``{'name': 'Janice', 'age': 6}`` that you probably expected.

Until recently, Python dictionaries held their entries in random order. For
accessing items that's not a problem, because you look them up based on the key
anyway. But having the lines of a YAML file reorganised in a random order can
make the file really hard to read! So YAtiML reads the file into a special
ordered dictionary, which preserves the order. That way, you can save it again
later without making a big mess. Otherwise it works just like a plain Python
``dict``, so you can do ``doc['name']`` and so on as usual.

Checking the input
------------------

In the example above, we didn't specify any constraints on what the input should
look like. This is often inconvenient. For example, if we have an age limit of
12 on our drawing contest and want to write a program that reads in the YAML
file for each submission, checks the age and then prints out the names of any
kids that are too old, then we really need to have both the name and the age in
each file, or it's not going to work.

The example code above will happily read any input. If there's a list of numbers
in the file, then ``doc`` will hold a list of numbers instead of a ``dict``, and
your program will probably crash somewhere with an error ``TypeError: list
indices must be integers or slices, not str`` and then you get to figure out
what went wrong and where.

It would be much better if we could check that our input is a really a
dictionary with keys ``name`` and ``age``. We could do that by hand, after
reading, but with YAtiML there's a better way to do it. We're going to make a
Python class that shows what the YAML should look like:

.. literalinclude:: examples/untyped_class.py
  :caption: ``docs/examples/untyped_class.py``
  :language: python


The main new bit of this example is the ``Submission`` class:

.. code-block:: python

  class Submission:
      def __init__(self, name, age):
          self.name = name
          self.age = age


This creates a Python class named ``Submission``. If you've never seen one, a
class is basically a group of variables, in this case ``name`` and ``age``.
Classes also have an *init function* with the special name ``__init__`` which is
used to create a variable containing an object holding those variables. So here
we have a class named ``Submission``. It can be used like this:

.. code-block:: python

  submission = Submission('Janice', 6)
  print(submission.name)    # prints Janice
  print(submission.age)     # prints 6

  submission.age = 7
  print(submission.age)     # prints 7


Now, we can pass this class to YAtiML when we ask it to create our load
function:

.. code-block:: python

  load = yatiml.load_function(Submission)


YAtiML will now create a load function for us that expects to read in a
dictionary containing keys ``name`` and ``age``.

We use the load function as before, and it will read the YAML file and convert
it into a ``Submission`` object. We can check that we really got one using
``type()``, and inspect the name and age of our contestant.

Of course, we got exactly the input we expected, so in this case everything went
fine. What if there's an error? Then you get an error message.

.. admonition:: Exercise

  Change the input in various ways in the previous example, and see what error
  messages you get when you try to load the incorrect input.


Checking types
--------------

If you have played around a little bit with the previous example, then you may
have noticed that there's a certain kind of problem that is not detected when
you load the YAML input into a ``Submission`` object, and that is that the
values for ``name`` and ``age`` may not be of the right type. For example,
someone could write their age as ``six`` instead of as ``6``, and you would
suddenly have a string where you expected a number. That would almost certainly
mess up the ``submission.age <= 12`` in your age check!

So it would be better if we could make sure that the inputs are of the right
type too, and give an error on loading if they are not. Here's how to do that:

.. literalinclude:: examples/typed_class.py
  :caption: ``docs/examples/untyped_class.py``
  :language: python


This example is almost the same as the previous one, except that the
``__init__`` function of our ``Submission`` class now has some *type
annotations*: instead of ``name`` it says ``name: str`` and instead of ``age``
it says ``age: int``. That is all it takes to make sure that any values given
for those keys in the YAML file are checked. (There's also ``-> None`` at the
end, which specifies that the function does not return anything. YAtiML ignores
this bit, and so can you if you want to.)

.. admonition:: Exercise

  Try changing the input to use values of a different type and see what happens.


``int`` and ``str`` are standard Python types, and adding them to the function
parameters as in the example is standard Python. For decimal numbers, you can
use ``float`` and for truth values (e.g. true, false, yes, no) the type
``bool``.

Lists and dicts are also supported, but they require some special types from the
standard Python ``typing`` package. For example, to allow multiple contestants
to make a drawing together, we could allow a list of strings for the ``name``
field, and a dictionary mapping each name to the corresponding age for ``age``.
That would look like this.

.. literalinclude:: examples/collaborative_submissions.py
  :caption: ``docs/examples/collaborative_submissions.py``
  :language: python


For dates you can use ``date`` from the ``datetime`` package, and if you need to
read the location of a file from a YAML file then you can use ``Path`` from
Python's ``pathlib``. If you want to explicitly accept any kind of YAML, then
you can use ``Any`` from ``typing``, which is the same as not specifying a type
at all like we did in the beginning.

Finally, ``Union`` from ``typing`` makes it possible to accept multiple
different types. Try this for example:

.. literalinclude:: examples/custom_class.py
  :caption: ``docs/examples/custom_class.py``
  :language: python


and see what YAML inputs it will accept.

Default values
--------------

One of the issues you will run into when implementing a complex YAML-based
format by hand, is default values. For example in a configuration file, it is
often much easier if the users can completely omit any options for which a
default value suffices. If you have nested optional structures (e.g. users are
allowed to omit an entire dictionary if its attributes have all been omitted),
then processing the data becomes a tedious set of nested ifs. In YAtiML,
default values are easy: since ``__init__`` parameters map to attributes, all
you have to do is declare a parameter with a default value:

.. literalinclude:: examples/default_values.py
  :caption: ``docs/examples/default_values.py``
  :language: python

Here we have added the tool that was used as an argument with a default value.
If the YAML file contains a key ``tool`` with a string value, that value will be
passed to the ``__init__`` method. If the key ``tool`` exists, but the value is
not of type string, a ``RecognitionError`` is raised. If the key is missing, the
default value is used.

Note that in this case, the ``tool`` attribute is optional in the YAML file, but
not in the class: every object of type ``Submission`` has to have a value for
``tool`` that is not ``None``. This allows you to conveniently skip the check,
which gets rid of those nested ifs if you have nested optional entries in your
YAML file.

However, you may want to make the attribute optional in the class as well, and
perhaps set ``None`` as the default value. That is done like this:

.. literalinclude:: examples/optional_attribute.py
  :caption: ``docs/examples/optional_attribute.py``
  :language: python

Now the value of a ``Submission`` object's ``tool`` attribute can be ``None``,
and it will be if that attribute is omitted in the YAML mapping. Note that this
definition is entirely standard Python 3, there is nothing YAtiML-specific in
it.

Saving to YAML
--------------

There is more to be said about loading YAML files with YAtiML, but let's first
have a look at saving objects back to YAML, or dumping as PyYAML and ruamel.yaml
call it. The code for this is a mirror image of the loading code:

.. literalinclude:: examples/saving.py
  :caption: ``docs/examples/saving.py``
  :language: python

And as expected, it outputs:

.. code-block:: none

  name: Youssou
  age: 7
  tool: pencils

YAtiML expects a public attribute with the same name for each parameter in the
``__init__`` method to exist, and will use its value in saving. This can be
overridden, see :ref:`hiding-attributes` below.

Note that the attributes are in the order of the parameters of the ``__init__``
method. YAtiML always outputs attributes in this order, even if the object was
read in with YAtiML from a YAML file and originally had a different order. While
it would be nice to do full round-trip formatting of the input YAML, support for
this in the ruamel.yaml library used by YAtiML is still developing, so for now
this is what YAtiML does.

:meth:`yatiml.dumps_function` creates a function that converts objects to a
string. If you want to write the output to a file directly, you can use
:meth:`yatiml.dump_function` instead to create a function that can do that.

As an example of the advantage of using YAtiML, saving a Submission document
with PyYAML or ruamel.yaml gives this:

.. code-block:: none

  !!python/object:__main__.Submission {age: 7, name: Youssou, tool: pencils}

which is not nearly as nice to read or write. (To be fair, ruamel.yaml can do a
bit nicer than this with its RoundTripDumper, which YAtiML uses, but the tag
with the exclamation marks remains.)

Saving to JSON
--------------

YAML is a superset of JSON, so YAtiML can read JSON files. If you want to save
JSON as well, then you can use :meth:`yatiml.dumps_json_function` instead:

.. code-block:: python

  # Create dumper
  dumps = yatiml.dumps_json_function(Submission)
Advanced features
=================

Class hierarchies
-----------------

One of the main features of object oriented design is inheritance. If your
objects can be categorised in classes and subclasses, then Python lets you code
them like that, and YAtiML can read and write them.

For example, let's add a description of the drawing to our Submission, in the
form of a list of the shapes that it consists of. We'll content ourselves with
a somewhat crude representation consisting of circles and squares.

.. literalinclude:: examples/class_hierarchy.py
  :caption: ``docs/examples/class_hierarchy.py``
  :language: python

Here, we have defined a class ``Shape``, and have added a list of Shapes as an
attribute to ``Submission``. Each shape has a location, its center, which is a
list of coordinates. Classes ``Circle`` and ``Square`` inherit from
``Shape``, and have some additional attributes. All the classe are passed when
creating the load function, and that's important, because only those classes
will be considered by YAtiML.

YAtiML will automatically recognize which subclass matches the object actually
specified in the list from the attributes that it has. If more than one subclass
matches, it will give an error message stating that the file being read is
ambiguous. If both a parent class and its child class match, YAtiML will prefer
the child class, and not consider it ambiguous.

Note that the child classes include the parent's class's ``center`` attribute in
their ``__init__``, and pass it on using ``super()``. This is required, as
otherwise YAtiML won't accept the ``center`` attribute for a subclass. Another
design option here would be to automatically merge the named attributes along
the inheritance path, and allow using a ``**kwargs`` on ``__init__`` to forward
additional attributes to the parent classes. The more explicit option is more
typing, but it also makes it easier to see what's going on when reading the
code, and that's very important for code maintainability. So that's what YAtiML
does.

Enumerations
------------

Enumerations, or enums, are types that are defined by listing a set of possible
values. In Python 3, they are made by creating a class that inherits from
``enum.Enum``. YAML does not have enumerations, but strings work fine provided
that you have something like YAtiML to check that the string that the user put
in actually matches one of the values of the enum type, and return the correct
value from the enum class. Here's how to add some colour to the drawings.

.. literalinclude:: examples/enums.py
  :caption: ``docs/examples/enums.py``
  :language: python

Note that the labels that YAtiML looks for are the names of the enum members,
not their values. In many existing standards, enums map to numerical values, or
if you're making something new, it's often convenient to use the values for
something else. The names are usually what you want to write though, so they're
probably easier for the users to write in the YAML file too. If you want
something else though, you can season your enumerations. See below for a general
explanation of seasoning, or look at :ref:`seasoning_enumerations` in the
Recipes section for some examples.

User-Defined Strings
--------------------

When defining file formats, you often find yourself with a need to define a
string with constraints. For example, Dutch postal codes consist of four digits,
followed by two uppercase letters. If you use a generic string type for a postal
code, you may end up accepting invalid values. A better solution is to define a
custom string type with a built-in constraint. In Python 3, you can do this
by deriving a class either from ``str`` or from ``collections.UserString``. The
latter is easier, so that's what we'll use in this example. Let's add the town
that our participant lives in to our YAML format, but insist that it be
capitalised.

.. literalinclude:: examples/user_defined_string.py
  :caption: ``docs/examples/user_defined_string.py``
  :language: python

Python's ``UserString`` provides an attribute ``data`` containing the actual
string, so all we need to do is test that for validity in the constructor. If
you spell the town using only lowercase letters, you'll get:

.. code-block:: none

  ValueError: Invalid TitleCaseString 'piedmont': Each word must start with a
  capital letter

Note that you can't make a TitleCaseString object containing 'piedmont' from
Python either, so the object model and the YAML format are consistent.

Python's UserString class tries very hard to look like a string by overloading
various special methods. Most of the time that's fine, but sometimes you have a
class that's really not much like a string on the Python side, but still should
be written to YAML as a string. In this case, you can add :class:`yatiml.String`
as a base class. YAtiML will then expect a string on the YAML side, call
``__init__`` with that string as the sole argument, and when dumping use
``str(obj)`` to obtain the string representation to write to the YAML file
(the result is then passed to ``_yatiml_sweeten()`` if you have it, so you can
still modify it if desired). Like classes derived from ``str`` and
``UserString``, such classes can be used as keys for dictionaries, but be sure
to implement ``__hash__()`` and ``__eq__()`` to make that work on the Python
side.

Seasoning your YAML
-------------------

For users who are manually typing YAML files, it is usually nice to have some
flexibility. For programmers processing the data read from such a file, it is
very convenient if everything is rigidly defined, so that they do not have to
take into account all sorts of corner cases. YAtiML helps you bridge this gap
with its support for seasoning.

In programming languages, small features that make the language easier to type,
but which do not add any real functionality are known as `syntactic sugar`. With
YAtiML, you can add a bit of extra processing halfway through the dumping
process to format your object in a nicer way. YAtiML calls this `sweetening`.
When loading, you can convert back to the single representation that matches
your class definition by `savourising`, savoury being the opposite of sweet.
Together, sweetening and savourising are referred to as `seasoning`.

Let's do another example. Having ages either as strings or as ints is not very
convenient if you want to check which age category to file a submission under.
So let's add a savourising function to convert strings to int on loading:

.. literalinclude:: examples/savorizing.py
  :caption: ``docs/examples/savorizing.py``
  :language: python

We have added a new ``_yatiml_savorize()`` class method to our Submission class.
This method will be called by YAtiML after the YAML text has been parsed, but
before our Submission object has been generated. This method is passed the
`node` representing the mapping that will become the object. The node is of
type :class:`yatiml.Node`, which in turn is a wrapper for an internal
ruamel.yaml object. Note that this method needs to be a classmethod, since
there is no object yet to call it with.

The :class:`yatiml.Node` class has a number of methods that you can use to
manipulate the node. In this case, we first check if there is an ``age``
attribute at all, and if so, whether it has a string as its value. This is
needed, because we are operating on the freshly-parsed YAML input, before any
type checks have taken place. In other words, that node may contain anything.
Next, we get the attribute's value, and then try to convert it to an int and set
it as the new value. If a string value was used that we do not know how to
convert, we raise a :class:`yatiml.SeasoningError`, which is the appropriate way
to signal an error during execution of ``_yatiml_savorize()``.

(At this point I should apologise for the language mix-up; the code uses
North-American spelling because it's rare to use British spelling in code and so
it would confuse everyone, while the documentation uses British spelling because
it's what its author is used to.)

When saving a Submission, we may want to apply the opposite transformation, and
convert some ints back to strings. That can be done with a ``_yatiml_sweeten``
classmethod:

.. literalinclude:: examples/sweetening.py
  :caption: ``docs/examples/sweetening.py``
  :language: python

The ``_yatiml_sweeten()`` method has the same signature as
``_yatiml_savorize()`` but is called when dumping rather than when loading. It
gives you access to the YAML node that has been produced from a Submission
object before it is written out to the YAML output. Here, we use the same
functions as before to convert some of the int values back to strings. Since we
converted all the strings to ints on loading above, we can assume that the value
is indeed an int, and we do not have to check.

Indeed, if we run this example, we get:

.. code-block:: none

  name: Youssou
  age: seven
  tool: pencils

However, there is still an issue. We have now used the seasoning functionality
of YAtiML to give the user the freedom to write ages either as words or as
numbers, while always giving the programmer ints to work with. However, the
programmer could still accidentally put a string into the age field when
constructing a Submission directly in the code, as the type annotation allows
it. This would then crash the ``_yatiml_sweeten()`` method when trying to dump
the object.

The solution, of course, is to change the type on the ``age`` attribute of
``__init__`` to ``int``. Unfortunately, this breaks loading. If you try to run
the savourise example above with ``age`` as type ``int``, then you will get

.. code-block:: none

  yatiml.exceptions.RecognitionError:   in "<unicode string>", line 1, column 1:
      name: Janice
      ^ (line: 1)
  Type mismatch, expected a Submission

The reason we get the error above is that by default, YAtiML recognises objects
of custom classes by their attributes, checking both names and types.
With the type of the ``age`` attribute now defined as ``int``, a mapping
containing an ``age`` with a string value is now no longer recognised as a
Submission object. A potential solution would be to apply seasoning before
trying to recognise, but to know how to savorise a mapping we need to know which
type it is or should be, and for that we need to recognise it. The way to fix
this is to override the default recognition function with our own, and make that
recognise both ``int`` and ``str`` values for ``age``.

Customising recognition
-----------------------

Customising the recognition function is done by adding a
``_yatiml_recognize()`` method to your class, like this:

.. literalinclude:: examples/custom_recognition.py
  :caption: ``docs/examples/custom_recognition.py``
  :language: python

This is again a classmethod, with a single argument of type
:class:`yatiml.UnknownNode` representing the node. Like
:class:`yatiml.Node`, :class:`yatiml.UnknownNode` wraps a YAML node, but this
class has helper functions intended for writing recognition functions.  Here,
we use :meth:`require_attribute` to list the required attributes and their
types. Since ``tool`` is optional, it is not required, and not listed. The
``age`` attribute is specified with the Union type we used before. Now, any
mapping that is in a place where we expect a Submission will be recognised as a
Submission, as long as it has a ``name`` attribute with a string value, and an
``age`` attribute that is either a string or an integer. If ``age`` is a
string, the ``_yatiml_savorize()`` method will convert it to an int, after
which a Submission object can be constructed without violating the type
constraint in the ``__init__()`` method.

In fact, the ``_yatiml_recognize()`` method here could be even simpler. In
every place in our document where a Submission can occur (namely the root),
only a Submission can occur. The Submission class does not have ancestors,
and it is never part of a Union. So there is never any doubt as to how to treat
the mapping, and in fact, the following will also work:

.. code-block:: python

  @classmethod
  def _yatiml_recognize(cls, node: yatiml.UnknownNode) -> None:
      pass

Now, if you try to read a document with, say, a float argument to ``age``, it
will be recognised as a Submission, the ``_yatiml_savorize()`` method will do
nothing with it, and you'll get an error message at the type check just before a
Submission is constructed.

This makes it clear that recognition is not a type check. Instead, its job is to
distinguish between different possible types in places where the class hierarchy
leaves some leeway to put in objects of different classes. If there is no such
leeway, the recognition stage does not need to do anything. If there is some
leeway, it just needs to do the minimum to exclude other possibilities.

However, since data models tend to evolve, it is usually a good idea to do a
full check anyway, so that if this class ends up being used in a Union, or if
you or someone else adds derived classes later, things will still work correctly
and there won't be any unnecessary ambiguity errors for the users.

Speaking of derived classes, note that while ``_yatiml_recognize()`` is
inherited by derived classes like any other Python method, YAtiML will only use
it for the class on which it is defined; derived classes will use automatic
recognition unless they have their own ``_yatiml_recognize()``. The same goes
for ``_yatiml_savorize()`` and `` _yatiml_sweeten()``.

Extra attributes
----------------

By default, YAtiML will match a mapping in a YAML file exactly: each required
attribute must be there, and any extraneous attributes give an error. However,
you may want to give your users the option of adding additional attributes. The
logical way for YAtiML to support this would be through having a ``**kwargs``
attribute to the ``__init__`` method, but unfortunately this would lose the
ordering information, since ``**kwargs`` is a plain unordered dict (although
this is in the process of changing in newer versions of Python). Also, there
wouldn't be an obvious way of saving such extra attributes again.

So, instead, extra attributes are sent to a ``_yatiml_extra`` parameter of type
``OrderedDict`` on ``__init__``, if there is one. You put this value into a
``_yatiml_extra`` attribute, whose contents YAtiML will then dump appended to
the normal attributes. If you want to be able to add extra attributes when
constructing an object using keyword arguments, then you can add a ``**kwargs``
parameter as well, and put any key-value pairs in it into ``self._yatiml_extra``
in your favourite order yourself.

Here is an example:

.. literalinclude:: examples/extra_attributes.py
  :caption: ``docs/examples/extra_attributes.py``
  :language: python

In this example, we use the ``tool`` attribute again, but with this code, we
could add any attribute, and it would show up in ``_yatiml_extra`` with no
errors generated.

Note that any explicit YAML tags on any mapping values of the extra attributes
or anywhere beneath them in the YAML tree will be stripped, so that this tree
will consist of plain lists and dicts. This is to avoid unexpected
user-controlled object construction, for safety reasons. These tags are
currently not added back on saving either, so it's good if the extra data does
not rely on them, better if it does not have any.

.. _hiding-attributes:

Hiding attributes
-----------------

By default, YAtiML assumes that your classes have a public attribute
corresponding to each parameter of their ``__init__`` method. If this
arrangement does not work for you, then you can override it by creating a
``_yatiml_attributes()`` method. This is `not` a classmethod, but an ordinary
method, because it is used for saving a particular instance of your class, to
which it needs access. If your custom class has a ``_yatiml_attributes()``
method defined, YAtiML will call that method instead of looking for public
attributes.  It should return an ``OrderedDict`` with names and values of the
attributes.

So far, we have been printing the values of public attributes to see the results
of our work. It would be better encapsulation to use private attributes instead,
with a ``__str__`` method to help printing. With ``_yatiml_attributes()``, that
can be done:

.. literalinclude:: examples/private_attributes.py
  :caption: ``docs/examples/private_attributes.py``
  :language: python

Further reading
---------------

You've reached the end of this tutorial, which means that you have seen all the
major features that YAtiML has. If you haven't already started, now is the time
to start making your awn YAML-based file format. You may want to have a look at
the :doc:`api`, and if you get stuck, there is the :doc:`problem_solving`
section to help you out.
.. YAtiML documentation master file, created by
   sphinx-quickstart on Thu Jun 21 11:07:11 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to YAtiML
=================

YAML-based file formats can be very handy, as YAML is easy to write by humans,
and parsing support for it is widely available. Just read your YAML file into a
document structure (a tree of nested dicts and lists), and manipulate that in
your code.

As long as that YAML file contains exactly what you expect, that works fine.
But if it contains a mistake, then you're likely to crash the program with a
cryptic error message, or worse (especially if the YAML file was loaded from the
Internet) it may do something unexpected.

To avoid that, you can validate your YAML using various schema checkers. You
write a description of what your YAML file must look like, then feed that to a
library which checks the incoming file against the description. That gives you a
better error message, but it's a lot of work.

YAtiML takes a different approach. Instead of a schema, you write a Python
class. You probably already know how to do that, so no need to learn anything.
YAtiML then generates loading and dumping functions for you, which convert
between YAML and Python objects. If needed, you can add some extra code to make
the YAML look nicer or implement special features.

YAtiML supports Python 3.6 and later.

If you use YAtiML for scientific work, we ask that you cite it. You can
`download a citation in various formats at the Research Software Directory
<https://www.research-software.nl/software/yatiml>`_.


Documentation Overview
======================

.. toctree::
   :maxdepth: 2

   why
   installation
   basic_tutorial
   advanced_features
   recipes
   problem_solving
   api


Development
===========

.. toctree::
  :maxdepth: 2

  development


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
API Reference
=============

Loading / Saving
----------------

.. autofunction:: yatiml.load_function
.. autofunction:: yatiml.dumps_function
.. autofunction:: yatiml.dumps_json_function
.. autofunction:: yatiml.dump_function
.. autofunction:: yatiml.dump_json_function

Seasoning
---------

.. autoclass:: yatiml.UnknownNode
    :no-special-members:

.. autoclass:: yatiml.Node
    :no-special-members:

Errors
------

.. autoclass:: yatiml.RecognitionError
.. autoclass:: yatiml.SeasoningError

Miscellaneous
-------------

.. autoclass:: yatiml.bool_union_fix
.. autoclass:: yatiml.logger
.. autoclass:: yatiml.String

Deprecated
----------

.. autoclass:: yatiml.Dumper
    :exclude-members: emit, emit_json

.. autofunction:: yatiml.add_to_dumper

.. autoclass:: yatiml.Loader
    :exclude-members: get_node, get_single_node

.. autofunction:: yatiml.add_to_loader
.. autofunction:: yatiml.set_document_type
