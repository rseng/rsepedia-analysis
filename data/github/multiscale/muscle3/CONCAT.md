############################
Contributing guidelines
############################

We welcome any kind of contribution to our software, from simple comment or question to a full fledged `pull request <https://help.github.com/articles/about-pull-requests/>`_. Please read and follow our `Code of Conduct <CODE_OF_CONDUCT.rst>`_.

A contribution can be one of the following cases:

#. you have a question;
#. you think you may have found a bug (including unexpected behavior);
#. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

You have a question
*******************

#. use the search functionality `here <https://github.com/yatiml/yatiml/issues>`_ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/yatiml/yatiml/issues>`_ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the `SHA hashcode <https://help.github.com/articles/autolinked-references-and-urls/#commit-shas>`_ of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
#. apply relevant labels to the newly created issue.

You want to make some kind of change to the code base
*****************************************************

#. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
#. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
#. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions `here <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`_ and `here <https://help.github.com/articles/syncing-a-fork/>`_);
#. make sure the existing tests still work by running ``python setup.py test``;
#. add your own tests (if necessary);
#. update or expand the documentation;
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the YAtiML repository on GitHub;
#. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`_.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
###########
Change Log
###########

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.

0.4.0
*****

Incompatible changes
--------------------

* `compute_elements` are now called `components` in .ymmsl files

Improved
--------

* Use latest OpenSSL library when installing it automatically

Fixed
-----

* Handling of non-contiguous and F-order numpy arrays
* C++ memory usage for large dicts/lists now more reasonable
* Improved shutdown when Python submodel crashes
* Logging warning message


0.3.2
*****

Improved
--------

* Accessing settings from C++ now more flexible
* Python produces more detailed logs to aid in debugging
* Improved pkg-config set-up
* Improved build system output to help find problems
* Documentation on logging in Python
* Protobuf dependency build now more compatible

Fixed
-----

* C++ list/dict building functions
* C++ use-after-free when receiving grids

Thanks
------

* Pavel for testing and reporting issues
* Dongwei for testing and reporting issues


0.3.1
*****

Added
-----

* Support for sending and receiving multidimensional grids/arrays
* Support for Python 3.8

Improved
--------

* Python 3.5.1 support
* Build compatibility on more operating systems

Thanks
------

* Olivier for testing, reporting and fixing build issues
* Pavel for testing and reporting build issues
* Hamid for testing and reporting build issues
* Ben for testing and reporting build issues


0.3.0
*****

Incompatible changes
--------------------

* Data::key() now returns std::string instead of DataConstRef.
* Data::value() now return Data rather than DataConstRef

Added
-----

* Support for Fortran, including MPI

Improved
--------

* Fixes to examples
* Small documentation improvements
* Improved compatibility with other packages using gRPC


Thanks
------

* Pavel for reporting documentation/examples issues
* Derek for testing on Eagle
* Dongwei for reporting the gRPC issue


0.2.0
*****

Added
-----

* Support for C++
* Support for MPI in C++

Improved
--------

* Cluster/HPC networking

Incompatible Changes
-------

* Fatal logic errors now throw instead of exiting, so that you have a chance
  to shut down the model cleanly before exiting.
* Instance.exit_error() was replaced by Instance.error_shutdown(), which no
  longer exits the process, it just shuts down the Instance.
* Central MUSCLE 3-managed settings are called settings everywhere now, not
  parameters. As a result, the API has changed in several places.


0.1.0
*****

Initial release of MUSCLE 3.

Added
-----
* Coupling different submodel instances
* Spatial and temporal scale separation and overlap
* Settings management
* Combining features
* Python support
* Initial distributed execution capability
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
.. image:: https://github.com/multiscale/muscle3/raw/develop/docs/source/muscle3_logo_readme.png
    :alt: MUSCLE 3

.. image:: https://readthedocs.org/projects/muscle3/badge/?version=develop
    :target: https://muscle3.readthedocs.io/en/develop/?badge=develop
    :alt: Documentation Build Status

.. image:: https://github.com/multiscale/muscle3/workflows/continuous_integration/badge.svg?branch=develop
    :target: https://github.com/multiscale/muscle3/actions
    :alt: Build Status

.. image:: https://api.codacy.com/project/badge/Coverage/ea0c833cf1ce4e13840c6498dfe27ff8
    :target: https://www.codacy.com/app/LourensVeen/muscle3
    :alt: Code Coverage

.. image:: https://api.codacy.com/project/badge/Grade/ea0c833cf1ce4e13840c6498dfe27ff8
    :target: https://www.codacy.com/app/LourensVeen/muscle3
    :alt: Codacy Grade

.. image:: https://requires.io/github/multiscale/muscle3/requirements.svg?branch=develop
     :target: https://requires.io/github/multiscale/muscle3/requirements/?branch=develop
     :alt: Requirements Status

.. image:: https://zenodo.org/badge/122876985.svg
   :target: https://zenodo.org/badge/latestdoi/122876985

.. image:: https://img.shields.io/badge/rsd-muscle3-00a3e3.svg
   :target: https://www.research-software.nl/software/muscle3

|

MUSCLE 3 is the third incarnation of the MUSCLE Multiscale Coupling Library and
Environment. It is developed by the e-MUSC project of the University of
Amsterdam and the Netherlands eScience Center.

Browse to `the MUSCLE 3 documentation`_ to get started.


Legal
=====

MUSCLE 3 is Copyright 2018-2020 University of Amsterdam and Netherlands eScience
Center. It is licensed under the Apache License 2.0.


Contributing
============

We welcome contributions to MUSCLE 3! For small additions or fixes, please
submit a pull request, and we'll review your changes and apply them. If you
want to add larger features, please make an issue and describe what you want to
do, so that we can coordinate and avoid double work. For more details, see the
`documentation on contributing`_.

The MUSCLE 3 community fosters a welcoming atmosphere, and participation is open
to everyone. In order to ensure that it stays that way, we have a
`code of conduct`_, and we will enforce it if necessary. Please treat your
fellow human beings with respect.


.. _`the MUSCLE 3 documentation`: https://muscle3.readthedocs.io
.. _`documentation on contributing`: http://muscle3.readthedocs.io/en/latest/contributing.html
.. _`code of conduct`: https://github.com/multiscale/muscle3/blob/develop/CODE_OF_CONDUCT.md
========================
Contributing to MUSCLE 3
========================

MUSCLE 3 is developed by a community of contributors, and we welcome anyone and
everyone to help us build it. Maybe you've found a mistake in the code (it
happens, we're only human!) that you would like to fix, want to add a new
feature, or improve the documentation. That's great!

However, because MUSCLE 3 is worked on by different people simultaneously, we do
need some system to keep track of all the changes, so that it doesn't become an
unintelligible mess. We use the Git version control system for this, together
with the Github hosting platform.

Describing all of Git and Github here would be too much (both have their own
documentation), but we do describe the process of contributing a change to
MUSCLE 3 in some detail. We're assuming that you already have Git installed, and
have a Github account that you're logged in with. Once that's done, please
proceed as below.


Make an Issue
=============

Issues are found in a tab at the top of the repository home page. Please check
to see that the bug you want to fix or the feature you want to add does not
already have an issue dedicated to it. If it does, feel free to add to the
discussion. If not, please make a new issue.

If you have discovered an error in MUSCLE 3, please describe

* What you were trying to achieve
* What you did to achieve this
* What result you expected
* What happened instead

With that information, it should be possible for someone else to reproduce the
problem, or to check that a proposed fix has really solved the issue.

If you would like a new feature (or would like to add it yourself), please
describe

* The context in which you're using MUSCLE 3
* What you are trying to achieve
* How the new feature would help you do that
* What the new feature should do and how it should work

If you want to fix the bug or implement the feature yourself, you'll have to set
up a development environment.


Get a local repository
======================

First, you'll need a public Git repository that you have write access to, and which
contains MUSCLE 3. If you have write access to the main MUSCLE 3 repository,
then you're set already. If not, click the Fork button at the top right of the
main MUSCLE 3 repository page. This will make a copy of the repository that you
can write to, under your own account.

Next, you need a copy of this public repository on your local computer, so that
you can easily edit files. This is done by cloning the repository::

  git clone git@github.com:/multiscale/muscle3.git

for the original repository, of if you've forked it, probably something like::

  git clone git@github.com:/UserName/muscle3.git

You can use the ``Clone or download`` button to copy the location into your
clipboard.

This will create a local directory named ``muscle3`` containing a Git
repository that you can work in.


Install the tools
=================

Inside this directory, first check out the ``develop`` branch::

  cd muscle3
  git checkout develop

Next, you'll need a virtual environment to install the development tool in::

  virtualenv -p python3 ~/Envs/muscle3

The path is arbitrary, a virtualenv is just a directory with software in it, so
you can put it anywhere you want. An ``Envs`` subdirectory in your home
directory is somewhat standard though.  Do note that once created, a virtualenv
cannot be moved around without breaking (but it's easy to make a new one in a
new location instead).

Then, activate the environment, and install the development tools
into it::

  source ~/Envs/muscle3/bin/activate
  pip install -e .[dev]

and run the test suite, just to make sure that you have everything in working
order::

  make test


Make changes
============

Changes should be made on a fresh branch, dedicated to the issue you are
addressing. That way, multiple people can work at the same time, and multiple
issues can be addressed simultaneously, without everyone getting confused. So to
start, make a new branch, starting from the `develop` branch::

  git checkout develop
  git checkout -b issue-123

where, instead of 123, you put the number of the issue you're addressing.

For a simple bug, the whole issue can probably be fixed with a single change. (A
single change may affect different functions or files, but those changes should
be related, so that it doesn't make sense to do one and not the other.)

First, make a test that detects the bug you're going to fix, and make sure it
fails by running the test suite::

  make test

Then, fix the bug, and run the tests again, now to make sure that they pass.
When everything is in order, you can commit the change::

  git add changed_file.py some_dir/other_changed_file.py
  git commit -m "Made x configurable"
  git push

The last command pushes your changes back to the server, so that other
developers can see them. Do this after every commit! It makes it much easier to
collaborate.

Note that commits should be as small as possible. Please do not hack away for
several days and then commit a whole giant bunch of changes in a single commit,
as that makes it impossible to figure out later what was changed when, and which
change introduced that bug we are trying to find.

If you are solving a more complex problem, then a single commit will not be
enough, and you'll have to make a series of them, repeating the above steps,
until the issue is solved. Starting with a test is often the best way of going
about adding a new feature as well. You'll find that you'll need to think about
what your new feature should do and how it should work to create the test(s),
and once you've done that implementing it is a lot easier!

One last note: **Never copy-paste code from another program!**. It's fine to
have external dependencies (although we do try to limit them, to try to keep
installation simple), but those should be kept separate. Copy-pasting code leads
to complicated legal issues that we would really like to avoid. So please, only
contribute code that you wrote yourself. Thanks!


Make a pull request
===================

Once you've made all the changes needed to resolve the issue, the next step is
to make a pull request. Your changes so far are on a branch, either in the main
repository, or in a fork of the main repository. A pull request is a request to
the maintainers of MUSCLE 3 to take the changes on your branch, and incorporate
them into the main version of the software.

To make a pull request, make sure that you have committed and pushed all your
changes, and that the tests pass. Then, go to the Github homepage of your fork,
if you have one, or the main MUSCLE 3 repository. If you've just pushed, then
Github will show a "Compare & pull request" button. Otherwise, look up your
branch using the top left drop-down button, and then click the "New pull
request" button next to it.

This gives you a page describing your pull request. You will want to request a
merge from your issue branch, to the develop branch in the main MUSCLE 3
repository. Add a description of the changes you've made, and click "Create pull
request", and you're all set.


Interact
========

Like issues, pull requests on Github are a kind of discussion forum, in which
the proposed changes can be discussed. We may ask you to make some improvements
before we accept your pull request. While the pull request is open, any
additional commits pushed to your public branch will automatically show up
there.

Once we're all satisfied with the change, the pull request will be accepted, and
your code will become part of MUSCLE 3. Thank you!
.. _api-docs-fortran:

API Documentation for Fortran
=============================

This page provides full documentation for the Fortran API of MUSCLE 3.

A note on types
---------------

Fortran variables have a type and a kind. The type specifies a broad category,
i.e. whether it is a logical value, text, or an integer or floating point
number. The kind further subdivides the numbers into kinds with different
precision or size. Traditional Fortran does not standardise which kinds are
available, and so this may vary from compiler to compiler and from machine to
machine. When Fortran was invented, number representations were much less
standardised than today, and so the language designers left it to the
implementation to avoid making the language impossible to implement on
particular hardware.

Sincle the libmuscle API uses integers and floating point numbers all over the
place, we need to be able to declare these kinds. Also, the type and/or kind
sometimes shows up in a function name, so we need some consistent naming to
avoid confusion.

This is complicated somewhat by the lack of guarantees on which types will be
available. However, these days, all hardware that is relevant to MUSCLE 3 uses
two's complement integers with sizes of 8, 16, 32 and 64 bits, and IEEE 754
floating point numbers using 32 and 64-bit precision. So we can standardise on
these without excluding anyone.

In many compilers, integer types of size n bytes can be written ``integer*n``.
This is not standard however, so it's not guaranteed to work. In Fortran 2003,
a standardised way of selecting a kind with at least a given precision became
available in the form of the ``selected_int_kind`` and ``selected_real_kind``
functions. These work, but they're a pain to write out all the time.

So, for each of the above mentioned integer and floating point types, libmuscle
has a unique name, which is used in function names, and a parameter using the
same name which defines the corresponding kind. You can use these to declare
variables which you use with the libmuscle API, or you can write the same
type/kind in another way that your compiler supports and it should work as well.
The names and types are as follows:

+---------+-----------------+----------------------------------------+
|  Name   | Kind parameter  | Description                            |
+=========+=================+========================================+
|  int1   | LIBMUSCLE_int1  | 8-bit signed 2's complement integer.   |
+---------+-----------------+----------------------------------------+
|  int2   | LIBMUSCLE_int2  | 16-bit signed 2's complement integer.  |
+---------+-----------------+----------------------------------------+
|  int4   | LIBMUSCLE_int4  | 32-bit signed 2's complement integer.  |
+---------+-----------------+----------------------------------------+
|  int8   | LIBMUSCLE_int8  | 64-bit signed 2's complement integer.  |
+---------+-----------------+----------------------------------------+
| integer |                 | Default integer.                       |
+---------+-----------------+----------------------------------------+
|  size   | LIBMUSCLE_size  | Used for array indiced and sizes.      |
+---------+-----------------+----------------------------------------+
|  real4  | LIBMUSCLE_real4 | 32-bit IEEE 754 floating point.        |
+---------+-----------------+----------------------------------------+
|  real8  | LIBMUSCLE_real8 | 64-bit IEEE 754 floating point.        |
+---------+-----------------+----------------------------------------+

As said above, if you are used to ``integer*n`` and ``real*n`` style types, then
you can use those, using e.g. ``integer*2`` where the API expects a
``LIBMUSCLE_int2``. For ``LIBMUSCLE_size``, ``integer*8`` will probably work,
otherwise try ``integer*4``.


Namespace LIBMUSCLE
-------------------

.. f:module:: libmuscle
    :synopsis: Fortran module for libmuscle

.. f:currentmodule:: libmuscle

LIBMUSCLE_Data
``````````````
.. f:type:: LIBMUSCLE_Data

    Represents a libmuscle Data object. This is an opaque object that may be
    returned from and passed to libmuscle functions, but does not contain any
    directly accessible members.

    With respect to creation, assignment and copying, :f:type:`LIBMUSCLE_Data`
    objects act like Python objects: basic values (logicals, strings, integers
    and reals) get copied, while for lists, dictionaries, and byte arrays, the
    variable contains a reference which gets copied.

.. f:type:: LIBMUSCLE_DataConstRef

    Represents a read-only reference to a libmuscle Data object. This is an
    opaque object like :f:type:`LIBMUSCLE_Data`. DataConstRef objects work
    exactly the same as Data objects, except that the names of the corresponding
    functions start with ``LIBMUSCLE_DataConstRef``, and that only creation,
    freeing, and non-modifying operations are supported. Since these functions
    are otherwise identical, they are not documented separately here. If you
    want to know how to use say ``LIBMUSCLE_DataConstRef_as_int8``, look
    up :f:func:`LIBMUSCLE_Data_as_int8`.

    This class is mainly used to represent received messages, which should not
    be modified.

.. f:function:: LIBMUSCLE_Data_create()

    Creates a Data object representing nil.

    Nil is a special "no data" value, like nullptr in C++ or None in Python.

    :r obj: The new Data object
    :rtype obj: LIBMUSCLE_Data

.. f:function:: LIBMUSCLE_Data_create(value)

    Creates a Data object representing a logical value.

    :p logical value: The value to represent.
    :r obj: The new Data object
    :rtype obj: LIBMUSCLE_Data

.. f:function:: LIBMUSCLE_Data_create(value)

    Creates a Data object representing a character value.

    :p character value: The value to represent.
    :r obj: The new Data object
    :rtype obj: LIBMUSCLE_Data

.. f:function:: LIBMUSCLE_Data_create(value)

    Creates a Data object representing an integer value. Supported kinds are
    the default ``integer`` kind, and ``LIBMUSCLE_intN`` with N set to 1, 2, 4
    or 8. See the note at the top of this page for more on integer types in
    libmuscle.

    Note that while libmuscle supports unsigned integers, these don't exist
    in Fortran. They will be mapped to the corresponding signed type, which
    may cause a silent overflow if the number is out of the signed type's
    range.

    :p integer value: The value to represent (see above).
    :r obj: The new Data object
    :rtype obj: LIBMUSCLE_Data

.. f:function:: LIBMUSCLE_Data_create(value)

    Creates a Data object representing a real value. Supported kinds are
    ``LIBMUSCLE_real4`` and ``LIBMUSCLE_real8``. See the note at the top of this
    page for more on real types in libmuscle.

    :p real value: The value to represent.
    :r obj: The new Data object
    :rtype obj: LIBMUSCLE_Data

.. f:function:: LIBMUSCLE_Data_create(value)

    Creates a Data object containing a :f:type:`YMMSL_Settings` object.

    :p YMMSL_Settings value: The Settings value to represent.
    :r obj: The new Data object.
    :rtype obj: LIBMUSCLE_Data

.. f:function:: LIBMUSCLE_Data_create_grid(data_array)

    Creates a Data object containing a grid (array).

    The argument must be an array of type ``logical`` and default kind,
    an array of type ``integer`` and kind ``LIBMUSCLE_int4`` or
    ``LIBMUSCLE_int8``, or an array of type ``real`` and kind
    ``LIBMUSCLE_real4`` or kind ``LIBMUSCLE_real8``.

    Grids created with this function have no index names.

    :p array data_array: The array of data to represent.
    :r obj: The new Data object.
    :rtype obj: LIBMUSCLE_Data

.. f:function:: LIBMUSCLE_Data_create_grid(data_array, index_name, ...)

    Creates a Data object containing a grid (array).

    The ``data_array`` argument must be an array of type ``logical`` and
    default kind, an array of type ``integer`` and kind ``LIBMUSCLE_int4``
    or ``LIBMUSCLE_int8``, or an array of type ``real`` and kind
    ``LIBMUSCLE_real4`` or kind ``LIBMUSCLE_real8``.

    If an ``n``-dimensional array is passed as the first argument, then
    there must be ``n`` additional arguments of type ``character``, giving
    the names of the indexes in order. For instance, if your 2D array
    represents a table and you index it ``data_array(row, column)`` then
    ``"row"`` and ``"column"`` would be reasonable index names here. Note
    that MUSCLE 3 does not use these names, they are here to make it
    easier to understand the message on the receiver side, or if it is
    saved and analysed later.

    :p array data_array: The array of data to represent.
    :p character index_name: The names of the grid's indexes.
    :r obj: The new Data object.
    :rtype obj: LIBMUSCLE_Data

.. f:function:: LIBMUSCLE_Data_create_dict()

    Creates a Data object containing an empty dictionary.

    :r obj: The new Data object
    :rtype obj: LIBMUSCLE_Data

.. f:function:: LIBMUSCLE_Data_create_list()

    Creates a Data object containing an empty list.

    :r obj: The new Data object
    :rtype obj: LIBMUSCLE_Data

.. f:function:: LIBMUSCLE_Data_create_nils(size)

    Creates a Data object containing a list of ``size`` nil values.

    :p integer size: The number of nil values to put into the list
        (kind=LIBMUSCLE_size).
    :r obj: The new Data object
    :rtype obj: LIBMUSCLE_Data

.. f:function:: LIBMUSCLE_Data_create_byte_array(size)

    Creates a Data object containing a byte array of the given number of bytes.

    :p integer size: The number of bytes to allocate for the array
            (LIBMUSCLE_size).
    :r obj: The new Data object.
    :rtype obj: LIBMUSCLE_Data

.. f:function:: LIBMUSCLE_Data_create_byte_array(buf)

    Creates a Data object referring to the given data.

    The buffer passed will not be copied! This creates a Data object that refers
    to your buffer, and you need to make sure that that buffer exists for as
    long as the Data object (and/or any copies of it) is used.

    Example:

    .. code-block:: fortran

        type(LIBMUSCLE_Data) :: d1
        character(len=1), dimension(1024) :: bytes
        character(len=1), dimension(:), allocatable :: buf

        ! Create some data
        do i = 1, 1024
            bytes(i) = achar(mod(i, 256))
        end do

        ! Create a Data object referring to it
        d1 = LIBMUSCLE_Data_create_byte_array(bytes)
        ! Now d1 contains a byte array of size 1024

        ! Extract the data again, into a new buffer
        allocate(buf(LIBMUSCLE_Data_size(d1)))
        call LIBMUSCLE_Data_as_byte_array(d1, buf)
        ! Now, ichar(buf(i)) equals mod(i, 256)

        ! Clean up the buffer and the Data object
        deallocate(buf)
        call LIBMUSCLE_Data_free(d1)

    :p character buf: An array of characters to refer to.
    :r obj: The new Data object
    :rtype obj: LIBMUSCLE_Data

.. f:subroutine:: LIBMUSCLE_Data_free(self)

    Frees a Data object.

    This frees the resources associated with the given Data object. Do not use
    the object for anything after calling this, because it will be invalid.

    :p LIBMUSCLE_Data self: The Data object to free.

.. f:subroutine:: LIBMUSCLE_Data_set(self, value)

    Assigns the value of Data object ``value`` to Data object ``self``. Both
    ``value`` and ``self`` must have a value (nil is okay, but an uninitialised
    :f:type:`LIBMUSCLE_Data` is not). If ``value`` holds a basic type, then the
    value will be copied into ``self``, overwriting any previous value in
    ``self``. If ``value`` holds a list, dict, or byte array, then ``self`` will
    end up referring to the same object as ``value``.

    If you haven't created ``self`` yet, then it's shorter to use
    :f:func:`LIBMUSCLE_Data_create`.

    This is the equivalent of ``self = value`` in C++ or Python.

    :p LIBMUSCLE_Data self: The Data object to assign to.
    :p LIBMUSCLE_Data value: The Data object to assign from.

.. f:subroutine:: LIBMUSCLE_Data_set(self, value)

    Assigns the value of ``value`` to ``self``. ``self`` must be an initialised
    :f:type:`LIBMUSCLE_Data` object. If your target :f:type:`LIBMUSCLE_Data`
    object does not exist yet, use :f:func:`LIBMUSCLE_Data_create` instead.

    This is the equivalent of ``self = value`` in C++ or Python.

    Value may be of types ``logical``, ``character``, ``integer`` or ``real``.
    Integer kinds may be those representing 8-bit, 16-bit, 32-bit and 64-bit
    values, real kinds may be 32-bit single and 64-bit double precision.

    :p LIBMUSCLE_Data self: The Data object to assign to.
    :p ``see_above`` value: The value to assign from.

.. f:function:: LIBMUSCLE_Data_is_a_logical(self)

    Determine whether the Data object contains a logical (boolean) value.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains a logical value.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_character(self)

    Determine whether the Data object contains a character (string) value.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains a character value.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_int(self)

    Determine whether the Data object contains an integer value.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains an integer value.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_int1(self)

    Determine whether the Data object contains an integer value.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains an integer value.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_int2(self)

    Determine whether the Data object contains an integer value.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains an integer value.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_int4(self)

    Determine whether the Data object contains an integer value.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains an integer value.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_real4(self)

    Determine whether the Data object contains a single precision floating
    point value.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains a single precision float value.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_real8(self)

    Determine whether the Data object contains a double precision floating
    point value.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains a double precision float value.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_settings(self)

    Determine whether the Data object contains a Settings value.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains a Settings value.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_nil(self)

    Determine whether the Data object contains a nil value.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains a nil value.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_dict(self)

    Determine whether the Data object contains a dictionary value.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains a dictionary.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_list(self)

    Determine whether the Data object contains a list value.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains a list.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_grid_of_logical(self)

    Determine whether the Data object contains a grid of logical values.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains a grid of logical.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_grid_of_int4(self)

    Determine whether the Data object contains a grid of integer values
    of kind ``LIBMUSCLE_int4``.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains a grid of int4.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_grid_of_int8(self)

    Determine whether the Data object contains a grid of integer values
    of kind ``LIBMUSCLE_int8``.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains a grid of int8.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_grid_of_real4(self)

    Determine whether the Data object contains a grid of real values
    of kind ``LIBMUSCLE_real4``.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains a grid of real4.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_grid_of_real8(self)

    Determine whether the Data object contains a grid of real8 values
    of kind ``LIBMUSCLE_real8``.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains a grid of real8.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_is_a_byte_array(self)

    Determine whether the Data object contains a byte array value.

    :p LIBMUSCLE_Data self: The Data object to inspect.
    :r is: True if the object contains a byte array.
    :rtype is: logical

.. f:function:: LIBMUSCLE_Data_size(self)

    Returns the size of a list (number of items), dict (number of key/value
    pairs), or byte array (number of bytes).

    :p LIBMUSCLE_Data self: The Data object to get the size of.
    :r size: The size of the object.
    :rtype size: integer (kind=LIBMUSCLE_size)

.. f:function:: LIBMUSCLE_Data_as_logical(self, err_code, err_msg)

    Access a logical value.

    You can use :f:func:`LIBMUSCLE_Data_is_a_logical` to ascertain that the Data
    object contains a logical value.

    If the Data object does not contain a logical (boolean) value, then an error
    message will be printed and execution will be halted.

    Alternatively, you can pass an argument for ``err_code``, or for both
    ``err_code`` and ``err_msg``, to catch the error.

    If ``err_code`` equals ``LIBMUSCLE_success`` after the call, then the
    returned value is the logical value held in this Data object. If it equals
    ``LIBMUSCLE_runtime_error``, the Data value did not contain a logical
    (boolean) value. If you passed an ``err_msg`` argument as well, then the
    passed variable will contain an appropriate error message in case of error,
    and needs to be deallocated (using ``deallocate()``) when you're done with
    it.

    Example:

    .. code-block:: fortran

        type(LIBMUSCLE_Data) :: mydata
        integer :: err_code
        logical :: val
        character(len=:), allocatable :: str, err_msg

        ! Create data object containing a logical value
        mydata = LIBMUSCLE_Data_create(.true.)
        ! Retrieve the value
        val = LIBMUSCLE_Data_as_logical(mydata)
        ! val equals .true. here
        ! Attempt to (incorrectly) retrieve a character value
        str = LIBMUSCLE_Data_as_character(mydata, err_code, err_msg)
        if (err_code .ne. LIBMUSCLE_success) then
            print *, err_msg
            ! Need to free the memory if an error message was returned
            deallocate(err_msg)
        end if
        ! Free the data object
        call LIBMUSCLE_Data_free(mydata)

    :p LIBMUSCLE_Data self: The Data object to get a logical value out of.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value.
    :rtype value: logical

.. f:function:: LIBMUSCLE_Data_as_character(self, err_code, err_msg)

    Access a character (string) value.

    You can use :f:func:`LIBMUSCLE_Data_is_a_character` to ascertain that the
    Data object contains a character value.

    If the Data object does not contain a character (string) value, then an
    error message will be printed and execution will be halted.

    Alternatively, you can pass an argument for ``err_code``, or for both
    ``err_code`` and ``err_msg``, to catch the error.

    If ``err_code`` equals ``LIBMUSCLE_success`` after the call, then the
    returned value is the character value held in this Data object. If it equals
    ``LIBMUSCLE_runtime_error``, the Data value did not contain a character
    (string) value. If you passed an ``err_msg`` argument as well, then the
    passed variable will contain an appropriate error message in case of error,
    and needs to be deallocated (using ``deallocate()``) when you're done with
    it.

    Note that the result variable will be allocated (unless an error occurs),
    and must be deallocated when you're done with the resulting string, or
    you'll have a memory leak.

    Example:

    .. code-block:: fortran

        type(LIBMUSCLE_Data) :: mydata
        character(len=:), allocatable :: str

        ! Create a data object containing a character value
        mydata = LIBMUSCLE_Data_create('Example')
        ! Retrieve the value
        str = LIBMUSCLE_Data_as_character(mydata)
        ! Free the retrieved copy of the character
        deallocate(str)
        ! Free the data object
        call LIBMUSCLE_Data_free(mydat

    See :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    :p LIBMUSCLE_Data self The Data object to get a character out of.
    :p integer err_code: An error code output (optional)
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value.
    :rtype value: character(len=:), allocatable

.. f:function:: LIBMUSCLE_Data_as_int(self, err_code, err_msg)

    Access an integer value of default kind.

    You can use :f:func:`LIBMUSCLE_Data_is_a_int` to ascertain that the Data
    object contains a default integer value.

    If the Data object does not contain an integer value, then an error
    message will be printed and execution will be halted.

    Alternatively, you can pass an argument for ``err_code``, or for both
    ``err_code`` and ``err_msg``, to catch the error.

    If ``err_code`` equals ``LIBMUSCLE_success`` after the call, then the
    returned value is the integer value held in this Data object. If it equals
    ``LIBMUSCLE_runtime_error``, the Data value did not contain an integer
    value. If you passed an ``err_msg`` argument as well, then the passed
    variable will contain an appropriate error message in case of error, and
    needs to be deallocated (using ``deallocate()``) when you're done with it.

    Example:

    .. code-block:: fortran

        type(LIBMUSCLE_Data) :: mydata
        integer :: number

        ! Create a data object containing an integer value
        mydata = LIBMUSCLE_Data_create(42424242)
        ! Retrieve the value
        number = LIBMUSCLE_Data_as_int(mydata)
        ! Free the data object
        call LIBMUSCLE_Data_free(mydata)

    See :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    :p LIBMUSCLE_Data self The Data object to get an integer value out of.
    :p integer err_code: An error code output (optional)
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value.
    :rtype value: integer

.. f:function:: LIBMUSCLE_Data_as_int1(self, err_code, err_msg)

    Access an int value that fits in 8 bits.

    You can use :f:func:`LIBMUSCLE_Data_is_a_int1` to ascertain that the Data
    object contains an ``int1`` value.

    If the Data object does not contain a ``LIBMUSCLE_int1`` (integer with
    ``selected_int_kind(2)``) value, then an error message will be printed and
    execution will be halted.

    Alternatively, you can pass an argument for ``err_code``, or for both
    ``err_code`` and ``err_msg``, to catch the error.

    If ``err_code`` equals ``LIBMUSCLE_success`` after the call, then the
    returned value is the integer value held in this Data object. If it equals
    ``LIBMUSCLE_runtime_error``, the Data value did not contain an integer
    value. If you passed an ``err_msg`` argument as well, then the passed
    variable will contain an appropriate error message in case of error, and
    needs to be deallocated (using ``deallocate()``) when you're done with it.

    Example:

    .. code-block:: fortran

        type(LIBMUSCLE_Data) :: mydata
        integer(kind=LIBMUSCLE_int1) :: number

        ! Create a data object containing an integer value
        mydata = LIBMUSCLE_Data_create(42_LIBMUSCLE_int1)
        ! Retrieve the value
        number = LIBMUSCLE_Data_as_int1(mydata)
        ! Free the data object
        call LIBMUSCLE_Data_free(mydata)

    See :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    :p LIBMUSCLE_Data self The Data object to get an integer value out of.
    :p integer err_code: An error code output (optional)
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value.
    :rtype value: integer(kind=LIBMUSCLE_int1)

.. f:function:: LIBMUSCLE_Data_as_int2(self, err_code, err_msg)

    Access an int value that fits in 16 bits.

    You can use :f:func:`LIBMUSCLE_Data_is_a_int2` to ascertain that the Data
    object contains an integer value.

    If the Data object does not contain a ``LIBMUSCLE_int2`` (integer with
    ``selected_int_kind(4)``) value, then an error message will be printed and
    execution will be halted.

    Alternatively, you can pass an argument for ``err_code``, or for both
    ``err_code`` and ``err_msg``, to catch the error.

    If ``err_code`` equals ``LIBMUSCLE_success`` after the call, then the
    returned value is the integer value held in this Data object. If it equals
    ``LIBMUSCLE_runtime_error``, the Data value did not contain an integer
    value. If you passed an ``err_msg`` argument as well, then the passed
    variable will contain an appropriate error message in case of error, and
    needs to be deallocated (using ``deallocate()``) when you're done with it.

    Example:

    .. code-block:: fortran

        type(LIBMUSCLE_Data) :: mydata
        integer(kind=LIBMUSCLE_int2) :: number

        ! Create a data object containing an integer value
        mydata = LIBMUSCLE_Data_create(4242_LIBMUSCLE_int2)
        ! Retrieve the value
        number = LIBMUSCLE_Data_as_int2(mydata)
        ! Free the data object
        call LIBMUSCLE_Data_free(mydata)

    See :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    :p LIBMUSCLE_Data self The Data object to get an integer value out of.
    :p integer err_code: An error code output (optional)
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value.
    :rtype value: integer(kind=LIBMUSCLE_int2)

.. f:function:: LIBMUSCLE_Data_as_int4(self, err_code, err_msg)

    Access an integer value.

    You can use :f:func:`LIBMUSCLE_Data_is_a_int4` to ascertain that the Data
    object contains an integer value of kind ``LIBMUSCLE_int4``.

    If the Data object does not contain an int (integer) value, then an error
    message will be printed and execution will be halted.

    Alternatively, you can pass an argument for ``err_code``, or for both
    ``err_code`` and ``err_msg``, to catch the error.

    If ``err_code`` equals ``LIBMUSCLE_success`` after the call, then the
    returned value is the integer value held in this Data object. If it equals
    ``LIBMUSCLE_runtime_error``, the Data value did not contain an integer
    value. If you passed an ``err_msg`` argument as well, then the passed
    variable will contain an appropriate error message in case of error, and
    needs to be deallocated (using ``deallocate()``) when you're done with it.

    Example:

    .. code-block:: fortran

        type(LIBMUSCLE_Data) :: mydata
        integer(LIBMUSCLE_int4) :: number

        ! Create a data object containing an integer value
        mydata = LIBMUSCLE_Data_create(42424242_LIBMUSCLE_int4)
        ! Retrieve the value
        number = LIBMUSCLE_Data_as_int4(mydata)
        ! Free the data object
        call LIBMUSCLE_Data_free(mydata)

    See :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    :p LIBMUSCLE_Data self The Data object to get an integer value out of.
    :p integer err_code: An error code output (optional)
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value.
    :rtype value: integer(LIBMUSCLE_int4)

.. f:function:: LIBMUSCLE_Data_as_int8(self, err_code, err_msg)

    Access an integer value of kind ``LIBMUSCLE_int8``..

    You can use :f:func:`LIBMUSCLE_Data_is_a_int8` to ascertain that the Data
    object contains a 64-bit integer value.

    If the Data object does not contain a ``LIBMUSCLE_int8`` (integer with
    ``selected_int_kind(18)``) value, then an error message will be printed and
    execution will be halted.

    Alternatively, you can pass an argument for ``err_code``, or for both
    ``err_code`` and ``err_msg``, to catch the error.

    If ``err_code`` equals ``LIBMUSCLE_success`` after the call, then the
    returned value is the integer value held in this Data object. If it equals
    ``LIBMUSCLE_runtime_error``, the Data value did not contain an integer
    value. If you passed an ``err_msg`` argument as well, then the passed
    variable will contain an appropriate error message in case of error, and
    needs to be deallocated (using ``deallocate()``) when you're done with it.

    Example:

    .. code-block:: fortran

        type(LIBMUSCLE_Data) :: mydata
        integer(kind=LIBMUSCLE_int8) :: number

        ! Create a data object containing an integer value
        mydata = LIBMUSCLE_Data_create(123456789123456789_LIBMUSCLE_int8)
        ! Retrieve the value
        number = LIBMUSCLE_Data_as_int8(mydata)
        ! Free the data object
        call LIBMUSCLE_Data_free(mydata)

    See :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    :p LIBMUSCLE_Data self: The Data object to get an integer value out of.
    :p integer err_code: An error code output (optional)
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value.
    :rtype value: integer(kind=LIBMUSCLE_int8)

.. f:function:: LIBMUSCLE_Data_as_real4(self, err_code, err_msg)

    Access a single-precision (4 byte) real value.

    You can use :f:func:`LIBMUSCLE_Data_is_a_real4` to ascertain that the Data
    object contains a single-precision real value.

    If the Data object does not contain a ``LIBMUSCLE_real4`` (real with
    ``selected_real_kind(6)``) value, then an error message will be printed and
    execution will be halted.

    Alternatively, you can pass an argument for ``err_code``, or for both
    ``err_code`` and ``err_msg``, to catch the error.

    If ``err_code`` equals ``LIBMUSCLE_success`` after the call, then the
    returned value is the real value held in this Data object. If it equals
    ``LIBMUSCLE_runtime_error``, the Data value did not contain a real
    value. If you passed an ``err_msg`` argument as well, then the passed
    variable will contain an appropriate error message in case of error, and
    needs to be deallocated (using ``deallocate()``) when you're done with it.

    Example:

    .. code-block:: fortran

        type(LIBMUSCLE_Data) :: mydata
        real(kind=LIBMUSCLE_real4) :: number

        ! Create a data object containing a real value
        mydata = LIBMUSCLE_Data_create(42.0)
        ! Retrieve the value
        number = LIBMUSCLE_Data_as_real4(mydata)
        ! Free the data object
        call LIBMUSCLE_Data_free(mydata)

    See :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    :p LIBMUSCLE_Data self: The Data object to get a single-precision real
            value out of.
    :p integer err_code: An error code output (optional)
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value.
    :rtype value: real(kind=LIBMUSCLE_real4)

.. f:function:: LIBMUSCLE_Data_as_real8(self, err_code, err_msg)

    Access a double-precision (8 byte) real value.

    You can use :f:func:`LIBMUSCLE_Data_is_a_real8` to ascertain that the Data
    object contains a double-precision real value.

    If the Data object does not contain a ``LIBMUSCLE_real8`` (real with
    ``selected_real_kind(15)``) value, then an error message will be printed and
    execution will be halted.

    Alternatively, you can pass an argument for ``err_code``, or for both
    ``err_code`` and ``err_msg``, to catch the error.

    If ``err_code`` equals ``LIBMUSCLE_success`` after the call, then the
    returned value is the real value held in this Data object. If it equals
    ``LIBMUSCLE_runtime_error``, the Data value did not contain a real
    value. If you passed an ``err_msg`` argument as well, then the passed
    variable will contain an appropriate error message in case of error, and
    needs to be deallocated (using ``deallocate()``) when you're done with it.

    Example:

    .. code-block:: fortran

        type(LIBMUSCLE_Data) :: mydata
        real(kind=LIBMUSCLE_real8) :: number

        ! Create a data object containing a real value
        mydata = LIBMUSCLE_Data_create(42.0d0)
        ! Retrieve the value
        number = LIBMUSCLE_Data_as_real8(mydata)
        ! Free the data object
        call LIBMUSCLE_Data_free(mydata)

    See :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    :p LIBMUSCLE_Data self: The Data object to get a double-precision real value
            out of.
    :p integer err_code: An error code output (optional)
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value.
    :rtype value: real(kind=LIBMUSCLE_real8)

.. f:function:: LIBMUSCLE_Data_as_settings(self, err_code, err_msg)

    Access a :f:type:`YMMSL_Settings` value.

    You can use :f:func:`LIBMUSCLE_Data_is_a_settings` to ascertain that the
    Data object contains a Settings-type value.

    If the Data object does not contain a :f:type:`YMMSL_Settings` value, then
    an error message will be printed and execution will be halted.

    Alternatively, you can pass an argument for ``err_code``, or for both
    ``err_code`` and ``err_msg``, to catch the error.

    If ``err_code`` equals ``LIBMUSCLE_success`` after the call, then the
    returned value is the value held in this Data object. If it equals
    ``LIBMUSCLE_runtime_error``, the Data value did not contain a Settings
    value. If you passed an ``err_msg`` argument as well, then the passed
    variable will contain an appropriate error message in case of error, and
    needs to be deallocated (using ``deallocate()``) when you're done with it.

    Example:

    .. code-block:: fortran

        type(LIBMUSCLE_Data) :: mydata
        type(YMMSL_Settings) :: settings1, settings2

        ! Create a Settings object
        settings1 = YMMSL_Settings_create()
        ! Create a data object containing the Settings value
        mydata = LIBMUSCLE_Data_create(settings1)
        ! Retrieve the value
        settings2 = LIBMUSCLE_Data_as_settings(mydata)
        ! Free the data object and the settings objects
        call YMMSL_Settings_free(settings1)
        call YMMSL_Settings_free(settings2)
        call LIBMUSCLE_Data_free(mydata)

    See :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    :p LIBMUSCLE_Data self: The Data object to get a Settings value out of.
    :p integer err_code: An error code output (optional)
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value.
    :rtype value: type(YMMSL_Settings)

.. f:subroutine:: LIBMUSCLE_Data_as_byte_array(self, buf, err_code, err_msg)

    Access a byte array value by copying it into ``buf``.

    You can use :f:func:`LIBMUSCLE_Data_is_a_byte_array` to ascertain that the
    Data object contains a byte array value. You can use
    :f:func:`LIBMUSCLE_Data_size` to get the number of bytes stored.

    If the Data object does not contain a byte array (character array) value,
    then an error message will be printed and execution will be halted.

    Alternatively, you can pass an argument for ``err_code``, or for both
    ``err_code`` and ``err_msg``, to catch the error.

    If ``err_code`` equals ``LIBMUSCLE_success`` after the call, then the
    contents of this variable will have been copied into ``buf``. If it equals
    ``LIBMUSCLE_runtime_error``, the Data value did not contain a byte array
    value. If you passed an ``err_msg`` argument as well, then the passed
    variable will contain an appropriate error message in case of error, and
    needs to be deallocated (using ``deallocate()``) when you're done with it.

    See :f:func:`LIBMUSCLE_Data_create_byte_array` for an example of
    creating and extracting byte array values. See
    :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    :p LIBMUSCLE_Data self: The Data object to get a byte array out of.
    :p character buf: A buffer large enough to hold the contents of the data
            object.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).

.. f:function:: LIBMUSCLE_Data_get_item(self, i, err_code, err_msg)

    Access an item in a list.

    This function is only valid for Data objects containing a list. You
    can use :f:func:`LIBMUSCLE_Data_is_a_list` to check whether that is the
    case.

    This returns a :f:type:`LIBMUSCLE_Data` object containing the value at the
    given index in the list object. If ``self`` does not contain a list, the
    result will be invalid, and ``err_code`` will be set to
    ``LIBMUSCLE_runtime_error``. If ``i`` is negative, zero, or larger than the
    number of items in the list (see :f:func:`LIBMUSCLE_Data_size`),
    ``err_code`` will be set to ``LIBMUSCLE_out_of_range``, and the result will
    be invalid.

    As with any returned :f:type:`LIBMUSCLE_Data` object, the result needs to be
    freed via :f:func:`LIBMUSCLE_Data_free` once you're done with it. Setting
    the value of the returned object will update the list, but it's easier and
    safer to use :f:func:`LIBMUSCLE_Data_set_item` instead.

    Example:

    .. code-block:: fortran

        type(LIBMUSCLE_Data) :: d1, d2

        d1 = LIBMUSCLE_Data_create_nils(10_LIBMUSCLE_size)
        d2 = LIBMUSCLE_Data_get_item(d1, 5_LIBMUSCLE_size)
        ! LIBMUSCLE_Data_is_nil(d2) returns .true. here
        call LIBMUSCLE_Data_free(d2)
        call LIBMUSCLE_Data_free(d1)

    See :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    :p LIBMUSCLE_Data self: The Data object to get an item out of.
    :p integer i: The index to get the value at, in range [1..size]
            (kind=LIBMUSCLE_size)
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value at the corresponding index.
    :rtype value: LIBMUSCLE_Data

.. f:function:: LIBMUSCLE_Data_get_item(self, key, err_code, err_msg)

    Access an item in a dictionary.

    This function is only valid for Data objects containing a dictionary. You
    can use :f:func:`LIBMUSCLE_Data_is_a_dict` to check whether that is the
    case.

    This returns a :f:type:`LIBMUSCLE_Data` object containing the value
    associated with the given key in the dictionary object. If ``self`` does not
    contain a dictionary, the result will be invalid, and ``err_code`` will be
    set to ``LIBMUSCLE_runtime_error``. If ``key`` does not exist in this
    dictionary, ``err_code`` will be set to ``LIBMUSCLE_out_of_range``, and the
    result will be invalid.

    As with any returned :f:type:`LIBMUSCLE_Data` object, the result needs to be
    freed via :f:func:`LIBMUSCLE_Data_free` once you're done with it. Note that
    the returned object will be invalidated if a new key is added to the
    dictionary.  Assigning to the returned object will update the dictionary,
    but it's easier and safer to use :f:func:`LIBMUSCLE_Data_set_item` instead.

    Example:

    .. code-block:: fortran

        type(LIBMUSCLE_Data) :: d1, d2, d3
        character(len=:), allocatable :: s1

        d1 = LIBMUSCLE_Data_create_dict()
        call LIBMUSCLE_Data_set_item(d1, 'key1', 'value1')
        d3 = LIBMUSCLE_Data_get_item(d1, 'key1')
        s1 = LIBMUSCLE_Data_as_character(d3)
        print *, s1     ! prints 'value1'
        deallocate(s1)
        call LIBMUSCLE_Data_free(d3)
        call LIBMUSCLE_Data_free(d2)
        call LIBMUSCLE_Data_free(d1)

    See :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    :p LIBMUSCLE_Data self: The Data object to get an item out of.
    :p character key: The key to get the value for.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value corresponding to the selected key.
    :rtype value: LIBMUSCLE_Data

.. f:subroutine:: LIBMUSCLE_Data_set_item(self, i, value, err_code, err_msg)

    Set an item in a list.

    This function is only valid for Data objects containing a list. You can
    use :f:func:`LIBMUSCLE_Data_is_a_list` to check whether that is the case.

    This subroutine sets the ``i``'th value in the list to ``value``. If a value
    is already stored at this position, then it will be replaced. If the Data
    object does not contain a list, ``err_code`` will be set to
    ``LIBMUSCLE_runtime_error``. If the position ``i`` is zero, negative, or
    larger than the size of the list, ``err_code`` will be set to
    ``LIBMUSCLE_out_of_range``.

    ``value`` may be of type logical, character, integer, real, or Data. See
    :f:func:`LIBMUSCLE_Data_get_item` for an example. See
    :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    :p LIBMUSCLE_Data self: The Data object to set an item value on.
    :p integer i: The position to set the value for, in range [1..size]
            (kind=LIBMUSCLE_size).
    :p see_above value: The value to set.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional)

.. f:subroutine:: LIBMUSCLE_Data_set_item(self, key, value, err_code, err_msg)

    Set an item in a dictionary.

    This function is only valid for Data objects containing a dictionary. You
    can use :f:func:`LIBMUSCLE_Data_is_a_dict` to check whether that is the
    case.

    This subroutine sets the value stored under ``key`` to ``value``. If a value
    is already stored under this key, then it will be replaced. If the Data
    object does not contain a dictionary, ``err_code`` will be set to
    ``LIBMUSCLE_runtime_error``.

    ``value`` may be of type logical, character, integer, real, or Data. See
    :f:func:`LIBMUSCLE_Data_get_item` for an example. See
    :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    :p LIBMUSCLE_Data self: The Data object to set an item value on.
    :p character key: The key to set the value for.
    :p see_above value: The value to set.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional)

.. f:function:: LIBMUSCLE_Data_key(self, i, err_code, err_msg)

    Get the i'th key in the dictionary.

    This function is only valid for Data objects containing a dictionary. You
    can use :f:func:`LIBMUSCLE_Data_is_a_dict` to check whether that is the
    case.

    The indices range from 1 to the number of items in the dictionary
    (inclusive), as usual in Fortran. Use :f:func:`LIBMUSCLE_Data_size` to get
    the number of items. Note that changes to the dictionary (e.g. inserting a
    new key) may change the order in which the key-value pairs are retrieved by
    this function. It's best to not change the dictionary while iterating
    through it.

    As always when a character value is returned by MUSCLE, the variable it ends
    up in must be allocatable, and must be deallocated after use.

    The corresponding value may be obtained via
    :f:func:`LIBMUSCLE_Data_value(i)`.

    .. code-block:: fortran

        type(LIBMUSCLE_Data) :: d1, val
        character(len=:), allocatable :: key, cval
        integer (kind=LIBMUSCLE_size) :: i
        integer intval

        d1 = LIBMUSCLE_Data_create_dict()
        call LIBMUSCLE_Data_set_item(d1, 'key1', 'value1')
        call LIBMUSCLE_Data_set_item(d1, 'key2', 'value2')

        do i = 1, LIBMUSCLE_Data_size(d1)
            key = LIBMUSCLE_Data_key(d1, i)
            val = LIBMUSCLE_Data_value(d1, i)
            cval = LIBMUSCLE_Data_as_character(val)
            print '(a8, a8)', key, cval
            deallocate(key)
            deallocate(cval)
            LIBMUSCLE_Data_free(val)
        end do

        call LIBMUSCLE_Data_free(d1)

    :p LIBMUSCLE_Data self: The Data object to get a key for.
    :p integer i: The index of the key to retrieve (LIBMUSCLE_size)
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r key: The key at the given index.
    :rtype key: character (allocatable)

.. f:function:: LIBMUSCLE_Data_value(self, i, err_code, err_msg)

    Get the i'th value in the dictionary.

    This function is only valid for Data objects containing a dictionary. You
    can use :f:func:`LIBMUSCLE_Data_is_a_dict` to check whether that is the
    case.

    The indices range from 1 to the number of items in the dictionary
    (inclusive), as usual in Fortran. Use :f:func:`LIBMUSCLE_Data_size` to get
    the number of items. Note that changes to the dictionary (e.g. inserting a
    new key) may change the order in which the key-value pairs are retrieved by
    this function. It's best to not change the dictionary while iterating
    through it.

    The corresponding key may be obtained via :f:func:`LIBMUSCLE_Data_key`. See
    there for an example as well.

    :p LIBMUSCLE_Data self: The Data object to get a value for.
    :p integer i: The index of the key to retrieve (LIBMUSCLE_size)
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value at the given index
    :rtype value: LIBMUSCLE_Data

.. f:function:: LIBMUSCLE_Data_num_dims(self, err_code, err_msg)

    Get the number of dimensions of a grid-valued Data object.

    This function is only valid for Data objects containing a grid. You
    can use :f:func:`LIBMUSCLE_Data_is_a_grid_of_logical` and similar
    functions to check that it is a grid. If the Data object does not
    contain a grid, ``err_code`` will be set to
    ``LIBMUSCLE_runtime_error``.

    :p LIBMUSCLE_Data self: The Data object to get the number of
            dimensions for.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r num_dims: The number of dimensions
    :rtype num_dims: integer (LIBMUSCLE_size)

.. f:subroutine:: LIBMUSCLE_Data_shape(self, shp, err_code, err_msg)

    Get the shape of the array of a grid-valued Data object.

    The array passed to receive the shape must be one-dimensional, and
    at least of length `n`, where `n` is the number of dimensions of
    the grid.

    This function is only valid for Data objects containing a grid. You
    can use :f:func:`LIBMUSCLE_Data_is_a_grid_of_logical` and similar
    functions to check that it is a grid. If the Data object does not
    contain a grid, ``err_code`` will be set to
    ``LIBMUSCLE_runtime_error``.

    :p LIBMUSCLE_Data self: The data object to get the shape of.
    :p integer shp: A 1D array of integer (LIBMUSCLE_size) to put the
            shape into (intent (out)).
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).

.. f:function:: LIBMUSCLE_Data_has_indexes(self, err_code, err_msg)

    Check whether a grid has index names.

    Returns ``.true.`` if the grid has index names, ``.false.`` otherwise.

    This function is only valid for Data objects containing a grid. You
    can use :f:func:`LIBMUSCLE_Data_is_a_grid_of_logical` and similar
    functions to check that it is a grid. If the Data object does not
    contain a grid, ``err_code`` will be set to
    ``LIBMUSCLE_runtime_error``.

    :p LIBMUSCLE_Data self: The data object to check for indexes.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r has_indexes: Whether there are indexes.
    :rtype has_indexes: logical

.. f:function:: LIBMUSCLE_Data_index(self, i, err_code, err_msg)

    Return the name of the i'th index.

    The value of ``i`` ranges from 1 to the number of dimensions.

    This function is only valid for Data objects containing a grid. You
    can use :f:func:`LIBMUSCLE_Data_is_a_grid_of_logical` and similar
    functions to check that it is a grid. If the Data object does not
    contain a grid, ``err_code`` will be set to
    ``LIBMUSCLE_runtime_error``. If the index is zero, negative, or
    larger than the number of dimensions, ``err_code`` will be set to
    ``LIBMUSCLE_out_of_range``.

    :p LIBMUSCLE_Data self: The data object to get the index of.
    :p integer i: The index of the index to get the name of (LIBMUSCLE_size).
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r index_name: The name of the index (allocatable)
    :rtype index_name: logical



LIBMUSCLE_Message
`````````````````
.. f:type:: LIBMUSCLE_Message

    Represents a message to be sent or that has been received.

    Messages have four attributes: A timestamp, an optional next timestamp,
    some data, and optional overlay settings.

    The timestamp holds the simulation time (not wallclock time) to which the
    data in this message corresponds. If a next timestamp is set, then that
    represents the simulation time of the next timestap of the model that
    generated the message. The data is the content of the message, and is
    model-specific. Overlay settings may be passed along with a message, and
    will be overlaid onto the receiving model's settings; this is normally only
    used by special simulation components.

.. f:function:: LIBMUSCLE_Message_create(timestamp, data)

    Create a new Message object.

    :p LIBMUSCLE_real8 timestamp: The simulated time to which the data in this
            message applies.
    :p LIBMUSCLE_Data data: An object to send or that was received

.. f:function:: LIBMUSCLE_Message_create(timestamp, next_timestamp, data)

    Create a new Message object.

    :p LIBMUSCLE_real8 timestamp: The simulated time to which the data in this
            message applies.
    :p LIBMUSCLE_real8 next_timestamp: Simulation time of the next message to
            be transmitted.
    :p LIBMUSCLE_Data data: An object to send or that was received

.. f:function:: LIBMUSCLE_Message_create(timestamp, data, settings)

    Create a new Message object.

    :p LIBMUSCLE_real8 timestamp: The simulated time to which the data in this
            message applies.
    :p LIBMUSCLE_Data data: An object to send or that was received
    :p YMMSL_Settings settings: Overlay settings to send or that were received.

.. f:function:: LIBMUSCLE_Message_create(timestamp, next_timestamp, data, settings)

    Create a new Message object.

    :p LIBMUSCLE_real8 timestamp: The simulated time to which the data in this
            message applies.
    :p LIBMUSCLE_real8 next_timestamp: Simulation time of the next message to
            be transmitted.
    :p LIBMUSCLE_Data data: An object to send or that was received
    :p YMMSL_Settings settings: Overlay settings to send or that were received.

.. f:function:: LIBMUSCLE_Message_timestamp(self)

    Returns the timestamp associated with the message.

    This is the simulated time to which the data in this message applies.

    :r timestamp: The timestamp.
    :rtype timestamp: LIBMUSCLE_real8

.. f:subroutine:: LIBMUSCLE_Message_set_timestamp(self, timestamp)

    Sets the timestamp associated with this message.

    This should be the simulated time to which the data in this message applies.

    :p LIBMUSCLE_Message self: The Message object to modify.
    :p timestamp: The timestamp to set.

.. f:function:: LIBMUSCLE_Message_has_next_timestamp(self)

    Returns whether the message has a next timestamp set.

    :p LIBMUSCLE_Message self: The Message object to inspect.
    :r has_next_timestamp: ``.true.`` if there's a next timestamp.
    :rtype has_next_timestamp: logical

.. f:function:: LIBMUSCLE_Message_next_timestamp(self)

    Returns the message's next timestamp.

    Only call if :f:func:`LIBMUSCLE_Message_has_next_timestamp` returns
    ``.true.``.

    :p LIBMUSCLE_Message self: The Message object to inspect.
    :r next_timestamp: The next timestamp for this message.
    :rtype next_timestamp: LIBMUSCLE_real8

.. f:subroutine:: LIBMUSCLE_Message_set_next_timestamp(self, timestamp)

    Sets the next timestamp associated with this message.

    This should be the simulated time of the next timestep of the model.

    :p LIBMUSCLE_Message self: The Message object to modify.
    :p timestamp: The timestamp to set.

.. f:subroutine:: LIBMUSCLE_Message_unset_next_timestamp(self)

    Unsets the next timestamp associated with this message.

    After calling this, :f:func:`LIBMUSCLE_Message_has_next_timestamp` will
    return ``.false.``.

    :p LIBMUSCLE_Message self: The Message object to modify.

.. f:function:: LIBMUSCLE_Message_get_data(self)

    Returns the data contained in the message.

    :p LIBMUSCLE_Message self: The Message object to inspect.
    :r data: The data contained in this message
    :rtype data: LIBMUSCLE_DataConstRef

.. f:subroutine:: LIBMUSCLE_Message_set_data(self, data)

    Sets the data contained by the message.

    Note that this will not transfer ownership of the data object, you still
    need to free it.

    :p LIBMUSCLE_Message self: The Message object to modify.
    :p LIBMUSCLE_Data data: The data object to take the value from.

.. f:subroutine:: LIBMUSCLE_Message_set_data(self, data)

    Sets the data contained by the message.

    Note that this will not transfer ownership of the data object, you still
    need to free it.

    :p LIBMUSCLE_Message self: The Message object to modify.
    :p LIBMUSCLE_DataConstRef data: The data object to take the value from.

.. f:function:: LIBMUSCLE_Message_has_settings(self)

    Returns whether the message has an associated Settings object.

    :p LIBMUSCLE_Message self: The Message object to inspect.
    :r has: ``.true.`` iff the message has settings.
    :rtype has: LIBMUSCLE_DataConstRef

.. f:function:: LIBMUSCLE_Message_get_settings(self)

    Returns the message's associated Settings object.

    Only call if :f:func:`LIBMUSCLE_Message_has_settings` returns ``.true.``.

    :p LIBMUSCLE_Message self: The Message object to inspect.
    :r settings: (A copy of) the associated settings object.
    :rtype settings: YMMSL_Settings

.. f:subroutine:: LIBMUSCLE_Message_set_settings(self, settings)

    Sets the message's associated Settings object.

    If the message has settings already, then they will be replaced by the new
    settings. After calling this, :f:func:`LIBMUSCLE_Message_has_settings` will
    return ``.true.``.

    :p LIBMUSCLE_Message self: The Message object to modify.
    :p YMMSL_Settings settings: The new settings.

.. f:subroutine:: LIBMUSCLE_Message_unset_settings(self)

    Removes any associated settings object from the message.

    This may be called whether the message currently has associated settings or
    not. After calling this function, :f:func:`LIBMUSCLE_Message_has_settings`
    will return ``.false.``.

    :p LIBMUSCLE_Message self: The Message object to modify.

LIBMUSCLE_PortsDescription
``````````````````````````
.. f:type:: LIBMUSCLE_PortsDescription

    Describes the ports of a compute element.

    This data structure is passed to libmuscle to describe how a
    compute element connects to the outside world, or it can be
    obtained from libmuscle in order to find out how a compute
    element with flexible ports was used.

    A PortsDescription contains a list of port names for each
    :f:type:`YMMSL_Operator`.

.. f:function:: LIBMUSCLE_PortsDescription_create()

    Create a PortsDescription containing no port names.

    :r ports_description: A new PortsDescription object.
    :rtype ports_description: LIBMUSCLE_PortsDescription

.. f:subroutine:: LIBMUSCLE_PortsDescription_free(self)

    Frees a PortsDescription object.

    This deallocates all resources associated with the object, and
    should be called for every PortsDescription object when you're
    done using it.

    :p LIBMUSCLE_PortsDescription self: The object to free.

.. f:subroutine:: LIBMUSCLE_PortsDescription_add(self, operator, port)

    Add a port name to a PortsDescription object.

    :p LIBMUSCLE_PortsDescription self: The object to modify.
    :p YMMSL_Operator operator: The operator under which to put the port.
    :p character port: The name of the port to add.

.. f:function:: LIBMUSCLE_PortsDescription_num_ports(self, operator)

    Get the number of ports in this object for the given operator.

    :p LIBMUSCLE_PortsDescription self: The object to inspect.
    :p YMMSL_Operator operator: A chosen operator.
    :r num_ports: The number of ports for this operator (of kind
            LIBMUSCLE_size).
    :rtype num_ports: integer

.. f:function:: LIBMUSCLE_PortsDescription_get(self, operator, i)

    Get the i'th port name for the given operator.

    Parameter ``i`` must be in the range 1..num_ports inclusive, where
    num_ports is the result of calling
    :f:func:`LIBMUSCLE_PortsDescription_num_ports` with the same
    object and operator.

    :p LIBMUSCLE_PortsDescription self: The object to inspect.
    :p YMMSL_Operator operator: A chosen operator.
    :r port_name: The name of the given port.
    :rtype port_name: character

LIBMUSCLE_Instance
``````````````````````````
.. f:type:: LIBMUSCLE_Instance

    The Instance class represents a compute element instance in a
    MUSCLE3 simulation. This class provides a low-level send/receive API for the
    instance to use.

.. f:function:: LIBMUSCLE_Instance_create()

    Create a new Instance object with ports from the configuration.

    For MPI-based compute elements, this will have libmuscle_mpi use a duplicate
    of ``MPI_COMM_WORLD`` to communicate, and the designated root process will
    be that with rank 0.

    This object must be freed when you're done with it using
    :f:func:`LIBMUSCLE_Instance_free`.

    :r instance: The newly created instance object.
    :rtype instance: LIBMUSCLE_Instance

.. f:function:: LIBMUSCLE_Instance_create(ports)

    Create a new Instance object with the given ports.

    For MPI-based compute elements, this will have libmuscle_mpi use a duplicate
    of ``MPI_COMM_WORLD`` to communicate, and the designated root process will
    be that with rank 0.

    This object must be freed when you're done with it using
    :f:func:`LIBMUSCLE_Instance_free`.

    :p LIBMUSCLE_PortsDescription ports: The ports of the new instance.
    :r instance: The newly created instance object.
    :rtype instance: LIBMUSCLE_Instance

.. f:function:: LIBMUSCLE_Instance_create(communicator, root)

    Create a new Instance object for MPI with ports from the configuration.

    For MPI-based compute elements, an MPI communicator and a root rank may be
    passed. The communicator must contain all processes in this instance, and
    ``root`` must be the rank of one of them. MUSCLE will create a duplicate of
    this communicator for its own use. Creating a :f:type:`LIBMUSCLE_Instance`
    for an MPI compute element is a collective operation, so it must be done in
    all processes simultaneously, with the same communicator and the same root.

    This object must be freed when you're done with it using
    :f:func:`LIBMUSCLE_Instance_free`.

    :p integer communicator: MPI communicator to use (optional, default
            MPI_COMM_WORLD).
    :p integer root: Rank of the root process (optional, default 0).
    :r instance: The newly created instance object.
    :rtype instance: LIBMUSCLE_Instance

.. f:function:: LIBMUSCLE_Instance_create(ports, communicator, root)

    Create a new Instance object for MPI with the given ports.

    For MPI-based compute elements, an MPI communicator and a root rank may be
    passed. The communicator must contain all processes in this instance, and
    ``root`` must be the rank of one of them. MUSCLE will create a duplicate of
    this communicator for its own use. Creating a :f:type:`LIBMUSCLE_Instance`
    for an MPI compute element is a collective operation, so it must be done in
    all processes simultaneously, with the same communicator and the same root.

    This object must be freed when you're done with it using
    :f:func:`LIBMUSCLE_Instance_free`.

    :p LIBMUSCLE_PortsDescription ports: The ports of the new instance.
    :p integer communicator: MPI communicator to use (optional, default
            MPI_COMM_WORLD).
    :p integer root: Rank of the root process (optional, default 0).
    :r instance: The newly created instance object.
    :rtype instance: LIBMUSCLE_Instance

.. f:subroutine:: LIBMUSCLE_Instance_free(self)

    Free resources associated with the given Instance object.

    :p LIBMUSCLE_Instance self: The object to free.

.. f:function:: LIBMUSCLE_Instance_reuse_instance(self)

    Checks whether to reuse this instance.

    This method must be called at the beginning of the reuse loop, i.e. before
    the F_INIT operator, and its return value should decide whether to enter
    that loop again.

    MPI-based compute elements must execute the reuse loop in each process in
    parallel, and call this function at the top of the reuse loop in each
    process.

    :p LIBMUSCLE_Instance self: The object to check for reuse.
    :r reuse: Whether to enter the reuse loop another time.
    :rtype reuse: logical

.. f:function:: LIBMUSCLE_Instance_reuse_instance(self, apply_overlay)

    Checks whether to reuse this instance.

    This method must be called at the beginning of the reuse loop, i.e. before
    the F_INIT operator, and its return value should decide whether to enter
    that loop again.

    This version of this function lets you choose whether to apply the received
    settings overlay or to return it with the message. If you're going to use
    :f:func:`LIBMUSCLE_Instance_receive_with_settings` on your F_INIT ports, set
    this to ``.false.``. If you don't know what that means, just call
    ``LIBMUSCLE_Instance_reuse_instance()`` with no arguments and all will be
    fine. If it turns out that you did need to specify ``.false.`` here, MUSCLE
    3 will tell you in an error message, and you can add it.

    :p LIBMUSCLE_Instance self: The object to check for reuse.
    :p logical apply_overlay: Whether to apply the received settings overlay.
    :r reuse: Whether to enter the reuse loop another time.
    :rtype reuse: logical

.. f:subroutine:: LIBMUSCLE_Instance_error_shutdown(self, message)

    Logs an error and shuts down the Instance.

    If you detect that something is wrong and want to stop this instance, then
    you should call this function to shut down the instance before stopping the
    program. This makes debugging easier.

    MPI-based compute elements may either call this function in all processes,
    or only in the root process (as passed to the constructor).

    :p LIBMUSCLE_Instance self: The instance to shut down.
    :p character message: An error message describing the problem encountered.

.. f:function:: LIBMUSCLE_Instance_is_setting_a_character(self, name, err_code, err_msg)

    Returns whether the setting is of type character.

    If no setting with the given name exists, ``err_code`` will be set to
    ``LIBMUSCLE_out_of_range``. See :f:func:`LIBMUSCLE_Data_as_logical` for an
    example of error handling.

    MPI-based compute elements may call this function at any time within the
    reuse loop, in any or all processes, simultaneously or not.

    :p LIBMUSCLE_Instance self: The instance to get the setting from.
    :p character name: The name of the setting to inspect.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r is_character: .true. if the setting is of character type.
    :rtype is_character: logical

.. f:function:: LIBMUSCLE_Instance_is_setting_a_int8(self, name, err_code, err_msg)

    Returns whether the setting is of type integer.

    If no setting with the given name exists, ``err_code`` will be set to
    ``LIBMUSCLE_out_of_range``. See :f:func:`LIBMUSCLE_Data_as_logical` for an
    example of error handling.

    MPI-based compute elements may call this function at any time within the
    reuse loop, in any or all processes, simultaneously or not.

    :p LIBMUSCLE_Instance self: The instance to get the setting from.
    :p character name: The name of the setting to inspect.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r is_int8: .true. if the setting is of integer type.
    :rtype is_int8: logical

.. f:function:: LIBMUSCLE_Instance_is_setting_a_real8(self, name, err_code, err_msg)

    Returns whether the setting is of type real.

    If no setting with the given name exists, ``err_code`` will be set to
    ``LIBMUSCLE_out_of_range``. See :f:func:`LIBMUSCLE_Data_as_logical` for an
    example of error handling.

    MPI-based compute elements may call this function at any time within the
    reuse loop, in any or all processes, simultaneously or not.

    :p LIBMUSCLE_Instance self: The instance to get the setting from.
    :p character name: The name of the setting to inspect.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r is_real8: .true. if the setting is of real type.
    :rtype is_real8: logical

.. f:function:: LIBMUSCLE_Instance_is_setting_a_logical(self, name, err_code, err_msg)

    Returns whether the setting is of type logical.

    If no setting with the given name exists, ``err_code`` will be set to
    ``LIBMUSCLE_out_of_range``. See :f:func:`LIBMUSCLE_Data_as_logical` for an
    example of error handling.

    MPI-based compute elements may call this function at any time within the
    reuse loop, in any or all processes, simultaneously or not.

    :p LIBMUSCLE_Instance self: The instance to get the setting from.
    :p character name: The name of the setting to inspect.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r is_logical: .true. if the setting is of logical type.
    :rtype is_logical: logical

.. f:function:: LIBMUSCLE_Instance_is_setting_a_real8array(self, name, err_code, err_msg)

    Returns whether the setting is a 1D array of real.

    If no setting with the given name exists, ``err_code`` will be set to
    ``LIBMUSCLE_out_of_range``. See :f:func:`LIBMUSCLE_Data_as_logical` for an
    example of error handling.

    MPI-based compute elements may call this function at any time within the
    reuse loop, in any or all processes, simultaneously or not.

    :p LIBMUSCLE_Instance self: The instance to get the setting from.
    :p character name: The name of the setting to inspect.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r is_real8array: .true. if the setting is of real8array type.
    :rtype is_real8array: logical

.. f:function:: LIBMUSCLE_Instance_is_setting_a_real8array2(self, name, err_code, err_msg)

    Returns whether the setting is a 2D array of real.

    If no setting with the given name exists, ``err_code`` will be set to
    ``LIBMUSCLE_out_of_range``. See :f:func:`LIBMUSCLE_Data_as_logical` for an
    example of error handling.

    MPI-based compute elements may call this function at any time within the
    reuse loop, in any or all processes, simultaneously or not.

    :p LIBMUSCLE_Instance self: The instance to get the setting from.
    :p character name: The name of the setting to inspect.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r is_real8array2: .true. if the setting is of real8array2 type.
    :rtype is_real8array2: logical

.. f:function:: LIBMUSCLE_Instance_get_setting_as_character(self, name, err_code, err_msg)

    Returns the value of a character-valued model setting.

    If no setting with the given name exists, ``err_code`` will be set to
    ``LIBMUSCLE_out_of_range``. If the value is not of type character, then
    ``err_code`` will be set to ``LIBMUSCLE_bad_cast``. See
    :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    MPI-based compute elements may call this function at any time within the
    reuse loop, in any or all processes, simultaneously or not.

    :p LIBMUSCLE_Instance self: The instance to get the setting from.
    :p character name: The name of the setting to retrieve.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The setting's value.
    :rtype value: character

.. f:function:: LIBMUSCLE_Instance_get_setting_as_int8(self, name, err_code, err_msg)

    Returns the value of an integer-valued model setting.

    If no setting with the given name exists, ``err_code`` will be set to
    ``LIBMUSCLE_out_of_range``. If the value is not of type integer, then
    ``err_code`` will be set to ``LIBMUSCLE_bad_cast``. See
    :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    MPI-based compute elements may call this function at any time within the
    reuse loop, in any or all processes, simultaneously or not.

    :p LIBMUSCLE_Instance self: The instance to get the setting from.
    :p character name: The name of the setting to retrieve.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The setting's value (kind=LIBMUSCLE_int8).
    :rtype value: integer

.. f:function:: LIBMUSCLE_Instance_get_setting_as_real8(self, name, err_code, err_msg)

    Returns the value of a real-valued model setting.

    If no setting with the given name exists, ``err_code`` will be set to
    ``LIBMUSCLE_out_of_range``. If the value is not of type integer, then
    ``err_code`` will be set to ``LIBMUSCLE_bad_cast``. See
    :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    MPI-based compute elements may call this function at any time within the
    reuse loop, in any or all processes, simultaneously or not.

    :p LIBMUSCLE_Instance self: The instance to get the setting from.
    :p character name: The name of the setting to retrieve.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The setting's value (kind=LIBMUSCLE_real8).
    :rtype value: real

.. f:function:: LIBMUSCLE_Instance_get_setting_as_logical(self, name, err_code, err_msg)

    Returns the value of a logical-valued model setting.

    If no setting with the given name exists, ``err_code`` will be set to
    ``LIBMUSCLE_out_of_range``. If the value is not of type integer, then
    ``err_code`` will be set to ``LIBMUSCLE_bad_cast``. See
    :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    MPI-based compute elements may call this function at any time within the
    reuse loop, in any or all processes, simultaneously or not.

    :p LIBMUSCLE_Instance self: The instance to get the setting from.
    :p character name: The name of the setting to retrieve.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The setting's value.
    :rtype value: logical

.. f:subroutine:: LIBMUSCLE_Instance_get_setting_as_real8array(self, name, value, err_code, err_msg)

    Returns the value of an array-of-real8-valued model setting.

    Note that there is currently no way to get the size of the array in advance.
    This feature is intended to be used for small fixed arrays, in which case
    the size will be known in advance to the programmer.

    If no setting with the given name exists, ``err_code`` will be set to
    ``LIBMUSCLE_out_of_range``. If the value is not a 1D array of real8, then
    ``err_code`` will be set to ``LIBMUSCLE_bad_cast``. See
    :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    MPI-based compute elements may call this function at any time within the
    reuse loop, in any or all processes, simultaneously or not.

    :p LIBMUSCLE_Instance self: The instance to get the setting from.
    :p character name: The name of the setting to retrieve.
    :p LIBMUSCLE_real8: The returned value (out, dimension(:))
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).

.. f:subroutine:: LIBMUSCLE_Instance_get_setting_as_real8array2(self, name, value, err_code, err_msg)

    Returns the value of a 2D-array-of-real8-valued model setting.

    Note that there is currently no way to get the size of the array in advance.
    This feature is intended to be used for small fixed arrays, in which case
    the size will be known in advance to the programmer.

    If no setting with the given name exists, ``err_code`` will be set to
    ``LIBMUSCLE_out_of_range``. If the value is not a 2D array of real8, then
    ``err_code`` will be set to ``LIBMUSCLE_bad_cast``. See
    :f:func:`LIBMUSCLE_Data_as_logical` for an example of error handling.

    MPI-based compute elements may call this function at any time within the
    reuse loop, in any or all processes, simultaneously or not.

    :p LIBMUSCLE_Instance self: The instance to get the setting from.
    :p character name: The name of the setting to retrieve.
    :p LIBMUSCLE_real8: The returned value (out, dimension(:,:))
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).

.. f:function:: LIBMUSCLE_Instance_list_ports(self)

    Returns a description of the ports of this instance.

    MPI-based compute elements may call this function only in the root process.

    :p LIBMUSCLE_Instance self: The instance whose ports to describe.
    :r ports: A description of the ports, organised by operator.
    :rtype ports: LIBMUSCLE_PortsDescription

.. f:function:: LIBMUSCLE_Instance_is_connected(self, port)

    Returns whether the given port is connected.

    MPI-based compute elements may call this function only in the root process.

    :p LIBMUSCLE_Instance self: The instance to inspect.
    :p character port: The name of the port to inspect
    :r connected: ``.true.`` if the port is connected.
    :rtype connected: logical

.. f:function:: LIBMUSCLE_Instance_is_vector_port(self, port)

    Returns whether the given port is a vector port.

    MPI-based compute elements may call this function only in the root process.

    :p LIBMUSCLE_Instance self: The instance to inspect.
    :p character port: The name of the port to inspect
    :r connected: ``.true.`` if the port is a vector port.
    :rtype connected: logical

.. f:function:: LIBMUSCLE_Instance_is_resizable(self, port, err_code, err_msg)

    Returns whether the port is resizable.

    This function must only be called on vector ports. If the port is a scalar
    port, ``err_code`` will be set to ``LIBMUSCLE_runtime_error`` and the return
    value will be invalid.

    MPI-based compute elements may call this function only in the root process.

    :p LIBMUSCLE_Instance self: The instance to inspect.
    :p character port: The name of the port to inspect
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r connected: ``.true.`` if the port is a vector port.
    :rtype connected: logical

.. f:function:: LIBMUSCLE_Instance_get_port_length(self, port, err_code, err_msg)

    Returns the current length of a vector port.

    This function must only be called on vector ports. If the port is a scalar
    port, ``err_code`` will be set to ``LIBMUSCLE_runtime_error`` and the return
    value will be invalid.

    MPI-based compute elements may call this function only in the root process.

    :p LIBMUSCLE_Instance self: The instance to inspect.
    :p character port: The name of the port to inspect
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r connected: ``.true.`` if the port is a vector port.
    :rtype connected: logical

.. f:subroutine:: LIBMUSCLE_Instance_set_port_length(self, port, length)

    Sets the current length of a vector port.

    This function must only be called on resizable vector ports. If the port is
    a scalar port or a non-resizable vector port, ``err_code`` will be set to
    ``LIBMUSCLE_runtime_error`` and the return value will be invalid.

    MPI-based compute elements may call this function only in the root process.

    :p LIBMUSCLE_Instance self: The instance to change a port on.
    :p character port: The name of the port to modify.
    :p integer length: The new length of the port.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).

.. f:subroutine:: LIBMUSCLE_Instance_send(self, port_name, message)

    Send a message to the outside world via a port.

    Sending is non-blocking, a copy of the message will be made and stored until
    the receiver is ready to receive it.

    MPI-based compute elements may call this function either in all processes,
    or only in the root process. In both cases, the message given by the root
    process will be sent, the others ignored. You may want to do a gather
    operation first to collect all the information that is to be sent in the
    root process.

    :p LIBMUSCLE_Instance self: The instance to send a message from.
    :p character port_name: The name of the port to send on.
    :p LIBMUSCLE_Message message: The message to send.

.. f:subroutine:: LIBMUSCLE_Instance_send(self, port_name, message, slot)

    Send a message to the outside world via a slot on a port.

    Sending is non-blocking, a copy of the message will be made and stored until
    the receiver is ready to receive it.

    MPI-based compute elements may call this function either in all processes,
    or only in the root process. In both cases, the message given by the root
    process will be sent, the others ignored. You may want to do a gather
    operation first to collect all the information that is to be sent in the
    root process.

    :p LIBMUSCLE_Instance self: The instance to send a message from.
    :p character port_name: The name of the port to send on.
    :p int slot: The slot to send on.
    :p LIBMUSCLE_Message message: The message to send.

.. f:function:: LIBMUSCLE_Instance_receive(self, port_name, err_code, err_msg)

    Receive a message from the outside world via a port.

    Receiving is a blocking operation. This function will contact the sender,
    wait for a message to be available, and receive and return it. Note that you
    must free the returned :f:type:`LIBMUSCLE_Message` object when you're done
    with it.

    If the port you are receiving on is not connected, then err_code will be set
    to ``LIBMUSCLE_runtime_error``.

    MPI-based compute elements must call this function in all processes
    simultaneously. The received message will be returned in the root process,
    all other processes will receive a dummy message. It is therefore up to the
    model code to scatter or broadcast the received message to the non-root
    processes, if necessary.

    :p LIBMUSCLE_Instance self: The instance to receive a message for.
    :p character port_name: The name of the port to receive on.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r message: The received message.
    :rtype message: LIBMUSCLE_Message

.. f:function:: LIBMUSCLE_Instance_receive(self, port_name, default_msg, err_code, err_msg)

    Receive a message from the outside world via a port.

    Receiving is a blocking operation. This function will contact the sender,
    wait for a message to be available, and receive and return it. Note that you
    must free the returned :f:type:`LIBMUSCLE_Message` object when you're done
    with it.

    If the port you are receiving on is not connected, then a copy of
    ``default_msg`` will be returned.

    MPI-based compute elements must call this function in all processes
    simultaneously. The received message will be returned in the root process,
    all other processes will receive a dummy message. It is therefore up to the
    model code to scatter or broadcast the received message to the non-root
    processes, if necessary.

    :p LIBMUSCLE_Instance self: The instance to receive a message for.
    :p character port_name: The name of the port to receive on.
    :p LIBMUSCLE_Message default_msg: A default message in case the port is not
            connected.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r message: The received message.
    :rtype message: LIBMUSCLE_Message

.. f:function:: LIBMUSCLE_Instance_receive_on_slot(self, port_name, slot, err_code, err_msg)

    Receive a message from the outside world via a slot on a vector port.

    Receiving is a blocking operation. This function will contact the sender,
    wait for a message to be available, and receive and return it. Note that you
    must free the returned :f:type:`LIBMUSCLE_Message` object when you're done
    with it.

    If the port you are receiving on is not connected, then err_code will be set
    to ``LIBMUSCLE_runtime_error``.

    MPI-based compute elements must call this function in all processes
    simultaneously. The received message will be returned in the root process,
    all other processes will receive a dummy message. It is therefore up to the
    model code to scatter or broadcast the received message to the non-root
    processes, if necessary.

    MPI-based compute elements must call this function in all processes
    simultaneously. The received message will be returned in the root process,
    all other processes will receive a dummy message. It is therefore up to the
    model code to scatter or broadcast the received message to the non-root
    processes, if necessary.

    :p LIBMUSCLE_Instance self: The instance to receive a message for.
    :p character port_name: The name of the vector port to receive on.
    :p integer slot: The slot to receive on.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r message: The received message.
    :rtype message: LIBMUSCLE_Message

.. f:function:: LIBMUSCLE_Instance_receive_on_slot(self, port_name, slot, default_msg, err_code, err_msg)

    Receive a message from the outside world via a slot on a vector port.

    Receiving is a blocking operation. This function will contact the sender,
    wait for a message to be available, and receive and return it. Note that you
    must free the returned :f:type:`LIBMUSCLE_Message` object when you're done
    with it.

    If the port you are receiving on is not connected, then a copy of
    ``default_msg`` will be returned.

    MPI-based compute elements must call this function in all processes
    simultaneously. The received message will be returned in the root process,
    all other processes will receive a dummy message. It is therefore up to the
    model code to scatter or broadcast the received message to the non-root
    processes, if necessary.

    :p LIBMUSCLE_Instance self: The instance to receive a message for.
    :p character port_name: The name of the port to receive on.
    :p integer slot: The slot to receive on.
    :p LIBMUSCLE_Message default_msg: A default message in case the port is not
            connected.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r message: The received message.
    :rtype message: LIBMUSCLE_Message

.. f:function:: LIBMUSCLE_Instance_receive_with_settings(self, port_name, err_code, err_msg)

    Receive a message with attached settings overlay.

    Receiving is a blocking operation. This function will contact the sender,
    wait for a message to be available, and receive and return it. Note that you
    must free the returned :f:type:`LIBMUSCLE_Message` object when you're done
    with it.

    If the port you are receiving on is not connected, then err_code will be set
    to ``LIBMUSCLE_runtime_error``.

    MPI-based compute elements must call this function in all processes
    simultaneously. The received message will be returned in the root process,
    all other processes will receive a dummy message. It is therefore up to the
    model code to scatter or broadcast the received message to the non-root
    processes, if necessary.

    :p LIBMUSCLE_Instance self: The instance to receive a message for.
    :p character port_name: The name of the port to receive on.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r message: The received message.
    :rtype message: LIBMUSCLE_Message

.. f:function:: LIBMUSCLE_Instance_receive_with_settings(self, port_name, default_msg, err_code, err_msg)

    Receive a message with attached settings overlay.

    Receiving is a blocking operation. This function will contact the sender,
    wait for a message to be available, and receive and return it. Note that you
    must free the returned :f:type:`LIBMUSCLE_Message` object when you're done
    with it.

    If the port you are receiving on is not connected, then a copy of
    ``default_msg`` will be returned.

    MPI-based compute elements must call this function in all processes
    simultaneously. The received message will be returned in the root process,
    all other processes will receive a dummy message. It is therefore up to the
    model code to scatter or broadcast the received message to the non-root
    processes, if necessary.

    :p LIBMUSCLE_Instance self: The instance to receive a message for.
    :p character port_name: The name of the port to receive on.
    :p LIBMUSCLE_Message default_msg: A default message in case the port is not
            connected.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r message: The received message.
    :rtype message: LIBMUSCLE_Message

.. f:function:: LIBMUSCLE_Instance_receive_with_settings_on_slot(self, port_name, slot, err_code, err_msg)

    Receive a message with attached settings overlay.

    Receiving is a blocking operation. This function will contact the sender,
    wait for a message to be available, and receive and return it. Note that you
    must free the returned :f:type:`LIBMUSCLE_Message` object when you're done
    with it.

    If the port you are receiving on is not connected, then err_code will be set
    to ``LIBMUSCLE_runtime_error``.

    MPI-based compute elements must call this function in all processes
    simultaneously. The received message will be returned in the root process,
    all other processes will receive a dummy message. It is therefore up to the
    model code to scatter or broadcast the received message to the non-root
    processes, if necessary.

    :p LIBMUSCLE_Instance self: The instance to receive a message for.
    :p character port_name: The name of the vector port to receive on.
    :p integer slot: The slot to receive on.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r message: The received message.
    :rtype message: LIBMUSCLE_Message

.. f:function:: LIBMUSCLE_Instance_receive_with_settings_on_slot(self, port_name, slot, default_msg, err_code, err_msg)

    Receive a message with attached settings overlay.

    Receiving is a blocking operation. This function will contact the sender,
    wait for a message to be available, and receive and return it. Note that you
    must free the returned :f:type:`LIBMUSCLE_Message` object when you're done
    with it.

    If the port you are receiving on is not connected, then a copy of
    ``default_msg`` will be returned.

    MPI-based compute elements must call this function in all processes
    simultaneously. The received message will be returned in the root process,
    all other processes will receive a dummy message. It is therefore up to the
    model code to scatter or broadcast the received message to the non-root
    processes, if necessary.

    :p LIBMUSCLE_Instance self: The instance to receive a message for.
    :p character port_name: The name of the port to receive on.
    :p integer slot: The slot to receive on.
    :p LIBMUSCLE_Message default_msg: A default message in case the port is not
            connected.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r message: The received message.
    :rtype message: LIBMUSCLE_Message

Namespace YMMSL
---------------

.. f:module:: ymmsl
    :synopsis: Fortran module for ymmsl

.. f:currentmodule:: ymmsl

YMMSL_Operator
``````````````

YMMSL operators are represented in Python by integer constants of kind
``YMMSL_Operator``. The following values are available:

+----------+-----------------------+
| Operator | Constant              |
+==========+=======================+
|   None   | YMMSL_Operator_NONE   |
+----------+-----------------------+
|  f_init  | YMMSL_Operator_F_INIT |
+----------+-----------------------+
|   O_i    | YMMSL_Operator_O_I    |
+----------+-----------------------+
|    S     | YMMSL_Operator_S      |
+----------+-----------------------+
|    B     | YMMSL_Operator_B      |
+----------+-----------------------+
|   O_f    | YMMSL_Operator_O_F    |
+----------+-----------------------+

YMMSL_Settings
``````````````
.. f:type:: YMMSL_Settings

    Represents a libmuscle Settings object. These are used to send and receive
    Settings objects to other compute elements. This is an opaque object that
    may be returned from and passed to libmuscle functions, but does not contain
    any directly accessible members.

    A Settings object is a dictionary-like object which is indexed by a string,
    and whose values can be strings, logicals, 8-byte integers, 8-byte real
    numbers, and one- and two-dimensional arrays of 8-byte real numbers.

.. f:function:: YMMSL_Settings_create()

    Creates an empty Settings object.

    :r obj: The new Settings object
    :rtype obj: YMMSL_Settings

.. f:subroutine:: YMMSL_Settings_free(self)

    Frees a Settings object.

    This frees the resources associated with the given Settings object. Do not
    use the object for anything after calling this, because it will be invalid.

    :p YMMSL_Settings self: The Settings object to free.

.. f:function:: YMMSL_Settings_equals(self, other)

    Compares two Settings objects for equality.

    This returns ``.true.`` if and only if the two :f:type:`YMMSL_Settings`
    objects contain the same keys and values.

    :p YMMSL_Settings self: The object to compare.
    :p YMMSL_Settings other: The object to compare to.
    :r equal: ``.true.`` if the objects are equal.
    :rtype equal: logical

.. f:function:: YMMSL_Settings_size(self)

    Returns the number of settings in this object.

    :p YMMSL_Settings self: The object to inspect.
    :r count: The number of key-value pairs in this Settings object
            (kind=YMMSL_size).
    :rtype count: integer

.. f:function:: YMMSL_Settings_empty(self)

    Returns ``.true.`` if and only if the Settings object has no items.

    :p YMMSL_Settings self: The object to inspect.
    :r empty: Whether the object is empty.
    :rtype empty: logical

.. f:function:: YMMSL_Settings_contains(self, key)

    Returns ``.true.`` if the Settings object contains the given key.

    :p YMMSL_Settings self: The object to inspect.
    :r contains: Whether the given key exists in this Settings object.
    :rtype contains: logical

.. f:function:: YMMSL_Settings_is_a_character(self, key, err_code, err_msg)

    Return whether a value is of type character.

    If the given key does not exist, then ``err_code`` will be set to
    ``YMMSL_out_of_bounds`` and the result will be invalid.

    :p YMMSL_Settings self: The Settings object to inspect.
    :p character key: The name of the setting to check.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r is: ``.true.`` if the value is of type character.
    :rtype is: logical

.. f:function:: YMMSL_Settings_is_a_logical(self, key, err_code, err_msg)

    Return whether a value is of type logical.

    If the given key does not exist, then ``err_code`` will be set to
    ``YMMSL_out_of_bounds`` and the result will be invalid.

    :p YMMSL_Settings self: The Settings object to inspect.
    :p character key: The name of the setting to check.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r is: ``.true.`` if the value is of type logical.
    :rtype is: logical

.. f:function:: YMMSL_Settings_is_a_int4(self, key, err_code, err_msg)

    Return whether a value is of type ``YMMSL_int4``.

    This returns ``.true.`` if the value is an integer and fits in an int4.

    If the given key does not exist, then ``err_code`` will be set to
    ``YMMSL_out_of_bounds`` and the result will be invalid.

    :p YMMSL_Settings self: The Settings object to inspect.
    :p character key: The name of the setting to check.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r is: ``.true.`` if the value is of type ``YMMSL_int4``.
    :rtype is: logical

.. f:function:: YMMSL_Settings_is_a_int8(self, key, err_code, err_msg)

    Return whether a value is of type ``YMMSL_int8``.

    If the given key does not exist, then ``err_code`` will be set to
    ``YMMSL_out_of_bounds`` and the result will be invalid.

    :p YMMSL_Settings self: The Settings object to inspect.
    :p character key: The name of the setting to check.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r is: ``.true.`` if the value is of type ``YMMSL_int8``.
    :rtype is: logical

.. f:function:: YMMSL_Settings_is_a_real8(self, key, err_code, err_msg)

    Return whether a value is of type ``YMMSL_real8``.

    This will also return ``.true.`` if the value is an integer, even if
    converting it would lose precision.

    If the given key does not exist, then ``err_code`` will be set to
    ``YMMSL_out_of_bounds`` and the result will be invalid.

    :p YMMSL_Settings self: The Settings object to inspect.
    :p character key: The name of the setting to check.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r is: ``.true.`` if the value is of type ``YMMSL_real8``.
    :rtype is: logical

.. f:function:: YMMSL_Settings_is_a_real8array(self, key, err_code, err_msg)

    Return whether a value is a 1D array of type ``YMMSL_real8``.

    If the given key does not exist, then ``err_code`` will be set to
    ``YMMSL_out_of_bounds`` and the result will be invalid.

    :p YMMSL_Settings self: The Settings object to inspect.
    :p character key: The name of the setting to check.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r is: ``.true.`` if the value is a 1D array of type ``YMMSL_real8``.
    :rtype is: logical

.. f:function:: YMMSL_Settings_is_a_real8array2(self, key, err_code, err_msg)

    Return whether a value is a 2D array of type ``YMMSL_real8``.

    If the given key does not exist, then ``err_code`` will be set to
    ``YMMSL_out_of_bounds`` and the result will be invalid.

    :p YMMSL_Settings self: The Settings object to inspect.
    :p character key: The name of the setting to check.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r is: ``.true.`` if the value is a 2D array of type ``YMMSL_real8``.
    :rtype is: logical

.. f:subroutine:: YMMSL_Settings_set(self, key, value)

    Sets a setting to the given value.

    If no setting with the given key exists, one is added, if one does,
    it is overwritten.

    ``value`` may be a character (string), a logical, a 4-byte integer (e.g.
    ``YMMSL_int4``), an 8-byte integer (e.g.  ``YMMSL_int8``), an 8-byte real
    number (``YMMSL_real8``), or a one- or two-dimensional arrays of 8-byte real
    numbers.

    :p YMMSL_Settings self: The Settings object to modify.
    :p character key: The name of the setting.
    :p see_above value: The value to set the setting to.

.. f:function:: YMMSL_Settings_get_as_character(self, key, err_code, err_msg)

    Return the value of a character-typed setting.

    If this setting is not currently set to a character-typed value,
    then ``err_code`` will be set to ``YMMSL_bad_cast`` and the result
    will be invalid.

    :p character key: The name of the setting to get.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value at the given index
    :rtype value: character

.. f:function:: YMMSL_Settings_get_as_logical(self, key, err_code, err_msg)

    Return the value of a logical-typed setting.

    If this setting is not currently set to a logical-typed value,
    then ``err_code`` will be set to ``YMMSL_bad_cast`` and the result
    will be invalid.

    :p character key: The name of the setting to get.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value at the given index
    :rtype value: logical

.. f:function:: YMMSL_Settings_get_as_int4(self, key, err_code, err_msg)

    Return the value of an integer-typed setting.

    If this setting is not currently set to a integer-typed value, or the
    value is out of range for an int4, then ``err_code`` will be set to
    ``YMMSL_bad_cast`` and the result will be invalid.

    :p character key: The name of the setting to get.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value at the given index (YMMSL_int4)
    :rtype value: integer

.. f:function:: YMMSL_Settings_get_as_int8(self, key, err_code, err_msg)

    Return the value of an integer-typed setting.

    If this setting is not currently set to a integer-typed value,
    then ``err_code`` will be set to ``YMMSL_bad_cast`` and the result
    will be invalid.

    :p character key: The name of the setting to get.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value at the given index (YMMSL_int8)
    :rtype value: integer

.. f:function:: YMMSL_Settings_get_as_real8(self, key, err_code, err_msg)

    Return the value of a real-typed setting.

    This will also work if the setting is integer-typed in which case it
    will be converted, with possible loss of precision.

    If this setting is not currently set to a real- or integer-typed value,
    then ``err_code`` will be set to ``YMMSL_bad_cast`` and the result
    will be invalid.

    :p character key: The name of the setting to get.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value at the given index (YMMSL_real8)
    :rtype value: real

.. f:subroutine:: YMMSL_Settings_get_as_real8array(self, key, value, err_code, err_msg)

    Return the value of a setting that is a 1D array of reals.

    If this setting is not currently set to a 1D array of reals-typed value,
    then ``err_code`` will be set to ``YMMSL_bad_cast`` and ``value``
    will be invalid.

    :p character key: The name of the setting to get.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value at the given index (dimension(:))
    :rtype value: YMMSL_real8

.. f:subroutine:: YMMSL_Settings_get_as_real8array2(self, key, value, err_code, err_msg)

    Return the value of a setting that is a 2D array of reals.

    If this setting is not currently set to a 2D array of reals-typed value,
    then ``err_code`` will be set to ``YMMSL_bad_cast`` and ``value``
    will be invalid.

    :p character key: The name of the setting to get.
    :p integer err_code: An error code output (optional).
    :p character err_msg: An error message output (allocatable, optional).
    :r value: The value at the given index (dimension(:,:))
    :rtype value: YMMSL_real8

.. f:function:: YMMSL_Settings_erase(self, key)

    Remove a setting from the Settings object.

    :p YMMSL_Settings self: The Settings object to modify.
    :p character key: The name of the setting to remove.
    :r removed: The number of settings removed (0 or 1; YMMSL_size)
    :rtype removed: integer

.. f:subroutine:: YMMSL_Settings_clear(self)

    Remove all settings from the Settings object.

    After calling this subroutine, the Settings object will be empty.

    :p YMMSL_Settings self: The Settings object to modify.

.. f:function:: YMMSL_Settings_key(self, i, err_code, err_msg)

    Get the i'th key in this Settings object.

    Note that any changes to the Settings object may change the order of the
    keys, so this is only stable if the Settings object is not changed.

    Parameter ``i`` must be in the range [1..N], where N is the number of items
    in the Settings object (see :f:func:`YMMSL_Settings_size`).

    :p YMMSL_Settings self: The Settings object to get a key of.
    :p YMMSL_size i: The index of the key to retrieve.
    :r key: The name of the i'th key (allocatable)
    :rtype key: character
.. _development:

Releasing
***********

MUSCLE 3 uses Git on GitHub for version management, using the `Git Flow`_
branching model. Making a release involves quite a few steps, so they're listed
here to help make the process more reliable; this information is really only
useful for the maintainers.

Make a release branch
---------------------

To start the release process, make a release branch

.. code-block:: bash

  git checkout -b release-x.y.z develop

MUSCLE 3 uses `Semantic Versioning`_, so name the new version accordingly.

Update the changelog
--------------------

Each release should have an entry in the CHANGELOG.rst describing the new
features and fixed problems. Use the git log to get a list of the changes, and
switch to the development branch:

.. code-block:: bash

  git log <your favourite options>

and then edit CHANGELOG.rst and commit.

.. code-block:: bash

  git add CHANGELOG.rst
  git commit -m 'Add version x.y.z to the change log'

Update version
--------------

Next, the version should be updated. There is a single version tag in the
``VERSION`` file in the root of the repository. On the development branch, the
version should be set to ``x.y.z-dev``, where ``x.y.z`` is the next expected
version (it's fine if that changes later, e.g. because you end up releasing
2.0.0 rather than 1.4.0).  On the release branch, it should be set to the number
of this release of course.

Check documentation
-------------------

Next, we should build the documentation to ensure that the new version number
shows up:

.. code-block:: bash

  make docs

It may give some warnings about missing references, that's a known issue and
normally harmless. Next, point your web browser to
``docs/build/html/index.html`` and verify that the documentation built
correctly. In particular, the new version number should be in the browser's
title bar as well as in the blue box on the top left of the page.

Run tests
---------

Before we make a commit, the tests should be run, and this is a good idea anyway
if we're making a release. So run ``make test`` and check that everything is in
order.

Commit the version update
-------------------------

This is the usual Git poem:

.. code-block:: bash

  git add VERSION
  git commit -m 'Set release version to x.y.z'
  git push --set-upstream origin release-x.y.z

This will trigger the Continuous Integration, so check that that's not giving
any errors while we're at it.

Fix badges
----------

The badges in the README.rst normally point to the development branch versions
of everything. For the master branch, they should point to the master version.
Note that for Codacy there is only one badge, so no change is needed.

.. code-block:: bash

  # edit README.rst
  git add README.rst
  git commit -m 'Update badges to point to master'
  git push

Merge into the master branch
----------------------------

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


Make a GitHub release
---------------------

In order to get a DOI for this release, we need to make a release on GitHub. Go
to the `MUSCLE 3 GitHub repository`_ and click 'Releases'. Select 'Draft a new
release', select the x.y.z. tag that we just uploaded, and use 'Release x.y.z'
as the title. Then copy-paste the description from the change log, convert it
from ReStructuredText to MarkDown and maybe prepend some text if there is
something interesting to mention. Optionally select 'This is a pre-release' if
it's not a final version, then publish it.

Build and release to PyPI
-------------------------

Finally, the new version needs to be built and uploaded to PyPI, so that people
can start using it. To build, use:

.. code-block:: bash

  rm -r ./build
  python3 setup.py sdist bdist_wheel

Note that we remove ``./build``, which is the build directory setuptools uses,
to ensure that we're doing a clean build, I've seen some weird mixes of versions
on occasion so it's better to be safe than sorry.

We can then check to see if everything is okay using

.. code-block:: bash

  twine check dist/muscle3-x.y.z*

and if all seems well, we can upload to PyPI:

.. code-block:: bash

  twine upload dist/muscle3-x.y.z*

Announce release
----------------

Announce the release in the usual places, so that people know it exists. There
should be a short release message listing new features and fixed bugs, and don't
forget to thank everyone who contributed!

Merge the release branch back into develop
------------------------------------------

The above concludes the release, but we need to do one more thing to be able to
continue developing. The release branch contains some changes to the change log
that we want to have back on the develop branch. So we'll merge it back in:

.. code-block:: bash

  git checkout develop
  git merge --no-commit release-x.y.z


We use --no-commit to give ourselves a chance to edit the changes before
committing them. Make sure that README.rst is taken from the develop side,
CHANGELOG.rst comes from the release branch, and VERSION is given a new number,
probably x.y.{z+1}-dev unless you have big plans. When done, commit the merge
and continue developing.


.. _`Git Flow`: http://nvie.com/posts/a-successful-git-branching-model/
.. _`Semantic Versioning`: http://www.semver.org
.. _`MUSCLE 3 GitHub repository`: https://github.com/multiscale/muscle3
Uncertainty Quantification
==========================

Only in very rare cases are the inputs and parameters of a model known exactly.
Usually, inputs and parameter values come from measurements done in the field or
in the lab, and those measurements are never exact. So every model is subject to
uncertainty. Uncertainty Quantification (UQ) provides techniques for analysing
th uncertainty in models, parameters and results. This section shows an example
of a black-box quasi-Monte Carlo uncertainty quantification of the
reaction-diffusion model we created before, implemented using MUSCLE 3.

Simulation design
-----------------

Let's start with an overview in the form of a gMMSL diagram of the set-up.

.. figure:: reaction_diffusion_qmc.png
  :scale: 40 %
  :align: center

  gMMSL diagram of the qMC reaction-diffusion simulation set-up.


At the bottom, we have the macro-micro model we saw earlier. The reaction model
is exactly identical to the previous example. For the diffusion model we have
removed the visualisation and instead send the final state on an output port
``final_state_out`` associated with the O_F operator. (It would probably have
been good to do that to begin with, using a separate visualisation or
save-to-disk component attached via a conduit.) The only other change here is
that there are now ten copies each of the ``macro`` and ``micro`` components.
There will also be ten instances of the conduits between ``macro`` and
``micro``, which MUSCLE 3 will create automatically. The instances are wired up
one-on-one, so that ``macro[i]`` connects to ``micro[i]``.

There are two new compute elements: the ``qmc`` element, which implements the
quasi-Monte Carlo UQ algorithm, and the ``rr`` element, which distributes the
ensemble evenly over a given number of model instances.

The quasi-Monte Carlo component is not a submodel, since it does not actually
model any real-world system. In MUSCLE/MMSL terminology it is a type of
component called a Proxy. We'll not get into the theoretical differences here,
but just have a practical look at how it works.

We will assume that the model parameters ``k`` and ``d`` of the reaction and
diffusion models are uncertain, and have a uniform distribution ranging from
-10% to +10% of their original values. The ``qmc`` component will generate a
quasi-random sample of the resulting box in the input parameter space, resulting
in a set of parameter values {(k, d)}. It writes these parameters to its
``parameters_out`` port.

The name of this port in the diagram ends in a pair of square brackets, which
designates the port a *vector port*. A vector port is a port that is used to
talk to multiple instances of a connected compute element. It has a number of
slots on which messages can be sent (or received, if it is a receiving port).
The number of slots is known as the *length* of a vector port. Since we are
going to run an ensemble, and we want to be able to run the ensemble members in
parallel, having a set of concurrent instances is just what we need. Of course,
we will then receive results from this same set as well, so we'll receive on a
vector port too, in this case ``states_in``.

The ``qmc.parameters_out`` port will have a length equal to the size of the
sample set it produces. This number may well be larger than the number of
instances of the model we can accomodate given the size of our compute
facilities. The ``rr`` component solves this issue by taking the messages it
receives on its ``front_in`` vector port, and distributing them evenly amongst
the slots of its ``back_out`` vector port in a round-robin fashion (hence the
name). The length of ``front_in`` matches that of ``qmc.parameters_out``, while
the length of ``back_out`` matches the number of instances of ``macro``. The
returning final states will be mapped back accordingly.

Next, ``rr`` needs to be connected to the model. This presents a problem: ``rr``
sends out sets of parameters, but ``macro`` has no F_INIT port to receive them.
It gets its parameter values, via MUSCLE 3, from the central configuration. So
we need a trick, and MUSCLE 3 provides one in the form of the
``muscle_settings_in`` port. This is a special F_INIT port that each MUSCLE 3
compute element automatically has. It can be connected to a port on a component
that sends ``Settings`` objects. MUSCLE 3 will automatically receive these
messages, and overlay the received settings on top of the base settings from the
central configuration. When the receiving submodel then asks MUSCLE 3 for a
setting, MUSCLE 3 will look at the overlay settings first. If the setting is
not found there, then MUSCLE 3 will fall back to the central base configuration.

So, now our diffusion model will ask for the value of ``d`` as it did before,
but this time, it will actually come from the values sent by ``qmc``, and it
will be different for each instance of the diffusion model. In a way, each
instance lives in its own universe, with its own universal constants. It is
unaware of this however; to the diffusion model instance everything looks just
the same as in the previous example when it was the only one around.

The one remaining piece of the puzzle is that the universe will be automatically
extended: overlay parameters that were sent from the ``rr`` component to a
macro-model instance will automatically be passed on to the corresponding
micro-model instance. In this way, each ``k`` parameter value arrives at the
reaction model that needs it.

Implementation
--------------

The full Python program implementing this looks like this:

.. literalinclude:: examples/python/reaction_diffusion_qmc.py
  :caption: ``docs/source/examples/python/reaction_diffusion_qmc.py``
  :language: python


The reaction model's implementation is completely unchanged, and the diffusion
model has only been made a bit more modular by removing the visualisation and
sending the final state on an additional O_F port, as mentioned above. The
``qmc``
and ``rr`` components are new however, and deserve a detailed explanation.

The quasi-Monte Carlo element
-----------------------------

The purpose of the qMC component is to generate a set of quasi-random pairs of
parameters ``(k, d)``, and pass these to the individual members of the ensemble
of models we will run. It receives back the final state of the simulation of
each ensemble member, and then calculates the ensemble mean final state.

.. code-block:: python

  instance = Instance({
          Operator.O_I: ['parameters_out[]'],
          Operator.S: ['states_in[]']})


The qMC component does not require any input messages to initialise, and it does
not produce a final result as a message (it makes a plot instead). It just sends
out a single set of parameter pairs, and receives back a single set of final
states. In terms of communication with the outside world, it therefore works
similarly to a submodel with only O_I and S ports and exactly one state update
step, so that each of those operators is run once. So that is how we describe it
to MUSCLE 3. (It may be better to have some Proxy-specific operator names here
in the future, opinions welcome!) We will send parameter sets on vector port
``parameters_out``, and receive final states on vector port ``states_in``.

Next, we enter the reuse loop as before, except that we pass ``False`` as an
argument, which will be explained shortly. We read and check settings in
F_INIT as before, and calculate a set of parameter values using a quasi-random
Sobol sequence.

.. code-block:: python

  # configure output port
  if not instance.is_resizable('parameters_out'):
      instance.error_shutdown('This component needs a resizable'
                          ' parameters_out port, but it is connected to'
                          ' something that cannot be resized. Maybe try'
                          ' adding a load balancer.')
      exit(1)

  instance.set_port_length('parameters_out', n_samples)


Next, we need to configure our output vector port. Since we will have
``n_samples`` sets of parameters, we will resize the port to that length. We do
need to check whether the port is resizable first, because that may or may not
be the case depending on what is attached to it. In order to preserve modularity
and reusability of individual compute elements, MUSCLE 3 tries to tell the
compute element as little as possible about what is on the other side of a port,
but you can ask it whether the port has a fixed size or not.

So that is what we do here, and we generate an error message if the port's
length is fixed. We use the function
:meth:`libmuscle.Instance.error_shutdown()` for this.  This function will log
the error, and then tell the rest of the simulation that we are shutting down.
Doing this instead of raising an exception reduces the chance that any part of
the simulation will sit around waiting forever for a message we will never
send, and the log will show the origin of the problem for easier debugging. In
this case, the port will be resizable and it will work as intended.

.. code-block:: python

  for sample in range(n_samples):
      uq_parameters = Settings({
          'd': ds[sample],
          'k': ks[sample]})
      msg = Message(0.0, None, uq_parameters)
      instance.send('parameters_out', msg, sample)

Since we only run our O_I and S once, we do not have a state update loop that
advances the simulation time. We do however need to loop through all of our
samples, and send a message on ``parameters_out`` for each one.

First, we create a ``Settings`` object that contains the parameters we are going
to send. This is the same ``Settings`` class we used before to define the model
configuration. Here, we only have ``d`` and ``k`` in there, since those are the
only ones we need to set per ensemble member; the rest is identical and defined
in the central configuration.

Next, we create a :class:`libmuscle.Message` object to send. Since our models
will start at time 0, we'll set that as the timestamp, and since we're only
running them once each, the next timestamp is ``None``. For the data, we send
the ``Settings`` object. (MUSCLE 3 contains special support for sending
``Settings`` objects, since being objects they're not normally
MessagePack-serialisable.)

We then send our message as normal, except that we pass an extra argument, the
*slot number*. Vector ports connect to sets of instances, and the slot number
selects the exact instance to send the message to. They're zero-based, so if the
length of the vector port is 10, then the valid slot range is [0..9]. We resized
our port to the number of samples before, so each sample has a corresponding
slot for us to send on.

.. code-block:: python

  for sample in range(n_samples):
      msg = instance.receive_with_settings('states_in', sample)


When the reaction-diffusion models are done, they will send their final states
back to us, so we need to receive those now. This is effectively our S operator.
For each sample, we receive the result on our ``states_in`` port, passing the
sample number as the slot to receive on. We're using a slightly different
receive function here. Rather than :meth:`libmuscle.Instance.receive`, we call
:meth:`libmuscle.Instance.receive_with_settings`.  The difference has to do
with the settings overlays.

Recall that each compute element instance has a settings overlay, which can be
set through the ``muscle_settings_in`` port and is automatically propagated to
corresponding instances of other compute elements. This propagation is done by
sending the overlay along with any messages that are sent.

Since the ``macro`` and ``micro`` instances are connected one-on-one, each pair
lives in its own universe with its own settings, unaware of the other pairs and
their settings.  Once the reaction-diffusion simulation is done however, the
``macro`` instances all send their result to a single ``qmc`` instance (via
``rr``, which is transparent in this respect). Thus, the universes meet, and the
question arises what the settings overlay for ``qmc`` should look like.

As there is no general answer to this, MUSCLE 3 cannot automatically propagate
the overlay from the different ``macro`` instances to ``qmc``, and it will give
an error message if you try to have it do this by receiving as usual with
:meth:`libmuscle.Instance.receive`. The
:meth:`libmuscle.Instance.receive_with_settings` function solves this problem,
and is in a way the counterpart of sending a message to ``muscle_settings_in``.
It will not try to merge the incoming settings overlay into the overlay for
``qmc``, but simply return it as the ``settings`` attribute of the received
message. It is then up to the receiver to decide what to do with it.

In this case, we ignore the settings, concatenate all the received states
together, plot them, and calculate the mean final state. You will probably want
to do a more complex analysis (for which the settings may be very useful!), or
save the raw data to disk for later processing.

The round-robin load balancer
-----------------------------

The round-robin load balancer has the job to sit between the ``qmc`` element and
the macro model, and distribute the many parameter sets ``qmc`` produces over a
limited number of macro model instances. It has a front side, which connects to
``qmc``, and a back side, which connects to ``macro``.

.. code-block:: python

  while instance.reuse_instance(False):
      # F_INIT
      started = 0     # number started and index of next to start
      done = 0        # number done and index of next to return

      num_calls = instance.get_port_length('front_in')
      num_workers = instance.get_port_length('back_out')

      instance.set_port_length('front_out', num_calls)
      while done < num_calls:
          while started - done < num_workers and started < num_calls:
              msg = instance.receive_with_settings('front_in', started)
              instance.send('back_out', msg, started % num_workers)
              started += 1
          msg = instance.receive_with_settings('back_in', done % num_workers)
          instance.send('front_out', msg, done)
          done += 1


In order to distribute the messages correctly, we first need to determine the
number of messages we'll receive from ``qmc``, as well as how many ``macro``
instances we have. We can determine this from the lengths of the corresponding
vector ports. Since ``front_out`` is connected to another vector port, it does
not have an intrinsic size, and we need to set the size explicitly to match
``front_in``.

Next, we process messages, reading them from ``front_in`` and forwarding them to
``back_out``, always sending each ``macro`` instance one message at a time, and
making sure that for each message we send on ``back_out``, we receive one on
``back_in``. (We could actually send multiple messages on the same slot before
receiving a result, they'll be queued up and processed in order.)

We use :meth:`libmuscle.Instance.receive_with_settings` everywhere, in order to
correctly pass on any settings overlays. Since we are using
:meth:`libmuscle.Instance.receive_with_settings` on an F_INIT port, we passed
``False`` to :meth:`libmuscle.Instance.reuse_instance`. It is a technical
requirement of MUSCLE 3 to do this, and MUSCLE will give an error message if
you call :meth:`libmuscle.Instance.receive_with_settings` without having passed
``False`` to :meth:`libmuscle.Instance.reuse_instance`. (There's just no other
way to implement this, or rather, all other options can lead to potentially
difficult-to-debug situations, while this can be checked and a clear error
message shown if it goes wrong. So we chose this as the preferable option.)

Discussion
----------

There are a few interesting things to remark here. First, the original model
code is essentially unchanged. All the changes needed to do Uncertainty
Quantification have been at the level of the couplings between the models.
MUSCLE encourages modularity and reuse of models, and this example shows the
benefits of this approach.

Second, we can see the beginnings of a library of reusable components. The
round-robin load balancer shown here is completely model-agnostic, and in a
future version of MUSCLE 3 will become a built-in standard component for general
use. There are other such components that can be made, such as the duplication
mapper that MUSCLE 2 already has. With a small extension to MUSCLE, the ``qmc``
component here can also be made generic and model-agnostic. While some
helper components will likely remain model-specific (such as scale bridges and
data converters), we expect that modeling complex systems and performing UQ on
them can be made significantly simpler using this approach.

The load balancer component described here uses a round-robin algorithm to
distribute work. This works well in this case, with all model runs taking
approximately the same amount of compute time. In general however, a more
flexible algorithm is desirable, which would require another small extension to
MUSCLE 3. We plan to add this in a future version.

Examples in C++ and Fortran
---------------------------

MUSCLE 3 comes with C++ and Fortran versions of this Uncertainty Quantification
use case. See the C++ or Fortran section for how to run the examples, and the
``reaction_diffusion_mc_cpp.sh`` and ``reaction_diffusion_mc_fortran.sh``
scripts respectively. The source code for the components may be found in
``docs/source/examples/cpp`` and ``docs/source/examples/fortran``.
=================
Development tools
=================

MUSCLE 3 is a fairly complex piece of software, being a coupling library
designed to link together independent bits of software written in different
languages. To help us develop MUSCLE 3, we use a variety of tools that check our
work and build software and documentation where required. This page describes
the development tools we use, what they do, and how they are configured.


Building and installation
=========================

The main language for MUSCLE 3 is Python 3. To help install and publish the
Python parts of MUSCLE 3, we use `setuptools`_. Setuptools reads information
from the file `setup.py` in the root of the repository. This file contains a
description of MUSCLE 3 which is used on PyPI (the standard online Python
package repository), information used by PIP when installing MUSCLE 3, and
information about required libraries.

The dependencies in `setup.py` are specified with an exact version number,
rather than an at-least version number. This way, we can be sure that no
incompatible changes are introduced that break MUSCLE 3. However, it also means
that the dependencies get out of date, as new versions are released. To keep
up-to-date on this, we use `Requires.io`_, a web service that scans our
`setup.py` and checks our versions against the latest releases. It produces a
button, which we embed in our README.rst, that shows whether we are up-to-date.

We will most likely have outdated dependencies quite often, as rapid updates are
a fact of life in today's Open Source world. This is not an immediate problem, as
the older versions will still be available and MUSCLE 3 will work with them, but
we should update regularly to stay current and avoid potential security issues.


Quality Control
===============

MUSCLE 3 is infrastructure used for building infrastructure (models) for doing
scientific experiments. As such, it is important that it is reliable, easy to
use, and easy to develop. Furthermore, it will frequently be used by
inexperienced programmers, who will look at MUSCLE 3 as an example of how to
write good code. It is therefore important that MUSCLE 3 is a well-implemented,
high quality program. To make sure that this happens, we use static checking
(linting), unit tests, and integration tests.

Static checking
---------------

The quality control system uses several tools, which are called from the central
Makefile. In our standard development environment, the command `make test` in
the top directory will run them all.

Python has a standard coding standard known as PEP8, which specifies how Python
code should be formatted. Consistently writing our code according to PEP8 makes
it easier to read, and makes it easier for new developers to get started. We use
the `PEP8 plug-in`_ for Pytest (see below) to check the formatting.

Although Python is a dynamically-typed language, it supports type annotations,
and using those for static type checking. While writing code with type
annotations sometimes takes a bit more thinking than doing without them, it
increases the reliability of the program a lot, as type checking can detect
mistakes that unit tests don't. So we use type annotations, and use `mypy`_ to
check them as part of the set of quality control software.

Finally, we use `Codacy`_, which is an online tool that runs various static
checking tools against our repository. It is run automatically for the main
branches, as well as for every pull request. Codacy tends to give false
positives fairly often, but the settings can be tuned and specific warnings can
be set to be ignored. So having some Codacy warnings on a pull request is not
necessarily a problem, but they should be reviewed and either fixed or disabled
as appropriate.

Tests
-----

MUSCLE 3 has a suite of tests, both unit tests (which test small bits of the
code one piece at a time) and integration tests (which test whether it all goes
together correctly). We run these tests, measure test coverage, and have
continuous integration just in case we forget to run the tests manually.

The main test framework for the Python part of the code is `PyTest`_. This
Python library provides a library with functions that make it easier to write
test cases, and a runner that automatically detects tests and runs them. Tests
are located in `test/` subdirectories in the code, and in the `integration_test`
directory in the root of the repository. In various places in the code, files
named `conftest.py` are located, which contain test fixtures. These files are
picked up automatically by PyTest. Having such a file in the root directory
causes PyTest to add that path to the PYTHONPATH, so that Python can find the
module under test if you import it from the test case. This doesn't work for the
integration tests, so there's a separate `include_libmuscle.py` there that sets
the path correctly.

We use a code coverage plug-in for PyTest, which measures which lines of the
code are executed during testing, and more importantly, which lines aren't.
Ideally, all code is covered, but in practice there will be things that are
difficult to test and not important enough to warrant testing. 90% coverage is a
nice (but sometimes challenging) target to aim for. This plug-in uses the
Python coverage library, which is configured in `.coveragerc`, where e.g.
generated code is excluded from testing.

Tests can be run locally using `make test`, but are also run in the cloud via
`Travis-CI`_ on every push to the Github repository. After Travis has run the
tests, it reports the code coverage to `Codacy`_, which integrates it with the
linting results into its dashboard. It also reports back to GitHub, adding an
icon to e.g. a pull request to signal whether the tests have passed or not.


Documentation
=============

The MUSCLE 3 documentation lives in the repository together with the code. This
is a much preferable arrangement to a separate wiki, as previous MUSCLEs had,
because wikis tend to get lost.

The source for the documentation is in `docs/source`, and it is written in
`ReStructuredText`_. It is processed into HTML by the `Sphinx`_ program, which
uses the `sphinx-autodoc`_ plug-in for extracting docstrings from the Python
source code and converting that into API documentation. Sphinx is configured
via `docs/source/conf.py`. This file contains a hook to automatically run
`sphinx-autodoc` whenever Sphinx is run, so that is does not have to be run
separately first.

The `ReadTheDocs`_ online service is configured to automatically render the
documentation and display it online. This uses Sphinx, but does not support
`sphinx-autodoc`, which is why we use the aforementioned hook.


.. _`setuptools`: https://setuptools.readthedocs.io
.. _`Requires.io`: https://requires.io/
.. _`PEP8 plug-in`: https://pypi.python.org/pypi/pytest-pep8
.. _`mypy`: https://mypy.readthedocs.io
.. _`Codacy`: https://support.codacy.com
.. _`PyTest`: https://pytest.org
.. _`Travis-CI`: https://docs.travis-ci.com
.. _`ReStructuredText`: http://docutils.sourceforge.net/rst.html
.. _`Sphinx`: http://www.sphinx-doc.org
.. _`sphinx-autodoc`: http://www.sphinx-doc.org/en/master/ext/autodoc.html
.. _`ReadTheDocs`: https://readthedocs.org
=================
MUSCLE 3 Road Map
=================

1. Implement logging server side

    * MMP logging call
    * Manager Logger component
    * Unit tests and Continuous Integration
    * Documentation build
    * Linting

2. Implement logging client side

    * LibMuscle logging API
    * Initial MMP Client
    
==========================
MUSCLE 3 Use Case Road Map
==========================

These use cases are probably far too few, but I hope Mark and Lourens can add more of them, or change them around to make the progression more natural?

1. MUSCLE2 App with Two coupled submodels, running each within the same work station.
2. MUSCLE2 App with Three coupled submodels, running each within the same work station.
3. MUSCLE2 App with Two coupled submodels, running on different nodes (using TCP communication)
4. MUSCLE2 App with Three coupled submodels, running on different nodes (using TCP communication)
5. MUSCLE2 App with Two coupled submodels, running on different nodes (using MPI communication)
6. MUSCLE2 App with Three coupled submodels, running on different sites (using TCP/MPWide communication)
MUSCLE and MPI
==============

MPI, the Message Passing Interface, is the most widely used communications
middleware for High-Performance Computing. It's very low-level, but also usually
very well optimised, and perhaps more importantly, the standard. Many existing
HPC applications use MPI for communication, either directly or through a
library.

Chances are, therefore, that you will run into a submodel that uses MPI and that
you want to connect to MUSCLE 3. MUSCLE 3 supports this (only in C++ and
Fortran for now), but as always with MPI, some extra care needs to be taken
because things run in parallel.

MPI follows the Single-Program Multiple Data paradigm, which means that there
are many copies of your program running simultaneously doing the same thing,
except that some of the data they're processing are different for each copy, and
every once in a while the copies communicate.

The MUSCLE parts of the simulation follow this: each process makes an Instance,
each process runs through the reuse loop, and each process calls the send and
receive functions at the same time. Two exceptions are that configuring ports
may only be done in the root process, and getting settings can be done anytime
and anywhere.

This section shows how to use MUSCLE 3 with MPI from C++ and Fortran, based on
the same reaction-diffusion model given in the C++ and Fortran sections, but now
with a reaction model that uses MPI. It will help if you've read the Python
tutorial and the C++ or Fortran section first, and some experience with MPI or
having browsed through a tutorial will help.

`The source code for the example in this section is here
<https://github.com/multiscale/muscle3/tree/master/docs/source/examples/cpp>`_.
You can also go to ``docs/source/examples/cpp`` in the source directory (that's
easiest, and has some handy scripts and a Makefile), or copy-paste the code from
here. For Fortran, there is also `an online version
<https://github.com/multiscale/muscle3/tree/master/docs/source/examples/fortran>`_
and a version in your local source tree at ``docs/source/examples/fortran``.

Building and running the example
--------------------------------

If you've just built and installed the C++ version of libmuscle, then you're all
set to build the examples. Go into the ``docs/source/examples`` subdirectory of
the source directory, and type ``MUSCLE3_HOME=<PREFIX> MUSCLE_ENABLE_MPI=1 make
cpp python`` and the examples will be compiled and linked against your new
libmuscle.  For Fortran, use ``MUSCLE3_HOME=<PREFIX> MUSCLE_ENABLE_MPI=1 make
fortran python`` instead.

In order to run the examples, you'll need the Python version of MUSCLE
3 installed, because you need the MUSCLE Manager, which is why ``python`` needs
to be added to the Make command. The build system will set up a virtual
environment in ``docs/source/examples/python/build/venv`` with everything needed
to run the examples.

You can then run the examples using the provided scripts in
``docs/source/examples``:

.. code-block:: bash

  ~/muscle3_source/docs/source/examples$ MUSCLE3_HOME=~/muscle3 ./reaction_diffusion_mpi_cpp.sh
  ~/muscle3_source/docs/source/examples$ MUSCLE3_HOME=~/muscle3 ./reaction_diffusion_mpi_fortran.sh


Note that these examples are compiled with slightly different options compared
to the non-MPI examples, and you should use those as well when compiling your
own MPI-based submodels. See the Installing section for details.


MPI Reaction model
------------------

Here is the C++-and-MPI version of the reaction model:

.. literalinclude:: examples/cpp/reaction_mpi.cpp
  :caption: ``docs/source/examples/cpp/reaction_mpi.cpp``
  :language: cpp

And this is the Fortran-and-MPI version of the reaction model:

.. literalinclude:: examples/fortran/reaction_mpi.f03
  :caption: ``docs/source/examples/fortran/reaction_mpi.f03``
  :language: fortran


We'll go through it top to bottom, one piece at a time, focusing on the
differences with the plain C++ and Fortran versions.

Headers and modules
```````````````````

.. code-block:: cpp
    :caption: C++

    #include <mpi.h>

    #include <libmuscle/libmuscle.hpp>
    #include <ymmsl/ymmsl.hpp>


In order to use MPI, we need to include the MPI header. The headers for
libmuscle are unchanged. Remember that we're compiling with
``MUSCLE_ENABLE_MPI`` defined however, and this changes the API slightly as
shown below.

.. code-block:: fortran
    :caption: Fortran

    use mpi
    use ymmsl
    use libmuscle_mpi


In Fortran, we use the ``mpi`` module to be able to make MPI calls, and
``libmuscle_mpi`` to get the MPI-enabled version of the MUSCLE 3 API.


.. code-block:: cpp
    :caption: C++

    void reaction(int argc, char * argv[]) {
        const int root_rank = 0;
        int rank, num_ranks;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);


.. code-block:: fortran
    :caption: Fortran

    integer, parameter :: root_rank = 0
    integer :: rank, num_ranks, ierr

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, num_ranks, ierr)


In MPI, all parallel processes are equal, but from time to time something must
be done from, to, or in only one of the processes, and so one process needs to
be declared special. This is called the root process, and it's usually the one
with rank zero. That's what we do here as well, and then we ask MPI for the rank
of the current process, and the total number of processes.

Creating an Instance
````````````````````

.. code-block:: cpp
    :caption: C++

    Instance instance(argc, argv, {
            {Operator::F_INIT, {"initial_state"}},  // list of double
            {Operator::O_F, {"final_state"}}},      // list of double
            MPI_COMM_WORLD, root_rank);


.. code-block:: fortran
    :caption: Fortran

    ports = LIBMUSCLE_PortsDescription_create()
    call LIBMUSCLE_PortsDescription_add(ports, YMMSL_Operator_F_INIT, 'initial_state')
    call LIBMUSCLE_PortsDescription_add(ports, YMMSL_Operator_O_F, 'final_state')
    instance = LIBMUSCLE_Instance_create(ports, MPI_COMM_WORLD, root_rank)
    call LIBMUSCLE_PortsDescription_free(ports)


When MUSCLE 3 is used in MPI mode, the Instance constructor takes two extra
arguments: an MPI communicator, and the rank of the root process. These default
to MPI_COMM_WORLD and 0, respectively. MUSCLE 3 will create a copy of the given
communicator for internal use; the one you pass must contain all processes in
your submodel. MUSCLE 3 will do all communication with other submodels from the
process with the given root rank, and that process is special when sending and
receiving as described below.

Reuse loop and settings
```````````````````````

.. code-block:: cpp
    :caption: C++

    while (instance.reuse_instance()) {
        // F_INIT
        double t_max = instance.get_setting_as<double>("t_max");
        double dt = instance.get_setting_as<double>("dt");
        double k = instance.get_setting_as<double>("k");


.. code-block:: fortran
    :caption: Fortran

    do while (LIBMUSCLE_Instance_reuse_instance(instance))
        ! F_INIT
        t_max = LIBMUSCLE_Instance_get_setting_as_real8(instance, 't_max')
        dt = LIBMUSCLE_Instance_get_setting_as_real8(instance, 'dt')
        k = LIBMUSCLE_Instance_get_setting_as_real8(instance, 'k')


This part is unchanged from the non-MPI versions. This means that all processes
enter the reuse loop together. With MPI enabled, ``reuse_instance()`` is
effectively a collective operation, so it must be called in all processes
simultaneously.

Settings may be obtained at any point and in any MPI process. Getting a setting
value does not require any MPI activity, so it can be done in any way you like
as long as it's within the reuse loop.

Distributed state
`````````````````

.. code-block:: cpp
    :caption: C++

    std::vector<double> U, U_all;
    int U_size;
    double t_cur, t_end;


.. code-block:: Fortran
    :caption: Fortran

    real (selected_real_kind(15)) :: t_cur, t_max, t_end, dt, k
    integer :: i, U_size, U_all_size
    real (selected_real_kind(15)), dimension(:), allocatable :: U, U_all


The state of the reaction model is a vector of doubles. In this MPI version,
we'll be dividing that state among the processes, so that they can each process
a part of it. The ``U`` variable exists in each process and contains its piece
of the total state. ``U_all`` contains the complete state in the root process,
and it's empty in the rest of the processes. ``U_size`` is the size of ``U``,
and then we have ``t_cur`` and ``t_end`` in all processes, containing
respectively the current simulation time, and the time at which to end the
current run.

Receiving messages
``````````````````

Next, it's time time to receive the initial state:

.. code-block:: cpp
    :caption: C++

    auto msg = instance.receive("initial_state");


.. code-block:: fortran
    :caption: Fortran

    rmsg = LIBMUSCLE_Instance_receive(instance, 'initial_state')


The ``receive()`` function is another collective operation, so it must be called
in all processes. All processes will block until a message is received. The
message will be returned in the root process (as designated when creating the
Instance), in all other processes, an empty dummy message is returned.

This may seem a bit odd (why not just receive only in the root process, or
return the received message in all processes), but it has a good reason. When
MPI needs to wait for something, it goes into a loop that continuously checks
whether the event has occurred (a spinloop). This keeps the CPU core it's
running on fully loaded. If you received only in the root process, and then
used a broadcast operation to distribute data, the root process would block
(without using CPU) on the network receive call until the other model was done,
while all other processes would spin frantically, keeping almost all your cores
occupied.

Since macro-micro models alternate execution, it's often nice if you can put
them on the same cores, but that only works if the waiting model doesn't
load them when it's not running. MUSCLE solves this problem using a TCP-based
barrier. You call ``receive()`` in all processes, and they will all block on a
network receive call, without using CPU, until a message is received and they
can continue.

The message is then sent only to the root, because MUSCLE does not know whether
it needs to be broadcast to all processes, or distributed somehow. You'll have
to do that part yourself.


.. code-block:: cpp
    :caption: C++

    if (rank == root_rank) {
        DataConstRef data(msg.data());
        U_all.resize(data.size());
        for (int i = 0; i < data.size(); ++i)
            U_all[i] = data[i].as<double>();

        t_cur = msg.timestamp();
        t_end = t_cur + t_max;

        U_size = U_all.size() / num_ranks;
        if (U_size * num_ranks != U_all.size()) {
            instance.error_shutdown("State does not divide evenly");
            throw std::runtime_error("State does not divide evenly");
        }
    }


.. code-block:: fortran
    :caption: Fortran

    if (rank == root_rank) then
        rdata = LIBMUSCLE_Message_get_data(rmsg)
        U_all_size = LIBMUSCLE_DataConstRef_size(rdata)
        allocate (U_all(U_all_size))
        do i = 1, U_all_size
            item = LIBMUSCLE_DataConstRef_get_item(rdata, int(i, LIBMUSCLE_size))
            U_all(i) = LIBMUSCLE_DataConstRef_as_real8(item)
            call LIBMUSCLE_DataConstRef_free(item)
        end do
        call LIBMUSCLE_DataConstRef_free(rdata)

        t_cur = LIBMUSCLE_Message_timestamp(rmsg)
        t_end = LIBMUSCLE_Message_timestamp(rmsg) + t_max
        call LIBMUSCLE_Message_free(rmsg)

        U_size = U_all_size / num_ranks
        if (U_size * num_ranks /= U_all_size) then
            call LIBMUSCLE_Instance_error_shutdown(instance, 'State does not divide evenly')
            print *, 'State does not divide evenly'
            stop
        end if
    end if


Here, we start in the root process only by unpacking the received list into
``U_all``, and by setting the current and final times. We also calculate the
size of the per-process portion of the state, and check that it divides evenly.
If it doesn't, we tell MUSCLE that we've encountered an error (you can do this
either from the root process, or from all processes simultaneously), and raise
an exception. This example comples configured with a 100-cell long grid and two
MPI processes, so it'll work well.

Next, we distribute the received information among the processes:

.. code-block:: cpp
    :caption: C++

    MPI_Bcast(&U_size, 1, MPI_INT, root_rank, MPI_COMM_WORLD);
    U.resize(U_size);
    MPI_Scatter(U_all.data(), U_size, MPI_DOUBLE,
                U.data(), U_size, MPI_DOUBLE,
                root_rank, MPI_COMM_WORLD);

    MPI_Bcast(&t_cur, 1, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);
    MPI_Bcast(&t_end, 1, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);


.. code-block:: fortran
    :caption: Fortran

    call MPI_Bcast(U_size, 1, MPI_INT, root_rank, MPI_COMM_WORLD, ierr)
    allocate (U(U_size))
    call MPI_Scatter(U_all, U_size, MPI_DOUBLE,  &
                     U, U_size, MPI_DOUBLE,      &
                     root_rank, MPI_COMM_WORLD, ierr)

    call MPI_Bcast(t_cur, 1, MPI_DOUBLE, root_rank, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(t_end, 1, MPI_DOUBLE, root_rank, MPI_COMM_WORLD, ierr)


This is fairly standard MPI code. First, we broadcast the size of the
per-process state. This is necessary, because it's derived from the received
message, and we only received it in the root process. Next, we make some space
in ``U`` for the new state information, and then we scatter ``U_all`` to the
``U`` vectors in the processes. The current time and end time are also derived
from the received message, so they need to be broadcast as well.

Next is the state update loop, which is completely unchanged. Each process
processes its part of the state, and since the reaction in a grid cell is not
affected by anything outside of that grid cell, we don't need to do any
communication. We could of course if this were e.g. the diffusion model, and we
needed to send data to the neighbours. That's not MUSCLE-related though, so we
refer to an MPI tutorial for more information on how to do that.

Finally, once we're done iterating, we need to send out the final state:

.. code-block:: cpp
    :caption: C++

    // O_F
    MPI_Gather(U.data(), U_size, MPI_DOUBLE,
               U_all.data(), U_size, MPI_DOUBLE,
               root_rank, MPI_COMM_WORLD);

    if (rank == 0) {
        auto result = Data::nils(U_all.size());
        for (int i = 0; i < U_all.size(); ++i)
            result[i] = U_all[i];
        instance.send("final_state", Message(t_cur, result));
    }


.. code-block:: fortran
    :caption: Fortran

    ! O_F
    call MPI_Gather(U, U_size, MPI_DOUBLE,       &
                    U_all, U_size, MPI_DOUBLE,   &
                    root_rank, MPI_COMM_WORLD, ierr)

    if (rank == root_rank) then
        sdata = LIBMUSCLE_Data_create_nils(int(U_all_size, LIBMUSCLE_size))
        do i = 1, U_all_size
            call LIBMUSCLE_Data_set_item(sdata, int(i, LIBMUSCLE_size), U_all(i))
        end do

        smsg = LIBMUSCLE_Message_create(t_cur, sdata)
        call LIBMUSCLE_Instance_send(instance, 'final_state', smsg)

        call LIBMUSCLE_Message_free(smsg)
        call LIBMUSCLE_Data_free(sdata)
        deallocate (U_all)
    end if


The ``send()`` function can be called either from the root process only, or from
all processes. In the latter case, it will simply do something only in the root
process, and return immediately in all other processes. This sometimes gives
cleaner code though, which is why it's an option.

In this case, we're going to do the sending only in the root process. We first
use a gather operation in all processes to collect the data from the local ``U``
variables into the root's ``U_all``. Then, we convert it into a ``Data`` object
and send it as before, but only in the root process.


MPI requires initialisation and finalisation, which for C++ we do in the main
function:

.. code-block:: cpp

    int main(int argc, char * argv[]) {
        MPI_Init(&argc, &argv);
        reaction(argc, argv);
        MPI_Finalize();
        return EXIT_SUCCESS;
    }


In the Fortran version of the example, it's in the main program.

Note that ``MPI_Init()`` must have been called before a MUSCLE ``Instance`` is
created, since the ``Instance`` constructor will make MPI calls. The ordering of
``MPI_Finalize()`` is less strict, but for symmetry we free the ``Instance``
first in the Fortran version. In C++, this is done automatically at the end of
the ``reaction()`` function.

Distributed execution
=====================

In the previous section, we created a simple macro-micro multiscale model with
MUSCLE 3, and ran it as a single Python script. This section briefly explains
how to go from there to a distributed simulation, possibly on multiple nodes.

Note that distributed simulations are not as easy to use as we would like them
to be yet. All the necessary pieces are in place however. Below, we explain how
they work. If you want to run the example below using a premade script, then you
can go to the ``docs/examples/`` directory in the source code, and set up and
run the Python example like this:

.. code-block:: bash

  docs/examples$ make python
  docs/examples$ ./reaction_diffusion_python.sh


Below, we explain how it works and how to run it by hand.


Making separate instances
-------------------------

Previously, we started a simulation by starting a single Python script that
contained all of the implementations and also the model configuration. For
distributed running, we want to run the instances separately, so our script
needs to be split up. Doing that is very simple: just copy the model function
(and anything it depends on) to a new file, and add a main clause that runs it.
For the reaction model of the previous example, that looks like this:

.. literalinclude:: examples/python/reaction.py
  :caption: ``docs/source/examples/python/reaction.py``
  :language: python


Note that the code is exactly the same, we've just removed everything related to
the diffusion model and to the coupling between them. We can do the same to the
diffusion model:

.. literalinclude:: examples/python/diffusion.py
  :caption: ``docs/source/examples/python/diffusion.py``
  :language: python


Again, it's exactly the same code as before, just split off into a separate
file.

yMMSL files
-----------

In a distributed set-up, the manager and each instance run as separate programs,
communicating via the network. To make this work, the manager is started first,
with a yMMSL file that describes the simulation. The yMMSL file corresponding to
our previous example looks like this:

.. literalinclude:: examples/reaction_diffusion.ymmsl
  :caption: ``docs/source/examples/reaction_diffusion.ymmsl``
  :language: yaml


As you can see, this looks very much like the object representation. You can
load a yMMSL file to objects using ``ymmsl.load`` and save it back using
``ymmsl.save``.

Starting the manager
--------------------

Every MUSCLE 3 simulation needs a manager. The manager helps the instances find
each other, and distributes settings to them. To start the manager, run

.. code-block:: bash

  (venv)$ muscle_manager reaction_diffusion.ymmsl


in an environment in which you've installed MUSCLE 3. The manager will start,
and will print its location on the standard output. This is the network address
on which it listens for connections from instances.

Starting instances
------------------

Starting an instance works just like starting any other Python script:

.. code-block:: bash

  (venv)$ python ./reaction.py --muscle-instance=micro


In the original script, we made a dictionary that mapped compute element names
to Python functions, which MUSCLE 3 used to start the submodel instances. Here,
we need to pass the same information. Instead of a function name, we start a
particular script, and we pass the instance name using the
``--muscle-instance=`` argument. Note that the argument parser is not very
flexible, so it needs to be written exactly like that, and of course the
instance name needs to match the compute element name in the yMMSL file passed
to the manager.

When the instance starts, it will register with the manager. Since we didn't
pass a location using the ``--muscle-manager=<host:port>`` option, MUSCLE 3
tried to connect to the default location (localhost:9000), which the manager
was conveniently listening on. The ``micro`` instance registered, and if you
have a look at the ``muscle3_manager.log`` file produced by the manager (in the
directory in which you started it), you should find a record of this there.

Now, we need to start the second submodel in the same way:

.. code-block:: bash

  (venv)$ python ./diffusion.py --muscle-instance=macro


and then the simulation will run. Once it's done, the instances disconnect and
the manager shuts itself down automatically.
Current status
==============

The current release includes most of the necessary features to run multiscale
simulations in both non-distributed and distributed mode, in Python 3 and C++.
This makes this release suitable for prototyping and learning to work with
MUSCLE 3 and multiscale modeling in general, and for production work using
Python and C++ models.

In particular, MUSCLE 3 currently provides the following features:

Coupling different submodel instances
-------------------------------------

MUSCLE connects simultaneously running instances of different submodels together
based on an MMSL description of how to do so. Single instances can be connected
to each other, or to sets of instances via so-called vector ports. MUSCLE allows
peer instances to find each other, and then transmit messages between themselves
in a peer-to-peer fashion. Communication between the models therefore does not
suffer from a bottleneck created by a central component.

Spatial and temporal scale separation and overlap
-------------------------------------------------

For a macro-micro coupled model with time scale separation, the micro-model
needs to be run once for each time step of the macro-model. With MUSCLE, models
can be run repeatedly to make this possible. If the macro-model does not have a
fixed number of timesteps (e.g. because it runs until convergence, and it is not
known in advance when that occurs), then the number of required runs of the
micro-model is unknown on beforehand as well. MUSCLE will automatically
coordinate the model's execution in order to ensure that exactly the right
number of runs happen and no deadlocks occur.

For time scale overlap, the models run simultaneously and exchange information
via MUSCLE while running, but interpolation may be needed to synchronise the
data that is exchanged. Messages passed between models via MUSCLE carry
simulation-time timestamps that allow you to implement a suitable interpolation
component.

For spatial scale separation, usually many instances of the micro-model need to
be run to cover the macro-scale domain, and they need to be linked to the
macro-model. This too can be described with the MMSL and implemented using
MUSCLE.

For adjacent spatial scales, information needs to be exchanged between the
micro-model instances. This effectively equivalent to domain decomposition,
which MUSCLE 3 cannot (yet) do. In practice, compute-intensive models usually do
this internally anyway, so this is not a major limitation, although it would be
nice to have.

MUSCLE also supports time-domain adjacency, which is to say a process that
occurs after another one. In this way, it is semantically a superset of a
scientific workflow system.

Settings management
--------------------

Most models have some, or many, model parameters and other settings that
determine how the simulation proceeds. In a large coupled simulation, managing
these settings is tedious, as they need to be distributed over various
configuration files (or worse, source files), which then need to be put into the
correct location for models to pick them up. With MUSCLE, settings can be
defined in the yMMSL file that defines the simulation to run, and queried in
each submodel. This simplifies deployment and collecting metadata, and ensures
that settings that are used in different models are always the same (unless
explicitly set to different values).

For Uncertainty Quantification and Sensitivity Analysis, an ensemble is usually
run, with the values of some parameters varying amongst the ensemble members.
With MUSCLE 3, components can be introduced into the simulation that connect to
a set of instances and inject different values for some of the settings into
each one. This makes it possible to implement advanced semi-intrusive methods
for efficient UQ and SA of multiscale models.

Combining features
------------------

The most powerful feature of the MMSL, and MUSCLE 3, is the ability to
arbitrarily combine the above features in a single model. For instance, if your
macro-micro model exhibits both time scale and spatial scale separation, then
you can have a set of instances that each are reused many times. You can then
make a set of each of these (i.e. a set of macro-model instances, and a set of
sets of micro-model instances) and attach a UQ sampler to the macro-model in
order to run an Uncertainty Quantification. The overlaid parameters will then
be propagated automatically by MUSCLE 3 to the correct micro-model instances.

Other scenarios could include combining two macro-micro set-ups into a
macro-meso-micro model, or connecting a macro-model to an ensemble of stochastic
micro-models via a component that replicates the message from the macro-model to
each micro-model instance, and another component that calculates the ensemble
mean of the results and passes it back to the macro-model.

This feature is important because MUSCLE 3 is a tool for scientific research.
While the MUSCLE 3 developers cannot predict what you will invent in the future,
its flexibility greatly increases the probability that you will be able to use
MUSCLE 3 to implement it.


Ongoing development
-------------------

Language support
````````````````
Currently, MUSCLE 3 supports Python 3, C++ and Fortran. While that covers quite
some existing models in computational science, support for other languages is
still needed. This mainly entails porting libmuscle, as the manager is a
separate program, and writing tools to manipulate yMMSL files is probably best
done in Python.

Next on the list in terms of language support are C and possibly Java, then
perhaps other languages. If you feel strongly about support for a specific
language, please `make an issue on GitHub
<https://github.com/multiscale/muscle3/issues>`_.

Distributed execution
`````````````````````
Distributed execution is now officially supported, although not yet widely
tested. If you want to experiment on your laptop or your supercomputer,
please have a look at the `Distributed execution`_ section of the manual. We are
planning to use `QCG PilotJob <https://github.com/vecma-project/QCG-PilotJob>`_
to make running simulations inside of cluster allocations easier.

In order to avoid writing job scripts and staging files manually, some kind of
launcher is needed. We are collaborating with the `VECMA project
<https://www.vecma.eu>`_ on getting MUSCLE 3 support into `FabSim3
<https://fabsim3.readthedocs.io/en/latest/>`_.

Dynamic instantiation
`````````````````````
Some simulations require varying amounts of compute resources over the course of
a run. For instance, the time taken by a set of micro-models may depend on the
state of the macro-model, which changes over the course of the simulation. In
these cases, it may be necessary to change the number of instances of the
micro-model during the simulation. MUSCLE 3 does not yet support this. Some kind
of pilot job framework, such as `QCG PilotJob
<https://github.com/vecma-project/QCG-PilotJob>`_ will have to be integrated
with MUSCLE 3 to make this work.

It should be noted that in most cases, the resources have been allocated to the
user by a scheduler, and whether or not they are used does not affect the cost
of the simulation in core hours. Energy use is of course affected, and in some
cases wall-clock time can be reduced by redistributing the available resources
over fewer needed instances. Also, when running on a cloud, it may be possible
to return resources when they are no longer needed, and avoid paying for them.

The simple case of repeated instantiation of a micro-model is taken care of by
the model reuse facility that MUSCLE 3 does offer; dynamic instantiation is not
needed for this.

Profiling
`````````
MUSCLE 3 contains a partial implementation of a simple profiler, which can
measure the amount of time it takes to send messages between the instances.
While measurements are taken and the information is sent to the manager, it is
not yet saved to disk for further processing and not yet supported in C++. This
should be a simple addition.

Tutorial with Python
====================

In this section, we'll look at a few examples of how to use MUSCLE 3 to create a
multiscale simulation in Python.

`The source code for these examples is here
<https://github.com/multiscale/muscle3/tree/master/docs/source/examples/python>`_.
You can also clone the repository or download the source (see the C++ part of
the Installation section) and go to ``docs/source/examples/python``, or
copy-paste the code from here.

The easiest way to get set up is to create a virtualenv and then install MUSCLE
3 and the additional requirements inside it. Running ``make`` in the
``examples/python`` directory will make a virtual environment in
``examples/python/build/venv``:

.. code-block:: bash

  examples/python$ make

  # or by hand:
  examples/python$ python3 -m venv build/venv              # create venv
  examples/python$ . build/venv/bin/activate               # activate it
  examples/python$ pip3 install -r requirements.txt        # install dependencies


If you get an error message saying amongst others ``error: invalid command
'bdist_wheel'``, try running ``pip3 install wheel`` and then ``pip3 install -r
requirements.txt`` again to fix it. Or if you have administrator rights,
``apt-get install python3-wheel`` will also work, and fix this issue for the
whole system.

You can the run the example described below by activating the virtual
environment, and then running the file ``reaction_diffusion.py``:

.. code-block:: bash

  examples/python$ . build/venv/bin/activate
  exmaples/python$ python3 reaction_diffusion.py


Our first example is a reaction-diffusion model on a 1D grid. It consists of a
reaction model coupled to a diffusion model in a macro-micro fashion, with the
diffusion model the macro-model and the reaction model the micro-model. In a
macro-micro model with timescale separation, the micro model does a full run (of
many timesteps) for every timestep of the macro-model. Thus, the macro-model
effectively calls the micro-model before each state update, sending it the
current state and using its output as an input for the state update operation.

Here's how to implement that with MUSCLE 3. (A detailed explanation follows
below the code.)

.. literalinclude:: examples/python/reaction_diffusion.py
  :caption: ``docs/source/examples/python/reaction_diffusion.py``
  :language: python

Let's take it step by step.


Importing headers
-----------------

.. code-block:: python

  from collections import OrderedDict
  import logging
  import os

  import numpy as np

  from libmuscle import Grid, Instance, Message
  from libmuscle.runner import run_simulation
  from ymmsl import (ComputeElement, Conduit, Configuration, Model, Operator,
                     Settings)


As usual, we begin by importing some required libraries. OrderedDict and logging
come from the Python standard library, and NumPy provides matrix math
functionality. From libmuscle, we import :class:`libmuscle.Instance`, which will
represent a model instance in the larger simulation, and
:class:`libmuscle.Message`, which represents a MUSCLE message as passed between
instances.  The :func:`libmuscle.run_simulation` function allows us to run a
complete simulation from a single Python file, as opposed to having to start the
different instances and the manager separately.

In order to describe our model, we use a number of definitions from
`ymmsl-python`. More about those below.

A simple submodel
-----------------

Next is the reaction model. It is defined as a single Python function that takes
no arguments and returns nothing. (Note that this code uses type annotations to
show the types of function arguments and return values. These are ignored by the
Python interpreter, and also by MUSCLE 3, so you don't have to use them if you
don't want to.)

The first step in a MUSCLE model is to create an :class:`libmuscle.Instance`
object:

.. code-block:: python

  def reaction() -> None:
      instance = Instance({
              Operator.F_INIT: ['initial_state'],       # list of float
              Operator.O_F: ['final_state']})           # list of float
      # ...

The constructor takes a single argument, a dictionary that maps operators to
lists of ports. Ports are used by submodels to exchange messages with the
outside world. Ports can be input ports or output ports, but not both at the
same time. In this model, we have a single input port named ``initial_state``
that will receive an initial state at the beginning of the model run, and a
single output port named ``final_state`` that sends the final state of the
reaction simulation to the rest of the simulation. Note that the type of data
that is sent is documented in a comment. This is obviously not required, but it
makes life a lot easier if someone else needs to use the code or you haven't
looked at it for a while, so it's highly recommended. Available types are
described below.

The operators will become clear shortly.

The reuse loop
--------------

.. code-block:: python

  while instance.reuse_instance():
      # F_INIT
      # State update loop
      # O_F

Now that we have an :class:`libmuscle.Instance` object, we can start the *reuse
loop*. In multiscale simulations, submodels often have to run multiple times,
for instance because they're used as a micro-model or because they are part of
an ensemble that cannot be completely parallelised. In order to accomplish
this, the entire model is wrapped into a loop. Exactly when this loop ends
depends on the behaviour of the whole model, and is not easy to determine, but
fortunately MUSCLE will do that for us if we call the
:meth:`libmuscle.Instance.reuse_instance` method.

Initialisation: Settings and receiving messages
-----------------------------------------------

.. code-block:: python

  while instance.reuse_instance():
      # F_INIT
      t_max = instance.get_setting('t_max', 'float')
      dt = instance.get_setting('dt', 'float')
      k = instance.get_setting('k', 'float')

      msg = instance.receive('initial_state')
      U = msg.data.array.copy()

      t_cur = msg.timestamp


Next is the first part of the model, in which the state is initialised. Pretty
much every model has such an initialisation phase at the beginning. In MUSCLE's
abstract idea of what a submodel does, the Submodel Execution Loop, this part of
the code is referred to as the F_INIT operator. Here, we ask MUSCLE for the
values of some settings that we're going to need later, in particular the total
time to simulate, the time step to use, and the model parameter ``k``. The
second argument, which specifies the expected type, is optional. If it's given,
MUSCLE will check that the user specified a value of the correct type, and if
not raise an exception. Note that getting settings needs to happen *within* the
reuse loop; doing it before can lead to incorrect results.

After getting our settings, we receive the initial state on the
``initial_state`` port. Note that we have declared that port above, and declared
it to be associated with the F_INIT operator. During F_INIT, messages can only
be received, not sent, so that declaration makes ``initial_state`` a receiving
port.

The message that we receive contains several bits of information. Here, we are
interested in the ``data`` attribute, which we assume to be a grid of floats
containing our initial state, which we'll call ``U``. The ``msg.data`` attribute
holds an object of type :class:`libmuscle.Grid`, which holds a read-only NumPy
array and optionally a list of index names. Here, we take the array and make a
copy of it for ``U``, so that we can modify ``U`` in our upcoming state update.
Without calling ``.copy()``, ``U`` would end up pointing to the same read-only
array, and we'd get an error message if we tried to modify it.

Finally, we'll initialise our simulation time to the time at which that state is
valid, which is contained in the ``timestamp`` attribute. This is a
double-precision float containing the number of simulated (not wall-clock)
seconds since the whole simulation started.

The state update loop
---------------------

.. code-block:: python

  # F_INIT

  while t_cur + dt < msg.timestamp + t_max:
      # O_I

      # S
      U += k * U * dt
      t_cur += dt

  # O_F

Having initialised the model, it is now time for the state update loop. This
code should look familiar: we loop until we reach our maximum time (now
relative to when we started), and on each iteration update the state according
to the model equation. This update is called operator S in the Submodel
Execution Loop, and in this model, it is determined entirely by the current
state.  Since no information from outside is needed, we do not receive any
messages, and in our :class:`libmuscle.Instance` declaration above, we did not
declare any ports associated with ``Operator.S``.

The Submodel Execution Loop specifies another operator within the state update
loop, which is called O_I and comes before S. This operator provides for
observations of intermediate states. In other words, here is where you can send
a message to the outside world with (part of) the current state. In this case,
the O_I operator is empty; we're not sending anything. While this makes the
model a bit less reusable (it won't work as a macro-model like this), it is
perfectly legal. MUSCLE tries hard to let you break the rules unless doing so
would break the model, and in that case it tries to give a helpful error
message.

Sending the final result
------------------------

.. code-block:: python

  # O_F
  instance.send('final_state', Message(t_cur, None, Grid(U, ['x'])))


After the update loop is done, the model has arrived at its final state. We
finish up by sharing that final state with the outside world, by sending it on
the ``final_state`` port. The part of the code after the state update loop (but
still within the reuse loop) is known as the O_F operator in the Submodel
Execution Loop, so that is where we declared this port to live in our
:class:`libmuscle.Instance` declaration above.

To send a message, we specify the port on which to send (which must match the
declaration by name and operator), and a Message object containing the current
simulation time and the current state, converted to a plain Python list. The
optional second parameter is a second timestamp, which will be discussed below,
and is set to ``None`` here.

MUSCLE 3 uses `MessagePack <https://msgpack.org>`_ to encode messages between
models. MessagePack is a binary encoding format which can be thought of as a
binary version of JSON. That means that the message can be an integer, float,
bool, string, or a list or dictionary containing such, and MUSCLE 3 also
supports NumPy arrays and byte arrays. Like with JSON, these can be nested, so
you can send a dictionary containing lists of floats for example.

MessagePack is self-describing, so you can inspect the received message to find
out what you were sent. In most cases, all data you send on a particular
port will be of the same type, and you will expect data in a particular format
when receiving. We intentionally chose a self-describing format to keep you from
having to change a definition and recompile bindings every time you change what
you send. That does have a downside in that it makes it more difficult to see
what should be sent to someone elses submodel. So please make sure that you
document, for each port, which data type you're expecting or sending! Your
future colleagues (and possibly your future self) will thank you.

MessagePack is an extensible format, and since sending grids is very common in
these kinds of models MUSCLE 3 supports sending NumPy arrays directly. Here, we
wrap our array U into a :class:`libmuscle.Grid` object, so that we can add the
name of the dimensions. In this case there's only one, and ``x`` is not very
descriptive, so we could have also passed ``U`` directly, in which case MUSCLE
would have sent a :class:`libmuscle.Grid` without index names automatically.

Note that grids (and NumPy arrays) are much more efficient than lists and
dictionaries. If you have a lot of data to send, then you should use those as
much as possible. For example, a set of agents is best sent as a dictionary of
1D NumPy arrays, with one array/grid for each attribute. A list of dictionaries
will be quite a bit slower and use a lot more memory, but it should work up to
a million or so objects.

Finally, if you want to use your own encoding, you can just send a ``bytes``
object, which will be transmitted as-is, with minimal overhead.

This concludes our reaction model. As you can see, submodels that are used with
MUSCLE 3 very much look like they would without. With one additional variable,
one extra loop, and a few send and receive statements you're ready to connect
your model into a larger simulation.


Message timestamps
------------------

In order to make a coupled simulation, we need at least two models. The second
model is the diffusion model. Its overall structure is the same as for the
reaction model, but it has a few additional features. The first of these is an
O_I operator.

.. code-block:: python

  # O_I
  t_next = t_cur + dt
  if t_next + dt > t_max:
      t_next = None
  cur_state_msg = Message(t_cur, t_next, Grid(U, ['x']))
  instance.send('state_out', cur_state_msg)


Since the diffusion model is the macro-submodel in this model, it needs to send
its state to the outside world on every timestep. This is done in the O_I
operator. The message simply contains the state, converted to a
:class:`libmuscle.Grid`, and it is sent on the ``state_out`` port, which was
declared for the O_I operator when we made the :class:`libmuscle.Instance` for
this model. The message is sent with the current simulation time, and a second
timestamp that gives the simulation time for the next message that will be sent
on this port. Since our time steps are fixed, this is easy to calculate. We do
need to take care to send ``None`` if this is the final message on this port
however, since there won't be another message in that case.

This deserves a bit more explanation. First, MUSCLE 3 does not use the
timestamps that are attached to the messages for anything, and in this
particular case, always sending ``None`` for the second timestamp will work
fine.  These timestamps are necessary if two submodels with timescale overlap
need to be connected together. In this case, both models will run concurrently,
and they will exchange messages between their O_I and S operators.

If their time scales (step sizes) do not match exactly, then the number of
messages sent by one submodel does not match the number of messages that the
other submodel is trying to receive, and one will end up waiting for the other,
bringing the whole simulation to a halt. To avoid this, an intermediate
component needs to be inserted, which interpolates the messages. In order to
send and receive in the correct order, this component needs to know when the
next message will be sent, and that is when the second timestamp is required.

So, to make your submodel more generically usable, it's good to set the second
timestamp. But perhaps you're trying to connect an existing codebase that uses
varying timestep sizes, and it's not easy to get it to tell you how big the
next timestep will be. In that case, if you're not doing time scale overlap,
just put ``None`` there and move on to the next problem, it'll work just fine.

Receiving messages with a default
---------------------------------

.. code-block:: python

  # S
  msg = instance.receive('state_in', default=cur_state_msg)
  if msg.timestamp > t_cur + dt:
      logger.warning('Received a message from the future!')
  np.copyto(U, msg.data.array)


The diffusion model being the macro-model, it will need to receive input for its
state update from the micro-model, which it does by calling
:meth:`libmuscle.Instance.receive`. This receive is a bit special in that we
are passing a default message. The default message is returned if this port is
not connected. We are cleverly passing the message containing our current
state, so that if this port is not connected, the model continues from its
current state. Since MUSCLE 3 will simply ignore a send command on a
disconnected port, this makes it possible to run the diffusion model without a
micro-model attached.

Of course, a sensible default value will not be available in all cases, but if
there is one, using it like this is a good idea.

Next, we check that we didn't receive a message from the future. This can happen
if the micro-model runs for longer than our macro-model timestep. In this case,
the number of steps is fixed, so that this warning will never be emitted.
However, if the micro-model runs until it detects convergence, then it can
happen that it runs for longer than the timestep of the macro-model, and that
would indicate that there is no timescale separation anymore. In that case, the
result could be wrong, and a warning is appropriate.

Then, we copy the received data into our state array. The received
:class:`libmuscle.Grid` object contains a read-only array, and just writing ``U
= msg.data.array`` would discard the state we have in ``U``, and instead make
the variable ``U`` refer to the received read-only array. We would then get an
error message when trying to update the state. The ``np.copyto`` function
instead copies the (read-only) contents of ``msg.data.array`` into the existing
(writable) array referred to by ``U``. This way, we can do our state update, and
it's also a bit more efficient than reallocating memory all the time.

The S operator here calls the ``laplacian()`` function (not shown). There is no
requirement for a submodel to be a single function, you can split it up, call
library functions, and so on. There has to be a top-level function however if
you want to run more than one submodel in a single Python program. Also, you
cannot share any data with other components, other than by sending and receiving
messages. In particular, you can't use global variables to communicate between
models. This is intentional, because it doesn't work if you're running as
separate programs on different computers.


Connecting it all together
--------------------------

With both models defined, we now need to instruct MUSCLE 3 on how to connect
them together. We do this by creating an object of type
``ymmsl.Configuration``, which contains all the information needed to run the
simulation. It often helps to draw a diagram first:

.. figure:: reaction_diffusion.png
  :scale: 40 %
  :align: center

  gMMSL diagram of the reaction-diffusion model.


This is a gMMSL diagram of the reaction-diffusion model. It shows that there are
two compute elements named ``macro`` and ``micro``. A conduit connects port
``state_out`` on ``macro`` to ``state_in`` on ``micro``. The symbols at the ends
of the conduit show the operators that the ports belong to, O_I for
``macro.state_out`` and F_INIT for ``micro.state_in``. Another conduit connects
port ``micro.final_state`` (O_F) to ``macro.state_in`` (S).

Note that there's a mnemonic here: Operators O_I and S, which are within the
state update loop, have a circular symbol, while F_INIT and O_F use a diamond
shape. Also, filled symbols designate ports on which messages are sent, while
open symbols designate receiving ports. We can therefore see that for each state
update, ``macro`` will send on ``state_out``, after which ``micro`` will do a
full run and send its final result as input for the next state update of
``macro``.

Since diagrams aren't valid Python, we need an alternative way of describing
this model in our code. For this, we use the ``ymmsl`` library to create a set
of objects that form a description.

.. code-block:: python

    elements = [
            ComputeElement('macro', 'diffusion'),
            ComputeElement('micro', 'reaction')]


First, we describe the two compute elements in this model. Compute elements can
be submodels, or helper components that convert data, control the simulation, or
otherwise implement required non-model functionality. In this simple example, we
only have two submodels: one named ``macro`` and one named ``micro``. Macro is
implemented by an implementation named ``diffusion``, while micro is implemented
by an implementation named ``reaction``.

The name of a compute element is used by MUSCLE as an address for communication
between the models. The implementation name is intended for use by a launcher,
which would start the corresponding program to create an instance of a compute
element. It is these instances that form the actual running simulation.

.. code-block:: python

    conduits = [
            Conduit('macro.state_out', 'micro.initial_state'),
            Conduit('micro.final_state', 'macro.state_in')]

    model = Model('reaction_diffusion', elements, conduits)


Next, we need to connect the compute elements together. This is done by
conduits, which have a sender and a receiver. Here, we connect sending port
``state_out`` on compute element ``macro`` to receiving port ``initial_state``
on compute element ``micro``. Note that these ports are actually defined in the
implementations, and not in this configuration file, and they are referred to
here.

The compute elements and the conduits together form a ``Model``, which has a
name and those two sets of components.

.. code-block:: python

    settings = Settings(OrderedDict([
                ('micro.t_max', 2.469136e-6),
                ('micro.dt', 2.469136e-8),
                ('macro.t_max', 1.234568e-4),
                ('macro.dt', 2.469136e-6),
                ('x_max', 1.01),
                ('dx', 0.01),
                ('k', -4.05e4),     # reaction parameter
                ('d', 4.05e-2)      # diffusion parameter
                ]))

    configuration = Configuration(model, settings)


Finally, we define the settings for our simulation. We are using an
``OrderedDict`` instead of a normal (unordered) Python dictionary because having
these in a logical order is really really useful if you're trying to edit them
and quickly find things. For complex models, you'll have many settings, and
having them ordered allows you to group them logically. Of course, the order in
this source code will not change, but if you save this to YAML, you want to
preserve the order, and for that OrderedDict is required. (Unofficially in
Python 3.6, and officially in Python 3.7, ``dict`` is now ordered, but for
compatibility with older versions, we'll use an OrderedDict here.)

Note that there are two duplicated names between the two models: ``t_max`` and
``dt``. With MUSCLE 3, you can create a global setting with that name to set it
to the same value everywhere (handy for settings that are used by multiple
compute elements), or you can prepend the name of the compute element to set the
value for a specific one. Specific settings go before generic ones, so if you
specify both ``micro.t_max`` and ``t_max``, then the ``micro`` compute element
will use ``micro.t_max`` and all others will use ``t_max``.

The model and the settings are combined into a Configuration, which contains
everything needed to run the simulation. A Configuration object is also the
Python-side equivalent to a yMMSL YAML file.


Launching the simulation
------------------------

Launching a simulation is strictly speaking outside of the scope of MUSCLE 3,
which primarily does coordination (helping instances find each other) and
communication (sending messages between them). However, for testing,
experimentation, learning, and small scale production use, having some means of
running a whole simulation from a single Python file is very nice. So MUSCLE 3
has a small facility for this in the form of the
:func:`libmuscle.run_simulation` function:

.. code-block:: python

  implementations = {'diffusion': diffusion, 'reaction': reaction}
  run_simulation(configuration, implementations)


Here, we make a dictionary that maps implementation names to the corresponding
Python functions. We then pass the configuration and the implementations to
:func:`libmuscle.run_simulation`, which will start the muscle_manager and the
implementations, and then wait for the simulation to finish.

Note that this will actually fork off a separate process for the manager and for
each implementation, so that this is a real distributed simulation, albeit on a
single machine. That makes this mode of operation a more realistic test case
for a real distributed run on an HPC machine, smoothing the transition to larger
compute resources.


Log output
----------

If you run the script, e.g. using

.. code-block:: bash

  (venv) python$ python3 reaction_diffusion.py

it will pop up a plot showing the state of the simulated system over time. If
you are on a machine without graphics support, then you will get an error
message if you run the above, saying something like ``couldn't connect to
display``. In that case, try the below command to disable graphical output:

.. code-block:: bash

  (venv) python$ DONTPLOT=1 python3 reaction_diffusion.py

You will also find three log files in this directory: ``muscle3_manager.log``,
``muscle3.macro.log`` and ``muscle3.micro.log``. These contain log output for
the manager and the submodel instances respectively. You can log messages in the
usual Python way in your models, and MUSCLE 3 will automatically take care of
writing them to the log file. Any messages at level ``WARNING`` or higher will
be sent to the manager log as well. This helps give an overview of what went
wrong in a single place in case of errors.

Note that by default, Python (and MUSCLE 3) only logs messages at level
``WARNING`` or higher. So if you add a statement like

.. code-block:: python

  logging.info('Some useful information')

then you'll find that it will not show up in either of the log files, because
``INFO`` is a lower level than ``WARNING``. To get more or less local log output
(e.g. in ``muscle3.macro.log``), you can use one of these commands:

.. code-block:: python

  # least output
  logging.getLogger().setLevel(logging.CRITICAL)
  logging.getLogger().setLevel(logging.ERROR)
  logging.getLogger().setLevel(logging.WARNING)
  logging.getLogger().setLevel(logging.INFO)
  logging.getLogger().setLevel(logging.DEBUG)
  # most output

The minimum log level for sending a message to the central log file
(``muscle_manager.log``) is set by the ``muscle_remote_log_level`` setting. For
example, if you change the example to read

.. code-block:: python

    settings = Settings(OrderedDict([
                ('muscle_remote_log_level', 'DEBUG'),
                ('micro.t_max', 2.469136e-6),
                ('micro.dt', 2.469136e-8),
                ('macro.t_max', 1.234568e-4),
                ('macro.dt', 2.469136e-6),
                ('x_max', 1.01),
                ('dx', 0.01),
                ('k', -4.05e4),     # reaction parameter
                ('d', 4.05e-2)      # diffusion parameter
                ]))

then all log messages from all submodels will be sent to the manager log. This
does slow down the simulation if you have many log statements, so it's good to
use this to debug and to learn, but do increase the level again for production
runs.
.. muscle3 documentation master file, created by
   sphinx-quickstart on Wed Apr 19 10:25:17 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MUSCLE 3's documentation!
====================================================

MUSCLE 3 is the third incarnation of the Multiscale Coupling Library and
Environment.

With MUSCLE 3, you can connect multiple simulation models together into a
multiscale simulation. Simulation models can be as simple as a single Python
file, or as complex as a combination of multiple separate simulation codes
written in C++ or Fortran, and running on an HPC machine.

If you use MUSCLE 3 for scientific work, please `cite
the version of the MUSCLE 3 software you used
<https://www.research-software.nl/software/muscle3>`_ and the following paper:

Veen L.E., Hoekstra A.G. (2020) Easing Multiscale Model Design and Coupling
with MUSCLE 3. In: Krzhizhanovskaya V. et al. (eds) Computational Science –
ICCS 2020. ICCS 2020.  Lecture Notes in Computer Science, vol 12142. Springer,
Cham.  `<https://doi.org/10.1007/978-3-030-50433-5_33>`_


.. toctree::
   :maxdepth: 4
   :caption: Contents:

   introduction
   current_status
   installing
   tutorial
   distributed_execution
   cplusplus
   fortran
   mpi
   uncertainty_quantification
   python_api
   cpp_api
   fortran_api

   contributing
   devtools
   releasing


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
API Documentation for C++
=========================

This page provides full documentation for the C++ API of MUSCLE 3.

Note that in a few places, classes are referred to as
``libmuscle::impl::<class>``. This is a bug in the documentation rendering
process, the class is actually available as ``libmuscle::<class>`` and
should be used as such.

Namespace libmuscle
-------------------

.. doxygenclass:: libmuscle::impl::Data
.. doxygenclass:: libmuscle::impl::DataConstRef
.. doxygenclass:: libmuscle::impl::Instance
.. doxygenclass:: libmuscle::impl::Message
.. doxygentypedef:: libmuscle::impl::PortsDescription


Namespace ymmsl
---------------

.. doxygenfunction:: ymmsl::impl::allows_sending
.. doxygenfunction:: ymmsl::impl::allows_receiving

.. doxygenclass:: ymmsl::impl::Conduit
.. doxygenclass:: ymmsl::impl::Identifier
.. doxygenfunction:: ymmsl::impl::operator<<(std::ostream&, Identifier const&)
.. doxygenenum:: ymmsl::impl::Operator
.. doxygenstruct:: ymmsl::impl::Port
.. doxygenclass:: ymmsl::impl::Reference
.. doxygenfunction:: ymmsl::impl::operator<<(std::ostream&, Reference const&)
.. doxygenclass:: ymmsl::impl::ReferencePart
.. doxygenclass:: ymmsl::impl::Settings
.. doxygenfunction:: ymmsl::impl::operator<<(std::ostream&, ymmsl::impl::Settings const&)
.. doxygenclass:: ymmsl::impl::SettingValue
.. doxygenfunction:: ymmsl::impl::operator<<(std::ostream&, ymmsl::impl::SettingValue const&)

Introduction
============

Welcome to MUSCLE 3!

MUSCLE 3 is the third incarnation of the Multiscale Coupling Library and
Environment, and the successor to MUSCLE 2. Its purpose is to make creating
coupled multiscale simulations easy, and to then enable efficient Uncertainty
Quantification of such models using advanced semi-intrusive algorithms.

MUSCLE 3 uses the Multiscale Modelling and Simulation Language (MMSL) to
describe the structure of a multiscale model. MMSL can be expressed in the form
of a diagram (gMMSL; not yet implemented) or as a YAML file (yMMSL; this is
convenient both for people and for software). The MMSL lets one describe which
*compute elements* (submodels, scale bridges, data converters, UQ components,
etc.) a multiscale model consist of, how many instances of each we need, and how
they are wired together.

MUSCLE 3 is intended to scale from your laptop to the exascale. At the low end,
it supports a non-distributed but parallel mode in which an entire multiscale
simulation, including all compute element implementations and the MMSL
configuration, is in a single (short) Python file. Next is a distributed mode
where the manager and compute element instances are started on multiple nodes
in a cluster, and communicate directly with one another. Beyond that, additional
components and optimisations are envisioned that would allow scaling to huge
machines or combinations of multiple machines. Our goal is to make the
transitions between these modes as smooth as possible, so that extra compute
power can be added gradually, as needed.

MUSCLE 3 consist of three components: **libmuscle**, **ymmsl-python**, and
**muscle_manager**. **Libmuscle** is a library that is used to connect compute
element implementations (e.g. new or existing models) to a coupled simulation.
It offers functions for sending and receiving messages to and from the outside
world, determining whether to restart the model at the end of its run, and
various utilities like centralised settings (for convenience and UQ),
centralised logging (to make debugging easier) and basic profiling.

**ymmsl-python** is a Python library that contains class definitions to
represent an MMSL model description. It supports loading a yMMSL file and
converting it to objects, and saving objects to a yMMSL file. Input is validated
and a clear error message is provided in case of errors, and output is formatted
for good readability, courtesy of the `YAtiML <https://yatiml.readthedocs.io>`_
library.  yMMSL files can also contain settings, so that all the information
needed for a simulation run can be described in a single file.

**muscle_manager** is the central run-time component of MUSCLE 3. It is started
together with the compute element intances, and provides a central coordination
point that the instances use to find each other. The manager also collects log
messages from the individual instances to aid in debugging, and some profiling
information to aid in scheduling and performance optimisation.

API Documentation for Python
============================

.. automodule:: libmuscle
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: libmuscle.runner
   :members: run_simulation
MUSCLE and C++
==============

This section shows how to use MUSCLE 3 from C++, based on the same
reaction-diffusion model as given in the Python tutorial.
`The source code for the examples in this section is here
<https://github.com/multiscale/muscle3/tree/master/docs/source/examples/cpp>`_.
You can also go to ``docs/source/examples/cpp`` in the source directory (that's
easiest, and has some handy scripts and a Makefile), or copy-paste the code from
here.

Building and running the examples
---------------------------------

If you've just built and installed the C++ version of libmuscle, then you're all
set to build the examples. To do that, go into the ``docs/source/examples``
subdirectory, and run Make, passing the installed location of MUSCLE 3:

.. code-block:: bash

  ~/muscle3_source$ cd docs/source/examples
  ~/muscle3_source/docs/source/examples$ MUSCLE3_HOME=~/muscle3 make cpp


We also need the Python version of MUSCLE 3 installed, because the MUSCLE
Manager comes with that, and we need it to run the simulation. To set that up,
do:

.. code-block:: bash

  ~/muscle3_source/docs/source/examples$ make python


You can then run the examples using the provided scripts in
``docs/source/examples``:

.. code-block:: bash

  ~/muscle3_source/docs/source/examples$ MUSCLE3_HOME=~/muscle3 ./reaction_diffusion_cpp.sh


Log output
----------

When you run a MUSCLE 3 simulation, the manager will produce a
``muscle3_manager.log`` file with log messages collected from the submodels. You
can set the log level for C++ instances by adding a ``muscle_remote_log_level``
setting to the yMMSL file, with a value of ``DEBUG``, ``INFO``, ``WARNING``,
``ERROR`` or ``CRITICAL``. The default is ``WARNING``, which is pretty quiet,
``INFO`` will give more output, but also slow things down a bit because those
messages have to be sent to the manager.

The running script will redirect standard out and standard error for each
instance to a file, so you can find any output produced by the submodels there.


Reaction-diffusion in C++
-------------------------

This page assumes that you've at least read the Python tutorial; we'll do the
same thing again but now in C++.

The C++ version of libmuscle doesn't support running multiple compute elements
in a single program, so we're going to do the distributed version of the
reaction-diffusion example in C++. Here's the basic reaction model, without
MUSCLE support, in C++:

Reaction model
``````````````

.. code-block:: cpp
    :caption: A simple reaction model on a 1D grid, in C++.

    void reaction() {
        // F_INIT
        double t_max = 2.469136e-07;
        double dt = 2.469136e-08;
        double x_max = 1.0;
        double dx = 0.01;
        double k = -40500.0;

        std::vector<double> U(lrint(x_max / dx));
        U[25] = 2.0;
        U[50] = 2.0;
        U[75] = 2.0;

        double t_cur = 0.0;
        while (t_cur + dt < t_max) {
            // O_I

            // S
            for (double & u : U)
                u += k * u * dt;
            t_cur += dt;
        }

        // O_F
    }


This is one of the simplest possible computational models: we set some
parameters, create a state variable ``U``, initialise the state, and then update
the state in a loop until we reach the final simulation time.

The MUSCLE version in C++ looks quite similar to the Python version:

.. literalinclude:: examples/cpp/reaction.cpp
  :caption: ``docs/source/examples/cpp/reaction.cpp``
  :language: cpp


We'll go through it top to bottom, one piece at a time.

Headers
```````

.. code-block:: cpp

    #include <cstdlib>

    #include <libmuscle/libmuscle.hpp>
    #include <ymmsl/ymmsl.hpp>


    using libmuscle::Data;
    using libmuscle::DataConstRef;
    using libmuscle::Instance;
    using libmuscle::Message;
    using ymmsl::Operator;


Here we include some headers, and import the classes we'll need into our root
namespace. ``cstdlib`` is there for the ``EXIT_SUCCESS`` constant at the bottom,
and not needed for MUSCLE. The C++ API mirrors the Python API where applicable,
so it has ``libmuscle`` and ``ymmsl`` like in Python. However, in C++ there is
no separate yMMSL library. Instead the classes needed are included with
libmuscle. This means that there is no support for manipulating yMMSL files from
C++.

There are two classes here that don't have a Python equivalent, ``Data`` and
``DataConstRef``. More on that when we use them below.

Creating an Instance
````````````````````

.. code-block:: cpp

    void reaction(int argc, char * argv[]) {
        Instance instance(argc, argv, {
                {Operator::F_INIT, {"initial_state"}},  // list of double
                {Operator::O_F, {"final_state"}}});     // list of double


The ``reaction`` function has changed a bit: it now takes the command line
parameters as arguments. Like the Python version, the C++ libmuscle takes some
information from the command line, but in C++ this needs to be passed
explicitly to the ``Instance`` object. Otherwise, constructing the ``Instance``
is the same: it is passed a dictionary mapping operators to a list of port
descriptions (see ``PortsDescription`` in the API documentation). Like in the
Python version, we'll be passing a list of doubles between the models.

Reuse loop and settings
```````````````````````

.. code-block:: cpp

    while (instance.reuse_instance()) {
        // F_INIT
        double t_max = instance.get_setting_as<double>("t_max");
        double dt = instance.get_setting_as<double>("dt");
        double k = instance.get_setting_as<double>("k");


As in Python, we have a reuse loop, and except for the syntax differences
between the languages, the code is exactly the same. Getting settings works as
shown; with C++ being a statically typed language, we need to specify the type
we're expecting. As in Python, if the actual type is different, an exception
will be thrown.

Supported types for settings are ``std::string``, ``bool``, ``int64_t``,
``double``, ``std::vector<double>`` and ``std::vector<std::vector<double>>``.
There is also a ``get_setting()`` member function, which returns a
``SettingValue`` object. A ``SettingValue`` can contain a value of any of the
above types, and it can be queried to see what is in it.

Receiving messages and DataConstRef
```````````````````````````````````

Instead of initialising our state from a constant, we're going to receive a
message on our ``F_INIT`` port with the initial state:

.. code-block:: cpp

    auto msg = instance.receive("initial_state");
    auto data_ptr = msg.data().elements<double>();
    std::vector<double> U(data_ptr, data_ptr + msg.data().size());


Calling the ``receive`` method on an instance yields an object of type
``libmuscle::Message``. In Python, this is a very simple class with several
public members. In C++, these members are wrapped in accessors, so while the API
is conceptually the same, the syntax is slightly different. Here, we're
interested in the data that was sent to us, so we use ``msg.data()`` to access
it.

This returns an object of type ``libmuscle::DataConstRef``. This type is new in
the C++ version of the library, and it exists because C++ is a statically typed
language, while in a MUSCLE simulation the type of data that is sent in a
message can vary from message to message. So we need some kind of class that can
contain data of many different types, and that's what ``Data`` and
``DataConstRef`` are for. Here are some key properties of ``DataConstRef``:

- ``DataConstRef`` is like a reference in that it points to an actual data
  object of some type. If you copy it into a second ``DataConstRef``, both
  ``DataConstRef`` variables will point to the same object. Like in Python or in
  Java.

- ``DataConstRef`` is constant, you cannot change the data through a
  ``DataConstRef``.

- Memory management for ``DataConstRef`` is automatic, there's no need to
  ``delete`` or ``free`` anything. Objects will be deleted from memory when all
  ``Data`` and ``DataConstRef`` objects pointing to them are gone. Like in
  Python or in Java.

- ''DataConstRef'' variables can be inspected to see what type of data they
  contain, and the data can be extracted into a C++ variable of the correct
  type.

In the code here, we don't bother with a check. Instead, we blindly assume that
we've been sent a 1D grid of doubles. If there is no grid, an exception will
be thrown and our program will halt. That's okay, because it means that there's
something wrong somewhere that we need to fix. MUSCLE is designed to let you get
away with being a bit sloppy as long as things actually go right, but it will
check for problems and let you know if something goes wrong.

A grid in MUSCLE 3 is an n-dimensional array of numbers or booleans. It has a
shape (size in each dimension), a size (total number of elements), a storage
order which says in which order the elements are arranged in memory, and
optionally names for the indexes.

Here, we expect a 1D grid of doubles, which is just a vector of numbers. Storage
order is irrelevant then. We'll use the standard C++ ``std::vector`` class to
contain our state. It has a constructor that takes a pointer to an array of
objects of the appropriate type and the number of them, and copies those objects
into itself. We use the :cpp:func:`DataConstRef::elements` to get a pointer to
the elements, and :cpp:func:`DataConstRef::size` to get the number of elements,
and that's all we need to create our state vector ``U``.

Note that MUSCLE 3 will in this case check that we have the right type of
elements (doubles), but it cannot know and therefore will not check that we're
expecting a 1D grid. We could check that by hand by checking the length of
``msg.data().shape()`` to see if it's 1. See the diffusion model below for an
example.


.. code-block:: cpp

        double t_cur = msg.timestamp();
        while (t_cur + dt < msg.timestamp() + t_max) {
            // O_I

            // S
            for (double & u : U)
                u += k * u * dt;
            t_cur += dt;
        }

The main loop of the model is almost unchanged, we just get the starting time
from the received message's timestamp so as to synchronise with the overall
simulation time correctly, and the stopping time is adjusted accordingly as
well.

Sending messages and Data
`````````````````````````

.. code-block:: cpp

        // O_F
        auto result = Data::grid(U.data(), {U.size()}, {"x"});
        instance.send("final_state", Message(t_cur, result));


Having computed our final state, we will send it to the outside world on the
``final_state`` port. In this case, we need to send a 1D grid of doubles, which
we first need to wrap up into a ``Data`` object. A ``Data`` object works just
like a ``DataConstRef``, except that it isn't constant, and can thus be
modified. (It is in fact a reference, like ``DataConstRef``, despite the name,
and it has automatic memory management as well.)

Here, we create a ``Data`` containing a grid. We pass it a pointer to the
numbers in ``U``, and the shape of our grid, which is a 1-element array because
the grid is 1-dimensional. The third argument contains the names of the indexes,
which helps to avoid confusion over which index is x and which is y (or
row/column, or latitude/longitude, or... you get the idea). You are not
required to add these (and MUSCLE 3 doesn't use them), but you or someone else
using your code will be very grateful you did at some point in the future.

With our data item constructed, we can send a ``Message`` containing the current
timestamp and the data to the ``final_state`` port. Note that there are
different ways of creating a ``Message``, depending on whether you set the next
timestamp, or are using an explicit settings overlay. See the API documentation
for details.

``Data`` is quite versatile, and makes it easier to send data of various
types between submodels. Here are some other examples of creating ``Data``
objects containing different kinds of data:

.. code-block:: cpp

    // create a Data containing a nil value (None in Python)
    Data d1;

    // strings, booleans
    Data d2("String data");
    Data d3(true);

    // various kinds of numbers
    Data d4(0);         // int
    Data d5(10u);       // unsigned int
    Data d6(100l);      // long int
    Data d7(3.141592f); // float
    Data d8(1.4142);    // double

    // constant list
    auto d9 = Data::list("String data", true, 0, 10u, 1.4142);

    // dictionary
    auto d10 = Data::dict(
            "x", x,
            "y", y,
            "note", "Keys must be strings",
            "why_not", Data::list(10, 3.141592, "that's fine"),
            "existing_object", d4
            );

    // grid
    // +-----------------+
    // | 1.0 | 2.0 | 3.0 |
    // +-----------------+
    // | 4.0 | 5.0 | 6.0 |
    // +-----------------+
    std::vector<double> array = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    auto d11 = Data::grid(array.data(), {2, 3}, {"row", "col"});

    // byte array
    std::vector<char> bytes;
    auto d12 = Data::byte_array(bytes.data(), bytes.size());

    // simple types are created automatically when making a Message
    instance.send("output_port", Message(t_cur, 1.0));

    std::string text("A text message");
    instance.send("output_port", Message(t_cur, text));

    // but you have to name dicts, lists, grids and byte arrays
    instance.send("output_port", Message(t_cur, Data::list("a", 10)));


As you can see, sending complex data types with MUSCLE is almost as easy in C++
as it is in Python. There is however a bit of a cost to this ease-of-use in the
form of extra memory usage. If you're sending large amounts of data, then it is
best to use as many grids as you can, as those are much more efficient than
dictionaries and lists. In particular, a set of objects (agents for example) is
much more efficiently sent as a dictionary-of-grids (with one 1D-grid for each
attribute) than as a list of dictionaries. (The latter will work up to a million
objects or so, but it will be slower.)

.. code-block:: cpp

    int main(int argc, char * argv[]) {
        reaction(argc, argv);
        return EXIT_SUCCESS;
    }


We finish the example with a simple ``main()`` function that runs the model,
passing the command line arguments.

Diffusion model
```````````````

If you've studied the above carefully, and have seen the Python version of the
diffusion model, then you should now be able to understand the C++ diffusion
model below:

.. literalinclude:: examples/cpp/diffusion.cpp
  :caption: ``docs/source/examples/cpp/diffusion.cpp``
  :language: cpp


In the examples directory, you will find a handy script called
``reaction_diffusion_cpp.sh``, which runs the C++ reaction-diffusion model
locally. This script launches the MUSCLE Manager and an instance of each
submodel, then waits for the simulation to complete. At the top of this page
there are instructions on how to run the examples.

Also in this directory are C++ versions of a Monte Carlo sampler and a
round-robin load balancer, plus a script to launch this extended example. See
the `Uncertainty Quantification` section to learn what those are for.

MUSCLE and Fortran
==================

This section shows how to use MUSCLE 3 from Fortran, based on the same
reaction-diffusion model given in the Python tutorial.
`The source code for the examples in this section is here
<https://github.com/multiscale/muscle3/tree/master/docs/source/examples/fortran>`_.
You can also go to ``docs/source/examples/fortran`` in the source directory
(that's easiest, and has some handy scripts and a Makefile), or copy-paste the
code from here.

Building and running the examples
---------------------------------

If you've just built and installed the C++ version of libmuscle (which includes
the Fortran bindings), then you're all set to build the examples. To do that, go
into the ``docs/source/examples`` subdirectory, and run Make:

.. code-block:: bash

  ~/muscle3_source$ cd docs/source/examples
  ~/muscle3_source/docs/source/examples$ MUSCLE3_HOME=~/muscle3 make fortran


We also need the Python version of MUSCLE 3 installed, because the MUSCLE
Manager comes with that, and we need it to run the simulation. To set that up,
do:

.. code-block:: bash

  ~/muscle3_source/docs/source/examples$ make python


You can then run the examples using the provided scripts in
``docs/source/examples``:

.. code-block:: bash

  ~/muscle3_source/docs/source/examples$ MUSCLE3_HOME=~/muscle3 ./reaction_diffusion_fortran.sh


Log output
----------

When you run a MUSCLE 3 simulation, the manager will produce a
``muscle3_manager.log`` file with log messages collected from the submodels. You
can set the log level for Fortran instances by adding a
``muscle_remote_log_level`` setting to the yMMSL file, with a value of
``DEBUG``, ``INFO``, ``WARNING``, ``ERROR`` or ``CRITICAL``. The default is
``WARNING``, which is pretty quiet, ``INFO`` will give more output, but also
slow things down a bit because those messages have to be sent to the manager.

The running script will redirect standard out and standard error for each
instance to a file, so you can find any output produced by the submodels there.


Reaction-diffusion in Fortran
-----------------------------

This page assumes that you've at least read the Python tutorial; we'll do the
same thing again but now in Fortran.

The Fortran version of libmuscle doesn't support running multiple compute
elements in a single program, so we're going to do the distributed version of
the reaction-diffusion example in Fortran. Here's the basic reaction model,
without MUSCLE support, in Fortran:

Reaction model
``````````````

.. code-block:: fortran
    :caption: A simple reaction model on a 1D grid, in Fortran.

    program reaction
        implicit none

        real (selected_real_kind(15)) :: t_cur, t_max, dt, x_max, dx, k
        integer :: i, U_size
        real (selected_real_kind(15)), dimension(:), allocatable :: U


        ! F_INIT
        t_max = 2.469136e-07
        dt = 2.469136e-08
        x_max = 1.0
        dx = 0.01
        k = -40500.0

        U_size = int(x_max / dx)
        allocate (U(U_size))
        U = 1e-20
        U(26) = 2.0
        U(51) = 2.0
        U(76) = 2.0

        t_cur = 0.0
        do while (t_cur + dt < t_max)
            ! O_I

            ! S
            U = k * U * dt
            t_cur = t_cur + dt
        end do

        ! O_F

        deallocate (U)

    end program reaction


This is one of the simplest possible computational models: we set some
parameters, create a state variable ``U``, initialise the state, and then update
the state in a loop until we reach the final simulation time.

The MUSCLE version in Fortran looks quite similar to the Python version:

.. literalinclude:: examples/fortran/reaction.f03
  :caption: ``docs/source/examples/fortran/reaction.f03``
  :language: fortran


We'll go through it top to bottom, one piece at a time.

Modules
```````

.. code-block:: fortran

    use ymmsl
    use libmuscle
    implicit none


Here we tell Fortran that we'll be using the ymmsl and libmuscle modules. These
mirror the corresponding Python packages. Like in C++, yMMSL support is limited
to what is needed to implement compute elements. Loading, manipulating and
saving yMMSL documents is better done in Python.


Variables
`````````

.. code-block:: fortran

    type(LIBMUSCLE_PortsDescription) :: ports
    type(LIBMUSCLE_Instance) :: instance

    type(LIBMUSCLE_Message) :: rmsg
    type(LIBMUSCLE_DataConstRef) :: rdata, item


    type(LIBMUSCLE_Message) :: smsg
    type(LIBMUSCLE_Data) :: sdata

    real (selected_real_kind(15)) :: t_cur, t_max, dt, k
    real (selected_real_kind(15)), dimension(:), allocatable :: U


Next, it's time to declare our variables. We have a few extra variables here
compared to the non-MUSCLE version. The :f:type:`LIBMUSCLE_PortsDescription` is
used to describe the ports we'll use to send and receive to MUSCLE 3. We need a
:f:type:`LIBMUSCLE_Instance` as the main interface to MUSCLE. We'll receive
messages into ``rmsg``, a :f:type:`LIBMUSCLE_Message` variable, and extract the
data from them into ``rdata`` and ``item``. Since received data cannot be
modified, this variable is of type :f:type:`LIBMUSCLE_DataConstRef`. For sending
messages we have ``smsg``, and we'll create a :f:type:`LIBMUSCLE_Data` object to
put our data into before sending. Note that the names are all prefixed with
``LIBMUSCLE``, to make sure that they don't collide with any other library you
may want to use.

Eagle-eyed readers will have noticed that ``dx`` and ``x_max`` are missing. That
is because we'll derive the size of the state vector of the model (``U``) from
the state we receive, rather than from the configuration.

Creating an Instance
````````````````````

.. code-block:: fortran

    ports = LIBMUSCLE_PortsDescription_create()
    call LIBMUSCLE_PortsDescription_add(ports, YMMSL_Operator_F_INIT, 'initial_state')
    call LIBMUSCLE_PortsDescription_add(ports, YMMSL_Operator_O_F, 'final_state')
    instance = LIBMUSCLE_Instance_create(ports)
    call LIBMUSCLE_PortsDescription_free(ports)


In order to talk to the rest of the MUSCLE simulation, we need a
:f:type:`LIBMUSCLE_Instance` object, and we need to pass a description of the
ports we'll use when we make that. So here, we first create a
:f:type:`LIBMUSCLE_PortsDescription`, and add some ports to the ``F_INIT`` and
``O_F`` operators. We can then create the Instance object, passing the ports.
Finally, we need to free the PortsDescription object. Since Fortran does not
have automatic memory management, you will have to be careful to always free any
object you create, or that is returned by a MUSCLE function, after you're done
using it. So we free the PortsDescription here, but we don't free the Instance,
since we need it in the following part of the code.


Reuse loop and settings
```````````````````````

.. code-block:: fortran

    do while (LIBMUSCLE_Instance_reuse_instance(instance))
        ! F_INIT
        t_max = LIBMUSCLE_Instance_get_setting_as_real8(instance, 't_max')
        dt = LIBMUSCLE_Instance_get_setting_as_real8(instance, 'dt')
        k = LIBMUSCLE_Instance_get_setting_as_real8(instance, 'k')


As in Python, we have a reuse loop, and except for the syntax differences
between the languages, the code is exactly the same. Getting settings works as
shown; with Fortran being a statically typed language, we need to specify the
type we're expecting by calling a different function. If the actual type is
different, an error will be printed and the program will be halted. It is
possible to detect and handle errors instead, see below for that.

Getting settings is done via the ``LIBMUSCLE_Instance_get_setting_as_<type>()``
functions. Supported types are ``character``, ``logical``, ``int8`` (a 64-bit
integer), ``real8`` (a 64-bit double precision real number), ``real8array`` (a
1D array of double precision reals) and ``real8array2`` (a 2D array of double
precision reals). See the :f:func:`API documentation
<LIBMUSCLE_Instance_get_setting_as_character>` for details.


Receiving messages and DataConstRef
```````````````````````````````````

Instead of initialising our state from a constant, we're going to receive a
message on our ``F_INIT`` port with the initial state:

.. code-block:: fortran

        rmsg = LIBMUSCLE_Instance_receive(instance, 'initial_state')
        rdata = LIBMUSCLE_Message_get_data(rmsg)
        allocate (U(LIBMUSCLE_DataConstRef_size(rdata)))
        call LIBMUSCLE_DataConstRef_elements(rdata, U)
        call LIBMUSCLE_DataConstRef_free(rdata)


Calling the :f:func:`LIBMUSCLE_Instance_receive` function with an instance and a
port name yields an object of type :f:type:`LIBMUSCLE_Message` containing the
received message. As always when using MUSCLE from Fortran, we have to remember
to free the returned Message object when we are done with it. That's not yet
the case though, because we still need to get the received data out. In Python,
Message is a very simple class with several public members. In Fortran, objects
are always accessed and manipulated through LIBMUSCLE functions, in this case
:f:func:`LIBMUSCLE_Message_get_data`, which we use to get a ``DataConstRef``
object. Again, since MUSCLE gave us an object, we have to remember to free it
when we're done.

First though, we'll use the data to initialise our state. We are going to assume
that we'll receive a 1D grid of numbers, just like in the equivalent examples in
the other supported languages. :f:func:`LIBMUSCLE_DataConstRef_size` will return
the size of the array that we received (or print an error and stop if we
received something else that doesn't have a size, see below for handling errors
differently). We use that to allocate the ``U`` array to the same size as the
received state, and then we copy the elements of the array into ``U`` using
:f:func:`LIBMUSCLE_DataConstRef_elements`. We can then free the ``rdata``
object, since we don't need it any longer. We'll hold on to the ``rmsg`` a bit
longer, since we'll need it later.

If all this freeing objects is starting to sound tedious, that's because it is,
and it's why more modern languages like Python and C++ do this for you.
Unfortunately, for Fortran, it has to be done manually.

Note that indices for the received array start at 1, as usual in Fortran. MUSCLE
3 follows the language in which you're using it and automatically translates, so
if this grid was sent from Python or C++, then received item 1 corresponds to
sent item 0.

If you have a DataConstRef object, then you can check which kind of value it
contains, e.g. using :f:func:`LIBMUSCLE_DataConstRef_is_a_grid_of_real8`. Here,
we don't bother with a check. Instead, we blindly assume that we've been sent a
1D grid of doubles. If that's not the case, an error will be printed and our
program will halt. That's okay, because it means that there's something wrong
somewhere that we need to fix. MUSCLE is designed to let you get away with being
a bit sloppy as long as things actually go right, but it will check for problems
and let you know if something goes wrong. If you want to make a submodel or
component that can handle different kinds of messages, then these inspection
functions will help you do so however.


.. code-block:: fortran

        t_cur = LIBMUSCLE_Message_timestamp(rmsg)
        t_max = LIBMUSCLE_Message_timestamp(rmsg) + t_max
        call LIBMUSCLE_Message_free(rmsg)

        do while (t_cur + dt < t_max)
            ! O_I

            ! S
            U = k * U * dt
            t_cur = t_cur + dt
        end do


The main loop of the model is almost unchanged, we just get the starting time
from the received message's timestamp so as to synchronise with the overall
simulation time correctly. The stopping time is adjusted accordingly as
well. Note that ``t_cur`` and ``t_max`` are numbers, not objects, so they do not
have to be freed. ``rmsg`` however does, and since we have all the information
we need from the received message, we free it here.


Sending messages and Data
`````````````````````````

.. code-block:: fortran

        ! O_F
        sdata = LIBMUSCLE_Data_create_grid(U, 'x')
        smsg = LIBMUSCLE_Message_create(t_cur, sdata)
        call LIBMUSCLE_Instance_send(instance, 'final_state', smsg)
        call LIBMUSCLE_Message_free(smsg)
        call LIBMUSCLE_Data_free(sdata)
        deallocate (U)
    end do

    call LIBMUSCLE_Instance_free(instance)


Having computed our final state, we will send it to the outside world on the
``final_state`` port. In this case, we need to send a vector of doubles, which
we first need to wrap up into a ``Data`` object. A ``Data`` object works just
like a ``DataConstRef``, except that it isn't constant, and can thus be
modified. We do this by creating a ``Data`` object containing a grid value, with
the array ``U`` and the index name ``x``.

Having put our data into the ``Data`` object is a grid, we can then put the grid
into a new :f:type:`LIBMUSCLE_Message` object, and call the
:f:func:`LIBMUSCLE_Instance_send` function to send it. Finally, we free the
Message and Data objects that we created, and deallocate the state as in the
original non-MUSCLE version.

That concludes the reuse loop. When we're done running the model, the reuse loop
will finish, and we can free our Instance object before we quit.

:f:type:`LIBMUSCLE_Data` is quite versatile, and makes it easier to send data of
various types between submodels. Here are some other examples of creating
:f:type:`LIBMUSCLE_Data` objects containing different kinds of data (note that
freeing objects is omitted here for brevity, of course you have to do that in
your model!):

.. code-block:: fortran

    type(LIBMUSCLE_Data) :: d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12
    integer (LIBMUSCLE_int8), dimension(2, 3) :: ar
    character(len=1), dimension(1024) :: bytes

    ! create a Data containing a nil value (None in Python)
    d1 = LIBMUSCLE_Data_create()

    ! character strings, logicals
    d2 = LIBMUSCLE_Data_create('String data')
    d3 = LIBMUSCLE_Data_create(.true.)

    ! various kinds of numbers
    d4 = LIBMUSCLE_Data_create(0)                   ! integer
    d5 = LIBMUSCLE_Data_create(100_LIBMUSCLE_int8)  ! 64-bit integer
    d6 = LIBMUSCLE_Data_create(3.141592_LIBMUSCLE_real4)    ! 32-bit real
    d7 = LIBMUSCLE_Data_create(1.4142_LIBMUSCLE_real8)      ! 64-bit real

    ! constant list
    d8 = LIBMUSCLE_Data_create_list()                   ! empty list
    d9 = LIBMUSCLE_Data_create_nils(4_LIBMUSCLE_size)   ! list of length 4
    call LIBMUSCLE_Data_set_item(d9, 1, 'String data')
    call LIBMUSCLE_Data_set_item(d9, 2, .true.)
    call LIBMUSCLE_Data_set_item(d9, 3, 0)
    call LIBMUSCLE_Data_set_item(d9, 4, 1.4142d0)

    ! dictionary
    d10 = LIBMUSCLE_Data_create_dict()
    call LIBMUSCLE_Data_set_item(d10, 'x', x)
    call LIBMUSCLE_Data_set_item(d10, 'y', y)
    call LIBMUSCLE_Data_set_item(d10, 'note', 'Keys must be strings')
    call LIBMUSCLE_Data_set_item(d10, 'nest all you want', d9)

    ! grid
    ar = reshape(spread((/1_LIBMUSCLE_int8/), 1, 6), (/2, 3/))
    d12 = LIBMUSCLE_Data_create_grid(ar)

    ! byte array
    d12 = LIBMUSCLE_Data_create_byte_array(bytes)


As you can see, sending complex data types with MUSCLE is a bit more difficult
in Fortran than in Python, but it is not too burdensome.

If you want to send large amounts of data, then use grids (arrays) as much as
possible. Lists and dicts are more flexible, but they are also slower and take
up much more memory. In particular, a set of objects (agents for example) is
better sent as a dict of arrays (with one 1D array for each attribute) than as a
list of dicts. (The latter will probably work up to 1 million objects or so, but
it will still be slower.)

Handling errors
```````````````

In some cases when you're asking MUSCLE 3 to do something, things can go wrong.
For example, receiving a message from a port does not work if the port is not
connected (unless you passed a default message, see
:f:func:`LIBMUSCLE_Instance_receive`). You cannot extract a real number from a
:f:type:`LIBMUSCLE_DataConstRef` which contains an integer, or obtain a setting
as a character string if the user set it to a logical value. There are two ways
of dealing with these situations: checking in advance if it's going to work and
acting accordingly, or just trying and checking whether you were successful. If
you do neither (which is easiest), MUSCLE 3 will print an error message and stop
your program if an error occurs (which is actually usually what you want).

How to check in advance whether an operation can work depends on what you want
to do. If you want to receive a message, you could use
:f:func:`LIBMUSCLE_Instance_is_connected` to check that the port is connected
before trying to receive. For settings that can be of different types, you can
use ``LIBMUSCLE_Instance_is_setting_a_<type>``, and for data there is
``LIBMUSCLE_DataConstRef_is_a_<type>``.

This does not work in all cases however, for example if you have received a
dictionary that may or may not have a particular key. In that case, (and
anywhere else you think it's more convenient) it's better to just try, and
handle any errors if they occur. This is done by adding one or two extra
arguments to the function call that you expect may fail, like this:

.. code-block:: fortran

    type(LIBMUSCLE_DataConstRef) :: rdata, item
    integer :: err_code
    character(len=:), allocatable :: err_msg

    logical :: value

    ! Receiving rdata omitted

    ! If there is a key 'key' and its value is a logical, set variable
    ! value to that, otherwise set it to .true..
    item = LIBMUSCLE_DataConstRef_get_item(rdata, 'key', err_code, err_msg)
    if (err_code /= LIBMUSCLE_success) then
        print *, err_msg
        ! Need to deallocate the message
        deallocate(err_msg)
        value = .true.
    else
        value = LIBMUSCLE_DataConstRef_as_logical(item, err_code)
        if (err_code /= LIBMUSCLE_success) then
            value = .true.
        end if
    end if


Note that if an error occurs, an error message will be set and the variable
must be deallocated after you're done using it. If no error occurs, the variable
will remain unset, and no deallocation is needed. Passing a variable for an
error message is optional, you can also only request the error code, as shown.
You cannot get only a message, and no code.

Valid error code values are ``LIBMUSCLE_success`` (no error),
``LIBMUSCLE_domain_error`` (incorrect input to a function),
``LIBMUSCLE_out_of_range`` (invalid index), ``LIBMUSCLE_logic_error`` (you tried
to do something that doesn't make sense), ``LIBMUSCLE_runtime_error`` (general
error, e.g.  you're trying to get the length of a scalar port, which does not
make sense), ``LIBMUSCLE_bad_cast`` (value doesn't have the type you're trying
to retrieve it as). See the :ref:`api-docs-fortran` for the function you are
calling to see which error code it will return when exactly.

(Implementation note: which error is returned when is a bit messy still, and
will be cleaned up in a future version.)


Diffusion model
```````````````

If you've studied the above carefully, and have seen the Python version of the
diffusion model, then you should now be able to understand the Fortran diffusion
model below:

.. literalinclude:: examples/fortran/diffusion.f03
  :caption: ``docs/source/examples/fortran/diffusion.f03``
  :language: fortran


In the examples directory, you will find a handy script called
``reaction_diffusion_fortran.sh``, which runs the Fortran reaction-diffusion
model locally. This script launches the MUSCLE Manager and an instance of each
submodel, then waits for the simulation to complete. At the top of this page
there are instructions on how to run the examples.

Also in this directory are Fortran versions of a Monte Carlo sampler and a
round-robin load balancer, plus a script to launch this extended example. See
the `Uncertainty Quantification` section to learn what those are for.

MUSCLE 3 Examples
=================

This directory contains examples for MUSCLE 3. Once you've `installed MUSCLE 3
<https://muscle3.readthedocs.io/en/latest/installing.html>`_, you can build them
by running Make from this directory, like this:

.. code-block: bash

    examples$ MUSCLE3_HOME=/path/to/muscle3 make


This will build the C++ and Fortran examples, and create a virtualenv for
running the Python examples and the manager. You can also build for each
language separately, using

.. code-block: bash

    examples$ MUSCLE3_HOME=/path/to/muscle3 make python
    examples$ MUSCLE3_HOME=/path/to/muscle3 make cpp
    examples$ MUSCLE3_HOME=/path/to/muscle3 make fortran


Note that you need to run the Python build in order to be able to run the C++ or
Fortran examples, as it does the set-up required to run the manager, which you
always need to run a MUSCLE simulation.

Once you've built the examples, you can run them using the shell scripts in this
directory, e.g.

.. code-block: bash

    examples$ MUSCLE3_HOME=/path/to/muscle3 ./reaction_diffusion_cpp.sh


You can also export the ``MUSCLE3_HOME`` variable to avoid having to type it all
the time:

.. code-block: bash

    examples$ export MUSCLE3_HOME=/path/to/muscle3
    examples$ ./reaction_diffusion_cpp.sh
======================
LibMuscle architecture
======================

LibMuscle is the MUSCLE 3 client library. Kernels link with it to get access to
functionality that lets them be part of a MUSCLE 3 simulation. As kernels can be
written in many languages, a version of LibMuscle will have to be made for each
such language.

LibMuscle has two core responsibilities: it manages configuration information,
and it takes care of communication with other instances. It can also take care
of execution of the kernel, allowing users to create kernels by simply
implementing a set of functions. Another ancillary task is providing centralised
logging.

The user uses the LibMuscle API to talk to it. It has a number of sub-APIs:

* Coordination API
* Configuration API
* Communication API
* Execution API
* Logging API

These API components contain their own logic, and use functionality
implemented by three main LibMuscle components:

* MMP client
* Configuration Store
* Communication Engine


----------------
Coordination API
----------------

The Coordination API consists of a single function, which
initialises MUSCLE 3, contacts the manager, and performs the
instance side of the MUSCLE initialisation process. This function
must always be called first, and should take the form of a
constructor of a central MUSCLE 3 object that all other communcation
goes through.

The init function is passed (if necessary) the command line
arguments, and reads and removes the MUSCLE 3-specific ones. These
tell it where the Manager is, which it uses to initialise the MMP
Client. It is also passed a description of the kernel, which it
passes to the Manager via the MMP Client on registration.

The Coordination API uses the MMP Client to connect to the Manager,
and stores the information obtained via the Manager in the
Configuration Store and the Communication Engine.


-----------------
Configuration API
-----------------

Using the Configuration API, the client can test whether a given
parameter is set, which type it has (in statically typed languages),
and request its value. The Configuration API obtains this
information from the Configuration Store.


-----------------
Communication API
-----------------

The Communication API is used to send and receive messages. It
contains functions for sending and receiving messages, as well as
support for packing complex data structures for transport (using
MessagePack). It uses the Communication Engine to do the actual
work.


-------------
Execution API
-------------

The Execution API implements the Submodel Execution Loop, as well as
execution loops for other types of kernels. It allows the user to
define their model as a class with a set of member functions, one
for each operator. The library then takes care of calling these
functions in the correct order, and passing them messages. See
`LittleMuscle`_'s `submodel.py`_ for a prototype.


-----------
Logging API
-----------

The Logging API is used by the client to log information about its
progress, errors, etc. If something goes wrong in a distributed
system, it is often difficult to find out what happened and why.
Having log output easily available in a single location makes this
much easier, so MUSCLE 3 provides this convenience. Furthermore, the
logging system can be used to obtain basic performance information,
which can be used by the planner to allocate resources more
efficiently.


----------
MMP Client
----------

The MMP Client communicates with the Manager. Its code is partially
generated from the gRPC specification of the MUSCLE Manager
Protocol. It provides functions to the other modules for registering
the instance, getting peers, getting configuration information, and
sending log messages, and implements these by communicating with the
Manager. It is initialised with the network location of the Manager.


-------------------
Configuration Store
-------------------

The Configuration Store contains information about the configuration
of this instance. Like in the Manager, the basic container of
configuration information is a dictionary structure that maps the
name of a model parameter (a string) to a value, which is a string,
integer, floating point number, or nil (None, null).

Unlike the Manager, the Configuration Store contains only
configuration information for its own instance, but it maintains a
stack of these dictionaries. Lookup of a parameter value starts at
the dictionary on the top of the stack, and falls through if the
parameter is not found. A parameter can be set to nil explicitly
(thus overriding any value lower in the stack), or if it is not
present in any of the dictionaries, then it will be nil by default.

When an instance calls the Coordination API (see below) to register
itself, it sends a dictionary with parameters and their default
values. This becomes the bottom layer of its stack.


--------------------
Communication Engine
--------------------

The Communication Engine component of LibMuscle takes care of
communication with other instances. It keeps a list of endpoints,
with for each sending endpoint the location of the peer on the other
side of the conduit connected to the endpoint. This information is
passed to it by the Coordination API as part of starting up. It
sends and receives messages on these end points (as directed by the
Communication API), but does not know anything about the content of
the messages.

Actual communication is done by a TCP submodule, which knows how to
set up a connection to a peer given its network location, and
performs actual communication. In the future, we envision submodules
for other protocols as well, e.g. MPI.


.. _`submodel.py`: https://github.com/multiscale/littlemuscle/blob/develop/littlemuscle/submodel.py
.. _`LittleMuscle`: https://github.com/multiscale/littlemuscle

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   introduction.rst
   manager.rst
   libmuscle.rst
=========================
MUSCLE 3 Technical Design
=========================

------------
Introduction
------------

Computational science involves the development of simulation models of complex
natural systems. Often, several processes in such a system are relevant to the
research question being posed, and in many cases these processes act at
different scales in time and space. To model such a system, several submodels,
one for each process, are combined together into a multi-scale model.

To implement a multiscale model in software, some kind of subsystem is required
that is responsible for *coupling* the models. This subsystem ensures that the
submodels can exchange information when needed, that the model configuration
and parameters are available to the submodels, and in case of dynamic use of
compute resources, that the correct set of submodels is running when needed.
Thus, the coupling library provides coordination and communication in addition
to the computing done by the submodels.

MUSCLE 3 is the third major version of the Multiscale Coupling Library and
Environment. MUSCLE is grounded in a solid theoretical basis provided by the
Multiscale Modelling Language (MML), which is being extended to include a
mechanism for running large, dynamically sized sets of submodels. Using MUSCLE
3, submodels written in different programming languages can be coupled
together, with data exchanged via a network to facilitate distributed
execution. As such, MUSCLE 3 is a core component of the Multiscale Modelling
and Simulation Framework (MMSF) for high-performance distributed computing.


-------------------------------
Anatomy of a MUSCLE3 simulation
-------------------------------

A MUSCLE 3 simulation run consists of a set of kernels (which include submodels
as well as other components), which are separate programs running independently
and exchanging information while doing so, and a manager that coordinates their
activity. The simulation starts with a start-up phase, during which the kernels
are connected to each other and settings are exchanged. Once all kernels are in
place and able to communicate, the execution phase starts. During the execution
phase, the kernels operate independently, exchanging information when needed.
The manager remains active to handle dynamic instantiation of new kernels, as
well as logging.

For execution across multiple high-performance computing (HPC) resources, a
separate program is usually required that forwards messages between the
clusters, due to constraints imposed by the network architecture (i.e.
firewalls that are in the way). This is the MUSCLE Transport Overlay (MTO).
Details about the MTO will be added to this document later.

MUSCLE 3 provides four components for building a simulation. First is the
central manager, which keeps track of where kernels are executing and informs
them of each other's location, so that communication channels can be set up. It
also provides a central facility for kernel configuration, it does centralised
logging, and it starts extra kernels on demand.

LibMuscle is a library that submodels and other kernels link with, and which
provides them with all the functionality needed to be part of a MUSCLE 3
simulation. LibMuscle will be available in several programming languages, so
that existing kernels can be used as-is, and so that users are free to choose a
suitable language.

Besides these two pieces of software, MUSCLE 3 comprises two communication
protocols. The first is the MUSCLE Manager Protocol (MMP), which is used by
kernels to communicate with the manager. The second is the MUSCLE Data Protocol
(MDP), which is used for exchange of data between kernels.


----------------------
Simulation walkthrough
----------------------

In this section, a MUSCLE 3 simulation is described in a bit more detail. As
described above, a MUSCLE simulation consists of a start-up phase followed by
an execution phase. We assume that the simulation is initiated by a launcher,
some external program that starts up the manager and an initial set of kernel
instances on one or more machines.

On start-up, the manager begins by reading a configuration file that describes
the model. This file is in the yMML format, a YAML-based representation of the
Multiscale Modelling Language. It contains a description of the submodels, the
connection topology between them, and a set of parameter values to run the
various submodels with. After having read the configuration file, the manager
waits for the instances to contact it.

When the instances start up, they contact the manager (having been passed its
network location via the command line) and register themselves, by providing
their name, a description of their communication endpoints, and the network
location they listen on). After registration, each instance will contact the
manager again to request the network locations of the peers with which it
should communicate. The manager passes the instance this information, which is
a combination of the contents of the yMML file and the information sent to it
by the other instances. Next, the instance will request a set of parameter
values from the manager, which returns the ones read from the yMML file. After
this, the kernel instance will start its execution.

During execution, an instance will normally execute in a loop, sending and
receiving messages before, during, and after loop execution. It does this in
terms of a set of endpoints. An endpoint can be either a sending endpoint, or a
receiving endpoint. These endpoints are a property of each kernel. Information
about which other instance each endpoint is connected to, is managed by
LibMuscle, and not seen or used by the user's code.

In some cases, dynamic instantiation of a kernel may be needed. This is
achieved by including a dynamic conduit in the model description. A dynamic
conduit connects a sending endpoint of a running instance, to a receiving
endpoint of an instance that is to be started on demand. When the sending
instance sends a message on the sending endpoint (which must happen through
LibMuscle), LibMuscle will contact the manager and ask it to provide a suitable
receiving submodel instance. The manager looks up what kind of instance to
provide, creates or reuses an instance, and returns the network location of the
new receiving instance. LibMuscle then sends it the message.

During execution, instances generate log messages that can be used for checking
progress and for debugging. If the log level of the message is at or above the
threshold set for local logging, it will be written to a log file on the
machine where the instance is running. If the log level is at or above the
threshold for central logging, then it will be sent to the manager for
recording in the central log file.


-----------------
Model description
-----------------

The description of a model comprises a description of its kernels, the
instances of those kernels that the running model is made of, and the
conduits that carry information between the instances. This information is read
by the manager from a configuration file in yMML format. In this section, the
model description is described in more detail.


Kernels
-------

A kernel is (in this context) a computer program that uses MUSCLE to
communicate with the outside world.  As a kernel executes, it runs through a
succession of *Operators*. The Multiscale Modelling Language specifies three
kinds of kernels: Submodels, Mappers, and Proxies, each with their own set of
operators and execution order. MUSCLE supports these specifically, but also
provides a generic interface that can be used if some other kind of kernel is
required.

A Submodel kernel has five operators, which are run through according to the
Submodel Execution Loop (SEL). In the SEL, a submodel, after startup, first
receives one or more messages that determine its initial state (f_init). Then,
it sends one or more messages describing (parts of) this initial state (O_i).
Next, it performs a state update, possibly receiving messages that provide
additional input (S). Next, it updates its boundary conditions, again possibly
receiving messages (B). Then, it loops back to sending messages describing the
new state.  After the boundary condition update, the submodel may leave the
loop if it is done, send one final set of messages describing its final state
(O_f), and end this present run. If needed, the submodel (as with any kernel)
may then be restarted from the beginning.

A Mapper has only a single operator, M, a set of receiving endpoints, and a set
of sending endpoints. A mapper with a single receiving endpoint and several
sending endpoints is a fan-out mapper; a mapper with several receiving
endpoints and a single sending endpoint is a fan-in mapper. A mapper first
receives one message on each of its receiving endpoints, and then sends one
message on each of its sending endpoints. The mapper may then be restarted, and
receive-and-send again. A mapper must send a message on each of its sending
endpoints whenever it receives messages on each of its receiving endpoints.

The Proxy is a new kind of kernel in MUSCLE 3 with a single operator, P, and
four sets of conduits. It will be described in more detail later.


Endpoints
---------

Endpoints are used by kernels to communicate with the outside world. They have
a unique (for that kernel) name, an associated operator, and a data type.
Operators may come with restrictions on endpoints, e.g. the submodel's S
operator may only receive messages. A sending endpoint with associated operator
S is therefore invalid.


Instances
---------

An instance is a process, a kernel running on some computer. Models may contain
one or more (and even very many) instances of the same kernel, for example with
each instance calculating some part of the spatial domain. The model
description contains a list of all the instances comprising the model. Each
instance definition specifies which kernel it is an instance of, and it
contains configuration information for the kernel.

The configuration consists of three parts: the space and time scales of the
modelled process, MUSCLE 3 built-in settings, and kernel parameters.

Scales in MUSCLE are defined by their grain (step or cell size), and extent
(total size). For a kernel operating on a grid, the space scales specify the
grid size. For kernels with a repeated solving step, the time scale specifies
the size of the time step and the overall duration of the simulation. Spatial
and temporal scales should be chosen with care, dependent on the spatial and
temporal characteristics of the modelled phenomenon.

MUSCLE 3 built-in settings are settings that are used by LibMuscle, not by the
user-written kernel code. These include e.g. configuration of the logging
subsystem.

Finally, kernel parameters are defined by the maker of the kernel, and can be
any kind of model parameters or configuration.


Conduits
--------

Conduits connect instances, allowing them to exchange messages. More
specifically, a conduit connects a sending endpoint on a given instance to a
receiving endpoint on another instance.
====================
Manager architecture
====================

The MUSCLE 3 manager consists of several modules, each with a specific
responsibility. Together, they perform all the tasks required of the manager.
These components are:

* the MMP Server, which listens for connections from kernels and responds to
  them,
* the Logger, which gets log messages and writes them to a log file,
* the Instance Registry, which stores information about instances,
* the Configuration Store, which stores kernel configuration information,
* the Topology Store, which stores information on connections between instances,
* LibMML, which reads yMML files and extracts information about kernels,
  instances, and conduits,
* Main, which creates the other components and performs start-up tasks.


----------
MMP Server
----------

The MMP Server is a gRPC server. It listens for requests from instances, and
services them by calling on functionality in the Logger, Instance Registry,
Configuration Store and Topology Store. These requests include:

* Log Requests
* Registration Requests
* Configuration Requests
* Peer Requests
* Deregistration Requests
* Dynamic Peer Requests (will be described later)

A Log Request is a request to log a message. Instances send these for log
messages with a high enough log level (e.g. critical errors, log messages used
for performance tracking). Log messages are structured (see below at the
description of the Logger).

A Registration Request is how an instance makes itself known to the manager.
The instance sends its name, a list of its endpoints, and the network location
it is listening on (a string). The server sends this information to the
Instance Registry, and returns nothing.

A Configuration Request is sent by an instance. It sends its instance name, and
receives a dictionary of parameters and their values (see below under
Configuration Store).

A Peer Request is sent by an instance after registering itself. It sends its
instance name, and receives either an error message, if its peers have not yet
registered themselves, or a list of records containing

* sender endpoint name
* slot number
* target network location
* receiver name
* receiver endpoint name

with one such record for each endpoint the instance registered in its
Registration Request. These are all strings, with the exception of the slot
number, which is an integer that is always equal to 0 (this will change when
vector conduits are added, to be described later). The MMP Server obtains this
information from the Instance Registry (to get the sending endpoints of the
requesting instance), the Topology Store (to get the outgoing conduit and the
name of the receiving instance), and again from the Instance Registry (to get
the network location of the receiving instance).

A Deregistration Request is sent by an instance when it shuts down. It contains
the instance's name, and causes the MMP Server to remove the instance from the
Instance Registry.


------
Logger
------

The logger component takes log messages, which consist of

* instance name
* operator
* time stamp
* log level
* message text

and writes them to a log file in text format. Initially, this can just be a
fixed muscle3.log in the current directory, with configurability of the
location and other niceties added later. The log level should be configurable
though. In the future, we will want the ability to write (instance name,
operator, time stamp) for a given log level to a second file in
machine-readable format (e.g. JSON).


-----------------
Instance Registry
-----------------

The Instance Registry stores information about the currently running instances.
It is a simple in-memory database. Its user is the MMP Server, which adds,
requests and removes entries. The Instance Registry stores the following
information for each instance:

* instance name
* network location
* list of endpoints

where each endpoint has a

* name
* operator
* data type.

New instances can be added to the Instance Registry, instances can be deleted
given their name, and instances can be retrieved by name.


-------------------
Configuration Store
-------------------

The Configuration Store contains information about the configuration of the
model that is being run. The basic container of configuration information is a
dictionary structure that maps the name of a model parameter (a string) to a
value, which is a string, integer, floating point number, or nil (None, null).

The Configuration Store maintains such a dictionary for each instance. This
dictionary is read from the yMML file on start-up by the Main component and
passed to the Configuration Store. Later, the MMP Server will request the
dictionary for a given instance in order to service a request from that
instance.


--------------
Topology Store
--------------

The Topology Store stores a list of conduits. Each conduit has the sending
instance's name, the sending endpoint name, the receiving instance's name, and
the receiving endpoint name. These are also read from the yMML file. This
information is used by the MMP Server to service Peer Requests.


------
LibMML
------

LibMML is a Python library for reading and writing yMML files, which are
YAML-based representations of MML graphs. The yMML file format still needs to
be designed; it will be based on xMML, with some extensions to support the
ComPat patterns.

LibMML will be based on the Python YAML library, with a shim layer to add some
syntactic sugar to the file format, to make it easier to type by hand. LibMML
is not a part of MUSCLE 3, but will be used by it.

yMML files will contain information on kernels, instances, the connections
between them, and the parameters to use them with. Information on how to deploy
the simulation will probably be kept separate, although the two are of course
related.


----
Main
----

The Main component of the MUSCLE 3 Manager performs start-up of the manager. It
parses the command line, setting the log level as given, and reading the yMML
file location. It then constructs all the other components, loads the yMML file
and stores the information into the Configuration Store and Topology Store. It
then starts the MMP Server and passes it control.
