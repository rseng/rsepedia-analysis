<div align="center">
  <img src="https://raw.githubusercontent.com/zarr-developers/community/master/logos/logo2.png"><br>
</div>

# Zarr

<table>
<tr>
  <td>Latest Release</td>
  <td>
    <a href="https://pypi.org/project/zarr/">
    <img src="https://badge.fury.io/py/zarr.svg" alt="latest release" />
    </a>
  </td>
</tr>
  <td></td>
  <td>
    <a href="https://anaconda.org/anaconda/zarr/">
    <img src="https://anaconda.org/conda-forge/zarr/badges/version.svg" alt="latest release" />
    </a>
</td>
</tr>
<tr>
  <td>Package Status</td>
  <td>
		<a href="https://pypi.org/project/zarr/">
		<img src="https://img.shields.io/pypi/status/zarr.svg" alt="status" />
		</a>
  </td>
</tr>
<tr>
  <td>License</td>
  <td>
    <a href="https://github.com/zarr-developers/zarr-python/blob/master/LICENSE">
    <img src="https://img.shields.io/pypi/l/zarr.svg" alt="license" />
    </a>
</td>
</tr>
<tr>
  <td>Build Status</td>
  <td>
    <a href="https://travis-ci.org/zarr-developers/zarr-python">
    <img src="https://travis-ci.org/zarr-developers/zarr-python.svg?branch=master" alt="travis build status" />
    </a>
  </td>
</tr>
<tr>
  <td>Coverage</td>
  <td>
    <a href="https://codecov.io/gh/zarr-developers/zarr-python">
    <img src="https://codecov.io/gh/zarr-developers/zarr-python/branch/master/graph/badge.svg"/ alt="coverage">
    </a>
  </td>
</tr>
<tr>
  <td>Downloads</td>
  <td>
    <a href="https://zarr.readthedocs.io">
    <img src="https://pepy.tech/badge/zarr" alt="pypi downloads" />
    </a>
  </td>
</tr>
<tr>
	<td>Gitter</td>
	<td>
		<a href="https://gitter.im/zarr-developers/community">
		<img src="https://badges.gitter.im/zarr-developers/community.svg" />
		</a>
	</td>
</tr>
<tr>
	<td>Citation</td>
	<td>
		<a href="https://doi.org/10.5281/zenodo.3773450">
			<img src="https://zenodo.org/badge/DOI/10.5281/zenodo.3773450.svg" alt="DOI">
		</a>
	</td>
</tr>

</table>

## What is it?

Zarr is a Python package providing an implementation of compressed, chunked, N-dimensional arrays, designed for use in parallel computing. See the [documentation](https://zarr.readthedocs.io) for more information.

## Main Features

- [**Create**](https://zarr.readthedocs.io/en/stable/tutorial.html#creating-an-array) N-dimensional arrays with any NumPy `dtype`.
- [**Chunk arrays**](https://zarr.readthedocs.io/en/stable/tutorial.html#chunk-optimizations) along any dimension.
- [**Compress**](https://zarr.readthedocs.io/en/stable/tutorial.html#compressors) and/or filter chunks using any NumCodecs codec.
- [**Store arrays**](https://zarr.readthedocs.io/en/stable/tutorial.html#tutorial-storage) in memory, on disk, inside a zip file, on S3, etc...
- [**Read**](https://zarr.readthedocs.io/en/stable/tutorial.html#reading-and-writing-data) an array [**concurrently**](https://zarr.readthedocs.io/en/stable/tutorial.html#parallel-computing-and-synchronization) from multiple threads or processes.
- Write to an array concurrently from multiple threads or processes.
- Organize arrays into hierarchies via [**groups**](https://zarr.readthedocs.io/en/stable/tutorial.html#groups).

## Where to get it

Zarr can be installed from PyPI using `pip`:

```bash
pip install zarr
```

or via `conda`:

```bash
conda install -c conda-forge zarr
```

For more details, including how to install from source, see the [installation documentation](https://zarr.readthedocs.io/en/stable/#installation).
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at zarr.conduct@gmail.com. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
For bug reports, please follow the template below. For enhancement proposals, feel free
to use whatever template makes sense (major new features should be discussed in the
Zarr specifications repository https://github.com/zarr-developers/zarr-specs).

#### Minimal, reproducible code sample, a copy-pastable example if possible

```python
# Your code here

```

#### Problem description

Explain why the current behavior is a problem, what the expected output/behaviour 
is, and why the expected output/behaviour is a better solution.

#### Version and installation information

Please provide the following:

* Value of ``zarr.__version__``
* Value of ``numcodecs.__version__``
* Version of Python interpreter
* Operating system (Linux/Windows/Mac)
* How Zarr was installed (e.g., "using pip into virtual environment", or "using conda")

Also, if you think it might be relevant, please provide the output from ``pip freeze`` or
``conda env export`` depending on which was used to install Zarr.
Contributing
============

Please see the [project documentation](http://zarr.readthedocs.io/en/stable/contributing.html) for information about contributing to Zarr.

[Description of PR]

TODO:
* [ ] Add unit tests and/or doctests in docstrings
* [ ] Add docstrings and API docs for any new/modified user-facing classes and functions
* [ ] New/modified features documented in docs/tutorial.rst
* [ ] Changes documented in docs/release.rst
* [ ] GitHub Actions have all passed
* [ ] Test coverage is 100% (Codecov passes)
Contributing to Zarr
====================

Zarr is a community maintained project. We welcome contributions in the form of bug
reports, bug fixes, documentation, enhancement proposals and more. This page provides
information on how best to contribute.

Asking for help
---------------

If you have a question about how to use Zarr, please post your question on
StackOverflow using the `"zarr" tag <https://stackoverflow.com/questions/tagged/zarr>`_.
If you don't get a response within a day or two, feel free to raise a `GitHub issue
<https://github.com/zarr-developers/zarr-python/issues/new>`_ including a link to your StackOverflow
question. We will try to respond to questions as quickly as possible, but please bear
in mind that there may be periods where we have limited time to answer questions
due to other commitments.

Bug reports
-----------

If you find a bug, please raise a `GitHub issue
<https://github.com/zarr-developers/zarr-python/issues/new>`_. Please include the following items in
a bug report:

1. A minimal, self-contained snippet of Python code reproducing the problem. You can
   format the code nicely using markdown, e.g.::


    ```python
    import zarr
    g = zarr.group()
    # etc.
    ```

2. An explanation of why the current behaviour is wrong/not desired, and what you
   expect instead.

3. Information about the version of Zarr, along with versions of dependencies and the
   Python interpreter, and installation information. The version of Zarr can be obtained
   from the ``zarr.__version__`` property. Please also state how Zarr was installed,
   e.g., "installed via pip into a virtual environment", or "installed using conda".
   Information about other packages installed can be obtained by executing ``pip freeze``
   (if using pip to install packages) or ``conda env export`` (if using conda to install
   packages) from the operating system command prompt. The version of the Python
   interpreter can be obtained by running a Python interactive session, e.g.::

    $ python
    Python 3.6.1 (default, Mar 22 2017, 06:17:05)
    [GCC 6.3.0 20170321] on linux

Enhancement proposals
---------------------

If you have an idea about a new feature or some other improvement to Zarr, please raise a
`GitHub issue <https://github.com/zarr-developers/zarr-python/issues/new>`_ first to discuss.

We very much welcome ideas and suggestions for how to improve Zarr, but please bear in
mind that we are likely to be conservative in accepting proposals for new features. The
reasons for this are that we would like to keep the Zarr code base lean and focused on
a core set of functionalities, and available time for development, review and maintenance
of new features is limited. But if you have a great idea, please don't let that stop
you from posting it on GitHub, just please don't be offended if we respond cautiously.

Contributing code and/or documentation
--------------------------------------

Forking the repository
~~~~~~~~~~~~~~~~~~~~~~

The Zarr source code is hosted on GitHub at the following location:

* `https://github.com/zarr-developers/zarr-python <https://github.com/zarr-developers/zarr-python>`_

You will need your own fork to work on the code. Go to the link above and hit
the "Fork" button. Then clone your fork to your local machine::

    $ git clone git@github.com:your-user-name/zarr.git
    $ cd zarr
    $ git remote add upstream git@github.com:zarr-developers/zarr-python.git

Creating a development environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To work with the Zarr source code, it is recommended to set up a Python virtual
environment and install all Zarr dependencies using the same versions as are used by
the core developers and continuous integration services. Assuming you have a Python
3 interpreter already installed, and have also installed the virtualenv package, and
you have cloned the Zarr source code and your current working directory is the root of
the repository, you can do something like the following::

    $ mkdir -p ~/pyenv/zarr-dev
    $ python -m venv ~/pyenv/zarr-dev
    $ source ~/pyenv/zarr-dev/bin/activate
    $ pip install -r requirements_dev_minimal.txt -r requirements_dev_numpy.txt
    $ pip install -e .

To verify that your development environment is working, you can run the unit tests::

    $ pytest -v zarr

Creating a branch
~~~~~~~~~~~~~~~~~

Before you do any new work or submit a pull request, please open an issue on GitHub to
report the bug or propose the feature you'd like to add.

It's best to synchronize your fork with the upstream repository, then create a
new, separate branch for each piece of work you want to do. E.g.::

    git checkout master
    git fetch upstream
    git rebase upstream/master
    git push
    git checkout -b shiny-new-feature
    git push -u origin shiny-new-feature

This changes your working directory to the 'shiny-new-feature' branch. Keep any changes in
this branch specific to one bug or feature so it is clear what the branch brings to
Zarr.

To update this branch with latest code from Zarr, you can retrieve the changes from
the master branch and perform a rebase::

    git fetch upstream
    git rebase upstream/master

This will replay your commits on top of the latest Zarr git master. If this leads to
merge conflicts, these need to be resolved before submitting a pull request.
Alternatively, you can merge the changes in from upstream/master instead of rebasing,
which can be simpler::

    git fetch upstream
    git merge upstream/master

Again, any conflicts need to be resolved before submitting a pull request.

Running the test suite
~~~~~~~~~~~~~~~~~~~~~~

Zarr includes a suite of unit tests, as well as doctests included in
function and class docstrings and in the tutorial and storage
spec. The simplest way to run the unit tests is to activate your
development environment (see `creating a development environment`_ above)
and invoke::

    $ pytest -v zarr

Some tests require optional dependencies to be installed, otherwise
the tests will be skipped. To install all optional dependencies, run::

    $ pip install -r requirements_dev_optional.txt

To also run the doctests within docstrings (requires optional
dependencies to be installed), run::

    $ pytest -v --doctest-plus zarr

To run the doctests within the tutorial and storage spec (requires
optional dependencies to be installed), run::

    $ python -m doctest -o NORMALIZE_WHITESPACE -o ELLIPSIS docs/tutorial.rst docs/spec/v2.rst

Note that some tests also require storage services to be running
locally. To run the Azure Blob Service storage tests, run an Azure
storage emulator (e.g., azurite) and set the environment variable
``ZARR_TEST_ABS=1``. If you're using Docker to run azurite, start the service with::

    docker run --rm -p 10000:10000 mcr.microsoft.com/azure-storage/azurite azurite-blob --loose --blobHost 0.0.0.0 

To run the Mongo DB storage tests, run a Mongo
server locally and set the environment variable ``ZARR_TEST_MONGO=1``.
To run the Redis storage tests, run a Redis server locally on port
6379 and set the environment variable ``ZARR_TEST_REDIS=1``.

All tests are automatically run via GitHub Actions for every pull
request and must pass before code can be accepted. Test coverage is
also collected automatically via the Codecov service, and total
coverage over all builds must be 100% (although individual builds
may be lower due to Python 2/3 or other differences).

Code standards
~~~~~~~~~~~~~~

All code must conform to the PEP8 standard. Regarding line length, lines up to 100
characters are allowed, although please try to keep under 90 wherever possible.
Conformance can be checked by running::

    $ flake8 --max-line-length=100 zarr

Test coverage
~~~~~~~~~~~~~

Zarr maintains 100% test coverage under the latest Python stable release (currently
Python 3.8). Both unit tests and docstring doctests are included when computing
coverage. Running ``tox -e py38`` will automatically run the test suite with coverage
and produce a coverage report. This should be 100% before code can be accepted into the
main code base.

When submitting a pull request, coverage will also be collected across all supported
Python versions via the Codecov service, and will be reported back within the pull
request. Codecov coverage must also be 100% before code can be accepted.

Documentation
~~~~~~~~~~~~~

Docstrings for user-facing classes and functions should follow the
`numpydoc
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
standard, including sections for Parameters and Examples. All examples
should run and pass as doctests under Python 3.8. To run doctests,
activate your development environment, install optional requirements,
and run::

    $ pytest -v --doctest-plus zarr

Zarr uses Sphinx for documentation, hosted on readthedocs.org. Documentation is
written in the RestructuredText markup language (.rst files) in the ``docs`` folder.
The documentation consists both of prose and API documentation. All user-facing classes
and functions should be included in the API documentation, under the ``docs/api``
folder. Any new features or important usage information should be included in the
tutorial (``docs/tutorial.rst``). Any changes should also be included in the release
notes (``docs/release.rst``).

The documentation can be built locally by running::

    $ tox -e docs

The resulting built documentation will be available in the ``.tox/docs/tmp/html`` folder.

Development best practices, policies and procedures
---------------------------------------------------

The following information is mainly for core developers, but may also be of interest to
contributors.

Merging pull requests
~~~~~~~~~~~~~~~~~~~~~

Pull requests submitted by an external contributor should be reviewed and approved by at least
one core developers before being merged. Ideally, pull requests submitted by a core developer
should be reviewed and approved by at least one other core developers before being merged.

Pull requests should not be merged until all CI checks have passed (GitHub Actions
Codecov) against code that has had the latest master merged in.

Compatibility and versioning policies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because Zarr is a data storage library, there are two types of compatibility to
consider: API compatibility and data format compatibility.

API compatibility
"""""""""""""""""

All functions, classes and methods that are included in the API
documentation (files under ``docs/api/*.rst``) are considered as part of the Zarr **public API**,
except if they have been documented as an experimental feature, in which case they are part of
the **experimental API**.

Any change to the public API that does **not** break existing third party
code importing Zarr, or cause third party code to behave in a different way, is a
**backwards-compatible API change**. For example, adding a new function, class or method is usually
a backwards-compatible change. However, removing a function, class or method; removing an argument
to a function or method; adding a required argument to a function or method; or changing the
behaviour of a function or method, are examples of **backwards-incompatible API changes**.

If a release contains no changes to the public API (e.g., contains only bug fixes or
other maintenance work), then the micro version number should be incremented (e.g.,
2.2.0 -> 2.2.1). If a release contains public API changes, but all changes are
backwards-compatible, then the minor version number should be incremented
(e.g., 2.2.1 -> 2.3.0). If a release contains any backwards-incompatible public API changes,
the major version number should be incremented (e.g., 2.3.0 -> 3.0.0).

Backwards-incompatible changes to the experimental API can be included in a minor release,
although this should be minimised if possible. I.e., it would be preferable to save up
backwards-incompatible changes to the experimental API to be included in a major release, and to
stabilise those features at the same time (i.e., move from experimental to public API), rather than
frequently tinkering with the experimental API in minor releases.

Data format compatibility
"""""""""""""""""""""""""

The data format used by Zarr is defined by a specification document, which should be
platform-independent and contain sufficient detail to construct an interoperable
software library to read and/or write Zarr data using any programming language. The
latest version of the specification document is available from the :ref:`spec` page.

Here, **data format compatibility** means that all software libraries that implement a
particular version of the Zarr storage specification are interoperable, in the sense
that data written by any one library can be read by all others. It is obviously
desirable to maintain data format compatibility wherever possible. However, if a change
is needed to the storage specification, and that change would break data format
compatibility in any way, then the storage specification version number should be
incremented (e.g., 2 -> 3).

The versioning of the Zarr software library is related to the versioning of the storage
specification as follows. A particular version of the Zarr library will
implement a particular version of the storage specification. For example, Zarr version
2.2.0 implements the Zarr storage specification version 2. If a release of the Zarr
library implements a different version of the storage specification, then the major
version number of the Zarr library should be incremented. E.g., if Zarr version 2.2.0
implements the storage spec version 2, and the next release of the Zarr library
implements storage spec version 3, then the next library release should have version
number 3.0.0. Note however that the major version number of the Zarr library may not
always correspond to the spec version number. For example, Zarr versions 2.x, 3.x, and
4.x might all implement the same version of the storage spec and thus maintain data
format compatibility, although they will not maintain API compatibility. The version number
of the storage specification that is currently implemented is stored under the
``zarr.meta.ZARR_FORMAT`` variable.

Note that the Zarr test suite includes a data fixture and tests to try and ensure that
data format compatibility is not accidentally broken. See the
:func:`test_format_compatibility` function in the :mod:`zarr.tests.test_storage` module
for details.

When to make a release
~~~~~~~~~~~~~~~~~~~~~~

Ideally, any bug fixes that don't change the public API should be released as soon as
possible. It is fine for a micro release to contain only a single bug fix.

When to make a minor release is at the discretion of the core developers. There are no
hard-and-fast rules, e.g., it is fine to make a minor release to make a single new
feature available; equally, it is fine to make a minor release that includes a number of
changes.

Major releases obviously need to be given careful consideration, and should be done as
infrequently as possible, as they will break existing code and/or affect data
compatibility in some way.

Release procedure
~~~~~~~~~~~~~~~~~

.. note:: 

   Most of the release process is now handled by github workflow which should
   automatically push a release to PyPI if a tag is pushed. 

Checkout and update the master branch::

    $ git checkout master
    $ git pull

Verify all tests pass on all supported Python versions, and docs build::

    $ tox

Tag the version (where "X.X.X" stands for the version number, e.g., "2.2.0")::

    $ version=X.X.X
    $ git tag -a v$version -m v$version
    $ git push origin v$version

Create a GitHub release in order to generate the Zenodo DOI and
review the automatically generated zarr-feedstock PR.
.. _tutorial:

Tutorial
========

Zarr provides classes and functions for working with N-dimensional arrays that
behave like NumPy arrays but whose data is divided into chunks and each chunk is
compressed. If you are already familiar with HDF5 then Zarr arrays provide
similar functionality, but with some additional flexibility.

.. _tutorial_create:

Creating an array
-----------------

Zarr has several functions for creating arrays. For example::

    >>> import zarr
    >>> z = zarr.zeros((10000, 10000), chunks=(1000, 1000), dtype='i4')
    >>> z
    <zarr.core.Array (10000, 10000) int32>

The code above creates a 2-dimensional array of 32-bit integers with 10000 rows
and 10000 columns, divided into chunks where each chunk has 1000 rows and 1000
columns (and so there will be 100 chunks in total).

For a complete list of array creation routines see the :mod:`zarr.creation`
module documentation.

.. _tutorial_array:

Reading and writing data
------------------------

Zarr arrays support a similar interface to NumPy arrays for reading and writing
data. For example, the entire array can be filled with a scalar value::

    >>> z[:] = 42

Regions of the array can also be written to, e.g.::

    >>> import numpy as np
    >>> z[0, :] = np.arange(10000)
    >>> z[:, 0] = np.arange(10000)

The contents of the array can be retrieved by slicing, which will load the
requested region into memory as a NumPy array, e.g.::

    >>> z[0, 0]
    0
    >>> z[-1, -1]
    42
    >>> z[0, :]
    array([   0,    1,    2, ..., 9997, 9998, 9999], dtype=int32)
    >>> z[:, 0]
    array([   0,    1,    2, ..., 9997, 9998, 9999], dtype=int32)
    >>> z[:]
    array([[   0,    1,    2, ..., 9997, 9998, 9999],
           [   1,   42,   42, ...,   42,   42,   42],
           [   2,   42,   42, ...,   42,   42,   42],
           ...,
           [9997,   42,   42, ...,   42,   42,   42],
           [9998,   42,   42, ...,   42,   42,   42],
           [9999,   42,   42, ...,   42,   42,   42]], dtype=int32)

.. _tutorial_persist:

Persistent arrays
-----------------

In the examples above, compressed data for each chunk of the array was stored in
main memory. Zarr arrays can also be stored on a file system, enabling
persistence of data between sessions. For example::

    >>> z1 = zarr.open('data/example.zarr', mode='w', shape=(10000, 10000),
    ...                chunks=(1000, 1000), dtype='i4')

The array above will store its configuration metadata and all compressed chunk
data in a directory called 'data/example.zarr' relative to the current working
directory. The :func:`zarr.convenience.open` function provides a convenient way
to create a new persistent array or continue working with an existing
array. Note that although the function is called "open", there is no need to
close an array: data are automatically flushed to disk, and files are
automatically closed whenever an array is modified.

Persistent arrays support the same interface for reading and writing data,
e.g.::

    >>> z1[:] = 42
    >>> z1[0, :] = np.arange(10000)
    >>> z1[:, 0] = np.arange(10000)

Check that the data have been written and can be read again::

    >>> z2 = zarr.open('data/example.zarr', mode='r')
    >>> np.all(z1[:] == z2[:])
    True

If you are just looking for a fast and convenient way to save NumPy arrays to
disk then load back into memory later, the functions
:func:`zarr.convenience.save` and :func:`zarr.convenience.load` may be
useful. E.g.::

    >>> a = np.arange(10)
    >>> zarr.save('data/example.zarr', a)
    >>> zarr.load('data/example.zarr')
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

Please note that there are a number of other options for persistent array
storage, see the section on :ref:`tutorial_storage` below.

.. _tutorial_resize:

Resizing and appending
----------------------

A Zarr array can be resized, which means that any of its dimensions can be
increased or decreased in length. For example::

    >>> z = zarr.zeros(shape=(10000, 10000), chunks=(1000, 1000))
    >>> z[:] = 42
    >>> z.resize(20000, 10000)
    >>> z.shape
    (20000, 10000)

Note that when an array is resized, the underlying data are not rearranged in
any way. If one or more dimensions are shrunk, any chunks falling outside the
new array shape will be deleted from the underlying store.

For convenience, Zarr arrays also provide an ``append()`` method, which can be
used to append data to any axis. E.g.::

    >>> a = np.arange(10000000, dtype='i4').reshape(10000, 1000)
    >>> z = zarr.array(a, chunks=(1000, 100))
    >>> z.shape
    (10000, 1000)
    >>> z.append(a)
    (20000, 1000)
    >>> z.append(np.vstack([a, a]), axis=1)
    (20000, 2000)
    >>> z.shape
    (20000, 2000)

.. _tutorial_compress:

Compressors
-----------

A number of different compressors can be used with Zarr. A separate package
called NumCodecs_ is available which provides a common interface to various
compressor libraries including Blosc, Zstandard, LZ4, Zlib, BZ2 and
LZMA. Different compressors can be provided via the ``compressor`` keyword
argument accepted by all array creation functions. For example::

    >>> from numcodecs import Blosc
    >>> compressor = Blosc(cname='zstd', clevel=3, shuffle=Blosc.BITSHUFFLE)
    >>> data = np.arange(100000000, dtype='i4').reshape(10000, 10000)
    >>> z = zarr.array(data, chunks=(1000, 1000), compressor=compressor)
    >>> z.compressor
    Blosc(cname='zstd', clevel=3, shuffle=BITSHUFFLE, blocksize=0)

This array above will use Blosc as the primary compressor, using the Zstandard
algorithm (compression level 3) internally within Blosc, and with the
bit-shuffle filter applied.

When using a compressor, it can be useful to get some diagnostics on the
compression ratio. Zarr arrays provide a ``info`` property which can be used to
print some diagnostics, e.g.::

    >>> z.info
    Type               : zarr.core.Array
    Data type          : int32
    Shape              : (10000, 10000)
    Chunk shape        : (1000, 1000)
    Order              : C
    Read-only          : False
    Compressor         : Blosc(cname='zstd', clevel=3, shuffle=BITSHUFFLE,
                       : blocksize=0)
    Store type         : zarr.storage.KVStore
    No. bytes          : 400000000 (381.5M)
    No. bytes stored   : 3379344 (3.2M)
    Storage ratio      : 118.4
    Chunks initialized : 100/100

If you don't specify a compressor, by default Zarr uses the Blosc
compressor. Blosc is generally very fast and can be configured in a variety of
ways to improve the compression ratio for different types of data. Blosc is in
fact a "meta-compressor", which means that it can use a number of different
compression algorithms internally to compress the data. Blosc also provides
highly optimized implementations of byte- and bit-shuffle filters, which can
improve compression ratios for some data. A list of the internal compression
libraries available within Blosc can be obtained via::

    >>> from numcodecs import blosc
    >>> blosc.list_compressors()
    ['blosclz', 'lz4', 'lz4hc', 'snappy', 'zlib', 'zstd']

In addition to Blosc, other compression libraries can also be used. For example,
here is an array using Zstandard compression, level 1::

    >>> from numcodecs import Zstd
    >>> z = zarr.array(np.arange(100000000, dtype='i4').reshape(10000, 10000),
    ...                chunks=(1000, 1000), compressor=Zstd(level=1))
    >>> z.compressor
    Zstd(level=1)

Here is an example using LZMA with a custom filter pipeline including LZMA's
built-in delta filter::

    >>> import lzma
    >>> lzma_filters = [dict(id=lzma.FILTER_DELTA, dist=4),
    ...                 dict(id=lzma.FILTER_LZMA2, preset=1)]
    >>> from numcodecs import LZMA
    >>> compressor = LZMA(filters=lzma_filters)
    >>> z = zarr.array(np.arange(100000000, dtype='i4').reshape(10000, 10000),
    ...                chunks=(1000, 1000), compressor=compressor)
    >>> z.compressor
    LZMA(format=1, check=-1, preset=None, filters=[{'dist': 4, 'id': 3}, {'id': 33, 'preset': 1}])

The default compressor can be changed by setting the value of the
``zarr.storage.default_compressor`` variable, e.g.::

    >>> import zarr.storage
    >>> from numcodecs import Zstd, Blosc
    >>> # switch to using Zstandard
    ... zarr.storage.default_compressor = Zstd(level=1)
    >>> z = zarr.zeros(100000000, chunks=1000000)
    >>> z.compressor
    Zstd(level=1)
    >>> # switch back to Blosc defaults
    ... zarr.storage.default_compressor = Blosc()

To disable compression, set ``compressor=None`` when creating an array, e.g.::

    >>> z = zarr.zeros(100000000, chunks=1000000, compressor=None)
    >>> z.compressor is None
    True

.. _tutorial_filters:

Filters
-------

In some cases, compression can be improved by transforming the data in some
way. For example, if nearby values tend to be correlated, then shuffling the
bytes within each numerical value or storing the difference between adjacent
values may increase compression ratio. Some compressors provide built-in filters
that apply transformations to the data prior to compression. For example, the
Blosc compressor has built-in implementations of byte- and bit-shuffle filters,
and the LZMA compressor has a built-in implementation of a delta
filter. However, to provide additional flexibility for implementing and using
filters in combination with different compressors, Zarr also provides a
mechanism for configuring filters outside of the primary compressor.

Here is an example using a delta filter with the Blosc compressor::

    >>> from numcodecs import Blosc, Delta
    >>> filters = [Delta(dtype='i4')]
    >>> compressor = Blosc(cname='zstd', clevel=1, shuffle=Blosc.SHUFFLE)
    >>> data = np.arange(100000000, dtype='i4').reshape(10000, 10000)
    >>> z = zarr.array(data, chunks=(1000, 1000), filters=filters, compressor=compressor)
    >>> z.info
    Type               : zarr.core.Array
    Data type          : int32
    Shape              : (10000, 10000)
    Chunk shape        : (1000, 1000)
    Order              : C
    Read-only          : False
    Filter [0]         : Delta(dtype='<i4')
    Compressor         : Blosc(cname='zstd', clevel=1, shuffle=SHUFFLE, blocksize=0)
    Store type         : zarr.storage.KVStore
    No. bytes          : 400000000 (381.5M)
    No. bytes stored   : 1290562 (1.2M)
    Storage ratio      : 309.9
    Chunks initialized : 100/100

For more information about available filter codecs, see the `Numcodecs
<http://numcodecs.readthedocs.io/>`_ documentation.

.. _tutorial_groups:

Groups
------

Zarr supports hierarchical organization of arrays via groups. As with arrays,
groups can be stored in memory, on disk, or via other storage systems that
support a similar interface.

To create a group, use the :func:`zarr.group` function::

    >>> root = zarr.group()
    >>> root
    <zarr.hierarchy.Group '/'>

Groups have a similar API to the Group class from `h5py
<http://www.h5py.org/>`_.  For example, groups can contain other groups::

    >>> foo = root.create_group('foo')
    >>> bar = foo.create_group('bar')

Groups can also contain arrays, e.g.::

    >>> z1 = bar.zeros('baz', shape=(10000, 10000), chunks=(1000, 1000), dtype='i4')
    >>> z1
    <zarr.core.Array '/foo/bar/baz' (10000, 10000) int32>

Arrays are known as "datasets" in HDF5 terminology. For compatibility with h5py,
Zarr groups also implement the ``create_dataset()`` and ``require_dataset()``
methods, e.g.::

    >>> z = bar.create_dataset('quux', shape=(10000, 10000), chunks=(1000, 1000), dtype='i4')
    >>> z
    <zarr.core.Array '/foo/bar/quux' (10000, 10000) int32>

Members of a group can be accessed via the suffix notation, e.g.::

    >>> root['foo']
    <zarr.hierarchy.Group '/foo'>

The '/' character can be used to access multiple levels of the hierarchy in one
call, e.g.::

    >>> root['foo/bar']
    <zarr.hierarchy.Group '/foo/bar'>
    >>> root['foo/bar/baz']
    <zarr.core.Array '/foo/bar/baz' (10000, 10000) int32>

The :func:`zarr.hierarchy.Group.tree` method can be used to print a tree
representation of the hierarchy, e.g.::

    >>> root.tree()
    /
     └── foo
         └── bar
             ├── baz (10000, 10000) int32
             └── quux (10000, 10000) int32

The :func:`zarr.convenience.open` function provides a convenient way to create or
re-open a group stored in a directory on the file-system, with sub-groups stored in
sub-directories, e.g.::

    >>> root = zarr.open('data/group.zarr', mode='w')
    >>> root
    <zarr.hierarchy.Group '/'>
    >>> z = root.zeros('foo/bar/baz', shape=(10000, 10000), chunks=(1000, 1000), dtype='i4')
    >>> z
    <zarr.core.Array '/foo/bar/baz' (10000, 10000) int32>

Groups can be used as context managers (in a ``with`` statement).
If the underlying store has a ``close`` method, it will be called on exit.

For more information on groups see the :mod:`zarr.hierarchy` and
:mod:`zarr.convenience` API docs.

.. _tutorial_diagnostics:

Array and group diagnostics
---------------------------

Diagnostic information about arrays and groups is available via the ``info``
property. E.g.::

    >>> root = zarr.group()
    >>> foo = root.create_group('foo')
    >>> bar = foo.zeros('bar', shape=1000000, chunks=100000, dtype='i8')
    >>> bar[:] = 42
    >>> baz = foo.zeros('baz', shape=(1000, 1000), chunks=(100, 100), dtype='f4')
    >>> baz[:] = 4.2
    >>> root.info
    Name        : /
    Type        : zarr.hierarchy.Group
    Read-only   : False
    Store type  : zarr.storage.MemoryStore
    No. members : 1
    No. arrays  : 0
    No. groups  : 1
    Groups      : foo

    >>> foo.info
    Name        : /foo
    Type        : zarr.hierarchy.Group
    Read-only   : False
    Store type  : zarr.storage.MemoryStore
    No. members : 2
    No. arrays  : 2
    No. groups  : 0
    Arrays      : bar, baz

    >>> bar.info
    Name               : /foo/bar
    Type               : zarr.core.Array
    Data type          : int64
    Shape              : (1000000,)
    Chunk shape        : (100000,)
    Order              : C
    Read-only          : False
    Compressor         : Blosc(cname='lz4', clevel=5, shuffle=SHUFFLE, blocksize=0)
    Store type         : zarr.storage.MemoryStore
    No. bytes          : 8000000 (7.6M)
    No. bytes stored   : 33240 (32.5K)
    Storage ratio      : 240.7
    Chunks initialized : 10/10

    >>> baz.info
    Name               : /foo/baz
    Type               : zarr.core.Array
    Data type          : float32
    Shape              : (1000, 1000)
    Chunk shape        : (100, 100)
    Order              : C
    Read-only          : False
    Compressor         : Blosc(cname='lz4', clevel=5, shuffle=SHUFFLE, blocksize=0)
    Store type         : zarr.storage.MemoryStore
    No. bytes          : 4000000 (3.8M)
    No. bytes stored   : 23943 (23.4K)
    Storage ratio      : 167.1
    Chunks initialized : 100/100

Groups also have the :func:`zarr.hierarchy.Group.tree` method, e.g.::

    >>> root.tree()
    /
     └── foo
         ├── bar (1000000,) int64
         └── baz (1000, 1000) float32

If you're using Zarr within a Jupyter notebook (requires
`ipytree <https://github.com/QuantStack/ipytree>`_), calling ``tree()`` will generate an
interactive tree representation, see the `repr_tree.ipynb notebook
<http://nbviewer.jupyter.org/github/zarr-developers/zarr-python/blob/master/notebooks/repr_tree.ipynb>`_
for more examples.

.. _tutorial_attrs:

User attributes
---------------

Zarr arrays and groups support custom key/value attributes, which can be useful for
storing application-specific metadata. For example::

    >>> root = zarr.group()
    >>> root.attrs['foo'] = 'bar'
    >>> z = root.zeros('zzz', shape=(10000, 10000))
    >>> z.attrs['baz'] = 42
    >>> z.attrs['qux'] = [1, 4, 7, 12]
    >>> sorted(root.attrs)
    ['foo']
    >>> 'foo' in root.attrs
    True
    >>> root.attrs['foo']
    'bar'
    >>> sorted(z.attrs)
    ['baz', 'qux']
    >>> z.attrs['baz']
    42
    >>> z.attrs['qux']
    [1, 4, 7, 12]

Internally Zarr uses JSON to store array attributes, so attribute values must be
JSON serializable.

.. _tutorial_indexing:

Advanced indexing
-----------------

As of version 2.2, Zarr arrays support several methods for advanced or "fancy"
indexing, which enable a subset of data items to be extracted or updated in an
array without loading the entire array into memory.

Note that although this functionality is similar to some of the advanced
indexing capabilities available on NumPy arrays and on h5py datasets, **the Zarr
API for advanced indexing is different from both NumPy and h5py**, so please
read this section carefully.  For a complete description of the indexing API,
see the documentation for the :class:`zarr.core.Array` class.

Indexing with coordinate arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Items from a Zarr array can be extracted by providing an integer array of
coordinates. E.g.::

    >>> z = zarr.array(np.arange(10))
    >>> z[:]
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    >>> z.get_coordinate_selection([1, 4])
    array([1, 4])

Coordinate arrays can also be used to update data, e.g.::

    >>> z.set_coordinate_selection([1, 4], [-1, -2])
    >>> z[:]
    array([ 0, -1,  2,  3, -2,  5,  6,  7,  8,  9])

For multidimensional arrays, coordinates must be provided for each dimension,
e.g.::

    >>> z = zarr.array(np.arange(15).reshape(3, 5))
    >>> z[:]
    array([[ 0,  1,  2,  3,  4],
           [ 5,  6,  7,  8,  9],
           [10, 11, 12, 13, 14]])
    >>> z.get_coordinate_selection(([0, 2], [1, 3]))
    array([ 1, 13])
    >>> z.set_coordinate_selection(([0, 2], [1, 3]), [-1, -2])
    >>> z[:]
    array([[ 0, -1,  2,  3,  4],
           [ 5,  6,  7,  8,  9],
           [10, 11, 12, -2, 14]])

For convenience, coordinate indexing is also available via the ``vindex``
property, as well as the square bracket operator, e.g.::

    >>> z.vindex[[0, 2], [1, 3]]
    array([-1, -2])
    >>> z.vindex[[0, 2], [1, 3]] = [-3, -4]
    >>> z[:]
    array([[ 0, -3,  2,  3,  4],
           [ 5,  6,  7,  8,  9],
           [10, 11, 12, -4, 14]])
    >>> z[[0, 2], [1, 3]]
    array([-3, -4])

When the indexing arrays have different shapes, they are broadcast together.
That is, the following two calls are equivalent::

    >>> z[1, [1, 3]]
    array([5, 7])
    >>> z[[1, 1], [1, 3]]
    array([5, 7])

Indexing with a mask array
~~~~~~~~~~~~~~~~~~~~~~~~~~

Items can also be extracted by providing a Boolean mask. E.g.::

    >>> z = zarr.array(np.arange(10))
    >>> z[:]
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    >>> sel = np.zeros_like(z, dtype=bool)
    >>> sel[1] = True
    >>> sel[4] = True
    >>> z.get_mask_selection(sel)
    array([1, 4])
    >>> z.set_mask_selection(sel, [-1, -2])
    >>> z[:]
    array([ 0, -1,  2,  3, -2,  5,  6,  7,  8,  9])

Here's a multidimensional example::

    >>> z = zarr.array(np.arange(15).reshape(3, 5))
    >>> z[:]
    array([[ 0,  1,  2,  3,  4],
           [ 5,  6,  7,  8,  9],
           [10, 11, 12, 13, 14]])
    >>> sel = np.zeros_like(z, dtype=bool)
    >>> sel[0, 1] = True
    >>> sel[2, 3] = True
    >>> z.get_mask_selection(sel)
    array([ 1, 13])
    >>> z.set_mask_selection(sel, [-1, -2])
    >>> z[:]
    array([[ 0, -1,  2,  3,  4],
           [ 5,  6,  7,  8,  9],
           [10, 11, 12, -2, 14]])

For convenience, mask indexing is also available via the ``vindex`` property,
e.g.::

    >>> z.vindex[sel]
    array([-1, -2])
    >>> z.vindex[sel] = [-3, -4]
    >>> z[:]
    array([[ 0, -3,  2,  3,  4],
           [ 5,  6,  7,  8,  9],
           [10, 11, 12, -4, 14]])

Mask indexing is conceptually the same as coordinate indexing, and is
implemented internally via the same machinery. Both styles of indexing allow
selecting arbitrary items from an array, also known as point selection.

Orthogonal indexing
~~~~~~~~~~~~~~~~~~~

Zarr arrays also support methods for orthogonal indexing, which allows
selections to be made along each dimension of an array independently. For
example, this allows selecting a subset of rows and/or columns from a
2-dimensional array. E.g.::

    >>> z = zarr.array(np.arange(15).reshape(3, 5))
    >>> z[:]
    array([[ 0,  1,  2,  3,  4],
           [ 5,  6,  7,  8,  9],
           [10, 11, 12, 13, 14]])
    >>> z.get_orthogonal_selection(([0, 2], slice(None)))  # select first and third rows
    array([[ 0,  1,  2,  3,  4],
           [10, 11, 12, 13, 14]])
    >>> z.get_orthogonal_selection((slice(None), [1, 3]))  # select second and fourth columns
    array([[ 1,  3],
           [ 6,  8],
           [11, 13]])
    >>> z.get_orthogonal_selection(([0, 2], [1, 3]))       # select rows [0, 2] and columns [1, 4]
    array([[ 1,  3],
           [11, 13]])

Data can also be modified, e.g.::

    >>> z.set_orthogonal_selection(([0, 2], [1, 3]), [[-1, -2], [-3, -4]])
    >>> z[:]
    array([[ 0, -1,  2, -2,  4],
           [ 5,  6,  7,  8,  9],
           [10, -3, 12, -4, 14]])

For convenience, the orthogonal indexing functionality is also available via the
``oindex`` property, e.g.::

    >>> z = zarr.array(np.arange(15).reshape(3, 5))
    >>> z.oindex[[0, 2], :]  # select first and third rows
    array([[ 0,  1,  2,  3,  4],
           [10, 11, 12, 13, 14]])
    >>> z.oindex[:, [1, 3]]  # select second and fourth columns
    array([[ 1,  3],
           [ 6,  8],
           [11, 13]])
    >>> z.oindex[[0, 2], [1, 3]]  # select rows [0, 2] and columns [1, 4]
    array([[ 1,  3],
           [11, 13]])
    >>> z.oindex[[0, 2], [1, 3]] = [[-1, -2], [-3, -4]]
    >>> z[:]
    array([[ 0, -1,  2, -2,  4],
           [ 5,  6,  7,  8,  9],
           [10, -3, 12, -4, 14]])

Any combination of integer, slice, 1D integer array and/or 1D Boolean array can
be used for orthogonal indexing.

Indexing fields in structured arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All selection methods support a ``fields`` parameter which allows retrieving or
replacing data for a specific field in an array with a structured dtype. E.g.::

    >>> a = np.array([(b'aaa', 1, 4.2),
    ...               (b'bbb', 2, 8.4),
    ...               (b'ccc', 3, 12.6)],
    ...              dtype=[('foo', 'S3'), ('bar', 'i4'), ('baz', 'f8')])
    >>> z = zarr.array(a)
    >>> z['foo']
    array([b'aaa', b'bbb', b'ccc'],
          dtype='|S3')
    >>> z['baz']
    array([  4.2,   8.4,  12.6])
    >>> z.get_basic_selection(slice(0, 2), fields='bar')
    array([1, 2], dtype=int32)
    >>> z.get_coordinate_selection([0, 2], fields=['foo', 'baz'])
    array([(b'aaa',   4.2), (b'ccc',  12.6)],
          dtype=[('foo', 'S3'), ('baz', '<f8')])

.. _tutorial_storage:

Storage alternatives
--------------------

Zarr can use any object that implements the ``MutableMapping`` interface from
the :mod:`collections` module in the Python standard library as the store for a
group or an array.

Some pre-defined storage classes are provided in the :mod:`zarr.storage`
module. For example, the :class:`zarr.storage.DirectoryStore` class provides a
``MutableMapping`` interface to a directory on the local file system. This is
used under the hood by the :func:`zarr.convenience.open` function. In other words,
the following code::

    >>> z = zarr.open('data/example.zarr', mode='w', shape=1000000, dtype='i4')

...is short-hand for::

    >>> store = zarr.DirectoryStore('data/example.zarr')
    >>> z = zarr.create(store=store, overwrite=True, shape=1000000, dtype='i4')

...and the following code::

    >>> root = zarr.open('data/example.zarr', mode='w')

...is short-hand for::

    >>> store = zarr.DirectoryStore('data/example.zarr')
    >>> root = zarr.group(store=store, overwrite=True)

Any other compatible storage class could be used in place of
:class:`zarr.storage.DirectoryStore` in the code examples above. For example,
here is an array stored directly into a Zip file, via the
:class:`zarr.storage.ZipStore` class::

    >>> store = zarr.ZipStore('data/example.zip', mode='w')
    >>> root = zarr.group(store=store)
    >>> z = root.zeros('foo/bar', shape=(1000, 1000), chunks=(100, 100), dtype='i4')
    >>> z[:] = 42
    >>> store.close()

Re-open and check that data have been written::

    >>> store = zarr.ZipStore('data/example.zip', mode='r')
    >>> root = zarr.group(store=store)
    >>> z = root['foo/bar']
    >>> z[:]
    array([[42, 42, 42, ..., 42, 42, 42],
           [42, 42, 42, ..., 42, 42, 42],
           [42, 42, 42, ..., 42, 42, 42],
           ...,
           [42, 42, 42, ..., 42, 42, 42],
           [42, 42, 42, ..., 42, 42, 42],
           [42, 42, 42, ..., 42, 42, 42]], dtype=int32)
    >>> store.close()

Note that there are some limitations on how Zip files can be used, because items
within a Zip file cannot be updated in place. This means that data in the array
should only be written once and write operations should be aligned with chunk
boundaries. Note also that the ``close()`` method must be called after writing
any data to the store, otherwise essential records will not be written to the
underlying zip file.

Another storage alternative is the :class:`zarr.storage.DBMStore` class, added
in Zarr version 2.2. This class allows any DBM-style database to be used for
storing an array or group. Here is an example using a Berkeley DB B-tree
database for storage (requires `bsddb3
<https://www.jcea.es/programacion/pybsddb.htm>`_ to be installed)::

    >>> import bsddb3
    >>> store = zarr.DBMStore('data/example.bdb', open=bsddb3.btopen)
    >>> root = zarr.group(store=store, overwrite=True)
    >>> z = root.zeros('foo/bar', shape=(1000, 1000), chunks=(100, 100), dtype='i4')
    >>> z[:] = 42
    >>> store.close()

Also added in Zarr version 2.2 is the :class:`zarr.storage.LMDBStore` class which
enables the lightning memory-mapped database (LMDB) to be used for storing an array or
group (requires `lmdb <http://lmdb.readthedocs.io/>`_ to be installed)::

    >>> store = zarr.LMDBStore('data/example.lmdb')
    >>> root = zarr.group(store=store, overwrite=True)
    >>> z = root.zeros('foo/bar', shape=(1000, 1000), chunks=(100, 100), dtype='i4')
    >>> z[:] = 42
    >>> store.close()

In Zarr version 2.3 is the :class:`zarr.storage.SQLiteStore` class which
enables the SQLite database to be used for storing an array or group (requires
Python is built with SQLite support)::

    >>> store = zarr.SQLiteStore('data/example.sqldb')
    >>> root = zarr.group(store=store, overwrite=True)
    >>> z = root.zeros('foo/bar', shape=(1000, 1000), chunks=(100, 100), dtype='i4')
    >>> z[:] = 42
    >>> store.close()

Also added in Zarr version 2.3 are two storage classes for interfacing with server-client
databases. The :class:`zarr.storage.RedisStore` class interfaces `Redis <https://redis.io/>`_
(an in memory data structure store), and the :class:`zarr.storage.MongoDB` class interfaces
with `MongoDB <https://www.mongodb.com/>`_ (an object oriented NoSQL database). These stores
respectively require the `redis-py <https://redis-py.readthedocs.io>`_ and
`pymongo <https://api.mongodb.com/python/current/>`_ packages to be installed. 

For compatibility with the `N5 <https://github.com/saalfeldlab/n5>`_ data format, Zarr also provides
an N5 backend (this is currently an experimental feature). Similar to the zip storage class, an
:class:`zarr.n5.N5Store` can be instantiated directly::

    >>> store = zarr.N5Store('data/example.n5')
    >>> root = zarr.group(store=store)
    >>> z = root.zeros('foo/bar', shape=(1000, 1000), chunks=(100, 100), dtype='i4')
    >>> z[:] = 42

For convenience, the N5 backend will automatically be chosen when the filename
ends with `.n5`::

    >>> root = zarr.open('data/example.n5', mode='w')

Distributed/cloud storage
~~~~~~~~~~~~~~~~~~~~~~~~~

It is also possible to use distributed storage systems. The Dask project has
implementations of the ``MutableMapping`` interface for Amazon S3 (`S3Map
<http://s3fs.readthedocs.io/en/latest/api.html#s3fs.mapping.S3Map>`_), Hadoop
Distributed File System (`HDFSMap
<http://hdfs3.readthedocs.io/en/latest/api.html#hdfs3.mapping.HDFSMap>`_) and
Google Cloud Storage (`GCSMap
<http://gcsfs.readthedocs.io/en/latest/api.html#gcsfs.mapping.GCSMap>`_), which
can be used with Zarr.

Here is an example using S3Map to read an array created previously::

    >>> import s3fs
    >>> import zarr
    >>> s3 = s3fs.S3FileSystem(anon=True, client_kwargs=dict(region_name='eu-west-2'))
    >>> store = s3fs.S3Map(root='zarr-demo/store', s3=s3, check=False)
    >>> root = zarr.group(store=store)
    >>> z = root['foo/bar/baz']
    >>> z
    <zarr.core.Array '/foo/bar/baz' (21,) |S1>
    >>> z.info
    Name               : /foo/bar/baz
    Type               : zarr.core.Array
    Data type          : |S1
    Shape              : (21,)
    Chunk shape        : (7,)
    Order              : C
    Read-only          : False
    Compressor         : Blosc(cname='lz4', clevel=5, shuffle=SHUFFLE, blocksize=0)
    Store type         : zarr.storage.KVStore
    No. bytes          : 21
    No. bytes stored   : 382
    Storage ratio      : 0.1
    Chunks initialized : 3/3
    >>> z[:]
    array([b'H', b'e', b'l', b'l', b'o', b' ', b'f', b'r', b'o', b'm', b' ',
           b't', b'h', b'e', b' ', b'c', b'l', b'o', b'u', b'd', b'!'],
          dtype='|S1')
    >>> z[:].tobytes()
    b'Hello from the cloud!'

Zarr now also has a builtin storage backend for Azure Blob Storage.
The class is :class:`zarr.storage.ABSStore` (requires
`azure-storage-blob <https://docs.microsoft.com/en-us/azure/storage/blobs/storage-quickstart-blobs-python>`_
to be installed)::

    >>> import azure.storage.blob
    >>> container_client = azure.storage.blob.ContainerClient(...)  # doctest: +SKIP
    >>> store = zarr.ABSStore(client=container_client, prefix='zarr-testing')  # doctest: +SKIP
    >>> root = zarr.group(store=store, overwrite=True)  # doctest: +SKIP
    >>> z = root.zeros('foo/bar', shape=(1000, 1000), chunks=(100, 100), dtype='i4')  # doctest: +SKIP
    >>> z[:] = 42  # doctest: +SKIP

When using an actual storage account, provide ``account_name`` and
``account_key`` arguments to :class:`zarr.storage.ABSStore`, the
above client is just testing against the emulator. Please also note
that this is an experimental feature.

Note that retrieving data from a remote service via the network can be significantly
slower than retrieving data from a local file system, and will depend on network latency
and bandwidth between the client and server systems. If you are experiencing poor
performance, there are several things you can try. One option is to increase the array
chunk size, which will reduce the number of chunks and thus reduce the number of network
round-trips required to retrieve data for an array (and thus reduce the impact of network
latency). Another option is to try to increase the compression ratio by changing
compression options or trying a different compressor (which will reduce the impact of
limited network bandwidth).

As of version 2.2, Zarr also provides the :class:`zarr.storage.LRUStoreCache`
which can be used to implement a local in-memory cache layer over a remote
store. E.g.::

    >>> s3 = s3fs.S3FileSystem(anon=True, client_kwargs=dict(region_name='eu-west-2'))
    >>> store = s3fs.S3Map(root='zarr-demo/store', s3=s3, check=False)
    >>> cache = zarr.LRUStoreCache(store, max_size=2**28)
    >>> root = zarr.group(store=cache)
    >>> z = root['foo/bar/baz']
    >>> from timeit import timeit
    >>> # first data access is relatively slow, retrieved from store
    ... timeit('print(z[:].tobytes())', number=1, globals=globals())  # doctest: +SKIP
    b'Hello from the cloud!'
    0.1081731989979744
    >>> # second data access is faster, uses cache
    ... timeit('print(z[:].tobytes())', number=1, globals=globals())  # doctest: +SKIP
    b'Hello from the cloud!'
    0.0009490990014455747

If you are still experiencing poor performance with distributed/cloud storage,
please raise an issue on the GitHub issue tracker with any profiling data you
can provide, as there may be opportunities to optimise further either within
Zarr or within the mapping interface to the storage.

IO with ``fsspec``
~~~~~~~~~~~~~~~~~~

As of version 2.5, zarr supports passing URLs directly to `fsspec`_,
and having it create the "mapping" instance automatically. This means, that
for all of the backend storage implementations `supported by fsspec`_,
you can skip importing and configuring the storage explicitly.
For example::

    >>> g = zarr.open_group("s3://zarr-demo/store", storage_options={'anon': True})   # doctest: +SKIP
    >>> g['foo/bar/baz'][:].tobytes()  # doctest: +SKIP
    b'Hello from the cloud!'

The provision of the protocol specifier "s3://" will select the correct backend.
Notice the kwargs ``storage_options``, used to pass parameters to that backend.

As of version 2.6, write mode and complex URLs are also supported, such as::

    >>> g = zarr.open_group("simplecache::s3://zarr-demo/store",
    ...                     storage_options={"s3": {'anon': True}})  # doctest: +SKIP
    >>> g['foo/bar/baz'][:].tobytes()  # downloads target file  # doctest: +SKIP
    b'Hello from the cloud!'
    >>> g['foo/bar/baz'][:].tobytes()  # uses cached file  # doctest: +SKIP
    b'Hello from the cloud!'

The second invocation here will be much faster. Note that the ``storage_options``
have become more complex here, to account for the two parts of the supplied
URL.

.. _fsspec: https://filesystem-spec.readthedocs.io/en/latest/

.. _supported by fsspec: https://filesystem-spec.readthedocs.io/en/latest/api.html#built-in-implementations

.. _tutorial_copy:

Consolidating metadata
~~~~~~~~~~~~~~~~~~~~~~

Since there is a significant overhead for every connection to a cloud object
store such as S3, the pattern described in the previous section may incur
significant latency while scanning the metadata of the array hierarchy, even
though each individual metadata object is small.  For cases such as these, once
the data are static and can be regarded as read-only, at least for the
metadata/structure of the array hierarchy, the many metadata objects can be
consolidated into a single one via
:func:`zarr.convenience.consolidate_metadata`. Doing this can greatly increase
the speed of reading the array metadata, e.g.::

   >>> zarr.consolidate_metadata(store)  # doctest: +SKIP

This creates a special key with a copy of all of the metadata from all of the
metadata objects in the store.

Later, to open a Zarr store with consolidated metadata, use
:func:`zarr.convenience.open_consolidated`, e.g.::

   >>> root = zarr.open_consolidated(store)  # doctest: +SKIP

This uses the special key to read all of the metadata in a single call to the
backend storage.

Note that, the hierarchy could still be opened in the normal way and altered,
causing the consolidated metadata to become out of sync with the real state of
the array hierarchy. In this case,
:func:`zarr.convenience.consolidate_metadata` would need to be called again.

To protect against consolidated metadata accidentally getting out of sync, the
root group returned by :func:`zarr.convenience.open_consolidated` is read-only
for the metadata, meaning that no new groups or arrays can be created, and
arrays cannot be resized. However, data values with arrays can still be updated.

Copying/migrating data
----------------------

If you have some data in an HDF5 file and would like to copy some or all of it
into a Zarr group, or vice-versa, the :func:`zarr.convenience.copy` and
:func:`zarr.convenience.copy_all` functions can be used. Here's an example
copying a group named 'foo' from an HDF5 file to a Zarr group::

    >>> import h5py
    >>> import zarr
    >>> import numpy as np
    >>> source = h5py.File('data/example.h5', mode='w')
    >>> foo = source.create_group('foo')
    >>> baz = foo.create_dataset('bar/baz', data=np.arange(100), chunks=(50,))
    >>> spam = source.create_dataset('spam', data=np.arange(100, 200), chunks=(30,))
    >>> zarr.tree(source)
    /
     ├── foo
     │   └── bar
     │       └── baz (100,) int64
     └── spam (100,) int64
    >>> dest = zarr.open_group('data/example.zarr', mode='w')
    >>> from sys import stdout
    >>> zarr.copy(source['foo'], dest, log=stdout)
    copy /foo
    copy /foo/bar
    copy /foo/bar/baz (100,) int64
    all done: 3 copied, 0 skipped, 800 bytes copied
    (3, 0, 800)
    >>> dest.tree()  # N.B., no spam
    /
     └── foo
         └── bar
             └── baz (100,) int64
    >>> source.close()

If rather than copying a single group or array you would like to copy all
groups and arrays, use :func:`zarr.convenience.copy_all`, e.g.::

    >>> source = h5py.File('data/example.h5', mode='r')
    >>> dest = zarr.open_group('data/example2.zarr', mode='w')
    >>> zarr.copy_all(source, dest, log=stdout)
    copy /foo
    copy /foo/bar
    copy /foo/bar/baz (100,) int64
    copy /spam (100,) int64
    all done: 4 copied, 0 skipped, 1,600 bytes copied
    (4, 0, 1600)
    >>> dest.tree()
    /
     ├── foo
     │   └── bar
     │       └── baz (100,) int64
     └── spam (100,) int64

If you need to copy data between two Zarr groups, the
:func:`zarr.convenience.copy` and :func:`zarr.convenience.copy_all` functions can
be used and provide the most flexibility. However, if you want to copy data
in the most efficient way possible, without changing any configuration options,
the :func:`zarr.convenience.copy_store` function can be used. This function
copies data directly between the underlying stores, without any decompression or
re-compression, and so should be faster. E.g.::

    >>> import zarr
    >>> import numpy as np
    >>> store1 = zarr.DirectoryStore('data/example.zarr')
    >>> root = zarr.group(store1, overwrite=True)
    >>> baz = root.create_dataset('foo/bar/baz', data=np.arange(100), chunks=(50,))
    >>> spam = root.create_dataset('spam', data=np.arange(100, 200), chunks=(30,))
    >>> root.tree()
    /
     ├── foo
     │   └── bar
     │       └── baz (100,) int64
     └── spam (100,) int64
    >>> from sys import stdout
    >>> store2 = zarr.ZipStore('data/example.zip', mode='w')
    >>> zarr.copy_store(store1, store2, log=stdout)
    copy .zgroup
    copy foo/.zgroup
    copy foo/bar/.zgroup
    copy foo/bar/baz/.zarray
    copy foo/bar/baz/0
    copy foo/bar/baz/1
    copy spam/.zarray
    copy spam/0
    copy spam/1
    copy spam/2
    copy spam/3
    all done: 11 copied, 0 skipped, 1,138 bytes copied
    (11, 0, 1138)
    >>> new_root = zarr.group(store2)
    >>> new_root.tree()
    /
     ├── foo
     │   └── bar
     │       └── baz (100,) int64
     └── spam (100,) int64
    >>> new_root['foo/bar/baz'][:]
    array([ 0,  1,  2,  ..., 97, 98, 99])
    >>> store2.close()  # zip stores need to be closed

.. _tutorial_strings:

String arrays
-------------

There are several options for storing arrays of strings.

If your strings are all ASCII strings, and you know the maximum length of the string in
your array, then you can use an array with a fixed-length bytes dtype. E.g.::

    >>> z = zarr.zeros(10, dtype='S6')
    >>> z
    <zarr.core.Array (10,) |S6>
    >>> z[0] = b'Hello'
    >>> z[1] = b'world!'
    >>> z[:]
    array([b'Hello', b'world!', b'', b'', b'', b'', b'', b'', b'', b''],
          dtype='|S6')

A fixed-length unicode dtype is also available, e.g.::

    >>> greetings = ['¡Hola mundo!', 'Hej Världen!', 'Servus Woid!', 'Hei maailma!',
    ...              'Xin chào thế giới', 'Njatjeta Botë!', 'Γεια σου κόσμε!',
    ...              'こんにちは世界', '世界，你好！', 'Helló, világ!', 'Zdravo svete!',
    ...              'เฮลโลเวิลด์']
    >>> text_data = greetings * 10000
    >>> z = zarr.array(text_data, dtype='U20')
    >>> z
    <zarr.core.Array (120000,) <U20>
    >>> z[:]
    array(['¡Hola mundo!', 'Hej Världen!', 'Servus Woid!', ...,
           'Helló, világ!', 'Zdravo svete!', 'เฮลโลเวิลด์'],
          dtype='<U20')

For variable-length strings, the ``object`` dtype can be used, but a codec must be
provided to encode the data (see also :ref:`tutorial_objects` below). At the time of
writing there are four codecs available that can encode variable length string
objects: :class:`numcodecs.VLenUTF8`, :class:`numcodecs.JSON`, :class:`numcodecs.MsgPack`.
and :class:`numcodecs.Pickle`. E.g. using ``VLenUTF8``::

    >>> import numcodecs
    >>> z = zarr.array(text_data, dtype=object, object_codec=numcodecs.VLenUTF8())
    >>> z
    <zarr.core.Array (120000,) object>
    >>> z.filters
    [VLenUTF8()]
    >>> z[:]
    array(['¡Hola mundo!', 'Hej Världen!', 'Servus Woid!', ...,
           'Helló, világ!', 'Zdravo svete!', 'เฮลโลเวิลด์'], dtype=object)

As a convenience, ``dtype=str`` (or ``dtype=unicode`` on Python 2.7) can be used, which
is a short-hand for ``dtype=object, object_codec=numcodecs.VLenUTF8()``, e.g.::

    >>> z = zarr.array(text_data, dtype=str)
    >>> z
    <zarr.core.Array (120000,) object>
    >>> z.filters
    [VLenUTF8()]
    >>> z[:]
    array(['¡Hola mundo!', 'Hej Världen!', 'Servus Woid!', ...,
           'Helló, világ!', 'Zdravo svete!', 'เฮลโลเวิลด์'], dtype=object)

Variable-length byte strings are also supported via ``dtype=object``. Again an
``object_codec`` is required, which can be one of :class:`numcodecs.VLenBytes` or
:class:`numcodecs.Pickle`. For convenience, ``dtype=bytes`` (or ``dtype=str`` on Python
2.7) can be used as a short-hand for ``dtype=object, object_codec=numcodecs.VLenBytes()``,
e.g.::

    >>> bytes_data = [g.encode('utf-8') for g in greetings] * 10000
    >>> z = zarr.array(bytes_data, dtype=bytes)
    >>> z
    <zarr.core.Array (120000,) object>
    >>> z.filters
    [VLenBytes()]
    >>> z[:]
    array([b'\xc2\xa1Hola mundo!', b'Hej V\xc3\xa4rlden!', b'Servus Woid!',
           ..., b'Hell\xc3\xb3, vil\xc3\xa1g!', b'Zdravo svete!',
           b'\xe0\xb9\x80\xe0\xb8\xae\xe0\xb8\xa5\xe0\xb9\x82\xe0\xb8\xa5\xe0\xb9\x80\xe0\xb8\xa7\xe0\xb8\xb4\xe0\xb8\xa5\xe0\xb8\x94\xe0\xb9\x8c'], dtype=object)

If you know ahead of time all the possible string values that can occur, you could
also use the :class:`numcodecs.Categorize` codec to encode each unique string value as an
integer. E.g.::

    >>> categorize = numcodecs.Categorize(greetings, dtype=object)
    >>> z = zarr.array(text_data, dtype=object, object_codec=categorize)
    >>> z
    <zarr.core.Array (120000,) object>
    >>> z.filters
    [Categorize(dtype='|O', astype='|u1', labels=['¡Hola mundo!', 'Hej Världen!', 'Servus Woid!', ...])]
    >>> z[:]
    array(['¡Hola mundo!', 'Hej Världen!', 'Servus Woid!', ...,
           'Helló, világ!', 'Zdravo svete!', 'เฮลโลเวิลด์'], dtype=object)


.. _tutorial_objects:

Object arrays
-------------

Zarr supports arrays with an "object" dtype. This allows arrays to contain any type of
object, such as variable length unicode strings, or variable length arrays of numbers, or
other possibilities. When creating an object array, a codec must be provided via the
``object_codec`` argument. This codec handles encoding (serialization) of Python objects.
The best codec to use will depend on what type of objects are present in the array.

At the time of writing there are three codecs available that can serve as a general
purpose object codec and support encoding of a mixture of object types:
:class:`numcodecs.JSON`, :class:`numcodecs.MsgPack`. and :class:`numcodecs.Pickle`.

For example, using the JSON codec::

    >>> z = zarr.empty(5, dtype=object, object_codec=numcodecs.JSON())
    >>> z[0] = 42
    >>> z[1] = 'foo'
    >>> z[2] = ['bar', 'baz', 'qux']
    >>> z[3] = {'a': 1, 'b': 2.2}
    >>> z[:]
    array([42, 'foo', list(['bar', 'baz', 'qux']), {'a': 1, 'b': 2.2}, None], dtype=object)

Not all codecs support encoding of all object types. The
:class:`numcodecs.Pickle` codec is the most flexible, supporting encoding any type
of Python object. However, if you are sharing data with anyone other than yourself, then
Pickle is not recommended as it is a potential security risk. This is because malicious
code can be embedded within pickled data. The JSON and MsgPack codecs do not have any
security issues and support encoding of unicode strings, lists and dictionaries.
MsgPack is usually faster for both encoding and decoding.

Ragged arrays
~~~~~~~~~~~~~

If you need to store an array of arrays, where each member array can be of any length
and stores the same primitive type (a.k.a. a ragged array), the
:class:`numcodecs.VLenArray` codec can be used, e.g.::

    >>> z = zarr.empty(4, dtype=object, object_codec=numcodecs.VLenArray(int))
    >>> z
    <zarr.core.Array (4,) object>
    >>> z.filters
    [VLenArray(dtype='<i8')]
    >>> z[0] = np.array([1, 3, 5])
    >>> z[1] = np.array([4])
    >>> z[2] = np.array([7, 9, 14])
    >>> z[:]
    array([array([1, 3, 5]), array([4]), array([ 7,  9, 14]),
           array([], dtype=int64)], dtype=object)

As a convenience, ``dtype='array:T'`` can be used as a short-hand for
``dtype=object, object_codec=numcodecs.VLenArray('T')``, where 'T' can be any NumPy
primitive dtype such as 'i4' or 'f8'. E.g.::

    >>> z = zarr.empty(4, dtype='array:i8')
    >>> z
    <zarr.core.Array (4,) object>
    >>> z.filters
    [VLenArray(dtype='<i8')]
    >>> z[0] = np.array([1, 3, 5])
    >>> z[1] = np.array([4])
    >>> z[2] = np.array([7, 9, 14])
    >>> z[:]
    array([array([1, 3, 5]), array([4]), array([ 7,  9, 14]),
           array([], dtype=int64)], dtype=object)

.. _tutorial_chunks:

Chunk optimizations
-------------------

.. _tutorial_chunks_shape:

Chunk size and shape
~~~~~~~~~~~~~~~~~~~~

In general, chunks of at least 1 megabyte (1M) uncompressed size seem to provide
better performance, at least when using the Blosc compression library.

The optimal chunk shape will depend on how you want to access the data. E.g.,
for a 2-dimensional array, if you only ever take slices along the first
dimension, then chunk across the second dimenson. If you know you want to chunk
across an entire dimension you can use ``None`` or ``-1`` within the ``chunks``
argument, e.g.::

    >>> z1 = zarr.zeros((10000, 10000), chunks=(100, None), dtype='i4')
    >>> z1.chunks
    (100, 10000)

Alternatively, if you only ever take slices along the second dimension, then
chunk across the first dimension, e.g.::

    >>> z2 = zarr.zeros((10000, 10000), chunks=(None, 100), dtype='i4')
    >>> z2.chunks
    (10000, 100)

If you require reasonable performance for both access patterns then you need to
find a compromise, e.g.::

    >>> z3 = zarr.zeros((10000, 10000), chunks=(1000, 1000), dtype='i4')
    >>> z3.chunks
    (1000, 1000)

If you are feeling lazy, you can let Zarr guess a chunk shape for your data by
providing ``chunks=True``, although please note that the algorithm for guessing
a chunk shape is based on simple heuristics and may be far from optimal. E.g.::

    >>> z4 = zarr.zeros((10000, 10000), chunks=True, dtype='i4')
    >>> z4.chunks
    (625, 625)

If you know you are always going to be loading the entire array into memory, you
can turn off chunks by providing ``chunks=False``, in which case there will be
one single chunk for the array::

    >>> z5 = zarr.zeros((10000, 10000), chunks=False, dtype='i4')
    >>> z5.chunks
    (10000, 10000)

.. _tutorial_chunks_order:

Chunk memory layout
~~~~~~~~~~~~~~~~~~~

The order of bytes **within each chunk** of an array can be changed via the
``order`` keyword argument, to use either C or Fortran layout. For
multi-dimensional arrays, these two layouts may provide different compression
ratios, depending on the correlation structure within the data. E.g.::

    >>> a = np.arange(100000000, dtype='i4').reshape(10000, 10000).T
    >>> c = zarr.array(a, chunks=(1000, 1000))
    >>> c.info
    Type               : zarr.core.Array
    Data type          : int32
    Shape              : (10000, 10000)
    Chunk shape        : (1000, 1000)
    Order              : C
    Read-only          : False
    Compressor         : Blosc(cname='lz4', clevel=5, shuffle=SHUFFLE, blocksize=0)
    Store type         : zarr.storage.KVStore
    No. bytes          : 400000000 (381.5M)
    No. bytes stored   : 6696010 (6.4M)
    Storage ratio      : 59.7
    Chunks initialized : 100/100
    >>> f = zarr.array(a, chunks=(1000, 1000), order='F')
    >>> f.info
    Type               : zarr.core.Array
    Data type          : int32
    Shape              : (10000, 10000)
    Chunk shape        : (1000, 1000)
    Order              : F
    Read-only          : False
    Compressor         : Blosc(cname='lz4', clevel=5, shuffle=SHUFFLE, blocksize=0)
    Store type         : zarr.storage.KVStore
    No. bytes          : 400000000 (381.5M)
    No. bytes stored   : 4684636 (4.5M)
    Storage ratio      : 85.4
    Chunks initialized : 100/100

In the above example, Fortran order gives a better compression ratio. This is an
artificial example but illustrates the general point that changing the order of
bytes within chunks of an array may improve the compression ratio, depending on
the structure of the data, the compression algorithm used, and which compression
filters (e.g., byte-shuffle) have been applied.

.. _tutorial_rechunking:

Changing chunk shapes (rechunking)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes you are not free to choose the initial chunking of your input data, or
you might have data saved with chunking which is not optimal for the analysis you
have planned. In such cases it can be advantageous to re-chunk the data. For small
datasets, or when the mismatch between input and output chunks is small
such that only a few chunks of the input dataset need to be read to create each
chunk in the output array, it is sufficient to simply copy the data to a new array
with the desired chunking, e.g. ::

    >>> a = zarr.zeros((10000, 10000), chunks=(100,100), dtype='uint16', store='a.zarr')
    >>> b = zarr.array(a, chunks=(100, 200), store='b.zarr')

If the chunk shapes mismatch, however, a simple copy can lead to non-optimal data
access patterns and incur a substantial performance hit when using
file based stores. One of the most pathological examples is
switching from column-based chunking to row-based chunking e.g. ::

    >>> a = zarr.zeros((10000,10000), chunks=(10000, 1), dtype='uint16, store='a.zarr')
    >>> b = zarr.array(a, chunks=(1,10000), store='b.zarr')

which will require every chunk in the input data set to be repeatedly read when creating
each output chunk. If the entire array will fit within memory, this is simply resolved
by forcing the entire input array into memory as a numpy array before converting
back to zarr with the desired chunking. ::

    >>> a = zarr.zeros((10000,10000), chunks=(10000, 1), dtype='uint16, store='a.zarr')
    >>> b = a[...]
    >>> c = zarr.array(b, chunks=(1,10000), store='c.zarr')

For data sets which have mismatched chunks and which do not fit in memory, a
more sophisticated approach to rechunking, such as offered by the
`rechunker <https://github.com/pangeo-data/rechunker>`_ package and discussed
`here <https://medium.com/pangeo/rechunker-the-missing-link-for-chunked-array-analytics-5b2359e9dc11>`_
may offer a substantial improvement in performance.

.. _tutorial_sync:

Parallel computing and synchronization
--------------------------------------

Zarr arrays have been designed for use as the source or sink for data in
parallel computations. By data source we mean that multiple concurrent read
operations may occur. By data sink we mean that multiple concurrent write
operations may occur, with each writer updating a different region of the
array. Zarr arrays have **not** been designed for situations where multiple
readers and writers are concurrently operating on the same array.

Both multi-threaded and multi-process parallelism are possible. The bottleneck
for most storage and retrieval operations is compression/decompression, and the
Python global interpreter lock (GIL) is released wherever possible during these
operations, so Zarr will generally not block other Python threads from running.

When using a Zarr array as a data sink, some synchronization (locking) may be
required to avoid data loss, depending on how data are being updated. If each
worker in a parallel computation is writing to a separate region of the array,
and if region boundaries are perfectly aligned with chunk boundaries, then no
synchronization is required. However, if region and chunk boundaries are not
perfectly aligned, then synchronization is required to avoid two workers
attempting to modify the same chunk at the same time, which could result in data
loss.

To give a simple example, consider a 1-dimensional array of length 60, ``z``,
divided into three chunks of 20 elements each. If three workers are running and
each attempts to write to a 20 element region (i.e., ``z[0:20]``, ``z[20:40]``
and ``z[40:60]``) then each worker will be writing to a separate chunk and no
synchronization is required. However, if two workers are running and each
attempts to write to a 30 element region (i.e., ``z[0:30]`` and ``z[30:60]``)
then it is possible both workers will attempt to modify the middle chunk at the
same time, and synchronization is required to prevent data loss.

Zarr provides support for chunk-level synchronization. E.g., create an array
with thread synchronization::

    >>> z = zarr.zeros((10000, 10000), chunks=(1000, 1000), dtype='i4',
    ...                 synchronizer=zarr.ThreadSynchronizer())
    >>> z
    <zarr.core.Array (10000, 10000) int32>

This array is safe to read or write within a multi-threaded program.

Zarr also provides support for process synchronization via file locking,
provided that all processes have access to a shared file system, and provided
that the underlying file system supports file locking (which is not the case for
some networked file systems). E.g.::

    >>> synchronizer = zarr.ProcessSynchronizer('data/example.sync')
    >>> z = zarr.open_array('data/example', mode='w', shape=(10000, 10000),
    ...                     chunks=(1000, 1000), dtype='i4',
    ...                     synchronizer=synchronizer)
    >>> z
    <zarr.core.Array (10000, 10000) int32>

This array is safe to read or write from multiple processes.

When using multiple processes to parallelize reads or writes on arrays using the Blosc
compression library, it may be necessary to set ``numcodecs.blosc.use_threads = False``,
as otherwise Blosc may share incorrect global state amongst processes causing programs
to hang. See also the section on :ref:`tutorial_tips_blosc` below.

Please note that support for parallel computing is an area of ongoing research
and development. If you are using Zarr for parallel computing, we welcome
feedback, experience, discussion, ideas and advice, particularly about issues
related to data integrity and performance.

.. _tutorial_pickle:

Pickle support
--------------

Zarr arrays and groups can be pickled, as long as the underlying store object can be
pickled. Instances of any of the storage classes provided in the :mod:`zarr.storage`
module can be pickled, as can the built-in ``dict`` class which can also be used for
storage.

Note that if an array or group is backed by an in-memory store like a ``dict`` or
:class:`zarr.storage.MemoryStore`, then when it is pickled all of the store data will be
included in the pickled data. However, if an array or group is backed by a persistent
store like a :class:`zarr.storage.DirectoryStore`, :class:`zarr.storage.ZipStore` or
:class:`zarr.storage.DBMStore` then the store data **are not** pickled. The only thing
that is pickled is the necessary parameters to allow the store to re-open any
underlying files or databases upon being unpickled.

E.g., pickle/unpickle an in-memory array::

    >>> import pickle
    >>> z1 = zarr.array(np.arange(100000))
    >>> s = pickle.dumps(z1)
    >>> len(s) > 5000  # relatively large because data have been pickled
    True
    >>> z2 = pickle.loads(s)
    >>> z1 == z2
    True
    >>> np.all(z1[:] == z2[:])
    True

E.g., pickle/unpickle an array stored on disk::

    >>> z3 = zarr.open('data/walnuts.zarr', mode='w', shape=100000, dtype='i8')
    >>> z3[:] = np.arange(100000)
    >>> s = pickle.dumps(z3)
    >>> len(s) < 200  # small because no data have been pickled
    True
    >>> z4 = pickle.loads(s)
    >>> z3 == z4
    True
    >>> np.all(z3[:] == z4[:])
    True

.. _tutorial_datetime:

Datetimes and timedeltas
------------------------

NumPy's ``datetime64`` ('M8') and ``timedelta64`` ('m8') dtypes are supported for Zarr
arrays, as long as the units are specified. E.g.::

    >>> z = zarr.array(['2007-07-13', '2006-01-13', '2010-08-13'], dtype='M8[D]')
    >>> z
    <zarr.core.Array (3,) datetime64[D]>
    >>> z[:]
    array(['2007-07-13', '2006-01-13', '2010-08-13'], dtype='datetime64[D]')
    >>> z[0]
    numpy.datetime64('2007-07-13')
    >>> z[0] = '1999-12-31'
    >>> z[:]
    array(['1999-12-31', '2006-01-13', '2010-08-13'], dtype='datetime64[D]')

.. _tutorial_tips:

Usage tips
----------

.. _tutorial_tips_copy:

Copying large arrays
~~~~~~~~~~~~~~~~~~~~

Data can be copied between large arrays without needing much memory, e.g.::

    >>> z1 = zarr.empty((10000, 10000), chunks=(1000, 1000), dtype='i4')
    >>> z1[:] = 42
    >>> z2 = zarr.empty_like(z1)
    >>> z2[:] = z1

Internally the example above works chunk-by-chunk, extracting only the data from
``z1`` required to fill each chunk in ``z2``. The source of the data (``z1``)
could equally be an h5py Dataset.

.. _tutorial_tips_blosc:

Configuring Blosc
~~~~~~~~~~~~~~~~~

The Blosc compressor is able to use multiple threads internally to accelerate
compression and decompression. By default, Blosc uses up to 8
internal threads. The number of Blosc threads can be changed to increase or
decrease this number, e.g.::

    >>> from numcodecs import blosc
    >>> blosc.set_nthreads(2)  # doctest: +SKIP
    8

When a Zarr array is being used within a multi-threaded program, Zarr
automatically switches to using Blosc in a single-threaded
"contextual" mode. This is generally better as it allows multiple
program threads to use Blosc simultaneously and prevents CPU thrashing
from too many active threads. If you want to manually override this
behaviour, set the value of the ``blosc.use_threads`` variable to
``True`` (Blosc always uses multiple internal threads) or ``False``
(Blosc always runs in single-threaded contextual mode). To re-enable
automatic switching, set ``blosc.use_threads`` to ``None``.

Please note that if Zarr is being used within a multi-process program, Blosc may not
be safe to use in multi-threaded mode and may cause the program to hang. If using Blosc
in a multi-process program then it is recommended to set ``blosc.use_threads = False``.
.. zarr documentation master file, created by
   sphinx-quickstart on Mon May  2 21:40:09 2016.

Zarr
====

Zarr is a format for the storage of chunked, compressed, N-dimensional arrays.
These documents describe the Zarr format and its Python implementation.

Highlights
----------

* Create N-dimensional arrays with any NumPy dtype.
* Chunk arrays along any dimension.
* Compress and/or filter chunks using any NumCodecs_ codec.
* Store arrays in memory, on disk, inside a Zip file, on S3, ...
* Read an array concurrently from multiple threads or processes.
* Write to an array concurrently from multiple threads or processes.
* Organize arrays into hierarchies via groups.

Status
------

Zarr is still a young project. Feedback and bug reports are very welcome, please get in touch via
the `GitHub issue tracker <https://github.com/zarr-developers/zarr-python/issues>`_. See
:doc:`contributing` for further information about contributing to Zarr.

Installation
------------

Zarr depends on NumPy. It is generally best to `install NumPy
<https://numpy.org/doc/stable/user/install.html>`_ first using whatever method is most
appropriate for you operating system and Python distribution. Other dependencies should be
installed automatically if using one of the installation methods below.

Install Zarr from PyPI::

    $ pip install zarr

Alternatively, install Zarr via conda::

    $ conda install -c conda-forge zarr

To install the latest development version of Zarr, you can use pip with the
latest GitHub master::

    $ pip install git+https://github.com/zarr-developers/zarr-python.git

To work with Zarr source code in development, install from GitHub::

    $ git clone --recursive https://github.com/zarr-developers/zarr-python.git
    $ cd zarr-python
    $ python setup.py install

To verify that Zarr has been fully installed, run the test suite::

    $ pip install pytest
    $ python -m pytest -v --pyargs zarr

Contents
--------

.. toctree::
    :maxdepth: 2

    tutorial
    api
    spec
    release
    contributing

Projects using Zarr
-------------------

If you are using Zarr, we would `love to hear about it
<https://github.com/zarr-developers/zarr-python/issues/228>`_.

Acknowledgments
---------------

The following people have contributed to the development of Zarr by contributing code,
documentation, code reviews, comments and/or ideas:

* :user:`Francesc Alted <FrancescAlted>`
* :user:`Martin Durant <martindurant>`
* :user:`Stephan Hoyer <shoyer>`
* :user:`John Kirkham <jakirkham>`
* :user:`Alistair Miles <alimanfoo>`
* :user:`Mamy Ratsimbazafy <mratsim>`
* :user:`Matthew Rocklin <mrocklin>`
* :user:`Vincent Schut <vincentschut>`
* :user:`Anthony Scopatz <scopatz>`
* :user:`Prakhar Goel <newt0311>`

Zarr is inspired by `HDF5 <https://www.hdfgroup.org/HDF5/>`_, `h5py
<http://www.h5py.org/>`_ and `bcolz <http://bcolz.blosc.org/>`_.

Development of Zarr is supported by the
`MRC Centre for Genomics and Global Health <http://www.cggh.org>`_.

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _NumCodecs: http://numcodecs.readthedocs.io/
.. _spec:

Specifications
==============

.. toctree::
    :maxdepth: 3

    spec/v1
    spec/v2
Release notes
=============

.. _unreleased:

Unreleased
----------

Enhancements
~~~~~~~~~~~~

* Allow to assign array ``fill_values`` and update metadata accordingly. :issue:`662`

* array indexing with [] (getitem and setitem) now supports fancy indexing.
  By :user:`Juan Nunez-Iglesias <jni>`; :issue:`725`.

* write_empty_chunks=False deletes chunks consisting of only fill_value.
  By :user:`Davis Bennett <d-v-b>`; :issue:`738`.

* Move metadata handling to a class.
  By :user:`Greggory Lee <grlee77>`; :issue:`839`.

* Create a Base store class for Zarr Store.
  By :user:`Greggory Lee <grlee77>`; :issue:`789`.

.. _release_2.10.3:

2.10.3
------

Bug fixes
~~~~~~~~~

* N5 keywords now emit UserWarning instead of raising a ValueError.
  By :user:`Boaz Mohar <boazmohar>`; :issue:`860`.

* blocks_to_decompress not used in read_part function.
  By :user:`Boaz Mohar <boazmohar>`; :issue:`861`.

* defines blocksize for array, updates hexdigest values.
  By :user:`Andrew Fulton <andrewfulton9>`; :issue:`867`.

* Fix test failure on Debian and conda-forge builds.
  By :user:`Josh Moore <joshmoore>`; :issue:`871`.

.. _release_2.10.2:

2.10.2
------

Bug fixes
~~~~~~~~~

* Fix NestedDirectoryStore datasets without dimension_separator metadata.
  By :user:`Josh Moore <joshmoore>`; :issue:`850`.

.. _release_2.10.1:

2.10.1
------

Bug fixes
~~~~~~~~~

* Fix regression by setting normalize_keys=False in fsstore constructor.
  By :user:`Davis Bennett <d-v-b>`; :issue:`842`.

.. _release_2.10.0:

2.10.0
------

Enhancements
~~~~~~~~~~~~

* Add N5FSStore.
  By :user:`Davis Bennett <d-v-b>`; :issue:`793`.

Bug fixes
~~~~~~~~~

* Ignore None dim_separators in save_array.
  By :user:`Josh Moore <joshmoore>`; :issue:`831`.

.. _release_2.9.5:

2.9.5
-----

Bug fixes
~~~~~~~~~

* Fix FSStore.listdir behavior for nested directories.
  By :user:`Greggory Lee <grlee77>`; :issue:`802`.

.. _release_2.9.4:

2.9.4
-----

Bug fixes
~~~~~~~~~

* Fix structured arrays that contain objects
  By :user: `Attila Bergou <abergou>`; :issue: `806`

.. _release_2.9.3:

2.9.3
-----

Maintenance
~~~~~~~~~~~

* Mark the fact that some tests that require ``fsspec``, without compromising the code coverage score.
  By :user:`Ben Williams <benjaminhwilliams>`; :issue:`823`.

* Only inspect alternate node type if desired isn't present.
  By :user:`Trevor Manz <manzt>`; :issue:`696`.

.. _release_2.9.2:

2.9.2
-----

Maintenance
~~~~~~~~~~~

* Correct conda-forge deployment of Zarr by fixing some Zarr tests.
  By :user:`Ben Williams <benjaminhwilliams>`; :issue:`821`.

.. _release_2.9.1:

2.9.1
-----

Maintenance
~~~~~~~~~~~

* Correct conda-forge deployment of Zarr.
  By :user:`Josh Moore <joshmoore>`; :issue:`XXX`.

.. _release_2.9.0:

2.9.0
-----

This release of Zarr Python is the first release of Zarr to not support Python 3.6.

Enhancements
~~~~~~~~~~~~

* Update ABSStore for compatibility with newer `azure.storage.blob`.
  By :user:`Tom Augspurger <TomAugspurger>`; :issue:`759`.

* Pathlib support.
  By :user:`Chris Barnes <clbarnes>`; :issue:`768`.

Documentation
~~~~~~~~~~~~~

* Clarify that arbitrary key/value pairs are OK for attributes.
  By :user:`Stephan Hoyer <shoyer>`; :issue:`751`.

* Clarify how to manually convert a DirectoryStore to a ZipStore.
  By :user:`pmav99 <pmav99>`; :issue:`763`.

Bug fixes
~~~~~~~~~

* Fix dimension_separator support.
  By :user:`Josh Moore <joshmoore>`; :issue:`775`.

* Extract ABSStore to zarr._storage.absstore.
  By :user:`Josh Moore <joshmoore>`; :issue:`781`.

* avoid NumPy 1.21.0 due to https://github.com/numpy/numpy/issues/19325
  By :user:`Greggory Lee <grlee77>`; :issue:`791`.

Maintenance
~~~~~~~~~~~

* Drop 3.6 builds.
  By :user:`Josh Moore <joshmoore>`; :issue:`774`, :issue:`778`.

* Fix build with Sphinx 4.
  By :user:`Elliott Sales de Andrade <QuLogic>`; :issue:`799`.

* TST: add missing assert in test_hexdigest.
  By :user:`Greggory Lee <grlee77>`; :issue:`801`.

.. _release_2.8.3:

2.8.3
-----

Bug fixes
~~~~~~~~~

* FSStore: default to normalize_keys=False
  By :user:`Josh Moore <joshmoore>`; :issue:`755`.
* ABSStore: compatibility with ``azure.storage.python>=12`` 
  By :user:`Tom Augspurger <tomaugspurger>`; :issue:`618`


.. _release_2.8.2:

2.8.2
-----

Documentation
~~~~~~~~~~~~~

* Add section on rechunking to tutorial
  By :user:`David Baddeley <David-Baddeley>`; :issue:`730`.

Bug fixes
~~~~~~~~~

* Expand FSStore tests and fix implementation issues
  By :user:`Davis Bennett <d-v-b>`; :issue:`709`.

Maintenance
~~~~~~~~~~~

* Updated ipytree warning for jlab3
  By :user:`Ian Hunt-Isaak <ianhi>`; :issue:`721`.

* b170a48a - (issue-728, copy-nested) Updated ipytree warning for jlab3 (#721) (3 weeks ago) <Ian Hunt-Isaak>
* Activate dependabot
  By :user:`Josh Moore <joshmoore>`; :issue:`734`.

* Update Python classifiers (Zarr is stable!)
  By :user:`Josh Moore <joshmoore>`; :issue:`731`.

.. _release_2.8.1:

2.8.1
-----

Bug fixes
~~~~~~~~~

* raise an error if create_dataset's dimension_separator is inconsistent
  By :user:`Gregory R. Lee <grlee77>`; :issue:`724`.

.. _release_2.8.0:

2.8.0
-----

V2 Specification Update
~~~~~~~~~~~~~~~~~~~~~~~

* Introduce optional dimension_separator .zarray key for nested chunks.
  By :user:`Josh Moore <joshmoore>`; :issue:`715`, :issue:`716`.

.. _release_2.7.1:

2.7.1
-----

Bug fixes
~~~~~~~~~

* Update Array to respect FSStore's key_separator  (#718)
  By :user:`Gregory R. Lee <grlee77>`; :issue:`718`.

.. _release_2.7.0:

2.7.0
-----

Enhancements
~~~~~~~~~~~~

* Start stop for iterator (`islice()`)
  By :user:`Sebastian Grill <yetyetanotherusername>`; :issue:`621`.

* Add capability to partially read and decompress chunks
  By :user:`Andrew Fulton <andrewfulton9>`; :issue:`667`.

Bug fixes
~~~~~~~~~

* Make DirectoryStore __setitem__ resilient against antivirus file locking
  By :user:`Eric Younkin <ericgyounkin>`; :issue:`698`.

* Compare test data's content generally
  By :user:`John Kirkham <jakirkham>`; :issue:`436`.

* Fix dtype usage in zarr/meta.py
  By :user:`Josh Moore <joshmoore>`; :issue:`700`.

* Fix FSStore key_seperator usage
  By :user:`Josh Moore <joshmoore>`; :issue:`669`.

* Simplify text handling in DB Store
  By :user:`John Kirkham <jakirkham>`; :issue:`670`.

* GitHub Actions migration
  By :user:`Matthias Bussonnier <Carreau>`;
  :issue:`641`, :issue:`671`, :issue:`674`, :issue:`676`, :issue:`677`, :issue:`678`,
  :issue:`679`, :issue:`680`, :issue:`682`, :issue:`684`, :issue:`685`, :issue:`686`,
  :issue:`687`, :issue:`695`, :issue:`706`.

.. _release_2.6.1:

2.6.1
-----

* Minor build fix
  By :user:`Matthias Bussonnier <Carreau>`; :issue:`666`.

.. _release_2.6.0:

2.6.0
-----

This release of Zarr Python is the first release of Zarr to not support Python 3.5.

* End Python 3.5 support.
  By :user:`Chris Barnes <clbarnes>`; :issue:`602`.

* Fix ``open_group/open_array`` to allow opening of read-only store with
  ``mode='r'`` :issue:`269`

* Add `Array` tests for FSStore.
  By :user:`Andrew Fulton <andrewfulton9>`; :issue: `644`.

* fix a bug in which ``attrs`` would not be copied on the root when using ``copy_all``; :issue:`613`

* Fix ``FileNotFoundError``  with dask/s3fs :issue:`649`

* Fix flaky fixture in test_storage.py :issue:`652`

* Fix FSStore getitems fails with arrays that have a 0 length shape dimension :issue:`644`

* Use async to fetch/write result concurrently when possible. :issue:`536`, See `this comment
  <https://github.com/zarr-developers/zarr-python/issues/536#issuecomment-721253094>`_ for some performance analysis
  showing order of magnitude faster response in some benchmark.

See `this link <https://github.com/zarr-developers/zarr-python/milestone/11?closed=1>` for the full list of closed and
merged PR tagged with the 2.6 milestone.

* Add ability to partially read and decompress arrays, see :issue:`667`. It is
  only available to chunks stored using fsspec and using Blosc as a compressor.

  For certain analysis case when only a small portion of chunks is needed it can
  be advantageous to only access and decompress part of the chunks. Doing
  partial read and decompression add high latency to many of the operation so
  should be used only when the subset of the data is small compared to the full
  chunks and is stored contiguously (that is to say either last dimensions for C
  layout, firsts for F). Pass ``partial_decompress=True`` as argument when
  creating an ``Array``, or when using ``open_array``. No option exists yet to
  apply partial read and decompress on a per-operation basis.

.. _release_2.5.0:

2.5.0
-----

This release will be the last to support Python 3.5, next version of Zarr will be Python 3.6+.

* `DirectoryStore` now uses `os.scandir`, which should make listing large store
  faster, :issue:`563`
  
* Remove a few remaining Python 2-isms.
  By :user:`Poruri Sai Rahul <rahulporuri>`; :issue:`393`.

* Fix minor bug in `N5Store`.
  By :user:`gsakkis`, :issue:`550`.

* Improve error message in Jupyter when trying to use the ``ipytree`` widget
  without ``ipytree`` installed.
  By :user:`Zain Patel <mzjp2>`; :issue:`537`

* Add typing information to many of the core functions :issue:`589`

* Explicitly close stores during testing.
  By :user:`Elliott Sales de Andrade <QuLogic>`; :issue:`442`

* Many of the convenience functions to emit errors (``err_*`` from
  ``zarr.errors``  have been replaced by ``ValueError`` subclasses. The corresponding
  ``err_*`` function have been removed. :issue:`590`, :issue:`614`)

* Improve consistency of terminology regarding arrays and datasets in the 
  documentation.
  By :user:`Josh Moore <joshmoore>`; :issue:`571`.

* Added support for generic URL opening by ``fsspec``, where the URLs have the
  form "protocol://[server]/path" or can be chained URls with "::" separators.
  The additional argument ``storage_options`` is passed to the backend, see
  the ``fsspec`` docs.
  By :user:`Martin Durant <martindurant>`; :issue:`546`

* Added support for fetching multiple items via ``getitems`` method of a
  store, if it exists. This allows for concurrent fetching of data blocks
  from stores that implement this; presently HTTP, S3, GCS. Currently only
  applies to reading.
  By :user:`Martin Durant <martindurant>`; :issue:`606`

* Efficient iteration expanded with option to pass start and stop index via
  ``array.islice``.
  By :user:`Sebastian Grill <yetyetanotherusername>`, :issue:`615`.

.. _release_2.4.0:

2.4.0
-----

Enhancements
~~~~~~~~~~~~

* Add key normalization option for ``DirectoryStore``, ``NestedDirectoryStore``,
  ``TempStore``, and ``N5Store``.
  By :user:`James Bourbeau <jrbourbeau>`; :issue:`459`.

* Add ``recurse`` keyword to ``Group.array_keys`` and ``Group.arrays`` methods.
  By :user:`James Bourbeau <jrbourbeau>`; :issue:`458`.

* Use uniform chunking for all dimensions when specifying ``chunks`` as an integer.
  Also adds support for specifying ``-1`` to chunk across an entire dimension.
  By :user:`James Bourbeau <jrbourbeau>`; :issue:`456`.

* Rename ``DictStore`` to ``MemoryStore``.
  By :user:`James Bourbeau <jrbourbeau>`; :issue:`455`.

* Rewrite ``.tree()`` pretty representation to use ``ipytree``.
  Allows it to work in both the Jupyter Notebook and JupyterLab.
  By :user:`John Kirkham <jakirkham>`; :issue:`450`.

* Do not rename Blosc parameters in n5 backend and add `blocksize` parameter,
  compatible with n5-blosc. By :user:`axtimwalde`, :issue:`485`.

* Update ``DirectoryStore`` to create files with more permissive permissions.
  By :user:`Eduardo Gonzalez <eddienko>` and :user:`James Bourbeau <jrbourbeau>`; :issue:`493`

* Use ``math.ceil`` for scalars.
  By :user:`John Kirkham <jakirkham>`; :issue:`500`.

* Ensure contiguous data using ``astype``.
  By :user:`John Kirkham <jakirkham>`; :issue:`513`.

* Refactor out ``_tofile``/``_fromfile`` from ``DirectoryStore``.
  By :user:`John Kirkham <jakirkham>`; :issue:`503`.

* Add ``__enter__``/``__exit__`` methods to ``Group`` for ``h5py.File`` compatibility.
  By :user:`Chris Barnes <clbarnes>`; :issue:`509`.

Bug fixes
~~~~~~~~~

* Fix Sqlite Store Wrong Modification.
  By :user:`Tommy Tran <potter420>`; :issue:`440`.

* Add intermediate step (using ``zipfile.ZipInfo`` object) to write
  inside ``ZipStore`` to solve too restrictive permission issue.
  By :user:`Raphael Dussin <raphaeldussin>`; :issue:`505`.

* Fix '/' prepend bug in ``ABSStore``.
  By :user:`Shikhar Goenka <shikharsg>`; :issue:`525`.

Documentation
~~~~~~~~~~~~~
* Fix hyperlink in ``README.md``.
  By :user:`Anderson Banihirwe <andersy005>`; :issue:`531`.

* Replace "nuimber" with "number".
  By :user:`John Kirkham <jakirkham>`; :issue:`512`.

* Fix azure link rendering in tutorial.
  By :user:`James Bourbeau <jrbourbeau>`; :issue:`507`.

* Update ``README`` file to be more detailed.
  By :user:`Zain Patel <mzjp2>`; :issue:`495`.

* Import blosc from numcodecs in tutorial.
  By :user:`James Bourbeau <jrbourbeau>`; :issue:`491`.

* Adds logo to docs.
  By :user:`James Bourbeau <jrbourbeau>`; :issue:`462`.

* Fix N5 link in tutorial.
  By :user:`James Bourbeau <jrbourbeau>`; :issue:`480`.

* Fix typo in code snippet.
  By :user:`Joe Jevnik <llllllllll>`; :issue:`461`.

* Fix URLs to point to zarr-python
  By :user:`John Kirkham <jakirkham>`; :issue:`453`.

Maintenance
~~~~~~~~~~~

* Add documentation build to CI.
  By :user:`James Bourbeau <jrbourbeau>`; :issue:`516`.

* Use ``ensure_ndarray`` in a few more places.
  By :user:`John Kirkham <jakirkham>`; :issue:`506`.

* Support Python 3.8.
  By :user:`John Kirkham <jakirkham>`; :issue:`499`.

* Require Numcodecs 0.6.4+ to use text handling functionality from it.
  By :user:`John Kirkham <jakirkham>`; :issue:`497`.

* Updates tests to use ``pytest.importorskip``.
  By :user:`James Bourbeau <jrbourbeau>`; :issue:`492`

* Removed support for Python 2.
  By :user:`jhamman`; :issue:`393`, :issue:`470`.

* Upgrade dependencies in the test matrices and resolve a
  compatibility issue with testing against the Azure Storage
  Emulator. By :user:`alimanfoo`; :issue:`468`, :issue:`467`.

* Use ``unittest.mock`` on Python 3.
  By :user:`Elliott Sales de Andrade <QuLogic>`; :issue:`426`.

* Drop ``decode`` from ``ConsolidatedMetadataStore``.
  By :user:`John Kirkham <jakirkham>`; :issue:`452`.


.. _release_2.3.2:

2.3.2
-----

Enhancements
~~~~~~~~~~~~

* Use ``scandir`` in ``DirectoryStore``'s ``getsize`` method.
  By :user:`John Kirkham <jakirkham>`; :issue:`431`.

Bug fixes
~~~~~~~~~

* Add and use utility functions to simplify reading and writing JSON.
  By :user:`John Kirkham <jakirkham>`; :issue:`429`, :issue:`430`.

* Fix ``collections``'s ``DeprecationWarning``\ s.
  By :user:`John Kirkham <jakirkham>`; :issue:`432`.

* Fix tests on big endian machines.
  By :user:`Elliott Sales de Andrade <QuLogic>`; :issue:`427`.


.. _release_2.3.1:

2.3.1
-----

Bug fixes
~~~~~~~~~

* Makes ``azure-storage-blob`` optional for testing.
  By :user:`John Kirkham <jakirkham>`; :issue:`419`, :issue:`420`.


.. _release_2.3.0:

2.3.0
-----

Enhancements
~~~~~~~~~~~~

* New storage backend, backed by Azure Blob Storage, class :class:`zarr.storage.ABSStore`.
  All data is stored as block blobs. By :user:`Shikhar Goenka <shikarsg>`,
  :user:`Tim Crone <tjcrone>` and :user:`Zain Patel <mzjp2>`; :issue:`345`.

* Add "consolidated" metadata as an experimental feature: use
  :func:`zarr.convenience.consolidate_metadata` to copy all metadata from the various
  metadata keys within a dataset hierarchy under a single key, and
  :func:`zarr.convenience.open_consolidated` to use this single key. This can greatly
  cut down the number of calls to the storage backend, and so remove a lot of overhead
  for reading remote data.
  By :user:`Martin Durant <martindurant>`, :user:`Alistair Miles <alimanfoo>`,
  :user:`Ryan Abernathey <rabernat>`, :issue:`268`, :issue:`332`, :issue:`338`.

* Support has been added for structured arrays with sub-array shape and/or nested fields. By
  :user:`Tarik Onalan <onalant>`, :issue:`111`, :issue:`296`.

* Adds the SQLite-backed :class:`zarr.storage.SQLiteStore` class enabling an
  SQLite database to be used as the backing store for an array or group.
  By :user:`John Kirkham <jakirkham>`, :issue:`368`, :issue:`365`.

* Efficient iteration over arrays by decompressing chunkwise.
  By :user:`Jerome Kelleher <jeromekelleher>`, :issue:`398`, :issue:`399`.

* Adds the Redis-backed :class:`zarr.storage.RedisStore` class enabling a
  Redis database to be used as the backing store for an array or group.
  By :user:`Joe Hamman <jhamman>`, :issue:`299`, :issue:`372`.

* Adds the MongoDB-backed :class:`zarr.storage.MongoDBStore` class enabling a
  MongoDB database to be used as the backing store for an array or group.
  By :user:`Noah D Brenowitz <nbren12>`, :user:`Joe Hamman <jhamman>`,
  :issue:`299`, :issue:`372`, :issue:`401`.

* **New storage class for N5 containers**. The :class:`zarr.n5.N5Store` has been
  added, which uses :class:`zarr.storage.NestedDirectoryStore` to support
  reading and writing from and to N5 containers.
  By :user:`Jan Funke <funkey>` and :user:`John Kirkham <jakirkham>`.

Bug fixes
~~~~~~~~~

* The implementation of the :class:`zarr.storage.DirectoryStore` class has been modified to
  ensure that writes are atomic and there are no race conditions where a chunk might appear
  transiently missing during a write operation. By :user:`sbalmer <sbalmer>`, :issue:`327`,
  :issue:`263`.

* Avoid raising in :class:`zarr.storage.DirectoryStore`'s ``__setitem__`` when file already exists.
  By :user:`Justin Swaney <jmswaney>`, :issue:`272`, :issue:`318`.

* The required version of the `Numcodecs`_ package has been upgraded
  to 0.6.2, which has enabled some code simplification and fixes a failing test involving
  msgpack encoding. By :user:`John Kirkham <jakirkham>`, :issue:`361`, :issue:`360`, :issue:`352`,
  :issue:`355`, :issue:`324`.

* Failing tests related to pickling/unpickling have been fixed. By :user:`Ryan Williams <ryan-williams>`,
  :issue:`273`, :issue:`308`.

* Corrects handling of ``NaT`` in ``datetime64`` and ``timedelta64`` in various
  compressors (by :user:`John Kirkham <jakirkham>`; :issue:`344`).

* Ensure ``DictStore`` contains only ``bytes`` to facilitate comparisons and protect against writes.
  By :user:`John Kirkham <jakirkham>`, :issue:`350`.

* Test and fix an issue (w.r.t. fill values) when storing complex data to ``Array``.
  By :user:`John Kirkham <jakirkham>`, :issue:`363`.

* Always use a ``tuple`` when indexing a NumPy ``ndarray``.
  By :user:`John Kirkham <jakirkham>`, :issue:`376`.

* Ensure when ``Array`` uses a ``dict``-based chunk store that it only contains
  ``bytes`` to facilitate comparisons and protect against writes. Drop the copy
  for the no filter/compressor case as this handles that case.
  By :user:`John Kirkham <jakirkham>`, :issue:`359`.

Maintenance
~~~~~~~~~~~

* Simplify directory creation and removal in ``DirectoryStore.rename``.
  By :user:`John Kirkham <jakirkham>`, :issue:`249`.

* CI and test environments have been upgraded to include Python 3.7, drop Python 3.4, and
  upgrade all pinned package requirements. :user:`Alistair Miles <alimanfoo>`, :issue:`308`.

* Start using pyup.io to maintain dependencies.
  :user:`Alistair Miles <alimanfoo>`, :issue:`326`.

* Configure flake8 line limit generally.
  :user:`John Kirkham <jakirkham>`, :issue:`335`.

* Add missing coverage pragmas.
  :user:`John Kirkham <jakirkham>`, :issue:`343`, :issue:`355`.

* Fix missing backslash in docs.
  :user:`John Kirkham <jakirkham>`, :issue:`254`, :issue:`353`.

* Include tests for stores' ``popitem`` and ``pop`` methods.
  By :user:`John Kirkham <jakirkham>`, :issue:`378`, :issue:`380`.

* Include tests for different compressors, endianness, and attributes.
  By :user:`John Kirkham <jakirkham>`, :issue:`378`, :issue:`380`.

* Test validity of stores' contents.
  By :user:`John Kirkham <jakirkham>`, :issue:`359`, :issue:`408`.


.. _release_2.2.0:

2.2.0
-----

Enhancements
~~~~~~~~~~~~

* **Advanced indexing**. The ``Array`` class has several new methods and
  properties that enable a selection of items in an array to be retrieved or
  updated. See the :ref:`tutorial_indexing` tutorial section for more
  information. There is also a `notebook
  <https://github.com/zarr-developers/zarr-python/blob/master/notebooks/advanced_indexing.ipynb>`_
  with extended examples and performance benchmarks. :issue:`78`, :issue:`89`,
  :issue:`112`, :issue:`172`.

* **New package for compressor and filter codecs**. The classes previously
  defined in the :mod:`zarr.codecs` module have been factored out into a
  separate package called `Numcodecs`_. The `Numcodecs`_ package also includes
  several new codec classes not previously available in Zarr, including
  compressor codecs for Zstd and LZ4. This change is backwards-compatible with
  existing code, as all codec classes defined by Numcodecs are imported into the
  :mod:`zarr.codecs` namespace. However, it is recommended to import codecs from
  the new package, see the tutorial sections on :ref:`tutorial_compress` and
  :ref:`tutorial_filters` for examples. With contributions by
  :user:`John Kirkham <jakirkham>`; :issue:`74`, :issue:`102`, :issue:`120`,
  :issue:`123`, :issue:`139`.

* **New storage class for DBM-style databases**. The
  :class:`zarr.storage.DBMStore` class enables any DBM-style database such as gdbm,
  ndbm or Berkeley DB, to be used as the backing store for an array or group. See the
  tutorial section on :ref:`tutorial_storage` for some examples. :issue:`133`,
  :issue:`186`.

* **New storage class for LMDB databases**. The :class:`zarr.storage.LMDBStore` class
  enables an LMDB "Lightning" database to be used as the backing store for an array or
  group. :issue:`192`.

* **New storage class using a nested directory structure for chunk files**. The
  :class:`zarr.storage.NestedDirectoryStore` has been added, which is similar to
  the existing :class:`zarr.storage.DirectoryStore` class but nests chunk files
  for multidimensional arrays into sub-directories. :issue:`155`, :issue:`177`.

* **New tree() method for printing hierarchies**. The ``Group`` class has a new
  :func:`zarr.hierarchy.Group.tree` method which enables a tree representation of
  a group hierarchy to be printed. Also provides an interactive tree
  representation when used within a Jupyter notebook. See the
  :ref:`tutorial_diagnostics` tutorial section for examples. By
  :user:`John Kirkham <jakirkham>`; :issue:`82`, :issue:`140`, :issue:`184`.

* **Visitor API**. The ``Group`` class now implements the h5py visitor API, see
  docs for the :func:`zarr.hierarchy.Group.visit`,
  :func:`zarr.hierarchy.Group.visititems` and
  :func:`zarr.hierarchy.Group.visitvalues` methods. By
  :user:`John Kirkham <jakirkham>`, :issue:`92`, :issue:`122`.

* **Viewing an array as a different dtype**. The ``Array`` class has a new
  :func:`zarr.core.Array.astype` method, which is a convenience that enables an
  array to be viewed as a different dtype. By :user:`John Kirkham <jakirkham>`,
  :issue:`94`, :issue:`96`.

* **New open(), save(), load() convenience functions**. The function
  :func:`zarr.convenience.open` provides a convenient way to open a persistent
  array or group, using either a ``DirectoryStore`` or ``ZipStore`` as the backing
  store. The functions :func:`zarr.convenience.save` and
  :func:`zarr.convenience.load` are also available and provide a convenient way to
  save an entire NumPy array to disk and load back into memory later. See the
  tutorial section :ref:`tutorial_persist` for examples. :issue:`104`,
  :issue:`105`, :issue:`141`, :issue:`181`.

* **IPython completions**. The ``Group`` class now implements ``__dir__()`` and
  ``_ipython_key_completions_()`` which enables tab-completion for group members
  to be used in any IPython interactive environment. :issue:`170`.

* **New info property; changes to __repr__**. The ``Group`` and
  ``Array`` classes have a new ``info`` property which can be used to print
  diagnostic information, including compression ratio where available. See the
  tutorial section on :ref:`tutorial_diagnostics` for examples. The string
  representation (``__repr__``) of these classes has been simplified to ensure
  it is cheap and quick to compute in all circumstances. :issue:`83`,
  :issue:`115`, :issue:`132`, :issue:`148`.

* **Chunk options**. When creating an array, ``chunks=False`` can be specified,
  which will result in an array with a single chunk only. Alternatively,
  ``chunks=True`` will trigger an automatic chunk shape guess. See
  :ref:`tutorial_chunks` for more on the ``chunks`` parameter. :issue:`106`,
  :issue:`107`, :issue:`183`.

* **Zero-dimensional arrays** and are now supported; by
  :user:`Prakhar Goel <newt0311>`, :issue:`154`, :issue:`161`.

* **Arrays with one or more zero-length dimensions** are now fully supported; by
  :user:`Prakhar Goel <newt0311>`, :issue:`150`, :issue:`154`, :issue:`160`.

* **The .zattrs key is now optional** and will now only be created when the first
  custom attribute is set; :issue:`121`, :issue:`200`.

* **New Group.move() method** supports moving a sub-group or array to a different
  location within the same hierarchy. By :user:`John Kirkham <jakirkham>`,
  :issue:`191`, :issue:`193`, :issue:`196`.

* **ZipStore is now thread-safe**; :issue:`194`, :issue:`192`.

* **New Array.hexdigest() method** computes an ``Array``'s hash with ``hashlib``.
  By :user:`John Kirkham <jakirkham>`, :issue:`98`, :issue:`203`.

* **Improved support for object arrays**. In previous versions of Zarr,
  creating an array with ``dtype=object`` was possible but could under certain
  circumstances lead to unexpected errors and/or segmentation faults. To make it easier
  to properly configure an object array, a new ``object_codec`` parameter has been
  added to array creation functions. See the tutorial section on :ref:`tutorial_objects`
  for more information and examples. Also, runtime checks have been added in both Zarr
  and Numcodecs so that segmentation faults are no longer possible, even with a badly
  configured array. This API change is backwards compatible and previous code that created
  an object array and provided an object codec via the ``filters`` parameter will
  continue to work, however a warning will be raised to encourage use of the
  ``object_codec`` parameter. :issue:`208`, :issue:`212`.

* **Added support for datetime64 and timedelta64 data types**;
  :issue:`85`, :issue:`215`.

* **Array and group attributes are now cached by default** to improve performance with
  slow stores, e.g., stores accessing data via the network; :issue:`220`, :issue:`218`,
  :issue:`204`.

* **New LRUStoreCache class**. The class :class:`zarr.storage.LRUStoreCache` has been
  added and provides a means to locally cache data in memory from a store that may be
  slow, e.g., a store that retrieves data from a remote server via the network;
  :issue:`223`.

* **New copy functions**. The new functions :func:`zarr.convenience.copy` and
  :func:`zarr.convenience.copy_all` provide a way to copy groups and/or arrays
  between HDF5 and Zarr, or between two Zarr groups. The
  :func:`zarr.convenience.copy_store` provides a more efficient way to copy
  data directly between two Zarr stores. :issue:`87`, :issue:`113`,
  :issue:`137`, :issue:`217`.

Bug fixes
~~~~~~~~~

* Fixed bug where ``read_only`` keyword argument was ignored when creating an
  array; :issue:`151`, :issue:`179`.

* Fixed bugs when using a ``ZipStore`` opened in 'w' mode; :issue:`158`,
  :issue:`182`.

* Fill values can now be provided for fixed-length string arrays; :issue:`165`,
  :issue:`176`.

* Fixed a bug where the number of chunks initialized could be counted
  incorrectly; :issue:`97`, :issue:`174`.

* Fixed a bug related to the use of an ellipsis (...) in indexing statements;
  :issue:`93`, :issue:`168`, :issue:`172`.

* Fixed a bug preventing use of other integer types for indexing; :issue:`143`,
  :issue:`147`.

Documentation
~~~~~~~~~~~~~

* Some changes have been made to the :ref:`spec_v2` document to clarify
  ambiguities and add some missing information. These changes do not break compatibility
  with any of the material as previously implemented, and so the changes have been made
  in-place in the document without incrementing the document version number. See the
  section on :ref:`spec_v2_changes` in the specification document for more information.
* A new :ref:`tutorial_indexing` section has been added to the tutorial.
* A new :ref:`tutorial_strings` section has been added to the tutorial
  (:issue:`135`, :issue:`175`).
* The :ref:`tutorial_chunks` tutorial section has been reorganised and updated.
* The :ref:`tutorial_persist` and :ref:`tutorial_storage` tutorial sections have
  been updated with new examples (:issue:`100`, :issue:`101`, :issue:`103`).
* A new tutorial section on :ref:`tutorial_pickle` has been added (:issue:`91`).
* A new tutorial section on :ref:`tutorial_datetime` has been added.
* A new tutorial section on :ref:`tutorial_diagnostics` has been added.
* The tutorial sections on :ref:`tutorial_sync` and :ref:`tutorial_tips_blosc` have been
  updated to provide information about how to avoid program hangs when using the Blosc
  compressor with multiple processes (:issue:`199`, :issue:`201`).

Maintenance
~~~~~~~~~~~

* A data fixture has been included in the test suite to ensure data format
  compatibility is maintained; :issue:`83`, :issue:`146`.
* The test suite has been migrated from nosetests to pytest; :issue:`189`, :issue:`225`.
* Various continuous integration updates and improvements; :issue:`118`, :issue:`124`,
  :issue:`125`, :issue:`126`, :issue:`109`, :issue:`114`, :issue:`171`.
* Bump numcodecs dependency to 0.5.3, completely remove nose dependency, :issue:`237`.
* Fix compatibility issues with NumPy 1.14 regarding fill values for structured arrays,
  :issue:`222`, :issue:`238`, :issue:`239`.

Acknowledgments
~~~~~~~~~~~~~~~

Code was contributed to this release by :user:`Alistair Miles <alimanfoo>`, :user:`John
Kirkham <jakirkham>` and :user:`Prakhar Goel <newt0311>`.

Documentation was contributed to this release by :user:`Mamy Ratsimbazafy <mratsim>`
and :user:`Charles Noyes <CSNoyes>`.

Thank you to :user:`John Kirkham <jakirkham>`, :user:`Stephan Hoyer <shoyer>`,
:user:`Francesc Alted <FrancescAlted>`, and :user:`Matthew Rocklin <mrocklin>` for code
reviews and/or comments on pull requests.

.. _release_2.1.4:

2.1.4
-----

* Resolved an issue where calling ``hasattr`` on a ``Group`` object erroneously
  returned a ``KeyError``. By :user:`Vincent Schut <vincentschut>`; :issue:`88`,
  :issue:`95`.

.. _release_2.1.3:

2.1.3
-----

* Resolved an issue with :func:`zarr.creation.array` where dtype was given as
  None (:issue:`80`).

.. _release_2.1.2:

2.1.2
-----

* Resolved an issue when no compression is used and chunks are stored in memory
  (:issue:`79`).

.. _release_2.1.1:

2.1.1
-----

Various minor improvements, including: ``Group`` objects support member access
via dot notation (``__getattr__``); fixed metadata caching for ``Array.shape``
property and derivatives; added ``Array.ndim`` property; fixed
``Array.__array__`` method arguments; fixed bug in pickling ``Array`` state;
fixed bug in pickling ``ThreadSynchronizer``.

.. _release_2.1.0:

2.1.0
-----

* Group objects now support member deletion via ``del`` statement
  (:issue:`65`).
* Added :class:`zarr.storage.TempStore` class for convenience to provide
  storage via a temporary directory
  (:issue:`59`).
* Fixed performance issues with :class:`zarr.storage.ZipStore` class
  (:issue:`66`).
* The Blosc extension has been modified to return bytes instead of array
  objects from compress and decompress function calls. This should
  improve compatibility and also provides a small performance increase for
  compressing high compression ratio data
  (:issue:`55`).
* Added ``overwrite`` keyword argument to array and group creation methods
  on the :class:`zarr.hierarchy.Group` class
  (:issue:`71`).
* Added ``cache_metadata`` keyword argument to array creation methods.
* The functions :func:`zarr.creation.open_array` and
  :func:`zarr.hierarchy.open_group` now accept any store as first argument
  (:issue:`56`).

.. _release_2.0.1:

2.0.1
-----

The bundled Blosc library has been upgraded to version 1.11.1.

.. _release_2.0.0:

2.0.0
-----

Hierarchies
~~~~~~~~~~~

Support has been added for organizing arrays into hierarchies via groups. See
the tutorial section on :ref:`tutorial_groups` and the :mod:`zarr.hierarchy`
API docs for more information.

Filters
~~~~~~~

Support has been added for configuring filters to preprocess chunk data prior
to compression. See the tutorial section on :ref:`tutorial_filters` and the
:mod:`zarr.codecs` API docs for more information.

Other changes
~~~~~~~~~~~~~

To accommodate support for hierarchies and filters, the Zarr metadata format
has been modified. See the :ref:`spec_v2` for more information. To migrate an
array stored using Zarr version 1.x, use the :func:`zarr.storage.migrate_1to2`
function.

The bundled Blosc library has been upgraded to version 1.11.0.

Acknowledgments
~~~~~~~~~~~~~~~

Thanks to :user:`Matthew Rocklin <mrocklin>`, :user:`Stephan Hoyer <shoyer>` and
:user:`Francesc Alted <FrancescAlted>` for contributions and comments.

.. _release_1.1.0:

1.1.0
-----

* The bundled Blosc library has been upgraded to version 1.10.0. The 'zstd'
  internal compression library is now available within Blosc. See the tutorial
  section on :ref:`tutorial_compress` for an example.
* When using the Blosc compressor, the default internal compression library
  is now 'lz4'.
* The default number of internal threads for the Blosc compressor has been
  increased to a maximum of 8 (previously 4).
* Added convenience functions :func:`zarr.blosc.list_compressors` and
  :func:`zarr.blosc.get_nthreads`.

.. _release_1.0.0:

1.0.0
-----

This release includes a complete re-organization of the code base. The
major version number has been bumped to indicate that there have been
backwards-incompatible changes to the API and the on-disk storage
format. However, Zarr is still in an early stage of development, so
please do not take the version number as an indicator of maturity.

Storage
~~~~~~~

The main motivation for re-organizing the code was to create an
abstraction layer between the core array logic and data storage (:issue:`21`).
In this release, any
object that implements the ``MutableMapping`` interface can be used as
an array store. See the tutorial sections on :ref:`tutorial_persist`
and :ref:`tutorial_storage`, the :ref:`spec_v1`, and the
:mod:`zarr.storage` module documentation for more information.

Please note also that the file organization and file name conventions
used when storing a Zarr array in a directory on the file system have
changed. Persistent Zarr arrays created using previous versions of the
software will not be compatible with this version. See the
:mod:`zarr.storage` API docs and the :ref:`spec_v1` for more
information.

Compression
~~~~~~~~~~~

An abstraction layer has also been created between the core array
logic and the code for compressing and decompressing array
chunks. This release still bundles the c-blosc library and uses Blosc
as the default compressor, however other compressors including zlib,
BZ2 and LZMA are also now supported via the Python standard
library. New compressors can also be dynamically registered for use
with Zarr. See the tutorial sections on :ref:`tutorial_compress` and
:ref:`tutorial_tips_blosc`, the :ref:`spec_v1`, and the
:mod:`zarr.compressors` module documentation for more information.

Synchronization
~~~~~~~~~~~~~~~

The synchronization code has also been refactored to create a layer of
abstraction, enabling Zarr arrays to be used in parallel computations
with a number of alternative synchronization methods. For more
information see the tutorial section on :ref:`tutorial_sync` and the
:mod:`zarr.sync` module documentation.

Changes to the Blosc extension
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NumPy is no longer a build dependency for the :mod:`zarr.blosc` Cython
extension, so setup.py will run even if NumPy is not already
installed, and should automatically install NumPy as a runtime
dependency. Manual installation of NumPy prior to installing Zarr is
still recommended, however, as the automatic installation of NumPy may
fail or be sub-optimal on some platforms.

Some optimizations have been made within the :mod:`zarr.blosc`
extension to avoid unnecessary memory copies, giving a ~10-20%
performance improvement for multi-threaded compression operations.

The :mod:`zarr.blosc` extension now automatically detects whether it
is running within a single-threaded or multi-threaded program and
adapts its internal behaviour accordingly (:issue:`27`). There is no need for
the user to make any API calls to switch Blosc between contextual and
non-contextual (global lock) mode. See also the tutorial section on
:ref:`tutorial_tips_blosc`.

Other changes
~~~~~~~~~~~~~

The internal code for managing chunks has been rewritten to be more
efficient. Now no state is maintained for chunks outside of the array
store, meaning that chunks do not carry any extra memory overhead not
accounted for by the store. This negates the need for the "lazy"
option present in the previous release, and this has been removed.

The memory layout within chunks can now be set as either "C"
(row-major) or "F" (column-major), which can help to provide better
compression for some data (:issue:`7`). See the tutorial
section on :ref:`tutorial_chunks_order` for more information.

A bug has been fixed within the ``__getitem__`` and ``__setitem__``
machinery for slicing arrays, to properly handle getting and setting
partial slices.

Acknowledgments
~~~~~~~~~~~~~~~

Thanks to :user:`Matthew Rocklin <mrocklin>`, :user:`Stephan Hoyer <shoyer>`,
:user:`Francesc Alted <FrancescAlted>`, :user:`Anthony Scopatz <scopatz>` and
:user:`Martin Durant <martindurant>` for contributions and comments.

.. _release_0.4.0:

0.4.0
-----

See `v0.4.0 release notes on GitHub
<https://github.com/zarr-developers/zarr-python/releases/tag/v0.4.0>`_.

.. _release_0.3.0:

0.3.0
-----

See `v0.3.0 release notes on GitHub
<https://github.com/zarr-developers/zarr-python/releases/tag/v0.3.0>`_.

.. _Numcodecs: http://numcodecs.readthedocs.io/
API reference
=============

.. toctree::
    :maxdepth: 3

    api/creation
    api/core
    api/hierarchy
    api/storage
    api/n5
    api/convenience
    api/codecs
    api/attrs
    api/sync
Zarr - scalable storage of tensor data for use in parallel and distributed computing
====================================================================================

SciPy 2019 submission.


Short summary
-------------

Many scientific problems involve computing over large N-dimensional
typed arrays of data, and reading or writing data is often the major
bottleneck limiting speed or scalability. The Zarr project is
developing a simple, scalable approach to storage of such data in a
way that is compatible with a range of approaches to distributed and
parallel computing. We describe the Zarr protocol and data storage
format, and the current state of implementations for various
programming languages including Python. We also describe current uses
of Zarr in malaria genomics, the Human Cell Atlas, and the Pangeo
project.


Abstract
--------

Background
~~~~~~~~~~

Across a broad range of scientific disciplines, data are naturally
represented and stored as N-dimensional typed arrays, also known as
tensors. The volume of data being generated is outstripping our
ability to analyse it, and scientific communities are looking for ways
to leverage modern multi-core CPUs and distributed computing
platforms, including cloud computing. Retrieval and storage of data is
often the major bottleneck, and new approaches to data storage are
needed to accelerate distributed computations and enable them to scale
on a variety of platforms.

Methods
~~~~~~~

We have designed a new storage format and protocol for tensor data
[1_], and have released an open source Python implementation [2_,
3_]. Our approach builds on data storage concepts from HDF5 [4_],
particularly chunking and compression, and hierarchical organisation
of datasets. Key design goals include: a simple protocol and format
that can be implemented in other programming languages; support for
multiple concurrent readers or writers; support for a variety of
parallel computing environments, from multi-threaded execution on a
single CPU to multi-process execution across a multi-node cluster;
pluggable storage subsystem with support for file systems, key-value
databases and cloud object stores; pluggable encoding subsystem with
support for a variety of modern compressors.

Results
~~~~~~~

We illustrate the use of Zarr with examples from several scientific
domains. Zarr is being used within the Pangeo project [5_], which is
building a community platform for big data geoscience. The Pangeo
community have converted a number of existing climate modelling and
satellite observation datasets to Zarr [6_], and have demonstrated
their use in computations using HPC and cloud computing
environments. Within the MalariaGEN project [7_], Zarr is used to
store genome variation data from next-generation sequencing of natural
populations of malaria parasites and mosquitoes [8_] and these data
are used as input to analyses of the evolution of these organisms in
response to selective pressure from anti-malarial drugs and
insecticides. Zarr is being used within the Human Cell Atlas (HCA)
project [9_], which is building a reference atlas of healthy human
cell types. This project hopes to leverage this information to better
understand the dysregulation of cellular states that underly human
disease. The Human Cell Atlas uses Zarr as the output data format
because it enables the project to easily generate matrices containing
user-selected subsets of cells.

Conclusions
~~~~~~~~~~~

Zarr is generating interest across a range of scientific domains, and
work is ongoing to establish a community process to support further
development of the specifications and implementations in other
programming languages [10_, 11_, 12_] and building interoperability
with a similar project called N5 [13_]. Other packages within the
PyData ecosystem, notably Dask [14_], Xarray [15_] and Intake [16_],
have added capability to read and write Zarr, and together these
packages provide a compelling solution for large scale data science
using Python [17_]. Zarr has recently been presented in several
venues, including a webinar for the ESIP Federation tech dive series
[18_], and a talk at the AGU Fall Meeting 2018 [19_].


References
~~~~~~~~~~

.. _1: https://zarr.readthedocs.io/en/stable/spec/v2.html
.. _2: https://github.com/zarr-developers/zarr-python
.. _3: https://github.com/zarr-developers/numcodecs
.. _4: https://www.hdfgroup.org/solutions/hdf5/
.. _5: https://pangeo.io/
.. _6: https://pangeo.io/catalog.html
.. _7: https://www.malariagen.net/
.. _8: http://alimanfoo.github.io/2016/09/21/genotype-compression-benchmark.html
.. _9: https://www.humancellatlas.org/
.. _10: https://github.com/constantinpape/z5
.. _11: https://github.com/lasersonlab/ndarray.scala
.. _12: https://github.com/meggart/ZarrNative.jl
.. _13: https://github.com/saalfeldlab/n5
.. _14: http://docs.dask.org/en/latest/array-creation.html
.. _15: http://xarray.pydata.org/en/stable/io.html
.. _16: https://github.com/ContinuumIO/intake-xarray
.. _17: http://matthewrocklin.com/blog/work/2018/01/22/pangeo-2
.. _18: http://wiki.esipfed.org/index.php/Interoperability_and_Technology/Tech_Dive_Webinar_Series#8_March.2C_2018:_.22Zarr:_A_simple.2C_open.2C_scalable_solution_for_big_NetCDF.2FHDF_data_on_the_Cloud.22:_Alistair_Miles.2C_University_of_Oxford.
.. _19: https://agu.confex.com/agu/fm18/meetingapp.cgi/Paper/390015


Authors
-------

Project contributors are listed in alphabetical order by surname.

* `Ryan Abernathey <https://github.com/rabernat>`_, Columbia University 
* `Stephan Balmer <https://github.com/sbalmer>`_, Meteotest
* `Ambrose Carr <https://github.com/ambrosejcarr>`_, Chan Zuckerberg Initiative
* `Tim Crone <https://github.com/tjcrone>`_, Columbia University
* `Martin Durant <https://github.com/martindurant>`_, Anaconda, inc.
* `Jan Funke <https://github.com/funkey>`_, HHMI Janelia
* `Darren Gallagher <https://github.com/dazzag24>`_, Satavia
* `Fabian Gans <https://github.com/meggart>`_, Max Planck Institute for Biogeochemistry
* `Shikhar Goenka <https://github.com/shikharsg>`_, Satavia
* `Joe Hamman <https://github.com/jhamman>`_, NCAR
* `Stephan Hoyer <https://github.com/shoyer>`_, Google
* `Jerome Kelleher <https://github.com/jeromekelleher>`_, University of Oxford
* `John Kirkham <https://github.com/jakirkham>`_, HHMI Janelia
* `Alistair Miles <https://github.com/alimanfoo>`_, University of Oxford
* `Josh Moore <https://github.com/joshmoore>`_, University of Dundee
* `Charles Noyes <https://github.com/CSNoyes>`_, University of Southern California
* `Tarik Onalan <https://github.com/onalant>`_
* `Constantin Pape <https://github.com/constantinpape>`_, University of Heidelberg
* `Zain Patel <https://github.com/mzjp2>`_, University of Cambridge
* `Matthew Rocklin <https://github.com/mrocklin>`_, NVIDIA
* `Stephan Saafeld <https://github.com/axtimwalde>`_, HHMI Janelia
* `Vincent Schut <https://github.com/vincentschut>`_, Satelligence
* `Justin Swaney <https://github.com/jmswaney>`_, MIT
* `Ryan Williams <https://github.com/ryan-williams>`_, Chan Zuckerberg Initiative
The Attributes class (``zarr.attrs``)
=====================================
.. module:: zarr.attrs

.. autoclass:: Attributes

    .. automethod:: __getitem__
    .. automethod:: __setitem__
    .. automethod:: __delitem__
    .. automethod:: __iter__
    .. automethod:: __len__
    .. automethod:: keys
    .. automethod:: asdict
    .. automethod:: put
    .. automethod:: update
    .. automethod:: refresh
Storage (``zarr.storage``)
==========================
.. automodule:: zarr.storage

.. autoclass:: MemoryStore
.. autoclass:: DirectoryStore
.. autoclass:: TempStore
.. autoclass:: NestedDirectoryStore
.. autoclass:: ZipStore

    .. automethod:: close
    .. automethod:: flush

.. autoclass:: DBMStore

    .. automethod:: close
    .. automethod:: flush

.. autoclass:: LMDBStore

    .. automethod:: close
    .. automethod:: flush

.. autoclass:: SQLiteStore

    .. automethod:: close

.. autoclass:: MongoDBStore
.. autoclass:: RedisStore
.. autoclass:: LRUStoreCache

    .. automethod:: invalidate
    .. automethod:: invalidate_values
    .. automethod:: invalidate_keys

.. autoclass:: ABSStore

.. autoclass:: FSStore

.. autoclass:: ConsolidatedMetadataStore

.. autofunction:: init_array
.. autofunction:: init_group
.. autofunction:: contains_array
.. autofunction:: contains_group
.. autofunction:: listdir
.. autofunction:: rmdir
.. autofunction:: getsize
.. autofunction:: rename
.. autofunction:: migrate_1to2
Compressors and filters (``zarr.codecs``)
=========================================
.. module:: zarr.codecs

This module contains compressor and filter classes for use with Zarr. Please note that this module
is provided for backwards compatibility with previous versions of Zarr. From Zarr version 2.2
onwards, all codec classes have been moved to a separate package called Numcodecs_. The two
packages (Zarr and Numcodecs_) are designed to be used together. For example, a Numcodecs_ codec
class can be used as a compressor for a Zarr array::

    >>> import zarr
    >>> from numcodecs import Blosc
    >>> z = zarr.zeros(1000000, compressor=Blosc(cname='zstd', clevel=1, shuffle=Blosc.SHUFFLE))

Codec classes can also be used as filters. See the tutorial section on :ref:`tutorial_filters`
for more information.

Please note that it is also relatively straightforward to define and register custom codec
classes. See the Numcodecs `codec API <http://numcodecs.readthedocs.io/en/latest/abc.html>`_ and
`codec registry <http://numcodecs.readthedocs.io/en/latest/registry.html>`_ documentation for more
information.

.. _Numcodecs: http://numcodecs.readthedocs.io/
Convenience functions (``zarr.convenience``)
============================================
.. automodule:: zarr.convenience
.. autofunction:: open
.. autofunction:: save
.. autofunction:: load
.. autofunction:: save_array
.. autofunction:: save_group
.. autofunction:: copy
.. autofunction:: copy_all
.. autofunction:: copy_store
.. autofunction:: tree
.. autofunction:: consolidate_metadata
.. autofunction:: open_consolidated
The Array class (``zarr.core``)
===============================
.. module:: zarr.core

.. autoclass:: Array

    .. automethod:: __getitem__
    .. automethod:: __setitem__
    .. automethod:: get_basic_selection
    .. automethod:: set_basic_selection
    .. automethod:: get_mask_selection
    .. automethod:: set_mask_selection
    .. automethod:: get_coordinate_selection
    .. automethod:: set_coordinate_selection
    .. automethod:: get_orthogonal_selection
    .. automethod:: set_orthogonal_selection
    .. automethod:: digest
    .. automethod:: hexdigest
    .. automethod:: resize
    .. automethod:: append
    .. automethod:: view
    .. automethod:: astype
Synchronization (``zarr.sync``)
===============================
.. module:: zarr.sync

.. autoclass:: ThreadSynchronizer
.. autoclass:: ProcessSynchronizer
Groups (``zarr.hierarchy``)
===========================
.. module:: zarr.hierarchy

.. autofunction:: group
.. autofunction:: open_group

.. autoclass:: Group

    .. automethod:: __len__
    .. automethod:: __iter__
    .. automethod:: __contains__
    .. automethod:: __getitem__
    .. automethod:: __enter__
    .. automethod:: __exit__
    .. automethod:: group_keys
    .. automethod:: groups
    .. automethod:: array_keys
    .. automethod:: arrays
    .. automethod:: visit
    .. automethod:: visitkeys
    .. automethod:: visitvalues
    .. automethod:: visititems
    .. automethod:: tree
    .. automethod:: create_group
    .. automethod:: require_group
    .. automethod:: create_groups
    .. automethod:: require_groups
    .. automethod:: create_dataset
    .. automethod:: require_dataset
    .. automethod:: create
    .. automethod:: empty
    .. automethod:: zeros
    .. automethod:: ones
    .. automethod:: full
    .. automethod:: array
    .. automethod:: empty_like
    .. automethod:: zeros_like
    .. automethod:: ones_like
    .. automethod:: full_like
    .. automethod:: moveArray creation (``zarr.creation``)
==================================
.. module:: zarr.creation
.. autofunction:: create
.. autofunction:: empty
.. autofunction:: zeros
.. autofunction:: ones
.. autofunction:: full
.. autofunction:: array
.. autofunction:: open_array
.. autofunction:: empty_like
.. autofunction:: zeros_like
.. autofunction:: ones_like
.. autofunction:: full_like
.. autofunction:: open_like
N5 (``zarr.n5``)
================
.. automodule:: zarr.n5

.. autoclass:: N5Store
.. _spec_v1:

Zarr storage specification version 1
====================================

This document provides a technical specification of the protocol and
format used for storing a Zarr array. The key words "MUST", "MUST
NOT", "REQUIRED", "SHALL", "SHALL NOT", "SHOULD", "SHOULD NOT",
"RECOMMENDED", "MAY", and "OPTIONAL" in this document are to be
interpreted as described in `RFC 2119
<https://www.ietf.org/rfc/rfc2119.txt>`_.

Status
------

This specification is deprecated. See :ref:`spec` for the latest version.

Storage
-------

A Zarr array can be stored in any storage system that provides a
key/value interface, where a key is an ASCII string and a value is an
arbitrary sequence of bytes, and the supported operations are read
(get the sequence of bytes associated with a given key), write (set
the sequence of bytes associated with a given key) and delete (remove
a key/value pair).

For example, a directory in a file system can provide this interface,
where keys are file names, values are file contents, and files can be
read, written or deleted via the operating system. Equally, an S3
bucket can provide this interface, where keys are resource names,
values are resource contents, and resources can be read, written or
deleted via HTTP.

Below an "array store" refers to any system implementing this
interface.

Metadata
--------

Each array requires essential configuration metadata to be stored,
enabling correct interpretation of the stored data. This metadata is
encoded using JSON and stored as the value of the 'meta' key within an
array store.

The metadata resource is a JSON object. The following keys MUST be
present within the object:

zarr_format
    An integer defining the version of the storage specification to which the
    array store adheres.
shape
    A list of integers defining the length of each dimension of the array.
chunks
    A list of integers defining the length of each dimension of a chunk of the
    array. Note that all chunks within a Zarr array have the same shape.
dtype
    A string or list defining a valid data type for the array. See also
    the subsection below on data type encoding.
compression
    A string identifying the primary compression library used to compress
    each chunk of the array.
compression_opts
    An integer, string or dictionary providing options to the primary
    compression library.
fill_value
    A scalar value providing the default value to use for uninitialized
    portions of the array.
order
    Either 'C' or 'F', defining the layout of bytes within each chunk of the
    array. 'C' means row-major order, i.e., the last dimension varies fastest;
    'F' means column-major order, i.e., the first dimension varies fastest.

Other keys MAY be present within the metadata object however they MUST
NOT alter the interpretation of the required fields defined above.

For example, the JSON object below defines a 2-dimensional array of
64-bit little-endian floating point numbers with 10000 rows and 10000
columns, divided into chunks of 1000 rows and 1000 columns (so there
will be 100 chunks in total arranged in a 10 by 10 grid). Within each
chunk the data are laid out in C contiguous order, and each chunk is
compressed using the Blosc compression library::

    {
        "chunks": [
            1000,
            1000
        ],
        "compression": "blosc",
        "compression_opts": {
            "clevel": 5,
            "cname": "lz4",
            "shuffle": 1
        },
        "dtype": "<f8",
        "fill_value": null,
        "order": "C",
        "shape": [
            10000,
            10000
        ],
        "zarr_format": 1
    }

Data type encoding
~~~~~~~~~~~~~~~~~~

Simple data types are encoded within the array metadata resource as a
string, following the `NumPy array protocol type string (typestr)
format
<numpy:arrays.interface>`_. The
format consists of 3 parts: a character describing the byteorder of
the data (``<``: little-endian, ``>``: big-endian, ``|``:
not-relevant), a character code giving the basic type of the array,
and an integer providing the number of bytes the type uses. The byte
order MUST be specified. E.g., ``"<f8"``, ``">i4"``, ``"|b1"`` and
``"|S12"`` are valid data types.

Structure data types (i.e., with multiple named fields) are encoded as
a list of two-element lists, following `NumPy array protocol type
descriptions (descr)
<numpy:arrays.interface>`_.
For example, the JSON list ``[["r", "|u1"], ["g", "|u1"], ["b",
"|u1"]]`` defines a data type composed of three single-byte unsigned
integers labelled 'r', 'g' and 'b'.

Chunks
------

Each chunk of the array is compressed by passing the raw bytes for the
chunk through the primary compression library to obtain a new sequence
of bytes comprising the compressed chunk data. No header is added to
the compressed bytes or any other modification made. The internal
structure of the compressed bytes will depend on which primary
compressor was used. For example, the `Blosc compressor
<https://github.com/Blosc/c-blosc/blob/master/README_HEADER.rst>`_
produces a sequence of bytes that begins with a 16-byte header
followed by compressed data.

The compressed sequence of bytes for each chunk is stored under a key
formed from the index of the chunk within the grid of chunks
representing the array. To form a string key for a chunk, the indices
are converted to strings and concatenated with the period character
('.') separating each index. For example, given an array with shape
(10000, 10000) and chunk shape (1000, 1000) there will be 100 chunks
laid out in a 10 by 10 grid. The chunk with indices (0, 0) provides
data for rows 0-1000 and columns 0-1000 and is stored under the key
'0.0'; the chunk with indices (2, 4) provides data for rows 2000-3000
and columns 4000-5000 and is stored under the key '2.4'; etc.

There is no need for all chunks to be present within an array
store. If a chunk is not present then it is considered to be in an
uninitialized state.  An uninitialized chunk MUST be treated as if it
was uniformly filled with the value of the 'fill_value' field in the
array metadata. If the 'fill_value' field is ``null`` then the
contents of the chunk are undefined.

Note that all chunks in an array have the same shape. If the length of
any array dimension is not exactly divisible by the length of the
corresponding chunk dimension then some chunks will overhang the edge
of the array. The contents of any chunk region falling outside the
array are undefined.

Attributes
----------

Each array can also be associated with custom attributes, which are
simple key/value items with application-specific meaning. Custom
attributes are encoded as a JSON object and stored under the 'attrs'
key within an array store. Even if the attributes are empty, the
'attrs' key MUST be present within an array store.

For example, the JSON object below encodes three attributes named
'foo', 'bar' and 'baz'::

    {
        "foo": 42,
        "bar": "apples",
        "baz": [1, 2, 3, 4]
    }

Example
-------

Below is an example of storing a Zarr array, using a directory on the
local file system as storage.

Initialize the store::

    >>> import zarr
    >>> store = zarr.DirectoryStore('example.zarr')
    >>> zarr.init_store(store, shape=(20, 20), chunks=(10, 10),
    ...                 dtype='i4', fill_value=42, compression='zlib',
    ...                 compression_opts=1, overwrite=True)

No chunks are initialized yet, so only the 'meta' and 'attrs' keys
have been set::

    >>> import os
    >>> sorted(os.listdir('example.zarr'))
    ['attrs', 'meta']

Inspect the array metadata::

    >>> print(open('example.zarr/meta').read())
    {
        "chunks": [
            10,
            10
        ],
        "compression": "zlib",
        "compression_opts": 1,
        "dtype": "<i4",
        "fill_value": 42,
        "order": "C",
        "shape": [
            20,
            20
        ],
        "zarr_format": 1
    }

Inspect the array attributes::

    >>> print(open('example.zarr/attrs').read())
    {}

Set some data::

    >>> z = zarr.Array(store)
    >>> z[0:10, 0:10] = 1
    >>> sorted(os.listdir('example.zarr'))
    ['0.0', 'attrs', 'meta']

Set some more data::

    >>> z[0:10, 10:20] = 2
    >>> z[10:20, :] = 3
    >>> sorted(os.listdir('example.zarr'))
    ['0.0', '0.1', '1.0', '1.1', 'attrs', 'meta']

Manually decompress a single chunk for illustration::

    >>> import zlib
    >>> b = zlib.decompress(open('example.zarr/0.0', 'rb').read())
    >>> import numpy as np
    >>> a = np.frombuffer(b, dtype='<i4')
    >>> a
    array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1], dtype=int32)

Modify the array attributes::

    >>> z.attrs['foo'] = 42
    >>> z.attrs['bar'] = 'apples'
    >>> z.attrs['baz'] = [1, 2, 3, 4]
    >>> print(open('example.zarr/attrs').read())
    {
        "bar": "apples",
        "baz": [
            1,
            2,
            3,
            4
        ],
        "foo": 42
    }
.. _spec_v2:

Zarr storage specification version 2
====================================

This document provides a technical specification of the protocol and format
used for storing Zarr arrays. The key words "MUST", "MUST NOT", "REQUIRED",
"SHALL", "SHALL NOT", "SHOULD", "SHOULD NOT", "RECOMMENDED", "MAY", and
"OPTIONAL" in this document are to be interpreted as described in `RFC 2119
<https://www.ietf.org/rfc/rfc2119.txt>`_.

Status
------

This specification is the latest version. See :ref:`spec` for previous
versions.

.. _spec_v2_storage:

Storage
-------

A Zarr array can be stored in any storage system that provides a key/value
interface, where a key is an ASCII string and a value is an arbitrary sequence
of bytes, and the supported operations are read (get the sequence of bytes
associated with a given key), write (set the sequence of bytes associated with
a given key) and delete (remove a key/value pair).

For example, a directory in a file system can provide this interface, where
keys are file names, values are file contents, and files can be read, written
or deleted via the operating system. Equally, an S3 bucket can provide this
interface, where keys are resource names, values are resource contents, and
resources can be read, written or deleted via HTTP.

Below an "array store" refers to any system implementing this interface.

.. _spec_v2_array:

Arrays
------

.. _spec_v2_array_metadata:

Metadata
~~~~~~~~

Each array requires essential configuration metadata to be stored, enabling
correct interpretation of the stored data. This metadata is encoded using JSON
and stored as the value of the ".zarray" key within an array store.

The metadata resource is a JSON object. The following keys MUST be present
within the object:

zarr_format
    An integer defining the version of the storage specification to which the
    array store adheres.
shape
    A list of integers defining the length of each dimension of the array.
chunks
    A list of integers defining the length of each dimension of a chunk of the
    array. Note that all chunks within a Zarr array have the same shape.
dtype
    A string or list defining a valid data type for the array. See also
    the subsection below on data type encoding.
compressor
    A JSON object identifying the primary compression codec and providing
    configuration parameters, or ``null`` if no compressor is to be used.
    The object MUST contain an ``"id"`` key identifying the codec to be used.
fill_value
    A scalar value providing the default value to use for uninitialized
    portions of the array, or ``null`` if no fill_value is to be used.
order
    Either "C" or "F", defining the layout of bytes within each chunk of the
    array. "C" means row-major order, i.e., the last dimension varies fastest;
    "F" means column-major order, i.e., the first dimension varies fastest.
filters
    A list of JSON objects providing codec configurations, or ``null`` if no
    filters are to be applied. Each codec configuration object MUST contain a
    ``"id"`` key identifying the codec to be used.

The following keys MAY be present within the object:

dimension_separator
    If present, either the string ``"."`` or ``"/""`` defining the separator placed
    between the dimensions of a chunk. If the value is not set, then the
    default MUST be assumed to be ``"."``, leading to chunk keys of the form "0.0".
    Arrays defined with ``"/"`` as the dimension separator can be considered to have
    nested, or hierarchical, keys of the form "0/0" that SHOULD where possible
    produce a directory-like structure.

Other keys SHOULD NOT be present within the metadata object and SHOULD be
ignored by implementations.

For example, the JSON object below defines a 2-dimensional array of 64-bit
little-endian floating point numbers with 10000 rows and 10000 columns, divided
into chunks of 1000 rows and 1000 columns (so there will be 100 chunks in total
arranged in a 10 by 10 grid). Within each chunk the data are laid out in C
contiguous order. Each chunk is encoded using a delta filter and compressed
using the Blosc compression library prior to storage::

    {
        "chunks": [
            1000,
            1000
        ],
        "compressor": {
            "id": "blosc",
            "cname": "lz4",
            "clevel": 5,
            "shuffle": 1
        },
        "dtype": "<f8",
        "fill_value": "NaN",
        "filters": [
            {"id": "delta", "dtype": "<f8", "astype": "<f4"}
        ],
        "order": "C",
        "shape": [
            10000,
            10000
        ],
        "zarr_format": 2
    }

.. _spec_v2_array_dtype:

Data type encoding
~~~~~~~~~~~~~~~~~~

Simple data types are encoded within the array metadata as a string,
following the :ref:`NumPy array protocol type string (typestr) format
<numpy:arrays.interface>`. The format
consists of 3 parts:

* One character describing the byteorder of the data (``"<"``: little-endian;
  ``">"``: big-endian; ``"|"``: not-relevant)
* One character code giving the basic type of the array (``"b"``: Boolean (integer
  type where all values are only True or False); ``"i"``: integer; ``"u"``: unsigned
  integer; ``"f"``: floating point; ``"c"``: complex floating point; ``"m"``: timedelta;
  ``"M"``: datetime; ``"S"``: string (fixed-length sequence of char); ``"U"``: unicode
  (fixed-length sequence of Py_UNICODE); ``"V"``: other (void * – each item is a
  fixed-size chunk of memory))
* An integer specifying the number of bytes the type uses.

The byte order MUST be specified. E.g., ``"<f8"``, ``">i4"``, ``"|b1"`` and
``"|S12"`` are valid data type encodings.

For datetime64 ("M") and timedelta64 ("m") data types, these MUST also include the
units within square brackets. A list of valid units and their definitions are given in
the :ref:`NumPy documentation on Datetimes and Timedeltas
<numpy:arrays.dtypes.dateunits>`.
For example, ``"<M8[ns]"`` specifies a datetime64 data type with nanosecond time units.

Structured data types (i.e., with multiple named fields) are encoded
as a list of lists, following :ref:`NumPy array protocol type descriptions
(descr)
<numpy:arrays.interface>`. Each
sub-list has the form ``[fieldname, datatype, shape]`` where ``shape``
is optional. ``fieldname`` is a string, ``datatype`` is a string
specifying a simple data type (see above), and ``shape`` is a list of
integers specifying subarray shape. For example, the JSON list below
defines a data type composed of three single-byte unsigned integer
fields named "r", "g" and "b"::

    [["r", "|u1"], ["g", "|u1"], ["b", "|u1"]]

For example, the JSON list below defines a data type composed of three
fields named "x", "y" and "z", where "x" and "y" each contain 32-bit
floats, and each item in "z" is a 2 by 2 array of floats::

    [["x", "<f4"], ["y", "<f4"], ["z", "<f4", [2, 2]]]

Structured data types may also be nested, e.g., the following JSON
list defines a data type with two fields "foo" and "bar", where "bar"
has two sub-fields "baz" and "qux"::

    [["foo", "<f4"], ["bar", [["baz", "<f4"], ["qux", "<i4"]]]] 
    
.. _spec_v2_array_fill_value:

Fill value encoding
~~~~~~~~~~~~~~~~~~~

For simple floating point data types, the following table MUST be used to
encode values of the "fill_value" field:

=================  ===============
Value              JSON encoding
=================  ===============
Not a Number       ``"NaN"``
Positive Infinity  ``"Infinity"``
Negative Infinity  ``"-Infinity"``
=================  ===============

If an array has a fixed length byte string data type (e.g., ``"|S12"``), or a
structured data type, and if the fill value is not null, then the fill value
MUST be encoded as an ASCII string using the standard Base64 alphabet.

.. _spec_v2_array_chunks:

Chunks
~~~~~~

Each chunk of the array is compressed by passing the raw bytes for the chunk
through the primary compression library to obtain a new sequence of bytes
comprising the compressed chunk data. No header is added to the compressed
bytes or any other modification made. The internal structure of the compressed
bytes will depend on which primary compressor was used. For example, the `Blosc
compressor <https://github.com/Blosc/c-blosc/blob/master/README_HEADER.rst>`_
produces a sequence of bytes that begins with a 16-byte header followed by
compressed data.

The compressed sequence of bytes for each chunk is stored under a key formed
from the index of the chunk within the grid of chunks representing the array.
To form a string key for a chunk, the indices are converted to strings and
concatenated with the period character (".") separating each index. For
example, given an array with shape (10000, 10000) and chunk shape (1000, 1000)
there will be 100 chunks laid out in a 10 by 10 grid. The chunk with indices
(0, 0) provides data for rows 0-1000 and columns 0-1000 and is stored under the
key "0.0"; the chunk with indices (2, 4) provides data for rows 2000-3000 and
columns 4000-5000 and is stored under the key "2.4"; etc.

There is no need for all chunks to be present within an array store. If a chunk
is not present then it is considered to be in an uninitialized state.  An
uninitialized chunk MUST be treated as if it was uniformly filled with the value
of the "fill_value" field in the array metadata. If the "fill_value" field is
``null`` then the contents of the chunk are undefined.

Note that all chunks in an array have the same shape. If the length of any
array dimension is not exactly divisible by the length of the corresponding
chunk dimension then some chunks will overhang the edge of the array. The
contents of any chunk region falling outside the array are undefined.

.. _spec_v2_array_filters:

Filters
~~~~~~~

Optionally a sequence of one or more filters can be used to transform chunk
data prior to compression. When storing data, filters are applied in the order
specified in array metadata to encode data, then the encoded data are passed to
the primary compressor. When retrieving data, stored chunk data are
decompressed by the primary compressor then decoded using filters in the
reverse order.

.. _spec_v2_hierarchy:

Hierarchies
-----------

.. _spec_v2_hierarchy_paths:

Logical storage paths
~~~~~~~~~~~~~~~~~~~~~

Multiple arrays can be stored in the same array store by associating each array
with a different logical path. A logical path is simply an ASCII string. The
logical path is used to form a prefix for keys used by the array. For example,
if an array is stored at logical path "foo/bar" then the array metadata will be
stored under the key "foo/bar/.zarray", the user-defined attributes will be
stored under the key "foo/bar/.zattrs", and the chunks will be stored under
keys like "foo/bar/0.0", "foo/bar/0.1", etc.

To ensure consistent behaviour across different storage systems, logical paths
MUST be normalized as follows:

* Replace all backward slash characters ("\\\\") with forward slash characters
  ("/")
* Strip any leading "/" characters
* Strip any trailing "/" characters
* Collapse any sequence of more than one "/" character into a single "/"
  character

The key prefix is then obtained by appending a single "/" character to the
normalized logical path.

After normalization, if splitting a logical path by the "/" character results
in any path segment equal to the string "." or the string ".." then an error
MUST be raised.

N.B., how the underlying array store processes requests to store values under
keys containing the "/" character is entirely up to the store implementation
and is not constrained by this specification. E.g., an array store could simply
treat all keys as opaque ASCII strings; equally, an array store could map
logical paths onto some kind of hierarchical storage (e.g., directories on a
file system).

.. _spec_v2_hierarchy_groups:

Groups
~~~~~~

Arrays can be organized into groups which can also contain other groups. A
group is created by storing group metadata under the ".zgroup" key under some
logical path. E.g., a group exists at the root of an array store if the
".zgroup" key exists in the store, and a group exists at logical path "foo/bar"
if the "foo/bar/.zgroup" key exists in the store.

If the user requests a group to be created under some logical path, then groups
MUST also be created at all ancestor paths. E.g., if the user requests group
creation at path "foo/bar" then groups MUST be created at path "foo" and the
root of the store, if they don't already exist.

If the user requests an array to be created under some logical path, then
groups MUST also be created at all ancestor paths. E.g., if the user requests
array creation at path "foo/bar/baz" then groups must be created at path
"foo/bar", path "foo", and the root of the store, if they don't already exist.

The group metadata resource is a JSON object. The following keys MUST be present
within the object:

zarr_format
    An integer defining the version of the storage specification to which the
    array store adheres.

Other keys MUST NOT be present within the metadata object.

The members of a group are arrays and groups stored under logical paths that
are direct children of the parent group's logical path. E.g., if groups exist
under the logical paths "foo" and "foo/bar" and an array exists at logical path
"foo/baz" then the members of the group at path "foo" are the group at path
"foo/bar" and the array at path "foo/baz".

.. _spec_v2_attrs:

Attributes
----------

An array or group can be associated with custom attributes, which are arbitrary
key/value pairs with application-specific meaning. Custom attributes are encoded
as a JSON object and stored under the ".zattrs" key within an array store. The
".zattrs" key does not have to be present, and if it is absent the attributes
should be treated as empty.

For example, the JSON object below encodes three attributes named
"foo", "bar" and "baz"::

    {
        "foo": 42,
        "bar": "apples",
        "baz": [1, 2, 3, 4]
    }

.. _spec_v2_examples:

Examples
--------

Storing a single array
~~~~~~~~~~~~~~~~~~~~~~

Below is an example of storing a Zarr array, using a directory on the
local file system as storage.

Create an array::

    >>> import zarr
    >>> store = zarr.DirectoryStore('data/example.zarr')
    >>> a = zarr.create(shape=(20, 20), chunks=(10, 10), dtype='i4',
    ...                 fill_value=42, compressor=zarr.Zlib(level=1),
    ...                 store=store, overwrite=True)

No chunks are initialized yet, so only the ".zarray" and ".zattrs" keys
have been set in the store::

    >>> import os
    >>> sorted(os.listdir('data/example.zarr'))
    ['.zarray']

Inspect the array metadata::

    >>> print(open('data/example.zarr/.zarray').read())
    {
        "chunks": [
            10,
            10
        ],
        "compressor": {
            "id": "zlib",
            "level": 1
        },
        "dtype": "<i4",
        "fill_value": 42,
        "filters": null,
        "order": "C",
        "shape": [
            20,
            20
        ],
        "zarr_format": 2
    }

Chunks are initialized on demand. E.g., set some data::

    >>> a[0:10, 0:10] = 1
    >>> sorted(os.listdir('data/example.zarr'))
    ['.zarray', '0.0']

Set some more data::

    >>> a[0:10, 10:20] = 2
    >>> a[10:20, :] = 3
    >>> sorted(os.listdir('data/example.zarr'))
    ['.zarray', '0.0', '0.1', '1.0', '1.1']

Manually decompress a single chunk for illustration::

    >>> import zlib
    >>> buf = zlib.decompress(open('data/example.zarr/0.0', 'rb').read())
    >>> import numpy as np
    >>> chunk = np.frombuffer(buf, dtype='<i4')
    >>> chunk
    array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
           1, 1, 1, 1, 1, 1, 1, 1], dtype=int32)

Modify the array attributes::

    >>> a.attrs['foo'] = 42
    >>> a.attrs['bar'] = 'apples'
    >>> a.attrs['baz'] = [1, 2, 3, 4]
    >>> sorted(os.listdir('data/example.zarr'))
    ['.zarray', '.zattrs', '0.0', '0.1', '1.0', '1.1']
    >>> print(open('data/example.zarr/.zattrs').read())
    {
        "bar": "apples",
        "baz": [
            1,
            2,
            3,
            4
        ],
        "foo": 42
    }

Storing multiple arrays in a hierarchy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below is an example of storing multiple Zarr arrays organized into a group
hierarchy, using a directory on the local file system as storage. This storage
implementation maps logical paths onto directory paths on the file system,
however this is an implementation choice and is not required.

Setup the store::

    >>> import zarr
    >>> store = zarr.DirectoryStore('data/group.zarr')

Create the root group::

    >>> root_grp = zarr.group(store, overwrite=True)

The metadata resource for the root group has been created::

    >>> import os
    >>> sorted(os.listdir('data/group.zarr'))
    ['.zgroup']

Inspect the group metadata::

    >>> print(open('data/group.zarr/.zgroup').read())
    {
        "zarr_format": 2
    }

Create a sub-group::

    >>> sub_grp = root_grp.create_group('foo')

What has been stored::

    >>> sorted(os.listdir('data/group.zarr'))
    ['.zgroup', 'foo']
    >>> sorted(os.listdir('data/group.zarr/foo'))
    ['.zgroup']

Create an array within the sub-group::

    >>> a = sub_grp.create_dataset('bar', shape=(20, 20), chunks=(10, 10))
    >>> a[:] = 42

Set a custom attributes::

    >>> a.attrs['comment'] = 'answer to life, the universe and everything'

What has been stored::

    >>> sorted(os.listdir('data/group.zarr'))
    ['.zgroup', 'foo']
    >>> sorted(os.listdir('data/group.zarr/foo'))
    ['.zgroup', 'bar']
    >>> sorted(os.listdir('data/group.zarr/foo/bar'))
    ['.zarray', '.zattrs', '0.0', '0.1', '1.0', '1.1']

Here is the same example using a Zip file as storage::

    >>> store = zarr.ZipStore('data/group.zip', mode='w')
    >>> root_grp = zarr.group(store)
    >>> sub_grp = root_grp.create_group('foo')
    >>> a = sub_grp.create_dataset('bar', shape=(20, 20), chunks=(10, 10))
    >>> a[:] = 42
    >>> a.attrs['comment'] = 'answer to life, the universe and everything'
    >>> store.close()

What has been stored::

    >>> import zipfile
    >>> zf = zipfile.ZipFile('data/group.zip', mode='r')
    >>> for name in sorted(zf.namelist()):
    ...     print(name)
    .zgroup
    foo/.zgroup
    foo/bar/.zarray
    foo/bar/.zattrs
    foo/bar/0.0
    foo/bar/0.1
    foo/bar/1.0
    foo/bar/1.1

.. _spec_v2_changes:

Changes
-------

Version 2 clarifications
~~~~~~~~~~~~~~~~~~~~~~~~

The following changes have been made to the version 2 specification since it was
initially published to clarify ambiguities and add some missing information.

* The specification now describes how bytes fill values should be encoded and
  decoded for arrays with a fixed-length byte string data type (:issue:`165`,
  :issue:`176`).

* The specification now clarifies that units must be specified for datetime64 and
  timedelta64 data types (:issue:`85`, :issue:`215`).

* The specification now clarifies that the '.zattrs' key does not have to be present for
  either arrays or groups, and if absent then custom attributes should be treated as
  empty.

* The specification now describes how structured datatypes with
  subarray shapes and/or with nested structured data types are encoded
  in array metadata (:issue:`111`, :issue:`296`).

* Clarified the key/value pairs of custom attributes as "arbitrary" rather than
  "simple".

Changes from version 1 to version 2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following changes were made between version 1 and version 2 of this specification:

* Added support for storing multiple arrays in the same store and organising
  arrays into hierarchies using groups.
* Array metadata is now stored under the ".zarray" key instead of the "meta"
  key.
* Custom attributes are now stored under the ".zattrs" key instead of the
  "attrs" key.
* Added support for filters.
* Changed encoding of "fill_value" field within array metadata.
* Changed encoding of compressor information within array metadata to be
  consistent with representation of filter information.
