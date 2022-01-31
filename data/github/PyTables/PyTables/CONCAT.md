# Contribute to PyTables

PyTables is actively seeking additional maintainers and developers. If you are interested in becoming on, check out our developer resources on the PyTables website: http://www.pytables.org/development.html#pytables-development

## Documentation

In general, the barrier to entry is lower when making a docs contribution. We encourage you to check out our open docs issues: https://github.com/PyTables/PyTables/labels/documentation

## Stay up to date

If you want to stay up to date with PyTables development and communicate with the development team, consider joining our developer mailing list: https://groups.google.com/g/pytables-dev

## Something else? 

Do you have an idea for a contribution but you are not sure if it fits into any of the aforementioned options? Please open an issue and propose it! We are happy to help support contributions of all types.
=======================================
 Release notes for PyTables 3.7 series
=======================================

:Author: PyTables Developers
:Contact: pytables-dev@googlegroups.com

.. py:currentmodule:: tables


Changes from 3.7.0 to 3.7.1
===========================

* TBW


Changes from 3.6.1 to 3.7.0
===========================

Improvements
------------
- Compatibility with Python 3.10, numpy 1.21 and HDF5 1.12.
- Support for Python 3.5 has been dropped (:issue:`840` and :issue:`850`).
- Windows: Significantly faster `import tables` PR #781.
  Thanks to Christoph Gohlke.
- Internal C-Blosc sources updated to 1.21.1 (:issue:`931`).
  Note that, starting from C-Blosc 1.19 does not include the Snappy codec
  sources anymore, so Snappy will be not available if you compile from
  included sources; other packages (like conda or wheels),
  may (or may not) include it.
- Stop using appveyor and deprecated ci-helpers (closes :issue:`827`).
- Switch to `git submodule` for the management of vendored c-blosc sources.
- CI moved to GitHub Actions (GHA).
- Drop Travis-CI.
- Improved code formatting and notation consistency (:issue:`873`,
  :issue:`868`, :issue:`865` thanks to Miroslav Šedivý).
- Improve the use of modern Python including :mod:`pathlib`, f-strings
  (:issue:`859`, :issue:`855`, :issue:`839` and :issue:`818`
  thanks to Miroslav Šedivý).
- Several improvements to wheels generation in CI
  (thanks to Andreas Motl @amotl and Matthias @xmatthias).
- Simplified management of version information.
- Drop dependency on the deprecated distutils.
- Modernize the setup script and add support for PEP517 (:issue:`907`).

Bugfixes
--------
- Fix `pkg-config` (`setup.py`) for Python 3.9 on Debian.
  Thanks to Marco Sulla PR #792.
- Fix ROFileNode fails to return the `fileno()` (:issue:`633`).
- Do not flush read only files (:issue:`915` thanks to @lrepiton).

Other changes
-------------
- Drop the deprecated `hdf5Version` and `File.open_count`.
- the :func:`get_tables_version` and :func:`get_hdf5_version` functions are
  now deprecated please use the coresponding :data:`tables.__version__` and
  :data:`tables.hdf5_version` instead.
===========================================
 PyTables: hierarchical datasets in Python
===========================================

.. image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/PyTables/PyTables
   :target: https://gitter.im/PyTables/PyTables

.. image:: https://github.com/PyTables/PyTables/workflows/CI/badge.svg
   :target: https://github.com/PyTables/PyTables/actions?query=workflow%3ACI

.. image:: https://img.shields.io/pypi/v/tables.svg
  :target: https://pypi.org/project/tables/

.. image:: https://img.shields.io/pypi/pyversions/tables.svg
  :target: https://pypi.org/project/tables/

.. image:: https://img.shields.io/pypi/l/tables
  :target: https://github.com/PyTables/PyTables/


:URL: http://www.pytables.org/


PyTables is a package for managing hierarchical datasets and designed
to efficiently cope with extremely large amounts of data.

It is built on top of the HDF5 library and the NumPy package. It
features an object-oriented interface that, combined with C extensions
for the performance-critical parts of the code (generated using
Cython), makes it a fast, yet extremely easy to use tool for
interactively save and retrieve very large amounts of data. One
important feature of PyTables is that it optimizes memory and disk
resources so that they take much less space (between a factor 3 to 5,
and more if the data is compressible) than other solutions, like for
example, relational or object oriented databases.

State-of-the-art compression
----------------------------

PyTables comes with out-of-box support for the `Blosc compressor
<http://www.blosc.org>`_.  This allows for extremely high compression
speed, while keeping decent compression ratios.  By doing so, I/O can
be accelerated by a large extent, and you may end achieving higher
performance than the bandwidth provided by your I/O subsystem.  See
the `Tuning The Chunksize section of the Optimization Tips chapter
<http://www.pytables.org/usersguide/optimization.html#fine-tuning-the-chunksize>`_
of user documentation for some benchmarks.

Not a RDBMS replacement
-----------------------

PyTables is not designed to work as a relational database replacement,
but rather as a teammate. If you want to work with large datasets of
multidimensional data (for example, for multidimensional analysis), or
just provide a categorized structure for some portions of your
cluttered RDBS, then give PyTables a try. It works well for storing
data from data acquisition systems (DAS), simulation software, network
data monitoring systems (for example, traffic measurements of IP
packets on routers), or as a centralized repository for system logs,
to name only a few possible uses.

Tables
------

A table is defined as a collection of records whose values are stored
in fixed-length fields. All records have the same structure and all
values in each field have the same data type. The terms "fixed-length"
and strict "data types" seems to be quite a strange requirement for an
interpreted language like Python, but they serve a useful function if
the goal is to save very large quantities of data (such as is
generated by many scientific applications, for example) in an
efficient manner that reduces demand on CPU time and I/O.

Arrays
------

There are other useful objects like arrays, enlargeable arrays or
variable length arrays that can cope with different missions on your
project.

Easy to use
-----------

One of the principal objectives of PyTables is to be user-friendly.
In addition, many different iterators have been implemented so as to
enable the interactive work to be as productive as possible.

Platforms
---------

We are using Linux on top of Intel32 and Intel64 boxes as the main
development platforms, but PyTables should be easy to compile/install
on other UNIX or Windows machines.

Compiling
---------

To compile PyTables you will need, at least, a recent version of HDF5
(C flavor) library, the Zlib compression library and the NumPy and
Numexpr packages. Besides, it comes with support for the Blosc, LZO
and bzip2 compressor libraries. Blosc is mandatory, but PyTables comes
with Blosc sources so, although it is recommended to have Blosc
installed in your system, you don't absolutely need to install it
separately.  LZO and bzip2 compression libraries are, however,
optional.

Installation
------------

1. Make sure you have HDF5 version 1.8.4 or above.

   On OSX you can install HDF5 using `Homebrew <http://brew.sh>`_::

       $ brew install hdf5

   On debian bases distributions::

       $ sudo apt-get install libhdf5-serial-dev

   If you have the HDF5 library in some non-standard location (that
   is, where the compiler and the linker can't find it) you can use
   the environment variable `HDF5_DIR` to specify its location. See
   `the manual
   <http://www.pytables.org/usersguide/installation.html>`_ for more
   details.

3. For stability (and performance too) reasons, it is strongly
   recommended that you install the C-Blosc library separately,
   although you might want PyTables to use its internal C-Blosc
   sources.

3. Optionally, consider to install the LZO compression library and/or
   the bzip2 compression library.

4. Install!::

       $ python3 -m pip install tables

5. To run the test suite run::

       $ python3 -m tables.tests.test_all

   If there is some test that does not pass, please send the
   complete output for tests back to us.


**Enjoy data!** -- The PyTables Team

.. Local Variables:
.. mode: text
.. coding: utf-8
.. fill-column: 70
.. End:
=====================
Blosc filter for HDF5
=====================

:Travis CI: |travis|
:And...: |powered|

.. |travis| image:: https://travis-ci.org/Blosc/hdf5.png?branch=master
        :target: https://travis-ci.org/Blosc/hdf5

.. |powered| image:: http://b.repl.ca/v1/Powered--By-Blosc-blue.png
        :target: https://blosc.org

This is an example of filter for HDF5 that uses the Blosc compressor.

You need to be a bit careful before using this filter because you
should not activate the shuffle right in HDF5, but rather from Blosc
itself.  This is because Blosc uses an SIMD shuffle internally which
is much faster.


Using the Blosc filter from HDF5
================================

In order to register Blosc into your HDF5 application, you only need
to call a function in blosc_filter.h, with the following signature:

    int register_blosc(char **version, char **date)

Calling this will register the filter with the HDF5 library and will
return info about the Blosc release in `**version` and `**date`
char pointers.

A non-negative return value indicates success.  If the registration
fails, an error is pushed onto the current error stack and a negative
value is returned.

An example C program ('src/example.c') is included which demonstrates
the proper use of the filter.

This filter has been tested against HDF5 versions 1.6.5 through
1.8.10.  It is released under the MIT license (see LICENSE.txt for
details).


Compiling
=========

The filter consists of a single 'src/blosc_filter.c' source file and
'src/blosc_filter.h' header, which will need the Blosc library
installed to work.


As an HDF5 plugin
=================

Also, you can use blosc as an HDF5 plugin; see 'src/blosc_plugin.c' for
details.


Acknowledgments
===============

See THANKS.rst.


----

  **Enjoy data!**
====================
PyTables Development
====================

If you want to follow the development of PyTables and take part in it,
you may have a look at the PyTables project pages on
`GitHub <https://github.com>`_.

The source code for PyTables may be found at the `GitHub project site`_.
You can get a copy of the latest version of the source code (under
development) from the master branch of the project repository using git::

    git clone --recursive git@github.com:PyTables/PyTables.git

Also, be sure to subscribe to the `Users' Mailing List`_ and/or the
`Developers' Mailing List`_.

.. _`GitHub project site`: https://github.com/PyTables
.. _`Users' Mailing List`: https://groups.google.com/group/pytables-users
.. _`Developers' Mailing List`: https://groups.google.com/group/pytables-dev

Other resources for developers:

* `GitHub project site`_
* :ref:`library_reference`
* `Git Repository browser <https://github.com/PyTables/PyTables>`_
* `Issue tracker <https://github.com/PyTables/PyTables/issues>`_
* `Developers wiki <https://github.com/PyTables/PyTables/wiki>`_
* `Users' Mailing List`_
* `Developers' Mailing List`_
* Continuous Integration:

  - `GitHub Actions (GHA) <https://github.com/PyTables/PyTables/actions>`_

.. todo:: improve this section
======================
PyTables Release Notes
======================

Migration
----------

.. toctree::
    :maxdepth: 1

    Migrating from 2.x to 3.x <MIGRATING_TO_3.x>
    Migrating from 1.x to 2.x <MIGRATING_TO_2.x>

PyTables
--------

.. toctree::
    :maxdepth: 1

    release-notes/RELEASE_NOTES_v3.7.x
    release-notes/RELEASE_NOTES_v3.6.x
    release-notes/RELEASE_NOTES_v3.5.x
    release-notes/RELEASE_NOTES_v3.4.x
    release-notes/RELEASE_NOTES_v3.3.x
    release-notes/RELEASE_NOTES_v3.2.x
    release-notes/RELEASE_NOTES_v3.1.x
    release-notes/RELEASE_NOTES_v3.0.x
    release-notes/RELEASE_NOTES_v2.4.x
    release-notes/RELEASE_NOTES_v2.3.x
    release-notes/RELEASE_NOTES_v2.2.x
    release-notes/RELEASE_NOTES_v2.1.x
    release-notes/RELEASE_NOTES_v2.0.x
    release-notes/RELEASE_NOTES_v1.4
    release-notes/RELEASE_NOTES_v1.3.3
    release-notes/RELEASE_NOTES_v1.3.2
    release-notes/RELEASE_NOTES_v1.3.1
    release-notes/RELEASE_NOTES_v1.3
    release-notes/RELEASE_NOTES_v1.2.3
    release-notes/RELEASE_NOTES_v1.2.2
    release-notes/RELEASE_NOTES_v1.2.1
    release-notes/RELEASE_NOTES_v1.2
    release-notes/RELEASE_NOTES_v1.1.1
    release-notes/RELEASE_NOTES_v1.1
    release-notes/RELEASE_NOTES_v1.0
    release-notes/RELEASE_NOTES_v0.9.1
    release-notes/RELEASE_NOTES_v0.9
    release-notes/RELEASE_NOTES_v0.8
    release-notes/RELEASE_NOTES_v0.7.2
    release-notes/RELEASE_NOTES_v0.7.1


PyTables Pro
------------

.. toctree::
    :maxdepth: 1

    release-notes/RELEASE_NOTES_v2.2.x-pro
    release-notes/RELEASE_NOTES_v2.1.x-pro
    release-notes/RELEASE_NOTES_v2.0.x-pro


Release timeline
----------------

=============== =========== ==========
PyTables        3.5.2       2019-05-31
PyTables        3.5.1       2019-03-14
PyTables        3.5.0       2019-03-13
PyTables        3.4.4       2018-06-11
PyTables        3.4.3       2018-04-17
PyTables        3.4.2       2017-04-19
PyTables        3.4.1       2017-04-12
PyTables        3.4.0       2017-04-11
PyTables        3.3.0       2016-09-12
PyTables        3.2.3.1     2016-07-05
PyTables        3.2.3       2016-07-04
PyTables        3.2.2       2015-09-22
PyTables        3.2.1.1     2015-08-31
PyTables        3.2.1       2015-08-04
PyTables        3.2.0       2015-05-06
PyTables        3.2.0rc2    2015-05-01
PyTables        3.2.0rc1    2015-04-21
PyTables        3.1.1       2014-03-25
PyTables        3.1.0       2014-02-05
PyTables        3.1.0rc2    2014-01-22
PyTables        3.1.0rc1    2014-01-17
PyTables        3.0         2013-06-01
PyTables        3.0rc3      2013-05-29
PyTables        3.0rc2      2013-05-17
PyTables        3.0rc1      2013-05-10
PyTables        3.0b1       2013-04-27
PyTables        2.4         2012-07-20
PyTables        2.4rc1      2012-07-16
PyTables        2.4b1       2012-07-07
PyTables        2.3.1       2011-10-28
PyTables        2.3         2011-09-23
PyTables        2.3rc1      2011-09-11
PyTables        2.2.1       2010-11-05
PyTables        2.2.1rc1    2010-11-03
Pytables        2.2         2010-07-01
PyTables        2.2rc2      2010-06-17
PyTables        2.2rc1      2010-05-20
PyTables        2.1.2       2009-09-14
PyTables        2.1.1       2009-03-13
PyTables Pro    2.1.1       2009-03-13
PyTables        2.1         2008-12-19
PyTables        2.1rc2      2008-11-18
PyTables        2.1rc1      2008-10-31
PyTables        2.0.4       2008-07-05
PyTables Pro    2.0.4       2008-07-05
PyTables        2.0.3       2008-03-07
PyTables Pro    2.0.2.1     2007-12-24
PyTables Pro    2.0.1       2007-09-20
PyTables        2.0.1       2007-09-20
PyTables Pro    2.0         2007-07-12
PyTables        2.0         2007-07-12
PyTables        2.0rc2      2007-05-28
PyTables        2.0rc1      2007-04-26
PyTables        1.4         2006-12-21
PyTables        1.3.3       2006-08-24
PyTables        1.3.2       2006-06-20
PyTables        1.3.1       2006-05-02
PyTables        1.3         2006-04-01
PyTables        1.2.3       2006-02-23
PyTables        1.2.2       2006-02-16
PyTables        1.2.1       2005-12-21
PyTables        1.2         2005-11-22
PyTables        1.1.1       2005-09-13
PyTables        1.1         2005-07-14
PyTables        1.0         2005-05-12
PyTables        0.9.1       2004-12-02
PyTables-       0.9         2004-11-08
PyTables        0.8.1       2004-07-13
PyTables        0.8         2004-03-03
PyTables        0.7.2       2003-09-22
PyTables        0.7         2003-07-31
PyTables        0.5.1       2003-05-14
PyTables        0.5         2003-05-10
Pytables        0.4         2003-03-19
=============== =========== ==========
:author: FrancescAlted
:date: 2011-06-13 08:40:20

.. py:currentmodule:: tables

===
FAQ
===

General questions
=================

What is PyTables?
-----------------

PyTables is a package for managing hierarchical datasets designed to
efficiently cope with extremely large amounts of data.

It is built on top of the HDF5_  library, the `Python language`_ and the
NumPy_ package.
It features an object-oriented interface that, combined with C extensions
for the performance-critical parts of the code, makes it a fast yet
extremely easy-to-use tool for interactively storing and retrieving very
large amounts of data.


What are PyTables' licensing terms?
-----------------------------------

PyTables is free for both commercial and non-commercial use, under the terms
of the `BSD 3-Clause License <http://opensource.org/licenses/BSD-3-Clause>`_.


I'm having problems. How can I get support?
-------------------------------------------

The most common and efficient way is to subscribe (remember you *need* to
subscribe prior to send messages) to the PyTables `users mailing list`_, and
send there a brief description of your issue and, if possible, a short script
that can reproduce it.
Hopefully, someone on the list will be able to help you.
It is also a good idea to check out the `archives of the user's list`_ (you may
want to check the `Gmane archives`_ instead) so as to see if the answer to your
question has already been dealt with.


Why HDF5?
---------

HDF5_ is the underlying C library and file format that enables PyTables to
efficiently deal with the data.  It has been chosen for the following reasons:

* Designed to efficiently manage very large datasets.
* Lets you organize datasets hierarchically.
* Very flexible and well tested in scientific environments.
* Good maintenance and improvement rate.
* Technical excellence (`R&D 100 Award`_).
* **It's Open Source software**


Why Python?
-----------

1. Python is interactive.

   People familiar with data processing understand how powerful command line
   interfaces are for exploring mathematical relationships and scientific data
   sets.  Python provides an interactive environment with the added benefit of
   a full featured programming language behind it.

2. Python is productive for beginners and experts alike.

   PyTables is targeted at engineers, scientists, system analysts, financial
   analysts, and others who consider programming a necessary evil.  Any time
   spent learning a language or tracking down bugs is time spent not solving
   their real problem.  Python has a short learning curve and most people can
   do real and useful work with it in a day of learning.  Its clean syntax and
   interactive nature facilitate this.

3. Python is data-handling friendly.

   Python comes with nice idioms that make the access to data much easier:
   general slicing (i.e. ``data[start:stop:step]``), list comprehensions,
   iterators, generators ... are constructs that make the interaction with your
   data very easy.


Why NumPy?
----------

NumPy_ is a Python package to efficiently deal with large datasets
**in-memory**, providing containers for homogeneous data, heterogeneous data,
and string arrays.
PyTables uses these NumPy containers as *in-memory buffers* to push the I/O
bandwith towards the platform limits.


Where can PyTables be applied?
==============================

In all the scenarios where one needs to deal with large datasets:

* Industrial applications

  - Data acquisition in real time
  - Quality control
  - Fast data processing

* Scientific applications

  - Meteorology, oceanography
  - Numerical simulations
  - Medicine (biological sensors, general data gathering & processing)

* Information systems

  - System log monitoring & consolidation
  - Tracing of routing data
  - Alert systems in security


Is PyTables safe?
-----------------

Well, first of all, let me state that PyTables does not support transactional
features yet (we don't even know if we will ever be motivated to implement
this!), so there is always the risk that you can lose your data in case of an
unexpected event while writing (like a power outage, system shutdowns ...).
Having said that, if your typical scenarios are *write once, read many*, then
the use of PyTables is perfectly safe, even for dealing extremely large amounts
of data.


Can PyTables be used in concurrent access scenarios?
----------------------------------------------------

It depends. Concurrent reads are no problem at all. However, whenever a process
(or thread) is trying to write, then problems will start to appear.  First,
PyTables doesn't support locking at any level, so several process writing
concurrently to the same PyTables file will probably end up corrupting it, so
don't do this!  Even having only one process writing and the others reading is
a hairy thing, because the reading processes might be reading incomplete data
from a concurrent data writing operation.

The solution would be to lock the file while writing and unlock it after a
flush over the file has been performed.  Also, in order to avoid cache (HDF5_,
PyTables) problems with read apps, you would need to re-open your files
whenever you are going to issue a read operation.  If a re-opening operation is
unacceptable in terms of speed, you may want to do all your I/O operations in
one single process (or thread) and communicate the results via sockets,
:class:`Queue.Queue` objects (in case of using threads), or whatever, with the
client process/thread.

The `examples` directory contains two scripts demonstrating methods of
accessing a PyTables file from multiple processes.

The first, *multiprocess_access_queues.py*, uses a
:class:`multiprocessing.Queue` object to transfer read and write requests from
multiple *DataProcessor* processes to a single process responsible for all
access to the PyTables file.  The results of read requests are then transferred
back to the originating processes using other :class:`Queue` objects.

The second example script, *multiprocess_access_benchmarks.py*, demonstrates
and benchmarks four methods of transferring PyTables array data between
processes.  The four methods are:

 * Using :class:`multiprocessing.Pipe` from the Python standard library.
 * Using a memory mapped file that is shared between two processes.  The NumPy
   array associated with the file is passed as the *out* argument to the
   :meth:`tables.Array.read` method.
 * Using a Unix domain socket.  Note that this example uses the 'abstract
   namespace' and will only work under Linux.
 * Using an IPv4 socket.

See also the discussion in :issue:`790`.


What kind of containers does PyTables implement?
------------------------------------------------

PyTables does support a series of data containers that address specific needs
of the user. Below is a brief description of them:

::class:`Table`:
    Lets you deal with heterogeneous datasets. Allows compression. Enlargeable.
    Supports nested types. Good performance for read/writing data.
::class:`Array`:
    Provides quick and dirty array handling. Not compression allowed.
    Not enlargeable. Can be used only with relatively small datasets (i.e.
    those that fit in memory). It provides the fastest I/O speed.
::class:`CArray`:
    Provides compressed array support. Not enlargeable. Good speed when
    reading/writing.
::class:`EArray`:
    Most general array support. Compressible and enlargeable. It is pretty
    fast at extending, and very good at reading.
::class:`VLArray`:
    Supports collections of homogeneous data with a variable number of entries.
    Compressible and enlargeable. I/O is not very fast.
::class:`Group`:
    The structural component.
    A hierarchically-addressable container for HDF5 nodes (each of these
    containers, including Group, are nodes), similar to a directory in a
    UNIX filesystem.

Please refer to the  :doc:`usersguide/libref` for more specific information.


Cool! I'd like to see some examples of use.
-------------------------------------------

Sure. Go to the HowToUse section to find simple examples that will help you
getting started.


Can you show me some screenshots?
---------------------------------

Well, PyTables is not a graphical library by itself.  However, you may want to
check out ViTables_, a GUI tool to browse and edit PyTables & HDF5_ files.


Is PyTables a replacement for a relational database?
----------------------------------------------------

No, by no means. PyTables lacks many features that are standard in most
relational databases.  In particular, it does not have support for
relationships (beyond the hierarchical one, of course) between datasets and it
does not have transactional features.  PyTables is more focused on speed and
dealing with really large datasets, than implementing the above features.  In
that sense, PyTables can be best viewed as a *teammate* of a relational
database.

For example, if you have very large tables in your existing relational
database, they will take lots of space on disk, potentially reducing the
performance of the relational engine.  In such a case, you can move those huge
tables out of your existing relational database to PyTables, and let your
relational engine do what it does best (i.e.  manage relatively small or medium
datasets with potentially complex relationships), and use PyTables for what it
has been designed for (i.e. manage large amounts of data which are loosely
related).


How can PyTables be fast if it is written in an interpreted language like Python?
---------------------------------------------------------------------------------

Actually, all of the critical I/O code in PyTables is a thin layer of code on
top of HDF5_, which is a very efficient C library. Cython_ is used as the
*glue* language to generate "wrappers" around HDF5 calls so that they can be
used in Python.  Also, the use of an efficient numerical package such as NumPy_
makes the most costly operations effectively run at C speed.  Finally,
time-critical loops are usually implemented in Cython_ (which, if used
properly, allows to generate code that runs at almost pure C speeds).


If it is designed to deal with very large datasets, then PyTables should consume a lot of memory, shouldn't it?
---------------------------------------------------------------------------------------------------------------

Well, you already know that PyTables sits on top of HDF5, Python and NumPy_,
and if we add its own logic (~7500 lines of code in Python, ~3000 in Cython and
~4000 in C), then we should conclude that PyTables isn't effectively a paradigm
of lightness.

Having said that, PyTables (as HDF5_ itself) tries very hard to optimize the
memory consumption by implementing a series of features like dynamic
determination of buffer sizes, *Least Recently Used* cache for keeping unused
nodes out of memory, and extensive use of compact NumPy_ data containers.
Moreover, PyTables is in a relatively mature state and most memory leaks have
been already addressed and fixed.

Just to give you an idea of what you can expect, a PyTables program can deal
with a table with around 30 columns and 1 million entries using as low as 13 MB
of memory (on a 32-bit platform).  All in all, it is not that much, is it?.


Why was PyTables born?
----------------------

Because, back in August 2002, one of its authors (`Francesc Alted`_) had a need
to save lots of hierarchical data in an efficient way for later post-processing
it.  After trying out several approaches, he found that they presented distinct
inconveniences.  For example, working with file sizes larger than, say, 100 MB,
was rather painful with ZODB (it took lots of memory with the version available
by that time).

The netCDF3_ interface provided by `Scientific Python`_ was great, but it did
not allow to structure the hierarchically; besides, netCDF3_ only supports
homogeneous datasets, not heterogeneous ones (i.e. tables). (As an aside,
netCDF4_ overcomes many of the limitations of netCDF3_, although curiously
enough, it is based on top of HDF5_, the library chosen as the base for
PyTables from the very beginning.)

So, he decided to give HDF5_ a try, start doing his own wrappings to it and
voilà, this is how the first public release of PyTables (0.1) saw the light in
October 2002, three months after his itch started to eat him ;-).


How does PyTables compare with the h5py project?
------------------------------------------------

Well, they are similar in that both packages are Python interfaces to the HDF5_
library, but there are some important differences to be noted.  h5py_ is an
attempt to map the HDF5_ feature set to NumPy_ as closely as possible.  In
addition, it also provides access to nearly all of the HDF5_ C API.

Instead, PyTables builds up an additional abstraction layer on top of HDF5_ and
NumPy_ where it implements things like an enhanced type system, an :ref:`engine
for enabling complex queries <searchOptim>`, an efficient computational
kernel, advanced indexing capabilities or an undo/redo feature, to name
just a few.  This additional layer also allows PyTables to be relatively
independent of its underlying libraries (and their possible limitations).  For
example, PyTables can support HDF5_ data types like `enumerated` or `time` that
are available in the HDF5_ library but not in the NumPy_ package; or even
perform powerful complex queries that are not implemented directly in neither
HDF5_ nor NumPy_.

Furthermore, PyTables also tries hard to be a high performance interface to
HDF5/NumPy, implementing niceties like internal LRU caches for nodes and other
data and metadata, :ref:`automatic computation of optimal chunk sizes
<chunksizeFineTune>` for the datasets, a variety of compressors, ranging from
slow but efficient (bzip2_) to extremely fast ones (Blosc_) in addition to the
standard `zlib`_.  Another difference is that PyTables makes use of numexpr_ so
as to accelerate internal computations (for example, in evaluating complex
queries) to a maximum.

For contrasting with other opinions, you may want to check the PyTables/h5py
comparison in a similar entry of the `FAQ of h5py`_.


I've found a bug.  What do I do?
--------------------------------

The PyTables development team works hard to make this eventuality as rare as
possible, but, as in any software made by human beings, bugs do occur.  If you
find any bug, please tell us by file a bug report in the `issue tracker`_ on
GitHub_.


Is it possible to get involved in PyTables development?
-------------------------------------------------------

Indeed. We are keen for more people to help out contributing code, unit tests,
documentation, and helping out maintaining this wiki. Drop us a mail on the
`users mailing list` and tell us in which area do you want to work.


How can I cite PyTables?
------------------------

The recommended way to cite PyTables in a paper or a presentation is as
following:

* Author: Francesc Alted, Ivan Vilata and others
* Title: PyTables: Hierarchical Datasets in Python
* Year: 2002 -
* URL: http://www.pytables.org

Here's an example of a BibTeX entry::

    @Misc{,
      author =    {PyTables Developers Team},
      title =     {{PyTables}: Hierarchical Datasets in {Python}},
      year =      {2002--},
      url = "http://www.pytables.org/"
    }


PyTables 2.x issues
===================

I'm having problems migrating my apps from PyTables 1.x into PyTables 2.x. Please, help!
----------------------------------------------------------------------------------------

Sure.  However, you should first check out the :doc:`MIGRATING_TO_2.x`
document.
It should provide hints to the most frequently asked questions on this regard.


For combined searches like `table.where('(x<5) & (x>3)')`, why was a `&` operator chosen instead of an `and`?
-------------------------------------------------------------------------------------------------------------

Search expressions are in fact Python expressions written as strings, and they
are evaluated as such.  This has the advantage of not having to learn a new
syntax, but it also implies some limitations with logical `and` and `or`
operators, namely that they can not be overloaded in Python.  Thus, it is
impossible right now to get an element-wise operation out of an expression like
`'array1 and array2'`.  That's why one has to choose some other operator, being
`&` and `|` the most similar to their C counterparts `&&` and `||`, which
aren't available in Python either.

You should be careful about expressions like `'x<5 & x>3'` and others like `'3
< x < 5'` which ''won't work as expected'', because of the different operator
precedence and the absence of an overloaded logical `and` operator.  More on
this in the appendix about condition syntax in the `HDF5 manual`_.

There are quite a few packages affected by those limitations including NumPy_
themselves and SQLObject_, and there have been quite longish discussions about
adding the possibility of overloading logical operators to Python (see `PEP
335`_ and `this thread`__ for more details).

__ https://mail.python.org/pipermail/python-dev/2004-September/048763.html


I can not select rows using in-kernel queries with a condition that involves an UInt64Col. Why?
-----------------------------------------------------------------------------------------------

This turns out to be a limitation of the numexpr_ package.  Internally,
numexpr_ uses a limited set of types for doing calculations, and unsigned
integers are always upcasted to the immediate signed integer that can fit the
information.  The problem here is that there is not a (standard) signed integer
that can be used to keep the information of a 64-bit unsigned integer.

So, your best bet right now is to avoid `uint64` types if you can.  If you
absolutely need `uint64`, the only way for doing selections with this is
through regular Python selections.  For example, if your table has a `colM`
column which is declared as an `UInt64Col`, then you can still filter its
values with::

    [row['colN'] for row in table if row['colM'] < X]


However, this approach will generally lead to slow speed (specially on Win32
platforms, where the values will be converted to Python `long` values).


I'm already using PyTables 2.x but I'm still getting numarray objects instead of NumPy ones!
--------------------------------------------------------------------------------------------

This is most probably due to the fact that you are using a file created with
PyTables 1.x series.  By default, PyTables 1.x was setting an HDF5 attribute
`FLAVOR` with the value `'numarray'` to all leaves.  Now, PyTables 2.x sees
this attribute and obediently converts the internal object (truly a NumPy
object) into a `numarray` one.  For PyTables 2.x files the `FLAVOR` attribute
will only be saved when explicitly set via the `leaf.flavor` property (or when
passing data to an :class:`Array` or :class:`Table` at creation time), so you
will be able to distinguish default flavors from user-set ones by checking the
existence of the `FLAVOR` attribute.

Meanwhile, if you don't want to receive `numarray` objects when reading old
files, you have several possibilities:

* Remove the flavor for your datasets by hand::

     for leaf in h5file.walkNodes(classname='Leaf'):
         del leaf.flavor

* Use the :program:'ptrepack` utility with the flag `--upgrade-flavors`
  so as to convert all flavors in old files to the default (effectively by
  removing the `FLAVOR` attribute).
* Remove the `numarray` (and/or `Numeric`) package from your system.
  Then PyTables 2.x will return you pure NumPy objects (it can't be
  otherwise!).


Installation issues
===================

Windows
-------

Error when importing tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You have installed the binary installer for Windows and, when importing the
*tables* package you are getting an error like::

    The command in "0x6714a822" refers to memory in "0x012011a0". The
    procedure "written" could not be executed.
    Click to ok to terminate.
    Click to abort to debug the program.

This problem can be due to a series of reasons, but the most probable one is
that you have a version of a DLL library that is needed by PyTables and it is
not at the correct version.  Please, double-check the versions of the required
libraries for PyTables and install newer versions, if needed. In most cases,
this solves the issue.

In case you continue getting problems, there are situations where other
programs do install libraries in the PATH that are **optional** to PyTables
(for example BZIP2 or LZO), but that they will be used if they are found in
your system (i.e. anywhere in your :envvar:`PATH`).  So, if you find any of
these libraries in your PATH, upgrade it to the latest version available (you
don't need to re-install PyTables).


-----


.. target-notes::

.. _HDF5: http://www.hdfgroup.org/HDF5
.. _`Python language`: http://www.python.org
.. _NumPy: http://www.numpy.org
.. _`users mailing list`: https://groups.google.com/group/pytables-users
.. _`archives of the user's list`: https://sourceforge.net/p/pytables/mailman/pytables-users/
.. _`Gmane archives`: http://www.mail-archive.com/pytables-users@lists.sourceforge.net/
.. _`R&D 100 Award`: http://www.hdfgroup.org/HDF5/RD100-2002/
.. _ViTables: http://vitables.org
.. _Cython: http://www.cython.org
.. _`Francesc Alted`: https://github.com/FrancescAlted
.. _netCDF3: http://www.unidata.ucar.edu/software/netcdf
.. _`Scientific Python`: http://dirac.cnrs-orleans.fr/ScientificPython.html
.. _netCDF4: http://www.unidata.ucar.edu/software/netcdf
.. _OPeNDAP: http://opendap.org
.. _`PyTables Manual`: http://www.pytables.org/usersguide/index.html
.. _h5py: http://www.h5py.org
.. _bzip2: http://www.bzip.org
.. _Blosc: https://www.blosc.org
.. _`zlib`: http://zlib.net
.. _numexpr: https://github.com/pydata/numexpr
.. _`FAQ of h5py`: http://docs.h5py.org/en/latest/faq.html#what-s-the-difference-between-h5py-and-pytables
.. _`issue tracker`: https://github.com/PyTables/PyTables/issues
.. _GitHub: https://github.com
.. _`HDF5 manual`: https://portal.hdfgroup.org/display/HDF5/Datatypes
.. _SQLObject: http://sqlobject.org
.. _`PEP 335`: http://www.python.org/dev/peps/pep-0335


.. todo:: fix links that point to wiki pages

==============
Other Material
==============

Videos
======

These are the videos of a series dedicated to introduce the main features of
PyTables in a visual and easy to grasp manner.
More videos will be made available with the time:

* `HDF5 is for Lovers, SciPy 2012 Tutorial <http://www.youtube.com/watch?v=Nzx0HAd3FiI>`_:
  a beginner's introduction to PyTables and HDF5.


Presentations
=============

Here are the slides of some presentations about PyTables that you may find
useful:

* HDF5 is for Lovers, SciPy 2012 Tutorial, July 2012, Austin, TX, USA,
  `slides (pdf) <https://raw.github.com/scopatz/scipy2012/master/hdf5/scopatz_scipy2012_hdf5.pdf>`_,
  `video <http://www.youtube.com/watch?v=Nzx0HAd3FiI>`_,
  `exercises <https://github.com/scopatz/scipy2012/tree/master/hdf5/exer>`_,
  `solutions <https://github.com/scopatz/scipy2012/tree/master/hdf5/sol>`_, and
  `repository <https://github.com/scopatz/scipy2012>`_.
* `An on-disk binary data container <http://www.pytables.org/docs/PUG-Austin-2012-v3.pdf>`_.
  Talk given at the `Austin Python Meetup <http://www.meetup.com/austinpython>`_,
  Austin, TX, USA (May 2012).
* `Large Data Analysis with Python <http://www.pytables.org/docs/LargeDataAnalysis.pdf>`_.
  Seminar given at the `German Neuroinformatics Node <http://www.g-node.org>`_,
  Munich, Germany (November 2010).
* `Starving CPUs (and coping with that in PyTables)
  <http://www.pytables.org/docs/StarvingCPUs-PyTablesUsages.pdf>`_.
  Seminar given at `FOM Institute for Plasma Physics Rijnhuizen <http://www.rijnhuizen.nl/>`_,
  The Netherlands (September 2009).
* `On The Data Access Issue (or Why Modern CPUs Are Starving)
  <http://www.pytables.org/docs/StarvingCPUs.pdf>`_.
  Keynote presented at `EuroSciPy 2009 <https://www.euroscipy.org/>`_ conference
  in Leipzig, Germany (July 2009).
* `An Overview of Future Improvements to OPSI
  <http://www.pytables.org/docs/THG-2007-PlansForNewOPSI.pdf>`_.
  Informal talk given at the `THG headquarters <http://www.hdfgroup.org>`_ in
  Urbana-Champaign, Illinois, USA (October 2007).
* `Finding Needles in a Huge DataStack
  <http://www.pytables.org/docs/FindingNeedles.pdf>`_.
  Talk given at the **EuroPython 2006 Conference**, held at CERN, Genève,
  Switzerland (July 2006).
* `Presentation given at the "HDF Workshop 2005"
  <http://www.pytables.org/docs/HDF_IX_Workshop.pdf>`_, held at San Francisco,
  USA (December 2005).
* `I <http://www.pytables.org/docs/taller-sf1-color.pdf>`_ and
  `II <http://www.pytables.org/docs/taller-sf2-color.pdf>`_ **Workshop in Free
  Software and Scientific Computing** given at the Universitat Jaume I,
  Castelló, Spain (October 2004). In Catalan.
* `Presentation given at the "SciPy Workshop 2004"
  <http://www.pytables.org/docs/SciPy04.pdf>`_, held at Caltech, Pasadena,
  USA (September 2004).
* `Slides <http://www.pytables.org/docs/EuroPython2003.pdf>`_ of presentation
  given at **EuroPython Conference** in Charleroi, Belgium (June 2003).
* `Presentation for the "iParty5" <http://www.pytables.org/docs/iparty2003.pdf>`_
  held at Castelló, Spain (May 2003). In Spanish.
* `Talk <http://www.pytables.org/docs/pycon2003.pdf>`_ on PyTables given at
  the **PyCon 2003 Convention** held at Washington, USA (March 203).


Reports
=======

* White Paper on `OPSI indexes <http://www.pytables.org/docs/OPSI-indexes.pdf>`_,
  explaining the powerful new indexing engine in PyTables Pro.
* `Performance study <http://www.pytables.org/docs/NewObjectTreeCache.pdf>`_
  on how the new object tree cache introduced in PyTables 1.2 can accelerate
  the opening of files with a large number of objects, while being quite less
  memory hungry.
* `Paper version <http://www.pytables.org/docs/pycon2003-paper.pdf>`_ of the
  presentation at PyCon2003.



Other sources for examples
==========================

The examples presented above show just a little amount of the full capabilities
of PyTables.
Please check out the documentation and the :file:`examples/` directory in the
source package for more examples.

==================================
Migrating from PyTables 2.x to 3.x
==================================

:Author: Antonio Valentino
:Author: Anthony Scopatz
:Author: Thomas Provoost

This document describes the major changes in PyTables in going from the
2.x to 3.x series and what you need to know when migrating downstream
code bases.

Python 3 at Last!
=================

The PyTables 3.x series now ships with full compatibility for Python 3.1+.
Additionally, we plan on maintaining compatibility with Python 2.7 for the
foreseeable future.  Python 2.6 is no longer supported but
may work in most cases.  Note that the entire 3.x series now relies on
numexpr v2.1+, which itself is the first version of numexpr support both
Python 2 & 3.

Numeric, Numarray, NetCDF3, & HDF5 1.6 No More!
===============================================

PyTables no longer supports numeric and numarray. Please use numpy instead.
Additionally, the ``tables.netcdf3`` module has been removed. Please refer
to the `netcdf4-python`_ project for further support. Lastly, the older
HDF5 1.6 API is no longer supported.  Please upgrade to HDF5 1.8+.

Unicode all the strings!
========================

In Python 3, all strings are natively in Unicode. This introduces some 
difficulties, as the native HDF5 string format is not Unicode-compatible. 
To minimize explicit conversion troubles when writing, especially :ref:`when 
creating data sets from existing Python objects <create-signatures>`, string 
objects are implicitly cast to non-Unicode for HDF5 storage. To make you
aware of this, a warning is raised when this happens.

This is certainly no true Unicode compatibility, but mainly for convenience 
with the pure-Unicode Python 3 string type. Any string that is not castable 
as ascii upon creation of your data set, will hence still raise an error. 
For true Unicode support, look into the ``VLUnicodeAtom`` class.

Major API Changes
=================

The PyTables developers, `by popular demand`_, have taken this opportunity
that a major version number upgrade affords to implement significant API
changes.  We have tried to do this in such a way that will not immediately
break most existing code, though in some breakages may still occur.

PEP 8 Compliance
****************
The PyTables 3.x series now follows `PEP 8`_ coding standard.  This makes
using PyTables more idiomatic with surrounding Python code that also adheres
to this standard.  The primary way that the 2.x series was *not* PEP 8
compliant was with respect to variable naming conventions.  Approximately
:ref:`450 API variables <api-name-changes>` were identified and updated for
PyTables 3.x.

To ease migration, PyTables ships with a new ``pt2to3`` command line tool.
This tool will run over a file and replace any instances of the old variable
names with the 3.x version of the name.  This tool covers the overwhelming
majority of cases was used to transition the PyTables code base itself!  However,
it may also accidentally also pick up variable names in 3rd party codes that
have *exactly* the same name as a PyTables' variable.  This is because ``pt2to3``
was implemented using regular expressions rather than a fancier AST-based
method. By using regexes, ``pt2to3`` works on Python and Cython code.


``pt2to3`` **help:**

.. code-block:: bash

    usage: pt2to3 [-h] [-r] [-p] [-o OUTPUT] [-i] filename

    PyTables 2.x -> 3.x API transition tool This tool displays to standard out, so
    it is common to pipe this to another file: $ pt2to3 oldfile.py > newfile.py

    positional arguments:
      filename              path to input file.

    optional arguments:
      -h, --help            show this help message and exit
      -r, --reverse         reverts changes, going from 3.x -> 2.x.
      -p, --no-ignore-previous
                            ignores previous_api() calls.
      -o OUTPUT             output file to write to.
      -i, --inplace         overwrites the file in-place.

Note that ``pt2to3`` only works on a single file, not a directory.  However,
a simple BASH script may be written to run ``pt2to3`` over an entire directory
and all sub-directories:

.. code-block:: bash

    #!/bin/bash
    for f in $(find .)
    do
        echo $f
        pt2to3 $f > temp.txt
        mv temp.txt $f
    done

.. note::

    :program:`pt2to3` uses the :mod:`argparse` module that is part of the
    Python standard library since Python 2.7.
    Users of Python 2.6 should install :mod:`argparse` separately
    (e.g. via :program:`pip`).

The old APIs and variable names will continue to be supported for the short term,
where possible.  (The major backwards incompatible changes come from the renaming
of some function and method arguments and keyword arguments.)  Using the 2.x APIs
in the 3.x series, however, will issue warnings.  The following is the release
plan for the warning types:

* 3.0 - PendingDeprecationWarning
* 3.1 - DeprecationWarning
* >=3.2 - Remove warnings, previous_api(), and _past.py; keep pt2to3,

The current plan is to maintain the old APIs for at least 2 years, though this
is subject to change.

.. _create-signatures:

Consistent ``create_xxx()`` Signatures
***************************************

Also by popular demand, it is now possible to create all data sets (``Array``,
``CArray``, ``EArray``, ``VLArray``, and ``Table``) from existing Python objects.
Constructors for these classes now accept either of the following keyword arguments:

* an ``obj`` to initialize with data
* or both ``atom`` and ``shape`` to initialize an empty structure, if possible.

These keyword arguments are also now part of the function signature for the
corresponding ``create_xxx()`` methods on the ``File`` class.  These would be called
as follows::

    # All create methods will support the following
    create_xxx(where, name, obj=obj)

    # All non-variable length arrays support the following:
    create_xxx(where, name, atom=atom, shape=shape)

Using ``obj`` or ``atom`` and ``shape`` are mutually exclusive. Previously only
``Array`` could be created with an existing Python object using the ``object``
keyword argument.


.. _api-name-changes:

API Name Changes
****************

The following tables shows the old 2.x names that have been update to their
new values in the new 3.x series.  Please use the ``pt2to3`` tool to convert
between these.

================================ ================================
**2.x Name**                     **3.x Name**
================================ ================================
AtomFromHDF5Type                 atom_from_hdf5_type
AtomToHDF5Type                   atom_to_hdf5_type
BoolTypeNextAfter                bool_type_next_after
HDF5ClassToString                hdf5_class_to_string
HDF5ToNPExtType                  hdf5_to_np_ext_type
HDF5ToNPNestedType               hdf5_to_np_nested_type
IObuf                            iobuf
IObufcpy                         iobufcpy
IntTypeNextAfter                 int_type_next_after
NPExtPrefixesToPTKinds           npext_prefixes_to_ptkinds
PTSpecialKinds                   pt_special_kinds
PTTypeToHDF5                     pttype_to_hdf5
StringNextAfter                  string_next_after
__allowedInitKwArgs              __allowed_init_kwargs
__getRootGroup                   __get_root_group
__next__inKernel                 __next__inkernel
_actionLogName                   _action_log_name
_actionLogParent                 _action_log_parent
_actionLogPath                   _action_log_path
_addRowsToIndex                  _add_rows_to_index
_appendZeros                     _append_zeros
_autoIndex                       _autoindex
_byteShape                       _byte_shape
_c_classId                       _c_classid
_c_shadowNameRE                  _c_shadow_name_re
_cacheDescriptionData            _cache_description_data
_checkAndSetPair                 _check_and_set_pair
_checkAttributes                 _check_attributes
_checkBase                       _checkbase
_checkColumn                     _check_column
_checkGroup                      _check_group
_checkNotClosed                  _check_not_closed
_checkOpen                       _check_open
_checkShape                      _check_shape
_checkShapeAppend                _check_shape_append
_checkUndoEnabled                _check_undo_enabled
_checkWritable                   _check_writable
_check_sortby_CSI                _check_sortby_csi
_closeFile                       _close_file
_codeToOp                        _code_to_op
_column__createIndex             _column__create_index
_compileCondition                _compile_condition
_conditionCache                  _condition_cache
_convertTime64                   _convert_time64
_convertTime64_                  _convert_time64_
_convertTypes                    _convert_types
_createArray                     _create_array
_createCArray                    _create_carray
_createMark                      _create_mark
_createPath                      _create_path
_createTable                     _create_table
_createTransaction               _create_transaction
_createTransactionGroup          _create_transaction_group
_disableIndexingInQueries        _disable_indexing_in_queries
_doReIndex                       _do_reindex
_emptyArrayCache                 _empty_array_cache
_enableIndexingInQueries         _enable_indexing_in_queries
_enabledIndexingInQueries        _enabled_indexing_in_queries
_exprvarsCache                   _exprvars_cache
_f_copyChildren                  _f_copy_children
_f_delAttr                       _f_delattr
_f_getAttr                       _f_getattr
_f_getChild                      _f_get_child
_f_isVisible                     _f_isvisible
_f_iterNodes                     _f_iter_nodes
_f_listNodes                     _f_list_nodes
_f_setAttr                       _f_setattr
_f_walkGroups                    _f_walk_groups
_f_walkNodes                     _f_walknodes
_fancySelection                  _fancy_selection
_fillCol                         _fill_col
_flushBufferedRows               _flush_buffered_rows
_flushFile                       _flush_file
_flushModRows                    _flush_mod_rows
_g_addChildrenNames              _g_add_children_names
_g_checkGroup                    _g_check_group
_g_checkHasChild                 _g_check_has_child
_g_checkName                     _g_check_name
_g_checkNotContains              _g_check_not_contains
_g_checkOpen                     _g_check_open
_g_closeDescendents              _g_close_descendents
_g_closeGroup                    _g_close_group
_g_copyAsChild                   _g_copy_as_child
_g_copyChildren                  _g_copy_children
_g_copyRows                      _g_copy_rows
_g_copyRows_optim                _g_copy_rows_optim
_g_copyWithStats                 _g_copy_with_stats
_g_createHardLink                _g_create_hard_link
_g_delAndLog                     _g_del_and_log
_g_delLocation                   _g_del_location
_g_flushGroup                    _g_flush_group
_g_getAttr                       _g_getattr
_g_getChildGroupClass            _g_get_child_group_class
_g_getChildLeafClass             _g_get_child_leaf_class
_g_getGChildAttr                 _g_get_gchild_attr
_g_getLChildAttr                 _g_get_lchild_attr
_g_getLinkClass                  _g_get_link_class
_g_listAttr                      _g_list_attr
_g_listGroup                     _g_list_group
_g_loadChild                     _g_load_child
_g_logAdd                        _g_log_add
_g_logCreate                     _g_log_create
_g_logMove                       _g_log_move
_g_maybeRemove                   _g_maybe_remove
_g_moveNode                      _g_move_node
_g_postInitHook                  _g_post_init_hook
_g_postReviveHook                _g_post_revive_hook
_g_preKillHook                   _g_pre_kill_hook
_g_propIndexes                   _g_prop_indexes
_g_readCoords                    _g_read_coords
_g_readSelection                 _g_read_selection
_g_readSlice                     _g_read_slice
_g_readSortedSlice               _g_read_sorted_slice
_g_refNode                       _g_refnode
_g_removeAndLog                  _g_remove_and_log
_g_setAttr                       _g_setattr
_g_setLocation                   _g_set_location
_g_setNestedNamesDescr           _g_set_nested_names_descr
_g_setPathNames                  _g_set_path_names
_g_unrefNode                     _g_unrefnode
_g_updateDependent               _g_update_dependent
_g_updateLocation                _g_update_location
_g_updateNodeLocation            _g_update_node_location
_g_updateTableLocation           _g_update_table_location
_g_widthWarning                  _g_width_warning
_g_writeCoords                   _g_write_coords
_g_writeSelection                _g_write_selection
_g_writeSlice                    _g_write_slice
_getColumnInstance               _get_column_instance
_getConditionKey                 _get_condition_key
_getContainer                    _get_container
_getEnumMap                      _get_enum_map
_getFileId                       _get_file_id
_getFinalAction                  _get_final_action
_getInfo                         _get_info
_getLinkClass                    _get_link_class
_getMarkID                       _get_mark_id
_getNode                         _get_node
_getOrCreatePath                 _get_or_create_path
_getTypeColNames                 _get_type_col_names
_getUnsavedNrows                 _get_unsaved_nrows
_getValueFromContainer           _get_value_from_container
_hiddenNameRE                    _hidden_name_re
_hiddenPathRE                    _hidden_path_re
_indexNameOf                     _index_name_of
_indexNameOf_                    _index_name_of_
_indexPathnameOf                 _index_pathname_of
_indexPathnameOfColumn           _index_pathname_of_column
_indexPathnameOfColumn_          _index_pathname_of_column_
_indexPathnameOf_                _index_pathname_of_
_initLoop                        _init_loop
_initSortedSlice                 _init_sorted_slice
_isWritable                      _iswritable
_is_CSI                          _is_csi
_killNode                        _killnode
_lineChunkSize                   _line_chunksize
_lineSeparator                   _line_separator
_markColumnsAsDirty              _mark_columns_as_dirty
_newBuffer                       _new_buffer
_notReadableError                _not_readable_error
_npSizeType                      _npsizetype
_nxTypeFromNPType                _nxtype_from_nptype
_opToCode                        _op_to_code
_openArray                       _open_array
_openUnImplemented               _open_unimplemented
_pointSelection                  _point_selection
_processRange                    _process_range
_processRangeRead                _process_range_read
_pythonIdRE                      _python_id_re
_reIndex                         _reindex
_readArray                       _read_array
_readCoordinates                 _read_coordinates
_readCoords                      _read_coords
_readIndexSlice                  _read_index_slice
_readSelection                   _read_selection
_readSlice                       _read_slice
_readSortedSlice                 _read_sorted_slice
_refNode                         _refnode
_requiredExprVars                _required_expr_vars
_reservedIdRE                    _reserved_id_re
_reviveNode                      _revivenode
_saveBufferedRows                _save_buffered_rows
_searchBin                       _search_bin
_searchBinNA_b                   _search_bin_na_b
_searchBinNA_d                   _search_bin_na_d
_searchBinNA_e                   _search_bin_na_e
_searchBinNA_f                   _search_bin_na_f
_searchBinNA_g                   _search_bin_na_g
_searchBinNA_i                   _search_bin_na_i
_searchBinNA_ll                  _search_bin_na_ll
_searchBinNA_s                   _search_bin_na_s
_searchBinNA_ub                  _search_bin_na_ub
_searchBinNA_ui                  _search_bin_na_ui
_searchBinNA_ull                 _search_bin_na_ull
_searchBinNA_us                  _search_bin_na_us
_setAttributes                   _set_attributes
_setColumnIndexing               _set_column_indexing
_shadowName                      _shadow_name
_shadowParent                    _shadow_parent
_shadowPath                      _shadow_path
_sizeToShape                     _size_to_shape
_tableColumnPathnameOfIndex      _table_column_pathname_of_index
_tableFile                       _table_file
_tablePath                       _table_path
_table__autoIndex                _table__autoindex
_table__getautoIndex             _table__getautoindex
_table__setautoIndex             _table__setautoindex
_table__whereIndexed             _table__where_indexed
_transGroupName                  _trans_group_name
_transGroupParent                _trans_group_parent
_transGroupPath                  _trans_group_path
_transName                       _trans_name
_transParent                     _trans_parent
_transPath                       _trans_path
_transVersion                    _trans_version
_unrefNode                       _unrefnode
_updateNodeLocations             _update_node_locations
_useIndex                        _use_index
_vShape                          _vshape
_vType                           _vtype
_v__nodeFile                     _v__nodefile
_v__nodePath                     _v__nodepath
_v_colObjects                    _v_colobjects
_v_maxGroupWidth                 _v_max_group_width
_v_maxTreeDepth                  _v_maxtreedepth
_v_nestedDescr                   _v_nested_descr
_v_nestedFormats                 _v_nested_formats
_v_nestedNames                   _v_nested_names
_v_objectID                      _v_objectid
_whereCondition                  _where_condition
_writeCoords                     _write_coords
_writeSelection                  _write_selection
_writeSlice                      _write_slice
appendLastRow                    append_last_row
attrFromShadow                   attr_from_shadow
attrToShadow                     attr_to_shadow
autoIndex                        autoindex
bufcoordsData                    bufcoords_data
calcChunksize                    calc_chunksize
checkFileAccess                  check_file_access
checkNameValidity                check_name_validity
childName                        childname
chunkmapData                     chunkmap_data
classIdDict                      class_id_dict
className                        classname
classNameDict                    class_name_dict
containerRef                     containerref
convertToNPAtom                  convert_to_np_atom
convertToNPAtom2                 convert_to_np_atom2
copyChildren                     copy_children
copyClass                        copyclass
copyFile                         copy_file
copyLeaf                         copy_leaf
copyNode                         copy_node
copyNodeAttrs                    copy_node_attrs
countLoggedInstances             count_logged_instances
createArray                      create_array
createCArray                     create_carray
createCSIndex                    create_csindex
createEArray                     create_earray
createExternalLink               create_external_link
createGroup                      create_group
createHardLink                   create_hard_link
createIndex                      create_index
createIndexesDescr               create_indexes_descr
createIndexesTable               create_indexes_table
createNestedType                 create_nested_type
createSoftLink                   create_soft_link
createTable                      create_table
createVLArray                    create_vlarray
defaultAutoIndex                 default_auto_index
defaultIndexFilters              default_index_filters
delAttr                          del_attr
delAttrs                         _del_attrs
delNodeAttr                      del_node_attr
detectNumberOfCores              detect_number_of_cores
disableUndo                      disable_undo
dumpGroup                        dump_group
dumpLeaf                         dump_leaf
dumpLoggedInstances              dump_logged_instances
enableUndo                       enable_undo
enumFromHDF5                     enum_from_hdf5
enumToHDF5                       enum_to_hdf5
fetchLoggedInstances             fetch_logged_instances
flushRowsToIndex                 flush_rows_to_index
getAttr                          get_attr
getAttrs                         _get_attrs
getClassByName                   get_class_by_name
getColsInOrder                   get_cols_in_order
getCurrentMark                   get_current_mark
getEnum                          get_enum
getFilters                       get_filters
getHDF5Version                   get_hdf5_version
getIndices                       get_indices
getLRUbounds                     get_lru_bounds
getLRUsorted                     get_lru_sorted
getLookupRange                   get_lookup_range
getNestedField                   get_nested_field
getNestedFieldCache              get_nested_field_cache
getNestedType                    get_nested_type
getNode                          get_node
getNodeAttr                      get_node_attr
getPyTablesVersion               get_pytables_version
getTypeEnum                      get_type_enum
getWhereList                     get_where_list
hdf5Extension                    hdf5extension
hdf5Version                      hdf5_version
indexChunk                       indexchunk
indexValid                       indexvalid
indexValidData                   index_valid_data
indexValues                      indexvalues
indexValuesData                  index_values_data
indexesExtension                 indexesextension
infType                          inftype
infinityF                        infinityf
infinityMap                      infinitymap
initRead                         initread
isHDF5File                       is_hdf5_file
isPyTablesFile                   is_pytables_file
isUndoEnabled                    is_undo_enabled
isVisible                        isvisible
isVisibleName                    isvisiblename
isVisibleNode                    is_visible_node
isVisiblePath                    isvisiblepath
is_CSI                           is_csi
iterNodes                        iter_nodes
iterseqMaxElements               iterseq_max_elements
joinPath                         join_path
joinPaths                        join_paths
linkExtension                    linkextension
listLoggedInstances              list_logged_instances
listNodes                        list_nodes
loadEnum                         load_enum
logInstanceCreation              log_instance_creation
lrucacheExtension                lrucacheextension
metaIsDescription                MetaIsDescription
modifyColumn                     modify_column
modifyColumns                    modify_columns
modifyCoordinates                modify_coordinates
modifyRows                       modify_rows
moveFromShadow                   move_from_shadow
moveNode                         move_node
moveToShadow                     move_to_shadow
newNode                          new_node
newSet                           newset
newdstGroup                      newdst_group
objectID                         object_id
oldPathname                      oldpathname
openFile                         open_file
openNode                         open_node
parentNode                       parentnode
parentPath                       parentpath
reIndex                          reindex
reIndexDirty                     reindex_dirty
readCoordinates                  read_coordinates
readIndices                      read_indices
readSlice                        read_slice
readSorted                       read_sorted
readWhere                        read_where
read_sliceLR                     read_slice_lr
recreateIndexes                  recreate_indexes
redoAddAttr                      redo_add_attr
redoCreate                       redo_create
redoDelAttr                      redo_del_attr
redoMove                         redo_move
redoRemove                       redo_remove
removeIndex                      remove_index
removeNode                       remove_node
removeRows                       remove_rows
renameNode                       rename_node
rootUEP                          root_uep
searchLastRow                    search_last_row
setAttr                          set_attr
setAttrs                         _set_attrs
setBloscMaxThreads               set_blosc_max_threads
setInputsRange                   set_inputs_range
setNodeAttr                      set_node_attr
setOutput                        set_output
setOutputRange                   set_output_range
silenceHDF5Messages              silence_hdf5_messages
splitPath                        split_path
tableExtension                   tableextension
undoAddAttr                      undo_add_attr
undoCreate                       undo_create
undoDelAttr                      undo_del_attr
undoMove                         undo_move
undoRemove                       undo_remove
utilsExtension                   utilsextension
walkGroups                       walk_groups
walkNodes                        walk_nodes
whereAppend                      append_where
whereCond                        wherecond
whichClass                       which_class
whichLibVersion                  which_lib_version
willQueryUseIndexing             will_query_use_indexing
================================ ================================

----

  **Enjoy data!**

  -- The PyTables Developers


.. Local Variables:
.. mode: rst
.. coding: utf-8
.. fill-column: 78
.. End:


.. _by popular demand: http://sourceforge.net/mailarchive/message.php?msg_id=29584752

.. _PEP 8: http://www.python.org/dev/peps/pep-0008/

.. _netcdf4-python: http://code.google.com/p/netcdf4-python/
================
Project pointers
================

* `Project Home Page <http://www.pytables.org>`_
* `GitHub Project Page <https://github.com/PyTables>`_
* `Online HTML Documentation <http://www.pytables.org>`_
* `Download area <http://sourceforge.net/projects/pytables/files/pytables>`_
* `Git Repository browser <https://github.com/PyTables/PyTables>`_
* `Users Mailing List <https://groups.google.com/group/pytables-users>`_
* `Announce Mailing List <https://lists.sourceforge.net/lists/listinfo/pytables-announce>`_
* `Developers Mailing List <https://groups.google.com/group/pytables-dev>`_
* Continuous Integration:

  - `GitHub Actions (GHA) <https://github.com/PyTables/PyTables/actions>`_

* `Project page on PyPi <https://pypi.python.org/pypi/tables>`_
* `Project Page on SourceForge.net <http://sourceforge.net/projects/pytables>`_
  (needs update)
* `Project page on Launchpad <https://launchpad.net/pytables>`_
  (going to be closed)
========================
PyTables Governance Team
========================

The PyTables team includes:

* Francesc Alted
* Ivan Vilata
* Scott Prater
* Vicent Mas
* Tom Hedley
* `Antonio Valentino`_
* Jeffrey Whitaker
* `Josh Moore`_
* `Anthony Scopatz`_
* `Andrea Bedini`_
* `Tom Kooij`_
* `Javier Sancho`_

.. _Anthony Scopatz: https://github.com/scopatz
.. _Antonio Valentino: https://github.com/avalentino
.. _Josh Moore: https://github.com/joshmoore
.. _Andrea Bedini: https://github.com/andreabedini
.. _Tom Kooij: https://github.com/tomkooij
.. _Javier Sancho: https://en.jsancho.org/
==================================
Migrating from PyTables 1.x to 2.x
==================================

:Author: Francesc Alted i Abad
:Contact: faltet@pytables.com
:Author: Ivan Vilata i Balaguer
:Contact: ivan@selidor.net


Next are described a series of issues that you must have in mind when
migrating from PyTables 1.x to PyTables 2.x series.


New type system
===============

In PyTables 2.x all the data types for leaves are described through a couple
of classes:

- ``Atom``: Describes homogeneous types of the atomic components in ``*Array``
   objects (``Array``, ``CArray``, ``EArray`` and ``VLArray``).

- ``Description``: Describes (possibly nested) heterogeneous types in
  ``Table`` objects.

So, in order to upgrade to the new type system, you must perform the next
replacements:

- ``*Array.stype`` --> ``*Array.atom.type`` (PyTables type)
- ``*Array.type`` --> ``*Array.atom.dtype`` (NumPy type)
- ``*Array.itemsize`` --> ``*Array.atom.itemsize`` (the size of the item)

Furthermore, the PyTables types (previously called "string types") have
changed to better adapt to NumPy conventions.  The next changes have been
applied:

- PyTables types are now written in lower case, so 'Type' becomes 'type'.  For
  example, 'Int64' becomes now 'int64'.

- 'CharType' --> 'string'

- 'Complex32', 'Complex64' --> 'complex64', 'complex128'.  Note that the
  numeric part of a 'complex' type refers now to the *size in bits* of the
  type and not to the precision, as before.

See Appendix I of the Users' Manual on supported data types for more
information on the new PyTables types.


Important changes in ``Atom`` specification
===========================================

- The ``dtype`` argument of ``EnumAtom`` and ``EnumCol`` constructors
  has been replaced by the ``base`` argument, which can take a
  full-blown atom, although it accepts bare PyTables types as well.
  This is a *mandatory* argument now.

- ``vlstring`` pseudo-atoms used in ``VLArray`` nodes do no longer imply UTF-8
  (nor any other) encoding, they only store and load *raw strings of bytes*.
  All encoding and decoding is left to the user.  Be warned that reading old
  files may yield raw UTF-8 encoded strings, which may be converted back to
  Unicode in this way::

      unistr = vlarray[index].decode('utf-8')

  If you need to work with variable-length Unicode strings, you may want to
  use the new ``vlunicode`` pseudo-atom, which fully supports Unicode strings
  with no encoding hassles.

- Finally, ``Atom`` and ``Col`` are now abstract classes, so you can't use
  them to create atoms or column definitions of an arbitrary type.  If you
  know the particular type you need, use the proper subclass; otherwise, use
  the ``Atom.from_*()`` or ``Col.from_*()`` factory methods.  See the section
  on declarative classes in the reference.

  You are also advised to avoid using the inheritance of atoms to check for
  their kind or type; for that purpose, use their ``kind`` and ``type``
  attributes.


New query system
================

- In-kernel conditions, since they are based now in Numexpr, must be written
  *as strings*.  For example, a condition that in 1.x was stated as::

      result = [row['col2'] for row in table.where(table.cols.col1 == 1)]

  now should read::

      result = [row['col2'] for row in table.where('col1 == 1')]

  That means that complex selections are possible now::

      result = [ row['col2'] for row in
                 table.where('(col1 == 1) & (col3**4 > 1)') ]

- For the same reason, conditions for indexed columns must be written as
  strings as well.


New indexing system
===================

The indexing system has been totally rewritten from scratch for PyTables 2.0
Pro Edition.  The new indexing system has been included into PyTables with
release 2.3.  Due to this, your existing indexes created with PyTables 1.x
will be useless, and although you will be able to continue using the actual
data in files, you won't be able to take advantage of any improvement in
speed.

You will be offered the possibility to automatically re-create the indexes
in PyTables 1.x format to the new 2.0 format by using the ``ptrepack``
utility.


New meanings for atom shape and ``*Array`` shape argument
=========================================================

With PyTables 1.x, the atom shape was used for different goals depending on
the context it was used.  For example, in ``createEArray()``, the shape of the
atom was used to specify the *dataset shape* of the object on disk, while in
``CArray`` the same atom shape was used to specify the *chunk shape* of the
dataset on disk.  Moreover, for ``VLArray`` objects, the very same atom shape
specified the *type shape* of the data type.  As you see, all of these was
quite a mess.

Starting with PyTables 2.x, an ``Atom`` only specifies properties of the data
type (à la ``VLArray`` in 1.x).  This lets the door open for specifying
multidimensional data types (that can be part of another layer of
multidimensional datasets) in a consistent way along all the ``*Array``
objects in PyTables.

As a consequence of this, ``File.createCArray()`` and ``File.createVLArray()``
methods have received new parameters in order to make possible to specify the
shapes of the datasets as well as chunk sizes (in fact, it is possible now to
specify the latter for all the chunked leaves, see below).  Please have this
in mind during the migration process.

Another consequence is that, now that the meaning of the atom shape is clearly
defined, it has been chosen as the main object to describe homogeneous data
types in PyTables.  See the Users' Manual for more info on this.


New argument ``chunkshape`` of chunked leaves
=============================================

It is possible now to specify the chunk shape for all the chunked leaves in
PyTables (all except ``Array``).  With PyTables 1.x this value was
automatically calculated so as to achieve decent results in most of the
situations.  However, the user may be interested in specifying its own chunk
shape based on her own needs (although this should be done only by advanced
users).

Of course, if this parameter is not specified, a sensible default is
calculated for the size of the leave (which is recommended).

A new attribute called ``chunkshape`` has been added to all leaves.  It is
read-only (you can't change the size of chunks once you have created a leaf),
but it can be useful for inspection by advanced users.


New flavor specification
========================

As of 2.x, flavors can *only* be set through the ``flavor`` attribute of
leaves, and they are *persistent*, so changing a flavor requires that the file
be writable.

Flavors can no longer be set through ``File.create*()`` methods, nor the
``flavor`` argument previously found in some ``Table`` methods, nor through
``Atom`` constructors or the ``_v_flavor`` attribute of descriptions.


System attributes can be deleted now
====================================

The protection against removing system attributes (like ``FILTERS``,
``FLAVOR`` or ``CLASS``, to name only a few) has been completely removed.  It
is now the responsibility of the user to make a proper use of this freedom.
With this, users can get rid of all proprietary PyTables attributes if they
want to (for example, for making a file to look more like an HDF5 native one).


Byteorder issues
================

Now, all the data coming from reads and internal buffers is always converted
on-the-fly, if needed, to the *native* byteorder.  This represents a big
advantage in terms of speed when operating with objects coming from files that
have been created in machines with a byte ordering different from native.

Besides, all leaf constructors have received a new ``byteorder`` parameter
that allows specifying the byteorder of data on disk.  In particular, a
``_v_byteorder`` entry in a Table description is no longer honored and you
should use the aforementioned ``byteorder`` parameter.


Tunable internal buffer sizes
=============================

You can change the size of the internal buffers for I/O purposes of PyTables
by changing the value of the new public attribute ``nrowsinbuf`` that is
present in all leaves.  By default, this contains a sensible value so as to
achieve a good balance between speed and memory consumption.  Be careful when
changing it, if you don't want to get unwanted results (very slow I/O, huge
memory consumption...).


Changes to module names
=======================

If your application is directly accessing modules under the ``tables``
package, you need to know that *the names of all modules are now all in
lowercase*.  This allows one to tell apart the ``tables.Array`` *class* from
the ``tables.array`` *module* (which was also called ``tables.Array`` before).
This includes subpackages like ``tables.nodes.FileNode``.

On top of that, more-or-less independent modules have also been renamed and
some of them grouped into subpackages.  The most important are:

- The ``tables.netcdf3`` subpackage replaces the old ``tables.NetCDF`` module.
- The ``tables.nra`` subpackage replaces the old ``nestedrecords.py`` with the
  implementation of the ``NestedRecArray`` class.

Also, the ``tables.misc`` package includes utility modules which do not depend
on PyTables.


Other changes
=============

- ``Filters.complib`` is ``None`` for filter properties created with
  ``complevel=0`` (i.e. disabled compression, which is the default).
- 'non-relevant' --> 'irrelevant' (applied to byteorders)
- ``Table.colstypes`` --> ``Table.coltypes``
- ``Table.coltypes`` --> ``Table.coldtypes``
- Added ``Table.coldescr``, dictionary of the ``Col`` descriptions.
- ``Table.colshapes`` has disappeared.  You can get it this way::

       colshapes = dict( (name, col.shape)
                         for (name, col) in table.coldescr.iteritems() )

- ``Table.colitemsizes`` has disappeared.  You can get it this way::

       colitemsizes = dict( (name, col.itemsize)
                            for (name, col) in table.coldescr.iteritems() )

- ``Description._v_totalsize`` --> ``Description._v_itemsize``
- ``Description._v_itemsizes`` and ``Description._v_totalsizes`` have
  disappeared.

- ``Leaf._v_chunksize`` --> ``Leaf.chunkshape``


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: rst
.. coding: utf-8
.. fill-column: 78
.. End:
===================================
Welcome to PyTables' documentation!
===================================

PyTables is a package for managing hierarchical datasets and designed
to efficiently and easily cope with extremely large amounts of data.
You can download PyTables and use it for free. You can access documentation,
some examples of use and presentations here.

PyTables is built on top of the HDF5 library, using the Python language
and the NumPy package. It features an object-oriented interface that,
combined with C extensions for the performance-critical parts of the
code (generated using Cython), makes it a fast, yet extremely easy to
use tool for interactively browsing, processing and searching very large
amounts of data. One important feature of PyTables is that it optimizes
memory and disk resources so that data takes much less space (specially 
if on-flight compression is used) than other solutions such as relational
object oriented databases.

You can also find more information by reading the PyTables :doc:`FAQ`.

PyTables development is a continuing effort and we are always looking for
more developers, testers, and users.  If you are interested in being
involved with this project, please contact us via `github`_ or the
`mailing list`_.

.. image:: images/NumFocusSponsoredStamp.png
   :alt: NumFocus Sponsored Stamp
   :align: center
   :width: 300
   :target: http://www.numfocus.org

Since August 2015, PyTables is a `NumFOCUS project`_, which means that
your donations are fiscally sponsored under the NumFOCUS umbrella.  Please
consider `donating to NumFOCUS`_.


--------
Contents
--------

.. toctree::
    :maxdepth: 1

    User’s Guide <usersguide/index>
    Cookbook <cookbook/index>
    FAQ
    other_material
    Migrating from 2.x to 3.x <MIGRATING_TO_3.x>
    downloads
    Release Notes <release_notes>
    project_pointers
    Development <development>
    Development Team <dev_team>


=============
Helpful Links
=============

* :ref:`genindex`
* :ref:`search`


.. _github: https://github.com/PyTables/PyTables
.. _`mailing list`: https://groups.google.com/group/pytables-users
.. _`NumFOCUS project`: http://www.numfocus.org/open-source-projects.html
.. _`donating to NumFOCUS`: https://numfocus.salsalabs.org/donate-to-pytables/index.html
=========
Downloads
=========

Stable Versions
---------------

The stable versions of PyTables can be downloaded from the file `download
area`_ on SourceForge.net.  The full distribution contains a copy of this
documentation in HTML.  The documentation in both HTML and PDF formats can
also be downloaded separately from the same URL.

A *pure source* version of the package (mainly intended for developers and
packagers) is available on the `tags page`_ on GitHub.  It contains all files
under SCM but not the (generated) files, HTML doc and *cythonized* C
extensions, so it is smaller that the standard package (about 3.5MB).

Windows binaries can be obtained from many different distributions, like
`Python(x,y)`_, ActiveState_, or Enthought_.
In addition, Christoph Gohlke normally does an excellent job by providing
binaries for many interesting software on his
`website <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_.

You may be interested to install the latest released stable version::

    $ pip install tables

Or, you may prefer to install the stable version in Git repository
using :program:`pip`. For example, for the stable 3.1 series, you can do::

    $ pip install --install-option='--prefix=<PREFIX>' \
    -e git+https://github.com/PyTables/PyTables.git@v.3.1#egg=tables

.. _`download area`: http://sourceforge.net/projects/pytables/files/pytables
.. _`tags page`: https://github.com/PyTables/PyTables/tags
.. _`Python(x,y)`: http://code.google.com/p/pythonxy
.. _ActiveState: http://www.activestate.com/activepython
.. _Enthought: https://www.enthought.com/products/epd


Bleeding Edge Versions
----------------------

The latest, coolest, and possibly buggiest ;-) sources can be obtained from
the new github repository:

https://github.com/PyTables/PyTables

A
`snapshot <https://github.com/PyTables/PyTables/archive/refs/heads/master.zip>`_
of the code in development is also available on the `GitHub project page`_.

.. _`GitHub project page`: https://github.com/PyTables/PyTables

.. _filenode_usersguide:

filenode - simulating a filesystem with PyTables
================================================

.. currentmodule:: tables.nodes.filenode

What is filenode?
-----------------
filenode is a module which enables you to create a PyTables database of nodes
which can be used like regular opened files in Python. In other words, you
can store a file in a PyTables database, and read and write it as you would
do with any other file in Python. Used in conjunction with PyTables
hierarchical database organization, you can have your database turned into an
open, extensible, efficient, high capacity, portable and metadata-rich
filesystem for data exchange with other systems (including backup purposes).

Between the main features of filenode, one can list:

- *Open:* Since it relies on PyTables, which in turn, sits over HDF5 (see
  :ref:`[HDGG1] <HDFG1>`), a standard hierarchical data format from NCSA.

- *Extensible:* You can define new types of nodes, and their instances will
  be safely preserved (as are normal groups, leafs and attributes) by
  PyTables applications having no knowledge of their types. Moreover, the set
  of possible attributes for a node is not fixed, so you can define your own
  node attributes.

- *Efficient:* Thanks to PyTables' proven extreme efficiency on handling huge
  amounts of data. filenode can make use of PyTables' on-the-fly compression
  and decompression of data.

- *High capacity:* Since PyTables and HDF5 are designed for massive data
  storage (they use 64-bit addressing even where the platform does not
  support it natively).

- *Portable:* Since the HDF5 format has an architecture-neutral design, and
  the HDF5 libraries and PyTables are known to run under a variety of
  platforms. Besides that, a PyTables database fits into a single file, which
  poses no trouble for transportation.

- *Metadata-rich:* Since PyTables can store arbitrary key-value pairs (even
  Python objects!) for every database node. Metadata may include authorship,
  keywords, MIME types and encodings, ownership information, access control
  lists (ACL), decoding functions and anything you can imagine!


Finding a filenode node
-----------------------
filenode nodes can be recognized because they have a NODE_TYPE system
attribute with a 'file' value. It is recommended that you use the
:meth:`File.get_node_attr` method of tables.File class to get the NODE_TYPE
attribute independently of the nature (group or leaf) of the node, so you do
not need to care about.


filenode - simulating files inside PyTables
-------------------------------------------
The filenode module is part of the nodes sub-package of PyTables. The
recommended way to import the module is::

    >>> from tables.nodes import filenode

However, filenode exports very few symbols, so you can import * for
interactive usage. In fact, you will most probably only use the NodeType
constant and the new_node() and open_node() calls.

The NodeType constant contains the value that the NODE_TYPE system attribute
of a node file is expected to contain ('file', as we have seen).
Although this is not expected to change, you should use filenode.NodeType
instead of the literal 'file' when possible.

new_node() and open_node() are the equivalent to the Python file() call (alias
open()) for ordinary files. Their arguments differ from that of file(), but
this is the only point where you will note the difference between working
with a node file and working with an ordinary file.

For this little tutorial, we will assume that we have a PyTables database
opened for writing. Also, if you are somewhat lazy at typing sentences, the
code that we are going to explain is included in the examples/filenodes1.py
file.

You can create a brand new file with these sentences::

    >>> import tables
    >>> h5file = tables.open_file('fnode.h5', 'w')


Creating a new file node
~~~~~~~~~~~~~~~~~~~~~~~~
Creation of a new file node is achieved with the new_node() call. You must
tell it in which PyTables file you want to create it, where in the PyTables
hierarchy you want to create the node and which will be its name. The
PyTables file is the first argument to new_node(); it will be also called the
'host PyTables file'. The other two arguments must be given as keyword
arguments where and name, respectively.
As a result of the call, a brand new appendable and readable file node object
is returned.

So let us create a new node file in the previously opened h5file PyTables
file, named 'fnode_test' and placed right under the root of the database
hierarchy. This is that command::

    >>> fnode = filenode.new_node(h5file, where='/', name='fnode_test')

That is basically all you need to create a file node. Simple, isn't it? From
that point on, you can use fnode as any opened Python file (i.e. you can
write data, read data, lines of text and so on).

new_node() accepts some more keyword arguments. You can give a title to your
file with the title argument. You can use PyTables' compression features with
the filters argument. If you know beforehand the size that your file will
have, you can give its final file size in bytes to the expectedsize argument
so that the PyTables library would be able to optimize the data access.

new_node() creates a PyTables node where it is told to. To prove it, we will
try to get the NODE_TYPE attribute from the newly created node::

    >>> print(h5file.get_node_attr('/fnode_test', 'NODE_TYPE'))
    file


Using a file node
~~~~~~~~~~~~~~~~~
As stated above, you can use the new node file as any other opened file. Let
us try to write some text in and read it::

    >>> print("This is a test text line.", file=fnode)
    >>> print("And this is another one.", file=fnode)
    >>> print(file=fnode)
    >>> fnode.write("Of course, file methods can also be used.")
    >>>
    >>> fnode.seek(0)  # Go back to the beginning of file.
    >>>
    >>> for line in fnode:
    ...     print(repr(line))
    'This is a test text line.\\n'
    'And this is another one.\\n'
    '\\n'
    'Of course, file methods can also be used.'

This was run on a Unix system, so newlines are expressed as '\n'. In fact,
you can override the line separator for a file by setting its line_separator
property to any string you want.

While using a file node, you should take care of closing it *before* you
close the PyTables host file.
Because of the way PyTables works, your data it will not be at a risk, but
every operation you execute after closing the host file will fail with a
ValueError. To close a file node, simply delete it or call its close()
method::

    >>> fnode.close()
    >>> print(fnode.closed)
    True


Opening an existing file node
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you have a file node that you created using new_node(), you can open it
later by calling open_node(). Its arguments are similar to that of file() or
open(): the first argument is the PyTables node that you want to open (i.e. a
node with a NODE_TYPE attribute having a 'file' value), and the second
argument is a mode string indicating how to open the file. Contrary to
file(), open_node() can not be used to create a new file node.

File nodes can be opened in read-only mode ('r') or in read-and-append mode
('a+'). Reading from a file node is allowed in both modes, but appending is
only allowed in the second one. Just like Python files do, writing data to an
appendable file places it after the file pointer if it is on or beyond the
end of the file, or otherwise after the existing data. Let us see an
example::

    >>> node = h5file.root.fnode_test
    >>> fnode = filenode.open_node(node, 'a+')
    >>> print(repr(fnode.readline()))
    'This is a test text line.\\n'
    >>> print(fnode.tell())
    26
    >>> print("This is a new line.", file=fnode)
    >>> print(repr(fnode.readline()))
    ''

Of course, the data append process places the pointer at the end of the file,
so the last readline() call hit EOF. Let us seek to the beginning of the file
to see the whole contents of our file::

    >>> fnode.seek(0)
    >>> for line in fnode:
    ...   print(repr(line))
    'This is a test text line.\\n'
    'And this is another one.\\n'
    '\\n'
    'Of course, file methods can also be used.This is a new line.\\n'

As you can check, the last string we wrote was correctly appended at the end
of the file, instead of overwriting the second line, where the file pointer
was positioned by the time of the appending.


Adding metadata to a file node
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You can associate arbitrary metadata to any open node file, regardless of its
mode, as long as the host PyTables file is writable. Of course, you could use
the set_node_attr() method of tables.File to do it directly on the proper node,
but filenode offers a much more comfortable way to do it. filenode objects
have an attrs property which gives you direct access to their corresponding
AttributeSet object.

For instance, let us see how to associate MIME type metadata to our file
node::

    >>> fnode.attrs.content_type = 'text/plain; charset=us-ascii'

As simple as A-B-C. You can put nearly anything in an attribute, which opens
the way to authorship, keywords, permissions and more. Moreover, there is not
a fixed list of attributes.
However, you should avoid names in all caps or starting with '_', since
PyTables and filenode may use them internally. Some valid examples::

    >>> fnode.attrs.author = "Ivan Vilata i Balaguer"
    >>> fnode.attrs.creation_date = '2004-10-20T13:25:25+0200'
    >>> fnode.attrs.keywords_en = ["FileNode", "test", "metadata"]
    >>> fnode.attrs.keywords_ca = ["FileNode", "prova", "metadades"]
    >>> fnode.attrs.owner = 'ivan'
    >>> fnode.attrs.acl = {'ivan': 'rw', '@users': 'r'}

You can check that these attributes get stored by running the ptdump command
on the host PyTables file.

.. code-block:: bash

    $ ptdump -a fnode.h5:/fnode_test
    /fnode_test (EArray(113,)) ''
    /fnode_test.attrs (AttributeSet), 14 attributes:
    [CLASS := 'EARRAY',
    EXTDIM := 0,
    FLAVOR := 'numpy',
    NODE_TYPE := 'file',
    NODE_TYPE_VERSION := 2,
    TITLE := '',
    VERSION := '1.2',
    acl := {'ivan': 'rw', '@users': 'r'},
    author := 'Ivan Vilata i Balaguer',
    content_type := 'text/plain; charset=us-ascii',
    creation_date := '2004-10-20T13:25:25+0200',
    keywords_ca := ['FileNode', 'prova', 'metadades'],
    keywords_en := ['FileNode', 'test', 'metadata'],
    owner := 'ivan']

Note that filenode makes no assumptions about the meaning of your metadata,
so its handling is entirely left to your needs and imagination.


Complementary notes
-------------------
You can use file nodes and PyTables groups to mimic a filesystem with files
and directories. Since you can store nearly anything you want as file
metadata, this enables you to use a PyTables file as a portable compressed
backup, even between radically different platforms. Take this with a grain of
salt, since node files are restricted in their naming (only valid Python
identifiers are valid); however, remember that you can use node titles and
metadata to overcome this limitation. Also, you may need to devise some
strategy to represent special files such as devices, sockets and such (not
necessarily using filenode).

We are eager to hear your opinion about filenode and its potential uses.
Suggestions to improve filenode and create other node types are also welcome.
Do not hesitate to contact us!


Current limitations
-------------------
filenode is still a young piece of software, so it lacks some functionality.
This is a list of known current limitations:

#. Node files can only be opened for read-only or read and append mode. This
   should be enhanced in the future.
#. Near future?
#. Only binary I/O is supported currently (read/write strings of bytes)
#. There is no universal newline support yet. The only new-line character
   used at the moment is ``\n``. This is likely to be improved in a near
   future.
#. Sparse files (files with lots of zeros) are not treated specially; if you
   want them to take less space, you should be better off using compression.

These limitations still make filenode entirely adequate to work with most
binary and text files. Of course, suggestions and patches are welcome.

See :ref:`filenode_classes` for detailed documentation on the filenode
interface.

Bibliography
============

.. _HDFG1:

:ref:`[HDFG1] <HDFG1>`
    The HDF Group. What is HDF5?. Concise description about HDF5 capabilities
    and its differences from earlier versions (HDF4).
    `<https://www.hdfgroup.org/solutions/hdf5>`_.

.. _HDFG2:

:ref:`[HDFG2] <HDFG2>`
    The HDF Group. Introduction to HDF5. Introduction to the HDF5 data model
    and programming model. `<http://www.hdfgroup.org/HDF5/doc/H5.intro.html>`_.

.. _HDFG3:

:ref:`[HDFG3] <HDFG3>`
    The HDF Group. The HDF5 table programming model. Examples on using HDF5
    tables with the C API. `<http://www.hdfgroup.org/HDF5/Tutor/h5table.html>`_.

.. _MERTZ:

:ref:`[MERTZ] <MERTZ>`
    David Mertz. Objectify. On the 'Pythonic' treatment of XML documents as
    objects(II). Article describing XML Objectify, a Python module that
    allows working with XML documents as Python objects.
    Some of the ideas presented here are used in PyTables.
    `<http://gnosis.cx/publish/programming/xml_matters_2.html>`_.

.. _CYTHON:

:ref:`[CYTHON] <CYTHON>`
    Stefan Behnel, Robert Bradshaw, Dag Sverre Seljebotn, and Greg Ewing.
    Cython. A language that makes writing C extensions for the Python
    language as easy as Python itself. `<http://www.cython.org>`_.

.. _NUMPY:

:ref:`[NUMPY] <NUMPY>`
    Travis Oliphant and et al. NumPy. Scientific Computing with Numerical
    Python. The latest and most powerful re-implementation of Numeric to
    date.
    It implements all the features that can be found in Numeric and numarray,
    plus a bunch of new others. In general, it is more efficient as well.
    `<http://www.numpy.org>`_.

.. _NUMEXPR:

:ref:`[NUMEXPR] <NUMEXPR>`
    David Cooke, Francesc Alted, and et al. Numexpr. Fast evaluation of array
    expressions by using a vector-based virtual machine.
    It is an enhaced computing kernel that is generally faster (between 1x
    and 10x, depending on the kind of operations) than NumPy at evaluating
    complex array expressions. `<http://code.google.com/p/numexpr>`_.

.. _ZLIB:

:ref:`[ZLIB] <ZLIB>`
    JeanLoup Gailly and Mark Adler. zlib. A Massively Spiffy Yet Delicately
    Unobtrusive Compression Library. A standard library for compression
    purposes. `<https://zlib.net>`_.

.. _LZO:

:ref:`[LZO] <LZO>`
    Markus F Oberhumer. LZO. A data compression library which is suitable for
    data de-/compression in real-time. It offers pretty fast compression and
    decompression with reasonable compression ratio.
    `<http://www.oberhumer.com/opensource/>`_.

.. _BZIP2:

:ref:`[BZIP2] <BZIP2>`
    Julian Seward. bzip2. A high performance lossless compressor.
    It offers very high compression ratios within reasonable times.
    `<http://www.bzip.org/>`_.

.. _BLOSC:

:ref:`[BLOSC] <BLOSC>`
    Francesc Alted. Blosc. A blocking, shuffling and loss-less compression
    library.  A compressor designed to transmit data from memory to CPU
    (and back) faster than a plain memcpy().
    `<http://www.blosc.org/>`_.

.. _GNUWIN32:

:ref:`[GNUWIN32] <GNUWIN32>`
    Alexis Wilke, Jerry S., Kees Zeelenberg, and Mathias Michaelis.
    GnuWin32. GNU (and other) tools ported to Win32.
    GnuWin32 provides native Win32-versions of GNU tools, or tools with a
    similar open source licence.
    `<http://gnuwin32.sourceforge.net/>`_.

.. _PSYCO:

:ref:`[PSYCO] <PSYCO>`
    Armin Rigo. Psyco. A Python specializing compiler.
    Run existing Python software faster, with no change in your source.
    `<http://psyco.sourceforge.net>`_.

.. _SCIPY1:

:ref:`[SCIPY1] <SCIPY1>`
    Konrad Hinsen. Scientific Python. Collection of Python modules useful for
    scientific computing.
    `<http://dirac.cnrs-orleans.fr/ScientificPython>`_.

.. _SCIPY2:

:ref:`[SCIPY2] <SCIPY2>`
    Eric Jones, Travis Oliphant, Pearu Peterson, and et al. SciPy.
    Scientific tools for Python. SciPy supplements the popular Numeric module,
    gathering a variety of high level science and engineering modules
    together as a single package.
    `<http://www.scipy.org>`_.

.. _OPTIM:

:ref:`[OPTIM] <OPTIM>`
    Francesc Alted and Ivan Vilata. Optimization of file openings in PyTables.
    This document explores the savings of the opening process in terms of
    both CPU time and memory, due to the adoption of a LRU cache for the
    nodes in the object tree.
    `<http://www.pytables.org/docs/NewObjectTreeCache.pdf>`_.

.. _OPSI:

:ref:`[OPSI] <OPSI>`
    Francesc Alted and Ivan Vilata. OPSI: The indexing system of PyTables 2
    Professional Edition. Exhaustive description and benchmarks about the
    indexing engine that comes with PyTables Pro.
    `<http://www.pytables.org/docs/OPSI-indexes.pdf>`_.

.. _VITABLES:

:ref:`[VITABLES] <VITABLES>`
    Vicent Mas. ViTables. A GUI for PyTables/HDF5 files.
    It is a graphical tool for browsing and editing files in both PyTables
    and HDF5 formats.
    `<http://vitables.org>`_.

.. _GIT:

:ref:`[GIT] <GIT>`
    Git is a free and open source, distributed version control system designed
    to handle everything from small to very large projects with speed and
    efficiency `<http://git-scm.com>`_.

.. _SPHINX:

:ref:`[SPHINX] <SPHINX>`
    Sphinx is a tool that makes it easy to create intelligent and beautiful
    documentation, written by Georg Brandl and licensed under the BSD license
    `<http://sphinx-doc.org>`_.

.. |Kuepper| unicode:: K U+00FC pper .. Kuepper

.. todo:: remove the above substitution. It is no more needed with sphinx
          1.0.8
.. _library_reference:

Library Reference
=================
PyTables implements several classes to represent the different nodes in the
object tree. They are named File, Group, Leaf, Table, Array, CArray, EArray,
VLArray and UnImplemented. Another one allows the user to complement the
information on these different objects; its name is AttributeSet. Finally,
another important class called IsDescription allows to build a Table record
description by declaring a subclass of it. Many other classes are defined in
PyTables, but they can be regarded as helpers whose goal is mainly to declare
the *data type properties* of the different first class objects and will be
described at the end of this chapter as well.

An important function, called open_file is responsible to create, open or append
to files. In addition, a few utility functions are defined to guess if the user
supplied file is a *PyTables* or *HDF5* file. These are called is_pytables_file()
and is_hdf5_file(), respectively. There exists also a function called
which_lib_version() that informs about the versions of the underlying C libraries
(for example, HDF5 or Zlib) and another called print_versions() that prints all
the versions of the software that PyTables relies on. Finally, test() lets you
run the complete test suite from a Python console interactively.

.. toctree::
    :maxdepth: 2

    libref/top_level
    libref/file_class
    libref/hierarchy_classes
    libref/structured_storage
    libref/homogenous_storage
    libref/link_classes
    libref/declarative_classes
    libref/helper_classes
    libref/expr_class
    libref/filenode_classes
Utilities
=========
PyTables comes with a couple of utilities that make the life easier to the
user. One is called ptdump and lets you see the contents of a PyTables file
(or generic HDF5 file, if supported). The other one is named ptrepack that
allows to (recursively) copy sub-hierarchies of objects present in a file
into another one, changing, if desired, some of the filters applied to the
leaves during the copy process.

Normally, these utilities will be installed somewhere in your PATH during the
process of installation of the PyTables package, so that you can invoke them
from any place in your file system after the installation has successfully
finished.


ptdump
------
As has been said before, ptdump utility allows you look into the contents of
your PyTables files. It lets you see not only the data but also the metadata
(that is, the *structure* and additional information in the form of
*attributes*).

Usage
~~~~~
For instructions on how to use it, just pass the -h flag to the command:

.. code-block:: bash

    $ ptdump -h

to see the message usage:

.. code-block:: bash

    usage: ptdump [-h] [-v] [-d] [-a] [-s] [-c] [-i] [-R RANGE]
                  filename[:nodepath]

    The ptdump utility allows you look into the contents of your PyTables files.
    It lets you see not only the data but also the metadata (that is, the
    *structure* and additional information in the form of *attributes*).

    positional arguments:
      filename[:nodepath]   name of the HDF5 file to dump

    optional arguments:
      -h, --help            show this help message and exit
      -v, --verbose         dump more metainformation on nodes
      -d, --dump            dump data information on leaves
      -a, --showattrs       show attributes in nodes (only useful when -v or -d
                            are active)
      -s, --sort            sort output by node name
      -c, --colinfo         show info of columns in tables (only useful when -v or
                            -d are active)
      -i, --idxinfo         show info of indexed columns (only useful when -v or
                            -d are active)
      -R RANGE, --range RANGE
                            select a RANGE of rows (in the form "start,stop,step")
                            during the copy of *all* the leaves. Default values
                            are "None,None,1", which means a copy of all the rows.

Read on for a brief introduction to this utility.


A small tutorial on ptdump
~~~~~~~~~~~~~~~~~~~~~~~~~~
Let's suppose that we want to know only the *structure* of a file. In order
to do that, just don't pass any flag, just the file as parameter.

.. code-block:: bash

    $ ptdump vlarray1.h5
    / (RootGroup) ''
    /vlarray1 (VLArray(3,), shuffle, zlib(1)) 'ragged array of ints'
    /vlarray2 (VLArray(3,), shuffle, zlib(1)) 'ragged array of strings'

we can see that the file contains just a leaf object called vlarray1, that is
an instance of VLArray, has 4 rows, and two filters has been used in order to
create it: shuffle and zlib (with a compression level of 1).

Let's say we want more meta-information. Just add the -v (verbose) flag:

.. code-block:: bash

    $ ptdump -v vlarray1.h5
    / (RootGroup) ''
    /vlarray1 (VLArray(3,), shuffle, zlib(1)) 'ragged array of ints'
      atom = Int32Atom(shape=(), dflt=0)
      byteorder = 'little'
      nrows = 3
      flavor = 'numpy'
    /vlarray2 (VLArray(3,), shuffle, zlib(1)) 'ragged array of strings'
      atom = StringAtom(itemsize=2, shape=(), dflt='')
      byteorder = 'irrelevant'
      nrows = 3
      flavor = 'python'

so we can see more info about the atoms that are the components of the
vlarray1 dataset, i.e. they are scalars of type Int32 and with NumPy
*flavor*.

If we want information about the attributes on the nodes, we must add the -a
flag:

.. code-block:: bash

    $ ptdump -va vlarray1.h5
    / (RootGroup) ''
      /._v_attrs (AttributeSet), 4 attributes:
       [CLASS := 'GROUP',
        PYTABLES_FORMAT_VERSION := '2.0',
        TITLE := '',
        VERSION := '1.0']
    /vlarray1 (VLArray(3,), shuffle, zlib(1)) 'ragged array of ints'
      atom = Int32Atom(shape=(), dflt=0)
      byteorder = 'little'
      nrows = 3
      flavor = 'numpy'
      /vlarray1._v_attrs (AttributeSet), 3 attributes:
       [CLASS := 'VLARRAY',
        TITLE := 'ragged array of ints',
        VERSION := '1.3']
    /vlarray2 (VLArray(3,), shuffle, zlib(1)) 'ragged array of strings'
      atom = StringAtom(itemsize=2, shape=(), dflt='')
      byteorder = 'irrelevant'
      nrows = 3
      flavor = 'python'
      /vlarray2._v_attrs (AttributeSet), 4 attributes:
       [CLASS := 'VLARRAY',
        FLAVOR := 'python',
        TITLE := 'ragged array of strings',
        VERSION := '1.3']


Let's have a look at the real data:

.. code-block:: bash

    $ ptdump -d vlarray1.h5
    / (RootGroup) ''
    /vlarray1 (VLArray(3,), shuffle, zlib(1)) 'ragged array of ints'
      Data dump:
    [0] [5 6]
    [1] [5 6 7]
    [2] [5 6 9 8]
    /vlarray2 (VLArray(3,), shuffle, zlib(1)) 'ragged array of strings'
      Data dump:
    [0] ['5', '66']
    [1] ['5', '6', '77']
    [2] ['5', '6', '9', '88']

We see here a data dump of the 4 rows in vlarray1 object, in the form of a
list. Because the object is a VLA, we see a different number of integers on
each row.

Say that we are interested only on a specific *row range* of the /vlarray1
object:

.. code-block:: bash

    ptdump -R2,3 -d vlarray1.h5:/vlarray1
    /vlarray1 (VLArray(3,), shuffle, zlib(1)) 'ragged array of ints'
      Data dump:
    [2] [5 6 9 8]

Here, we have specified the range of rows between 2 and 4 (the upper limit
excluded, as usual in Python). See how we have selected only the /vlarray1
object for doing the dump (vlarray1.h5:/vlarray1).

Finally, you can mix several information at once:

.. code-block:: bash

    $ ptdump -R2,3 -vad vlarray1.h5:/vlarray1
    /vlarray1 (VLArray(3,), shuffle, zlib(1)) 'ragged array of ints'
      atom = Int32Atom(shape=(), dflt=0)
      byteorder = 'little'
      nrows = 3
      flavor = 'numpy'
      /vlarray1._v_attrs (AttributeSet), 3 attributes:
       [CLASS := 'VLARRAY',
        TITLE := 'ragged array of ints',
        VERSION := '1.3']
      Data dump:
    [2] [5 6 9 8]


.. _ptrepackDescr:

ptrepack
--------
This utility is a very powerful one and lets you copy any leaf, group or
complete subtree into another file. During the copy process you are allowed
to change the filter properties if you want so. Also, in the case of
duplicated pathnames, you can decide if you want to overwrite already
existing nodes on the destination file. Generally speaking, ptrepack can be
useful in may situations, like replicating a subtree in another file, change
the filters in objects and see how affect this to the compression degree or
I/O performance, consolidating specific data in repositories or even
*importing* generic HDF5 files and create true PyTables counterparts.


Usage
~~~~~
For instructions on how to use it, just pass the -h flag to the command:

.. code-block:: bash

    $ ptrepack -h

to see the message usage:

.. code-block:: bash

    usage: ptrepack [-h] [-v] [-o] [-R RANGE] [--non-recursive]
                    [--dest-title TITLE] [--dont-create-sysattrs]
                    [--dont-copy-userattrs] [--overwrite-nodes]
                    [--complevel COMPLEVEL]
                    [--complib {zlib,lzo,bzip2,blosc,blosc:blosclz,blosc:lz4,blosc:lz4hc,blosc:snappy,blosc:zlib,blosc:zstd}]
                    [--shuffle {0,1}] [--bitshuffle {0,1}] [--fletcher32 {0,1}]
                    [--keep-source-filters] [--chunkshape CHUNKSHAPE]
                    [--upgrade-flavors] [--dont-regenerate-old-indexes]
                    [--sortby COLUMN] [--checkCSI] [--propindexes]
                    sourcefile:sourcegroup destfile:destgroup

    This utility is very powerful and lets you copy any leaf, group or complete
    subtree into another file. During the copy process you are allowed to change
    the filter properties if you want so. Also, in the case of duplicated
    pathnames, you can decide if you want to overwrite already existing nodes on
    the destination file. Generally speaking, ptrepack can be useful in may
    situations, like replicating a subtree in another file, change the filters in
    objects and see how affect this to the compression degree or I/O performance,
    consolidating specific data in repositories or even *importing* generic HDF5
    files and create true PyTables counterparts.

    positional arguments:
      sourcefile:sourcegroup
                            source file/group
      destfile:destgroup    destination file/group

    optional arguments:
      -h, --help            show this help message and exit
      -v, --verbose         show verbose information
      -o, --overwrite       overwrite destination file
      -R RANGE, --range RANGE
                            select a RANGE of rows (in the form "start,stop,step")
                            during the copy of *all* the leaves. Default values
                            are "None,None,1", which means a copy of all the rows.
      --non-recursive       do not do a recursive copy. Default is to do it
      --dest-title TITLE    title for the new file (if not specified, the source
                            is copied)
      --dont-create-sysattrs
                            do not create sys attrs (default is to do it)
      --dont-copy-userattrs
                            do not copy the user attrs (default is to do it)
      --overwrite-nodes     overwrite destination nodes if they exist. Default is
                            to not overwrite them
      --complevel COMPLEVEL
                            set a compression level (0 for no compression, which
                            is the default)
      --complib {zlib,lzo,bzip2,blosc,blosc:blosclz,blosc:lz4,blosc:lz4hc,blosc:snappy,blosc:zlib,blosc:zstd}
                            set the compression library to be used during the
                            copy. Defaults to zlib
      --shuffle {0,1}       activate or not the shuffle filter (default is active
                            if complevel > 0)
      --bitshuffle {0,1}    activate or not the bitshuffle filter (not active by
                            default)
      --fletcher32 {0,1}    whether to activate or not the fletcher32 filter (not
                            active by default)
      --keep-source-filters
                            use the original filters in source files. The default
                            is not doing that if any of --complevel, --complib,
                            --shuffle --bitshuffle or --fletcher32 option is
                            specified
      --chunkshape CHUNKSHAPE
                            set a chunkshape. Possible options are: "keep" |
                            "auto" | int | tuple. A value of "auto" computes a
                            sensible value for the chunkshape of the leaves
                            copied. The default is to "keep" the original value
      --upgrade-flavors     when repacking PyTables 1.x or PyTables 2.x files, the
                            flavor of leaves will be unset. With this, such a
                            leaves will be serialized as objects with the internal
                            flavor ('numpy' for 3.x series)
      --dont-regenerate-old-indexes
                            disable regenerating old indexes. The default is to
                            regenerate old indexes as they are found
      --sortby COLUMN       do a table copy sorted by the index in "column". For
                            reversing the order, use a negative value in the
                            "step" part of "RANGE" (see "-r" flag). Only applies
                            to table objects
      --checkCSI            force the check for a CSI index for the --sortby
                            column
      --propindexes         propagate the indexes existing in original tables. The
                            default is to not propagate them. Only applies to
                            table objects
      --dont-allow-padding  remove the possible padding in compound types in
                            source files. The default is to propagate it. Only
                            applies to table objects


Read on for a brief introduction to this utility.

A small tutorial on ptrepack
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Imagine that we have ended the tutorial 1 (see the output of
examples/tutorial1-1.py), and we want to copy our reduced data (i.e. those
datasets that hangs from the /column group) to another file. First, let's
remember the content of the examples/tutorial1.h5:

.. code-block:: bash

    $ ptdump tutorial1.h5
    / (RootGroup) 'Test file'
    /columns (Group) 'Pressure and Name'
    /columns/name (Array(3,)) 'Name column selection'
    /columns/pressure (Array(3,)) 'Pressure column selection'
    /detector (Group) 'Detector information'
    /detector/readout (Table(10,)) 'Readout example'

Now, copy the /columns to other non-existing file. That's easy:

.. code-block:: bash

    $ ptrepack tutorial1.h5:/columns reduced.h5

That's all. Let's see the contents of the newly created reduced.h5 file:

.. code-block:: bash

    $ ptdump reduced.h5
    / (RootGroup) ''
    /name (Array(3,)) 'Name column selection'
    /pressure (Array(3,)) 'Pressure column selection'

so, you have copied the children of /columns group into the *root* of the
reduced.h5 file.

Now, you suddenly realized that what you intended to do was to copy all the
hierarchy, the group /columns itself included. You can do that by just
specifying the destination group:

.. code-block:: bash

    $ ptrepack tutorial1.h5:/columns reduced.h5:/columns
    $ ptdump reduced.h5
    / (RootGroup) ''
    /name (Array(3,)) 'Name column selection'
    /pressure (Array(3,)) 'Pressure column selection'
    /columns (Group) ''
    /columns/name (Array(3,)) 'Name column selection'
    /columns/pressure (Array(3,)) 'Pressure column selection'

OK. Much better. But you want to get rid of the existing nodes on the new
file. You can achieve this by adding the -o flag:

.. code-block:: bash

    $ ptrepack -o tutorial1.h5:/columns reduced.h5:/columns
    $ ptdump reduced.h5
    / (RootGroup) ''
    /columns (Group) ''
    /columns/name (Array(3,)) 'Name column selection'
    /columns/pressure (Array(3,)) 'Pressure column selection'

where you can see how the old contents of the reduced.h5 file has been
overwritten.

You can copy just one single node in the repacking operation and change its
name in destination:

.. code-block:: bash

    $ ptrepack tutorial1.h5:/detector/readout reduced.h5:/rawdata
    $ ptdump reduced.h5
    / (RootGroup) ''
    /rawdata (Table(10,)) 'Readout example'
    /columns (Group) ''
    /columns/name (Array(3,)) 'Name column selection'
    /columns/pressure (Array(3,)) 'Pressure column selection'

where the /detector/readout has been copied to /rawdata in destination.

We can change the filter properties as well:

.. code-block:: bash

    $ ptrepack --complevel=1 tutorial1.h5:/detector/readout reduced.h5:/rawdata
    Problems doing the copy from 'tutorial1.h5:/detector/readout' to 'reduced.h5:/rawdata'
    The error was --> tables.exceptions.NodeError: destination group \``/\`` already has a node named \``rawdata``; you may want to use the \``overwrite`` argument
    The destination file looks like:
    / (RootGroup) ''
    /rawdata (Table(10,)) 'Readout example'
    /columns (Group) ''
    /columns/name (Array(3,)) 'Name column selection'
    /columns/pressure (Array(3,)) 'Pressure column selection'
    Traceback (most recent call last):
      File "utils/ptrepack", line 3, in ?
        main()
      File ".../tables/scripts/ptrepack.py", line 349, in main
        stats = stats, start = start, stop = stop, step = step)
      File ".../tables/scripts/ptrepack.py", line 107, in copy_leaf
        raise RuntimeError, "Please check that the node names are not
        duplicated in destination, and if so, add the --overwrite-nodes flag
        if desired."
    RuntimeError: Please check that the node names are not duplicated in
    destination, and if so, add the --overwrite-nodes flag if desired.

Ooops! We ran into problems: we forgot that the /rawdata pathname already
existed in destination file. Let's add the --overwrite-nodes, as the verbose
error suggested:

.. code-block:: bash

    $ ptrepack --overwrite-nodes --complevel=1 tutorial1.h5:/detector/readout
    reduced.h5:/rawdata
    $ ptdump reduced.h5
    / (RootGroup) ''
    /rawdata (Table(10,), shuffle, zlib(1)) 'Readout example'
    /columns (Group) ''
    /columns/name (Array(3,)) 'Name column selection'
    /columns/pressure (Array(3,)) 'Pressure column selection'

you can check how the filter properties has been changed for the /rawdata
table. Check as the other nodes still exists.

Finally, let's copy a *slice* of the readout table in origin to destination,
under a new group called /slices and with the name, for example, aslice:

.. code-block:: bash

    $ ptrepack -R1,8,3 tutorial1.h5:/detector/readout reduced.h5:/slices/aslice
    $ ptdump reduced.h5
    / (RootGroup) ''
    /rawdata (Table(10,), shuffle, zlib(1)) 'Readout example'
    /columns (Group) ''
    /columns/name (Array(3,)) 'Name column selection'
    /columns/pressure (Array(3,)) 'Pressure column selection'
    /slices (Group) ''
    /slices/aslice (Table(3,)) 'Readout example'

note how only 3 rows of the original readout table has been copied to the new
aslice destination. Note as well how the previously nonexistent slices group
has been created in the same operation.



pt2to3
------

The PyTables 3.x series now follows `PEP 8`_ coding standard.  This makes
using PyTables more idiomatic with surrounding Python code that also adheres
to this standard.  The primary way that the 2.x series was *not* PEP 8
compliant was with respect to variable naming conventions.  Approximately 450
API variables were identified and updated for PyTables 3.x.

To ease migration, PyTables ships with a new ``pt2to3`` command line tool.
This tool will run over a file and replace any instances of the old variable
names with the 3.x version of the name.  This tool covers the overwhelming
majority of cases was used to transition the PyTables code base itself!  However,
it may also accidentally also pick up variable names in 3rd party codes that
have *exactly* the same name as a PyTables' variable.  This is because ``pt2to3``
was implemented using regular expressions rather than a fancier AST-based
method. By using regexes, ``pt2to3`` works on Python and Cython code.


``pt2to3`` **help:**

.. code-block:: bash

    usage: pt2to3 [-h] [-r] [-p] [-o OUTPUT] [-i] filename

    PyTables 2.x -> 3.x API transition tool This tool displays to standard out, so
    it is common to pipe this to another file: $ pt2to3 oldfile.py > newfile.py

    positional arguments:
      filename              path to input file.

    optional arguments:
      -h, --help            show this help message and exit
      -r, --reverse         reverts changes, going from 3.x -> 2.x.
      -p, --no-ignore-previous
                            ignores previous_api() calls.
      -o OUTPUT             output file to write to.
      -i, --inplace         overwrites the file in-place.

Note that ``pt2to3`` only works on a single file, not a a directory.  However,
a simple BASH script may be written to run ``pt2to3`` over an entire directory
and all sub-directories:

.. code-block:: bash

    #!/bin/bash
    for f in $(find .)
    do
        echo $f
        pt2to3 $f > temp.txt
        mv temp.txt $f
    done

.. _PEP 8: http://www.python.org/dev/peps/pep-0008/
Tutorials
=========
.. epigraph::

    Seràs la clau que obre tots els panys,
    seràs la llum, la llum il.limitada,
    seràs confí on l'aurora comença,
    seràs forment, escala il.luminada!

    -- Lyrics: Vicent Andrés i Estellés. Music: Ovidi Montllor, Toti Soler, M'aclame a tu


This chapter consists of a series of simple yet comprehensive
tutorials that will enable you to understand PyTables' main features. If
you would like more information about some particular instance variable,
global function, or method, look at the doc strings or go to the library
reference in :ref:`library_reference`. If you are reading this in PDF or HTML
formats, follow the corresponding hyperlink near each newly introduced
entity.

Please note that throughout this document the terms *column* and *field*
will be used interchangeably, as will the terms *row* and *record*.

.. currentmodule:: tables

Getting started
---------------
In this section, we will see how to define our own records in Python and save
collections of them (i.e. a *table*) into a file. Then we will select some of
the data in the table using Python cuts and create NumPy arrays to store this
selection as separate objects in a tree.

In *examples/tutorial1-1.py* you will find the working version of all the
code in this section. Nonetheless, this tutorial series has been written to
allow you reproduce it in a Python interactive console. I encourage you to do
parallel testing and inspect the created objects (variables, docs, children
objects, etc.) during the course of the tutorial!


Importing tables objects
~~~~~~~~~~~~~~~~~~~~~~~~
Before starting you need to import the public objects in the tables package.
You normally do that by executing::

    >>> import tables

This is the recommended way to import tables if you don't want to pollute
your namespace. However, PyTables has a contained set of first-level
primitives, so you may consider using the alternative::

    >>> from tables import *

If you are going to work with NumPy arrays (and normally, you will) you will
also need to import functions from the numpy package. So most PyTables
programs begin with::

    >>> import tables   # but in this tutorial we use "from tables import \*"
    >>> import numpy as np


Declaring a Column Descriptor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now, imagine that we have a particle detector and we want to create a table
object in order to save data retrieved from it. You need first to define the
table, the number of columns it has, what kind of object is contained in each
column, and so on.

Our particle detector has a TDC (Time to Digital Converter) counter with a
dynamic range of 8 bits and an ADC (Analogical to Digital Converter) with a
range of 16 bits. For these values, we will define 2 fields in our record
object called TDCcount and ADCcount. We also want to save the grid position
in which the particle has been detected, so we will add two new fields called
grid_i and grid_j. Our instrumentation also can obtain the pressure and
energy of the particle. The resolution of the pressure-gauge allows us to use
a single-precision float to store pressure readings, while the energy value
will need a double-precision float. Finally, to track the particle we want to
assign it a name to identify the kind of the particle it is and a unique
numeric identifier. So we will add two more fields: name will be a string of
up to 16 characters, and idnumber will be an integer of 64 bits (to allow us
to store records for extremely large numbers of particles).

Having determined our columns and their types, we can now declare a new
Particle class that will contain all this information::

    >>> from tables import *
    >>> class Particle(IsDescription):
    ...     name      = StringCol(16)   # 16-character String
    ...     idnumber  = Int64Col()      # Signed 64-bit integer
    ...     ADCcount  = UInt16Col()     # Unsigned short integer
    ...     TDCcount  = UInt8Col()      # unsigned byte
    ...     grid_i    = Int32Col()      # 32-bit integer
    ...     grid_j    = Int32Col()      # 32-bit integer
    ...     pressure  = Float32Col()    # float  (single-precision)
    ...     energy    = Float64Col()    # double (double-precision)
    >>>

This definition class is self-explanatory. Basically, you declare a class
variable for each field you need. As its value you assign an instance of the
appropriate Col subclass, according to the kind of column defined (the data
type, the length, the shape, etc). See the :ref:`ColClassDescr` for a
complete description of these subclasses. See also :ref:`datatypes` for a
list of data types supported by the Col constructor.

From now on, we can use Particle instances as a descriptor for our detector
data table. We will see later on how to pass this object to construct the
table. But first, we must create a file where all the actual data pushed into
our table will be saved.


Creating a PyTables file from scratch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the top-level :func:`open_file` function to create a PyTables file::

    >>> h5file = open_file("tutorial1.h5", mode="w", title="Test file")

:func:`open_file` is one of the objects imported by the
```from tables import *``` statement. Here, we are saying that we want to
create a new file in the current working directory called "tutorial1.h5" in
"w"rite mode and with a descriptive title string ("Test file").
This function attempts to open the file, and if successful, returns the File
(see :ref:`FileClassDescr`) object instance h5file. The root of the object
tree is specified in the instance's root attribute.


Creating a new group
~~~~~~~~~~~~~~~~~~~~
Now, to better organize our data, we will create a group called *detector*
that branches from the root node. We will save our particle data table in
this group::

    >>> group = h5file.create_group("/", 'detector', 'Detector information')

Here, we have taken the File instance h5file and invoked its
:meth:`File.create_group` method to create a  new group called *detector*
branching from "*/*" (another way to refer to the h5file.root object we
mentioned above). This will create a new Group (see :ref:`GroupClassDescr`)
object instance that will be assigned to the variable group.


Creating a new table
~~~~~~~~~~~~~~~~~~~~
Let's now create a Table (see :ref:`TableClassDescr`) object as a branch off
the newly-created group. We do that by calling the :meth:`File.create_table`
method of the h5file object::

    >>> table = h5file.create_table(group, 'readout', Particle, "Readout example")

We create the Table instance under group. We assign this table the node name
"*readout*". The Particle class declared before is the *description*
parameter (to define the columns of the table) and finally we set
"*Readout example*" as the Table title. With all this information, a new
Table instance is created and assigned to the variable *table*.

If you are curious about how the object tree looks right now, simply print
the File instance variable *h5file*, and examine the output::

    >>> print(h5file)
    tutorial1.h5 (File) 'Test file'
    Last modif.: 'Wed Mar  7 11:06:12 2007'
    Object Tree:
    / (RootGroup) 'Test file'
    /detector (Group) 'Detector information'
    /detector/readout (Table(0,)) 'Readout example'

As you can see, a dump of the object tree is displayed. It's easy to see the
Group and Table objects we have just created. If you want more information,
just type the variable containing the File instance::

    >>> h5file
    File(filename='tutorial1.h5', title='Test file', mode='w', root_uep='/', filters=Filters(complevel=0, shuffle=False, bitshuffle=False, fletcher32=False))
    / (RootGroup) 'Test file'
    /detector (Group) 'Detector information'
    /detector/readout (Table(0,)) 'Readout example'
    description := {
      "ADCcount": UInt16Col(shape=(), dflt=0, pos=0),
      "TDCcount": UInt8Col(shape=(), dflt=0, pos=1),
      "energy": Float64Col(shape=(), dflt=0.0, pos=2),
      "grid_i": Int32Col(shape=(), dflt=0, pos=3),
      "grid_j": Int32Col(shape=(), dflt=0, pos=4),
      "idnumber": Int64Col(shape=(), dflt=0, pos=5),
      "name": StringCol(itemsize=16, shape=(), dflt='', pos=6),
      "pressure": Float32Col(shape=(), dflt=0.0, pos=7)}
      byteorder := 'little'
      chunkshape := (87,)

More detailed information is displayed about each object in the tree. Note
how Particle, our table descriptor class, is printed as part of the *readout*
table description information. In general, you can obtain much more
information about the objects and their children by just printing them. That
introspection capability is very useful, and I recommend that you use it
extensively.

The time has come to fill this table with some values. First we will get a
pointer to the Row (see :ref:`RowClassDescr`) instance of this table
instance::

    >>> particle = table.row

The row attribute of table points to the Row instance that will be used to
write data rows into the table. We write data simply by assigning the Row
instance the values for each row as if it were a dictionary (although it is
actually an *extension class*), using the column names as keys.

Below is an example of how to write rows::

    >>> for i in range(10):
    ...     particle['name']  = f'Particle: {i:6d}'
    ...     particle['TDCcount'] = i % 256
    ...     particle['ADCcount'] = (i * 256) % (1 << 16)
    ...     particle['grid_i'] = i
    ...     particle['grid_j'] = 10 - i
    ...     particle['pressure'] = float(i*i)
    ...     particle['energy'] = float(particle['pressure'] ** 4)
    ...     particle['idnumber'] = i * (2 ** 34)
    ...     # Insert a new particle record
    ...     particle.append()
    >>>

This code should be easy to understand. The lines inside the loop just assign
values to the different columns in the Row instance particle (see
:ref:`RowClassDescr`). A call to its append() method writes this information
to the table I/O buffer.

After we have processed all our data, we should flush the table's I/O buffer
if we want to write all this data to disk. We achieve that by calling the
table.flush() method::

    >>> table.flush()

Remember, flushing a table is a *very important* step as it will not only
help to maintain the integrity of your file, but also will free valuable
memory resources (i.e. internal buffers) that your program may need for other
things.


.. _readingAndSelectingUsage:

Reading (and selecting) data in a table
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ok. We have our data on disk, and now we need to access it and select from
specific columns the values we are interested in. See the example below::

    >>> table = h5file.root.detector.readout
    >>> pressure = [x['pressure'] for x in table.iterrows() if x['TDCcount'] > 3 and 20 <= x['pressure'] < 50]
    >>> pressure
    [25.0, 36.0, 49.0]

The first line creates a "shortcut" to the *readout* table deeper on the
object tree. As you can see, we use the *natural naming* schema to access it.
We also could have used the h5file.get_node() method, as we will do later on.

You will recognize the last two lines as a Python list comprehension.
It loops over the rows in *table* as they are provided by the
:meth:`Table.iterrows` iterator. The iterator returns values until all the
data in table is exhausted. These rows are filtered using the expression::

    x['TDCcount'] > 3 and 20 <= x['pressure'] < 50

So, we are selecting the values of the pressure column from filtered records
to create the final list and assign it to pressure variable.

We could have used a normal for loop to accomplish the same purpose, but I
find comprehension syntax to be more compact and elegant.

PyTables do offer other, more powerful ways of performing selections which
may be more suitable if you have very large tables or if you need very high
query speeds. They are called *in-kernel* and *indexed* queries, and you can
use them through :meth:`Table.where` and other related methods.

Let's use an in-kernel selection to query the name column for the same set of
cuts::

    >>> names = [ x['name'] for x in table.where("""(TDCcount > 3) & (20 <= pressure) & (pressure < 50)""") ]
    >>> names
    ['Particle:      5', 'Particle:      6', 'Particle:      7']

In-kernel and indexed queries are not only much faster, but as you can see,
they also look more compact, and are among the greatest features for
PyTables, so be sure that you use them a lot. See :ref:`condition_syntax` and
:ref:`searchOptim` for more information on in-kernel and indexed selections.

.. note::

    A special care should be taken when the query condition includes
    string literals.  Indeed Python 2 string literals are string of
    bytes while Python 3 strings are unicode objects.

    With reference to the above definition of :class:`Particle` it has to be
    noted that the type of the "name" column do not change depending on the
    Python version used (of course).
    It always corresponds to strings of bytes.

    Any condition involving the "name" column should be written using the
    appropriate type for string literals in order to avoid
    :exc:`TypeError`\ s.

    Suppose one wants to get rows corresponding to specific particle names.

    The code below will work fine in Python 2 but will fail with a
    :exc:`TypeError` in Python 3::

        >>> condition = '(name == "Particle:      5") | (name == "Particle:      7")'
        >>> for record in table.where(condition):  # TypeError in Python3
        ...     # do something with "record"

    The reason is that in Python 3 "condition" implies a comparison
    between a string of bytes ("name" column contents) and an unicode
    literals.

    The correct way to write the condition is::

        >>> condition = '(name == b"Particle:      5") | (name == b"Particle:      7")'

That's enough about selections for now. The next section will show you how to
save these selected results to a file.


Creating new array objects
~~~~~~~~~~~~~~~~~~~~~~~~~~
In order to separate the selected data from the mass of detector data, we
will create a new group columns branching off the root group. Afterwards,
under this group, we will create two arrays that will contain the selected
data. First, we create the group::

    >>> gcolumns = h5file.create_group(h5file.root, "columns", "Pressure and Name")

Note that this time we have specified the first parameter using *natural
naming* (h5file.root) instead of with an absolute path string ("/").

Now, create the first of the two Array objects we've just mentioned::

    >>> h5file.create_array(gcolumns, 'pressure', np.array(pressure), "Pressure column selection")
    /columns/pressure (Array(3,)) 'Pressure column selection'
      atom := Float64Atom(shape=(), dflt=0.0)
      maindim := 0
      flavor := 'numpy'
      byteorder := 'little'
      chunkshape := None

We already know the first two parameters of the :meth:`File.create_array`
methods (these are the same as the first two in create_table): they are the
parent group *where* Array will be created and the Array instance *name*.
The third parameter is the *object* we want to save to disk. In this case, it
is a NumPy array that is built from the selection list we created before.
The fourth parameter is the *title*.

Now, we will save the second array. It contains the list of strings we
selected before: we save this object as-is, with no further conversion::

    >>> h5file.create_array(gcolumns, 'name', names, "Name column selection")
    /columns/name (Array(3,)) 'Name column selection'
      atom := StringAtom(itemsize=16, shape=(), dflt='')
      maindim := 0
      flavor := 'python'
      byteorder := 'irrelevant'
      chunkshape := None

As you can see, :meth:`File.create_array` accepts *names* (which is a regular
Python list) as an *object* parameter. Actually, it accepts a variety of
different regular objects (see :func:`create_array`) as parameters. The flavor
attribute (see the output above) saves the original kind of object that was
saved. Based on this *flavor*, PyTables will be able to retrieve exactly the
same object from disk later on.

Note that in these examples, the create_array method returns an Array instance
that is not assigned to any variable. Don't worry, this is intentional to
show the kind of object we have created by displaying its representation. The
Array objects have been attached to the object tree and saved to disk, as you
can see if you print the complete object tree::

    >>> print(h5file)
    tutorial1.h5 (File) 'Test file'
    Last modif.: 'Wed Mar  7 19:40:44 2007'
    Object Tree:
    / (RootGroup) 'Test file'
    /columns (Group) 'Pressure and Name'
    /columns/name (Array(3,)) 'Name column selection'
    /columns/pressure (Array(3,)) 'Pressure column selection'
    /detector (Group) 'Detector information'
    /detector/readout (Table(10,)) 'Readout example'


Closing the file and looking at its content
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To finish this first tutorial, we use the close method of the h5file File
object to close the file before exiting Python::

    >>> h5file.close()
    >>> ^D
    $

You have now created your first PyTables file with a table and two arrays.
You can examine it with any generic HDF5 tool, such as h5dump or h5ls. Here is
what the tutorial1.h5 looks like when read with the h5ls program.

.. code-block:: bash

    $ h5ls -rd tutorial1.h5
    /columns                 Group
    /columns/name            Dataset {3}
        Data:
            (0) "Particle:      5", "Particle:      6", "Particle:      7"
    /columns/pressure        Dataset {3}
        Data:
            (0) 25, 36, 49
    /detector                Group
    /detector/readout        Dataset {10/Inf}
        Data:
            (0) {0, 0, 0, 0, 10, 0, "Particle:      0", 0},
            (1) {256, 1, 1, 1, 9, 17179869184, "Particle:      1", 1},
            (2) {512, 2, 256, 2, 8, 34359738368, "Particle:      2", 4},
            (3) {768, 3, 6561, 3, 7, 51539607552, "Particle:      3", 9},
            (4) {1024, 4, 65536, 4, 6, 68719476736, "Particle:      4", 16},
            (5) {1280, 5, 390625, 5, 5, 85899345920, "Particle:      5", 25},
            (6) {1536, 6, 1679616, 6, 4, 103079215104, "Particle:      6", 36},
            (7) {1792, 7, 5764801, 7, 3, 120259084288, "Particle:      7", 49},
            (8) {2048, 8, 16777216, 8, 2, 137438953472, "Particle:      8", 64},
            (9) {2304, 9, 43046721, 9, 1, 154618822656, "Particle:      9", 81}

Here's the output as displayed by the "ptdump" PyTables utility (located in
utils/ directory).

.. code-block:: bash

    $ ptdump tutorial1.h5
    / (RootGroup) 'Test file'
    /columns (Group) 'Pressure and Name'
    /columns/name (Array(3,)) 'Name column selection'
    /columns/pressure (Array(3,)) 'Pressure column selection'
    /detector (Group) 'Detector information'
    /detector/readout (Table(10,)) 'Readout example'

You can pass the `-v` or `-d` options to ptdump if you want
more verbosity. Try them out!

Also, in :ref:`Figure 1 <tutorial1-1-tableview>`, you can admire how the
tutorial1.h5 looks like using the `ViTables <http://vitables.org>`_ graphical
interface.

.. _tutorial1-1-tableview:

.. figure:: images/tutorial1-1-tableview.png
    :align: center

    **Figure 1. The initial version of the data file for tutorial 1, with a
    view of the data objects.**


Browsing the *object tree*
--------------------------
In this section, we will learn how to browse the tree and retrieve data and
also meta-information about the actual data.

In *examples/tutorial1-2.py* you will find the working version of all the
code in this section. As before, you are encouraged to use a python shell and
inspect the object tree during the course of the tutorial.


Traversing the object tree
~~~~~~~~~~~~~~~~~~~~~~~~~~
Let's start by opening the file we created in last tutorial section::

    >>> h5file = open_file("tutorial1.h5", "a")

This time, we have opened the file in "a"ppend mode. We use this mode to add
more information to the file.

PyTables, following the Python tradition, offers powerful introspection
capabilities, i.e. you can easily ask information about any component of the
object tree as well as search the tree.

To start with, you can get a preliminary overview of the object tree by
simply printing the existing File instance::

    >>> print(h5file)
    tutorial1.h5 (File) 'Test file'
    Last modif.: 'Wed Mar  7 19:50:57 2007'
    Object Tree:
    / (RootGroup) 'Test file'
    /columns (Group) 'Pressure and Name'
    /columns/name (Array(3,)) 'Name column selection'
    /columns/pressure (Array(3,)) 'Pressure column selection'
    /detector (Group) 'Detector information'
    /detector/readout (Table(10,)) 'Readout example'

It looks like all of our objects are there. Now let's make use of the File
iterator to see how to list all the nodes in the object tree::

    >>> for node in h5file:
    ...     print(node)
    / (RootGroup) 'Test file'
    /columns (Group) 'Pressure and Name'
    /detector (Group) 'Detector information'
    /columns/name (Array(3,)) 'Name column selection'
    /columns/pressure (Array(3,)) 'Pressure column selection'
    /detector/readout (Table(10,)) 'Readout example'

We can use the :meth:`File.walk_groups` method of the File class to list only
the *groups* on tree::

    >>> for group in h5file.walk_groups():
    ...     print(group)
    / (RootGroup) 'Test file'
    /columns (Group) 'Pressure and Name'
    /detector (Group) 'Detector information'

Note that :meth:`File.walk_groups` actually returns an *iterator*, not a list
of objects. Using this iterator with the list_nodes() method is a powerful
combination. Let's see an example listing of all the arrays in the tree::

    >>> for group in h5file.walk_groups("/"):
    ...     for array in h5file.list_nodes(group, classname='Array'):
    ...         print(array)
    /columns/name (Array(3,)) 'Name column selection'
    /columns/pressure (Array(3,)) 'Pressure column selection'

:meth:`File.list_nodes` returns a list containing all the nodes hanging off a
specific Group. If the *classname* keyword is specified, the method will
filter out all instances which are not descendants of the class. We have
asked for only Array instances. There exist also an iterator counterpart
called :meth:`File.iter_nodes` that might be handy is some situations, like
for example when dealing with groups with a large number of nodes behind it.

We can combine both calls by using the :meth:`File.walk_nodes` special method
of the File object. For example::

    >>> for array in h5file.walk_nodes("/", "Array"):
    ...     print(array)
    /columns/name (Array(3,)) 'Name column selection'
    /columns/pressure (Array(3,)) 'Pressure column selection'

This is a nice shortcut when working interactively.

Finally, we will list all the Leaf, i.e. Table and Array instances (see
:ref:`LeafClassDescr` for detailed information on Leaf class), in the
/detector group. Note that only one instance of the Table class (i.e.
readout) will be selected in this group (as should be the case)::

    >>> for leaf in h5file.root.detector._f_walknodes('Leaf'):
    ...     print(leaf)
    /detector/readout (Table(10,)) 'Readout example'

We have used a call to the :meth:`Group._f_walknodes` method, using the
*natural naming* path specification.

Of course you can do more sophisticated node selections using these powerful
methods. But first, let's take a look at some important PyTables object
instance variables.


Setting and getting user attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PyTables provides an easy and concise way to complement the meaning of your
node objects on the tree by using the AttributeSet class (see
:ref:`AttributeSetClassDescr`). You can access this object through the
standard attribute attrs in Leaf nodes and _v_attrs in Group nodes.

For example, let's imagine that we want to save the date indicating when the
data in /detector/readout table has been acquired, as well as the temperature
during the gathering process::

    >>> table = h5file.root.detector.readout
    >>> table.attrs.gath_date = "Wed, 06/12/2003 18:33"
    >>> table.attrs.temperature = 18.4
    >>> table.attrs.temp_scale = "Celsius"

Now, let's set a somewhat more complex attribute in the /detector group::

    >>> detector = h5file.root.detector
    >>> detector._v_attrs.stuff = [5, (2.3, 4.5), "Integer and tuple"]

Note how the AttributeSet instance is accessed with the _v_attrs attribute
because detector is a Group node. In general, you can save any standard
Python data structure as an attribute node. See :ref:`AttributeSetClassDescr`
for a more detailed explanation of how they are serialized for export to
disk.

Retrieving the attributes is equally simple::

    >>> table.attrs.gath_date
    'Wed, 06/12/2003 18:33'
    >>> table.attrs.temperature
    18.399999999999999
    >>> table.attrs.temp_scale
    'Celsius'
    >>> detector._v_attrs.stuff
    [5, (2.2999999999999998, 4.5), 'Integer and tuple']

You can probably guess how to delete attributes::

    >>> del table.attrs.gath_date

If you want to examine the current user attribute set of /detector/table, you
can print its representation (try hitting the TAB key twice if you are on a
Unix Python console with the rlcompleter module active)::

    >>> table.attrs
    /detector/readout._v_attrs (AttributeSet), 23 attributes:
        [CLASS := 'TABLE',
         FIELD_0_FILL := 0,
         FIELD_0_NAME := 'ADCcount',
         FIELD_1_FILL := 0,
         FIELD_1_NAME := 'TDCcount',
         FIELD_2_FILL := 0.0,
         FIELD_2_NAME := 'energy',
         FIELD_3_FILL := 0,
         FIELD_3_NAME := 'grid_i',
         FIELD_4_FILL := 0,
         FIELD_4_NAME := 'grid_j',
         FIELD_5_FILL := 0,
         FIELD_5_NAME := 'idnumber',
         FIELD_6_FILL := '',
         FIELD_6_NAME := 'name',
         FIELD_7_FILL := 0.0,
         FIELD_7_NAME := 'pressure',
         FLAVOR := 'numpy',
         NROWS := 10,
         TITLE := 'Readout example',
         VERSION := '2.6',
         temp_scale := 'Celsius',
         temperature := 18.399999999999999]

We've got all the attributes (including the *system* attributes). You can get
a list of *all* attributes or only the *user* or *system* attributes with the
_f_list() method::

    >>> print(table.attrs._f_list("all"))
    ['CLASS', 'FIELD_0_FILL', 'FIELD_0_NAME', 'FIELD_1_FILL', 'FIELD_1_NAME',
    'FIELD_2_FILL', 'FIELD_2_NAME', 'FIELD_3_FILL', 'FIELD_3_NAME', 'FIELD_4_FILL',
    'FIELD_4_NAME', 'FIELD_5_FILL', 'FIELD_5_NAME', 'FIELD_6_FILL', 'FIELD_6_NAME',
    'FIELD_7_FILL', 'FIELD_7_NAME', 'FLAVOR', 'NROWS', 'TITLE', 'VERSION',
    'temp_scale', 'temperature']
    >>> print(table.attrs._f_list("user"))
    ['temp_scale', 'temperature']
    >>> print(table.attrs._f_list("sys"))
    ['CLASS', 'FIELD_0_FILL', 'FIELD_0_NAME', 'FIELD_1_FILL', 'FIELD_1_NAME',
    'FIELD_2_FILL', 'FIELD_2_NAME', 'FIELD_3_FILL', 'FIELD_3_NAME', 'FIELD_4_FILL',
    'FIELD_4_NAME', 'FIELD_5_FILL', 'FIELD_5_NAME', 'FIELD_6_FILL', 'FIELD_6_NAME',
    'FIELD_7_FILL', 'FIELD_7_NAME', 'FLAVOR', 'NROWS', 'TITLE', 'VERSION']

You can also rename attributes::

    >>> table.attrs._f_rename("temp_scale","tempScale")
    >>> print(table.attrs._f_list())
    ['tempScale', 'temperature']

And, from PyTables 2.0 on, you are allowed also to set, delete or rename
system attributes::

    >>> table.attrs._f_rename("VERSION", "version")
    >>> table.attrs.VERSION
    Traceback (most recent call last):
        File "<stdin>", line 1, in <module>
        File "tables/attributeset.py", line 222, in __getattr__
            (name, self._v__nodepath)
    AttributeError: Attribute 'VERSION' does not exist in node: '/detector/readout'
    >>> table.attrs.version
    '2.6'

*Caveat emptor:* you must be careful when modifying system attributes because
you may end fooling PyTables and ultimately getting unwanted behaviour. Use
this only if you know what are you doing.

So, given the caveat above, we will proceed to restore the original name of
VERSION attribute::

    >>> table.attrs._f_rename("version", "VERSION")
    >>> table.attrs.VERSION
    '2.6'

Ok, that's better. If you would terminate your session now, you would be able
to use the h5ls command to read the /detector/readout attributes from the
file written to disk.

.. code-block:: bash

    $ h5ls -vr tutorial1.h5/detector/readout
    Opened "tutorial1.h5" with sec2 driver.
    /detector/readout        Dataset {10/Inf}
        Attribute: CLASS     scalar
            Type:      6-byte null-terminated ASCII string
            Data:  "TABLE"
        Attribute: VERSION   scalar
            Type:      4-byte null-terminated ASCII string
            Data:  "2.6"
        Attribute: TITLE     scalar
            Type:      16-byte null-terminated ASCII string
            Data:  "Readout example"
        Attribute: NROWS     scalar
            Type:      native long long
            Data:  10
        Attribute: FIELD_0_NAME scalar
            Type:      9-byte null-terminated ASCII string
            Data:  "ADCcount"
        Attribute: FIELD_1_NAME scalar
            Type:      9-byte null-terminated ASCII string
            Data:  "TDCcount"
        Attribute: FIELD_2_NAME scalar
            Type:      7-byte null-terminated ASCII string
            Data:  "energy"
        Attribute: FIELD_3_NAME scalar
            Type:      7-byte null-terminated ASCII string
            Data:  "grid_i"
        Attribute: FIELD_4_NAME scalar
            Type:      7-byte null-terminated ASCII string
            Data:  "grid_j"
        Attribute: FIELD_5_NAME scalar
            Type:      9-byte null-terminated ASCII string
            Data:  "idnumber"
        Attribute: FIELD_6_NAME scalar
            Type:      5-byte null-terminated ASCII string
            Data:  "name"
        Attribute: FIELD_7_NAME scalar
            Type:      9-byte null-terminated ASCII string
            Data:  "pressure"
        Attribute: FLAVOR    scalar
            Type:      5-byte null-terminated ASCII string
            Data:  "numpy"
        Attribute: tempScale scalar
            Type:      7-byte null-terminated ASCII string
            Data:  "Celsius"
        Attribute: temperature scalar
            Type:      native double
            Data:  18.4
        Location:  0:1:0:1952
        Links:     1
        Modified:  2006-12-11 10:35:13 CET
        Chunks:    {85} 3995 bytes
        Storage:   470 logical bytes, 3995 allocated bytes, 11.76% utilization
        Type:      struct {
                        "ADCcount"         +0    native unsigned short
                        "TDCcount"         +2    native unsigned char
                        "energy"           +3    native double
                        "grid_i"           +11   native int
                        "grid_j"           +15   native int
                        "idnumber"         +19   native long long
                        "name"             +27   16-byte null-terminated ASCII string
                        "pressure"         +43   native float
                    } 47 bytes

Attributes are a useful mechanism to add persistent (meta) information to
your data.


Getting object metadata
~~~~~~~~~~~~~~~~~~~~~~~
Each object in PyTables has *metadata* information about the data in the
file. Normally this *meta-information* is accessible through the node
instance variables. Let's take a look at some examples::

    >>> print("Object:", table)
    Object: /detector/readout (Table(10,)) 'Readout example'
    >>> print("Table name:", table.name)
    Table name: readout
    >>> print("Table title:", table.title)
    Table title: Readout example
    >>> print("Number of rows in table:", table.nrows)
    Number of rows in table: 10
    >>> print("Table variable names with their type and shape:")
    Table variable names with their type and shape:
    >>> for name in table.colnames:
    ...     print(name, ':= %s, %s' % (table.coldtypes[name], table.coldtypes[name].shape))
    ADCcount := uint16, ()
    TDCcount := uint8, ()
    energy := float64, ()
    grid_i := int32, ()
    grid_j := int32, ()
    idnumber := int64, ()
    name := |S16, ()
    pressure := float32, ()

Here, the name, title, nrows, colnames and coldtypes attributes (see
:class:`Table` for a complete attribute list) of the Table object gives us
quite a bit of information about the table data.

You can interactively retrieve general information about the public objects
in PyTables by asking for help::

    >>> help(table)
    Help on Table in module tables.table:
    class Table(tableextension.Table, tables.leaf.Leaf)
    |  This class represents heterogeneous datasets in an HDF5 file.
    |
    |  Tables are leaves (see the `Leaf` class) whose data consists of a
    |  unidimensional sequence of *rows*, where each row contains one or
    |  more *fields*.  Fields have an associated unique *name* and
    |  *position*, with the first field having position 0.  All rows have
    |  the same fields, which are arranged in *columns*.
    [snip]
    |
    |  Instance variables
    |  ------------------
    |
    |  The following instance variables are provided in addition to those
    |  in `Leaf`.  Please note that there are several `col` dictionaries
    |  to ease retrieving information about a column directly by its path
    |  name, avoiding the need to walk through `Table.description` or
    |  `Table.cols`.
    |
    |  autoindex
    |      Automatically keep column indexes up to date?
    |
    |      Setting this value states whether existing indexes should be
    |      automatically updated after an append operation or recomputed
    |      after an index-invalidating operation (i.e. removal and
    |      modification of rows).  The default is true.
    [snip]
    |  rowsize
    |      The size in bytes of each row in the table.
    |
    |  Public methods -- reading
    |  -------------------------
    |
    |  * col(name)
    |  * iterrows([start][, stop][, step])
    |  * itersequence(sequence)
    * itersorted(sortby[, checkCSI][, start][, stop][, step])
    |  * read([start][, stop][, step][, field][, coords])
    |  * read_coordinates(coords[, field])
    * read_sorted(sortby[, checkCSI][, field,][, start][, stop][, step])
    |  * __getitem__(key)
    |  * __iter__()
    |
    |  Public methods -- writing
    |  -------------------------
    |
    |  * append(rows)
    |  * modify_column([start][, stop][, step][, column][, colname])
    [snip]

Try getting help with other object docs by yourself::

    >>> help(h5file)
    >>> help(table.remove_rows)

To examine metadata in the */columns/pressure* Array object::

    >>> pressureObject = h5file.get_node("/columns", "pressure")
    >>> print("Info on the object:", repr(pressureObject))
    Info on the object: /columns/pressure (Array(3,)) 'Pressure column selection'
      atom := Float64Atom(shape=(), dflt=0.0)
      maindim := 0
      flavor := 'numpy'
      byteorder := 'little'
      chunkshape := None
    >>> print("  shape: ==>", pressureObject.shape)
      shape: ==> (3,)
    >>> print("  title: ==>", pressureObject.title)
      title: ==> Pressure column selection
    >>> print("  atom: ==>", pressureObject.atom)
      atom: ==> Float64Atom(shape=(), dflt=0.0)

Observe that we have used the :meth:`File.get_node` method of the File class
to access a node in the tree, instead of the natural naming method. Both are
useful, and depending on the context you will prefer one or the other.
:meth:`File.get_node` has the advantage that it can get a node from the
pathname string (as in this example) and can also act as a filter to show
only nodes in a particular location that are instances of class *classname*.
In general, however, I consider natural naming to be more elegant and easier
to use, especially if you are using the name completion capability present in
interactive console. Try this powerful combination of natural naming and
completion capabilities present in most Python consoles, and see how pleasant
it is to browse the object tree (well, as pleasant as such an activity can
be).

If you look at the type attribute of the pressureObject object, you can
verify that it is a "*float64*" array. By looking at its shape attribute, you
can deduce that the array on disk is unidimensional and has 3 elements.
See :class:`Array` or the internal doc strings for the complete Array
attribute list.


Reading data from Array objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Once you have found the desired Array, use the read() method of the Array
object to retrieve its data::

    >>> pressureArray = pressureObject.read()
    >>> pressureArray
    array([ 25.,  36.,  49.])
    >>> print("pressureArray is an object of type:", type(pressureArray))
    pressureArray is an object of type: <type 'numpy.ndarray'>
    >>> nameArray = h5file.root.columns.name.read()
    >>> print("nameArray is an object of type:", type(nameArray))
    nameArray is an object of type: <type 'list'>
    >>>
    >>> print("Data on arrays nameArray and pressureArray:")
    Data on arrays nameArray and pressureArray:
    >>> for i in range(pressureObject.shape[0]):
    ...     print(nameArray[i], "-->", pressureArray[i])
    Particle:      5 --> 25.0
    Particle:      6 --> 36.0
    Particle:      7 --> 49.0

You can see that the :meth:`Array.read` method returns an authentic NumPy
object for the pressureObject instance by looking at the output of the type()
call. A read() of the nameArray object instance returns a native Python list
(of strings). The type of the object saved is stored as an HDF5 attribute
(named FLAVOR) for objects on disk. This attribute is then read as Array
meta-information (accessible through in the Array.attrs.FLAVOR variable),
enabling the read array to be converted into the original object. This
provides a means to save a large variety of objects as arrays with the
guarantee that you will be able to later recover them in their original form.
See :meth:`File.create_array` for a complete list of supported objects for the
Array object class.


Commiting data to tables and arrays
-----------------------------------
We have seen how to create tables and arrays and how to browse both data and
metadata in the object tree. Let's examine more closely now one of the most
powerful capabilities of PyTables, namely, how to modify already created
tables and arrays [1]_


Appending data to an existing table
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Now, let's have a look at how we can add records to an existing table on
disk. Let's use our well-known *readout* Table object and append some new
values to it::

    >>> table = h5file.root.detector.readout
    >>> particle = table.row
    >>> for i in range(10, 15):
    ...     particle['name']  = f'Particle: {i:6d}'
    ...     particle['TDCcount'] = i % 256
    ...     particle['ADCcount'] = (i * 256) % (1 << 16)
    ...     particle['grid_i'] = i
    ...     particle['grid_j'] = 10 - i
    ...     particle['pressure'] = float(i*i)
    ...     particle['energy'] = float(particle['pressure'] ** 4)
    ...     particle['idnumber'] = i * (2 ** 34)
    ...     particle.append()
    >>> table.flush()

It's the same method we used to fill a new table. PyTables knows that this
table is on disk, and when you add new records, they are appended to the end
of the table [2]_.

If you look carefully at the code you will see that we have used the
table.row attribute to create a table row and fill it with the new values.
Each time that its append() method is called, the actual row is committed to
the output buffer and the row pointer is incremented to point to the next
table record. When the buffer is full, the data is saved on disk, and the
buffer is reused again for the next cycle.

*Caveat emptor*: Do not forget to always call the flush() method after a
write operation, or else your tables will not be updated!

Let's have a look at some rows in the modified table and verify that our new
data has been appended::

    >>> for r in table.iterrows():
    ...     print("%-16s | %11.1f | %11.4g | %6d | %6d | %8d \|" % \\
    ...         (r['name'], r['pressure'], r['energy'], r['grid_i'], r['grid_j'],
    ...         r['TDCcount']))
    Particle:      0 |         0.0 |           0 |      0 |     10 |        0 |
    Particle:      1 |         1.0 |           1 |      1 |      9 |        1 |
    Particle:      2 |         4.0 |         256 |      2 |      8 |        2 |
    Particle:      3 |         9.0 |        6561 |      3 |      7 |        3 |
    Particle:      4 |        16.0 |   6.554e+04 |      4 |      6 |        4 |
    Particle:      5 |        25.0 |   3.906e+05 |      5 |      5 |        5 |
    Particle:      6 |        36.0 |    1.68e+06 |      6 |      4 |        6 |
    Particle:      7 |        49.0 |   5.765e+06 |      7 |      3 |        7 |
    Particle:      8 |        64.0 |   1.678e+07 |      8 |      2 |        8 |
    Particle:      9 |        81.0 |   4.305e+07 |      9 |      1 |        9 |
    Particle:     10 |       100.0 |       1e+08 |     10 |      0 |       10 |
    Particle:     11 |       121.0 |   2.144e+08 |     11 |     -1 |       11 |
    Particle:     12 |       144.0 |     4.3e+08 |     12 |     -2 |       12 |
    Particle:     13 |       169.0 |   8.157e+08 |     13 |     -3 |       13 |
    Particle:     14 |       196.0 |   1.476e+09 |     14 |     -4 |       14 |


Modifying data in tables
~~~~~~~~~~~~~~~~~~~~~~~~
Ok, until now, we've been only reading and writing (appending) values to our
tables. But there are times that you need to modify your data once you have
saved it on disk (this is specially true when you need to modify the real
world data to adapt your goals ;).
Let's see how we can modify the values that were saved in our existing tables.
We will start modifying single cells in the first row of the Particle table::

    >>> print("Before modif-->", table[0])
    Before modif--> (0, 0, 0.0, 0, 10, 0L, 'Particle:      0', 0.0)
    >>> table.cols.TDCcount[0] = 1
    >>> print("After modifying first row of ADCcount-->", table[0])
    After modifying first row of ADCcount--> (0, 1, 0.0, 0, 10, 0L, 'Particle:      0', 0.0)
    >>> table.cols.energy[0] = 2
    >>> print("After modifying first row of energy-->", table[0])
    After modifying first row of energy--> (0, 1, 2.0, 0, 10, 0L, 'Particle:      0', 0.0)

We can modify complete ranges of columns as well::

    >>> table.cols.TDCcount[2:5] = [2,3,4]
    >>> print("After modifying slice [2:5] of TDCcount-->", table[0:5])
    After modifying slice [2:5] of TDCcount-->
    [(0, 1, 2.0, 0, 10, 0L, 'Particle:      0', 0.0)
     (256, 1, 1.0, 1, 9, 17179869184L, 'Particle:      1', 1.0)
     (512, 2, 256.0, 2, 8, 34359738368L, 'Particle:      2', 4.0)
     (768, 3, 6561.0, 3, 7, 51539607552L, 'Particle:      3', 9.0)
     (1024, 4, 65536.0, 4, 6, 68719476736L, 'Particle:      4', 16.0)]
    >>> table.cols.energy[1:9:3] = [2,3,4]
    >>> print("After modifying slice [1:9:3] of energy-->", table[0:9])
    After modifying slice [1:9:3] of energy-->
    [(0, 1, 2.0, 0, 10, 0L, 'Particle:      0', 0.0)
     (256, 1, 2.0, 1, 9, 17179869184L, 'Particle:      1', 1.0)
     (512, 2, 256.0, 2, 8, 34359738368L, 'Particle:      2', 4.0)
     (768, 3, 6561.0, 3, 7, 51539607552L, 'Particle:      3', 9.0)
     (1024, 4, 3.0, 4, 6, 68719476736L, 'Particle:      4', 16.0)
     (1280, 5, 390625.0, 5, 5, 85899345920L, 'Particle:      5', 25.0)
     (1536, 6, 1679616.0, 6, 4, 103079215104L, 'Particle:      6', 36.0)
     (1792, 7, 4.0, 7, 3, 120259084288L, 'Particle:      7', 49.0)
     (2048, 8, 16777216.0, 8, 2, 137438953472L, 'Particle:      8', 64.0)]

Check that the values have been correctly modified!

.. hint::

    remember that column TDCcount is the second one, and that energy is the
    third. Look for more info on modifying columns in
    :meth:`Column.__setitem__`.

PyTables also lets you modify complete sets of rows at the same time. As a
demonstration of these capability, see the next example::

    >>> table.modify_rows(start=1, step=3,
    ...                 rows=[(1, 2, 3.0, 4, 5, 6L, 'Particle:   None', 8.0),
    ...                       (2, 4, 6.0, 8, 10, 12L, 'Particle: None*2', 16.0)])
    2
    >>> print("After modifying the complete third row-->", table[0:5])
    After modifying the complete third row-->
    [(0, 1, 2.0, 0, 10, 0L, 'Particle:      0', 0.0)
     (1, 2, 3.0, 4, 5, 6L, 'Particle:   None', 8.0)
     (512, 2, 256.0, 2, 8, 34359738368L, 'Particle:      2', 4.0)
     (768, 3, 6561.0, 3, 7, 51539607552L, 'Particle:      3', 9.0)
     (2, 4, 6.0, 8, 10, 12L, 'Particle: None*2', 16.0)]

As you can see, the modify_rows() call has modified the rows second and fifth,
and it returned the number of modified rows.

Apart of :meth:`Table.modify_rows`, there exists another method, called
:meth:`Table.modify_column` to modify specific columns as well.

Finally, there is another way of modifying tables that is generally more
handy than the one described above. This new way uses the :meth:`Row.update` method of
the Row instance that is attached to every table, so it
is meant to be used in table iterators. Take a look at the following example::

    >>> for row in table.where('TDCcount <= 2'):
    ...     row['energy'] = row['TDCcount']*2
    ...     row.update()
    >>> print("After modifying energy column (where TDCcount <=2)-->", table[0:4])
    After modifying energy column (where TDCcount <=2)-->
    [(0, 1, 2.0, 0, 10, 0L, 'Particle:      0', 0.0)
     (1, 2, 4.0, 4, 5, 6L, 'Particle:   None', 8.0)
     (512, 2, 4.0, 2, 8, 34359738368L, 'Particle:      2', 4.0)
     (768, 3, 6561.0, 3, 7, 51539607552L, 'Particle:      3', 9.0)]

.. note::

    The authors find this way of updating tables (i.e. using Row.update())
    to be both convenient and efficient. Please make sure to use it
    extensively.

*Caveat emptor*: Currently, :meth:`Row.update` will not work (the table will
not be updated) if the loop is broken with ``break`` statement. A possible
workaround consists in manually flushing the row internal buffer by calling
``row._flushModRows()`` just before the ``break`` statement.

Modifying data in arrays
~~~~~~~~~~~~~~~~~~~~~~~~

We are now going to see how to modify data in array objects.
The basic way to do this is through the use of :meth:`Array.__setitem__`
special method. Let's see how to modify data on the pressureObject array::

    >>> pressureObject = h5file.root.columns.pressure
    >>> print("Before modif-->", pressureObject[:])
    Before modif--> [ 25.  36.  49.]
    >>> pressureObject[0] = 2
    >>> print("First modif-->", pressureObject[:])
    First modif--> [  2.  36.  49.]
    >>> pressureObject[1:3] = [2.1, 3.5]
    >>> print("Second modif-->", pressureObject[:])
    Second modif--> [ 2.   2.1  3.5]
    >>> pressureObject[::2] = [1,2]
    >>> print("Third modif-->", pressureObject[:])
    Third modif--> [ 1.   2.1  2. ]

So, in general, you can use any combination of (multidimensional) extended
slicing.

With the sole exception that you cannot use negative values for step to refer
to indexes that you want to modify. See :meth:`Array.__getitem__` for more
examples on how to use extended slicing in PyTables objects.

Similarly, with an array of strings::

    >>> nameObject = h5file.root.columns.name
    >>> print("Before modif-->", nameObject[:])
    Before modif--> ['Particle:      5', 'Particle:      6', 'Particle:      7']
    >>> nameObject[0] = 'Particle:   None'
    >>> print("First modif-->", nameObject[:])
    First modif--> ['Particle:   None', 'Particle:      6', 'Particle:      7']
    >>> nameObject[1:3] = ['Particle:      0', 'Particle:      1']
    >>> print("Second modif-->", nameObject[:])
    Second modif--> ['Particle:   None', 'Particle:      0', 'Particle:      1']
    >>> nameObject[::2] = ['Particle:     -3', 'Particle:     -5']
    >>> print("Third modif-->", nameObject[:])
    Third modif--> ['Particle:     -3', 'Particle:      0', 'Particle:     -5']


And finally... how to delete rows from a table
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We'll finish this tutorial by deleting some rows from the table we have.
Suppose that we want to delete the 5th to 9th rows (inclusive)::

    >>> table.remove_rows(5,10)
    5

:meth:`Table.remove_rows` deletes the rows in the range (start, stop). It
returns the number of rows effectively removed.

We have reached the end of this first tutorial. Don't forget to close the
file when you finish::

    >>> h5file.close()
    >>> ^D
    $

In :ref:`Figure 2 <tutorial1-2-tableview>` you can see a graphical view of the
PyTables file with the datasets we have just created. In
:ref:`tutorial1-general` are displayed the general properties of the table
/detector/readout.

.. _tutorial1-2-tableview:

.. figure:: images/tutorial1-2-tableview.png
    :align: center

    **Figure 2. The final version of the data file for tutorial 1.**


.. _tutorial1-general:

.. figure:: images/tutorial1-general.png
    :align: center

    **Figure 3. General properties of the /detector/readout table.**


.. _secondExample:

Multidimensional table cells and automatic sanity checks
--------------------------------------------------------
Now it's time for a more real-life example (i.e. with errors in the code). We
will create two groups that branch directly from the root node, Particles and
Events. Then, we will put three tables in each group. In Particles we will
put tables based on the Particle descriptor and in Events, the tables based
the Event descriptor.

Afterwards, we will provision the tables with a number of records. Finally,
we will read the newly-created table /Events/TEvent3 and select some values
from it, using a comprehension list.

Look at the next script (you can find it in :file:`examples/tutorial2.py`).
It appears to do all of the above, but it contains some small bugs. Note that
this Particle class is not directly related to the one defined in
last tutorial; this class is simpler (note, however, the *multidimensional*
columns called pressure and temperature).

We also introduce a new manner to describe a Table as a structured NumPy
dtype (or even as a dictionary), as you can see in the Event description. See
:meth:`File.create_table` for different kinds of descriptor objects that
can be passed to this method::

    import tables as tb
    import numpy as np

    # Describe a particle record
    class Particle(tb.IsDescription):
        name        = tb.StringCol(itemsize=16)  # 16-character string
        lati        = tb.Int32Col()              # integer
        longi       = tb.Int32Col()              # integer
        pressure    = tb.Float32Col(shape=(2,3)) # array of floats (single-precision)
        temperature = tb.Float64Col(shape=(2,3)) # array of doubles (double-precision)

    # Native NumPy dtype instances are also accepted
    Event = np.dtype([
        ("name"     , "S16"),
        ("TDCcount" , np.uint8),
        ("ADCcount" , np.uint16),
        ("xcoord"   , np.float32),
        ("ycoord"   , np.float32)
        ])

    # And dictionaries too (this defines the same structure as above)
    # Event = {
    #     "name"     : tb.StringCol(itemsize=16),
    #     "TDCcount" : tb.UInt8Col(),
    #     "ADCcount" : tb.UInt16Col(),
    #     "xcoord"   : tb.Float32Col(),
    #     "ycoord"   : tb.Float32Col(),
    #     }

    # Open a file in "w"rite mode
    fileh = tb.open_file("tutorial2.h5", mode="w")

    # Get the HDF5 root group
    root = fileh.root

    # Create the groups:
    for groupname in ("Particles", "Events"):
        group = fileh.create_group(root, groupname)

    # Now, create and fill the tables in Particles group
    gparticles = root.Particles

    # Create 3 new tables
    for tablename in ("TParticle1", "TParticle2", "TParticle3"):
        # Create a table
        table = fileh.create_table("/Particles", tablename, Particle, "Particles: "+tablename)

        # Get the record object associated with the table:
        particle = table.row

        # Fill the table with 257 particles
        for i in range(257):
            # First, assign the values to the Particle record
            particle['name'] = f'Particle: {i:6d}'
            particle['lati'] = i
            particle['longi'] = 10 - i

            ########### Detectable errors start here. Play with them!
            particle['pressure'] = i * np.arange(2 * 4).reshape(2, 4)  # Incorrect
            #particle['pressure'] = i * np.arange(2 * 3).reshape(2, 3) # Correct
            ########### End of errors

            particle['temperature'] = i ** 2     # Broadcasting

            # This injects the Record values
            particle.append()

        # Flush the table buffers
        table.flush()

    # Now, go for Events:
    for tablename in ("TEvent1", "TEvent2", "TEvent3"):
        # Create a table in Events group
        table = fileh.create_table(root.Events, tablename, Event, "Events: "+tablename)

        # Get the record object associated with the table:
        event = table.row

        # Fill the table with 257 events
        for i in range(257):
            # First, assign the values to the Event record
            event['name']  = f'Event: {i:6d}'
            event['TDCcount'] = i % (1<<8)   # Correct range

            ########### Detectable errors start here. Play with them!
            event['xcoor'] = float(i ** 2)     # Wrong spelling
            #event['xcoord'] = float(i ** 2)   # Correct spelling
            event['ADCcount'] = "sss"        # Wrong type
            #event['ADCcount'] = i * 2       # Correct type
            ########### End of errors

            event['ycoord'] = float(i) ** 4

            # This injects the Record values
            event.append()

        # Flush the buffers
        table.flush()

    # Read the records from table "/Events/TEvent3" and select some
    table = root.Events.TEvent3
    e = [ p['TDCcount'] for p in table if p['ADCcount'] < 20 and 4 <= p['TDCcount'] < 15 ]
    print(f"Last record ==> {p}")
    print("Selected values ==> {e}")
    print("Total selected records ==> {len(e)}")

    # Finally, close the file (this also will flush all the remaining buffers!)
    fileh.close()


Shape checking
~~~~~~~~~~~~~~
If you look at the code carefully, you'll see that it won't work. You will
get the following error.

.. code-block:: bash

    $ python3 tutorial2.py
    Traceback (most recent call last):
      File "tutorial2.py", line 60, in <module>
        particle['pressure'] = array(i * arange(2 * 3)).reshape((2, 4))  # Incorrect
    ValueError: total size of new array must be unchanged
    Closing remaining open files: tutorial2.h5... done

This error indicates that you are trying to assign an array with an
incompatible shape to a table cell. Looking at the source, we see that we
were trying to assign an array of shape (2,4) to a pressure element, which
was defined with the shape (2,3).

In general, these kinds of operations are forbidden, with one valid
exception: when you assign a *scalar* value to a multidimensional column
cell, all the cell elements are populated with the value of the scalar.
For example::

    particle['temperature'] = i ** 2    # Broadcasting

The value i**2 is assigned to all the elements of the temperature table cell.
This capability is provided by the NumPy package and is known as
*broadcasting*.


Field name checking
~~~~~~~~~~~~~~~~~~~
After fixing the previous error and rerunning the program, we encounter
another error.

.. code-block:: bash

    $ python3 tutorial2.py
    Traceback (most recent call last):
      File "tutorial2.py", line 73, in ?
        event['xcoor'] = float(i ** 2)     # Wrong spelling
      File "tableextension.pyx", line 1094, in tableextension.Row.__setitem__
      File "tableextension.pyx", line 127, in tableextension.get_nested_field_cache
      File "utilsextension.pyx", line 331, in utilsextension.get_nested_field
    KeyError: 'no such column: xcoor'

This error indicates that we are attempting to assign a value to a
non-existent field in the *event* table object. By looking carefully at the
Event class attributes, we see that we misspelled the xcoord field (we wrote
xcoor instead). This is unusual behavior for Python, as normally when you
assign a value to a non-existent instance variable, Python creates a new
variable with that name. Such a feature can be dangerous when dealing with an
object that contains a fixed list of field names. PyTables checks that the
field exists and raises a KeyError if the check fails.


Data type checking
~~~~~~~~~~~~~~~~~~
Finally, the last issue which we will find here is a TypeError exception.

.. code-block:: bash

    $ python3 tutorial2.py
    Traceback (most recent call last):
      File "tutorial2.py", line 75, in ?
        event['ADCcount'] = "sss"          # Wrong type
      File "tableextension.pyx", line 1111, in tableextension.Row.__setitem__
    TypeError: invalid type (<type 'str'>) for column ``ADCcount``

And, if we change the affected line to read::

    event.ADCcount = i * 2        # Correct type

we will see that the script ends well.

You can see the structure created with this (corrected) script in
:ref:`Figure 4 <tutorial2-tableview>`.
In particular, note the multidimensional column cells in table
/Particles/TParticle2.

.. _tutorial2-tableview:

.. figure:: images/tutorial2-tableview.png
    :align: center

    **Figure 4. Table hierarchy for tutorial 2.**


.. _LinksTutorial:

Using links for more convenient access to nodes
-----------------------------------------------
Links are special nodes that can be used to create additional paths to your
existing nodes.  PyTables supports three kinds of links: hard links, soft
links (aka symbolic links) and external links.

Hard links let the user create additional paths to access another node in the
same file, and once created, they are indistinguishable from the referred
node object, except that they have different paths in the object tree.  For
example, if the referred node is, say, a Table object, then the new hard link
will become a Table object itself.  From this point on, you will be able to
access the same Table object from two different paths: the original one and
the new hard link path.  If you delete one path to the table, you will be
able to reach it via the other path.

Soft links are similar to hard links, but they keep their own personality.
When you create a soft link to another node, you will get a new SoftLink
object that *refers* to that node.  However, in order to access the referred
node, you need to *dereference* it.

Finally, external links are like soft links, with the difference that these
are meant to point to nodes in *external* files instead of nodes in the same
file.  They are represented by the ExternalLink class and, like soft links,
you need to dereference them in order to get access to the pointed node.


Interactive example
~~~~~~~~~~~~~~~~~~~
Now we are going to learn how to deal with links. You can find the code used
in this section in :file:`examples/links.py`.

First, let's create a file with some group structure::

    >>> import tables as tb
    >>> f1 = tb.open_file('links1.h5', 'w')
    >>> g1 = f1.create_group('/', 'g1')
    >>> g2 = f1.create_group(g1, 'g2')

Now, we will put some datasets on the /g1 and /g1/g2 groups::

    >>> a1 = f1.create_carray(g1, 'a1', tb.Int64Atom(), shape=(10000,))
    >>> t1 = f1.create_table(g2, 't1', {'f1': tb.IntCol(), 'f2': tb.FloatCol()})

We can start the party now.  We are going to create a new group, say /gl,
where we will put our links and will start creating one hard link too::

    >>> gl = f1.create_group('/', 'gl')
    >>> ht = f1.create_hard_link(gl, 'ht', '/g1/g2/t1')  # ht points to t1
    >>> print(f"``{ht}`` is a hard link to: ``{t1}``")
    ``/gl/ht (Table(0,)) `` is a hard link to: ``/g1/g2/t1 (Table(0,)) ``

You can see how we've created a hard link in /gl/ht which is pointing to the
existing table in /g1/g2/t1.  Have look at how the hard link is represented;
it looks like a table, and actually, it is a *real* table.  We have two
different paths to access that table, the original /g1/g2/t1 and the new one
/gl/ht.  If we remove the original path we still can reach the table by using
the new path::

    >>> t1.remove()
    >>> print(f"table continues to be accessible in: ``{f1.get_node('/gl/ht')}``")
    table continues to be accessible in: ``/gl/ht (Table(0,)) ``

So far so good. Now, let's create a couple of soft links::

    >>> la1 = f1.create_soft_link(gl, 'la1', '/g1/a1')  # la1 points to a1
    >>> print(f"``{la1}`` is a soft link to: ``{la1.target}``")
    ``/gl/la1 (SoftLink) -> /g1/a1`` is a soft link to: ``/g1/a1``
    >>> lt = f1.create_soft_link(gl, 'lt', '/g1/g2/t1')  # lt points to t1
    >>> print(f"``{lt}`` is a soft link to: ``{lt.target}``")
    ``/gl/lt (SoftLink) -> /g1/g2/t1 (dangling)`` is a soft link to: ``/g1/g2/t1``

Okay, we see how the first link /gl/la1 points to the array /g1/a1.  Notice
how the link prints as a SoftLink, and how the referred node is stored in the
target instance attribute.  The second link (/gt/lt) pointing to /g1/g2/t1
also has been created successfully, but by better inspecting the string
representation of it, we see that is labeled as '(dangling)'.  Why is this?
Well, you should remember that we recently removed the /g1/g2/t1 path to
access table t1.  When printing it, the object knows that it points to
*nowhere* and reports this.
This is a nice way to quickly know whether a soft link points to an exiting
node or not.

So, let's re-create the removed path to t1 table::

    >>> t1 = f1.create_hard_link('/g1/g2', 't1', '/gl/ht')
    >>> print(f"``{lt}`` is not dangling anymore")
    ``/gl/lt (SoftLink) -> /g1/g2/t1`` is not dangling anymore

and the soft link is pointing to an existing node now.

Of course, for soft links to serve any actual purpose we need a way to get
the pointed node.  It happens that soft links are callable, and that's the
way to get the referred nodes back::

    >>> plt = lt()
    >>> print(f"dereferred lt node: ``{plt}``")
    dereferred lt node: ``/g1/g2/t1 (Table(0,)) ``
    >>> pla1 = la1()
    >>> print(f"dereferred la1 node: ``{pla1}``")
    dereferred la1 node: ``/g1/a1 (CArray(10000,)) ``

Now, plt is a Python reference to the t1 table while pla1 refers to the a1
array.  Easy, uh?

Let's suppose now that a1 is an array whose access speed is critical for our
application.  One possible solution is to move the entire file into a faster
disk, say, a solid state disk so that access latencies can be reduced quite a
lot.  However, it happens that our file is too big to fit into our shiny new
(although small in capacity) SSD disk.  A solution is to copy just the a1
array into a separate file that would fit into our SSD disk.  However, our
application would be able to handle two files instead of only one, adding
significantly more complexity, which is not a good thing.

External links to the rescue!  As we've already said, external links are like
soft links, but they are designed to link objects in external files.  Back to
our problem, let's copy the a1 array into a different file::

    >>> f2 = tb.open_file('links2.h5', 'w')
    >>> new_a1 = a1.copy(f2.root, 'a1')
    >>> f2.close()  # close the other file

And now, we can remove the existing soft link and create the external link in
its place::

    >>> la1.remove()
    >>> la1 = f1.create_external_link(gl, 'la1', 'links2.h5:/a1')
    >>> print(f"``{la1}`` is an external link to: ``{la1.target}``")
    ``/gl/la1 (ExternalLink) -> links2.h5:/a1`` is an external link to: ``links2.h5:/a1``

Let's try dereferring it::

    >>> new_a1 = la1()  # dereferrencing la1 returns a1 in links2.h5
    >>> print(f"dereferred la1 node:  ``{new_a1}``")
    dereferred la1 node:  ``/a1 (CArray(10000,)) ``

Well, it seems like we can access the external node.  But just to make sure
that the node is in the other file::

    >>> print("new_a1 file:", new_a1._v_file.filename)
    new_a1 file: links2.h5

Okay, the node is definitely in the external file.  So, you won't have to
worry about your application: it will work exactly the same no matter the
link is internal (soft) or external.

Finally, here it is a dump of the objects in the final file, just to get a
better idea of what we ended with::

    >>> f1.close()
    >>> exit()
    $ ptdump links1.h5
    / (RootGroup) ''
    /g1 (Group) ''
    /g1/a1 (CArray(10000,)) ''
    /gl (Group) ''
    /gl/ht (Table(0,)) ''
    /gl/la1 (ExternalLink) -> links2.h5:/a1
    /gl/lt (SoftLink) -> /g1/g2/t1
    /g1/g2 (Group) ''
    /g1/g2/t1 (Table(0,)) ''

This ends this tutorial.  I hope it helped you to appreciate how useful links
can be.  I'm sure you will find other ways in which you can use links that
better fit your own needs.


Exercising the Undo/Redo feature
--------------------------------
PyTables has integrated support for undoing and/or redoing actions. This
functionality lets you put marks in specific places of your hierarchy
manipulation operations, so that you can make your HDF5 file pop back
(*undo*) to a specific mark (for example for inspecting how your hierarchy
looked at that point). You can also go forward to a more recent marker
(*redo*). You can even do jumps to the marker you want using just one
instruction as we will see shortly.

You can undo/redo all the operations that are related to object tree
management, like creating, deleting, moving or renaming nodes (or complete
sub-hierarchies) inside a given object tree. You can also undo/redo
operations (i.e. creation, deletion or modification) of persistent node
attributes. However, actions which include *internal* modifications of datasets
(that includes Table.append, Table.modify_rows or Table.remove_rows, among others) 
cannot be currently undone/redone.

This capability can be useful in many situations, for example when doing
simulations with multiple branches. When you have to choose a path to follow
in such a situation, you can put a mark there and, if the simulation is not
going well, you can go back to that mark and start another path. Another
possible application is defining coarse-grained operations which operate in a
transactional-like way, i.e. which return the database to its previous state
if the operation finds some kind of problem while running. You can probably
devise many other scenarios where the Undo/Redo feature can be useful to you
[3]_.


A basic example
~~~~~~~~~~~~~~~
In this section, we are going to show the basic behavior of the Undo/Redo
feature. You can find the code used in this example in
:file:`examples/tutorial3-1.py`. A somewhat more complex example will be
explained in the next section.

First, let's create a file::

    >>> import tables
    >>> fileh = tables.open_file("tutorial3-1.h5", "w", title="Undo/Redo demo 1")

And now, activate the Undo/Redo feature with the method
:meth:`File.enable_undo` of File::

    >>> fileh.enable_undo()

From now on, all our actions will be logged internally by PyTables. Now, we
are going to create a node (in this case an Array object)::

    >>> one = fileh.create_array('/', 'anarray', [3,4], "An array")

Now, mark this point::

    >>> fileh.mark()
    1

We have marked the current point in the sequence of actions.
In addition, the mark() method has returned the identifier assigned to this
new mark, that is 1 (mark #0 is reserved for the implicit mark at the
beginning of the action log). In the next section we will see that you can
also assign a *name* to a mark (see :meth:`File.mark` for more info on
mark()).
Now, we are going to create another array::

    >>> another = fileh.create_array('/', 'anotherarray', [4,5], "Another array")

Right. Now, we can start doing funny things. Let's say that we want to pop
back to the previous mark (that whose value was 1, do you remember?). Let's
introduce the undo() method (see :meth:`File.undo`)::

    >>> fileh.undo()

Fine, what do you think happened? Well, let's have a look at the object
tree::

    >>> print(fileh)
    tutorial3-1.h5 (File) 'Undo/Redo demo 1'
    Last modif.: 'Tue Mar 13 11:43:55 2007'
    Object Tree:
    / (RootGroup) 'Undo/Redo demo 1'
    /anarray (Array(2,)) 'An array'

What happened with the /anotherarray node we've just created? You've guessed it,
it has disappeared because it was created *after* the mark 1. If you are
curious enough you may ask - where has it gone? Well, it has not been
deleted completely; it has just been moved into a special, hidden group of
PyTables that renders it invisible and waiting for a chance to be reborn.

Now, unwind once more, and look at the object tree::

    >>> fileh.undo()
    >>> print(fileh)
    tutorial3-1.h5 (File) 'Undo/Redo demo 1'
    Last modif.: 'Tue Mar 13 11:43:55 2007'
    Object Tree:
    / (RootGroup) 'Undo/Redo demo 1'

Oops, /anarray has disappeared as well!
Don't worry, it will visit us again very shortly. So, you might be somewhat lost
right now; at which mark are we? Let's ask the :meth:`File.get_current_mark`
method in the file handler::

    >>> print(fileh.get_current_mark())
    0

So we are at mark #0, remember? Mark #0 is an implicit mark that is created
when you start the log of actions when calling File.enable_undo(). Fine, but
you are missing your too-young-to-die arrays. What can we do about that?
:meth:`File.redo` to the rescue::

    >>> fileh.redo()
    >>> print(fileh)
    tutorial3-1.h5 (File) 'Undo/Redo demo 1'
    Last modif.: 'Tue Mar 13 11:43:55 2007'
    Object Tree:
    / (RootGroup) 'Undo/Redo demo 1'
    /anarray (Array(2,)) 'An array'

Great! The /anarray array has come into life again. Just check that it is
alive and well::

    >>> fileh.root.anarray.read()
    [3, 4]
    >>> fileh.root.anarray.title
    'An array'

Well, it looks pretty similar to its past life;
what's more, it is exactly the same object::

    >>> fileh.root.anarray is one
    True

It was just moved to the hidden group and back again, that's all!
That's kind of fun, so we are going to do the same with /anotherarray::

    >>> fileh.redo()
    >>> print(fileh)
    tutorial3-1.h5 (File) 'Undo/Redo demo 1'
    Last modif.: 'Tue Mar 13 11:43:55 2007'
    Object Tree:
    / (RootGroup) 'Undo/Redo demo 1'
    /anarray (Array(2,)) 'An array'
    /anotherarray (Array(2,)) 'Another array'

Welcome back, /anotherarray! Just a couple of sanity checks::

    >>> assert fileh.root.anotherarray.read() == [4,5]
    >>> assert fileh.root.anotherarray.title == "Another array"
    >>> fileh.root.anotherarray is another
    True

Nice, you managed to turn your data back to life.
Congratulations! But wait, do not forget to close your action log when you
don't need this feature anymore::

    >>> fileh.disable_undo()

That will allow you to continue working with your data without actually
requiring PyTables to keep track of all your actions, and more importantly,
allowing your objects to die completely if they have to, not requiring to
keep them anywhere, hence saving process time and space in your database
file.


A more complete example
~~~~~~~~~~~~~~~~~~~~~~~
Now, time for a somewhat more sophisticated demonstration of the Undo/Redo
feature. In it, several marks will be set in different parts of the code flow
and we will see how to jump between these marks with just one method call.
You can find the code used in this example in :file:`examples/tutorial3-2.py`

Let's introduce the first part of the code::

    import tables

    # Create an HDF5 file
    fileh = tables.open_file('tutorial3-2.h5', 'w', title='Undo/Redo demo 2')

            #'-**-**-**-**-**-**- enable undo/redo log  -**-**-**-**-**-**-**-'
    fileh.enable_undo()

    # Start undoable operations
    fileh.create_array('/', 'otherarray1', [3,4], 'Another array 1')
    fileh.create_group('/', 'agroup', 'Group 1')

    # Create a 'first' mark
    fileh.mark('first')
    fileh.create_array('/agroup', 'otherarray2', [4,5], 'Another array 2')
    fileh.create_group('/agroup', 'agroup2', 'Group 2')

    # Create a 'second' mark
    fileh.mark('second')
    fileh.create_array('/agroup/agroup2', 'otherarray3', [5,6], 'Another array 3')

    # Create a 'third' mark
    fileh.mark('third')
    fileh.create_array('/', 'otherarray4', [6,7], 'Another array 4')
    fileh.create_array('/agroup', 'otherarray5', [7,8], 'Another array 5')

You can see how we have set several marks interspersed in the code flow,
representing different states of the database. Also, note that we have
assigned *names* to these marks, namely 'first', 'second' and 'third'.

Now, start doing some jumps back and forth in the states of the database::

    # Now go to mark 'first'
    fileh.goto('first')
    assert '/otherarray1' in fileh
    assert '/agroup' in fileh
    assert '/agroup/agroup2' not in fileh
    assert '/agroup/otherarray2' not in fileh
    assert '/agroup/agroup2/otherarray3' not in fileh
    assert '/otherarray4' not in fileh
    assert '/agroup/otherarray5' not in fileh

    # Go to mark 'third'
    fileh.goto('third')
    assert '/otherarray1' in fileh
    assert '/agroup' in fileh
    assert '/agroup/agroup2' in fileh
    assert '/agroup/otherarray2' in fileh
    assert '/agroup/agroup2/otherarray3' in fileh
    assert '/otherarray4' not in fileh
    assert '/agroup/otherarray5' not in fileh

    # Now go to mark 'second'
    fileh.goto('second')
    assert '/otherarray1' in fileh
    assert '/agroup' in fileh
    assert '/agroup/agroup2' in fileh
    assert '/agroup/otherarray2' in fileh
    assert '/agroup/agroup2/otherarray3' not in fileh
    assert '/otherarray4' not in fileh
    assert '/agroup/otherarray5' not in fileh

Well, the code above shows how easy is to jump to a certain mark in the
database by using the :meth:`File.goto` method.

There are also a couple of implicit marks for going to the beginning or the
end of the saved states: 0 and -1. Going to mark #0 means go to the beginning
of the saved actions, that is, when method fileh.enable_undo() was called.
Going to mark #-1 means go to the last recorded action, that is the last
action in the code flow.

Let's see what happens when going to the end of the action log::

    # Go to the end
    fileh.goto(-1)
    assert '/otherarray1' in fileh
    assert '/agroup' in fileh
    assert '/agroup/agroup2' in fileh
    assert '/agroup/otherarray2' in fileh
    assert '/agroup/agroup2/otherarray3' in fileh
    assert '/otherarray4' in fileh
    assert '/agroup/otherarray5' in fileh

    # Check that objects have come back to life in a sane state
    assert fileh.root.otherarray1.read() == [3,4]
    assert fileh.root.agroup.otherarray2.read() == [4,5]
    assert fileh.root.agroup.agroup2.otherarray3.read() == [5,6]
    assert fileh.root.otherarray4.read() == [6,7]
    assert fileh.root.agroup.otherarray5.read() == [7,8]

Try going to the beginning of the action log yourself (remember, mark #0)
and check contents of the object tree.

We have nearly finished this demonstration. As always, do not forget to close
the action log as well as the database::

    #'-**-**-**-**-**-**- disable undo/redo log  -**-**-**-**-**-**-**-'
    fileh.disable_undo()
    # Close the file
    fileh.close()

You might want to check other examples on Undo/Redo feature that appear in
:file:`examples/undo-redo.py`.


Using enumerated types
----------------------
PyTables includes support for handling enumerated types. Those types are
defined by providing an exhaustive *set* or *list* of possible, named values
for a variable of that type. Enumerated variables of the same type are
usually compared between them for equality and sometimes for order, but are
not usually operated upon.

Enumerated values have an associated *name* and *concrete value*. Every name
is unique and so are concrete values. An enumerated variable always takes the
concrete value, not its name. Usually, the concrete value is not used
directly, and frequently it is entirely irrelevant. For the same reason, an
enumerated variable is not usually compared with concrete values out of its
enumerated type. For that kind of use, standard variables and constants are
more adequate.

PyTables provides the Enum (see :ref:`EnumClassDescr`) class to provide
support for enumerated types. Each instance of Enum is an enumerated type (or
*enumeration*). For example, let us create an enumeration of colors

All these examples can be found in :file:`examples/play-with-enums.py`::

    >>> import tables
    >>> colorList = ['red', 'green', 'blue', 'white', 'black']
    >>> colors = tables.Enum(colorList)

Here we used a simple list giving the names of enumerated values, but we left
the choice of concrete values up to the Enum class. Let us see the enumerated
pairs to check those values::

    >>> print("Colors:", [v for v in colors])
    Colors: [('blue', 2), ('black', 4), ('white', 3), ('green', 1), ('red', 0)]

Names have been given automatic integer concrete values. We can iterate over
values in an enumeration, but we will usually be more interested in
accessing single values. We can get the concrete value associated with a name
by accessing it as an attribute or as an item (the later can be useful for
names not resembling Python identifiers)::

    >>> print("Value of 'red' and 'white':", (colors.red, colors.white))
    Value of 'red' and 'white': (0, 3)
    >>> print("Value of 'yellow':", colors.yellow)
    Value of 'yellow':
    Traceback (most recent call last):
      File "<stdin>", line 1, in ?
      File ".../tables/misc/enum.py", line 230, in __getattr__
        raise AttributeError(\*ke.args)
    AttributeError: no enumerated value with that name: 'yellow'
    >>>
    >>> print("Value of 'red' and 'white':", (colors['red'], colors['white']))
    Value of 'red' and 'white': (0, 3)
    >>> print("Value of 'yellow':", colors['yellow'])
    Value of 'yellow':
    Traceback (most recent call last):
      File "<stdin>", line 1, in ?
      File ".../tables/misc/enum.py", line 189, in __getitem__
        raise KeyError("no enumerated value with that name: %r" % (name,))
    KeyError: "no enumerated value with that name: 'yellow'"

See how accessing a value that is not in the enumeration raises the
appropriate exception. We can also do the opposite and get the name
that matches a concrete value by using the __call__() method of Enum::

    >>> print(f"Name of value {colors.red}:", colors(colors.red))
    Name of value 0: red
    >>> print("Name of value 1234:", colors(1234))
    Name of value 1234:
    Traceback (most recent call last):
      File "<stdin>", line 1, in ?
      File ".../tables/misc/enum.py", line 320, in __call__
        raise ValueError(
    ValueError: no enumerated value with that concrete value: 1234

You can see what we made as using the enumerated type to *convert* a concrete
value into a name in the enumeration. Of course, values out of the
enumeration can not be converted.


Enumerated columns
~~~~~~~~~~~~~~~~~~
Columns of an enumerated type can be declared by using the EnumCol (see
:ref:`ColClassDescr`) class. To see how this works, let us open a new
PyTables file and create a table to collect the simulated results of a
probabilistic experiment. In it, we have a bag full of colored balls; we take
a ball out and annotate the time of extraction and the color of the ball::

    >>> h5f = tables.open_file('enum.h5', 'w')
    >>> class BallExt(tables.IsDescription):
    ...     ballTime = tables.Time32Col()
    ...     ballColor = tables.EnumCol(colors, 'black', base='uint8')
    >>> tbl = h5f.create_table('/', 'extractions', BallExt, title="Random ball extractions")
    >>>

We declared the ballColor column to be of the enumerated type colors, with a
default value of black. We also stated that we are going to store concrete
values as unsigned 8-bit integer values [4]_.

Let us use some random values to fill the table::

    >>> import time
    >>> import random
    >>> now = time.time()
    >>> row = tbl.row
    >>> for i in range(10):
    ...     row['ballTime'] = now + i
    ...     row['ballColor'] = colors[random.choice(colorList)]  # take note of this
    ...     row.append()
    >>>

Notice how we used the __getitem__() call of colors to get the concrete value
to store in ballColor. This way of appending values to a table automatically checks
for the validity on enumerated values. For instance::

    >>> row['ballTime'] = now + 42
    >>> row['ballColor'] = 1234
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "tableextension.pyx", line 1086, in tableextension.Row.__setitem__
      File ".../tables/misc/enum.py", line 320, in __call__
        "no enumerated value with that concrete value: %r" % (value,))
    ValueError: no enumerated value with that concrete value: 1234

Note that this check is performed *only* by row.append() and not in other
methods such as tbl.append() or tbl.modify_rows(). Now, after flushing the
table we can see the result of insertions::

    >>> tbl.flush()
    >>> for r in tbl:
    ...     ballTime = r['ballTime']
    ...     ballColor = colors(r['ballColor'])  # notice this
    ...     print("Ball extracted on %d is of color %s." % (ballTime, ballColor))
    Ball extracted on 1173785568 is of color green.
    Ball extracted on 1173785569 is of color black.
    Ball extracted on 1173785570 is of color white.
    Ball extracted on 1173785571 is of color black.
    Ball extracted on 1173785572 is of color black.
    Ball extracted on 1173785573 is of color red.
    Ball extracted on 1173785574 is of color green.
    Ball extracted on 1173785575 is of color red.
    Ball extracted on 1173785576 is of color white.
    Ball extracted on 1173785577 is of color white.

As a final note, you may be wondering how to access an enumeration associated 
with ballColor once the file is closed and reopened. You can call tbl.get_enum('ballColor')
(see :meth:`Table.get_enum`) to get the enumeration back.


Enumerated arrays
~~~~~~~~~~~~~~~~~
EArray and VLArray leaves can also be declared to store enumerated values by
means of the EnumAtom (see :ref:`AtomClassDescr`) class, which works very
much like EnumCol for tables. Also, Array leaves can be used to open native
HDF enumerated arrays.

Let us create a sample EArray containing ranges of working days as
bidimensional values::

    >>> workingDays = {'Mon': 1, 'Tue': 2, 'Wed': 3, 'Thu': 4, 'Fri': 5}
    >>> dayRange = tables.EnumAtom(workingDays, 'Mon', base='uint16')
    >>> earr = h5f.create_earray('/', 'days', dayRange, (0, 2), title="Working day ranges")
    >>> earr.flavor = 'python'

Nothing unusual, except for two details. Firstly, we use
a *dictionary* instead of a list to explicitly set concrete values in the
enumeration. Secondly, there is no explicit Enum instance created!
Instead, the dictionary is passed as the first argument to the constructor of
EnumAtom. If the constructor receives a list or a dictionary instead of an
enumeration, it automatically builds the enumeration from it.

Now let us feed some data to the array::

    >>> wdays = earr.get_enum()
    >>> earr.append([(wdays.Mon, wdays.Fri), (wdays.Wed, wdays.Fri)])
    >>> earr.append([(wdays.Mon, 1234)])

Please note that, since we had no explicit Enum instance, we were forced to
use get_enum() (see :ref:`EArrayMethodsDescr`) to get it from the array (we
could also have used dayRange.enum).  Also note that we were able to append
an invalid value (1234). Array methods do not check the validity of
enumerated values.

Finally, we will print contents of the array::

    >>> for (d1, d2) in earr:
    ...     print(f"From {wdays(d1) to {wdays(d2) ({d2 - d1 + 1} days).")
    From Mon to Fri (5 days).
    From Wed to Fri (3 days).
    Traceback (most recent call last):
      File "<stdin>", line 2, in <module>
      File ".../tables/misc/enum.py", line 320, in __call__
        "no enumerated value with that concrete value: %r" % (value,))
    ValueError: no enumerated value with that concrete value: 1234

That was an example of operating on concrete values. It also showed how the
value-to-name conversion failed because of the value not belonging to the
enumeration.

Now we will close the file, and this little tutorial on enumerated types is
done::

    >>> h5f.close()


Nested structures in tables
----------------------------------------
PyTables supports handling of nested structures (or, in other words, nested datatypes)
in table objects, allowing you to define nested columns of arbitrary depth. 

Why is that useful? Suppose your data has a certain structure on the column level, 
which you would like to represent in the data model. You can do so by creating 
nested subclasses of IsDescription. The benefit is the ability to group
and retrieve data more easily. Example below may be a bit silly, but it will serve 
as an illustration of the concept::

    import tables as tb

    class Info(tb.IsDescription):
        """A sub-structure of NestedDescr"""
        _v_pos = 2   # The position in the whole structure
        name = tb.StringCol(10)
        value = tb.Float64Col(pos=0)

    colors = tb.Enum(['red', 'green', 'blue'])

    class NestedDescr(tb.IsDescription):
        """A description that has several nested columns"""
        color = tb.EnumCol(colors, 'red', base='uint32')
        info1 = tb.Info()

        class info2(tb.IsDescription):
            _v_pos = 1
            name = tb.StringCol(10)
            value = tb.Float64Col(pos=0)

            class info3(tb.IsDescription):
                x = tb.Float64Col(dflt=1)
                y = tb.UInt8Col(dflt=1)

NestedDescr is the root class with two *substructures* in it: info1 and info2. Note info1 is an instance of class Info which is defined prior to NestedDescr. info2 is declared within NestedDescr. Also, there is a third substructure, info3, that is in turn declared within substructure info2. You can define positions of substructures in the containing object by declaring special class attribute _v_pos.


Creating nested tables
~~~~~~~~~~~~~~~~~~~~~~
Now that we have defined our nested structure, let's create a *nested* table,
that is a table with columns which contain subcolumns::

    >>> fileh = tb.open_file("nested-tut.h5", "w")
    >>> table = fileh.create_table(fileh.root, 'table', NestedDescr)

Done! Now, to populate the table with values, assign a value to each field. Referencing nested fields can be accomplished by providing a full path. Follow the structure defined earlier - use '/' to access each sublevel of the structure, similar to how you would access a subdirectory on a Unix filesystem::

    >>> row = table.row
    >>> for i in range(10):
    ...     row['color'] = colors[['red', 'green', 'blue'][i % 3]]
    ...     row['info1/name'] = f"name1-{i}"
    ...     row['info2/name'] = f"name2-{i}"
    ...     row['info2/info3/y'] =  i
    ...     # Remaining fields will be filled with defaults
    ...     row.append()
    >>> table.flush()
    >>> table.nrows
    10

As demonstrated above, substructure's field can be accessed by specifying its full path as defined in the table hierarchy.


Reading nested tables
~~~~~~~~~~~~~~~~~~~~~
Now, what happens if we want to read the table? What kind of data container
would we get? Let's find out::

    >>> nra = table[::4]
    >>> nra
    array([(((1.0, 0), 'name2-0', 0.0), ('name1-0', 0.0), 0L),
           (((1.0, 4), 'name2-4', 0.0), ('name1-4', 0.0), 1L),
           (((1.0, 8), 'name2-8', 0.0), ('name1-8', 0.0), 2L)],
           dtype=[('info2', [('info3', [('x', '>f8'), ('y', '\|u1')]),
                  ('name', '\|S10'), ('value', '>f8')]),
                  ('info1', [('name', '\|S10'), ('value', '>f8')]),
                  ('color', '>u4')])

What we've got is a NumPy array with a *compound, nested datatype*, i.e. its dtype is
a list of name-datatype tuples. For every fourth row in the table (note [::4]), we get one resulting row,
giving a total of three.

You can make use of the above object in many different ways.
For example, you can use it to append new data to an existing table object::

    >>> table.append(nra)
    >>> table.nrows
    13

Or to create new tables::

    >>> table2 = fileh.create_table(fileh.root, 'table2', nra)
    >>> table2[:]
    array([(((1.0, 0), 'name2-0', 0.0), ('name1-0', 0.0), 0L),
           (((1.0, 4), 'name2-4', 0.0), ('name1-4', 0.0), 1L),
           (((1.0, 8), 'name2-8', 0.0), ('name1-8', 0.0), 2L)],
           dtype=[('info2', [('info3', [('x', '<f8'), ('y', '\|u1')]),
                  ('name', '\|S10'), ('value', '<f8')]),
                  ('info1', [('name', '\|S10'), ('value', '<f8')]),
                  ('color', '<u4')])

Finally, you can select nested values that meet an arbitrary condition::

    >>> names = [ x['info2/name'] for x in table if x['color'] == colors.red ]
    >>> names
    ['name2-0', 'name2-3', 'name2-6', 'name2-9', 'name2-0']

Note that row accessor does not provide natural naming feature, so
you have to specify an absolute path to your desired column in order to
access it.


Using Cols accessor
~~~~~~~~~~~~~~~~~~~
We can use cols attribute object (see :ref:`ColsClassDescr`) of the table
to conveniently access data stored in a substructure::

    >>> table.cols.info2[1:5]
    array([((1.0, 1), 'name2-1', 0.0), ((1.0, 2), 'name2-2', 0.0),
           ((1.0, 3), 'name2-3', 0.0), ((1.0, 4), 'name2-4', 0.0)],
           dtype=[('info3', [('x', '<f8'), ('y', '\|u1')]), ('name', '\|S10'),
                  ('value', '<f8')])

Here, we have made use of the cols accessor to access *info2*
substructure and a slice operation to retrieve a subset of data we
were interested in; you probably recognized natural naming approach
used here as well. We can continue and ask for data in *info3* substructure::

    >>> table.cols.info2.info3[1:5]
    array([(1.0, 1), (1.0, 2), (1.0, 3), (1.0, 4)],
           dtype=[('x', '<f8'), ('y', '\|u1')])

You can also use the _f_col method to get a handler for the column::

    >>> table.cols._f_col('info2')
    /table.cols.info2 (Cols), 3 columns
      info3 (Cols(), Description)
      name (Column(), \|S10)
      value (Column(), float64)

Here, you've got another Cols object handler because *info2* was a nested
column. If you select a non-nested column, you will get a regular Column
instance::

    >>> table.cols._f_col('info2/info3/y')
    /table.cols.info2.info3.y (Column(), uint8, idx=None)

To summarize, cols accessor is a very handy and powerful tool to access data
in nested tables. Don't hesitate to use it, especially when doing interactive work.


Accessing meta-information of nested tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tables have a *description* attribute, which returns an instance of
the Description class (see :ref:`DescriptionClassDescr`) with table meta-data.
It can be helpful in understanding the table structure, including nested columns::

    >>> table.description
    {
      "info2": {
        "info3": {
          "x": Float64Col(shape=(), dflt=1.0, pos=0),
          "y": UInt8Col(shape=(), dflt=1, pos=1)},
        "name": StringCol(itemsize=10, shape=(), dflt='', pos=1),
        "value": Float64Col(shape=(), dflt=0.0, pos=2)},
      "info1": {
        "name": StringCol(itemsize=10, shape=(), dflt='', pos=0),
        "value": Float64Col(shape=(), dflt=0.0, pos=1)},
      "color": EnumCol(enum=Enum({'blue': 2, 'green': 1, 'red': 0}), dflt='red',
                       base=UInt32Atom(shape=(), dflt=0), shape=(), pos=2)}

As you can see, it provides very useful information on both the format and
the structure of columns in your table.

You can use natural naming approach with the description attribute to gain access to subcolumn's
meta-data::

    >>> table.description.info1
    {"name": StringCol(itemsize=10, shape=(), dflt='', pos=0),
     "value": Float64Col(shape=(), dflt=0.0, pos=1)}
    >>> table.description.info2.info3
    {"x": Float64Col(shape=(), dflt=1.0, pos=0),
     "y": UInt8Col(shape=(), dflt=1, pos=1)}

_v_nested_names attribute provides names of columns as well as structures embedded in them::

    >>> table.description._v_nested_names
    [('info2', [('info3', ['x', 'y']), 'name', 'value']),
     ('info1', ['name', 'value']), 'color']
    >>> table.description.info1._v_nested_names
    ['name', 'value']

Regardless of which level of the structure is accessed, the output of _v_nested_names contains
the same kind of information. This is because even for nested structures, a Description
object is returned.

It is possible to create arrays that immitate nested table-like structure with _v_nested_descr attribute::

    >>> import numpy
    >>> table.description._v_nested_descr
    [('info2', [('info3', [('x', '()f8'), ('y', '()u1')]), ('name', '()S10'),
     ('value', '()f8')]), ('info1', [('name', '()S10'), ('value', '()f8')]),
     ('color', '()u4')]
    >>> numpy.rec.array(None, shape=0,
                        dtype=table.description._v_nested_descr)
    recarray([],
          dtype=[('info2', [('info3', [('x', '>f8'), ('y', '|u1')]),
                 ('name', '|S10'), ('value', '>f8')]),
                 ('info1', [('name', '|S10'), ('value', '>f8')]),
                 ('color', '>u4')])
    >>> numpy.rec.array(None, shape=0,
                        dtype=table.description.info2._v_nested_descr)
    recarray([],
          dtype=[('info3', [('x', '>f8'), ('y', '|u1')]), ('name', '|S10'),
                 ('value', '>f8')])

Last but not least, there is a special iterator of the Description class: _f_walk,
which returns different columns of the table::

    >>> for coldescr in table.description._f_walk():
    ...     print(f"column--> {coldescr}")
    column--> Description([('info2', [('info3', [('x', '()f8'), ('y', '()u1')]),
                           ('name', '()S10'), ('value', '()f8')]),
                           ('info1', [('name', '()S10'), ('value', '()f8')]),
                           ('color', '()u4')])
    column--> EnumCol(enum=Enum({'blue': 2, 'green': 1, 'red': 0}), dflt='red',
                                base=UInt32Atom(shape=(), dflt=0), shape=(), pos=2)
    column--> Description([('info3', [('x', '()f8'), ('y', '()u1')]), ('name', '()S10'),
                           ('value', '()f8')])
    column--> StringCol(itemsize=10, shape=(), dflt='', pos=1)
    column--> Float64Col(shape=(), dflt=0.0, pos=2)
    column--> Description([('name', '()S10'), ('value', '()f8')])
    column--> StringCol(itemsize=10, shape=(), dflt='', pos=0)
    column--> Float64Col(shape=(), dflt=0.0, pos=1)
    column--> Description([('x', '()f8'), ('y', '()u1')])
    column--> Float64Col(shape=(), dflt=1.0, pos=0)
    column--> UInt8Col(shape=(), dflt=1, pos=1)

See the :ref:`DescriptionClassDescr` for a complete listing of attributes
and methods of the Description object.

Well, this is the end of this tutorial. As always, remember to close
your files::

    >>> fileh.close()

Finally, you may want to have a look at your resulting data file.

.. code-block:: bash

    $ ptdump -d nested-tut.h5
    / (RootGroup) ''
    /table (Table(13,)) ''
      Data dump:
    [0] (((1.0, 0), 'name2-0', 0.0), ('name1-0', 0.0), 0L)
    [1] (((1.0, 1), 'name2-1', 0.0), ('name1-1', 0.0), 1L)
    [2] (((1.0, 2), 'name2-2', 0.0), ('name1-2', 0.0), 2L)
    [3] (((1.0, 3), 'name2-3', 0.0), ('name1-3', 0.0), 0L)
    [4] (((1.0, 4), 'name2-4', 0.0), ('name1-4', 0.0), 1L)
    [5] (((1.0, 5), 'name2-5', 0.0), ('name1-5', 0.0), 2L)
    [6] (((1.0, 6), 'name2-6', 0.0), ('name1-6', 0.0), 0L)
    [7] (((1.0, 7), 'name2-7', 0.0), ('name1-7', 0.0), 1L)
    [8] (((1.0, 8), 'name2-8', 0.0), ('name1-8', 0.0), 2L)
    [9] (((1.0, 9), 'name2-9', 0.0), ('name1-9', 0.0), 0L)
    [10] (((1.0, 0), 'name2-0', 0.0), ('name1-0', 0.0), 0L)
    [11] (((1.0, 4), 'name2-4', 0.0), ('name1-4', 0.0), 1L)
    [12] (((1.0, 8), 'name2-8', 0.0), ('name1-8', 0.0), 2L)
    /table2 (Table(3,)) ''
      Data dump:
    [0] (((1.0, 0), 'name2-0', 0.0), ('name1-0', 0.0), 0L)
    [1] (((1.0, 4), 'name2-4', 0.0), ('name1-4', 0.0), 1L)
    [2] (((1.0, 8), 'name2-8', 0.0), ('name1-8', 0.0), 2L)

Most of the code in this section is also available in
:file:`examples/nested-tut.py`.

PyTables provides a comprehensive set of tools to work with nested structures and to address your classification needs.
Try to avoid nesting your data too deeply as it may lead to very long and convoluted series of lists, tuples, and description objects, which can be hard to read and understand.


Other examples in PyTables distribution
---------------------------------------
Feel free to examine the rest of examples in directory :file:`examples/`, and
try to understand them. We have written several practical sample scripts to
give you an idea of the PyTables capabilities, its way of dealing with HDF5
objects, and how it can be used in the real world.

------------

.. [1] Appending data to arrays is also supported, but you need to create
       special objects called EArray (see :ref:`EArrayClassDescr` for more
       info).

.. [2] Note that you can append not only scalar values to tables,
       but also fully multidimensional array objects.

.. [3] You can even *hide* nodes temporarily. Can you think of a way to do it?

.. [4] In fact, only integer values are supported right now, but
       this may change in the future.


.. _datatypes:

Supported data types in PyTables
================================

All PyTables datasets can handle the complete set of data types supported by
the NumPy (see :ref:`[NUMPY] <NUMPY>`) package in Python.
The data types for table fields can be set via instances of the Col class and
its descendants (see :ref:`ColClassDescr`), while the data type of array
elements can be set through the use of the Atom class and its descendants
(see :ref:`AtomClassDescr`).

PyTables uses ordinary strings to represent its *types*, with most of them
matching the names of NumPy scalar types. Usually, a PyTables type consists
of two parts: a *kind* and a *precision* in bits.
The precision may be omitted in types with just one supported precision (like
bool) or with a non-fixed size (like string).

There are eight kinds of types supported by PyTables:

- bool: Boolean (true/false) types.
  Supported precisions: 8 (default) bits.

- int: Signed integer types.
  Supported precisions: 8, 16, 32 (default) and 64 bits.

- uint: Unsigned integer types.
  Supported precisions: 8, 16, 32 (default) and 64 bits.

- float: Floating point types.
  Supported precisions: 16, 32, 64 (default) bits and extended precision
  floating point (see
  :ref:`note on floating point types<floating-point-note>`).

- complex: Complex number types.
  Supported precisions: 64 (32+32), 128 (64+64, default) bits and extended
  precision complex (see
  :ref:`note on floating point types<floating-point-note>`).

- string: Raw string types.
  Supported precisions: 8-bit positive multiples.

- time: Data/time types.
  Supported precisions: 32 and 64 (default) bits.

- enum: Enumerated types.
  Precision depends on base type.

.. _floating-point-note:
.. note:: Floating point types.

   The half precision floating point data type (float16) and extended
   precision ones (fload96, float128, complex192, complex256) are only
   available if numpy_ supports them on the host platform.

   Also, in order to use the half precision floating point type (float16)
   it is required numpy_ >= 1.6.0.

.. _numpy: http://www.numpy.org

The time and enum kinds area little bit special, since they represent HDF5
types which have no direct Python counterpart, though atoms of these kinds
have a more-or-less equivalent NumPy data type.

There are two types of time: 4-byte signed integer (time32) and 8-byte double
precision floating point (time64). Both of them reflect the number of seconds
since the Unix epoch, i.e. Jan 1 00:00:00 UTC 1970. They are stored in memory
as NumPy's int32 and float64, respectively, and in the HDF5 file using the
H5T_TIME class. Integer times are stored on disk as such, while floating
point times are split into two signed integer values representing seconds and
microseconds (beware: smaller decimals will be lost!).

PyTables also supports HDF5 H5T_ENUM *enumerations* (restricted sets of
unique name and unique value pairs). The NumPy representation of an
enumerated value (an Enum, see :ref:`EnumClassDescr`) depends on the concrete
*base type* used to store the enumeration in the HDF5 file.
Currently, only scalar integer values (both signed and unsigned) are
supported in enumerations. This restriction may be lifted when HDF5 supports
other kinds on enumerated values.

Here you have a quick reference to the complete set of supported data types:

.. table:: **Data types supported for array elements and tables columns in
             PyTables.**

    ================== ========================== ====================== =============== ==================
    Type Code          Description                C Type                 Size (in bytes) Python Counterpart
    ================== ========================== ====================== =============== ==================
    bool               boolean                    unsigned char          1               bool
    int8               8-bit integer              signed char            1               int
    uint8              8-bit unsigned integer     unsigned char          1               int
    int16              16-bit integer             short                  2               int
    uint16             16-bit unsigned integer    unsigned short         2               int
    int32              integer                    int                    4               int
    uint32             unsigned integer           unsigned int           4               long
    int64              64-bit integer             long long              8               long
    uint64             unsigned 64-bit integer    unsigned long long     8               long
    float16 [1]_       half-precision float       -                      2               -
    float32            single-precision float     float                  4               float
    float64            double-precision float     double                 8               float
    float96 [1]_ [2]_  extended precision float   -                      12              -
    float128 [1]_ [2]_ extended precision float   -                      16              -
    complex64          single-precision complex   struct {float r, i;}   8               complex
    complex128         double-precision complex   struct {double r, i;}  16              complex
    complex192 [1]_    extended precision complex -                      24              -
    complex256 [1]_    extended precision complex -                      32              -
    string             arbitrary length string    char[]                 *               str
    time32             integer time               POSIX's time_t         4               int
    time64             floating point time        POSIX's struct timeval 8               float
    enum               enumerated value           enum                   -               -
    ================== ========================== ====================== =============== ==================

.. rubric:: Footnotes

.. [1] see the above :ref:`note on floating point types <floating-point-note>`.
.. [2] currently in numpy_. "float96" and "float128" are equivalent of
       "longdouble" i.e. 80 bit extended precision floating point.
Installation
============

.. epigraph::

    Make things as simple as possible, but not any simpler.

    -- Albert Einstein


The Python Distutils are used to build and install PyTables, so it is fairly
simple to get the application up and running. If you want to install the
package from sources you can go on reading to the next section.

However, if you want to go straight to binaries that 'just work' for the main
platforms (Linux, Mac OSX and Windows), you might want to use the excellent
Anaconda_, ActivePython_, Canopy_ distributions.  PyTables usually distributes its own
Windows binaries too; go :ref:`binaryInstallationDescr` for instructions.
Finally `Christoph Gohlke`_ also maintains an excellent suite of a variety of
binary packages for Windows at his site.

.. _Anaconda: https://store.continuum.io/cshop/anaconda/
.. _Canopy: https://www.enthought.com/products/canopy/
.. _ActivePython: https://www.activestate.com/activepython/downloads
.. _`Christoph Gohlke`: http://www.lfd.uci.edu/~gohlke/pythonlibs/


Installation from source
------------------------

These instructions are for both Unix/MacOS X and Windows systems. If you are
using Windows, it is assumed that you have a recent version of MS Visual C++
compiler installed.
A GCC compiler is assumed for Unix, but other compilers should work as well.

Extensions in PyTables have been developed in Cython (see
:ref:`[CYTHON] <CYTHON>`) and the C language. You can rebuild everything from
scratch if you have Cython installed, but this is not necessary, as the Cython
compiled source is included in the source distribution.

To compile PyTables you will need a recent version of Python, the HDF5 (C
flavor) library from http://www.hdfgroup.org, and the NumPy (see
:ref:`[NUMPY] <NUMPY>`) and Numexpr (see :ref:`[NUMEXPR] <NUMEXPR>`)
packages.


Prerequisites
~~~~~~~~~~~~~

First, make sure that you have

* Python_ >= 3.6 (PyTables-3.5 was the last release with Python 2.7 support)
* HDF5_ >= 1.8.4 (>=1.8.15 is strongly recommended)
* NumPy_ >= 1.19.0
* Numexpr_ >= 2.6.2
* Cython_ >= 0.29.21
* c-blosc_ >= 1.4.1 (sources are bundled with PyTables sources but the user can
  use an external version of sources using the :envvar:`BLOSC_DIR` environment
  variable or the `--blosc` flag of the :file:`setup.py`)

installed (for testing purposes, we are using HDF5_ 1.8.15, NumPy_ 1.10.2
and Numexpr_ 2.5.2 currently). If you don't, fetch and install them before
proceeding.

.. _Python: http://www.python.org
.. _HDF5: http://www.hdfgroup.org/HDF5
.. _NumPy: http://www.numpy.org
.. _Numexpr: http://code.google.com/p/numexpr
.. _Cython: http://www.cython.org
.. _c-blosc: http://blosc.org

Compile and install these packages (but see :ref:`prerequisitesBinInst` for
instructions on how to install pre-compiled binaries if you are not willing
to compile the prerequisites on Windows systems).

For compression (and possibly improved performance), you will need to install
the Zlib (see :ref:`[ZLIB] <ZLIB>`), which is also required by HDF5 as well.
You may also optionally install the excellent LZO compression library (see
:ref:`[LZO] <LZO>` and :ref:`compressionIssues`). The high-performance bzip2
compression library can also be used with PyTables (see
:ref:`[BZIP2] <BZIP2>`).

The Blosc (see :ref:`[BLOSC] <BLOSC>`) compression library is embedded
in PyTables, so this will be used in case it is not found in the
system.  So, in case the installer warns about not finding it, do not
worry too much ;)

**Unix**

    setup.py will detect HDF5, Blosc, LZO, or bzip2 libraries and include
    files under :file:`/usr` or :file:`/usr/local`; this will cover most
    manual installations as well as installations from packages.  If setup.py
    can not find libhdf5, libhdf5 (or liblzo, or libbz2 that you may wish to
    use) or if you have several versions of a library installed and want to
    use a particular one, then you can set the path to the resource in the
    environment, by setting the values of the :envvar:`HDF5_DIR`,
    :envvar:`LZO_DIR`, :envvar:`BZIP2_DIR` or :envvar:`BLOSC_DIR` environment
    variables to the path to the particular resource. You may also specify the
    locations of the resource root directories on the setup.py command line.
    For example::

        --hdf5=/stuff/hdf5-1.8.12
        --blosc=/stuff/blosc-1.8.1
        --lzo=/stuff/lzo-2.02
        --bzip2=/stuff/bzip2-1.0.5

    If your HDF5 library was built as a shared library not in the runtime load
    path, then you can specify the additional linker flags needed to find the
    shared library on the command line as well. For example::

        --lflags="-Xlinker -rpath -Xlinker /stuff/hdf5-1.8.12/lib"

    You may also want to try setting the :envvar:`LD_LIBRARY_PATH`
    environment variable to point to the directory where the shared libraries
    can be found. Check your compiler and linker documentation as well as the
    Python Distutils documentation for the correct syntax or environment
    variable names.
    It is also possible to link with specific libraries by setting the
    :envvar:`LIBS` environment variable::

        LIBS="hdf5-1.8.12 nsl"

    Starting from PyTables 3.2 can also query the *pkg-config* database to
    find the required packages. If available, pkg-config is used by default
    unless explicitly disabled.

    To suppress the use of *pkg-config*::

      $ python3 setup.py build --use-pkgconfig=FALSE

    or use the :envvar:`USE-PKGCONFIG` environment variable::

      $ env USE_PKGCONFIG=FALSE python3 setup.py build

**Windows**

    You can get ready-to-use Windows binaries and other development files for
    most of the following libraries from the GnuWin32 project (see
    :ref:`[GNUWIN32] <GNUWIN32>`).  In case you cannot find the LZO binaries
    in the GnuWin32 repository, you can find them at
    http://sourceforge.net/projects/pytables/files/lzo-win.
    Once you have installed the prerequisites, setup.py needs to know where
    the necessary library *stub* (.lib) and *header* (.h) files are installed.
    You can set the path to the include and dll directories for the HDF5
    (mandatory) and LZO, BZIP2, BLOSC (optional) libraries in the environment,
    by setting the values of the :envvar:`HDF5_DIR`, :envvar:`LZO_DIR`,
    :envvar:`BZIP2_DIR` or :envvar:`BLOSC_DIR` environment variables to the
    path to the particular resource.  For example::

        set HDF5_DIR=c:\\stuff\\hdf5-1.8.5-32bit-VS2008-IVF101\\release
        set BLOSC_DIR=c:\\Program Files (x86)\\Blosc
        set LZO_DIR=c:\\Program Files (x86)\\GnuWin32
        set BZIP2_DIR=c:\\Program Files (x86)\\GnuWin32

    You may also specify the locations of the resource root directories on the
    setup.py command line.
    For example::

        --hdf5=c:\\stuff\\hdf5-1.8.5-32bit-VS2008-IVF101\\release
        --blosc=c:\\Program Files (x86)\\Blosc
        --lzo=c:\\Program Files (x86)\\GnuWin32
        --bzip2=c:\\Program Files (x86)\\GnuWin32

**Conda**

    Pre-built packages for PyTables are available in the anaconda (default)
    channel::

        conda install pytables

    The most recent version is usually available in the conda-forge
    channel::

        conda config --add channels conda-forge
        conda install pytables

    The HDF5 libraries and other helper packages are automatically found in
    a conda environment. During installation setup.py uses the `CONDA_PREFIX`
    environment variable to detect a conda environment. If detected it will
    try to find all packages within this environment. PyTables needs at least
    the hdf5 package::

        conda install hdf5
        python3 setup.py install

    It is still possible to override package locations using the
    :envvar:`HDF5_DIR`, :envvar:`LZO_DIR`, :envvar:`BZIP2_DIR` or
    :envvar:`BLOSC_DIR` environment variables.

    When inside a conda environment *pkg-config* will not work. To disable
    using the conda environment and fall back to *pkg-config* use `--no-conda`::

          python3 setup.py install --no-conda

    When the `--use-pkgconfig` flag is used, `--no-conda` is assumed.

**Development version (Unix)**

    Installation of the development version is very similar to installation
    from a source package (described above).  There are two main differences:

    #. sources have to be downloaded from the `PyTables source repository`_
       hosted on GitHub_. Git (see :ref:`[GIT] <GIT>`) is used as VCS.
       The following command create a local copy of latest development version
       sources::

        $ git clone --recursive https://github.com/PyTables/PyTables.git

    #. sources in the git repository do not include pre-built documentation
       and pre-generated C code of Cython extension modules.  To be able to
       generate them, both Cython (see :ref:`[CYTHON] <CYTHON>`) and
       sphinx >= 1.0.7 (see :ref:`[SPHINX] <SPHINX>`) are mandatory
       prerequisites.

.. _`PyTables source repository`: https://github.com/PyTables/PyTables
.. _GitHub: https://github.com


PyTables package installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you have installed the HDF5 library and the NumPy and Numexpr packages,
you can proceed with the PyTables package itself.

#. Run this command from the main PyTables distribution directory, including
   any extra command line arguments as discussed above::

      $ python3 setup.py build

   If the HDF5 installation is in a custom path, e.g. $HOME/hdf5-1.8.15pre7,
   one of the following commands can be used::

      $ python3 setup.py build --hdf5=$HOME/hdf5-1.8.15pre7

   .. note::

       AVX2 support is detected automatically for your machine and, if found,
       it is enabled by default.  In some situations you may want to disable
       AVX2 explicitly (maybe your binaries have to be exported and run on
       machines that do not have AVX2 support).  In that case, define the
       DISABLE_AVX2 environment variable::

          $ DISABLE_AVX2=True python3 setup.py build  # for bash and its variants

#. To run the test suite, execute any of these commands.

   **Unix**
      In the sh shell and its variants::

        $ cd build/lib.linux-x86_64-3.3
        $ env PYTHONPATH=. python3 tables/tests/test_all.py

      or, if you prefer::

        $ cd build/lib.linux-x86_64-3.3
        $ env PYTHONPATH=. python3 -c "import tables; tables.test()"

      .. note::

          the syntax used above overrides original contents of the
          :envvar:`PYTHONPATH` environment variable.
          If this is not the desired behaviour and the user just wants to add
          some path before existing ones, then the safest syntax to use is
          the following::

            $ env PYTHONPATH=.${PYTHONPATH:+:$PYTHONPATH} python3 tables/tests/test_all.py

          Please refer to your :program:`sh` documentation for details.

   **Windows**

      Open the command prompt (cmd.exe or command.com) and type::

        > cd build\\lib.linux-x86_64-2.7
        > set PYTHONPATH=.;%PYTHONPATH%
        > python3 tables\\tests\\test_all.py

      or::

        > cd build\\lib.linux-x86_64-2.7
        > set PYTHONPATH=.;%PYTHONPATH%
        > python3 -c "import tables; tables.test()"

   Both commands do the same thing, but the latter still works on an already
   installed PyTables (so, there is no need to set the :envvar:`PYTHONPATH`
   variable for this case).
   However, before installation, the former is recommended because it is
   more flexible, as you can see below.
   If you would like to see verbose output from the tests simply add the
   `-v` flag and/or the word verbose to the first of the command lines
   above. You can also run only the tests in a particular test module.
   For example, to execute just the test_types test suite, you only have to
   specify it::

      # change to backslashes for win
      $ python3 tables/tests/test_types.py -v

   You have other options to pass to the :file:`test_all.py` driver::

      # change to backslashes for win
      $ python3 tables/tests/test_all.py --heavy

   The command above runs every test in the test unit. Beware, it can take a
   lot of time, CPU and memory resources to complete::

      # change to backslashes for win
      $ python3 tables/tests/test_all.py --print-versions

   The command above shows the versions for all the packages that PyTables
   relies on. Please be sure to include this when reporting bugs::

      # only under Linux 2.6.x
      $ python3 tables/tests/test_all.py --show-memory

   The command above prints out the evolution of the memory consumption after
   each test module completion. It's useful for locating memory leaks in
   PyTables (or packages behind it). Only valid for Linux 2.6.x kernels.
   And last, but not least, in case a test fails, please run the failing test
   module again and enable the verbose output::

      $ python3 tables/tests/test_<module>.py -v verbose

   and, very important, obtain your PyTables version information by using the
   `--print-versions` flag (see above) and send back both outputs to
   developers so that we may continue improving PyTables.
   If you run into problems because Python can not load the HDF5 library or
   other shared libraries.

   **Unix**

      Try setting the LD_LIBRARY_PATH or equivalent environment variable to
      point to the directory where the missing libraries can be found.

   **Windows**

      Put the DLL libraries (hdf5dll.dll and, optionally, lzo1.dll,
      bzip2.dll or blosc.dll) in a directory listed in your
      :envvar:`PATH` environment variable. The setup.py installation
      program will print out a warning to that effect if the libraries
      can not be found.

#. To install the entire PyTables Python package, change back to the root
   distribution directory and run the following command (make sure you have
   sufficient permissions to write to the directories where the PyTables files
   will be installed)::

      $ python3 setup.py install

   Again if one needs to point to libraries installed in custom paths, then
   specific setup.py options can be used::

      $ python3 setup.py install --hdf5=/hdf5/custom/path

   or::

      $ env HDF5_DIR=/hdf5/custom/path python3 setup.py install

   Of course, you will need super-user privileges if you want to install
   PyTables on a system-protected area. You can select, though, a different
   place to install the package using the `--prefix` flag::

      $ python3 setup.py install --prefix="/home/myuser/mystuff"

   Have in mind, however, that if you use the `--prefix` flag to
   install in a non-standard place, you should properly setup your
   :envvar:`PYTHONPATH` environment variable, so that the Python interpreter
   would be able to find your new PyTables installation.
   You have more installation options available in the Distutils package.
   Issue a::

      $ python3 setup.py install --help

   for more information on that subject.

That's it! Now you can skip to the next chapter to learn how to use PyTables.


Installation with :program:`pip`
--------------------------------

Many users find it useful to use the :program:`pip` program (or similar ones)
to install python packages.

As explained in previous sections the user should in any case ensure that all
dependencies listed in the `Prerequisites`_ section are correctly installed.

The simplest way to install PyTables using :program:`pip` is the following::

  $ python3 -m pip install tables

The following example shows how to install the latest stable version of
PyTables in the user folder when a older version of the package is already
installed at system level::

  $ python3 -m pip install --user --upgrade tables

The `--user` option tells to the :program:`pip` tool to install the package in
the user folder (``$HOME/.local`` on GNU/Linux and Unix systems), while the
`--upgrade` option forces the installation of the latest version even if an
older version of the package is already installed.

Additional options for the setup.py script can be specified using them
`--install-option`::

  $ python3 -m pip install --install-option='--hdf5=/custom/path/to/hdf5' tables

or::

  $ env HDF5_DIR=/custom/path/to/hdf5 python3 -m pip install tables

The :program:`pip` tool can also be used to install packages from a source
tar-ball::

  $ python3 -m pip install tables-3.0.0.tar.gz

To install the development version of PyTables from the *develop* branch of
the main :program:`git` :ref:`[GIT] <GIT>` repository the command is the
following::

  $ python3 -m pip install git+https://github.com/PyTables/PyTables.git@develop#egg=tables

A similar command can be used to install a specific tagged version::

  $ python3 -m pip install git+https://github.com/PyTables/PyTables.git@v.2.4.0#egg=tables

Of course the `pip` can be used to install only python packages.
Other dependencies like the HDF5 library of compression libraries have to
be installed by the user.

.. note::

   Recent versions of Debian_ and Ubuntu_ the HDF5 library is installed in
   with a very peculiar layout that allows to have both the serial and MPI
   versions installed at the same time.

   PyTables >= 3.2 natively supports the new layout via *pkg-config* (that
   is expected to be installed on the system at build time).

   If *pkg-config* is not available or PyTables is older than version 3.2,
   then the following command can be used::

     $ env CPPFLAGS=-I/usr/include/hdf5/serial \
     LDFLAGS=-L/usr/lib/x86_64-linux-gnu/hdf5/serial python3 setup.py install

   or::

     $ env CPPFLAGS=-I/usr/include/hdf5/serial \
     LDFLAGS=-L/usr/lib/x86_64-linux-gnu/hdf5/serial python3 -m pip install tables

.. _Debian: https://www.debian.org
.. _Ubuntu: http://www.ubuntu.com


.. _binaryInstallationDescr:

Binary installation (Windows)
-----------------------------

This section is intended for installing precompiled binaries on Windows
platforms. Binaries are distribution in wheel format, which can be downloaded
and installed using pip as described above. You may also find it useful for
instructions on how to install *binary prerequisites* even if you want to
compile PyTables itself on Windows.

.. _prerequisitesBinInst:

Windows prerequisites
~~~~~~~~~~~~~~~~~~~~~

First, make sure that you have Python 3, NumPy 1.8.0 and Numexpr 2.5.2 or
higher installed.

To enable compression with the optional LZO library (see the
:ref:`compressionIssues` for hints about how it may be used to improve
performance), fetch and install the LZO from
http://sourceforge.net/projects/pytables/files/lzo-win (choose v1.x for
Windows 32-bit and v2.x for Windows 64-bit).
Normally, you will only need to fetch that package and copy the included
lzo1.dll/lzo2.dll file in a directory in the PATH environment variable
(for example C:\\WINDOWS\\SYSTEM) or
python_installation_path\\Lib\\site-packages\\tables (the last directory may
not exist yet, so if you want to install the DLL there, you should do so
*after* installing the PyTables package), so that it can be found by the
PyTables extensions.

Please note that PyTables has internal machinery for dealing with uninstalled
optional compression libraries, so, you don't need to install the LZO or bzip2
dynamic libraries if you don't want to.


PyTables package installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On PyPI wheels for 32 and 64-bit versions of Windows and are usually provided. They
are automatically found and installed using pip::

    $ python3 -m pip install tables

If a matching wheel cannot be found for your installation, third party built wheels
can be found e.g. at the `Unofficial Windows Binaries for Python Extension Packages
<http://www.lfd.uci.edu/~gohlke/pythonlibs/#pytables>`_ page. Download the wheel
matching the version of python and either the 32 or 64-bit version and install
using pip::

    # python 3.6 64-bit:
    $ python3 -m pip install tables-3.6.1-2-cp36-cp36m-win_amd64.whl

You can (and *you should*) test your installation by running the next
commands::

    >>> import tables
    >>> tables.test()

on your favorite python shell. If all the tests pass (possibly with a few
warnings, related to the potential unavailability of LZO lib) you already have
a working, well-tested copy of PyTables installed! If any test fails, please
copy the output of the error messages as well as the output of::

    >>> tables.print_versions()

and mail them to the developers so that the problem can be fixed in future
releases.

You can proceed now to the next chapter to see how to use PyTables.
PyTables File Format
====================
PyTables has a powerful capability to deal with native HDF5 files created
with another tools. However, there are situations were you may want to create
truly native PyTables files with those tools while retaining fully
compatibility with PyTables format. That is perfectly possible, and in this
appendix is presented the format that you should endow to your own-generated
files in order to get a fully PyTables compatible file.

We are going to describe the *2.0 version of PyTables file format*
(introduced in PyTables version 2.0). As time goes by, some changes might be
introduced (and documented here) in order to cope with new necessities.
However, the changes will be carefully pondered so as to ensure backward
compatibility whenever is possible.

A PyTables file is composed with arbitrarily large amounts of HDF5 groups
(Groups in PyTables naming scheme) and datasets (Leaves in PyTables naming
scheme). For groups, the only requirements are that they must have some
*system attributes* available. By convention, system attributes in PyTables
are written in upper case, and user attributes in lower case but this is not
enforced by the software. In the case of datasets, besides the mandatory
system attributes, some conditions are further needed in their storage
layout, as well as in the datatypes used in there, as we will see shortly.

As a final remark, you can use any filter as you want to create a PyTables
file, provided that the filter is a standard one in HDF5, like *zlib*,
*shuffle* or *szip* (although the last one can not be used from within
PyTables to create a new file, datasets compressed with szip can be read,
because it is the HDF5 library which do the decompression transparently).


.. currentmodule:: tables

Mandatory attributes for a File
-------------------------------
The File object is, in fact, an special HDF5 *group* structure that is *root*
for the rest of the objects on the object tree. The next attributes are
mandatory for the HDF5 *root group* structure in PyTables files:

* *CLASS*: This attribute should always be set to 'GROUP' for group
  structures.
* *PYTABLES_FORMAT_VERSION*: It represents the internal format version, and
  currently should be set to the '2.0' string.
* *TITLE*: A string where the user can put some description on what is this
  group used for.
* *VERSION*: Should contains the string '1.0'.


Mandatory attributes for a Group
--------------------------------
The next attributes are mandatory for *group* structures:

* *CLASS*: This attribute should always be set to 'GROUP' for group structures.
* *TITLE*: A string where the user can put some description on what is this
  group used for.
* *VERSION*: Should contains the string '1.0'.


Optional attributes for a Group
-------------------------------
The next attributes are optional for *group* structures:

* *FILTERS*: When present, this attribute contains the filter properties (a
  Filters instance, see section :ref:`FiltersClassDescr`) that may be
  inherited by leaves or groups created immediately under this group. This is
  a packed 64-bit integer structure, where

  - *byte 0* (the least-significant byte) is the compression level
    (complevel).
  - *byte 1* is the compression library used (complib): 0 when irrelevant, 1
    for Zlib, 2 for LZO and 3 for Bzip2.
  - *byte 2* indicates which parameterless filters are enabled (shuffle and
    fletcher32): bit 0 is for *Shuffle* while bit 1 is for*Fletcher32*.
  - other bytes are reserved for future use.


Mandatory attributes, storage layout and supported data types for Leaves
------------------------------------------------------------------------
This depends on the kind of Leaf. The format for each type follows.


.. _TableFormatDescr:

Table format
~~~~~~~~~~~~

Mandatory attributes
^^^^^^^^^^^^^^^^^^^^
The next attributes are mandatory for *table* structures:

* *CLASS*: Must be set to 'TABLE'.
* *TITLE*: A string where the user can put some description on what is this
  dataset used for.
* *VERSION*: Should contain the string '2.6'.
* *FIELD_X_NAME*: It contains the names of the different fields. The X means
  the number of the field, zero-based (beware, order do matter). You should
  add as many attributes of this kind as fields you have in your records.
* *FIELD_X_FILL*: It contains the default values of the different fields. All
  the datatypes are supported natively, except for complex types that are
  currently serialized using Pickle.  The X means the number of the field,
  zero-based (beware, order do matter). You should add as many attributes of
  this kind as fields you have in your records.  These fields are meant for
  saving the default values persistently and their existence is optional.
* *NROWS*: This should contain the number of *compound* data type entries in
  the dataset. It must be an *int* data type.


Storage Layout
^^^^^^^^^^^^^^
A Table has a *dataspace* with a *1-dimensional chunked* layout.

Datatypes supported
^^^^^^^^^^^^^^^^^^^
The datatype of the elements (rows) of Table must be the H5T_COMPOUND
*compound* data type, and each of these compound components must be built
with only the next HDF5 data types *classes*:

* *H5T_BITFIELD*: This class is used to represent the Bool type. Such a type
  must be build using a H5T_NATIVE_B8 datatype, followed by a HDF5
  H5Tset_precision call to set its precision to be just 1 bit.
* *H5T_INTEGER*: This includes the next data types:
    * *H5T_NATIVE_SCHAR*: This represents a *signed char* C type, but it is
      effectively used to represent an Int8 type.
    * *H5T_NATIVE_UCHAR*:  This represents an *unsigned char* C type, but it
      is effectively used to represent an UInt8 type.
    * *H5T_NATIVE_SHORT*: This represents a *short* C type, and it is
      effectively used to represent an Int16 type.
    * *H5T_NATIVE_USHORT*: This represents an *unsigned short* C type, and it
      is effectively used to represent an UInt16 type.
    * *H5T_NATIVE_INT*: This represents an *int* C type, and it is
      effectively used to represent an Int32 type.
    * *H5T_NATIVE_UINT*: This represents an *unsigned int* C type, and it is
      effectively used to represent an UInt32 type.
    * *H5T_NATIVE_LONG*: This represents a *long* C type, and it is
      effectively used to represent an Int32 or an Int64, depending on
      whether you are running a 32-bit or 64-bit architecture.
    * *H5T_NATIVE_ULONG*: This represents an *unsigned long* C type, and it
      is effectively used to represent an UInt32 or an UInt64, depending on
      whether you are running a 32-bit or 64-bit architecture.
    * *H5T_NATIVE_LLONG*: This represents a *long long* C type (__int64, if
      you are using a Windows system) and it is effectively used to represent
      an Int64 type.
    * *H5T_NATIVE_ULLONG*: This represents an *unsigned long long* C type
      (beware: this type does not have a correspondence on Windows systems)
      and it is effectively used to represent an UInt64 type.
* *H5T_FLOAT*: This includes the next datatypes:
    * *H5T_NATIVE_FLOAT*: This represents a *float* C type and it is
      effectively used to represent an Float32 type.
    * *H5T_NATIVE_DOUBLE*: This represents a *double* C type and it is
      effectively used to represent an Float64 type.
* *H5T_TIME*: This includes the next datatypes:
    * *H5T_UNIX_D32*: This represents a POSIX *time_t* C type and it is
      effectively used to represent a 'Time32' aliasing type, which
      corresponds to an Int32 type.
    * *H5T_UNIX_D64*: This represents a POSIX *struct timeval* C type and it
      is effectively used to represent a 'Time64' aliasing type, which
      corresponds to a Float64 type.
* *H5T_STRING*: The datatype used to describe strings in PyTables is H5T_C_S1
  (i.e. a *string* C type) followed with a call to the HDF5 H5Tset_size()
  function to set their length.
* *H5T_ARRAY*: This allows the construction of homogeneous, multidimensional
  arrays, so that you can include such objects in compound records. The types
  supported as elements of H5T_ARRAY data types are the ones described above.
  Currently, PyTables does not support nested H5T_ARRAY types.
* *H5T_COMPOUND*: This allows the support for datatypes that are compounds of
  compounds (this is also known as *nested types* along this manual).

  This support can also be used for defining complex numbers. Its format is
  described below:

  The H5T_COMPOUND type class contains two members. Both members must have
  the H5T_FLOAT atomic datatype class. The name of the first member should be
  "r" and represents the real part. The name of the second member should be
  "i" and represents the imaginary part. The *precision* property of both of
  the H5T_FLOAT members must be either 32 significant bits (e.g.
  H5T_NATIVE_FLOAT) or 64 significant bits (e.g. H5T_NATIVE_DOUBLE). They
  represent Complex32 and Complex64 types respectively.


Array format
~~~~~~~~~~~~

Mandatory attributes
^^^^^^^^^^^^^^^^^^^^
The next attributes are mandatory for *array* structures:

* *CLASS*: Must be set to 'ARRAY'.
* *TITLE*: A string where the user can put some description on what is this
  dataset used for.
* *VERSION*: Should contain the string '2.3'.


Storage Layout
^^^^^^^^^^^^^^
An Array has a *dataspace* with a *N-dimensional contiguous* layout (if you
prefer a *chunked* layout see EArray below).


Datatypes supported
^^^^^^^^^^^^^^^^^^^
The elements of Array must have either HDF5 *atomic* data types or a
*compound* data type representing a complex number. The atomic data types can
currently be one of the next HDF5 data type *classes*: H5T_BITFIELD,
H5T_INTEGER, H5T_FLOAT and H5T_STRING. The H5T_TIME class is also supported
for reading existing Array objects, but not for creating them. See the Table
format description in :ref:`TableFormatDescr` for more info about these
types.

In addition to the HDF5 atomic data types, the Array format supports complex
numbers with the H5T_COMPOUND data type class.
See the Table format description in :ref:`TableFormatDescr` for more info
about this special type.

You should note that H5T_ARRAY class datatypes are not allowed in Array
objects.


CArray format
~~~~~~~~~~~~~

Mandatory attributes
^^^^^^^^^^^^^^^^^^^^
The next attributes are mandatory for *CArray* structures:

* *CLASS*: Must be set to 'CARRAY'.
* *TITLE*: A string where the user can put some description on what is this
  dataset used for.
* *VERSION*: Should contain the string '1.0'.


Storage Layout
^^^^^^^^^^^^^^
An CArray has a *dataspace* with a *N-dimensional chunked* layout.

Datatypes supported
^^^^^^^^^^^^^^^^^^^
The elements of CArray must have either HDF5 *atomic* data types or a
*compound* data type representing a complex number. The atomic data types can
currently be one of the next HDF5 data type *classes*: H5T_BITFIELD,
H5T_INTEGER, H5T_FLOAT and H5T_STRING. The H5T_TIME class is also supported
for reading existing CArray objects, but not for creating them. See the Table
format description in :ref:`TableFormatDescr` for more info about these
types.

In addition to the HDF5 atomic data types, the CArray format supports complex
numbers with the H5T_COMPOUND data type class.
See the Table format description in :ref:`TableFormatDescr` for more info
about this special type.

You should note that H5T_ARRAY class datatypes are not allowed yet in Array
objects.


EArray format
~~~~~~~~~~~~~

Mandatory attributes
^^^^^^^^^^^^^^^^^^^^
The next attributes are mandatory for *earray* structures:

* *CLASS*: Must be set to 'EARRAY'.
* *EXTDIM*: (*Integer*) Must be set to the extendable dimension. Only one
  extendable dimension is supported right now.
* *TITLE*: A string where the user can put some description on what is this
  dataset used for.
* *VERSION*: Should contain the string '1.3'.


Storage Layout
^^^^^^^^^^^^^^
An EArray has a *dataspace* with a *N-dimensional chunked* layout.


Datatypes supported
^^^^^^^^^^^^^^^^^^^
The elements of EArray are allowed to have the same data types as for the
elements in the Array format. They can be one of the HDF5 *atomic* data type
*classes*: H5T_BITFIELD, H5T_INTEGER, H5T_FLOAT, H5T_TIME or H5T_STRING, see
the Table format description in :ref:`TableFormatDescr` for more info about
these types. They can also be a H5T_COMPOUND datatype representing a complex
number, see the Table format description in :ref:`TableFormatDescr`.

You should note that H5T_ARRAY class data types are not allowed in EArray
objects.


.. _VLArrayFormatDescr:

VLArray format
~~~~~~~~~~~~~~

Mandatory attributes
^^^^^^^^^^^^^^^^^^^^
The next attributes are mandatory for *vlarray* structures:

* *CLASS*: Must be set to 'VLARRAY'.
* *PSEUDOATOM*: This is used so as to specify the kind of pseudo-atom (see
  :ref:`VLArrayFormatDescr`) for the VLArray. It can take the values
  'vlstring', 'vlunicode' or 'object'. If your atom is not a pseudo-atom then
  you should not specify it.
* *TITLE*: A string where the user can put some description on what is this
  dataset used for.
* *VERSION*: Should contain the string '1.3'.


Storage Layout
^^^^^^^^^^^^^^
An VLArray has a *dataspace* with a *1-dimensional chunked* layout.


Data types supported
^^^^^^^^^^^^^^^^^^^^
The data type of the elements (rows) of VLArray objects must be the H5T_VLEN
*variable-length* (or VL for short) datatype, and the base datatype specified
for the VL datatype can be of any *atomic* HDF5 datatype that is listed in
the Table format description :ref:`TableFormatDescr`.  That includes the
classes:

- H5T_BITFIELD
- H5T_INTEGER
- H5T_FLOAT
- H5T_TIME
- H5T_STRING
- H5T_ARRAY

They can also be a H5T_COMPOUND data type representing a complex number, see
the Table format description in :ref:`TableFormatDescr` for a detailed
description.

You should note that this does not include another VL datatype, or a compound
datatype that does not fit the description of a complex number. Note as well
that, for object and vlstring pseudo-atoms, the base for the VL datatype is
always a H5T_NATIVE_UCHAR (H5T_NATIVE_UINT for vlunicode). That means that
the complete row entry in the dataset has to be used in order to fully
serialize the object or the variable length string.


Optional attributes for Leaves
------------------------------
The next attributes are optional for *leaves*:

* *FLAVOR*: This is meant to provide the information about the kind of object
  kept in the Leaf, i.e. when the dataset is read, it will be converted to
  the indicated flavor.
  It can take one the next string values:

    * *"numpy"*: Read data (structures arrays, arrays, records, scalars) will
      be returned as NumPy objects.
    * *"python"*: Read data will be returned as Python lists, tuples, or
      scalars.

.. _condition_syntax:

Condition Syntax
================
.. currentmodule:: tables

Conditions in PyTables are used in methods related with in-kernel and indexed
searches such as :meth:`Table.where` or :meth:`Table.read_where`.
They are interpreted using Numexpr, a powerful package for achieving C-speed
computation of array operations (see :ref:`[NUMEXPR] <NUMEXPR>`).

A condition on a table is just a *string* containing a Python expression
involving *at least one column*, and maybe some constants and external
variables, all combined with algebraic operators and functions. The result of
a valid condition is always a *boolean array* of the same length as the
table, where the *i*-th element is true if the value of the expression on the
*i*-th row of the table evaluates to true

That is the reason why multidimensional fields in a table are not supported
in conditions, since the truth value of each resulting multidimensional
boolean value is not obvious.
Usually, a method using a condition will only consider the rows where the
boolean result is true.

For instance, the condition 'sqrt(x*x + y*y) < 1' applied on a table with x
and y columns consisting of floating point numbers results in a boolean array
where the *i*-th element is true if (unsurprisingly) the value of the square
root of the sum of squares of x and y is less than 1.
The sqrt() function works element-wise, the 1 constant is adequately
broadcast to an array of ones of the length of the table for evaluation, and
the *less than* operator makes the result a valid boolean array. A condition
like 'mycolumn' alone will not usually be valid, unless mycolumn is itself a
column of scalar, boolean values.

In the previous conditions, mycolumn, x and y are examples of *variables*
which are associated with columns.
Methods supporting conditions do usually provide their own ways of binding
variable names to columns and other values. You can read the documentation of
:meth:`Table.where` for more information on that. Also, please note that the
names None, True and False, besides the names of functions (see below) *can
not be overridden*, but you can always define other new names for the objects
you intend to use.

Values in a condition may have the following types:

- 8-bit boolean (bool).

- 32-bit signed integer (int).

- 64-bit signed integer (long).

- 32-bit, single-precision floating point number (float or float32).

- 64-bit, double-precision floating point number (double or float64).

- 2x64-bit, double-precision complex number (complex).

- Raw string of bytes (str).

Nevertheless, if the type passed is not among the above ones, it will be
silently upcasted, so you don't need to worry too much about passing
supported types, except for the Unsigned 64 bits integer, that cannot be
upcasted to any of the supported types.

However, the types in PyTables conditions are somewhat stricter than those of
Python. For instance, the *only* valid constants for booleans are True and
False, and they are *never* automatically cast to integers. The type
strengthening also affects the availability of operators and functions.
Beyond that, the usual type inference rules apply.

Conditions support the set of operators listed below:

- Logical operators: &, \|, ~.

- Comparison operators: <, <=, ==, !=, >=, >.

- Unary arithmetic operators: -.

- Binary arithmetic operators: +, -, \*, /, \**, %.

Types do not support all operators. Boolean values only support logical and
strict (in)equality comparison operators, while strings only support
comparisons, numbers do not work with logical operators, and complex
comparisons can only check for strict (in)equality. Unsupported operations
(including invalid castings) raise NotImplementedError exceptions.

You may have noticed the special meaning of the usually bitwise operators &,
| and ~. Because of the way Python handles the short-circuiting of logical
operators and the truth values of their operands, conditions must use the
bitwise operator equivalents instead.
This is not difficult to remember, but you must be careful because bitwise
operators have a *higher precedence* than logical operators. For instance,
'a and b == c' (*a is true AND b is equal to c*) is *not* equivalent to
'a & b == c' (*a AND b is equal to c)*. The safest way to avoid confusions is
to *use parentheses* around logical operators, like this: 'a & (b == c)'.
Another effect of short-circuiting is that expressions like '0 < x < 1' will
*not* work as expected; you should use '(0 < x) & (x < 1)'.

All of this may be solved if Python supported overloadable boolean operators
(see PEP 335) or some kind of non-shortcircuiting boolean operators (like C's
&&, || and !).

You can also use the following functions in conditions:

- where(bool, number1, number2):
  number - number1 if the bool condition is true, number2 otherwise.

- {sin,cos,tan}(float|complex):
  float|complex - trigonometric sine, cosine or tangent.

- {arcsin,arccos,arctan}(float|complex):
  float|complex - trigonometric inverse sine, cosine or tangent.

- arctan2(float1, float2):
  float - trigonometric inverse tangent of float1/float2.

- {sinh,cosh,tanh}(float|complex):
  float|complex - hyperbolic sine, cosine or tangent.

- {arcsinh,arccosh,arctanh}(float|complex):
  float|complex - hyperbolic inverse sine, cosine or tangent.

- {log,log10,log1p}(float|complex):
  float|complex - natural, base-10 and log(1+x) logarithms.

- {exp,expm1}(float|complex):
  float|complex - exponential and exponential minus one.

- sqrt(float|complex): float|complex - square root.

- abs(float|complex): float|complex - absolute value.

- {real,imag}(complex):
  float - real or imaginary part of complex.

- complex(float, float):
  complex - complex from real and imaginary parts.

=====================
PyTables User's Guide
=====================
-------------------------------
Hierarchical datasets in Python
-------------------------------

:Authors: **Francesc Alted, Ivan Vilata, Scott Prater, Vicent Mas, Tom Hedley,
          Antonio Valentino, Jeffrey Whitaker, Anthony Scopatz, Josh Moore**
:Copyright: |copy| 2002, 2003, 2004 - Francesc Alted

            |copy| 2005, 2006, 2007 - Cárabos Coop. V.

            |copy| 2008, 2009, 2010 - Francesc Alted

            |copy| 2011–2021 - PyTables maintainers

--------
Contents
--------
.. toctree::
    :maxdepth: 1

    introduction
    installation
    tutorials
    libref
    optimization
    filenode
    datatypes
    condition_syntax
    parameter_files
    utilities
    file_format
    bibliography

--------------------------------------------------------
Copyright Notice and Statement for PyTables User's Guide
--------------------------------------------------------
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

a. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

b. Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the
distribution.

c. Neither the name of Francesc Alted nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

.. |copy|   unicode:: U+000A9 .. COPYRIGHT SIGN
.. _parameter_files:

PyTables parameter files
========================

.. currentmodule:: tables.parameters

PyTables issues warnings when certain limits are exceeded.  Those limits are
not intrinsic limitations of the underlying software, but rather are
proactive measures to avoid large resource consumptions.  The default limits
should be enough for most of cases, and users should try to respect them.
However, in some situations, it can be convenient to increase (or decrease)
these limits.

Also, and in order to get maximum performance, PyTables implements a series
of sophisticated features, like I/O buffers or different kind of caches (for
nodes, chunks and other internal metadata).  These features comes with a
default set of parameters that ensures a decent performance in most of
situations.  But, as there is always a need for every case, it is handy to
have the possibility to fine-tune some of these parameters.

Because of these reasons, PyTables implements a couple of ways to change the
values of these parameters.  All the *tunable* parameters live in the
:file:`tables/parameters.py`.  The user can choose to change them in the
parameter files themselves for a global and persistent change.  Moreover, if
he wants a finer control, he can pass any of these parameters directly to the
:func:`tables.open_file` function, and the new parameters will only take
effect in the corresponding file (the defaults will continue to be in the
parameter files).

A description of all of the tunable parameters follows.  As the defaults
stated here may change from release to release, please check with your actual
parameter files so as to know your actual default values.

.. warning::

    Changing the next parameters may have a very bad effect in the resource
    consumption and performance of your PyTables scripts.

    Please be careful when touching these!


Tunable parameters in parameters.py
-----------------------------------

Recommended maximum values
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autodata:: MAX_COLUMNS

.. autodata:: MAX_NODE_ATTRS

.. autodata:: MAX_GROUP_WIDTH

.. autodata:: MAX_TREE_DEPTH

.. autodata:: MAX_UNDO_PATH_LENGTH


Cache limits
~~~~~~~~~~~~
.. autodata:: CHUNK_CACHE_NELMTS

.. autodata:: CHUNK_CACHE_PREEMPT

.. autodata:: CHUNK_CACHE_SIZE

.. autodata:: COND_CACHE_SLOTS

.. autodata:: METADATA_CACHE_SIZE

.. autodata:: NODE_CACHE_SLOTS


Parameters for the different internal caches
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autodata:: BOUNDS_MAX_SIZE

.. autodata:: BOUNDS_MAX_SLOTS

.. autodata:: ITERSEQ_MAX_ELEMENTS

.. autodata:: ITERSEQ_MAX_SIZE

.. autodata:: ITERSEQ_MAX_SLOTS

.. autodata:: LIMBOUNDS_MAX_SIZE

.. autodata:: LIMBOUNDS_MAX_SLOTS

.. autodata:: TABLE_MAX_SIZE

.. autodata:: SORTED_MAX_SIZE

.. autodata:: SORTEDLR_MAX_SIZE

.. autodata:: SORTEDLR_MAX_SLOTS


Parameters for general cache behaviour
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. warning::

    The next parameters will not take any effect if passed to the open_file()
    function, so they can only be changed in a *global* way.  You can change
    them in the file, but this is strongly discouraged unless you know well
    what you are doing.

.. autodata:: DISABLE_EVERY_CYCLES

.. autodata:: ENABLE_EVERY_CYCLES

.. autodata:: LOWEST_HIT_RATIO


Parameters for the I/O buffer in Leaf objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autodata:: IO_BUFFER_SIZE

.. autodata:: BUFFER_TIMES


Miscellaneous
~~~~~~~~~~~~~

.. autodata:: EXPECTED_ROWS_EARRAY

.. autodata:: EXPECTED_ROWS_TABLE

.. autodata:: PYTABLES_SYS_ATTRS

.. autodata:: MAX_NUMEXPR_THREADS

.. autodata:: MAX_BLOSC_THREADS

.. autodata:: USER_BLOCK_SIZE

.. autodata:: ALLOW_PADDING


HDF5 driver management
~~~~~~~~~~~~~~~~~~~~~~

.. autodata:: DRIVER

.. autodata:: DRIVER_DIRECT_ALIGNMENT

.. autodata:: DRIVER_DIRECT_BLOCK_SIZE

.. autodata:: DRIVER_DIRECT_CBUF_SIZE

.. autodata:: DRIVER_CORE_INCREMENT

.. autodata:: DRIVER_CORE_BACKING_STORE

.. autodata:: DRIVER_CORE_IMAGE

.. autodata:: DRIVER_SPLIT_META_EXT

.. autodata:: DRIVER_SPLIT_RAW_EXT
Introduction
============
.. epigraph::

    La sabiduría no vale la pena si no es posible servirse de ella para
    inventar una nueva manera de preparar los garbanzos.

    [Wisdom isn't worth anything if you can't use it to come up with a new
    way to cook garbanzos.]

    -- Gabriel García Márquez, A wise Catalan in *"Cien años de soledad"*


The goal of PyTables is to enable the end user to manipulate easily data
*tables* and *array* objects in a hierarchical structure. The foundation of
the underlying hierarchical data organization is the excellent HDF5 library
(see :ref:`[HDGF1] <HDFG1>`).

It should be noted that this package is not intended to serve as a complete
wrapper for the entire HDF5 API, but only to provide a flexible, *very
pythonic* tool to deal with (arbitrarily) large amounts of data (typically
bigger than available memory) in tables and arrays organized in a
hierarchical and persistent disk storage structure.

A table is defined as a collection of records whose values are stored in
*fixed-length* fields. All records have the same structure and all values in
each field have the same *data type*. The terms *fixed-length* and strict
*data types* may seem to be a strange requirement for an interpreted language
like Python, but they serve a useful function if the goal is to save very
large quantities of data (such as is generated by many data acquisition
systems, Internet services or scientific applications, for example) in an
efficient manner that reduces demand on CPU time and I/O.

In order to emulate in Python records mapped to HDF5 C structs PyTables
implements a special class so as to easily define all its fields and other
properties. PyTables also provides a powerful interface to mine data in
tables. Records in tables are also known in the HDF5 naming scheme as
*compound* data types.

For example, you can define arbitrary tables in Python simply by declaring a
class with named fields and type information, such as in the following
example::

    class Particle(IsDescription):
        name      = StringCol(16)   # 16-character String
        idnumber  = Int64Col()      # signed 64-bit integer
        ADCcount  = UInt16Col()     # unsigned short integer
        TDCcount  = UInt8Col()      # unsigned byte
        grid_i    = Int32Col()      # integer
        grid_j    = Int32Col()      # integer

        # A sub-structure (nested data-type)
        class Properties(IsDescription):
            pressure = Float32Col(shape=(2,3)) # 2-D float array (single-precision)
            energy   = Float64Col(shape=(2,3,4)) # 3-D float array (double-precision)

You then pass this class to the table constructor, fill its rows with your
values, and save (arbitrarily large) collections of them to a file for
persistent storage. After that, the data can be retrieved and post-processed
quite easily with PyTables or even with another HDF5 application (in C,
Fortran, Java or whatever language that provides a library to interface with
HDF5).

Other important entities in PyTables are *array* objects, which are analogous
to tables with the difference that all of their components are homogeneous.
They come in different flavors, like *generic* (they provide a quick and fast
way to deal with for numerical arrays), *enlargeable* (arrays can be extended
along a single dimension) and *variable length* (each row in the array can
have a different number of elements).

The next section describes the most interesting capabilities of PyTables.


Main Features
-------------
PyTables takes advantage of the object orientation and introspection
capabilities offered by Python, the powerful data management features of
HDF5, and NumPy's flexibility and Numexpr's high-performance manipulation of
large sets of objects organized in a grid-like fashion to provide these
features:

- *Support for table entities:* You can tailor your data adding or deleting
  records in your tables. Large numbers of rows (up to 2**63, much more than
  will fit into memory) are supported as well.

- *Multidimensional and nested table cells:* You can declare a column to
  consist of values having any number of dimensions besides scalars, which is
  the only dimensionality allowed by the majority of relational databases.
  You can even declare columns that are made of other columns (of different
  types).

- *Indexing support for columns of tables:*
  Very useful if you have large tables and you want to quickly look up for
  values in columns satisfying some criteria.

- *Support for numerical arrays:*
  NumPy (see :ref:`[NUMPY] <NUMPY>`) arrays can be used as a useful
  complement of tables to store homogeneous data.

- *Enlargeable arrays:* You can add new
  elements to existing arrays on disk in any dimension you want (but only
  one). Besides, you are able to access just a slice of your datasets by
  using the powerful extended slicing mechanism, without need to load all
  your complete dataset in memory.

- *Variable length arrays:* The number of elements in these arrays can vary
  from row to row. This provides a lot of flexibility when dealing with
  complex data.

- *Supports a hierarchical data model:*
  Allows the user to clearly structure all data. PyTables builds up an
  *object tree* in memory that replicates the underlying file data structure.
  Access to objects in the file is achieved by walking through and
  manipulating this object tree.
  Besides, this object tree is built in a lazy way, for efficiency purposes.

- *User defined metadata:* Besides
  supporting system metadata (like the number of rows of a table, shape,
  flavor, etc.) the user may specify arbitrary metadata (as for example, room
  temperature, or protocol for IP traffic that was collected) that complement
  the meaning of actual data.

- *Ability to read/modify generic HDF5 files:* PyTables can access a wide
  range of objects in generic HDF5 files, like compound type datasets (that
  can be mapped to Table objects), homogeneous datasets (that can be mapped
  to Array objects) or variable length record datasets (that can be mapped to
  VLArray objects). Besides, if a dataset is not supported, it will be mapped
  to a special UnImplemented class (see :ref:`UnImplementedClassDescr`), that
  will let the user see that the data is there, although it will be
  unreachable (still, you will be able to access the attributes and some
  metadata in the dataset). With that, PyTables probably can access and
  *modify* most of the HDF5 files out there.

- *Data compression:* Supports data compression (using the *Zlib*, *LZO*,
  *bzip2* and *Blosc* compression libraries) out of the box. This is
  important when you have repetitive data patterns and don't want to spend
  time searching for an optimized way to store them (saving you time spent
  analyzing your data organization).

- *High performance I/O:* On modern systems storing large amounts of data,
  tables and array objects can be read and written at a speed only limited by
  the performance of the underlying I/O subsystem. Moreover, if your data is
  compressible, even that limit is surmountable!

- *Support of files bigger than 2 GB:*
  PyTables automatically inherits this capability from the underlying HDF5
  library (assuming your platform supports the C long long integer, or, on
  Windows, __int64).

- *Architecture-independent:* PyTables has been carefully coded (as HDF5
  itself) with little-endian/big-endian byte ordering issues in mind. So, you
  can write a file on a big-endian machine (like a Sparc or MIPS) and read it
  on other little-endian machine (like an Intel or Alpha) without problems.
  In addition, it has been tested successfully with 64 bit platforms
  (Intel-64, AMD-64, PowerPC-G5, MIPS, UltraSparc) using code generated with
  64 bit aware compilers.


.. _ObjectTreeSection:

The Object Tree
---------------
The hierarchical model of the underlying HDF5 library allows PyTables to
manage tables and arrays in a tree-like structure. In order to achieve this,
an *object tree* entity is *dynamically* created imitating the HDF5 structure
on disk. The HDF5 objects are read by walking through this object tree. You
can get a good picture of what kind of data is kept in the object by
examining the *metadata* nodes.

The different nodes in the object tree are instances of PyTables classes.
There are several types of classes, but the most important ones are the Node,
Group and Leaf classes. All nodes in a PyTables tree are instances of the
Node class. The Group and Leaf classes are descendants of Node. Group
instances (referred to as *groups* from now on) are a grouping structure
containing instances of zero or more groups or leaves, together with
supplementary metadata. Leaf instances (referred to as *leaves*) are
containers for actual data and can not contain further groups or leaves. The
Table, Array, CArray, EArray, VLArray and UnImplemented classes are
descendants of Leaf, and inherit all its properties.

Working with groups and leaves is similar in many ways to working with
directories and files on a Unix filesystem, i.e. a node (file or directory)
is always a *child* of one and only one group (directory), its *parent group*
[1]_.
Inside of that group, the node is accessed by its *name*. As is the case with
Unix directories and files, objects in the object tree are often referenced
by giving their full (absolute) path names. In PyTables this full path can be
specified either as string (such as '/subgroup2/table3', using / as a
parent/child separator) or as a complete object path written in a format
known as the *natural name* schema (such as file.root.subgroup2.table3).

Support for *natural naming* is a key aspect of PyTables. It means that the
names of instance variables of the node objects are the same as the names of
its children [2]_. This is very *Pythonic* and intuitive in many cases. Check
the tutorial :ref:`readingAndSelectingUsage` for usage examples.

You should also be aware that not all the data present in a file is loaded
into the object tree. The *metadata* (i.e. special data that describes the
structure of the actual data) is loaded only when the user want to access to
it (see later). Moreover, the actual data is not read until she request it
(by calling a method on a particular node). Using the object tree (the
metadata) you can retrieve information about the objects on disk such as
table names, titles, column names, data types in columns, numbers of rows,
or, in the case of arrays, their shapes, typecodes, etc. You can also search
through the tree for specific kinds of data then read it and process it. In a
certain sense, you can think of PyTables as a tool that applies the same
introspection capabilities of Python objects to large amounts of data in
persistent storage.

It is worth noting that PyTables sports a *metadata cache system* that loads
nodes *lazily* (i.e. on-demand), and unloads nodes that have not been used
for some time (following a *Least Recently Used* schema). It is important to
stress out that the nodes enter the cache after they have been unreferenced
(in the sense of Python reference counting), and that they can be revived (by
referencing them again) directly from the cache without performing the
de-serialization process from disk. This feature allows dealing with files
with large hierarchies very quickly and with low memory consumption, while
retaining all the powerful browsing capabilities of the previous
implementation of the object tree. See :ref:`[OPTIM] <OPTIM>` for more facts
about the advantages introduced by this new metadata cache system.

To better understand the dynamic nature of this object tree entity, let's
start with a sample PyTables script (which you can find in
examples/objecttree.py) to create an HDF5 file::

    import tables as tb

    class Particle(tb.IsDescription):
        identity = tb.StringCol(itemsize=22, dflt=" ", pos=0)  # character String
        idnumber = tb.Int16Col(dflt=1, pos = 1)  # short integer
        speed    = tb.Float32Col(dflt=1, pos = 2)  # single-precision

    # Open a file in "w"rite mode
    fileh = tb.open_file("objecttree.h5", mode = "w")

    # Get the HDF5 root group
    root = fileh.root

    # Create the groups
    group1 = fileh.create_group(root, "group1")
    group2 = fileh.create_group(root, "group2")

    # Now, create an array in root group
    array1 = fileh.create_array(root, "array1", ["string", "array"], "String array")

    # Create 2 new tables in group1
    table1 = fileh.create_table(group1, "table1", Particle)
    table2 = fileh.create_table("/group2", "table2", Particle)

    # Create the last table in group2
    array2 = fileh.create_array("/group1", "array2", [1,2,3,4])

    # Now, fill the tables
    for table in (table1, table2):
        # Get the record object associated with the table:
        row = table.row

        # Fill the table with 10 records
        for i in range(10):
            # First, assign the values to the Particle record
            row['identity']  = f'This is particle: {i:2d}'
            row['idnumber'] = i
            row['speed']  = i * 2.

            # This injects the Record values
            row.append()

        # Flush the table buffers
        table.flush()

    # Finally, close the file (this also will flush all the remaining buffers!)
    fileh.close()

This small program creates a simple HDF5 file called objecttree.h5 with the
structure that appears in :ref:`Figure 1 <objecttree-h5>` [3]_.
When the file is created, the metadata in the object tree is updated in
memory while the actual data is saved to disk. When you close the file the
object tree is no longer available. However, when you reopen this file the
object tree will be reconstructed in memory from the metadata on disk (this
is done in a lazy way, in order to load only the objects that are required by
the user), allowing you to work with it in exactly the same way as when you
originally created it.

.. _objecttree-h5:

.. figure:: images/objecttree-h5.png
    :align: center

    **Figure 1: An HDF5 example with 2 subgroups, 2 tables and 1 array.**

In :ref:`Figure2 <objecttree>`, you can see an example of the object tree
created when the above objecttree.h5 file is read (in fact, such an object
tree is always created when reading any supported generic HDF5 file).
It is worthwhile to take your time to understand it [4]_.
It will help you understand the relationships of in-memory PyTables objects.

.. _objecttree:

.. figure:: images/objecttree.png
    :width: 100%
    :align: center

    **Figure 2: A PyTables object tree example.**

---------------------------

.. [1] PyTables does not support hard links - for the moment.

.. [2] I got this simple but powerful idea from the excellent Objectify
       module by David Mertz (see :ref:`[MERTZ] <MERTZ>`).

.. [3] We have used ViTables (see :ref:`[VITABLES] <VITABLES>`) in order to
       create this snapshot.

.. [4] Bear in mind, however, that this diagram is *not* a standard UML class
       diagram; it is rather meant to show the connections between the
       PyTables objects and some of its most important attributes and
       methods.

Optimization tips
=================
.. epigraph::

    ... durch planmässiges Tattonieren.

    [... through systematic, palpable experimentation.]

    -- Johann Karl Friedrich Gauss [asked how he came upon his theorems]

.. currentmodule:: tables

On this chapter, you will get deeper knowledge of PyTables internals.
PyTables has many tunable features so that you can improve the performance of
your application.  If you are planning to deal with really large data, you
should read carefully this section in order to learn how to get an important
efficiency boost for your code.  But if your datasets are small (say, up to
10 MB) or your number of nodes is contained (up to 1000), you should not
worry about that as the default parameters in PyTables are already tuned for
those sizes (although you may want to adjust them further anyway).  At any
rate, reading this chapter will help you in your life with PyTables.


Understanding chunking
----------------------
The underlying HDF5 library that is used by PyTables allows for certain
datasets (the so-called *chunked* datasets) to take the data in bunches of a
certain length, named *chunks*, and write them on disk as a whole, i.e. the
HDF5 library treats chunks as atomic objects and disk I/O is always made in
terms of complete chunks.  This allows data filters to be defined by the
application to perform tasks such as compression, encryption, check-summing,
etc. on entire chunks.

HDF5 keeps a B-tree in memory that is used to map chunk structures on disk.
The more chunks that are allocated for a dataset the larger the B-tree.
Large B-trees take memory and cause file storage overhead as well as more
disk I/O and higher contention for the metadata cache.  Consequently, it's
important to balance between memory and I/O overhead (small B-trees) and time
to access data (big B-trees).

In the next couple of sections, you will discover how to inform PyTables
about the expected size of your datasets for allowing a sensible computation
of the chunk sizes.  Also, you will be presented some experiments so that you
can get a feeling on the consequences of manually specifying the chunk size.
Although doing this latter is only reserved to experienced people, these
benchmarks may allow you to understand more deeply the chunk size
implications and let you quickly start with the fine-tuning of this important
parameter.


.. _expectedRowsOptim:

Informing PyTables about expected number of rows in tables or arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PyTables can determine a sensible chunk size to your dataset size if you
help it by providing an estimation of the final number of rows for an
extensible leaf [1]_.  You should provide this information at leaf creation
time by passing this value to the ``expectedrows`` argument of the
:meth:`File.create_table` method or :meth:`File.create_earray` method (see
:ref:`EArrayClassDescr`).

When your leaf size is bigger than 10 MB (take this figure only as a
reference, not strictly), by providing this guess you will be optimizing the
access to your data.  When the table or array size is larger than, say 100MB,
you are *strongly* suggested to provide such a guess; failing to do that may
cause your application to do very slow I/O operations and to demand *huge*
amounts of memory. You have been warned!


.. _chunksizeFineTune:

Fine-tuning the chunksize
~~~~~~~~~~~~~~~~~~~~~~~~~
.. warning::

    This section is mostly meant for experts.  If you are a beginner, you
    must know that setting manually the chunksize is a potentially dangerous
    action.

Most of the time, informing PyTables about the extent of your dataset is
enough.  However, for more sophisticated applications, when one has special
requirements for doing the I/O or when dealing with really large datasets,
you should really understand the implications of the chunk size in order to
be able to find the best value for your own application.

You can specify the chunksize for every chunked dataset in PyTables by
passing the chunkshape argument to the corresponding constructors. It is
important to point out that chunkshape is not exactly the same thing than a
chunksize; in fact, the chunksize of a dataset can be computed multiplying
all the dimensions of the chunkshape among them and multiplying the outcome
by the size of the atom.

We are going to describe a series of experiments where an EArray of 15 GB is
written with different chunksizes, and then it is accessed in both sequential
(i.e. first element 0, then element 1 and so on and so forth until the data
is exhausted) and random mode (i.e. single elements are read randomly all
through the dataset). These benchmarks have been carried out with
PyTables 2.1 on a machine with an Intel Core2 processor @ 3 GHz and a RAID-0
made of two SATA disks spinning at 7200 RPM, and using GNU/Linux with an XFS
filesystem.  The script used for the benchmarks is available in
bench/optimal-chunksize.py.

In figures :ref:`Figure 1 <createTime-chunksize>`,
:ref:`Figure 2 <fileSizes-chunksize>`, :ref:`Figure 3 <seqTime-chunksize>`
and :ref:`Figure 4 <randomTime-chunksize>`, you can see how the chunksize
affects different aspects, like creation time, file sizes, sequential read
time and random read time.  So, if you properly inform PyTables about the
extent of your datasets, you will get an automatic chunksize value (256 KB in
this case) that is pretty optimal for most of uses.  However, if what you
want is, for example, optimize the creation time when using the
Zlib compressor, you may want to reduce the chunksize to 32 KB (see
:ref:`Figure 1 <createTime-chunksize>`). Or, if your goal is to optimize the
sequential access time for an dataset compressed with Blosc, you may want to
increase the chunksize to 512 KB (see :ref:`Figure 3 <seqTime-chunksize>`).

You will notice that, by manually specifying the chunksize of a leave you
will not normally get a drastic increase in performance, but at least, you
have the opportunity to fine-tune such an important parameter for improve
performance.

.. _createTime-chunksize:

.. figure:: images/create-chunksize-15GB.png
    :align: center

    **Figure 1. Creation time per element for a 15 GB EArray and different
    chunksizes.**


.. _fileSizes-chunksize:

.. figure:: images/filesizes-chunksize-15GB.png
    :align: center

    **Figure 2. File sizes for a 15 GB EArray and different chunksizes.**

.. _seqTime-chunksize:


.. figure:: images/seq-chunksize-15GB.png
    :align: center

    **Figure 3. Sequential access time per element for a 15 GB EArray and
    different chunksizes.**


.. _randomTime-chunksize:

.. figure:: images/random-chunksize-15GB.png
    :align: center

    **Figure 4. Random access time per element for a 15 GB EArray and
    different chunksizes.**


Finally, it is worth noting that adjusting the chunksize can be specially
important if you want to access your dataset by blocks of certain dimensions.
In this case, it is normally a good idea to set your chunkshape to be the
same than these dimensions; you only have to be careful to not end with a too
small or too large chunksize.  As always, experimenting prior to pass your
application into production is your best ally.


.. _searchOptim:

Accelerating your searches
--------------------------

.. note::

    Many of the explanations and plots in this section and the forthcoming
    ones still need to be updated to include Blosc (see
    :ref:`[BLOSC] <BLOSC>`), the new and powerful compressor added in
    PyTables 2.2 series.  You should expect it to be the fastest compressor
    among all the described here, and its use is strongly recommended
    whenever you need extreme speed and not a very high compression ratio.

Searching in tables is one of the most common and time consuming operations
that a typical user faces in the process of mining through his data.  Being
able to perform queries as fast as possible will allow more opportunities for
finding the desired information quicker and also allows to deal with larger
datasets.

PyTables offers many sort of techniques so as to speed-up the search process
as much as possible and, in order to give you hints to use them based, a
series of benchmarks have been designed and carried out.  All the results
presented in this section have been obtained with synthetic, random data and
using PyTables 2.1.  Also, the tests have been conducted on a machine with an
Intel Core2 (64-bit) @ 3 GHz processor with RAID-0 disk storage (made of four
spinning disks @ 7200 RPM), using GNU/Linux with an XFS filesystem.  The
script used for the benchmarks is available in bench/indexed_search.py.
As your data, queries and platform may be totally different for your case,
take this just as a guide because your mileage may vary (and will vary).

In order to be able to play with tables with a number of rows as large as
possible, the record size has been chosen to be rather small (24 bytes). Here
it is its definition::

    class Record(tables.IsDescription):
        col1 = tables.Int32Col()
        col2 = tables.Int32Col()
        col3 = tables.Float64Col()
        col4 = tables.Float64Col()

In the next sections, we will be optimizing the times for a relatively
complex query like this::

    result = [row['col2'] for row in table if (
              ((row['col4'] >= lim1 and row['col4'] < lim2) or
              ((row['col2'] > lim3 and row['col2'] < lim4])) and
              ((row['col1']+3.1*row['col2']+row['col3']*row['col4']) > lim5)
              )]

(for future reference, we will call this sort of queries *regular* queries).
So, if you want to see how to greatly improve the time taken to run queries
like this, keep reading.


.. _inkernelSearch:

In-kernel searches
~~~~~~~~~~~~~~~~~~
PyTables provides a way to accelerate data selections inside of a single
table, through the use of the :ref:`TableMethods_querying` iterator and
related query methods. This mode of selecting data is called *in-kernel*.
Let's see an example of an *in-kernel* query based on the *regular* one
mentioned above::

    result = [row['col2'] for row in table.where(
                '''(((col4 >= lim1) & (col4 < lim2)) |
                   ((col2 > lim3) & (col2 < lim4)) &
                   ((col1+3.1*col2+col3*col4) > lim5))''')]

This simple change of mode selection can improve search times quite a lot and
actually make PyTables very competitive when compared against typical
relational databases as you can see in :ref:`Figure 5 <sequentialTimes-10m>`
and :ref:`Figure 6 <sequentialTimes-1g>`.

.. _sequentialTimes-10m:

.. figure:: images/Q7-10m-noidx.png
    :align: center

    **Figure 5. Times for non-indexed complex queries in a small table with
    10 millions of rows: the data fits in memory.**

By looking at :ref:`Figure 5 <sequentialTimes-10m>` you can see how in the
case that table data fits easily in memory, in-kernel searches on
uncompressed tables are generally much faster (10x) than standard queries as
well as PostgreSQL (5x).  Regarding compression, we can see how Zlib
compressor actually slows down the performance of in-kernel queries by a
factor 3.5x; however, it remains faster than PostgreSQL (40%).
On his hand, LZO compressor only decreases the performance by a 75% with
respect to uncompressed in-kernel queries and is still a lot faster than
PostgreSQL (3x).  Finally, one can observe that, for low selectivity queries
(large number of hits), PostgreSQL performance degrades quite steadily, while
in PyTables this slow down rate is significantly smaller.  The reason of this
behaviour is not entirely clear to the authors, but the fact is clearly
reproducible in our benchmarks.

But, why in-kernel queries are so fast when compared with regular ones?.
The answer is that in regular selection mode the data for all the rows in
table has to be brought into Python space so as to evaluate the condition and
decide if the corresponding field should be added to the result list.  On the
contrary, in the in-kernel mode, the condition is passed to the PyTables
kernel (hence the name), written in C, and evaluated there at full C speed
(with the help of the integrated Numexpr package, see
:ref:`[NUMEXPR] <NUMEXPR>`), so that the only values that are brought to
Python space are the rows that fulfilled the condition.  Hence, for
selections that only have a relatively small number of hits (compared with
the total amount of rows), the savings are very large.  It is also
interesting to note the fact that, although for queries with a large number
of hits the speed-up is not as high, it is still very important.

On the other hand, when the table is too large to fit in memory (see
:ref:`Figure 6 <sequentialTimes-1g>`), the difference in speed between
regular and in-kernel is not so important, but still significant (2x).  Also,
and curiously enough, large tables compressed with Zlib offers slightly
better performance (around 20%) than uncompressed ones; this is because the
additional CPU spent by the uncompressor is compensated by the savings in
terms of net I/O (one has to read less actual data from disk).  However, when
using the extremely fast LZO compressor, it gives a clear advantage over
Zlib, and is up to 2.5x faster than not using compression at all.  The reason
is that LZO decompression speed is much faster than Zlib, and that allows
PyTables to read the data at full disk speed (i.e. the bottleneck is in the
I/O subsystem, not in the CPU).  In this case the compression rate is around
2.5x, and this is why the data can be read 2.5x faster.  So, in general,
using the LZO compressor is the best way to ensure best reading/querying
performance for out-of-core datasets (more about how compression affects
performance in :ref:`compressionIssues`).

.. _sequentialTimes-1g:

.. figure:: images/Q8-1g-noidx.png
    :align: center

    **Figure 6. Times for non-indexed complex queries in a large table with 1
    billion of rows: the data does not fit in memory.**

Furthermore, you can mix the *in-kernel* and *regular* selection modes for
evaluating arbitrarily complex conditions making use of external functions.
Look at this example::

    result = [ row['var2']
               for row in table.where('(var3 == "foo") & (var1 <= 20)')
               if your_function(row['var2']) ]

Here, we use an *in-kernel* selection to choose rows according to the values
of the var3 and var1 fields.  Then, we apply a *regular* selection to
complete the query. Of course, when you mix the *in-kernel* and *regular*
selection modes you should pass the most restrictive condition to the
*in-kernel* part, i.e. to the where() iterator.  In situations where it is
not clear which is the most restrictive condition, you might want to
experiment a bit in order to find the best combination.

However, since in-kernel condition strings allow rich expressions allowing
the coexistence of multiple columns, variables, arithmetic operations and
many typical functions, it is unlikely that you will be forced to use
external regular selections in conditions of small to medium complexity.
See :ref:`condition_syntax` for more information on in-kernel condition
syntax.


Indexed searches
~~~~~~~~~~~~~~~~
When you need more speed than *in-kernel* selections can offer you, PyTables
offers a third selection method, the so-called *indexed* mode (based on the
highly efficient OPSI indexing engine ).  In this mode, you have to decide
which column(s) you are going to apply your selections over, and index them.
Indexing is just a kind of sorting operation over a column, so that searches
along such a column (or columns) will look at this sorted information by
using a *binary search* which is much faster than the *sequential search*
described in the previous section.

You can index the columns you want by calling the :meth:`Column.create_index`
method on an already created table.  For example::

    indexrows = table.cols.var1.create_index()
    indexrows = table.cols.var2.create_index()
    indexrows = table.cols.var3.create_index()

will create indexes for all var1, var2 and var3 columns.

After you have indexed a series of columns, the PyTables query optimizer will
try hard to discover the usable indexes in a potentially complex expression.
However, there are still places where it cannot determine that an index can
be used. See below for examples where the optimizer can safely determine if
an index, or series of indexes, can be used or not.

Example conditions where an index can be used:

- var1 >= "foo" (var1 is used)

- var1 >= mystr (var1 is used)

- (var1 >= "foo") & (var4 > 0.0) (var1 is used)

- ("bar" <= var1) & (var1 < "foo") (var1 is used)

- (("bar" <= var1) & (var1 < "foo")) & (var4 > 0.0) (var1 is used)

- (var1 >= "foo") & (var3 > 10) (var1 and var3 are used)

- (var1 >= "foo") | (var3 > 10) (var1 and var3 are used)

- ~(var1 >= "foo") | ~(var3 > 10) (var1 and var3 are used)

Example conditions where an index can *not* be used:

- var4 > 0.0 (var4 is not indexed)

- var1 != 0.0 (range has two pieces)

- ~(("bar" <= var1) & (var1 < "foo")) & (var4 > 0.0) (negation of a complex boolean expression)

.. note:: From PyTables 2.3 on, several indexes can be used in a single query.

.. note::

    If you want to know for sure whether a particular query will use indexing
    or not (without actually running it), you are advised to use the
    :meth:`Table.will_query_use_indexing` method.

One important aspect of the new indexing in PyTables (>= 2.3) is that it has
been designed from the ground up with the goal of being capable to
effectively manage very large tables.  To this goal, it sports a wide
spectrum of different quality levels (also called optimization levels) for
its indexes so that the user can choose the best one that suits her needs
(more or less size, more or less performance).

In :ref:`Figure 7 <createIndexTimes>`, you can see that the times to index
columns in tables can be really short.  In particular, the time to index a
column with 1 billion rows (1 Gigarow) with the lowest optimization level is
less than 4 minutes while indexing the same column with full optimization (so
as to get a completely sorted index or CSI) requires around 1 hour.  These
are rather competitive figures compared with a relational database (in this
case, PostgreSQL 8.3.1, which takes around 1.5 hours for getting the index
done).  This is because PyTables is geared towards read-only or append-only
tables and takes advantage of this fact to optimize the indexes properly.  On
the contrary, most relational databases have to deliver decent performance in
other scenarios as well (specially updates and deletions), and this fact
leads not only to slower index creation times, but also to indexes taking
much more space on disk, as you can see in :ref:`Figure 8 <indexSizes>`.

.. _createIndexTimes:

.. figure:: images/create-index-time-int32-float64.png
    :align: center

    **Figure 7. Times for indexing an Int32 and Float64 column.**


.. _indexSizes:

.. figure:: images/indexes-sizes2.png
    :align: center

    **Figure 8. Sizes for an index of a Float64 column with 1 billion of rows.**


The user can select the index quality by passing the desired optlevel and
kind arguments to the :meth:`Column.create_index` method.  We can see in
figures :ref:`Figure 7 <createIndexTimes>` and :ref:`Figure 8 <indexSizes>`
how the different optimization levels affects index time creation and index
sizes.

So, which is the effect of the different optimization levels in terms of
query times?  You can see that in :ref:`Figure 9 <queryTimes-indexed-optlevels>`.

.. _queryTimes-indexed-optlevels:

.. figure:: images/Q8-1g-idx-optlevels.png
    :align: center

    **Figure 9. Times for complex queries with a cold cache (mean of 5 first
    random queries) for different optimization levels. Benchmark made on a machine with Intel Core2 (64-bit) @ 3 GHz processor with RAID-0 disk storage.**

Of course, compression also has an effect when doing indexed queries,
although not very noticeable, as can be seen in
:ref:`Figure 10 <queryTimes-indexed-compress>`.
As you can see, the difference between using no compression and using Zlib or
LZO is very little, although LZO achieves relatively better performance
generally speaking.

.. _queryTimes-indexed-compress:

.. figure:: images/Q8-1g-idx-compress.png
    :align: center

    **Figure 10. Times for complex queries with a cold cache (mean of 5 first
    random queries) for different compressors.**

You can find a more complete description and benchmarks about OPSI, the
indexing system of PyTables (>= 2.3) in :ref:`[OPSI] <OPSI>`.


Indexing and Solid State Disks (SSD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Lately, the long promised Solid State Disks (SSD for brevity) with decent
capacities and affordable prices have finally hit the market and will
probably stay in coexistence with the traditional spinning disks for the
foreseeable future (separately or forming *hybrid* systems).  SSD have many
advantages over spinning disks, like much less power consumption and better
throughput.  But of paramount importance, specially in the context of
accelerating indexed queries, is its very reduced latency during disk seeks,
which is typically 100x better than traditional disks.
Such a huge improvement has to have a clear impact in reducing the query
times, specially when the selectivity is high (i.e. the number of hits is
small).

In order to offer an estimate on the performance improvement we can expect
when using a low-latency SSD instead of traditional spinning disks, the
benchmark in the previous section has been repeated, but this time using a
single SSD disk instead of the four spinning disks in RAID-0.  The result can
be seen in :ref:`Figure 11 <queryTimes-indexed-SSD>`.  There one can see how
a query in a table of 1 billion of rows with 100 hits took just 1 tenth of
second when using a SSD, instead of 1 second that needed the RAID made of
spinning disks.  This factor of 10x of speed-up for high-selectivity queries
is nothing to sneeze at, and should be kept in mind when really high
performance in queries is needed.  It is also interesting that using
compression with LZO does have a clear advantage over when no compression is
done.

.. _queryTimes-indexed-SSD:

.. figure:: images/Q8-1g-idx-SSD.png
    :align: center

    **Figure 11. Times for complex queries with a cold cache (mean of 5 first
    random queries) for different disk storage (SSD vs spinning disks).**

Finally, we should remark that SSD can't compete with traditional spinning
disks in terms of capacity as they can only provide, for a similar cost,
between 1/10th and 1/50th of the size of traditional disks.  It is here where
the compression capabilities of PyTables can be very helpful because both
tables and indexes can be compressed and the final space can be reduced by
typically 2x to 5x (4x to 10x when compared with traditional relational
databases).
Best of all, as already mentioned, performance is not degraded when
compression is used, but actually *improved*.
So, by using PyTables and SSD you can query larger datasets that otherwise
would require spinning disks when using other databases

In fact, we were unable to run the PostgreSQL benchmark in this case because
the space needed exceeded the capacity of our SSD., while allowing
improvements in the speed of indexed queries between 2x (for medium to low
selectivity queries) and 10x (for high selectivity queries).


Achieving ultimate speed: sorted tables and beyond
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

    Sorting a large table is a costly operation.  The next procedure should
    only be performed when your dataset is mainly read-only and meant to be
    queried many times.

When querying large tables, most of the query time is spent in locating the
interesting rows to be read from disk.  In some occasions, you may have
queries whose result depends *mainly* of one single column (a query with only
one single condition is the trivial example), so we can guess that sorting
the table by this column would lead to locate the interesting rows in a much
more efficient way (because they would be mostly *contiguous*).  We are going
to confirm this guess.

For the case of the query that we have been using in the previous sections::

    result = [row['col2'] for row in table.where(
                '''(((col4 >= lim1) & (col4 < lim2)) |
                   ((col2 > lim3) & (col2 < lim4)) &
                   ((col1+3.1*col2+col3*col4) > lim5))''')]

it is possible to determine, by analysing the data distribution and the query
limits, that col4 is such a *main column*.  So, by ordering the table by the
col4 column (for example, by specifying setting the column to sort by in the
sortby parameter in the :meth:`Table.copy` method and re-indexing col2 and
col4 afterwards, we should get much faster performance for our query.  This
is effectively demonstrated in :ref:`Figure 12 <queryTimes-indexed-sorted>`,
where one can see how queries with a low to medium (up to 10000) number of
hits can be done in around 1 tenth of second for a RAID-0 setup and in around
1 hundredth of second for a SSD disk.  This represents up to more that 100x
improvement in speed with respect to the times with unsorted tables.  On the
other hand, when the number of hits is large (> 1 million), the query times
grow almost linearly, showing a near-perfect scalability for both RAID-0 and
SSD setups (the sequential access to disk becomes the bottleneck in this
case).

.. _queryTimes-indexed-sorted:

.. figure:: images/Q8-1g-idx-sorted.png
    :align: center

    **Figure 12. Times for complex queries with a cold cache (mean of 5 first
    random queries) for unsorted and sorted tables.**

Even though we have shown many ways to improve query times that should
fulfill the needs of most of people, for those needing more, you can for sure
discover new optimization opportunities.  For example, querying against
sorted tables is limited mainly by sequential access to data on disk and data
compression capability, so you may want to read :ref:`chunksizeFineTune`, for
ways on improving this aspect.
Reading the other sections of this chapter will help in finding new roads for
increasing the performance as well.  You know, the limit for stopping the
optimization process is basically your imagination (but, most plausibly, your
available time ;-).


.. _compressionIssues:

Compression issues
------------------
One of the beauties of PyTables is that it supports compression on tables and
arrays [2]_, although it is not used by default. Compression of big amounts
of data might be a bit controversial feature, because it has a legend of
being a very big consumer of CPU time resources. However, if you are willing
to check if compression can help not only by reducing your dataset file size
but *also* by improving I/O efficiency, specially when dealing with very
large datasets, keep reading.


A study on supported compression libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The compression library used by default is the *Zlib* (see
:ref:`[ZLIB] <ZLIB>`). Since HDF5 *requires* it, you can safely use it and
expect that your HDF5 files will be readable on any other platform that has
HDF5 libraries installed. Zlib provides good compression ratio, although
somewhat slow, and reasonably fast decompression.  Because of that, it is a
good candidate to be used for compressing you data.

However, in some situations it is critical to have a *very good decompression
speed* (at the expense of lower compression ratios or more CPU wasted on
compression, as we will see soon). In others, the emphasis is put in
achieving the *maximum compression ratios*, no matter which reading speed
will result. This is why support for two additional compressors has been
added to PyTables: LZO (see :ref:`[LZO] <LZO>`) and bzip2 (see
:ref:`[BZIP2] <BZIP2>`). Following the author of LZO (and checked by the
author of this section, as you will see soon), LZO offers pretty fast
compression and extremely fast decompression. In fact, LZO is so fast when
compressing/decompressing that it may well happen (that depends on your data,
of course) that writing or reading a compressed dataset is sometimes faster
than if it is not compressed at all (specially when dealing with extremely
large datasets). This fact is very important, specially if you have to deal
with very large amounts of data. Regarding bzip2, it has a reputation of
achieving excellent compression ratios, but at the price of spending much
more CPU time, which results in very low compression/decompression speeds.

Be aware that the LZO and bzip2 support in PyTables is not standard on HDF5,
so if you are going to use your PyTables files in other contexts different
from PyTables you will not be able to read them. Still, see the
:ref:`ptrepackDescr` (where the ptrepack utility is described) to find a way
to free your files from LZO or bzip2 dependencies, so that you can use these
compressors locally with the warranty that you can replace them with Zlib (or
even remove compression completely) if you want to use these files with other
HDF5 tools or platforms afterwards.

In order to allow you to grasp what amount of compression can be achieved,
and how this affects performance, a series of experiments has been carried
out. All the results presented in this section (and in the next one) have
been obtained with synthetic data and using PyTables 1.3. Also, the tests
have been conducted on a IBM OpenPower 720 (e-series) with a PowerPC G5 at
1.65 GHz and a hard disk spinning at 15K RPM. As your data and platform may
be totally different for your case, take this just as a guide because your
mileage may vary. Finally, and to be able to play with tables with a number
of rows as large as possible, the record size has been chosen to be small (16
bytes). Here is its definition::

    class Bench(IsDescription):
        var1 = StringCol(length=4)
        var2 = IntCol()
        var3 = FloatCol()

With this setup, you can look at the compression ratios that can be achieved
in :ref:`Figure 13 <comprTblComparison>`. As you can see, LZO is the
compressor that performs worse in this sense, but, curiously enough, there is
not much difference between Zlib and bzip2.

.. _comprTblComparison:

.. figure:: images/compressed-recordsize.png
    :align: center

    **Figure 13. Comparison between different compression libraries.**

Also, PyTables lets you select different compression levels for Zlib and
bzip2, although you may get a bit disappointed by the small improvement that
these compressors show when dealing with a combination of numbers and strings
as in our example. As a reference, see plot
:ref:`Figure 14 <comprZlibComparison>` for a comparison of the compression
achieved by selecting different levels of Zlib.  Very oddly, the best
compression ratio corresponds to level 1 (!).  See later for an explanation
and more figures on this subject.

.. _comprZlibComparison:

.. figure:: images/compressed-recordsize-zlib.png
    :align: center

    **Figure 14. Comparison between different compression levels of Zlib.**

Have also a look at :ref:`Figure 15 <comprWriteComparison>`. It shows how the
speed of writing rows evolves as the size (number of rows) of the table
grows. Even though in these graphs the size of one single row is 16 bytes,
you can most probably extrapolate these figures to other row sizes.

.. _comprWriteComparison:

.. figure:: images/compressed-writing.png
    :align: center

    **Figure 15. Writing tables with several compressors.**

In :ref:`Figure 16 <comprReadNoCacheComparison>` you can see how compression
affects the reading performance. In fact, what you see in the plot is an
*in-kernel selection* speed, but provided that this operation is very fast
(see :ref:`inkernelSearch`), we can accept it as an actual read test.
Compared with the reference line without compression, the general trend here
is that LZO does not affect too much the reading performance (and in some
points it is actually better), Zlib makes speed drop to a half, while bzip2
is performing very slow (up to 8x slower).

Also, in the same :ref:`Figure 16 <comprReadNoCacheComparison>` you can
notice some strange peaks in the speed that we might be tempted to attribute
to libraries on which PyTables relies (HDF5, compressors...), or to PyTables
itself.
However, :ref:`Figure 17 <comprReadCacheComparison>` reveals that, if we put
the file in the filesystem cache (by reading it several times before, for
example), the evolution of the performance is much smoother. So, the most
probable explanation would be that such peaks are a consequence of the
underlying OS filesystem, rather than a flaw in PyTables (or any other
library behind it). Another consequence that can be derived from the
aforementioned plot is that LZO decompression performance is much better than
Zlib, allowing an improvement in overall speed of more than 2x, and perhaps
more important, the read performance for really large datasets (i.e. when
they do not fit in the OS filesystem cache) can be actually *better* than not
using compression at all. Finally, one can see that reading performance is
very badly affected when bzip2 is used (it is 10x slower than LZO and 4x than
Zlib), but this was somewhat expected anyway.

.. _comprReadNoCacheComparison:

.. figure:: images/compressed-select-nocache.png
    :align: center

    **Figure 16. Selecting values in tables with several compressors.
    The file is not in the OS cache.**


.. _comprReadCacheComparison:

.. figure:: images/compressed-select-cache.png
    :align: center

    **Figure 17. Selecting values in tables with several compressors.
    The file is in the OS cache.**

So, generally speaking and looking at the experiments above, you can expect
that LZO will be the fastest in both compressing and decompressing, but the
one that achieves the worse compression ratio (although that may be just OK
for many situations, specially when used with shuffling - see
:ref:`ShufflingOptim`).  bzip2 is the slowest, by large, in both compressing
and decompressing, and besides, it does not achieve any better compression
ratio than Zlib. Zlib represents a balance between them: it's somewhat slow
compressing (2x) and decompressing (3x) than LZO, but it normally achieves
better compression ratios.

Finally, by looking at the plots :ref:`Figure 18 <comprWriteZlibComparison>`,
:ref:`Figure 19 <comprReadZlibComparison>`, and the aforementioned
:ref:`Figure 14 <comprZlibComparison>` you can see why the recommended
compression level to use for all compression libraries is 1.  This is the
lowest level of compression, but as the size of the underlying HDF5 chunk
size is normally rather small compared with the size of compression buffers,
there is not much point in increasing the latter (i.e. increasing the
compression level).  Nonetheless, in some situations (like for example, in
extremely large tables or arrays, where the computed chunk size can be rather
large) you may want to check, on your own, how the different compression
levels do actually affect your application.

You can select the compression library and level by setting the complib and
complevel keywords in the Filters class (see :ref:`FiltersClassDescr`). A
compression level of 0 will completely disable compression (the default), 1
is the less memory and CPU time demanding level, while 9 is the maximum level
and the most memory demanding and CPU intensive. Finally, have in mind that
LZO is not accepting a compression level right now, so, when using LZO, 0
means that compression is not active, and any other value means that LZO is
active.

So, in conclusion, if your ultimate goal is writing and reading as fast as
possible, choose LZO. If you want to reduce as much as possible your data,
while retaining acceptable read speed, choose Zlib. Finally, if portability
is important for you, Zlib is your best bet. So, when you want to use bzip2?
Well, looking at the results, it is difficult to recommend its use in
general, but you may want to experiment with it in those cases where you know
that it is well suited for your data pattern (for example, for dealing with
repetitive string datasets).

.. _comprWriteZlibComparison:

.. figure:: images/compressed-writing-zlib.png
    :align: center

    **Figure 18. Writing in tables with different levels of compression.**

.. _comprReadZlibComparison:

.. figure:: images/compressed-select-cache-zlib.png
    :align: center

    **Figure 19. Selecting values in tables with different levels of
    compression. The file is in the OS cache.**


.. _ShufflingOptim:

Shuffling (or how to make the compression process more effective)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The HDF5 library provides an interesting filter that can leverage the results
of your favorite compressor. Its name is *shuffle*, and because it can
greatly benefit compression and it does not take many CPU resources (see
below for a justification), it is active *by default* in PyTables whenever
compression is activated (independently of the chosen compressor). It is
deactivated when compression is off (which is the default, as you already
should know). Of course, you can deactivate it if you want, but this is not
recommended.

.. note::

   Since PyTables 3.3, a new *bitshuffle* filter for Blosc compressor
   has been added.  Contrarily to *shuffle* that shuffles bytes,
   *bitshuffle* shuffles the chunk data at bit level which **could**
   improve compression ratios at the expense of some speed penalty.
   Look at the :ref:`FiltersClassDescr` documentation on how to
   activate bitshuffle and experiment with it so as to decide if it
   can be useful for you.


So, how does this mysterious filter exactly work? From the HDF5 reference
manual::

    "The shuffle filter de-interlaces a block of data by reordering the
    bytes. All the bytes from one consistent byte position of each data
    element are placed together in one block; all bytes from a second
    consistent byte position of each data element are placed together a
    second block; etc. For example, given three data elements of a 4-byte
    datatype stored as 012301230123, shuffling will re-order data as
    000111222333. This can be a valuable step in an effective compression
    algorithm because the bytes in each byte position are often closely
    related to each other and putting them together can increase the
    compression ratio."

In :ref:`Figure 20 <comprShuffleComparison>` you can see a benchmark that
shows how the *shuffle* filter can help the different libraries in
compressing data. In this experiment, shuffle has made LZO compress almost 3x
more (!), while Zlib and bzip2 are seeing improvements of 2x. Once again, the
data for this experiment is synthetic, and *shuffle* seems to do a great work
with it, but in general, the results will vary in each case [3]_.

.. _comprShuffleComparison:

.. figure:: images/compressed-recordsize-shuffle.png
    :align: center

    **Figure 20. Comparison between different compression libraries with and
    without the shuffle filter.**

At any rate, the most remarkable fact about the *shuffle* filter is the
relatively high level of compression that compressor filters can achieve when
used in combination with it. A curious thing to note is that the Bzip2
compression rate does not seem very much improved (less than a 40%), and what
is more striking, Bzip2+shuffle does compress quite *less* than Zlib+shuffle
or LZO+shuffle combinations, which is kind of unexpected. The thing that
seems clear is that Bzip2 is not very good at compressing patterns that
result of shuffle application. As always, you may want to experiment with
your own data before widely applying the Bzip2+shuffle combination in order
to avoid surprises.

Now, how does shuffling affect performance? Well, if you look at plots
:ref:`Figure 21 <comprWriteShuffleComparison>`,
:ref:`Figure 22 <comprReadNoCacheShuffleComparison>` and
:ref:`Figure 23 <comprReadCacheShuffleComparison>`, you will get a somewhat
unexpected (but pleasant) surprise. Roughly, *shuffle* makes the writing
process (shuffling+compressing) faster (approximately a 15% for LZO, 30% for
Bzip2 and a 80% for Zlib), which is an interesting result by itself.
But perhaps more exciting is the fact that the reading process
(unshuffling+decompressing) is also accelerated by a similar extent (a 20%
for LZO, 60% for Zlib and a 75% for Bzip2, roughly).

.. _comprWriteShuffleComparison:

.. figure:: images/compressed-writing-shuffle.png
    :align: center

    **Figure 21. Writing with different compression libraries with and
    without the shuffle filter.**


.. _comprReadNoCacheShuffleComparison:

.. figure:: images/compressed-select-nocache-shuffle-only.png
    :align: center

    **Figure 22. Reading with different compression libraries with the
    shuffle filter. The file is not in OS cache.**



.. _comprReadCacheShuffleComparison:

.. figure:: images/compressed-select-cache-shuffle.png
    :align: center

    **Figure 23. Reading with different compression libraries with and
    without the shuffle filter. The file is in OS cache.**

You may wonder why introducing another filter in the write/read pipelines
does effectively accelerate the throughput. Well, maybe data elements are
more similar or related column-wise than row-wise, i.e. contiguous elements
in the same column are more alike, so shuffling makes the job of the
compressor easier (faster) and more effective (greater ratios). As a side
effect, compressed chunks do fit better in the CPU cache (at least, the
chunks are smaller!) so that the process of unshuffle/decompress can make a
better use of the cache (i.e. reducing the number of CPU cache faults).

So, given the potential gains (faster writing and reading, but specially
much improved compression level), it is a good thing to have such a filter
enabled by default in the battle for discovering redundancy when you want to
compress your data, just as PyTables does.


Using Psyco
-----------
Psyco (see :ref:`[PSYCO] <PSYCO>`) is a kind of specialized compiler for
Python that typically accelerates Python applications with no change in
source code. You can think of Psyco as a kind of just-in-time (JIT) compiler,
a little bit like Java's, that emits machine code on the fly instead of
interpreting your Python program step by step. The result is that your
unmodified Python programs run faster.

Psyco is very easy to install and use, so in most scenarios it is worth to
give it a try. However, it only runs on Intel 386 architectures, so if you
are using other architectures, you are out of luck (and, moreover, it seems
that there are no plans to support other platforms).  Besides, with the
addition of flexible (and very fast) in-kernel queries (by the way, they
cannot be optimized at all by Psyco), the use of Psyco will only help in
rather few scenarios.  In fact, the only important situation that you might
benefit right now from using Psyco (I mean, in PyTables contexts) is for
speeding-up the write speed in tables when using the Row interface (see
:ref:`RowClassDescr`).  But again, this latter case can also be accelerated
by using the :meth:`Table.append` method and building your own buffers [4]_.

As an example, imagine that you have a small script that reads and selects
data over a series of datasets, like this::

    def read_file(filename):
        "Select data from all the tables in filename"
        fileh = open_file(filename, mode = "r")
        result = []
        for table in fileh("/", 'Table'):
            result = [p['var3'] for p in table if p['var2'] <= 20]
        fileh.close()
        return result

    if __name__=="__main__":
        print(read_file("myfile.h5"))

In order to accelerate this piece of code, you can rewrite your main program
to look like::

    if __name__=="__main__":
        import psyco
        psyco.bind(read_file)
        print(read_file("myfile.h5"))

That's all!  From now on, each time that you execute your Python script,
Psyco will deploy its sophisticated algorithms so as to accelerate your
calculations.

You can see in the graphs :ref:`Figure 24 <psycoWriteComparison>` and
:ref:`Figure 25 <psycoReadComparison>` how much I/O speed improvement you can
get by using Psyco. By looking at this figures you can get an idea if these
improvements are of your interest or not. In general, if you are not going to
use compression you will take advantage of Psyco if your tables are medium
sized (from a thousand to a million rows), and this advantage will disappear
progressively when the number of rows grows well over one million. However if
you use compression, you will probably see improvements even beyond this
limit (see :ref:`compressionIssues`).
As always, there is no substitute for experimentation with your own dataset.

.. _psycoWriteComparison:

.. figure:: images/write-medium-psyco-nopsyco-comparison.png
    :align: center

    **Figure 24. Writing tables with/without Psyco.**


.. _psycoReadComparison:

.. figure:: images/read-medium-psyco-nopsyco-comparison.png
    :align: center

    **Figure 25. Reading tables with/without Psyco.**


.. _LRUOptim:

Getting the most from the node LRU cache
----------------------------------------
One limitation of the initial versions of PyTables was that they needed to
load all nodes in a file completely before being ready to deal with them,
making the opening times for files with a lot of nodes very high and
unacceptable in many cases.

Starting from PyTables 1.2 on, a new lazy node loading schema was setup that
avoids loading all the nodes of the *object tree* in memory. In addition, a
new LRU cache was introduced in order to accelerate the access to already
visited nodes. This cache (one per file) is responsible for keeping up the
most recently visited nodes in memory and discard the least recent used ones.
This represents a big advantage over the old schema, not only in terms of
memory usage (as there is no need to load *every* node in memory), but it
also adds very convenient optimizations for working interactively like, for
example, speeding-up the opening times of files with lots of nodes, allowing
to open almost any kind of file in typically less than one tenth of second
(compare this with the more than 10 seconds for files with more than 10000
nodes in PyTables pre-1.2 era) as well as optimizing the access to frequently
visited nodes. See for more info on the advantages (and also drawbacks) of
this approach.

One thing that deserves some discussion is the election of the parameter that
sets the maximum amount of nodes to be kept in memory at any time.
As PyTables is meant to be deployed in machines that can have potentially low
memory, the default for it is quite conservative (you can look at its actual
value in the :data:`parameters.NODE_CACHE_SLOTS` parameter in module
:file:`tables/parameters.py`). However, if you usually need to deal with
files that have many more nodes than the maximum default, and you have a lot
of free memory in your system, then you may want to experiment in order to
see which is the appropriate value of :data:`parameters.NODE_CACHE_SLOTS` that
fits better your needs.

As an example, look at the next code::

    def browse_tables(filename):
        fileh = open_file(filename,'a')
        group = fileh.root.newgroup
        for j in range(10):
            for tt in fileh.walk_nodes(group, "Table"):
                title = tt.attrs.TITLE
                for row in tt:
                    pass
        fileh.close()

We will be running the code above against a couple of files having a
``/newgroup`` containing 100 tables and 1000 tables respectively.  In addition,
this benchmark is run twice for two different values of the LRU cache size,
specifically 256 and 1024. You can see the results in
:ref:`table <optimization_table_1>`.

.. _optimization_table_1:

.. only:: not latex

    .. table:: **Retrieval speed and memory consumption depending on the number of nodes in LRU cache.**

        ====================== =========== === ======= ==== ==== === ======= ==== ====
        Number:                                   100 nodes             1000 nodes
        ---------------------------------- --------------------- ---------------------
        Mem & Speed                        Memory (MB) Time (ms) Memory (MB) Time (ms)
        ---------------------------------- ----------- --------- ----------- ---------
        Node is coming from... Cache size  256 1024    256  1024 256 1024    256  1024
        ====================== =========== === ======= ==== ==== === ======= ==== ====
        Disk                               14  14      1.24 1.24 51  66      1.33 1.31
        Cache                              14  14      0.53 0.52 65  73      1.35 0.68
        ====================== =========== === ======= ==== ==== === ======= ==== ====

.. raw:: latex

    \begin{threeparttable}
    \capstart\caption{Retrieval speed and memory consumption depending on the number of nodes in LRU cache.}

    \begin{tabulary}{\linewidth}{|l|l|r|r|r|r|r|r|r|r|}
    \hline
    \multicolumn{2}{|l|}{\textbf{Number:}} & \multicolumn{4}{|c|}{\textbf{100 nodes}} & \multicolumn{4}{|c|}{\textbf{1000 nodes}} \\
    \hline
    \multicolumn{2}{|l|}{\textbf{Mem and Speed}} & \multicolumn{2}{|c|}{\textbf{Memory (MB)}} & \multicolumn{2}{|c|}{\textbf{Time (ms)}}  & \multicolumn{2}{|c|}{\textbf{Memory (MB)}} & \multicolumn{2}{|c|}{\textbf{Time (ms)}}\\
    \hline
    \textbf{Node is coming from...} & \textbf{Cache size} & \textbf{256} & \textbf{1024} & \textbf{256} & \textbf{1024} & \textbf{256} & \textbf{1024} & \textbf{256} & \textbf{1024}\\
    \hline
    Disk   &  & 14 & 14 & 1.24 & 1.24 & 51 & 66 & 1.33 & 1.31 \\
    Cache  &  & 14 & 14 & 0.53 & 0.52 & 65 & 73 & 1.35 & 0.68 \\
    \hline
    \end{tabulary}

    \end{threeparttable}


From the data in :ref:`table <optimization_table_1>`, one can see that when
the number of objects that you are dealing with does fit in cache, you will
get better access times to them. Also, incrementing the node cache size
effectively consumes more memory *only* if the total nodes exceeds the slots
in cache; otherwise the memory consumption remains the same. It is also worth
noting that incrementing the node cache size in the case you want to fit all
your nodes in cache does not take much more memory than being too
conservative. On the other hand, it might happen that the speed-up that you
can achieve by allocating more slots in your cache is not worth the amount of
memory used.

Also worth noting is that if you have a lot of memory available and
performance is absolutely critical, you may want to try out a negative value
for :data:`parameters.NODE_CACHE_SLOTS`.  This will cause that all the touched
nodes will be kept in an internal dictionary and this is the faster way to
load/retrieve nodes.
However, and in order to avoid a large memory consumption, the user will be
warned when the number of loaded nodes will reach the ``-NODE_CACHE_SLOTS``
value.

Finally, a value of zero in :data:`parameters.NODE_CACHE_SLOTS` means that
any cache mechanism is disabled.

At any rate, if you feel that this issue is important for you, there is no
replacement for setting your own experiments up in order to proceed to
fine-tune the :data:`parameters.NODE_CACHE_SLOTS` parameter.

.. note::

    PyTables >= 2.3 sports an optimized LRU cache node written in C, so
    you should expect significantly faster LRU cache operations when
    working with it.


.. note::

    Numerical results reported in :ref:`table <optimization_table_1>` have been
    obtained with PyTables < 3.1. In PyTables 3.1 the node cache mechanism has
    been completely redesigned so while all comments above are still valid,
    numerical values could be a little bit different from the ones reported in
    :ref:`table <optimization_table_1>`.


Compacting your PyTables files
------------------------------
Let's suppose that you have a file where you have made a lot of row deletions
on one or more tables, or deleted many leaves or even entire subtrees. These
operations might leave *holes* (i.e. space that is not used anymore) in your
files that may potentially affect not only the size of the files but, more
importantly, the performance of I/O. This is because when you delete a lot of
rows in a table, the space is not automatically recovered on the fly.
In addition, if you add many more rows to a table than specified in the
expectedrows keyword at creation time this may affect performance as well, as
explained in :ref:`expectedRowsOptim`.

In order to cope with these issues, you should be aware that PyTables
includes a handy utility called ptrepack which can be very useful not only to
compact *fragmented* files, but also to adjust some internal parameters in
order to use better buffer and chunk sizes for optimum I/O speed.
Please check the :ref:`ptrepackDescr` for a brief tutorial on its use.

Another thing that you might want to use ptrepack for is changing the
compression filters or compression levels on your existing data for different
goals, like checking how this can affect both final size and I/O performance,
or getting rid of the optional compressors like LZO or bzip2 in your existing
files, in case you want to use them with generic HDF5 tools that do not have
support for these filters.

--------------

.. [1] CArray nodes, though not
       extensible, are chunked and have their optimum chunk size
       automatically computed at creation time, since their final shape is known.

.. [2] Except for Array objects.

.. [3] Some users reported that the typical improvement with real
       data is between a factor 1.5x and 2.5x over the already compressed
       datasets.

.. [4] So, there is not much point in using Psyco
       with recent versions of PyTables anymore.

:orphan:

=====================
PyTables User's Guide
=====================

.. raw:: latex

    \listoffigures
    \listoftables
    \clearpage


:Authors:   Francesc Alted, Ivan Vilata, Scott Prater, Vicent Mas, Tom Hedley,
            Antonio Valentino, Jeffrey Whitaker, Anthony Scopatz, Josh Moore
:Copyright: |copy| 2002, 2003, 2004 - Francesc Alted

            |copy| 2005, 2006, 2007 - Cárabos Coop. V.

            |copy| 2008, 2009, 2010 - Francesc Alted

            |copy| 2011–2021 - PyTables maintainers
:Date:      |today|
:Version:   |version|
:Home Page: http://www.pytables.org

.. raw:: latex

    \clearpage


.. rubric:: Copyright Notice and Statement for PyTables User's Guide

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

a. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

b. Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the
distribution.

c. Neither the name of Francesc Alted nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

.. |copy|   unicode:: U+000A9 .. COPYRIGHT SIGN


-------------------------
The PyTables Core Library
-------------------------

.. toctree::
    :maxdepth: 1

    introduction
    installation
    tutorials
    libref
    optimization


---------------------
Complementary modules
---------------------


.. toctree::
    :maxdepth: 1

    filenode


----------
Appendixes
----------

.. todo:: check why "latex_appendices" config option doesn't work with parts
.. todo:: try to use raw latex \appendix

.. toctree::
    :maxdepth: 1

    datatypes
    condition_syntax
    parameter_files
    utilities
    file_format


.. raw:: latex

    \bookmarksetup{startatroot}
    \addtocontents{toc}{\bigskip}


.. toctree::
    :maxdepth: 1

    bibliography
.. currentmodule:: tables

Top-level variables and functions
=================================

Global variables
----------------

.. autodata:: __version__

.. autodata:: hdf5_version


Global functions
----------------
.. autofunction:: copy_file

.. autofunction:: is_hdf5_file

.. autofunction:: is_pytables_file

.. autofunction:: open_file

.. autofunction:: set_blosc_max_threads

.. autofunction:: print_versions

.. autofunction:: restrict_flavors

.. autofunction:: split_type

.. autofunction:: test

.. autofunction:: which_lib_version

.. autofunction:: silence_hdf5_messages
.. currentmodule:: tables

General purpose expression evaluator class
==========================================

The Expr class
--------------
.. autoclass:: Expr

..  These are defined in the class docstring.
    Expr instance variables
    ~~~~~~~~~~~~~~~~~~~~~~~
    .. autoattribute:: Expr.append_mode
    .. autoattribute:: Expr.maindim
    .. autoattribute:: Expr.names
    .. autoattribute:: Expr.out
    .. autoattribute:: Expr.o_start
    .. autoattribute:: Expr.o_stop
    .. autoattribute:: Expr.o_step
    .. autoattribute:: Expr.shape
    .. autoattribute:: Expr.values


Expr methods
~~~~~~~~~~~~
.. automethod:: Expr.eval

.. automethod:: Expr.set_inputs_range

.. automethod:: Expr.set_output

.. automethod:: Expr.set_output_range


Expr special methods
~~~~~~~~~~~~~~~~~~~~
.. automethod:: Expr.__iter__
.. currentmodule:: tables

Link classes
============

.. _LinkClassDescr:

The Link class
--------------
.. autoclass:: tables.link.Link

..  These are defined in the class docstring
    .. autoattribute:: tables.link.Link.target

Link instance variables
~~~~~~~~~~~~~~~~~~~~~~~
.. autoattribute:: tables.link.Link._v_attrs


Link methods
~~~~~~~~~~~~
The following methods are useful for copying, moving, renaming and removing
links.

.. automethod:: tables.link.Link.copy

.. automethod:: tables.link.Link.move

.. automethod:: tables.link.Link.remove

.. automethod:: tables.link.Link.rename


.. _SoftLinkClassDescr:

The SoftLink class
------------------
.. autoclass:: tables.link.SoftLink


SoftLink special methods
~~~~~~~~~~~~~~~~~~~~~~~~
The following methods are specific for dereferencing and representing soft
links.

.. automethod:: tables.link.SoftLink.__call__

.. automethod:: tables.link.SoftLink.__str__


The ExternalLink class
----------------------
.. autoclass:: tables.link.ExternalLink

..  This is defined in the class docstring
    ExternalLink instance variables
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    .. autoattribute:: tables.link.ExternalLink.extfile


ExternalLink methods
~~~~~~~~~~~~~~~~~~~~
.. automethod:: tables.link.ExternalLink.umount


ExternalLink special methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following methods are specific for dereferencing and representing
external links.

.. automethod:: tables.link.ExternalLink.__call__

.. automethod:: tables.link.ExternalLink.__str__
.. currentmodule:: tables

Hierarchy definition classes
============================


.. _NodeClassDescr:

The Node class
--------------
.. autoclass:: Node

.. These are defined in class docstring
    .. autoattribute:: Node._v_depth
    .. autoattribute:: Node._v_file
    .. autoattribute:: Node._v_name
    .. autoattribute:: Node._v_pathname
    .. autoattribute:: Node._v_objectid (location independent)

Node instance variables - location dependent
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoattribute:: Node._v_parent


Node instance variables - location independent
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoattribute:: Node._v_attrs

.. autoattribute:: Node._v_isopen


Node instance variables - attribute shorthands
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoattribute:: Node._v_title


Node methods - hierarchy manipulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: Node._f_close

.. automethod:: Node._f_copy

.. automethod:: Node._f_isvisible

.. automethod:: Node._f_move

.. automethod:: Node._f_remove

.. automethod:: Node._f_rename


Node methods - attribute handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: Node._f_delattr

.. automethod:: Node._f_getattr

.. automethod:: Node._f_setattr


.. _GroupClassDescr:

The Group class
---------------
.. autoclass:: Group

..  These are defined in the class docstring
    Group instance variables
    ~~~~~~~~~~~~~~~~~~~~~~~~
    The following instance variables are provided in addition to those in Node
    (see :ref:`NodeClassDescr`):

    .. autoattribute:: Group._v_children
    .. autoattribute:: Group._v_groups
    .. autoattribute:: Group._v_hidden
    .. autoattribute:: Group._v_leaves
    .. autoattribute:: Group._v_links
    .. autoattribute:: Group._v_unknown

Group properties
~~~~~~~~~~~~~~~~
.. autoattribute:: Group._v_nchildren

.. autoattribute:: Group._v_filters


Group methods
~~~~~~~~~~~~~

.. important::

    *Caveat:* The following methods are documented for completeness, and they
    can be used without any problem. However, you should use the high-level
    counterpart methods in the File class (see :ref:`FileClassDescr`, because
    they are most used in documentation and examples, and are a bit more
    powerful than those exposed here.

The following methods are provided in addition to those in
Node (see :ref:`NodeClassDescr`):


.. automethod:: Group._f_close

.. automethod:: Group._f_copy

.. automethod:: Group._f_copy_children

.. automethod:: Group._f_get_child

.. automethod:: Group._f_iter_nodes

.. automethod:: Group._f_list_nodes

.. automethod:: Group._f_walk_groups

.. automethod:: Group._f_walknodes


Group special methods
~~~~~~~~~~~~~~~~~~~~~
Following are described the methods that automatically trigger actions when a
Group instance is accessed in a special way.

This class defines the :meth:`__setattr__`, :meth:`__getattr__` and
:meth:`__delattr__` methods, and they set, get and delete *ordinary Python
attributes* as normally intended. In addition to that, :meth:`__getattr__`
allows getting *child nodes* by their name for the sake of easy interaction
on the command line, as long as there is no Python attribute with the same
name. Groups also allow the interactive completion (when using readline) of
the names of child nodes. For instance::

    # get a Python attribute
    nchild = group._v_nchildren

    # Add a Table child called 'table' under 'group'.
    h5file.create_table(group, 'table', my_description)
    table = group.table          # get the table child instance
    group.table = 'foo'          # set a Python attribute

    # (PyTables warns you here about using the name of a child node.)
    foo = group.table            # get a Python attribute
    del group.table              # delete a Python attribute
    table = group.table          # get the table child instance again

.. automethod:: Group.__contains__

.. automethod:: Group.__delattr__

.. automethod:: Group.__getattr__

.. automethod:: Group.__iter__

.. automethod:: Group.__repr__

.. automethod:: Group.__setattr__

.. automethod:: Group.__str__


.. _LeafClassDescr:

The Leaf class
--------------
.. autoclass:: Leaf

..  These are defined in the class docstring
    .. _LeafInstanceVariables:

    Leaf instance variables
    ~~~~~~~~~~~~~~~~~~~~~~~
    These instance variables are provided in addition to those in Node
    (see :ref:`NodeClassDescr`):

    .. autoattribute:: Leaf.byteorder
    .. autoattribute:: Leaf.dtype
    .. autoattribute:: Leaf.extdim
    .. autoattribute:: Leaf.nrows
    .. autoattribute:: Leaf.nrowsinbuf
    .. autoattribute:: Leaf.shape


Leaf properties
~~~~~~~~~~~~~~~
.. autoattribute:: Leaf.chunkshape

.. autoattribute:: Leaf.ndim

.. autoattribute:: Leaf.filters

.. autoattribute:: Leaf.maindim

.. autoattribute:: Leaf.flavor

.. attribute:: Leaf.size_in_memory

    The size of this leaf's data in bytes when it is fully loaded into
    memory.

.. autoattribute:: Leaf.size_on_disk


Leaf instance variables - aliases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following are just easier-to-write aliases to their Node (see
:ref:`NodeClassDescr`) counterparts (indicated between parentheses):

.. autoattribute:: Leaf.attrs

.. autoattribute:: Leaf.name

.. autoattribute:: Leaf.object_id

.. autoattribute:: Leaf.title


Leaf methods
~~~~~~~~~~~~
.. automethod:: Leaf.close

.. automethod:: Leaf.copy

.. automethod:: Leaf.flush

.. automethod:: Leaf.isvisible

.. automethod:: Leaf.move

.. automethod:: Leaf.rename

.. automethod:: Leaf.remove

.. automethod:: Leaf.get_attr

.. automethod:: Leaf.set_attr

.. automethod:: Leaf.del_attr

.. automethod:: Leaf.truncate

.. automethod:: Leaf.__len__

.. automethod:: Leaf._f_close
.. currentmodule:: tables

Structured storage classes
==========================

.. _TableClassDescr:

The Table class
---------------
.. autoclass:: Table

.. These are defined in the class docstring
    .. _TableInstanceVariablesDescr:

    Table instance variables
    ~~~~~~~~~~~~~~~~~~~~~~~~
    The following instance variables are provided in addition to those in Leaf
    (see :ref:`LeafClassDescr`).  Please note that there are several col*
    dictionaries to ease retrieving information about a column directly by its
    path name, avoiding the need to walk through Table.description or
    :attr:`Table.cols`.

    .. autoattribute:: Table.coldescrs
    .. autoattribute:: Table.coldflts
    .. autoattribute:: Table.coldtypes
    .. autoattribute:: Table.colindexed
    .. autoattribute:: Table.colinstances
    .. autoattribute:: Table.colnames
    .. autoattribute:: Table.colpathnames
    .. autoattribute:: Table.cols
    .. autoattribute:: Table.coltypes
    .. autoattribute:: Table.description
    .. autoattribute:: Table.extdim
    .. autoattribute:: Table.indexed
    .. autoattribute:: Table.nrows


Table properties
~~~~~~~~~~~~~~~~
.. autoattribute:: Table.autoindex

.. autoattribute:: Table.colindexes

.. autoattribute:: Table.indexedcolpathnames

.. autoattribute:: Table.row

.. autoattribute:: Table.rowsize


Table methods - reading
~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: Table.col

.. automethod:: Table.iterrows

.. automethod:: Table.itersequence

.. automethod:: Table.itersorted

.. automethod:: Table.read

.. automethod:: Table.read_coordinates

.. automethod:: Table.read_sorted

.. automethod:: Table.__getitem__

.. automethod:: Table.__iter__


Table methods - writing
~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: Table.append

.. automethod:: Table.modify_column

.. automethod:: Table.modify_columns

.. automethod:: Table.modify_coordinates

.. automethod:: Table.modify_rows

.. automethod:: Table.remove_rows

.. automethod:: Table.remove_row

.. automethod:: Table.__setitem__


.. _TableMethods_querying:

Table methods - querying
~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: Table.get_where_list

.. automethod:: Table.read_where

.. automethod:: Table.where

.. automethod:: Table.append_where

.. automethod:: Table.will_query_use_indexing


Table methods - other
~~~~~~~~~~~~~~~~~~~~~
.. automethod:: Table.copy

.. automethod:: Table.flush_rows_to_index

.. automethod:: Table.get_enum

.. automethod:: Table.reindex

.. automethod:: Table.reindex_dirty


.. _DescriptionClassDescr:

The Description class
~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: Description

..  These are defined in the class docstring
    Description instance variables
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    .. autoattribute:: Description._v_col_objects
    .. autoattribute:: Description._v_dflts
    .. autoattribute:: Description._v_dtype
    .. autoattribute:: Description._v_dtypes
    .. autoattribute:: Description._v_is_nested
    .. autoattribute:: Description._v_itemsize
    .. autoattribute:: Description._v_name
    .. autoattribute:: Description._v_names
    .. autoattribute:: Description._v_nested_descr
    .. autoattribute:: Description._v_nested_formats
    .. autoattribute:: Description._v_nestedlvl
    .. autoattribute:: Description._v_nested_names
    .. autoattribute:: Description._v_pathname
    .. autoattribute:: Description._v_pathnames
    .. autoattribute:: Description._v_types
    .. autoattribute:: Description._v_offsets


Description methods
^^^^^^^^^^^^^^^^^^^
.. automethod:: Description._f_walk


.. _RowClassDescr:

The Row class
~~~~~~~~~~~~~
.. autoclass:: tables.tableextension.Row

..  These are defined in the class docstring
    Row instance variables
    ^^^^^^^^^^^^^^^^^^^^^^
    .. autoattribute:: tables.tableextension.Row.nrow


Row methods
^^^^^^^^^^^
.. automethod:: tables.tableextension.Row.append

.. automethod:: tables.tableextension.Row.fetch_all_fields

.. automethod:: tables.tableextension.Row.update


.. _RowSpecialMethods:

Row special methods
^^^^^^^^^^^^^^^^^^^
.. automethod:: tables.tableextension.Row.__contains__

.. automethod:: tables.tableextension.Row.__getitem__

.. automethod:: tables.tableextension.Row.__setitem__


.. _ColsClassDescr:

The Cols class
~~~~~~~~~~~~~~
.. autoclass:: Cols

..  These are defined in the class docstring
    Cols instance variables
    ^^^^^^^^^^^^^^^^^^^^^^^
    .. autoattribute:: Cols._v_colnames
    .. autoattribute:: Cols._v_colpathnames
    .. autoattribute:: Cols._v_desc


Cols properties
^^^^^^^^^^^^^^^
.. autoattribute:: Cols._v_table


Cols methods
^^^^^^^^^^^^
.. automethod:: Cols._f_col

.. automethod:: Cols.__getitem__

.. automethod:: Cols.__len__

.. automethod:: Cols.__setitem__


.. _ColumnClassDescr:

The Column class
~~~~~~~~~~~~~~~~
.. autoclass:: Column

.. These are defined in the class docstring

    .. autoattribute:: Column.descr
    .. autoattribute:: Column.name
    .. autoattribute:: Column.pathname

Column instance variables
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autoattribute:: Column.dtype

.. autoattribute:: Column.index

.. autoattribute:: Column.is_indexed

.. autoattribute:: Column.maindim

.. autoattribute:: Column.shape

.. autoattribute:: Column.table

.. autoattribute:: Column.type


Column methods
^^^^^^^^^^^^^^
.. automethod:: Column.create_index

.. automethod:: Column.create_csindex

.. automethod:: Column.reindex

.. automethod:: Column.reindex_dirty

.. automethod:: Column.remove_index


Column special methods
^^^^^^^^^^^^^^^^^^^^^^
.. automethod:: Column.__getitem__

.. automethod:: Column.__len__

.. automethod:: Column.__setitem__
.. currentmodule:: tables


Helper classes
==============
This section describes some classes that do not fit in any other
section and that mainly serve for ancillary purposes.


.. _FiltersClassDescr:

The Filters class
-----------------
.. autoclass:: Filters

..  These are defined in the class docstring.
    Filters instance variables
    ^^^^^^^^^^^^^^^^^^^^^^^^^^
    .. autoattribute:: Filters.bitshuffle
    .. autoattribute:: Filters.fletcher32
    .. autoattribute:: Filters.complevel
    .. autoattribute:: Filters.complib
    .. autoattribute:: Filters.shuffle


Filters methods
~~~~~~~~~~~~~~~
.. automethod:: Filters.copy


.. _IndexClassDescr:

The Index class
---------------
.. autoclass:: tables.index.Index

..  This is defined in the class docstring
    .. autoattribute:: tables.index.Index.nelements

Index instance variables
~~~~~~~~~~~~~~~~~~~~~~~~
.. autoattribute:: tables.index.Index.column

.. autoattribute:: tables.index.Index.dirty

.. autoattribute:: tables.index.Index.filters

.. autoattribute:: tables.index.Index.is_csi

.. attribute:: tables.index.Index.nelements

    The number of currently indexed rows for this column.


Index methods
~~~~~~~~~~~~~
.. automethod:: tables.index.Index.read_sorted

.. automethod:: tables.index.Index.read_indices


Index special methods
~~~~~~~~~~~~~~~~~~~~~
.. automethod:: tables.index.Index.__getitem__


The IndexArray class
--------------------

.. autoclass:: tables.indexes.IndexArray
    :members:


.. _EnumClassDescr:

The Enum class
--------------
.. autoclass:: tables.misc.enum.Enum


Enum special methods
~~~~~~~~~~~~~~~~~~~~
.. automethod:: Enum.__call__

.. automethod:: Enum.__contains__

.. automethod:: Enum.__eq__

.. automethod:: Enum.__getattr__

.. automethod:: Enum.__getitem__

.. automethod:: Enum.__iter__

.. automethod:: Enum.__len__

.. automethod:: Enum.__repr__


.. _UnImplementedClassDescr:

The UnImplemented class
-----------------------
.. autoclass:: UnImplemented
    :members:


The Unknown class
-----------------
.. autoclass:: Unknown
    :members:


.. _ExceptionsDescr:

Exceptions module
-----------------
In the :mod:`exceptions` module exceptions and warnings that are specific
to PyTables are declared.

.. autoexception:: HDF5ExtError
    :members:

.. autoexception:: ClosedNodeError

.. autoexception:: ClosedFileError

.. autoexception:: FileModeError

.. autoexception:: NodeError

.. autoexception:: NoSuchNodeError

.. autoexception:: UndoRedoError

.. autoexception:: UndoRedoWarning

.. autoexception:: NaturalNameWarning

.. autoexception:: PerformanceWarning

.. autoexception:: FlavorError

.. autoexception:: FlavorWarning

.. autoexception:: FiltersWarning

.. autoexception:: OldIndexWarning

.. autoexception:: DataTypeWarning

.. autoexception:: ExperimentalFeatureWarning
.. currentmodule:: tables.nodes.filenode

.. _filenode_classes:

Filenode Module
===============

.. automodule:: tables.nodes.filenode


Module constants
----------------

.. autodata:: NodeType

.. autodata:: NodeTypeVersions


Module functions
----------------

.. autofunction:: new_node

.. autofunction:: open_node

.. autofunction:: read_from_filenode

.. autofunction:: save_to_filenode


The RawPyTablesIO base class
----------------------------

.. autoclass:: RawPyTablesIO


RawPyTablesIO attributes
~~~~~~~~~~~~~~~~~~~~~~~~

.. autoattribute:: RawPyTablesIO.mode


RawPyTablesIO methods
~~~~~~~~~~~~~~~~~~~~~

.. automethod:: RawPyTablesIO.tell

.. automethod:: RawPyTablesIO.seek

.. automethod:: RawPyTablesIO.seekable

.. automethod:: RawPyTablesIO.fileno

.. automethod:: RawPyTablesIO.close

.. automethod:: RawPyTablesIO.flush

.. automethod:: RawPyTablesIO.truncate

.. automethod:: RawPyTablesIO.readable

.. automethod:: RawPyTablesIO.writable

.. automethod:: RawPyTablesIO.readinto

.. automethod:: RawPyTablesIO.readline

.. automethod:: RawPyTablesIO.write


The ROFileNode class
--------------------

.. autoclass:: ROFileNode


ROFileNode attributes
~~~~~~~~~~~~~~~~~~~~~

.. autoattribute:: ROFileNode.attrs


ROFileNode methods
~~~~~~~~~~~~~~~~~~

.. automethod:: ROFileNode.flush

.. automethod:: ROFileNode.read

.. automethod:: ROFileNode.readline

.. automethod:: ROFileNode.readlines

.. automethod:: ROFileNode.close

.. automethod:: ROFileNode.seek

.. automethod:: ROFileNode.tell

.. automethod:: ROFileNode.readable

.. automethod:: ROFileNode.writable

.. automethod:: ROFileNode.seekable

.. automethod:: ROFileNode.fileno


The RAFileNode class
--------------------

.. autoclass:: RAFileNode


RAFileNode attributes
~~~~~~~~~~~~~~~~~~~~~

.. autoattribute:: RAFileNode.attrs


RAFileNode methods
~~~~~~~~~~~~~~~~~~

.. automethod:: RAFileNode.flush

.. automethod:: RAFileNode.read

.. automethod:: RAFileNode.readline

.. automethod:: RAFileNode.readlines

.. automethod:: RAFileNode.truncate

.. automethod:: RAFileNode.write

.. automethod:: RAFileNode.writelines

.. automethod:: RAFileNode.close

.. automethod:: RAFileNode.seek

.. automethod:: RAFileNode.tell

.. automethod:: RAFileNode.readable

.. automethod:: RAFileNode.writable

.. automethod:: RAFileNode.seekable

.. automethod:: RAFileNode.fileno

.. currentmodule:: tables

Declarative classes
===================
In this section a series of classes that are meant to
*declare* datatypes that are required for creating
primary PyTables datasets are described.


.. _AtomClassDescr:

The Atom class and its descendants
----------------------------------
.. autoclass:: Atom

..  These are defined in the class docstring
    Atom instance variables
    ^^^^^^^^^^^^^^^^^^^^^^^
    .. autoattribute:: Atom.dflt
    .. autoattribute:: Atom.dtype
    .. autoattribute:: Atom.itemsize
    .. autoattribute:: Atom.kind
    .. autoattribute:: Atom.shape
    .. autoattribute:: Atom.type


Atom properties
~~~~~~~~~~~~~~~
.. autoattribute:: Atom.ndim

.. autoattribute:: Atom.recarrtype

.. autoattribute:: Atom.size


Atom methods
~~~~~~~~~~~~
.. automethod:: Atom.copy


Atom factory methods
~~~~~~~~~~~~~~~~~~~~
.. automethod:: Atom.from_dtype

.. automethod:: Atom.from_kind

.. automethod:: Atom.from_sctype

.. automethod:: Atom.from_type


Atom Sub-classes
~~~~~~~~~~~~~~~~
.. autoclass:: StringAtom
    :members:

.. autoclass:: BoolAtom
    :members:

.. autoclass:: IntAtom
    :members:

.. autoclass:: Int8Atom
    :members:

.. autoclass:: Int16Atom
    :members:

.. autoclass:: Int32Atom
    :members:

.. autoclass:: Int64Atom
    :members:

.. autoclass:: UIntAtom
    :members:

.. autoclass:: UInt8Atom
    :members:

.. autoclass:: UInt16Atom
    :members:

.. autoclass:: UInt32Atom
    :members:

.. autoclass:: UInt64Atom
    :members:

.. autoclass:: FloatAtom
    :members:

.. autoclass:: Float32Atom
    :members:

.. autoclass:: Float64Atom
    :members:

.. autoclass:: ComplexAtom
    :members:

.. autoclass:: Time32Atom
    :members:

.. autoclass:: Time64Atom
    :members:

.. autoclass:: EnumAtom
    :members:


Pseudo atoms
~~~~~~~~~~~~
Now, there come three special classes, ObjectAtom, VLStringAtom and
VLUnicodeAtom, that actually do not descend from Atom, but which goal is so
similar that they should be described here. Pseudo-atoms can only be used with
VLArray datasets (see :ref:`VLArrayClassDescr`), and they do not support
multidimensional values, nor multiple values per row.

They can be recognised because they also have kind, type and shape attributes,
but no size, itemsize or dflt ones. Instead, they have a base atom which
defines the elements used for storage.

See :file:`examples/vlarray1.py` and :file:`examples/vlarray2.py` for further
examples on VLArray datasets, including object serialization and string
management.


ObjectAtom
^^^^^^^^^^
.. autoclass:: ObjectAtom
    :members:


.. _VLStringAtom:

VLStringAtom
^^^^^^^^^^^^
.. autoclass:: VLStringAtom
    :members:


.. _VLUnicodeAtom:

VLUnicodeAtom
^^^^^^^^^^^^^
.. autoclass:: VLUnicodeAtom
    :members:


.. _ColClassDescr:

The Col class and its descendants
---------------------------------
.. autoclass:: Col

..
    Col instance variables
    ^^^^^^^^^^^^^^^^^^^^^^
   .. autoattribute:: _v_pos


Col instance variables
~~~~~~~~~~~~~~~~~~~~~~
In addition to the variables that they inherit from the Atom class, Col
instances have the following attributes.

.. attribute:: Col._v_pos

    The *relative* position of this column with regard to its column
    siblings.


Col factory methods
~~~~~~~~~~~~~~~~~~~
.. automethod:: Col.from_atom


Col sub-classes
~~~~~~~~~~~~~~~
.. autoclass:: StringCol
    :members:

.. autoclass:: BoolCol
    :members:

.. autoclass:: IntCol
    :members:

.. autoclass:: Int8Col
    :members:

.. autoclass:: Int16Col
    :members:

.. autoclass:: Int32Col
    :members:

.. autoclass:: Int64Col
    :members:

.. autoclass:: UIntCol
    :members:

.. autoclass:: UInt8Col
    :members:

.. autoclass:: UInt16Col
    :members:

.. autoclass:: UInt32Col
    :members:

.. autoclass:: UInt64Col
    :members:

.. autoclass:: Float32Col
    :members:

.. autoclass:: Float64Col
    :members:

.. autoclass:: ComplexCol
    :members:

.. autoclass:: TimeCol
    :members:

.. autoclass:: Time32Col
    :members:

.. autoclass:: Time64Col
    :members:

.. autoclass:: EnumCol
    :members:


.. _IsDescriptionClassDescr:

The IsDescription class
-----------------------
.. autoclass:: IsDescription


Description helper functions
----------------------------
.. autofunction:: tables.description.descr_from_dtype

.. autofunction:: tables.description.dtype_from_descr


.. _AttributeSetClassDescr:

The AttributeSet class
----------------------
.. autoclass:: tables.attributeset.AttributeSet

..  These are defined in the class docstring
    AttributeSet attributes
    ~~~~~~~~~~~~~~~~~~~~~~~
    .. autoattribute:: tables.attributeset.AttributeSet._v_attrnames
    .. autoattribute:: tables.attributeset.AttributeSet._v_attrnamessys
    .. autoattribute:: tables.attributeset.AttributeSet._v_attrnamesuser
    .. autoattribute:: tables.attributeset.AttributeSet._v_unimplemented

AttributeSet properties
~~~~~~~~~~~~~~~~~~~~~~~
.. autoattribute:: tables.attributeset.AttributeSet._v_node


AttributeSet methods
~~~~~~~~~~~~~~~~~~~~
.. automethod:: tables.attributeset.AttributeSet._f_copy

.. automethod:: tables.attributeset.AttributeSet._f_list

.. automethod:: tables.attributeset.AttributeSet._f_rename

.. automethod:: tables.attributeset.AttributeSet.__contains__
.. currentmodule:: tables

File manipulation class
=======================

.. _FileClassDescr:

The File Class
--------------
.. autoclass:: File

..  These are defined in the class docstring.
    This is necessary because attributes created in a class's
    __init__ method can't be documented with autoattribute.
    See Sphinx bug #904.
    https://bitbucket.org/birkenfeld/sphinx/issue/904

    Attributes
    ~~~~~~~~~~
    .. autoattribute:: File.filename
    .. autoattribute:: File.format_version
    .. autoattribute:: File.isopen
    .. autoattribute:: File.mode
    .. autoattribute:: File.root
    .. autoattribute:: File.root_uep


File properties
~~~~~~~~~~~~~~~
.. autoattribute:: File.title

.. autoattribute:: File.filters


File methods - file handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: File.close

.. automethod:: File.copy_file

.. automethod:: File.flush

.. automethod:: File.fileno

.. automethod:: File.__enter__

.. automethod:: File.__exit__

.. automethod:: File.__str__

.. automethod:: File.__repr__

.. automethod:: File.get_file_image

.. automethod:: File.get_filesize

.. automethod:: File.get_userblock_size


File methods - hierarchy manipulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: File.copy_children

.. automethod:: File.copy_node

.. automethod:: File.create_array

.. automethod:: File.create_carray

.. automethod:: File.create_earray

.. automethod:: File.create_external_link

.. automethod:: File.create_group

.. automethod:: File.create_hard_link

.. automethod:: File.create_soft_link

.. automethod:: File.create_table

.. automethod:: File.create_vlarray

.. automethod:: File.move_node

.. automethod:: File.remove_node

.. automethod:: File.rename_node


File methods - tree traversal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: File.get_node

.. automethod:: File.is_visible_node

.. automethod:: File.iter_nodes

.. automethod:: File.list_nodes

.. automethod:: File.walk_groups

.. automethod:: File.walk_nodes

.. automethod:: File.__contains__

.. automethod:: File.__iter__


File methods - Undo/Redo support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: File.disable_undo

.. automethod:: File.enable_undo

.. automethod:: File.get_current_mark

.. automethod:: File.goto

.. automethod:: File.is_undo_enabled

.. automethod:: File.mark

.. automethod:: File.redo

.. automethod:: File.undo


File methods - attribute handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automethod:: File.copy_node_attrs

.. automethod:: File.del_node_attr

.. automethod:: File.get_node_attr

.. automethod:: File.set_node_attr
.. currentmodule:: tables

Homogenous storage classes
==========================

.. _ArrayClassDescr:

The Array class
---------------
.. autoclass:: Array


Array instance variables
~~~~~~~~~~~~~~~~~~~~~~~~
.. attribute:: Array.atom

    An Atom (see :ref:`AtomClassDescr`) instance representing the *type*
    and *shape* of the atomic objects to be saved.

.. autoattribute:: Array.rowsize

.. attribute:: Array.nrow

    On iterators, this is the index of the current row.

.. autoattribute:: Array.nrows


Array methods
~~~~~~~~~~~~~
.. automethod:: Array.get_enum

.. automethod:: Array.iterrows

.. automethod:: Array.__next__

.. automethod:: Array.read


Array special methods
~~~~~~~~~~~~~~~~~~~~~
The following methods automatically trigger actions when an :class:`Array`
instance is accessed in a special way (e.g. ``array[2:3,...,::2]`` will be
equivalent to a call to
``array.__getitem__((slice(2, 3, None), Ellipsis, slice(None, None, 2))))``.

.. automethod:: Array.__getitem__

.. automethod:: Array.__iter__

.. automethod:: Array.__setitem__


.. _CArrayClassDescr:

The CArray class
----------------
.. autoclass:: CArray


.. _EArrayClassDescr:

The EArray class
----------------
.. autoclass:: EArray


.. _EArrayMethodsDescr:

EArray methods
~~~~~~~~~~~~~~

.. automethod:: EArray.append


.. _VLArrayClassDescr:

The VLArray class
-----------------
.. autoclass:: VLArray

..  These are defined in the class docstring
    VLArray instance variables
    ~~~~~~~~~~~~~~~~~~~~~~~~~~
    .. autoattribute:: VLArray.atom
    .. autoattribute:: VLArray.flavor
    .. autoattribute:: VLArray.nrow
    .. autoattribute:: VLArray.nrows
    .. autoattribute:: VLArray.extdim
    .. autoattribute:: VLArray.nrows


VLArray properties
~~~~~~~~~~~~~~~~~~
.. autoattribute:: VLArray.size_on_disk

.. autoattribute:: VLArray.size_in_memory


VLArray methods
~~~~~~~~~~~~~~~
.. automethod:: VLArray.append

.. automethod:: VLArray.get_enum

.. automethod:: VLArray.iterrows

.. automethod:: VLArray.__next__

.. automethod:: VLArray.read

.. automethod:: VLArray.get_row_size


VLArray special methods
~~~~~~~~~~~~~~~~~~~~~~~
The following methods automatically trigger actions when a :class:`VLArray`
instance is accessed in a special way (e.g., vlarray[2:5] will be equivalent
to a call to vlarray.__getitem__(slice(2, 5, None)).

.. automethod:: VLArray.__getitem__

.. automethod:: VLArray.__iter__

.. automethod:: VLArray.__setitem__
=======================================
 Release notes for PyTables 3.5 series
=======================================

:Author: PyTables Developers
:Contact: pytables-dev@googlegroups.com

.. py:currentmodule:: tables


Changes from 3.5.1 to 3.5.2
===========================

- Fixed compatibility with python 3.8: Fixed `Dictionary keys changed during
  iteration` RuntimeError while moving/renaming a node.
  Thanks to Christoph Gohlke for reporting and Miro Hrončok for help with
  building PyTables for python 3.8alpha (cython compatibility).
  see :issue:`733` and PR #737.
- Fixed a bug in offset calculations producing floats instead of ints
  affecting python 3. See PR #736. Thanks to Brad Montgomery.


Changes from 3.5.0 to 3.5.1
===========================

- Maintenance release to fix how PyPi repo is handling wheel versions.


Changes from 3.4.4 to 3.5.0
===========================

Improvements
------------
 - When copying data from native HDF5 files with padding in compound types,
   the padding is not removed now by default.  This allows for better
   compatibility with existing HDF5 applications that expect the padding
   to stay there.
   Also, when the `description` is a NumPy struct array with padding, this
   is honored now.  The previous behaviour (i.e. getting rid of paddings) can
   be replicated by passing the new `allow_padding` parameter when opening
   a file.  For some examples, see the new `examples/tables-with-padding.py`
   and `examples/attrs-with-padding.py`.  For details on the implementation
   see :issue:`720`.
 - Added a new flag `--dont-allow-padding` in `ptrepack` utility so as to
   replicate the previous behaviour of removing padding during file copies.
   The default is to honor the original padding in copies.
 - Improve compatibility with numpy 1.16.
 - Improve detection of the LZO2 library at build time.
 - Suppress several warnings.
 - Add AVX2 support for Windows.  See PR #716.  Thanks to Robert McLeod.
What's new in PyTables 0.9
==========================

On this release you will find a series of quite
exciting new features, being the most important the indexing
capabilities, in-kernel selections, support for complex datatypes and
the possibility to modify values in both tables *and* arrays (yeah,
finally :).

New features:
-------------

- Indexing of columns in tables. That allow to make data selections on
  tables up to 500 times faster than standard selections (for
  ex. doing a selection along an indexed column of 100 million of rows
  takes less than 1 second on a modern CPU). Perhaps the most
  interesting thing about the indexing algorithm implemented by
  PyTables is that the time taken to index grows *linearly* with the
  length of the data, so, making the indexation process to be
  *scalable* (quite differently to many relational databases). This
  means that it can index, in a relatively quick way, arbitrarily
  large table columns (for ex. indexing a column of 100 million of rows
  takes just 100 seconds, i.e. at a rate of 1 Mrow/sec). See more
  detailed info about that in http://www.pytables.org/docs/SciPy04.pdf.

- In-kernel selections. This feature allow to make data selections on
  tables up to 5 times faster than standard selections (i.e. pre-0.9
  selections), without a need to create an index. As a hint of how
  fast these selections can be, they are up to 10 times faster than a
  traditional relational database. Again, see
  http://www.pytables.org/docs/SciPy04.pdf for some experiments on that
  matter.

- Support of complex datatypes for all the data objects (i.e. Table,
  Array, EArray and VLArray). With that, the complete set of datatypes
  of Numeric and numarray packages are supported. Thanks to Tom Hedley
  for providing the patches for Array, EArray and VLArray objects, as
  well as updating the User's Manual and adding unit tests for the new
  functionality.

- Modification of values. You can modify Table, Array, EArray and
  VLArray values. See Table.modifyRows, Table.modifyColumns() and the
  newly introduced __setitem__() method for Table, Array, EArray and
  VLArray entities in the Library Reference of User's Manual.

- A new sub-package called "nodes" is there. On it, there will be
  included different modules to make more easy working with different
  entities (like images, files, ...). The first module that has been
  added to this sub-package is "FileNode", whose mission is to enable
  the creation of a database of nodes which can be used like regular
  opened files in Python.  In other words, you can store a set of
  files in a PyTables database, and read and write it as you would do
  with any other file in Python. Thanks to Ivan Vilata i Balaguer for
  contributing this.

Improvements:
-------------

- New __len__(self) methods added in Arrays, Tables and Columns. This,
  in combination with __getitem__(self,key) allows to better emulate
  sequences.

- Better capabilities to import generic HDF5 files. In particular,
  Table objects (in the HDF5_HL naming schema) with "holes" in their
  compound type definition are supported. That allows to read certain
  files produced by NASA (thanks to Stephen Walton for reporting this).

- Much improved test units. More than 2000 different tests has been
  implemented which accounts for more than 13000 loc (this represents
  twice of the PyTables library code itself (!)).

Backward-incompatible API changes:
----------------------------------

- The __call__ special method has been removed from objects File,
  Group, Table, Array, EArray and VLArray. Now, you should use
  walkNodes() in File and Group and iterrows in Table, Array, EArray
  and VLArray to get the same functionality. This would provide better
  compatibility with IPython as well.

'nctoh5', a new importing utility:

- Jeff Whitaker has contributed a script to easily convert NetCDF
  files into HDF5 files using Scientific Python and PyTables. It has
  been included and documented as a new utility.

Bug fixes:
----------

- A call to File.flush() now invoke a call to H5Fflush() so to
  effectively flushing all the file contents to disk. Thanks to Shack
  Toms for reporting this and providing a patch.

- SF #1054683: Security hole in utils.checkNameValidity(). Reported in
  2004-10-26 by ivilata

- SF #1049297: Suggestion: new method File.delAttrNode(). Reported in
  2004-10-18 by ivilata

- SF #1049285: Leak in AttributeSet.__delattr__(). Reported in
  2004-10-18 by ivilata

- SF #1014298: Wrong method call in examples/tutorial1-2.py. Reported
  in 2004-08-23 by ivilata

- SF #1013202: Cryptic error appending to EArray on RO file. Reported
  in 2004-08-21 by ivilata

- SF #991715: Table.read(field="var1", flavor="List") fails. Reported
  in 2004-07-15 by falted

- SF #988547: Wrong file type assumption in File.__new__. Reported in
  2004-07-10 by ivilata


Bon profit!,

-- Francesc Altet
falted@pytables.org

==============================
 What's new in PyTables 1.3.2
==============================


:Author: Francesc Altet
:Contact: faltet@carabos.com
:Author: Ivan Vilata i Balaguer
:Contact: ivilata@carabos.com


This document details the modifications to PyTables since version 1.2.  Its
main purpose is help you ensure that your programs will be runnable when you
switch from PyTables 1.2 to PyTables 1.3.2.


API additions
=============

- The ``Table.Cols`` accessor has received a new ``__setitem__()`` method that
  allows doing things like::

      table.cols[4] = record
      table.cols.x[4:1000:2] = array   # homogeneous column
      table.cols.Info[4:1000:2] = recarray   # nested column


Backward-incompatible changes
=============================

- None


Deprecated features
===================

- None


API refinements
===============

- ``Table.itersequence()`` has changed the default value for the ``sort``
  parameter.  It is now false by default, as it is not clear if this actually
  accelerates the iterator, so it is better to let the user do the proper
  checks (if interested).


Bug fixes (affecting API)
=========================

- None


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: text
.. coding: utf-8
.. fill-column: 78
.. End:
What's new in PyTables 0.7.2
----------------------------

This is a mainly a maintenance release, where the next issues has
been addressed:

- Fixed a nasty memory leak located on the C libraries (It was
  occurring during attribute writes). After that, the memory
  consumption when using large object trees has dropped quite
  a bit. However, there remains some small leaks that has been
  tracked down to the underlying numarray library. These leaks
  has been reported, and hopefully they should be fixed more
  sooner than later.

- Table buffers are built dinamically now, so if Tables are
  not accessed for reading or writing this memory will not be
  booked. This will help to reduce the memory consumption.

- The opening of files with lots of nodes has been optimized
  between a factor 2 and 3. For example, a file with 10 groups
  and 3000 tables that takes 9.3 seconds to open in 0.7.1, now
  takes only 2.8 seconds.

- The Table.read() method has been refactored and optimized
  and some parts of its code has been moved to Pyrex. In
  particular, in the special case of step=1, up to a factor 5
  of speedup (reaching 160 MB/s on a Pentium4 @ 2 GHz) when
  reading table contents can be achieved now.


Enjoy!,

-- Francesc Alted
falted@openlc.org

==============================
 What's new in PyTables 1.2.3
==============================


:Author: Francesc Altet
:Contact: faltet@carabos.com
:Author: Ivan Vilata i Balaguer
:Contact: ivilata@carabos.com


This document details the modifications to PyTables since version 1.2.  Its
main purpose is help you ensure that your programs will be runnable when you
switch from PyTables 1.2 to PyTables 1.2.3.


API additions
=============

- None

Backward-incompatible changes
=============================

- None

Deprecated features
===================

- None


API refinements
===============

- None


Bug fixes (affecting API)
=========================

- None


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: text
.. coding: utf-8
.. fill-column: 78
.. End:
=======================================
 Release notes for PyTables 2.2 series
=======================================

:Author: Francesc Alted i Abad
:Contact: faltet@pytables.org


Changes from 2.2.1rc1 to 2.2.1
==============================

- The `Row` accessor implements a new `__contains__` special method that
  allows doing things like::

    for row in table:
        if item in row:
            print "Value found in row", row.nrow
            break

  Closes #309.

- PyTables is more friendly with easy_install and pip now, as all the
  Python dependencies should be installed automatically.  Closes #298.


Changes from 2.2 to 2.2.1rc1
============================

- When using `ObjectAtom` objects in `VLArrays` the ``HIGHEST_PROTOCOL``
  is used for pickling objects.  For NumPy arrays, this simple change
  leads to space savings up to 3x and time improvements up to 30x.
  Closes #301.

- tables.Expr can perform operations on scalars now.  Thanks to Gaëtan
  de Menten for providing a patch for this.  Closes #287.

- Fixed a problem with indexes larger than 32-bit on leaf objects on
  32-bit machines.  Fixes #283.

- Merged in Blosc 1.1.2 for fixing a problem with large datatypes and
  subprocess issues.  Closes #288 and #295.

- Due to the adoption of Blosc 1.1.2, the pthreads-win32 library
  dependency is dropped on Windows platforms.

- Fixed a problem with tables.Expr and operands with vary large
  rowsizes. Closes #300.

- ``leaf[numpy.array[scalar]]`` idiom returns a NumPy array instead of
  an scalar.  This has been done for compatibility with NumPy.  Closes
  #303.

- Optimization for `Table.copy()` so that ``FIELD_*`` attrs are not
  overwritten during the copy.  This can lead to speed-ups up to 100x
  for short tables that have hundreds of columns.  Closes #304.

- For external links, its relative paths are resolved now with respect
  to the directory of the main HDF5 file, rather than with respect to
  the current directory.  Closes #306.

- ``Expr.setInputsRange()`` and ``Expr.setOutputRange()`` do support
  ``numpy.integer`` types now.  Closes #285.

- Column names in tables can start with '__' now.  Closes #291.

- Unicode empty strings are supported now as attributes.  Addresses #307.

- Cython 0.13 and higher is supported now.  Fixes #293.

- PyTables should be more 'easy_install'-able now.  Addresses #298.


Changes from 2.2rc2 to 2.2 (final)
==================================

- Updated Blosc to 1.0 (final).

- Filter ID of Blosc changed from wrong 32010 to reserved 32001.  This
  will prevent PyTables 2.2 (final) to read files created with Blosc and
  PyTables 2.2 pre-final.  `ptrepack` can be used to retrieve those
  files, if necessary.  More info in ticket #281.

- Recent benchmarks suggest a new parametrization is better in most
  scenarios:

  * The default chunksize has been doubled for every dataset size.  This
    works better in most of scenarios, specially with the new Blosc
    compressor.

  * The HDF5 CHUNK_CACHE_SIZE parameter has been raised to 2 MB in order
    to better adapt to the chunksize increase.  This provides better hit
    ratio (at the cost of consuming more memory).

  Some plots have been added to the User's Manual (chapter 5) showing
  how the new parametrization works.


Changes from 2.2rc1 to 2.2rc2
=============================

- A new version of Blosc (0.9.5) is included.  This version is now
  considered to be stable and apt for production.  Thanks for all
  PyTables users that have contributed to find and report bugs.

- Added a new `IO_BUFFER_SIZE` parameter to ``tables/parameters.py``
  that allows to set the internal PyTables' buffer for doing I/O.  This
  replaces `CHUNKTIMES` but it is more general because it affects to all
  `Leaf` objects and also the `tables.Expr` module (and not only tables
  as before).

- `BUFFERTIMES` parameter in ``tables/parameters.py`` has been
  renamed to `BUFFER_TIMES` which is more consistent with other
  parameter names.

- On Windows platforms, the path to the tables module is now appended to
  sys.path and the PATH environment variable. That way DLLs and PYDs in
  the tables directory are to be found now.  Thanks to Christoph Gohlke
  for the hint.

- A replacement for barriers for Mac OSX, or other systems not
  implementing them, has been carried out.  This allows to compile
  PyTables on such platforms.  Fixes #278

- Fixed a couple of warts that raise compatibility warnings with
  forthcoming Python 2.7.

-  HDF5 1.8.5 is used in Windows binaries.

Changes from 2.2b3 to 2.2rc1
============================

- Numexpr is not included anymore in PyTables and has become a requisite
  instead.  This is because Numexpr already has decent enough installers
  and is available in the PyPI repository also, so it should be easy for
  users to fulfill this dependency.

- When using a Numexpr package that is turbo-loaded with Intel's
  VML/MKL, the parameter `MAX_THREADS` will control the number of
  threads that VML can use during computations.  For a finer control,
  the `numexpr.set_vml_num_threads()` can always be used.

- Cython is used now instead of Pyrex for Pyrex extensions.

- Updated to 0.9 version of Blosc compressor.  This version can make use
  of threads so as to accelerate the compression/decompression process.
  In order to change the maximum number of threads that Blosc can use (2
  by default), you can modify the `MAX_THREADS` variable in
  ``tables/parameters.py`` or make use of the new `setBloscMaxThreads()`
  global function.

- Reopening already opened files is supported now, provided that there is
  not incompatibility among intended usages (for example, you cannot
  reopen in append mode an already opened file in read-only mode).

- Option ``--print-versions`` for ``test_all.py`` script is now
  preferred over the deprecated ``--show-versions``.  This is more
  consistent with the existing `print_versions()` function.

- Fixed a bug that, under some circumstances, prevented the use of table
  iterators in `itertool.groupby()`.  Now, you can safely do things
  like::

    sel_rows = table.where('(row_id >= 3)')
    for group_id, grouped_rows in itertools.groupby(sel_rows, f_group):
        group_mean = average([row['row_id'] for row in grouped_rows])

  Fixes #264.

- Copies of `Array` objects with multidimensional atoms (coming from
  native HDF5 files) work correctly now (i.e. the copy holds the atom
  dimensionality).  Fixes #275.

- The `tables.openFile()` function does not try anymore to open/close
  the file in order to guess whether it is a HDF5 or PyTables one before
  opening it definitely.  This allows the `fcntl.flock()` and
  `fcntl.lockf()` Python functions to work correctly now (that's useful
  for arbitrating access to the file by different processes).  Thanks to
  Dag Sverre Seljebotn and Ivan Vilata for their suggestions on hunting
  this one!  Fixes #185.

- The estimation of the chunksize when using multidimensional atoms in
  EArray/Carray was wrong because it did not take in account the shape
  of the atom.  Thanks to Ralf Juengling for reporting.  Fixes #273.

- Non-contiguous arrays can now safely be saved as attributes.  Before,
  if arrays were not contiguous, incorrect data was saved in attr.
  Fixes #270.

- EXTDIM attribute for CArray/EArray now saves the correct extendable
  dimension, instead of rubbish.  This does not affected functionality,
  because extendable dimension was retrieved directly from shape
  information, but it was providing misleading information to the user.
  Fixes #268.

API changes
-----------

- Now, `Table.Cols.__len__()` returns the number of top level columns
  instead of the number of rows in table.  This is more consistent in
  that `Table.Cols` is an accessor for *columns*.  Fixes #276.


Changes from 2.2b2 to 2.2b3
===========================

- Blosc compressor has been added as an additional filter, in addition
  to the existing Zlib, LZO and bzip2.  This new compressor is meant for
  fast compression and extremely fast decompression.  Fixes #265.

- In `File.copyFile()` method, `copyuserattrs` was set to false as
  default.  This was inconsistent with other methods where the default
  value for `copyuserattrs` is true.  The default for this is true now.
  Closes #261.

- `tables.copyFile` and `File.copyFile` recognize now the parameters
  present in ``tables/parameters.py``.  Fixes #262.

- Backported fix for issue #25 in Numexpr (OP_NEG_LL treats the argument
  as an int, not a long long).  Thanks to David Cooke for this.

- CHUNK_CACHE_NELMTS in `tables/parameters.py` set to a prime number as
  Neil Fortner suggested.

- Workaround for a problem in Python 2.6.4 (and probably other versions
  too) for pickling strings like "0" or "0.".  Fixes #253.


Changes from 2.2b1 to 2.2b2
===========================

Enhancements
------------

- Support for HDF5 hard links, soft links and external links (when
  PyTables is compiled against HDF5 1.8.x series).  A new tutorial about
  its usage has been added to the 'Tutorials' chapter of User's Manual.
  Closes #239 and #247.

- Added support for setting HDF5 chunk cache parameters in file
  opening/creating time.  'CHUNK_CACHE_NELMTS', 'CHUNK_CACHE_PREEMPT'
  and 'CHUNK_CACHE_SIZE' are the new parameters.  See "PyTables'
  parameter files" appendix in User's Manual for more info.  Closes
  #221.

- New `Unknown` class added so that objects that HDF5 identifies as
  ``H5G_UNKNOWN`` can be mapped to it and continue operations
  gracefully.

- Added flag `--dont-create-sysattrs` to ``ptrepack`` so as to not
  create sys attrs (default is to do it).

- Support for native compound types in attributes.  This allows for
  better compatibility with HDF5 files.  Closes #208.

- Support for native NumPy dtype in the description parameter of
  `File.createTable()`.  Closes #238.


Bugs fixed
----------

- Added missing `_c_classId` attribute to the `UnImplemented` class.
  ``ptrepack`` no longer chokes while copying `Unimplemented` classes.

- The ``FIELD_*`` sys attrs are no longer copied when the
  ``PYTABLES_SYS_ATTRS`` parameter is set to false.

- `File.createTable()` no longer segfaults if description=None.  Closes
  #248.

- Workaround for avoiding a Python issue causing a segfault when saving
  and then retrieving a string attribute with values "0" or "0.".
  Closes #253.


API changes
-----------

- `Row.__contains__()` disabled because it has little sense to query for
  a key in Row, and the correct way should be to query for it in
  `Table.colnames` or `Table.colpathnames` better.  Closes #241.

- [Semantic change] To avoid a common pitfall when asking for the string
  representation of a `Row` class, `Row.__str__()` has been redefined.
  Now, it prints something like::

      >>> for row in table:
      ...     print row
      ...
      /newgroup/table.row (Row), pointing to row #0
      /newgroup/table.row (Row), pointing to row #1
      /newgroup/table.row (Row), pointing to row #2

  instead of::

      >>> for row in table:
      ...     print row
      ...
      ('Particle:      0', 0, 10, 0.0, 0.0)
      ('Particle:      1', 1, 9, 1.0, 1.0)
      ('Particle:      2', 2, 8, 4.0, 4.0)

  Use `print row[:]` idiom if you want to reproduce the old behaviour.
  Closes #252.


Other changes
-------------

- After some improvements in both HDF5 and PyTables, the limit before
  emitting a `PerformanceWarning` on the number of children in a group
  has been raised from 4096 to 16384.


Changes from 2.1.1 to 2.2b1
===========================

Enhancements
------------

- Added `Expr`, a class for evaluating expressions containing
  array-like objects.  It can evaluate expressions (like '3*a+4*b')
  that operate on arbitrary large arrays while optimizing the
  resources (basically main memory and CPU cache memory) required to
  perform them.  It is similar to the Numexpr package, but in addition
  to NumPy objects, it also accepts disk-based homogeneous arrays,
  like the `Array`, `CArray`, `EArray` and `Column` PyTables objects.

- Added support for NumPy's extended slicing in all `Leaf` objects.
  With that, you can do the next sort of selections::

      array1 = array[4]                       # simple selection
      array2 = array[4:1000:2]                # slice selection
      array3 = array[1, ..., ::2, 1:4, 4:]    # general slice selection
      array4 = array[1, [1,5,10], ..., -1]    # fancy selection
      array5 = array[np.where(array[:] > 4)]  # point selection
      array6 = array[array[:] > 4]            # boolean selection

  Thanks to Andrew Collette for implementing this for h5py, from which
  it has been backported.  Closes #198 and #209.

- Numexpr updated to 1.3.1.  This can lead to up a 25% improvement of
  the time for both in-kernel and indexed queries for unaligned
  tables.

- HDF5 1.8.3 supported.


Bugs fixed
----------

- Fixed problems when modifying multidimensional columns in Table
  objects.  Closes #228.

- Row attribute is no longer stalled after a table move or rename.
  Fixes #224.

- Array.__getitem__(scalar) returns a NumPy scalar now, instead of a
  0-dim NumPy array.  This should not be noticed by normal users,
  unless they check for the type of returned value.  Fixes #222.


API changes
-----------

- Added a `dtype` attribute for all leaves.  This is the NumPy
  ``dtype`` that most closely matches the leaf type.  This allows for
  a quick-and-dirty check of leaf types.  Closes #230.

- Added a `shape` attribute for `Column` objects.  This is formed by
  concatenating the length of the column and the shape of its type.
  Also, the representation of columns has changed an now includes the
  length of the column as the leading dimension.  Closes #231.

- Added a new `maindim` attribute for `Column` which has the 0 value
  (the leading dimension).  This allows for a better similarity with
  other \*Array objects.

- In order to be consistent and allow the extended slicing to happen
  in `VLArray` objects too, `VLArray.__setitem__()` is not able to
  partially modify rows based on the second dimension passed as key.
  If this is tried, an `IndexError` is raised now.  Closes #210.

- The `forceCSI` flag has been replaced by `checkCSI` in the next
  `Table` methods: `copy()`, `readSorted()` and `itersorted()`.  The
  change reflects the fact that a re-index operation cannot be
  triggered from these methods anymore.  The rational for the change
  is that an indexing operation is a potentially very expensive
  operation that should be carried out explicitly instead of being
  triggered by methods that should not be in charge of this task.
  Closes #216.


Backward incompatible changes
-----------------------------

- After the introduction of the `shape` attribute for `Column`
  objects, the shape information for multidimensional columns has been
  removed from the `dtype` attribute (it is set to the base type of
  the column now).  Closes #232.


  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: rst
.. coding: utf-8
.. fill-column: 72
.. End:
PyTables 0.7.1 is out!
----------------------

This is a mainly a bug-fixing release, where the next problems has
been addressed:

- Fixed several memory leaks. After that, the memory
  consumption when using large object trees has dropped
  sensibly. However, there remains some small leaks, but
  hopefully they are not very important unless you use *huge*
  object trees.

- Fixed a bug that make the __getitem__ special method in
  table to fail when the stop parameter in a extended slice
  was not specified. That is, table[10:] now correctly returns
  table[10:table.nrows+1], and not table[10:11].

- The removeRows() method in Table did not update the NROWS
  attribute in Table objects, giving place to errors after
  doing further updating operations (removing or adding more
  rows) in the same table. This has been fixed now.

Apart of these fixes, a new lazy reading algorithm for attributes has
been activated by default. With that, the opening of objects with
large hierarchies has been improved by 60% (you can obtain another
additional 10% if using python 2.3 instead of python 2.2).  The
documentation has been updated as well, specially a more detailed
instructions on the compression (zlib) libraries installation.

Also, a stress test has been conducted in order to see if PyTables can
*really* work not only with large data tables, but also with large
object trees. In it, it has been generated and checked a file with
more than 1 TB of size and more than 100 thousand tables on it!.
=======================================
 Release notes for PyTables 2.3 series
=======================================

:Author: PyTables maintainers
:Contact: pytables@googlemail.com


Changes from 2.3 to 2.3.1
=========================

- Fixed a bug that prevented to read scalar datasets of UnImplemented types
  (closes :issue:`111`). Thanks to Kamil Kisiel.

- Fixed a bug in `setup.py` that caused installation of PyTables 2.3 to fail
  on hosts with multiple python versions installed (closes :issue:`113`).
  Thanks to sbinet.


Changes from 2.2.1 to 2.3
=========================

Features coming from (now liberated) PyTables Pro
-------------------------------------------------

- OPSI is a powerful and innovative indexing engine allowing PyTables to
  perform fast queries on arbitrarily large tables. Moreover, it offers a wide
  range of optimization levels for its indexes so that the user can choose the
  best one that suits her needs (more or less size, more or less performance).
  Indexation code also takes advantage of the vectorization capabilities of the
  NumPy and Numexpr packages to ensure really short indexing and search times.

- A fine-tuned LRU cache for both metadata (nodes) and regular data that lets
  you achieve maximum speed for intensive object tree browsing during data
  reads and queries. It complements the already efficient cache present in
  HDF5, although this is more geared towards high-level structures that are
  specific to PyTables and that are critical for achieving very high
  performance.

Other changes
-------------

- Indexes with no elements are now evaluated as non-CSI ones.  Closes
  #312.

- Numexpr presence is tested now in setup.py, provided that user is not
  using setuptools (i.e. ``easy_install`` or ``pip`` tools).  When using
  setuptools, numexpr continues to be a requisite (and Cython too).
  Closes #298.

- Cython is enforced now during compilation time.  Also, it is not
  required when running tests.

- Repeatedly closing a file that has been reopened several times is
  supported now.  Closes #318.

- The number of times a file has been currently reopened is available
  now in the new `File.open_count` read-only attribute.

- The entire documentation set has been converted to sphinx (close
  :issue:`85` and :issue:`86`) that now also has an index
  (closes :issue`39`).

- The entire test suite has been updated to use unittest specific
  assertions (closes :issue:`66`).

- PyTables has been tested against the latest version of numpy (v. 1.6.1
  and 2.0dev) and Cython (v, 0.15) packages. Closes :issue:`84`.

- The setup.py script has been improved to better detect runtimes
  (closes :issue:`73`).

Deprecations
------------

Support for some old packages and related features has been deprecated
and will be removed in future versions:

- Numeric (closes :issue:`76`)
- numarray (closes :issue`76` and :issue:`75`)
- HDF5 1.6.x (closes :issue`96`)

At the API level the following are now deprecated:

- the tables.is_pro constant is deprecated because PyTables Pro
  has been released under an open source license.
- the netcdf3 sub-package (closes :issue:`67`)
- the nra sub-package


  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: rst
.. coding: utf-8
.. fill-column: 72
.. End:
============================
 What's new in PyTables 1.1
============================


:Author: Francesc Altet
:Contact: faltet@carabos.com
:Author: Ivan Vilata i Balaguer
:Contact: ivilata@carabos.com


This document details the modifications to PyTables since version 1.0.  Its
main purpose is help you ensure that your programs will be runnable when you
switch from PyTables 1.0 to PyTables 1.1.


API additions
=============

- something...

Backward-incompatible changes
=============================

- ``Table.read()`` raises a ``KeyError`` instead of a ``ValueError`` when a
  nonexistent field name is specified, for consistency with other methods.
  The same goes for the ``col()`` method.

- ``File.__contains__()`` returns a true value when it is asked for an existent
  node, be it visible or not.  This is more consistent with
  ``Group.__contains__()``.


API refinements
===============

- Using ``table.cols['colname']`` is deprecated.  The usage of
  ``table.cols._f_col('colname')`` (with the new ``Cols._f_col()`` method) is
  preferred.

Bug fixes (affecting API)
=========================

- something...


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: text
.. coding: utf-8
.. fill-column: 78
.. End:
=======================================
 Release notes for PyTables 3.4 series
=======================================

Changes from 3.4.3 to 3.4.4
===========================

Improvements
------------
 - Environment variable to control the use of embedded libraries.
   Thanks to avalentino.
 - Include citation in repository.
   :issue:`690`. Thanks to katrinleinweber.

Bugs fixed
----------
 - Fixed import error with numexpr 2.6.5.dev0
   :issue:`685`. Thanks to cgohlke.
 - Fixed linter warnings.
   Thanks to avalentino.
 - Fixed for re.split() is version detection.
   :issue:`687`. Thanks to mingwandroid.
 - Fixed test failures with Python 2.7 and NumPy 1.14.3
   :issue:`688` & :issue:`689`. Thanks to oleksandr-pavlyk.


Changes from 3.4.2 to 3.4.3
===========================

Improvements
------------
 - On interactive python sessions, group/attribute  `__dir__()` method
   autocompletes children that are named as valid python identifiers.
   :issue:`624` & :issue:`625` thanks to ankostis.
 - Implement `Group.__getitem__()` to have groups act as python-containers,
   so code like this works: ``hfile.root['some child']``.
   :issue:`628` thanks to ankostis.
 - Enable building with Intel compiler (icc/icpc).
   Thanks to rohit-jamuar.
 - PEP 519 support, using new `os.fspath` method.
   Thanks to mruffalo.
 - Optional disable recording of ctime (metadata creation time) when
   creating datasets that makes possible to get bitwise identical output
   from repeated runs.
   Thanks to alex-cobb.
 - Prevent from reading all rows for each coord in a VLArray when
   indexing using a list .
   Thanks to igormq.
 - Internal Blosc version updated to 1.14.3

Bugs fixed
----------
 - Fixed division by zero when using `_convert_time64()` with an empty
   nparr array.
   :issue:`653`. Thanks to alobbs.
 - Fixed deprecation warnings with numpy 1.14.
   Thanks to oleksandr-pavlyk.
 - Skip DLL check when running from a frozen app.
   :issue:`675`. Thanks to jwiggins.
 - Fixed behaviour with slices out of range.
   :issue:`651`. Thanks to jackdbd.


Changes from 3.4.1 to 3.4.2
===========================

Improvements
------------
 - setup.py detects conda env and uses installed conda (hdf5, bzip2, lzo
   and/or blosc) packages when building from source.

Bugs fixed
----------
 - Linux wheels now built against built-in blosc.
 - Fixed windows absolute paths in ptrepack, ptdump, ptree.
   :issue:`616`. Thanks to oscar6echo.


Changes from 3.4.0 to 3.4.1
===========================

Bugs fixed
----------
 - Fixed bug in ptrepack


Changes from 3.3.0 to 3.4.0
===========================

Improvements
------------
 - Support for HDF5 v1.10.x (see :issue:`582`)
 - Fix compatibility with the upcoming Python 2.7.13, 3.5.3 and 3.6 versions.
   See also :issue:`590`. Thanks to Yaroslav Halchenko
 - Internal Blosc version updated to 1.11.3
 - Gracefully handle cpuinfo failure. (PR #578)
   Thanks to Zbigniew Jędrzejewski-Szmek
 - Update internal py-cpuinfo to 3.3.0. Thanks to Gustavo Serra Scalet.

Bugs fixed
----------
 - Fix conversion of python 2 `long` type to `six.integer_types` in atom.py.
   See also :issue:`598`. Thanks to Kyle Keppler for reporting.
 - Fix important bug in bitshuffle filter in internal Blosc on big-endian
   machines. See also :issue:`583`.
 - Fix allow for long type in nextafter. (PR #587) Thanks to Yaroslav Halchenko.
 - Fix unicode bug in group and tables names. :issue:`514`
=======================================
 Release notes for PyTables 2.1 series
=======================================

:Author: Francesc Alted i Abad
:Contact: faltet@pytables.org


Changes from 2.1.1 to 2.1.2
===========================

Bug fixes
---------

- Solved problems with Table.modifyColumn() when the column(s) is
  multidimensional. Fixes #228.

- The row attribute of a table seems stalled after a table move or
  rename.  Fixes #224.

- Fixed a problem with ``len(array)`` in 32-bit platforms when array
  is large enough (> 2**31).

- Added missing `_c_classId` attribute to the `UnImplemented` class.
  ``ptrepack`` no longer chokes while copying `Unimplemented` classes.

- The ``FIELD_*`` sys attrs are no longer copied when the
  ``PYTABLES_SYS_ATTRS`` parameter is set to false.

- The ``FILTERS`` attribute is not added anymore when
  ``PYTABLES_SYS_ATTR`` parameter is set to false.

- Disable the printing of Unicode characters that cannot be printed on
  win32 platform. Fixes #235.

Other changes
-------------

- When retrieving a row of a 1-dimensional array, a 0-dim array was
  returned instead of a numpy scalar.  Now, an actual numpy scalar is
  returned.  Closes #222.

- LZO and bzip2 filters adapted to an API fix introduced in HDF5
  1.8.3.  Closes #225.

- Unsupported HDF5 types in attributes are no longer transferred
  during copies. A new `_v_unimplemented` list have been added in
  `AttributeSet` class so as to keep track of such attributes.  Closes
  #240.

- LZO binaries have disappeared from the GnuWin32 repository.  Until
  they come eventually back, they have been put at
  http://www.pytables.org/download/lzo-win. This has been documented
  in the install chapter.


Changes from 2.1 to 2.1.1
=========================

Bug fixes
---------

- Fixed a memory leak when a lot of queries were made.  Closes #203
  and #207.

- The chunkshape="auto" parameter value of `Leaf.copy()` is honored
  now, even when the (start, stop, step) parameters are specified.
  Closes #204.

- Due to a flaw in its design, the `File` class was not able to be
  subclassed.  This has been fixed.  Closes #205.

- Default values were not correctly retrieved when opening already
  created CArray/EArray objects.  Fixed.  Closes #212.

- Fixed a problem with the installation of the ``nctoh5`` script that
  prevented it from being executed.  Closes #215.

- [Pro] The ``iterseq`` cache ignored non-indexed conditions, giving
  wrong results when those appeared in condition expressions.  This
  has been fixed.  Closes #206.

Other changes
-------------

- `openFile()`, `isHDF5File()` and `isPyTablesFile()` functions accept
  Unicode filenames now.  Closes #202 and #214.

- When creating large type sizes (exceeding 64 KB), HDF5 complained
  and refused to do so.  The HDF5 team has logged the issue as a bug,
  but meanwhile it has been implemented a workaround in PyTables that
  allows to create such large datatypes for situations that does not
  require defaults other than zero.  Addresses #211.

- In order to be consistent with how are stored the other data types,
  Unicode attributes are retrieved now as NumPy scalars instead of
  Python Unicode strings or NumPy arrays.  For the moment, I've fixed
  this through pickling the Unicode strings.  In the future, when HDF5
  1.8.x series would be a requirement, that should be done via a HDF5
  native Unicode type.  Closes #213.


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: rst
.. coding: utf-8
.. fill-column: 72
.. End:
==============================
 What's new in PyTables 1.3.3
==============================


:Author: Francesc Altet
:Contact: faltet@carabos.com
:Author: Ivan Vilata i Balaguer
:Contact: ivilata@carabos.com


This document details the modifications to PyTables since version 1.2.  Its
main purpose is help you ensure that your programs will be runnable when you
switch from PyTables 1.2 to PyTables 1.3.3.


API additions
=============

- None

Backward-incompatible changes
=============================

- None


Deprecated features
===================

- None


API refinements
===============

- None


Bug fixes (affecting API)
=========================

- None


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: text
.. coding: utf-8
.. fill-column: 78
.. End:
==============================
 What's new in PyTables 1.1.1
==============================


:Author: Francesc Altet
:Contact: faltet@carabos.com
:Author: Ivan Vilata i Balaguer
:Contact: ivilata@carabos.com


This document details the modifications to PyTables since version 1.0.
Its main purpose is help you ensure that your programs will be runnable
when you switch from PyTables 1.0 to PyTables 1.1.1.


API additions
=============

- None

Backward-incompatible changes
=============================

- ``Table.read()`` raises a ``KeyError`` instead of a ``ValueError``
  when a nonexistent field name is specified, for consistency with other
  methods.  The same goes for the ``col()`` method.

- ``File.__contains__()`` returns a true value when it is asked for an
  existent node, be it visible or not.  This is more consistent with
  ``Group.__contains__()``.


API refinements
===============

- Using ``table.cols['colname']`` is deprecated.  The usage of
  ``table.cols._f_col('colname')`` (with the new ``Cols._f_col()``
  method) is preferred.

Bug fixes (affecting API)
=========================

- None


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: text
.. coding: utf-8
.. fill-column: 72
.. End:
.. include:: ../../../RELEASE_NOTES.rst
=======================================
 Release notes for PyTables 2.4 series
=======================================

:Author: PyTables maintainers
:Contact: pytables@googlemail.com

.. py:currentmodule:: tables


Changes from 2.3.1 to 2.4
=========================

New features
------------

- Improved HDF5 error logging management:

  * added a new function, :func:`silenceHDF5Messages`, for suppressing
    (and re-enabling) HDF5 messages.  By default HDF5 error logging is now
    suppressed. Closes :issue:`87`.
  * now all HDF5 error messages and trace-backs are trapped and attached to
    the :exc:`exceptions.HDF5ExtError` exception instances.
    Closes :issue:`120`.

- Added support for the float16 data type.  It is only available if numpy_
  provides it as well (i.e. numpy_ >= 1.6).  See :issue:`51`.

- Leaf nodes now have attributes for retrieving the size of data in memory
  and on disk.  Data on disk can be compressed, so the new attributes make it
  easy to compute the data compression ration.
  Thanks to Josh Ayers (close :issue:`141`).

- The maximum number of threads for Blosc_ and Numexpr_ is now handled using
  the :data:`parameters.MAX_BLOSC_THREADS` and
  :data:`parameters.MAX_NUMEXPR_THREADS` parameters respectively.
  This allows a more fine grained configuration capability.
  Closes :issue:`142`.

- `ndim` (read-only) attribute added to :class:`Leaf`, :class:`Atom` and
  :class:`Col` objects (closes :issue:`126`).

- Added read support for variable length string attributes (non scalar
  attributes are converted into numpy_ arrays with 'O8' type).
  See :issue:`54`.


Other improvements
------------------

- Dropped support for HDF5 1.6.x. Now PyTables uses the HDF5 1.8 API
  (closes :issue:`105`).

- Blosc_ updated to v. 1.1.3.

- The Blosc_ compression library is now automatically disabled on platforms
  that do not support unaligned memory access (see also
  https://github.com/FrancescAlted/blosc/issues/3 and
  http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=661286).

- Improved bzip2 detection on Windows (:issue:`116`).  Thanks to cgohlke.

- For Windows, the setup.py script now has the ability to automatically find
  the HDF5_DIR in the system PATH.  Thanks to Mark (mwiebe).

- Improved multi-arch support in GNU/Linux platforms (closes :issue:`124`)
  Thanks to Julian Taylor and Picca Frederic-Emmanuel.

- Use new style syntax for exception raising. Closes :issue:`93`.

- Fixed most of the warnings related to py3k compatibility (see :issue:`92`).

- Fixed pyflakes_ warnings (closes :issue:`102`).

- Cython_ extensions updated to use new constructs (closes :issue:`100`).

- Reduced the number of build warnings (closes :issue:`101`).

- Removed the old lrucache module. It is no more needed after the merge with
  PyTables Pro (closes :issue:`118`).

- Added explicit (import time) testing for hdf5dll.dll on Windows to improve
  diagnostics (closes :issue:`146`).  Thanks to Mark (mwiebe).


Documentation improvements
--------------------------

- new cookbook section (contents have been coming from the PyTables wiki
  on http://www.pytables.org)

- complete rework of the library reference.  Now the entire chapter is
  generated from docstrings using the sphinx autodoc extension.
  A big thank you to Josh Ayers.  Closes :issue:`148`.

- new sphinx theme based on the cloud template


Bugs fixed
----------

- Fixed a segfault on platforms that do not support unaligned memory access
  (closes: :issue:`134`).  Thanks to Julian Taylor.

- Fixed broken inheritance in :class:`IsDescription` classes (thanks to
  Andrea Bedini).  Closes :issue:`65`.

- Fixed table descriptions copy method (closes :issue:`131`).

- Fixed open failures handling (closes :issue:`158`).
  Errors that happen when one tries to open an invalid HDF5 file (e.g. an
  empty file) are now detected earlier by PyTables and a proper exception
  (:exc:`exceptions.HDF5ExtError`) is raised.
  Also, in case of open failures, invalid file descriptors are no more cached.
  Before is fix it was not possible to completely close the bad file and reopen
  the same path, even if a valid file was created in the meanwhile.
  Thanks to Daniele for reporting and for the useful test code.

- Fixed support to rich structured  numpy.dtype in
  :func:`description.descr_from_dtype`.   Closes :issue:`160`.

- Fixed sorting of nested tables that caused AttributeError.
  Closes :issue:`156` and :issue:`157`.  Thanks to Uwe Mayer.

- Fixed flavor deregistration (closes :issue:`163`)


Deprecations
------------

- The :data:`parameters.MAX_THREADS` configuration parameter is now
  deprecated.  Please use :data:`parameters.MAX_BLOSC_THREADS` and
  :data:`parameters.MAX_NUMEXPR_THREADS` instead.
  See :issue:`142`.

- Since the support for HDF5 1.6.x has been dropped, the *warn16incompat*
  argument of the :meth:`File.createExternalLink` method and the
  :exc:`exceptions.Incompat16Warning` exception class are now deprecated.

.. _pyflakes: https://launchpad.net/pyflakes
.. _numpy: http://www.numpy.org
.. _Blosc: https://github.com/FrancescAlted/blosc
.. _Numexpr: http://code.google.com/p/numexpr
.. _Cython: http://www.cython.org


  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: rst
.. coding: utf-8
.. fill-column: 72
.. End:
============================
 What's new in PyTables 1.4
============================


:Author: Francesc Altet
:Contact: faltet@carabos.com
:Author: Ivan Vilata i Balaguer
:Contact: ivilata@carabos.com


This document details the modifications to PyTables since version 1.3.  Its
main purpose is help you ensure that your programs will be runnable when you
switch from PyTables 1.3 to PyTables 1.4.


API additions
=============

- The ``Table.getWhereList()`` method has got a new ``sort`` parameter.  The
  default now is to get the list of parameters unsorted.  Set ``sort`` to True
  to get the old behaviour.  We've done this to avoid unnecessary ordering of
  potentially large sets of coordinates.

- Node creation, copying and moving operations have received a new optional
  `createparents` argument.  When true, the necessary groups in the target
  path that don't exist at the time of running the operation are automatically
  created, so that the target group of the operation always exists.


Backward-incompatible changes
=============================

- None


Deprecated features
===================

- None


API refinements
===============

- ``Description._v_walk()`` has been renamed to ``_f_walk()``, since it is a
  public method, not a value.

- ``Table.removeIndex()`` now accepts a column name in addition to an
  ``Index`` instance (the later is deprecated).  This avoids the user having
  to retrieve the needed ``Index`` object.


Bug fixes (affecting API)
=========================

- None


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: text
.. coding: utf-8
.. fill-column: 78
.. End:
What's new in PyTables 0.8
----------------------------

On this release, many enhancements has been added and some bugs has
been fixed. Here is the (non-exhaustive) list:

- The new VLArray class enables you to store large lists of rows
  containing variable numbers of elements. The elements can
  be scalars or fully multidimensional objects, in the PyTables
  tradition. This class supports two special objects as rows:
  Unicode strings (UTF-8 codification is used internally) and
  generic Python objects (through the use of cPickle).

- The new EArray class allows you to enlarge already existing
  multidimensional homogeneous data objects. Consider it
  an extension of the already existing Array class, but
  with more functionality. Online compression or other filters
  can be applied to EArray instances, for example.

  Another nice feature of EA's is their support for fully
  multidimensional data selection with extended slices.  You
  can write "earray[1,2:3,...,4:200]", for example, to get the
  desired dataset slice from the disk. This is implemented
  using the powerful selection capabilities of the HDF5
  library, which results in very highly efficient I/O
  operations. The same functionality has been added to Array
  objects as well.

- New UnImplemented class. If a dataset contains unsupported
  datatypes, it will be associated with an UnImplemented
  instance, then inserted into to the object tree as usual.
  This allows you to continue to work with supported objects
  while retaining access to attributes of unsupported
  datasets.  This has changed from previous versions, where a
  RuntimeError occurred when an unsupported object was
  encountered.

  The combination of the new UnImplemented class with the
  support for new datatypes will enable PyTables to greatly
  increase the number of types of native HDF5 files that can
  be read and modified.

- Boolean support has been added for all the Leaf objects.

- The Table class has now an append() method that allows you
  to save large buffers of data in one go (i.e. bypassing the
  Row accessor). This can greatly improve data gathering
  speed.

- The standard HDF5 shuffle filter (to further enhance the
      compression level) is supported.

- The standard HDF5 fletcher32 checksum filter is supported.

- As the supported number of filters is growing (and may be
  further increased in the future), a Filters() class has been
  introduced to handle filters more easily.  In order to add
  support for this class, it was necessary to make a change in
  the createTable() method that is not backwards compatible:
  the "compress" and "complib" parameters are deprecated now
  and the "filters" parameter should be used in their
  place. You will be able to continue using the old parameters
  (only a Deprecation warning will be issued) for the next few
  releases, but you should migrate to the new version as soon
  as possible. In general, you can easily migrate old code by
  substituting code in its place::

    table = fileh.createTable(group, 'table', Test, '', complevel, complib)

  should be replaced by::

    table = fileh.createTable(group, 'table', Test, '',
                              Filters(complevel, complib))

- A copy() method that supports slicing and modification of
  filtering capabilities has been added for all the Leaf
  objects. See the User's Manual for more information.

- A couple of new methods, namely copyFile() and copyChilds(),
  have been added to File class, to permit easy replication
  of complete hierarchies or sub-hierarchies, even to
  other files. You can change filters during the copy
  process as well.

- Two new utilities has been added: ptdump and
  ptrepack. The utility ptdump allows the user to examine
  the contents of PyTables files (both metadata and actual
  data). The powerful ptrepack utility lets you
  selectively copy (portions of) hierarchies to specific
  locations in other files. It can be also used as an
  importer for generic HDF5 files.

- The meaning of the stop parameter in read() methods has
  changed. Now a value of 'None' means the last row, and a
  value of 0 (zero) means the first row. This is more
  consistent with the range() function in python and the
  __getitem__() special method in numarray.

- The method Table.removeRows() is no longer limited by table
  size.  You can now delete rows regardless of the size of the
  table.

- The "numarray" value has been added to the flavor parameter
  in the Table.read() method for completeness.

- The attributes (.attr instance variable) are Python
  properties now. Access to their values is no longer
  lazy, i.e. you will be able to see both system or user
  attributes from the command line using the tab-completion
  capability of your python console (if enabled).

- Documentation has been greatly improved to explain all the
  new functionality. In particular, the internal format of
  PyTables is now fully described. You can now build
  "native" PyTables files using any generic HDF5 software
  by just duplicating their format.

- Many new tests have been added, not only to check new
  functionality but also to more stringently check
  existing functionality. There are more than 800 different
  tests now (and the number is increasing :).

- PyTables has a new record in the data size that fits in one
  single file: more than 5 TB (yeah, more than 5000 GB), that
  accounts for 11 GB compressed, has been created on an AMD
  Opteron machine running Linux-64 (the 64 bits version of the
  Linux kernel). See the gory details in:
  http://pytables.sf.net/html/HowFast.html.

- New platforms supported: PyTables has been compiled and tested
  under Linux32 (Intel), Linux64 (AMD Opteron and Alpha), Win32
  (Intel), MacOSX (PowerPC), FreeBSD (Intel), Solaris (6, 7, 8
  and 9 with UltraSparc), IRIX64 (IRIX 6.5 with R12000) and it
  probably works in many more architectures. In particular,
  release 0.8 is the first one that provides a relatively clean
  porting to 64-bit platforms.

- As always, some bugs have been solved (especially bugs that
  occur when deleting and/or overwriting attributes).

- And last, but definitely not least, a new donations section
  has been added to the PyTables web site
  (http://sourceforge.net/projects/pytables, then follow the
  "Donations" tag). If you like PyTables and want this effort
  to continue, please, donate!

Enjoy!,

-- Francesc Alted
falted@pytables.org

============================
 What's new in PyTables 1.2
============================


:Author: Francesc Altet
:Contact: faltet@carabos.com
:Author: Ivan Vilata i Balaguer
:Contact: ivilata@carabos.com


This document details the modifications to PyTables since version 1.1.  Its
main purpose is help you ensure that your programs will be runnable when you
switch from PyTables 1.1 to PyTables 1.2.


API additions
=============

- The user is now allowed to set arbitrary Python (non-persistent) attributes
  on any instance of ``Node``.  If the name matches that of a child node, the
  later will no longer be accessible via natural naming, but it will still be
  available via ``File.getNode()``, ``Group._f_getChild()`` and the group
  children dictionaries.

  Of course, this allows the user to overwrite internal (``^_[cfgv]_``)
  PyTables variables, but this is the way most Python packages work.

- The new ``Group._f_getChild()`` method allows to get a child node (be it
  visible or not) by its name.  This should be more intuitive that using
  ``getattr()`` or using the group children dictionaries.

- The new ``File.isVisibleNode()``, ``Node._f_isVisible()`` and
  ``Leaf.isVisible()`` methods tell whether a node is visible or not, i.e. if
  the node will appear in listing operations such as ``Group._f_listNodes()``.


Backward-incompatible changes
=============================

- ``File.objects``, ``File.groups`` and ``File.leaves`` can no longer be used
  to iterate over all the nodes in the file.  However, they still may be used
  to access any node by its path.

- ``File.__contains__()`` returns a true value when it is asked for an
  existent node, be it visible or not.  This is more consistent with
  ``Group.__contains__()``.

- Using ``Group.__delattr__()`` to remove a child is no longer supported.
  Please use ``Group._f_remove()`` instead.

- The ``indexprops`` attribute is now present on all ``Table`` instances, be
  they indexed or not.  In the last case, it is ``None``.

- Table.getWhereList() now has flavor parameter equal to "NumArray" by
  default, which is more consistent with other methods. Before, flavor
  defaulted to "List".

- The ``extVersion`` variable does no longer exist.  It did not make much
  sense either, since the canonical version of the whole PyTables package is
  that of ``__version__``.

- The ``Row.nrow()`` has been converted into a property, so you have to
  replace any call to ``Row.nrow()`` into ``Row.nrow``.


Deprecated features
===================

- The ``objects``, ``groups`` and ``leaves`` mappings in ``File`` are retained
  only for compatibility purposes.  Using ``File.getNode()`` is recommended to
  access nodes, ``File.__contains__()`` to check for node existence, and
  ``File.walkNodes()`` for iteration purposes.  Using ``isinstance()`` and
  ``*isVisible*()`` methods is the preferred way of checking node type and
  visibility.

  Please note that the aforementioned mappings use the named methods
  internally, so the former have no special performance gains over the later.


API refinements
===============

- The ``isHDF5File()`` and ``isPyTablesFile()`` functions know how to handle
  nonexistent or unreadable files.  An ``IOError`` is raised in those cases.


Bug fixes (affecting API)
=========================

- None


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: text
.. coding: utf-8
.. fill-column: 78
.. End:
=======================================
 Release notes for PyTables 2.2 series
=======================================

:Author: Francesc Alted i Abad
:Contact: faltet@pytables.org


Changes from 2.2.1rc1 to 2.2.1
==============================

- The `Row` accessor implements a new `__contains__` special method that
  allows doing things like::

    for row in table:
        if item in row:
            print "Value found in row", row.nrow
            break

  Closes #309.

- PyTables is more friendly with easy_install and pip now, as all the
  Python dependencies should be installed automatically.  Closes #298.


Changes from 2.2 to 2.2.1rc1
============================

- When using `ObjectAtom` objects in `VLArrays` the ``HIGHEST_PROTOCOL``
  is used for pickling objects.  For NumPy arrays, this simple change
  leads to space savings up to 3x and time improvements up to 30x.
  Closes #301.

- The `Row` accessor implements a new `__contains__` special method that
  allows doing things like::

    for row in table:
        if item in row:
            print "Value found in row", row.nrow
            break

  Closes #309.

- tables.Expr can perform operations on scalars now.  Thanks to Gaëtan
  de Menten for providing a patch for this.  Closes #287.

- Fixed a problem with indexes larger than 32-bit on leaf objects on
  32-bit machines.  Fixes #283.

- Merged in Blosc 1.1.2 for fixing a problem with large datatypes and
  subprocess issues.  Closes #288 and #295.

- Due to the adoption of Blosc 1.1.2, the pthreads-win32 library
  dependency is dropped on Windows platforms.

- Fixed a problem with tables.Expr and operands with vary large
  rowsizes. Closes #300.

- ``leaf[numpy.array[scalar]]`` idiom returns a NumPy array instead of
  an scalar.  This has been done for compatibility with NumPy.  Closes
  #303.

- Optimization for `Table.copy()` so that ``FIELD_*`` attrs are not
  overwritten during the copy.  This can lead to speed-ups up to 100x
  for short tables that have hundreds of columns.  Closes #304.

- For external links, its relative paths are resolved now with respect
  to the directory of the main HDF5 file, rather than with respect to
  the current directory.  Closes #306.

- ``Expr.setInputsRange()`` and ``Expr.setOutputRange()`` do support
  ``numpy.integer`` types now.  Closes #285.

- Column names in tables can start with '__' now.  Closes #291.

- Unicode empty strings are supported now as attributes.  Addresses #307.

- Cython 0.13 and higher is supported now.  Fixes #293.

- PyTables should be more 'easy_install'-able now.  Addresses #298.


Changes from 2.2rc2 to 2.2 (final)
==================================

- Updated Blosc to 1.0 (final).

- Filter ID of Blosc changed from wrong 32010 to reserved 32001.  This
  will prevent PyTables 2.2 (final) to read files created with Blosc and
  PyTables 2.2 pre-final.  `ptrepack` can be used to retrieve those
  files, if necessary.  More info in ticket #281.

- Recent benchmarks suggest a new parametrization is better in most
  scenarios:

  * The default chunksize has been doubled for every dataset size.  This
    works better in most of scenarios, specially with the new Blosc
    compressor.

  * The HDF5 CHUNK_CACHE_SIZE parameter has been raised to 2 MB in order
    to better adapt to the chunksize increase.  This provides better hit
    ratio (at the cost of consuming more memory).

  Some plots have been added to the User's Manual (chapter 5) showing
  how the new parametrization works.


Changes from 2.2rc1 to 2.2rc2
=============================

- A new version of Blosc (0.9.5) is included.  This version is now
  considered to be stable and apt for production.  Thanks for all
  PyTables users that have contributed to find and report bugs.

- Added a new `IO_BUFFER_SIZE` parameter to ``tables/parameters.py``
  that allows to set the internal PyTables' buffer for doing I/O.  This
  replaces `CHUNKTIMES` but it is more general because it affects to all
  `Leaf` objects and also the `tables.Expr` module (and not only tables
  as before).

- `BUFFERTIMES` parameter in ``tables/parameters.py`` has been
  renamed to `BUFFER_TIMES` which is more consistent with other
  parameter names.

- On Windows platforms, the path to the tables module is now appended to
  sys.path and the PATH environment variable. That way DLLs and PYDs in
  the tables directory are to be found now.  Thanks to Christoph Gohlke
  for the hint.

- A replacement for barriers for Mac OSX, or other systems not
  implementing them, has been carried out.  This allows to compile
  PyTables on such platforms.  Fixes #278

- Fixed a couple of warts that raise compatibility warnings with
  forthcoming Python 2.7.

-  HDF5 1.8.5 is used in Windows binaries.

Changes from 2.2b3 to 2.2rc1
============================

- Numexpr is not included anymore in PyTables and has become a requisite
  instead.  This is because Numexpr already has decent enough installers
  and is available in the PyPI repository also, so it should be easy for
  users to fulfill this dependency.

- When using a Numexpr package that is turbo-loaded with Intel's
  VML/MKL, the parameter `MAX_THREADS` will control the number of
  threads that VML can use during computations.  For a finer control,
  the `numexpr.set_vml_num_threads()` can always be used.

- Cython is used now instead of Pyrex for Pyrex extensions.

- Updated to 0.9 version of Blosc compressor.  This version can make use
  of threads so as to accelerate the compression/decompression process.
  In order to change the maximum number of threads that Blosc can use (2
  by default), you can modify the `MAX_THREADS` variable in
  ``tables/parameters.py`` or make use of the new `setBloscMaxThreads()`
  global function.

- Reopening already opened files is supported now, provided that there is
  not incompatibility among intended usages (for example, you cannot
  reopen in append mode an already opened file in read-only mode).

- Option ``--print-versions`` for ``test_all.py`` script is now
  preferred over the deprecated ``--show-versions``.  This is more
  consistent with the existing `print_versions()` function.

- Fixed a bug that, under some circumstances, prevented the use of table
  iterators in `itertool.groupby()`.  Now, you can safely do things
  like::

    sel_rows = table.where('(row_id >= 3)')
    for group_id, grouped_rows in itertools.groupby(sel_rows, f_group):
        group_mean = average([row['row_id'] for row in grouped_rows])

  Fixes #264.

- Copies of `Array` objects with multidimensional atoms (coming from
  native HDF5 files) work correctly now (i.e. the copy holds the atom
  dimensionality).  Fixes #275.

- The `tables.openFile()` function does not try anymore to open/close
  the file in order to guess whether it is a HDF5 or PyTables one before
  opening it definitely.  This allows the `fcntl.flock()` and
  `fcntl.lockf()` Python functions to work correctly now (that's useful
  for arbitrating access to the file by different processes).  Thanks to
  Dag Sverre Seljebotn and Ivan Vilata for their suggestions on hunting
  this one!  Fixes #185.

- The estimation of the chunksize when using multidimensional atoms in
  EArray/Carray was wrong because it did not take in account the shape
  of the atom.  Thanks to Ralf Juengling for reporting.  Fixes #273.

- Non-contiguous arrays can now safely be saved as attributes.  Before,
  if arrays were not contiguous, incorrect data was saved in attr.
  Fixes #270.

- EXTDIM attribute for CArray/EArray now saves the correct extendable
  dimension, instead of rubbish.  This does not affected functionality,
  because extendable dimension was retrieved directly from shape
  information, but it was providing misleading information to the user.
  Fixes #268.

API changes
-----------

- Now, `Table.Cols.__len__()` returns the number of top level columns
  instead of the number of rows in table.  This is more consistent in
  that `Table.Cols` is an accessor for *columns*.  Fixes #276.


Changes from 2.2b2 to 2.2b3
===========================

- Blosc compressor has been added as an additional filter, in addition
  to the existing Zlib, LZO and bzip2.  This new compressor is meant for
  fast compression and extremely fast decompression.  Fixes #265.

- In `File.copyFile()` method, `copyuserattrs` was set to false as
  default.  This was inconsistent with other methods where the default
  value for `copyuserattrs` is true.  The default for this is true now.
  Closes #261.

- `tables.copyFile` and `File.copyFile` recognize now the parameters
  present in ``tables/parameters.py``.  Fixes #262.

- Backported fix for issue #25 in Numexpr (OP_NEG_LL treats the argument
  as an int, not a long long).  Thanks to David Cooke for this.

- CHUNK_CACHE_NELMTS in `tables/parameters.py` set to a prime number as
  Neil Fortner suggested.

- Workaround for a problem in Python 2.6.4 (and probably other versions
  too) for pickling strings like "0" or "0.".  Fixes #253.


Changes from 2.2b1 to 2.2b2
===========================

Enhancements
------------

- Support for HDF5 hard links, soft links and external links (when
  PyTables is compiled against HDF5 1.8.x series).  A new tutorial about
  its usage has been added to the 'Tutorials' chapter of User's Manual.
  Closes #239 and #247.

- Added support for setting HDF5 chunk cache parameters in file
  opening/creating time.  'CHUNK_CACHE_NELMTS', 'CHUNK_CACHE_PREEMPT'
  and 'CHUNK_CACHE_SIZE' are the new parameters.  See "PyTables'
  parameter files" appendix in User's Manual for more info.  Closes
  #221.

- New `Unknown` class added so that objects that HDF5 identifies as
  ``H5G_UNKNOWN`` can be mapped to it and continue operations
  gracefully.

- Optimization in the indexed queries when the resulting rows increase
  monotonically.  From 3x (for medium-size query results) and 10x (for very
  large query results) speed-ups can be expected.

- Added flag `--dont-create-sysattrs` to ``ptrepack`` so as to not
  create sys attrs (default is to do it).

- Support for native compound types in attributes.  This allows for
  better compatibility with HDF5 files.  Closes #208.

- Support for native NumPy dtype in the description parameter of
  `File.createTable()`.  Closes #238.


Bugs fixed
----------

- Added missing `_c_classId` attribute to the `UnImplemented` class.
  ``ptrepack`` no longer chokes while copying `Unimplemented` classes.

- The ``FIELD_*`` sys attrs are no longer copied when the
  ``PYTABLES_SYS_ATTRS`` parameter is set to false.

- `File.createTable()` no longer segfaults if description=None.  Closes
  #248.

- Workaround for avoiding a Python issue causing a segfault when saving
  and then retrieving a string attribute with values "0" or "0.".
  Closes #253.


API changes
-----------

- `Row.__contains__()` disabled because it has little sense to query for
  a key in Row, and the correct way should be to query for it in
  `Table.colnames` or `Table.colpathnames` better.  Closes #241.

- [Semantic change] To avoid a common pitfall when asking for the string
  representation of a `Row` class, `Row.__str__()` has been redefined.
  Now, it prints something like::

      >>> for row in table:
      ...     print row
      ...
      /newgroup/table.row (Row), pointing to row #0
      /newgroup/table.row (Row), pointing to row #1
      /newgroup/table.row (Row), pointing to row #2

  instead of::

      >>> for row in table:
      ...     print row
      ...
      ('Particle:      0', 0, 10, 0.0, 0.0)
      ('Particle:      1', 1, 9, 1.0, 1.0)
      ('Particle:      2', 2, 8, 4.0, 4.0)

  Use `print row[:]` idiom if you want to reproduce the old behaviour.
  Closes #252.


Other changes
-------------

- After some improvements in both HDF5 and PyTables, the limit before
  emitting a `PerformanceWarning` on the number of children in a group
  has been raised from 4096 to 16384.


Changes from 2.1.1 to 2.2b1
===========================

Enhancements
------------

- Added `Expr`, a class for evaluating expressions containing
  array-like objects.  It can evaluate expressions (like '3*a+4*b')
  that operate on arbitrary large arrays while optimizing the
  resources (basically main memory and CPU cache memory) required to
  perform them.  It is similar to the Numexpr package, but in addition
  to NumPy objects, it also accepts disk-based homogeneous arrays,
  like the `Array`, `CArray`, `EArray` and `Column` PyTables objects.

- Added support for NumPy's extended slicing in all `Leaf` objects.
  With that, you can do the next sort of selections::

      array1 = array[4]                       # simple selection
      array2 = array[4:1000:2]                # slice selection
      array3 = array[1, ..., ::2, 1:4, 4:]    # general slice selection
      array4 = array[1, [1,5,10], ..., -1]    # fancy selection
      array5 = array[np.where(array[:] > 4)]  # point selection
      array6 = array[array[:] > 4]            # boolean selection

  Thanks to Andrew Collette for implementing this for h5py, from which
  it has been backported.  Closes #198 and #209.

- Numexpr updated to 1.3.1.  This can lead to up a 25% improvement of
  the time for both in-kernel and indexed queries for unaligned
  tables.

- HDF5 1.8.3 supported.


Bugs fixed
----------

- Fixed problems when modifying multidimensional columns in Table
  objects.  Closes #228.

- Row attribute is no longer stalled after a table move or rename.
  Fixes #224.

- Array.__getitem__(scalar) returns a NumPy scalar now, instead of a
  0-dim NumPy array.  This should not be noticed by normal users,
  unless they check for the type of returned value.  Fixes #222.


API changes
-----------

- Added a `dtype` attribute for all leaves.  This is the NumPy
  ``dtype`` that most closely matches the leaf type.  This allows for
  a quick-and-dirty check of leaf types.  Closes #230.

- Added a `shape` attribute for `Column` objects.  This is formed by
  concatenating the length of the column and the shape of its type.
  Also, the representation of columns has changed an now includes the
  length of the column as the leading dimension.  Closes #231.

- Added a new `maindim` attribute for `Column` which has the 0 value
  (the leading dimension).  This allows for a better similarity with
  other \*Array objects.

- In order to be consistent and allow the extended slicing to happen
  in `VLArray` objects too, `VLArray.__setitem__()` is not able to
  partially modify rows based on the second dimension passed as key.
  If this is tried, an `IndexError` is raised now.  Closes #210.

- The `forceCSI` flag has been replaced by `checkCSI` in the next
  `Table` methods: `copy()`, `readSorted()` and `itersorted()`.  The
  change reflects the fact that a re-index operation cannot be
  triggered from these methods anymore.  The rational for the change
  is that an indexing operation is a potentially very expensive
  operation that should be carried out explicitly instead of being
  triggered by methods that should not be in charge of this task.
  Closes #216.


Backward incompatible changes
-----------------------------

- After the introduction of the `shape` attribute for `Column`
  objects, the shape information for multidimensional columns has been
  removed from the `dtype` attribute (it is set to the base type of
  the column now).  Closes #232.


  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: rst
.. coding: utf-8
.. fill-column: 72
.. End:
=======================================
 Release notes for PyTables 3.3 series
=======================================

Changes from 3.2.3.1 to 3.3.0
=============================

Improvements
------------

- Single codebase Python 2 and 3 support (PR #493).
- Internal Blosc version updated to 1.11.1 (closes :issue:`541`)
- Full BitShuffle support for new Blosc versions (>= 1.8).
- It is now possible to remove all rows from a table.
- It is now possible to read reference types by dereferencing them as
  numpy array of objects (closes :issue:`518` and :issue:`519`).
  Thanks to Ehsan Azar
- Get rid of the `-native` compile flag (closes :issue:`503`)
- The default number of threads to run numexpr (MAX_NUMEXPR_THREADS)
  internally has been raised from 2 to 4.  This is because we are in
  2016 and 4 core configurations are becoming common.
- In order to avoid locking issues when using PyTables concurrently in
  several process, MAX_BLOSC_THREADS has been set to 1 by default.  If
  you are running PyTables in one single process, you may want to
  experiment if higher values (like 2 or 4) bring better performance for
  you.

Bugs fixed
----------

- On python 3 try 'latin1' encoding before 'bytes' encoding during unpickling
  of node attributes pickled on python 2. Better fix for :issue:`560`.
- Fixed Windows 32 and 64-bit builds.
What's new in PyTables 0.9.1
----------------------------

This release is mainly a maintenance version. In it, some bugs has
been fixed and a few improvements has been made. One important thing
is that chunk sizes in EArrays has been re-tuned to get much better
performance. Besides, it has been tested against the latest Python 2.4
and all unit tests seems to pass fine.

More in detail:

Improvements:

- The chunksize computation for EArrays has been re-tuned to allow the
  compression rations that were usual before 0.9 release.

- New --unpackshort and --quantize flags has been added to nctoh5
  script. --unpackshort unpack short integer variables to float
  variables using scale_factor and add_offset netCDF variable
  attributes. --quantize quantize data to improve compression using
  least_significant_digit netCDF variable attribute (not active by
  default).  See https://www.ogc.org/standards/netcdf
  for further explanation of what this attribute means. Thanks to Jeff
  Whitaker for providing this.

- Table.itersequence has received a new parameter called "sort". This
  allows to disable the sorting of the sequence in case the user wants
  so.

Backward-incompatible changes:

- Now, the AttributeSet class throw an AttributeError on __getattr__
  for nonexistent attributes in it. Formerly, the routine returned
  None, which is pretty much against convention in Python and breaks
  the built-in hasattr() function. Thanks to Robert Nemec for noting
  this and offering a patch.

- VLArray.read() has changed its behaviour. Now, it always returns a
  list, as stated in documentation, even when the number of elements
  to return is 0 or 1. This is much more consistent when representing
  the actual number of elements on a certain VLArray row.

API additions:

- A Row.getTable() has been added. It is an accessor for the associated
  Table object.

- A File.copyAttrs() has been added. It allows copying attributes from
  one leaf to other. Properly speaking, this was already there, but not
  documented :-/

Bug fixes:

- Now, the copy of hierarchies works even when there are scalar Arrays
  (i.e. Arrays which shape is ()) on it. Thanks to Robert Nemec for
  providing a patch.

- Solved a memory leak regarding the Filters instance associated with
  the File object, that was not released after closing the file. Now,
  there are no known leaks on PyTables itself.

- Improved security of nodes name checking. Closes #1074335


Enjoy data!,

-- Francesc Altet
falted@pytables.org

=======================================
 Release notes for PyTables 3.6 series
=======================================

:Author: PyTables Developers
:Contact: pytables-dev@googlegroups.com

.. py:currentmodule:: tables


Changes from 3.6.0 to 3.6.1
===========================

Maintenance release to fix packaging issues. No new features or bugfixes.


Changes from 3.5.3 to 3.6.0
===========================

PyTables 3.6 no longer supports Python 2.7 see PR #747.

Improvements
------------
- Full python 3.8 support.
- On Windows PyTables wheels on PyPI are linked to `pytables_hdf5.dll` instead
  of `hdf5.dll` to prevent collisions with other packages/wheels that also
  vendor `hdf5.dll`.
  This should prevent problems that arise when a different version of a dll
  is imported that the version to which the program was linked to.
  This problem is known as "DLL Hell".
  With the renaming of the HDF5 DLL to `pytables_hdf5.dll` these problems
  should be solved.

Bugfixes
--------
- Bugfix for HDF5 files/types with padding. For details see :issue:`734`.
- More fixes for python 3.8 compatibility: Replace deprecated time.clock
  with time.perf_counter
  Thanks to Sergio Pascual (sergiopasra). see :issue:`744` and PR #745.
- Improvements in tests as well as clean up from dropping Python 2.7 support.
  Thanks to Seth Troisi (sethtroisi).
==============================
 What's new in PyTables 1.2.2
==============================


:Author: Francesc Altet
:Contact: faltet@carabos.com
:Author: Ivan Vilata i Balaguer
:Contact: ivilata@carabos.com


This document details the modifications to PyTables since version 1.2.  Its
main purpose is help you ensure that your programs will be runnable when you
switch from PyTables 1.2 to PyTables 1.2.2.


API additions
=============

- None

Backward-incompatible changes
=============================

- None

Deprecated features
===================

- None


API refinements
===============

- None


Bug fixes (affecting API)
=========================

- None


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: text
.. coding: utf-8
.. fill-column: 78
.. End:
==============================
 What's new in PyTables 1.3.1
==============================


:Author: Francesc Altet
:Contact: faltet@carabos.com
:Author: Ivan Vilata i Balaguer
:Contact: ivilata@carabos.com


This document details the modifications to PyTables since version 1.2.  Its
main purpose is help you ensure that your programs will be runnable when you
switch from PyTables 1.2 to PyTables 1.3.1.


API additions
=============

- The Table.Cols accessor has received a new __setitem__() method that
  allows doing things like:

            table.cols[4] = record
            table.cols.x[4:1000:2] = array   # homogeneous column
            table.cols.Info[4:1000:2] = recarray   # nested column


Backward-incompatible changes
=============================

- None


Deprecated features
===================

- None


API refinements
===============

- Table.itersequence has changed the default value for 'sort' parameter. It is
  now False by default, as it is not clear if this actually accelerates the
  iterator, so it is better to let to the user doing the proper checks (if he
  is interested at all).


Bug fixes (affecting API)
=========================

- None


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: text
.. coding: utf-8
.. fill-column: 78
.. End:
=======================================
 Release notes for PyTables 3.0 series
=======================================

:Author: PyTables Developers
:Contact: pytables@googlemail.com

.. py:currentmodule:: tables


Changes from 2.4 to 3.0
=======================

New features
------------

- Since this release PyTables provides full support to Python_ 3
  (closes :issue:`188`).

- The entire code base is now more compliant with coding style guidelines
  describe in the PEP8_ (closes :issue:`103` and :issue:`224`).
  See `API changes`_ for more details.

- Basic support for HDF5 drivers.  Now it is possible to open/create an
  HDF5 file using one of the SEC2, DIRECT, LOG, WINDOWS, STDIO or CORE
  drivers.  Users can also set the main driver parameters (closes
  :issue:`166`).
  Thanks to Michal Slonina.

- Basic support for in-memory image files.  An HDF5 file can be set from or
  copied into a memory buffer (thanks to Michal Slonina).  This feature is
  only available if PyTables is built against HDF5 1.8.9 or newer.
  Closes :issue:`165` and :issue:`173`.

- New :meth:`File.get_filesize` method for retrieving the HDF5 file size.

- Implemented methods to get/set the user block size in a HDF5 file
  (closes :issue:`123`)

- Improved support for PyInstaller_.  Now it is easier to pack frozen
  applications that use the PyTables package (closes: :issue:`177`).
  Thanks to Stuart Mentzer and Christoph Gohlke.

- All read methods now have an optional *out* argument that allows to pass a
  pre-allocated array to store data (closes :issue:`192`)

- Added support for the floating point data types with extended precision
  (Float96, Float128, Complex192 and Complex256).  This feature is only
  available if numpy_ provides it as well.
  Closes :issue:`51` and :issue:`214`.  Many thanks to Andrea Bedini.

- Consistent ``create_xxx()`` signatures.  Now it is possible to create all
  data sets :class:`Array`, :class:`CArray`, :class:`EArray`,
  :class:`VLArray`, and :class:`Table` from existing Python objects (closes
  :issue:`61` and :issue:`249`).  See also the `API changes`_ section.

- Complete rewrite of the :mod:`nodes.filenode` module. Now it is fully
  compliant with the interfaces defined in the standard :mod:`io` module.
  Only non-buffered binary I/O is supported currently.
  See also the `API changes`_ section.  Closes :issue:`244`.

- New :program:`pt2to3` tool is provided to help users to port their
  applications to the new API (see `API changes`_ section).


Improvements
------------

- Improved runtime checks on dynamic loading of libraries: meaningful error
  messages are generated in case of failure.
  Also, now PyTables no more alters the system PATH.
  Closes :issue:`178` and :issue:`179` (thanks to Christoph Gohlke).

- Improved list of search paths for libraries as suggested by Nicholaus
  Halecky (see :issue:`219`).

- Removed deprecated Cython_ include (.pxi) files. Contents of
  :file:`convtypetables.pxi` have been moved in :file:`utilsextension.pyx`.
  Closes :issue:`217`.

- The internal Blosc_ library has been upgraded to version 1.2.3.

- Pre-load the bzip2_ library on windows (closes :issue:`205`)

- The :meth:`File.get_node` method now accepts unicode paths
  (closes :issue:`203`)

- Improved compatibility with Cython_ 0.19 (see :issue:`220` and
  :issue:`221`)

- Improved compatibility with numexpr_ 2.1 (see also :issue:`199` and
  :issue:`241`)

- Improved compatibility with development versions of numpy_
  (see :issue:`193`)

- Packaging: since this release the standard tar-ball package no more includes
  the PDF version of the "PyTables User Guide", so it is a little bit smaller
  now.  The complete and pre-build version of the documentation both in HTML
  and PDF format is available on the file `download area`_ on SourceForge.net.
  Closes: :issue:`172`.

- Now PyTables also uses `Travis-CI`_ as continuous integration service.
  All branches and all pull requests are automatically tested with different
  Python_ versions.  Closes :issue:`212`.


Other changes
-------------

- PyTables now requires Python 2.6 or newer.

- Minimum supported version of Numexpr_ is now 2.0.


API changes
-----------

The entire PyTables API as been made more PEP8_ compliant (see :issue:`224`).

This means that many methods, attributes, module global variables and also
keyword parameters have been renamed to be compliant with PEP8_ style
guidelines (e.g. the ``tables.hdf5Version`` constant has been renamed into
``tables.hdf5_version``).

We made the best effort to maintain compatibility to the old API for existing
applications.  In most cases, the old 2.x API is still available and usable
even if it is now deprecated (see the Deprecations_ section).

The only important backwards incompatible API changes are for names of
function/methods arguments.  All uses of keyword arguments should be
checked and fixed to use the new naming convention.

The new :program:`pt2to3` tool can be used to port PyTables based applications
to the new API.

Many deprecated features and support for obsolete modules has been dropped:

- The deprecated :data:`is_pro` module constant has been removed

- The nra module and support for the obsolete numarray module has been removed.
  The *numarray* flavor is no more supported as well (closes :issue:`107`).

- Support for the obsolete Numeric module has been removed.
  The *numeric* flavor is no longer available (closes :issue:`108`).

- The tables.netcdf3 module has been removed (closes :issue:`68`).

- The deprecated :exc:`exceptions.Incompat16Warning` exception has been
  removed

- The :meth:`File.create_external_link` method no longer has a keyword
  parameter named *warn16incompat*.  It was deprecated in PyTables 2.4.

Moreover:

- The :meth:`File.create_array`, :meth:`File.create_carray`,
  :meth:`File.create_earray`, :meth:`File.create_vlarray`, and
  :meth:`File.create_table` methods of the :class:`File` objects gained a
  new (optional) keyword argument named ``obj``.  It can be used to initialize
  the newly created dataset with an existing Python object, though normally
  these are numpy_ arrays.

  The *atom*/*descriptor* and *shape* parameters are now optional if the
  *obj* argument is provided.

- The :mod:`nodes.filenode` has been completely rewritten to be fully
  compliant with the interfaces defined in the :mod:`io` module.

  The FileNode classes currently implemented are intended for binary I/O.

  Main changes:

  * the FileNode base class is no more available,
  * the new version of :class:`nodes.filenode.ROFileNode` and
    :class:`nodes.filenode.RAFileNode` objects no more expose the *offset*
    attribute (the *seek* and *tell* methods can be used instead),
  * the *lineSeparator* property is no more available and the ``\n``
    character is always used as line separator.

- The `__version__` module constants has been removed from almost all the
  modules (it was not used after the switch to Git).  Of course the package
  level constant (:data:`tables.__version__`) still remains.
  Closes :issue:`112`.

- The :func:`lrange` has been dropped in favor of xrange (:issue:`181`)

- The :data:`parameters.MAX_THREADS` configuration parameter has been dropped
  in favor of :data:`parameters.MAX_BLOSC_THREADS` and
  :data:`parameters.MAX_NUMEXPR_THREADS` (closes :issue:`147`).

- The :func:`conditions.compile_condition` function no more has a *copycols*
  argument, it was no more necessary since Numexpr_ 1.3.1.
  Closes :issue:`117`.

- The *expectedsizeinMB* parameter of the :meth:`File.create_vlarray` and of
  the :meth:`VLArrsy.__init__` methods has been replaced by *expectedrows*.
  See also (:issue:`35`).

- The :meth:`Table.whereAppend` method has been renamed into
  :meth:`Table.append_where` (closes :issue:`248`).

Please refer to the :doc:`../MIGRATING_TO_3.x` document for more details about
API changes and for some useful hint about the migration process from the 2.X
API to the new one.


Other possibly incompatible changes
-----------------------------------

- All methods of the :class:`Table` class that take *start*, *stop* and
  *step* parameters (including :meth:`Table.read`, :meth:`Table.where`,
  :meth:`Table.iterrows`, etc) have been redesigned to have a consistent
  behaviour.  The meaning of the *start*, *stop* and *step* and their default
  values now always work exactly like in the standard :class:`slice` objects.
  Closes :issue:`44` and :issue:`255`.

- Unicode attributes are not stored in the HDF5 file as pickled string.
  They are now saved on the HDF5 file as UTF-8 encoded strings.

  Although this does not introduce any API breakage, files produced are
  different (for unicode attributes) from the ones produced by earlier
  versions of PyTables.

- System attributes are now stored in the HDF5 file using the character set
  that reflects the native string behaviour: ASCII for Python 2 and UTF8 for
  Python 3.  In any case, system attributes are represented as Python string.

- The :meth:`iterrows` method of :class:`*Array` and :class:`Table` as well
  as the :meth:`Table.itersorted` now behave like functions in the standard
  :mod:`itertools` module.
  If the *start* parameter is provided and *stop* is None then the
  array/table is iterated from *start* to the last line.
  In PyTables < 3.0 only one element was returned.


Deprecations
------------

- As described in `API changes`_, all functions, methods and attribute names
  that was not compliant with the PEP8_ guidelines have been changed.
  Old names are still available but they are deprecated.

- The use of upper-case keyword arguments in the :func:`open_file` function
  and the :class:`File` class initializer is now deprecated.  All parameters
  defined in the :file:`tables/parameters.py` module can still be passed as
  keyword argument to the :func:`open_file` function just using a lower-case
  version of the parameter name.


Bugs fixed
----------

- Better check access on closed files (closes :issue:`62`)

- Fix for :meth:`File.renameNode` where in certain cases
  :meth:`File._g_updateLocation` was wrongly called (closes :issue:`208`).
  Thanks to Michka Popoff.

- Fixed ptdump failure on data with nested columns (closes :issue:`213`).
  Thanks to Alexander Ford.

- Fixed an error in :func:`open_file` when *filename* is a :class:`numpy.str_`
  (closes :issue:`204`)

- Fixed :issue:`119`, :issue:`230` and :issue:`232`, where an index on
  :class:`Time64Col` (only, :class:`Time32Col` was ok) hides the data on
  selection from a Tables. Thanks to Jeff Reback.

- Fixed ``tables.tests.test_nestedtypes.ColsTestCase.test_00a_repr`` test
  method.  Now the ``repr`` of cols on big-endian platforms is correctly
  handled  (closes :issue:`237`).

- Fixes bug with completely sorted indexes where *nrowsinbuf* must be equal
  to or greater than the *chunksize* (thanks to Thadeus Burgess).
  Closes :issue:`206` and :issue:`238`.

- Fixed an issue of the :meth:`Table.itersorted` with reverse iteration
  (closes :issue:`252` and :issue:`253`).


.. _Python: http://www.python.org
.. _PEP8: http://www.python.org/dev/peps/pep-0008
.. _PyInstaller: http://www.pyinstaller.org
.. _Blosc: https://github.com/FrancescAlted/blosc
.. _bzip2: http://www.bzip.org
.. _Cython: http://www.cython.org
.. _Numexpr: http://code.google.com/p/numexpr
.. _numpy: http://www.numpy.org
.. _`download area`: http://sourceforge.net/projects/pytables/files/pytables
.. _`Travis-CI`: https://travis-ci.org


  **Enjoy data!**

  -- The PyTables Developers


.. Local Variables:
.. mode: rst
.. coding: utf-8
.. fill-column: 72
.. End:
Changes from 3.1.0 to 3.1.1
===========================

Bugs fixed
----------

- Fixed a critical bug that caused an exception at import time.
  The error was triggered when a bug in long-double detection is detected
  in the HDF5 library (see :issue:`275`) and numpy_ does not expose
  `float96` or `float128`. Closes :issue:`344`.
- The internal Blosc_ library has been updated to version 1.3.5.
  This fixes a false buffer overrun condition that made c-blosc to fail,
  even if the problem was not real.


Improvements
------------

- Do not create a temporary array when the *obj* parameter is not specified
  in :meth:`File.create_array` (thanks to Francesc).
  Closes :issue:`337` and :issue:`339`).
- Added two new utility functions
  (:func:`tables.nodes.filenode.read_from_filenode` and
  :func:`tables.nodes.filenode.save_to_filenode`) for the direct copy from
  filesystem to filenode and vice versa (closes :issue:`342`).
  Thanks to Andreas Hilboll.
- Removed the :file:`examples/nested-iter.py` considered no longer useful.
  Closes :issue:`343`.
- Better detection of the `-msse2` compiler flag.


Changes from 3.0 to 3.1.0
=========================

New features
------------

- Now PyTables is able to save/restore the default value of :class:`EnumAtom`
  types (closes :issue:`234`).
- Implemented support for the H5FD_SPLIT driver (closes :issue:`288`,
  :issue:`289` and :issue:`295`). Many thanks to simleo.
- New quantization filter: the filter truncates floating point data to a
  specified precision before writing to disk. This can significantly improve
  the performance of compressors (closes :issue:`261`).
  Thanks to Andreas Hilboll.
- Added new :meth:`VLArray.get_row_size` method to :class:`VLArray` for
  querying the number of atoms of a :class:`VLArray` row.
  Closes :issue:`24` and :issue:`315`.
- The internal Blosc_ library has been updated to version 1.3.2.
  All new features introduced in the Blosc_ 1.3.x series, and in particular
  the ability to leverage different compressors within Blosc_ (see the `Blosc
  Release Notes`_), are now available in PyTables via the blosc filter
  (closes: :issue:`324`). A big thank you to Francesc.


Improvements
------------

- The node caching mechanism has been completely redesigned to be simpler and
  less dependent from specific behaviours of the ``__del__`` method.
  Now PyTables is compatible with the forthcoming Python 3.4.
  Closes :issue:`306`.
- PyTables no longer uses shared/cached file handlers. This change somewhat
  improves support for concurrent reading allowing the user to safely open the
  same file in different threads for reading (requires HDF5 >= 1.8.7).
  More details about this change can be found in the `Backward incompatible
  changes`_ section.
  See also :issue:`130`, :issue:`129` :issue:`292` and :issue:`216`.
- PyTables is now able to detect and use external installations of the Blosc_
  library (closes :issue:`104`).  If Blosc_ is not found in the system, and the
  user do not specify a custom installation directory, then it is used an internal
  copy of the Blosc_ source code.
- Automatically disable extended float support if a buggy version of HDF5
  is detected (see also `Issues with H5T_NATIVE_LDOUBLE`_).
  See also :issue:`275`, :issue:`290` and :issue:`300`.
- Documented an unexpected behaviour with string literals in query conditions
  on Python 3 (closes :issue:`265`)
- The deprecated :mod:`getopt` module has been dropped in favour of
  :mod:`argparse` in all command line utilities (close :issue:`251`)
- Improved the installation section of the :doc:`../usersguide/index`.

  * instructions for installing PyTables via pip_ have been added.
  * added a reference to the Anaconda_, Canopy_ and `Christoph Gohlke suites`_
    (closes :issue:`291`)

- Enabled `Travis-CI`_ builds for Python_ 3.3
- :meth:`Tables.read_coordinates` now also works with boolean indices input.
  Closes :issue:`287` and :issue:`298`.
- Improved compatibility with numpy_ >= 1.8 (see :issue:`259`)
- The code of the benchmark programs (bench directory) has been updated.
  Closes :issue:`114`.
- Fixed some warning related to non-unicode file names (the Windows bytes API
  has been deprecated in Python 3.4)


Bugs fixed
----------

- Fixed detection of platforms supporting Blosc_
- Fixed a crash that occurred when one attempts to write a numpy_ array to
  an :class:`Atom` (closes :issue:`209` and :issue:`296`)
- Prevent creation of a table with no columns (closes :issue:`18` and
  :issue:`299`)
- Fixed a memory leak that occured when iterating over
  :class:`CArray`/:class:`EArray` objects (closes :issue:`308`,
  see also :issue:`309`).
  Many thanks to Alistair Muldal.
- Make NaN types sort to the end. Closes :issue:`282` and :issue:`313`
- Fixed selection on float columns when NaNs are present (closes :issue:`327`
  and :issue:`330`)
- Fix computation of the buffer size for iterations on rows.
  The buffers size was overestimated resulting in a :exc:`MemoryError`
  in some cases.
  Closes :issue:`316`. Thamks to bbudescu.
- Better check of file open mode. Closes :issue:`318`.
- The Blosc filter now works correctly together with fletcher32.
  Closes :issue:`21`.
- Close the file handle before trying to delete the corresponding file.
  Fixes a test failure on Windows.
- Use integer division for computing indices (fixes some warning on Windows)


Deprecations
------------

Following the plan for the complete transition to the new (PEP8_ compliant)
API, all calls to the old API will raise a :exc:`DeprecationWarning`.

The new API has been introduced in PyTables 3.0 and is backward incompatible.
In order to guarantee a smoother transition the old API is still usable even
if it is now deprecated.

The plan for the complete transition to the new API is outlined in
:issue:`224`.


Backward incompatible changes
-----------------------------

In PyTables <= 3.0 file handles (objects that are returned by the
:func:`open_file` function) were stored in an internal registry and re-used
when possible.

Two subsequent attempts to open the same file (with compatible open mode)
returned the same file handle in PyTables <= 3.0::

    In [1]: import tables
    In [2]: print(tables.__version__)
    3.0.0
    In [3]: a = tables.open_file('test.h5', 'a')
    In [4]: b = tables.open_file('test.h5', 'a')
    In [5]: a is b
    Out[5]: True

All this is an implementation detail, it happened under the hood and the user
had no control over the process.

This kind of behaviour was considered a feature since it can speed up opening
of files in case of repeated opens and it also avoids any potential problem
related to multiple opens, a practice that the HDF5 developers recommend to
avoid (see also H5Fopen_ reference page).

The trick, of course, is that files are not opened multiple times at HDF5
level, rather an open file is referenced several times.

The big drawback of this approach is that there are really few chances to use
PyTables safely in a multi thread program.  Several bug reports have been
filed regarding this topic.

After long discussions about the possibility to actually achieve concurrent I/O
and about patterns that should be used for the I/O in concurrent programs
PyTables developers decided to remove the *black magic under the hood* and
allow the users to implement the patterns they want.

Starting from PyTables 3.1 file handles are no more re-used (*shared*) and
each call to the :func:`open_file` function returns a new file handle::

    In [1]: import tables
    In [2]: print tables.__version__
    3.1.0
    In [3]: a = tables.open_file('test.h5', 'a')
    In [4]: b = tables.open_file('test.h5', 'a')
    In [5]: a is b
    Out[5]: False

It is important to stress that the new implementation still has an internal
registry (implementation detail) and it is still **not thread safe**.
Just now a smart enough developer should be able to use PyTables in a
muti-thread program without too much headaches.

The new implementation behaves differently from the previous one, although the
API has not been changed.  Now users should pay more attention when they open a
file multiple times (as recommended in the `HDF5 reference`__ ) and they
should take care of using them in an appropriate way.

__ H5Fopen_

Please note that the :attr:`File.open_count` property was originally intended
to keep track of the number of references to the same file handle.
In PyTables >= 3.1, despite of the name, it maintains the same semantics, just
now its value should never be higher that 1.

.. note::

    HDF5 versions lower than 1.8.7 are not fully compatible with PyTables 3.1.
    A partial support to HDF5 < 1.8.7 is still provided but in that case
    multiple file opens are not allowed at all (even in read-only mode).


.. _pip: http://www.pip-installer.org
.. _Anaconda: https://store.continuum.io/cshop/anaconda
.. _Canopy: https://www.enthought.com/products/canopy
.. _`Christoph Gohlke suites`: http://www.lfd.uci.edu/~gohlke/pythonlibs
.. _`Issues with H5T_NATIVE_LDOUBLE`: https://forum.hdfgroup.org/t/issues-with-h5t-native-ldouble/2991
.. _Python: http://www.python.org
.. _Blosc: http://www.blosc.org
.. _numpy: http://www.numpy.org
.. _`Travis-CI`: https://travis-ci.org
.. _PEP8: http://www.python.org/dev/peps/pep-0008
.. _`Blosc Release Notes`: https://github.com/FrancescAlted/blosc/wiki/Release-notes
.. _H5Fopen: https://portal.hdfgroup.org/display/HDF5/Files
=======================================
 Release notes for PyTables 3.2 series
=======================================

:Author: PyTables Developers
:Contact: pytables-dev@googlegroups.com

.. py:currentmodule:: tables


Changes from 3.2.3 to 3.2.3.1
=============================

Fixed issues with pip install.


Changes from 3.2.2 to 3.2.3
===========================

Improvements
------------

- It is now possible to use HDF5 with the new shared library naming scheme
  (>= 1.8.10, hdf5.dll instead of hdf5dll.dll) on Windows (:issue:`540`).
  Thanks to Tadeu Manoel.
- Now :program: `ptdump` sorts output by node name and does not print a
  backtrace if file cannot be opened.
  Thanks to Zbigniew Jędrzejewski-Szmek.


Bugs fixed
----------

- Only run `tables.tests.test_basics.UnicodeFilename` if the filesystem
  encoding is utf-8. Closes :issue:`485`.
- Add lib64 to posix search path. (closes :issue:`507`)
  Thanks to Mehdi Sadeghi.
- Ensure cache entries are removed if fewer than 10 (closes :issue:`529`).
  Thanks to Graham Jones.
- Fix segmentation fault in a number of test cases that use
  :class:`index.Index` (closes :issue:`532` and :issue:`533`).
  Thanks to Diane Trout.
- Fixed the evaluation of transcendental functions when numexpr is
  compiled with VML support (closes :issue:`534`, PR #536).
  Thanks to Tom Kooij.
- Make sure that index classes use buffersizes that are a multiple
  of chunkshape[0] (closes :issue:`538`, PR #538).
  Thanks to Tom Kooij.
- Ensure benchmark paths exist before benchmarks are executed (PR #544).
  Thanks to rohitjamuar.

Other changes
-------------

- Minimum Cython_ version is now v0.21


.. _Cython: http://cython.org


Changes from 3.2.1.1 to 3.2.2
=============================

Bug fixed
---------

- Fix AssertionError in Row.__init_loop. See :issue:`477`.
- Fix issues with Cython 0.23. See :issue:`481`.
- Only run `tables.tests.test_basics.UnicodeFilename` if the filesystem
  encoding is utf-8. Closes :issue:`485`.
- Fix missing PyErr_Clear. See :issue:`486`.
- Fix the C type of some numpy attributes. See :issue:`494`.
- Cast selection indices to integer. See :issue:`496`.
- Fix indexesextension._keysort_string. Closes :issue:`497` and :issue:`498`.


Changes from 3.2.1 to 3.2.1.1
=============================

- Fix permission on distributed source distribution

Other changes
-------------

- Minimum Cython_ version is now v0.21


.. _Cython: http://cython.org


Changes from 3.2.0 to 3.2.1
===========================

Bug fixed
---------

- Fix indexesextension._keysort. Fixes :issue:`455`. Thanks to Andrew Lin.


Changes from 3.1.1 to 3.2.0
===========================

Improvements
------------

- The `nrowsinbuf` is better computed now for EArray/CArray having
  a small `chunkshape` in the main dimension.  Fixes #285.

- PyTables should be installable very friendly via pip, including NumPy
  being installed automatically in the unlikely case it is not yet
  installed in the system.  Thanks to Andrea Bedini.

- setup.py has been largely simplified and now it requires *setuptools*.
  Although we think this is a good step, please keep us informed this is
  breaking some installation in a very bad manner.

- setup.py now is able to used *pkg-config*, if available, to locate required
  libraries (hdf5, bzip2, etc.). The use of *pkg-config* can be controlled
  via setup.py command line flags or via environment variables.
  Please refer to the installation guide (in the *User Manual*) for details.
  Closes :issue:`442`.

- It is now possible to create a new node whose parent is a softlink to another
  group (see :issue:`422`). Thanks to Alistair Muldal.

- :class:`link.SoftLink` objects no longer need to be explicitly dereferenced.
  Methods and attributes of the linked object are now automatically accessed
  when the user acts on a soft-link (see :issue:`399`).
  Thanks to Alistair Muldal.

- Now :program:`ptrepack` recognizes hardlinks and replicates them in the
  output (*repacked*) file. This saves disk space and makes repacked files
  more conformal to the original one. Closes :issue:`380`.

- New :program:`pttree` script for printing HDF5 file contents as a pretty
  ASCII tree (closes :issue:`400`). Thanks to Alistair Muldal.

- The internal Blosc library has been downgraded to version 1.4.4.  This
  is in order to still allow using multiple threads *inside* Blosc, even
  on multithreaded applications (see :issue:`411`, :issue:`412`,
  :issue:`437` and :issue:`448`).

- The :func:`print_versions` function now also reports the version of
  compression libraries used by Blosc.

- Now the :file:`setup.py` tries to use the '-march=native' C flag by
  default. In falls back on '-msse2' if '-march=native' is not supported
  by the compiler. Closes :issue:`379`.

- Fixed a spurious unicode comparison warning (closes :issue:`372` and
  :issue:`373`).

- Improved handling of empty string attributes. In previous versions of
  PyTables empty string were stored as scalar HDF5 attributes having size 1
  and value '\0' (an empty null terminated string).
  Now empty string are stored as HDF5 attributes having zero size

- Added a new cookbook recipe and a couple of examples for simple threading
  with PyTables.

- The redundant :func:`utilsextension.get_indices` function has been
  eliminated (replaced by :meth:`slice.indices`). Closes :issue:`195`.

- Allow negative indices in point selection (closes :issue:`360`)

- Index wasn't being used if it claimed there were no results.
  Closes :issue:`351` (see also :issue:`353`)

- Atoms and Col types are no longer generated dynamically so now it is easier
  for IDEs and static analysis tool to handle them (closes :issue:`345`)

- The keysort functions in idx-opt.c have been cythonised using fused types.
  The perfomance is mostly unchanged, but the code is much more simpler now.
  Thanks to Andrea Bedini.

- Small unit tests re-factoring:

  * :func:`print_versions` and :func:`tests.common.print_heavy` functions
     moved to the :mod:`tests.common` module

  * always use :func:`print_versions` when test modules are called as scripts

  * use the unittest2_ package in Python 2.6.x

  * removed internal machinery used to replicate unittest2_ features

  * always use :class:`tests.common.PyTablesTestCase` as base class for all
    test cases

  * code of the old :func:`tasts.common.cleanup` function has been moved to
    :meth:`tests.common.PyTablesTestCase.tearDown` method

  * new implementation of :meth:`tests.common.PyTablesTestCase.assertWarns`
    compatible with the one provided by the standard :mod:`unittest` module
    in Python >= 3.2

  * use :meth:`tests.common.PyTablesTestCase.assertWarns` as context manager
    when appropriate

  * use the :func:`unittest.skipIf` decorator when appropriate

  * new :class:tests.comon.TestFileMixin: class


.. _unittest2: https://pypi.python.org/pypi/unittest2


Bugs fixed
----------

- Fixed compatibility problems with numpy 1.9 and 1.10-dev
  (closes :issue:`362` and :issue:`366`)

- Fixed compatibility with Cython >= 0.20 (closes :issue:`386` and
  :issue:`387`)

- Fixed support for unicode node names in LRU cache (only Python 2 was
  affected). Closes :issue:`367` and :issue:`369`.

- Fixed support for unicode node titles (only Python 2 was affected).
  Closes :issue:`370` and :issue:`374`.

- Fixed a bug that caused the silent truncation of unicode attributes
  containing the '\0' character. Closes :issue:`371`.

- Fixed :func:`descr_from_dtype` to work as expected with complex types.
  Closes :issue:`381`.

- Fixed the :class:`tests.test_basics.ThreadingTestCase` test case.
  Closes :issue:`359`.

- Fix incomplete results when performing the same query twice and exhausting
  the second iterator before the first. The first one writes incomplete
  results to *seqcache* (:issue:`353`)

- Fix false results potentially going to *seqcache* if
  :meth:`tableextension.Row.update` is used during iteration
  (see :issue:`353`)

- Fix :meth:`Column.create_csindex` when there's NaNs

- Fixed handling of unicode file names on windows (closes :issue:`389`)

- No longer not modify :data:`sys.argv` at import time (closes :issue:`405`)

- Fixed a performance issue on NFS (closes :issue:`402`)

- Fixed a nasty problem affecting results of indexed queries.
  Closes :issue:`319` and probably :issue:`419` too.

- Fixed another problem affecting results of indexed queries too.
  Closes :issue:`441`.

- Replaced "len(xrange(start, stop, step))" -> "len(xrange(0, stop -
  start, step))" to fix issues with large row counts with Python 2.x.
  Fixes #447.


Other changes
-------------

- Cython is not a hard dependency anymore (although developers will need it
  so as to generated the C extension code).

- The number of threads used by default for numexpr and Blosc operation that
  was set to the number of available cores have been reduced to 2.  This is
  a much more reasonable setting for not creating too much overhead.


  **Enjoy data!**

  -- The PyTables Developers


.. Local Variables:
.. mode: rst
.. coding: utf-8
.. fill-column: 72
.. End:
============================
 What's new in PyTables 1.0
============================


:Author: Francesc Altet
:Contact: faltet@carabos.com
:Author: Ivan Vilata i Balaguer
:Contact: ivilata@carabos.com


This document details the modifications to PyTables since version 0.9.1.  Its
main purpose is help you ensure that your programs will be runnable when you
switch from PyTables 0.9.1 to PyTables 1.0.


API additions
=============

- The new ``Table.col()`` method can be used to get a column from a table as a
  ``NumArray`` or ``CharArray`` object.  This is preferred over the syntax
  ``table['colname']``.

- The new ``Table.readCoordinates()`` method reads a set of rows given their
  indexes into an in-memory object.

- The new ``Table.readAppend()`` method Append rows fulfilling the condition
  to a destination table.

Backward-incompatible changes
=============================

- Trying to open a nonexistent file or a file of unknown type raises
  ``IOError`` instead of ``RuntimeError``.  Using an invalid mode raises
  ``ValueError`` instead of ``RuntimeError``.

- Getting a child node from a closed group raises ``ValueError`` instead of
  ``RuntimeError``.

- Running an action on the wrong type of node now (i.e. using
  ``file.listNodes()`` on a leaf) raises a ``TypeError`` instead of a
  ``NodeError``.

- Removing a non-existing child now raises a ``NoSuchNodeError``, instead of
  doing nothing.

- Removing a non-empty child group using ``del group.child`` fails with a
  ``NodeError`` instead of recursively doing the removal.  This is because of
  the potential damage it may cause when used inadvertently.  If a recursive
  behavior is needed, use the ``_f_remove()`` method of the child node.

- The `recursive` flag of ``Group._f_walkNodes()`` is ``True`` by default now.
  Before it was ``False``.

- Now, deleting and getting a non-existing attribute raises an
  ``AttributeError`` instead of a ``RuntimeError``.

- Swapped last two arguments of ``File.copyAttrs()`` to match the other
  methods.  Please use ``File.copyNodeAttrs()`` anyway.

- Failing to infer the size of a string column raises ``ValueError`` instead
  of ``RuntimeError``.

- Excessive table column name length and number of columns now raise
  ``ValueError`` instead of ``IndexError`` and ``NameError``.

- Excessive table row length now raises ``ValueError`` instead of
  ``RuntimeError``.

- ``table[integer]`` returns a ``numarray.records.Record`` object instead of a
  tuple.  This was the original behavior before PyTables 0.9 and proved to be
  more consistent than the last one (tables do not have an explicit ordering
  of columns).

- Specifying a nonexistent column in ``Table.read()`` raises a ``ValueError``
  instead of a ``LookupError``.

- When ``start >= stop`` an empty iterator is returned by ``Table.iterrows()``
  instead of an empty ``RecArray``.  Thanks to Ashley Walsh for noting this.

- The interface of ``isHDF5File()`` and ``isPyTablesFile()`` file has been
  unified so that they both return true or false values on success and raise
  ``HDF5ExtError`` or errors.  The true value in ``isPyTablesFile()`` is the
  format version string of the file.

- ``Table.whereIndexed()`` and ``Table.whereInRange()`` are now *private*
  methods, since the ``Table.where()`` method is able to choose the most
  adequate option.

- The global variables ``ExtVersion`` and ``HDF5Version`` have been renamed to
  ``extVersion`` and ``hdf5Version``, respectively.

- ``whichLibVersion()`` returns ``None`` on querying unavailable libraries,
  and raises ``ValueError`` on unknown ones.

The following modifications, though being (strictly speaking) modifications of
the API, will most probably not cause compatibility problems (but your mileage
may vary):

- The default values for ``name`` and ``classname`` arguments in
  ``File.getNode()`` are now ``None``, although the empty string is still
  allowed for backwards compatibility.  File hierarchy manipulation and
  attribute handling operations using those arguments have changed to reflect
  this.

- Copy operations (``Group._f_copyChildren()``, ``File.copyChildren()``,
  ``File.copyNode()``...) do no longer return a tuple with the new node and
  statistics.  Instead, they only return the new node, and statistics are
  collected via an optional keyword argument.

- The ``copyFile()`` function in ``File.py`` has changed its signature from::

      copyFile(srcfilename=None, dstfilename=None, title=None, filters=None,
               copyuserattrs=True, overwrite=False, stats=None)

  to::

      copyFile(srcfilename, dstfilename, overwrite=False, **kwargs)

  Thus, the function allows the same options as ``File.copyFile()``.

- The ``File.copyFile()`` method has changed its signature from::

      copyFile(self, dstfilename=None, title=None, filters=None,
               copyuserattrs=1, overwrite=0, stats=None):

  to::

      copyFile(self, dstfilename, overwrite=False, **kwargs)

  This enables this method to pass on arbitrary flags and options supported by
  copying methods of inner nodes in the hierarchy.

- The ``File.copyChildren()`` method has changed its signature from::

      copyChildren(self, wheresrc, wheredst, recursive=False, filters=None,
                   copyuserattrs=True, start=0, stop=None, step=1,
                   overwrite=False, stats=None)

  to::

      copyChildren(self, srcgroup, dstgroup, overwrite=False, recursive=False,
                   **kwargs):

  Thus, the function allows the same options as ``Group._f_copyChildren()``.

- The ``Group._f_copyChildren()`` method has changed its signature from::

      _f_copyChildren(self, where, recursive=False, filters=None,
                      copyuserattrs=True, start=0, stop=None, step=1,
                      overwrite=False, stats=None)

  to::

      _f_copyChildren(self, dstgroup, overwrite=False, recursive=False,
                      **kwargs)

  This enables this method to pass on arbitrary flags and options supported by
  copying methods of inner nodes in the group.

- Renamed ``srcFilename`` and ``dstFilename`` arguments in ``copyFile()`` and
  ``File.copyFile()`` to ``srcfilename`` and ``dstfilename``, respectively.
  Renamed ``whereSrc`` and ``whereDst`` arguments in ``File.copyChildren()``
  to ``wheresrc`` and ``wheredst``, respectively.  Renamed ``dstNode``
  argument in ``File.copyAttrs()`` to ``dstnode``.  Tose arguments should be
  easier to type in interactive sessions (although 99% of the time it is not
  necessary to specify them).

- Renamed ``object`` argument in ``EArray.append()`` to ``sequence``.

- The ``rows`` argument in ``Table.append()`` is now compulsory.

- The ``start`` argument in ``Table.removeRows()`` is now compulsory.


API refinements
===============

- The ``isHDF5()`` function has been deprecated in favor of ``isHDF5File()``.

- Node attribute-handling methods in ``File`` have been renamed for a better
  coherence and understanding of their purpose:

  * ``getAttrNode()`` is now called ``getNodeAttr()``
  * ``setAttrNode()`` is now called ``setNodeAttr()``
  * ``delAttrNode()`` is now called ``delNodeAttr()``
  * ``copyAttrs()`` is now called ``copyNodeAttrs()``

  They keep their respective signatures, and the old versions still exist for
  backwards compatibility, though they issue a ``DeprecationWarning``.

- Using ``VLArray.append()`` with multiple arguments is now deprecated for its
  ambiguity.  You should put the arguments in a single sequence object (list,
  tuple, array...) and pass it as the only argument.

- Using ``table['colname']`` is deprecated.  Using ``table.col('colname')``
  (with the new ``col()`` method) is preferred.


Bug fixes (affecting API)
=========================

- ``Table.iterrows()`` returns an empty iterator when no rows are selected,
  instead of returning ``None``.


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: text
.. coding: utf-8
.. fill-column: 78
.. End:
==============================
 What's new in PyTables 1.2.1
==============================


:Author: Francesc Altet
:Contact: faltet@carabos.com
:Author: Ivan Vilata i Balaguer
:Contact: ivilata@carabos.com


This document details the modifications to PyTables since version 1.2.  Its
main purpose is help you ensure that your programs will be runnable when you
switch from PyTables 1.2 to PyTables 1.2.1.


API additions
=============

- None

Backward-incompatible changes
=============================

- None

Deprecated features
===================

- None


API refinements
===============

- None


Bug fixes (affecting API)
=========================

- None


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: text
.. coding: utf-8
.. fill-column: 78
.. End:
============================
 What's new in PyTables 1.3
============================


:Author: Francesc Altet
:Contact: faltet@carabos.com
:Author: Ivan Vilata i Balaguer
:Contact: ivilata@carabos.com


This document details the modifications to PyTables since version 1.2.  Its
main purpose is help you ensure that your programs will be runnable when you
switch from PyTables 1.2 to PyTables 1.3.


API additions
=============

- The Table.Cols accessor has received a new __setitem__() method that
  allows doing things like:

            table.cols[4] = record
            table.cols.x[4:1000:2] = array   # homogeneous column
            table.cols.Info[4:1000:2] = recarray   # nested column


Backward-incompatible changes
=============================

- None


Deprecated features
===================

- None


API refinements
===============

- Table.itersequence has changed the default value for 'sort' parameter. It is
  now False by default, as it is not clear if this actually accelerates the
  iterator, so it is better to let to the user doing the proper checks (if he
  is interested at all).


Bug fixes (affecting API)
=========================

- None


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: text
.. coding: utf-8
.. fill-column: 78
.. End:
=======================================
 Release notes for PyTables 2.0 series
=======================================

:Author: Francesc Alted i Abad
:Contact: faltet@pytables.com
:Author: Ivan Vilata i Balaguer
:Contact: ivan@selidor.net


Changes from 2.0.3 to 2.0.4
===========================

- Selections in tables works now in threaded environments.  The problem was in
  the Numexpr package -- the solution has been reported to the upstream
  authors too.  Fixes #164.

- PyTables had problems importing native HDF5 files with gaps in nested
  compound types.  This has been solved.  Fixes #173.

- In order to prevent a bug existing in HDF5 1.6 series, the
  ``EArray.truncate()`` method refused to accept a 0 as parameter
  (i.e. truncate an existing EArray to have zero rows did not work).  As this
  has been fixed in the recent HDF5 1.8 series, this limitation has been
  removed (but only if the user has one of these installed).  Fixes #171.

- Small fixes for allowing the test suite to pass when using the new NumPy
  1.1.  However, it remains a small issue with the way the new NumPy
  represents complex numbers.  I'm not fixing that in the PyTables suite, as
  there are chances that this can be fixed in NumPy itself (see ticket #841).


Changes from 2.0.2 to 2.0.3
===========================

- Replaced the algorithm for computing chunksizes by another that is
  more general and useful for a larger range of expected dataset sizes.
  The outcome of the new calculation is the same than before for
  dataset sizes <= 100 GB. For datasets between 100 GB <= size < 10
  TB, larger values are returned. For sizes >= 10 TB a maximum value
  of 1 MB is always returned.

- Fixed a problem when updating multidimensional cells using the
  Row.update() method in the middle of table iterators .  Fixes #149.

- Added support for the latest 1.8.0 version of the HDF5 library.
  Fixes ticket #127.

- PyTables compiles now against latest versions of Pyrex (0.9.6.4).  For the
  first time, the extensions do compile without warnings!  Fixes #159.

- Numexpr module has been put in sync with the version in SciPy sandbox.

- Added a couple of warnings in User's Guide so as to tell the user that it is
  not safe to use methods that can change the number of rows of a table in the
  middle of a row iterator. Fixes #153.


Changes from 2.0.1 to 2.0.2
===========================

- Added ``__enter__()`` and ``__exit__()`` methods to ``File``; fixes #113.
  With this, and if using Python 2.5 you can do things like:

    with tables.openFile("test.h5") as h5file:
        ...

- Carefully preserve type when converting NumPy scalar to numarray; fixes
  #125.

- Fixed a nasty bug that appeared when moving or renaming groups due to a bad
  interaction between ``Group._g_updateChildrenLocation()`` and the LRU cache.
  Solves #126.

- Return 0 when no rows are given to ``Table.modifyRows()``; fixes #128.

- Added an informative message when the ``nctoh5`` utility is run without the
  NetCDF interface of ScientificPython being installed.

- Now, a default representation of closed nodes is provided; fixes #129.


Changes from 2.0 to 2.0.1
=========================

- The ``coords`` argument of ``Table.readCoords()`` was not checked
  for contiguousness, raising fatal errors when it was discontiguous.
  This has been fixed.

- There is an inconsistency in the way used to specify the atom shape
  in ``Atom`` constructors.  When the shape is specified as
  ``shape=()`` it means a scalar atom and when it is specified as
  ``shape=N`` it means an atom with ``shape=(N,)``.  But when the
  shape is specified as ``shape=1`` (i.e. in the default case) then a
  scalar atom is obtained instead of an atom with ``shape=(1,)``.
  This is inconsistent and not the behavior that NumPy exhibits.

  Changing this will require a migration path which includes
  deprecating the old behaviour if we want to make the change happen
  before a new major version.  The proposed path is:

   1. In PyTables 2.0.1, we are changing the default value of the
      ``shape`` argument to ``()``, and issue a ``DeprecationWarning``
      when someone uses ``shape=1`` stating that, for the time being,
      it is equivalent to ``()``, but in near future versions it will
      become equivalent to ``(1,)``, and recommending the user to pass
      ``shape=()`` if a scalar is desired.

   2. In PyTables 2.1, we will remove the previous warning and take
      ``shape=N`` to mean ``shape=(N,)`` for any value of N.

  See ticket #96 for more info.

- The info about the ``chunkshape`` attribute of a leaf is now printed
  in the ``__repr__()`` of chunked leaves (all except ``Array``).

- After some scrupulous benchmarking job, the size of the I/O buffer
  for ``Table`` objects has been reduced to the minimum that allows
  maximum performance.  This represents more than 10x of reduction in
  size for that buffer, which will benefit those programs dealing with
  many tables simultaneously (#109).

- In the ``ptrepack`` utility, when ``--complevel`` and ``--shuffle``
  were specified at the same time, the 'shuffle' filter was always set
  to 'off'.  This has been fixed (#104).

- An ugly bug related with the integrated Numexpr not being aware of
  all the variations of data arrangements in recarray objects has been
  fixed (#103).  We should stress that the bug only affected the
  Numexpr version integrated in PyTables, and *not* the original one.

- When passing a record array to a table at creation time, its real
  length is now used instead of the default value for
  ``expectedrows``.  This allows for better performance (#97).

- Added some workarounds so that NumPy scalars can be successfully
  converted to numarray objects.  Fixes #98.

- PyTables is now able to access table rows beyond 2**31 in 32-bit
  Python.  The problem was a limitation of ``xrange`` and we have
  replaced it by a new ``lrange`` class written in Pyrex.  Moreover,
  ``lrange`` has been made publicly accessible as a safe 64-bit
  replacement for ``xrange`` for 32-bit platforms users.  Fixes #99.

- If a group and a table are created in a function, and the table is
  accessed through the group, the table can be flushed now.  Fixes
  #94.

- It is now possible to directly assign a field in a nested record of
  a table using the natural naming notation (#93).


Changes from 2.0rc2 to 2.0
==========================

- Added support for recognizing native HDF5 files with datasets compressed
  with szip compressor.

- Fixed a problem when asking for the string representation (str()) of closed
  files. Fixes ticket #79.

- Do not take LZO as available when its initialisation fails.

- Fixed a glitch in ptrepack utility. When the user wants a copy of a group,
  and a group is *to be created* in destination, the attributes of the
  original group *are* copied. If it is *not to be created*, the attributes
  will *not be* copied. I think this should be what the user would expect most
  of the times.

- Fixed the check for creating intermediate groups in ptrepack utility.
  Solves ticket #83.

- Before, when reading a dataset with an unknown CLASS id, a warning was
  issued and the dataset mapped to ``UnImplemented``. This closed the door to
  have the opportunity to try to recognize the dataset and map it to a
  supported CLASS. Now, when a CLASS attribute is not recognized, an attempt
  to recognize its associated dataset is made. If it is recognized, the
  matching class is associated with the dataset. If it is not recognized, then
  a warning is issued and the dataset becomes mapped to ``UnImplemented``.

- Always pass verbose and heavy values in the common test module to test().
  Fixes ticket #85.

- Now, the ``verbose`` and ``--heavy`` flag passed to test_all.py are honored.

- All the DLL's of dependencies are included now in Windows binaries.  This
  should allow for better portability of the binaries.

- Fixed the description of Node._v_objectID that was misleading.


Changes from 2.0rc1 to 2.0rc2
=============================

- The "Optimization tips" chapter of the User's Guide has been completely
  updated to adapt to PyTables 2.0 series.  In particular, new benchmarks on
  the much improved indexed queries have been included; you will see that
  PyTables indexing is competitive (and sometimes much faster) than that of
  traditional relational databases.  With this, the manual should be fairly
  finished for 2.0 final release.

- Large refactoring done on the ``Row`` class.  The most important change is
  that ``Table.row`` is now a single object.  This allows to reuse the same
  ``Row`` instance even after ``Table.flush()`` calls, which can be convenient
  in many situations.

- I/O buffers unified in the ``Row`` class.  That allows for bigger savings in
  memory space whenever the ``Row`` extension is used.

- Improved speed (up to a 70%) with unaligned column operations (a quite
  common scenario when dealing with ``Table`` objects) through the integrated
  Numexpr.  In-kernel searches take advantage of this optimization.

- Added ``VLUnicodeAtom`` for storing variable-length Unicode strings in
  ``VLArray`` objects regardless of encoding.  Closes ticket #51.

- Added support for ``time`` datatypes to be portable between big-endian and
  low-endian architectures.  This feature is not currently supported natively
  by the HDF5 library, so the support for such conversion has been added in
  PyTables itself.  Fixes #72.

- Added slice arguments to ``Table.readWhere()`` and ``Table.getWhereList()``.
  Although API changes are frozen, this may still be seen as an inconsistency
  with other query methods.  The patch is backwards-compatible anyway.

- Added missing overwrite argument to ``File.renameNode()`` and
  ``Node._f_rename()``.  Fixes ticket #66.

- Calling ``tables.test()`` no longer exits the interpreter session.  Fixes
  ticket #67.

- Fix comparing strings where one is a prefix of the other in integrated
  Numexpr.  Fixes ticket #76.

- Added a check for avoiding an ugly HDF5 message when copying a file over
  itself (for both ``copyFile()`` and ``File.copyFile()``).  Fixes ticket #73.

- Corrected the appendix E, were it was said that PyTables doesn't support
  compounds of compounds (it does since version 1.2!).


Changes from 2.0b2 to 2.0rc1
============================

- The API Reference section of the User's Manual (and the matching docstrings)
  has been completely reviewed, expanded and corrected.  This process has
  unveiled some errors and inconsistencies which have also been fixed.

- Fixed ``VLArray.__getitem__()`` to behave as expected in Python when using
  slices, instead of following the semantics of PyTables' ``read()`` methods
  (e.g. reading just one element when no stop is provided).  Fixes ticket #50.

- Removed implicit UTF-8 encoding from ``VLArray`` data using ``vlstring``
  atoms.  Now a variable-length string is stored as is, which lets users use
  any encoding of their choice, or none of them.  A ``vlunicode`` atom will
  probably be added to the next release so as to fix ticket #51.

- Allow non-sequence objects to be passed to ``VLArray.append()`` when using
  an ``object`` atom.  This was already possible in 1.x but stopped working
  when the old append syntax was dropped in 2.0.  Fixes ticket #63.

- Changed ``Cols.__len__()`` to return the number of rows of the table or
  nested column (instead of the number of fields), like its counterparts in
  ``Table`` and ``Column``.

- Python scalars cached in ``AttributeSet`` instances are now kept as NumPy
  objects instead of Python ones, because they do become NumPy objects when
  retrieved from disk.  Fixes ticket #59.

- Avoid HDF5 error when appending an empty array to a ``Table`` (ticket #57)
  or ``EArray`` (ticket #49) dataset.

- Fix wrong implementation of the top-level ``table.description._v_dflts``
  map, which was also including the pathnames of columns inside nested
  columns.  Fixes ticket #45.

- Optimized the access to unaligned arrays in Numexpr between a 30% and a 70%.

- Fixed a die-hard bug that caused the loading of groups while closing a file.
  This only showed with certain usage patterns of the LRU cache (e.g. the one
  caused by ``ManyNodesTestCase`` in ``test_indexes.py`` under Pro).

- Avoid copious warnings about unused functions and variables when compiling
  Numexpr.

- Several fixes to Numexpr expressions with all constant values.  Fixed
  tickets #53, #54, #55, #58.  Reported bugs to mainstream developers.

- Solved an issue when trying to open one of the included test files in append
  mode on a system-wide installation by a normal user with no write privileges
  on it.  The file isn't being modified anyway, so the test is skipped then.

- Added a new benchmark to compare the I/O speed of ``Array`` and ``EArray``
  objects with that of ``cPickle``.

- The old ``Row.__call__()`` is no longer available as a public method.  It
  was not documented, anyway.  Fixes ticket #46.

- ``Cols._f_close()`` is no longer public.  Fixes ticket #47.

- ``Attributes._f_close()`` is no longer public.  Fixes ticket #52.

- The undocumented ``Description.classdict`` attribute has been completely
  removed.  Fixes ticket #44.


Changes from 2.0b1 to 2.0b2
===========================

- A very exhaustive overhauling of the User's Manual is in process.  The
  chapters 1 (Introduction), 2 (Installation), 3 (Tutorials) have been
  completed (and hopefully, the lines of code are easier to copy&paste now),
  while chapter 4 (API Reference) has been done up to (and including) the
  Table class.  During this tedious (but critical in a library) overhauling
  work, we have tried hard to synchronize the text in the User's Guide with
  that which appears on the docstrings.

- Removed the ``recursive`` argument in ``Group._f_walkNodes()``.  Using it
  with a false value was redundant with ``Group._f_iterNodes()``.  Fixes
  ticket #42.

- Removed the ``coords`` argument from ``Table.read()``.  It was undocumented
  and redundant with ``Table.readCoordinates()``.  Fixes ticket #41.

- Fixed the signature of ``Group.__iter__()`` (by removing its parameters).

- Added new ``Table.coldescrs`` and ``Table.description._v_itemsize``
  attributes.

- Added a couple of new attributes for leaves:

  * ``nrowsinbuf``: the number of rows that fit in the internal buffers.
  * ``chunkshape``: the chunk size for chunked datasets.

- Fixed setuptools so that making an egg out of the PyTables 2 package is
  possible now.

- Added a new ``tables.restrict_flavors()`` function allowing to restrict
  available flavors to a given set.  This can be useful e.g. if you want to
  force PyTables to get NumPy data out of an old, ``numarray``-flavored
  PyTables file even if the ``numarray`` package is installed.

- Fixed a bug which caused filters of unavailable compression libraries to be
  loaded as using the default Zlib library, after issuing a warning.  Added a
  new ``FiltersWarning`` and a ``Filters.copy()``.


Important changes from 1.4.x to 2.0
===================================

API additions
-------------

- ``Column.createIndex()`` has received a couple of new parameters:
  ``optlevel`` and ``filters``.  The first one sets the desired quality level
  of the index, while the second one allows the user to specify the filters
  for the index.

- ``Table.indexprops`` has been split into ``Table.indexFilters`` and
  ``Table.autoIndex``.  The later groups the functionality of the old ``auto``
  and ``reindex``.

- The new ``Table.colpathnames`` is a sequence which contains the full
  pathnames of all bottom-level columns in a table.  This can be used to walk
  all ``Column`` objects in a table when used with ``Table.colinstances``.

- The new ``Table.colinstances`` dictionary maps column pathnames to their
  associated ``Column`` or ``Cols`` object for simple or nested columns,
  respectively.  This is similar to ``Table.cols._f_col()``, but faster.

- ``Row`` has received a new ``Row.fetch_all_fields()`` method in order to
  return all the fields in the current row.  This returns a NumPy void scalar
  for each call.

- New ``tables.test(verbose=False, heavy=False)`` high level function for
  interactively running the complete test suite from the Python console.

- Added a ``tables.print_versions()`` for easily getting the versions for all
  the software on which PyTables relies on.


Backward-incompatible changes
-----------------------------

- You can no longer mark a column for indexing in a ``Col`` declaration.  The
  only way of creating an index for a column is to invoke the
  ``createIndex()`` method of the proper column object *after the table has
  been created*.

- Now the ``Table.colnames`` attribute is just a list of the names of
  top-level columns in a table.  You can still get something similar to the
  old structure by using ``Table.description._v_nestedNames``.  See also the
  new ``Table.colpathnames`` attribute.

- The ``File.objects``, ``File.leaves`` and ``File.groups`` dictionaries have
  been removed.  If you still need this functionality, please use the
  ``File.getNode()`` and ``File.walkNodes()`` instead.

- ``Table.removeIndex()`` is no longer available; to remove an index on a
  column, one must use the ``removeIndex()`` method of the associated
  ``Column`` instance.

- ``Column.dirty`` is no longer available.  If you want to check
  column index dirtiness, use ``Column.index.dirty``.

- ``complib`` and ``complevel`` parameters have been removed from
  ``File.createTable()``, ``File.createEArray()``, ``File.createCArray()`` and
  ``File.createVLArray()``.  They were already deprecated in PyTables 1.x.

- The ``shape`` and ``atom`` parameters have been swapped in
  ``File.createCArray()``.  This has been done to be consistent with
  ``Atom()`` definitions (i.e. type comes before and shape after).

Deprecated features
-------------------

- ``Node._v_rootgroup`` has been removed.  Please use ``node._v_file.root``
  instead.

- The ``Node._f_isOpen()`` and ``Leaf.isOpen()`` methods have been removed.
  Please use the ``Node._v_isopen`` attribute instead (it is much faster).

- The ``File.getAttrNode()``, ``File.setAttrNode()`` and
  ``File.delAttrNode()`` methods have been removed.  Please use
  ``File.getNodeAttr()``, ``File.setNodeAttr()`` and ``File.delNodeAttr()``
  instead.

- ``File.copyAttrs()`` has been removed.  Please use ``File.copyNodeAttrs()``
  instead.

- The ``table[colname]`` idiom is no longer supported.  You can use
  ``table.cols._f_col(column)`` for doing the same.

API refinements
---------------

- ``File.createEArray()`` received a new ``shape`` parameter.  This allows to
  not have to use the shape of the atom so as to set the shape of the
  underlying dataset on disk.

- All the leaf constructors have received a new ``chunkshape`` parameter that
  allows specifying the chunk sizes of datasets on disk.

- All ``File.create*()`` factories for ``Leaf`` nodes have received a new
  ``byteorder`` parameter that allows the user to specify the byteorder in
  which data will be written to disk (data in memory is now always handled in
  *native* order).


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: rst
.. coding: utf-8
.. fill-column: 78
.. End:
=======================================
 Release notes for PyTables 2.1 series
=======================================

:Author: Francesc Alted i Abad
:Contact: faltet@pytables.org


Changes from 2.1.1 to 2.1.2
===========================

Bug fixes
---------

- Solved problems with Table.modifyColumn() when the column(s) is
  multidimensional. Fixes #228.

- The row attribute of a table seems stalled after a table move or
  rename.  Fixes #224.

- Fixed a problem with ``len(array)`` in 32-bit platforms when array
  is large enough (> 2**31).

- Added missing `_c_classId` attribute to the `UnImplemented` class.
  ``ptrepack`` no longer chokes while copying `Unimplemented` classes.

- The ``FIELD_*`` sys attrs are no longer copied when the
  ``PYTABLES_SYS_ATTRS`` parameter is set to false.

- The ``FILTERS`` attribute is not added anymore when
  ``PYTABLES_SYS_ATTR`` parameter is set to false.

- Disable the printing of Unicode characters that cannot be printed on
  win32 platform. Fixes #235.

Other changes
-------------

- When retrieving a row of a 1-dimensional array, a 0-dim array was
  returned instead of a numpy scalar.  Now, an actual numpy scalar is
  returned.  Closes #222.

- LZO and bzip2 filters adapted to an API fix introduced in HDF5
  1.8.3.  Closes #225.

- Unsupported HDF5 types in attributes are no longer transferred
  during copies. A new `_v_unimplemented` list have been added in
  `AttributeSet` class so as to keep track of such attributes.  Closes
  #240.

- LZO binaries have disappeared from the GnuWin32 repository.  Until
  they come eventually back, they have been put at
  http://www.pytables.org/download/lzo-win. This has been documented
  in the install chapter.


Changes from 2.1 to 2.1.1
=========================

Bug fixes
---------

- Fixed a memory leak when a lot of queries were made.  Closes #203
  and #207.

- The chunkshape="auto" parameter value of `Leaf.copy()` is honored
  now, even when the (start, stop, step) parameters are specified.
  Closes #204.

- Due to a flaw in its design, the `File` class was not able to be
  subclassed.  This has been fixed.  Closes #205.

- Default values were not correctly retrieved when opening already
  created CArray/EArray objects.  Fixed.  Closes #212.

- Fixed a problem with the installation of the ``nctoh5`` script that
  prevented it from being executed.  Closes #215.

- [Pro] The ``iterseq`` cache ignored non-indexed conditions, giving
  wrong results when those appeared in condition expressions.  This
  has been fixed.  Closes #206.

Other changes
-------------

- `openFile()`, `isHDF5File()` and `isPyTablesFile()` functions accept
  Unicode filenames now.  Closes #202 and #214.

- When creating large type sizes (exceeding 64 KB), HDF5 complained
  and refused to do so.  The HDF5 team has logged the issue as a bug,
  but meanwhile it has been implemented a workaround in PyTables that
  allows to create such large datatypes for situations that does not
  require defaults other than zero.  Addresses #211.

- In order to be consistent with how are stored the other data types,
  Unicode attributes are retrieved now as NumPy scalars instead of
  Python Unicode strings or NumPy arrays.  For the moment, I've fixed
  this through pickling the Unicode strings.  In the future, when HDF5
  1.8.x series would be a requirement, that should be done via a HDF5
  native Unicode type.  Closes #213.


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: rst
.. coding: utf-8
.. fill-column: 72
.. End:
===========================================
 Release notes for PyTables Pro 2.0 series
===========================================

:Author: Francesc Alted i Abad
:Contact: faltet@pytables.com
:Author: Ivan Vilata i Balaguer
:Contact: ivan@selidor.net


Changes from 2.0.3 to 2.0.4
===========================

- Selections in tables works now in threaded environments.  The problem was in
  the Numexpr package -- the solution has been reported to the upstream
  authors too.  Fixes #164.

- PyTables had problems importing native HDF5 files with gaps in nested
  compound types.  This has been solved.  Fixes #173.

- In order to prevent a bug existing in HDF5 1.6 series, the
  ``EArray.truncate()`` method refused to accept a 0 as parameter
  (i.e. truncate an existing EArray to have zero rows did not work).  As this
  has been fixed in the recent HDF5 1.8 series, this limitation has been
  removed (but only if the user has one of these installed).  Fixes #171.

- Small fixes for allowing the test suite to pass when using the new NumPy
  1.1.  However, it remains a small issue with the way the new NumPy
  represents complex numbers.  I'm not fixing that in the PyTables suite, as
  there are chances that this can be fixed in NumPy itself (see ticket #841).


Changes from 2.0.2.1 to 2.0.3
=============================

- Replaced the algorithm for computing chunksizes by another that is
  more general and useful for a larger range of expected dataset
  sizes. The outcome of the new calculation is the same than before
  for dataset sizes <= 100 GB. For datasets between 100 GB <= size <
  10 TB, larger values are returned. For sizes >= 10 TB a maximum
  value of 1 MB is always returned.

- Added support for the latest 1.8.0 version of the HDF5 library.
  Fixes ticket #127.

- PyTables compiles now against latest versions of Pyrex (0.9.6.4).  For the
  first time, the extensions do compile without warnings!  Fixes #159.

- Numexpr module has been put in sync with the version in SciPy sandbox.

- Added a couple of warnings in User's Guide so as to tell the user that it is
  not safe to use methods that can change the number of rows of a table in the
  middle of a row iterator. Fixes #153.

- Fixed a problem when updating multidimensional cells using the
  Row.update() method in the middle of table iterators .  Fixes #149.

- Fixed a problem when using 64-bit indexes in 32-bit platforms.
  Solves ticket #148.

- Table.indexFilters is working now as documented. However, as per ticket
  #155, its use is now deprecated (will be removed in 2.1). Fixes #155.


Changes from 2.0.2 to 2.0.2.1
=============================

- Optimization added for avoid to unnecessarily update index columns
  that have not been modified in table update operations.  Fixes #139.


Changes from 2.0.1 to 2.0.2
===========================

- Fixed a critical bug that returned wrong results when doing repetitive
  queries affecting the last row part of indices.  Fixes #60 of the private
  Trac of Carabos.

- Added ``__enter__()`` and ``__exit__()`` methods to ``File``; fixes #113.
  With this, and if using Python 2.5 you can do things like:

    with tables.openFile("test.h5") as h5file:
        ...

- Carefully preserve type when converting NumPy scalar to numarray; fixes
  #125.

- Fixed a nasty bug that appeared when moving or renaming groups due to a bad
  interaction between ``Group._g_updateChildrenLocation()`` and the LRU cache.
  Solves #126.

- Return 0 when no rows are given to ``Table.modifyRows()``; fixes #128.

- Added an informative message when the ``nctoh5`` utility is run without the
  NetCDF interface of ScientificPython bening installed.

- Now, a default representation of closed nodes is provided; fixes #129.


Changes from 2.0 to 2.0.1
=========================

- The ``coords`` argument of ``Table.readCoords()`` was not checked
  for contiguousness, raising fatal errors when it was discontiguous.
  This has been fixed.

- There is an inconsistency in the way used to specify the atom shape
  in ``Atom`` constructors.  When the shape is specified as
  ``shape=()`` it means a scalar atom and when it is specified as
  ``shape=N`` it means an atom with ``shape=(N,)``.  But when the
  shape is specified as ``shape=1`` (i.e. in the default case) then a
  scalar atom is obtained instead of an atom with ``shape=(1,)``.
  This is inconsistent and not the behavior that NumPy exhibits.

  Changing this will require a migration path which includes
  deprecating the old behaviour if we want to make the change happen
  before a new major version.  The proposed path is:

   1. In PyTables 2.0.1, we are changing the default value of the
      ``shape`` argument to ``()``, and issue a ``DeprecationWarning``
      when someone uses ``shape=1`` stating that, for the time being,
      it is equivalent to ``()``, but in near future versions it will
      become equivalent to ``(1,)``, and recommending the user to pass
      ``shape=()`` if a scalar is desired.

   2. In PyTables 2.1, we will remove the previous warning and take
      ``shape=N`` to mean ``shape=(N,)`` for any value of N.

  See ticket #96 for more info.

- The info about the ``chunkshape`` attribute of a leaf is now printed
  in the ``__repr__()`` of chunked leaves (all except ``Array``).

- After some scrupulous benchmarking job, the size of the I/O buffer
  for ``Table`` objects has been reduced to the minimum that allows
  maximum performance.  This represents more than 10x of reduction in
  size for that buffer, which will benefit those programs dealing with
  many tables simultaneously (#109).

- In the ``ptrepack`` utility, when ``--complevel`` and ``--shuffle``
  were specified at the same time, the 'shuffle' filter was always set
  to 'off'.  This has been fixed (#104).

- An ugly bug related with the integrated Numexpr not being aware of
  all the variations of data arrangements in recarray objects has been
  fixed (#103).  We should stress that the bug only affected the
  Numexpr version integrated in PyTables, and *not* the original one.

- When passing a record array to a table at creation time, its real
  length is now used instead of the default value for
  ``expectedrows``.  This allows for better performance (#97).

- Added some workarounds so that NumPy scalars can be successfully
  converted to numarray objects.  Fixes #98.

- PyTables is now able to access table rows beyond 2**31 in 32-bit
  Python.  The problem was a limitation of ``xrange`` and we have
  replaced it by a new ``lrange`` class written in Pyrex.  Moreover,
  ``lrange`` has been made publicly accessible as a safe 64-bit
  replacement for ``xrange`` for 32-bit platforms users.  Fixes #99.

- If a group and a table are created in a function, and the table is
  accessed through the group, the table can be flushed now.  Fixes
  #94.

- It is now possible to directly assign a field in a nested record of
  a table using the natural naming notation (#93).


Changes from 2.0rc2 to 2.0
==========================

- Added support for recognizing native HDF5 files with datasets compressed
  with szip compressor.

- Fixed a problem when asking for the string representation (str()) of closed
  files. Fixes ticket #79.

- Do not take LZO as available when its initialisation fails.

- Fixed a glitch in ptrepack utility. When the user wants a copy of a group,
  and a group is *to be created* in destination, the attributes of the
  original group *are* copied. If it is *not to be created*, the attributes
  will *not be* copied. I think this should be what the user would expect most
  of the times.

- Fixed the check for creating intermediate groups in ptrepack utility.
  Solves ticket #83.

- Before, when reading a dataset with an unknown CLASS id, a warning was
  issued and the dataset mapped to ``UnImplemented``. This closed the door to
  have the opportunity to try to recognize the dataset and map it to a
  supported CLASS. Now, when a CLASS attribute is not recognized, an attempt
  to recognize its associated dataset is made. If it is recognized, the
  matching class is associated with the dataset. If it is not recognized, then
  a warning is issued and the dataset becomes mapped to ``UnImplemented``.

- Always pass verbose and heavy values in the common test module to test().
  Fixes ticket #85.

- Now, the ``verbose`` and ``--heavy`` flag passed to test_all.py are honored.

- All the DLL's of dependencies are included now in Windows binaries.  This
  should allow for better portability of the binaries.

- Fixed the description of Node._v_objectID that was misleading.


Changes from 2.0rc1 to 2.0rc2
=============================

- The "Optimization tips" chapter of the User's Guide has been completely
  updated to adapt to PyTables 2.0 series.  In particular, new benchmarks on
  the much improved indexed queries have been included; you will see that
  PyTables indexing is competitive (and sometimes much faster) than that of
  traditional relational databases.  With this, the manual should be fairly
  finished for 2.0 final release.

- Large refactoring done on the ``Row`` class.  The most important change is
  that ``Table.row`` is now a single object.  This allows to reuse the same
  ``Row`` instance even after ``Table.flush()`` calls, which can be convenient
  in many situations.

- I/O buffers unified in the ``Row`` class.  That allows for bigger savings in
  memory space whenever the ``Row`` extension is used.

- Improved speed (up to a 70%) with unaligned column operations (a quite
  common scenario when dealing with ``Table`` objects) through the integrated
  Numexpr.  In-kernel searches take advantage of this optimization.

- Added ``VLUnicodeAtom`` for storing variable-length Unicode strings in
  ``VLArray`` objects regardless of encoding.  Closes ticket #51.

- Added support for ``time`` datatypes to be portable between big-endian and
  low-endian architectures.  This feature is not currently supported natively
  by the HDF5 library, so the support for such conversion has been added in
  PyTables itself.  Fixes #72.

- Added slice arguments to ``Table.readWhere()`` and ``Table.getWhereList()``.
  Although API changes are frozen, this may still be seen as an inconsistency
  with other query methods.  The patch is backwards-compatible anyway.

- Added missing overwrite argument to ``File.renameNode()`` and
  ``Node._f_rename()``.  Fixes ticket #66.

- Calling ``tables.test()`` no longer exits the interpreter session.  Fixes
  ticket #67.

- Fix comparing strings where one is a prefix of the other in integrated
  Numexpr.  Fixes ticket #76.

- Added a check for avoiding an ugly HDF5 message when copying a file over
  itself (for both ``copyFile()`` and ``File.copyFile()``).  Fixes ticket #73.

- Corrected the appendix E, were it was said that PyTables doesn't support
  compounds of compounds (it does since version 1.2!).


Changes from 2.0b2 to 2.0rc1
============================

- The ``lastrow`` argument of ``Table.flushRowsToIndex()`` is no longer
  public.  It was not documented, anyway.  Fixes ticket #43.

- Added a ``memlevel`` argument to ``Cols.createIndex()`` which allows the
  user to control the amount of memory required for creating an index.

- Added ``blocksizes`` and ``opts`` arguments to ``Cols.createIndex()``, which
  allow the user to control the sizes of index datasets, and to specify
  different optimization levels for each index dataset, respectively.  These
  are very low-level options meant only for experienced users.  Normal users
  should stick to the higher-level ``memlevel`` and ``optlevel``.

- Query tests have been tuned to exhaustively check the new parametrization of
  indexes.

- A new algorithm has been implemented that better reduces the entropy of
  indexes.

- The API Reference section of the User's Manual (and the matching docstrings)
  has been completely reviewed, expanded and corrected.  This process has
  unveiled some errors and inconsistencies which have also been fixed.

- Fixed ``VLArray.__getitem__()`` to behave as expected in Python when using
  slices, instead of following the semantics of PyTables' ``read()`` methods
  (e.g. reading just one element when no stop is provided).  Fixes ticket #50.

- Removed implicit UTF-8 encoding from ``VLArray`` data using ``vlstring``
  atoms.  Now a variable-length string is stored as is, which lets users use
  any encoding of their choice, or none of them.  A ``vlunicode`` atom will
  probably be added to the next release so as to fix ticket #51.

- Allow non-sequence objects to be passed to ``VLArray.append()`` when using
  an ``object`` atom.  This was already possible in 1.x but stopped working
  when the old append syntax was dropped in 2.0.  Fixes ticket #63.

- Changed ``Cols.__len__()`` to return the number of rows of the table or
  nested column (instead of the number of fields), like its counterparts in
  ``Table`` and ``Column``.

- Python scalars cached in ``AttributeSet`` instances are now kept as NumPy
  objects instead of Python ones, because they do become NumPy objects when
  retrieved from disk.  Fixes ticket #59.

- Avoid HDF5 error when appending an empty array to a ``Table`` (ticket #57)
  or ``EArray`` (ticket #49) dataset.

- Fix wrong implementation of the top-level ``table.description._v_dflts``
  map, which was also including the pathnames of columns inside nested
  columns.  Fixes ticket #45.

- Optimized the access to unaligned arrays in Numexpr between a 30% and a 70%.

- Fixed a die-hard bug that caused the loading of groups while closing a file.
  This only showed with certain usage patterns of the LRU cache (e.g. the one
  caused by ``ManyNodesTestCase`` in ``test_indexes.py`` under Pro).

- Avoid copious warnings about unused functions and variables when compiling
  Numexpr.

- Several fixes to Numexpr expressions with all constant values.  Fixed
  tickets #53, #54, #55, #58.  Reported bugs to mainstream developers.

- Solved an issue when trying to open one of the included test files in append
  mode on a system-wide installation by a normal user with no write privileges
  on it.  The file isn't being modified anyway, so the test is skipped then.

- Added a new benchmark to compare the I/O speed of ``Array`` and ``EArray``
  objects with that of ``cPickle``.

- The old ``Row.__call__()`` is no longer available as a public method.  It
  was not documented, anyway.  Fixes ticket #46.

- ``Cols._f_close()`` is no longer public.  Fixes ticket #47.

- ``Attributes._f_close()`` is no longer public.  Fixes ticket #52.

- The undocumented ``Description.classdict`` attribute has been completely
  removed.  Fixes ticket #44.


Changes from 2.0b1 to 2.0b2
===========================

- A very exhaustive overhauling of the User's Manual is in process.  The
  chapters 1 (Introduction), 2 (Installation), 3 (Tutorials) have been
  completed (and hopefully, the lines of code are easier to copy&paste now),
  while chapter 4 (API Reference) has been done up to (and including) the
  Table class.  During this tedious (but critical in a library) overhauling
  work, we have tried hard to synchronize the text in the User's Guide with
  that which appears on the docstrings.

- Removed the ``recursive`` argument in ``Group._f_walkNodes()``.  Using it
  with a false value was redundant with ``Group._f_iterNodes()``.  Fixes
  ticket #42.

- Removed the ``coords`` argument from ``Table.read()``.  It was undocumented
  and redundant with ``Table.readCoordinates()``.  Fixes ticket #41.

- Fixed the signature of ``Group.__iter__()`` (by removing its parameters).

- Added new ``Table.coldescrs`` and ``Table.description._v_itemsize``
  attributes.

- Added a couple of new attributes for leaves:

  * ``nrowsinbuf``: the number of rows that fit in the internal buffers.
  * ``chunkshape``: the chunk size for chunked datasets.

- Fixed setuptools so that making an egg out of the PyTables 2 package is
  possible now.

- Added a new ``tables.restrict_flavors()`` function allowing to restrict
  available flavors to a given set.  This can be useful e.g. if you want to
  force PyTables to get NumPy data out of an old, ``numarray``-flavored
  PyTables file even if the ``numarray`` package is installed.

- Fixed a bug which caused filters of unavailable compression libraries to be
  loaded as using the default Zlib library, after issuing a warning.  Added a
  new ``FiltersWarning`` and a ``Filters.copy()``.


Changes from 1.4.x to 2.0b1
===========================

API additions
-------------

- ``Column.createIndex()`` has received a couple of new parameters:
  ``optlevel`` and ``filters``.  The first one sets the desired quality level
  of the index, while the second one allows the user to specify the filters
  for the index.

- ``Table.indexprops`` has been split into ``Table.indexFilters`` and
  ``Table.autoIndex``.  The later groups the functionality of the old ``auto``
  and ``reindex``.

- The new ``Table.colpathnames`` is a sequence which contains the full
  pathnames of all bottom-level columns in a table.  This can be used to walk
  all ``Column`` objects in a table when used with ``Table.colinstances``.

- The new ``Table.colinstances`` dictionary maps column pathnames to their
  associated ``Column`` or ``Cols`` object for simple or nested columns,
  respectively.  This is similar to ``Table.cols._f_col()``, but faster.

- ``Row`` has received a new ``Row.fetch_all_fields()`` method in order to
  return all the fields in the current row.  This returns a NumPy void scalar
  for each call.

- New ``tables.test(verbose=False, heavy=False)`` high level function for
  interactively running the complete test suite from the Python console.

- Added a ``tables.print_versions()`` for easily getting the versions for all
  the software on which PyTables relies on.

Backward-incompatible changes
-----------------------------

- You can no longer mark a column for indexing in a ``Col`` declaration.  The
  only way of creating an index for a column is to invoke the
  ``createIndex()`` method of the proper column object *after the table has
  been created*.

- Now the ``Table.colnames`` attribute is just a list of the names of
  top-level columns in a table.  You can still get something similar to the
  old structure by using ``Table.description._v_nestedNames``.  See also the
  new ``Table.colpathnames`` attribute.

- The ``File.objects``, ``File.leaves`` and ``File.groups`` dictionaries have
  been removed.  If you still need this functionality, please use the
  ``File.getNode()`` and ``File.walkNodes()`` instead.

- ``Table.removeIndex()`` is no longer available; to remove an index on a
  column, one must use the ``removeIndex()`` method of the associated
  ``Column`` instance.

- ``Column.dirty`` is no longer available.  If you want to check
  column index dirtiness, use ``Column.index.dirty``.

- ``complib`` and ``complevel`` parameters have been removed from
  ``File.createTable()``, ``File.createEArray()``, ``File.createCArray()`` and
  ``File.createVLArray()``.  They were already deprecated in PyTables 1.x.

- The ``shape`` and ``atom`` parameters have been swapped in
  ``File.createCArray()``.  This has been done to be consistent with
  ``Atom()`` definitions (i.e. type comes before and shape after).

Deprecated features
-------------------

- ``Node._v_rootgroup`` has been removed.  Please use ``node._v_file.root``
  instead.

- The ``Node._f_isOpen()`` and ``Leaf.isOpen()`` methods have been removed.
  Please use the ``Node._v_isopen`` attribute instead (it is much faster).

- The ``File.getAttrNode()``, ``File.setAttrNode()`` and
  ``File.delAttrNode()`` methods have been removed.  Please use
  ``File.getNodeAttr()``, ``File.setNodeAttr()`` and ``File.delNodeAttr()``
  instead.

- ``File.copyAttrs()`` has been removed.  Please use ``File.copyNodeAttrs()``
  instead.

- The ``table[colname]`` idiom is no longer supported.  You can use
  ``table.cols._f_col(column)`` for doing the same.

API refinements
---------------

- ``File.createEArray()`` received a new ``shape`` parameter.  This allows to
  not have to use the shape of the atom so as to set the shape of the
  underlying dataset on disk.

- All the leaf constructors have received a new ``chunkshape`` parameter that
  allows specifying the chunk sizes of datasets on disk.

- All ``File.create*()`` factories for ``Leaf`` nodes have received a new
  ``byteorder`` parameter that allows the user to specify the byteorder in
  which data will be written to disk (data in memory is now always handled in
  *native* order).


----

  **Enjoy data!**

  -- The PyTables Team


.. Local Variables:
.. mode: rst
.. coding: utf-8
.. fill-column: 78
.. End:
:author: valhallasw
:date: 2012-06-18 10:15:15

===================
Hints for SQL users
===================

This page is intended to be **a guide to new PyTables for users who are used
to writing SQL code** to access their relational databases.
It will cover the most usual SQL statements.
If you are missing a particular statement or usage example, you can ask at the
`PyTables users' list`_ for it.
If you know some examples yourself, you can also write them here!

This page is under development: you can come back frequently to check for new
examples.
Also, this is no replacement for the `User's Guide`_;
if you don't read the manual, you'll be missing lots of features not available
in relational databases!

Examples in Python assume that you have imported the PyTables package like
this::

    import tables

.. .. contents:: Table Of Contents


Creating a new database
=======================

RDBMs happen to have several syntaxes for creating a database.
A usual syntax is::

    CREATE DATABASE database_name

In PyTables, each database goes to a different HDF5_ file (much like
SQLite_ or MS Access).
To create a new HDF5_ file, you use the :func:`tables.open_file` function with
the ``'w'`` mode (which deletes the database if it already exists), like this::

    h5f = tables.open_file('database_name.h5', 'w')

In this way you get the ``h5f`` PyTables file handle (an instance of the
:class:`tables.File` class), which is a concept similar to a *database
connection*, and a new :file:`database_name.h5` file is created in the current
directory (you can use full paths here).
You can close the handle (like you close the connection) with::

    h5f.close()

This is important for PyTables to dump pending changes to the database.
In case you forget to do it, PyTables closes all open database handles for
you when you exit your program or interactive session, but it is always safer
to close your files explicitly.
If you want to use the database after closing it, you just call
:func:`open_file` again, but using the ``'r+'`` or ``'r'`` modes, depending on
whether you do or don't need to modify the database, respectively.

You may use several PyTables databases simultaneously in a program, so you
must be explicit on which database you want to act upon (by using its handle).

A note on concurrency under PyTables
------------------------------------

Unlike most RDBMs, PyTables is not intended to serve concurrent accesses to a
database.
It has no protections whatsoever against corruption for different (or even the
same) programs accessing the same database.
Opening several handles to the same database in read-only mode is safe, though.


Creating a table
================

PyTables supports some other *datasets* besides tables, and they're not
arranged in a flat namespace, but rather into a *hierarchical* one (see an
introduction to the _ref:`object tree <ObjectTreeSection>`);
however, due to the nature of these recipes, we'll limit ourselves to tables
in the *root group*.
The basic syntax for table creation under SQL is::

    CREATE TABLE table_name (
        column_name1 column_type1,
        column_name2 column_type2,
        ...
        column_nameN column_typeN
    )


Table descriptions
------------------

In PyTables, one first *describes* the structure of a table.
PyTables allows you to *reuse a description* for creating several tables with
the same structure, just by using the description object (``description_name``
below) or getting it from a created table.
This is specially useful for creating temporary tables holding query results.

You can create a table description using a dictionary::

    description_name = {
        'column_name1': column_type1,
        'column_name2': column_type2,
        'column_name3': column_type3,
        ...
        'column_nameN': column_typeN
    }

or a subclass of :class:`tables.IsDescription`::

    class description_name(tables.IsDescription):
        column_name1 = column_type1
        column_name2 = column_type2
        column_name3 = column_type3
        ...
        column_nameN = column_typeN

Please note that dictionaries are the only way of describing structures with
names which cannot be Python identifiers.
Also, if an explicit order is desired for columns, it must be specified through
the column type declarations (see below), since dictionary keys and class
attributes aren't ordered.
Otherwise, columns are ordered in alphabetic increasing order.
It is important to note that PyTables doesn't have a concept of primary or
foreign keys, so relationships between tables are left to the user.


Column type declarations
------------------------

PyTables supports lots of types (including nested and multidimensional
columns).
Non-nested columns are declared through instances of :class:`tables.Col`
subclasses (which you can also reuse).
These are some correspondences with SQL:

==================== ==========================
SQL type declaration PyTables type declaration
==================== ==========================
INTEGER(digits)      tables.IntCol(itemsize)
REAL                 tables.FloatCol()
VARCHAR(length)      tables.StringCol(itemsize)
DATE                 tables.Time32Col()
TIMESTAMP            tables.Time64Col()
==================== ==========================

See a complete description of :ref:`PyTables types <datatypes>`.
Note that some types admit different *item sizes*, which are specified in
bytes.
For types with a limited set of supported item sizes, you may also use specific
subclasses which are named after the type and its *precision*, e.g. ``Int32Col``
for 4-byte (32 bit) item size.

Cells in a PyTables' table always have a value of the cell type, so there is
no ``NULL``.
Instead, cells take a *default value* (zero or empty) which can be changed in
the type declaration, like this: ``col_name = StringCol(10, dflt='nothing')``
(``col_name`` takes the value ``'nothing'`` if unset).
The declaration also allows you to set *column order* via the ``pos`` argument,
like this::

    class ParticleDescription(tables.IsDescription):
        name = tables.StringCol(10, pos=1)
        x = tables.FloatCol(pos=2)
        y = tables.FloatCol(pos=3)
        temperature = tables.FloatCol(pos=4)


Using a description
===================

Once you have a table description ``description_name`` and a writeable file
handle ``h5f``, creating a table with that description is as easy as::

    tbl = h5f.create_table('/', 'table_name', description_name)

PyTables is very object-oriented, and database is usually done through
methods of :class:`tables.File`.
The first argument indicates the *path* where the table will be created,
i.e. the root path (HDF5 uses Unix-like paths).
The :meth:`tables.File.create_table` method has many options e.g. for setting
a table title or compression properties. What you get back is an instance of
:class:`tables.Table`, a handle for accessing the data in that table.

As with files, table handles can also be closed with ``tbl.close()``.
If you want to access an already created table, you can use::

    tbl = h5f.get_node('/', 'table_name')

(PyTables uses the concept of *node* for datasets -tables and others- and
groups in the object tree) or, using *natural naming*::

    tbl = h5f.root.table_name

Once you have created a table, you can access (and reuse) its description by
accessing the ``description`` attribute of its handle.


Creating an index
=================

RDBMs use to allow named indexes on any set of columns (or all of them) in a
table, using a syntax like::

    CREATE INDEX index_name
    ON table_name (column_name1, column_name2, column_name3...)

and

    DROP INDEX index_name

Indexing is supported in the versions of PyTables >= 2.3 (and in PyTablesPro).
However, indexes don't have names and they are bound to single columns.
Following the object-oriented philosophy of PyTables, index creation is a
method (:meth:`tables.Column.create_index`) of a :class:`tables.Column` object
of a table, which you can access trough its ``cols`` accessor.

::
    tbl.cols.column_name.create_index()

For dropping an index on a column::

    tbl.cols.column_name.remove_index()


Altering a table
================

The first case of table alteration is renaming::

    ALTER TABLE old_name RENAME TO new_name

This is accomplished in !PyTables with::

    h5f.rename_node('/', name='old_name', newname='new_name')

or through the table handle::

    tbl.rename('new_name')

A handle to a table is still usable after renaming.
The second alteration, namely column addition, is currently not supported in
PyTables.


Dropping a table
================

In SQL you can remove a table using::

    DROP TABLE table_name

In PyTables, tables are removed as other nodes, using the
:meth:`tables.File.remove_node` method::

    h5f.remove_node('/', 'table_name')

or through the table handle::

    tbl.remove()

When you remove a table, its associated indexes are automatically removed.


Inserting data
==============

In SQL you can insert data one row at a time (fetching from a selection will
be covered later) using a syntax like::

    INSERT INTO table_name (column_name1, column_name2...)
    VALUES (value1, value2...)

In PyTables, rows in a table form a *sequence*, so data isn't *inserted* into
a set, but rather *appended* to the end of the sequence.
This also implies that identical rows may exist in a table (but they have a
different *row number*).
There are two ways of appending rows: one at a time or in a block.
The first one is conceptually similar to the SQL case::

    tbl.row['column_name1'] = value1
    tbl.row['column_name2'] = value2
    ...
    tbl.row.append()

The ``tbl.row`` accessor represents a *new row* in the table.
You just set the values you want to set (the others take the default value
from their column declarations - see above) and the effectively append the
new row.
This code is usually enclosed in some kind of loop, like::

    row = tbl.row
    while some_condition:
        row['column_name1'] = value1
        ...
        row.append()

For appending a block of rows in a single shot, :meth:`tables.Table.append`
is more adequate.
You just pass a NumPy_ record array or Python sequence with elements which
match the expected columns.
For example, given the ``tbl`` handle for a table with the ``ParticleDescription``
structure described above::

    rows = [
        ('foo', 0.0, 0.0, 150.0),
        ('bar', 0.5, 0.0, 100.0),
        ('foo', 1.0, 1.0,  25.0)
    ]
    tbl.append(rows)

    # Using a NumPy container.
    import numpy
    rows = numpy.rec.array(rows)
    tbl.append(rows)


A note on transactions
----------------------

PyTables doesn't support transactions nor checkpointing or rolling back (there
is undo support for operations performed on the object tree, but this is
unrelated).
Changes to the database are optimised for maximum performance and reasonable
memory requirements, which means that you can't tell whether e.g.
``tbl.append()`` has actually committed all, some or no data to disk when it ends.

However, you can *force* PyTables to commit changes to disk using the ``flush()``
method of table and file handles::

    tbl.flush()  # flush data in the table
    h5f.flush()  # flush all pending data

Closing a table or a database actually flushes it, but it is recommended that
you explicitly flush frequently (specially with tables).


Updating data
=============

We're now looking for alternatives to the SQL ``UPDATE`` statement::

    UPDATE table_name
    SET column_name1 = expression1, column_name2 = expression2...
    [WHERE condition]

There are different ways of approaching this, depending on your needs.
If you aren't using a condition, then the ``SET`` clause updates all rows,
something you can do in PyTables by iterating over the table::

    for row in tbl:
        row['column_name1'] = expression1
        row['column_name2'] = expression2
        ...
        row.update()

Don't forget to call ``update()`` or no value will be changed!
Also, since the used iterator allows you to read values from the current row,
you can implement a simple *conditional update*, like this::

    for row in tbl:
        if condition on row['column_name1'], row['column_name2']...:
            row['column_name1'] = expression1
            row['column_name2'] = expression2
            ...
            row.update()

There are substantially more efficient ways of locating rows fulfilling a
condition.
Given the main PyTables usage scenarios, querying and modifying data are
quite decoupled operations, so we will have a look at querying later and
assume that you already know the set of rows you want to update.

If the set happens to be a slice of the table, you may use the
:`meth:`tables.Table.modify_rows` method or its equivalent
:meth:`tables.Table.__setitem__` notation::

    rows = [
        ('foo', 0.0, 0.0, 150.0),
        ('bar', 0.5, 0.0, 100.0),
        ('foo', 1.0, 1.0,  25.0)
    ]
    tbl.modifyRows(start=6, stop=13, step=3, rows=rows)
    tbl[6:13:3] = rows  # this is the same

If you just want to update some columns in the slice, use the
:meth:`tables.Table.modify_columns` or :meth:`tables.Table.modify_column`
methods::

    cols = [
        [150.0, 100.0, 25.0]
    ]
    # These are all equivalent.
    tbl.modify_columns(start=6, stop=13, step=3, columns=cols, names=['temperature'])
    tbl.modify_column(start=6, stop=13, step=3, column=cols[0], colname='temperature')
    tbl.cols.temperature[6:13:3] = cols[0]

The last line shows an example of using the ``cols`` accessor to get to the
desired :class:`tables.Column` of the table using natural naming and apply
``setitem`` on it.

If the set happens to be an array of sparse coordinates, you can also use
PyTables' extended slice notation::

    rows = [
        ('foo', 0.0, 0.0, 150.0),
        ('bar', 0.5, 0.0, 100.0),
        ('foo', 1.0, 1.0,  25.0)
    ]
    rownos = [2, 735, 371913476]
    tbl[rownos] = rows


instead of the traditional::

    for row_id, datum in zip(rownos, rows):
         tbl[row_id] = datum

Since you are modifying table data in all cases, you should also remember to
``flush()`` the table when you're done.


Deleting data
=============

Rows are deleted from a table with the following SQL syntax::

    DELETE FROM table_name
    [WHERE condition]

:meth:`tables.Table.remove_rows` is the method used for deleting rows in
PyTables.
However, it is very simple (only contiguous blocks of rows can be deleted) and
quite inefficient, and one should consider whether *dumping filtered data from
one table into another* isn't a much more convenient approach.
This is a far more optimized operation under PyTables which will be covered
later.

Anyway, using ``remove_row()`` or ``remove_rows()`` is quite straightforward::

    tbl.remove_row(12)  # delete one single row (12)
    tbl.remove_rows(12, 20)  # delete all rows from 12 to 19 (included)
    tbl.remove_rows(0, tbl.nrows)  # delete all rows unconditionally
    tbl.remove_rows(-4, tbl.nrows)  # delete the last 4 rows


Reading data
============

The most basic syntax in SQL for reading rows in a table without using a
condition is::

    SELECT (column_name1, column_name2... | *) FROM table_name

Which reads all rows (though maybe not all columns) from a table.
In PyTables there are two ways of retrieving data: *iteratively* or *at once*.
You'll notice some similarities with how we appended and updated data above,
since this dichotomy is widespread here.

For a clearer separation with conditional queries (covered further below),
and since the concept of *row number* doesn't exist in relational databases,
we'll be including here the cases where you want to read a **known** *slice*
or *sequence* of rows, besides the case of reading *all* rows.


Iterating over rows
-------------------

This is similar to using the ``fetchone()`` method of a DB ``cursor`` in a
`Python DBAPI`_-compliant package, i.e. you *iterate* over the list of wanted
rows, getting one *row handle* at a time.
In this case, the handle is an instance of the :class:`tables.Row` class,
which allows access to individual columns as items accessed by key (so there
is no special way of selecting columns: you just use the ones you want
whenever you want).

This way of reading rows is recommended when you want to perform operations
on individual rows in a simple manner, and specially if you want to process
a lot of rows in the table (i.e. when loading them all at once would take too
much memory).
Iterators are also handy for using with the ``itertools`` Python module for
grouping, sorting and other operations.

For iterating over *all* rows, use plain iteration or the
:meth:`tables.Table.iterrows` method::

    for row in tbl:  # or tbl.iterrows()
        do something with row['column_name1'], row['column_name2']...

For iterating over a *slice* of rows, use the
:meth:`tables.Table.iterrows|Table.iterrows` method::

    for row in tbl.iterrows(start=6, stop=13, step=3):
        do something with row['column_name1'], row['column_name2']...

For iterating over a *sequence* of rows, use the
:meth:`tables.Table.itersequence` method::

    for row in tbl.itersequence([6, 7, 9, 11]):
        do something with row['column_name1'], row['column_name2']...

Reading rows at once
--------------------

In contrast with iteration, you can fetch all desired rows into a single
*container* in memory (usually an efficient NumPy_ record-array) in a single
operation, like the ``fetchall()`` or ``fetchmany()`` methods of a DBAPI ``cursor``.
This is specially useful when you want to transfer the read data to another
component in your program, avoiding loops to construct your own containers.
However, you should be careful about the amount of data you are fetching into
memory, since it can be quite large (and even exceed its physical capacity).

You can choose between the ``Table.read*()`` methods or the
:meth:`tables.Table.__getitem__` syntax for this kind of reads.
The ``read*()`` methods offer you the chance to choose a single column to read
via their ``field`` argument (which isn't still as powerful as the SQL ``SELECT``
column spec).

For reading *all* rows, use ``[:]`` or the :meth:`tables.Table.read` method::

    rows = tbl.read()
    rows = tbl[:]  # equivalent

For reading a *slice* of rows, use ``[slice]`` or the
:meth:`tables.Table.read|Table.read` method::

    rows = tbl.read(start=6, stop=13, step=3)
    rows = tbl[6:13:3]  # equivalent

For reading a *sequence* of rows, use the :meth:`tables.Table.read_coordinates`
method::

    rows = tbl.read_coordinates([6, 7, 9, 11])

Please note that you can add a ``field='column_name'`` argument to ``read*()``
methods in order to get only the given column instead of them all.


Selecting data
==============

When you want to read a subset of rows which match a given condition from a
table you use a syntax like this in SQL::

    SELECT column_specification FROM table_name
    WHERE condition

The ``condition`` is an expression yielding a boolean value based on a
combination of column names and constants with functions and operators.
If the condition holds true for a given row, the ``column_specification`` is
applied on it and the resulting row is added to the result.

In PyTables, you may filter rows using two approaches: the first one is
achieved through standard Python comparisons (similar to what we used for
conditional update), like this::

    for row in tbl:
        if condition on row['column_name1'], row['column_name2']...:
            do something with row

This is easy for newcomers, but not very efficient. That's why PyTables offers
another approach: **in-kernel** searches, which are much more efficient than
standard searches, and can take advantage of indexing (under PyTables >= 2.3).

In-kernel searches are used through the *where methods* in ``Table``, which are
passed a *condition string* describing the condition in a Python-like syntax.
For instance, with the ``ParticleDescription`` we defined above, we may specify
a condition for selecting particles at most 1 unit apart from the origin with
a temperature under 100 with a condition string like this one::

    '(sqrt(x**2 + y**2) <= 1) & (temperature < 100)'

Where ``x``, ``y`` and ``temperature`` are the names of columns in the table.
The operators and functions you may use in a condition string are described
in the :ref:`appendix on condition syntax <condition_syntax>` in the
`User's Guide`_.


Iterating over selected rows
----------------------------

You can iterate over the rows in a table which fulfill a condition (a la DBAPI
``fetchone()``) by using the :meth:`tables.Table.where` method, which is very
similar to the :meth:`tables.Table.iterrows` one discussed above, and which
can be used in the same circumstances (i.e. performing operations on individual
rows or having results exceeding available memory).

Here is an example of using ``where()`` with the previous example condition::

    for row in tbl.where('(sqrt(x**2 + y**2) <= 1) & (temperature < 100)'):
        do something with row['name'], row['x']...


Reading selected rows at once
-----------------------------

Like the aforementioned :meth:`tables.Table.read`,
:meth:`tables.Table.read_where` gets all the rows fulfilling the given
condition and packs them in a single container (a la DBAPI ``fetchmany()``).
The same warning applies: be careful on how many rows you expect to retrieve,
or you may run out of memory!

Here is an example of using ``read_where()`` with the previous example
condition::

    rows = tbl.read_where('(sqrt(x**2 + y**2) <= 1) & (temperature < 100)')

Please note that both :meth:`tables.Table.where` and
:meth:`tables.Table.read_where` can also take slicing arguments.


Getting the coordinates of selected rows
----------------------------------------

There is yet another method for querying tables:
:meth:`tables.Table.get_where_list`.
It returns just a sequence of the numbers of the rows which fulfil the given
condition.
You may pass that sequence to :meth:`tables.Table.read_coordinates`, e.g. to
retrieve data from a different table where rows with the same number as the
queried one refer to the same first-class object or entity.


A note on table joins
---------------------

You may have noticed that queries in PyTables only cover one table.
In fact, there is no way of directly performing a join between two tables in
PyTables (remember that it's not a relational database).
You may however work around this limitation depending on your case:

* If one table is an *extension* of another (i.e. it contains additional
  columns for the same entities), your best bet is to arrange rows of the
  same entity so that they are placed in the same positions in both tables.
  For instance, if ``tbl1`` and ``tbl2`` follow this rule, you may do something
  like this to emulate a natural join::

    for row1 in tbl1.where('condition'):
        row2 = tbl2[row1.nrow]
        if condition on row2['column_name1'], row2['column_name2']...:
            do something with row1 and row2...

   (Note that ``row1`` is a ``Row`` instance and ``row2`` is a record of the current
   flavor.)

* If rows in both tables are linked by a common value (e.g. acting as an
  identifier), you'll need to split your condition in one for the first table
  and one for the second table, and then nest your queries, placing the most
  restrictive one first. For instance::

    SELECT clients.name, bills.item_id FROM clients, bills
    WHERE clients.id = bills.client_id and clients.age > 50 and bills.price > 200

  could be written as::

    for client in clients.where('age > 50'):
        # Note that the following query is different for each client.
        for bill in bills.where('(client_id == %r) & (price > 200)' % client['id']):
            do something with client['name'] and bill['item_id']

  In this example, indexing the ``client_id`` column of ``bills`` could speed up
  the inner query quite a lot.
  Also, you could avoid parsing the inner condition each time by using
  *condition variables*::

    for client in clients.where('age > 50'):
        for bill in bills.where('(client_id == cid) & (price > 200)', {'cid': client['id']}):
            do something with client['name'] and bill['item_id']


Summary of row selection methods
================================

+----------------------+-----------------+---------------------+-----------------------+-------------------------+
|                      | **All rows**    | **Range of rows**   | **Sequence of rows**  | **Condition**           |
+----------------------+-----------------+---------------------+-----------------------+-------------------------+
| **Iterative access** | ``__iter__()``, | ``iterrows(range)`` | ``itersequence()``    | ``where(condition)``    |
|                      | ``iterrows()``  |                     |                       |                         |
+----------------------+-----------------+---------------------+-----------------------+-------------------------+
| **Block access**     | ``[:]``,        | ``[range]``,        | ``readCoordinates()`` |``read_where(condition)``|
|                      | ``read()``      | ``read(range)``     |                       |                         |
+----------------------+-----------------+---------------------+-----------------------+-------------------------+


Sorting the results of a selection
==================================

*Do you feel like writing this section? Your contribution is welcome!*


Grouping the results of a selection
===================================

By making use of the :func:`itertools.groupby` utility, you can group results
by field::

    group = {} # dictionary to put results grouped by 'pressure'
    def pressure_selector(row):
        return row['pressure']
    for pressure, rows_grouped_by_pressure in itertools.groupby(mytable, pressure_selector):
        group[pressure] = sum((r['energy'] + r['ADCcount'] for r in rows_grouped_by_pressure))

However, :func:`itertools.groupby` assumes the incoming array is sorted by the
grouping field.
If not, there are multiple groups with the same grouper returned.
In the example, mytable thus has to be sorted on pressure, or the last line
should be changed to::

    group[pressure] += sum((r['energy'] + r['ADCcount'] for r in rows_grouped_by_pressure))


-----


.. target-notes::

.. _`PyTables users' list`: https://lists.sourceforge.net/lists/listinfo/pytables-users
.. _`User's Guide`: https://www.pytables.org/usersguide
.. _HDF5: http://www.hdfgroup.org/HDF5
.. _SQLite: http://www.sqlite.org
.. _NumPy: http://www.numpy.org
.. _`Python DBAPI`: http://www.python.org/dev/peps/pep-0249
====================
In-memory HDF5 files
====================

The HDF5 library provides functions to allow an application to work with a
file in memory for faster reads and writes. File contents are kept in memory
until the file is closed.  At closing, the memory version of the file can be
written back to disk or abandoned.


Open an existing file in memory
===============================

Assuming the :file:`sample.h5` exists in the current folder, it is possible to
open it in memory simply using the CORE driver at opening time.

The HDF5 driver that one intend to use to open/create a file can be specified
using the *driver* keyword argument of the :func:`tables.open_file` function::

    >>> import tables
    >>> h5file = tables.open_file("sample.h5", driver="H5FD_CORE")

The content of the :file`sample.h5` is opened for reading. It is loaded into
memory and all reading operations are performed without disk I/O overhead.

.. note::

    the initial loading of the entire file into memory can be time expensive
    depending on the size of the opened file and on the performances of the
    disk subsystem.

.. seealso::

    general information about HDF5 drivers can be found in the `Alternate
    File Storage Layouts and Low-level File Drivers`__ section of the `HDF5
    User's Guide`_.

__ `HDF5 drivers`_


Creating a new file in memory
=============================

Creating a new file in memory is as simple as creating a regular file, just
one needs to specify to use the CORE driver::

    >>> import tables
    >>> h5file = tables.open_file("new_sample.h5", "w", driver="H5FD_CORE")
    >>> import numpy
    >>> a = h5file.create_array(h5file.root, "array", numpy.zeros((300, 300)))
    >>> h5file.close()


Backing store
=============

In the previous example contents of the in-memory `h5file` are automatically
saved to disk when the file descriptor is closed, so a new
:file:`new_sample.h5` file is created and all data are transferred to disk.

Again this can be time a time expensive action depending on the amount of
data in the HDF5 file and depending on how fast the disk I/O is.

Saving data to disk is the default behavior for the CORE driver in PyTables.

This feature can be controlled using the *driver_core_backing_store*
parameter of the :func:`tables.open_file` function.  Setting it to `False`
disables the backing store feature and all changes in the working `h5file`
are lost after closing::

    >>> h5file = tables.open_file("new_sample.h5", "w", driver="H5FD_CORE",
    ...                           driver_core_backing_store=0)

Please note that the *driver_core_backing_store* disables saving of data, not
loading.
In the following example the :file:`sample.h5` file is opened in-memory in
append mode.  All data in the existing :file:`sample.h5` file are loaded into
memory and contents can be actually modified by the user::

    >>> import tables
    >>> h5file = tables.open_file("sample.h5", "a", driver="H5FD_CORE",
                                  driver_core_backing_store=0)
    >>> import numpy
    >>> h5file.create_array(h5file.root, "new_array", numpy.arange(20),
                            title="New array")
    >>> array2 = h5file.root.array2
    >>> print(array2)
    /array2 (Array(20,)) 'New array'
    >>> h5file.close()

Modifications are lost when the `h5file` descriptor is closed.


Memory images of HDF5 files
===========================

It is possible to get a memory image of an HDF5 file (see
`HDF5 File Image Operations`_).  This feature is only available if PyTables
is build against version 1.8.9 or newer of the HDF5 library.

In particular getting a memory image of an HDF5 file is possible only if the
file has been opened with one of the following drivers: SEC2 (the default
one), STDIO or CORE.

An example of how to get an image::

    >>> import tables
    >>> h5file = tables.open_file("sample.h5")
    >>> image = h5file.get_file_image()
    >>> h5file.close()

The memory ìmage of the :file:`sample.h5` file is copied into the `ìmage`
string (of bytes).

.. note::

    the `ìmage` string contains all data stored in the HDF5 file so, of
    course, it can be quite large.

The `ìmage` string can be passed around and can also be used to initialize a
new HDF5 file descriptor::

    >>> import tables
    >>> h5file = tables.open_file("in-memory-sample.h5", driver="H5FD_CORE",
                                  driver_core_image=image,
                                  driver_core_backing_store=0)
    >>> print(h5file.root.array)
    /array (Array(300, 300)) 'Array'
    >>> h5file.setNodeAttr(h5file.root, "description", "In memory file example")



-----


.. target-notes::

.. _`HDF5 drivers`: http://www.hdfgroup.org/HDF5/doc/UG/08_TheFile.html#Drivers
.. _`HDF5 User's Guide`: https://portal.hdfgroup.org/display/HDF5/HDF5+User+Guides
.. _`HDF5 File Image Operations`: http://www.hdfgroup.org/HDF5/doc/Advanced/FileImageOperations/HDF5FileImageOperations.pdf
:author: localhost
:date: 2008-04-21 11:12:44

.. todo:: update to use new SW versions


Installing PyTables when you're not root
========================================

By `Koen van de Sande <http://www.tibed.net>`_.

.. warning:: contents of this recipe may be outdated.

This guide describes how to install PyTables and its dependencies on Linux or
other \*nix systems when your user account is not root.
Installing the HDF5_ shared libraries and Python extension
NumPy requires some non-trivial steps to work.
We describe all steps needed.
They only assumption is that you have Python 3.6 or higher and a C/C++
compiler (gcc) installed.


Installing HDF5
---------------

* First go to or make a temporary folder where we can download and compile
  software.
  We'll assume you're in this temporary folder in the rest of this section.
* Download `hdf5-1.12.1.tar.gz` from https://www.hdfgroup.org/downloads/hdf5
* Extract the archive to the current folder::

    tar xzvf hdf5-1.12.1.tar.gz

* Go to the extracted HDF5 folder::

    cd hdf5-1.12.1

* Run the configure script::

    ./configure

* Run make::

    make install

* We've now compiled HDF5_ into the `hdf5` folder inside the source tree.
  We'll need to move this to its final location.
  For this guide, we'll make a `software` folder inside your home directory
  to store installed libraries::

    mkdir ~/software

* Move the files to the right location::

    mv hdf5 ~/software/


Installing NumPy
----------------

* From the `NumPy page on PyPI <https://pypi.org/project/numpy/>`_
  download NumPy 1.21.5 (at time of writing) to our temporary folder.
* Extract the archive::

    tar xzvf numpy-1.21.5.tar.gz

* Go to the NumPy folder::

    cd numpy-1.21.5
* Build and install the Python module into our software folder::

    python3 setup.py install --home=~/software


Python wrapper script
---------------------

We've installed all dependencies of PyTables.
We need to create a wrapper script for Python to let PyTables actually find
all these dependencies.
Had we installed them as root, they'd be trivial to find, but now we need to
help a bit.

* Create a script with the following contents (I've called this script `p` on
  my machine)::

    #!/bin/bash
    export PYTHONPATH=~/software/lib/python
    export HDF5_DIR=~/software/hdf5
    export LD_LIBRARY_PATH=~/software/lib/python/tables:~/software/hdf5/lib
    python3 $*

* Make the script executable::

    chmod 755 p

* Place the script somewhere on your path (for example, inside a folder
  called `bin` inside your home dir, which is normally added to the path
  automatically).
  If you do not add this script to your path, you'll have to replace `p` in
  scripts below by the full path (and name of) your script, e.g.
  `~/pytablespython.sh` if you called it `pytablespython.sh` and put it in
  your home dir.
* Test your Python wrapper script::

    p

* It should now start Python. And you should be able to import `numpy`
  without errors::

    >>> import numpy


.. note::

    you could do this differently by defining these environment settings
    somewhere in your startup scripts, but this wrapper script approach is
    cleaner.


Installing PyTables
-------------------

* From the `PyPI page <https://pypi.org/project/tables/>`_
  download PyTables 3.7.0 (at time of writing) to our temporary folder.
* Extract the archive::

    tar xzvf pytables-3.7.0.tar.gz

* Go to the PyTables folder::

    cd pytables-3.7.0

* Install PyTables using our wrapper script::

    p setup.py install --home=~/software


Running Python with PyTables support
------------------------------------

* Use your Python wrapper script to start Python::

    p

* You can now import `tables` without errors::

    >>> import tables
    >>> tables.__version__
    '3.7.0'


Concluding remarks
------------------

* It is safe to remove the temporary folder we have used in this guide,
  there are no dependencies on it.
* This guide was written for and tested with HDF5 1.12.1, PyTables 3.7.6 and
  Numpy 1.21.5.


Enjoy working with PyTables!

*Koen*


-----


.. target-notes::

.. _HDF5: http://www.hdfgroup.org/HDF5
:author: KennethArnold
:date: 2009-07-14 21:51:07


================================
Using your own custom data types
================================

You can make your own data types by subclassing Table (or other PyTables types,
such as :class:`tables.Leaf`).
This can be useful for storing a specialized type of data or presenting a
customized API.

Submitted by Kevin R. Thornton.

::

    import numpy as np

    import tables
    from tables import File, Table
    from tables.file import _checkfilters

    from tables.parameters import EXPECTED_ROWS_TABLE


    class DerivedFromTable(Table):
        _c_classId = 'DerivedFromTable'

        def __init__(self, parentNode, name, description=None,
                     title="", filters=None,
                     expectedrows=EXPECTED_ROWS_TABLE,
                     chunkshape=None, byteorder=None, _log=True):
            super().__init__(parentNode, name,
                             description=description, title=title,
                             filters=filters,
                             expectedrows=expectedrows,
                             chunkshape=chunkshape, byteorder=byteorder,
                             _log=_log)

        def read(self, start=None, stop=None, step=None, field=None):
            print("HERE!")
            data = Table.read(self, start=start, stop=stop, step=step,
                              field=field)
            return data


    def createDerivedFromTable(self, where, name, data, title="",
                               filters=None, expectedrows=10000,
                               chunkshape=None, byteorder=None,
                               createparents=False):
        parentNode = self._get_or_create_path(where, createparents)

        _checkfilters(filters)
        return DerivedFromTable(parentNode, name, data,
                                title=title, filters=filters,
                                expectedrows=expectedrows,
                                chunkshape=chunkshape, byteorder=byteorder)


    File.createDerivedFromTable = createDerivedFromTable


    if __name__ == '__main__':
        x = np.random.rand(100).reshape(50,2)
        x.dtype = [('x',float), ('y',float)]
        h5file = tables.open_file('tester.hdf5', 'w')
        mtab = h5file.createDerivedFromTable(h5file.root, 'random', x)

        h5file.flush()
        print(type(mtab))
        mtab_read = mtab.read()
        h5file.close()
        h5file = tables.open_file('tester.hdf5', 'r')
        mtab = h5file.root.random

        print(type(mtab))
        mtab_read2 = mtab.read()
        print(np.array_equal(mtab_read, mtab_read2))


There is an issue that the DerivedFromTable read function will not be called
when the file is re-opened. The notion that the H5 file contains a derived
object gets lost. The output shows that the read function is only called before
the function is closed:

::
        <class '__main__.DerivedFromTable'>
        HERE!
        <class 'tables.table.Table'>
        True
        Closing remaining open files:tester.hdf5...done


I ran into this because I wanted a custom read that returned a more complex
object implemented in C++. Using pybind11, I'm easily able to write to a
Table via a record array. I was hoping that I could read back in, construct the
correct C++-based type, and return it. The example seems to suggest that this
is not possible.
:author: AskJakobsen
:date: 2012-06-01 09:31:14

========================
Tailoring `atexit` hooks
========================

In some situations you may want to tailor the typical messages that PyTables
outputs::

    Closing remaining open files: /tmp/prova.h5... done

The responsible of this behaviour is the :func:`tables.file.close_open_files`
function that is being registered via :func:`atexit.register` Python function.
Although you can't de-register already registered cleanup functions, you can
register new ones to tailor the existing behaviour.
For example, if you  register this function::

    def my_close_open_files(verbose):
        from packaging import Version

        open_files = tables.file._open_files

        are_open_files = len(open_files) > 0

        if verbose and are_open_files:
            sys.stderr.write("Closing remaining open files:")

        if Version(tables.__version__) >= Version("3.1.0"):
            # make a copy of the open_files.handlers container for the iteration
            handlers = list(open_files.handlers)
        else:
            # for older versions of pytables, setup the handlers list from the
            # keys
            keys = open_files.keys()
            handlers = []
            for key in keys:
                handlers.append(open_files[key])

        for fileh in handlers:
            if verbose:
                sys.stderr.write("%s..." % fileh.filename)

            fileh.close()

            if verbose:
                sys.stderr.write("done")

        if verbose and are_open_files:
            sys.stderr.write("\n")

    import sys, atexit
    atexit.register(my_close_open_files, False)

then, you won't get the closing messages anymore because the new registered
function is executed before the existing one.
If you want the messages back again, just set the verbose parameter to true.

You can also use the `atexit` hooks to perform other cleanup functions as well.

:author: localhost
:date: 2008-04-21 11:12:45

.. todo:: update the code example to numpy

=============================================================
How to integrate PyTables in your application by using py2exe
=============================================================

This document shortly describes how to build an executable when using PyTables.
Py2exe_ is a third party product that converts python scripts into standalone
windows application/programs.
For more information about py2exe please visit http://www.py2exe.org.

To be able to use py2exe you have to download and install it.
Please follow the instructions at http://www.py2exe.org.

Let’s assume that you have written a python script as in the attachment
:download:`py2exe_howto/pytables_test.py`

.. literalinclude:: py2exe_howto/pytables_test.py
   :linenos:

To wrap this script into an executable you have to create a setup script and a
configuration script in your program directory.

The setup script will look like this::

    from setuptools import setup
    import py2exe
    setup(console=['pytables_test.py'])

The configuration script (:file:`setup.cfg`) specifies which modules to be
included and excluded::

    [py2exe]
    excludes= Tkconstants,Tkinter,tcl
    includes= encodings.*, tables.*, numpy.*

As you can see I have included everything from tables (tables.*) and numpy
(numpy.*).

Now you are ready to build the executable file (:file:`pytable_test.exe`).
During the build process a subfolder called *dist* will be created.
This folder contains everything needed for your program.
All dependencies (dll's and such stuff) will be copied into this folder.
When you distribute your application you have to distribute all files and
folders inside the *dist* folder.

Below you can see how to start the build process (`python setup.py py2exe`)::

    c:pytables_test> python3 setup.py py2exe
    ...
    BUILDING EXECUTABLE
    ...

After the build process I enter the *dist* folder and start
:file:`pytables_test.exe`.

::

    c:pytables_test> cd dist

    c:pytables_testdist> pytables_test.exe
    tutorial.h5 (File) 'Test file'
    Last modif.: 'Tue Apr 04 23:09:17 2006'
    Object Tree:
    / (RootGroup) 'Test file'
    /detector (Group) 'Detector information'
    /detector/readout (Table(0,)) 'Readout example'

    [25.0, 36.0, 49.0]

DONE!


-----


.. target-notes::

.. _py2exe: http://www.py2exe.org
=================
PyTables Cookbook
=================

--------
Contents
--------

.. toctree::
    :maxdepth: 1

    hints_for_sql_users
    PyTables & py2exe Howto (by Tommy Edvardsen) <py2exe_howto>
    How to install PyTables when you're not root (by Koen van de Sande) <no_root_install>
    tailoring_atexit_hooks
    custom_data_types
    simple_table
    inmemory_hdf5_files
    threading
:author: FrancescAlted
:date: 2010-04-20 16:44:41


===================================================
SimpleTable: simple wrapper around the Table object
===================================================

Here it is yet another example on how to inherit from the :class:`tables.Table`
object so as to build an easy-to-use Table object.
Thanks to Brent Pedersen for this one (taken from
https://pypi.python.org/pypi/simpletable).

::

    """

    SimpleTable: simple wrapper around pytables hdf5
    ------------------------------------------------------------------------------

    Example Usage::

      >>> from simpletable import SimpleTable
      >>> import tables

      # define the table as a subclass of simple table.
      >>> class ATable(SimpleTable):
      ...     x = tables.Float32Col()
      ...     y = tables.Float32Col()
      ...     name = tables.StringCol(16)

      # instantiate with: args: filename, tablename
      >>> tbl = ATable('test_docs.h5', 'atable1')

      # insert as with pytables:
      >>> row = tbl.row
      >>> for i in range(50):
      ...    row['x'], row['y'] = i, i * 10
      ...    row['name'] = "name_%i" % i
      ...    row.append()
      >>> tbl.flush()

      # there is also insert_many() method() with takes an iterable
      # of dicts with keys matching the colunns (x, y, name) in this
      # case.

      # query the data (query() alias of tables' readWhere()
      >>> tbl.query('(x > 4) & (y < 70)') #doctest: +NORMALIZE_WHITESPACE
      array([('name_5', 5.0, 50.0), ('name_6', 6.0, 60.0)],
            dtype=[('name', '|S16'), ('x', '<f4'), ('y', '<f4')])

    """

    import tables
    _filter = tables.Filters(complib="lzo", complevel=1, shuffle=True)

    class SimpleTable(tables.Table):
        def __init__(self, file_name, table_name, description=None,
                     group_name='default', mode='a', title="", filters=_filter,
                     expectedrows=512000):

            f = tables.openFile(file_name, mode)
            self.uservars = None

            if group_name is None: group_name = 'default'
            parentNode = f._getOrCreatePath('/' + group_name, True)

            if table_name in parentNode: # existing table
                description = None
            elif description is None: # pull the description from the attrs
                description = dict(self._get_description())

            tables.Table.__init__(self, parentNode, table_name,
                           description=description, title=title,
                           filters=filters,
                           expectedrows=expectedrows,
                           _log=False)
            self._c_classId = self.__class__.__name__

        def _get_description(self):
            # pull the description from the attrs
            for attr_name in dir(self):
                if attr_name[0] == '_': continue
                try:
                    attr = getattr(self, attr_name)
                except:
                    continue
                if isinstance(attr, tables.Atom):
                    yield attr_name, attr

        def insert_many(self, data_generator, attr=False):
            row = self.row
            cols = self.colnames
            if not attr:
                for d in data_generator:
                    for c in cols:
                        row[c] = d[c]
                    row.append()
            else:
                for d in data_generator:
                    for c in cols:
                        row[c] = getattr(d, c)
                    row.append()
            self.flush()

        query = tables.Table.readWhere

    # convience sublcass that i use a lot.
    class BlastTable(SimpleTable):
          query      = tables.StringCol(5)
          subject    = tables.StringCol(5)

          pctid      = tables.Float32Col()
          hitlen     = tables.UInt16Col()
          nmismatch  = tables.UInt16Col()
          ngaps      = tables.UInt16Col()

          qstart     = tables.UInt32Col()
          qstop      = tables.UInt32Col()
          sstart     = tables.UInt32Col()
          sstop      = tables.UInt32Col()

          evalue     = tables.Float64Col()
          score      = tables.Float32Col()


    if __name__ == '__main__':
        import doctest
        doctest.testmod()
        import os
        os.unlink('test_docs.h5')

=========
Threading
=========

.. py:currentmodule:: tables


Background
==========

Several bug reports have been filed in the past by the users regarding
problems related to the impossibility to use PyTables in multi-thread
programs.

The problem was mainly related to an internal registry that forced the
sharing of HDF5 file handles across multiple threads.

In PyTables 3.1.0 the code for file handles management has been completely
redesigned (see the *Backward incompatible changes* section in 
:doc:`../release-notes/RELEASE_NOTES_v3.1.x`) to be more simple and
transparent and to allow the use of PyTables in multi-thread programs.

Citing the :doc:`../release-notes/RELEASE_NOTES_v3.1.x`::

    It is important to stress that the new implementation still has an
    internal registry (implementation detail) and it is still
    **not thread safe**.
    Just now a smart enough developer should be able to use PyTables in a
    muti-thread program without too much headaches.


A common schema for concurrency
===============================

Although it is probably not the most efficient or elegant solution to solve
a certain class of problems, many users seems to like the possibility to
load a portion of data and process it inside a *thread function* using
multiple threads to process the entire dataset.

Each thread is responsible of:

* opening the (same) HDF5 file for reading,
* load data from it and
* close the HDF5 file itself

Each file handle is of exclusive use of the thread that opened it and
file handles are never shared across threads.

In order to do it in a safe way with PyTables some care should be used
during the phase of opening and closing HDF5 files in order ensure the
correct behaviour of the internal machinery used to manage HDF5 file handles.


Very simple solution
====================

A very simple solution for this kind of scenario is to use a
:class:`threading.Lock` around part of the code that are considered critical
e.g. the :func:`open_file` function and the :meth:`File.close` method::

    import threading
    
    lock = threading.Lock()

    def synchronized_open_file(*args, **kwargs):
        with lock:
            return tb.open_file(*args, **kwargs)

    def synchronized_close_file(self, *args, **kwargs):
        with lock:
            return self.close(*args, **kwargs)


The :func:`synchronized_open_file` and :func:`synchronized_close_file` can
be used in the *thread function* to open and close the HDF5 file::

    import numpy as np
    import tables as tb

    def run(filename, path, inqueue, outqueue):
        try:
            yslice = inqueue.get()
            h5file = synchronized_open_file(filename, mode='r')
            h5array = h5file.get_node(path)
            data = h5array[yslice, ...]
            psum = np.sum(data)
        except Exception as e:
            outqueue.put(e)
        else:
            outqueue.put(psum)
        finally:
            synchronized_close_file(h5file)


Finally the main function of the program:

* instantiates the input and output :class:`queue.Queue`,
* starts all threads, 
* sends the processing requests on the input :class:`queue.Queue`
* collects results reading from the output :class:`queue.Queue`
* performs finalization actions (:meth:`threading.Thread.join`)

.. code-block:: python

    import os
    import queue
    import threading

    import numpy as np
    import tables as tb

    SIZE = 100
    NTHREADS = 5
    FILENAME = 'simple_threading.h5'
    H5PATH = '/array'

    def create_test_file(filename):
        data = np.random.rand(SIZE, SIZE)

        with tb.open_file(filename, 'w') as h5file:
            h5file.create_array('/', 'array', title="Test Array", obj=data)

    def chunk_generator(data_size, nchunks):
        chunk_size = int(np.ceil(data_size / nchunks))
        for start in range(0, data_size, chunk_size):
            yield slice(start, start + chunk_size)

    def main():
        # generate the test data
        if not os.path.exists(FILENAME):
            create_test_file(FILENAME)

        threads = []
        inqueue = queue.Queue()
        outqueue = queue.Queue()

        # start all threads
        for i in range(NTHREADS):
            thread = threading.Thread(
                target=run, args=(FILENAME, H5PATH, inqueue, outqueue))
            thread.start()
            threads.append(thread)

        # push requests in the input queue
        for yslice in chunk_generator(SIZE, len(threads)):
            inqueue.put(yslice)

        # collect results
        try:
            mean_ = 0.

            for i in range(len(threads)):
                out = outqueue.get()
                if isinstance(out, Exception):
                    raise out
                else:
                    mean_ += out

            mean_ /= SIZE * SIZE

        finally:
            for thread in threads:
                thread.join()

        # print results
        print('Mean: {}'.format(mean_))

    if __name__ == '__main__':
        main()

The program in the example computes the mean value of a potentially huge
dataset splinting the computation across :data:`NTHREADS` (5 in this case)
threads.

The complete and working code of this example (Python 3 is required) can be
found in the :file:`examples` directory:
:download:`simple_threading.py <../../../examples/simple_threading.py>`.

The approach presented in this section is very simple and readable but has
the **drawback** that the user code have to be modified to replace
:func:`open_file` and :meth:`File.close` calls with their safe version
(:func:`synchronized_open_file` and :func:`synchronized_close_file`).

Also, the solution shown in the example does not cover the entire PyTables
API (e.g. although not recommended HDF5 files can be opened using the
:class:`File` constructor) and makes it impossible to use *pythonic*
constructs like the *with* statement::

    with tb.open_file(filename) as h5file:
        do_something(h5file)


Monkey-patching PyTables
========================

An alternative implementation with respect to the `Very simple solution`_
presented in the previous section consists in monkey-patching the PyTables
package to replace some of its components with a more thread-safe version of
themselves::

    import threading

    import tables as tb
    import tables.file as _tables_file

    class ThreadsafeFileRegistry(_tables_file._FileRegistry):
        lock = threading.RLock()

        @property
        def handlers(self):
            return self._handlers.copy()

        def add(self, handler):
            with self.lock:
                return super().add(handler)

        def remove(self, handler):
            with self.lock:
                return super().remove(handler)

        def close_all(self):
            with self.lock:
                return super().close_all(handler)

    class ThreadsafeFile(_tables_file.File):
        def __init__(self, *args, **kargs):
            with ThreadsafeFileRegistry.lock:
                super().__init__(*args, **kargs)

        def close(self):
            with ThreadsafeFileRegistry.lock:
                super().close()

    @functools.wraps(tb.open_file)
    def synchronized_open_file(*args, **kwargs):
        with ThreadsafeFileRegistry.lock:
            return _tables_file._original_open_file(*args, **kwargs)

    # monkey patch the tables package
    _tables_file._original_open_file = _tables_file.open_file
    _tables_file.open_file = synchronized_open_file
    tb.open_file = synchronized_open_file

    _tables_file._original_File = _tables_file.File
    _tables_file.File = ThreadsafeFile
    tb.File = ThreadsafeFile

    _tables_file._open_files = ThreadsafeFileRegistry()


At this point PyTables can be used transparently in the example program presented
in the previous section.
In particular the standard PyTables API (including *with* statements) can be
used in the *thread function*::

    def run(filename, path, inqueue, outqueue):
        try:
            yslice = inqueue.get()
            with tb.open_file(filename, mode='r') as h5file:
                h5array = h5file.get_node(path)
                data = h5array[yslice, ...]
            psum = np.sum(data)
        except Exception as e:
            outqueue.put(e)
        else:
            outqueue.put(psum)


The complete code of this version of the example can be found in the
:file:`examples` folder:
:download:`simple_threading.py <../../../examples/threading_monkeypatch.py>`.
Python 3 is required.


