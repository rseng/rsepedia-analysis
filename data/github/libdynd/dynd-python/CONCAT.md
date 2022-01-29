Running The Python Tests
========================

The Python tests are written using Python's built in unittest module,
and can either be run by executing the scripts or with nose.  The two
modules dynd.nd and dynd.ndt have separate tests. To run a test script
directly, simply execute the test .py file.

Examples:

    D:\Develop\dynd-python>python dynd\ndt\test\test_type_basics.py
    D:\Develop\dynd-python>python dynd\nd\test\test_array_basics.py

To execute the full test suite for both modules, run the following
(not from the root project directory):

    D:\Develop>python -c "import dynd; dynd.test()"
    Running unit tests for the DyND Python bindings
    Python version: 3.3.2 |Continuum Analytics, Inc.| (default, May 17 2013, 11:32:27) [MSC v.1500 64 bit (AMD64)]
    Python prefix: C:\Anaconda\envs\py33
    DyND-Python module: C:\Anaconda\envs\py33\lib\site-packages\dynd
    DyND-Python version: 0.3.1-8-gd0620eb-dirty
    DyND-Python git sha1: d0620ebad73a9ef840e40c9649ba97bff9444e0c
    LibDyND version: 0.3.1-2-g25d6d4e
    LibDyND git sha1: 25d6d4ef6ba4b1213359b0da80897c04c8d7474a
    NumPy version: 1.7.1
    ..........................................................................................................
    ----------------------------------------------------------------------
    Ran 106 tests in 0.348s

    OK

To generate Jenkins-compatible XML output, control the verbosity of the
output, etc, you may pass some extra parameters to the test function.

    D:\Develop>python
    Python 3.3.2 |Continuum Analytics, Inc.| (default, May 17 2013, 11:32:27) [MSC v.1500 64 bit (AMD64)] on win32
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import dynd
    >>> help(dynd.test)
    Help on function test in module dynd:

    test(verbosity=1, xunitfile=None, exit=False)
        Runs the full DyND test suite, outputing
        the results of the tests to  sys.stdout.

        Parameters
        ----------
        verbosity : int, optional
            Value 0 prints very little, 1 prints a little bit,
            and 2 prints the test names while testing.
        xunitfile : string, optional
            If provided, writes the test results to an xunit
            style xml file. This is useful for running the tests
            in a CI server such as Jenkins.
        exit : bool, optional
            If True, the function will call sys.exit with an
            error code after the tests are finished.
DyND-Python
===========

TravisCI: [![Build Status](https://api.travis-ci.org/libdynd/dynd-python.svg?branch=master)](https://travis-ci.org/libdynd/dynd-python)
AppVeyor: [![Build Status](https://ci.appveyor.com/api/projects/status/cv2bnq3oghe4nqnj/branch/master?svg=true)](https://ci.appveyor.com/project/libdynd/dynd-python/branch/master)

DyND-Python, a component of [the Blaze project](http://blaze.pydata.org/),
is the Python exposure of [the DyND dynamic multi-dimensional array library](http://libdynd.org).

To discuss the development of this library, subscribe to the
[LibDyND Development List](https://groups.google.com/forum/#!forum/libdynd-dev).

Python versions 2.7, 3.4 and 3.5 are supported.

http://libdynd.org

https://github.com/libdynd/libdynd

https://github.com/libdynd/dynd-python

Trying Out DyND
---------------

The easiest way to try it out is through the Anaconda
Python distribution. The latest release of Anaconda includes
a version of DyND.

http://continuum.io/downloads

For trying the latest updates, there is also an automated
build configured which tracks the latest git master. When
all the tests pass, it uploads conda packages to the anaconda.org
channel "dynd/channel/dev". To get these versions, you
can run the following command.

```
conda install dynd-python --channel dynd/channel/dev
```

It may work best to install development versions of DyND into an
environment instead of the main Anaconda directory.
You can do this with a command like:

```
C:\>conda create -n dynd-env python=3.3 dynd-python --channel dynd/channel/dev
```

Developing DyND
---------------

See the [build and install instructions](BUILD_INSTALL.md) for details on
building the software for environments not supported by Anaconda, or
if you would like to modify or contribute to the project.
STEP BY STEP BUILD AND INSTALL
==============================

1. Check the C++ compiler version.

  Ensure you have a suitable C++14 compiler. On Windows, Visual
Studio 2015 is the minimum supported compiler. On Mac OS X, clang is the
recommended compiler. On Linux, gcc 4.9 and later, and
clang 3.4 and later have been tested.

2. Get the prerequisites.
  * CMake >= 2.8.11
  * Python 2.7, 3.4, or 3.5
  * Cython >= 0.24
  * NumPy >= 1.7.1
  * git (for cloning the github repositories)
  * Nose (Only for generating xunit .xml output when running tests)

BUILD FOR DEVELOPMENT
---------------------

Development of DyND should be done using the build configuration combining
libdynd and the DyND Python bindings in a single combined build in a
development mode. With this configuration, a single build will update both
libdynd and dynd-python in a way that both the C++ and Python tests can
be run with no further installation steps.

This combined build configuration works both with make-style builds and
MSVC solution files. The following instructions are the same for both
Windows and Linux/OS X, but on a unix platform run `make` instead of
loading the `dynd-python.sln` file at the end.

  ```
  C:\>git clone --recursive https://github.com/libdynd/dynd-python
  Cloning into 'dynd-python'...
  <...>
  C:\>cd dynd-python
  C:\dynd-python>git clone --recursive https://github.com/libdynd/libdynd
  Cloning into 'libdynd'...
  <...>
  C:\dynd-python>python setup.py develop
  <...>
  C:\dynd-python>cd build-dev
  C:\dynd-python\build-dev>start dynd-python.sln
  ```

BUILD FOR INSTALLATION
----------------

3. Get the source code.

  Check out the dynd-python and libdynd source code. The following commands
should work equivalently on Windows and Unix-like operating systems.

  ```
  ~ $ git clone --recursive https://github.com/libdynd/libdynd
  Cloning into libdynd...
  ~ $ git clone --recursive https://github.com/libdynd/dynd-python
  Cloning into dynd-python...
  ```

4. Build and install libdynd

  ```
  ~ $ cd libdynd
  ~/libdynd $ mkdir build
  ~/libdynd $ cd build
  ~/libdynd/build $ cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
  <...>
  ~/libdynd/build $ make
  <...>
  ~/libdynd/build $ sudo make install
  ```

  If you want to control where the libdynd shared object is
installed, add `-DCMAKE_INSTALL_PREFIX=<prefix>`
to the cmake command.

4. Build and install dynd-python

  ```
  ~ $ cd dynd-python
  ~/dynd-python $ python setup.py install
  ```

  If you want to control where the Python module goes, add
`-DPYTHON_PACKAGE_INSTALL_PREFIX=<site-pkg-dir>`
to the cmake command.


ALTERNATIVE COMPILERS
=====================

  If you want to build with a different compiler, for
  example to use the static analyzer in clang, you can
  customize the compiler at this step.

  ```
  ~/libdynd $ mkdir build-analyze
  ~/libdynd $ cd build-analyze
  ~/libdynd/build-analyze $ export CCC_CC=clang
  ~/libdynd/build-analyze $ export CCC_CXX=clang++
  ~/libdynd/build-analyze $ cmake -DCMAKE_CXX_COMPILER=c++-analyzer -DCMAKE_C_COMPILER=ccc-analyzer ..
  ```

CONFIGURATION OPTIONS
=====================

These are some options which can be configured by calling
CMake with an argument like "-DCMAKE_BUILD_TYPE=Release".

CMAKE_BUILD_TYPE
    Which kind of build, such as Release, RelWithDebInfo, Debug.
CMAKE_INSTALL_PREFIX
    The prefix for installing shared libraries such as
    libdynd.so.

# PyObject to DyND Array and Back

For seamless interoperability with Python, DyND needs extremely solid conversion between Python objects and nd::arrays. What exactly that means isn't completely obvious, as NumPy's experience shows, in having some excellent behavior but sometimes doing unexpected things.

## Requirements

* Assignment needs to be pluggable
  * Dynamically registered DyND types should be able to specify Python type conversions
  * Newly defined Python types should be able to specify conversions
* Automatic type deductions should be unsurprising, or fail if no reasonable choice (e.g. avoid copying NumPy 'object' array behaviour)
* Trivial cases with Python native types should be very fast (e.g. `nd.array(3.14)`, `nd.array([1,2,3,4])`)
* Nested generators and iterators should work intuitively. This means conversion with deduction should make only one pass through the PyObject dimensional structure, because processing an iterator consumes it.

## The Types Framing the Discussion

All DyND conversions are defined in terms of a source type and destination type, e.g. `(float32) -> int64`. In DyND, the `nd::array` object can be placed inside another `nd::array` via the `array[Any]` type. This allows signatures for general conversions to/from DyND to be expressed the same way as specific conversions for types like `float32`. The Python bindings define a type `pyobject` which holds a single reference to a Python object, i.e. a `PyObject *` pointer in the case of the CPython interpreter. 

Combining these give us the two general conversion signatures, `(pyobject) -> array[Any]` and `(array[Any]) -> pyobject`.

## Example Conversions

Instead of beginning this with general rules for conversion, let's go through a bunch of specific examples and use our intuition to define what we think the conversion should be. Then we can take the collection of examples and produce a system which matches it. Iterating the system through finding examples that behave unintuitively, adding them to this list, and tweaking the system to account for them should hopefully converge to a good end result.

### `(pyobject) -> array[Any]`

#### Scalar cases

The following is from existing tests as well as some additional examples

```python
# Python object => DyND type
True                                           => "bool"
10                                             => "int32"
-2200000000                                    => "int64"
5.125                                          => "float64"
5.125 - 2.5j                                   => "complex[float64]"
'abcdef'                                       => "string"
u'abcdef'                                      => "string"
b'abcdef'                                      => "bytes"          # On Python 3
```

#### Array cases

```python
[]                                             => "0 * int32"      # Current choice for empty dynamic list
[[], [], []]                                   => "3 * 0 * int32"
[1, 2, 3]                                      => "3 * int32"
[True, False]                                  => "2 * bool"
[1, True]                                      => "2 * int32"      # Current integer behaviour
[10000000000, 1, False]                        => "3 * int64"
[10000000000, 3.25, 2, False]                  => "4 * float64"
[3.25j, 3.25, 1, 2, True]                      => "5 * complex[float64]"
[str(x) + 'test' for x in range(10)]           => "10 * string"
[u'test', 'test2']                             => "2 * string"
[b'x'*x for x in range(10)]                    => "10 * bytes"    # On Python 3
[[True, 2, 3], [4, 5, 6.5], [1, 2, 3]]         => "3 * 3 * float64"
[[1], [2, 3, 4], [5, 6]]                       => "3 * var * int32"
[[True, False], [False, 2, 3], [-10000000000], [True, 10, 3.125, 5.5j]] => "4 * var * complex[float64]"
[[], [False, 2, 3]]                            => "2 * var * int32"
[[], [[]], [[[1, 3]]]]                         => "3 * var * var * 2 * int32"
```

#### Iterator protocol cases

Generators and other iterators can present tricky cases.

```python
(x for x in [])                                => "0 * int32" # Equivalent to []
iter([...]) for all list cases above           => "..."       # Same as the list version but as iterator
```

#### Error cases

The following should throw errors (or we should at least think hard about how "smart" DyND should be).

```python
[[1], [[2]]]       # The 1 and 2 values imply different numbers of dimensions
[1, "test"]        # Implicitly converting numbers to strings not allowed
[b"test", "test"]  # (Python 3) Implicitly converting bytes to strings or vice versa is not allowed

```

#### Questions

* What should the DyND type be for an empty dynamic list? This question applies across all such contexts, and should be the same for JSON parsing, Python object conversion, etc. Current choice is `0 * int32`.
* What should the deduced integer type be?
  * `int32` if it fits, `int64` if not, `int128` if still not, `bigint` if still not.
  * `int64` if it fits, `int128` if not, `bigint` if still not
  * `bigint` (SSO-based arbitrary-sized integer)
* Do we like Python's implicit bool to int conversion? Is it an option for us to raise an error for it instead?

### `(array[Any]) -> pyobject`

# Developer Documentation

This folder contains design and development documentation intended
for developers of the DyND Python bindings.
﻿Vagrant Images
==============

This directory contains a few vagrant images which
can be used to build and develop on the latest DyND.
To get started, you'll first want to install
[Vagrant](http://www.vagrantup.com/downloads.html)
and [VirtualBox](https://www.virtualbox.org/wiki/Downloads).

These images work on Linux, Mac OS X, and Windows, and
are an easy way to get a consistent, isolated environment
for Blaze development and experimentation.

Starting a Vagrant Image
------------------------

To provision and start one of these images, use
the following commands:

```
$ cd trusty64-py34
$ make start
```

To connect with ssh, use:

```
$ make login
```

On Windows,
[PuTTY](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html)
is a good way to ssh into these images. If you run
`vagrant ssh-config`, you will see output similar to
the following:

```
Host default
  HostName 127.0.0.1
  User vagrant
  Port 2222
  UserKnownHostsFile /dev/null
  StrictHostKeyChecking no
  PasswordAuthentication no
  IdentityFile C:/Users/[name]/.vagrant.d/insecure_private_key
  IdentitiesOnly yes
  LogLevel FATAL
```

You can use the PuTTYgen tool to convert the IdentityFile
private key from ssh format to the PuTTY .ppk format.
Then create a PuTTY configuration using the settings
and converted .ppk file.

Updating VirtualBox Guest Additions
-----------------------------------

You may find it useful to install the vagrant-vbguest
plugin, which automatically keeps the VirtualBox guest
additions up to date. When these are out of date,
features such as access to the `/vagrant` folder
may not work reliably. The command is:

```
$ vagrant plugin install vagrant-vbguest
```
