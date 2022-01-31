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
DyND GFuncs
===========

NOTE: GFuncs are not currently operational due to refactoring
of the core object representation.

GFuncs are the array-programming primitive functions created and
used by blaze-local. They are currently in a very preliminary form,
but as they are they provide a preview of what is to come.

Lazy Evaluation
---------------

GFuncs and many other operations, like simple arithmetic, are
evaluated in a lazy function. This means that an expression
tree is built up, and the system has an opportunity to analyze
the complete expression DAG before evaluating the final values
when they are requested.

A simple illustration of this is to create an expression,
then evaluate it before and after modifying one of its inputs.

.. code-block:: python

    >>> a = nd.array([1,2])
    >>> b = a + 3
    >>> b
    nd.array([4, 5], int32)

    >>> a[0].val_assign(10)
    >>> b
    nd.array([13, 5], int32)

The evaluation engine is currently very primitive, and does not
support evaluating expressions with complicated expression graphs.
To evaluate such expressions, use the ``vals()`` method to create
intermediate results.

.. code-block:: python

    >>> (a + 3) * 2
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "_pydnd.pyx", line 289, in _pydnd.w_ndarray.__repr__ (D:\Develop\blaze\build\_pydnd.cxx:4170)
    RuntimeError: evaluating this expression graph is not yet supported:
    ("elementwise_binary_kernel",
      [SNIP]
    )

    >>> (a + 3).vals() * 2
    nd.array([8, 10], int32)

Using the Builtin GFuncs
------------------------

Blaze-local has a few gfuncs built in, at this stage primarily to
demonstrate the system. There is no implicit type promotion, so
any calls to the functions much have dtypes precisely matching
a function signature within the gfunc.

.. code-block:: python

    >>> nd.sin(3.1)
    nd.array(0.0415807, float64)

    >>> nd.sin(3)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "_pydnd.pyx", line 374, in _pydnd.w_elwise_gfunc.__call__ (D:\Develop\blaze\build\_pydnd.cxx:5735)
    RuntimeError: sin: could not find a gfunc kernel matching input argument types (int32)

    >>> nd.sum([1,3,5,9])
    nd.array(18, int32)

Creating Elementwise GFuncs
---------------------------

GFuncs are created using the constructors in the gfunc namespace. Once created, kernel
functions can be added, which presently must be provided as ctypes function pointer.
The ``nd.elwise_kernels`` namespace has a number of elementwise kernels to play with.

To demonstrate gfunc creation, we make a gfunc with a few different kernels that behave
differently. Normally one would want there to be a system to how overloads are done,
just as in C++, Java, or C#, this is only for demonstration purposes.

.. code-block:: python

    >>> myfunc = nd.gfunc.elwise('myfunc')
    >>> myfunc.add_kernel(nd.elwise_kernels.add_int32)
    >>> myfunc.add_kernel(nd.elwise_kernels.multiply_float64)
    >>> myfunc.add_kernel(nd.elwise_kernels.square_float64)

Now we can call the function if we provide operands of the right type.

.. code-block:: python

    >>> myfunc(3.0)
    nd.array(9, float64)

    >>> myfunc(1,2)
    nd.array(3, int32)

    >>> myfunc(1.0, 2.0)
    nd.array(2, float64)

    >>> myfunc(1, 2.0)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "_pydnd.pyx", line 374, in _pydnd.w_elwise_gfunc.__call__ (D:\Develop\blaze\build\_pydnd.cxx:5735)
    RuntimeError: myfunc: could not find a gfunc kernel matching input argument types (int32, float64)

Creating Elementwise Reduction GFuncs
-------------------------------------

Blaze-local supports simple reductions, like ``sum`` and ``min`` as well, through
elementwise reduction gfuncs. Presently, these require binary kernels whose output
is the same type as the two inputs. This will expand to other types in the future,
for example if an identity is provided, the output type could be different from the
input array.

We can make an example reduction that does ``sum`` or ``product`` depending on
the input dtype.

.. code-block:: python

    >>> myred = nd.gfunc.elwise_reduce('myred')
    >>> myred.add_kernel(nd.elwise_kernels.add_int32, associative=True, commutative=True, identity=0)
    >>> myred.add_kernel(nd.elwise_kernels.multiply_float64, associative=True, commutative=True, identity=1)

    >>> myred([1,2,3,4])
    nd.array(10, int32)
    >>> myred([1.,2.,3.,4.])
    nd.array(24, float64)

Groupby Reductions
------------------

Blaze-local has a simple ``nd.groupby`` function which, when combined with elementwise
reductions, can be used for groupby reductions. Here's a simple example.

.. code-block:: python

    >>> data = np.array([0, 1, 2, 3, 4, 5, 6, 7])
    >>> by = np.array(['a', 'a', 'c', 'a', 'b', 'c', 'a', 'd'])
    >>> groups = nd.factor_categorical_dtype(by)
    >>> gb = nd.groupby(data, by, groups)

    >>> print(groups)
    categorical<fixedstring<ascii,1>, ["a", "b", "c", "d"]>

    >>> print("max:     ", nd.max(gb, axis=1))
    ('max:     ', nd.array([6, 4, 5, 7], int32))

    >>> print("min:     ", nd.min(gb, axis=1))
    ('min:     ', nd.array([0, 4, 2, 7], int32))

    >>> print("sum:     ", nd.sum(gb, axis=1))
    ('sum:     ', nd.array([10, 4, 7, 7], int32))

    >>> print("product: ", nd.product(gb, axis=1))
    ('product: ', nd.array([0, 4, 10, 7], int32))

DyND Types
==========

DyND has a preliminary set of types, similar to the ones
in NumPy. The string representation of types is based on
the `DataShape grammar <https://github.com/ContinuumIO/datashape/blob/master/docs/source/grammar.rst>`_.

This means types are printed using angle brackets [] for the
parameters determining the nature of the type.

Primitive Types
---------------

To start are the primitive numeric types with a fixed size.

=================== =============== =======================================
Type                 Size (bytes)    Notes
=================== =============== =======================================
ndt.bool             1               True (value 0) or False (value 1)
ndt.int8             1               8 bit 2's complement signed integer
ndt.int16            2               16 bit 2's complement signed integer
ndt.int32            4               32 bit 2's complement signed integer
ndt.int64            8               64 bit 2's complement signed integer
ndt.uint8            1               8 bit unsigned integer
ndt.uint16           2               16 bit unsigned integer
ndt.uint32           4               32 bit unsigned integer
ndt.uint64           8               64 bit unsigned integer
ndt.float32          4               32-bit IEEE floating point
ndt.float64          8               64-bit IEEE floating point
ndt.complex_float32  8               complex made of 32-bit IEEE floats
ndt.complex_float64  16              complex made of 64-bit IEEE floats
=================== =============== =======================================

Bytes Type
----------

There is a fixed size bytes type, which has a configurable size
and alignment. For example, one could make a bytes type with
size 8 and alignment 4 to match the storage properties of complex
float32.

.. code-block:: python

    >>> ndt.make_fixedbytes(8, 4)
    ndt.type('fixedbytes[8,4]')

Unaligned Type
--------------

Any memory containing data for a type must always obey the type's
alignment. This is different from Numpy, where data may have different
alignment, and NumPy checks whether there is alignment and
inserts buffering where necessary. To deal with unaligned data,
a type which adapts it into aligned storage must be used.

.. code-block:: python

    >>> ndt.make_unaligned(ndt.float64)
    ndt.type('unaligned[float64]')

Byteswap Type
-------------

When data has the wrong endianness for the platform, it must be
byteswapped before it can be used. The `byteswap` adapter type
plays this role.

.. code-block:: python

    >>> ndt.make_byteswap(ndt.int32)
    ndt.type('byteswap[int32]')

Convert Type
------------

Conversion between types is handled by the `convert` type. The
error mode used for the conversion can be controlled with an extra
parameter.

.. code-block:: python

    >>> ndt.make_convert(ndt.int32, ndt.float64)
    ndt.type('convert[to=int32, from=float64]')

    >>> ndt.make_convert(ndt.int32, ndt.float64, errmode='overflow')
    ndt.type('convert[to=int32, from=float64, errmode=overflow]')


String Types
------------

There are presently two different kinds of string types. There are
fixed size strings, similar to those in NumPy, where the string data
is null-terminated within a buffer of particular size.

There are also variable size strings, whose data is in a separate
memory block, reference counted at the memory block level (blockref).
These strings cannot be resized in place, but are designed to be used
in a functional style where results are created once, not repeatedly
mutated in place.

Also planned are strings with full dynamic behavior, where each string
manages its own memory allocation. These are not implemented yet.

To create string types, use the `ndt.make_string` and
`ndt.make_fixedstring`. There is a default string type
`ndt.string`, which is UTF-8.

.. code-block:: python

    >>> ndt.string
    ndt.string

    >>> ndt.make_string('ascii')
    ndt.type('string['ascii']')

    >>> ndt.make_fixedstring(16, 'utf_32')
    ndt.type('string[16,'utf32']')

When creating nd::array objects from Python lists, blockref strings
are used by default.

.. code-block:: python

    >>> nd.array(['abcdefg', u'안녕', u'Testing'])
    nd.array(["abcdefg", "\uc548\ub155", "Testing"], type="strided * string")

Categorical Type
----------------

There is a preliminary categorical type, used by the `nd.groupby`
function.

.. code-block:: python

    >>> groups = nd.array(['a', 'b', 'c'],
                     dtype=ndt.make_fixedstring(1, 'ascii'))
    >>> ndt.make_categorical(groups)
    ndt.type('categorical[string[1,'ascii'], ["a", "b", "c"]]')

Pointer Type
------------

This type presently exists to help with `ctypes` function pointer
interoperability, but eventually will behave in a blockref manner,
similar to the blockref string type.

.. code-block:: python

    >>> ndt.make_pointer(ndt.complex_float32)
    ndt.type('pointer[complex[float32]]')

Array Types
-----------

There are a few different array types, with different properties
with respect to the flexibility of size and memory layout.

.. code-block:: python

    >>> ndt.strided * ndt.float32
    ndt.type('strided * float32')

    >>> ndt.fixed[10] * ndt.bool
    ndt.type('10 * bool')

    >>> ndt.cfixed[12] * ndt.int32
    ndt.type('cfixed[12] * int32')
    
    >>> ndt.strided * ndt.var * '{x : int32, y : float32}'
    ndt.type('strided * var * {x : int32, y : float32}')

========================================
Converting Python Objects To DyND Arrays
========================================

In order to provide a user friendly boundary between Python objects and
DyND arrays, the Python exposure of LibDyND does some fairly extensive
analysis of input objects passed to the ``nd.array`` constructor. This
document provides an overview of the different ways this occurs.

There are a few ways provided to construct arrays from Python
objects, generally through the ``nd.array`` constructor. There
is also the ``nd.asarray`` function which is not as fully fleshed
out, which retains a view where possible such as when a buffer
interface is present.::

    >>> # Array from an object - fully autodetected
    >>> nd.array([1, 2, 3])
    nd.array([1, 2, 3], type="strided * int32")

    >>> # Array from an object - dtype provided
    >>> nd.array([1, 2, 3], dtype=ndt.int16)
    nd.array([1, 2, 3], type="strided * int16")

    >>> # Array from an object - full type provided
    >>> nd.array([1, 2, 3], type='strided * int16')
    nd.array([1, 2, 3], type="strided * int16")

The entry point for this is in the array class's ``__cinit__`` function
in Cython.

https://github.com/libdynd/dynd-python/blob/master/src/_pydynd.pyx#L971

This calls the ``array_init_from_pyobject`` function, which has two overloads
defined depending on whether a type is provided or not.

https://github.com/libdynd/dynd-python/blob/master/src/array_functions.cpp#L234

These call ``array_from_py``, similarly with and without type information::

https://github.com/libdynd/dynd-python/blob/master/src/array_from_py.cpp#L368

https://github.com/libdynd/dynd-python/blob/master/src/array_from_py.cpp#L592

Constructing With A Known Full Type
===================================

The simplest case is when the type of the desired array is fully known,
for example::

    >>> nd.array([1, 2, 3], type="strided * int16")
    nd.array([1, 2, 3], type="strided * int16")

To be able to do this conversion, DyND needs to do these steps:

* Traverse the object to determine its shape (3,)
* Create an empty ``nd.array`` using the shape information and type
* Copy the Python object into the array

This code is here:

https://github.com/libdynd/dynd-python/blob/master/src/array_from_py.cpp#L685

where it first checks whether the shape is required, and deduces
the shape, then calls ``array_nodim_broadcast_assign_from_py`` which
does an assignment from the pyobject without allowing new dimensions
to be broadcast. This function is defined in this file:

https://github.com/libdynd/dynd-python/blob/master/src/array_assign_from_py.cpp

Constructing with a Known DType
===============================

A little bit more complicated is the case to construct the array
when the dtype is known, but dimensions should be deduced from the
PyObject::

    >>> nd.array([1, 2, 3], dtype=ndt.int16)
    nd.array([1, 2, 3], type="strided * int16")

    >>> nd.array([['test', 1], ['two', 3]],
                 dtype='{name: string, value: int32}')
    nd.array([["test", 1], ["two", 3]], type="strided * {name : string, value : int32}")

To be able to do this, DyND needs to:

* Deduce the number of dimensions and shape of the object,
  using the dtype as a hint. This can be tricky, especially
  allowing for structure initialization from sequences.
* Create an empty ``nd.array`` using the shape information and type
* Copy the Python object into the array

This code is here:

https://github.com/libdynd/dynd-python/blob/master/src/array_from_py.cpp#L595

The code which deduces the shape and number of dimensions, then reconciles
it with the provided dtype is here:

https://github.com/libdynd/dynd-python/blob/master/src/array_from_py.cpp#L626

And finally the object is copied with the same function as in the
full type case, ``array_nodim_broadcast_assign_from_py``.

Constructing Fully Automatically
================================

When the PyObject isn't an instantly recognizable type, it gets to here:

https://github.com/libdynd/dynd-python/blob/master/src/array_from_py.cpp#L536

where it calls ``array_from_py_dynamic``. This function dynamically
updates the type, promoting both the dtype and dimension types
as it goes. This function is here:

https://github.com/libdynd/dynd-python/blob/master/src/array_from_py_dynamic.cpp

This function does need some more work, it does not support numpy
scalar types and arrays intermixed with iterators, for example.

DyND Array
==========

In this preview release of DyND, there is an `nd.array` object
which hold a multi-dimensional array of elements of a specific
`ndt.type`.

One of the important features of this library is interoperability with
NumPy that is as seamless as possible. When a dtype is compatible across
both systems, arrays from NumPy and DyND will implicitly move
between the libraries.

Importing DyND
--------------

The standard way to import DyND is with the command::

    >>> from dynd import nd, ndt

Constructing From Python Scalars
--------------------------------

Currently, scalars are always zero-dimensional arrays in the system.
This is likely to change in a future preview release.

There's a difference with NumPy in the default treatment of integers.
In NumPy, the C/C++ ``long`` type is used, which is 32-bits on 32-bit
platforms and on 64-bit Windows, and which is 64-bits on Linux and OS X.
In DyND if the Python integer fits in 32-bits, it will always
be initialized by default as a 32-bit integer.

.. code-block:: python

    >>> nd.array(0)
    nd.array(0, type="int32")

    >>> nd.array(3.14)
    nd.array(3.14, type="float64")

    >>> nd.array(1 + 1j)
    nd.array((1,1), type="complex[float64]")

Strings default to `blockref` strings, which are variable-sized strings.
The support for them is still preliminary, but some basic functionality
like converting between different unicode encodings is implemented.

.. code-block:: python

    >>> nd.array('testing')
    nd.array("testing", type="string")

Constructing from Python Lists
------------------------------

Similar to in NumPy, arrays can be constructed from lists of
objects. This code does not try to be as clever as NumPy, and
will fail if something is inconsistent in the input data.

.. code-block:: python

    >>> nd.array([True, False])
    nd.array([true, false], type="strided * bool")

    >>> nd.array([1,2,3])
    nd.array([1, 2, 3], type="strided * int32")

    >>> nd.array([[1.0,0],[0,1.0]])
    nd.array([[1, 0], [0, 1]], type="strided * strided * float64")

    >>> nd.array(["testing", "one", u"two", "three"])
    nd.array(["testing", "one", "two", "three"], type="strided * string")

Converting to Python Types
--------------------------

To convert back into native Python objects, there is an ``nd.as_py(a)``
function.

.. code-block:: python

    >> x = nd.array([True, False])
    >> nd.as_py(x)
    [True, False]

    >> x = nd.array("testing")
    >> nd.as_py(x)
    u'testing'

Constructing from NumPy Scalars
-------------------------------

Numpy scalars are also supported as input, and the dtype is preserved
in the conversion.

.. code-block:: python

    >>> import numpy as np

    >>> x = np.bool_(False)
    >>> nd.array(x)
    nd.array(false, type="bool")

    >>> x = np.int16(1000)
    >>> nd.array(x)
    nd.array(1000, type="int16")

    >>> x = np.complex128(3.1)
    >>> nd.array(x)
    nd.array((3.1,0), type="complex[float64]")

Constructing from NumPy Arrays
------------------------------

When the dtype is supported by DyND, NumPy arrays can
be converted into DyND arrays. When using the `nd.array` constructor,
a new array is created, but there is also `nd.asarray` which creates
a view when possible, and `nd.view` which always creates a view and
raises if that is not possible.

.. code-block:: python

    >>> x = np.arange(6.).reshape(3,2)
    >>> nd.array(x)
    nd.array([[0, 1], [2, 3], [4, 5]], type="strided * strided * float64")
    >>> nd.asarray(x)
    nd.array([[0, 1], [2, 3], [4, 5]], type="strided * strided * float64")
    >>> nd.view(x)
    nd.array([[0, 1], [2, 3], [4, 5]], type="strided * strided * float64")

    >>> x = np.array(['testing', 'one', 'two', 'three'])
    >>> nd.asarray(x)
    nd.array(["testing", "one", "two", "three"], type="strided * string[7,'utf32']")


Converting to NumPy Arrays
--------------------------

To support naturally feeding data into NumPy operations, the
NumPy array interface is used via the C struct PyArrayInterface.
This means NumPy operations will work on arrays with compatible
dtypes.

.. code-block:: python

    >>> x = nd.array([1, 2, 3.5])
    >>> np.square(x)
    array([  1.  ,   4.  ,  12.25])

There are some cases where a DyND array will not seamlessly convert
into a NumPy array. The behavior of NumPy is usually to ignore errors
that occur, and switch to an "object" array instead.

.. code-block:: python

    >>> x = nd.array([1, 2, 3]).ucast(ndt.float32)
    >>> x
    nd.array([1, 2, 3], type="strided * convert[to=float32, from=int32]")
    >>> np.array(x)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    TypeError: expected a readable buffer object
    >>> nd.as_numpy(x)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "_pydynd.pyx", line 1406, in _pydynd.as_numpy (_pydynd.cxx:9568)
    RuntimeError: cannot view dynd array with dtype strided_dim<convert<to=float32,
    from=int32>> as numpy without making a copy
    >>> nd.as_numpy(x, allow_copy=True)
    array([ 1.,  2.,  3.], dtype=float32)
        
.. DyND-Python documentation master file, created by
   sphinx-quickstart on Fri Aug 17 17:12:56 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the DyND Python bindings documentation!
===================================

.. warning::

   DyND is in a preview release stage, things are expected
   to change significantly during pre-1.0 development.

This is some preliminary documentation for DyND in Python, meant to provide
a simple reference to some of the key features in the library. Not everything
in the library is covered, and many things you might expect to work are
not implemented yet.

In any examples, it is assumed that ``from dynd import nd, ndt`` has been
executed prior to the sample Python code. In examples demonstrating
Numpy interoperability, ``import numpy as np`` is also assumed.

Contents:

.. toctree::
    :maxdepth: 2

    type
    array
    gfunc
    numpy-compat
    debugging
    dev
    dev-frompython

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

================================
DyND Python Bindings Development
================================

The DyND library, in its initial conception, is trying to
work well in two worlds, as both a C++ library and
a Python module. In both languages, DyND should feel
like a native first-class project, interacting with
the features of each language in the way heavy users
of the language will expect intuitively.
CPython is the present target for the Python bindings.

In Python, this means that the library should
exercise fine-grained control over how it interacts
with all the native Python types and interfaces.
NumPy accomplishes this by working with the CPython
API directly. For DyND, we did not want to go this
route because it requires a large amount of boilerplate
code.

Cython + CPython API for Bindings
---------------------------------

The approach chosen for DyND is to use Cython for
its ability to generate the CPython API boilerplate
for extension classes and other constructs, combined
with directly using the CPython API from C++ for
the rest. A minimal set of support classes and functions,
such as an RAII PyObject* class, are used to reduce
potential errors in accessing the API.

A number of techniques are used to make things integrate
smoothly and perform well. Let's go through them one
by one.

See http://docs.cython.org/src/userguide/wrapping_CPlusPlus.html
for the Cython documentation about wrapping C++.

Only Thin Wrappers in Cython
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This was a choice made after trying to do a few things in
Cython, and running into troubles with Cython's support
of C++ features like pointers, references, and operator
overloading. The Cython .pyx file mostly consists of
class definitions, docstrings, and calls to C++ functions.

Inline Functions
~~~~~~~~~~~~~~~~

For functionality we would like to inline into
the Cython functions, but would like to write in
C++ because Cython doesn't support the appropriate
syntax, or C++ is easier, we can use inline functions.

An example where this is necessary is to do
an assignment statement. The Cython developers debated
where to draw the line on this feature, whether
Cython should be more Python-like or C++-like, and
the result is that some code must be in C++ to access
C++ features.

Placement New/Delete
~~~~~~~~~~~~~~~~~~~~

The way the Cython documentation recommends that C++
objects be wrapped is to allocate them on the heap
and store a pointer in the Python object.

http://docs.cython.org/src/userguide/wrapping_CPlusPlus.html#create-cython-wrapper-class

Most objects in DyND behave in a smart pointer fashion,
so storing them with a heap allocation would add an
extra inefficient memory allocation and pointer indirection.

Cython does not allow C++ objects to be placed directly
in its data, only POD C structs. We can still accomplish
what we desire by manually calling the placement new
and destructor. In `placement_wrappers.hpp`, a set
of inline functions is defined to do placement new,
delete, get a reference to the C++ value, and
assign the C++ value.

This slightly obfuscates the Cython code, but the
names have been chosen in a way which is working
well in practice. Here's how the dtype class wrapper
starts::

    cdef class w_type:
        cdef ndt_type_placement_wrapper v

        def __cinit__(self, rep=None):
            placement_new(self.v)
            if rep is not None:
                SET(self.v, make_dtype_from_pyobject(rep))
        def __dealloc__(self):
            placement_delete(self.v)

Accessing Type Objects from C++
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For proper integration, its necessary to have access to
each wrapper class's type object, with corresponding
`<classname>_Check` function.

Cython has a mechanism to automatically generate a
header file, but in a test of generating a header for
just one of the wrapped classes, it generated a 1300 line
header, and still required calling a special `import_dynd`
function from the C++ side. This mechanism is designed for
allowing third party libraries to access the classes, but even
for that purpose a handwritten header seems much more appropriate,
as it would be human readable.

There is just one piece of information the C++ code needs to learn
from Cython, the TypeObject instance. Ideally one could declare an
exported variable `PyTypeObject *WArray_Type` in the Cython code,
and reference it with a one-liner in a header file. The Cython syntax
doesn't appear to support this, so the alternative being used is to
declare a C++ function::

    PyTypeObject *pydynd::WArray_Type;

    void pydynd::init_w_array_typeobject(PyObject *type)
    {
        WArray_Type = (PyTypeObject *)type;
    }

and call it from Cython at the outer scope, which is executed
on initialization::

    init_w_array_typeobject(w_array)
    
The full implementation of this in the C++ header is then::

    extern PyTypeObject *WArray_Type;
    inline bool WArray_CheckExact(PyObject *obj) {
        return Py_TYPE(obj) == WArray_Type;
    }
    inline bool WArray_Check(PyObject *obj) {
        return PyObject_TypeCheck(obj, WArray_Type);
    }
    struct WArray {
      PyObject_HEAD;
      // This is array_placement_wrapper in Cython-land
      dynd::array v;
    };
    void init_w_array_typeobject(PyObject *type);

There is a special consideration that must be made when constructing
the Cython classes from C++, which is that calling the
`WArray_Type->tp_alloc` method does not call the Cython
`__cinit__` function. This leads to the following wrapper code::

    inline PyObject *wrap_array(const dynd::array& n) {
        WArray *result = (WArray *)WArray_Type->tp_alloc(WArray_Type, 0);
        if (!result) {
            throw std::runtime_error("");
        }
        // Calling tp_alloc doesn't call Cython's __cinit__, so do the placement new here
        pydynd::placement_new(reinterpret_cast<pydynd::array_placement_wrapper &>(result->v));
        result->v = n;
        return (PyObject *)result;
    }

Translating C++ Exceptions to Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

http://docs.cython.org/src/userguide/wrapping_CPlusPlus.html#exceptions

Cython supports an automatic mechanism for translating
C++ exceptions into Python exceptions. The default
way to handle this translation is when declaring
functions imported from header files, to add 'except +'
to the end of the definition, as follows::

    void pydynd::translate_exception()
    {
        try {
            if (PyErr_Occurred())
                ; // let the latest Python exn pass through and ignore the current one
            else
                throw;
        } catch (const dynd::broadcast_error& exn) {
            PyErr_SetString(BroadcastException, exn.message());
        } catch (const dynd::too_many_indices& exn) {
            PyErr_SetString(PyExc_IndexError, exn.message());
        ...
        } catch (const std::exception& exn) {
            PyErr_SetString(PyExc_RuntimeError, exn.what());
        }
    }

The naked `throw` reraises the exception caught by the Cython code,
and uses an appropriate PyErr_SetString or PyErr_SetObject
to translate the exception. It appears that this is
conformant C++, as it is rethrowing the exception while within
scope of another catch statement, even though that statement
is within another function.

The only problem encountered is on Mac OS X, on an older version of
clang, where catching subclasses don't appear to work, and explicit
catches of every single possible exception was required. The solution
at the time was to switch to using g++ 4.2. This may have been caused
by a mismatch in libc++ vs. libstdc++ or something similar, but was
not tracked down.

Defining Custom Python Exceptions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The standard Python exceptions do not cover all the cases needed
by DyND, so we need to define some additional exception types
These new exceptions are defined in Cython, and their `TypeObject`
is passed to C++ in the same way as others. Here is the Cython
code for the `BroadcastError` class::

    # Exceptions to convert from C++
    class BroadcastError(Exception):
        pass

    # Register all the exception objects with the exception translator
    set_broadcast_exception(BroadcastError)

and the corresponding C++ code::

    PyObject *BroadcastException = NULL;

    void pydynd::set_broadcast_exception(PyObject *e)
    {
        BroadcastException = e;
    }

Accessing CTypes Structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~

CTypes doesn't define a C API for accessing its objects, so
to provide integration with CTypes requires some additional
work. This is done by defining a struct with all the needed
`PyTypeObject` instances from CTypes, and initializing them
using the CPython API at startup.::

    /**
     * Struct with data about the _ctypes module.
     */
    struct ctypes_info {
        // The _ctypes module (for C-implementation details)
        PyObject *_ctypes;
        // These match the corresponding names within _ctypes.c
        PyObject *PyCData_Type;
        PyObject *PyCStructType_Type;
        PyObject *UnionType_Type;
        PyObject *PyCPointerType_Type;
        PyObject *PyCArrayType_Type;
        PyObject *PyCSimpleType_Type;
        PyObject *PyCFuncPtrType_Type;
    };

    extern ctypes_info ctypes;

    /**
     * Should be called at module initialization, this
     * stores some internal information about the ctypes
     * classes for later.
     */
    void init_ctypes_interop();

PEP 3118 / Python Buffer Protocol
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For images and matrices, the buffer protocol introduced
in Python 2.6 is a great way to expose and consume regular
array data. This mechanism supports communicating multi-dimensional
arrays between C modules. For any module, like DyND, which exposes
images, matrices, or other multi-dimensional strided data,
supporting this is mandatory to interoperate properly with NumPy,
Cython, and other Python numerical libraries.

There are still some rough edges in the specification and
implementation. In the official Python 3.3 documentation, the
buffer protocol refers to the `struct` module for the specification
of the `format` string, but the `struct` module doesn't
include some of the additions proposed in PEP 3118, such as
complex numbers. Because most programs are still using this
protocol for low-level communication of arrays, supporting
16-bit floating point, complex, and other types specified
in PEP 3118 is possible without requiring the Python `struct`
module to support everything.

Debugging DyND
==============

One of the simplest tools to get started debugging and
understanding how DyND works is the ``nd.debug_repr(a)``
method which exists on most objects.

Here are a few examples to show what it prints.

Debug Printing ND.Array
-----------------------

.. code-block:: python

    >>> nd.debug_repr(nd.array(True))
    ------ array
     address: 0000000000415380
     refcount: 1
     type:
      pointer: 0000000000000001
      type: bool
     arrmeta:
      flags: 5 (read_access immutable )
     data:
       pointer: 00000000004153B0
       reference: 0000000000000000 (embedded in array memory)
    ------

    >>> nd.debug_repr(nd.array("testing"))
    ------ array
     address: 00000000003FAE90
     refcount: 1
     type:
      pointer: 0000000000415380
      type: string
     arrmeta:
      flags: 5 (read_access immutable )
      type-specific arrmeta:
       string arrmeta
        ------ NULL memory block
     data:
       pointer: 00000000003FAEC8
       reference: 0000000000000000 (embedded in array memory)
    ------

    >>> nd.debug_repr(nd.array([1,2,3,4,5]))
    ------ array
     address: 00000000003FAE90
     refcount: 1
     type:
      pointer: 000007FEE194E8F0
      type: strided * int32
     arrmeta:
      flags: 5 (read_access immutable )
      type-specific arrmeta:
       strided_dim arrmeta
        stride: 4
        size: 5
     data:
       pointer: 00000000003FAED0
       reference: 0000000000000000 (embedded in array memory)
    ------


Debug Printing GFuncs
---------------------

This functionality was disabled during code refactoring,
and will resurface at some point in the future.

.. code-block:: python

    >>> nd.maximum.debug_dump()

    >>> nd.sum.debug_dump()

===================
NumPy Compatibility
===================

While DyND doesn't aim to be a drop-in replacement for NumPy, having
a significant amount of compatibility is useful to help project
experiment with it. This document aims to track compatibility choices
that are made, both to be compatible and to have a different interface.

Type System
-----------

In NumPy, dtype and multi-dimensional array shape are separate properties
of the array.  DyND uses a different type system, developed within the
Blaze project and called datashape. Dimensions are a part of the types,
and the shape is typically encoded in the type as well.

Basic Constructors
------------------

The constructors ``nd.empty``, ``nd.zeros``, and ``nd.ones``, are very
similar to the NumPy functions of the same name. There is compatibility
in the sense that passing a shape as a tuple followed by a dtype works
in both. DyND is more flexible about the input parameters though.

The DyND constructors require a dtype be specified, they do not default
to float64 as in NumPy.

The ``nd.full`` constructor has an incompatible signature with NumPy.

Array Attributes
----------------

DyND arrays expose ``ndim``, ``shape``, ``strides``, and ``dtype``.
