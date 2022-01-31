
We welcome people who want to make contributions to Numba, big or small!
Even simple documentation improvements are encouraged.

# Asking questions

Numba has a [discourse forum](https://numba.discourse.group/) for longer/more
involved questions and an IRC channel on
[gitter.im](https://gitter.im/numba/numba) for quick questions and interactive
help.

# Ways to help:

There's lots of ways to help improve Numba, some of these require creating code
changes, see **contributing patches** below.

## Quick things:

* Answer a question asked on [discourse](https://numba.discourse.group/) or
  [gitter.im](https://gitter.im/numba/numba).
* Review a page of documentation, check it makes sense, that it's clear and
  still relevant, that the examples are present, good and working. Fix anything
  that needs updating in a pull request.
* Make a file that is not `flake8` compliant meet the standard, a list of all
  failing files is in the `exclude` section of the [`.flake8` config](https://github.com/numba/numba/blob/master/.flake8),
  then create a pull request with the change.

## More involved things:

* Review a pull request, you don't need to be a compiler engineer to do an
  initial review of a pull request. It's incredibly helpful to have pull
  requests go through a review to just make sure the code change is well formed,
  documented, efficient and clear. Further, if the code is fixing a bug, making
  sure that tests are present demonstrating it is fixed! Look out for PRs with
  the [`needs initial review`](https://github.com/numba/numba/labels/needs%20initial%20review)
  label.
* Work on fixing or implementing something in the code base, there are a lot of
  [`good first issue's`](https://github.com/numba/numba/labels/good%20first%20issue)
  and [`good second issue's`](https://github.com/numba/numba/labels/good%20first%20issue).
  For implementing new features/functionality, the extension API is the best
  thing to use and a guide to using `@overload` in particular is
  [here](https://numba.pydata.org/numba-doc/dev/extending/overloading-guide.html)
  and the API documentation is [here](https://numba.pydata.org/numba-doc/latest/extending/high-level.html#implementing-functions).

## Contributing patches

Please fork the Numba repository on Github, and create a new branch
containing your work.  When you are done, open a pull request.

# Further reading

Please read the [contributing guide](
https://numba.pydata.org/numba-doc/dev/developer/contributing.html).
<!--

Thanks for wanting to contribute to Numba :)

First, if you need some help or want to chat to the core developers, please
visit https://gitter.im/numba/numba for real time chat or post to the Numba
forum https://numba.discourse.group/.

Here's some guidelines to help the review process go smoothly.

0. Please write a description in this text box of the changes that are being
   made.

1. Please ensure that you have written units tests for the changes made/features
   added.

2. If you are closing an issue please use one of the automatic closing words as
   noted here: https://help.github.com/articles/closing-issues-using-keywords/

3. If your pull request is not ready for review but you want to make use of the
   continuous integration testing facilities here, please click the arrow besides
   "Create Pull Request" and choose "Create Draft Pull Request".
   When it's ready for review, you can click the button "ready to review" near
   the end of the pull request
   (besides "This pull request is still a work in progress".)
   The maintainers will then be automatically notified to review it.

4. Once review has taken place please do not add features or make changes out of
   the scope of those requested by the reviewer (doing this just add delays as
   already reviewed code ends up having to be re-reviewed/it is hard to tell
   what is new etc!). Further, please do not rebase your branch on master/force
   push/rewrite history, doing any of these causes the context of any comments
   made by reviewers to be lost. If conflicts occur against master they should
   be resolved by merging master into the branch used for making the pull
   request.

Many thanks in advance for your cooperation!

-->
---
name: Feature Request
about: Tell us about something in the Python language/NumPy you'd like Numba to support. Not for asking general questions - see below.

---

---

<!--

Thanks for opening an issue! To help the Numba team handle your information
efficiently, please first ensure that there is no other issue present that
already describes the issue you have
(search at https://github.com/numba/numba/issues?&q=is%3Aissue).

-->

## Feature request

<!--

Please include details of the feature you would like to see, why you would
like to see it/the use case.

-->
---
name: First Release Candidate Checklist (maintainer only)
about: Checklist template for the first release of every series
title: Numba X.Y.Zrc1 Checklist (FIXME)
labels: task

---


## Numba X.Y.Z

* [ ] Merge to master.
    - [ ] "remaining Pull-Requests from milestone".
* [ ] Review deprecation schedule and notices. Make PRs if need be.
* [ ] Merge change log changes.
    - [ ] "PR with changelog entries".
* [ ] Create X.Y release branch.
* [ ] Pin llvmlite to `>=0.A.0rc1,<0.A+1.0`.
* [ ] Pin NumPy if needed
* [ ] Pin tbb if needed
* [ ] Annotated tag X.Y.Zrc1 on release branch.
* [ ] Build and upload conda packages on buildfarm (check "upload").
* [ ] Build wheels (`$PYTHON_VERSIONS`) on the buildfarm.
* [ ] Verify packages uploaded to Anaconda Cloud and move to `numba/label/main`.
* [ ] Upload wheels and sdist to PyPI (upload from `ci_artifacts`).
* [ ] Verify wheels for all platforms arrived on PyPi.
* [ ] Initialize and verify ReadTheDocs build.
* [ ] Clean up `ci_artifacts`.
* [ ] Send RC announcement email / post announcement to discourse group.
* [ ] Post link to Twitter.

### Post Release:

* [ ] Tag X.Y+1.0dev0 to start new development cycle on `master`.
* [ ] Update llvmlite dependency spec to match next version via PR to `master`.
* [ ] Update release checklist template with any additional bullet points that
      may have arisen during the release.
* [ ] Close milestone (and then close this release issue).
---
name: Bug Report
about: Report a bug. Not for asking general questions - see below.

---

<!--

Thanks for opening an issue! To help the Numba team handle your information
efficiently, please first ensure that there is no other issue present that
already describes the issue you have
(search at https://github.com/numba/numba/issues?&q=is%3Aissue).

-->

## Reporting a bug

<!--

Before submitting a bug report please ensure that you can check off these boxes:

-->

- [ ] I have tried using the latest released version of Numba (most recent is
 visible in the change log (https://github.com/numba/numba/blob/master/CHANGE_LOG).
- [ ] I have included a self contained code sample to reproduce the problem.
  i.e. it's possible to run as 'python bug.py'.

<!--

Please include details of the bug here, including, if applicable, what you
expected to happen!

-->
---
name: Subsequent Release Candidate Checklist (maintainer only)
about: Checklist template for all subsequent releases (RC 2-N, FINAL and PATCH) of every series
title: Numba X.Y.Zrc1 Checklist (FIXME)
labels: task

---


## numba X.Y.Z

* [ ] Cherry-pick items from the X.Y.Z milestone into a PR.
* [ ] Approve change log modifications and cherry-pick.
* [ ] Merge change log modifications and cherry-picks to X.Y release branch.
  * [ ] https://github.com/numba/numba/pull/XXXX
* [ ] Review, merge and check execution of release notebook. (FINAL ONLY)
* [ ] Annotated tag X.Y.Z on release branch (no `v` prefix).
* [ ] Build and upload conda packages on buildfarm (check `upload`).
* [ ] Verify packages uploaded to Anaconda Cloud and move to
  `numba/label/main`.
* [ ] Build wheels (`$PYTHON_VERSIONS`) on the buildfarm.
* [ ] Upload wheels and sdist to PyPI (upload from `ci_artifacts`).
* [ ] Verify wheels for all platforms arrived on PyPi.
* [ ] Verify ReadTheDocs build.
* [ ] Send RC/FINAL announcement email / post announcement to discourse group.
* [ ] Post link to Twitter.
* [ ] Post link to python-announce-list@python.org.

### Post release

* [ ] Clean up `ci_artifacts` by moving files to subdirectories
* [ ] Update release checklist template.
* [ ] Ping Anaconda Distro team to trigger a build for `defaults` (FINAL ONLY).
* [ ] Create a release on Github at https://github.com/numba/numba/releases (FINAL ONLY).
* [ ] Close milestone (and then close this release issue).
# DAG Roadmap

This directory includes a representation of the Numba roadmap in the form of a
DAG.  We have done this to enable a highly granular display of enhancements to
Numba that also shows the relationships between these tasks. Many tasks have
prerequisites, and we've found that issue trackers, Kanban boards, and
time-bucketed roadmap documentation all fail to represent this information in
different ways.

## Requirements

```
conda install jinja2 python-graphviz pyyaml
```

## Usage

```
./render.py -o dagmap.html dagmap.yaml
```

The generated HTML file will look for `jquery.graphviz.svg.js` in the same
directory.

## Updating the DAG

Copy one of the existing tasks and edit:
  * `label`: text appears on the node.  Embed `\n` for line breaks.
  * `id`: Referenced to indicate a dependency
  * `description`: Shown in the tooltip.  Automatically word-wrapped.
  * `depends_on`: Optional list of task IDs which this task depends on.

The `style` section of the file is not used yet.

## Notes

The HTML rendering of the graph is based on a slightly modified version of
(jquery.graphviz.svg)[https://github.com/mountainstorm/jquery.graphviz.svg/].
Its license is:
```
Copyright (c) 2015 Mountainstorm
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```*****
Numba
*****

.. image:: https://badges.gitter.im/numba/numba.svg
   :target: https://gitter.im/numba/numba?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
   :alt: Gitter

.. image:: https://img.shields.io/badge/discuss-on%20discourse-blue
   :target: https://numba.discourse.group/
   :alt: Discourse

.. image:: https://zenodo.org/badge/3659275.svg
   :target: https://zenodo.org/badge/latestdoi/3659275
   :alt: Zenodo DOI
   
.. image:: https://img.shields.io/pypi/v/numba.svg
   :target: https://pypi.python.org/pypi/numba/
   :alt: PyPI

A Just-In-Time Compiler for Numerical Functions in Python
#########################################################

Numba is an open source, NumPy-aware optimizing compiler for Python sponsored
by Anaconda, Inc.  It uses the LLVM compiler project to generate machine code
from Python syntax.

Numba can compile a large subset of numerically-focused Python, including many
NumPy functions.  Additionally, Numba has support for automatic
parallelization of loops, generation of GPU-accelerated code, and creation of
ufuncs and C callbacks.

For more information about Numba, see the Numba homepage:
https://numba.pydata.org

Supported Platforms
===================

* Operating systems and CPUs:

  - Linux: x86 (32-bit), x86_64, ppc64le (POWER8 and 9), ARMv7 (32-bit),
    ARMv8 (64-bit).
  - Windows: x86, x86_64.
  - macOS: x86_64, (M1/Arm64, unofficial support only).
  - \*BSD: (unofficial support only).

* (Optional) Accelerators and GPUs:

  * NVIDIA GPUs (Kepler architecture or later) via CUDA driver on Linux and
    Windows.

Dependencies
============

* Python versions: 3.7-3.10
* llvmlite 0.39.*
* NumPy >=1.18 (can build with 1.11 for ABI compatibility).

Optionally:

* SciPy >=1.0.0 (for ``numpy.linalg`` support).


Installing
==========

The easiest way to install Numba and get updates is by using the Anaconda
Distribution: https://www.anaconda.com/download

::

   $ conda install numba

For more options, see the Installation Guide:
https://numba.readthedocs.io/en/stable/user/installing.html

Documentation
=============

https://numba.readthedocs.io/en/stable/index.html


Contact
=======

Numba has a discourse forum for discussions:

* https://numba.discourse.group



Continuous Integration
======================

.. image:: https://dev.azure.com/numba/numba/_apis/build/status/numba.numba?branchName=master
    :target: https://dev.azure.com/numba/numba/_build/latest?definitionId=1?branchName=master
    :alt: Azure Pipelines
======================
Release Notes
======================

.. include:: ../../CHANGE_LOG

Glossary
========

.. glossary::

   ahead-of-time compilation
   AOT compilation
   AOT
      Compilation of a function in a separate step before running the
      program code, producing an on-disk binary object which can be distributed
      independently.  This is the traditional kind of compilation known
      in languages such as C, C++ or Fortran.

   bytecode
   Python bytecode
      The original form in which Python functions are executed.  Python
      bytecode describes a stack-machine executing abstract (untyped)
      operations using operands from both the function stack and the
      execution environment (e.g. global variables).

   compile-time constant
      An expression whose value Numba can infer and freeze at compile-time.
      Global variables and closure variables are compile-time constants.

   just-in-time compilation
   JIT compilation
   JIT
      Compilation of a function at execution time, as opposed to
      :term:`ahead-of-time compilation`.

   JIT function
      Shorthand for "a function :term:`JIT-compiled <JIT>` with Numba using
      the :ref:`@jit <jit>` decorator."

   loop-lifting
   loop-jitting
      A feature of compilation in :term:`object mode` where a loop can be
      automatically extracted and compiled in :term:`nopython mode`.  This
      allows functions with operations unsupported in nopython mode to see
      significant performance improvements if they contain loops with only
      nopython-supported operations.

   lowering
      The act of translating :term:`Numba IR` into LLVM IR.  The term
      "lowering" stems from the fact that LLVM IR is low-level and
      machine-specific while Numba IR is high-level and abstract.

   NPM
   nopython mode
      A Numba compilation mode that generates code that does not access the
      Python C API.  This compilation mode produces the highest performance
      code, but requires that the native types of all values in the function
      can be :term:`inferred <type inference>`.  Unless otherwise instructed,
      the ``@jit`` decorator will automatically fall back to :term:`object
      mode` if nopython mode cannot be used.

   Numba IR
   Numba intermediate representation
      A representation of a piece of Python code which is more amenable
      to analysis and transformations than the original Python
      :term:`bytecode`.

   object mode
      A Numba compilation mode that generates code that handles all values
      as Python objects and uses the Python C API to perform all operations
      on those objects.  Code compiled in object mode will often run
      no faster than Python interpreted code, unless the Numba compiler can
      take advantage of :term:`loop-jitting`.

   ``OptionalType``
     An ``OptionalType`` is effectively a type union of a ``type`` and ``None``.
     They typically occur in practice due to a variable being set to ``None``
     and then in a branch the variable being set to some other value. It's
     often not possible at compile time to determine if the branch will execute
     so to permit :term:`type inference` to complete, the type of the variable
     becomes the union of a ``type`` (from the value) and ``None``,
     i.e. ``OptionalType(type)``.

   type inference
      The process by which Numba determines the specialized types of all
      values within a function being compiled.  Type inference can fail
      if arguments or globals have Python types unknown to Numba, or if
      functions are used that are not recognized by Numba.  Successful
      type inference is a prerequisite for compilation in
      :term:`nopython mode`.

   typing
      The act of running :term:`type inference` on a value or operation.

   ufunc
      A NumPy `universal function <http://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_.
      Numba can create new compiled ufuncs with
      the :ref:`@vectorize <vectorize>` decorator.

   reflection
      In numba, when a mutable container is passed as argument to a nopython
      function from the Python interpreter, the container object and all its
      contained elements are converted into nopython values.  To match the
      semantics of Python, any mutation on the container inside the nopython
      function must be visible in the Python interpreter.  To do so, Numba
      must update the container and its elements and convert them back into
      Python objects during the transition back into the interpreter.

      Not to be confused with Python's "reflection" in the context of binary
      operators (see https://docs.python.org/3.5/reference/datamodel.html).
.. Numba documentation master file, created by
   sphinx-quickstart on Tue Dec 30 11:55:40 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Numba documentation
===================

This is the Numba documentation.  Unless you are already acquainted
with Numba, we suggest you start with the :doc:`User manual <user/index>`.


.. toctree::
   :caption: For all users
   :maxdepth: 2

   user/index.rst
   reference/index.rst


.. toctree::
   :caption: For CUDA users
   :maxdepth: 2

   cuda/index.rst
   cuda-reference/index.rst


.. toctree::
   :caption: For advanced users & developers
   :maxdepth: 2

   extending/index.rst
   developer/index.rst
   proposals/index.rst
   glossary.rst
   release-notes.rst
===================
NBEP 3: JIT Classes
===================

:Author: Siu Kwan Lam
:Date: Dec 2015
:Status: Draft

Introduction
============

Numba does not yet support user-defined classes.
Classes provide useful abstraction and promote modularity when used
right.  In the simplest sense, a class specifies the set of data and
operations as attributes and methods, respectively.
A class instance is an instantiation of that class.
This proposal will focus on supporting this simple usecase of classes--with
just attributes and methods.  Other features, such as class methods, static
methods, and inheritance are deferred to another proposal, but we believe
these features can be easily implemented given the foundation described here.


Proposal: jit-classes
=====================

A JIT-classes is more restricted than a Python class.
We will focus on the following operations on a class and its instance:

* Instantiation: create an instance of a class using the class object as the
  constructor: ``cls(*args, **kwargs)``
* Destruction: remove resources allocated during instantiation and release
  all references to other objects.
* Attribute access: loading and storing attributes using ``instance.attr``
  syntax.
* Method access: loading methods using ``instance.method`` syntax.

With these operations, a class object (not the instance) does not need to be
materialize. Using the class object as a constructor is fully resolved (a
runtime implementation is picked) during the typing phase in the compiler.
This means **a class object will not be first class**.  On the other hand,
implementating a first-class class object will require an
"interface" type, or the type of class.

The instantiation of a class will allocate resources for storing the data
attributes.  This is described in the "Storage model" section.  Methods are
never stored in the instance.  They are information attached to the class.
Since a class object only exists in the type domain, the methods will also be
fully resolved at the typing phase.  Again, numba do not have first-class
function value and each function type maps uniquely to each function
implementation (this needs to be changed to support function value as argument).

A class instance can contain other NRT reference-counted object as attributes.
To properly clean up an instance, a destructor is called when the reference
count of the instance is dropped to zero.  This is described in the
"Reference count and descructor" section.

Storage model
~~~~~~~~~~~~~

For compatibility with C, attributes are stored in a simple plain-old-data
structure.  Each attribute are stored in a user-defined order in a padded
(for proper alignment), contiguous memory region. An instance that contains
three fields of int32, float32, complex64 will be compatible with the following
C structure::

    struct {
        int32     field0;
        float32   field1;
        complex64 field2;
    };

This will also be comptabile with an aligned numpy structure dtype.


Methods
~~~~~~~

Methods are regular function that can be bounded to an instance.
They can be compiled as regular function by numba.
The operation ``getattr(instance, name)`` (getting an attribute ``name`` from
``instance``) binds the instance to the requested method at runtime.


The special ``__init__`` method is also handled like regular functions.


``__del__`` is not supported at this time.


Reference count and destructor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An instance of jit-class is reference-counted by NRT. Since it may contain
other NRT tracked object, it must call a destructor when its reference count
dropped to zero.  The destructor will decrement the reference count of all
attributes by one.

At this time, there is no support for user defined ``__del__`` method.

Proper cleanup for cyclic reference is not handled at this time.
Cycles will cause memory leak.

Type inference
~~~~~~~~~~~~~~

So far we have not described the type of the attributes or the methods.
Type information is necessary to materailize the instance (e.g. allocate the
storage).  The simplest way is to let user provide the type of each attributes
as well as the ordering; for instance::

    dct = OrderedDict()
    dct['x'] = int32
    dct['y'] = float32

Allowing user to supply an ordered dictionary will provide the name, ordering
and types of the attributes.  However, this statically typed semantic is not as
flexible as the Python semantic which behaves like a generic class.

Inferring the type of attributes is difficult.  In a previous attempt to
implement JIT classes, the ``__init__`` method is specialized to capture
the type stored into the attributes.  Since the method can contain arbitrary
logic, the problem can become a dependent typing problem if types are assigned
conditionally depending on the value. (Very few languages implement dependent
typing and those that does are mostly theorem provers.)

Example: typing function using an OrderedDict
---------------------------------------------

.. code-block:: python

    spec = OrderedDict()
    spec['x'] = numba.int32
    spec['y'] = numba.float32

    @jitclass(spec)
    class Vec(object):
        def __init__(self, x, y):
            self.x = x
            self.y = y

        def add(self, dx, dy):
            self.x += dx
            self.y += dy

Example: typing function using a list of 2-tuples
-------------------------------------------------

.. code-block:: python

    spec = [('x', numba.int32),
            ('y', numba.float32)]

    @jitclass(spec)
    class Vec(object):
        ...

Creating multiple jitclasses from a single class object
-------------------------------------------------------

The `jitclass(spec)` decorator creates a new jitclass type even when applied to
the same class object and the same type specification.

.. code-block:: python

    class Vec(object):
      ...

    Vec1 = jitclass(spec)(Vec)
    Vec2 = jitclass(spec)(Vec)
    # Vec1 and Vec2 are two different jitclass types

Usage from the Interpreter
~~~~~~~~~~~~~~~~~~~~~~~~~~

When constructing a new instance of a jitclass, a "box" is created that wraps
the underlying jitclass instance from numba.  Attributes and methods are
accessible from the interpreter.  The actual implementation will be in numba
compiled code.  Any Python object is converted to its native
representation for consumption in numba.  Similarly, the returned value is
converted to its Python representation.  As a result, there may be overhead in
manipulating jitclass instances in the interpreter.  This overhead is minimal
and should be easily amortized by more efficient computation in the compiled
methods.

Support for property, staticmethod and classmethod
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The use of ``property`` is accepted for getter and setter only.  Deleter is not
supported.

The use of ``staticmethod`` is not supported.

The use of ``classmethod`` is not supported.

Inheritance
~~~~~~~~~~~

Class inhertance is not considered in this proposal.  The only accepted base
class for a jitclass is `object`.

Supported targets
~~~~~~~~~~~~~~~~~~

Only the CPU target (including the parallel target) is supported.
GPUs (e.g. CUDA and HSA) targets are supported via an immutable version of the
jitclass instance, which will be described in a separate NBEP.


Other properties
~~~~~~~~~~~~~~~~

Given:

.. code-block:: python

    spec = [('x', numba.int32),
            ('y', numba.float32)]

    @jitclass(spec)
    class Vec(object):
        ...

* ``isinstance(Vec(1, 2), Vec)`` is True.
* ``type(Vec(1, 2))`` may not be ``Vec``.

Future enhancements
~~~~~~~~~~~~~~~~~~~

This proposal has only described the basic semantic and functionality of a
jitclass.  Additional features will be described in future enhancement
proposals.
========================
NBEP 6: Typing Recursion
========================

:Author: Siu Kwan Lam
:Date: Sept 2016
:Status: Draft

Introduction
============

This document proposes an enhancement to the type inference algorithm to
support recursion without explicitly annotating the function signature.
As a result, the proposal enables numba to type-infer both self-recursive and
mutual-recursive functions under some limitations.  In practice, these
limitions can be easily overcome by specifying a compilation order.


The Current State
=================

Recursion support in numba is currently limited to self-recursion with explicit
type annotation for the function.  This limitation comes from the inability to
determine the return type of a recursive call.  This is because the callee is
either the current function (for self-recursion) or a parent function
(mutual-recursion) and its type inference process has been suspended while waiting for
the function-type of its callee.  This results in the formation of a cyclic
dependency.  For example, given a function ``foo()`` that calls ``bar()``,
which in turns call ``foo()``::

    def foo(x):
        if x > 0:
            return bar(x)
        else:
            return 1

    def bar(x):
        return foo(x - 1)


The type inferrence process of ``foo()`` depends on that of ``bar()``,
which depends on ``foo()``.  Therefore ``foo()`` depends on itself and the type
inference algorithm cannot terminate.


The Solution
============

The proposed solution has two components:

1. The introduction of a compile-time *callstack* that tracks the compiling functions.
2. The allowance of a partial type inference on functions by leveraging the return type
   on non-recursive control-flow paths.

The compile-time callstack stores typing information of the functions being
compiled.  Like an ordinary callstack, it pushes a new record every time a
function is "called".  Since this occurs at compile-time, a "call" triggers
a compilation of the callee.

To detect recursion, the compile-time callstack is searched bottom-up
(stack grows downward) for a record that matches the callee.
As the record contains a reference to the type inference state,
the type inference process can be resumed to determine the return type.

Recall that the type inference process cannot be resumed normally because of the cyclic
dependency of the return type.  In practice, we can assume that a useful
program must have a terminating condition, a path that does not recurse.  So,
the type inference process can make an initial guess for the return-type at the recursive
call by using the return-type determined by the non-recursive paths.  This
allows type information to propagate on the recursive paths to generate the
final return type, which is used to refine the type information by the
subsequent iteration in the type inference process.

The following figure illustrates the compile-time callstack when the compiler
reaches the recursive call to ``foo()`` from ``bar()``:

.. image:: recursion_callstack.svg
    :width: 400px

At this time, the type inference process of ``foo()`` is suspended and that of ``bar()``
is active.  The compiler can see that the callee is already compiling by
searching the callstack.  Knowing that it is a recursive call, the compiler
can resume the type-inference on ``foo()`` by ignoring the paths that contain
recursive calls.  This means only the ``else`` branch is considered and we can
easily tell that ``foo()`` returns an ``int`` in this case.  The compiler will
then set the initial return type of ``foo()`` and ``bar()`` to ``int``.  The
subsequent type propagation can use this information to complete the type
inference of both functions, unifying the return-type of all returning paths.


Limitations
===========

For the proposed type inference algorithm to terminate, it assumes that
at least one of the control path leads to a return-statement without undertaking
a recursive call.  Should this not be the case, the algorithm will raise an
exception indicating a potential runaway recursion.

For example::

    @jit
    def first(x):
        # The recursing call must have a path that is non-recursing.
        if x > 0:
            return second(x)
        else:
            return 1

    @jit
    def second(x):
        return third(x)

    @jit
    def third(x):
        return first(x - 1)


The ``first()`` function must be the compiled first for the type inference algorithm to
complete successfully.  Compiling any other function first will lead to a failure
in type inference.  The type inference algorithm will treat it as a runaway
recursion due to the lack of a non-recursive exit in the recursive callee.

For example, compiling ``second()`` first will move the recursive call to
``first()``.  When the compiler tries to resume the type inference process of
``second()``, it will fail to find a non-recursive path.

This is a small limitation and can be overcome easily by code restructuring or
precompiling in a specific order.

.. _nbep-7:

===============================================
NBEP 7: CUDA External Memory Management Plugins
===============================================

:Author: Graham Markall, NVIDIA
:Contributors: Thomson Comer, Peter Entschev, Leo Fang, John Kirkham, Keith Kraus
:Date: March 2020
:Status: Final

Background and goals
--------------------

The :ref:`CUDA Array Interface <cuda-array-interface>` enables sharing of data
between different Python libraries that access CUDA devices. However, each
library manages its own memory distinctly from the others. For example:


* `Numba <https://numba.pydata.org/>`_ internally manages memory for the creation
  of device and mapped host arrays.
* `The RAPIDS libraries <https://rapids.ai/>`_ (cuDF, cuML, etc.) use the `Rapids
  Memory Manager <https://github.com/rapidsai/rmm>`_ for allocating device
  memory.
* `CuPy <https://cupy.chainer.org/>`_ includes a `memory pool
  implementation <https://docs-cupy.chainer.org/en/stable/reference/memory.html>`_
  for both device and pinned memory.

The goal of this NBEP is to describe a plugin interface that enables Numba's
internal memory management to be replaced with an external memory manager by the
user. When the plugin interface is in use, Numba no longer directly allocates or
frees any memory when creating arrays, but instead requests allocations and
frees through the external manager.

Requirements
------------

Provide an *External Memory Manager (EMM)* interface in Numba.


* When the EMM is in use, Numba will make all memory allocation using the EMM.
  It will never directly call functions such as ``CuMemAlloc``\ , ``cuMemFree``\ , etc.
* When not using an *External Memory Manager (EMM)*\ , Numba's present behaviour
  is unchanged (at the time of writing, the current version is the 0.48
  release).

If an EMM is to be used, it will entirely replace Numba's internal memory
management for the duration of program execution. An interface for setting the
memory manager will be provided.

Device vs. Host memory
^^^^^^^^^^^^^^^^^^^^^^^

An EMM will always take responsibility for the management of device memory.
However, not all CUDA memory management libraries also support managing host
memory, so a facility for Numba to continue the management of host memory 
whilst ceding control of device memory to the EMM will be provided.

Deallocation strategies
^^^^^^^^^^^^^^^^^^^^^^^

Numba's internal memory management uses a :ref:`deallocation strategy
<deallocation-behavior>` designed to increase efficiency by deferring
deallocations until a significant quantity are pending. It also provides a
mechanism for preventing deallocations entirely during critical sections, using
the :func:`~numba.cuda.defer_cleanup` context manager.


* When the EMM is not in use, the deallocation strategy and operation of
  ``defer_cleanup`` remain unchanged.
* When the EMM is in use, the deallocation strategy is implemented by the EMM,
  and Numba's internal deallocation mechanism is not used. For example:

  * A similar strategy to Numba's could be implemented by the EMM, or
  * Deallocated memory might immediately be returned to a memory pool.

* The ``defer_cleanup`` context manager may behave differently with an EMM - an
  EMM should be accompanied by documentation of the behaviour of the
  ``defer_cleanup`` context manager when it is in use.

  * For example, a pool allocator could always immediately return memory to a
    pool even when the context manager is in use, but could choose
    not to free empty pools until ``defer_cleanup`` is not in use.

Management of other objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to memory, Numba manages the allocation and deallocation of
:ref:`events <events>`, :ref:`streams <streams>`, and modules (a module is a
compiled object, which is generated from ``@cuda.jit``\ -ted functions). The
management of streams, events, and modules should be unchanged by the presence
or absence of an EMM.

Asynchronous allocation / deallocation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An asynchronous memory manager might provide the facility for an allocation or
free to take a CUDA stream and execute asynchronously. For freeing, this is
unlikely to cause issues since it operates at a layer beneath Python, but for
allocations this could be problematic if the user tries to then launch a kernel
on the default stream from this asynchronous memory allocation.

The interface described in this proposal will not be required to support
asynchronous allocation and deallocation, and as such these use cases will not
be considered further. However, nothing in this proposal should preclude the
straightforward addition of asynchronous operations in future versions of the
interface.

Non-requirements
^^^^^^^^^^^^^^^^

In order to minimise complexity and constrain this proposal to a reasonable
scope, the following will not be supported:


* Using different memory manager implementations for different contexts. All
  contexts will use the same memory manager implementation - either the Numba
  internal implementation or an external implementation.
* Changing the memory manager once execution has begun. It is not practical to
  change the memory manager and retain all allocations. Cleaning up the entire
  state and then changing to a different memory allocator (rather than starting
  a new process) appears to be a rather niche use case.
* Any changes to the ``__cuda_array_interface__`` to further define its semantics,
  e.g. for acquiring / releasing memory as discussed in `Numba Issue
  #4886 <https://github.com/numba/numba/issues/4886>`_ - these are independent,
  and can be addressed as part of separate proposals.
* Managed memory / UVM is not supported. At present Numba does not support UVM -
  see `Numba Issue #4362 <https://github.com/numba/numba/issues/4362>`_ for
  discussion of support.

Interface for Plugin developers
-------------------------------

New classes and functions will be added to ``numba.cuda.cudadrv.driver``:

* ``BaseCUDAMemoryManager`` and ``HostOnlyCUDAMemoryManager``\ : base classes for
  EMM plugin implementations.
* ``set_memory_manager``: a method for registering an external memory manager with
  Numba.

These will be exposed through the public API, in the ``numba.cuda`` module.
Additionally, some classes that are already part of the `driver` module will be
exposed as part of the public API:

* ``MemoryPointer``: used to encapsulate information about a pointer to device
  memory.
* ``MappedMemory``: used to hold information about host memory that is mapped into
  the device address space (a subclass of ``MemoryPointer``\ ).
* ``PinnedMemory``: used to hold information about host memory that is pinned (a
  subclass of ``mviewbuf.MemAlloc``\ , a class internal to Numba).

As an alternative to calling the ``set_memory_manager`` function, an environment
variable can be used to set the memory manager. The value of the environment
variable should be the name of the module containing the memory manager in its
global scope, named ``_numba_memory_manager``\ :

.. code-block::

   export NUMBA_CUDA_MEMORY_MANAGER="<module>"

When this variable is set, Numba will automatically use the memory manager from
the specified module. Calls to ``set_memory_manager`` will issue a warning, but
otherwise be ignored.

Plugin Base Classes
^^^^^^^^^^^^^^^^^^^

An EMM plugin is implemented by inheriting from the ``BaseCUDAMemoryManager``
class, which is defined as:

.. code-block:: python

   class BaseCUDAMemoryManager(object, metaclass=ABCMeta):
       @abstractmethod
       def memalloc(self, size):
           """
           Allocate on-device memory in the current context. Arguments:

           - `size`: Size of allocation in bytes

           Returns: a `MemoryPointer` to the allocated memory.
           """

       @abstractmethod
       def memhostalloc(self, size, mapped, portable, wc):
           """
           Allocate pinned host memory. Arguments:

           - `size`: Size of the allocation in bytes
           - `mapped`: Whether the allocated memory should be mapped into the CUDA
                       address space.
           - `portable`: Whether the memory will be considered pinned by all
                         contexts, and not just the calling context.
           - `wc`: Whether to allocate the memory as write-combined.

           Returns a `MappedMemory` or `PinnedMemory` instance that owns the
           allocated memory, depending on whether the region was mapped into
           device memory.
           """

       @abstractmethod
       def mempin(self, owner, pointer, size, mapped):
           """
           Pin a region of host memory that is already allocated. Arguments:

           - `owner`: An object owning the memory - e.g. a `DeviceNDArray`.
           - `pointer`: The pointer to the beginning of the region to pin.
           - `size`: The size of the region to pin.
           - `mapped`: Whether the region should also be mapped into device memory.

           Returns a `MappedMemory` or `PinnedMemory` instance that refers to the
           allocated memory, depending on whether the region was mapped into device
           memory.
           """

       @abstractmethod
       def initialize(self):
           """
           Perform any initialization required for the EMM plugin to be ready to
           use.
           """

       @abstractmethod
       def get_memory_info(self):
           """
           Returns (free, total) memory in bytes in the context
           """

       @abstractmethod
       def get_ipc_handle(self, memory):
           """
           Return an `IpcHandle` from a GPU allocation. Arguments:

           - `memory`: A `MemoryPointer` for which the IPC handle should be created.
           """

       @abstractmethod
       def reset(self):
           """
           Clear up all memory allocated in this context.
           """

       @abstractmethod
       def defer_cleanup(self):
           """
           Returns a context manager that ensures the implementation of deferred
           cleanup whilst it is active.
           """

       @property
       @abstractmethod
       def interface_version(self):
           """
           Returns an integer specifying the version of the EMM Plugin interface
           supported by the plugin implementation. Should always return 1 for
           implementations described in this proposal.
           """

All of the methods of an EMM plugin are called from within Numba - they never
need to be invoked directly by a Numba user.

The ``initialize`` method is called by Numba prior to any memory allocations
being requested. This gives the EMM an opportunity to initialize any data
structures, etc., that it needs for its normal operations. The method may be
called multiple times during the lifetime of the program - subsequent calls
should not invalidate or reset the state of the EMM.

The ``memalloc``\ , ``memhostalloc``\ , and ``mempin`` methods are called when Numba
requires an allocation of device or host memory, or pinning of host memory.
Device memory should always be allocated in the current context.

``get_ipc_handle`` is called when an IPC handle for an array is required. Note
that there is no method for closing an IPC handle - this is because the
``IpcHandle`` object constructed by ``get_ipc_handle`` contains a ``close()`` method
as part of its definition in Numba, which closes the handle by calling
``cuIpcCloseMemHandle``. It is expected that this is sufficient for general use
cases, so no facility for customising the closing of IPC handles is provided by
the EMM Plugin interface.

``get_memory_info`` may be called at any time after ``initialize``.

``reset`` is called as part of resetting a context. Numba does not normally call
reset spontaneously, but it may be called at the behest of the user. Calls to
``reset`` may even occur before ``initialize`` is called, so the plugin should be
robust against this occurrence.

``defer_cleanup`` is called when the ``numba.cuda.defer_cleanup`` context manager
is used from user code.

``interface_version`` is called by Numba when the memory manager is set, to
ensure that the version of the interface implemented by the plugin is
compatible with the version of Numba in use.

Representing pointers
^^^^^^^^^^^^^^^^^^^^^

Device Memory
~~~~~~~~~~~~~

The ``MemoryPointer`` class is used to represent a pointer to memory. Whilst there
are various details of its implementation, the only aspect relevant to EMM
plugin development is its initialization. The ``__init__`` method has the
following interface:

.. code-block:: python

   class MemoryPointer:
       def __init__(self, context, pointer, size, owner=None, finalizer=None):


* ``context``\ : The context in which the pointer was allocated.
* ``pointer``\ : A ``ctypes`` pointer (e.g. ``ctypes.c_uint64``\ ) holding the address of
  the memory.
* ``size``\ : The size of the allocation in bytes.
* ``owner``\ : The owner is sometimes set by the internals of the class, or used for
  Numba's internal memory management, but need not be provided by the writer of
  an EMM plugin - the default of ``None`` should always suffice.
* ``finalizer``\ : A method that is called when the last reference to the
  ``MemoryPointer`` object is released. Usually this will make a call to the
  external memory management library to inform it that the memory is no longer
  required, and that it could potentially be freed (though the EMM is not
  required to free it immediately).

Host Memory
~~~~~~~~~~~

Memory mapped into the CUDA address space (which is created when the
``memhostalloc`` or ``mempin`` methods are called with ``mapped=True``\ ) is managed
using the ``MappedMemory`` class:

.. code-block:: python

   class MappedMemory(AutoFreePointer):
       def __init__(self, context, pointer, size, owner, finalizer=None):


* ``context``\ : The context in which the pointer was allocated.
* ``pointer``\ : A ``ctypes`` pointer (e.g. ``ctypes.c_void_p``\ ) holding the address of
  the allocated memory.
* ``size``\ : The size of the allocated memory in bytes.
* ``owner``\ : A Python object that owns the memory, e.g. a ``DeviceNDArray``
  instance.
* ``finalizer``\ : A method that is called when the last reference to the
  ``MappedMemory`` object is released. For example, this method could call
  ``cuMemFreeHost`` on the pointer to deallocate the memory immediately.

Note that the inheritance from ``AutoFreePointer`` is an implementation detail and
need not concern the developer of an EMM plugin - ``MemoryPointer`` is higher in
the MRO of ``MappedMemory``.

Memory that is only in the host address space and has been pinned is represented
with the ``PinnedMemory`` class:

.. code-block:: python

   class PinnedMemory(mviewbuf.MemAlloc):
       def __init__(self, context, pointer, size, owner, finalizer=None):


* ``context``\ : The context in which the pointer was allocated.
* ``pointer``\ : A ``ctypes`` pointer (e.g. ``ctypes.c_void_p``\ ) holding the address of
  the pinned memory.
* ``size``\ : The size of the pinned region in bytes.
* ``owner``\ : A Python object that owns the memory, e.g. a ``DeviceNDArray``
  instance.
* ``finalizer``\ : A method that is called when the last reference to the
  ``PinnedMemory`` object is released. This method could e.g. call
  ``cuMemHostUnregister`` on the pointer to unpin the memory immediately.

Providing device memory management only
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some external memory managers will support management of on-device memory but
not host memory. To make it easy to implement an EMM plugin using one of these
managers, Numba will provide a memory manager class with implementations of the
``memhostalloc`` and ``mempin`` methods. An abridged definition of this class
follows:

.. code-block:: python

   class HostOnlyCUDAMemoryManager(BaseCUDAMemoryManager):
       # Unimplemented methods:
       #
       # - memalloc
       # - get_memory_info

       def memhostalloc(self, size, mapped, portable, wc):
           # Implemented.

       def mempin(self, owner, pointer, size, mapped):
           # Implemented.

       def initialize(self):
           # Implemented.
           #
           # Must be called by any subclass when its initialize() method is
           # called.

       def reset(self):
           # Implemented.
           #
           # Must be called by any subclass when its reset() method is
           # called.

       def defer_cleanup(self):
           # Implemented.
           #
           # Must be called by any subclass when its defer_cleanup() method is
           # called.

A class can subclass the ``HostOnlyCUDAMemoryManager`` and then it only needs to
add implementations of methods for on-device memory. Any subclass must observe
the following rules:


* If the subclass implements ``__init__``\ , then it must also call
  ``HostOnlyCUDAMemoryManager.__init__``\ , as this is used to initialize some of
  its data structures (\ ``self.allocations`` and ``self.deallocations``\ ).
* The subclass must implement ``memalloc`` and ``get_memory_info``.
* The ``initialize`` and ``reset`` methods perform initialisation of structures
  used by the ``HostOnlyCUDAMemoryManager``.

  * If the subclass has nothing to do on initialisation (possibly) or reset
    (unlikely) then it need not implement these methods. 
  * However, if it does implement these methods then it must also call the
    methods from ``HostOnlyCUDAMemoryManager`` in its own implementations.

* Similarly if ``defer_cleanup`` is implemented, it should enter the context
  provided by ``HostOnlyCUDAManager.defer_cleanup()`` prior to ``yield``\ ing (or in
  the ``__enter__`` method) and release it prior to exiting (or in the ``__exit__``
  method).

Import order
^^^^^^^^^^^^

The order in which Numba and the library implementing an EMM Plugin should not
matter. For example, if ``rmm`` were to implement and register an EMM Plugin,
then:

.. code-block:: python

   from numba import cuda
   import rmm

and

.. code-block:: python

   import rmm
   from numba import cuda

are equivalent - this is because Numba does not initialize CUDA or allocate any
memory until the first call to a CUDA function - neither instantiating and
registering an EMM plugin, nor importing ``numba.cuda`` causes a call to a CUDA
function.

Numba as a Dependency
^^^^^^^^^^^^^^^^^^^^^

Adding the implementation of an EMM Plugin to a library naturally makes Numba a
dependency of the library where it may not have been previously. In order to
make the dependency optional, if this is desired, one might conditionally
instantiate and register the EMM Plugin like:

.. code-block:: python

   try:
       import numba
       from mylib.numba_utils import MyNumbaMemoryManager
       numba.cuda.cudadrv.driver.set_memory_manager(MyNumbaMemoryManager)
   except:
       print("Numba not importable - not registering EMM Plugin")

so that ``mylib.numba_utils``\ , which contains the implementation of the EMM
Plugin, is only imported if Numba is already present. If Numba is not available,
then ``mylib.numba_utils`` (which necessarily imports ``numba``\ ), will never be
imported.

It is recommended that any library with an EMM Plugin includes at least some
environments with Numba for testing with the EMM Plugin in use, as well as some
environments without Numba, to avoid introducing an accidental Numba dependency.

Example implementation - A RAPIDS Memory Manager (RMM) Plugin
-------------------------------------------------------------

An implementation of an EMM plugin within the `Rapids Memory Manager
(RMM) <https://github.com/rapidsai/rmm>`_ is sketched out in this section. This is
intended to show an overview of the implementation in order to support the
descriptions above and to illustrate how the plugin interface can be used -
different choices may be made for a production-ready implementation.

The plugin implementation consists of additions to `python/rmm/rmm.py
<https://github.com/rapidsai/rmm/blob/d5831ac5ebb5408ee83f63b7c7d03d8870ecb361/python/rmm/rmm.py>`_:

.. code-block:: python

   # New imports:
   from contextlib import context_manager
   # RMM already has Numba as a dependency, so these imports need not be guarded
   # by a check for the presence of numba.
   from numba.cuda import (HostOnlyCUDAMemoryManager, MemoryPointer, IpcHandle,
                           set_memory_manager)


   # New class implementing the EMM Plugin:
   class RMMNumbaManager(HostOnlyCUDAMemoryManager):
       def memalloc(self, size):
           # Allocates device memory using RMM functions. The finalizer for the
           # allocated memory calls back to RMM to free the memory.
           addr = librmm.rmm_alloc(bytesize, 0)
           ctx = cuda.current_context()
           ptr = ctypes.c_uint64(int(addr))
           finalizer = _make_finalizer(addr, stream)
           return MemoryPointer(ctx, ptr, size, finalizer=finalizer)

      def get_ipc_handle(self, memory):
           """ 
           Get an IPC handle for the memory with offset modified by the RMM memory
           pool.
           """
           # This implementation provides a functional implementation and illustrates
           # what get_ipc_handle needs to do, but it is not a very "clean"
           # implementation, and it relies on borrowing bits of Numba internals to
           # initialise ipchandle. 
           #
           # A more polished implementation might make use of additional functions in
           # the RMM C++ layer for initialising IPC handles, and not use any Numba
           # internals.
           ipchandle = (ctypes.c_byte * 64)()  # IPC handle is 64 bytes
           cuda.cudadrv.memory.driver_funcs.cuIpcGetMemHandle(
               ctypes.byref(ipchandle),
               memory.owner.handle,
           )
           source_info = cuda.current_context().device.get_device_identity()
           ptr = memory.device_ctypes_pointer.value
           offset = librmm.rmm_getallocationoffset(ptr, 0)
           return IpcHandle(memory, ipchandle, memory.size, source_info,
                            offset=offset)

       def get_memory_info(self):
           # Returns a tuple of (free, total) using RMM functionality.
           return get_info() # Function defined in rmm.py

       def initialize(self):
           # Nothing required to initialize RMM here, but this method is added
           # to illustrate that the super() method should also be called.
           super().initialize()

       @contextmanager
       def defer_cleanup(self):
           # Does nothing to defer cleanup - a full implementation may choose to
           # implement a different policy.
           with super().defer_cleanup():
               yield

       @property
       def interface_version(self):
           # As required by the specification
           return 1

   # The existing _make_finalizer function is used by RMMNumbaManager:
   def _make_finalizer(handle, stream):
       """
       Factory to make the finalizer function.
       We need to bind *handle* and *stream* into the actual finalizer, which
       takes no args.
       """

       def finalizer():
           """
           Invoked when the MemoryPointer is freed
           """
           librmm.rmm_free(handle, stream)

       return finalizer

   # Utility function register `RMMNumbaManager` as an EMM:
   def use_rmm_for_numba():
       set_memory_manager(RMMNumbaManager)

   # To support `NUMBA_CUDA_MEMORY_MANAGER=rmm`:
   _numba_memory_manager = RMMNumbaManager

Example usage
^^^^^^^^^^^^^

A simple example that configures Numba to use RMM for memory management and
creates a device array is as follows:

.. code-block:: python

   # example.py
   import rmm 
   import numpy as np

   from numba import cuda

   rmm.use_rmm_for_numba()

   a = np.zeros(10)
   d_a = cuda.to_device(a)
   del(d_a)
   print(rmm.csv_log())

Running this should result in output similar to the following:

.. code-block::

   Event Type,Device ID,Address,Stream,Size (bytes),Free Memory,Total Memory,Current Allocs,Start,End,Elapsed,Location
   Alloc,0,0x7fae06600000,0,80,0,0,1,1.10549,1.1074,0.00191666,<path>/numba/numba/cuda/cudadrv/driver.py:683
   Free,0,0x7fae06600000,0,0,0,0,0,1.10798,1.10921,0.00122238,<path>/numba/numba/utils.py:678

Note that there is some scope for improvement in RMM for detecting the line
number at which the allocation / free occurred, but this is outside the scope of
the example in this proposal.

Setting the memory manager through the environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rather than calling ``rmm.use_rmm_for_numba()`` in the example above, the memory
manager could also be set to use RMM globally with an environment variable, so
the Python interpreter is invoked to run the example as:

.. code-block::

   NUMBA_CUDA_MEMORY_MANAGER="rmm.RMMNumbaManager" python example.py

Numba internal changes
----------------------

This section is intended primarily for Numba developers - those with an interest
in the external interface for implementing EMM plugins may choose to skip over
this section.

Current model / implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At present, memory management is implemented in the
:class:`~numba.cuda.cudadrv.driver.Context` class. It maintains lists of
allocations and deallocations:

* ``allocations`` is a ``numba.core.utils.UniqueDict``, created at context
  creation time.
* ``deallocations`` is an instance of the ``_PendingDeallocs`` class, and is created
  when ``Context.prepare_for_use()`` is called.

These are used to track allocations and deallocations of:

* Device memory
* Pinned memory
* Mapped memory
* Streams
* Events
* Modules

The ``_PendingDeallocs`` class implements the deferred deallocation strategy -
cleanup functions (such as ``cuMemFree``\ ) for the items above are added to its
list of pending deallocations by the finalizers of objects representing
allocations. These finalizers are run when the objects owning them are
garbage-collected by the Python interpreter. When the addition of a new
cleanup function to the deallocation list causes the number or size of pending
deallocations to exceed a configured ratio, the ``_PendingDeallocs`` object runs
deallocators for all items it knows about and then clears its internal pending
list.

See :ref:`deallocation-behavior` for more details of this implementation.

Proposed changes
^^^^^^^^^^^^^^^^

This section outlines the major changes that will be made to support the EMM
plugin interface - there will be various small changes to other parts of Numba
that will be required in order to adapt to these changes; an exhaustive list of
these is not provided.

Context changes
~~~~~~~~~~~~~~~

The ``numba.cuda.cudadrv.driver.Context`` class will no longer directly allocate
and free memory. Instead, the context will hold a reference to a memory manager
instance, and its memory allocation methods will call into the memory manager,
e.g.:

.. code-block:: python

       def memalloc(self, size):
           return self.memory_manager.memalloc(size)

       def memhostalloc(self, size, mapped=False, portable=False, wc=False):
           return self.memory_manager.memhostalloc(size, mapped, portable, wc)

       def mempin(self, owner, pointer, size, mapped=False):
           if mapped and not self.device.CAN_MAP_HOST_MEMORY:
               raise CudaDriverError("%s cannot map host memory" % self.device)
           return self.memory_manager.mempin(owner, pointer, size, mapped)

       def prepare_for_use(self):
           self.memory_manager.initialize()

       def get_memory_info(self):
           self.memory_manager.get_memory_info()

       def get_ipc_handle(self, memory):
           return self.memory_manager.get_ipc_handle(memory)

       def reset(self):
           # ... Already-extant reset logic, plus:
           self._memory_manager.reset()

The ``memory_manager`` member is initialised when the context is created.

The ``memunpin`` method (not shown above but currently exists in the ``Context``
class) has never been implemented - it presently raises a ``NotImplementedError``.
This method arguably un-needed - pinned memory is immediately unpinned by its
finalizer, and unpinning before a finalizer runs would invalidate the state of
``PinnedMemory`` objects for which references are still held. It is proposed that
this is removed when making the other changes to the ``Context`` class.

The ``Context`` class will still instantiate ``self.allocations`` and
``self.deallocations`` as before - these will still be used by the context to
manage the allocations and deallocations of events, streams, and modules, which
are not handled by the EMM plugin.

New components of the ``driver`` module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* ``BaseCUDAMemoryManager``\ : An abstract class, as defined in the plugin interface
  above.
* ``HostOnlyCUDAMemoryManager``\ : A subclass of ``BaseCUDAMemoryManager``\ , with the
  logic from ``Context.memhostalloc`` and ``Context.mempin`` moved into it. This
  class will also create its own ``allocations`` and ``deallocations`` members,
  similarly to how the ``Context`` class creates them. These are used to manage
  the allocations and deallocations of pinned and mapped host memory.
* ``NumbaCUDAMemoryManager``\ : A subclass of ``HostOnlyCUDAMemoryManager``\ , which
  also contains an implementation of ``memalloc`` based on that presently existing
  in the ``Context`` class. This is the default memory manager, and its use
  preserves the behaviour of Numba prior to the addition of the EMM plugin
  interface - that is, all memory allocation and deallocation for Numba arrays
  is handled within Numba.

  * This class shares the ``allocations`` and ``deallocations`` members with its
    parent class ``HostOnlyCUDAMemoryManager``\ , and it uses these for the
    management of device memory that it allocates.

* The ``set_memory_manager`` function, which sets a global pointing to the memory
  manager class. This global initially holds ``NumbaCUDAMemoryManager`` (the
  default).

Staged IPC
~~~~~~~~~~

Staged IPC should not take ownership of the memory that it allocates. When the
default internal memory manager is in use, the memory allocated for the staging
array is already owned. When an EMM plugin is in use, it is not legitimate to
take ownership of the memory.

This change can be made by applying the following small patch, which has been
tested to have no effect on the CUDA test suite:

.. code-block:: diff

   diff --git a/numba/cuda/cudadrv/driver.py b/numba/cuda/cudadrv/driver.py
   index 7832955..f2c1352 100644
   --- a/numba/cuda/cudadrv/driver.py
   +++ b/numba/cuda/cudadrv/driver.py
   @@ -922,7 +922,11 @@ class _StagedIpcImpl(object):
            with cuda.gpus[srcdev.id]:
                impl.close()

   -        return newmem.own()
   +        return newmem

Testing
~~~~~~~

Alongside the addition of appropriate tests for new functionality, there will be
some refactoring of existing tests required, but these changes are not
substantial. Tests of the deallocation strategy (e.g. ``TestDeallocation``\ ,
``TestDeferCleanup``\ ) will need to be modified to ensure that they are
examining the correct set of deallocations. When an EMM plugin is in use, they
will need to be skipped.

Prototyping / experimental implementation
-----------------------------------------

Some prototype / experimental implementations have been produced to guide the
designs presented in this document. The current implementations can be found in:


* Numba branch: https://github.com/gmarkall/numba/tree/grm-numba-nbep-7.
* RMM branch: https://github.com/gmarkall/rmm/tree/grm-numba-nbep-7.
* CuPy implementation:
  https://github.com/gmarkall/nbep-7/blob/master/nbep7/cupy_mempool.py - uses
  an unmodified CuPy.

  * See `CuPy memory management
    docs <https://docs-cupy.chainer.org/en/stable/reference/memory.html>`_.

Current implementation status
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

RMM Plugin
~~~~~~~~~~

For a minimal example, a simple allocation and free using RMM works as expected.
For the example code (similar to the RMM example above):

.. code-block:: python

   import rmm
   import numpy as np

   from numba import cuda

   rmm.use_rmm_for_numba()

   a = np.zeros(10)
   d_a = cuda.to_device(a)
   del(d_a)
   print(rmm.csv_log())

We see the following output:

.. code-block::

   Event Type,Device ID,Address,Stream,Size (bytes),Free Memory,Total Memory,Current Allocs,Start,End,Elapsed,Location
   Alloc,0,0x7f96c7400000,0,80,0,0,1,1.13396,1.13576,0.00180059,<path>/numba/numba/cuda/cudadrv/driver.py:686
   Free,0,0x7f96c7400000,0,0,0,0,0,1.13628,1.13723,0.000956004,<path>/numba/numba/utils.py:678

This output is similar to the expected output from the example usage presented
above (though note that the pointer addresses and timestamps vary compared to
the example), and provides some validation of the example use case.

CuPy Plugin
~~~~~~~~~~~

.. code-block:: python

   from nbep7.cupy_mempool import use_cupy_mm_for_numba
   import numpy as np

   from numba import cuda

   use_cupy_mm_for_numba()

   a = np.zeros(10)
   d_a = cuda.to_device(a)
   del(d_a)

The prototype CuPy plugin has somewhat primitive logging, so we see the output:

.. code-block::

   Allocated 80 bytes at 7f004d400000
   Freeing 80 bytes at 7f004d400000

Numba CUDA Unit tests
^^^^^^^^^^^^^^^^^^^^^

As well as providing correct execution of a simple example, all relevant Numba
CUDA unit tests also pass with the prototype branch, for both the internal memory
manager and the RMM EMM Plugin.

RMM
~~~

The unit test suite can be run with the RMM EMM Plugin with:

.. code-block::

   NUMBA_CUDA_MEMORY_MANAGER=rmm python -m numba.runtests numba.cuda.tests

A summary of the unit test suite output is:

.. code-block::

   Ran 564 tests in 142.211s

   OK (skipped=11)

When running with the built-in Numba memory management, the output is:

.. code-block::

   Ran 564 tests in 133.396s

   OK (skipped=5)

i.e. the changes for using an external memory manager do not break the built-in
Numba memory management. There are an additional 6 skipped tests, from:


* ``TestDeallocation``\ : skipped as it specifically tests Numba's internal
  deallocation strategy.
* ``TestDeferCleanup``\ : skipped as it specifically tests Numba's implementation of
  deferred cleanup.
* ``TestCudaArrayInterface.test_ownership``\ : skipped as Numba does not own memory
  when an EMM Plugin is used, but ownership is assumed by this test case.

CuPy
~~~~

The test suite can be run with the CuPy plugin using:

.. code-block::

   NUMBA_CUDA_MEMORY_MANAGER=nbep7.cupy_mempool python -m numba.runtests numba.cuda.tests

This plugin implementation is presently more primitive than the RMM
implementation, and results in some errors with the unit test suite:

.. code-block::

   Ran 564 tests in 111.699s

   FAILED (errors=8, skipped=11)

The 8 errors are due to a lack of implementation of ``get_ipc_handle`` in the
CuPy EMM Plugin implementation. It is expected that this implementation will be
re-visited and completed so that CuPy can be used stably as an allocator for
Numba in the future.
============================
NBEP 4: Defining C callbacks
============================

:Author: Antoine Pitrou
:Date: April 2016
:Status: Draft


Interfacing with some native libraries (for example written in C
or C++) can necessitate writing native callbacks to provide business logic
to the library.  Some Python-facing libraries may also provide the
alternative of passing a ctypes-wrapped native callback instead of a
Python callback for better performance.  A simple example is the
``scipy.integrate`` package where the user passes the function to be
integrated as a callback.

Users of those libraries may want to benefit from the performance advantage
of running purely native code, while writing their code in Python.
This proposal outlines a scheme to provide such a functionality in
Numba.


Basic usage
===========

We propose adding a new decorator, ``@cfunc``, importable from the main
package.  This decorator allows defining a callback as in the following
example::

   from numba import cfunc
   from numba.types import float64

   # A callback with the C signature `double(double)`

   @cfunc(float64(float64), nopython=True)
   def integrand(x):
       return 1 / x


The ``@cfunc`` decorator returns a "C function" object holding the
resources necessary to run the given compiled function (for example its
LLVM module).  This object has several attributes and methods:

* the ``ctypes`` attribute is a ctypes function object representing
  the native function.

* the ``address`` attribute is the address of the native function code, as
  an integer (note this can also be computed from the ``ctypes`` attribute).

* the ``native_name`` attribute is the symbol under which the function
  can be looked up inside the current process.

* the ``inspect_llvm()`` method returns the IR for the LLVM module
  in which the function is compiled.  It is expected that the ``native_name``
  attribute corresponds to the function's name in the LLVM IR.

The general signature of the decorator is ``cfunc(signature, **options)``.

The ``signature`` must specify the argument types and return type of the
function using Numba types.  In contrary to ``@jit``, the return type cannot
be omitted.

The ``options`` are keyword-only parameters specifying compilation options.
We are expecting that the standard ``@jit`` options (``nopython``,
``forceobj``, ``cache``) can be made to work with ``@cfunc``.


Calling from Numba-compiled functions
-------------------------------------

While the intended use is to pass a callback's address to foreign C
code expecting a function pointer, it should be made possible to call
the C callback from a Numba-compiled function.


Passing array data
==================

Native platform ABIs as used by C or C++ don't have the notion of a shaped
array as in Numpy.  One common solution is to pass a raw data pointer and
one or several size arguments (depending on dimensionality).  Numba must
provide a way to rebuild an array view of this data inside the callback.

::

   from numba import cfunc, carray
   from numba.types import float64, CPointer, void, intp

   # A callback with the C signature `void(double *, double *, size_t)`

   @cfunc(void(CPointer(float64), CPointer(float64), intp))
   def invert(in_ptr, out_ptr, n):
       in_ = carray(in_ptr, (n,))
       out = carray(out_ptr, (n,))
       for i in range(n):
           out[i] = 1 / in_[i]


The ``carray`` function takes ``(pointer, shape, dtype)`` arguments
(``dtype`` being optional) and returns a C-layout array view over the
data *pointer*, with the given *shape* and *dtype*.  *pointer* must
be a ctypes pointer object (not a Python integer).  The array's
dimensionality corresponds to the *shape* tuple's length.  If *dtype*
is not given, the array's dtype corresponds to the *pointer*'s pointee
type.

The ``farray`` function is similar except that it returns a F-layout
array view.


Error handling
==============

There is no standard mechanism in C for error reporting.  Unfortunately,
Numba currently doesn't handle ``try..except`` blocks, which makes it more
difficult for the user to implement the required error reporting scheme.
The current stance of this proposal is to let users guard against invalid
arguments where necessary, and do whatever is required to inform the caller
of the error.

Based on user feedback, we can later add support for some error reporting
schemes, such as returning an integer error code depending on whether an
exception was raised, or setting ``errno``.


Deferred topics
===============

Ahead-of-Time compilation
-------------------------

This proposal doesn't make any provision for AOT compilation of C callbacks.
It would probably necessitate a separate API (a new method on the
``numba.pycc.CC`` object), and the implementation would require exposing
a subset of the C function object's functionality from the compiled C
extension module.

Opaque data pointers
--------------------

Some libraries allow passing an opaque data pointer (``void *``) to a
user-provided callback, to provide any required context for execution
of the callback.  Taking advantage of this functionality would require
adding specific support in Numba, for example the ability to do generic
conversion from ``types.voidptr`` and to take the address of a
Python-facing ``jitclass`` instance.
========================
NBEP 2: Extension points
========================

:Author: Antoine Pitrou
:Date: July 2015
:Status: Draft


Implementing new types or functions in Numba requires hooking into
various mechanisms along the compilation chain (and potentially
outside of it).  This document aims, first, at examining the
current ways of doing so and, second, at making proposals to make
extending easier.

If some of the proposals are implemented, we should first strive
to use and exercise them internally, before exposing the APIs to the
public.

.. note::
   This document doesn't cover CUDA or any other non-CPU backend.


High-level API
==============

There is currently no high-level API, making some use cases more
complicated than they should be.

Proposed changes
----------------

Dedicated module
''''''''''''''''

We propose the addition of a ``numba.extending`` module exposing the main
APIs useful for extending Numba.

Implementing a function
'''''''''''''''''''''''

We propose the addition of a ``@overload`` decorator allowing the
implementation of a given function for use in :term:`nopython mode`.
The overloading function has the same formal signature as the implemented
function, and receives the actual argument types.  It should return a
Python function implementing the overloaded function for the given types.

The following example implements :func:`numpy.where` with
this approach.

.. literalinclude:: np-where-override.py

It is also possible to implement functions already known to Numba, to
support additional types.  The following example implements the
built-in function :func:`len` for tuples with this approach::

   @overload(len)
   def tuple_len(x):
      if isinstance(x, types.BaseTuple):
         # The tuple length is known at compile-time, so simply reify it
         # as a constant.
         n = len(x)
         def len_impl(x):
            return n
         return len_impl


Implementing an attribute
'''''''''''''''''''''''''

We propose the addition of a ``@overload_attribute`` decorator allowing
the implementation of an attribute getter for use in :term:`nopython mode`.

The following example implements the ``.nbytes`` attribute on Numpy arrays::

   @overload_attribute(types.Array, 'nbytes')
   def array_nbytes(arr):
      def get(arr):
          return arr.size * arr.itemsize
      return get

.. note::
   The overload_attribute() signature allows for expansion to also define
   setters and deleters, by letting the decorated function return a
   ``getter, setter, deleter`` tuple instead of a single ``getter``.


Implementing a method
'''''''''''''''''''''

We propose the addition of a ``@overload_method`` decorator allowing the
implementation of an instance method for use in :term:`nopython mode`.

The following example implements the ``.take()`` method on Numpy arrays::

   @overload_method(types.Array, 'take')
   def array_take(arr, indices):
      if isinstance(indices, types.Array):
          def take_impl(arr, indices):
              n = indices.shape[0]
              res = np.empty(n, arr.dtype)
              for i in range(n):
                  res[i] = arr[indices[i]]
              return res
          return take_impl


Exposing a structure member
'''''''''''''''''''''''''''

We propose the addition of a ``make_attribute_wrapper()`` function exposing
an internal field as a visible read-only attribute, for those types backed
by a ``StructModel`` data model.

For example, assuming ``PdIndexType`` is the Numba type of pandas indices,
here is how to expose the underlying Numpy array as a ``._data`` attribute::

   @register_model(PdIndexType)
   class PdIndexModel(models.StructModel):
       def __init__(self, dmm, fe_type):
           members = [
               ('values', fe_type.as_array),
               ]
           models.StructModel.__init__(self, dmm, fe_type, members)

   make_attribute_wrapper(PdIndexType, 'values', '_data')


Typing
======

Numba types
-----------

Numba's standard types are declared in :mod:`numba.types`.  To declare
a new type, one subclasses the base :class:`Type` class or one of its
existing abstract subclasses, and implements the required functionality.

Proposed changes
''''''''''''''''

No change required.


Type inference on values
------------------------

Values of a new type need to be type-inferred if they can appear as
function arguments or constants.  The core machinery is in
:mod:`numba.typing.typeof`.

In the common case where some Python class or classes map exclusively
to the new type, one can extend a generic function to dispatch on said
classes, e.g.::

   from numba.typing.typeof import typeof_impl

   @typeof_impl(MyClass)
   def _typeof_myclass(val, c):
      if "some condition":
         return MyType(...)

The ``typeof_impl`` specialization must return a Numba type instance,
or None if the value failed typing.

(when one controls the class being type-inferred, an alternative
to ``typeof_impl`` is to define a ``_numba_type_`` property on the class)

In the rarer case where the new type can denote various Python classes
that are impossible to enumerate, one must insert a manual check in the
fallback implementation of the ``typeof_impl`` generic function.

Proposed changes
''''''''''''''''

Allow people to define a generic hook without monkeypatching the
fallback implementation.


Fast path for type inference on function arguments
--------------------------------------------------

Optionally, one may want to allow a new type to participate in the
fast type resolution (written in C code) to minimize function call
overhead when a JIT-compiled function is called with the new type.
One must then insert the required checks and implementation in
the ``_typeof.c`` file, presumably inside the ``compute_fingerprint()``
function.

Proposed changes
''''''''''''''''

None.  Adding generic hooks to C code embedded in a C Python extension
is too delicate a change.


Type inference on operations
----------------------------

Values resulting from various operations (function calls, operators, etc.)
are typed using a set of helpers called "templates".  One can define a
new template by subclass one of the existing base classes and implement
the desired inference mechanism.  The template is explicitly registered
with the type inference machinery using a decorator.

The :class:`ConcreteTemplate` base class allows one to define inference as
a set of supported signatures for a given operation.  The following example
types the modulo operator::

   @builtin
   class BinOpMod(ConcreteTemplate):
       key = "%"
       cases = [signature(op, op, op)
                for op in sorted(types.signed_domain)]
       cases += [signature(op, op, op)
                 for op in sorted(types.unsigned_domain)]
       cases += [signature(op, op, op) for op in sorted(types.real_domain)]

(note that type *instances* are used in the signatures, severely
limiting the amount of genericity that can be expressed)

The :class:`AbstractTemplate` base class allows to define inference
programmatically, giving it full flexibility.  Here is a simplistic
example of how tuple indexing (i.e. the ``__getitem__`` operator) can
be expressed::

   @builtin
   class GetItemUniTuple(AbstractTemplate):
       key = "getitem"

       def generic(self, args, kws):
           tup, idx = args
           if isinstance(tup, types.UniTuple) and isinstance(idx, types.Integer):
               return signature(tup.dtype, tup, idx)


The :class:`AttributeTemplate` base class allows to type the attributes
and methods of a given type.  Here is an example, typing the ``.real``
and ``.imag`` attributes of complex numbers::

   @builtin_attr
   class ComplexAttribute(AttributeTemplate):
       key = types.Complex

       def resolve_real(self, ty):
           return ty.underlying_float

       def resolve_imag(self, ty):
           return ty.underlying_float

.. note::
   :class:`AttributeTemplate` only works for getting attributes.  Setting
   an attribute's value is hardcoded in :mod:`numba.typeinfer`.

The :class:`CallableTemplate` base class offers an easier way to parse
flexible function signatures, by letting one define a callable that has
the same definition as the function being typed.  For example, here is how
one could hypothetically type Python's ``sorted`` function if Numba supported
lists::

   @builtin
   class Sorted(CallableTemplate):
       key = sorted

       def generic(self):
           def typer(iterable, key=None, reverse=None):
               if reverse is not None and not isinstance(reverse, types.Boolean):
                   return
               if key is not None and not isinstance(key, types.Callable):
                   return
               if not isinstance(iterable, types.Iterable):
                   return
               return types.List(iterable.iterator_type.yield_type)

           return typer

(note you can return just the function's return type instead of the
full signature)

Proposed changes
''''''''''''''''

Naming of the various decorators is quite vague and confusing.  We propose
renaming ``@builtin`` to ``@infer``, ``@builtin_attr`` to ``@infer_getattr``
and ``builtin_global`` to ``infer_global``.

The two-step declaration for global values is a bit verbose, we propose
simplifying it by allowing the use of ``infer_global`` as a decorator::

   @infer_global(len)
   class Len(AbstractTemplate):
       key = len

       def generic(self, args, kws):
           assert not kws
           (val,) = args
           if isinstance(val, (types.Buffer, types.BaseTuple)):
               return signature(types.intp, val)

The class-based API can feel clumsy, we can add a functional API for
some of the template kinds:

.. code-block:: python

   @type_callable(sorted)
   def type_sorted(context):
       def typer(iterable, key=None, reverse=None):
           # [same function as above]

       return typer


Code generation
===============

Concrete representation of values of a Numba type
-------------------------------------------------

Any concrete Numba type must be able to be represented in LLVM form
(for variable storage, argument passing, etc.).  One defines that
representation by implementing a datamodel class and registering it
with a decorator.  Datamodel classes for standard types are defined
in :mod:`numba.datamodel.models`.

Proposed changes
''''''''''''''''

No change required.

Conversion between types
------------------------

Implicit conversion between Numba types is currently implemented as a
monolithic sequence of choices and type checks in the
:meth:`BaseContext.cast` method.  To add a new implicit conversion, one
appends a type-specific check in that method.

Boolean evaluation is a special case of implicit conversion (the
destination type being :class:`types.Boolean`).

.. note::
   Explicit conversion is seen as a regular operation, e.g. a constructor
   call.

Proposed changes
''''''''''''''''

Add a generic function for implicit conversion, with multiple dispatch
based on the source and destination types.  Here is an example showing
how to write a float-to-integer conversion::

   @lower_cast(types.Float, types.Integer)
   def float_to_integer(context, builder, fromty, toty, val):
       lty = context.get_value_type(toty)
       if toty.signed:
           return builder.fptosi(val, lty)
       else:
           return builder.fptoui(val, lty)


Implementation of an operation
------------------------------

Other operations are implemented and registered using a set of generic
functions and decorators.  For example, here is how lookup for a the ``.ndim``
attribute on Numpy arrays is implemented::

   @builtin_attr
   @impl_attribute(types.Kind(types.Array), "ndim", types.intp)
   def array_ndim(context, builder, typ, value):
       return context.get_constant(types.intp, typ.ndim)

And here is how calling ``len()`` on a tuple value is implemented::

   @builtin
   @implement(types.len_type, types.Kind(types.BaseTuple))
   def tuple_len(context, builder, sig, args):
       tupty, = sig.args
       retty = sig.return_type
       return context.get_constant(retty, len(tupty.types))

Proposed changes
''''''''''''''''

Review and streamine the API.  Drop the requirement to write
``types.Kind(...)`` explicitly.  Remove the separate ``@implement``
decorator and rename ``@builtin`` to ``@lower_builtin``, ``@builtin_attr``
to ``@lower_getattr``, etc.

Add decorators to implement ``setattr()`` operations, named
``@lower_setattr`` and ``@lower_setattr_generic``.


Conversion from / to Python objects
-----------------------------------

Some types need to be converted from or to Python objects, if they can
be passed as function arguments or returned from a function.  The
corresponding boxing and unboxing operations are implemented using
a generic function.  The implementations for standard Numba types
are in :mod:`numba.targets.boxing`.  For example, here is the boxing
implementation for a boolean value::

   @box(types.Boolean)
   def box_bool(c, typ, val):
       longval = c.builder.zext(val, c.pyapi.long)
       return c.pyapi.bool_from_long(longval)

Proposed changes
''''''''''''''''

Change the implementation signature from ``(c, typ, val)`` to
``(typ, val, c)``, to match the one chosen for the ``typeof_impl``
generic function.
======================
NBEP 5: Type Inference
======================

:Author: Siu Kwan Lam
:Date: Sept 2016
:Status: Draft


This document describes the current type inference implementation in numba.


Introduction
============

Numba uses type information to ensure that every variable in the user code can
be correctly lowered (translated into a low-level representation).  The type of
a variable describes the set of valid operations and available attributes.
Resolving this information during compilation avoids the overhead of type
checking and dispatching at runtime.  However, Python is dynamically typed and
the user does not declare variable types.  Since type information is absent,
we use type inference to reconstruct the missing information.


Numba Type Semantic
===================

Type inference operates on :term:`Numba IR`, a mostly static-single-assignment (SSA)
encoding of the Python bytecode.  Conceptually, all intermediate values in the
Python code are explicitly assigned to a variable in the IR.  Numba enforces
that each IR variable to have one type only.  A user variable (from the Python
source code) can be mapped to multiple variables in the IR.  They are *versions*
of a variable.  Each time a user variable is assigned to, a new version is
created.  From that point, all subsequent references will use the new version.
The user variable *evolves* as the function logic updates its type.  Merge
points (e.g. subsequent block to an if-else, the loop body, etc..) in the control
flow need extra care. At each merge point, a new version is implicitly created
to merge the different variable versions from the incoming paths.
The merging of the variable versions may translate into an implicit cast.

Numba uses function overloading to emulate Python duck-typing.  The type of a
function can contain multiple call signatures that accept different argument
types and yield different return types.  The process to decide the best
signature for an overloaded function is called *overload resolution*.
Numba partially implements the C++ overload resolution scheme
(`ISOCPP`_ 13.3 Overload Resolution).  The scheme uses a "best fit" algorithm by
ranking each argument symmetrically.  The five possible rankings in increasing
order of penalty are:

* *Exact*: the expected type is the same as the actual type.
* *Promotion*: the actual type can be upcast to the expected type by extending
  the precision without changing the behavior.
* *Safe conversion*: the actual type can be cast to the expected type by changing
  the type without losing information.
* *Unsafe conversion*: the actual type can be cast to the expected type by
  changing the type or downcasting the type even if it is imprecise.
* *No match*: no valid operation can convert the actual type to the expected type.

It is possible to have an ambiguous resolution.  For example, a function with
signatures ``(int16, int32)`` and ``(int32, int16)`` can become ambiguous if
presented with the argument types ``(int32, int32)``, because demoting either
argument to ``int16`` is equally "fit".  Fortunately, numba can usually resolve
such ambiguity by compiling a new version with the exact signature
``(int32, int32)``.  When compilation is disabled and there are multiple
signatures with equal fit, an exception is raised.

Type Inference
==============

The type inference in numba has three important components---type
variable, constraint network, and typing context.

* The *typing context* provides all the type information and typing related
  operations, including the logic for type unification, and the logic for typing
  of global and constant values.  It defines the semantic of the language that
  can be compiled by numba.

* A *type variable* holds the type of each variable (in the Numba IR).
  Conceptually, it is initialized to the universal type and, as it is re-assigned,
  it stores a common type by unifying the new type with the existing type.  The
  common type must be able to represent values of the new type and the existing
  type.  Type conversion is applied as necessary and precision loss is
  accepted for usability reason.

* The *constraint network* is a dependency graph built from the IR.  Each
  node represents an operation in the Numba IR and updates at least one type
  variable.  There may be cycles due to loops in user code.

The type inference process starts by seeding the argument types.  These initial
types are propagated in the constraint network, which eventually fills all the
type variables.  Due to cycles in the network, the process repeats until all
type variables converge or it fails with undecidable types.

Type unification always returns a more "general" (quoted because unsafe conversion
is allowed) type.  Types will converge to the least "general" type that
can represent all possible values that the variable can hold.  Since unification
will never move down the type hierarchy and there is a single top type, the
universal type---``object``, the type inference is guaranteed to converge.

A failure in type inference can be caused by two reasons.  The first reason is user
error due to incorrect use of a type.  This type of error will also trigger an
exception in regular python execution.  The second reason is due to the use of an
unsupported feature, but the code is otherwise valid in regular python
execution.  Upon an error, the type inference will set all types to the object
type.  As a result, numba will fallback to *object-mode*.

Since functions can be overloaded, the type inference needs to decide the
type signature used at each call site.  The overload resolution is applied to
all known overload versions of the callee function described in *call-templates*.
A call-template can either be concrete or abstract.  A concrete call-template
defines a fixed list of all possible signatures.  An abstract call-template
defines the logic to compute the accepted signature and it is used to implement
generic functions.

Numba-compiled functions are generic functions due to their ability to compile
new versions.  When it sees a new set of argument types, it triggers type
inference to validate and determine the return type. When there are nested calls
for numba-compiled functions, each call-site triggers type inference.
This poses a problem to recursive functions because the type inference will also
be triggered recursively.  Currently, simple single recursion is supported if
the signature is user-annotated by the user, which avoids unbound recursion in
type inference that will never terminate.

.. _ISOCPP: http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2014/n4296.pdf===========================
Numba Enhancement Proposals
===========================

Numba Enhancement Proposals (not really abbreviated "NEPs", since "NEP"
is already taken by the Numpy project) describe proposed changes to Numba.
They are modeled on Python Enhancement Proposals (PEPs) and Numpy Enhancement
Proposals, and are typically written up when important changes
(behavioural changes, feature additions...) to Numba are proposed.

This page provides an overview of all proposals, making only a distinction
between the ones that have been implemented and those that have not been
implemented.

Implemented proposals
---------------------

.. toctree::
   :maxdepth: 1

   integer-typing.rst
   external-memory-management.rst

Other proposals
---------------

.. toctree::
   :maxdepth: 1

   extension-points.rst
   jit-classes.rst
   cfunc.rst
   type-inference.rst
   typing_recursion.rst

.. _nbep-1:

=================================
NBEP 1: Changes in integer typing
=================================

:Author: Antoine Pitrou
:Date: July 2015
:Status: Final


Current semantics
=================

Type inference of integers in Numba currently has some subtleties
and some corner cases.  The simple case is when some variable has an obvious
Numba type (for example because it is the result of a constructor call to a
Numpy scalar type such as ``np.int64``). That case suffers no ambiguity.

The less simple case is when a variable doesn't bear such explicit
information.  This can happen because it is inferred from a built-in Python
``int`` value, or from an arithmetic operation between two integers, or
other cases yet.  Then Numba has a number of rules to infer the resulting
Numba type, especially its signedness and bitwidth.

Currently, the generic case could be summarized as: *start small,
grow bigger as required*.  Concretely:

1. Each constant or pseudo-constant is inferred using the *smallest signed
   integer type* that can correctly represent it (or, possibly, ``uint64``
   for positive integers between ``2**63`` and ``2**64 - 1``).
2. The result of an operation is typed so as to ensure safe representation
   in the face of overflow and other magnitude increases (for example,
   ``int32 + int32`` would be typed ``int64``).
3. As an exception, a Python ``int`` used as function argument is always
   typed ``intp``, a pointer-size integer.  This is to avoid the proliferation
   of compiled specializations, as otherwise various integer bitwidths
   in input arguments may produce multiple signatures.

.. note::
   The second rule above (the "respect magnitude increases" rule)
   reproduces Numpy's behaviour with arithmetic on scalar values.
   Numba, however, has different implementation and performance constraints
   than Numpy scalars.

   It is worth nothing, by the way, that Numpy arrays do not implement
   said rule (i.e. ``array(int32) + array(int32)`` is typed ``array(int32)``,
   not ``array(int64)``).  Probably because this makes performance more
   controllable.

This has several non-obvious side-effects:

1. It is difficult to predict the precise type of a value inside a function,
   after several operations.  The basic operands in an expression tree
   may for example be ``int8`` but the end result may be ``int64``.  Whether
   this is desirable or not is an open question; it is good for correctness,
   but potentially bad for performance.

2. In trying to follow the correctness over predictability rule, some values
   can actually leave the integer realm.  For example, ``int64 + uint64``
   is typed ``float64`` in order to avoid magnitude losses (but incidentally
   will lose precision on large integer values...), again following Numpy's
   semantics for scalars.  This is usually not intended by the user.

3. More complicated scenarios can produce unexpected errors at the type unification
   stage.  An example is at `Github issue 1299 <https://github.com/numba/numba/issues/1299>`_,
   the gist of which is reproduced here::

      @jit(nopython=True)
      def f():
          variable = 0
          for i in range(1):
              variable = variable + 1
          return np.arange(variable)

   At the time of this writing, this fails compiling, on a 64-bit system,
   with the error::

      numba.errors.TypingError: Failed at nopython (nopython frontend)
      Can't unify types of variable '$48.4': $48.4 := {array(int32, 1d, C), array(int64, 1d, C)}

   People expert with Numba's type unification system can understand why.
   But the user is caught in mystery.


Proposal: predictable width-conserving typing
=============================================

We propose to turn the current typing philosophy on its head.  Instead
of "*start small and grow as required*", we propose "*start big and keep
the width unchanged*".

Concretely:

1. The typing of Python ``int`` values used as function arguments doesn't
   change, as it works satisfyingly and doesn't surprise the user.

2. The typing of integer *constants* (and pseudo-constants) changes to match
   the typing of integer arguments.  That is, every non-explicitly typed
   integer constant is typed ``intp``, the pointer-sized integer; except for
   the rare cases where ``int64`` (on 32-bit systems) or ``uint64`` is
   required.

3. Operations on integers promote bitwidth to ``intp``, if smaller, otherwise
   they don't promote.  For example, on a 32-bit machine, ``int8 + int8``
   is typed ``int32``, as is ``int32 + int32``.  However, ``int64 + int64``
   is typed ``int64``.

4. Furthermore, mixed operations between signed and unsigned fall back to
   signed, while following the same bitwidth rule.  For example, on a
   32-bit machine, ``int8 + uint16`` is typed ``int32``, as is
   ``uint32 + int32``.


Proposal impact
===============

Semantics
---------

With this proposal, the semantics become clearer.  Regardless of whether
the arguments and constants of a function were explicitly typed or not,
the results of various expressions at any point in the function have
easily predictable types.

When using built-in Python ``int``, the user gets acceptable magnitude
(32 or 64 bits depending on the system's bitness), and the type remains
the same across all computations.

When explicitly using smaller bitwidths, intermediate results don't
suffer from magnitude loss, since their bitwidth is promoted to ``intp``.

There is also less potential for annoyances with the type unification
system as demonstrated above.  The user would have to force several
different types to be faced with such an error.

One potential cause for concern is the discrepancy with Numpy's scalar
semantics; but at the same time this brings Numba scalar semantics closer
to array semantics (both Numba's and Numpy's), which seems a desirable
outcome as well.

It is worth pointing out that some sources of integer numbers, such
as the ``range()`` built-in, always yield 32-bit integers or larger.
This proposal could be an opportunity to standardize them on ``intp``.

Performance
-----------

Except in trivial cases, it seems unlikely that the current "best fit"
behaviour for integer constants really brings a performance benefit.  After
all, most integers in Numba code would either be stored in arrays (with
well-known types, chosen by the user) or be used as indices, where a ``int8``
is highly unlikely to fare better than a ``intp`` (actually, it may be worse,
if LLVM isn't able to optimize away the required sign-extension).

As a side note, the default use of ``intp`` rather than ``int64``
ensures that 32-bit systems won't suffer from poor arithmetic performance.

Implementation
--------------

Optimistically, this proposal may simplify some Numba internals a bit.
Or, at least, it doesn't threaten to make them significantly more complicated.

Limitations
-----------

This proposal doesn't really solve the combination of signed and unsigned
integers.  It is geared mostly at solving the bitwidth issues, which are
a somewhat common cause of pain for users.  Unsigned integers are in
practice very uncommon in Numba-compiled code, except when explicitly
asked for, and therefore much less of a pain point.

On the bitwidth front, 32-bit systems could still show discrepancies based
on the values of constants: if a constant is too large to fit in 32 bits,
it is typed ``int64``, which propagates through other computations.
This would be a reminiscence of the current behaviour, but rarer and much
more controlled still.

Long-term horizon
-----------------

While we believe this proposal makes Numba's behaviour more regular and more
predictable, it also pulls it further from general compatibility with pure
Python semantics, where users can assume arbitrary-precision integers without
any truncation issues.

Contributing to Numba
=====================

We welcome people who want to make contributions to Numba, big or small!
Even simple documentation improvements are encouraged.  If you have
questions, don't hesitate to ask them (see below).


Communication
-------------

Real-time Chat
''''''''''''''

Numba uses Gitter for public real-time chat.  To help improve the
signal-to-noise ratio, we have two channels:

* `numba/numba <https://gitter.im/numba/numba>`_: General Numba discussion,
  questions, and debugging help.
* `numba/numba-dev <https://gitter.im/numba/numba-dev>`_: Discussion of PRs,
  planning, release coordination, etc.

Both channels are public, but we may ask that discussions on numba-dev move to
the numba channel.  This is simply to ensure that numba-dev is easy for core
developers to keep up with.

Note that the Github issue tracker is the best place to report bugs.  Bug
reports in chat are difficult to track and likely to be lost.

Forum
.....

Numba uses Discourse as a forum for longer running threads such as design
discussions and roadmap planning. There are various categories available and it
can be reached at: `numba.discourse.group <https://numba.discourse.group/>`_.

Weekly Meetings
'''''''''''''''

The core Numba developers have a weekly video conference to discuss roadmap,
feature planning, and outstanding issues.  These meetings are entirely public,
details are posted on
`numba.discourse.group Announcements <https://numba.discourse.group/c/announcements/>`_
and everyone is welcome to join the discussion. Minutes will be taken and will
be posted to the
`Numba wiki <https://github.com/numba/numba/wiki/Meeting-Minutes>`_.

.. _report-numba-bugs:

Bug tracker
''''''''''''

We use the `Github issue tracker <https://github.com/numba/numba/issues>`_
to track both bug reports and feature requests.  If you report an issue,
please include specifics:

* what you are trying to do;
* which operating system you have and which version of Numba you are running;
* how Numba is misbehaving, e.g. the full error traceback, or the unexpected
  results you are getting;
* as far as possible, a code snippet that allows full reproduction of your
  problem.

Getting set up
--------------

If you want to contribute, we recommend you fork our `Github repository
<https://github.com/numba/numba>`_, then create a branch representing
your work.  When your work is ready, you should submit it as a pull
request from the Github interface.

If you want, you can submit a pull request even when you haven't finished
working.  This can be useful to gather feedback, or to stress your changes
against the :ref:`continuous integration <continuous_integration_testing>`
platform.  In this case, please prepend ``[WIP]`` to your pull request's title.

.. _buildenv:

Build environment
'''''''''''''''''

Numba has a number of dependencies (mostly `NumPy <http://www.numpy.org/>`_
and `llvmlite <https://github.com/numba/llvmlite>`_) with non-trivial build
instructions.  Unless you want to build those dependencies yourself, we
recommend you use `conda <http://conda.pydata.org/miniconda.html>`_ to
create a dedicated development environment and install precompiled versions
of those dependencies there.

First add the Anaconda Cloud ``numba`` channel so as to get development builds
of the llvmlite library::

   $ conda config --add channels numba

Then create an environment with the right dependencies::

   $ conda create -n numbaenv python=3.8 llvmlite numpy scipy jinja2 cffi

.. note::
   This installs an environment based on Python 3.8, but you can of course
   choose another version supported by Numba.  To test additional features,
   you may also need to install ``tbb`` and/or ``llvm-openmp`` and
   ``intel-openmp``.

To activate the environment for the current shell session::

   $ conda activate numbaenv

.. note::
   These instructions are for a standard Linux shell.  You may need to
   adapt them for other platforms.

Once the environment is activated, you have a dedicated Python with the
required dependencies::

    $ python
    Python 3.8.5 (default, Sep  4 2020, 07:30:14) 
    [GCC 7.3.0] :: Anaconda, Inc. on linux
    Type "help", "copyright", "credits" or "license" for more information.

    >>> import llvmlite
    >>> llvmlite.__version__
    '0.35.0'


Building Numba
''''''''''''''

For a convenient development workflow, we recommend you build Numba inside
its source checkout::

   $ git clone git://github.com/numba/numba.git
   $ cd numba
   $ python setup.py build_ext --inplace

This assumes you have a working C compiler and runtime on your development
system.  You will have to run this command again whenever you modify
C files inside the Numba source tree.

The ``build_ext`` command in Numba's setup also accepts the following
arguments:

- ``--noopt``: This disables optimization when compiling Numba's CPython
  extensions, which makes debugging them much easier. Recommended in
  conjunction with the standard ``build_ext`` option ``--debug``.
- ``--werror``: Compiles Numba's CPython extensions with the ``-Werror`` flag.
- ``--wall``: Compiles Numba's CPython extensions with the ``-Wall`` flag.

Note that Numba's CI and the conda recipe for Linux build with the ``--werror``
and ``--wall`` flags, so any contributions that change the CPython extensions
should be tested with these flags too.

Running tests
'''''''''''''

Numba is validated using a test suite comprised of various kind of tests
(unit tests, functional tests). The test suite is written using the
standard :py:mod:`unittest` framework.

The tests can be executed via ``python -m numba.runtests``.  If you are
running Numba from a source checkout, you can type ``./runtests.py``
as a shortcut.  Various options are supported to influence test running
and reporting.  Pass ``-h`` or ``--help`` to get a glimpse at those options.
Examples:

* to list all available tests::

    $ python -m numba.runtests -l

* to list tests from a specific (sub-)suite::

    $ python -m numba.runtests -l numba.tests.test_usecases

* to run those tests::

    $ python -m numba.runtests numba.tests.test_usecases

* to run all tests in parallel, using multiple sub-processes::

    $ python -m numba.runtests -m

* For a detailed list of all options::

    $ python -m numba.runtests -h

The numba test suite can take a long time to complete.  When you want to avoid
the long wait,  it is useful to focus on the failing tests first with the
following test runner options:

* The ``--failed-first`` option is added to capture the list of failed tests
  and to re-execute them first::

    $ python -m numba.runtests --failed-first -m -v -b

* The ``--last-failed`` option is used with ``--failed-first`` to execute
  the previously failed tests only::

    $ python -m numba.runtests --last-failed -m -v -b

When debugging, it is useful to turn on logging.  Numba logs using the
standard ``logging`` module.  One can use the standard ways (i.e.
``logging.basicConfig``) to configure the logging behavior.  To enable logging
in the test runner, there is a ``--log`` flag for convenience::

    $ python -m numba.runtests --log

To enable :ref:`runtime type-checking <type_anno_check>`, set the environment
variable ``NUMBA_USE_TYPEGUARD=1`` and use `runtests.py` from the source root
instead. For example::

    $ NUMBA_USE_TYPEGUARD=1 python runtests.py


Development rules
-----------------

Code reviews
''''''''''''

Any non-trivial change should go through a code review by one or several of
the core developers.  The recommended process is to submit a pull request
on github.

A code review should try to assess the following criteria:

* general design and correctness
* code structure and maintainability
* coding conventions
* docstrings, comments
* test coverage

Coding conventions
''''''''''''''''''

All Python code should follow :pep:`8`.  Our C code doesn't have a
well-defined coding style (would it be nice to follow :pep:`7`?).
Code and documentation should generally fit within 80 columns, for
maximum readability with all existing tools (such as code review UIs).

Numba uses `Flake8 <http://flake8.pycqa.org/en/latest/>`_ to ensure a consistent
Python code format throughout the project. ``flake8`` can be installed
with ``pip`` or ``conda`` and then run from the root of the Numba repository::

    flake8 numba

Optionally, you may wish to setup `pre-commit hooks <https://pre-commit.com/>`_
to automatically run ``flake8`` when you make a git commit. This can be
done by installing ``pre-commit``::

    pip install pre-commit

and then running::

    pre-commit install

from the root of the Numba repository. Now ``flake8`` will be run each time
you commit changes. You can skip this check with ``git commit --no-verify``.

Numba has started the process of using `type hints <https://www.python.org/dev/peps/pep-0484/>`_ in its code base. This
will be a gradual process of extending the number of files that use type hints, as well as going from voluntary to
mandatory type hints for new features. `Mypy <http://mypy-lang.org/>`_ is used for automated static checking.

At the moment, only certain files are checked by mypy. The list can be found in ``mypy.ini``. When making changes to
those files, it is necessary to add the required type hints such that mypy tests will pass. Only in exceptional
circumstances should ``type: ignore`` comments be used.

If you are contributing a new feature, we encourage you to use type hints, even if the file is not currently in the
checklist. If you want to contribute type hints to enable a new file to be in the checklist, please add the file to the
``files`` variable in ``mypy.ini``, and decide what level of compliance you are targetting. Level 3 is basic static
checks, while levels 2 and 1 represent stricter checking. The levels are described in details in ``mypy.ini``.

There is potential for confusion between the Numba module ``typing`` and Python built-in module ``typing`` used for type
hints, as well as between Numba types---such as ``Dict`` or ``Literal``---and ``typing`` types of the same name.
To mitigate the risk of confusion we use a naming convention by which objects of the built-in ``typing`` module are
imported with an ``pt`` prefix. For example, ``typing.Dict`` is imported as ``from typing import Dict as ptDict``.

Stability
'''''''''

The repository's ``master`` branch is expected to be stable at all times.
This translates into the fact that the test suite passes without errors
on all supported platforms (see below).  This also means that a pull request
also needs to pass the test suite before it is merged in.

.. _platform_support:

Platform support
''''''''''''''''

Every commit to the master branch is automatically tested on all of the
platforms Numba supports. This includes ARMv8, POWER8, and NVIDIA GPUs.
The build system however is internal to Anaconda, so we also use
`Azure <https://dev.azure.com/numba/numba/_build>`_ to provide public continuous
integration information for as many combinations as can be supported by the
service.  Azure CI automatically tests all pull requests on Windows, OS X and
Linux, as well as a sampling of different Python and NumPy versions. If you see
problems on platforms you are unfamiliar with, feel free to ask for help in your
pull request. The Numba core developers can help diagnose cross-platform
compatibility issues. Also see the :ref:`continuous integration
<continuous_integration_testing>` section on how public CI is implemented.

.. _continuous_integration_testing:

Continuous integration testing
''''''''''''''''''''''''''''''

The Numba test suite causes CI systems a lot of grief:

#. It's huge, 9000+ tests.
#. In part because of 1. and that compilers are pretty involved, the test suite
   takes a long time to run.
#. There's sections of the test suite that are deliberately designed to stress
   systems almost to the point of failure (tests which concurrently compile and
   execute with threads and fork processes etc).
#. The combination of things that Numba has to test well exceeds the capacity of
   any public CI system, (Python versions x NumPy versions x Operating systems
   x Architectures x feature libraries (e.g. SVML) x threading backends
   (e.g. OpenMP, TBB)) and then there's CUDA too and all its version
   variants.

As a result of the above, public CI is implemented as follows:

#. The combination of OS x Python x NumPy x Various Features in the testing
   matrix is designed to give a good indicative result for whether "this pull
   request is probably ok".
#. When public CI runs it:

   #. Looks for files that contain tests that have been altered by the proposed
      change and runs these on the whole testing matrix.
   #. Runs a subset of the test suite on each part of the testing matrix. i.e.
      slice the test suite up by the number of combinations in the testing
      matrix and each combination runs one chunk. This is done for speed,
      because public CI cannot cope with the load else.

If a Pull Request (PR) changes CUDA code or will affect the CUDA target, it
needs to be run on `gpuCI <https://gpuci.gpuopenanalytics.com/job/numba/>`_.
This can be triggered by one of the Numba maintainers commenting ``run gpuCI
tests`` on the PR discussion. This runs the CUDA testsuite with various CUDA
toolkit versions on Linux, to provide some initial confidence in the
correctness of the changes with respect to CUDA. Following approval, the PR
will also be run on Numba's build farm to test other configurations with CUDA
(including Windows, which is not tested by gpuCI).

If the PR is not CUDA-related but makes changes to something that the core
developers consider risky, then it will also be run on the Numba farm just to
make sure. The Numba project's private build and test farm will actually
exercise all the applicable tests on all the combinations noted above on real
hardware!


.. _type_anno_check:

Type annotation and runtime type checking
'''''''''''''''''''''''''''''''''''''''''

Numba is slowly gaining type annotations. To facilitate the review of pull
requests that are incrementally adding type annotations, the test suite uses
`typeguard`_ to perform runtime type checking. This helps verify the validity
of type annotations.

To enable runtime type checking in the test suite, users can use
`runtests.py`_ in the source root as the test runner and set environment
variable ``NUMBA_USE_TYPEGUARD=1``. For example::

    $ NUMBA_USE_TYPEGUARD=1 python runtests.py numba.tests

Things that help with pull requests
'''''''''''''''''''''''''''''''''''

Even with the mitigating design above public CI can get overloaded which causes
a backlog of builds. It's therefore really helpful when opening pull requests if
you can limit the frequency of pushing changes. Ideally, please squash commits
to reduce the number of patches and/or push as infrequently as possible. Also,
once a pull request review has started, please don't rebase/force push/squash
or do anything that rewrites history of the reviewed code as GitHub cannot track
this and it makes it very hard for reviewers to see what has changed.

The core developers thank everyone for their cooperation with the above!

Why is my pull request/issue seemingly being ignored?
'''''''''''''''''''''''''''''''''''''''''''''''''''''

Numba is an open source project and like many similar projects it has limited
resources. As a result, it is unfortunately necessary for the core developers to
associate a priority with issues/pull requests (PR). A great way to move your
issue/PR up the priority queue is to help out somewhere else in the project so
as to free up core developer time. Examples of ways to help:

* Perform an initial review on a PR. This often doesn't require compiler
  engineering knowledge and just involves checking that the proposed patch is of
  good quality, fixes the problem/implements the feature, is well tested and
  documented.
* Debug an issue, there are numerous issues which `"need triage" <https://github.com/numba/numba/issues?q=is%3Aissue+is%3Aopen+label%3Aneedtriage>`_
  which essentially involves debugging the reported problem. Even if you cannot
  get right to the bottom of a problem, leaving notes about what was discovered
  for someone else is also helpful.
* Answer questions/provide help for users on `discourse <https://numba.discourse.group/>`_
  and/or `gitter.im <https://gitter.im/numba/numba>`_.

The core developers thank everyone for their understanding with the above!

Documentation
-------------

The Numba documentation is split over two repositories:

* This documentation is in the ``docs`` directory inside the
  `Numba repository <https://github.com/numba/numba>`_.

* The `Numba homepage <https://numba.pydata.org>`_ has its sources in a
  separate repository at https://github.com/numba/numba-webpage


Main documentation
''''''''''''''''''

This documentation is under the ``docs`` directory of the `Numba repository`_.
It is built with `Sphinx <http://sphinx-doc.org/>`_ and
`numpydoc <https://numpydoc.readthedocs.io/>`_, which are available using
conda or pip; i.e. ``conda install sphinx numpydoc``.

To build the documentation, you need the bootstrap theme::

   $ pip install sphinx_bootstrap_theme

You can edit the source files under ``docs/source/``, after which you can
build and check the documentation::

   $ make html
   $ open _build/html/index.html

Core developers can upload this documentation to the Numba website
at https://numba.pydata.org by using the ``gh-pages.py`` script under ``docs``::

   $ python gh-pages.py version  # version can be 'dev' or '0.16' etc

then verify the repository under the ``gh-pages`` directory and use
``git push``.

Web site homepage
'''''''''''''''''

The Numba homepage on https://numba.pydata.org can be fetched from here:
https://github.com/numba/numba-webpage

After pushing documentation to a new version, core developers will want to
update the website.  Some notable files:

* ``index.rst``       # Update main page
* ``_templates/sidebar_versions.html``    # Update sidebar links
* ``doc.rst``         # Update after adding a new version for numba docs
* ``download.rst``    # Updata after uploading new numba version to pypi

After updating run::

   $ make html

and check out ``_build/html/index.html``.  To push updates to the Web site::

   $ python _scripts/gh-pages.py

then verify the repository under the ``gh-pages`` directory.  Make sure the
``CNAME`` file is present and contains a single line for ``numba.pydata.org``.
Finally, use ``git push`` to update the website.


.. _typeguard: https://typeguard.readthedocs.io/en/latest/
.. _runtests.py: https://github.com/numba/numba/blob/master/runtests.py

.. _arch-generators:

===================
Notes on generators
===================

Numba recently gained support for compiling generator functions.  This
document explains some of the implementation choices.


Terminology
===========

For clarity, we distinguish between *generator functions* and
*generators*.  A generator function is a function containing one or
several ``yield`` statements.  A generator (sometimes also called "generator
iterator") is the return value of a generator function; it resumes
execution inside its frame each time :py:func:`next` is called.

A *yield point* is the place where a ``yield`` statement is called.
A *resumption point* is the place just after a *yield point* where execution
is resumed when :py:func:`next` is called again.


Function analysis
=================

Suppose we have the following simple generator function::

   def gen(x, y):
       yield x + y
       yield x - y

Here is its CPython bytecode, as printed out using :py:func:`dis.dis`::

  7           0 LOAD_FAST                0 (x)
              3 LOAD_FAST                1 (y)
              6 BINARY_ADD
              7 YIELD_VALUE
              8 POP_TOP

  8           9 LOAD_FAST                0 (x)
             12 LOAD_FAST                1 (y)
             15 BINARY_SUBTRACT
             16 YIELD_VALUE
             17 POP_TOP
             18 LOAD_CONST               0 (None)
             21 RETURN_VALUE

When compiling this function with :envvar:`NUMBA_DUMP_IR` set to 1, the
following information is printed out::

   ----------------------------------IR DUMP: gen----------------------------------
   label 0:
       x = arg(0, name=x)                       ['x']
       y = arg(1, name=y)                       ['y']
       $0.3 = x + y                             ['$0.3', 'x', 'y']
       $0.4 = yield $0.3                        ['$0.3', '$0.4']
       del $0.4                                 []
       del $0.3                                 []
       $0.7 = x - y                             ['$0.7', 'x', 'y']
       del y                                    []
       del x                                    []
       $0.8 = yield $0.7                        ['$0.7', '$0.8']
       del $0.8                                 []
       del $0.7                                 []
       $const0.9 = const(NoneType, None)        ['$const0.9']
       $0.10 = cast(value=$const0.9)            ['$0.10', '$const0.9']
       del $const0.9                            []
       return $0.10                             ['$0.10']
   ------------------------------GENERATOR INFO: gen-------------------------------
   generator state variables: ['$0.3', '$0.7', 'x', 'y']
   yield point #1: live variables = ['x', 'y'], weak live variables = ['$0.3']
   yield point #2: live variables = [], weak live variables = ['$0.7']


What does it mean? The first part is the Numba IR, as already seen in
:ref:`arch_generate_numba_ir`.  We can see the two yield points (``yield $0.3``
and ``yield $0.7``).

The second part shows generator-specific information.  To understand it
we have to understand what suspending and resuming a generator means.

When suspending a generator, we are not merely returning a value to the
caller (the operand of the ``yield`` statement).  We also have to save the
generator's *current state* in order to resume execution.  In trivial use
cases, perhaps the CPU's register values or stack slots would be preserved
until the next call to next().  However, any non-trivial case will hopelessly
clobber those values, so we have to save them in a well-defined place.

What are the values we need to save?  Well, in the context of the Numba
Intermediate Representation, we must save all *live variables* at each
yield point.  These live variables are computed thanks to the control
flow graph.

Once live variables are saved and the generator is suspended, resuming
the generator simply involves the inverse operation: the live variables
are restored from the saved generator state.

.. note::
   It is the same analysis which helps insert Numba ``del`` instructions
   where appropriate.

Let's go over the generator info again::

   generator state variables: ['$0.3', '$0.7', 'x', 'y']
   yield point #1: live variables = ['x', 'y'], weak live variables = ['$0.3']
   yield point #2: live variables = [], weak live variables = ['$0.7']

Numba has computed the union of all live variables (denoted as "state
variables").  This will help define the layout of the :ref:`generator
structure <generator-structure>`.  Also, for each yield point, we have
computed two sets of variables:

* the *live variables* are the variables which are used by code following
  the resumption point (i.e. after the ``yield`` statement)

* the *weak live variables* are variables which are del'ed immediately
  after the resumption point; they have to be saved in :term:`object mode`,
  to ensure proper reference cleanup


.. _generator-structure:

The generator structure
=======================

Layout
------

Function analysis helps us gather enough information to define the
layout of the generator structure, which will store the entire execution
state of a generator.  Here is a sketch of the generator structure's layout,
in pseudo-code::

   struct gen_struct_t {
      int32_t resume_index;
      struct gen_args_t {
         arg_0_t arg0;
         arg_1_t arg1;
         ...
         arg_N_t argN;
      }
      struct gen_state_t {
         state_0_t state_var0;
         state_1_t state_var1;
         ...
         state_N_t state_varN;
      }
   }

Let's describe those fields in order.

* The first member, the *resume index*, is an integer telling the generator
  at which resumption point execution must resume.  By convention, it can
  have two special values: 0 means execution must start at the beginning of
  the generator (i.e. the first time :py:func:`next` is called); -1 means
  the generator is exhausted and resumption must immediately raise
  StopIteration.  Other values indicate the yield point's index starting from 1
  (corresponding to the indices shown in the generator info above).

* The second member, the *arguments structure* is read-only after it is first
  initialized.  It stores the values of the arguments the generator function
  was called with.  In our example, these are the values of ``x`` and ``y``.

* The third member, the *state structure*, stores the live variables as
  computed above.

Concretely, our example's generator structure (assuming the generator
function is called with floating-point numbers) is then::

   struct gen_struct_t {
      int32_t resume_index;
      struct gen_args_t {
         double arg0;
         double arg1;
      }
      struct gen_state_t {
         double $0.3;
         double $0.7;
         double x;
         double y;
      }
   }

Note that here, saving ``x`` and ``y`` is redundant: Numba isn't able to
recognize that the state variables ``x`` and ``y`` have the same value
as ``arg0`` and ``arg1``.

Allocation
----------

How does Numba ensure the generator structure is preserved long enough?
There are two cases:

* When a Numba-compiled generator function is called from a Numba-compiled
  function, the structure is allocated on the stack by the callee.  In this
  case, generator instantiation is practically costless.

* When a Numba-compiled generator function is called from regular Python
  code, a CPython-compatible wrapper is instantiated that has the right
  amount of allocated space to store the structure, and whose
  :c:member:`~PyTypeObject.tp_iternext` slot is a wrapper around the
  generator's native code.


Compiling to native code
========================

When compiling a generator function, three native functions are actually
generated by Numba:

* An initialization function.  This is the function corresponding
  to the generator function itself: it receives the function arguments and
  stores them inside the generator structure (which is passed by pointer).
  It also initialized the *resume index* to 0, indicating that the generator
  hasn't started yet.

* A next() function.  This is the function called to resume execution
  inside the generator.  Its single argument is a pointer to the generator
  structure and it returns the next yielded value (or a special exit code
  is used if the generator is exhausted, for quick checking when called
  from Numba-compiled functions).

* An optional finalizer.  In object mode, this function ensures that all
  live variables stored in the generator state are decref'ed, even if the
  generator is destroyed without having been exhausted.

The next() function
-------------------

The next() function is the least straight-forward of the three native
functions.  It starts with a trampoline which dispatches execution to the
right resume point depending on the *resume index* stored in the generator
structure.  Here is how the function start may look like in our example:

.. code-block:: llvm

   define i32 @"__main__.gen.next"(
      double* nocapture %retptr,
      { i8*, i32 }** nocapture readnone %excinfo,
      i8* nocapture readnone %env,
      { i32, { double, double }, { double, double, double, double } }* nocapture %arg.gen)
   {
     entry:
        %gen.resume_index = getelementptr { i32, { double, double }, { double, double, double, double } }* %arg.gen, i64 0, i32 0
        %.47 = load i32* %gen.resume_index, align 4
        switch i32 %.47, label %stop_iteration [
          i32 0, label %B0
          i32 1, label %generator_resume1
          i32 2, label %generator_resume2
        ]

     ; rest of the function snipped

(uninteresting stuff trimmed from the LLVM IR to make it more readable)

We recognize the pointer to the generator structure in ``%arg.gen``.
The trampoline switch has three targets (one for each *resume index* 0, 1
and 2), and a fallback target label named ``stop_iteration``.  Label ``B0``
represents the function's start, ``generator_resume1`` (resp.
``generator_resume2``) is the resumption point after the first
(resp. second) yield point.

After generation by LLVM, the whole native assembly code for this function
may look like this (on x86-64):

.. code-block:: asm

           .globl  __main__.gen.next
           .align  16, 0x90
   __main__.gen.next:
           movl    (%rcx), %eax
           cmpl    $2, %eax
           je      .LBB1_5
           cmpl    $1, %eax
           jne     .LBB1_2
           movsd   40(%rcx), %xmm0
           subsd   48(%rcx), %xmm0
           movl    $2, (%rcx)
           movsd   %xmm0, (%rdi)
           xorl    %eax, %eax
           retq
   .LBB1_5:
           movl    $-1, (%rcx)
           jmp     .LBB1_6
   .LBB1_2:
           testl   %eax, %eax
           jne     .LBB1_6
           movsd   8(%rcx), %xmm0
           movsd   16(%rcx), %xmm1
           movaps  %xmm0, %xmm2
           addsd   %xmm1, %xmm2
           movsd   %xmm1, 48(%rcx)
           movsd   %xmm0, 40(%rcx)
           movl    $1, (%rcx)
           movsd   %xmm2, (%rdi)
           xorl    %eax, %eax
           retq
   .LBB1_6:
           movl    $-3, %eax
           retq

Note the function returns 0 to indicate a value is yielded, -3 to indicate
StopIteration. ``%rcx`` points to the start of the generator structure,
where the resume index is stored.
A Map of the Numba Repository
=============================

The Numba repository is quite large, and due to age has functionality spread
around many locations.  To help orient developers, this document will try to
summarize where different categories of functionality can be found.


Support Files
-------------

Build and Packaging
'''''''''''''''''''

- :ghfile:`setup.py` - Standard Python distutils/setuptools script
- :ghfile:`MANIFEST.in` - Distutils packaging instructions
- :ghfile:`requirements.txt` - Pip package requirements, not used by conda
- :ghfile:`versioneer.py` - Handles automatic setting of version in
  installed package from git tags
- :ghfile:`.flake8` - Preferences for code formatting.  Files should be
  fixed and removed from the exception list as time allows.
- :ghfile:`.pre-commit-config.yaml` - Configuration file for pre-commit hooks.
- :ghfile:`.readthedocs.yml` - Configuration file for Read the Docs.
- :ghfile:`buildscripts/condarecipe.local` - Conda build recipe
- :ghfile:`buildscripts/condarecipe_clone_icc_rt` - Recipe to build a
  standalone icc_rt package.


Continuous Integration
''''''''''''''''''''''
- :ghfile:`azure-pipelines.yml` - Azure Pipelines CI config (active:
  Win/Mac/Linux)
- :ghfile:`buildscripts/azure/` - Azure Pipeline configuration for specific
  platforms
- :ghfile:`buildscripts/appveyor/` - Appveyor build scripts
- :ghfile:`buildscripts/incremental/` - Generic scripts for building Numba
  on various CI systems
- :ghfile:`codecov.yml` - Codecov.io coverage reporting


Documentation
'''''''''''''
- :ghfile:`LICENSE` - License for Numba
- :ghfile:`LICENSES.third-party` - License for third party code vendored
  into Numba
- :ghfile:`README.rst` - README for repo, also uploaded to PyPI
- :ghfile:`CONTRIBUTING.md` - Documentation on how to contribute to project
  (out of date, should be updated to point to Sphinx docs)
- :ghfile:`CHANGE_LOG` - History of Numba releases, also directly embedded
  into Sphinx documentation
- :ghfile:`docs/` - Documentation source
- :ghfile:`docs/_templates/` - Directory for templates (to override defaults
  with Sphinx theme)
- :ghfile:`docs/Makefile` - Used to build Sphinx docs with ``make``
- :ghfile:`docs/source` - ReST source for Numba documentation
- :ghfile:`docs/_static/` - Static CSS and image assets for Numba docs
- :ghfile:`docs/gh-pages.py` - Utility script to update Numba docs (stored
  as gh-pages)
- :ghfile:`docs/make.bat` - Not used (remove?)
- :ghfile:`docs/requirements.txt` - Pip package requirements for building docs
  with Read the Docs.
- :ghfile:`numba/scripts/generate_lower_listing.py` - Dump all registered
  implementations decorated with ``@lower*`` for reference
  documentation.  Currently misses implementations from the higher
  level extension API.


Numba Source Code
-----------------

Numba ships with both the source code and tests in one package.

- :ghfile:`numba/` - all of the source code and tests


Public API
''''''''''

These define aspects of the public Numba interface.

- :ghfile:`numba/core/decorators.py` - User-facing decorators for compiling
  regular functions on the CPU
- :ghfile:`numba/core/extending.py` - Public decorators for extending Numba
  (``overload``, ``intrinsic``, etc)
  - :ghfile:`numba/experimental/structref.py` - Public API for defining a mutable struct
- :ghfile:`numba/core/ccallback.py` - ``@cfunc`` decorator for compiling
  functions to a fixed C signature.  Used to make callbacks.
- :ghfile:`numba/np/ufunc/decorators.py` - ufunc/gufunc compilation
  decorators
- :ghfile:`numba/core/config.py` - Numba global config options and environment
  variable handling
- :ghfile:`numba/core/annotations` - Gathering and printing type annotations of
  Numba IR
- :ghfile:`numba/core/annotations/pretty_annotate.py` - Code highlighting of
  Numba functions and types (both ANSI terminal and HTML)
- :ghfile:`numba/core/event.py` - A simple event system for applications to
  listen to specific compiler events.


Dispatching
'''''''''''

- :ghfile:`numba/core/dispatcher.py` - Dispatcher objects are compiled functions
  produced by ``@jit``.  A dispatcher has different implementations
  for different type signatures.
- :ghfile:`numba/_dispatcher.cpp` - C++ dispatcher implementation (for speed on
  common data types)
- :ghfile:`numba/core/retarget.py` - Support for dispatcher objects to switch
  target via a specific with-context.


Compiler Pipeline
'''''''''''''''''

- :ghfile:`numba/core/compiler.py` - Compiler pipelines and flags
- :ghfile:`numba/core/errors.py` - Numba exception and warning classes
- :ghfile:`numba/core/ir.py` - Numba IR data structure objects
- :ghfile:`numba/core/bytecode.py` - Bytecode parsing and function identity (??)
- :ghfile:`numba/core/interpreter.py` - Translate Python interpreter bytecode to
  Numba IR
- :ghfile:`numba/core/analysis.py` - Utility functions to analyze Numba IR
  (variable lifetime, prune branches, etc)
- :ghfile:`numba/core/dataflow.py` - Dataflow analysis for Python bytecode (used
  in analysis.py)
- :ghfile:`numba/core/controlflow.py` - Control flow analysis of Numba IR and
  Python bytecode
- :ghfile:`numba/core/typeinfer.py` - Type inference algorithm
- :ghfile:`numba/core/transforms.py` - Numba IR transformations
- :ghfile:`numba/core/rewrites` - Rewrite passes used by compiler
- :ghfile:`numba/core/rewrites/__init__.py` - Loads all rewrite passes so they
  are put into the registry
- :ghfile:`numba/core/rewrites/registry.py` - Registry object for collecting
  rewrite passes
- :ghfile:`numba/core/rewrites/ir_print.py` - Write print() calls into special
  print nodes in the IR
- :ghfile:`numba/core/rewrites/static_raise.py` - Converts exceptions with
  static arguments into a special form that can be lowered
- :ghfile:`numba/core/rewrites/static_getitem.py` - Rewrites getitem and setitem
  with constant arguments to allow type inference
- :ghfile:`numba/core/rewrites/static_binop.py` - Rewrites binary operations
  (specifically ``**``) with constant arguments so faster code can be
  generated
- :ghfile:`numba/core/inline_closurecall.py` - Inlines body of closure functions
  to call site.  Support for array comprehensions, reduction inlining,
  and stencil inlining.
- :ghfile:`numba/core/postproc.py` - Postprocessor for Numba IR that computes
  variable lifetime, inserts del operations, and handles generators
- :ghfile:`numba/core/lowering.py` - General implementation of lowering Numba IR
  to LLVM
  :ghfile:`numba/core/environment.py` - Runtime environment object
- :ghfile:`numba/core/withcontexts.py` - General scaffolding for implementing
  context managers in nopython mode, and the objectmode context
  manager
- :ghfile:`numba/core/pylowering.py` - Lowering of Numba IR in object mode
- :ghfile:`numba/core/pythonapi.py` - LLVM IR code generation to interface with
  CPython API
- :ghfile:`numba/core/targetconfig.py` - Utils for target configurations such
  as compiler flags.


Type Management
'''''''''''''''

- :ghfile:`numba/core/typeconv/` - Implementation of type casting and type
  signature matching in both C++ and Python
- :ghfile:`numba/capsulethunk.h` - Used by typeconv
- :ghfile:`numba/core/types/` - definition of the Numba type hierarchy, used
  everywhere in compiler to select implementations
- :ghfile:`numba/core/consts.py` - Constant inference (used to make constant
  values available during codegen when possible)
- :ghfile:`numba/core/datamodel` - LLVM IR representations of data types in
  different contexts
- :ghfile:`numba/core/datamodel/models.py` - Models for most standard types
- :ghfile:`numba/core/datamodel/registry.py` - Decorator to register new data
  models
- :ghfile:`numba/core/datamodel/packer.py` - Pack typed values into a data
  structure
- :ghfile:`numba/core/datamodel/testing.py` - Data model tests (this should
  move??)
- :ghfile:`numba/core/datamodel/manager.py` - Map types to data models


Compiled Extensions
'''''''''''''''''''

Numba uses a small amount of compiled C/C++ code for core
functionality, like dispatching and type matching where performance
matters, and it is more convenient to encapsulate direct interaction
with CPython APIs.

- :ghfile:`numba/_arraystruct.h` - Struct for holding NumPy array
  attributes.  Used in helperlib and the Numba Runtime.
- :ghfile:`numba/_helperlib.c` - C functions required by Numba compiled code
  at runtime.  Linked into ahead-of-time compiled modules
- :ghfile:`numba/_helpermod.c` - Python extension module with pointers to
  functions from ``_helperlib.c`` and ``_npymath_exports.c``
- :ghfile:`numba/_npymath_exports.c` - Export function pointer table to
  NumPy C math functions
- :ghfile:`numba/_dynfuncmod.c` - Python extension module exporting
  _dynfunc.c functionality
- :ghfile:`numba/_dynfunc.c` - C level Environment and Closure objects (keep
  in sync with numba/target/base.py)
- :ghfile:`numba/mathnames.h` - Macros for defining names of math functions
- :ghfile:`numba/_pymodule.h` - C macros for Python 2/3 portable naming of C
  API functions
- :ghfile:`numba/mviewbuf.c` - Handles Python memoryviews
- :ghfile:`numba/_typeof.{h,c}` - C implementation of type fingerprinting,
  used by dispatcher
- :ghfile:`numba/_numba_common.h` - Portable C macro for marking symbols
  that can be shared between object files, but not outside the
  library.



Misc Support
''''''''''''

- :ghfile:`numba/_version.py` - Updated by versioneer
- :ghfile:`numba/core/runtime` - Language runtime.  Currently manages
  reference-counted memory allocated on the heap by Numba-compiled
  functions
- :ghfile:`numba/core/ir_utils.py` - Utility functions for working with Numba IR
  data structures
- :ghfile:`numba/core/cgutils.py` - Utility functions for generating common code
  patterns in LLVM IR
- :ghfile:`numba/core/utils.py` - Python 2 backports of Python 3 functionality
  (also imports local copy of ``six``)
- :ghfile:`numba/core/overload_glue.py` - Functions for wrapping split typing
  and lowering API use cases into overloads.
- :ghfile:`numba/misc/appdirs.py` - Vendored package for determining application
  config directories on every platform
- :ghfile:`numba/core/compiler_lock.py` - Global compiler lock because Numba's
  usage of LLVM is not thread-safe
- :ghfile:`numba/misc/special.py` - Python stub implementations of special Numba
  functions (prange, gdb*)
- :ghfile:`numba/core/itanium_mangler.py` - Python implementation of Itanium C++
  name mangling
- :ghfile:`numba/misc/findlib.py` - Helper function for locating shared
  libraries on all platforms
- :ghfile:`numba/core/debuginfo.py` - Helper functions to construct LLVM IR
  debug
  info
- :ghfile:`numba/core/unsafe/refcount.py` - Read reference count of object
- :ghfile:`numba/core/unsafe/eh.py` - Exception handling helpers
- :ghfile:`numba/core/unsafe/nrt.py` - Numba runtime (NRT) helpers
- :ghfile:`numba/cpython/unsafe/tuple.py` - Replace a value in a tuple slot
- :ghfile:`numba/np/unsafe/ndarray.py` - NumPy array helpers
- :ghfile:`numba/core/unsafe/bytes.py` - Copying and dereferencing data from
  void pointers
- :ghfile:`numba/misc/dummyarray.py` - Used by GPU backends to hold array
  information on the host, but not the data.
- :ghfile:`numba/core/callwrapper.py` - Handles argument unboxing and releasing
  the GIL when moving from Python to nopython mode
- :ghfile:`numba/np/numpy_support.py` - Helper functions for working with NumPy
  and translating Numba types to and from NumPy dtypes.
- :ghfile:`numba/core/tracing.py` - Decorator for tracing Python calls and
  emitting log messages
- :ghfile:`numba/core/funcdesc.py` - Classes for describing function metadata
  (used in the compiler)
- :ghfile:`numba/core/sigutils.py` - Helper functions for parsing and
  normalizing Numba type signatures
- :ghfile:`numba/core/serialize.py` - Support for pickling compiled functions
- :ghfile:`numba/core/caching.py` - Disk cache for compiled functions
- :ghfile:`numba/np/npdatetime.py` - Helper functions for implementing NumPy
  datetime64 support
- :ghfile:`numba/misc/llvm_pass_timings.py` - Helper to record timings of
  LLVM passes.
- :ghfile:`numba/cloudpickle` - Vendored cloudpickle subpackage

Core Python Data Types
''''''''''''''''''''''

- :ghfile:`numba/_hashtable.{h,c}` - Adaptation of the Python 3.7 hash table
  implementation
- :ghfile:`numba/cext/dictobject.{h,c}` - C level implementation of typed
  dictionary
- :ghfile:`numba/typed/dictobject.py` - Nopython mode wrapper for typed
  dictionary
- :ghfile:`numba/cext/listobject.{h,c}` - C level implementation of typed list
- :ghfile:`numba/typed/listobject.py` - Nopython mode wrapper for typed list
- :ghfile:`numba/typed/typedobjectutils.py` - Common utilities for typed
  dictionary and list
- :ghfile:`numba/cpython/unicode.py` - Unicode strings (Python 3.5 and later)
- :ghfile:`numba/typed` - Python interfaces to statically typed containers
- :ghfile:`numba/typed/typeddict.py` - Python interface to typed dictionary
- :ghfile:`numba/typed/typedlist.py` - Python interface to typed list
- :ghfile:`numba/experimental/jitclass` - Implementation of experimental JIT
  compilation of Python classes
- :ghfile:`numba/core/generators.py` - Support for lowering Python generators


Math
''''

- :ghfile:`numba/_random.c` - Reimplementation of NumPy / CPython random
  number generator
- :ghfile:`numba/_lapack.c` - Wrappers for calling BLAS and LAPACK functions
  (requires SciPy)


ParallelAccelerator
'''''''''''''''''''

Code transformation passes that extract parallelizable code from
a function and convert it into multithreaded gufunc calls.

- :ghfile:`numba/parfors/parfor.py` - General ParallelAccelerator
- :ghfile:`numba/parfors/parfor_lowering.py` - gufunc lowering for
  ParallelAccelerator
- :ghfile:`numba/parfors/array_analysis.py` - Array analysis passes used in
  ParallelAccelerator


Stencil
'''''''

Implementation of ``@stencil``:

- :ghfile:`numba/stencils/stencil.py` - Stencil function decorator (implemented
  without ParallelAccelerator)
- :ghfile:`numba/stencils/stencilparfor.py` - ParallelAccelerator implementation
  of stencil


Debugging Support
'''''''''''''''''

- :ghfile:`numba/misc/gdb_hook.py` - Hooks to jump into GDB from nopython
  mode
- :ghfile:`numba/misc/cmdlang.gdb` - Commands to setup GDB for setting
  explicit breakpoints from Python


Type Signatures (CPU)
'''''''''''''''''''''

Some (usually older) Numba supported functionality separates the
declaration of allowed type signatures from the definition of
implementations.  This package contains registries of type signatures
that must be matched during type inference.

- :ghfile:`numba/core/typing` - Type signature module
- :ghfile:`numba/core/typing/templates.py` - Base classes for type signature
  templates
- :ghfile:`numba/core/typing/cmathdecl.py` - Python complex math (``cmath``)
  module
- :ghfile:`numba/core/typing/bufproto.py` - Interpreting objects supporting the
  buffer protocol
- :ghfile:`numba/core/typing/mathdecl.py` - Python ``math`` module
- :ghfile:`numba/core/typing/listdecl.py` - Python lists
- :ghfile:`numba/core/typing/builtins.py` - Python builtin global functions and
  operators
- :ghfile:`numba/core/typing/randomdecl.py` - Python and NumPy ``random``
  modules
- :ghfile:`numba/core/typing/setdecl.py` - Python sets
- :ghfile:`numba/core/typing/npydecl.py` - NumPy ndarray (and operators), NumPy
  functions
- :ghfile:`numba/core/typing/arraydecl.py` - Python ``array`` module
- :ghfile:`numba/core/typing/context.py` - Implementation of typing context
  (class that collects methods used in type inference)
- :ghfile:`numba/core/typing/collections.py` - Generic container operations and
  namedtuples
- :ghfile:`numba/core/typing/ctypes_utils.py` - Typing ctypes-wrapped function
  pointers
- :ghfile:`numba/core/typing/enumdecl.py` - Enum types
- :ghfile:`numba/core/typing/cffi_utils.py` - Typing of CFFI objects
- :ghfile:`numba/core/typing/typeof.py` - Implementation of typeof operations
  (maps Python object to Numba type)
- :ghfile:`numba/core/typing/asnumbatype.py` - Implementation of
  ``as_numba_type`` operations (maps Python types to Numba type)
- :ghfile:`numba/core/typing/npdatetime.py` - Datetime dtype support for NumPy
  arrays


Target Implementations (CPU)
''''''''''''''''''''''''''''

Implementations of Python / NumPy functions and some data models.
These modules are responsible for generating LLVM IR during lowering.
Note that some of these modules do not have counterparts in the typing
package because newer Numba extension APIs (like overload) allow
typing and implementation to be specified together.

- :ghfile:`numba/core/cpu.py` - Context for code gen on CPU
- :ghfile:`numba/core/base.py` - Base class for all target contexts
- :ghfile:`numba/core/codegen.py` - Driver for code generation
- :ghfile:`numba/core/boxing.py` - Boxing and unboxing for most data
  types
- :ghfile:`numba/core/intrinsics.py` - Utilities for converting LLVM
  intrinsics to other math calls
- :ghfile:`numba/core/callconv.py` - Implements different calling
  conventions for Numba-compiled functions
- :ghfile:`numba/core/options.py` - Container for options that control
  lowering
- :ghfile:`numba/core/optional.py` - Special type representing value or
  ``None``
- :ghfile:`numba/core/registry.py` - Registry object for collecting
  implementations for a specific target
- :ghfile:`numba/core/imputils.py` - Helper functions for lowering
- :ghfile:`numba/core/externals.py` - Registers external C functions
  needed to link generated code
- :ghfile:`numba/core/fastmathpass.py` - Rewrite pass to add fastmath
  attributes to function call sites and binary operations
- :ghfile:`numba/core/removerefctpass.py` - Rewrite pass to remove
  unnecessary incref/decref pairs
- :ghfile:`numba/core/descriptors.py` - empty base class for all target
  descriptors (is this needed?)
- :ghfile:`numba/cpython/builtins.py` - Python builtin functions and
  operators
- :ghfile:`numba/cpython/cmathimpl.py` - Python complex math module
- :ghfile:`numba/cpython/enumimpl.py` - Enum objects
- :ghfile:`numba/cpython/hashing.py` - Hashing algorithms
- :ghfile:`numba/cpython/heapq.py` - Python ``heapq`` module
- :ghfile:`numba/cpython/iterators.py` - Iterable data types and iterators
- :ghfile:`numba/cpython/listobj.py` - Python lists
- :ghfile:`numba/cpython/mathimpl.py` - Python ``math`` module
- :ghfile:`numba/cpython/numbers.py` - Numeric values (int, float, etc)
- :ghfile:`numba/cpython/printimpl.py` - Print function
- :ghfile:`numba/cpython/randomimpl.py` - Python and NumPy ``random``
  modules
- :ghfile:`numba/cpython/rangeobj.py` - Python `range` objects
- :ghfile:`numba/cpython/slicing.py` - Slice objects, and index calculations
  used in slicing
- :ghfile:`numba/cpython/setobj.py` - Python set type
- :ghfile:`numba/cpython/tupleobj.py` - Tuples (statically typed as
  immutable struct)
- :ghfile:`numba/misc/cffiimpl.py` - CFFI functions
- :ghfile:`numba/misc/quicksort.py` - Quicksort implementation used with
  list and array objects
- :ghfile:`numba/misc/mergesort.py` - Mergesort implementation used with
  array objects
- :ghfile:`numba/np/arraymath.py` - Math operations on arrays (both
  Python and NumPy)
- :ghfile:`numba/np/arrayobj.py` - Array operations (both NumPy and
  buffer protocol)
- :ghfile:`numba/np/linalg.py` - NumPy linear algebra operations
- :ghfile:`numba/np/npdatetime.py` - NumPy datetime operations
- :ghfile:`numba/np/npyfuncs.py` - Kernels used in generating some
  NumPy ufuncs
- :ghfile:`numba/np/npyimpl.py` - Implementations of most NumPy ufuncs
- :ghfile:`numba/np/polynomial.py` - ``numpy.roots`` function
- :ghfile:`numba/np/ufunc_db.py` - Big table mapping types to ufunc
  implementations


Ufunc Compiler and Runtime
''''''''''''''''''''''''''

- :ghfile:`numba/np/ufunc` - ufunc compiler implementation
- :ghfile:`numba/np/ufunc/_internal.{h,c}` - Python extension module with
  helper functions that use CPython & NumPy C API
- :ghfile:`numba/np/ufunc/_ufunc.c` - Used by `_internal.c`
- :ghfile:`numba/np/ufunc/deviceufunc.py` - Custom ufunc dispatch for
  non-CPU targets
- :ghfile:`numba/np/ufunc/gufunc_scheduler.{h,cpp}` - Schedule work chunks
  to threads
- :ghfile:`numba/np/ufunc/dufunc.py` - Special ufunc that can compile new
  implementations at call time
- :ghfile:`numba/np/ufunc/ufuncbuilder.py` - Top-level orchestration of
  ufunc/gufunc compiler pipeline
- :ghfile:`numba/np/ufunc/sigparse.py` - Parser for generalized ufunc
  indexing signatures
- :ghfile:`numba/np/ufunc/parallel.py` - Codegen for ``parallel`` target
- :ghfile:`numba/np/ufunc/array_exprs.py` - Rewrite pass for turning array
  expressions in regular functions into ufuncs
- :ghfile:`numba/np/ufunc/wrappers.py` - Wrap scalar function kernel with
  loops
- :ghfile:`numba/np/ufunc/workqueue.{h,c}` - Threading backend based on
  pthreads/Windows threads and queues
- :ghfile:`numba/np/ufunc/omppool.cpp` - Threading backend based on OpenMP
- :ghfile:`numba/np/ufunc/tbbpool.cpp` - Threading backend based on TBB



Unit Tests (CPU)
''''''''''''''''

CPU unit tests (GPU target unit tests listed in later sections

- :ghfile:`runtests.py` - Convenience script that launches test runner and
  turns on full compiler tracebacks
- :ghfile:`run_coverage.py` - Runs test suite with coverage tracking enabled
- :ghfile:`.coveragerc` - Coverage.py configuration
- :ghfile:`numba/runtests.py` - Entry point to unittest runner
- :ghfile:`numba/testing/_runtests.py` - Implementation of custom test runner
  command line interface
- :ghfile:`numba/tests/test_*` - Test cases
- :ghfile:`numba/tests/*_usecases.py` - Python functions compiled by some
  unit tests
- :ghfile:`numba/tests/support.py` - Helper functions for testing and
  special TestCase implementation
- :ghfile:`numba/tests/dummy_module.py` - Module used in
  ``test_dispatcher.py``
- :ghfile:`numba/tests/npyufunc` - ufunc / gufunc compiler tests
- :ghfile:`numba/testing` - Support code for testing
- :ghfile:`numba/testing/loader.py` - Find tests on disk
- :ghfile:`numba/testing/notebook.py` - Support for testing notebooks
- :ghfile:`numba/testing/main.py` - Numba test runner


Command Line Utilities
''''''''''''''''''''''
- :ghfile:`bin/numba` - Command line stub, delegates to main in
  ``numba_entry.py``
- :ghfile:`numba/misc/numba_entry.py` - Main function for ``numba`` command line
  tool
- :ghfile:`numba/pycc` - Ahead of time compilation of functions to shared
  library extension
- :ghfile:`numba/pycc/__init__.py` - Main function for ``pycc`` command line
  tool
- :ghfile:`numba/pycc/cc.py` - User-facing API for tagging functions to
  compile ahead of time
- :ghfile:`numba/pycc/compiler.py` - Compiler pipeline for creating
  standalone Python extension modules
- :ghfile:`numba/pycc/llvm_types.py` - Aliases to LLVM data types used by
  ``compiler.py``
- :ghfile:`numba/pycc/pycc` - Stub to call main function.  Is this still
  used?
- :ghfile:`numba/pycc/modulemixin.c` - C file compiled into every compiled
  extension.  Pulls in C source from Numba core that is needed to make
  extension standalone
- :ghfile:`numba/pycc/platform.py` - Portable interface to platform-specific
  compiler toolchains
- :ghfile:`numba/pycc/decorators.py` - Deprecated decorators for tagging
  functions to compile.  Use ``cc.py`` instead.


CUDA GPU Target
'''''''''''''''

Note that the CUDA target does reuse some parts of the CPU target.

- :ghfile:`numba/cuda/` - The implementation of the CUDA (NVIDIA GPU) target
  and associated unit tests
- :ghfile:`numba/cuda/decorators.py` - Compiler decorators for CUDA kernels
  and device functions
- :ghfile:`numba/cuda/dispatcher.py` - Dispatcher for CUDA JIT functions
- :ghfile:`numba/cuda/printimpl.py` - Special implementation of device printing
- :ghfile:`numba/cuda/libdevice.py` - Registers libdevice functions
- :ghfile:`numba/cuda/kernels/` - Custom kernels for reduction and transpose
- :ghfile:`numba/cuda/device_init.py` - Initializes the CUDA target when
  imported
- :ghfile:`numba/cuda/compiler.py` - Compiler pipeline for CUDA target
- :ghfile:`numba/cuda/intrinsic_wrapper.py` - CUDA device intrinsics
  (shuffle, ballot, etc)
- :ghfile:`numba/cuda/initialize.py` - Defered initialization of the CUDA
  device and subsystem.  Called only when user imports ``numba.cuda``
- :ghfile:`numba/cuda/simulator_init.py` - Initalizes the CUDA simulator
  subsystem (only when user requests it with env var)
- :ghfile:`numba/cuda/random.py` - Implementation of random number generator
- :ghfile:`numba/cuda/api.py` - User facing APIs imported into ``numba.cuda.*``
- :ghfile:`numba/cuda/stubs.py` - Python placeholders for functions that
  only can be used in GPU device code
- :ghfile:`numba/cuda/simulator/` - Simulate execution of CUDA kernels in
  Python interpreter
- :ghfile:`numba/cuda/vectorizers.py` - Subclasses of ufunc/gufunc compilers
  for CUDA
- :ghfile:`numba/cuda/args.py` - Management of kernel arguments, including
  host<->device transfers
- :ghfile:`numba/cuda/target.py` - Typing and target contexts for GPU
- :ghfile:`numba/cuda/cudamath.py` - Type signatures for math functions in
  CUDA Python
- :ghfile:`numba/cuda/errors.py` - Validation of kernel launch configuration
- :ghfile:`numba/cuda/nvvmutils.py` - Helper functions for generating
  NVVM-specific IR
- :ghfile:`numba/cuda/testing.py` - Support code for creating CUDA unit
  tests and capturing standard out
- :ghfile:`numba/cuda/cudadecl.py` - Type signatures of CUDA API (threadIdx,
  blockIdx, atomics) in Python on GPU
- :ghfile:`numba/cuda/cudaimpl.py` - Implementations of CUDA API functions
  on GPU
- :ghfile:`numba/cuda/codegen.py` - Code generator object for CUDA target
- :ghfile:`numba/cuda/cudadrv/` - Wrapper around CUDA driver API
- :ghfile:`numba/cuda/tests/` - CUDA unit tests, skipped when CUDA is not
  detected
- :ghfile:`numba/cuda/tests/cudasim/` - Tests of CUDA simulator
- :ghfile:`numba/cuda/tests/nocuda/` - Tests for NVVM functionality when
  CUDA not present
- :ghfile:`numba/cuda/tests/cudapy/` - Tests of compiling Python functions
  for GPU
- :ghfile:`numba/cuda/tests/cudadrv/` - Tests of Python wrapper around CUDA
  API

.. Copyright (c) 2017 Intel Corporation
   SPDX-License-Identifier: BSD-2-Clause

.. _arch-stencil:

=================
Notes on stencils 
=================

Numba provides the :ref:`@stencil decorator <numba-stencil>` to
represent stencil computations.  This document explains how this
feature is implemented in the several different modes available in
Numba.  Currently, calls to the stencil from non-jitted code is
supported as well as calls from jitted code, either with or without
the :ref:`parallel=True <parallel_jit_option>` option.

The stencil decorator
=====================

The stencil decorator itself just returns a ``StencilFunc`` object.
This object encapsulates the original stencil kernel function
as specified in the program and the options passed to the
stencil decorator.  Also of note is that after the first compilation
of the stencil, the computed neighborhood of the stencil is
stored in the ``StencilFunc`` object in the ``neighborhood`` attribute.

Handling the three modes
========================

As mentioned above, Numba supports the calling of stencils
from inside or outside a ``@jit`` compiled function, with or
without the :ref:`parallel=True <parallel_jit_option>` option.

Outside jit context
-------------------

``StencilFunc`` overrides the ``__call__`` method so that calls
to ``StencilFunc`` objects execute the stencil::

    def __call__(self, *args, **kwargs):                                        
        result = kwargs.get('out')
                                                                                
        new_stencil_func = self._stencil_wrapper(result, None, *args)           
                                                                                
        if result is None:                                                      
            return new_stencil_func.entry_point(*args)                          
        else:                                                                   
            return new_stencil_func.entry_point(*args, result)                  

First, the presence of the optional :ref:`out <stencil-function-out>`
parameter is checked.  If it is present then the output array is
stored in ``result``.  Then, the call to ``_stencil_wrapper``
generates the stencil function given the result and argument types
and finally the generated stencil function is executed and its result
returned.

Jit without ``parallel=True``
-----------------------------

When constructed, a ``StencilFunc`` inserts itself into the typing
context's set of user functions and provides the ``_type_me``
callback.  In this way, the standard Numba compiler is able to
determine the output type and signature of a ``StencilFunc``.
Each ``StencilFunc`` maintains a cache of previously seen combinations
of input argument types and keyword types.  If previously seen,
the ``StencilFunc`` returns the computed signature.  If not previously
computed, the ``StencilFunc`` computes the return type of the stencil
by running the Numba compiler frontend on the stencil kernel and
then performing type inference on the :term:`Numba IR` (IR) to get the scalar
return type of the kernel.  From that, a Numpy array type is constructed
whose element type matches that scalar return type.

After computing the signature of the stencil for a previously
unseen combination of input and keyword types, the ``StencilFunc``
then :ref:`creates the stencil function <arch-stencil-create-function>` itself.
``StencilFunc`` then installs the new stencil function's definition
in the target context so that jitted code is able to call it.

Thus, in this mode, the generated stencil function is a stand-alone
function called like a normal function from within jitted code.

Jit with ``parallel=True``
--------------------------

When calling a ``StencilFunc`` from a jitted context with ``parallel=True``,
a separate stencil function as generated by :ref:`arch-stencil-create-function`
is not used.  Instead, `parfors` (:ref:`parallel-accelerator`) are
created within the current function that implements the stencil.
This code again starts with the stencil kernel and does a similar kernel
size computation but then rather than standard Python looping syntax,
corresponding `parfors` are created so that the execution of the stencil
will take place in parallel.

The stencil to `parfor` translations can also be selectively disabled
by setting ``parallel={'stencil': False}``, among other sub-options
described in :ref:`parallel-accelerator`.

.. _arch-stencil-create-function:

Creating the stencil function
=============================

Conceptually, a stencil function is created from the user-specified
stencil kernel by adding looping code around the kernel, transforming
the relative kernel indices into absolute array indices based on the
loop indices, and replacing the kernel's ``return`` statement with
a statement to assign the computed value into the output array.

To accomplish this transformation, first, a copy of the stencil 
kernel IR is created so that subsequent modifications of the IR
for different stencil signatures will not effect each other.

Then, an approach similar to how GUFunc's are created for `parfors`
is employed.  In a text buffer, a Python function is created with
a unique name.  The input array parameter is added to the function
definition and if the ``out`` argument type is present then an
``out`` parameter is added to the stencil function definition.
If the ``out`` argument is not present then first an output array
is created with ``numpy.zeros`` having the same shape as the
input array.

The kernel is then analyzed to compute the stencil size and the
shape of the boundary (or the ``neighborhood`` stencil decorator
argument is used for this purpose if present).
Then, one ``for`` loop for each dimension of the input array is
added to the stencil function definition.  The range of each
loop is controlled by the stencil kernel size previously computed
so that the boundary of the output image is not modified but instead
left as is.   The body of the innermost ``for`` loop is a single
``sentinel`` statement that is easily recognized in the IR.
A call to ``exec`` with the text buffer is used to force the
stencil function into existence and an ``eval`` is used to get
access to the corresponding function on which ``run_frontend`` is
used to get the stencil function IR.

Various renaming and relabeling is performed on the stencil function
IR and the kernel IR so that the two can be combined without conflict.
The relative indices in the kernel IR (i.e., ``getitem`` calls) are
replaced with expressions where the corresponding loop index variables
are added to the relative indices.  The ``return`` statement in the
kernel IR is replaced with a ``setitem`` for the corresponding element
in the output array.
The stencil function IR is then scanned for the sentinel and the
sentinel replaced with the modified kernel IR.

Next, ``compile_ir`` is used to compile the combined stencil function
IR.  The resulting compile result is cached in the ``StencilFunc`` so that
other calls to the same stencil do not need to undertake this process
again.

Exceptions raised
=================

Various checks are performed during stencil compilation to make sure
that user-specified options do not conflict with each other or with
other runtime parameters.  For example, if the user has manually
specified a ``neighborhood`` to the stencil decorator, the length of
that neighborhood must match the dimensionality of the input array.
If this is not the case, a ``ValueError`` is raised.

If the neighborhood has not been specified then it must be inferred
and a requirement to infer the kernel is that all indices are constant
integers.  If they are not, a ``ValueError`` is raised indicating that
kernel indices may not be non-constant.

Finally, the stencil implementation detects the output array type
by running Numba type inference on the stencil kernel.  If the
return type of this kernel does not match the type of the value
passed to the ``cval`` stencil decorator option then a ``ValueError``
is raised.

=======================
Polymorphic dispatching
=======================

Functions compiled using :func:`~numba.jit` or :func:`~numba.vectorize`
are open-ended: they can be called with many different input types and
have to select (possibly compile on-the-fly) the right low-level
specialization.  We hereby explain how this mechanism is implemented.


Requirements
============

JIT-compiled functions can take several arguments and each of them is
taken into account when selecting a specialization.  Thus it is a
form of multiple dispatch, more complex than single dispatch.

Each argument weighs in the selection based on its :ref:`Numba type
<numba-types>`.  Numba types are often more granular than Python types:
for example, Numba types Numpy arrays differently depending on their
dimensionality and their layout (C-contiguous, etc.).

Once a Numba type is inferred for each argument, a specialization must
be chosen amongst the available ones; or, if not suitable specialization
is found, a new one must be compiled.  This is not a trivial decision:
there can be multiple specializations compatible with a given concrete
signature (for example, say a two-argument function has compiled
specializations for ``(float64, float64)`` and ``(complex64, complex64)``,
and it is called with ``(float32, float32)``).

Therefore, there are two crucial steps in the dispatch mechanism:

1. infer the Numba types of the concrete arguments
2. select the best available specialization (or choose to compile a new one)
   for the inferred Numba types

Compile-time vs. run-time
-------------------------

This document discusses dispatching when it is done at runtime, i.e.
when a JIT-compiled function is called from pure Python.  In that context,
performance is important.  To stay in the realm of normal function call
overhead in Python, the overhead of dispatching should stay under a
microsecond.  Of course, *the faster the better*...

When a JIT-compiled function is called from another JIT-compiled
function (in :term:`nopython mode`), the polymorphism is resolved at
compile-time, using a non-performance critical mechanism, bearing zero
runtime performance overhead.

.. note::
   In practice, the performance-critical parts described here are coded in C.


Type resolution
===============

The first step is therefore to infer, at call-time, a Numba type for each
of the function's concrete arguments.  Given the finer granularity of
Numba types compared to Python types, one cannot simply lookup an object's
class and key a dictionary with it to obtain the corresponding Numba type.

Instead, there is a machinery to inspect the object and, based on its
Python type, query various properties to infer the appropriate Numba
type.  This can be more or less complex: for example, a Python ``int``
argument will always infer to a Numba ``intp`` (a pointer-sized integer),
but a Python ``tuple`` argument can infer to multiple Numba types (depending
on the tuple's size and the concrete type of each of its elements).

The Numba type system is high-level and written in pure Python; there is
a pure Python machinery, based on a generic function, to do said inference
(in :mod:`numba.typing.typeof`).  That machinery is used for compile-time
inference, e.g. on constants.  Unfortunately, it is too slow for run-time
value-based dispatching.  It is only used as a fallback for rarely used
(or difficult to infer) types, and exhibits multiple-microsecond overhead.

Typecodes
---------

The Numba type system is really too high-level to be manipulated efficiently
from C code.  Therefore, the C dispatching layer uses another representation
based on integer typecodes.  Each Numba type gets a unique integer typecode
when constructed; also, an interning system ensure no two instances of same
type are created.  The dispatching layer is therefore able to *eschew*
the overhead of the Numba type system by working with simple integer
typecodes, amenable to well-known optimizations (fast hash tables, etc.).

The goal of the type resolution step becomes: infer a Numba *typecode*
for each of the function's concrete arguments.  Ideally, it doesn't deal
with Numba types anymore...

Hard-coded fast paths
---------------------

While eschewing the abstraction and object-orientation overhead of the type
system, the integer typecodes still have the same conceptual complexity.
Therefore, an important technique to speed up inference is to first go
through checks for the most important types, and hard-code a fast resolution
for each of them.

Several types benefit from such an optimization, notably:

* basic Python scalars (``bool``, ``int``, ``float``, ``complex``);
* basic Numpy scalars (the various kinds of integer, floating-point,
  complex numbers);
* Numpy arrays of certain dimensionalities and basic element types.

Each of those fast paths ideally uses a hard-coded result value or a direct
table lookup after a few simple checks.

However, we can't apply that technique to all argument types; there would
be an explosion of ad-hoc internal caches, and it would become difficult to
maintain.  Besides, the recursive application of hard-coded fast paths
would not necessarily combine into a low overhead (in the nested tuple
case, for example).

Fingerprint-based typecode cache
--------------------------------

For non-so-trivial types (imagine a tuple, or a Numpy ``datetime64`` array,
for example), the hard-coded fast paths don't match.  Another mechanism
then kicks in, more generic.

The principle here is to examine each argument value, as the pure Python
machinery would do, and to describe its Numba type unambiguously.  The
difference is that *we don't actually compute a Numba type*.  Instead, we
compute a simple bytestring, a low-level possible denotation of that
Numba type: a *fingerprint*.  The fingerprint format is designed to be
short and extremely simple to compute from C code (in practice, it has
a bytecode-like format).

Once the fingerprint is computed, it is looked up in a cache mapping
fingerprints to typecodes.  The cache is a hash table, and the lookup
is fast thanks to the fingerprints being generally very short (rarely
more than 20 bytes).

If the cache lookup fails, the typecode must first be computed using the
slow pure Python machinery.  Luckily, this would only happen once: on
subsequent calls, the cached typecode would be returned for the given
fingerprint.

In rare cases, a fingerprint cannot be computed efficiently.  This is
the case for some types which cannot be easily inspected from C: for
example ``cffi`` function pointers.  Then, the slow Pure Python machinery
is invoked at each function call with such an argument.

.. note::
   Two fingerprints may denote a single Numba type.  This does not make
   the mechanism incorrect; it only creates more cache entries.


Summary
-------

Type resolution of a function argument involves the following mechanisms
in order:

* Try a few hard-coded fast paths, for common simple types.
* If the above failed, compute a fingerprint for the argument and lookup
  its typecode in a cache.
* If all the above failed, invoke the pure Python machinery which will
  determine a Numba type for the argument (and look up its typecode).


Specialization selection
========================

At the previous step, an integer typecode has been determined for each
concrete argument to the JIT-compiled function.  Now it remains to match
that concrete signature against each of the available specializations for
the function.  There can be three outcomes:

* There is a satisfying best match: the corresponding specialization
  is then invoked (it will handle argument unboxing and other details).
* There is a tie between two or more "best matches": an exception is raised,
  refusing to solve the ambiguity.
* There is no satisfying match: a new specialization is compiled tailored
  for the concrete argument types that were inferred.

The selection works by looping over all available specializations, and
computing the compatibility of each concrete argument type with the
corresponding type in the specialization's intended signature.  Specifically,
we are interested in:

1. Whether the concrete argument type is allowed to convert implicitly to
   the specialization's argument type;
2. If so, at what semantic (user-visible) cost the conversion comes.

Implicit conversion rules
-------------------------

There are five possible kinds of implicit conversion from a source type
to a destination type (note this is an asymmetric relationship):

1. *exact match*: the two types are identical; this is the ideal case,
   since the specialization would behave exactly as intended;
2. *same-kind promotion*: the two types belong to the same "kind" (for
   example ``int32`` and ``int64`` are two integer types), and the source
   type can be converted losslessly to the destination type (e.g. from
   ``int32`` to ``int64``, but not the reverse);
3. *safe conversion*: the two types belong to different kinds, but the
   source type can be reasonably converted to the destination type
   (e.g. from ``int32`` to ``float64``, but not the reverse);
4. *unsafe conversion*: a conversion is available from the source type
   to the destination type, but it may lose precision, magnitude, or
   another desirable quality.
5. *no conversion*: there is no correct or reasonably efficient way to
   convert between the two types (for example between an ``int64`` and a
   ``datetime64``, or a C-contiguous array and a Fortran-contiguous array).

When a specialization is examined, the latter two cases eliminate it from
the final choice: i.e. when at least one argument has *no conversion* or
only an *unsafe conversion* to the signature's argument type.

.. note::
   However, if the function is compiled with explicit signatures
   in the :func:`~numba.jit` call (and therefore it is not allowed to compile
   new specializations), *unsafe conversion* is allowed.

Candidates and best match
-------------------------

If a specialization is not eliminated by the rule above, it enters the
list of *candidates* for the final choice.  Those candidates are ranked
by an ordered 4-uple of integers: ``(number of unsafe conversions,
number of safe conversions, number of same-kind promotions, number of
exact matches)`` (note the sum of the tuple's elements is equal to the
number of arguments).  The best match is then the #1 result in sorted
ascending order, thereby preferring exact matches over promotions,
promotions over safe conversions, safe conversions over unsafe conversions.

Implementation
--------------

The above-described mechanism works on integer typecodes, not on Numba
types.  It uses an internal hash table storing the possible conversion
kind for each pair of compatible types.  The internal hash table is in part
built at startup (for built-in trivial types such as ``int32``, ``int64``
etc.), in part filled dynamically (for arbitrarily complex types such
as array types: for example to allow using a C-contiguous 2D array where
a function expects a non-contiguous 2D array).

Summary
-------

Selecting the right specialization involves the following steps:

* Examine each available specialization and match it against the concrete
  argument types.
* Eliminate any specialization where at least one argument doesn't offer
  sufficient compatibility.
* If there are remaining candidates, choose the best one in terms of
  preserving the types' semantics.


Miscellaneous
=============

Some `benchmarks of dispatch performance
<https://github.com/numba/numba-benchmark/blob/master/benchmarks/bench_dispatch.py>`_
exist in the `Numba benchmarks <https://github.com/numba/numba-benchmark>`_
repository.

Some unit tests of specific aspects of the machinery are available
in :mod:`numba.tests.test_typeinfer` and :mod:`numba.tests.test_typeof`.
Higher-level dispatching tests are in :mod:`numba.tests.test_dispatcher`.
.. _arch-numba-runtime:

======================
Notes on Numba Runtime
======================


The *Numba Runtime (NRT)* provides the language runtime to the *nopython mode*
Python subset.  NRT is a standalone C library with a Python binding.  This
allows :term:`NPM` runtime feature to be used without the GIL.  Currently, the 
only language feature implemented in NRT is memory management.


Memory Management
=================

NRT implements memory management for :term:`NPM` code.  It uses *atomic 
reference count* for threadsafe, deterministic memory management.  NRT maintains 
a separate ``MemInfo`` structure for storing information about each allocation.

Cooperating with CPython
------------------------

For NRT to cooperate with CPython, the NRT python binding provides adaptors for
converting python objects that export a memory region.  When such an
object is used as an argument to a :term:`NPM` function, a new ``MemInfo`` is 
created and it acquires a reference to the Python object.  When a :term:`NPM` 
value is returned to the Python interpreter, the associated ``MemInfo`` 
(if any) is checked.  If the ``MemInfo`` references a Python object, the 
underlying Python object is released and returned instead.  Otherwise, the 
``MemInfo`` is wrapped in a Python object and returned. Additional process 
maybe required depending on the type.

The current implementation supports Numpy array and any buffer-exporting types.


Compiler-side Cooperation
-------------------------

NRT reference counting requires the compiler to emit incref/decref operations
according to the usage.  When the reference count drops to zero, the compiler
must call the destructor routine in NRT.


.. _nrt-refct-opt-pass:

Optimizations
-------------

The compiler is allowed to emit incref/decref operations naively.  It relies
on an optimization pass to remove redundant reference count operations.

A new optimization pass is implemented in version 0.52.0 to remove reference
count operations that fall into the following four categories of control-flow
structure---per basic-block, diamond, fanout, fanout+raise. See the documentation
for :envvar:`NUMBA_LLVM_REFPRUNE_FLAGS` for their descriptions.

The old optimization pass runs at block level to avoid control flow analysis.
It depends on LLVM function optimization pass to simplify the control flow,
stack-to-register, and simplify instructions.  It works by matching and
removing incref and decref pairs within each block.  The old pass can be
enabled by setting :envvar:`NUMBA_LLVM_REFPRUNE_PASS` to `0`.

Important assumptions
---------------------

Both the old (pre-0.52.0) and the new (post-0.52.0) optimization passes assume
that the only function that can consume a reference is ``NRT_decref``.
It is important that there are no other functions that will consume references.
Since the passes operate on LLVM IR, the "functions" here are referring to any
callee in a LLVM call instruction.

To summarize, all functions exposed to the refcount optimization pass
**must not** consume counted references unless done so via ``NRT_decref``.


Quirks of the old optimization pass
-----------------------------------

Since the pre-0.52.0 `refcount optimization pass <nrt-refct-opt-pass_>`_
requires the LLVM function optimization pass, the pass works on the LLVM IR as
text. The optimized IR is then materialized again as a new LLVM in-memory
bitcode object.


Debugging Leaks
---------------

To debug reference leaks in NRT MemInfo, each MemInfo python object has a
``.refcount`` attribute for inspection.  To get the MemInfo from a ndarray
allocated by NRT, use the ``.base`` attribute.

To debug memory leaks in NRT, the ``numba.core.runtime.rtsys`` defines
``.get_allocation_stats()``.  It returns a namedtuple containing the
number of allocation and deallocation since the start of the program.
Checking that the allocation and deallocation counters are matching is the
simplest way to know if the NRT is leaking.


Debugging Leaks in C
--------------------

The start of `numba/core/runtime/nrt.h
<https://github.com/numba/numba/blob/master/numba/core/runtime/nrt.h>`_
has these lines:

.. code-block:: C

  /* Debugging facilities - enabled at compile-time */
  /* #undef NDEBUG */
  #if 0
  #   define NRT_Debug(X) X
  #else
  #   define NRT_Debug(X) if (0) { X; }
  #endif

Undefining NDEBUG (uncomment the ``#undef NDEBUG`` line) enables the assertion
check in NRT.

Enabling the NRT_Debug (replace ``#if 0`` with ``#if 1``) turns on
debug print inside NRT.


Recursion Support
=================

During the compilation of a pair of mutually recursive functions, one of the
functions will contain unresolved symbol references since the compiler handles
one function at a time.  The memory for the unresolved symbols is allocated and
initialized to the address of the *unresolved symbol abort* function
(``nrt_unresolved_abort``) just before the machine code is
generated by LLVM.  These symbols are tracked and resolved as new functions are
compiled.  If a bug prevents the resolution of these symbols,
the abort function will be called, raising a ``RuntimeError`` exception.

The *unresolved symbol abort* function is defined in the NRT with a zero-argument
signature. The caller is safe to call it with arbitrary number of
arguments.  Therefore, it is safe to be used inplace of the intended callee.

Using the NRT from C code
=========================

Externally compiled C code should use the ``NRT_api_functions`` struct as a
function table to access the NRT API. The struct is defined in
:ghfile:`numba/core/runtime/nrt_external.h`. Users can use the utility function
``numba.extending.include_path()`` to determine the include directory for
Numba provided C headers.

.. literalinclude:: ../../../numba/core/runtime/nrt_external.h
  :language: C
  :caption: `numba/core/runtime/nrt_external.h`

Inside Numba compiled code, the ``numba.core.unsafe.nrt.NRT_get_api()``
intrinsic can be used to obtain a pointer to the ``NRT_api_functions``.

Here is an example that uses the ``nrt_external.h``:

.. code-block:: C

  #include <stdio.h>
  #include "numba/core/runtime/nrt_external.h"

  void my_dtor(void *ptr) {
      free(ptr);
  }

  NRT_MemInfo* my_allocate(NRT_api_functions *nrt) {
      /* heap allocate some memory */
      void * data = malloc(10);
      /* wrap the allocated memory; yield a new reference */
      NRT_MemInfo *mi = nrt->manage_memory(data, my_dtor);
      /* acquire reference */
      nrt->acquire(mi);
      /* release reference */
      nrt->release(mi);
      return mi;
  }

It is important to ensure that the NRT is initialized prior to making calls to
it, calling ``numba.core.runtime.nrt.rtsys.initialize(context)`` from Python
will have the desired effect. Similarly the code snippet:

.. code-block:: Python

  from numba.core.registry import cpu_target # Get the CPU target singleton
  cpu_target.target_context # Access the target_context property to initialize

will achieve the same specifically for Numba's CPU target (the default). Failure
to initialize the NRT will result in access violations as function pointers for
various internal atomic operations will be missing in the ``NRT_MemSys`` struct.

Future Plan
===========

The plan for NRT is to make a standalone shared library that can be linked to
Numba compiled code, including use within the Python interpreter and without
the Python interpreter.  To make that work, we will be doing some refactoring:

* numba :term:`NPM` code references statically compiled code in "helperlib.c".
  Those functions should be moved to NRT.
.. _arch-pipeline:

========================
Customizing the Compiler
========================

.. warning:: The custom pipeline feature is for expert use only.  Modifying
             the compiler behavior can invalidate internal assumptions in the
             numba source code.


For library developers looking for a way to extend or modify the compiler
behavior, you can do so by defining a custom compiler by inheriting from
``numba.compiler.CompilerBase``.  The default Numba compiler is defined
as ``numba.compiler.Compiler``, implementing the ``.define_pipelines()``
method, which adds the *nopython-mode*, *object-mode* and *interpreted-mode*
pipelines. For convenience these three pipelines are defined in
``numba.compiler.DefaultPassBuilder`` by the methods:

* ``.define_nopython_pipeline()``
* ``.define_objectmode_pipeline()``
* ``.define_interpreted_pipeline()``

respectively.

To use a custom subclass of ``CompilerBase``, supply it as the
``pipeline_class`` keyword argument to the ``@jit`` and ``@generated_jit``
decorators.  By doing so, the effect of the custom pipeline is limited to the
function being decorated.

Implementing a compiler pass
----------------------------

Numba makes it possible to implement a new compiler pass and does so through the
use of an API similar to that of LLVM. The following demonstrates the basic
process involved.


Compiler pass classes
#####################

All passes must inherit from ``numba.compiler_machinery.CompilerPass``, commonly
used subclasses are:

* ``numba.compiler_machinery.FunctionPass`` for describing a pass that operates
  on a function-at-once level and may mutate the IR state.
* ``numba.compiler_machinery.AnalysisPass`` for describing a pass that performs
  analysis only.
* ``numba.compiler_machinery.LoweringPass`` for describing a pass that performs
  lowering only.

In this example a new compiler pass will be implemented that will rewrite all
``ir.Const(x)`` nodes, where ``x`` is a subclass of ``numbers.Number``, such
that the value of x is incremented by one. There is no use for this pass other
than to serve as a pedagogical vehicle!

The ``numba.compiler_machinery.FunctionPass`` is appropriate for the suggested
pass behavior and so is the base class of the new pass. Further, a ``run_pass``
method is defined to do the work (this method is abstract, all compiler passes
must implement it).

First the new class:

.. literalinclude:: compiler_pass_example.py
   :language: python
   :dedent: 4
   :start-after: magictoken.ex_compiler_pass.begin
   :end-before: magictoken.ex_compiler_pass.end


Note also that the class must be registered with Numba's compiler machinery
using ``@register_pass``. This in part is to allow the declaration of whether
the pass mutates the control flow graph and whether it is an analysis only pass.

Next, define a new compiler based on the existing
``numba.compiler.CompilerBase``. The compiler pipeline is defined through the
use of an existing pipeline and the new pass declared above is added to be run
after the ``IRProcessing`` pass.


.. literalinclude:: compiler_pass_example.py
   :language: python
   :dedent: 4
   :start-after: magictoken.ex_compiler_defn.begin
   :end-before: magictoken.ex_compiler_defn.end

Finally update the ``@njit`` decorator at the call site to make use of the newly
defined compilation pipeline.

.. literalinclude:: compiler_pass_example.py
   :language: python
   :dedent: 4
   :start-after: magictoken.ex_compiler_call.begin
   :end-before: magictoken.ex_compiler_call.end

Debugging compiler passes
-------------------------

Observing IR Changes
####################

It is often useful to be able to see the changes a pass makes to the IR. Numba
conveniently permits this through the use of the environment variable
:envvar:`NUMBA_DEBUG_PRINT_AFTER`. In the case of the above pass, running the
example code with ``NUMBA_DEBUG_PRINT_AFTER="ir_processing,consts_add_one"``
gives:


.. code-block:: none
    :emphasize-lines: 4, 7, 24, 27

    ----------------------------nopython: ir_processing-----------------------------
    label 0:
        x = arg(0, name=x)                       ['x']
        $const0.1 = const(int, 10)               ['$const0.1']
        a = $const0.1                            ['$const0.1', 'a']
        del $const0.1                            []
        $const0.2 = const(float, 20.2)           ['$const0.2']
        b = $const0.2                            ['$const0.2', 'b']
        del $const0.2                            []
        $0.5 = x + a                             ['$0.5', 'a', 'x']
        del x                                    []
        del a                                    []
        $0.7 = $0.5 + b                          ['$0.5', '$0.7', 'b']
        del b                                    []
        del $0.5                                 []
        c = $0.7                                 ['$0.7', 'c']
        del $0.7                                 []
        $0.9 = cast(value=c)                     ['$0.9', 'c']
        del c                                    []
        return $0.9                              ['$0.9']
    ----------------------------nopython: consts_add_one----------------------------
    label 0:
        x = arg(0, name=x)                       ['x']
        $const0.1 = const(int, 11)               ['$const0.1']
        a = $const0.1                            ['$const0.1', 'a']
        del $const0.1                            []
        $const0.2 = const(float, 21.2)           ['$const0.2']
        b = $const0.2                            ['$const0.2', 'b']
        del $const0.2                            []
        $0.5 = x + a                             ['$0.5', 'a', 'x']
        del x                                    []
        del a                                    []
        $0.7 = $0.5 + b                          ['$0.5', '$0.7', 'b']
        del b                                    []
        del $0.5                                 []
        c = $0.7                                 ['$0.7', 'c']
        del $0.7                                 []
        $0.9 = cast(value=c)                     ['$0.9', 'c']
        del c                                    []
        return $0.9                              ['$0.9']

Note the change in the values in the ``const`` nodes.

Pass execution times
####################

Numba has built-in support for timing all compiler passes, the execution times
are stored in the metadata associated with a compilation result. This
demonstrates one way of accessing this information based on the previously
defined function, ``foo``:

.. literalinclude:: compiler_pass_example.py
   :language: python
   :dedent: 4
   :start-after: magictoken.ex_compiler_timings.begin
   :end-before: magictoken.ex_compiler_timings.end

the output of which is, for example::

    pass_timings(init=1.914000677061267e-06, run=4.308700044930447e-05, finalize=1.7400006981915794e-06)

this displaying the pass initialization, run and finalization times in seconds.

================
Notes on Hashing
================

Numba supports the built-in :func:`hash` and does so by simply calling the
:func:`__hash__` member function on the supplied argument. This makes it
trivial to add hash support for new types as all that is required is the
application of the extension API :func:`overload_method` decorator to overload
a function for computing the hash value for the new type registered to the
type's :func:`__hash__` method. For example::

    from numba.extending import overload_method

    @overload_method(myType, '__hash__')
    def myType_hash_overload(obj):
        # implementation details


The Implementation
==================

The implementation of the Numba hashing functions strictly follows that of
Python 3. The only exception to this is that for hashing Unicode and bytes (for
content longer than ``sys.hash_info.cutoff``) the only supported algorithm is
``siphash24`` (default in CPython 3). As a result Numba will match Python 3
hash values for all supported types under the default conditions described.

Unicode hash cache differences
------------------------------

Both Numba and CPython Unicode string internal representations have a ``hash``
member for the purposes of caching the string's hash value. This member is
always checked ahead of computing a hash value with the view of simply providing
a value from cache as it is considerably cheaper to do so. The Numba Unicode
string hash caching implementation behaves in a similar way to that of
CPython's. The only notable behavioral change (and its only impact is a minor
potential change in performance) is that Numba always computes and caches the
hash for Unicode strings created in ``nopython mode`` at the time they are boxed
for reuse in Python, this is too eager in some cases in comparison to CPython
which may delay hashing a new Unicode string depending on creation method. It
should also be noted that Numba copies in the ``hash`` member of the CPython
internal representation for Unicode strings when unboxing them to its own
representation so as to not recompute the hash of a string that already has a
hash value associated with it.

The accommodation of ``PYTHONHASHSEED``
---------------------------------------

The ``PYTHONHASHSEED`` environment variable can be used to seed the CPython
hashing algorithms for e.g. the purposes of reproduciblity. The Numba hashing
implementation directly reads the CPython hashing algorithms' internal state and
as a result the influence of ``PYTHONHASHSEED`` is replicated in Numba's
hashing implementations.
.. _live variable analysis:

======================
Live Variable Analysis
======================

(Related issue https://github.com/numba/numba/pull/1611)

Numba uses reference-counting for garbage collection, a technique that
requires cooperation by the compiler.  The Numba IR encodes the location
where a decref must be inserted.  These locations are determined by live
variable analysis.  The corresponding source code is the ``_insert_var_dels()``
method in https://github.com/numba/numba/blob/master/numba/interpreter.py.


In Python semantic, once a variable is defined inside a function, it is alive
until the variable is explicitly deleted or the function scope is ended.
However, Numba analyzes the code to determine the minimum bound of the lifetime
of each variable by its definition and usages during compilation.
As soon as a variable is unreachable, a ``del`` instruction is inserted at the
closest basic-block (either at the start of the next block(s) or at the
end of the current block).  This means variables can be released earlier than in
regular Python code.

The behavior of the live variable analysis affects memory usage of the compiled
code.  Internally, Numba does not differentiate temporary variables and user
variables.  Since each operation generates at least one temporary variable,
a function can accumulate a high number of temporary variables if they are not
released as soon as possible.
Our generator implementation can benefit from early releasing of variables,
which reduces the size of the state to suspend at each yield point.


Notes on behavior of the live variable analysis
================================================


Variable deleted before definition
-----------------------------------

(Related issue: https://github.com/numba/numba/pull/1738)

When a variable lifetime is confined within the loop body (its definition and
usage does not escape the loop body), like:

.. code-block:: python

    def f(arr):
      # BB 0
      res = 0
      # BB 1
      for i in (0, 1):
          # BB 2
          t = arr[i]
          if t[i] > 1:
              # BB 3
              res += t[i]
      # BB 4
      return res


Variable ``t`` is never referenced outside of the loop.
A ``del`` instruction is emitted for ``t`` at the head of the loop (BB 1)
before a variable is defined.  The reason is obvious once we know the control
flow graph::

             +------------------------------> BB4
             |
             |
    BB 0 --> BB 1  -->  BB 2 ---> BB 3
             ^          |          |
             |          V          V
             +---------------------+


Variable ``t`` is defined in BB 1.  In BB 2, the evaluation of
``t[i] > 1`` uses ``t``, which is the last use if execution takes the false
branch and goto BB 1.  In BB 3, ``t`` is only used in ``res += t[i]``, which is
the last use if execution takes the true branch.  Because BB 3, an outgoing
branch of BB 2 uses ``t``, ``t`` must be deleted at the common predecessor.
The closest point is BB 1, which does not have ``t`` defined from the incoming
edge of BB 0.

Alternatively, if ``t`` is deleted at BB 4, we will still have to delete the
variable before its definition because BB4 can be executed without executing
the loop body (BB 2 and BB 3), where the variable is defined.

.. _architecture:

==================
Numba architecture
==================

Introduction
============

Numba is a compiler for Python bytecode with optional type-specialization.

Suppose you enter a function like this into the standard Python interpreter
(henceforward referred to as "CPython")::

    def add(a, b):
        return a + b

The interpreter will immediately parse the function and convert it into a
bytecode representation that describes how the CPython interpreter should
execute the function at a low level.  For the example above, it looks
something like this::

    >>> import dis
    >>> dis.dis(add)
    2           0 LOAD_FAST                0 (a)
                3 LOAD_FAST                1 (b)
                6 BINARY_ADD
                7 RETURN_VALUE


CPython uses a stack-based interpreter (much like an HP calculator), so the
code first pushes two local variables onto the stack.  The ``BINARY_ADD``
opcode pops the top two arguments off the stack and makes a Python C API
function call that is equivalent to calling ``a.__add__(b)``.  The result is
then pushed onto the top of the interpreter stack.  Finally, the
``RETURN_VALUE`` opcode returns value on the top of the stack as the result of
the function call.

Numba can take this bytecode and compile it to machine code that performs the
same operations as the CPython interpreter, treating ``a`` and ``b`` as
generic Python objects.  The full semantics of Python are preserved, and the
compiled function can be used with any kind of objects that have the add
operator defined.  When a Numba function is compiled this way, we say that it
has been compiled in :term:`object mode`, because the code still manipulates
Python objects.

Numba code compiled in object mode is not much faster than executing the
original Python function in the CPython interpreter.  However, if we
specialize the function to only run with certain data types, Numba can
generate much shorter and more efficient code that manipulates the data
natively without any calls into the Python C API.  When code has been compiled
for specific data types so that the function body no longer relies on the
Python runtime, we say the function has been compiled in :term:`nopython mode`.
Numeric code compiled in nopython mode can be hundreds of times faster
than the original Python.


Compiler architecture
=====================

Like many compilers, Numba can be conceptually divided into a
*frontend* and a *backend*.

The Numba *frontend* comprises the stages which analyze the Python bytecode,
translate it to :term:`Numba IR` and perform various transformations and
analysis steps on the IR.  One of the key steps is :term:`type inference`.
The frontend must succeed in typing all variables unambiguously in order
for the backend to generate code in :term:`nopython mode`, because the
backend uses type information to match appropriate code generators with
the values they operate on.

The Numba *backend* walks the Numba IR resulting from the frontend analyses
and exploits the type information deduced by the type inference phase to
produce the right LLVM code for each encountered operation.  After LLVM
code is produced, the LLVM library is asked to optimize it and generate
native processor code for the final, native function.

There are other pieces besides the compiler frontend and backend, such
as the caching machinery for JIT functions.  Those pieces are not considered
in this document.


Contexts
========

Numba is quite flexible, allowing it to generate code for different hardware
architectures like CPUs and GPUs.  In order to support these different
applications, Numba uses a *typing context* and a *target context*.

A *typing context* is used in the compiler frontend to perform type inference
on operations and values in the function.  Similar typing contexts could be
used for many architectures because for nearly all cases, typing inference
is hardware-independent.  However, Numba currently has a different typing
context for each target.

A *target context* is used to generate the specific instruction sequence
required to operate on the Numba types identified during type inference.
Target contexts are architecture-specific and are flexible in defining
the execution model and available Python APIs.  For example, Numba has a "cpu"
and a "cuda" context for those two kinds of architecture, and a "parallel"
context which produces multithreaded CPU code.


Compiler stages
===============

The :func:`~numba.jit` decorator in Numba ultimately calls
``numba.compiler.compile_extra()`` which compiles the Python function in a
multi-stage process, described below.

Stage 1: Analyze bytecode
-------------------------

At the start of compilation, the function bytecode is passed to an instance of
the Numba interpreter (``numba.interpreter``).  The interpreter object
analyzes the bytecode to find the control flow graph (``numba.controlflow``).
The control flow graph (CFG) describes the ways that execution can move from one
block to the next inside the function as a result of loops and branches.

The data flow analysis (``numba.dataflow``) takes the control flow graph and
traces how values get pushed and popped off the Python interpreter stack for
different code paths.  This is important to understand the lifetimes of
variables on the stack, which are needed in Stage 2.

If you set the environment variable ``NUMBA_DUMP_CFG`` to 1, Numba will dump
the results of the control flow graph analysis to the screen.  Our ``add()``
example is pretty boring, since there is only one statement block::

    CFG adjacency lists:
    {0: []}
    CFG dominators:
    {0: set([0])}
    CFG post-dominators:
    {0: set([0])}
    CFG back edges: []
    CFG loops:
    {}
    CFG node-to-loops:
    {0: []}

A function with more complex flow control will have a more interesting
control flow graph.  This function::

    def doloops(n):
        acc = 0
        for i in range(n):
            acc += 1
            if n == 10:
                break
        return acc

compiles to this bytecode::

      9           0 LOAD_CONST               1 (0)
                  3 STORE_FAST               1 (acc)

     10           6 SETUP_LOOP              46 (to 55)
                  9 LOAD_GLOBAL              0 (range)
                 12 LOAD_FAST                0 (n)
                 15 CALL_FUNCTION            1
                 18 GET_ITER
            >>   19 FOR_ITER                32 (to 54)
                 22 STORE_FAST               2 (i)

     11          25 LOAD_FAST                1 (acc)
                 28 LOAD_CONST               2 (1)
                 31 INPLACE_ADD
                 32 STORE_FAST               1 (acc)

     12          35 LOAD_FAST                0 (n)
                 38 LOAD_CONST               3 (10)
                 41 COMPARE_OP               2 (==)
                 44 POP_JUMP_IF_FALSE       19

     13          47 BREAK_LOOP
                 48 JUMP_ABSOLUTE           19
                 51 JUMP_ABSOLUTE           19
            >>   54 POP_BLOCK

     14     >>   55 LOAD_FAST                1 (acc)
                 58 RETURN_VALUE

The corresponding CFG for this bytecode is::

    CFG adjacency lists:
    {0: [6], 6: [19], 19: [54, 22], 22: [19, 47], 47: [55], 54: [55], 55: []}
    CFG dominators:
    {0: set([0]),
     6: set([0, 6]),
     19: set([0, 6, 19]),
     22: set([0, 6, 19, 22]),
     47: set([0, 6, 19, 22, 47]),
     54: set([0, 6, 19, 54]),
     55: set([0, 6, 19, 55])}
    CFG post-dominators:
    {0: set([0, 6, 19, 55]),
     6: set([6, 19, 55]),
     19: set([19, 55]),
     22: set([22, 55]),
     47: set([47, 55]),
     54: set([54, 55]),
     55: set([55])}
    CFG back edges: [(22, 19)]
    CFG loops:
    {19: Loop(entries=set([6]), exits=set([54, 47]), header=19, body=set([19, 22]))}
    CFG node-to-loops:
    {0: [], 6: [], 19: [19], 22: [19], 47: [], 54: [], 55: []}

The numbers in the CFG refer to the bytecode offsets shown just to the left
of the opcode names above.

.. _arch_generate_numba_ir:

Stage 2: Generate the Numba IR
------------------------------

Once the control flow and data analyses are complete, the Numba interpreter
can step through the bytecode and translate it into an Numba-internal
intermediate representation.  This translation process changes the function
from a stack machine representation (used by the Python interpreter) to a
register machine representation (used by LLVM).

Although the IR is stored in memory as a tree of objects, it can be serialized
to a string for debugging.  If you set the environment variable
``NUMBA_DUMP_IR`` equal to 1, the Numba IR will be dumped to the screen.  For
the ``add()`` function described above, the Numba IR looks like::

   label 0:
       a = arg(0, name=a)                       ['a']
       b = arg(1, name=b)                       ['b']
       $0.3 = a + b                             ['$0.3', 'a', 'b']
       del b                                    []
       del a                                    []
       $0.4 = cast(value=$0.3)                  ['$0.3', '$0.4']
       del $0.3                                 []
       return $0.4                              ['$0.4']

The ``del`` instructions are produced by :ref:`live variable analysis`.
Those instructions ensure references are not leaked.
In :term:`nopython mode`, some objects are tracked by the Numba runtime and
some are not.  For tracked objects, a dereference operation is emitted;
otherwise, the instruction is an no-op.
In :term:`object mode` each variable contains an owned reference to a PyObject.


.. _`rewrite-untyped-ir`:

Stage 3: Rewrite untyped IR
---------------------------

Before running type inference, it may be desired to run certain
transformations on the Numba IR.  One such example is to detect ``raise``
statements which have an implicitly constant argument, so as to
support them in :term:`nopython mode`.  Let's say you compile the
following function with Numba::

   def f(x):
      if x == 0:
         raise ValueError("x cannot be zero")

If you set the :envvar:`NUMBA_DUMP_IR` environment variable to ``1``,
you'll see the IR being rewritten before the type inference phase::

   REWRITING:
       del $0.3                                 []
       $12.1 = global(ValueError: <class 'ValueError'>) ['$12.1']
       $const12.2 = const(str, x cannot be zero) ['$const12.2']
       $12.3 = call $12.1($const12.2)           ['$12.1', '$12.3', '$const12.2']
       del $const12.2                           []
       del $12.1                                []
       raise $12.3                              ['$12.3']
   ____________________________________________________________
       del $0.3                                 []
       $12.1 = global(ValueError: <class 'ValueError'>) ['$12.1']
       $const12.2 = const(str, x cannot be zero) ['$const12.2']
       $12.3 = call $12.1($const12.2)           ['$12.1', '$12.3', '$const12.2']
       del $const12.2                           []
       del $12.1                                []
       raise <class 'ValueError'>('x cannot be zero') []


.. _arch_type_inference:

Stage 4: Infer types
--------------------

Now that the Numba IR has been generated, type analysis can be performed.  The
types of the function arguments can be taken either from the explicit function
signature given in the ``@jit`` decorator (such as ``@jit('float64(float64,
float64)')``), or they can be taken from the types of the actual function
arguments if compilation is happening when the function is first called.

The type inference engine is found in ``numba.typeinfer``.  Its job is to
assign a type to every intermediate variable in the Numba IR.  The result of
this pass can be seen by setting the :envvar:`NUMBA_DUMP_ANNOTATION`
environment variable to 1:

.. code-block:: python

   -----------------------------------ANNOTATION-----------------------------------
   # File: archex.py
   # --- LINE 4 ---

   @jit(nopython=True)

   # --- LINE 5 ---

   def add(a, b):

       # --- LINE 6 ---
       # label 0
       #   a = arg(0, name=a)  :: int64
       #   b = arg(1, name=b)  :: int64
       #   $0.3 = a + b  :: int64
       #   del b
       #   del a
       #   $0.4 = cast(value=$0.3)  :: int64
       #   del $0.3
       #   return $0.4

       return a + b


If type inference fails to find a consistent type assignment for all the
intermediate variables, it will label every variable as type ``pyobject`` and
fall back to object mode.  Type inference can fail when unsupported Python
types, language features, or functions are used in the function body.


.. _`rewrite-typed-ir`:

Stage 5a: Rewrite typed IR
--------------------------

This pass's purpose is to perform any high-level optimizations that still
require, or could at least benefit from, Numba IR type information.

One example of a problem domain that isn't as easily optimized once
lowered is the domain of multidimensional array operations.  When
Numba lowers an array operation, Numba treats the operation like a
full ufunc kernel.  During lowering a single array operation, Numba
generates an inline broadcasting loop that creates a new result array.
Then Numba generates an application loop that applies the operator
over the array inputs.  Recognizing and rewriting these loops once
they are lowered into LLVM is hard, if not impossible.

An example pair of optimizations in the domain of array operators is
loop fusion and shortcut deforestation.  When the optimizer
recognizes that the output of one array operator is being fed into
another array operator, and only to that array operator, it can fuse
the two loops into a single loop.  The optimizer can further eliminate
the temporary array allocated for the initial operation by directly
feeding the result of the first operation into the second, skipping
the store and load to the intermediate array.  This elimination is
known as shortcut deforestation.  Numba currently uses the rewrite
pass to implement these array optimizations.  For more information,
please consult the ":ref:`case-study-array-expressions`" subsection,
later in this document.

One can see the result of rewriting by setting the
:envvar:`NUMBA_DUMP_IR` environment variable to a non-zero value (such
as 1).  The following example shows the output of the rewrite pass as
it recognizes an array expression consisting of a multiply and add,
and outputs a fused kernel as a special operator, :func:`arrayexpr`::

  ______________________________________________________________________
  REWRITING:
  a0 = arg(0, name=a0)                     ['a0']
  a1 = arg(1, name=a1)                     ['a1']
  a2 = arg(2, name=a2)                     ['a2']
  $0.3 = a0 * a1                           ['$0.3', 'a0', 'a1']
  del a1                                   []
  del a0                                   []
  $0.5 = $0.3 + a2                         ['$0.3', '$0.5', 'a2']
  del a2                                   []
  del $0.3                                 []
  $0.6 = cast(value=$0.5)                  ['$0.5', '$0.6']
  del $0.5                                 []
  return $0.6                              ['$0.6']
  ____________________________________________________________
  a0 = arg(0, name=a0)                     ['a0']
  a1 = arg(1, name=a1)                     ['a1']
  a2 = arg(2, name=a2)                     ['a2']
  $0.5 = arrayexpr(ty=array(float64, 1d, C), expr=('+', [('*', [Var(a0, test.py (14)), Var(a1, test.py (14))]), Var(a2, test.py (14))])) ['$0.5', 'a0', 'a1', 'a2']
  del a0                                   []
  del a1                                   []
  del a2                                   []
  $0.6 = cast(value=$0.5)                  ['$0.5', '$0.6']
  del $0.5                                 []
  return $0.6                              ['$0.6']
  ______________________________________________________________________

Following this rewrite, Numba lowers the array expression into a new
ufunc-like function that is inlined into a single loop that only
allocates a single result array.


.. _`parallel-accelerator`:

Stage 5b: Perform Automatic Parallelization
-------------------------------------------

This pass is only performed if the ``parallel`` option in the :func:`~numba.jit`
decorator is set to ``True``.  This pass finds parallelism implicit in the
semantics of operations in the Numba IR and replaces those operations
with explicitly parallel representations of those operations using a
special `parfor` operator.  Then, optimizations are performed to maximize
the number of parfors that are adjacent to each other such that they can
then be fused together into one parfor that takes only one pass over the
data and will thus typically have better cache performance.  Finally,
during lowering, these parfor operators are converted to a form similar
to guvectorize to implement the actual parallelism.

The automatic parallelization pass has a number of sub-passes, many of
which are controllable using a dictionary of options passed via the
``parallel`` keyword argument to :func:`~numba.jit`::

   { 'comprehension': True/False,  # parallel comprehension
     'prange':        True/False,  # parallel for-loop
     'numpy':         True/False,  # parallel numpy calls
     'reduction':     True/False,  # parallel reduce calls
     'setitem':       True/False,  # parallel setitem
     'stencil':       True/False,  # parallel stencils
     'fusion':        True/False,  # enable fusion or not
   }

The default is set to `True` for all of them. The sub-passes are
described in more detail in the following paragraphs.

#. CFG Simplification
    Sometimes Numba IR will contain chains of blocks containing no loops which
    are merged in this sub-pass into single blocks.  This sub-pass simplifies
    subsequent analysis of the IR.

#. Numpy canonicalization
    Some Numpy operations can be written as operations on Numpy objects (e.g.
    ``arr.sum()``), or as calls to Numpy taking those objects (e.g.
    ``numpy.sum(arr)``).  This sub-pass converts all such operations to the
    latter form for cleaner subsequent analysis.

#. Array analysis
    A critical requirement for later parfor fusion is that parfors have
    identical iteration spaces and these iteration spaces typically correspond
    to the sizes of the dimensions of Numpy arrays.  In this sub-pass, the IR is
    analyzed to determine equivalence classes for the dimensions of Numpy
    arrays.  Consider the example, ``a = b + 1``, where ``a`` and ``b`` are both
    Numpy arrays.  Here, we know that each dimension of ``a`` must have the same
    equivalence class as the corresponding dimension of ``b``.  Typically,
    routines rich in Numpy operations will enable equivalence classes to be
    fully known for all arrays created within a function.

    Array analysis will also reason about size equivalence for slice selection,
    and boolean array masking (one dimensional only). For example, it is able to
    infer that ``a[1 : n-1]`` is of the same size as ``b[0 : n-2]``.

    Array analysis may also insert safety assumptions to ensure pre-conditions
    related to array sizes are met before an operation can be parallelized.
    For example, ``np.dot(X, w)`` between a 2-D matrix ``X`` and a 1-D vector ``w``
    requires that the second dimension of ``X`` is of the same size as ``w``.
    Usually this kind of runtime check is automatically inserted, but if array
    analysis can infer such equivalence, it will skip them.

    Users can even help array analysis by turning implicit knowledge about
    array sizes into explicit assertions. For example, in the code below:

    .. code-block:: python

       @numba.njit(parallel=True)
       def logistic_regression(Y, X, w, iterations):
           assert(X.shape == (Y.shape[0], w.shape[0]))
           for i in range(iterations):
               w -= np.dot(((1.0 / (1.0 + np.exp(-Y * np.dot(X, w))) - 1.0) * Y), X)
           return w

    Making the explicit assertion helps eliminate all bounds checks in the
    rest of the function.

#. ``prange()`` to parfor
    The use of prange (:ref:`numba-prange`) in a for loop is an explicit
    indication from the programmer that all iterations of the for loop can
    execute in parallel.  In this sub-pass, we analyze the CFG to locate loops
    and to convert those loops controlled by a prange object to the explicit
    `parfor` operator.  Each explicit parfor operator consists of:

    a. A list of loop nest information that describes the iteration space of the
       parfor.  Each entry in the loop nest list contains an indexing variable,
       the start of the range, the end of the range, and the step value for each
       iteration.
    #. An initialization (init) block which contains instructions to be executed
       one time before the parfor begins executing.
    #. A loop body comprising a set of basic blocks that correspond to the body
       of the loop and compute one point in the iteration space.
    #. The index variables used for each dimension of the iteration space.

    For parfor `pranges`, the loop nest is a single entry where the start,
    stop, and step fields come from the specified `prange`.  The init block is
    empty for `prange` parfors and the loop body is the set of blocks in the
    loop minus the loop header.

    With parallelization on, array comprehensions (:ref:`pysupported-comprehension`)
    will also be translated to prange so as to run in parallel. This behavior
    be disabled by setting ``parallel={'comprehension': False}``.

    Likewise, the overall `prange` to `parfor` translation can be disabled by
    setting ``parallel={'prange': False}``, in which case `prange` is treated the
    same as `range`.

#. Numpy to parfor
    In this sub-pass, Numpy functions such as ``ones``, ``zeros``, ``dot``, most
    of the random number generating functions, arrayexprs (from Section
    :ref:`rewrite-typed-ir`), and Numpy reductions are converted to parfors.
    Generally, this conversion creates the loop nest list, whose length is equal
    to the number of dimensions of the left-hand side of the assignment
    instruction in the IR.  The number and size of the dimensions of the
    left-hand-side array is taken from the array analysis information generated
    in sub-pass 3 above.  An instruction to create the result Numpy array is
    generated and stored in the new parfor's init block.  A basic block is
    created for the loop body and an instruction is generated and added to the
    end of that block to store the result of the computation into the array at
    the current point in the iteration space.  The result stored into the array
    depends on the operation that is being converted.  For example, for ``ones``,
    the value stored is a constant 1.  For calls to generate a random array, the
    value comes from a call to the same random number function but with the size
    parameter dropped and therefore returning a scalar.  For arrayexpr operators,
    the arrayexpr tree is converted to Numba IR and the value at the root of that
    expression tree is used to write into the output array. The translation from
    Numpy functions and arrayexpr operators to `parfor` can be disabled by
    setting ``parallel={'numpy': False}``.

    For reductions, the loop nest list is similarly created using the array
    analysis information for the array being reduced.  In the init block, the
    initial value is assigned to the reduction variable.  The loop body consists
    of a single block in which the next value in the iteration space is fetched
    and the reduction operation is applied to that value and the current
    reduction value and the result stored back into the reduction value.
    The translation of reduction functions to `parfor` can be disabled by
    setting ``parallel={'reduction': False}``.

    Setting the :envvar:`NUMBA_DEBUG_ARRAY_OPT_STATS` environment variable to
    1 will show some statistics about parfor conversions in general.

#. Setitem to parfor
    Setting a range of array elements using a slice or boolean array selection
    can also run in parallel.  Statement such as ``A[P] = B[Q]``
    (or a simpler case ``A[P] = c``, where ``c`` is a scalar) is translated to
    `parfor` if one of the following conditions is met:

     a. ``P`` and ``Q`` are slices or multi-dimensional selector involving
        scalar and slices, and ``A[P]`` and ``B[Q]`` are considered size
        equivalent by array analysis. Only 2-value slice/range is supported,
        3-value with a step will not be translated to `parfor`.
     #. ``P`` and ``Q`` are the same boolean array.

    This translation can be disabled by setting ``parallel={'setitem': False}``.

#. Simplification
    Performs a copy propagation and dead code elimination pass.

#. Fusion
    This sub-pass first processes each basic block and does a reordering of the
    instructions within the block with the goal of pushing parfors lower in the
    block and lifting non-parfors towards the start of the block.  In practice,
    this approach does a good job of getting parfors adjacent to each other in
    the IR, which enables more parfors to then be fused.  During parfor fusion,
    each basic block is repeatedly scanned until no further fusion is possible.
    During this scan, each set of adjacent instructions are considered.
    Adjacent instructions are fused together if:

    a. they are both parfors
    #. the parfors' loop nests are the same size and the array equivalence
       classes for each dimension of the loop nests are the same, and
    #. the first parfor does not create a reduction variable used by the
       second parfor.

    The two parfors are fused together by adding the second parfor's init block
    to the first's, merging the two parfors' loop bodies together and replacing
    the instances of the second parfor's loop index variables in the second
    parfor's body with the loop index variables for the first parfor.
    Fusion can be disabled by setting ``parallel={'fusion': False}``.

    Setting the :envvar:`NUMBA_DEBUG_ARRAY_OPT_STATS` environment variable to
    1 will show some statistics about parfor fusions.

#. Push call objects and compute parfor parameters
    In the lowering phase described in Section :ref:`lowering`, each parfor
    becomes a separate function executed in parallel in ``guvectorize``
    (:ref:`guvectorize`) style.  Since parfors may use variables defined
    previously in a function, when those parfors become separate functions,
    those variables must be passed to the parfor function as parameters.  In
    this sub-pass, a use-def scan is made over each parfor body and liveness
    information is used to determine which variables are used but not defined by
    the parfor.  That list of variables is stored here in the parfor for use
    during lowering.  Function variables are a special case in this process
    since function variables cannot be passed to functions compiled in nopython
    mode.  Instead, for function variables, this sub-pass pushes the assignment
    instruction to the function variable into the parfor body so that those do
    not need to be passed as parameters.

    To see the intermediate IR between the above sub-passes and other debugging
    information, set the :envvar:`NUMBA_DEBUG_ARRAY_OPT` environment variable to
    1. For the example in Section :ref:`rewrite-typed-ir`, the following IR with
    a parfor is generated during this stage::

     ______________________________________________________________________
     label 0:
         a0 = arg(0, name=a0)                     ['a0']
         a0_sh_attr0.0 = getattr(attr=shape, value=a0) ['a0', 'a0_sh_attr0.0']
         $consta00.1 = const(int, 0)              ['$consta00.1']
         a0size0.2 = static_getitem(value=a0_sh_attr0.0, index_var=$consta00.1, index=0) ['$consta00.1', 'a0_sh_attr0.0', 'a0size0.2']
         a1 = arg(1, name=a1)                     ['a1']
         a1_sh_attr0.3 = getattr(attr=shape, value=a1) ['a1', 'a1_sh_attr0.3']
         $consta10.4 = const(int, 0)              ['$consta10.4']
         a1size0.5 = static_getitem(value=a1_sh_attr0.3, index_var=$consta10.4, index=0) ['$consta10.4', 'a1_sh_attr0.3', 'a1size0.5']
         a2 = arg(2, name=a2)                     ['a2']
         a2_sh_attr0.6 = getattr(attr=shape, value=a2) ['a2', 'a2_sh_attr0.6']
         $consta20.7 = const(int, 0)              ['$consta20.7']
         a2size0.8 = static_getitem(value=a2_sh_attr0.6, index_var=$consta20.7, index=0) ['$consta20.7', 'a2_sh_attr0.6', 'a2size0.8']
     ---begin parfor 0---
     index_var =  parfor_index.9
     LoopNest(index_variable=parfor_index.9, range=0,a0size0.2,1 correlation=5)
     init block:
         $np_g_var.10 = global(np: <module 'numpy' from '/usr/local/lib/python3.5/dist-packages/numpy/__init__.py'>) ['$np_g_var.10']
         $empty_attr_attr.11 = getattr(attr=empty, value=$np_g_var.10) ['$empty_attr_attr.11', '$np_g_var.10']
         $np_typ_var.12 = getattr(attr=float64, value=$np_g_var.10) ['$np_g_var.10', '$np_typ_var.12']
         $0.5 = call $empty_attr_attr.11(a0size0.2, $np_typ_var.12, kws=(), func=$empty_attr_attr.11, vararg=None, args=[Var(a0size0.2, test2.py (7)), Var($np_typ_var.12, test2.py (7))]) ['$0.5', '$empty_attr_attr.11', '$np_typ_var.12', 'a0size0.2']
     label 1:
         $arg_out_var.15 = getitem(value=a0, index=parfor_index.9) ['$arg_out_var.15', 'a0', 'parfor_index.9']
         $arg_out_var.16 = getitem(value=a1, index=parfor_index.9) ['$arg_out_var.16', 'a1', 'parfor_index.9']
         $arg_out_var.14 = $arg_out_var.15 * $arg_out_var.16 ['$arg_out_var.14', '$arg_out_var.15', '$arg_out_var.16']
         $arg_out_var.17 = getitem(value=a2, index=parfor_index.9) ['$arg_out_var.17', 'a2', 'parfor_index.9']
         $expr_out_var.13 = $arg_out_var.14 + $arg_out_var.17 ['$arg_out_var.14', '$arg_out_var.17', '$expr_out_var.13']
         $0.5[parfor_index.9] = $expr_out_var.13  ['$0.5', '$expr_out_var.13', 'parfor_index.9']
     ----end parfor 0----
         $0.6 = cast(value=$0.5)                  ['$0.5', '$0.6']
         return $0.6                              ['$0.6']
     ______________________________________________________________________

  .. _`lowering`:

Stage 6a: Generate nopython LLVM IR
-----------------------------------

If type inference succeeds in finding a Numba type for every intermediate
variable, then Numba can (potentially) generate specialized native code.  This
process is called :term:`lowering`.  The Numba IR tree is translated into
LLVM IR by using helper classes from `llvmlite <http://llvmlite.pydata.org/>`_.
The machine-generated LLVM IR can seem unnecessarily verbose, but the LLVM
toolchain is able to optimize it quite easily into compact, efficient code.

The basic lowering algorithm is generic, but the specifics of how particular
Numba IR nodes are translated to LLVM instructions is handled by the
target context selected for compilation.  The default target context is
the "cpu" context, defined in ``numba.targets.cpu``.

The LLVM IR can be displayed by setting the :envvar:`NUMBA_DUMP_LLVM` environment
variable to 1.  For the "cpu" context, our ``add()`` example would look like:

.. code-block:: llvm

   define i32 @"__main__.add$1.int64.int64"(i64* %"retptr",
                                            {i8*, i32}** %"excinfo",
                                            i8* %"env",
                                            i64 %"arg.a", i64 %"arg.b")
   {
      entry:
        %"a" = alloca i64
        %"b" = alloca i64
        %"$0.3" = alloca i64
        %"$0.4" = alloca i64
        br label %"B0"
      B0:
        store i64 %"arg.a", i64* %"a"
        store i64 %"arg.b", i64* %"b"
        %".8" = load i64* %"a"
        %".9" = load i64* %"b"
        %".10" = add i64 %".8", %".9"
        store i64 %".10", i64* %"$0.3"
        %".12" = load i64* %"$0.3"
        store i64 %".12", i64* %"$0.4"
        %".14" = load i64* %"$0.4"
        store i64 %".14", i64* %"retptr"
        ret i32 0
   }

The post-optimization LLVM IR can be output by setting
:envvar:`NUMBA_DUMP_OPTIMIZED` to 1.  The optimizer shortens the code
generated above quite significantly:

.. code-block:: llvm

   define i32 @"__main__.add$1.int64.int64"(i64* nocapture %retptr,
                                            { i8*, i32 }** nocapture readnone %excinfo,
                                            i8* nocapture readnone %env,
                                            i64 %arg.a, i64 %arg.b)
   {
      entry:
        %.10 = add i64 %arg.b, %arg.a
        store i64 %.10, i64* %retptr, align 8
        ret i32 0
   }

If created during :ref:`parallel-accelerator`, parfor operations are
lowered in the following manner.  First, instructions in the parfor's init
block are lowered into the existing function using the normal lowering code.
Second, the loop body of the parfor is turned into a separate GUFunc.
Third, code is emitted for the current function to call the parallel GUFunc.

To create a GUFunc from the parfor body, the signature of the GUFunc is
created by taking the parfor parameters as identified in step 9 of
Stage :ref:`parallel-accelerator` and adding to that a special `schedule`
parameter, across which the GUFunc will be parallelized.  The schedule
parameter is in effect a static schedule mapping portions of the parfor
iteration space to Numba threads and so the length of the schedule
array is the same as the number of configured Numba threads.  To make
this process easier and somewhat less dependent on changes to Numba IR,
this stage creates a Python function as text that contains the parameters
to the GUFunc and iteration code that takes the current schedule entry
and loops through the specified portion of the iteration space.  In the
body of that loop, a special sentinel is inserted for subsequent easy
location.  This code that handles the processing of the iteration space
is then ``eval``'ed into existence and the Numba compiler's run_frontend
function is called to generate IR.  That IR is scanned to locate the
sentinel and the sentinel is replaced with the loop body of the parfor.
Then, the process of creating the parallel GUFunc is completed by
compiling this merged IR with the Numba compiler's ``compile_ir`` function.

To call the parallel GUFunc, the static schedule must be created.
Code is inserted to call a function named ``do_scheduling.``  This function
is called with the size of each of the parfor's dimensions and the number
`N` of configured Numba threads (:envvar:`NUMBA_NUM_THREADS`).
The ``do_scheduling`` function will divide
the iteration space into N approximately equal sized regions (linear for
1D, rectangular for 2D, or hyperrectangles for 3+D) and the resulting
schedule is passed to the parallel GUFunc.  The number of threads
dedicated to a given dimension of the full iteration space is roughly
proportional to the ratio of the size of the given dimension to the sum
of the sizes of all the dimensions of the iteration space.

Parallel reductions are not natively provided by GUFuncs but the parfor
lowering strategy allows us to use GUFuncs in a way that reductions can
be performed in parallel.  To accomplish this, for each reduction variable
computed by a parfor, the parallel GUFunc and the code that calls it are
modified to make the scalar reduction variable into an array of reduction
variables whose length is equal to the number of Numba threads.  In addition,
the GUFunc still contains a scalar version of the reduction variable that
is updated by the parfor body during each iteration.  One time at the
end of the GUFunc this local reduction variable is copied into the
reduction array.  In this way, false sharing of the reduction array is
prevented.  Code is also inserted into the main
function after the parallel GUFunc has returned that does a reduction
across this smaller reduction array and this final reduction value is
then stored into the original scalar reduction variable.

The GUFunc corresponding to the example from Section :ref:`parallel-accelerator`
can be seen below::

  ______________________________________________________________________
  label 0:
      sched.29 = arg(0, name=sched)            ['sched.29']
      a0 = arg(1, name=a0)                     ['a0']
      a1 = arg(2, name=a1)                     ['a1']
      a2 = arg(3, name=a2)                     ['a2']
      _0_5 = arg(4, name=_0_5)                 ['_0_5']
      $3.1.24 = global(range: <class 'range'>) ['$3.1.24']
      $const3.3.21 = const(int, 0)             ['$const3.3.21']
      $3.4.23 = getitem(value=sched.29, index=$const3.3.21) ['$3.4.23', '$const3.3.21', 'sched.29']
      $const3.6.28 = const(int, 1)             ['$const3.6.28']
      $3.7.27 = getitem(value=sched.29, index=$const3.6.28) ['$3.7.27', '$const3.6.28', 'sched.29']
      $const3.8.32 = const(int, 1)             ['$const3.8.32']
      $3.9.31 = $3.7.27 + $const3.8.32         ['$3.7.27', '$3.9.31', '$const3.8.32']
      $3.10.36 = call $3.1.24($3.4.23, $3.9.31, kws=[], func=$3.1.24, vararg=None, args=[Var($3.4.23, <string> (2)), Var($3.9.31, <string> (2))]) ['$3.1.24', '$3.10.36', '$3.4.23', '$3.9.31']
      $3.11.30 = getiter(value=$3.10.36)       ['$3.10.36', '$3.11.30']
      jump 1                                   []
  label 1:
      $28.2.35 = iternext(value=$3.11.30)      ['$28.2.35', '$3.11.30']
      $28.3.25 = pair_first(value=$28.2.35)    ['$28.2.35', '$28.3.25']
      $28.4.40 = pair_second(value=$28.2.35)   ['$28.2.35', '$28.4.40']
      branch $28.4.40, 2, 3                    ['$28.4.40']
  label 2:
      $arg_out_var.15 = getitem(value=a0, index=$28.3.25) ['$28.3.25', '$arg_out_var.15', 'a0']
      $arg_out_var.16 = getitem(value=a1, index=$28.3.25) ['$28.3.25', '$arg_out_var.16', 'a1']
      $arg_out_var.14 = $arg_out_var.15 * $arg_out_var.16 ['$arg_out_var.14', '$arg_out_var.15', '$arg_out_var.16']
      $arg_out_var.17 = getitem(value=a2, index=$28.3.25) ['$28.3.25', '$arg_out_var.17', 'a2']
      $expr_out_var.13 = $arg_out_var.14 + $arg_out_var.17 ['$arg_out_var.14', '$arg_out_var.17', '$expr_out_var.13']
      _0_5[$28.3.25] = $expr_out_var.13        ['$28.3.25', '$expr_out_var.13', '_0_5']
      jump 1                                   []
  label 3:
      $const44.1.33 = const(NoneType, None)    ['$const44.1.33']
      $44.2.39 = cast(value=$const44.1.33)     ['$44.2.39', '$const44.1.33']
      return $44.2.39                          ['$44.2.39']
  ______________________________________________________________________


Stage 6b: Generate object mode LLVM IR
--------------------------------------

If type inference fails to find Numba types for all values inside a function,
the function will be compiled in object mode.  The generated LLVM will be
significantly longer, as the compiled code will need to make calls to the
`Python C API <https://docs.python.org/3/c-api/>`_ to perform basically all
operations.  The optimized LLVM for our example ``add()`` function is:

.. code-block:: llvm

   @PyExc_SystemError = external global i8
   @".const.Numba_internal_error:_object_mode_function_called_without_an_environment" = internal constant [73 x i8] c"Numba internal error: object mode function called without an environment\00"
   @".const.name_'a'_is_not_defined" = internal constant [24 x i8] c"name 'a' is not defined\00"
   @PyExc_NameError = external global i8
   @".const.name_'b'_is_not_defined" = internal constant [24 x i8] c"name 'b' is not defined\00"

   define i32 @"__main__.add$1.pyobject.pyobject"(i8** nocapture %retptr, { i8*, i32 }** nocapture readnone %excinfo, i8* readnone %env, i8* %arg.a, i8* %arg.b) {
   entry:
     %.6 = icmp eq i8* %env, null
     br i1 %.6, label %entry.if, label %entry.endif, !prof !0

   entry.if:                                         ; preds = %entry
     tail call void @PyErr_SetString(i8* @PyExc_SystemError, i8* getelementptr inbounds ([73 x i8]* @".const.Numba_internal_error:_object_mode_function_called_without_an_environment", i64 0, i64 0))
     ret i32 -1

   entry.endif:                                      ; preds = %entry
     tail call void @Py_IncRef(i8* %arg.a)
     tail call void @Py_IncRef(i8* %arg.b)
     %.21 = icmp eq i8* %arg.a, null
     br i1 %.21, label %B0.if, label %B0.endif, !prof !0

   B0.if:                                            ; preds = %entry.endif
     tail call void @PyErr_SetString(i8* @PyExc_NameError, i8* getelementptr inbounds ([24 x i8]* @".const.name_'a'_is_not_defined", i64 0, i64 0))
     tail call void @Py_DecRef(i8* null)
     tail call void @Py_DecRef(i8* %arg.b)
     ret i32 -1

   B0.endif:                                         ; preds = %entry.endif
     %.30 = icmp eq i8* %arg.b, null
     br i1 %.30, label %B0.endif1, label %B0.endif1.1, !prof !0

   B0.endif1:                                        ; preds = %B0.endif
     tail call void @PyErr_SetString(i8* @PyExc_NameError, i8* getelementptr inbounds ([24 x i8]* @".const.name_'b'_is_not_defined", i64 0, i64 0))
     tail call void @Py_DecRef(i8* %arg.a)
     tail call void @Py_DecRef(i8* null)
     ret i32 -1

   B0.endif1.1:                                      ; preds = %B0.endif
     %.38 = tail call i8* @PyNumber_Add(i8* %arg.a, i8* %arg.b)
     %.39 = icmp eq i8* %.38, null
     br i1 %.39, label %B0.endif1.1.if, label %B0.endif1.1.endif, !prof !0

   B0.endif1.1.if:                                   ; preds = %B0.endif1.1
     tail call void @Py_DecRef(i8* %arg.a)
     tail call void @Py_DecRef(i8* %arg.b)
     ret i32 -1

   B0.endif1.1.endif:                                ; preds = %B0.endif1.1
     tail call void @Py_DecRef(i8* %arg.b)
     tail call void @Py_DecRef(i8* %arg.a)
     tail call void @Py_IncRef(i8* %.38)
     tail call void @Py_DecRef(i8* %.38)
     store i8* %.38, i8** %retptr, align 8
     ret i32 0
   }

   declare void @PyErr_SetString(i8*, i8*)

   declare void @Py_IncRef(i8*)

   declare void @Py_DecRef(i8*)

   declare i8* @PyNumber_Add(i8*, i8*)


The careful reader might notice several unnecessary calls to ``Py_IncRef``
and ``Py_DecRef`` in the generated code.  Currently Numba isn't able to
optimize those away.

Object mode compilation will also attempt to identify loops which can be
extracted and statically-typed for "nopython" compilation.  This process is
called *loop-lifting*, and results in the creation of a hidden nopython mode
function just containing the loop which is then called from the original
function.  Loop-lifting helps improve the performance of functions that
need to access uncompilable code (such as I/O or plotting code) but still
contain a time-intensive section of compilable code.

Stage 7: Compile LLVM IR to machine code
----------------------------------------

In both :term:`object mode` and :term:`nopython mode`, the generated LLVM IR
is compiled by the LLVM JIT compiler and the machine code is loaded into
memory.  A Python wrapper is also created (defined in
``numba.dispatcher.Dispatcher``) which can do the dynamic dispatch to the
correct version of the compiled function if multiple type specializations
were generated (for example, for both ``float32`` and ``float64`` versions
of the same function).

The machine assembly code generated by LLVM can be dumped to the screen by
setting the :envvar:`NUMBA_DUMP_ASSEMBLY` environment variable to 1:

.. code-block:: gas

           .globl  __main__.add$1.int64.int64
           .align  16, 0x90
           .type   __main__.add$1.int64.int64,@function
   __main__.add$1.int64.int64:
           addq    %r8, %rcx
           movq    %rcx, (%rdi)
           xorl    %eax, %eax
           retq

The assembly output will also include the generated wrapper function that
translates the Python arguments to native data types.
Event API
=========

.. automodule:: numba.core.event
    :members:==========================
Notes on Target Extensions
==========================

.. warning:: All features and APIs described in this page are in-development and
             may change at any time without deprecation notices being issued.


Inheriting compiler flags from the caller
=========================================

Compiler flags, i.e. options such as ``fastmath``, ``nrt`` in
``@jit(nrt=True, fastmath=True))`` are specified per-function but their
effects are not well-defined---some flags affect the entire callgraph, some
flags affect only the current function. Sometimes it is necessary for callees
to inherit flags from the caller; for example the ``fastmath`` flag should be
infectious.

To address the problem, the following are needed:

1. Better definitions for the semantics of compiler flags. Preferably, all flags should
   limit their effect to the current function. (TODO)
2. Allow compiler flags to be inherited from the caller. (Done)
3. Consider compiler flags in function resolution. (TODO)

:class:`numba.core.targetconfig.ConfigStack` is used to propagate the compiler flags
throughout the compiler. At the start of the compilation, the flags are pushed
into the ``ConfigStack``, which maintains a thread-local stack for the
compilation. Thus, callees can check the flags in the caller.

.. autoclass:: numba.core.targetconfig.ConfigStack
    :members:

Compiler flags
--------------

`Compiler flags`_ are defined as a subclass of ``TargetConfig``:

.. _Compiler flags: https://github.com/numba/numba/blob/7e8538140ce3f8d01a5273a39233b5481d8b20b1/numba/core/compiler.py#L39

.. autoclass:: numba.core.targetconfig.TargetConfig
    :members:


These are internal compiler flags and they are different from the user-facing
options used in the jit decorators.

Internally, `the user-facing options are mapped to the internal compiler flags <https://github.com/numba/numba/blob/7e8538140ce3f8d01a5273a39233b5481d8b20b1/numba/core/options.py#L72>`_
by :class:`numba.core.options.TargetOptions`. Each target can override the
default compiler flags and control the flag inheritance in
``TargetOptions.finalize``. `The CPU target overrides it.
<https://github.com/numba/numba/blob/7e8538140ce3f8d01a5273a39233b5481d8b20b1/numba/core/cpu.py#L259>`_

.. autoclass:: numba.core.options.TargetOptions
    :members: finalize


In :meth:`numba.core.options.TargetOptions.finalize`,
use :meth:`numba.core.targetconfig.TargetConfig.inherit_if_not_set`
to request a compiler flag from the caller if it is not set for the current
function.
.. _developer-caching:

================
Notes on Caching
================

Numba supports caching of compiled functions into the filesystem for future
use of the same functions.


The Implementation
==================

Caching is done by saving the compiled *object code*, the ELF object of the
executable code.  By using the *object code*, cached functions have minimal
overhead because no compilation is needed. The cached data is saved under the
cache directory (see :envvar:`NUMBA_CACHE_DIR`). The index of the cache is
stored in a ``.nbi`` file, with one index per function, and it lists all the
overloaded signatures compiled for the function. The *object code* is stored in
files with an ``.nbc`` extension, one file per overload. The data in both files
is serialized with :mod:`pickle`.

.. note:: On Python <=3.7, Numba extends ``pickle`` using the pure-Python
          pickler. To use the faster C Pickler, install ``pickle5``
          from ``pip``. ``pickle5`` backports Python 3.8 pickler features.


Requirements for Cacheability
-----------------------------

Developers should note the requirements of a function to permit it to be cached
to ensure that the features they are working on are compatible with caching.

Requirements for cacheable function:

- The LLVM module must be *self-contained*, meaning that it cannot rely on
  other compiled units without linking to them.
- The only allowed external symbols are from the
  :ref:`NRT <arch-numba-runtime>` or other common symbols from system libraries
  (i.e. libc and libm).

Debugging note:

- Look for the usage of ``inttoptr`` in the LLVM IR or
  ``target_context.add_dynamic_add()`` in the lowering code in Python.
  They indicate potential usage of runtime address. Not all uses are
  problematic and some are necessary. Only the conversion of constant integers
  into pointers will affect caching.
- Misuse of dynamic address or dynamic symbols will likely result in a
  segfault.
- Linking order matters because unused symbols are dropped after linking.
  Linking should start from the leaf nodes of the dependency graph.


Features Compatible with Caching
--------------------------------

The following features are explicitly verified to work with caching.

- ufuncs and gufuncs for the ``cpu`` and ``parallel`` target
- parallel accelerator features (i.e. ``parallel=True``)


Caching Limitations
-------------------

This is a list of known limitation of the cache:

- Cache invalidation fails to recognize changes in symbols defined in a
  different file.
- Global variables are treated as constants. The cache will remember the value
  in the global variable used at compilation. On cache load, the cached
  function will not rebind to the new value of the global variable.


.. _cache-sharing:

Cache Sharing
-------------

It is safe to share and reuse the contents in the cache directory on a
different machine. The cache remembers the CPU model and the available
CPU features during compilation. If the CPU model and the CPU features do
not match exactly, the cache contents will not be considered.
(Also see :envvar:`NUMBA_CPU_NAME`)

If the cache directory is shared on a network filesystem, concurrent
read/write of the cache is safe only if file replacement operation is atomic
for the filesystem. Numba always writes to a unique temporary file first, it
then replaces the target cache file path with the temporary file. Numba is
tolerant against lost cache files and lost cache entries.

.. _cache-clearing:

Cache Clearing
--------------

The cache is invalidated when the corresponding source file is modified.
However, it is necessary sometimes to clear the cache directory manually.
For instance, changes in the compiler will not be recognized because the source
files are not modified.

To clear the cache, the cache directory can be simply removed.

Removing the cache directory when a Numba application is running may cause an
``OSError`` exception to be raised at the compilation site.

Related Environment Variables
-----------------------------

See :ref:`env-vars for caching <numba-envvars-caching>`.

==================
Environment Object
==================

The Environment object (Env) is used to maintain references to python objects
that are needed to support compiled functions for both object-mode and
nopython-mode.

In nopython-mode, the Env is used for:

* Storing pyobjects for reconstruction from native values,
  such as:
  
  * for printing native values of NumPy arrays;
  * for returning or yielding native values back to the interpreter.

In object-mode, the Env is used for:

* storing constant values referenced in the code.
* storing a reference to the function's global dictionary to load global
  values.


The Implementation
==================

The Env is implemented in two parts.  In ``_dynfunc.c``, the Env is defined
as ``EnvironmentObject`` as a Python C-extension type.  In ``lowering.py``,
the `EnvironmentObject`` (exported as ``_dynfunc.Environment``) is extended
to support necessary operations needed at lowering.


Serialization
-------------

The Env supports being pickled.  Compilation cache files and ahead-of-time
compiled modules serialize all the used Envs for recreation at the runtime.

Usage
-----

At the start of the lowering for a function or a generator, an Env is created.
Throughout the compilation, the Env is mutated to attach additional
information.  The compiled code references an Env via a global variable in
the emitted LLVM IR.  The global variable is zero-initialized with "common"
linkage, which is the default linkage for C global values.  The use of this
linkage allows multiple definitions of the global variable to be merged into
a single definition when the modules are linked together.  The name of the
global variable is computed from the name of the function
(see ``FunctionDescriptor.env_name`` and ``.get_env_name()`` of the target
context).

The Env is initialized when the compiled-function is loaded. The JIT engine
finds the address of the associated global variable for the Env and stores the
address of the Env into it.  For cached functions, the same process applies.
For ahead-of-time compiled functions, the module initializer in the generated
library is responsible for initializing the global variables of all the Envs
in the module.
=========================================
Notes on Numba's threading implementation
=========================================

The execution of the work presented by the Numba ``parallel`` targets is
undertaken by the Numba threading layer. Practically, the "threading layer"
is a Numba built-in library that can perform the required concurrent execution.
At the time of writing there are three threading layers available, each
implemented via a different lower level native threading library. More
information on the threading layers and appropriate selection of a threading
layer for a given application/system can be found in the
:ref:`threading layer documentation <numba-threading-layer>`.

The pertinent information to note for the following sections is that the
function in the threading library that performs the parallel execution is the
``parallel_for`` function. The job of this function is to both orchestrate and
execute the parallel tasks.

The relevant source files referenced in this document are

- ``numba/np/ufunc/tbbpool.cpp``
- ``numba/np/ufunc/omppool.cpp``
- ``numba/np/ufunc/workqueue.c``

  These files contain the TBB, OpenMP, and workqueue threadpool
  implementations, respectively. Each includes the functions
  ``set_num_threads()``, ``get_num_threads()``, and ``get_thread_id()``, as
  well as the relevant logic for thread masking in their respective
  schedulers. Note that the basic thread local variable logic is duplicated in
  each of these files, and not shared between them.

- ``numba/np/ufunc/parallel.py``

  This file contains the Python and JIT compatible wrappers for
  ``set_num_threads()``, ``get_num_threads()``, and ``get_thread_id()``, as
  well as the code that loads the above libraries into Python and launches the
  threadpool.

- ``numba/parfors/parfor_lowering.py``

  This file contains the main logic for generating code for the parallel
  backend. The thread mask is accessed in this file in the code that generates
  scheduler code, and passed to the relevant backend scheduler function (see
  below).

Thread masking
--------------

As part of its design, Numba never launches new threads beyond the threads
that are launched initially with ``numba.np.ufunc.parallel._launch_threads()``
when the first parallel execution is run. This is due to the way threads were
already implemented in Numba prior to thread masking being implemented. This
restriction was kept to keep the design simple, although it could be removed
in the future. Consequently, it's possible to programmatically set the number
of threads, but only to less than or equal to the total number that have
already been launched. This is done by "masking" out unused threads, causing
them to do no work. For example, on a 16 core machine, if the user were to
call ``set_num_threads(4)``, Numba would always have 16 threads present, but
12 of them would sit idle for parallel computations. A further call to
``set_num_threads(16)`` would cause those same threads to do work in later
computations.

:ref:`Thread masking <numba-threading-layer-thread-masking>` was added to make
it possible for a user to programmatically alter the number of threads
performing work in the threading layer. Thread masking proved challenging to
implement as it required the development of a programming model that is suitable
for users, easy to reason about, and could be implemented safely, with
consistent behavior across the various threading layers.

Programming model
~~~~~~~~~~~~~~~~~

The programming model chosen is similar to that found in OpenMP. The reasons
for this choice were that it is familiar to a lot of users, restricted in
scope and also simple. The number of threads in use is specified by calling
``set_num_threads`` and the number of threads in use can be queried by calling
``get_num_threads``.These two functions are synonymous with their OpenMP
counterparts (with the above restriction that the mask must be less than or
equal to the number of launched threads). The execution semantics are also
similar to OpenMP in that once a parallel region is launched, altering the
thread mask has no impact on the currently executing region, but will have an
impact on parallel regions executed subsequently.

The Implementation
~~~~~~~~~~~~~~~~~~

So as to place no further restrictions on user code other than those that
already existed in the threading layer libraries, careful consideration of the
design of thread masking was required. The "thread mask" cannot be stored in a
global value as concurrent use of the threading layer may result in classic
forms of race conditions on the value itself. Numerous designs were discussed
involving various types of mutex on such a global value, all of which were
eventually broken through thought experiment alone. It eventually transpired
that, following some OpenMP implementations, the "thread mask" is best
implemented as a ``thread local``. This means each thread that executes a Numba
parallel function will have a thread local storage (TLS) slot that contains the
value of the thread mask to use when scheduling threads in the ``parallel_for``
function.

The above notion of TLS use for a thread mask is relatively easy to implement,
``get_num_threads`` and ``set_num_threads`` simply need to address the TLS slot
in a given threading layer. This also means that the execution schedule for a
parallel region can be derived from a run time call to ``get_num_threads``. This
is achieved via a well known and relatively easy to implement pattern of a ``C``
library function registration and wrapping it in the internal Numba
implementation.

In addition to satisfying the original upfront thread masking requirements, a
few more complicated scenarios needed consideration as follows.

Nested parallelism
******************

In all threading layers a "main thread" will invoke the ``parallel_for``
function and then in the parallel region, depending on the threading layer,
some number of additional threads will assist in doing the actual work.
If the work contains a call to another parallel function (i.e. nested
parallelism) it is necessary for the thread making the call to know what the
"thread mask" of the main thread is so that it can propagate it into the
``parallel_for`` call it makes when executing the nested parallel function.
The implementation of this behavior is threading layer specific but the general
principle is for the "main thread" to always "send" the value of the thread mask
from its TLS slot to all threads in the threading layer that are active in the
parallel region. These active threads then update their TLS slots with this
value prior to performing any work. The net result of this implementation detail
is that:

* thread masks correctly propagate into nested functions
* it's still possible for each thread in a parallel region to safely have a
  different mask with which to call nested functions, if it's not set explicitly
  then the inherited mask from the "main thread" is used
* threading layers which have dynamic scheduling with threads potentially
  joining and leaving the active pool during a ``parallel_for`` execution are
  successfully accommodated
* any "main thread" thread mask is entirely decoupled from the in-flux nature
  of the thread masks of the threads in the active thread pool

Python threads independently invoking parallel functions
********************************************************

The threading layer launch sequence is heavily guarded to ensure that the
launch is both thread and process safe and run once per process. In a system
with numerous Python ``threading`` module threads all using Numba, the first
thread through the launch sequence will get its thread mask set appropriately,
but no further threads can run the launch sequence. This means that other
threads will need their initial thread mask set some other way. This is
achieved when ``get_num_threads`` is called and no thread mask is present, in
this case the thread mask will be set to the default. In the implementation,
"no thread mask is present" is represented by the value ``-1`` and the "default
thread mask" (unset) is represented by the value ``0``. The implementation also
immediately calls ``set_num_threads(NUMBA_NUM_THREADS)`` after doing this, so
if either ``-1`` or ``0`` is encountered as a result from ``get_num_threads()`` it
indicates a bug in the above processes.

OS ``fork()`` calls
*******************

The use of TLS was also in part driven by the Linux (the most popular
platform for Numba use by far) having a ``fork(2, 3P)`` call that will do TLS
propagation into child processes, see ``clone(2)``\ 's ``CLONE_SETTLS``.

Thread ID
*********

A private ``get_thread_id()`` function was added to each threading backend,
which returns a unique ID for each thread. This can be accessed from Python by
``numba.np.ufunc.parallel._get_thread_id()`` (it can also be used inside a
JIT compiled function). The thread ID function is useful for testing that the
thread masking behavior is correct, but it should not be used outside of the
tests. For example, one can call ``set_num_threads(4)`` and then collect all
unique ``_get_thread_id()``\ s in a parallel region to verify that only 4
threads are run.

Caveats
~~~~~~~

Some caveats to be aware of when testing thread masking:

- The TBB backend may choose to schedule fewer than the given mask number of
  threads. Thus a test such as the one described above may return fewer than 4
  unique threads.

- The workqueue backend is not threadsafe, so attempts to do multithreading
  nested parallelism with it may result in deadlocks or other undefined
  behavior. The workqueue backend will raise a SIGABRT signal if it detects
  nested parallelism.

- Certain backends may reuse the main thread for computation, but this
  behavior shouldn't be relied upon (for instance, if propagating exceptions).

Use in Code Generation
~~~~~~~~~~~~~~~~~~~~~~

The general pattern for using ``get_num_threads`` in code generation is

.. code:: python

   import llvmlite.llvmpy.core as lc

   get_num_threads = cgutils.get_or_insert_function(builder.module
       lc.Type.function(lc.Type.int(types.intp.bitwidth), []),
       name="get_num_threads")

   num_threads = builder.call(get_num_threads, [])

   with cgutils.if_unlikely(builder, builder.icmp_signed('<=', num_threads,
                                                 num_threads.type(0))):
       cgutils.printf(builder, "num_threads: %d\n", num_threads)
       context.call_conv.return_user_exc(builder, RuntimeError,
                                                 ("Invalid number of threads. "
                                                  "This likely indicates a bug in Numba.",))

   # Pass num_threads through to the appropriate backend function here

See the code in ``numba/parfors/parfor_lowering.py``.

The guard against ``num_threads`` being <= 0 is not strictly necessary, but it
can protect against accidentally incorrect behavior in case the thread masking
logic contains a bug.

The ``num_threads`` variable should be passed through to the appropriate
backend function, such as ``do_scheduling`` or ``parallel_for``. If it's used
in some way other than passing it through to the backend function, the above
considerations should be taken into account to ensure the use of the
``num_threads`` variable is safe. It would probably be better to keep such
logic in the threading backends, rather than trying to do it in code
generation.
=====================
Numba Project Roadmap
=====================

.. note::
    This page was last revised in *December 2018*.

This roadmap is for informational purposes only.  Priorities and resources
change, so we may choose to reorder or abandon things on this list.
Additionally, the further out items are, the less concrete they will be.  If
you have an interest in working on one of these items, please open an issue
where we can discuss the design and approach first.

Short Term: 2019H1
==================

* Container improvements:

  * Numba dictionary support
  * Refactor lists to follow new container best practices.
    See the discussion in `issue 3546 <https://github.com/numba/numba/issues/3546#issuecomment-443008201>`_.

* Deprecate Python 2.7 support
* Improve caching:

  * Full support for functions compiled with ParallelAccelerator
  * Safe caching of generated functions (eval of strings)
  * Expire cache when any function in call chain (even in other files) changes
  * Process for distributing pre-populated cache

* Continue to improve usability and debugging:

  * Trap more unsupported features earlier in pipeline (especially things that parfors can’t handle)
  * Error messages
  * Diagnostic tools for debugging and understanding performance
  * Better on-boarding for new users and contributors (revise docs, more examples)

* Begin refactoring existing features that cause common bug reports:

  * Enhance description of interfaces provided by Numba functions to give more type information
  * Convert older Numba function implementations to use public extension mechanisms
  * More unit testing and modularization of ParallelAccelerator passes

Medium Term: 2019H2
===================

* Unify dispatch of regular functions, ufuncs, and gufuncs
* Declare Numba 1.0 with stable interfaces
* Continue to improve usability and debugging (see above)
* Continue refactoring Numba internals to solve common bug reports (see above)
* JIT class review and improvement
* Improve compilation speed
* Improve memory management of Numba-allocated memory
* Better support for writing code transformation passes
* Make caching and parallel execution features opt-out instead of opt-in

  * add heuristic to determine if parfor passes will be beneficial

Long Term: 2020 and beyond
==========================

* Unify GPU backends (share more code and interfaces)
* Improve ahead of time compilation (for low powered devices)
* Improve cross language connections (C++, JVM?, Julia?, R?)

  * Call Numba from other languages,
  * Call from Numba into other languages

* Better support for "hybrid" CPU/GPU/TPU/etc programming
* Partial / deferred compilation of functions
* Foster integration of Numba into core PyData packages:

  * scipy/scikit-learn/scikit-image/pandas

* More support for efforts to put Numba into other applications (databases, etc) for compiling user-defined functions
* More support for usage of Numba as a “compiler toolkit” to create custom compilers (like HPAT, automatic differentiation of functions, etc)
* Investigate AST-based Numba frontend in addition to existing bytecode-based frontend
=================
Notes on Inlining
=================

There are occasions where it is useful to be able to inline a function at its
call site, at the Numba IR level of representation. The decorators such as
:func:`numba.jit`, :func:`numba.extending.overload` and
:func:`register_jitable` support the keyword argument ``inline``, to facilitate
this behaviour.

When attempting to inline at this level, it is important to understand what
purpose this serves and what effect this will have. In contrast to the inlining
performed by LLVM, which is aimed at improving performance, the main reason to
inline at the Numba IR level is to allow type inference to cross function
boundaries.

As an example, consider the following snippet:

.. code:: python

    from numba import njit


    @njit
    def bar(a):
        a.append(10)


    @njit
    def foo():
        z = []
        bar(z)


    foo()

This will fail to compile and run, because the type of ``z`` can not be inferred
as it will only be refined within ``bar``. If we now add ``inline=True`` to the
decorator for ``bar`` the snippet will compile and run. This is because inlining
the call to ``a.append(10)`` will mean that ``z`` will be refined to hold integers
and so type inference will succeed.

So, to recap, inlining at the Numba IR level is unlikely to have a performance
benefit. Whereas inlining at the LLVM level stands a better chance.

The ``inline`` keyword argument can be one of three values:

* The string ``'never'``, this is the default and results in the function not
  being inlined under any circumstances.
* The string ``'always'``, this results in the function being inlined at all
  call sites.
* A python function that takes three arguments. The first argument is always the
  ``ir.Expr`` node that is the ``call`` requesting the inline, this is present
  to allow the function to make call contextually aware decisions. The second
  and third arguments are:

  * In the case of an untyped inline, i.e. that which occurs when using the
    :func:`numba.jit` family of decorators, both arguments are
    ``numba.ir.FunctionIR`` instances. The second argument corresponding to the
    IR of the caller, the third argument corresponding to the IR of the callee.

  * In the case of a typed inline, i.e. that which occurs when using
    :func:`numba.extending.overload`, both arguments are instances of a
    ``namedtuple`` with fields (corresponding to their standard use in the
    compiler internals):

    * ``func_ir`` - the function's Numba IR.
    * ``typemap`` - the function's type map.
    * ``calltypes`` - the call types of any calls in the function.
    * ``signature`` - the function's signature.

    The second argument holds the information from the caller, the third holds
    the information from the callee.

  In all cases the function should return True to inline and return False to not
  inline, this essentially permitting custom inlining rules (typical use might
  be cost models).
* Recursive functions with ``inline='always'`` will result in a non-terminating
  compilation. If you wish to avoid this, supply a function to limit the
  recursion depth (see below).

.. note:: No guarantee is made about the order in which functions are assessed
          for inlining or about the order in which they are inlined.


Example using :func:`numba.jit`
===============================

An example of using all three options to ``inline`` in the :func:`numba.njit`
decorator:

.. literalinclude:: inline_example.py

which produces the following when executed (with a print of the IR after the
legalization pass, enabled via the environment variable
``NUMBA_DEBUG_PRINT_AFTER="ir_legalization"``):

.. code-block:: none
    :emphasize-lines: 2, 3, 9, 16, 17, 21, 22, 26, 35

    label 0:
        $0.1 = global(never_inline: CPUDispatcher(<function never_inline at 0x7f890ccf9048>)) ['$0.1']
        $0.2 = call $0.1(func=$0.1, args=[], kws=(), vararg=None) ['$0.1', '$0.2']
        del $0.1                                 []
        a = $0.2                                 ['$0.2', 'a']
        del $0.2                                 []
        $0.3 = global(always_inline: CPUDispatcher(<function always_inline at 0x7f890ccf9598>)) ['$0.3']
        del $0.3                                 []
        $const0.1.0 = const(int, 200)            ['$const0.1.0']
        $0.2.1 = $const0.1.0                     ['$0.2.1', '$const0.1.0']
        del $const0.1.0                          []
        $0.4 = $0.2.1                            ['$0.2.1', '$0.4']
        del $0.2.1                               []
        b = $0.4                                 ['$0.4', 'b']
        del $0.4                                 []
        $0.5 = global(maybe_inline1: CPUDispatcher(<function maybe_inline1 at 0x7f890ccf9ae8>)) ['$0.5']
        $0.6 = call $0.5(func=$0.5, args=[], kws=(), vararg=None) ['$0.5', '$0.6']
        del $0.5                                 []
        d = $0.6                                 ['$0.6', 'd']
        del $0.6                                 []
        $const0.7 = const(int, 13)               ['$const0.7']
        magic_const = $const0.7                  ['$const0.7', 'magic_const']
        del $const0.7                            []
        $0.8 = global(maybe_inline1: CPUDispatcher(<function maybe_inline1 at 0x7f890ccf9ae8>)) ['$0.8']
        del $0.8                                 []
        $const0.1.2 = const(int, 300)            ['$const0.1.2']
        $0.2.3 = $const0.1.2                     ['$0.2.3', '$const0.1.2']
        del $const0.1.2                          []
        $0.9 = $0.2.3                            ['$0.2.3', '$0.9']
        del $0.2.3                               []
        e = $0.9                                 ['$0.9', 'e']
        del $0.9                                 []
        $0.10 = global(maybe_inline2: CPUDispatcher(<function maybe_inline2 at 0x7f890ccf9b70>)) ['$0.10']
        del $0.10                                []
        $const0.1.4 = const(int, 37)             ['$const0.1.4']
        $0.2.5 = $const0.1.4                     ['$0.2.5', '$const0.1.4']
        del $const0.1.4                          []
        $0.11 = $0.2.5                           ['$0.11', '$0.2.5']
        del $0.2.5                               []
        c = $0.11                                ['$0.11', 'c']
        del $0.11                                []
        $0.14 = a + b                            ['$0.14', 'a', 'b']
        del b                                    []
        del a                                    []
        $0.16 = $0.14 + c                        ['$0.14', '$0.16', 'c']
        del c                                    []
        del $0.14                                []
        $0.18 = $0.16 + d                        ['$0.16', '$0.18', 'd']
        del d                                    []
        del $0.16                                []
        $0.20 = $0.18 + e                        ['$0.18', '$0.20', 'e']
        del e                                    []
        del $0.18                                []
        $0.22 = $0.20 + magic_const              ['$0.20', '$0.22', 'magic_const']
        del magic_const                          []
        del $0.20                                []
        $0.23 = cast(value=$0.22)                ['$0.22', '$0.23']
        del $0.22                                []
        return $0.23                             ['$0.23']


Things to note in the above:

1. The call to the function ``never_inline`` remains as a call.
2. The ``always_inline`` function has been inlined, note its
   ``const(int, 200)`` in the caller body.
3. There is a call to ``maybe_inline1`` before the ``const(int, 13)``
   declaration, the cost model prevented this from being inlined.
4. After the ``const(int, 13)`` the subsequent call to ``maybe_inline1`` has
   been inlined as shown by the ``const(int, 300)`` in the caller body.
5. The function ``maybe_inline2`` has been inlined as demonstrated by
   ``const(int, 37)`` in the caller body.
6. That dead code elimination has not been performed and as a result there are
   superfluous statements present in the IR.


Example using :func:`numba.extending.overload`
==============================================

An example of using inlining with the  :func:`numba.extending.overload`
decorator. It is most interesting to note that if a function is supplied as the
argument to ``inline`` a lot more information is available via the supplied
function arguments for use in decision making. Also that different
``@overload`` s can have different inlining behaviours, with multiple ways to
achieve this:

.. literalinclude:: inline_overload_example.py

which produces the following when executed (with a print of the IR after the
legalization pass, enabled via the environment variable
``NUMBA_DEBUG_PRINT_AFTER="ir_legalization"``):

.. code-block:: none
    :emphasize-lines: 2, 3, 4, 5, 6, 15, 16, 17, 18, 19, 20, 21, 22, 28, 29, 30

    label 0:
        $const0.2 = const(tuple, (1, 2, 3))      ['$const0.2']
        x.0 = $const0.2                          ['$const0.2', 'x.0']
        del $const0.2                            []
        $const0.2.2 = const(int, 0)              ['$const0.2.2']
        $0.3.3 = getitem(value=x.0, index=$const0.2.2) ['$0.3.3', '$const0.2.2', 'x.0']
        del x.0                                  []
        del $const0.2.2                          []
        $0.4.4 = $0.3.3                          ['$0.3.3', '$0.4.4']
        del $0.3.3                               []
        $0.3 = $0.4.4                            ['$0.3', '$0.4.4']
        del $0.4.4                               []
        a = $0.3                                 ['$0.3', 'a']
        del $0.3                                 []
        $const0.5 = const(int, 100)              ['$const0.5']
        x.5 = $const0.5                          ['$const0.5', 'x.5']
        del $const0.5                            []
        $const0.2.7 = const(int, 1)              ['$const0.2.7']
        $0.3.8 = x.5 + $const0.2.7               ['$0.3.8', '$const0.2.7', 'x.5']
        del x.5                                  []
        del $const0.2.7                          []
        $0.4.9 = $0.3.8                          ['$0.3.8', '$0.4.9']
        del $0.3.8                               []
        $0.6 = $0.4.9                            ['$0.4.9', '$0.6']
        del $0.4.9                               []
        b = $0.6                                 ['$0.6', 'b']
        del $0.6                                 []
        $0.7 = global(bar: <function bar at 0x7f6c3710d268>) ['$0.7']
        $const0.8 = const(complex, 300j)         ['$const0.8']
        $0.9 = call $0.7($const0.8, func=$0.7, args=[Var($const0.8, inline_overload_example.py (56))], kws=(), vararg=None) ['$0.7', '$0.9', '$const0.8']
        del $const0.8                            []
        del $0.7                                 []
        c = $0.9                                 ['$0.9', 'c']
        del $0.9                                 []
        $0.12 = a + b                            ['$0.12', 'a', 'b']
        del b                                    []
        del a                                    []
        $0.14 = $0.12 + c                        ['$0.12', '$0.14', 'c']
        del c                                    []
        del $0.12                                []
        $0.15 = cast(value=$0.14)                ['$0.14', '$0.15']
        del $0.14                                []
        return $0.15                             ['$0.15']

Things to note in the above:

1. The first highlighted section is the always inlined overload for the
   ``UniTuple`` argument type.
2. The second highlighted section is the overload for the ``Number`` argument
   type that has been inlined as the cost model function decided to do so as the
   argument was an ``Integer`` type instance.
3. The third highlighted section is the overload for the ``Number`` argument
   type that has not inlined as the cost model function decided to reject it as
   the argument was an ``Complex`` type instance.
4. That dead code elimination has not been performed and as a result there are
   superfluous statements present in the IR.

Using a function to limit the inlining depth of a recursive function
====================================================================

When using recursive inlines, you can terminate the compilation by using
a cost model.

.. code:: python

    from numba import njit
    import numpy as np

    class CostModel(object):
        def __init__(self, max_inlines):
            self._count = 0
            self._max_inlines = max_inlines

        def __call__(self, expr, caller, callee):
            ret = self._count < self._max_inlines
            self._count += 1
            return ret

    @njit(inline=CostModel(3))
    def factorial(n):
        if n <= 0:
            return 1
        return n * factorial(n - 1)

    factorial(5)

.. _developer-manual:

Developer Manual
================

.. toctree::
   :maxdepth: 2

   contributing.rst
   release.rst
   repomap.rst
   architecture.rst
   dispatching.rst
   generators.rst
   numba-runtime.rst
   rewrites.rst
   live_variable_analysis.rst
   listings.rst
   stencil.rst
   custom_pipeline.rst
   inlining.rst
   environment.rst
   hashing.rst
   caching.rst
   threading_implementation.rst
   literal.rst
   llvm_timings.rst
   debugging.rst
   event_api.rst
   target_extension.rst
   roadmap.rst
Listings
========

This shows listings from compiler internal registries  (e.g. lowering
definitions).  The information is provided as developer reference.
When possible, links to source code are provided via github links.

New style listings
------------------

The following listings are generated from ``numba.help.inspector.write_listings()``. Users can run ``python -m numba.help.inspector --format=rst <package>`` to recreate the the documentation.

.. toctree::
   :maxdepth: 2

   autogen_builtins_listing.rst
   autogen_math_listing.rst
   autogen_cmath_listing.rst
   autogen_numpy_listing.rst


Old style listings
------------------

.. toctree::
   :maxdepth: 2

   autogen_lower_listing.rst

.. _developer-debugging:

==================
Notes on Debugging
==================

This section describes techniques that can be useful in debugging the
compilation and execution of generated code.

.. seealso::
   :ref:`debugging-jit-compiled-code`


Memcheck
--------

Memcheck_ is a memory error detector implemented using Valgrind_. It is useful
for detecting memory errors in compiled code, particularly out-of-bounds
accesses and use-after-free errors. Buggy or miscompiled native code can
generate these kinds of errors. The `Memcheck documentation
<https://valgrind.org/docs/manual/mc-manual.html>`_ explains its usage; here, we
discuss only the specifics of using it with Numba.

.. _Memcheck: https://valgrind.org/docs/manual/mc-manual.html
.. _Valgrind: https://valgrind.org/

The Python interpreter and some of the libraries used by Numba can generate
false positives with Memcheck - see `this section of the manual
<https://valgrind.org/docs/manual/mc-manual.html#mc-manual.machine>`_ for more
information on why false positives occur. The false positives can make it
difficult to determine when an actual error has occurred, so it is helpful to
suppress known false positives. This can be done by supplying a suppressions
file, which instructs Memcheck to ignore errors that match the suppressions
defined in it.

The CPython source distribution includes a suppressions file, in the file
``Misc/valgrind-python.supp``. Using this file prevents a lot of spurious errors
generated by Python's memory allocation implementation. Additionally, the Numba
repository includes a suppressions file in ``contrib/valgrind-numba.supp``.

.. note:: It is important to use the suppressions files from the versions of the
   Python interpreter and Numba that you are using - these files evolve over
   time, so non-current versions can fail to suppress some errors, or
   erroneously suppress actual errors.

To run the Python interpreter under Memcheck with both suppressions
files, it is invoked with the following command::

   valgrind --tool=memcheck \
            --suppressions=${CPYTHON_SRC_DIR}/Misc/valgrind-python.supp \
            --suppressions=${NUMBA_SRC_DIR}/contrib/valgrind-numba.supp \
            python ${PYTHON_ARGS}

where ``${CPYTHON_SRC_DIR}`` is set to the location of the CPython source
distribution, ``${NUMBA_SRC_DIR}`` is the location of the Numba source dir, and
``${PYTHON_ARGS}`` are the arguments to the Python interpreter.

If there are errors, then messages describing them will be printed to standard
error. An example of an error is::

   ==77113==    at 0x24169A: PyLong_FromLong (longobject.c:251)
   ==77113==    by 0x241881: striter_next (bytesobject.c:3084)
   ==77113==    by 0x2D3C95: _PyEval_EvalFrameDefault (ceval.c:2809)
   ==77113==    by 0x21B499: _PyEval_EvalCodeWithName (ceval.c:3930)
   ==77113==    by 0x26B436: _PyFunction_FastCallKeywords (call.c:433)
   ==77113==    by 0x2D3605: call_function (ceval.c:4616)
   ==77113==    by 0x2D3605: _PyEval_EvalFrameDefault (ceval.c:3124)
   ==77113==    by 0x21B977: _PyEval_EvalCodeWithName (ceval.c:3930)
   ==77113==    by 0x21C2A4: _PyFunction_FastCallDict (call.c:376)
   ==77113==    by 0x2D5129: do_call_core (ceval.c:4645)
   ==77113==    by 0x2D5129: _PyEval_EvalFrameDefault (ceval.c:3191)
   ==77113==    by 0x21B499: _PyEval_EvalCodeWithName (ceval.c:3930)
   ==77113==    by 0x26B436: _PyFunction_FastCallKeywords (call.c:433)
   ==77113==    by 0x2D46DA: call_function (ceval.c:4616)
   ==77113==    by 0x2D46DA: _PyEval_EvalFrameDefault (ceval.c:3139)
   ==77113==
   ==77113== Use of uninitialised value of size 8

The traceback provided only outlines the C call stack, which can make it
difficult to determine what the Python interpreter was doing at the time of the
error. One can learn more about the state of the stack by looking at the
backtrace in the `GNU Debugger (GDB) <https://www.gnu.org/software/gdb/>`_.
Launch ``valgrind`` with an additional argument, ``--vgdb-error=0`` and attach
to the process using GDB as instructed by the output. Once an error is
encountered, GDB will stop at the error and the stack can be inspected.

GDB does provide support for backtracing through the Python stack, but this
requires symbols which may not be easily available in your Python distribution.
In this case, it is still possible to determine some information about what was
happening in Python, but this depends on examining the backtrace closely. For
example, in a backtrace corresponding to the above error, we see items such as:

.. code-block::

   #18 0x00000000002722da in slot_tp_call (
       self=<_wrap_impl(_callable=<_wrap_missing_loc(func=<function at remote
       0x1cf66c20>) at remote 0x1d200bd0>, _imp=<function at remote 0x1d0e7440>,
       _context=<CUDATargetContext(address_size=64,
       typing_context=<CUDATypingContext(_registries={<Registry(functions=[<type
       at remote 0x65be5e0>, <type at remote 0x65be9d0>, <type at remote
       0x65bedc0>, <type at remote 0x65bf1b0>, <type at remote 0x8b78000>, <type
       at remote 0x8b783f0>, <type at remote 0x8b787e0>, <type at remote
       0x8b78bd0>, <type at remote 0x8b78fc0>, <type at remote 0x8b793b0>, <type
       at remote 0x8b797a0>, <type at remote 0x8b79b90>, <type at remote
       0x8b79f80>, <type at remote 0x8b7a370>, <type at remote 0x8b7a760>, <type
       at remote 0x8b7ab50>, <type at remote 0x8b7af40>, <type at remote
       0x8b7b330>, <type at remote 0x8b7b720>, <type at remote 0x8b7bf00>, <type
       at remote 0x8b7c2f0>, <type at remote 0x8b7c6e0>], attributes=[<type at
       remote 0x8b7cad0>, <type at remote 0x8b7cec0>, <type at remote
       0x8b7d2b0>, <type at remote 0x8b7d6a0>, <type at remote 0x8b7da90>,
       <t...(truncated),
       args=(<Builder(_block=<Block(parent=<Function(parent=<Module(context=<Context(scope=<NameScope(_useset={''},
       _basenamemap={}) at remote 0xbb5ae10>, identified_types={}) at remote
       0xbb5add0>, name='cuconstRecAlign$7',
       data_layout='e-p:64:64:64-i1:8:8-i8:8:8-i16:16:16-i32:32:32-i64:64:64-f32:32:32-f64:64:64-v16:16:16-v32:32:32-v64:64:64-v128:128:128-n16:32:64',
       scope=<NameScope(_useset={'',
       '_ZN08NumbaEnv5numba4cuda5tests6cudapy13test_constmem19cuconstRecAlign$247E5ArrayIdLi1E1C7mutable7alignedE5ArrayIdLi1E1C7mutable7alignedE5ArrayIdLi1E1C7mutable7alignedE5ArrayIdLi1E1C7mutable7alignedE5ArrayIdLi1E1C7mutable7alignedE',
       '_ZN5numba4cuda5tests6cudapy13test_constmem19cuconstRecAlign$247E5ArrayIdLi1E1C7mutable7alignedE5ArrayIdLi1E1C7mutable7alignedE5ArrayIdLi1E1C7mutable7alignedE5ArrayIdLi1E1C7mutable7alignedE5ArrayIdLi1E1C7mutable7alignedE'},
       _basenamemap={}) at remote 0x1d27bf10>, triple='nvptx64-nvidia-cuda',
       globals={'_ZN08NumbaEnv5numba4cuda5tests6cudapy13test_constmem19cuconstRecAlign$247E5ArrayIdLi1E1C7mutable7ali...(truncated),
       kwds=0x0)

We can see some of the arguments, in particular the names of the compiled functions, e.g::

   _ZN5numba4cuda5tests6cudapy13test_constmem19cuconstRecAlign$247E5ArrayIdLi1E1C7mutable7alignedE5ArrayIdLi1E1C7mutable7alignedE5ArrayIdLi1E1C7mutable7alignedE5ArrayIdLi1E1C7mutable7alignedE5ArrayIdLi1E1C7mutable7alignedE

We can run this through ``c++filt`` to see a more human-readable representation::

   numba::cuda::tests::cudapy::test_constmem::cuconstRecAlign$247(
     Array<double, 1, C, mutable, aligned>,
     Array<double, 1, C, mutable, aligned>,
     Array<double, 1, C, mutable, aligned>,
     Array<double, 1, C, mutable, aligned>,
     Array<double, 1, C, mutable, aligned>)

which is the fully qualified name of a jitted function and the types with which
it was called.

Numba Release Process
=====================

The goal of the Numba release process -- from a high level perspective -- is to
publish source and binary artifacts that correspond to a given version
number. This usually involves a sequence of individual tasks that must be
performed in the correct order and with diligence. Numba and llvmlite are
commonly released in lockstep since there is usually a one-to-one mapping
between a Numba version and a corresponding llvmlite version.

This section contains various notes and templates that can be used to create a
Numba release checklist on the Numba Github issue tracker. This is an aid for
the maintainers during the release process and helps to ensure that all tasks
are completed in the correct order and that no tasks are accidentally omitted.

If new or additional items do appear during release, please do remember to add
them to the checklist templates. Also note that the release process itself is
always a work in progress. This means that some of the information here may be
outdated. If you notice this please do remember to submit a pull-request to
update this document.

All release checklists are available as Gitub issue templates. To create a new
release checklist simply open a new issue and select the correct template.


Primary Release Candidate Checklist
-----------------------------------

This is for the first/primary release candidate for minor release i.e. the
first release of every series. It is special, because during this release, the
release branch will have to be created. Release candidate indexing begins at 1.

.. literalinclude:: ../../../.github/ISSUE_TEMPLATE/first_rc_checklist.md
    :language: md
    :lines: 9-

`Open a primary release checklist <https://github.com/numba/numba/issues/new?template=first_rc_checklist.md>`_.

Subsequent Release Candidates, Final Releases and Patch Releases
----------------------------------------------------------------

Releases subsequent to the first release in a series usually involves a series
of cherry-picks, the recipe is therefore slightly different.

.. literalinclude:: ../../../.github/ISSUE_TEMPLATE/sub_rc_checklist.md
    :language: md
    :lines: 9-

`Open a subsequent release checklist <https://github.com/numba/numba/issues/new?template=sub_rc_checklist.md>`_.
.. _developer-literally:

======================
Notes on Literal Types
======================

.. note:: This document describes an advanced feature designed to overcome
          some limitations of the compilation mechanism relating to types.

Some features need to specialize based on the literal value during
compliation to produce type stable code necessary for successful compilation in
Numba. This can be achieved by propagating the literal value through the type
system. Numba recognizes inline literal values as :class:`numba.types.Literal`.
For example::

    def foo(x):
        a = 123
        return bar(x, a)

Numba will infer the type of ``a`` as ``Literal[int](123)``. The definition of
``bar()`` can subsequently specialize its implementation knowing that the
second argument is an ``int`` with the value ``123``.

``Literal`` Type
----------------

Classes and methods related to the ``Literal`` type.

.. autoclass:: numba.types.Literal

.. autofunction:: numba.types.literal

.. autofunction:: numba.types.unliteral

.. autofunction:: numba.types.maybe_literal

Specifying for Literal Typing
-----------------------------

To specify a value as a ``Literal`` type in code scheduled for JIT compilation,
use the following function:

.. autofunction:: numba.literally

Code Example
~~~~~~~~~~~~

.. literalinclude:: ../../../numba/tests/doc_examples/test_literally_usage.py
   :language: python
   :caption: from ``test_literally_usage`` of ``numba/tests/doc_examples/test_literally_usage.py``
   :start-after: magictoken.ex_literally_usage.begin
   :end-before: magictoken.ex_literally_usage.end
   :dedent: 4
   :linenos:


Internal Details
~~~~~~~~~~~~~~~~

Internally, the compiler raises a ``ForceLiteralArgs`` exception to signal
the dispatcher to wrap specified arguments using the ``Literal`` type.

.. autoclass:: numba.errors.ForceLiteralArg
    :members: __init__, combine, __or__


Inside Extensions
-----------------

``@overload`` extensions can use ``literally`` inside the implementation body
like in normal jit-code.

Explicit handling of literal requirements is possible through use of the
following:

.. autoclass:: numba.extending.SentryLiteralArgs
    :members:

.. autoclass:: numba.extending.BoundLiteralArgs
    :members:

.. autofunction:: numba.extending.sentry_literal_args
=====================================================
Using the Numba Rewrite Pass for Fun and Optimization
=====================================================

Overview
========

This section introduces intermediate representation (IR) rewrites, and
how they can be used to implement optimizations.

As discussed earlier in ":ref:`rewrite-typed-ir`", rewriting the Numba
IR allows us to perform optimizations that would be much more
difficult to perform at the lower LLVM level.  Similar to the Numba
type and lowering subsystems, the rewrite subsystem is user
extensible.  This extensibility affords Numba the possibility of
supporting a wide variety of domain-specific optimizations (DSO's).

The remaining subsections detail the mechanics of implementing a
rewrite, registering a rewrite with the rewrite registry, and provide
examples of adding new rewrites, as well as internals of the array
expression optimization pass.  We conclude by reviewing some use cases
exposed in the examples, as well as reviewing any points where
developers should take care.


Rewriting Passes
================

Rewriting passes have a simple :func:`~Rewrite.match` and
:func:`~Rewrite.apply` interface.  The division between matching and
rewriting follows how one would define a term rewrite in a declarative
domain-specific languages (DSL's).  In such DSL's, one may write a
rewrite as follows::

  <match> => <replacement>


The ``<match>`` and ``<replacement>`` symbols represent IR term
expressions, where the left-hand side presents a pattern to match, and
the right-hand side an IR term constructor to build upon matching.
Whenever the rewrite matches an IR pattern, any free variables in the
left-hand side are bound within a custom environment.  When applied,
the rewrite uses the pattern matching environment to bind any free
variables in the right-hand side.

As Python is not commonly used in a declarative capacity, Numba uses
object state to handle the transfer of information between the
matching and application steps.


The :class:`Rewrite` Base Class
-------------------------------

.. class:: Rewrite

   The :class:`Rewrite` class simply defines an abstract base class
   for Numba rewrites.  Developers should define rewrites as
   subclasses of this base type, overloading the
   :func:`~Rewrite.match` and :func:`~Rewrite.apply` methods.

   .. attribute:: pipeline

       The pipeline attribute contains the
       :class:`numba.compiler.Pipeline` instance that is currently
       compiling the function under consideration for rewriting.

   .. method:: __init__(self, pipeline, *args, **kws)

       The base constructor for rewrites simply stashes its arguments
       into attributes of the same name.  Unless being used in
       debugging or testing, rewrites should only be constructed by
       the :class:`RewriteRegistry` in the
       :func:`RewriteRegistry.apply` method, and the construction
       interface should remain stable (though the pipeline will
       commonly contain just about everything there is to know).

   .. method:: match(self, block, typemap, callmap)

      The :func:`~Rewrite.match` method takes four arguments other
      than *self*:

      * *func_ir*: This is an instance of :class:`numba.ir.FunctionIR` for the
        function being rewritten.

      * *block*: This is an instance of :class:`numba.ir.Block`.  The
        matching method should iterate over the instructions contained
        in the :attr:`numba.ir.Block.body` member.

      * *typemap*: This is a Python :class:`dict` instance mapping
        from symbol names in the IR, represented as strings, to Numba
        types.

      * *callmap*: This is another :class:`dict` instance mapping from
        calls, represented as :class:`numba.ir.Expr` instances, to
        their corresponding call site type signatures, represented as
        a :class:`numba.typing.templates.Signature` instance.

      The :func:`~Rewrite.match` method should return a :class:`bool`
      result.  A :obj:`True` result should indicate that one or more
      matches were found, and the :func:`~Rewrite.apply` method will
      return a new replacement :class:`numba.ir.Block` instance.  A
      :obj:`False` result should indicate that no matches were found, and
      subsequent calls to :func:`~Rewrite.apply` will return undefined
      or invalid results.

   .. method:: apply(self)

      The :func:`~Rewrite.apply` method should only be invoked
      following a successful call to :func:`~Rewrite.match`.  This
      method takes no additional parameters other than *self*, and
      should return a replacement :class:`numba.ir.Block` instance.

      As mentioned above, the behavior of calling
      :func:`~Rewrite.apply` is undefined unless
      :func:`~Rewrite.match` has already been called and returned
      :obj:`True`.


Subclassing :class:`Rewrite`
----------------------------

Before going into the expectations for the overloaded methods any
:class:`Rewrite` subclass must have, let's step back a minute to
review what is taking place here.  By providing an extensible
compiler, Numba opens itself to user-defined code generators which may
be incomplete, or worse, incorrect.  When a code generator goes awry,
it can cause abnormal program behavior or early termination.
User-defined rewrites add a new level of complexity because they must
not only generate correct code, but the code they generate should
ensure that the compiler does not get stuck in a match/apply loop.
Non-termination by the compiler will directly lead to non-termination
of user function calls.

There are several ways to help ensure that a rewrite terminates:

* *Typing*: A rewrite should generally attempt to decompose composite
  types, and avoid composing new types.  If the rewrite is matching a
  specific type, changing expression types to a lower-level type will
  ensure they will no long match after the rewrite is applied.

* *Special instructions*: A rewrite may synthesize custom operators or
  use special functions in the target IR.  This technique again
  generates code that is no longer within the domain of the original
  match, and the rewrite will terminate.

In the ":ref:`case-study-array-expressions`" subsection, below, we'll
see how the array expression rewriter uses both of these techniques.


Overloading :func:`Rewrite.match`
---------------------------------

Every rewrite developer should seek to have their implementation of
:func:`~Rewrite.match` return a :obj:`False` value as quickly as
possible.  Numba is a just-in-time compiler, and adding compilation
time ultimately adds to the user's run time.  When a rewrite returns
:obj:`False` for a given block, the registry will no longer process that
block with that rewrite, and the compiler is that much closer to
proceeding to lowering.

This need for timeliness has to be balanced against collecting the
necessary information to make a match for a rewrite.  Rewrite
developers should be comfortable adding dynamic attributes to their
subclasses, and then having these new attributes guide construction of
the replacement basic block.


Overloading :func:`Rewrite.apply`
-----------------------------------

The :func:`~Rewrite.apply` method should return a replacement
:class:`numba.ir.Block` instance to replace the basic block that
contained a match for the rewrite.  As mentioned above, the IR built
by :func:`~Rewrite.apply` methods should preserve the semantics of the
user's code, but also seek to avoid generating another match for the
same rewrite or set of rewrites.


The Rewrite Registry
====================

When you want to include a rewrite in the rewrite pass, you should
register it with the rewrite registry.  The :mod:`numba.rewrites`
module provides both the abstract base class and a class decorator for
hooking into the Numba rewrite subsystem.  The following illustrates a
stub definition of a new rewrite::

  from numba import rewrites

  @rewrites.register_rewrite
  class MyRewrite(rewrites.Rewrite):

      def match(self, block, typemap, calltypes):
          raise NotImplementedError("FIXME")

      def apply(self):
          raise NotImplementedError("FIXME")


Developers should note that using the class decorator as shown above
will register a rewrite at import time.  It is the developer's
responsibility to ensure their extensions are loaded before
compilation starts.


.. _`case-study-array-expressions`:

Case study: Array Expressions
=============================

This subsection looks at the array expression rewriter in more depth.
The array expression rewriter, and most of its support functionality,
are found in the :mod:`numba.npyufunc.array_exprs` module.  The
rewriting pass itself is implemented in the :class:`RewriteArrayExprs`
class.  In addition to the rewriter, the
:mod:`~numba.npyufunc.array_exprs` module includes a function for
lowering array expressions,
:func:`~numba.npyufunc.array_exprs._lower_array_expr`.  The overall
optimization process is as follows:

* :func:`RewriteArrayExprs.match`: The rewrite pass looks for two or
  more array operations that form an array expression.

* :func:`RewriteArrayExprs.apply`: Once an array expression is found,
  the rewriter replaces the individual array operations with a new
  kind of IR expression, the ``arrayexpr``.

* :func:`numba.npyufunc.array_exprs._lower_array_expr`: During
  lowering, the code generator calls
  :func:`~numba.npyufunc.array_exprs._lower_array_expr` whenever it
  finds an ``arrayexpr`` IR expression.

More details on each step of the optimization are given below.


The :func:`RewriteArrayExprs.match` method
------------------------------------------

The array expression optimization pass starts by looking for array
operations, including calls to supported :class:`~numpy.ufunc`\'s and
user-defined :class:`~numba.DUFunc`\'s.  Numba IR follows the
conventions of a static single assignment (SSA) language, meaning that
the search for array operators begins with looking for assignment
instructions.

When the rewriting pass calls the :func:`RewriteArrayExprs.match`
method, it first checks to see if it can trivially reject the basic
block.  If the method determines the block to be a candidate for
matching, it sets up the following state variables in the rewrite
object:

* *crnt_block*: The current basic block being matched.

* *typemap*: The *typemap* for the function being matched.

* *matches*: A list of variable names that reference array expressions.

* *array_assigns*: A map from assignment variable names to the actual
  assignment instructions that define the given variable.

* *const_assigns*: A map from assignment variable names to the
  constant valued expression that defines the constant variable.

At this point, the match method iterates over the assignment
instructions in the input basic block.  For each assignment
instruction, the matcher looks for one of two things:

* Array operations: If the right-hand side of the assignment
  instruction is an expression, and the result of that expression is
  an array type, the matcher checks to see if the expression is either
  a known array operation, or a call to a universal function.  If an
  array operator is found, the matcher stores the left-hand variable
  name and the whole instruction in the *array_assigns* member.
  Finally, the matcher tests to see if any operands of the array
  operation have also been identified as targets of other array
  operations.  If one or more operands are also targets of array
  operations, then the matcher will also append the left-hand side
  variable name to the *matches* member.

* Constants: Constants (even scalars) can be operands to array
  operations.  Without worrying about the constant being apart of an
  array expression, the matcher stores constant names and values in
  the *const_assigns* member.

The end of the matching method simply checks for a non-empty *matches*
list, returning :obj:`True` if there were one or more matches, and
:obj:`False` when *matches* is empty.


The :func:`RewriteArrayExprs.apply` method
------------------------------------------

When one or matching array expressions are found by
:func:`RewriteArrayExprs.match`, the rewriting pass will call
:func:`RewriteArrayExprs.apply`.  The apply method works in two
passes.  The first pass iterates over the matches found, and builds a
map from instructions in the old basic block to new instructions in
the new basic block.  The second pass iterates over the instructions
in the old basic block, copying instructions that are not changed by
the rewrite, and replacing or deleting instructions that were
identified by the first pass.

The :func:`RewriteArrayExprs._handle_matches` implements the first
pass of the code generation portion of the rewrite.  For each match,
this method builds a special IR expression that contains an expression
tree for the array expression.  To compute the leaves of the
expression tree, the :func:`~RewriteArrayExprs._handle_matches` method
iterates over the operands of the identified root operation.  If the
operand is another array operation, it is translated into an
expression sub-tree.  If the operand is a constant,
:func:`~RewriteArrayExprs._handle_matches` copies the constant value.
Otherwise, the operand is marked as being used by an array expression.
As the method builds array expression nodes, it builds a map from old
instructions to new instructions (*replace_map*), as well as sets of
variables that may have moved (*used_vars*), and variables that should
be removed altogether (*dead_vars*).  These three data structures are
returned back to the calling :func:`RewriteArrayExprs.apply` method.

The remaining part of the :func:`RewriteArrayExprs.apply` method
iterates over the instructions in the old basic block.  For each
instruction, this method either replaces, deletes, or duplicates that
instruction based on the results of
:func:`RewriteArrayExprs._handle_matches`.  The following list
describes how the optimization handles individual instructions:

* When an instruction is an assignment,
  :func:`~RewriteArrayExprs.apply` checks to see if it is in the
  replacement instruction map.  When an assignment instruction is found
  in the instruction map, :func:`~RewriteArrayExprs.apply` must then
  check to see if the replacement instruction is also in the replacement
  map.  The optimizer continues this check until it either arrives at a
  :obj:`None` value or an instruction that isn't in the replacement map.
  Instructions that have a replacement that is :obj:`None` are deleted.
  Instructions that have a non-:obj:`None` replacement are replaced.
  Assignment instructions not in the replacement map are appended to the
  new basic block with no changes made.

* When the instruction is a delete instruction, the rewrite checks to
  see if it deletes a variable that may still be used by a later array
  expression, or if it deletes a dead variable.  Delete instructions for
  used variables are added to a map of deferred delete instructions that
  :func:`~RewriteArrayExprs.apply` uses to move them past any uses of
  that variable.  The loop copies delete instructions for non-dead
  variables, and ignores delete instructions for dead variables
  (effectively removing them from the basic block).

* All other instructions are appended to the new basic block.

Finally, the :func:`~RewriteArrayExprs.apply` method returns the new
basic block for lowering.


The :func:`~numba.npyufunc.array_exprs._lower_array_expr` function
------------------------------------------------------------------

If we left things at just the rewrite, then the lowering stage of the
compiler would fail, complaining it doesn't know how to lower
``arrayexpr`` operations.  We start by hooking a lowering function
into the target context whenever the :class:`RewriteArrayExprs` class
is instantiated by the compiler.  This hook causes the lowering pass to
call :func:`~numba.npyufunc.array_exprs._lower_array_expr` whenever it
encounters an ``arrayexr`` operator.

This function has two steps:

* Synthesize a Python function that implements the array expression:
  This new Python function essentially behaves like a Numpy
  :class:`~numpy.ufunc`, returning the result of the expression on
  scalar values in the broadcasted array arguments.  The lowering
  function accomplishes this by translating from the array expression
  tree into a Python AST.

* Compile the synthetic Python function into a kernel:  At this point,
  the lowering function relies on existing code for lowering ufunc and
  DUFunc kernels, calling
  :func:`numba.targets.numpyimpl.numpy_ufunc_kernel` after defining
  how to lower calls to the synthetic function.

The end result is similar to loop lifting in Numba's object mode.


Conclusions and Caveats
=======================

We have seen how to implement rewrites in Numba, starting with the
interface, and ending with an actual optimization.  The key points of
this section are:

* When writing a good plug-in, the matcher should try to get a
  go/no-go result as soon as possible.

* The rewrite application portion can be more computationally
  expensive, but should still generate code that won't cause infinite
  loops in the compiler.

* We use object state to communicate any results of matching to the
  rewrite application pass.
.. _developer-llvm-timings:

====================
Notes on timing LLVM
====================


Getting LLVM Pass Timings
-------------------------

The dispatcher stores LLVM pass timings in the dispatcher object metadata under
the ``llvm_pass_timings`` key when :envvar:`NUMBA_LLVM_PASS_TIMINGS` is
enabled or ``numba.config.LLVM_PASS_TIMINGS`` is set to truthy.
The timings information contains details on how much time
has been spent in each pass. The pass timings are also grouped by their purpose.
For example, there will be pass timings for function-level pre-optimizations,
module-level optimizations, and object code generation.


Code Example
~~~~~~~~~~~~

.. literalinclude:: ../../../numba/tests/doc_examples/test_llvm_pass_timings.py
   :language: python
   :caption: from ``test_pass_timings`` of ``numba/tests/doc_examples/test_llvm_pass_timings.py``
   :start-after: magictoken.ex_llvm_pass_timings.begin
   :end-before: magictoken.ex_llvm_pass_timings.end
   :dedent: 16
   :linenos:

Example output:

.. code-block:: text

  Printing pass timings for JITCodeLibrary('DocsLLVMPassTimings.test_pass_timings.<locals>.foo')
  Total time: 0.0376
  == #0 Function passes on '_ZN5numba5tests12doc_examples22test_llvm_pass_timings19DocsLLVMPassTimings17test_pass_timings12$3clocals$3e7foo$241Ex'
  Percent: 4.8%
  Total 0.0018s
  Top timings:
    0.0015s ( 81.6%) SROA #3
    0.0002s (  9.3%) Early CSE #2
    0.0001s (  4.0%) Simplify the CFG #9
    0.0000s (  1.5%) Prune NRT refops #4
    0.0000s (  1.1%) Post-Dominator Tree Construction #5
  == #1 Function passes on '_ZN7cpython5numba5tests12doc_examples22test_llvm_pass_timings19DocsLLVMPassTimings17test_pass_timings12$3clocals$3e7foo$241Ex'
  Percent: 0.8%
  Total 0.0003s
  Top timings:
    0.0001s ( 30.4%) Simplify the CFG #10
    0.0001s ( 24.1%) Early CSE #3
    0.0001s ( 17.8%) SROA #4
    0.0000s (  8.8%) Prune NRT refops #5
    0.0000s (  5.6%) Post-Dominator Tree Construction #6
  == #2 Function passes on 'cfunc._ZN5numba5tests12doc_examples22test_llvm_pass_timings19DocsLLVMPassTimings17test_pass_timings12$3clocals$3e7foo$241Ex'
  Percent: 0.5%
  Total 0.0002s
  Top timings:
    0.0001s ( 27.7%) Early CSE #4
    0.0001s ( 26.8%) Simplify the CFG #11
    0.0000s ( 13.8%) Prune NRT refops #6
    0.0000s (  7.4%) Post-Dominator Tree Construction #7
    0.0000s (  6.7%) Dominator Tree Construction #29
  == #3 Module passes (cheap optimization for refprune)
  Percent: 3.7%
  Total 0.0014s
  Top timings:
    0.0007s ( 52.0%) Combine redundant instructions
    0.0001s (  5.4%) Function Integration/Inlining
    0.0001s (  4.9%) Prune NRT refops #2
    0.0001s (  4.8%) Natural Loop Information
    0.0001s (  4.6%) Post-Dominator Tree Construction #2
  == #4 Module passes (full optimization)
  Percent: 43.9%
  Total 0.0165s
  Top timings:
    0.0032s ( 19.5%) Combine redundant instructions #9
    0.0022s ( 13.5%) Combine redundant instructions #7
    0.0010s (  6.1%) Induction Variable Simplification
    0.0008s (  4.8%) Unroll loops #2
    0.0007s (  4.5%) Loop Vectorization
  == #5 Finalize object
  Percent: 46.3%
  Total 0.0174s
  Top timings:
    0.0060s ( 34.6%) X86 DAG->DAG Instruction Selection #2
    0.0019s ( 11.0%) Greedy Register Allocator #2
    0.0013s (  7.4%) Machine Instruction Scheduler #2
    0.0012s (  7.1%) Loop Strength Reduction
    0.0004s (  2.3%) Induction Variable Users


API for custom analysis
~~~~~~~~~~~~~~~~~~~~~~~

It is possible to get more details then the summary text in the above example.
The pass timings are stored in a
:class:`numba.misc.llvm_pass_timings.PassTimingsCollection`, which contains
methods for accessing individual record for each pass.

.. autoclass:: numba.misc.llvm_pass_timings.PassTimingsCollection
    :members: get_total_time, list_longest_first, summary, __getitem__, __len__

.. autoclass:: numba.misc.llvm_pass_timings.ProcessedPassTimings
    :members: get_raw_data, get_total_time, list_records, list_top, summary

.. autoclass:: numba.misc.llvm_pass_timings.PassTimingRecord

.. _low-level-extending:

Low-level extension API
=======================

This extension API is available through the :mod:`numba.extending` module.
It allows you to hook directly into the Numba compilation chain.  As such,
it distinguished between several compilation phases:

* The :term:`typing` phase deduces the types of variables in a compiled
  function by looking at the operations performed.

* The :term:`lowering` phase converts high-level Python operations into
  low-level LLVM code.  This phase exploits the typing information derived
  by the typing phase.

* *Boxing* and *unboxing* convert Python objects into native values, and
  vice-versa.  They occur at the boundaries of calling a Numba function
  from the Python interpreter.


Typing
------

.. XXX the API described here can be insufficient for some use cases.
   Should we describe the whole templates menagerie?

Type inference -- or simply *typing* -- is the process of assigning
Numba types to all values involved in a function, so as to enable
efficient code generation.  Broadly speaking, typing comes in two flavours:
typing plain Python *values* (e.g. function arguments or global variables)
and typing *operations* (or *functions*) on known value types.

.. decorator:: typeof_impl.register(cls)

   Register the decorated function as typing Python values of class *cls*.
   The decorated function will be called with the signature ``(val, c)``
   where *val* is the Python value being typed and *c* is a context
   object.


.. decorator:: type_callable(func)

   Register the decorated function as typing the callable *func*.
   *func* can be either an actual Python callable or a string denoting
   a operation internally known to Numba (for example ``'getitem'``).
   The decorated function is called with a single *context* argument
   and must return a typer function.  The typer function should have
   the same signature as the function being typed, and it is called
   with the Numba *types* of the function arguments; it should return
   either the Numba type of the function's return value, or ``None``
   if inference failed.

.. function:: as_numba_type.register(py_type, numba_type)

   Register that the Python type *py_type* corresponds with the Numba type
   *numba_type*.  This can be used to register a new type or overwrite the
   existing default (e.g. to treat ``float`` as ``numba.float32`` instead of
   ``numba.float64``).

.. decorator:: as_numba_type.register

   Register the decorated function as a type inference function used by
   ``as_numba_type`` when trying to infer the Numba type of a Python type.
   The decorated function is called with a single *py_type* argument
   and returns either a corresponding Numba type, or None if it cannot infer
   that *py_type*.


Lowering
--------

The following decorators all take a type specification of some kind.
A type specification is usually a type class (such as ``types.Float``)
or a specific type instance (such as ``types.float64``).  Some values
have a special meaning:

* ``types.Any`` matches any type; this allows doing your own dispatching
  inside the implementation

* ``types.VarArg(<some type>)`` matches any number of arguments of the
  given type; it can only appear as the last type specification when
  describing a function's arguments.

A *context* argument in the following APIs is a target context providing
various utility methods for code generation (such as creating a constant,
converting from a type to another, looking up the implementation of a
specific function, etc.).  A *builder* argument is a
:class:`llvmlite.ir.IRBuilder` instance for the LLVM code being generated.

A *signature* is an object specifying the concrete type of an operation.
The ``args`` attribute of the signature is a tuple of the argument types.
The ``return_type`` attribute of the signature is the type that the
operation should return.

.. note::
   Numba always reasons on Numba types, but the values being passed
   around during lowering are LLVM values: they don't hold the required
   type information, which is why Numba types are passed explicitly too.

   LLVM has its own, very low-level type system: you can access the LLVM
   type of a value by looking up its ``.type`` attribute.


Native operations
'''''''''''''''''

.. decorator:: lower_builtin(func, typespec, ...)

   Register the decorated function as implementing the callable *func*
   for the arguments described by the given Numba *typespecs*.
   As with :func:`type_callable`, *func* can be either an actual Python
   callable or a string denoting a operation internally known to Numba
   (for example ``'getitem'``).

   The decorated function is called with four arguments
   ``(context, builder, sig, args)``.  ``sig`` is the concrete signature
   the callable is being invoked with.  ``args`` is a tuple of the values
   of the arguments the callable is being invoked with; each value in
   ``args`` corresponds to a type in ``sig.args``.  The function
   must return a value compatible with the type ``sig.return_type``.

.. decorator:: lower_getattr(typespec, name)

   Register the decorated function as implementing the attribute *name*
   of the given *typespec*.  The decorated function is called with four
   arguments ``(context, builder, typ, value)``.  *typ* is the concrete
   type the attribute is being looked up on.  *value* is the value the
   attribute is being looked up on.

.. decorator:: lower_getattr_generic(typespec)

   Register the decorated function as a fallback for attribute lookup
   on a given *typespec*.  Any attribute that does not have a corresponding
   :func:`lower_getattr` declaration will go through
   :func:`lower_getattr_generic`.  The decorated function is called with
   five arguments ``(context, builder, typ, value, name)``.  *typ*
   and *value* are as in :func:`lower_getattr`.  *name* is the name
   of the attribute being looked up.

.. decorator:: lower_cast(fromspec, tospec)

   Register the decorated function as converting from types described by
   *fromspec* to types described by *tospec*.  The decorated function
   is called with five arguments ``(context, builder, fromty, toty, value)``.
   *fromty* and *toty* are the concrete types being converted from and to,
   respectively.  *value* is the value being converted.  The function
   must return a value compatible with the type ``toty``.


Constants
'''''''''

.. decorator:: lower_constant(typespec)

   Register the decorated function as implementing the creation of
   constants for the Numba *typespec*.  The decorated function
   is called with four arguments ``(context, builder, ty, pyval)``.
   *ty* is the concrete type to create a constant for.  *pyval*
   is the Python value to convert into a LLVM constant.
   The function must return a value compatible with the type ``ty``.


Boxing and unboxing
'''''''''''''''''''

In these functions, *c* is a convenience object with several attributes:

* its ``context`` attribute is a target context as above
* its ``builder`` attribute is a :class:`llvmlite.ir.IRBuilder` as above
* its ``pyapi`` attribute is an object giving access to a subset of the
  `Python interpreter's C API <https://docs.python.org/3/c-api/index.html>`_

An object, as opposed to a native value, is a ``PyObject *`` pointer.
Such pointers can be produced or processed by the methods in the ``pyapi``
object.

.. decorator:: box(typespec)

   Register the decorated function as boxing values matching the *typespec*.
   The decorated function is called with three arguments ``(typ, val, c)``.
   *typ* is the concrete type being boxed.  *val* is the value being
   boxed.  The function should return a Python object, or NULL to signal
   an error.

.. decorator:: unbox(typespec)

   Register the decorated function as unboxing values matching the *typespec*.
   The decorated function is called with three arguments ``(typ, obj, c)``.
   *typ* is the concrete type being unboxed.  *obj* is the Python object
   (a ``PyObject *`` pointer, in C terms) being unboxed.  The function
   should return a ``NativeValue`` object giving the unboxing result value
   and an optional error bit.

.. _overloading-guide:

==============================
A guide to using ``@overload``
==============================


As mentioned in the :ref:`high-level extension API <high-level-extending>`, you
can use the ``@overload`` decorator to create a Numba implementation of a
function that can be used in :term:`nopython mode` functions. A common use case
is to re-implement NumPy functions so that they can be called in ``@jit``
decorated code. This section discusses how and when to use the ``@overload``
decorator and what contributing such a function to the Numba code base might
entail. This should help you get started when needing to use the ``@overload``
decorator or when attempting to contribute new functions to Numba itself.

The ``@overload`` decorator and it's variants are useful when you have a
third-party library that you do not control and you wish to provide Numba
compatible implementations for specific functions from that library.

Concrete Example
================

Let's assume that you are working on a minimization algorithm that makes use of
|scipy.linalg.norm|_ to find different vector norms and the `frobenius
norm <https://en.wikipedia.org/wiki/Frobenius_inner_product>`_ for matrices.
You know that only integer and real numbers will be involved. (While this may
sound like an artificial example, especially because a Numba implementation of
``numpy.linalg.norm`` exists, it is largely pedagogical and serves to
illustrate how and when to use ``@overload``).

.. |scipy.linalg.norm| replace:: ``scipy.linalg.norm``
.. _scipy.linalg.norm: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.norm.html

The skeleton might look something like this::

    def algorithm():
        # setup
        v = ...
        while True:
            # take a step
            d = scipy.linalg.norm(v)
            if d < tolerance:
                break

Now, let's further assume, that you have heard of Numba and you now wish to use
it to accelerate your function. However, after adding the
``jit(nopython=True)``
decorator, Numba complains that ``scipy.linalg.norm`` isn't supported. From
looking at the documentation, you realize that a norm is probably fairly easy
to implement using NumPy. A good starting point is the following template.

.. literalinclude:: template.py

After some deliberation and tinkering, you end up with the following code:

.. literalinclude:: mynorm.py

As you can see, the implementation only supports what you need right now:

* Only supports integer and floating-point types
* All vector norms
* Only the Frobenius norm for matrices
* Code sharing between vector and matrix implementations using
  ``@register_jitable``.
* Norms are implemented using NumPy syntax. (This is possible because
  Numba is very aware of NumPy and many functions are supported.)

So what actually happens here? The ``overload`` decorator registers a suitable
implementation for ``scipy.linalg.norm`` in case a call to this is encountered
in code that is being JIT-compiled, for example when you decorate your
``algorithm`` function with ``@jit(nopython=True)``. In that case, the function
``jit_norm`` will be called with the currently encountered types and will then
return either ``_oneD_norm_x`` in the vector case and ``_two_D_norm_2``.

You can download the example code here: :download:`mynorm.py </extending/mynorm.py>`

Implementing ``@overload`` for NumPy functions
==============================================

Numba supports NumPy through the provision of ``@jit`` compatible
re-implementations of NumPy functions. In such cases ``@overload`` is a very
convenient option for writing such implementations, however there are a few
additional things to watch out for.

* The Numba implementation should match the NumPy implementation as closely as
  feasible with respect to accepted types, arguments, raised exceptions and
  algorithmic complexity (Big-O / Landau order).

* When implementing supported argument types, bear in mind that, due to
  duck typing, NumPy does tend to accept a multitude of argument types beyond
  NumPy arrays such as scalar, list, tuple, set, iterator, generator etc.
  You will need to account for that during type inference and subsequently as
  part of the tests.

* A NumPy function may return a scalar, array or a data structure
  which matches one of its inputs, you need to be aware of type
  unification problems and dispatch to appropriate implementations. For
  example, |np.corrcoef|_ may return an array or a scalar depending on its
  inputs.

.. |np.corrcoef| replace:: ``np.corrcoef``
.. _np.corrcoef: https://docs.scipy.org/doc/numpy/reference/generated/numpy.corrcoef.html

* If you are implementing a new function, you should always update the
  `documentation
  <https://numba.pydata.org/numba-doc/latest/reference/numpysupported.html>`_.
  The sources can be found in ``docs/source/reference/numpysupported.rst``. Be
  sure to mention any limitations that your implementation has, e.g. no support
  for the ``axis`` keyword.

* When writing tests for the functionality itself, it's useful to include
  handling of non-finite values, arrays with different shapes and layouts,
  complex inputs, scalar inputs, inputs with types for which support is not
  documented (e.g. a function which the NumPy docs say requires a float or int
  input might also 'work' if given a bool or complex input).

* When writing tests for exceptions, for example if adding tests to
  ``numba/tests/test_np_functions.py``, you may encounter the following error
  message:

  .. code::

        ======================================================================
        FAIL: test_foo (numba.tests.test_np_functions.TestNPFunctions)
        ----------------------------------------------------------------------
        Traceback (most recent call last):
        File "<path>/numba/numba/tests/support.py", line 645, in tearDown
            self.memory_leak_teardown()
        File "<path>/numba/numba/tests/support.py", line 619, in memory_leak_teardown
            self.assert_no_memory_leak()
        File "<path>/numba/numba/tests/support.py", line 628, in assert_no_memory_leak
            self.assertEqual(total_alloc, total_free)
        AssertionError: 36 != 35

  This occurs because raising exceptions from jitted code leads to reference
  leaks. Ideally, you will place all exception testing in a separate test
  method and then add a call in each test to ``self.disable_leak_check()`` to
  disable the leak-check (inherit from ``numba.tests.support.TestCase`` to make
  that available).

* For many of the functions that are available in NumPy, there are
  corresponding methods defined on the NumPy ``ndarray`` type. For example, the
  function ``repeat`` is available as a NumPy module level function and a
  member function on the ``ndarray`` class.

  .. code:: python

        import numpy as np
        a = np.arange(10)
        # function
        np.repeat(a, 10)
        # method
        a.repeat(10)

  Once you have written the function implementation, you can easily use
  ``@overload_method`` and reuse it. Just be sure to check that NumPy doesn't
  diverge in the implementations of its function/method.

  As an example, the ``repeat`` function/method:

  .. code:: python

        @extending.overload_method(types.Array, 'repeat')
        def array_repeat(a, repeats):
            def array_repeat_impl(a, repeat):
                # np.repeat has already been overloaded
                return np.repeat(a, repeat)

            return array_repeat_impl

* If you need to create ancillary functions, for example to re-use a small
  utility function or to split your implementation across functions for the
  sake of readability, you can make use of the ``@register_jitable`` decorator.
  This will make those functions available from within your ``@jit`` and
  ``@overload`` decorated functions.

* The Numba continuous integration (CI) set up tests a wide variety of NumPy
  versions, you'll sometimes be alerted to a change in behaviour from some
  previous NumPy version. If you can find supporting evidence in the NumPy
  change log / repository, then you'll need to decide whether to create
  branches and attempt to replicate the logic across versions, or use a version
  gate (with associated wording in the documentation) to advertise that Numba
  replicates NumPy from some particular version onwards.

* You can look at the Numba source code for inspiration, many of the overloaded
  NumPy functions and methods are in ``numba/targets/arrayobj.py``. Below, you
  will find a list of implementations to look at that are well implemented in
  terms of accepted types and test coverage.

  * ``np.repeat``

Example: an interval type
=========================

We will extend the Numba frontend to support a class that it does not
currently support so as to allow:

* Passing an instance of the class to a Numba function
* Accessing attributes of the class in a Numba function
* Constructing and returning a new instance of the class from a Numba function

(all the above in :term:`nopython mode`)

We will mix APIs from the :ref:`high-level extension API <high-level-extending>`
and the :ref:`low-level extension API <low-level-extending>`, depending on what is
available for a given task.

The starting point for our example is the following pure Python class::

   class Interval(object):
       """
       A half-open interval on the real number line.
       """
       def __init__(self, lo, hi):
           self.lo = lo
           self.hi = hi

       def __repr__(self):
           return 'Interval(%f, %f)' % (self.lo, self.hi)

       @property
       def width(self):
           return self.hi - self.lo


Extending the typing layer
""""""""""""""""""""""""""

Creating a new Numba type
-------------------------

As the ``Interval`` class is not known to Numba, we must create a new Numba
type to represent instances of it.  Numba does not deal with Python types
directly: it has its own type system that allows a different level of
granularity as well as various meta-information not available with regular
Python types.

We first create a type class ``IntervalType`` and, since we don't need the
type to be parametric, we instantiate a single type instance ``interval_type``::

   from numba import types

   class IntervalType(types.Type):
       def __init__(self):
           super(IntervalType, self).__init__(name='Interval')

   interval_type = IntervalType()


Type inference for Python values
--------------------------------

In itself, creating a Numba type doesn't do anything.  We must teach Numba
how to infer some Python values as instances of that type.  In this example,
it is trivial: any instance of the ``Interval`` class should be treated as
belonging to the type ``interval_type``::

   from numba.extending import typeof_impl

   @typeof_impl.register(Interval)
   def typeof_index(val, c):
       return interval_type

Function arguments and global values will thusly be recognized as belonging
to ``interval_type`` whenever they are instances of ``Interval``.


Type inference for Python annotations
-------------------------------------

While ``typeof`` is used to infer the Numba type of Python objects,
``as_numba_type`` is used to infer the Numba type of Python types.  For simple
cases, we can simply register that the Python type ``Interval`` corresponds with
the Numba type ``interval_type``::

   from numba.extending import as_numba_type

   as_numba_type.register(Interval, interval_type)

Note that ``as_numba_type`` is only used to infer types from type annotations at
compile time.  The ``typeof`` registry above is used to infer the type of
objects at runtime.


Type inference for operations
-----------------------------

We want to be able to construct interval objects from Numba functions, so
we must teach Numba to recognize the two-argument ``Interval(lo, hi)``
constructor.  The arguments should be floating-point numbers::

   from numba.extending import type_callable

   @type_callable(Interval)
   def type_interval(context):
       def typer(lo, hi):
           if isinstance(lo, types.Float) and isinstance(hi, types.Float):
               return interval_type
       return typer


The :func:`type_callable` decorator specifies that the decorated function
should be invoked when running type inference for the given callable object
(here the ``Interval`` class itself).  The decorated function must simply
return a typer function that will be called with the argument types.  The
reason for this seemingly convoluted setup is for the typer function to have
*exactly* the same signature as the typed callable.  This allows handling
keyword arguments correctly.

The *context* argument received by the decorated function is useful in
more sophisticated cases where computing the callable's return type
requires resolving other types.


Extending the lowering layer
""""""""""""""""""""""""""""

We have finished teaching Numba about our type inference additions.
We must now teach Numba how to actually generate code and data for
the new operations.


Defining the data model for native intervals
--------------------------------------------

As a general rule, :term:`nopython mode` does not work on Python objects
as they are generated by the CPython interpreter.  The representations
used by the interpreter are far too inefficient for fast native code.
Each type supported in :term:`nopython mode` therefore has to define
a tailored native representation, also called a *data model*.

A common case of data model is an immutable struct-like data model, that
is akin to a C ``struct``.  Our interval datatype conveniently falls in
that category, and here is a possible data model for it::

   from numba.extending import models, register_model

   @register_model(IntervalType)
   class IntervalModel(models.StructModel):
       def __init__(self, dmm, fe_type):
           members = [
               ('lo', types.float64),
               ('hi', types.float64),
               ]
           models.StructModel.__init__(self, dmm, fe_type, members)


This instructs Numba that values of type ``IntervalType`` (or any instance
thereof) are represented as a structure of two fields ``lo`` and ``hi``,
each of them a double-precision floating-point number (``types.float64``).

.. note::
   Mutable types need more sophisticated data models to be able to
   persist their values after modification.  They typically cannot be
   stored and passed on the stack or in registers like immutable types do.


Exposing data model attributes
------------------------------

We want the data model attributes ``lo`` and ``hi`` to be exposed under
the same names for use in Numba functions.  Numba provides a convenience
function to do exactly that::

   from numba.extending import make_attribute_wrapper

   make_attribute_wrapper(IntervalType, 'lo', 'lo')
   make_attribute_wrapper(IntervalType, 'hi', 'hi')

This will expose the attributes in read-only mode.  As mentioned above,
writable attributes don't fit in this model.


Exposing a property
-------------------

As the ``width`` property is computed rather than stored in the structure,
we cannot simply expose it like we did for ``lo`` and ``hi``.  We have to
re-implement it explicitly::

   from numba.extending import overload_attribute

   @overload_attribute(IntervalType, "width")
   def get_width(interval):
       def getter(interval):
           return interval.hi - interval.lo
       return getter

You might ask why we didn't need to expose a type inference hook for this
attribute? The answer is that ``@overload_attribute`` is part of the
high-level API: it combines type inference and code generation in a
single API.


Implementing the constructor
----------------------------

Now we want to implement the two-argument ``Interval`` constructor::

   from numba.extending import lower_builtin
   from numba.core import cgutils

   @lower_builtin(Interval, types.Float, types.Float)
   def impl_interval(context, builder, sig, args):
       typ = sig.return_type
       lo, hi = args
       interval = cgutils.create_struct_proxy(typ)(context, builder)
       interval.lo = lo
       interval.hi = hi
       return interval._getvalue()


There is a bit more going on here.  ``@lower_builtin`` decorates the
implementation of the given callable or operation (here the ``Interval``
constructor) for some specific argument types.  This allows defining
type-specific implementations of a given operation, which is important
for heavily overloaded functions such as :func:`len`.

``types.Float`` is the class of all floating-point types (``types.float64``
is an instance of ``types.Float``).  It is generally more future-proof
to match argument types on their class rather than on specific instances
(however, when *returning* a type -- chiefly during the type inference
phase --, you must usually return a type instance).

``cgutils.create_struct_proxy()`` and ``interval._getvalue()`` are a bit
of boilerplate due to how Numba passes values around.  Values are passed
as instances of :class:`llvmlite.ir.Value`, which can be too limited:
LLVM structure values especially are quite low-level.  A struct proxy
is a temporary wrapper around a LLVM structure value allowing to easily
get or set members of the structure. The ``_getvalue()`` call simply
gets the LLVM value out of the wrapper.


Boxing and unboxing
-------------------

If you try to use an ``Interval`` instance at this point, you'll certainly
get the error *"cannot convert Interval to native value"*.  This is because
Numba doesn't yet know how to make a native interval value from a Python
``Interval`` instance.  Let's teach it how to do it::

   from numba.extending import unbox, NativeValue

   @unbox(IntervalType)
   def unbox_interval(typ, obj, c):
       """
       Convert a Interval object to a native interval structure.
       """
       lo_obj = c.pyapi.object_getattr_string(obj, "lo")
       hi_obj = c.pyapi.object_getattr_string(obj, "hi")
       interval = cgutils.create_struct_proxy(typ)(c.context, c.builder)
       interval.lo = c.pyapi.float_as_double(lo_obj)
       interval.hi = c.pyapi.float_as_double(hi_obj)
       c.pyapi.decref(lo_obj)
       c.pyapi.decref(hi_obj)
       is_error = cgutils.is_not_null(c.builder, c.pyapi.err_occurred())
       return NativeValue(interval._getvalue(), is_error=is_error)

*Unbox* is the other name for "convert a Python object to a native value"
(it fits the idea of a Python object as a sophisticated box containing
a simple native value).  The function returns a ``NativeValue`` object
which gives its caller access to the computed native value, the error bit
and possibly other information.

The snippet above makes abundant use of the ``c.pyapi`` object, which
gives access to a subset of the
`Python interpreter's C API <https://docs.python.org/3/c-api/index.html>`_.
Note the use of ``c.pyapi.err_occurred()`` to detect any errors that
may have happened when unboxing the object (try passing ``Interval('a', 'b')``
for example).

We also want to do the reverse operation, called *boxing*, so as to return
interval values from Numba functions::

   from numba.extending import box

   @box(IntervalType)
   def box_interval(typ, val, c):
       """
       Convert a native interval structure to an Interval object.
       """
       interval = cgutils.create_struct_proxy(typ)(c.context, c.builder, value=val)
       lo_obj = c.pyapi.float_from_double(interval.lo)
       hi_obj = c.pyapi.float_from_double(interval.hi)
       class_obj = c.pyapi.unserialize(c.pyapi.serialize_object(Interval))
       res = c.pyapi.call_function_objargs(class_obj, (lo_obj, hi_obj))
       c.pyapi.decref(lo_obj)
       c.pyapi.decref(hi_obj)
       c.pyapi.decref(class_obj)
       return res


Using it
""""""""

:term:`nopython mode` functions are now able to make use of Interval objects
and the various operations you have defined on them.  You can try for
example the following functions::

   from numba import jit

   @jit(nopython=True)
   def inside_interval(interval, x):
       return interval.lo <= x < interval.hi

   @jit(nopython=True)
   def interval_width(interval):
       return interval.width

   @jit(nopython=True)
   def sum_intervals(i, j):
       return Interval(i.lo + j.lo, i.hi + j.hi)


Conclusion
""""""""""

We have shown how to do the following tasks:

* Define a new Numba type class by subclassing the ``Type`` class
* Define a singleton Numba type instance for a non-parametric type
* Teach Numba how to infer the Numba type of Python values of a certain class,
  using ``typeof_impl.register``
* Teach Numba how to infer the Numba type of the Python type itself, using
  ``as_numba_type.register``
* Define the data model for a Numba type using ``StructModel``
  and ``register_model``
* Implementing a boxing function for a Numba type using the ``@box`` decorator
* Implementing an unboxing function for a Numba type using the ``@unbox`` decorator
  and the ``NativeValue`` class
* Type and implement a callable using the ``@type_callable`` and
  ``@lower_builtin`` decorators
* Expose a read-only structure attribute using the ``make_attribute_wrapper``
  convenience function
* Implement a read-only property using the ``@overload_attribute`` decorator
Registering Extensions with Entry Points
========================================

Often, third party packages will have a user-facing API as well as define
extensions to the Numba compiler.  In those situations, the new types and
overloads can registered with Numba when the package is imported by the user.
However, there are situations where a Numba extension would not normally be
imported directly by the user, but must still be registered with the Numba
compiler.  An example of this is the `numba-scipy
<https://github.com/numba/numba-scipy>`_ package, which adds support for some
SciPy functions to Numba.  The end user does not need to ``import
numba_scipy`` to enable compiler support for SciPy, the extension only needs
to be installed in the Python environment.

Numba discovers extensions using the `entry points
<https://setuptools.readthedocs.io/en/latest/setuptools.html#dynamic-discovery-of-services-and-plugins>`_
feature of ``setuptools``.  This allows a Python package to register an
initializer function that will be called before ``numba`` compiles for the
first time.  The delay ensures that the cost of importing extensions is
deferred until it is necessary.


Adding Support for the "Init" Entry Point
-----------------------------------------

A package can register an initialization function with Numba by adding the
``entry_points`` argument to the ``setup()`` function call in ``setup.py``:

.. code-block:: python

    setup(
        ...,
        entry_points={
            "numba_extensions": [
                "init = numba_scipy:_init_extension",
            ],
        },
        ...
    )

Numba currently only looks for the ``init`` entry point in the
``numba_extensions`` group.  The entry point should be a function (any name,
as long as it matches what is listed in ``setup.py``) that takes no arguments,
and the return value is ignored.  This function should register types,
overloads, or call other Numba extension APIs.  The order of initialization of
extensions is undefined.

Testing your Entry Point
------------------------

Numba loads all entry points when the first function is compiled. To test your
entry point, it is not sufficient to just ``import numba``; you have to define
and run a small function, like this:

.. code-block:: python

    import numba; numba.njit(lambda x: x + 1)(123)

It is not necessary to import your module: entry points are identified by the
``entry_points.txt`` file in your library's ``*.egg-info`` directory.

The ``setup.py build`` command does not create eggs, but ``setup.py sdist``
(for testing in a local directory) and ``setup.py install`` do. All entry points
registered in eggs that are on the Python path are loaded. Be sure to check for
stale ``entry_points.txt`` when debugging.

Extending Numba
===============

.. module:: numba.extending

This chapter describes how to extend Numba to make it recognize and support
additional operations, functions or types.  Numba provides two categories
of APIs to this end:

* The high-level APIs provide abstracted entry points which are sufficient
  for simple uses.  They require little knowledge of Numba's internal
  compilation chain.

* The low-level APIs reflect Numba's internal compilation chain and allow
  flexible interaction with its various layers, but require more effort
  and experience with Numba internals.

It may be helpful for readers of this chapter to also read some of the
documents in the :doc:`developer manual <../developer/index>`, especially
the :doc:`architecture document <../developer/architecture>`.


.. toctree::
   high-level.rst
   low-level.rst
   interval-example.rst
   overloading-guide.rst
   entrypoints.rst


.. _high-level-extending:

High-level extension API
========================

This extension API is exposed through the :mod:`numba.extending` module.

To aid debugging extensions to Numba, it's recommended to set the following
environment variable::

    NUMBA_CAPTURED_ERRORS="new_style"

this makes it easy to differentiate between errors in implementation and
acceptable errors that can take part in e.g. type inference. For more
information see :envvar:`NUMBA_CAPTURED_ERRORS`.

Implementing functions
----------------------

The ``@overload`` decorator allows you to implement arbitrary functions
for use in :term:`nopython mode` functions.  The function decorated with
``@overload`` is called at compile-time with the *types* of the function's
runtime arguments.  It should return a callable representing the
*implementation* of the function for the given types.  The returned
implementation is compiled by Numba as if it were a normal function
decorated with ``@jit``.  Additional options to ``@jit`` can be passed as
dictionary using the ``jit_options`` argument.

For example, let's pretend Numba doesn't support the :func:`len` function
on tuples yet.  Here is how to implement it using ``@overload``::

   from numba import types
   from numba.extending import overload

   @overload(len)
   def tuple_len(seq):
      if isinstance(seq, types.BaseTuple):
          n = len(seq)
          def len_impl(seq):
              return n
          return len_impl


You might wonder, what happens if :func:`len()` is called with something
else than a tuple? If a function decorated with ``@overload`` doesn't
return anything (i.e. returns None), other definitions are tried until
one succeeds.  Therefore, multiple libraries may overload :func:`len()`
for different types without conflicting with each other.

Implementing methods
--------------------

The ``@overload_method`` decorator similarly allows implementing a
method on a type well-known to Numba.

.. autofunction:: numba.core.extending.overload_method

Implementing classmethods
-------------------------

The ``@overload_classmethod`` decorator similarly allows implementing a
classmethod on a type well-known to Numba.

.. autofunction:: numba.core.extending.overload_classmethod


Implementing attributes
-----------------------

The ``@overload_attribute`` decorator allows implementing a data
attribute (or property) on a type.  Only reading the attribute is
possible; writable attributes are only supported through the
:ref:`low-level API <low-level-extending>`.

The following example implements the :attr:`~numpy.ndarray.nbytes` attribute
on Numpy arrays::

   @overload_attribute(types.Array, 'nbytes')
   def array_nbytes(arr):
      def get(arr):
          return arr.size * arr.itemsize
      return get

.. _cython-support:

Importing Cython Functions
--------------------------

The function ``get_cython_function_address`` obtains the address of a
C function in a Cython extension module. The address can be used to
access the C function via a :func:`ctypes.CFUNCTYPE` callback, thus
allowing use of the C function inside a Numba jitted function. For
example, suppose that you have the file ``foo.pyx``::

   from libc.math cimport exp

   cdef api double myexp(double x):
       return exp(x)

You can access ``myexp`` from Numba in the following way::

   import ctypes
   from numba.extending import get_cython_function_address

   addr = get_cython_function_address("foo", "myexp")
   functype = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
   myexp = functype(addr)

The function ``myexp`` can now be used inside jitted functions, for
example::

   @njit
   def double_myexp(x):
       return 2*myexp(x)

One caveat is that if your function uses Cython's fused types, then
the function's name will be mangled. To find out the mangled name of
your function you can check the extension module's ``__pyx_capi__``
attribute.

Implementing intrinsics
-----------------------

The ``@intrinsic`` decorator is used for marking a function *func* as typing and
implementing the function in ``nopython`` mode using the
`llvmlite IRBuilder API <http://llvmlite.pydata.org/en/latest/user-guide/ir/ir-builder.html>`_.
This is an escape hatch for expert users to build custom LLVM IR that will be
inlined into the caller, there is no safety net!

The first argument to *func* is the typing context.  The rest of the arguments
corresponds to the type of arguments of the decorated function. These arguments
are also used as the formal argument of the decorated function.  If *func* has
the signature ``foo(typing_context, arg0, arg1)``, the decorated function will
have the signature ``foo(arg0, arg1)``.

The return values of *func* should be a 2-tuple of expected type signature, and
a code-generation function that will passed to
:func:`~numba.targets.imputils.lower_builtin`. For an unsupported operation,
return ``None``.

Here is an example that cast any integer to a byte pointer::

    from numba import types
    from numba.extending import intrinsic

    @intrinsic
    def cast_int_to_byte_ptr(typingctx, src):
        # check for accepted types
        if isinstance(src, types.Integer):
            # create the expected type signature
            result_type = types.CPointer(types.uint8)
            sig = result_type(types.uintp)
            # defines the custom code generation
            def codegen(context, builder, signature, args):
                # llvm IRBuilder code here
                [src] = args
                rtype = signature.return_type
                llrtype = context.get_value_type(rtype)
                return builder.inttoptr(src, llrtype)
            return sig, codegen

it may be used as follows::

    from numba import njit

    @njit('void(int64)')
    def foo(x):
        y = cast_int_to_byte_ptr(x)

    foo.inspect_types()

and the output of ``.inspect_types()`` demonstrates the cast (note the
``uint8*``)::

    def foo(x):

        #   x = arg(0, name=x)  :: int64
        #   $0.1 = global(cast_int_to_byte_ptr: <intrinsic cast_int_to_byte_ptr>)  :: Function(<intrinsic cast_int_to_byte_ptr>)
        #   $0.3 = call $0.1(x, func=$0.1, args=[Var(x, check_intrin.py (24))], kws=(), vararg=None)  :: (uint64,) -> uint8*
        #   del x
        #   del $0.1
        #   y = $0.3  :: uint8*
        #   del y
        #   del $0.3
        #   $const0.4 = const(NoneType, None)  :: none
        #   $0.5 = cast(value=$const0.4)  :: none
        #   del $const0.4
        #   return $0.5

        y = cast_int_to_byte_ptr(x)


Implementing mutable structures
-------------------------------

.. warning:: This is an experimental feature, the API may change without warning.

The ``numba.experimental.structref`` module provides utilities for defining
mutable pass-by-reference structures, a ``StructRef``. The following example
demonstrates how to define a basic mutable structure:

Defining a StructRef
''''''''''''''''''''

.. literalinclude:: ../../../numba/tests/doc_examples/test_structref_usage.py
   :language: python
   :caption: from ``numba/tests/doc_examples/test_structref_usage.py``
   :start-after: magictoken.ex_structref_type_definition.begin
   :end-before: magictoken.ex_structref_type_definition.end
   :dedent: 0
   :linenos:

The following demonstrates using the above mutable struct definition:

.. literalinclude:: ../../../numba/tests/doc_examples/test_structref_usage.py
   :language: python
   :caption: from ``test_type_definition`` of ``numba/tests/doc_examples/test_structref_usage.py``
   :start-after: magictoken.ex_structref_type_definition_test.begin
   :end-before: magictoken.ex_structref_type_definition_test.end
   :dedent: 8
   :linenos:


Defining a method on StructRef
''''''''''''''''''''''''''''''

Methods and attributes can be attached using ``@overload_*`` as shown in the
previous sections.

The following demonstrates the use of ``@overload_method`` to insert a
method for instances of ``MyStructType``:

.. literalinclude:: ../../../numba/tests/doc_examples/test_structref_usage.py
   :language: python
   :caption: from ``test_overload_method`` of ``numba/tests/doc_examples/test_structref_usage.py``
   :start-after: magictoken.ex_structref_method.begin
   :end-before: magictoken.ex_structref_method.end
   :dedent: 8
   :linenos:


``numba.experimental.structref`` API Reference
''''''''''''''''''''''''''''''''''''''''''''''

.. automodule:: numba.experimental.structref
    :members:

Determining if a function is already wrapped by a ``jit`` family decorator
--------------------------------------------------------------------------

The following function is provided for this purpose.

.. automethod:: numba.extending.is_jitted
Memory Management
=================

.. autofunction:: numba.cuda.to_device
.. autofunction:: numba.cuda.device_array
.. autofunction:: numba.cuda.device_array_like
.. autofunction:: numba.cuda.pinned_array
.. autofunction:: numba.cuda.pinned_array_like
.. autofunction:: numba.cuda.mapped_array
.. autofunction:: numba.cuda.mapped_array_like
.. autofunction:: numba.cuda.managed_array
.. autofunction:: numba.cuda.pinned
.. autofunction:: numba.cuda.mapped

Device Objects
--------------

.. autoclass:: numba.cuda.cudadrv.devicearray.DeviceNDArray
   :members: copy_to_device, copy_to_host, is_c_contiguous, is_f_contiguous,
              ravel, reshape, split
.. autoclass:: numba.cuda.cudadrv.devicearray.DeviceRecord
   :members: copy_to_device, copy_to_host
.. autoclass:: numba.cuda.cudadrv.devicearray.MappedNDArray
   :members: copy_to_device, copy_to_host, split
CUDA Host API
=============

Device Management
-----------------

Device detection and enquiry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following functions are available for querying the available hardware:

.. autofunction:: numba.cuda.is_available

.. autofunction:: numba.cuda.detect

Context management
~~~~~~~~~~~~~~~~~~

CUDA Python functions execute within a CUDA context. Each CUDA device in a
system has an associated CUDA context, and Numba presently allows only one context
per thread. For further details on CUDA Contexts, refer to the `CUDA Driver API
Documentation on Context Management
<http://docs.nvidia.com/cuda/cuda-driver-api/group__CUDA__CTX.html>`_ and the
`CUDA C Programming Guide Context Documentation
<http://docs.nvidia.com/cuda/cuda-c-programming-guide/#context>`_. CUDA Contexts
are instances of the :class:`~numba.cuda.cudadrv.driver.Context` class:

.. autoclass:: numba.cuda.cudadrv.driver.Context
   :members: reset, get_memory_info, push, pop

The following functions can be used to get or select the context:

.. autofunction:: numba.cuda.current_context
.. autofunction:: numba.cuda.require_context

The following functions affect the current context:

.. autofunction:: numba.cuda.synchronize
.. autofunction:: numba.cuda.close

Device management
~~~~~~~~~~~~~~~~~

Numba maintains a list of supported CUDA-capable devices:

.. attribute:: numba.cuda.gpus

   An indexable list of supported CUDA devices. This list is indexed by integer
   device ID.

Alternatively, the current device can be obtained:

.. function:: numba.cuda.gpus.current

   Return the currently-selected device.

Getting a device through :attr:`numba.cuda.gpus` always provides an instance of
:class:`numba.cuda.cudadrv.devices._DeviceContextManager`, which acts as a
context manager for the selected device:

.. autoclass:: numba.cuda.cudadrv.devices._DeviceContextManager

One may also select a context and device or get the current device using the
following three functions:

.. autofunction:: numba.cuda.select_device
.. autofunction:: numba.cuda.get_current_device
.. autofunction:: numba.cuda.list_devices

The :class:`numba.cuda.cudadrv.driver.Device` class can be used to enquire about
the functionality of the selected device:

.. class:: numba.cuda.cudadrv.driver.Device

   The device associated with a particular context.

   .. attribute:: compute_capability

      A tuple, *(major, minor)* indicating the supported compute capability.

   .. attribute:: id

      The integer ID of the device.

   .. attribute:: name

      The name of the device (e.g. "GeForce GTX 970").

   .. attribute:: uuid

      The UUID of the device (e.g. "GPU-e6489c45-5b68-3b03-bab7-0e7c8e809643").

   .. method:: reset

      Delete the context for the device. This will destroy all memory
      allocations, events, and streams created within the context.


Compilation
-----------

Numba provides an entry point for compiling a Python function to PTX without
invoking any of the driver API. This can be useful for:

- Generating PTX that is to be inlined into other PTX code (e.g. from outside
  the Numba / Python ecosystem).
- Generating code when there is no device present.
- Generating code prior to a fork without initializing CUDA.

.. note:: It is the user's responsibility to manage any ABI issues arising from
   the use of compilation to PTX.

.. autofunction:: numba.cuda.compile_ptx


The environment variable ``NUMBA_CUDA_DEFAULT_PTX_CC`` can be set to control
the default compute capability targeted by ``compile_ptx`` - see
:ref:`numba-envvars-gpu-support`. If PTX for the compute capability of the
current device is required, the ``compile_ptx_for_current_device`` function can
be used:

.. autofunction:: numba.cuda.compile_ptx_for_current_device



Measurement
-----------

.. _cuda-profiling:

Profiling
~~~~~~~~~

The NVidia Visual Profiler can be used directly on executing CUDA Python code -
it is not a requirement to insert calls to these functions into user code.
However, these functions can be used to allow profiling to be performed
selectively on specific portions of the code. For further information on
profiling, see the `NVidia Profiler User's Guide
<https://docs.nvidia.com/cuda/profiler-users-guide/>`_.

.. autofunction:: numba.cuda.profile_start
.. autofunction:: numba.cuda.profile_stop
.. autofunction:: numba.cuda.profiling


.. _events:

Events
~~~~~~

Events can be used to monitor the progress of execution and to record the
timestamps of specific points being reached. Event creation returns immediately,
and the created event can be queried to determine if it has been reached. For
further information, see the `CUDA C Programming Guide Events section
<http://docs.nvidia.com/cuda/cuda-c-programming-guide/#events>`_.

The following functions are used for creating and measuring the time between
events:

.. autofunction:: numba.cuda.event
.. autofunction:: numba.cuda.event_elapsed_time

Events are instances of the :class:`numba.cuda.cudadrv.driver.Event` class:

.. autoclass:: numba.cuda.cudadrv.driver.Event
   :members: query, record, synchronize, wait


.. _streams:

Stream Management
-----------------

Streams allow concurrency of execution on a single device within a given
context. Queued work items in the same stream execute sequentially, but work
items in different streams may execute concurrently. Most operations involving a
CUDA device can be performed asynchronously using streams, including data
transfers and kernel execution. For further details on streams, see the `CUDA C
Programming Guide Streams section
<http://docs.nvidia.com/cuda/cuda-c-programming-guide/#streams>`_.

Numba defaults to using the legacy default stream as the default stream. The
per-thread default stream can be made the default stream by setting the
environment variable ``NUMBA_CUDA_PER_THREAD_DEFAULT_STREAM`` to ``1`` (see the
:ref:`CUDA Environment Variables section <numba-envvars-gpu-support>`).
Regardless of this setting, the objects representing the legacy and per-thread
default streams can be constructed using the functions below.

Streams are instances of :class:`numba.cuda.cudadrv.driver.Stream`:

.. autoclass:: numba.cuda.cudadrv.driver.Stream
   :members: synchronize, auto_synchronize, add_callback, async_done

To create a new stream:

.. autofunction:: numba.cuda.stream

To get the default stream:

.. autofunction:: numba.cuda.default_stream

To get the default stream with an explicit choice of whether it is the legacy
or per-thread default stream:

.. autofunction:: numba.cuda.legacy_default_stream

.. autofunction:: numba.cuda.per_thread_default_stream

To construct a Numba ``Stream`` object using a stream allocated elsewhere, the
``external_stream`` function is provided. Note that the lifetime of external
streams must be managed by the user - Numba will not deallocate an external
stream, and the stream must remain valid whilst the Numba ``Stream`` object is
in use.

.. autofunction:: numba.cuda.external_stream


Runtime
-------

Numba generally uses the Driver API, but it provides a simple wrapper to the
Runtime API so that the version of the runtime in use can be queried. This is
accessed through ``cuda.runtime``, which is an instance of the
:class:`numba.cuda.cudadrv.runtime.Runtime` class:

.. autoclass:: numba.cuda.cudadrv.runtime.Runtime
   :members: get_version, is_supported_version, supported_versions

Whether the current runtime is officially supported and tested with the current
version of Numba can also be queried:

.. autofunction:: numba.cuda.is_supported_version
Libdevice functions
===================

All wrapped libdevice functions are listed in this section. All functions in
libdevice are wrapped, with the exception of ``__nv_nan`` and ``__nv_nanf``.
These functions return a representation of a quiet NaN, but the argument they
take (a pointer to an object specifying the representation) is undocumented, and
follows an unusual form compared to the rest of libdevice - it is not an output
like every other pointer argument. If a NaN is required, one can be obtained in
CUDA Python by other means, e.g. ``math.nan``.

Wrapped functions
-----------------

.. automodule:: numba.cuda.libdevice
   :members:
CUDA Kernel API
===============

Kernel declaration
------------------

The ``@cuda.jit`` decorator is used to create a CUDA dispatcher object that can
be configured and launched:

.. autofunction:: numba.cuda.jit


Dispatcher objects
------------------

The usual syntax for configuring a Dispatcher with a launch configuration uses
subscripting, with the arguments being as in the following:

.. code-block:: python

   # func is some function decorated with @cuda.jit
   func[griddim, blockdim, stream, sharedmem]


The ``griddim`` and ``blockdim`` arguments specify the size of the grid and
thread blocks, and may be either integers or tuples of length up to 3. The
``stream`` parameter is an optional stream on which the kernel will be launched,
and the ``sharedmem`` parameter specifies the size of dynamic shared memory in
bytes.

Subscripting the Dispatcher returns a configuration object that can be called
with the kernel arguments:

.. code-block:: python

   configured = func[griddim, blockdim, stream, sharedmem]
   configured(x, y, z)


However, it is more idiomatic to configure and call the kernel within a single
statement:

.. code-block:: python

   func[griddim, blockdim, stream, sharedmem](x, y, z)

This is similar to launch configuration in CUDA C/C++:

.. code-block:: cuda

   func<<<griddim, blockdim, sharedmem, stream>>>(x, y, z)

.. note:: The order of ``stream`` and ``sharedmem`` are reversed in Numba
   compared to in CUDA C/C++.

Dispatcher objects also provide several utility methods for inspection and
creating a specialized instance:

.. autoclass:: numba.cuda.compiler.Dispatcher
   :members: inspect_asm, inspect_llvm, inspect_sass, inspect_types,
             get_regs_per_thread, specialize, specialized, extensions, forall


Intrinsic Attributes and Functions
----------------------------------

The remainder of the attributes and functions in this section may only be called
from within a CUDA Kernel.

Thread Indexing
~~~~~~~~~~~~~~~

.. attribute:: numba.cuda.threadIdx

    The thread indices in the current thread block, accessed through the
    attributes ``x``, ``y``, and ``z``. Each index is an integer spanning the
    range from 0 inclusive to the corresponding value of the attribute in
    :attr:`numba.cuda.blockDim` exclusive.

.. attribute:: numba.cuda.blockIdx

    The block indices in the grid of thread blocks, accessed through the
    attributes ``x``, ``y``, and ``z``. Each index is an integer spanning the
    range from 0 inclusive to the corresponding value of the attribute in
    :attr:`numba.cuda.gridDim` exclusive.

.. attribute:: numba.cuda.blockDim

    The shape of a block of threads, as declared when instantiating the
    kernel.  This value is the same for all threads in a given kernel, even
    if they belong to different blocks (i.e. each block is "full").

.. attribute:: numba.cuda.gridDim

    The shape of the grid of blocks, accessed through the attributes ``x``,
    ``y``, and ``z``.

.. attribute:: numba.cuda.laneid

    The thread index in the current warp, as an integer spanning the range
    from 0 inclusive to the :attr:`numba.cuda.warpsize` exclusive.

.. attribute:: numba.cuda.warpsize

    The size in threads of a warp on the GPU. Currently this is always 32.

.. function:: numba.cuda.grid(ndim)

   Return the absolute position of the current thread in the entire
   grid of blocks.  *ndim* should correspond to the number of dimensions
   declared when instantiating the kernel.  If *ndim* is 1, a single integer
   is returned.  If *ndim* is 2 or 3, a tuple of the given number of
   integers is returned.

   Computation of the first integer is as follows::

      cuda.threadIdx.x + cuda.blockIdx.x * cuda.blockDim.x

   and is similar for the other two indices, but using the ``y`` and ``z``
   attributes.

.. function:: numba.cuda.gridsize(ndim)

   Return the absolute size (or shape) in threads of the entire grid of
   blocks. *ndim* should correspond to the number of dimensions declared when
   instantiating the kernel.

   Computation of the first integer is as follows::

       cuda.blockDim.x * cuda.gridDim.x

   and is similar for the other two indices, but using the ``y`` and ``z``
   attributes.

Memory Management
~~~~~~~~~~~~~~~~~

.. function:: numba.cuda.shared.array(shape, dtype)

   Creates an array in the local memory space of the CUDA kernel with
   the given ``shape`` and ``dtype``.

   Returns an array with its content uninitialized.

   .. note:: All threads in the same thread block sees the same array.

.. function:: numba.cuda.local.array(shape, dtype)

   Creates an array in the local memory space of the CUDA kernel with the
   given ``shape`` and ``dtype``.

   Returns an array with its content uninitialized.

   .. note:: Each thread sees a unique array.

.. function:: numba.cuda.const.array_like(ary)

   Copies the ``ary`` into constant memory space on the CUDA kernel at compile
   time.

   Returns an array like the ``ary`` argument.

   .. note:: All threads and blocks see the same array.

Synchronization and Atomic Operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: numba.cuda.atomic.add(array, idx, value)

    Perform ``array[idx] += value``. Support int32, int64, float32 and
    float64 only. The ``idx`` argument can be an integer or a tuple of integer
    indices for indexing into multiple dimensional arrays. The number of element
    in ``idx`` must match the number of dimension of ``array``.

    Returns the value of ``array[idx]`` before the storing the new value.
    Behaves like an atomic load.

.. function:: numba.cuda.atomic.sub(array, idx, value)

    Perform ``array[idx] -= value``. Supports int32, int64, float32 and
    float64 only. The ``idx`` argument can be an integer or a tuple of integer
    indices for indexing into multi-dimensional arrays. The number of elements
    in ``idx`` must match the number of dimensions of ``array``.

    Returns the value of ``array[idx]`` before the storing the new value.
    Behaves like an atomic load.

.. function:: numba.cuda.atomic.and_(array, idx, value)

    Perform ``array[idx] &= value``. Supports int32, uint32, int64,
    and uint64 only. The ``idx`` argument can be an integer or a tuple of
    integer indices for indexing into multi-dimensional arrays. The number
    of elements in ``idx`` must match the number of dimensions of ``array``.

    Returns the value of ``array[idx]`` before the storing the new value.
    Behaves like an atomic load.

.. function:: numba.cuda.atomic.or_(array, idx, value)

    Perform ``array[idx] |= value``. Supports int32, uint32, int64,
    and uint64 only. The ``idx`` argument can be an integer or a tuple of
    integer indices for indexing into multi-dimensional arrays. The number
    of elements in ``idx`` must match the number of dimensions of ``array``.

    Returns the value of ``array[idx]`` before the storing the new value.
    Behaves like an atomic load.

.. function:: numba.cuda.atomic.xor(array, idx, value)

    Perform ``array[idx] ^= value``. Supports int32, uint32, int64,
    and uint64 only. The ``idx`` argument can be an integer or a tuple of
    integer indices for indexing into multi-dimensional arrays. The number
    of elements in ``idx`` must match the number of dimensions of ``array``.

    Returns the value of ``array[idx]`` before the storing the new value.
    Behaves like an atomic load.

.. function:: numba.cuda.atomic.exch(array, idx, value)

    Perform ``array[idx] = value``. Supports int32, uint32, int64,
    and uint64 only. The ``idx`` argument can be an integer or a tuple of
    integer indices for indexing into multi-dimensional arrays. The number
    of elements in ``idx`` must match the number of dimensions of ``array``.

    Returns the value of ``array[idx]`` before the storing the new value.
    Behaves like an atomic load.

.. function:: numba.cuda.atomic.inc(array, idx, value)

    Perform ``array[idx] = (0 if array[idx] >= value else array[idx] + 1)``.
    Supports uint32, and uint64 only. The ``idx`` argument can be an integer
    or a tuple of integer indices for indexing into multi-dimensional arrays.
    The number of elements in ``idx`` must match the number of dimensions of
    ``array``.

    Returns the value of ``array[idx]`` before the storing the new value.
    Behaves like an atomic load.

.. function:: numba.cuda.atomic.dec(array, idx, value)

    Perform ``array[idx] =
    (value if (array[idx] == 0) or (array[idx] > value) else array[idx] - 1)``.
    Supports uint32, and uint64 only. The ``idx`` argument can be an integer
    or a tuple of integer indices for indexing into multi-dimensional arrays.
    The number of elements in ``idx`` must match the number of dimensions of
    ``array``.

    Returns the value of ``array[idx]`` before the storing the new value.
    Behaves like an atomic load.

.. function:: numba.cuda.atomic.max(array, idx, value)

    Perform ``array[idx] = max(array[idx], value)``. Support int32, int64,
    float32 and float64 only. The ``idx`` argument can be an integer or a
    tuple of integer indices for indexing into multiple dimensional arrays.
    The number of element in ``idx`` must match the number of dimension of
    ``array``.

    Returns the value of ``array[idx]`` before the storing the new value.
    Behaves like an atomic load.


.. function:: numba.cuda.syncthreads

    Synchronize all threads in the same thread block.  This function implements
    the same pattern as barriers in traditional multi-threaded programming: this
    function waits until all threads in the block call it, at which point it
    returns control to all its callers.

.. function:: numba.cuda.syncthreads_count(predicate)

    An extension to :attr:`numba.cuda.syncthreads` where the return value is a count
    of the threads where ``predicate`` is true.

.. function:: numba.cuda.syncthreads_and(predicate)

    An extension to :attr:`numba.cuda.syncthreads` where 1 is returned if ``predicate`` is
    true for all threads or 0 otherwise.

.. function:: numba.cuda.syncthreads_or(predicate)

    An extension to :attr:`numba.cuda.syncthreads` where 1 is returned if ``predicate`` is
    true for any thread or 0 otherwise.

    .. warning:: All syncthreads functions must be called by every thread in the
                 thread-block. Falling to do so may result in undefined behavior.


Cooperative Groups
~~~~~~~~~~~~~~~~~~

.. function:: numba.cuda.cg.this_grid()

   Get the current grid group.

   :return: The current grid group
   :rtype: numba.cuda.cg.GridGroup

.. class:: numba.cuda.cg.GridGroup

   A grid group. Users should not construct a GridGroup directly - instead, get
   the current grid group using :func:`cg.this_grid() <numba.cuda.cg.this_grid>`.

   .. method:: sync()

      Synchronize the current grid group.


Memory Fences
~~~~~~~~~~~~~

The memory fences are used to guarantee the effect of memory operations
are visible by other threads within the same thread-block, the same GPU device,
and the same system (across GPUs on global memory). Memory loads and stores
are guaranteed to not move across the memory fences by optimization passes.

.. warning:: The memory fences are considered to be advanced API and most
             usercases should use the thread barrier (e.g. ``syncthreads()``).



.. function:: numba.cuda.threadfence

   A memory fence at device level (within the GPU).

.. function:: numba.cuda.threadfence_block

   A memory fence at thread block level.

.. function:: numba.cuda.threadfence_system


   A memory fence at system level (across GPUs).

Warp Intrinsics
~~~~~~~~~~~~~~~

The argument ``membermask`` is a 32 bit integer mask with each bit
corresponding to a thread in the warp, with 1 meaning the thread is in the
subset of threads within the function call. The ``membermask`` must be all 1 if
the GPU compute capability is below 7.x.

.. function:: numba.cuda.syncwarp(membermask)

   Synchronize a masked subset of the threads in a warp.

.. function:: numba.cuda.all_sync(membermask, predicate)

    If the ``predicate`` is true for all threads in the masked warp, then
    a non-zero value is returned, otherwise 0 is returned.

.. function:: numba.cuda.any_sync(membermask, predicate)

    If the ``predicate`` is true for any thread in the masked warp, then
    a non-zero value is returned, otherwise 0 is returned.

.. function:: numba.cuda.eq_sync(membermask, predicate)

    If the boolean ``predicate`` is the same for all threads in the masked warp,
    then a non-zero value is returned, otherwise 0 is returned.

.. function:: numba.cuda.ballot_sync(membermask, predicate)

    Returns a mask of all threads in the warp whose ``predicate`` is true,
    and are within the given mask.

.. function:: numba.cuda.shfl_sync(membermask, value, src_lane)

    Shuffles ``value`` across the masked warp and returns the ``value``
    from ``src_lane``. If this is outside the warp, then the
    given ``value`` is returned.

.. function:: numba.cuda.shfl_up_sync(membermask, value, delta)

    Shuffles ``value`` across the masked warp and returns the ``value``
    from ``laneid - delta``. If this is outside the warp, then the
    given ``value`` is returned.

.. function:: numba.cuda.shfl_down_sync(membermask, value, delta)

    Shuffles ``value`` across the masked warp and returns the ``value``
    from ``laneid + delta``. If this is outside the warp, then the
    given ``value`` is returned.

.. function:: numba.cuda.shfl_xor_sync(membermask, value, lane_mask)

    Shuffles ``value`` across the masked warp and returns the ``value``
    from ``laneid ^ lane_mask``.

.. function:: numba.cuda.match_any_sync(membermask, value, lane_mask)

    Returns a mask of threads that have same ``value`` as the given ``value``
    from within the masked warp.

.. function:: numba.cuda.match_all_sync(membermask, value, lane_mask)

    Returns a tuple of (mask, pred), where mask is a mask of threads that have
    same ``value`` as the given ``value`` from within the masked warp, if they
    all have the same value, otherwise it is 0. And pred is a boolean of whether
    or not all threads in the mask warp have the same warp.

.. function:: numba.cuda.activemask()

    Returns a 32-bit integer mask of all currently active threads in the
    calling warp. The Nth bit is set if the Nth lane in the warp is active when
    activemask() is called. Inactive threads are represented by 0 bits in the
    returned mask. Threads which have exited the kernel are always marked as
    inactive.

.. function:: numba.cuda.lanemask_lt()

    Returns a 32-bit integer mask of all lanes (including inactive ones) with
    ID less than the current lane.


Integer Intrinsics
~~~~~~~~~~~~~~~~~~

A subset of the CUDA Math API's integer intrinsics are available. For further
documentation, including semantics, please refer to the `CUDA Toolkit
documentation
<https://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__INTRINSIC__INT.html>`_.


.. function:: numba.cuda.popc(x)

   Returns the number of bits set in ``x``.

.. function:: numba.cuda.brev(x)

   Returns the reverse of the bit pattern of ``x``. For example, ``0b10110110``
   becomes ``0b01101101``.

.. function:: numba.cuda.clz(x)

   Returns the number of leading zeros in ``x``.

.. function:: numba.cuda.ffs(x)

   Returns the position of the first (least significant) bit set to 1 in ``x``,
   where the least significant bit position is 1. ``ffs(0)`` returns 0.


Floating Point Intrinsics
~~~~~~~~~~~~~~~~~~~~~~~~~

A subset of the CUDA Math API's floating point intrinsics are available. For further
documentation, including semantics, please refer to the `single
<https://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__SINGLE.html>`_ and
`double <https://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__DOUBLE.html>`_
precision parts of the CUDA Toolkit documentation.


.. function:: numba.cuda.fma

   Perform the fused multiply-add operation. Named after the ``fma`` and ``fmaf`` in
   the C api, but maps to the ``fma.rn.f32`` and ``fma.rn.f64`` (round-to-nearest-even)
   PTX instructions.

.. function:: numba.cuda.cbrt (x)

   Perform the cube root operation, x ** (1/3). Named after the functions
   ``cbrt`` and ``cbrtf`` in the C api. Supports float32, and float64 arguments
   only.

16-bit Floating Point Intrinsics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following functions are used to operate on 16-bit floating point operands.
These functions return a 16-bit floating point result.


.. function:: numba.cuda.fp16.hfma (a, b, c)

   Perform the fused multiply-add operation ``(a * b) + c`` on 16-bit
   floating point arguments in round to nearest mode. Maps to the ``fma.rn.f16``
   PTX instruction.

   Returns the 16-bit floating point result of the fused multiply-add.

.. function:: numba.cuda.fp16.hadd (a, b)

   Perform the add operation ``a + b`` on 16-bit floating point arguments in
   round to nearest mode. Maps to the ``add.f16`` PTX instruction.

   Returns the 16-bit floating point result of the addition.

.. function:: numba.cuda.fp16.hsub (a, b)

   Perform the subtract operation ``a - b`` on 16-bit floating point arguments in
   round to nearest mode. Maps to the ``sub.f16`` PTX instruction.

   Returns the 16-bit floating point result of the subtraction.

.. function:: numba.cuda.fp16.hmul (a, b)

   Perform the multiply operation ``a * b`` on 16-bit floating point arguments in
   round to nearest mode. Maps to the ``mul.f16`` PTX instruction.

   Returns the 16-bit floating point result of the multiplication.

.. function:: numba.cuda.fp16.hneg (a)

   Perform the negation operation ``-a`` on the 16-bit floating point argument.
   Maps to the ``neg.f16`` PTX instruction.

   Returns the 16-bit floating point result of the negation.

.. function:: numba.cuda.fp16.habs (a)

   Perform the absolute value operation ``|a|`` on the 16-bit floating point argument.
   Maps to the ``abs.f16`` PTX instruction.

   Returns the 16-bit floating point result of the absolute value operation.

Control Flow Instructions
~~~~~~~~~~~~~~~~~~~~~~~~~

A subset of the CUDA's control flow instructions are directly available as
intrinsics. Avoiding branches is a key way to improve CUDA performance, and
using these intrinsics mean you don't have to rely on the ``nvcc`` optimizer
identifying and removing branches. For further documentation, including
semantics, please refer to the `relevant CUDA Toolkit documentation
<https://docs.nvidia.com/cuda/parallel-thread-execution/index.html#comparison-and-selection-instructions>`_.


.. function:: numba.cuda.selp

    Select between two expressions, depending on the value of the first
    argument. Similar to LLVM's ``select`` instruction.


Timer Intrinsics
~~~~~~~~~~~~~~~~

.. function:: numba.cuda.nanosleep(ns)

    Suspends the thread for a sleep duration approximately close to the delay
    ``ns``, specified in nanoseconds.
CUDA Python Reference
=====================

.. toctree::

   host.rst
   kernel.rst
   memory.rst
   libdevice.rst
.. _deprecation:

===================
Deprecation Notices
===================

This section contains information about deprecation of behaviours, features and
APIs that have become undesirable/obsolete. Any information about the schedule
for their deprecation and reasoning behind the changes, along with examples, is
provided. However, first is a small section on how to suppress deprecation
warnings that may be raised from Numba so as to prevent warnings propagating
into code that is consuming Numba.

Suppressing Deprecation warnings
================================
All Numba deprecations are issued via ``NumbaDeprecationWarning`` or
``NumbaPendingDeprecationWarning`` s, to suppress the reporting of
these the following code snippet can be used::

    from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
    import warnings

    warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
    warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

The ``action`` used above is ``'ignore'``, other actions are available, see
`The Warnings Filter <https://docs.python.org/3/library/warnings.html#the-warnings-filter>`_
documentation for more information.

.. note:: It is **strongly recommended** that applications and libraries which
          choose to suppress these warnings should pin their Numba dependency
          to a suitable version because their users will no longer be aware of
          the coming incompatibility.

Deprecation of reflection for List and Set types
================================================
Reflection (:term:`reflection`) is the jargon used in Numba to describe the
process of ensuring that changes made by compiled code to arguments that are
mutable Python container data types are visible in the Python interpreter when
the compiled function returns. Numba has for some time supported reflection of
``list`` and ``set`` data types and it is support for this reflection that
is scheduled for deprecation with view to replace with a better implementation.

Reason for deprecation
----------------------
First recall that for Numba to be able to compile a function in ``nopython``
mode all the variables must have a concrete type ascertained through type
inference. In simple cases, it is clear how to reflect changes to containers
inside ``nopython`` mode back to the original Python containers. However,
reflecting changes to complex data structures with nested container types (for
example, lists of lists of integers) quickly becomes impossible to do
efficiently and consistently. After a number of years of experience with this
problem, it is clear that providing this behaviour is both fraught with
difficulty and often leads to code which does not have good performance (all
reflected data has to go through special APIs to convert the data to native
formats at call time and then back to CPython formats at return time). As a
result of this, the sheer number of reported problems in the issue tracker, and
how well a new approach that was taken with ``typed.Dict`` (typed dictionaries)
has gone, the core developers have decided to deprecate the noted ``reflection``
behaviour.


Example(s) of the impact
------------------------

At present only a warning of the upcoming change is issued. In future code such
as::

  from numba import njit

  @njit
  def foo(x):
      x.append(10)

  a = [1, 2, 3]
  foo(a)

will require adjustment to use a ``typed.List`` instance, this typed container
is synonymous to the :ref:`feature-typed-dict`. An example of translating the
above is::

    from numba import njit
    from numba.typed import List

    @njit
    def foo(x):
        x.append(10)

    a = [1, 2, 3]
    typed_a = List()
    [typed_a.append(x) for x in a]
    foo(typed_a)

For more information about ``typed.List`` see :ref:`feature-typed-list`. Further
usability enhancements for this feature were made in the 0.47.0 release
cycle.

Schedule
--------
This feature will be removed with respect to this schedule:

* Pending-deprecation warnings will be issued in version 0.44.0
* Prominent notice will be given for a minimum of two releases prior to full
  removal.

Recommendations
---------------
Projects that need/rely on the deprecated behaviour should pin their dependency
on Numba to a version prior to removal of this behaviour, or consider following
replacement instructions that will be issued outlining how to adjust to the
change.

Expected Replacement
--------------------
As noted above ``typed.List`` will be used to permit similar functionality to
reflection in the case of ``list`` s, a ``typed.Set`` will provide the
equivalent for ``set`` (not implemented yet!). The advantages to this approach
are:

* That the containers are typed means type inference has to work less hard.
* Nested containers (containers of containers of ...) are more easily
  supported.
* Performance penalties currently incurred translating data to/from native
  formats are largely avoided.
* Numba's ``typed.Dict`` will be able to use these containers as values.


Deprecation of :term:`object mode` `fall-back` behaviour when using ``@jit``
============================================================================
The ``numba.jit`` decorator has for a long time followed the behaviour of first
attempting to compile the decorated function in :term:`nopython mode` and should
this compilation fail it will `fall-back` and try again to compile but this time
in :term:`object mode`. It it this `fall-back` behaviour which is being
deprecated, the result of which will be that ``numba.jit`` will by default
compile in :term:`nopython mode` and :term:`object mode` compilation will
become `opt-in` only.


Reason for deprecation
----------------------
The `fall-back` has repeatedly caused confusion for users as seemingly innocuous
changes in user code can lead to drastic performance changes as code which may
have once compiled in :term:`nopython mode` mode may silently switch to
compiling in :term:`object mode` e.g::

    from numba import jit

    @jit
    def foo():
        l = []
        for x in range(10):
            l.append(x)
        return l

    foo()

    assert foo.nopython_signatures # this was compiled in nopython mode

    @jit
    def bar():
        l = []
        for x in range(10):
            l.append(x)
        return reversed(l) # innocuous change, but no reversed support in nopython mode

    bar()

    assert not bar.nopython_signatures # this was not compiled in nopython mode

Another reason to remove the `fall-back` is that it is confusing for the
compiler engineers developing Numba as it causes internal state problems that
are really hard to debug and it makes manipulating the compiler pipelines
incredibly challenging.

Further, it has long been considered best practice that the
:term:`nopython mode` keyword argument in the ``numba.jit`` decorator is set to
``True`` and that any user effort spent should go into making code work in this
mode as there's very little gain if it does not. The result is that, as Numba
has evolved, the amount of use :term:`object mode` gets in practice and its
general utility has decreased. It can be noted that there are some minor
improvements available through the notion of :term:`loop-lifting`, the cases of
this being used in practice are, however, rare and often a legacy from use of
less-recent Numba whereby such behaviour was better accommodated/the use of
``@jit`` with `fall-back` was recommended.


Example(s) of the impact
------------------------
At present a warning of the upcoming change is issued if ``@jit`` decorated code
uses the `fall-back` compilation path. In future code such as::

    @jit
    def bar():
        l = []
        for x in range(10):
            l.append(x)
        return reversed(l)

    bar()

will simply not compile, a ``TypingError`` would be raised.

Schedule
--------
This feature will be removed with respect to this schedule:

* Deprecation warnings will be issued in version 0.44.0
* Prominent notice will be given for a minimum of two releases prior to full
  removal.

Recommendations
---------------
Projects that need/rely on the deprecated behaviour should pin their dependency
on Numba to a version prior to removal of this behaviour. Alternatively, to
accommodate the scheduled deprecations, users with code compiled at present with
``@jit`` can supply the ``nopython=True`` keyword argument, if the code
continues to compile then the code is already ready for this change. If the code
does not compile, continue using the ``@jit`` decorator without
``nopython=True`` and profile the performance of the function. Then remove the
decorator and again check the performance of the function. If there is no
benefit to having the ``@jit`` decorator present consider removing it! If there
is benefit to having the ``@jit`` decorator present, then to be future proof
supply the keyword argument ``forceobj=True`` to ensure the function is always
compiled in :term:`object mode`.


.. _deprecation-strict-strides:


Deprecation of the ``inspect_ptx()`` method
===========================================

The undocumented ``inspect_ptx()`` method of functions decorated with
``@cuda.jit(device=True)`` is sometimes used to compile a Python function to
PTX for use outside of Numba. An interface for this specific purpose is
provided in the :func:`compile_ptx() <numba.cuda.compile_ptx>` function.
``inspect_ptx()`` has one or two longstanding issues and presents a maintenance
burden for upcoming changes in the CUDA target, so it is deprecated and will be
removed in favor of the use of :func:`compile_ptx() <numba.cuda.compile_ptx>`.

Recommendations
---------------

Replace any code that compiles device functions to PTX using the following
pattern:

.. code-block:: python

    @cuda.jit(signature, device=True)
    def func(args):
        ...

    ptx_code = func.inspect_ptx(nvvm_options=nvvm_options).decode()

with:

.. code-block:: python

    def func(args):
        ...

    ptx_code, return_type = compile_ptx(func, signature, device=True, nvvm_options=nvvm_options)

Schedule
--------

- In Numba 0.54: ``inspect_ptx()`` was deprecated.
- In Numba 0.55: ``inspect_ptx()`` was removed.


Deprecation of eager compilation of CUDA device functions
=========================================================

In future versions of Numba, the ``device`` kwarg to the ``@cuda.jit`` decorator
will be obviated, and whether a device function or global kernel is compiled will
be inferred from the context. With respect to kernel / device functions and lazy
/ eager compilation, four cases were handled:

1. ``device=True``, eager compilation with a signature provided
2. ``device=False``, eager compilation with a signature provided
3. ``device=True``, lazy compilation with no signature
4. ``device=False``, lazy compilation with no signature

The latter two cases can be differentiated without the ``device`` kwarg, because
it can be inferred from the calling context - if the call is from the host, then
a global kernel should be compiled, and if the call is from a kernel or another
device function, then a device function should be compiled.

The first two cases cannot be differentiated in the absence of the ``device``
kwarg - without it, it will not be clear from a signature alone whether a device
function or global kernel should be compiled. In order to resolve this, device
functions will no longer be eagerly compiled. When a signature is provided to a
device function, it will only be used to enforce the types of arguments that
the function accepts.

.. note::

   In previous releases this notice stated that support for providing
   signatures to device functions would be removed completely - however, this
   precludes the common use case of enforcing the types that can be passed to a
   device function (and the automatic insertion of casts that it implies) so
   this notice has been updated to retain support for passing signatures.


Schedule
--------

- In Numba 0.54: Eager compilation of device functions will be deprecated.
- In Numba 0.55: Eager compilation of device functions will be unsupported and
  the provision of signatures for device functions will only enforce casting.


Deprecation of ``numba.core.base.BaseContext.add_user_function()``
==================================================================

``add_user_function()``  offers the same functionality as
``insert_user_function()``, only with a check that the function has already
been inserted at least once.  It is now deprecated as it is no longer used
internally and it is expected that it is not used externally.

Recommendations
---------------

Replace any uses of ``add_user_function()`` with ``insert_user_function()``.

Schedule
--------

- In Numba 0.55: ``add_user_function()`` will be deprecated.
- In Numba 0.56: ``add_user_function()`` will be removed.


Deprecation of CUDA Toolkits < 10.2 and devices with CC < 5.3
=============================================================

Support for:

- Devices with Compute Capability < 5.3, and
- CUDA toolkits less than 10.2

is deprecated and will be removed in future.

Recommendations
---------------

- For devices of Compute Capability 3.0 - 5.2, Numba 0.55.1 or earlier will be
  required.
- CUDA toolkit 10.2 or later (ideally 11.2 or later) should be installed.

Schedule
--------

- In Numba 0.55.1: support for CC < 5.3 and CUDA toolkits < 10.2 are deprecated.
- In Numba 0.56: support for CC < 5.3 and CUDA toolkits < 10.2 will be removed.
Just-in-Time compilation
========================


JIT functions
-------------

.. _jit-decorator:

.. decorator:: numba.jit(signature=None, nopython=False, nogil=False, cache=False, forceobj=False, parallel=False, error_model='python', fastmath=False, locals={}, boundscheck=False)

   Compile the decorated function on-the-fly to produce efficient machine
   code.  All parameters are optional.

   If present, the *signature* is either a single signature or a list of
   signatures representing the expected :ref:`numba-types` of function
   arguments and return values.  Each signature can be given in several
   forms:

   * A tuple of :ref:`numba-types` arguments (for example
     ``(numba.int32, numba.double)``) representing the types of the
     function's arguments; Numba will then infer an appropriate return
     type from the arguments.
   * A call signature using :ref:`numba-types`, specifying both return
     type and argument types. This can be given in intuitive form
     (for example ``numba.void(numba.int32, numba.double)``).
   * A string representation of one of the above, for example
     ``"void(int32, double)"``.  All type names used in the string are assumed
     to be defined in the ``numba.types`` module.

   *nopython* and *nogil* are boolean flags.  *locals* is a mapping of
   local variable names to :ref:`numba-types`.

   This decorator has several modes of operation:

   * If one or more signatures are given in *signature*, a specialization is
     compiled for each of them.  Calling the decorated function will then try
     to choose the best matching signature, and raise a :class:`TypeError` if
     no appropriate conversion is available for the function arguments.  If
     converting succeeds, the compiled machine code is executed with the
     converted arguments and the return value is converted back according to
     the signature.

   * If no *signature* is given, the decorated function implements
     lazy compilation.  Each call to the decorated function will try to
     re-use an existing specialization if it exists (for example, a call
     with two integer arguments may re-use a specialization for argument
     types ``(numba.int64, numba.int64)``).  If no suitable specialization
     exists, a new specialization is compiled on-the-fly, stored for later
     use, and executed with the converted arguments.

   If true, *nopython* forces the function to be compiled in :term:`nopython
   mode`. If not possible, compilation will raise an error.

   If true, *forceobj* forces the function to be compiled in :term:`object
   mode`.  Since object mode is slower than nopython mode, this is mostly
   useful for testing purposes.

   If true, *nogil* tries to release the :py:term:`global interpreter lock`
   inside the compiled function.  The GIL will only be released if Numba can
   compile the function in :term:`nopython mode`, otherwise a compilation
   warning will be printed.

   .. _jit-decorator-cache:

   If true, *cache* enables a file-based cache to shorten compilation times
   when the function was already compiled in a previous invocation.
   The cache is maintained in the ``__pycache__`` subdirectory of
   the directory containing the source file; if the current user is not
   allowed to write to it, though, it falls back to a platform-specific
   user-wide cache directory (such as ``$HOME/.cache/numba`` on Unix
   platforms).

   .. _jit-decorator-parallel:

   If true, *parallel* enables the automatic parallelization of a number of
   common Numpy constructs as well as the fusion of adjacent parallel
   operations to maximize cache locality.

   The *error_model* option controls the divide-by-zero behavior.
   Setting it to 'python' causes divide-by-zero to raise exception like CPython.
   Setting it to 'numpy' causes divide-by-zero to set the result to *+/-inf* or
   *nan*.

   Not all functions can be cached, since some functionality cannot be
   always persisted to disk.  When a function cannot be cached, a
   warning is emitted.

   .. _jit-decorator-fastmath:

   If true, *fastmath* enables the use of otherwise unsafe floating point
   transforms as described in the
   `LLVM documentation <https://llvm.org/docs/LangRef.html#fast-math-flags>`_.
   Further, if :ref:`Intel SVML <intel-svml>` is installed faster but less
   accurate versions of some math intrinsics are used (answers to within
   ``4 ULP``).

   .. _jit-decorator-boundscheck:

   If True, ``boundscheck`` enables bounds checking for array indices. Out of
   bounds accesses will raise IndexError. The default is to not do bounds
   checking. If bounds checking is disabled, out of bounds accesses can
   produce garbage results or segfaults. However, enabling bounds checking
   will slow down typical functions, so it is recommended to only use this
   flag for debugging. You can also set the `NUMBA_BOUNDSCHECK` environment
   variable to 0 or 1 to globally override this flag.

   The *locals* dictionary may be used to force the :ref:`numba-types`
   of particular local variables, for example if you want to force the
   use of single precision floats at some point.  In general, we recommend
   you let Numba's compiler infer the types of local variables by itself.

   Here is an example with two signatures::

      @jit(["int32(int32)", "float32(float32)"], nopython=True)
      def f(x): ...

   Not putting any parentheses after the decorator is equivalent to calling
   the decorator without any arguments, i.e.::

      @jit
      def f(x): ...

   is equivalent to::

      @jit()
      def f(x): ...

   The decorator returns a :class:`Dispatcher` object.

   .. note::
      If no *signature* is given, compilation errors will be raised when
      the actual compilation occurs, i.e. when the function is first called
      with some given argument types.

   .. note::
      Compilation can be influenced by some dedicated :ref:`numba-envvars`.


Generated JIT functions
-----------------------

.. decorator:: numba.generated_jit(nopython=False, nogil=False, cache=False, forceobj=False, locals={})

   Like the :func:`~numba.jit` decorator, but calls the decorated function at
   compile-time, passing the *types* of the function's arguments.
   The decorated function must return a callable which will be compiled as
   the function's implementation for those types, allowing flexible kinds of
   specialization.

   The :func:`~numba.generated_jit` decorator returns a :class:`Dispatcher` object.


Dispatcher objects
------------------

.. class:: Dispatcher

   The class of objects created by calling :func:`~numba.jit` or
   :func:`~numba.generated_jit`.  You shouldn't try to create such an object
   in any other way.  Calling a Dispatcher object calls the compiled
   specialization for the arguments with which it is called, letting it
   act as an accelerated replacement for the Python function which was compiled.

   In addition, Dispatcher objects have the following methods and attributes:

   .. attribute:: py_func

      The pure Python function which was compiled.

   .. method:: inspect_types(file=None, pretty=False)

      Print out a listing of the function source code annotated line-by-line
      with the corresponding Numba IR, and the inferred types of the various
      variables.  If *file* is specified, printing is done to that file
      object, otherwise to sys.stdout. If *pretty* is set to True then colored
      ANSI will be produced in a terminal and HTML in a notebook.

      .. seealso:: :ref:`architecture`

   .. method:: inspect_llvm(signature=None)

      Return a dictionary keying compiled function signatures to the human
      readable LLVM IR generated for the function.  If the signature
      keyword is specified a string corresponding to that individual
      signature is returned.

   .. method:: inspect_asm(signature=None)

      Return a dictionary keying compiled function signatures to the
      human-readable native assembly code for the function.  If the
      signature keyword is specified a string corresponding to that
      individual signature is returned.

   .. method:: inspect_cfg(signature=None, show_wrapped)

      Return a dictionary keying compiled function signatures to the
      control-flow graph objects for the function.  If the signature keyword is
      specified a string corresponding to that individual signature is returned.

      The control-flow graph objects can be stringified (``str`` or ``repr``)
      to get the textual representation of the graph in DOT format.  Or, use
      its ``.display(filename=None, view=False)`` method to plot the graph.
      The *filename* option can be set to a specific path for the rendered
      output to write to.  If *view* option is True, the plot is opened by
      the system default application for the image format (PDF). In IPython
      notebook, the returned object can be plot inlined.

      Usage::

        @jit
        def foo():
          ...

        # opens the CFG in system default application
        foo.inspect_cfg(foo.signatures[0]).display(view=True)


   .. method:: inspect_disasm_cfg(signature=None)

      Return a dictionary keying compiled function signatures to the
      control-flow graph of the disassembly of the underlying compiled ``ELF``
      object.  If the signature keyword is specified a control-flow graph
      corresponding to that individual signature is returned. This function is
      execution environment aware and will produce SVG output in Jupyter
      notebooks and ASCII in terminals.

      Example::

        @njit
        def foo(x):
            if x < 3:
                return x + 1
            return x + 2

        foo(10)

        print(foo.inspect_disasm_cfg(signature=foo.signatures[0]))

      Gives::

        [0x08000040]>  # method.__main__.foo_241_long_long (int64_t arg1, int64_t arg3);
         ─────────────────────────────────────────────────────────────────────┐
        │  0x8000040                                                          │
        │ ; arg3 ; [02] -r-x section size 279 named .text                     │
        │   ;-- section..text:                                                │
        │   ;-- .text:                                                        │
        │   ;-- __main__::foo$241(long long):                                 │
        │   ;-- rip:                                                          │
        │ 25: method.__main__.foo_241_long_long (int64_t arg1, int64_t arg3); │
        │ ; arg int64_t arg1 @ rdi                                            │
        │ ; arg int64_t arg3 @ rdx                                            │
        │ ; 2                                                                 │
        │ cmp rdx, 2                                                          │
        │ jg 0x800004f                                                        │
        └─────────────────────────────────────────────────────────────────────┘
                f t
                │ │
                │ └──────────────────────────────┐
                └──┐                             │
                   │                             │
            ┌─────────────────────────┐   ┌─────────────────────────┐
            │  0x8000046              │   │  0x800004f              │
            │ ; arg3                  │   │ ; arg3                  │
            │ inc rdx                 │   │ add rdx, 2              │
            │ ; arg3                  │   │ ; arg3                  │
            │ mov qword [rdi], rdx    │   │ mov qword [rdi], rdx    │
            │ xor eax, eax            │   │ xor eax, eax            │
            │ ret                     │   │ ret                     │
            └─────────────────────────┘   └─────────────────────────┘

   .. method:: recompile()

      Recompile all existing signatures.  This can be useful for example if
      a global or closure variable was frozen by your function and its value
      in Python has changed.  Since compiling isn't cheap, this is mainly
      for testing and interactive use.

   .. method:: parallel_diagnostics(signature=None, level=1)

      Print parallel diagnostic information for the given signature. If no
      signature is present it is printed for all known signatures. ``level`` is
      used to adjust the verbosity, ``level=1`` (default) is minimum verbosity,
      levels 2, 3, and 4 provide increasing levels of verbosity.

   .. method:: get_metadata(signature=None)

      Obtain the compilation metadata for a given signature. This is useful for
      developers of Numba and Numba extensions.


Vectorized functions (ufuncs and DUFuncs)
-----------------------------------------

.. decorator:: numba.vectorize(*, signatures=[], identity=None, nopython=True, target='cpu', forceobj=False, cache=False, locals={})

   Compile the decorated function and wrap it either as a `Numpy
   ufunc`_ or a Numba :class:`~numba.DUFunc`.  The optional
   *nopython*, *forceobj* and *locals* arguments have the same meaning
   as in :func:`numba.jit`.

   *signatures* is an optional list of signatures expressed in the
   same form as in the :func:`numba.jit` *signature* argument.  If
   *signatures* is non-empty, then the decorator will compile the user
   Python function into a Numpy ufunc.  If no *signatures* are given,
   then the decorator will wrap the user Python function in a
   :class:`~numba.DUFunc` instance, which will compile the user
   function at call time whenever Numpy can not find a matching loop
   for the input arguments.  *signatures* is required if *target* is
   ``"parallel"``.

   *identity* is the identity (or unit) value of the function being
   implemented.  Possible values are 0, 1, None, and the string
   ``"reorderable"``.  The default is None.  Both None and
   ``"reorderable"`` mean the function has no identity value;
   ``"reorderable"`` additionally specifies that reductions along multiple
   axes can be reordered.

   If there are several *signatures*, they must be ordered from the more
   specific to the least specific.  Otherwise, Numpy's type-based
   dispatching may not work as expected.  For example, the following is
   wrong::

      @vectorize(["float64(float64)", "float32(float32)"])
      def f(x): ...

   as running it over a single-precision array will choose the ``float64``
   version of the compiled function, leading to much less efficient
   execution.  The correct invocation is::

      @vectorize(["float32(float32)", "float64(float64)"])
      def f(x): ...

   *target* is a string for backend target; Available values are "cpu",
   "parallel", and "cuda".  To use a multithreaded version, change the
   target to "parallel" (which requires signatures to be specified)::

      @vectorize(["float64(float64)", "float32(float32)"], target='parallel')
      def f(x): ...

   For the CUDA target, use "cuda"::

      @vectorize(["float64(float64)", "float32(float32)"], target='cuda')
      def f(x): ...

   The compiled function can be cached to reduce future compilation time.
   It is enabled by setting *cache* to True. Only the "cpu" and "parallel"
   targets support caching.


.. decorator:: numba.guvectorize(signatures, layout, *, identity=None, nopython=True, target='cpu', forceobj=False, cache=False, locals={})

   Generalized version of :func:`numba.vectorize`.  While
   :func:`numba.vectorize` will produce a simple ufunc whose core
   functionality (the function you are decorating) operates on scalar
   operands and returns a scalar value, :func:`numba.guvectorize`
   allows you to create a `Numpy ufunc`_ whose core function takes array
   arguments of various dimensions.

   The additional argument *layout* is a string specifying, in symbolic
   form, the dimensionality and size relationship of the argument types
   and return types.  For example, a matrix multiplication will have
   a layout string of ``"(m,n),(n,p)->(m,p)"``.  Its definition might
   be (function body omitted)::

      @guvectorize(["void(float64[:,:], float64[:,:], float64[:,:])"],
                   "(m,n),(n,p)->(m,p)")
      def f(a, b, result):
          """Fill-in *result* matrix such as result := a * b"""
          ...

   If one of the arguments should be a scalar, the corresponding layout
   specification is ``()`` and the argument will really be given to
   you as a zero-dimension array (you have to dereference it to get the
   scalar value).  For example, a :ref:`one-dimension moving average <example-movemean>`
   with a parameterable window width may have a layout string of ``"(n),()->(n)"``.

   Note that any output will be given to you preallocated as an additional
   function argument: your code has to fill it with the appropriate values
   for the function you are implementing.

   If your function doesn't take an output array, you should omit the "arrow"
   in the layout string (e.g. ``"(n),(n)"``). When doing this, it is important
   to be aware that changes to the input arrays cannot always be relied on to be
   visible outside the execution of the ufunc, as NumPy may pass in temporary
   arrays as inputs (for example, if a cast is required).

   .. seealso::
      Specification of the `layout string <https://numpy.org/doc/stable/reference/c-api/generalized-ufuncs.html#details-of-signature>`_
      as supported by Numpy.  Note that Numpy uses the term "signature",
      which we unfortunately use for something else.

   The compiled function can be cached to reduce future compilation time.
   It is enabled by setting *cache* to True. Only the "cpu" and "parallel"
   targets support caching.

.. _Numpy ufunc: http://docs.scipy.org/doc/numpy/reference/ufuncs.html

.. class:: numba.DUFunc

   The class of objects created by calling :func:`numba.vectorize`
   with no signatures.

   DUFunc instances should behave similarly to Numpy
   :class:`~numpy.ufunc` objects with one important difference:
   call-time loop generation.  When calling a ufunc, Numpy looks at
   the existing loops registered for that ufunc, and will raise a
   :class:`~python.TypeError` if it cannot find a loop that it cannot
   safely cast the inputs to suit.  When calling a DUFunc, Numba
   delegates the call to Numpy.  If the Numpy ufunc call fails, then
   Numba attempts to build a new loop for the given input types, and
   calls the ufunc again.  If this second call attempt fails or a
   compilation error occurs, then DUFunc passes along the exception to
   the caller.

   .. seealso::

      The ":ref:`dynamic-universal-functions`" section in the user's
      guide demonstrates the call-time behavior of
      :class:`~numba.DUFunc`, and discusses the impact of call order
      on how Numba generates the underlying :class:`~numpy.ufunc`.

   .. attribute:: ufunc

      The actual Numpy :class:`~numpy.ufunc` object being built by the
      :class:`~numba.DUFunc` instance.  Note that the
      :class:`~numba.DUFunc` object maintains several important data
      structures required for proper ufunc functionality (specifically
      the dynamically compiled loops).  Users should not pass the
      :class:`~numpy.ufunc` value around without ensuring the
      underlying :class:`~numba.DUFunc` will not be garbage collected.

   .. attribute:: nin

      The number of DUFunc (ufunc) inputs.  See `ufunc.nin`_.

   .. attribute:: nout

      The number of DUFunc outputs.  See `ufunc.nout`_.

   .. attribute:: nargs

      The total number of possible DUFunc arguments (should be
      :attr:`~numba.DUFunc.nin` + :attr:`~numba.DUFunc.nout`).
      See `ufunc.nargs`_.

   .. attribute:: ntypes

      The number of input types supported by the DUFunc.  See
      `ufunc.ntypes`_.

   .. attribute:: types

      A list of the supported types given as strings.  See
      `ufunc.types`_.

   .. attribute:: identity

      The identity value when using the ufunc as a reduction.  See
      `ufunc.identity`_.

   .. method:: reduce(A, *, axis, dtype, out, keepdims)

      Reduces *A*\'s dimension by one by applying the DUFunc along one
      axis.  See `ufunc.reduce`_.

   .. method:: accumulate(A, *, axis, dtype, out)

      Accumulate the result of applying the operator to all elements.
      See `ufunc.accumulate`_.

   .. method:: reduceat(A, indices, *, axis, dtype, out)

      Performs a (local) reduce with specified slices over a single
      axis.  See `ufunc.reduceat`_.

   .. method:: outer(A, B)

      Apply the ufunc to all pairs (*a*, *b*) with *a* in *A*, and *b*
      in *B*.  See `ufunc.outer`_.

   .. method:: at(A, indices, *, B)

      Performs unbuffered in place operation on operand *A* for
      elements specified by *indices*.  If you are using Numpy 1.7 or
      earlier, this method will not be present.  See `ufunc.at`_.


.. note::
   Vectorized functions can, in rare circumstances, show
   :ref:`unexpected warnings or errors <ufunc-fpu-errors>`.


.. _`ufunc.nin`: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ufunc.nin.html#numpy.ufunc.nin

.. _`ufunc.nout`: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ufunc.nout.html#numpy.ufunc.nout

.. _`ufunc.nargs`: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ufunc.nargs.html#numpy.ufunc.nargs

.. _`ufunc.ntypes`: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ufunc.ntypes.html#numpy.ufunc.ntypes

.. _`ufunc.types`: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ufunc.types.html#numpy.ufunc.types

.. _`ufunc.identity`: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ufunc.identity.html#numpy.ufunc.identity

.. _`ufunc.reduce`: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ufunc.reduce.html#numpy.ufunc.reduce

.. _`ufunc.accumulate`: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ufunc.accumulate.html#numpy.ufunc.accumulate

.. _`ufunc.reduceat`: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ufunc.reduceat.html#numpy.ufunc.reduceat

.. _`ufunc.outer`: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ufunc.outer.html#numpy.ufunc.outer

.. _`ufunc.at`: http://docs.scipy.org/doc/numpy/reference/generated/numpy.ufunc.at.html#numpy.ufunc.at


C callbacks
-----------

.. decorator:: numba.cfunc(signature, nopython=False, cache=False, locals={})

   Compile the decorated function on-the-fly to produce efficient machine
   code.  The compiled code is wrapped in a thin C callback that makes it
   callable using the natural C ABI.

   The *signature* is a single signature representing the signature of the
   C callback.  It must have the same form as in :func:`~numba.jit`.
   The decorator does not check that the types in the signature have
   a well-defined representation in C.

   *nopython* and *cache* are boolean flags.  *locals* is a mapping of
   local variable names to :ref:`numba-types`.  They all have the same
   meaning as in :func:`~numba.jit`.

   The decorator returns a :class:`CFunc` object.

   .. note::
      C callbacks currently do not support :term:`object mode`.


.. class:: CFunc

   The class of objects created by :func:`~numba.cfunc`.  :class:`CFunc`
   objects expose the following attributes and methods:

   .. attribute:: address

      The address of the compiled C callback, as an integer.

   .. attribute:: cffi

      A `cffi`_ function pointer instance, to be passed as an argument to
      `cffi`_-wrapped functions.  The pointer's type is ``void *``, so
      only minimal type checking will happen when passing it to `cffi`_.

   .. attribute:: ctypes

      A :mod:`ctypes` callback instance, as if it were created using
      :func:`ctypes.CFUNCTYPE`.

   .. attribute:: native_name

      The name of the compiled C callback.

   .. method:: inspect_llvm()

      Return the human-readable LLVM IR generated for the C callback.
      :attr:`native_name` is the name under which this callback is defined
      in the IR.


.. _cffi: https://cffi.readthedocs.org/

=========
Utilities
=========

Dealing with pointers
=====================

These functions can be called from pure Python as well as in
:term:`nopython mode`.


.. function:: numba.carray(ptr, shape, dtype=None)

   Return a Numpy array view over the data pointed to by *ptr* with the
   given *shape*, in C order.  If *dtype* is given, it is used as the array's
   dtype, otherwise the array's dtype is inferred from *ptr*'s type.
   As the returned array is a view, not a copy, writing to it will modify
   the original data.

   *ptr* should be a ctypes pointer object (either a typed pointer
   as created using :func:`~ctypes.POINTER`, or a :class:`~ctypes.c_void_p`).

   *shape* should be an integer or a tuple of integers.

   *dtype* should be a Numpy dtype or scalar class (i.e. both
   ``np.dtype('int8')`` and ``np.int8`` are accepted).


.. function:: numba.farray(ptr, shape, dtype=None)

   Same as :func:`~numba.carray`, but the data is assumed to be laid out
   in Fortran order, and the array view is constructed accordingly.


.. _numpy-support:

========================
Supported NumPy features
========================

One objective of Numba is having a seamless integration with `NumPy`_.
NumPy arrays provide an efficient storage method for homogeneous sets of
data.  NumPy dtypes provide type information useful when compiling, and
the regular, structured storage of potentially large amounts of data
in memory provides an ideal memory layout for code generation.  Numba
excels at generating code that executes on top of NumPy arrays.

NumPy support in Numba comes in many forms:

* Numba understands calls to NumPy `ufuncs`_ and is able to generate
  equivalent native code for many of them.

* NumPy arrays are directly supported in Numba.  Access to Numpy arrays
  is very efficient, as indexing is lowered to direct memory accesses
  when possible.

* Numba is able to generate `ufuncs`_ and `gufuncs`_. This means that it
  is possible to implement ufuncs and gufuncs within Python, getting
  speeds comparable to that of ufuncs/gufuncs implemented in C extension
  modules using the NumPy C API.

.. _NumPy: http://www.numpy.org/
.. _ufuncs: http://docs.scipy.org/doc/numpy/reference/ufuncs.html
.. _gufuncs: http://docs.scipy.org/doc/numpy/reference/c-api.generalized-ufuncs.html

The following sections focus on the Numpy features supported in
:term:`nopython mode`, unless otherwise stated.


Scalar types
============

Numba supports the following Numpy scalar types:

* **Integers**: all integers of either signedness, and any width up to 64 bits
* **Booleans**
* **Real numbers:** single-precision (32-bit) and double-precision (64-bit) reals
* **Complex numbers:** single-precision (2x32-bit) and double-precision (2x64-bit) complex numbers
* **Datetimes and timestamps:** of any unit
* **Character sequences** (but no operations are available on them)
* **Structured scalars:** structured scalars made of any of the types above and arrays of the types above

The following scalar types and features are not supported:

* **Arbitrary Python objects**
* **Half-precision and extended-precision** real and complex numbers
* **Nested structured scalars** the fields of structured scalars may not contain other structured scalars

The operations supported on NumPy scalars are almost the same as on the
equivalent built-in types such as ``int`` or ``float``.  You can use a type's
constructor to convert from a different type or width. In addition you can use
the ``view(np.<dtype>)`` method to bitcast all ``int`` and ``float`` types
within the same width. However, you must define the scalar using a NumPy
constructor within a jitted function. For example, the following will work:

.. code:: pycon

    >>> import numpy as np
    >>> from numba import njit
    >>> @njit
    ... def bitcast():
    ...     i = np.int64(-1)
    ...     print(i.view(np.uint64))
    ...
    >>> bitcast()
    18446744073709551615


Whereas the following will not work:


.. code:: pycon

    >>> import numpy as np
    >>> from numba import njit
    >>> @njit
    ... def bitcast(i):
    ...     print(i.view(np.uint64))
    ...
    >>> bitcast(np.int64(-1))
    ---------------------------------------------------------------------------
    TypingError                               Traceback (most recent call last)
        ...
    TypingError: Failed in nopython mode pipeline (step: ensure IR is legal prior to lowering)
    'view' can only be called on NumPy dtypes, try wrapping the variable with 'np.<dtype>()'

    File "<ipython-input-3-fc40aaab84c4>", line 3:
    def bitcast(i):
        print(i.view(np.uint64))

Structured scalars support attribute getting and setting, as well as
member lookup using constant strings. Strings stored in a local or global tuple
are considered constant strings and can be used for member lookup.



.. literalinclude:: ../../../numba/tests/doc_examples/test_rec_array.py
   :language: python
   :start-after: magictoken.ex_rec_arr_const_index.begin
   :end-before: magictoken.ex_rec_arr_const_index.end
   :dedent: 8

It is also possible to use local or global tuples together with ``literal_unroll``:

.. literalinclude:: ../../../numba/tests/doc_examples/test_rec_array.py
   :language: python
   :start-after: magictoken.ex_rec_arr_lit_unroll_index.begin
   :end-before: magictoken.ex_rec_arr_lit_unroll_index.end
   :dedent: 8


Record subtyping
----------------
.. warning::
   This is an experimental feature.

Numba allows `width subtyping <https://en.wikipedia.org/wiki/Subtyping#Record_types>`_ of structured scalars.
For example, ``dtype([('a', 'f8'), ('b', 'i8')])`` will be considered a subtype of ``dtype([('a', 'f8')]``, because
the second is a strict subset of the first, i.e. field ``a`` is of the same type and is in the same position in both
types. The subtyping relationship will matter in cases where compilation for a certain input is not allowed, but the
input is a subtype of another, allowed type.

.. code-block:: python

    import numpy as np
    from numba import njit, typeof
    from numba.core import types
    record1 = np.array([1], dtype=[('a', 'f8')])[0]
    record2 = np.array([(2,3)], dtype=[('a', 'f8'), ('b', 'f8')])[0]

    @njit(types.float64(typeof(record1)))
    def foo(rec):
        return rec['a']

    foo(record1)
    foo(record2)

Without subtyping the last line would fail. With subtyping, no new compilation will be triggered, but the
compiled function for ``record1`` will be used for ``record2``.

.. seealso::
   `Numpy scalars <http://docs.scipy.org/doc/numpy/reference/arrays.scalars.html>`_
   reference.


Array types
===========

`Numpy arrays <http://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html>`_
of any of the scalar types above are supported, regardless of the shape
or layout.

Array access
------------

Arrays support normal iteration.  Full basic indexing and slicing is
supported.  A subset of advanced indexing is also supported: only one
advanced index is allowed, and it has to be a one-dimensional array
(it can be combined with an arbitrary number of basic indices as well).

.. seealso::
   `Numpy indexing <http://docs.scipy.org/doc/numpy/reference/arrays.indexing.html>`_
   reference.


.. _structured-array-access:

Structured array access
-----------------------

Numba presently supports accessing fields of individual elements in structured
arrays by attribute as well as by getting and setting. This goes slightly
beyond the NumPy API, which only allows accessing fields by getting and
setting. For example:

.. code:: python

   from numba import njit
   import numpy as np

   record_type = np.dtype([("ival", np.int32), ("fval", np.float64)], align=True)

   def f(rec):
       value = 2.5
       rec[0].ival = int(value)
       rec[0].fval = value
       return rec

   arr = np.ones(1, dtype=record_type)

   cfunc = njit(f)

   # Works
   print(cfunc(arr))

   # Does not work
   print(f(arr))

The above code results in the output:

.. code:: none

   [(2, 2.5)]
   Traceback (most recent call last):
     File "repro.py", line 22, in <module>
       print(f(arr))
     File "repro.py", line 9, in f
       rec[0].ival = int(value)
   AttributeError: 'numpy.void' object has no attribute 'ival'

The Numba-compiled version of the function executes, but the pure Python
version raises an error because of the unsupported use of attribute access.

.. note::
   This behavior will eventually be deprecated and removed.

Attributes
----------

The following attributes of Numpy arrays are supported:

* :attr:`~numpy.ndarray.dtype`
* :attr:`~numpy.ndarray.flags`
* :attr:`~numpy.ndarray.flat`
* :attr:`~numpy.ndarray.itemsize`
* :attr:`~numpy.ndarray.ndim`
* :attr:`~numpy.ndarray.shape`
* :attr:`~numpy.ndarray.size`
* :attr:`~numpy.ndarray.strides`
* :attr:`~numpy.ndarray.T`
* :attr:`~numpy.ndarray.real`
* :attr:`~numpy.ndarray.imag`

The ``flags`` object
''''''''''''''''''''

The object returned by the :attr:`~numpy.ndarray.flags` attribute supports
the ``contiguous``, ``c_contiguous`` and ``f_contiguous`` attributes.

The ``flat`` object
'''''''''''''''''''

The object returned by the :attr:`~numpy.ndarray.flat` attribute supports
iteration and indexing, but be careful: indexing is very slow on
non-C-contiguous arrays.

The ``real`` and ``imag`` attributes
''''''''''''''''''''''''''''''''''''

Numpy supports these attributes regardless of the dtype but Numba chooses to
limit their support to avoid potential user error.  For numeric dtypes,
Numba follows Numpy's behavior.  The :attr:`~numpy.ndarray.real` attribute
returns a view of the real part of the complex array and it behaves as an identity
function for other numeric dtypes.  The :attr:`~numpy.ndarray.imag` attribute
returns a view of the imaginary part of the complex array and it returns a zero
array with the same shape and dtype for other numeric dtypes.  For non-numeric
dtypes, including all structured/record dtypes, using these attributes will
result in a compile-time (`TypingError`) error.  This behavior differs from
Numpy's but it is chosen to avoid the potential confusion with field names that
overlap these attributes.

Calculation
-----------

The following methods of Numpy arrays are supported in their basic form
(without any optional arguments):

* :meth:`~numpy.ndarray.all`
* :meth:`~numpy.ndarray.any`
* :meth:`~numpy.ndarray.clip`
* :meth:`~numpy.ndarray.conj`
* :meth:`~numpy.ndarray.conjugate`
* :meth:`~numpy.ndarray.cumprod`
* :meth:`~numpy.ndarray.cumsum`
* :meth:`~numpy.ndarray.max`
* :meth:`~numpy.ndarray.mean`
* :meth:`~numpy.ndarray.min`
* :meth:`~numpy.ndarray.nonzero`
* :meth:`~numpy.ndarray.prod`
* :meth:`~numpy.ndarray.std`
* :meth:`~numpy.ndarray.take`
* :meth:`~numpy.ndarray.var`

The corresponding top-level Numpy functions (such as :func:`numpy.prod`)
are similarly supported.

Other methods
-------------

The following methods of Numpy arrays are supported:

* :meth:`~numpy.ndarray.argmax` (``axis`` keyword argument supported).
* :meth:`~numpy.ndarray.argmin` (``axis`` keyword argument supported).
* :meth:`~numpy.ndarray.argsort` (``kind`` key word argument supported for
  values ``'quicksort'`` and ``'mergesort'``)
* :meth:`~numpy.ndarray.astype` (only the 1-argument form)
* :meth:`~numpy.ndarray.copy` (without arguments)
* :meth:`~numpy.ndarray.dot` (only the 1-argument form)
* :meth:`~numpy.ndarray.flatten` (no order argument; 'C' order only)
* :meth:`~numpy.ndarray.item` (without arguments)
* :meth:`~numpy.ndarray.itemset` (only the 1-argument form)
* :meth:`~numpy.ndarray.ptp` (without arguments)
* :meth:`~numpy.ndarray.ravel` (no order argument; 'C' order only)
* :meth:`~numpy.ndarray.repeat` (no axis argument)
* :meth:`~numpy.ndarray.reshape` (only the 1-argument form)
* :meth:`~numpy.ndarray.sort` (without arguments)
* :meth:`~numpy.ndarray.sum` (with or without the ``axis`` and/or ``dtype``
  arguments.)

  * ``axis`` only supports ``integer`` values.
  * If the ``axis`` argument is a compile-time constant, all valid values
    are supported.
    An out-of-range value will result in a ``LoweringError`` at compile-time.
  * If the ``axis`` argument is not a compile-time constant, only values
    from 0 to 3 are supported.
    An out-of-range value will result in a runtime exception.
  * All numeric ``dtypes`` are supported in the ``dtype`` parameter.
    ``timedelta`` arrays can be used as input arrays but ``timedelta`` is not
    supported as ``dtype`` parameter.
  * When a ``dtype`` is given, it determines the type of the internal
    accumulator. When it is not, the selection is made automatically based on
    the input array's ``dtype``, mostly following the same rules as NumPy.
    However, on 64-bit Windows, Numba uses a 64-bit accumulator for integer
    inputs (``int64`` for ``int32`` inputs and ``uint64`` for ``uint32``
    inputs), while NumPy would use a 32-bit accumulator in those cases.


* :meth:`~numpy.ndarray.transpose`
* :meth:`~numpy.ndarray.view` (only the 1-argument form)
* :meth:`~numpy.ndarray.__contains__` 

Where applicable, the corresponding top-level NumPy functions (such as
:func:`numpy.argmax`) are similarly supported.

.. warning::
   Sorting may be slightly slower than Numpy's implementation.


Functions
=========

Linear algebra
--------------

Basic linear algebra is supported on 1-D and 2-D contiguous arrays of
floating-point and complex numbers:

* :func:`numpy.dot`
* :func:`numpy.kron` ('C' and 'F' order only)
* :func:`numpy.outer`
* :func:`numpy.trace` (only the first argument).
* :func:`numpy.vdot`
* On Python 3.5 and above, the matrix multiplication operator from
  :pep:`465` (i.e. ``a @ b`` where ``a`` and ``b`` are 1-D or 2-D arrays).
* :func:`numpy.linalg.cholesky`
* :func:`numpy.linalg.cond` (only non string values in ``p``).
* :func:`numpy.linalg.det`
* :func:`numpy.linalg.eig` (only running with data that does not cause a domain
  change is supported e.g. real input -> real
  output, complex input -> complex output).
* :func:`numpy.linalg.eigh` (only the first argument).
* :func:`numpy.linalg.eigvals` (only running with data that does not cause a
  domain change is supported e.g. real input -> real output,
  complex input -> complex output).
* :func:`numpy.linalg.eigvalsh` (only the first argument).
* :func:`numpy.linalg.inv`
* :func:`numpy.linalg.lstsq`
* :func:`numpy.linalg.matrix_power`
* :func:`numpy.linalg.matrix_rank`
* :func:`numpy.linalg.norm` (only the 2 first arguments and only non string
  values in ``ord``).
* :func:`numpy.linalg.pinv`
* :func:`numpy.linalg.qr` (only the first argument).
* :func:`numpy.linalg.slogdet`
* :func:`numpy.linalg.solve`
* :func:`numpy.linalg.svd` (only the 2 first arguments).

.. note::
   The implementation of these functions needs SciPy to be installed.

Reductions
----------

The following reduction functions are supported:

* :func:`numpy.diff` (only the 2 first arguments)
* :func:`numpy.median` (only the first argument)
* :func:`numpy.nancumprod` (only the first argument)
* :func:`numpy.nancumsum` (only the first argument)
* :func:`numpy.nanmax` (only the first argument)
* :func:`numpy.nanmean` (only the first argument)
* :func:`numpy.nanmedian` (only the first argument)
* :func:`numpy.nanmin` (only the first argument)
* :func:`numpy.nanpercentile` (only the 2 first arguments, complex dtypes
  unsupported)
* :func:`numpy.nanquantile` (only the 2 first arguments, complex dtypes
  unsupported)
* :func:`numpy.nanprod` (only the first argument)
* :func:`numpy.nanstd` (only the first argument)
* :func:`numpy.nansum` (only the first argument)
* :func:`numpy.nanvar` (only the first argument)
* :func:`numpy.percentile` (only the 2 first arguments, complex dtypes
  unsupported)
* :func:`numpy.quantile` (only the 2 first arguments, complex dtypes
  unsupported)

Other functions
---------------

The following top-level functions are supported:

* :func:`numpy.append`
* :func:`numpy.arange`
* :func:`numpy.argsort` (``kind`` key word argument supported for values
  ``'quicksort'`` and ``'mergesort'``)
* :func:`numpy.argwhere`
* :func:`numpy.array` (only the 2 first arguments)
* :func:`numpy.array_equal`
* :func:`numpy.array_split`
* :func:`numpy.asarray` (only the 2 first arguments)
* :func:`numpy.asarray_chkfinite` (only the 2 first arguments)
* :func:`numpy.asfarray`
* :func:`numpy.asfortranarray` (only the first argument)
* :func:`numpy.atleast_1d`
* :func:`numpy.atleast_2d`
* :func:`numpy.atleast_3d`
* :func:`numpy.bartlett`
* :func:`numpy.bincount`
* :func:`numpy.blackman`
* :func:`numpy.broadcast_to` (only the 2 first arguments)
* :func:`numpy.column_stack`
* :func:`numpy.concatenate`
* :func:`numpy.convolve` (only the 2 first arguments)
* :func:`numpy.copy` (only the first argument)
* :func:`numpy.corrcoef` (only the 3 first arguments, requires SciPy)
* :func:`numpy.correlate` (only the 2 first arguments)
* :func:`numpy.count_nonzero` (axis only supports scalar values)
* :func:`numpy.cov` (only the 5 first arguments)
* :func:`numpy.cross` (only the 2 first arguments; at least one of the input
  arrays should have ``shape[-1] == 3``)

  * If ``shape[-1] == 2`` for both inputs, please replace your
    :func:`numpy.cross` call with :func:`numba.np.extensions.cross2d`.

* :func:`numpy.delete` (only the 2 first arguments)
* :func:`numpy.diag`
* :func:`numpy.digitize`
* :func:`numpy.dstack`
* :func:`numpy.dtype` (only the first argument)
* :func:`numpy.ediff1d`
* :func:`numpy.empty` (only the 2 first arguments)
* :func:`numpy.empty_like` (only the 2 first arguments)
* :func:`numpy.expand_dims`
* :func:`numpy.extract`
* :func:`numpy.eye`
* :func:`numpy.fill_diagonal`
* :func:`numpy.flatten` (no order argument; 'C' order only)
* :func:`numpy.flatnonzero`
* :func:`numpy.flip` (no axis argument)
* :func:`numpy.fliplr`
* :func:`numpy.flipud`
* :func:`numpy.frombuffer` (only the 2 first arguments)
* :func:`numpy.full` (only the 3 first arguments)
* :func:`numpy.full_like` (only the 3 first arguments)
* :func:`numpy.hamming`
* :func:`numpy.hanning`
* :func:`numpy.histogram` (only the 3 first arguments)
* :func:`numpy.hstack`
* :func:`numpy.identity`
* :func:`numpy.kaiser`
* :func:`numpy.iscomplex`
* :func:`numpy.iscomplexobj`
* :func:`numpy.isneginf`
* :func:`numpy.isposinf`
* :func:`numpy.isreal`
* :func:`numpy.isrealobj`
* :func:`numpy.isscalar`
* :func:`numpy.interp` (only the 3 first arguments)
* :func:`numpy.intersect1d` (only first 2 arguments, ar1 and ar2)
* :func:`numpy.linspace` (only the 3-argument form)
* :func:`numpy.logspace` (only the 3 first arguments)
* :class:`numpy.ndenumerate`
* :class:`numpy.ndindex`
* :class:`numpy.nditer` (only the first argument)
* :func:`numpy.ones` (only the 2 first arguments)
* :func:`numpy.ones_like` (only the 2 first arguments)
* :func:`numpy.partition` (only the 2 first arguments)
* :func:`numpy.ptp` (only the first argument)
* :func:`numpy.ravel` (no order argument; 'C' order only)
* :func:`numpy.repeat` (no axis argument)
* :func:`numpy.reshape` (no order argument; 'C' order only)
* :func:`numpy.roll` (only the 2 first arguments; second argument ``shift``
  must be an integer)
* :func:`numpy.roots`
* :func:`numpy.rot90` (only the 2 first arguments)
* :func:`numpy.round_`
* :func:`numpy.searchsorted` (only the 3 first arguments)
* :func:`numpy.select` (only using homogeneous lists or tuples for the first
  two arguments, condlist and choicelist). Additionally, these two arguments
  can only contain arrays (unlike Numpy that also accepts tuples).
* :func:`numpy.shape`
* :func:`numpy.sinc`
* :func:`numpy.sort` (no optional arguments)
* :func:`numpy.split`
* :func:`numpy.stack`
* :func:`numpy.swapaxes`
* :func:`numpy.take` (only the 2 first arguments)
* :func:`numpy.take_along_axis` (the axis argument must be a literal value)
* :func:`numpy.transpose`
* :func:`numpy.trapz` (only the 3 first arguments)
* :func:`numpy.tri` (only the 3 first arguments; third argument ``k`` must be an integer)
* :func:`numpy.tril` (second argument ``k`` must be an integer)
* :func:`numpy.tril_indices` (all arguments must be integer)
* :func:`numpy.tril_indices_from` (second argument ``k`` must be an integer)
* :func:`numpy.triu` (second argument ``k`` must be an integer)
* :func:`numpy.triu_indices` (all arguments must be integer)
* :func:`numpy.triu_indices_from` (second argument ``k`` must be an integer)
* :func:`numpy.unique` (only the first argument)
* :func:`numpy.vander`
* :func:`numpy.vstack`
* :func:`numpy.where`
* :func:`numpy.zeros` (only the 2 first arguments)
* :func:`numpy.zeros_like` (only the 2 first arguments)

The following constructors are supported, both with a numeric input (to
construct a scalar) or a sequence (to construct an array):

* :class:`numpy.bool_`
* :class:`numpy.complex64`
* :class:`numpy.complex128`
* :class:`numpy.float32`
* :class:`numpy.float64`
* :class:`numpy.int8`
* :class:`numpy.int16`
* :class:`numpy.int32`
* :class:`numpy.int64`
* :class:`numpy.intc`
* :class:`numpy.intp`
* :class:`numpy.uint8`
* :class:`numpy.uint16`
* :class:`numpy.uint32`
* :class:`numpy.uint64`
* :class:`numpy.uintc`
* :class:`numpy.uintp`

The following machine parameter classes are supported, with all purely numerical
attributes:

* :class:`numpy.iinfo`
* :class:`numpy.finfo` (``machar`` attribute not supported)
* :class:`numpy.MachAr` (with no arguments to the constructor)


Literal arrays
--------------

.. XXX should this part of the user's guide?

Neither Python nor Numba has actual array literals, but you can construct
arbitrary arrays by calling :func:`numpy.array` on a nested tuple::

   a = numpy.array(((a, b, c), (d, e, f)))

(nested lists are not yet supported by Numba)


Modules
=======

.. _numpy-random:

``random``
----------

Numba supports top-level functions from the
`numpy.random <http://docs.scipy.org/doc/numpy/reference/routines.random.html>`_
module, but does not allow you to create individual RandomState instances.
The same algorithms are used as for :ref:`the standard
random module <pysupported-random>` (and therefore the same notes apply),
but with an independent internal state: seeding or drawing numbers from
one generator won't affect the other.

The following functions are supported.

Initialization
''''''''''''''

* :func:`numpy.random.seed`: with an integer argument only

.. warning::
   Calling :func:`numpy.random.seed` from interpreted code (including from :term:`object mode`
   code) will seed the NumPy random generator, not the Numba random generator.
   To seed the Numba random generator, see the example below.

.. code-block:: python

  from numba import njit
  import numpy as np

  @njit
  def seed(a):
      np.random.seed(a)

  @njit
  def rand():
      return np.random.rand()


  # Incorrect seeding
  np.random.seed(1234)
  print(rand())

  np.random.seed(1234)
  print(rand())

  # Correct seeding
  seed(1234)
  print(rand())

  seed(1234)
  print(rand())




Simple random data
''''''''''''''''''

* :func:`numpy.random.rand`
* :func:`numpy.random.randint` (only the first two arguments)
* :func:`numpy.random.randn`
* :func:`numpy.random.random`
* :func:`numpy.random.random_sample`
* :func:`numpy.random.ranf`
* :func:`numpy.random.sample`

Permutations
''''''''''''

* :func:`numpy.random.choice`: the optional *p* argument (probabilities
  array) is not supported
* :func:`numpy.random.permutation`
* :func:`numpy.random.shuffle`: the sequence argument must be a one-dimension
  Numpy array or buffer-providing object (such as a :class:`bytearray`
  or :class:`array.array`)

Distributions
'''''''''''''

The following functions support all arguments. 

* :func:`numpy.random.beta`
* :func:`numpy.random.binomial`
* :func:`numpy.random.chisquare`
* :func:`numpy.random.dirichlet`
* :func:`numpy.random.exponential`
* :func:`numpy.random.f`
* :func:`numpy.random.gamma`
* :func:`numpy.random.geometric`
* :func:`numpy.random.gumbel`
* :func:`numpy.random.hypergeometric`
* :func:`numpy.random.laplace`
* :func:`numpy.random.logistic`
* :func:`numpy.random.lognormal`
* :func:`numpy.random.logseries`
* :func:`numpy.random.multinomial`
* :func:`numpy.random.negative_binomial`
* :func:`numpy.random.normal`
* :func:`numpy.random.pareto`
* :func:`numpy.random.poisson`
* :func:`numpy.random.power`
* :func:`numpy.random.rayleigh`
* :func:`numpy.random.standard_cauchy`
* :func:`numpy.random.standard_exponential`
* :func:`numpy.random.standard_gamma`
* :func:`numpy.random.standard_normal`
* :func:`numpy.random.standard_t`
* :func:`numpy.random.triangular`
* :func:`numpy.random.uniform`
* :func:`numpy.random.vonmises`
* :func:`numpy.random.wald`
* :func:`numpy.random.weibull`
* :func:`numpy.random.zipf`

.. note::
   Calling :func:`numpy.random.seed` from non-Numba code (or from
   :term:`object mode` code) will seed the Numpy random generator, not the
   Numba random generator.

.. note::
   Since version 0.28.0, the generator is thread-safe and fork-safe.  Each
   thread and each process will produce independent streams of random numbers.


``stride_tricks``
-----------------

The following function from the :mod:`numpy.lib.stride_tricks` module
is supported:

* :func:`~numpy.lib.stride_tricks.as_strided` (the *strides* argument
  is mandatory, the *subok* argument is not supported)

.. _supported_ufuncs:

Standard ufuncs
===============

One objective of Numba is having all the
`standard ufuncs in NumPy <http://docs.scipy.org/doc/numpy/reference/ufuncs.html#available-ufuncs>`_
understood by Numba.  When a supported ufunc is found when compiling a
function, Numba maps the ufunc to equivalent native code.  This allows the
use of those ufuncs in Numba code that gets compiled in :term:`nopython mode`.

Limitations
-----------

Right now, only a selection of the standard ufuncs work in :term:`nopython mode`.
Following is a list of the different standard ufuncs that Numba is aware of,
sorted in the same way as in the NumPy documentation.


Math operations
---------------

==============  =============  ===============
    UFUNC                  MODE
--------------  ------------------------------
    name         object mode    nopython mode
==============  =============  ===============
 add                 Yes          Yes
 subtract            Yes          Yes
 multiply            Yes          Yes
 divide              Yes          Yes
 logaddexp           Yes          Yes
 logaddexp2          Yes          Yes
 true_divide         Yes          Yes
 floor_divide        Yes          Yes
 negative            Yes          Yes
 power               Yes          Yes
 float_power         Yes          Yes
 remainder           Yes          Yes
 mod                 Yes          Yes
 fmod                Yes          Yes
 divmod (*)          Yes          Yes
 abs                 Yes          Yes
 absolute            Yes          Yes
 fabs                Yes          Yes
 rint                Yes          Yes
 sign                Yes          Yes
 conj                Yes          Yes
 exp                 Yes          Yes
 exp2                Yes          Yes
 log                 Yes          Yes
 log2                Yes          Yes
 log10               Yes          Yes
 expm1               Yes          Yes
 log1p               Yes          Yes
 sqrt                Yes          Yes
 square              Yes          Yes
 cbrt                Yes          Yes
 reciprocal          Yes          Yes
 conjugate           Yes          Yes
 gcd                 Yes          Yes
 lcm                 Yes          Yes
==============  =============  ===============

(\*) not supported on timedelta types

Trigonometric functions
-----------------------

==============  =============  ===============
    UFUNC                  MODE
--------------  ------------------------------
    name         object mode    nopython mode
==============  =============  ===============
 sin                 Yes          Yes
 cos                 Yes          Yes
 tan                 Yes          Yes
 arcsin              Yes          Yes
 arccos              Yes          Yes
 arctan              Yes          Yes
 arctan2             Yes          Yes
 hypot               Yes          Yes
 sinh                Yes          Yes
 cosh                Yes          Yes
 tanh                Yes          Yes
 arcsinh             Yes          Yes
 arccosh             Yes          Yes
 arctanh             Yes          Yes
 deg2rad             Yes          Yes
 rad2deg             Yes          Yes
 degrees             Yes          Yes
 radians             Yes          Yes
==============  =============  ===============


Bit-twiddling functions
-----------------------

==============  =============  ===============
    UFUNC                  MODE
--------------  ------------------------------
    name         object mode    nopython mode
==============  =============  ===============
 bitwise_and         Yes          Yes
 bitwise_or          Yes          Yes
 bitwise_xor         Yes          Yes
 bitwise_not         Yes          Yes
 invert              Yes          Yes
 left_shift          Yes          Yes
 right_shift         Yes          Yes
==============  =============  ===============


Comparison functions
--------------------

==============  =============  ===============
    UFUNC                  MODE
--------------  ------------------------------
    name         object mode    nopython mode
==============  =============  ===============
 greater             Yes          Yes
 greater_equal       Yes          Yes
 less                Yes          Yes
 less_equal          Yes          Yes
 not_equal           Yes          Yes
 equal               Yes          Yes
 logical_and         Yes          Yes
 logical_or          Yes          Yes
 logical_xor         Yes          Yes
 logical_not         Yes          Yes
 maximum             Yes          Yes
 minimum             Yes          Yes
 fmax                Yes          Yes
 fmin                Yes          Yes
==============  =============  ===============


Floating functions
------------------

==============  =============  ===============
    UFUNC                  MODE
--------------  ------------------------------
    name         object mode    nopython mode
==============  =============  ===============
 isfinite            Yes          Yes
 isinf               Yes          Yes
 isnan               Yes          Yes
 signbit             Yes          Yes
 copysign            Yes          Yes
 nextafter           Yes          Yes
 modf                Yes          No
 ldexp               Yes (*)      Yes
 frexp               Yes          No
 floor               Yes          Yes
 ceil                Yes          Yes
 trunc               Yes          Yes
 spacing             Yes          Yes
==============  =============  ===============

(\*) not supported on windows 32 bit


Datetime functions
------------------

==============  =============  ===============
    UFUNC                  MODE
--------------  ------------------------------
    name         object mode    nopython mode
==============  =============  ===============
 isnat            Yes          Yes
==============  =============  ===============
.. _pysupported:

=========================
Supported Python features
=========================

Apart from the :ref:`pysupported-language` part below, which applies to both
:term:`object mode` and :term:`nopython mode`, this page only lists the
features supported in :term:`nopython mode`.

.. warning::
    Numba behavior differs from Python semantics in some situations.  We
    strongly advise reviewing :ref:`pysemantics` to become familiar with these
    differences.


.. _pysupported-language:

Language
========

Constructs
----------

Numba strives to support as much of the Python language as possible, but
some language features are not available inside Numba-compiled functions.
Below is a quick reference for the support level of Python constructs.


**Supported** constructs:

- conditional branch: ``if .. elif .. else``
- loops: ``while``, ``for .. in``, ``break``, ``continue``
- basic generator: ``yield``
- assertion: ``assert``

**Partially supported** constructs:

- exceptions: ``try .. except``, ``raise``, ``else`` and ``finally``
  (See details in this :ref:`section <pysupported-exception-handling>`)

- context manager:
  ``with`` (only support :ref:`numba.objmode() <with_objmode>`)

- list comprehension (see details in this
  :ref:`section <pysupported-comprehension>`)

**Unsupported** constructs:

- async features: ``async with``, ``async for`` and ``async def``
- class definition: ``class`` (except for :ref:`@jitclass <jitclass>`)
- set, dict and generator comprehensions
- generator delegation: ``yield from``

Functions
---------

Function calls
''''''''''''''

Numba supports function calls using positional and named arguments, as well
as arguments with default values and ``*args`` (note the argument for
``*args`` can only be a tuple, not a list).  Explicit ``**kwargs`` are
not supported.

Function calls to locally defined inner functions are supported as long as
they can be fully inlined.

Functions as arguments
''''''''''''''''''''''

Functions can be passed as argument into another function.  But, they cannot
be returned. For example:

.. code-block:: python

  from numba import jit

  @jit
  def add1(x):
      return x + 1

  @jit
  def bar(fn, x):
      return fn(x)

  @jit
  def foo(x):
      return bar(add1, x)

  # Passing add1 within numba compiled code.
  print(foo(1))
  # Passing add1 into bar from interpreted code
  print(bar(add1, 1))

.. note:: Numba does not handle function objects as real objects.  Once a
          function is assigned to a variable, the variable cannot be
          re-assigned to a different function.


Inner function and closure
'''''''''''''''''''''''''''

Numba now supports inner functions as long as they are non-recursive
and only called locally, but not passed as argument or returned as
result. The use of closure variables (variables defined in outer scopes)
within an inner function is also supported.

Recursive calls
'''''''''''''''

Most recursive call patterns are supported.  The only restriction is that the
recursive callee must have a control-flow path that returns without recursing.
Numba is able to type-infer recursive functions without specifying the function
type signature (which is required in numba 0.28 and earlier).
Recursive calls can even call into a different overload of the function.

.. XXX add reference to NBEP

Generators
----------

Numba supports generator functions and is able to compile them in
:term:`object mode` and :term:`nopython mode`.  The returned generator
can be used both from Numba-compiled code and from regular Python code.

Coroutine features of generators are not supported (i.e. the
:meth:`generator.send`, :meth:`generator.throw`, :meth:`generator.close`
methods).

.. _pysupported-exception-handling:

Exception handling
------------------

``raise`` statement
'''''''''''''''''''

The ``raise`` statement is only supported in the following forms:

* ``raise SomeException``
* ``raise SomeException(<arguments>)``: in :term:`nopython mode`, constructor
  arguments must be :term:`compile-time constants <compile-time constant>`

It is currently unsupported to re-raise an exception created in compiled code.

``try .. except``
'''''''''''''''''

The ``try .. except`` construct is partially supported. The following forms
of are supported:

* the *bare* except that captures all exceptions:

  .. code-block:: python

    try:
        ...
    except:
        ...

* using exactly the ``Exception`` class in the ``except`` clause:

  .. code-block:: python

    try:
      ...
    except Exception:
      ...

  This will match any exception that is a subclass of ``Exception`` as
  expected. Currently, instances of ``Exception`` and it's subclasses are the
  only kind of exception that can be raised in compiled code.

.. warning:: Numba currently masks signals like ``KeyboardInterrupt`` and
  ``SystemExit``. These signaling exceptions are ignored during the execution of
  Numba compiled code. The Python interpreter will handle them as soon as
  the control is returned to it.

Currently, exception objects are not materialized inside compiled functions.
As a result, it is not possible to store an exception object into a user
variable or to re-raise an exception. With this limitation, the only realistic
use-case would look like:

.. code-block:: python

  try:
     do_work()
  except Exception:
     handle_error_case()
     return error_code

``try .. except .. else .. finally``
''''''''''''''''''''''''''''''''''''

The ``else`` block and the ``finally`` block of a ``try .. except`` are
supported:

  .. code-block:: python

    >>> @jit(nopython=True)
    ... def foo():
    ...     try:
    ...         print('main block')
    ...     except Exception:
    ...         print('handler block')
    ...     else:
    ...         print('else block')
    ...     finally:
    ...         print('final block')
    ...
    >>> foo()
    main block
    else block
    final block

The ``try .. finally`` construct without the ``except`` clause is also
supported.

.. _pysupported-builtin-types:

Built-in types
==============

int, bool
---------

Arithmetic operations as well as truth values are supported.

The following attributes and methods are supported:

* ``.conjugate()``
* ``.real``
* ``.imag``

float, complex
--------------

Arithmetic operations as well as truth values are supported.

The following attributes and methods are supported:

* ``.conjugate()``
* ``.real``
* ``.imag``

str
---

Numba supports (Unicode) strings in Python 3.  Strings can be passed into
:term:`nopython mode` as arguments, as well as constructed and returned from
:term:`nopython mode`. As in Python, slices (even of length 1) return a new,
reference counted string.  Optimized code paths for efficiently accessing
single characters may be introduced in the future.

The in-memory representation is the same as was introduced in Python 3.4, with
each string having a tag to indicate whether the string is using a 1, 2, or 4
byte character width in memory.  When strings of different encodings are
combined (as in concatenation), the resulting string automatically uses the
larger character width of the two input strings.  String slices also use the
same character width as the original string, even if the slice could be
represented with a narrower character width.  (These details are invisible to
the user, of course.)

The following constructors, functions, attributes and methods are currently
supported:

* ``str(int)``
* ``len()``
* ``+`` (concatenation of strings)
* ``*`` (repetition of strings)
* ``in``, ``.contains()``
* ``==``, ``<``, ``<=``, ``>``, ``>=`` (comparison)
* ``.capitalize()``
* ``.casefold()``
* ``.center()``
* ``.count()``
* ``.endswith()``
* ``.endswith()``
* ``.expandtabs()``
* ``.find()``
* ``.index()``
* ``.isalnum()``
* ``.isalpha()``
* ``.isdecimal()``
* ``.isdigit()``
* ``.isidentifier()``
* ``.islower()``
* ``.isnumeric()``
* ``.isprintable()``
* ``.isspace()``
* ``.istitle()``
* ``.isupper()``
* ``.join()``
* ``.ljust()``
* ``.lower()``
* ``.lstrip()``
* ``.partition()``
* ``.replace()``
* ``.rfind()``
* ``.rindex()``
* ``.rjust()``
* ``.rpartition()``
* ``.rsplit()``
* ``.rstrip()``
* ``.split()``
* ``.splitlines()``
* ``.startswith()``
* ``.strip()``
* ``.swapcase()``
* ``.title()``
* ``.upper()``
* ``.zfill()``

Regular string literals (e.g. ``"ABC"``) as well as f-strings without format specs
(e.g. ``"ABC_{a+1}"``)
that only use string and integer variables (types with ``str()`` overload)
are supported in :term:`nopython mode`.

Additional operations as well as support for Python 2 strings / Python 3 bytes
will be added in a future version of Numba.  Python 2 Unicode objects will
likely never be supported.

.. warning::
    The performance of some operations is known to be slower than the CPython
    implementation. These include substring search (``in``, ``.contains()``
    and ``find()``) and string creation (like ``.split()``).  Improving the
    string performance is an ongoing task, but the speed of CPython is
    unlikely to be surpassed for basic string operation in isolation.
    Numba is most successfully used for larger algorithms that happen to
    involve strings, where basic string operations are not the bottleneck.


tuple
-----

Tuple support is categorised into two categories based on the contents of a
tuple. The first category is homogeneous tuples, these are tuples where the type
of all the values in the tuple are the same, the second is heterogeneous tuples,
these are tuples where the types of the values are different.

.. note::

    The ``tuple()`` constructor itself is NOT supported.

homogeneous tuples
------------------

An example of a homogeneous tuple:

.. code-block:: python

    homogeneous_tuple = (1, 2, 3, 4)

The following operations are supported on homogeneous tuples:

* Tuple construction.
* Tuple unpacking.
* Comparison between tuples.
* Iteration and indexing.
* Addition (concatenation) between tuples.
* Slicing tuples with a constant slice.
* The index method on tuples.

heterogeneous tuples
--------------------

An example of a heterogeneous tuple:

.. code-block:: python

    heterogeneous_tuple = (1, 2j, 3.0, "a")

The following operations are supported on heterogeneous tuples:

* Comparison between tuples.
* Indexing using an index value that is a compile time constant
  e.g. ``mytuple[7]``, where ``7`` is evidently a constant.
* Iteration over a tuple (requires experimental :func:`literal_unroll` feature,
  see below).

.. warning::
   The following feature (:func:`literal_unroll`) is experimental and was added
   in version 0.47.

To permit iteration over a heterogeneous tuple the special function
:func:`numba.literal_unroll` must be used. This function has no effect other
than to act as a token to permit the use of this feature. Example use:

.. code-block:: python

    from numba import njit, literal_unroll

    @njit
    def foo():
        heterogeneous_tuple = (1, 2j, 3.0, "a")
        for i in literal_unroll(heterogeneous_tuple):
            print(i)

.. warning::
    The following restrictions apply to the use of :func:`literal_unroll`:

    * :func:`literal_unroll` can only be used on tuples and constant lists of
      compile time constants, e.g. ``[1, 2j, 3, "a"]`` and the list not being
      mutated.
    * The only supported use pattern for :func:`literal_unroll` is loop
      iteration.
    * Only one :func:`literal_unroll` call is permitted per loop nest (i.e.
      nested heterogeneous tuple iteration loops are forbidden).
    * The usual type inference/stability rules still apply.

A more involved use of :func:`literal_unroll` might be type specific dispatch,
recall that string and integer literal values are considered their own type,
for example:

.. code-block:: python

    from numba import njit, types, literal_unroll
    from numba.extending import overload

    def dt(x):
        # dummy function to overload
        pass

    @overload(dt, inline='always')
    def ol_dt(li):
        if isinstance(li, types.StringLiteral):
            value = li.literal_value
            if value == "apple":
                def impl(li):
                    return 1
            elif value == "orange":
                def impl(li):
                    return 2
            elif value == "banana":
                def impl(li):
                    return 3
            return impl
        elif isinstance(li, types.IntegerLiteral):
            value = li.literal_value
            if value == 0xca11ab1e:
                def impl(li):
                    # capture the dispatcher literal value
                    return 0x5ca1ab1e + value
                return impl

    @njit
    def foo():
        acc = 0
        for t in literal_unroll(('apple', 'orange', 'banana', 3390155550)):
            acc += dt(t)
        return acc

    print(foo())


list
----


.. warning::
    As of version 0.45.x the internal implementation for the list datatype in
    Numba is changing. Until recently, only a single implementation of the list
    datatype was available, the so-called *reflected-list* (see below).
    However, it was scheduled for deprecation from version 0.44.0 onwards due
    to its limitations. As of version 0.45.0 a new implementation, the
    so-called *typed-list* (see below), is available as an experimental
    feature. For more information, please see: :ref:`deprecation`.

Creating and returning lists from JIT-compiled functions is supported,
as well as all methods and operations.  Lists must be strictly homogeneous:
Numba will reject any list containing objects of different types, even if
the types are compatible (for example, ``[1, 2.5]`` is rejected as it
contains a :class:`int` and a :class:`float`).

For example, to create a list of arrays::

  In [1]: from numba import njit

  In [2]: import numpy as np

  In [3]: @njit
    ...: def foo(x):
    ...:     lst = []
    ...:     for i in range(x):
    ...:         lst.append(np.arange(i))
    ...:     return lst
    ...:

  In [4]: foo(4)
  Out[4]: [array([], dtype=int64), array([0]), array([0, 1]), array([0, 1, 2])]


.. _feature-reflected-list:

List Reflection
'''''''''''''''

In nopython mode, Numba does not operate on Python objects.  ``list`` are
compiled into an internal representation.  Any ``list`` arguments must be
converted into this representation on the way in to nopython mode and their
contained elements must be restored in the original Python objects via a
process called :term:`reflection`.  Reflection is required to maintain the same
semantics as found in regular Python code.  However, the reflection process
can be expensive for large lists and it is not supported for lists that contain
reflected data types.  Users cannot use list-of-list as an argument because
of this limitation.

.. note::
   When passing a list into a JIT-compiled function, any modifications
   made to the list will not be visible to the Python interpreter until
   the function returns.  (A limitation of the reflection process.)

.. warning::
   List sorting currently uses a quicksort algorithm, which has different
   performance characterics than the algorithm used by Python.

.. _feature-list-initial-value:

Initial Values
''''''''''''''
.. warning::
  This is an experimental feature!

Lists that:

* Are constructed using the square braces syntax
* Have values of a literal type

will have their initial value stored in the ``.initial_value`` property on the
type so as to permit inspection of these values at compile time. If required,
to force value based dispatch the :ref:`literally <developer-literally>`
function will accept such a list.

Example:

.. literalinclude:: ../../../numba/tests/doc_examples/test_literal_container_usage.py
   :language: python
   :caption: from ``test_ex_initial_value_list_compile_time_consts`` of ``numba/tests/doc_examples/test_literal_container_usage.py``
   :start-after: magictoken.test_ex_initial_value_list_compile_time_consts.begin
   :end-before: magictoken.test_ex_initial_value_list_compile_time_consts.end
   :dedent: 12
   :linenos:

.. _feature-typed-list:

Typed List
''''''''''

.. note::
  ``numba.typed.List`` is an experimental feature, if you encounter any bugs in
  functionality or suffer from unexpectedly bad performance, please report
  this, ideally by opening an issue on the Numba issue tracker.

As of version 0.45.0 a new implementation of the list data type is available,
the so-called *typed-list*. This is compiled library backed, type-homogeneous
list data type that is an improvement over the *reflected-list* mentioned
above.  Additionally, lists can now be arbitrarily nested. Since the
implementation is considered experimental, you will need to import it
explicitly from the `numba.typed` module::

    In [1]: from numba.typed import List

    In [2]: from numba import njit

    In [3]: @njit
    ...: def foo(l):
    ...:     l.append(23)
    ...:     return l
    ...:

    In [4]: mylist = List()

    In [5]: mylist.append(1)

    In [6]: foo(mylist)
    Out[6]: ListType[int64]([1, 23])


.. note::
    As the typed-list stabilizes it will fully replace the reflected-list and the
    constructors `[]` and `list()` will create a typed-list instead of a
    reflected one.


Here's an example using ``List()`` to create ``numba.typed.List`` inside a
jit-compiled function and letting the compiler infer the item type:

.. literalinclude:: ../../../numba/tests/doc_examples/test_typed_list_usage.py
   :language: python
   :caption: from ``ex_inferred_list_jit`` of ``numba/tests/doc_examples/test_typed_list_usage.py``
   :start-after: magictoken.ex_inferred_list_jit.begin
   :end-before: magictoken.ex_inferred_list_jit.end
   :dedent: 12
   :linenos:

Here's an example of using ``List()`` to create a ``numba.typed.List`` outside of
a jit-compiled function and then using it as an argument to a jit-compiled
function:

.. literalinclude:: ../../../numba/tests/doc_examples/test_typed_list_usage.py
   :language: python
   :caption: from ``ex_inferred_list`` of ``numba/tests/doc_examples/test_typed_list_usage.py``
   :start-after: magictoken.ex_inferred_list.begin
   :end-before: magictoken.ex_inferred_list.end
   :dedent: 12
   :linenos:

Finally, here's an example of using a nested `List()`:

.. literalinclude:: ../../../numba/tests/doc_examples/test_typed_list_usage.py
   :language: python
   :caption: from ``ex_nested_list`` of ``numba/tests/doc_examples/test_typed_list_usage.py``
   :start-after: magictoken.ex_nested_list.begin
   :end-before: magictoken.ex_nested_list.end
   :dedent: 12
   :linenos:

.. _feature-literal-list:

Literal List
''''''''''''

.. warning::
  This is an experimental feature!

Numba supports the use of literal lists containing any values, for example::

  l = ['a', 1, 2j, np.zeros(5,)]

the predominant use of these lists is for use as a configuration object.
The lists appear as a ``LiteralList`` type which inherits from ``Literal``, as a
result the literal values of the list items are available at compile time.
For example:

.. literalinclude:: ../../../numba/tests/doc_examples/test_literal_container_usage.py
   :language: python
   :caption: from ``test_ex_literal_list`` of ``numba/tests/doc_examples/test_literal_container_usage.py``
   :start-after: magictoken.test_ex_literal_list.begin
   :end-before: magictoken.test_ex_literal_list.end
   :dedent: 12
   :linenos:

Important things to note about these kinds of lists:

#. They are immutable, use of mutating methods e.g. ``.pop()`` will result in
   compilation failure. Read-only static access and read only methods are
   supported e.g. ``len()``.
#. Dynamic access of items is not possible, e.g. ``some_list[x]``, for a
   value ``x`` which is not a compile time constant. This is because it's
   impossible to statically determine the type of the item being accessed.
#. Inside the compiler, these lists are actually just tuples with some extra
   things added to make them look like they are lists.
#. They cannot be returned to the interpreter from a compiled function.

.. _pysupported-comprehension:

List comprehension
''''''''''''''''''

Numba supports list comprehension.  For example::


  In [1]: from numba import njit

  In [2]: @njit
    ...: def foo(x):
    ...:     return [[i for i in range(n)] for n in range(x)]
    ...:

  In [3]: foo(3)
  Out[3]: [[], [0], [0, 1]]


.. note::
  Prior to version 0.39.0, Numba did not support the creation of nested lists.


Numba also supports "array comprehension" that is a list comprehension
followed immediately by a call to :func:`numpy.array`. The following
is an example that produces a 2D Numpy array::

    from numba import jit
    import numpy as np

    @jit(nopython=True)
    def f(n):
      return np.array([ [ x * y for x in range(n) ] for y in range(n) ])

In this case, Numba is able to optimize the program to allocate and
initialize the result array directly without allocating intermediate
list objects.  Therefore, the nesting of list comprehension here is
not a problem since a multi-dimensional array is being created here
instead of a nested list.

Additionally, Numba supports parallel array comprehension when combined
with the :ref:`parallel_jit_option` option on CPUs.

set
---

All methods and operations on sets are supported in JIT-compiled functions.

Sets must be strictly homogeneous: Numba will reject any set containing
objects of different types, even if the types are compatible (for example,
``{1, 2.5}`` is rejected as it contains a :class:`int` and a :class:`float`).
The use of reference counted types, e.g. strings, in sets is unsupported.

.. note::
   When passing a set into a JIT-compiled function, any modifications
   made to the set will not be visible to the Python interpreter until
   the function returns.

.. _feature-typed-dict:

Typed Dict
----------

.. warning::
  ``numba.typed.Dict`` is an experimental feature.  The API may change
  in the future releases.

.. note::
  ``dict()`` was not supported in versions prior to 0.44.  Currently, calling
  ``dict()`` translates to calling ``numba.typed.Dict()``.

Numba only supports the use of ``dict()`` without any arguments.  Such use is
semantically equivalent to ``{}`` and ``numba.typed.Dict()``.  It will create
an instance of ``numba.typed.Dict`` where the key-value types will be later
inferred by usage.

Numba does not fully support the Python ``dict`` because it is an untyped
container that can have any Python types as members. To generate efficient
machine code, Numba needs the keys and the values of the dictionary to have
fixed types, declared in advance. To achieve this, Numba has a typed dictionary,
``numba.typed.Dict``, for which the type-inference mechanism must be able to
infer the key-value types by use, or the user must explicitly declare the
key-value type using the ``Dict.empty()`` constructor method.
This typed dictionary has the same API as the Python ``dict``,  it implements
the ``collections.MutableMapping`` interface and is usable in both interpreted
Python code and JIT-compiled Numba functions.
Because the typed dictionary stores keys and values in Numba's native,
unboxed data layout, passing a Numba dictionary into nopython mode has very low
overhead. However, this means that using a typed dictionary from the Python
interpreter is slower than a regular dictionary because Numba has to box and
unbox key and value objects when getting or setting items.

An important difference of the typed dictionary in comparison to Python's
``dict`` is that **implicit casting** occurs when a key or value is stored.
As a result the *setitem* operation may fail should the type-casting fail.

It should be noted that the Numba typed dictionary is implemented using the same
algorithm as the CPython 3.7 dictionary. As a consequence, the typed dictionary
is ordered and has the same collision resolution as the CPython implementation.

Further to the above in relation to type specification, there are limitations
placed on the types that can be used as keys and/or values in the typed
dictionary, most notably the Numba ``Set`` and ``List`` types are currently
unsupported. Acceptable key/value types include but are not limited to: unicode
strings, arrays (value only), scalars, tuples. It is expected that these
limitations will be relaxed as Numba continues to improve.

Here's an example of using ``dict()`` and ``{}`` to create ``numba.typed.Dict``
instances and letting the compiler infer the key-value types:

.. literalinclude:: ../../../numba/tests/doc_examples/test_typed_dict_usage.py
   :language: python
   :caption: from ``test_ex_inferred_dict_njit`` of ``numba/tests/doc_examples/test_typed_dict_usage.py``
   :start-after: magictoken.ex_inferred_dict_njit.begin
   :end-before: magictoken.ex_inferred_dict_njit.end
   :dedent: 12
   :linenos:

Here's an example of creating a ``numba.typed.Dict`` instance from interpreted
code and using the dictionary in jit code:

.. literalinclude:: ../../../numba/tests/doc_examples/test_typed_dict_usage.py
   :language: python
   :caption: from ``test_ex_typed_dict_from_cpython`` of ``numba/tests/doc_examples/test_typed_dict_usage.py``
   :start-after: magictoken.ex_typed_dict_from_cpython.begin
   :end-before: magictoken.ex_typed_dict_from_cpython.end
   :dedent: 12
   :linenos:

Here's an example of creating a ``numba.typed.Dict`` instance from jit code and
using the dictionary in interpreted code:

.. literalinclude:: ../../../numba/tests/doc_examples/test_typed_dict_usage.py
   :language: python
   :caption: from ``test_ex_typed_dict_njit`` of ``numba/tests/doc_examples/test_typed_dict_usage.py``
   :start-after: magictoken.ex_typed_dict_njit.begin
   :end-before: magictoken.ex_typed_dict_njit.end
   :dedent: 12
   :linenos:

It should be noted that ``numba.typed.Dict`` is not thread-safe.
Specifically, functions which modify a dictionary from multiple
threads will potentially corrupt memory, causing a
range of possible failures. However, the dictionary can be safely read from
multiple threads as long as the contents of the dictionary do not
change during the parallel access.

Dictionary comprehension
''''''''''''''''''''''''

Numba supports dictionary comprehension under the assumption that a
``numba.typed.Dict`` instance can be created from the comprehension.  For
example::

  In [1]: from numba import njit

  In [2]: @njit
     ...: def foo(n):
     ...:     return {i: i**2 for i in range(n)}
     ...:

  In [3]: foo(3)
  Out[3]: DictType[int64,int64]<iv=None>({0: 0, 1: 1, 2: 4})

.. _feature-dict-initial-value:

Initial Values
''''''''''''''
.. warning::
  This is an experimental feature!

Typed dictionaries that:

* Are constructed using the curly braces syntax
* Have literal string keys
* Have values of a literal type

will have their initial value stored in the ``.initial_value`` property on the
type so as to permit inspection of these values at compile time. If required,
to force value based dispatch the :ref:`literally <developer-literally>`
function will accept a typed dictionary.

Example:

.. literalinclude:: ../../../numba/tests/doc_examples/test_literal_container_usage.py
   :language: python
   :caption: from ``test_ex_initial_value_dict_compile_time_consts`` of ``numba/tests/doc_examples/test_literal_container_usage.py``
   :start-after: magictoken.test_ex_initial_value_dict_compile_time_consts.begin
   :end-before: magictoken.test_ex_initial_value_dict_compile_time_consts.end
   :dedent: 12
   :linenos:

.. _feature-literal-str-key-dict:

Heterogeneous Literal String Key Dictionary
-------------------------------------------

.. warning::
  This is an experimental feature!

Numba supports the use of statically declared string key to any value
dictionaries, for example::

  d = {'a': 1, 'b': 'data', 'c': 2j}

the predominant use of these dictionaries is to orchestrate advanced compilation
dispatch or as a container for use as a configuration object. The dictionaries
appear as a ``LiteralStrKeyDict`` type which inherits from ``Literal``, as a
result the literal values of the keys and the types of the items are available
at compile time. For example:

.. literalinclude:: ../../../numba/tests/doc_examples/test_literal_container_usage.py
   :language: python
   :caption: from ``test_ex_literal_dict_compile_time_consts`` of ``numba/tests/doc_examples/test_literal_container_usage.py``
   :start-after: magictoken.test_ex_literal_dict_compile_time_consts.begin
   :end-before: magictoken.test_ex_literal_dict_compile_time_consts.end
   :dedent: 12
   :linenos:

Important things to note about these kinds of dictionaries:

#. They are immutable, use of mutating methods e.g. ``.pop()`` will result in
   compilation failure. Read-only static access and read only methods are
   supported e.g. ``len()``.
#. Dynamic access of items is not possible, e.g. ``some_dictionary[x]``, for a
   value ``x`` which is not a compile time constant. This is because it's
   impossible statically determine the type of the item being accessed.
#. Inside the compiler, these dictionaries are actually just named tuples with
   some extra things added to make them look like they are dictionaries.
#. They cannot be returned to the interpreter from a compiled function.
#. The ``.keys()``, ``.values()`` and ``.items()`` methods all functionally
   operate but return tuples opposed to iterables.

None
----

The None value is supported for identity testing (when using an
:class:`~numba.optional` type).


bytes, bytearray, memoryview
----------------------------

The :class:`bytearray` type and, on Python 3, the :class:`bytes` type
support indexing, iteration and retrieving the len().

The :class:`memoryview` type supports indexing, slicing, iteration,
retrieving the len(), and also the following attributes:

* :attr:`~memoryview.contiguous`
* :attr:`~memoryview.c_contiguous`
* :attr:`~memoryview.f_contiguous`
* :attr:`~memoryview.itemsize`
* :attr:`~memoryview.nbytes`
* :attr:`~memoryview.ndim`
* :attr:`~memoryview.readonly`
* :attr:`~memoryview.shape`
* :attr:`~memoryview.strides`


Built-in functions
==================

The following built-in functions are supported:

.. warning::
  Support for ``isinstance`` is an experimental feature. This feature is
  automatically enabled by simply using ``isinstance`` in JIT compiled code.

* :func:`abs`
* :class:`bool`
* :func:`chr`
* :class:`complex`
* :func:`divmod`
* :func:`enumerate`
* :func:`filter`
* :class:`float`
* :func:`hash` (see :ref:`pysupported-hashing` below)
* :class:`int`: only the one-argument form
* :func:`iter`: only the one-argument form
* :func:`isinstance` (experimental support only)
* :func:`len`
* :func:`min`
* :func:`map`
* :func:`max`
* :func:`next`: only the one-argument form
* :func:`ord`
* :func:`print`: only numbers and strings; no ``file`` or ``sep`` argument
* :class:`range`: The only permitted use of range is as a callable function
  (cannot pass range as an argument to a jitted function or return a range from
  a jitted function).
* :func:`round`
* :func:`sorted`: the ``key`` argument is not supported
* :func:`sum`
* :func:`type`: only the one-argument form, and only on some types
  (e.g. numbers and named tuples)
* :func:`zip`

.. _pysupported-hashing:

Hashing
-------

The :func:`hash` built-in is supported and produces hash values for all
supported hashable types with the following Python version specific behavior:

Under Python 3, hash values computed by Numba will exactly match those computed
in CPython under the condition that the :attr:`sys.hash_info.algorithm` is
``siphash24`` (default).

The ``PYTHONHASHSEED`` environment variable influences the hashing behavior in
precisely the manner described in the CPython documentation.


Standard library modules
========================

``array``
---------

Limited support for the :class:`array.array` type is provided through
the buffer protocol.  Indexing, iteration and taking the len() is supported.
All type codes are supported except for ``"u"``.

``cmath``
---------

The following functions from the :mod:`cmath` module are supported:

* :func:`cmath.acos`
* :func:`cmath.acosh`
* :func:`cmath.asin`
* :func:`cmath.asinh`
* :func:`cmath.atan`
* :func:`cmath.atanh`
* :func:`cmath.cos`
* :func:`cmath.cosh`
* :func:`cmath.exp`
* :func:`cmath.isfinite`
* :func:`cmath.isinf`
* :func:`cmath.isnan`
* :func:`cmath.log`
* :func:`cmath.log10`
* :func:`cmath.phase`
* :func:`cmath.polar`
* :func:`cmath.rect`
* :func:`cmath.sin`
* :func:`cmath.sinh`
* :func:`cmath.sqrt`
* :func:`cmath.tan`
* :func:`cmath.tanh`

``collections``
---------------

Named tuple classes, as returned by :func:`collections.namedtuple`, are
supported in the same way regular tuples are supported.  Attribute access
and named parameters in the constructor are also supported.

Creating a named tuple class inside Numba code is *not* supported; the class
must be created at the global level.

.. _ctypes-support:

``ctypes``
----------

Numba is able to call ctypes-declared functions with the following argument
and return types:

* :class:`ctypes.c_int8`
* :class:`ctypes.c_int16`
* :class:`ctypes.c_int32`
* :class:`ctypes.c_int64`
* :class:`ctypes.c_uint8`
* :class:`ctypes.c_uint16`
* :class:`ctypes.c_uint32`
* :class:`ctypes.c_uint64`
* :class:`ctypes.c_float`
* :class:`ctypes.c_double`
* :class:`ctypes.c_void_p`

``enum``
--------

Both :class:`enum.Enum` and :class:`enum.IntEnum` subclasses are supported.

``math``
--------

The following functions from the :mod:`math` module are supported:

* :func:`math.acos`
* :func:`math.acosh`
* :func:`math.asin`
* :func:`math.asinh`
* :func:`math.atan`
* :func:`math.atan2`
* :func:`math.atanh`
* :func:`math.ceil`
* :func:`math.copysign`
* :func:`math.cos`
* :func:`math.cosh`
* :func:`math.degrees`
* :func:`math.erf`
* :func:`math.erfc`
* :func:`math.exp`
* :func:`math.expm1`
* :func:`math.fabs`
* :func:`math.floor`
* :func:`math.frexp`
* :func:`math.gamma`
* :func:`math.gcd`
* :func:`math.hypot`
* :func:`math.isfinite`
* :func:`math.isinf`
* :func:`math.isnan`
* :func:`math.ldexp`
* :func:`math.lgamma`
* :func:`math.log`
* :func:`math.log10`
* :func:`math.log1p`
* :func:`math.pow`
* :func:`math.radians`
* :func:`math.sin`
* :func:`math.sinh`
* :func:`math.sqrt`
* :func:`math.tan`
* :func:`math.tanh`
* :func:`math.trunc`

``operator``
------------

The following functions from the :mod:`operator` module are supported:

* :func:`operator.add`
* :func:`operator.and_`
* :func:`operator.eq`
* :func:`operator.floordiv`
* :func:`operator.ge`
* :func:`operator.gt`
* :func:`operator.iadd`
* :func:`operator.iand`
* :func:`operator.ifloordiv`
* :func:`operator.ilshift`
* :func:`operator.imatmul` (Python 3.5 and above)
* :func:`operator.imod`
* :func:`operator.imul`
* :func:`operator.invert`
* :func:`operator.ior`
* :func:`operator.ipow`
* :func:`operator.irshift`
* :func:`operator.isub`
* :func:`operator.itruediv`
* :func:`operator.ixor`
* :func:`operator.le`
* :func:`operator.lshift`
* :func:`operator.lt`
* :func:`operator.matmul` (Python 3.5 and above)
* :func:`operator.mod`
* :func:`operator.mul`
* :func:`operator.ne`
* :func:`operator.neg`
* :func:`operator.not_`
* :func:`operator.or_`
* :func:`operator.pos`
* :func:`operator.pow`
* :func:`operator.rshift`
* :func:`operator.sub`
* :func:`operator.truediv`
* :func:`operator.xor`

``functools``
-------------

The :func:`functools.reduce` function is supported but the `initializer`
argument is required.

.. _pysupported-random:

``random``
----------

Numba supports top-level functions from the :mod:`random` module, but does
not allow you to create individual Random instances.  A Mersenne-Twister
generator is used, with a dedicated internal state.  It is initialized at
startup with entropy drawn from the operating system.

* :func:`random.betavariate`
* :func:`random.expovariate`
* :func:`random.gammavariate`
* :func:`random.gauss`
* :func:`random.getrandbits`: number of bits must not be greater than 64
* :func:`random.lognormvariate`
* :func:`random.normalvariate`
* :func:`random.paretovariate`
* :func:`random.randint`
* :func:`random.random`
* :func:`random.randrange`
* :func:`random.seed`: with an integer argument only
* :func:`random.shuffle`: the sequence argument must be a one-dimension
  Numpy array or buffer-providing object (such as a :class:`bytearray`
  or :class:`array.array`); the second (optional) argument is not supported
* :func:`random.uniform`
* :func:`random.triangular`
* :func:`random.vonmisesvariate`
* :func:`random.weibullvariate`

.. warning::
   Calling :func:`random.seed` from non-Numba code (or from :term:`object mode`
   code) will seed the Python random generator, not the Numba random generator.
   To seed the Numba random generator, see the example below.

.. code-block:: python

  from numba import njit
  import random

  @njit
  def seed(a):
      random.seed(a)

  @njit
  def rand():
      return random.random()


  # Incorrect seeding
  random.seed(1234)
  print(rand())

  random.seed(1234)
  print(rand())

  # Correct seeding
  seed(1234)
  print(rand())

  seed(1234)
  print(rand())


.. note::
   Since version 0.28.0, the generator is thread-safe and fork-safe.  Each
   thread and each process will produce independent streams of random numbers.

.. seealso::
   Numba also supports most additional distributions from the :ref:`Numpy
   random module <numpy-random>`.

``heapq``
---------

The following functions from the :mod:`heapq` module are supported:

* :func:`heapq.heapify`
* :func:`heapq.heappop`
* :func:`heapq.heappush`
* :func:`heapq.heappushpop`
* :func:`heapq.heapreplace`
* :func:`heapq.nlargest` : first two arguments only
* :func:`heapq.nsmallest` : first two arguments only

Note: the heap must be seeded with at least one value to allow its type to be
inferred; heap items are assumed to be homogeneous in type.


Third-party modules
===================

.. I put this here as there's only one module (apart from Numpy), otherwise
   it should be a separate page.

.. _cffi-support:

``cffi``
--------

Similarly to ctypes, Numba is able to call into `cffi`_-declared external
functions, using the following C types and any derived pointer types:

* :c:type:`char`
* :c:type:`short`
* :c:type:`int`
* :c:type:`long`
* :c:type:`long long`
* :c:type:`unsigned char`
* :c:type:`unsigned short`
* :c:type:`unsigned int`
* :c:type:`unsigned long`
* :c:type:`unsigned long long`
* :c:type:`int8_t`
* :c:type:`uint8_t`
* :c:type:`int16_t`
* :c:type:`uint16_t`
* :c:type:`int32_t`
* :c:type:`uint32_t`
* :c:type:`int64_t`
* :c:type:`uint64_t`
* :c:type:`float`
* :c:type:`double`
* :c:type:`ssize_t`
* :c:type:`size_t`
* :c:type:`void`

The ``from_buffer()`` method of ``cffi.FFI`` and ``CompiledFFI`` objects is
supported for passing Numpy arrays and other buffer-like objects.  Only
*contiguous* arguments are accepted.  The argument to ``from_buffer()``
is converted to a raw pointer of the appropriate C type (for example a
``double *`` for a ``float64`` array).

Additional type mappings for the conversion from a buffer to the appropriate C
type may be registered with Numba. This may include struct types, though it is
only permitted to call functions that accept pointers to structs - passing a
struct by value is unsupported. For registering a mapping, use:

.. function:: numba.core.typing.cffi_utils.register_type(cffi_type, numba_type)

Out-of-line cffi modules must be registered with Numba prior to the use of any
of their functions from within Numba-compiled functions:

.. function:: numba.core.typing.cffi_utils.register_module(mod)

   Register the cffi out-of-line module ``mod`` with Numba.

Inline cffi modules require no registration.

.. _cffi: https://cffi.readthedocs.org/
.. _numba-types:

====================
Types and signatures
====================

Rationale
=========

As an optimizing compiler, Numba needs to decide on the type of each
variable to generate efficient machine code.  Python's standard types
are not precise enough for that, so we had to develop our own fine-grained
type system.

You will encounter Numba types mainly when trying to inspect the results
of Numba's type inference, for :ref:`debugging <numba-envvars>` or
:ref:`educational <architecture>` purposes.  However, you need to use
types explicitly if compiling code :ref:`ahead-of-time <pycc>`.


Signatures
==========

A signature specifies the type of a function.  Exactly which kind
of signature is allowed depends on the context (:term:`AOT` or :term:`JIT`
compilation), but signatures always involve some representation of Numba
types to specify the concrete types for the function's arguments and,
if required, the function's return type.

An example function signature would be the string ``"f8(i4, i4)"``
(or the equivalent ``"float64(int32, int32)"``) which specifies a
function taking two 32-bit integers and returning a double-precision float.


Basic types
===========

The most basic types can be expressed through simple expressions.  The
symbols below refer to attributes of the main ``numba`` module (so if
you read "boolean", it means that symbol can be accessed as ``numba.boolean``).
Many types are available both as a canonical name and a shorthand alias,
following Numpy's conventions.

Numbers
-------

The following table contains the elementary numeric types currently defined
by Numba and their aliases.

===================     =========        ===================================
Type name(s)            Shorthand        Comments
===================     =========        ===================================
boolean                 b1               represented as a byte
uint8, byte             u1               8-bit unsigned byte
uint16                  u2               16-bit unsigned integer
uint32                  u4               32-bit unsigned integer
uint64                  u8               64-bit unsigned integer

int8, char              i1               8-bit signed byte
int16                   i2               16-bit signed integer
int32                   i4               32-bit signed integer
int64                   i8               64-bit signed integer

intc                    --               C int-sized integer
uintc                   --               C int-sized unsigned integer
intp                    --               pointer-sized integer
uintp                   --               pointer-sized unsigned integer
ssize_t                 --               C ssize_t
size_t                  --               C size_t

float32                 f4               single-precision floating-point number
float64, double         f8               double-precision floating-point number

complex64               c8               single-precision complex number
complex128              c16              double-precision complex number
===================     =========        ===================================

Arrays
------

The easy way to declare :class:`~numba.types.Array` types is to subscript an
elementary type according to the number of dimensions. For example a 
1-dimension single-precision array::

   >>> numba.float32[:]
   array(float32, 1d, A)

or a 3-dimension array of the same underlying type::

   >>> numba.float32[:, :, :]
   array(float32, 3d, A)

This syntax defines array types with no particular layout (producing code
that accepts both non-contiguous and contiguous arrays), but you can
specify a particular contiguity by using the ``::1`` index either at
the beginning or the end of the index specification::

   >>> numba.float32[::1]
   array(float32, 1d, C)
   >>> numba.float32[:, :, ::1]
   array(float32, 3d, C)
   >>> numba.float32[::1, :, :]
   array(float32, 3d, F)

Functions
---------

.. warning::
   The feature of considering functions as first-class type objects is
   under development.

Functions are often considered as certain transformations of
input arguments to output values. Within Numba :term:`JIT` compiled
functions, the functions can also be considered as objects, that is,
functions can be passed around as arguments or return values, or used
as items in sequences, in addition to being callable.

First-class function support is enabled for all Numba :term:`JIT`
compiled functions and Numba ``cfunc`` compiled functions except when:

- using a non-CPU compiler,
- the compiled function is a Python generator,
- the compiled function has Omitted arguments,
- or the compiled function returns Optional value.

To disable first-class function support, use ``no_cfunc_wrapper=True``
decorator option.

For instance, consider an example where the Numba :term:`JIT` compiled
function applies user-specified functions as a composition to an input
argument::

    >>> @numba.njit
    ... def composition(funcs, x):
    ...     r = x
    ...     for f in funcs[::-1]:
    ...         r = f(r)
    ...     return r
    ...
    >>> @numba.cfunc("double(double)")
    ... def a(x):
    ...     return x + 1.0
    ...
    >>> @numba.njit
    ... def b(x):
    ...     return x * x
    ...
    >>> composition((a, b), 0.5), 0.5 ** 2 + 1
    (1.25, 1.25)
    >>> composition((b, a, b, b, a), 0.5), b(a(b(b(a(0.5)))))
    (36.75390625, 36.75390625)

Here, ``cfunc`` compiled functions ``a`` and ``b`` are considered as
first-class function objects because these are passed in to the Numba
:term:`JIT` compiled function ``composition`` as arguments, that is, the
``composition`` is :term:`JIT` compiled independently from its argument function
objects (that are collected in the input argument ``funcs``).

Currently, first-class function objects can be Numba ``cfunc`` compiled
functions, :term:`JIT` compiled functions, and objects that implement the
Wrapper Address Protocol (WAP, see below) with the following restrictions:

========================   ============   ==============   ===========
Context                    JIT compiled   cfunc compiled   WAP objects
========================   ============   ==============   ===========
Can be used as arguments   yes            yes              yes
Can be called              yes            yes              yes
Can be used as items       yes\*          yes              yes
Can be returned            yes            yes              yes
Namespace scoping          yes            yes              yes
Automatic overload         yes            no               no
========================   ============   ==============   ===========

\* at least one of the items in a sequence of first-class function objects must
have a precise type.


Wrapper Address Protocol - WAP
++++++++++++++++++++++++++++++

Wrapper Address Protocol provides an API for making any Python object
a first-class function for Numba :term:`JIT` compiled functions. This assumes
that the Python object represents a compiled function that can be
called via its memory address (function pointer value) from Numba :term:`JIT`
compiled functions. The so-called WAP objects must define the
following two methods:

.. method:: __wrapper_address__(self) -> int

            Return the memory address of a first-class function. This
            method is used when a Numba :term:`JIT` compiled function tries to
            call the given WAP instance.

.. method:: signature(self) -> numba.typing.Signature

            Return the signature of the given first-class
            function. This method is used when passing in the given
            WAP instance to a Numba :term:`JIT` compiled function.

In addition, the WAP object may implement the ``__call__``
method. This is necessary when calling WAP objects from Numba
:term:`JIT` compiled functions in :term:`object mode`.

As an example, let us call the standard math library function ``cos``
within a Numba :term:`JIT` compiled function. The memory address of ``cos`` can
be established after loading the math library and using the ``ctypes``
package::

    >>> import numba, ctypes, ctypes.util, math
    >>> libm = ctypes.cdll.LoadLibrary(ctypes.util.find_library('m'))
    >>> class LibMCos(numba.types.WrapperAddressProtocol):
    ...     def __wrapper_address__(self):
    ...         return ctypes.cast(libm.cos, ctypes.c_voidp).value
    ...     def signature(self):
    ...         return numba.float64(numba.float64)
    ...
    >>> @numba.njit
    ... def foo(f, x):
    ...     return f(x)
    ...
    >>> foo(LibMCos(), 0.0)
    1.0
    >>> foo(LibMCos(), 0.5), math.cos(0.5)
    (0.8775825618903728, 0.8775825618903728)

Miscellaneous Types
-------------------

There are some non-numerical types that do not fit into the other categories.

===================   =================================================
Type name(s)          Comments
===================   =================================================
pyobject              generic Python object
voidptr               raw pointer, no operations can be performed on it
===================   =================================================

Advanced types
==============

For more advanced declarations, you have to explicitly call helper
functions or classes provided by Numba.

.. warning::
   The APIs documented here are not guaranteed to be stable.  Unless
   necessary, it is recommended to let Numba infer argument types by using
   the :ref:`signature-less variant of @jit <jit-lazy>`.

.. A word of note: I only documented those types that can be genuinely
   useful to users, i.e. types that can be passed as parameters to a JIT
   function.  Other types such as tuple are only usable in type inference.


Inference
---------

.. function:: numba.typeof(value)

   Create a Numba type accurately describing the given Python *value*.
   ``ValueError`` is raised if the value isn't supported in
   :term:`nopython mode`.

   ::

      >>> numba.typeof(np.empty(3))
      array(float64, 1d, C)
      >>> numba.typeof((1, 2.0))
      (int64, float64)
      >>> numba.typeof([0])
      reflected list(int64)


Numpy scalars
-------------

Instead of using :func:`~numba.typeof`, non-trivial scalars such as
structured types can also be constructed programmatically.

.. function:: numba.from_dtype(dtype)

   Create a Numba type corresponding to the given Numpy *dtype*::

      >>> struct_dtype = np.dtype([('row', np.float64), ('col', np.float64)])
      >>> ty = numba.from_dtype(struct_dtype)
      >>> ty
      Record([('row', '<f8'), ('col', '<f8')])
      >>> ty[:, :]
      unaligned array(Record([('row', '<f8'), ('col', '<f8')]), 2d, A)

.. class:: numba.types.NPDatetime(unit)

   Create a Numba type for Numpy datetimes of the given *unit*.  *unit*
   should be a string amongst the codes recognized by Numpy (e.g.
   ``Y``, ``M``, ``D``, etc.).

.. class:: numba.types.NPTimedelta(unit)

   Create a Numba type for Numpy timedeltas of the given *unit*.  *unit*
   should be a string amongst the codes recognized by Numpy (e.g.
   ``Y``, ``M``, ``D``, etc.).

   .. seealso::
      Numpy `datetime units <http://docs.scipy.org/doc/numpy/reference/arrays.datetime.html#datetime-units>`_.


Arrays
------

.. class:: numba.types.Array(dtype, ndim, layout)

   Create an array type.  *dtype* should be a Numba type.  *ndim* is the
   number of dimensions of the array (a positive integer).  *layout*
   is a string giving the layout of the array: ``A`` means any layout, ``C``
   means C-contiguous and ``F`` means Fortran-contiguous.


Optional types
--------------

.. class:: numba.optional(typ)

   Create an optional type based on the underlying Numba type *typ*.
   The optional type will allow any value of either *typ* or :const:`None`.

   ::

      >>> @jit((optional(intp),))
      ... def f(x):
      ...     return x is not None
      ...
      >>> f(0)
      True
      >>> f(None)
      False


Type annotations
-----------------

.. function:: numba.extending.as_numba_type(py_type)

   Create a Numba type corresponding to the given Python *type annotation*.
   ``TypingError`` is raised if the type annotation can't be mapped to a Numba
   type.  This function is meant to be used at statically compile time to
   evaluate Python type annotations.  For runtime checking of Python objects
   see ``typeof`` above.

   For any numba type, ``as_numba_type(nb_type) == nb_type``.

      >>> numba.extending.as_numba_type(int)
      int64
      >>> import typing  # the Python library, not the Numba one
      >>> numba.extending.as_numba_type(typing.List[float])
      ListType[float64]
      >>> numba.extending.as_numba_type(numba.int32)
      int32

   ``as_numba_type`` is automatically updated to include any ``@jitclass``.

      >>> @jitclass
      ... class Counter:
      ...     x: int
      ...
      ...     def __init__(self):
      ...         self.x = 0
      ...
      ...     def inc(self):
      ...         old_val = self.x
      ...         self.x += 1
      ...         return old_val
      ...
      >>> numba.extending.as_numba_type(Counter)
      instance.jitclass.Counter#11bad4278<x:int64>

   Currently ``as_numba_type`` is only used to infer fields for ``@jitclass``.
.. _aot-compilation:

Ahead-of-Time compilation
=========================

.. currentmodule:: numba.pycc

.. class:: CC(extension_name, source_module=None)

   An object used to generate compiled extensions from Numba-compiled
   Python functions.  *extension_name* is the name of the extension
   to be generated.  *source_module* is the Python module
   containing the functions; if ``None``, it is inferred by examining
   the call stack.

   :class:`CC` instances have the following attributes and methods:

   .. attribute:: name

      (read-only attribute) The name of the extension module to be generated.

   .. attribute:: output_dir

      (read-write attribute) The directory the extension module will be
      written into.  By default it is the directory the *source_module* is
      located in.

   .. attribute:: output_file

      (read-write attribute) The name of the file the extension module will
      be written to.  By default this follows the Python naming convention
      for the current platform.

   .. attribute:: target_cpu

      (read-write attribute) The name of the CPU model to generate code for.
      This will select the appropriate instruction set extensions.  By
      default, a generic CPU is selected in order to produce portable code.

      Recognized names for this attribute depend on the current architecture
      and LLVM version.  If you have LLVM installed, ``llc -mcpu=help``
      will give you a list.  Examples on x86-64 are ``"ivybridge"``,
      ``"haswell"``, ``"skylake"`` or ``"broadwell"``.  You can also give
      the value ``"host"`` which will select the current host CPU.

   .. attribute:: verbose

      (read-write attribute) If true, print out information while
      compiling the extension.  False by default.

   .. decorator:: export(exported_name, sig)

      Mark the decorated function for compilation with the signature *sig*.
      The compiled function will be exposed as *exported_name* in the
      generated extension module.

      All exported names within a given :class:`CC` instance must be
      distinct, otherwise an exception is raised.

   .. method:: compile()

      Compile all exported functions and generate the extension module
      as specified by :attr:`output_dir` and :attr:`output_file`.

   .. method:: distutils_extension(**kwargs)

      Return a :py:class:`distutils.core.Extension` instance allowing
      to integrate generation of the extension module in a conventional
      ``setup.py``-driven build process.  The optional *kwargs* let you
      pass optional parameters to the :py:class:`~distutils.core.Extension`
      constructor.

      In this mode of operation, it is not necessary to call :meth:`compile`
      yourself.  Also, :attr:`output_dir` and :attr:`output_file` will be
      ignored.

.. _numba-envvars:

Environment variables
=====================

.. note:: This section relates to environment variables that impact Numba's
          runtime, for compile time environment variables see
          :ref:`numba-source-install-env_vars`.

Numba allows its behaviour to be changed through the use of environment
variables. Unless otherwise mentioned, those variables have integer values and
default to zero.

For convenience, Numba also supports the use of a configuration file to persist
configuration settings. Note: To use this feature ``pyyaml`` must be installed.

The configuration file must be named ``.numba_config.yaml`` and be present in
the directory from which the Python interpreter is invoked. The configuration
file, if present, is read for configuration settings before the environment
variables are searched. This means that the environment variable settings will
override the settings obtained from a configuration file (the configuration file
is for setting permanent preferences whereas the environment variables are for
ephemeral preferences).

The format of the configuration file is a dictionary in ``YAML`` format that
maps the environment variables below (without the ``NUMBA_`` prefix) to a
desired value. For example, to permanently switch on developer mode
(``NUMBA_DEVELOPER_MODE`` environment variable) and control flow graph printing
(``NUMBA_DUMP_CFG`` environment variable), create a configuration file with the
contents::

    developer_mode: 1
    dump_cfg: 1

This can be especially useful in the case of wanting to use a set color scheme
based on terminal background color. For example, if the terminal background
color is black, the ``dark_bg`` color scheme would be well suited and can be set
for permanent use by adding::

    color_scheme: dark_bg

Jit flags
---------

These variables globally override flags to the :func:`~numba.jit` decorator.

.. envvar:: NUMBA_BOUNDSCHECK

   If set to 0 or 1, globally disable or enable bounds checking, respectively.
   The default if the variable is not set or set to an empty string is to use
   the ``boundscheck`` flag passed to the :func:`~numba.jit` decorator for a
   given function. See the documentation of :ref:`@jit
   <jit-decorator-boundscheck>` for more information.

   Note, due to limitations in numba, the bounds checking currently produces
   exception messages that do not match those from NumPy. If you set
   ``NUMBA_FULL_TRACEBACKS=1``, the full exception message with the axis,
   index, and shape information will be printed to the terminal.

Debugging
---------

These variables influence what is printed out during compilation of
:term:`JIT functions <JIT function>`.

.. envvar:: NUMBA_DEVELOPER_MODE

    If set to non-zero, developer mode produces full tracebacks and disables
    help instructions. Default is zero.

.. envvar:: NUMBA_FULL_TRACEBACKS

    If set to non-zero, enable full tracebacks when an exception occurs.
    Defaults to the value set by `NUMBA_DEVELOPER_MODE`.

.. envvar:: NUMBA_SHOW_HELP

    If set to non-zero, show resources for getting help. Default is zero.

.. envvar:: NUMBA_CAPTURED_ERRORS

    Alters the way in which Numba captures and handles exceptions that do not
    inherit from ``numba.core.errors.NumbaError`` during compilation (e.g.
    standard Python exceptions). This does not impact runtime exception
    handling. Valid values are:

    - ``"old_style"`` (default): this is the exception handling behaviour that
      is present in Numba versions <= 0.54.x. Numba will capture and wrap all
      errors occuring in compilation and depending on the compilation phase they
      will likely materialize as part of the message in a ``TypingError`` or a
      ``LoweringError``.
    - ``"new_style"`` this will treat any exception that does not inherit from
      ``numba.core.errors.NumbaError`` **and** is raised during compilation as a
      "hard error", i.e. the exception will propagate and compilation will halt.
      The purpose of this new style is to differentiate between intentionally
      raised exceptions and those which occur due to mistakes. For example, if
      an ``AttributeError`` occurs in the typing of an ``@overload`` function,
      under this new behaviour it is assumed that this a mistake in the
      implementation and compilation will halt due to this exception. This
      behaviour will eventually become the default.

.. envvar:: NUMBA_DISABLE_ERROR_MESSAGE_HIGHLIGHTING

    If set to non-zero error message highlighting is disabled. This is useful
    for running the test suite on CI systems.

.. envvar:: NUMBA_COLOR_SCHEME

   Alters the color scheme used in error reporting (requires the ``colorama``
   package to be installed to work). Valid values are:

   - ``no_color`` No color added, just bold font weighting.
   - ``dark_bg`` Suitable for terminals with a dark background.
   - ``light_bg`` Suitable for terminals with a light background.
   - ``blue_bg`` Suitable for terminals with a blue background.
   - ``jupyter_nb`` Suitable for use in Jupyter Notebooks.

   *Default value:* ``no_color``. The type of the value is ``string``.

.. envvar:: NUMBA_HIGHLIGHT_DUMPS

   If set to non-zero and ``pygments`` is installed, syntax highlighting is
   applied to Numba IR, LLVM IR and assembly dumps. Default is zero.

.. envvar:: NUMBA_DISABLE_PERFORMANCE_WARNINGS

   If set to non-zero the issuing of performance warnings is disabled. Default
   is zero.

.. envvar:: NUMBA_DEBUG

   If set to non-zero, print out all possible debugging information during
   function compilation.  Finer-grained control can be obtained using other
   variables below.

.. envvar:: NUMBA_DEBUG_FRONTEND

   If set to non-zero, print out debugging information during operation
   of the compiler frontend, up to and including generation of the Numba
   Intermediate Representation.

.. envvar:: NUMBA_DEBUGINFO

   If set to non-zero, enable debug for the full application by setting
   the default value of the ``debug`` option in ``jit``. Beware that
   enabling debug info significantly increases the memory consumption
   for each compiled function.
   Default value equals to the value of `NUMBA_ENABLE_PROFILING`.

.. envvar:: NUMBA_EXTEND_VARIABLE_LIFETIMES

    If set to non-zero, extend the lifetime of variables to the end of the block
    in which their lifetime ends. This is particularly useful in conjunction
    with :envvar:`NUMBA_DEBUGINFO` as it helps with introspection of values.
    Default is zero.

.. envvar:: NUMBA_GDB_BINARY

   Set the ``gdb`` binary for use in Numba's ``gdb`` support, this takes the
   form  of a path and full name of the binary, for example:
   ``/path/from/root/to/binary/name_of_gdb_binary`` This is to permit
   the use of a ``gdb`` from a non-default location with a non-default name. If
   not set ``gdb`` is assumed to reside at ``/usr/bin/gdb``.

.. envvar:: NUMBA_DEBUG_TYPEINFER

   If set to non-zero, print out debugging information about type inference.

.. envvar:: NUMBA_ENABLE_PROFILING

   Enables JIT events of LLVM in order to support profiling of jitted functions.
   This option is automatically enabled under certain profilers.

.. envvar:: NUMBA_TRACE

   If set to non-zero, trace certain function calls (function entry and exit
   events, including arguments and return values).

.. envvar:: NUMBA_DUMP_BYTECODE

   If set to non-zero, print out the Python :py:term:`bytecode` of
   compiled functions.

.. envvar:: NUMBA_DUMP_CFG

   If set to non-zero, print out information about the Control Flow Graph
   of compiled functions.

.. envvar:: NUMBA_DUMP_IR

   If set to non-zero, print out the Numba Intermediate Representation
   of compiled functions.


.. envvar:: NUMBA_DUMP_SSA

   If set to non-zero, print out the Numba Intermediate Representation of
   compiled functions after conversion to Static Single Assignment (SSA) form.

.. envvar:: NUMBA_DEBUG_PRINT_AFTER

   Dump the Numba IR after declared pass(es). This is useful for debugging IR
   changes made by given passes. Accepted values are:

   * Any pass name (as given by the ``.name()`` method on the class)
   * Multiple pass names as a comma separated list, i.e. ``"foo_pass,bar_pass"``
   * The token ``"all"``, which will print after all passes.

   The default value is ``"none"`` so as to prevent output.

.. envvar:: NUMBA_DUMP_ANNOTATION

   If set to non-zero, print out types annotations for compiled functions.

.. envvar:: NUMBA_DUMP_LLVM

   Dump the unoptimized LLVM assembly source of compiled functions.
   Unoptimized code is usually very verbose; therefore,
   :envvar:`NUMBA_DUMP_OPTIMIZED` is recommended instead.

.. envvar:: NUMBA_DUMP_FUNC_OPT

   Dump the LLVM assembly source after the LLVM "function optimization"
   pass, but before the "module optimization" pass.  This is useful mostly
   when developing Numba itself, otherwise use :envvar:`NUMBA_DUMP_OPTIMIZED`.

.. envvar:: NUMBA_DUMP_OPTIMIZED

   Dump the LLVM assembly source of compiled functions after all
   optimization passes.  The output includes the raw function as well as
   its CPython-compatible wrapper (whose name begins with ``wrapper.``).
   Note that the function is often inlined inside the wrapper, as well.

.. envvar:: NUMBA_DEBUG_ARRAY_OPT

   Dump debugging information related to the processing associated with
   the ``parallel=True`` jit decorator option.

.. envvar:: NUMBA_DEBUG_ARRAY_OPT_RUNTIME

   Dump debugging information related to the runtime scheduler associated
   with the ``parallel=True`` jit decorator option.

.. envvar:: NUMBA_DEBUG_ARRAY_OPT_STATS

   Dump statistics about how many operators/calls are converted to
   parallel for-loops and how many are fused together, which are associated
   with the ``parallel=True`` jit decorator option.

.. envvar:: NUMBA_PARALLEL_DIAGNOSTICS

   If set to an integer value between 1 and 4 (inclusive) diagnostic information
   about parallel transforms undertaken by Numba will be written to STDOUT. The
   higher the value set the more detailed the information produced.

.. envvar:: NUMBA_DUMP_ASSEMBLY

   Dump the native assembly code of compiled functions.

.. envvar:: NUMBA_LLVM_PASS_TIMINGS

    Set to ``1`` to enable recording of pass timings in LLVM;
    e.g. ``NUMBA_LLVM_PASS_TIMINGS=1``.
    See :ref:`developer-llvm-timings`.

    *Default value*: ``0`` (Off)

.. seealso::
   :ref:`numba-troubleshooting` and :ref:`architecture`.


Compilation options
-------------------

.. envvar:: NUMBA_OPT

   The optimization level; this option is passed straight to LLVM.

   *Default value:* 3

.. envvar:: NUMBA_LOOP_VECTORIZE

   If set to non-zero, enable LLVM loop vectorization.

   *Default value:* 1 (except on 32-bit Windows)

.. envvar:: NUMBA_SLP_VECTORIZE

   If set to non-zero, enable LLVM superword-level parallelism vectorization.

   *Default value:* 1

.. envvar:: NUMBA_ENABLE_AVX

   If set to non-zero, enable AVX optimizations in LLVM.  This is disabled
   by default on Sandy Bridge and Ivy Bridge architectures as it can sometimes
   result in slower code on those platforms.

.. envvar:: NUMBA_DISABLE_INTEL_SVML

    If set to non-zero and Intel SVML is available, the use of SVML will be
    disabled.

.. envvar:: NUMBA_DISABLE_JIT

   Disable JIT compilation entirely.  The :func:`~numba.jit` decorator acts
   as if it performs no operation, and the invocation of decorated functions
   calls the original Python function instead of a compiled version.  This
   can be useful if you want to run the Python debugger over your code.

.. envvar:: NUMBA_CPU_NAME
.. envvar:: NUMBA_CPU_FEATURES

    Override CPU and CPU features detection.
    By setting ``NUMBA_CPU_NAME=generic``, a generic CPU model is picked
    for the CPU architecture and the feature list (``NUMBA_CPU_FEATURES``)
    defaults to empty.  CPU features must be listed with the format
    ``+feature1,-feature2`` where ``+`` indicates enable and ``-`` indicates
    disable. For example, ``+sse,+sse2,-avx,-avx2`` enables SSE and SSE2, and
    disables AVX and AVX2.

    These settings are passed to LLVM for configuring the compilation target.
    To get a list of available options, use the ``llc`` commandline tool
    from LLVM, for example::

        llc -march=x86 -mattr=help


    .. tip:: To force all caching functions (``@jit(cache=True)``) to emit
        portable code (portable within the same architecture and OS),
        simply set ``NUMBA_CPU_NAME=generic``.

.. envvar:: NUMBA_FUNCTION_CACHE_SIZE

    Override the size of the function cache for retaining recently
    deserialized functions in memory.  In systems like
    `Dask <http://dask.pydata.org>`_, it is common for functions to be deserialized
    multiple times.  Numba will cache functions as long as there is a
    reference somewhere in the interpreter.  This cache size variable controls
    how many functions that are no longer referenced will also be retained,
    just in case they show up in the future.  The implementation of this is
    not a true LRU, but the large size of the cache should be sufficient for
    most situations.

    Note: this is unrelated to the compilation cache.

    *Default value:* 128

.. envvar:: NUMBA_LLVM_REFPRUNE_PASS

    Turns on the LLVM pass level reference-count pruning pass and disables the
    regex based implementation in Numba.

    *Default value:* 1 (On)

.. envvar:: NUMBA_LLVM_REFPRUNE_FLAGS

    When ``NUMBA_LLVM_REFPRUNE_PASS`` is on, this allows configuration
    of subpasses in the reference-count pruning LLVM pass.

    Valid values are any combinations of the below separated by `,`
    (case-insensitive):

    - ``all``: enable all subpasses.
    - ``per_bb``: enable per-basic-block level pruning, which is same as the
      old regex based implementation.
    - ``diamond``: enable inter-basic-block pruning that is a diamond shape
      pattern, i.e. a single-entry single-exit CFG subgraph where has an incref
      in the entry and a corresponding decref in the exit.
    - ``fanout``: enable inter-basic-block pruning that has a fanout pattern,
      i.e. a single-entry multiple-exit CFG subgraph where the entry has an
      incref and every exit has a corresponding decref.
    - ``fanout_raise``: same as ``fanout`` but allow subgraph exit nodes to be
      raising an exception and not have a corresponding decref.

    For example, ``all`` is the same as
    ``per_bb, diamond, fanout, fanout_raise``

    *Default value:* "all"


.. _numba-envvars-caching:

Caching options
---------------

Options for the compilation cache.

.. envvar:: NUMBA_DEBUG_CACHE

   If set to non-zero, print out information about operation of the
   :ref:`JIT compilation cache <jit-cache>`.

.. envvar:: NUMBA_CACHE_DIR

    Override the location of the cache directory. If defined, this should be
    a valid directory path.

    If not defined, Numba picks the cache directory in the following order:

    1. In-tree cache. Put the cache next to the corresponding source file under
       a ``__pycache__`` directory following how ``.pyc`` files are stored.
    2. User-wide cache. Put the cache in the user's application directory using
       ``appdirs.user_cache_dir`` from the
       `Appdirs package <https://github.com/ActiveState/appdirs>`_.
    3. IPython cache. Put the cache in an IPython specific application
       directory.
       Stores are made under the ``numba_cache`` in the directory returned by
       ``IPython.paths.get_ipython_cache_dir()``.

    Also see :ref:`docs on cache sharing <cache-sharing>` and
    :ref:`docs on cache clearing <cache-clearing>`


.. _numba-envvars-gpu-support:

GPU support
-----------

.. envvar:: NUMBA_DISABLE_CUDA

   If set to non-zero, disable CUDA support.

.. envvar:: NUMBA_FORCE_CUDA_CC

   If set, force the CUDA compute capability to the given version (a
   string of the type ``major.minor``), regardless of attached devices.

.. envvar:: NUMBA_CUDA_DEFAULT_PTX_CC

   The default compute capability (a string of the type ``major.minor``) to
   target when compiling to PTX using ``cuda.compile_ptx``. The default is
   5.2, which is the lowest non-deprecated compute capability in the most
   recent version of the CUDA toolkit supported (10.2 at present).

.. envvar:: NUMBA_ENABLE_CUDASIM

   If set, don't compile and execute code for the GPU, but use the CUDA
   Simulator instead. For debugging purposes.


.. envvar:: NUMBA_CUDA_ARRAY_INTERFACE_SYNC

   Whether to synchronize on streams provided by objects imported using the CUDA
   Array Interface. This defaults to 1. If set to 0, then no synchronization
   takes place, and the user of Numba (and other CUDA libraries) is responsible
   for ensuring correctness with respect to synchronization on streams.

.. envvar:: NUMBA_CUDA_LOG_LEVEL

   For debugging purposes. If no other logging is configured, the value of this
   variable is the logging level for CUDA API calls. The default value is
   ``CRITICAL`` - to trace all API calls on standard error, set this to
   ``DEBUG``.

.. envvar:: NUMBA_CUDA_LOG_API_ARGS

   By default the CUDA API call logs only give the names of functions called.
   Setting this variable to 1 also includes the values of arguments to Driver
   API calls in the logs.

.. envvar:: NUMBA_CUDA_DRIVER

   Path of the directory in which the CUDA driver libraries are to be found.
   Normally this should not need to be set as Numba can locate the driver in
   standard locations. However, this variable can be used if the driver is in a
   non-standard location.

.. envvar:: NUMBA_CUDA_LOG_SIZE

   Buffer size for logs produced by CUDA driver API operations. This defaults
   to 1024 and should not normally need to be modified - however, if an error
   in an API call produces a large amount of output that appears to be
   truncated (perhaps due to multiple long function names, for example) then
   this variable can be used to increase the buffer size and view the full
   error message.

.. envvar:: NUMBA_CUDA_VERBOSE_JIT_LOG

   Whether the CUDA driver should produce verbose log messages. Defaults to 1,
   indicating that verbose messaging is enabled. This should not need to be
   modified under normal circumstances.

.. envvar:: NUMBA_CUDA_PER_THREAD_DEFAULT_STREAM

   When set to 1, the default stream is the per-thread default stream. When set
   to 0, the default stream is the legacy default stream. This defaults to 0,
   for the legacy default stream. See `Stream Synchronization Behavior
   <https://docs.nvidia.com/cuda/cuda-runtime-api/stream-sync-behavior.html>`_
   for an explanation of the legacy and per-thread default streams.

   This variable only takes effect when using Numba's internal CUDA bindings;
   when using the NVIDIA bindings, use the environment variable
   ``CUDA_PYTHON_CUDA_PER_THREAD_DEFAULT_STREAM`` instead.

   .. seealso::

      The `Default Stream section
      <https://nvidia.github.io/cuda-python/release/11.6.0-notes.html#default-stream>`_
      in the NVIDIA Bindings documentation.

.. envvar:: NUMBA_CUDA_LOW_OCCUPANCY_WARNINGS

   Enable warnings if the grid size is too small relative to the number of
   streaming multiprocessors (SM). This option is on by default (default value is 1).

   The heuristic checked is whether ``gridsize < 2 * (number of SMs)``. NOTE: The absence of
   a warning does not imply a good gridsize relative to the number of SMs. Disabling
   this warning will reduce the number of CUDA API calls (during JIT compilation), as the
   heuristic needs to check the number of SMs available on the device in the
   current context.

.. envvar:: NUMBA_CUDA_WARN_ON_IMPLICIT_COPY

   Enable warnings if a kernel is launched with host memory which forces a copy to and
   from the device. This option is on by default (default value is 1).

.. envvar:: NUMBA_CUDA_USE_NVIDIA_BINDING

   When set to 1, Numba will attempt to use the `NVIDIA CUDA Python binding
   <https://nvidia.github.io/cuda-python/>`_ to make calls to the driver API
   instead of using its own ctypes binding. This defaults to 0 (off), as the
   NVIDIA binding is currently missing support for Per-Thread Default
   Streams and the profiler APIs.

Threading Control
-----------------

.. envvar:: NUMBA_NUM_THREADS

   If set, the number of threads in the thread pool for the parallel CPU target
   will take this value. Must be greater than zero. This value is independent
   of ``OMP_NUM_THREADS`` and ``MKL_NUM_THREADS``.

   *Default value:* The number of CPU cores on the system as determined at run
   time. This can be accessed via :obj:`numba.config.NUMBA_DEFAULT_NUM_THREADS`.

   See also the section on :ref:`setting_the_number_of_threads` for
   information on how to set the number of threads at runtime.

.. envvar:: NUMBA_THREADING_LAYER

   This environment variable controls the library used for concurrent execution
   for the CPU parallel targets (``@vectorize(target='parallel')``,
   ``@guvectorize(target='parallel')``  and ``@njit(parallel=True)``). The
   variable type is string and by default is ``default`` which will select a
   threading layer based on what is available in the runtime. The valid values
   are (for more information about these see
   :ref:`the threading layer documentation <numba-threading-layer>`):

   * ``default`` - select a threading layer based on what is available in the
     current runtime.
   * ``safe`` - select a threading layer that is both fork and thread safe
     (requires the TBB package).
   * ``forksafe`` - select a threading layer that is fork safe.
   * ``threadsafe`` - select a threading layer that is thread safe.
   * ``tbb`` - A threading layer backed by Intel TBB.
   * ``omp`` - A threading layer backed by OpenMP.
   * ``workqueue`` - A simple built-in work-sharing task scheduler.

.. envvar:: NUMBA_THREADING_LAYER_PRIORITY

   This environment variable controls the order in which the libraries used for
   concurrent execution, for the CPU parallel targets
   (``@vectorize(target='parallel')``, ``@guvectorize(target='parallel')``
   and ``@njit(parallel=True)``), are prioritized for use. The variable type is
   string and by default is ``tbb omp workqueue``, with the priority taken based
   on position from the left of the string, left most being the highest. Valid
   values are any permutation of the three choices (for more information about
   these see :ref:`the threading layer documentation <numba-threading-layer>`.)

Reference Manual
================

.. toctree::

   types.rst
   jit-compilation.rst
   aot-compilation.rst
   utils.rst
   envvars.rst
   pysupported.rst
   numpysupported.rst
   pysemantics.rst
   fpsemantics.rst
   deprecation.rst

Floating-point pitfalls
=======================

Precision and accuracy
----------------------

For some operations, Numba may use a different algorithm than Python or
Numpy.  The results may not be bit-by-bit compatible.  The difference
should generally be small and within reasonable expectations.  However,
small accumulated differences might produce large differences at the end,
especially if a divergent function is involved.

Math library implementations
''''''''''''''''''''''''''''

Numba supports a variety of platforms and operating systems, each of which
has its own math library implementation (referred to as ``libm`` from here
in).  The majority of math functions included in ``libm`` have specific
requirements as set out by the IEEE 754 standard (like ``sin()``, ``exp()``
etc.), but each implementation may have bugs.  Thus, on some platforms
Numba has to exercise special care in order to workaround known ``libm``
issues.

Another typical problem is when an operating system's ``libm`` function
set is incomplete and needs to be supplemented by additional functions.
These are provided with reference to the IEEE 754 and C99 standards
and are often implemented in Numba in a manner similar to equivalent
CPython functions.

Linear algebra
''''''''''''''

Numpy forces some linear algebra operations to run in double-precision mode
even when a ``float32`` input is given.  Numba will always observe
the input's precision, and invoke single-precision linear algebra routines
when all inputs are ``float32`` or ``complex64``.

The implementations of the ``numpy.linalg`` routines in Numba only support the
floating point types that are used in the LAPACK functions that provide
the underlying core functionality. As a result  only ``float32``, ``float64``,
``complex64`` and ``complex128`` types are supported. If a user has e.g. an
``int32`` type, an appropriate type conversion must be performed to a
floating point type prior to its use in these routines. The reason for this
decision is to essentially avoid having to replicate type conversion choices
made in Numpy and to also encourage the user to choose the optimal floating
point type for the operation they are undertaking.


Mixed-types operations
''''''''''''''''''''''

Numpy will most often return a ``float64`` as a result of a computation
with mixed integer and floating-point operands (a typical example is the
power operator ``**``).  Numba by contrast will select the highest precision
amongst the floating-point operands, so for example ``float32 ** int32``
will return a ``float32``, regardless of the input values.  This makes
performance characteristics easier to predict, but you should explicitly
cast the input to ``float64`` if you need the extra precision.


.. _ufunc-fpu-errors:

Warnings and errors
-------------------

When calling a :term:`ufunc` created with :func:`~numba.vectorize`,
Numpy will determine whether an error occurred by examining the FPU
error word.  It may then print out a warning or raise an exception
(such as ``RuntimeWarning: divide by zero encountered``),
depending on the current error handling settings.

Depending on how LLVM optimized the ufunc's code, however, some spurious
warnings or errors may appear.  If you get caught by this issue, we
recommend you call :func:`numpy.seterr` to change Numpy's error handling
settings, or the :class:`numpy.errstate` context manager to switch them
temporarily::

   with np.errstate(all='ignore'):
       x = my_ufunc(y)

.. _pysemantics:

Deviations from Python Semantics
================================

Bounds Checking
---------------

By default, instead of causing an :class:`IndexError`, accessing an
out-of-bound index of an array in a Numba-compiled function will return
invalid values or lead to an access violation error (it's reading from
invalid memory locations). Bounds checking can be enabled on a specific
function via the :ref:`boundscheck <jit-decorator-boundscheck>`
option of the jit decorator. Additionally, the :envvar:`NUMBA_BOUNDSCHECK`
can be set to 0 or 1 to globally override this flag.

.. note::
  Bounds checking will slow down typical functions so it is recommended to only
  use this flag for debugging purposes.

Exceptions and Memory Allocation
--------------------------------

Due to limitations in the current compiler when handling exceptions, memory
allocated (almost always NumPy arrays) within a function that raises an
exception will **leak**.  This is a known issue that will be fixed, but in the
meantime, it is best to do memory allocation outside of functions that can
also raise exceptions.

Integer width
-------------

While Python has arbitrary-sized integers, integers in Numba-compiled
functions get a fixed size through :term:`type inference` (usually,
the size of a machine integer).  This means that arithmetic
operations can wrapround or produce undefined results or overflow.

Type inference can be overridden by an explicit type specification,
if fine-grained control of integer width is desired.

.. seealso::
   :ref:`Enhancement proposal 1: Changes in integer typing <nbep-1>`


Boolean inversion
-----------------

Calling the bitwise complement operator (the ``~`` operator) on a Python
boolean returns an integer, while the same operator on a Numpy boolean
returns another boolean::

   >>> ~True
   -2
   >>> ~np.bool_(True)
   False

Numba follows the Numpy semantics.


Global and closure variables
----------------------------

In :term:`nopython mode`, global and closure variables are *frozen* by
Numba: a Numba-compiled function sees the value of those variables at the
time the function was compiled.  Also, it is not possible to change their
values from the function.

Numba **may or may not** copy global variables referenced inside a compiled
function.  Small global arrays are copied for potential compiler optimization
with immutability assumption.  However, large global arrays are not copied to
conserve memory.  The definition of "small" and "large" may change.


Zero initialization of variables
--------------------------------

Numba does not track variable liveness at runtime. For simplicity of
implementation, all variables are zero-initialized. Example::

    from numba import njit

    @njit
    def foo():
        for i in range(0):
            pass
        print(i) # will print 0 and not raise UnboundLocalError

    foo()
=================
Memory management
=================

.. _cuda-device-memory:

Data transfer
=============

Even though Numba can automatically transfer NumPy arrays to the device,
it can only do so conservatively by always transferring device memory back to
the host when a kernel finishes. To avoid the unnecessary transfer for
read-only arrays, you can use the following APIs to manually control the
transfer:

.. autofunction:: numba.cuda.device_array
   :noindex:
.. autofunction:: numba.cuda.device_array_like
   :noindex:
.. autofunction:: numba.cuda.to_device
   :noindex:

In addition to the device arrays, Numba can consume any object that implements
:ref:`cuda array interface <cuda-array-interface>`.  These objects also can be
manually converted into a Numba device array by creating a view of the GPU
buffer using the following APIs:

.. autofunction:: numba.cuda.as_cuda_array
  :noindex:
.. autofunction:: numba.cuda.is_cuda_array
  :noindex:


Device arrays
-------------

Device array references have the following methods.  These methods are to be
called in host code, not within CUDA-jitted functions.

.. autoclass:: numba.cuda.cudadrv.devicearray.DeviceNDArray
    :members: copy_to_host, is_c_contiguous, is_f_contiguous, ravel, reshape
    :noindex:


.. note:: DeviceNDArray defines the :ref:`cuda array interface <cuda-array-interface>`.


Pinned memory
=============

.. autofunction:: numba.cuda.pinned
   :noindex:
.. autofunction:: numba.cuda.pinned_array
   :noindex:
.. autofunction:: numba.cuda.pinned_array_like
   :noindex:


Mapped memory
=============

.. autofunction:: numba.cuda.mapped
   :noindex:
.. autofunction:: numba.cuda.mapped_array
   :noindex:
.. autofunction:: numba.cuda.mapped_array_like
   :noindex:



Managed memory
==============

.. autofunction:: numba.cuda.managed_array
   :noindex:


Streams
=======

Streams can be passed to functions that accept them (e.g. copies between the
host and device) and into kernel launch configurations so that the operations
are executed asynchronously.

.. autofunction:: numba.cuda.stream
   :noindex:

.. autofunction:: numba.cuda.default_stream
   :noindex:

.. autofunction:: numba.cuda.legacy_default_stream
   :noindex:

.. autofunction:: numba.cuda.per_thread_default_stream
   :noindex:

.. autofunction:: numba.cuda.external_stream
   :noindex:

CUDA streams have the following methods:

.. autoclass:: numba.cuda.cudadrv.driver.Stream
    :members: synchronize, auto_synchronize
    :noindex:

.. _cuda-shared-memory:

Shared memory and thread synchronization
========================================

A limited amount of shared memory can be allocated on the device to speed
up access to data, when necessary.  That memory will be shared (i.e. both
readable and writable) amongst all threads belonging to a given block
and has faster access times than regular device memory.  It also allows
threads to cooperate on a given solution.  You can think of it as a
manually-managed data cache.

The memory is allocated once for the duration of the kernel, unlike
traditional dynamic memory management.

.. function:: numba.cuda.shared.array(shape, type)
   :noindex:

   Allocate a shared array of the given *shape* and *type* on the device.
   This function must be called on the device (i.e. from a kernel or
   device function). *shape* is either an integer or a tuple of integers
   representing the array's dimensions and must be a simple constant
   expression. *type* is a :ref:`Numba type <numba-types>` of the elements
   needing to be stored in the array.

   The returned array-like object can be read and written to like any normal
   device array (e.g. through indexing).

   A common pattern is to have each thread populate one element in the
   shared array and then wait for all threads to finish using :func:`.syncthreads`.


.. function:: numba.cuda.syncthreads()
   :noindex:

   Synchronize all threads in the same thread block.  This function
   implements the same pattern as `barriers <http://en.wikipedia.org/wiki/Barrier_%28computer_science%29>`_
   in traditional multi-threaded programming: this function waits
   until all threads in the block call it, at which point it returns
   control to all its callers.

.. seealso::
   :ref:`Matrix multiplication example <cuda-matmul>`.

.. _cuda-local-memory:

Local memory
============

Local memory is an area of memory private to each thread.  Using local
memory helps allocate some scratchpad area when scalar local variables
are not enough.  The memory is allocated once for the duration of the kernel,
unlike traditional dynamic memory management.

.. function:: numba.cuda.local.array(shape, type)
   :noindex:

   Allocate a local array of the given *shape* and *type* on the device.
   *shape* is either an integer or a tuple of integers representing the array's
   dimensions and must be a simple constant expression. *type* is a :ref:`Numba
   type <numba-types>` of the elements needing to be stored in the array. The
   array is private to the current thread. An array-like object is returned
   which can be read and written to like any standard array (e.g. through
   indexing).

   .. seealso:: The Local Memory section of `Device Memory Accesses
      <https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#device-memory-accesses>`_
      in the CUDA programming guide.

Constant memory
===============

Constant memory is an area of memory that is read only, cached and off-chip, it
is accessible by all threads and is host allocated. A method of
creating an array in constant memory is through the use of:

.. function:: numba.cuda.const.array_like(arr)
   :noindex:

   Allocate and make accessible an array in constant memory based on array-like
   *arr*.


.. _deallocation-behavior:

Deallocation Behavior
=====================

This section describes the deallocation behaviour of Numba's internal memory
management. If an External Memory Management Plugin is in use (see
:ref:`cuda-emm-plugin`), then deallocation behaviour may differ; you may refer to the
documentation for the EMM Plugin to understand its deallocation behaviour.

Deallocation of all CUDA resources are tracked on a per-context basis.
When the last reference to a device memory is dropped, the underlying memory
is scheduled to be deallocated.  The deallocation does not occur immediately.
It is added to a queue of pending deallocations.  This design has two benefits:

1. Resource deallocation API may cause the device to synchronize; thus, breaking
   any asynchronous execution.  Deferring the deallocation could avoid latency
   in performance critical code section.
2. Some deallocation errors may cause all the remaining deallocations to fail.
   Continued deallocation errors can cause critical errors at the CUDA driver
   level.  In some cases, this could mean a segmentation fault in the CUDA
   driver. In the worst case, this could cause the system GUI to freeze and
   could only recover with a system reset.  When an error occurs during a
   deallocation, the remaining pending deallocations are cancelled.  Any
   deallocation error will be reported.  When the process is terminated, the
   CUDA driver is able to release all allocated resources by the terminated
   process.

The deallocation queue is flushed automatically as soon as the following events
occur:

- An allocation failed due to out-of-memory error.  Allocation is retried after
  flushing all deallocations.
- The deallocation queue has reached its maximum size, which is default to 10.
  User can override by setting the environment variable
  `NUMBA_CUDA_MAX_PENDING_DEALLOCS_COUNT`.  For example,
  `NUMBA_CUDA_MAX_PENDING_DEALLOCS_COUNT=20`, increases the limit to 20.
- The maximum accumulated byte size of resources that are pending deallocation
  is reached.  This is default to 20% of the device memory capacity.
  User can override by setting the environment variable
  `NUMBA_CUDA_MAX_PENDING_DEALLOCS_RATIO`. For example,
  `NUMBA_CUDA_MAX_PENDING_DEALLOCS_RATIO=0.5` sets the limit to 50% of the
  capacity.

Sometimes, it is desired to defer resource deallocation until a code section
ends.  Most often, users want to avoid any implicit synchronization due to
deallocation.  This can be done by using the following context manager:

.. autofunction:: numba.cuda.defer_cleanup
GPU Reduction
==============

Writing a reduction algorithm for CUDA GPU can be tricky.  Numba provides a
``@reduce`` decorator for converting a simple binary operation into a reduction
kernel. An example follows::

    import numpy
    from numba import cuda

    @cuda.reduce
    def sum_reduce(a, b):
        return a + b

    A = (numpy.arange(1234, dtype=numpy.float64)) + 1
    expect = A.sum()      # numpy sum reduction
    got = sum_reduce(A)   # cuda sum reduction
    assert expect == got

Lambda functions can also be used here::

    sum_reduce = cuda.reduce(lambda a, b: a + b)

The Reduce class
----------------

The ``reduce`` decorator creates an instance of the ``Reduce`` class.
Currently, ``reduce`` is an alias to ``Reduce``, but this behavior is not
guaranteed.

.. autoclass:: numba.cuda.Reduce
   :members: __init__, __call__
   :member-order: bysource
.. _cuda-array-interface:

================================
CUDA Array Interface (Version 3)
================================

The *CUDA Array Interface* (or CAI) is created for interoperability between
different implementations of CUDA array-like objects in various projects. The
idea is borrowed from the `NumPy array interface`_.


.. note::
    Currently, we only define the Python-side interface.  In the future, we may
    add a C-side interface for efficient exchange of the information in
    compiled code.


Python Interface Specification
==============================

.. note:: Experimental feature.  Specification may change.

The ``__cuda_array_interface__`` attribute returns a dictionary (``dict``)
that must contain the following entries:

- **shape**: ``(integer, ...)``

  A tuple of ``int`` (or ``long``) representing the size of each dimension.

- **typestr**: ``str``

  The type string.  This has the same definition as ``typestr`` in the
  `numpy array interface`_.

- **data**: ``(integer, boolean)``

  The **data** is a 2-tuple.  The first element is the data pointer
  as a Python ``int`` (or ``long``).  The data must be device-accessible.
  For zero-size arrays, use ``0`` here.
  The second element is the read-only flag as a Python ``bool``.

  Because the user of the interface may or may not be in the same context,
  the most common case is to use ``cuPointerGetAttribute`` with
  ``CU_POINTER_ATTRIBUTE_DEVICE_POINTER`` in the CUDA driver API (or the
  equivalent CUDA Runtime API) to retrieve a device pointer that
  is usable in the currently active context.

- **version**: ``integer``

  An integer for the version of the interface being exported.
  The current version is *3*.


The following are optional entries:

- **strides**: ``None`` or ``(integer, ...)``

  If **strides** is not given, or it is ``None``, the array is in
  C-contiguous layout. Otherwise, a tuple of ``int`` (or ``long``) is explicitly
  given for representing the number of bytes to skip to access the next
  element at each dimension.

- **descr**

  This is for describing more complicated types.  This follows the same
  specification as in the `numpy array interface`_.

- **mask**: ``None`` or object exposing the ``__cuda_array_interface__``

  If ``None`` then all values in **data** are valid. All elements of the mask
  array should be interpreted only as true or not true indicating which
  elements of this array are valid. This has the same definition as ``mask``
  in the `numpy array interface`_.

  .. note:: Numba does not currently support working with masked CUDA arrays
            and will raise a ``NotImplementedError`` exception if one is passed
            to a GPU function.

- **stream**: ``None`` or ``integer``

  An optional stream upon which synchronization must take place at the point of
  consumption, either by synchronizing on the stream or enqueuing operations on
  the data on the given stream. Integer values in this entry are as follows:

  - ``0``: This is disallowed as it would be ambiguous between ``None`` and the
    default stream, and also between the legacy and per-thread default streams.
    Any use case where ``0`` might be given should either use ``None``, ``1``,
    or ``2`` instead for clarity.
  - ``1``: The legacy default stream.
  - ``2``: The per-thread default stream.
  - Any other integer: a ``cudaStream_t`` represented as a Python integer.

  When ``None``, no synchronization is required. See the
  :ref:`cuda-array-interface-synchronization` section below for further details.

  In a future revision of the interface, this entry may be expanded (or another
  entry added) so that an event to synchronize on can be specified instead of a
  stream.


.. _cuda-array-interface-synchronization:

Synchronization
---------------

Definitions
~~~~~~~~~~~

When discussing synchronization, the following definitions are used:

- *Producer*: The library / object on which ``__cuda_array_interface__`` is
  accessed.
- *Consumer*: The library / function that accesses the
  ``__cuda_array_interface__`` of the Producer.
- *User Code*: Code that induces a Producer and Consumer to share data through
  the CAI.
- *User*: The person writing or maintaining the User Code. The User may
  implement User Code without knowledge of the CAI, since the CAI accesses can
  be hidden from their view.

In the following example:

.. code-block:: python

   import cupy
   from numba import cuda

   @cuda.jit
   def add(x, y, out):
       start = cuda.grid(1)
       stride = cuda.gridsize(1)
       for i in range(start, x.shape[0], stride):
           out[i] = x[i] + y[i]

   a = cupy.arange(10)
   b = a * 2
   out = cupy.zeros_like(a)

   add[1, 32](a, b, out)

When the ``add`` kernel is launched:

- ``a``, ``b``, ``out`` are Producers.
- The ``add`` kernel is the Consumer.
- The User Code is specifically ``add[1, 32](a, b, out)``.
- The author of the code is the User.


Design Motivations
~~~~~~~~~~~~~~~~~~

Elements of the CAI design related to synchronization seek to fulfill these
requirements:

1. Producers and Consumers that exchange data through the CAI must be able to do
   so without data races.
2. Requirement 1 should be met without requiring the user to be
   aware of any particulars of the CAI - in other words, exchanging data between
   Producers and Consumers that operate on data asynchronously should be correct
   by default.

   - An exception to this requirement is made for Producers and Consumers that
     explicitly document that the User is required to take additional steps to
     ensure correctness with respect to synchronization. In this case, Users
     are required to understand the details of the CUDA Array Interface, and
     the Producer/Consumer library documentation must specify the steps that
     Users are required to take.

     Use of this exception should be avoided where possible, as it is provided
     for libraries that cannot implement the synchronization semantics without
     the involvement of the User - for example, those interfacing with
     third-party libraries oblivious to the CUDA Array Interface.

3. Where the User is aware of the particulars of the CAI and implementation
   details of the Producer and Consumer, they should be able to, at their
   discretion, override some of the synchronization semantics of the interface
   to reduce the synchronization overhead. Overriding synchronization semantics
   implies that:

   - The CAI design, and the design and implementation of the Producer and
     Consumer do not specify or guarantee correctness with respect to data
     races.
   - Instead, the User is responsible for ensuring correctness with respect to
     data races.


Interface Requirements
~~~~~~~~~~~~~~~~~~~~~~

The ``stream`` entry enables Producers and Consumers to avoid hazards when
exchanging data. Expected behaviour of the Consumer is as follows:

* When ``stream`` is not present or is ``None``:

  - No synchronization is required on the part of the Consumer.
  - The Consumer may enqueue operations on the underlying data immediately on
    any stream.

* When ``stream`` is an integer, its value indicates the stream on which the
  Producer may have in-progress operations on the data, and which the Consumer
  is expected to either:

  - Synchronize on before accessing the data, or
  - Enqueue operations in when accessing the data.

  The Consumer can choose which mechanism to use, with the following
  considerations:

  - If the Consumer synchronizes on the provided stream prior to accessing the
    data, then it must ensure that no computation can take place in the provided
    stream until its operations in its own choice of stream have taken place.
    This could be achieved by either:

    - Placing a wait on an event in the provided stream that occurs once all
      of the Consumer's operations on the data are completed, or
    - Avoiding returning control to the user code until after its operations
      on its own stream have completed.

  - If the consumer chooses to only enqueue operations on the data in the
    provided stream, then it may return control to the User code immediately
    after enqueueing its work, as the work will all be serialized on the
    exported array's stream. This is sufficient to ensure correctness even if
    the User code were to induce the Producer to subsequently start enqueueing
    more work on the same stream.

* If the User has set the Consumer to ignore CAI synchronization semantics, the
  Consumer may assume it can operate on the data immediately in any stream with
  no further synchronization, even if the ``stream`` member has an integer
  value.


When exporting an array through the CAI, Producers must ensure that:

* If there is work on the data enqueued in one or more streams, then
  synchronization on the provided ``stream`` is sufficient to ensure
  synchronization with all pending work.

  - If the Producer has no enqueued work, or work only enqueued on the stream
    identified by ``stream``, then this condition is met.
  - If the Producer has enqueued work on the data on multiple streams, then it
    must enqueue events on those streams that follow the enqueued work, and
    then wait on those events in the provided ``stream``. For example:

    1. Work is enqueued by the Producer on streams ``7``, ``9``, and ``15``.
    2. Events are then enqueued on each of streams ``7``, ``9``, and ``15``.
    3. Producer then tells stream ``3`` to wait on the events from Step 2, and
       the ``stream`` entry is set to ``3``.

* If there is no work enqueued on the data, then the ``stream`` entry may be
  either ``None``, or not provided.

Optionally, to facilitate the User relaxing conformance to synchronization
semantics:

* Producers may provide a configuration option to always set ``stream`` to
  ``None``.
* Consumers may provide a configuration option to ignore the value of ``stream``
  and act as if it were ``None`` or not provided.  This elides synchronization
  on the Producer-provided streams, and allows enqueuing work on streams other
  than that provided by the Producer.

These options should not be set by default in either a Producer or a Consumer.
The CAI specification does not prescribe the exact mechanism by which these
options are set, or related options that Producers or Consumers might provide
to allow the user further control over synchronization behavior.


Synchronization in Numba
~~~~~~~~~~~~~~~~~~~~~~~~

Numba is neither strictly a Producer nor a Consumer - it may be used to
implement either by a User. In order to facilitate the correct implementation of
synchronization semantics, Numba exhibits the following behaviors related to
synchronization of the interface:

- When Numba acts as a Consumer (for example when an array-like object is passed
  to a kernel launch): If ``stream`` is an integer, then Numba will immediately
  synchronize on the provided ``stream``. A Numba :class:`Device Array
  <numba.cuda.cudadrv.devicearray.DeviceNDArray>` created from an array-like
  object has its *default stream* set to the provided stream.

- When Numba acts as a Producer (when the ``__cuda_array_interface__`` property
  of a Numba CUDA Array is accessed): If the exported CUDA Array has a
  *default stream*, then it is given as the ``stream`` entry. Otherwise,
  ``stream`` is set to ``None``.

.. note:: In Numba's terminology, an array's *default stream* is a property
          specifying the stream that Numba will enqueue asynchronous
          transfers in if no other stream is provided as an argument to the
          function invoking the transfer. It is not the same as the `Default
          Stream
          <https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#default-stream>`_
          in normal CUDA terminology.

Numba's synchronization behavior results in the following intended
consequences:

- Exchanging data either as a Producer or a Consumer will be correct without
  the need for any further action from the User, provided that the other side
  of the interaction also follows the CAI synchronization semantics.
- The User is expected to either:

  - Avoid launching kernels or other operations on streams that
    are not the default stream for their parameters, or
  - When launching operations on a stream that is not the default stream for
    a given parameter, they should then insert an event into the stream that
    they are operating in, and wait on that event in the default stream for
    the parameter. For an example of this, :ref:`see below
    <example-multi-streams>`.

The User may override Numba's synchronization behavior by setting the
environment variable ``NUMBA_CUDA_ARRAY_INTERFACE_SYNC`` or the config variable
``CUDA_ARRAY_INTERFACE_SYNC`` to ``0`` (see :ref:`GPU Support Environment
Variables <numba-envvars-gpu-support>`).  When set, Numba will not synchronize
on the streams of imported arrays, and it is the responsibility of the user to
ensure correctness with respect to stream synchronization. Synchronization when
creating a Numba CUDA Array from an object exporting the CUDA Array Interface
may also be elided by passing ``sync=False`` when creating the Numba CUDA
Array with :func:`numba.cuda.as_cuda_array` or
:func:`numba.cuda.from_cuda_array_interface`.

There is scope for Numba's synchronization implementation to be optimized in
the future, by eliding synchronizations when a kernel or driver API operation
(e.g.  a memcopy or memset) is launched on the same stream as an imported
array.


.. _example-multi-streams:

An example launching on an array's non-default stream
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This example shows how to ensure that a Consumer can safely consume an array
with a default stream when it is passed to a kernel launched in a different
stream.

First we need to import Numba and a consumer library (a fictitious library named
``other_cai_library`` for this example):

.. code-block:: python

   from numba import cuda, int32, void
   import other_cai_library

Now we'll define a kernel - this initializes the elements of the array, setting
each entry to its index:

.. code-block:: python

   @cuda.jit(void, int32[::1])
   def initialize_array(x):
       i = cuda.grid(1)
       if i < len(x):
           x[i] = i

Next we will create two streams:

.. code-block:: python

   array_stream = cuda.stream()
   kernel_stream = cuda.stream()

Then create an array with one of the streams as its default stream:

.. code-block:: python

   N = 16384
   x = cuda.device_array(N, stream=array_stream)

Now we launch the kernel in the other stream:

.. code-block:: python

   nthreads = 256
   nblocks = N // nthreads

   initialize_array[nthreads, nblocks, kernel_stream](x)

If we were to pass ``x`` to a Consumer now, there is a risk that it may operate on
it in ``array_stream`` whilst the kernel is still running in ``kernel_stream``.
To prevent operations in ``array_stream`` starting before the kernel launch is
finished, we create an event and wait on it:

.. code-block:: python

   # Create event
   evt = cuda.event()
   # Record the event after the kernel launch in kernel_stream
   evt.record(kernel_stream)
   # Wait for the event in array_stream
   evt.wait(array_stream)

It is now safe for ``other_cai_library`` to consume ``x``:

.. code-block:: python

   other_cai_library.consume(x)


Lifetime management
-------------------

Data
~~~~

Obtaining the value of the ``__cuda_array_interface__`` property of any object
has no effect on the lifetime of the object from which it was created. In
particular, note that the interface has no slot for the owner of the data.

The User code must preserve the lifetime of the object owning the data for as
long as the Consumer might use it.


Streams
~~~~~~~

Like data, CUDA streams also have a finite lifetime. It is therefore required
that a Producer exporting data on the interface with an associated stream
ensures that the exported stream's lifetime is equal to or surpasses the
lifetime of the object from which the interface was exported.


Lifetime management in Numba
----------------------------

Producing Arrays
~~~~~~~~~~~~~~~~

Numba takes no steps to maintain the lifetime of an object from which the
interface is exported - it is the user's responsibility to ensure that the
underlying object is kept alive for the duration that the exported interface
might be used.

The lifetime of any Numba-managed stream exported on the interface is guaranteed
to equal or surpass the lifetime of the underlying object, because the
underlying object holds a reference to the stream.

.. note:: Numba-managed streams are those created with
          ``cuda.default_stream()``, ``cuda.legacy_default_stream()``, or
          ``cuda.per_thread_default_stream()``. Streams not managed by Numba
          are created from an external stream with ``cuda.external_stream()``.


Consuming Arrays
~~~~~~~~~~~~~~~~

Numba provides two mechanisms for creating device arrays from objects exporting
the CUDA Array Interface. Which to use depends on whether the created device
array should maintain the life of the object from which it is created:

- ``as_cuda_array``: This creates a device array that holds a reference to the
  owning object. As long as a reference to the device array is held, its
  underlying data will also be kept alive, even if all other references to the
  original owning object have been dropped.
- ``from_cuda_array_interface``: This creates a device array with no reference
  to the owning object by default. The owning object, or some other object to
  be considered the owner can be passed in the ``owner`` parameter.

The interfaces of these functions are:

.. automethod:: numba.cuda.as_cuda_array

.. automethod:: numba.cuda.from_cuda_array_interface


Pointer Attributes
------------------

Additional information about the data pointer can be retrieved using
``cuPointerGetAttribute`` or ``cudaPointerGetAttributes``.  Such information
include:

- the CUDA context that owns the pointer;
- is the pointer host-accessible?
- is the pointer a managed memory?


.. _numpy array interface: https://docs.scipy.org/doc/numpy-1.13.0/reference/arrays.interface.html#__array_interface__


Differences with CUDA Array Interface (Version 0)
-------------------------------------------------

Version 0 of the CUDA Array Interface did not have the optional **mask**
attribute to support masked arrays.


Differences with CUDA Array Interface (Version 1)
-------------------------------------------------

Versions 0 and 1 of the CUDA Array Interface neither clarified the
**strides** attribute for C-contiguous arrays nor specified the treatment for
zero-size arrays.


Differences with CUDA Array Interface (Version 2)
-------------------------------------------------

Prior versions of the CUDA Array Interface made no statement about
synchronization.


Interoperability
----------------

The following Python libraries have adopted the CUDA Array Interface:

- Numba
- `CuPy <https://docs-cupy.chainer.org/en/stable/reference/interoperability.html>`_
- `PyTorch <https://pytorch.org>`_
- `PyArrow <https://arrow.apache.org/docs/python/generated/pyarrow.cuda.Context.html#pyarrow.cuda.Context.buffer_from_object>`_
- `mpi4py <https://mpi4py.readthedocs.io/en/latest/overview.html#support-for-cuda-aware-mpi>`_
- `ArrayViews <https://github.com/xnd-project/arrayviews>`_
- `JAX <https://jax.readthedocs.io/en/latest/index.html>`_
- `PyCUDA <https://documen.tician.de/pycuda/tutorial.html#interoperability-with-other-libraries-using-the-cuda-array-interface>`_
- `DALI: the NVIDIA Data Loading Library <https://github.com/NVIDIA/DALI>`_ :

    - `TensorGPU objects
      <https://docs.nvidia.com/deeplearning/dali/user-guide/docs/data_types.html#nvidia.dali.backend.TensorGPU>`_
      expose the CUDA Array Interface.
    - `The External Source operator
      <https://docs.nvidia.com/deeplearning/dali/user-guide/docs/supported_ops.html#nvidia.dali.fn.external_source>`_
      consumes objects exporting the CUDA Array Interface.
- The RAPIDS stack:

    - `cuDF <https://rapidsai.github.io/projects/cudf/en/0.11.0/10min-cudf-cupy.html>`_
    - `cuML <https://docs.rapids.ai/api/cuml/nightly/>`_
    - `cuSignal <https://github.com/rapidsai/cusignal>`_
    - `RMM <https://docs.rapids.ai/api/rmm/stable/>`_

If your project is not on this list, please feel free to report it on the `Numba issue tracker <https://github.com/numba/numba/issues>`_.

========
Examples
========

.. _cuda-matmul:

Matrix multiplication
=====================
First, import the modules needed for this example:

.. literalinclude:: ../../../numba/cuda/tests/doc_examples/test_matmul.py
   :language: python
   :caption: from ``test_ex_matmul`` in ``numba/cuda/tests/doc_examples/test_matmul.py``
   :start-after: magictoken.ex_import.begin
   :end-before: magictoken.ex_import.end
   :dedent: 8
   :linenos:

Here is a naïve implementation of matrix multiplication using a CUDA kernel:

.. literalinclude:: ../../../numba/cuda/tests/doc_examples/test_matmul.py
   :language: python
   :caption: from ``test_ex_matmul`` in ``numba/cuda/tests/doc_examples/test_matmul.py``
   :start-after: magictoken.ex_matmul.begin
   :end-before: magictoken.ex_matmul.end
   :dedent: 8
   :linenos:

An example usage of this function is as follows:

.. literalinclude:: ../../../numba/cuda/tests/doc_examples/test_matmul.py
   :language: python
   :caption: from ``test_ex_matmul`` in ``numba/cuda/tests/doc_examples/test_matmul.py``
   :start-after: magictoken.ex_run_matmul.begin
   :end-before: magictoken.ex_run_matmul.end
   :dedent: 8
   :linenos:

This implementation is straightforward and intuitive but performs poorly,
because the same matrix elements will be loaded multiple times from device
memory, which is slow (some devices may have transparent data caches, but
they may not be large enough to hold the entire inputs at once).

It will be faster if we use a blocked algorithm to reduce accesses to the
device memory.  CUDA provides a fast :ref:`shared memory <cuda-shared-memory>`
for threads in a block to cooperatively compute on a task.  The following
implements a faster version of the square matrix multiplication using shared
memory:

.. literalinclude:: ../../../numba/cuda/tests/doc_examples/test_matmul.py
   :language: python
   :caption: from ``test_ex_matmul`` in ``numba/cuda/tests/doc_examples/test_matmul.py``
   :start-after: magictoken.ex_fast_matmul.begin
   :end-before: magictoken.ex_fast_matmul.end
   :dedent: 8
   :linenos:


Because the shared memory is a limited resource, the code preloads a small
block at a time from the input arrays.  Then, it calls
:func:`~numba.cuda.syncthreads` to wait until all threads have finished
preloading and before doing the computation on the shared memory.
It synchronizes again after the computation to ensure all threads
have finished with the data in shared memory before overwriting it
in the next loop iteration.

An example usage of the ``fast_matmul`` function is as follows:

.. literalinclude:: ../../../numba/cuda/tests/doc_examples/test_matmul.py
   :language: python
   :caption: from ``test_ex_matmul`` in ``numba/cuda/tests/doc_examples/test_matmul.py``
   :start-after: magictoken.ex_run_fast_matmul.begin
   :end-before: magictoken.ex_run_fast_matmul.end
   :dedent: 8
   :linenos:


This passes a :ref:`CUDA memory check test <debugging-cuda-python-code>`, which
can help with debugging. Running the code above produces the following output:

.. code-block:: none

    $ python fast_matmul.py
    [[ 6.  6.  6.  6.]
    [22. 22. 22. 22.]
    [38. 38. 38. 38.]
    [54. 54. 54. 54.]]
    [[ 6.  6.  6.  6.]
    [22. 22. 22. 22.]
    [38. 38. 38. 38.]
    [54. 54. 54. 54.]]

.. note:: For high performance matrix multiplication in CUDA, see also the `CuPy implementation <https://docs.cupy.dev/en/stable/reference/generated/cupy.matmul.html>`_.

The approach outlined here generalizes to non-square matrix multiplication as
follows by adjusting the ``blockspergrid`` variable:

Again, here is an example usage:

.. literalinclude:: ../../../numba/cuda/tests/doc_examples/test_matmul.py
   :language: python
   :caption: from ``test_ex_matmul`` in ``numba/cuda/tests/doc_examples/test_matmul.py``
   :start-after: magictoken.ex_run_nonsquare.begin
   :end-before: magictoken.ex_run_nonsquare.end
   :dedent: 8
   :linenos:

and the corresponding output:

.. code-block:: none

  $ python nonsquare_matmul.py
  [[ 253.  253.  253.  253.  253.  253.  253.]
  [ 782.  782.  782.  782.  782.  782.  782.]
  [1311. 1311. 1311. 1311. 1311. 1311. 1311.]
  [1840. 1840. 1840. 1840. 1840. 1840. 1840.]
  [2369. 2369. 2369. 2369. 2369. 2369. 2369.]]
  [[ 253.  253.  253.  253.  253.  253.  253.]
  [ 782.  782.  782.  782.  782.  782.  782.]
  [1311. 1311. 1311. 1311. 1311. 1311. 1311.]
  [1840. 1840. 1840. 1840. 1840. 1840. 1840.]
  [2369. 2369. 2369. 2369. 2369. 2369. 2369.]]

Device management
=================

For multi-GPU machines, users may want to select which GPU to use.
By default the CUDA driver selects the fastest GPU as the device 0,
which is the default device used by Numba.

The features introduced on this page are generally not of interest
unless working with systems hosting/offering more than one CUDA-capable GPU.

Device Selection
----------------

If at all required, device selection must be done before any CUDA feature is
used.

::

    from numba import cuda
    cuda.select_device(0)

The device can be closed by:

::

    cuda.close()

Users can then create a new context with another device.

::

    cuda.select_device(1)  # assuming we have 2 GPUs


.. function:: numba.cuda.select_device(device_id)
   :noindex:

   Create a new CUDA context for the selected *device_id*.  *device_id*
   should be the number of the device (starting from 0; the device order
   is determined by the CUDA libraries).  The context is associated with
   the current thread.  Numba currently allows only one context per thread.

   If successful, this function returns a device instance.

   .. XXX document device instances?


.. function:: numba.cuda.close
   :noindex:

   Explicitly close all contexts in the current thread.

   .. note::
      Compiled functions are associated with the CUDA context.
      This makes it not very useful to close and create new devices, though it
      is certainly useful for choosing which device to use when the machine
      has multiple GPUs.

The Device List
===============

The Device List is a list of all the GPUs in the system, and can be indexed to
obtain a context manager that ensures execution on the selected GPU.

.. attribute:: numba.cuda.gpus
   :noindex:
.. attribute:: numba.cuda.cudadrv.devices.gpus

:py:data:`numba.cuda.gpus` is an instance of the ``_DeviceList`` class, from
which the current GPU context can also be retrieved:

.. autoclass:: numba.cuda.cudadrv.devices._DeviceList
    :members: current
    :noindex:


Device UUIDs
============

The UUID of a device (equal to that returned by ``nvidia-smi -L``) is available
in the :attr:`uuid <numba.cuda.cudadrv.driver.Device.uuid>` attribute of a CUDA
device object.

For example, to obtain the UUID of the current device:

.. code-block:: python

   dev = cuda.current_context().device
   # prints e.g. "GPU-e6489c45-5b68-3b03-bab7-0e7c8e809643"
   print(dev.uuid)


.. _cuda-fast-math:

CUDA Fast Math
==============

As noted in :ref:`fast-math`, for certain classes of applications that utilize
floating point, strict IEEE-754 conformance is not required. For this subset of
applications, performance speedups may be possible.

The CUDA target implements :ref:`fast-math` behavior with two differences.

* First, the ``fastmath`` argument to the :func:`@jit decorator
  <numba.cuda.jit>` is limited to the values ``True`` and ``False``.
  When ``True``, the following optimizations are enabled:

  - Flushing of denormals to zero.
  - Use of a fast approximation to the square root function.
  - Use of a fast approximation to the division operation.
  - Contraction of multiply and add operations into single fused multiply-add
    operations.

  See the `documentation for nvvmCompileProgram <https://docs.nvidia.com/cuda/libnvvm-api/group__compilation.html#group__compilation_1g76ac1e23f5d0e2240e78be0e63450346>`_ for more details of these optimizations.

* Secondly, calls to a subset of math module functions on ``float32`` operands
  will be implemented using fast approximate implementations from the libdevice
  library.

  - :func:`math.cos`: Implemented using `__nv_fast_cosf <https://docs.nvidia.com/cuda/libdevice-users-guide/__nv_fast_cosf.html>`_.
  - :func:`math.sin`: Implemented using `__nv_fast_sinf <https://docs.nvidia.com/cuda/libdevice-users-guide/__nv_fast_sinf.html>`_.
  - :func:`math.tan`: Implemented using `__nv_fast_tanf <https://docs.nvidia.com/cuda/libdevice-users-guide/__nv_fast_tanf.html>`_.
  - :func:`math.exp`: Implemented using `__nv_fast_expf <https://docs.nvidia.com/cuda/libdevice-users-guide/__nv_fast_expf.html>`_.
  - :func:`math.log2`: Implemented using `__nv_fast_log2f <https://docs.nvidia.com/cuda/libdevice-users-guide/__nv_fast_log2f.html>`_.
  - :func:`math.log10`: Implemented using `__nv_fast_log10f <https://docs.nvidia.com/cuda/libdevice-users-guide/__nv_fast_log10f.html>`_.
  - :func:`math.log`: Implemented using `__nv_fast_logf <https://docs.nvidia.com/cuda/libdevice-users-guide/__nv_fast_logf.html>`_.
  - :func:`math.pow`: Implemented using `__nv_fast_powf <https://docs.nvidia.com/cuda/libdevice-users-guide/__nv_fast_powf.html>`_.

.. _cuda-random:

Random Number Generation
========================

Numba provides a random number generation algorithm that can be executed on
the GPU.  Due to technical issues with how NVIDIA implemented cuRAND, however,
Numba's GPU random number generator is not based on cuRAND.  Instead, Numba's
GPU RNG is an implementation of the `xoroshiro128+ algorithm
<http://xoroshiro.di.unimi.it/>`_. The xoroshiro128+ algorithm has a period of
``2**128 - 1``, which is shorter than the period of the XORWOW algorithm
used by default in cuRAND, but xoroshiro128+ still passes the BigCrush tests
of random number generator quality.

When using any RNG on the GPU, it is important to make sure that each thread
has its own RNG state, and they have been initialized to produce non-overlapping
sequences.  The  numba.cuda.random module provides a host function to do this,
as well as CUDA device functions to obtain uniformly or normally distributed
random numbers.

.. note:: Numba (like cuRAND) uses the
    `Box-Muller transform <https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform>`
    to generate normally distributed random numbers from a uniform generator.
    However, Box-Muller generates pairs of random numbers, and the current
    implementation only returns one of them.  As a result, generating normally
    distributed values is half the speed of uniformly distributed values.

.. automodule:: numba.cuda.random
    :members: create_xoroshiro128p_states, init_xoroshiro128p_states, xoroshiro128p_uniform_float32, xoroshiro128p_uniform_float64, xoroshiro128p_normal_float32, xoroshiro128p_normal_float64
    :noindex:

A simple example
''''''''''''''''

Here is a sample program that uses the random number generator::

    from __future__ import print_function, absolute_import

    from numba import cuda
    from numba.cuda.random import create_xoroshiro128p_states, xoroshiro128p_uniform_float32
    import numpy as np

    @cuda.jit
    def compute_pi(rng_states, iterations, out):
        """Find the maximum value in values and store in result[0]"""
        thread_id = cuda.grid(1)

        # Compute pi by drawing random (x, y) points and finding what
        # fraction lie inside a unit circle
        inside = 0
        for i in range(iterations):
            x = xoroshiro128p_uniform_float32(rng_states, thread_id)
            y = xoroshiro128p_uniform_float32(rng_states, thread_id)
            if x**2 + y**2 <= 1.0:
                inside += 1

        out[thread_id] = 4.0 * inside / iterations

    threads_per_block = 64
    blocks = 24
    rng_states = create_xoroshiro128p_states(threads_per_block * blocks, seed=1)
    out = np.zeros(threads_per_block * blocks, dtype=np.float32)

    compute_pi[blocks, threads_per_block](rng_states, 10000, out)
    print('pi:', out.mean())

An example of managing RNG state size and using a 3D grid
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The number of RNG states scales with the number of threads using the RNG, so it
is often better to use strided loops in conjunction with the RNG in order to
keep the state size manageable.

In the following example, which initializes a large 3D array with random
numbers, using one thread per output element would result in 453,617,100 RNG
states.  This would take a long time to initialize and poorly utilize the GPU.
Instead, it uses a fixed size 3D grid with a total of 2,097,152 (``(16 ** 3) *
(8 ** 3)``) threads striding over the output array. The 3D thread indices
``startx``, ``starty``, and ``startz``  are linearized into a 1D index,
``tid``, to index into the 2,097,152 RNG states.


.. literalinclude:: ../../../numba/cuda/tests/doc_examples/test_random.py
   :language: python
   :caption: from ``test_ex_3d_grid of ``numba/cuda/tests/doc_example/test_random.py``
   :start-after: magictoken.ex_3d_grid.begin
   :end-before: magictoken.ex_3d_grid.end
   :dedent: 8
   :linenos:

Supported Atomic Operations
===========================

Numba provides access to some of the atomic operations supported in CUDA. Those
that are presently implemented are as follows:

.. automodule:: numba.cuda
    :members: atomic
    :noindex:

Example
'''''''

The following code demonstrates the use of :class:`numba.cuda.atomic.max` to
find the maximum value in an array. Note that this is not the most efficient way
of finding a maximum in this case, but that it serves as an example::

    from numba import cuda
    import numpy as np

    @cuda.jit
    def max_example(result, values):
        """Find the maximum value in values and store in result[0]"""
        tid = cuda.threadIdx.x
        bid = cuda.blockIdx.x
        bdim = cuda.blockDim.x
        i = (bid * bdim) + tid
        cuda.atomic.max(result, 0, values[i])


    arr = np.random.rand(16384)
    result = np.zeros(1, dtype=np.float64)

    max_example[256,64](result, arr)
    print(result[0]) # Found using cuda.atomic.max
    print(max(arr))  # Print max(arr) for comparison (should be equal!)


Multiple dimension arrays are supported by using a tuple of ints for the index::


    @cuda.jit
    def max_example_3d(result, values):
        """
        Find the maximum value in values and store in result[0].
        Both result and values are 3d arrays.
        """
        i, j, k = cuda.grid(3)
        # Atomically store to result[0,1,2] from values[i, j, k]
        cuda.atomic.max(result, (0, 1, 2), values[i, j, k])

    arr = np.random.rand(1000).reshape(10,10,10)
    result = np.zeros((3, 3, 3), dtype=np.float64)
    max_example_3d[(2, 2, 2), (5, 5, 5)](result, arr)
    print(result[0, 1, 2], '==', np.max(arr))



====================
Writing CUDA Kernels
====================

Introduction
============

CUDA has an execution model unlike the traditional sequential model used
for programming CPUs.  In CUDA, the code you write will be executed by
multiple threads at once (often hundreds or thousands).  Your solution will
be modeled by defining a thread hierarchy of *grid*, *blocks* and *threads*.

Numba's CUDA support exposes facilities to declare and manage this
hierarchy of threads.  The facilities are largely similar to those
exposed by NVidia's CUDA C language.

Numba also exposes three kinds of GPU memory: global :ref:`device memory
<cuda-device-memory>` (the large, relatively slow
off-chip memory that's connected to the GPU itself), on-chip
:ref:`shared memory <cuda-shared-memory>` and :ref:`local memory <cuda-local-memory>`.
For all but the simplest algorithms, it is important that you carefully
consider how to use and access memory in order to minimize bandwidth
requirements and contention.


Kernel declaration
==================

A *kernel function* is a GPU function that is meant to be called from CPU
code (*).  It gives it two fundamental characteristics:

* kernels cannot explicitly return a value; all result data must be written
  to an array passed to the function (if computing a scalar, you will
  probably pass a one-element array);

* kernels explicitly declare their thread hierarchy when called: i.e.
  the number of thread blocks and the number of threads per block
  (note that while a kernel is compiled once, it can be called multiple
  times with different block sizes or grid sizes).

At first sight, writing a CUDA kernel with Numba looks very much like
writing a :term:`JIT function` for the CPU::

    @cuda.jit
    def increment_by_one(an_array):
        """
        Increment all array elements by one.
        """
        # code elided here; read further for different implementations

(*) Note: newer CUDA devices support device-side kernel launching; this feature
is called *dynamic parallelism* but Numba does not support it currently)


.. _cuda-kernel-invocation:

Kernel invocation
=================

A kernel is typically launched in the following way::

    threadsperblock = 32
    blockspergrid = (an_array.size + (threadsperblock - 1)) // threadsperblock
    increment_by_one[blockspergrid, threadsperblock](an_array)

We notice two steps here:

* Instantiate the kernel proper, by specifying a number of blocks
  (or "blocks per grid"), and a number of threads per block.  The product
  of the two will give the total number of threads launched.  Kernel
  instantiation is done by taking the compiled kernel function
  (here ``increment_by_one``) and indexing it with a tuple of integers.

* Running the kernel, by passing it the input array (and any separate
  output arrays if necessary). Kernels run asynchronously: launches queue their
  execution on the device and then return immediately.  You can use
  :func:`cuda.synchronize() <numba.cuda.synchronize>` to wait for all previous
  kernel launches to finish executing.

.. note:: Passing an array that resides in host memory will implicitly cause a
   copy back to the host, which will be synchronous. In this case, the kernel
   launch will not return until the data is copied back, and therefore appears
   to execute synchronously.

Choosing the block size
-----------------------

It might seem curious to have a two-level hierarchy when declaring the
number of threads needed by a kernel.  The block size (i.e. number of
threads per block) is often crucial:

* On the software side, the block size determines how many threads
  share a given area of :ref:`shared memory <cuda-shared-memory>`.

* On the hardware side, the block size must be large enough for full
  occupation of execution units; recommendations can be found in the
  `CUDA C Programming Guide`_.

Multi-dimensional blocks and grids
----------------------------------

To help deal with multi-dimensional arrays, CUDA allows you to specify
multi-dimensional blocks and grids.  In the example above, you could
make ``blockspergrid`` and ``threadsperblock`` tuples of one, two
or three integers.  Compared to 1D declarations of equivalent sizes,
this doesn't change anything to the efficiency or behaviour of generated
code, but can help you write your algorithms in a more natural way.


Thread positioning
==================

When running a kernel, the kernel function's code is executed by every
thread once.  It therefore has to know which thread it is in, in order
to know which array element(s) it is responsible for (complex algorithms
may define more complex responsibilities, but the underlying principle
is the same).

One way is for the thread to determine its position in the grid and block
and manually compute the corresponding array position::

    @cuda.jit
    def increment_by_one(an_array):
        # Thread id in a 1D block
        tx = cuda.threadIdx.x
        # Block id in a 1D grid
        ty = cuda.blockIdx.x
        # Block width, i.e. number of threads per block
        bw = cuda.blockDim.x
        # Compute flattened index inside the array
        pos = tx + ty * bw
        if pos < an_array.size:  # Check array boundaries
            an_array[pos] += 1

.. note:: Unless you are sure the block size and grid size is a divisor
   of your array size, you **must** check boundaries as shown above.

:attr:`.threadIdx`, :attr:`.blockIdx`, :attr:`.blockDim` and :attr:`.gridDim`
are special objects provided by the CUDA backend for the sole purpose of
knowing the geometry of the thread hierarchy and the position of the
current thread within that geometry.

These objects can be 1D, 2D or 3D, depending on how the kernel was
:ref:`invoked <cuda-kernel-invocation>`.  To access the value at each
dimension, use the ``x``, ``y`` and ``z`` attributes of these objects,
respectively.

.. attribute:: numba.cuda.threadIdx
   :noindex:

   The thread indices in the current thread block.  For 1D blocks, the index
   (given by the ``x`` attribute) is an integer spanning the range from 0
   inclusive to :attr:`numba.cuda.blockDim` exclusive.  A similar rule
   exists for each dimension when more than one dimension is used.

.. attribute:: numba.cuda.blockDim
   :noindex:

   The shape of the block of threads, as declared when instantiating the
   kernel.  This value is the same for all threads in a given kernel, even
   if they belong to different blocks (i.e. each block is "full").

.. attribute:: numba.cuda.blockIdx
   :noindex:

   The block indices in the grid of threads launched a kernel.  For a 1D grid,
   the index (given by the ``x`` attribute) is an integer spanning the range
   from 0 inclusive to :attr:`numba.cuda.gridDim` exclusive.  A similar rule
   exists for each dimension when more than one dimension is used.

.. attribute:: numba.cuda.gridDim
   :noindex:

   The shape of the grid of blocks, i.e. the total number of blocks launched
   by this kernel invocation, as declared when instantiating the kernel.

Absolute positions
------------------

Simple algorithms will tend to always use thread indices in the
same way as shown in the example above.  Numba provides additional facilities
to automate such calculations:

.. function:: numba.cuda.grid(ndim)
   :noindex:

   Return the absolute position of the current thread in the entire
   grid of blocks.  *ndim* should correspond to the number of dimensions
   declared when instantiating the kernel.  If *ndim* is 1, a single integer
   is returned.  If *ndim* is 2 or 3, a tuple of the given number of
   integers is returned.

.. function:: numba.cuda.gridsize(ndim)
   :noindex:

   Return the absolute size (or shape) in threads of the entire grid of
   blocks.  *ndim* has the same meaning as in :func:`.grid` above.

With these functions, the incrementation example can become::

    @cuda.jit
    def increment_by_one(an_array):
        pos = cuda.grid(1)
        if pos < an_array.size:
            an_array[pos] += 1

The same example for a 2D array and grid of threads would be::

    @cuda.jit
    def increment_a_2D_array(an_array):
        x, y = cuda.grid(2)
        if x < an_array.shape[0] and y < an_array.shape[1]:
           an_array[x, y] += 1

Note the grid computation when instantiating the kernel must still be
done manually, for example::

    threadsperblock = (16, 16)
    blockspergrid_x = math.ceil(an_array.shape[0] / threadsperblock[0])
    blockspergrid_y = math.ceil(an_array.shape[1] / threadsperblock[1])
    blockspergrid = (blockspergrid_x, blockspergrid_y)
    increment_a_2D_array[blockspergrid, threadsperblock](an_array)


Further Reading
----------------

Please refer to the the `CUDA C Programming Guide`_ for a detailed discussion
of CUDA programming.


.. _CUDA C Programming Guide: http://docs.nvidia.com/cuda/cuda-c-programming-guide
========================================
Supported Python features in CUDA Python
========================================

This page lists the Python features supported in the CUDA Python.  This includes
all kernel and device functions compiled with ``@cuda.jit`` and other higher
level Numba decorators that targets the CUDA GPU.

Language
========

Execution Model
---------------

CUDA Python maps directly to the *single-instruction multiple-thread*
execution (SIMT) model of CUDA.  Each instruction is implicitly
executed by multiple threads in parallel.  With this execution model, array
expressions are less useful because we don't want multiple threads to perform
the same task.  Instead, we want threads to perform a task in a cooperative
fashion.

For details please consult the
`CUDA Programming Guide
<http://docs.nvidia.com/cuda/cuda-c-programming-guide/#programming-model>`_.

Constructs
----------

The following Python constructs are not supported:

* Exception handling (``try .. except``, ``try .. finally``)
* Context management (the ``with`` statement)
* Comprehensions (either list, dict, set or generator comprehensions)
* Generator (any ``yield`` statements)

The ``raise`` statement is supported.

The ``assert`` statement is supported, but only has an effect when
``debug=True`` is passed to the :func:`numba.cuda.jit` decorator. This is
similar to the behavior of the ``assert`` keyword in CUDA C/C++, which is
ignored unless compiling with device debug turned on.


Printing of strings, integers, and floats is supported, but printing is an
asynchronous operation - in order to ensure that all output is printed after a
kernel launch, it is necessary to call :func:`numba.cuda.synchronize`. Eliding
the call to ``synchronize`` is acceptable, but output from a kernel may appear
during other later driver operations (e.g. subsequent kernel launches, memory
transfers, etc.), or fail to appear before the program execution completes. Up
to 32 arguments may be passed to the ``print`` function - if more are passed
then a format string will be emitted instead and a warning will be produced.
This is due to a general limitation in CUDA printing, as outlined in the
`section on limitations in printing
<https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#limitations>`_
in the CUDA C++ Programming Guide.

Built-in types
===============

The following built-in types support are inherited from CPU nopython mode.

* int
* float
* complex
* bool
* None
* tuple

See :ref:`nopython built-in types <pysupported-builtin-types>`.

There is also some very limited support for character sequences (bytes and
unicode strings) used in NumPy arrays. Note that this support can only be used
with CUDA 11.2 onwards.

Built-in functions
==================

The following built-in functions are supported:

* :func:`abs`
* :class:`bool`
* :class:`complex`
* :func:`enumerate`
* :class:`float`
* :class:`int`: only the one-argument form
* :func:`len`
* :func:`min`: only the multiple-argument form
* :func:`max`: only the multiple-argument form
* :func:`pow`
* :class:`range`
* :func:`round`
* :func:`zip`


Standard library modules
========================


``cmath``
---------

The following functions from the :mod:`cmath` module are supported:

* :func:`cmath.acos`
* :func:`cmath.acosh`
* :func:`cmath.asin`
* :func:`cmath.asinh`
* :func:`cmath.atan`
* :func:`cmath.atanh`
* :func:`cmath.cos`
* :func:`cmath.cosh`
* :func:`cmath.exp`
* :func:`cmath.isfinite`
* :func:`cmath.isinf`
* :func:`cmath.isnan`
* :func:`cmath.log`
* :func:`cmath.log10`
* :func:`cmath.phase`
* :func:`cmath.polar`
* :func:`cmath.rect`
* :func:`cmath.sin`
* :func:`cmath.sinh`
* :func:`cmath.sqrt`
* :func:`cmath.tan`
* :func:`cmath.tanh`

``math``
--------

The following functions from the :mod:`math` module are supported:

* :func:`math.acos`
* :func:`math.asin`
* :func:`math.atan`
* :func:`math.acosh`
* :func:`math.asinh`
* :func:`math.atanh`
* :func:`math.cos`
* :func:`math.sin`
* :func:`math.tan`
* :func:`math.hypot`
* :func:`math.cosh`
* :func:`math.sinh`
* :func:`math.tanh`
* :func:`math.atan2`
* :func:`math.erf`
* :func:`math.erfc`
* :func:`math.exp`
* :func:`math.expm1`
* :func:`math.fabs`
* :func:`math.frexp`
* :func:`math.ldexp`
* :func:`math.gamma`
* :func:`math.lgamma`
* :func:`math.log`
* :func:`math.log2`
* :func:`math.log10`
* :func:`math.log1p`
* :func:`math.sqrt`
* :func:`math.remainder`: Python 3.7+
* :func:`math.pow`
* :func:`math.ceil`
* :func:`math.floor`
* :func:`math.copysign`
* :func:`math.fmod`
* :func:`math.modf`
* :func:`math.isnan`
* :func:`math.isinf`
* :func:`math.isfinite`


``operator``
------------

The following functions from the :mod:`operator` module are supported:

* :func:`operator.add`
* :func:`operator.and_`
* :func:`operator.eq`
* :func:`operator.floordiv`
* :func:`operator.ge`
* :func:`operator.gt`
* :func:`operator.iadd`
* :func:`operator.iand`
* :func:`operator.ifloordiv`
* :func:`operator.ilshift`
* :func:`operator.imod`
* :func:`operator.imul`
* :func:`operator.invert`
* :func:`operator.ior`
* :func:`operator.ipow`
* :func:`operator.irshift`
* :func:`operator.isub`
* :func:`operator.itruediv`
* :func:`operator.ixor`
* :func:`operator.le`
* :func:`operator.lshift`
* :func:`operator.lt`
* :func:`operator.mod`
* :func:`operator.mul`
* :func:`operator.ne`
* :func:`operator.neg`
* :func:`operator.not_`
* :func:`operator.or_`
* :func:`operator.pos`
* :func:`operator.pow`
* :func:`operator.rshift`
* :func:`operator.sub`
* :func:`operator.truediv`
* :func:`operator.xor`


Numpy support
=============

Due to the CUDA programming model, dynamic memory allocation inside a kernel is
inefficient and is often not needed.  Numba disallows any memory allocating features.
This disables a large number of NumPy APIs.  For best performance, users should write
code such that each thread is dealing with a single element at a time.

Supported numpy features:

* accessing `ndarray` attributes `.shape`, `.strides`, `.ndim`, `.size`, etc..
* scalar ufuncs that have equivalents in the `math` module; i.e. ``np.sin(x[0])``, where x is a 1D array.
* indexing and slicing works.

Unsupported numpy features:

* array creation APIs.
* array methods.
* functions that returns a new array.
CUDA Ufuncs and Generalized Ufuncs
==================================

This page describes the CUDA ufunc-like object.

To support the programming pattern of CUDA programs, CUDA Vectorize and
GUVectorize cannot produce a conventional ufunc.  Instead, a ufunc-like
object is returned.  This object is a close analog but not fully
compatible with a regular NumPy ufunc.  The CUDA ufunc adds support for
passing intra-device arrays (already on the GPU device) to reduce
traffic over the PCI-express bus.  It also accepts a `stream` keyword
for launching in asynchronous mode.

Example: Basic Example
------------------------

::

    import math
    from numba import vectorize, cuda
    import numpy as np

    @vectorize(['float32(float32, float32, float32)',
                'float64(float64, float64, float64)'],
               target='cuda')
    def cu_discriminant(a, b, c):
        return math.sqrt(b ** 2 - 4 * a * c)

    N = 10000
    dtype = np.float32

    # prepare the input
    A = np.array(np.random.sample(N), dtype=dtype)
    B = np.array(np.random.sample(N) + 10, dtype=dtype)
    C = np.array(np.random.sample(N), dtype=dtype)

    D = cu_discriminant(A, B, C)

    print(D)  # print result

Example: Calling Device Functions
----------------------------------

All CUDA ufunc kernels have the ability to call other CUDA device functions::

    from numba import vectorize, cuda

    # define a device function
    @cuda.jit('float32(float32, float32, float32)', device=True, inline=True)
    def cu_device_fn(x, y, z):
        return x ** y / z

    # define a ufunc that calls our device function
    @vectorize(['float32(float32, float32, float32)'], target='cuda')
    def cu_ufunc(x, y, z):
        return cu_device_fn(x, y, z)


Generalized CUDA ufuncs
-----------------------

Generalized ufuncs may be executed on the GPU using CUDA, analogous to
the CUDA ufunc functionality.  This may be accomplished as follows::

    from numba import guvectorize

    @guvectorize(['void(float32[:,:], float32[:,:], float32[:,:])'], 
                 '(m,n),(n,p)->(m,p)', target='cuda')
    def matmulcore(A, B, C):
        ...

There are times when the gufunc kernel uses too many of a GPU's
resources, which can cause the kernel launch to fail.  The user can
explicitly control the maximum size of the thread block by setting
the `max_blocksize` attribute on the compiled gufunc object.

::

    from numba import guvectorize

    @guvectorize(..., target='cuda')
    def very_complex_kernel(A, B, C):
        ...

    very_complex_kernel.max_blocksize = 32  # limits to 32 threads per block

.. comment

    Example: A Chunk at a Time
    ---------------------------

    Partitioning your data into chunks allows computation and memory transfer
    to be overlapped.  This can increase the throughput of your ufunc and
    enables your ufunc to operate on data that is larger than the memory
    capacity of your GPU.  For example:

    ::

        import math
        from numba import vectorize, cuda
        import numpy as np

        # the ufunc kernel
        def discriminant(a, b, c):
            return math.sqrt(b ** 2 - 4 * a * c)

        cu_discriminant = vectorize(['float32(float32, float32, float32)',
                                     'float64(float64, float64, float64)'],
                                    target='cuda')(discriminant)

        N = int(1e+8)
        dtype = np.float32

        # prepare the input
        A = np.array(np.random.sample(N), dtype=dtype)
        B = np.array(np.random.sample(N) + 10, dtype=dtype)
        C = np.array(np.random.sample(N), dtype=dtype)
        D = np.empty(A.shape, dtype=A.dtype)

        # create a CUDA stream
        stream = cuda.stream()

        chunksize = 1e+6
        chunkcount = N // chunksize

        # partition numpy arrays into chunks
        # no copying is performed
        sA = np.split(A, chunkcount)
        sB = np.split(B, chunkcount)
        sC = np.split(C, chunkcount)
        sD = np.split(D, chunkcount)

        device_ptrs = []

        with stream.auto_synchronize():
            # every operation in this context with be launched asynchronously
            # by using the CUDA stream

            # for each chunk
            for a, b, c, d in zip(sA, sB, sC, sD):
                # transfer to device
                dA = cuda.to_device(a, stream)
                dB = cuda.to_device(b, stream)
                dC = cuda.to_device(c, stream)
                dD = cuda.to_device(d, stream, copy=False) # no copying
                # launch kernel
                cu_discriminant(dA, dB, dC, out=dD, stream=stream)
                # retrieve result
                dD.copy_to_host(d, stream)
                # store device pointers to prevent them from freeing before
                # the kernel is scheduled
                device_ptrs.extend([dA, dB, dC, dD])

        # data is ready at this point inside D
==================
Cooperative Groups
==================

Supported features
------------------

Numba's Cooperative Groups support presently provides grid groups and grid
synchronization, along with cooperative kernel launches.

Cooperative groups are supported on Linux, and Windows for devices in `TCC
mode
<https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#tesla-compute-cluster-mode-for-windows>`_.
Cooperative Groups also require the CUDA Device Runtime library, ``cudadevrt``,
to be available - for conda default channel-installed CUDA toolkit packages, it
is only available in versions 10.2 onwards. System-installed toolkits (e.g. from
NVIDIA distribution packages or runfiles) all include ``cudadevrt``.

Using Grid Groups
-----------------

To get the current grid group, use the :meth:`cg.this_grid()
<numba.cuda.cg.this_grid>` function:

.. code-block:: python

   g = cuda.cg.this_grid()

Synchronizing the grid is done with the :meth:`sync()
<numba.cuda.cg.GridGroup.sync>` method of the grid group:

.. code-block:: python

   g.sync()


Cooperative Launches
--------------------

Unlike the CUDA C/C++ API, a cooperative launch is invoked using the same syntax
as a normal kernel launch - Numba automatically determines whether a cooperative
launch is required based on whether a grid group is synchronized in the kernel.

The grid size limit for a cooperative launch is more restrictive than for a
normal launch - the grid must be no larger than the maximum number of active
blocks on the device on which it is launched. To get maximum grid size for a
cooperative launch of a kernel with a given block size and dynamic shared
memory requirement, use the ``max_cooperative_grid_blocks()`` method of kernel
overloads:

.. automethod:: numba.cuda.compiler._Kernel.max_cooperative_grid_blocks

This can be used to ensure that the kernel is launched with no more than the
maximum number of blocks. Exceeding the maximum number of blocks for the
cooperative launch will result in a ``CUDA_ERROR_COOPERATIVE_LAUNCH_TOO_LARGE``
error. 


Applications and Example
------------------------

Grid group synchronization can be used to implement a global barrier across all
threads in the grid - applications of this include a global reduction to a
single value, or looping over rows of a large matrix sequentially using the
entire grid to operate on column elements in parallel.

In the following example, rows are written sequentially by the grid. Each thread
in the grid reads a value from the previous row written by it's *opposite*
thread. A grid sync is needed to ensure that threads in the grid don't run ahead
of threads in other blocks, or fail to see updates from their opposite thread.

First we'll define our kernel:

.. literalinclude:: ../../../numba/cuda/tests/doc_examples/test_cg.py
   :language: python
   :caption: from ``test_grid_sync`` of ``numba/cuda/tests/doc_example/test_cg.py``
   :start-after: magictoken.ex_grid_sync_kernel.begin
   :end-before: magictoken.ex_grid_sync_kernel.end
   :dedent: 8
   :linenos:

Then create some empty input data and determine the grid and block sizes:

.. literalinclude:: ../../../numba/cuda/tests/doc_examples/test_cg.py
   :language: python
   :caption: from ``test_grid_sync`` of ``numba/cuda/tests/doc_example/test_cg.py``
   :start-after: magictoken.ex_grid_sync_data.begin
   :end-before: magictoken.ex_grid_sync_data.end
   :dedent: 8
   :linenos:

Finally we launch the kernel and print the result:

.. literalinclude:: ../../../numba/cuda/tests/doc_examples/test_cg.py
   :language: python
   :caption: from ``test_grid_sync`` of ``numba/cuda/tests/doc_example/test_cg.py``
   :start-after: magictoken.ex_grid_sync_launch.begin
   :end-before: magictoken.ex_grid_sync_launch.end
   :dedent: 8
   :linenos:


The maximum grid size for ``sequential_rows`` can be enquired using:


.. code-block:: python

   overload = sequential_rows.overloads[(int32[:,::1],)
   max_blocks = overload.max_cooperative_grid_blocks(blockdim)
   print(max_blocks)
   # 1152 (e.g. on Quadro RTX 8000 with Numba 0.52.1 and CUDA 11.0)

.. _cudafaq:

=================================================
CUDA Frequently Asked Questions
=================================================

nvprof reports "No kernels were profiled"
-----------------------------------------

When using the ``nvprof`` tool to profile Numba jitted code for the CUDA
target, the output contains ``No kernels were profiled`` but there are clearly
running kernels present, what is going on?

This is quite likely due to the profiling data not being flushed on program
exit, see the `NVIDIA CUDA documentation
<http://docs.nvidia.com/cuda/profiler-users-guide/#flush-profile-data>`_ for
details. To fix this simply add a call to ``numba.cuda.profile_stop()`` prior
to the exit point in your program (or wherever you want to stop profiling).
For more on CUDA profiling support in Numba, see :ref:`cuda-profiling`.

Writing Device Functions
========================

CUDA device functions can only be invoked from within the device (by a kernel
or another device function).  To define a device function::

    from numba import cuda

    @cuda.jit(device=True)
    def a_device_function(a, b):
        return a + b

Unlike a kernel function, a device function can return a value like normal
functions.
========
Overview
========

Numba supports CUDA GPU programming by directly compiling a restricted subset
of Python code into CUDA kernels and device functions following the CUDA
execution model.  Kernels written in Numba appear to have direct access
to NumPy arrays.  NumPy arrays are transferred between the CPU and the
GPU automatically.


Terminology
===========

Several important terms in the topic of CUDA programming are listed here:

- *host*: the CPU
- *device*: the GPU
- *host memory*: the system main memory
- *device memory*: onboard memory on a GPU card
- *kernels*: a GPU function launched by the host and executed on the device
- *device function*: a GPU function executed on the device which can only be
  called from the device (i.e. from a kernel or another device function)


Programming model
=================

Most CUDA programming facilities exposed by Numba map directly to the CUDA
C language offered by NVidia.  Therefore, it is recommended you read the
official `CUDA C programming guide <http://docs.nvidia.com/cuda/cuda-c-programming-guide>`_.


Requirements
============

Supported GPUs
--------------

Numba supports CUDA-enabled GPUs with Compute Capability 3.0 or greater.
Support for devices with Compute Capability less than 5.3 is deprecated, and
will be removed in the next Numba release (0.56).

Devices with Compute Capability 5.3 or greater include (but are not limited to):

- Embedded platforms: NVIDIA Jetson Nano, TX1, TX2, Xavier NX, AGX Xavier.
- Desktop / Server GPUs: All GPUs with Pascal microarchitecture or later. E.g.
  GTX 10 / 16 series, RTX 20 / 30 series, Quadro P / V / RTX series, RTX A series.
- Laptop GPUs: All GPUs with Pascal microarchitecture or later. E.g. MX series,
  Quadro P / T series (mobile), RTX 20 / 30 series (mobile), RTX A series (mobile).

Software
--------

Numba aims to support CUDA Toolkit versions released within the last 3 years.
An NVIDIA driver sufficient for the toolkit version is also required.
Presently:

* 9.2 is the minimum required toolkit version.
* Support for versions less than 10.2 is deprecated, and will be removed in the
  next Numba release (0.56).
* 11.2 or later is recommended, as it uses an NVVM version based on LLVM 7 (as
  opposed to 3.4 in earlier releases).

CUDA is supported on 64-bit Linux and Windows.

If you are using Conda, you can install the CUDA toolkit with::

   $ conda install cudatoolkit

If you are not using Conda or if you want to use a different version of CUDA
toolkit, the following describes how Numba searches for a CUDA toolkit
installation.

.. _cuda-bindings:

CUDA Bindings
~~~~~~~~~~~~~

Numba supports interacting with the CUDA Driver API via the `NVIDIA CUDA Python
bindings <https://nvidia.github.io/cuda-python/>`_ and its own ctypes-based
bindings. Functionality is equivalent between the two bindings. The
ctypes-based bindings are presently the default, but the NVIDIA bindings will
be used by default (if they are available in the environment) in a future Numba
release.

You can install the NVIDIA bindings with::

   $ conda install nvidia::cuda-python

if you are using Conda, or::

   $ pip install cuda-python

if you are using pip.

The use of the NVIDIA bindings is enabled by setting the environment variable
:envvar:`NUMBA_CUDA_USE_NVIDIA_BINDING` to ``"1"``.

.. _cudatoolkit-lookup:

Setting CUDA Installation Path
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Numba searches for a CUDA toolkit installation in the following order:

1. Conda installed `cudatoolkit` package.
2. Environment variable ``CUDA_HOME``, which points to the directory of the
   installed CUDA toolkit (i.e. ``/home/user/cuda-10``)
3. System-wide installation at exactly ``/usr/local/cuda`` on Linux platforms.
   Versioned installation paths (i.e. ``/usr/local/cuda-10.0``) are intentionally
   ignored.  Users can use ``CUDA_HOME`` to select specific versions.

In addition to the CUDA toolkit libraries, which can be installed by conda into
an environment or installed system-wide by the `CUDA SDK installer
<(https://developer.nvidia.com/cuda-downloads)>`_, the CUDA target in Numba
also requires an up-to-date NVIDIA graphics driver.  Updated graphics drivers
are also installed by the CUDA SDK installer, so there is no need to do both.
Note that on macOS, the CUDA SDK must be installed to get the required driver,
and the driver is only supported on macOS prior to 10.14 (Mojave).  If the
``libcuda`` library is in a non-standard location, users can set environment
variable ``NUMBA_CUDA_DRIVER`` to the file path (not the directory path) of the
shared library file.


Missing CUDA Features
=====================

Numba does not implement all features of CUDA, yet.  Some missing features
are listed below:

* dynamic parallelism
* texture memory

.. _cuda-index:

Numba for CUDA GPUs
===================

.. toctree::

   overview.rst
   kernels.rst
   memory.rst
   device-functions.rst
   cudapysupported.rst
   fastmath.rst
   intrinsics.rst
   cooperative_groups.rst
   random.rst
   device-management.rst
   examples.rst
   simulator.rst
   reduction.rst
   ufunc.rst
   ipc.rst
   cuda_array_interface.rst
   external-memory.rst
   bindings.rst
   faq.rst
CUDA Bindings
=============

Numba supports two bindings to the CUDA Driver APIs: its own internal bindings
based on ctypes, and the official `NVIDIA CUDA Python bindings
<https://nvidia.github.io/cuda-python/>`_. Functionality is equivalent between
the two bindings.

The internal bindings are used by default. If the NVIDIA bindings are installed,
then they can be used by setting the environment variable
``NUMBA_CUDA_USE_NVIDIA_BINDING`` to ``1`` prior to the import of Numba. Once
Numba has been imported, the selected binding cannot be changed.


Per-Thread Default Streams
--------------------------

Responsibility for handling Per-Thread Default Streams (PTDS) is delegated to
the NVIDIA bindings when they are in use. To use PTDS with the NVIDIA bindings,
set the environment variable ``CUDA_PYTHON_CUDA_PER_THREAD_DEFAULT_STREAM`` to
``1`` instead of Numba's environmnent variable
:envvar:`NUMBA_CUDA_PER_THREAD_DEFAULT_STREAM`.

.. seealso::

   The `Default Stream section
   <https://nvidia.github.io/cuda-python/release/11.6.0-notes.html#default-stream>`_
   in the NVIDIA Bindings documentation.


Roadmap
-------

In Numba 0.56, the NVIDIA Bindings will be used by default, if they are
installed.

In future versions of Numba:

- The internal bindings will be deprecated.
- The internal bindings will be removed.

At present, no specific release is planned for the deprecation or removal of
the internal bindings.
===================
Sharing CUDA Memory
===================

.. _cuda-ipc-memory:

Sharing between process
=======================

Export device array to another process
--------------------------------------

A device array can be shared with another process in the same machine using
the CUDA IPC API.  To do so, use the ``.get_ipc_handle()`` method on the device
array to get a ``IpcArrayHandle`` object, which can be transferred to another
process.


.. automethod:: numba.cuda.cudadrv.devicearray.DeviceNDArray.get_ipc_handle
    :noindex:

.. autoclass:: numba.cuda.cudadrv.devicearray.IpcArrayHandle
    :members: open, close


Import IPC memory from another process
--------------------------------------

The following function is used to open IPC handle from another process
as a device array.

.. automethod:: numba.cuda.open_ipc_array
.. _cuda-emm-plugin:

=================================================
External Memory Management (EMM) Plugin interface
=================================================

The :ref:`CUDA Array Interface <cuda-array-interface>` enables sharing of data
between different Python libraries that access CUDA devices. However, each
library manages its own memory distinctly from the others. For example:

- By default, Numba allocates memory on CUDA devices by interacting with the
  CUDA driver API to call functions such as ``cuMemAlloc`` and ``cuMemFree``,
  which is suitable for many use cases.
- The RAPIDS libraries (cuDF, cuML, etc.) use the `RAPIDS Memory Manager (RMM)
  <https://github.com/rapidsai/rmm>`_ for allocating device memory.
- `CuPy <https://cupy.chainer.org/>`_ includes a `memory pool implementation
  <https://docs-cupy.chainer.org/en/stable/reference/memory.html>`_ for both
  device and pinned memory.

When multiple CUDA-aware libraries are used together, it may be preferable for
Numba to defer to another library for memory management. The EMM Plugin
interface facilitates this, by enabling Numba to use another CUDA-aware library
for all allocations and deallocations.

An EMM Plugin is used to facilitate the use of an external library for memory
management. An EMM Plugin can be a part of an external library, or could be
implemented as a separate library.


Overview of External Memory Management
======================================

When an EMM Plugin is in use (see :ref:`setting-emm-plugin`), Numba will make
memory allocations and deallocations through the Plugin. It will never directly call
functions such as ``cuMemAlloc``, ``cuMemFree``, etc.

EMM Plugins always take responsibility for the management of device memory.
However, not all CUDA-aware libraries also support managing host memory, so a
facility for Numba to continue the management of host memory whilst ceding
control of device memory to the EMM is provided (see
:ref:`host-only-cuda-memory-manager`).


Effects on Deallocation Strategies
----------------------------------

Numba's internal :ref:`deallocation-behavior` is designed to increase efficiency
by deferring deallocations until a significant quantity are pending. It also
provides a mechanism for preventing deallocations entirely during critical
sections, using the :func:`~numba.cuda.defer_cleanup` context manager.

When an EMM Plugin is in use, the deallocation strategy is implemented by the
EMM, and Numba's internal deallocation mechanism is not used. The EMM
Plugin could implement:
  
- A similar strategy to the Numba deallocation behaviour, or
- Something more appropriate to the plugin - for example, deallocated memory
  might immediately be returned to a memory pool.

The ``defer_cleanup`` context manager may behave differently with an EMM Plugin
- an EMM Plugin should be accompanied by documentation of the behaviour of the
``defer_cleanup`` context manager when it is in use. For example, a pool
allocator could always immediately return memory to a pool even when the
context manager is in use, but could choose not to free empty pools until
``defer_cleanup`` is not in use.


Management of other objects
---------------------------

In addition to memory, Numba manages the allocation and deallocation of
:ref:`events <events>`, :ref:`streams <streams>`, and modules (a module is a
compiled object, which is generated from ``@cuda.jit``\ -ted functions). The
management of events, streams, and modules is unchanged by the use of an EMM
Plugin.


Asynchronous allocation and deallocation
----------------------------------------

The present EMM Plugin interface does not provide support for asynchronous
allocation and deallocation. This may be added to a future version of the
interface.


Implementing an EMM Plugin
==========================

An EMM Plugin is implemented by deriving from
:class:`~numba.cuda.BaseCUDAMemoryManager`. A summary of considerations for the
implementation follows:

- Numba instantiates one instance of the EMM Plugin class per context. The
  context that owns an EMM Plugin object is accessible through ``self.context``,
  if required.
- The EMM Plugin is transparent to any code that uses Numba - all its methods
  are invoked by Numba, and never need to be called by code that uses Numba.
- The allocation methods ``memalloc``, ``memhostalloc``, and ``mempin``, should
  use the underlying library to allocate and/or pin device or host memory, and
  construct an instance of a :ref:`memory pointer <memory-pointers>`
  representing the memory to return back to Numba. These methods are always
  called when the current CUDA context is the context that owns the EMM Plugin
  instance.
- The ``initialize`` method is called by Numba prior to the first use of the EMM
  Plugin object for a context. This method should do anything required to
  prepare the underlying library for allocations in the current context. This
  method may be called multiple times, and must not invalidate previous state
  when it is called.
- The ``reset`` method is called when all allocations in the context are to be
  cleaned up. It may be called even prior to ``initialize``, and an EMM Plugin
  implementation needs to guard against this.
- To support inter-GPU communication, the ``get_ipc_handle`` method should
  provide an :class:`~numba.cuda.IpcHandle` for a given
  :class:`~numba.cuda.MemoryPointer` instance. This method is part of the EMM
  interface (rather than being handled within Numba) because the base address of
  the allocation is only known by the underlying library. Closing an IPC handle
  is handled internally within Numba.
- It is optional to provide memory info from the ``get_memory_info`` method, which
  provides a count of the total and free memory on the device for the context.
  It is preferrable to implement the method, but this may not be practical for
  all allocators. If memory info is not provided, this method should raise a
  :class:`RuntimeError`.
- The ``defer_cleanup`` method should return a context manager that ensures that
  expensive cleanup operations are avoided whilst it is active. The nuances of
  this will vary between plugins, so the plugin documentation should include an
  explanation of how deferring cleanup affects deallocations, and performance in
  general.
- The ``interface_version`` property is used to ensure that the plugin version
  matches the interface provided by the version of Numba. At present, this
  should always be 1.

Full documentation for the base class follows:

.. autoclass:: numba.cuda.BaseCUDAMemoryManager
   :members: memalloc, memhostalloc, mempin, initialize, get_ipc_handle,
             get_memory_info, reset, defer_cleanup, interface_version
   :member-order: bysource


.. _host-only-cuda-memory-manager:

The Host-Only CUDA Memory Manager
---------------------------------

Some external memory managers will support management of on-device memory but
not host memory. For implementing EMM Plugins using one of these memory
managers, a partial implementation of a plugin that implements host-side
allocation and pinning is provided. To use it, derive from
:class:`~numba.cuda.HostOnlyCUDAMemoryManager` instead of
:class:`~numba.cuda.BaseCUDAMemoryManager`. Guidelines for using this class
are:

- The host-only memory manager implements ``memhostalloc`` and ``mempin`` - the
  EMM Plugin should still implement ``memalloc``.
- If ``reset`` is overridden, it must also call ``super().reset()`` to allow the
  host allocations to be cleaned up.
- If ``defer_cleanup`` is overridden, it must hold an active context manager
  from ``super().defer_cleanup()`` to ensure that host-side cleanup is also
  deferred.

Documentation for the methods of :class:`~numba.cuda.HostOnlyCUDAMemoryManager`
follows:

.. autoclass:: numba.cuda.HostOnlyCUDAMemoryManager
   :members: memhostalloc, mempin, reset, defer_cleanup
   :member-order: bysource


The IPC Handle Mixin
--------------------

An implementation of the ``get_ipc_handle()`` function is is provided in the
``GetIpcHandleMixin`` class. This uses the driver API to determine the base
address of an allocation for opening an IPC handle. If this implementation is
appropriate for an EMM plugin, it can be added by mixing in the
``GetIpcHandleMixin`` class:

.. autoclass:: numba.cuda.GetIpcHandleMixin
   :members: get_ipc_handle


Classes and structures of returned objects
==========================================

This section provides an overview of the classes and structures that need to be
constructed by an EMM Plugin.

.. _memory-pointers:

Memory Pointers
---------------

EMM Plugins should construct memory pointer instances that represent their
allocations, for return to Numba. The appropriate memory pointer class to use in
each method is:

- :class:`~numba.cuda.MemoryPointer`: returned from ``memalloc``
- :class:`~numba.cuda.MappedMemory`: returned from ``memhostalloc`` or
  ``mempin`` when the host memory is mapped into the device memory space.
- :class:`~numba.cuda.PinnedMemory`: return from ``memhostalloc`` or ``mempin``
  when the host memory is not mapped into the device memory space.

Memory pointers can take a finalizer, which is a function that is called when
the buffer is no longer needed. Usually the finalizer will make a call to the
memory management library (either internal to Numba, or external if allocated
by an EMM Plugin) to inform it that the memory is no longer required, and that
it could potentially be freed and/or unpinned. The memory manager may choose to
defer actually cleaning up the memory to any later time after the finalizer
runs - it is not required to free the buffer immediately.

Documentation for the memory pointer classes follows.

.. autoclass:: numba.cuda.MemoryPointer

The ``AutoFreePointer`` class need not be used directly, but is documented here
as it is subclassed by :class:`numba.cuda.MappedMemory`:

.. autoclass:: numba.cuda.cudadrv.driver.AutoFreePointer

.. autoclass:: numba.cuda.MappedMemory

.. autoclass:: numba.cuda.PinnedMemory


Memory Info
-----------

If an implementation of
:meth:`~numba.cuda.BaseCUDAMemoryManager.get_memory_info` is to provide a
result, then it should return an instance of the ``MemoryInfo`` named tuple:

.. autoclass:: numba.cuda.MemoryInfo


IPC
---

An instance of ``IpcHandle`` is required to be returned from an implementation
of :meth:`~numba.cuda.BaseCUDAMemoryManager.get_ipc_handle`:

.. autoclass:: numba.cuda.IpcHandle

Guidance for constructing an IPC handle in the context of implementing an EMM
Plugin:

- The ``memory`` parameter passed to the ``get_ipc_handle`` method of an EMM
  Plugin can be passed as the ``base`` parameter.
- A suitable type for the ``handle`` can be constructed as ``ctypes.c_byte *
  64``. The data for ``handle`` must be populated using a method for obtaining a
  CUDA IPC handle appropriate to the underlying library.
- ``size`` should match the size of the original allocation, which can be
  obtained with ``memory.size`` in ``get_ipc_handle``.
- An appropriate value for ``source_info`` can be created by calling
  ``self.context.device.get_device_identity()``.
- If the underlying memory does not point to the base of an allocation returned
  by the CUDA driver or runtime API (e.g. if a pool allocator is in use) then
  the ``offset`` from the base must be provided.


.. _setting-emm-plugin:

Setting the EMM Plugin
======================

By default, Numba uses its internal memory management - if an EMM Plugin is to
be used, it must be configured. There are two mechanisms for configuring the use
of an EMM Plugin: an environment variable, and a function.


Environment variable
--------------------

A module name can be provided in the environment variable,
``NUMBA_CUDA_MEMORY_MANAGER``. If this environment variable is set, Numba will
attempt to import the module, and and use its ``_numba_memory_manager`` global
variable as the memory manager class. This is primarily useful for running the
Numba test suite with an EMM Plugin, e.g.:

.. code::

   $ NUMBA_CUDA_MEMORY_MANAGER=rmm python -m numba.runtests numba.cuda.tests


Function
--------

The :func:`~numba.cuda.set_memory_manager` function can be used to set the
memory manager at runtime. This should be called prior to the initialization of
any contexts, as EMM Plugin instances are instantiated along with contexts.

.. autofunction:: numba.cuda.set_memory_manager


Resetting the memory manager
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is recommended that the memory manager is set once prior to using any CUDA
functionality, and left unchanged for the remainder of execution. It is possible
to set the memory manager multiple times, noting the following:

* At the time of their creation, contexts are bound to an instance of a memory
  manager for their lifetime.
* Changing the memory manager will have no effect on existing contexts - only
  contexts created after the memory manager was updated will use instances of
  the new memory manager.
* :func:`numba.cuda.close` can be used to destroy contexts after setting the
  memory manager so that they get re-created with the new memory manager.

  - This will invalidate any arrays, streams, events, and modules owned by the
    context.
  - Attempting to use invalid arrays, streams, or events will likely fail with
    an exception being raised due to a ``CUDA_ERROR_INVALID_CONTEXT`` or
    ``CUDA_ERROR_CONTEXT_IS_DESTROYED`` return code from a Driver API function.
  - Attempting to use an invalid module will result in similar, or in some
    cases a segmentation fault / access violation.

.. note:: The invalidation of modules means that all functions compiled with
          ``@cuda.jit`` prior to context destruction will need to be
          redefined, as the code underlying them will also have been unloaded
          from the GPU.

.. _simulator:

=================================================
Debugging CUDA Python with the the CUDA Simulator
=================================================

Numba includes a CUDA Simulator that implements most of the semantics in CUDA
Python using the Python interpreter and some additional Python code. This can
be used to debug CUDA Python code, either by adding print statements to your
code, or by using the debugger to step through the execution of an individual
thread.

The simulator deliberately allows running non-CUDA code like starting a debugger 
and printing arbitrary expressions for debugging purposes. Therefore, it is
best to start from code that compiles for the CUDA target, and then move over to
the simulator to investigate issues.

Execution of kernels is performed by the simulator one block at a time. One
thread is spawned for each thread in the block, and scheduling of the execution
of these threads is left up to the operating system.

Using the simulator
===================

The simulator is enabled by setting the environment variable
:envvar:`NUMBA_ENABLE_CUDASIM` to 1 prior to importing Numba. CUDA Python code 
may then be executed as normal. The easiest way to use the debugger inside a
kernel is to only stop a single thread, otherwise the interaction with the
debugger is difficult to handle. For example, the kernel below will stop in
the thread ``<<<(3,0,0), (1, 0, 0)>>>``::

    @cuda.jit
    def vec_add(A, B, out):
        x = cuda.threadIdx.x
        bx = cuda.blockIdx.x
        bdx = cuda.blockDim.x
        if x == 1 and bx == 3:
            from pdb import set_trace; set_trace()
        i = bx * bdx + x
        out[i] = A[i] + B[i]

when invoked with a one-dimensional grid and one-dimensional blocks.

Supported features
==================

The simulator aims to provide as complete a simulation of execution on a real
GPU as possible - in particular, the following are supported:

* Atomic operations
* Constant memory
* Local memory
* Shared memory: declarations of shared memory arrays must be on separate source
  lines, since the simulator uses source line information to keep track of
  allocations of shared memory across threads.
* Mapped arrays.
* Host and device memory operations: copying and setting memory.
* :func:`.syncthreads` is supported - however, in the case where divergent
  threads enter different :func:`.syncthreads` calls, the launch will not fail,
  but unexpected behaviour will occur. A future version of the simulator may
  detect this condition.
* The stream API is supported, but all operations occur sequentially and
  synchronously, unlike on a real device. Synchronising on a stream is therefore
  a no-op.
* The event API is also supported, but provides no meaningful timing
  information.
* Data transfer to and from the GPU - in particular, creating array objects with
  :func:`.device_array` and :func:`.device_array_like`. The APIs for pinned memory
  :func:`.pinned` and :func:`.pinned_array` are also supported, but no pinning
  takes place.
* The driver API implementation of the list of GPU contexts (``cuda.gpus`` and
  ``cuda.cudadrv.devices.gpus``) is supported, and reports a single GPU context.
  This context can be closed and reset as the real one would.
* The :func:`.detect` function is supported, and reports one device called
  `SIMULATOR`.
* Cooperative grids: A cooperative kernel can be launched, but with only one
  block - the simulator always returns ``1`` from a kernel overload's
  :meth:`~numba.cuda.compiler._Kernel.max_cooperative_grid_blocks` method.

Some limitations of the simulator include:

* It does not perform type checking/type inference. If any argument types to a
  jitted function are incorrect, or if the specification of the type of any
  local variables are incorrect, this will not be detected by the simulator.
* Only one GPU is simulated.
* Multithreaded accesses to a single GPU are not supported, and will result in
  unexpected behaviour.
* Most of the driver API is unimplemented.
* It is not possible to link PTX code with CUDA Python functions.
* Warps and warp-level operations are not yet implemented.
* Because the simulator executes kernels using the Python interpreter,
  structured array access by attribute that works with the hardware target may
  fail in the simulator - see :ref:`structured-array-access`.
* Operations directly against device arrays are only partially supported, that
  is, testing equality, less than, greater than, and basic mathematical 
  operations are supported, but many other operations, such as the in-place 
  operators and bit operators are not.
* The :func:`ffs() <numba.cuda.ffs>` function only works correctly for values
  that can be represented using 32-bit integers.

Obviously, the speed of the simulator is also much lower than that of a real
device. It may be necessary to reduce the size of input data and the size of the
CUDA grid in order to make debugging with the simulator tractable.

Installation
============

Compatibility
-------------

Numba is compatible with Python 3.7--3.10, and Numpy versions 1.18 or later.

Our supported platforms are:

* Linux x86 (32-bit and 64-bit)
* Linux ppcle64 (POWER8, POWER9)
* Windows 7 and later (32-bit and 64-bit)
* OS X 10.9 and later (64-bit and unofficial support on M1/Arm64)
* \*BSD (unofficial support only)
* NVIDIA GPUs of compute capability 5.3 and later

  * Compute capabilities 3.0 - 5.2 are supported, but deprecated.
* ARMv7 (32-bit little-endian, such as Raspberry Pi 2 and 3)
* ARMv8 (64-bit little-endian, such as the NVIDIA Jetson)

:ref:`numba-parallel` is only available on 64-bit platforms.

Installing using conda on x86/x86_64/POWER Platforms
----------------------------------------------------

The easiest way to install Numba and get updates is by using ``conda``,
a cross-platform package manager and software distribution maintained
by Anaconda, Inc.  You can either use `Anaconda
<https://www.anaconda.com/download>`_ to get the full stack in one download,
or `Miniconda <https://conda.io/miniconda.html>`_ which will install
the minimum packages required for a conda environment.

Once you have conda installed, just type::

    $ conda install numba

or::

    $ conda update numba

Note that Numba, like Anaconda, only supports PPC in 64-bit little-endian mode.

To enable CUDA GPU support for Numba, install the latest `graphics drivers from
NVIDIA <https://www.nvidia.com/Download/index.aspx>`_ for your platform.
(Note that the open source Nouveau drivers shipped by default with many Linux
distributions do not support CUDA.)  Then install the ``cudatoolkit`` package::

    $ conda install cudatoolkit

You do not need to install the CUDA SDK from NVIDIA.


Installing using pip on x86/x86_64 Platforms
--------------------------------------------

Binary wheels for Windows, Mac, and Linux are also available from `PyPI
<https://pypi.org/project/numba/>`_.  You can install Numba using ``pip``::

    $ pip install numba

This will download all of the needed dependencies as well.  You do not need to
have LLVM installed to use Numba (in fact, Numba will ignore all LLVM
versions installed on the system) as the required components are bundled into
the llvmlite wheel.

To use CUDA with Numba installed by `pip`, you need to install the `CUDA SDK
<https://developer.nvidia.com/cuda-downloads>`_ from NVIDIA.  Please refer to
:ref:`cudatoolkit-lookup` for details. Numba can also detect CUDA libraries
installed system-wide on Linux.


.. _numba-install-armv7:

Installing on Linux ARMv7 Platforms
-----------------------------------

`Berryconda <https://github.com/jjhelmus/berryconda>`_ is a
conda-based Python distribution for the Raspberry Pi.  We are now uploading
packages to the ``numba`` channel on Anaconda Cloud for 32-bit little-endian,
ARMv7-based boards, which currently includes the Raspberry Pi 2 and 3,
but not the Pi 1 or Zero.  These can be installed using conda from the
``numba`` channel::

    $ conda install -c numba numba

Berryconda and Numba may work on other Linux-based ARMv7 systems, but this has
not been tested.


Installing on Linux ARMv8 (AArch64) Platforms
---------------------------------------------

We build and test conda packages on the `NVIDIA Jetson TX2
<https://www.nvidia.com/en-us/autonomous-machines/embedded-systems-dev-kits-modules/>`_,
but they are likely to work for other AArch64 platforms.  (Note that while the
Raspberry Pi CPU is 64-bit, Raspbian runs it in 32-bit mode, so look at
:ref:`numba-install-armv7` instead.)

Conda-forge support for AArch64 is still quite experimental and packages are limited,
but it does work enough for Numba to build and pass tests.  To set up the environment:

* Install `miniforge <https://github.com/conda-forge/miniforge>`_.
  This will create a minimal conda environment.

* Then you can install Numba from the ``numba`` channel::

    $ conda install -c numba numba

On CUDA-enabled systems, like the Jetson, the CUDA toolkit should be
automatically detected in the environment.

.. _numba-source-install-instructions:

Installing from source
----------------------

Installing Numba from source is fairly straightforward (similar to other
Python packages), but installing `llvmlite
<https://github.com/numba/llvmlite>`_ can be quite challenging due to the need
for a special LLVM build.  If you are building from source for the purposes of
Numba development, see :ref:`buildenv` for details on how to create a Numba
development environment with conda.

If you are building Numba from source for other reasons, first follow the
`llvmlite installation guide <https://llvmlite.readthedocs.io/en/latest/admin-guide/install.html>`_.
Once that is completed, you can download the latest Numba source code from
`Github <https://github.com/numba/numba>`_::

    $ git clone git://github.com/numba/numba.git

Source archives of the latest release can also be found on
`PyPI <https://pypi.org/project/numba/>`_.  In addition to ``llvmlite``, you will also need:

* A C compiler compatible with your Python installation.  If you are using
  Anaconda, you can use the following conda packages:

  * Linux ``x86``: ``gcc_linux-32`` and ``gxx_linux-32``
  * Linux ``x86_64``: ``gcc_linux-64`` and ``gxx_linux-64``
  * Linux ``POWER``: ``gcc_linux-ppc64le`` and ``gxx_linux-ppc64le``
  * Linux ``ARM``: no conda packages, use the system compiler
  * Mac OSX: ``clang_osx-64`` and ``clangxx_osx-64`` or the system compiler at
    ``/usr/bin/clang`` (Mojave onwards)
  * Windows: a version of Visual Studio appropriate for the Python version in
    use

* `NumPy <http://www.numpy.org/>`_

Then you can build and install Numba from the top level of the source tree::

    $ python setup.py install

.. _numba-source-install-env_vars:

Build time environment variables and configuration of optional components
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below are environment variables that are applicable to altering how Numba would
otherwise build by default along with information on configuration options.

.. envvar:: NUMBA_DISABLE_OPENMP (default: not set)

  To disable compilation of the OpenMP threading backend set this environment
  variable to a non-empty string when building. If not set (default):

  * For Linux and Windows it is necessary to provide OpenMP C headers and
    runtime  libraries compatible with the compiler tool chain mentioned above,
    and for these to be accessible to the compiler via standard flags.
  * For OSX the conda packages ``llvm-openmp`` and ``intel-openmp`` provide
    suitable C headers and libraries. If the compilation requirements are not
    met the OpenMP threading backend will not be compiled

.. envvar:: NUMBA_DISABLE_TBB (default: not set)

  To disable the compilation of the TBB threading backend set this environment
  variable to a non-empty string when building. If not set (default) the TBB C
  headers and libraries must be available at compile time. If building with
  ``conda build`` this requirement can be met by installing the ``tbb-devel``
  package. If not building with ``conda build`` the requirement can be met via a
  system installation of TBB or through the use of the ``TBBROOT`` environment
  variable to provide the location of the TBB installation. For more
  information about setting ``TBBROOT`` see the `Intel documentation <https://software.intel.com/content/www/us/en/develop/documentation/advisor-user-guide/top/appendix/adding-parallelism-to-your-program/adding-the-parallel-framework-to-your-build-environment/defining-the-tbbroot-environment-variable.html>`_.

.. _numba-source-install-check:

Dependency List
---------------

Numba has numerous required and optional dependencies which additionally may
vary with target operating system and hardware. The following lists them all
(as of July 2020).

* Required build time:

  * ``setuptools``
  * ``numpy``
  * ``llvmlite``
  * Compiler toolchain mentioned above

* Required run time:

  * ``setuptools``
  * ``numpy``
  * ``llvmlite``

* Optional build time:

  See :ref:`numba-source-install-env_vars` for more details about additional
  options for the configuration and specification of these optional components.

  * ``llvm-openmp`` (OSX) - provides headers for compiling OpenMP support into
    Numba's threading backend
  * ``intel-openmp`` (OSX) - provides OpenMP library support for Numba's
    threading backend.
  * ``tbb-devel`` - provides TBB headers/libraries for compiling TBB support
    into Numba's threading backend (version >= 2021 required).

* Optional runtime are:

  * ``scipy`` - provides cython bindings used in Numba's ``np.linalg.*``
    support
  * ``tbb`` - provides the TBB runtime libraries used by Numba's TBB threading
    backend (version >= 2021 required).
  * ``jinja2`` - for "pretty" type annotation output (HTML) via the ``numba``
    CLI
  * ``cffi`` - permits use of CFFI bindings in Numba compiled functions
  * ``intel-openmp`` - (OSX) provides OpenMP library support for Numba's OpenMP
    threading backend
  * ``ipython`` - if in use, caching will use IPython's cache
    directories/caching still works
  * ``pyyaml`` - permits the use of a ``.numba_config.yaml``
    file for storing per project configuration options
  * ``colorama`` - makes error message highlighting work
  * ``icc_rt`` - (numba channel) allows Numba to use Intel SVML for extra
    performance
  * ``pygments`` - for "pretty" type annotation
  * ``gdb`` as an executable on the ``$PATH`` - if you would like to use the gdb
    support
  * Compiler toolchain mentioned above, if you would like to use ``pycc`` for
    Ahead-of-Time (AOT) compilation
  * ``r2pipe`` - required for assembly CFG inspection.
  * ``radare2`` as an executable on the ``$PATH`` - required for assembly CFG
    inspection. `See here <https://github.com/radareorg/radare2>`_ for
    information on obtaining and installing.
  * ``graphviz`` - for some CFG inspection functionality.
  * ``pickle5`` - provides Python 3.8 pickling features for faster pickling in
    Python 3.7.
  * ``typeguard`` - used by ``runtests.py`` for
    :ref:`runtime type-checking <type_anno_check>`.
  * ``cuda-python`` - The NVIDIA CUDA Python bindings. See :ref:`cuda-bindings`.
    Numba requires Version 11.6 or greater.

* To build the documentation:

  * ``sphinx``
  * ``pygments``
  * ``sphinx_rtd_theme``
  * ``numpydoc``
  * ``make`` as an executable on the ``$PATH``

Checking your installation
--------------------------

You should be able to import Numba from the Python prompt::

    $ python
    Python 3.10.2 | packaged by conda-forge | (main, Jan 14 2022, 08:02:09) [GCC 9.4.0] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import numba
    >>> numba.__version__
    '0.55.1'

You can also try executing the ``numba --sysinfo`` (or ``numba -s`` for short)
command to report information about your system capabilities. See :ref:`cli` for
further information.

::

    $ numba -s
    System info:
    --------------------------------------------------------------------------------
    __Time Stamp__
    Report started (local time)                   : 2022-01-18 10:35:08.981319

    __Hardware Information__
    Machine                                       : x86_64
    CPU Name                                      : skylake-avx512
    CPU Count                                     : 12
    CPU Features                                  :
    64bit adx aes avx avx2 avx512bw avx512cd avx512dq avx512f avx512vl bmi bmi2
    clflushopt clwb cmov cx16 cx8 f16c fma fsgsbase fxsr invpcid lzcnt mmx
    movbe pclmul pku popcnt prfchw rdrnd rdseed rtm sahf sse sse2 sse3 sse4.1
    sse4.2 ssse3 xsave xsavec xsaveopt xsaves

    __OS Information__
    Platform Name                                 : Linux-5.4.0-94-generic-x86_64-with-glibc2.31
    Platform Release                              : 5.4.0-94-generic
    OS Name                                       : Linux
    OS Version                                    : #106-Ubuntu SMP Thu Jan 6 23:58:14 UTC 2022

    __Python Information__
    Python Compiler                               : GCC 9.4.0
    Python Implementation                         : CPython
    Python Version                                : 3.10.2
    Python Locale                                 : en_GB.UTF-8

    __LLVM information__
    LLVM Version                                  : 11.1.0

    __CUDA Information__
    Found 1 CUDA devices
    id 0      b'Quadro RTX 8000'                              [SUPPORTED]
                          Compute Capability: 7.5
                               PCI Device ID: 0
                                  PCI Bus ID: 21
                                        UUID: GPU-e6489c45-5b68-3b03-bab7-0e7c8e809643
                                    Watchdog: Enabled
                 FP32/FP64 Performance Ratio: 32

(output truncated due to length)
.. Copyright (c) 2017 Intel Corporation
   SPDX-License-Identifier: BSD-2-Clause

.. _numba-stencil:

================================
Using the ``@stencil`` decorator
================================

Stencils are a common computational pattern in which array elements
are updated according to some fixed pattern called the stencil kernel.
Numba provides the ``@stencil`` decorator so that users may
easily specify a stencil kernel and Numba then generates the looping
code necessary to apply that kernel to some input array.  Thus, the
stencil decorator allows clearer, more concise code and in conjunction
with :ref:`the parallel jit option <parallel_jit_option>` enables higher
performance through parallelization of the stencil execution.


Basic usage
===========

An example use of the ``@stencil`` decorator::

   from numba import stencil

   @stencil
   def kernel1(a):
       return 0.25 * (a[0, 1] + a[1, 0] + a[0, -1] + a[-1, 0])

The stencil kernel is specified by what looks like a standard Python
function definition but there are different semantics with
respect to array indexing.
Stencils produce an output array of the same size and shape as the
input array although depending on the kernel definition may have a
different type.
Conceptually, the stencil kernel is run once for each element in the
output array.  The return value from the stencil kernel is the value
written into the output array for that particular element.

The parameter ``a`` represents the input array over which the
kernel is applied.
Indexing into this array takes place with respect to the current element
of the output array being processed.  For example, if element ``(x, y)``
is being processed then ``a[0, 0]`` in the stencil kernel corresponds to
``a[x + 0, y + 0]`` in the input array.  Similarly, ``a[-1, 1]`` in the stencil
kernel corresponds to ``a[x - 1, y + 1]`` in the input array.

Depending on the specified kernel, the kernel may not be applicable to the
borders of the output array as this may cause the input array to be
accessed out-of-bounds.  The way in which the stencil decorator handles
this situation is dependent upon which :ref:`stencil-mode` is selected.
The default mode is for the stencil decorator to set the border elements
of the output array to zero.

To invoke a stencil on an input array, call the stencil as if it were
a regular function and pass the input array as the argument. For example, using
the kernel defined above::

   >>> import numpy as np
   >>> input_arr = np.arange(100).reshape((10, 10))
   array([[ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9],
          [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
          [20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
          [30, 31, 32, 33, 34, 35, 36, 37, 38, 39],
          [40, 41, 42, 43, 44, 45, 46, 47, 48, 49],
          [50, 51, 52, 53, 54, 55, 56, 57, 58, 59],
          [60, 61, 62, 63, 64, 65, 66, 67, 68, 69],
          [70, 71, 72, 73, 74, 75, 76, 77, 78, 79],
          [80, 81, 82, 83, 84, 85, 86, 87, 88, 89],
          [90, 91, 92, 93, 94, 95, 96, 97, 98, 99]])
   >>> output_arr = kernel1(input_arr)
   array([[  0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.],
          [  0.,  11.,  12.,  13.,  14.,  15.,  16.,  17.,  18.,   0.],
          [  0.,  21.,  22.,  23.,  24.,  25.,  26.,  27.,  28.,   0.],
          [  0.,  31.,  32.,  33.,  34.,  35.,  36.,  37.,  38.,   0.],
          [  0.,  41.,  42.,  43.,  44.,  45.,  46.,  47.,  48.,   0.],
          [  0.,  51.,  52.,  53.,  54.,  55.,  56.,  57.,  58.,   0.],
          [  0.,  61.,  62.,  63.,  64.,  65.,  66.,  67.,  68.,   0.],
          [  0.,  71.,  72.,  73.,  74.,  75.,  76.,  77.,  78.,   0.],
          [  0.,  81.,  82.,  83.,  84.,  85.,  86.,  87.,  88.,   0.],
          [  0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.]])
   >>> input_arr.dtype
   dtype('int64')
   >>> output_arr.dtype
   dtype('float64')

Note that the stencil decorator has determined that the output type
of the specified stencil kernel is ``float64`` and has thus created the
output array as ``float64`` while the input array is of type ``int64``.

Stencil Parameters
==================

Stencil kernel definitions may take any number of arguments with
the following provisions.  The first argument must be an array.
The size and shape of the output array will be the same as that of the
first argument.  Additional arguments may either be scalars or
arrays.  For array arguments, those arrays must be at least as large
as the first argument (array) in each dimension.  Array indexing is relative for
all such input array arguments.

.. _stencil-kernel-shape-inference:

Kernel shape inference and border handling
==========================================

In the above example and in most cases, the array indexing in the
stencil kernel will exclusively use ``Integer`` literals.
In such cases, the stencil decorator is able to analyze the stencil
kernel to determine its size.  In the above example, the stencil
decorator determines that the kernel is ``3 x 3`` in shape since indices
``-1`` to ``1`` are used for both the first and second dimensions.  Note that
the stencil decorator also correctly handles non-symmetric and
non-square stencil kernels.

Based on the size of the stencil kernel, the stencil decorator is
able to compute the size of the border in the output array.  If
applying the kernel to some element of input array would cause
an index to be out-of-bounds then that element belongs to the border
of the output array.  In the above example, points ``-1`` and ``+1`` are
accessed in each dimension and thus the output array has a border
of size one in all dimensions.

The parallel mode is able to infer kernel indices as constants from
simple expressions if possible. For example::

    @njit(parallel=True)
    def stencil_test(A):
        c = 2
        B = stencil(
            lambda a, c: 0.3 * (a[-c+1] + a[0] + a[c-1]))(A, c)
        return B


Stencil decorator options
=========================

.. note::
   The stencil decorator may be augmented in the future to provide additional
   mechanisms for border handling. At present, only one behaviour is
   implemented, ``"constant"`` (see ``func_or_mode`` below for details).

.. _stencil-neighborhood:

``neighborhood``
----------------

Sometimes it may be inconvenient to write the stencil kernel
exclusively with ``Integer`` literals.  For example, let us say we
would like to compute the trailing 30-day moving average of a
time series of data.  One could write
``(a[-29] + a[-28] + ... + a[-1] + a[0]) / 30`` but the stencil
decorator offers a more concise form using the ``neighborhood``
option::

   @stencil(neighborhood = ((-29, 0),))
   def kernel2(a):
       cumul = 0
       for i in range(-29, 1):
           cumul += a[i]
       return cumul / 30

The neighborhood option is a tuple of tuples.  The outer tuple's
length is equal to the number of dimensions of the input array.
The inner tuple's lengths are always two because
each element of the inner tuple corresponds to minimum and
maximum index offsets used in the corresponding dimension.

If a user specifies a neighborhood but the kernel accesses elements outside the
specified neighborhood, **the behavior is undefined.**

.. _stencil-mode:

``func_or_mode``
----------------

The optional ``func_or_mode`` parameter controls how the border of the output array
is handled.  Currently, there is only one supported value, ``"constant"``.
In ``constant`` mode, the stencil kernel is not applied in cases where
the kernel would access elements outside the valid range of the input
array.  In such cases, those elements in the output array are assigned
to a constant value, as specified by the ``cval`` parameter.

``cval``
--------

The optional cval parameter defaults to zero but can be set to any
desired value, which is then used for the border of the output array
if the ``func_or_mode`` parameter is set to ``constant``.  The cval parameter is
ignored in all other modes.  The type of the cval parameter must match
the return type of the stencil kernel.  If the user wishes the output
array to be constructed from a particular type then they should ensure
that the stencil kernel returns that type.

``standard_indexing``
---------------------

By default, all array accesses in a stencil kernel are processed as
relative indices as described above.  However, sometimes it may be
advantageous to pass an auxiliary array (e.g. an array of weights)
to a stencil kernel and have that array use standard Python indexing
rather than relative indexing.  For this purpose, there is the
stencil decorator option ``standard_indexing`` whose value is a
collection of strings whose names match those parameters to the
stencil function that are to be accessed with standard Python indexing
rather than relative indexing::

    @stencil(standard_indexing=("b",))
    def kernel3(a, b):
        return a[-1] * b[0] + a[0] + b[1]

``StencilFunc``
===============

The stencil decorator returns a callable object of type ``StencilFunc``.
``StencilFunc`` objects contains a number of attributes but the only one of
potential interest to users is the ``neighborhood`` attribute.
If the ``neighborhood`` option was passed to the stencil decorator then
the provided neighborhood is stored in this attribute.  Else, upon
first execution or compilation, the system calculates the neighborhood
as described above and then stores the computed neighborhood into this
attribute.  A user may then inspect the attribute if they wish to verify
that the calculated neighborhood is correct.

Stencil invocation options
==========================

Internally, the stencil decorator transforms the specified stencil
kernel into a regular Python function.  This function will have the
same parameters as specified in the stencil kernel definition but will
also include the following optional parameter.

.. _stencil-function-out:

``out``
-------

The optional ``out`` parameter is added to every stencil function
generated by Numba.  If specified, the ``out`` parameter tells
Numba that the user is providing their own pre-allocated array
to be used for the output of the stencil.  In this case, the
stencil function will not allocate its own output array.
Users should assure that the return type of the stencil kernel can
be safely cast to the element-type of the user-specified output array
following the `Numpy ufunc casting rules`_.

.. _`Numpy ufunc casting rules`: http://docs.scipy.org/doc/numpy/reference/ufuncs.html#casting-rules

An example usage is shown below::

   >>> import numpy as np
   >>> input_arr = np.arange(100).reshape((10, 10))
   >>> output_arr = np.full(input_arr.shape, 0.0)
   >>> kernel1(input_arr, out=output_arr)
.. _performance-tips:

Performance Tips
================

This is a short guide to features present in Numba that can help with obtaining
the best performance from code. Two examples are used, both are entirely
contrived and exist purely for pedagogical reasons to motivate discussion.
The first is the computation of the trigonometric identity
``cos(x)^2 + sin(x)^2``, the second is a simple element wise square root of a
vector with reduction over summation. All performance numbers are indicative
only and unless otherwise stated were taken from running on an Intel ``i7-4790``
CPU (4 hardware threads) with an input of ``np.arange(1.e7)``.

.. note::
   A reasonably effective approach to achieving high performance code is to
   profile the code running with real data and use that to guide performance
   tuning. The information presented here is to demonstrate features, not to act
   as canonical guidance!

No Python mode vs Object mode
-----------------------------

A common pattern is to decorate functions with ``@jit`` as this is the most
flexible decorator offered by Numba. ``@jit`` essentially encompasses two modes
of compilation, first it will try and compile the decorated function in no
Python mode, if this fails it will try again to compile the function using
object mode. Whilst the use of looplifting in object mode can enable some
performance increase, getting functions to compile under no python mode is
really the key to good performance. To make it such that only no python mode is
used and if compilation fails an exception is raised the decorators ``@njit``
and ``@jit(nopython=True)`` can be used (the first is an alias of the
second for convenience).

Loops
-----
Whilst NumPy has developed a strong idiom around the use of vector operations,
Numba is perfectly happy with loops too. For users familiar with C or Fortran,
writing Python in this style will work fine in Numba (after all, LLVM gets a
lot of use in compiling C lineage languages). For example::

    @njit
    def ident_np(x):
        return np.cos(x) ** 2 + np.sin(x) ** 2

    @njit
    def ident_loops(x):
        r = np.empty_like(x)
        n = len(x)
        for i in range(n):
            r[i] = np.cos(x[i]) ** 2 + np.sin(x[i]) ** 2
        return r

The above run at almost identical speeds when decorated with ``@njit``, without
the decorator the vectorized function is a couple of orders of magnitude faster.

+-----------------+-------+----------------+
| Function Name   | @njit | Execution time |
+=================+=======+================+
| ``ident_np``    | No    |     0.581s     |
+-----------------+-------+----------------+
| ``ident_np``    | Yes   |     0.659s     |
+-----------------+-------+----------------+
| ``ident_loops`` | No    |     25.2s      |
+-----------------+-------+----------------+
| ``ident_loops`` | Yes   |     0.670s     |
+-----------------+-------+----------------+

.. _fast-math:

Fastmath
--------
In certain classes of applications strict IEEE 754 compliance is less
important. As a result it is possible to relax some numerical rigour with
view of gaining additional performance. The way to achieve this behaviour in
Numba is through the use of the ``fastmath`` keyword argument::

    @njit(fastmath=False)
    def do_sum(A):
        acc = 0.
        # without fastmath, this loop must accumulate in strict order
        for x in A:
            acc += np.sqrt(x)
        return acc

    @njit(fastmath=True)
    def do_sum_fast(A):
        acc = 0.
        # with fastmath, the reduction can be vectorized as floating point
        # reassociation is permitted.
        for x in A:
            acc += np.sqrt(x)
        return acc


+-----------------+-----------------+
| Function Name   | Execution time  |
+=================+=================+
| ``do_sum``      |      35.2 ms    |
+-----------------+-----------------+
| ``do_sum_fast`` |      17.8 ms    |
+-----------------+-----------------+

In some cases you may wish to opt-in to only a subset of possible fast-math
optimizations. This can be done by supplying a set of `LLVM fast-math flags
<https://llvm.org/docs/LangRef.html#fast-math-flags>`_ to ``fastmath``.::

    def add_assoc(x, y):
        return (x - y) + y

    print(njit(fastmath=False)(add_assoc)(0, np.inf)) # nan
    print(njit(fastmath=True) (add_assoc)(0, np.inf)) # 0.0
    print(njit(fastmath={'reassoc', 'nsz'})(add_assoc)(0, np.inf)) # 0.0
    print(njit(fastmath={'reassoc'})       (add_assoc)(0, np.inf)) # nan
    print(njit(fastmath={'nsz'})           (add_assoc)(0, np.inf)) # nan


Parallel=True
-------------
If code contains operations that are parallelisable (:ref:`and supported
<numba-parallel-supported>`) Numba can compile a version that will run in
parallel on multiple native threads (no GIL!). This parallelisation is performed
automatically and is enabled by simply adding the ``parallel`` keyword
argument::

    @njit(parallel=True)
    def ident_parallel(x):
        return np.cos(x) ** 2 + np.sin(x) ** 2


Executions times are as follows:

+--------------------+-----------------+
| Function Name      | Execution time  |
+====================+=================+
| ``ident_parallel`` | 112 ms          |
+--------------------+-----------------+


The execution speed of this function with ``parallel=True`` present is
approximately 5x that of the NumPy equivalent and 6x that of standard
``@njit``.


Numba parallel execution also has support for explicit parallel loop
declaration similar to that in OpenMP. To indicate that a loop should be
executed in parallel the ``numba.prange`` function should be used, this function
behaves like Python ``range`` and if ``parallel=True`` is not set it acts
simply as an alias of ``range``. Loops induced with ``prange`` can be used for
embarrassingly parallel computation and also reductions.

Revisiting the reduce over sum example, assuming it is safe for the sum to be
accumulated out of order, the loop in ``n`` can be parallelised through the use
of ``prange``. Further, the ``fastmath=True`` keyword argument can be added
without concern in this case as the assumption that out of order execution is
valid has already been made through the use of ``parallel=True`` (as each thread
computes a partial sum).
::

    @njit(parallel=True)
    def do_sum_parallel(A):
        # each thread can accumulate its own partial sum, and then a cross
        # thread reduction is performed to obtain the result to return
        n = len(A)
        acc = 0.
        for i in prange(n):
            acc += np.sqrt(A[i])
        return acc

    @njit(parallel=True, fastmath=True)
    def do_sum_parallel_fast(A):
        n = len(A)
        acc = 0.
        for i in prange(n):
            acc += np.sqrt(A[i])
        return acc


Execution times are as follows, ``fastmath`` again improves performance.

+-------------------------+-----------------+
| Function Name           | Execution time  |
+=========================+=================+
| ``do_sum_parallel``     |      9.81 ms    |
+-------------------------+-----------------+
| ``do_sum_parallel_fast``|      5.37 ms    |
+-------------------------+-----------------+

.. _intel-svml:

Intel SVML
----------

Intel provides a short vector math library (SVML) that contains a large number
of optimised transcendental functions available for use as compiler
intrinsics. If the ``icc_rt`` package is present in the environment (or the SVML
libraries are simply locatable!) then Numba automatically configures the LLVM
back end to use the SVML intrinsic functions where ever possible. SVML provides
both high and low accuracy versions of each intrinsic and the version that is
used is determined through the use of the ``fastmath`` keyword. The default is
to use high accuracy which is accurate to within ``1 ULP``, however if 
``fastmath`` is set to ``True`` then the lower accuracy versions of the
intrinsics are used (answers to within ``4 ULP``).


First obtain SVML, using conda for example::

    conda install -c numba icc_rt

Rerunning the identity function example ``ident_np`` from above with various
combinations of options to ``@njit`` and with/without SVML yields the following
performance results (input size ``np.arange(1.e8)``). For reference, with just
NumPy the function executed in ``5.84s``:

+-----------------------------------+--------+-------------------+
| ``@njit`` kwargs                  |  SVML  | Execution time    |
+===================================+========+===================+
| ``None``                          | No     | 5.95s             |
+-----------------------------------+--------+-------------------+
| ``None``                          | Yes    | 2.26s             |
+-----------------------------------+--------+-------------------+
| ``fastmath=True``                 | No     | 5.97s             |
+-----------------------------------+--------+-------------------+
| ``fastmath=True``                 | Yes    | 1.8s              |
+-----------------------------------+--------+-------------------+
| ``parallel=True``                 | No     | 1.36s             |
+-----------------------------------+--------+-------------------+
| ``parallel=True``                 | Yes    | 0.624s            |
+-----------------------------------+--------+-------------------+
| ``parallel=True, fastmath=True``  | No     | 1.32s             |
+-----------------------------------+--------+-------------------+
| ``parallel=True, fastmath=True``  | Yes    | 0.576s            |
+-----------------------------------+--------+-------------------+

It is evident that SVML significantly increases the performance of this
function. The impact of ``fastmath`` in the case of SVML not being present is
zero, this is expected as there is nothing in the original function that would
benefit from relaxing numerical strictness.

Linear algebra
--------------
Numba supports most of ``numpy.linalg`` in no Python mode. The internal
implementation relies on a LAPACK and BLAS library to do the numerical work
and it obtains the bindings for the necessary functions from SciPy. Therefore,
to achieve good performance in ``numpy.linalg`` functions with Numba it is
necessary to use a SciPy built against a well optimised LAPACK/BLAS library.
In the case of the Anaconda distribution SciPy is built against Intel's MKL
which is highly optimised and as a result Numba makes use of this performance.
========
Examples
========


Mandelbrot
----------

.. literalinclude:: ../../../numba/tests/doc_examples/test_examples.py
   :language: python
   :caption: from ``test_mandelbrot`` of ``numba/tests/doc_examples/test_examples.py``
   :start-after: magictoken.ex_mandelbrot.begin
   :end-before: magictoken.ex_mandelbrot.end
   :dedent: 12
   :linenos:

.. _example-movemean:

Moving average
--------------

.. literalinclude:: ../../../numba/tests/doc_examples/test_examples.py
   :language: python
   :caption: from ``test_moving_average`` of ``numba/tests/doc_examples/test_examples.py``
   :start-after: magictoken.ex_moving_average.begin
   :end-before: magictoken.ex_moving_average.end
   :dedent: 12
   :linenos:

Multi-threading
---------------

The code below showcases the potential performance improvement when
using the :ref:`nogil <jit-nogil>` feature.  For example, on a 4-core machine,
the following results were printed::

   numpy (1 thread)       145 ms
   numba (1 thread)       128 ms
   numba (4 threads)       35 ms

.. note::
   If preferred it's possible to use the standard `concurrent.futures
   <https://docs.python.org/3/library/concurrent.futures.html>`_ module
   rather than spawn threads and dispatch tasks by hand.

.. literalinclude:: ../../../numba/tests/doc_examples/test_examples.py
   :language: python
   :caption: from ``test_no_gil`` of ``numba/tests/doc_examples/test_examples.py``
   :start-after: magictoken.ex_no_gil.begin
   :end-before: magictoken.ex_no_gil.end
   :dedent: 12
   :linenos:
.. _cfunc:

====================================
Creating C callbacks with ``@cfunc``
====================================

Interfacing with some native libraries (for example written in C or C++)
can necessitate writing native callbacks to provide business logic to the
library.  The :func:`numba.cfunc` decorator creates a compiled function
callable from foreign C code, using the signature of your choice.


Basic usage
===========

The ``@cfunc`` decorator has a similar usage to ``@jit``, but with an
important difference: passing a single signature is mandatory.
It determines the visible signature of the C callback::

   from numba import cfunc

   @cfunc("float64(float64, float64)")
   def add(x, y):
       return x + y


The C function object exposes the address of the compiled C callback as
the :attr:`~CFunc.address` attribute, so that you can pass it to any
foreign C or C++ library.  It also exposes a :mod:`ctypes` callback
object pointing to that callback; that object is also callable from
Python, making it easy to check the compiled code::

   @cfunc("float64(float64, float64)")
   def add(x, y):
       return x + y

   print(add.ctypes(4.0, 5.0))  # prints "9.0"


Example
=======

In this example, we are going to be using the ``scipy.integrate.quad``
function.  That function accepts either a regular Python callback or
a C callback wrapped in a :mod:`ctypes` callback object.

Let's define a pure Python integrand and compile it as a
C callback::

   >>> import numpy as np
   >>> from numba import cfunc
   >>> def integrand(t):
           return np.exp(-t) / t**2
      ...:
   >>> nb_integrand = cfunc("float64(float64)")(integrand)

We can pass the ``nb_integrand`` object's :mod:`ctypes` callback to
``scipy.integrate.quad`` and check that the results are the same as with
the pure Python function::

   >>> import scipy.integrate as si
   >>> def do_integrate(func):
           """
           Integrate the given function from 1.0 to +inf.
           """
           return si.quad(func, 1, np.inf)
      ...:
   >>> do_integrate(integrand)
   (0.14849550677592208, 3.8736750296130505e-10)
   >>> do_integrate(nb_integrand.ctypes)
   (0.14849550677592208, 3.8736750296130505e-10)


Using the compiled callback, the integration function does not invoke the
Python interpreter each time it evaluates the integrand.  In our case, the
integration is made 18 times faster::

   >>> %timeit do_integrate(integrand)
   1000 loops, best of 3: 242 µs per loop
   >>> %timeit do_integrate(nb_integrand.ctypes)
   100000 loops, best of 3: 13.5 µs per loop


Dealing with pointers and array memory
======================================

A less trivial use case of C callbacks involves doing operation on some
array of data passed by the caller.  As C doesn't have a high-level
abstraction similar to Numpy arrays, the C callback's signature will pass
low-level pointer and size arguments.  Nevertheless, the Python code for
the callback will expect to exploit the power and expressiveness of Numpy
arrays.

In the following example, the C callback is expected to operate on 2-d arrays,
with the signature ``void(double *input, double *output, int m, int n)``.
You can implement such a callback thusly::

   from numba import cfunc, types, carray

   c_sig = types.void(types.CPointer(types.double),
                      types.CPointer(types.double),
                      types.intc, types.intc)

   @cfunc(c_sig)
   def my_callback(in_, out, m, n):
       in_array = carray(in_, (m, n))
       out_array = carray(out, (m, n))
       for i in range(m):
           for j in range(n):
               out_array[i, j] = 2 * in_array[i, j]


The :func:`numba.carray` function takes as input a data pointer and a shape
and returns an array view of the given shape over that data.  The data is
assumed to be laid out in C order.  If the data is laid out in Fortran order,
:func:`numba.farray` should be used instead.


Handling C structures
=====================


With CFFI
---------

For applications that have a lot of state, it is useful to pass data in C
structures.  To simplify the interoperability with C code, numba can convert
a ``cffi`` type into a numba ``Record`` type using
``numba.core.typing.cffi_utils.map_type``::

   from numba.core.typing import cffi_utils

   nbtype = cffi_utils.map_type(cffi_type, use_record_dtype=True)

.. note:: **use_record_dtype=True** is needed otherwise pointers to C
    structures are returned as void pointers.

.. note:: From v0.49 the ``numba.cffi_support`` module has been phased out
    in favour of ``numba.core.typing.cffi_utils``


For example::

   from cffi import FFI

   src = """

   /* Define the C struct */
   typedef struct my_struct {
      int    i1;
      float  f2;
      double d3;
      float  af4[7]; // arrays are supported
   } my_struct;

   /* Define a callback function */
   typedef double (*my_func)(my_struct*, size_t);
   """

   ffi = FFI()
   ffi.cdef(src)

   # Get the function signature from *my_func*
   sig = cffi_utils.map_type(ffi.typeof('my_func'), use_record_dtype=True)

   # Make the cfunc
   from numba import cfunc, carray

   @cfunc(sig)
   def foo(ptr, n):
      base = carray(ptr, n)  # view pointer as an array of my_struct
      tmp = 0
      for i in range(n):
         tmp += base[i].i1 * base[i].f2 / base[i].d3
         tmp += base[i].af4.sum()  # nested arrays are like normal numpy array
      return tmp


With ``numba.types.Record.make_c_struct``
-----------------------------------------

The ``numba.types.Record`` type can be created manually to follow a
C-structure's layout.  To do that, use ``Record.make_c_struct``, for example::

   my_struct = types.Record.make_c_struct([
      # Provides a sequence of 2-tuples i.e. (name:str, type:Type)
      ('i1', types.int32),
      ('f2', types.float32),
      ('d3', types.float64),
      ('af4', types.NestedArray(dtype=types.float32, shape=(7,))),
   ])

Due to ABI limitations, structures should be passed as pointers
using ``types.CPointer(my_struct)`` as the argument type.  Inside the ``cfunc``
body, the ``my_struct*`` can be accessed with ``carray``.

Full example
------------

See full example in ``examples/notebooks/Accessing C Struct Data.ipynb``.


Signature specification
=======================

The explicit ``@cfunc`` signature can use any :ref:`Numba types <numba-types>`,
but only a subset of them make sense for a C callback.  You should
generally limit yourself to scalar types (such as ``int8`` or ``float64``)
,pointers to them (for example ``types.CPointer(types.int8)``), or pointers
to ``Record`` type.


Compilation options
===================

A number of keyword-only arguments can be passed to the ``@cfunc``
decorator: ``nopython`` and ``cache``.  Their meaning is similar to those
in the ``@jit`` decorator.

============================
Compiling code ahead of time
============================

.. _pycc:

While Numba's main use case is :term:`Just-in-Time compilation`, it also
provides a facility for :term:`Ahead-of-Time compilation` (AOT).


Overview
========

Benefits
--------

#. AOT compilation produces a compiled extension module which does not depend
   on Numba: you can distribute the module on machines which do not have
   Numba installed (but Numpy is required).

#. There is no compilation overhead at runtime (but see the
   ``@jit`` :ref:`cache <jit-cache>` option), nor any overhead of importing
   Numba.

.. seealso::
   Compiled extension modules are discussed in the
   `Python packaging user guide <https://packaging.python.org/en/latest/extensions/>`_.


Limitations
-----------

#. AOT compilation only allows for regular functions, not :term:`ufuncs <ufunc>`.

#. You have to specify function signatures explicitly.

#. Each exported function can have only one signature (but you can export
   several different signatures under different names).

#. Exported functions do not check the types of the arguments that are passed
   to them; the caller is expected to provide arguments of the correct type.

#. AOT compilation produces generic code for your CPU's architectural family
   (for example "x86-64"), while JIT compilation produces code optimized
   for your particular CPU model.


Usage
=====

Standalone example
------------------

::

   from numba.pycc import CC

   cc = CC('my_module')
   # Uncomment the following line to print out the compilation steps
   #cc.verbose = True

   @cc.export('multf', 'f8(f8, f8)')
   @cc.export('multi', 'i4(i4, i4)')
   def mult(a, b):
       return a * b

   @cc.export('square', 'f8(f8)')
   def square(a):
       return a ** 2

   if __name__ == "__main__":
       cc.compile()


If you run this Python script, it will generate an extension module named
``my_module``.  Depending on your platform, the actual filename may be
``my_module.so``, ``my_module.pyd``, ``my_module.cpython-34m.so``, etc.

The generated module has three functions: ``multf``, ``multi`` and ``square``.
``multi`` operates on 32-bit integers (``i4``), while ``multf`` and ``square``
operate on double-precision floats (``f8``)::

   >>> import my_module
   >>> my_module.multi(3, 4)
   12
   >>> my_module.square(1.414)
   1.9993959999999997


Distutils integration
---------------------

You can also integrate the compilation step for your extension modules
in your ``setup.py`` script, using distutils or setuptools::

   from distutils.core import setup

   from source_module import cc

   setup(...,
         ext_modules=[cc.distutils_extension()])


The ``source_module`` above is the module defining the ``cc`` object.
Extensions compiled like this will be automatically included in the
build files for your Python project, so you can distribute them inside
binary packages such as wheels or Conda packages. Note that in the case of
using conda, the compilers used for AOT need to be those that are available
in the Anaconda distribution.


Signature syntax
----------------

The syntax for exported signatures is the same as in the ``@jit``
decorator.  You can read more about it in the :ref:`types <numba-types>`
reference.

Here is an example of exporting an implementation of the second-order
centered difference on a 1d array::

   @cc.export('centdiff_1d', 'f8[:](f8[:], f8)')
   def centdiff_1d(u, dx):
       D = np.empty_like(u)
       D[0] = 0
       D[-1] = 0
       for i in range(1, len(D) - 1):
           D[i] = (u[i+1] - 2 * u[i] + u[i-1]) / dx**2
       return D

.. (example from http://nbviewer.ipython.org/gist/ketch/ae87a94f4ef0793d5d52)

You can also omit the return type, which will then be inferred by Numba::

   @cc.export('centdiff_1d', '(f8[:], f8)')
   def centdiff_1d(u, dx):
       # Same code as above
       ...

.. _jit-module:

============================================
Automatic module jitting with ``jit_module``
============================================

A common usage pattern is to have an entire module containing user-defined
functions that all need to be jitted. One option to accomplish this is to
manually apply the ``@jit`` decorator to each function definition. This approach
works and is great in many cases. However, for large modules with many functions,
manually ``jit``-wrapping each function definition can be tedious. For these
situations, Numba provides another option, the ``jit_module`` function, to
automatically replace functions declared in a module with their ``jit``-wrapped
equivalents.

It's important to note the conditions under which ``jit_module`` will *not*
impact a function:

1. Functions which have already been wrapped with a Numba decorator (e.g.
   ``jit``, ``vectorize``, ``cfunc``, etc.) are not impacted by ``jit_module``.

2. Functions which are declared outside the module from which ``jit_module``
   is called are not automatically ``jit``-wrapped.

3. Function declarations which occur logically after calling ``jit_module``
   are not impacted.

All other functions in a module will have the ``@jit`` decorator automatically
applied to them. See the following section for an example use case.

.. note:: This feature is for use by module authors. ``jit_module`` should not
    be called outside the context of a module containing functions to be jitted.


Example usage
=============

Let's assume we have a Python module we've created, ``mymodule.py`` (shown
below), which contains several functions. Some of these functions are defined
in ``mymodule.py`` while others are imported from other modules. We wish to have
all the functions which are defined in ``mymodule.py`` jitted using
``jit_module``.

.. _jit-module-usage:

.. code-block:: python

   # mymodule.py

   from numba import jit, jit_module

   def inc(x):
      return x + 1
   
   def add(x, y):
      return x + y
   
   import numpy as np
   # Use NumPy's mean function
   mean = np.mean
   
   @jit(nogil=True)
   def mul(a, b):
      return a * b
   
   jit_module(nopython=True, error_model="numpy")

   def div(a, b):
       return a / b

There are several things to note in the above example:

- Both the ``inc`` and ``add`` functions will be replaced with their
  ``jit``-wrapped equivalents with :ref:`compilation options <jit-options>`
  ``nopython=True`` and ``error_model="numpy"``.

- The ``mean`` function, because it's defined *outside* of ``mymodule.py`` in
  NumPy, will not be modified.

- ``mul`` will not be modified because it has been manually decorated with
  ``jit``.

- ``div`` will not be automatically ``jit``-wrapped because it is declared
  after ``jit_module`` is called.

When the above module is imported, we have:

.. code-block:: python

   >>> import mymodule
   >>> mymodule.inc
   CPUDispatcher(<function inc at 0x1032f86a8>)
   >>> mymodule.mean
   <function mean at 0x1096b8950>


API
===
.. warning:: This feature is experimental. The supported features may change
    with or without notice.

.. autofunction:: numba.jit_module

.. _cli:

Command line interface
======================

Numba is a Python package, usually you ``import numba`` from Python and use the
Python application programming interface (API). However, Numba also ships with a
command line interface (CLI), i.e. a tool ``numba`` that is installed when you
install Numba.

Currently, the only purpose of the CLI is to allow you to quickly show some
information about your system and installation, or to quickly get some debugging
information for a Python script using Numba.

.. _cli_usage:

Usage
-----

To use the Numba CLI from the terminal, use ``numba`` followed by the options
and arguments like ``--help`` or ``-s``, as explained below.

Sometimes it can happen that you get a "command not found" error when you type
``numba``, because your ``PATH`` isn't configured properly. In that case you can
use the equivalent command ``python -m numba``. If that still gives "command
not found", try to ``import numba`` as suggested here:
:ref:`numba-source-install-check`.

The two versions ``numba`` and ``python -m numba`` are the same. The first is
shorter to type, but if you get a "command not found" error because your
``PATH`` doesn't contain the location where ``numba`` is installed, having the
``python -m numba`` variant is useful.

To use the Numba CLI from IPython or Jupyter, use ``!numba``, i.e. prefix the
command with an exclamation mark. This is a general IPython/Jupyter feature to
execute shell commands, it is not available in the regular ``python`` terminal.

.. _cli_help:

Help
----

To see all available options, use ``numba --help``::

    $ numba --help
    usage: numba [-h] [--annotate] [--dump-llvm] [--dump-optimized]
                 [--dump-assembly] [--annotate-html ANNOTATE_HTML] [-s]
                 [--sys-json SYS_JSON]
                 [filename]

    positional arguments:
    filename              Python source filename

    optional arguments:
    -h, --help            show this help message and exit
    --annotate            Annotate source
    --dump-llvm           Print generated llvm assembly
    --dump-optimized      Dump the optimized llvm assembly
    --dump-assembly       Dump the LLVM generated assembly
    --annotate-html ANNOTATE_HTML
                            Output source annotation as html
    -s, --sysinfo         Output system information for bug reporting
    --sys-json SYS_JSON   Saves the system info dict as a json file


.. _cli_sysinfo:

System information
------------------

The ``numba -s`` (or the equivalent ``numba --sysinfo``) command prints a lot of
information about your system and your Numba installation and relevant
dependencies.

Remember: you can use ``!numba -s`` with an exclamation mark to see this
information from IPython or Jupyter.

Example output::

    $ numba -s

    System info:
    --------------------------------------------------------------------------------
    __Time Stamp__
    2019-05-07 14:15:39.733994

    __Hardware Information__
    Machine                                       : x86_64
    CPU Name                                      : haswell
    CPU count                                     : 8
    CPU Features                                  : 
    aes avx avx2 bmi bmi2 cmov cx16 f16c fma fsgsbase invpcid lzcnt mmx movbe pclmul
    popcnt rdrnd sahf sse sse2 sse3 sse4.1 sse4.2 ssse3 xsave xsaveopt

    __OS Information__
    Platform                                      : Darwin-18.5.0-x86_64-i386-64bit
    Release                                       : 18.5.0
    System Name                                   : Darwin
    Version                                       : Darwin Kernel Version 18.5.0: Mon Mar 11 20:40:32 PDT 2019; root:xnu-4903.251.3~3/RELEASE_X86_64
    OS specific info                              : 10.14.4   x86_64

    __Python Information__
    Python Compiler                               : Clang 4.0.1 (tags/RELEASE_401/final)
    Python Implementation                         : CPython
    Python Version                                : 3.7.3
    Python Locale                                 : en_US UTF-8

    __LLVM information__
    LLVM version                                  : 7.0.0

    __CUDA Information__
    CUDA driver library cannot be found or no CUDA enabled devices are present.
    Error class: <class 'numba.cuda.cudadrv.error.CudaSupportError'>

    __SVML Information__
    SVML state, config.USING_SVML                 : False
    SVML library found and loaded                 : False
    llvmlite using SVML patched LLVM              : True
    SVML operational                              : False

    __Threading Layer Information__
    TBB Threading layer available                 : False
    +--> Disabled due to                          : Unknown import problem.
    OpenMP Threading layer available              : False
    +--> Disabled due to                          : Unknown import problem.
    Workqueue Threading layer available           : True

    __Numba Environment Variable Information__
    None set.

    __Conda Information__
    conda_build_version                           : 3.17.8
    conda_env_version                             : 4.6.14
    platform                                      : osx-64
    python_version                                : 3.7.3.final.0
    root_writable                                 : True

    __Current Conda Env__
    (output truncated due to length)

.. _cli_debug:

Debugging
---------

As shown in the help output above, the ``numba`` command includes options that
can help you to debug Numba compiled code.

To try it out, create an example script called ``myscript.py``::

    import numba

    @numba.jit
    def f(x):
        return 2 * x
    
    f(42)

and then execute one of the following commands::

    $ numba myscript.py --annotate
    $ numba myscript.py --annotate-html myscript.html
    $ numba myscript.py --dump-llvm
    $ numba myscript.py --dump-optimized
    $ numba myscript.py --dump-assembly

.. _numba-troubleshooting:

========================
Troubleshooting and tips
========================

.. _what-to-compile:

What to compile
===============

The general recommendation is that you should only try to compile the
critical paths in your code.  If you have a piece of performance-critical
computational code amongst some higher-level code, you may factor out
the performance-critical code in a separate function and compile the
separate function with Numba.  Letting Numba focus on that small piece
of performance-critical code has several advantages:

* it reduces the risk of hitting unsupported features;
* it reduces the compilation times;
* it allows you to evolve the higher-level code which is outside of the
  compiled function much easier.

.. _code-doesnt-compile:

My code doesn't compile
=======================

There can be various reasons why Numba cannot compile your code, and raises
an error instead.  One common reason is that your code relies on an
unsupported Python feature, especially in :term:`nopython mode`.
Please see the list of :ref:`pysupported`.  If you find something that
is listed there and still fails compiling, please
:ref:`report a bug <report-numba-bugs>`.

When Numba tries to compile your code it first tries to work out the types of
all the variables in use, this is so it can generate a type specific
implementation of your code that can be compiled down to machine code. A common
reason for Numba failing to compile (especially in :term:`nopython mode`) is a
type inference failure, essentially Numba cannot work out what the type of all
the variables in your code should be.

For example, let's consider this trivial function::

    @jit(nopython=True)
    def f(x, y):
        return x + y

If you call it with two numbers, Numba is able to infer the types properly::

    >>> f(1, 2)
        3

If however you call it with a tuple and a number, Numba is unable to say
what the result of adding a tuple and number is, and therefore compilation
errors out::

    >>> f(1, (2,))
    Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    File "<path>/numba/numba/dispatcher.py", line 339, in _compile_for_args
        reraise(type(e), e, None)
    File "<path>/numba/numba/six.py", line 658, in reraise
        raise value.with_traceback(tb)
    numba.errors.TypingError: Failed at nopython (nopython frontend)
    Invalid use of + with parameters (int64, tuple(int64 x 1))
    Known signatures:
    * (int64, int64) -> int64
    * (int64, uint64) -> int64
    * (uint64, int64) -> int64
    * (uint64, uint64) -> uint64
    * (float32, float32) -> float32
    * (float64, float64) -> float64
    * (complex64, complex64) -> complex64
    * (complex128, complex128) -> complex128
    * (uint16,) -> uint64
    * (uint8,) -> uint64
    * (uint64,) -> uint64
    * (uint32,) -> uint64
    * (int16,) -> int64
    * (int64,) -> int64
    * (int8,) -> int64
    * (int32,) -> int64
    * (float32,) -> float32
    * (float64,) -> float64
    * (complex64,) -> complex64
    * (complex128,) -> complex128
    * parameterized
    [1] During: typing of intrinsic-call at <stdin> (3)

    File "<stdin>", line 3:

The error message helps you find out what went wrong:
"Invalid use of + with parameters (int64, tuple(int64 x 1))" is to be
interpreted as "Numba encountered an addition of variables typed as integer
and 1-tuple of integer, respectively, and doesn't know about any such
operation".

Note that if you allow object mode::

    @jit
    def g(x, y):
        return x + y

compilation will succeed and the compiled function will raise at runtime as
Python would do::

   >>> g(1, (2,))
   Traceback (most recent call last):
     File "<stdin>", line 1, in <module>
   TypeError: unsupported operand type(s) for +: 'int' and 'tuple'


My code has a type unification problem
======================================

Another common reason for Numba not being able to compile your code is that it
cannot statically determine the return type of a function. The most likely
cause of this is the return type depending on a value that is available only at
runtime. Again, this is most often problematic when using
:term:`nopython mode`. The concept of type unification is simply trying to find
a type in which two variables could safely be represented. For example a 64 bit
float and a 64 bit complex number could both be represented in a 128 bit complex
number.

As an example of type unification failure, this function has a return type that
is determined at runtime based on the value of `x`::

    In [1]: from numba import jit

    In [2]: @jit(nopython=True)
    ...: def f(x):
    ...:     if x > 10:
    ...:         return (1,)
    ...:     else:
    ...:         return 1
    ...:

    In [3]: f(10)

Trying to execute this function, errors out as follows::

    TypingError: Failed at nopython (nopython frontend)
    Can't unify return type from the following types: tuple(int64 x 1), int64
    Return of: IR name '$8.2', type '(int64 x 1)', location:
    File "<ipython-input-2-51ef1cc64bea>", line 4:
    def f(x):
        <source elided>
        if x > 10:
            return (1,)
            ^
    Return of: IR name '$12.2', type 'int64', location:
    File "<ipython-input-2-51ef1cc64bea>", line 6:
    def f(x):
        <source elided>
        else:
            return 1

The error message "Can't unify return type from the following types:
tuple(int64 x 1), int64" should be read as "Numba cannot find a type that
can safely represent a 1-tuple of integer and an integer".

.. _code-has-untyped-list:

My code has an untyped list problem
===================================

As :ref:`noted previously <code-doesnt-compile>` the first part of Numba
compiling your code involves working out what the types of all the variables
are. In the case of lists, a list must contain items that are of the same type
or can be empty if the type can be inferred from some later operation. What is
not possible is to have a list which is defined as empty and has no inferable
type (i.e. an untyped list).

For example, this is using a list of a known type::

    from numba import jit
    @jit(nopython=True)
    def f():
        return [1, 2, 3] # this list is defined on construction with `int` type

This is using an empty list, but the type can be inferred::

    from numba import jit
    @jit(nopython=True)
    def f(x):
        tmp = [] # defined empty
        for i in range(x):
            tmp.append(i) # list type can be inferred from the type of `i`
        return tmp

This is using an empty list and the type cannot be inferred::

    from numba import jit
    @jit(nopython=True)
    def f(x):
        tmp = [] # defined empty
        return (tmp, x) # ERROR: the type of `tmp` is unknown

Whilst slightly contrived, if you need an empty list and the type cannot be
inferred but you know what type you want the list to be, this "trick" can be
used to instruct the typing mechanism::

    from numba import jit
    import numpy as np
    @jit(nopython=True)
    def f(x):
        # define empty list, but instruct that the type is np.complex64
        tmp = [np.complex64(x) for x in range(0)]
        return (tmp, x) # the type of `tmp` is known, but it is still empty

The compiled code is too slow
=============================

The most common reason for slowness of a compiled JIT function is that
compiling in :term:`nopython mode` has failed and the Numba compiler has
fallen back to :term:`object mode`.  :term:`object mode` currently provides
little to no speedup compared to regular Python interpretation, and its
main point is to allow an internal optimization known as
:term:`loop-lifting`: this optimization will allow to compile inner
loops in :term:`nopython mode` regardless of what code surrounds those
inner loops.

To find out if type inference succeeded on your function, you can use
the :meth:`~Dispatcher.inspect_types` method on the compiled function.

For example, let's take the following function::

   @jit
   def f(a, b):
       s = a + float(b)
       return s

When called with numbers, this function should be fast as Numba is able
to convert number types to floating-point numbers.  Let's see::

   >>> f(1, 2)
   3.0
   >>> f.inspect_types()
   f (int64, int64)
   --------------------------------------------------------------------------------
   # --- LINE 7 ---

   @jit

   # --- LINE 8 ---

   def f(a, b):

       # --- LINE 9 ---
       # label 0
       #   a.1 = a  :: int64
       #   del a
       #   b.1 = b  :: int64
       #   del b
       #   $0.2 = global(float: <class 'float'>)  :: Function(<class 'float'>)
       #   $0.4 = call $0.2(b.1, )  :: (int64,) -> float64
       #   del b.1
       #   del $0.2
       #   $0.5 = a.1 + $0.4  :: float64
       #   del a.1
       #   del $0.4
       #   s = $0.5  :: float64
       #   del $0.5

       s = a + float(b)

       # --- LINE 10 ---
       #   $0.7 = cast(value=s)  :: float64
       #   del s
       #   return $0.7

       return s

Without trying to understand too much of the Numba intermediate representation,
it is still visible that all variables and temporary values have had their
types inferred properly: for example *a* has the type ``int64``, *$0.5* has
the type ``float64``, etc.

However, if *b* is passed as a string, compilation will fall back on object
mode as the float() constructor with a string is currently not supported
by Numba::

   >>> f(1, "2")
   3.0
   >>> f.inspect_types()
   [... snip annotations for other signatures, see above ...]
   ================================================================================
   f (int64, str)
   --------------------------------------------------------------------------------
   # --- LINE 7 ---

   @jit

   # --- LINE 8 ---

   def f(a, b):

       # --- LINE 9 ---
       # label 0
       #   a.1 = a  :: pyobject
       #   del a
       #   b.1 = b  :: pyobject
       #   del b
       #   $0.2 = global(float: <class 'float'>)  :: pyobject
       #   $0.4 = call $0.2(b.1, )  :: pyobject
       #   del b.1
       #   del $0.2
       #   $0.5 = a.1 + $0.4  :: pyobject
       #   del a.1
       #   del $0.4
       #   s = $0.5  :: pyobject
       #   del $0.5

       s = a + float(b)

       # --- LINE 10 ---
       #   $0.7 = cast(value=s)  :: pyobject
       #   del s
       #   return $0.7

       return s

Here we see that all variables end up typed as ``pyobject``.  This means
that the function was compiled in object mode and values are passed
around as generic Python objects, without Numba trying to look into them
to reason about their raw values.  This is a situation you want to avoid
when caring about the speed of your code.

If a function fails to compile in ``nopython`` mode warnings will be emitted
with explanation as to why compilation failed. For example with the ``f()``
function above (slightly edited for documentation purposes)::

    >>> f(1, 2)
    3.0
    >>> f(1, "2")
    example.py:7: NumbaWarning:
    Compilation is falling back to object mode WITH looplifting enabled because Function "f" failed type inference due to: Invalid use of Function(<class 'float'>) with argument(s) of type(s): (unicode_type)
    * parameterized
    In definition 0:
        TypeError: float() only support for numbers
        raised from <path>/numba/typing/builtins.py:880
    In definition 1:
        TypeError: float() only support for numbers
        raised from <path>/numba/typing/builtins.py:880
    This error is usually caused by passing an argument of a type that is unsupported by the named function.
    [1] During: resolving callee type: Function(<class 'float'>)
    [2] During: typing of call at example.py (9)


    File "example.py", line 9:
    def f(a, b):
        s = a + float(b)
        ^

    <path>/numba/compiler.py:722: NumbaWarning: Function "f" was compiled in object mode without forceobj=True.

    File "example.py", line 8:
    @jit
    def f(a, b):
    ^

    3.0


Disabling JIT compilation
=========================

In order to debug code, it is possible to disable JIT compilation, which makes
the ``jit`` decorator (and the ``njit`` decorator) act as if
they perform no operation, and the invocation of decorated functions calls the
original Python function instead of a compiled version. This can be toggled by
setting the :envvar:`NUMBA_DISABLE_JIT` enviroment variable to ``1``.

When this mode is enabled, the ``vectorize`` and ``guvectorize`` decorators will
still result in compilation of a ufunc, as there is no straightforward pure
Python implementation of these functions.


.. _debugging-jit-compiled-code:

Debugging JIT compiled code with GDB
====================================

Setting the ``debug`` keyword argument in the ``jit`` decorator
(e.g. ``@jit(debug=True)``) enables the emission of debug info in the jitted
code.  To debug, GDB version 7.0 or above is required.  Currently, the following
debug info is available:

* Function name will be shown in the backtrace along with type information and
  values (if available).
* Source location (filename and line number) is available.  For example,
  users can set a break point by the absolute filename and line number;
  e.g. ``break /path/to/myfile.py:6``.
* Arguments to the current function can be show with ``info args``
* Local variables in the current function can be shown with ``info locals``.
* The type of variables can be shown with ``whatis myvar``.
* The value of variables can be shown with ``print myvar`` or ``display myvar``.

  * Simple numeric types, i.e. int, float and double, are shown in their
    native representation.
  * Other types are shown as a structure based on Numba's memory model
    representation of the type.

Further, the Numba ``gdb`` printing extension can be loaded into ``gdb`` (if the
``gdb`` has Python support) to permit the printing of variables as they would be
in native Python. The extension does this by reinterpreting Numba's memory model
representations as Python types. Information about the ``gdb`` installation that
Numba is using, including the path to load the ``gdb`` printing extension, can
be displayed by using the ``numba -g`` command. For best results ensure that the
Python that ``gdb`` is using has a NumPy module accessible. An example output
of the ``gdb`` information follows:

.. code-block:: none
  :emphasize-lines: 1

    $ numba -g
    GDB info:
    --------------------------------------------------------------------------------
    Binary location                               : <some path>/gdb
    Print extension location                      : <some python path>/numba/misc/gdb_print_extension.py
    Python version                                : 3.8
    NumPy version                                 : 1.20.0
    Numba printing extension supported            : True

    To load the Numba gdb printing extension, execute the following from the gdb prompt:

    source <some python path>/numba/misc/gdb_print_extension.py

    --------------------------------------------------------------------------------

Known issues:

* Stepping depends heavily on optimization level. At full optimization
  (equivalent to O3), most of the variables are optimized out. It is often
  beneficial to use the jit option ``_dbg_optnone=True`` 
  or the environment variable :envvar:`NUMBA_OPT` to adjust the 
  optimization level and the jit option ``_dbg_extend_lifetimes=True`` 
  (which is on by default if ``debug=True``) or
  :envvar:`NUMBA_EXTEND_VARIABLE_LIFETIMES` to extend
  the lifetime of variables to the end of their scope so as to get a debugging
  experience closer to the semantics of Python execution.

* Memory consumption increases significantly with debug info enabled.
  The compiler emits extra information (`DWARF <http://www.dwarfstd.org/>`_)
  along with the instructions.  The emitted object code can be 2x bigger with
  debug info.

Internal details:

* Since Python semantics allow variables to bind to value of different types,
  Numba internally creates multiple versions of the variable for each type.
  So for code like::

    x = 1         # type int
    x = 2.3       # type float
    x = (1, 2, 3) # type 3-tuple of int

  Each assignments will store to a different variable name.  In the debugger,
  the variables will be ``x``, ``x$1`` and ``x$2``.  (In the Numba IR, they are
  ``x``, ``x.1`` and ``x.2``.)

* When debug is enabled, inlining of functions at LLVM IR level is disabled.

JIT options for debug
---------------------

* ``debug`` (bool). Set to ``True`` to enable debug info. Defaults to ``False``.
* ``_dbg_optnone`` (bool). Set to ``True`` to disable all LLVM optimization passes 
  on the function. Defaults to ``False``. See :envvar:`NUMBA_OPT` for a global setting
  to disable optimization.
* ``_dbg_extend_lifetimes`` (bool). Set to ``True`` to extend the lifetime of
  objects such that they more closely follow the semantics of Python.
  Automatically set to ``True`` when 
  ``debug=True``; otherwise, defaults to ``False``. Users can explicitly set this option 
  to ``False`` to retain the normal execution semantics of compiled code.
  See :envvar:`NUMBA_EXTEND_VARIABLE_LIFETIMES` for a global option to extend object 
  lifetimes.

Example debug usage
-------------------

The python source:

.. code-block:: python
  :linenos:

  from numba import njit

  @njit(debug=True)
  def foo(a):
      b = a + 1
      c = a * 2.34
      d = (a, b, c)
      print(a, b, c, d)

  r = foo(123)
  print(r)

In the terminal:

.. code-block:: none
  :emphasize-lines: 1, 3, 7, 12, 14, 16, 20, 22, 26, 28, 30, 32, 34, 36

    $ NUMBA_OPT=0 NUMBA_EXTEND_VARIABLE_LIFETIMES=1 gdb -q python
    Reading symbols from python...
    (gdb) break test1.py:5
    No source file named test1.py.
    Make breakpoint pending on future shared library load? (y or [n]) y
    Breakpoint 1 (test1.py:5) pending.
    (gdb) run test1.py
    Starting program: <path>/bin/python test1.py
    ...
    Breakpoint 1, __main__::foo_241[abi:c8tJTC_2fWgEeGLSgydRTQUgiqKEZ6gEoDvQJmaQIA](long long) (a=123) at test1.py:5
    5           b = a + 1
    (gdb) info args
    a = 123
    (gdb) n
    6           c = a * 2.34
    (gdb) info locals
    b = 124
    c = 0
    d = {f0 = 0, f1 = 0, f2 = 0}
    (gdb) n
    7           d = (a, b, c)
    (gdb) info locals
    b = 124
    c = 287.81999999999999
    d = {f0 = 0, f1 = 0, f2 = 0}
    (gdb) whatis b
    type = int64
    (gdb) whatis d
    type = Tuple(int64, int64, float64) ({i64, i64, double})
    (gdb) n
    8           print(a, b, c, d)
    (gdb) print b
    $1 = 124
    (gdb) print d
    $2 = {f0 = 123, f1 = 124, f2 = 287.81999999999999}
    (gdb) bt
    #0  __main__::foo_241[abi:c8tJTC_2fWgEeGLSgydRTQUgiqKEZ6gEoDvQJmaQIA](long long) (a=123) at test1.py:8
    #1  0x00007ffff06439fa in cpython::__main__::foo_241[abi:c8tJTC_2fWgEeGLSgydRTQUgiqKEZ6gEoDvQJmaQIA](long long) ()


Another example follows that makes use of the Numba ``gdb`` printing extension
mentioned above, note the change in the print format once the extension is
loaded with ``source`` :

The Python source:

.. code-block:: python
  :linenos:

    from numba import njit
    import numpy as np

    @njit(debug=True)
    def foo(n):
        x = np.arange(n)
        y = (x[0], x[-1])
        return x, y

    foo(4)

In the terminal:

.. code-block:: none
  :emphasize-lines: 1, 3, 4, 7, 12, 14, 16, 17, 20

    $ NUMBA_OPT=0 NUMBA_EXTEND_VARIABLE_LIFETIMES=1 gdb -q python
    Reading symbols from python...
    (gdb) set breakpoint pending on
    (gdb) break test2.py:8
    No source file named test2.py.
    Breakpoint 1 (test2.py:8) pending.
    (gdb) run test2.py
    Starting program: <path>/bin/python test2.py
    ...
    Breakpoint 1, __main__::foo_241[abi:c8tJTC_2fWgEeGLSgydRTQUgiqKEZ6gEoDvQJmaQIA](long long) (n=4) at test2.py:8
    8           return x, y
    (gdb) print x
    $1 = {meminfo = 0x55555688f470 "\001", parent = 0x0, nitems = 4, itemsize = 8, data = 0x55555688f4a0, shape = {4}, strides = {8}}
    (gdb) print y
    $2 = {0, 3}
    (gdb) source numba/misc/gdb_print_extension.py
    (gdb) print x
    $3 =
    [0 1 2 3]
    (gdb) print y
    $4 = (0, 3)



Globally override debug setting
-------------------------------

It is possible to enable debug for the full application by setting environment
variable ``NUMBA_DEBUGINFO=1``.  This sets the default value of the ``debug``
option in ``jit``.  Debug can be turned off on individual functions by setting
``debug=False``.

Beware that enabling debug info significantly increases the memory consumption
for each compiled function.  For large application, this may cause out-of-memory
error.

Using Numba's direct ``gdb`` bindings in ``nopython``  mode
===========================================================

Numba (version 0.42.0 and later) has some additional functions relating to
``gdb`` support for CPUs that make it easier to debug programs. All the ``gdb``
related functions described in the following work in the same manner
irrespective of whether they are called from the standard CPython interpreter or
code compiled in either :term:`nopython mode` or :term:`object mode`.

.. note:: This feature is experimental!

.. warning:: This feature does unexpected things if used from Jupyter or
             alongside the ``pdb`` module. It's behaviour is harmless, just hard
             to predict!

Set up
------
Numba's ``gdb`` related functions make use of a ``gdb`` binary, the location and
name of this binary can be configured via the :envvar:`NUMBA_GDB_BINARY`
environment variable if desired.

.. note:: Numba's ``gdb`` support requires the ability for ``gdb`` to attach to
          another process. On some systems (notably Ubuntu Linux) default
          security restrictions placed on ``ptrace`` prevent this from being
          possible. This restriction is enforced at the system level by the
          Linux security module `Yama`. Documentation for this module and the
          security implications of making changes to its behaviour can be found
          in the `Linux Kernel documentation <https://www.kernel.org/doc/Documentation/admin-guide/LSM/Yama.rst>`_.
          The `Ubuntu Linux security documentation <https://wiki.ubuntu.com/Security/Features#ptrace>`_
          discusses how to adjust the behaviour of `Yama` on with regards to
          ``ptrace_scope`` so as to permit the required behaviour.

Basic ``gdb`` support
---------------------

.. warning:: Calling :func:`numba.gdb` and/or :func:`numba.gdb_init` more than
             once in the same program is not advisable, unexpected things may
             happen. If multiple breakpoints are desired within a program,
             launch ``gdb`` once via :func:`numba.gdb` or :func:`numba.gdb_init`
             and then use :func:`numba.gdb_breakpoint` to register additional
             breakpoint locations.

The most simple function for adding ``gdb`` support is :func:`numba.gdb`, which,
at the call location, will:

* launch ``gdb`` and attach it to the running process.
* create a breakpoint at the site of the :func:`numba.gdb()` function call, the
  attached ``gdb`` will pause execution here awaiting user input.

use of this functionality is best motivated by example, continuing with the
example used above:

.. code-block:: python
  :linenos:

  from numba import njit, gdb

  @njit(debug=True)
  def foo(a):
      b = a + 1
      gdb() # instruct Numba to attach gdb at this location and pause execution
      c = a * 2.34
      d = (a, b, c)
      print(a, b, c, d)

  r= foo(123)
  print(r)

In the terminal (``...`` on a line by itself indicates output that is not
presented for brevity):

.. code-block:: none
    :emphasize-lines: 1, 4, 8, 13, 24, 26, 28, 30, 32, 37

    $ NUMBA_OPT=0 NUMBA_EXTEND_VARIABLE_LIFETIMES=1 python demo_gdb.py
    ...
    Breakpoint 1, 0x00007fb75238d830 in numba_gdb_breakpoint () from numba/_helperlib.cpython-39-x86_64-linux-gnu.so
    (gdb) s
    Single stepping until exit from function numba_gdb_breakpoint,
    which has no line number information.
    0x00007fb75233e1cf in numba::misc::gdb_hook::hook_gdb::_3clocals_3e::impl_242[abi:c8tJTIeFCjyCbUFRqqOAK_2f6h0phxApMogijRBAA_3d](StarArgTuple) ()
    (gdb) s
    Single stepping until exit from function _ZN5numba4misc8gdb_hook8hook_gdb12_3clocals_3e8impl_242B44c8tJTIeFCjyCbUFRqqOAK_2f6h0phxApMogijRBAA_3dE12StarArgTuple,
    which has no line number information.
    __main__::foo_241[abi:c8tJTC_2fWgEeGLSgydRTQUgiqKEZ6gEoDvQJmaQIA](long long) (a=123) at demo_gdb.py:7
    7           c = a * 2.34
    (gdb) l
    2
    3       @njit(debug=True)
    4       def foo(a):
    5           b = a + 1
    6           gdb() # instruct Numba to attach gdb at this location and pause execution
    7           c = a * 2.34
    8           d = (a, b, c)
    9           print(a, b, c, d)
    10
    11      r= foo(123)
    (gdb) p a
    $1 = 123
    (gdb) p b
    $2 = 124
    (gdb) p c
    $3 = 0
    (gdb) b 9
    Breakpoint 2 at 0x7fb73d1f7287: file demo_gdb.py, line 9.
    (gdb) c
    Continuing.

    Breakpoint 2, __main__::foo_241[abi:c8tJTC_2fWgEeGLSgydRTQUgiqKEZ6gEoDvQJmaQIA](long long) (a=123) at demo_gdb.py:9
    9           print(a, b, c, d)
    (gdb) info locals
    b = 124
    c = 287.81999999999999
    d = {f0 = 123, f1 = 124, f2 = 287.81999999999999}


It can be seen in the above example that execution of the code is paused at the
location of the ``gdb()`` function call at end of the ``numba_gdb_breakpoint``
function (this is the Numba internal symbol registered as breakpoint with
``gdb``). Issuing a ``step`` twice at this point moves to the stack frame of the
compiled Python source. From there, it can be seen that the variables ``a`` and
``b`` have been evaluated but ``c`` has not, as demonstrated by printing their
values, this is precisely as expected given the location of the ``gdb()`` call.
Issuing a ``break`` on line 9 and then continuing execution leads to the
evaluation of line ``7``. The variable ``c`` is assigned a value as a result of
the execution and this can be seen in output of ``info locals`` when the
breakpoint is hit.

Running with ``gdb`` enabled
----------------------------

The functionality provided by :func:`numba.gdb` (launch and attach ``gdb`` to
the executing process and pause on a breakpoint) is also available as two
separate functions:

* :func:`numba.gdb_init` this function injects code at the call site to launch
  and attach ``gdb`` to the executing process but does not pause execution.
* :func:`numba.gdb_breakpoint` this function injects code at the call site that
  will call the special ``numba_gdb_breakpoint`` function that is registered as
  a breakpoint in Numba's ``gdb`` support. This is demonstrated in the next
  section.

This functionality enables more complex debugging capabilities. Again, motivated
by example, debugging a 'segfault' (memory access violation signalling
``SIGSEGV``):

.. code-block:: python
  :linenos:

    from numba import njit, gdb_init
    import numpy as np

    # NOTE debug=True switches bounds-checking on, but for the purposes of this
    # example it is explicitly turned off so that the out of bounds index is
    # not caught!
    @njit(debug=True, boundscheck=False)
    def foo(a, index):
        gdb_init() # instruct Numba to attach gdb at this location, but not to pause execution
        b = a + 1
        c = a * 2.34
        d = c[index] # access an address that is a) invalid b) out of the page
        print(a, b, c, d)

    bad_index = int(1e9) # this index is invalid
    z = np.arange(10)
    r = foo(z, bad_index)
    print(r)

In the terminal (``...`` on a line by itself indicates output that is not
presented for brevity):

.. code-block:: none
    :emphasize-lines: 1, 6, 8, 10, 12

    $ NUMBA_OPT=0 python demo_gdb_segfault.py
    ...
    Program received signal SIGSEGV, Segmentation fault.
    0x00007f5a4ca655eb in __main__::foo_241[abi:c8tJTC_2fWgEeGLSgydRTQUgiqKEZ6gEoDvQJmaQIA](Array<long long, 1, C, mutable, aligned>, long long) (a=..., index=1000000000) at demo_gdb_segfault.py:12
    12          d = c[index] # access an address that is a) invalid b) out of the page
    (gdb) p index
    $1 = 1000000000
    (gdb) p c
    $2 = {meminfo = 0x5586cfb95830 "\001", parent = 0x0, nitems = 10, itemsize = 8, data = 0x5586cfb95860, shape = {10}, strides = {8}}
    (gdb) whatis c
    type = array(float64, 1d, C) ({i8*, i8*, i64, i64, double*, [1 x i64], [1 x i64]})
    (gdb) p c.nitems
    $3 = 10

In the ``gdb`` output it can be noted that a ``SIGSEGV`` signal was caught, and
the line in which the access violation occurred is printed.

Continuing the example as a debugging session demonstration, first ``index``
can be printed, and it is evidently 1e9. Printing ``c`` shows that it is a
structure, so the type needs looking up and it can be seen that is it an
``array(float64, 1d, C)`` type. Given the segfault came from an invalid access
it would be informative to check the number of items in the array and compare
that to the index requested. Inspecting the ``nitems`` member of the structure
``c`` shows 10 items. It's therefore clear that the segfault comes from an
invalid access of index ``1000000000`` in an array containing ``10`` items.

Adding breakpoints to code
--------------------------

The next example demonstrates using multiple breakpoints that are defined
through the invocation of the :func:`numba.gdb_breakpoint` function:

.. code-block:: python
  :linenos:

  from numba import njit, gdb_init, gdb_breakpoint

  @njit(debug=True)
  def foo(a):
      gdb_init() # instruct Numba to attach gdb at this location
      b = a + 1
      gdb_breakpoint() # instruct gdb to break at this location
      c = a * 2.34
      d = (a, b, c)
      gdb_breakpoint() # and to break again at this location
      print(a, b, c, d)

  r= foo(123)
  print(r)

In the terminal (``...`` on a line by itself indicates output that is not
presented for brevity):

.. code-block:: none
    :emphasize-lines: 1, 4, 9, 20, 22, 24, 29, 31

    $ NUMBA_OPT=0 python demo_gdb_breakpoints.py
    ...
    Breakpoint 1, 0x00007fb65bb4c830 in numba_gdb_breakpoint () from numba/_helperlib.cpython-39-x86_64-linux-gnu.so
    (gdb) step
    Single stepping until exit from function numba_gdb_breakpoint,
    which has no line number information.
    __main__::foo_241[abi:c8tJTC_2fWgEeGLSgydRTQUgiqKEZ6gEoDvQJmaQIA](long long) (a=123) at demo_gdb_breakpoints.py:8
    8           c = a * 2.34
    (gdb) l
    3       @njit(debug=True)
    4       def foo(a):
    5           gdb_init() # instruct Numba to attach gdb at this location
    6           b = a + 1
    7           gdb_breakpoint() # instruct gdb to break at this location
    8           c = a * 2.34
    9           d = (a, b, c)
    10          gdb_breakpoint() # and to break again at this location
    11          print(a, b, c, d)
    12
    (gdb) p b
    $1 = 124
    (gdb) p c
    $2 = 0
    (gdb) c
    Continuing.

    Breakpoint 1, 0x00007fb65bb4c830 in numba_gdb_breakpoint ()
    from numba/_helperlib.cpython-39-x86_64-linux-gnu.so
    (gdb) step
    11          print(a, b, c, d)
    (gdb) p c
    $3 = 287.81999999999999

From the ``gdb`` output it can be seen that execution paused at line 8 as a
breakpoint was hit, and after a ``continue`` was issued, it broke again at line
11 where the next breakpoint was hit.

Debugging in parallel regions
-----------------------------

The follow example is quite involved, it executes with ``gdb`` instrumentation
from the outset as per the example above, but it also uses threads and makes use
of the breakpoint functionality. Further, the last iteration of the parallel
section calls the function ``work``, which is actually just a binding to
``glibc``'s ``free(3)`` in this case, but could equally be some involved
function that is presenting a segfault for unknown reasons.

.. code-block:: python
  :linenos:

    from numba import njit, prange, gdb_init, gdb_breakpoint
    import ctypes

    def get_free():
        lib = ctypes.cdll.LoadLibrary('libc.so.6')
        free_binding = lib.free
        free_binding.argtypes = [ctypes.c_void_p,]
        free_binding.restype = None
        return free_binding

    work = get_free()

    @njit(debug=True, parallel=True)
    def foo():
        gdb_init() # instruct Numba to attach gdb at this location, but not to pause execution
        counter = 0
        n = 9
        for i in prange(n):
            if i > 3 and i < 8: # iterations 4, 5, 6, 7 will break here
                gdb_breakpoint()

            if i == 8: # last iteration segfaults
                work(0xBADADD)

            counter += 1
        return counter

    r = foo()
    print(r)

In the terminal (``...`` on a line by itself indicates output that is not
presented for brevity), note the setting of ``NUMBA_NUM_THREADS`` to 4 to ensure
that there are 4 threads running in the parallel section:

.. code-block:: none
    :emphasize-lines: 1, 19, 29, 44, 50, 56, 62, 69

    $ NUMBA_NUM_THREADS=4 NUMBA_OPT=0 python demo_gdb_threads.py
    Attaching to PID: 21462
    ...
    Attaching to process 21462
    [New LWP 21467]
    [New LWP 21468]
    [New LWP 21469]
    [New LWP 21470]
    [Thread debugging using libthread_db enabled]
    Using host libthread_db library "/lib64/libthread_db.so.1".
    0x00007f59ec31756d in nanosleep () at ../sysdeps/unix/syscall-template.S:81
    81      T_PSEUDO (SYSCALL_SYMBOL, SYSCALL_NAME, SYSCALL_NARGS)
    Breakpoint 1 at 0x7f59d631e8f0: file numba/_helperlib.c, line 1090.
    Continuing.
    [Switching to Thread 0x7f59d1fd1700 (LWP 21470)]

    Thread 5 "python" hit Breakpoint 1, numba_gdb_breakpoint () at numba/_helperlib.c:1090
    1090    }
    (gdb) info threads
    Id   Target Id         Frame
    1    Thread 0x7f59eca2f740 (LWP 21462) "python" pthread_cond_wait@@GLIBC_2.3.2 ()
        at ../nptl/sysdeps/unix/sysv/linux/x86_64/pthread_cond_wait.S:185
    2    Thread 0x7f59d37d4700 (LWP 21467) "python" pthread_cond_wait@@GLIBC_2.3.2 ()
        at ../nptl/sysdeps/unix/sysv/linux/x86_64/pthread_cond_wait.S:185
    3    Thread 0x7f59d2fd3700 (LWP 21468) "python" pthread_cond_wait@@GLIBC_2.3.2 ()
        at ../nptl/sysdeps/unix/sysv/linux/x86_64/pthread_cond_wait.S:185
    4    Thread 0x7f59d27d2700 (LWP 21469) "python" numba_gdb_breakpoint () at numba/_helperlib.c:1090
    * 5    Thread 0x7f59d1fd1700 (LWP 21470) "python" numba_gdb_breakpoint () at numba/_helperlib.c:1090
    (gdb) thread apply 2-5 info locals

    Thread 2 (Thread 0x7f59d37d4700 (LWP 21467)):
    No locals.

    Thread 3 (Thread 0x7f59d2fd3700 (LWP 21468)):
    No locals.

    Thread 4 (Thread 0x7f59d27d2700 (LWP 21469)):
    No locals.

    Thread 5 (Thread 0x7f59d1fd1700 (LWP 21470)):
    sched$35 = '\000' <repeats 55 times>
    counter__arr = '\000' <repeats 16 times>, "\001\000\000\000\000\000\000\000\b\000\000\000\000\000\000\000\370B]\"hU\000\000\001", '\000' <repeats 14 times>
    counter = 0
    (gdb) continue
    Continuing.
    [Switching to Thread 0x7f59d27d2700 (LWP 21469)]

    Thread 4 "python" hit Breakpoint 1, numba_gdb_breakpoint () at numba/_helperlib.c:1090
    1090    }
    (gdb) continue
    Continuing.
    [Switching to Thread 0x7f59d1fd1700 (LWP 21470)]

    Thread 5 "python" hit Breakpoint 1, numba_gdb_breakpoint () at numba/_helperlib.c:1090
    1090    }
    (gdb) continue
    Continuing.
    [Switching to Thread 0x7f59d27d2700 (LWP 21469)]

    Thread 4 "python" hit Breakpoint 1, numba_gdb_breakpoint () at numba/_helperlib.c:1090
    1090    }
    (gdb) continue
    Continuing.

    Thread 5 "python" received signal SIGSEGV, Segmentation fault.
    [Switching to Thread 0x7f59d1fd1700 (LWP 21470)]
    __GI___libc_free (mem=0xbadadd) at malloc.c:2935
    2935      if (chunk_is_mmapped(p))                       /* release mmapped memory. */
    (gdb) bt
    #0  __GI___libc_free (mem=0xbadadd) at malloc.c:2935
    #1  0x00007f59d37ded84 in $3cdynamic$3e::__numba_parfor_gufunc__0x7ffff80a61ae3e31$244(Array<unsigned long long, 1, C, mutable, aligned>, Array<long long, 1, C, mutable, aligned>) () at <string>:24
    #2  0x00007f59d17ce326 in __gufunc__._ZN13$3cdynamic$3e45__numba_parfor_gufunc__0x7ffff80a61ae3e31$244E5ArrayIyLi1E1C7mutable7alignedE5ArrayIxLi1E1C7mutable7alignedE ()
    #3  0x00007f59d37d7320 in thread_worker ()
    from <path>/numba/numba/npyufunc/workqueue.cpython-37m-x86_64-linux-gnu.so
    #4  0x00007f59ec626e25 in start_thread (arg=0x7f59d1fd1700) at pthread_create.c:308
    #5  0x00007f59ec350bad in clone () at ../sysdeps/unix/sysv/linux/x86_64/clone.S:113

In the output it can be seen that there are 4 threads launched and that they all
break at the breakpoint, further that ``Thread 5`` receives a signal ``SIGSEGV``
and that back tracing shows that it came from ``__GI___libc_free`` with the
invalid address in ``mem``, as expected.

Using the ``gdb`` command language
----------------------------------
Both the :func:`numba.gdb` and :func:`numba.gdb_init` functions accept unlimited
string arguments which will be passed directly to ``gdb`` as command line
arguments when it initializes, this makes it easy to set breakpoints on other
functions and perform repeated debugging tasks without having to manually type
them every time. For example, this code runs with ``gdb`` attached and sets a
breakpoint on ``_dgesdd`` (say for example the arguments passed to the LAPACK's
double precision divide and conqueror SVD function need debugging).

.. code-block:: python
  :linenos:

    from numba import njit, gdb
    import numpy as np

    @njit(debug=True)
    def foo(a):
        # instruct Numba to attach gdb at this location and on launch, switch
        # breakpoint pending on , and then set a breakpoint on the function
        # _dgesdd, continue execution, and once the breakpoint is hit, backtrace
        gdb('-ex', 'set breakpoint pending on',
            '-ex', 'b dgesdd_',
            '-ex','c',
            '-ex','bt')
        b = a + 10
        u, s, vh = np.linalg.svd(b)
        return s # just return singular values

    z = np.arange(70.).reshape(10, 7)
    r = foo(z)
    print(r)

In the terminal (``...`` on a line by itself indicates output that is not
presented for brevity), note that no interaction is required to break and
backtrace:

.. code-block:: none
    :emphasize-lines: 1

    $ NUMBA_OPT=0 python demo_gdb_args.py
    Attaching to PID: 22300
    GNU gdb (GDB) Red Hat Enterprise Linux 8.0.1-36.el7
    ...
    Attaching to process 22300
    Reading symbols from <py_env>/bin/python3.7...done.
    0x00007f652305a550 in __nanosleep_nocancel () at ../sysdeps/unix/syscall-template.S:81
    81      T_PSEUDO (SYSCALL_SYMBOL, SYSCALL_NAME, SYSCALL_NARGS)
    Breakpoint 1 at 0x7f650d0618f0: file numba/_helperlib.c, line 1090.
    Continuing.

    Breakpoint 1, numba_gdb_breakpoint () at numba/_helperlib.c:1090
    1090    }
    Breakpoint 2 at 0x7f65102322e0 (2 locations)
    Continuing.

    Breakpoint 2, 0x00007f65182be5f0 in mkl_lapack.dgesdd_ ()
    from <py_env>/lib/python3.7/site-packages/numpy/core/../../../../libmkl_rt.so
    #0  0x00007f65182be5f0 in mkl_lapack.dgesdd_ ()
    from <py_env>/lib/python3.7/site-packages/numpy/core/../../../../libmkl_rt.so
    #1  0x00007f650d065b71 in numba_raw_rgesdd (kind=kind@entry=100 'd', jobz=<optimized out>, jobz@entry=65 'A', m=m@entry=10,
        n=n@entry=7, a=a@entry=0x561c6fbb20c0, lda=lda@entry=10, s=0x561c6facf3a0, u=0x561c6fb680e0, ldu=10, vt=0x561c6fd375c0,
        ldvt=7, work=0x7fff4c926c30, lwork=-1, iwork=0x7fff4c926c40, info=0x7fff4c926c20) at numba/_lapack.c:1277
    #2  0x00007f650d06768f in numba_ez_rgesdd (ldvt=7, vt=0x561c6fd375c0, ldu=10, u=0x561c6fb680e0, s=0x561c6facf3a0, lda=10,
        a=0x561c6fbb20c0, n=7, m=10, jobz=65 'A', kind=<optimized out>) at numba/_lapack.c:1307
    #3  numba_ez_gesdd (kind=<optimized out>, jobz=<optimized out>, m=10, n=7, a=0x561c6fbb20c0, lda=10, s=0x561c6facf3a0,
        u=0x561c6fb680e0, ldu=10, vt=0x561c6fd375c0, ldvt=7) at numba/_lapack.c:1477
    #4  0x00007f650a3147a3 in numba::targets::linalg::svd_impl::$3clocals$3e::svd_impl$243(Array<double, 2, C, mutable, aligned>, omitted$28default$3d1$29) ()
    #5  0x00007f650a1c0489 in __main__::foo$241(Array<double, 2, C, mutable, aligned>) () at demo_gdb_args.py:15
    #6  0x00007f650a1c2110 in cpython::__main__::foo$241(Array<double, 2, C, mutable, aligned>) ()
    #7  0x00007f650cd096a4 in call_cfunc ()
    from <path>/numba/numba/_dispatcher.cpython-37m-x86_64-linux-gnu.so
    ...


How does the ``gdb`` binding work?
----------------------------------
For advanced users and debuggers of Numba applications it's important to know
some of the internal implementation details of the outlined ``gdb`` bindings.
The :func:`numba.gdb` and :func:`numba.gdb_init` functions work by injecting the
following into the function's LLVM IR:

* At the call site of the function first inject a call to ``getpid(3)`` to get
  the PID of the executing process and store this for use later, then inject a
  ``fork(3)`` call:

  * In the parent:

    * Inject a call ``sleep(3)`` (hence the pause whilst ``gdb`` loads).
    * Inject a call to the ``numba_gdb_breakpoint`` function (only
      :func:`numba.gdb` does this).

  * In the child:

    * Inject a call to ``execl(3)`` with the arguments
      ``numba.config.GDB_BINARY``, the ``attach`` command and the PID recorded
      earlier. Numba has a special ``gdb`` command file that contains
      instructions to break on the symbol ``numba_gdb_breakpoint`` and then
      ``finish``, this is to make sure that the program stops on the
      breakpoint but the frame it stops in is the compiled Python frame (or
      one ``step`` away from, depending on optimisation). This command file is
      also added to the arguments and finally and any user specified arguments
      are added.

At the call site of a :func:`numba.gdb_breakpoint` a call is injected to the
special ``numba_gdb_breakpoint`` symbol, which is already registered and
instrumented as a place to break and ``finish`` immediately.

As a result of this, a e.g. :func:`numba.gdb` call will cause a fork in the
program, the parent will sleep whilst the child launches ``gdb`` and attaches it
to the parent and tells the parent to continue. The launched ``gdb`` has the
``numba_gdb_breakpoint`` symbol registered as a breakpoint and when the parent
continues and stops sleeping it will immediately call ``numba_gdb_breakpoint``
on which the child will break. Additional :func:`numba.gdb_breakpoint` calls
create calls to the registered breakpoint hence the program will also break at
these locations.

.. _debugging-cuda-python-code:

Debugging CUDA Python code
==========================

Using the simulator
-------------------

CUDA Python code can be run in the Python interpreter using the CUDA Simulator,
allowing it to be debugged with the Python debugger or with print statements. To
enable the CUDA simulator, set the environment variable
:envvar:`NUMBA_ENABLE_CUDASIM` to 1. For more information on the CUDA Simulator,
see :ref:`the CUDA Simulator documentation <simulator>`.


Debug Info
----------

By setting the ``debug`` argument to ``cuda.jit`` to ``True``
(``@cuda.jit(debug=True)``), Numba will emit source location in the compiled
CUDA code.  Unlike the CPU target, only filename and line information are
available, but no variable type information is emitted.  The information
is sufficient to debug memory error with
`cuda-memcheck <http://docs.nvidia.com/cuda/cuda-memcheck/index.html>`_.

For example, given the following cuda python code:

.. code-block:: python
  :linenos:

  import numpy as np
  from numba import cuda

  @cuda.jit(debug=True)
  def foo(arr):
      arr[cuda.threadIdx.x] = 1

  arr = np.arange(30)
  foo[1, 32](arr)   # more threads than array elements

We can use ``cuda-memcheck`` to find the memory error:

.. code-block:: none

  $ cuda-memcheck python chk_cuda_debug.py
  ========= CUDA-MEMCHECK
  ========= Invalid __global__ write of size 8
  =========     at 0x00000148 in /home/user/chk_cuda_debug.py:6:cudapy::__main__::foo$241(Array<__int64, int=1, C, mutable, aligned>)
  =========     by thread (31,0,0) in block (0,0,0)
  =========     Address 0x500a600f8 is out of bounds
  ...
  =========
  ========= Invalid __global__ write of size 8
  =========     at 0x00000148 in /home/user/chk_cuda_debug.py:6:cudapy::__main__::foo$241(Array<__int64, int=1, C, mutable, aligned>)
  =========     by thread (30,0,0) in block (0,0,0)
  =========     Address 0x500a600f0 is out of bounds
  ...
.. _jit:

===================================
Compiling Python code with ``@jit``
===================================

Numba provides several utilities for code generation, but its central
feature is the :func:`numba.jit` decorator.  Using this decorator, you can mark
a function for optimization by Numba's JIT compiler.  Various invocation
modes trigger differing compilation options and behaviours.


Basic usage
===========

.. _jit-lazy:

Lazy compilation
----------------

The recommended way to use the ``@jit`` decorator is to let Numba decide
when and how to optimize::

   from numba import jit

   @jit
   def f(x, y):
       # A somewhat trivial example
       return x + y

In this mode, compilation will be deferred until the first function
execution.  Numba will infer the argument types at call time, and generate
optimized code based on this information.  Numba will also be able to
compile separate specializations depending on the input types.  For example,
calling the ``f()`` function above with integer or complex numbers will
generate different code paths::

   >>> f(1, 2)
   3
   >>> f(1j, 2)
   (2+1j)

Eager compilation
-----------------

You can also tell Numba the function signature you are expecting.  The
function ``f()`` would now look like::

   from numba import jit, int32

   @jit(int32(int32, int32))
   def f(x, y):
       # A somewhat trivial example
       return x + y

``int32(int32, int32)`` is the function's signature.  In this case, the
corresponding specialization will be compiled by the ``@jit`` decorator,
and no other specialization will be allowed. This is useful if you want
fine-grained control over types chosen by the compiler (for example,
to use single-precision floats).

If you omit the return type, e.g. by writing ``(int32, int32)`` instead of
``int32(int32, int32)``, Numba will try to infer it for you.  Function
signatures can also be strings, and you can pass several of them as a list;
see the :func:`numba.jit` documentation for more details.

Of course, the compiled function gives the expected results::

   >>> f(1,2)
   3

and if we specified ``int32`` as return type, the higher-order bits get
discarded::

   >>> f(2**31, 2**31 + 1)
   1


Calling and inlining other functions
====================================

Numba-compiled functions can call other compiled functions.  The function
calls may even be inlined in the native code, depending on optimizer
heuristics.  For example::

   @jit
   def square(x):
       return x ** 2

   @jit
   def hypot(x, y):
       return math.sqrt(square(x) + square(y))

The ``@jit`` decorator *must* be added to any such library function,
otherwise Numba may generate much slower code.


Signature specifications
========================

Explicit ``@jit`` signatures can use a number of types.  Here are some
common ones:

* ``void`` is the return type of functions returning nothing (which
  actually return :const:`None` when called from Python)
* ``intp`` and ``uintp`` are pointer-sized integers (signed and unsigned,
  respectively)
* ``intc`` and ``uintc`` are equivalent to C ``int`` and ``unsigned int``
  integer types
* ``int8``, ``uint8``, ``int16``, ``uint16``, ``int32``, ``uint32``,
  ``int64``, ``uint64`` are fixed-width integers of the corresponding bit
  width (signed and unsigned)
* ``float32`` and ``float64`` are single- and double-precision floating-point
  numbers, respectively
* ``complex64`` and ``complex128`` are single- and double-precision complex
  numbers, respectively
* array types can be specified by indexing any numeric type, e.g. ``float32[:]``
  for a one-dimensional single-precision array or ``int8[:,:]`` for a
  two-dimensional array of 8-bit integers.


.. _jit-options:

Compilation options
===================

A number of keyword-only arguments can be passed to the ``@jit`` decorator.

.. _jit-nopython:

``nopython``
------------

Numba has two compilation modes: :term:`nopython mode` and
:term:`object mode`.  The former produces much faster code, but has
limitations that can force Numba to fall back to the latter.  To prevent
Numba from falling back, and instead raise an error, pass ``nopython=True``.

::

   @jit(nopython=True)
   def f(x, y):
       return x + y

.. seealso:: :ref:`numba-troubleshooting`

.. _jit-nogil:

``nogil``
---------

Whenever Numba optimizes Python code to native code that only works on
native types and variables (rather than Python objects), it is not necessary
anymore to hold Python's :py:term:`global interpreter lock` (GIL).
Numba will release the GIL when entering such a compiled function if you
passed ``nogil=True``.

::

   @jit(nogil=True)
   def f(x, y):
       return x + y

Code running with the GIL released runs concurrently with other
threads executing Python or Numba code (either the same compiled function,
or another one), allowing you to take advantage of multi-core systems.
This will not be possible if the function is compiled in :term:`object mode`.

When using ``nogil=True``, you'll have to be wary of the usual pitfalls
of multi-threaded programming (consistency, synchronization, race conditions,
etc.).

.. _jit-cache:

``cache``
---------

To avoid compilation times each time you invoke a Python program,
you can instruct Numba to write the result of function compilation into
a file-based cache.  This is done by passing ``cache=True``::

   @jit(cache=True)
   def f(x, y):
       return x + y

.. _parallel_jit_option:

``parallel``
------------

Enables automatic parallelization (and related optimizations) for those
operations in the function known to have parallel semantics.  For a list of
supported operations, see :ref:`numba-parallel`.  This feature is enabled by
passing ``parallel=True`` and must be used in conjunction with
``nopython=True``::

   @jit(nopython=True, parallel=True)
   def f(x, y):
       return x + y

.. seealso:: :ref:`numba-parallel`
.. _numba-threading-layer:

The Threading Layers
====================

This section is about the Numba threading layer, this is the library that is
used internally to perform the parallel execution that occurs through the use of
the ``parallel`` targets for CPUs, namely:

* The use of the ``parallel=True`` kwarg in ``@jit`` and ``@njit``.
* The use of the ``target='parallel'`` kwarg in ``@vectorize`` and
  ``@guvectorize``.

.. note::
    If a code base does not use the ``threading`` or ``multiprocessing``
    modules (or any other sort of parallelism) the defaults for the threading
    layer that ship with Numba will work well, no further action is required!


Which threading layers are available?
-------------------------------------
There are three threading layers available and they are named as follows:

* ``tbb`` - A threading layer backed by Intel TBB.
* ``omp`` - A threading layer backed by OpenMP.
* ``workqueue`` -A simple built-in work-sharing task scheduler.

In practice, the only threading layer guaranteed to be present is ``workqueue``.
The ``omp`` layer requires the presence of a suitable OpenMP runtime library.
The ``tbb`` layer requires the presence of Intel's TBB libraries, these can be
obtained via the conda command::

    $ conda install tbb

If you installed Numba with ``pip``, TBB can be enabled by running::

    $ pip install tbb

Due to compatibility issues with manylinux1 and other portability concerns,
the OpenMP threading layer is disabled in the Numba binary wheels on PyPI.

.. note::
    The default manner in which Numba searches for and loads a threading layer
    is tolerant of missing libraries, incompatible runtimes etc.


.. _numba-threading-layer-setting-mech:

Setting the threading layer
---------------------------


The threading layer is set via the environment variable
``NUMBA_THREADING_LAYER`` or through assignment to
``numba.config.THREADING_LAYER``. If the programmatic approach to setting the
threading layer is used it must occur logically before any Numba based
compilation for a parallel target has occurred. There are two approaches to
choosing a threading layer, the first is by selecting a threading layer that is
safe under various forms of parallel execution, the second is through explicit
selection via the threading layer name (e.g. ``tbb``).

Setting the threading layer selection priority
----------------------------------------------

By default the threading layers are searched in the order of ``'tbb'``,
``'omp'``, then ``'workqueue'``. To change this search order whilst
maintaining the selection of a threading layer based on availability, the
environment variable :envvar:`NUMBA_THREADING_LAYER_PRIORITY` can be used.

Note that it can also be set via
:py:data:`numba.config.THREADING_LAYER_PRIORITY`.
Similar to :py:data:`numba.config.THREADING_LAYER`,
it must occur logically before any Numba based
compilation for a parallel target has occurred.

For example, to instruct Numba to choose ``omp`` first if available,
then ``tbb`` and so on, set the environment variable as
``NUMBA_THREADING_LAYER_PRIORITY="omp tbb workqueue"``.
Or programmatically,
``numba.config.THREADING_LAYER_PRIORITY = ["omp", "tbb", "workqueue"]``.

Selecting a threading layer for safe parallel execution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Parallel execution is fundamentally derived from core Python libraries in four
forms (the first three also apply to code using parallel execution via other
means!):

* ``threads`` from the ``threading`` module.
* ``spawn`` ing processes from the ``multiprocessing`` module via ``spawn``
  (default on Windows, only available in Python 3.4+ on Unix)
* ``fork`` ing processes from the ``multiprocessing`` module via ``fork``
  (default on Unix).
* ``fork`` ing processes from the ``multiprocessing`` module through the use of
  a ``forkserver`` (only available in Python 3 on Unix). Essentially a new
  process is spawned and then forks are made from this new process on request.

Any library in use with these forms of parallelism must exhibit safe behaviour
under the given paradigm. As a result, the threading layer selection methods
are designed to provide a way to choose a threading layer library that is safe
for a given paradigm in an easy, cross platform and environment tolerant manner.
The options that can be supplied to the
:ref:`setting mechanisms <numba-threading-layer-setting-mech>` are as
follows:

* ``default`` provides no specific safety guarantee and is the default.
* ``safe`` is both fork and thread safe, this requires the ``tbb`` package
  (Intel TBB libraries) to be installed.
* ``forksafe`` provides a fork safe library.
* ``threadsafe`` provides a thread safe library.

To discover the threading layer that was selected, the function
``numba.threading_layer()`` may be called after parallel execution. For example,
on a Linux machine with no TBB installed::

    from numba import config, njit, threading_layer
    import numpy as np

    # set the threading layer before any parallel target compilation
    config.THREADING_LAYER = 'threadsafe'

    @njit(parallel=True)
    def foo(a, b):
        return a + b

    x = np.arange(10.)
    y = x.copy()

    # this will force the compilation of the function, select a threading layer
    # and then execute in parallel
    foo(x, y)

    # demonstrate the threading layer chosen
    print("Threading layer chosen: %s" % threading_layer())

which produces::

    Threading layer chosen: omp

and this makes sense as GNU OpenMP, as present on Linux, is thread safe.

Selecting a named threading layer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Advanced users may wish to select a specific threading layer for their use case,
this is done by directly supplying the threading layer name to the
:ref:`setting mechanisms <numba-threading-layer-setting-mech>`. The options
and requirements are as follows:

+----------------------+-----------+-------------------------------------------+
| Threading Layer Name | Platform  | Requirements                              |
+======================+===========+===========================================+
| ``tbb``              | All       | The ``tbb`` package (``$ conda install    |
|                      |           | tbb``)                                    |
+----------------------+-----------+-------------------------------------------+
| ``omp``              | Linux     | GNU OpenMP libraries (very likely this    |
|                      |           | will already exist)                       |
|                      |           |                                           |
|                      | Windows   | MS OpenMP libraries (very likely this will|
|                      |           | already exist)                            |
|                      |           |                                           |
|                      | OSX       | The ``intel-openmp`` package (``$ conda   |
|                      |           | install intel-openmp``)                   |
+----------------------+-----------+-------------------------------------------+
| ``workqueue``        | All       | None                                      |
+----------------------+-----------+-------------------------------------------+

Should the threading layer not load correctly Numba will detect this and provide
a hint about how to resolve the problem. It should also be noted that the Numba
diagnostic command ``numba -s`` has a section
``__Threading Layer Information__`` that reports on the availability of
threading layers in the current environment.


Extra notes
-----------
The threading layers have fairly complex interactions with CPython internals and
system level libraries, some additional things to note:

* The installation of Intel's TBB libraries vastly widens the options available
  in the threading layer selection process.
* On Linux, the ``omp`` threading layer is not fork safe due to the GNU OpenMP
  runtime library (``libgomp``) not being fork safe. If a fork occurs in a
  program that is using the ``omp`` threading layer, a detection mechanism is
  present that will try and gracefully terminate the forked child and print an
  error message to ``STDERR``.
* On systems with the ``fork(2)`` system call available, if the TBB backed
  threading layer is in use and a ``fork`` call is made from a thread other than
  the thread that launched TBB (typically the main thread) then this results in
  undefined behaviour and a warning will be displayed on ``STDERR``. As
  ``spawn`` is essentially ``fork`` followed by ``exec`` it is safe to ``spawn``
  from a non-main thread, but as this cannot be differentiated from just a
  ``fork`` call the warning message will still be displayed.
* On OSX, the ``intel-openmp`` package is required to enable the OpenMP based
  threading layer.

.. _setting_the_number_of_threads:

Setting the Number of Threads
-----------------------------

The number of threads used by numba is based on the number of CPU cores
available (see :obj:`numba.config.NUMBA_DEFAULT_NUM_THREADS`), but it can be
overridden with the :envvar:`NUMBA_NUM_THREADS` environment variable.

The total number of threads that numba launches is in the variable
:obj:`numba.config.NUMBA_NUM_THREADS`.

For some use cases, it may be desirable to set the number of threads to a
lower value, so that numba can be used with higher level parallelism.

The number of threads can be set dynamically at runtime using
:func:`numba.set_num_threads`. Note that :func:`~.set_num_threads` only allows
setting the number of threads to a smaller value than
:obj:`~.NUMBA_NUM_THREADS`. Numba always launches
:obj:`numba.config.NUMBA_NUM_THREADS` threads, but :func:`~.set_num_threads`
causes it to mask out unused threads so they aren't used in computations.

The current number of threads used by numba can be accessed with
:func:`numba.get_num_threads`. Both functions work inside of a jitted
function.

.. _numba-threading-layer-thread-masking:

Example of Limiting the Number of Threads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example, suppose the machine we are running on has 8 cores (so
:obj:`numba.config.NUMBA_NUM_THREADS` would be ``8``). Suppose we want to run
some code with ``@njit(parallel=True)``, but we also want to run our code
concurrently in 4 different processes. With the default number of threads,
each Python process would run 8 threads, for a total in 4*8 = 32 threads,
which is oversubscription for our 8 cores. We should rather limit each process
to 2 threads, so that the total will be 4*2 = 8, which matches our number of
physical cores.

There are two ways to do this. One is to set the :envvar:`NUMBA_NUM_THREADS`
environment variable to ``2``.

.. code:: bash

   $ NUMBA_NUM_THREADS=2 python ourcode.py

However, there are two downsides to this approach:

1. :envvar:`NUMBA_NUM_THREADS` must be set before Numba is imported, and
   ideally before Python is launched. As soon as Numba is imported the
   environment variable is read and that number of threads is locked in as the
   number of threads Numba launches.

2. If we want to later increase the number of threads used by the process, we
   cannot. :envvar:`NUMBA_NUM_THREADS` sets the *maximum* number of threads
   that are launched for a process. Calling :func:`~.set_num_threads()` with a
   value greater than :obj:`numba.config.NUMBA_NUM_THREADS` results in an
   error.

The advantage of this approach is that we can do it from outside of the
process without changing the code.

Another approach is to use the :func:`numba.set_num_threads` function in our code

.. code:: python

   from numba import njit, set_num_threads

   @njit(parallel=True)
   def func():
       ...

   set_num_threads(2)
   func()

If we call ``set_num_threads(2)`` before executing our parallel code, it has
the same effect as calling the process with ``NUMBA_NUM_THREADS=2``, in that
the parallel code will only execute on 2 threads. However, we can later call
``set_num_threads(8)`` to increase the number of threads back to the default
size. And we do not have to worry about setting it before Numba gets imported.
It only needs to be called before the parallel function is run.

API Reference
~~~~~~~~~~~~~

.. py:data:: numba.config.NUMBA_NUM_THREADS

   The total (maximum) number of threads launched by numba.

   Defaults to :obj:`numba.config.NUMBA_DEFAULT_NUM_THREADS`, but can be
   overridden with the :envvar:`NUMBA_NUM_THREADS` environment variable.

.. py:data:: numba.config.NUMBA_DEFAULT_NUM_THREADS

   The number of usable CPU cores on the system (as determined by
   ``len(os.sched_getaffinity(0))``, if supported by the OS, or
   ``multiprocessing.cpu_count()`` if not).
   This is the default value for :obj:`numba.config.NUMBA_NUM_THREADS` unless
   the :envvar:`NUMBA_NUM_THREADS` environment variable is set.

.. autofunction:: numba.set_num_threads

.. autofunction:: numba.get_num_threads

Talks and Tutorials
===================

.. note:: This is a selection of talks and tutorials that have been given by members of
    the Numba team as well as Numba users.  If you know of a Numba-related talk
    that should be included on this list, please `open an issue <https://github.com/numba/numba/issues>`_.

Talks on Numba
--------------

* AnacondaCON 2018 - Accelerating Scientific Workloads with Numba - Siu Kwan Lam (`Video <https://www.youtube.com/watch?v=6oXedk2tGfk>`__)
* `DIANA-HEP Meeting, 23 April 2018 <https://indico.cern.ch/event/709711/>`__ - Overview of Numba - Stan Seibert

Talks on Applications of Numba
------------------------------

* GPU Technology Conference 2016 - Accelerating a Spectral Algorithm for Plasma Physics with Python/Numba on GPU - Manuel Kirchen & Rémi Lehe (`Slides <http://on-demand.gputechconf.com/gtc/2016/presentation/s6353-manuel-kirchen-spectral-algorithm-plasma-physics.pdf>`__)
* `DIANA-HEP Meeting, 23 April 2018 <https://indico.cern.ch/event/709711/>`_ - Use of Numba in XENONnT - Chris Tunnell
* `DIANA-HEP Meeting, 23 April 2018 <https://indico.cern.ch/event/709711/>`_ - Extending Numba for HEP data types - Jim Pivarski
* STAC Summit, Nov 1 2017 - Scaling High-Performance Python with Minimal Effort - Ehsan Totoni (`Video <https://stacresearch.com/STAC-Summit-1-Nov-2017-Intel-Totoni>`__, `Slides <https://stacresearch.com/system/files/resource/files/STAC-Summit-1-Nov-2017-Intel-Totoni.pdf>`__)
* SciPy 2018 - UMAP: Uniform Manifold Approximation and Projection for Dimensional Reduction - Leland McInnes (`Video <https://www.youtube.com/watch?v=nq6iPZVUxZU>`__, `Github <https://github.com/lmcinnes/umap>`__)
* PyData Berlin 2018 - Extending Pandas using Apache Arrow and Numba - Uwe L. Korn (`Video <https://www.youtube.com/watch?v=tvmX8YAFK80>`__, `Blog <https://uwekorn.com/2018/08/03/use-numba-to-work-with-apache-arrow-in-pure-python.html>`__)
* FOSDEM 2019 - Extending Numba - Joris Geessels  (`Video, Slides & Examples <https://fosdem.org/2019/schedule/event/python_extending_numba/>`__)
* PyCon India 2019 - Real World Numba: Taking the Path of Least Resistance - Ankit Mahato (`Video <https://www.youtube.com/watch?v=rhbegsr8stc>`__)
* SciPy 2019 - How to Accelerate an Existing Codebase with Numba - Siu Kwan Lam & Stanley Seibert (`Video <https://www.youtube.com/watch?v=-4tD8kNHdXs>`__)
* SciPy 2019 - Real World Numba: Creating a Skeleton Analysis Library - Juan Nunez-Iglesias (`Video <https://www.youtube.com/watch?v=0pUPNMglnaE>`__)
* SciPy 2019 - Fast Gradient Boosting Decision Trees with PyGBM and Numba - Nicholas Hug (`Video <https://www.youtube.com/watch?v=cLpIh8Aiy2w>`__)
* PyCon Sweden 2020 - Accelerating Scientific Computing using Numba - Ankit Mahato (`Video <https://www.youtube.com/watch?v=d_21Q0UoWrQ>`__)

Tutorials
---------

* SciPy 2017 - Numba: Tell those C++ Bullies to Get Lost - Gil Forsyth & Lorena Barba (`Video <https://www.youtube.com/watch?v=1AwG0T4gaO0>`__, `Notebooks <https://github.com/gforsyth/numba_tutorial_scipy2017>`__)
* GPU Technology Conference 2018 - GPU Computing in Python with Numba - Stan Seibert (`Notebooks <https://github.com/ContinuumIO/gtc2018-numba>`__)
* PyData Amsterdam 2019 - Create CUDA kernels from Python using Numba and CuPy - Valentin Haenel (`Video <https://www.youtube.com/watch?v=CQDsT81GyS8>`__)
.. _jitclass:

===========================================
Compiling Python classes with ``@jitclass``
===========================================

.. note::

  This is a early version of jitclass support. Not all compiling features are
  exposed or implemented, yet.


Numba supports code generation for classes via the
:func:`numba.experimental.jitclass` decorator.  A class can be marked for
optimization using this decorator along with a specification of the types of
each field.  We call the resulting class object a *jitclass*.  All methods of a
jitclass are compiled into nopython functions.  The data of a jitclass instance
is allocated on the heap as a C-compatible structure so that any compiled
functions can have direct access to the underlying data, bypassing the
interpreter.


Basic usage
===========

Here's an example of a jitclass:

.. literalinclude:: ../../../numba/tests/doc_examples/test_jitclass.py
   :language: python
   :start-after: magictoken.ex_jitclass.begin
   :end-before: magictoken.ex_jitclass.end
   :dedent: 8

In the above example, a ``spec`` is provided as a list of 2-tuples.  The tuples
contain the name of the field and the Numba type of the field.  Alternatively,
user can use a dictionary (an ``OrderedDict`` preferably for stable field
ordering), which maps field names to types.

The definition of the class requires at least a ``__init__`` method for
initializing each defined fields.  Uninitialized fields contains garbage data.
Methods and properties (getters and setters only) can be defined.  They will be
automatically compiled.


Inferred class member types from type annotations with ``as_numba_type``
========================================================================

Fields of a ``jitclass`` can also be inferred from Python type annotations.

.. literalinclude:: ../../../numba/tests/doc_examples/test_jitclass.py
   :language: python
   :start-after: magictoken.ex_jitclass_type_hints.begin
   :end-before: magictoken.ex_jitclass_type_hints.end
   :dedent: 8

Any type annotations on the class will be used to extend the spec if that field
is not already present.  The Numba type corresponding to the given Python type
is inferred using ``as_numba_type``.  For example, if we have the class

.. code-block:: python

    @jitclass([("w", int32), ("y", float64[:])])
    class Foo:
        w: int
        x: float
        y: np.ndarray
        z: SomeOtherType

        def __init__(self, w: int, x: float, y: np.ndarray, z: SomeOtherType):
            ...

then the full spec used for ``Foo`` will be:

* ``"w": int32`` (specified in the ``spec``)
* ``"x": float64`` (added from type annotation)
* ``"y": array(float64, 1d, A)`` (specified in the ``spec``)
* ``"z": numba.as_numba_type(SomeOtherType)`` (added from type annotation)

Here ``SomeOtherType`` could be any supported Python type (e.g.
``bool``, ``typing.Dict[int, typing.Tuple[float, float]]``, or another
``jitclass``).

Note that only type annotations on the class will be used to infer spec
elements.  Method type annotations (e.g. those of ``__init__`` above) are
ignored.

Numba requires knowing the dtype and rank of numpy arrays, which cannot
currently be expressed with type annotations. Because of this, numpy arrays need
to be included in the ``spec`` explicitly.


Specifying ``numba.typed`` containers as class members explicitly
=================================================================

The following patterns demonstrate how to specify a ``numba.typed.Dict`` or
``numba.typed.List`` explicitly as part of the ``spec`` passed to ``jitclass``.

First, using explicit Numba types and explicit construction.

.. code-block:: python

    from numba import jitclass, types, typed

    # key and value types
    kv_ty = (types.int64, types.unicode_type)

    # A container class with:
    # * member 'd' holding a typed dictionary of int64 -> unicode string (kv_ty)
    # * member 'l' holding a typed list of float64
    @jitclass([('d', types.DictType(*kv_ty)),
               ('l', types.ListType(types.float64))])
    class ContainerHolder(object):
        def __init__(self):
            # initialize the containers
            self.d = typed.Dict.empty(*kv_ty)
            self.l = typed.List.empty_list(types.float64)

    container = ContainerHolder()
    container.d[1] = "apple"
    container.d[2] = "orange"
    container.l.append(123.)
    container.l.append(456.)
    print(container.d) # {1: apple, 2: orange}
    print(container.l) # [123.0, 456.0]

Another useful pattern is to use the ``numba.typed`` container attribute
``_numba_type_`` to find the type of a container, this can be accessed directly
from an instance of the container in the Python interpreter. The same
information can be obtained by calling :func:`numba.typeof` on the instance. For
example:

.. code-block:: python

    from numba import jitclass, typed, typeof

    d = typed.Dict()
    d[1] = "apple"
    d[2] = "orange"
    l = typed.List()
    l.append(123.)
    l.append(456.)


    @jitclass([('d', typeof(d)), ('l', typeof(l))])
    class ContainerInstHolder(object):
        def __init__(self, dict_inst, list_inst):
            self.d = dict_inst
            self.l = list_inst

    container = ContainerInstHolder(d, l)
    print(container.d) # {1: apple, 2: orange}
    print(container.l) # [123.0, 456.0]

It is worth noting that the instance of the container in a ``jitclass`` must be
initialized before use, for example, this will cause an invalid memory access
as ``self.d`` is written to without ``d`` being initialized as a ``type.Dict``
instance of the type specified.

.. code-block:: python

    from numba import jitclass, types

    dict_ty = types.DictType(types.int64, types.unicode_type)

    @jitclass([('d', dict_ty)])
    class NotInitialisingContainer(object):
        def __init__(self):
            self.d[10] = "apple" # this is invalid, `d` is not initialized

    NotInitialisingContainer() # segmentation fault/memory access violation


Support operations
==================

The following operations of jitclasses work in both the interpreter and Numba
compiled functions:

* calling the jitclass class object to construct a new instance
  (e.g. ``mybag = Bag(123)``);
* read/write access to attributes and properties (e.g. ``mybag.value``);
* calling methods (e.g. ``mybag.increment(3)``);
* calling static methods as instance attributes (e.g. ``mybag.add(1, 1)``);
* calling static methods as class attributes (e.g. ``Bag.add(1, 2)``);

Using jitclasses in Numba compiled function is more efficient.
Short methods can be inlined (at the discretion of LLVM inliner).
Attributes access are simply reading from a C structure.
Using jitclasses from the interpreter has the same overhead of calling any
Numba compiled function from the interpreter.  Arguments and return values
must be unboxed or boxed between Python objects and native representation.
Values encapsulated by a jitclass does not get boxed into Python object when
the jitclass instance is handed to the interpreter.  It is during attribute
access to the field values that they are boxed.
Calling static methods as class attributes is only supported outside of the
class definition (i.e. code cannot call ``Bag.add()`` from within another method
of ``Bag``).


Limitations
===========

* A jitclass class object is treated as a function (the constructor) inside
  a Numba compiled function.
* ``isinstance()`` only works in the interpreter.
* Manipulating jitclass instances in the interpreter is not optimized, yet.
* Support for jitclasses are available on CPU only.
  (Note: Support for GPU devices is planned for a future release.)


The decorator: ``@jitclass``
============================

.. autofunction:: numba.experimental.jitclass
.. _generated-jit:

================================================
Flexible specializations with ``@generated_jit``
================================================


While the :func:`~numba.jit` decorator is useful for many situations,
sometimes you want to write a function that has different implementations
depending on its input types.  The :func:`~numba.generated_jit` decorator
allows the user to control the selection of a specialization at compile-time,
while fully retaining runtime execution speed of a JIT function.


Example
=======

Suppose you want to write a function which returns whether a given value
is a "missing" value according to certain conventions.  For the sake of
the example, let's adopt the following definition:

- for floating-point arguments, a missing value is a ``NaN``
- for Numpy datetime64 and timedelta64 arguments, a missing value is a ``NaT``
- other types don't have the concept of a missing value.

That compile-time logic is easily implemented using the
:func:`~numba.generated_jit` decorator::

   import numpy as np

   from numba import generated_jit, types

   @generated_jit(nopython=True)
   def is_missing(x):
       """
       Return True if the value is missing, False otherwise.
       """
       if isinstance(x, types.Float):
           return lambda x: np.isnan(x)
       elif isinstance(x, (types.NPDatetime, types.NPTimedelta)):
           # The corresponding Not-a-Time value
           missing = x('NaT')
           return lambda x: x == missing
       else:
           return lambda x: False


There are several things to note here:

* The decorated function is called with the :ref:`Numba types <numba-types>`
  of the arguments, not their values.

* The decorated function doesn't actually compute a result, it returns
  a callable implementing the actual definition of the function for the
  given types.

* It is possible to pre-compute some data at compile-time (the ``missing``
  variable above) to have them reused inside the compiled implementation.

* The function definitions use the same names for arguments as in the
  decorated function, this is required to ensure passing arguments by
  name works as expected.


Compilation options
===================

The :func:`~numba.generated_jit` decorator supports the same keyword-only
arguments as the :func:`~numba.jit` decorator, for example the ``nopython``
and ``cache`` options.


==========================
Frequently Asked Questions
==========================

Installation
============

Numba could not be imported
---------------------------

If you are seeing an exception on importing Numba with an error message
that starts with::

    ImportError: Numba could not be imported.

here are some common issues and things to try to fix it.

#. Your installation has more than one version of Numba a given environment.

   Common ways this occurs include:

   * Installing Numba with conda and then installing again with pip.
   * Installing Numba with pip and then updating to a new version with pip (pip
     re-installations don't seem to always clean up very well).

   To fix this the best approach is to create an entirely new environment and
   install a single version of Numba in that environment using a package manager
   of your choice.

#. Your installation has Numba for Python version X but you are running with
   Python version Y.

   This occurs due to a variety of Python environment mix-up/mismatch problems.
   The most common mismatch comes from installing Numba into the
   site-packages/environment of one version of Python by using a base or
   system installation of Python that is a different version, this typically
   happens through the use of the "wrong" ``pip`` binary. This will obviously
   cause problems as the C-Extensions on which Numba relies are bound to
   specific Python versions. A way to check if this likely the problem is to
   see if the path to the ``python`` binary at::

       python -c 'import sys; print(sys.executable)'

   matches the path to your installation tool and/or matches the reported
   installation location and if the Python versions match up across all of
   these. Note that Python version ``X.Y.A`` is compatible with ``X.Y.B``.

   To fix this the best approach is to create an entirely new environment and
   ensure that the installation tool used to install Numba is the one from that
   environment/the Python versions at install and run time match.

#. Your core system libraries are too old.

   This is a somewhat rare occurrence, but there are occasions when a very old
   (typically out of support) version of Linux is in use it doesn't have a
   ``glibc`` library with sufficiently new versioned symbols for Numba's shared
   libraries to resolve against. The fix for this is to update your OS system
   libraries/update your OS.

#. You are using an IDE e.g. Spyder.

   There are some unknown issues in relation to installing Numba via IDEs, but
   it would appear that these are likely variations of 1. or 2. with the same
   suggested fixes. Also, try installation from outside of the IDE with the
   command line.


If you have an installation problem which is not one of the above problems,
please do ask on `numba.discourse.group <https://numba.discourse.group/>`_ and
if possible include the path where Numba is installed and also the output of::

    python -c 'import sys; print(sys.executable)'


Programming
===========

Can I pass a function as an argument to a jitted function?
----------------------------------------------------------

As of Numba 0.39, you can, so long as the function argument has also been
JIT-compiled::

   @jit(nopython=True)
   def f(g, x):
       return g(x) + g(-x)

   result = f(jitted_g_function, 1)

However, dispatching with arguments that are functions has extra overhead.
If this matters for your application, you can also use a factory function to 
capture the function argument in a closure::

   def make_f(g):
       # Note: a new f() is created each time make_f() is called!
       @jit(nopython=True)
       def f(x):
           return g(x) + g(-x)
       return f

   f = make_f(jitted_g_function)
   result = f(1)

Improving the dispatch performance of functions in Numba is an ongoing task.

Numba doesn't seem to care when I modify a global variable
----------------------------------------------------------

Numba considers global variables as compile-time constants.  If you want
your jitted function to update itself when you have modified a global
variable's value, one solution is to recompile it using the
:meth:`~Dispatcher.recompile` method.  This is a relatively slow operation,
though, so you may instead decide to rearchitect your code and turn the
global variable into a function argument.

Can I debug a jitted function?
------------------------------

Calling into :mod:`pdb` or other such high-level facilities is currently not
supported from Numba-compiled code.  However, you can temporarily disable
compilation by setting the :envvar:`NUMBA_DISABLE_JIT` environment
variable.

How can I create a Fortran-ordered array?
-----------------------------------------

Numba currently doesn't support the ``order`` argument to most Numpy
functions such as :func:`numpy.empty` (because of limitations in the
:term:`type inference` algorithm).  You can work around this issue by
creating a C-ordered array and then transposing it.  For example::

   a = np.empty((3, 5), order='F')
   b = np.zeros(some_shape, order='F')

can be rewritten as::

   a = np.empty((5, 3)).T
   b = np.zeros(some_shape[::-1]).T

How can I increase integer width?
---------------------------------

By default, Numba will generally use machine integer width for integer
variables.  On a 32-bit machine, you may sometimes need the magnitude of
64-bit integers instead.  You can simply initialize relevant variables as
``np.int64`` (for example ``np.int64(0)`` instead of ``0``).  It will
propagate to all computations involving those variables.

.. _parallel_faqs:

How can I tell if ``parallel=True`` worked?
-------------------------------------------

If the ``parallel=True`` transformations failed for a function
decorated as such, a warning will be displayed. See also
:ref:`numba-parallel-diagnostics` for information about parallel diagnostics.

Performance
===========

Does Numba inline functions?
----------------------------

Numba gives enough information to LLVM so that functions short enough
can be inlined.  This only works in :term:`nopython mode`.

Does Numba vectorize array computations (SIMD)?
-----------------------------------------------

Numba doesn't implement such optimizations by itself, but it lets LLVM
apply them.

Why has my loop not vectorized?
-------------------------------

Numba enables the loop-vectorize optimization in LLVM by default.
While it is a powerful optimization, not all loops are applicable.
Sometimes, loop-vectorization may fail due to subtle details like memory access
pattern. To see additional diagnostic information from LLVM,
add the following lines:

.. code-block:: python

    import llvmlite.binding as llvm
    llvm.set_option('', '--debug-only=loop-vectorize')

This tells LLVM to print debug information from the **loop-vectorize**
pass to stderr.  Each function entry looks like:

.. code-block:: text

    LV: Checking a loop in "<low-level symbol name>" from <function name>
    LV: Loop hints: force=? width=0 unroll=0
    ...
    LV: Vectorization is possible but not beneficial.
    LV: Interleaving is not beneficial.

Each function entry is separated by an empty line.  The reason for rejecting
the vectorization is usually at the end of the entry.  In the example above,
LLVM rejected the vectorization because doing so will not speedup the loop.
In this case, it can be due to memory access pattern.  For instance, the
array being looped over may not be in contiguous layout.

When memory access pattern is non-trivial such that it cannot determine the
access memory region, LLVM may reject with the following message:

.. code-block:: text

    LV: Can't vectorize due to memory conflicts

Another common reason is:

.. code-block:: text

    LV: Not vectorizing: loop did not meet vectorization requirements.

In this case, vectorization is rejected because the vectorized code may behave
differently.  This is a case to try turning on ``fastmath=True`` to allow
fastmath instructions.

Why are the ``typed`` containers slower when used from the interpreter?
-----------------------------------------------------------------------

The Numba ``typed`` containers found in ``numba.typed`` e.g.
``numba.typed.List`` store their data in an efficient form for access from JIT
compiled code. When these containers are used from the CPython interpreter, the
data involved has to be converted from/to the container format. This process is
relatively costly and as a result impacts performance. In JIT compiled code no
such penalty exists and so operations on the containers are much quicker and
often faster than the pure Python equivalent.

Does Numba automatically parallelize code?
------------------------------------------

It can, in some cases:

* Ufuncs and gufuncs with the ``target="parallel"`` option will run on multiple threads.
* The ``parallel=True`` option to ``@jit`` will attempt to optimize array
  operations and run them in parallel.  It also adds support for ``prange()`` to
  explicitly parallelize a loop.

You can also manually run computations on multiple threads yourself and use
the ``nogil=True`` option (see :ref:`releasing the GIL <jit-nogil>`).  Numba
can also target parallel execution on GPU architectures using its CUDA and HSA
backends.


Can Numba speed up short-running functions?
-------------------------------------------

Not significantly.  New users sometimes expect to JIT-compile such
functions::

   def f(x, y):
       return x + y

and get a significant speedup over the Python interpreter.  But there isn't
much Numba can improve here: most of the time is probably spent in CPython's
function call mechanism, rather than the function itself.  As a rule of
thumb, if a function takes less than 10 µs to execute: leave it.

The exception is that you *should* JIT-compile that function if it is called
from another jitted function.

There is a delay when JIT-compiling a complicated function, how can I improve it?
---------------------------------------------------------------------------------

Try to pass ``cache=True`` to the ``@jit`` decorator.  It will keep the
compiled version on disk for later use.

A more radical alternative is :ref:`ahead-of-time compilation <pycc>`.


GPU Programming
===============

How do I work around the ``CUDA intialized before forking`` error?
------------------------------------------------------------------

On Linux, the ``multiprocessing`` module in the Python standard library
defaults to using the ``fork`` method for creating new processes.  Because of
the way process forking duplicates state between the parent and child
processes, CUDA will not work correctly in the child process if the CUDA
runtime was initialized *prior* to the fork.  Numba detects this and raises a
``CudaDriverError`` with the message ``CUDA initialized before forking``.

One approach to avoid this error is to make all calls to ``numba.cuda``
functions inside the child processes or after the process pool is created.
However, this is not always possible, as you might want to query the number of
available GPUs before starting the process pool.  In Python 3, you can change
the process start method, as described in the `multiprocessing documentation
<https://docs.python.org/3.9/library/multiprocessing.html#contexts-and-start-methods>`_.
Switching from ``fork`` to ``spawn`` or ``forkserver`` will avoid the CUDA
initalization issue, although the child processes will not inherit any global
variables from their parent.


Integration with other utilities
================================

Can I "freeze" an application which uses Numba?
-----------------------------------------------

If you're using PyInstaller or a similar utility to freeze an application,
you may encounter issues with llvmlite.  llvmlite needs a non-Python DLL
for its working, but it won't be automatically detected by freezing utilities.
You have to inform the freezing utility of the DLL's location: it will
usually be named ``llvmlite/binding/libllvmlite.so`` or
``llvmlite/binding/llvmlite.dll``, depending on your system.

I get errors when running a script twice under Spyder
-----------------------------------------------------

When you run a script in a console under Spyder, Spyder first tries to
reload existing modules.  This doesn't work well with Numba, and can
produce errors like ``TypeError: No matching definition for argument type(s)``.

There is a fix in the Spyder preferences. Open the "Preferences" window,
select "Console", then "Advanced Settings", click the "Set UMR excluded
modules" button, and add ``numba`` inside the text box that pops up.

To see the setting take effect, be sure to restart the IPython console or
kernel.

.. _llvm-locale-bug:

Why does Numba complain about the current locale?
-------------------------------------------------

If you get an error message such as the following::

   RuntimeError: Failed at nopython (nopython mode backend)
   LLVM will produce incorrect floating-point code in the current locale

it means you have hit a LLVM bug which causes incorrect handling of
floating-point constants.  This is known to happen with certain third-party
libraries such as the Qt backend to matplotlib.

To work around the bug, you need to force back the locale to its default
value, for example::

   import locale
   locale.setlocale(locale.LC_NUMERIC, 'C')

How do I get Numba development builds?
--------------------------------------

Pre-release versions of Numba can be installed with conda::

    $ conda install -c numba/label/dev numba


Miscellaneous
=============

Where does the project name "Numba" come from?
----------------------------------------------

"Numba" is a combination of "NumPy" and "Mamba". Mambas are some of the fastest
snakes in the world, and Numba makes your Python code fast.

How do I reference/cite/acknowledge Numba in other work?
--------------------------------------------------------
For academic use, the best option is to cite our ACM Proceedings: `Numba: a
LLVM-based Python JIT compiler.
<http://dl.acm.org/citation.cfm?id=2833162&dl=ACM&coll=DL>`_ You can also find
`the sources on github <https://github.com/numba/Numba-SC15-Paper>`_, including
`a pre-print pdf
<https://github.com/numba/Numba-SC15-Paper/raw/master/numba_sc15.pdf>`_, in case
you don't have access to the ACM site but would like to read the paper.

Other related papers
~~~~~~~~~~~~~~~~~~~~
A paper describing ParallelAccelerator technology, that is activated when the
``parallel=True`` jit option is used, can be found `here
<http://drops.dagstuhl.de/opus/volltexte/2017/7269/pdf/LIPIcs-ECOOP-2017-4.pdf>`_.

How do I write a minimal working reproducer for a problem with Numba?
---------------------------------------------------------------------

A minimal working reproducer for Numba should include:

1. The source code of the function(s) that reproduce the problem.
2. Some example data and a demonstration of calling the reproducing code with
   that data. As Numba compiles based on type information, unless your problem
   is numerical, it's fine to just provide dummy data of the right type, e.g.
   use ``numpy.ones`` of the correct ``dtype``/size/shape for arrays.
3. Ideally put 1. and 2. into a script with all the correct imports. Make sure
   your script actually executes and reproduces the problem before submitting
   it! The target is to make it so that the script can just be copied directly
   from the `issue tracker <https://github.com/numba/numba/issues>`_ and run by
   someone else such that they can see the same problem as you are having.

Having made a reproducer, now remove every part of the code that does not
contribute directly to reproducing the problem to create a "minimal" reproducer.
This means removing imports that aren't used, removing variables that aren't
used or have no effect, removing lines of code which have no effect, reducing
the complexity of expressions, and shrinking input data to the minimal amount
required to trigger the problem.

Doing the above really helps out the Numba issue triage process and will enable
a faster response to your problem!

`Suggested further reading
<http://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports>`_ on
writing minimal working reproducers.

Overview
========

Numba is a compiler for Python array and numerical functions that gives 
you the power to speed up your applications with high performance
functions written directly in Python.

Numba generates optimized machine code from pure Python code using
the `LLVM compiler infrastructure <http://llvm.org/>`_.  With a few simple
annotations, array-oriented and math-heavy Python code can be
just-in-time optimized to performance similar as C, C++ and Fortran, without
having to switch languages or Python interpreters.

Numba's main features are:

* :ref:`on-the-fly code generation <jit>` (at import time or runtime, at the
  user's preference)
* native code generation for the CPU (default) and
  :doc:`GPU hardware <../cuda/index>`
* integration with the Python scientific software stack (thanks to Numpy)

Here is how a Numba-optimized function, taking a Numpy array as argument,
might look like::

   @numba.jit
   def sum2d(arr):
       M, N = arr.shape
       result = 0.0
       for i in range(M):
           for j in range(N):
               result += arr[i,j]
       return result

============================================================
Callback into the Python Interpreter from within JIT'ed code
============================================================

There are rare but real cases when a nopython-mode function needs to callback
into the Python interpreter to invoke code that cannot be compiled by Numba.
Such cases include:

- logging progress for long running JIT'ed functions;
- use data structures that are not currently supported by Numba;
- debugging inside JIT'ed code using the Python debugger.

When Numba callbacks into the Python interpreter, the following has to happen:

- acquire the GIL;
- convert values in native representation back into Python objects;
- call-back into the Python interpreter;
- convert returned values from the Python-code into native representation;
- release the GIL.

These steps can be expensive.  Users **should not** rely on the feature
described here on performance-critical paths.


.. _with_objmode:

The ``objmode`` context-manager
===============================

.. warning:: This feature can be easily mis-used.  Users should first consider
    alternative approaches to achieve their intended goal before using
    this feature.

.. autofunction:: numba.objmode

User Manual
===========

.. toctree::

   5minguide.rst
   overview.rst
   installing.rst
   jit.rst
   generated-jit.rst
   vectorize.rst
   jitclass.rst
   cfunc.rst
   pycc.rst
   parallel.rst
   stencil.rst
   withobjmode.rst
   jit-module.rst
   performance-tips.rst
   threading-layer.rst
   cli.rst
   troubleshoot.rst
   faq.rst
   examples.rst
   talks.rst
.. Copyright (c) 2017 Intel Corporation
   SPDX-License-Identifier: BSD-2-Clause

.. _numba-parallel:

=======================================
Automatic parallelization with ``@jit``
=======================================

Setting the :ref:`parallel_jit_option` option for :func:`~numba.jit` enables
a Numba transformation pass that attempts to automatically parallelize and
perform other optimizations on (part of) a function. At the moment, this
feature only works on CPUs.

Some operations inside a user defined function, e.g. adding a scalar value to
an array, are known to have parallel semantics.  A user program may contain
many such operations and while each operation could be parallelized
individually, such an approach often has lackluster performance due to poor
cache behavior.  Instead, with auto-parallelization, Numba attempts to
identify such operations in a user program, and fuse adjacent ones together,
to form one or more kernels that are automatically run in parallel.
The process is fully automated without modifications to the user program,
which is in contrast to Numba's :func:`~numba.vectorize` or
:func:`~numba.guvectorize` mechanism, where manual effort is required
to create parallel kernels.

.. _numba-parallel-supported:

Supported Operations
====================

In this section, we give a list of all the array operations that have
parallel semantics and for which we attempt to parallelize.

#. All numba array operations that are supported by :ref:`case-study-array-expressions`,
   which include common arithmetic functions between Numpy arrays, and between
   arrays and scalars, as well as Numpy ufuncs. They are often called
   `element-wise` or `point-wise` array operations:

    * unary operators: ``+`` ``-`` ``~``
    * binary operators: ``+`` ``-`` ``*`` ``/`` ``/?`` ``%`` ``|`` ``>>`` ``^`` ``<<`` ``&`` ``**`` ``//``
    * comparison operators: ``==`` ``!=`` ``<`` ``<=`` ``>`` ``>=``
    * :ref:`Numpy ufuncs <supported_ufuncs>` that are supported in :term:`nopython mode`.
    * User defined :class:`~numba.DUFunc` through :func:`~numba.vectorize`.

#. Numpy reduction functions ``sum``, ``prod``, ``min``, ``max``, ``argmin``,
   and ``argmax``. Also, array math functions ``mean``, ``var``, and ``std``.

#. Numpy array creation functions ``zeros``, ``ones``, ``arange``, ``linspace``,
   and several random functions (rand, randn, ranf, random_sample, sample,
   random, standard_normal, chisquare, weibull, power, geometric, exponential,
   poisson, rayleigh, normal, uniform, beta, binomial, f, gamma, lognormal,
   laplace, randint, triangular).

#. Numpy ``dot`` function between a matrix and a vector, or two vectors.
   In all other cases, Numba's default implementation is used.

#. Multi-dimensional arrays are also supported for the above operations
   when operands have matching dimension and size. The full semantics of
   Numpy broadcast between arrays with mixed dimensionality or size is
   not supported, nor is the reduction across a selected dimension.

#. Array assignment in which the target is an array selection using a slice
   or a boolean array, and the value being assigned is either a scalar or
   another selection where the slice range or bitarray are inferred to be
   compatible.

#. The ``reduce`` operator of ``functools`` is supported for specifying parallel
   reductions on 1D Numpy arrays but the initial value argument is mandatory.

.. _numba-prange:

Explicit Parallel Loops
========================

Another feature of the code transformation pass (when ``parallel=True``) is
support for explicit parallel loops. One can use Numba's ``prange`` instead of
``range`` to specify that a loop can be parallelized. The user is required to
make sure that the loop does not have cross iteration dependencies except for
supported reductions.

A reduction is inferred automatically if a variable is updated by a binary
function/operator using its previous value in the loop body. The initial value
of the reduction is inferred automatically for the ``+=``, ``-=``,  ``*=``,
and ``/=`` operators.
For other functions/operators, the reduction variable should hold the identity
value right before entering the ``prange`` loop.  Reductions in this manner
are supported for scalars and for arrays of arbitrary dimensions.

The example below demonstrates a parallel loop with a
reduction (``A`` is a one-dimensional Numpy array)::

    from numba import njit, prange

    @njit(parallel=True)
    def prange_test(A):
        s = 0
        # Without "parallel=True" in the jit-decorator
        # the prange statement is equivalent to range
        for i in prange(A.shape[0]):
            s += A[i]
        return s

The following example demonstrates a product reduction on a two-dimensional array::

    from numba import njit, prange
    import numpy as np

    @njit(parallel=True)
    def two_d_array_reduction_prod(n):
        shp = (13, 17)
        result1 = 2 * np.ones(shp, np.int_)
        tmp = 2 * np.ones_like(result1)

        for i in prange(n):
            result1 *= tmp

        return result1

.. note:: When using Python's ``range`` to induce a loop, Numba types the
          induction variable as a signed integer. This is also the case for
          Numba's ``prange`` when ``parallel=False``. However, for
          ``parallel=True``, if the range is identifiable as strictly positive,
          the type of the induction variable  will be ``uint64``. The impact of
          a ``uint64`` induction variable is often most noticable when
          undertaking operations involving it and a signed integer. Under
          Numba's type coercion rules, such a case will commonly result in the
          operation producing a floating point result type.


Care should be taken, however, when reducing into slices or elements of an array 
if the elements specified by the slice or index are written to simultaneously by 
multiple parallel threads. The compiler may not detect such cases and then a race condition
would occur.

The following example demonstrates such a case where a race condition in the execution of the 
parallel for-loop results in an incorrect return value::

    from numba import njit, prange
    import numpy as np

    @njit(parallel=True)
    def prange_wrong_result(x):
        n = x.shape[0]
        y = np.zeros(4)
        for i in prange(n):
            # accumulating into the same element of `y` from different
            # parallel iterations of the loop results in a race condition
            y[:] += x[i]

        return y

as does the following example where the accumulating element is explicitly specified::

    from numba import njit, prange
    import numpy as np

    @njit(parallel=True)
    def prange_wrong_result(x):
        n = x.shape[0]
        y = np.zeros(4)
        for i in prange(n):
            # accumulating into the same element of `y` from different
            # parallel iterations of the loop results in a race condition
            y[i % 4] += x[i]

        return y

whereas performing a whole array reduction is fine::

   from numba import njit, prange
   import numpy as np

   @njit(parallel=True)
   def prange_ok_result_whole_arr(x):
       n = x.shape[0]
       y = np.zeros(4)
       for i in prange(n):
           y += x[i]
       return y

as is creating a slice reference outside of the parallel reduction loop::

   from numba import njit, prange
   import numpy as np

   @njit(parallel=True)
   def prange_ok_result_outer_slice(x):
       n = x.shape[0]
       y = np.zeros(4)
       z = y[:]
       for i in prange(n):
           z += x[i]
       return y

Examples
========

In this section, we give an example of how this feature helps
parallelize Logistic Regression::

    @numba.jit(nopython=True, parallel=True)
    def logistic_regression(Y, X, w, iterations):
        for i in range(iterations):
            w -= np.dot(((1.0 / (1.0 + np.exp(-Y * np.dot(X, w))) - 1.0) * Y), X)
        return w

We will not discuss details of the algorithm, but instead focus on how
this program behaves with auto-parallelization:

1. Input ``Y`` is a vector of size ``N``, ``X`` is an ``N x D`` matrix,
   and ``w`` is a vector of size ``D``.

2. The function body is an iterative loop that updates variable ``w``.
   The loop body consists of a sequence of vector and matrix operations.

3. The inner ``dot`` operation produces a vector of size ``N``, followed by a
   sequence of arithmetic operations either between a scalar and vector of
   size ``N``, or two vectors both of size ``N``.

4. The outer ``dot`` produces a vector of size ``D``, followed by an inplace
   array subtraction on variable ``w``.

5. With auto-parallelization, all operations that produce array of size
   ``N`` are fused together to become a single parallel kernel. This includes
   the inner ``dot`` operation and all point-wise array operations following it.

6. The outer ``dot`` operation produces a result array of different dimension,
   and is not fused with the above kernel.

Here, the only thing required to take advantage of parallel hardware is to set
the :ref:`parallel_jit_option` option for :func:`~numba.jit`, with no
modifications to the ``logistic_regression`` function itself.  If we were to
give an equivalence parallel implementation using :func:`~numba.guvectorize`,
it would require a pervasive change that rewrites the code to extract kernel
computation that can be parallelized, which was both tedious and challenging.

Unsupported Operations
======================

This section contains a non-exhaustive list of commonly encountered but
currently unsupported features:

#. **Mutating a list is not threadsafe**

   Concurrent write operations on container types (i.e. lists, sets and
   dictionaries) in a ``prange`` parallel region are not threadsafe e.g.::

    @njit(parallel=True)
    def invalid():
        z = []
        for i in prange(10000):
            z.append(i)
        return z

   It is highly likely that the above will result in corruption or an access
   violation as containers require thread-safety under mutation but this feature
   is not implemented.

#. **Induction variables are not associated with thread ID**

   The use of the induction variable induced by a ``prange`` based loop in
   conjunction with ``get_num_threads`` as a method of ensuring safe writes into
   a pre-sized container is not valid e.g.::

    @njit(parallel=True)
    def invalid():
        n = get_num_threads()
        z = [0 for _ in range(n)]
        for i in prange(100):
            z[i % n] += i
        return z

   The above can on occasion appear to work, but it does so by luck. There's no
   guarantee about which indexes are assigned to which executing threads or the
   order in which the loop iterations execute.

.. _numba-parallel-diagnostics:

Diagnostics
===========

.. note:: At present not all parallel transforms and functions can be tracked
          through the code generation process. Occasionally diagnostics about
          some loops or transforms may be missing.

The :ref:`parallel_jit_option` option for :func:`~numba.jit` can produce
diagnostic information about the transforms undertaken in automatically
parallelizing the decorated code. This information can be accessed in two ways,
the first is by setting the environment variable
:envvar:`NUMBA_PARALLEL_DIAGNOSTICS`, the second is by calling
:meth:`~Dispatcher.parallel_diagnostics`, both methods give the same information
and print to ``STDOUT``. The level of verbosity in the diagnostic information is
controlled by an integer argument of value between 1 and 4 inclusive, 1 being
the least verbose and 4 the most. For example::

    @njit(parallel=True)
    def test(x):
        n = x.shape[0]
        a = np.sin(x)
        b = np.cos(a * a)
        acc = 0
        for i in prange(n - 2):
            for j in prange(n - 1):
                acc += b[i] + b[j + 1]
        return acc

    test(np.arange(10))

    test.parallel_diagnostics(level=4)

produces::

    ================================================================================
    ======= Parallel Accelerator Optimizing:  Function test, example.py (4)  =======
    ================================================================================


    Parallel loop listing for  Function test, example.py (4)
    --------------------------------------|loop #ID
    @njit(parallel=True)                  |
    def test(x):                          |
        n = x.shape[0]                    |
        a = np.sin(x)---------------------| #0
        b = np.cos(a * a)-----------------| #1
        acc = 0                           |
        for i in prange(n - 2):-----------| #3
            for j in prange(n - 1):-------| #2
                acc += b[i] + b[j + 1]    |
        return acc                        |
    --------------------------------- Fusing loops ---------------------------------
    Attempting fusion of parallel loops (combines loops with similar properties)...
    Trying to fuse loops #0 and #1:
        - fusion succeeded: parallel for-loop #1 is fused into for-loop #0.
    Trying to fuse loops #0 and #3:
        - fusion failed: loop dimension mismatched in axis 0. slice(0, x_size0.1, 1)
    != slice(0, $40.4, 1)
    ----------------------------- Before Optimization ------------------------------
    Parallel region 0:
    +--0 (parallel)
    +--1 (parallel)


    Parallel region 1:
    +--3 (parallel)
    +--2 (parallel)


    --------------------------------------------------------------------------------
    ------------------------------ After Optimization ------------------------------
    Parallel region 0:
    +--0 (parallel, fused with loop(s): 1)


    Parallel region 1:
    +--3 (parallel)
    +--2 (serial)



    Parallel region 0 (loop #0) had 1 loop(s) fused.

    Parallel region 1 (loop #3) had 0 loop(s) fused and 1 loop(s) serialized as part
    of the larger parallel loop (#3).
    --------------------------------------------------------------------------------
    --------------------------------------------------------------------------------

    ---------------------------Loop invariant code motion---------------------------

    Instruction hoisting:
    loop #0:
    Failed to hoist the following:
        dependency: $arg_out_var.10 = getitem(value=x, index=$parfor__index_5.99)
        dependency: $0.6.11 = getattr(value=$0.5, attr=sin)
        dependency: $expr_out_var.9 = call $0.6.11($arg_out_var.10, func=$0.6.11, args=[Var($arg_out_var.10, example.py (7))], kws=(), vararg=None)
        dependency: $arg_out_var.17 = $expr_out_var.9 * $expr_out_var.9
        dependency: $0.10.20 = getattr(value=$0.9, attr=cos)
        dependency: $expr_out_var.16 = call $0.10.20($arg_out_var.17, func=$0.10.20, args=[Var($arg_out_var.17, example.py (8))], kws=(), vararg=None)
    loop #3:
    Has the following hoisted:
        $const58.3 = const(int, 1)
        $58.4 = _n_23 - $const58.3
    --------------------------------------------------------------------------------



To aid users unfamiliar with the transforms undertaken when the
:ref:`parallel_jit_option` option is used, and to assist in the understanding of
the subsequent sections, the following definitions are provided:

* Loop fusion
    `Loop fusion <https://en.wikipedia.org/wiki/Loop_fission_and_fusion>`_ is a
    technique whereby loops with equivalent bounds may be combined under certain
    conditions to produce a loop with a larger body (aiming to improve data
    locality).

* Loop serialization
    Loop serialization occurs when any number of ``prange`` driven loops are
    present inside another ``prange`` driven loop. In this case the outermost
    of all the ``prange`` loops executes in parallel and any inner ``prange``
    loops (nested or otherwise) are treated as standard ``range`` based loops.
    Essentially, nested parallelism does not occur.

* Loop invariant code motion
    `Loop invariant code motion
    <https://en.wikipedia.org/wiki/Loop-invariant_code_motion>`_ is an
    optimization technique that analyses a loop to look for statements that can
    be moved outside the loop body without changing the result of executing the
    loop, these statements are then "hoisted" out of the loop to save repeated
    computation.

* Allocation hoisting
    Allocation hoisting is a specialized case of loop invariant code motion that
    is possible due to the design of some common NumPy allocation methods.
    Explanation of this technique is best driven by an example:

    .. code-block:: python

        @njit(parallel=True)
        def test(n):
            for i in prange(n):
                temp = np.zeros((50, 50)) # <--- Allocate a temporary array with np.zeros()
                for j in range(50):
                    temp[j, j] = i

            # ...do something with temp

    internally, this is transformed to approximately the following:

    .. code-block:: python

        @njit(parallel=True)
        def test(n):
            for i in prange(n):
                temp = np.empty((50, 50)) # <--- np.zeros() is rewritten as np.empty()
                temp[:] = 0               # <--- and then a zero initialisation
                for j in range(50):
                    temp[j, j] = i

            # ...do something with temp

    then after hoisting:

    .. code-block:: python

        @njit(parallel=True)
        def test(n):
            temp = np.empty((50, 50)) # <--- allocation is hoisted as a loop invariant as `np.empty` is considered pure
            for i in prange(n):
                temp[:] = 0           # <--- this remains as assignment is a side effect
                for j in range(50):
                    temp[j, j] = i

            # ...do something with temp

    it can be seen that the ``np.zeros`` allocation is split into an allocation
    and an assignment, and then the allocation is hoisted out of the loop in
    ``i``, this producing more efficient code as the allocation only occurs
    once.

The parallel diagnostics report sections
----------------------------------------

The report is split into the following sections:

#. Code annotation
    This is the first section and contains the source code of the decorated
    function with loops that have parallel semantics identified and enumerated.
    The ``loop #ID`` column on the right of the source code lines up with
    identified parallel loops. From the example, ``#0`` is ``np.sin``, ``#1``
    is ``np.cos`` and ``#2`` and ``#3`` are ``prange()``:

    .. code-block:: python

        Parallel loop listing for  Function test, example.py (4)
        --------------------------------------|loop #ID
        @njit(parallel=True)                  |
        def test(x):                          |
            n = x.shape[0]                    |
            a = np.sin(x)---------------------| #0
            b = np.cos(a * a)-----------------| #1
            acc = 0                           |
            for i in prange(n - 2):-----------| #3
                for j in prange(n - 1):-------| #2
                    acc += b[i] + b[j + 1]    |
            return acc                        |

    It is worth noting that the loop IDs are enumerated in the order they are
    discovered which is not necessarily the same order as present in the source.
    Further, it should also be noted that the parallel transforms use a static
    counter for loop ID indexing. As a consequence it is possible for the loop
    ID index to not start at 0 due to use of the same counter for internal
    optimizations/transforms taking place that are invisible to the user.

#. Fusing loops
    This section describes the attempts made at fusing discovered
    loops noting which succeeded and which failed. In the case of failure to
    fuse a reason is given (e.g. dependency on other data). From the example:

    .. code-block:: text

        --------------------------------- Fusing loops ---------------------------------
        Attempting fusion of parallel loops (combines loops with similar properties)...
        Trying to fuse loops #0 and #1:
            - fusion succeeded: parallel for-loop #1 is fused into for-loop #0.
        Trying to fuse loops #0 and #3:
            - fusion failed: loop dimension mismatched in axis 0. slice(0, x_size0.1, 1)
        != slice(0, $40.4, 1)

    It can be seen that fusion of loops ``#0`` and ``#1`` was attempted and this
    succeeded (both are based on the same dimensions of ``x``). Following the
    successful fusion of ``#0`` and ``#1``, fusion was attempted between ``#0``
    (now including the fused ``#1`` loop) and ``#3``. This fusion failed because
    there is a loop dimension mismatch, ``#0`` is size ``x.shape`` whereas
    ``#3`` is size ``x.shape[0] - 2``.

#. Before Optimization
    This section shows the structure of the parallel regions in the code before
    any optimization has taken place, but with loops associated with their final
    parallel region (this is to make before/after optimization output directly
    comparable). Multiple parallel regions may exist if there are loops which
    cannot be fused, in this case code within each region will execute in
    parallel, but each parallel region will run sequentially. From the example:

    .. code-block:: text

        Parallel region 0:
        +--0 (parallel)
        +--1 (parallel)


        Parallel region 1:
        +--3 (parallel)
        +--2 (parallel)

    As alluded to by the `Fusing loops` section, there are necessarily two
    parallel regions in the code. The first contains loops ``#0`` and ``#1``,
    the second contains ``#3`` and ``#2``, all loops are marked ``parallel`` as
    no optimization has taken place yet.

#. After Optimization
    This section shows the structure of the parallel regions in the code after
    optimization has taken place. Again, parallel regions are enumerated with
    their corresponding loops but this time loops which are fused or serialized
    are noted and a summary is presented. From the example:

    .. code-block:: text

        Parallel region 0:
        +--0 (parallel, fused with loop(s): 1)


        Parallel region 1:
        +--3 (parallel)
           +--2 (serial)

        Parallel region 0 (loop #0) had 1 loop(s) fused.

        Parallel region 1 (loop #3) had 0 loop(s) fused and 1 loop(s) serialized as part
        of the larger parallel loop (#3).


    It can be noted that parallel region 0 contains loop ``#0`` and, as seen in
    the `fusing loops` section, loop ``#1`` is fused into loop ``#0``. It can
    also be noted that parallel region 1 contains loop ``#3`` and that loop
    ``#2`` (the inner ``prange()``) has been serialized for execution in the
    body of loop ``#3``.

#. Loop invariant code motion
    This section shows for each loop, after optimization has occurred:

    * the instructions that failed to be hoisted and the reason for failure
      (dependency/impure).
    * the instructions that were hoisted.
    * any allocation hoisting that may have occurred.

    From the example:

    .. code-block:: text

        Instruction hoisting:
        loop #0:
        Failed to hoist the following:
            dependency: $arg_out_var.10 = getitem(value=x, index=$parfor__index_5.99)
            dependency: $0.6.11 = getattr(value=$0.5, attr=sin)
            dependency: $expr_out_var.9 = call $0.6.11($arg_out_var.10, func=$0.6.11, args=[Var($arg_out_var.10, example.py (7))], kws=(), vararg=None)
            dependency: $arg_out_var.17 = $expr_out_var.9 * $expr_out_var.9
            dependency: $0.10.20 = getattr(value=$0.9, attr=cos)
            dependency: $expr_out_var.16 = call $0.10.20($arg_out_var.17, func=$0.10.20, args=[Var($arg_out_var.17, example.py (8))], kws=(), vararg=None)
        loop #3:
        Has the following hoisted:
            $const58.3 = const(int, 1)
            $58.4 = _n_23 - $const58.3

    The first thing to note is that this information is for advanced users as it
    refers to the :term:`Numba IR` of the function being transformed. As an
    example, the expression ``a * a`` in the example source partly translates to
    the expression ``$arg_out_var.17 = $expr_out_var.9 * $expr_out_var.9`` in
    the IR, this clearly cannot be hoisted out of ``loop #0`` because it is not
    loop invariant! Whereas in ``loop #3``, the expression
    ``$const58.3 = const(int, 1)`` comes from the source ``b[j + 1]``, the
    number ``1`` is clearly a constant and so can be hoisted out of the loop.

.. seealso:: :ref:`parallel_jit_option`, :ref:`Parallel FAQs <parallel_FAQs>`
==================================
Creating NumPy universal functions
==================================

There are two types of universal functions:

* Those which operate on scalars, these are "universal functions" or *ufuncs*
  (see ``@vectorize`` below).
* Those which operate on higher dimensional arrays and scalars, these are
  "generalized universal functions" or *gufuncs* (``@guvectorize`` below).

.. _vectorize:

The ``@vectorize`` decorator
============================

Numba's vectorize allows Python functions taking scalar input arguments to
be used as NumPy `ufuncs`_.  Creating a traditional NumPy ufunc is
not the most straightforward process and involves writing some C code.
Numba makes this easy.  Using the :func:`~numba.vectorize` decorator, Numba
can compile a pure Python function into a ufunc that operates over NumPy
arrays as fast as traditional ufuncs written in C.

.. _ufuncs: http://docs.scipy.org/doc/numpy/reference/ufuncs.html

Using :func:`~numba.vectorize`, you write your function as operating over
input scalars, rather than arrays.  Numba will generate the surrounding
loop (or *kernel*) allowing efficient iteration over the actual inputs.

The :func:`~numba.vectorize` decorator has two modes of operation:

* Eager, or decoration-time, compilation: If you pass one or more type
  signatures to the decorator, you will be building a Numpy universal
  function (ufunc).  The rest of this subsection describes building
  ufuncs using decoration-time compilation.

* Lazy, or call-time, compilation: When not given any signatures, the
  decorator will give you a Numba dynamic universal function
  (:class:`~numba.DUFunc`) that dynamically compiles a new kernel when
  called with a previously unsupported input type.  A later
  subsection, ":ref:`dynamic-universal-functions`", describes this mode in
  more depth.

As described above, if you pass a list of signatures to the
:func:`~numba.vectorize` decorator, your function will be compiled
into a Numpy ufunc.  In the basic case, only one signature will be
passed::

   from numba import vectorize, float64

   @vectorize([float64(float64, float64)])
   def f(x, y):
       return x + y

If you pass several signatures, beware that you have to pass most specific
signatures before least specific ones (e.g., single-precision floats
before double-precision floats), otherwise type-based dispatching will not work
as expected::

   @vectorize([int32(int32, int32),
               int64(int64, int64),
               float32(float32, float32),
               float64(float64, float64)])
   def f(x, y):
       return x + y

The function will work as expected over the specified array types::

   >>> a = np.arange(6)
   >>> f(a, a)
   array([ 0,  2,  4,  6,  8, 10])
   >>> a = np.linspace(0, 1, 6)
   >>> f(a, a)
   array([ 0. ,  0.4,  0.8,  1.2,  1.6,  2. ])

but it will fail working on other types::

   >>> a = np.linspace(0, 1+1j, 6)
   >>> f(a, a)
   Traceback (most recent call last):
     File "<stdin>", line 1, in <module>
   TypeError: ufunc 'ufunc' not supported for the input types, and the inputs could not be safely coerced to any supported types according to the casting rule ''safe''


You might ask yourself, "why would I go through this instead of compiling
a simple iteration loop using the :ref:`@jit <jit>` decorator?".  The
answer is that NumPy ufuncs automatically get other features such as
reduction, accumulation or broadcasting.  Using the example above::

   >>> a = np.arange(12).reshape(3, 4)
   >>> a
   array([[ 0,  1,  2,  3],
          [ 4,  5,  6,  7],
          [ 8,  9, 10, 11]])
   >>> f.reduce(a, axis=0)
   array([12, 15, 18, 21])
   >>> f.reduce(a, axis=1)
   array([ 6, 22, 38])
   >>> f.accumulate(a)
   array([[ 0,  1,  2,  3],
          [ 4,  6,  8, 10],
          [12, 15, 18, 21]])
   >>> f.accumulate(a, axis=1)
   array([[ 0,  1,  3,  6],
          [ 4,  9, 15, 22],
          [ 8, 17, 27, 38]])

.. seealso::
   `Standard features of ufuncs <http://docs.scipy.org/doc/numpy/reference/ufuncs.html#ufunc>`_ (NumPy documentation).

.. note::
   Only the broadcasting features of ufuncs are supported in compiled code.

The :func:`~numba.vectorize` decorator supports multiple ufunc targets:

=================       ===============================================================
Target                    Description
=================       ===============================================================
cpu                     Single-threaded CPU


parallel                Multi-core CPU


cuda                    CUDA GPU

                        .. NOTE:: This creates an *ufunc-like* object.
			  See `documentation for CUDA ufunc <../cuda/ufunc.html>`_ for detail.
=================       ===============================================================

A general guideline is to choose different targets for different data sizes
and algorithms.
The "cpu" target works well for small data sizes (approx. less than 1KB) and low
compute intensity algorithms. It has the least amount of overhead.
The "parallel" target works well for medium data sizes (approx. less than 1MB).
Threading adds a small delay.
The "cuda" target works well for big data sizes (approx. greater than 1MB) and
high compute intensity algorithms.  Transferring memory to and from the GPU adds
significant overhead.


.. _guvectorize:

The ``@guvectorize`` decorator
==============================

While :func:`~numba.vectorize` allows you to write ufuncs that work on one
element at a time, the :func:`~numba.guvectorize` decorator takes the concept
one step further and allows you to write ufuncs that will work on an
arbitrary number of elements of input arrays, and take and return arrays of
differing dimensions.  The typical example is a running median or a
convolution filter.

Contrary to :func:`~numba.vectorize` functions, :func:`~numba.guvectorize`
functions don't return their result value: they take it as an array
argument, which must be filled in by the function.  This is because the
array is actually allocated by NumPy's dispatch mechanism, which calls into
the Numba-generated code.

Similar to :func:`~numba.vectorize` decorator, :func:`~numba.guvectorize`
also has two modes of operation: Eager, or decoration-time compilation and
lazy, or call-time compilation.


Here is a very simple example::

   @guvectorize([(int64[:], int64, int64[:])], '(n),()->(n)')
   def g(x, y, res):
       for i in range(x.shape[0]):
           res[i] = x[i] + y

The underlying Python function simply adds a given scalar (``y``) to all
elements of a 1-dimension array.  What's more interesting is the declaration.
There are two things there:

* the declaration of input and output *layouts*, in symbolic form:
  ``(n),()->(n)`` tells NumPy that the function takes a *n*-element one-dimension
  array, a scalar (symbolically denoted by the empty tuple ``()``) and
  returns a *n*-element one-dimension array;

* the list of supported concrete *signatures* as per ``@vectorize``; here,
  as in the above example, we demonstrate ``int64`` arrays.

.. note::
   1D array type can also receive scalar arguments (those with shape ``()``).
   In the above example, the second argument also could be declared as
   ``int64[:]``.  In that case, the value must be read by ``y[0]``.

We can now check what the compiled ufunc does, over a simple example::

   >>> a = np.arange(5)
   >>> a
   array([0, 1, 2, 3, 4])
   >>> g(a, 2)
   array([2, 3, 4, 5, 6])

The nice thing is that NumPy will automatically dispatch over more
complicated inputs, depending on their shapes::

   >>> a = np.arange(6).reshape(2, 3)
   >>> a
   array([[0, 1, 2],
          [3, 4, 5]])
   >>> g(a, 10)
   array([[10, 11, 12],
          [13, 14, 15]])
   >>> g(a, np.array([10, 20]))
   array([[10, 11, 12],
          [23, 24, 25]])


.. note::
   Both :func:`~numba.vectorize` and :func:`~numba.guvectorize` support
   passing ``nopython=True`` :ref:`as in the @jit decorator <jit-nopython>`.
   Use it to ensure the generated code does not fallback to
   :term:`object mode`.

.. _overwriting-input-values:

Overwriting input values
------------------------

In most cases, writing to inputs may also appear to work - however, this
behaviour cannot be relied on. Consider the following example function::

   @guvectorize([(float64[:], float64[:])], '()->()')
   def init_values(invals, outvals):
       invals[0] = 6.5
       outvals[0] = 4.2

Calling the `init_values` function with an array of `float64` type results in
visible changes to the input::

   >>> invals = np.zeros(shape=(3, 3), dtype=np.float64)
   >>> outvals = init_values(invals)
   >>> invals
   array([[6.5, 6.5, 6.5],
          [6.5, 6.5, 6.5],
          [6.5, 6.5, 6.5]])
   >>> outvals
   array([[4.2, 4.2, 4.2],
       [4.2, 4.2, 4.2],
       [4.2, 4.2, 4.2]])

This works because NumPy can pass the input data directly into the `init_values`
function as the data `dtype` matches that of the declared argument.  However, it
may also create and pass in a temporary array, in which case changes to the
input are lost. For example, this can occur when casting is required. To
demonstrate, we can  use an array of `float32` with the `init_values` function::

   >>> invals = np.zeros(shape=(3, 3), dtype=np.float32)
   >>> outvals = init_values(invals)
   >>> invals
   array([[0., 0., 0.],
          [0., 0., 0.],
          [0., 0., 0.]], dtype=float32)

In this case, there is no change to the `invals` array because the temporary
casted array was mutated instead.

.. _dynamic-universal-functions:

Dynamic universal functions
===========================

As described above, if you do not pass any signatures to the
:func:`~numba.vectorize` decorator, your Python function will be used
to build a dynamic universal function, or :class:`~numba.DUFunc`.  For
example::

   from numba import vectorize

   @vectorize
   def f(x, y):
       return x * y

The resulting :func:`f` is a :class:`~numba.DUFunc` instance that
starts with no supported input types.  As you make calls to :func:`f`,
Numba generates new kernels whenever you pass a previously unsupported
input type.  Given the example above, the following set of interpreter
interactions illustrate how dynamic compilation works::

   >>> f
   <numba._DUFunc 'f'>
   >>> f.ufunc
   <ufunc 'f'>
   >>> f.ufunc.types
   []

The example above shows that :class:`~numba.DUFunc` instances are not
ufuncs.  Rather than subclass ufunc's, :class:`~numba.DUFunc`
instances work by keeping a :attr:`~numba.DUFunc.ufunc` member, and
then delegating ufunc property reads and method calls to this member
(also known as type aggregation).  When we look at the initial types
supported by the ufunc, we can verify there are none.

Let's try to make a call to :func:`f`::

   >>> f(3,4)
   12
   >>> f.types   # shorthand for f.ufunc.types
   ['ll->l']

If this was a normal Numpy ufunc, we would have seen an exception
complaining that the ufunc couldn't handle the input types.  When we
call :func:`f` with integer arguments, not only do we receive an
answer, but we can verify that Numba created a loop supporting C
:code:`long` integers.

We can add additional loops by calling :func:`f` with different inputs::

   >>> f(1.,2.)
   2.0
   >>> f.types
   ['ll->l', 'dd->d']

We can now verify that Numba added a second loop for dealing with
floating-point inputs, :code:`"dd->d"`.

If we mix input types to :func:`f`, we can verify that `Numpy ufunc
casting rules`_ are still in effect::

   >>> f(1,2.)
   2.0
   >>> f.types
   ['ll->l', 'dd->d']

.. _`Numpy ufunc casting rules`: http://docs.scipy.org/doc/numpy/reference/ufuncs.html#casting-rules

This example demonstrates that calling :func:`f` with mixed types
caused Numpy to select the floating-point loop, and cast the integer
argument to a floating-point value.  Thus, Numba did not create a
special :code:`"dl->d"` kernel.

This :class:`~numba.DUFunc` behavior leads us to a point similar to
the warning given above in "`The @vectorize decorator`_" subsection,
but instead of signature declaration order in the decorator, call
order matters.  If we had passed in floating-point arguments first,
any calls with integer arguments would be cast to double-precision
floating-point values.  For example::

   >>> @vectorize
   ... def g(a, b): return a / b
   ...
   >>> g(2.,3.)
   0.66666666666666663
   >>> g(2,3)
   0.66666666666666663
   >>> g.types
   ['dd->d']

If you require precise support for various type signatures, you should
specify them in the :func:`~numba.vectorize` decorator, and not rely
on dynamic compilation.

Dynamic generalized universal functions
=======================================

Similar to a dynamic universal function, if you do not specify any types to
the :func:`~numba.guvectorize` decorator, your Python function will be used
to build a dynamic generalized universal function, or :class:`~numba.GUFunc`.
For example::

   from numba import guvectorize

   @guvectorize('(n),()->(n)')
   def g(x, y, res):
       for i in range(x.shape[0]):
           res[i] = x[i] + y

We can verify the resulting function :func:`g` is a :class:`~numba.GUFunc`
instance that starts with no supported input types. For instance::

   >>> g
   <numba._GUFunc 'g'>
   >>> g.ufunc
   <ufunc 'g'>
   >>> g.ufunc.types
   []

Similar to a :class:`~numba.DUFunc`, as one make calls to :func:`g()`,
numba generates new kernels for previously unsupported input types. The
following set of interpreter interactions will illustrate how dynamic
compilation works for a :class:`~numba.GUFunc`::

   >>> x = np.arange(5, dtype=np.int64)
   >>> y = 10
   >>> res = np.zeros_like(x)
   >>> g(x, y, res)
   >>> res
   array([5, 6, 7, 8, 9])
   >>> g.types
   ['ll->l']

If this was a normal :func:`guvectorize` function, we would have seen an
exception complaining that the ufunc could not handle the given input types.
When we call :func:`g()` with the input arguments, numba creates a new loop
for the input types.

We can add additional loops by calling :func:`g` with new arguments::

   >>> x = np.arange(5, dtype=np.double)
   >>> y = 2.2
   >>> res = np.zeros_like(x)
   >>> g(x, y, res)

We can now verify that Numba added a second loop for dealing with
floating-point inputs, :code:`"dd->d"`.

   >>> g.types  # shorthand for g.ufunc.types
   ['ll->l', 'dd->d']

One can also verify that Numpy ufunc casting rules are working as expected::

   >>> x = np.arange(5, dtype=np.int64)
   >>> y = 2.2
   >>> res = np.zeros_like(x)
   >>> g(x, y, res)
   >>> res

If you need precise support for various type signatures, you should not rely on dynamic
compilation and instead, specify the types them as first
argument in the :func:`~numba.guvectorize` decorator.
.. _numba-5_mins:

A ~5 minute guide to Numba
==========================

Numba is a just-in-time compiler for Python that works best on code that uses
NumPy arrays and functions, and loops. The most common way to use Numba is
through its collection of decorators that can be applied to your functions to
instruct Numba to compile them. When a call is made to a Numba-decorated
function it is compiled to machine code "just-in-time" for execution and all or
part of your code can subsequently run at native machine code speed!

Out of the box Numba works with the following:

* OS: Windows (32 and 64 bit), OSX, Linux (32 and 64 bit). Unofficial support on
  \*BSD.
* Architecture: x86, x86_64, ppc64le, armv7l, armv8l (aarch64). Unofficial
  support on M1/Arm64.
* GPUs: Nvidia CUDA.
* CPython
* NumPy 1.18 - latest

How do I get it?
----------------
Numba is available as a `conda <https://conda.io/docs/>`_ package for the
`Anaconda Python distribution <https://www.anaconda.com/>`_::

  $ conda install numba

Numba also has wheels available::

  $ pip install numba

Numba can also be
:ref:`compiled from source <numba-source-install-instructions>`, although we do
not recommend it for first-time Numba users.

Numba is often used as a core package so its dependencies are kept to an
absolute minimum, however, extra packages can be installed as follows to provide
additional functionality:

* ``scipy`` - enables support for compiling ``numpy.linalg`` functions.
* ``colorama`` - enables support for color highlighting in backtraces/error
  messages.
* ``pyyaml`` - enables configuration of Numba via a YAML config file.
* ``icc_rt`` - allows the use of the Intel SVML (high performance short vector
  math library, x86_64 only). Installation instructions are in the
  :ref:`performance tips <intel-svml>`.

Will Numba work for my code?
----------------------------
This depends on what your code looks like, if your code is numerically
orientated (does a lot of math), uses NumPy a lot and/or has a lot of loops,
then Numba is often a good choice. In these examples we'll apply the most
fundamental of Numba's JIT decorators, ``@jit``, to try and speed up some
functions to demonstrate what works well and what does not.

Numba works well on code that looks like this::

    from numba import jit
    import numpy as np

    x = np.arange(100).reshape(10, 10)

    @jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
    def go_fast(a): # Function is compiled to machine code when called the first time
        trace = 0.0
        for i in range(a.shape[0]):   # Numba likes loops
            trace += np.tanh(a[i, i]) # Numba likes NumPy functions
        return a + trace              # Numba likes NumPy broadcasting

    print(go_fast(x))


It won't work very well, if at all, on code that looks like this::

    from numba import jit
    import pandas as pd

    x = {'a': [1, 2, 3], 'b': [20, 30, 40]}

    @jit
    def use_pandas(a): # Function will not benefit from Numba jit
        df = pd.DataFrame.from_dict(a) # Numba doesn't know about pd.DataFrame
        df += 1                        # Numba doesn't understand what this is
        return df.cov()                # or this!

    print(use_pandas(x))

Note that Pandas is not understood by Numba and as a result Numba would simply
run this code via the interpreter but with the added cost of the Numba internal
overheads!

What is ``nopython`` mode?
--------------------------
The Numba ``@jit`` decorator fundamentally operates in two compilation modes,
``nopython`` mode and ``object`` mode. In the ``go_fast`` example above,
``nopython=True`` is set in the ``@jit`` decorator; this is instructing Numba to
operate in ``nopython`` mode. The behaviour of the ``nopython`` compilation mode
is to essentially compile the decorated function so that it will run entirely
without the involvement of the Python interpreter. This is the recommended and
best-practice way to use the Numba ``jit`` decorator as it leads to the best
performance.

Should the compilation in ``nopython`` mode fail, Numba can compile using
``object mode``. This is a fall back mode for the ``@jit`` decorator if
``nopython=True`` is not set (as seen in the ``use_pandas`` example above). In
this mode Numba will identify loops that it can compile and compile those into
functions that run in machine code, and it will run the rest of the code in the
interpreter. For best performance avoid using this mode!

How to measure the performance of Numba?
----------------------------------------
First, recall that Numba has to compile your function for the argument types
given before it executes the machine code version of your function. This takes
time. However, once the compilation has taken place Numba caches the machine
code version of your function for the particular types of arguments presented.
If it is called again with the same types, it can reuse the cached version
instead of having to compile again.

A really common mistake when measuring performance is to not account for the
above behaviour and to time code once with a simple timer that includes the
time taken to compile your function in the execution time.

For example::

    from numba import jit
    import numpy as np
    import time

    x = np.arange(100).reshape(10, 10)

    @jit(nopython=True)
    def go_fast(a): # Function is compiled and runs in machine code
        trace = 0.0
        for i in range(a.shape[0]):
            trace += np.tanh(a[i, i])
        return a + trace

    # DO NOT REPORT THIS... COMPILATION TIME IS INCLUDED IN THE EXECUTION TIME!
    start = time.time()
    go_fast(x)
    end = time.time()
    print("Elapsed (with compilation) = %s" % (end - start))

    # NOW THE FUNCTION IS COMPILED, RE-TIME IT EXECUTING FROM CACHE
    start = time.time()
    go_fast(x)
    end = time.time()
    print("Elapsed (after compilation) = %s" % (end - start))

This, for example prints::

    Elapsed (with compilation) = 0.33030009269714355
    Elapsed (after compilation) = 6.67572021484375e-06

A good way to measure the impact Numba JIT has on your code is to time execution
using the `timeit <https://docs.python.org/3/library/timeit.html>`_ module
functions; these measure multiple iterations of execution and, as a result,
can be made to accommodate for the compilation time in the first execution.

As a side note, if compilation time is an issue, Numba JIT supports
:ref:`on-disk caching <jit-decorator-cache>` of compiled functions and also has
an :ref:`Ahead-Of-Time <aot-compilation>` compilation mode.

How fast is it?
---------------
Assuming Numba can operate in ``nopython`` mode, or at least compile some loops,
it will target compilation to your specific CPU. Speed up varies depending on
application but can be one to two orders of magnitude. Numba has a
:ref:`performance guide <performance-tips>` that covers common options for
gaining extra performance.

How does Numba work?
--------------------
Numba reads the Python bytecode for a decorated function and combines this with
information about the types of the input arguments to the function. It analyzes
and optimizes your code, and finally uses the LLVM compiler library to generate
a machine code version of your function, tailored to your CPU capabilities. This
compiled version is then used every time your function is called.

Other things of interest:
-------------------------
Numba has quite a few decorators, we've seen ``@jit``, but there's
also:

* ``@njit`` - this is an alias for ``@jit(nopython=True)`` as it is so commonly
  used!
* ``@vectorize`` - produces NumPy ``ufunc`` s (with all the ``ufunc`` methods
  supported). :ref:`Docs are here <vectorize>`.
* ``@guvectorize`` - produces NumPy generalized ``ufunc`` s.
  :ref:`Docs are here <guvectorize>`.
* ``@stencil`` - declare a function as a kernel for a stencil like operation.
  :ref:`Docs are here <numba-stencil>`.
* ``@jitclass`` - for jit aware classes. :ref:`Docs are here <jitclass>`.
* ``@cfunc`` - declare a function for use as a native call back (to be called
  from C/C++ etc). :ref:`Docs are here <cfunc>`.
* ``@overload`` - register your own implementation of a function for use in
  nopython mode, e.g. ``@overload(scipy.special.j0)``.
  :ref:`Docs are here <high-level-extending>`.

Extra options available in some decorators:

* ``parallel = True`` - :ref:`enable <jit-decorator-parallel>` the
  :ref:`automatic parallelization <numba-parallel>` of the function.
* ``fastmath = True`` - enable :ref:`fast-math <jit-decorator-fastmath>`
  behaviour for the function.

ctypes/cffi/cython interoperability:

* ``cffi`` - The calling of :ref:`CFFI  <cffi-support>` functions is supported
  in ``nopython`` mode.
* ``ctypes`` - The calling of :ref:`ctypes  <ctypes-support>` wrapped
  functions is supported in ``nopython`` mode.
* Cython exported functions :ref:`are callable <cython-support>`.

GPU targets:
~~~~~~~~~~~~

Numba can target `Nvidia CUDA <https://developer.nvidia.com/cuda-zone>`_ GPUs.
You can write a kernel in pure Python and have Numba handle the computation and
data movement (or do this explicitly). Click for Numba documentation on
:ref:`CUDA <cuda-index>`.
