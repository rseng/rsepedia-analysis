

## Removed
 
 * `noodles.file` unused module
 
 
## Changed

 * Added `SerPath` class to serial namespace


## Fixed
 * `Fail` class serialization
---
title: Noodles - parallel programming in Python
---

[![Travis](https://travis-ci.org/NLeSC/noodles.svg?branch=master)](https://travis-ci.org/NLeSC/noodles)
[![Zenodo DOI](https://zenodo.org/badge/45391130.svg)](https://zenodo.org/badge/latestdoi/45391130)
[![Code coverage](https://codecov.io/gh/NLeSC/noodles/branch/master/graph/badge.svg)](https://codecov.io/gh/NLeSC/noodles)
[![Documentation](https://readthedocs.org/projects/noodles/badge/?version=latest)](https://noodles.readthedocs.io/en/latest/?badge=latest)

::: {.splash}
* Write readable code
* Parallelise with a dash of Noodle sauce!
* Scale your applications from laptop to HPC using Xenon
    + [Learn more about Xenon](https://xenon-middleware.github.io/xenon)
* Read our [documentation](https://noodles.rtfd.io/), including tutorials on:
    + [Creating parallel programs](https://noodles.readthedocs.io/en/latest/poetry_tutorial.html)
    + [Circumventing the global interpreter lock](https://noodles.readthedocs.io/en/latest/prime_numbers.html)
    + [Handling errors in a meaningful way](https://noodles.readthedocs.io/en/latest/errors.html)
    + [Serialising your data](https://noodles.readthedocs.io/en/latest/serialisation.html)
    + [Functional programming and flow control](https://noodles.readthedocs.io/en/latest/control_your_flow.html)
:::

# What is Noodles?

Noodles is a task-based parallel programming model in Python that offers the same intuitive interface when running complex workflows on your laptop or on large computer clusters.

# Installation
To install the latest version from PyPI:

```
pip install noodles
```

To enable the Xenon backend for remote job execution,

```
pip install noodles[xenon]
```

This requires a Java Runtime to be installed, you may check this by running

```
java --version
```

which should print the version of the currently installed JRE.


# Documentation
All the latest documentation is available on [Read the Docs](https://noodles.rtfd.io/).

# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

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

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at j.hidding@esciencecenter.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
Thank you for being interested to contribute to Noodles!

## Issues / Pull requests
If you find a bug or unexpected behaviour in Noodles you are welcome to open an issue.
We'll do our best to address the issue if it is within our capacity to do so.

Pull requests are certainly welcome. Please first open an issue outlining the bug or feature request that is being addressed.

All contributions will be integrated per the Apache License agreement (see LICENSE); contributing authors will so be attributed.
# Howto run

From this folder say:
    > PYTHONPATH=".." python3 example<n>.py

Noodles build system: BOIL
==========================

This example shows how to enhance your daily workflow with Noodles! The program
reads a file `boil.ini` from the current directory and starts to compile your
favourite C/C++ project. Dependencies for each source file are being checked
using `gcc -MM` and files are only compiled when the dependcies are newer.
The user can give the `-j` parameter to determine how many jobs should run
simultaneously.

If you look at the BOIL source code, you may notice that very little of the
parallelism of this code creeps into the basic logic of the program.


Testing BOIL
~~~~~~~~~~~~

To test BOIL we provided a small C++ program. There is a main program and
a unittest which can be compiled separately using the given configuration
in `boil.ini`.==============================================
SOBA: Utility for non-directed graph execution
==============================================

Noodles gives a way to execute directed acyclic graphs. There are use cases where
the graph just gives us mutual exclusion of jobs, because they are writing to the
same output location (either memory or disk). In this case we want to do the
scheduling of jobs dynamically, with the exclusion information added as meta-data.

For this example, we have the graph being fed as a JSON file. We add a worker broker
to the pool to deal with just this problem.
Development documentation
=========================
.. automodule:: noodles
    :members:

Internal Specs
--------------
.. automodule:: noodles.workflow
    :members:

Promised object
---------------
.. automodule:: noodles.interface
    :members:

Runners
-------
.. automodule:: noodles.run.scheduler
    :members:

.. automodule:: noodles.run.hybrid
    :members:

Serialisation
-------------
.. automodule:: noodles.serial
    :members:

.. automodule:: noodles.serial.registry
    :members:

Worker executable
-----------------
.. automodule:: noodles.worker
    :members:

Streams
-------
.. automodule:: noodles.lib
    :members:
Implementation
==============

.. toctree::
    :maxdepth: 2

    development
    scheduler
.. highlight:: python
    :linenothreshold: 5


Real World Tutorial 2: Boil, a make tool.
=========================================

.. toctree::

    boil_source

This tutorial should teach you the basics of using Noodles to parallelise a Python program. We will see some of the quirks that come with programming in a strict functional environment.

"make -j"
~~~~~~~~~

The first of all workflow engines is called 'make'. It builds a tree of interdepending jobs and executes them in order. The '-j' option allows the user to specify a number of jobs to be run simultaneously to speed up execution. The syntax of make is notoriously terse, partly due to it's long heritage (from 1976). This example shows how we can write a script that compiles a C/C++ program using GCC or CLANG, then how we can parallelise it using Noodles. We have even included a small C++ program to play with!

Functions
~~~~~~~~~

Compiling a C program is a two-stage process. First we need to compile each source file to an object file, then link these object files to an executable. The :py:func:`compile_source` function looks like this::

    @schedule
    def compile_source(source_file, object_file, config):
        p = subprocess.run(
            [config['cc'], '-c'] + config['cflags'].split() \
                + [source_file, '-o', object_file],
            stderr=subprocess.PIPE, universal_newlines=True)
        p.check_returncode()

        return object_file

It takes a source file, an object file and a configuration object. This configuration contains all the information on which compiler to use, with which flags, and so on. If the compilation was successful, the name of the object file is returned. This last bit is crucial if we want to include this function in a workflow.

.. NOTE:: Each dependency (in this case linking to an excecutable requires compilation first) should be represented by return values of one function ending up as arguments to another function.

The function for linking object files to an executable looks very similar::

    @schedule
    def link(object_files, config):
        p = subprocess.run(
            [config['cc']] + object_files + ['-o', config['target']] \
                + config['ldflags'].split(),
            stderr=subprocess.PIPE, universal_newlines=True)
        p.check_returncode()

        return config['target']

In this case the function takes a list of object file names and the same configuration object that we saw before. Again, this function returns the name of the target executable file. The caller of this function already knows the name of the target file, but we need it to track dependencies between function calls.

Since both the :py:func:`link` and the :py:func:`compile_source` functions do actual work that we'd like to see being done in a concurrent environment, they need to be decorated with the :py:func:`schedule` decorator.

One of the great perks of using Make, is that it skips building any files that are already up-to-date with the source code. If our little build script is to compete with such efficiency, we should do the same!

::

    def is_out_of_date(f, deps):
        if not os.path.exists(f):
            return True

        f_stat = os.stat(f)

        for d in deps:
            d_stat = os.stat(d)

            if d_stat.st_mtime_ns > f_stat.st_mtime_ns:
                return True

        return False

This function takes a file `f` and a list of files `deps` and checks the modification dates of all of the files in `deps` against that of `f`.

One of the *quirks*, that we will need to deal with, is that some *logic* in a program needs to have knowledge of the actual objects that are computed, not just the possibility of such an object in the future. When designing programs for Noodles, you need to be aware that such logic can only be performed *inside* the functions. Say we have a condition under which one or the other action needs to be taken, and this condition depends on the outcome of a previous element in the workflow. The actual Python `if` statement evaluating this condition needs to be inside a scheduled function. One way around this, is to write a wrapper::

    @schedule
    def cond(predicate, when_true, when_false):
        if predicate:
            return when_true
        else:
            return when_false

However, there is a big problem with this approach! The Noodles engine sees two arguments to the `cond` function that it wants to evaluate before heading into the call to `cond`. Both arguments will be evaluated before we can even decide which of the two we should use! We will present a solution to this problem at a later stage, but for now we will have to work our way around this by embedding the conditional logic in a more specific function.


In this case we have the function `is_out_of_date` that determines whether we need to recompile a file or leave it as it is. The second stage, linking the object files to an executable, only needs to happen if any of the object files is younger than the executable. However these object files are part of the logic inside the workflow! The conditional execution of the linker needs to be called by another scheduled function::

    @schedule
    def get_target(obj_files, config):
        if is_out_of_date(config['target'], obj_files):
            return link(obj_files, config)
        else:
            return message("target is up-to-date.")

Since we need the answer to :py:func:`is_out_of_date` in the present, the :py:func:`is_out_of_date` function cannot be a scheduled function. Python doesn't know the truth value of a :py:class:`PromisedObject`. The `message` function is not a special function; it just prints a message and returns a value (optional second argument). We still need to optionalise the compilation step. Since all of the information needed to decide whether to compile or not is already present, we can make this a normal Python function; there is no need to schedule anything (even though everything would still work if we did)::

    def get_object_file(src_dir, src_file, config):
        obj_path = object_filename(src_dir, src_file, config)
        src_path = os.path.join(src_dir, src_file)

        deps = dependencies(src_path, config)
        if is_out_of_date(obj_path, deps):
            return compile_source(src_path, obj_path, config)
        else:
            return obj_path


The `object_filename` is a little helper function creating a sensible name for the object file; also it makes sure that the directory in which the object file is placed exists. `dependencies` Runs the compiler with '-MM' flags to obtain the header dependencies of the C-file.

We are now ready to put these functions together in a workflow!

::

    def make_target(config):
        dirs = [config['srcdir']] + [
            os.path.normpath(os.path.join(config['srcdir'], d))
            for d in config['modules'].split()
        ]

        files = chain(*(
            find_files(d, config['ext'])
            for d in dirs)
        )

        obj_files = noodles.gather(*[
            get_object_file(src_dir, src_file, config)
            for src_dir, src_file in files
            ])

        return get_target(obj_files, config)

Let's go through this step-by-step. The `make_target` function takes one argument, the config object. We obtain from the configuration the directories to search for source files. We then search these directories for any files with the correct file extension, stored in `config['ext']`. The variable `files` now contains a list of pairs, each pair having a directory and file name. So far we have not yet used any Noodles code.

Next we pass each source file through the `get_object_file` function in a list comprehension. The resulting list contains both :py:class:`PromisedObject` and string values; strings for all the object files that are already up-to-date. To pass this list to the linking stage we have to make sure that Noodles understands that the list is something that is being promised. If we were to pass it as is, Noodles just sees a list as an argument to `get_target` and doesn't look any deeper.

.. NOTE:: Every `PromisedObject` has to be passed as an argument to a scheduled function in order to be evaluated. To pass a list to a scheduled function, we have to convert the list of promises into a promise of a list.

The function :py:func:`gather` solves this little problem; it's definition is very simple::

    @schedule
    def gather(*args):
        return list(args)

Now that the variable `obj_files` is a :py:class:`PromisedObject`, we can pass it to `get_target`, giving us the final workflow. Running this workflow can be as simple as ``run_single(wf)`` or ``run_parallel(wf, n_threads=4)``.

Friendly output and error handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The code as defined above will run, but if the compiler gives error messages it will crash in a very ugly manner. Noodles has some features that will make our fledgeling Make utility much prettier. We can decorate our functions further with information on how to notify the user of things happening::

    @schedule_hint(display="Compiling {source_file} ... ",
                   confirm=True)
    def compile_source(source_file, object_file, config):
        pass

The :py:func:`schedule_hint` decorator attaches hints to the scheduled function. These hints can be used in any fashion we like, depending on the workers that we use to evaluate the workflow. In this case we use the :py:func:`run_logging` worker, with the :py:class:`SimpleDisplay` class to take care of screen output::

    with SimpleDisplay(error_filter) as display:
        noodles.run_logging(work, args.n_threads, display)

The `error_filter` function determines what errors are expected and how we print the exceptions that are caught. In our case we expect exceptions of type :py:exc:`subprocess.CalledProcessError` in the case of a compiler error. Otherwise the exception is unexpected and should be treated as a bug in Boil!

::

    def error_filter(xcptn):
        if isinstance(xcptn, subprocess.CalledProcessError):
            return xcptn.stderr
        else:
            return None

The `display` tag in the hinst tells the `display` object to print a text when the job is started.
The `confirm` flag in the hints tells the `display` object that a function is error-prone and to draw a green checkmark on success and a red X in case of failure.

Conclusion
~~~~~~~~~~

You should now be able to fully understand the sourcecode of Boil! Try it out on the sample code provided:

.. code-block:: bash

    > ./boil --help
    
.. code-block:: none

    usage: boil [-h] [-j N_THREADS] [-dumb] target

    Compile software. Configuration is in 'boil.ini'.

    positional arguments:
      target        target to build, give 'list' to list targets.

    optional arguments:
      -h, --help    show this help message and exit
      -j N_THREADS  number of threads to run simultaneously.
      -dumb         print info without special term codes.

.. code-block:: bash

    > cat boil.ini

.. code-block:: ini

    [generic]
    objdir = obj
    ldflags = -lm
    cflags = -g -std=c++11 -O2 -fdiagnostics-color -Wpedantic
    cc = g++
    ext = .cc

    [main]
    srcdir = main
    target = hello
    modules = ../src

    [test]
    srcdir = test
    target = unittests
    modules = ../src

    [clean]
    command = rm -rf ${generic:objdir} ${test:target} ${main:target}


.. code-block:: bash
    
    > ./boil main -j4

.. code-block:: none

    ╭─(Building target hello)
    │   Compiling src/mandel.cc ...                   (✔)
    │   Compiling src/common.cc ...                   (✔)
    │   Compiling src/render.cc ...                   (✔)
    │   Compiling src/iterate.cc ...                  (✔)
    │   Compiling src/julia.cc ...                    (✔)
    │   Compiling main/main.cc ...                    (✔)
    │   Linking ...                                   (✔)
    ╰─(success)

Eating noodles (user docs)
==========================

The primary goal of the *noodles* library is to ease the construction and execution of *computational workflows* using the Python language. This library is meant for scientists who want to perform complex compute-intensive tasks on parallel/distributed infrastructures in a readable, scalable and sustainable/reproducible? manner.
A workflow is commonly modelled as a *directed acyclic graph* (DAG or simply graph) in which the computations are represented as nodes whereas the dependencies between them are represented as directed edges (indicating data transport).

A first example
---------------

Let's look at a small example of creating a diamond workflow, which consists of simple (arithmetic) functions:

.. code:: python

    from noodles import run_single
    from noodles.tutorial import (add, sub, mul)

    u = add(5, 4)
    v = sub(u, 3)
    w = sub(u, 2)
    x = mul(v, w)

    answer = run_single(x)

    print("The answer is {0}.".format(answer))

That allmost looks like normal Python! The only difference is the :py:func:`~noodles.run_single` statement at the end of this program.
The catch is that none of the computation is actually done until the :py:func:`~noodles.run_single` statement has been given.
The variables ``u``, ``v``, ``w``, and ``x`` only represent the *promise* of a value.
The functions that we imported are wrapped, such that they construct the directed acyclic graph of the computation in stead of just computing the result immediately.
This DAG then looks like this:

.. figure:: _static/images/dag1.png
    :alt: the diamond workflow DAG
    :align: center
    :figwidth: 50%

    The diamond workflow.

Running this program will first evaluate the result to ``add(5, 4)``.
The resulting value is then inserted into the empty slots in the depending nodes.
Each time a node has no empty slots left, it is scheduled for evaluation.
At the end, the program should print:

::

    The answer is 42.

At this point it is good to know what the module ``noodles.tutorial`` looks like.
It looks very simple.
However, you should be aware of what happens behind the curtains, to understand the limitations of this approach.

.. code:: python

    from noodles import schedule

    @schedule
    def add(a, b):
      """Adding up numbers! (is very uplifting)"""
      return a + b

    @schedule
    def sub(a, b):
      """Subtracting numbers."""
      return a - b

    @schedule
    def mul(a, b):
      """Multiplying numbers."""
      return a * b

    ...

The :py:func:`@schedule <noodles.schedule>` decorators takes care that the functions actually return *promises* instead of values.
Such a :py:class:`~noodles.interface.PromisedObject` is a placeholder for the expected result.
It stores the workflow graph that is needed to compute the promise.
When another :py:func:`schedule <noodles.schedule>`-decorated function is called with a promise, the graphs of the dependencies are merged to create a new workflow graph.

.. NOTE:: The promised object can be of any type and can be used as a normal object.
          You access attributes and functions of the object that is promised as you normally would.
          Be aware, however, that it is important to program in a functional way, so changing the attributes of a promised object is not a good idea.
          Instead, return a copy of the object with the changed values.


Doing things parallel
~~~~~~~~~~~~~~~~~~~~~

Using the Noodles approach it becomes very easy to paralellise computations. Let's look at a second example.

.. code:: python

    from noodles import (gather, run_parallel)
    from noodles.tutorial import (add, sub, mul, accumulate)


    def my_func(a, b, c):
        d = add(a, b)
        return mul(d, c)


    u = add(1, 1)
    v = sub(3, u)
    w = [my_func(i, v, u) for i in range(6)]
    x = accumulate(gather(*w))

    answer = run_parallel(x, n_threads=4)

    print("The answer is {0}, again.".format(answer))

This time the workflow graph will look a bit more complicated.

.. figure:: _static/images/dag2.png
    :alt: the workflow graph of the second example
    :align: center
    :figwidth: 100%

    The workflow graph of the second example.

Here we see how a user can define normal python functions and use them to build a larger workflow.
Furthermore, we introduce a new bit of magic: the :py:func:`gather <noodles.gather>` function.
When you build a list of computations using a list-comprehension like above, you essentially store a *list of promises* in variable ``w``.
However, schedule-decorated functions cannot easily see which arguments contain promised values, such as ``w``, and which arguments are plain Python.
The :py:func:`gather <noodles.gather>` function converts the list of promises into a promise of a list, making it clear to the scheduled function this argument is a promise.
The :py:func:`gather <noodles.gather>` function is defined as follows:

.. code:: python

    @schedule
    def gather(*lst):
        return lst

By unpacking the list (by doing ``gather(*w)``) in the call to gather, each item in ``w`` becomes a dependency of the ``gather`` node in this workflow, as we can see in the figure above.

To make use of the parallelism in this workflow, we run it with :py:func:`~noodles.run_parallel`.
This runner function creates a specified number of threads, each taking jobs from the Noodles scheduler and returning results.

Running workflows
-----------------

Noodles ships with a few ready-made functions that run the workflow for you, depending on the amount of work that needs to be done.

:py:func:`~noodles.run_single`, local single thread
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Runs your workflow in the same thread as the caller.
This function is mainly for testing.
When running workflows you almost always want to use one of the other functions.

:py:func:`~noodles.run_parallel`, local multi-thread
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Runs your workflow in parallel using any number of threads.
Usually, specifying the number of cores in your CPU will give optimal performance for this runner.

.. NOTE:: If you are very **very** certain that your workflow will never need to scale to cluster-computing, this runner is more lenient on the kinds of Python that is supported, because function arguments are not converted to and from JSON. Think of nested functions, lambda forms, generators, etc.

:py:func:`~noodles.run_process`, local multi-process
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Starts a second process to run jobs. This is usefull for testing the JSON compatability of your workflow on your own machine.

Xenon
~~~~~
Xenon_ is a Java library offering a uniform interface to all manners of job schedulers. Running a job on your local machine is as easy as submitting it to SLURM or Torque on your groceries supercomputer. To talk to Xenon from Python we use pyxenon_.

Using the Xenon runner, there are two modes of operation: *batch* and *online*. In online mode, jobs are streamed to the worker and results read back. If your laptop crashes while an online computation is running, that is to say, the connection is broken, the worker dies and you may lose results. Getting the online mode to be more robust is one of the aims for upcomming releases.

The Xenon runner needs a way to setup the virtualenv on the remote side, so a worker script needs to be specified. We have included a bash-script ``worker.sh`` that should work in the simplest cases.

.. code-block:: bash

    #!/bin/bash

    # run in the directory where the script is located
    cd "$(dirname "${BASH_SOURCE[0]}")"

    # activate the virtualenv that is given as first argument
    # invoking this script.
    if [ -e $1/bin/activate ]; then
    	source $1/bin/activate;
    fi

    # start the worker with the rest of the arguments.
    # stderr is written to a file.
    python -m noodles.worker ${@:2} 2> errlog

    # close the virtualenv.
    if [ -z ${VIRTUAL_ENV+x} ]; then
    	deactivate;
    fi

If you need to setup some more aspects of the environment, load modules, set variables etc., modify this script and put it in the directory where you want to run the jobs. Specify this directory in the Python script.

.. code:: python

    from noodles import (
        serial, gather)
    from noodles.run.xenon import (
        XenonConfig, RemoteJobConfig, XenonKeeper, run_xenon_prov)
    from noodles.display import (
        NCDisplay)

    from noodles.tutorial import add, accumulate

    if __name__ == "__main__":
        a = [add(i, j) for i in range(5) for j in range(5)]
        b = accumulate(gather(*a))

        # XenonKeeper is the root Xenon object that gives access
        # to the Xenon Java library
        with XenonKeeper() as Xe:
            # We recommend loging in on your compute resource
            # through private/public key pairs. This prevents
            # passwords ending up as ASCII in your source files.
            certificate = Xe.credentials.newCertificateCredential(
                'ssh', os.environ['HOME'] + '/.ssh/id_rsa', '<username>', '', None)

            # Configure Xenon to access your favourite super computer.
            xenon_config = XenonConfig(
                jobs_scheme='slurm',
                location='login.super-duper-computer.darpa.net',
                credential=certificate
            )

            # Specify how to submit jobs.
            job_config = RemoteJobConfig(
                registry=serial.base,
                prefix='<path-to-virtualenv>',
                working_dir='<project-path>',
                time_out=5000
            )


            # Run jobs with NCurses based console feedback
            with NCDisplay() as display:
                result = run_xenon_prov(
                    b, Xe, "cache.json", 2, xenon_config, job_config,
                    display=display)

        print("This test is working {0}%!".format(result))


Hybrid mode
~~~~~~~~~~~
We may have a situation where a workflow consists of some very heavy *compute* jobs and a lot of smaller jobs that do some bookkeeping. If we were to schedule all the menial jobs to a SLURM queue we actually slow down the computation through the overhead of job submission. The Noodles cook may provide the schedule functions with hints on the type of job the function represents. Depending on these hints we may dispatch the job to a remote worker or keep it on the local machine.

We provide an example on how to use the hybrid worker in the source.

If you really need to, it is not too complicated to develop your own job runner based on some of these examples. Elsewhere in this documentation we elaborate on the architecture and interaction between runners and the scheduler, see: :ref:`noodles-scheduler`.

.. _Xenon: http://nlesc.github.io/Xenon/
.. _pyxenon: http://github.com/NLeSC/pyxenon
Tutorials
=========

.. toctree::
    :maxdepth: 2

    first_steps
    poetry_tutorial
    boil_tutorial
    prime_numbers
    errors
    serialisation
    control_your_flow
.. highlight:: python
    :linenothreshold: 5

.. _noodles-scheduler:

The Noodles Scheduler
=====================

The Noodles scheduler is completely separated from the worker infrastructure. The scheduler accepts a single worker as an argument. This worker provides the scheduler with two coroutines. One acts as a generator of results, the other as a sink for jobs (the scheduler calls the ``send()`` method on it).

Both jobs and results are accompanied by a unique key to identify the associated job. The scheduler loops over the results as follows (more or less):

.. code-block:: python

    for (key, result) in source:
        """process result"""
        ...

        for node in workflow.nodes:
            if node.ready():
                sink.send((node.key, node.job))

Local workers
-------------

The single worker
~~~~~~~~~~~~~~~~~

It is the responsibility of the worker to keep a queue where so desired. A single result may trigger many new nodes to be ready for evaluation. This means that either the jobs or the results must be buffered in a queue. In the simplest case we have a single worker in the same thread as the scheduler.

.. figure:: _static/images/sd-single.svg
    :alt: sequence diagram, single thread
    :align: center

    Sequence diagram for a single threaded execution model.

The worker code looks like this:

::

    from noodles.coroutines import (IOQueue, Connection)
    from noodles.run_common import run_job

    def single_worker():
        """Sets up a single worker co-routine."""
        jobs = IOQueue()

        def get_result():
            source = jobs.source()

            for key, job in source:
                yield (key, run_job(job))

        return Connection(get_result, jobs.sink)

The ``IOQueue`` class wraps a standard Python queue. It provides a ``sink`` member pushing elements onto the queue, and a ``source`` member yielding elements from the queue, calling ``Queue.task_done()`` when the coroutine regains control.
The ``Connection`` class packs a coroutine source (a generator) and a sink. Together these objects provide a plug-board interface for the scheduler and a hierarchy of workers.

Now, when the scheduler calls ``sink.send(...)``, the job is pushed onto the queue that is created in ``single_worker()``. When the scheduler iterates over the results, ``get_result()`` feeds it results that it computes itself (through ``run_job``).

The Python queue is thread-safe. We may call ``jobs.source()`` in a different thread in another worker. This worker then safely pulls jobs from the same queue.

The Threaded worker
~~~~~~~~~~~~~~~~~~~

To have several workers run in tandem we need to keep a result queue in addition to the job queue. In the next sequence diagram we see how any number of threads are completely decoupled from the thread that manages the scheduling.

.. figure:: _static/images/sd-threaded.svg
    :alt: sequence diagram, multiple threads
    :align: center

    Sequence diagram where the actual job execution is deferred to one or more additional threads.

In Python source this looks as follows:

::

    def threaded_worker(n_threads):
        """Sets up a number of threads, each polling for jobs."""
        job_q = IOQueue()
        result_q = IOQueue()

        worker_connection = QueueConnection(job_q, result_q)
        scheduler_connection = QueueConnection(result_q, job_q)

        def worker(source, sink):
            for key, job in source:
                sink.send((key, run_job(job)))

        for i in range(n_threads):
            t = threading.Thread(
                target=worker,
                args=worker_connection.setup())

            t.daemon = True
            t.start()

        return scheduler_connection


The Hybrid worker
~~~~~~~~~~~~~~~~~

.. figure:: _static/images/sd-hybrid.svg
    :alt: sequence diagram, hybrid model
    :align: center

    Sequence diagram where the jobs get dispatched, each to a worker selected by a dispatcher.

Remote workers
--------------

Xenon
~~~~~

Fireworks
~~~~~~~~~
Cooking of Noodles (library docs)
=================================

The cooking of good Noodles can be tricky. We try to make it as easy as possible, but to write good Noodles you need to settle in a *functional style* of programming. The functions you design cannot write to some global state, or modify its arguments and expect these modifications to persist throughout the program. This is not a restriction of Noodles itself, this is a fundamental principle that applies to all possible frameworks for parallel and distributed programming. So get used to it!

Every function call in Noodles (that is, calls to scheduled function) can be visualised as a node in a call graph. You should be able to draw this graph conceptually when designing the program. Luckily there is (almost) always a way to write down non-functional code in a functional way.

.. NOTE:: Golden Rule: if you modify something, return it.


Call by value
-------------

Suppose we have the following program

::

    from noodles import (schedule, run_single)

    @schedule
    def double(x):
        return x['value'] * 2

    @schedule
    def add(x, y):
        return x + y

    a = {'value': 4}
    b = double(a)
    a['value'] = 5
    c = double(a)
    d = add(b, c)

    print(run_single(d))

If this were undecorated Python, the answer would be 18. However, the computation of this answer depends on the time-dependency of the Python interpreter. In Python, dictionaries are passed by reference. The promised object `b` then contains a reference to the dictionary in `a`. If we then change the value in this dictionary, the call producing the value of `b` is retroactively changed to double the value 5 instead of 4.

If Noodles is to evaluate this program correctly it needs to :py:func:`deepcopy` every argument to a scheduled function. There is another way to have the same semantics produce a correct result. This is by making `a` a promised object in the first place. The third solution is to teach your user *functional programming*.
Deep copying function arguments can result in a significant performance penalty on the side of the job scheduler. In most applications that we target this is not the bottle neck.

Since we aim for the maximum ease of use for the end-user, we chose to enable call-by-value by default.


Monads (sort of)
----------------

We still have ways to do object oriented programming and assignments. The :py:class:`PromisedObject` class has several magic methods overloaded to translate to functional equivalents.

Member assignment
~~~~~~~~~~~~~~~~~

Especially member assignment is treated in a particular way. Suppose ``a`` is a :py:class:`PromisedObject`, then the statement

::

    a.b = 3

is (conceptually) transformed into

::

    a = _setattr(a, 'b', 3)

where :py:func:`_setattr` is a scheduled function. The :py:class:`PromisedObject` contains a representation of the complete workflow representing the computation to get to the value of `a`. In member assignment, this workflow is replaced with the new workflow containing this last instruction.

This is not a recommended way of programming. Every assignment results in a nested function call. The `statefulness` of the program is then implemented in the composition of functions, similar to how other functional languages do it using `monads`. It results in sequential code that will not parallelise so well.

Other magic methods
~~~~~~~~~~~~~~~~~~~

Next to member assignment, we also (obviously) support member reference, method function call and object function call (with `__call__`).


Storable
--------



Serialisation
-------------
.. Noodles documentation master file, created by
   sphinx-quickstart on Wed Nov 11 13:52:27 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Noodles's documentation!
===================================

Introduction
------------
Often, a computer program can be sped up by executing parts of its code *in
parallel* (simultaneously), as opposed to *synchronously* (one part after
another).

A simple example may be where you assign two variables, as follows ``a = 2 * i``
and ``b = 3 * i``. Either statement is only dependent on ``i``, but whether you
assign ``a`` before ``b`` or vice versa, does not matter for how your program
works. Whenever this is the case, there is potential to speed up a program,
because the assignment of ``a`` and ``b`` could be done in parallel, using
multiple cores on your computer's CPU. Obviously, for simple assignments like
``a = 2 * i``, there is not much time to be gained, but what if ``a`` is the
result of a time-consuming function, e.g. ``a = very_difficult_function(i)``?
And what if your program makes many calls to that function, e.g. ``list_of_a =
[very_difficult_function(i) for i in list_of_i]``? The potential speed-up could
be tremendous.

So, parallel execution of computer programs is great for improving performance,
but how do you tell the computer which parts should be executed in parallel, and
which parts should be executed synchronously? How do you identify the order in
which to execute each part, since the optimal order may be different from the
order in which the parts appear in your program. These questions quickly become
nearly impossible to answer as your program grows and changes during
development. Because of this, many developers accept the slow execution of their
program only because it saves them from the headaches associated with keeping
track of which parts of their program depend on which other parts.

Enter Noodles.

Noodles is a Python package that can automatically construct a *callgraph*
for a given Python program, listing exactly which parts depend on which parts.
Moreover, Noodles can subsequently use the callgraph to execute code in parallel
on your local machine using multiple cores. If you so choose, you can even
configure Noodles such that it will execute the code remotely, for example on a
big compute node in a cluster computer.

Copyright & Licence
-------------------

Noodles 0.3.0 is copyright by the *Netherlands eScience Center (NLeSC)* and released under the Apache v2 License.

See http://www.esciencecenter.nl for more information on the NLeSC.

Installation
------------

.. WARNING:: We don't support Python versions lower than 3.5.

The core of Noodles runs on **Python 3.5** and above. To run Noodles on your own machine, no extra dependencies are required. It is advised to install Noodles in a virtualenv. If you want support for `Xenon`_, install `pyxenon`_ too.

.. code-block:: bash

    # create the virtualenv
    virtualenv -p python3 <venv-dir>
    . <venv-dir>/bin/activate

    # install noodles
    pip install noodles

Noodles has several optional dependencies. To be able to use the Xenon job scheduler, install Noodles with::

    pip install noodles[xenon]

The provenance/caching feature needs TinyDB installed::

    pip install noodles[prov]

To be able to run the unit tests::

    pip install noodles[test]

Documentation Contents
======================

.. toctree::
    :maxdepth: 2

    Introduction <self>
    eating
    cooking
    tutorials
    implementation


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _Xenon: http://nlesc.github.io/Xenon/
.. _pyxenon: http://github.com/NLeSC/pyxenon
.. _`generating SSH keys`: https://help.github.com/articles/generating-ssh-keys/
.. _`decorators`: https://www.thecodeship.com/patterns/guide-to-python-function-decorators/
.. highlight:: python
    :linenothreshold: 5

Boil: the source
~~~~~~~~~~~~~~~~

.. literalinclude:: ../../examples/boil/boil
    :linenos:
